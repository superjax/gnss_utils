#include "gnss_utils/gps_sat.h"

using namespace Eigen;

namespace gnss_utils
{

const double GPSSat::MAXDTOE = 7200.0; // max time difference to GPS Toe (s)
const double GPSSat::FREQL1 = 1.57542e9;
const double GPSSat::LAMBDA_L1 = GPSSat::C_LIGHT / GPSSat::FREQL1;

GPSSat::GPSSat(int id, int idx) :
  Satellite(id, idx)
{
    memset(&eph_, 0, sizeof(eph_t));
}

GPSSat::GPSSat(const eph_t& eph, int idx) :
  Satellite(eph.sat, idx)
{
    addEphemeris(eph);
}


void GPSSat::addEphemeris(const eph_t &_eph)
{
    assert( _eph.sat == id_);
//    ASSERT((eph_.size() > 0) ? eph.toe >= eph_.back().toe : true,
//           "tried to push ephemeris out of order");

    assert(_eph.toe.tow_sec <= DateTime::SECONDS_IN_WEEK);
    assert(_eph.toe.week <= 1000000);
    eph_ = _eph;
    t.week = 0;
    t.tow_sec = 0;

//    if (eph_.size() > 0 && (eph.toe == eph_.back().toe))
//    {
//        eph_.back() = eph;
//        if (eph_.size() == 1)
//        {
//            eph = &eph_[0];
//        }
//    }
//    else
//    {
//        eph_.push_back(eph);
//    }

//    if (!eph)
//    {
//        eph = &eph_[0];
//    }
}

//#define DBG(x) printf(#x": %6.6f\n", x); std::cout << std::flush;
void GPSSat::computeMeas(const GTime& rec_time, const Vector3d& rec_pos, const Vector3d& rec_vel, const Vector2d& clk_bias, Vector3d& z )
{
    update(rec_time);
    double range = (pos - rec_pos).norm();
    double sagnac = OMEGA_EARTH * (pos.x()*rec_pos.y() - pos.y()*rec_pos.x())/C_LIGHT;
    range += sagnac;
    double tau = range / C_LIGHT;  // Time it took for observation to get here

    // extrapolate satellite position backwards in time
    Vector3d sat_pos_ts = pos - vel * tau;

    // Earth rotation correction. The change in velocity can be neglected.
    Vector3d earth_rot = sat_pos_ts.cross(e_z * OMEGA_EARTH * tau);
    sat_pos_ts += earth_rot;

    // Re-calculate the line-of-sight vector with the adjusted position
    Vector3d los = sat_pos_ts - rec_pos;
    range = los.norm();
    sagnac = OMEGA_EARTH * (sat_pos_ts.x()*rec_pos.y() - sat_pos_ts.y()*rec_pos.x())/C_LIGHT;
    range += sagnac;

    // adjust range by the satellite clock offset
    z(0) = range + C_LIGHT * (clk_bias(0) - clk(0));

    // compute relative velocity between receiver and satellite, adjusted by the clock drift rate and coriolis
    z(1) = (vel - rec_vel).dot(los/range)
            + GPSSat::OMEGA_EARTH / GPSSat::C_LIGHT * (vel[1]*rec_pos[0] + pos[1]*rec_vel[0] - vel[0]*rec_pos[1] - pos[0]*rec_vel[1]) /// TODO check this math
            + C_LIGHT * (clk_bias(1) - clk(1));

    // Don't calculate lla if we are at (0, 0, 0)
    if (rec_pos.norm() > 0)
    {
        // Compute Azimuth and Elevation to satellite
        Vector2d az_el;
        los2AzEl(rec_pos, los, az_el);
        Vector3d lla = WGS84::ecef2lla(rec_pos);

        // Compute and incorporate ionospheric delay
        double ion_delay = ionDelay(rec_time, lla, az_el);
        double trop_delay = tropDelay(rec_time, lla, az_el);
        z(0) += ion_delay + trop_delay;
    }

    z(2) = z(0) / LAMBDA_L1;

    return;
}

void GPSSat::update(const GTime &time)
{
    if (time == t)
        return;

    double dt = selectEphemeris(time);

    if (dt > MAXDTOE)
        return;

    // https://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
    double n0 = std::sqrt(GM_EARTH/(eph_.A*eph_.A*eph_.A));
    double mkdot = n0 + eph_.deln;
    double mk = eph_.M0 + mkdot*dt;

    int i = 0;
    double ek_prev;
    double ek = mk;
    while (std::abs(ek - ek_prev) > 1e-13 && i < 30)
    {
        ek_prev = ek;
        ek -= (ek - eph_.e * std::sin(ek) - mk) / (1.0 - eph_.e*std::cos(ek));
        i++;
    }
    double sek = std::sin(ek);
    double cek = std::cos(ek);


    double ekdot = mkdot/(1.0 - eph_.e * cek);

    double tak = std::atan2(std::sqrt(1.0-eph_.e*eph_.e) * sek, cek - eph_.e);
    double takdot = sek*ekdot*(1.0+eph_.e*std::cos(tak)) / (std::sin(tak)*(1.0-eph_.e*cek));

    double phik = tak + eph_.omg;
    double sphik2 = std::sin(2.0 * phik);
    double cphik2 = std::cos(2.0 * phik);
    double corr_u = eph_.cus * sphik2 + eph_.cuc * cphik2;
    double corr_r = eph_.crs * sphik2 + eph_.crc * cphik2;
    double corr_i = eph_.cis * sphik2 + eph_.cic * cphik2;
    double uk = phik + corr_u;
    double rk = eph_.A*(1.0 - eph_.e*cek) + corr_r;
    double ik = eph_.i0 + eph_.idot*dt + corr_i;

    double s2uk = std::sin(2.0*uk);
    double c2uk = std::cos(2.0*uk);

    double ukdot = takdot + 2.0 * (eph_.cus * c2uk - eph_.cuc*s2uk) * takdot;
    double rkdot = eph_.A * eph_.e * sek * mkdot / (1.0 - eph_.e * cek) + 2.0 * (eph_.crs * c2uk - eph_.crc * s2uk) * takdot;
    double ikdot = eph_.idot + (eph_.cis * c2uk - eph_.cic * s2uk) * 2.0 * takdot;

    double cuk = std::cos(uk);
    double suk = std::sin(uk);

    double xpk = rk * cuk;
    double ypk = rk * suk;

    double xpkdot = rkdot * cuk - ypk * ukdot;
    double ypkdot = rkdot * suk + xpk * ukdot;

    double omegak = eph_.OMG0 + (eph_.OMGd - OMEGA_EARTH) * dt - OMEGA_EARTH * eph_.toes;
    double omegakdot = eph_.OMGd - OMEGA_EARTH;

    double cwk = std::cos(omegak);
    double swk = std::sin(omegak);
    double cik = std::cos(ik);
    double sik = std::sin(ik);

    pos.x() = xpk * cwk - ypk * swk * cik;
    pos.y() = xpk * swk + ypk * cwk * cik;
    pos.z() = ypk * sik;

    vel.x() = ( xpkdot - ypk*cik*omegakdot )*cwk
            - ( xpk*omegakdot + ypkdot*cik - ypk*sik*ikdot )*swk;
    vel.y() = ( xpkdot - ypk*cik*omegakdot )*swk
            + ( xpk*omegakdot + ypkdot*cik - ypk*sik*ikdot )*cwk;
    vel.z() = ypkdot*sik + ypk*cik*ikdot;

    dt = (time - eph_.toc).toSec();
    double dts = eph_.f0 + eph_.f1*dt + eph_.f2*dt*dt;

    // Correct for relativistic effects on the satellite clock
    dts -= 2.0*std::sqrt(GM_EARTH * eph_.A) * eph_.e * sek/(C_LIGHT * C_LIGHT);

    clk(0) = dts; // satellite clock bias
    clk(1) = eph_.f1 + eph_.f2*dt; // satellite drift rate
    t = time;
}


void GPSSat::computePVT(const GTime& time, const Ref<Vector3d> &_pos, const Ref<Vector3d> &_vel, const Ref<Vector2d>& _clock)
{
    update(time);

    // const-cast hackery to get around Ref
    // https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
    Ref<Vector3d> pos_ref = const_cast<Ref<Vector3d>&>(_pos);
    Ref<Vector3d> vel_ref = const_cast<Ref<Vector3d>&>(_vel);
    Ref<Vector2d> clock_ref = const_cast<Ref<Vector2d>&>(_clock);

    pos_ref = pos;
    vel_ref = vel;
    clock_ref = clk;
    return;
}

void GPSSat::readFromRawFile(const std::string &filename)
{
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    const int size = sizeof(eph_t);
    char buf[size];
    eph_t* eph = (eph_t*)buf;

    if (!file.is_open())
    {
        throw std::runtime_error(std::string("unable to open ") + filename);
    }

    while (!file.eof())
    {
        file.read(buf, size);
        if (!file)
            break;
        if (eph->sat == id_)
        {
            eph->toe = GTime::fromUTC(eph->toe.week, eph->toe.tow_sec);
            eph->toc = GTime::fromUTC(eph->toc.week, eph->toc.tow_sec);
            eph->ttr = GTime::fromUTC(eph->ttr.week, eph->ttr.tow_sec);
            addEphemeris(*eph);
        }
    }
}

double GPSSat::selectEphemeris(const GTime &time) const
{
  return (time - eph_.toe).toSec();
}

}

