#include <Eigen/Core>

#include "GNSS_utils/satellite.h"
#include "GNSS_utils/wgs84.h"

using namespace Eigen;

const double Satellite::GM_EARTH = 3.986005e14;
const double Satellite::OMEGA_EARTH = 7.2921151467e-5;
const double Satellite::PI = 3.1415926535898;
const double Satellite::C_LIGHT = 299792458.0;
const double Satellite::MAXDTOE = 7200.0; // max time difference to GPS Toe (s)
const double Satellite::FREQL1 = 1.57542e9;
const double Satellite::LAMBDA_L1 = Satellite::C_LIGHT / Satellite::FREQL1;



Satellite::Satellite(int id, int idx)
{
    id_ = id;
    idx_ = idx;
    memset(&eph_, 0, sizeof(eph_t));
    t.week = 0;
    t.tow_sec = 0;
}

Satellite::Satellite(const eph_t& eph, int idx)
{
    id_ = eph.sat;
    idx_ = idx;
    addEphemeris(eph);
    t.week = 0;
    t.tow_sec = 0;
}


void Satellite::addEphemeris(const eph_t &_eph)
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
void Satellite::computeMeasurement(const GTime& rec_time, const Vector3d& rec_pos, const Vector3d& receiver_vel, const Vector2d& clk_bias, Vector3d& z )
{
    Vector3d sat_pos_tr, sat_vel_tr;
    Vector2d sat_clk;
    computePositionVelocityClock(rec_time, sat_pos_tr, sat_vel_tr, sat_clk);
    double range = (sat_pos_tr - rec_pos).norm();
    double sagnac = OMEGA_EARTH * (sat_pos_tr.x()*rec_pos.y() - sat_pos_tr.y()*rec_pos.x())/C_LIGHT;
    range += sagnac;
    double tau = range / C_LIGHT;  // Time it took for observation to get here

    // extrapolate satellite position backwards in time
    Vector3d sat_pos_ts = sat_pos_tr - sat_vel_tr * tau;

    // Earth rotation correction. The change in velocity can be neglected.
    Vector3d earth_rot = sat_pos_ts.cross(e_z * OMEGA_EARTH * tau);
    sat_pos_ts += earth_rot;

    // Re-calculate the line-of-sight vector with the adjusted position
    Vector3d los = sat_pos_ts - rec_pos;
    range = los.norm();
    sagnac = OMEGA_EARTH * (sat_pos_ts.x()*rec_pos.y() - sat_pos_ts.y()*rec_pos.x())/C_LIGHT;
    range += sagnac;

    // adjust range by the satellite clock offset
    z(0) = range + C_LIGHT * (clk_bias(0) - sat_clk(0));

    // compute relative velocity between receiver and satellite, adjusted by the clock drift rate
    z(1) = (vel - receiver_vel).dot(los/range)
            + OMEGA_EARTH / C_LIGHT * e_z.dot(sat_vel.cross(receiver_vel))
            + C_LIGHT * (clk_bias(1) - sat_clk(1));
    printf("coriolis prange = %f\n", OMEGA_EARTH / C_LIGHT * e_z.dot(vel.cross(receiver_vel)));

    // Don't calculate lla if we are at (0, 0, 0)
    if (rec_pos.norm() > 0)
    {
        // Compute Azimuth and Elevation to satellite
        Vector2d az_el;
        los2azimuthElevation(rec_pos, los, az_el);
        Vector3d lla = WGS84::ecef2lla(rec_pos);

        // Compute and incorporate ionospheric delay
        double ion_delay = ionosphericDelay(rec_time, lla, az_el);
        double trop_delay = troposphericDelay(rec_time, lla, az_el);
        z(0) += ion_delay + trop_delay;
    }

    z(2) = z(0) / LAMBDA_L1;

    return;
}

Vector2d Satellite::los2azimuthElevation(const Vector3d &receiver_pos_ecef, const Vector3d &los_ecef) const
{
    Vector2d az_el;
    los2azimuthElevation(receiver_pos_ecef, los_ecef, az_el);
    return az_el;
}

Vector2d Satellite::azimuthElevation(const GTime& t, const Vector3d &rec_pos_ecef)
{
  update(t);
  Vector3d los_ecef = pos - rec_pos_ecef;
  Vector2d az_el;
  los2azimuthElevation(rec_pos_ecef, los_ecef, az_el);
  return az_el;
}

void Satellite::los2azimuthElevation(const Vector3d& receiver_pos_ecef, const Vector3d& los_ecef, Vector2d& az_el) const
{
    xform::Xformd x_e2n = WGS84::x_ecef2ned(receiver_pos_ecef);
    Vector3d los_ned = x_e2n.q().rotp(los_ecef.normalized());
    quat::Quatd q_los = quat::Quatd::from_two_unit_vectors(e_x, los_ned);
    az_el(0) = q_los.yaw();
    az_el(1) = q_los.pitch();
}

double Satellite::troposphericDelay(const GTime &t, const Vector3d &pos, const Vector2d &azel) const
{
    // Saastamoninen Troposphere Model
    const double temp0=15.0; /* temparature at sea level */
    double humi = 0.7; // relative humidity

    if (pos[2] < -100.0 || 1E4<pos[2] || azel[1] <= 0)
        return 0.0;

    /* standard atmosphere */
    double hgt = pos[2]<0.0 ? 0.0 : pos[2];

    double pres=1013.25*std::pow(1.0-2.2557E-5*hgt,5.2568);
    double temp=temp0-6.5E-3*hgt+273.16;
    double e=6.108*humi*std::exp((17.15*temp-4684.0)/(temp-38.45));

    double z=PI/2.0-azel[1];
    double trph=0.0022768*pres/(1.0-0.00266*std::cos(2.0*pos[0])-0.00028*hgt/1E3)/std::cos(z);
    double trpw=0.002277*(1255.0/temp+0.05)*e/std::cos(z);
    return trph+trpw;
}


double Satellite::ionosphericDelay(const GTime& gtime, const Vector3d& lla, const Vector2d& az_el) const
{
    // Klobuchar Algorithm:
    // https://gssc.esa.int/navipedia/index.php/Klobuchar_Ionospheric_Model

    const double ion[]={ /* 2004/1/1 */
                         0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
                         0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
                       };

    // Earth-Centered Angle (Elevation in Semicircles)
    double psi = 0.0137 / (az_el(1)/M_PI + 0.11) - 0.022;

    // latitude of the ionosphere pierce point (IPP)
    double phi_I = lla(0)/M_PI + psi * std::cos(az_el(0));
    phi_I = phi_I > 0.416 ? 0.416 : phi_I < -0.416 ? -0.416 : phi_I;

    // longitude of IPP
    double lambda_I = lla(1)/M_PI + (psi * std::sin(az_el(0))/cos(phi_I*M_PI));

    // geomagnetic latitude of the IPP
    double phi_m = phi_I + 0.064 * std::cos((lambda_I - 1.617) * M_PI);

    // local time at hte IPP
    double t= 43200*lambda_I + gtime.tow_sec;
    t -= std::floor(t/86400.0) * 86400.0;

    // Amplitude of Ionospheric Delay
    double Amp = ion[0] + phi_m * (ion[1] + phi_m * (ion[2] + phi_m * ion[3]));
    Amp = Amp < 0 ? 0 : Amp;

    // Period of Ionospheric Delay
    double Per = ion[4] + phi_m * (ion[5] + phi_m * (ion[6] + phi_m * ion[7]));
    Per = Per < 72000 ? 72000 : Per;

    // Phase of Ionospheric Delay
    double X_I = 2.0 * M_PI * (t - 50400.0) / Per;

    // Slant Factor
    double F = 1.0 + 16.0 * pow((0.53 - az_el(1)/M_PI), 3.0);


    // Compute Ionospheric Time Delay (meters)
    if (std::abs(X_I) <= 1.57)
    {
        double X2 = X_I*X_I;
        double X4 = X2*X2;
        return C_LIGHT * F * (5e-9 + Amp * (1.0 - X2/2.0 + X4/24.0));
    }
    else
        return C_LIGHT * F * 5e-9;
}

double Satellite::selectEphemeris(const GTime &time) const
{
    // find the closest ephemeris
    return (time - eph_.toe).toSec();
//    assert(eph_.size() > 0);
//    while (1)
//    {
//        dt = (time - eph.toe).toSec();
//        if (dt < 0 && (closest_eph_idx_ > 0))
//        {
//            double prev_dt = (time - eph_[closest_eph_idx_-1].toe).toSec();
//            if (std::abs(prev_dt) < std::abs(dt))
//            {
//                dt = prev_dt;
//                closest_eph_idx_ -= 1;
//                eph = &eph_[closest_eph_idx_];
//            }
//            else
//                break;
//        }

//        else if (dt > 0 && (closest_eph_idx_ < eph_.size()-1))
//        {
//            double next_dt = (time - eph_[closest_eph_idx_+1].toe).toSec();
//            if (std::abs(next_dt) < std::abs(dt))
//            {
//                dt = next_dt;
//                closest_eph_idx_ += 1;
//                eph = &eph_[closest_eph_idx_];
//            }
//            else
//                break;
//        }
//        else
//            break;
//    }
//    return dt;
}

void Satellite::update(const GTime &time)
{
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


void Satellite::computePositionVelocityClock(const GTime& time, const Ref<Vector3d> &_pos, const Ref<Vector3d> &_vel, const Ref<Vector2d>& _clock)
{
    if (time != t)
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

void Satellite::readFromRawFile(std::string filename)
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

Obs::Obs()
{
    sat_idx = -1;
}
bool Obs::operator <(const Obs& other)
{
    return sat_idx < other.sat_idx;
}
