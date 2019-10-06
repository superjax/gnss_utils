#include <Eigen/Core>

#include "gnss_utils/satellite.h"
#include "gnss_utils/wgs84.h"

using namespace Eigen;

namespace gnss_utils
{

const double Satellite::MAXDTOE = 7200.0; // max time difference to Toe (s)
const double Satellite::GM_EARTH = 3.986005e14;
const double Satellite::OMEGA_EARTH = 7.2921151467e-5;
const double Satellite::PI = 3.1415926535898;
const double Satellite::C_LIGHT = 299792458.0;

Satellite::Satellite(int id, int idx)
{
    id_ = id;
    idx_ = idx;
    t.week = 0;
    t.tow_sec = 0;
}

Vector2d Satellite::los2AzEl(const Vector3d &receiver_pos_ecef, const Vector3d &los_ecef) const
{
    Vector2d az_el;
    los2AzEl(receiver_pos_ecef, los_ecef, az_el);
    return az_el;
}

Vector2d Satellite::azEl(const GTime& t, const Vector3d &rec_pos_ecef)
{
  update(t);
  Vector3d los_ecef = pos - rec_pos_ecef;
  Vector2d az_el;
  los2AzEl(rec_pos_ecef, los_ecef, az_el);
  return az_el;
}

void Satellite::los2AzEl(const Vector3d& receiver_pos_ecef, const Vector3d& los_ecef, Vector2d& az_el) const
{
    xform::Xformd x_e2n = WGS84::x_ecef2ned(receiver_pos_ecef);
    Vector3d los_ned = x_e2n.q().rotp(los_ecef.normalized());
    quat::Quatd q_los = quat::Quatd::from_two_unit_vectors(e_x, los_ned);
    az_el(0) = q_los.yaw();
    az_el(1) = q_los.pitch();
}

double Satellite::tropDelay(const GTime &t, const Vector3d &pos, const Vector2d &azel) const
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


double Satellite::ionDelay(const GTime& gtime, const Vector3d& lla, const Vector2d& az_el) const
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

}
