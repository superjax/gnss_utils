#pragma once
#include <stdint.h>
#include <vector>
#include <fstream>

#include <Eigen/Core>

#ifndef NDEBUG
#define ASSERT(condition, ...) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << printf(__VA_ARGS__) << std::endl; \
            assert(condition); \
        } \
    } while (false)
#else
#   define ASSERT(...)
#endif

#include "GNSS_utils/datetime.h"
#include "GNSS_utils/gtime.h"

typedef struct {
  int32_t sat; // satellite number
  int32_t iode; // IODE Issue of Data, Ephemeris (ephemeris version)
  int32_t iodc; // IODC Issue of Data, Clock (clock version)
  int32_t sva; // SV accuracy (URA index) IRN-IS-200H p.97
  int32_t svh; // SV health GPS/QZS (0:ok)
  int32_t week; // GPS/QZS: gps week, GAL: galileo week (00=Invalid, 01 = P Code ON, 11 = C/A code ON, 11 = Invalid) // GPS/QZS: code on L2 * GAL/CMP: data sources
  int32_t code;  // GPS/QZS: L2 P data flag (indicates that the NAV data stream was commanded OFF on the P-code of the in-phase component of the L2 channel) *  CMP: nav type
  int32_t flag;
  GTime toe; // Toe
  GTime toc; // clock data reference time (s) (20.3.4.5)
  GTime ttr; // T_trans
  double A; // Semi-Major Axis m
  double e; // Eccentricity (no units)
  double i0; // Inclination Angle at Reference Time (rad)
  double OMG0; // Longitude of Ascending Node of Orbit Plane at Weekly Epoch (rad)
  double omg; // Argument of Perigee (rad)
  double M0; // Mean Anomaly at Reference Time (rad)
  double deln; // Mean Motion Difference From Computed Value (rad)
  double OMGd; // Rate of Right Ascension (rad/s)
  double idot; // Rate of Inclination Angle (rad/s)
  double crc; // Amplitude of the Cosine Harmonic Correction Term to the Orbit Radius
  double crs; // Amplitude of the Sine Harmonic Correction Term to the Orbit Radius (m)
  double cuc; // Amplitude of the Cosine Harmonic Correction Term to the Argument of Latitude (rad)
  double cus; // Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude (rad)
  double cic; // Amplitude of the Cosine Harmonic Correction Term to the Angle of Inclination (rad)
  double cis; // Amplitude of the Sine Harmonic Correction Term to the Angle of Inclination (rad)
  double toes; // Reference Time Ephemeris in week (s)
  double fit; // fit interval (h) (0: 4 hours, 1:greater than 4 hours)
  double f0; // SV clock parameters - af0
  double f1; // SV clock parameters - af1
  double f2; // SV clock parameters - af2 * * GPS/QZS:tgd[0]=TGD (IRN-IS-200H p.103) // group delay parameter * GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E * CMP    :tgd[0]=BGD1,tgd[1]=BGD
  double tgd[4];
  double Adot; // Adot for CNAV
  double ndot; // ndot for CNAV
} eph_t;

using namespace Eigen;

class Satellite
{
public:
    static const double GM_EARTH;
    static const double OMEGA_EARTH;
    static const double PI;
    static const double C_LIGHT;
    static const double MAXDTOE;
    static const double FREQL1;
    static const double LAMBDA_L1;

    Satellite(int id, int idx);
    Satellite(const eph_t& eph, int idx);
    void update(const GTime &time);
    void computePositionVelocityClock(const GTime &g, const Ref<Vector3d> &pos, const Ref<Vector3d> &vel, const Ref<Vector2d> &clock);
    void computeMeasurement(const GTime& rec_time, const Vector3d& receiver_pos, const Vector3d &receiver_vel, const Vector2d &clk_bias, Vector3d &z);
    void los2azimuthElevation(const Vector3d& receiver_pos_ecef, const Vector3d& los_ecef, Vector2d& az_el) const;
    Vector2d los2azimuthElevation(const Vector3d& receiver_pos_ecef, const Vector3d& los_ecef) const;
    double ionosphericDelay(const GTime &t, const Vector3d& lla, const Vector2d& az_el) const;
    double troposphericDelay(const GTime &t, const Vector3d& lla, const Vector2d& az_el) const;
    double selectEphemeris(const GTime& time) const;
    void readFromRawFile(std::string filename);
    void addEphemeris(const eph_t& eph_);
    Vector2d azimuthElevation(const GTime &t, const Vector3d& rec_pos_ecef);

    int id_;
    int idx_;
    eph_t eph_ = { 0 };
    GTime t;
    Vector3d pos;
    Vector3d vel;
    Vector2d clk;
};
