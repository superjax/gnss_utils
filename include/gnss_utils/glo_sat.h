#pragma once


#include "gnss_utils/satellite.h"

namespace gnss_utils
{

typedef struct
{
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    GTime toe;        /* epoch of epherides (gpst) */
    GTime tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

class GloSat : public Satellite
{
public:
    static const double TSTEP;
    static const double ERREPH_GLO;
    static const double OMGE_GLO;
    static const double MU_GLO;
    static const double J2_GLO;
    static const double RE_GLO;

    GloSat(int id, int idx);
    GloSat(const geph_t& _geph, int idx);

    void update(const GTime &time) override;


    void computePVT(const GTime &g,
                    const Eigen::Ref<Eigen::Vector3d> &_pos,
                    const Eigen::Ref<Eigen::Vector3d> &_vel,
                    const Eigen::Ref<Eigen::Vector2d> &_clock) override;

    void computeMeas(const GTime& rec_time,
                     const Eigen::Vector3d& receiver_pos,
                     const Eigen::Vector3d &rec_vel,
                     const Eigen::Vector2d &clk_bias,
                     Eigen::Vector3d &z) override;
    double selectEphemeris(const GTime &time) const;
    void readFromRawFile(const std::string& filename) override;

    void addEphemeris(const geph_t& _eph);

    static Vector6d orbit(const Vector6d& x, const Eigen::Vector3d& acc)
    {
      auto p = x.head<3>();
      auto v = x.tail<3>();
      Vector6d xdot;
      auto pdot = xdot.head<3>();
      auto vdot = xdot.tail<3>();
      const double r2 = p.transpose() * p;
      const double r3 = r2 * std::sqrt(r2);
      const double r5 = r2 * r3;
      const double omg2 = GloSat::OMGE_GLO * GloSat::OMGE_GLO;
      const double re2 = GloSat::RE_GLO * GloSat::RE_GLO;

      double a = 1.5 * GloSat::J2_GLO * GloSat::MU_GLO * re2 / r5; //3/2*J2*mu*Ae^2/r^5
      printf("a = %f\n", a);
      double b = 5.0 * p.z() * p.z() / r2;
      printf("b = %f\n", b);
      double c = -GloSat::MU_GLO / r3 - a * (1.0 - b);
      printf("c = %f\n", c);
      pdot = v;
      vdot.x() = (c + omg2)  * p.x() + 2.0*GloSat::OMGE_GLO * v.y() + acc.x();
      vdot.y() = (c + omg2)  * p.y() - 2.0*GloSat::OMGE_GLO * v.x() + acc.y();
      vdot.z() = (c - 2.0*a) * p.z() + acc.z();
      return xdot;
    }
    geph_t geph_ = { 0 };
};

}
