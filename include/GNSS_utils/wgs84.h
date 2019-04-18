#pragma once

#include <Eigen/Core>
#include "geometry/xform.h"
#include "GNSS_utils/gtime.h"
#include "GNSS_utils/satellite.h"

using namespace Eigen;
using namespace xform;


struct WSG84
{
    typedef std::vector<Vector3d, aligned_allocator<Vector3d>> VecVec3;
    static constexpr double A = 6378137.0;       // WGS-84 Earth semimajor axis (m)
    static constexpr double B = 6356752.314245;  // Derived Earth semiminor axis (m)
    static constexpr double F = (A - B) / A;     // Ellipsoid Flatness
    static constexpr double F_INV = 1.0 / F;     // Inverse flattening
    static constexpr double A2 = A * A;
    static constexpr double B2 = B * B;
    static constexpr double E2 = F * (2 - F);    // Square of Eccentricity

    static Vector3d ecef2lla(const Vector3d& ecef)
    {
        Vector3d lla;
        ecef2lla(ecef, lla);
        return lla;
    }

    static void ecef2lla(const Vector3d& ecef, Vector3d& lla)
    {
        static const double e2 = F * (2.0 - F);

        double r2 = ecef.x()*ecef.x() + ecef.y()*ecef.y();
        double z=ecef.z();
        double v;
        double zk;
        do
        {
            zk = z;
            double sinp = z / std::sqrt(r2 + z*z);
            v = A / std::sqrt(1.0 - e2*sinp*sinp);
            z = ecef.z() + v*e2*sinp;
        }
        while (std::abs(z - zk) >= 1e-4);

        lla.x() = r2 > 1e-12 ? std::atan(z / std::sqrt(r2)) : (ecef.z() > 0.0 ? M_PI/2.0 : -M_PI/2.0);
        lla.y() = r2 > 1e-12 ? std::atan2(ecef.y(), ecef.x()) : 0.0;
        lla.z() = std::sqrt(r2+z*z) - v;
    }

    static Vector3d lla2ecef(const Vector3d& lla)
    {
        Vector3d ecef;
        lla2ecef(lla, ecef);
        return ecef;
    }

    static void lla2ecef(const Vector3d& lla, Vector3d& ecef)
    {
        double sinp=sin(lla[0]);
        double cosp=cos(lla[0]);
        double sinl=sin(lla[1]);
        double cosl=cos(lla[1]);
        double e2=F*(2.0-F);
        double v=A/sqrt(1.0-e2*sinp*sinp);

        ecef[0]=(v+lla[2])*cosp*cosl;
        ecef[1]=(v+lla[2])*cosp*sinl;
        ecef[2]=(v*(1.0-e2)+lla[2])*sinp;
    }

    static void x_ecef2ned(const Vector3d& ecef, xform::Xformd& X_e2n)
    {
        X_e2n.q() = q_e2n(ecef2lla(ecef));
        X_e2n.t() = ecef;
    }

    static xform::Xformd x_ecef2ned(const Vector3d& ecef)
    {
        xform::Xformd X_e2n;
        x_ecef2ned(ecef, X_e2n);
        return X_e2n;
    }

    static Vector3d ned2ecef(const xform::Xformd x_e2n, const Vector3d& ned)
    {
        return x_e2n.transforma(ned);
    }

    static void ned2ecef(const xform::Xformd x_e2n, const Vector3d& ned, Vector3d& ecef)
    {
        ecef = x_e2n.transforma(ned);
    }

    static Vector3d ecef2ned(const xform::Xformd x_e2n, const Vector3d& ecef)
    {
        return x_e2n.transformp(ecef);
    }

    static void ecef2ned(const xform::Xformd x_e2n, const Vector3d& ecef, Vector3d& ned)
    {
        ned = x_e2n.transformp(ecef);
    }

    static void lla2ned(const Vector3d& lla0, const Vector3d& lla, Vector3d& ned)
    {
        xform::Xformd x_e2n;
        x_e2n.q() = q_e2n(lla0);
        x_e2n.t() = lla2ecef(lla0);
        ecef2ned(x_e2n, lla2ecef(lla), ned);
    }

    static Vector3d lla2ned(const Vector3d& lla0, const Vector3d& lla)
    {
        Vector3d ned;
        lla2ned(lla0, lla, ned);
        return ned;
    }

    static void ned2lla(const Vector3d& lla0, const Vector3d& ned, Vector3d&lla)
    {
        xform::Xformd x_e2n;
        x_e2n.q() = q_e2n(lla0);
        x_e2n.t() = lla2ecef(lla0);
        ecef2lla(ned2ecef(x_e2n, ned), lla);
    }

    static Vector3d ned2lla(const Vector3d& lla0, const Vector3d& ned)
    {
        Vector3d lla;
        ned2lla(lla0, ned, lla);
        return lla;
    }

    static quat::Quatd q_e2n(const Vector3d& lla)
    {
        Quatd q1, q2;
        q1 = quat::Quatd::from_axis_angle(e_z, lla(1));
        q2 = quat::Quatd::from_axis_angle(e_y, -M_PI/2.0 - lla(0));
        return q1 * q2;
    }


    static bool pointPositioning(const GTime &t, const VecVec3 &z, std::vector<Satellite> &sats, Vector3d &xhat)
    {
      const int nsat = sats.size();
      MatrixXd A, b;
      A.resize(nsat, 4);
      b.resize(nsat, 1);
      Matrix<double, 4, 1> dx;
      GTime that = t;
      ColPivHouseholderQR<MatrixXd> solver;

      int iter = 0;
      do
      {
        iter++;
        int i = 0;
        for (Satellite sat : sats)
        {
          Vector3d sat_pos, sat_vel;
          Vector2d sat_clk_bias;
          sat.computePositionVelocityClock(t, sat_pos, sat_vel, sat_clk_bias);

          Vector3d zhat ;
          sat.computeMeasurement(t, xhat, Vector3d::Zero(), Vector2d::Zero(), zhat);
          b(i) = z[i](0) - zhat(0);

          A.block<1,3>(i,0) = (xhat - sat_pos).normalized().transpose();
          A(i,3) = Satellite::C_LIGHT;
          i++;
        }

        solver.compute(A);
        dx = solver.solve(b);

        xhat += dx.topRows<3>();
        that += dx(3);
      } while (dx.norm() > 1e-4 && iter < 10);
      return iter < 10;
    }
};

inline void printLla(const Vector3d& lla)
{
    std::cout << lla(0)*180.0/M_PI << ", " << lla(1)*180.0/M_PI << ", " << lla(2);
}


