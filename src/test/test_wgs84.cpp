#include <gtest/gtest.h>

#include "gnss_utils/wgs84.h"
#include "gnss_utils/test_common.h"

using namespace Eigen;

TEST (Gnss, lla2ecef)
{
    Vector3d lla = {40.246184 * DEG2RAD , -111.647769 * DEG2RAD, 1387.997511}; // BYU Campus
    Vector3d ecef_known = {-1798810.23, -4532232.54, 4099784.74};

    Vector3d ecef_calc = WGS84::lla2ecef(lla);

    ASSERT_MAT_NEAR(ecef_known, ecef_calc, 1e-2);
}

TEST (Gnss, ecef2lla)
{
    Vector3d ecef = {-1798810.23, -4532232.54, 4099784.74};
    Vector3d lla_known = {40.246184 * DEG2RAD , -111.647769 * DEG2RAD, 1387.998309};

    Vector3d lla_calc = WGS84::ecef2lla(ecef);

    ASSERT_MAT_NEAR(lla_known, lla_calc, 1e-6);
}

TEST (Gnss, ecef2lla2ecef)
{
    Vector3d ecef = {-1798810.23, -4532232.54, 4099784.74};
    Vector3d lla_known = {40.246184 * DEG2RAD , -111.647769 * DEG2RAD, 1387.998309};

    Vector3d lla_calc = WGS84::ecef2lla(ecef);
    Vector3d ecef_calc = WGS84::lla2ecef(lla_calc);

    ASSERT_MAT_NEAR(ecef_calc, ecef, 1e-6);
}

TEST (Gnss, x_ned2ecef)
{
    Vector3d lla0 = {40.247082 * DEG2RAD, -111.647776 * DEG2RAD, 1387.998309};
    Vector3d ecef0 = WGS84::lla2ecef(lla0);

    Vector3d ned1 = {-54.976484, 1.276565, 0.000237};
    Vector3d lla1 = {40.246587 * DEG2RAD, -111.647761 * DEG2RAD, 1387.998309};
    Vector3d ecef1 = WGS84::lla2ecef(lla1);

    xform::Xformd x_e2n = WGS84::x_ecef2ned(ecef0);
    Vector3d ecef_hat = x_e2n.transforma(ned1);
    Vector3d ned1_hat = x_e2n.transformp(ecef1);

    EXPECT_MAT_NEAR(ecef_hat, ecef1, 1e-6);
    EXPECT_MAT_NEAR(ned1_hat, ned1, 1e-3);
}

TEST (Gnss, ecef2ned_check_axes)
{
    Vector3d lla0 = {40.247082 * DEG2RAD, -111.647776 * DEG2RAD, 1387.998309};
    Vector3d ecef0 = WGS84::lla2ecef(lla0);
    xform::Xformd x_e2n = WGS84::x_ecef2ned(ecef0);

    double sp = std::sin(lla0(0));
    double cp = std::cos(lla0(0));
    double sl = std::sin(lla0(1));
    double cl = std::cos(lla0(1));

    Matrix3d R; // https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates
    R << -sp*cl, -sl, -cp*cl,
         -sp*sl,  cl, -cp*sl,
          cp,     0,  -sp;
    EXPECT_MAT_NEAR(x_e2n.q().R(), R.transpose(), 1e-8);

    Vector3d E_r_N_E = 1.0*ecef0;
    E_r_N_E /= E_r_N_E.stableNorm();
    Vector3d E_z_N = x_e2n.q().rota(e_z);
    EXPECT_MAT_NEAR(-1.0 * E_z_N, E_r_N_E, 3e-3);
}

TEST (Gnss, ned2ecef)
{
    Vector3d lla0 = {40.247082 * DEG2RAD, -111.647776 * DEG2RAD, 1387.998309};
    Vector3d ecef0 = WGS84::lla2ecef(lla0);

    Vector3d ned1 = {-54.976484, 1.276565, 0.000237};
    Vector3d lla1 = {40.246587 * DEG2RAD, -111.647761 * DEG2RAD, 1387.998309};
    Vector3d ecef1 = WGS84::lla2ecef(lla1);

    xform::Xformd x_e2n = WGS84::x_ecef2ned(ecef0);
    Vector3d ecef_hat = WGS84::ned2ecef(x_e2n, ned1);
    Vector3d ned1_hat = WGS84::ecef2ned(x_e2n, ecef1);

    EXPECT_MAT_NEAR(ecef_hat, ecef1, 1e-6);
    EXPECT_MAT_NEAR(ned1_hat, ned1, 1e-6);
}

TEST (Gnss, lla2ned)
{
    Vector3d lla0 = {40.247082 * DEG2RAD, -111.647776 * DEG2RAD, 1387.998309};
    Vector3d lla1 = {40.246587 * DEG2RAD, -111.647761 * DEG2RAD, 1387.998309};
    Vector3d ned1 = {-54.976484, 1.276565, 0.000237};

    Vector3d ned1_hat = WGS84::lla2ned(lla0, lla1);
    Vector3d lla1_hat = WGS84::ned2lla(lla0, ned1);

    EXPECT_MAT_NEAR(lla1_hat, lla1, 1e-6);
    EXPECT_MAT_NEAR(ned1_hat, ned1, 1e-6);
}
