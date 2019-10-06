#include <gtest/gtest.h>

#include "gnss_utils/glo_sat.h"
#include "gnss_utils/test_common.h"
#include "gnss_utils/wgs84.h"
#include "gnss_utils/datetime.h"

using namespace Eigen;
using namespace gnss_utils;

//class TestGloSat: public ::testing::Test
//{
//  GTime time;
//  geph_t geph_;
//  GloSat sat;
//};

TEST (GloSat, ReadFromFile)
{
  std::vector<int> sat_ids = {42, 43};
  for (int i = 0; i < sat_ids.size(); i++)
  {
    GloSat sat(sat_ids[i],i);
    sat.readFromRawFile(GNSS_UTILS_DIR"/sample/geph.dat");
    ASSERT_GT(sat.geph_.toe.week, 0);
  }
}

TEST (GloSat, ReadFromFileCheckTime)
{
  GTime log_start = GTime::fromUTC(1561988124,0.550);
  std::vector<int> sat_ids = {42, 43};
  for (int i = 0; i < sat_ids.size(); i++)
  {
    GloSat sat(sat_ids[i],i);
    sat.readFromRawFile(GNSS_UTILS_DIR"/sample/geph.dat");
    EXPECT_LE(std::abs((log_start - sat.geph_.toe).toSec()), Satellite::MAXDTOE);
  }
}

TEST (GloSat, OrbitEquation)
{
  GloSat sat(42, 0);
  Vector6d x;
  x.head<3>() = Map<Vector3d>(sat.geph_.pos);
  x.tail<3>() = Map<Vector3d>(sat.geph_.vel);

  GloSat::orbit(x, Map<Vector3d>(sat.geph_.acc));
  Vector3d x_oracle;
  double clk;
  double var;
  geph2pos(sat.geph_.toe, &sat.geph_, x_oracle, &clk, &var);

  EXPECT_MAT_EQ(x_oracle, x.head<3>());
}

TEST (GloSat, ReadFromFileCheckPVT)
{
  GTime log_start = GTime::fromUTC(1561988124,0.550);
  std::vector<int> sat_ids = {42, 43};
  for (int i = 0; i < sat_ids.size(); i++)
  {
    GloSat sat(sat_ids[i],i);
    sat.readFromRawFile(GNSS_UTILS_DIR"/sample/geph.dat");

    Vector3d pos, vel;
    Vector3d oracle_pos, oracle_pos_p, oracle_pos_m;
    double oracle_clock, oracle_clock_p, oracle_clock_m;
    double var;
    Vector2d clock;
    sat.computePVT(log_start, pos, vel, clock);
    geph2pos(log_start, &sat.geph_, oracle_pos, &oracle_clock, &var);

    // numerically differentiate for velocity
    double eps = 1e-3;
    GTime tp = log_start+eps;
    GTime tm = log_start-eps;
    geph2pos(tp, &sat.geph_, oracle_pos_p, &oracle_clock_p, &var);
    geph2pos(tm, &sat.geph_, oracle_pos_m, &oracle_clock_m, &var);
    Vector3d oracle_vel = (oracle_pos_p - oracle_pos_m) / (tp - tm).toSec();
    double oracle_clock_rate = (oracle_clock_p - oracle_clock_m) / (tp-tm).toSec();

    EXPECT_MAT_NEAR(pos, oracle_pos, 1e-8);
    EXPECT_MAT_NEAR(vel, oracle_vel, 1e-4);
    EXPECT_NEAR(clock(0), oracle_clock, 1e-16);
    EXPECT_NEAR(clock(1), oracle_clock_rate, 1e-10);
  }
}

TEST (GloSat, ReadFromFileCheckAzEl)
{
  std::vector<int> sat_ids     = {42, 43};
  // Truth Data gathered from https://in-the-sky.org/satmap_radar.php?year=2018&month=11&day=5
  // It's not especially accurate, but it's relatively close
  std::vector<double> sat_dist = {25051000, 25277000, 22050000, 22132000, 20560000, 20874000, 23237000, 21193000, 20854000};
  std::vector<double> sat_az   = {-62, -122, 104, -81, -14, -89, -58, 166, 40};
  std::vector<double> sat_el   = {6, 5, 38, 33, 78, 54, 26, 48, 60};
  std::vector<GPSSat> GPSSats;

  GTime log_start = GTime::fromUTC(1561988124,0.550);
  std::cout << DateTime(log_start) << std::endl;

  Vector3d rec_pos {-1798904.13, -4532227.1 ,  4099781.95};

  for (int i = 0; i < sat_ids.size(); i++)
  {
    GPSSat sat(sat_ids[i], i);
    sat.readFromRawFile(GNSS_UTILS_DIR"/sample/eph.dat");

    Vector3d pos, vel;
    Vector2d clock, az_el;
    sat.computePVT(log_start, pos, vel, clock);
    Vector3d los_ecef = pos - rec_pos;
    sat.los2AzEl(rec_pos, los_ecef, az_el);
    std::cout << "id: " << sat_ids[i] << std::endl;
    std::cout << "azel: " << az_el.transpose() << std::endl;
    std::cout << "dist: " << los_ecef.norm() << std::endl;
    EXPECT_NEAR((pos-rec_pos).norm(), sat_dist[i], 10000);
    EXPECT_NEAR(az_el(0)*RAD2DEG, sat_az[i], 1);
    EXPECT_NEAR(az_el(1)*RAD2DEG, sat_el[i], 1);
  }
}
