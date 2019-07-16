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
    EXPECT_LE(std::abs((log_start - sat.geph_.toe).toSec()), 30);
  }
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
