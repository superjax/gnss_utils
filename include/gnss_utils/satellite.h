#pragma once

#include <Eigen/Core>

#include "gnss_utils/datetime.h"
#include "gnss_utils/gtime.h"

#include "geometry/xform.h"

namespace gnss_utils
{

class Satellite
{

public:
  static const double GM_EARTH;
  static const double OMEGA_EARTH;
  static const double PI;
  static const double C_LIGHT;

  GTime t;
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Vector2d clk;
  int id_;
  int idx_;

  Satellite(int id, int idx);
  virtual void update(const GTime& time) {}
  virtual void computePVT(const GTime &g, const Eigen::Ref<Eigen::Vector3d> &pos,
                          const Eigen::Ref<Eigen::Vector3d> &vel, const Eigen::Ref<Eigen::Vector2d> &clock) {}
  virtual void computeMeas(const GTime& rec_time, const Eigen::Vector3d& receiver_pos,
                   const Eigen::Vector3d &rec_vel, const Eigen::Vector2d &clk_bias,
                           Eigen::Vector3d &z){}
  virtual double selectEphemeris(const GTime& time) const {}
  virtual void readFromRawFile(const std::string& filename) {}

  void los2AzEl(const Eigen::Vector3d& receiver_pos_ecef,
                const Eigen::Vector3d& los_ecef,
                Eigen::Vector2d& az_el) const;
  Eigen::Vector2d los2AzEl(const Eigen::Vector3d& receiver_pos_ecef,
                           const Eigen::Vector3d& los_ecef) const;
  double ionDelay(const GTime &t, const Eigen::Vector3d& lla, const Eigen::Vector2d& az_el) const;
  double tropDelay(const GTime &t, const Eigen::Vector3d& lla, const Eigen::Vector2d& az_el) const;

  Eigen::Vector2d azEl(const GTime &t, const Eigen::Vector3d& rec_pos_ecef);
};

}
