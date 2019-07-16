#include <fstream>
#include <functional>

#include "gnss_utils/glo_sat.h"

using namespace Eigen;

namespace gnss_utils
{

const double GloSat::TSTEP = 60.0;
const double GloSat::ERREPH_GLO = 5.0;
const double GloSat::OMGE_GLO = 7.292115E-5;
const double GloSat::MU_GLO = 3.9860044E14;
const double GloSat::J2_GLO = 1.0826257E-3;
const double GloSat::RE_GLO = 6378136.0;

static Eigen::VectorXd RK4(std::function<Eigen::VectorXd(const Eigen::VectorXd& x, const Eigen::VectorXd& u)> f,
                double dt, const Eigen::VectorXd& x0, const Eigen::VectorXd& u)
{
  VectorXd k1 = f(x0, u);
  VectorXd k2 = f(x0+k1*dt/2.0, u);
  VectorXd k3 = f(x0+k2*dt/2.0, u);
  VectorXd k4 = f(x0+k3*dt, u);
  return x0 + (k1 + 2.0*k2 + 2.0*k3 + k4) * (dt / 6.0);
}

GloSat::GloSat(int id, int idx) :
  Satellite(id, idx)
{
  memset(&geph_, 0, sizeof(geph_t));
}

GloSat::GloSat(const geph_t &_geph, int idx) :
  Satellite(_geph.sat, idx)
{
  geph_ = _geph;
}

void GloSat::readFromRawFile(const std::string &filename)
{
  std::ifstream file(filename, std::ios::in | std::ios::binary);
  const int size = sizeof(geph_t);
  char buf[size];
  geph_t* geph = (geph_t*)buf;

  if (!file.is_open())
  {
      throw std::runtime_error(std::string("unable to open ") + filename);
  }

  while (!file.eof())
  {
      file.read(buf, size);
      if (!file)
          break;
      if (geph->sat == id_)
      {
          addEphemeris(*geph);
      }
  }
}

void GloSat::addEphemeris(const geph_t &_eph)
{
  geph_ = _eph;
}


#define dbg(x) printf("%s:%d %s=%f\n",__FILE__,__LINE__,#x,x)
void GloSat::update(const GTime &time)
{
  if (time == t)
    return;

  double dt = selectEphemeris(time);
  dbg(dt);

  // Clock biases
  clk(0) = -geph_.taun + geph_.gamn*dt;
  clk(1) = geph_.gamn;

  Vector6d x;
  auto p = x.head<3>();
  auto v = x.tail<3>();
  p = Map<Vector3d>(geph_.pos);
  v = Map<Vector3d>(geph_.vel);

  double timestep = dt < 0.0 ? -GloSat::TSTEP : GloSat::TSTEP;
  while (std::abs(dt) > 1e-9)
  {
    if (std::abs(dt) < GloSat::TSTEP)
      timestep = dt;

    x = RK4(GloSat::orbit, timestep, x, Map<const Vector3d>(geph_.acc));

    dt -= timestep;
  }

  pos = x.head<3>();
  vel = x.tail<3>();
}

void GloSat::computePVT(const GTime &g, const Ref<Vector3d> &_pos, const Ref<Vector3d> &_vel,
                        const Ref<Vector2d> &_clock)
{
  update(g);

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

void GloSat::computeMeas(const GTime &rec_time,
                         const Eigen::Vector3d &receiver_pos,
                         const Eigen::Vector3d &rec_vel,
                         const Eigen::Vector2d &clk_bias,
                         Eigen::Vector3d &z)
{

}

double GloSat::selectEphemeris(const GTime &time) const
{
  return (time - geph_.toe).toSec();
}


}
