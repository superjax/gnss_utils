#include <ros/ros.h>

#include "GNSS_utils/positioning_demo.h"
#include "GNSS_utils/wgs84.h"

PP_Demo::PP_Demo() :
    nh_(), nh_private_("~")
{
    obs_sub_ = nh_.subscribe("gps/obs", 1, &PP_Demo::obsVecCB, this);
    eph_sub_ = nh_.subscribe("gps/eph", 1, &PP_Demo::ephCB, this);
    gps_sub_ = nh_.subscribe("gps", 1, &PP_Demo::gpsCB, this);

    lla_pub_ = nh_.advertise<geometry_msgs::Vector3>("lla", 1);
    vel_pub_ = nh_.advertise<geometry_msgs::Vector3>("vel", 1);

    std::string log_name;
    nh_private_.param("min_sat_elevation", min_satellite_elevation_, 20*M_PI/180.0);
    nh_private_.param("log_name", log_name, std::string("/tmp/positioning_demo.bin"));
    log_.open(log_name);
    rec_pos_prev_.setZero();
    time_prev_.tow_sec = 0;
    time_prev_.week = 0;
    start_time_.tow_sec = 0;
    start_time_.week = 0;
}

void PP_Demo::gpsCB(const inertial_sense::GPSConstPtr &gps)
{
    ref_lla_ << gps->latitude, gps->longitude, gps->altitude;
    ref_vel_ecef_ << gps->velEcef.x, gps->velEcef.y, gps->velEcef.z;
}

void PP_Demo::obsVecCB(const inertial_sense::GNSSObsVecConstPtr &msg)
{
    obsVec new_obs;
    new_obs.reserve(msg->obs.size());

    if (start_time_.week == 0 && start_time_.tow_sec == 0)
        start_time_ = GTime::fromUTC(msg->header.stamp.sec, msg->header.stamp.nsec*1e-9);


    for (const auto& o : msg->obs)
    {
        Obs new_o;
        new_o.t = GTime::fromTime(o.time.time, o.time.sec);
        new_o.sat = o.sat;
        new_o.rcv = o.rcv;
        new_o.SNR = o.SNR;
        new_o.LLI = o.LLI;
        new_o.code = o.code;
        new_o.z << o.P, Satellite::LAMBDA_L1*o.D, o.L;
        new_obs.push_back(new_o);
    }
    filterObs(new_obs);
    if (obs_.size() > 7)
    {
        calcPosition();
        logPosition();
    }
}

void PP_Demo::ephCB(const inertial_sense::GNSSEphemerisConstPtr &eph)
{
    eph_t new_eph;
    new_eph.sat = eph->sat;
    new_eph.iode = eph->iode;
    new_eph.iodc = eph->iodc;
    new_eph.sva = eph->sva;
    new_eph.svh = eph->svh;
    new_eph.week = eph->week;
    new_eph.code = eph->code;
    new_eph.flag = eph->flag;
    new_eph.toe = GTime::fromTime(eph->toe.time, eph->toe.sec);
    new_eph.toc = GTime::fromTime(eph->toc.time, eph->toc.sec);
    new_eph.ttr = GTime::fromTime(eph->ttr.time, eph->ttr.sec);
    new_eph.A = eph->A;
    new_eph.e = eph->e;
    new_eph.i0 = eph->i0;
    new_eph.OMG0 = eph->OMG0;
    new_eph.omg = eph->omg;
    new_eph.M0 = eph->M0;
    new_eph.deln = eph->deln;
    new_eph.OMGd = eph->OMGd;
    new_eph.idot = eph->idot;
    new_eph.crc = eph->crc;
    new_eph.crs = eph->crs;
    new_eph.cuc = eph->cuc;
    new_eph.cus = eph->cus;
    new_eph.cic = eph->cic;
    new_eph.cis = eph->cis;
    new_eph.toes = eph->toes;
    new_eph.fit = eph->fit;
    new_eph.f0 = eph->f0;
    new_eph.f1 = eph->f1;
    new_eph.f2 = eph->f2;
    new_eph.tgd[0] = eph->tgd[0];
    new_eph.tgd[1] = eph->tgd[1];
    new_eph.tgd[2] = eph->tgd[2];
    new_eph.tgd[3] = eph->tgd[3];
    new_eph.Adot = eph->Adot;
    new_eph.ndot = eph->ndot;

    // Is this a new satellite?
    auto s = sats_.begin();
    while (s != sats_.end())
    {
        if (s->id_ == new_eph.sat)
            break;
        s++;
    }
    bool new_sat = (s == sats_.end());

    // Is this satellite high enough in the sky to care about?
    Satellite sat(new_eph, sats_.size());
    bool high_enough = true;
    if (rec_pos_prev_.x() > 0)
        high_enough = sat.azimuthElevation(time_prev_, rec_pos_prev_)(1) > min_satellite_elevation_;

    if (new_sat)
    {
        if (high_enough)
            sats_.push_back(sat);
    }
    else
    {
        if (high_enough)
            s->addEphemeris(new_eph);
        else
            sats_.erase(s);
    }
    refreshSatIdx();
}

void PP_Demo::refreshSatIdx()
{
    for (int i = 0; i < sats_.size(); i++)
    {
        sats_[i].idx_ = i;
    }
}

void PP_Demo::calcPosition()
{
    if (obs_.size() < 7)
    {
        ROS_INFO_THROTTLE(1, "Waiting for Ephemeris");
        return;
    }

    Vec8 xhat;
    xhat.setZero();
    xhat.topRows<3>() = rec_pos_prev_;
    lsPositioning(obs_[0].t, obs_, sats_, xhat);

    lla_ = WGS84::ecef2lla(xhat.topRows<3>());
    vel_ned_ = WGS84::x_ecef2ned(xhat.topRows<3>()).rotp(xhat.segment<3>(3));

    lla_.x() *= 180.0/M_PI;
    lla_.y() *= 180.0/M_PI;

    geometry_msgs::Vector3 lla_msg;
    lla_msg.x = lla_.x();
    lla_msg.y = lla_.y();
    lla_msg.z = lla_.z();
    lla_pub_.publish(lla_msg);
    geometry_msgs::Vector3 vel_msg;
    vel_msg.x = vel_ned_.x();
    vel_msg.y = vel_ned_.y();
    vel_msg.z = vel_ned_.z();
    vel_pub_.publish(vel_msg);
}

void PP_Demo::logPosition()
{
    log_.log((double)(time_prev_ - start_time_).toSec());
    Vector3d ref_vel_ned = WGS84::x_ecef2ned(WGS84::lla2ecef(lla_)).rotp(ref_vel_ecef_);
    std::cout << "lla: " << lla_.transpose() << ", ref: " << ref_lla_.transpose() << ",\t ";
    std::cout << "vel: " << vel_ned_.transpose() << ", ref: " << ref_vel_ned.transpose() << std::endl;
    log_.logVectors(lla_, ref_lla_, vel_ned_, ref_vel_ned);
}

void PP_Demo::filterObs(const obsVec &rawObs)
{
    obs_.clear();
    GTime time;
    time.tow_sec = 0.0;

    // Only add observations we have the satellites for
    // and make sure they all have the same time stamp
    for (auto& o : rawObs)
    {
        assert(time.tow_sec == 0.0 || time == o.t);
        time = o.t;
        int idx = getSatIdx(o.sat);
        if (idx >= 0)
        {
            obs_.push_back(o);
            obs_.back().sat_idx = idx;
        }
    }
    if (obs_.size() > 1)
        time_prev_ = obs_[0].t;
}

int PP_Demo::getSatIdx(int sat_id) const
{
    for (auto& s : sats_)
    {
        if (s.id_ == sat_id)
            return s.idx_;
    }
    return -1;
}

void PP_Demo::lsPositioning(const GTime &t, const obsVec &obs, satVec &sats, Vec8 &xhat)
{
    const int nobs = obs.size();
    MatrixXd A, b;
    A.setZero(nobs*2, 8);
    b.setZero(nobs*2, 1);
    Vector8d dx;
    ColPivHouseholderQR<MatrixXd> solver;

    int iter = 0;
    do
    {
        int i = 0;
        for (auto&& o : obs)
        {
            Satellite& sat(sats[o.sat_idx]);
            Vector3d sat_pos, sat_vel;
            Vector2d sat_clk_bias;
            auto phat = xhat.segment<3>(0);
            auto vhat = xhat.segment<3>(3);
            auto that = xhat.segment<2>(6);
            GTime tnew = t + that(0);
            sat.computePositionVelocityClock(tnew, sat_pos, sat_vel, sat_clk_bias);

            Vector3d zhat ;
            sat.computeMeasurement(tnew, phat, vhat, that, zhat);
            b(2*i) = o.z(0) - zhat(0);
            b(2*i + 1) = o.z(1) - zhat(1);

            Vector3d e_i = (sat_pos - phat).normalized();
            A.block<1,3>(2*i,0) = -e_i.transpose();
            A(2*i,6) = Satellite::C_LIGHT;
            A.block<1,3>(2*i+1,3) = -e_i.transpose();
            A(2*i+1,7) = Satellite::C_LIGHT;

            i++;
        }


        solver.compute(A);
        dx = solver.solve(b);

        xhat += dx;
    } while (dx.norm() > 1e-4 && ++iter < 10);
}



int main(int argc, char** argv)
{
    ros::init(argc, argv, "gnss_utils_demo");
    PP_Demo thing;
    ros::spin();
}
