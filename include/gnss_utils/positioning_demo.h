#pragma once

#include <ros/ros.h>

#include <inertial_sense/GNSSObsVec.h>
#include <inertial_sense/GNSSEphemeris.h>
#include <inertial_sense/GlonassEphemeris.h>
#include <inertial_sense/GPS.h>
#include <geometry_msgs/Vector3.h>
#include <gnss_utils/SatInfo.h>


#include "gnss_utils/logger.h"
#include "gnss_utils/gps_sat.h"
#include "gnss_utils/glo_sat.h"
#include "gnss_utils/obs.h"
#include "gnss_utils/gtime.h"
#include "gnss_utils/datetime.h"

class PP_Demo
{
public:
    typedef std::vector<gnss_utils::Obs, Eigen::aligned_allocator<gnss_utils::Obs>> obsVec;
    typedef std::vector<gnss_utils::GPSSat, Eigen::aligned_allocator<gnss_utils::GPSSat>> satVec;
    typedef Eigen::Matrix<double, 8, 1> Vec8;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    PP_Demo();
    ~PP_Demo();
    void obsVecCB(const inertial_sense::GNSSObsVecConstPtr& msg);
    void ephCB(const inertial_sense::GNSSEphemerisConstPtr &eph);
    void gephCB(const inertial_sense::GlonassEphemerisConstPtr& geph);
    void gpsCB(const inertial_sense::GPSConstPtr& gps);

    void calcPosition();
    void logPosition();
    void refreshSatIdx();
    int getSatIdx(int sat_id) const;
    void filterObs(const obsVec& rawObs);
    void lsPositioning(const gnss_utils::GTime &t, const obsVec &obs, satVec &sats, Vec8 &xhat);

    satVec sats_;
    obsVec obs_;

    std::ofstream out_;

    ros::Subscriber obs_sub_;
    ros::Subscriber eph_sub_;
    ros::Subscriber geph_sub_;
    ros::Subscriber gps_sub_;
    std::vector<ros::Publisher> sat_pos_pub_;
    ros::Publisher lla_pub_;
    ros::Publisher vel_pub_;
    ros::NodeHandle nh_, nh_private_;

    Eigen::Vector3d rec_pos_prev_;
    gnss_utils::GTime time_prev_;
    double min_satellite_elevation_;

    Eigen::Vector3d ref_lla_, lla_;
    Eigen::Vector3d ref_vel_ecef_, vel_ned_;
    Logger log_;
    gnss_utils::GTime start_time_;

    int sat_to_track = -1;
    double prange_prev;
    gnss_utils::GTime t_prev;
};
