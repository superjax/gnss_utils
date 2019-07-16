#pragma once

#include <cstdint>
#include <Eigen/Core>

#include "gnss_utils/gtime.h"

namespace gnss_utils
{

struct Obs
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    GTime t;
    uint8_t sat_idx; // index in sats_ SatVec
    uint8_t sat;
    uint8_t rcv;
    uint8_t SNR;
    uint8_t LLI; // loss-of-lock indicator
    uint8_t code;
    uint8_t qualL; // carrier phase cov
    uint8_t qualP; // psuedorange cov
    Eigen::Vector3d z; // [prange, doppler, cphase]

    Obs();
    bool operator < (const Obs& other);
};

}
