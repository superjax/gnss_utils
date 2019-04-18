#pragma once
#include <gtest/gtest.h>

#include <Eigen/Core>

#include "GNSS_utils/satellite.h"
#include "GNSS_utils/gtime.h"
#include "GNSS_utils/datetime.h"

static const double RAD2DEG = 180.0 / M_PI;
static const double DEG2RAD = M_PI / 180.0;


#define ASSERT_MAT_NEAR(mat1, mat2, tol)\
do{  \
    ASSERT_EQ((mat1).rows(), (mat2).rows()); \
    ASSERT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat1).rows(); r++) { \
        for (int c = 0; c < (mat1).cols(); c++) {\
            ASSERT_NEAR((mat1)(r,c), (mat2)(r,c), tol); \
        }\
    } \
} while(0)

#define EXPECT_MAT_NEAR(mat1, mat2, tol)\
do{  \
    EXPECT_EQ((mat1).rows(), (mat2).rows()); \
    EXPECT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat2).rows(); r++) { \
        for (int c = 0; c < (mat2).cols(); c++) {\
            EXPECT_NEAR((mat1)(r,c), (mat2)(r,c), tol); \
        }\
    } \
} while(0)


#define ASSERT_MAT_EQ(mat1, mat2)\
do{  \
    ASSERT_EQ((mat1).rows(), (mat2).rows()); \
    ASSERT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat1).rows(); r++) { \
        for (int c = 0; c < (mat1).cols(); c++) {\
            ASSERT_FLOAT_EQ((mat1)(r,c), (mat2)(r,c)); \
        }\
    } \
} while(0)

#define EXPECT_MAT_EQ(mat1, mat2)\
do{  \
    EXPECT_EQ((mat1).rows(), (mat2).rows()); \
    EXPECT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat1).rows(); r++) { \
        for (int c = 0; c < (mat1).cols(); c++) {\
            EXPECT_FLOAT_EQ((mat1)(r,c), (mat2)(r,c)); \
        }\
    } \
} while(0)

#define EXPECT_MAT_NE(mat1, mat2)\
do{  \
    EXPECT_EQ((mat1).rows(), (mat2).rows()); \
    EXPECT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat1).rows(); r++) { \
        for (int c = 0; c < (mat1).cols(); c++) {\
            EXPECT_NE((mat1)(r,c), (mat2)(r,c)); \
        }\
    }\
} while(0)


#define ASSERT_MAT_NE(mat1, mat2)\
do{  \
    ASSERT_EQ((mat1).rows(), (mat2).rows()); \
    ASSERT_EQ((mat1).cols(), (mat2).cols()); \
    for (int r = 0; r < (mat1).rows(); r++) { \
        for (int c = 0; c < (mat1).cols(); c++) {\
            ASSERT_NE((mat1)(r,c), (mat2)(r,c)); \
        }\
    } \
} while(0)


typedef struct
{
    GTime g;
    double range; // pseudorange
    double rate;
    double d; // geometric distance
    double azel[2];
    double iono_delay;
} range_t;
typedef struct
{
    int enable;
    int vflg;
    double alpha0,alpha1,alpha2,alpha3;
    double beta0,beta1,beta2,beta3;
} ionoutc_t;
void eph2pos(const GTime& t, const eph_t *eph, Eigen::Vector3d& pos, double *dts);
double ionmodel(const GTime& t, const double *pos, const double *azel);
double ionosphericDelay(const ionoutc_t *ionoutc, GTime g, double *llh, double *azel);
void computeRange(range_t *rho, Satellite &eph, ionoutc_t *ionoutc, GTime g, Vector3d& xyz);
