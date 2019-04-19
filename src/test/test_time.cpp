#include <gtest/gtest.h>

#include <gnss_utils/gtime.h>
#include <gnss_utils/datetime.h>

using namespace gnss_utils;

TEST (Time, FromDateTimeKnown)
{
    DateTime dtime;
    dtime.year = 2019;
    dtime.month = 1;
    dtime.day = 13;
    dtime.hour = 8;
    dtime.minute = 57;
    dtime.second = 10;


    GTime gtime_known;
    gtime_known.week = 2036;
    gtime_known.tow_sec = 32248.0;

    GTime gtime_new = dtime;

    ASSERT_EQ(gtime_new.week, gtime_known.week);
    ASSERT_FLOAT_EQ(gtime_new.tow_sec, gtime_known.tow_sec);
}

TEST (Time, FromGTimeKnown)
{
    GTime gtime;
    gtime.week = 2036;
    gtime.tow_sec = 32248.0;

    DateTime dtime;
    dtime.year = 2019;
    dtime.month = 1;
    dtime.day = 13;
    dtime.hour = 8;
    dtime.minute = 57;
    dtime.second = 10;

    DateTime dtime_new = gtime;

    ASSERT_EQ(dtime.year, dtime_new.year);
    ASSERT_EQ(dtime.month, dtime_new.month);
    ASSERT_EQ(dtime.day, dtime_new.day);
    ASSERT_EQ(dtime.hour, dtime_new.hour);
    ASSERT_EQ(dtime.minute, dtime_new.minute);
    ASSERT_FLOAT_EQ(dtime.second, dtime_new.second);
}
