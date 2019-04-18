#include <cmath>

#include "GNSS_utils/datetime.h"
#include "GNSS_utils/gtime.h"

DateTime::DateTime()
{}

DateTime::DateTime(const GTime& g)
{
    (*this) = g.toDate();
}

DateTime& DateTime::operator= (const GTime& g)
{
    (*this) = g.toDate();
    return (*this);
}

GTime DateTime::toGTime() const
{
    GTime g;
    constexpr int doy[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
    int y;
    int d;
    int num_leap_days;

    y = year - 1980;

    // Compute the number of leap days since Jan 5/Jan 6, 1980.
    num_leap_days = y/4 + 1;
    if ((y % 4) == 0 && month <= 2)
        num_leap_days--;

    // Compute the number of days elapsed since Jan 5/Jan 6, 1980.
    d = y*365 + doy[month-1] + day + num_leap_days - 6;

    // Convert time to GPS weeks and seconds.
    g.week = d / 7;
    g.tow_sec = static_cast<double>((d % 7)*SECONDS_IN_DAY + hour*SECONDS_IN_HOUR + minute*SECONDS_IN_MINUTE + second + LEAP_SECONDS);
    return g;
}

std::ostream& operator<<(std::ostream& os, const DateTime& dt)
{
    os << dt.month << "/" << dt.day << "/" << dt.year << " "
       << dt.hour << ":" << dt.minute << ":" << dt.second;
    return os;
}
