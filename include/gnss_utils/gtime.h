#pragma once
#include <stdint.h>
#include <iostream>

#include "gnss_utils/datetime.h"

namespace gnss_utils
{

class DateTime;
class GTime
{
public:
    int64_t week;
    double tow_sec;

    GTime();
    GTime(int64_t week, double tow_sec);
    GTime(const DateTime& t);
    GTime& operator= (const DateTime& t);
    GTime operator -(const GTime& gt2) const;
    GTime operator -(const double& sec) const;
    GTime operator +(const GTime& gt2) const;
    GTime operator +(const double& tow_sec) const;
    GTime& operator +=(const double& tow_sec);
    DateTime toDate() const;
    bool operator >(const GTime& g2) const;
    bool operator >=(const GTime& g2) const;
    bool operator <(const GTime& g2) const;
    bool operator <=(const GTime& g2) const;
    bool operator ==(const GTime& g2) const;
    bool operator !=(const GTime& g2) const;
    double toSec() const;
    static GTime fromUTC(int time_sec, double subsec);
    static GTime fromTime(int time, double sec);
};

inline GTime operator+ (const double& sec, const GTime& t)
{
    return t + sec;
}

inline std::ostream & operator << (std::ostream &out, const GTime &t)
{
    out << "[ " << t.week;
    out << ", " << t.tow_sec << " ]";
    return out;
}

}
