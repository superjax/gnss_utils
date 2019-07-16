#include "gnss_utils/obs.h"

namespace gnss_utils
{

Obs::Obs()
{
    sat_idx = -1;
}
bool Obs::operator <(const Obs& other)
{
    return sat_idx < other.sat_idx;
}

}
