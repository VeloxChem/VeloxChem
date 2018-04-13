//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SystemClock.hpp"

#include <sstream>

CSystemClock::CSystemClock()

    : _refTime(std::chrono::system_clock::now())
{

}

CSystemClock::CSystemClock(const CSystemClock& source)

    : _refTime(source._refTime)
{

}

CSystemClock::CSystemClock(CSystemClock&& source) noexcept

    : _refTime(std::move(source._refTime))
{

}

CSystemClock::~CSystemClock()
{

}

CSystemClock& CSystemClock::operator=(const CSystemClock& source)
{
    if (this == &source) return *this;

    _refTime = source._refTime;

    return *this;
}

CSystemClock& CSystemClock::operator=(CSystemClock&& source) noexcept
{
    if (this == &source) return *this;

    _refTime = std::move(source._refTime);

    return *this;
}

void CSystemClock::restart()
{
    _refTime = std::chrono::system_clock::now();
}

std::string CSystemClock::getStartDate() const
{
    auto mtime = _getLocalTime(std::chrono::system_clock::to_time_t(_refTime));

    return _date(&mtime);
}

std::string CSystemClock::getCurrentDate() const
{
    auto curtime = std::chrono::system_clock::now();

    auto mtime = _getLocalTime(std::chrono::system_clock::to_time_t(curtime));

    return _date(&mtime);
}

std::string CSystemClock::getElapsedTime() const
{
    auto tdiff = std::chrono::system_clock::now() - _refTime;

    auto thour = std::chrono::duration_cast<std::chrono::hours>(tdiff);

    auto mdiff = tdiff % std::chrono::hours(1);

    auto tmin = std::chrono::duration_cast<std::chrono::minutes>(mdiff);

    auto sdiff = tdiff % std::chrono::minutes(1);

    auto tsec = std::chrono::duration_cast<std::chrono::seconds>(sdiff);

    std::stringstream stream;

    stream << thour.count() << " h. " << tmin.count() << " min. ";

    stream << tsec.count() << " sec.";

    return stream.str();
}

double CSystemClock::getElapsedTimeInSeconds() const
{
    std::chrono::duration<double> tsec = std::chrono::system_clock::now()
                                       - _refTime;

    return tsec.count();
}

std::tm CSystemClock::_getLocalTime(const std::time_t& time) const
{
    std::tm ctime;

    localtime_r(&time, &ctime);

    return ctime;
}

std::string CSystemClock::_date(const std::tm* time) const
{
    char cstring[81];

    if (std::strftime(cstring, sizeof(cstring)-1, "%T on %F", time) <= 0)
    {
        *cstring = '\0';
    }

    return std::string(cstring);
}
