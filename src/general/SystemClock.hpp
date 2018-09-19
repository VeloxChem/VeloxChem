//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef SystemClock_hpp
#define SystemClock_hpp

#include <chrono>
#include <string>
#include <ctime>

/**
 Class CSystemClock stores time information and provides set of methods for
 manipulating time and dates.
 
 @author Z. Rinkevicius
 */
class CSystemClock
{
    /**
     The reference time clock.
     */
    std::chrono::system_clock::time_point _refTime; // reference time clock

    /**
     Converts std::time_t data structure to std::tm data structure, which
     contains local time zone information.

     @param time the std::time_t data structure.
     @return the std::tm data structure.
     */
    std::tm _getLocalTime(const std::time_t& time) const;

    /**
     Converts std:tm data structure to string in format "hours:minutes:seconds
     on year-month-day".

     @param time the std::tm data structure.
     @return the "hours:minutes:seconds on year-month-day" string.
     */
    std::string _date(const std::tm* time) const;

public:

    /**
     Creates a system clock object and initializes reference clock.
     */
    CSystemClock();

    /**
     Creates a system clock object by copying other system clock object.
     
     @param source the system clock object.
     */
    CSystemClock(const CSystemClock& source);

    /**
     Creates a system clock object by by moving other system clock object.
     
     @param source the system clock object.
     */
    CSystemClock(CSystemClock&& source) noexcept;

    /**
     Destroys a system clock object.
     */
    ~CSystemClock();

    /**
     Assigns a system clock object by copying other system clock object.
     
     @param source the system clock object.
     */
    CSystemClock& operator=(const CSystemClock& source);

    /**
     Assigns a system clock object by moving other system clock object.
     
     @param source the system clock object.
     */
    CSystemClock& operator=(CSystemClock&& source) noexcept;

    /**
     Resets reference clock in system clock object.
     */
    void restart();

    /**
     Converts reference clock time to date string in format
     "hours:minutes:seconds on year-month-day".

     @return the reference time as formatted string.
     */
    std::string getStartDate() const;

    /**
     Converts current clock time to date string in format
     "hours:minutes:seconds on year-month-day".

     @return the reference time as formatted string.
     */
    std::string getCurrentDate() const;

    /**
     Converts elapsed time (current clock time - reference clock time) to string
     in format "hours h. minutes min. seconds sec.".

     @return the elapsed time as formatted string.
     */
    std::string getElapsedTime() const;

    /**
     Converts elapsed time (current clock time - reference clock time) to
     seconds.

     @return the elapsed time in seconds.
     */
    double getElapsedTimeInSeconds() const;
};

#endif /* SystemClock_hpp */
