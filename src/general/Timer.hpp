//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef Timer_hpp
#define Timer_hpp

// to avoid MSVC error
// "error: capturing a structured binding is not yet supported in OpenMP"
#ifndef _MSC_VER
#include <chrono>
#endif

#include <string>

/**
 Class CTimer implements timer.
 */
class CTimer
{
    /**
     The time duration (accumulated).
     */
#ifndef _MSC_VER
    std::chrono::steady_clock::duration _duration;
#else
    double _duration;
#endif

    /**
     Whether the timer has started.
     */
    bool _started;

    /**
     The time point when the timer started.
     */
#ifndef _MSC_VER
    std::chrono::steady_clock::time_point _startTime;
#else
    double _startTime;
#endif

   public:
    /**
     Creates a timer object.
     */
    CTimer();

    /**
     Resets the timer.
     */
    void reset();

    /**
     Starts the timer.
     */
    void start();

    /**
     Stops the timer.
     */
    void stop();

    /**
     Gets elapsed time (accumulated).

     @return the elapsed time as a string.
     */
    std::string getElapsedTime() const;
};

#endif /* Timer_hpp */
