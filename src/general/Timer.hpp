//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include <chrono>
#include <string>

/**
 Class CTimer implements timer.

 @author X. Li
 */
class CTimer
{
    /**
     The time duration (accumulated).
     */
    std::chrono::steady_clock::duration _duration;

    /**
     Whether the timer has started.
     */
    bool _started;

    /**
     The time point when the timer started.
     */
    std::chrono::steady_clock::time_point _startTime;

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
