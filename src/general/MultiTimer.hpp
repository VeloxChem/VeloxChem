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

#ifndef MultiTimer_hpp
#define MultiTimer_hpp

#include <string>
#include <vector>

#include "Timer.hpp"

/**
 Class CMultiTimer implements timer with multiple timing.

 @author X. Li
 */
class CMultiTimer
{
    /**
     The lables of the timers.
     */
    std::vector<std::string> _labels;

    /**
     The timers.
     */
    std::vector<CTimer> _timers;

   public:
    /**
     Creates a timer object.
     */
    CMultiTimer();

    /**
     Clears the multi-timer.
     */
    void clear();

    /**
     Starts a timer with a label.
     */
    void start(const std::string& label);

    /**
     Stops a timer with a label.
     */
    void stop(const std::string& label);

    /**
     Gets a summary of elapsed time (accumulated).

     @return the summary as a string.
     */
    std::string getSummary() const;
};

#endif /* MultiTimer_hpp */
