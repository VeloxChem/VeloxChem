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

#include "MultiTimer.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>

CMultiTimer::CMultiTimer()
{
}

void
CMultiTimer::clear()
{
    _labels.clear();

    _timers.clear();
}

void
CMultiTimer::start(const std::string& label)
{
    auto label_iter = std::find(_labels.begin(), _labels.end(), label);

    if (label_iter != _labels.end())
    {
        // found label

        int64_t index = label_iter - _labels.begin();

        _timers[index].start();
    }
    else
    {
        // add new label

        _labels.push_back(label);

        _timers.push_back(CTimer());

        _timers.back().start();
    }
}

void
CMultiTimer::stop(const std::string& label)
{
    auto label_iter = std::find(_labels.begin(), _labels.end(), label);

    if (label_iter != _labels.end())
    {
        // found label

        int64_t index = label_iter - _labels.begin();

        _timers[index].stop();
    }
}

std::string
CMultiTimer::getSummary() const
{
    std::stringstream ss;

    for (size_t timer_id = 0; timer_id < _timers.size(); timer_id++)
    {
        const auto& label = _labels[timer_id];

        const auto& timer = _timers[timer_id];

        ss << std::left << std::setw(25) << label << ": " << timer.getElapsedTime() << "\n";
    }

    return ss.str();
}
