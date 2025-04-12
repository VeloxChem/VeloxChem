//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
