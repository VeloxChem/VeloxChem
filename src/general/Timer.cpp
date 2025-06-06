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

#include "Timer.hpp"

#include <iomanip>
#include <sstream>

CTimer::CTimer()
{
    reset();
}

void
CTimer::reset()
{
#ifndef _MSC_VER
    auto t0 = std::chrono::steady_clock::now();

    _duration = t0 - t0;
#else
    _duration = 0.0;
#endif

    _started = false;
}

void
CTimer::start()
{
    if (!_started)
    {
#ifndef _MSC_VER
        _startTime = std::chrono::steady_clock::now();
#else
        _startTime = 0.0;
#endif

        _started = true;
    }
}

void
CTimer::stop()
{
    if (_started)
    {
#ifndef _MSC_VER
        _duration += std::chrono::steady_clock::now() - _startTime;
#endif

        _started = false;
    }
}

std::string
CTimer::getElapsedTime() const
{
    std::stringstream ss;

#ifndef _MSC_VER
    auto duration_double = std::chrono::duration_cast<std::chrono::duration<double>>(_duration);

    ss << std::fixed << std::setw(15) << std::setprecision(3) << duration_double.count() << " sec";
#else
    ss << "         N/A";
#endif

    return ss.str();
}
