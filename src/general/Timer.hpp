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
