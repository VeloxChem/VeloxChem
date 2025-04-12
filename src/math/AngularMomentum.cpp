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

#include "AngularMomentum.hpp"

#include <array>

namespace angmom {  // angmom namespace

auto
getStringOfAngularMomentum(const int64_t angmom, const int64_t component) -> std::string
{
    if (angmom == 0) return std::string("s  ");

    if (angmom == 1)
    {
        const std::array<std::string, 3> clist({"p-1", "p0 ", "p+1"});

        return clist[component];
    }

    if (angmom == 2)
    {
        const std::array<std::string, 5> clist({"d-2", "d-1", "d0 ", "d+1", "d+2"});

        return clist[component];
    }

    if (angmom == 3)
    {
        const std::array<std::string, 7> clist({"f-3", "f-2", "f-1", "f0 ", "f+1", "f+2", "f+3"});

        return clist[component];
    }

    if (angmom == 4)
    {
        const std::array<std::string, 9> clist({"g-4", "g-3", "g-2", "g-1", "g0 ", "g+1", "g+2", "g+3", "g+4"});

        return clist[component];
    }

    if (angmom == 5)
    {
        const std::array<std::string, 11> clist({"h-5", "h-4", "h-3", "h-2", "h-1", "h0 ", "h+1", "h+2", "h+3", "h+4", "h+5"});

        return clist[component];
    }

    if (angmom == 6)
    {
        const std::array<std::string, 13> clist({"i-6", "i-5", "i-4", "i-3", "i-2", "i-1", "i0 ", "i+1", "i+2", "i+3", "i+4", "i+5", "i+6"});

        return clist[component];
    }

    return std::string();
}

}  // namespace angmom
