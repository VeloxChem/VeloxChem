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

#ifndef XCTimingReport_hpp
#define XCTimingReport_hpp

#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>

#include "MultiTimer.hpp"

/// @brief The phases of an exchange correlation integration, reported when asked.
///
/// @note The integrators have carried named timers for a long while and printed
/// none of them: the calls which would have were commented out, one per function,
/// and would have written a block for every thread. This turns them on behind
/// VLX_XC_PROFILE and adds the threads together, which is what `getTimings` was put
/// there for.
///
/// @note The quadrature is the whole of a pure functional's Fock build once the
/// Coulomb matrix is fitted -- 86 to 92 per cent of it, as the RI-J tables record --
/// so where its own time goes is worth being able to ask.
namespace xcprof {  // xcprof namespace

/// @brief Whether the phases were asked for.
/// @note Read once. The check is on a path which runs every iteration.
inline auto
wanted() -> bool
{
    static const bool asked = (std::getenv("VLX_XC_PROFILE") != nullptr);

    return asked;
}

/// @brief Writes the phases of one integration and their share of it.
/// @param what Which integration it was, named in the heading.
/// @param timer The timer of the serial part, which holds the total.
/// @param omptimers One timer for each thread, whose labels are added together.
/// @param boxes The boxes of the grid, which are the tasks the threads are handed.
/// A grid with fewer boxes than there are threads cannot keep them busy, and the
/// phases will then fall short of the total by whatever the idle ones did not do.
/// @note The threads are summed and divided by their number, so a phase reads as
/// the wall time it would take if the threads shared it evenly. A phase which is
/// unbalanced therefore reads low, and the total is the honest number to check the
/// parts against.
auto report(const std::string             &what,
            const CMultiTimer             &timer,
            const std::vector<CMultiTimer> &omptimers,
            const size_t                   boxes) -> void;

/// @brief Writes what the prescreening left, and what it saved.
/// @param naos The basis functions of the molecule.
/// @param boxes The boxes of the grid.
/// @param points The grid points, summed over the boxes.
/// @param ao_sum The surviving functions, summed over the boxes.
/// @param ao_max The most any one box kept.
/// @param point_ao The surviving functions weighted by the points of their box.
/// @param point_ao_sq The same weighted by the square, which is the work the two
/// matrix phases really do.
/// @note The screening exists to hold the surviving count steady as the molecule
/// grows, which would make the quadrature linear in it. Whether it does is the
/// question these numbers answer: the work is points times the square of what
/// survives, so the points weighted figures are the ones to read and the plain
/// average is there only to show how far apart the two are.
auto report_blocks(const size_t naos,
                   const size_t boxes,
                   const size_t points,
                   const size_t ao_sum,
                   const size_t ao_max,
                   const size_t point_ao,
                   const size_t point_ao_sq) -> void;

}  // namespace xcprof

#endif /* XCTimingReport_hpp */
