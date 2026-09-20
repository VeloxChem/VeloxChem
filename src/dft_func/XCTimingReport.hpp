//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

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
namespace xcprof {

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
auto report(const std::string    &what,
            const CMultiTimer    &timer,
            const std::vector<CMultiTimer> &omptimers,
            const size_t          boxes) -> void;

}  // namespace xcprof

#endif /* XCTimingReport_hpp */
