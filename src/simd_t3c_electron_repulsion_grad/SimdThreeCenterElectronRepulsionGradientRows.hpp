//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdThreeCenterElectronRepulsionGradientRows_hpp
#define SimdThreeCenterElectronRepulsionGradientRows_hpp

#include <cstddef>

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief The rows the buffer of a combination needs to differentiate the first
/// center on bra side.
/// @note The two generated tables declare one function of one name in one
/// namespace, so a unit which included both would not compile. Each is wrapped
/// in a unit of its own and named for the center it differentiates, which leaves
/// the generated files untouched and a regeneration a straight copy.
auto geom_100_buffer_rows(const int a_angular_momentum, const int b_angular_momentum, const int c_angular_momentum)
    -> size_t;

/// @brief The rows needed to differentiate the second center on bra side.
/// @note The two differ: differentiating one center or the other raises a
/// different shell by one, so the tables are not the same and the arena of a
/// driver which calls both takes the larger.
auto geom_010_buffer_rows(const int a_angular_momentum, const int b_angular_momentum, const int c_angular_momentum)
    -> size_t;

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGradientRows_hpp */
