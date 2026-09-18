//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdThreeCenterElectronRepulsionGradientRows.hpp"

#include "SimdThreeCenterElectronRepulsionGeom010BufferRows.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

// NOTE: this wrapper keeps a translation unit which uses the table from having
// to include the header which declares it. That mattered while every set
// declared its table under one name and a unit which saw two of them was ill
// formed; the tables are named apart now, and this stays because callers use
// its name rather than because it is still load bearing.

auto
geom_010_buffer_rows(const int a_angular_momentum, const int b_angular_momentum, const int c_angular_momentum) -> size_t
{
    return simdt3ceri::number_of_geom_010_buffer_rows(a_angular_momentum, b_angular_momentum, c_angular_momentum);
}

}  // namespace simdt3cerigrad
