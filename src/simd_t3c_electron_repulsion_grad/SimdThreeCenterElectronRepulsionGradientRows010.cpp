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

auto
geom_010_buffer_rows(const int a_angular_momentum, const int b_angular_momentum, const int c_angular_momentum) -> size_t
{
    return simdt3ceri::number_of_buffer_rows(a_angular_momentum, b_angular_momentum, c_angular_momentum);
}

}  // namespace simdt3cerigrad
