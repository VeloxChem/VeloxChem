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

#ifndef Prescreener_hpp
#define Prescreener_hpp

#include <array>
#include <cstdint>
#include <tuple>
#include <vector>

#include "GtoBlock.hpp"

namespace prescr {  // prescr namespace

/**
 Gets grid box dimension.

 @param gridBlockPosition the displacement of grid points in this box.
 @param nGridPoints the number of grid points in this box.
 @param xcoords the X coordinates of grid points.
 @param ycoords the Y coordinates of grid points.
 @param zcoords the Z coordinates of grid points.
 @return grid box dimension as (xmin, ymin, zmin, xmax, ymax, zmax).
 */
auto getGridBoxDimension(const int gridBlockPosition,
                         const int nGridPoints,
                         const double* xcoords,
                         const double* ycoords,
                         const double* zcoords) -> std::array<double, 6>;

/**
 Prescreens GTOs block for a grid box.

 @param gtoBlock the GTO block.
 @param gtoDeriv the level of GTO derivative.
 @param gtoThreshold the screening threshold for GTO.
 @param boxDimension the dimension of the grid box.
 @return the mask indices for CGTOs and the pre-screened indices for AOs.
 */
auto preScreenGtoBlock(const CGtoBlock& gtoBlock, const int gtoDeriv, const double gtoThreshold, const std::array<double, 6>& boxDimension)
    -> std::tuple<std::vector<int>, std::vector<int>>;

}  // namespace prescr

#endif /* Prescreener_hpp */
