//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

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
auto getGridBoxDimension(const int64_t gridBlockPosition,
                         const int64_t nGridPoints,
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
auto preScreenGtoBlock(const CGtoBlock& gtoBlock, const int64_t gtoDeriv, const double gtoThreshold, const std::array<double, 6>& boxDimension)
    -> std::tuple<std::vector<int64_t>, std::vector<int64_t>>;

}  // namespace prescr

#endif /* Prescreener_hpp */
