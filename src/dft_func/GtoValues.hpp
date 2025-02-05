//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef GtoValues_hpp
#define GtoValues_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

/**
 Computes LDA GTO values on given grid.

 @param gto_block the GTOs block.
 @param grid_coords_x the vector of Cartesian X coordinates of grid.
 @param grid_coords_y the vector of Cartesian Y coordinates of grid.
 @param grid_coords_z the vector of Cartesian Z coordinates of grid.
 @param gtos_mask the mask for GTOs (1 evaluate, 0 skip).
 */
auto get_gto_values_for_lda(const CGtoBlock&            gto_block,
                        const std::vector<double>&  grid_coords_x,
                        const std::vector<double>&  grid_coords_y,
                        const std::vector<double>&  grid_coords_z,
                        const std::vector<int>& gtos_mask) -> CMatrix;

/**
 Computes GGA GTO values on given grid.

 @param gto_block the GTOs block.
 @param grid_coords_x the vector of Cartesian X coordinates of grid.
 @param grid_coords_y the vector of Cartesian Y coordinates of grid.
 @param grid_coords_z the vector of Cartesian Z coordinates of grid.
 @param gtos_mask the mask for GTOs (1 evaluate, 0 skip).
 */
auto get_gto_values_for_gga(const CGtoBlock&            gto_block,
                        const std::vector<double>&  grid_coords_x,
                        const std::vector<double>&  grid_coords_y,
                        const std::vector<double>&  grid_coords_z,
                        const std::vector<int>& gtos_mask) -> CMatrix;

auto get_gto_values_for_mgga(const CGtoBlock&            gto_block,
                             const std::vector<double>&  grid_coords_x,
                             const std::vector<double>&  grid_coords_y,
                             const std::vector<double>&  grid_coords_z,
                             const std::vector<int>& gtos_mask) -> CMatrix;

auto get_gto_values_for_3rd_order(const CGtoBlock&            gto_block,
                                  const std::vector<double>&  grid_coords_x,
                                  const std::vector<double>&  grid_coords_y,
                                  const std::vector<double>&  grid_coords_z,
                                  const std::vector<int>& gtos_mask) -> CMatrix;

}  // namespace gtoval

#endif /* GtoValues_hpp */
