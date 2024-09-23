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

#ifndef GtoValuesRecP_hpp
#define GtoValuesRecP_hpp

#include <cstddef>
#include <vector>

#include "GtoBlock.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

/// @brief Computes LDA basis function values on given grid for S type basis functions.
/// @param gto_block The basis functions block.
/// @param grid_coords_x The vector of Cartesian X coordinates of grid.
/// @param grid_coords_y The vector of Cartesian Y coordinates of grid.
/// @param grid_coords_z The vector of Cartesian Z coordinates of grid.
/// @param gtos_mask The mask for basis functions (1 evaluate, 0 skip).
/// @return The matrix with LDA basis function values on grid points.
auto get_lda_values_rec_p(const CGtoBlock&            gto_block,
                          const std::vector<double>&  grid_coords_x,
                          const std::vector<double>&  grid_coords_y,
                          const std::vector<double>&  grid_coords_z,
                          const std::vector<int>& gtos_mask) -> CMatrix;

/// @brief Computes GGA basis function values on given grid for S type basis functions.
/// @param gto_block The basis functions block.
/// @param grid_coords_x The vector of Cartesian X coordinates of grid.
/// @param grid_coords_y The vector of Cartesian Y coordinates of grid.
/// @param grid_coords_z The vector of Cartesian Z coordinates of grid.
/// @param gtos_mask The mask for basis functions (1 evaluate, 0 skip).
/// @return The matrix with GGA basis function values on grid points.
auto get_gga_values_rec_p(const CGtoBlock&            gto_block,
                          const std::vector<double>&  grid_coords_x,
                          const std::vector<double>&  grid_coords_y,
                          const std::vector<double>&  grid_coords_z,
                          const std::vector<int>& gtos_mask) -> CMatrix;

}  // namespace gtoval

#endif /* GtoValuesRecP_hpp */
