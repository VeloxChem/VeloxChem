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

#ifndef GtoValuesRecD_hpp
#define GtoValuesRecD_hpp

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
auto get_lda_values_rec_d(const CGtoBlock&            gto_block,
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
auto get_gga_values_rec_d(const CGtoBlock&            gto_block,
                          const std::vector<double>&  grid_coords_x,
                          const std::vector<double>&  grid_coords_y,
                          const std::vector<double>&  grid_coords_z,
                          const std::vector<int>& gtos_mask) -> CMatrix;

auto get_mgga_values_rec_d(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix;

auto get_3rd_order_values_rec_d(const CGtoBlock&            gto_block,
                                const std::vector<double>&  grid_coords_x,
                                const std::vector<double>&  grid_coords_y,
                                const std::vector<double>&  grid_coords_z,
                                const std::vector<int>&     gtos_mask) -> CMatrix;

}  // namespace gtoval

#endif /* GtoValuesRecD_hpp */
