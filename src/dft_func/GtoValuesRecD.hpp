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

}  // namespace gtoval

#endif /* GtoValuesRecD_hpp */
