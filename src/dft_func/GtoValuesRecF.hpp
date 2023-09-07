#ifndef GtoValuesRecF_hpp
#define GtoValuesRecF_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

/**
 Computes GTO values on given grid for F type GTOs.

 @param gto_block the GTOs block.
 @param n_points the number of grid points.
 @param grid_coords_x the pointer to Cartesian X coordinates of grid.
 @param grid_coords_y the pointer to Cartesian Y coordinates of grid.
 @param grid_coords_z the pointer to Cartesian Z coordinates of grid.
 @param gtos_mask the mask for GTOs (1 evaluate, 0 skip).
 */
auto getLdaValuesRecF(const CGtoBlock&            gto_block,
                      const int64_t               n_points,
                      const double*               grid_coords_x,
                      const double*               grid_coords_y,
                      const double*               grid_coords_z,
                      const std::vector<int64_t>& gtos_mask) -> CMatrix;

}  // namespace gtoval

#endif /* GtoValuesRecF_hpp */
