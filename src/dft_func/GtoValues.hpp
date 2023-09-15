#ifndef GtoValues_hpp
#define GtoValues_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

/**
 Computes GTO values on given grid.

 @param gto_block the GTOs block.
 @param grid_coords_x the vector of Cartesian X coordinates of grid.
 @param grid_coords_y the vector of Cartesian Y coordinates of grid.
 @param grid_coords_z the vector of Cartesian Z coordinates of grid.
 @param gtos_mask the mask for GTOs (1 evaluate, 0 skip).
 */
auto getGtoValuesForLda(const CGtoBlock&            gto_block,
                        const std::vector<double>&  grid_coords_x,
                        const std::vector<double>&  grid_coords_y,
                        const std::vector<double>&  grid_coords_z,
                        const std::vector<int64_t>& gtos_mask) -> CMatrix;

}  // namespace gtoval

#endif /* GtoValues_hpp */
