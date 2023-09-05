#ifndef GtoValuesRecP_hpp
#define GtoValuesRecP_hpp

#include <cstdint>
#include <vector>

#include "Matrix.hpp"
#include "GtoBlock.hpp"

namespace dft {  // dft namespace

/**
 Computes GTO values on given grid for P type GTOs.

 @param gto_block the GTOs block.
 @param grid_coords_x the vector of Cartesian X coordinates of grid.
 @param grid_coords_y the vector of Cartesian Y coordinates of grid.
 @param grid_coords_z the vector of Cartesian Z coordinates of grid.
 @param gtos_mask the mask for GTOs (1 evaluate, 0 skip).
 */
auto getLDAValuesRecP(const CGtoBlock&            gto_block,
                      const std::vector<double>&  grid_coords_x,
                      const std::vector<double>&  grid_coords_y,
                      const std::vector<double>&  grid_coords_z,
                      const std::vector<int64_t>& gtos_mask) -> CMatrix;

}  // namespace dft

#endif /* GtoValuesRecP_hpp */
