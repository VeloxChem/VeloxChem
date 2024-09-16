#include "GtoValues.hpp"

#include "GtoValuesRecD.hpp"
#include "GtoValuesRecP.hpp"
#include "GtoValuesRecS.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

auto
get_gto_values_for_lda(const CGtoBlock&            gto_block,
                       const std::vector<double>&  grid_coords_x,
                       const std::vector<double>&  grid_coords_y,
                       const std::vector<double>&  grid_coords_z,
                       const std::vector<int>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.angular_momentum();

    if (gto_ang == 0)
    {
        return gtoval::get_lda_values_rec_s(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::get_lda_values_rec_p(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::get_lda_values_rec_d(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

auto
get_gto_values_for_gga(const CGtoBlock&            gto_block,
                       const std::vector<double>&  grid_coords_x,
                       const std::vector<double>&  grid_coords_y,
                       const std::vector<double>&  grid_coords_z,
                       const std::vector<int>&    gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.angular_momentum();

    if (gto_ang == 0)
    {
        return gtoval::get_gga_values_rec_s(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::get_gga_values_rec_p(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::get_gga_values_rec_d(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

}  // namespace gtoval
