#include "GtoValues.hpp"

#include "GtoValuesRecD.hpp"
#include "GtoValuesRecF.hpp"
#include "GtoValuesRecP.hpp"
#include "GtoValuesRecS.hpp"
#include "Matrix.hpp"

namespace gtoval {  // gtoval namespace

auto
getGtoValuesForLda(const CGtoBlock&            gto_block,
                   const int64_t               n_points,
                   const double*               grid_coords_x,
                   const double*               grid_coords_y,
                   const double*               grid_coords_z,
                   const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.getAngularMomentum();

    if (gto_ang == 0)
    {
        return gtoval::getLdaValuesRecS(gto_block, n_points, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::getLdaValuesRecP(gto_block, n_points, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::getLdaValuesRecD(gto_block, n_points, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::getLdaValuesRecF(gto_block, n_points, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

}  // namespace gtoval
