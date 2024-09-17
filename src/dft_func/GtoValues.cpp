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

#include "GtoValues.hpp"

#include "GtoValuesRecD.hpp"
#include "GtoValuesRecF.hpp"
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
        return gtoval::getLdaValuesRecS(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::getLdaValuesRecP(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::getLdaValuesRecD(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::getLdaValuesRecF(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

auto
get_gto_values_for_gga(const CGtoBlock&            gto_block,
                   const std::vector<double>&  grid_coords_x,
                   const std::vector<double>&  grid_coords_y,
                   const std::vector<double>&  grid_coords_z,
                   const std::vector<int>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.angular_momentum();

    if (gto_ang == 0)
    {
        return gtoval::getGgaValuesRecS(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::getGgaValuesRecP(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::getGgaValuesRecD(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    return CMatrix();
}

}  // namespace gtoval
