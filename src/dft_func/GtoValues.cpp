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

#include "ErrorHandler.hpp"
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

    std::string errangmom("get_gto_values_for_lda: Only implemented up to d-orbitals");

    errors::assertMsgCritical(false, errangmom);

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

    std::string errangmom("get_gto_values_for_gga: Only implemented up to d-orbitals");

    errors::assertMsgCritical(false, errangmom);

    return CMatrix();
}

}  // namespace gtoval
