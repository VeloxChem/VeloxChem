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

#include "GtoValues.hpp"

#include "ErrorHandler.hpp"
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
    else if (gto_ang == 3)
    {
        return gtoval::get_lda_values_rec_f(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    std::string errangmom("get_gto_values_for_lda: Only implemented up to f-orbitals");

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
    else if (gto_ang == 3)
    {
        return gtoval::get_gga_values_rec_f(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    std::string errangmom("get_gto_values_for_gga: Only implemented up to f-orbitals");

    errors::assertMsgCritical(false, errangmom);

    return CMatrix();
}

auto
get_gto_values_for_mgga(const CGtoBlock&            gto_block,
                   const std::vector<double>&  grid_coords_x,
                   const std::vector<double>&  grid_coords_y,
                   const std::vector<double>&  grid_coords_z,
                   const std::vector<int>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.angular_momentum();

    if (gto_ang == 0)
    {
        return gtoval::get_mgga_values_rec_s(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::get_mgga_values_rec_p(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::get_mgga_values_rec_d(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::get_mgga_values_rec_f(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    std::string errangmom("get_gto_values_for_mgga: Only implemented up to f-orbitals");

    errors::assertMsgCritical(false, errangmom);

    return CMatrix();
}

auto
get_gto_values_for_3rd_order(const CGtoBlock&            gto_block,
                             const std::vector<double>&  grid_coords_x,
                             const std::vector<double>&  grid_coords_y,
                             const std::vector<double>&  grid_coords_z,
                             const std::vector<int>& gtos_mask) -> CMatrix
{
    auto gto_ang = gto_block.angular_momentum();

    if (gto_ang == 0)
    {
        return gtoval::get_3rd_order_values_rec_s(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 1)
    {
        return gtoval::get_3rd_order_values_rec_p(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 2)
    {
        return gtoval::get_3rd_order_values_rec_d(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }
    else if (gto_ang == 3)
    {
        return gtoval::get_3rd_order_values_rec_f(gto_block, grid_coords_x, grid_coords_y, grid_coords_z, gtos_mask);
    }

    std::string errangmom("get_gto_values_for_3rd_order: Only implemented up to f-orbitals");

    errors::assertMsgCritical(false, errangmom);

    return CMatrix();
}

}  // namespace gtoval
