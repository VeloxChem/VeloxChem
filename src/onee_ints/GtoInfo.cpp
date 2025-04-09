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

#include "GtoInfo.hpp"

#include "MathFunc.hpp"

namespace gtoinfo {  // gtoinfo namespace

auto
updatePrimitiveInfoForS(double* s_prim_info, int* s_prim_aoinds, const int s_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int s_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 0)
        {
            const auto ncgtos = gto_block.number_of_basis_functions();
            const auto npgtos = gto_block.number_of_primitives();

            const auto gto_exps   = gto_block.exponents();
            const auto gto_norms  = gto_block.normalization_factors();
            const auto gto_coords = gto_block.coordinates();

            const auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

            for (int j = 0; j < npgtos; j++)
            {
                for (int i = 0; i < ncgtos; i++)
                {
                    auto gto_coords_i = gto_coords[i].coordinates();

                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = gto_exps[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 1] = gto_norms[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 2] = gto_coords_i[0];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 3] = gto_coords_i[1];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 4] = gto_coords_i[2];

                    s_prim_aoinds[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = gto_ao_inds[i];
                }
            }

            s_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
updatePrimitiveInfoForP(double* p_prim_info, int* p_prim_aoinds, const int p_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int p_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 1)
        {
            const auto ncgtos = gto_block.number_of_basis_functions();
            const auto npgtos = gto_block.number_of_primitives();

            const auto gto_exps     = gto_block.exponents();
            const auto gto_norms    = gto_block.normalization_factors();
            const auto gto_coords = gto_block.coordinates();

            const auto gto_ao_inds_cart = gto_block.getAtomicOrbitalsIndexesForCartesian();

            for (int j = 0; j < npgtos; j++)
            {
                for (int i = 0; i < ncgtos; i++)
                {
                    auto gto_coords_i = gto_coords[i].coordinates();

                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 0] = gto_exps[i + j * ncgtos];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 1] = gto_norms[i + j * ncgtos];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 2] = gto_coords_i[0];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 3] = gto_coords_i[1];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 4] = gto_coords_i[2];

                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 0] = gto_ao_inds_cart[i + ncgtos * 0];
                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 1] = gto_ao_inds_cart[i + ncgtos * 1];
                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 2] = gto_ao_inds_cart[i + ncgtos * 2];
                }
            }

            p_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
updatePrimitiveInfoForD(double* d_prim_info, int* d_prim_aoinds, const int d_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int d_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 2)
        {
            const auto ncgtos = gto_block.number_of_basis_functions();
            const auto npgtos = gto_block.number_of_primitives();

            const auto gto_exps     = gto_block.exponents();
            const auto gto_norms    = gto_block.normalization_factors();
            const auto gto_coords = gto_block.coordinates();

            const auto gto_ao_inds_cart = gto_block.getAtomicOrbitalsIndexesForCartesian();

            for (int j = 0; j < npgtos; j++)
            {
                for (int i = 0; i < ncgtos; i++)
                {
                    auto gto_coords_i = gto_coords[i].coordinates();

                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 0] = gto_exps[i + j * ncgtos];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 1] = gto_norms[i + j * ncgtos];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 2] = gto_coords_i[0];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 3] = gto_coords_i[1];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 4] = gto_coords_i[2];

                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 0] = gto_ao_inds_cart[i + ncgtos * 0];
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 1] = gto_ao_inds_cart[i + ncgtos * 1];
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 2] = gto_ao_inds_cart[i + ncgtos * 2];
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 3] = gto_ao_inds_cart[i + ncgtos * 3];
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 4] = gto_ao_inds_cart[i + ncgtos * 4];
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 5] = gto_ao_inds_cart[i + ncgtos * 5];
                }
            }

            d_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
updatePrimitiveInfoForF(double* f_prim_info, int* f_prim_aoinds, const int f_prim_count, const std::vector<CGtoBlock>& gto_blocks, const int ncgtos_d) -> void
{
    int f_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 3)
        {
            const auto ncgtos = gto_block.number_of_basis_functions();
            const auto npgtos = gto_block.number_of_primitives();

            const auto gto_exps     = gto_block.exponents();
            const auto gto_norms    = gto_block.normalization_factors();
            const auto gto_coords = gto_block.coordinates();

            const auto gto_ao_inds_cart = gto_block.getAtomicOrbitalsIndexesForCartesian(ncgtos_d);

            for (int j = 0; j < npgtos; j++)
            {
                for (int i = 0; i < ncgtos; i++)
                {
                    auto gto_coords_i = gto_coords[i].coordinates();

                    f_prim_info[f_prim_offset + i + j * ncgtos + f_prim_count * 0] = gto_exps[i + j * ncgtos];
                    f_prim_info[f_prim_offset + i + j * ncgtos + f_prim_count * 1] = gto_norms[i + j * ncgtos];
                    f_prim_info[f_prim_offset + i + j * ncgtos + f_prim_count * 2] = gto_coords_i[0];
                    f_prim_info[f_prim_offset + i + j * ncgtos + f_prim_count * 3] = gto_coords_i[1];
                    f_prim_info[f_prim_offset + i + j * ncgtos + f_prim_count * 4] = gto_coords_i[2];

                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 0] = gto_ao_inds_cart[i + ncgtos * 0];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 1] = gto_ao_inds_cart[i + ncgtos * 1];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 2] = gto_ao_inds_cart[i + ncgtos * 2];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 3] = gto_ao_inds_cart[i + ncgtos * 3];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 4] = gto_ao_inds_cart[i + ncgtos * 4];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 5] = gto_ao_inds_cart[i + ncgtos * 5];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 6] = gto_ao_inds_cart[i + ncgtos * 6];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 7] = gto_ao_inds_cart[i + ncgtos * 7];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 8] = gto_ao_inds_cart[i + ncgtos * 8];
                    f_prim_aoinds[f_prim_offset + i + j * ncgtos + f_prim_count * 9] = gto_ao_inds_cart[i + ncgtos * 9];
                }
            }

            f_prim_offset += npgtos * ncgtos;
        }
    }
}

}  // namespace gtoinfo
