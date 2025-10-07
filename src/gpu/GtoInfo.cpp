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
updatePrimitiveInfoForS(double* s_prim_info, int32_t* s_prim_aoinds, const int64_t s_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int64_t s_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0)
        {
            const auto ncgtos = gto_block.getNumberOfBasisFunctions();
            const auto npgtos = gto_block.getNumberOfPrimitives();

            const auto gto_exps   = gto_block.getExponents();
            const auto gto_norms  = gto_block.getNormalizationFactors();
            const auto gto_coords = gto_block.getCoordinates();

            const auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

            for (int64_t j = 0; j < npgtos; j++)
            {
                for (int64_t i = 0; i < ncgtos; i++)
                {
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = gto_exps[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 1] = gto_norms[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 2] = gto_coords[i][0];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 3] = gto_coords[i][1];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 4] = gto_coords[i][2];

                    s_prim_aoinds[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = static_cast<int32_t>(gto_ao_inds[i]);
                }
            }

            s_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
updatePrimitiveInfoForP(double* p_prim_info, int32_t* p_prim_aoinds, const int64_t p_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int64_t p_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 1)
        {
            const auto ncgtos = gto_block.getNumberOfBasisFunctions();
            const auto npgtos = gto_block.getNumberOfPrimitives();

            const auto gto_exps   = gto_block.getExponents();
            const auto gto_norms  = gto_block.getNormalizationFactors();
            const auto gto_coords = gto_block.getCoordinates();

            const auto gto_ao_inds_cart = gto_block.getAtomicOrbitalsIndexesForCartesian();

            for (int64_t j = 0; j < npgtos; j++)
            {
                for (int64_t i = 0; i < ncgtos; i++)
                {
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 0] = gto_exps[i + j * ncgtos];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 1] = gto_norms[i + j * ncgtos];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 2] = gto_coords[i][0];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 3] = gto_coords[i][1];
                    p_prim_info[p_prim_offset + i + j * ncgtos + p_prim_count * 4] = gto_coords[i][2];

                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 0] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 0]);
                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 1] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 1]);
                    p_prim_aoinds[p_prim_offset + i + j * ncgtos + p_prim_count * 2] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 2]);
                }
            }

            p_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
updatePrimitiveInfoForD(double* d_prim_info, int32_t* d_prim_aoinds, const int64_t d_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void
{
    int64_t d_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 2)
        {
            const auto ncgtos = gto_block.getNumberOfBasisFunctions();
            const auto npgtos = gto_block.getNumberOfPrimitives();

            const auto gto_exps   = gto_block.getExponents();
            const auto gto_norms  = gto_block.getNormalizationFactors();
            const auto gto_coords = gto_block.getCoordinates();

            const auto gto_ao_inds_cart = gto_block.getAtomicOrbitalsIndexesForCartesian();

            for (int64_t j = 0; j < npgtos; j++)
            {
                for (int64_t i = 0; i < ncgtos; i++)
                {
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 0] = gto_exps[i + j * ncgtos];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 1] = gto_norms[i + j * ncgtos];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 2] = gto_coords[i][0];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 3] = gto_coords[i][1];
                    d_prim_info[d_prim_offset + i + j * ncgtos + d_prim_count * 4] = gto_coords[i][2];

                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 0] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 0]);
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 1] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 1]);
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 2] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 2]);
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 3] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 3]);
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 4] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 4]);
                    d_prim_aoinds[d_prim_offset + i + j * ncgtos + d_prim_count * 5] = static_cast<int32_t>(gto_ao_inds_cart[i + ncgtos * 5]);
                }
            }

            d_prim_offset += npgtos * ncgtos;
        }
    }
}

auto
getGtoInfo(const CGtoBlock gto_block) -> std::vector<double>
{
    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    std::vector<double> gto_info(5 * ncgtos * npgtos);

    auto gto_info_ptr = gto_info.data();

    for (int64_t i = 0; i < ncgtos; i++)
    {
        const auto r_x = gto_coords[i][0];
        const auto r_y = gto_coords[i][1];
        const auto r_z = gto_coords[i][2];

        for (int64_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_exps[j * ncgtos + i];
            const auto fnorm = gto_norms[j * ncgtos + i];

            gto_info_ptr[i + j * ncgtos + npgtos * ncgtos * 0] = fexp;
            gto_info_ptr[i + j * ncgtos + npgtos * ncgtos * 1] = fnorm;
            gto_info_ptr[i + j * ncgtos + npgtos * ncgtos * 2] = r_x;
            gto_info_ptr[i + j * ncgtos + npgtos * ncgtos * 3] = r_y;
            gto_info_ptr[i + j * ncgtos + npgtos * ncgtos * 4] = r_z;
        }
    }

    return gto_info;
}

auto
getGtoInfo(const CGtoBlock gto_block, const std::vector<int64_t>& gtos_mask) -> std::vector<double>
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // set up data on host and device

    std::vector<double> gto_info(5 * nrows * npgtos);

    auto gto_info_ptr = gto_info.data();

    for (int64_t i = 0, irow = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            const auto r_x = gto_coords[i][0];
            const auto r_y = gto_coords[i][1];
            const auto r_z = gto_coords[i][2];

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto fexp  = gto_exps[j * ncgtos + i];
                const auto fnorm = gto_norms[j * ncgtos + i];

                gto_info_ptr[irow + j * nrows + npgtos * nrows * 0] = fexp;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 1] = fnorm;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 2] = r_x;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 3] = r_y;
                gto_info_ptr[irow + j * nrows + npgtos * nrows * 4] = r_z;
            }

            irow++;
        }
    }

    return gto_info;
}

}  // namespace gtoinfo
