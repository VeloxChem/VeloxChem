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

#include "GpuRuntime.hpp"


#include <cstdint>

#include "GtoValuesLdaGPU.hpp"

namespace gpu {  // gpu namespace

__global__ void
gtoValuesLdaRecS(double*        gto_values,
                 const uint32_t row_offset,
                 const double*  gto_info,
                 const double*  grid_x,
                 const double*  grid_y,
                 const double*  grid_z,
                 const uint32_t grid_offset,
                 const uint32_t nrows,
                 const uint32_t npgtos,
                 const uint32_t ncols)
{
    const uint32_t g = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t i = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        double s0 = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            s0 += fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
        }

        gto_values[g + (i + row_offset) * ncols] = s0;
    }
}

__global__ void
gtoValuesLdaRecP(double*        gto_values_p3,
                 const uint32_t row_offset,
                 const double*  gto_info,
                 const double*  grid_x,
                 const double*  grid_y,
                 const double*  grid_z,
                 const uint32_t grid_offset,
                 const uint32_t nrows,
                 const uint32_t npgtos,
                 const uint32_t ncols)
{
    const uint32_t g = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t i = blockDim.y * blockIdx.y + threadIdx.y;

    if ((i < nrows) && (g < ncols))
    {
        double px = 0.0;
        double py = 0.0;
        double pz = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            const auto fss = fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            px += gr_x * fss;
            py += gr_y * fss;
            pz += gr_z * fss;
        }

        // p-1: py
        // p_0: pz
        // p+1: px

        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 0] = py;
        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 1] = pz;
        gto_values_p3[g + (i + row_offset) * ncols + nrows * ncols * 2] = px;
    }
}

__global__ void
gtoValuesLdaRecD(double*        gto_values_d5,
                 const uint32_t row_offset,
                 const double   f2_3,
                 const double*  gto_info,
                 const double*  grid_x,
                 const double*  grid_y,
                 const double*  grid_z,
                 const uint32_t grid_offset,
                 const uint32_t nrows,
                 const uint32_t npgtos,
                 const uint32_t ncols)
{
    const uint32_t i = blockDim.y * blockIdx.y + threadIdx.y;
    const uint32_t g = blockDim.x * blockIdx.x + threadIdx.x;

    if ((i < nrows) && (g < ncols))
    {
        double dxx = 0.0;
        double dxy = 0.0;
        double dxz = 0.0;
        double dyy = 0.0;
        double dyz = 0.0;
        double dzz = 0.0;

        const auto g_x = grid_x[g + grid_offset];
        const auto g_y = grid_y[g + grid_offset];
        const auto g_z = grid_z[g + grid_offset];

        for (uint32_t j = 0; j < npgtos; j++)
        {
            const auto fexp  = gto_info[i + j * nrows + npgtos * nrows * 0];
            const auto fnorm = gto_info[i + j * nrows + npgtos * nrows * 1];
            const auto r_x   = gto_info[i + j * nrows + npgtos * nrows * 2];
            const auto r_y   = gto_info[i + j * nrows + npgtos * nrows * 3];
            const auto r_z   = gto_info[i + j * nrows + npgtos * nrows * 4];

            const auto gr_x = g_x - r_x;
            const auto gr_y = g_y - r_y;
            const auto gr_z = g_z - r_z;

            const auto fss = fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            dxx += gr_x * gr_x * fss;
            dxy += gr_x * gr_y * fss;
            dxz += gr_x * gr_z * fss;
            dyy += gr_y * gr_y * fss;
            dyz += gr_y * gr_z * fss;
            dzz += gr_z * gr_z * fss;
        }

        // d-2: dxy * f2_3
        // d-1: dyz * f2_3
        // d_0: dzz * 2.0 - dxx - dyy
        // d+1: dxz * f2_3
        // d+2: (dxx - dyy) * 0.5 * f2_3

        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 0] = dxy * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 1] = dyz * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 2] = (dzz * 2.0 - dxx - dyy);
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 3] = dxz * f2_3;
        gto_values_d5[g + (i + row_offset) * ncols + nrows * ncols * 4] = (dxx - dyy) * 0.5 * f2_3;
    }
}

}  // namespace gpu
