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

#include <hip/hip_runtime.h>

#include <cstdint>

#include "GtoValuesMggaGPU.hpp"

namespace gpu {  // gpu namespace

__global__ void
gtoValuesMggaRecS(double*        gto_values_0,
                  double*        gto_values_x,
                  double*        gto_values_y,
                  double*        gto_values_z,
                  double*        gto_values_xx,
                  double*        gto_values_xy,
                  double*        gto_values_xz,
                  double*        gto_values_yy,
                  double*        gto_values_yz,
                  double*        gto_values_zz,
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
        double s_0_0 = 0.0;

        double s_x_0 = 0.0;
        double s_y_0 = 0.0;
        double s_z_0 = 0.0;

        double s_xx_0 = 0.0;
        double s_xy_0 = 0.0;
        double s_xz_0 = 0.0;
        double s_yy_0 = 0.0;
        double s_yz_0 = 0.0;
        double s_zz_0 = 0.0;

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

            const auto f00 = fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto fg0 = -2.0 * fexp;
            const auto fg1 = f00 * fg0;
            const auto fg2 = fg1 * fg0;

            s_0_0 += f00;

            s_x_0 += gr_x * fg1;
            s_y_0 += gr_y * fg1;
            s_z_0 += gr_z * fg1;

            s_xx_0 += fg2 * gr_x * gr_x + fg1;
            s_xy_0 += fg2 * gr_x * gr_y;
            s_xz_0 += fg2 * gr_x * gr_z;
            s_yy_0 += fg2 * gr_y * gr_y + fg1;
            s_yz_0 += fg2 * gr_y * gr_z;
            s_zz_0 += fg2 * gr_z * gr_z + fg1;
        }

        gto_values_0[g + (i + row_offset) * ncols] = s_0_0;

        gto_values_x[g + (i + row_offset) * ncols] = s_x_0;
        gto_values_y[g + (i + row_offset) * ncols] = s_y_0;
        gto_values_z[g + (i + row_offset) * ncols] = s_z_0;

        gto_values_xx[g + (i + row_offset) * ncols] = s_xx_0;
        gto_values_xy[g + (i + row_offset) * ncols] = s_xy_0;
        gto_values_xz[g + (i + row_offset) * ncols] = s_xz_0;
        gto_values_yy[g + (i + row_offset) * ncols] = s_yy_0;
        gto_values_yz[g + (i + row_offset) * ncols] = s_yz_0;
        gto_values_zz[g + (i + row_offset) * ncols] = s_zz_0;
    }
}

__global__ void
gtoValuesMggaRecP(double*        gto_values_p3_0,
                  double*        gto_values_p3_x,
                  double*        gto_values_p3_y,
                  double*        gto_values_p3_z,
                  double*        gto_values_p3_xx,
                  double*        gto_values_p3_xy,
                  double*        gto_values_p3_xz,
                  double*        gto_values_p3_yy,
                  double*        gto_values_p3_yz,
                  double*        gto_values_p3_zz,
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
        double p_0_x = 0.0;

        double p_x_x = 0.0;
        double p_y_x = 0.0;
        double p_z_x = 0.0;

        double p_xx_x = 0.0;
        double p_xy_x = 0.0;
        double p_xz_x = 0.0;
        double p_yy_x = 0.0;
        double p_yz_x = 0.0;
        double p_zz_x = 0.0;

        double p_0_y = 0.0;

        double p_x_y = 0.0;
        double p_y_y = 0.0;
        double p_z_y = 0.0;

        double p_xx_y = 0.0;
        double p_xy_y = 0.0;
        double p_xz_y = 0.0;
        double p_yy_y = 0.0;
        double p_yz_y = 0.0;
        double p_zz_y = 0.0;

        double p_0_z = 0.0;

        double p_x_z = 0.0;
        double p_y_z = 0.0;
        double p_z_z = 0.0;

        double p_xx_z = 0.0;
        double p_xy_z = 0.0;
        double p_xz_z = 0.0;
        double p_yy_z = 0.0;
        double p_yz_z = 0.0;
        double p_zz_z = 0.0;

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

            const auto f00 = fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto fg0 = -2.0 * fexp;

            const auto fg1 = f00 * fg0;

            const auto fg2 = fg1 * fg0;

            p_0_x += f00 * gr_x;

            p_x_x += f00 * (1.0 + gr_x * gr_x * fg0);
            p_y_x += f00 * gr_x * gr_y * fg0;
            p_z_x += f00 * gr_x * gr_z * fg0;

            p_xx_x += fg2 * gr_x * gr_x * gr_x + 3.0 * fg1 * gr_x;
            p_xy_x += fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
            p_xz_x += fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
            p_yy_x += fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
            p_yz_x += fg2 * gr_x * gr_y * gr_z;
            p_zz_x += fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;

            p_0_y += f00 * gr_y;

            p_x_y += f00 * gr_y * gr_x * fg0;
            p_y_y += f00 * (1.0 + gr_y * gr_y * fg0);
            p_z_y += f00 * gr_y * gr_z * fg0;

            p_xx_y +=  + fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
            p_xy_y +=  + fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
            p_xz_y +=  + fg2 * gr_x * gr_y * gr_z;
            p_yy_y +=  + fg2 * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_y;
            p_yz_y +=  + fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
            p_zz_y +=  + fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;

            p_0_z += f00 * gr_z;

            p_x_z += f00 * gr_z * gr_x * fg0;
            p_y_z += f00 * gr_z * gr_y * fg0;
            p_z_z += f00 * (1.0 + gr_z * gr_z * fg0);

            p_xx_z +=  + fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
            p_xy_z +=  + fg2 * gr_x * gr_y * gr_z;
            p_xz_z +=  + fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;
            p_yy_z +=  + fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
            p_yz_z +=  + fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;
            p_zz_z +=  + fg2 * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_z;
        }

        // p-1: py
        // p_0: pz
        // p+1: px

        gto_values_p3_0[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_0_y;
        gto_values_p3_0[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_0_z;
        gto_values_p3_0[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_0_x;

        gto_values_p3_x[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_x_y;
        gto_values_p3_x[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_x_z;
        gto_values_p3_x[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_x_x;

        gto_values_p3_y[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_y_y;
        gto_values_p3_y[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_y_z;
        gto_values_p3_y[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_y_x;

        gto_values_p3_z[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_z_y;
        gto_values_p3_z[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_z_z;
        gto_values_p3_z[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_z_x;

        gto_values_p3_xx[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_xx_y;
        gto_values_p3_xx[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_xx_z;
        gto_values_p3_xx[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_xx_x;

        gto_values_p3_xy[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_xy_y;
        gto_values_p3_xy[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_xy_z;
        gto_values_p3_xy[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_xy_x;

        gto_values_p3_xz[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_xz_y;
        gto_values_p3_xz[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_xz_z;
        gto_values_p3_xz[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_xz_x;

        gto_values_p3_yy[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_yy_y;
        gto_values_p3_yy[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_yy_z;
        gto_values_p3_yy[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_yy_x;

        gto_values_p3_yz[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_yz_y;
        gto_values_p3_yz[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_yz_z;
        gto_values_p3_yz[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_yz_x;

        gto_values_p3_zz[g + (i + row_offset) * ncols + nrows * ncols * 0] = p_zz_y;
        gto_values_p3_zz[g + (i + row_offset) * ncols + nrows * ncols * 1] = p_zz_z;
        gto_values_p3_zz[g + (i + row_offset) * ncols + nrows * ncols * 2] = p_zz_x;
    }
}

__global__ void
gtoValuesMggaRecD(double*        gto_values_d5_0,
                  double*        gto_values_d5_x,
                  double*        gto_values_d5_y,
                  double*        gto_values_d5_z,
                  double*        gto_values_d5_xx,
                  double*        gto_values_d5_xy,
                  double*        gto_values_d5_xz,
                  double*        gto_values_d5_yy,
                  double*        gto_values_d5_yz,
                  double*        gto_values_d5_zz,
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
        double d_0_xx = 0.0;
        double d_0_xy = 0.0;
        double d_0_xz = 0.0;
        double d_0_yy = 0.0;
        double d_0_yz = 0.0;
        double d_0_zz = 0.0;

        double d_x_xx = 0.0;
        double d_x_xy = 0.0;
        double d_x_xz = 0.0;
        double d_x_yy = 0.0;
        double d_x_yz = 0.0;
        double d_x_zz = 0.0;

        double d_y_xx = 0.0;
        double d_y_xy = 0.0;
        double d_y_xz = 0.0;
        double d_y_yy = 0.0;
        double d_y_yz = 0.0;
        double d_y_zz = 0.0;

        double d_z_xx = 0.0;
        double d_z_xy = 0.0;
        double d_z_xz = 0.0;
        double d_z_yy = 0.0;
        double d_z_yz = 0.0;
        double d_z_zz = 0.0;

        double d_xx_xx = 0.0;
        double d_xx_xy = 0.0;
        double d_xx_xz = 0.0;
        double d_xx_yy = 0.0;
        double d_xx_yz = 0.0;
        double d_xx_zz = 0.0;

        double d_xy_xx = 0.0;
        double d_xy_xy = 0.0;
        double d_xy_xz = 0.0;
        double d_xy_yy = 0.0;
        double d_xy_yz = 0.0;
        double d_xy_zz = 0.0;

        double d_xz_xx = 0.0;
        double d_xz_xy = 0.0;
        double d_xz_xz = 0.0;
        double d_xz_yy = 0.0;
        double d_xz_yz = 0.0;
        double d_xz_zz = 0.0;

        double d_yy_xx = 0.0;
        double d_yy_xy = 0.0;
        double d_yy_xz = 0.0;
        double d_yy_yy = 0.0;
        double d_yy_yz = 0.0;
        double d_yy_zz = 0.0;

        double d_yz_xx = 0.0;
        double d_yz_xy = 0.0;
        double d_yz_xz = 0.0;
        double d_yz_yy = 0.0;
        double d_yz_yz = 0.0;
        double d_yz_zz = 0.0;

        double d_zz_xx = 0.0;
        double d_zz_xy = 0.0;
        double d_zz_xz = 0.0;
        double d_zz_yy = 0.0;
        double d_zz_yz = 0.0;
        double d_zz_zz = 0.0;

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

            const auto f00 = fnorm * exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto fg0 = -2.0 * fexp;

            const auto f0x = gr_x * f00;
            const auto f0y = gr_y * f00;
            const auto f0z = gr_z * f00;

            const auto fg1 = f00 * fg0;

            const auto fg2 = fg1 * fg0;

            d_0_xx += f0x * gr_x;

            d_x_xx += f0x * (2.0 + fg0 * gr_x * gr_x);
            d_y_xx += f0x * fg0 * gr_x * gr_y;
            d_z_xx += f0x * fg0 * gr_x * gr_z;

            d_xx_xx += fg2 * gr_x * gr_x * gr_x * gr_x + 5.0 * fg1 * gr_x * gr_x + 2.0 * f00;
            d_xy_xx += fg2 * gr_x * gr_x * gr_x * gr_y + 2.0 * fg1 * gr_x * gr_y;
            d_xz_xx += fg2 * gr_x * gr_x * gr_x * gr_z + 2.0 * fg1 * gr_x * gr_z;
            d_yy_xx += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x;
            d_yz_xx += fg2 * gr_x * gr_x * gr_y * gr_z;
            d_zz_xx += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x;

            d_0_xy += f0x * gr_y;

            d_x_xy += f0y * (1.0 + fg0 * gr_x * gr_x);
            d_y_xy += f0x * (1.0 + fg0 * gr_y * gr_y);
            d_z_xy += f0x * fg0 * gr_y * gr_z;

            d_xx_xy += fg2 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg1 * gr_x * gr_y;
            d_xy_xy += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x + fg1 * gr_y * gr_y + f00;
            d_xz_xy += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
            d_yy_xy += fg2 * gr_x * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y;
            d_yz_xy += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
            d_zz_xy += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;

            d_0_xz += f0x * gr_z;

            d_x_xz += f0z * (1.0 + fg0 * gr_x * gr_x);
            d_y_xz += f0x * fg0 * gr_z * gr_y;
            d_z_xz += f0x * (1.0 + fg0 * gr_z * gr_z);

            d_xx_xz += fg2 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg1 * gr_x * gr_z;
            d_xy_xz += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
            d_xz_xz += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x + fg1 * gr_z * gr_z + f00;
            d_yy_xz += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
            d_yz_xz += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;
            d_zz_xz += fg2 * gr_x * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z;

            d_0_yy += f0y * gr_y;

            d_x_yy += f0y * fg0 * gr_y * gr_x;
            d_y_yy += f0y * (2.0 + fg0 * gr_y * gr_y);
            d_z_yy += f0y * fg0 * gr_y * gr_z;

            d_xx_yy += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_y * gr_y;
            d_xy_yy += fg2 * gr_x * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_x * gr_y;
            d_xz_yy += fg2 * gr_x * gr_y * gr_y * gr_z;
            d_yy_yy += fg2 * gr_y * gr_y * gr_y * gr_y + 5.0 * fg1 * gr_y * gr_y + 2.0 * f00;
            d_yz_yy += fg2 * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
            d_zz_yy += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y;

            d_0_yz += f0y * gr_z;

            d_x_yz += f0y * fg0 * gr_z * gr_x;
            d_y_yz += f0z * (1.0 + fg0 * gr_y * gr_y);
            d_z_yz += f0y * (1.0 + fg0 * gr_z * gr_z);

            d_xx_yz += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
            d_xy_yz += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
            d_xz_yz += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;
            d_yy_yz += fg2 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_z;
            d_yz_yz += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y + fg1 * gr_z * gr_z + f00;
            d_zz_yz += fg2 * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z;

            d_0_zz += f0z * gr_z;

            d_x_zz += f0z * fg0 * gr_z * gr_x;
            d_y_zz += f0z * fg0 * gr_z * gr_y;
            d_z_zz += f0z * (2.0 + fg0 * gr_z * gr_z);

            d_xx_zz += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_z * gr_z;
            d_xy_zz += fg2 * gr_x * gr_y * gr_z * gr_z;
            d_xz_zz += fg2 * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_z;
            d_yy_zz += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_z * gr_z;
            d_yz_zz += fg2 * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_z;
            d_zz_zz += fg2 * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_z * gr_z + 2.0 * f00;
        }

        // d-2: dxy * f2_3
        // d-1: dyz * f2_3
        // d_0: dzz * 2.0 - dxx - dyy
        // d+1: dxz * f2_3
        // d+2: (dxx - dyy) * 0.5 * f2_3

        gto_values_d5_0[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_0_xy * f2_3;
        gto_values_d5_0[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_0_yz * f2_3;
        gto_values_d5_0[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_0_zz * 2.0 - d_0_xx - d_0_yy);
        gto_values_d5_0[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_0_xz * f2_3;
        gto_values_d5_0[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_0_xx - d_0_yy) * 0.5 * f2_3;

        gto_values_d5_x[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_x_xy * f2_3;
        gto_values_d5_x[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_x_yz * f2_3;
        gto_values_d5_x[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_x_zz * 2.0 - d_x_xx - d_x_yy);
        gto_values_d5_x[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_x_xz * f2_3;
        gto_values_d5_x[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_x_xx - d_x_yy) * 0.5 * f2_3;

        gto_values_d5_y[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_y_xy * f2_3;
        gto_values_d5_y[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_y_yz * f2_3;
        gto_values_d5_y[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_y_zz * 2.0 - d_y_xx - d_y_yy);
        gto_values_d5_y[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_y_xz * f2_3;
        gto_values_d5_y[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_y_xx - d_y_yy) * 0.5 * f2_3;

        gto_values_d5_z[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_z_xy * f2_3;
        gto_values_d5_z[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_z_yz * f2_3;
        gto_values_d5_z[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_z_zz * 2.0 - d_z_xx - d_z_yy);
        gto_values_d5_z[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_z_xz * f2_3;
        gto_values_d5_z[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_z_xx - d_z_yy) * 0.5 * f2_3;

        gto_values_d5_xx[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_xx_xy * f2_3;
        gto_values_d5_xx[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_xx_yz * f2_3;
        gto_values_d5_xx[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_xx_zz * 2.0 - d_xx_xx - d_xx_yy);
        gto_values_d5_xx[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_xx_xz * f2_3;
        gto_values_d5_xx[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_xx_xx - d_xx_yy) * 0.5 * f2_3;

        gto_values_d5_xy[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_xy_xy * f2_3;
        gto_values_d5_xy[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_xy_yz * f2_3;
        gto_values_d5_xy[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_xy_zz * 2.0 - d_xy_xx - d_xy_yy);
        gto_values_d5_xy[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_xy_xz * f2_3;
        gto_values_d5_xy[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_xy_xx - d_xy_yy) * 0.5 * f2_3;

        gto_values_d5_xz[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_xz_xy * f2_3;
        gto_values_d5_xz[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_xz_yz * f2_3;
        gto_values_d5_xz[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_xz_zz * 2.0 - d_xz_xx - d_xz_yy);
        gto_values_d5_xz[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_xz_xz * f2_3;
        gto_values_d5_xz[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_xz_xx - d_xz_yy) * 0.5 * f2_3;

        gto_values_d5_yy[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_yy_xy * f2_3;
        gto_values_d5_yy[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_yy_yz * f2_3;
        gto_values_d5_yy[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_yy_zz * 2.0 - d_yy_xx - d_yy_yy);
        gto_values_d5_yy[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_yy_xz * f2_3;
        gto_values_d5_yy[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_yy_xx - d_yy_yy) * 0.5 * f2_3;

        gto_values_d5_yz[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_yz_xy * f2_3;
        gto_values_d5_yz[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_yz_yz * f2_3;
        gto_values_d5_yz[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_yz_zz * 2.0 - d_yz_xx - d_yz_yy);
        gto_values_d5_yz[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_yz_xz * f2_3;
        gto_values_d5_yz[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_yz_xx - d_yz_yy) * 0.5 * f2_3;

        gto_values_d5_zz[g + (i + row_offset) * ncols + nrows * ncols * 0] = d_zz_xy * f2_3;
        gto_values_d5_zz[g + (i + row_offset) * ncols + nrows * ncols * 1] = d_zz_yz * f2_3;
        gto_values_d5_zz[g + (i + row_offset) * ncols + nrows * ncols * 2] = (d_zz_zz * 2.0 - d_zz_xx - d_zz_yy);
        gto_values_d5_zz[g + (i + row_offset) * ncols + nrows * ncols * 3] = d_zz_xz * f2_3;
        gto_values_d5_zz[g + (i + row_offset) * ncols + nrows * ncols * 4] = (d_zz_xx - d_zz_yy) * 0.5 * f2_3;
    }
}

}  // namespace gpu
