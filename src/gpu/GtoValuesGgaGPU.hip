//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include <hip/hip_runtime.h>

#include <cstdint>

#include "GtoValuesGgaGPU.hpp"

namespace gpu {  // gpu namespace

__global__ void
gtoValuesGgaRecS(double*        gto_values_0,
                 double*        gto_values_x,
                 double*        gto_values_y,
                 double*        gto_values_z,
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

            const auto g0 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto g1 = -2.0 * fexp * g0;

            s_0_0 += g0;
            s_x_0 += gr_x * g1;
            s_y_0 += gr_y * g1;
            s_z_0 += gr_z * g1;
        }

        gto_values_0[g + (i + row_offset) * ncols] = s_0_0;
        gto_values_x[g + (i + row_offset) * ncols] = s_x_0;
        gto_values_y[g + (i + row_offset) * ncols] = s_y_0;
        gto_values_z[g + (i + row_offset) * ncols] = s_z_0;
    }
}

__global__ void
gtoValuesGgaRecP(double*        gto_values_p3_0,
                 double*        gto_values_p3_x,
                 double*        gto_values_p3_y,
                 double*        gto_values_p3_z,
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
        double p_0_y = 0.0;
        double p_0_z = 0.0;

        double p_x_x = 0.0;
        double p_x_y = 0.0;
        double p_x_z = 0.0;

        double p_y_x = 0.0;
        double p_y_y = 0.0;
        double p_y_z = 0.0;

        double p_z_x = 0.0;
        double p_z_y = 0.0;
        double p_z_z = 0.0;

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

            const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto fg0 = -2.0 * fexp;

            p_0_x += f00 * gr_x;
            p_x_x += f00 * (1.0 + gr_x * gr_x * fg0);
            p_y_x += f00 * gr_x * gr_y * fg0;
            p_z_x += f00 * gr_x * gr_z * fg0;

            p_0_y += f00 * gr_y;
            p_x_y += f00 * gr_y * gr_x * fg0;
            p_y_y += f00 * (1.0 + gr_y * gr_y * fg0);
            p_z_y += f00 * gr_y * gr_z * fg0;

            p_0_z += f00 * gr_z;
            p_x_z += f00 * gr_z * gr_x * fg0;
            p_y_z += f00 * gr_z * gr_y * fg0;
            p_z_z += f00 * (1.0 + gr_z * gr_z * fg0);
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
    }
}

__global__ void
gtoValuesGgaRecD(double*        gto_values_d5_0,
                 double*        gto_values_d5_x,
                 double*        gto_values_d5_y,
                 double*        gto_values_d5_z,
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

            const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

            const auto fg0 = -2.0 * fexp;

            const auto f0x = gr_x * f00;
            const auto f0y = gr_y * f00;
            const auto f0z = gr_z * f00;

            d_0_xx += f0x * gr_x;
            d_x_xx += f0x * (2.0 + fg0 * gr_x * gr_x);
            d_y_xx += f0x * fg0 * gr_x * gr_y;
            d_z_xx += f0x * fg0 * gr_x * gr_z;

            d_0_xy += f0x * gr_y;
            d_x_xy += f0y * (1.0 + fg0 * gr_x * gr_x);
            d_y_xy += f0x * (1.0 + fg0 * gr_y * gr_y);
            d_z_xy += f0x * fg0 * gr_y * gr_z;

            d_0_xz += f0x * gr_z;
            d_x_xz += f0z * (1.0 + fg0 * gr_x * gr_x);
            d_y_xz += f0x * fg0 * gr_z * gr_y;
            d_z_xz += f0x * (1.0 + fg0 * gr_z * gr_z);

            d_0_yy += f0y * gr_y;
            d_x_yy += f0y * fg0 * gr_y * gr_x;
            d_y_yy += f0y * (2.0 + fg0 * gr_y * gr_y);
            d_z_yy += f0y * fg0 * gr_y * gr_z;

            d_0_yz += f0y * gr_z;
            d_x_yz += f0y * fg0 * gr_z * gr_x;
            d_y_yz += f0z * (1.0 + fg0 * gr_y * gr_y);
            d_z_yz += f0y * (1.0 + fg0 * gr_z * gr_z);

            d_0_zz += f0z * gr_z;
            d_x_zz += f0z * fg0 * gr_z * gr_x;
            d_y_zz += f0z * fg0 * gr_z * gr_y;
            d_z_zz += f0z * (2.0 + fg0 * gr_z * gr_z);
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
    }
}

}  // namespace gpu
