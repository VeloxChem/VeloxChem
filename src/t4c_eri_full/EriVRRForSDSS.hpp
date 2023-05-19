//                                                                              
//                           VELOXCHEM 1.0-RC2                                  
//         ----------------------------------------------------                 
//                     An Electronic Structure Code                             
//                                                                              
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.         
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

#include <cstdint>

#include "Buffer.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostVRRForSDSS_V0(      BufferHostXY<T>&      intsBufferSDSS,
                      const BufferHostX<int32_t>& intsIndexesSDSS0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSPSS0,
                      const BufferHostX<int32_t>& intsIndexesSPSS0,
                      const BufferHostXY<T>&      intsBufferSPSS1,
                      const BufferHostX<int32_t>& intsIndexesSPSS1,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fz_0 = osFactorsBraZeta;

    auto frz2_0 = osFactorsBraRhoZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SDSS]^(0) integral components

    t_0_zz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(0));

    t_0_yz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(1));

    t_0_yy_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(2));

    t_0_xz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(3));

    t_0_xy_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(4));

    t_0_xx_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(5));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SPSS]^(0) integral components

    t_0_z_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(0));

    t_0_y_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(1));

    t_0_x_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(2));

    // set up [SPSS]^(1) integral components

    t_0_z_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(0));

    t_0_y_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(1));

    t_0_x_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_xy_0_0_0, t_0_xz_0_0_0, t_0_y_0_0_0, t_0_y_0_0_1,\
                               t_0_yy_0_0_0, t_0_yz_0_0_0, t_0_z_0_0_0, t_0_z_0_0_1,\
                               t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] += rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yz_0_0_0[i] += rpb_y[i] * t_0_z_0_0_0[i] + rwp_y[i] * t_0_z_0_0_1[i];

            t_0_yy_0_0_0[i] += rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xz_0_0_0[i] += rpb_x[i] * t_0_z_0_0_0[i] + rwp_x[i] * t_0_z_0_0_1[i];

            t_0_xy_0_0_0[i] += rpb_x[i] * t_0_y_0_0_0[i] + rwp_x[i] * t_0_y_0_0_1[i];

            t_0_xx_0_0_0[i] += rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_xy_0_0_0, t_0_xz_0_0_0, t_0_y_0_0_0, t_0_y_0_0_1,\
                               t_0_yy_0_0_0, t_0_yz_0_0_0, t_0_z_0_0_0, t_0_z_0_0_1,\
                               t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] = rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yz_0_0_0[i] = rpb_y[i] * t_0_z_0_0_0[i] + rwp_y[i] * t_0_z_0_0_1[i];

            t_0_yy_0_0_0[i] = rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xz_0_0_0[i] = rpb_x[i] * t_0_z_0_0_0[i] + rwp_x[i] * t_0_z_0_0_1[i];

            t_0_xy_0_0_0[i] = rpb_x[i] * t_0_y_0_0_0[i] + rwp_x[i] * t_0_y_0_0_1[i];

            t_0_xx_0_0_0[i] = rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSDSS_V1(      BufferHostXY<T>&      intsBufferSDSS,
                      const BufferHostX<int32_t>& intsIndexesSDSS0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSPSS0,
                      const BufferHostX<int32_t>& intsIndexesSPSS0,
                      const BufferHostXY<T>&      intsBufferSPSS1,
                      const BufferHostX<int32_t>& intsIndexesSPSS1,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fz_0 = osFactorsBraZeta;

    auto frz2_0 = osFactorsBraRhoZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SDSS]^(0) integral components

    t_0_zz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(0));

    t_0_yz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(1));

    t_0_yy_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(2));

    t_0_xx_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(3));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SPSS]^(0) integral components

    t_0_z_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(0));

    t_0_y_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(1));

    t_0_x_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(2));

    // set up [SPSS]^(1) integral components

    t_0_z_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(0));

    t_0_y_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(1));

    t_0_x_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_y_0_0_0, t_0_y_0_0_1, t_0_yy_0_0_0, t_0_yz_0_0_0,\
                               t_0_z_0_0_0, t_0_z_0_0_1, t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] += rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yz_0_0_0[i] += rpb_y[i] * t_0_z_0_0_0[i] + rwp_y[i] * t_0_z_0_0_1[i];

            t_0_yy_0_0_0[i] += rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xx_0_0_0[i] += rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_y_0_0_0, t_0_y_0_0_1, t_0_yy_0_0_0, t_0_yz_0_0_0,\
                               t_0_z_0_0_0, t_0_z_0_0_1, t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] = rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yz_0_0_0[i] = rpb_y[i] * t_0_z_0_0_0[i] + rwp_y[i] * t_0_z_0_0_1[i];

            t_0_yy_0_0_0[i] = rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xx_0_0_0[i] = rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSDSS_V2(      BufferHostXY<T>&      intsBufferSDSS,
                      const BufferHostX<int32_t>& intsIndexesSDSS0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSPSS0,
                      const BufferHostX<int32_t>& intsIndexesSPSS0,
                      const BufferHostXY<T>&      intsBufferSPSS1,
                      const BufferHostX<int32_t>& intsIndexesSPSS1,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fz_0 = osFactorsBraZeta;

    auto frz2_0 = osFactorsBraRhoZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SDSS]^(0) integral components

    t_0_zz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(0));

    t_0_yy_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(1));

    t_0_xx_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(2));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SPSS]^(0) integral components

    t_0_z_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(0));

    t_0_y_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(1));

    t_0_x_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(2));

    // set up [SPSS]^(1) integral components

    t_0_z_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(0));

    t_0_y_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(1));

    t_0_x_0_0_1 = intsBufferSPSS1.data(intsIndexesSPSS1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_y_0_0_0, t_0_y_0_0_1, t_0_yy_0_0_0, t_0_z_0_0_0,\
                               t_0_z_0_0_1, t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] += rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yy_0_0_0[i] += rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xx_0_0_0[i] += rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_x_0_0_0, t_0_x_0_0_1, t_0_xx_0_0_0,\
                               t_0_y_0_0_0, t_0_y_0_0_1, t_0_yy_0_0_0, t_0_z_0_0_0,\
                               t_0_z_0_0_1, t_0_zz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_0_0[i] = rpb_z[i] * t_0_z_0_0_0[i] + rwp_z[i] * t_0_z_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_yy_0_0_0[i] = rpb_y[i] * t_0_y_0_0_0[i] + rwp_y[i] * t_0_y_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];

            t_0_xx_0_0_0[i] = rpb_x[i] * t_0_x_0_0_0[i] + rwp_x[i] * t_0_x_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_0_1[i];
        }
    }
}


} // derirec namespace
