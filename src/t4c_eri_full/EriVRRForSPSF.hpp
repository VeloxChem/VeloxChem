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
compHostVRRForSPSF_V0(      BufferHostXY<T>&      intsBufferSPSF,
                      const BufferHostX<int32_t>& intsIndexesSPSF0,
                      const BufferHostXY<T>&      intsBufferSSSD1,
                      const BufferHostX<int32_t>& intsIndexesSSSD1,
                      const BufferHostXY<T>&      intsBufferSSSF0,
                      const BufferHostX<int32_t>& intsIndexesSSSF0,
                      const BufferHostXY<T>&      intsBufferSSSF1,
                      const BufferHostX<int32_t>& intsIndexesSSSF1,
                      const T*                    osFactorsZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SPSF]^(0) integral components

    t_0_z_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(0));

    t_0_z_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(1));

    t_0_z_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(2));

    t_0_z_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(3));

    t_0_z_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(4));

    t_0_z_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(5));

    t_0_z_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(6));

    t_0_z_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(7));

    t_0_z_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(8));

    t_0_z_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(9));

    t_0_y_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(10));

    t_0_y_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(11));

    t_0_y_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(12));

    t_0_y_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(13));

    t_0_y_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(14));

    t_0_y_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(15));

    t_0_y_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(16));

    t_0_y_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(17));

    t_0_y_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(18));

    t_0_y_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(19));

    t_0_x_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(20));

    t_0_x_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(21));

    t_0_x_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(22));

    t_0_x_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(23));

    t_0_x_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(24));

    t_0_x_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(25));

    t_0_x_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(26));

    t_0_x_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(27));

    t_0_x_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(28));

    t_0_x_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(29));

    // set up [SSSD]^(1) integral components

    t_0_0_0_zz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(0));

    t_0_0_0_yz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(1));

    t_0_0_0_yy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(2));

    t_0_0_0_xz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(3));

    t_0_0_0_xy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(4));

    t_0_0_0_xx_1 = intsBufferSSSD1.data(intsIndexesSSSD1(5));

    // set up [SSSF]^(0) integral components

    t_0_0_0_zzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(0));

    t_0_0_0_yzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(1));

    t_0_0_0_yyz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(2));

    t_0_0_0_yyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(3));

    t_0_0_0_xzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(4));

    t_0_0_0_xyz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(5));

    t_0_0_0_xyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(6));

    t_0_0_0_xxz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(7));

    t_0_0_0_xxy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(8));

    t_0_0_0_xxx_0 = intsBufferSSSF0.data(intsIndexesSSSF0(9));

    // set up [SSSF]^(1) integral components

    t_0_0_0_zzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(0));

    t_0_0_0_yzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(1));

    t_0_0_0_yyz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(2));

    t_0_0_0_yyy_1 = intsBufferSSSF1.data(intsIndexesSSSF1(3));

    t_0_0_0_xzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(4));

    t_0_0_0_xyz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(5));

    t_0_0_0_xyy_1 = intsBufferSSSF1.data(intsIndexesSSSF1(6));

    t_0_0_0_xxz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(7));

    t_0_0_0_xxy_1 = intsBufferSSSF1.data(intsIndexesSSSF1(8));

    t_0_0_0_xxx_1 = intsBufferSSSF1.data(intsIndexesSSSF1(9));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxx_1, t_0_0_0_xxy_0, t_0_0_0_xxy_1,\
                               t_0_0_0_xxz_0, t_0_0_0_xxz_1, t_0_0_0_xy_1, t_0_0_0_xyy_0,\
                               t_0_0_0_xyy_1, t_0_0_0_xyz_0, t_0_0_0_xyz_1, t_0_0_0_xz_1,\
                               t_0_0_0_xzz_0, t_0_0_0_xzz_1, t_0_0_0_yy_1, t_0_0_0_yyy_0,\
                               t_0_0_0_yyy_1, t_0_0_0_yyz_0, t_0_0_0_yyz_1, t_0_0_0_yz_1,\
                               t_0_0_0_yzz_0, t_0_0_0_yzz_1, t_0_0_0_zz_1, t_0_0_0_zzz_0,\
                               t_0_0_0_zzz_1, t_0_x_0_xxx_0, t_0_x_0_xxy_0, t_0_x_0_xxz_0,\
                               t_0_x_0_xyy_0, t_0_x_0_xyz_0, t_0_x_0_xzz_0, t_0_x_0_yyy_0,\
                               t_0_x_0_yyz_0, t_0_x_0_yzz_0, t_0_x_0_zzz_0, t_0_y_0_xxx_0,\
                               t_0_y_0_xxy_0, t_0_y_0_xxz_0, t_0_y_0_xyy_0, t_0_y_0_xyz_0,\
                               t_0_y_0_xzz_0, t_0_y_0_yyy_0, t_0_y_0_yyz_0, t_0_y_0_yzz_0,\
                               t_0_y_0_zzz_0, t_0_z_0_xxx_0, t_0_z_0_xxy_0, t_0_z_0_xxz_0,\
                               t_0_z_0_xyy_0, t_0_z_0_xyz_0, t_0_z_0_xzz_0, t_0_z_0_yyy_0,\
                               t_0_z_0_yyz_0, t_0_z_0_yzz_0, t_0_z_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzz_0[i] += rpb_z[i] * t_0_0_0_zzz_0[i] + rwp_z[i] * t_0_0_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_z_0_yzz_0[i] += rpb_z[i] * t_0_0_0_yzz_0[i] + rwp_z[i] * t_0_0_0_yzz_1[i] + fze_0[i] * t_0_0_0_yz_1[i];

            t_0_z_0_yyz_0[i] += rpb_z[i] * t_0_0_0_yyz_0[i] + rwp_z[i] * t_0_0_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_z_0_yyy_0[i] += rpb_z[i] * t_0_0_0_yyy_0[i] + rwp_z[i] * t_0_0_0_yyy_1[i];

            t_0_z_0_xzz_0[i] += rpb_z[i] * t_0_0_0_xzz_0[i] + rwp_z[i] * t_0_0_0_xzz_1[i] + fze_0[i] * t_0_0_0_xz_1[i];

            t_0_z_0_xyz_0[i] += rpb_z[i] * t_0_0_0_xyz_0[i] + rwp_z[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xy_1[i];

            t_0_z_0_xyy_0[i] += rpb_z[i] * t_0_0_0_xyy_0[i] + rwp_z[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxz_0[i] += rpb_z[i] * t_0_0_0_xxz_0[i] + rwp_z[i] * t_0_0_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xx_1[i];

            t_0_z_0_xxy_0[i] += rpb_z[i] * t_0_0_0_xxy_0[i] + rwp_z[i] * t_0_0_0_xxy_1[i];

            t_0_z_0_xxx_0[i] += rpb_z[i] * t_0_0_0_xxx_0[i] + rwp_z[i] * t_0_0_0_xxx_1[i];

            t_0_y_0_zzz_0[i] += rpb_y[i] * t_0_0_0_zzz_0[i] + rwp_y[i] * t_0_0_0_zzz_1[i];

            t_0_y_0_yzz_0[i] += rpb_y[i] * t_0_0_0_yzz_0[i] + rwp_y[i] * t_0_0_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_y_0_yyz_0[i] += rpb_y[i] * t_0_0_0_yyz_0[i] + rwp_y[i] * t_0_0_0_yyz_1[i] + fze_0[i] * t_0_0_0_yz_1[i];

            t_0_y_0_yyy_0[i] += rpb_y[i] * t_0_0_0_yyy_0[i] + rwp_y[i] * t_0_0_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_y_0_xzz_0[i] += rpb_y[i] * t_0_0_0_xzz_0[i] + rwp_y[i] * t_0_0_0_xzz_1[i];

            t_0_y_0_xyz_0[i] += rpb_y[i] * t_0_0_0_xyz_0[i] + rwp_y[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xz_1[i];

            t_0_y_0_xyy_0[i] += rpb_y[i] * t_0_0_0_xyy_0[i] + rwp_y[i] * t_0_0_0_xyy_1[i] + fze_0[i] * t_0_0_0_xy_1[i];

            t_0_y_0_xxz_0[i] += rpb_y[i] * t_0_0_0_xxz_0[i] + rwp_y[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_xxy_0[i] += rpb_y[i] * t_0_0_0_xxy_0[i] + rwp_y[i] * t_0_0_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xx_1[i];

            t_0_y_0_xxx_0[i] += rpb_y[i] * t_0_0_0_xxx_0[i] + rwp_y[i] * t_0_0_0_xxx_1[i];

            t_0_x_0_zzz_0[i] += rpb_x[i] * t_0_0_0_zzz_0[i] + rwp_x[i] * t_0_0_0_zzz_1[i];

            t_0_x_0_yzz_0[i] += rpb_x[i] * t_0_0_0_yzz_0[i] + rwp_x[i] * t_0_0_0_yzz_1[i];

            t_0_x_0_yyz_0[i] += rpb_x[i] * t_0_0_0_yyz_0[i] + rwp_x[i] * t_0_0_0_yyz_1[i];

            t_0_x_0_yyy_0[i] += rpb_x[i] * t_0_0_0_yyy_0[i] + rwp_x[i] * t_0_0_0_yyy_1[i];

            t_0_x_0_xzz_0[i] += rpb_x[i] * t_0_0_0_xzz_0[i] + rwp_x[i] * t_0_0_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_x_0_xyz_0[i] += rpb_x[i] * t_0_0_0_xyz_0[i] + rwp_x[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yz_1[i];

            t_0_x_0_xyy_0[i] += rpb_x[i] * t_0_0_0_xyy_0[i] + rwp_x[i] * t_0_0_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_x_0_xxz_0[i] += rpb_x[i] * t_0_0_0_xxz_0[i] + rwp_x[i] * t_0_0_0_xxz_1[i] + fze_0[i] * t_0_0_0_xz_1[i];

            t_0_x_0_xxy_0[i] += rpb_x[i] * t_0_0_0_xxy_0[i] + rwp_x[i] * t_0_0_0_xxy_1[i] + fze_0[i] * t_0_0_0_xy_1[i];

            t_0_x_0_xxx_0[i] += rpb_x[i] * t_0_0_0_xxx_0[i] + rwp_x[i] * t_0_0_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxx_1, t_0_0_0_xxy_0, t_0_0_0_xxy_1,\
                               t_0_0_0_xxz_0, t_0_0_0_xxz_1, t_0_0_0_xy_1, t_0_0_0_xyy_0,\
                               t_0_0_0_xyy_1, t_0_0_0_xyz_0, t_0_0_0_xyz_1, t_0_0_0_xz_1,\
                               t_0_0_0_xzz_0, t_0_0_0_xzz_1, t_0_0_0_yy_1, t_0_0_0_yyy_0,\
                               t_0_0_0_yyy_1, t_0_0_0_yyz_0, t_0_0_0_yyz_1, t_0_0_0_yz_1,\
                               t_0_0_0_yzz_0, t_0_0_0_yzz_1, t_0_0_0_zz_1, t_0_0_0_zzz_0,\
                               t_0_0_0_zzz_1, t_0_x_0_xxx_0, t_0_x_0_xxy_0, t_0_x_0_xxz_0,\
                               t_0_x_0_xyy_0, t_0_x_0_xyz_0, t_0_x_0_xzz_0, t_0_x_0_yyy_0,\
                               t_0_x_0_yyz_0, t_0_x_0_yzz_0, t_0_x_0_zzz_0, t_0_y_0_xxx_0,\
                               t_0_y_0_xxy_0, t_0_y_0_xxz_0, t_0_y_0_xyy_0, t_0_y_0_xyz_0,\
                               t_0_y_0_xzz_0, t_0_y_0_yyy_0, t_0_y_0_yyz_0, t_0_y_0_yzz_0,\
                               t_0_y_0_zzz_0, t_0_z_0_xxx_0, t_0_z_0_xxy_0, t_0_z_0_xxz_0,\
                               t_0_z_0_xyy_0, t_0_z_0_xyz_0, t_0_z_0_xzz_0, t_0_z_0_yyy_0,\
                               t_0_z_0_yyz_0, t_0_z_0_yzz_0, t_0_z_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzz_0[i] = rpb_z[i] * t_0_0_0_zzz_0[i] + rwp_z[i] * t_0_0_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_z_0_yzz_0[i] = rpb_z[i] * t_0_0_0_yzz_0[i] + rwp_z[i] * t_0_0_0_yzz_1[i] + fze_0[i] * t_0_0_0_yz_1[i];

            t_0_z_0_yyz_0[i] = rpb_z[i] * t_0_0_0_yyz_0[i] + rwp_z[i] * t_0_0_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_z_0_yyy_0[i] = rpb_z[i] * t_0_0_0_yyy_0[i] + rwp_z[i] * t_0_0_0_yyy_1[i];

            t_0_z_0_xzz_0[i] = rpb_z[i] * t_0_0_0_xzz_0[i] + rwp_z[i] * t_0_0_0_xzz_1[i] + fze_0[i] * t_0_0_0_xz_1[i];

            t_0_z_0_xyz_0[i] = rpb_z[i] * t_0_0_0_xyz_0[i] + rwp_z[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xy_1[i];

            t_0_z_0_xyy_0[i] = rpb_z[i] * t_0_0_0_xyy_0[i] + rwp_z[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxz_0[i] = rpb_z[i] * t_0_0_0_xxz_0[i] + rwp_z[i] * t_0_0_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xx_1[i];

            t_0_z_0_xxy_0[i] = rpb_z[i] * t_0_0_0_xxy_0[i] + rwp_z[i] * t_0_0_0_xxy_1[i];

            t_0_z_0_xxx_0[i] = rpb_z[i] * t_0_0_0_xxx_0[i] + rwp_z[i] * t_0_0_0_xxx_1[i];

            t_0_y_0_zzz_0[i] = rpb_y[i] * t_0_0_0_zzz_0[i] + rwp_y[i] * t_0_0_0_zzz_1[i];

            t_0_y_0_yzz_0[i] = rpb_y[i] * t_0_0_0_yzz_0[i] + rwp_y[i] * t_0_0_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_y_0_yyz_0[i] = rpb_y[i] * t_0_0_0_yyz_0[i] + rwp_y[i] * t_0_0_0_yyz_1[i] + fze_0[i] * t_0_0_0_yz_1[i];

            t_0_y_0_yyy_0[i] = rpb_y[i] * t_0_0_0_yyy_0[i] + rwp_y[i] * t_0_0_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_y_0_xzz_0[i] = rpb_y[i] * t_0_0_0_xzz_0[i] + rwp_y[i] * t_0_0_0_xzz_1[i];

            t_0_y_0_xyz_0[i] = rpb_y[i] * t_0_0_0_xyz_0[i] + rwp_y[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xz_1[i];

            t_0_y_0_xyy_0[i] = rpb_y[i] * t_0_0_0_xyy_0[i] + rwp_y[i] * t_0_0_0_xyy_1[i] + fze_0[i] * t_0_0_0_xy_1[i];

            t_0_y_0_xxz_0[i] = rpb_y[i] * t_0_0_0_xxz_0[i] + rwp_y[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_xxy_0[i] = rpb_y[i] * t_0_0_0_xxy_0[i] + rwp_y[i] * t_0_0_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xx_1[i];

            t_0_y_0_xxx_0[i] = rpb_y[i] * t_0_0_0_xxx_0[i] + rwp_y[i] * t_0_0_0_xxx_1[i];

            t_0_x_0_zzz_0[i] = rpb_x[i] * t_0_0_0_zzz_0[i] + rwp_x[i] * t_0_0_0_zzz_1[i];

            t_0_x_0_yzz_0[i] = rpb_x[i] * t_0_0_0_yzz_0[i] + rwp_x[i] * t_0_0_0_yzz_1[i];

            t_0_x_0_yyz_0[i] = rpb_x[i] * t_0_0_0_yyz_0[i] + rwp_x[i] * t_0_0_0_yyz_1[i];

            t_0_x_0_yyy_0[i] = rpb_x[i] * t_0_0_0_yyy_0[i] + rwp_x[i] * t_0_0_0_yyy_1[i];

            t_0_x_0_xzz_0[i] = rpb_x[i] * t_0_0_0_xzz_0[i] + rwp_x[i] * t_0_0_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_zz_1[i];

            t_0_x_0_xyz_0[i] = rpb_x[i] * t_0_0_0_xyz_0[i] + rwp_x[i] * t_0_0_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yz_1[i];

            t_0_x_0_xyy_0[i] = rpb_x[i] * t_0_0_0_xyy_0[i] + rwp_x[i] * t_0_0_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yy_1[i];

            t_0_x_0_xxz_0[i] = rpb_x[i] * t_0_0_0_xxz_0[i] + rwp_x[i] * t_0_0_0_xxz_1[i] + fze_0[i] * t_0_0_0_xz_1[i];

            t_0_x_0_xxy_0[i] = rpb_x[i] * t_0_0_0_xxy_0[i] + rwp_x[i] * t_0_0_0_xxy_1[i] + fze_0[i] * t_0_0_0_xy_1[i];

            t_0_x_0_xxx_0[i] = rpb_x[i] * t_0_0_0_xxx_0[i] + rwp_x[i] * t_0_0_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xx_1[i];
        }
    }
}


} // derirec namespace
