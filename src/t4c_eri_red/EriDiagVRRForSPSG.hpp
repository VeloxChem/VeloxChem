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
compHostVRRForSPSG_V0(      BufferHostXY<T>&      intsBufferSPSG,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSSSF1,
                      const BufferHostX<int32_t>& intsIndexesSSSF1,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
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

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_z_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_z_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_z_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_y_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_y_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(15));

    t_0_y_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(16));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(17));

    t_0_y_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(18));

    t_0_x_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(19));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(20));

    t_0_x_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(21));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(22));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(23));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(24));

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

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyy_0, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxzz_0, t_0_y_0_xxxy_0, t_0_y_0_xxyy_0, t_0_y_0_xxyz_0,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyz_0, t_0_y_0_xyzz_0, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyz_0, t_0_y_0_yyzz_0, t_0_z_0_xxxz_0, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyyz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0,\
                               t_0_z_0_yyyz_0, t_0_z_0_yyzz_0, t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] += rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] += rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] += rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_yyyz_0[i] += rpb_z[i] * t_0_0_0_yyyz_0[i] + rwp_z[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_z_0_xzzz_0[i] += rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] += rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xyyz_0[i] += rpb_z[i] * t_0_0_0_xyyz_0[i] + rwp_z[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxzz_0[i] += rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_z_0_xxyz_0[i] += rpb_z[i] * t_0_0_0_xxyz_0[i] + rwp_z[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_z_0_xxxz_0[i] += rpb_z[i] * t_0_0_0_xxxz_0[i] + rwp_z[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxx_1[i];

            t_0_y_0_yyzz_0[i] += rpb_y[i] * t_0_0_0_yyzz_0[i] + rwp_y[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_y_0_yyyz_0[i] += rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] += rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyzz_0[i] += rpb_y[i] * t_0_0_0_xyzz_0[i] + rwp_y[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_y_0_xyyz_0[i] += rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] += rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyz_0[i] += rpb_y[i] * t_0_0_0_xxyz_0[i] + rwp_y[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_xxyy_0[i] += rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_y_0_xxxy_0[i] += rpb_y[i] * t_0_0_0_xxxy_0[i] + rwp_y[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxx_1[i];

            t_0_x_0_xxzz_0[i] += rpb_x[i] * t_0_0_0_xxzz_0[i] + rwp_x[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_x_0_xxyz_0[i] += rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxyy_0[i] += rpb_x[i] * t_0_0_0_xxyy_0[i] + rwp_x[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_x_0_xxxz_0[i] += rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] += rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] += rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyy_0, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxzz_0, t_0_y_0_xxxy_0, t_0_y_0_xxyy_0, t_0_y_0_xxyz_0,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyz_0, t_0_y_0_xyzz_0, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyz_0, t_0_y_0_yyzz_0, t_0_z_0_xxxz_0, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyyz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0,\
                               t_0_z_0_yyyz_0, t_0_z_0_yyzz_0, t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] = rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] = rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] = rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_yyyz_0[i] = rpb_z[i] * t_0_0_0_yyyz_0[i] + rwp_z[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_z_0_xzzz_0[i] = rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] = rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xyyz_0[i] = rpb_z[i] * t_0_0_0_xyyz_0[i] + rwp_z[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxzz_0[i] = rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_z_0_xxyz_0[i] = rpb_z[i] * t_0_0_0_xxyz_0[i] + rwp_z[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_z_0_xxxz_0[i] = rpb_z[i] * t_0_0_0_xxxz_0[i] + rwp_z[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxx_1[i];

            t_0_y_0_yyzz_0[i] = rpb_y[i] * t_0_0_0_yyzz_0[i] + rwp_y[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_y_0_yyyz_0[i] = rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] = rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyzz_0[i] = rpb_y[i] * t_0_0_0_xyzz_0[i] + rwp_y[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_y_0_xyyz_0[i] = rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] = rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyz_0[i] = rpb_y[i] * t_0_0_0_xxyz_0[i] + rwp_y[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_xxyy_0[i] = rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_y_0_xxxy_0[i] = rpb_y[i] * t_0_0_0_xxxy_0[i] + rwp_y[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxx_1[i];

            t_0_x_0_xxzz_0[i] = rpb_x[i] * t_0_0_0_xxzz_0[i] + rwp_x[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_x_0_xxyz_0[i] = rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxyy_0[i] = rpb_x[i] * t_0_0_0_xxyy_0[i] + rwp_x[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_x_0_xxxz_0[i] = rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] = rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] = rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSPSG_V1(      BufferHostXY<T>&      intsBufferSPSG,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSSSF1,
                      const BufferHostX<int32_t>& intsIndexesSSSF1,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
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

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_z_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_y_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_x_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(15));

    t_0_x_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(16));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(17));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(18));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(19));

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

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyy_0, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxzz_0, t_0_y_0_xxyy_0, t_0_y_0_xyyy_0, t_0_y_0_xyyz_0,\
                               t_0_y_0_yyyy_0, t_0_y_0_yyyz_0, t_0_y_0_yyzz_0, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyyz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0,\
                               t_0_z_0_yyzz_0, t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] += rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] += rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] += rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_xzzz_0[i] += rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] += rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xyyz_0[i] += rpb_z[i] * t_0_0_0_xyyz_0[i] + rwp_z[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxzz_0[i] += rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_z_0_xxyz_0[i] += rpb_z[i] * t_0_0_0_xxyz_0[i] + rwp_z[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_y_0_yyzz_0[i] += rpb_y[i] * t_0_0_0_yyzz_0[i] + rwp_y[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_y_0_yyyz_0[i] += rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] += rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyyz_0[i] += rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] += rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyy_0[i] += rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxzz_0[i] += rpb_x[i] * t_0_0_0_xxzz_0[i] + rwp_x[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_x_0_xxyz_0[i] += rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxyy_0[i] += rpb_x[i] * t_0_0_0_xxyy_0[i] + rwp_x[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_x_0_xxxz_0[i] += rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] += rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] += rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyy_0, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxzz_0, t_0_y_0_xxyy_0, t_0_y_0_xyyy_0, t_0_y_0_xyyz_0,\
                               t_0_y_0_yyyy_0, t_0_y_0_yyyz_0, t_0_y_0_yyzz_0, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyyz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0,\
                               t_0_z_0_yyzz_0, t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] = rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] = rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] = rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_xzzz_0[i] = rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] = rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xyyz_0[i] = rpb_z[i] * t_0_0_0_xyyz_0[i] + rwp_z[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_z_0_xxzz_0[i] = rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_z_0_xxyz_0[i] = rpb_z[i] * t_0_0_0_xxyz_0[i] + rwp_z[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_y_0_yyzz_0[i] = rpb_y[i] * t_0_0_0_yyzz_0[i] + rwp_y[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_y_0_yyyz_0[i] = rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] = rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyyz_0[i] = rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] = rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyy_0[i] = rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxzz_0[i] = rpb_x[i] * t_0_0_0_xxzz_0[i] + rwp_x[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_x_0_xxyz_0[i] = rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxyy_0[i] = rpb_x[i] * t_0_0_0_xxyy_0[i] + rwp_x[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_x_0_xxxz_0[i] = rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] = rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] = rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSPSG_V2(      BufferHostXY<T>&      intsBufferSPSG,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSSSF1,
                      const BufferHostX<int32_t>& intsIndexesSSSF1,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
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

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

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

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up scaling factors

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyz_0, t_0_y_0_xxyy_0,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyz_0, t_0_y_0_yyyy_0, t_0_y_0_yyyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0, t_0_z_0_yyzz_0,\
                               t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] += rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] += rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] += rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_xzzz_0[i] += rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] += rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xxzz_0[i] += rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_yyyz_0[i] += rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] += rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyyz_0[i] += rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] += rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyy_0[i] += rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxyz_0[i] += rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxxz_0[i] += rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] += rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] += rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1,\
                               t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxy_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyy_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyz_1,\
                               t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyy_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1,\
                               t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1,\
                               t_0_0_0_zzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxz_0, t_0_x_0_xxyz_0, t_0_y_0_xxyy_0,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyz_0, t_0_y_0_yyyy_0, t_0_y_0_yyyz_0,\
                               t_0_z_0_xxzz_0, t_0_z_0_xyzz_0, t_0_z_0_xzzz_0, t_0_z_0_yyzz_0,\
                               t_0_z_0_yzzz_0, t_0_z_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zzzz_0[i] = rpb_z[i] * t_0_0_0_zzzz_0[i] + rwp_z[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_0_0_zzz_1[i];

            t_0_z_0_yzzz_0[i] = rpb_z[i] * t_0_0_0_yzzz_0[i] + rwp_z[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yzz_1[i];

            t_0_z_0_yyzz_0[i] = rpb_z[i] * t_0_0_0_yyzz_0[i] + rwp_z[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_z_0_xzzz_0[i] = rpb_z[i] * t_0_0_0_xzzz_0[i] + rwp_z[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xzz_1[i];

            t_0_z_0_xyzz_0[i] = rpb_z[i] * t_0_0_0_xyzz_0[i] + rwp_z[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_z_0_xxzz_0[i] = rpb_z[i] * t_0_0_0_xxzz_0[i] + rwp_z[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_y_0_yyyz_0[i] = rpb_y[i] * t_0_0_0_yyyz_0[i] + rwp_y[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_yyz_1[i];

            t_0_y_0_yyyy_0[i] = rpb_y[i] * t_0_0_0_yyyy_0[i] + rwp_y[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_0_0_yyy_1[i];

            t_0_y_0_xyyz_0[i] = rpb_y[i] * t_0_0_0_xyyz_0[i] + rwp_y[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_y_0_xyyy_0[i] = rpb_y[i] * t_0_0_0_xyyy_0[i] + rwp_y[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xyy_1[i];

            t_0_y_0_xxyy_0[i] = rpb_y[i] * t_0_0_0_xxyy_0[i] + rwp_y[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxyz_0[i] = rpb_x[i] * t_0_0_0_xxyz_0[i] + rwp_x[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_0_0_xyz_1[i];

            t_0_x_0_xxxz_0[i] = rpb_x[i] * t_0_0_0_xxxz_0[i] + rwp_x[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxz_1[i];

            t_0_x_0_xxxy_0[i] = rpb_x[i] * t_0_0_0_xxxy_0[i] + rwp_x[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_0_0_xxy_1[i];

            t_0_x_0_xxxx_0[i] = rpb_x[i] * t_0_0_0_xxxx_0[i] + rwp_x[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_0_0_xxx_1[i];
        }
    }
}


} // derirec namespace
