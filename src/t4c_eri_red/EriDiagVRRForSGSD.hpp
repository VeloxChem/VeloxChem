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
compHostVRRForSGSD_V0(      BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD0,
                      const BufferHostXY<T>&      intsBufferSDSD0,
                      const BufferHostX<int32_t>& intsIndexesSDSD0,
                      const BufferHostXY<T>&      intsBufferSDSD1,
                      const BufferHostX<int32_t>& intsIndexesSDSD1,
                      const BufferHostXY<T>&      intsBufferSFSP1,
                      const BufferHostX<int32_t>& intsIndexesSFSP1,
                      const BufferHostXY<T>&      intsBufferSFSD0,
                      const BufferHostX<int32_t>& intsIndexesSFSD0,
                      const BufferHostXY<T>&      intsBufferSFSD1,
                      const BufferHostX<int32_t>& intsIndexesSFSD1,
                      const T*                    osFactorsZeta,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

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

    // set up [SGSD]^(0) integral components

    t_0_zzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(0));

    t_0_yzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(1));

    t_0_yzzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(2));

    t_0_yyzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(3));

    t_0_yyzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(4));

    t_0_yyzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(5));

    t_0_yyyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(6));

    t_0_yyyz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(7));

    t_0_yyyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(8));

    t_0_xzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(9));

    t_0_xzzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(10));

    t_0_xyzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(11));

    t_0_xyzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(12));

    t_0_xyzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(13));

    t_0_xyzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(14));

    t_0_xyyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(15));

    t_0_xyyz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(16));

    t_0_xyyz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(17));

    t_0_xyyz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(18));

    t_0_xyyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(19));

    t_0_xyyy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(20));

    t_0_xxzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(21));

    t_0_xxzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(22));

    t_0_xxzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(23));

    t_0_xxyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(24));

    t_0_xxyz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(25));

    t_0_xxyz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(26));

    t_0_xxyz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(27));

    t_0_xxyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(28));

    t_0_xxyy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(29));

    t_0_xxyy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(30));

    t_0_xxxz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(31));

    t_0_xxxz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(32));

    t_0_xxxy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(33));

    t_0_xxxy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(34));

    t_0_xxxx_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(35));

    // set up [SDSD]^(0) integral components

    t_0_zz_0_zz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(0));

    t_0_zz_0_yz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(1));

    t_0_zz_0_xz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(2));

    t_0_yy_0_yy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(3));

    t_0_yy_0_xy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(4));

    t_0_xx_0_xx_0 = intsBufferSDSD0.data(intsIndexesSDSD0(5));

    // set up [SDSD]^(1) integral components

    t_0_zz_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(0));

    t_0_zz_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(1));

    t_0_zz_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(2));

    t_0_yy_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(3));

    t_0_yy_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(4));

    t_0_xx_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(5));

    // set up [SFSP]^(1) integral components

    t_0_zzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(0));

    t_0_yzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(1));

    t_0_yzz_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(2));

    t_0_yyz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(3));

    t_0_yyy_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(4));

    t_0_xzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(5));

    t_0_xyy_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(6));

    t_0_xxz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(7));

    t_0_xxx_0_x_1 = intsBufferSFSP1.data(intsIndexesSFSP1(8));

    // set up [SFSD]^(0) integral components

    t_0_zzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(0));

    t_0_zzz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(1));

    t_0_zzz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(2));

    t_0_yzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(3));

    t_0_yzz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(4));

    t_0_yzz_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(5));

    t_0_yyz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(6));

    t_0_yyz_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(7));

    t_0_yyz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(8));

    t_0_yyy_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(9));

    t_0_yyy_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(10));

    t_0_yyy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(11));

    t_0_xzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(12));

    t_0_xzz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(13));

    t_0_xyy_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(14));

    t_0_xyy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(15));

    t_0_xxz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(16));

    t_0_xxz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(17));

    t_0_xxz_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(18));

    t_0_xxy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(19));

    t_0_xxy_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(20));

    t_0_xxx_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(21));

    t_0_xxx_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(22));

    t_0_xxx_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(23));

    // set up [SFSD]^(1) integral components

    t_0_zzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(0));

    t_0_zzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(1));

    t_0_zzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(2));

    t_0_yzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(3));

    t_0_yzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(4));

    t_0_yzz_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(5));

    t_0_yyz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(6));

    t_0_yyz_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(7));

    t_0_yyz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(8));

    t_0_yyy_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(9));

    t_0_yyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(10));

    t_0_yyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(11));

    t_0_xzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(12));

    t_0_xzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(13));

    t_0_xyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(14));

    t_0_xyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(15));

    t_0_xxz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(16));

    t_0_xxz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(17));

    t_0_xxz_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(18));

    t_0_xxy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(19));

    t_0_xxy_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(20));

    t_0_xxx_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(21));

    t_0_xxx_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(22));

    t_0_xxx_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(23));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xxx_0_x_1, t_0_xxx_0_xx_0,\
                               t_0_xxx_0_xx_1, t_0_xxx_0_xy_0, t_0_xxx_0_xy_1, t_0_xxx_0_xz_0,\
                               t_0_xxx_0_xz_1, t_0_xxxx_0_xx_0, t_0_xxxy_0_xx_0, t_0_xxxy_0_xy_0,\
                               t_0_xxxz_0_xx_0, t_0_xxxz_0_xz_0, t_0_xxy_0_xx_0, t_0_xxy_0_xx_1,\
                               t_0_xxy_0_xy_0, t_0_xxy_0_xy_1, t_0_xxyy_0_xx_0, t_0_xxyy_0_xy_0,\
                               t_0_xxyy_0_yy_0, t_0_xxyz_0_xx_0, t_0_xxyz_0_xy_0, t_0_xxyz_0_xz_0,\
                               t_0_xxyz_0_yz_0, t_0_xxz_0_xx_0, t_0_xxz_0_xx_1, t_0_xxz_0_xz_0,\
                               t_0_xxz_0_xz_1, t_0_xxz_0_yz_0, t_0_xxz_0_yz_1, t_0_xxz_0_z_1,\
                               t_0_xxzz_0_xx_0, t_0_xxzz_0_xz_0, t_0_xxzz_0_zz_0, t_0_xyy_0_xy_0,\
                               t_0_xyy_0_xy_1, t_0_xyy_0_y_1, t_0_xyy_0_yy_0, t_0_xyy_0_yy_1,\
                               t_0_xyyy_0_xy_0, t_0_xyyy_0_yy_0, t_0_xyyz_0_xy_0, t_0_xyyz_0_xz_0,\
                               t_0_xyyz_0_yy_0, t_0_xyyz_0_yz_0, t_0_xyzz_0_xy_0, t_0_xyzz_0_xz_0,\
                               t_0_xyzz_0_yz_0, t_0_xyzz_0_zz_0, t_0_xzz_0_xz_0, t_0_xzz_0_xz_1,\
                               t_0_xzz_0_z_1, t_0_xzz_0_zz_0, t_0_xzz_0_zz_1, t_0_xzzz_0_xz_0,\
                               t_0_xzzz_0_zz_0, t_0_yy_0_xy_0, t_0_yy_0_xy_1, t_0_yy_0_yy_0,\
                               t_0_yy_0_yy_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1, t_0_yyy_0_y_1,\
                               t_0_yyy_0_yy_0, t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1,\
                               t_0_yyyy_0_yy_0, t_0_yyyz_0_yy_0, t_0_yyyz_0_yz_0, t_0_yyz_0_xz_0,\
                               t_0_yyz_0_xz_1, t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyz_0_yz_0,\
                               t_0_yyz_0_yz_1, t_0_yyz_0_z_1, t_0_yyzz_0_yy_0, t_0_yyzz_0_yz_0,\
                               t_0_yyzz_0_zz_0, t_0_yzz_0_xy_0, t_0_yzz_0_xy_1, t_0_yzz_0_y_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_z_1, t_0_yzz_0_zz_0,\
                               t_0_yzz_0_zz_1, t_0_yzzz_0_yz_0, t_0_yzzz_0_zz_0, t_0_zz_0_xz_0,\
                               t_0_zz_0_xz_1, t_0_zz_0_yz_0, t_0_zz_0_yz_1, t_0_zz_0_zz_0,\
                               t_0_zz_0_zz_1, t_0_zzz_0_xz_0, t_0_zzz_0_xz_1, t_0_zzz_0_yz_0,\
                               t_0_zzz_0_yz_1, t_0_zzz_0_z_1, t_0_zzz_0_zz_0, t_0_zzz_0_zz_1,\
                               t_0_zzzz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zz_0[i] += rpb_z[i] * t_0_zzz_0_zz_0[i] + rwp_z[i] * t_0_zzz_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zz_1[i] + fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yzzz_0_zz_0[i] += rpb_y[i] * t_0_zzz_0_zz_0[i] + rwp_y[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yz_0[i] += rpb_y[i] * t_0_zzz_0_yz_0[i] + rwp_y[i] * t_0_zzz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yyzz_0_zz_0[i] += rpb_y[i] * t_0_yzz_0_zz_0[i] + rwp_y[i] * t_0_yzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_yyzz_0_yz_0[i] += rpb_y[i] * t_0_yzz_0_yz_0[i] + rwp_y[i] * t_0_yzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_z_1[i];

            t_0_yyzz_0_yy_0[i] += rpb_z[i] * t_0_yyz_0_yy_0[i] + rwp_z[i] * t_0_yyz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyyz_0_yz_0[i] += rpb_z[i] * t_0_yyy_0_yz_0[i] + rwp_z[i] * t_0_yyy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyz_0_yy_0[i] += rpb_z[i] * t_0_yyy_0_yy_0[i] + rwp_z[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyy_0_yy_0[i] += rpb_y[i] * t_0_yyy_0_yy_0[i] + rwp_y[i] * t_0_yyy_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yy_1[i] + fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xzzz_0_zz_0[i] += rpb_x[i] * t_0_zzz_0_zz_0[i] + rwp_x[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_xz_0[i] += rpb_x[i] * t_0_zzz_0_xz_0[i] + rwp_x[i] * t_0_zzz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_xyzz_0_zz_0[i] += rpb_x[i] * t_0_yzz_0_zz_0[i] + rwp_x[i] * t_0_yzz_0_zz_1[i];

            t_0_xyzz_0_yz_0[i] += rpb_x[i] * t_0_yzz_0_yz_0[i] + rwp_x[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_xz_0[i] += rpb_y[i] * t_0_xzz_0_xz_0[i] + rwp_y[i] * t_0_xzz_0_xz_1[i];

            t_0_xyzz_0_xy_0[i] += rpb_x[i] * t_0_yzz_0_xy_0[i] + rwp_x[i] * t_0_yzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_y_1[i];

            t_0_xyyz_0_yz_0[i] += rpb_x[i] * t_0_yyz_0_yz_0[i] + rwp_x[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_yy_0[i] += rpb_x[i] * t_0_yyz_0_yy_0[i] + rwp_x[i] * t_0_yyz_0_yy_1[i];

            t_0_xyyz_0_xz_0[i] += rpb_x[i] * t_0_yyz_0_xz_0[i] + rwp_x[i] * t_0_yyz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_z_1[i];

            t_0_xyyz_0_xy_0[i] += rpb_z[i] * t_0_xyy_0_xy_0[i] + rwp_z[i] * t_0_xyy_0_xy_1[i];

            t_0_xyyy_0_yy_0[i] += rpb_x[i] * t_0_yyy_0_yy_0[i] + rwp_x[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xy_0[i] += rpb_x[i] * t_0_yyy_0_xy_0[i] + rwp_x[i] * t_0_yyy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xxzz_0_zz_0[i] += rpb_x[i] * t_0_xzz_0_zz_0[i] + rwp_x[i] * t_0_xzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_xxzz_0_xz_0[i] += rpb_x[i] * t_0_xzz_0_xz_0[i] + rwp_x[i] * t_0_xzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_z_1[i];

            t_0_xxzz_0_xx_0[i] += rpb_z[i] * t_0_xxz_0_xx_0[i] + rwp_z[i] * t_0_xxz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxyz_0_yz_0[i] += rpb_y[i] * t_0_xxz_0_yz_0[i] + rwp_y[i] * t_0_xxz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_z_1[i];

            t_0_xxyz_0_xz_0[i] += rpb_y[i] * t_0_xxz_0_xz_0[i] + rwp_y[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xy_0[i] += rpb_z[i] * t_0_xxy_0_xy_0[i] + rwp_z[i] * t_0_xxy_0_xy_1[i];

            t_0_xxyz_0_xx_0[i] += rpb_y[i] * t_0_xxz_0_xx_0[i] + rwp_y[i] * t_0_xxz_0_xx_1[i];

            t_0_xxyy_0_yy_0[i] += rpb_x[i] * t_0_xyy_0_yy_0[i] + rwp_x[i] * t_0_xyy_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_xxyy_0_xy_0[i] += rpb_x[i] * t_0_xyy_0_xy_0[i] + rwp_x[i] * t_0_xyy_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_y_1[i];

            t_0_xxyy_0_xx_0[i] += rpb_y[i] * t_0_xxy_0_xx_0[i] + rwp_y[i] * t_0_xxy_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxxz_0_xz_0[i] += rpb_z[i] * t_0_xxx_0_xz_0[i] + rwp_z[i] * t_0_xxx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxz_0_xx_0[i] += rpb_z[i] * t_0_xxx_0_xx_0[i] + rwp_z[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_xy_0[i] += rpb_y[i] * t_0_xxx_0_xy_0[i] + rwp_y[i] * t_0_xxx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxy_0_xx_0[i] += rpb_y[i] * t_0_xxx_0_xx_0[i] + rwp_y[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxx_0_xx_0[i] += rpb_x[i] * t_0_xxx_0_xx_0[i] + rwp_x[i] * t_0_xxx_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xx_1[i] + fze_0[i] * t_0_xxx_0_x_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xxx_0_x_1, t_0_xxx_0_xx_0,\
                               t_0_xxx_0_xx_1, t_0_xxx_0_xy_0, t_0_xxx_0_xy_1, t_0_xxx_0_xz_0,\
                               t_0_xxx_0_xz_1, t_0_xxxx_0_xx_0, t_0_xxxy_0_xx_0, t_0_xxxy_0_xy_0,\
                               t_0_xxxz_0_xx_0, t_0_xxxz_0_xz_0, t_0_xxy_0_xx_0, t_0_xxy_0_xx_1,\
                               t_0_xxy_0_xy_0, t_0_xxy_0_xy_1, t_0_xxyy_0_xx_0, t_0_xxyy_0_xy_0,\
                               t_0_xxyy_0_yy_0, t_0_xxyz_0_xx_0, t_0_xxyz_0_xy_0, t_0_xxyz_0_xz_0,\
                               t_0_xxyz_0_yz_0, t_0_xxz_0_xx_0, t_0_xxz_0_xx_1, t_0_xxz_0_xz_0,\
                               t_0_xxz_0_xz_1, t_0_xxz_0_yz_0, t_0_xxz_0_yz_1, t_0_xxz_0_z_1,\
                               t_0_xxzz_0_xx_0, t_0_xxzz_0_xz_0, t_0_xxzz_0_zz_0, t_0_xyy_0_xy_0,\
                               t_0_xyy_0_xy_1, t_0_xyy_0_y_1, t_0_xyy_0_yy_0, t_0_xyy_0_yy_1,\
                               t_0_xyyy_0_xy_0, t_0_xyyy_0_yy_0, t_0_xyyz_0_xy_0, t_0_xyyz_0_xz_0,\
                               t_0_xyyz_0_yy_0, t_0_xyyz_0_yz_0, t_0_xyzz_0_xy_0, t_0_xyzz_0_xz_0,\
                               t_0_xyzz_0_yz_0, t_0_xyzz_0_zz_0, t_0_xzz_0_xz_0, t_0_xzz_0_xz_1,\
                               t_0_xzz_0_z_1, t_0_xzz_0_zz_0, t_0_xzz_0_zz_1, t_0_xzzz_0_xz_0,\
                               t_0_xzzz_0_zz_0, t_0_yy_0_xy_0, t_0_yy_0_xy_1, t_0_yy_0_yy_0,\
                               t_0_yy_0_yy_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1, t_0_yyy_0_y_1,\
                               t_0_yyy_0_yy_0, t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1,\
                               t_0_yyyy_0_yy_0, t_0_yyyz_0_yy_0, t_0_yyyz_0_yz_0, t_0_yyz_0_xz_0,\
                               t_0_yyz_0_xz_1, t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyz_0_yz_0,\
                               t_0_yyz_0_yz_1, t_0_yyz_0_z_1, t_0_yyzz_0_yy_0, t_0_yyzz_0_yz_0,\
                               t_0_yyzz_0_zz_0, t_0_yzz_0_xy_0, t_0_yzz_0_xy_1, t_0_yzz_0_y_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_z_1, t_0_yzz_0_zz_0,\
                               t_0_yzz_0_zz_1, t_0_yzzz_0_yz_0, t_0_yzzz_0_zz_0, t_0_zz_0_xz_0,\
                               t_0_zz_0_xz_1, t_0_zz_0_yz_0, t_0_zz_0_yz_1, t_0_zz_0_zz_0,\
                               t_0_zz_0_zz_1, t_0_zzz_0_xz_0, t_0_zzz_0_xz_1, t_0_zzz_0_yz_0,\
                               t_0_zzz_0_yz_1, t_0_zzz_0_z_1, t_0_zzz_0_zz_0, t_0_zzz_0_zz_1,\
                               t_0_zzzz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zz_0[i] = rpb_z[i] * t_0_zzz_0_zz_0[i] + rwp_z[i] * t_0_zzz_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zz_1[i] + fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yzzz_0_zz_0[i] = rpb_y[i] * t_0_zzz_0_zz_0[i] + rwp_y[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yz_0[i] = rpb_y[i] * t_0_zzz_0_yz_0[i] + rwp_y[i] * t_0_zzz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yyzz_0_zz_0[i] = rpb_y[i] * t_0_yzz_0_zz_0[i] + rwp_y[i] * t_0_yzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_yyzz_0_yz_0[i] = rpb_y[i] * t_0_yzz_0_yz_0[i] + rwp_y[i] * t_0_yzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_z_1[i];

            t_0_yyzz_0_yy_0[i] = rpb_z[i] * t_0_yyz_0_yy_0[i] + rwp_z[i] * t_0_yyz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyyz_0_yz_0[i] = rpb_z[i] * t_0_yyy_0_yz_0[i] + rwp_z[i] * t_0_yyy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyz_0_yy_0[i] = rpb_z[i] * t_0_yyy_0_yy_0[i] + rwp_z[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyy_0_yy_0[i] = rpb_y[i] * t_0_yyy_0_yy_0[i] + rwp_y[i] * t_0_yyy_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yy_1[i] + fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xzzz_0_zz_0[i] = rpb_x[i] * t_0_zzz_0_zz_0[i] + rwp_x[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_xz_0[i] = rpb_x[i] * t_0_zzz_0_xz_0[i] + rwp_x[i] * t_0_zzz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_xyzz_0_zz_0[i] = rpb_x[i] * t_0_yzz_0_zz_0[i] + rwp_x[i] * t_0_yzz_0_zz_1[i];

            t_0_xyzz_0_yz_0[i] = rpb_x[i] * t_0_yzz_0_yz_0[i] + rwp_x[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_xz_0[i] = rpb_y[i] * t_0_xzz_0_xz_0[i] + rwp_y[i] * t_0_xzz_0_xz_1[i];

            t_0_xyzz_0_xy_0[i] = rpb_x[i] * t_0_yzz_0_xy_0[i] + rwp_x[i] * t_0_yzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_y_1[i];

            t_0_xyyz_0_yz_0[i] = rpb_x[i] * t_0_yyz_0_yz_0[i] + rwp_x[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_yy_0[i] = rpb_x[i] * t_0_yyz_0_yy_0[i] + rwp_x[i] * t_0_yyz_0_yy_1[i];

            t_0_xyyz_0_xz_0[i] = rpb_x[i] * t_0_yyz_0_xz_0[i] + rwp_x[i] * t_0_yyz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_z_1[i];

            t_0_xyyz_0_xy_0[i] = rpb_z[i] * t_0_xyy_0_xy_0[i] + rwp_z[i] * t_0_xyy_0_xy_1[i];

            t_0_xyyy_0_yy_0[i] = rpb_x[i] * t_0_yyy_0_yy_0[i] + rwp_x[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xy_0[i] = rpb_x[i] * t_0_yyy_0_xy_0[i] + rwp_x[i] * t_0_yyy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xxzz_0_zz_0[i] = rpb_x[i] * t_0_xzz_0_zz_0[i] + rwp_x[i] * t_0_xzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_xxzz_0_xz_0[i] = rpb_x[i] * t_0_xzz_0_xz_0[i] + rwp_x[i] * t_0_xzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_z_1[i];

            t_0_xxzz_0_xx_0[i] = rpb_z[i] * t_0_xxz_0_xx_0[i] + rwp_z[i] * t_0_xxz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxyz_0_yz_0[i] = rpb_y[i] * t_0_xxz_0_yz_0[i] + rwp_y[i] * t_0_xxz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_z_1[i];

            t_0_xxyz_0_xz_0[i] = rpb_y[i] * t_0_xxz_0_xz_0[i] + rwp_y[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xy_0[i] = rpb_z[i] * t_0_xxy_0_xy_0[i] + rwp_z[i] * t_0_xxy_0_xy_1[i];

            t_0_xxyz_0_xx_0[i] = rpb_y[i] * t_0_xxz_0_xx_0[i] + rwp_y[i] * t_0_xxz_0_xx_1[i];

            t_0_xxyy_0_yy_0[i] = rpb_x[i] * t_0_xyy_0_yy_0[i] + rwp_x[i] * t_0_xyy_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_xxyy_0_xy_0[i] = rpb_x[i] * t_0_xyy_0_xy_0[i] + rwp_x[i] * t_0_xyy_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_y_1[i];

            t_0_xxyy_0_xx_0[i] = rpb_y[i] * t_0_xxy_0_xx_0[i] + rwp_y[i] * t_0_xxy_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxxz_0_xz_0[i] = rpb_z[i] * t_0_xxx_0_xz_0[i] + rwp_z[i] * t_0_xxx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxz_0_xx_0[i] = rpb_z[i] * t_0_xxx_0_xx_0[i] + rwp_z[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_xy_0[i] = rpb_y[i] * t_0_xxx_0_xy_0[i] + rwp_y[i] * t_0_xxx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxy_0_xx_0[i] = rpb_y[i] * t_0_xxx_0_xx_0[i] + rwp_y[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxx_0_xx_0[i] = rpb_x[i] * t_0_xxx_0_xx_0[i] + rwp_x[i] * t_0_xxx_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xx_1[i] + fze_0[i] * t_0_xxx_0_x_1[i];
        }
    }
}


} // derirec namespace
