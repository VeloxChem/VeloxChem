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
compHostVRRForSGSG_V0(      BufferHostXY<T>&      intsBufferSGSG,
                      const BufferHostX<int32_t>& intsIndexesSGSG0,
                      const BufferHostXY<T>&      intsBufferSDSG0,
                      const BufferHostX<int32_t>& intsIndexesSDSG0,
                      const BufferHostXY<T>&      intsBufferSDSG1,
                      const BufferHostX<int32_t>& intsIndexesSDSG1,
                      const BufferHostXY<T>&      intsBufferSFSF1,
                      const BufferHostX<int32_t>& intsIndexesSFSF1,
                      const BufferHostXY<T>&      intsBufferSFSG0,
                      const BufferHostX<int32_t>& intsIndexesSFSG0,
                      const BufferHostXY<T>&      intsBufferSFSG1,
                      const BufferHostX<int32_t>& intsIndexesSFSG1,
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

    // set up [SGSG]^(0) integral components

    t_0_zzzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(0));

    t_0_yzzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(1));

    t_0_yyzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(2));

    t_0_yyyz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(3));

    t_0_yyyy_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(4));

    t_0_xzzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(5));

    t_0_xyzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(6));

    t_0_xyyz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(7));

    t_0_xyyy_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(8));

    t_0_xxzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(9));

    t_0_xxyz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(10));

    t_0_xxyy_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(11));

    t_0_xxxz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(12));

    t_0_xxxy_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(13));

    t_0_xxxx_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(14));

    // set up [SDSG]^(0) integral components

    t_0_zz_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(0));

    t_0_zz_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(1));

    t_0_zz_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(2));

    t_0_yy_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(3));

    t_0_yy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(4));

    t_0_xx_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(5));

    // set up [SDSG]^(1) integral components

    t_0_zz_0_zzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(0));

    t_0_zz_0_yyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(1));

    t_0_zz_0_xxzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(2));

    t_0_yy_0_yyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(3));

    t_0_yy_0_xxyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(4));

    t_0_xx_0_xxxx_1 = intsBufferSDSG1.data(intsIndexesSDSG1(5));

    // set up [SFSF]^(1) integral components

    t_0_zzz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(0));

    t_0_yzz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(1));

    t_0_yyz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(2));

    t_0_yyy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(3));

    t_0_xzz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(4));

    t_0_xyy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(5));

    t_0_xxz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(6));

    t_0_xxx_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(7));

    // set up [SFSG]^(0) integral components

    t_0_zzz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(0));

    t_0_zzz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(1));

    t_0_zzz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(2));

    t_0_yzz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(3));

    t_0_yzz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(4));

    t_0_yyz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(5));

    t_0_yyy_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(6));

    t_0_yyy_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(7));

    t_0_yyy_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(8));

    t_0_xzz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(9));

    t_0_xyy_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(10));

    t_0_xxz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(11));

    t_0_xxx_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(12));

    t_0_xxx_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(13));

    t_0_xxx_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(14));

    // set up [SFSG]^(1) integral components

    t_0_zzz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(0));

    t_0_zzz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(1));

    t_0_zzz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(2));

    t_0_yzz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(3));

    t_0_yzz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(4));

    t_0_yyz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(5));

    t_0_yyy_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(6));

    t_0_yyy_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(7));

    t_0_yyy_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(8));

    t_0_xzz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(9));

    t_0_xyy_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(10));

    t_0_xxz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(11));

    t_0_xxx_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(12));

    t_0_xxx_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(13));

    t_0_xxx_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(14));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1, t_0_xxx_0_xxx_1,\
                               t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1, t_0_xxx_0_xxxy_0,\
                               t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0, t_0_xxx_0_xxxz_1,\
                               t_0_xxxx_0_xxxx_0, t_0_xxxy_0_xxxy_0, t_0_xxxz_0_xxxz_0,\
                               t_0_xxyy_0_xxyy_0, t_0_xxyz_0_xxyz_0, t_0_xxz_0_xxyz_0,\
                               t_0_xxz_0_xxyz_1, t_0_xxz_0_xxz_1, t_0_xxzz_0_xxzz_0,\
                               t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1, t_0_xyy_0_xyy_1,\
                               t_0_xyyy_0_xyyy_0, t_0_xyyz_0_xyyz_0, t_0_xyzz_0_xyzz_0,\
                               t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1, t_0_xzz_0_xzz_1,\
                               t_0_xzzz_0_xzzz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_yyy_1, t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1,\
                               t_0_yyy_0_yyyz_0, t_0_yyy_0_yyyz_1, t_0_yyyy_0_yyyy_0,\
                               t_0_yyyz_0_yyyz_0, t_0_yyz_0_xyyz_0, t_0_yyz_0_xyyz_1,\
                               t_0_yyz_0_yyz_1, t_0_yyzz_0_yyzz_0, t_0_yzz_0_xyzz_0,\
                               t_0_yzz_0_xyzz_1, t_0_yzz_0_yyzz_0, t_0_yzz_0_yyzz_1,\
                               t_0_yzz_0_yzz_1, t_0_yzzz_0_yzzz_0, t_0_zz_0_xxzz_0,\
                               t_0_zz_0_xxzz_1, t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1, t_0_zz_0_zzzz_0,\
                               t_0_zz_0_zzzz_1, t_0_zzz_0_xzzz_0, t_0_zzz_0_xzzz_1,\
                               t_0_zzz_0_yzzz_0, t_0_zzz_0_yzzz_1, t_0_zzz_0_zzz_1,\
                               t_0_zzz_0_zzzz_0, t_0_zzz_0_zzzz_1, t_0_zzzz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzzz_0[i] += rpb_z[i] * t_0_zzz_0_zzzz_0[i] + rwp_z[i] * t_0_zzz_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yzzz_0[i] += rpb_y[i] * t_0_zzz_0_yzzz_0[i] + rwp_y[i] * t_0_zzz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yyzz_0_yyzz_0[i] += rpb_y[i] * t_0_yzz_0_yyzz_0[i] + rwp_y[i] * t_0_yzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_yyyz_0_yyyz_0[i] += rpb_z[i] * t_0_yyy_0_yyyz_0[i] + rwp_z[i] * t_0_yyy_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyy_0_yyyy_0[i] += rpb_y[i] * t_0_yyy_0_yyyy_0[i] + rwp_y[i] * t_0_yyy_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xzzz_0_xzzz_0[i] += rpb_x[i] * t_0_zzz_0_xzzz_0[i] + rwp_x[i] * t_0_zzz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_xyzz_0_xyzz_0[i] += rpb_x[i] * t_0_yzz_0_xyzz_0[i] + rwp_x[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyyz_0_xyyz_0[i] += rpb_x[i] * t_0_yyz_0_xyyz_0[i] + rwp_x[i] * t_0_yyz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyy_0_xyyy_0[i] += rpb_x[i] * t_0_yyy_0_xyyy_0[i] + rwp_x[i] * t_0_yyy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xxzz_0_xxzz_0[i] += rpb_x[i] * t_0_xzz_0_xxzz_0[i] + rwp_x[i] * t_0_xzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_xzz_0_xzz_1[i];

            t_0_xxyz_0_xxyz_0[i] += rpb_y[i] * t_0_xxz_0_xxyz_0[i] + rwp_y[i] * t_0_xxz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xxz_1[i];

            t_0_xxyy_0_xxyy_0[i] += rpb_x[i] * t_0_xyy_0_xxyy_0[i] + rwp_x[i] * t_0_xyy_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_xyy_0_xyy_1[i];

            t_0_xxxz_0_xxxz_0[i] += rpb_z[i] * t_0_xxx_0_xxxz_0[i] + rwp_z[i] * t_0_xxx_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_xxxy_0[i] += rpb_y[i] * t_0_xxx_0_xxxy_0[i] + rwp_y[i] * t_0_xxx_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxx_0_xxxx_0[i] += rpb_x[i] * t_0_xxx_0_xxxx_0[i] + rwp_x[i] * t_0_xxx_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1, t_0_xxx_0_xxx_1,\
                               t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1, t_0_xxx_0_xxxy_0,\
                               t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0, t_0_xxx_0_xxxz_1,\
                               t_0_xxxx_0_xxxx_0, t_0_xxxy_0_xxxy_0, t_0_xxxz_0_xxxz_0,\
                               t_0_xxyy_0_xxyy_0, t_0_xxyz_0_xxyz_0, t_0_xxz_0_xxyz_0,\
                               t_0_xxz_0_xxyz_1, t_0_xxz_0_xxz_1, t_0_xxzz_0_xxzz_0,\
                               t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1, t_0_xyy_0_xyy_1,\
                               t_0_xyyy_0_xyyy_0, t_0_xyyz_0_xyyz_0, t_0_xyzz_0_xyzz_0,\
                               t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1, t_0_xzz_0_xzz_1,\
                               t_0_xzzz_0_xzzz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_yyy_1, t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1,\
                               t_0_yyy_0_yyyz_0, t_0_yyy_0_yyyz_1, t_0_yyyy_0_yyyy_0,\
                               t_0_yyyz_0_yyyz_0, t_0_yyz_0_xyyz_0, t_0_yyz_0_xyyz_1,\
                               t_0_yyz_0_yyz_1, t_0_yyzz_0_yyzz_0, t_0_yzz_0_xyzz_0,\
                               t_0_yzz_0_xyzz_1, t_0_yzz_0_yyzz_0, t_0_yzz_0_yyzz_1,\
                               t_0_yzz_0_yzz_1, t_0_yzzz_0_yzzz_0, t_0_zz_0_xxzz_0,\
                               t_0_zz_0_xxzz_1, t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1, t_0_zz_0_zzzz_0,\
                               t_0_zz_0_zzzz_1, t_0_zzz_0_xzzz_0, t_0_zzz_0_xzzz_1,\
                               t_0_zzz_0_yzzz_0, t_0_zzz_0_yzzz_1, t_0_zzz_0_zzz_1,\
                               t_0_zzz_0_zzzz_0, t_0_zzz_0_zzzz_1, t_0_zzzz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzzz_0[i] = rpb_z[i] * t_0_zzz_0_zzzz_0[i] + rwp_z[i] * t_0_zzz_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yzzz_0[i] = rpb_y[i] * t_0_zzz_0_yzzz_0[i] + rwp_y[i] * t_0_zzz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yyzz_0_yyzz_0[i] = rpb_y[i] * t_0_yzz_0_yyzz_0[i] + rwp_y[i] * t_0_yzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_yyyz_0_yyyz_0[i] = rpb_z[i] * t_0_yyy_0_yyyz_0[i] + rwp_z[i] * t_0_yyy_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyy_0_yyyy_0[i] = rpb_y[i] * t_0_yyy_0_yyyy_0[i] + rwp_y[i] * t_0_yyy_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xzzz_0_xzzz_0[i] = rpb_x[i] * t_0_zzz_0_xzzz_0[i] + rwp_x[i] * t_0_zzz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_xyzz_0_xyzz_0[i] = rpb_x[i] * t_0_yzz_0_xyzz_0[i] + rwp_x[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyyz_0_xyyz_0[i] = rpb_x[i] * t_0_yyz_0_xyyz_0[i] + rwp_x[i] * t_0_yyz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyy_0_xyyy_0[i] = rpb_x[i] * t_0_yyy_0_xyyy_0[i] + rwp_x[i] * t_0_yyy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xxzz_0_xxzz_0[i] = rpb_x[i] * t_0_xzz_0_xxzz_0[i] + rwp_x[i] * t_0_xzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_xzz_0_xzz_1[i];

            t_0_xxyz_0_xxyz_0[i] = rpb_y[i] * t_0_xxz_0_xxyz_0[i] + rwp_y[i] * t_0_xxz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xxz_1[i];

            t_0_xxyy_0_xxyy_0[i] = rpb_x[i] * t_0_xyy_0_xxyy_0[i] + rwp_x[i] * t_0_xyy_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_xyy_0_xyy_1[i];

            t_0_xxxz_0_xxxz_0[i] = rpb_z[i] * t_0_xxx_0_xxxz_0[i] + rwp_z[i] * t_0_xxx_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_xxxy_0[i] = rpb_y[i] * t_0_xxx_0_xxxy_0[i] + rwp_y[i] * t_0_xxx_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxx_0_xxxx_0[i] = rpb_x[i] * t_0_xxx_0_xxxx_0[i] + rwp_x[i] * t_0_xxx_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];
        }
    }
}


} // derirec namespace
