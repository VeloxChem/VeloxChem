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
compHostVRRForSSSG_V0(      BufferHostXY<T>&      intsBufferSSSG,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSD0,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const BufferHostXY<T>&      intsBufferSSSD1,
                      const BufferHostX<int32_t>& intsIndexesSSSD1,
                      const BufferHostXY<T>&      intsBufferSSSF0,
                      const BufferHostX<int32_t>& intsIndexesSSSF0,
                      const BufferHostXY<T>&      intsBufferSSSF1,
                      const BufferHostX<int32_t>& intsIndexesSSSF1,
                      const T*                    osFactorsKetZeta,
                      const BufferHostMY<T, 3>&   rDistancesQD,
                      const BufferHostMY<T, 3>&   rDistancesWQ,
                      const T*                    osFactorsKetRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fe_0 = osFactorsKetZeta;

    auto fre2_0 = osFactorsKetRhoZeta;

    // set up R(QD) distances

    auto rqd_z = rDistancesQD.data(2);

    auto rqd_y = rDistancesQD.data(1);

    auto rqd_x = rDistancesQD.data(0);

    // set up R(WQ) distances

    auto rwq_z = rDistancesWQ.data(2);

    auto rwq_y = rDistancesWQ.data(1);

    auto rwq_x = rDistancesWQ.data(0);

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

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    // set up [SSSD]^(1) integral components

    t_0_0_0_zz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(0));

    t_0_0_0_yy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(1));

    t_0_0_0_xx_1 = intsBufferSSSD1.data(intsIndexesSSSD1(2));

    // set up [SSSF]^(0) integral components

    t_0_0_0_zzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(0));

    t_0_0_0_yzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(1));

    t_0_0_0_yyz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(2));

    t_0_0_0_yyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(3));

    t_0_0_0_xzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(4));

    t_0_0_0_xyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(5));

    t_0_0_0_xxz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(6));

    t_0_0_0_xxx_0 = intsBufferSSSF0.data(intsIndexesSSSF0(7));

    // set up [SSSF]^(1) integral components

    t_0_0_0_zzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(0));

    t_0_0_0_yzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(1));

    t_0_0_0_yyz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(2));

    t_0_0_0_yyy_1 = intsBufferSSSF1.data(intsIndexesSSSF1(3));

    t_0_0_0_xzz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(4));

    t_0_0_0_xyy_1 = intsBufferSSSF1.data(intsIndexesSSSF1(5));

    t_0_0_0_xxz_1 = intsBufferSSSF1.data(intsIndexesSSSF1(6));

    t_0_0_0_xxx_1 = intsBufferSSSF1.data(intsIndexesSSSF1(7));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_xx_0, t_0_0_0_xx_1, t_0_0_0_xxx_0, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxy_0, t_0_0_0_xxxz_0, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyz_0, t_0_0_0_xxz_0, t_0_0_0_xxz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xyy_0, t_0_0_0_xyy_1, t_0_0_0_xyyy_0, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyzz_0, t_0_0_0_xzz_0, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_yy_0, t_0_0_0_yy_1, t_0_0_0_yyy_0, t_0_0_0_yyy_1,\
                               t_0_0_0_yyyy_0, t_0_0_0_yyyz_0, t_0_0_0_yyz_0, t_0_0_0_yyz_1,\
                               t_0_0_0_yyzz_0, t_0_0_0_yzz_0, t_0_0_0_yzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0, t_0_0_0_zzz_1,\
                               t_0_0_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzzz_0[i] += rqd_z[i] * t_0_0_0_zzz_0[i] + rwq_z[i] * t_0_0_0_zzz_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yzzz_0[i] += rqd_y[i] * t_0_0_0_zzz_0[i] + rwq_y[i] * t_0_0_0_zzz_1[i];

            t_0_0_0_yyzz_0[i] += rqd_y[i] * t_0_0_0_yzz_0[i] + rwq_y[i] * t_0_0_0_yzz_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyyz_0[i] += rqd_z[i] * t_0_0_0_yyy_0[i] + rwq_z[i] * t_0_0_0_yyy_1[i];

            t_0_0_0_yyyy_0[i] += rqd_y[i] * t_0_0_0_yyy_0[i] + rwq_y[i] * t_0_0_0_yyy_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_yy_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xzzz_0[i] += rqd_x[i] * t_0_0_0_zzz_0[i] + rwq_x[i] * t_0_0_0_zzz_1[i];

            t_0_0_0_xyzz_0[i] += rqd_x[i] * t_0_0_0_yzz_0[i] + rwq_x[i] * t_0_0_0_yzz_1[i];

            t_0_0_0_xyyz_0[i] += rqd_x[i] * t_0_0_0_yyz_0[i] + rwq_x[i] * t_0_0_0_yyz_1[i];

            t_0_0_0_xyyy_0[i] += rqd_x[i] * t_0_0_0_yyy_0[i] + rwq_x[i] * t_0_0_0_yyy_1[i];

            t_0_0_0_xxzz_0[i] += rqd_x[i] * t_0_0_0_xzz_0[i] + rwq_x[i] * t_0_0_0_xzz_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xxyz_0[i] += rqd_y[i] * t_0_0_0_xxz_0[i] + rwq_y[i] * t_0_0_0_xxz_1[i];

            t_0_0_0_xxyy_0[i] += rqd_x[i] * t_0_0_0_xyy_0[i] + rwq_x[i] * t_0_0_0_xyy_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_yy_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxxz_0[i] += rqd_z[i] * t_0_0_0_xxx_0[i] + rwq_z[i] * t_0_0_0_xxx_1[i];

            t_0_0_0_xxxy_0[i] += rqd_y[i] * t_0_0_0_xxx_0[i] + rwq_y[i] * t_0_0_0_xxx_1[i];

            t_0_0_0_xxxx_0[i] += rqd_x[i] * t_0_0_0_xxx_0[i] + rwq_x[i] * t_0_0_0_xxx_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_xx_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_xx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_xx_0, t_0_0_0_xx_1, t_0_0_0_xxx_0, t_0_0_0_xxx_1,\
                               t_0_0_0_xxxx_0, t_0_0_0_xxxy_0, t_0_0_0_xxxz_0, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyz_0, t_0_0_0_xxz_0, t_0_0_0_xxz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xyy_0, t_0_0_0_xyy_1, t_0_0_0_xyyy_0, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyzz_0, t_0_0_0_xzz_0, t_0_0_0_xzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_yy_0, t_0_0_0_yy_1, t_0_0_0_yyy_0, t_0_0_0_yyy_1,\
                               t_0_0_0_yyyy_0, t_0_0_0_yyyz_0, t_0_0_0_yyz_0, t_0_0_0_yyz_1,\
                               t_0_0_0_yyzz_0, t_0_0_0_yzz_0, t_0_0_0_yzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0, t_0_0_0_zzz_1,\
                               t_0_0_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzzz_0[i] = rqd_z[i] * t_0_0_0_zzz_0[i] + rwq_z[i] * t_0_0_0_zzz_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yzzz_0[i] = rqd_y[i] * t_0_0_0_zzz_0[i] + rwq_y[i] * t_0_0_0_zzz_1[i];

            t_0_0_0_yyzz_0[i] = rqd_y[i] * t_0_0_0_yzz_0[i] + rwq_y[i] * t_0_0_0_yzz_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyyz_0[i] = rqd_z[i] * t_0_0_0_yyy_0[i] + rwq_z[i] * t_0_0_0_yyy_1[i];

            t_0_0_0_yyyy_0[i] = rqd_y[i] * t_0_0_0_yyy_0[i] + rwq_y[i] * t_0_0_0_yyy_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_yy_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xzzz_0[i] = rqd_x[i] * t_0_0_0_zzz_0[i] + rwq_x[i] * t_0_0_0_zzz_1[i];

            t_0_0_0_xyzz_0[i] = rqd_x[i] * t_0_0_0_yzz_0[i] + rwq_x[i] * t_0_0_0_yzz_1[i];

            t_0_0_0_xyyz_0[i] = rqd_x[i] * t_0_0_0_yyz_0[i] + rwq_x[i] * t_0_0_0_yyz_1[i];

            t_0_0_0_xyyy_0[i] = rqd_x[i] * t_0_0_0_yyy_0[i] + rwq_x[i] * t_0_0_0_yyy_1[i];

            t_0_0_0_xxzz_0[i] = rqd_x[i] * t_0_0_0_xzz_0[i] + rwq_x[i] * t_0_0_0_xzz_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_zz_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xxyz_0[i] = rqd_y[i] * t_0_0_0_xxz_0[i] + rwq_y[i] * t_0_0_0_xxz_1[i];

            t_0_0_0_xxyy_0[i] = rqd_x[i] * t_0_0_0_xyy_0[i] + rwq_x[i] * t_0_0_0_xyy_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_yy_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxxz_0[i] = rqd_z[i] * t_0_0_0_xxx_0[i] + rwq_z[i] * t_0_0_0_xxx_1[i];

            t_0_0_0_xxxy_0[i] = rqd_y[i] * t_0_0_0_xxx_0[i] + rwq_y[i] * t_0_0_0_xxx_1[i];

            t_0_0_0_xxxx_0[i] = rqd_x[i] * t_0_0_0_xxx_0[i] + rwq_x[i] * t_0_0_0_xxx_1[i] + fact_3_2 * fe_0[i] * t_0_0_0_xx_0[i] - fact_3_2 * fre2_0[i] * t_0_0_0_xx_1[i];
        }
    }
}


} // derirec namespace
