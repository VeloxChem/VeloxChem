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
compHostVRRForSSSF_V0(      BufferHostXY<T>&      intsBufferSSSF,
                      const BufferHostX<int32_t>& intsIndexesSSSF0,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
                      const BufferHostXY<T>&      intsBufferSSSD0,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const BufferHostXY<T>&      intsBufferSSSD1,
                      const BufferHostX<int32_t>& intsIndexesSSSD1,
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

    // set up [SSSP]^(0) integral components

    t_0_0_0_z_0 = intsBufferSSSP0.data(intsIndexesSSSP0(0));

    t_0_0_0_y_0 = intsBufferSSSP0.data(intsIndexesSSSP0(1));

    t_0_0_0_x_0 = intsBufferSSSP0.data(intsIndexesSSSP0(2));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(3));

    // set up [SSSD]^(1) integral components

    t_0_0_0_zz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(0));

    t_0_0_0_yz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(1));

    t_0_0_0_yy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(2));

    t_0_0_0_xx_1 = intsBufferSSSD1.data(intsIndexesSSSD1(3));

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxy_0, t_0_0_0_xxz_0, t_0_0_0_xyy_0,\
                               t_0_0_0_xyz_0, t_0_0_0_xzz_0, t_0_0_0_y_0, t_0_0_0_y_1,\
                               t_0_0_0_yy_0, t_0_0_0_yy_1, t_0_0_0_yyy_0, t_0_0_0_yyz_0,\
                               t_0_0_0_yz_0, t_0_0_0_yz_1, t_0_0_0_yzz_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzz_0[i] += rqd_z[i] * t_0_0_0_zz_0[i] + rwq_z[i] * t_0_0_0_zz_1[i] + fe_0[i] * t_0_0_0_z_0[i] - fre2_0[i] * t_0_0_0_z_1[i];

            t_0_0_0_yzz_0[i] += rqd_y[i] * t_0_0_0_zz_0[i] + rwq_y[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyz_0[i] += rqd_z[i] * t_0_0_0_yy_0[i] + rwq_z[i] * t_0_0_0_yy_1[i];

            t_0_0_0_yyy_0[i] += rqd_y[i] * t_0_0_0_yy_0[i] + rwq_y[i] * t_0_0_0_yy_1[i] + fe_0[i] * t_0_0_0_y_0[i] - fre2_0[i] * t_0_0_0_y_1[i];

            t_0_0_0_xzz_0[i] += rqd_x[i] * t_0_0_0_zz_0[i] + rwq_x[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xyz_0[i] += rqd_x[i] * t_0_0_0_yz_0[i] + rwq_x[i] * t_0_0_0_yz_1[i];

            t_0_0_0_xyy_0[i] += rqd_x[i] * t_0_0_0_yy_0[i] + rwq_x[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxz_0[i] += rqd_z[i] * t_0_0_0_xx_0[i] + rwq_z[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxy_0[i] += rqd_y[i] * t_0_0_0_xx_0[i] + rwq_y[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxx_0[i] += rqd_x[i] * t_0_0_0_xx_0[i] + rwq_x[i] * t_0_0_0_xx_1[i] + fe_0[i] * t_0_0_0_x_0[i] - fre2_0[i] * t_0_0_0_x_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxy_0, t_0_0_0_xxz_0, t_0_0_0_xyy_0,\
                               t_0_0_0_xyz_0, t_0_0_0_xzz_0, t_0_0_0_y_0, t_0_0_0_y_1,\
                               t_0_0_0_yy_0, t_0_0_0_yy_1, t_0_0_0_yyy_0, t_0_0_0_yyz_0,\
                               t_0_0_0_yz_0, t_0_0_0_yz_1, t_0_0_0_yzz_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzz_0[i] = rqd_z[i] * t_0_0_0_zz_0[i] + rwq_z[i] * t_0_0_0_zz_1[i] + fe_0[i] * t_0_0_0_z_0[i] - fre2_0[i] * t_0_0_0_z_1[i];

            t_0_0_0_yzz_0[i] = rqd_y[i] * t_0_0_0_zz_0[i] + rwq_y[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyz_0[i] = rqd_z[i] * t_0_0_0_yy_0[i] + rwq_z[i] * t_0_0_0_yy_1[i];

            t_0_0_0_yyy_0[i] = rqd_y[i] * t_0_0_0_yy_0[i] + rwq_y[i] * t_0_0_0_yy_1[i] + fe_0[i] * t_0_0_0_y_0[i] - fre2_0[i] * t_0_0_0_y_1[i];

            t_0_0_0_xzz_0[i] = rqd_x[i] * t_0_0_0_zz_0[i] + rwq_x[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xyz_0[i] = rqd_x[i] * t_0_0_0_yz_0[i] + rwq_x[i] * t_0_0_0_yz_1[i];

            t_0_0_0_xyy_0[i] = rqd_x[i] * t_0_0_0_yy_0[i] + rwq_x[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxz_0[i] = rqd_z[i] * t_0_0_0_xx_0[i] + rwq_z[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxy_0[i] = rqd_y[i] * t_0_0_0_xx_0[i] + rwq_y[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxx_0[i] = rqd_x[i] * t_0_0_0_xx_0[i] + rwq_x[i] * t_0_0_0_xx_1[i] + fe_0[i] * t_0_0_0_x_0[i] - fre2_0[i] * t_0_0_0_x_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSSSF_V1(      BufferHostXY<T>&      intsBufferSSSF,
                      const BufferHostX<int32_t>& intsIndexesSSSF0,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
                      const BufferHostXY<T>&      intsBufferSSSD0,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const BufferHostXY<T>&      intsBufferSSSD1,
                      const BufferHostX<int32_t>& intsIndexesSSSD1,
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

    // set up [SSSF]^(0) integral components

    t_0_0_0_zzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(0));

    t_0_0_0_yzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(1));

    t_0_0_0_yyz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(2));

    t_0_0_0_yyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(3));

    t_0_0_0_xzz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(4));

    t_0_0_0_xyy_0 = intsBufferSSSF0.data(intsIndexesSSSF0(5));

    t_0_0_0_xxz_0 = intsBufferSSSF0.data(intsIndexesSSSF0(6));

    t_0_0_0_xxx_0 = intsBufferSSSF0.data(intsIndexesSSSF0(7));

    // set up [SSSP]^(0) integral components

    t_0_0_0_z_0 = intsBufferSSSP0.data(intsIndexesSSSP0(0));

    t_0_0_0_y_0 = intsBufferSSSP0.data(intsIndexesSSSP0(1));

    t_0_0_0_x_0 = intsBufferSSSP0.data(intsIndexesSSSP0(2));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    // set up [SSSD]^(1) integral components

    t_0_0_0_zz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(0));

    t_0_0_0_yy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(1));

    t_0_0_0_xx_1 = intsBufferSSSD1.data(intsIndexesSSSD1(2));

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxz_0, t_0_0_0_xyy_0, t_0_0_0_xzz_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_yy_1,\
                               t_0_0_0_yyy_0, t_0_0_0_yyz_0, t_0_0_0_yzz_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzz_0[i] += rqd_z[i] * t_0_0_0_zz_0[i] + rwq_z[i] * t_0_0_0_zz_1[i] + fe_0[i] * t_0_0_0_z_0[i] - fre2_0[i] * t_0_0_0_z_1[i];

            t_0_0_0_yzz_0[i] += rqd_y[i] * t_0_0_0_zz_0[i] + rwq_y[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyz_0[i] += rqd_z[i] * t_0_0_0_yy_0[i] + rwq_z[i] * t_0_0_0_yy_1[i];

            t_0_0_0_yyy_0[i] += rqd_y[i] * t_0_0_0_yy_0[i] + rwq_y[i] * t_0_0_0_yy_1[i] + fe_0[i] * t_0_0_0_y_0[i] - fre2_0[i] * t_0_0_0_y_1[i];

            t_0_0_0_xzz_0[i] += rqd_x[i] * t_0_0_0_zz_0[i] + rwq_x[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xyy_0[i] += rqd_x[i] * t_0_0_0_yy_0[i] + rwq_x[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxz_0[i] += rqd_z[i] * t_0_0_0_xx_0[i] + rwq_z[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxx_0[i] += rqd_x[i] * t_0_0_0_xx_0[i] + rwq_x[i] * t_0_0_0_xx_1[i] + fe_0[i] * t_0_0_0_x_0[i] - fre2_0[i] * t_0_0_0_x_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0, t_0_0_0_xx_1,\
                               t_0_0_0_xxx_0, t_0_0_0_xxz_0, t_0_0_0_xyy_0, t_0_0_0_xzz_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_yy_1,\
                               t_0_0_0_yyy_0, t_0_0_0_yyz_0, t_0_0_0_yzz_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_0_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zzz_0[i] = rqd_z[i] * t_0_0_0_zz_0[i] + rwq_z[i] * t_0_0_0_zz_1[i] + fe_0[i] * t_0_0_0_z_0[i] - fre2_0[i] * t_0_0_0_z_1[i];

            t_0_0_0_yzz_0[i] = rqd_y[i] * t_0_0_0_zz_0[i] + rwq_y[i] * t_0_0_0_zz_1[i];

            t_0_0_0_yyz_0[i] = rqd_z[i] * t_0_0_0_yy_0[i] + rwq_z[i] * t_0_0_0_yy_1[i];

            t_0_0_0_yyy_0[i] = rqd_y[i] * t_0_0_0_yy_0[i] + rwq_y[i] * t_0_0_0_yy_1[i] + fe_0[i] * t_0_0_0_y_0[i] - fre2_0[i] * t_0_0_0_y_1[i];

            t_0_0_0_xzz_0[i] = rqd_x[i] * t_0_0_0_zz_0[i] + rwq_x[i] * t_0_0_0_zz_1[i];

            t_0_0_0_xyy_0[i] = rqd_x[i] * t_0_0_0_yy_0[i] + rwq_x[i] * t_0_0_0_yy_1[i];

            t_0_0_0_xxz_0[i] = rqd_z[i] * t_0_0_0_xx_0[i] + rwq_z[i] * t_0_0_0_xx_1[i];

            t_0_0_0_xxx_0[i] = rqd_x[i] * t_0_0_0_xx_0[i] + rwq_x[i] * t_0_0_0_xx_1[i] + fe_0[i] * t_0_0_0_x_0[i] - fre2_0[i] * t_0_0_0_x_1[i];
        }
    }
}


} // derirec namespace
