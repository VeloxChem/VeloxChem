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
compHostVRRForSSSD_V0(      BufferHostXY<T>&      intsBufferSSSD,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
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

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    t_0_0_0_xz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(3));

    t_0_0_0_xy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(4));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(5));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SSSP]^(0) integral components

    t_0_0_0_z_0 = intsBufferSSSP0.data(intsIndexesSSSP0(0));

    t_0_0_0_y_0 = intsBufferSSSP0.data(intsIndexesSSSP0(1));

    t_0_0_0_x_0 = intsBufferSSSP0.data(intsIndexesSSSP0(2));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_xy_0, t_0_0_0_xz_0, t_0_0_0_y_0, t_0_0_0_y_1,\
                               t_0_0_0_yy_0, t_0_0_0_yz_0, t_0_0_0_z_0, t_0_0_0_z_1,\
                               t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] += rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yz_0[i] += rqd_y[i] * t_0_0_0_z_0[i] + rwq_y[i] * t_0_0_0_z_1[i];

            t_0_0_0_yy_0[i] += rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xz_0[i] += rqd_x[i] * t_0_0_0_z_0[i] + rwq_x[i] * t_0_0_0_z_1[i];

            t_0_0_0_xy_0[i] += rqd_x[i] * t_0_0_0_y_0[i] + rwq_x[i] * t_0_0_0_y_1[i];

            t_0_0_0_xx_0[i] += rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_xy_0, t_0_0_0_xz_0, t_0_0_0_y_0, t_0_0_0_y_1,\
                               t_0_0_0_yy_0, t_0_0_0_yz_0, t_0_0_0_z_0, t_0_0_0_z_1,\
                               t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] = rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yz_0[i] = rqd_y[i] * t_0_0_0_z_0[i] + rwq_y[i] * t_0_0_0_z_1[i];

            t_0_0_0_yy_0[i] = rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xz_0[i] = rqd_x[i] * t_0_0_0_z_0[i] + rwq_x[i] * t_0_0_0_z_1[i];

            t_0_0_0_xy_0[i] = rqd_x[i] * t_0_0_0_y_0[i] + rwq_x[i] * t_0_0_0_y_1[i];

            t_0_0_0_xx_0[i] = rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSSSD_V1(      BufferHostXY<T>&      intsBufferSSSD,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
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

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(3));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SSSP]^(0) integral components

    t_0_0_0_z_0 = intsBufferSSSP0.data(intsIndexesSSSP0(0));

    t_0_0_0_y_0 = intsBufferSSSP0.data(intsIndexesSSSP0(1));

    t_0_0_0_x_0 = intsBufferSSSP0.data(intsIndexesSSSP0(2));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_yz_0,\
                               t_0_0_0_z_0, t_0_0_0_z_1, t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] += rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yz_0[i] += rqd_y[i] * t_0_0_0_z_0[i] + rwq_y[i] * t_0_0_0_z_1[i];

            t_0_0_0_yy_0[i] += rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xx_0[i] += rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_yz_0,\
                               t_0_0_0_z_0, t_0_0_0_z_1, t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] = rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yz_0[i] = rqd_y[i] * t_0_0_0_z_0[i] + rwq_y[i] * t_0_0_0_z_1[i];

            t_0_0_0_yy_0[i] = rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xx_0[i] = rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSSSD_V2(      BufferHostXY<T>&      intsBufferSSSD,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
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

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    // set up [SSSP]^(0) integral components

    t_0_0_0_z_0 = intsBufferSSSP0.data(intsIndexesSSSP0(0));

    t_0_0_0_y_0 = intsBufferSSSP0.data(intsIndexesSSSP0(1));

    t_0_0_0_x_0 = intsBufferSSSP0.data(intsIndexesSSSP0(2));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] += rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yy_0[i] += rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xx_0[i] += rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fe_0, fre2_0, rqd_x, rqd_y, rqd_z, rwq_x, rwq_y, rwq_z,\
                               t_0_0_0_0_0, t_0_0_0_0_1, t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_xx_0,\
                               t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_yy_0, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_0_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_0_0_zz_0[i] = rqd_z[i] * t_0_0_0_z_0[i] + rwq_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_yy_0[i] = rqd_y[i] * t_0_0_0_y_0[i] + rwq_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];

            t_0_0_0_xx_0[i] = rqd_x[i] * t_0_0_0_x_0[i] + rwq_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fe_0[i] * t_0_0_0_0_0[i] - fact_1_2 * fre2_0[i] * t_0_0_0_0_1[i];
        }
    }
}


} // derirec namespace
