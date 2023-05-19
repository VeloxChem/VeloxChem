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
compHostVRRForSPSD_V0(      BufferHostXY<T>&      intsBufferSPSD,
                      const BufferHostX<int32_t>& intsIndexesSPSD0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
                      const BufferHostXY<T>&      intsBufferSSSD0,
                      const BufferHostX<int32_t>& intsIndexesSSSD0,
                      const BufferHostXY<T>&      intsBufferSSSD1,
                      const BufferHostX<int32_t>& intsIndexesSSSD1,
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

    // set up [SPSD]^(0) integral components

    t_0_z_0_zz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(0));

    t_0_z_0_yz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(1));

    t_0_z_0_yy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(2));

    t_0_z_0_xz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(3));

    t_0_z_0_xy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(4));

    t_0_z_0_xx_0 = intsBufferSPSD0.data(intsIndexesSPSD0(5));

    t_0_y_0_zz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(6));

    t_0_y_0_yz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(7));

    t_0_y_0_yy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(8));

    t_0_y_0_xz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(9));

    t_0_y_0_xy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(10));

    t_0_y_0_xx_0 = intsBufferSPSD0.data(intsIndexesSPSD0(11));

    t_0_x_0_zz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(12));

    t_0_x_0_yz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(13));

    t_0_x_0_yy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(14));

    t_0_x_0_xz_0 = intsBufferSPSD0.data(intsIndexesSPSD0(15));

    t_0_x_0_xy_0 = intsBufferSPSD0.data(intsIndexesSPSD0(16));

    t_0_x_0_xx_0 = intsBufferSPSD0.data(intsIndexesSPSD0(17));

    // set up [SSSP]^(1) integral components

    t_0_0_0_z_1 = intsBufferSSSP1.data(intsIndexesSSSP1(0));

    t_0_0_0_y_1 = intsBufferSSSP1.data(intsIndexesSSSP1(1));

    t_0_0_0_x_1 = intsBufferSSSP1.data(intsIndexesSSSP1(2));

    // set up [SSSD]^(0) integral components

    t_0_0_0_zz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(0));

    t_0_0_0_yz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(1));

    t_0_0_0_yy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(2));

    t_0_0_0_xz_0 = intsBufferSSSD0.data(intsIndexesSSSD0(3));

    t_0_0_0_xy_0 = intsBufferSSSD0.data(intsIndexesSSSD0(4));

    t_0_0_0_xx_0 = intsBufferSSSD0.data(intsIndexesSSSD0(5));

    // set up [SSSD]^(1) integral components

    t_0_0_0_zz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(0));

    t_0_0_0_yz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(1));

    t_0_0_0_yy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(2));

    t_0_0_0_xz_1 = intsBufferSSSD1.data(intsIndexesSSSD1(3));

    t_0_0_0_xy_1 = intsBufferSSSD1.data(intsIndexesSSSD1(4));

    t_0_0_0_xx_1 = intsBufferSSSD1.data(intsIndexesSSSD1(5));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_x_1,\
                               t_0_0_0_xx_0, t_0_0_0_xx_1, t_0_0_0_xy_0, t_0_0_0_xy_1,\
                               t_0_0_0_xz_0, t_0_0_0_xz_1, t_0_0_0_y_1, t_0_0_0_yy_0,\
                               t_0_0_0_yy_1, t_0_0_0_yz_0, t_0_0_0_yz_1, t_0_0_0_z_1,\
                               t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_x_0_xx_0, t_0_x_0_xy_0,\
                               t_0_x_0_xz_0, t_0_x_0_yy_0, t_0_x_0_yz_0, t_0_x_0_zz_0,\
                               t_0_y_0_xx_0, t_0_y_0_xy_0, t_0_y_0_xz_0, t_0_y_0_yy_0,\
                               t_0_y_0_yz_0, t_0_y_0_zz_0, t_0_z_0_xx_0, t_0_z_0_xy_0,\
                               t_0_z_0_xz_0, t_0_z_0_yy_0, t_0_z_0_yz_0, t_0_z_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zz_0[i] += rpb_z[i] * t_0_0_0_zz_0[i] + rwp_z[i] * t_0_0_0_zz_1[i] + fze_0[i] * t_0_0_0_z_1[i];

            t_0_z_0_yz_0[i] += rpb_z[i] * t_0_0_0_yz_0[i] + rwp_z[i] * t_0_0_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_y_1[i];

            t_0_z_0_yy_0[i] += rpb_z[i] * t_0_0_0_yy_0[i] + rwp_z[i] * t_0_0_0_yy_1[i];

            t_0_z_0_xz_0[i] += rpb_z[i] * t_0_0_0_xz_0[i] + rwp_z[i] * t_0_0_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_x_1[i];

            t_0_z_0_xy_0[i] += rpb_z[i] * t_0_0_0_xy_0[i] + rwp_z[i] * t_0_0_0_xy_1[i];

            t_0_z_0_xx_0[i] += rpb_z[i] * t_0_0_0_xx_0[i] + rwp_z[i] * t_0_0_0_xx_1[i];

            t_0_y_0_zz_0[i] += rpb_y[i] * t_0_0_0_zz_0[i] + rwp_y[i] * t_0_0_0_zz_1[i];

            t_0_y_0_yz_0[i] += rpb_y[i] * t_0_0_0_yz_0[i] + rwp_y[i] * t_0_0_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_z_1[i];

            t_0_y_0_yy_0[i] += rpb_y[i] * t_0_0_0_yy_0[i] + rwp_y[i] * t_0_0_0_yy_1[i] + fze_0[i] * t_0_0_0_y_1[i];

            t_0_y_0_xz_0[i] += rpb_y[i] * t_0_0_0_xz_0[i] + rwp_y[i] * t_0_0_0_xz_1[i];

            t_0_y_0_xy_0[i] += rpb_y[i] * t_0_0_0_xy_0[i] + rwp_y[i] * t_0_0_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_x_1[i];

            t_0_y_0_xx_0[i] += rpb_y[i] * t_0_0_0_xx_0[i] + rwp_y[i] * t_0_0_0_xx_1[i];

            t_0_x_0_zz_0[i] += rpb_x[i] * t_0_0_0_zz_0[i] + rwp_x[i] * t_0_0_0_zz_1[i];

            t_0_x_0_yz_0[i] += rpb_x[i] * t_0_0_0_yz_0[i] + rwp_x[i] * t_0_0_0_yz_1[i];

            t_0_x_0_yy_0[i] += rpb_x[i] * t_0_0_0_yy_0[i] + rwp_x[i] * t_0_0_0_yy_1[i];

            t_0_x_0_xz_0[i] += rpb_x[i] * t_0_0_0_xz_0[i] + rwp_x[i] * t_0_0_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_z_1[i];

            t_0_x_0_xy_0[i] += rpb_x[i] * t_0_0_0_xy_0[i] + rwp_x[i] * t_0_0_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_y_1[i];

            t_0_x_0_xx_0[i] += rpb_x[i] * t_0_0_0_xx_0[i] + rwp_x[i] * t_0_0_0_xx_1[i] + fze_0[i] * t_0_0_0_x_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_x_1,\
                               t_0_0_0_xx_0, t_0_0_0_xx_1, t_0_0_0_xy_0, t_0_0_0_xy_1,\
                               t_0_0_0_xz_0, t_0_0_0_xz_1, t_0_0_0_y_1, t_0_0_0_yy_0,\
                               t_0_0_0_yy_1, t_0_0_0_yz_0, t_0_0_0_yz_1, t_0_0_0_z_1,\
                               t_0_0_0_zz_0, t_0_0_0_zz_1, t_0_x_0_xx_0, t_0_x_0_xy_0,\
                               t_0_x_0_xz_0, t_0_x_0_yy_0, t_0_x_0_yz_0, t_0_x_0_zz_0,\
                               t_0_y_0_xx_0, t_0_y_0_xy_0, t_0_y_0_xz_0, t_0_y_0_yy_0,\
                               t_0_y_0_yz_0, t_0_y_0_zz_0, t_0_z_0_xx_0, t_0_z_0_xy_0,\
                               t_0_z_0_xz_0, t_0_z_0_yy_0, t_0_z_0_yz_0, t_0_z_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_zz_0[i] = rpb_z[i] * t_0_0_0_zz_0[i] + rwp_z[i] * t_0_0_0_zz_1[i] + fze_0[i] * t_0_0_0_z_1[i];

            t_0_z_0_yz_0[i] = rpb_z[i] * t_0_0_0_yz_0[i] + rwp_z[i] * t_0_0_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_y_1[i];

            t_0_z_0_yy_0[i] = rpb_z[i] * t_0_0_0_yy_0[i] + rwp_z[i] * t_0_0_0_yy_1[i];

            t_0_z_0_xz_0[i] = rpb_z[i] * t_0_0_0_xz_0[i] + rwp_z[i] * t_0_0_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_x_1[i];

            t_0_z_0_xy_0[i] = rpb_z[i] * t_0_0_0_xy_0[i] + rwp_z[i] * t_0_0_0_xy_1[i];

            t_0_z_0_xx_0[i] = rpb_z[i] * t_0_0_0_xx_0[i] + rwp_z[i] * t_0_0_0_xx_1[i];

            t_0_y_0_zz_0[i] = rpb_y[i] * t_0_0_0_zz_0[i] + rwp_y[i] * t_0_0_0_zz_1[i];

            t_0_y_0_yz_0[i] = rpb_y[i] * t_0_0_0_yz_0[i] + rwp_y[i] * t_0_0_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_z_1[i];

            t_0_y_0_yy_0[i] = rpb_y[i] * t_0_0_0_yy_0[i] + rwp_y[i] * t_0_0_0_yy_1[i] + fze_0[i] * t_0_0_0_y_1[i];

            t_0_y_0_xz_0[i] = rpb_y[i] * t_0_0_0_xz_0[i] + rwp_y[i] * t_0_0_0_xz_1[i];

            t_0_y_0_xy_0[i] = rpb_y[i] * t_0_0_0_xy_0[i] + rwp_y[i] * t_0_0_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_x_1[i];

            t_0_y_0_xx_0[i] = rpb_y[i] * t_0_0_0_xx_0[i] + rwp_y[i] * t_0_0_0_xx_1[i];

            t_0_x_0_zz_0[i] = rpb_x[i] * t_0_0_0_zz_0[i] + rwp_x[i] * t_0_0_0_zz_1[i];

            t_0_x_0_yz_0[i] = rpb_x[i] * t_0_0_0_yz_0[i] + rwp_x[i] * t_0_0_0_yz_1[i];

            t_0_x_0_yy_0[i] = rpb_x[i] * t_0_0_0_yy_0[i] + rwp_x[i] * t_0_0_0_yy_1[i];

            t_0_x_0_xz_0[i] = rpb_x[i] * t_0_0_0_xz_0[i] + rwp_x[i] * t_0_0_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_z_1[i];

            t_0_x_0_xy_0[i] = rpb_x[i] * t_0_0_0_xy_0[i] + rwp_x[i] * t_0_0_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_y_1[i];

            t_0_x_0_xx_0[i] = rpb_x[i] * t_0_0_0_xx_0[i] + rwp_x[i] * t_0_0_0_xx_1[i] + fze_0[i] * t_0_0_0_x_1[i];
        }
    }
}


} // derirec namespace
