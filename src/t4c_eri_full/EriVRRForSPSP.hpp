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
compHostVRRForSPSP_V0(      BufferHostXY<T>&      intsBufferSPSP,
                      const BufferHostX<int32_t>& intsIndexesSPSP0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostXY<T>&      intsBufferSSSP0,
                      const BufferHostX<int32_t>& intsIndexesSSSP0,
                      const BufferHostXY<T>&      intsBufferSSSP1,
                      const BufferHostX<int32_t>& intsIndexesSSSP1,
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

    // set up [SPSP]^(0) integral components

    t_0_z_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(0));

    t_0_z_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(1));

    t_0_z_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(2));

    t_0_y_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(3));

    t_0_y_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(4));

    t_0_y_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(5));

    t_0_x_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(6));

    t_0_x_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(7));

    t_0_x_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(8));

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
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_0_1,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_x_0_x_0, t_0_x_0_y_0, t_0_x_0_z_0, t_0_y_0_x_0,\
                               t_0_y_0_y_0, t_0_y_0_z_0, t_0_z_0_x_0, t_0_z_0_y_0, t_0_z_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_z_0[i] += rpb_z[i] * t_0_0_0_z_0[i] + rwp_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];

            t_0_z_0_y_0[i] += rpb_z[i] * t_0_0_0_y_0[i] + rwp_z[i] * t_0_0_0_y_1[i];

            t_0_z_0_x_0[i] += rpb_z[i] * t_0_0_0_x_0[i] + rwp_z[i] * t_0_0_0_x_1[i];

            t_0_y_0_z_0[i] += rpb_y[i] * t_0_0_0_z_0[i] + rwp_y[i] * t_0_0_0_z_1[i];

            t_0_y_0_y_0[i] += rpb_y[i] * t_0_0_0_y_0[i] + rwp_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];

            t_0_y_0_x_0[i] += rpb_y[i] * t_0_0_0_x_0[i] + rwp_y[i] * t_0_0_0_x_1[i];

            t_0_x_0_z_0[i] += rpb_x[i] * t_0_0_0_z_0[i] + rwp_x[i] * t_0_0_0_z_1[i];

            t_0_x_0_y_0[i] += rpb_x[i] * t_0_0_0_y_0[i] + rwp_x[i] * t_0_0_0_y_1[i];

            t_0_x_0_x_0[i] += rpb_x[i] * t_0_0_0_x_0[i] + rwp_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_0_1,\
                               t_0_0_0_x_0, t_0_0_0_x_1, t_0_0_0_y_0, t_0_0_0_y_1, t_0_0_0_z_0,\
                               t_0_0_0_z_1, t_0_x_0_x_0, t_0_x_0_y_0, t_0_x_0_z_0, t_0_y_0_x_0,\
                               t_0_y_0_y_0, t_0_y_0_z_0, t_0_z_0_x_0, t_0_z_0_y_0, t_0_z_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_z_0[i] = rpb_z[i] * t_0_0_0_z_0[i] + rwp_z[i] * t_0_0_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];

            t_0_z_0_y_0[i] = rpb_z[i] * t_0_0_0_y_0[i] + rwp_z[i] * t_0_0_0_y_1[i];

            t_0_z_0_x_0[i] = rpb_z[i] * t_0_0_0_x_0[i] + rwp_z[i] * t_0_0_0_x_1[i];

            t_0_y_0_z_0[i] = rpb_y[i] * t_0_0_0_z_0[i] + rwp_y[i] * t_0_0_0_z_1[i];

            t_0_y_0_y_0[i] = rpb_y[i] * t_0_0_0_y_0[i] + rwp_y[i] * t_0_0_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];

            t_0_y_0_x_0[i] = rpb_y[i] * t_0_0_0_x_0[i] + rwp_y[i] * t_0_0_0_x_1[i];

            t_0_x_0_z_0[i] = rpb_x[i] * t_0_0_0_z_0[i] + rwp_x[i] * t_0_0_0_z_1[i];

            t_0_x_0_y_0[i] = rpb_x[i] * t_0_0_0_y_0[i] + rwp_x[i] * t_0_0_0_y_1[i];

            t_0_x_0_x_0[i] = rpb_x[i] * t_0_0_0_x_0[i] + rwp_x[i] * t_0_0_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_0_0_0_1[i];
        }
    }
}


} // derirec namespace
