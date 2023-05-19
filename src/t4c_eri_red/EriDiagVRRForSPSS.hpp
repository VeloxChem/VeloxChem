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
compHostVRRForSPSS_V0(      BufferHostXY<T>&      intsBufferSPSS,
                      const BufferHostX<int32_t>& intsIndexesSPSS0,
                      const T*                    intsBufferSSSS0,
                      const T*                    intsBufferSSSS1,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SPSS]^(0) integral components

    t_0_z_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(0));

    t_0_y_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(1));

    t_0_x_0_0_0 = intsBufferSPSS0.data(intsIndexesSPSS0(2));

    // set up [SSSS]^(0) integral components

    t_0_0_0_0_0 = intsBufferSSSS0;

    // set up [SSSS]^(1) integral components

    t_0_0_0_0_1 = intsBufferSSSS1;

    if (useSummation)
    {
        #pragma omp simd align(rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_0_0,\
                               t_0_0_0_0_1, t_0_x_0_0_0, t_0_y_0_0_0, t_0_z_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_0_0[i] += rpb_z[i] * t_0_0_0_0_0[i] + rwp_z[i] * t_0_0_0_0_1[i];

            t_0_y_0_0_0[i] += rpb_y[i] * t_0_0_0_0_0[i] + rwp_y[i] * t_0_0_0_0_1[i];

            t_0_x_0_0_0[i] += rpb_x[i] * t_0_0_0_0_0[i] + rwp_x[i] * t_0_0_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z, t_0_0_0_0_0,\
                               t_0_0_0_0_1, t_0_x_0_0_0, t_0_y_0_0_0, t_0_z_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_z_0_0_0[i] = rpb_z[i] * t_0_0_0_0_0[i] + rwp_z[i] * t_0_0_0_0_1[i];

            t_0_y_0_0_0[i] = rpb_y[i] * t_0_0_0_0_0[i] + rwp_y[i] * t_0_0_0_0_1[i];

            t_0_x_0_0_0[i] = rpb_x[i] * t_0_0_0_0_0[i] + rwp_x[i] * t_0_0_0_0_1[i];
        }
    }
}


} // derirec namespace
