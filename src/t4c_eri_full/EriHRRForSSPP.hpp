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
compHostHRRForSSPP_V0(      BufferHostXY<T>&      intsBufferSSPP,
                      const BufferHostX<int32_t>& intsIndexesSSPP,
                      const BufferHostXY<T>&      intsBufferSSSP,
                      const BufferHostX<int32_t>& intsIndexesSSSP,
                      const BufferHostXY<T>&      intsBufferSSSD,
                      const BufferHostX<int32_t>& intsIndexesSSSD,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SSPP) integral components

    t_0_0_z_z = intsBufferSSPP.data(intsIndexesSSPP(0));

    t_0_0_z_y = intsBufferSSPP.data(intsIndexesSSPP(1));

    t_0_0_z_x = intsBufferSSPP.data(intsIndexesSSPP(2));

    t_0_0_y_z = intsBufferSSPP.data(intsIndexesSSPP(3));

    t_0_0_y_y = intsBufferSSPP.data(intsIndexesSSPP(4));

    t_0_0_y_x = intsBufferSSPP.data(intsIndexesSSPP(5));

    t_0_0_x_z = intsBufferSSPP.data(intsIndexesSSPP(6));

    t_0_0_x_y = intsBufferSSPP.data(intsIndexesSSPP(7));

    t_0_0_x_x = intsBufferSSPP.data(intsIndexesSSPP(8));

    // set up (SSSP) integral components

    t_0_0_0_z = intsBufferSSSP.data(intsIndexesSSSP(0));

    t_0_0_0_y = intsBufferSSSP.data(intsIndexesSSSP(1));

    t_0_0_0_x = intsBufferSSSP.data(intsIndexesSSSP(2));

    // set up (SSSD) integral components

    t_0_0_0_zz = intsBufferSSSD.data(intsIndexesSSSD(0));

    t_0_0_0_yz = intsBufferSSSD.data(intsIndexesSSSD(1));

    t_0_0_0_yy = intsBufferSSSD.data(intsIndexesSSSD(2));

    t_0_0_0_xz = intsBufferSSSD.data(intsIndexesSSSD(3));

    t_0_0_0_xy = intsBufferSSSD.data(intsIndexesSSSD(4));

    t_0_0_0_xx = intsBufferSSSD.data(intsIndexesSSSD(5));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_0_0_x, t_0_0_0_xx, t_0_0_0_xy, t_0_0_0_xz,\
                           t_0_0_0_y, t_0_0_0_yy, t_0_0_0_yz, t_0_0_0_z, t_0_0_0_zz,\
                           t_0_0_x_x, t_0_0_x_y, t_0_0_x_z, t_0_0_y_x, t_0_0_y_y, t_0_0_y_z,\
                           t_0_0_z_x, t_0_0_z_y, t_0_0_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_0_z_z[i] = t_0_0_0_zz[i] - rcd_z[i] * t_0_0_0_z[i];

        t_0_0_z_y[i] = t_0_0_0_yz[i] - rcd_z[i] * t_0_0_0_y[i];

        t_0_0_z_x[i] = t_0_0_0_xz[i] - rcd_z[i] * t_0_0_0_x[i];

        t_0_0_y_z[i] = t_0_0_0_yz[i] - rcd_y[i] * t_0_0_0_z[i];

        t_0_0_y_y[i] = t_0_0_0_yy[i] - rcd_y[i] * t_0_0_0_y[i];

        t_0_0_y_x[i] = t_0_0_0_xy[i] - rcd_y[i] * t_0_0_0_x[i];

        t_0_0_x_z[i] = t_0_0_0_xz[i] - rcd_x[i] * t_0_0_0_z[i];

        t_0_0_x_y[i] = t_0_0_0_xy[i] - rcd_x[i] * t_0_0_0_y[i];

        t_0_0_x_x[i] = t_0_0_0_xx[i] - rcd_x[i] * t_0_0_0_x[i];
    }
}


} // derirec namespace
