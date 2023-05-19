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
compHostHRRForPPPP_V0(      BufferHostXY<T>&      intsBufferPPPP,
                      const BufferHostX<int32_t>& intsIndexesPPPP,
                      const BufferHostXY<T>&      intsBufferSPPP,
                      const BufferHostX<int32_t>& intsIndexesSPPP,
                      const BufferHostXY<T>&      intsBufferSDPP,
                      const BufferHostX<int32_t>& intsIndexesSDPP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPPP) integral components

    t_z_z_z_z = intsBufferPPPP.data(intsIndexesPPPP(0));

    t_z_y_z_y = intsBufferPPPP.data(intsIndexesPPPP(1));

    t_z_x_z_x = intsBufferPPPP.data(intsIndexesPPPP(2));

    t_y_z_y_z = intsBufferPPPP.data(intsIndexesPPPP(3));

    t_y_y_y_y = intsBufferPPPP.data(intsIndexesPPPP(4));

    t_y_x_y_x = intsBufferPPPP.data(intsIndexesPPPP(5));

    t_x_z_x_z = intsBufferPPPP.data(intsIndexesPPPP(6));

    t_x_y_x_y = intsBufferPPPP.data(intsIndexesPPPP(7));

    t_x_x_x_x = intsBufferPPPP.data(intsIndexesPPPP(8));

    // set up (SPPP) integral components

    t_0_z_z_z = intsBufferSPPP.data(intsIndexesSPPP(0));

    t_0_z_y_z = intsBufferSPPP.data(intsIndexesSPPP(1));

    t_0_z_x_z = intsBufferSPPP.data(intsIndexesSPPP(2));

    t_0_y_z_y = intsBufferSPPP.data(intsIndexesSPPP(3));

    t_0_y_y_y = intsBufferSPPP.data(intsIndexesSPPP(4));

    t_0_y_x_y = intsBufferSPPP.data(intsIndexesSPPP(5));

    t_0_x_z_x = intsBufferSPPP.data(intsIndexesSPPP(6));

    t_0_x_y_x = intsBufferSPPP.data(intsIndexesSPPP(7));

    t_0_x_x_x = intsBufferSPPP.data(intsIndexesSPPP(8));

    // set up (SDPP) integral components

    t_0_zz_z_z = intsBufferSDPP.data(intsIndexesSDPP(0));

    t_0_yz_z_y = intsBufferSDPP.data(intsIndexesSDPP(1));

    t_0_yz_y_z = intsBufferSDPP.data(intsIndexesSDPP(2));

    t_0_yy_y_y = intsBufferSDPP.data(intsIndexesSDPP(3));

    t_0_xz_z_x = intsBufferSDPP.data(intsIndexesSDPP(4));

    t_0_xz_x_z = intsBufferSDPP.data(intsIndexesSDPP(5));

    t_0_xy_y_x = intsBufferSDPP.data(intsIndexesSDPP(6));

    t_0_xy_x_y = intsBufferSDPP.data(intsIndexesSDPP(7));

    t_0_xx_x_x = intsBufferSDPP.data(intsIndexesSDPP(8));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_x_x_x, t_0_x_y_x, t_0_x_z_x, t_0_xx_x_x,\
                           t_0_xy_x_y, t_0_xy_y_x, t_0_xz_x_z, t_0_xz_z_x, t_0_y_x_y,\
                           t_0_y_y_y, t_0_y_z_y, t_0_yy_y_y, t_0_yz_y_z, t_0_yz_z_y,\
                           t_0_z_x_z, t_0_z_y_z, t_0_z_z_z, t_0_zz_z_z, t_x_x_x_x, t_x_y_x_y,\
                           t_x_z_x_z, t_y_x_y_x, t_y_y_y_y, t_y_z_y_z, t_z_x_z_x, t_z_y_z_y,\
                           t_z_z_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_z_z[i] = t_0_zz_z_z[i] - rab_z[i] * t_0_z_z_z[i];

        t_z_y_z_y[i] = t_0_yz_z_y[i] - rab_z[i] * t_0_y_z_y[i];

        t_z_x_z_x[i] = t_0_xz_z_x[i] - rab_z[i] * t_0_x_z_x[i];

        t_y_z_y_z[i] = t_0_yz_y_z[i] - rab_y[i] * t_0_z_y_z[i];

        t_y_y_y_y[i] = t_0_yy_y_y[i] - rab_y[i] * t_0_y_y_y[i];

        t_y_x_y_x[i] = t_0_xy_y_x[i] - rab_y[i] * t_0_x_y_x[i];

        t_x_z_x_z[i] = t_0_xz_x_z[i] - rab_x[i] * t_0_z_x_z[i];

        t_x_y_x_y[i] = t_0_xy_x_y[i] - rab_x[i] * t_0_y_x_y[i];

        t_x_x_x_x[i] = t_0_xx_x_x[i] - rab_x[i] * t_0_x_x_x[i];
    }
}


} // derirec namespace
