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
compHostHRRForPPSS_V0(      BufferHostXY<T>&      intsBufferPPSS,
                      const BufferHostX<int32_t>& intsIndexesPPSS,
                      const BufferHostXY<T>&      intsBufferSPSS,
                      const BufferHostX<int32_t>& intsIndexesSPSS,
                      const BufferHostXY<T>&      intsBufferSDSS,
                      const BufferHostX<int32_t>& intsIndexesSDSS,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPSS) integral components

    t_z_z_0_0 = intsBufferPPSS.data(intsIndexesPPSS(0));

    t_z_y_0_0 = intsBufferPPSS.data(intsIndexesPPSS(1));

    t_z_x_0_0 = intsBufferPPSS.data(intsIndexesPPSS(2));

    t_y_z_0_0 = intsBufferPPSS.data(intsIndexesPPSS(3));

    t_y_y_0_0 = intsBufferPPSS.data(intsIndexesPPSS(4));

    t_y_x_0_0 = intsBufferPPSS.data(intsIndexesPPSS(5));

    t_x_z_0_0 = intsBufferPPSS.data(intsIndexesPPSS(6));

    t_x_y_0_0 = intsBufferPPSS.data(intsIndexesPPSS(7));

    t_x_x_0_0 = intsBufferPPSS.data(intsIndexesPPSS(8));

    // set up (SPSS) integral components

    t_0_z_0_0 = intsBufferSPSS.data(intsIndexesSPSS(0));

    t_0_y_0_0 = intsBufferSPSS.data(intsIndexesSPSS(1));

    t_0_x_0_0 = intsBufferSPSS.data(intsIndexesSPSS(2));

    // set up (SDSS) integral components

    t_0_zz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(0));

    t_0_yz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(1));

    t_0_yy_0_0 = intsBufferSDSS.data(intsIndexesSDSS(2));

    t_0_xz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(3));

    t_0_xy_0_0 = intsBufferSDSS.data(intsIndexesSDSS(4));

    t_0_xx_0_0 = intsBufferSDSS.data(intsIndexesSDSS(5));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_x_0_0, t_0_xx_0_0, t_0_xy_0_0, t_0_xz_0_0,\
                           t_0_y_0_0, t_0_yy_0_0, t_0_yz_0_0, t_0_z_0_0, t_0_zz_0_0,\
                           t_x_x_0_0, t_x_y_0_0, t_x_z_0_0, t_y_x_0_0, t_y_y_0_0, t_y_z_0_0,\
                           t_z_x_0_0, t_z_y_0_0, t_z_z_0_0 : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_0_0[i] = t_0_zz_0_0[i] - rab_z[i] * t_0_z_0_0[i];

        t_z_y_0_0[i] = t_0_yz_0_0[i] - rab_z[i] * t_0_y_0_0[i];

        t_z_x_0_0[i] = t_0_xz_0_0[i] - rab_z[i] * t_0_x_0_0[i];

        t_y_z_0_0[i] = t_0_yz_0_0[i] - rab_y[i] * t_0_z_0_0[i];

        t_y_y_0_0[i] = t_0_yy_0_0[i] - rab_y[i] * t_0_y_0_0[i];

        t_y_x_0_0[i] = t_0_xy_0_0[i] - rab_y[i] * t_0_x_0_0[i];

        t_x_z_0_0[i] = t_0_xz_0_0[i] - rab_x[i] * t_0_z_0_0[i];

        t_x_y_0_0[i] = t_0_xy_0_0[i] - rab_x[i] * t_0_y_0_0[i];

        t_x_x_0_0[i] = t_0_xx_0_0[i] - rab_x[i] * t_0_x_0_0[i];
    }
}


} // derirec namespace
