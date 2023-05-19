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
compHostHRRForPPSP_V0(      BufferHostXY<T>&      intsBufferPPSP,
                      const BufferHostX<int32_t>& intsIndexesPPSP,
                      const BufferHostXY<T>&      intsBufferSPSP,
                      const BufferHostX<int32_t>& intsIndexesSPSP,
                      const BufferHostXY<T>&      intsBufferSDSP,
                      const BufferHostX<int32_t>& intsIndexesSDSP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPSP) integral components

    t_z_z_0_z = intsBufferPPSP.data(intsIndexesPPSP(0));

    t_z_z_0_y = intsBufferPPSP.data(intsIndexesPPSP(1));

    t_z_z_0_x = intsBufferPPSP.data(intsIndexesPPSP(2));

    t_z_y_0_z = intsBufferPPSP.data(intsIndexesPPSP(3));

    t_z_y_0_y = intsBufferPPSP.data(intsIndexesPPSP(4));

    t_z_y_0_x = intsBufferPPSP.data(intsIndexesPPSP(5));

    t_z_x_0_z = intsBufferPPSP.data(intsIndexesPPSP(6));

    t_z_x_0_y = intsBufferPPSP.data(intsIndexesPPSP(7));

    t_z_x_0_x = intsBufferPPSP.data(intsIndexesPPSP(8));

    t_y_z_0_z = intsBufferPPSP.data(intsIndexesPPSP(9));

    t_y_z_0_y = intsBufferPPSP.data(intsIndexesPPSP(10));

    t_y_z_0_x = intsBufferPPSP.data(intsIndexesPPSP(11));

    t_y_y_0_z = intsBufferPPSP.data(intsIndexesPPSP(12));

    t_y_y_0_y = intsBufferPPSP.data(intsIndexesPPSP(13));

    t_y_y_0_x = intsBufferPPSP.data(intsIndexesPPSP(14));

    t_y_x_0_z = intsBufferPPSP.data(intsIndexesPPSP(15));

    t_y_x_0_y = intsBufferPPSP.data(intsIndexesPPSP(16));

    t_y_x_0_x = intsBufferPPSP.data(intsIndexesPPSP(17));

    t_x_z_0_z = intsBufferPPSP.data(intsIndexesPPSP(18));

    t_x_z_0_y = intsBufferPPSP.data(intsIndexesPPSP(19));

    t_x_z_0_x = intsBufferPPSP.data(intsIndexesPPSP(20));

    t_x_y_0_z = intsBufferPPSP.data(intsIndexesPPSP(21));

    t_x_y_0_y = intsBufferPPSP.data(intsIndexesPPSP(22));

    t_x_y_0_x = intsBufferPPSP.data(intsIndexesPPSP(23));

    t_x_x_0_z = intsBufferPPSP.data(intsIndexesPPSP(24));

    t_x_x_0_y = intsBufferPPSP.data(intsIndexesPPSP(25));

    t_x_x_0_x = intsBufferPPSP.data(intsIndexesPPSP(26));

    // set up (SPSP) integral components

    t_0_z_0_z = intsBufferSPSP.data(intsIndexesSPSP(0));

    t_0_z_0_y = intsBufferSPSP.data(intsIndexesSPSP(1));

    t_0_z_0_x = intsBufferSPSP.data(intsIndexesSPSP(2));

    t_0_y_0_z = intsBufferSPSP.data(intsIndexesSPSP(3));

    t_0_y_0_y = intsBufferSPSP.data(intsIndexesSPSP(4));

    t_0_y_0_x = intsBufferSPSP.data(intsIndexesSPSP(5));

    t_0_x_0_z = intsBufferSPSP.data(intsIndexesSPSP(6));

    t_0_x_0_y = intsBufferSPSP.data(intsIndexesSPSP(7));

    t_0_x_0_x = intsBufferSPSP.data(intsIndexesSPSP(8));

    // set up (SDSP) integral components

    t_0_zz_0_z = intsBufferSDSP.data(intsIndexesSDSP(0));

    t_0_zz_0_y = intsBufferSDSP.data(intsIndexesSDSP(1));

    t_0_zz_0_x = intsBufferSDSP.data(intsIndexesSDSP(2));

    t_0_yz_0_z = intsBufferSDSP.data(intsIndexesSDSP(3));

    t_0_yz_0_y = intsBufferSDSP.data(intsIndexesSDSP(4));

    t_0_yz_0_x = intsBufferSDSP.data(intsIndexesSDSP(5));

    t_0_yy_0_z = intsBufferSDSP.data(intsIndexesSDSP(6));

    t_0_yy_0_y = intsBufferSDSP.data(intsIndexesSDSP(7));

    t_0_yy_0_x = intsBufferSDSP.data(intsIndexesSDSP(8));

    t_0_xz_0_z = intsBufferSDSP.data(intsIndexesSDSP(9));

    t_0_xz_0_y = intsBufferSDSP.data(intsIndexesSDSP(10));

    t_0_xz_0_x = intsBufferSDSP.data(intsIndexesSDSP(11));

    t_0_xy_0_z = intsBufferSDSP.data(intsIndexesSDSP(12));

    t_0_xy_0_y = intsBufferSDSP.data(intsIndexesSDSP(13));

    t_0_xy_0_x = intsBufferSDSP.data(intsIndexesSDSP(14));

    t_0_xx_0_z = intsBufferSDSP.data(intsIndexesSDSP(15));

    t_0_xx_0_y = intsBufferSDSP.data(intsIndexesSDSP(16));

    t_0_xx_0_x = intsBufferSDSP.data(intsIndexesSDSP(17));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_x_0_x, t_0_x_0_y, t_0_x_0_z, t_0_xx_0_x,\
                           t_0_xx_0_y, t_0_xx_0_z, t_0_xy_0_x, t_0_xy_0_y, t_0_xy_0_z,\
                           t_0_xz_0_x, t_0_xz_0_y, t_0_xz_0_z, t_0_y_0_x, t_0_y_0_y,\
                           t_0_y_0_z, t_0_yy_0_x, t_0_yy_0_y, t_0_yy_0_z, t_0_yz_0_x,\
                           t_0_yz_0_y, t_0_yz_0_z, t_0_z_0_x, t_0_z_0_y, t_0_z_0_z,\
                           t_0_zz_0_x, t_0_zz_0_y, t_0_zz_0_z, t_x_x_0_x, t_x_x_0_y,\
                           t_x_x_0_z, t_x_y_0_x, t_x_y_0_y, t_x_y_0_z, t_x_z_0_x, t_x_z_0_y,\
                           t_x_z_0_z, t_y_x_0_x, t_y_x_0_y, t_y_x_0_z, t_y_y_0_x, t_y_y_0_y,\
                           t_y_y_0_z, t_y_z_0_x, t_y_z_0_y, t_y_z_0_z, t_z_x_0_x, t_z_x_0_y,\
                           t_z_x_0_z, t_z_y_0_x, t_z_y_0_y, t_z_y_0_z, t_z_z_0_x, t_z_z_0_y,\
                           t_z_z_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_0_z[i] = t_0_zz_0_z[i] - rab_z[i] * t_0_z_0_z[i];

        t_z_z_0_y[i] = t_0_zz_0_y[i] - rab_z[i] * t_0_z_0_y[i];

        t_z_z_0_x[i] = t_0_zz_0_x[i] - rab_z[i] * t_0_z_0_x[i];

        t_z_y_0_z[i] = t_0_yz_0_z[i] - rab_z[i] * t_0_y_0_z[i];

        t_z_y_0_y[i] = t_0_yz_0_y[i] - rab_z[i] * t_0_y_0_y[i];

        t_z_y_0_x[i] = t_0_yz_0_x[i] - rab_z[i] * t_0_y_0_x[i];

        t_z_x_0_z[i] = t_0_xz_0_z[i] - rab_z[i] * t_0_x_0_z[i];

        t_z_x_0_y[i] = t_0_xz_0_y[i] - rab_z[i] * t_0_x_0_y[i];

        t_z_x_0_x[i] = t_0_xz_0_x[i] - rab_z[i] * t_0_x_0_x[i];

        t_y_z_0_z[i] = t_0_yz_0_z[i] - rab_y[i] * t_0_z_0_z[i];

        t_y_z_0_y[i] = t_0_yz_0_y[i] - rab_y[i] * t_0_z_0_y[i];

        t_y_z_0_x[i] = t_0_yz_0_x[i] - rab_y[i] * t_0_z_0_x[i];

        t_y_y_0_z[i] = t_0_yy_0_z[i] - rab_y[i] * t_0_y_0_z[i];

        t_y_y_0_y[i] = t_0_yy_0_y[i] - rab_y[i] * t_0_y_0_y[i];

        t_y_y_0_x[i] = t_0_yy_0_x[i] - rab_y[i] * t_0_y_0_x[i];

        t_y_x_0_z[i] = t_0_xy_0_z[i] - rab_y[i] * t_0_x_0_z[i];

        t_y_x_0_y[i] = t_0_xy_0_y[i] - rab_y[i] * t_0_x_0_y[i];

        t_y_x_0_x[i] = t_0_xy_0_x[i] - rab_y[i] * t_0_x_0_x[i];

        t_x_z_0_z[i] = t_0_xz_0_z[i] - rab_x[i] * t_0_z_0_z[i];

        t_x_z_0_y[i] = t_0_xz_0_y[i] - rab_x[i] * t_0_z_0_y[i];

        t_x_z_0_x[i] = t_0_xz_0_x[i] - rab_x[i] * t_0_z_0_x[i];

        t_x_y_0_z[i] = t_0_xy_0_z[i] - rab_x[i] * t_0_y_0_z[i];

        t_x_y_0_y[i] = t_0_xy_0_y[i] - rab_x[i] * t_0_y_0_y[i];

        t_x_y_0_x[i] = t_0_xy_0_x[i] - rab_x[i] * t_0_y_0_x[i];

        t_x_x_0_z[i] = t_0_xx_0_z[i] - rab_x[i] * t_0_x_0_z[i];

        t_x_x_0_y[i] = t_0_xx_0_y[i] - rab_x[i] * t_0_x_0_y[i];

        t_x_x_0_x[i] = t_0_xx_0_x[i] - rab_x[i] * t_0_x_0_x[i];
    }
}


} // derirec namespace
