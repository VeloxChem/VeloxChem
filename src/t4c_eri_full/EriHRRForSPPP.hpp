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
compHostHRRForSPPP_V0(      BufferHostXY<T>&      intsBufferSPPP,
                      const BufferHostX<int32_t>& intsIndexesSPPP,
                      const BufferHostXY<T>&      intsBufferSPSP,
                      const BufferHostX<int32_t>& intsIndexesSPSP,
                      const BufferHostXY<T>&      intsBufferSPSD,
                      const BufferHostX<int32_t>& intsIndexesSPSD,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SPPP) integral components

    t_0_z_z_z = intsBufferSPPP.data(intsIndexesSPPP(0));

    t_0_z_z_y = intsBufferSPPP.data(intsIndexesSPPP(1));

    t_0_z_z_x = intsBufferSPPP.data(intsIndexesSPPP(2));

    t_0_z_y_z = intsBufferSPPP.data(intsIndexesSPPP(3));

    t_0_z_y_y = intsBufferSPPP.data(intsIndexesSPPP(4));

    t_0_z_y_x = intsBufferSPPP.data(intsIndexesSPPP(5));

    t_0_z_x_z = intsBufferSPPP.data(intsIndexesSPPP(6));

    t_0_z_x_y = intsBufferSPPP.data(intsIndexesSPPP(7));

    t_0_z_x_x = intsBufferSPPP.data(intsIndexesSPPP(8));

    t_0_y_z_z = intsBufferSPPP.data(intsIndexesSPPP(9));

    t_0_y_z_y = intsBufferSPPP.data(intsIndexesSPPP(10));

    t_0_y_z_x = intsBufferSPPP.data(intsIndexesSPPP(11));

    t_0_y_y_z = intsBufferSPPP.data(intsIndexesSPPP(12));

    t_0_y_y_y = intsBufferSPPP.data(intsIndexesSPPP(13));

    t_0_y_y_x = intsBufferSPPP.data(intsIndexesSPPP(14));

    t_0_y_x_z = intsBufferSPPP.data(intsIndexesSPPP(15));

    t_0_y_x_y = intsBufferSPPP.data(intsIndexesSPPP(16));

    t_0_y_x_x = intsBufferSPPP.data(intsIndexesSPPP(17));

    t_0_x_z_z = intsBufferSPPP.data(intsIndexesSPPP(18));

    t_0_x_z_y = intsBufferSPPP.data(intsIndexesSPPP(19));

    t_0_x_z_x = intsBufferSPPP.data(intsIndexesSPPP(20));

    t_0_x_y_z = intsBufferSPPP.data(intsIndexesSPPP(21));

    t_0_x_y_y = intsBufferSPPP.data(intsIndexesSPPP(22));

    t_0_x_y_x = intsBufferSPPP.data(intsIndexesSPPP(23));

    t_0_x_x_z = intsBufferSPPP.data(intsIndexesSPPP(24));

    t_0_x_x_y = intsBufferSPPP.data(intsIndexesSPPP(25));

    t_0_x_x_x = intsBufferSPPP.data(intsIndexesSPPP(26));

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

    // set up (SPSD) integral components

    t_0_z_0_zz = intsBufferSPSD.data(intsIndexesSPSD(0));

    t_0_z_0_yz = intsBufferSPSD.data(intsIndexesSPSD(1));

    t_0_z_0_yy = intsBufferSPSD.data(intsIndexesSPSD(2));

    t_0_z_0_xz = intsBufferSPSD.data(intsIndexesSPSD(3));

    t_0_z_0_xy = intsBufferSPSD.data(intsIndexesSPSD(4));

    t_0_z_0_xx = intsBufferSPSD.data(intsIndexesSPSD(5));

    t_0_y_0_zz = intsBufferSPSD.data(intsIndexesSPSD(6));

    t_0_y_0_yz = intsBufferSPSD.data(intsIndexesSPSD(7));

    t_0_y_0_yy = intsBufferSPSD.data(intsIndexesSPSD(8));

    t_0_y_0_xz = intsBufferSPSD.data(intsIndexesSPSD(9));

    t_0_y_0_xy = intsBufferSPSD.data(intsIndexesSPSD(10));

    t_0_y_0_xx = intsBufferSPSD.data(intsIndexesSPSD(11));

    t_0_x_0_zz = intsBufferSPSD.data(intsIndexesSPSD(12));

    t_0_x_0_yz = intsBufferSPSD.data(intsIndexesSPSD(13));

    t_0_x_0_yy = intsBufferSPSD.data(intsIndexesSPSD(14));

    t_0_x_0_xz = intsBufferSPSD.data(intsIndexesSPSD(15));

    t_0_x_0_xy = intsBufferSPSD.data(intsIndexesSPSD(16));

    t_0_x_0_xx = intsBufferSPSD.data(intsIndexesSPSD(17));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_x_0_x, t_0_x_0_xx, t_0_x_0_xy, t_0_x_0_xz,\
                           t_0_x_0_y, t_0_x_0_yy, t_0_x_0_yz, t_0_x_0_z, t_0_x_0_zz,\
                           t_0_x_x_x, t_0_x_x_y, t_0_x_x_z, t_0_x_y_x, t_0_x_y_y, t_0_x_y_z,\
                           t_0_x_z_x, t_0_x_z_y, t_0_x_z_z, t_0_y_0_x, t_0_y_0_xx, t_0_y_0_xy,\
                           t_0_y_0_xz, t_0_y_0_y, t_0_y_0_yy, t_0_y_0_yz, t_0_y_0_z,\
                           t_0_y_0_zz, t_0_y_x_x, t_0_y_x_y, t_0_y_x_z, t_0_y_y_x, t_0_y_y_y,\
                           t_0_y_y_z, t_0_y_z_x, t_0_y_z_y, t_0_y_z_z, t_0_z_0_x, t_0_z_0_xx,\
                           t_0_z_0_xy, t_0_z_0_xz, t_0_z_0_y, t_0_z_0_yy, t_0_z_0_yz,\
                           t_0_z_0_z, t_0_z_0_zz, t_0_z_x_x, t_0_z_x_y, t_0_z_x_z, t_0_z_y_x,\
                           t_0_z_y_y, t_0_z_y_z, t_0_z_z_x, t_0_z_z_y, t_0_z_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_z_z_z[i] = t_0_z_0_zz[i] - rcd_z[i] * t_0_z_0_z[i];

        t_0_z_z_y[i] = t_0_z_0_yz[i] - rcd_z[i] * t_0_z_0_y[i];

        t_0_z_z_x[i] = t_0_z_0_xz[i] - rcd_z[i] * t_0_z_0_x[i];

        t_0_z_y_z[i] = t_0_z_0_yz[i] - rcd_y[i] * t_0_z_0_z[i];

        t_0_z_y_y[i] = t_0_z_0_yy[i] - rcd_y[i] * t_0_z_0_y[i];

        t_0_z_y_x[i] = t_0_z_0_xy[i] - rcd_y[i] * t_0_z_0_x[i];

        t_0_z_x_z[i] = t_0_z_0_xz[i] - rcd_x[i] * t_0_z_0_z[i];

        t_0_z_x_y[i] = t_0_z_0_xy[i] - rcd_x[i] * t_0_z_0_y[i];

        t_0_z_x_x[i] = t_0_z_0_xx[i] - rcd_x[i] * t_0_z_0_x[i];

        t_0_y_z_z[i] = t_0_y_0_zz[i] - rcd_z[i] * t_0_y_0_z[i];

        t_0_y_z_y[i] = t_0_y_0_yz[i] - rcd_z[i] * t_0_y_0_y[i];

        t_0_y_z_x[i] = t_0_y_0_xz[i] - rcd_z[i] * t_0_y_0_x[i];

        t_0_y_y_z[i] = t_0_y_0_yz[i] - rcd_y[i] * t_0_y_0_z[i];

        t_0_y_y_y[i] = t_0_y_0_yy[i] - rcd_y[i] * t_0_y_0_y[i];

        t_0_y_y_x[i] = t_0_y_0_xy[i] - rcd_y[i] * t_0_y_0_x[i];

        t_0_y_x_z[i] = t_0_y_0_xz[i] - rcd_x[i] * t_0_y_0_z[i];

        t_0_y_x_y[i] = t_0_y_0_xy[i] - rcd_x[i] * t_0_y_0_y[i];

        t_0_y_x_x[i] = t_0_y_0_xx[i] - rcd_x[i] * t_0_y_0_x[i];

        t_0_x_z_z[i] = t_0_x_0_zz[i] - rcd_z[i] * t_0_x_0_z[i];

        t_0_x_z_y[i] = t_0_x_0_yz[i] - rcd_z[i] * t_0_x_0_y[i];

        t_0_x_z_x[i] = t_0_x_0_xz[i] - rcd_z[i] * t_0_x_0_x[i];

        t_0_x_y_z[i] = t_0_x_0_yz[i] - rcd_y[i] * t_0_x_0_z[i];

        t_0_x_y_y[i] = t_0_x_0_yy[i] - rcd_y[i] * t_0_x_0_y[i];

        t_0_x_y_x[i] = t_0_x_0_xy[i] - rcd_y[i] * t_0_x_0_x[i];

        t_0_x_x_z[i] = t_0_x_0_xz[i] - rcd_x[i] * t_0_x_0_z[i];

        t_0_x_x_y[i] = t_0_x_0_xy[i] - rcd_x[i] * t_0_x_0_y[i];

        t_0_x_x_x[i] = t_0_x_0_xx[i] - rcd_x[i] * t_0_x_0_x[i];
    }
}


} // derirec namespace
