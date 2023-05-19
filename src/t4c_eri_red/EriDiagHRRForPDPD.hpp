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
compHostHRRForPDPD_V0(      BufferHostXY<T>&      intsBufferPDPD,
                      const BufferHostX<int32_t>& intsIndexesPDPD,
                      const BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PDPD) integral components

    t_z_zz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(0));

    t_z_yz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(1));

    t_z_yy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(2));

    t_z_xz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(3));

    t_z_xy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(4));

    t_z_xx_z_xx = intsBufferPDPD.data(intsIndexesPDPD(5));

    t_y_zz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(6));

    t_y_yz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(7));

    t_y_yy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(8));

    t_y_xz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(9));

    t_y_xy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(10));

    t_y_xx_y_xx = intsBufferPDPD.data(intsIndexesPDPD(11));

    t_x_zz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(12));

    t_x_yz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(13));

    t_x_yy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(14));

    t_x_xz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(15));

    t_x_xy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(16));

    t_x_xx_x_xx = intsBufferPDPD.data(intsIndexesPDPD(17));

    // set up (SDPD) integral components

    t_0_zz_z_zz = intsBufferSDPD.data(intsIndexesSDPD(0));

    t_0_zz_y_zz = intsBufferSDPD.data(intsIndexesSDPD(1));

    t_0_zz_x_zz = intsBufferSDPD.data(intsIndexesSDPD(2));

    t_0_yz_z_yz = intsBufferSDPD.data(intsIndexesSDPD(3));

    t_0_yz_y_yz = intsBufferSDPD.data(intsIndexesSDPD(4));

    t_0_yz_x_yz = intsBufferSDPD.data(intsIndexesSDPD(5));

    t_0_yy_z_yy = intsBufferSDPD.data(intsIndexesSDPD(6));

    t_0_yy_y_yy = intsBufferSDPD.data(intsIndexesSDPD(7));

    t_0_yy_x_yy = intsBufferSDPD.data(intsIndexesSDPD(8));

    t_0_xz_z_xz = intsBufferSDPD.data(intsIndexesSDPD(9));

    t_0_xz_y_xz = intsBufferSDPD.data(intsIndexesSDPD(10));

    t_0_xz_x_xz = intsBufferSDPD.data(intsIndexesSDPD(11));

    t_0_xy_z_xy = intsBufferSDPD.data(intsIndexesSDPD(12));

    t_0_xy_y_xy = intsBufferSDPD.data(intsIndexesSDPD(13));

    t_0_xy_x_xy = intsBufferSDPD.data(intsIndexesSDPD(14));

    t_0_xx_z_xx = intsBufferSDPD.data(intsIndexesSDPD(15));

    t_0_xx_y_xx = intsBufferSDPD.data(intsIndexesSDPD(16));

    t_0_xx_x_xx = intsBufferSDPD.data(intsIndexesSDPD(17));

    // set up (SFPD) integral components

    t_0_zzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(0));

    t_0_yzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(1));

    t_0_yzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(2));

    t_0_yyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(3));

    t_0_yyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(4));

    t_0_yyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(5));

    t_0_xzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(6));

    t_0_xzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(7));

    t_0_xyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(8));

    t_0_xyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(9));

    t_0_xyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(10));

    t_0_xyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(11));

    t_0_xyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(12));

    t_0_xxz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(13));

    t_0_xxz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(14));

    t_0_xxy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(15));

    t_0_xxy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(16));

    t_0_xxx_x_xx = intsBufferSFPD.data(intsIndexesSFPD(17));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_xx_x_xx, t_0_xx_y_xx, t_0_xx_z_xx,\
                           t_0_xxx_x_xx, t_0_xxy_x_xy, t_0_xxy_y_xx, t_0_xxz_x_xz, t_0_xxz_z_xx,\
                           t_0_xy_x_xy, t_0_xy_y_xy, t_0_xy_z_xy, t_0_xyy_x_yy, t_0_xyy_y_xy,\
                           t_0_xyz_x_yz, t_0_xyz_y_xz, t_0_xyz_z_xy, t_0_xz_x_xz, t_0_xz_y_xz,\
                           t_0_xz_z_xz, t_0_xzz_x_zz, t_0_xzz_z_xz, t_0_yy_x_yy, t_0_yy_y_yy,\
                           t_0_yy_z_yy, t_0_yyy_y_yy, t_0_yyz_y_yz, t_0_yyz_z_yy, t_0_yz_x_yz,\
                           t_0_yz_y_yz, t_0_yz_z_yz, t_0_yzz_y_zz, t_0_yzz_z_yz, t_0_zz_x_zz,\
                           t_0_zz_y_zz, t_0_zz_z_zz, t_0_zzz_z_zz, t_x_xx_x_xx, t_x_xy_x_xy,\
                           t_x_xz_x_xz, t_x_yy_x_yy, t_x_yz_x_yz, t_x_zz_x_zz, t_y_xx_y_xx,\
                           t_y_xy_y_xy, t_y_xz_y_xz, t_y_yy_y_yy, t_y_yz_y_yz, t_y_zz_y_zz,\
                           t_z_xx_z_xx, t_z_xy_z_xy, t_z_xz_z_xz, t_z_yy_z_yy, t_z_yz_z_yz,\
                           t_z_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_z_zz[i] = t_0_zzz_z_zz[i] - rab_z[i] * t_0_zz_z_zz[i];

        t_z_yz_z_yz[i] = t_0_yzz_z_yz[i] - rab_z[i] * t_0_yz_z_yz[i];

        t_z_yy_z_yy[i] = t_0_yyz_z_yy[i] - rab_z[i] * t_0_yy_z_yy[i];

        t_z_xz_z_xz[i] = t_0_xzz_z_xz[i] - rab_z[i] * t_0_xz_z_xz[i];

        t_z_xy_z_xy[i] = t_0_xyz_z_xy[i] - rab_z[i] * t_0_xy_z_xy[i];

        t_z_xx_z_xx[i] = t_0_xxz_z_xx[i] - rab_z[i] * t_0_xx_z_xx[i];

        t_y_zz_y_zz[i] = t_0_yzz_y_zz[i] - rab_y[i] * t_0_zz_y_zz[i];

        t_y_yz_y_yz[i] = t_0_yyz_y_yz[i] - rab_y[i] * t_0_yz_y_yz[i];

        t_y_yy_y_yy[i] = t_0_yyy_y_yy[i] - rab_y[i] * t_0_yy_y_yy[i];

        t_y_xz_y_xz[i] = t_0_xyz_y_xz[i] - rab_y[i] * t_0_xz_y_xz[i];

        t_y_xy_y_xy[i] = t_0_xyy_y_xy[i] - rab_y[i] * t_0_xy_y_xy[i];

        t_y_xx_y_xx[i] = t_0_xxy_y_xx[i] - rab_y[i] * t_0_xx_y_xx[i];

        t_x_zz_x_zz[i] = t_0_xzz_x_zz[i] - rab_x[i] * t_0_zz_x_zz[i];

        t_x_yz_x_yz[i] = t_0_xyz_x_yz[i] - rab_x[i] * t_0_yz_x_yz[i];

        t_x_yy_x_yy[i] = t_0_xyy_x_yy[i] - rab_x[i] * t_0_yy_x_yy[i];

        t_x_xz_x_xz[i] = t_0_xxz_x_xz[i] - rab_x[i] * t_0_xz_x_xz[i];

        t_x_xy_x_xy[i] = t_0_xxy_x_xy[i] - rab_x[i] * t_0_xy_x_xy[i];

        t_x_xx_x_xx[i] = t_0_xxx_x_xx[i] - rab_x[i] * t_0_xx_x_xx[i];
    }
}


} // derirec namespace
