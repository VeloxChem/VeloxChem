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
compHostHRRForSDPD_V0(      BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostXY<T>&      intsBufferSDSD,
                      const BufferHostX<int32_t>& intsIndexesSDSD,
                      const BufferHostXY<T>&      intsBufferSDSF,
                      const BufferHostX<int32_t>& intsIndexesSDSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SDSD) integral components

    t_0_zz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(0));

    t_0_yz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(1));

    t_0_yy_0_yy = intsBufferSDSD.data(intsIndexesSDSD(2));

    t_0_xz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(3));

    t_0_xy_0_xy = intsBufferSDSD.data(intsIndexesSDSD(4));

    t_0_xx_0_xx = intsBufferSDSD.data(intsIndexesSDSD(5));

    // set up (SDSF) integral components

    t_0_zz_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(0));

    t_0_zz_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(1));

    t_0_zz_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(2));

    t_0_yz_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(3));

    t_0_yz_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(4));

    t_0_yz_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(5));

    t_0_yy_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(6));

    t_0_yy_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(7));

    t_0_yy_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(8));

    t_0_xz_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(9));

    t_0_xz_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(10));

    t_0_xz_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(11));

    t_0_xy_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(12));

    t_0_xy_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(13));

    t_0_xy_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(14));

    t_0_xx_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(15));

    t_0_xx_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(16));

    t_0_xx_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(17));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_0_xx, t_0_xx_0_xxx, t_0_xx_0_xxy,\
                           t_0_xx_0_xxz, t_0_xx_x_xx, t_0_xx_y_xx, t_0_xx_z_xx, t_0_xy_0_xxy,\
                           t_0_xy_0_xy, t_0_xy_0_xyy, t_0_xy_0_xyz, t_0_xy_x_xy, t_0_xy_y_xy,\
                           t_0_xy_z_xy, t_0_xz_0_xxz, t_0_xz_0_xyz, t_0_xz_0_xz, t_0_xz_0_xzz,\
                           t_0_xz_x_xz, t_0_xz_y_xz, t_0_xz_z_xz, t_0_yy_0_xyy, t_0_yy_0_yy,\
                           t_0_yy_0_yyy, t_0_yy_0_yyz, t_0_yy_x_yy, t_0_yy_y_yy, t_0_yy_z_yy,\
                           t_0_yz_0_xyz, t_0_yz_0_yyz, t_0_yz_0_yz, t_0_yz_0_yzz, t_0_yz_x_yz,\
                           t_0_yz_y_yz, t_0_yz_z_yz, t_0_zz_0_xzz, t_0_zz_0_yzz, t_0_zz_0_zz,\
                           t_0_zz_0_zzz, t_0_zz_x_zz, t_0_zz_y_zz, t_0_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_z_zz[i] = t_0_zz_0_zzz[i] - rcd_z[i] * t_0_zz_0_zz[i];

        t_0_zz_y_zz[i] = t_0_zz_0_yzz[i] - rcd_y[i] * t_0_zz_0_zz[i];

        t_0_zz_x_zz[i] = t_0_zz_0_xzz[i] - rcd_x[i] * t_0_zz_0_zz[i];

        t_0_yz_z_yz[i] = t_0_yz_0_yzz[i] - rcd_z[i] * t_0_yz_0_yz[i];

        t_0_yz_y_yz[i] = t_0_yz_0_yyz[i] - rcd_y[i] * t_0_yz_0_yz[i];

        t_0_yz_x_yz[i] = t_0_yz_0_xyz[i] - rcd_x[i] * t_0_yz_0_yz[i];

        t_0_yy_z_yy[i] = t_0_yy_0_yyz[i] - rcd_z[i] * t_0_yy_0_yy[i];

        t_0_yy_y_yy[i] = t_0_yy_0_yyy[i] - rcd_y[i] * t_0_yy_0_yy[i];

        t_0_yy_x_yy[i] = t_0_yy_0_xyy[i] - rcd_x[i] * t_0_yy_0_yy[i];

        t_0_xz_z_xz[i] = t_0_xz_0_xzz[i] - rcd_z[i] * t_0_xz_0_xz[i];

        t_0_xz_y_xz[i] = t_0_xz_0_xyz[i] - rcd_y[i] * t_0_xz_0_xz[i];

        t_0_xz_x_xz[i] = t_0_xz_0_xxz[i] - rcd_x[i] * t_0_xz_0_xz[i];

        t_0_xy_z_xy[i] = t_0_xy_0_xyz[i] - rcd_z[i] * t_0_xy_0_xy[i];

        t_0_xy_y_xy[i] = t_0_xy_0_xyy[i] - rcd_y[i] * t_0_xy_0_xy[i];

        t_0_xy_x_xy[i] = t_0_xy_0_xxy[i] - rcd_x[i] * t_0_xy_0_xy[i];

        t_0_xx_z_xx[i] = t_0_xx_0_xxz[i] - rcd_z[i] * t_0_xx_0_xx[i];

        t_0_xx_y_xx[i] = t_0_xx_0_xxy[i] - rcd_y[i] * t_0_xx_0_xx[i];

        t_0_xx_x_xx[i] = t_0_xx_0_xxx[i] - rcd_x[i] * t_0_xx_0_xx[i];
    }
}


} // derirec namespace
