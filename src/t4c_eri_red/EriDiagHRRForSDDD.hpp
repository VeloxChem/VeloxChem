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
compHostHRRForSDDD_V0(      BufferHostXY<T>&      intsBufferSDDD,
                      const BufferHostX<int32_t>& intsIndexesSDDD,
                      const BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostXY<T>&      intsBufferSDPF,
                      const BufferHostX<int32_t>& intsIndexesSDPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SDDD) integral components

    t_0_zz_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(0));

    t_0_zz_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(1));

    t_0_zz_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(2));

    t_0_zz_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(3));

    t_0_zz_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(4));

    t_0_zz_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(5));

    t_0_yz_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(6));

    t_0_yz_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(7));

    t_0_yz_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(8));

    t_0_yz_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(9));

    t_0_yz_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(10));

    t_0_yz_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(11));

    t_0_yy_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(12));

    t_0_yy_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(13));

    t_0_yy_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(14));

    t_0_yy_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(15));

    t_0_yy_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(16));

    t_0_yy_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(17));

    t_0_xz_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(18));

    t_0_xz_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(19));

    t_0_xz_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(20));

    t_0_xz_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(21));

    t_0_xz_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(22));

    t_0_xz_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(23));

    t_0_xy_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(24));

    t_0_xy_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(25));

    t_0_xy_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(26));

    t_0_xy_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(27));

    t_0_xy_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(28));

    t_0_xy_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(29));

    t_0_xx_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(30));

    t_0_xx_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(31));

    t_0_xx_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(32));

    t_0_xx_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(33));

    t_0_xx_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(34));

    t_0_xx_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(35));

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

    // set up (SDPF) integral components

    t_0_zz_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(0));

    t_0_zz_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(1));

    t_0_zz_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(2));

    t_0_zz_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(3));

    t_0_zz_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(4));

    t_0_zz_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(5));

    t_0_yz_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(6));

    t_0_yz_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(7));

    t_0_yz_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(8));

    t_0_yz_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(9));

    t_0_yz_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(10));

    t_0_yz_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(11));

    t_0_yy_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(12));

    t_0_yy_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(13));

    t_0_yy_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(14));

    t_0_yy_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(15));

    t_0_yy_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(16));

    t_0_yy_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(17));

    t_0_xz_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(18));

    t_0_xz_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(19));

    t_0_xz_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(20));

    t_0_xz_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(21));

    t_0_xz_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(22));

    t_0_xz_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(23));

    t_0_xy_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(24));

    t_0_xy_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(25));

    t_0_xy_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(26));

    t_0_xy_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(27));

    t_0_xy_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(28));

    t_0_xy_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(29));

    t_0_xx_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(30));

    t_0_xx_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(31));

    t_0_xx_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(32));

    t_0_xx_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(33));

    t_0_xx_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(34));

    t_0_xx_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(35));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_x_xx, t_0_xx_x_xxx, t_0_xx_x_xxy,\
                           t_0_xx_x_xxz, t_0_xx_xx_xx, t_0_xx_xy_xx, t_0_xx_xz_xx, t_0_xx_y_xx,\
                           t_0_xx_y_xxy, t_0_xx_y_xxz, t_0_xx_yy_xx, t_0_xx_yz_xx, t_0_xx_z_xx,\
                           t_0_xx_z_xxz, t_0_xx_zz_xx, t_0_xy_x_xxy, t_0_xy_x_xy, t_0_xy_x_xyy,\
                           t_0_xy_x_xyz, t_0_xy_xx_xy, t_0_xy_xy_xy, t_0_xy_xz_xy, t_0_xy_y_xy,\
                           t_0_xy_y_xyy, t_0_xy_y_xyz, t_0_xy_yy_xy, t_0_xy_yz_xy, t_0_xy_z_xy,\
                           t_0_xy_z_xyz, t_0_xy_zz_xy, t_0_xz_x_xxz, t_0_xz_x_xyz, t_0_xz_x_xz,\
                           t_0_xz_x_xzz, t_0_xz_xx_xz, t_0_xz_xy_xz, t_0_xz_xz_xz, t_0_xz_y_xyz,\
                           t_0_xz_y_xz, t_0_xz_y_xzz, t_0_xz_yy_xz, t_0_xz_yz_xz, t_0_xz_z_xz,\
                           t_0_xz_z_xzz, t_0_xz_zz_xz, t_0_yy_x_xyy, t_0_yy_x_yy, t_0_yy_x_yyy,\
                           t_0_yy_x_yyz, t_0_yy_xx_yy, t_0_yy_xy_yy, t_0_yy_xz_yy, t_0_yy_y_yy,\
                           t_0_yy_y_yyy, t_0_yy_y_yyz, t_0_yy_yy_yy, t_0_yy_yz_yy, t_0_yy_z_yy,\
                           t_0_yy_z_yyz, t_0_yy_zz_yy, t_0_yz_x_xyz, t_0_yz_x_yyz, t_0_yz_x_yz,\
                           t_0_yz_x_yzz, t_0_yz_xx_yz, t_0_yz_xy_yz, t_0_yz_xz_yz, t_0_yz_y_yyz,\
                           t_0_yz_y_yz, t_0_yz_y_yzz, t_0_yz_yy_yz, t_0_yz_yz_yz, t_0_yz_z_yz,\
                           t_0_yz_z_yzz, t_0_yz_zz_yz, t_0_zz_x_xzz, t_0_zz_x_yzz, t_0_zz_x_zz,\
                           t_0_zz_x_zzz, t_0_zz_xx_zz, t_0_zz_xy_zz, t_0_zz_xz_zz, t_0_zz_y_yzz,\
                           t_0_zz_y_zz, t_0_zz_y_zzz, t_0_zz_yy_zz, t_0_zz_yz_zz, t_0_zz_z_zz,\
                           t_0_zz_z_zzz, t_0_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_zz_zz[i] = t_0_zz_z_zzz[i] - rcd_z[i] * t_0_zz_z_zz[i];

        t_0_zz_yz_zz[i] = t_0_zz_y_zzz[i] - rcd_z[i] * t_0_zz_y_zz[i];

        t_0_zz_yy_zz[i] = t_0_zz_y_yzz[i] - rcd_y[i] * t_0_zz_y_zz[i];

        t_0_zz_xz_zz[i] = t_0_zz_x_zzz[i] - rcd_z[i] * t_0_zz_x_zz[i];

        t_0_zz_xy_zz[i] = t_0_zz_x_yzz[i] - rcd_y[i] * t_0_zz_x_zz[i];

        t_0_zz_xx_zz[i] = t_0_zz_x_xzz[i] - rcd_x[i] * t_0_zz_x_zz[i];

        t_0_yz_zz_yz[i] = t_0_yz_z_yzz[i] - rcd_z[i] * t_0_yz_z_yz[i];

        t_0_yz_yz_yz[i] = t_0_yz_y_yzz[i] - rcd_z[i] * t_0_yz_y_yz[i];

        t_0_yz_yy_yz[i] = t_0_yz_y_yyz[i] - rcd_y[i] * t_0_yz_y_yz[i];

        t_0_yz_xz_yz[i] = t_0_yz_x_yzz[i] - rcd_z[i] * t_0_yz_x_yz[i];

        t_0_yz_xy_yz[i] = t_0_yz_x_yyz[i] - rcd_y[i] * t_0_yz_x_yz[i];

        t_0_yz_xx_yz[i] = t_0_yz_x_xyz[i] - rcd_x[i] * t_0_yz_x_yz[i];

        t_0_yy_zz_yy[i] = t_0_yy_z_yyz[i] - rcd_z[i] * t_0_yy_z_yy[i];

        t_0_yy_yz_yy[i] = t_0_yy_y_yyz[i] - rcd_z[i] * t_0_yy_y_yy[i];

        t_0_yy_yy_yy[i] = t_0_yy_y_yyy[i] - rcd_y[i] * t_0_yy_y_yy[i];

        t_0_yy_xz_yy[i] = t_0_yy_x_yyz[i] - rcd_z[i] * t_0_yy_x_yy[i];

        t_0_yy_xy_yy[i] = t_0_yy_x_yyy[i] - rcd_y[i] * t_0_yy_x_yy[i];

        t_0_yy_xx_yy[i] = t_0_yy_x_xyy[i] - rcd_x[i] * t_0_yy_x_yy[i];

        t_0_xz_zz_xz[i] = t_0_xz_z_xzz[i] - rcd_z[i] * t_0_xz_z_xz[i];

        t_0_xz_yz_xz[i] = t_0_xz_y_xzz[i] - rcd_z[i] * t_0_xz_y_xz[i];

        t_0_xz_yy_xz[i] = t_0_xz_y_xyz[i] - rcd_y[i] * t_0_xz_y_xz[i];

        t_0_xz_xz_xz[i] = t_0_xz_x_xzz[i] - rcd_z[i] * t_0_xz_x_xz[i];

        t_0_xz_xy_xz[i] = t_0_xz_x_xyz[i] - rcd_y[i] * t_0_xz_x_xz[i];

        t_0_xz_xx_xz[i] = t_0_xz_x_xxz[i] - rcd_x[i] * t_0_xz_x_xz[i];

        t_0_xy_zz_xy[i] = t_0_xy_z_xyz[i] - rcd_z[i] * t_0_xy_z_xy[i];

        t_0_xy_yz_xy[i] = t_0_xy_y_xyz[i] - rcd_z[i] * t_0_xy_y_xy[i];

        t_0_xy_yy_xy[i] = t_0_xy_y_xyy[i] - rcd_y[i] * t_0_xy_y_xy[i];

        t_0_xy_xz_xy[i] = t_0_xy_x_xyz[i] - rcd_z[i] * t_0_xy_x_xy[i];

        t_0_xy_xy_xy[i] = t_0_xy_x_xyy[i] - rcd_y[i] * t_0_xy_x_xy[i];

        t_0_xy_xx_xy[i] = t_0_xy_x_xxy[i] - rcd_x[i] * t_0_xy_x_xy[i];

        t_0_xx_zz_xx[i] = t_0_xx_z_xxz[i] - rcd_z[i] * t_0_xx_z_xx[i];

        t_0_xx_yz_xx[i] = t_0_xx_y_xxz[i] - rcd_z[i] * t_0_xx_y_xx[i];

        t_0_xx_yy_xx[i] = t_0_xx_y_xxy[i] - rcd_y[i] * t_0_xx_y_xx[i];

        t_0_xx_xz_xx[i] = t_0_xx_x_xxz[i] - rcd_z[i] * t_0_xx_x_xx[i];

        t_0_xx_xy_xx[i] = t_0_xx_x_xxy[i] - rcd_y[i] * t_0_xx_x_xx[i];

        t_0_xx_xx_xx[i] = t_0_xx_x_xxx[i] - rcd_x[i] * t_0_xx_x_xx[i];
    }
}


} // derirec namespace
