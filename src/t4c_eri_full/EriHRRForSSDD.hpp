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
compHostHRRForSSDD_V0(      BufferHostXY<T>&      intsBufferSSDD,
                      const BufferHostX<int32_t>& intsIndexesSSDD,
                      const BufferHostXY<T>&      intsBufferSSPD,
                      const BufferHostX<int32_t>& intsIndexesSSPD,
                      const BufferHostXY<T>&      intsBufferSSPF,
                      const BufferHostX<int32_t>& intsIndexesSSPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SSDD) integral components

    t_0_0_zz_zz = intsBufferSSDD.data(intsIndexesSSDD(0));

    t_0_0_zz_yz = intsBufferSSDD.data(intsIndexesSSDD(1));

    t_0_0_zz_yy = intsBufferSSDD.data(intsIndexesSSDD(2));

    t_0_0_zz_xz = intsBufferSSDD.data(intsIndexesSSDD(3));

    t_0_0_zz_xy = intsBufferSSDD.data(intsIndexesSSDD(4));

    t_0_0_zz_xx = intsBufferSSDD.data(intsIndexesSSDD(5));

    t_0_0_yz_zz = intsBufferSSDD.data(intsIndexesSSDD(6));

    t_0_0_yz_yz = intsBufferSSDD.data(intsIndexesSSDD(7));

    t_0_0_yz_yy = intsBufferSSDD.data(intsIndexesSSDD(8));

    t_0_0_yz_xz = intsBufferSSDD.data(intsIndexesSSDD(9));

    t_0_0_yz_xy = intsBufferSSDD.data(intsIndexesSSDD(10));

    t_0_0_yz_xx = intsBufferSSDD.data(intsIndexesSSDD(11));

    t_0_0_yy_zz = intsBufferSSDD.data(intsIndexesSSDD(12));

    t_0_0_yy_yz = intsBufferSSDD.data(intsIndexesSSDD(13));

    t_0_0_yy_yy = intsBufferSSDD.data(intsIndexesSSDD(14));

    t_0_0_yy_xz = intsBufferSSDD.data(intsIndexesSSDD(15));

    t_0_0_yy_xy = intsBufferSSDD.data(intsIndexesSSDD(16));

    t_0_0_yy_xx = intsBufferSSDD.data(intsIndexesSSDD(17));

    t_0_0_xz_zz = intsBufferSSDD.data(intsIndexesSSDD(18));

    t_0_0_xz_yz = intsBufferSSDD.data(intsIndexesSSDD(19));

    t_0_0_xz_yy = intsBufferSSDD.data(intsIndexesSSDD(20));

    t_0_0_xz_xz = intsBufferSSDD.data(intsIndexesSSDD(21));

    t_0_0_xz_xy = intsBufferSSDD.data(intsIndexesSSDD(22));

    t_0_0_xz_xx = intsBufferSSDD.data(intsIndexesSSDD(23));

    t_0_0_xy_zz = intsBufferSSDD.data(intsIndexesSSDD(24));

    t_0_0_xy_yz = intsBufferSSDD.data(intsIndexesSSDD(25));

    t_0_0_xy_yy = intsBufferSSDD.data(intsIndexesSSDD(26));

    t_0_0_xy_xz = intsBufferSSDD.data(intsIndexesSSDD(27));

    t_0_0_xy_xy = intsBufferSSDD.data(intsIndexesSSDD(28));

    t_0_0_xy_xx = intsBufferSSDD.data(intsIndexesSSDD(29));

    t_0_0_xx_zz = intsBufferSSDD.data(intsIndexesSSDD(30));

    t_0_0_xx_yz = intsBufferSSDD.data(intsIndexesSSDD(31));

    t_0_0_xx_yy = intsBufferSSDD.data(intsIndexesSSDD(32));

    t_0_0_xx_xz = intsBufferSSDD.data(intsIndexesSSDD(33));

    t_0_0_xx_xy = intsBufferSSDD.data(intsIndexesSSDD(34));

    t_0_0_xx_xx = intsBufferSSDD.data(intsIndexesSSDD(35));

    // set up (SSPD) integral components

    t_0_0_z_zz = intsBufferSSPD.data(intsIndexesSSPD(0));

    t_0_0_z_yz = intsBufferSSPD.data(intsIndexesSSPD(1));

    t_0_0_z_yy = intsBufferSSPD.data(intsIndexesSSPD(2));

    t_0_0_z_xz = intsBufferSSPD.data(intsIndexesSSPD(3));

    t_0_0_z_xy = intsBufferSSPD.data(intsIndexesSSPD(4));

    t_0_0_z_xx = intsBufferSSPD.data(intsIndexesSSPD(5));

    t_0_0_y_zz = intsBufferSSPD.data(intsIndexesSSPD(6));

    t_0_0_y_yz = intsBufferSSPD.data(intsIndexesSSPD(7));

    t_0_0_y_yy = intsBufferSSPD.data(intsIndexesSSPD(8));

    t_0_0_y_xz = intsBufferSSPD.data(intsIndexesSSPD(9));

    t_0_0_y_xy = intsBufferSSPD.data(intsIndexesSSPD(10));

    t_0_0_y_xx = intsBufferSSPD.data(intsIndexesSSPD(11));

    t_0_0_x_zz = intsBufferSSPD.data(intsIndexesSSPD(12));

    t_0_0_x_yz = intsBufferSSPD.data(intsIndexesSSPD(13));

    t_0_0_x_yy = intsBufferSSPD.data(intsIndexesSSPD(14));

    t_0_0_x_xz = intsBufferSSPD.data(intsIndexesSSPD(15));

    t_0_0_x_xy = intsBufferSSPD.data(intsIndexesSSPD(16));

    t_0_0_x_xx = intsBufferSSPD.data(intsIndexesSSPD(17));

    // set up (SSPF) integral components

    t_0_0_z_zzz = intsBufferSSPF.data(intsIndexesSSPF(0));

    t_0_0_z_yzz = intsBufferSSPF.data(intsIndexesSSPF(1));

    t_0_0_z_yyz = intsBufferSSPF.data(intsIndexesSSPF(2));

    t_0_0_z_xzz = intsBufferSSPF.data(intsIndexesSSPF(3));

    t_0_0_z_xyz = intsBufferSSPF.data(intsIndexesSSPF(4));

    t_0_0_z_xxz = intsBufferSSPF.data(intsIndexesSSPF(5));

    t_0_0_y_zzz = intsBufferSSPF.data(intsIndexesSSPF(6));

    t_0_0_y_yzz = intsBufferSSPF.data(intsIndexesSSPF(7));

    t_0_0_y_yyz = intsBufferSSPF.data(intsIndexesSSPF(8));

    t_0_0_y_yyy = intsBufferSSPF.data(intsIndexesSSPF(9));

    t_0_0_y_xzz = intsBufferSSPF.data(intsIndexesSSPF(10));

    t_0_0_y_xyz = intsBufferSSPF.data(intsIndexesSSPF(11));

    t_0_0_y_xyy = intsBufferSSPF.data(intsIndexesSSPF(12));

    t_0_0_y_xxz = intsBufferSSPF.data(intsIndexesSSPF(13));

    t_0_0_y_xxy = intsBufferSSPF.data(intsIndexesSSPF(14));

    t_0_0_x_zzz = intsBufferSSPF.data(intsIndexesSSPF(15));

    t_0_0_x_yzz = intsBufferSSPF.data(intsIndexesSSPF(16));

    t_0_0_x_yyz = intsBufferSSPF.data(intsIndexesSSPF(17));

    t_0_0_x_yyy = intsBufferSSPF.data(intsIndexesSSPF(18));

    t_0_0_x_xzz = intsBufferSSPF.data(intsIndexesSSPF(19));

    t_0_0_x_xyz = intsBufferSSPF.data(intsIndexesSSPF(20));

    t_0_0_x_xyy = intsBufferSSPF.data(intsIndexesSSPF(21));

    t_0_0_x_xxz = intsBufferSSPF.data(intsIndexesSSPF(22));

    t_0_0_x_xxy = intsBufferSSPF.data(intsIndexesSSPF(23));

    t_0_0_x_xxx = intsBufferSSPF.data(intsIndexesSSPF(24));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_0_x_xx, t_0_0_x_xxx, t_0_0_x_xxy,\
                           t_0_0_x_xxz, t_0_0_x_xy, t_0_0_x_xyy, t_0_0_x_xyz, t_0_0_x_xz,\
                           t_0_0_x_xzz, t_0_0_x_yy, t_0_0_x_yyy, t_0_0_x_yyz, t_0_0_x_yz,\
                           t_0_0_x_yzz, t_0_0_x_zz, t_0_0_x_zzz, t_0_0_xx_xx, t_0_0_xx_xy,\
                           t_0_0_xx_xz, t_0_0_xx_yy, t_0_0_xx_yz, t_0_0_xx_zz, t_0_0_xy_xx,\
                           t_0_0_xy_xy, t_0_0_xy_xz, t_0_0_xy_yy, t_0_0_xy_yz, t_0_0_xy_zz,\
                           t_0_0_xz_xx, t_0_0_xz_xy, t_0_0_xz_xz, t_0_0_xz_yy, t_0_0_xz_yz,\
                           t_0_0_xz_zz, t_0_0_y_xx, t_0_0_y_xxy, t_0_0_y_xxz, t_0_0_y_xy,\
                           t_0_0_y_xyy, t_0_0_y_xyz, t_0_0_y_xz, t_0_0_y_xzz, t_0_0_y_yy,\
                           t_0_0_y_yyy, t_0_0_y_yyz, t_0_0_y_yz, t_0_0_y_yzz, t_0_0_y_zz,\
                           t_0_0_y_zzz, t_0_0_yy_xx, t_0_0_yy_xy, t_0_0_yy_xz, t_0_0_yy_yy,\
                           t_0_0_yy_yz, t_0_0_yy_zz, t_0_0_yz_xx, t_0_0_yz_xy, t_0_0_yz_xz,\
                           t_0_0_yz_yy, t_0_0_yz_yz, t_0_0_yz_zz, t_0_0_z_xx, t_0_0_z_xxz,\
                           t_0_0_z_xy, t_0_0_z_xyz, t_0_0_z_xz, t_0_0_z_xzz, t_0_0_z_yy,\
                           t_0_0_z_yyz, t_0_0_z_yz, t_0_0_z_yzz, t_0_0_z_zz, t_0_0_z_zzz,\
                           t_0_0_zz_xx, t_0_0_zz_xy, t_0_0_zz_xz, t_0_0_zz_yy, t_0_0_zz_yz,\
                           t_0_0_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_0_zz_zz[i] = t_0_0_z_zzz[i] - rcd_z[i] * t_0_0_z_zz[i];

        t_0_0_zz_yz[i] = t_0_0_z_yzz[i] - rcd_z[i] * t_0_0_z_yz[i];

        t_0_0_zz_yy[i] = t_0_0_z_yyz[i] - rcd_z[i] * t_0_0_z_yy[i];

        t_0_0_zz_xz[i] = t_0_0_z_xzz[i] - rcd_z[i] * t_0_0_z_xz[i];

        t_0_0_zz_xy[i] = t_0_0_z_xyz[i] - rcd_z[i] * t_0_0_z_xy[i];

        t_0_0_zz_xx[i] = t_0_0_z_xxz[i] - rcd_z[i] * t_0_0_z_xx[i];

        t_0_0_yz_zz[i] = t_0_0_y_zzz[i] - rcd_z[i] * t_0_0_y_zz[i];

        t_0_0_yz_yz[i] = t_0_0_y_yzz[i] - rcd_z[i] * t_0_0_y_yz[i];

        t_0_0_yz_yy[i] = t_0_0_y_yyz[i] - rcd_z[i] * t_0_0_y_yy[i];

        t_0_0_yz_xz[i] = t_0_0_y_xzz[i] - rcd_z[i] * t_0_0_y_xz[i];

        t_0_0_yz_xy[i] = t_0_0_y_xyz[i] - rcd_z[i] * t_0_0_y_xy[i];

        t_0_0_yz_xx[i] = t_0_0_y_xxz[i] - rcd_z[i] * t_0_0_y_xx[i];

        t_0_0_yy_zz[i] = t_0_0_y_yzz[i] - rcd_y[i] * t_0_0_y_zz[i];

        t_0_0_yy_yz[i] = t_0_0_y_yyz[i] - rcd_y[i] * t_0_0_y_yz[i];

        t_0_0_yy_yy[i] = t_0_0_y_yyy[i] - rcd_y[i] * t_0_0_y_yy[i];

        t_0_0_yy_xz[i] = t_0_0_y_xyz[i] - rcd_y[i] * t_0_0_y_xz[i];

        t_0_0_yy_xy[i] = t_0_0_y_xyy[i] - rcd_y[i] * t_0_0_y_xy[i];

        t_0_0_yy_xx[i] = t_0_0_y_xxy[i] - rcd_y[i] * t_0_0_y_xx[i];

        t_0_0_xz_zz[i] = t_0_0_x_zzz[i] - rcd_z[i] * t_0_0_x_zz[i];

        t_0_0_xz_yz[i] = t_0_0_x_yzz[i] - rcd_z[i] * t_0_0_x_yz[i];

        t_0_0_xz_yy[i] = t_0_0_x_yyz[i] - rcd_z[i] * t_0_0_x_yy[i];

        t_0_0_xz_xz[i] = t_0_0_x_xzz[i] - rcd_z[i] * t_0_0_x_xz[i];

        t_0_0_xz_xy[i] = t_0_0_x_xyz[i] - rcd_z[i] * t_0_0_x_xy[i];

        t_0_0_xz_xx[i] = t_0_0_x_xxz[i] - rcd_z[i] * t_0_0_x_xx[i];

        t_0_0_xy_zz[i] = t_0_0_x_yzz[i] - rcd_y[i] * t_0_0_x_zz[i];

        t_0_0_xy_yz[i] = t_0_0_x_yyz[i] - rcd_y[i] * t_0_0_x_yz[i];

        t_0_0_xy_yy[i] = t_0_0_x_yyy[i] - rcd_y[i] * t_0_0_x_yy[i];

        t_0_0_xy_xz[i] = t_0_0_x_xyz[i] - rcd_y[i] * t_0_0_x_xz[i];

        t_0_0_xy_xy[i] = t_0_0_x_xyy[i] - rcd_y[i] * t_0_0_x_xy[i];

        t_0_0_xy_xx[i] = t_0_0_x_xxy[i] - rcd_y[i] * t_0_0_x_xx[i];

        t_0_0_xx_zz[i] = t_0_0_x_xzz[i] - rcd_x[i] * t_0_0_x_zz[i];

        t_0_0_xx_yz[i] = t_0_0_x_xyz[i] - rcd_x[i] * t_0_0_x_yz[i];

        t_0_0_xx_yy[i] = t_0_0_x_xyy[i] - rcd_x[i] * t_0_0_x_yy[i];

        t_0_0_xx_xz[i] = t_0_0_x_xxz[i] - rcd_x[i] * t_0_0_x_xz[i];

        t_0_0_xx_xy[i] = t_0_0_x_xxy[i] - rcd_x[i] * t_0_0_x_xy[i];

        t_0_0_xx_xx[i] = t_0_0_x_xxx[i] - rcd_x[i] * t_0_0_x_xx[i];
    }
}


} // derirec namespace
