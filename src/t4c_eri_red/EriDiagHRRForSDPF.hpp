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
compHostHRRForSDPF_V0(      BufferHostXY<T>&      intsBufferSDPF,
                      const BufferHostX<int32_t>& intsIndexesSDPF,
                      const BufferHostXY<T>&      intsBufferSDSF,
                      const BufferHostX<int32_t>& intsIndexesSDSF,
                      const BufferHostXY<T>&      intsBufferSDSG,
                      const BufferHostX<int32_t>& intsIndexesSDSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SDSG) integral components

    t_0_zz_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(0));

    t_0_zz_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(1));

    t_0_zz_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(2));

    t_0_zz_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(3));

    t_0_zz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(4));

    t_0_zz_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(5));

    t_0_yz_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(6));

    t_0_yz_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(7));

    t_0_yz_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(8));

    t_0_yz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(9));

    t_0_yz_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(10));

    t_0_yz_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(11));

    t_0_yy_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(12));

    t_0_yy_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(13));

    t_0_yy_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(14));

    t_0_yy_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(15));

    t_0_yy_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(16));

    t_0_yy_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(17));

    t_0_xz_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(18));

    t_0_xz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(19));

    t_0_xz_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(20));

    t_0_xz_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(21));

    t_0_xz_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(22));

    t_0_xz_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(23));

    t_0_xy_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(24));

    t_0_xy_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(25));

    t_0_xy_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(26));

    t_0_xy_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(27));

    t_0_xy_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(28));

    t_0_xy_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(29));

    t_0_xx_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(30));

    t_0_xx_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(31));

    t_0_xx_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(32));

    t_0_xx_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(33));

    t_0_xx_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(34));

    t_0_xx_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(35));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_0_xxx, t_0_xx_0_xxxx, t_0_xx_0_xxxy,\
                           t_0_xx_0_xxxz, t_0_xx_0_xxy, t_0_xx_0_xxyy, t_0_xx_0_xxyz,\
                           t_0_xx_0_xxz, t_0_xx_0_xxzz, t_0_xx_x_xxx, t_0_xx_x_xxy,\
                           t_0_xx_x_xxz, t_0_xx_y_xxy, t_0_xx_y_xxz, t_0_xx_z_xxz, t_0_xy_0_xxxy,\
                           t_0_xy_0_xxy, t_0_xy_0_xxyy, t_0_xy_0_xxyz, t_0_xy_0_xyy,\
                           t_0_xy_0_xyyy, t_0_xy_0_xyyz, t_0_xy_0_xyz, t_0_xy_0_xyzz,\
                           t_0_xy_x_xxy, t_0_xy_x_xyy, t_0_xy_x_xyz, t_0_xy_y_xyy, t_0_xy_y_xyz,\
                           t_0_xy_z_xyz, t_0_xz_0_xxxz, t_0_xz_0_xxyz, t_0_xz_0_xxz,\
                           t_0_xz_0_xxzz, t_0_xz_0_xyyz, t_0_xz_0_xyz, t_0_xz_0_xyzz,\
                           t_0_xz_0_xzz, t_0_xz_0_xzzz, t_0_xz_x_xxz, t_0_xz_x_xyz,\
                           t_0_xz_x_xzz, t_0_xz_y_xyz, t_0_xz_y_xzz, t_0_xz_z_xzz, t_0_yy_0_xxyy,\
                           t_0_yy_0_xyy, t_0_yy_0_xyyy, t_0_yy_0_xyyz, t_0_yy_0_yyy,\
                           t_0_yy_0_yyyy, t_0_yy_0_yyyz, t_0_yy_0_yyz, t_0_yy_0_yyzz,\
                           t_0_yy_x_xyy, t_0_yy_x_yyy, t_0_yy_x_yyz, t_0_yy_y_yyy, t_0_yy_y_yyz,\
                           t_0_yy_z_yyz, t_0_yz_0_xxyz, t_0_yz_0_xyyz, t_0_yz_0_xyz,\
                           t_0_yz_0_xyzz, t_0_yz_0_yyyz, t_0_yz_0_yyz, t_0_yz_0_yyzz,\
                           t_0_yz_0_yzz, t_0_yz_0_yzzz, t_0_yz_x_xyz, t_0_yz_x_yyz,\
                           t_0_yz_x_yzz, t_0_yz_y_yyz, t_0_yz_y_yzz, t_0_yz_z_yzz, t_0_zz_0_xxzz,\
                           t_0_zz_0_xyzz, t_0_zz_0_xzz, t_0_zz_0_xzzz, t_0_zz_0_yyzz,\
                           t_0_zz_0_yzz, t_0_zz_0_yzzz, t_0_zz_0_zzz, t_0_zz_0_zzzz,\
                           t_0_zz_x_xzz, t_0_zz_x_yzz, t_0_zz_x_zzz, t_0_zz_y_yzz, t_0_zz_y_zzz,\
                           t_0_zz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_z_zzz[i] = t_0_zz_0_zzzz[i] - rcd_z[i] * t_0_zz_0_zzz[i];

        t_0_zz_y_zzz[i] = t_0_zz_0_yzzz[i] - rcd_y[i] * t_0_zz_0_zzz[i];

        t_0_zz_y_yzz[i] = t_0_zz_0_yyzz[i] - rcd_y[i] * t_0_zz_0_yzz[i];

        t_0_zz_x_zzz[i] = t_0_zz_0_xzzz[i] - rcd_x[i] * t_0_zz_0_zzz[i];

        t_0_zz_x_yzz[i] = t_0_zz_0_xyzz[i] - rcd_x[i] * t_0_zz_0_yzz[i];

        t_0_zz_x_xzz[i] = t_0_zz_0_xxzz[i] - rcd_x[i] * t_0_zz_0_xzz[i];

        t_0_yz_z_yzz[i] = t_0_yz_0_yzzz[i] - rcd_z[i] * t_0_yz_0_yzz[i];

        t_0_yz_y_yzz[i] = t_0_yz_0_yyzz[i] - rcd_y[i] * t_0_yz_0_yzz[i];

        t_0_yz_y_yyz[i] = t_0_yz_0_yyyz[i] - rcd_y[i] * t_0_yz_0_yyz[i];

        t_0_yz_x_yzz[i] = t_0_yz_0_xyzz[i] - rcd_x[i] * t_0_yz_0_yzz[i];

        t_0_yz_x_yyz[i] = t_0_yz_0_xyyz[i] - rcd_x[i] * t_0_yz_0_yyz[i];

        t_0_yz_x_xyz[i] = t_0_yz_0_xxyz[i] - rcd_x[i] * t_0_yz_0_xyz[i];

        t_0_yy_z_yyz[i] = t_0_yy_0_yyzz[i] - rcd_z[i] * t_0_yy_0_yyz[i];

        t_0_yy_y_yyz[i] = t_0_yy_0_yyyz[i] - rcd_y[i] * t_0_yy_0_yyz[i];

        t_0_yy_y_yyy[i] = t_0_yy_0_yyyy[i] - rcd_y[i] * t_0_yy_0_yyy[i];

        t_0_yy_x_yyz[i] = t_0_yy_0_xyyz[i] - rcd_x[i] * t_0_yy_0_yyz[i];

        t_0_yy_x_yyy[i] = t_0_yy_0_xyyy[i] - rcd_x[i] * t_0_yy_0_yyy[i];

        t_0_yy_x_xyy[i] = t_0_yy_0_xxyy[i] - rcd_x[i] * t_0_yy_0_xyy[i];

        t_0_xz_z_xzz[i] = t_0_xz_0_xzzz[i] - rcd_z[i] * t_0_xz_0_xzz[i];

        t_0_xz_y_xzz[i] = t_0_xz_0_xyzz[i] - rcd_y[i] * t_0_xz_0_xzz[i];

        t_0_xz_y_xyz[i] = t_0_xz_0_xyyz[i] - rcd_y[i] * t_0_xz_0_xyz[i];

        t_0_xz_x_xzz[i] = t_0_xz_0_xxzz[i] - rcd_x[i] * t_0_xz_0_xzz[i];

        t_0_xz_x_xyz[i] = t_0_xz_0_xxyz[i] - rcd_x[i] * t_0_xz_0_xyz[i];

        t_0_xz_x_xxz[i] = t_0_xz_0_xxxz[i] - rcd_x[i] * t_0_xz_0_xxz[i];

        t_0_xy_z_xyz[i] = t_0_xy_0_xyzz[i] - rcd_z[i] * t_0_xy_0_xyz[i];

        t_0_xy_y_xyz[i] = t_0_xy_0_xyyz[i] - rcd_y[i] * t_0_xy_0_xyz[i];

        t_0_xy_y_xyy[i] = t_0_xy_0_xyyy[i] - rcd_y[i] * t_0_xy_0_xyy[i];

        t_0_xy_x_xyz[i] = t_0_xy_0_xxyz[i] - rcd_x[i] * t_0_xy_0_xyz[i];

        t_0_xy_x_xyy[i] = t_0_xy_0_xxyy[i] - rcd_x[i] * t_0_xy_0_xyy[i];

        t_0_xy_x_xxy[i] = t_0_xy_0_xxxy[i] - rcd_x[i] * t_0_xy_0_xxy[i];

        t_0_xx_z_xxz[i] = t_0_xx_0_xxzz[i] - rcd_z[i] * t_0_xx_0_xxz[i];

        t_0_xx_y_xxz[i] = t_0_xx_0_xxyz[i] - rcd_y[i] * t_0_xx_0_xxz[i];

        t_0_xx_y_xxy[i] = t_0_xx_0_xxyy[i] - rcd_y[i] * t_0_xx_0_xxy[i];

        t_0_xx_x_xxz[i] = t_0_xx_0_xxxz[i] - rcd_x[i] * t_0_xx_0_xxz[i];

        t_0_xx_x_xxy[i] = t_0_xx_0_xxxy[i] - rcd_x[i] * t_0_xx_0_xxy[i];

        t_0_xx_x_xxx[i] = t_0_xx_0_xxxx[i] - rcd_x[i] * t_0_xx_0_xxx[i];
    }
}


} // derirec namespace
