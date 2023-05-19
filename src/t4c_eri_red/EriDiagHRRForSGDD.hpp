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
compHostHRRForSGDD_V0(      BufferHostXY<T>&      intsBufferSGDD,
                      const BufferHostX<int32_t>& intsIndexesSGDD,
                      const BufferHostXY<T>&      intsBufferSGPD,
                      const BufferHostX<int32_t>& intsIndexesSGPD,
                      const BufferHostXY<T>&      intsBufferSGPF,
                      const BufferHostX<int32_t>& intsIndexesSGPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SGDD) integral components

    t_0_zzzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(0));

    t_0_yzzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(1));

    t_0_yzzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(2));

    t_0_yyzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(3));

    t_0_yyzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(4));

    t_0_yyzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(5));

    t_0_yyyz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(6));

    t_0_yyyz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(7));

    t_0_yyyy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(8));

    t_0_xzzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(9));

    t_0_xzzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(10));

    t_0_xyzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(11));

    t_0_xyzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(12));

    t_0_xyzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(13));

    t_0_xyzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(14));

    t_0_xyyz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(15));

    t_0_xyyz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(16));

    t_0_xyyz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(17));

    t_0_xyyz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(18));

    t_0_xyyy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(19));

    t_0_xyyy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(20));

    t_0_xxzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(21));

    t_0_xxzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(22));

    t_0_xxzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(23));

    t_0_xxyz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(24));

    t_0_xxyz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(25));

    t_0_xxyz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(26));

    t_0_xxyz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(27));

    t_0_xxyy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(28));

    t_0_xxyy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(29));

    t_0_xxyy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(30));

    t_0_xxxz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(31));

    t_0_xxxz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(32));

    t_0_xxxy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(33));

    t_0_xxxy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(34));

    t_0_xxxx_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(35));

    // set up (SGPD) integral components

    t_0_zzzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(0));

    t_0_yzzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(1));

    t_0_yzzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(2));

    t_0_yyzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(3));

    t_0_yyzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(4));

    t_0_yyzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(5));

    t_0_yyyz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(6));

    t_0_yyyz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(7));

    t_0_yyyy_y_yy = intsBufferSGPD.data(intsIndexesSGPD(8));

    t_0_xzzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(9));

    t_0_xzzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(10));

    t_0_xyzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(11));

    t_0_xyzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(12));

    t_0_xyzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(13));

    t_0_xyzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(14));

    t_0_xyyz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(15));

    t_0_xyyz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(16));

    t_0_xyyz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(17));

    t_0_xyyz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(18));

    t_0_xyyy_y_xy = intsBufferSGPD.data(intsIndexesSGPD(19));

    t_0_xyyy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(20));

    t_0_xxzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(21));

    t_0_xxzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(22));

    t_0_xxzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(23));

    t_0_xxyz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(24));

    t_0_xxyz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(25));

    t_0_xxyz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(26));

    t_0_xxyz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(27));

    t_0_xxyy_y_xx = intsBufferSGPD.data(intsIndexesSGPD(28));

    t_0_xxyy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(29));

    t_0_xxyy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(30));

    t_0_xxxz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(31));

    t_0_xxxz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(32));

    t_0_xxxy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(33));

    t_0_xxxy_x_xx = intsBufferSGPD.data(intsIndexesSGPD(34));

    t_0_xxxx_x_xx = intsBufferSGPD.data(intsIndexesSGPD(35));

    // set up (SGPF) integral components

    t_0_zzzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(0));

    t_0_yzzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(1));

    t_0_yzzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(2));

    t_0_yyzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(3));

    t_0_yyzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(4));

    t_0_yyyz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(5));

    t_0_yyyy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(6));

    t_0_xzzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(7));

    t_0_xzzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(8));

    t_0_xyzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(9));

    t_0_xyzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(10));

    t_0_xyzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(11));

    t_0_xyyz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(12));

    t_0_xyyz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(13));

    t_0_xyyy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(14));

    t_0_xyyy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(15));

    t_0_xxzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(16));

    t_0_xxzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(17));

    t_0_xxyz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(18));

    t_0_xxyz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(19));

    t_0_xxyy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(20));

    t_0_xxyy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(21));

    t_0_xxxz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(22));

    t_0_xxxy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(23));

    t_0_xxxx_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(24));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_x_xx, t_0_xxxx_x_xxx, t_0_xxxx_xx_xx,\
                           t_0_xxxy_x_xx, t_0_xxxy_x_xxy, t_0_xxxy_x_xy, t_0_xxxy_xx_xy,\
                           t_0_xxxy_xy_xx, t_0_xxxz_x_xx, t_0_xxxz_x_xxz, t_0_xxxz_x_xz,\
                           t_0_xxxz_xx_xz, t_0_xxxz_xz_xx, t_0_xxyy_x_xy, t_0_xxyy_x_xyy,\
                           t_0_xxyy_x_yy, t_0_xxyy_xx_yy, t_0_xxyy_xy_xy, t_0_xxyy_y_xx,\
                           t_0_xxyy_y_xxy, t_0_xxyy_yy_xx, t_0_xxyz_x_xy, t_0_xxyz_x_xyz,\
                           t_0_xxyz_x_xz, t_0_xxyz_x_yz, t_0_xxyz_xx_yz, t_0_xxyz_xy_xz,\
                           t_0_xxyz_xz_xy, t_0_xxyz_y_xx, t_0_xxyz_y_xxz, t_0_xxyz_yz_xx,\
                           t_0_xxzz_x_xz, t_0_xxzz_x_xzz, t_0_xxzz_x_zz, t_0_xxzz_xx_zz,\
                           t_0_xxzz_xz_xz, t_0_xxzz_z_xx, t_0_xxzz_z_xxz, t_0_xxzz_zz_xx,\
                           t_0_xyyy_x_yy, t_0_xyyy_x_yyy, t_0_xyyy_xy_yy, t_0_xyyy_y_xy,\
                           t_0_xyyy_y_xyy, t_0_xyyy_yy_xy, t_0_xyyz_x_yy, t_0_xyyz_x_yyz,\
                           t_0_xyyz_x_yz, t_0_xyyz_xy_yz, t_0_xyyz_xz_yy, t_0_xyyz_y_xy,\
                           t_0_xyyz_y_xyz, t_0_xyyz_y_xz, t_0_xyyz_yy_xz, t_0_xyyz_yz_xy,\
                           t_0_xyzz_x_yz, t_0_xyzz_x_yzz, t_0_xyzz_x_zz, t_0_xyzz_xy_zz,\
                           t_0_xyzz_xz_yz, t_0_xyzz_y_xz, t_0_xyzz_y_xzz, t_0_xyzz_yz_xz,\
                           t_0_xyzz_z_xy, t_0_xyzz_z_xyz, t_0_xyzz_zz_xy, t_0_xzzz_x_zz,\
                           t_0_xzzz_x_zzz, t_0_xzzz_xz_zz, t_0_xzzz_z_xz, t_0_xzzz_z_xzz,\
                           t_0_xzzz_zz_xz, t_0_yyyy_y_yy, t_0_yyyy_y_yyy, t_0_yyyy_yy_yy,\
                           t_0_yyyz_y_yy, t_0_yyyz_y_yyz, t_0_yyyz_y_yz, t_0_yyyz_yy_yz,\
                           t_0_yyyz_yz_yy, t_0_yyzz_y_yz, t_0_yyzz_y_yzz, t_0_yyzz_y_zz,\
                           t_0_yyzz_yy_zz, t_0_yyzz_yz_yz, t_0_yyzz_z_yy, t_0_yyzz_z_yyz,\
                           t_0_yyzz_zz_yy, t_0_yzzz_y_zz, t_0_yzzz_y_zzz, t_0_yzzz_yz_zz,\
                           t_0_yzzz_z_yz, t_0_yzzz_z_yzz, t_0_yzzz_zz_yz, t_0_zzzz_z_zz,\
                           t_0_zzzz_z_zzz, t_0_zzzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_zz_zz[i] = t_0_zzzz_z_zzz[i] - rcd_z[i] * t_0_zzzz_z_zz[i];

        t_0_yzzz_zz_yz[i] = t_0_yzzz_z_yzz[i] - rcd_z[i] * t_0_yzzz_z_yz[i];

        t_0_yzzz_yz_zz[i] = t_0_yzzz_y_zzz[i] - rcd_z[i] * t_0_yzzz_y_zz[i];

        t_0_yyzz_zz_yy[i] = t_0_yyzz_z_yyz[i] - rcd_z[i] * t_0_yyzz_z_yy[i];

        t_0_yyzz_yz_yz[i] = t_0_yyzz_y_yzz[i] - rcd_z[i] * t_0_yyzz_y_yz[i];

        t_0_yyzz_yy_zz[i] = t_0_yyzz_y_yzz[i] - rcd_y[i] * t_0_yyzz_y_zz[i];

        t_0_yyyz_yz_yy[i] = t_0_yyyz_y_yyz[i] - rcd_z[i] * t_0_yyyz_y_yy[i];

        t_0_yyyz_yy_yz[i] = t_0_yyyz_y_yyz[i] - rcd_y[i] * t_0_yyyz_y_yz[i];

        t_0_yyyy_yy_yy[i] = t_0_yyyy_y_yyy[i] - rcd_y[i] * t_0_yyyy_y_yy[i];

        t_0_xzzz_zz_xz[i] = t_0_xzzz_z_xzz[i] - rcd_z[i] * t_0_xzzz_z_xz[i];

        t_0_xzzz_xz_zz[i] = t_0_xzzz_x_zzz[i] - rcd_z[i] * t_0_xzzz_x_zz[i];

        t_0_xyzz_zz_xy[i] = t_0_xyzz_z_xyz[i] - rcd_z[i] * t_0_xyzz_z_xy[i];

        t_0_xyzz_yz_xz[i] = t_0_xyzz_y_xzz[i] - rcd_z[i] * t_0_xyzz_y_xz[i];

        t_0_xyzz_xz_yz[i] = t_0_xyzz_x_yzz[i] - rcd_z[i] * t_0_xyzz_x_yz[i];

        t_0_xyzz_xy_zz[i] = t_0_xyzz_x_yzz[i] - rcd_y[i] * t_0_xyzz_x_zz[i];

        t_0_xyyz_yz_xy[i] = t_0_xyyz_y_xyz[i] - rcd_z[i] * t_0_xyyz_y_xy[i];

        t_0_xyyz_yy_xz[i] = t_0_xyyz_y_xyz[i] - rcd_y[i] * t_0_xyyz_y_xz[i];

        t_0_xyyz_xz_yy[i] = t_0_xyyz_x_yyz[i] - rcd_z[i] * t_0_xyyz_x_yy[i];

        t_0_xyyz_xy_yz[i] = t_0_xyyz_x_yyz[i] - rcd_y[i] * t_0_xyyz_x_yz[i];

        t_0_xyyy_yy_xy[i] = t_0_xyyy_y_xyy[i] - rcd_y[i] * t_0_xyyy_y_xy[i];

        t_0_xyyy_xy_yy[i] = t_0_xyyy_x_yyy[i] - rcd_y[i] * t_0_xyyy_x_yy[i];

        t_0_xxzz_zz_xx[i] = t_0_xxzz_z_xxz[i] - rcd_z[i] * t_0_xxzz_z_xx[i];

        t_0_xxzz_xz_xz[i] = t_0_xxzz_x_xzz[i] - rcd_z[i] * t_0_xxzz_x_xz[i];

        t_0_xxzz_xx_zz[i] = t_0_xxzz_x_xzz[i] - rcd_x[i] * t_0_xxzz_x_zz[i];

        t_0_xxyz_yz_xx[i] = t_0_xxyz_y_xxz[i] - rcd_z[i] * t_0_xxyz_y_xx[i];

        t_0_xxyz_xz_xy[i] = t_0_xxyz_x_xyz[i] - rcd_z[i] * t_0_xxyz_x_xy[i];

        t_0_xxyz_xy_xz[i] = t_0_xxyz_x_xyz[i] - rcd_y[i] * t_0_xxyz_x_xz[i];

        t_0_xxyz_xx_yz[i] = t_0_xxyz_x_xyz[i] - rcd_x[i] * t_0_xxyz_x_yz[i];

        t_0_xxyy_yy_xx[i] = t_0_xxyy_y_xxy[i] - rcd_y[i] * t_0_xxyy_y_xx[i];

        t_0_xxyy_xy_xy[i] = t_0_xxyy_x_xyy[i] - rcd_y[i] * t_0_xxyy_x_xy[i];

        t_0_xxyy_xx_yy[i] = t_0_xxyy_x_xyy[i] - rcd_x[i] * t_0_xxyy_x_yy[i];

        t_0_xxxz_xz_xx[i] = t_0_xxxz_x_xxz[i] - rcd_z[i] * t_0_xxxz_x_xx[i];

        t_0_xxxz_xx_xz[i] = t_0_xxxz_x_xxz[i] - rcd_x[i] * t_0_xxxz_x_xz[i];

        t_0_xxxy_xy_xx[i] = t_0_xxxy_x_xxy[i] - rcd_y[i] * t_0_xxxy_x_xx[i];

        t_0_xxxy_xx_xy[i] = t_0_xxxy_x_xxy[i] - rcd_x[i] * t_0_xxxy_x_xy[i];

        t_0_xxxx_xx_xx[i] = t_0_xxxx_x_xxx[i] - rcd_x[i] * t_0_xxxx_x_xx[i];
    }
}


} // derirec namespace
