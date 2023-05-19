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
compHostHRRForSGPD_V0(      BufferHostXY<T>&      intsBufferSGPD,
                      const BufferHostX<int32_t>& intsIndexesSGPD,
                      const BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD,
                      const BufferHostXY<T>&      intsBufferSGSF,
                      const BufferHostX<int32_t>& intsIndexesSGSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SGSD) integral components

    t_0_zzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(0));

    t_0_yzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(1));

    t_0_yzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(2));

    t_0_yyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(3));

    t_0_yyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(4));

    t_0_yyzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(5));

    t_0_yyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(6));

    t_0_yyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(7));

    t_0_yyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(8));

    t_0_xzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(9));

    t_0_xzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(10));

    t_0_xyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(11));

    t_0_xyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(12));

    t_0_xyzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(13));

    t_0_xyzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(14));

    t_0_xyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(15));

    t_0_xyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(16));

    t_0_xyyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(17));

    t_0_xyyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(18));

    t_0_xyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(19));

    t_0_xyyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(20));

    t_0_xxzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(21));

    t_0_xxzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(22));

    t_0_xxzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(23));

    t_0_xxyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(24));

    t_0_xxyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(25));

    t_0_xxyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(26));

    t_0_xxyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(27));

    t_0_xxyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(28));

    t_0_xxyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(29));

    t_0_xxyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(30));

    t_0_xxxz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(31));

    t_0_xxxz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(32));

    t_0_xxxy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(33));

    t_0_xxxy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(34));

    t_0_xxxx_0_xx = intsBufferSGSD.data(intsIndexesSGSD(35));

    // set up (SGSF) integral components

    t_0_zzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(0));

    t_0_yzzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(1));

    t_0_yyzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(2));

    t_0_yyzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(3));

    t_0_yyyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(4));

    t_0_yyyz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(5));

    t_0_yyyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(6));

    t_0_xzzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(7));

    t_0_xyzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(8));

    t_0_xyzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(9));

    t_0_xyyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(10));

    t_0_xyyz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(11));

    t_0_xyyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(12));

    t_0_xxzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(13));

    t_0_xxzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(14));

    t_0_xxyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(15));

    t_0_xxyz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(16));

    t_0_xxyz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(17));

    t_0_xxyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(18));

    t_0_xxyy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(19));

    t_0_xxxz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(20));

    t_0_xxxz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(21));

    t_0_xxxy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(22));

    t_0_xxxy_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(23));

    t_0_xxxx_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(24));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_0_xx, t_0_xxxx_0_xxx, t_0_xxxx_x_xx,\
                           t_0_xxxy_0_xx, t_0_xxxy_0_xxx, t_0_xxxy_0_xxy, t_0_xxxy_0_xy,\
                           t_0_xxxy_x_xx, t_0_xxxy_x_xy, t_0_xxxz_0_xx, t_0_xxxz_0_xxx,\
                           t_0_xxxz_0_xxz, t_0_xxxz_0_xz, t_0_xxxz_x_xx, t_0_xxxz_x_xz,\
                           t_0_xxyy_0_xx, t_0_xxyy_0_xxy, t_0_xxyy_0_xy, t_0_xxyy_0_xyy,\
                           t_0_xxyy_0_yy, t_0_xxyy_x_xy, t_0_xxyy_x_yy, t_0_xxyy_y_xx,\
                           t_0_xxyz_0_xx, t_0_xxyz_0_xxy, t_0_xxyz_0_xxz, t_0_xxyz_0_xy,\
                           t_0_xxyz_0_xyz, t_0_xxyz_0_xz, t_0_xxyz_0_yz, t_0_xxyz_x_xy,\
                           t_0_xxyz_x_xz, t_0_xxyz_x_yz, t_0_xxyz_y_xx, t_0_xxzz_0_xx,\
                           t_0_xxzz_0_xxz, t_0_xxzz_0_xz, t_0_xxzz_0_xzz, t_0_xxzz_0_zz,\
                           t_0_xxzz_x_xz, t_0_xxzz_x_zz, t_0_xxzz_z_xx, t_0_xyyy_0_xy,\
                           t_0_xyyy_0_xyy, t_0_xyyy_0_yy, t_0_xyyy_x_yy, t_0_xyyy_y_xy,\
                           t_0_xyyz_0_xy, t_0_xyyz_0_xyy, t_0_xyyz_0_xyz, t_0_xyyz_0_xz,\
                           t_0_xyyz_0_yy, t_0_xyyz_0_yz, t_0_xyyz_x_yy, t_0_xyyz_x_yz,\
                           t_0_xyyz_y_xy, t_0_xyyz_y_xz, t_0_xyzz_0_xy, t_0_xyzz_0_xyz,\
                           t_0_xyzz_0_xz, t_0_xyzz_0_xzz, t_0_xyzz_0_yz, t_0_xyzz_0_zz,\
                           t_0_xyzz_x_yz, t_0_xyzz_x_zz, t_0_xyzz_y_xz, t_0_xyzz_z_xy,\
                           t_0_xzzz_0_xz, t_0_xzzz_0_xzz, t_0_xzzz_0_zz, t_0_xzzz_x_zz,\
                           t_0_xzzz_z_xz, t_0_yyyy_0_yy, t_0_yyyy_0_yyy, t_0_yyyy_y_yy,\
                           t_0_yyyz_0_yy, t_0_yyyz_0_yyy, t_0_yyyz_0_yyz, t_0_yyyz_0_yz,\
                           t_0_yyyz_y_yy, t_0_yyyz_y_yz, t_0_yyzz_0_yy, t_0_yyzz_0_yyz,\
                           t_0_yyzz_0_yz, t_0_yyzz_0_yzz, t_0_yyzz_0_zz, t_0_yyzz_y_yz,\
                           t_0_yyzz_y_zz, t_0_yyzz_z_yy, t_0_yzzz_0_yz, t_0_yzzz_0_yzz,\
                           t_0_yzzz_0_zz, t_0_yzzz_y_zz, t_0_yzzz_z_yz, t_0_zzzz_0_zz,\
                           t_0_zzzz_0_zzz, t_0_zzzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_z_zz[i] = t_0_zzzz_0_zzz[i] - rcd_z[i] * t_0_zzzz_0_zz[i];

        t_0_yzzz_z_yz[i] = t_0_yzzz_0_yzz[i] - rcd_z[i] * t_0_yzzz_0_yz[i];

        t_0_yzzz_y_zz[i] = t_0_yzzz_0_yzz[i] - rcd_y[i] * t_0_yzzz_0_zz[i];

        t_0_yyzz_z_yy[i] = t_0_yyzz_0_yyz[i] - rcd_z[i] * t_0_yyzz_0_yy[i];

        t_0_yyzz_y_zz[i] = t_0_yyzz_0_yzz[i] - rcd_y[i] * t_0_yyzz_0_zz[i];

        t_0_yyzz_y_yz[i] = t_0_yyzz_0_yyz[i] - rcd_y[i] * t_0_yyzz_0_yz[i];

        t_0_yyyz_y_yz[i] = t_0_yyyz_0_yyz[i] - rcd_y[i] * t_0_yyyz_0_yz[i];

        t_0_yyyz_y_yy[i] = t_0_yyyz_0_yyy[i] - rcd_y[i] * t_0_yyyz_0_yy[i];

        t_0_yyyy_y_yy[i] = t_0_yyyy_0_yyy[i] - rcd_y[i] * t_0_yyyy_0_yy[i];

        t_0_xzzz_z_xz[i] = t_0_xzzz_0_xzz[i] - rcd_z[i] * t_0_xzzz_0_xz[i];

        t_0_xzzz_x_zz[i] = t_0_xzzz_0_xzz[i] - rcd_x[i] * t_0_xzzz_0_zz[i];

        t_0_xyzz_z_xy[i] = t_0_xyzz_0_xyz[i] - rcd_z[i] * t_0_xyzz_0_xy[i];

        t_0_xyzz_y_xz[i] = t_0_xyzz_0_xyz[i] - rcd_y[i] * t_0_xyzz_0_xz[i];

        t_0_xyzz_x_zz[i] = t_0_xyzz_0_xzz[i] - rcd_x[i] * t_0_xyzz_0_zz[i];

        t_0_xyzz_x_yz[i] = t_0_xyzz_0_xyz[i] - rcd_x[i] * t_0_xyzz_0_yz[i];

        t_0_xyyz_y_xz[i] = t_0_xyyz_0_xyz[i] - rcd_y[i] * t_0_xyyz_0_xz[i];

        t_0_xyyz_y_xy[i] = t_0_xyyz_0_xyy[i] - rcd_y[i] * t_0_xyyz_0_xy[i];

        t_0_xyyz_x_yz[i] = t_0_xyyz_0_xyz[i] - rcd_x[i] * t_0_xyyz_0_yz[i];

        t_0_xyyz_x_yy[i] = t_0_xyyz_0_xyy[i] - rcd_x[i] * t_0_xyyz_0_yy[i];

        t_0_xyyy_y_xy[i] = t_0_xyyy_0_xyy[i] - rcd_y[i] * t_0_xyyy_0_xy[i];

        t_0_xyyy_x_yy[i] = t_0_xyyy_0_xyy[i] - rcd_x[i] * t_0_xyyy_0_yy[i];

        t_0_xxzz_z_xx[i] = t_0_xxzz_0_xxz[i] - rcd_z[i] * t_0_xxzz_0_xx[i];

        t_0_xxzz_x_zz[i] = t_0_xxzz_0_xzz[i] - rcd_x[i] * t_0_xxzz_0_zz[i];

        t_0_xxzz_x_xz[i] = t_0_xxzz_0_xxz[i] - rcd_x[i] * t_0_xxzz_0_xz[i];

        t_0_xxyz_y_xx[i] = t_0_xxyz_0_xxy[i] - rcd_y[i] * t_0_xxyz_0_xx[i];

        t_0_xxyz_x_yz[i] = t_0_xxyz_0_xyz[i] - rcd_x[i] * t_0_xxyz_0_yz[i];

        t_0_xxyz_x_xz[i] = t_0_xxyz_0_xxz[i] - rcd_x[i] * t_0_xxyz_0_xz[i];

        t_0_xxyz_x_xy[i] = t_0_xxyz_0_xxy[i] - rcd_x[i] * t_0_xxyz_0_xy[i];

        t_0_xxyy_y_xx[i] = t_0_xxyy_0_xxy[i] - rcd_y[i] * t_0_xxyy_0_xx[i];

        t_0_xxyy_x_yy[i] = t_0_xxyy_0_xyy[i] - rcd_x[i] * t_0_xxyy_0_yy[i];

        t_0_xxyy_x_xy[i] = t_0_xxyy_0_xxy[i] - rcd_x[i] * t_0_xxyy_0_xy[i];

        t_0_xxxz_x_xz[i] = t_0_xxxz_0_xxz[i] - rcd_x[i] * t_0_xxxz_0_xz[i];

        t_0_xxxz_x_xx[i] = t_0_xxxz_0_xxx[i] - rcd_x[i] * t_0_xxxz_0_xx[i];

        t_0_xxxy_x_xy[i] = t_0_xxxy_0_xxy[i] - rcd_x[i] * t_0_xxxy_0_xy[i];

        t_0_xxxy_x_xx[i] = t_0_xxxy_0_xxx[i] - rcd_x[i] * t_0_xxxy_0_xx[i];

        t_0_xxxx_x_xx[i] = t_0_xxxx_0_xxx[i] - rcd_x[i] * t_0_xxxx_0_xx[i];
    }
}


} // derirec namespace
