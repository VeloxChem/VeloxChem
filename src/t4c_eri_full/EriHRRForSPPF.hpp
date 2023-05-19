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
compHostHRRForSPPF_V0(      BufferHostXY<T>&      intsBufferSPPF,
                      const BufferHostX<int32_t>& intsIndexesSPPF,
                      const BufferHostXY<T>&      intsBufferSPSF,
                      const BufferHostX<int32_t>& intsIndexesSPSF,
                      const BufferHostXY<T>&      intsBufferSPSG,
                      const BufferHostX<int32_t>& intsIndexesSPSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SPPF) integral components

    t_0_z_z_zzz = intsBufferSPPF.data(intsIndexesSPPF(0));

    t_0_z_z_yzz = intsBufferSPPF.data(intsIndexesSPPF(1));

    t_0_z_z_yyz = intsBufferSPPF.data(intsIndexesSPPF(2));

    t_0_z_z_xzz = intsBufferSPPF.data(intsIndexesSPPF(3));

    t_0_z_z_xyz = intsBufferSPPF.data(intsIndexesSPPF(4));

    t_0_z_z_xxz = intsBufferSPPF.data(intsIndexesSPPF(5));

    t_0_z_y_zzz = intsBufferSPPF.data(intsIndexesSPPF(6));

    t_0_z_y_yzz = intsBufferSPPF.data(intsIndexesSPPF(7));

    t_0_z_y_yyz = intsBufferSPPF.data(intsIndexesSPPF(8));

    t_0_z_y_yyy = intsBufferSPPF.data(intsIndexesSPPF(9));

    t_0_z_y_xzz = intsBufferSPPF.data(intsIndexesSPPF(10));

    t_0_z_y_xyz = intsBufferSPPF.data(intsIndexesSPPF(11));

    t_0_z_y_xyy = intsBufferSPPF.data(intsIndexesSPPF(12));

    t_0_z_y_xxz = intsBufferSPPF.data(intsIndexesSPPF(13));

    t_0_z_y_xxy = intsBufferSPPF.data(intsIndexesSPPF(14));

    t_0_z_x_zzz = intsBufferSPPF.data(intsIndexesSPPF(15));

    t_0_z_x_yzz = intsBufferSPPF.data(intsIndexesSPPF(16));

    t_0_z_x_yyz = intsBufferSPPF.data(intsIndexesSPPF(17));

    t_0_z_x_yyy = intsBufferSPPF.data(intsIndexesSPPF(18));

    t_0_z_x_xzz = intsBufferSPPF.data(intsIndexesSPPF(19));

    t_0_z_x_xyz = intsBufferSPPF.data(intsIndexesSPPF(20));

    t_0_z_x_xyy = intsBufferSPPF.data(intsIndexesSPPF(21));

    t_0_z_x_xxz = intsBufferSPPF.data(intsIndexesSPPF(22));

    t_0_z_x_xxy = intsBufferSPPF.data(intsIndexesSPPF(23));

    t_0_z_x_xxx = intsBufferSPPF.data(intsIndexesSPPF(24));

    t_0_y_z_zzz = intsBufferSPPF.data(intsIndexesSPPF(25));

    t_0_y_z_yzz = intsBufferSPPF.data(intsIndexesSPPF(26));

    t_0_y_z_yyz = intsBufferSPPF.data(intsIndexesSPPF(27));

    t_0_y_z_xzz = intsBufferSPPF.data(intsIndexesSPPF(28));

    t_0_y_z_xyz = intsBufferSPPF.data(intsIndexesSPPF(29));

    t_0_y_z_xxz = intsBufferSPPF.data(intsIndexesSPPF(30));

    t_0_y_y_zzz = intsBufferSPPF.data(intsIndexesSPPF(31));

    t_0_y_y_yzz = intsBufferSPPF.data(intsIndexesSPPF(32));

    t_0_y_y_yyz = intsBufferSPPF.data(intsIndexesSPPF(33));

    t_0_y_y_yyy = intsBufferSPPF.data(intsIndexesSPPF(34));

    t_0_y_y_xzz = intsBufferSPPF.data(intsIndexesSPPF(35));

    t_0_y_y_xyz = intsBufferSPPF.data(intsIndexesSPPF(36));

    t_0_y_y_xyy = intsBufferSPPF.data(intsIndexesSPPF(37));

    t_0_y_y_xxz = intsBufferSPPF.data(intsIndexesSPPF(38));

    t_0_y_y_xxy = intsBufferSPPF.data(intsIndexesSPPF(39));

    t_0_y_x_zzz = intsBufferSPPF.data(intsIndexesSPPF(40));

    t_0_y_x_yzz = intsBufferSPPF.data(intsIndexesSPPF(41));

    t_0_y_x_yyz = intsBufferSPPF.data(intsIndexesSPPF(42));

    t_0_y_x_yyy = intsBufferSPPF.data(intsIndexesSPPF(43));

    t_0_y_x_xzz = intsBufferSPPF.data(intsIndexesSPPF(44));

    t_0_y_x_xyz = intsBufferSPPF.data(intsIndexesSPPF(45));

    t_0_y_x_xyy = intsBufferSPPF.data(intsIndexesSPPF(46));

    t_0_y_x_xxz = intsBufferSPPF.data(intsIndexesSPPF(47));

    t_0_y_x_xxy = intsBufferSPPF.data(intsIndexesSPPF(48));

    t_0_y_x_xxx = intsBufferSPPF.data(intsIndexesSPPF(49));

    t_0_x_z_zzz = intsBufferSPPF.data(intsIndexesSPPF(50));

    t_0_x_z_yzz = intsBufferSPPF.data(intsIndexesSPPF(51));

    t_0_x_z_yyz = intsBufferSPPF.data(intsIndexesSPPF(52));

    t_0_x_z_xzz = intsBufferSPPF.data(intsIndexesSPPF(53));

    t_0_x_z_xyz = intsBufferSPPF.data(intsIndexesSPPF(54));

    t_0_x_z_xxz = intsBufferSPPF.data(intsIndexesSPPF(55));

    t_0_x_y_zzz = intsBufferSPPF.data(intsIndexesSPPF(56));

    t_0_x_y_yzz = intsBufferSPPF.data(intsIndexesSPPF(57));

    t_0_x_y_yyz = intsBufferSPPF.data(intsIndexesSPPF(58));

    t_0_x_y_yyy = intsBufferSPPF.data(intsIndexesSPPF(59));

    t_0_x_y_xzz = intsBufferSPPF.data(intsIndexesSPPF(60));

    t_0_x_y_xyz = intsBufferSPPF.data(intsIndexesSPPF(61));

    t_0_x_y_xyy = intsBufferSPPF.data(intsIndexesSPPF(62));

    t_0_x_y_xxz = intsBufferSPPF.data(intsIndexesSPPF(63));

    t_0_x_y_xxy = intsBufferSPPF.data(intsIndexesSPPF(64));

    t_0_x_x_zzz = intsBufferSPPF.data(intsIndexesSPPF(65));

    t_0_x_x_yzz = intsBufferSPPF.data(intsIndexesSPPF(66));

    t_0_x_x_yyz = intsBufferSPPF.data(intsIndexesSPPF(67));

    t_0_x_x_yyy = intsBufferSPPF.data(intsIndexesSPPF(68));

    t_0_x_x_xzz = intsBufferSPPF.data(intsIndexesSPPF(69));

    t_0_x_x_xyz = intsBufferSPPF.data(intsIndexesSPPF(70));

    t_0_x_x_xyy = intsBufferSPPF.data(intsIndexesSPPF(71));

    t_0_x_x_xxz = intsBufferSPPF.data(intsIndexesSPPF(72));

    t_0_x_x_xxy = intsBufferSPPF.data(intsIndexesSPPF(73));

    t_0_x_x_xxx = intsBufferSPPF.data(intsIndexesSPPF(74));

    // set up (SPSF) integral components

    t_0_z_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(0));

    t_0_z_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(1));

    t_0_z_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(2));

    t_0_z_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(3));

    t_0_z_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(4));

    t_0_z_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(5));

    t_0_z_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(6));

    t_0_z_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(7));

    t_0_z_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(8));

    t_0_z_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(9));

    t_0_y_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(10));

    t_0_y_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(11));

    t_0_y_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(12));

    t_0_y_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(13));

    t_0_y_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(14));

    t_0_y_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(15));

    t_0_y_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(16));

    t_0_y_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(17));

    t_0_y_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(18));

    t_0_y_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(19));

    t_0_x_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(20));

    t_0_x_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(21));

    t_0_x_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(22));

    t_0_x_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(23));

    t_0_x_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(24));

    t_0_x_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(25));

    t_0_x_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(26));

    t_0_x_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(27));

    t_0_x_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(28));

    t_0_x_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(29));

    // set up (SPSG) integral components

    t_0_z_0_zzzz = intsBufferSPSG.data(intsIndexesSPSG(0));

    t_0_z_0_yzzz = intsBufferSPSG.data(intsIndexesSPSG(1));

    t_0_z_0_yyzz = intsBufferSPSG.data(intsIndexesSPSG(2));

    t_0_z_0_yyyz = intsBufferSPSG.data(intsIndexesSPSG(3));

    t_0_z_0_yyyy = intsBufferSPSG.data(intsIndexesSPSG(4));

    t_0_z_0_xzzz = intsBufferSPSG.data(intsIndexesSPSG(5));

    t_0_z_0_xyzz = intsBufferSPSG.data(intsIndexesSPSG(6));

    t_0_z_0_xyyz = intsBufferSPSG.data(intsIndexesSPSG(7));

    t_0_z_0_xyyy = intsBufferSPSG.data(intsIndexesSPSG(8));

    t_0_z_0_xxzz = intsBufferSPSG.data(intsIndexesSPSG(9));

    t_0_z_0_xxyz = intsBufferSPSG.data(intsIndexesSPSG(10));

    t_0_z_0_xxyy = intsBufferSPSG.data(intsIndexesSPSG(11));

    t_0_z_0_xxxz = intsBufferSPSG.data(intsIndexesSPSG(12));

    t_0_z_0_xxxy = intsBufferSPSG.data(intsIndexesSPSG(13));

    t_0_z_0_xxxx = intsBufferSPSG.data(intsIndexesSPSG(14));

    t_0_y_0_zzzz = intsBufferSPSG.data(intsIndexesSPSG(15));

    t_0_y_0_yzzz = intsBufferSPSG.data(intsIndexesSPSG(16));

    t_0_y_0_yyzz = intsBufferSPSG.data(intsIndexesSPSG(17));

    t_0_y_0_yyyz = intsBufferSPSG.data(intsIndexesSPSG(18));

    t_0_y_0_yyyy = intsBufferSPSG.data(intsIndexesSPSG(19));

    t_0_y_0_xzzz = intsBufferSPSG.data(intsIndexesSPSG(20));

    t_0_y_0_xyzz = intsBufferSPSG.data(intsIndexesSPSG(21));

    t_0_y_0_xyyz = intsBufferSPSG.data(intsIndexesSPSG(22));

    t_0_y_0_xyyy = intsBufferSPSG.data(intsIndexesSPSG(23));

    t_0_y_0_xxzz = intsBufferSPSG.data(intsIndexesSPSG(24));

    t_0_y_0_xxyz = intsBufferSPSG.data(intsIndexesSPSG(25));

    t_0_y_0_xxyy = intsBufferSPSG.data(intsIndexesSPSG(26));

    t_0_y_0_xxxz = intsBufferSPSG.data(intsIndexesSPSG(27));

    t_0_y_0_xxxy = intsBufferSPSG.data(intsIndexesSPSG(28));

    t_0_y_0_xxxx = intsBufferSPSG.data(intsIndexesSPSG(29));

    t_0_x_0_zzzz = intsBufferSPSG.data(intsIndexesSPSG(30));

    t_0_x_0_yzzz = intsBufferSPSG.data(intsIndexesSPSG(31));

    t_0_x_0_yyzz = intsBufferSPSG.data(intsIndexesSPSG(32));

    t_0_x_0_yyyz = intsBufferSPSG.data(intsIndexesSPSG(33));

    t_0_x_0_yyyy = intsBufferSPSG.data(intsIndexesSPSG(34));

    t_0_x_0_xzzz = intsBufferSPSG.data(intsIndexesSPSG(35));

    t_0_x_0_xyzz = intsBufferSPSG.data(intsIndexesSPSG(36));

    t_0_x_0_xyyz = intsBufferSPSG.data(intsIndexesSPSG(37));

    t_0_x_0_xyyy = intsBufferSPSG.data(intsIndexesSPSG(38));

    t_0_x_0_xxzz = intsBufferSPSG.data(intsIndexesSPSG(39));

    t_0_x_0_xxyz = intsBufferSPSG.data(intsIndexesSPSG(40));

    t_0_x_0_xxyy = intsBufferSPSG.data(intsIndexesSPSG(41));

    t_0_x_0_xxxz = intsBufferSPSG.data(intsIndexesSPSG(42));

    t_0_x_0_xxxy = intsBufferSPSG.data(intsIndexesSPSG(43));

    t_0_x_0_xxxx = intsBufferSPSG.data(intsIndexesSPSG(44));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_y_0_xxz, t_0_y_0_xxzz, t_0_y_0_xyz,\
                           t_0_y_0_xyzz, t_0_y_0_xzz, t_0_y_0_xzzz, t_0_y_0_yyy, t_0_y_0_yyyy,\
                           t_0_y_0_yyyz, t_0_y_0_yyz, t_0_y_0_yyzz, t_0_y_0_yzz, t_0_y_0_yzzz,\
                           t_0_y_0_zzz, t_0_y_0_zzzz, t_0_y_y_xzz, t_0_y_y_yyy, t_0_y_y_yyz,\
                           t_0_y_y_yzz, t_0_y_y_zzz, t_0_y_z_xxz, t_0_y_z_xyz, t_0_y_z_xzz,\
                           t_0_y_z_yyz, t_0_y_z_yzz, t_0_y_z_zzz, t_0_z_0_xxx, t_0_z_0_xxxx,\
                           t_0_z_0_xxxy, t_0_z_0_xxxz, t_0_z_0_xxy, t_0_z_0_xxyy, t_0_z_0_xxyz,\
                           t_0_z_0_xxz, t_0_z_0_xxzz, t_0_z_0_xyy, t_0_z_0_xyyy, t_0_z_0_xyyz,\
                           t_0_z_0_xyz, t_0_z_0_xyzz, t_0_z_0_xzz, t_0_z_0_xzzz, t_0_z_0_yyy,\
                           t_0_z_0_yyyy, t_0_z_0_yyyz, t_0_z_0_yyz, t_0_z_0_yyzz, t_0_z_0_yzz,\
                           t_0_z_0_yzzz, t_0_z_0_zzz, t_0_z_0_zzzz, t_0_z_x_xxx, t_0_z_x_xxy,\
                           t_0_z_x_xxz, t_0_z_x_xyy, t_0_z_x_xyz, t_0_z_x_xzz, t_0_z_x_yyy,\
                           t_0_z_x_yyz, t_0_z_x_yzz, t_0_z_x_zzz, t_0_z_y_xxy, t_0_z_y_xxz,\
                           t_0_z_y_xyy, t_0_z_y_xyz, t_0_z_y_xzz, t_0_z_y_yyy, t_0_z_y_yyz,\
                           t_0_z_y_yzz, t_0_z_y_zzz, t_0_z_z_xxz, t_0_z_z_xyz, t_0_z_z_xzz,\
                           t_0_z_z_yyz, t_0_z_z_yzz, t_0_z_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_z_z_zzz[i] = t_0_z_0_zzzz[i] - rcd_z[i] * t_0_z_0_zzz[i];

        t_0_z_z_yzz[i] = t_0_z_0_yzzz[i] - rcd_z[i] * t_0_z_0_yzz[i];

        t_0_z_z_yyz[i] = t_0_z_0_yyzz[i] - rcd_z[i] * t_0_z_0_yyz[i];

        t_0_z_z_xzz[i] = t_0_z_0_xzzz[i] - rcd_z[i] * t_0_z_0_xzz[i];

        t_0_z_z_xyz[i] = t_0_z_0_xyzz[i] - rcd_z[i] * t_0_z_0_xyz[i];

        t_0_z_z_xxz[i] = t_0_z_0_xxzz[i] - rcd_z[i] * t_0_z_0_xxz[i];

        t_0_z_y_zzz[i] = t_0_z_0_yzzz[i] - rcd_y[i] * t_0_z_0_zzz[i];

        t_0_z_y_yzz[i] = t_0_z_0_yyzz[i] - rcd_y[i] * t_0_z_0_yzz[i];

        t_0_z_y_yyz[i] = t_0_z_0_yyyz[i] - rcd_y[i] * t_0_z_0_yyz[i];

        t_0_z_y_yyy[i] = t_0_z_0_yyyy[i] - rcd_y[i] * t_0_z_0_yyy[i];

        t_0_z_y_xzz[i] = t_0_z_0_xyzz[i] - rcd_y[i] * t_0_z_0_xzz[i];

        t_0_z_y_xyz[i] = t_0_z_0_xyyz[i] - rcd_y[i] * t_0_z_0_xyz[i];

        t_0_z_y_xyy[i] = t_0_z_0_xyyy[i] - rcd_y[i] * t_0_z_0_xyy[i];

        t_0_z_y_xxz[i] = t_0_z_0_xxyz[i] - rcd_y[i] * t_0_z_0_xxz[i];

        t_0_z_y_xxy[i] = t_0_z_0_xxyy[i] - rcd_y[i] * t_0_z_0_xxy[i];

        t_0_z_x_zzz[i] = t_0_z_0_xzzz[i] - rcd_x[i] * t_0_z_0_zzz[i];

        t_0_z_x_yzz[i] = t_0_z_0_xyzz[i] - rcd_x[i] * t_0_z_0_yzz[i];

        t_0_z_x_yyz[i] = t_0_z_0_xyyz[i] - rcd_x[i] * t_0_z_0_yyz[i];

        t_0_z_x_yyy[i] = t_0_z_0_xyyy[i] - rcd_x[i] * t_0_z_0_yyy[i];

        t_0_z_x_xzz[i] = t_0_z_0_xxzz[i] - rcd_x[i] * t_0_z_0_xzz[i];

        t_0_z_x_xyz[i] = t_0_z_0_xxyz[i] - rcd_x[i] * t_0_z_0_xyz[i];

        t_0_z_x_xyy[i] = t_0_z_0_xxyy[i] - rcd_x[i] * t_0_z_0_xyy[i];

        t_0_z_x_xxz[i] = t_0_z_0_xxxz[i] - rcd_x[i] * t_0_z_0_xxz[i];

        t_0_z_x_xxy[i] = t_0_z_0_xxxy[i] - rcd_x[i] * t_0_z_0_xxy[i];

        t_0_z_x_xxx[i] = t_0_z_0_xxxx[i] - rcd_x[i] * t_0_z_0_xxx[i];

        t_0_y_z_zzz[i] = t_0_y_0_zzzz[i] - rcd_z[i] * t_0_y_0_zzz[i];

        t_0_y_z_yzz[i] = t_0_y_0_yzzz[i] - rcd_z[i] * t_0_y_0_yzz[i];

        t_0_y_z_yyz[i] = t_0_y_0_yyzz[i] - rcd_z[i] * t_0_y_0_yyz[i];

        t_0_y_z_xzz[i] = t_0_y_0_xzzz[i] - rcd_z[i] * t_0_y_0_xzz[i];

        t_0_y_z_xyz[i] = t_0_y_0_xyzz[i] - rcd_z[i] * t_0_y_0_xyz[i];

        t_0_y_z_xxz[i] = t_0_y_0_xxzz[i] - rcd_z[i] * t_0_y_0_xxz[i];

        t_0_y_y_zzz[i] = t_0_y_0_yzzz[i] - rcd_y[i] * t_0_y_0_zzz[i];

        t_0_y_y_yzz[i] = t_0_y_0_yyzz[i] - rcd_y[i] * t_0_y_0_yzz[i];

        t_0_y_y_yyz[i] = t_0_y_0_yyyz[i] - rcd_y[i] * t_0_y_0_yyz[i];

        t_0_y_y_yyy[i] = t_0_y_0_yyyy[i] - rcd_y[i] * t_0_y_0_yyy[i];

        t_0_y_y_xzz[i] = t_0_y_0_xyzz[i] - rcd_y[i] * t_0_y_0_xzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_x_0_xxy, t_0_x_0_xxyy, t_0_x_0_xxyz,\
                           t_0_x_0_xxz, t_0_x_0_xxzz, t_0_x_0_xyy, t_0_x_0_xyyy, t_0_x_0_xyyz,\
                           t_0_x_0_xyz, t_0_x_0_xyzz, t_0_x_0_xzz, t_0_x_0_xzzz, t_0_x_0_yyy,\
                           t_0_x_0_yyyy, t_0_x_0_yyyz, t_0_x_0_yyz, t_0_x_0_yyzz, t_0_x_0_yzz,\
                           t_0_x_0_yzzz, t_0_x_0_zzz, t_0_x_0_zzzz, t_0_x_x_xyy, t_0_x_x_xyz,\
                           t_0_x_x_xzz, t_0_x_x_yyy, t_0_x_x_yyz, t_0_x_x_yzz, t_0_x_x_zzz,\
                           t_0_x_y_xxy, t_0_x_y_xxz, t_0_x_y_xyy, t_0_x_y_xyz, t_0_x_y_xzz,\
                           t_0_x_y_yyy, t_0_x_y_yyz, t_0_x_y_yzz, t_0_x_y_zzz, t_0_x_z_xxz,\
                           t_0_x_z_xyz, t_0_x_z_xzz, t_0_x_z_yyz, t_0_x_z_yzz, t_0_x_z_zzz,\
                           t_0_y_0_xxx, t_0_y_0_xxxx, t_0_y_0_xxxy, t_0_y_0_xxxz, t_0_y_0_xxy,\
                           t_0_y_0_xxyy, t_0_y_0_xxyz, t_0_y_0_xxz, t_0_y_0_xxzz, t_0_y_0_xyy,\
                           t_0_y_0_xyyy, t_0_y_0_xyyz, t_0_y_0_xyz, t_0_y_0_xyzz, t_0_y_0_xzz,\
                           t_0_y_0_xzzz, t_0_y_0_yyy, t_0_y_0_yyz, t_0_y_0_yzz, t_0_y_0_zzz,\
                           t_0_y_x_xxx, t_0_y_x_xxy, t_0_y_x_xxz, t_0_y_x_xyy, t_0_y_x_xyz,\
                           t_0_y_x_xzz, t_0_y_x_yyy, t_0_y_x_yyz, t_0_y_x_yzz, t_0_y_x_zzz,\
                           t_0_y_y_xxy, t_0_y_y_xxz, t_0_y_y_xyy, t_0_y_y_xyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_y_y_xyz[i] = t_0_y_0_xyyz[i] - rcd_y[i] * t_0_y_0_xyz[i];

        t_0_y_y_xyy[i] = t_0_y_0_xyyy[i] - rcd_y[i] * t_0_y_0_xyy[i];

        t_0_y_y_xxz[i] = t_0_y_0_xxyz[i] - rcd_y[i] * t_0_y_0_xxz[i];

        t_0_y_y_xxy[i] = t_0_y_0_xxyy[i] - rcd_y[i] * t_0_y_0_xxy[i];

        t_0_y_x_zzz[i] = t_0_y_0_xzzz[i] - rcd_x[i] * t_0_y_0_zzz[i];

        t_0_y_x_yzz[i] = t_0_y_0_xyzz[i] - rcd_x[i] * t_0_y_0_yzz[i];

        t_0_y_x_yyz[i] = t_0_y_0_xyyz[i] - rcd_x[i] * t_0_y_0_yyz[i];

        t_0_y_x_yyy[i] = t_0_y_0_xyyy[i] - rcd_x[i] * t_0_y_0_yyy[i];

        t_0_y_x_xzz[i] = t_0_y_0_xxzz[i] - rcd_x[i] * t_0_y_0_xzz[i];

        t_0_y_x_xyz[i] = t_0_y_0_xxyz[i] - rcd_x[i] * t_0_y_0_xyz[i];

        t_0_y_x_xyy[i] = t_0_y_0_xxyy[i] - rcd_x[i] * t_0_y_0_xyy[i];

        t_0_y_x_xxz[i] = t_0_y_0_xxxz[i] - rcd_x[i] * t_0_y_0_xxz[i];

        t_0_y_x_xxy[i] = t_0_y_0_xxxy[i] - rcd_x[i] * t_0_y_0_xxy[i];

        t_0_y_x_xxx[i] = t_0_y_0_xxxx[i] - rcd_x[i] * t_0_y_0_xxx[i];

        t_0_x_z_zzz[i] = t_0_x_0_zzzz[i] - rcd_z[i] * t_0_x_0_zzz[i];

        t_0_x_z_yzz[i] = t_0_x_0_yzzz[i] - rcd_z[i] * t_0_x_0_yzz[i];

        t_0_x_z_yyz[i] = t_0_x_0_yyzz[i] - rcd_z[i] * t_0_x_0_yyz[i];

        t_0_x_z_xzz[i] = t_0_x_0_xzzz[i] - rcd_z[i] * t_0_x_0_xzz[i];

        t_0_x_z_xyz[i] = t_0_x_0_xyzz[i] - rcd_z[i] * t_0_x_0_xyz[i];

        t_0_x_z_xxz[i] = t_0_x_0_xxzz[i] - rcd_z[i] * t_0_x_0_xxz[i];

        t_0_x_y_zzz[i] = t_0_x_0_yzzz[i] - rcd_y[i] * t_0_x_0_zzz[i];

        t_0_x_y_yzz[i] = t_0_x_0_yyzz[i] - rcd_y[i] * t_0_x_0_yzz[i];

        t_0_x_y_yyz[i] = t_0_x_0_yyyz[i] - rcd_y[i] * t_0_x_0_yyz[i];

        t_0_x_y_yyy[i] = t_0_x_0_yyyy[i] - rcd_y[i] * t_0_x_0_yyy[i];

        t_0_x_y_xzz[i] = t_0_x_0_xyzz[i] - rcd_y[i] * t_0_x_0_xzz[i];

        t_0_x_y_xyz[i] = t_0_x_0_xyyz[i] - rcd_y[i] * t_0_x_0_xyz[i];

        t_0_x_y_xyy[i] = t_0_x_0_xyyy[i] - rcd_y[i] * t_0_x_0_xyy[i];

        t_0_x_y_xxz[i] = t_0_x_0_xxyz[i] - rcd_y[i] * t_0_x_0_xxz[i];

        t_0_x_y_xxy[i] = t_0_x_0_xxyy[i] - rcd_y[i] * t_0_x_0_xxy[i];

        t_0_x_x_zzz[i] = t_0_x_0_xzzz[i] - rcd_x[i] * t_0_x_0_zzz[i];

        t_0_x_x_yzz[i] = t_0_x_0_xyzz[i] - rcd_x[i] * t_0_x_0_yzz[i];

        t_0_x_x_yyz[i] = t_0_x_0_xyyz[i] - rcd_x[i] * t_0_x_0_yyz[i];

        t_0_x_x_yyy[i] = t_0_x_0_xyyy[i] - rcd_x[i] * t_0_x_0_yyy[i];

        t_0_x_x_xzz[i] = t_0_x_0_xxzz[i] - rcd_x[i] * t_0_x_0_xzz[i];

        t_0_x_x_xyz[i] = t_0_x_0_xxyz[i] - rcd_x[i] * t_0_x_0_xyz[i];

        t_0_x_x_xyy[i] = t_0_x_0_xxyy[i] - rcd_x[i] * t_0_x_0_xyy[i];
    }

    #pragma omp simd align(rcd_x, t_0_x_0_xxx, t_0_x_0_xxxx, t_0_x_0_xxxy, t_0_x_0_xxxz,\
                           t_0_x_0_xxy, t_0_x_0_xxz, t_0_x_x_xxx, t_0_x_x_xxy, t_0_x_x_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_x_x_xxz[i] = t_0_x_0_xxxz[i] - rcd_x[i] * t_0_x_0_xxz[i];

        t_0_x_x_xxy[i] = t_0_x_0_xxxy[i] - rcd_x[i] * t_0_x_0_xxy[i];

        t_0_x_x_xxx[i] = t_0_x_0_xxxx[i] - rcd_x[i] * t_0_x_0_xxx[i];
    }
}


} // derirec namespace
