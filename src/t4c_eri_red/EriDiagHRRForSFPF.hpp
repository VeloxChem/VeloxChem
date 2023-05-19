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
compHostHRRForSFPF_V0(      BufferHostXY<T>&      intsBufferSFPF,
                      const BufferHostX<int32_t>& intsIndexesSFPF,
                      const BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF,
                      const BufferHostXY<T>&      intsBufferSFSG,
                      const BufferHostX<int32_t>& intsIndexesSFSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPF) integral components

    t_0_zzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(0));

    t_0_zzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(1));

    t_0_zzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(2));

    t_0_yzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(3));

    t_0_yzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(4));

    t_0_yzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(5));

    t_0_yzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(6));

    t_0_yyz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(7));

    t_0_yyz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(8));

    t_0_yyz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(9));

    t_0_yyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(10));

    t_0_yyy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(11));

    t_0_yyy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(12));

    t_0_yyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(13));

    t_0_xzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(14));

    t_0_xzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(15));

    t_0_xzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(16));

    t_0_xzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(17));

    t_0_xzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(18));

    t_0_xyz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(19));

    t_0_xyz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(20));

    t_0_xyz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(21));

    t_0_xyz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(22));

    t_0_xyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(23));

    t_0_xyz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(24));

    t_0_xyy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(25));

    t_0_xyy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(26));

    t_0_xyy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(27));

    t_0_xyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(28));

    t_0_xyy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(29));

    t_0_xxz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(30));

    t_0_xxz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(31));

    t_0_xxz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(32));

    t_0_xxz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(33));

    t_0_xxz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(34));

    t_0_xxy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(35));

    t_0_xxy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(36));

    t_0_xxy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(37));

    t_0_xxy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(38));

    t_0_xxy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(39));

    t_0_xxx_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(40));

    t_0_xxx_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(41));

    t_0_xxx_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(42));

    // set up (SFSF) integral components

    t_0_zzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(0));

    t_0_yzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(1));

    t_0_yzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(2));

    t_0_yyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(3));

    t_0_yyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(4));

    t_0_yyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(5));

    t_0_yyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(6));

    t_0_xzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(7));

    t_0_xzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(8));

    t_0_xzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(9));

    t_0_xyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(10));

    t_0_xyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(11));

    t_0_xyz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(12));

    t_0_xyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(13));

    t_0_xyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(14));

    t_0_xyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(15));

    t_0_xyy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(16));

    t_0_xyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(17));

    t_0_xxz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(18));

    t_0_xxz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(19));

    t_0_xxz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(20));

    t_0_xxy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(21));

    t_0_xxy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(22));

    t_0_xxy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(23));

    t_0_xxy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(24));

    t_0_xxx_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(25));

    t_0_xxx_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(26));

    t_0_xxx_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(27));

    // set up (SFSG) integral components

    t_0_zzz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(0));

    t_0_zzz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(1));

    t_0_zzz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(2));

    t_0_yzz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(3));

    t_0_yzz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(4));

    t_0_yzz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(5));

    t_0_yyz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(6));

    t_0_yyz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(7));

    t_0_yyz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(8));

    t_0_yyy_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(9));

    t_0_yyy_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(10));

    t_0_yyy_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(11));

    t_0_xzz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(12));

    t_0_xzz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(13));

    t_0_xzz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(14));

    t_0_xyz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(15));

    t_0_xyz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(16));

    t_0_xyz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(17));

    t_0_xyy_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(18));

    t_0_xyy_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(19));

    t_0_xyy_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(20));

    t_0_xxz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(21));

    t_0_xxz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(22));

    t_0_xxz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(23));

    t_0_xxy_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(24));

    t_0_xxy_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(25));

    t_0_xxy_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(26));

    t_0_xxx_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(27));

    t_0_xxx_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(28));

    t_0_xxx_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(29));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxy_0_xxyz, t_0_xxy_0_xxz, t_0_xxy_y_xxz,\
                           t_0_xxz_0_xxxz, t_0_xxz_0_xxyz, t_0_xxz_0_xxz, t_0_xxz_0_xxzz,\
                           t_0_xxz_0_xyz, t_0_xxz_0_xzz, t_0_xxz_x_xxz, t_0_xxz_x_xyz,\
                           t_0_xxz_x_xzz, t_0_xxz_y_xxz, t_0_xxz_z_xxz, t_0_xyy_0_xxyy,\
                           t_0_xyy_0_xyy, t_0_xyy_0_xyyy, t_0_xyy_0_xyyz, t_0_xyy_0_xyz,\
                           t_0_xyy_0_yyy, t_0_xyy_0_yyz, t_0_xyy_x_xyy, t_0_xyy_x_yyy,\
                           t_0_xyy_x_yyz, t_0_xyy_y_xyy, t_0_xyy_y_xyz, t_0_xyz_0_xxyz,\
                           t_0_xyz_0_xyyz, t_0_xyz_0_xyz, t_0_xyz_0_xyzz, t_0_xyz_0_xzz,\
                           t_0_xyz_0_yyz, t_0_xyz_0_yzz, t_0_xyz_x_xyz, t_0_xyz_x_yyz,\
                           t_0_xyz_x_yzz, t_0_xyz_y_xyz, t_0_xyz_y_xzz, t_0_xyz_z_xyz,\
                           t_0_xzz_0_xxzz, t_0_xzz_0_xyzz, t_0_xzz_0_xzz, t_0_xzz_0_xzzz,\
                           t_0_xzz_0_yzz, t_0_xzz_0_zzz, t_0_xzz_x_xzz, t_0_xzz_x_yzz,\
                           t_0_xzz_x_zzz, t_0_xzz_y_xzz, t_0_xzz_z_xzz, t_0_yyy_0_xyyy,\
                           t_0_yyy_0_yyy, t_0_yyy_0_yyyy, t_0_yyy_0_yyyz, t_0_yyy_0_yyz,\
                           t_0_yyy_x_yyy, t_0_yyy_y_yyy, t_0_yyy_y_yyz, t_0_yyz_0_xyyz,\
                           t_0_yyz_0_yyyz, t_0_yyz_0_yyz, t_0_yyz_0_yyzz, t_0_yyz_0_yzz,\
                           t_0_yyz_x_yyz, t_0_yyz_y_yyz, t_0_yyz_y_yzz, t_0_yyz_z_yyz,\
                           t_0_yzz_0_xyzz, t_0_yzz_0_yyzz, t_0_yzz_0_yzz, t_0_yzz_0_yzzz,\
                           t_0_yzz_0_zzz, t_0_yzz_x_yzz, t_0_yzz_y_yzz, t_0_yzz_y_zzz,\
                           t_0_yzz_z_yzz, t_0_zzz_0_xzzz, t_0_zzz_0_yzzz, t_0_zzz_0_zzz,\
                           t_0_zzz_0_zzzz, t_0_zzz_x_zzz, t_0_zzz_y_zzz, t_0_zzz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_zzz[i] = t_0_zzz_0_zzzz[i] - rcd_z[i] * t_0_zzz_0_zzz[i];

        t_0_zzz_y_zzz[i] = t_0_zzz_0_yzzz[i] - rcd_y[i] * t_0_zzz_0_zzz[i];

        t_0_zzz_x_zzz[i] = t_0_zzz_0_xzzz[i] - rcd_x[i] * t_0_zzz_0_zzz[i];

        t_0_yzz_z_yzz[i] = t_0_yzz_0_yzzz[i] - rcd_z[i] * t_0_yzz_0_yzz[i];

        t_0_yzz_y_zzz[i] = t_0_yzz_0_yzzz[i] - rcd_y[i] * t_0_yzz_0_zzz[i];

        t_0_yzz_y_yzz[i] = t_0_yzz_0_yyzz[i] - rcd_y[i] * t_0_yzz_0_yzz[i];

        t_0_yzz_x_yzz[i] = t_0_yzz_0_xyzz[i] - rcd_x[i] * t_0_yzz_0_yzz[i];

        t_0_yyz_z_yyz[i] = t_0_yyz_0_yyzz[i] - rcd_z[i] * t_0_yyz_0_yyz[i];

        t_0_yyz_y_yzz[i] = t_0_yyz_0_yyzz[i] - rcd_y[i] * t_0_yyz_0_yzz[i];

        t_0_yyz_y_yyz[i] = t_0_yyz_0_yyyz[i] - rcd_y[i] * t_0_yyz_0_yyz[i];

        t_0_yyz_x_yyz[i] = t_0_yyz_0_xyyz[i] - rcd_x[i] * t_0_yyz_0_yyz[i];

        t_0_yyy_y_yyz[i] = t_0_yyy_0_yyyz[i] - rcd_y[i] * t_0_yyy_0_yyz[i];

        t_0_yyy_y_yyy[i] = t_0_yyy_0_yyyy[i] - rcd_y[i] * t_0_yyy_0_yyy[i];

        t_0_yyy_x_yyy[i] = t_0_yyy_0_xyyy[i] - rcd_x[i] * t_0_yyy_0_yyy[i];

        t_0_xzz_z_xzz[i] = t_0_xzz_0_xzzz[i] - rcd_z[i] * t_0_xzz_0_xzz[i];

        t_0_xzz_y_xzz[i] = t_0_xzz_0_xyzz[i] - rcd_y[i] * t_0_xzz_0_xzz[i];

        t_0_xzz_x_zzz[i] = t_0_xzz_0_xzzz[i] - rcd_x[i] * t_0_xzz_0_zzz[i];

        t_0_xzz_x_yzz[i] = t_0_xzz_0_xyzz[i] - rcd_x[i] * t_0_xzz_0_yzz[i];

        t_0_xzz_x_xzz[i] = t_0_xzz_0_xxzz[i] - rcd_x[i] * t_0_xzz_0_xzz[i];

        t_0_xyz_z_xyz[i] = t_0_xyz_0_xyzz[i] - rcd_z[i] * t_0_xyz_0_xyz[i];

        t_0_xyz_y_xzz[i] = t_0_xyz_0_xyzz[i] - rcd_y[i] * t_0_xyz_0_xzz[i];

        t_0_xyz_y_xyz[i] = t_0_xyz_0_xyyz[i] - rcd_y[i] * t_0_xyz_0_xyz[i];

        t_0_xyz_x_yzz[i] = t_0_xyz_0_xyzz[i] - rcd_x[i] * t_0_xyz_0_yzz[i];

        t_0_xyz_x_yyz[i] = t_0_xyz_0_xyyz[i] - rcd_x[i] * t_0_xyz_0_yyz[i];

        t_0_xyz_x_xyz[i] = t_0_xyz_0_xxyz[i] - rcd_x[i] * t_0_xyz_0_xyz[i];

        t_0_xyy_y_xyz[i] = t_0_xyy_0_xyyz[i] - rcd_y[i] * t_0_xyy_0_xyz[i];

        t_0_xyy_y_xyy[i] = t_0_xyy_0_xyyy[i] - rcd_y[i] * t_0_xyy_0_xyy[i];

        t_0_xyy_x_yyz[i] = t_0_xyy_0_xyyz[i] - rcd_x[i] * t_0_xyy_0_yyz[i];

        t_0_xyy_x_yyy[i] = t_0_xyy_0_xyyy[i] - rcd_x[i] * t_0_xyy_0_yyy[i];

        t_0_xyy_x_xyy[i] = t_0_xyy_0_xxyy[i] - rcd_x[i] * t_0_xyy_0_xyy[i];

        t_0_xxz_z_xxz[i] = t_0_xxz_0_xxzz[i] - rcd_z[i] * t_0_xxz_0_xxz[i];

        t_0_xxz_y_xxz[i] = t_0_xxz_0_xxyz[i] - rcd_y[i] * t_0_xxz_0_xxz[i];

        t_0_xxz_x_xzz[i] = t_0_xxz_0_xxzz[i] - rcd_x[i] * t_0_xxz_0_xzz[i];

        t_0_xxz_x_xyz[i] = t_0_xxz_0_xxyz[i] - rcd_x[i] * t_0_xxz_0_xyz[i];

        t_0_xxz_x_xxz[i] = t_0_xxz_0_xxxz[i] - rcd_x[i] * t_0_xxz_0_xxz[i];

        t_0_xxy_y_xxz[i] = t_0_xxy_0_xxyz[i] - rcd_y[i] * t_0_xxy_0_xxz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, t_0_xxx_0_xxx, t_0_xxx_0_xxxx, t_0_xxx_0_xxxy,\
                           t_0_xxx_0_xxxz, t_0_xxx_0_xxy, t_0_xxx_0_xxz, t_0_xxx_x_xxx,\
                           t_0_xxx_x_xxy, t_0_xxx_x_xxz, t_0_xxy_0_xxxy, t_0_xxy_0_xxy,\
                           t_0_xxy_0_xxyy, t_0_xxy_0_xxyz, t_0_xxy_0_xyy, t_0_xxy_0_xyz,\
                           t_0_xxy_x_xxy, t_0_xxy_x_xyy, t_0_xxy_x_xyz, t_0_xxy_y_xxy : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxy_y_xxy[i] = t_0_xxy_0_xxyy[i] - rcd_y[i] * t_0_xxy_0_xxy[i];

        t_0_xxy_x_xyz[i] = t_0_xxy_0_xxyz[i] - rcd_x[i] * t_0_xxy_0_xyz[i];

        t_0_xxy_x_xyy[i] = t_0_xxy_0_xxyy[i] - rcd_x[i] * t_0_xxy_0_xyy[i];

        t_0_xxy_x_xxy[i] = t_0_xxy_0_xxxy[i] - rcd_x[i] * t_0_xxy_0_xxy[i];

        t_0_xxx_x_xxz[i] = t_0_xxx_0_xxxz[i] - rcd_x[i] * t_0_xxx_0_xxz[i];

        t_0_xxx_x_xxy[i] = t_0_xxx_0_xxxy[i] - rcd_x[i] * t_0_xxx_0_xxy[i];

        t_0_xxx_x_xxx[i] = t_0_xxx_0_xxxx[i] - rcd_x[i] * t_0_xxx_0_xxx[i];
    }
}


} // derirec namespace
