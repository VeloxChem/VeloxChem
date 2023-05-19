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
compHostHRRForSFDD_V0(      BufferHostXY<T>&      intsBufferSFDD,
                      const BufferHostX<int32_t>& intsIndexesSFDD,
                      const BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostXY<T>&      intsBufferSFPF,
                      const BufferHostX<int32_t>& intsIndexesSFPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFDD) integral components

    t_0_zzz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(0));

    t_0_zzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(1));

    t_0_zzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(2));

    t_0_yzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(3));

    t_0_yzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(4));

    t_0_yzz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(5));

    t_0_yzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(6));

    t_0_yzz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(7));

    t_0_yzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(8));

    t_0_yyz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(9));

    t_0_yyz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(10));

    t_0_yyz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(11));

    t_0_yyz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(12));

    t_0_yyz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(13));

    t_0_yyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(14));

    t_0_yyy_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(15));

    t_0_yyy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(16));

    t_0_yyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(17));

    t_0_xzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(18));

    t_0_xzz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(19));

    t_0_xzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(20));

    t_0_xzz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(21));

    t_0_xzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(22));

    t_0_xzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(23));

    t_0_xyz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(24));

    t_0_xyz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(25));

    t_0_xyz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(26));

    t_0_xyz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(27));

    t_0_xyz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(28));

    t_0_xyz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(29));

    t_0_xyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(30));

    t_0_xyz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(31));

    t_0_xyz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(32));

    t_0_xyy_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(33));

    t_0_xyy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(34));

    t_0_xyy_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(35));

    t_0_xyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(36));

    t_0_xyy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(37));

    t_0_xyy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(38));

    t_0_xxz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(39));

    t_0_xxz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(40));

    t_0_xxz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(41));

    t_0_xxz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(42));

    t_0_xxz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(43));

    t_0_xxz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(44));

    t_0_xxy_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(45));

    t_0_xxy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(46));

    t_0_xxy_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(47));

    t_0_xxy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(48));

    t_0_xxy_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(49));

    t_0_xxy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(50));

    t_0_xxx_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(51));

    t_0_xxx_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(52));

    t_0_xxx_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(53));

    // set up (SFPD) integral components

    t_0_zzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(0));

    t_0_zzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(1));

    t_0_zzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(2));

    t_0_yzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(3));

    t_0_yzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(4));

    t_0_yzz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(5));

    t_0_yzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(6));

    t_0_yzz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(7));

    t_0_yyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(8));

    t_0_yyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(9));

    t_0_yyz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(10));

    t_0_yyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(11));

    t_0_yyz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(12));

    t_0_yyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(13));

    t_0_yyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(14));

    t_0_xzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(15));

    t_0_xzz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(16));

    t_0_xzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(17));

    t_0_xzz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(18));

    t_0_xyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(19));

    t_0_xyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(20));

    t_0_xyz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(21));

    t_0_xyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(22));

    t_0_xyz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(23));

    t_0_xyz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(24));

    t_0_xyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(25));

    t_0_xyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(26));

    t_0_xyy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(27));

    t_0_xxz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(28));

    t_0_xxz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(29));

    t_0_xxz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(30));

    t_0_xxz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(31));

    t_0_xxy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(32));

    t_0_xxy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(33));

    t_0_xxy_x_xx = intsBufferSFPD.data(intsIndexesSFPD(34));

    t_0_xxx_x_xx = intsBufferSFPD.data(intsIndexesSFPD(35));

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

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyy_x_yy, t_0_xyy_x_yyz, t_0_xyy_xz_yy,\
                           t_0_xyy_y_xy, t_0_xyy_y_xyy, t_0_xyy_y_xyz, t_0_xyy_yy_xy,\
                           t_0_xyy_yz_xy, t_0_xyz_x_xy, t_0_xyz_x_xyz, t_0_xyz_x_xz,\
                           t_0_xyz_x_yyz, t_0_xyz_x_yz, t_0_xyz_x_yzz, t_0_xyz_xx_yz,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yz, t_0_xyz_xz_xy, t_0_xyz_xz_yz,\
                           t_0_xyz_y_xy, t_0_xyz_y_xyz, t_0_xyz_y_xz, t_0_xyz_y_xzz,\
                           t_0_xyz_yy_xz, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_z_xy,\
                           t_0_xyz_z_xyz, t_0_xyz_zz_xy, t_0_xzz_x_xz, t_0_xzz_x_xzz,\
                           t_0_xzz_x_yzz, t_0_xzz_x_zz, t_0_xzz_x_zzz, t_0_xzz_xx_zz,\
                           t_0_xzz_xy_zz, t_0_xzz_xz_xz, t_0_xzz_xz_zz, t_0_xzz_y_xz,\
                           t_0_xzz_y_xzz, t_0_xzz_yz_xz, t_0_xzz_z_xz, t_0_xzz_z_xzz,\
                           t_0_xzz_zz_xz, t_0_yyy_x_yy, t_0_yyy_x_yyy, t_0_yyy_xy_yy,\
                           t_0_yyy_y_yy, t_0_yyy_y_yyy, t_0_yyy_y_yyz, t_0_yyy_yy_yy,\
                           t_0_yyy_yz_yy, t_0_yyz_x_yy, t_0_yyz_x_yyz, t_0_yyz_x_yz,\
                           t_0_yyz_xy_yz, t_0_yyz_xz_yy, t_0_yyz_y_yy, t_0_yyz_y_yyz,\
                           t_0_yyz_y_yz, t_0_yyz_y_yzz, t_0_yyz_yy_yz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_z_yy, t_0_yyz_z_yyz, t_0_yyz_zz_yy,\
                           t_0_yzz_x_yz, t_0_yzz_x_yzz, t_0_yzz_x_zz, t_0_yzz_xy_zz,\
                           t_0_yzz_xz_yz, t_0_yzz_y_yz, t_0_yzz_y_yzz, t_0_yzz_y_zz,\
                           t_0_yzz_y_zzz, t_0_yzz_yy_zz, t_0_yzz_yz_yz, t_0_yzz_yz_zz,\
                           t_0_yzz_z_yz, t_0_yzz_z_yzz, t_0_yzz_zz_yz, t_0_zzz_x_zz,\
                           t_0_zzz_x_zzz, t_0_zzz_xz_zz, t_0_zzz_y_zz, t_0_zzz_y_zzz,\
                           t_0_zzz_yz_zz, t_0_zzz_z_zz, t_0_zzz_z_zzz, t_0_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_zz_zz[i] = t_0_zzz_z_zzz[i] - rcd_z[i] * t_0_zzz_z_zz[i];

        t_0_zzz_yz_zz[i] = t_0_zzz_y_zzz[i] - rcd_z[i] * t_0_zzz_y_zz[i];

        t_0_zzz_xz_zz[i] = t_0_zzz_x_zzz[i] - rcd_z[i] * t_0_zzz_x_zz[i];

        t_0_yzz_zz_yz[i] = t_0_yzz_z_yzz[i] - rcd_z[i] * t_0_yzz_z_yz[i];

        t_0_yzz_yz_zz[i] = t_0_yzz_y_zzz[i] - rcd_z[i] * t_0_yzz_y_zz[i];

        t_0_yzz_yz_yz[i] = t_0_yzz_y_yzz[i] - rcd_z[i] * t_0_yzz_y_yz[i];

        t_0_yzz_yy_zz[i] = t_0_yzz_y_yzz[i] - rcd_y[i] * t_0_yzz_y_zz[i];

        t_0_yzz_xz_yz[i] = t_0_yzz_x_yzz[i] - rcd_z[i] * t_0_yzz_x_yz[i];

        t_0_yzz_xy_zz[i] = t_0_yzz_x_yzz[i] - rcd_y[i] * t_0_yzz_x_zz[i];

        t_0_yyz_zz_yy[i] = t_0_yyz_z_yyz[i] - rcd_z[i] * t_0_yyz_z_yy[i];

        t_0_yyz_yz_yz[i] = t_0_yyz_y_yzz[i] - rcd_z[i] * t_0_yyz_y_yz[i];

        t_0_yyz_yz_yy[i] = t_0_yyz_y_yyz[i] - rcd_z[i] * t_0_yyz_y_yy[i];

        t_0_yyz_yy_yz[i] = t_0_yyz_y_yyz[i] - rcd_y[i] * t_0_yyz_y_yz[i];

        t_0_yyz_xz_yy[i] = t_0_yyz_x_yyz[i] - rcd_z[i] * t_0_yyz_x_yy[i];

        t_0_yyz_xy_yz[i] = t_0_yyz_x_yyz[i] - rcd_y[i] * t_0_yyz_x_yz[i];

        t_0_yyy_yz_yy[i] = t_0_yyy_y_yyz[i] - rcd_z[i] * t_0_yyy_y_yy[i];

        t_0_yyy_yy_yy[i] = t_0_yyy_y_yyy[i] - rcd_y[i] * t_0_yyy_y_yy[i];

        t_0_yyy_xy_yy[i] = t_0_yyy_x_yyy[i] - rcd_y[i] * t_0_yyy_x_yy[i];

        t_0_xzz_zz_xz[i] = t_0_xzz_z_xzz[i] - rcd_z[i] * t_0_xzz_z_xz[i];

        t_0_xzz_yz_xz[i] = t_0_xzz_y_xzz[i] - rcd_z[i] * t_0_xzz_y_xz[i];

        t_0_xzz_xz_zz[i] = t_0_xzz_x_zzz[i] - rcd_z[i] * t_0_xzz_x_zz[i];

        t_0_xzz_xz_xz[i] = t_0_xzz_x_xzz[i] - rcd_z[i] * t_0_xzz_x_xz[i];

        t_0_xzz_xy_zz[i] = t_0_xzz_x_yzz[i] - rcd_y[i] * t_0_xzz_x_zz[i];

        t_0_xzz_xx_zz[i] = t_0_xzz_x_xzz[i] - rcd_x[i] * t_0_xzz_x_zz[i];

        t_0_xyz_zz_xy[i] = t_0_xyz_z_xyz[i] - rcd_z[i] * t_0_xyz_z_xy[i];

        t_0_xyz_yz_xz[i] = t_0_xyz_y_xzz[i] - rcd_z[i] * t_0_xyz_y_xz[i];

        t_0_xyz_yz_xy[i] = t_0_xyz_y_xyz[i] - rcd_z[i] * t_0_xyz_y_xy[i];

        t_0_xyz_yy_xz[i] = t_0_xyz_y_xyz[i] - rcd_y[i] * t_0_xyz_y_xz[i];

        t_0_xyz_xz_yz[i] = t_0_xyz_x_yzz[i] - rcd_z[i] * t_0_xyz_x_yz[i];

        t_0_xyz_xz_xy[i] = t_0_xyz_x_xyz[i] - rcd_z[i] * t_0_xyz_x_xy[i];

        t_0_xyz_xy_yz[i] = t_0_xyz_x_yyz[i] - rcd_y[i] * t_0_xyz_x_yz[i];

        t_0_xyz_xy_xz[i] = t_0_xyz_x_xyz[i] - rcd_y[i] * t_0_xyz_x_xz[i];

        t_0_xyz_xx_yz[i] = t_0_xyz_x_xyz[i] - rcd_x[i] * t_0_xyz_x_yz[i];

        t_0_xyy_yz_xy[i] = t_0_xyy_y_xyz[i] - rcd_z[i] * t_0_xyy_y_xy[i];

        t_0_xyy_yy_xy[i] = t_0_xyy_y_xyy[i] - rcd_y[i] * t_0_xyy_y_xy[i];

        t_0_xyy_xz_yy[i] = t_0_xyy_x_yyz[i] - rcd_z[i] * t_0_xyy_x_yy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_x_xx, t_0_xxx_x_xxx, t_0_xxx_x_xxy,\
                           t_0_xxx_x_xxz, t_0_xxx_xx_xx, t_0_xxx_xy_xx, t_0_xxx_xz_xx,\
                           t_0_xxy_x_xx, t_0_xxy_x_xxy, t_0_xxy_x_xy, t_0_xxy_x_xyy,\
                           t_0_xxy_x_xyz, t_0_xxy_xx_xy, t_0_xxy_xy_xx, t_0_xxy_xy_xy,\
                           t_0_xxy_xz_xy, t_0_xxy_y_xx, t_0_xxy_y_xxy, t_0_xxy_y_xxz,\
                           t_0_xxy_yy_xx, t_0_xxy_yz_xx, t_0_xxz_x_xx, t_0_xxz_x_xxz,\
                           t_0_xxz_x_xyz, t_0_xxz_x_xz, t_0_xxz_x_xzz, t_0_xxz_xx_xz,\
                           t_0_xxz_xy_xz, t_0_xxz_xz_xx, t_0_xxz_xz_xz, t_0_xxz_y_xx,\
                           t_0_xxz_y_xxz, t_0_xxz_yz_xx, t_0_xxz_z_xx, t_0_xxz_z_xxz,\
                           t_0_xxz_zz_xx, t_0_xyy_x_xy, t_0_xyy_x_xyy, t_0_xyy_x_yy,\
                           t_0_xyy_x_yyy, t_0_xyy_xx_yy, t_0_xyy_xy_xy, t_0_xyy_xy_yy : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyy_xy_yy[i] = t_0_xyy_x_yyy[i] - rcd_y[i] * t_0_xyy_x_yy[i];

        t_0_xyy_xy_xy[i] = t_0_xyy_x_xyy[i] - rcd_y[i] * t_0_xyy_x_xy[i];

        t_0_xyy_xx_yy[i] = t_0_xyy_x_xyy[i] - rcd_x[i] * t_0_xyy_x_yy[i];

        t_0_xxz_zz_xx[i] = t_0_xxz_z_xxz[i] - rcd_z[i] * t_0_xxz_z_xx[i];

        t_0_xxz_yz_xx[i] = t_0_xxz_y_xxz[i] - rcd_z[i] * t_0_xxz_y_xx[i];

        t_0_xxz_xz_xz[i] = t_0_xxz_x_xzz[i] - rcd_z[i] * t_0_xxz_x_xz[i];

        t_0_xxz_xz_xx[i] = t_0_xxz_x_xxz[i] - rcd_z[i] * t_0_xxz_x_xx[i];

        t_0_xxz_xy_xz[i] = t_0_xxz_x_xyz[i] - rcd_y[i] * t_0_xxz_x_xz[i];

        t_0_xxz_xx_xz[i] = t_0_xxz_x_xxz[i] - rcd_x[i] * t_0_xxz_x_xz[i];

        t_0_xxy_yz_xx[i] = t_0_xxy_y_xxz[i] - rcd_z[i] * t_0_xxy_y_xx[i];

        t_0_xxy_yy_xx[i] = t_0_xxy_y_xxy[i] - rcd_y[i] * t_0_xxy_y_xx[i];

        t_0_xxy_xz_xy[i] = t_0_xxy_x_xyz[i] - rcd_z[i] * t_0_xxy_x_xy[i];

        t_0_xxy_xy_xy[i] = t_0_xxy_x_xyy[i] - rcd_y[i] * t_0_xxy_x_xy[i];

        t_0_xxy_xy_xx[i] = t_0_xxy_x_xxy[i] - rcd_y[i] * t_0_xxy_x_xx[i];

        t_0_xxy_xx_xy[i] = t_0_xxy_x_xxy[i] - rcd_x[i] * t_0_xxy_x_xy[i];

        t_0_xxx_xz_xx[i] = t_0_xxx_x_xxz[i] - rcd_z[i] * t_0_xxx_x_xx[i];

        t_0_xxx_xy_xx[i] = t_0_xxx_x_xxy[i] - rcd_y[i] * t_0_xxx_x_xx[i];

        t_0_xxx_xx_xx[i] = t_0_xxx_x_xxx[i] - rcd_x[i] * t_0_xxx_x_xx[i];
    }
}


} // derirec namespace
