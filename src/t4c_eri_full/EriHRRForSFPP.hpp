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
compHostHRRForSFPP_V0(      BufferHostXY<T>&      intsBufferSFPP,
                      const BufferHostX<int32_t>& intsIndexesSFPP,
                      const BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPP) integral components

    t_0_zzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(0));

    t_0_zzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(1));

    t_0_zzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(2));

    t_0_zzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(3));

    t_0_zzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(4));

    t_0_zzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(5));

    t_0_zzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(6));

    t_0_zzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(7));

    t_0_zzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(8));

    t_0_yzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(9));

    t_0_yzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(10));

    t_0_yzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(11));

    t_0_yzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(12));

    t_0_yzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(13));

    t_0_yzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(14));

    t_0_yzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(15));

    t_0_yzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(16));

    t_0_yzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(17));

    t_0_yyz_z_z = intsBufferSFPP.data(intsIndexesSFPP(18));

    t_0_yyz_z_y = intsBufferSFPP.data(intsIndexesSFPP(19));

    t_0_yyz_z_x = intsBufferSFPP.data(intsIndexesSFPP(20));

    t_0_yyz_y_z = intsBufferSFPP.data(intsIndexesSFPP(21));

    t_0_yyz_y_y = intsBufferSFPP.data(intsIndexesSFPP(22));

    t_0_yyz_y_x = intsBufferSFPP.data(intsIndexesSFPP(23));

    t_0_yyz_x_z = intsBufferSFPP.data(intsIndexesSFPP(24));

    t_0_yyz_x_y = intsBufferSFPP.data(intsIndexesSFPP(25));

    t_0_yyz_x_x = intsBufferSFPP.data(intsIndexesSFPP(26));

    t_0_yyy_z_z = intsBufferSFPP.data(intsIndexesSFPP(27));

    t_0_yyy_z_y = intsBufferSFPP.data(intsIndexesSFPP(28));

    t_0_yyy_z_x = intsBufferSFPP.data(intsIndexesSFPP(29));

    t_0_yyy_y_z = intsBufferSFPP.data(intsIndexesSFPP(30));

    t_0_yyy_y_y = intsBufferSFPP.data(intsIndexesSFPP(31));

    t_0_yyy_y_x = intsBufferSFPP.data(intsIndexesSFPP(32));

    t_0_yyy_x_z = intsBufferSFPP.data(intsIndexesSFPP(33));

    t_0_yyy_x_y = intsBufferSFPP.data(intsIndexesSFPP(34));

    t_0_yyy_x_x = intsBufferSFPP.data(intsIndexesSFPP(35));

    t_0_xzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(36));

    t_0_xzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(37));

    t_0_xzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(38));

    t_0_xzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(39));

    t_0_xzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(40));

    t_0_xzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(41));

    t_0_xzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(42));

    t_0_xzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(43));

    t_0_xzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(44));

    t_0_xyz_z_z = intsBufferSFPP.data(intsIndexesSFPP(45));

    t_0_xyz_z_y = intsBufferSFPP.data(intsIndexesSFPP(46));

    t_0_xyz_z_x = intsBufferSFPP.data(intsIndexesSFPP(47));

    t_0_xyz_y_z = intsBufferSFPP.data(intsIndexesSFPP(48));

    t_0_xyz_y_y = intsBufferSFPP.data(intsIndexesSFPP(49));

    t_0_xyz_y_x = intsBufferSFPP.data(intsIndexesSFPP(50));

    t_0_xyz_x_z = intsBufferSFPP.data(intsIndexesSFPP(51));

    t_0_xyz_x_y = intsBufferSFPP.data(intsIndexesSFPP(52));

    t_0_xyz_x_x = intsBufferSFPP.data(intsIndexesSFPP(53));

    t_0_xyy_z_z = intsBufferSFPP.data(intsIndexesSFPP(54));

    t_0_xyy_z_y = intsBufferSFPP.data(intsIndexesSFPP(55));

    t_0_xyy_z_x = intsBufferSFPP.data(intsIndexesSFPP(56));

    t_0_xyy_y_z = intsBufferSFPP.data(intsIndexesSFPP(57));

    t_0_xyy_y_y = intsBufferSFPP.data(intsIndexesSFPP(58));

    t_0_xyy_y_x = intsBufferSFPP.data(intsIndexesSFPP(59));

    t_0_xyy_x_z = intsBufferSFPP.data(intsIndexesSFPP(60));

    t_0_xyy_x_y = intsBufferSFPP.data(intsIndexesSFPP(61));

    t_0_xyy_x_x = intsBufferSFPP.data(intsIndexesSFPP(62));

    t_0_xxz_z_z = intsBufferSFPP.data(intsIndexesSFPP(63));

    t_0_xxz_z_y = intsBufferSFPP.data(intsIndexesSFPP(64));

    t_0_xxz_z_x = intsBufferSFPP.data(intsIndexesSFPP(65));

    t_0_xxz_y_z = intsBufferSFPP.data(intsIndexesSFPP(66));

    t_0_xxz_y_y = intsBufferSFPP.data(intsIndexesSFPP(67));

    t_0_xxz_y_x = intsBufferSFPP.data(intsIndexesSFPP(68));

    t_0_xxz_x_z = intsBufferSFPP.data(intsIndexesSFPP(69));

    t_0_xxz_x_y = intsBufferSFPP.data(intsIndexesSFPP(70));

    t_0_xxz_x_x = intsBufferSFPP.data(intsIndexesSFPP(71));

    t_0_xxy_z_z = intsBufferSFPP.data(intsIndexesSFPP(72));

    t_0_xxy_z_y = intsBufferSFPP.data(intsIndexesSFPP(73));

    t_0_xxy_z_x = intsBufferSFPP.data(intsIndexesSFPP(74));

    t_0_xxy_y_z = intsBufferSFPP.data(intsIndexesSFPP(75));

    t_0_xxy_y_y = intsBufferSFPP.data(intsIndexesSFPP(76));

    t_0_xxy_y_x = intsBufferSFPP.data(intsIndexesSFPP(77));

    t_0_xxy_x_z = intsBufferSFPP.data(intsIndexesSFPP(78));

    t_0_xxy_x_y = intsBufferSFPP.data(intsIndexesSFPP(79));

    t_0_xxy_x_x = intsBufferSFPP.data(intsIndexesSFPP(80));

    t_0_xxx_z_z = intsBufferSFPP.data(intsIndexesSFPP(81));

    t_0_xxx_z_y = intsBufferSFPP.data(intsIndexesSFPP(82));

    t_0_xxx_z_x = intsBufferSFPP.data(intsIndexesSFPP(83));

    t_0_xxx_y_z = intsBufferSFPP.data(intsIndexesSFPP(84));

    t_0_xxx_y_y = intsBufferSFPP.data(intsIndexesSFPP(85));

    t_0_xxx_y_x = intsBufferSFPP.data(intsIndexesSFPP(86));

    t_0_xxx_x_z = intsBufferSFPP.data(intsIndexesSFPP(87));

    t_0_xxx_x_y = intsBufferSFPP.data(intsIndexesSFPP(88));

    t_0_xxx_x_x = intsBufferSFPP.data(intsIndexesSFPP(89));

    // set up (SFSP) integral components

    t_0_zzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(0));

    t_0_zzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(1));

    t_0_zzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(2));

    t_0_yzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(3));

    t_0_yzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(4));

    t_0_yzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(5));

    t_0_yyz_0_z = intsBufferSFSP.data(intsIndexesSFSP(6));

    t_0_yyz_0_y = intsBufferSFSP.data(intsIndexesSFSP(7));

    t_0_yyz_0_x = intsBufferSFSP.data(intsIndexesSFSP(8));

    t_0_yyy_0_z = intsBufferSFSP.data(intsIndexesSFSP(9));

    t_0_yyy_0_y = intsBufferSFSP.data(intsIndexesSFSP(10));

    t_0_yyy_0_x = intsBufferSFSP.data(intsIndexesSFSP(11));

    t_0_xzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(12));

    t_0_xzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(13));

    t_0_xzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(14));

    t_0_xyz_0_z = intsBufferSFSP.data(intsIndexesSFSP(15));

    t_0_xyz_0_y = intsBufferSFSP.data(intsIndexesSFSP(16));

    t_0_xyz_0_x = intsBufferSFSP.data(intsIndexesSFSP(17));

    t_0_xyy_0_z = intsBufferSFSP.data(intsIndexesSFSP(18));

    t_0_xyy_0_y = intsBufferSFSP.data(intsIndexesSFSP(19));

    t_0_xyy_0_x = intsBufferSFSP.data(intsIndexesSFSP(20));

    t_0_xxz_0_z = intsBufferSFSP.data(intsIndexesSFSP(21));

    t_0_xxz_0_y = intsBufferSFSP.data(intsIndexesSFSP(22));

    t_0_xxz_0_x = intsBufferSFSP.data(intsIndexesSFSP(23));

    t_0_xxy_0_z = intsBufferSFSP.data(intsIndexesSFSP(24));

    t_0_xxy_0_y = intsBufferSFSP.data(intsIndexesSFSP(25));

    t_0_xxy_0_x = intsBufferSFSP.data(intsIndexesSFSP(26));

    t_0_xxx_0_z = intsBufferSFSP.data(intsIndexesSFSP(27));

    t_0_xxx_0_y = intsBufferSFSP.data(intsIndexesSFSP(28));

    t_0_xxx_0_x = intsBufferSFSP.data(intsIndexesSFSP(29));

    // set up (SFSD) integral components

    t_0_zzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(0));

    t_0_zzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(1));

    t_0_zzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(2));

    t_0_zzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(3));

    t_0_zzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(4));

    t_0_zzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(5));

    t_0_yzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(6));

    t_0_yzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(7));

    t_0_yzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(8));

    t_0_yzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(9));

    t_0_yzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(10));

    t_0_yzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(11));

    t_0_yyz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(12));

    t_0_yyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(13));

    t_0_yyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(14));

    t_0_yyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(15));

    t_0_yyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(16));

    t_0_yyz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(17));

    t_0_yyy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(18));

    t_0_yyy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(19));

    t_0_yyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(20));

    t_0_yyy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(21));

    t_0_yyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(22));

    t_0_yyy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(23));

    t_0_xzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(24));

    t_0_xzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(25));

    t_0_xzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(26));

    t_0_xzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(27));

    t_0_xzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(28));

    t_0_xzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(29));

    t_0_xyz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(30));

    t_0_xyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(31));

    t_0_xyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(32));

    t_0_xyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(33));

    t_0_xyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(34));

    t_0_xyz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(35));

    t_0_xyy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(36));

    t_0_xyy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(37));

    t_0_xyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(38));

    t_0_xyy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(39));

    t_0_xyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(40));

    t_0_xyy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(41));

    t_0_xxz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(42));

    t_0_xxz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(43));

    t_0_xxz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(44));

    t_0_xxz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(45));

    t_0_xxz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(46));

    t_0_xxz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(47));

    t_0_xxy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(48));

    t_0_xxy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(49));

    t_0_xxy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(50));

    t_0_xxy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(51));

    t_0_xxy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(52));

    t_0_xxy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(53));

    t_0_xxx_0_zz = intsBufferSFSD.data(intsIndexesSFSD(54));

    t_0_xxx_0_yz = intsBufferSFSD.data(intsIndexesSFSD(55));

    t_0_xxx_0_yy = intsBufferSFSD.data(intsIndexesSFSD(56));

    t_0_xxx_0_xz = intsBufferSFSD.data(intsIndexesSFSD(57));

    t_0_xxx_0_xy = intsBufferSFSD.data(intsIndexesSFSD(58));

    t_0_xxx_0_xx = intsBufferSFSD.data(intsIndexesSFSD(59));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyy_0_x, t_0_yyy_0_xx, t_0_yyy_0_xy,\
                           t_0_yyy_0_xz, t_0_yyy_0_y, t_0_yyy_0_yy, t_0_yyy_0_yz, t_0_yyy_0_z,\
                           t_0_yyy_0_zz, t_0_yyy_x_x, t_0_yyy_x_y, t_0_yyy_x_z, t_0_yyy_y_x,\
                           t_0_yyy_y_y, t_0_yyy_y_z, t_0_yyy_z_x, t_0_yyy_z_y, t_0_yyy_z_z,\
                           t_0_yyz_0_x, t_0_yyz_0_xx, t_0_yyz_0_xy, t_0_yyz_0_xz, t_0_yyz_0_y,\
                           t_0_yyz_0_yy, t_0_yyz_0_yz, t_0_yyz_0_z, t_0_yyz_0_zz, t_0_yyz_x_x,\
                           t_0_yyz_x_y, t_0_yyz_x_z, t_0_yyz_y_x, t_0_yyz_y_y, t_0_yyz_y_z,\
                           t_0_yyz_z_x, t_0_yyz_z_y, t_0_yyz_z_z, t_0_yzz_0_x, t_0_yzz_0_xx,\
                           t_0_yzz_0_xy, t_0_yzz_0_xz, t_0_yzz_0_y, t_0_yzz_0_yy, t_0_yzz_0_yz,\
                           t_0_yzz_0_z, t_0_yzz_0_zz, t_0_yzz_x_x, t_0_yzz_x_y, t_0_yzz_x_z,\
                           t_0_yzz_y_x, t_0_yzz_y_y, t_0_yzz_y_z, t_0_yzz_z_x, t_0_yzz_z_y,\
                           t_0_yzz_z_z, t_0_zzz_0_x, t_0_zzz_0_xx, t_0_zzz_0_xy, t_0_zzz_0_xz,\
                           t_0_zzz_0_y, t_0_zzz_0_yy, t_0_zzz_0_yz, t_0_zzz_0_z, t_0_zzz_0_zz,\
                           t_0_zzz_x_x, t_0_zzz_x_y, t_0_zzz_x_z, t_0_zzz_y_x, t_0_zzz_y_y,\
                           t_0_zzz_y_z, t_0_zzz_z_x, t_0_zzz_z_y, t_0_zzz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_z[i] = t_0_zzz_0_zz[i] - rcd_z[i] * t_0_zzz_0_z[i];

        t_0_zzz_z_y[i] = t_0_zzz_0_yz[i] - rcd_z[i] * t_0_zzz_0_y[i];

        t_0_zzz_z_x[i] = t_0_zzz_0_xz[i] - rcd_z[i] * t_0_zzz_0_x[i];

        t_0_zzz_y_z[i] = t_0_zzz_0_yz[i] - rcd_y[i] * t_0_zzz_0_z[i];

        t_0_zzz_y_y[i] = t_0_zzz_0_yy[i] - rcd_y[i] * t_0_zzz_0_y[i];

        t_0_zzz_y_x[i] = t_0_zzz_0_xy[i] - rcd_y[i] * t_0_zzz_0_x[i];

        t_0_zzz_x_z[i] = t_0_zzz_0_xz[i] - rcd_x[i] * t_0_zzz_0_z[i];

        t_0_zzz_x_y[i] = t_0_zzz_0_xy[i] - rcd_x[i] * t_0_zzz_0_y[i];

        t_0_zzz_x_x[i] = t_0_zzz_0_xx[i] - rcd_x[i] * t_0_zzz_0_x[i];

        t_0_yzz_z_z[i] = t_0_yzz_0_zz[i] - rcd_z[i] * t_0_yzz_0_z[i];

        t_0_yzz_z_y[i] = t_0_yzz_0_yz[i] - rcd_z[i] * t_0_yzz_0_y[i];

        t_0_yzz_z_x[i] = t_0_yzz_0_xz[i] - rcd_z[i] * t_0_yzz_0_x[i];

        t_0_yzz_y_z[i] = t_0_yzz_0_yz[i] - rcd_y[i] * t_0_yzz_0_z[i];

        t_0_yzz_y_y[i] = t_0_yzz_0_yy[i] - rcd_y[i] * t_0_yzz_0_y[i];

        t_0_yzz_y_x[i] = t_0_yzz_0_xy[i] - rcd_y[i] * t_0_yzz_0_x[i];

        t_0_yzz_x_z[i] = t_0_yzz_0_xz[i] - rcd_x[i] * t_0_yzz_0_z[i];

        t_0_yzz_x_y[i] = t_0_yzz_0_xy[i] - rcd_x[i] * t_0_yzz_0_y[i];

        t_0_yzz_x_x[i] = t_0_yzz_0_xx[i] - rcd_x[i] * t_0_yzz_0_x[i];

        t_0_yyz_z_z[i] = t_0_yyz_0_zz[i] - rcd_z[i] * t_0_yyz_0_z[i];

        t_0_yyz_z_y[i] = t_0_yyz_0_yz[i] - rcd_z[i] * t_0_yyz_0_y[i];

        t_0_yyz_z_x[i] = t_0_yyz_0_xz[i] - rcd_z[i] * t_0_yyz_0_x[i];

        t_0_yyz_y_z[i] = t_0_yyz_0_yz[i] - rcd_y[i] * t_0_yyz_0_z[i];

        t_0_yyz_y_y[i] = t_0_yyz_0_yy[i] - rcd_y[i] * t_0_yyz_0_y[i];

        t_0_yyz_y_x[i] = t_0_yyz_0_xy[i] - rcd_y[i] * t_0_yyz_0_x[i];

        t_0_yyz_x_z[i] = t_0_yyz_0_xz[i] - rcd_x[i] * t_0_yyz_0_z[i];

        t_0_yyz_x_y[i] = t_0_yyz_0_xy[i] - rcd_x[i] * t_0_yyz_0_y[i];

        t_0_yyz_x_x[i] = t_0_yyz_0_xx[i] - rcd_x[i] * t_0_yyz_0_x[i];

        t_0_yyy_z_z[i] = t_0_yyy_0_zz[i] - rcd_z[i] * t_0_yyy_0_z[i];

        t_0_yyy_z_y[i] = t_0_yyy_0_yz[i] - rcd_z[i] * t_0_yyy_0_y[i];

        t_0_yyy_z_x[i] = t_0_yyy_0_xz[i] - rcd_z[i] * t_0_yyy_0_x[i];

        t_0_yyy_y_z[i] = t_0_yyy_0_yz[i] - rcd_y[i] * t_0_yyy_0_z[i];

        t_0_yyy_y_y[i] = t_0_yyy_0_yy[i] - rcd_y[i] * t_0_yyy_0_y[i];

        t_0_yyy_y_x[i] = t_0_yyy_0_xy[i] - rcd_y[i] * t_0_yyy_0_x[i];

        t_0_yyy_x_z[i] = t_0_yyy_0_xz[i] - rcd_x[i] * t_0_yyy_0_z[i];

        t_0_yyy_x_y[i] = t_0_yyy_0_xy[i] - rcd_x[i] * t_0_yyy_0_y[i];

        t_0_yyy_x_x[i] = t_0_yyy_0_xx[i] - rcd_x[i] * t_0_yyy_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxz_0_x, t_0_xxz_0_xx, t_0_xxz_0_xy,\
                           t_0_xxz_0_xz, t_0_xxz_0_y, t_0_xxz_0_yy, t_0_xxz_0_yz, t_0_xxz_0_z,\
                           t_0_xxz_0_zz, t_0_xxz_x_x, t_0_xxz_x_y, t_0_xxz_x_z, t_0_xxz_y_x,\
                           t_0_xxz_y_y, t_0_xxz_y_z, t_0_xxz_z_x, t_0_xxz_z_y, t_0_xxz_z_z,\
                           t_0_xyy_0_x, t_0_xyy_0_xx, t_0_xyy_0_xy, t_0_xyy_0_xz, t_0_xyy_0_y,\
                           t_0_xyy_0_yy, t_0_xyy_0_yz, t_0_xyy_0_z, t_0_xyy_0_zz, t_0_xyy_x_x,\
                           t_0_xyy_x_y, t_0_xyy_x_z, t_0_xyy_y_x, t_0_xyy_y_y, t_0_xyy_y_z,\
                           t_0_xyy_z_x, t_0_xyy_z_y, t_0_xyy_z_z, t_0_xyz_0_x, t_0_xyz_0_xx,\
                           t_0_xyz_0_xy, t_0_xyz_0_xz, t_0_xyz_0_y, t_0_xyz_0_yy, t_0_xyz_0_yz,\
                           t_0_xyz_0_z, t_0_xyz_0_zz, t_0_xyz_x_x, t_0_xyz_x_y, t_0_xyz_x_z,\
                           t_0_xyz_y_x, t_0_xyz_y_y, t_0_xyz_y_z, t_0_xyz_z_x, t_0_xyz_z_y,\
                           t_0_xyz_z_z, t_0_xzz_0_x, t_0_xzz_0_xx, t_0_xzz_0_xy, t_0_xzz_0_xz,\
                           t_0_xzz_0_y, t_0_xzz_0_yy, t_0_xzz_0_yz, t_0_xzz_0_z, t_0_xzz_0_zz,\
                           t_0_xzz_x_x, t_0_xzz_x_y, t_0_xzz_x_z, t_0_xzz_y_x, t_0_xzz_y_y,\
                           t_0_xzz_y_z, t_0_xzz_z_x, t_0_xzz_z_y, t_0_xzz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzz_z_z[i] = t_0_xzz_0_zz[i] - rcd_z[i] * t_0_xzz_0_z[i];

        t_0_xzz_z_y[i] = t_0_xzz_0_yz[i] - rcd_z[i] * t_0_xzz_0_y[i];

        t_0_xzz_z_x[i] = t_0_xzz_0_xz[i] - rcd_z[i] * t_0_xzz_0_x[i];

        t_0_xzz_y_z[i] = t_0_xzz_0_yz[i] - rcd_y[i] * t_0_xzz_0_z[i];

        t_0_xzz_y_y[i] = t_0_xzz_0_yy[i] - rcd_y[i] * t_0_xzz_0_y[i];

        t_0_xzz_y_x[i] = t_0_xzz_0_xy[i] - rcd_y[i] * t_0_xzz_0_x[i];

        t_0_xzz_x_z[i] = t_0_xzz_0_xz[i] - rcd_x[i] * t_0_xzz_0_z[i];

        t_0_xzz_x_y[i] = t_0_xzz_0_xy[i] - rcd_x[i] * t_0_xzz_0_y[i];

        t_0_xzz_x_x[i] = t_0_xzz_0_xx[i] - rcd_x[i] * t_0_xzz_0_x[i];

        t_0_xyz_z_z[i] = t_0_xyz_0_zz[i] - rcd_z[i] * t_0_xyz_0_z[i];

        t_0_xyz_z_y[i] = t_0_xyz_0_yz[i] - rcd_z[i] * t_0_xyz_0_y[i];

        t_0_xyz_z_x[i] = t_0_xyz_0_xz[i] - rcd_z[i] * t_0_xyz_0_x[i];

        t_0_xyz_y_z[i] = t_0_xyz_0_yz[i] - rcd_y[i] * t_0_xyz_0_z[i];

        t_0_xyz_y_y[i] = t_0_xyz_0_yy[i] - rcd_y[i] * t_0_xyz_0_y[i];

        t_0_xyz_y_x[i] = t_0_xyz_0_xy[i] - rcd_y[i] * t_0_xyz_0_x[i];

        t_0_xyz_x_z[i] = t_0_xyz_0_xz[i] - rcd_x[i] * t_0_xyz_0_z[i];

        t_0_xyz_x_y[i] = t_0_xyz_0_xy[i] - rcd_x[i] * t_0_xyz_0_y[i];

        t_0_xyz_x_x[i] = t_0_xyz_0_xx[i] - rcd_x[i] * t_0_xyz_0_x[i];

        t_0_xyy_z_z[i] = t_0_xyy_0_zz[i] - rcd_z[i] * t_0_xyy_0_z[i];

        t_0_xyy_z_y[i] = t_0_xyy_0_yz[i] - rcd_z[i] * t_0_xyy_0_y[i];

        t_0_xyy_z_x[i] = t_0_xyy_0_xz[i] - rcd_z[i] * t_0_xyy_0_x[i];

        t_0_xyy_y_z[i] = t_0_xyy_0_yz[i] - rcd_y[i] * t_0_xyy_0_z[i];

        t_0_xyy_y_y[i] = t_0_xyy_0_yy[i] - rcd_y[i] * t_0_xyy_0_y[i];

        t_0_xyy_y_x[i] = t_0_xyy_0_xy[i] - rcd_y[i] * t_0_xyy_0_x[i];

        t_0_xyy_x_z[i] = t_0_xyy_0_xz[i] - rcd_x[i] * t_0_xyy_0_z[i];

        t_0_xyy_x_y[i] = t_0_xyy_0_xy[i] - rcd_x[i] * t_0_xyy_0_y[i];

        t_0_xyy_x_x[i] = t_0_xyy_0_xx[i] - rcd_x[i] * t_0_xyy_0_x[i];

        t_0_xxz_z_z[i] = t_0_xxz_0_zz[i] - rcd_z[i] * t_0_xxz_0_z[i];

        t_0_xxz_z_y[i] = t_0_xxz_0_yz[i] - rcd_z[i] * t_0_xxz_0_y[i];

        t_0_xxz_z_x[i] = t_0_xxz_0_xz[i] - rcd_z[i] * t_0_xxz_0_x[i];

        t_0_xxz_y_z[i] = t_0_xxz_0_yz[i] - rcd_y[i] * t_0_xxz_0_z[i];

        t_0_xxz_y_y[i] = t_0_xxz_0_yy[i] - rcd_y[i] * t_0_xxz_0_y[i];

        t_0_xxz_y_x[i] = t_0_xxz_0_xy[i] - rcd_y[i] * t_0_xxz_0_x[i];

        t_0_xxz_x_z[i] = t_0_xxz_0_xz[i] - rcd_x[i] * t_0_xxz_0_z[i];

        t_0_xxz_x_y[i] = t_0_xxz_0_xy[i] - rcd_x[i] * t_0_xxz_0_y[i];

        t_0_xxz_x_x[i] = t_0_xxz_0_xx[i] - rcd_x[i] * t_0_xxz_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_0_x, t_0_xxx_0_xx, t_0_xxx_0_xy,\
                           t_0_xxx_0_xz, t_0_xxx_0_y, t_0_xxx_0_yy, t_0_xxx_0_yz, t_0_xxx_0_z,\
                           t_0_xxx_0_zz, t_0_xxx_x_x, t_0_xxx_x_y, t_0_xxx_x_z, t_0_xxx_y_x,\
                           t_0_xxx_y_y, t_0_xxx_y_z, t_0_xxx_z_x, t_0_xxx_z_y, t_0_xxx_z_z,\
                           t_0_xxy_0_x, t_0_xxy_0_xx, t_0_xxy_0_xy, t_0_xxy_0_xz, t_0_xxy_0_y,\
                           t_0_xxy_0_yy, t_0_xxy_0_yz, t_0_xxy_0_z, t_0_xxy_0_zz, t_0_xxy_x_x,\
                           t_0_xxy_x_y, t_0_xxy_x_z, t_0_xxy_y_x, t_0_xxy_y_y, t_0_xxy_y_z,\
                           t_0_xxy_z_x, t_0_xxy_z_y, t_0_xxy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxy_z_z[i] = t_0_xxy_0_zz[i] - rcd_z[i] * t_0_xxy_0_z[i];

        t_0_xxy_z_y[i] = t_0_xxy_0_yz[i] - rcd_z[i] * t_0_xxy_0_y[i];

        t_0_xxy_z_x[i] = t_0_xxy_0_xz[i] - rcd_z[i] * t_0_xxy_0_x[i];

        t_0_xxy_y_z[i] = t_0_xxy_0_yz[i] - rcd_y[i] * t_0_xxy_0_z[i];

        t_0_xxy_y_y[i] = t_0_xxy_0_yy[i] - rcd_y[i] * t_0_xxy_0_y[i];

        t_0_xxy_y_x[i] = t_0_xxy_0_xy[i] - rcd_y[i] * t_0_xxy_0_x[i];

        t_0_xxy_x_z[i] = t_0_xxy_0_xz[i] - rcd_x[i] * t_0_xxy_0_z[i];

        t_0_xxy_x_y[i] = t_0_xxy_0_xy[i] - rcd_x[i] * t_0_xxy_0_y[i];

        t_0_xxy_x_x[i] = t_0_xxy_0_xx[i] - rcd_x[i] * t_0_xxy_0_x[i];

        t_0_xxx_z_z[i] = t_0_xxx_0_zz[i] - rcd_z[i] * t_0_xxx_0_z[i];

        t_0_xxx_z_y[i] = t_0_xxx_0_yz[i] - rcd_z[i] * t_0_xxx_0_y[i];

        t_0_xxx_z_x[i] = t_0_xxx_0_xz[i] - rcd_z[i] * t_0_xxx_0_x[i];

        t_0_xxx_y_z[i] = t_0_xxx_0_yz[i] - rcd_y[i] * t_0_xxx_0_z[i];

        t_0_xxx_y_y[i] = t_0_xxx_0_yy[i] - rcd_y[i] * t_0_xxx_0_y[i];

        t_0_xxx_y_x[i] = t_0_xxx_0_xy[i] - rcd_y[i] * t_0_xxx_0_x[i];

        t_0_xxx_x_z[i] = t_0_xxx_0_xz[i] - rcd_x[i] * t_0_xxx_0_z[i];

        t_0_xxx_x_y[i] = t_0_xxx_0_xy[i] - rcd_x[i] * t_0_xxx_0_y[i];

        t_0_xxx_x_x[i] = t_0_xxx_0_xx[i] - rcd_x[i] * t_0_xxx_0_x[i];
    }
}


} // derirec namespace
