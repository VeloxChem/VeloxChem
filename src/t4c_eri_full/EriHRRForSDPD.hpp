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

    t_0_zz_z_yz = intsBufferSDPD.data(intsIndexesSDPD(1));

    t_0_zz_z_yy = intsBufferSDPD.data(intsIndexesSDPD(2));

    t_0_zz_z_xz = intsBufferSDPD.data(intsIndexesSDPD(3));

    t_0_zz_z_xy = intsBufferSDPD.data(intsIndexesSDPD(4));

    t_0_zz_z_xx = intsBufferSDPD.data(intsIndexesSDPD(5));

    t_0_zz_y_zz = intsBufferSDPD.data(intsIndexesSDPD(6));

    t_0_zz_y_yz = intsBufferSDPD.data(intsIndexesSDPD(7));

    t_0_zz_y_yy = intsBufferSDPD.data(intsIndexesSDPD(8));

    t_0_zz_y_xz = intsBufferSDPD.data(intsIndexesSDPD(9));

    t_0_zz_y_xy = intsBufferSDPD.data(intsIndexesSDPD(10));

    t_0_zz_y_xx = intsBufferSDPD.data(intsIndexesSDPD(11));

    t_0_zz_x_zz = intsBufferSDPD.data(intsIndexesSDPD(12));

    t_0_zz_x_yz = intsBufferSDPD.data(intsIndexesSDPD(13));

    t_0_zz_x_yy = intsBufferSDPD.data(intsIndexesSDPD(14));

    t_0_zz_x_xz = intsBufferSDPD.data(intsIndexesSDPD(15));

    t_0_zz_x_xy = intsBufferSDPD.data(intsIndexesSDPD(16));

    t_0_zz_x_xx = intsBufferSDPD.data(intsIndexesSDPD(17));

    t_0_yz_z_zz = intsBufferSDPD.data(intsIndexesSDPD(18));

    t_0_yz_z_yz = intsBufferSDPD.data(intsIndexesSDPD(19));

    t_0_yz_z_yy = intsBufferSDPD.data(intsIndexesSDPD(20));

    t_0_yz_z_xz = intsBufferSDPD.data(intsIndexesSDPD(21));

    t_0_yz_z_xy = intsBufferSDPD.data(intsIndexesSDPD(22));

    t_0_yz_z_xx = intsBufferSDPD.data(intsIndexesSDPD(23));

    t_0_yz_y_zz = intsBufferSDPD.data(intsIndexesSDPD(24));

    t_0_yz_y_yz = intsBufferSDPD.data(intsIndexesSDPD(25));

    t_0_yz_y_yy = intsBufferSDPD.data(intsIndexesSDPD(26));

    t_0_yz_y_xz = intsBufferSDPD.data(intsIndexesSDPD(27));

    t_0_yz_y_xy = intsBufferSDPD.data(intsIndexesSDPD(28));

    t_0_yz_y_xx = intsBufferSDPD.data(intsIndexesSDPD(29));

    t_0_yz_x_zz = intsBufferSDPD.data(intsIndexesSDPD(30));

    t_0_yz_x_yz = intsBufferSDPD.data(intsIndexesSDPD(31));

    t_0_yz_x_yy = intsBufferSDPD.data(intsIndexesSDPD(32));

    t_0_yz_x_xz = intsBufferSDPD.data(intsIndexesSDPD(33));

    t_0_yz_x_xy = intsBufferSDPD.data(intsIndexesSDPD(34));

    t_0_yz_x_xx = intsBufferSDPD.data(intsIndexesSDPD(35));

    t_0_yy_z_zz = intsBufferSDPD.data(intsIndexesSDPD(36));

    t_0_yy_z_yz = intsBufferSDPD.data(intsIndexesSDPD(37));

    t_0_yy_z_yy = intsBufferSDPD.data(intsIndexesSDPD(38));

    t_0_yy_z_xz = intsBufferSDPD.data(intsIndexesSDPD(39));

    t_0_yy_z_xy = intsBufferSDPD.data(intsIndexesSDPD(40));

    t_0_yy_z_xx = intsBufferSDPD.data(intsIndexesSDPD(41));

    t_0_yy_y_zz = intsBufferSDPD.data(intsIndexesSDPD(42));

    t_0_yy_y_yz = intsBufferSDPD.data(intsIndexesSDPD(43));

    t_0_yy_y_yy = intsBufferSDPD.data(intsIndexesSDPD(44));

    t_0_yy_y_xz = intsBufferSDPD.data(intsIndexesSDPD(45));

    t_0_yy_y_xy = intsBufferSDPD.data(intsIndexesSDPD(46));

    t_0_yy_y_xx = intsBufferSDPD.data(intsIndexesSDPD(47));

    t_0_yy_x_zz = intsBufferSDPD.data(intsIndexesSDPD(48));

    t_0_yy_x_yz = intsBufferSDPD.data(intsIndexesSDPD(49));

    t_0_yy_x_yy = intsBufferSDPD.data(intsIndexesSDPD(50));

    t_0_yy_x_xz = intsBufferSDPD.data(intsIndexesSDPD(51));

    t_0_yy_x_xy = intsBufferSDPD.data(intsIndexesSDPD(52));

    t_0_yy_x_xx = intsBufferSDPD.data(intsIndexesSDPD(53));

    t_0_xz_z_zz = intsBufferSDPD.data(intsIndexesSDPD(54));

    t_0_xz_z_yz = intsBufferSDPD.data(intsIndexesSDPD(55));

    t_0_xz_z_yy = intsBufferSDPD.data(intsIndexesSDPD(56));

    t_0_xz_z_xz = intsBufferSDPD.data(intsIndexesSDPD(57));

    t_0_xz_z_xy = intsBufferSDPD.data(intsIndexesSDPD(58));

    t_0_xz_z_xx = intsBufferSDPD.data(intsIndexesSDPD(59));

    t_0_xz_y_zz = intsBufferSDPD.data(intsIndexesSDPD(60));

    t_0_xz_y_yz = intsBufferSDPD.data(intsIndexesSDPD(61));

    t_0_xz_y_yy = intsBufferSDPD.data(intsIndexesSDPD(62));

    t_0_xz_y_xz = intsBufferSDPD.data(intsIndexesSDPD(63));

    t_0_xz_y_xy = intsBufferSDPD.data(intsIndexesSDPD(64));

    t_0_xz_y_xx = intsBufferSDPD.data(intsIndexesSDPD(65));

    t_0_xz_x_zz = intsBufferSDPD.data(intsIndexesSDPD(66));

    t_0_xz_x_yz = intsBufferSDPD.data(intsIndexesSDPD(67));

    t_0_xz_x_yy = intsBufferSDPD.data(intsIndexesSDPD(68));

    t_0_xz_x_xz = intsBufferSDPD.data(intsIndexesSDPD(69));

    t_0_xz_x_xy = intsBufferSDPD.data(intsIndexesSDPD(70));

    t_0_xz_x_xx = intsBufferSDPD.data(intsIndexesSDPD(71));

    t_0_xy_z_zz = intsBufferSDPD.data(intsIndexesSDPD(72));

    t_0_xy_z_yz = intsBufferSDPD.data(intsIndexesSDPD(73));

    t_0_xy_z_yy = intsBufferSDPD.data(intsIndexesSDPD(74));

    t_0_xy_z_xz = intsBufferSDPD.data(intsIndexesSDPD(75));

    t_0_xy_z_xy = intsBufferSDPD.data(intsIndexesSDPD(76));

    t_0_xy_z_xx = intsBufferSDPD.data(intsIndexesSDPD(77));

    t_0_xy_y_zz = intsBufferSDPD.data(intsIndexesSDPD(78));

    t_0_xy_y_yz = intsBufferSDPD.data(intsIndexesSDPD(79));

    t_0_xy_y_yy = intsBufferSDPD.data(intsIndexesSDPD(80));

    t_0_xy_y_xz = intsBufferSDPD.data(intsIndexesSDPD(81));

    t_0_xy_y_xy = intsBufferSDPD.data(intsIndexesSDPD(82));

    t_0_xy_y_xx = intsBufferSDPD.data(intsIndexesSDPD(83));

    t_0_xy_x_zz = intsBufferSDPD.data(intsIndexesSDPD(84));

    t_0_xy_x_yz = intsBufferSDPD.data(intsIndexesSDPD(85));

    t_0_xy_x_yy = intsBufferSDPD.data(intsIndexesSDPD(86));

    t_0_xy_x_xz = intsBufferSDPD.data(intsIndexesSDPD(87));

    t_0_xy_x_xy = intsBufferSDPD.data(intsIndexesSDPD(88));

    t_0_xy_x_xx = intsBufferSDPD.data(intsIndexesSDPD(89));

    t_0_xx_z_zz = intsBufferSDPD.data(intsIndexesSDPD(90));

    t_0_xx_z_yz = intsBufferSDPD.data(intsIndexesSDPD(91));

    t_0_xx_z_yy = intsBufferSDPD.data(intsIndexesSDPD(92));

    t_0_xx_z_xz = intsBufferSDPD.data(intsIndexesSDPD(93));

    t_0_xx_z_xy = intsBufferSDPD.data(intsIndexesSDPD(94));

    t_0_xx_z_xx = intsBufferSDPD.data(intsIndexesSDPD(95));

    t_0_xx_y_zz = intsBufferSDPD.data(intsIndexesSDPD(96));

    t_0_xx_y_yz = intsBufferSDPD.data(intsIndexesSDPD(97));

    t_0_xx_y_yy = intsBufferSDPD.data(intsIndexesSDPD(98));

    t_0_xx_y_xz = intsBufferSDPD.data(intsIndexesSDPD(99));

    t_0_xx_y_xy = intsBufferSDPD.data(intsIndexesSDPD(100));

    t_0_xx_y_xx = intsBufferSDPD.data(intsIndexesSDPD(101));

    t_0_xx_x_zz = intsBufferSDPD.data(intsIndexesSDPD(102));

    t_0_xx_x_yz = intsBufferSDPD.data(intsIndexesSDPD(103));

    t_0_xx_x_yy = intsBufferSDPD.data(intsIndexesSDPD(104));

    t_0_xx_x_xz = intsBufferSDPD.data(intsIndexesSDPD(105));

    t_0_xx_x_xy = intsBufferSDPD.data(intsIndexesSDPD(106));

    t_0_xx_x_xx = intsBufferSDPD.data(intsIndexesSDPD(107));

    // set up (SDSD) integral components

    t_0_zz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(0));

    t_0_zz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(1));

    t_0_zz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(2));

    t_0_zz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(3));

    t_0_zz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(4));

    t_0_zz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(5));

    t_0_yz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(6));

    t_0_yz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(7));

    t_0_yz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(8));

    t_0_yz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(9));

    t_0_yz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(10));

    t_0_yz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(11));

    t_0_yy_0_zz = intsBufferSDSD.data(intsIndexesSDSD(12));

    t_0_yy_0_yz = intsBufferSDSD.data(intsIndexesSDSD(13));

    t_0_yy_0_yy = intsBufferSDSD.data(intsIndexesSDSD(14));

    t_0_yy_0_xz = intsBufferSDSD.data(intsIndexesSDSD(15));

    t_0_yy_0_xy = intsBufferSDSD.data(intsIndexesSDSD(16));

    t_0_yy_0_xx = intsBufferSDSD.data(intsIndexesSDSD(17));

    t_0_xz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(18));

    t_0_xz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(19));

    t_0_xz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(20));

    t_0_xz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(21));

    t_0_xz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(22));

    t_0_xz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(23));

    t_0_xy_0_zz = intsBufferSDSD.data(intsIndexesSDSD(24));

    t_0_xy_0_yz = intsBufferSDSD.data(intsIndexesSDSD(25));

    t_0_xy_0_yy = intsBufferSDSD.data(intsIndexesSDSD(26));

    t_0_xy_0_xz = intsBufferSDSD.data(intsIndexesSDSD(27));

    t_0_xy_0_xy = intsBufferSDSD.data(intsIndexesSDSD(28));

    t_0_xy_0_xx = intsBufferSDSD.data(intsIndexesSDSD(29));

    t_0_xx_0_zz = intsBufferSDSD.data(intsIndexesSDSD(30));

    t_0_xx_0_yz = intsBufferSDSD.data(intsIndexesSDSD(31));

    t_0_xx_0_yy = intsBufferSDSD.data(intsIndexesSDSD(32));

    t_0_xx_0_xz = intsBufferSDSD.data(intsIndexesSDSD(33));

    t_0_xx_0_xy = intsBufferSDSD.data(intsIndexesSDSD(34));

    t_0_xx_0_xx = intsBufferSDSD.data(intsIndexesSDSD(35));

    // set up (SDSF) integral components

    t_0_zz_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(0));

    t_0_zz_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(1));

    t_0_zz_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(2));

    t_0_zz_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(3));

    t_0_zz_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(4));

    t_0_zz_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(5));

    t_0_zz_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(6));

    t_0_zz_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(7));

    t_0_zz_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(8));

    t_0_zz_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(9));

    t_0_yz_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(10));

    t_0_yz_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(11));

    t_0_yz_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(12));

    t_0_yz_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(13));

    t_0_yz_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(14));

    t_0_yz_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(15));

    t_0_yz_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(16));

    t_0_yz_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(17));

    t_0_yz_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(18));

    t_0_yz_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(19));

    t_0_yy_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(20));

    t_0_yy_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(21));

    t_0_yy_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(22));

    t_0_yy_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(23));

    t_0_yy_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(24));

    t_0_yy_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(25));

    t_0_yy_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(26));

    t_0_yy_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(27));

    t_0_yy_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(28));

    t_0_yy_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(29));

    t_0_xz_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(30));

    t_0_xz_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(31));

    t_0_xz_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(32));

    t_0_xz_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(33));

    t_0_xz_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(34));

    t_0_xz_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(35));

    t_0_xz_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(36));

    t_0_xz_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(37));

    t_0_xz_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(38));

    t_0_xz_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(39));

    t_0_xy_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(40));

    t_0_xy_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(41));

    t_0_xy_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(42));

    t_0_xy_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(43));

    t_0_xy_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(44));

    t_0_xy_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(45));

    t_0_xy_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(46));

    t_0_xy_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(47));

    t_0_xy_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(48));

    t_0_xy_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(49));

    t_0_xx_0_zzz = intsBufferSDSF.data(intsIndexesSDSF(50));

    t_0_xx_0_yzz = intsBufferSDSF.data(intsIndexesSDSF(51));

    t_0_xx_0_yyz = intsBufferSDSF.data(intsIndexesSDSF(52));

    t_0_xx_0_yyy = intsBufferSDSF.data(intsIndexesSDSF(53));

    t_0_xx_0_xzz = intsBufferSDSF.data(intsIndexesSDSF(54));

    t_0_xx_0_xyz = intsBufferSDSF.data(intsIndexesSDSF(55));

    t_0_xx_0_xyy = intsBufferSDSF.data(intsIndexesSDSF(56));

    t_0_xx_0_xxz = intsBufferSDSF.data(intsIndexesSDSF(57));

    t_0_xx_0_xxy = intsBufferSDSF.data(intsIndexesSDSF(58));

    t_0_xx_0_xxx = intsBufferSDSF.data(intsIndexesSDSF(59));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yz_0_xx, t_0_yz_0_xxx, t_0_yz_0_xxy,\
                           t_0_yz_0_xxz, t_0_yz_0_xy, t_0_yz_0_xyy, t_0_yz_0_xyz, t_0_yz_0_xz,\
                           t_0_yz_0_xzz, t_0_yz_0_yy, t_0_yz_0_yyy, t_0_yz_0_yyz, t_0_yz_0_yz,\
                           t_0_yz_0_yzz, t_0_yz_0_zz, t_0_yz_0_zzz, t_0_yz_x_xx, t_0_yz_x_xy,\
                           t_0_yz_x_xz, t_0_yz_x_yy, t_0_yz_x_yz, t_0_yz_x_zz, t_0_yz_y_xx,\
                           t_0_yz_y_xy, t_0_yz_y_xz, t_0_yz_y_yy, t_0_yz_y_yz, t_0_yz_y_zz,\
                           t_0_yz_z_xx, t_0_yz_z_xy, t_0_yz_z_xz, t_0_yz_z_yy, t_0_yz_z_yz,\
                           t_0_yz_z_zz, t_0_zz_0_xx, t_0_zz_0_xxx, t_0_zz_0_xxy, t_0_zz_0_xxz,\
                           t_0_zz_0_xy, t_0_zz_0_xyy, t_0_zz_0_xyz, t_0_zz_0_xz, t_0_zz_0_xzz,\
                           t_0_zz_0_yy, t_0_zz_0_yyy, t_0_zz_0_yyz, t_0_zz_0_yz, t_0_zz_0_yzz,\
                           t_0_zz_0_zz, t_0_zz_0_zzz, t_0_zz_x_xx, t_0_zz_x_xy, t_0_zz_x_xz,\
                           t_0_zz_x_yy, t_0_zz_x_yz, t_0_zz_x_zz, t_0_zz_y_xx, t_0_zz_y_xy,\
                           t_0_zz_y_xz, t_0_zz_y_yy, t_0_zz_y_yz, t_0_zz_y_zz, t_0_zz_z_xx,\
                           t_0_zz_z_xy, t_0_zz_z_xz, t_0_zz_z_yy, t_0_zz_z_yz, t_0_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_z_zz[i] = t_0_zz_0_zzz[i] - rcd_z[i] * t_0_zz_0_zz[i];

        t_0_zz_z_yz[i] = t_0_zz_0_yzz[i] - rcd_z[i] * t_0_zz_0_yz[i];

        t_0_zz_z_yy[i] = t_0_zz_0_yyz[i] - rcd_z[i] * t_0_zz_0_yy[i];

        t_0_zz_z_xz[i] = t_0_zz_0_xzz[i] - rcd_z[i] * t_0_zz_0_xz[i];

        t_0_zz_z_xy[i] = t_0_zz_0_xyz[i] - rcd_z[i] * t_0_zz_0_xy[i];

        t_0_zz_z_xx[i] = t_0_zz_0_xxz[i] - rcd_z[i] * t_0_zz_0_xx[i];

        t_0_zz_y_zz[i] = t_0_zz_0_yzz[i] - rcd_y[i] * t_0_zz_0_zz[i];

        t_0_zz_y_yz[i] = t_0_zz_0_yyz[i] - rcd_y[i] * t_0_zz_0_yz[i];

        t_0_zz_y_yy[i] = t_0_zz_0_yyy[i] - rcd_y[i] * t_0_zz_0_yy[i];

        t_0_zz_y_xz[i] = t_0_zz_0_xyz[i] - rcd_y[i] * t_0_zz_0_xz[i];

        t_0_zz_y_xy[i] = t_0_zz_0_xyy[i] - rcd_y[i] * t_0_zz_0_xy[i];

        t_0_zz_y_xx[i] = t_0_zz_0_xxy[i] - rcd_y[i] * t_0_zz_0_xx[i];

        t_0_zz_x_zz[i] = t_0_zz_0_xzz[i] - rcd_x[i] * t_0_zz_0_zz[i];

        t_0_zz_x_yz[i] = t_0_zz_0_xyz[i] - rcd_x[i] * t_0_zz_0_yz[i];

        t_0_zz_x_yy[i] = t_0_zz_0_xyy[i] - rcd_x[i] * t_0_zz_0_yy[i];

        t_0_zz_x_xz[i] = t_0_zz_0_xxz[i] - rcd_x[i] * t_0_zz_0_xz[i];

        t_0_zz_x_xy[i] = t_0_zz_0_xxy[i] - rcd_x[i] * t_0_zz_0_xy[i];

        t_0_zz_x_xx[i] = t_0_zz_0_xxx[i] - rcd_x[i] * t_0_zz_0_xx[i];

        t_0_yz_z_zz[i] = t_0_yz_0_zzz[i] - rcd_z[i] * t_0_yz_0_zz[i];

        t_0_yz_z_yz[i] = t_0_yz_0_yzz[i] - rcd_z[i] * t_0_yz_0_yz[i];

        t_0_yz_z_yy[i] = t_0_yz_0_yyz[i] - rcd_z[i] * t_0_yz_0_yy[i];

        t_0_yz_z_xz[i] = t_0_yz_0_xzz[i] - rcd_z[i] * t_0_yz_0_xz[i];

        t_0_yz_z_xy[i] = t_0_yz_0_xyz[i] - rcd_z[i] * t_0_yz_0_xy[i];

        t_0_yz_z_xx[i] = t_0_yz_0_xxz[i] - rcd_z[i] * t_0_yz_0_xx[i];

        t_0_yz_y_zz[i] = t_0_yz_0_yzz[i] - rcd_y[i] * t_0_yz_0_zz[i];

        t_0_yz_y_yz[i] = t_0_yz_0_yyz[i] - rcd_y[i] * t_0_yz_0_yz[i];

        t_0_yz_y_yy[i] = t_0_yz_0_yyy[i] - rcd_y[i] * t_0_yz_0_yy[i];

        t_0_yz_y_xz[i] = t_0_yz_0_xyz[i] - rcd_y[i] * t_0_yz_0_xz[i];

        t_0_yz_y_xy[i] = t_0_yz_0_xyy[i] - rcd_y[i] * t_0_yz_0_xy[i];

        t_0_yz_y_xx[i] = t_0_yz_0_xxy[i] - rcd_y[i] * t_0_yz_0_xx[i];

        t_0_yz_x_zz[i] = t_0_yz_0_xzz[i] - rcd_x[i] * t_0_yz_0_zz[i];

        t_0_yz_x_yz[i] = t_0_yz_0_xyz[i] - rcd_x[i] * t_0_yz_0_yz[i];

        t_0_yz_x_yy[i] = t_0_yz_0_xyy[i] - rcd_x[i] * t_0_yz_0_yy[i];

        t_0_yz_x_xz[i] = t_0_yz_0_xxz[i] - rcd_x[i] * t_0_yz_0_xz[i];

        t_0_yz_x_xy[i] = t_0_yz_0_xxy[i] - rcd_x[i] * t_0_yz_0_xy[i];

        t_0_yz_x_xx[i] = t_0_yz_0_xxx[i] - rcd_x[i] * t_0_yz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xz_0_xx, t_0_xz_0_xxx, t_0_xz_0_xxy,\
                           t_0_xz_0_xxz, t_0_xz_0_xy, t_0_xz_0_xyy, t_0_xz_0_xyz, t_0_xz_0_xz,\
                           t_0_xz_0_xzz, t_0_xz_0_yy, t_0_xz_0_yyy, t_0_xz_0_yyz, t_0_xz_0_yz,\
                           t_0_xz_0_yzz, t_0_xz_0_zz, t_0_xz_0_zzz, t_0_xz_x_xx, t_0_xz_x_xy,\
                           t_0_xz_x_xz, t_0_xz_x_yy, t_0_xz_x_yz, t_0_xz_x_zz, t_0_xz_y_xx,\
                           t_0_xz_y_xy, t_0_xz_y_xz, t_0_xz_y_yy, t_0_xz_y_yz, t_0_xz_y_zz,\
                           t_0_xz_z_xx, t_0_xz_z_xy, t_0_xz_z_xz, t_0_xz_z_yy, t_0_xz_z_yz,\
                           t_0_xz_z_zz, t_0_yy_0_xx, t_0_yy_0_xxx, t_0_yy_0_xxy, t_0_yy_0_xxz,\
                           t_0_yy_0_xy, t_0_yy_0_xyy, t_0_yy_0_xyz, t_0_yy_0_xz, t_0_yy_0_xzz,\
                           t_0_yy_0_yy, t_0_yy_0_yyy, t_0_yy_0_yyz, t_0_yy_0_yz, t_0_yy_0_yzz,\
                           t_0_yy_0_zz, t_0_yy_0_zzz, t_0_yy_x_xx, t_0_yy_x_xy, t_0_yy_x_xz,\
                           t_0_yy_x_yy, t_0_yy_x_yz, t_0_yy_x_zz, t_0_yy_y_xx, t_0_yy_y_xy,\
                           t_0_yy_y_xz, t_0_yy_y_yy, t_0_yy_y_yz, t_0_yy_y_zz, t_0_yy_z_xx,\
                           t_0_yy_z_xy, t_0_yy_z_xz, t_0_yy_z_yy, t_0_yy_z_yz, t_0_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yy_z_zz[i] = t_0_yy_0_zzz[i] - rcd_z[i] * t_0_yy_0_zz[i];

        t_0_yy_z_yz[i] = t_0_yy_0_yzz[i] - rcd_z[i] * t_0_yy_0_yz[i];

        t_0_yy_z_yy[i] = t_0_yy_0_yyz[i] - rcd_z[i] * t_0_yy_0_yy[i];

        t_0_yy_z_xz[i] = t_0_yy_0_xzz[i] - rcd_z[i] * t_0_yy_0_xz[i];

        t_0_yy_z_xy[i] = t_0_yy_0_xyz[i] - rcd_z[i] * t_0_yy_0_xy[i];

        t_0_yy_z_xx[i] = t_0_yy_0_xxz[i] - rcd_z[i] * t_0_yy_0_xx[i];

        t_0_yy_y_zz[i] = t_0_yy_0_yzz[i] - rcd_y[i] * t_0_yy_0_zz[i];

        t_0_yy_y_yz[i] = t_0_yy_0_yyz[i] - rcd_y[i] * t_0_yy_0_yz[i];

        t_0_yy_y_yy[i] = t_0_yy_0_yyy[i] - rcd_y[i] * t_0_yy_0_yy[i];

        t_0_yy_y_xz[i] = t_0_yy_0_xyz[i] - rcd_y[i] * t_0_yy_0_xz[i];

        t_0_yy_y_xy[i] = t_0_yy_0_xyy[i] - rcd_y[i] * t_0_yy_0_xy[i];

        t_0_yy_y_xx[i] = t_0_yy_0_xxy[i] - rcd_y[i] * t_0_yy_0_xx[i];

        t_0_yy_x_zz[i] = t_0_yy_0_xzz[i] - rcd_x[i] * t_0_yy_0_zz[i];

        t_0_yy_x_yz[i] = t_0_yy_0_xyz[i] - rcd_x[i] * t_0_yy_0_yz[i];

        t_0_yy_x_yy[i] = t_0_yy_0_xyy[i] - rcd_x[i] * t_0_yy_0_yy[i];

        t_0_yy_x_xz[i] = t_0_yy_0_xxz[i] - rcd_x[i] * t_0_yy_0_xz[i];

        t_0_yy_x_xy[i] = t_0_yy_0_xxy[i] - rcd_x[i] * t_0_yy_0_xy[i];

        t_0_yy_x_xx[i] = t_0_yy_0_xxx[i] - rcd_x[i] * t_0_yy_0_xx[i];

        t_0_xz_z_zz[i] = t_0_xz_0_zzz[i] - rcd_z[i] * t_0_xz_0_zz[i];

        t_0_xz_z_yz[i] = t_0_xz_0_yzz[i] - rcd_z[i] * t_0_xz_0_yz[i];

        t_0_xz_z_yy[i] = t_0_xz_0_yyz[i] - rcd_z[i] * t_0_xz_0_yy[i];

        t_0_xz_z_xz[i] = t_0_xz_0_xzz[i] - rcd_z[i] * t_0_xz_0_xz[i];

        t_0_xz_z_xy[i] = t_0_xz_0_xyz[i] - rcd_z[i] * t_0_xz_0_xy[i];

        t_0_xz_z_xx[i] = t_0_xz_0_xxz[i] - rcd_z[i] * t_0_xz_0_xx[i];

        t_0_xz_y_zz[i] = t_0_xz_0_yzz[i] - rcd_y[i] * t_0_xz_0_zz[i];

        t_0_xz_y_yz[i] = t_0_xz_0_yyz[i] - rcd_y[i] * t_0_xz_0_yz[i];

        t_0_xz_y_yy[i] = t_0_xz_0_yyy[i] - rcd_y[i] * t_0_xz_0_yy[i];

        t_0_xz_y_xz[i] = t_0_xz_0_xyz[i] - rcd_y[i] * t_0_xz_0_xz[i];

        t_0_xz_y_xy[i] = t_0_xz_0_xyy[i] - rcd_y[i] * t_0_xz_0_xy[i];

        t_0_xz_y_xx[i] = t_0_xz_0_xxy[i] - rcd_y[i] * t_0_xz_0_xx[i];

        t_0_xz_x_zz[i] = t_0_xz_0_xzz[i] - rcd_x[i] * t_0_xz_0_zz[i];

        t_0_xz_x_yz[i] = t_0_xz_0_xyz[i] - rcd_x[i] * t_0_xz_0_yz[i];

        t_0_xz_x_yy[i] = t_0_xz_0_xyy[i] - rcd_x[i] * t_0_xz_0_yy[i];

        t_0_xz_x_xz[i] = t_0_xz_0_xxz[i] - rcd_x[i] * t_0_xz_0_xz[i];

        t_0_xz_x_xy[i] = t_0_xz_0_xxy[i] - rcd_x[i] * t_0_xz_0_xy[i];

        t_0_xz_x_xx[i] = t_0_xz_0_xxx[i] - rcd_x[i] * t_0_xz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_0_xx, t_0_xx_0_xxx, t_0_xx_0_xxy,\
                           t_0_xx_0_xxz, t_0_xx_0_xy, t_0_xx_0_xyy, t_0_xx_0_xyz, t_0_xx_0_xz,\
                           t_0_xx_0_xzz, t_0_xx_0_yy, t_0_xx_0_yyy, t_0_xx_0_yyz, t_0_xx_0_yz,\
                           t_0_xx_0_yzz, t_0_xx_0_zz, t_0_xx_0_zzz, t_0_xx_x_xx, t_0_xx_x_xy,\
                           t_0_xx_x_xz, t_0_xx_x_yy, t_0_xx_x_yz, t_0_xx_x_zz, t_0_xx_y_xx,\
                           t_0_xx_y_xy, t_0_xx_y_xz, t_0_xx_y_yy, t_0_xx_y_yz, t_0_xx_y_zz,\
                           t_0_xx_z_xx, t_0_xx_z_xy, t_0_xx_z_xz, t_0_xx_z_yy, t_0_xx_z_yz,\
                           t_0_xx_z_zz, t_0_xy_0_xx, t_0_xy_0_xxx, t_0_xy_0_xxy, t_0_xy_0_xxz,\
                           t_0_xy_0_xy, t_0_xy_0_xyy, t_0_xy_0_xyz, t_0_xy_0_xz, t_0_xy_0_xzz,\
                           t_0_xy_0_yy, t_0_xy_0_yyy, t_0_xy_0_yyz, t_0_xy_0_yz, t_0_xy_0_yzz,\
                           t_0_xy_0_zz, t_0_xy_0_zzz, t_0_xy_x_xx, t_0_xy_x_xy, t_0_xy_x_xz,\
                           t_0_xy_x_yy, t_0_xy_x_yz, t_0_xy_x_zz, t_0_xy_y_xx, t_0_xy_y_xy,\
                           t_0_xy_y_xz, t_0_xy_y_yy, t_0_xy_y_yz, t_0_xy_y_zz, t_0_xy_z_xx,\
                           t_0_xy_z_xy, t_0_xy_z_xz, t_0_xy_z_yy, t_0_xy_z_yz, t_0_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xy_z_zz[i] = t_0_xy_0_zzz[i] - rcd_z[i] * t_0_xy_0_zz[i];

        t_0_xy_z_yz[i] = t_0_xy_0_yzz[i] - rcd_z[i] * t_0_xy_0_yz[i];

        t_0_xy_z_yy[i] = t_0_xy_0_yyz[i] - rcd_z[i] * t_0_xy_0_yy[i];

        t_0_xy_z_xz[i] = t_0_xy_0_xzz[i] - rcd_z[i] * t_0_xy_0_xz[i];

        t_0_xy_z_xy[i] = t_0_xy_0_xyz[i] - rcd_z[i] * t_0_xy_0_xy[i];

        t_0_xy_z_xx[i] = t_0_xy_0_xxz[i] - rcd_z[i] * t_0_xy_0_xx[i];

        t_0_xy_y_zz[i] = t_0_xy_0_yzz[i] - rcd_y[i] * t_0_xy_0_zz[i];

        t_0_xy_y_yz[i] = t_0_xy_0_yyz[i] - rcd_y[i] * t_0_xy_0_yz[i];

        t_0_xy_y_yy[i] = t_0_xy_0_yyy[i] - rcd_y[i] * t_0_xy_0_yy[i];

        t_0_xy_y_xz[i] = t_0_xy_0_xyz[i] - rcd_y[i] * t_0_xy_0_xz[i];

        t_0_xy_y_xy[i] = t_0_xy_0_xyy[i] - rcd_y[i] * t_0_xy_0_xy[i];

        t_0_xy_y_xx[i] = t_0_xy_0_xxy[i] - rcd_y[i] * t_0_xy_0_xx[i];

        t_0_xy_x_zz[i] = t_0_xy_0_xzz[i] - rcd_x[i] * t_0_xy_0_zz[i];

        t_0_xy_x_yz[i] = t_0_xy_0_xyz[i] - rcd_x[i] * t_0_xy_0_yz[i];

        t_0_xy_x_yy[i] = t_0_xy_0_xyy[i] - rcd_x[i] * t_0_xy_0_yy[i];

        t_0_xy_x_xz[i] = t_0_xy_0_xxz[i] - rcd_x[i] * t_0_xy_0_xz[i];

        t_0_xy_x_xy[i] = t_0_xy_0_xxy[i] - rcd_x[i] * t_0_xy_0_xy[i];

        t_0_xy_x_xx[i] = t_0_xy_0_xxx[i] - rcd_x[i] * t_0_xy_0_xx[i];

        t_0_xx_z_zz[i] = t_0_xx_0_zzz[i] - rcd_z[i] * t_0_xx_0_zz[i];

        t_0_xx_z_yz[i] = t_0_xx_0_yzz[i] - rcd_z[i] * t_0_xx_0_yz[i];

        t_0_xx_z_yy[i] = t_0_xx_0_yyz[i] - rcd_z[i] * t_0_xx_0_yy[i];

        t_0_xx_z_xz[i] = t_0_xx_0_xzz[i] - rcd_z[i] * t_0_xx_0_xz[i];

        t_0_xx_z_xy[i] = t_0_xx_0_xyz[i] - rcd_z[i] * t_0_xx_0_xy[i];

        t_0_xx_z_xx[i] = t_0_xx_0_xxz[i] - rcd_z[i] * t_0_xx_0_xx[i];

        t_0_xx_y_zz[i] = t_0_xx_0_yzz[i] - rcd_y[i] * t_0_xx_0_zz[i];

        t_0_xx_y_yz[i] = t_0_xx_0_yyz[i] - rcd_y[i] * t_0_xx_0_yz[i];

        t_0_xx_y_yy[i] = t_0_xx_0_yyy[i] - rcd_y[i] * t_0_xx_0_yy[i];

        t_0_xx_y_xz[i] = t_0_xx_0_xyz[i] - rcd_y[i] * t_0_xx_0_xz[i];

        t_0_xx_y_xy[i] = t_0_xx_0_xyy[i] - rcd_y[i] * t_0_xx_0_xy[i];

        t_0_xx_y_xx[i] = t_0_xx_0_xxy[i] - rcd_y[i] * t_0_xx_0_xx[i];

        t_0_xx_x_zz[i] = t_0_xx_0_xzz[i] - rcd_x[i] * t_0_xx_0_zz[i];

        t_0_xx_x_yz[i] = t_0_xx_0_xyz[i] - rcd_x[i] * t_0_xx_0_yz[i];

        t_0_xx_x_yy[i] = t_0_xx_0_xyy[i] - rcd_x[i] * t_0_xx_0_yy[i];

        t_0_xx_x_xz[i] = t_0_xx_0_xxz[i] - rcd_x[i] * t_0_xx_0_xz[i];

        t_0_xx_x_xy[i] = t_0_xx_0_xxy[i] - rcd_x[i] * t_0_xx_0_xy[i];

        t_0_xx_x_xx[i] = t_0_xx_0_xxx[i] - rcd_x[i] * t_0_xx_0_xx[i];
    }
}


} // derirec namespace
