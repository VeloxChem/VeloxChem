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
compHostHRRForSPDD_V0(      BufferHostXY<T>&      intsBufferSPDD,
                      const BufferHostX<int32_t>& intsIndexesSPDD,
                      const BufferHostXY<T>&      intsBufferSPPD,
                      const BufferHostX<int32_t>& intsIndexesSPPD,
                      const BufferHostXY<T>&      intsBufferSPPF,
                      const BufferHostX<int32_t>& intsIndexesSPPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SPDD) integral components

    t_0_z_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(0));

    t_0_z_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(1));

    t_0_z_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(2));

    t_0_z_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(3));

    t_0_z_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(4));

    t_0_z_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(5));

    t_0_z_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(6));

    t_0_z_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(7));

    t_0_z_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(8));

    t_0_z_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(9));

    t_0_z_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(10));

    t_0_z_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(11));

    t_0_z_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(12));

    t_0_z_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(13));

    t_0_z_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(14));

    t_0_z_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(15));

    t_0_z_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(16));

    t_0_z_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(17));

    t_0_z_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(18));

    t_0_z_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(19));

    t_0_z_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(20));

    t_0_z_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(21));

    t_0_z_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(22));

    t_0_z_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(23));

    t_0_z_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(24));

    t_0_z_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(25));

    t_0_z_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(26));

    t_0_z_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(27));

    t_0_z_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(28));

    t_0_z_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(29));

    t_0_z_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(30));

    t_0_z_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(31));

    t_0_z_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(32));

    t_0_z_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(33));

    t_0_z_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(34));

    t_0_z_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(35));

    t_0_y_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(36));

    t_0_y_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(37));

    t_0_y_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(38));

    t_0_y_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(39));

    t_0_y_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(40));

    t_0_y_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(41));

    t_0_y_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(42));

    t_0_y_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(43));

    t_0_y_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(44));

    t_0_y_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(45));

    t_0_y_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(46));

    t_0_y_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(47));

    t_0_y_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(48));

    t_0_y_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(49));

    t_0_y_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(50));

    t_0_y_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(51));

    t_0_y_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(52));

    t_0_y_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(53));

    t_0_y_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(54));

    t_0_y_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(55));

    t_0_y_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(56));

    t_0_y_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(57));

    t_0_y_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(58));

    t_0_y_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(59));

    t_0_y_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(60));

    t_0_y_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(61));

    t_0_y_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(62));

    t_0_y_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(63));

    t_0_y_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(64));

    t_0_y_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(65));

    t_0_y_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(66));

    t_0_y_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(67));

    t_0_y_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(68));

    t_0_y_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(69));

    t_0_y_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(70));

    t_0_y_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(71));

    t_0_x_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(72));

    t_0_x_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(73));

    t_0_x_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(74));

    t_0_x_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(75));

    t_0_x_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(76));

    t_0_x_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(77));

    t_0_x_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(78));

    t_0_x_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(79));

    t_0_x_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(80));

    t_0_x_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(81));

    t_0_x_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(82));

    t_0_x_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(83));

    t_0_x_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(84));

    t_0_x_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(85));

    t_0_x_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(86));

    t_0_x_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(87));

    t_0_x_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(88));

    t_0_x_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(89));

    t_0_x_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(90));

    t_0_x_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(91));

    t_0_x_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(92));

    t_0_x_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(93));

    t_0_x_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(94));

    t_0_x_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(95));

    t_0_x_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(96));

    t_0_x_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(97));

    t_0_x_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(98));

    t_0_x_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(99));

    t_0_x_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(100));

    t_0_x_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(101));

    t_0_x_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(102));

    t_0_x_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(103));

    t_0_x_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(104));

    t_0_x_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(105));

    t_0_x_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(106));

    t_0_x_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(107));

    // set up (SPPD) integral components

    t_0_z_z_zz = intsBufferSPPD.data(intsIndexesSPPD(0));

    t_0_z_z_yz = intsBufferSPPD.data(intsIndexesSPPD(1));

    t_0_z_z_yy = intsBufferSPPD.data(intsIndexesSPPD(2));

    t_0_z_z_xz = intsBufferSPPD.data(intsIndexesSPPD(3));

    t_0_z_z_xy = intsBufferSPPD.data(intsIndexesSPPD(4));

    t_0_z_z_xx = intsBufferSPPD.data(intsIndexesSPPD(5));

    t_0_z_y_zz = intsBufferSPPD.data(intsIndexesSPPD(6));

    t_0_z_y_yz = intsBufferSPPD.data(intsIndexesSPPD(7));

    t_0_z_y_yy = intsBufferSPPD.data(intsIndexesSPPD(8));

    t_0_z_y_xz = intsBufferSPPD.data(intsIndexesSPPD(9));

    t_0_z_y_xy = intsBufferSPPD.data(intsIndexesSPPD(10));

    t_0_z_y_xx = intsBufferSPPD.data(intsIndexesSPPD(11));

    t_0_z_x_zz = intsBufferSPPD.data(intsIndexesSPPD(12));

    t_0_z_x_yz = intsBufferSPPD.data(intsIndexesSPPD(13));

    t_0_z_x_yy = intsBufferSPPD.data(intsIndexesSPPD(14));

    t_0_z_x_xz = intsBufferSPPD.data(intsIndexesSPPD(15));

    t_0_z_x_xy = intsBufferSPPD.data(intsIndexesSPPD(16));

    t_0_z_x_xx = intsBufferSPPD.data(intsIndexesSPPD(17));

    t_0_y_z_zz = intsBufferSPPD.data(intsIndexesSPPD(18));

    t_0_y_z_yz = intsBufferSPPD.data(intsIndexesSPPD(19));

    t_0_y_z_yy = intsBufferSPPD.data(intsIndexesSPPD(20));

    t_0_y_z_xz = intsBufferSPPD.data(intsIndexesSPPD(21));

    t_0_y_z_xy = intsBufferSPPD.data(intsIndexesSPPD(22));

    t_0_y_z_xx = intsBufferSPPD.data(intsIndexesSPPD(23));

    t_0_y_y_zz = intsBufferSPPD.data(intsIndexesSPPD(24));

    t_0_y_y_yz = intsBufferSPPD.data(intsIndexesSPPD(25));

    t_0_y_y_yy = intsBufferSPPD.data(intsIndexesSPPD(26));

    t_0_y_y_xz = intsBufferSPPD.data(intsIndexesSPPD(27));

    t_0_y_y_xy = intsBufferSPPD.data(intsIndexesSPPD(28));

    t_0_y_y_xx = intsBufferSPPD.data(intsIndexesSPPD(29));

    t_0_y_x_zz = intsBufferSPPD.data(intsIndexesSPPD(30));

    t_0_y_x_yz = intsBufferSPPD.data(intsIndexesSPPD(31));

    t_0_y_x_yy = intsBufferSPPD.data(intsIndexesSPPD(32));

    t_0_y_x_xz = intsBufferSPPD.data(intsIndexesSPPD(33));

    t_0_y_x_xy = intsBufferSPPD.data(intsIndexesSPPD(34));

    t_0_y_x_xx = intsBufferSPPD.data(intsIndexesSPPD(35));

    t_0_x_z_zz = intsBufferSPPD.data(intsIndexesSPPD(36));

    t_0_x_z_yz = intsBufferSPPD.data(intsIndexesSPPD(37));

    t_0_x_z_yy = intsBufferSPPD.data(intsIndexesSPPD(38));

    t_0_x_z_xz = intsBufferSPPD.data(intsIndexesSPPD(39));

    t_0_x_z_xy = intsBufferSPPD.data(intsIndexesSPPD(40));

    t_0_x_z_xx = intsBufferSPPD.data(intsIndexesSPPD(41));

    t_0_x_y_zz = intsBufferSPPD.data(intsIndexesSPPD(42));

    t_0_x_y_yz = intsBufferSPPD.data(intsIndexesSPPD(43));

    t_0_x_y_yy = intsBufferSPPD.data(intsIndexesSPPD(44));

    t_0_x_y_xz = intsBufferSPPD.data(intsIndexesSPPD(45));

    t_0_x_y_xy = intsBufferSPPD.data(intsIndexesSPPD(46));

    t_0_x_y_xx = intsBufferSPPD.data(intsIndexesSPPD(47));

    t_0_x_x_zz = intsBufferSPPD.data(intsIndexesSPPD(48));

    t_0_x_x_yz = intsBufferSPPD.data(intsIndexesSPPD(49));

    t_0_x_x_yy = intsBufferSPPD.data(intsIndexesSPPD(50));

    t_0_x_x_xz = intsBufferSPPD.data(intsIndexesSPPD(51));

    t_0_x_x_xy = intsBufferSPPD.data(intsIndexesSPPD(52));

    t_0_x_x_xx = intsBufferSPPD.data(intsIndexesSPPD(53));

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

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_z_x_xx, t_0_z_x_xxx, t_0_z_x_xxy,\
                           t_0_z_x_xxz, t_0_z_x_xy, t_0_z_x_xyy, t_0_z_x_xyz, t_0_z_x_xz,\
                           t_0_z_x_xzz, t_0_z_x_yy, t_0_z_x_yyy, t_0_z_x_yyz, t_0_z_x_yz,\
                           t_0_z_x_yzz, t_0_z_x_zz, t_0_z_x_zzz, t_0_z_xx_xx, t_0_z_xx_xy,\
                           t_0_z_xx_xz, t_0_z_xx_yy, t_0_z_xx_yz, t_0_z_xx_zz, t_0_z_xy_xx,\
                           t_0_z_xy_xy, t_0_z_xy_xz, t_0_z_xy_yy, t_0_z_xy_yz, t_0_z_xy_zz,\
                           t_0_z_xz_xx, t_0_z_xz_xy, t_0_z_xz_xz, t_0_z_xz_yy, t_0_z_xz_yz,\
                           t_0_z_xz_zz, t_0_z_y_xx, t_0_z_y_xxy, t_0_z_y_xxz, t_0_z_y_xy,\
                           t_0_z_y_xyy, t_0_z_y_xyz, t_0_z_y_xz, t_0_z_y_xzz, t_0_z_y_yy,\
                           t_0_z_y_yyy, t_0_z_y_yyz, t_0_z_y_yz, t_0_z_y_yzz, t_0_z_y_zz,\
                           t_0_z_y_zzz, t_0_z_yy_xx, t_0_z_yy_xy, t_0_z_yy_xz, t_0_z_yy_yy,\
                           t_0_z_yy_yz, t_0_z_yy_zz, t_0_z_yz_xx, t_0_z_yz_xy, t_0_z_yz_xz,\
                           t_0_z_yz_yy, t_0_z_yz_yz, t_0_z_yz_zz, t_0_z_z_xx, t_0_z_z_xxz,\
                           t_0_z_z_xy, t_0_z_z_xyz, t_0_z_z_xz, t_0_z_z_xzz, t_0_z_z_yy,\
                           t_0_z_z_yyz, t_0_z_z_yz, t_0_z_z_yzz, t_0_z_z_zz, t_0_z_z_zzz,\
                           t_0_z_zz_xx, t_0_z_zz_xy, t_0_z_zz_xz, t_0_z_zz_yy, t_0_z_zz_yz,\
                           t_0_z_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_z_zz_zz[i] = t_0_z_z_zzz[i] - rcd_z[i] * t_0_z_z_zz[i];

        t_0_z_zz_yz[i] = t_0_z_z_yzz[i] - rcd_z[i] * t_0_z_z_yz[i];

        t_0_z_zz_yy[i] = t_0_z_z_yyz[i] - rcd_z[i] * t_0_z_z_yy[i];

        t_0_z_zz_xz[i] = t_0_z_z_xzz[i] - rcd_z[i] * t_0_z_z_xz[i];

        t_0_z_zz_xy[i] = t_0_z_z_xyz[i] - rcd_z[i] * t_0_z_z_xy[i];

        t_0_z_zz_xx[i] = t_0_z_z_xxz[i] - rcd_z[i] * t_0_z_z_xx[i];

        t_0_z_yz_zz[i] = t_0_z_y_zzz[i] - rcd_z[i] * t_0_z_y_zz[i];

        t_0_z_yz_yz[i] = t_0_z_y_yzz[i] - rcd_z[i] * t_0_z_y_yz[i];

        t_0_z_yz_yy[i] = t_0_z_y_yyz[i] - rcd_z[i] * t_0_z_y_yy[i];

        t_0_z_yz_xz[i] = t_0_z_y_xzz[i] - rcd_z[i] * t_0_z_y_xz[i];

        t_0_z_yz_xy[i] = t_0_z_y_xyz[i] - rcd_z[i] * t_0_z_y_xy[i];

        t_0_z_yz_xx[i] = t_0_z_y_xxz[i] - rcd_z[i] * t_0_z_y_xx[i];

        t_0_z_yy_zz[i] = t_0_z_y_yzz[i] - rcd_y[i] * t_0_z_y_zz[i];

        t_0_z_yy_yz[i] = t_0_z_y_yyz[i] - rcd_y[i] * t_0_z_y_yz[i];

        t_0_z_yy_yy[i] = t_0_z_y_yyy[i] - rcd_y[i] * t_0_z_y_yy[i];

        t_0_z_yy_xz[i] = t_0_z_y_xyz[i] - rcd_y[i] * t_0_z_y_xz[i];

        t_0_z_yy_xy[i] = t_0_z_y_xyy[i] - rcd_y[i] * t_0_z_y_xy[i];

        t_0_z_yy_xx[i] = t_0_z_y_xxy[i] - rcd_y[i] * t_0_z_y_xx[i];

        t_0_z_xz_zz[i] = t_0_z_x_zzz[i] - rcd_z[i] * t_0_z_x_zz[i];

        t_0_z_xz_yz[i] = t_0_z_x_yzz[i] - rcd_z[i] * t_0_z_x_yz[i];

        t_0_z_xz_yy[i] = t_0_z_x_yyz[i] - rcd_z[i] * t_0_z_x_yy[i];

        t_0_z_xz_xz[i] = t_0_z_x_xzz[i] - rcd_z[i] * t_0_z_x_xz[i];

        t_0_z_xz_xy[i] = t_0_z_x_xyz[i] - rcd_z[i] * t_0_z_x_xy[i];

        t_0_z_xz_xx[i] = t_0_z_x_xxz[i] - rcd_z[i] * t_0_z_x_xx[i];

        t_0_z_xy_zz[i] = t_0_z_x_yzz[i] - rcd_y[i] * t_0_z_x_zz[i];

        t_0_z_xy_yz[i] = t_0_z_x_yyz[i] - rcd_y[i] * t_0_z_x_yz[i];

        t_0_z_xy_yy[i] = t_0_z_x_yyy[i] - rcd_y[i] * t_0_z_x_yy[i];

        t_0_z_xy_xz[i] = t_0_z_x_xyz[i] - rcd_y[i] * t_0_z_x_xz[i];

        t_0_z_xy_xy[i] = t_0_z_x_xyy[i] - rcd_y[i] * t_0_z_x_xy[i];

        t_0_z_xy_xx[i] = t_0_z_x_xxy[i] - rcd_y[i] * t_0_z_x_xx[i];

        t_0_z_xx_zz[i] = t_0_z_x_xzz[i] - rcd_x[i] * t_0_z_x_zz[i];

        t_0_z_xx_yz[i] = t_0_z_x_xyz[i] - rcd_x[i] * t_0_z_x_yz[i];

        t_0_z_xx_yy[i] = t_0_z_x_xyy[i] - rcd_x[i] * t_0_z_x_yy[i];

        t_0_z_xx_xz[i] = t_0_z_x_xxz[i] - rcd_x[i] * t_0_z_x_xz[i];

        t_0_z_xx_xy[i] = t_0_z_x_xxy[i] - rcd_x[i] * t_0_z_x_xy[i];

        t_0_z_xx_xx[i] = t_0_z_x_xxx[i] - rcd_x[i] * t_0_z_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_y_x_xx, t_0_y_x_xxx, t_0_y_x_xxy,\
                           t_0_y_x_xxz, t_0_y_x_xy, t_0_y_x_xyy, t_0_y_x_xyz, t_0_y_x_xz,\
                           t_0_y_x_xzz, t_0_y_x_yy, t_0_y_x_yyy, t_0_y_x_yyz, t_0_y_x_yz,\
                           t_0_y_x_yzz, t_0_y_x_zz, t_0_y_x_zzz, t_0_y_xx_xx, t_0_y_xx_xy,\
                           t_0_y_xx_xz, t_0_y_xx_yy, t_0_y_xx_yz, t_0_y_xx_zz, t_0_y_xy_xx,\
                           t_0_y_xy_xy, t_0_y_xy_xz, t_0_y_xy_yy, t_0_y_xy_yz, t_0_y_xy_zz,\
                           t_0_y_xz_xx, t_0_y_xz_xy, t_0_y_xz_xz, t_0_y_xz_yy, t_0_y_xz_yz,\
                           t_0_y_xz_zz, t_0_y_y_xx, t_0_y_y_xxy, t_0_y_y_xxz, t_0_y_y_xy,\
                           t_0_y_y_xyy, t_0_y_y_xyz, t_0_y_y_xz, t_0_y_y_xzz, t_0_y_y_yy,\
                           t_0_y_y_yyy, t_0_y_y_yyz, t_0_y_y_yz, t_0_y_y_yzz, t_0_y_y_zz,\
                           t_0_y_y_zzz, t_0_y_yy_xx, t_0_y_yy_xy, t_0_y_yy_xz, t_0_y_yy_yy,\
                           t_0_y_yy_yz, t_0_y_yy_zz, t_0_y_yz_xx, t_0_y_yz_xy, t_0_y_yz_xz,\
                           t_0_y_yz_yy, t_0_y_yz_yz, t_0_y_yz_zz, t_0_y_z_xx, t_0_y_z_xxz,\
                           t_0_y_z_xy, t_0_y_z_xyz, t_0_y_z_xz, t_0_y_z_xzz, t_0_y_z_yy,\
                           t_0_y_z_yyz, t_0_y_z_yz, t_0_y_z_yzz, t_0_y_z_zz, t_0_y_z_zzz,\
                           t_0_y_zz_xx, t_0_y_zz_xy, t_0_y_zz_xz, t_0_y_zz_yy, t_0_y_zz_yz,\
                           t_0_y_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_y_zz_zz[i] = t_0_y_z_zzz[i] - rcd_z[i] * t_0_y_z_zz[i];

        t_0_y_zz_yz[i] = t_0_y_z_yzz[i] - rcd_z[i] * t_0_y_z_yz[i];

        t_0_y_zz_yy[i] = t_0_y_z_yyz[i] - rcd_z[i] * t_0_y_z_yy[i];

        t_0_y_zz_xz[i] = t_0_y_z_xzz[i] - rcd_z[i] * t_0_y_z_xz[i];

        t_0_y_zz_xy[i] = t_0_y_z_xyz[i] - rcd_z[i] * t_0_y_z_xy[i];

        t_0_y_zz_xx[i] = t_0_y_z_xxz[i] - rcd_z[i] * t_0_y_z_xx[i];

        t_0_y_yz_zz[i] = t_0_y_y_zzz[i] - rcd_z[i] * t_0_y_y_zz[i];

        t_0_y_yz_yz[i] = t_0_y_y_yzz[i] - rcd_z[i] * t_0_y_y_yz[i];

        t_0_y_yz_yy[i] = t_0_y_y_yyz[i] - rcd_z[i] * t_0_y_y_yy[i];

        t_0_y_yz_xz[i] = t_0_y_y_xzz[i] - rcd_z[i] * t_0_y_y_xz[i];

        t_0_y_yz_xy[i] = t_0_y_y_xyz[i] - rcd_z[i] * t_0_y_y_xy[i];

        t_0_y_yz_xx[i] = t_0_y_y_xxz[i] - rcd_z[i] * t_0_y_y_xx[i];

        t_0_y_yy_zz[i] = t_0_y_y_yzz[i] - rcd_y[i] * t_0_y_y_zz[i];

        t_0_y_yy_yz[i] = t_0_y_y_yyz[i] - rcd_y[i] * t_0_y_y_yz[i];

        t_0_y_yy_yy[i] = t_0_y_y_yyy[i] - rcd_y[i] * t_0_y_y_yy[i];

        t_0_y_yy_xz[i] = t_0_y_y_xyz[i] - rcd_y[i] * t_0_y_y_xz[i];

        t_0_y_yy_xy[i] = t_0_y_y_xyy[i] - rcd_y[i] * t_0_y_y_xy[i];

        t_0_y_yy_xx[i] = t_0_y_y_xxy[i] - rcd_y[i] * t_0_y_y_xx[i];

        t_0_y_xz_zz[i] = t_0_y_x_zzz[i] - rcd_z[i] * t_0_y_x_zz[i];

        t_0_y_xz_yz[i] = t_0_y_x_yzz[i] - rcd_z[i] * t_0_y_x_yz[i];

        t_0_y_xz_yy[i] = t_0_y_x_yyz[i] - rcd_z[i] * t_0_y_x_yy[i];

        t_0_y_xz_xz[i] = t_0_y_x_xzz[i] - rcd_z[i] * t_0_y_x_xz[i];

        t_0_y_xz_xy[i] = t_0_y_x_xyz[i] - rcd_z[i] * t_0_y_x_xy[i];

        t_0_y_xz_xx[i] = t_0_y_x_xxz[i] - rcd_z[i] * t_0_y_x_xx[i];

        t_0_y_xy_zz[i] = t_0_y_x_yzz[i] - rcd_y[i] * t_0_y_x_zz[i];

        t_0_y_xy_yz[i] = t_0_y_x_yyz[i] - rcd_y[i] * t_0_y_x_yz[i];

        t_0_y_xy_yy[i] = t_0_y_x_yyy[i] - rcd_y[i] * t_0_y_x_yy[i];

        t_0_y_xy_xz[i] = t_0_y_x_xyz[i] - rcd_y[i] * t_0_y_x_xz[i];

        t_0_y_xy_xy[i] = t_0_y_x_xyy[i] - rcd_y[i] * t_0_y_x_xy[i];

        t_0_y_xy_xx[i] = t_0_y_x_xxy[i] - rcd_y[i] * t_0_y_x_xx[i];

        t_0_y_xx_zz[i] = t_0_y_x_xzz[i] - rcd_x[i] * t_0_y_x_zz[i];

        t_0_y_xx_yz[i] = t_0_y_x_xyz[i] - rcd_x[i] * t_0_y_x_yz[i];

        t_0_y_xx_yy[i] = t_0_y_x_xyy[i] - rcd_x[i] * t_0_y_x_yy[i];

        t_0_y_xx_xz[i] = t_0_y_x_xxz[i] - rcd_x[i] * t_0_y_x_xz[i];

        t_0_y_xx_xy[i] = t_0_y_x_xxy[i] - rcd_x[i] * t_0_y_x_xy[i];

        t_0_y_xx_xx[i] = t_0_y_x_xxx[i] - rcd_x[i] * t_0_y_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_x_x_xx, t_0_x_x_xxx, t_0_x_x_xxy,\
                           t_0_x_x_xxz, t_0_x_x_xy, t_0_x_x_xyy, t_0_x_x_xyz, t_0_x_x_xz,\
                           t_0_x_x_xzz, t_0_x_x_yy, t_0_x_x_yyy, t_0_x_x_yyz, t_0_x_x_yz,\
                           t_0_x_x_yzz, t_0_x_x_zz, t_0_x_x_zzz, t_0_x_xx_xx, t_0_x_xx_xy,\
                           t_0_x_xx_xz, t_0_x_xx_yy, t_0_x_xx_yz, t_0_x_xx_zz, t_0_x_xy_xx,\
                           t_0_x_xy_xy, t_0_x_xy_xz, t_0_x_xy_yy, t_0_x_xy_yz, t_0_x_xy_zz,\
                           t_0_x_xz_xx, t_0_x_xz_xy, t_0_x_xz_xz, t_0_x_xz_yy, t_0_x_xz_yz,\
                           t_0_x_xz_zz, t_0_x_y_xx, t_0_x_y_xxy, t_0_x_y_xxz, t_0_x_y_xy,\
                           t_0_x_y_xyy, t_0_x_y_xyz, t_0_x_y_xz, t_0_x_y_xzz, t_0_x_y_yy,\
                           t_0_x_y_yyy, t_0_x_y_yyz, t_0_x_y_yz, t_0_x_y_yzz, t_0_x_y_zz,\
                           t_0_x_y_zzz, t_0_x_yy_xx, t_0_x_yy_xy, t_0_x_yy_xz, t_0_x_yy_yy,\
                           t_0_x_yy_yz, t_0_x_yy_zz, t_0_x_yz_xx, t_0_x_yz_xy, t_0_x_yz_xz,\
                           t_0_x_yz_yy, t_0_x_yz_yz, t_0_x_yz_zz, t_0_x_z_xx, t_0_x_z_xxz,\
                           t_0_x_z_xy, t_0_x_z_xyz, t_0_x_z_xz, t_0_x_z_xzz, t_0_x_z_yy,\
                           t_0_x_z_yyz, t_0_x_z_yz, t_0_x_z_yzz, t_0_x_z_zz, t_0_x_z_zzz,\
                           t_0_x_zz_xx, t_0_x_zz_xy, t_0_x_zz_xz, t_0_x_zz_yy, t_0_x_zz_yz,\
                           t_0_x_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_x_zz_zz[i] = t_0_x_z_zzz[i] - rcd_z[i] * t_0_x_z_zz[i];

        t_0_x_zz_yz[i] = t_0_x_z_yzz[i] - rcd_z[i] * t_0_x_z_yz[i];

        t_0_x_zz_yy[i] = t_0_x_z_yyz[i] - rcd_z[i] * t_0_x_z_yy[i];

        t_0_x_zz_xz[i] = t_0_x_z_xzz[i] - rcd_z[i] * t_0_x_z_xz[i];

        t_0_x_zz_xy[i] = t_0_x_z_xyz[i] - rcd_z[i] * t_0_x_z_xy[i];

        t_0_x_zz_xx[i] = t_0_x_z_xxz[i] - rcd_z[i] * t_0_x_z_xx[i];

        t_0_x_yz_zz[i] = t_0_x_y_zzz[i] - rcd_z[i] * t_0_x_y_zz[i];

        t_0_x_yz_yz[i] = t_0_x_y_yzz[i] - rcd_z[i] * t_0_x_y_yz[i];

        t_0_x_yz_yy[i] = t_0_x_y_yyz[i] - rcd_z[i] * t_0_x_y_yy[i];

        t_0_x_yz_xz[i] = t_0_x_y_xzz[i] - rcd_z[i] * t_0_x_y_xz[i];

        t_0_x_yz_xy[i] = t_0_x_y_xyz[i] - rcd_z[i] * t_0_x_y_xy[i];

        t_0_x_yz_xx[i] = t_0_x_y_xxz[i] - rcd_z[i] * t_0_x_y_xx[i];

        t_0_x_yy_zz[i] = t_0_x_y_yzz[i] - rcd_y[i] * t_0_x_y_zz[i];

        t_0_x_yy_yz[i] = t_0_x_y_yyz[i] - rcd_y[i] * t_0_x_y_yz[i];

        t_0_x_yy_yy[i] = t_0_x_y_yyy[i] - rcd_y[i] * t_0_x_y_yy[i];

        t_0_x_yy_xz[i] = t_0_x_y_xyz[i] - rcd_y[i] * t_0_x_y_xz[i];

        t_0_x_yy_xy[i] = t_0_x_y_xyy[i] - rcd_y[i] * t_0_x_y_xy[i];

        t_0_x_yy_xx[i] = t_0_x_y_xxy[i] - rcd_y[i] * t_0_x_y_xx[i];

        t_0_x_xz_zz[i] = t_0_x_x_zzz[i] - rcd_z[i] * t_0_x_x_zz[i];

        t_0_x_xz_yz[i] = t_0_x_x_yzz[i] - rcd_z[i] * t_0_x_x_yz[i];

        t_0_x_xz_yy[i] = t_0_x_x_yyz[i] - rcd_z[i] * t_0_x_x_yy[i];

        t_0_x_xz_xz[i] = t_0_x_x_xzz[i] - rcd_z[i] * t_0_x_x_xz[i];

        t_0_x_xz_xy[i] = t_0_x_x_xyz[i] - rcd_z[i] * t_0_x_x_xy[i];

        t_0_x_xz_xx[i] = t_0_x_x_xxz[i] - rcd_z[i] * t_0_x_x_xx[i];

        t_0_x_xy_zz[i] = t_0_x_x_yzz[i] - rcd_y[i] * t_0_x_x_zz[i];

        t_0_x_xy_yz[i] = t_0_x_x_yyz[i] - rcd_y[i] * t_0_x_x_yz[i];

        t_0_x_xy_yy[i] = t_0_x_x_yyy[i] - rcd_y[i] * t_0_x_x_yy[i];

        t_0_x_xy_xz[i] = t_0_x_x_xyz[i] - rcd_y[i] * t_0_x_x_xz[i];

        t_0_x_xy_xy[i] = t_0_x_x_xyy[i] - rcd_y[i] * t_0_x_x_xy[i];

        t_0_x_xy_xx[i] = t_0_x_x_xxy[i] - rcd_y[i] * t_0_x_x_xx[i];

        t_0_x_xx_zz[i] = t_0_x_x_xzz[i] - rcd_x[i] * t_0_x_x_zz[i];

        t_0_x_xx_yz[i] = t_0_x_x_xyz[i] - rcd_x[i] * t_0_x_x_yz[i];

        t_0_x_xx_yy[i] = t_0_x_x_xyy[i] - rcd_x[i] * t_0_x_x_yy[i];

        t_0_x_xx_xz[i] = t_0_x_x_xxz[i] - rcd_x[i] * t_0_x_x_xz[i];

        t_0_x_xx_xy[i] = t_0_x_x_xxy[i] - rcd_x[i] * t_0_x_x_xy[i];

        t_0_x_xx_xx[i] = t_0_x_x_xxx[i] - rcd_x[i] * t_0_x_x_xx[i];
    }
}


} // derirec namespace
