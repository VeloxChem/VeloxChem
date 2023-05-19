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

    t_0_zz_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(1));

    t_0_zz_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(2));

    t_0_zz_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(3));

    t_0_zz_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(4));

    t_0_zz_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(5));

    t_0_zz_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(6));

    t_0_zz_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(7));

    t_0_zz_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(8));

    t_0_zz_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(9));

    t_0_zz_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(10));

    t_0_zz_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(11));

    t_0_zz_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(12));

    t_0_zz_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(13));

    t_0_zz_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(14));

    t_0_zz_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(15));

    t_0_zz_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(16));

    t_0_zz_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(17));

    t_0_zz_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(18));

    t_0_zz_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(19));

    t_0_zz_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(20));

    t_0_zz_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(21));

    t_0_zz_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(22));

    t_0_zz_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(23));

    t_0_zz_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(24));

    t_0_yz_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(25));

    t_0_yz_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(26));

    t_0_yz_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(27));

    t_0_yz_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(28));

    t_0_yz_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(29));

    t_0_yz_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(30));

    t_0_yz_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(31));

    t_0_yz_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(32));

    t_0_yz_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(33));

    t_0_yz_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(34));

    t_0_yz_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(35));

    t_0_yz_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(36));

    t_0_yz_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(37));

    t_0_yz_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(38));

    t_0_yz_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(39));

    t_0_yz_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(40));

    t_0_yz_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(41));

    t_0_yz_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(42));

    t_0_yz_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(43));

    t_0_yz_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(44));

    t_0_yz_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(45));

    t_0_yz_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(46));

    t_0_yz_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(47));

    t_0_yz_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(48));

    t_0_yz_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(49));

    t_0_yy_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(50));

    t_0_yy_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(51));

    t_0_yy_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(52));

    t_0_yy_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(53));

    t_0_yy_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(54));

    t_0_yy_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(55));

    t_0_yy_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(56));

    t_0_yy_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(57));

    t_0_yy_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(58));

    t_0_yy_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(59));

    t_0_yy_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(60));

    t_0_yy_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(61));

    t_0_yy_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(62));

    t_0_yy_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(63));

    t_0_yy_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(64));

    t_0_yy_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(65));

    t_0_yy_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(66));

    t_0_yy_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(67));

    t_0_yy_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(68));

    t_0_yy_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(69));

    t_0_yy_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(70));

    t_0_yy_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(71));

    t_0_yy_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(72));

    t_0_yy_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(73));

    t_0_yy_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(74));

    t_0_xz_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(75));

    t_0_xz_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(76));

    t_0_xz_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(77));

    t_0_xz_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(78));

    t_0_xz_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(79));

    t_0_xz_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(80));

    t_0_xz_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(81));

    t_0_xz_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(82));

    t_0_xz_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(83));

    t_0_xz_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(84));

    t_0_xz_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(85));

    t_0_xz_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(86));

    t_0_xz_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(87));

    t_0_xz_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(88));

    t_0_xz_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(89));

    t_0_xz_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(90));

    t_0_xz_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(91));

    t_0_xz_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(92));

    t_0_xz_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(93));

    t_0_xz_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(94));

    t_0_xz_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(95));

    t_0_xz_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(96));

    t_0_xz_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(97));

    t_0_xz_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(98));

    t_0_xz_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(99));

    t_0_xy_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(100));

    t_0_xy_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(101));

    t_0_xy_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(102));

    t_0_xy_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(103));

    t_0_xy_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(104));

    t_0_xy_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(105));

    t_0_xy_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(106));

    t_0_xy_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(107));

    t_0_xy_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(108));

    t_0_xy_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(109));

    t_0_xy_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(110));

    t_0_xy_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(111));

    t_0_xy_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(112));

    t_0_xy_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(113));

    t_0_xy_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(114));

    t_0_xy_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(115));

    t_0_xy_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(116));

    t_0_xy_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(117));

    t_0_xy_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(118));

    t_0_xy_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(119));

    t_0_xy_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(120));

    t_0_xy_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(121));

    t_0_xy_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(122));

    t_0_xy_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(123));

    t_0_xy_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(124));

    t_0_xx_z_zzz = intsBufferSDPF.data(intsIndexesSDPF(125));

    t_0_xx_z_yzz = intsBufferSDPF.data(intsIndexesSDPF(126));

    t_0_xx_z_yyz = intsBufferSDPF.data(intsIndexesSDPF(127));

    t_0_xx_z_xzz = intsBufferSDPF.data(intsIndexesSDPF(128));

    t_0_xx_z_xyz = intsBufferSDPF.data(intsIndexesSDPF(129));

    t_0_xx_z_xxz = intsBufferSDPF.data(intsIndexesSDPF(130));

    t_0_xx_y_zzz = intsBufferSDPF.data(intsIndexesSDPF(131));

    t_0_xx_y_yzz = intsBufferSDPF.data(intsIndexesSDPF(132));

    t_0_xx_y_yyz = intsBufferSDPF.data(intsIndexesSDPF(133));

    t_0_xx_y_yyy = intsBufferSDPF.data(intsIndexesSDPF(134));

    t_0_xx_y_xzz = intsBufferSDPF.data(intsIndexesSDPF(135));

    t_0_xx_y_xyz = intsBufferSDPF.data(intsIndexesSDPF(136));

    t_0_xx_y_xyy = intsBufferSDPF.data(intsIndexesSDPF(137));

    t_0_xx_y_xxz = intsBufferSDPF.data(intsIndexesSDPF(138));

    t_0_xx_y_xxy = intsBufferSDPF.data(intsIndexesSDPF(139));

    t_0_xx_x_zzz = intsBufferSDPF.data(intsIndexesSDPF(140));

    t_0_xx_x_yzz = intsBufferSDPF.data(intsIndexesSDPF(141));

    t_0_xx_x_yyz = intsBufferSDPF.data(intsIndexesSDPF(142));

    t_0_xx_x_yyy = intsBufferSDPF.data(intsIndexesSDPF(143));

    t_0_xx_x_xzz = intsBufferSDPF.data(intsIndexesSDPF(144));

    t_0_xx_x_xyz = intsBufferSDPF.data(intsIndexesSDPF(145));

    t_0_xx_x_xyy = intsBufferSDPF.data(intsIndexesSDPF(146));

    t_0_xx_x_xxz = intsBufferSDPF.data(intsIndexesSDPF(147));

    t_0_xx_x_xxy = intsBufferSDPF.data(intsIndexesSDPF(148));

    t_0_xx_x_xxx = intsBufferSDPF.data(intsIndexesSDPF(149));

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

    // set up (SDSG) integral components

    t_0_zz_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(0));

    t_0_zz_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(1));

    t_0_zz_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(2));

    t_0_zz_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(3));

    t_0_zz_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(4));

    t_0_zz_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(5));

    t_0_zz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(6));

    t_0_zz_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(7));

    t_0_zz_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(8));

    t_0_zz_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(9));

    t_0_zz_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(10));

    t_0_zz_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(11));

    t_0_zz_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(12));

    t_0_zz_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(13));

    t_0_zz_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(14));

    t_0_yz_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(15));

    t_0_yz_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(16));

    t_0_yz_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(17));

    t_0_yz_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(18));

    t_0_yz_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(19));

    t_0_yz_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(20));

    t_0_yz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(21));

    t_0_yz_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(22));

    t_0_yz_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(23));

    t_0_yz_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(24));

    t_0_yz_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(25));

    t_0_yz_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(26));

    t_0_yz_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(27));

    t_0_yz_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(28));

    t_0_yz_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(29));

    t_0_yy_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(30));

    t_0_yy_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(31));

    t_0_yy_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(32));

    t_0_yy_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(33));

    t_0_yy_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(34));

    t_0_yy_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(35));

    t_0_yy_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(36));

    t_0_yy_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(37));

    t_0_yy_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(38));

    t_0_yy_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(39));

    t_0_yy_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(40));

    t_0_yy_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(41));

    t_0_yy_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(42));

    t_0_yy_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(43));

    t_0_yy_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(44));

    t_0_xz_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(45));

    t_0_xz_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(46));

    t_0_xz_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(47));

    t_0_xz_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(48));

    t_0_xz_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(49));

    t_0_xz_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(50));

    t_0_xz_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(51));

    t_0_xz_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(52));

    t_0_xz_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(53));

    t_0_xz_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(54));

    t_0_xz_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(55));

    t_0_xz_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(56));

    t_0_xz_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(57));

    t_0_xz_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(58));

    t_0_xz_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(59));

    t_0_xy_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(60));

    t_0_xy_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(61));

    t_0_xy_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(62));

    t_0_xy_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(63));

    t_0_xy_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(64));

    t_0_xy_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(65));

    t_0_xy_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(66));

    t_0_xy_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(67));

    t_0_xy_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(68));

    t_0_xy_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(69));

    t_0_xy_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(70));

    t_0_xy_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(71));

    t_0_xy_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(72));

    t_0_xy_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(73));

    t_0_xy_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(74));

    t_0_xx_0_zzzz = intsBufferSDSG.data(intsIndexesSDSG(75));

    t_0_xx_0_yzzz = intsBufferSDSG.data(intsIndexesSDSG(76));

    t_0_xx_0_yyzz = intsBufferSDSG.data(intsIndexesSDSG(77));

    t_0_xx_0_yyyz = intsBufferSDSG.data(intsIndexesSDSG(78));

    t_0_xx_0_yyyy = intsBufferSDSG.data(intsIndexesSDSG(79));

    t_0_xx_0_xzzz = intsBufferSDSG.data(intsIndexesSDSG(80));

    t_0_xx_0_xyzz = intsBufferSDSG.data(intsIndexesSDSG(81));

    t_0_xx_0_xyyz = intsBufferSDSG.data(intsIndexesSDSG(82));

    t_0_xx_0_xyyy = intsBufferSDSG.data(intsIndexesSDSG(83));

    t_0_xx_0_xxzz = intsBufferSDSG.data(intsIndexesSDSG(84));

    t_0_xx_0_xxyz = intsBufferSDSG.data(intsIndexesSDSG(85));

    t_0_xx_0_xxyy = intsBufferSDSG.data(intsIndexesSDSG(86));

    t_0_xx_0_xxxz = intsBufferSDSG.data(intsIndexesSDSG(87));

    t_0_xx_0_xxxy = intsBufferSDSG.data(intsIndexesSDSG(88));

    t_0_xx_0_xxxx = intsBufferSDSG.data(intsIndexesSDSG(89));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yz_0_xxz, t_0_yz_0_xxzz, t_0_yz_0_xyz,\
                           t_0_yz_0_xyzz, t_0_yz_0_xzz, t_0_yz_0_xzzz, t_0_yz_0_yyy,\
                           t_0_yz_0_yyyy, t_0_yz_0_yyyz, t_0_yz_0_yyz, t_0_yz_0_yyzz,\
                           t_0_yz_0_yzz, t_0_yz_0_yzzz, t_0_yz_0_zzz, t_0_yz_0_zzzz,\
                           t_0_yz_y_xzz, t_0_yz_y_yyy, t_0_yz_y_yyz, t_0_yz_y_yzz, t_0_yz_y_zzz,\
                           t_0_yz_z_xxz, t_0_yz_z_xyz, t_0_yz_z_xzz, t_0_yz_z_yyz, t_0_yz_z_yzz,\
                           t_0_yz_z_zzz, t_0_zz_0_xxx, t_0_zz_0_xxxx, t_0_zz_0_xxxy,\
                           t_0_zz_0_xxxz, t_0_zz_0_xxy, t_0_zz_0_xxyy, t_0_zz_0_xxyz,\
                           t_0_zz_0_xxz, t_0_zz_0_xxzz, t_0_zz_0_xyy, t_0_zz_0_xyyy,\
                           t_0_zz_0_xyyz, t_0_zz_0_xyz, t_0_zz_0_xyzz, t_0_zz_0_xzz,\
                           t_0_zz_0_xzzz, t_0_zz_0_yyy, t_0_zz_0_yyyy, t_0_zz_0_yyyz,\
                           t_0_zz_0_yyz, t_0_zz_0_yyzz, t_0_zz_0_yzz, t_0_zz_0_yzzz,\
                           t_0_zz_0_zzz, t_0_zz_0_zzzz, t_0_zz_x_xxx, t_0_zz_x_xxy,\
                           t_0_zz_x_xxz, t_0_zz_x_xyy, t_0_zz_x_xyz, t_0_zz_x_xzz, t_0_zz_x_yyy,\
                           t_0_zz_x_yyz, t_0_zz_x_yzz, t_0_zz_x_zzz, t_0_zz_y_xxy, t_0_zz_y_xxz,\
                           t_0_zz_y_xyy, t_0_zz_y_xyz, t_0_zz_y_xzz, t_0_zz_y_yyy, t_0_zz_y_yyz,\
                           t_0_zz_y_yzz, t_0_zz_y_zzz, t_0_zz_z_xxz, t_0_zz_z_xyz, t_0_zz_z_xzz,\
                           t_0_zz_z_yyz, t_0_zz_z_yzz, t_0_zz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_z_zzz[i] = t_0_zz_0_zzzz[i] - rcd_z[i] * t_0_zz_0_zzz[i];

        t_0_zz_z_yzz[i] = t_0_zz_0_yzzz[i] - rcd_z[i] * t_0_zz_0_yzz[i];

        t_0_zz_z_yyz[i] = t_0_zz_0_yyzz[i] - rcd_z[i] * t_0_zz_0_yyz[i];

        t_0_zz_z_xzz[i] = t_0_zz_0_xzzz[i] - rcd_z[i] * t_0_zz_0_xzz[i];

        t_0_zz_z_xyz[i] = t_0_zz_0_xyzz[i] - rcd_z[i] * t_0_zz_0_xyz[i];

        t_0_zz_z_xxz[i] = t_0_zz_0_xxzz[i] - rcd_z[i] * t_0_zz_0_xxz[i];

        t_0_zz_y_zzz[i] = t_0_zz_0_yzzz[i] - rcd_y[i] * t_0_zz_0_zzz[i];

        t_0_zz_y_yzz[i] = t_0_zz_0_yyzz[i] - rcd_y[i] * t_0_zz_0_yzz[i];

        t_0_zz_y_yyz[i] = t_0_zz_0_yyyz[i] - rcd_y[i] * t_0_zz_0_yyz[i];

        t_0_zz_y_yyy[i] = t_0_zz_0_yyyy[i] - rcd_y[i] * t_0_zz_0_yyy[i];

        t_0_zz_y_xzz[i] = t_0_zz_0_xyzz[i] - rcd_y[i] * t_0_zz_0_xzz[i];

        t_0_zz_y_xyz[i] = t_0_zz_0_xyyz[i] - rcd_y[i] * t_0_zz_0_xyz[i];

        t_0_zz_y_xyy[i] = t_0_zz_0_xyyy[i] - rcd_y[i] * t_0_zz_0_xyy[i];

        t_0_zz_y_xxz[i] = t_0_zz_0_xxyz[i] - rcd_y[i] * t_0_zz_0_xxz[i];

        t_0_zz_y_xxy[i] = t_0_zz_0_xxyy[i] - rcd_y[i] * t_0_zz_0_xxy[i];

        t_0_zz_x_zzz[i] = t_0_zz_0_xzzz[i] - rcd_x[i] * t_0_zz_0_zzz[i];

        t_0_zz_x_yzz[i] = t_0_zz_0_xyzz[i] - rcd_x[i] * t_0_zz_0_yzz[i];

        t_0_zz_x_yyz[i] = t_0_zz_0_xyyz[i] - rcd_x[i] * t_0_zz_0_yyz[i];

        t_0_zz_x_yyy[i] = t_0_zz_0_xyyy[i] - rcd_x[i] * t_0_zz_0_yyy[i];

        t_0_zz_x_xzz[i] = t_0_zz_0_xxzz[i] - rcd_x[i] * t_0_zz_0_xzz[i];

        t_0_zz_x_xyz[i] = t_0_zz_0_xxyz[i] - rcd_x[i] * t_0_zz_0_xyz[i];

        t_0_zz_x_xyy[i] = t_0_zz_0_xxyy[i] - rcd_x[i] * t_0_zz_0_xyy[i];

        t_0_zz_x_xxz[i] = t_0_zz_0_xxxz[i] - rcd_x[i] * t_0_zz_0_xxz[i];

        t_0_zz_x_xxy[i] = t_0_zz_0_xxxy[i] - rcd_x[i] * t_0_zz_0_xxy[i];

        t_0_zz_x_xxx[i] = t_0_zz_0_xxxx[i] - rcd_x[i] * t_0_zz_0_xxx[i];

        t_0_yz_z_zzz[i] = t_0_yz_0_zzzz[i] - rcd_z[i] * t_0_yz_0_zzz[i];

        t_0_yz_z_yzz[i] = t_0_yz_0_yzzz[i] - rcd_z[i] * t_0_yz_0_yzz[i];

        t_0_yz_z_yyz[i] = t_0_yz_0_yyzz[i] - rcd_z[i] * t_0_yz_0_yyz[i];

        t_0_yz_z_xzz[i] = t_0_yz_0_xzzz[i] - rcd_z[i] * t_0_yz_0_xzz[i];

        t_0_yz_z_xyz[i] = t_0_yz_0_xyzz[i] - rcd_z[i] * t_0_yz_0_xyz[i];

        t_0_yz_z_xxz[i] = t_0_yz_0_xxzz[i] - rcd_z[i] * t_0_yz_0_xxz[i];

        t_0_yz_y_zzz[i] = t_0_yz_0_yzzz[i] - rcd_y[i] * t_0_yz_0_zzz[i];

        t_0_yz_y_yzz[i] = t_0_yz_0_yyzz[i] - rcd_y[i] * t_0_yz_0_yzz[i];

        t_0_yz_y_yyz[i] = t_0_yz_0_yyyz[i] - rcd_y[i] * t_0_yz_0_yyz[i];

        t_0_yz_y_yyy[i] = t_0_yz_0_yyyy[i] - rcd_y[i] * t_0_yz_0_yyy[i];

        t_0_yz_y_xzz[i] = t_0_yz_0_xyzz[i] - rcd_y[i] * t_0_yz_0_xzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yy_0_xxy, t_0_yy_0_xxyy, t_0_yy_0_xxyz,\
                           t_0_yy_0_xxz, t_0_yy_0_xxzz, t_0_yy_0_xyy, t_0_yy_0_xyyy,\
                           t_0_yy_0_xyyz, t_0_yy_0_xyz, t_0_yy_0_xyzz, t_0_yy_0_xzz,\
                           t_0_yy_0_xzzz, t_0_yy_0_yyy, t_0_yy_0_yyyy, t_0_yy_0_yyyz,\
                           t_0_yy_0_yyz, t_0_yy_0_yyzz, t_0_yy_0_yzz, t_0_yy_0_yzzz,\
                           t_0_yy_0_zzz, t_0_yy_0_zzzz, t_0_yy_x_xyy, t_0_yy_x_xyz,\
                           t_0_yy_x_xzz, t_0_yy_x_yyy, t_0_yy_x_yyz, t_0_yy_x_yzz, t_0_yy_x_zzz,\
                           t_0_yy_y_xxy, t_0_yy_y_xxz, t_0_yy_y_xyy, t_0_yy_y_xyz, t_0_yy_y_xzz,\
                           t_0_yy_y_yyy, t_0_yy_y_yyz, t_0_yy_y_yzz, t_0_yy_y_zzz, t_0_yy_z_xxz,\
                           t_0_yy_z_xyz, t_0_yy_z_xzz, t_0_yy_z_yyz, t_0_yy_z_yzz, t_0_yy_z_zzz,\
                           t_0_yz_0_xxx, t_0_yz_0_xxxx, t_0_yz_0_xxxy, t_0_yz_0_xxxz,\
                           t_0_yz_0_xxy, t_0_yz_0_xxyy, t_0_yz_0_xxyz, t_0_yz_0_xxz,\
                           t_0_yz_0_xxzz, t_0_yz_0_xyy, t_0_yz_0_xyyy, t_0_yz_0_xyyz,\
                           t_0_yz_0_xyz, t_0_yz_0_xyzz, t_0_yz_0_xzz, t_0_yz_0_xzzz,\
                           t_0_yz_0_yyy, t_0_yz_0_yyz, t_0_yz_0_yzz, t_0_yz_0_zzz, t_0_yz_x_xxx,\
                           t_0_yz_x_xxy, t_0_yz_x_xxz, t_0_yz_x_xyy, t_0_yz_x_xyz, t_0_yz_x_xzz,\
                           t_0_yz_x_yyy, t_0_yz_x_yyz, t_0_yz_x_yzz, t_0_yz_x_zzz, t_0_yz_y_xxy,\
                           t_0_yz_y_xxz, t_0_yz_y_xyy, t_0_yz_y_xyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yz_y_xyz[i] = t_0_yz_0_xyyz[i] - rcd_y[i] * t_0_yz_0_xyz[i];

        t_0_yz_y_xyy[i] = t_0_yz_0_xyyy[i] - rcd_y[i] * t_0_yz_0_xyy[i];

        t_0_yz_y_xxz[i] = t_0_yz_0_xxyz[i] - rcd_y[i] * t_0_yz_0_xxz[i];

        t_0_yz_y_xxy[i] = t_0_yz_0_xxyy[i] - rcd_y[i] * t_0_yz_0_xxy[i];

        t_0_yz_x_zzz[i] = t_0_yz_0_xzzz[i] - rcd_x[i] * t_0_yz_0_zzz[i];

        t_0_yz_x_yzz[i] = t_0_yz_0_xyzz[i] - rcd_x[i] * t_0_yz_0_yzz[i];

        t_0_yz_x_yyz[i] = t_0_yz_0_xyyz[i] - rcd_x[i] * t_0_yz_0_yyz[i];

        t_0_yz_x_yyy[i] = t_0_yz_0_xyyy[i] - rcd_x[i] * t_0_yz_0_yyy[i];

        t_0_yz_x_xzz[i] = t_0_yz_0_xxzz[i] - rcd_x[i] * t_0_yz_0_xzz[i];

        t_0_yz_x_xyz[i] = t_0_yz_0_xxyz[i] - rcd_x[i] * t_0_yz_0_xyz[i];

        t_0_yz_x_xyy[i] = t_0_yz_0_xxyy[i] - rcd_x[i] * t_0_yz_0_xyy[i];

        t_0_yz_x_xxz[i] = t_0_yz_0_xxxz[i] - rcd_x[i] * t_0_yz_0_xxz[i];

        t_0_yz_x_xxy[i] = t_0_yz_0_xxxy[i] - rcd_x[i] * t_0_yz_0_xxy[i];

        t_0_yz_x_xxx[i] = t_0_yz_0_xxxx[i] - rcd_x[i] * t_0_yz_0_xxx[i];

        t_0_yy_z_zzz[i] = t_0_yy_0_zzzz[i] - rcd_z[i] * t_0_yy_0_zzz[i];

        t_0_yy_z_yzz[i] = t_0_yy_0_yzzz[i] - rcd_z[i] * t_0_yy_0_yzz[i];

        t_0_yy_z_yyz[i] = t_0_yy_0_yyzz[i] - rcd_z[i] * t_0_yy_0_yyz[i];

        t_0_yy_z_xzz[i] = t_0_yy_0_xzzz[i] - rcd_z[i] * t_0_yy_0_xzz[i];

        t_0_yy_z_xyz[i] = t_0_yy_0_xyzz[i] - rcd_z[i] * t_0_yy_0_xyz[i];

        t_0_yy_z_xxz[i] = t_0_yy_0_xxzz[i] - rcd_z[i] * t_0_yy_0_xxz[i];

        t_0_yy_y_zzz[i] = t_0_yy_0_yzzz[i] - rcd_y[i] * t_0_yy_0_zzz[i];

        t_0_yy_y_yzz[i] = t_0_yy_0_yyzz[i] - rcd_y[i] * t_0_yy_0_yzz[i];

        t_0_yy_y_yyz[i] = t_0_yy_0_yyyz[i] - rcd_y[i] * t_0_yy_0_yyz[i];

        t_0_yy_y_yyy[i] = t_0_yy_0_yyyy[i] - rcd_y[i] * t_0_yy_0_yyy[i];

        t_0_yy_y_xzz[i] = t_0_yy_0_xyzz[i] - rcd_y[i] * t_0_yy_0_xzz[i];

        t_0_yy_y_xyz[i] = t_0_yy_0_xyyz[i] - rcd_y[i] * t_0_yy_0_xyz[i];

        t_0_yy_y_xyy[i] = t_0_yy_0_xyyy[i] - rcd_y[i] * t_0_yy_0_xyy[i];

        t_0_yy_y_xxz[i] = t_0_yy_0_xxyz[i] - rcd_y[i] * t_0_yy_0_xxz[i];

        t_0_yy_y_xxy[i] = t_0_yy_0_xxyy[i] - rcd_y[i] * t_0_yy_0_xxy[i];

        t_0_yy_x_zzz[i] = t_0_yy_0_xzzz[i] - rcd_x[i] * t_0_yy_0_zzz[i];

        t_0_yy_x_yzz[i] = t_0_yy_0_xyzz[i] - rcd_x[i] * t_0_yy_0_yzz[i];

        t_0_yy_x_yyz[i] = t_0_yy_0_xyyz[i] - rcd_x[i] * t_0_yy_0_yyz[i];

        t_0_yy_x_yyy[i] = t_0_yy_0_xyyy[i] - rcd_x[i] * t_0_yy_0_yyy[i];

        t_0_yy_x_xzz[i] = t_0_yy_0_xxzz[i] - rcd_x[i] * t_0_yy_0_xzz[i];

        t_0_yy_x_xyz[i] = t_0_yy_0_xxyz[i] - rcd_x[i] * t_0_yy_0_xyz[i];

        t_0_yy_x_xyy[i] = t_0_yy_0_xxyy[i] - rcd_x[i] * t_0_yy_0_xyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xy_0_xxz, t_0_xy_0_xxzz, t_0_xy_0_xyz,\
                           t_0_xy_0_xyzz, t_0_xy_0_xzz, t_0_xy_0_xzzz, t_0_xy_0_yyz,\
                           t_0_xy_0_yyzz, t_0_xy_0_yzz, t_0_xy_0_yzzz, t_0_xy_0_zzz,\
                           t_0_xy_0_zzzz, t_0_xy_y_yzz, t_0_xy_y_zzz, t_0_xy_z_xxz,\
                           t_0_xy_z_xyz, t_0_xy_z_xzz, t_0_xy_z_yyz, t_0_xy_z_yzz, t_0_xy_z_zzz,\
                           t_0_xz_0_xxx, t_0_xz_0_xxxx, t_0_xz_0_xxxy, t_0_xz_0_xxxz,\
                           t_0_xz_0_xxy, t_0_xz_0_xxyy, t_0_xz_0_xxyz, t_0_xz_0_xxz,\
                           t_0_xz_0_xxzz, t_0_xz_0_xyy, t_0_xz_0_xyyy, t_0_xz_0_xyyz,\
                           t_0_xz_0_xyz, t_0_xz_0_xyzz, t_0_xz_0_xzz, t_0_xz_0_xzzz,\
                           t_0_xz_0_yyy, t_0_xz_0_yyyy, t_0_xz_0_yyyz, t_0_xz_0_yyz,\
                           t_0_xz_0_yyzz, t_0_xz_0_yzz, t_0_xz_0_yzzz, t_0_xz_0_zzz,\
                           t_0_xz_0_zzzz, t_0_xz_x_xxx, t_0_xz_x_xxy, t_0_xz_x_xxz,\
                           t_0_xz_x_xyy, t_0_xz_x_xyz, t_0_xz_x_xzz, t_0_xz_x_yyy, t_0_xz_x_yyz,\
                           t_0_xz_x_yzz, t_0_xz_x_zzz, t_0_xz_y_xxy, t_0_xz_y_xxz, t_0_xz_y_xyy,\
                           t_0_xz_y_xyz, t_0_xz_y_xzz, t_0_xz_y_yyy, t_0_xz_y_yyz, t_0_xz_y_yzz,\
                           t_0_xz_y_zzz, t_0_xz_z_xxz, t_0_xz_z_xyz, t_0_xz_z_xzz, t_0_xz_z_yyz,\
                           t_0_xz_z_yzz, t_0_xz_z_zzz, t_0_yy_0_xxx, t_0_yy_0_xxxx,\
                           t_0_yy_0_xxxy, t_0_yy_0_xxxz, t_0_yy_0_xxy, t_0_yy_0_xxz,\
                           t_0_yy_x_xxx, t_0_yy_x_xxy, t_0_yy_x_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yy_x_xxz[i] = t_0_yy_0_xxxz[i] - rcd_x[i] * t_0_yy_0_xxz[i];

        t_0_yy_x_xxy[i] = t_0_yy_0_xxxy[i] - rcd_x[i] * t_0_yy_0_xxy[i];

        t_0_yy_x_xxx[i] = t_0_yy_0_xxxx[i] - rcd_x[i] * t_0_yy_0_xxx[i];

        t_0_xz_z_zzz[i] = t_0_xz_0_zzzz[i] - rcd_z[i] * t_0_xz_0_zzz[i];

        t_0_xz_z_yzz[i] = t_0_xz_0_yzzz[i] - rcd_z[i] * t_0_xz_0_yzz[i];

        t_0_xz_z_yyz[i] = t_0_xz_0_yyzz[i] - rcd_z[i] * t_0_xz_0_yyz[i];

        t_0_xz_z_xzz[i] = t_0_xz_0_xzzz[i] - rcd_z[i] * t_0_xz_0_xzz[i];

        t_0_xz_z_xyz[i] = t_0_xz_0_xyzz[i] - rcd_z[i] * t_0_xz_0_xyz[i];

        t_0_xz_z_xxz[i] = t_0_xz_0_xxzz[i] - rcd_z[i] * t_0_xz_0_xxz[i];

        t_0_xz_y_zzz[i] = t_0_xz_0_yzzz[i] - rcd_y[i] * t_0_xz_0_zzz[i];

        t_0_xz_y_yzz[i] = t_0_xz_0_yyzz[i] - rcd_y[i] * t_0_xz_0_yzz[i];

        t_0_xz_y_yyz[i] = t_0_xz_0_yyyz[i] - rcd_y[i] * t_0_xz_0_yyz[i];

        t_0_xz_y_yyy[i] = t_0_xz_0_yyyy[i] - rcd_y[i] * t_0_xz_0_yyy[i];

        t_0_xz_y_xzz[i] = t_0_xz_0_xyzz[i] - rcd_y[i] * t_0_xz_0_xzz[i];

        t_0_xz_y_xyz[i] = t_0_xz_0_xyyz[i] - rcd_y[i] * t_0_xz_0_xyz[i];

        t_0_xz_y_xyy[i] = t_0_xz_0_xyyy[i] - rcd_y[i] * t_0_xz_0_xyy[i];

        t_0_xz_y_xxz[i] = t_0_xz_0_xxyz[i] - rcd_y[i] * t_0_xz_0_xxz[i];

        t_0_xz_y_xxy[i] = t_0_xz_0_xxyy[i] - rcd_y[i] * t_0_xz_0_xxy[i];

        t_0_xz_x_zzz[i] = t_0_xz_0_xzzz[i] - rcd_x[i] * t_0_xz_0_zzz[i];

        t_0_xz_x_yzz[i] = t_0_xz_0_xyzz[i] - rcd_x[i] * t_0_xz_0_yzz[i];

        t_0_xz_x_yyz[i] = t_0_xz_0_xyyz[i] - rcd_x[i] * t_0_xz_0_yyz[i];

        t_0_xz_x_yyy[i] = t_0_xz_0_xyyy[i] - rcd_x[i] * t_0_xz_0_yyy[i];

        t_0_xz_x_xzz[i] = t_0_xz_0_xxzz[i] - rcd_x[i] * t_0_xz_0_xzz[i];

        t_0_xz_x_xyz[i] = t_0_xz_0_xxyz[i] - rcd_x[i] * t_0_xz_0_xyz[i];

        t_0_xz_x_xyy[i] = t_0_xz_0_xxyy[i] - rcd_x[i] * t_0_xz_0_xyy[i];

        t_0_xz_x_xxz[i] = t_0_xz_0_xxxz[i] - rcd_x[i] * t_0_xz_0_xxz[i];

        t_0_xz_x_xxy[i] = t_0_xz_0_xxxy[i] - rcd_x[i] * t_0_xz_0_xxy[i];

        t_0_xz_x_xxx[i] = t_0_xz_0_xxxx[i] - rcd_x[i] * t_0_xz_0_xxx[i];

        t_0_xy_z_zzz[i] = t_0_xy_0_zzzz[i] - rcd_z[i] * t_0_xy_0_zzz[i];

        t_0_xy_z_yzz[i] = t_0_xy_0_yzzz[i] - rcd_z[i] * t_0_xy_0_yzz[i];

        t_0_xy_z_yyz[i] = t_0_xy_0_yyzz[i] - rcd_z[i] * t_0_xy_0_yyz[i];

        t_0_xy_z_xzz[i] = t_0_xy_0_xzzz[i] - rcd_z[i] * t_0_xy_0_xzz[i];

        t_0_xy_z_xyz[i] = t_0_xy_0_xyzz[i] - rcd_z[i] * t_0_xy_0_xyz[i];

        t_0_xy_z_xxz[i] = t_0_xy_0_xxzz[i] - rcd_z[i] * t_0_xy_0_xxz[i];

        t_0_xy_y_zzz[i] = t_0_xy_0_yzzz[i] - rcd_y[i] * t_0_xy_0_zzz[i];

        t_0_xy_y_yzz[i] = t_0_xy_0_yyzz[i] - rcd_y[i] * t_0_xy_0_yzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_0_xxy, t_0_xx_0_xxyy, t_0_xx_0_xxyz,\
                           t_0_xx_0_xxz, t_0_xx_0_xxzz, t_0_xx_0_xyy, t_0_xx_0_xyyy,\
                           t_0_xx_0_xyyz, t_0_xx_0_xyz, t_0_xx_0_xyzz, t_0_xx_0_xzz,\
                           t_0_xx_0_xzzz, t_0_xx_0_yyy, t_0_xx_0_yyyy, t_0_xx_0_yyyz,\
                           t_0_xx_0_yyz, t_0_xx_0_yyzz, t_0_xx_0_yzz, t_0_xx_0_yzzz,\
                           t_0_xx_0_zzz, t_0_xx_0_zzzz, t_0_xx_x_yyy, t_0_xx_x_yyz,\
                           t_0_xx_x_yzz, t_0_xx_x_zzz, t_0_xx_y_xxy, t_0_xx_y_xxz, t_0_xx_y_xyy,\
                           t_0_xx_y_xyz, t_0_xx_y_xzz, t_0_xx_y_yyy, t_0_xx_y_yyz, t_0_xx_y_yzz,\
                           t_0_xx_y_zzz, t_0_xx_z_xxz, t_0_xx_z_xyz, t_0_xx_z_xzz, t_0_xx_z_yyz,\
                           t_0_xx_z_yzz, t_0_xx_z_zzz, t_0_xy_0_xxx, t_0_xy_0_xxxx,\
                           t_0_xy_0_xxxy, t_0_xy_0_xxxz, t_0_xy_0_xxy, t_0_xy_0_xxyy,\
                           t_0_xy_0_xxyz, t_0_xy_0_xxz, t_0_xy_0_xxzz, t_0_xy_0_xyy,\
                           t_0_xy_0_xyyy, t_0_xy_0_xyyz, t_0_xy_0_xyz, t_0_xy_0_xyzz,\
                           t_0_xy_0_xzz, t_0_xy_0_xzzz, t_0_xy_0_yyy, t_0_xy_0_yyyy,\
                           t_0_xy_0_yyyz, t_0_xy_0_yyz, t_0_xy_0_yzz, t_0_xy_0_zzz,\
                           t_0_xy_x_xxx, t_0_xy_x_xxy, t_0_xy_x_xxz, t_0_xy_x_xyy, t_0_xy_x_xyz,\
                           t_0_xy_x_xzz, t_0_xy_x_yyy, t_0_xy_x_yyz, t_0_xy_x_yzz, t_0_xy_x_zzz,\
                           t_0_xy_y_xxy, t_0_xy_y_xxz, t_0_xy_y_xyy, t_0_xy_y_xyz, t_0_xy_y_xzz,\
                           t_0_xy_y_yyy, t_0_xy_y_yyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xy_y_yyz[i] = t_0_xy_0_yyyz[i] - rcd_y[i] * t_0_xy_0_yyz[i];

        t_0_xy_y_yyy[i] = t_0_xy_0_yyyy[i] - rcd_y[i] * t_0_xy_0_yyy[i];

        t_0_xy_y_xzz[i] = t_0_xy_0_xyzz[i] - rcd_y[i] * t_0_xy_0_xzz[i];

        t_0_xy_y_xyz[i] = t_0_xy_0_xyyz[i] - rcd_y[i] * t_0_xy_0_xyz[i];

        t_0_xy_y_xyy[i] = t_0_xy_0_xyyy[i] - rcd_y[i] * t_0_xy_0_xyy[i];

        t_0_xy_y_xxz[i] = t_0_xy_0_xxyz[i] - rcd_y[i] * t_0_xy_0_xxz[i];

        t_0_xy_y_xxy[i] = t_0_xy_0_xxyy[i] - rcd_y[i] * t_0_xy_0_xxy[i];

        t_0_xy_x_zzz[i] = t_0_xy_0_xzzz[i] - rcd_x[i] * t_0_xy_0_zzz[i];

        t_0_xy_x_yzz[i] = t_0_xy_0_xyzz[i] - rcd_x[i] * t_0_xy_0_yzz[i];

        t_0_xy_x_yyz[i] = t_0_xy_0_xyyz[i] - rcd_x[i] * t_0_xy_0_yyz[i];

        t_0_xy_x_yyy[i] = t_0_xy_0_xyyy[i] - rcd_x[i] * t_0_xy_0_yyy[i];

        t_0_xy_x_xzz[i] = t_0_xy_0_xxzz[i] - rcd_x[i] * t_0_xy_0_xzz[i];

        t_0_xy_x_xyz[i] = t_0_xy_0_xxyz[i] - rcd_x[i] * t_0_xy_0_xyz[i];

        t_0_xy_x_xyy[i] = t_0_xy_0_xxyy[i] - rcd_x[i] * t_0_xy_0_xyy[i];

        t_0_xy_x_xxz[i] = t_0_xy_0_xxxz[i] - rcd_x[i] * t_0_xy_0_xxz[i];

        t_0_xy_x_xxy[i] = t_0_xy_0_xxxy[i] - rcd_x[i] * t_0_xy_0_xxy[i];

        t_0_xy_x_xxx[i] = t_0_xy_0_xxxx[i] - rcd_x[i] * t_0_xy_0_xxx[i];

        t_0_xx_z_zzz[i] = t_0_xx_0_zzzz[i] - rcd_z[i] * t_0_xx_0_zzz[i];

        t_0_xx_z_yzz[i] = t_0_xx_0_yzzz[i] - rcd_z[i] * t_0_xx_0_yzz[i];

        t_0_xx_z_yyz[i] = t_0_xx_0_yyzz[i] - rcd_z[i] * t_0_xx_0_yyz[i];

        t_0_xx_z_xzz[i] = t_0_xx_0_xzzz[i] - rcd_z[i] * t_0_xx_0_xzz[i];

        t_0_xx_z_xyz[i] = t_0_xx_0_xyzz[i] - rcd_z[i] * t_0_xx_0_xyz[i];

        t_0_xx_z_xxz[i] = t_0_xx_0_xxzz[i] - rcd_z[i] * t_0_xx_0_xxz[i];

        t_0_xx_y_zzz[i] = t_0_xx_0_yzzz[i] - rcd_y[i] * t_0_xx_0_zzz[i];

        t_0_xx_y_yzz[i] = t_0_xx_0_yyzz[i] - rcd_y[i] * t_0_xx_0_yzz[i];

        t_0_xx_y_yyz[i] = t_0_xx_0_yyyz[i] - rcd_y[i] * t_0_xx_0_yyz[i];

        t_0_xx_y_yyy[i] = t_0_xx_0_yyyy[i] - rcd_y[i] * t_0_xx_0_yyy[i];

        t_0_xx_y_xzz[i] = t_0_xx_0_xyzz[i] - rcd_y[i] * t_0_xx_0_xzz[i];

        t_0_xx_y_xyz[i] = t_0_xx_0_xyyz[i] - rcd_y[i] * t_0_xx_0_xyz[i];

        t_0_xx_y_xyy[i] = t_0_xx_0_xyyy[i] - rcd_y[i] * t_0_xx_0_xyy[i];

        t_0_xx_y_xxz[i] = t_0_xx_0_xxyz[i] - rcd_y[i] * t_0_xx_0_xxz[i];

        t_0_xx_y_xxy[i] = t_0_xx_0_xxyy[i] - rcd_y[i] * t_0_xx_0_xxy[i];

        t_0_xx_x_zzz[i] = t_0_xx_0_xzzz[i] - rcd_x[i] * t_0_xx_0_zzz[i];

        t_0_xx_x_yzz[i] = t_0_xx_0_xyzz[i] - rcd_x[i] * t_0_xx_0_yzz[i];

        t_0_xx_x_yyz[i] = t_0_xx_0_xyyz[i] - rcd_x[i] * t_0_xx_0_yyz[i];

        t_0_xx_x_yyy[i] = t_0_xx_0_xyyy[i] - rcd_x[i] * t_0_xx_0_yyy[i];
    }

    #pragma omp simd align(rcd_x, t_0_xx_0_xxx, t_0_xx_0_xxxx, t_0_xx_0_xxxy, t_0_xx_0_xxxz,\
                           t_0_xx_0_xxy, t_0_xx_0_xxyy, t_0_xx_0_xxyz, t_0_xx_0_xxz,\
                           t_0_xx_0_xxzz, t_0_xx_0_xyy, t_0_xx_0_xyz, t_0_xx_0_xzz,\
                           t_0_xx_x_xxx, t_0_xx_x_xxy, t_0_xx_x_xxz, t_0_xx_x_xyy, t_0_xx_x_xyz,\
                           t_0_xx_x_xzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xx_x_xzz[i] = t_0_xx_0_xxzz[i] - rcd_x[i] * t_0_xx_0_xzz[i];

        t_0_xx_x_xyz[i] = t_0_xx_0_xxyz[i] - rcd_x[i] * t_0_xx_0_xyz[i];

        t_0_xx_x_xyy[i] = t_0_xx_0_xxyy[i] - rcd_x[i] * t_0_xx_0_xyy[i];

        t_0_xx_x_xxz[i] = t_0_xx_0_xxxz[i] - rcd_x[i] * t_0_xx_0_xxz[i];

        t_0_xx_x_xxy[i] = t_0_xx_0_xxxy[i] - rcd_x[i] * t_0_xx_0_xxy[i];

        t_0_xx_x_xxx[i] = t_0_xx_0_xxxx[i] - rcd_x[i] * t_0_xx_0_xxx[i];
    }
}


} // derirec namespace
