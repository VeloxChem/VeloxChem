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
compHostHRRForSFPD_V0(      BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPD) integral components

    t_0_zzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(0));

    t_0_zzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(1));

    t_0_zzz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(2));

    t_0_zzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(3));

    t_0_zzz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(4));

    t_0_zzz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(5));

    t_0_zzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(6));

    t_0_zzz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(7));

    t_0_zzz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(8));

    t_0_zzz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(9));

    t_0_zzz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(10));

    t_0_zzz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(11));

    t_0_zzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(12));

    t_0_zzz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(13));

    t_0_zzz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(14));

    t_0_zzz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(15));

    t_0_zzz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(16));

    t_0_zzz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(17));

    t_0_yzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(18));

    t_0_yzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(19));

    t_0_yzz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(20));

    t_0_yzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(21));

    t_0_yzz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(22));

    t_0_yzz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(23));

    t_0_yzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(24));

    t_0_yzz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(25));

    t_0_yzz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(26));

    t_0_yzz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(27));

    t_0_yzz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(28));

    t_0_yzz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(29));

    t_0_yzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(30));

    t_0_yzz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(31));

    t_0_yzz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(32));

    t_0_yzz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(33));

    t_0_yzz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(34));

    t_0_yzz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(35));

    t_0_yyz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(36));

    t_0_yyz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(37));

    t_0_yyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(38));

    t_0_yyz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(39));

    t_0_yyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(40));

    t_0_yyz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(41));

    t_0_yyz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(42));

    t_0_yyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(43));

    t_0_yyz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(44));

    t_0_yyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(45));

    t_0_yyz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(46));

    t_0_yyz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(47));

    t_0_yyz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(48));

    t_0_yyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(49));

    t_0_yyz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(50));

    t_0_yyz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(51));

    t_0_yyz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(52));

    t_0_yyz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(53));

    t_0_yyy_z_zz = intsBufferSFPD.data(intsIndexesSFPD(54));

    t_0_yyy_z_yz = intsBufferSFPD.data(intsIndexesSFPD(55));

    t_0_yyy_z_yy = intsBufferSFPD.data(intsIndexesSFPD(56));

    t_0_yyy_z_xz = intsBufferSFPD.data(intsIndexesSFPD(57));

    t_0_yyy_z_xy = intsBufferSFPD.data(intsIndexesSFPD(58));

    t_0_yyy_z_xx = intsBufferSFPD.data(intsIndexesSFPD(59));

    t_0_yyy_y_zz = intsBufferSFPD.data(intsIndexesSFPD(60));

    t_0_yyy_y_yz = intsBufferSFPD.data(intsIndexesSFPD(61));

    t_0_yyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(62));

    t_0_yyy_y_xz = intsBufferSFPD.data(intsIndexesSFPD(63));

    t_0_yyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(64));

    t_0_yyy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(65));

    t_0_yyy_x_zz = intsBufferSFPD.data(intsIndexesSFPD(66));

    t_0_yyy_x_yz = intsBufferSFPD.data(intsIndexesSFPD(67));

    t_0_yyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(68));

    t_0_yyy_x_xz = intsBufferSFPD.data(intsIndexesSFPD(69));

    t_0_yyy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(70));

    t_0_yyy_x_xx = intsBufferSFPD.data(intsIndexesSFPD(71));

    t_0_xzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(72));

    t_0_xzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(73));

    t_0_xzz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(74));

    t_0_xzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(75));

    t_0_xzz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(76));

    t_0_xzz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(77));

    t_0_xzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(78));

    t_0_xzz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(79));

    t_0_xzz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(80));

    t_0_xzz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(81));

    t_0_xzz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(82));

    t_0_xzz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(83));

    t_0_xzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(84));

    t_0_xzz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(85));

    t_0_xzz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(86));

    t_0_xzz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(87));

    t_0_xzz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(88));

    t_0_xzz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(89));

    t_0_xyz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(90));

    t_0_xyz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(91));

    t_0_xyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(92));

    t_0_xyz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(93));

    t_0_xyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(94));

    t_0_xyz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(95));

    t_0_xyz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(96));

    t_0_xyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(97));

    t_0_xyz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(98));

    t_0_xyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(99));

    t_0_xyz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(100));

    t_0_xyz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(101));

    t_0_xyz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(102));

    t_0_xyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(103));

    t_0_xyz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(104));

    t_0_xyz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(105));

    t_0_xyz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(106));

    t_0_xyz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(107));

    t_0_xyy_z_zz = intsBufferSFPD.data(intsIndexesSFPD(108));

    t_0_xyy_z_yz = intsBufferSFPD.data(intsIndexesSFPD(109));

    t_0_xyy_z_yy = intsBufferSFPD.data(intsIndexesSFPD(110));

    t_0_xyy_z_xz = intsBufferSFPD.data(intsIndexesSFPD(111));

    t_0_xyy_z_xy = intsBufferSFPD.data(intsIndexesSFPD(112));

    t_0_xyy_z_xx = intsBufferSFPD.data(intsIndexesSFPD(113));

    t_0_xyy_y_zz = intsBufferSFPD.data(intsIndexesSFPD(114));

    t_0_xyy_y_yz = intsBufferSFPD.data(intsIndexesSFPD(115));

    t_0_xyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(116));

    t_0_xyy_y_xz = intsBufferSFPD.data(intsIndexesSFPD(117));

    t_0_xyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(118));

    t_0_xyy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(119));

    t_0_xyy_x_zz = intsBufferSFPD.data(intsIndexesSFPD(120));

    t_0_xyy_x_yz = intsBufferSFPD.data(intsIndexesSFPD(121));

    t_0_xyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(122));

    t_0_xyy_x_xz = intsBufferSFPD.data(intsIndexesSFPD(123));

    t_0_xyy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(124));

    t_0_xyy_x_xx = intsBufferSFPD.data(intsIndexesSFPD(125));

    t_0_xxz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(126));

    t_0_xxz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(127));

    t_0_xxz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(128));

    t_0_xxz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(129));

    t_0_xxz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(130));

    t_0_xxz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(131));

    t_0_xxz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(132));

    t_0_xxz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(133));

    t_0_xxz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(134));

    t_0_xxz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(135));

    t_0_xxz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(136));

    t_0_xxz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(137));

    t_0_xxz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(138));

    t_0_xxz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(139));

    t_0_xxz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(140));

    t_0_xxz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(141));

    t_0_xxz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(142));

    t_0_xxz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(143));

    t_0_xxy_z_zz = intsBufferSFPD.data(intsIndexesSFPD(144));

    t_0_xxy_z_yz = intsBufferSFPD.data(intsIndexesSFPD(145));

    t_0_xxy_z_yy = intsBufferSFPD.data(intsIndexesSFPD(146));

    t_0_xxy_z_xz = intsBufferSFPD.data(intsIndexesSFPD(147));

    t_0_xxy_z_xy = intsBufferSFPD.data(intsIndexesSFPD(148));

    t_0_xxy_z_xx = intsBufferSFPD.data(intsIndexesSFPD(149));

    t_0_xxy_y_zz = intsBufferSFPD.data(intsIndexesSFPD(150));

    t_0_xxy_y_yz = intsBufferSFPD.data(intsIndexesSFPD(151));

    t_0_xxy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(152));

    t_0_xxy_y_xz = intsBufferSFPD.data(intsIndexesSFPD(153));

    t_0_xxy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(154));

    t_0_xxy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(155));

    t_0_xxy_x_zz = intsBufferSFPD.data(intsIndexesSFPD(156));

    t_0_xxy_x_yz = intsBufferSFPD.data(intsIndexesSFPD(157));

    t_0_xxy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(158));

    t_0_xxy_x_xz = intsBufferSFPD.data(intsIndexesSFPD(159));

    t_0_xxy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(160));

    t_0_xxy_x_xx = intsBufferSFPD.data(intsIndexesSFPD(161));

    t_0_xxx_z_zz = intsBufferSFPD.data(intsIndexesSFPD(162));

    t_0_xxx_z_yz = intsBufferSFPD.data(intsIndexesSFPD(163));

    t_0_xxx_z_yy = intsBufferSFPD.data(intsIndexesSFPD(164));

    t_0_xxx_z_xz = intsBufferSFPD.data(intsIndexesSFPD(165));

    t_0_xxx_z_xy = intsBufferSFPD.data(intsIndexesSFPD(166));

    t_0_xxx_z_xx = intsBufferSFPD.data(intsIndexesSFPD(167));

    t_0_xxx_y_zz = intsBufferSFPD.data(intsIndexesSFPD(168));

    t_0_xxx_y_yz = intsBufferSFPD.data(intsIndexesSFPD(169));

    t_0_xxx_y_yy = intsBufferSFPD.data(intsIndexesSFPD(170));

    t_0_xxx_y_xz = intsBufferSFPD.data(intsIndexesSFPD(171));

    t_0_xxx_y_xy = intsBufferSFPD.data(intsIndexesSFPD(172));

    t_0_xxx_y_xx = intsBufferSFPD.data(intsIndexesSFPD(173));

    t_0_xxx_x_zz = intsBufferSFPD.data(intsIndexesSFPD(174));

    t_0_xxx_x_yz = intsBufferSFPD.data(intsIndexesSFPD(175));

    t_0_xxx_x_yy = intsBufferSFPD.data(intsIndexesSFPD(176));

    t_0_xxx_x_xz = intsBufferSFPD.data(intsIndexesSFPD(177));

    t_0_xxx_x_xy = intsBufferSFPD.data(intsIndexesSFPD(178));

    t_0_xxx_x_xx = intsBufferSFPD.data(intsIndexesSFPD(179));

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

    // set up (SFSF) integral components

    t_0_zzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(0));

    t_0_zzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(1));

    t_0_zzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(2));

    t_0_zzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(3));

    t_0_zzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(4));

    t_0_zzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(5));

    t_0_zzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(6));

    t_0_zzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(7));

    t_0_zzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(8));

    t_0_zzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(9));

    t_0_yzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(10));

    t_0_yzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(11));

    t_0_yzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(12));

    t_0_yzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(13));

    t_0_yzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(14));

    t_0_yzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(15));

    t_0_yzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(16));

    t_0_yzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(17));

    t_0_yzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(18));

    t_0_yzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(19));

    t_0_yyz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(20));

    t_0_yyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(21));

    t_0_yyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(22));

    t_0_yyz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(23));

    t_0_yyz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(24));

    t_0_yyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(25));

    t_0_yyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(26));

    t_0_yyz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(27));

    t_0_yyz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(28));

    t_0_yyz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(29));

    t_0_yyy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(30));

    t_0_yyy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(31));

    t_0_yyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(32));

    t_0_yyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(33));

    t_0_yyy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(34));

    t_0_yyy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(35));

    t_0_yyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(36));

    t_0_yyy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(37));

    t_0_yyy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(38));

    t_0_yyy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(39));

    t_0_xzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(40));

    t_0_xzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(41));

    t_0_xzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(42));

    t_0_xzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(43));

    t_0_xzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(44));

    t_0_xzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(45));

    t_0_xzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(46));

    t_0_xzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(47));

    t_0_xzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(48));

    t_0_xzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(49));

    t_0_xyz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(50));

    t_0_xyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(51));

    t_0_xyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(52));

    t_0_xyz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(53));

    t_0_xyz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(54));

    t_0_xyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(55));

    t_0_xyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(56));

    t_0_xyz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(57));

    t_0_xyz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(58));

    t_0_xyz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(59));

    t_0_xyy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(60));

    t_0_xyy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(61));

    t_0_xyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(62));

    t_0_xyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(63));

    t_0_xyy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(64));

    t_0_xyy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(65));

    t_0_xyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(66));

    t_0_xyy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(67));

    t_0_xyy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(68));

    t_0_xyy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(69));

    t_0_xxz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(70));

    t_0_xxz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(71));

    t_0_xxz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(72));

    t_0_xxz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(73));

    t_0_xxz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(74));

    t_0_xxz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(75));

    t_0_xxz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(76));

    t_0_xxz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(77));

    t_0_xxz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(78));

    t_0_xxz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(79));

    t_0_xxy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(80));

    t_0_xxy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(81));

    t_0_xxy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(82));

    t_0_xxy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(83));

    t_0_xxy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(84));

    t_0_xxy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(85));

    t_0_xxy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(86));

    t_0_xxy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(87));

    t_0_xxy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(88));

    t_0_xxy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(89));

    t_0_xxx_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(90));

    t_0_xxx_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(91));

    t_0_xxx_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(92));

    t_0_xxx_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(93));

    t_0_xxx_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(94));

    t_0_xxx_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(95));

    t_0_xxx_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(96));

    t_0_xxx_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(97));

    t_0_xxx_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(98));

    t_0_xxx_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(99));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzz_0_xx, t_0_yzz_0_xxx, t_0_yzz_0_xxy,\
                           t_0_yzz_0_xxz, t_0_yzz_0_xy, t_0_yzz_0_xyy, t_0_yzz_0_xyz,\
                           t_0_yzz_0_xz, t_0_yzz_0_xzz, t_0_yzz_0_yy, t_0_yzz_0_yyy,\
                           t_0_yzz_0_yyz, t_0_yzz_0_yz, t_0_yzz_0_yzz, t_0_yzz_0_zz,\
                           t_0_yzz_0_zzz, t_0_yzz_x_xx, t_0_yzz_x_xy, t_0_yzz_x_xz,\
                           t_0_yzz_x_yy, t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_xx, t_0_yzz_y_xy,\
                           t_0_yzz_y_xz, t_0_yzz_y_yy, t_0_yzz_y_yz, t_0_yzz_y_zz, t_0_yzz_z_xx,\
                           t_0_yzz_z_xy, t_0_yzz_z_xz, t_0_yzz_z_yy, t_0_yzz_z_yz, t_0_yzz_z_zz,\
                           t_0_zzz_0_xx, t_0_zzz_0_xxx, t_0_zzz_0_xxy, t_0_zzz_0_xxz,\
                           t_0_zzz_0_xy, t_0_zzz_0_xyy, t_0_zzz_0_xyz, t_0_zzz_0_xz,\
                           t_0_zzz_0_xzz, t_0_zzz_0_yy, t_0_zzz_0_yyy, t_0_zzz_0_yyz,\
                           t_0_zzz_0_yz, t_0_zzz_0_yzz, t_0_zzz_0_zz, t_0_zzz_0_zzz,\
                           t_0_zzz_x_xx, t_0_zzz_x_xy, t_0_zzz_x_xz, t_0_zzz_x_yy, t_0_zzz_x_yz,\
                           t_0_zzz_x_zz, t_0_zzz_y_xx, t_0_zzz_y_xy, t_0_zzz_y_xz, t_0_zzz_y_yy,\
                           t_0_zzz_y_yz, t_0_zzz_y_zz, t_0_zzz_z_xx, t_0_zzz_z_xy, t_0_zzz_z_xz,\
                           t_0_zzz_z_yy, t_0_zzz_z_yz, t_0_zzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_zz[i] = t_0_zzz_0_zzz[i] - rcd_z[i] * t_0_zzz_0_zz[i];

        t_0_zzz_z_yz[i] = t_0_zzz_0_yzz[i] - rcd_z[i] * t_0_zzz_0_yz[i];

        t_0_zzz_z_yy[i] = t_0_zzz_0_yyz[i] - rcd_z[i] * t_0_zzz_0_yy[i];

        t_0_zzz_z_xz[i] = t_0_zzz_0_xzz[i] - rcd_z[i] * t_0_zzz_0_xz[i];

        t_0_zzz_z_xy[i] = t_0_zzz_0_xyz[i] - rcd_z[i] * t_0_zzz_0_xy[i];

        t_0_zzz_z_xx[i] = t_0_zzz_0_xxz[i] - rcd_z[i] * t_0_zzz_0_xx[i];

        t_0_zzz_y_zz[i] = t_0_zzz_0_yzz[i] - rcd_y[i] * t_0_zzz_0_zz[i];

        t_0_zzz_y_yz[i] = t_0_zzz_0_yyz[i] - rcd_y[i] * t_0_zzz_0_yz[i];

        t_0_zzz_y_yy[i] = t_0_zzz_0_yyy[i] - rcd_y[i] * t_0_zzz_0_yy[i];

        t_0_zzz_y_xz[i] = t_0_zzz_0_xyz[i] - rcd_y[i] * t_0_zzz_0_xz[i];

        t_0_zzz_y_xy[i] = t_0_zzz_0_xyy[i] - rcd_y[i] * t_0_zzz_0_xy[i];

        t_0_zzz_y_xx[i] = t_0_zzz_0_xxy[i] - rcd_y[i] * t_0_zzz_0_xx[i];

        t_0_zzz_x_zz[i] = t_0_zzz_0_xzz[i] - rcd_x[i] * t_0_zzz_0_zz[i];

        t_0_zzz_x_yz[i] = t_0_zzz_0_xyz[i] - rcd_x[i] * t_0_zzz_0_yz[i];

        t_0_zzz_x_yy[i] = t_0_zzz_0_xyy[i] - rcd_x[i] * t_0_zzz_0_yy[i];

        t_0_zzz_x_xz[i] = t_0_zzz_0_xxz[i] - rcd_x[i] * t_0_zzz_0_xz[i];

        t_0_zzz_x_xy[i] = t_0_zzz_0_xxy[i] - rcd_x[i] * t_0_zzz_0_xy[i];

        t_0_zzz_x_xx[i] = t_0_zzz_0_xxx[i] - rcd_x[i] * t_0_zzz_0_xx[i];

        t_0_yzz_z_zz[i] = t_0_yzz_0_zzz[i] - rcd_z[i] * t_0_yzz_0_zz[i];

        t_0_yzz_z_yz[i] = t_0_yzz_0_yzz[i] - rcd_z[i] * t_0_yzz_0_yz[i];

        t_0_yzz_z_yy[i] = t_0_yzz_0_yyz[i] - rcd_z[i] * t_0_yzz_0_yy[i];

        t_0_yzz_z_xz[i] = t_0_yzz_0_xzz[i] - rcd_z[i] * t_0_yzz_0_xz[i];

        t_0_yzz_z_xy[i] = t_0_yzz_0_xyz[i] - rcd_z[i] * t_0_yzz_0_xy[i];

        t_0_yzz_z_xx[i] = t_0_yzz_0_xxz[i] - rcd_z[i] * t_0_yzz_0_xx[i];

        t_0_yzz_y_zz[i] = t_0_yzz_0_yzz[i] - rcd_y[i] * t_0_yzz_0_zz[i];

        t_0_yzz_y_yz[i] = t_0_yzz_0_yyz[i] - rcd_y[i] * t_0_yzz_0_yz[i];

        t_0_yzz_y_yy[i] = t_0_yzz_0_yyy[i] - rcd_y[i] * t_0_yzz_0_yy[i];

        t_0_yzz_y_xz[i] = t_0_yzz_0_xyz[i] - rcd_y[i] * t_0_yzz_0_xz[i];

        t_0_yzz_y_xy[i] = t_0_yzz_0_xyy[i] - rcd_y[i] * t_0_yzz_0_xy[i];

        t_0_yzz_y_xx[i] = t_0_yzz_0_xxy[i] - rcd_y[i] * t_0_yzz_0_xx[i];

        t_0_yzz_x_zz[i] = t_0_yzz_0_xzz[i] - rcd_x[i] * t_0_yzz_0_zz[i];

        t_0_yzz_x_yz[i] = t_0_yzz_0_xyz[i] - rcd_x[i] * t_0_yzz_0_yz[i];

        t_0_yzz_x_yy[i] = t_0_yzz_0_xyy[i] - rcd_x[i] * t_0_yzz_0_yy[i];

        t_0_yzz_x_xz[i] = t_0_yzz_0_xxz[i] - rcd_x[i] * t_0_yzz_0_xz[i];

        t_0_yzz_x_xy[i] = t_0_yzz_0_xxy[i] - rcd_x[i] * t_0_yzz_0_xy[i];

        t_0_yzz_x_xx[i] = t_0_yzz_0_xxx[i] - rcd_x[i] * t_0_yzz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyy_0_xx, t_0_yyy_0_xxx, t_0_yyy_0_xxy,\
                           t_0_yyy_0_xxz, t_0_yyy_0_xy, t_0_yyy_0_xyy, t_0_yyy_0_xyz,\
                           t_0_yyy_0_xz, t_0_yyy_0_xzz, t_0_yyy_0_yy, t_0_yyy_0_yyy,\
                           t_0_yyy_0_yyz, t_0_yyy_0_yz, t_0_yyy_0_yzz, t_0_yyy_0_zz,\
                           t_0_yyy_0_zzz, t_0_yyy_x_xx, t_0_yyy_x_xy, t_0_yyy_x_xz,\
                           t_0_yyy_x_yy, t_0_yyy_x_yz, t_0_yyy_x_zz, t_0_yyy_y_xx, t_0_yyy_y_xy,\
                           t_0_yyy_y_xz, t_0_yyy_y_yy, t_0_yyy_y_yz, t_0_yyy_y_zz, t_0_yyy_z_xx,\
                           t_0_yyy_z_xy, t_0_yyy_z_xz, t_0_yyy_z_yy, t_0_yyy_z_yz, t_0_yyy_z_zz,\
                           t_0_yyz_0_xx, t_0_yyz_0_xxx, t_0_yyz_0_xxy, t_0_yyz_0_xxz,\
                           t_0_yyz_0_xy, t_0_yyz_0_xyy, t_0_yyz_0_xyz, t_0_yyz_0_xz,\
                           t_0_yyz_0_xzz, t_0_yyz_0_yy, t_0_yyz_0_yyy, t_0_yyz_0_yyz,\
                           t_0_yyz_0_yz, t_0_yyz_0_yzz, t_0_yyz_0_zz, t_0_yyz_0_zzz,\
                           t_0_yyz_x_xx, t_0_yyz_x_xy, t_0_yyz_x_xz, t_0_yyz_x_yy, t_0_yyz_x_yz,\
                           t_0_yyz_x_zz, t_0_yyz_y_xx, t_0_yyz_y_xy, t_0_yyz_y_xz, t_0_yyz_y_yy,\
                           t_0_yyz_y_yz, t_0_yyz_y_zz, t_0_yyz_z_xx, t_0_yyz_z_xy, t_0_yyz_z_xz,\
                           t_0_yyz_z_yy, t_0_yyz_z_yz, t_0_yyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyz_z_zz[i] = t_0_yyz_0_zzz[i] - rcd_z[i] * t_0_yyz_0_zz[i];

        t_0_yyz_z_yz[i] = t_0_yyz_0_yzz[i] - rcd_z[i] * t_0_yyz_0_yz[i];

        t_0_yyz_z_yy[i] = t_0_yyz_0_yyz[i] - rcd_z[i] * t_0_yyz_0_yy[i];

        t_0_yyz_z_xz[i] = t_0_yyz_0_xzz[i] - rcd_z[i] * t_0_yyz_0_xz[i];

        t_0_yyz_z_xy[i] = t_0_yyz_0_xyz[i] - rcd_z[i] * t_0_yyz_0_xy[i];

        t_0_yyz_z_xx[i] = t_0_yyz_0_xxz[i] - rcd_z[i] * t_0_yyz_0_xx[i];

        t_0_yyz_y_zz[i] = t_0_yyz_0_yzz[i] - rcd_y[i] * t_0_yyz_0_zz[i];

        t_0_yyz_y_yz[i] = t_0_yyz_0_yyz[i] - rcd_y[i] * t_0_yyz_0_yz[i];

        t_0_yyz_y_yy[i] = t_0_yyz_0_yyy[i] - rcd_y[i] * t_0_yyz_0_yy[i];

        t_0_yyz_y_xz[i] = t_0_yyz_0_xyz[i] - rcd_y[i] * t_0_yyz_0_xz[i];

        t_0_yyz_y_xy[i] = t_0_yyz_0_xyy[i] - rcd_y[i] * t_0_yyz_0_xy[i];

        t_0_yyz_y_xx[i] = t_0_yyz_0_xxy[i] - rcd_y[i] * t_0_yyz_0_xx[i];

        t_0_yyz_x_zz[i] = t_0_yyz_0_xzz[i] - rcd_x[i] * t_0_yyz_0_zz[i];

        t_0_yyz_x_yz[i] = t_0_yyz_0_xyz[i] - rcd_x[i] * t_0_yyz_0_yz[i];

        t_0_yyz_x_yy[i] = t_0_yyz_0_xyy[i] - rcd_x[i] * t_0_yyz_0_yy[i];

        t_0_yyz_x_xz[i] = t_0_yyz_0_xxz[i] - rcd_x[i] * t_0_yyz_0_xz[i];

        t_0_yyz_x_xy[i] = t_0_yyz_0_xxy[i] - rcd_x[i] * t_0_yyz_0_xy[i];

        t_0_yyz_x_xx[i] = t_0_yyz_0_xxx[i] - rcd_x[i] * t_0_yyz_0_xx[i];

        t_0_yyy_z_zz[i] = t_0_yyy_0_zzz[i] - rcd_z[i] * t_0_yyy_0_zz[i];

        t_0_yyy_z_yz[i] = t_0_yyy_0_yzz[i] - rcd_z[i] * t_0_yyy_0_yz[i];

        t_0_yyy_z_yy[i] = t_0_yyy_0_yyz[i] - rcd_z[i] * t_0_yyy_0_yy[i];

        t_0_yyy_z_xz[i] = t_0_yyy_0_xzz[i] - rcd_z[i] * t_0_yyy_0_xz[i];

        t_0_yyy_z_xy[i] = t_0_yyy_0_xyz[i] - rcd_z[i] * t_0_yyy_0_xy[i];

        t_0_yyy_z_xx[i] = t_0_yyy_0_xxz[i] - rcd_z[i] * t_0_yyy_0_xx[i];

        t_0_yyy_y_zz[i] = t_0_yyy_0_yzz[i] - rcd_y[i] * t_0_yyy_0_zz[i];

        t_0_yyy_y_yz[i] = t_0_yyy_0_yyz[i] - rcd_y[i] * t_0_yyy_0_yz[i];

        t_0_yyy_y_yy[i] = t_0_yyy_0_yyy[i] - rcd_y[i] * t_0_yyy_0_yy[i];

        t_0_yyy_y_xz[i] = t_0_yyy_0_xyz[i] - rcd_y[i] * t_0_yyy_0_xz[i];

        t_0_yyy_y_xy[i] = t_0_yyy_0_xyy[i] - rcd_y[i] * t_0_yyy_0_xy[i];

        t_0_yyy_y_xx[i] = t_0_yyy_0_xxy[i] - rcd_y[i] * t_0_yyy_0_xx[i];

        t_0_yyy_x_zz[i] = t_0_yyy_0_xzz[i] - rcd_x[i] * t_0_yyy_0_zz[i];

        t_0_yyy_x_yz[i] = t_0_yyy_0_xyz[i] - rcd_x[i] * t_0_yyy_0_yz[i];

        t_0_yyy_x_yy[i] = t_0_yyy_0_xyy[i] - rcd_x[i] * t_0_yyy_0_yy[i];

        t_0_yyy_x_xz[i] = t_0_yyy_0_xxz[i] - rcd_x[i] * t_0_yyy_0_xz[i];

        t_0_yyy_x_xy[i] = t_0_yyy_0_xxy[i] - rcd_x[i] * t_0_yyy_0_xy[i];

        t_0_yyy_x_xx[i] = t_0_yyy_0_xxx[i] - rcd_x[i] * t_0_yyy_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyz_0_xx, t_0_xyz_0_xxx, t_0_xyz_0_xxy,\
                           t_0_xyz_0_xxz, t_0_xyz_0_xy, t_0_xyz_0_xyy, t_0_xyz_0_xyz,\
                           t_0_xyz_0_xz, t_0_xyz_0_xzz, t_0_xyz_0_yy, t_0_xyz_0_yyy,\
                           t_0_xyz_0_yyz, t_0_xyz_0_yz, t_0_xyz_0_yzz, t_0_xyz_0_zz,\
                           t_0_xyz_0_zzz, t_0_xyz_x_xx, t_0_xyz_x_xy, t_0_xyz_x_xz,\
                           t_0_xyz_x_yy, t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx, t_0_xyz_y_xy,\
                           t_0_xyz_y_xz, t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz, t_0_xyz_z_xx,\
                           t_0_xyz_z_xy, t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz, t_0_xyz_z_zz,\
                           t_0_xzz_0_xx, t_0_xzz_0_xxx, t_0_xzz_0_xxy, t_0_xzz_0_xxz,\
                           t_0_xzz_0_xy, t_0_xzz_0_xyy, t_0_xzz_0_xyz, t_0_xzz_0_xz,\
                           t_0_xzz_0_xzz, t_0_xzz_0_yy, t_0_xzz_0_yyy, t_0_xzz_0_yyz,\
                           t_0_xzz_0_yz, t_0_xzz_0_yzz, t_0_xzz_0_zz, t_0_xzz_0_zzz,\
                           t_0_xzz_x_xx, t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy, t_0_xzz_x_yz,\
                           t_0_xzz_x_zz, t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz, t_0_xzz_y_yy,\
                           t_0_xzz_y_yz, t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy, t_0_xzz_z_xz,\
                           t_0_xzz_z_yy, t_0_xzz_z_yz, t_0_xzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzz_z_zz[i] = t_0_xzz_0_zzz[i] - rcd_z[i] * t_0_xzz_0_zz[i];

        t_0_xzz_z_yz[i] = t_0_xzz_0_yzz[i] - rcd_z[i] * t_0_xzz_0_yz[i];

        t_0_xzz_z_yy[i] = t_0_xzz_0_yyz[i] - rcd_z[i] * t_0_xzz_0_yy[i];

        t_0_xzz_z_xz[i] = t_0_xzz_0_xzz[i] - rcd_z[i] * t_0_xzz_0_xz[i];

        t_0_xzz_z_xy[i] = t_0_xzz_0_xyz[i] - rcd_z[i] * t_0_xzz_0_xy[i];

        t_0_xzz_z_xx[i] = t_0_xzz_0_xxz[i] - rcd_z[i] * t_0_xzz_0_xx[i];

        t_0_xzz_y_zz[i] = t_0_xzz_0_yzz[i] - rcd_y[i] * t_0_xzz_0_zz[i];

        t_0_xzz_y_yz[i] = t_0_xzz_0_yyz[i] - rcd_y[i] * t_0_xzz_0_yz[i];

        t_0_xzz_y_yy[i] = t_0_xzz_0_yyy[i] - rcd_y[i] * t_0_xzz_0_yy[i];

        t_0_xzz_y_xz[i] = t_0_xzz_0_xyz[i] - rcd_y[i] * t_0_xzz_0_xz[i];

        t_0_xzz_y_xy[i] = t_0_xzz_0_xyy[i] - rcd_y[i] * t_0_xzz_0_xy[i];

        t_0_xzz_y_xx[i] = t_0_xzz_0_xxy[i] - rcd_y[i] * t_0_xzz_0_xx[i];

        t_0_xzz_x_zz[i] = t_0_xzz_0_xzz[i] - rcd_x[i] * t_0_xzz_0_zz[i];

        t_0_xzz_x_yz[i] = t_0_xzz_0_xyz[i] - rcd_x[i] * t_0_xzz_0_yz[i];

        t_0_xzz_x_yy[i] = t_0_xzz_0_xyy[i] - rcd_x[i] * t_0_xzz_0_yy[i];

        t_0_xzz_x_xz[i] = t_0_xzz_0_xxz[i] - rcd_x[i] * t_0_xzz_0_xz[i];

        t_0_xzz_x_xy[i] = t_0_xzz_0_xxy[i] - rcd_x[i] * t_0_xzz_0_xy[i];

        t_0_xzz_x_xx[i] = t_0_xzz_0_xxx[i] - rcd_x[i] * t_0_xzz_0_xx[i];

        t_0_xyz_z_zz[i] = t_0_xyz_0_zzz[i] - rcd_z[i] * t_0_xyz_0_zz[i];

        t_0_xyz_z_yz[i] = t_0_xyz_0_yzz[i] - rcd_z[i] * t_0_xyz_0_yz[i];

        t_0_xyz_z_yy[i] = t_0_xyz_0_yyz[i] - rcd_z[i] * t_0_xyz_0_yy[i];

        t_0_xyz_z_xz[i] = t_0_xyz_0_xzz[i] - rcd_z[i] * t_0_xyz_0_xz[i];

        t_0_xyz_z_xy[i] = t_0_xyz_0_xyz[i] - rcd_z[i] * t_0_xyz_0_xy[i];

        t_0_xyz_z_xx[i] = t_0_xyz_0_xxz[i] - rcd_z[i] * t_0_xyz_0_xx[i];

        t_0_xyz_y_zz[i] = t_0_xyz_0_yzz[i] - rcd_y[i] * t_0_xyz_0_zz[i];

        t_0_xyz_y_yz[i] = t_0_xyz_0_yyz[i] - rcd_y[i] * t_0_xyz_0_yz[i];

        t_0_xyz_y_yy[i] = t_0_xyz_0_yyy[i] - rcd_y[i] * t_0_xyz_0_yy[i];

        t_0_xyz_y_xz[i] = t_0_xyz_0_xyz[i] - rcd_y[i] * t_0_xyz_0_xz[i];

        t_0_xyz_y_xy[i] = t_0_xyz_0_xyy[i] - rcd_y[i] * t_0_xyz_0_xy[i];

        t_0_xyz_y_xx[i] = t_0_xyz_0_xxy[i] - rcd_y[i] * t_0_xyz_0_xx[i];

        t_0_xyz_x_zz[i] = t_0_xyz_0_xzz[i] - rcd_x[i] * t_0_xyz_0_zz[i];

        t_0_xyz_x_yz[i] = t_0_xyz_0_xyz[i] - rcd_x[i] * t_0_xyz_0_yz[i];

        t_0_xyz_x_yy[i] = t_0_xyz_0_xyy[i] - rcd_x[i] * t_0_xyz_0_yy[i];

        t_0_xyz_x_xz[i] = t_0_xyz_0_xxz[i] - rcd_x[i] * t_0_xyz_0_xz[i];

        t_0_xyz_x_xy[i] = t_0_xyz_0_xxy[i] - rcd_x[i] * t_0_xyz_0_xy[i];

        t_0_xyz_x_xx[i] = t_0_xyz_0_xxx[i] - rcd_x[i] * t_0_xyz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxz_0_xx, t_0_xxz_0_xxx, t_0_xxz_0_xxy,\
                           t_0_xxz_0_xxz, t_0_xxz_0_xy, t_0_xxz_0_xyy, t_0_xxz_0_xyz,\
                           t_0_xxz_0_xz, t_0_xxz_0_xzz, t_0_xxz_0_yy, t_0_xxz_0_yyy,\
                           t_0_xxz_0_yyz, t_0_xxz_0_yz, t_0_xxz_0_yzz, t_0_xxz_0_zz,\
                           t_0_xxz_0_zzz, t_0_xxz_x_xx, t_0_xxz_x_xy, t_0_xxz_x_xz,\
                           t_0_xxz_x_yy, t_0_xxz_x_yz, t_0_xxz_x_zz, t_0_xxz_y_xx, t_0_xxz_y_xy,\
                           t_0_xxz_y_xz, t_0_xxz_y_yy, t_0_xxz_y_yz, t_0_xxz_y_zz, t_0_xxz_z_xx,\
                           t_0_xxz_z_xy, t_0_xxz_z_xz, t_0_xxz_z_yy, t_0_xxz_z_yz, t_0_xxz_z_zz,\
                           t_0_xyy_0_xx, t_0_xyy_0_xxx, t_0_xyy_0_xxy, t_0_xyy_0_xxz,\
                           t_0_xyy_0_xy, t_0_xyy_0_xyy, t_0_xyy_0_xyz, t_0_xyy_0_xz,\
                           t_0_xyy_0_xzz, t_0_xyy_0_yy, t_0_xyy_0_yyy, t_0_xyy_0_yyz,\
                           t_0_xyy_0_yz, t_0_xyy_0_yzz, t_0_xyy_0_zz, t_0_xyy_0_zzz,\
                           t_0_xyy_x_xx, t_0_xyy_x_xy, t_0_xyy_x_xz, t_0_xyy_x_yy, t_0_xyy_x_yz,\
                           t_0_xyy_x_zz, t_0_xyy_y_xx, t_0_xyy_y_xy, t_0_xyy_y_xz, t_0_xyy_y_yy,\
                           t_0_xyy_y_yz, t_0_xyy_y_zz, t_0_xyy_z_xx, t_0_xyy_z_xy, t_0_xyy_z_xz,\
                           t_0_xyy_z_yy, t_0_xyy_z_yz, t_0_xyy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyy_z_zz[i] = t_0_xyy_0_zzz[i] - rcd_z[i] * t_0_xyy_0_zz[i];

        t_0_xyy_z_yz[i] = t_0_xyy_0_yzz[i] - rcd_z[i] * t_0_xyy_0_yz[i];

        t_0_xyy_z_yy[i] = t_0_xyy_0_yyz[i] - rcd_z[i] * t_0_xyy_0_yy[i];

        t_0_xyy_z_xz[i] = t_0_xyy_0_xzz[i] - rcd_z[i] * t_0_xyy_0_xz[i];

        t_0_xyy_z_xy[i] = t_0_xyy_0_xyz[i] - rcd_z[i] * t_0_xyy_0_xy[i];

        t_0_xyy_z_xx[i] = t_0_xyy_0_xxz[i] - rcd_z[i] * t_0_xyy_0_xx[i];

        t_0_xyy_y_zz[i] = t_0_xyy_0_yzz[i] - rcd_y[i] * t_0_xyy_0_zz[i];

        t_0_xyy_y_yz[i] = t_0_xyy_0_yyz[i] - rcd_y[i] * t_0_xyy_0_yz[i];

        t_0_xyy_y_yy[i] = t_0_xyy_0_yyy[i] - rcd_y[i] * t_0_xyy_0_yy[i];

        t_0_xyy_y_xz[i] = t_0_xyy_0_xyz[i] - rcd_y[i] * t_0_xyy_0_xz[i];

        t_0_xyy_y_xy[i] = t_0_xyy_0_xyy[i] - rcd_y[i] * t_0_xyy_0_xy[i];

        t_0_xyy_y_xx[i] = t_0_xyy_0_xxy[i] - rcd_y[i] * t_0_xyy_0_xx[i];

        t_0_xyy_x_zz[i] = t_0_xyy_0_xzz[i] - rcd_x[i] * t_0_xyy_0_zz[i];

        t_0_xyy_x_yz[i] = t_0_xyy_0_xyz[i] - rcd_x[i] * t_0_xyy_0_yz[i];

        t_0_xyy_x_yy[i] = t_0_xyy_0_xyy[i] - rcd_x[i] * t_0_xyy_0_yy[i];

        t_0_xyy_x_xz[i] = t_0_xyy_0_xxz[i] - rcd_x[i] * t_0_xyy_0_xz[i];

        t_0_xyy_x_xy[i] = t_0_xyy_0_xxy[i] - rcd_x[i] * t_0_xyy_0_xy[i];

        t_0_xyy_x_xx[i] = t_0_xyy_0_xxx[i] - rcd_x[i] * t_0_xyy_0_xx[i];

        t_0_xxz_z_zz[i] = t_0_xxz_0_zzz[i] - rcd_z[i] * t_0_xxz_0_zz[i];

        t_0_xxz_z_yz[i] = t_0_xxz_0_yzz[i] - rcd_z[i] * t_0_xxz_0_yz[i];

        t_0_xxz_z_yy[i] = t_0_xxz_0_yyz[i] - rcd_z[i] * t_0_xxz_0_yy[i];

        t_0_xxz_z_xz[i] = t_0_xxz_0_xzz[i] - rcd_z[i] * t_0_xxz_0_xz[i];

        t_0_xxz_z_xy[i] = t_0_xxz_0_xyz[i] - rcd_z[i] * t_0_xxz_0_xy[i];

        t_0_xxz_z_xx[i] = t_0_xxz_0_xxz[i] - rcd_z[i] * t_0_xxz_0_xx[i];

        t_0_xxz_y_zz[i] = t_0_xxz_0_yzz[i] - rcd_y[i] * t_0_xxz_0_zz[i];

        t_0_xxz_y_yz[i] = t_0_xxz_0_yyz[i] - rcd_y[i] * t_0_xxz_0_yz[i];

        t_0_xxz_y_yy[i] = t_0_xxz_0_yyy[i] - rcd_y[i] * t_0_xxz_0_yy[i];

        t_0_xxz_y_xz[i] = t_0_xxz_0_xyz[i] - rcd_y[i] * t_0_xxz_0_xz[i];

        t_0_xxz_y_xy[i] = t_0_xxz_0_xyy[i] - rcd_y[i] * t_0_xxz_0_xy[i];

        t_0_xxz_y_xx[i] = t_0_xxz_0_xxy[i] - rcd_y[i] * t_0_xxz_0_xx[i];

        t_0_xxz_x_zz[i] = t_0_xxz_0_xzz[i] - rcd_x[i] * t_0_xxz_0_zz[i];

        t_0_xxz_x_yz[i] = t_0_xxz_0_xyz[i] - rcd_x[i] * t_0_xxz_0_yz[i];

        t_0_xxz_x_yy[i] = t_0_xxz_0_xyy[i] - rcd_x[i] * t_0_xxz_0_yy[i];

        t_0_xxz_x_xz[i] = t_0_xxz_0_xxz[i] - rcd_x[i] * t_0_xxz_0_xz[i];

        t_0_xxz_x_xy[i] = t_0_xxz_0_xxy[i] - rcd_x[i] * t_0_xxz_0_xy[i];

        t_0_xxz_x_xx[i] = t_0_xxz_0_xxx[i] - rcd_x[i] * t_0_xxz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_0_xx, t_0_xxx_0_xxx, t_0_xxx_0_xxy,\
                           t_0_xxx_0_xxz, t_0_xxx_0_xy, t_0_xxx_0_xyy, t_0_xxx_0_xyz,\
                           t_0_xxx_0_xz, t_0_xxx_0_xzz, t_0_xxx_0_yy, t_0_xxx_0_yyy,\
                           t_0_xxx_0_yyz, t_0_xxx_0_yz, t_0_xxx_0_yzz, t_0_xxx_0_zz,\
                           t_0_xxx_0_zzz, t_0_xxx_x_xx, t_0_xxx_x_xy, t_0_xxx_x_xz,\
                           t_0_xxx_x_yy, t_0_xxx_x_yz, t_0_xxx_x_zz, t_0_xxx_y_xx, t_0_xxx_y_xy,\
                           t_0_xxx_y_xz, t_0_xxx_y_yy, t_0_xxx_y_yz, t_0_xxx_y_zz, t_0_xxx_z_xx,\
                           t_0_xxx_z_xy, t_0_xxx_z_xz, t_0_xxx_z_yy, t_0_xxx_z_yz, t_0_xxx_z_zz,\
                           t_0_xxy_0_xx, t_0_xxy_0_xxx, t_0_xxy_0_xxy, t_0_xxy_0_xxz,\
                           t_0_xxy_0_xy, t_0_xxy_0_xyy, t_0_xxy_0_xyz, t_0_xxy_0_xz,\
                           t_0_xxy_0_xzz, t_0_xxy_0_yy, t_0_xxy_0_yyy, t_0_xxy_0_yyz,\
                           t_0_xxy_0_yz, t_0_xxy_0_yzz, t_0_xxy_0_zz, t_0_xxy_0_zzz,\
                           t_0_xxy_x_xx, t_0_xxy_x_xy, t_0_xxy_x_xz, t_0_xxy_x_yy, t_0_xxy_x_yz,\
                           t_0_xxy_x_zz, t_0_xxy_y_xx, t_0_xxy_y_xy, t_0_xxy_y_xz, t_0_xxy_y_yy,\
                           t_0_xxy_y_yz, t_0_xxy_y_zz, t_0_xxy_z_xx, t_0_xxy_z_xy, t_0_xxy_z_xz,\
                           t_0_xxy_z_yy, t_0_xxy_z_yz, t_0_xxy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxy_z_zz[i] = t_0_xxy_0_zzz[i] - rcd_z[i] * t_0_xxy_0_zz[i];

        t_0_xxy_z_yz[i] = t_0_xxy_0_yzz[i] - rcd_z[i] * t_0_xxy_0_yz[i];

        t_0_xxy_z_yy[i] = t_0_xxy_0_yyz[i] - rcd_z[i] * t_0_xxy_0_yy[i];

        t_0_xxy_z_xz[i] = t_0_xxy_0_xzz[i] - rcd_z[i] * t_0_xxy_0_xz[i];

        t_0_xxy_z_xy[i] = t_0_xxy_0_xyz[i] - rcd_z[i] * t_0_xxy_0_xy[i];

        t_0_xxy_z_xx[i] = t_0_xxy_0_xxz[i] - rcd_z[i] * t_0_xxy_0_xx[i];

        t_0_xxy_y_zz[i] = t_0_xxy_0_yzz[i] - rcd_y[i] * t_0_xxy_0_zz[i];

        t_0_xxy_y_yz[i] = t_0_xxy_0_yyz[i] - rcd_y[i] * t_0_xxy_0_yz[i];

        t_0_xxy_y_yy[i] = t_0_xxy_0_yyy[i] - rcd_y[i] * t_0_xxy_0_yy[i];

        t_0_xxy_y_xz[i] = t_0_xxy_0_xyz[i] - rcd_y[i] * t_0_xxy_0_xz[i];

        t_0_xxy_y_xy[i] = t_0_xxy_0_xyy[i] - rcd_y[i] * t_0_xxy_0_xy[i];

        t_0_xxy_y_xx[i] = t_0_xxy_0_xxy[i] - rcd_y[i] * t_0_xxy_0_xx[i];

        t_0_xxy_x_zz[i] = t_0_xxy_0_xzz[i] - rcd_x[i] * t_0_xxy_0_zz[i];

        t_0_xxy_x_yz[i] = t_0_xxy_0_xyz[i] - rcd_x[i] * t_0_xxy_0_yz[i];

        t_0_xxy_x_yy[i] = t_0_xxy_0_xyy[i] - rcd_x[i] * t_0_xxy_0_yy[i];

        t_0_xxy_x_xz[i] = t_0_xxy_0_xxz[i] - rcd_x[i] * t_0_xxy_0_xz[i];

        t_0_xxy_x_xy[i] = t_0_xxy_0_xxy[i] - rcd_x[i] * t_0_xxy_0_xy[i];

        t_0_xxy_x_xx[i] = t_0_xxy_0_xxx[i] - rcd_x[i] * t_0_xxy_0_xx[i];

        t_0_xxx_z_zz[i] = t_0_xxx_0_zzz[i] - rcd_z[i] * t_0_xxx_0_zz[i];

        t_0_xxx_z_yz[i] = t_0_xxx_0_yzz[i] - rcd_z[i] * t_0_xxx_0_yz[i];

        t_0_xxx_z_yy[i] = t_0_xxx_0_yyz[i] - rcd_z[i] * t_0_xxx_0_yy[i];

        t_0_xxx_z_xz[i] = t_0_xxx_0_xzz[i] - rcd_z[i] * t_0_xxx_0_xz[i];

        t_0_xxx_z_xy[i] = t_0_xxx_0_xyz[i] - rcd_z[i] * t_0_xxx_0_xy[i];

        t_0_xxx_z_xx[i] = t_0_xxx_0_xxz[i] - rcd_z[i] * t_0_xxx_0_xx[i];

        t_0_xxx_y_zz[i] = t_0_xxx_0_yzz[i] - rcd_y[i] * t_0_xxx_0_zz[i];

        t_0_xxx_y_yz[i] = t_0_xxx_0_yyz[i] - rcd_y[i] * t_0_xxx_0_yz[i];

        t_0_xxx_y_yy[i] = t_0_xxx_0_yyy[i] - rcd_y[i] * t_0_xxx_0_yy[i];

        t_0_xxx_y_xz[i] = t_0_xxx_0_xyz[i] - rcd_y[i] * t_0_xxx_0_xz[i];

        t_0_xxx_y_xy[i] = t_0_xxx_0_xyy[i] - rcd_y[i] * t_0_xxx_0_xy[i];

        t_0_xxx_y_xx[i] = t_0_xxx_0_xxy[i] - rcd_y[i] * t_0_xxx_0_xx[i];

        t_0_xxx_x_zz[i] = t_0_xxx_0_xzz[i] - rcd_x[i] * t_0_xxx_0_zz[i];

        t_0_xxx_x_yz[i] = t_0_xxx_0_xyz[i] - rcd_x[i] * t_0_xxx_0_yz[i];

        t_0_xxx_x_yy[i] = t_0_xxx_0_xyy[i] - rcd_x[i] * t_0_xxx_0_yy[i];

        t_0_xxx_x_xz[i] = t_0_xxx_0_xxz[i] - rcd_x[i] * t_0_xxx_0_xz[i];

        t_0_xxx_x_xy[i] = t_0_xxx_0_xxy[i] - rcd_x[i] * t_0_xxx_0_xy[i];

        t_0_xxx_x_xx[i] = t_0_xxx_0_xxx[i] - rcd_x[i] * t_0_xxx_0_xx[i];
    }
}


} // derirec namespace
