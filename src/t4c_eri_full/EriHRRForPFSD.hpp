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
compHostHRRForPFSD_V0(      BufferHostXY<T>&      intsBufferPFSD,
                      const BufferHostX<int32_t>& intsIndexesPFSD,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PFSD) integral components

    t_z_zzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(0));

    t_z_zzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(1));

    t_z_zzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(2));

    t_z_zzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(3));

    t_z_zzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(4));

    t_z_zzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(5));

    t_z_yzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(6));

    t_z_yzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(7));

    t_z_yzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(8));

    t_z_yzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(9));

    t_z_yzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(10));

    t_z_yzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(11));

    t_z_yyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(12));

    t_z_yyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(13));

    t_z_yyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(14));

    t_z_yyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(15));

    t_z_yyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(16));

    t_z_yyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(17));

    t_z_xzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(18));

    t_z_xzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(19));

    t_z_xzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(20));

    t_z_xzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(21));

    t_z_xzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(22));

    t_z_xzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(23));

    t_z_xyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(24));

    t_z_xyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(25));

    t_z_xyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(26));

    t_z_xyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(27));

    t_z_xyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(28));

    t_z_xyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(29));

    t_z_xxz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(30));

    t_z_xxz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(31));

    t_z_xxz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(32));

    t_z_xxz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(33));

    t_z_xxz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(34));

    t_z_xxz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(35));

    t_y_zzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(36));

    t_y_zzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(37));

    t_y_zzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(38));

    t_y_zzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(39));

    t_y_zzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(40));

    t_y_zzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(41));

    t_y_yzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(42));

    t_y_yzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(43));

    t_y_yzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(44));

    t_y_yzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(45));

    t_y_yzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(46));

    t_y_yzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(47));

    t_y_yyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(48));

    t_y_yyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(49));

    t_y_yyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(50));

    t_y_yyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(51));

    t_y_yyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(52));

    t_y_yyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(53));

    t_y_yyy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(54));

    t_y_yyy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(55));

    t_y_yyy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(56));

    t_y_yyy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(57));

    t_y_yyy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(58));

    t_y_yyy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(59));

    t_y_xzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(60));

    t_y_xzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(61));

    t_y_xzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(62));

    t_y_xzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(63));

    t_y_xzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(64));

    t_y_xzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(65));

    t_y_xyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(66));

    t_y_xyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(67));

    t_y_xyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(68));

    t_y_xyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(69));

    t_y_xyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(70));

    t_y_xyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(71));

    t_y_xyy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(72));

    t_y_xyy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(73));

    t_y_xyy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(74));

    t_y_xyy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(75));

    t_y_xyy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(76));

    t_y_xyy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(77));

    t_y_xxz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(78));

    t_y_xxz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(79));

    t_y_xxz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(80));

    t_y_xxz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(81));

    t_y_xxz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(82));

    t_y_xxz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(83));

    t_y_xxy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(84));

    t_y_xxy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(85));

    t_y_xxy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(86));

    t_y_xxy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(87));

    t_y_xxy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(88));

    t_y_xxy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(89));

    t_x_zzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(90));

    t_x_zzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(91));

    t_x_zzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(92));

    t_x_zzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(93));

    t_x_zzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(94));

    t_x_zzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(95));

    t_x_yzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(96));

    t_x_yzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(97));

    t_x_yzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(98));

    t_x_yzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(99));

    t_x_yzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(100));

    t_x_yzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(101));

    t_x_yyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(102));

    t_x_yyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(103));

    t_x_yyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(104));

    t_x_yyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(105));

    t_x_yyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(106));

    t_x_yyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(107));

    t_x_yyy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(108));

    t_x_yyy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(109));

    t_x_yyy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(110));

    t_x_yyy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(111));

    t_x_yyy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(112));

    t_x_yyy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(113));

    t_x_xzz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(114));

    t_x_xzz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(115));

    t_x_xzz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(116));

    t_x_xzz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(117));

    t_x_xzz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(118));

    t_x_xzz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(119));

    t_x_xyz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(120));

    t_x_xyz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(121));

    t_x_xyz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(122));

    t_x_xyz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(123));

    t_x_xyz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(124));

    t_x_xyz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(125));

    t_x_xyy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(126));

    t_x_xyy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(127));

    t_x_xyy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(128));

    t_x_xyy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(129));

    t_x_xyy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(130));

    t_x_xyy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(131));

    t_x_xxz_0_zz = intsBufferPFSD.data(intsIndexesPFSD(132));

    t_x_xxz_0_yz = intsBufferPFSD.data(intsIndexesPFSD(133));

    t_x_xxz_0_yy = intsBufferPFSD.data(intsIndexesPFSD(134));

    t_x_xxz_0_xz = intsBufferPFSD.data(intsIndexesPFSD(135));

    t_x_xxz_0_xy = intsBufferPFSD.data(intsIndexesPFSD(136));

    t_x_xxz_0_xx = intsBufferPFSD.data(intsIndexesPFSD(137));

    t_x_xxy_0_zz = intsBufferPFSD.data(intsIndexesPFSD(138));

    t_x_xxy_0_yz = intsBufferPFSD.data(intsIndexesPFSD(139));

    t_x_xxy_0_yy = intsBufferPFSD.data(intsIndexesPFSD(140));

    t_x_xxy_0_xz = intsBufferPFSD.data(intsIndexesPFSD(141));

    t_x_xxy_0_xy = intsBufferPFSD.data(intsIndexesPFSD(142));

    t_x_xxy_0_xx = intsBufferPFSD.data(intsIndexesPFSD(143));

    t_x_xxx_0_zz = intsBufferPFSD.data(intsIndexesPFSD(144));

    t_x_xxx_0_yz = intsBufferPFSD.data(intsIndexesPFSD(145));

    t_x_xxx_0_yy = intsBufferPFSD.data(intsIndexesPFSD(146));

    t_x_xxx_0_xz = intsBufferPFSD.data(intsIndexesPFSD(147));

    t_x_xxx_0_xy = intsBufferPFSD.data(intsIndexesPFSD(148));

    t_x_xxx_0_xx = intsBufferPFSD.data(intsIndexesPFSD(149));

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

    // set up (SGSD) integral components

    t_0_zzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(0));

    t_0_zzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(1));

    t_0_zzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(2));

    t_0_zzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(3));

    t_0_zzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(4));

    t_0_zzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(5));

    t_0_yzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(6));

    t_0_yzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(7));

    t_0_yzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(8));

    t_0_yzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(9));

    t_0_yzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(10));

    t_0_yzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(11));

    t_0_yyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(12));

    t_0_yyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(13));

    t_0_yyzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(14));

    t_0_yyzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(15));

    t_0_yyzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(16));

    t_0_yyzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(17));

    t_0_yyyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(18));

    t_0_yyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(19));

    t_0_yyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(20));

    t_0_yyyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(21));

    t_0_yyyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(22));

    t_0_yyyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(23));

    t_0_yyyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(24));

    t_0_yyyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(25));

    t_0_yyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(26));

    t_0_yyyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(27));

    t_0_yyyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(28));

    t_0_yyyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(29));

    t_0_xzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(30));

    t_0_xzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(31));

    t_0_xzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(32));

    t_0_xzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(33));

    t_0_xzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(34));

    t_0_xzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(35));

    t_0_xyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(36));

    t_0_xyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(37));

    t_0_xyzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(38));

    t_0_xyzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(39));

    t_0_xyzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(40));

    t_0_xyzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(41));

    t_0_xyyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(42));

    t_0_xyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(43));

    t_0_xyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(44));

    t_0_xyyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(45));

    t_0_xyyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(46));

    t_0_xyyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(47));

    t_0_xyyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(48));

    t_0_xyyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(49));

    t_0_xyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(50));

    t_0_xyyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(51));

    t_0_xyyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(52));

    t_0_xyyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(53));

    t_0_xxzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(54));

    t_0_xxzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(55));

    t_0_xxzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(56));

    t_0_xxzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(57));

    t_0_xxzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(58));

    t_0_xxzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(59));

    t_0_xxyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(60));

    t_0_xxyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(61));

    t_0_xxyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(62));

    t_0_xxyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(63));

    t_0_xxyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(64));

    t_0_xxyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(65));

    t_0_xxyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(66));

    t_0_xxyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(67));

    t_0_xxyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(68));

    t_0_xxyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(69));

    t_0_xxyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(70));

    t_0_xxyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(71));

    t_0_xxxz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(72));

    t_0_xxxz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(73));

    t_0_xxxz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(74));

    t_0_xxxz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(75));

    t_0_xxxz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(76));

    t_0_xxxz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(77));

    t_0_xxxy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(78));

    t_0_xxxy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(79));

    t_0_xxxy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(80));

    t_0_xxxy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(81));

    t_0_xxxy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(82));

    t_0_xxxy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(83));

    t_0_xxxx_0_zz = intsBufferSGSD.data(intsIndexesSGSD(84));

    t_0_xxxx_0_yz = intsBufferSGSD.data(intsIndexesSGSD(85));

    t_0_xxxx_0_yy = intsBufferSGSD.data(intsIndexesSGSD(86));

    t_0_xxxx_0_xz = intsBufferSGSD.data(intsIndexesSGSD(87));

    t_0_xxxx_0_xy = intsBufferSGSD.data(intsIndexesSGSD(88));

    t_0_xxxx_0_xx = intsBufferSGSD.data(intsIndexesSGSD(89));

    #pragma omp simd align(rab_z, t_0_xxz_0_xx, t_0_xxz_0_xy, t_0_xxz_0_xz, t_0_xxz_0_yy,\
                           t_0_xxz_0_yz, t_0_xxz_0_zz, t_0_xxzz_0_xx, t_0_xxzz_0_xy,\
                           t_0_xxzz_0_xz, t_0_xxzz_0_yy, t_0_xxzz_0_yz, t_0_xxzz_0_zz,\
                           t_0_xyz_0_xx, t_0_xyz_0_xy, t_0_xyz_0_xz, t_0_xyz_0_yy, t_0_xyz_0_yz,\
                           t_0_xyz_0_zz, t_0_xyzz_0_xx, t_0_xyzz_0_xy, t_0_xyzz_0_xz,\
                           t_0_xyzz_0_yy, t_0_xyzz_0_yz, t_0_xyzz_0_zz, t_0_xzz_0_xx,\
                           t_0_xzz_0_xy, t_0_xzz_0_xz, t_0_xzz_0_yy, t_0_xzz_0_yz, t_0_xzz_0_zz,\
                           t_0_xzzz_0_xx, t_0_xzzz_0_xy, t_0_xzzz_0_xz, t_0_xzzz_0_yy,\
                           t_0_xzzz_0_yz, t_0_xzzz_0_zz, t_0_yyz_0_xx, t_0_yyz_0_xy,\
                           t_0_yyz_0_xz, t_0_yyz_0_yy, t_0_yyz_0_yz, t_0_yyz_0_zz, t_0_yyzz_0_xx,\
                           t_0_yyzz_0_xy, t_0_yyzz_0_xz, t_0_yyzz_0_yy, t_0_yyzz_0_yz,\
                           t_0_yyzz_0_zz, t_0_yzz_0_xx, t_0_yzz_0_xy, t_0_yzz_0_xz,\
                           t_0_yzz_0_yy, t_0_yzz_0_yz, t_0_yzz_0_zz, t_0_yzzz_0_xx,\
                           t_0_yzzz_0_xy, t_0_yzzz_0_xz, t_0_yzzz_0_yy, t_0_yzzz_0_yz,\
                           t_0_yzzz_0_zz, t_0_zzz_0_xx, t_0_zzz_0_xy, t_0_zzz_0_xz,\
                           t_0_zzz_0_yy, t_0_zzz_0_yz, t_0_zzz_0_zz, t_0_zzzz_0_xx,\
                           t_0_zzzz_0_xy, t_0_zzzz_0_xz, t_0_zzzz_0_yy, t_0_zzzz_0_yz,\
                           t_0_zzzz_0_zz, t_z_xxz_0_xx, t_z_xxz_0_xy, t_z_xxz_0_xz,\
                           t_z_xxz_0_yy, t_z_xxz_0_yz, t_z_xxz_0_zz, t_z_xyz_0_xx, t_z_xyz_0_xy,\
                           t_z_xyz_0_xz, t_z_xyz_0_yy, t_z_xyz_0_yz, t_z_xyz_0_zz, t_z_xzz_0_xx,\
                           t_z_xzz_0_xy, t_z_xzz_0_xz, t_z_xzz_0_yy, t_z_xzz_0_yz, t_z_xzz_0_zz,\
                           t_z_yyz_0_xx, t_z_yyz_0_xy, t_z_yyz_0_xz, t_z_yyz_0_yy, t_z_yyz_0_yz,\
                           t_z_yyz_0_zz, t_z_yzz_0_xx, t_z_yzz_0_xy, t_z_yzz_0_xz, t_z_yzz_0_yy,\
                           t_z_yzz_0_yz, t_z_yzz_0_zz, t_z_zzz_0_xx, t_z_zzz_0_xy, t_z_zzz_0_xz,\
                           t_z_zzz_0_yy, t_z_zzz_0_yz, t_z_zzz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_0_zz[i] = t_0_zzzz_0_zz[i] - rab_z[i] * t_0_zzz_0_zz[i];

        t_z_zzz_0_yz[i] = t_0_zzzz_0_yz[i] - rab_z[i] * t_0_zzz_0_yz[i];

        t_z_zzz_0_yy[i] = t_0_zzzz_0_yy[i] - rab_z[i] * t_0_zzz_0_yy[i];

        t_z_zzz_0_xz[i] = t_0_zzzz_0_xz[i] - rab_z[i] * t_0_zzz_0_xz[i];

        t_z_zzz_0_xy[i] = t_0_zzzz_0_xy[i] - rab_z[i] * t_0_zzz_0_xy[i];

        t_z_zzz_0_xx[i] = t_0_zzzz_0_xx[i] - rab_z[i] * t_0_zzz_0_xx[i];

        t_z_yzz_0_zz[i] = t_0_yzzz_0_zz[i] - rab_z[i] * t_0_yzz_0_zz[i];

        t_z_yzz_0_yz[i] = t_0_yzzz_0_yz[i] - rab_z[i] * t_0_yzz_0_yz[i];

        t_z_yzz_0_yy[i] = t_0_yzzz_0_yy[i] - rab_z[i] * t_0_yzz_0_yy[i];

        t_z_yzz_0_xz[i] = t_0_yzzz_0_xz[i] - rab_z[i] * t_0_yzz_0_xz[i];

        t_z_yzz_0_xy[i] = t_0_yzzz_0_xy[i] - rab_z[i] * t_0_yzz_0_xy[i];

        t_z_yzz_0_xx[i] = t_0_yzzz_0_xx[i] - rab_z[i] * t_0_yzz_0_xx[i];

        t_z_yyz_0_zz[i] = t_0_yyzz_0_zz[i] - rab_z[i] * t_0_yyz_0_zz[i];

        t_z_yyz_0_yz[i] = t_0_yyzz_0_yz[i] - rab_z[i] * t_0_yyz_0_yz[i];

        t_z_yyz_0_yy[i] = t_0_yyzz_0_yy[i] - rab_z[i] * t_0_yyz_0_yy[i];

        t_z_yyz_0_xz[i] = t_0_yyzz_0_xz[i] - rab_z[i] * t_0_yyz_0_xz[i];

        t_z_yyz_0_xy[i] = t_0_yyzz_0_xy[i] - rab_z[i] * t_0_yyz_0_xy[i];

        t_z_yyz_0_xx[i] = t_0_yyzz_0_xx[i] - rab_z[i] * t_0_yyz_0_xx[i];

        t_z_xzz_0_zz[i] = t_0_xzzz_0_zz[i] - rab_z[i] * t_0_xzz_0_zz[i];

        t_z_xzz_0_yz[i] = t_0_xzzz_0_yz[i] - rab_z[i] * t_0_xzz_0_yz[i];

        t_z_xzz_0_yy[i] = t_0_xzzz_0_yy[i] - rab_z[i] * t_0_xzz_0_yy[i];

        t_z_xzz_0_xz[i] = t_0_xzzz_0_xz[i] - rab_z[i] * t_0_xzz_0_xz[i];

        t_z_xzz_0_xy[i] = t_0_xzzz_0_xy[i] - rab_z[i] * t_0_xzz_0_xy[i];

        t_z_xzz_0_xx[i] = t_0_xzzz_0_xx[i] - rab_z[i] * t_0_xzz_0_xx[i];

        t_z_xyz_0_zz[i] = t_0_xyzz_0_zz[i] - rab_z[i] * t_0_xyz_0_zz[i];

        t_z_xyz_0_yz[i] = t_0_xyzz_0_yz[i] - rab_z[i] * t_0_xyz_0_yz[i];

        t_z_xyz_0_yy[i] = t_0_xyzz_0_yy[i] - rab_z[i] * t_0_xyz_0_yy[i];

        t_z_xyz_0_xz[i] = t_0_xyzz_0_xz[i] - rab_z[i] * t_0_xyz_0_xz[i];

        t_z_xyz_0_xy[i] = t_0_xyzz_0_xy[i] - rab_z[i] * t_0_xyz_0_xy[i];

        t_z_xyz_0_xx[i] = t_0_xyzz_0_xx[i] - rab_z[i] * t_0_xyz_0_xx[i];

        t_z_xxz_0_zz[i] = t_0_xxzz_0_zz[i] - rab_z[i] * t_0_xxz_0_zz[i];

        t_z_xxz_0_yz[i] = t_0_xxzz_0_yz[i] - rab_z[i] * t_0_xxz_0_yz[i];

        t_z_xxz_0_yy[i] = t_0_xxzz_0_yy[i] - rab_z[i] * t_0_xxz_0_yy[i];

        t_z_xxz_0_xz[i] = t_0_xxzz_0_xz[i] - rab_z[i] * t_0_xxz_0_xz[i];

        t_z_xxz_0_xy[i] = t_0_xxzz_0_xy[i] - rab_z[i] * t_0_xxz_0_xy[i];

        t_z_xxz_0_xx[i] = t_0_xxzz_0_xx[i] - rab_z[i] * t_0_xxz_0_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyyz_0_xx, t_0_xyyz_0_xy, t_0_xyyz_0_xz, t_0_xyyz_0_yy,\
                           t_0_xyyz_0_yz, t_0_xyyz_0_zz, t_0_xyz_0_xx, t_0_xyz_0_xy,\
                           t_0_xyz_0_xz, t_0_xyz_0_yy, t_0_xyz_0_yz, t_0_xyz_0_zz, t_0_xyzz_0_xx,\
                           t_0_xyzz_0_xy, t_0_xyzz_0_xz, t_0_xyzz_0_yy, t_0_xyzz_0_yz,\
                           t_0_xyzz_0_zz, t_0_xzz_0_xx, t_0_xzz_0_xy, t_0_xzz_0_xz,\
                           t_0_xzz_0_yy, t_0_xzz_0_yz, t_0_xzz_0_zz, t_0_yyy_0_xx, t_0_yyy_0_xy,\
                           t_0_yyy_0_xz, t_0_yyy_0_yy, t_0_yyy_0_yz, t_0_yyy_0_zz, t_0_yyyy_0_xx,\
                           t_0_yyyy_0_xy, t_0_yyyy_0_xz, t_0_yyyy_0_yy, t_0_yyyy_0_yz,\
                           t_0_yyyy_0_zz, t_0_yyyz_0_xx, t_0_yyyz_0_xy, t_0_yyyz_0_xz,\
                           t_0_yyyz_0_yy, t_0_yyyz_0_yz, t_0_yyyz_0_zz, t_0_yyz_0_xx,\
                           t_0_yyz_0_xy, t_0_yyz_0_xz, t_0_yyz_0_yy, t_0_yyz_0_yz, t_0_yyz_0_zz,\
                           t_0_yyzz_0_xx, t_0_yyzz_0_xy, t_0_yyzz_0_xz, t_0_yyzz_0_yy,\
                           t_0_yyzz_0_yz, t_0_yyzz_0_zz, t_0_yzz_0_xx, t_0_yzz_0_xy,\
                           t_0_yzz_0_xz, t_0_yzz_0_yy, t_0_yzz_0_yz, t_0_yzz_0_zz, t_0_yzzz_0_xx,\
                           t_0_yzzz_0_xy, t_0_yzzz_0_xz, t_0_yzzz_0_yy, t_0_yzzz_0_yz,\
                           t_0_yzzz_0_zz, t_0_zzz_0_xx, t_0_zzz_0_xy, t_0_zzz_0_xz,\
                           t_0_zzz_0_yy, t_0_zzz_0_yz, t_0_zzz_0_zz, t_y_xyz_0_xx, t_y_xyz_0_xy,\
                           t_y_xyz_0_xz, t_y_xyz_0_yy, t_y_xyz_0_yz, t_y_xyz_0_zz, t_y_xzz_0_xx,\
                           t_y_xzz_0_xy, t_y_xzz_0_xz, t_y_xzz_0_yy, t_y_xzz_0_yz, t_y_xzz_0_zz,\
                           t_y_yyy_0_xx, t_y_yyy_0_xy, t_y_yyy_0_xz, t_y_yyy_0_yy, t_y_yyy_0_yz,\
                           t_y_yyy_0_zz, t_y_yyz_0_xx, t_y_yyz_0_xy, t_y_yyz_0_xz, t_y_yyz_0_yy,\
                           t_y_yyz_0_yz, t_y_yyz_0_zz, t_y_yzz_0_xx, t_y_yzz_0_xy, t_y_yzz_0_xz,\
                           t_y_yzz_0_yy, t_y_yzz_0_yz, t_y_yzz_0_zz, t_y_zzz_0_xx, t_y_zzz_0_xy,\
                           t_y_zzz_0_xz, t_y_zzz_0_yy, t_y_zzz_0_yz, t_y_zzz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zzz_0_zz[i] = t_0_yzzz_0_zz[i] - rab_y[i] * t_0_zzz_0_zz[i];

        t_y_zzz_0_yz[i] = t_0_yzzz_0_yz[i] - rab_y[i] * t_0_zzz_0_yz[i];

        t_y_zzz_0_yy[i] = t_0_yzzz_0_yy[i] - rab_y[i] * t_0_zzz_0_yy[i];

        t_y_zzz_0_xz[i] = t_0_yzzz_0_xz[i] - rab_y[i] * t_0_zzz_0_xz[i];

        t_y_zzz_0_xy[i] = t_0_yzzz_0_xy[i] - rab_y[i] * t_0_zzz_0_xy[i];

        t_y_zzz_0_xx[i] = t_0_yzzz_0_xx[i] - rab_y[i] * t_0_zzz_0_xx[i];

        t_y_yzz_0_zz[i] = t_0_yyzz_0_zz[i] - rab_y[i] * t_0_yzz_0_zz[i];

        t_y_yzz_0_yz[i] = t_0_yyzz_0_yz[i] - rab_y[i] * t_0_yzz_0_yz[i];

        t_y_yzz_0_yy[i] = t_0_yyzz_0_yy[i] - rab_y[i] * t_0_yzz_0_yy[i];

        t_y_yzz_0_xz[i] = t_0_yyzz_0_xz[i] - rab_y[i] * t_0_yzz_0_xz[i];

        t_y_yzz_0_xy[i] = t_0_yyzz_0_xy[i] - rab_y[i] * t_0_yzz_0_xy[i];

        t_y_yzz_0_xx[i] = t_0_yyzz_0_xx[i] - rab_y[i] * t_0_yzz_0_xx[i];

        t_y_yyz_0_zz[i] = t_0_yyyz_0_zz[i] - rab_y[i] * t_0_yyz_0_zz[i];

        t_y_yyz_0_yz[i] = t_0_yyyz_0_yz[i] - rab_y[i] * t_0_yyz_0_yz[i];

        t_y_yyz_0_yy[i] = t_0_yyyz_0_yy[i] - rab_y[i] * t_0_yyz_0_yy[i];

        t_y_yyz_0_xz[i] = t_0_yyyz_0_xz[i] - rab_y[i] * t_0_yyz_0_xz[i];

        t_y_yyz_0_xy[i] = t_0_yyyz_0_xy[i] - rab_y[i] * t_0_yyz_0_xy[i];

        t_y_yyz_0_xx[i] = t_0_yyyz_0_xx[i] - rab_y[i] * t_0_yyz_0_xx[i];

        t_y_yyy_0_zz[i] = t_0_yyyy_0_zz[i] - rab_y[i] * t_0_yyy_0_zz[i];

        t_y_yyy_0_yz[i] = t_0_yyyy_0_yz[i] - rab_y[i] * t_0_yyy_0_yz[i];

        t_y_yyy_0_yy[i] = t_0_yyyy_0_yy[i] - rab_y[i] * t_0_yyy_0_yy[i];

        t_y_yyy_0_xz[i] = t_0_yyyy_0_xz[i] - rab_y[i] * t_0_yyy_0_xz[i];

        t_y_yyy_0_xy[i] = t_0_yyyy_0_xy[i] - rab_y[i] * t_0_yyy_0_xy[i];

        t_y_yyy_0_xx[i] = t_0_yyyy_0_xx[i] - rab_y[i] * t_0_yyy_0_xx[i];

        t_y_xzz_0_zz[i] = t_0_xyzz_0_zz[i] - rab_y[i] * t_0_xzz_0_zz[i];

        t_y_xzz_0_yz[i] = t_0_xyzz_0_yz[i] - rab_y[i] * t_0_xzz_0_yz[i];

        t_y_xzz_0_yy[i] = t_0_xyzz_0_yy[i] - rab_y[i] * t_0_xzz_0_yy[i];

        t_y_xzz_0_xz[i] = t_0_xyzz_0_xz[i] - rab_y[i] * t_0_xzz_0_xz[i];

        t_y_xzz_0_xy[i] = t_0_xyzz_0_xy[i] - rab_y[i] * t_0_xzz_0_xy[i];

        t_y_xzz_0_xx[i] = t_0_xyzz_0_xx[i] - rab_y[i] * t_0_xzz_0_xx[i];

        t_y_xyz_0_zz[i] = t_0_xyyz_0_zz[i] - rab_y[i] * t_0_xyz_0_zz[i];

        t_y_xyz_0_yz[i] = t_0_xyyz_0_yz[i] - rab_y[i] * t_0_xyz_0_yz[i];

        t_y_xyz_0_yy[i] = t_0_xyyz_0_yy[i] - rab_y[i] * t_0_xyz_0_yy[i];

        t_y_xyz_0_xz[i] = t_0_xyyz_0_xz[i] - rab_y[i] * t_0_xyz_0_xz[i];

        t_y_xyz_0_xy[i] = t_0_xyyz_0_xy[i] - rab_y[i] * t_0_xyz_0_xy[i];

        t_y_xyz_0_xx[i] = t_0_xyyz_0_xx[i] - rab_y[i] * t_0_xyz_0_xx[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_0_xxy_0_xx, t_0_xxy_0_xy, t_0_xxy_0_xz, t_0_xxy_0_yy,\
                           t_0_xxy_0_yz, t_0_xxy_0_zz, t_0_xxyy_0_xx, t_0_xxyy_0_xy,\
                           t_0_xxyy_0_xz, t_0_xxyy_0_yy, t_0_xxyy_0_yz, t_0_xxyy_0_zz,\
                           t_0_xxyz_0_xx, t_0_xxyz_0_xy, t_0_xxyz_0_xz, t_0_xxyz_0_yy,\
                           t_0_xxyz_0_yz, t_0_xxyz_0_zz, t_0_xxz_0_xx, t_0_xxz_0_xy,\
                           t_0_xxz_0_xz, t_0_xxz_0_yy, t_0_xxz_0_yz, t_0_xxz_0_zz, t_0_xyy_0_xx,\
                           t_0_xyy_0_xy, t_0_xyy_0_xz, t_0_xyy_0_yy, t_0_xyy_0_yz, t_0_xyy_0_zz,\
                           t_0_xyyy_0_xx, t_0_xyyy_0_xy, t_0_xyyy_0_xz, t_0_xyyy_0_yy,\
                           t_0_xyyy_0_yz, t_0_xyyy_0_zz, t_0_xyyz_0_xx, t_0_xyyz_0_xy,\
                           t_0_xyyz_0_xz, t_0_xyyz_0_yy, t_0_xyyz_0_yz, t_0_xyyz_0_zz,\
                           t_0_xyzz_0_xx, t_0_xyzz_0_xy, t_0_xyzz_0_xz, t_0_xyzz_0_yy,\
                           t_0_xyzz_0_yz, t_0_xyzz_0_zz, t_0_xzzz_0_xx, t_0_xzzz_0_xy,\
                           t_0_xzzz_0_xz, t_0_xzzz_0_yy, t_0_xzzz_0_yz, t_0_xzzz_0_zz,\
                           t_0_yyz_0_xx, t_0_yyz_0_xy, t_0_yyz_0_xz, t_0_yyz_0_yy, t_0_yyz_0_yz,\
                           t_0_yyz_0_zz, t_0_yzz_0_xx, t_0_yzz_0_xy, t_0_yzz_0_xz, t_0_yzz_0_yy,\
                           t_0_yzz_0_yz, t_0_yzz_0_zz, t_0_zzz_0_xx, t_0_zzz_0_xy, t_0_zzz_0_xz,\
                           t_0_zzz_0_yy, t_0_zzz_0_yz, t_0_zzz_0_zz, t_x_yyz_0_xx, t_x_yyz_0_xy,\
                           t_x_yyz_0_xz, t_x_yyz_0_yy, t_x_yyz_0_yz, t_x_yyz_0_zz, t_x_yzz_0_xx,\
                           t_x_yzz_0_xy, t_x_yzz_0_xz, t_x_yzz_0_yy, t_x_yzz_0_yz, t_x_yzz_0_zz,\
                           t_x_zzz_0_xx, t_x_zzz_0_xy, t_x_zzz_0_xz, t_x_zzz_0_yy, t_x_zzz_0_yz,\
                           t_x_zzz_0_zz, t_y_xxy_0_xx, t_y_xxy_0_xy, t_y_xxy_0_xz, t_y_xxy_0_yy,\
                           t_y_xxy_0_yz, t_y_xxy_0_zz, t_y_xxz_0_xx, t_y_xxz_0_xy, t_y_xxz_0_xz,\
                           t_y_xxz_0_yy, t_y_xxz_0_yz, t_y_xxz_0_zz, t_y_xyy_0_xx, t_y_xyy_0_xy,\
                           t_y_xyy_0_xz, t_y_xyy_0_yy, t_y_xyy_0_yz, t_y_xyy_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyy_0_zz[i] = t_0_xyyy_0_zz[i] - rab_y[i] * t_0_xyy_0_zz[i];

        t_y_xyy_0_yz[i] = t_0_xyyy_0_yz[i] - rab_y[i] * t_0_xyy_0_yz[i];

        t_y_xyy_0_yy[i] = t_0_xyyy_0_yy[i] - rab_y[i] * t_0_xyy_0_yy[i];

        t_y_xyy_0_xz[i] = t_0_xyyy_0_xz[i] - rab_y[i] * t_0_xyy_0_xz[i];

        t_y_xyy_0_xy[i] = t_0_xyyy_0_xy[i] - rab_y[i] * t_0_xyy_0_xy[i];

        t_y_xyy_0_xx[i] = t_0_xyyy_0_xx[i] - rab_y[i] * t_0_xyy_0_xx[i];

        t_y_xxz_0_zz[i] = t_0_xxyz_0_zz[i] - rab_y[i] * t_0_xxz_0_zz[i];

        t_y_xxz_0_yz[i] = t_0_xxyz_0_yz[i] - rab_y[i] * t_0_xxz_0_yz[i];

        t_y_xxz_0_yy[i] = t_0_xxyz_0_yy[i] - rab_y[i] * t_0_xxz_0_yy[i];

        t_y_xxz_0_xz[i] = t_0_xxyz_0_xz[i] - rab_y[i] * t_0_xxz_0_xz[i];

        t_y_xxz_0_xy[i] = t_0_xxyz_0_xy[i] - rab_y[i] * t_0_xxz_0_xy[i];

        t_y_xxz_0_xx[i] = t_0_xxyz_0_xx[i] - rab_y[i] * t_0_xxz_0_xx[i];

        t_y_xxy_0_zz[i] = t_0_xxyy_0_zz[i] - rab_y[i] * t_0_xxy_0_zz[i];

        t_y_xxy_0_yz[i] = t_0_xxyy_0_yz[i] - rab_y[i] * t_0_xxy_0_yz[i];

        t_y_xxy_0_yy[i] = t_0_xxyy_0_yy[i] - rab_y[i] * t_0_xxy_0_yy[i];

        t_y_xxy_0_xz[i] = t_0_xxyy_0_xz[i] - rab_y[i] * t_0_xxy_0_xz[i];

        t_y_xxy_0_xy[i] = t_0_xxyy_0_xy[i] - rab_y[i] * t_0_xxy_0_xy[i];

        t_y_xxy_0_xx[i] = t_0_xxyy_0_xx[i] - rab_y[i] * t_0_xxy_0_xx[i];

        t_x_zzz_0_zz[i] = t_0_xzzz_0_zz[i] - rab_x[i] * t_0_zzz_0_zz[i];

        t_x_zzz_0_yz[i] = t_0_xzzz_0_yz[i] - rab_x[i] * t_0_zzz_0_yz[i];

        t_x_zzz_0_yy[i] = t_0_xzzz_0_yy[i] - rab_x[i] * t_0_zzz_0_yy[i];

        t_x_zzz_0_xz[i] = t_0_xzzz_0_xz[i] - rab_x[i] * t_0_zzz_0_xz[i];

        t_x_zzz_0_xy[i] = t_0_xzzz_0_xy[i] - rab_x[i] * t_0_zzz_0_xy[i];

        t_x_zzz_0_xx[i] = t_0_xzzz_0_xx[i] - rab_x[i] * t_0_zzz_0_xx[i];

        t_x_yzz_0_zz[i] = t_0_xyzz_0_zz[i] - rab_x[i] * t_0_yzz_0_zz[i];

        t_x_yzz_0_yz[i] = t_0_xyzz_0_yz[i] - rab_x[i] * t_0_yzz_0_yz[i];

        t_x_yzz_0_yy[i] = t_0_xyzz_0_yy[i] - rab_x[i] * t_0_yzz_0_yy[i];

        t_x_yzz_0_xz[i] = t_0_xyzz_0_xz[i] - rab_x[i] * t_0_yzz_0_xz[i];

        t_x_yzz_0_xy[i] = t_0_xyzz_0_xy[i] - rab_x[i] * t_0_yzz_0_xy[i];

        t_x_yzz_0_xx[i] = t_0_xyzz_0_xx[i] - rab_x[i] * t_0_yzz_0_xx[i];

        t_x_yyz_0_zz[i] = t_0_xyyz_0_zz[i] - rab_x[i] * t_0_yyz_0_zz[i];

        t_x_yyz_0_yz[i] = t_0_xyyz_0_yz[i] - rab_x[i] * t_0_yyz_0_yz[i];

        t_x_yyz_0_yy[i] = t_0_xyyz_0_yy[i] - rab_x[i] * t_0_yyz_0_yy[i];

        t_x_yyz_0_xz[i] = t_0_xyyz_0_xz[i] - rab_x[i] * t_0_yyz_0_xz[i];

        t_x_yyz_0_xy[i] = t_0_xyyz_0_xy[i] - rab_x[i] * t_0_yyz_0_xy[i];

        t_x_yyz_0_xx[i] = t_0_xyyz_0_xx[i] - rab_x[i] * t_0_yyz_0_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxxy_0_xx, t_0_xxxy_0_xy, t_0_xxxy_0_xz, t_0_xxxy_0_yy,\
                           t_0_xxxy_0_yz, t_0_xxxy_0_zz, t_0_xxxz_0_xx, t_0_xxxz_0_xy,\
                           t_0_xxxz_0_xz, t_0_xxxz_0_yy, t_0_xxxz_0_yz, t_0_xxxz_0_zz,\
                           t_0_xxy_0_xx, t_0_xxy_0_xy, t_0_xxy_0_xz, t_0_xxy_0_yy, t_0_xxy_0_yz,\
                           t_0_xxy_0_zz, t_0_xxyy_0_xx, t_0_xxyy_0_xy, t_0_xxyy_0_xz,\
                           t_0_xxyy_0_yy, t_0_xxyy_0_yz, t_0_xxyy_0_zz, t_0_xxyz_0_xx,\
                           t_0_xxyz_0_xy, t_0_xxyz_0_xz, t_0_xxyz_0_yy, t_0_xxyz_0_yz,\
                           t_0_xxyz_0_zz, t_0_xxz_0_xx, t_0_xxz_0_xy, t_0_xxz_0_xz,\
                           t_0_xxz_0_yy, t_0_xxz_0_yz, t_0_xxz_0_zz, t_0_xxzz_0_xx,\
                           t_0_xxzz_0_xy, t_0_xxzz_0_xz, t_0_xxzz_0_yy, t_0_xxzz_0_yz,\
                           t_0_xxzz_0_zz, t_0_xyy_0_xx, t_0_xyy_0_xy, t_0_xyy_0_xz,\
                           t_0_xyy_0_yy, t_0_xyy_0_yz, t_0_xyy_0_zz, t_0_xyyy_0_xx,\
                           t_0_xyyy_0_xy, t_0_xyyy_0_xz, t_0_xyyy_0_yy, t_0_xyyy_0_yz,\
                           t_0_xyyy_0_zz, t_0_xyz_0_xx, t_0_xyz_0_xy, t_0_xyz_0_xz,\
                           t_0_xyz_0_yy, t_0_xyz_0_yz, t_0_xyz_0_zz, t_0_xzz_0_xx, t_0_xzz_0_xy,\
                           t_0_xzz_0_xz, t_0_xzz_0_yy, t_0_xzz_0_yz, t_0_xzz_0_zz, t_0_yyy_0_xx,\
                           t_0_yyy_0_xy, t_0_yyy_0_xz, t_0_yyy_0_yy, t_0_yyy_0_yz, t_0_yyy_0_zz,\
                           t_x_xxy_0_xx, t_x_xxy_0_xy, t_x_xxy_0_xz, t_x_xxy_0_yy, t_x_xxy_0_yz,\
                           t_x_xxy_0_zz, t_x_xxz_0_xx, t_x_xxz_0_xy, t_x_xxz_0_xz, t_x_xxz_0_yy,\
                           t_x_xxz_0_yz, t_x_xxz_0_zz, t_x_xyy_0_xx, t_x_xyy_0_xy, t_x_xyy_0_xz,\
                           t_x_xyy_0_yy, t_x_xyy_0_yz, t_x_xyy_0_zz, t_x_xyz_0_xx, t_x_xyz_0_xy,\
                           t_x_xyz_0_xz, t_x_xyz_0_yy, t_x_xyz_0_yz, t_x_xyz_0_zz, t_x_xzz_0_xx,\
                           t_x_xzz_0_xy, t_x_xzz_0_xz, t_x_xzz_0_yy, t_x_xzz_0_yz, t_x_xzz_0_zz,\
                           t_x_yyy_0_xx, t_x_yyy_0_xy, t_x_yyy_0_xz, t_x_yyy_0_yy, t_x_yyy_0_yz,\
                           t_x_yyy_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yyy_0_zz[i] = t_0_xyyy_0_zz[i] - rab_x[i] * t_0_yyy_0_zz[i];

        t_x_yyy_0_yz[i] = t_0_xyyy_0_yz[i] - rab_x[i] * t_0_yyy_0_yz[i];

        t_x_yyy_0_yy[i] = t_0_xyyy_0_yy[i] - rab_x[i] * t_0_yyy_0_yy[i];

        t_x_yyy_0_xz[i] = t_0_xyyy_0_xz[i] - rab_x[i] * t_0_yyy_0_xz[i];

        t_x_yyy_0_xy[i] = t_0_xyyy_0_xy[i] - rab_x[i] * t_0_yyy_0_xy[i];

        t_x_yyy_0_xx[i] = t_0_xyyy_0_xx[i] - rab_x[i] * t_0_yyy_0_xx[i];

        t_x_xzz_0_zz[i] = t_0_xxzz_0_zz[i] - rab_x[i] * t_0_xzz_0_zz[i];

        t_x_xzz_0_yz[i] = t_0_xxzz_0_yz[i] - rab_x[i] * t_0_xzz_0_yz[i];

        t_x_xzz_0_yy[i] = t_0_xxzz_0_yy[i] - rab_x[i] * t_0_xzz_0_yy[i];

        t_x_xzz_0_xz[i] = t_0_xxzz_0_xz[i] - rab_x[i] * t_0_xzz_0_xz[i];

        t_x_xzz_0_xy[i] = t_0_xxzz_0_xy[i] - rab_x[i] * t_0_xzz_0_xy[i];

        t_x_xzz_0_xx[i] = t_0_xxzz_0_xx[i] - rab_x[i] * t_0_xzz_0_xx[i];

        t_x_xyz_0_zz[i] = t_0_xxyz_0_zz[i] - rab_x[i] * t_0_xyz_0_zz[i];

        t_x_xyz_0_yz[i] = t_0_xxyz_0_yz[i] - rab_x[i] * t_0_xyz_0_yz[i];

        t_x_xyz_0_yy[i] = t_0_xxyz_0_yy[i] - rab_x[i] * t_0_xyz_0_yy[i];

        t_x_xyz_0_xz[i] = t_0_xxyz_0_xz[i] - rab_x[i] * t_0_xyz_0_xz[i];

        t_x_xyz_0_xy[i] = t_0_xxyz_0_xy[i] - rab_x[i] * t_0_xyz_0_xy[i];

        t_x_xyz_0_xx[i] = t_0_xxyz_0_xx[i] - rab_x[i] * t_0_xyz_0_xx[i];

        t_x_xyy_0_zz[i] = t_0_xxyy_0_zz[i] - rab_x[i] * t_0_xyy_0_zz[i];

        t_x_xyy_0_yz[i] = t_0_xxyy_0_yz[i] - rab_x[i] * t_0_xyy_0_yz[i];

        t_x_xyy_0_yy[i] = t_0_xxyy_0_yy[i] - rab_x[i] * t_0_xyy_0_yy[i];

        t_x_xyy_0_xz[i] = t_0_xxyy_0_xz[i] - rab_x[i] * t_0_xyy_0_xz[i];

        t_x_xyy_0_xy[i] = t_0_xxyy_0_xy[i] - rab_x[i] * t_0_xyy_0_xy[i];

        t_x_xyy_0_xx[i] = t_0_xxyy_0_xx[i] - rab_x[i] * t_0_xyy_0_xx[i];

        t_x_xxz_0_zz[i] = t_0_xxxz_0_zz[i] - rab_x[i] * t_0_xxz_0_zz[i];

        t_x_xxz_0_yz[i] = t_0_xxxz_0_yz[i] - rab_x[i] * t_0_xxz_0_yz[i];

        t_x_xxz_0_yy[i] = t_0_xxxz_0_yy[i] - rab_x[i] * t_0_xxz_0_yy[i];

        t_x_xxz_0_xz[i] = t_0_xxxz_0_xz[i] - rab_x[i] * t_0_xxz_0_xz[i];

        t_x_xxz_0_xy[i] = t_0_xxxz_0_xy[i] - rab_x[i] * t_0_xxz_0_xy[i];

        t_x_xxz_0_xx[i] = t_0_xxxz_0_xx[i] - rab_x[i] * t_0_xxz_0_xx[i];

        t_x_xxy_0_zz[i] = t_0_xxxy_0_zz[i] - rab_x[i] * t_0_xxy_0_zz[i];

        t_x_xxy_0_yz[i] = t_0_xxxy_0_yz[i] - rab_x[i] * t_0_xxy_0_yz[i];

        t_x_xxy_0_yy[i] = t_0_xxxy_0_yy[i] - rab_x[i] * t_0_xxy_0_yy[i];

        t_x_xxy_0_xz[i] = t_0_xxxy_0_xz[i] - rab_x[i] * t_0_xxy_0_xz[i];

        t_x_xxy_0_xy[i] = t_0_xxxy_0_xy[i] - rab_x[i] * t_0_xxy_0_xy[i];

        t_x_xxy_0_xx[i] = t_0_xxxy_0_xx[i] - rab_x[i] * t_0_xxy_0_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxx_0_xx, t_0_xxx_0_xy, t_0_xxx_0_xz, t_0_xxx_0_yy,\
                           t_0_xxx_0_yz, t_0_xxx_0_zz, t_0_xxxx_0_xx, t_0_xxxx_0_xy,\
                           t_0_xxxx_0_xz, t_0_xxxx_0_yy, t_0_xxxx_0_yz, t_0_xxxx_0_zz,\
                           t_x_xxx_0_xx, t_x_xxx_0_xy, t_x_xxx_0_xz, t_x_xxx_0_yy, t_x_xxx_0_yz,\
                           t_x_xxx_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxx_0_zz[i] = t_0_xxxx_0_zz[i] - rab_x[i] * t_0_xxx_0_zz[i];

        t_x_xxx_0_yz[i] = t_0_xxxx_0_yz[i] - rab_x[i] * t_0_xxx_0_yz[i];

        t_x_xxx_0_yy[i] = t_0_xxxx_0_yy[i] - rab_x[i] * t_0_xxx_0_yy[i];

        t_x_xxx_0_xz[i] = t_0_xxxx_0_xz[i] - rab_x[i] * t_0_xxx_0_xz[i];

        t_x_xxx_0_xy[i] = t_0_xxxx_0_xy[i] - rab_x[i] * t_0_xxx_0_xy[i];

        t_x_xxx_0_xx[i] = t_0_xxxx_0_xx[i] - rab_x[i] * t_0_xxx_0_xx[i];
    }
}


} // derirec namespace
