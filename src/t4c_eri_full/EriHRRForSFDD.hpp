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

    t_0_zzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(1));

    t_0_zzz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(2));

    t_0_zzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(3));

    t_0_zzz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(4));

    t_0_zzz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(5));

    t_0_zzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(6));

    t_0_zzz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(7));

    t_0_zzz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(8));

    t_0_zzz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(9));

    t_0_zzz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(10));

    t_0_zzz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(11));

    t_0_zzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(12));

    t_0_zzz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(13));

    t_0_zzz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(14));

    t_0_zzz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(15));

    t_0_zzz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(16));

    t_0_zzz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(17));

    t_0_zzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(18));

    t_0_zzz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(19));

    t_0_zzz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(20));

    t_0_zzz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(21));

    t_0_zzz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(22));

    t_0_zzz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(23));

    t_0_zzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(24));

    t_0_zzz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(25));

    t_0_zzz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(26));

    t_0_zzz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(27));

    t_0_zzz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(28));

    t_0_zzz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(29));

    t_0_zzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(30));

    t_0_zzz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(31));

    t_0_zzz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(32));

    t_0_zzz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(33));

    t_0_zzz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(34));

    t_0_zzz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(35));

    t_0_yzz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(36));

    t_0_yzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(37));

    t_0_yzz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(38));

    t_0_yzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(39));

    t_0_yzz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(40));

    t_0_yzz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(41));

    t_0_yzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(42));

    t_0_yzz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(43));

    t_0_yzz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(44));

    t_0_yzz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(45));

    t_0_yzz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(46));

    t_0_yzz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(47));

    t_0_yzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(48));

    t_0_yzz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(49));

    t_0_yzz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(50));

    t_0_yzz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(51));

    t_0_yzz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(52));

    t_0_yzz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(53));

    t_0_yzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(54));

    t_0_yzz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(55));

    t_0_yzz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(56));

    t_0_yzz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(57));

    t_0_yzz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(58));

    t_0_yzz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(59));

    t_0_yzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(60));

    t_0_yzz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(61));

    t_0_yzz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(62));

    t_0_yzz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(63));

    t_0_yzz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(64));

    t_0_yzz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(65));

    t_0_yzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(66));

    t_0_yzz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(67));

    t_0_yzz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(68));

    t_0_yzz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(69));

    t_0_yzz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(70));

    t_0_yzz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(71));

    t_0_yyz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(72));

    t_0_yyz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(73));

    t_0_yyz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(74));

    t_0_yyz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(75));

    t_0_yyz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(76));

    t_0_yyz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(77));

    t_0_yyz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(78));

    t_0_yyz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(79));

    t_0_yyz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(80));

    t_0_yyz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(81));

    t_0_yyz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(82));

    t_0_yyz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(83));

    t_0_yyz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(84));

    t_0_yyz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(85));

    t_0_yyz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(86));

    t_0_yyz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(87));

    t_0_yyz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(88));

    t_0_yyz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(89));

    t_0_yyz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(90));

    t_0_yyz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(91));

    t_0_yyz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(92));

    t_0_yyz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(93));

    t_0_yyz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(94));

    t_0_yyz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(95));

    t_0_yyz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(96));

    t_0_yyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(97));

    t_0_yyz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(98));

    t_0_yyz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(99));

    t_0_yyz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(100));

    t_0_yyz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(101));

    t_0_yyz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(102));

    t_0_yyz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(103));

    t_0_yyz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(104));

    t_0_yyz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(105));

    t_0_yyz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(106));

    t_0_yyz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(107));

    t_0_yyy_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(108));

    t_0_yyy_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(109));

    t_0_yyy_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(110));

    t_0_yyy_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(111));

    t_0_yyy_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(112));

    t_0_yyy_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(113));

    t_0_yyy_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(114));

    t_0_yyy_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(115));

    t_0_yyy_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(116));

    t_0_yyy_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(117));

    t_0_yyy_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(118));

    t_0_yyy_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(119));

    t_0_yyy_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(120));

    t_0_yyy_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(121));

    t_0_yyy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(122));

    t_0_yyy_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(123));

    t_0_yyy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(124));

    t_0_yyy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(125));

    t_0_yyy_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(126));

    t_0_yyy_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(127));

    t_0_yyy_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(128));

    t_0_yyy_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(129));

    t_0_yyy_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(130));

    t_0_yyy_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(131));

    t_0_yyy_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(132));

    t_0_yyy_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(133));

    t_0_yyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(134));

    t_0_yyy_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(135));

    t_0_yyy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(136));

    t_0_yyy_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(137));

    t_0_yyy_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(138));

    t_0_yyy_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(139));

    t_0_yyy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(140));

    t_0_yyy_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(141));

    t_0_yyy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(142));

    t_0_yyy_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(143));

    t_0_xzz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(144));

    t_0_xzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(145));

    t_0_xzz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(146));

    t_0_xzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(147));

    t_0_xzz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(148));

    t_0_xzz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(149));

    t_0_xzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(150));

    t_0_xzz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(151));

    t_0_xzz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(152));

    t_0_xzz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(153));

    t_0_xzz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(154));

    t_0_xzz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(155));

    t_0_xzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(156));

    t_0_xzz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(157));

    t_0_xzz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(158));

    t_0_xzz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(159));

    t_0_xzz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(160));

    t_0_xzz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(161));

    t_0_xzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(162));

    t_0_xzz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(163));

    t_0_xzz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(164));

    t_0_xzz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(165));

    t_0_xzz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(166));

    t_0_xzz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(167));

    t_0_xzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(168));

    t_0_xzz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(169));

    t_0_xzz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(170));

    t_0_xzz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(171));

    t_0_xzz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(172));

    t_0_xzz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(173));

    t_0_xzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(174));

    t_0_xzz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(175));

    t_0_xzz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(176));

    t_0_xzz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(177));

    t_0_xzz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(178));

    t_0_xzz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(179));

    t_0_xyz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(180));

    t_0_xyz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(181));

    t_0_xyz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(182));

    t_0_xyz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(183));

    t_0_xyz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(184));

    t_0_xyz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(185));

    t_0_xyz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(186));

    t_0_xyz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(187));

    t_0_xyz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(188));

    t_0_xyz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(189));

    t_0_xyz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(190));

    t_0_xyz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(191));

    t_0_xyz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(192));

    t_0_xyz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(193));

    t_0_xyz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(194));

    t_0_xyz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(195));

    t_0_xyz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(196));

    t_0_xyz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(197));

    t_0_xyz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(198));

    t_0_xyz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(199));

    t_0_xyz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(200));

    t_0_xyz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(201));

    t_0_xyz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(202));

    t_0_xyz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(203));

    t_0_xyz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(204));

    t_0_xyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(205));

    t_0_xyz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(206));

    t_0_xyz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(207));

    t_0_xyz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(208));

    t_0_xyz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(209));

    t_0_xyz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(210));

    t_0_xyz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(211));

    t_0_xyz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(212));

    t_0_xyz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(213));

    t_0_xyz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(214));

    t_0_xyz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(215));

    t_0_xyy_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(216));

    t_0_xyy_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(217));

    t_0_xyy_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(218));

    t_0_xyy_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(219));

    t_0_xyy_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(220));

    t_0_xyy_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(221));

    t_0_xyy_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(222));

    t_0_xyy_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(223));

    t_0_xyy_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(224));

    t_0_xyy_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(225));

    t_0_xyy_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(226));

    t_0_xyy_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(227));

    t_0_xyy_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(228));

    t_0_xyy_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(229));

    t_0_xyy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(230));

    t_0_xyy_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(231));

    t_0_xyy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(232));

    t_0_xyy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(233));

    t_0_xyy_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(234));

    t_0_xyy_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(235));

    t_0_xyy_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(236));

    t_0_xyy_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(237));

    t_0_xyy_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(238));

    t_0_xyy_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(239));

    t_0_xyy_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(240));

    t_0_xyy_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(241));

    t_0_xyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(242));

    t_0_xyy_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(243));

    t_0_xyy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(244));

    t_0_xyy_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(245));

    t_0_xyy_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(246));

    t_0_xyy_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(247));

    t_0_xyy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(248));

    t_0_xyy_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(249));

    t_0_xyy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(250));

    t_0_xyy_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(251));

    t_0_xxz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(252));

    t_0_xxz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(253));

    t_0_xxz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(254));

    t_0_xxz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(255));

    t_0_xxz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(256));

    t_0_xxz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(257));

    t_0_xxz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(258));

    t_0_xxz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(259));

    t_0_xxz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(260));

    t_0_xxz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(261));

    t_0_xxz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(262));

    t_0_xxz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(263));

    t_0_xxz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(264));

    t_0_xxz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(265));

    t_0_xxz_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(266));

    t_0_xxz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(267));

    t_0_xxz_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(268));

    t_0_xxz_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(269));

    t_0_xxz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(270));

    t_0_xxz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(271));

    t_0_xxz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(272));

    t_0_xxz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(273));

    t_0_xxz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(274));

    t_0_xxz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(275));

    t_0_xxz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(276));

    t_0_xxz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(277));

    t_0_xxz_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(278));

    t_0_xxz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(279));

    t_0_xxz_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(280));

    t_0_xxz_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(281));

    t_0_xxz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(282));

    t_0_xxz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(283));

    t_0_xxz_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(284));

    t_0_xxz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(285));

    t_0_xxz_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(286));

    t_0_xxz_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(287));

    t_0_xxy_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(288));

    t_0_xxy_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(289));

    t_0_xxy_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(290));

    t_0_xxy_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(291));

    t_0_xxy_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(292));

    t_0_xxy_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(293));

    t_0_xxy_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(294));

    t_0_xxy_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(295));

    t_0_xxy_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(296));

    t_0_xxy_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(297));

    t_0_xxy_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(298));

    t_0_xxy_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(299));

    t_0_xxy_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(300));

    t_0_xxy_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(301));

    t_0_xxy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(302));

    t_0_xxy_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(303));

    t_0_xxy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(304));

    t_0_xxy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(305));

    t_0_xxy_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(306));

    t_0_xxy_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(307));

    t_0_xxy_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(308));

    t_0_xxy_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(309));

    t_0_xxy_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(310));

    t_0_xxy_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(311));

    t_0_xxy_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(312));

    t_0_xxy_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(313));

    t_0_xxy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(314));

    t_0_xxy_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(315));

    t_0_xxy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(316));

    t_0_xxy_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(317));

    t_0_xxy_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(318));

    t_0_xxy_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(319));

    t_0_xxy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(320));

    t_0_xxy_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(321));

    t_0_xxy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(322));

    t_0_xxy_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(323));

    t_0_xxx_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(324));

    t_0_xxx_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(325));

    t_0_xxx_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(326));

    t_0_xxx_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(327));

    t_0_xxx_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(328));

    t_0_xxx_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(329));

    t_0_xxx_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(330));

    t_0_xxx_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(331));

    t_0_xxx_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(332));

    t_0_xxx_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(333));

    t_0_xxx_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(334));

    t_0_xxx_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(335));

    t_0_xxx_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(336));

    t_0_xxx_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(337));

    t_0_xxx_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(338));

    t_0_xxx_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(339));

    t_0_xxx_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(340));

    t_0_xxx_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(341));

    t_0_xxx_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(342));

    t_0_xxx_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(343));

    t_0_xxx_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(344));

    t_0_xxx_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(345));

    t_0_xxx_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(346));

    t_0_xxx_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(347));

    t_0_xxx_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(348));

    t_0_xxx_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(349));

    t_0_xxx_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(350));

    t_0_xxx_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(351));

    t_0_xxx_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(352));

    t_0_xxx_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(353));

    t_0_xxx_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(354));

    t_0_xxx_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(355));

    t_0_xxx_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(356));

    t_0_xxx_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(357));

    t_0_xxx_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(358));

    t_0_xxx_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(359));

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

    // set up (SFPF) integral components

    t_0_zzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(0));

    t_0_zzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(1));

    t_0_zzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(2));

    t_0_zzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(3));

    t_0_zzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(4));

    t_0_zzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(5));

    t_0_zzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(6));

    t_0_zzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(7));

    t_0_zzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(8));

    t_0_zzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(9));

    t_0_zzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(10));

    t_0_zzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(11));

    t_0_zzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(12));

    t_0_zzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(13));

    t_0_zzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(14));

    t_0_zzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(15));

    t_0_zzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(16));

    t_0_zzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(17));

    t_0_zzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(18));

    t_0_zzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(19));

    t_0_zzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(20));

    t_0_zzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(21));

    t_0_zzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(22));

    t_0_zzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(23));

    t_0_zzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(24));

    t_0_yzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(25));

    t_0_yzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(26));

    t_0_yzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(27));

    t_0_yzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(28));

    t_0_yzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(29));

    t_0_yzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(30));

    t_0_yzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(31));

    t_0_yzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(32));

    t_0_yzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(33));

    t_0_yzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(34));

    t_0_yzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(35));

    t_0_yzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(36));

    t_0_yzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(37));

    t_0_yzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(38));

    t_0_yzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(39));

    t_0_yzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(40));

    t_0_yzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(41));

    t_0_yzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(42));

    t_0_yzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(43));

    t_0_yzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(44));

    t_0_yzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(45));

    t_0_yzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(46));

    t_0_yzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(47));

    t_0_yzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(48));

    t_0_yzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(49));

    t_0_yyz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(50));

    t_0_yyz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(51));

    t_0_yyz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(52));

    t_0_yyz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(53));

    t_0_yyz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(54));

    t_0_yyz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(55));

    t_0_yyz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(56));

    t_0_yyz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(57));

    t_0_yyz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(58));

    t_0_yyz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(59));

    t_0_yyz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(60));

    t_0_yyz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(61));

    t_0_yyz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(62));

    t_0_yyz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(63));

    t_0_yyz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(64));

    t_0_yyz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(65));

    t_0_yyz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(66));

    t_0_yyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(67));

    t_0_yyz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(68));

    t_0_yyz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(69));

    t_0_yyz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(70));

    t_0_yyz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(71));

    t_0_yyz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(72));

    t_0_yyz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(73));

    t_0_yyz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(74));

    t_0_yyy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(75));

    t_0_yyy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(76));

    t_0_yyy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(77));

    t_0_yyy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(78));

    t_0_yyy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(79));

    t_0_yyy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(80));

    t_0_yyy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(81));

    t_0_yyy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(82));

    t_0_yyy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(83));

    t_0_yyy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(84));

    t_0_yyy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(85));

    t_0_yyy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(86));

    t_0_yyy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(87));

    t_0_yyy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(88));

    t_0_yyy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(89));

    t_0_yyy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(90));

    t_0_yyy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(91));

    t_0_yyy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(92));

    t_0_yyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(93));

    t_0_yyy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(94));

    t_0_yyy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(95));

    t_0_yyy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(96));

    t_0_yyy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(97));

    t_0_yyy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(98));

    t_0_yyy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(99));

    t_0_xzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(100));

    t_0_xzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(101));

    t_0_xzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(102));

    t_0_xzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(103));

    t_0_xzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(104));

    t_0_xzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(105));

    t_0_xzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(106));

    t_0_xzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(107));

    t_0_xzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(108));

    t_0_xzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(109));

    t_0_xzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(110));

    t_0_xzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(111));

    t_0_xzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(112));

    t_0_xzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(113));

    t_0_xzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(114));

    t_0_xzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(115));

    t_0_xzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(116));

    t_0_xzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(117));

    t_0_xzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(118));

    t_0_xzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(119));

    t_0_xzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(120));

    t_0_xzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(121));

    t_0_xzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(122));

    t_0_xzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(123));

    t_0_xzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(124));

    t_0_xyz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(125));

    t_0_xyz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(126));

    t_0_xyz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(127));

    t_0_xyz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(128));

    t_0_xyz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(129));

    t_0_xyz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(130));

    t_0_xyz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(131));

    t_0_xyz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(132));

    t_0_xyz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(133));

    t_0_xyz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(134));

    t_0_xyz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(135));

    t_0_xyz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(136));

    t_0_xyz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(137));

    t_0_xyz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(138));

    t_0_xyz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(139));

    t_0_xyz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(140));

    t_0_xyz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(141));

    t_0_xyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(142));

    t_0_xyz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(143));

    t_0_xyz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(144));

    t_0_xyz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(145));

    t_0_xyz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(146));

    t_0_xyz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(147));

    t_0_xyz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(148));

    t_0_xyz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(149));

    t_0_xyy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(150));

    t_0_xyy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(151));

    t_0_xyy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(152));

    t_0_xyy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(153));

    t_0_xyy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(154));

    t_0_xyy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(155));

    t_0_xyy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(156));

    t_0_xyy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(157));

    t_0_xyy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(158));

    t_0_xyy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(159));

    t_0_xyy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(160));

    t_0_xyy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(161));

    t_0_xyy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(162));

    t_0_xyy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(163));

    t_0_xyy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(164));

    t_0_xyy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(165));

    t_0_xyy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(166));

    t_0_xyy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(167));

    t_0_xyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(168));

    t_0_xyy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(169));

    t_0_xyy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(170));

    t_0_xyy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(171));

    t_0_xyy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(172));

    t_0_xyy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(173));

    t_0_xyy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(174));

    t_0_xxz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(175));

    t_0_xxz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(176));

    t_0_xxz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(177));

    t_0_xxz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(178));

    t_0_xxz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(179));

    t_0_xxz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(180));

    t_0_xxz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(181));

    t_0_xxz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(182));

    t_0_xxz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(183));

    t_0_xxz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(184));

    t_0_xxz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(185));

    t_0_xxz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(186));

    t_0_xxz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(187));

    t_0_xxz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(188));

    t_0_xxz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(189));

    t_0_xxz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(190));

    t_0_xxz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(191));

    t_0_xxz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(192));

    t_0_xxz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(193));

    t_0_xxz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(194));

    t_0_xxz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(195));

    t_0_xxz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(196));

    t_0_xxz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(197));

    t_0_xxz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(198));

    t_0_xxz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(199));

    t_0_xxy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(200));

    t_0_xxy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(201));

    t_0_xxy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(202));

    t_0_xxy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(203));

    t_0_xxy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(204));

    t_0_xxy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(205));

    t_0_xxy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(206));

    t_0_xxy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(207));

    t_0_xxy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(208));

    t_0_xxy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(209));

    t_0_xxy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(210));

    t_0_xxy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(211));

    t_0_xxy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(212));

    t_0_xxy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(213));

    t_0_xxy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(214));

    t_0_xxy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(215));

    t_0_xxy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(216));

    t_0_xxy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(217));

    t_0_xxy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(218));

    t_0_xxy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(219));

    t_0_xxy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(220));

    t_0_xxy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(221));

    t_0_xxy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(222));

    t_0_xxy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(223));

    t_0_xxy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(224));

    t_0_xxx_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(225));

    t_0_xxx_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(226));

    t_0_xxx_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(227));

    t_0_xxx_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(228));

    t_0_xxx_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(229));

    t_0_xxx_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(230));

    t_0_xxx_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(231));

    t_0_xxx_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(232));

    t_0_xxx_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(233));

    t_0_xxx_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(234));

    t_0_xxx_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(235));

    t_0_xxx_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(236));

    t_0_xxx_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(237));

    t_0_xxx_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(238));

    t_0_xxx_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(239));

    t_0_xxx_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(240));

    t_0_xxx_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(241));

    t_0_xxx_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(242));

    t_0_xxx_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(243));

    t_0_xxx_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(244));

    t_0_xxx_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(245));

    t_0_xxx_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(246));

    t_0_xxx_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(247));

    t_0_xxx_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(248));

    t_0_xxx_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(249));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_zzz_x_xx, t_0_zzz_x_xxx, t_0_zzz_x_xxy,\
                           t_0_zzz_x_xxz, t_0_zzz_x_xy, t_0_zzz_x_xyy, t_0_zzz_x_xyz,\
                           t_0_zzz_x_xz, t_0_zzz_x_xzz, t_0_zzz_x_yy, t_0_zzz_x_yyy,\
                           t_0_zzz_x_yyz, t_0_zzz_x_yz, t_0_zzz_x_yzz, t_0_zzz_x_zz,\
                           t_0_zzz_x_zzz, t_0_zzz_xx_xx, t_0_zzz_xx_xy, t_0_zzz_xx_xz,\
                           t_0_zzz_xx_yy, t_0_zzz_xx_yz, t_0_zzz_xx_zz, t_0_zzz_xy_xx,\
                           t_0_zzz_xy_xy, t_0_zzz_xy_xz, t_0_zzz_xy_yy, t_0_zzz_xy_yz,\
                           t_0_zzz_xy_zz, t_0_zzz_xz_xx, t_0_zzz_xz_xy, t_0_zzz_xz_xz,\
                           t_0_zzz_xz_yy, t_0_zzz_xz_yz, t_0_zzz_xz_zz, t_0_zzz_y_xx,\
                           t_0_zzz_y_xxy, t_0_zzz_y_xxz, t_0_zzz_y_xy, t_0_zzz_y_xyy,\
                           t_0_zzz_y_xyz, t_0_zzz_y_xz, t_0_zzz_y_xzz, t_0_zzz_y_yy,\
                           t_0_zzz_y_yyy, t_0_zzz_y_yyz, t_0_zzz_y_yz, t_0_zzz_y_yzz,\
                           t_0_zzz_y_zz, t_0_zzz_y_zzz, t_0_zzz_yy_xx, t_0_zzz_yy_xy,\
                           t_0_zzz_yy_xz, t_0_zzz_yy_yy, t_0_zzz_yy_yz, t_0_zzz_yy_zz,\
                           t_0_zzz_yz_xx, t_0_zzz_yz_xy, t_0_zzz_yz_xz, t_0_zzz_yz_yy,\
                           t_0_zzz_yz_yz, t_0_zzz_yz_zz, t_0_zzz_z_xx, t_0_zzz_z_xxz,\
                           t_0_zzz_z_xy, t_0_zzz_z_xyz, t_0_zzz_z_xz, t_0_zzz_z_xzz,\
                           t_0_zzz_z_yy, t_0_zzz_z_yyz, t_0_zzz_z_yz, t_0_zzz_z_yzz,\
                           t_0_zzz_z_zz, t_0_zzz_z_zzz, t_0_zzz_zz_xx, t_0_zzz_zz_xy,\
                           t_0_zzz_zz_xz, t_0_zzz_zz_yy, t_0_zzz_zz_yz, t_0_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_zz_zz[i] = t_0_zzz_z_zzz[i] - rcd_z[i] * t_0_zzz_z_zz[i];

        t_0_zzz_zz_yz[i] = t_0_zzz_z_yzz[i] - rcd_z[i] * t_0_zzz_z_yz[i];

        t_0_zzz_zz_yy[i] = t_0_zzz_z_yyz[i] - rcd_z[i] * t_0_zzz_z_yy[i];

        t_0_zzz_zz_xz[i] = t_0_zzz_z_xzz[i] - rcd_z[i] * t_0_zzz_z_xz[i];

        t_0_zzz_zz_xy[i] = t_0_zzz_z_xyz[i] - rcd_z[i] * t_0_zzz_z_xy[i];

        t_0_zzz_zz_xx[i] = t_0_zzz_z_xxz[i] - rcd_z[i] * t_0_zzz_z_xx[i];

        t_0_zzz_yz_zz[i] = t_0_zzz_y_zzz[i] - rcd_z[i] * t_0_zzz_y_zz[i];

        t_0_zzz_yz_yz[i] = t_0_zzz_y_yzz[i] - rcd_z[i] * t_0_zzz_y_yz[i];

        t_0_zzz_yz_yy[i] = t_0_zzz_y_yyz[i] - rcd_z[i] * t_0_zzz_y_yy[i];

        t_0_zzz_yz_xz[i] = t_0_zzz_y_xzz[i] - rcd_z[i] * t_0_zzz_y_xz[i];

        t_0_zzz_yz_xy[i] = t_0_zzz_y_xyz[i] - rcd_z[i] * t_0_zzz_y_xy[i];

        t_0_zzz_yz_xx[i] = t_0_zzz_y_xxz[i] - rcd_z[i] * t_0_zzz_y_xx[i];

        t_0_zzz_yy_zz[i] = t_0_zzz_y_yzz[i] - rcd_y[i] * t_0_zzz_y_zz[i];

        t_0_zzz_yy_yz[i] = t_0_zzz_y_yyz[i] - rcd_y[i] * t_0_zzz_y_yz[i];

        t_0_zzz_yy_yy[i] = t_0_zzz_y_yyy[i] - rcd_y[i] * t_0_zzz_y_yy[i];

        t_0_zzz_yy_xz[i] = t_0_zzz_y_xyz[i] - rcd_y[i] * t_0_zzz_y_xz[i];

        t_0_zzz_yy_xy[i] = t_0_zzz_y_xyy[i] - rcd_y[i] * t_0_zzz_y_xy[i];

        t_0_zzz_yy_xx[i] = t_0_zzz_y_xxy[i] - rcd_y[i] * t_0_zzz_y_xx[i];

        t_0_zzz_xz_zz[i] = t_0_zzz_x_zzz[i] - rcd_z[i] * t_0_zzz_x_zz[i];

        t_0_zzz_xz_yz[i] = t_0_zzz_x_yzz[i] - rcd_z[i] * t_0_zzz_x_yz[i];

        t_0_zzz_xz_yy[i] = t_0_zzz_x_yyz[i] - rcd_z[i] * t_0_zzz_x_yy[i];

        t_0_zzz_xz_xz[i] = t_0_zzz_x_xzz[i] - rcd_z[i] * t_0_zzz_x_xz[i];

        t_0_zzz_xz_xy[i] = t_0_zzz_x_xyz[i] - rcd_z[i] * t_0_zzz_x_xy[i];

        t_0_zzz_xz_xx[i] = t_0_zzz_x_xxz[i] - rcd_z[i] * t_0_zzz_x_xx[i];

        t_0_zzz_xy_zz[i] = t_0_zzz_x_yzz[i] - rcd_y[i] * t_0_zzz_x_zz[i];

        t_0_zzz_xy_yz[i] = t_0_zzz_x_yyz[i] - rcd_y[i] * t_0_zzz_x_yz[i];

        t_0_zzz_xy_yy[i] = t_0_zzz_x_yyy[i] - rcd_y[i] * t_0_zzz_x_yy[i];

        t_0_zzz_xy_xz[i] = t_0_zzz_x_xyz[i] - rcd_y[i] * t_0_zzz_x_xz[i];

        t_0_zzz_xy_xy[i] = t_0_zzz_x_xyy[i] - rcd_y[i] * t_0_zzz_x_xy[i];

        t_0_zzz_xy_xx[i] = t_0_zzz_x_xxy[i] - rcd_y[i] * t_0_zzz_x_xx[i];

        t_0_zzz_xx_zz[i] = t_0_zzz_x_xzz[i] - rcd_x[i] * t_0_zzz_x_zz[i];

        t_0_zzz_xx_yz[i] = t_0_zzz_x_xyz[i] - rcd_x[i] * t_0_zzz_x_yz[i];

        t_0_zzz_xx_yy[i] = t_0_zzz_x_xyy[i] - rcd_x[i] * t_0_zzz_x_yy[i];

        t_0_zzz_xx_xz[i] = t_0_zzz_x_xxz[i] - rcd_x[i] * t_0_zzz_x_xz[i];

        t_0_zzz_xx_xy[i] = t_0_zzz_x_xxy[i] - rcd_x[i] * t_0_zzz_x_xy[i];

        t_0_zzz_xx_xx[i] = t_0_zzz_x_xxx[i] - rcd_x[i] * t_0_zzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzz_x_xx, t_0_yzz_x_xxx, t_0_yzz_x_xxy,\
                           t_0_yzz_x_xxz, t_0_yzz_x_xy, t_0_yzz_x_xyy, t_0_yzz_x_xyz,\
                           t_0_yzz_x_xz, t_0_yzz_x_xzz, t_0_yzz_x_yy, t_0_yzz_x_yyy,\
                           t_0_yzz_x_yyz, t_0_yzz_x_yz, t_0_yzz_x_yzz, t_0_yzz_x_zz,\
                           t_0_yzz_x_zzz, t_0_yzz_xx_xx, t_0_yzz_xx_xy, t_0_yzz_xx_xz,\
                           t_0_yzz_xx_yy, t_0_yzz_xx_yz, t_0_yzz_xx_zz, t_0_yzz_xy_xx,\
                           t_0_yzz_xy_xy, t_0_yzz_xy_xz, t_0_yzz_xy_yy, t_0_yzz_xy_yz,\
                           t_0_yzz_xy_zz, t_0_yzz_xz_xx, t_0_yzz_xz_xy, t_0_yzz_xz_xz,\
                           t_0_yzz_xz_yy, t_0_yzz_xz_yz, t_0_yzz_xz_zz, t_0_yzz_y_xx,\
                           t_0_yzz_y_xxy, t_0_yzz_y_xxz, t_0_yzz_y_xy, t_0_yzz_y_xyy,\
                           t_0_yzz_y_xyz, t_0_yzz_y_xz, t_0_yzz_y_xzz, t_0_yzz_y_yy,\
                           t_0_yzz_y_yyy, t_0_yzz_y_yyz, t_0_yzz_y_yz, t_0_yzz_y_yzz,\
                           t_0_yzz_y_zz, t_0_yzz_y_zzz, t_0_yzz_yy_xx, t_0_yzz_yy_xy,\
                           t_0_yzz_yy_xz, t_0_yzz_yy_yy, t_0_yzz_yy_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_xx, t_0_yzz_yz_xy, t_0_yzz_yz_xz, t_0_yzz_yz_yy,\
                           t_0_yzz_yz_yz, t_0_yzz_yz_zz, t_0_yzz_z_xx, t_0_yzz_z_xxz,\
                           t_0_yzz_z_xy, t_0_yzz_z_xyz, t_0_yzz_z_xz, t_0_yzz_z_xzz,\
                           t_0_yzz_z_yy, t_0_yzz_z_yyz, t_0_yzz_z_yz, t_0_yzz_z_yzz,\
                           t_0_yzz_z_zz, t_0_yzz_z_zzz, t_0_yzz_zz_xx, t_0_yzz_zz_xy,\
                           t_0_yzz_zz_xz, t_0_yzz_zz_yy, t_0_yzz_zz_yz, t_0_yzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yzz_zz_zz[i] = t_0_yzz_z_zzz[i] - rcd_z[i] * t_0_yzz_z_zz[i];

        t_0_yzz_zz_yz[i] = t_0_yzz_z_yzz[i] - rcd_z[i] * t_0_yzz_z_yz[i];

        t_0_yzz_zz_yy[i] = t_0_yzz_z_yyz[i] - rcd_z[i] * t_0_yzz_z_yy[i];

        t_0_yzz_zz_xz[i] = t_0_yzz_z_xzz[i] - rcd_z[i] * t_0_yzz_z_xz[i];

        t_0_yzz_zz_xy[i] = t_0_yzz_z_xyz[i] - rcd_z[i] * t_0_yzz_z_xy[i];

        t_0_yzz_zz_xx[i] = t_0_yzz_z_xxz[i] - rcd_z[i] * t_0_yzz_z_xx[i];

        t_0_yzz_yz_zz[i] = t_0_yzz_y_zzz[i] - rcd_z[i] * t_0_yzz_y_zz[i];

        t_0_yzz_yz_yz[i] = t_0_yzz_y_yzz[i] - rcd_z[i] * t_0_yzz_y_yz[i];

        t_0_yzz_yz_yy[i] = t_0_yzz_y_yyz[i] - rcd_z[i] * t_0_yzz_y_yy[i];

        t_0_yzz_yz_xz[i] = t_0_yzz_y_xzz[i] - rcd_z[i] * t_0_yzz_y_xz[i];

        t_0_yzz_yz_xy[i] = t_0_yzz_y_xyz[i] - rcd_z[i] * t_0_yzz_y_xy[i];

        t_0_yzz_yz_xx[i] = t_0_yzz_y_xxz[i] - rcd_z[i] * t_0_yzz_y_xx[i];

        t_0_yzz_yy_zz[i] = t_0_yzz_y_yzz[i] - rcd_y[i] * t_0_yzz_y_zz[i];

        t_0_yzz_yy_yz[i] = t_0_yzz_y_yyz[i] - rcd_y[i] * t_0_yzz_y_yz[i];

        t_0_yzz_yy_yy[i] = t_0_yzz_y_yyy[i] - rcd_y[i] * t_0_yzz_y_yy[i];

        t_0_yzz_yy_xz[i] = t_0_yzz_y_xyz[i] - rcd_y[i] * t_0_yzz_y_xz[i];

        t_0_yzz_yy_xy[i] = t_0_yzz_y_xyy[i] - rcd_y[i] * t_0_yzz_y_xy[i];

        t_0_yzz_yy_xx[i] = t_0_yzz_y_xxy[i] - rcd_y[i] * t_0_yzz_y_xx[i];

        t_0_yzz_xz_zz[i] = t_0_yzz_x_zzz[i] - rcd_z[i] * t_0_yzz_x_zz[i];

        t_0_yzz_xz_yz[i] = t_0_yzz_x_yzz[i] - rcd_z[i] * t_0_yzz_x_yz[i];

        t_0_yzz_xz_yy[i] = t_0_yzz_x_yyz[i] - rcd_z[i] * t_0_yzz_x_yy[i];

        t_0_yzz_xz_xz[i] = t_0_yzz_x_xzz[i] - rcd_z[i] * t_0_yzz_x_xz[i];

        t_0_yzz_xz_xy[i] = t_0_yzz_x_xyz[i] - rcd_z[i] * t_0_yzz_x_xy[i];

        t_0_yzz_xz_xx[i] = t_0_yzz_x_xxz[i] - rcd_z[i] * t_0_yzz_x_xx[i];

        t_0_yzz_xy_zz[i] = t_0_yzz_x_yzz[i] - rcd_y[i] * t_0_yzz_x_zz[i];

        t_0_yzz_xy_yz[i] = t_0_yzz_x_yyz[i] - rcd_y[i] * t_0_yzz_x_yz[i];

        t_0_yzz_xy_yy[i] = t_0_yzz_x_yyy[i] - rcd_y[i] * t_0_yzz_x_yy[i];

        t_0_yzz_xy_xz[i] = t_0_yzz_x_xyz[i] - rcd_y[i] * t_0_yzz_x_xz[i];

        t_0_yzz_xy_xy[i] = t_0_yzz_x_xyy[i] - rcd_y[i] * t_0_yzz_x_xy[i];

        t_0_yzz_xy_xx[i] = t_0_yzz_x_xxy[i] - rcd_y[i] * t_0_yzz_x_xx[i];

        t_0_yzz_xx_zz[i] = t_0_yzz_x_xzz[i] - rcd_x[i] * t_0_yzz_x_zz[i];

        t_0_yzz_xx_yz[i] = t_0_yzz_x_xyz[i] - rcd_x[i] * t_0_yzz_x_yz[i];

        t_0_yzz_xx_yy[i] = t_0_yzz_x_xyy[i] - rcd_x[i] * t_0_yzz_x_yy[i];

        t_0_yzz_xx_xz[i] = t_0_yzz_x_xxz[i] - rcd_x[i] * t_0_yzz_x_xz[i];

        t_0_yzz_xx_xy[i] = t_0_yzz_x_xxy[i] - rcd_x[i] * t_0_yzz_x_xy[i];

        t_0_yzz_xx_xx[i] = t_0_yzz_x_xxx[i] - rcd_x[i] * t_0_yzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyz_x_xx, t_0_yyz_x_xxx, t_0_yyz_x_xxy,\
                           t_0_yyz_x_xxz, t_0_yyz_x_xy, t_0_yyz_x_xyy, t_0_yyz_x_xyz,\
                           t_0_yyz_x_xz, t_0_yyz_x_xzz, t_0_yyz_x_yy, t_0_yyz_x_yyy,\
                           t_0_yyz_x_yyz, t_0_yyz_x_yz, t_0_yyz_x_yzz, t_0_yyz_x_zz,\
                           t_0_yyz_x_zzz, t_0_yyz_xx_xx, t_0_yyz_xx_xy, t_0_yyz_xx_xz,\
                           t_0_yyz_xx_yy, t_0_yyz_xx_yz, t_0_yyz_xx_zz, t_0_yyz_xy_xx,\
                           t_0_yyz_xy_xy, t_0_yyz_xy_xz, t_0_yyz_xy_yy, t_0_yyz_xy_yz,\
                           t_0_yyz_xy_zz, t_0_yyz_xz_xx, t_0_yyz_xz_xy, t_0_yyz_xz_xz,\
                           t_0_yyz_xz_yy, t_0_yyz_xz_yz, t_0_yyz_xz_zz, t_0_yyz_y_xx,\
                           t_0_yyz_y_xxy, t_0_yyz_y_xxz, t_0_yyz_y_xy, t_0_yyz_y_xyy,\
                           t_0_yyz_y_xyz, t_0_yyz_y_xz, t_0_yyz_y_xzz, t_0_yyz_y_yy,\
                           t_0_yyz_y_yyy, t_0_yyz_y_yyz, t_0_yyz_y_yz, t_0_yyz_y_yzz,\
                           t_0_yyz_y_zz, t_0_yyz_y_zzz, t_0_yyz_yy_xx, t_0_yyz_yy_xy,\
                           t_0_yyz_yy_xz, t_0_yyz_yy_yy, t_0_yyz_yy_yz, t_0_yyz_yy_zz,\
                           t_0_yyz_yz_xx, t_0_yyz_yz_xy, t_0_yyz_yz_xz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_yz_zz, t_0_yyz_z_xx, t_0_yyz_z_xxz,\
                           t_0_yyz_z_xy, t_0_yyz_z_xyz, t_0_yyz_z_xz, t_0_yyz_z_xzz,\
                           t_0_yyz_z_yy, t_0_yyz_z_yyz, t_0_yyz_z_yz, t_0_yyz_z_yzz,\
                           t_0_yyz_z_zz, t_0_yyz_z_zzz, t_0_yyz_zz_xx, t_0_yyz_zz_xy,\
                           t_0_yyz_zz_xz, t_0_yyz_zz_yy, t_0_yyz_zz_yz, t_0_yyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyz_zz_zz[i] = t_0_yyz_z_zzz[i] - rcd_z[i] * t_0_yyz_z_zz[i];

        t_0_yyz_zz_yz[i] = t_0_yyz_z_yzz[i] - rcd_z[i] * t_0_yyz_z_yz[i];

        t_0_yyz_zz_yy[i] = t_0_yyz_z_yyz[i] - rcd_z[i] * t_0_yyz_z_yy[i];

        t_0_yyz_zz_xz[i] = t_0_yyz_z_xzz[i] - rcd_z[i] * t_0_yyz_z_xz[i];

        t_0_yyz_zz_xy[i] = t_0_yyz_z_xyz[i] - rcd_z[i] * t_0_yyz_z_xy[i];

        t_0_yyz_zz_xx[i] = t_0_yyz_z_xxz[i] - rcd_z[i] * t_0_yyz_z_xx[i];

        t_0_yyz_yz_zz[i] = t_0_yyz_y_zzz[i] - rcd_z[i] * t_0_yyz_y_zz[i];

        t_0_yyz_yz_yz[i] = t_0_yyz_y_yzz[i] - rcd_z[i] * t_0_yyz_y_yz[i];

        t_0_yyz_yz_yy[i] = t_0_yyz_y_yyz[i] - rcd_z[i] * t_0_yyz_y_yy[i];

        t_0_yyz_yz_xz[i] = t_0_yyz_y_xzz[i] - rcd_z[i] * t_0_yyz_y_xz[i];

        t_0_yyz_yz_xy[i] = t_0_yyz_y_xyz[i] - rcd_z[i] * t_0_yyz_y_xy[i];

        t_0_yyz_yz_xx[i] = t_0_yyz_y_xxz[i] - rcd_z[i] * t_0_yyz_y_xx[i];

        t_0_yyz_yy_zz[i] = t_0_yyz_y_yzz[i] - rcd_y[i] * t_0_yyz_y_zz[i];

        t_0_yyz_yy_yz[i] = t_0_yyz_y_yyz[i] - rcd_y[i] * t_0_yyz_y_yz[i];

        t_0_yyz_yy_yy[i] = t_0_yyz_y_yyy[i] - rcd_y[i] * t_0_yyz_y_yy[i];

        t_0_yyz_yy_xz[i] = t_0_yyz_y_xyz[i] - rcd_y[i] * t_0_yyz_y_xz[i];

        t_0_yyz_yy_xy[i] = t_0_yyz_y_xyy[i] - rcd_y[i] * t_0_yyz_y_xy[i];

        t_0_yyz_yy_xx[i] = t_0_yyz_y_xxy[i] - rcd_y[i] * t_0_yyz_y_xx[i];

        t_0_yyz_xz_zz[i] = t_0_yyz_x_zzz[i] - rcd_z[i] * t_0_yyz_x_zz[i];

        t_0_yyz_xz_yz[i] = t_0_yyz_x_yzz[i] - rcd_z[i] * t_0_yyz_x_yz[i];

        t_0_yyz_xz_yy[i] = t_0_yyz_x_yyz[i] - rcd_z[i] * t_0_yyz_x_yy[i];

        t_0_yyz_xz_xz[i] = t_0_yyz_x_xzz[i] - rcd_z[i] * t_0_yyz_x_xz[i];

        t_0_yyz_xz_xy[i] = t_0_yyz_x_xyz[i] - rcd_z[i] * t_0_yyz_x_xy[i];

        t_0_yyz_xz_xx[i] = t_0_yyz_x_xxz[i] - rcd_z[i] * t_0_yyz_x_xx[i];

        t_0_yyz_xy_zz[i] = t_0_yyz_x_yzz[i] - rcd_y[i] * t_0_yyz_x_zz[i];

        t_0_yyz_xy_yz[i] = t_0_yyz_x_yyz[i] - rcd_y[i] * t_0_yyz_x_yz[i];

        t_0_yyz_xy_yy[i] = t_0_yyz_x_yyy[i] - rcd_y[i] * t_0_yyz_x_yy[i];

        t_0_yyz_xy_xz[i] = t_0_yyz_x_xyz[i] - rcd_y[i] * t_0_yyz_x_xz[i];

        t_0_yyz_xy_xy[i] = t_0_yyz_x_xyy[i] - rcd_y[i] * t_0_yyz_x_xy[i];

        t_0_yyz_xy_xx[i] = t_0_yyz_x_xxy[i] - rcd_y[i] * t_0_yyz_x_xx[i];

        t_0_yyz_xx_zz[i] = t_0_yyz_x_xzz[i] - rcd_x[i] * t_0_yyz_x_zz[i];

        t_0_yyz_xx_yz[i] = t_0_yyz_x_xyz[i] - rcd_x[i] * t_0_yyz_x_yz[i];

        t_0_yyz_xx_yy[i] = t_0_yyz_x_xyy[i] - rcd_x[i] * t_0_yyz_x_yy[i];

        t_0_yyz_xx_xz[i] = t_0_yyz_x_xxz[i] - rcd_x[i] * t_0_yyz_x_xz[i];

        t_0_yyz_xx_xy[i] = t_0_yyz_x_xxy[i] - rcd_x[i] * t_0_yyz_x_xy[i];

        t_0_yyz_xx_xx[i] = t_0_yyz_x_xxx[i] - rcd_x[i] * t_0_yyz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyy_x_xx, t_0_yyy_x_xxx, t_0_yyy_x_xxy,\
                           t_0_yyy_x_xxz, t_0_yyy_x_xy, t_0_yyy_x_xyy, t_0_yyy_x_xyz,\
                           t_0_yyy_x_xz, t_0_yyy_x_xzz, t_0_yyy_x_yy, t_0_yyy_x_yyy,\
                           t_0_yyy_x_yyz, t_0_yyy_x_yz, t_0_yyy_x_yzz, t_0_yyy_x_zz,\
                           t_0_yyy_x_zzz, t_0_yyy_xx_xx, t_0_yyy_xx_xy, t_0_yyy_xx_xz,\
                           t_0_yyy_xx_yy, t_0_yyy_xx_yz, t_0_yyy_xx_zz, t_0_yyy_xy_xx,\
                           t_0_yyy_xy_xy, t_0_yyy_xy_xz, t_0_yyy_xy_yy, t_0_yyy_xy_yz,\
                           t_0_yyy_xy_zz, t_0_yyy_xz_xx, t_0_yyy_xz_xy, t_0_yyy_xz_xz,\
                           t_0_yyy_xz_yy, t_0_yyy_xz_yz, t_0_yyy_xz_zz, t_0_yyy_y_xx,\
                           t_0_yyy_y_xxy, t_0_yyy_y_xxz, t_0_yyy_y_xy, t_0_yyy_y_xyy,\
                           t_0_yyy_y_xyz, t_0_yyy_y_xz, t_0_yyy_y_xzz, t_0_yyy_y_yy,\
                           t_0_yyy_y_yyy, t_0_yyy_y_yyz, t_0_yyy_y_yz, t_0_yyy_y_yzz,\
                           t_0_yyy_y_zz, t_0_yyy_y_zzz, t_0_yyy_yy_xx, t_0_yyy_yy_xy,\
                           t_0_yyy_yy_xz, t_0_yyy_yy_yy, t_0_yyy_yy_yz, t_0_yyy_yy_zz,\
                           t_0_yyy_yz_xx, t_0_yyy_yz_xy, t_0_yyy_yz_xz, t_0_yyy_yz_yy,\
                           t_0_yyy_yz_yz, t_0_yyy_yz_zz, t_0_yyy_z_xx, t_0_yyy_z_xxz,\
                           t_0_yyy_z_xy, t_0_yyy_z_xyz, t_0_yyy_z_xz, t_0_yyy_z_xzz,\
                           t_0_yyy_z_yy, t_0_yyy_z_yyz, t_0_yyy_z_yz, t_0_yyy_z_yzz,\
                           t_0_yyy_z_zz, t_0_yyy_z_zzz, t_0_yyy_zz_xx, t_0_yyy_zz_xy,\
                           t_0_yyy_zz_xz, t_0_yyy_zz_yy, t_0_yyy_zz_yz, t_0_yyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyy_zz_zz[i] = t_0_yyy_z_zzz[i] - rcd_z[i] * t_0_yyy_z_zz[i];

        t_0_yyy_zz_yz[i] = t_0_yyy_z_yzz[i] - rcd_z[i] * t_0_yyy_z_yz[i];

        t_0_yyy_zz_yy[i] = t_0_yyy_z_yyz[i] - rcd_z[i] * t_0_yyy_z_yy[i];

        t_0_yyy_zz_xz[i] = t_0_yyy_z_xzz[i] - rcd_z[i] * t_0_yyy_z_xz[i];

        t_0_yyy_zz_xy[i] = t_0_yyy_z_xyz[i] - rcd_z[i] * t_0_yyy_z_xy[i];

        t_0_yyy_zz_xx[i] = t_0_yyy_z_xxz[i] - rcd_z[i] * t_0_yyy_z_xx[i];

        t_0_yyy_yz_zz[i] = t_0_yyy_y_zzz[i] - rcd_z[i] * t_0_yyy_y_zz[i];

        t_0_yyy_yz_yz[i] = t_0_yyy_y_yzz[i] - rcd_z[i] * t_0_yyy_y_yz[i];

        t_0_yyy_yz_yy[i] = t_0_yyy_y_yyz[i] - rcd_z[i] * t_0_yyy_y_yy[i];

        t_0_yyy_yz_xz[i] = t_0_yyy_y_xzz[i] - rcd_z[i] * t_0_yyy_y_xz[i];

        t_0_yyy_yz_xy[i] = t_0_yyy_y_xyz[i] - rcd_z[i] * t_0_yyy_y_xy[i];

        t_0_yyy_yz_xx[i] = t_0_yyy_y_xxz[i] - rcd_z[i] * t_0_yyy_y_xx[i];

        t_0_yyy_yy_zz[i] = t_0_yyy_y_yzz[i] - rcd_y[i] * t_0_yyy_y_zz[i];

        t_0_yyy_yy_yz[i] = t_0_yyy_y_yyz[i] - rcd_y[i] * t_0_yyy_y_yz[i];

        t_0_yyy_yy_yy[i] = t_0_yyy_y_yyy[i] - rcd_y[i] * t_0_yyy_y_yy[i];

        t_0_yyy_yy_xz[i] = t_0_yyy_y_xyz[i] - rcd_y[i] * t_0_yyy_y_xz[i];

        t_0_yyy_yy_xy[i] = t_0_yyy_y_xyy[i] - rcd_y[i] * t_0_yyy_y_xy[i];

        t_0_yyy_yy_xx[i] = t_0_yyy_y_xxy[i] - rcd_y[i] * t_0_yyy_y_xx[i];

        t_0_yyy_xz_zz[i] = t_0_yyy_x_zzz[i] - rcd_z[i] * t_0_yyy_x_zz[i];

        t_0_yyy_xz_yz[i] = t_0_yyy_x_yzz[i] - rcd_z[i] * t_0_yyy_x_yz[i];

        t_0_yyy_xz_yy[i] = t_0_yyy_x_yyz[i] - rcd_z[i] * t_0_yyy_x_yy[i];

        t_0_yyy_xz_xz[i] = t_0_yyy_x_xzz[i] - rcd_z[i] * t_0_yyy_x_xz[i];

        t_0_yyy_xz_xy[i] = t_0_yyy_x_xyz[i] - rcd_z[i] * t_0_yyy_x_xy[i];

        t_0_yyy_xz_xx[i] = t_0_yyy_x_xxz[i] - rcd_z[i] * t_0_yyy_x_xx[i];

        t_0_yyy_xy_zz[i] = t_0_yyy_x_yzz[i] - rcd_y[i] * t_0_yyy_x_zz[i];

        t_0_yyy_xy_yz[i] = t_0_yyy_x_yyz[i] - rcd_y[i] * t_0_yyy_x_yz[i];

        t_0_yyy_xy_yy[i] = t_0_yyy_x_yyy[i] - rcd_y[i] * t_0_yyy_x_yy[i];

        t_0_yyy_xy_xz[i] = t_0_yyy_x_xyz[i] - rcd_y[i] * t_0_yyy_x_xz[i];

        t_0_yyy_xy_xy[i] = t_0_yyy_x_xyy[i] - rcd_y[i] * t_0_yyy_x_xy[i];

        t_0_yyy_xy_xx[i] = t_0_yyy_x_xxy[i] - rcd_y[i] * t_0_yyy_x_xx[i];

        t_0_yyy_xx_zz[i] = t_0_yyy_x_xzz[i] - rcd_x[i] * t_0_yyy_x_zz[i];

        t_0_yyy_xx_yz[i] = t_0_yyy_x_xyz[i] - rcd_x[i] * t_0_yyy_x_yz[i];

        t_0_yyy_xx_yy[i] = t_0_yyy_x_xyy[i] - rcd_x[i] * t_0_yyy_x_yy[i];

        t_0_yyy_xx_xz[i] = t_0_yyy_x_xxz[i] - rcd_x[i] * t_0_yyy_x_xz[i];

        t_0_yyy_xx_xy[i] = t_0_yyy_x_xxy[i] - rcd_x[i] * t_0_yyy_x_xy[i];

        t_0_yyy_xx_xx[i] = t_0_yyy_x_xxx[i] - rcd_x[i] * t_0_yyy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xzz_x_xx, t_0_xzz_x_xxx, t_0_xzz_x_xxy,\
                           t_0_xzz_x_xxz, t_0_xzz_x_xy, t_0_xzz_x_xyy, t_0_xzz_x_xyz,\
                           t_0_xzz_x_xz, t_0_xzz_x_xzz, t_0_xzz_x_yy, t_0_xzz_x_yyy,\
                           t_0_xzz_x_yyz, t_0_xzz_x_yz, t_0_xzz_x_yzz, t_0_xzz_x_zz,\
                           t_0_xzz_x_zzz, t_0_xzz_xx_xx, t_0_xzz_xx_xy, t_0_xzz_xx_xz,\
                           t_0_xzz_xx_yy, t_0_xzz_xx_yz, t_0_xzz_xx_zz, t_0_xzz_xy_xx,\
                           t_0_xzz_xy_xy, t_0_xzz_xy_xz, t_0_xzz_xy_yy, t_0_xzz_xy_yz,\
                           t_0_xzz_xy_zz, t_0_xzz_xz_xx, t_0_xzz_xz_xy, t_0_xzz_xz_xz,\
                           t_0_xzz_xz_yy, t_0_xzz_xz_yz, t_0_xzz_xz_zz, t_0_xzz_y_xx,\
                           t_0_xzz_y_xxy, t_0_xzz_y_xxz, t_0_xzz_y_xy, t_0_xzz_y_xyy,\
                           t_0_xzz_y_xyz, t_0_xzz_y_xz, t_0_xzz_y_xzz, t_0_xzz_y_yy,\
                           t_0_xzz_y_yyy, t_0_xzz_y_yyz, t_0_xzz_y_yz, t_0_xzz_y_yzz,\
                           t_0_xzz_y_zz, t_0_xzz_y_zzz, t_0_xzz_yy_xx, t_0_xzz_yy_xy,\
                           t_0_xzz_yy_xz, t_0_xzz_yy_yy, t_0_xzz_yy_yz, t_0_xzz_yy_zz,\
                           t_0_xzz_yz_xx, t_0_xzz_yz_xy, t_0_xzz_yz_xz, t_0_xzz_yz_yy,\
                           t_0_xzz_yz_yz, t_0_xzz_yz_zz, t_0_xzz_z_xx, t_0_xzz_z_xxz,\
                           t_0_xzz_z_xy, t_0_xzz_z_xyz, t_0_xzz_z_xz, t_0_xzz_z_xzz,\
                           t_0_xzz_z_yy, t_0_xzz_z_yyz, t_0_xzz_z_yz, t_0_xzz_z_yzz,\
                           t_0_xzz_z_zz, t_0_xzz_z_zzz, t_0_xzz_zz_xx, t_0_xzz_zz_xy,\
                           t_0_xzz_zz_xz, t_0_xzz_zz_yy, t_0_xzz_zz_yz, t_0_xzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzz_zz_zz[i] = t_0_xzz_z_zzz[i] - rcd_z[i] * t_0_xzz_z_zz[i];

        t_0_xzz_zz_yz[i] = t_0_xzz_z_yzz[i] - rcd_z[i] * t_0_xzz_z_yz[i];

        t_0_xzz_zz_yy[i] = t_0_xzz_z_yyz[i] - rcd_z[i] * t_0_xzz_z_yy[i];

        t_0_xzz_zz_xz[i] = t_0_xzz_z_xzz[i] - rcd_z[i] * t_0_xzz_z_xz[i];

        t_0_xzz_zz_xy[i] = t_0_xzz_z_xyz[i] - rcd_z[i] * t_0_xzz_z_xy[i];

        t_0_xzz_zz_xx[i] = t_0_xzz_z_xxz[i] - rcd_z[i] * t_0_xzz_z_xx[i];

        t_0_xzz_yz_zz[i] = t_0_xzz_y_zzz[i] - rcd_z[i] * t_0_xzz_y_zz[i];

        t_0_xzz_yz_yz[i] = t_0_xzz_y_yzz[i] - rcd_z[i] * t_0_xzz_y_yz[i];

        t_0_xzz_yz_yy[i] = t_0_xzz_y_yyz[i] - rcd_z[i] * t_0_xzz_y_yy[i];

        t_0_xzz_yz_xz[i] = t_0_xzz_y_xzz[i] - rcd_z[i] * t_0_xzz_y_xz[i];

        t_0_xzz_yz_xy[i] = t_0_xzz_y_xyz[i] - rcd_z[i] * t_0_xzz_y_xy[i];

        t_0_xzz_yz_xx[i] = t_0_xzz_y_xxz[i] - rcd_z[i] * t_0_xzz_y_xx[i];

        t_0_xzz_yy_zz[i] = t_0_xzz_y_yzz[i] - rcd_y[i] * t_0_xzz_y_zz[i];

        t_0_xzz_yy_yz[i] = t_0_xzz_y_yyz[i] - rcd_y[i] * t_0_xzz_y_yz[i];

        t_0_xzz_yy_yy[i] = t_0_xzz_y_yyy[i] - rcd_y[i] * t_0_xzz_y_yy[i];

        t_0_xzz_yy_xz[i] = t_0_xzz_y_xyz[i] - rcd_y[i] * t_0_xzz_y_xz[i];

        t_0_xzz_yy_xy[i] = t_0_xzz_y_xyy[i] - rcd_y[i] * t_0_xzz_y_xy[i];

        t_0_xzz_yy_xx[i] = t_0_xzz_y_xxy[i] - rcd_y[i] * t_0_xzz_y_xx[i];

        t_0_xzz_xz_zz[i] = t_0_xzz_x_zzz[i] - rcd_z[i] * t_0_xzz_x_zz[i];

        t_0_xzz_xz_yz[i] = t_0_xzz_x_yzz[i] - rcd_z[i] * t_0_xzz_x_yz[i];

        t_0_xzz_xz_yy[i] = t_0_xzz_x_yyz[i] - rcd_z[i] * t_0_xzz_x_yy[i];

        t_0_xzz_xz_xz[i] = t_0_xzz_x_xzz[i] - rcd_z[i] * t_0_xzz_x_xz[i];

        t_0_xzz_xz_xy[i] = t_0_xzz_x_xyz[i] - rcd_z[i] * t_0_xzz_x_xy[i];

        t_0_xzz_xz_xx[i] = t_0_xzz_x_xxz[i] - rcd_z[i] * t_0_xzz_x_xx[i];

        t_0_xzz_xy_zz[i] = t_0_xzz_x_yzz[i] - rcd_y[i] * t_0_xzz_x_zz[i];

        t_0_xzz_xy_yz[i] = t_0_xzz_x_yyz[i] - rcd_y[i] * t_0_xzz_x_yz[i];

        t_0_xzz_xy_yy[i] = t_0_xzz_x_yyy[i] - rcd_y[i] * t_0_xzz_x_yy[i];

        t_0_xzz_xy_xz[i] = t_0_xzz_x_xyz[i] - rcd_y[i] * t_0_xzz_x_xz[i];

        t_0_xzz_xy_xy[i] = t_0_xzz_x_xyy[i] - rcd_y[i] * t_0_xzz_x_xy[i];

        t_0_xzz_xy_xx[i] = t_0_xzz_x_xxy[i] - rcd_y[i] * t_0_xzz_x_xx[i];

        t_0_xzz_xx_zz[i] = t_0_xzz_x_xzz[i] - rcd_x[i] * t_0_xzz_x_zz[i];

        t_0_xzz_xx_yz[i] = t_0_xzz_x_xyz[i] - rcd_x[i] * t_0_xzz_x_yz[i];

        t_0_xzz_xx_yy[i] = t_0_xzz_x_xyy[i] - rcd_x[i] * t_0_xzz_x_yy[i];

        t_0_xzz_xx_xz[i] = t_0_xzz_x_xxz[i] - rcd_x[i] * t_0_xzz_x_xz[i];

        t_0_xzz_xx_xy[i] = t_0_xzz_x_xxy[i] - rcd_x[i] * t_0_xzz_x_xy[i];

        t_0_xzz_xx_xx[i] = t_0_xzz_x_xxx[i] - rcd_x[i] * t_0_xzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyz_x_xx, t_0_xyz_x_xxx, t_0_xyz_x_xxy,\
                           t_0_xyz_x_xxz, t_0_xyz_x_xy, t_0_xyz_x_xyy, t_0_xyz_x_xyz,\
                           t_0_xyz_x_xz, t_0_xyz_x_xzz, t_0_xyz_x_yy, t_0_xyz_x_yyy,\
                           t_0_xyz_x_yyz, t_0_xyz_x_yz, t_0_xyz_x_yzz, t_0_xyz_x_zz,\
                           t_0_xyz_x_zzz, t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz,\
                           t_0_xyz_xx_yy, t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx,\
                           t_0_xyz_xy_xy, t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz,\
                           t_0_xyz_xy_zz, t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz,\
                           t_0_xyz_xz_yy, t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_y_xx,\
                           t_0_xyz_y_xxy, t_0_xyz_y_xxz, t_0_xyz_y_xy, t_0_xyz_y_xyy,\
                           t_0_xyz_y_xyz, t_0_xyz_y_xz, t_0_xyz_y_xzz, t_0_xyz_y_yy,\
                           t_0_xyz_y_yyy, t_0_xyz_y_yyz, t_0_xyz_y_yz, t_0_xyz_y_yzz,\
                           t_0_xyz_y_zz, t_0_xyz_y_zzz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_z_xx, t_0_xyz_z_xxz,\
                           t_0_xyz_z_xy, t_0_xyz_z_xyz, t_0_xyz_z_xz, t_0_xyz_z_xzz,\
                           t_0_xyz_z_yy, t_0_xyz_z_yyz, t_0_xyz_z_yz, t_0_xyz_z_yzz,\
                           t_0_xyz_z_zz, t_0_xyz_z_zzz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyz_zz_zz[i] = t_0_xyz_z_zzz[i] - rcd_z[i] * t_0_xyz_z_zz[i];

        t_0_xyz_zz_yz[i] = t_0_xyz_z_yzz[i] - rcd_z[i] * t_0_xyz_z_yz[i];

        t_0_xyz_zz_yy[i] = t_0_xyz_z_yyz[i] - rcd_z[i] * t_0_xyz_z_yy[i];

        t_0_xyz_zz_xz[i] = t_0_xyz_z_xzz[i] - rcd_z[i] * t_0_xyz_z_xz[i];

        t_0_xyz_zz_xy[i] = t_0_xyz_z_xyz[i] - rcd_z[i] * t_0_xyz_z_xy[i];

        t_0_xyz_zz_xx[i] = t_0_xyz_z_xxz[i] - rcd_z[i] * t_0_xyz_z_xx[i];

        t_0_xyz_yz_zz[i] = t_0_xyz_y_zzz[i] - rcd_z[i] * t_0_xyz_y_zz[i];

        t_0_xyz_yz_yz[i] = t_0_xyz_y_yzz[i] - rcd_z[i] * t_0_xyz_y_yz[i];

        t_0_xyz_yz_yy[i] = t_0_xyz_y_yyz[i] - rcd_z[i] * t_0_xyz_y_yy[i];

        t_0_xyz_yz_xz[i] = t_0_xyz_y_xzz[i] - rcd_z[i] * t_0_xyz_y_xz[i];

        t_0_xyz_yz_xy[i] = t_0_xyz_y_xyz[i] - rcd_z[i] * t_0_xyz_y_xy[i];

        t_0_xyz_yz_xx[i] = t_0_xyz_y_xxz[i] - rcd_z[i] * t_0_xyz_y_xx[i];

        t_0_xyz_yy_zz[i] = t_0_xyz_y_yzz[i] - rcd_y[i] * t_0_xyz_y_zz[i];

        t_0_xyz_yy_yz[i] = t_0_xyz_y_yyz[i] - rcd_y[i] * t_0_xyz_y_yz[i];

        t_0_xyz_yy_yy[i] = t_0_xyz_y_yyy[i] - rcd_y[i] * t_0_xyz_y_yy[i];

        t_0_xyz_yy_xz[i] = t_0_xyz_y_xyz[i] - rcd_y[i] * t_0_xyz_y_xz[i];

        t_0_xyz_yy_xy[i] = t_0_xyz_y_xyy[i] - rcd_y[i] * t_0_xyz_y_xy[i];

        t_0_xyz_yy_xx[i] = t_0_xyz_y_xxy[i] - rcd_y[i] * t_0_xyz_y_xx[i];

        t_0_xyz_xz_zz[i] = t_0_xyz_x_zzz[i] - rcd_z[i] * t_0_xyz_x_zz[i];

        t_0_xyz_xz_yz[i] = t_0_xyz_x_yzz[i] - rcd_z[i] * t_0_xyz_x_yz[i];

        t_0_xyz_xz_yy[i] = t_0_xyz_x_yyz[i] - rcd_z[i] * t_0_xyz_x_yy[i];

        t_0_xyz_xz_xz[i] = t_0_xyz_x_xzz[i] - rcd_z[i] * t_0_xyz_x_xz[i];

        t_0_xyz_xz_xy[i] = t_0_xyz_x_xyz[i] - rcd_z[i] * t_0_xyz_x_xy[i];

        t_0_xyz_xz_xx[i] = t_0_xyz_x_xxz[i] - rcd_z[i] * t_0_xyz_x_xx[i];

        t_0_xyz_xy_zz[i] = t_0_xyz_x_yzz[i] - rcd_y[i] * t_0_xyz_x_zz[i];

        t_0_xyz_xy_yz[i] = t_0_xyz_x_yyz[i] - rcd_y[i] * t_0_xyz_x_yz[i];

        t_0_xyz_xy_yy[i] = t_0_xyz_x_yyy[i] - rcd_y[i] * t_0_xyz_x_yy[i];

        t_0_xyz_xy_xz[i] = t_0_xyz_x_xyz[i] - rcd_y[i] * t_0_xyz_x_xz[i];

        t_0_xyz_xy_xy[i] = t_0_xyz_x_xyy[i] - rcd_y[i] * t_0_xyz_x_xy[i];

        t_0_xyz_xy_xx[i] = t_0_xyz_x_xxy[i] - rcd_y[i] * t_0_xyz_x_xx[i];

        t_0_xyz_xx_zz[i] = t_0_xyz_x_xzz[i] - rcd_x[i] * t_0_xyz_x_zz[i];

        t_0_xyz_xx_yz[i] = t_0_xyz_x_xyz[i] - rcd_x[i] * t_0_xyz_x_yz[i];

        t_0_xyz_xx_yy[i] = t_0_xyz_x_xyy[i] - rcd_x[i] * t_0_xyz_x_yy[i];

        t_0_xyz_xx_xz[i] = t_0_xyz_x_xxz[i] - rcd_x[i] * t_0_xyz_x_xz[i];

        t_0_xyz_xx_xy[i] = t_0_xyz_x_xxy[i] - rcd_x[i] * t_0_xyz_x_xy[i];

        t_0_xyz_xx_xx[i] = t_0_xyz_x_xxx[i] - rcd_x[i] * t_0_xyz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyy_x_xx, t_0_xyy_x_xxx, t_0_xyy_x_xxy,\
                           t_0_xyy_x_xxz, t_0_xyy_x_xy, t_0_xyy_x_xyy, t_0_xyy_x_xyz,\
                           t_0_xyy_x_xz, t_0_xyy_x_xzz, t_0_xyy_x_yy, t_0_xyy_x_yyy,\
                           t_0_xyy_x_yyz, t_0_xyy_x_yz, t_0_xyy_x_yzz, t_0_xyy_x_zz,\
                           t_0_xyy_x_zzz, t_0_xyy_xx_xx, t_0_xyy_xx_xy, t_0_xyy_xx_xz,\
                           t_0_xyy_xx_yy, t_0_xyy_xx_yz, t_0_xyy_xx_zz, t_0_xyy_xy_xx,\
                           t_0_xyy_xy_xy, t_0_xyy_xy_xz, t_0_xyy_xy_yy, t_0_xyy_xy_yz,\
                           t_0_xyy_xy_zz, t_0_xyy_xz_xx, t_0_xyy_xz_xy, t_0_xyy_xz_xz,\
                           t_0_xyy_xz_yy, t_0_xyy_xz_yz, t_0_xyy_xz_zz, t_0_xyy_y_xx,\
                           t_0_xyy_y_xxy, t_0_xyy_y_xxz, t_0_xyy_y_xy, t_0_xyy_y_xyy,\
                           t_0_xyy_y_xyz, t_0_xyy_y_xz, t_0_xyy_y_xzz, t_0_xyy_y_yy,\
                           t_0_xyy_y_yyy, t_0_xyy_y_yyz, t_0_xyy_y_yz, t_0_xyy_y_yzz,\
                           t_0_xyy_y_zz, t_0_xyy_y_zzz, t_0_xyy_yy_xx, t_0_xyy_yy_xy,\
                           t_0_xyy_yy_xz, t_0_xyy_yy_yy, t_0_xyy_yy_yz, t_0_xyy_yy_zz,\
                           t_0_xyy_yz_xx, t_0_xyy_yz_xy, t_0_xyy_yz_xz, t_0_xyy_yz_yy,\
                           t_0_xyy_yz_yz, t_0_xyy_yz_zz, t_0_xyy_z_xx, t_0_xyy_z_xxz,\
                           t_0_xyy_z_xy, t_0_xyy_z_xyz, t_0_xyy_z_xz, t_0_xyy_z_xzz,\
                           t_0_xyy_z_yy, t_0_xyy_z_yyz, t_0_xyy_z_yz, t_0_xyy_z_yzz,\
                           t_0_xyy_z_zz, t_0_xyy_z_zzz, t_0_xyy_zz_xx, t_0_xyy_zz_xy,\
                           t_0_xyy_zz_xz, t_0_xyy_zz_yy, t_0_xyy_zz_yz, t_0_xyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyy_zz_zz[i] = t_0_xyy_z_zzz[i] - rcd_z[i] * t_0_xyy_z_zz[i];

        t_0_xyy_zz_yz[i] = t_0_xyy_z_yzz[i] - rcd_z[i] * t_0_xyy_z_yz[i];

        t_0_xyy_zz_yy[i] = t_0_xyy_z_yyz[i] - rcd_z[i] * t_0_xyy_z_yy[i];

        t_0_xyy_zz_xz[i] = t_0_xyy_z_xzz[i] - rcd_z[i] * t_0_xyy_z_xz[i];

        t_0_xyy_zz_xy[i] = t_0_xyy_z_xyz[i] - rcd_z[i] * t_0_xyy_z_xy[i];

        t_0_xyy_zz_xx[i] = t_0_xyy_z_xxz[i] - rcd_z[i] * t_0_xyy_z_xx[i];

        t_0_xyy_yz_zz[i] = t_0_xyy_y_zzz[i] - rcd_z[i] * t_0_xyy_y_zz[i];

        t_0_xyy_yz_yz[i] = t_0_xyy_y_yzz[i] - rcd_z[i] * t_0_xyy_y_yz[i];

        t_0_xyy_yz_yy[i] = t_0_xyy_y_yyz[i] - rcd_z[i] * t_0_xyy_y_yy[i];

        t_0_xyy_yz_xz[i] = t_0_xyy_y_xzz[i] - rcd_z[i] * t_0_xyy_y_xz[i];

        t_0_xyy_yz_xy[i] = t_0_xyy_y_xyz[i] - rcd_z[i] * t_0_xyy_y_xy[i];

        t_0_xyy_yz_xx[i] = t_0_xyy_y_xxz[i] - rcd_z[i] * t_0_xyy_y_xx[i];

        t_0_xyy_yy_zz[i] = t_0_xyy_y_yzz[i] - rcd_y[i] * t_0_xyy_y_zz[i];

        t_0_xyy_yy_yz[i] = t_0_xyy_y_yyz[i] - rcd_y[i] * t_0_xyy_y_yz[i];

        t_0_xyy_yy_yy[i] = t_0_xyy_y_yyy[i] - rcd_y[i] * t_0_xyy_y_yy[i];

        t_0_xyy_yy_xz[i] = t_0_xyy_y_xyz[i] - rcd_y[i] * t_0_xyy_y_xz[i];

        t_0_xyy_yy_xy[i] = t_0_xyy_y_xyy[i] - rcd_y[i] * t_0_xyy_y_xy[i];

        t_0_xyy_yy_xx[i] = t_0_xyy_y_xxy[i] - rcd_y[i] * t_0_xyy_y_xx[i];

        t_0_xyy_xz_zz[i] = t_0_xyy_x_zzz[i] - rcd_z[i] * t_0_xyy_x_zz[i];

        t_0_xyy_xz_yz[i] = t_0_xyy_x_yzz[i] - rcd_z[i] * t_0_xyy_x_yz[i];

        t_0_xyy_xz_yy[i] = t_0_xyy_x_yyz[i] - rcd_z[i] * t_0_xyy_x_yy[i];

        t_0_xyy_xz_xz[i] = t_0_xyy_x_xzz[i] - rcd_z[i] * t_0_xyy_x_xz[i];

        t_0_xyy_xz_xy[i] = t_0_xyy_x_xyz[i] - rcd_z[i] * t_0_xyy_x_xy[i];

        t_0_xyy_xz_xx[i] = t_0_xyy_x_xxz[i] - rcd_z[i] * t_0_xyy_x_xx[i];

        t_0_xyy_xy_zz[i] = t_0_xyy_x_yzz[i] - rcd_y[i] * t_0_xyy_x_zz[i];

        t_0_xyy_xy_yz[i] = t_0_xyy_x_yyz[i] - rcd_y[i] * t_0_xyy_x_yz[i];

        t_0_xyy_xy_yy[i] = t_0_xyy_x_yyy[i] - rcd_y[i] * t_0_xyy_x_yy[i];

        t_0_xyy_xy_xz[i] = t_0_xyy_x_xyz[i] - rcd_y[i] * t_0_xyy_x_xz[i];

        t_0_xyy_xy_xy[i] = t_0_xyy_x_xyy[i] - rcd_y[i] * t_0_xyy_x_xy[i];

        t_0_xyy_xy_xx[i] = t_0_xyy_x_xxy[i] - rcd_y[i] * t_0_xyy_x_xx[i];

        t_0_xyy_xx_zz[i] = t_0_xyy_x_xzz[i] - rcd_x[i] * t_0_xyy_x_zz[i];

        t_0_xyy_xx_yz[i] = t_0_xyy_x_xyz[i] - rcd_x[i] * t_0_xyy_x_yz[i];

        t_0_xyy_xx_yy[i] = t_0_xyy_x_xyy[i] - rcd_x[i] * t_0_xyy_x_yy[i];

        t_0_xyy_xx_xz[i] = t_0_xyy_x_xxz[i] - rcd_x[i] * t_0_xyy_x_xz[i];

        t_0_xyy_xx_xy[i] = t_0_xyy_x_xxy[i] - rcd_x[i] * t_0_xyy_x_xy[i];

        t_0_xyy_xx_xx[i] = t_0_xyy_x_xxx[i] - rcd_x[i] * t_0_xyy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxz_x_xx, t_0_xxz_x_xxx, t_0_xxz_x_xxy,\
                           t_0_xxz_x_xxz, t_0_xxz_x_xy, t_0_xxz_x_xyy, t_0_xxz_x_xyz,\
                           t_0_xxz_x_xz, t_0_xxz_x_xzz, t_0_xxz_x_yy, t_0_xxz_x_yyy,\
                           t_0_xxz_x_yyz, t_0_xxz_x_yz, t_0_xxz_x_yzz, t_0_xxz_x_zz,\
                           t_0_xxz_x_zzz, t_0_xxz_xx_xx, t_0_xxz_xx_xy, t_0_xxz_xx_xz,\
                           t_0_xxz_xx_yy, t_0_xxz_xx_yz, t_0_xxz_xx_zz, t_0_xxz_xy_xx,\
                           t_0_xxz_xy_xy, t_0_xxz_xy_xz, t_0_xxz_xy_yy, t_0_xxz_xy_yz,\
                           t_0_xxz_xy_zz, t_0_xxz_xz_xx, t_0_xxz_xz_xy, t_0_xxz_xz_xz,\
                           t_0_xxz_xz_yy, t_0_xxz_xz_yz, t_0_xxz_xz_zz, t_0_xxz_y_xx,\
                           t_0_xxz_y_xxy, t_0_xxz_y_xxz, t_0_xxz_y_xy, t_0_xxz_y_xyy,\
                           t_0_xxz_y_xyz, t_0_xxz_y_xz, t_0_xxz_y_xzz, t_0_xxz_y_yy,\
                           t_0_xxz_y_yyy, t_0_xxz_y_yyz, t_0_xxz_y_yz, t_0_xxz_y_yzz,\
                           t_0_xxz_y_zz, t_0_xxz_y_zzz, t_0_xxz_yy_xx, t_0_xxz_yy_xy,\
                           t_0_xxz_yy_xz, t_0_xxz_yy_yy, t_0_xxz_yy_yz, t_0_xxz_yy_zz,\
                           t_0_xxz_yz_xx, t_0_xxz_yz_xy, t_0_xxz_yz_xz, t_0_xxz_yz_yy,\
                           t_0_xxz_yz_yz, t_0_xxz_yz_zz, t_0_xxz_z_xx, t_0_xxz_z_xxz,\
                           t_0_xxz_z_xy, t_0_xxz_z_xyz, t_0_xxz_z_xz, t_0_xxz_z_xzz,\
                           t_0_xxz_z_yy, t_0_xxz_z_yyz, t_0_xxz_z_yz, t_0_xxz_z_yzz,\
                           t_0_xxz_z_zz, t_0_xxz_z_zzz, t_0_xxz_zz_xx, t_0_xxz_zz_xy,\
                           t_0_xxz_zz_xz, t_0_xxz_zz_yy, t_0_xxz_zz_yz, t_0_xxz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxz_zz_zz[i] = t_0_xxz_z_zzz[i] - rcd_z[i] * t_0_xxz_z_zz[i];

        t_0_xxz_zz_yz[i] = t_0_xxz_z_yzz[i] - rcd_z[i] * t_0_xxz_z_yz[i];

        t_0_xxz_zz_yy[i] = t_0_xxz_z_yyz[i] - rcd_z[i] * t_0_xxz_z_yy[i];

        t_0_xxz_zz_xz[i] = t_0_xxz_z_xzz[i] - rcd_z[i] * t_0_xxz_z_xz[i];

        t_0_xxz_zz_xy[i] = t_0_xxz_z_xyz[i] - rcd_z[i] * t_0_xxz_z_xy[i];

        t_0_xxz_zz_xx[i] = t_0_xxz_z_xxz[i] - rcd_z[i] * t_0_xxz_z_xx[i];

        t_0_xxz_yz_zz[i] = t_0_xxz_y_zzz[i] - rcd_z[i] * t_0_xxz_y_zz[i];

        t_0_xxz_yz_yz[i] = t_0_xxz_y_yzz[i] - rcd_z[i] * t_0_xxz_y_yz[i];

        t_0_xxz_yz_yy[i] = t_0_xxz_y_yyz[i] - rcd_z[i] * t_0_xxz_y_yy[i];

        t_0_xxz_yz_xz[i] = t_0_xxz_y_xzz[i] - rcd_z[i] * t_0_xxz_y_xz[i];

        t_0_xxz_yz_xy[i] = t_0_xxz_y_xyz[i] - rcd_z[i] * t_0_xxz_y_xy[i];

        t_0_xxz_yz_xx[i] = t_0_xxz_y_xxz[i] - rcd_z[i] * t_0_xxz_y_xx[i];

        t_0_xxz_yy_zz[i] = t_0_xxz_y_yzz[i] - rcd_y[i] * t_0_xxz_y_zz[i];

        t_0_xxz_yy_yz[i] = t_0_xxz_y_yyz[i] - rcd_y[i] * t_0_xxz_y_yz[i];

        t_0_xxz_yy_yy[i] = t_0_xxz_y_yyy[i] - rcd_y[i] * t_0_xxz_y_yy[i];

        t_0_xxz_yy_xz[i] = t_0_xxz_y_xyz[i] - rcd_y[i] * t_0_xxz_y_xz[i];

        t_0_xxz_yy_xy[i] = t_0_xxz_y_xyy[i] - rcd_y[i] * t_0_xxz_y_xy[i];

        t_0_xxz_yy_xx[i] = t_0_xxz_y_xxy[i] - rcd_y[i] * t_0_xxz_y_xx[i];

        t_0_xxz_xz_zz[i] = t_0_xxz_x_zzz[i] - rcd_z[i] * t_0_xxz_x_zz[i];

        t_0_xxz_xz_yz[i] = t_0_xxz_x_yzz[i] - rcd_z[i] * t_0_xxz_x_yz[i];

        t_0_xxz_xz_yy[i] = t_0_xxz_x_yyz[i] - rcd_z[i] * t_0_xxz_x_yy[i];

        t_0_xxz_xz_xz[i] = t_0_xxz_x_xzz[i] - rcd_z[i] * t_0_xxz_x_xz[i];

        t_0_xxz_xz_xy[i] = t_0_xxz_x_xyz[i] - rcd_z[i] * t_0_xxz_x_xy[i];

        t_0_xxz_xz_xx[i] = t_0_xxz_x_xxz[i] - rcd_z[i] * t_0_xxz_x_xx[i];

        t_0_xxz_xy_zz[i] = t_0_xxz_x_yzz[i] - rcd_y[i] * t_0_xxz_x_zz[i];

        t_0_xxz_xy_yz[i] = t_0_xxz_x_yyz[i] - rcd_y[i] * t_0_xxz_x_yz[i];

        t_0_xxz_xy_yy[i] = t_0_xxz_x_yyy[i] - rcd_y[i] * t_0_xxz_x_yy[i];

        t_0_xxz_xy_xz[i] = t_0_xxz_x_xyz[i] - rcd_y[i] * t_0_xxz_x_xz[i];

        t_0_xxz_xy_xy[i] = t_0_xxz_x_xyy[i] - rcd_y[i] * t_0_xxz_x_xy[i];

        t_0_xxz_xy_xx[i] = t_0_xxz_x_xxy[i] - rcd_y[i] * t_0_xxz_x_xx[i];

        t_0_xxz_xx_zz[i] = t_0_xxz_x_xzz[i] - rcd_x[i] * t_0_xxz_x_zz[i];

        t_0_xxz_xx_yz[i] = t_0_xxz_x_xyz[i] - rcd_x[i] * t_0_xxz_x_yz[i];

        t_0_xxz_xx_yy[i] = t_0_xxz_x_xyy[i] - rcd_x[i] * t_0_xxz_x_yy[i];

        t_0_xxz_xx_xz[i] = t_0_xxz_x_xxz[i] - rcd_x[i] * t_0_xxz_x_xz[i];

        t_0_xxz_xx_xy[i] = t_0_xxz_x_xxy[i] - rcd_x[i] * t_0_xxz_x_xy[i];

        t_0_xxz_xx_xx[i] = t_0_xxz_x_xxx[i] - rcd_x[i] * t_0_xxz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxy_x_xx, t_0_xxy_x_xxx, t_0_xxy_x_xxy,\
                           t_0_xxy_x_xxz, t_0_xxy_x_xy, t_0_xxy_x_xyy, t_0_xxy_x_xyz,\
                           t_0_xxy_x_xz, t_0_xxy_x_xzz, t_0_xxy_x_yy, t_0_xxy_x_yyy,\
                           t_0_xxy_x_yyz, t_0_xxy_x_yz, t_0_xxy_x_yzz, t_0_xxy_x_zz,\
                           t_0_xxy_x_zzz, t_0_xxy_xx_xx, t_0_xxy_xx_xy, t_0_xxy_xx_xz,\
                           t_0_xxy_xx_yy, t_0_xxy_xx_yz, t_0_xxy_xx_zz, t_0_xxy_xy_xx,\
                           t_0_xxy_xy_xy, t_0_xxy_xy_xz, t_0_xxy_xy_yy, t_0_xxy_xy_yz,\
                           t_0_xxy_xy_zz, t_0_xxy_xz_xx, t_0_xxy_xz_xy, t_0_xxy_xz_xz,\
                           t_0_xxy_xz_yy, t_0_xxy_xz_yz, t_0_xxy_xz_zz, t_0_xxy_y_xx,\
                           t_0_xxy_y_xxy, t_0_xxy_y_xxz, t_0_xxy_y_xy, t_0_xxy_y_xyy,\
                           t_0_xxy_y_xyz, t_0_xxy_y_xz, t_0_xxy_y_xzz, t_0_xxy_y_yy,\
                           t_0_xxy_y_yyy, t_0_xxy_y_yyz, t_0_xxy_y_yz, t_0_xxy_y_yzz,\
                           t_0_xxy_y_zz, t_0_xxy_y_zzz, t_0_xxy_yy_xx, t_0_xxy_yy_xy,\
                           t_0_xxy_yy_xz, t_0_xxy_yy_yy, t_0_xxy_yy_yz, t_0_xxy_yy_zz,\
                           t_0_xxy_yz_xx, t_0_xxy_yz_xy, t_0_xxy_yz_xz, t_0_xxy_yz_yy,\
                           t_0_xxy_yz_yz, t_0_xxy_yz_zz, t_0_xxy_z_xx, t_0_xxy_z_xxz,\
                           t_0_xxy_z_xy, t_0_xxy_z_xyz, t_0_xxy_z_xz, t_0_xxy_z_xzz,\
                           t_0_xxy_z_yy, t_0_xxy_z_yyz, t_0_xxy_z_yz, t_0_xxy_z_yzz,\
                           t_0_xxy_z_zz, t_0_xxy_z_zzz, t_0_xxy_zz_xx, t_0_xxy_zz_xy,\
                           t_0_xxy_zz_xz, t_0_xxy_zz_yy, t_0_xxy_zz_yz, t_0_xxy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxy_zz_zz[i] = t_0_xxy_z_zzz[i] - rcd_z[i] * t_0_xxy_z_zz[i];

        t_0_xxy_zz_yz[i] = t_0_xxy_z_yzz[i] - rcd_z[i] * t_0_xxy_z_yz[i];

        t_0_xxy_zz_yy[i] = t_0_xxy_z_yyz[i] - rcd_z[i] * t_0_xxy_z_yy[i];

        t_0_xxy_zz_xz[i] = t_0_xxy_z_xzz[i] - rcd_z[i] * t_0_xxy_z_xz[i];

        t_0_xxy_zz_xy[i] = t_0_xxy_z_xyz[i] - rcd_z[i] * t_0_xxy_z_xy[i];

        t_0_xxy_zz_xx[i] = t_0_xxy_z_xxz[i] - rcd_z[i] * t_0_xxy_z_xx[i];

        t_0_xxy_yz_zz[i] = t_0_xxy_y_zzz[i] - rcd_z[i] * t_0_xxy_y_zz[i];

        t_0_xxy_yz_yz[i] = t_0_xxy_y_yzz[i] - rcd_z[i] * t_0_xxy_y_yz[i];

        t_0_xxy_yz_yy[i] = t_0_xxy_y_yyz[i] - rcd_z[i] * t_0_xxy_y_yy[i];

        t_0_xxy_yz_xz[i] = t_0_xxy_y_xzz[i] - rcd_z[i] * t_0_xxy_y_xz[i];

        t_0_xxy_yz_xy[i] = t_0_xxy_y_xyz[i] - rcd_z[i] * t_0_xxy_y_xy[i];

        t_0_xxy_yz_xx[i] = t_0_xxy_y_xxz[i] - rcd_z[i] * t_0_xxy_y_xx[i];

        t_0_xxy_yy_zz[i] = t_0_xxy_y_yzz[i] - rcd_y[i] * t_0_xxy_y_zz[i];

        t_0_xxy_yy_yz[i] = t_0_xxy_y_yyz[i] - rcd_y[i] * t_0_xxy_y_yz[i];

        t_0_xxy_yy_yy[i] = t_0_xxy_y_yyy[i] - rcd_y[i] * t_0_xxy_y_yy[i];

        t_0_xxy_yy_xz[i] = t_0_xxy_y_xyz[i] - rcd_y[i] * t_0_xxy_y_xz[i];

        t_0_xxy_yy_xy[i] = t_0_xxy_y_xyy[i] - rcd_y[i] * t_0_xxy_y_xy[i];

        t_0_xxy_yy_xx[i] = t_0_xxy_y_xxy[i] - rcd_y[i] * t_0_xxy_y_xx[i];

        t_0_xxy_xz_zz[i] = t_0_xxy_x_zzz[i] - rcd_z[i] * t_0_xxy_x_zz[i];

        t_0_xxy_xz_yz[i] = t_0_xxy_x_yzz[i] - rcd_z[i] * t_0_xxy_x_yz[i];

        t_0_xxy_xz_yy[i] = t_0_xxy_x_yyz[i] - rcd_z[i] * t_0_xxy_x_yy[i];

        t_0_xxy_xz_xz[i] = t_0_xxy_x_xzz[i] - rcd_z[i] * t_0_xxy_x_xz[i];

        t_0_xxy_xz_xy[i] = t_0_xxy_x_xyz[i] - rcd_z[i] * t_0_xxy_x_xy[i];

        t_0_xxy_xz_xx[i] = t_0_xxy_x_xxz[i] - rcd_z[i] * t_0_xxy_x_xx[i];

        t_0_xxy_xy_zz[i] = t_0_xxy_x_yzz[i] - rcd_y[i] * t_0_xxy_x_zz[i];

        t_0_xxy_xy_yz[i] = t_0_xxy_x_yyz[i] - rcd_y[i] * t_0_xxy_x_yz[i];

        t_0_xxy_xy_yy[i] = t_0_xxy_x_yyy[i] - rcd_y[i] * t_0_xxy_x_yy[i];

        t_0_xxy_xy_xz[i] = t_0_xxy_x_xyz[i] - rcd_y[i] * t_0_xxy_x_xz[i];

        t_0_xxy_xy_xy[i] = t_0_xxy_x_xyy[i] - rcd_y[i] * t_0_xxy_x_xy[i];

        t_0_xxy_xy_xx[i] = t_0_xxy_x_xxy[i] - rcd_y[i] * t_0_xxy_x_xx[i];

        t_0_xxy_xx_zz[i] = t_0_xxy_x_xzz[i] - rcd_x[i] * t_0_xxy_x_zz[i];

        t_0_xxy_xx_yz[i] = t_0_xxy_x_xyz[i] - rcd_x[i] * t_0_xxy_x_yz[i];

        t_0_xxy_xx_yy[i] = t_0_xxy_x_xyy[i] - rcd_x[i] * t_0_xxy_x_yy[i];

        t_0_xxy_xx_xz[i] = t_0_xxy_x_xxz[i] - rcd_x[i] * t_0_xxy_x_xz[i];

        t_0_xxy_xx_xy[i] = t_0_xxy_x_xxy[i] - rcd_x[i] * t_0_xxy_x_xy[i];

        t_0_xxy_xx_xx[i] = t_0_xxy_x_xxx[i] - rcd_x[i] * t_0_xxy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_x_xx, t_0_xxx_x_xxx, t_0_xxx_x_xxy,\
                           t_0_xxx_x_xxz, t_0_xxx_x_xy, t_0_xxx_x_xyy, t_0_xxx_x_xyz,\
                           t_0_xxx_x_xz, t_0_xxx_x_xzz, t_0_xxx_x_yy, t_0_xxx_x_yyy,\
                           t_0_xxx_x_yyz, t_0_xxx_x_yz, t_0_xxx_x_yzz, t_0_xxx_x_zz,\
                           t_0_xxx_x_zzz, t_0_xxx_xx_xx, t_0_xxx_xx_xy, t_0_xxx_xx_xz,\
                           t_0_xxx_xx_yy, t_0_xxx_xx_yz, t_0_xxx_xx_zz, t_0_xxx_xy_xx,\
                           t_0_xxx_xy_xy, t_0_xxx_xy_xz, t_0_xxx_xy_yy, t_0_xxx_xy_yz,\
                           t_0_xxx_xy_zz, t_0_xxx_xz_xx, t_0_xxx_xz_xy, t_0_xxx_xz_xz,\
                           t_0_xxx_xz_yy, t_0_xxx_xz_yz, t_0_xxx_xz_zz, t_0_xxx_y_xx,\
                           t_0_xxx_y_xxy, t_0_xxx_y_xxz, t_0_xxx_y_xy, t_0_xxx_y_xyy,\
                           t_0_xxx_y_xyz, t_0_xxx_y_xz, t_0_xxx_y_xzz, t_0_xxx_y_yy,\
                           t_0_xxx_y_yyy, t_0_xxx_y_yyz, t_0_xxx_y_yz, t_0_xxx_y_yzz,\
                           t_0_xxx_y_zz, t_0_xxx_y_zzz, t_0_xxx_yy_xx, t_0_xxx_yy_xy,\
                           t_0_xxx_yy_xz, t_0_xxx_yy_yy, t_0_xxx_yy_yz, t_0_xxx_yy_zz,\
                           t_0_xxx_yz_xx, t_0_xxx_yz_xy, t_0_xxx_yz_xz, t_0_xxx_yz_yy,\
                           t_0_xxx_yz_yz, t_0_xxx_yz_zz, t_0_xxx_z_xx, t_0_xxx_z_xxz,\
                           t_0_xxx_z_xy, t_0_xxx_z_xyz, t_0_xxx_z_xz, t_0_xxx_z_xzz,\
                           t_0_xxx_z_yy, t_0_xxx_z_yyz, t_0_xxx_z_yz, t_0_xxx_z_yzz,\
                           t_0_xxx_z_zz, t_0_xxx_z_zzz, t_0_xxx_zz_xx, t_0_xxx_zz_xy,\
                           t_0_xxx_zz_xz, t_0_xxx_zz_yy, t_0_xxx_zz_yz, t_0_xxx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxx_zz_zz[i] = t_0_xxx_z_zzz[i] - rcd_z[i] * t_0_xxx_z_zz[i];

        t_0_xxx_zz_yz[i] = t_0_xxx_z_yzz[i] - rcd_z[i] * t_0_xxx_z_yz[i];

        t_0_xxx_zz_yy[i] = t_0_xxx_z_yyz[i] - rcd_z[i] * t_0_xxx_z_yy[i];

        t_0_xxx_zz_xz[i] = t_0_xxx_z_xzz[i] - rcd_z[i] * t_0_xxx_z_xz[i];

        t_0_xxx_zz_xy[i] = t_0_xxx_z_xyz[i] - rcd_z[i] * t_0_xxx_z_xy[i];

        t_0_xxx_zz_xx[i] = t_0_xxx_z_xxz[i] - rcd_z[i] * t_0_xxx_z_xx[i];

        t_0_xxx_yz_zz[i] = t_0_xxx_y_zzz[i] - rcd_z[i] * t_0_xxx_y_zz[i];

        t_0_xxx_yz_yz[i] = t_0_xxx_y_yzz[i] - rcd_z[i] * t_0_xxx_y_yz[i];

        t_0_xxx_yz_yy[i] = t_0_xxx_y_yyz[i] - rcd_z[i] * t_0_xxx_y_yy[i];

        t_0_xxx_yz_xz[i] = t_0_xxx_y_xzz[i] - rcd_z[i] * t_0_xxx_y_xz[i];

        t_0_xxx_yz_xy[i] = t_0_xxx_y_xyz[i] - rcd_z[i] * t_0_xxx_y_xy[i];

        t_0_xxx_yz_xx[i] = t_0_xxx_y_xxz[i] - rcd_z[i] * t_0_xxx_y_xx[i];

        t_0_xxx_yy_zz[i] = t_0_xxx_y_yzz[i] - rcd_y[i] * t_0_xxx_y_zz[i];

        t_0_xxx_yy_yz[i] = t_0_xxx_y_yyz[i] - rcd_y[i] * t_0_xxx_y_yz[i];

        t_0_xxx_yy_yy[i] = t_0_xxx_y_yyy[i] - rcd_y[i] * t_0_xxx_y_yy[i];

        t_0_xxx_yy_xz[i] = t_0_xxx_y_xyz[i] - rcd_y[i] * t_0_xxx_y_xz[i];

        t_0_xxx_yy_xy[i] = t_0_xxx_y_xyy[i] - rcd_y[i] * t_0_xxx_y_xy[i];

        t_0_xxx_yy_xx[i] = t_0_xxx_y_xxy[i] - rcd_y[i] * t_0_xxx_y_xx[i];

        t_0_xxx_xz_zz[i] = t_0_xxx_x_zzz[i] - rcd_z[i] * t_0_xxx_x_zz[i];

        t_0_xxx_xz_yz[i] = t_0_xxx_x_yzz[i] - rcd_z[i] * t_0_xxx_x_yz[i];

        t_0_xxx_xz_yy[i] = t_0_xxx_x_yyz[i] - rcd_z[i] * t_0_xxx_x_yy[i];

        t_0_xxx_xz_xz[i] = t_0_xxx_x_xzz[i] - rcd_z[i] * t_0_xxx_x_xz[i];

        t_0_xxx_xz_xy[i] = t_0_xxx_x_xyz[i] - rcd_z[i] * t_0_xxx_x_xy[i];

        t_0_xxx_xz_xx[i] = t_0_xxx_x_xxz[i] - rcd_z[i] * t_0_xxx_x_xx[i];

        t_0_xxx_xy_zz[i] = t_0_xxx_x_yzz[i] - rcd_y[i] * t_0_xxx_x_zz[i];

        t_0_xxx_xy_yz[i] = t_0_xxx_x_yyz[i] - rcd_y[i] * t_0_xxx_x_yz[i];

        t_0_xxx_xy_yy[i] = t_0_xxx_x_yyy[i] - rcd_y[i] * t_0_xxx_x_yy[i];

        t_0_xxx_xy_xz[i] = t_0_xxx_x_xyz[i] - rcd_y[i] * t_0_xxx_x_xz[i];

        t_0_xxx_xy_xy[i] = t_0_xxx_x_xyy[i] - rcd_y[i] * t_0_xxx_x_xy[i];

        t_0_xxx_xy_xx[i] = t_0_xxx_x_xxy[i] - rcd_y[i] * t_0_xxx_x_xx[i];

        t_0_xxx_xx_zz[i] = t_0_xxx_x_xzz[i] - rcd_x[i] * t_0_xxx_x_zz[i];

        t_0_xxx_xx_yz[i] = t_0_xxx_x_xyz[i] - rcd_x[i] * t_0_xxx_x_yz[i];

        t_0_xxx_xx_yy[i] = t_0_xxx_x_xyy[i] - rcd_x[i] * t_0_xxx_x_yy[i];

        t_0_xxx_xx_xz[i] = t_0_xxx_x_xxz[i] - rcd_x[i] * t_0_xxx_x_xz[i];

        t_0_xxx_xx_xy[i] = t_0_xxx_x_xxy[i] - rcd_x[i] * t_0_xxx_x_xy[i];

        t_0_xxx_xx_xx[i] = t_0_xxx_x_xxx[i] - rcd_x[i] * t_0_xxx_x_xx[i];
    }
}


} // derirec namespace
