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
compHostHRRForPFDD_V0(      BufferHostXY<T>&      intsBufferPFDD,
                      const BufferHostX<int32_t>& intsIndexesPFDD,
                      const BufferHostXY<T>&      intsBufferSFDD,
                      const BufferHostX<int32_t>& intsIndexesSFDD,
                      const BufferHostXY<T>&      intsBufferSGDD,
                      const BufferHostX<int32_t>& intsIndexesSGDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PFDD) integral components

    t_z_zzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(0));

    t_z_zzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(1));

    t_z_zzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(2));

    t_z_zzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(3));

    t_z_zzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(4));

    t_z_zzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(5));

    t_z_zzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(6));

    t_z_zzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(7));

    t_z_zzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(8));

    t_z_zzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(9));

    t_z_zzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(10));

    t_z_zzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(11));

    t_z_zzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(12));

    t_z_zzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(13));

    t_z_zzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(14));

    t_z_zzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(15));

    t_z_zzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(16));

    t_z_zzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(17));

    t_z_zzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(18));

    t_z_zzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(19));

    t_z_zzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(20));

    t_z_zzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(21));

    t_z_zzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(22));

    t_z_zzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(23));

    t_z_zzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(24));

    t_z_zzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(25));

    t_z_zzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(26));

    t_z_zzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(27));

    t_z_zzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(28));

    t_z_zzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(29));

    t_z_zzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(30));

    t_z_zzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(31));

    t_z_zzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(32));

    t_z_zzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(33));

    t_z_zzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(34));

    t_z_zzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(35));

    t_z_yzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(36));

    t_z_yzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(37));

    t_z_yzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(38));

    t_z_yzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(39));

    t_z_yzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(40));

    t_z_yzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(41));

    t_z_yzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(42));

    t_z_yzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(43));

    t_z_yzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(44));

    t_z_yzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(45));

    t_z_yzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(46));

    t_z_yzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(47));

    t_z_yzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(48));

    t_z_yzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(49));

    t_z_yzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(50));

    t_z_yzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(51));

    t_z_yzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(52));

    t_z_yzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(53));

    t_z_yzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(54));

    t_z_yzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(55));

    t_z_yzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(56));

    t_z_yzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(57));

    t_z_yzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(58));

    t_z_yzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(59));

    t_z_yzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(60));

    t_z_yzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(61));

    t_z_yzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(62));

    t_z_yzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(63));

    t_z_yzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(64));

    t_z_yzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(65));

    t_z_yzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(66));

    t_z_yzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(67));

    t_z_yzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(68));

    t_z_yzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(69));

    t_z_yzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(70));

    t_z_yzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(71));

    t_z_yyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(72));

    t_z_yyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(73));

    t_z_yyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(74));

    t_z_yyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(75));

    t_z_yyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(76));

    t_z_yyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(77));

    t_z_yyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(78));

    t_z_yyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(79));

    t_z_yyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(80));

    t_z_yyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(81));

    t_z_yyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(82));

    t_z_yyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(83));

    t_z_yyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(84));

    t_z_yyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(85));

    t_z_yyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(86));

    t_z_yyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(87));

    t_z_yyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(88));

    t_z_yyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(89));

    t_z_yyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(90));

    t_z_yyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(91));

    t_z_yyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(92));

    t_z_yyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(93));

    t_z_yyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(94));

    t_z_yyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(95));

    t_z_yyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(96));

    t_z_yyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(97));

    t_z_yyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(98));

    t_z_yyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(99));

    t_z_yyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(100));

    t_z_yyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(101));

    t_z_yyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(102));

    t_z_yyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(103));

    t_z_yyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(104));

    t_z_yyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(105));

    t_z_yyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(106));

    t_z_yyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(107));

    t_z_xzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(108));

    t_z_xzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(109));

    t_z_xzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(110));

    t_z_xzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(111));

    t_z_xzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(112));

    t_z_xzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(113));

    t_z_xzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(114));

    t_z_xzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(115));

    t_z_xzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(116));

    t_z_xzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(117));

    t_z_xzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(118));

    t_z_xzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(119));

    t_z_xzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(120));

    t_z_xzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(121));

    t_z_xzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(122));

    t_z_xzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(123));

    t_z_xzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(124));

    t_z_xzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(125));

    t_z_xzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(126));

    t_z_xzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(127));

    t_z_xzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(128));

    t_z_xzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(129));

    t_z_xzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(130));

    t_z_xzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(131));

    t_z_xzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(132));

    t_z_xzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(133));

    t_z_xzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(134));

    t_z_xzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(135));

    t_z_xzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(136));

    t_z_xzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(137));

    t_z_xzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(138));

    t_z_xzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(139));

    t_z_xzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(140));

    t_z_xzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(141));

    t_z_xzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(142));

    t_z_xzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(143));

    t_z_xyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(144));

    t_z_xyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(145));

    t_z_xyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(146));

    t_z_xyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(147));

    t_z_xyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(148));

    t_z_xyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(149));

    t_z_xyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(150));

    t_z_xyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(151));

    t_z_xyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(152));

    t_z_xyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(153));

    t_z_xyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(154));

    t_z_xyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(155));

    t_z_xyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(156));

    t_z_xyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(157));

    t_z_xyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(158));

    t_z_xyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(159));

    t_z_xyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(160));

    t_z_xyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(161));

    t_z_xyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(162));

    t_z_xyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(163));

    t_z_xyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(164));

    t_z_xyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(165));

    t_z_xyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(166));

    t_z_xyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(167));

    t_z_xyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(168));

    t_z_xyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(169));

    t_z_xyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(170));

    t_z_xyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(171));

    t_z_xyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(172));

    t_z_xyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(173));

    t_z_xyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(174));

    t_z_xyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(175));

    t_z_xyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(176));

    t_z_xyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(177));

    t_z_xyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(178));

    t_z_xyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(179));

    t_z_xxz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(180));

    t_z_xxz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(181));

    t_z_xxz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(182));

    t_z_xxz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(183));

    t_z_xxz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(184));

    t_z_xxz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(185));

    t_z_xxz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(186));

    t_z_xxz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(187));

    t_z_xxz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(188));

    t_z_xxz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(189));

    t_z_xxz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(190));

    t_z_xxz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(191));

    t_z_xxz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(192));

    t_z_xxz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(193));

    t_z_xxz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(194));

    t_z_xxz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(195));

    t_z_xxz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(196));

    t_z_xxz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(197));

    t_z_xxz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(198));

    t_z_xxz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(199));

    t_z_xxz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(200));

    t_z_xxz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(201));

    t_z_xxz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(202));

    t_z_xxz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(203));

    t_z_xxz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(204));

    t_z_xxz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(205));

    t_z_xxz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(206));

    t_z_xxz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(207));

    t_z_xxz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(208));

    t_z_xxz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(209));

    t_z_xxz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(210));

    t_z_xxz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(211));

    t_z_xxz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(212));

    t_z_xxz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(213));

    t_z_xxz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(214));

    t_z_xxz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(215));

    t_y_zzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(216));

    t_y_zzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(217));

    t_y_zzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(218));

    t_y_zzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(219));

    t_y_zzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(220));

    t_y_zzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(221));

    t_y_zzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(222));

    t_y_zzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(223));

    t_y_zzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(224));

    t_y_zzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(225));

    t_y_zzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(226));

    t_y_zzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(227));

    t_y_zzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(228));

    t_y_zzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(229));

    t_y_zzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(230));

    t_y_zzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(231));

    t_y_zzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(232));

    t_y_zzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(233));

    t_y_zzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(234));

    t_y_zzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(235));

    t_y_zzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(236));

    t_y_zzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(237));

    t_y_zzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(238));

    t_y_zzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(239));

    t_y_zzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(240));

    t_y_zzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(241));

    t_y_zzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(242));

    t_y_zzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(243));

    t_y_zzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(244));

    t_y_zzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(245));

    t_y_zzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(246));

    t_y_zzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(247));

    t_y_zzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(248));

    t_y_zzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(249));

    t_y_zzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(250));

    t_y_zzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(251));

    t_y_yzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(252));

    t_y_yzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(253));

    t_y_yzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(254));

    t_y_yzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(255));

    t_y_yzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(256));

    t_y_yzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(257));

    t_y_yzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(258));

    t_y_yzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(259));

    t_y_yzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(260));

    t_y_yzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(261));

    t_y_yzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(262));

    t_y_yzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(263));

    t_y_yzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(264));

    t_y_yzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(265));

    t_y_yzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(266));

    t_y_yzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(267));

    t_y_yzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(268));

    t_y_yzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(269));

    t_y_yzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(270));

    t_y_yzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(271));

    t_y_yzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(272));

    t_y_yzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(273));

    t_y_yzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(274));

    t_y_yzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(275));

    t_y_yzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(276));

    t_y_yzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(277));

    t_y_yzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(278));

    t_y_yzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(279));

    t_y_yzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(280));

    t_y_yzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(281));

    t_y_yzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(282));

    t_y_yzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(283));

    t_y_yzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(284));

    t_y_yzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(285));

    t_y_yzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(286));

    t_y_yzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(287));

    t_y_yyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(288));

    t_y_yyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(289));

    t_y_yyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(290));

    t_y_yyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(291));

    t_y_yyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(292));

    t_y_yyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(293));

    t_y_yyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(294));

    t_y_yyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(295));

    t_y_yyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(296));

    t_y_yyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(297));

    t_y_yyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(298));

    t_y_yyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(299));

    t_y_yyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(300));

    t_y_yyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(301));

    t_y_yyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(302));

    t_y_yyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(303));

    t_y_yyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(304));

    t_y_yyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(305));

    t_y_yyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(306));

    t_y_yyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(307));

    t_y_yyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(308));

    t_y_yyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(309));

    t_y_yyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(310));

    t_y_yyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(311));

    t_y_yyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(312));

    t_y_yyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(313));

    t_y_yyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(314));

    t_y_yyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(315));

    t_y_yyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(316));

    t_y_yyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(317));

    t_y_yyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(318));

    t_y_yyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(319));

    t_y_yyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(320));

    t_y_yyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(321));

    t_y_yyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(322));

    t_y_yyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(323));

    t_y_yyy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(324));

    t_y_yyy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(325));

    t_y_yyy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(326));

    t_y_yyy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(327));

    t_y_yyy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(328));

    t_y_yyy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(329));

    t_y_yyy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(330));

    t_y_yyy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(331));

    t_y_yyy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(332));

    t_y_yyy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(333));

    t_y_yyy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(334));

    t_y_yyy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(335));

    t_y_yyy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(336));

    t_y_yyy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(337));

    t_y_yyy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(338));

    t_y_yyy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(339));

    t_y_yyy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(340));

    t_y_yyy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(341));

    t_y_yyy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(342));

    t_y_yyy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(343));

    t_y_yyy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(344));

    t_y_yyy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(345));

    t_y_yyy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(346));

    t_y_yyy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(347));

    t_y_yyy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(348));

    t_y_yyy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(349));

    t_y_yyy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(350));

    t_y_yyy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(351));

    t_y_yyy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(352));

    t_y_yyy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(353));

    t_y_yyy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(354));

    t_y_yyy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(355));

    t_y_yyy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(356));

    t_y_yyy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(357));

    t_y_yyy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(358));

    t_y_yyy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(359));

    t_y_xzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(360));

    t_y_xzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(361));

    t_y_xzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(362));

    t_y_xzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(363));

    t_y_xzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(364));

    t_y_xzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(365));

    t_y_xzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(366));

    t_y_xzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(367));

    t_y_xzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(368));

    t_y_xzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(369));

    t_y_xzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(370));

    t_y_xzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(371));

    t_y_xzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(372));

    t_y_xzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(373));

    t_y_xzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(374));

    t_y_xzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(375));

    t_y_xzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(376));

    t_y_xzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(377));

    t_y_xzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(378));

    t_y_xzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(379));

    t_y_xzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(380));

    t_y_xzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(381));

    t_y_xzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(382));

    t_y_xzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(383));

    t_y_xzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(384));

    t_y_xzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(385));

    t_y_xzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(386));

    t_y_xzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(387));

    t_y_xzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(388));

    t_y_xzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(389));

    t_y_xzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(390));

    t_y_xzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(391));

    t_y_xzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(392));

    t_y_xzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(393));

    t_y_xzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(394));

    t_y_xzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(395));

    t_y_xyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(396));

    t_y_xyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(397));

    t_y_xyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(398));

    t_y_xyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(399));

    t_y_xyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(400));

    t_y_xyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(401));

    t_y_xyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(402));

    t_y_xyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(403));

    t_y_xyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(404));

    t_y_xyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(405));

    t_y_xyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(406));

    t_y_xyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(407));

    t_y_xyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(408));

    t_y_xyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(409));

    t_y_xyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(410));

    t_y_xyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(411));

    t_y_xyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(412));

    t_y_xyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(413));

    t_y_xyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(414));

    t_y_xyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(415));

    t_y_xyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(416));

    t_y_xyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(417));

    t_y_xyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(418));

    t_y_xyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(419));

    t_y_xyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(420));

    t_y_xyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(421));

    t_y_xyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(422));

    t_y_xyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(423));

    t_y_xyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(424));

    t_y_xyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(425));

    t_y_xyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(426));

    t_y_xyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(427));

    t_y_xyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(428));

    t_y_xyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(429));

    t_y_xyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(430));

    t_y_xyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(431));

    t_y_xyy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(432));

    t_y_xyy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(433));

    t_y_xyy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(434));

    t_y_xyy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(435));

    t_y_xyy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(436));

    t_y_xyy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(437));

    t_y_xyy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(438));

    t_y_xyy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(439));

    t_y_xyy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(440));

    t_y_xyy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(441));

    t_y_xyy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(442));

    t_y_xyy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(443));

    t_y_xyy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(444));

    t_y_xyy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(445));

    t_y_xyy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(446));

    t_y_xyy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(447));

    t_y_xyy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(448));

    t_y_xyy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(449));

    t_y_xyy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(450));

    t_y_xyy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(451));

    t_y_xyy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(452));

    t_y_xyy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(453));

    t_y_xyy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(454));

    t_y_xyy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(455));

    t_y_xyy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(456));

    t_y_xyy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(457));

    t_y_xyy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(458));

    t_y_xyy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(459));

    t_y_xyy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(460));

    t_y_xyy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(461));

    t_y_xyy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(462));

    t_y_xyy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(463));

    t_y_xyy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(464));

    t_y_xyy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(465));

    t_y_xyy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(466));

    t_y_xyy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(467));

    t_y_xxz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(468));

    t_y_xxz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(469));

    t_y_xxz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(470));

    t_y_xxz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(471));

    t_y_xxz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(472));

    t_y_xxz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(473));

    t_y_xxz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(474));

    t_y_xxz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(475));

    t_y_xxz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(476));

    t_y_xxz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(477));

    t_y_xxz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(478));

    t_y_xxz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(479));

    t_y_xxz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(480));

    t_y_xxz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(481));

    t_y_xxz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(482));

    t_y_xxz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(483));

    t_y_xxz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(484));

    t_y_xxz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(485));

    t_y_xxz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(486));

    t_y_xxz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(487));

    t_y_xxz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(488));

    t_y_xxz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(489));

    t_y_xxz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(490));

    t_y_xxz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(491));

    t_y_xxz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(492));

    t_y_xxz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(493));

    t_y_xxz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(494));

    t_y_xxz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(495));

    t_y_xxz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(496));

    t_y_xxz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(497));

    t_y_xxz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(498));

    t_y_xxz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(499));

    t_y_xxz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(500));

    t_y_xxz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(501));

    t_y_xxz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(502));

    t_y_xxz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(503));

    t_y_xxy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(504));

    t_y_xxy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(505));

    t_y_xxy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(506));

    t_y_xxy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(507));

    t_y_xxy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(508));

    t_y_xxy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(509));

    t_y_xxy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(510));

    t_y_xxy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(511));

    t_y_xxy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(512));

    t_y_xxy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(513));

    t_y_xxy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(514));

    t_y_xxy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(515));

    t_y_xxy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(516));

    t_y_xxy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(517));

    t_y_xxy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(518));

    t_y_xxy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(519));

    t_y_xxy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(520));

    t_y_xxy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(521));

    t_y_xxy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(522));

    t_y_xxy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(523));

    t_y_xxy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(524));

    t_y_xxy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(525));

    t_y_xxy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(526));

    t_y_xxy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(527));

    t_y_xxy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(528));

    t_y_xxy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(529));

    t_y_xxy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(530));

    t_y_xxy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(531));

    t_y_xxy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(532));

    t_y_xxy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(533));

    t_y_xxy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(534));

    t_y_xxy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(535));

    t_y_xxy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(536));

    t_y_xxy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(537));

    t_y_xxy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(538));

    t_y_xxy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(539));

    t_x_zzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(540));

    t_x_zzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(541));

    t_x_zzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(542));

    t_x_zzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(543));

    t_x_zzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(544));

    t_x_zzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(545));

    t_x_zzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(546));

    t_x_zzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(547));

    t_x_zzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(548));

    t_x_zzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(549));

    t_x_zzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(550));

    t_x_zzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(551));

    t_x_zzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(552));

    t_x_zzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(553));

    t_x_zzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(554));

    t_x_zzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(555));

    t_x_zzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(556));

    t_x_zzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(557));

    t_x_zzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(558));

    t_x_zzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(559));

    t_x_zzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(560));

    t_x_zzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(561));

    t_x_zzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(562));

    t_x_zzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(563));

    t_x_zzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(564));

    t_x_zzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(565));

    t_x_zzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(566));

    t_x_zzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(567));

    t_x_zzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(568));

    t_x_zzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(569));

    t_x_zzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(570));

    t_x_zzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(571));

    t_x_zzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(572));

    t_x_zzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(573));

    t_x_zzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(574));

    t_x_zzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(575));

    t_x_yzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(576));

    t_x_yzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(577));

    t_x_yzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(578));

    t_x_yzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(579));

    t_x_yzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(580));

    t_x_yzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(581));

    t_x_yzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(582));

    t_x_yzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(583));

    t_x_yzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(584));

    t_x_yzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(585));

    t_x_yzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(586));

    t_x_yzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(587));

    t_x_yzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(588));

    t_x_yzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(589));

    t_x_yzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(590));

    t_x_yzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(591));

    t_x_yzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(592));

    t_x_yzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(593));

    t_x_yzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(594));

    t_x_yzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(595));

    t_x_yzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(596));

    t_x_yzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(597));

    t_x_yzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(598));

    t_x_yzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(599));

    t_x_yzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(600));

    t_x_yzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(601));

    t_x_yzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(602));

    t_x_yzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(603));

    t_x_yzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(604));

    t_x_yzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(605));

    t_x_yzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(606));

    t_x_yzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(607));

    t_x_yzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(608));

    t_x_yzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(609));

    t_x_yzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(610));

    t_x_yzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(611));

    t_x_yyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(612));

    t_x_yyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(613));

    t_x_yyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(614));

    t_x_yyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(615));

    t_x_yyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(616));

    t_x_yyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(617));

    t_x_yyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(618));

    t_x_yyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(619));

    t_x_yyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(620));

    t_x_yyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(621));

    t_x_yyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(622));

    t_x_yyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(623));

    t_x_yyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(624));

    t_x_yyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(625));

    t_x_yyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(626));

    t_x_yyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(627));

    t_x_yyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(628));

    t_x_yyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(629));

    t_x_yyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(630));

    t_x_yyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(631));

    t_x_yyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(632));

    t_x_yyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(633));

    t_x_yyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(634));

    t_x_yyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(635));

    t_x_yyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(636));

    t_x_yyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(637));

    t_x_yyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(638));

    t_x_yyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(639));

    t_x_yyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(640));

    t_x_yyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(641));

    t_x_yyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(642));

    t_x_yyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(643));

    t_x_yyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(644));

    t_x_yyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(645));

    t_x_yyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(646));

    t_x_yyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(647));

    t_x_yyy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(648));

    t_x_yyy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(649));

    t_x_yyy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(650));

    t_x_yyy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(651));

    t_x_yyy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(652));

    t_x_yyy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(653));

    t_x_yyy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(654));

    t_x_yyy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(655));

    t_x_yyy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(656));

    t_x_yyy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(657));

    t_x_yyy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(658));

    t_x_yyy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(659));

    t_x_yyy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(660));

    t_x_yyy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(661));

    t_x_yyy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(662));

    t_x_yyy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(663));

    t_x_yyy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(664));

    t_x_yyy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(665));

    t_x_yyy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(666));

    t_x_yyy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(667));

    t_x_yyy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(668));

    t_x_yyy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(669));

    t_x_yyy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(670));

    t_x_yyy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(671));

    t_x_yyy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(672));

    t_x_yyy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(673));

    t_x_yyy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(674));

    t_x_yyy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(675));

    t_x_yyy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(676));

    t_x_yyy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(677));

    t_x_yyy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(678));

    t_x_yyy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(679));

    t_x_yyy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(680));

    t_x_yyy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(681));

    t_x_yyy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(682));

    t_x_yyy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(683));

    t_x_xzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(684));

    t_x_xzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(685));

    t_x_xzz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(686));

    t_x_xzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(687));

    t_x_xzz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(688));

    t_x_xzz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(689));

    t_x_xzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(690));

    t_x_xzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(691));

    t_x_xzz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(692));

    t_x_xzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(693));

    t_x_xzz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(694));

    t_x_xzz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(695));

    t_x_xzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(696));

    t_x_xzz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(697));

    t_x_xzz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(698));

    t_x_xzz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(699));

    t_x_xzz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(700));

    t_x_xzz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(701));

    t_x_xzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(702));

    t_x_xzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(703));

    t_x_xzz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(704));

    t_x_xzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(705));

    t_x_xzz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(706));

    t_x_xzz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(707));

    t_x_xzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(708));

    t_x_xzz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(709));

    t_x_xzz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(710));

    t_x_xzz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(711));

    t_x_xzz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(712));

    t_x_xzz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(713));

    t_x_xzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(714));

    t_x_xzz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(715));

    t_x_xzz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(716));

    t_x_xzz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(717));

    t_x_xzz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(718));

    t_x_xzz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(719));

    t_x_xyz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(720));

    t_x_xyz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(721));

    t_x_xyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(722));

    t_x_xyz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(723));

    t_x_xyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(724));

    t_x_xyz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(725));

    t_x_xyz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(726));

    t_x_xyz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(727));

    t_x_xyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(728));

    t_x_xyz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(729));

    t_x_xyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(730));

    t_x_xyz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(731));

    t_x_xyz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(732));

    t_x_xyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(733));

    t_x_xyz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(734));

    t_x_xyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(735));

    t_x_xyz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(736));

    t_x_xyz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(737));

    t_x_xyz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(738));

    t_x_xyz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(739));

    t_x_xyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(740));

    t_x_xyz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(741));

    t_x_xyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(742));

    t_x_xyz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(743));

    t_x_xyz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(744));

    t_x_xyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(745));

    t_x_xyz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(746));

    t_x_xyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(747));

    t_x_xyz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(748));

    t_x_xyz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(749));

    t_x_xyz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(750));

    t_x_xyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(751));

    t_x_xyz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(752));

    t_x_xyz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(753));

    t_x_xyz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(754));

    t_x_xyz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(755));

    t_x_xyy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(756));

    t_x_xyy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(757));

    t_x_xyy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(758));

    t_x_xyy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(759));

    t_x_xyy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(760));

    t_x_xyy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(761));

    t_x_xyy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(762));

    t_x_xyy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(763));

    t_x_xyy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(764));

    t_x_xyy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(765));

    t_x_xyy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(766));

    t_x_xyy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(767));

    t_x_xyy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(768));

    t_x_xyy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(769));

    t_x_xyy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(770));

    t_x_xyy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(771));

    t_x_xyy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(772));

    t_x_xyy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(773));

    t_x_xyy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(774));

    t_x_xyy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(775));

    t_x_xyy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(776));

    t_x_xyy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(777));

    t_x_xyy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(778));

    t_x_xyy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(779));

    t_x_xyy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(780));

    t_x_xyy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(781));

    t_x_xyy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(782));

    t_x_xyy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(783));

    t_x_xyy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(784));

    t_x_xyy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(785));

    t_x_xyy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(786));

    t_x_xyy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(787));

    t_x_xyy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(788));

    t_x_xyy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(789));

    t_x_xyy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(790));

    t_x_xyy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(791));

    t_x_xxz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(792));

    t_x_xxz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(793));

    t_x_xxz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(794));

    t_x_xxz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(795));

    t_x_xxz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(796));

    t_x_xxz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(797));

    t_x_xxz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(798));

    t_x_xxz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(799));

    t_x_xxz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(800));

    t_x_xxz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(801));

    t_x_xxz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(802));

    t_x_xxz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(803));

    t_x_xxz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(804));

    t_x_xxz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(805));

    t_x_xxz_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(806));

    t_x_xxz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(807));

    t_x_xxz_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(808));

    t_x_xxz_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(809));

    t_x_xxz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(810));

    t_x_xxz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(811));

    t_x_xxz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(812));

    t_x_xxz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(813));

    t_x_xxz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(814));

    t_x_xxz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(815));

    t_x_xxz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(816));

    t_x_xxz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(817));

    t_x_xxz_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(818));

    t_x_xxz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(819));

    t_x_xxz_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(820));

    t_x_xxz_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(821));

    t_x_xxz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(822));

    t_x_xxz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(823));

    t_x_xxz_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(824));

    t_x_xxz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(825));

    t_x_xxz_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(826));

    t_x_xxz_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(827));

    t_x_xxy_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(828));

    t_x_xxy_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(829));

    t_x_xxy_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(830));

    t_x_xxy_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(831));

    t_x_xxy_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(832));

    t_x_xxy_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(833));

    t_x_xxy_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(834));

    t_x_xxy_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(835));

    t_x_xxy_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(836));

    t_x_xxy_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(837));

    t_x_xxy_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(838));

    t_x_xxy_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(839));

    t_x_xxy_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(840));

    t_x_xxy_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(841));

    t_x_xxy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(842));

    t_x_xxy_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(843));

    t_x_xxy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(844));

    t_x_xxy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(845));

    t_x_xxy_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(846));

    t_x_xxy_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(847));

    t_x_xxy_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(848));

    t_x_xxy_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(849));

    t_x_xxy_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(850));

    t_x_xxy_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(851));

    t_x_xxy_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(852));

    t_x_xxy_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(853));

    t_x_xxy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(854));

    t_x_xxy_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(855));

    t_x_xxy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(856));

    t_x_xxy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(857));

    t_x_xxy_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(858));

    t_x_xxy_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(859));

    t_x_xxy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(860));

    t_x_xxy_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(861));

    t_x_xxy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(862));

    t_x_xxy_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(863));

    t_x_xxx_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(864));

    t_x_xxx_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(865));

    t_x_xxx_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(866));

    t_x_xxx_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(867));

    t_x_xxx_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(868));

    t_x_xxx_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(869));

    t_x_xxx_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(870));

    t_x_xxx_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(871));

    t_x_xxx_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(872));

    t_x_xxx_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(873));

    t_x_xxx_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(874));

    t_x_xxx_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(875));

    t_x_xxx_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(876));

    t_x_xxx_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(877));

    t_x_xxx_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(878));

    t_x_xxx_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(879));

    t_x_xxx_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(880));

    t_x_xxx_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(881));

    t_x_xxx_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(882));

    t_x_xxx_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(883));

    t_x_xxx_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(884));

    t_x_xxx_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(885));

    t_x_xxx_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(886));

    t_x_xxx_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(887));

    t_x_xxx_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(888));

    t_x_xxx_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(889));

    t_x_xxx_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(890));

    t_x_xxx_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(891));

    t_x_xxx_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(892));

    t_x_xxx_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(893));

    t_x_xxx_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(894));

    t_x_xxx_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(895));

    t_x_xxx_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(896));

    t_x_xxx_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(897));

    t_x_xxx_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(898));

    t_x_xxx_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(899));

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

    // set up (SGDD) integral components

    t_0_zzzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(0));

    t_0_zzzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(1));

    t_0_zzzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(2));

    t_0_zzzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(3));

    t_0_zzzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(4));

    t_0_zzzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(5));

    t_0_zzzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(6));

    t_0_zzzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(7));

    t_0_zzzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(8));

    t_0_zzzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(9));

    t_0_zzzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(10));

    t_0_zzzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(11));

    t_0_zzzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(12));

    t_0_zzzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(13));

    t_0_zzzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(14));

    t_0_zzzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(15));

    t_0_zzzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(16));

    t_0_zzzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(17));

    t_0_zzzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(18));

    t_0_zzzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(19));

    t_0_zzzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(20));

    t_0_zzzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(21));

    t_0_zzzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(22));

    t_0_zzzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(23));

    t_0_zzzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(24));

    t_0_zzzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(25));

    t_0_zzzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(26));

    t_0_zzzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(27));

    t_0_zzzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(28));

    t_0_zzzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(29));

    t_0_zzzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(30));

    t_0_zzzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(31));

    t_0_zzzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(32));

    t_0_zzzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(33));

    t_0_zzzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(34));

    t_0_zzzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(35));

    t_0_yzzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(36));

    t_0_yzzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(37));

    t_0_yzzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(38));

    t_0_yzzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(39));

    t_0_yzzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(40));

    t_0_yzzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(41));

    t_0_yzzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(42));

    t_0_yzzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(43));

    t_0_yzzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(44));

    t_0_yzzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(45));

    t_0_yzzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(46));

    t_0_yzzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(47));

    t_0_yzzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(48));

    t_0_yzzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(49));

    t_0_yzzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(50));

    t_0_yzzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(51));

    t_0_yzzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(52));

    t_0_yzzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(53));

    t_0_yzzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(54));

    t_0_yzzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(55));

    t_0_yzzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(56));

    t_0_yzzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(57));

    t_0_yzzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(58));

    t_0_yzzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(59));

    t_0_yzzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(60));

    t_0_yzzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(61));

    t_0_yzzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(62));

    t_0_yzzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(63));

    t_0_yzzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(64));

    t_0_yzzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(65));

    t_0_yzzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(66));

    t_0_yzzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(67));

    t_0_yzzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(68));

    t_0_yzzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(69));

    t_0_yzzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(70));

    t_0_yzzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(71));

    t_0_yyzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(72));

    t_0_yyzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(73));

    t_0_yyzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(74));

    t_0_yyzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(75));

    t_0_yyzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(76));

    t_0_yyzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(77));

    t_0_yyzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(78));

    t_0_yyzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(79));

    t_0_yyzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(80));

    t_0_yyzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(81));

    t_0_yyzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(82));

    t_0_yyzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(83));

    t_0_yyzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(84));

    t_0_yyzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(85));

    t_0_yyzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(86));

    t_0_yyzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(87));

    t_0_yyzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(88));

    t_0_yyzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(89));

    t_0_yyzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(90));

    t_0_yyzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(91));

    t_0_yyzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(92));

    t_0_yyzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(93));

    t_0_yyzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(94));

    t_0_yyzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(95));

    t_0_yyzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(96));

    t_0_yyzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(97));

    t_0_yyzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(98));

    t_0_yyzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(99));

    t_0_yyzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(100));

    t_0_yyzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(101));

    t_0_yyzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(102));

    t_0_yyzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(103));

    t_0_yyzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(104));

    t_0_yyzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(105));

    t_0_yyzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(106));

    t_0_yyzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(107));

    t_0_yyyz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(108));

    t_0_yyyz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(109));

    t_0_yyyz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(110));

    t_0_yyyz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(111));

    t_0_yyyz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(112));

    t_0_yyyz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(113));

    t_0_yyyz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(114));

    t_0_yyyz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(115));

    t_0_yyyz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(116));

    t_0_yyyz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(117));

    t_0_yyyz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(118));

    t_0_yyyz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(119));

    t_0_yyyz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(120));

    t_0_yyyz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(121));

    t_0_yyyz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(122));

    t_0_yyyz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(123));

    t_0_yyyz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(124));

    t_0_yyyz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(125));

    t_0_yyyz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(126));

    t_0_yyyz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(127));

    t_0_yyyz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(128));

    t_0_yyyz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(129));

    t_0_yyyz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(130));

    t_0_yyyz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(131));

    t_0_yyyz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(132));

    t_0_yyyz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(133));

    t_0_yyyz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(134));

    t_0_yyyz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(135));

    t_0_yyyz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(136));

    t_0_yyyz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(137));

    t_0_yyyz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(138));

    t_0_yyyz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(139));

    t_0_yyyz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(140));

    t_0_yyyz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(141));

    t_0_yyyz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(142));

    t_0_yyyz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(143));

    t_0_yyyy_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(144));

    t_0_yyyy_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(145));

    t_0_yyyy_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(146));

    t_0_yyyy_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(147));

    t_0_yyyy_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(148));

    t_0_yyyy_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(149));

    t_0_yyyy_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(150));

    t_0_yyyy_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(151));

    t_0_yyyy_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(152));

    t_0_yyyy_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(153));

    t_0_yyyy_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(154));

    t_0_yyyy_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(155));

    t_0_yyyy_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(156));

    t_0_yyyy_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(157));

    t_0_yyyy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(158));

    t_0_yyyy_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(159));

    t_0_yyyy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(160));

    t_0_yyyy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(161));

    t_0_yyyy_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(162));

    t_0_yyyy_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(163));

    t_0_yyyy_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(164));

    t_0_yyyy_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(165));

    t_0_yyyy_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(166));

    t_0_yyyy_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(167));

    t_0_yyyy_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(168));

    t_0_yyyy_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(169));

    t_0_yyyy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(170));

    t_0_yyyy_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(171));

    t_0_yyyy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(172));

    t_0_yyyy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(173));

    t_0_yyyy_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(174));

    t_0_yyyy_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(175));

    t_0_yyyy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(176));

    t_0_yyyy_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(177));

    t_0_yyyy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(178));

    t_0_yyyy_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(179));

    t_0_xzzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(180));

    t_0_xzzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(181));

    t_0_xzzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(182));

    t_0_xzzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(183));

    t_0_xzzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(184));

    t_0_xzzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(185));

    t_0_xzzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(186));

    t_0_xzzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(187));

    t_0_xzzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(188));

    t_0_xzzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(189));

    t_0_xzzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(190));

    t_0_xzzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(191));

    t_0_xzzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(192));

    t_0_xzzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(193));

    t_0_xzzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(194));

    t_0_xzzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(195));

    t_0_xzzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(196));

    t_0_xzzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(197));

    t_0_xzzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(198));

    t_0_xzzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(199));

    t_0_xzzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(200));

    t_0_xzzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(201));

    t_0_xzzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(202));

    t_0_xzzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(203));

    t_0_xzzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(204));

    t_0_xzzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(205));

    t_0_xzzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(206));

    t_0_xzzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(207));

    t_0_xzzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(208));

    t_0_xzzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(209));

    t_0_xzzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(210));

    t_0_xzzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(211));

    t_0_xzzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(212));

    t_0_xzzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(213));

    t_0_xzzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(214));

    t_0_xzzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(215));

    t_0_xyzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(216));

    t_0_xyzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(217));

    t_0_xyzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(218));

    t_0_xyzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(219));

    t_0_xyzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(220));

    t_0_xyzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(221));

    t_0_xyzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(222));

    t_0_xyzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(223));

    t_0_xyzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(224));

    t_0_xyzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(225));

    t_0_xyzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(226));

    t_0_xyzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(227));

    t_0_xyzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(228));

    t_0_xyzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(229));

    t_0_xyzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(230));

    t_0_xyzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(231));

    t_0_xyzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(232));

    t_0_xyzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(233));

    t_0_xyzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(234));

    t_0_xyzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(235));

    t_0_xyzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(236));

    t_0_xyzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(237));

    t_0_xyzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(238));

    t_0_xyzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(239));

    t_0_xyzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(240));

    t_0_xyzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(241));

    t_0_xyzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(242));

    t_0_xyzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(243));

    t_0_xyzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(244));

    t_0_xyzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(245));

    t_0_xyzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(246));

    t_0_xyzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(247));

    t_0_xyzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(248));

    t_0_xyzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(249));

    t_0_xyzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(250));

    t_0_xyzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(251));

    t_0_xyyz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(252));

    t_0_xyyz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(253));

    t_0_xyyz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(254));

    t_0_xyyz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(255));

    t_0_xyyz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(256));

    t_0_xyyz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(257));

    t_0_xyyz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(258));

    t_0_xyyz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(259));

    t_0_xyyz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(260));

    t_0_xyyz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(261));

    t_0_xyyz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(262));

    t_0_xyyz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(263));

    t_0_xyyz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(264));

    t_0_xyyz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(265));

    t_0_xyyz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(266));

    t_0_xyyz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(267));

    t_0_xyyz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(268));

    t_0_xyyz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(269));

    t_0_xyyz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(270));

    t_0_xyyz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(271));

    t_0_xyyz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(272));

    t_0_xyyz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(273));

    t_0_xyyz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(274));

    t_0_xyyz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(275));

    t_0_xyyz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(276));

    t_0_xyyz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(277));

    t_0_xyyz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(278));

    t_0_xyyz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(279));

    t_0_xyyz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(280));

    t_0_xyyz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(281));

    t_0_xyyz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(282));

    t_0_xyyz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(283));

    t_0_xyyz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(284));

    t_0_xyyz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(285));

    t_0_xyyz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(286));

    t_0_xyyz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(287));

    t_0_xyyy_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(288));

    t_0_xyyy_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(289));

    t_0_xyyy_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(290));

    t_0_xyyy_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(291));

    t_0_xyyy_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(292));

    t_0_xyyy_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(293));

    t_0_xyyy_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(294));

    t_0_xyyy_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(295));

    t_0_xyyy_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(296));

    t_0_xyyy_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(297));

    t_0_xyyy_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(298));

    t_0_xyyy_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(299));

    t_0_xyyy_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(300));

    t_0_xyyy_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(301));

    t_0_xyyy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(302));

    t_0_xyyy_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(303));

    t_0_xyyy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(304));

    t_0_xyyy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(305));

    t_0_xyyy_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(306));

    t_0_xyyy_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(307));

    t_0_xyyy_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(308));

    t_0_xyyy_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(309));

    t_0_xyyy_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(310));

    t_0_xyyy_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(311));

    t_0_xyyy_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(312));

    t_0_xyyy_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(313));

    t_0_xyyy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(314));

    t_0_xyyy_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(315));

    t_0_xyyy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(316));

    t_0_xyyy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(317));

    t_0_xyyy_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(318));

    t_0_xyyy_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(319));

    t_0_xyyy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(320));

    t_0_xyyy_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(321));

    t_0_xyyy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(322));

    t_0_xyyy_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(323));

    t_0_xxzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(324));

    t_0_xxzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(325));

    t_0_xxzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(326));

    t_0_xxzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(327));

    t_0_xxzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(328));

    t_0_xxzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(329));

    t_0_xxzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(330));

    t_0_xxzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(331));

    t_0_xxzz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(332));

    t_0_xxzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(333));

    t_0_xxzz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(334));

    t_0_xxzz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(335));

    t_0_xxzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(336));

    t_0_xxzz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(337));

    t_0_xxzz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(338));

    t_0_xxzz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(339));

    t_0_xxzz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(340));

    t_0_xxzz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(341));

    t_0_xxzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(342));

    t_0_xxzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(343));

    t_0_xxzz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(344));

    t_0_xxzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(345));

    t_0_xxzz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(346));

    t_0_xxzz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(347));

    t_0_xxzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(348));

    t_0_xxzz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(349));

    t_0_xxzz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(350));

    t_0_xxzz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(351));

    t_0_xxzz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(352));

    t_0_xxzz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(353));

    t_0_xxzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(354));

    t_0_xxzz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(355));

    t_0_xxzz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(356));

    t_0_xxzz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(357));

    t_0_xxzz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(358));

    t_0_xxzz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(359));

    t_0_xxyz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(360));

    t_0_xxyz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(361));

    t_0_xxyz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(362));

    t_0_xxyz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(363));

    t_0_xxyz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(364));

    t_0_xxyz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(365));

    t_0_xxyz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(366));

    t_0_xxyz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(367));

    t_0_xxyz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(368));

    t_0_xxyz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(369));

    t_0_xxyz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(370));

    t_0_xxyz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(371));

    t_0_xxyz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(372));

    t_0_xxyz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(373));

    t_0_xxyz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(374));

    t_0_xxyz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(375));

    t_0_xxyz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(376));

    t_0_xxyz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(377));

    t_0_xxyz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(378));

    t_0_xxyz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(379));

    t_0_xxyz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(380));

    t_0_xxyz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(381));

    t_0_xxyz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(382));

    t_0_xxyz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(383));

    t_0_xxyz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(384));

    t_0_xxyz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(385));

    t_0_xxyz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(386));

    t_0_xxyz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(387));

    t_0_xxyz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(388));

    t_0_xxyz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(389));

    t_0_xxyz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(390));

    t_0_xxyz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(391));

    t_0_xxyz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(392));

    t_0_xxyz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(393));

    t_0_xxyz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(394));

    t_0_xxyz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(395));

    t_0_xxyy_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(396));

    t_0_xxyy_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(397));

    t_0_xxyy_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(398));

    t_0_xxyy_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(399));

    t_0_xxyy_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(400));

    t_0_xxyy_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(401));

    t_0_xxyy_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(402));

    t_0_xxyy_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(403));

    t_0_xxyy_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(404));

    t_0_xxyy_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(405));

    t_0_xxyy_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(406));

    t_0_xxyy_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(407));

    t_0_xxyy_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(408));

    t_0_xxyy_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(409));

    t_0_xxyy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(410));

    t_0_xxyy_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(411));

    t_0_xxyy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(412));

    t_0_xxyy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(413));

    t_0_xxyy_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(414));

    t_0_xxyy_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(415));

    t_0_xxyy_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(416));

    t_0_xxyy_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(417));

    t_0_xxyy_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(418));

    t_0_xxyy_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(419));

    t_0_xxyy_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(420));

    t_0_xxyy_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(421));

    t_0_xxyy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(422));

    t_0_xxyy_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(423));

    t_0_xxyy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(424));

    t_0_xxyy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(425));

    t_0_xxyy_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(426));

    t_0_xxyy_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(427));

    t_0_xxyy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(428));

    t_0_xxyy_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(429));

    t_0_xxyy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(430));

    t_0_xxyy_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(431));

    t_0_xxxz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(432));

    t_0_xxxz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(433));

    t_0_xxxz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(434));

    t_0_xxxz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(435));

    t_0_xxxz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(436));

    t_0_xxxz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(437));

    t_0_xxxz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(438));

    t_0_xxxz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(439));

    t_0_xxxz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(440));

    t_0_xxxz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(441));

    t_0_xxxz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(442));

    t_0_xxxz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(443));

    t_0_xxxz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(444));

    t_0_xxxz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(445));

    t_0_xxxz_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(446));

    t_0_xxxz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(447));

    t_0_xxxz_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(448));

    t_0_xxxz_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(449));

    t_0_xxxz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(450));

    t_0_xxxz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(451));

    t_0_xxxz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(452));

    t_0_xxxz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(453));

    t_0_xxxz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(454));

    t_0_xxxz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(455));

    t_0_xxxz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(456));

    t_0_xxxz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(457));

    t_0_xxxz_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(458));

    t_0_xxxz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(459));

    t_0_xxxz_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(460));

    t_0_xxxz_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(461));

    t_0_xxxz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(462));

    t_0_xxxz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(463));

    t_0_xxxz_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(464));

    t_0_xxxz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(465));

    t_0_xxxz_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(466));

    t_0_xxxz_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(467));

    t_0_xxxy_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(468));

    t_0_xxxy_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(469));

    t_0_xxxy_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(470));

    t_0_xxxy_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(471));

    t_0_xxxy_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(472));

    t_0_xxxy_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(473));

    t_0_xxxy_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(474));

    t_0_xxxy_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(475));

    t_0_xxxy_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(476));

    t_0_xxxy_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(477));

    t_0_xxxy_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(478));

    t_0_xxxy_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(479));

    t_0_xxxy_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(480));

    t_0_xxxy_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(481));

    t_0_xxxy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(482));

    t_0_xxxy_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(483));

    t_0_xxxy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(484));

    t_0_xxxy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(485));

    t_0_xxxy_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(486));

    t_0_xxxy_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(487));

    t_0_xxxy_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(488));

    t_0_xxxy_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(489));

    t_0_xxxy_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(490));

    t_0_xxxy_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(491));

    t_0_xxxy_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(492));

    t_0_xxxy_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(493));

    t_0_xxxy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(494));

    t_0_xxxy_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(495));

    t_0_xxxy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(496));

    t_0_xxxy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(497));

    t_0_xxxy_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(498));

    t_0_xxxy_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(499));

    t_0_xxxy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(500));

    t_0_xxxy_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(501));

    t_0_xxxy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(502));

    t_0_xxxy_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(503));

    t_0_xxxx_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(504));

    t_0_xxxx_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(505));

    t_0_xxxx_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(506));

    t_0_xxxx_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(507));

    t_0_xxxx_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(508));

    t_0_xxxx_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(509));

    t_0_xxxx_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(510));

    t_0_xxxx_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(511));

    t_0_xxxx_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(512));

    t_0_xxxx_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(513));

    t_0_xxxx_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(514));

    t_0_xxxx_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(515));

    t_0_xxxx_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(516));

    t_0_xxxx_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(517));

    t_0_xxxx_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(518));

    t_0_xxxx_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(519));

    t_0_xxxx_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(520));

    t_0_xxxx_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(521));

    t_0_xxxx_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(522));

    t_0_xxxx_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(523));

    t_0_xxxx_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(524));

    t_0_xxxx_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(525));

    t_0_xxxx_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(526));

    t_0_xxxx_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(527));

    t_0_xxxx_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(528));

    t_0_xxxx_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(529));

    t_0_xxxx_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(530));

    t_0_xxxx_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(531));

    t_0_xxxx_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(532));

    t_0_xxxx_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(533));

    t_0_xxxx_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(534));

    t_0_xxxx_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(535));

    t_0_xxxx_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(536));

    t_0_xxxx_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(537));

    t_0_xxxx_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(538));

    t_0_xxxx_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(539));

    #pragma omp simd align(rab_z, t_0_zzz_xx_xx, t_0_zzz_xx_xy, t_0_zzz_xx_xz, t_0_zzz_xx_yy,\
                           t_0_zzz_xx_yz, t_0_zzz_xx_zz, t_0_zzz_xy_xx, t_0_zzz_xy_xy,\
                           t_0_zzz_xy_xz, t_0_zzz_xy_yy, t_0_zzz_xy_yz, t_0_zzz_xy_zz,\
                           t_0_zzz_xz_xx, t_0_zzz_xz_xy, t_0_zzz_xz_xz, t_0_zzz_xz_yy,\
                           t_0_zzz_xz_yz, t_0_zzz_xz_zz, t_0_zzz_yy_xx, t_0_zzz_yy_xy,\
                           t_0_zzz_yy_xz, t_0_zzz_yy_yy, t_0_zzz_yy_yz, t_0_zzz_yy_zz,\
                           t_0_zzz_yz_xx, t_0_zzz_yz_xy, t_0_zzz_yz_xz, t_0_zzz_yz_yy,\
                           t_0_zzz_yz_yz, t_0_zzz_yz_zz, t_0_zzz_zz_xx, t_0_zzz_zz_xy,\
                           t_0_zzz_zz_xz, t_0_zzz_zz_yy, t_0_zzz_zz_yz, t_0_zzz_zz_zz,\
                           t_0_zzzz_xx_xx, t_0_zzzz_xx_xy, t_0_zzzz_xx_xz, t_0_zzzz_xx_yy,\
                           t_0_zzzz_xx_yz, t_0_zzzz_xx_zz, t_0_zzzz_xy_xx, t_0_zzzz_xy_xy,\
                           t_0_zzzz_xy_xz, t_0_zzzz_xy_yy, t_0_zzzz_xy_yz, t_0_zzzz_xy_zz,\
                           t_0_zzzz_xz_xx, t_0_zzzz_xz_xy, t_0_zzzz_xz_xz, t_0_zzzz_xz_yy,\
                           t_0_zzzz_xz_yz, t_0_zzzz_xz_zz, t_0_zzzz_yy_xx, t_0_zzzz_yy_xy,\
                           t_0_zzzz_yy_xz, t_0_zzzz_yy_yy, t_0_zzzz_yy_yz, t_0_zzzz_yy_zz,\
                           t_0_zzzz_yz_xx, t_0_zzzz_yz_xy, t_0_zzzz_yz_xz, t_0_zzzz_yz_yy,\
                           t_0_zzzz_yz_yz, t_0_zzzz_yz_zz, t_0_zzzz_zz_xx, t_0_zzzz_zz_xy,\
                           t_0_zzzz_zz_xz, t_0_zzzz_zz_yy, t_0_zzzz_zz_yz, t_0_zzzz_zz_zz,\
                           t_z_zzz_xx_xx, t_z_zzz_xx_xy, t_z_zzz_xx_xz, t_z_zzz_xx_yy,\
                           t_z_zzz_xx_yz, t_z_zzz_xx_zz, t_z_zzz_xy_xx, t_z_zzz_xy_xy,\
                           t_z_zzz_xy_xz, t_z_zzz_xy_yy, t_z_zzz_xy_yz, t_z_zzz_xy_zz,\
                           t_z_zzz_xz_xx, t_z_zzz_xz_xy, t_z_zzz_xz_xz, t_z_zzz_xz_yy,\
                           t_z_zzz_xz_yz, t_z_zzz_xz_zz, t_z_zzz_yy_xx, t_z_zzz_yy_xy,\
                           t_z_zzz_yy_xz, t_z_zzz_yy_yy, t_z_zzz_yy_yz, t_z_zzz_yy_zz,\
                           t_z_zzz_yz_xx, t_z_zzz_yz_xy, t_z_zzz_yz_xz, t_z_zzz_yz_yy,\
                           t_z_zzz_yz_yz, t_z_zzz_yz_zz, t_z_zzz_zz_xx, t_z_zzz_zz_xy,\
                           t_z_zzz_zz_xz, t_z_zzz_zz_yy, t_z_zzz_zz_yz, t_z_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_zz_zz[i] = t_0_zzzz_zz_zz[i] - rab_z[i] * t_0_zzz_zz_zz[i];

        t_z_zzz_zz_yz[i] = t_0_zzzz_zz_yz[i] - rab_z[i] * t_0_zzz_zz_yz[i];

        t_z_zzz_zz_yy[i] = t_0_zzzz_zz_yy[i] - rab_z[i] * t_0_zzz_zz_yy[i];

        t_z_zzz_zz_xz[i] = t_0_zzzz_zz_xz[i] - rab_z[i] * t_0_zzz_zz_xz[i];

        t_z_zzz_zz_xy[i] = t_0_zzzz_zz_xy[i] - rab_z[i] * t_0_zzz_zz_xy[i];

        t_z_zzz_zz_xx[i] = t_0_zzzz_zz_xx[i] - rab_z[i] * t_0_zzz_zz_xx[i];

        t_z_zzz_yz_zz[i] = t_0_zzzz_yz_zz[i] - rab_z[i] * t_0_zzz_yz_zz[i];

        t_z_zzz_yz_yz[i] = t_0_zzzz_yz_yz[i] - rab_z[i] * t_0_zzz_yz_yz[i];

        t_z_zzz_yz_yy[i] = t_0_zzzz_yz_yy[i] - rab_z[i] * t_0_zzz_yz_yy[i];

        t_z_zzz_yz_xz[i] = t_0_zzzz_yz_xz[i] - rab_z[i] * t_0_zzz_yz_xz[i];

        t_z_zzz_yz_xy[i] = t_0_zzzz_yz_xy[i] - rab_z[i] * t_0_zzz_yz_xy[i];

        t_z_zzz_yz_xx[i] = t_0_zzzz_yz_xx[i] - rab_z[i] * t_0_zzz_yz_xx[i];

        t_z_zzz_yy_zz[i] = t_0_zzzz_yy_zz[i] - rab_z[i] * t_0_zzz_yy_zz[i];

        t_z_zzz_yy_yz[i] = t_0_zzzz_yy_yz[i] - rab_z[i] * t_0_zzz_yy_yz[i];

        t_z_zzz_yy_yy[i] = t_0_zzzz_yy_yy[i] - rab_z[i] * t_0_zzz_yy_yy[i];

        t_z_zzz_yy_xz[i] = t_0_zzzz_yy_xz[i] - rab_z[i] * t_0_zzz_yy_xz[i];

        t_z_zzz_yy_xy[i] = t_0_zzzz_yy_xy[i] - rab_z[i] * t_0_zzz_yy_xy[i];

        t_z_zzz_yy_xx[i] = t_0_zzzz_yy_xx[i] - rab_z[i] * t_0_zzz_yy_xx[i];

        t_z_zzz_xz_zz[i] = t_0_zzzz_xz_zz[i] - rab_z[i] * t_0_zzz_xz_zz[i];

        t_z_zzz_xz_yz[i] = t_0_zzzz_xz_yz[i] - rab_z[i] * t_0_zzz_xz_yz[i];

        t_z_zzz_xz_yy[i] = t_0_zzzz_xz_yy[i] - rab_z[i] * t_0_zzz_xz_yy[i];

        t_z_zzz_xz_xz[i] = t_0_zzzz_xz_xz[i] - rab_z[i] * t_0_zzz_xz_xz[i];

        t_z_zzz_xz_xy[i] = t_0_zzzz_xz_xy[i] - rab_z[i] * t_0_zzz_xz_xy[i];

        t_z_zzz_xz_xx[i] = t_0_zzzz_xz_xx[i] - rab_z[i] * t_0_zzz_xz_xx[i];

        t_z_zzz_xy_zz[i] = t_0_zzzz_xy_zz[i] - rab_z[i] * t_0_zzz_xy_zz[i];

        t_z_zzz_xy_yz[i] = t_0_zzzz_xy_yz[i] - rab_z[i] * t_0_zzz_xy_yz[i];

        t_z_zzz_xy_yy[i] = t_0_zzzz_xy_yy[i] - rab_z[i] * t_0_zzz_xy_yy[i];

        t_z_zzz_xy_xz[i] = t_0_zzzz_xy_xz[i] - rab_z[i] * t_0_zzz_xy_xz[i];

        t_z_zzz_xy_xy[i] = t_0_zzzz_xy_xy[i] - rab_z[i] * t_0_zzz_xy_xy[i];

        t_z_zzz_xy_xx[i] = t_0_zzzz_xy_xx[i] - rab_z[i] * t_0_zzz_xy_xx[i];

        t_z_zzz_xx_zz[i] = t_0_zzzz_xx_zz[i] - rab_z[i] * t_0_zzz_xx_zz[i];

        t_z_zzz_xx_yz[i] = t_0_zzzz_xx_yz[i] - rab_z[i] * t_0_zzz_xx_yz[i];

        t_z_zzz_xx_yy[i] = t_0_zzzz_xx_yy[i] - rab_z[i] * t_0_zzz_xx_yy[i];

        t_z_zzz_xx_xz[i] = t_0_zzzz_xx_xz[i] - rab_z[i] * t_0_zzz_xx_xz[i];

        t_z_zzz_xx_xy[i] = t_0_zzzz_xx_xy[i] - rab_z[i] * t_0_zzz_xx_xy[i];

        t_z_zzz_xx_xx[i] = t_0_zzzz_xx_xx[i] - rab_z[i] * t_0_zzz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_yzz_xx_xx, t_0_yzz_xx_xy, t_0_yzz_xx_xz, t_0_yzz_xx_yy,\
                           t_0_yzz_xx_yz, t_0_yzz_xx_zz, t_0_yzz_xy_xx, t_0_yzz_xy_xy,\
                           t_0_yzz_xy_xz, t_0_yzz_xy_yy, t_0_yzz_xy_yz, t_0_yzz_xy_zz,\
                           t_0_yzz_xz_xx, t_0_yzz_xz_xy, t_0_yzz_xz_xz, t_0_yzz_xz_yy,\
                           t_0_yzz_xz_yz, t_0_yzz_xz_zz, t_0_yzz_yy_xx, t_0_yzz_yy_xy,\
                           t_0_yzz_yy_xz, t_0_yzz_yy_yy, t_0_yzz_yy_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_xx, t_0_yzz_yz_xy, t_0_yzz_yz_xz, t_0_yzz_yz_yy,\
                           t_0_yzz_yz_yz, t_0_yzz_yz_zz, t_0_yzz_zz_xx, t_0_yzz_zz_xy,\
                           t_0_yzz_zz_xz, t_0_yzz_zz_yy, t_0_yzz_zz_yz, t_0_yzz_zz_zz,\
                           t_0_yzzz_xx_xx, t_0_yzzz_xx_xy, t_0_yzzz_xx_xz, t_0_yzzz_xx_yy,\
                           t_0_yzzz_xx_yz, t_0_yzzz_xx_zz, t_0_yzzz_xy_xx, t_0_yzzz_xy_xy,\
                           t_0_yzzz_xy_xz, t_0_yzzz_xy_yy, t_0_yzzz_xy_yz, t_0_yzzz_xy_zz,\
                           t_0_yzzz_xz_xx, t_0_yzzz_xz_xy, t_0_yzzz_xz_xz, t_0_yzzz_xz_yy,\
                           t_0_yzzz_xz_yz, t_0_yzzz_xz_zz, t_0_yzzz_yy_xx, t_0_yzzz_yy_xy,\
                           t_0_yzzz_yy_xz, t_0_yzzz_yy_yy, t_0_yzzz_yy_yz, t_0_yzzz_yy_zz,\
                           t_0_yzzz_yz_xx, t_0_yzzz_yz_xy, t_0_yzzz_yz_xz, t_0_yzzz_yz_yy,\
                           t_0_yzzz_yz_yz, t_0_yzzz_yz_zz, t_0_yzzz_zz_xx, t_0_yzzz_zz_xy,\
                           t_0_yzzz_zz_xz, t_0_yzzz_zz_yy, t_0_yzzz_zz_yz, t_0_yzzz_zz_zz,\
                           t_z_yzz_xx_xx, t_z_yzz_xx_xy, t_z_yzz_xx_xz, t_z_yzz_xx_yy,\
                           t_z_yzz_xx_yz, t_z_yzz_xx_zz, t_z_yzz_xy_xx, t_z_yzz_xy_xy,\
                           t_z_yzz_xy_xz, t_z_yzz_xy_yy, t_z_yzz_xy_yz, t_z_yzz_xy_zz,\
                           t_z_yzz_xz_xx, t_z_yzz_xz_xy, t_z_yzz_xz_xz, t_z_yzz_xz_yy,\
                           t_z_yzz_xz_yz, t_z_yzz_xz_zz, t_z_yzz_yy_xx, t_z_yzz_yy_xy,\
                           t_z_yzz_yy_xz, t_z_yzz_yy_yy, t_z_yzz_yy_yz, t_z_yzz_yy_zz,\
                           t_z_yzz_yz_xx, t_z_yzz_yz_xy, t_z_yzz_yz_xz, t_z_yzz_yz_yy,\
                           t_z_yzz_yz_yz, t_z_yzz_yz_zz, t_z_yzz_zz_xx, t_z_yzz_zz_xy,\
                           t_z_yzz_zz_xz, t_z_yzz_zz_yy, t_z_yzz_zz_yz, t_z_yzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yzz_zz_zz[i] = t_0_yzzz_zz_zz[i] - rab_z[i] * t_0_yzz_zz_zz[i];

        t_z_yzz_zz_yz[i] = t_0_yzzz_zz_yz[i] - rab_z[i] * t_0_yzz_zz_yz[i];

        t_z_yzz_zz_yy[i] = t_0_yzzz_zz_yy[i] - rab_z[i] * t_0_yzz_zz_yy[i];

        t_z_yzz_zz_xz[i] = t_0_yzzz_zz_xz[i] - rab_z[i] * t_0_yzz_zz_xz[i];

        t_z_yzz_zz_xy[i] = t_0_yzzz_zz_xy[i] - rab_z[i] * t_0_yzz_zz_xy[i];

        t_z_yzz_zz_xx[i] = t_0_yzzz_zz_xx[i] - rab_z[i] * t_0_yzz_zz_xx[i];

        t_z_yzz_yz_zz[i] = t_0_yzzz_yz_zz[i] - rab_z[i] * t_0_yzz_yz_zz[i];

        t_z_yzz_yz_yz[i] = t_0_yzzz_yz_yz[i] - rab_z[i] * t_0_yzz_yz_yz[i];

        t_z_yzz_yz_yy[i] = t_0_yzzz_yz_yy[i] - rab_z[i] * t_0_yzz_yz_yy[i];

        t_z_yzz_yz_xz[i] = t_0_yzzz_yz_xz[i] - rab_z[i] * t_0_yzz_yz_xz[i];

        t_z_yzz_yz_xy[i] = t_0_yzzz_yz_xy[i] - rab_z[i] * t_0_yzz_yz_xy[i];

        t_z_yzz_yz_xx[i] = t_0_yzzz_yz_xx[i] - rab_z[i] * t_0_yzz_yz_xx[i];

        t_z_yzz_yy_zz[i] = t_0_yzzz_yy_zz[i] - rab_z[i] * t_0_yzz_yy_zz[i];

        t_z_yzz_yy_yz[i] = t_0_yzzz_yy_yz[i] - rab_z[i] * t_0_yzz_yy_yz[i];

        t_z_yzz_yy_yy[i] = t_0_yzzz_yy_yy[i] - rab_z[i] * t_0_yzz_yy_yy[i];

        t_z_yzz_yy_xz[i] = t_0_yzzz_yy_xz[i] - rab_z[i] * t_0_yzz_yy_xz[i];

        t_z_yzz_yy_xy[i] = t_0_yzzz_yy_xy[i] - rab_z[i] * t_0_yzz_yy_xy[i];

        t_z_yzz_yy_xx[i] = t_0_yzzz_yy_xx[i] - rab_z[i] * t_0_yzz_yy_xx[i];

        t_z_yzz_xz_zz[i] = t_0_yzzz_xz_zz[i] - rab_z[i] * t_0_yzz_xz_zz[i];

        t_z_yzz_xz_yz[i] = t_0_yzzz_xz_yz[i] - rab_z[i] * t_0_yzz_xz_yz[i];

        t_z_yzz_xz_yy[i] = t_0_yzzz_xz_yy[i] - rab_z[i] * t_0_yzz_xz_yy[i];

        t_z_yzz_xz_xz[i] = t_0_yzzz_xz_xz[i] - rab_z[i] * t_0_yzz_xz_xz[i];

        t_z_yzz_xz_xy[i] = t_0_yzzz_xz_xy[i] - rab_z[i] * t_0_yzz_xz_xy[i];

        t_z_yzz_xz_xx[i] = t_0_yzzz_xz_xx[i] - rab_z[i] * t_0_yzz_xz_xx[i];

        t_z_yzz_xy_zz[i] = t_0_yzzz_xy_zz[i] - rab_z[i] * t_0_yzz_xy_zz[i];

        t_z_yzz_xy_yz[i] = t_0_yzzz_xy_yz[i] - rab_z[i] * t_0_yzz_xy_yz[i];

        t_z_yzz_xy_yy[i] = t_0_yzzz_xy_yy[i] - rab_z[i] * t_0_yzz_xy_yy[i];

        t_z_yzz_xy_xz[i] = t_0_yzzz_xy_xz[i] - rab_z[i] * t_0_yzz_xy_xz[i];

        t_z_yzz_xy_xy[i] = t_0_yzzz_xy_xy[i] - rab_z[i] * t_0_yzz_xy_xy[i];

        t_z_yzz_xy_xx[i] = t_0_yzzz_xy_xx[i] - rab_z[i] * t_0_yzz_xy_xx[i];

        t_z_yzz_xx_zz[i] = t_0_yzzz_xx_zz[i] - rab_z[i] * t_0_yzz_xx_zz[i];

        t_z_yzz_xx_yz[i] = t_0_yzzz_xx_yz[i] - rab_z[i] * t_0_yzz_xx_yz[i];

        t_z_yzz_xx_yy[i] = t_0_yzzz_xx_yy[i] - rab_z[i] * t_0_yzz_xx_yy[i];

        t_z_yzz_xx_xz[i] = t_0_yzzz_xx_xz[i] - rab_z[i] * t_0_yzz_xx_xz[i];

        t_z_yzz_xx_xy[i] = t_0_yzzz_xx_xy[i] - rab_z[i] * t_0_yzz_xx_xy[i];

        t_z_yzz_xx_xx[i] = t_0_yzzz_xx_xx[i] - rab_z[i] * t_0_yzz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_yyz_xx_xx, t_0_yyz_xx_xy, t_0_yyz_xx_xz, t_0_yyz_xx_yy,\
                           t_0_yyz_xx_yz, t_0_yyz_xx_zz, t_0_yyz_xy_xx, t_0_yyz_xy_xy,\
                           t_0_yyz_xy_xz, t_0_yyz_xy_yy, t_0_yyz_xy_yz, t_0_yyz_xy_zz,\
                           t_0_yyz_xz_xx, t_0_yyz_xz_xy, t_0_yyz_xz_xz, t_0_yyz_xz_yy,\
                           t_0_yyz_xz_yz, t_0_yyz_xz_zz, t_0_yyz_yy_xx, t_0_yyz_yy_xy,\
                           t_0_yyz_yy_xz, t_0_yyz_yy_yy, t_0_yyz_yy_yz, t_0_yyz_yy_zz,\
                           t_0_yyz_yz_xx, t_0_yyz_yz_xy, t_0_yyz_yz_xz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_yz_zz, t_0_yyz_zz_xx, t_0_yyz_zz_xy,\
                           t_0_yyz_zz_xz, t_0_yyz_zz_yy, t_0_yyz_zz_yz, t_0_yyz_zz_zz,\
                           t_0_yyzz_xx_xx, t_0_yyzz_xx_xy, t_0_yyzz_xx_xz, t_0_yyzz_xx_yy,\
                           t_0_yyzz_xx_yz, t_0_yyzz_xx_zz, t_0_yyzz_xy_xx, t_0_yyzz_xy_xy,\
                           t_0_yyzz_xy_xz, t_0_yyzz_xy_yy, t_0_yyzz_xy_yz, t_0_yyzz_xy_zz,\
                           t_0_yyzz_xz_xx, t_0_yyzz_xz_xy, t_0_yyzz_xz_xz, t_0_yyzz_xz_yy,\
                           t_0_yyzz_xz_yz, t_0_yyzz_xz_zz, t_0_yyzz_yy_xx, t_0_yyzz_yy_xy,\
                           t_0_yyzz_yy_xz, t_0_yyzz_yy_yy, t_0_yyzz_yy_yz, t_0_yyzz_yy_zz,\
                           t_0_yyzz_yz_xx, t_0_yyzz_yz_xy, t_0_yyzz_yz_xz, t_0_yyzz_yz_yy,\
                           t_0_yyzz_yz_yz, t_0_yyzz_yz_zz, t_0_yyzz_zz_xx, t_0_yyzz_zz_xy,\
                           t_0_yyzz_zz_xz, t_0_yyzz_zz_yy, t_0_yyzz_zz_yz, t_0_yyzz_zz_zz,\
                           t_z_yyz_xx_xx, t_z_yyz_xx_xy, t_z_yyz_xx_xz, t_z_yyz_xx_yy,\
                           t_z_yyz_xx_yz, t_z_yyz_xx_zz, t_z_yyz_xy_xx, t_z_yyz_xy_xy,\
                           t_z_yyz_xy_xz, t_z_yyz_xy_yy, t_z_yyz_xy_yz, t_z_yyz_xy_zz,\
                           t_z_yyz_xz_xx, t_z_yyz_xz_xy, t_z_yyz_xz_xz, t_z_yyz_xz_yy,\
                           t_z_yyz_xz_yz, t_z_yyz_xz_zz, t_z_yyz_yy_xx, t_z_yyz_yy_xy,\
                           t_z_yyz_yy_xz, t_z_yyz_yy_yy, t_z_yyz_yy_yz, t_z_yyz_yy_zz,\
                           t_z_yyz_yz_xx, t_z_yyz_yz_xy, t_z_yyz_yz_xz, t_z_yyz_yz_yy,\
                           t_z_yyz_yz_yz, t_z_yyz_yz_zz, t_z_yyz_zz_xx, t_z_yyz_zz_xy,\
                           t_z_yyz_zz_xz, t_z_yyz_zz_yy, t_z_yyz_zz_yz, t_z_yyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yyz_zz_zz[i] = t_0_yyzz_zz_zz[i] - rab_z[i] * t_0_yyz_zz_zz[i];

        t_z_yyz_zz_yz[i] = t_0_yyzz_zz_yz[i] - rab_z[i] * t_0_yyz_zz_yz[i];

        t_z_yyz_zz_yy[i] = t_0_yyzz_zz_yy[i] - rab_z[i] * t_0_yyz_zz_yy[i];

        t_z_yyz_zz_xz[i] = t_0_yyzz_zz_xz[i] - rab_z[i] * t_0_yyz_zz_xz[i];

        t_z_yyz_zz_xy[i] = t_0_yyzz_zz_xy[i] - rab_z[i] * t_0_yyz_zz_xy[i];

        t_z_yyz_zz_xx[i] = t_0_yyzz_zz_xx[i] - rab_z[i] * t_0_yyz_zz_xx[i];

        t_z_yyz_yz_zz[i] = t_0_yyzz_yz_zz[i] - rab_z[i] * t_0_yyz_yz_zz[i];

        t_z_yyz_yz_yz[i] = t_0_yyzz_yz_yz[i] - rab_z[i] * t_0_yyz_yz_yz[i];

        t_z_yyz_yz_yy[i] = t_0_yyzz_yz_yy[i] - rab_z[i] * t_0_yyz_yz_yy[i];

        t_z_yyz_yz_xz[i] = t_0_yyzz_yz_xz[i] - rab_z[i] * t_0_yyz_yz_xz[i];

        t_z_yyz_yz_xy[i] = t_0_yyzz_yz_xy[i] - rab_z[i] * t_0_yyz_yz_xy[i];

        t_z_yyz_yz_xx[i] = t_0_yyzz_yz_xx[i] - rab_z[i] * t_0_yyz_yz_xx[i];

        t_z_yyz_yy_zz[i] = t_0_yyzz_yy_zz[i] - rab_z[i] * t_0_yyz_yy_zz[i];

        t_z_yyz_yy_yz[i] = t_0_yyzz_yy_yz[i] - rab_z[i] * t_0_yyz_yy_yz[i];

        t_z_yyz_yy_yy[i] = t_0_yyzz_yy_yy[i] - rab_z[i] * t_0_yyz_yy_yy[i];

        t_z_yyz_yy_xz[i] = t_0_yyzz_yy_xz[i] - rab_z[i] * t_0_yyz_yy_xz[i];

        t_z_yyz_yy_xy[i] = t_0_yyzz_yy_xy[i] - rab_z[i] * t_0_yyz_yy_xy[i];

        t_z_yyz_yy_xx[i] = t_0_yyzz_yy_xx[i] - rab_z[i] * t_0_yyz_yy_xx[i];

        t_z_yyz_xz_zz[i] = t_0_yyzz_xz_zz[i] - rab_z[i] * t_0_yyz_xz_zz[i];

        t_z_yyz_xz_yz[i] = t_0_yyzz_xz_yz[i] - rab_z[i] * t_0_yyz_xz_yz[i];

        t_z_yyz_xz_yy[i] = t_0_yyzz_xz_yy[i] - rab_z[i] * t_0_yyz_xz_yy[i];

        t_z_yyz_xz_xz[i] = t_0_yyzz_xz_xz[i] - rab_z[i] * t_0_yyz_xz_xz[i];

        t_z_yyz_xz_xy[i] = t_0_yyzz_xz_xy[i] - rab_z[i] * t_0_yyz_xz_xy[i];

        t_z_yyz_xz_xx[i] = t_0_yyzz_xz_xx[i] - rab_z[i] * t_0_yyz_xz_xx[i];

        t_z_yyz_xy_zz[i] = t_0_yyzz_xy_zz[i] - rab_z[i] * t_0_yyz_xy_zz[i];

        t_z_yyz_xy_yz[i] = t_0_yyzz_xy_yz[i] - rab_z[i] * t_0_yyz_xy_yz[i];

        t_z_yyz_xy_yy[i] = t_0_yyzz_xy_yy[i] - rab_z[i] * t_0_yyz_xy_yy[i];

        t_z_yyz_xy_xz[i] = t_0_yyzz_xy_xz[i] - rab_z[i] * t_0_yyz_xy_xz[i];

        t_z_yyz_xy_xy[i] = t_0_yyzz_xy_xy[i] - rab_z[i] * t_0_yyz_xy_xy[i];

        t_z_yyz_xy_xx[i] = t_0_yyzz_xy_xx[i] - rab_z[i] * t_0_yyz_xy_xx[i];

        t_z_yyz_xx_zz[i] = t_0_yyzz_xx_zz[i] - rab_z[i] * t_0_yyz_xx_zz[i];

        t_z_yyz_xx_yz[i] = t_0_yyzz_xx_yz[i] - rab_z[i] * t_0_yyz_xx_yz[i];

        t_z_yyz_xx_yy[i] = t_0_yyzz_xx_yy[i] - rab_z[i] * t_0_yyz_xx_yy[i];

        t_z_yyz_xx_xz[i] = t_0_yyzz_xx_xz[i] - rab_z[i] * t_0_yyz_xx_xz[i];

        t_z_yyz_xx_xy[i] = t_0_yyzz_xx_xy[i] - rab_z[i] * t_0_yyz_xx_xy[i];

        t_z_yyz_xx_xx[i] = t_0_yyzz_xx_xx[i] - rab_z[i] * t_0_yyz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xzz_xx_xx, t_0_xzz_xx_xy, t_0_xzz_xx_xz, t_0_xzz_xx_yy,\
                           t_0_xzz_xx_yz, t_0_xzz_xx_zz, t_0_xzz_xy_xx, t_0_xzz_xy_xy,\
                           t_0_xzz_xy_xz, t_0_xzz_xy_yy, t_0_xzz_xy_yz, t_0_xzz_xy_zz,\
                           t_0_xzz_xz_xx, t_0_xzz_xz_xy, t_0_xzz_xz_xz, t_0_xzz_xz_yy,\
                           t_0_xzz_xz_yz, t_0_xzz_xz_zz, t_0_xzz_yy_xx, t_0_xzz_yy_xy,\
                           t_0_xzz_yy_xz, t_0_xzz_yy_yy, t_0_xzz_yy_yz, t_0_xzz_yy_zz,\
                           t_0_xzz_yz_xx, t_0_xzz_yz_xy, t_0_xzz_yz_xz, t_0_xzz_yz_yy,\
                           t_0_xzz_yz_yz, t_0_xzz_yz_zz, t_0_xzz_zz_xx, t_0_xzz_zz_xy,\
                           t_0_xzz_zz_xz, t_0_xzz_zz_yy, t_0_xzz_zz_yz, t_0_xzz_zz_zz,\
                           t_0_xzzz_xx_xx, t_0_xzzz_xx_xy, t_0_xzzz_xx_xz, t_0_xzzz_xx_yy,\
                           t_0_xzzz_xx_yz, t_0_xzzz_xx_zz, t_0_xzzz_xy_xx, t_0_xzzz_xy_xy,\
                           t_0_xzzz_xy_xz, t_0_xzzz_xy_yy, t_0_xzzz_xy_yz, t_0_xzzz_xy_zz,\
                           t_0_xzzz_xz_xx, t_0_xzzz_xz_xy, t_0_xzzz_xz_xz, t_0_xzzz_xz_yy,\
                           t_0_xzzz_xz_yz, t_0_xzzz_xz_zz, t_0_xzzz_yy_xx, t_0_xzzz_yy_xy,\
                           t_0_xzzz_yy_xz, t_0_xzzz_yy_yy, t_0_xzzz_yy_yz, t_0_xzzz_yy_zz,\
                           t_0_xzzz_yz_xx, t_0_xzzz_yz_xy, t_0_xzzz_yz_xz, t_0_xzzz_yz_yy,\
                           t_0_xzzz_yz_yz, t_0_xzzz_yz_zz, t_0_xzzz_zz_xx, t_0_xzzz_zz_xy,\
                           t_0_xzzz_zz_xz, t_0_xzzz_zz_yy, t_0_xzzz_zz_yz, t_0_xzzz_zz_zz,\
                           t_z_xzz_xx_xx, t_z_xzz_xx_xy, t_z_xzz_xx_xz, t_z_xzz_xx_yy,\
                           t_z_xzz_xx_yz, t_z_xzz_xx_zz, t_z_xzz_xy_xx, t_z_xzz_xy_xy,\
                           t_z_xzz_xy_xz, t_z_xzz_xy_yy, t_z_xzz_xy_yz, t_z_xzz_xy_zz,\
                           t_z_xzz_xz_xx, t_z_xzz_xz_xy, t_z_xzz_xz_xz, t_z_xzz_xz_yy,\
                           t_z_xzz_xz_yz, t_z_xzz_xz_zz, t_z_xzz_yy_xx, t_z_xzz_yy_xy,\
                           t_z_xzz_yy_xz, t_z_xzz_yy_yy, t_z_xzz_yy_yz, t_z_xzz_yy_zz,\
                           t_z_xzz_yz_xx, t_z_xzz_yz_xy, t_z_xzz_yz_xz, t_z_xzz_yz_yy,\
                           t_z_xzz_yz_yz, t_z_xzz_yz_zz, t_z_xzz_zz_xx, t_z_xzz_zz_xy,\
                           t_z_xzz_zz_xz, t_z_xzz_zz_yy, t_z_xzz_zz_yz, t_z_xzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xzz_zz_zz[i] = t_0_xzzz_zz_zz[i] - rab_z[i] * t_0_xzz_zz_zz[i];

        t_z_xzz_zz_yz[i] = t_0_xzzz_zz_yz[i] - rab_z[i] * t_0_xzz_zz_yz[i];

        t_z_xzz_zz_yy[i] = t_0_xzzz_zz_yy[i] - rab_z[i] * t_0_xzz_zz_yy[i];

        t_z_xzz_zz_xz[i] = t_0_xzzz_zz_xz[i] - rab_z[i] * t_0_xzz_zz_xz[i];

        t_z_xzz_zz_xy[i] = t_0_xzzz_zz_xy[i] - rab_z[i] * t_0_xzz_zz_xy[i];

        t_z_xzz_zz_xx[i] = t_0_xzzz_zz_xx[i] - rab_z[i] * t_0_xzz_zz_xx[i];

        t_z_xzz_yz_zz[i] = t_0_xzzz_yz_zz[i] - rab_z[i] * t_0_xzz_yz_zz[i];

        t_z_xzz_yz_yz[i] = t_0_xzzz_yz_yz[i] - rab_z[i] * t_0_xzz_yz_yz[i];

        t_z_xzz_yz_yy[i] = t_0_xzzz_yz_yy[i] - rab_z[i] * t_0_xzz_yz_yy[i];

        t_z_xzz_yz_xz[i] = t_0_xzzz_yz_xz[i] - rab_z[i] * t_0_xzz_yz_xz[i];

        t_z_xzz_yz_xy[i] = t_0_xzzz_yz_xy[i] - rab_z[i] * t_0_xzz_yz_xy[i];

        t_z_xzz_yz_xx[i] = t_0_xzzz_yz_xx[i] - rab_z[i] * t_0_xzz_yz_xx[i];

        t_z_xzz_yy_zz[i] = t_0_xzzz_yy_zz[i] - rab_z[i] * t_0_xzz_yy_zz[i];

        t_z_xzz_yy_yz[i] = t_0_xzzz_yy_yz[i] - rab_z[i] * t_0_xzz_yy_yz[i];

        t_z_xzz_yy_yy[i] = t_0_xzzz_yy_yy[i] - rab_z[i] * t_0_xzz_yy_yy[i];

        t_z_xzz_yy_xz[i] = t_0_xzzz_yy_xz[i] - rab_z[i] * t_0_xzz_yy_xz[i];

        t_z_xzz_yy_xy[i] = t_0_xzzz_yy_xy[i] - rab_z[i] * t_0_xzz_yy_xy[i];

        t_z_xzz_yy_xx[i] = t_0_xzzz_yy_xx[i] - rab_z[i] * t_0_xzz_yy_xx[i];

        t_z_xzz_xz_zz[i] = t_0_xzzz_xz_zz[i] - rab_z[i] * t_0_xzz_xz_zz[i];

        t_z_xzz_xz_yz[i] = t_0_xzzz_xz_yz[i] - rab_z[i] * t_0_xzz_xz_yz[i];

        t_z_xzz_xz_yy[i] = t_0_xzzz_xz_yy[i] - rab_z[i] * t_0_xzz_xz_yy[i];

        t_z_xzz_xz_xz[i] = t_0_xzzz_xz_xz[i] - rab_z[i] * t_0_xzz_xz_xz[i];

        t_z_xzz_xz_xy[i] = t_0_xzzz_xz_xy[i] - rab_z[i] * t_0_xzz_xz_xy[i];

        t_z_xzz_xz_xx[i] = t_0_xzzz_xz_xx[i] - rab_z[i] * t_0_xzz_xz_xx[i];

        t_z_xzz_xy_zz[i] = t_0_xzzz_xy_zz[i] - rab_z[i] * t_0_xzz_xy_zz[i];

        t_z_xzz_xy_yz[i] = t_0_xzzz_xy_yz[i] - rab_z[i] * t_0_xzz_xy_yz[i];

        t_z_xzz_xy_yy[i] = t_0_xzzz_xy_yy[i] - rab_z[i] * t_0_xzz_xy_yy[i];

        t_z_xzz_xy_xz[i] = t_0_xzzz_xy_xz[i] - rab_z[i] * t_0_xzz_xy_xz[i];

        t_z_xzz_xy_xy[i] = t_0_xzzz_xy_xy[i] - rab_z[i] * t_0_xzz_xy_xy[i];

        t_z_xzz_xy_xx[i] = t_0_xzzz_xy_xx[i] - rab_z[i] * t_0_xzz_xy_xx[i];

        t_z_xzz_xx_zz[i] = t_0_xzzz_xx_zz[i] - rab_z[i] * t_0_xzz_xx_zz[i];

        t_z_xzz_xx_yz[i] = t_0_xzzz_xx_yz[i] - rab_z[i] * t_0_xzz_xx_yz[i];

        t_z_xzz_xx_yy[i] = t_0_xzzz_xx_yy[i] - rab_z[i] * t_0_xzz_xx_yy[i];

        t_z_xzz_xx_xz[i] = t_0_xzzz_xx_xz[i] - rab_z[i] * t_0_xzz_xx_xz[i];

        t_z_xzz_xx_xy[i] = t_0_xzzz_xx_xy[i] - rab_z[i] * t_0_xzz_xx_xy[i];

        t_z_xzz_xx_xx[i] = t_0_xzzz_xx_xx[i] - rab_z[i] * t_0_xzz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz, t_0_xyz_xx_yy,\
                           t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx, t_0_xyz_xy_xy,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz, t_0_xyz_xy_zz,\
                           t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz, t_0_xyz_xz_yy,\
                           t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz,\
                           t_0_xyzz_xx_xx, t_0_xyzz_xx_xy, t_0_xyzz_xx_xz, t_0_xyzz_xx_yy,\
                           t_0_xyzz_xx_yz, t_0_xyzz_xx_zz, t_0_xyzz_xy_xx, t_0_xyzz_xy_xy,\
                           t_0_xyzz_xy_xz, t_0_xyzz_xy_yy, t_0_xyzz_xy_yz, t_0_xyzz_xy_zz,\
                           t_0_xyzz_xz_xx, t_0_xyzz_xz_xy, t_0_xyzz_xz_xz, t_0_xyzz_xz_yy,\
                           t_0_xyzz_xz_yz, t_0_xyzz_xz_zz, t_0_xyzz_yy_xx, t_0_xyzz_yy_xy,\
                           t_0_xyzz_yy_xz, t_0_xyzz_yy_yy, t_0_xyzz_yy_yz, t_0_xyzz_yy_zz,\
                           t_0_xyzz_yz_xx, t_0_xyzz_yz_xy, t_0_xyzz_yz_xz, t_0_xyzz_yz_yy,\
                           t_0_xyzz_yz_yz, t_0_xyzz_yz_zz, t_0_xyzz_zz_xx, t_0_xyzz_zz_xy,\
                           t_0_xyzz_zz_xz, t_0_xyzz_zz_yy, t_0_xyzz_zz_yz, t_0_xyzz_zz_zz,\
                           t_z_xyz_xx_xx, t_z_xyz_xx_xy, t_z_xyz_xx_xz, t_z_xyz_xx_yy,\
                           t_z_xyz_xx_yz, t_z_xyz_xx_zz, t_z_xyz_xy_xx, t_z_xyz_xy_xy,\
                           t_z_xyz_xy_xz, t_z_xyz_xy_yy, t_z_xyz_xy_yz, t_z_xyz_xy_zz,\
                           t_z_xyz_xz_xx, t_z_xyz_xz_xy, t_z_xyz_xz_xz, t_z_xyz_xz_yy,\
                           t_z_xyz_xz_yz, t_z_xyz_xz_zz, t_z_xyz_yy_xx, t_z_xyz_yy_xy,\
                           t_z_xyz_yy_xz, t_z_xyz_yy_yy, t_z_xyz_yy_yz, t_z_xyz_yy_zz,\
                           t_z_xyz_yz_xx, t_z_xyz_yz_xy, t_z_xyz_yz_xz, t_z_xyz_yz_yy,\
                           t_z_xyz_yz_yz, t_z_xyz_yz_zz, t_z_xyz_zz_xx, t_z_xyz_zz_xy,\
                           t_z_xyz_zz_xz, t_z_xyz_zz_yy, t_z_xyz_zz_yz, t_z_xyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xyz_zz_zz[i] = t_0_xyzz_zz_zz[i] - rab_z[i] * t_0_xyz_zz_zz[i];

        t_z_xyz_zz_yz[i] = t_0_xyzz_zz_yz[i] - rab_z[i] * t_0_xyz_zz_yz[i];

        t_z_xyz_zz_yy[i] = t_0_xyzz_zz_yy[i] - rab_z[i] * t_0_xyz_zz_yy[i];

        t_z_xyz_zz_xz[i] = t_0_xyzz_zz_xz[i] - rab_z[i] * t_0_xyz_zz_xz[i];

        t_z_xyz_zz_xy[i] = t_0_xyzz_zz_xy[i] - rab_z[i] * t_0_xyz_zz_xy[i];

        t_z_xyz_zz_xx[i] = t_0_xyzz_zz_xx[i] - rab_z[i] * t_0_xyz_zz_xx[i];

        t_z_xyz_yz_zz[i] = t_0_xyzz_yz_zz[i] - rab_z[i] * t_0_xyz_yz_zz[i];

        t_z_xyz_yz_yz[i] = t_0_xyzz_yz_yz[i] - rab_z[i] * t_0_xyz_yz_yz[i];

        t_z_xyz_yz_yy[i] = t_0_xyzz_yz_yy[i] - rab_z[i] * t_0_xyz_yz_yy[i];

        t_z_xyz_yz_xz[i] = t_0_xyzz_yz_xz[i] - rab_z[i] * t_0_xyz_yz_xz[i];

        t_z_xyz_yz_xy[i] = t_0_xyzz_yz_xy[i] - rab_z[i] * t_0_xyz_yz_xy[i];

        t_z_xyz_yz_xx[i] = t_0_xyzz_yz_xx[i] - rab_z[i] * t_0_xyz_yz_xx[i];

        t_z_xyz_yy_zz[i] = t_0_xyzz_yy_zz[i] - rab_z[i] * t_0_xyz_yy_zz[i];

        t_z_xyz_yy_yz[i] = t_0_xyzz_yy_yz[i] - rab_z[i] * t_0_xyz_yy_yz[i];

        t_z_xyz_yy_yy[i] = t_0_xyzz_yy_yy[i] - rab_z[i] * t_0_xyz_yy_yy[i];

        t_z_xyz_yy_xz[i] = t_0_xyzz_yy_xz[i] - rab_z[i] * t_0_xyz_yy_xz[i];

        t_z_xyz_yy_xy[i] = t_0_xyzz_yy_xy[i] - rab_z[i] * t_0_xyz_yy_xy[i];

        t_z_xyz_yy_xx[i] = t_0_xyzz_yy_xx[i] - rab_z[i] * t_0_xyz_yy_xx[i];

        t_z_xyz_xz_zz[i] = t_0_xyzz_xz_zz[i] - rab_z[i] * t_0_xyz_xz_zz[i];

        t_z_xyz_xz_yz[i] = t_0_xyzz_xz_yz[i] - rab_z[i] * t_0_xyz_xz_yz[i];

        t_z_xyz_xz_yy[i] = t_0_xyzz_xz_yy[i] - rab_z[i] * t_0_xyz_xz_yy[i];

        t_z_xyz_xz_xz[i] = t_0_xyzz_xz_xz[i] - rab_z[i] * t_0_xyz_xz_xz[i];

        t_z_xyz_xz_xy[i] = t_0_xyzz_xz_xy[i] - rab_z[i] * t_0_xyz_xz_xy[i];

        t_z_xyz_xz_xx[i] = t_0_xyzz_xz_xx[i] - rab_z[i] * t_0_xyz_xz_xx[i];

        t_z_xyz_xy_zz[i] = t_0_xyzz_xy_zz[i] - rab_z[i] * t_0_xyz_xy_zz[i];

        t_z_xyz_xy_yz[i] = t_0_xyzz_xy_yz[i] - rab_z[i] * t_0_xyz_xy_yz[i];

        t_z_xyz_xy_yy[i] = t_0_xyzz_xy_yy[i] - rab_z[i] * t_0_xyz_xy_yy[i];

        t_z_xyz_xy_xz[i] = t_0_xyzz_xy_xz[i] - rab_z[i] * t_0_xyz_xy_xz[i];

        t_z_xyz_xy_xy[i] = t_0_xyzz_xy_xy[i] - rab_z[i] * t_0_xyz_xy_xy[i];

        t_z_xyz_xy_xx[i] = t_0_xyzz_xy_xx[i] - rab_z[i] * t_0_xyz_xy_xx[i];

        t_z_xyz_xx_zz[i] = t_0_xyzz_xx_zz[i] - rab_z[i] * t_0_xyz_xx_zz[i];

        t_z_xyz_xx_yz[i] = t_0_xyzz_xx_yz[i] - rab_z[i] * t_0_xyz_xx_yz[i];

        t_z_xyz_xx_yy[i] = t_0_xyzz_xx_yy[i] - rab_z[i] * t_0_xyz_xx_yy[i];

        t_z_xyz_xx_xz[i] = t_0_xyzz_xx_xz[i] - rab_z[i] * t_0_xyz_xx_xz[i];

        t_z_xyz_xx_xy[i] = t_0_xyzz_xx_xy[i] - rab_z[i] * t_0_xyz_xx_xy[i];

        t_z_xyz_xx_xx[i] = t_0_xyzz_xx_xx[i] - rab_z[i] * t_0_xyz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xxz_xx_xx, t_0_xxz_xx_xy, t_0_xxz_xx_xz, t_0_xxz_xx_yy,\
                           t_0_xxz_xx_yz, t_0_xxz_xx_zz, t_0_xxz_xy_xx, t_0_xxz_xy_xy,\
                           t_0_xxz_xy_xz, t_0_xxz_xy_yy, t_0_xxz_xy_yz, t_0_xxz_xy_zz,\
                           t_0_xxz_xz_xx, t_0_xxz_xz_xy, t_0_xxz_xz_xz, t_0_xxz_xz_yy,\
                           t_0_xxz_xz_yz, t_0_xxz_xz_zz, t_0_xxz_yy_xx, t_0_xxz_yy_xy,\
                           t_0_xxz_yy_xz, t_0_xxz_yy_yy, t_0_xxz_yy_yz, t_0_xxz_yy_zz,\
                           t_0_xxz_yz_xx, t_0_xxz_yz_xy, t_0_xxz_yz_xz, t_0_xxz_yz_yy,\
                           t_0_xxz_yz_yz, t_0_xxz_yz_zz, t_0_xxz_zz_xx, t_0_xxz_zz_xy,\
                           t_0_xxz_zz_xz, t_0_xxz_zz_yy, t_0_xxz_zz_yz, t_0_xxz_zz_zz,\
                           t_0_xxzz_xx_xx, t_0_xxzz_xx_xy, t_0_xxzz_xx_xz, t_0_xxzz_xx_yy,\
                           t_0_xxzz_xx_yz, t_0_xxzz_xx_zz, t_0_xxzz_xy_xx, t_0_xxzz_xy_xy,\
                           t_0_xxzz_xy_xz, t_0_xxzz_xy_yy, t_0_xxzz_xy_yz, t_0_xxzz_xy_zz,\
                           t_0_xxzz_xz_xx, t_0_xxzz_xz_xy, t_0_xxzz_xz_xz, t_0_xxzz_xz_yy,\
                           t_0_xxzz_xz_yz, t_0_xxzz_xz_zz, t_0_xxzz_yy_xx, t_0_xxzz_yy_xy,\
                           t_0_xxzz_yy_xz, t_0_xxzz_yy_yy, t_0_xxzz_yy_yz, t_0_xxzz_yy_zz,\
                           t_0_xxzz_yz_xx, t_0_xxzz_yz_xy, t_0_xxzz_yz_xz, t_0_xxzz_yz_yy,\
                           t_0_xxzz_yz_yz, t_0_xxzz_yz_zz, t_0_xxzz_zz_xx, t_0_xxzz_zz_xy,\
                           t_0_xxzz_zz_xz, t_0_xxzz_zz_yy, t_0_xxzz_zz_yz, t_0_xxzz_zz_zz,\
                           t_z_xxz_xx_xx, t_z_xxz_xx_xy, t_z_xxz_xx_xz, t_z_xxz_xx_yy,\
                           t_z_xxz_xx_yz, t_z_xxz_xx_zz, t_z_xxz_xy_xx, t_z_xxz_xy_xy,\
                           t_z_xxz_xy_xz, t_z_xxz_xy_yy, t_z_xxz_xy_yz, t_z_xxz_xy_zz,\
                           t_z_xxz_xz_xx, t_z_xxz_xz_xy, t_z_xxz_xz_xz, t_z_xxz_xz_yy,\
                           t_z_xxz_xz_yz, t_z_xxz_xz_zz, t_z_xxz_yy_xx, t_z_xxz_yy_xy,\
                           t_z_xxz_yy_xz, t_z_xxz_yy_yy, t_z_xxz_yy_yz, t_z_xxz_yy_zz,\
                           t_z_xxz_yz_xx, t_z_xxz_yz_xy, t_z_xxz_yz_xz, t_z_xxz_yz_yy,\
                           t_z_xxz_yz_yz, t_z_xxz_yz_zz, t_z_xxz_zz_xx, t_z_xxz_zz_xy,\
                           t_z_xxz_zz_xz, t_z_xxz_zz_yy, t_z_xxz_zz_yz, t_z_xxz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xxz_zz_zz[i] = t_0_xxzz_zz_zz[i] - rab_z[i] * t_0_xxz_zz_zz[i];

        t_z_xxz_zz_yz[i] = t_0_xxzz_zz_yz[i] - rab_z[i] * t_0_xxz_zz_yz[i];

        t_z_xxz_zz_yy[i] = t_0_xxzz_zz_yy[i] - rab_z[i] * t_0_xxz_zz_yy[i];

        t_z_xxz_zz_xz[i] = t_0_xxzz_zz_xz[i] - rab_z[i] * t_0_xxz_zz_xz[i];

        t_z_xxz_zz_xy[i] = t_0_xxzz_zz_xy[i] - rab_z[i] * t_0_xxz_zz_xy[i];

        t_z_xxz_zz_xx[i] = t_0_xxzz_zz_xx[i] - rab_z[i] * t_0_xxz_zz_xx[i];

        t_z_xxz_yz_zz[i] = t_0_xxzz_yz_zz[i] - rab_z[i] * t_0_xxz_yz_zz[i];

        t_z_xxz_yz_yz[i] = t_0_xxzz_yz_yz[i] - rab_z[i] * t_0_xxz_yz_yz[i];

        t_z_xxz_yz_yy[i] = t_0_xxzz_yz_yy[i] - rab_z[i] * t_0_xxz_yz_yy[i];

        t_z_xxz_yz_xz[i] = t_0_xxzz_yz_xz[i] - rab_z[i] * t_0_xxz_yz_xz[i];

        t_z_xxz_yz_xy[i] = t_0_xxzz_yz_xy[i] - rab_z[i] * t_0_xxz_yz_xy[i];

        t_z_xxz_yz_xx[i] = t_0_xxzz_yz_xx[i] - rab_z[i] * t_0_xxz_yz_xx[i];

        t_z_xxz_yy_zz[i] = t_0_xxzz_yy_zz[i] - rab_z[i] * t_0_xxz_yy_zz[i];

        t_z_xxz_yy_yz[i] = t_0_xxzz_yy_yz[i] - rab_z[i] * t_0_xxz_yy_yz[i];

        t_z_xxz_yy_yy[i] = t_0_xxzz_yy_yy[i] - rab_z[i] * t_0_xxz_yy_yy[i];

        t_z_xxz_yy_xz[i] = t_0_xxzz_yy_xz[i] - rab_z[i] * t_0_xxz_yy_xz[i];

        t_z_xxz_yy_xy[i] = t_0_xxzz_yy_xy[i] - rab_z[i] * t_0_xxz_yy_xy[i];

        t_z_xxz_yy_xx[i] = t_0_xxzz_yy_xx[i] - rab_z[i] * t_0_xxz_yy_xx[i];

        t_z_xxz_xz_zz[i] = t_0_xxzz_xz_zz[i] - rab_z[i] * t_0_xxz_xz_zz[i];

        t_z_xxz_xz_yz[i] = t_0_xxzz_xz_yz[i] - rab_z[i] * t_0_xxz_xz_yz[i];

        t_z_xxz_xz_yy[i] = t_0_xxzz_xz_yy[i] - rab_z[i] * t_0_xxz_xz_yy[i];

        t_z_xxz_xz_xz[i] = t_0_xxzz_xz_xz[i] - rab_z[i] * t_0_xxz_xz_xz[i];

        t_z_xxz_xz_xy[i] = t_0_xxzz_xz_xy[i] - rab_z[i] * t_0_xxz_xz_xy[i];

        t_z_xxz_xz_xx[i] = t_0_xxzz_xz_xx[i] - rab_z[i] * t_0_xxz_xz_xx[i];

        t_z_xxz_xy_zz[i] = t_0_xxzz_xy_zz[i] - rab_z[i] * t_0_xxz_xy_zz[i];

        t_z_xxz_xy_yz[i] = t_0_xxzz_xy_yz[i] - rab_z[i] * t_0_xxz_xy_yz[i];

        t_z_xxz_xy_yy[i] = t_0_xxzz_xy_yy[i] - rab_z[i] * t_0_xxz_xy_yy[i];

        t_z_xxz_xy_xz[i] = t_0_xxzz_xy_xz[i] - rab_z[i] * t_0_xxz_xy_xz[i];

        t_z_xxz_xy_xy[i] = t_0_xxzz_xy_xy[i] - rab_z[i] * t_0_xxz_xy_xy[i];

        t_z_xxz_xy_xx[i] = t_0_xxzz_xy_xx[i] - rab_z[i] * t_0_xxz_xy_xx[i];

        t_z_xxz_xx_zz[i] = t_0_xxzz_xx_zz[i] - rab_z[i] * t_0_xxz_xx_zz[i];

        t_z_xxz_xx_yz[i] = t_0_xxzz_xx_yz[i] - rab_z[i] * t_0_xxz_xx_yz[i];

        t_z_xxz_xx_yy[i] = t_0_xxzz_xx_yy[i] - rab_z[i] * t_0_xxz_xx_yy[i];

        t_z_xxz_xx_xz[i] = t_0_xxzz_xx_xz[i] - rab_z[i] * t_0_xxz_xx_xz[i];

        t_z_xxz_xx_xy[i] = t_0_xxzz_xx_xy[i] - rab_z[i] * t_0_xxz_xx_xy[i];

        t_z_xxz_xx_xx[i] = t_0_xxzz_xx_xx[i] - rab_z[i] * t_0_xxz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yzzz_xx_xx, t_0_yzzz_xx_xy, t_0_yzzz_xx_xz, t_0_yzzz_xx_yy,\
                           t_0_yzzz_xx_yz, t_0_yzzz_xx_zz, t_0_yzzz_xy_xx, t_0_yzzz_xy_xy,\
                           t_0_yzzz_xy_xz, t_0_yzzz_xy_yy, t_0_yzzz_xy_yz, t_0_yzzz_xy_zz,\
                           t_0_yzzz_xz_xx, t_0_yzzz_xz_xy, t_0_yzzz_xz_xz, t_0_yzzz_xz_yy,\
                           t_0_yzzz_xz_yz, t_0_yzzz_xz_zz, t_0_yzzz_yy_xx, t_0_yzzz_yy_xy,\
                           t_0_yzzz_yy_xz, t_0_yzzz_yy_yy, t_0_yzzz_yy_yz, t_0_yzzz_yy_zz,\
                           t_0_yzzz_yz_xx, t_0_yzzz_yz_xy, t_0_yzzz_yz_xz, t_0_yzzz_yz_yy,\
                           t_0_yzzz_yz_yz, t_0_yzzz_yz_zz, t_0_yzzz_zz_xx, t_0_yzzz_zz_xy,\
                           t_0_yzzz_zz_xz, t_0_yzzz_zz_yy, t_0_yzzz_zz_yz, t_0_yzzz_zz_zz,\
                           t_0_zzz_xx_xx, t_0_zzz_xx_xy, t_0_zzz_xx_xz, t_0_zzz_xx_yy,\
                           t_0_zzz_xx_yz, t_0_zzz_xx_zz, t_0_zzz_xy_xx, t_0_zzz_xy_xy,\
                           t_0_zzz_xy_xz, t_0_zzz_xy_yy, t_0_zzz_xy_yz, t_0_zzz_xy_zz,\
                           t_0_zzz_xz_xx, t_0_zzz_xz_xy, t_0_zzz_xz_xz, t_0_zzz_xz_yy,\
                           t_0_zzz_xz_yz, t_0_zzz_xz_zz, t_0_zzz_yy_xx, t_0_zzz_yy_xy,\
                           t_0_zzz_yy_xz, t_0_zzz_yy_yy, t_0_zzz_yy_yz, t_0_zzz_yy_zz,\
                           t_0_zzz_yz_xx, t_0_zzz_yz_xy, t_0_zzz_yz_xz, t_0_zzz_yz_yy,\
                           t_0_zzz_yz_yz, t_0_zzz_yz_zz, t_0_zzz_zz_xx, t_0_zzz_zz_xy,\
                           t_0_zzz_zz_xz, t_0_zzz_zz_yy, t_0_zzz_zz_yz, t_0_zzz_zz_zz,\
                           t_y_zzz_xx_xx, t_y_zzz_xx_xy, t_y_zzz_xx_xz, t_y_zzz_xx_yy,\
                           t_y_zzz_xx_yz, t_y_zzz_xx_zz, t_y_zzz_xy_xx, t_y_zzz_xy_xy,\
                           t_y_zzz_xy_xz, t_y_zzz_xy_yy, t_y_zzz_xy_yz, t_y_zzz_xy_zz,\
                           t_y_zzz_xz_xx, t_y_zzz_xz_xy, t_y_zzz_xz_xz, t_y_zzz_xz_yy,\
                           t_y_zzz_xz_yz, t_y_zzz_xz_zz, t_y_zzz_yy_xx, t_y_zzz_yy_xy,\
                           t_y_zzz_yy_xz, t_y_zzz_yy_yy, t_y_zzz_yy_yz, t_y_zzz_yy_zz,\
                           t_y_zzz_yz_xx, t_y_zzz_yz_xy, t_y_zzz_yz_xz, t_y_zzz_yz_yy,\
                           t_y_zzz_yz_yz, t_y_zzz_yz_zz, t_y_zzz_zz_xx, t_y_zzz_zz_xy,\
                           t_y_zzz_zz_xz, t_y_zzz_zz_yy, t_y_zzz_zz_yz, t_y_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zzz_zz_zz[i] = t_0_yzzz_zz_zz[i] - rab_y[i] * t_0_zzz_zz_zz[i];

        t_y_zzz_zz_yz[i] = t_0_yzzz_zz_yz[i] - rab_y[i] * t_0_zzz_zz_yz[i];

        t_y_zzz_zz_yy[i] = t_0_yzzz_zz_yy[i] - rab_y[i] * t_0_zzz_zz_yy[i];

        t_y_zzz_zz_xz[i] = t_0_yzzz_zz_xz[i] - rab_y[i] * t_0_zzz_zz_xz[i];

        t_y_zzz_zz_xy[i] = t_0_yzzz_zz_xy[i] - rab_y[i] * t_0_zzz_zz_xy[i];

        t_y_zzz_zz_xx[i] = t_0_yzzz_zz_xx[i] - rab_y[i] * t_0_zzz_zz_xx[i];

        t_y_zzz_yz_zz[i] = t_0_yzzz_yz_zz[i] - rab_y[i] * t_0_zzz_yz_zz[i];

        t_y_zzz_yz_yz[i] = t_0_yzzz_yz_yz[i] - rab_y[i] * t_0_zzz_yz_yz[i];

        t_y_zzz_yz_yy[i] = t_0_yzzz_yz_yy[i] - rab_y[i] * t_0_zzz_yz_yy[i];

        t_y_zzz_yz_xz[i] = t_0_yzzz_yz_xz[i] - rab_y[i] * t_0_zzz_yz_xz[i];

        t_y_zzz_yz_xy[i] = t_0_yzzz_yz_xy[i] - rab_y[i] * t_0_zzz_yz_xy[i];

        t_y_zzz_yz_xx[i] = t_0_yzzz_yz_xx[i] - rab_y[i] * t_0_zzz_yz_xx[i];

        t_y_zzz_yy_zz[i] = t_0_yzzz_yy_zz[i] - rab_y[i] * t_0_zzz_yy_zz[i];

        t_y_zzz_yy_yz[i] = t_0_yzzz_yy_yz[i] - rab_y[i] * t_0_zzz_yy_yz[i];

        t_y_zzz_yy_yy[i] = t_0_yzzz_yy_yy[i] - rab_y[i] * t_0_zzz_yy_yy[i];

        t_y_zzz_yy_xz[i] = t_0_yzzz_yy_xz[i] - rab_y[i] * t_0_zzz_yy_xz[i];

        t_y_zzz_yy_xy[i] = t_0_yzzz_yy_xy[i] - rab_y[i] * t_0_zzz_yy_xy[i];

        t_y_zzz_yy_xx[i] = t_0_yzzz_yy_xx[i] - rab_y[i] * t_0_zzz_yy_xx[i];

        t_y_zzz_xz_zz[i] = t_0_yzzz_xz_zz[i] - rab_y[i] * t_0_zzz_xz_zz[i];

        t_y_zzz_xz_yz[i] = t_0_yzzz_xz_yz[i] - rab_y[i] * t_0_zzz_xz_yz[i];

        t_y_zzz_xz_yy[i] = t_0_yzzz_xz_yy[i] - rab_y[i] * t_0_zzz_xz_yy[i];

        t_y_zzz_xz_xz[i] = t_0_yzzz_xz_xz[i] - rab_y[i] * t_0_zzz_xz_xz[i];

        t_y_zzz_xz_xy[i] = t_0_yzzz_xz_xy[i] - rab_y[i] * t_0_zzz_xz_xy[i];

        t_y_zzz_xz_xx[i] = t_0_yzzz_xz_xx[i] - rab_y[i] * t_0_zzz_xz_xx[i];

        t_y_zzz_xy_zz[i] = t_0_yzzz_xy_zz[i] - rab_y[i] * t_0_zzz_xy_zz[i];

        t_y_zzz_xy_yz[i] = t_0_yzzz_xy_yz[i] - rab_y[i] * t_0_zzz_xy_yz[i];

        t_y_zzz_xy_yy[i] = t_0_yzzz_xy_yy[i] - rab_y[i] * t_0_zzz_xy_yy[i];

        t_y_zzz_xy_xz[i] = t_0_yzzz_xy_xz[i] - rab_y[i] * t_0_zzz_xy_xz[i];

        t_y_zzz_xy_xy[i] = t_0_yzzz_xy_xy[i] - rab_y[i] * t_0_zzz_xy_xy[i];

        t_y_zzz_xy_xx[i] = t_0_yzzz_xy_xx[i] - rab_y[i] * t_0_zzz_xy_xx[i];

        t_y_zzz_xx_zz[i] = t_0_yzzz_xx_zz[i] - rab_y[i] * t_0_zzz_xx_zz[i];

        t_y_zzz_xx_yz[i] = t_0_yzzz_xx_yz[i] - rab_y[i] * t_0_zzz_xx_yz[i];

        t_y_zzz_xx_yy[i] = t_0_yzzz_xx_yy[i] - rab_y[i] * t_0_zzz_xx_yy[i];

        t_y_zzz_xx_xz[i] = t_0_yzzz_xx_xz[i] - rab_y[i] * t_0_zzz_xx_xz[i];

        t_y_zzz_xx_xy[i] = t_0_yzzz_xx_xy[i] - rab_y[i] * t_0_zzz_xx_xy[i];

        t_y_zzz_xx_xx[i] = t_0_yzzz_xx_xx[i] - rab_y[i] * t_0_zzz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyzz_xx_xx, t_0_yyzz_xx_xy, t_0_yyzz_xx_xz, t_0_yyzz_xx_yy,\
                           t_0_yyzz_xx_yz, t_0_yyzz_xx_zz, t_0_yyzz_xy_xx, t_0_yyzz_xy_xy,\
                           t_0_yyzz_xy_xz, t_0_yyzz_xy_yy, t_0_yyzz_xy_yz, t_0_yyzz_xy_zz,\
                           t_0_yyzz_xz_xx, t_0_yyzz_xz_xy, t_0_yyzz_xz_xz, t_0_yyzz_xz_yy,\
                           t_0_yyzz_xz_yz, t_0_yyzz_xz_zz, t_0_yyzz_yy_xx, t_0_yyzz_yy_xy,\
                           t_0_yyzz_yy_xz, t_0_yyzz_yy_yy, t_0_yyzz_yy_yz, t_0_yyzz_yy_zz,\
                           t_0_yyzz_yz_xx, t_0_yyzz_yz_xy, t_0_yyzz_yz_xz, t_0_yyzz_yz_yy,\
                           t_0_yyzz_yz_yz, t_0_yyzz_yz_zz, t_0_yyzz_zz_xx, t_0_yyzz_zz_xy,\
                           t_0_yyzz_zz_xz, t_0_yyzz_zz_yy, t_0_yyzz_zz_yz, t_0_yyzz_zz_zz,\
                           t_0_yzz_xx_xx, t_0_yzz_xx_xy, t_0_yzz_xx_xz, t_0_yzz_xx_yy,\
                           t_0_yzz_xx_yz, t_0_yzz_xx_zz, t_0_yzz_xy_xx, t_0_yzz_xy_xy,\
                           t_0_yzz_xy_xz, t_0_yzz_xy_yy, t_0_yzz_xy_yz, t_0_yzz_xy_zz,\
                           t_0_yzz_xz_xx, t_0_yzz_xz_xy, t_0_yzz_xz_xz, t_0_yzz_xz_yy,\
                           t_0_yzz_xz_yz, t_0_yzz_xz_zz, t_0_yzz_yy_xx, t_0_yzz_yy_xy,\
                           t_0_yzz_yy_xz, t_0_yzz_yy_yy, t_0_yzz_yy_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_xx, t_0_yzz_yz_xy, t_0_yzz_yz_xz, t_0_yzz_yz_yy,\
                           t_0_yzz_yz_yz, t_0_yzz_yz_zz, t_0_yzz_zz_xx, t_0_yzz_zz_xy,\
                           t_0_yzz_zz_xz, t_0_yzz_zz_yy, t_0_yzz_zz_yz, t_0_yzz_zz_zz,\
                           t_y_yzz_xx_xx, t_y_yzz_xx_xy, t_y_yzz_xx_xz, t_y_yzz_xx_yy,\
                           t_y_yzz_xx_yz, t_y_yzz_xx_zz, t_y_yzz_xy_xx, t_y_yzz_xy_xy,\
                           t_y_yzz_xy_xz, t_y_yzz_xy_yy, t_y_yzz_xy_yz, t_y_yzz_xy_zz,\
                           t_y_yzz_xz_xx, t_y_yzz_xz_xy, t_y_yzz_xz_xz, t_y_yzz_xz_yy,\
                           t_y_yzz_xz_yz, t_y_yzz_xz_zz, t_y_yzz_yy_xx, t_y_yzz_yy_xy,\
                           t_y_yzz_yy_xz, t_y_yzz_yy_yy, t_y_yzz_yy_yz, t_y_yzz_yy_zz,\
                           t_y_yzz_yz_xx, t_y_yzz_yz_xy, t_y_yzz_yz_xz, t_y_yzz_yz_yy,\
                           t_y_yzz_yz_yz, t_y_yzz_yz_zz, t_y_yzz_zz_xx, t_y_yzz_zz_xy,\
                           t_y_yzz_zz_xz, t_y_yzz_zz_yy, t_y_yzz_zz_yz, t_y_yzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yzz_zz_zz[i] = t_0_yyzz_zz_zz[i] - rab_y[i] * t_0_yzz_zz_zz[i];

        t_y_yzz_zz_yz[i] = t_0_yyzz_zz_yz[i] - rab_y[i] * t_0_yzz_zz_yz[i];

        t_y_yzz_zz_yy[i] = t_0_yyzz_zz_yy[i] - rab_y[i] * t_0_yzz_zz_yy[i];

        t_y_yzz_zz_xz[i] = t_0_yyzz_zz_xz[i] - rab_y[i] * t_0_yzz_zz_xz[i];

        t_y_yzz_zz_xy[i] = t_0_yyzz_zz_xy[i] - rab_y[i] * t_0_yzz_zz_xy[i];

        t_y_yzz_zz_xx[i] = t_0_yyzz_zz_xx[i] - rab_y[i] * t_0_yzz_zz_xx[i];

        t_y_yzz_yz_zz[i] = t_0_yyzz_yz_zz[i] - rab_y[i] * t_0_yzz_yz_zz[i];

        t_y_yzz_yz_yz[i] = t_0_yyzz_yz_yz[i] - rab_y[i] * t_0_yzz_yz_yz[i];

        t_y_yzz_yz_yy[i] = t_0_yyzz_yz_yy[i] - rab_y[i] * t_0_yzz_yz_yy[i];

        t_y_yzz_yz_xz[i] = t_0_yyzz_yz_xz[i] - rab_y[i] * t_0_yzz_yz_xz[i];

        t_y_yzz_yz_xy[i] = t_0_yyzz_yz_xy[i] - rab_y[i] * t_0_yzz_yz_xy[i];

        t_y_yzz_yz_xx[i] = t_0_yyzz_yz_xx[i] - rab_y[i] * t_0_yzz_yz_xx[i];

        t_y_yzz_yy_zz[i] = t_0_yyzz_yy_zz[i] - rab_y[i] * t_0_yzz_yy_zz[i];

        t_y_yzz_yy_yz[i] = t_0_yyzz_yy_yz[i] - rab_y[i] * t_0_yzz_yy_yz[i];

        t_y_yzz_yy_yy[i] = t_0_yyzz_yy_yy[i] - rab_y[i] * t_0_yzz_yy_yy[i];

        t_y_yzz_yy_xz[i] = t_0_yyzz_yy_xz[i] - rab_y[i] * t_0_yzz_yy_xz[i];

        t_y_yzz_yy_xy[i] = t_0_yyzz_yy_xy[i] - rab_y[i] * t_0_yzz_yy_xy[i];

        t_y_yzz_yy_xx[i] = t_0_yyzz_yy_xx[i] - rab_y[i] * t_0_yzz_yy_xx[i];

        t_y_yzz_xz_zz[i] = t_0_yyzz_xz_zz[i] - rab_y[i] * t_0_yzz_xz_zz[i];

        t_y_yzz_xz_yz[i] = t_0_yyzz_xz_yz[i] - rab_y[i] * t_0_yzz_xz_yz[i];

        t_y_yzz_xz_yy[i] = t_0_yyzz_xz_yy[i] - rab_y[i] * t_0_yzz_xz_yy[i];

        t_y_yzz_xz_xz[i] = t_0_yyzz_xz_xz[i] - rab_y[i] * t_0_yzz_xz_xz[i];

        t_y_yzz_xz_xy[i] = t_0_yyzz_xz_xy[i] - rab_y[i] * t_0_yzz_xz_xy[i];

        t_y_yzz_xz_xx[i] = t_0_yyzz_xz_xx[i] - rab_y[i] * t_0_yzz_xz_xx[i];

        t_y_yzz_xy_zz[i] = t_0_yyzz_xy_zz[i] - rab_y[i] * t_0_yzz_xy_zz[i];

        t_y_yzz_xy_yz[i] = t_0_yyzz_xy_yz[i] - rab_y[i] * t_0_yzz_xy_yz[i];

        t_y_yzz_xy_yy[i] = t_0_yyzz_xy_yy[i] - rab_y[i] * t_0_yzz_xy_yy[i];

        t_y_yzz_xy_xz[i] = t_0_yyzz_xy_xz[i] - rab_y[i] * t_0_yzz_xy_xz[i];

        t_y_yzz_xy_xy[i] = t_0_yyzz_xy_xy[i] - rab_y[i] * t_0_yzz_xy_xy[i];

        t_y_yzz_xy_xx[i] = t_0_yyzz_xy_xx[i] - rab_y[i] * t_0_yzz_xy_xx[i];

        t_y_yzz_xx_zz[i] = t_0_yyzz_xx_zz[i] - rab_y[i] * t_0_yzz_xx_zz[i];

        t_y_yzz_xx_yz[i] = t_0_yyzz_xx_yz[i] - rab_y[i] * t_0_yzz_xx_yz[i];

        t_y_yzz_xx_yy[i] = t_0_yyzz_xx_yy[i] - rab_y[i] * t_0_yzz_xx_yy[i];

        t_y_yzz_xx_xz[i] = t_0_yyzz_xx_xz[i] - rab_y[i] * t_0_yzz_xx_xz[i];

        t_y_yzz_xx_xy[i] = t_0_yyzz_xx_xy[i] - rab_y[i] * t_0_yzz_xx_xy[i];

        t_y_yzz_xx_xx[i] = t_0_yyzz_xx_xx[i] - rab_y[i] * t_0_yzz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyyz_xx_xx, t_0_yyyz_xx_xy, t_0_yyyz_xx_xz, t_0_yyyz_xx_yy,\
                           t_0_yyyz_xx_yz, t_0_yyyz_xx_zz, t_0_yyyz_xy_xx, t_0_yyyz_xy_xy,\
                           t_0_yyyz_xy_xz, t_0_yyyz_xy_yy, t_0_yyyz_xy_yz, t_0_yyyz_xy_zz,\
                           t_0_yyyz_xz_xx, t_0_yyyz_xz_xy, t_0_yyyz_xz_xz, t_0_yyyz_xz_yy,\
                           t_0_yyyz_xz_yz, t_0_yyyz_xz_zz, t_0_yyyz_yy_xx, t_0_yyyz_yy_xy,\
                           t_0_yyyz_yy_xz, t_0_yyyz_yy_yy, t_0_yyyz_yy_yz, t_0_yyyz_yy_zz,\
                           t_0_yyyz_yz_xx, t_0_yyyz_yz_xy, t_0_yyyz_yz_xz, t_0_yyyz_yz_yy,\
                           t_0_yyyz_yz_yz, t_0_yyyz_yz_zz, t_0_yyyz_zz_xx, t_0_yyyz_zz_xy,\
                           t_0_yyyz_zz_xz, t_0_yyyz_zz_yy, t_0_yyyz_zz_yz, t_0_yyyz_zz_zz,\
                           t_0_yyz_xx_xx, t_0_yyz_xx_xy, t_0_yyz_xx_xz, t_0_yyz_xx_yy,\
                           t_0_yyz_xx_yz, t_0_yyz_xx_zz, t_0_yyz_xy_xx, t_0_yyz_xy_xy,\
                           t_0_yyz_xy_xz, t_0_yyz_xy_yy, t_0_yyz_xy_yz, t_0_yyz_xy_zz,\
                           t_0_yyz_xz_xx, t_0_yyz_xz_xy, t_0_yyz_xz_xz, t_0_yyz_xz_yy,\
                           t_0_yyz_xz_yz, t_0_yyz_xz_zz, t_0_yyz_yy_xx, t_0_yyz_yy_xy,\
                           t_0_yyz_yy_xz, t_0_yyz_yy_yy, t_0_yyz_yy_yz, t_0_yyz_yy_zz,\
                           t_0_yyz_yz_xx, t_0_yyz_yz_xy, t_0_yyz_yz_xz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_yz_zz, t_0_yyz_zz_xx, t_0_yyz_zz_xy,\
                           t_0_yyz_zz_xz, t_0_yyz_zz_yy, t_0_yyz_zz_yz, t_0_yyz_zz_zz,\
                           t_y_yyz_xx_xx, t_y_yyz_xx_xy, t_y_yyz_xx_xz, t_y_yyz_xx_yy,\
                           t_y_yyz_xx_yz, t_y_yyz_xx_zz, t_y_yyz_xy_xx, t_y_yyz_xy_xy,\
                           t_y_yyz_xy_xz, t_y_yyz_xy_yy, t_y_yyz_xy_yz, t_y_yyz_xy_zz,\
                           t_y_yyz_xz_xx, t_y_yyz_xz_xy, t_y_yyz_xz_xz, t_y_yyz_xz_yy,\
                           t_y_yyz_xz_yz, t_y_yyz_xz_zz, t_y_yyz_yy_xx, t_y_yyz_yy_xy,\
                           t_y_yyz_yy_xz, t_y_yyz_yy_yy, t_y_yyz_yy_yz, t_y_yyz_yy_zz,\
                           t_y_yyz_yz_xx, t_y_yyz_yz_xy, t_y_yyz_yz_xz, t_y_yyz_yz_yy,\
                           t_y_yyz_yz_yz, t_y_yyz_yz_zz, t_y_yyz_zz_xx, t_y_yyz_zz_xy,\
                           t_y_yyz_zz_xz, t_y_yyz_zz_yy, t_y_yyz_zz_yz, t_y_yyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yyz_zz_zz[i] = t_0_yyyz_zz_zz[i] - rab_y[i] * t_0_yyz_zz_zz[i];

        t_y_yyz_zz_yz[i] = t_0_yyyz_zz_yz[i] - rab_y[i] * t_0_yyz_zz_yz[i];

        t_y_yyz_zz_yy[i] = t_0_yyyz_zz_yy[i] - rab_y[i] * t_0_yyz_zz_yy[i];

        t_y_yyz_zz_xz[i] = t_0_yyyz_zz_xz[i] - rab_y[i] * t_0_yyz_zz_xz[i];

        t_y_yyz_zz_xy[i] = t_0_yyyz_zz_xy[i] - rab_y[i] * t_0_yyz_zz_xy[i];

        t_y_yyz_zz_xx[i] = t_0_yyyz_zz_xx[i] - rab_y[i] * t_0_yyz_zz_xx[i];

        t_y_yyz_yz_zz[i] = t_0_yyyz_yz_zz[i] - rab_y[i] * t_0_yyz_yz_zz[i];

        t_y_yyz_yz_yz[i] = t_0_yyyz_yz_yz[i] - rab_y[i] * t_0_yyz_yz_yz[i];

        t_y_yyz_yz_yy[i] = t_0_yyyz_yz_yy[i] - rab_y[i] * t_0_yyz_yz_yy[i];

        t_y_yyz_yz_xz[i] = t_0_yyyz_yz_xz[i] - rab_y[i] * t_0_yyz_yz_xz[i];

        t_y_yyz_yz_xy[i] = t_0_yyyz_yz_xy[i] - rab_y[i] * t_0_yyz_yz_xy[i];

        t_y_yyz_yz_xx[i] = t_0_yyyz_yz_xx[i] - rab_y[i] * t_0_yyz_yz_xx[i];

        t_y_yyz_yy_zz[i] = t_0_yyyz_yy_zz[i] - rab_y[i] * t_0_yyz_yy_zz[i];

        t_y_yyz_yy_yz[i] = t_0_yyyz_yy_yz[i] - rab_y[i] * t_0_yyz_yy_yz[i];

        t_y_yyz_yy_yy[i] = t_0_yyyz_yy_yy[i] - rab_y[i] * t_0_yyz_yy_yy[i];

        t_y_yyz_yy_xz[i] = t_0_yyyz_yy_xz[i] - rab_y[i] * t_0_yyz_yy_xz[i];

        t_y_yyz_yy_xy[i] = t_0_yyyz_yy_xy[i] - rab_y[i] * t_0_yyz_yy_xy[i];

        t_y_yyz_yy_xx[i] = t_0_yyyz_yy_xx[i] - rab_y[i] * t_0_yyz_yy_xx[i];

        t_y_yyz_xz_zz[i] = t_0_yyyz_xz_zz[i] - rab_y[i] * t_0_yyz_xz_zz[i];

        t_y_yyz_xz_yz[i] = t_0_yyyz_xz_yz[i] - rab_y[i] * t_0_yyz_xz_yz[i];

        t_y_yyz_xz_yy[i] = t_0_yyyz_xz_yy[i] - rab_y[i] * t_0_yyz_xz_yy[i];

        t_y_yyz_xz_xz[i] = t_0_yyyz_xz_xz[i] - rab_y[i] * t_0_yyz_xz_xz[i];

        t_y_yyz_xz_xy[i] = t_0_yyyz_xz_xy[i] - rab_y[i] * t_0_yyz_xz_xy[i];

        t_y_yyz_xz_xx[i] = t_0_yyyz_xz_xx[i] - rab_y[i] * t_0_yyz_xz_xx[i];

        t_y_yyz_xy_zz[i] = t_0_yyyz_xy_zz[i] - rab_y[i] * t_0_yyz_xy_zz[i];

        t_y_yyz_xy_yz[i] = t_0_yyyz_xy_yz[i] - rab_y[i] * t_0_yyz_xy_yz[i];

        t_y_yyz_xy_yy[i] = t_0_yyyz_xy_yy[i] - rab_y[i] * t_0_yyz_xy_yy[i];

        t_y_yyz_xy_xz[i] = t_0_yyyz_xy_xz[i] - rab_y[i] * t_0_yyz_xy_xz[i];

        t_y_yyz_xy_xy[i] = t_0_yyyz_xy_xy[i] - rab_y[i] * t_0_yyz_xy_xy[i];

        t_y_yyz_xy_xx[i] = t_0_yyyz_xy_xx[i] - rab_y[i] * t_0_yyz_xy_xx[i];

        t_y_yyz_xx_zz[i] = t_0_yyyz_xx_zz[i] - rab_y[i] * t_0_yyz_xx_zz[i];

        t_y_yyz_xx_yz[i] = t_0_yyyz_xx_yz[i] - rab_y[i] * t_0_yyz_xx_yz[i];

        t_y_yyz_xx_yy[i] = t_0_yyyz_xx_yy[i] - rab_y[i] * t_0_yyz_xx_yy[i];

        t_y_yyz_xx_xz[i] = t_0_yyyz_xx_xz[i] - rab_y[i] * t_0_yyz_xx_xz[i];

        t_y_yyz_xx_xy[i] = t_0_yyyz_xx_xy[i] - rab_y[i] * t_0_yyz_xx_xy[i];

        t_y_yyz_xx_xx[i] = t_0_yyyz_xx_xx[i] - rab_y[i] * t_0_yyz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyy_xx_xx, t_0_yyy_xx_xy, t_0_yyy_xx_xz, t_0_yyy_xx_yy,\
                           t_0_yyy_xx_yz, t_0_yyy_xx_zz, t_0_yyy_xy_xx, t_0_yyy_xy_xy,\
                           t_0_yyy_xy_xz, t_0_yyy_xy_yy, t_0_yyy_xy_yz, t_0_yyy_xy_zz,\
                           t_0_yyy_xz_xx, t_0_yyy_xz_xy, t_0_yyy_xz_xz, t_0_yyy_xz_yy,\
                           t_0_yyy_xz_yz, t_0_yyy_xz_zz, t_0_yyy_yy_xx, t_0_yyy_yy_xy,\
                           t_0_yyy_yy_xz, t_0_yyy_yy_yy, t_0_yyy_yy_yz, t_0_yyy_yy_zz,\
                           t_0_yyy_yz_xx, t_0_yyy_yz_xy, t_0_yyy_yz_xz, t_0_yyy_yz_yy,\
                           t_0_yyy_yz_yz, t_0_yyy_yz_zz, t_0_yyy_zz_xx, t_0_yyy_zz_xy,\
                           t_0_yyy_zz_xz, t_0_yyy_zz_yy, t_0_yyy_zz_yz, t_0_yyy_zz_zz,\
                           t_0_yyyy_xx_xx, t_0_yyyy_xx_xy, t_0_yyyy_xx_xz, t_0_yyyy_xx_yy,\
                           t_0_yyyy_xx_yz, t_0_yyyy_xx_zz, t_0_yyyy_xy_xx, t_0_yyyy_xy_xy,\
                           t_0_yyyy_xy_xz, t_0_yyyy_xy_yy, t_0_yyyy_xy_yz, t_0_yyyy_xy_zz,\
                           t_0_yyyy_xz_xx, t_0_yyyy_xz_xy, t_0_yyyy_xz_xz, t_0_yyyy_xz_yy,\
                           t_0_yyyy_xz_yz, t_0_yyyy_xz_zz, t_0_yyyy_yy_xx, t_0_yyyy_yy_xy,\
                           t_0_yyyy_yy_xz, t_0_yyyy_yy_yy, t_0_yyyy_yy_yz, t_0_yyyy_yy_zz,\
                           t_0_yyyy_yz_xx, t_0_yyyy_yz_xy, t_0_yyyy_yz_xz, t_0_yyyy_yz_yy,\
                           t_0_yyyy_yz_yz, t_0_yyyy_yz_zz, t_0_yyyy_zz_xx, t_0_yyyy_zz_xy,\
                           t_0_yyyy_zz_xz, t_0_yyyy_zz_yy, t_0_yyyy_zz_yz, t_0_yyyy_zz_zz,\
                           t_y_yyy_xx_xx, t_y_yyy_xx_xy, t_y_yyy_xx_xz, t_y_yyy_xx_yy,\
                           t_y_yyy_xx_yz, t_y_yyy_xx_zz, t_y_yyy_xy_xx, t_y_yyy_xy_xy,\
                           t_y_yyy_xy_xz, t_y_yyy_xy_yy, t_y_yyy_xy_yz, t_y_yyy_xy_zz,\
                           t_y_yyy_xz_xx, t_y_yyy_xz_xy, t_y_yyy_xz_xz, t_y_yyy_xz_yy,\
                           t_y_yyy_xz_yz, t_y_yyy_xz_zz, t_y_yyy_yy_xx, t_y_yyy_yy_xy,\
                           t_y_yyy_yy_xz, t_y_yyy_yy_yy, t_y_yyy_yy_yz, t_y_yyy_yy_zz,\
                           t_y_yyy_yz_xx, t_y_yyy_yz_xy, t_y_yyy_yz_xz, t_y_yyy_yz_yy,\
                           t_y_yyy_yz_yz, t_y_yyy_yz_zz, t_y_yyy_zz_xx, t_y_yyy_zz_xy,\
                           t_y_yyy_zz_xz, t_y_yyy_zz_yy, t_y_yyy_zz_yz, t_y_yyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yyy_zz_zz[i] = t_0_yyyy_zz_zz[i] - rab_y[i] * t_0_yyy_zz_zz[i];

        t_y_yyy_zz_yz[i] = t_0_yyyy_zz_yz[i] - rab_y[i] * t_0_yyy_zz_yz[i];

        t_y_yyy_zz_yy[i] = t_0_yyyy_zz_yy[i] - rab_y[i] * t_0_yyy_zz_yy[i];

        t_y_yyy_zz_xz[i] = t_0_yyyy_zz_xz[i] - rab_y[i] * t_0_yyy_zz_xz[i];

        t_y_yyy_zz_xy[i] = t_0_yyyy_zz_xy[i] - rab_y[i] * t_0_yyy_zz_xy[i];

        t_y_yyy_zz_xx[i] = t_0_yyyy_zz_xx[i] - rab_y[i] * t_0_yyy_zz_xx[i];

        t_y_yyy_yz_zz[i] = t_0_yyyy_yz_zz[i] - rab_y[i] * t_0_yyy_yz_zz[i];

        t_y_yyy_yz_yz[i] = t_0_yyyy_yz_yz[i] - rab_y[i] * t_0_yyy_yz_yz[i];

        t_y_yyy_yz_yy[i] = t_0_yyyy_yz_yy[i] - rab_y[i] * t_0_yyy_yz_yy[i];

        t_y_yyy_yz_xz[i] = t_0_yyyy_yz_xz[i] - rab_y[i] * t_0_yyy_yz_xz[i];

        t_y_yyy_yz_xy[i] = t_0_yyyy_yz_xy[i] - rab_y[i] * t_0_yyy_yz_xy[i];

        t_y_yyy_yz_xx[i] = t_0_yyyy_yz_xx[i] - rab_y[i] * t_0_yyy_yz_xx[i];

        t_y_yyy_yy_zz[i] = t_0_yyyy_yy_zz[i] - rab_y[i] * t_0_yyy_yy_zz[i];

        t_y_yyy_yy_yz[i] = t_0_yyyy_yy_yz[i] - rab_y[i] * t_0_yyy_yy_yz[i];

        t_y_yyy_yy_yy[i] = t_0_yyyy_yy_yy[i] - rab_y[i] * t_0_yyy_yy_yy[i];

        t_y_yyy_yy_xz[i] = t_0_yyyy_yy_xz[i] - rab_y[i] * t_0_yyy_yy_xz[i];

        t_y_yyy_yy_xy[i] = t_0_yyyy_yy_xy[i] - rab_y[i] * t_0_yyy_yy_xy[i];

        t_y_yyy_yy_xx[i] = t_0_yyyy_yy_xx[i] - rab_y[i] * t_0_yyy_yy_xx[i];

        t_y_yyy_xz_zz[i] = t_0_yyyy_xz_zz[i] - rab_y[i] * t_0_yyy_xz_zz[i];

        t_y_yyy_xz_yz[i] = t_0_yyyy_xz_yz[i] - rab_y[i] * t_0_yyy_xz_yz[i];

        t_y_yyy_xz_yy[i] = t_0_yyyy_xz_yy[i] - rab_y[i] * t_0_yyy_xz_yy[i];

        t_y_yyy_xz_xz[i] = t_0_yyyy_xz_xz[i] - rab_y[i] * t_0_yyy_xz_xz[i];

        t_y_yyy_xz_xy[i] = t_0_yyyy_xz_xy[i] - rab_y[i] * t_0_yyy_xz_xy[i];

        t_y_yyy_xz_xx[i] = t_0_yyyy_xz_xx[i] - rab_y[i] * t_0_yyy_xz_xx[i];

        t_y_yyy_xy_zz[i] = t_0_yyyy_xy_zz[i] - rab_y[i] * t_0_yyy_xy_zz[i];

        t_y_yyy_xy_yz[i] = t_0_yyyy_xy_yz[i] - rab_y[i] * t_0_yyy_xy_yz[i];

        t_y_yyy_xy_yy[i] = t_0_yyyy_xy_yy[i] - rab_y[i] * t_0_yyy_xy_yy[i];

        t_y_yyy_xy_xz[i] = t_0_yyyy_xy_xz[i] - rab_y[i] * t_0_yyy_xy_xz[i];

        t_y_yyy_xy_xy[i] = t_0_yyyy_xy_xy[i] - rab_y[i] * t_0_yyy_xy_xy[i];

        t_y_yyy_xy_xx[i] = t_0_yyyy_xy_xx[i] - rab_y[i] * t_0_yyy_xy_xx[i];

        t_y_yyy_xx_zz[i] = t_0_yyyy_xx_zz[i] - rab_y[i] * t_0_yyy_xx_zz[i];

        t_y_yyy_xx_yz[i] = t_0_yyyy_xx_yz[i] - rab_y[i] * t_0_yyy_xx_yz[i];

        t_y_yyy_xx_yy[i] = t_0_yyyy_xx_yy[i] - rab_y[i] * t_0_yyy_xx_yy[i];

        t_y_yyy_xx_xz[i] = t_0_yyyy_xx_xz[i] - rab_y[i] * t_0_yyy_xx_xz[i];

        t_y_yyy_xx_xy[i] = t_0_yyyy_xx_xy[i] - rab_y[i] * t_0_yyy_xx_xy[i];

        t_y_yyy_xx_xx[i] = t_0_yyyy_xx_xx[i] - rab_y[i] * t_0_yyy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyzz_xx_xx, t_0_xyzz_xx_xy, t_0_xyzz_xx_xz, t_0_xyzz_xx_yy,\
                           t_0_xyzz_xx_yz, t_0_xyzz_xx_zz, t_0_xyzz_xy_xx, t_0_xyzz_xy_xy,\
                           t_0_xyzz_xy_xz, t_0_xyzz_xy_yy, t_0_xyzz_xy_yz, t_0_xyzz_xy_zz,\
                           t_0_xyzz_xz_xx, t_0_xyzz_xz_xy, t_0_xyzz_xz_xz, t_0_xyzz_xz_yy,\
                           t_0_xyzz_xz_yz, t_0_xyzz_xz_zz, t_0_xyzz_yy_xx, t_0_xyzz_yy_xy,\
                           t_0_xyzz_yy_xz, t_0_xyzz_yy_yy, t_0_xyzz_yy_yz, t_0_xyzz_yy_zz,\
                           t_0_xyzz_yz_xx, t_0_xyzz_yz_xy, t_0_xyzz_yz_xz, t_0_xyzz_yz_yy,\
                           t_0_xyzz_yz_yz, t_0_xyzz_yz_zz, t_0_xyzz_zz_xx, t_0_xyzz_zz_xy,\
                           t_0_xyzz_zz_xz, t_0_xyzz_zz_yy, t_0_xyzz_zz_yz, t_0_xyzz_zz_zz,\
                           t_0_xzz_xx_xx, t_0_xzz_xx_xy, t_0_xzz_xx_xz, t_0_xzz_xx_yy,\
                           t_0_xzz_xx_yz, t_0_xzz_xx_zz, t_0_xzz_xy_xx, t_0_xzz_xy_xy,\
                           t_0_xzz_xy_xz, t_0_xzz_xy_yy, t_0_xzz_xy_yz, t_0_xzz_xy_zz,\
                           t_0_xzz_xz_xx, t_0_xzz_xz_xy, t_0_xzz_xz_xz, t_0_xzz_xz_yy,\
                           t_0_xzz_xz_yz, t_0_xzz_xz_zz, t_0_xzz_yy_xx, t_0_xzz_yy_xy,\
                           t_0_xzz_yy_xz, t_0_xzz_yy_yy, t_0_xzz_yy_yz, t_0_xzz_yy_zz,\
                           t_0_xzz_yz_xx, t_0_xzz_yz_xy, t_0_xzz_yz_xz, t_0_xzz_yz_yy,\
                           t_0_xzz_yz_yz, t_0_xzz_yz_zz, t_0_xzz_zz_xx, t_0_xzz_zz_xy,\
                           t_0_xzz_zz_xz, t_0_xzz_zz_yy, t_0_xzz_zz_yz, t_0_xzz_zz_zz,\
                           t_y_xzz_xx_xx, t_y_xzz_xx_xy, t_y_xzz_xx_xz, t_y_xzz_xx_yy,\
                           t_y_xzz_xx_yz, t_y_xzz_xx_zz, t_y_xzz_xy_xx, t_y_xzz_xy_xy,\
                           t_y_xzz_xy_xz, t_y_xzz_xy_yy, t_y_xzz_xy_yz, t_y_xzz_xy_zz,\
                           t_y_xzz_xz_xx, t_y_xzz_xz_xy, t_y_xzz_xz_xz, t_y_xzz_xz_yy,\
                           t_y_xzz_xz_yz, t_y_xzz_xz_zz, t_y_xzz_yy_xx, t_y_xzz_yy_xy,\
                           t_y_xzz_yy_xz, t_y_xzz_yy_yy, t_y_xzz_yy_yz, t_y_xzz_yy_zz,\
                           t_y_xzz_yz_xx, t_y_xzz_yz_xy, t_y_xzz_yz_xz, t_y_xzz_yz_yy,\
                           t_y_xzz_yz_yz, t_y_xzz_yz_zz, t_y_xzz_zz_xx, t_y_xzz_zz_xy,\
                           t_y_xzz_zz_xz, t_y_xzz_zz_yy, t_y_xzz_zz_yz, t_y_xzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xzz_zz_zz[i] = t_0_xyzz_zz_zz[i] - rab_y[i] * t_0_xzz_zz_zz[i];

        t_y_xzz_zz_yz[i] = t_0_xyzz_zz_yz[i] - rab_y[i] * t_0_xzz_zz_yz[i];

        t_y_xzz_zz_yy[i] = t_0_xyzz_zz_yy[i] - rab_y[i] * t_0_xzz_zz_yy[i];

        t_y_xzz_zz_xz[i] = t_0_xyzz_zz_xz[i] - rab_y[i] * t_0_xzz_zz_xz[i];

        t_y_xzz_zz_xy[i] = t_0_xyzz_zz_xy[i] - rab_y[i] * t_0_xzz_zz_xy[i];

        t_y_xzz_zz_xx[i] = t_0_xyzz_zz_xx[i] - rab_y[i] * t_0_xzz_zz_xx[i];

        t_y_xzz_yz_zz[i] = t_0_xyzz_yz_zz[i] - rab_y[i] * t_0_xzz_yz_zz[i];

        t_y_xzz_yz_yz[i] = t_0_xyzz_yz_yz[i] - rab_y[i] * t_0_xzz_yz_yz[i];

        t_y_xzz_yz_yy[i] = t_0_xyzz_yz_yy[i] - rab_y[i] * t_0_xzz_yz_yy[i];

        t_y_xzz_yz_xz[i] = t_0_xyzz_yz_xz[i] - rab_y[i] * t_0_xzz_yz_xz[i];

        t_y_xzz_yz_xy[i] = t_0_xyzz_yz_xy[i] - rab_y[i] * t_0_xzz_yz_xy[i];

        t_y_xzz_yz_xx[i] = t_0_xyzz_yz_xx[i] - rab_y[i] * t_0_xzz_yz_xx[i];

        t_y_xzz_yy_zz[i] = t_0_xyzz_yy_zz[i] - rab_y[i] * t_0_xzz_yy_zz[i];

        t_y_xzz_yy_yz[i] = t_0_xyzz_yy_yz[i] - rab_y[i] * t_0_xzz_yy_yz[i];

        t_y_xzz_yy_yy[i] = t_0_xyzz_yy_yy[i] - rab_y[i] * t_0_xzz_yy_yy[i];

        t_y_xzz_yy_xz[i] = t_0_xyzz_yy_xz[i] - rab_y[i] * t_0_xzz_yy_xz[i];

        t_y_xzz_yy_xy[i] = t_0_xyzz_yy_xy[i] - rab_y[i] * t_0_xzz_yy_xy[i];

        t_y_xzz_yy_xx[i] = t_0_xyzz_yy_xx[i] - rab_y[i] * t_0_xzz_yy_xx[i];

        t_y_xzz_xz_zz[i] = t_0_xyzz_xz_zz[i] - rab_y[i] * t_0_xzz_xz_zz[i];

        t_y_xzz_xz_yz[i] = t_0_xyzz_xz_yz[i] - rab_y[i] * t_0_xzz_xz_yz[i];

        t_y_xzz_xz_yy[i] = t_0_xyzz_xz_yy[i] - rab_y[i] * t_0_xzz_xz_yy[i];

        t_y_xzz_xz_xz[i] = t_0_xyzz_xz_xz[i] - rab_y[i] * t_0_xzz_xz_xz[i];

        t_y_xzz_xz_xy[i] = t_0_xyzz_xz_xy[i] - rab_y[i] * t_0_xzz_xz_xy[i];

        t_y_xzz_xz_xx[i] = t_0_xyzz_xz_xx[i] - rab_y[i] * t_0_xzz_xz_xx[i];

        t_y_xzz_xy_zz[i] = t_0_xyzz_xy_zz[i] - rab_y[i] * t_0_xzz_xy_zz[i];

        t_y_xzz_xy_yz[i] = t_0_xyzz_xy_yz[i] - rab_y[i] * t_0_xzz_xy_yz[i];

        t_y_xzz_xy_yy[i] = t_0_xyzz_xy_yy[i] - rab_y[i] * t_0_xzz_xy_yy[i];

        t_y_xzz_xy_xz[i] = t_0_xyzz_xy_xz[i] - rab_y[i] * t_0_xzz_xy_xz[i];

        t_y_xzz_xy_xy[i] = t_0_xyzz_xy_xy[i] - rab_y[i] * t_0_xzz_xy_xy[i];

        t_y_xzz_xy_xx[i] = t_0_xyzz_xy_xx[i] - rab_y[i] * t_0_xzz_xy_xx[i];

        t_y_xzz_xx_zz[i] = t_0_xyzz_xx_zz[i] - rab_y[i] * t_0_xzz_xx_zz[i];

        t_y_xzz_xx_yz[i] = t_0_xyzz_xx_yz[i] - rab_y[i] * t_0_xzz_xx_yz[i];

        t_y_xzz_xx_yy[i] = t_0_xyzz_xx_yy[i] - rab_y[i] * t_0_xzz_xx_yy[i];

        t_y_xzz_xx_xz[i] = t_0_xyzz_xx_xz[i] - rab_y[i] * t_0_xzz_xx_xz[i];

        t_y_xzz_xx_xy[i] = t_0_xyzz_xx_xy[i] - rab_y[i] * t_0_xzz_xx_xy[i];

        t_y_xzz_xx_xx[i] = t_0_xyzz_xx_xx[i] - rab_y[i] * t_0_xzz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyyz_xx_xx, t_0_xyyz_xx_xy, t_0_xyyz_xx_xz, t_0_xyyz_xx_yy,\
                           t_0_xyyz_xx_yz, t_0_xyyz_xx_zz, t_0_xyyz_xy_xx, t_0_xyyz_xy_xy,\
                           t_0_xyyz_xy_xz, t_0_xyyz_xy_yy, t_0_xyyz_xy_yz, t_0_xyyz_xy_zz,\
                           t_0_xyyz_xz_xx, t_0_xyyz_xz_xy, t_0_xyyz_xz_xz, t_0_xyyz_xz_yy,\
                           t_0_xyyz_xz_yz, t_0_xyyz_xz_zz, t_0_xyyz_yy_xx, t_0_xyyz_yy_xy,\
                           t_0_xyyz_yy_xz, t_0_xyyz_yy_yy, t_0_xyyz_yy_yz, t_0_xyyz_yy_zz,\
                           t_0_xyyz_yz_xx, t_0_xyyz_yz_xy, t_0_xyyz_yz_xz, t_0_xyyz_yz_yy,\
                           t_0_xyyz_yz_yz, t_0_xyyz_yz_zz, t_0_xyyz_zz_xx, t_0_xyyz_zz_xy,\
                           t_0_xyyz_zz_xz, t_0_xyyz_zz_yy, t_0_xyyz_zz_yz, t_0_xyyz_zz_zz,\
                           t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz, t_0_xyz_xx_yy,\
                           t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx, t_0_xyz_xy_xy,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz, t_0_xyz_xy_zz,\
                           t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz, t_0_xyz_xz_yy,\
                           t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz,\
                           t_y_xyz_xx_xx, t_y_xyz_xx_xy, t_y_xyz_xx_xz, t_y_xyz_xx_yy,\
                           t_y_xyz_xx_yz, t_y_xyz_xx_zz, t_y_xyz_xy_xx, t_y_xyz_xy_xy,\
                           t_y_xyz_xy_xz, t_y_xyz_xy_yy, t_y_xyz_xy_yz, t_y_xyz_xy_zz,\
                           t_y_xyz_xz_xx, t_y_xyz_xz_xy, t_y_xyz_xz_xz, t_y_xyz_xz_yy,\
                           t_y_xyz_xz_yz, t_y_xyz_xz_zz, t_y_xyz_yy_xx, t_y_xyz_yy_xy,\
                           t_y_xyz_yy_xz, t_y_xyz_yy_yy, t_y_xyz_yy_yz, t_y_xyz_yy_zz,\
                           t_y_xyz_yz_xx, t_y_xyz_yz_xy, t_y_xyz_yz_xz, t_y_xyz_yz_yy,\
                           t_y_xyz_yz_yz, t_y_xyz_yz_zz, t_y_xyz_zz_xx, t_y_xyz_zz_xy,\
                           t_y_xyz_zz_xz, t_y_xyz_zz_yy, t_y_xyz_zz_yz, t_y_xyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyz_zz_zz[i] = t_0_xyyz_zz_zz[i] - rab_y[i] * t_0_xyz_zz_zz[i];

        t_y_xyz_zz_yz[i] = t_0_xyyz_zz_yz[i] - rab_y[i] * t_0_xyz_zz_yz[i];

        t_y_xyz_zz_yy[i] = t_0_xyyz_zz_yy[i] - rab_y[i] * t_0_xyz_zz_yy[i];

        t_y_xyz_zz_xz[i] = t_0_xyyz_zz_xz[i] - rab_y[i] * t_0_xyz_zz_xz[i];

        t_y_xyz_zz_xy[i] = t_0_xyyz_zz_xy[i] - rab_y[i] * t_0_xyz_zz_xy[i];

        t_y_xyz_zz_xx[i] = t_0_xyyz_zz_xx[i] - rab_y[i] * t_0_xyz_zz_xx[i];

        t_y_xyz_yz_zz[i] = t_0_xyyz_yz_zz[i] - rab_y[i] * t_0_xyz_yz_zz[i];

        t_y_xyz_yz_yz[i] = t_0_xyyz_yz_yz[i] - rab_y[i] * t_0_xyz_yz_yz[i];

        t_y_xyz_yz_yy[i] = t_0_xyyz_yz_yy[i] - rab_y[i] * t_0_xyz_yz_yy[i];

        t_y_xyz_yz_xz[i] = t_0_xyyz_yz_xz[i] - rab_y[i] * t_0_xyz_yz_xz[i];

        t_y_xyz_yz_xy[i] = t_0_xyyz_yz_xy[i] - rab_y[i] * t_0_xyz_yz_xy[i];

        t_y_xyz_yz_xx[i] = t_0_xyyz_yz_xx[i] - rab_y[i] * t_0_xyz_yz_xx[i];

        t_y_xyz_yy_zz[i] = t_0_xyyz_yy_zz[i] - rab_y[i] * t_0_xyz_yy_zz[i];

        t_y_xyz_yy_yz[i] = t_0_xyyz_yy_yz[i] - rab_y[i] * t_0_xyz_yy_yz[i];

        t_y_xyz_yy_yy[i] = t_0_xyyz_yy_yy[i] - rab_y[i] * t_0_xyz_yy_yy[i];

        t_y_xyz_yy_xz[i] = t_0_xyyz_yy_xz[i] - rab_y[i] * t_0_xyz_yy_xz[i];

        t_y_xyz_yy_xy[i] = t_0_xyyz_yy_xy[i] - rab_y[i] * t_0_xyz_yy_xy[i];

        t_y_xyz_yy_xx[i] = t_0_xyyz_yy_xx[i] - rab_y[i] * t_0_xyz_yy_xx[i];

        t_y_xyz_xz_zz[i] = t_0_xyyz_xz_zz[i] - rab_y[i] * t_0_xyz_xz_zz[i];

        t_y_xyz_xz_yz[i] = t_0_xyyz_xz_yz[i] - rab_y[i] * t_0_xyz_xz_yz[i];

        t_y_xyz_xz_yy[i] = t_0_xyyz_xz_yy[i] - rab_y[i] * t_0_xyz_xz_yy[i];

        t_y_xyz_xz_xz[i] = t_0_xyyz_xz_xz[i] - rab_y[i] * t_0_xyz_xz_xz[i];

        t_y_xyz_xz_xy[i] = t_0_xyyz_xz_xy[i] - rab_y[i] * t_0_xyz_xz_xy[i];

        t_y_xyz_xz_xx[i] = t_0_xyyz_xz_xx[i] - rab_y[i] * t_0_xyz_xz_xx[i];

        t_y_xyz_xy_zz[i] = t_0_xyyz_xy_zz[i] - rab_y[i] * t_0_xyz_xy_zz[i];

        t_y_xyz_xy_yz[i] = t_0_xyyz_xy_yz[i] - rab_y[i] * t_0_xyz_xy_yz[i];

        t_y_xyz_xy_yy[i] = t_0_xyyz_xy_yy[i] - rab_y[i] * t_0_xyz_xy_yy[i];

        t_y_xyz_xy_xz[i] = t_0_xyyz_xy_xz[i] - rab_y[i] * t_0_xyz_xy_xz[i];

        t_y_xyz_xy_xy[i] = t_0_xyyz_xy_xy[i] - rab_y[i] * t_0_xyz_xy_xy[i];

        t_y_xyz_xy_xx[i] = t_0_xyyz_xy_xx[i] - rab_y[i] * t_0_xyz_xy_xx[i];

        t_y_xyz_xx_zz[i] = t_0_xyyz_xx_zz[i] - rab_y[i] * t_0_xyz_xx_zz[i];

        t_y_xyz_xx_yz[i] = t_0_xyyz_xx_yz[i] - rab_y[i] * t_0_xyz_xx_yz[i];

        t_y_xyz_xx_yy[i] = t_0_xyyz_xx_yy[i] - rab_y[i] * t_0_xyz_xx_yy[i];

        t_y_xyz_xx_xz[i] = t_0_xyyz_xx_xz[i] - rab_y[i] * t_0_xyz_xx_xz[i];

        t_y_xyz_xx_xy[i] = t_0_xyyz_xx_xy[i] - rab_y[i] * t_0_xyz_xx_xy[i];

        t_y_xyz_xx_xx[i] = t_0_xyyz_xx_xx[i] - rab_y[i] * t_0_xyz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyy_xx_xx, t_0_xyy_xx_xy, t_0_xyy_xx_xz, t_0_xyy_xx_yy,\
                           t_0_xyy_xx_yz, t_0_xyy_xx_zz, t_0_xyy_xy_xx, t_0_xyy_xy_xy,\
                           t_0_xyy_xy_xz, t_0_xyy_xy_yy, t_0_xyy_xy_yz, t_0_xyy_xy_zz,\
                           t_0_xyy_xz_xx, t_0_xyy_xz_xy, t_0_xyy_xz_xz, t_0_xyy_xz_yy,\
                           t_0_xyy_xz_yz, t_0_xyy_xz_zz, t_0_xyy_yy_xx, t_0_xyy_yy_xy,\
                           t_0_xyy_yy_xz, t_0_xyy_yy_yy, t_0_xyy_yy_yz, t_0_xyy_yy_zz,\
                           t_0_xyy_yz_xx, t_0_xyy_yz_xy, t_0_xyy_yz_xz, t_0_xyy_yz_yy,\
                           t_0_xyy_yz_yz, t_0_xyy_yz_zz, t_0_xyy_zz_xx, t_0_xyy_zz_xy,\
                           t_0_xyy_zz_xz, t_0_xyy_zz_yy, t_0_xyy_zz_yz, t_0_xyy_zz_zz,\
                           t_0_xyyy_xx_xx, t_0_xyyy_xx_xy, t_0_xyyy_xx_xz, t_0_xyyy_xx_yy,\
                           t_0_xyyy_xx_yz, t_0_xyyy_xx_zz, t_0_xyyy_xy_xx, t_0_xyyy_xy_xy,\
                           t_0_xyyy_xy_xz, t_0_xyyy_xy_yy, t_0_xyyy_xy_yz, t_0_xyyy_xy_zz,\
                           t_0_xyyy_xz_xx, t_0_xyyy_xz_xy, t_0_xyyy_xz_xz, t_0_xyyy_xz_yy,\
                           t_0_xyyy_xz_yz, t_0_xyyy_xz_zz, t_0_xyyy_yy_xx, t_0_xyyy_yy_xy,\
                           t_0_xyyy_yy_xz, t_0_xyyy_yy_yy, t_0_xyyy_yy_yz, t_0_xyyy_yy_zz,\
                           t_0_xyyy_yz_xx, t_0_xyyy_yz_xy, t_0_xyyy_yz_xz, t_0_xyyy_yz_yy,\
                           t_0_xyyy_yz_yz, t_0_xyyy_yz_zz, t_0_xyyy_zz_xx, t_0_xyyy_zz_xy,\
                           t_0_xyyy_zz_xz, t_0_xyyy_zz_yy, t_0_xyyy_zz_yz, t_0_xyyy_zz_zz,\
                           t_y_xyy_xx_xx, t_y_xyy_xx_xy, t_y_xyy_xx_xz, t_y_xyy_xx_yy,\
                           t_y_xyy_xx_yz, t_y_xyy_xx_zz, t_y_xyy_xy_xx, t_y_xyy_xy_xy,\
                           t_y_xyy_xy_xz, t_y_xyy_xy_yy, t_y_xyy_xy_yz, t_y_xyy_xy_zz,\
                           t_y_xyy_xz_xx, t_y_xyy_xz_xy, t_y_xyy_xz_xz, t_y_xyy_xz_yy,\
                           t_y_xyy_xz_yz, t_y_xyy_xz_zz, t_y_xyy_yy_xx, t_y_xyy_yy_xy,\
                           t_y_xyy_yy_xz, t_y_xyy_yy_yy, t_y_xyy_yy_yz, t_y_xyy_yy_zz,\
                           t_y_xyy_yz_xx, t_y_xyy_yz_xy, t_y_xyy_yz_xz, t_y_xyy_yz_yy,\
                           t_y_xyy_yz_yz, t_y_xyy_yz_zz, t_y_xyy_zz_xx, t_y_xyy_zz_xy,\
                           t_y_xyy_zz_xz, t_y_xyy_zz_yy, t_y_xyy_zz_yz, t_y_xyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyy_zz_zz[i] = t_0_xyyy_zz_zz[i] - rab_y[i] * t_0_xyy_zz_zz[i];

        t_y_xyy_zz_yz[i] = t_0_xyyy_zz_yz[i] - rab_y[i] * t_0_xyy_zz_yz[i];

        t_y_xyy_zz_yy[i] = t_0_xyyy_zz_yy[i] - rab_y[i] * t_0_xyy_zz_yy[i];

        t_y_xyy_zz_xz[i] = t_0_xyyy_zz_xz[i] - rab_y[i] * t_0_xyy_zz_xz[i];

        t_y_xyy_zz_xy[i] = t_0_xyyy_zz_xy[i] - rab_y[i] * t_0_xyy_zz_xy[i];

        t_y_xyy_zz_xx[i] = t_0_xyyy_zz_xx[i] - rab_y[i] * t_0_xyy_zz_xx[i];

        t_y_xyy_yz_zz[i] = t_0_xyyy_yz_zz[i] - rab_y[i] * t_0_xyy_yz_zz[i];

        t_y_xyy_yz_yz[i] = t_0_xyyy_yz_yz[i] - rab_y[i] * t_0_xyy_yz_yz[i];

        t_y_xyy_yz_yy[i] = t_0_xyyy_yz_yy[i] - rab_y[i] * t_0_xyy_yz_yy[i];

        t_y_xyy_yz_xz[i] = t_0_xyyy_yz_xz[i] - rab_y[i] * t_0_xyy_yz_xz[i];

        t_y_xyy_yz_xy[i] = t_0_xyyy_yz_xy[i] - rab_y[i] * t_0_xyy_yz_xy[i];

        t_y_xyy_yz_xx[i] = t_0_xyyy_yz_xx[i] - rab_y[i] * t_0_xyy_yz_xx[i];

        t_y_xyy_yy_zz[i] = t_0_xyyy_yy_zz[i] - rab_y[i] * t_0_xyy_yy_zz[i];

        t_y_xyy_yy_yz[i] = t_0_xyyy_yy_yz[i] - rab_y[i] * t_0_xyy_yy_yz[i];

        t_y_xyy_yy_yy[i] = t_0_xyyy_yy_yy[i] - rab_y[i] * t_0_xyy_yy_yy[i];

        t_y_xyy_yy_xz[i] = t_0_xyyy_yy_xz[i] - rab_y[i] * t_0_xyy_yy_xz[i];

        t_y_xyy_yy_xy[i] = t_0_xyyy_yy_xy[i] - rab_y[i] * t_0_xyy_yy_xy[i];

        t_y_xyy_yy_xx[i] = t_0_xyyy_yy_xx[i] - rab_y[i] * t_0_xyy_yy_xx[i];

        t_y_xyy_xz_zz[i] = t_0_xyyy_xz_zz[i] - rab_y[i] * t_0_xyy_xz_zz[i];

        t_y_xyy_xz_yz[i] = t_0_xyyy_xz_yz[i] - rab_y[i] * t_0_xyy_xz_yz[i];

        t_y_xyy_xz_yy[i] = t_0_xyyy_xz_yy[i] - rab_y[i] * t_0_xyy_xz_yy[i];

        t_y_xyy_xz_xz[i] = t_0_xyyy_xz_xz[i] - rab_y[i] * t_0_xyy_xz_xz[i];

        t_y_xyy_xz_xy[i] = t_0_xyyy_xz_xy[i] - rab_y[i] * t_0_xyy_xz_xy[i];

        t_y_xyy_xz_xx[i] = t_0_xyyy_xz_xx[i] - rab_y[i] * t_0_xyy_xz_xx[i];

        t_y_xyy_xy_zz[i] = t_0_xyyy_xy_zz[i] - rab_y[i] * t_0_xyy_xy_zz[i];

        t_y_xyy_xy_yz[i] = t_0_xyyy_xy_yz[i] - rab_y[i] * t_0_xyy_xy_yz[i];

        t_y_xyy_xy_yy[i] = t_0_xyyy_xy_yy[i] - rab_y[i] * t_0_xyy_xy_yy[i];

        t_y_xyy_xy_xz[i] = t_0_xyyy_xy_xz[i] - rab_y[i] * t_0_xyy_xy_xz[i];

        t_y_xyy_xy_xy[i] = t_0_xyyy_xy_xy[i] - rab_y[i] * t_0_xyy_xy_xy[i];

        t_y_xyy_xy_xx[i] = t_0_xyyy_xy_xx[i] - rab_y[i] * t_0_xyy_xy_xx[i];

        t_y_xyy_xx_zz[i] = t_0_xyyy_xx_zz[i] - rab_y[i] * t_0_xyy_xx_zz[i];

        t_y_xyy_xx_yz[i] = t_0_xyyy_xx_yz[i] - rab_y[i] * t_0_xyy_xx_yz[i];

        t_y_xyy_xx_yy[i] = t_0_xyyy_xx_yy[i] - rab_y[i] * t_0_xyy_xx_yy[i];

        t_y_xyy_xx_xz[i] = t_0_xyyy_xx_xz[i] - rab_y[i] * t_0_xyy_xx_xz[i];

        t_y_xyy_xx_xy[i] = t_0_xyyy_xx_xy[i] - rab_y[i] * t_0_xyy_xx_xy[i];

        t_y_xyy_xx_xx[i] = t_0_xyyy_xx_xx[i] - rab_y[i] * t_0_xyy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xxyz_xx_xx, t_0_xxyz_xx_xy, t_0_xxyz_xx_xz, t_0_xxyz_xx_yy,\
                           t_0_xxyz_xx_yz, t_0_xxyz_xx_zz, t_0_xxyz_xy_xx, t_0_xxyz_xy_xy,\
                           t_0_xxyz_xy_xz, t_0_xxyz_xy_yy, t_0_xxyz_xy_yz, t_0_xxyz_xy_zz,\
                           t_0_xxyz_xz_xx, t_0_xxyz_xz_xy, t_0_xxyz_xz_xz, t_0_xxyz_xz_yy,\
                           t_0_xxyz_xz_yz, t_0_xxyz_xz_zz, t_0_xxyz_yy_xx, t_0_xxyz_yy_xy,\
                           t_0_xxyz_yy_xz, t_0_xxyz_yy_yy, t_0_xxyz_yy_yz, t_0_xxyz_yy_zz,\
                           t_0_xxyz_yz_xx, t_0_xxyz_yz_xy, t_0_xxyz_yz_xz, t_0_xxyz_yz_yy,\
                           t_0_xxyz_yz_yz, t_0_xxyz_yz_zz, t_0_xxyz_zz_xx, t_0_xxyz_zz_xy,\
                           t_0_xxyz_zz_xz, t_0_xxyz_zz_yy, t_0_xxyz_zz_yz, t_0_xxyz_zz_zz,\
                           t_0_xxz_xx_xx, t_0_xxz_xx_xy, t_0_xxz_xx_xz, t_0_xxz_xx_yy,\
                           t_0_xxz_xx_yz, t_0_xxz_xx_zz, t_0_xxz_xy_xx, t_0_xxz_xy_xy,\
                           t_0_xxz_xy_xz, t_0_xxz_xy_yy, t_0_xxz_xy_yz, t_0_xxz_xy_zz,\
                           t_0_xxz_xz_xx, t_0_xxz_xz_xy, t_0_xxz_xz_xz, t_0_xxz_xz_yy,\
                           t_0_xxz_xz_yz, t_0_xxz_xz_zz, t_0_xxz_yy_xx, t_0_xxz_yy_xy,\
                           t_0_xxz_yy_xz, t_0_xxz_yy_yy, t_0_xxz_yy_yz, t_0_xxz_yy_zz,\
                           t_0_xxz_yz_xx, t_0_xxz_yz_xy, t_0_xxz_yz_xz, t_0_xxz_yz_yy,\
                           t_0_xxz_yz_yz, t_0_xxz_yz_zz, t_0_xxz_zz_xx, t_0_xxz_zz_xy,\
                           t_0_xxz_zz_xz, t_0_xxz_zz_yy, t_0_xxz_zz_yz, t_0_xxz_zz_zz,\
                           t_y_xxz_xx_xx, t_y_xxz_xx_xy, t_y_xxz_xx_xz, t_y_xxz_xx_yy,\
                           t_y_xxz_xx_yz, t_y_xxz_xx_zz, t_y_xxz_xy_xx, t_y_xxz_xy_xy,\
                           t_y_xxz_xy_xz, t_y_xxz_xy_yy, t_y_xxz_xy_yz, t_y_xxz_xy_zz,\
                           t_y_xxz_xz_xx, t_y_xxz_xz_xy, t_y_xxz_xz_xz, t_y_xxz_xz_yy,\
                           t_y_xxz_xz_yz, t_y_xxz_xz_zz, t_y_xxz_yy_xx, t_y_xxz_yy_xy,\
                           t_y_xxz_yy_xz, t_y_xxz_yy_yy, t_y_xxz_yy_yz, t_y_xxz_yy_zz,\
                           t_y_xxz_yz_xx, t_y_xxz_yz_xy, t_y_xxz_yz_xz, t_y_xxz_yz_yy,\
                           t_y_xxz_yz_yz, t_y_xxz_yz_zz, t_y_xxz_zz_xx, t_y_xxz_zz_xy,\
                           t_y_xxz_zz_xz, t_y_xxz_zz_yy, t_y_xxz_zz_yz, t_y_xxz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xxz_zz_zz[i] = t_0_xxyz_zz_zz[i] - rab_y[i] * t_0_xxz_zz_zz[i];

        t_y_xxz_zz_yz[i] = t_0_xxyz_zz_yz[i] - rab_y[i] * t_0_xxz_zz_yz[i];

        t_y_xxz_zz_yy[i] = t_0_xxyz_zz_yy[i] - rab_y[i] * t_0_xxz_zz_yy[i];

        t_y_xxz_zz_xz[i] = t_0_xxyz_zz_xz[i] - rab_y[i] * t_0_xxz_zz_xz[i];

        t_y_xxz_zz_xy[i] = t_0_xxyz_zz_xy[i] - rab_y[i] * t_0_xxz_zz_xy[i];

        t_y_xxz_zz_xx[i] = t_0_xxyz_zz_xx[i] - rab_y[i] * t_0_xxz_zz_xx[i];

        t_y_xxz_yz_zz[i] = t_0_xxyz_yz_zz[i] - rab_y[i] * t_0_xxz_yz_zz[i];

        t_y_xxz_yz_yz[i] = t_0_xxyz_yz_yz[i] - rab_y[i] * t_0_xxz_yz_yz[i];

        t_y_xxz_yz_yy[i] = t_0_xxyz_yz_yy[i] - rab_y[i] * t_0_xxz_yz_yy[i];

        t_y_xxz_yz_xz[i] = t_0_xxyz_yz_xz[i] - rab_y[i] * t_0_xxz_yz_xz[i];

        t_y_xxz_yz_xy[i] = t_0_xxyz_yz_xy[i] - rab_y[i] * t_0_xxz_yz_xy[i];

        t_y_xxz_yz_xx[i] = t_0_xxyz_yz_xx[i] - rab_y[i] * t_0_xxz_yz_xx[i];

        t_y_xxz_yy_zz[i] = t_0_xxyz_yy_zz[i] - rab_y[i] * t_0_xxz_yy_zz[i];

        t_y_xxz_yy_yz[i] = t_0_xxyz_yy_yz[i] - rab_y[i] * t_0_xxz_yy_yz[i];

        t_y_xxz_yy_yy[i] = t_0_xxyz_yy_yy[i] - rab_y[i] * t_0_xxz_yy_yy[i];

        t_y_xxz_yy_xz[i] = t_0_xxyz_yy_xz[i] - rab_y[i] * t_0_xxz_yy_xz[i];

        t_y_xxz_yy_xy[i] = t_0_xxyz_yy_xy[i] - rab_y[i] * t_0_xxz_yy_xy[i];

        t_y_xxz_yy_xx[i] = t_0_xxyz_yy_xx[i] - rab_y[i] * t_0_xxz_yy_xx[i];

        t_y_xxz_xz_zz[i] = t_0_xxyz_xz_zz[i] - rab_y[i] * t_0_xxz_xz_zz[i];

        t_y_xxz_xz_yz[i] = t_0_xxyz_xz_yz[i] - rab_y[i] * t_0_xxz_xz_yz[i];

        t_y_xxz_xz_yy[i] = t_0_xxyz_xz_yy[i] - rab_y[i] * t_0_xxz_xz_yy[i];

        t_y_xxz_xz_xz[i] = t_0_xxyz_xz_xz[i] - rab_y[i] * t_0_xxz_xz_xz[i];

        t_y_xxz_xz_xy[i] = t_0_xxyz_xz_xy[i] - rab_y[i] * t_0_xxz_xz_xy[i];

        t_y_xxz_xz_xx[i] = t_0_xxyz_xz_xx[i] - rab_y[i] * t_0_xxz_xz_xx[i];

        t_y_xxz_xy_zz[i] = t_0_xxyz_xy_zz[i] - rab_y[i] * t_0_xxz_xy_zz[i];

        t_y_xxz_xy_yz[i] = t_0_xxyz_xy_yz[i] - rab_y[i] * t_0_xxz_xy_yz[i];

        t_y_xxz_xy_yy[i] = t_0_xxyz_xy_yy[i] - rab_y[i] * t_0_xxz_xy_yy[i];

        t_y_xxz_xy_xz[i] = t_0_xxyz_xy_xz[i] - rab_y[i] * t_0_xxz_xy_xz[i];

        t_y_xxz_xy_xy[i] = t_0_xxyz_xy_xy[i] - rab_y[i] * t_0_xxz_xy_xy[i];

        t_y_xxz_xy_xx[i] = t_0_xxyz_xy_xx[i] - rab_y[i] * t_0_xxz_xy_xx[i];

        t_y_xxz_xx_zz[i] = t_0_xxyz_xx_zz[i] - rab_y[i] * t_0_xxz_xx_zz[i];

        t_y_xxz_xx_yz[i] = t_0_xxyz_xx_yz[i] - rab_y[i] * t_0_xxz_xx_yz[i];

        t_y_xxz_xx_yy[i] = t_0_xxyz_xx_yy[i] - rab_y[i] * t_0_xxz_xx_yy[i];

        t_y_xxz_xx_xz[i] = t_0_xxyz_xx_xz[i] - rab_y[i] * t_0_xxz_xx_xz[i];

        t_y_xxz_xx_xy[i] = t_0_xxyz_xx_xy[i] - rab_y[i] * t_0_xxz_xx_xy[i];

        t_y_xxz_xx_xx[i] = t_0_xxyz_xx_xx[i] - rab_y[i] * t_0_xxz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xxy_xx_xx, t_0_xxy_xx_xy, t_0_xxy_xx_xz, t_0_xxy_xx_yy,\
                           t_0_xxy_xx_yz, t_0_xxy_xx_zz, t_0_xxy_xy_xx, t_0_xxy_xy_xy,\
                           t_0_xxy_xy_xz, t_0_xxy_xy_yy, t_0_xxy_xy_yz, t_0_xxy_xy_zz,\
                           t_0_xxy_xz_xx, t_0_xxy_xz_xy, t_0_xxy_xz_xz, t_0_xxy_xz_yy,\
                           t_0_xxy_xz_yz, t_0_xxy_xz_zz, t_0_xxy_yy_xx, t_0_xxy_yy_xy,\
                           t_0_xxy_yy_xz, t_0_xxy_yy_yy, t_0_xxy_yy_yz, t_0_xxy_yy_zz,\
                           t_0_xxy_yz_xx, t_0_xxy_yz_xy, t_0_xxy_yz_xz, t_0_xxy_yz_yy,\
                           t_0_xxy_yz_yz, t_0_xxy_yz_zz, t_0_xxy_zz_xx, t_0_xxy_zz_xy,\
                           t_0_xxy_zz_xz, t_0_xxy_zz_yy, t_0_xxy_zz_yz, t_0_xxy_zz_zz,\
                           t_0_xxyy_xx_xx, t_0_xxyy_xx_xy, t_0_xxyy_xx_xz, t_0_xxyy_xx_yy,\
                           t_0_xxyy_xx_yz, t_0_xxyy_xx_zz, t_0_xxyy_xy_xx, t_0_xxyy_xy_xy,\
                           t_0_xxyy_xy_xz, t_0_xxyy_xy_yy, t_0_xxyy_xy_yz, t_0_xxyy_xy_zz,\
                           t_0_xxyy_xz_xx, t_0_xxyy_xz_xy, t_0_xxyy_xz_xz, t_0_xxyy_xz_yy,\
                           t_0_xxyy_xz_yz, t_0_xxyy_xz_zz, t_0_xxyy_yy_xx, t_0_xxyy_yy_xy,\
                           t_0_xxyy_yy_xz, t_0_xxyy_yy_yy, t_0_xxyy_yy_yz, t_0_xxyy_yy_zz,\
                           t_0_xxyy_yz_xx, t_0_xxyy_yz_xy, t_0_xxyy_yz_xz, t_0_xxyy_yz_yy,\
                           t_0_xxyy_yz_yz, t_0_xxyy_yz_zz, t_0_xxyy_zz_xx, t_0_xxyy_zz_xy,\
                           t_0_xxyy_zz_xz, t_0_xxyy_zz_yy, t_0_xxyy_zz_yz, t_0_xxyy_zz_zz,\
                           t_y_xxy_xx_xx, t_y_xxy_xx_xy, t_y_xxy_xx_xz, t_y_xxy_xx_yy,\
                           t_y_xxy_xx_yz, t_y_xxy_xx_zz, t_y_xxy_xy_xx, t_y_xxy_xy_xy,\
                           t_y_xxy_xy_xz, t_y_xxy_xy_yy, t_y_xxy_xy_yz, t_y_xxy_xy_zz,\
                           t_y_xxy_xz_xx, t_y_xxy_xz_xy, t_y_xxy_xz_xz, t_y_xxy_xz_yy,\
                           t_y_xxy_xz_yz, t_y_xxy_xz_zz, t_y_xxy_yy_xx, t_y_xxy_yy_xy,\
                           t_y_xxy_yy_xz, t_y_xxy_yy_yy, t_y_xxy_yy_yz, t_y_xxy_yy_zz,\
                           t_y_xxy_yz_xx, t_y_xxy_yz_xy, t_y_xxy_yz_xz, t_y_xxy_yz_yy,\
                           t_y_xxy_yz_yz, t_y_xxy_yz_zz, t_y_xxy_zz_xx, t_y_xxy_zz_xy,\
                           t_y_xxy_zz_xz, t_y_xxy_zz_yy, t_y_xxy_zz_yz, t_y_xxy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xxy_zz_zz[i] = t_0_xxyy_zz_zz[i] - rab_y[i] * t_0_xxy_zz_zz[i];

        t_y_xxy_zz_yz[i] = t_0_xxyy_zz_yz[i] - rab_y[i] * t_0_xxy_zz_yz[i];

        t_y_xxy_zz_yy[i] = t_0_xxyy_zz_yy[i] - rab_y[i] * t_0_xxy_zz_yy[i];

        t_y_xxy_zz_xz[i] = t_0_xxyy_zz_xz[i] - rab_y[i] * t_0_xxy_zz_xz[i];

        t_y_xxy_zz_xy[i] = t_0_xxyy_zz_xy[i] - rab_y[i] * t_0_xxy_zz_xy[i];

        t_y_xxy_zz_xx[i] = t_0_xxyy_zz_xx[i] - rab_y[i] * t_0_xxy_zz_xx[i];

        t_y_xxy_yz_zz[i] = t_0_xxyy_yz_zz[i] - rab_y[i] * t_0_xxy_yz_zz[i];

        t_y_xxy_yz_yz[i] = t_0_xxyy_yz_yz[i] - rab_y[i] * t_0_xxy_yz_yz[i];

        t_y_xxy_yz_yy[i] = t_0_xxyy_yz_yy[i] - rab_y[i] * t_0_xxy_yz_yy[i];

        t_y_xxy_yz_xz[i] = t_0_xxyy_yz_xz[i] - rab_y[i] * t_0_xxy_yz_xz[i];

        t_y_xxy_yz_xy[i] = t_0_xxyy_yz_xy[i] - rab_y[i] * t_0_xxy_yz_xy[i];

        t_y_xxy_yz_xx[i] = t_0_xxyy_yz_xx[i] - rab_y[i] * t_0_xxy_yz_xx[i];

        t_y_xxy_yy_zz[i] = t_0_xxyy_yy_zz[i] - rab_y[i] * t_0_xxy_yy_zz[i];

        t_y_xxy_yy_yz[i] = t_0_xxyy_yy_yz[i] - rab_y[i] * t_0_xxy_yy_yz[i];

        t_y_xxy_yy_yy[i] = t_0_xxyy_yy_yy[i] - rab_y[i] * t_0_xxy_yy_yy[i];

        t_y_xxy_yy_xz[i] = t_0_xxyy_yy_xz[i] - rab_y[i] * t_0_xxy_yy_xz[i];

        t_y_xxy_yy_xy[i] = t_0_xxyy_yy_xy[i] - rab_y[i] * t_0_xxy_yy_xy[i];

        t_y_xxy_yy_xx[i] = t_0_xxyy_yy_xx[i] - rab_y[i] * t_0_xxy_yy_xx[i];

        t_y_xxy_xz_zz[i] = t_0_xxyy_xz_zz[i] - rab_y[i] * t_0_xxy_xz_zz[i];

        t_y_xxy_xz_yz[i] = t_0_xxyy_xz_yz[i] - rab_y[i] * t_0_xxy_xz_yz[i];

        t_y_xxy_xz_yy[i] = t_0_xxyy_xz_yy[i] - rab_y[i] * t_0_xxy_xz_yy[i];

        t_y_xxy_xz_xz[i] = t_0_xxyy_xz_xz[i] - rab_y[i] * t_0_xxy_xz_xz[i];

        t_y_xxy_xz_xy[i] = t_0_xxyy_xz_xy[i] - rab_y[i] * t_0_xxy_xz_xy[i];

        t_y_xxy_xz_xx[i] = t_0_xxyy_xz_xx[i] - rab_y[i] * t_0_xxy_xz_xx[i];

        t_y_xxy_xy_zz[i] = t_0_xxyy_xy_zz[i] - rab_y[i] * t_0_xxy_xy_zz[i];

        t_y_xxy_xy_yz[i] = t_0_xxyy_xy_yz[i] - rab_y[i] * t_0_xxy_xy_yz[i];

        t_y_xxy_xy_yy[i] = t_0_xxyy_xy_yy[i] - rab_y[i] * t_0_xxy_xy_yy[i];

        t_y_xxy_xy_xz[i] = t_0_xxyy_xy_xz[i] - rab_y[i] * t_0_xxy_xy_xz[i];

        t_y_xxy_xy_xy[i] = t_0_xxyy_xy_xy[i] - rab_y[i] * t_0_xxy_xy_xy[i];

        t_y_xxy_xy_xx[i] = t_0_xxyy_xy_xx[i] - rab_y[i] * t_0_xxy_xy_xx[i];

        t_y_xxy_xx_zz[i] = t_0_xxyy_xx_zz[i] - rab_y[i] * t_0_xxy_xx_zz[i];

        t_y_xxy_xx_yz[i] = t_0_xxyy_xx_yz[i] - rab_y[i] * t_0_xxy_xx_yz[i];

        t_y_xxy_xx_yy[i] = t_0_xxyy_xx_yy[i] - rab_y[i] * t_0_xxy_xx_yy[i];

        t_y_xxy_xx_xz[i] = t_0_xxyy_xx_xz[i] - rab_y[i] * t_0_xxy_xx_xz[i];

        t_y_xxy_xx_xy[i] = t_0_xxyy_xx_xy[i] - rab_y[i] * t_0_xxy_xx_xy[i];

        t_y_xxy_xx_xx[i] = t_0_xxyy_xx_xx[i] - rab_y[i] * t_0_xxy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xzzz_xx_xx, t_0_xzzz_xx_xy, t_0_xzzz_xx_xz, t_0_xzzz_xx_yy,\
                           t_0_xzzz_xx_yz, t_0_xzzz_xx_zz, t_0_xzzz_xy_xx, t_0_xzzz_xy_xy,\
                           t_0_xzzz_xy_xz, t_0_xzzz_xy_yy, t_0_xzzz_xy_yz, t_0_xzzz_xy_zz,\
                           t_0_xzzz_xz_xx, t_0_xzzz_xz_xy, t_0_xzzz_xz_xz, t_0_xzzz_xz_yy,\
                           t_0_xzzz_xz_yz, t_0_xzzz_xz_zz, t_0_xzzz_yy_xx, t_0_xzzz_yy_xy,\
                           t_0_xzzz_yy_xz, t_0_xzzz_yy_yy, t_0_xzzz_yy_yz, t_0_xzzz_yy_zz,\
                           t_0_xzzz_yz_xx, t_0_xzzz_yz_xy, t_0_xzzz_yz_xz, t_0_xzzz_yz_yy,\
                           t_0_xzzz_yz_yz, t_0_xzzz_yz_zz, t_0_xzzz_zz_xx, t_0_xzzz_zz_xy,\
                           t_0_xzzz_zz_xz, t_0_xzzz_zz_yy, t_0_xzzz_zz_yz, t_0_xzzz_zz_zz,\
                           t_0_zzz_xx_xx, t_0_zzz_xx_xy, t_0_zzz_xx_xz, t_0_zzz_xx_yy,\
                           t_0_zzz_xx_yz, t_0_zzz_xx_zz, t_0_zzz_xy_xx, t_0_zzz_xy_xy,\
                           t_0_zzz_xy_xz, t_0_zzz_xy_yy, t_0_zzz_xy_yz, t_0_zzz_xy_zz,\
                           t_0_zzz_xz_xx, t_0_zzz_xz_xy, t_0_zzz_xz_xz, t_0_zzz_xz_yy,\
                           t_0_zzz_xz_yz, t_0_zzz_xz_zz, t_0_zzz_yy_xx, t_0_zzz_yy_xy,\
                           t_0_zzz_yy_xz, t_0_zzz_yy_yy, t_0_zzz_yy_yz, t_0_zzz_yy_zz,\
                           t_0_zzz_yz_xx, t_0_zzz_yz_xy, t_0_zzz_yz_xz, t_0_zzz_yz_yy,\
                           t_0_zzz_yz_yz, t_0_zzz_yz_zz, t_0_zzz_zz_xx, t_0_zzz_zz_xy,\
                           t_0_zzz_zz_xz, t_0_zzz_zz_yy, t_0_zzz_zz_yz, t_0_zzz_zz_zz,\
                           t_x_zzz_xx_xx, t_x_zzz_xx_xy, t_x_zzz_xx_xz, t_x_zzz_xx_yy,\
                           t_x_zzz_xx_yz, t_x_zzz_xx_zz, t_x_zzz_xy_xx, t_x_zzz_xy_xy,\
                           t_x_zzz_xy_xz, t_x_zzz_xy_yy, t_x_zzz_xy_yz, t_x_zzz_xy_zz,\
                           t_x_zzz_xz_xx, t_x_zzz_xz_xy, t_x_zzz_xz_xz, t_x_zzz_xz_yy,\
                           t_x_zzz_xz_yz, t_x_zzz_xz_zz, t_x_zzz_yy_xx, t_x_zzz_yy_xy,\
                           t_x_zzz_yy_xz, t_x_zzz_yy_yy, t_x_zzz_yy_yz, t_x_zzz_yy_zz,\
                           t_x_zzz_yz_xx, t_x_zzz_yz_xy, t_x_zzz_yz_xz, t_x_zzz_yz_yy,\
                           t_x_zzz_yz_yz, t_x_zzz_yz_zz, t_x_zzz_zz_xx, t_x_zzz_zz_xy,\
                           t_x_zzz_zz_xz, t_x_zzz_zz_yy, t_x_zzz_zz_yz, t_x_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zzz_zz_zz[i] = t_0_xzzz_zz_zz[i] - rab_x[i] * t_0_zzz_zz_zz[i];

        t_x_zzz_zz_yz[i] = t_0_xzzz_zz_yz[i] - rab_x[i] * t_0_zzz_zz_yz[i];

        t_x_zzz_zz_yy[i] = t_0_xzzz_zz_yy[i] - rab_x[i] * t_0_zzz_zz_yy[i];

        t_x_zzz_zz_xz[i] = t_0_xzzz_zz_xz[i] - rab_x[i] * t_0_zzz_zz_xz[i];

        t_x_zzz_zz_xy[i] = t_0_xzzz_zz_xy[i] - rab_x[i] * t_0_zzz_zz_xy[i];

        t_x_zzz_zz_xx[i] = t_0_xzzz_zz_xx[i] - rab_x[i] * t_0_zzz_zz_xx[i];

        t_x_zzz_yz_zz[i] = t_0_xzzz_yz_zz[i] - rab_x[i] * t_0_zzz_yz_zz[i];

        t_x_zzz_yz_yz[i] = t_0_xzzz_yz_yz[i] - rab_x[i] * t_0_zzz_yz_yz[i];

        t_x_zzz_yz_yy[i] = t_0_xzzz_yz_yy[i] - rab_x[i] * t_0_zzz_yz_yy[i];

        t_x_zzz_yz_xz[i] = t_0_xzzz_yz_xz[i] - rab_x[i] * t_0_zzz_yz_xz[i];

        t_x_zzz_yz_xy[i] = t_0_xzzz_yz_xy[i] - rab_x[i] * t_0_zzz_yz_xy[i];

        t_x_zzz_yz_xx[i] = t_0_xzzz_yz_xx[i] - rab_x[i] * t_0_zzz_yz_xx[i];

        t_x_zzz_yy_zz[i] = t_0_xzzz_yy_zz[i] - rab_x[i] * t_0_zzz_yy_zz[i];

        t_x_zzz_yy_yz[i] = t_0_xzzz_yy_yz[i] - rab_x[i] * t_0_zzz_yy_yz[i];

        t_x_zzz_yy_yy[i] = t_0_xzzz_yy_yy[i] - rab_x[i] * t_0_zzz_yy_yy[i];

        t_x_zzz_yy_xz[i] = t_0_xzzz_yy_xz[i] - rab_x[i] * t_0_zzz_yy_xz[i];

        t_x_zzz_yy_xy[i] = t_0_xzzz_yy_xy[i] - rab_x[i] * t_0_zzz_yy_xy[i];

        t_x_zzz_yy_xx[i] = t_0_xzzz_yy_xx[i] - rab_x[i] * t_0_zzz_yy_xx[i];

        t_x_zzz_xz_zz[i] = t_0_xzzz_xz_zz[i] - rab_x[i] * t_0_zzz_xz_zz[i];

        t_x_zzz_xz_yz[i] = t_0_xzzz_xz_yz[i] - rab_x[i] * t_0_zzz_xz_yz[i];

        t_x_zzz_xz_yy[i] = t_0_xzzz_xz_yy[i] - rab_x[i] * t_0_zzz_xz_yy[i];

        t_x_zzz_xz_xz[i] = t_0_xzzz_xz_xz[i] - rab_x[i] * t_0_zzz_xz_xz[i];

        t_x_zzz_xz_xy[i] = t_0_xzzz_xz_xy[i] - rab_x[i] * t_0_zzz_xz_xy[i];

        t_x_zzz_xz_xx[i] = t_0_xzzz_xz_xx[i] - rab_x[i] * t_0_zzz_xz_xx[i];

        t_x_zzz_xy_zz[i] = t_0_xzzz_xy_zz[i] - rab_x[i] * t_0_zzz_xy_zz[i];

        t_x_zzz_xy_yz[i] = t_0_xzzz_xy_yz[i] - rab_x[i] * t_0_zzz_xy_yz[i];

        t_x_zzz_xy_yy[i] = t_0_xzzz_xy_yy[i] - rab_x[i] * t_0_zzz_xy_yy[i];

        t_x_zzz_xy_xz[i] = t_0_xzzz_xy_xz[i] - rab_x[i] * t_0_zzz_xy_xz[i];

        t_x_zzz_xy_xy[i] = t_0_xzzz_xy_xy[i] - rab_x[i] * t_0_zzz_xy_xy[i];

        t_x_zzz_xy_xx[i] = t_0_xzzz_xy_xx[i] - rab_x[i] * t_0_zzz_xy_xx[i];

        t_x_zzz_xx_zz[i] = t_0_xzzz_xx_zz[i] - rab_x[i] * t_0_zzz_xx_zz[i];

        t_x_zzz_xx_yz[i] = t_0_xzzz_xx_yz[i] - rab_x[i] * t_0_zzz_xx_yz[i];

        t_x_zzz_xx_yy[i] = t_0_xzzz_xx_yy[i] - rab_x[i] * t_0_zzz_xx_yy[i];

        t_x_zzz_xx_xz[i] = t_0_xzzz_xx_xz[i] - rab_x[i] * t_0_zzz_xx_xz[i];

        t_x_zzz_xx_xy[i] = t_0_xzzz_xx_xy[i] - rab_x[i] * t_0_zzz_xx_xy[i];

        t_x_zzz_xx_xx[i] = t_0_xzzz_xx_xx[i] - rab_x[i] * t_0_zzz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyzz_xx_xx, t_0_xyzz_xx_xy, t_0_xyzz_xx_xz, t_0_xyzz_xx_yy,\
                           t_0_xyzz_xx_yz, t_0_xyzz_xx_zz, t_0_xyzz_xy_xx, t_0_xyzz_xy_xy,\
                           t_0_xyzz_xy_xz, t_0_xyzz_xy_yy, t_0_xyzz_xy_yz, t_0_xyzz_xy_zz,\
                           t_0_xyzz_xz_xx, t_0_xyzz_xz_xy, t_0_xyzz_xz_xz, t_0_xyzz_xz_yy,\
                           t_0_xyzz_xz_yz, t_0_xyzz_xz_zz, t_0_xyzz_yy_xx, t_0_xyzz_yy_xy,\
                           t_0_xyzz_yy_xz, t_0_xyzz_yy_yy, t_0_xyzz_yy_yz, t_0_xyzz_yy_zz,\
                           t_0_xyzz_yz_xx, t_0_xyzz_yz_xy, t_0_xyzz_yz_xz, t_0_xyzz_yz_yy,\
                           t_0_xyzz_yz_yz, t_0_xyzz_yz_zz, t_0_xyzz_zz_xx, t_0_xyzz_zz_xy,\
                           t_0_xyzz_zz_xz, t_0_xyzz_zz_yy, t_0_xyzz_zz_yz, t_0_xyzz_zz_zz,\
                           t_0_yzz_xx_xx, t_0_yzz_xx_xy, t_0_yzz_xx_xz, t_0_yzz_xx_yy,\
                           t_0_yzz_xx_yz, t_0_yzz_xx_zz, t_0_yzz_xy_xx, t_0_yzz_xy_xy,\
                           t_0_yzz_xy_xz, t_0_yzz_xy_yy, t_0_yzz_xy_yz, t_0_yzz_xy_zz,\
                           t_0_yzz_xz_xx, t_0_yzz_xz_xy, t_0_yzz_xz_xz, t_0_yzz_xz_yy,\
                           t_0_yzz_xz_yz, t_0_yzz_xz_zz, t_0_yzz_yy_xx, t_0_yzz_yy_xy,\
                           t_0_yzz_yy_xz, t_0_yzz_yy_yy, t_0_yzz_yy_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_xx, t_0_yzz_yz_xy, t_0_yzz_yz_xz, t_0_yzz_yz_yy,\
                           t_0_yzz_yz_yz, t_0_yzz_yz_zz, t_0_yzz_zz_xx, t_0_yzz_zz_xy,\
                           t_0_yzz_zz_xz, t_0_yzz_zz_yy, t_0_yzz_zz_yz, t_0_yzz_zz_zz,\
                           t_x_yzz_xx_xx, t_x_yzz_xx_xy, t_x_yzz_xx_xz, t_x_yzz_xx_yy,\
                           t_x_yzz_xx_yz, t_x_yzz_xx_zz, t_x_yzz_xy_xx, t_x_yzz_xy_xy,\
                           t_x_yzz_xy_xz, t_x_yzz_xy_yy, t_x_yzz_xy_yz, t_x_yzz_xy_zz,\
                           t_x_yzz_xz_xx, t_x_yzz_xz_xy, t_x_yzz_xz_xz, t_x_yzz_xz_yy,\
                           t_x_yzz_xz_yz, t_x_yzz_xz_zz, t_x_yzz_yy_xx, t_x_yzz_yy_xy,\
                           t_x_yzz_yy_xz, t_x_yzz_yy_yy, t_x_yzz_yy_yz, t_x_yzz_yy_zz,\
                           t_x_yzz_yz_xx, t_x_yzz_yz_xy, t_x_yzz_yz_xz, t_x_yzz_yz_yy,\
                           t_x_yzz_yz_yz, t_x_yzz_yz_zz, t_x_yzz_zz_xx, t_x_yzz_zz_xy,\
                           t_x_yzz_zz_xz, t_x_yzz_zz_yy, t_x_yzz_zz_yz, t_x_yzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yzz_zz_zz[i] = t_0_xyzz_zz_zz[i] - rab_x[i] * t_0_yzz_zz_zz[i];

        t_x_yzz_zz_yz[i] = t_0_xyzz_zz_yz[i] - rab_x[i] * t_0_yzz_zz_yz[i];

        t_x_yzz_zz_yy[i] = t_0_xyzz_zz_yy[i] - rab_x[i] * t_0_yzz_zz_yy[i];

        t_x_yzz_zz_xz[i] = t_0_xyzz_zz_xz[i] - rab_x[i] * t_0_yzz_zz_xz[i];

        t_x_yzz_zz_xy[i] = t_0_xyzz_zz_xy[i] - rab_x[i] * t_0_yzz_zz_xy[i];

        t_x_yzz_zz_xx[i] = t_0_xyzz_zz_xx[i] - rab_x[i] * t_0_yzz_zz_xx[i];

        t_x_yzz_yz_zz[i] = t_0_xyzz_yz_zz[i] - rab_x[i] * t_0_yzz_yz_zz[i];

        t_x_yzz_yz_yz[i] = t_0_xyzz_yz_yz[i] - rab_x[i] * t_0_yzz_yz_yz[i];

        t_x_yzz_yz_yy[i] = t_0_xyzz_yz_yy[i] - rab_x[i] * t_0_yzz_yz_yy[i];

        t_x_yzz_yz_xz[i] = t_0_xyzz_yz_xz[i] - rab_x[i] * t_0_yzz_yz_xz[i];

        t_x_yzz_yz_xy[i] = t_0_xyzz_yz_xy[i] - rab_x[i] * t_0_yzz_yz_xy[i];

        t_x_yzz_yz_xx[i] = t_0_xyzz_yz_xx[i] - rab_x[i] * t_0_yzz_yz_xx[i];

        t_x_yzz_yy_zz[i] = t_0_xyzz_yy_zz[i] - rab_x[i] * t_0_yzz_yy_zz[i];

        t_x_yzz_yy_yz[i] = t_0_xyzz_yy_yz[i] - rab_x[i] * t_0_yzz_yy_yz[i];

        t_x_yzz_yy_yy[i] = t_0_xyzz_yy_yy[i] - rab_x[i] * t_0_yzz_yy_yy[i];

        t_x_yzz_yy_xz[i] = t_0_xyzz_yy_xz[i] - rab_x[i] * t_0_yzz_yy_xz[i];

        t_x_yzz_yy_xy[i] = t_0_xyzz_yy_xy[i] - rab_x[i] * t_0_yzz_yy_xy[i];

        t_x_yzz_yy_xx[i] = t_0_xyzz_yy_xx[i] - rab_x[i] * t_0_yzz_yy_xx[i];

        t_x_yzz_xz_zz[i] = t_0_xyzz_xz_zz[i] - rab_x[i] * t_0_yzz_xz_zz[i];

        t_x_yzz_xz_yz[i] = t_0_xyzz_xz_yz[i] - rab_x[i] * t_0_yzz_xz_yz[i];

        t_x_yzz_xz_yy[i] = t_0_xyzz_xz_yy[i] - rab_x[i] * t_0_yzz_xz_yy[i];

        t_x_yzz_xz_xz[i] = t_0_xyzz_xz_xz[i] - rab_x[i] * t_0_yzz_xz_xz[i];

        t_x_yzz_xz_xy[i] = t_0_xyzz_xz_xy[i] - rab_x[i] * t_0_yzz_xz_xy[i];

        t_x_yzz_xz_xx[i] = t_0_xyzz_xz_xx[i] - rab_x[i] * t_0_yzz_xz_xx[i];

        t_x_yzz_xy_zz[i] = t_0_xyzz_xy_zz[i] - rab_x[i] * t_0_yzz_xy_zz[i];

        t_x_yzz_xy_yz[i] = t_0_xyzz_xy_yz[i] - rab_x[i] * t_0_yzz_xy_yz[i];

        t_x_yzz_xy_yy[i] = t_0_xyzz_xy_yy[i] - rab_x[i] * t_0_yzz_xy_yy[i];

        t_x_yzz_xy_xz[i] = t_0_xyzz_xy_xz[i] - rab_x[i] * t_0_yzz_xy_xz[i];

        t_x_yzz_xy_xy[i] = t_0_xyzz_xy_xy[i] - rab_x[i] * t_0_yzz_xy_xy[i];

        t_x_yzz_xy_xx[i] = t_0_xyzz_xy_xx[i] - rab_x[i] * t_0_yzz_xy_xx[i];

        t_x_yzz_xx_zz[i] = t_0_xyzz_xx_zz[i] - rab_x[i] * t_0_yzz_xx_zz[i];

        t_x_yzz_xx_yz[i] = t_0_xyzz_xx_yz[i] - rab_x[i] * t_0_yzz_xx_yz[i];

        t_x_yzz_xx_yy[i] = t_0_xyzz_xx_yy[i] - rab_x[i] * t_0_yzz_xx_yy[i];

        t_x_yzz_xx_xz[i] = t_0_xyzz_xx_xz[i] - rab_x[i] * t_0_yzz_xx_xz[i];

        t_x_yzz_xx_xy[i] = t_0_xyzz_xx_xy[i] - rab_x[i] * t_0_yzz_xx_xy[i];

        t_x_yzz_xx_xx[i] = t_0_xyzz_xx_xx[i] - rab_x[i] * t_0_yzz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyyz_xx_xx, t_0_xyyz_xx_xy, t_0_xyyz_xx_xz, t_0_xyyz_xx_yy,\
                           t_0_xyyz_xx_yz, t_0_xyyz_xx_zz, t_0_xyyz_xy_xx, t_0_xyyz_xy_xy,\
                           t_0_xyyz_xy_xz, t_0_xyyz_xy_yy, t_0_xyyz_xy_yz, t_0_xyyz_xy_zz,\
                           t_0_xyyz_xz_xx, t_0_xyyz_xz_xy, t_0_xyyz_xz_xz, t_0_xyyz_xz_yy,\
                           t_0_xyyz_xz_yz, t_0_xyyz_xz_zz, t_0_xyyz_yy_xx, t_0_xyyz_yy_xy,\
                           t_0_xyyz_yy_xz, t_0_xyyz_yy_yy, t_0_xyyz_yy_yz, t_0_xyyz_yy_zz,\
                           t_0_xyyz_yz_xx, t_0_xyyz_yz_xy, t_0_xyyz_yz_xz, t_0_xyyz_yz_yy,\
                           t_0_xyyz_yz_yz, t_0_xyyz_yz_zz, t_0_xyyz_zz_xx, t_0_xyyz_zz_xy,\
                           t_0_xyyz_zz_xz, t_0_xyyz_zz_yy, t_0_xyyz_zz_yz, t_0_xyyz_zz_zz,\
                           t_0_yyz_xx_xx, t_0_yyz_xx_xy, t_0_yyz_xx_xz, t_0_yyz_xx_yy,\
                           t_0_yyz_xx_yz, t_0_yyz_xx_zz, t_0_yyz_xy_xx, t_0_yyz_xy_xy,\
                           t_0_yyz_xy_xz, t_0_yyz_xy_yy, t_0_yyz_xy_yz, t_0_yyz_xy_zz,\
                           t_0_yyz_xz_xx, t_0_yyz_xz_xy, t_0_yyz_xz_xz, t_0_yyz_xz_yy,\
                           t_0_yyz_xz_yz, t_0_yyz_xz_zz, t_0_yyz_yy_xx, t_0_yyz_yy_xy,\
                           t_0_yyz_yy_xz, t_0_yyz_yy_yy, t_0_yyz_yy_yz, t_0_yyz_yy_zz,\
                           t_0_yyz_yz_xx, t_0_yyz_yz_xy, t_0_yyz_yz_xz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_yz_zz, t_0_yyz_zz_xx, t_0_yyz_zz_xy,\
                           t_0_yyz_zz_xz, t_0_yyz_zz_yy, t_0_yyz_zz_yz, t_0_yyz_zz_zz,\
                           t_x_yyz_xx_xx, t_x_yyz_xx_xy, t_x_yyz_xx_xz, t_x_yyz_xx_yy,\
                           t_x_yyz_xx_yz, t_x_yyz_xx_zz, t_x_yyz_xy_xx, t_x_yyz_xy_xy,\
                           t_x_yyz_xy_xz, t_x_yyz_xy_yy, t_x_yyz_xy_yz, t_x_yyz_xy_zz,\
                           t_x_yyz_xz_xx, t_x_yyz_xz_xy, t_x_yyz_xz_xz, t_x_yyz_xz_yy,\
                           t_x_yyz_xz_yz, t_x_yyz_xz_zz, t_x_yyz_yy_xx, t_x_yyz_yy_xy,\
                           t_x_yyz_yy_xz, t_x_yyz_yy_yy, t_x_yyz_yy_yz, t_x_yyz_yy_zz,\
                           t_x_yyz_yz_xx, t_x_yyz_yz_xy, t_x_yyz_yz_xz, t_x_yyz_yz_yy,\
                           t_x_yyz_yz_yz, t_x_yyz_yz_zz, t_x_yyz_zz_xx, t_x_yyz_zz_xy,\
                           t_x_yyz_zz_xz, t_x_yyz_zz_yy, t_x_yyz_zz_yz, t_x_yyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yyz_zz_zz[i] = t_0_xyyz_zz_zz[i] - rab_x[i] * t_0_yyz_zz_zz[i];

        t_x_yyz_zz_yz[i] = t_0_xyyz_zz_yz[i] - rab_x[i] * t_0_yyz_zz_yz[i];

        t_x_yyz_zz_yy[i] = t_0_xyyz_zz_yy[i] - rab_x[i] * t_0_yyz_zz_yy[i];

        t_x_yyz_zz_xz[i] = t_0_xyyz_zz_xz[i] - rab_x[i] * t_0_yyz_zz_xz[i];

        t_x_yyz_zz_xy[i] = t_0_xyyz_zz_xy[i] - rab_x[i] * t_0_yyz_zz_xy[i];

        t_x_yyz_zz_xx[i] = t_0_xyyz_zz_xx[i] - rab_x[i] * t_0_yyz_zz_xx[i];

        t_x_yyz_yz_zz[i] = t_0_xyyz_yz_zz[i] - rab_x[i] * t_0_yyz_yz_zz[i];

        t_x_yyz_yz_yz[i] = t_0_xyyz_yz_yz[i] - rab_x[i] * t_0_yyz_yz_yz[i];

        t_x_yyz_yz_yy[i] = t_0_xyyz_yz_yy[i] - rab_x[i] * t_0_yyz_yz_yy[i];

        t_x_yyz_yz_xz[i] = t_0_xyyz_yz_xz[i] - rab_x[i] * t_0_yyz_yz_xz[i];

        t_x_yyz_yz_xy[i] = t_0_xyyz_yz_xy[i] - rab_x[i] * t_0_yyz_yz_xy[i];

        t_x_yyz_yz_xx[i] = t_0_xyyz_yz_xx[i] - rab_x[i] * t_0_yyz_yz_xx[i];

        t_x_yyz_yy_zz[i] = t_0_xyyz_yy_zz[i] - rab_x[i] * t_0_yyz_yy_zz[i];

        t_x_yyz_yy_yz[i] = t_0_xyyz_yy_yz[i] - rab_x[i] * t_0_yyz_yy_yz[i];

        t_x_yyz_yy_yy[i] = t_0_xyyz_yy_yy[i] - rab_x[i] * t_0_yyz_yy_yy[i];

        t_x_yyz_yy_xz[i] = t_0_xyyz_yy_xz[i] - rab_x[i] * t_0_yyz_yy_xz[i];

        t_x_yyz_yy_xy[i] = t_0_xyyz_yy_xy[i] - rab_x[i] * t_0_yyz_yy_xy[i];

        t_x_yyz_yy_xx[i] = t_0_xyyz_yy_xx[i] - rab_x[i] * t_0_yyz_yy_xx[i];

        t_x_yyz_xz_zz[i] = t_0_xyyz_xz_zz[i] - rab_x[i] * t_0_yyz_xz_zz[i];

        t_x_yyz_xz_yz[i] = t_0_xyyz_xz_yz[i] - rab_x[i] * t_0_yyz_xz_yz[i];

        t_x_yyz_xz_yy[i] = t_0_xyyz_xz_yy[i] - rab_x[i] * t_0_yyz_xz_yy[i];

        t_x_yyz_xz_xz[i] = t_0_xyyz_xz_xz[i] - rab_x[i] * t_0_yyz_xz_xz[i];

        t_x_yyz_xz_xy[i] = t_0_xyyz_xz_xy[i] - rab_x[i] * t_0_yyz_xz_xy[i];

        t_x_yyz_xz_xx[i] = t_0_xyyz_xz_xx[i] - rab_x[i] * t_0_yyz_xz_xx[i];

        t_x_yyz_xy_zz[i] = t_0_xyyz_xy_zz[i] - rab_x[i] * t_0_yyz_xy_zz[i];

        t_x_yyz_xy_yz[i] = t_0_xyyz_xy_yz[i] - rab_x[i] * t_0_yyz_xy_yz[i];

        t_x_yyz_xy_yy[i] = t_0_xyyz_xy_yy[i] - rab_x[i] * t_0_yyz_xy_yy[i];

        t_x_yyz_xy_xz[i] = t_0_xyyz_xy_xz[i] - rab_x[i] * t_0_yyz_xy_xz[i];

        t_x_yyz_xy_xy[i] = t_0_xyyz_xy_xy[i] - rab_x[i] * t_0_yyz_xy_xy[i];

        t_x_yyz_xy_xx[i] = t_0_xyyz_xy_xx[i] - rab_x[i] * t_0_yyz_xy_xx[i];

        t_x_yyz_xx_zz[i] = t_0_xyyz_xx_zz[i] - rab_x[i] * t_0_yyz_xx_zz[i];

        t_x_yyz_xx_yz[i] = t_0_xyyz_xx_yz[i] - rab_x[i] * t_0_yyz_xx_yz[i];

        t_x_yyz_xx_yy[i] = t_0_xyyz_xx_yy[i] - rab_x[i] * t_0_yyz_xx_yy[i];

        t_x_yyz_xx_xz[i] = t_0_xyyz_xx_xz[i] - rab_x[i] * t_0_yyz_xx_xz[i];

        t_x_yyz_xx_xy[i] = t_0_xyyz_xx_xy[i] - rab_x[i] * t_0_yyz_xx_xy[i];

        t_x_yyz_xx_xx[i] = t_0_xyyz_xx_xx[i] - rab_x[i] * t_0_yyz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyyy_xx_xx, t_0_xyyy_xx_xy, t_0_xyyy_xx_xz, t_0_xyyy_xx_yy,\
                           t_0_xyyy_xx_yz, t_0_xyyy_xx_zz, t_0_xyyy_xy_xx, t_0_xyyy_xy_xy,\
                           t_0_xyyy_xy_xz, t_0_xyyy_xy_yy, t_0_xyyy_xy_yz, t_0_xyyy_xy_zz,\
                           t_0_xyyy_xz_xx, t_0_xyyy_xz_xy, t_0_xyyy_xz_xz, t_0_xyyy_xz_yy,\
                           t_0_xyyy_xz_yz, t_0_xyyy_xz_zz, t_0_xyyy_yy_xx, t_0_xyyy_yy_xy,\
                           t_0_xyyy_yy_xz, t_0_xyyy_yy_yy, t_0_xyyy_yy_yz, t_0_xyyy_yy_zz,\
                           t_0_xyyy_yz_xx, t_0_xyyy_yz_xy, t_0_xyyy_yz_xz, t_0_xyyy_yz_yy,\
                           t_0_xyyy_yz_yz, t_0_xyyy_yz_zz, t_0_xyyy_zz_xx, t_0_xyyy_zz_xy,\
                           t_0_xyyy_zz_xz, t_0_xyyy_zz_yy, t_0_xyyy_zz_yz, t_0_xyyy_zz_zz,\
                           t_0_yyy_xx_xx, t_0_yyy_xx_xy, t_0_yyy_xx_xz, t_0_yyy_xx_yy,\
                           t_0_yyy_xx_yz, t_0_yyy_xx_zz, t_0_yyy_xy_xx, t_0_yyy_xy_xy,\
                           t_0_yyy_xy_xz, t_0_yyy_xy_yy, t_0_yyy_xy_yz, t_0_yyy_xy_zz,\
                           t_0_yyy_xz_xx, t_0_yyy_xz_xy, t_0_yyy_xz_xz, t_0_yyy_xz_yy,\
                           t_0_yyy_xz_yz, t_0_yyy_xz_zz, t_0_yyy_yy_xx, t_0_yyy_yy_xy,\
                           t_0_yyy_yy_xz, t_0_yyy_yy_yy, t_0_yyy_yy_yz, t_0_yyy_yy_zz,\
                           t_0_yyy_yz_xx, t_0_yyy_yz_xy, t_0_yyy_yz_xz, t_0_yyy_yz_yy,\
                           t_0_yyy_yz_yz, t_0_yyy_yz_zz, t_0_yyy_zz_xx, t_0_yyy_zz_xy,\
                           t_0_yyy_zz_xz, t_0_yyy_zz_yy, t_0_yyy_zz_yz, t_0_yyy_zz_zz,\
                           t_x_yyy_xx_xx, t_x_yyy_xx_xy, t_x_yyy_xx_xz, t_x_yyy_xx_yy,\
                           t_x_yyy_xx_yz, t_x_yyy_xx_zz, t_x_yyy_xy_xx, t_x_yyy_xy_xy,\
                           t_x_yyy_xy_xz, t_x_yyy_xy_yy, t_x_yyy_xy_yz, t_x_yyy_xy_zz,\
                           t_x_yyy_xz_xx, t_x_yyy_xz_xy, t_x_yyy_xz_xz, t_x_yyy_xz_yy,\
                           t_x_yyy_xz_yz, t_x_yyy_xz_zz, t_x_yyy_yy_xx, t_x_yyy_yy_xy,\
                           t_x_yyy_yy_xz, t_x_yyy_yy_yy, t_x_yyy_yy_yz, t_x_yyy_yy_zz,\
                           t_x_yyy_yz_xx, t_x_yyy_yz_xy, t_x_yyy_yz_xz, t_x_yyy_yz_yy,\
                           t_x_yyy_yz_yz, t_x_yyy_yz_zz, t_x_yyy_zz_xx, t_x_yyy_zz_xy,\
                           t_x_yyy_zz_xz, t_x_yyy_zz_yy, t_x_yyy_zz_yz, t_x_yyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yyy_zz_zz[i] = t_0_xyyy_zz_zz[i] - rab_x[i] * t_0_yyy_zz_zz[i];

        t_x_yyy_zz_yz[i] = t_0_xyyy_zz_yz[i] - rab_x[i] * t_0_yyy_zz_yz[i];

        t_x_yyy_zz_yy[i] = t_0_xyyy_zz_yy[i] - rab_x[i] * t_0_yyy_zz_yy[i];

        t_x_yyy_zz_xz[i] = t_0_xyyy_zz_xz[i] - rab_x[i] * t_0_yyy_zz_xz[i];

        t_x_yyy_zz_xy[i] = t_0_xyyy_zz_xy[i] - rab_x[i] * t_0_yyy_zz_xy[i];

        t_x_yyy_zz_xx[i] = t_0_xyyy_zz_xx[i] - rab_x[i] * t_0_yyy_zz_xx[i];

        t_x_yyy_yz_zz[i] = t_0_xyyy_yz_zz[i] - rab_x[i] * t_0_yyy_yz_zz[i];

        t_x_yyy_yz_yz[i] = t_0_xyyy_yz_yz[i] - rab_x[i] * t_0_yyy_yz_yz[i];

        t_x_yyy_yz_yy[i] = t_0_xyyy_yz_yy[i] - rab_x[i] * t_0_yyy_yz_yy[i];

        t_x_yyy_yz_xz[i] = t_0_xyyy_yz_xz[i] - rab_x[i] * t_0_yyy_yz_xz[i];

        t_x_yyy_yz_xy[i] = t_0_xyyy_yz_xy[i] - rab_x[i] * t_0_yyy_yz_xy[i];

        t_x_yyy_yz_xx[i] = t_0_xyyy_yz_xx[i] - rab_x[i] * t_0_yyy_yz_xx[i];

        t_x_yyy_yy_zz[i] = t_0_xyyy_yy_zz[i] - rab_x[i] * t_0_yyy_yy_zz[i];

        t_x_yyy_yy_yz[i] = t_0_xyyy_yy_yz[i] - rab_x[i] * t_0_yyy_yy_yz[i];

        t_x_yyy_yy_yy[i] = t_0_xyyy_yy_yy[i] - rab_x[i] * t_0_yyy_yy_yy[i];

        t_x_yyy_yy_xz[i] = t_0_xyyy_yy_xz[i] - rab_x[i] * t_0_yyy_yy_xz[i];

        t_x_yyy_yy_xy[i] = t_0_xyyy_yy_xy[i] - rab_x[i] * t_0_yyy_yy_xy[i];

        t_x_yyy_yy_xx[i] = t_0_xyyy_yy_xx[i] - rab_x[i] * t_0_yyy_yy_xx[i];

        t_x_yyy_xz_zz[i] = t_0_xyyy_xz_zz[i] - rab_x[i] * t_0_yyy_xz_zz[i];

        t_x_yyy_xz_yz[i] = t_0_xyyy_xz_yz[i] - rab_x[i] * t_0_yyy_xz_yz[i];

        t_x_yyy_xz_yy[i] = t_0_xyyy_xz_yy[i] - rab_x[i] * t_0_yyy_xz_yy[i];

        t_x_yyy_xz_xz[i] = t_0_xyyy_xz_xz[i] - rab_x[i] * t_0_yyy_xz_xz[i];

        t_x_yyy_xz_xy[i] = t_0_xyyy_xz_xy[i] - rab_x[i] * t_0_yyy_xz_xy[i];

        t_x_yyy_xz_xx[i] = t_0_xyyy_xz_xx[i] - rab_x[i] * t_0_yyy_xz_xx[i];

        t_x_yyy_xy_zz[i] = t_0_xyyy_xy_zz[i] - rab_x[i] * t_0_yyy_xy_zz[i];

        t_x_yyy_xy_yz[i] = t_0_xyyy_xy_yz[i] - rab_x[i] * t_0_yyy_xy_yz[i];

        t_x_yyy_xy_yy[i] = t_0_xyyy_xy_yy[i] - rab_x[i] * t_0_yyy_xy_yy[i];

        t_x_yyy_xy_xz[i] = t_0_xyyy_xy_xz[i] - rab_x[i] * t_0_yyy_xy_xz[i];

        t_x_yyy_xy_xy[i] = t_0_xyyy_xy_xy[i] - rab_x[i] * t_0_yyy_xy_xy[i];

        t_x_yyy_xy_xx[i] = t_0_xyyy_xy_xx[i] - rab_x[i] * t_0_yyy_xy_xx[i];

        t_x_yyy_xx_zz[i] = t_0_xyyy_xx_zz[i] - rab_x[i] * t_0_yyy_xx_zz[i];

        t_x_yyy_xx_yz[i] = t_0_xyyy_xx_yz[i] - rab_x[i] * t_0_yyy_xx_yz[i];

        t_x_yyy_xx_yy[i] = t_0_xyyy_xx_yy[i] - rab_x[i] * t_0_yyy_xx_yy[i];

        t_x_yyy_xx_xz[i] = t_0_xyyy_xx_xz[i] - rab_x[i] * t_0_yyy_xx_xz[i];

        t_x_yyy_xx_xy[i] = t_0_xyyy_xx_xy[i] - rab_x[i] * t_0_yyy_xx_xy[i];

        t_x_yyy_xx_xx[i] = t_0_xyyy_xx_xx[i] - rab_x[i] * t_0_yyy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxzz_xx_xx, t_0_xxzz_xx_xy, t_0_xxzz_xx_xz, t_0_xxzz_xx_yy,\
                           t_0_xxzz_xx_yz, t_0_xxzz_xx_zz, t_0_xxzz_xy_xx, t_0_xxzz_xy_xy,\
                           t_0_xxzz_xy_xz, t_0_xxzz_xy_yy, t_0_xxzz_xy_yz, t_0_xxzz_xy_zz,\
                           t_0_xxzz_xz_xx, t_0_xxzz_xz_xy, t_0_xxzz_xz_xz, t_0_xxzz_xz_yy,\
                           t_0_xxzz_xz_yz, t_0_xxzz_xz_zz, t_0_xxzz_yy_xx, t_0_xxzz_yy_xy,\
                           t_0_xxzz_yy_xz, t_0_xxzz_yy_yy, t_0_xxzz_yy_yz, t_0_xxzz_yy_zz,\
                           t_0_xxzz_yz_xx, t_0_xxzz_yz_xy, t_0_xxzz_yz_xz, t_0_xxzz_yz_yy,\
                           t_0_xxzz_yz_yz, t_0_xxzz_yz_zz, t_0_xxzz_zz_xx, t_0_xxzz_zz_xy,\
                           t_0_xxzz_zz_xz, t_0_xxzz_zz_yy, t_0_xxzz_zz_yz, t_0_xxzz_zz_zz,\
                           t_0_xzz_xx_xx, t_0_xzz_xx_xy, t_0_xzz_xx_xz, t_0_xzz_xx_yy,\
                           t_0_xzz_xx_yz, t_0_xzz_xx_zz, t_0_xzz_xy_xx, t_0_xzz_xy_xy,\
                           t_0_xzz_xy_xz, t_0_xzz_xy_yy, t_0_xzz_xy_yz, t_0_xzz_xy_zz,\
                           t_0_xzz_xz_xx, t_0_xzz_xz_xy, t_0_xzz_xz_xz, t_0_xzz_xz_yy,\
                           t_0_xzz_xz_yz, t_0_xzz_xz_zz, t_0_xzz_yy_xx, t_0_xzz_yy_xy,\
                           t_0_xzz_yy_xz, t_0_xzz_yy_yy, t_0_xzz_yy_yz, t_0_xzz_yy_zz,\
                           t_0_xzz_yz_xx, t_0_xzz_yz_xy, t_0_xzz_yz_xz, t_0_xzz_yz_yy,\
                           t_0_xzz_yz_yz, t_0_xzz_yz_zz, t_0_xzz_zz_xx, t_0_xzz_zz_xy,\
                           t_0_xzz_zz_xz, t_0_xzz_zz_yy, t_0_xzz_zz_yz, t_0_xzz_zz_zz,\
                           t_x_xzz_xx_xx, t_x_xzz_xx_xy, t_x_xzz_xx_xz, t_x_xzz_xx_yy,\
                           t_x_xzz_xx_yz, t_x_xzz_xx_zz, t_x_xzz_xy_xx, t_x_xzz_xy_xy,\
                           t_x_xzz_xy_xz, t_x_xzz_xy_yy, t_x_xzz_xy_yz, t_x_xzz_xy_zz,\
                           t_x_xzz_xz_xx, t_x_xzz_xz_xy, t_x_xzz_xz_xz, t_x_xzz_xz_yy,\
                           t_x_xzz_xz_yz, t_x_xzz_xz_zz, t_x_xzz_yy_xx, t_x_xzz_yy_xy,\
                           t_x_xzz_yy_xz, t_x_xzz_yy_yy, t_x_xzz_yy_yz, t_x_xzz_yy_zz,\
                           t_x_xzz_yz_xx, t_x_xzz_yz_xy, t_x_xzz_yz_xz, t_x_xzz_yz_yy,\
                           t_x_xzz_yz_yz, t_x_xzz_yz_zz, t_x_xzz_zz_xx, t_x_xzz_zz_xy,\
                           t_x_xzz_zz_xz, t_x_xzz_zz_yy, t_x_xzz_zz_yz, t_x_xzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xzz_zz_zz[i] = t_0_xxzz_zz_zz[i] - rab_x[i] * t_0_xzz_zz_zz[i];

        t_x_xzz_zz_yz[i] = t_0_xxzz_zz_yz[i] - rab_x[i] * t_0_xzz_zz_yz[i];

        t_x_xzz_zz_yy[i] = t_0_xxzz_zz_yy[i] - rab_x[i] * t_0_xzz_zz_yy[i];

        t_x_xzz_zz_xz[i] = t_0_xxzz_zz_xz[i] - rab_x[i] * t_0_xzz_zz_xz[i];

        t_x_xzz_zz_xy[i] = t_0_xxzz_zz_xy[i] - rab_x[i] * t_0_xzz_zz_xy[i];

        t_x_xzz_zz_xx[i] = t_0_xxzz_zz_xx[i] - rab_x[i] * t_0_xzz_zz_xx[i];

        t_x_xzz_yz_zz[i] = t_0_xxzz_yz_zz[i] - rab_x[i] * t_0_xzz_yz_zz[i];

        t_x_xzz_yz_yz[i] = t_0_xxzz_yz_yz[i] - rab_x[i] * t_0_xzz_yz_yz[i];

        t_x_xzz_yz_yy[i] = t_0_xxzz_yz_yy[i] - rab_x[i] * t_0_xzz_yz_yy[i];

        t_x_xzz_yz_xz[i] = t_0_xxzz_yz_xz[i] - rab_x[i] * t_0_xzz_yz_xz[i];

        t_x_xzz_yz_xy[i] = t_0_xxzz_yz_xy[i] - rab_x[i] * t_0_xzz_yz_xy[i];

        t_x_xzz_yz_xx[i] = t_0_xxzz_yz_xx[i] - rab_x[i] * t_0_xzz_yz_xx[i];

        t_x_xzz_yy_zz[i] = t_0_xxzz_yy_zz[i] - rab_x[i] * t_0_xzz_yy_zz[i];

        t_x_xzz_yy_yz[i] = t_0_xxzz_yy_yz[i] - rab_x[i] * t_0_xzz_yy_yz[i];

        t_x_xzz_yy_yy[i] = t_0_xxzz_yy_yy[i] - rab_x[i] * t_0_xzz_yy_yy[i];

        t_x_xzz_yy_xz[i] = t_0_xxzz_yy_xz[i] - rab_x[i] * t_0_xzz_yy_xz[i];

        t_x_xzz_yy_xy[i] = t_0_xxzz_yy_xy[i] - rab_x[i] * t_0_xzz_yy_xy[i];

        t_x_xzz_yy_xx[i] = t_0_xxzz_yy_xx[i] - rab_x[i] * t_0_xzz_yy_xx[i];

        t_x_xzz_xz_zz[i] = t_0_xxzz_xz_zz[i] - rab_x[i] * t_0_xzz_xz_zz[i];

        t_x_xzz_xz_yz[i] = t_0_xxzz_xz_yz[i] - rab_x[i] * t_0_xzz_xz_yz[i];

        t_x_xzz_xz_yy[i] = t_0_xxzz_xz_yy[i] - rab_x[i] * t_0_xzz_xz_yy[i];

        t_x_xzz_xz_xz[i] = t_0_xxzz_xz_xz[i] - rab_x[i] * t_0_xzz_xz_xz[i];

        t_x_xzz_xz_xy[i] = t_0_xxzz_xz_xy[i] - rab_x[i] * t_0_xzz_xz_xy[i];

        t_x_xzz_xz_xx[i] = t_0_xxzz_xz_xx[i] - rab_x[i] * t_0_xzz_xz_xx[i];

        t_x_xzz_xy_zz[i] = t_0_xxzz_xy_zz[i] - rab_x[i] * t_0_xzz_xy_zz[i];

        t_x_xzz_xy_yz[i] = t_0_xxzz_xy_yz[i] - rab_x[i] * t_0_xzz_xy_yz[i];

        t_x_xzz_xy_yy[i] = t_0_xxzz_xy_yy[i] - rab_x[i] * t_0_xzz_xy_yy[i];

        t_x_xzz_xy_xz[i] = t_0_xxzz_xy_xz[i] - rab_x[i] * t_0_xzz_xy_xz[i];

        t_x_xzz_xy_xy[i] = t_0_xxzz_xy_xy[i] - rab_x[i] * t_0_xzz_xy_xy[i];

        t_x_xzz_xy_xx[i] = t_0_xxzz_xy_xx[i] - rab_x[i] * t_0_xzz_xy_xx[i];

        t_x_xzz_xx_zz[i] = t_0_xxzz_xx_zz[i] - rab_x[i] * t_0_xzz_xx_zz[i];

        t_x_xzz_xx_yz[i] = t_0_xxzz_xx_yz[i] - rab_x[i] * t_0_xzz_xx_yz[i];

        t_x_xzz_xx_yy[i] = t_0_xxzz_xx_yy[i] - rab_x[i] * t_0_xzz_xx_yy[i];

        t_x_xzz_xx_xz[i] = t_0_xxzz_xx_xz[i] - rab_x[i] * t_0_xzz_xx_xz[i];

        t_x_xzz_xx_xy[i] = t_0_xxzz_xx_xy[i] - rab_x[i] * t_0_xzz_xx_xy[i];

        t_x_xzz_xx_xx[i] = t_0_xxzz_xx_xx[i] - rab_x[i] * t_0_xzz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxyz_xx_xx, t_0_xxyz_xx_xy, t_0_xxyz_xx_xz, t_0_xxyz_xx_yy,\
                           t_0_xxyz_xx_yz, t_0_xxyz_xx_zz, t_0_xxyz_xy_xx, t_0_xxyz_xy_xy,\
                           t_0_xxyz_xy_xz, t_0_xxyz_xy_yy, t_0_xxyz_xy_yz, t_0_xxyz_xy_zz,\
                           t_0_xxyz_xz_xx, t_0_xxyz_xz_xy, t_0_xxyz_xz_xz, t_0_xxyz_xz_yy,\
                           t_0_xxyz_xz_yz, t_0_xxyz_xz_zz, t_0_xxyz_yy_xx, t_0_xxyz_yy_xy,\
                           t_0_xxyz_yy_xz, t_0_xxyz_yy_yy, t_0_xxyz_yy_yz, t_0_xxyz_yy_zz,\
                           t_0_xxyz_yz_xx, t_0_xxyz_yz_xy, t_0_xxyz_yz_xz, t_0_xxyz_yz_yy,\
                           t_0_xxyz_yz_yz, t_0_xxyz_yz_zz, t_0_xxyz_zz_xx, t_0_xxyz_zz_xy,\
                           t_0_xxyz_zz_xz, t_0_xxyz_zz_yy, t_0_xxyz_zz_yz, t_0_xxyz_zz_zz,\
                           t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz, t_0_xyz_xx_yy,\
                           t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx, t_0_xyz_xy_xy,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz, t_0_xyz_xy_zz,\
                           t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz, t_0_xyz_xz_yy,\
                           t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz,\
                           t_x_xyz_xx_xx, t_x_xyz_xx_xy, t_x_xyz_xx_xz, t_x_xyz_xx_yy,\
                           t_x_xyz_xx_yz, t_x_xyz_xx_zz, t_x_xyz_xy_xx, t_x_xyz_xy_xy,\
                           t_x_xyz_xy_xz, t_x_xyz_xy_yy, t_x_xyz_xy_yz, t_x_xyz_xy_zz,\
                           t_x_xyz_xz_xx, t_x_xyz_xz_xy, t_x_xyz_xz_xz, t_x_xyz_xz_yy,\
                           t_x_xyz_xz_yz, t_x_xyz_xz_zz, t_x_xyz_yy_xx, t_x_xyz_yy_xy,\
                           t_x_xyz_yy_xz, t_x_xyz_yy_yy, t_x_xyz_yy_yz, t_x_xyz_yy_zz,\
                           t_x_xyz_yz_xx, t_x_xyz_yz_xy, t_x_xyz_yz_xz, t_x_xyz_yz_yy,\
                           t_x_xyz_yz_yz, t_x_xyz_yz_zz, t_x_xyz_zz_xx, t_x_xyz_zz_xy,\
                           t_x_xyz_zz_xz, t_x_xyz_zz_yy, t_x_xyz_zz_yz, t_x_xyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xyz_zz_zz[i] = t_0_xxyz_zz_zz[i] - rab_x[i] * t_0_xyz_zz_zz[i];

        t_x_xyz_zz_yz[i] = t_0_xxyz_zz_yz[i] - rab_x[i] * t_0_xyz_zz_yz[i];

        t_x_xyz_zz_yy[i] = t_0_xxyz_zz_yy[i] - rab_x[i] * t_0_xyz_zz_yy[i];

        t_x_xyz_zz_xz[i] = t_0_xxyz_zz_xz[i] - rab_x[i] * t_0_xyz_zz_xz[i];

        t_x_xyz_zz_xy[i] = t_0_xxyz_zz_xy[i] - rab_x[i] * t_0_xyz_zz_xy[i];

        t_x_xyz_zz_xx[i] = t_0_xxyz_zz_xx[i] - rab_x[i] * t_0_xyz_zz_xx[i];

        t_x_xyz_yz_zz[i] = t_0_xxyz_yz_zz[i] - rab_x[i] * t_0_xyz_yz_zz[i];

        t_x_xyz_yz_yz[i] = t_0_xxyz_yz_yz[i] - rab_x[i] * t_0_xyz_yz_yz[i];

        t_x_xyz_yz_yy[i] = t_0_xxyz_yz_yy[i] - rab_x[i] * t_0_xyz_yz_yy[i];

        t_x_xyz_yz_xz[i] = t_0_xxyz_yz_xz[i] - rab_x[i] * t_0_xyz_yz_xz[i];

        t_x_xyz_yz_xy[i] = t_0_xxyz_yz_xy[i] - rab_x[i] * t_0_xyz_yz_xy[i];

        t_x_xyz_yz_xx[i] = t_0_xxyz_yz_xx[i] - rab_x[i] * t_0_xyz_yz_xx[i];

        t_x_xyz_yy_zz[i] = t_0_xxyz_yy_zz[i] - rab_x[i] * t_0_xyz_yy_zz[i];

        t_x_xyz_yy_yz[i] = t_0_xxyz_yy_yz[i] - rab_x[i] * t_0_xyz_yy_yz[i];

        t_x_xyz_yy_yy[i] = t_0_xxyz_yy_yy[i] - rab_x[i] * t_0_xyz_yy_yy[i];

        t_x_xyz_yy_xz[i] = t_0_xxyz_yy_xz[i] - rab_x[i] * t_0_xyz_yy_xz[i];

        t_x_xyz_yy_xy[i] = t_0_xxyz_yy_xy[i] - rab_x[i] * t_0_xyz_yy_xy[i];

        t_x_xyz_yy_xx[i] = t_0_xxyz_yy_xx[i] - rab_x[i] * t_0_xyz_yy_xx[i];

        t_x_xyz_xz_zz[i] = t_0_xxyz_xz_zz[i] - rab_x[i] * t_0_xyz_xz_zz[i];

        t_x_xyz_xz_yz[i] = t_0_xxyz_xz_yz[i] - rab_x[i] * t_0_xyz_xz_yz[i];

        t_x_xyz_xz_yy[i] = t_0_xxyz_xz_yy[i] - rab_x[i] * t_0_xyz_xz_yy[i];

        t_x_xyz_xz_xz[i] = t_0_xxyz_xz_xz[i] - rab_x[i] * t_0_xyz_xz_xz[i];

        t_x_xyz_xz_xy[i] = t_0_xxyz_xz_xy[i] - rab_x[i] * t_0_xyz_xz_xy[i];

        t_x_xyz_xz_xx[i] = t_0_xxyz_xz_xx[i] - rab_x[i] * t_0_xyz_xz_xx[i];

        t_x_xyz_xy_zz[i] = t_0_xxyz_xy_zz[i] - rab_x[i] * t_0_xyz_xy_zz[i];

        t_x_xyz_xy_yz[i] = t_0_xxyz_xy_yz[i] - rab_x[i] * t_0_xyz_xy_yz[i];

        t_x_xyz_xy_yy[i] = t_0_xxyz_xy_yy[i] - rab_x[i] * t_0_xyz_xy_yy[i];

        t_x_xyz_xy_xz[i] = t_0_xxyz_xy_xz[i] - rab_x[i] * t_0_xyz_xy_xz[i];

        t_x_xyz_xy_xy[i] = t_0_xxyz_xy_xy[i] - rab_x[i] * t_0_xyz_xy_xy[i];

        t_x_xyz_xy_xx[i] = t_0_xxyz_xy_xx[i] - rab_x[i] * t_0_xyz_xy_xx[i];

        t_x_xyz_xx_zz[i] = t_0_xxyz_xx_zz[i] - rab_x[i] * t_0_xyz_xx_zz[i];

        t_x_xyz_xx_yz[i] = t_0_xxyz_xx_yz[i] - rab_x[i] * t_0_xyz_xx_yz[i];

        t_x_xyz_xx_yy[i] = t_0_xxyz_xx_yy[i] - rab_x[i] * t_0_xyz_xx_yy[i];

        t_x_xyz_xx_xz[i] = t_0_xxyz_xx_xz[i] - rab_x[i] * t_0_xyz_xx_xz[i];

        t_x_xyz_xx_xy[i] = t_0_xxyz_xx_xy[i] - rab_x[i] * t_0_xyz_xx_xy[i];

        t_x_xyz_xx_xx[i] = t_0_xxyz_xx_xx[i] - rab_x[i] * t_0_xyz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxyy_xx_xx, t_0_xxyy_xx_xy, t_0_xxyy_xx_xz, t_0_xxyy_xx_yy,\
                           t_0_xxyy_xx_yz, t_0_xxyy_xx_zz, t_0_xxyy_xy_xx, t_0_xxyy_xy_xy,\
                           t_0_xxyy_xy_xz, t_0_xxyy_xy_yy, t_0_xxyy_xy_yz, t_0_xxyy_xy_zz,\
                           t_0_xxyy_xz_xx, t_0_xxyy_xz_xy, t_0_xxyy_xz_xz, t_0_xxyy_xz_yy,\
                           t_0_xxyy_xz_yz, t_0_xxyy_xz_zz, t_0_xxyy_yy_xx, t_0_xxyy_yy_xy,\
                           t_0_xxyy_yy_xz, t_0_xxyy_yy_yy, t_0_xxyy_yy_yz, t_0_xxyy_yy_zz,\
                           t_0_xxyy_yz_xx, t_0_xxyy_yz_xy, t_0_xxyy_yz_xz, t_0_xxyy_yz_yy,\
                           t_0_xxyy_yz_yz, t_0_xxyy_yz_zz, t_0_xxyy_zz_xx, t_0_xxyy_zz_xy,\
                           t_0_xxyy_zz_xz, t_0_xxyy_zz_yy, t_0_xxyy_zz_yz, t_0_xxyy_zz_zz,\
                           t_0_xyy_xx_xx, t_0_xyy_xx_xy, t_0_xyy_xx_xz, t_0_xyy_xx_yy,\
                           t_0_xyy_xx_yz, t_0_xyy_xx_zz, t_0_xyy_xy_xx, t_0_xyy_xy_xy,\
                           t_0_xyy_xy_xz, t_0_xyy_xy_yy, t_0_xyy_xy_yz, t_0_xyy_xy_zz,\
                           t_0_xyy_xz_xx, t_0_xyy_xz_xy, t_0_xyy_xz_xz, t_0_xyy_xz_yy,\
                           t_0_xyy_xz_yz, t_0_xyy_xz_zz, t_0_xyy_yy_xx, t_0_xyy_yy_xy,\
                           t_0_xyy_yy_xz, t_0_xyy_yy_yy, t_0_xyy_yy_yz, t_0_xyy_yy_zz,\
                           t_0_xyy_yz_xx, t_0_xyy_yz_xy, t_0_xyy_yz_xz, t_0_xyy_yz_yy,\
                           t_0_xyy_yz_yz, t_0_xyy_yz_zz, t_0_xyy_zz_xx, t_0_xyy_zz_xy,\
                           t_0_xyy_zz_xz, t_0_xyy_zz_yy, t_0_xyy_zz_yz, t_0_xyy_zz_zz,\
                           t_x_xyy_xx_xx, t_x_xyy_xx_xy, t_x_xyy_xx_xz, t_x_xyy_xx_yy,\
                           t_x_xyy_xx_yz, t_x_xyy_xx_zz, t_x_xyy_xy_xx, t_x_xyy_xy_xy,\
                           t_x_xyy_xy_xz, t_x_xyy_xy_yy, t_x_xyy_xy_yz, t_x_xyy_xy_zz,\
                           t_x_xyy_xz_xx, t_x_xyy_xz_xy, t_x_xyy_xz_xz, t_x_xyy_xz_yy,\
                           t_x_xyy_xz_yz, t_x_xyy_xz_zz, t_x_xyy_yy_xx, t_x_xyy_yy_xy,\
                           t_x_xyy_yy_xz, t_x_xyy_yy_yy, t_x_xyy_yy_yz, t_x_xyy_yy_zz,\
                           t_x_xyy_yz_xx, t_x_xyy_yz_xy, t_x_xyy_yz_xz, t_x_xyy_yz_yy,\
                           t_x_xyy_yz_yz, t_x_xyy_yz_zz, t_x_xyy_zz_xx, t_x_xyy_zz_xy,\
                           t_x_xyy_zz_xz, t_x_xyy_zz_yy, t_x_xyy_zz_yz, t_x_xyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xyy_zz_zz[i] = t_0_xxyy_zz_zz[i] - rab_x[i] * t_0_xyy_zz_zz[i];

        t_x_xyy_zz_yz[i] = t_0_xxyy_zz_yz[i] - rab_x[i] * t_0_xyy_zz_yz[i];

        t_x_xyy_zz_yy[i] = t_0_xxyy_zz_yy[i] - rab_x[i] * t_0_xyy_zz_yy[i];

        t_x_xyy_zz_xz[i] = t_0_xxyy_zz_xz[i] - rab_x[i] * t_0_xyy_zz_xz[i];

        t_x_xyy_zz_xy[i] = t_0_xxyy_zz_xy[i] - rab_x[i] * t_0_xyy_zz_xy[i];

        t_x_xyy_zz_xx[i] = t_0_xxyy_zz_xx[i] - rab_x[i] * t_0_xyy_zz_xx[i];

        t_x_xyy_yz_zz[i] = t_0_xxyy_yz_zz[i] - rab_x[i] * t_0_xyy_yz_zz[i];

        t_x_xyy_yz_yz[i] = t_0_xxyy_yz_yz[i] - rab_x[i] * t_0_xyy_yz_yz[i];

        t_x_xyy_yz_yy[i] = t_0_xxyy_yz_yy[i] - rab_x[i] * t_0_xyy_yz_yy[i];

        t_x_xyy_yz_xz[i] = t_0_xxyy_yz_xz[i] - rab_x[i] * t_0_xyy_yz_xz[i];

        t_x_xyy_yz_xy[i] = t_0_xxyy_yz_xy[i] - rab_x[i] * t_0_xyy_yz_xy[i];

        t_x_xyy_yz_xx[i] = t_0_xxyy_yz_xx[i] - rab_x[i] * t_0_xyy_yz_xx[i];

        t_x_xyy_yy_zz[i] = t_0_xxyy_yy_zz[i] - rab_x[i] * t_0_xyy_yy_zz[i];

        t_x_xyy_yy_yz[i] = t_0_xxyy_yy_yz[i] - rab_x[i] * t_0_xyy_yy_yz[i];

        t_x_xyy_yy_yy[i] = t_0_xxyy_yy_yy[i] - rab_x[i] * t_0_xyy_yy_yy[i];

        t_x_xyy_yy_xz[i] = t_0_xxyy_yy_xz[i] - rab_x[i] * t_0_xyy_yy_xz[i];

        t_x_xyy_yy_xy[i] = t_0_xxyy_yy_xy[i] - rab_x[i] * t_0_xyy_yy_xy[i];

        t_x_xyy_yy_xx[i] = t_0_xxyy_yy_xx[i] - rab_x[i] * t_0_xyy_yy_xx[i];

        t_x_xyy_xz_zz[i] = t_0_xxyy_xz_zz[i] - rab_x[i] * t_0_xyy_xz_zz[i];

        t_x_xyy_xz_yz[i] = t_0_xxyy_xz_yz[i] - rab_x[i] * t_0_xyy_xz_yz[i];

        t_x_xyy_xz_yy[i] = t_0_xxyy_xz_yy[i] - rab_x[i] * t_0_xyy_xz_yy[i];

        t_x_xyy_xz_xz[i] = t_0_xxyy_xz_xz[i] - rab_x[i] * t_0_xyy_xz_xz[i];

        t_x_xyy_xz_xy[i] = t_0_xxyy_xz_xy[i] - rab_x[i] * t_0_xyy_xz_xy[i];

        t_x_xyy_xz_xx[i] = t_0_xxyy_xz_xx[i] - rab_x[i] * t_0_xyy_xz_xx[i];

        t_x_xyy_xy_zz[i] = t_0_xxyy_xy_zz[i] - rab_x[i] * t_0_xyy_xy_zz[i];

        t_x_xyy_xy_yz[i] = t_0_xxyy_xy_yz[i] - rab_x[i] * t_0_xyy_xy_yz[i];

        t_x_xyy_xy_yy[i] = t_0_xxyy_xy_yy[i] - rab_x[i] * t_0_xyy_xy_yy[i];

        t_x_xyy_xy_xz[i] = t_0_xxyy_xy_xz[i] - rab_x[i] * t_0_xyy_xy_xz[i];

        t_x_xyy_xy_xy[i] = t_0_xxyy_xy_xy[i] - rab_x[i] * t_0_xyy_xy_xy[i];

        t_x_xyy_xy_xx[i] = t_0_xxyy_xy_xx[i] - rab_x[i] * t_0_xyy_xy_xx[i];

        t_x_xyy_xx_zz[i] = t_0_xxyy_xx_zz[i] - rab_x[i] * t_0_xyy_xx_zz[i];

        t_x_xyy_xx_yz[i] = t_0_xxyy_xx_yz[i] - rab_x[i] * t_0_xyy_xx_yz[i];

        t_x_xyy_xx_yy[i] = t_0_xxyy_xx_yy[i] - rab_x[i] * t_0_xyy_xx_yy[i];

        t_x_xyy_xx_xz[i] = t_0_xxyy_xx_xz[i] - rab_x[i] * t_0_xyy_xx_xz[i];

        t_x_xyy_xx_xy[i] = t_0_xxyy_xx_xy[i] - rab_x[i] * t_0_xyy_xx_xy[i];

        t_x_xyy_xx_xx[i] = t_0_xxyy_xx_xx[i] - rab_x[i] * t_0_xyy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxxz_xx_xx, t_0_xxxz_xx_xy, t_0_xxxz_xx_xz, t_0_xxxz_xx_yy,\
                           t_0_xxxz_xx_yz, t_0_xxxz_xx_zz, t_0_xxxz_xy_xx, t_0_xxxz_xy_xy,\
                           t_0_xxxz_xy_xz, t_0_xxxz_xy_yy, t_0_xxxz_xy_yz, t_0_xxxz_xy_zz,\
                           t_0_xxxz_xz_xx, t_0_xxxz_xz_xy, t_0_xxxz_xz_xz, t_0_xxxz_xz_yy,\
                           t_0_xxxz_xz_yz, t_0_xxxz_xz_zz, t_0_xxxz_yy_xx, t_0_xxxz_yy_xy,\
                           t_0_xxxz_yy_xz, t_0_xxxz_yy_yy, t_0_xxxz_yy_yz, t_0_xxxz_yy_zz,\
                           t_0_xxxz_yz_xx, t_0_xxxz_yz_xy, t_0_xxxz_yz_xz, t_0_xxxz_yz_yy,\
                           t_0_xxxz_yz_yz, t_0_xxxz_yz_zz, t_0_xxxz_zz_xx, t_0_xxxz_zz_xy,\
                           t_0_xxxz_zz_xz, t_0_xxxz_zz_yy, t_0_xxxz_zz_yz, t_0_xxxz_zz_zz,\
                           t_0_xxz_xx_xx, t_0_xxz_xx_xy, t_0_xxz_xx_xz, t_0_xxz_xx_yy,\
                           t_0_xxz_xx_yz, t_0_xxz_xx_zz, t_0_xxz_xy_xx, t_0_xxz_xy_xy,\
                           t_0_xxz_xy_xz, t_0_xxz_xy_yy, t_0_xxz_xy_yz, t_0_xxz_xy_zz,\
                           t_0_xxz_xz_xx, t_0_xxz_xz_xy, t_0_xxz_xz_xz, t_0_xxz_xz_yy,\
                           t_0_xxz_xz_yz, t_0_xxz_xz_zz, t_0_xxz_yy_xx, t_0_xxz_yy_xy,\
                           t_0_xxz_yy_xz, t_0_xxz_yy_yy, t_0_xxz_yy_yz, t_0_xxz_yy_zz,\
                           t_0_xxz_yz_xx, t_0_xxz_yz_xy, t_0_xxz_yz_xz, t_0_xxz_yz_yy,\
                           t_0_xxz_yz_yz, t_0_xxz_yz_zz, t_0_xxz_zz_xx, t_0_xxz_zz_xy,\
                           t_0_xxz_zz_xz, t_0_xxz_zz_yy, t_0_xxz_zz_yz, t_0_xxz_zz_zz,\
                           t_x_xxz_xx_xx, t_x_xxz_xx_xy, t_x_xxz_xx_xz, t_x_xxz_xx_yy,\
                           t_x_xxz_xx_yz, t_x_xxz_xx_zz, t_x_xxz_xy_xx, t_x_xxz_xy_xy,\
                           t_x_xxz_xy_xz, t_x_xxz_xy_yy, t_x_xxz_xy_yz, t_x_xxz_xy_zz,\
                           t_x_xxz_xz_xx, t_x_xxz_xz_xy, t_x_xxz_xz_xz, t_x_xxz_xz_yy,\
                           t_x_xxz_xz_yz, t_x_xxz_xz_zz, t_x_xxz_yy_xx, t_x_xxz_yy_xy,\
                           t_x_xxz_yy_xz, t_x_xxz_yy_yy, t_x_xxz_yy_yz, t_x_xxz_yy_zz,\
                           t_x_xxz_yz_xx, t_x_xxz_yz_xy, t_x_xxz_yz_xz, t_x_xxz_yz_yy,\
                           t_x_xxz_yz_yz, t_x_xxz_yz_zz, t_x_xxz_zz_xx, t_x_xxz_zz_xy,\
                           t_x_xxz_zz_xz, t_x_xxz_zz_yy, t_x_xxz_zz_yz, t_x_xxz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxz_zz_zz[i] = t_0_xxxz_zz_zz[i] - rab_x[i] * t_0_xxz_zz_zz[i];

        t_x_xxz_zz_yz[i] = t_0_xxxz_zz_yz[i] - rab_x[i] * t_0_xxz_zz_yz[i];

        t_x_xxz_zz_yy[i] = t_0_xxxz_zz_yy[i] - rab_x[i] * t_0_xxz_zz_yy[i];

        t_x_xxz_zz_xz[i] = t_0_xxxz_zz_xz[i] - rab_x[i] * t_0_xxz_zz_xz[i];

        t_x_xxz_zz_xy[i] = t_0_xxxz_zz_xy[i] - rab_x[i] * t_0_xxz_zz_xy[i];

        t_x_xxz_zz_xx[i] = t_0_xxxz_zz_xx[i] - rab_x[i] * t_0_xxz_zz_xx[i];

        t_x_xxz_yz_zz[i] = t_0_xxxz_yz_zz[i] - rab_x[i] * t_0_xxz_yz_zz[i];

        t_x_xxz_yz_yz[i] = t_0_xxxz_yz_yz[i] - rab_x[i] * t_0_xxz_yz_yz[i];

        t_x_xxz_yz_yy[i] = t_0_xxxz_yz_yy[i] - rab_x[i] * t_0_xxz_yz_yy[i];

        t_x_xxz_yz_xz[i] = t_0_xxxz_yz_xz[i] - rab_x[i] * t_0_xxz_yz_xz[i];

        t_x_xxz_yz_xy[i] = t_0_xxxz_yz_xy[i] - rab_x[i] * t_0_xxz_yz_xy[i];

        t_x_xxz_yz_xx[i] = t_0_xxxz_yz_xx[i] - rab_x[i] * t_0_xxz_yz_xx[i];

        t_x_xxz_yy_zz[i] = t_0_xxxz_yy_zz[i] - rab_x[i] * t_0_xxz_yy_zz[i];

        t_x_xxz_yy_yz[i] = t_0_xxxz_yy_yz[i] - rab_x[i] * t_0_xxz_yy_yz[i];

        t_x_xxz_yy_yy[i] = t_0_xxxz_yy_yy[i] - rab_x[i] * t_0_xxz_yy_yy[i];

        t_x_xxz_yy_xz[i] = t_0_xxxz_yy_xz[i] - rab_x[i] * t_0_xxz_yy_xz[i];

        t_x_xxz_yy_xy[i] = t_0_xxxz_yy_xy[i] - rab_x[i] * t_0_xxz_yy_xy[i];

        t_x_xxz_yy_xx[i] = t_0_xxxz_yy_xx[i] - rab_x[i] * t_0_xxz_yy_xx[i];

        t_x_xxz_xz_zz[i] = t_0_xxxz_xz_zz[i] - rab_x[i] * t_0_xxz_xz_zz[i];

        t_x_xxz_xz_yz[i] = t_0_xxxz_xz_yz[i] - rab_x[i] * t_0_xxz_xz_yz[i];

        t_x_xxz_xz_yy[i] = t_0_xxxz_xz_yy[i] - rab_x[i] * t_0_xxz_xz_yy[i];

        t_x_xxz_xz_xz[i] = t_0_xxxz_xz_xz[i] - rab_x[i] * t_0_xxz_xz_xz[i];

        t_x_xxz_xz_xy[i] = t_0_xxxz_xz_xy[i] - rab_x[i] * t_0_xxz_xz_xy[i];

        t_x_xxz_xz_xx[i] = t_0_xxxz_xz_xx[i] - rab_x[i] * t_0_xxz_xz_xx[i];

        t_x_xxz_xy_zz[i] = t_0_xxxz_xy_zz[i] - rab_x[i] * t_0_xxz_xy_zz[i];

        t_x_xxz_xy_yz[i] = t_0_xxxz_xy_yz[i] - rab_x[i] * t_0_xxz_xy_yz[i];

        t_x_xxz_xy_yy[i] = t_0_xxxz_xy_yy[i] - rab_x[i] * t_0_xxz_xy_yy[i];

        t_x_xxz_xy_xz[i] = t_0_xxxz_xy_xz[i] - rab_x[i] * t_0_xxz_xy_xz[i];

        t_x_xxz_xy_xy[i] = t_0_xxxz_xy_xy[i] - rab_x[i] * t_0_xxz_xy_xy[i];

        t_x_xxz_xy_xx[i] = t_0_xxxz_xy_xx[i] - rab_x[i] * t_0_xxz_xy_xx[i];

        t_x_xxz_xx_zz[i] = t_0_xxxz_xx_zz[i] - rab_x[i] * t_0_xxz_xx_zz[i];

        t_x_xxz_xx_yz[i] = t_0_xxxz_xx_yz[i] - rab_x[i] * t_0_xxz_xx_yz[i];

        t_x_xxz_xx_yy[i] = t_0_xxxz_xx_yy[i] - rab_x[i] * t_0_xxz_xx_yy[i];

        t_x_xxz_xx_xz[i] = t_0_xxxz_xx_xz[i] - rab_x[i] * t_0_xxz_xx_xz[i];

        t_x_xxz_xx_xy[i] = t_0_xxxz_xx_xy[i] - rab_x[i] * t_0_xxz_xx_xy[i];

        t_x_xxz_xx_xx[i] = t_0_xxxz_xx_xx[i] - rab_x[i] * t_0_xxz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxxy_xx_xx, t_0_xxxy_xx_xy, t_0_xxxy_xx_xz, t_0_xxxy_xx_yy,\
                           t_0_xxxy_xx_yz, t_0_xxxy_xx_zz, t_0_xxxy_xy_xx, t_0_xxxy_xy_xy,\
                           t_0_xxxy_xy_xz, t_0_xxxy_xy_yy, t_0_xxxy_xy_yz, t_0_xxxy_xy_zz,\
                           t_0_xxxy_xz_xx, t_0_xxxy_xz_xy, t_0_xxxy_xz_xz, t_0_xxxy_xz_yy,\
                           t_0_xxxy_xz_yz, t_0_xxxy_xz_zz, t_0_xxxy_yy_xx, t_0_xxxy_yy_xy,\
                           t_0_xxxy_yy_xz, t_0_xxxy_yy_yy, t_0_xxxy_yy_yz, t_0_xxxy_yy_zz,\
                           t_0_xxxy_yz_xx, t_0_xxxy_yz_xy, t_0_xxxy_yz_xz, t_0_xxxy_yz_yy,\
                           t_0_xxxy_yz_yz, t_0_xxxy_yz_zz, t_0_xxxy_zz_xx, t_0_xxxy_zz_xy,\
                           t_0_xxxy_zz_xz, t_0_xxxy_zz_yy, t_0_xxxy_zz_yz, t_0_xxxy_zz_zz,\
                           t_0_xxy_xx_xx, t_0_xxy_xx_xy, t_0_xxy_xx_xz, t_0_xxy_xx_yy,\
                           t_0_xxy_xx_yz, t_0_xxy_xx_zz, t_0_xxy_xy_xx, t_0_xxy_xy_xy,\
                           t_0_xxy_xy_xz, t_0_xxy_xy_yy, t_0_xxy_xy_yz, t_0_xxy_xy_zz,\
                           t_0_xxy_xz_xx, t_0_xxy_xz_xy, t_0_xxy_xz_xz, t_0_xxy_xz_yy,\
                           t_0_xxy_xz_yz, t_0_xxy_xz_zz, t_0_xxy_yy_xx, t_0_xxy_yy_xy,\
                           t_0_xxy_yy_xz, t_0_xxy_yy_yy, t_0_xxy_yy_yz, t_0_xxy_yy_zz,\
                           t_0_xxy_yz_xx, t_0_xxy_yz_xy, t_0_xxy_yz_xz, t_0_xxy_yz_yy,\
                           t_0_xxy_yz_yz, t_0_xxy_yz_zz, t_0_xxy_zz_xx, t_0_xxy_zz_xy,\
                           t_0_xxy_zz_xz, t_0_xxy_zz_yy, t_0_xxy_zz_yz, t_0_xxy_zz_zz,\
                           t_x_xxy_xx_xx, t_x_xxy_xx_xy, t_x_xxy_xx_xz, t_x_xxy_xx_yy,\
                           t_x_xxy_xx_yz, t_x_xxy_xx_zz, t_x_xxy_xy_xx, t_x_xxy_xy_xy,\
                           t_x_xxy_xy_xz, t_x_xxy_xy_yy, t_x_xxy_xy_yz, t_x_xxy_xy_zz,\
                           t_x_xxy_xz_xx, t_x_xxy_xz_xy, t_x_xxy_xz_xz, t_x_xxy_xz_yy,\
                           t_x_xxy_xz_yz, t_x_xxy_xz_zz, t_x_xxy_yy_xx, t_x_xxy_yy_xy,\
                           t_x_xxy_yy_xz, t_x_xxy_yy_yy, t_x_xxy_yy_yz, t_x_xxy_yy_zz,\
                           t_x_xxy_yz_xx, t_x_xxy_yz_xy, t_x_xxy_yz_xz, t_x_xxy_yz_yy,\
                           t_x_xxy_yz_yz, t_x_xxy_yz_zz, t_x_xxy_zz_xx, t_x_xxy_zz_xy,\
                           t_x_xxy_zz_xz, t_x_xxy_zz_yy, t_x_xxy_zz_yz, t_x_xxy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxy_zz_zz[i] = t_0_xxxy_zz_zz[i] - rab_x[i] * t_0_xxy_zz_zz[i];

        t_x_xxy_zz_yz[i] = t_0_xxxy_zz_yz[i] - rab_x[i] * t_0_xxy_zz_yz[i];

        t_x_xxy_zz_yy[i] = t_0_xxxy_zz_yy[i] - rab_x[i] * t_0_xxy_zz_yy[i];

        t_x_xxy_zz_xz[i] = t_0_xxxy_zz_xz[i] - rab_x[i] * t_0_xxy_zz_xz[i];

        t_x_xxy_zz_xy[i] = t_0_xxxy_zz_xy[i] - rab_x[i] * t_0_xxy_zz_xy[i];

        t_x_xxy_zz_xx[i] = t_0_xxxy_zz_xx[i] - rab_x[i] * t_0_xxy_zz_xx[i];

        t_x_xxy_yz_zz[i] = t_0_xxxy_yz_zz[i] - rab_x[i] * t_0_xxy_yz_zz[i];

        t_x_xxy_yz_yz[i] = t_0_xxxy_yz_yz[i] - rab_x[i] * t_0_xxy_yz_yz[i];

        t_x_xxy_yz_yy[i] = t_0_xxxy_yz_yy[i] - rab_x[i] * t_0_xxy_yz_yy[i];

        t_x_xxy_yz_xz[i] = t_0_xxxy_yz_xz[i] - rab_x[i] * t_0_xxy_yz_xz[i];

        t_x_xxy_yz_xy[i] = t_0_xxxy_yz_xy[i] - rab_x[i] * t_0_xxy_yz_xy[i];

        t_x_xxy_yz_xx[i] = t_0_xxxy_yz_xx[i] - rab_x[i] * t_0_xxy_yz_xx[i];

        t_x_xxy_yy_zz[i] = t_0_xxxy_yy_zz[i] - rab_x[i] * t_0_xxy_yy_zz[i];

        t_x_xxy_yy_yz[i] = t_0_xxxy_yy_yz[i] - rab_x[i] * t_0_xxy_yy_yz[i];

        t_x_xxy_yy_yy[i] = t_0_xxxy_yy_yy[i] - rab_x[i] * t_0_xxy_yy_yy[i];

        t_x_xxy_yy_xz[i] = t_0_xxxy_yy_xz[i] - rab_x[i] * t_0_xxy_yy_xz[i];

        t_x_xxy_yy_xy[i] = t_0_xxxy_yy_xy[i] - rab_x[i] * t_0_xxy_yy_xy[i];

        t_x_xxy_yy_xx[i] = t_0_xxxy_yy_xx[i] - rab_x[i] * t_0_xxy_yy_xx[i];

        t_x_xxy_xz_zz[i] = t_0_xxxy_xz_zz[i] - rab_x[i] * t_0_xxy_xz_zz[i];

        t_x_xxy_xz_yz[i] = t_0_xxxy_xz_yz[i] - rab_x[i] * t_0_xxy_xz_yz[i];

        t_x_xxy_xz_yy[i] = t_0_xxxy_xz_yy[i] - rab_x[i] * t_0_xxy_xz_yy[i];

        t_x_xxy_xz_xz[i] = t_0_xxxy_xz_xz[i] - rab_x[i] * t_0_xxy_xz_xz[i];

        t_x_xxy_xz_xy[i] = t_0_xxxy_xz_xy[i] - rab_x[i] * t_0_xxy_xz_xy[i];

        t_x_xxy_xz_xx[i] = t_0_xxxy_xz_xx[i] - rab_x[i] * t_0_xxy_xz_xx[i];

        t_x_xxy_xy_zz[i] = t_0_xxxy_xy_zz[i] - rab_x[i] * t_0_xxy_xy_zz[i];

        t_x_xxy_xy_yz[i] = t_0_xxxy_xy_yz[i] - rab_x[i] * t_0_xxy_xy_yz[i];

        t_x_xxy_xy_yy[i] = t_0_xxxy_xy_yy[i] - rab_x[i] * t_0_xxy_xy_yy[i];

        t_x_xxy_xy_xz[i] = t_0_xxxy_xy_xz[i] - rab_x[i] * t_0_xxy_xy_xz[i];

        t_x_xxy_xy_xy[i] = t_0_xxxy_xy_xy[i] - rab_x[i] * t_0_xxy_xy_xy[i];

        t_x_xxy_xy_xx[i] = t_0_xxxy_xy_xx[i] - rab_x[i] * t_0_xxy_xy_xx[i];

        t_x_xxy_xx_zz[i] = t_0_xxxy_xx_zz[i] - rab_x[i] * t_0_xxy_xx_zz[i];

        t_x_xxy_xx_yz[i] = t_0_xxxy_xx_yz[i] - rab_x[i] * t_0_xxy_xx_yz[i];

        t_x_xxy_xx_yy[i] = t_0_xxxy_xx_yy[i] - rab_x[i] * t_0_xxy_xx_yy[i];

        t_x_xxy_xx_xz[i] = t_0_xxxy_xx_xz[i] - rab_x[i] * t_0_xxy_xx_xz[i];

        t_x_xxy_xx_xy[i] = t_0_xxxy_xx_xy[i] - rab_x[i] * t_0_xxy_xx_xy[i];

        t_x_xxy_xx_xx[i] = t_0_xxxy_xx_xx[i] - rab_x[i] * t_0_xxy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxx_xx_xx, t_0_xxx_xx_xy, t_0_xxx_xx_xz, t_0_xxx_xx_yy,\
                           t_0_xxx_xx_yz, t_0_xxx_xx_zz, t_0_xxx_xy_xx, t_0_xxx_xy_xy,\
                           t_0_xxx_xy_xz, t_0_xxx_xy_yy, t_0_xxx_xy_yz, t_0_xxx_xy_zz,\
                           t_0_xxx_xz_xx, t_0_xxx_xz_xy, t_0_xxx_xz_xz, t_0_xxx_xz_yy,\
                           t_0_xxx_xz_yz, t_0_xxx_xz_zz, t_0_xxx_yy_xx, t_0_xxx_yy_xy,\
                           t_0_xxx_yy_xz, t_0_xxx_yy_yy, t_0_xxx_yy_yz, t_0_xxx_yy_zz,\
                           t_0_xxx_yz_xx, t_0_xxx_yz_xy, t_0_xxx_yz_xz, t_0_xxx_yz_yy,\
                           t_0_xxx_yz_yz, t_0_xxx_yz_zz, t_0_xxx_zz_xx, t_0_xxx_zz_xy,\
                           t_0_xxx_zz_xz, t_0_xxx_zz_yy, t_0_xxx_zz_yz, t_0_xxx_zz_zz,\
                           t_0_xxxx_xx_xx, t_0_xxxx_xx_xy, t_0_xxxx_xx_xz, t_0_xxxx_xx_yy,\
                           t_0_xxxx_xx_yz, t_0_xxxx_xx_zz, t_0_xxxx_xy_xx, t_0_xxxx_xy_xy,\
                           t_0_xxxx_xy_xz, t_0_xxxx_xy_yy, t_0_xxxx_xy_yz, t_0_xxxx_xy_zz,\
                           t_0_xxxx_xz_xx, t_0_xxxx_xz_xy, t_0_xxxx_xz_xz, t_0_xxxx_xz_yy,\
                           t_0_xxxx_xz_yz, t_0_xxxx_xz_zz, t_0_xxxx_yy_xx, t_0_xxxx_yy_xy,\
                           t_0_xxxx_yy_xz, t_0_xxxx_yy_yy, t_0_xxxx_yy_yz, t_0_xxxx_yy_zz,\
                           t_0_xxxx_yz_xx, t_0_xxxx_yz_xy, t_0_xxxx_yz_xz, t_0_xxxx_yz_yy,\
                           t_0_xxxx_yz_yz, t_0_xxxx_yz_zz, t_0_xxxx_zz_xx, t_0_xxxx_zz_xy,\
                           t_0_xxxx_zz_xz, t_0_xxxx_zz_yy, t_0_xxxx_zz_yz, t_0_xxxx_zz_zz,\
                           t_x_xxx_xx_xx, t_x_xxx_xx_xy, t_x_xxx_xx_xz, t_x_xxx_xx_yy,\
                           t_x_xxx_xx_yz, t_x_xxx_xx_zz, t_x_xxx_xy_xx, t_x_xxx_xy_xy,\
                           t_x_xxx_xy_xz, t_x_xxx_xy_yy, t_x_xxx_xy_yz, t_x_xxx_xy_zz,\
                           t_x_xxx_xz_xx, t_x_xxx_xz_xy, t_x_xxx_xz_xz, t_x_xxx_xz_yy,\
                           t_x_xxx_xz_yz, t_x_xxx_xz_zz, t_x_xxx_yy_xx, t_x_xxx_yy_xy,\
                           t_x_xxx_yy_xz, t_x_xxx_yy_yy, t_x_xxx_yy_yz, t_x_xxx_yy_zz,\
                           t_x_xxx_yz_xx, t_x_xxx_yz_xy, t_x_xxx_yz_xz, t_x_xxx_yz_yy,\
                           t_x_xxx_yz_yz, t_x_xxx_yz_zz, t_x_xxx_zz_xx, t_x_xxx_zz_xy,\
                           t_x_xxx_zz_xz, t_x_xxx_zz_yy, t_x_xxx_zz_yz, t_x_xxx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxx_zz_zz[i] = t_0_xxxx_zz_zz[i] - rab_x[i] * t_0_xxx_zz_zz[i];

        t_x_xxx_zz_yz[i] = t_0_xxxx_zz_yz[i] - rab_x[i] * t_0_xxx_zz_yz[i];

        t_x_xxx_zz_yy[i] = t_0_xxxx_zz_yy[i] - rab_x[i] * t_0_xxx_zz_yy[i];

        t_x_xxx_zz_xz[i] = t_0_xxxx_zz_xz[i] - rab_x[i] * t_0_xxx_zz_xz[i];

        t_x_xxx_zz_xy[i] = t_0_xxxx_zz_xy[i] - rab_x[i] * t_0_xxx_zz_xy[i];

        t_x_xxx_zz_xx[i] = t_0_xxxx_zz_xx[i] - rab_x[i] * t_0_xxx_zz_xx[i];

        t_x_xxx_yz_zz[i] = t_0_xxxx_yz_zz[i] - rab_x[i] * t_0_xxx_yz_zz[i];

        t_x_xxx_yz_yz[i] = t_0_xxxx_yz_yz[i] - rab_x[i] * t_0_xxx_yz_yz[i];

        t_x_xxx_yz_yy[i] = t_0_xxxx_yz_yy[i] - rab_x[i] * t_0_xxx_yz_yy[i];

        t_x_xxx_yz_xz[i] = t_0_xxxx_yz_xz[i] - rab_x[i] * t_0_xxx_yz_xz[i];

        t_x_xxx_yz_xy[i] = t_0_xxxx_yz_xy[i] - rab_x[i] * t_0_xxx_yz_xy[i];

        t_x_xxx_yz_xx[i] = t_0_xxxx_yz_xx[i] - rab_x[i] * t_0_xxx_yz_xx[i];

        t_x_xxx_yy_zz[i] = t_0_xxxx_yy_zz[i] - rab_x[i] * t_0_xxx_yy_zz[i];

        t_x_xxx_yy_yz[i] = t_0_xxxx_yy_yz[i] - rab_x[i] * t_0_xxx_yy_yz[i];

        t_x_xxx_yy_yy[i] = t_0_xxxx_yy_yy[i] - rab_x[i] * t_0_xxx_yy_yy[i];

        t_x_xxx_yy_xz[i] = t_0_xxxx_yy_xz[i] - rab_x[i] * t_0_xxx_yy_xz[i];

        t_x_xxx_yy_xy[i] = t_0_xxxx_yy_xy[i] - rab_x[i] * t_0_xxx_yy_xy[i];

        t_x_xxx_yy_xx[i] = t_0_xxxx_yy_xx[i] - rab_x[i] * t_0_xxx_yy_xx[i];

        t_x_xxx_xz_zz[i] = t_0_xxxx_xz_zz[i] - rab_x[i] * t_0_xxx_xz_zz[i];

        t_x_xxx_xz_yz[i] = t_0_xxxx_xz_yz[i] - rab_x[i] * t_0_xxx_xz_yz[i];

        t_x_xxx_xz_yy[i] = t_0_xxxx_xz_yy[i] - rab_x[i] * t_0_xxx_xz_yy[i];

        t_x_xxx_xz_xz[i] = t_0_xxxx_xz_xz[i] - rab_x[i] * t_0_xxx_xz_xz[i];

        t_x_xxx_xz_xy[i] = t_0_xxxx_xz_xy[i] - rab_x[i] * t_0_xxx_xz_xy[i];

        t_x_xxx_xz_xx[i] = t_0_xxxx_xz_xx[i] - rab_x[i] * t_0_xxx_xz_xx[i];

        t_x_xxx_xy_zz[i] = t_0_xxxx_xy_zz[i] - rab_x[i] * t_0_xxx_xy_zz[i];

        t_x_xxx_xy_yz[i] = t_0_xxxx_xy_yz[i] - rab_x[i] * t_0_xxx_xy_yz[i];

        t_x_xxx_xy_yy[i] = t_0_xxxx_xy_yy[i] - rab_x[i] * t_0_xxx_xy_yy[i];

        t_x_xxx_xy_xz[i] = t_0_xxxx_xy_xz[i] - rab_x[i] * t_0_xxx_xy_xz[i];

        t_x_xxx_xy_xy[i] = t_0_xxxx_xy_xy[i] - rab_x[i] * t_0_xxx_xy_xy[i];

        t_x_xxx_xy_xx[i] = t_0_xxxx_xy_xx[i] - rab_x[i] * t_0_xxx_xy_xx[i];

        t_x_xxx_xx_zz[i] = t_0_xxxx_xx_zz[i] - rab_x[i] * t_0_xxx_xx_zz[i];

        t_x_xxx_xx_yz[i] = t_0_xxxx_xx_yz[i] - rab_x[i] * t_0_xxx_xx_yz[i];

        t_x_xxx_xx_yy[i] = t_0_xxxx_xx_yy[i] - rab_x[i] * t_0_xxx_xx_yy[i];

        t_x_xxx_xx_xz[i] = t_0_xxxx_xx_xz[i] - rab_x[i] * t_0_xxx_xx_xz[i];

        t_x_xxx_xx_xy[i] = t_0_xxxx_xx_xy[i] - rab_x[i] * t_0_xxx_xx_xy[i];

        t_x_xxx_xx_xx[i] = t_0_xxxx_xx_xx[i] - rab_x[i] * t_0_xxx_xx_xx[i];
    }
}


} // derirec namespace
