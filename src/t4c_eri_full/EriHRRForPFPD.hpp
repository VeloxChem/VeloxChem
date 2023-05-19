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
compHostHRRForPFPD_V0(      BufferHostXY<T>&      intsBufferPFPD,
                      const BufferHostX<int32_t>& intsIndexesPFPD,
                      const BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostXY<T>&      intsBufferSGPD,
                      const BufferHostX<int32_t>& intsIndexesSGPD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PFPD) integral components

    t_z_zzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(0));

    t_z_zzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(1));

    t_z_zzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(2));

    t_z_zzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(3));

    t_z_zzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(4));

    t_z_zzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(5));

    t_z_zzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(6));

    t_z_zzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(7));

    t_z_zzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(8));

    t_z_zzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(9));

    t_z_zzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(10));

    t_z_zzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(11));

    t_z_zzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(12));

    t_z_zzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(13));

    t_z_zzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(14));

    t_z_zzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(15));

    t_z_zzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(16));

    t_z_zzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(17));

    t_z_yzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(18));

    t_z_yzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(19));

    t_z_yzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(20));

    t_z_yzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(21));

    t_z_yzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(22));

    t_z_yzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(23));

    t_z_yzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(24));

    t_z_yzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(25));

    t_z_yzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(26));

    t_z_yzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(27));

    t_z_yzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(28));

    t_z_yzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(29));

    t_z_yzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(30));

    t_z_yzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(31));

    t_z_yzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(32));

    t_z_yzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(33));

    t_z_yzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(34));

    t_z_yzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(35));

    t_z_yyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(36));

    t_z_yyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(37));

    t_z_yyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(38));

    t_z_yyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(39));

    t_z_yyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(40));

    t_z_yyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(41));

    t_z_yyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(42));

    t_z_yyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(43));

    t_z_yyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(44));

    t_z_yyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(45));

    t_z_yyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(46));

    t_z_yyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(47));

    t_z_yyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(48));

    t_z_yyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(49));

    t_z_yyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(50));

    t_z_yyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(51));

    t_z_yyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(52));

    t_z_yyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(53));

    t_z_xzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(54));

    t_z_xzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(55));

    t_z_xzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(56));

    t_z_xzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(57));

    t_z_xzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(58));

    t_z_xzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(59));

    t_z_xzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(60));

    t_z_xzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(61));

    t_z_xzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(62));

    t_z_xzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(63));

    t_z_xzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(64));

    t_z_xzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(65));

    t_z_xzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(66));

    t_z_xzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(67));

    t_z_xzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(68));

    t_z_xzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(69));

    t_z_xzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(70));

    t_z_xzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(71));

    t_z_xyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(72));

    t_z_xyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(73));

    t_z_xyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(74));

    t_z_xyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(75));

    t_z_xyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(76));

    t_z_xyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(77));

    t_z_xyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(78));

    t_z_xyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(79));

    t_z_xyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(80));

    t_z_xyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(81));

    t_z_xyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(82));

    t_z_xyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(83));

    t_z_xyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(84));

    t_z_xyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(85));

    t_z_xyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(86));

    t_z_xyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(87));

    t_z_xyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(88));

    t_z_xyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(89));

    t_z_xxz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(90));

    t_z_xxz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(91));

    t_z_xxz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(92));

    t_z_xxz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(93));

    t_z_xxz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(94));

    t_z_xxz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(95));

    t_z_xxz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(96));

    t_z_xxz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(97));

    t_z_xxz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(98));

    t_z_xxz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(99));

    t_z_xxz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(100));

    t_z_xxz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(101));

    t_z_xxz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(102));

    t_z_xxz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(103));

    t_z_xxz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(104));

    t_z_xxz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(105));

    t_z_xxz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(106));

    t_z_xxz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(107));

    t_y_zzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(108));

    t_y_zzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(109));

    t_y_zzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(110));

    t_y_zzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(111));

    t_y_zzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(112));

    t_y_zzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(113));

    t_y_zzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(114));

    t_y_zzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(115));

    t_y_zzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(116));

    t_y_zzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(117));

    t_y_zzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(118));

    t_y_zzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(119));

    t_y_zzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(120));

    t_y_zzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(121));

    t_y_zzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(122));

    t_y_zzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(123));

    t_y_zzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(124));

    t_y_zzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(125));

    t_y_yzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(126));

    t_y_yzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(127));

    t_y_yzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(128));

    t_y_yzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(129));

    t_y_yzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(130));

    t_y_yzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(131));

    t_y_yzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(132));

    t_y_yzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(133));

    t_y_yzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(134));

    t_y_yzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(135));

    t_y_yzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(136));

    t_y_yzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(137));

    t_y_yzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(138));

    t_y_yzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(139));

    t_y_yzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(140));

    t_y_yzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(141));

    t_y_yzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(142));

    t_y_yzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(143));

    t_y_yyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(144));

    t_y_yyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(145));

    t_y_yyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(146));

    t_y_yyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(147));

    t_y_yyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(148));

    t_y_yyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(149));

    t_y_yyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(150));

    t_y_yyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(151));

    t_y_yyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(152));

    t_y_yyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(153));

    t_y_yyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(154));

    t_y_yyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(155));

    t_y_yyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(156));

    t_y_yyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(157));

    t_y_yyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(158));

    t_y_yyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(159));

    t_y_yyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(160));

    t_y_yyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(161));

    t_y_yyy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(162));

    t_y_yyy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(163));

    t_y_yyy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(164));

    t_y_yyy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(165));

    t_y_yyy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(166));

    t_y_yyy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(167));

    t_y_yyy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(168));

    t_y_yyy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(169));

    t_y_yyy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(170));

    t_y_yyy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(171));

    t_y_yyy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(172));

    t_y_yyy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(173));

    t_y_yyy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(174));

    t_y_yyy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(175));

    t_y_yyy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(176));

    t_y_yyy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(177));

    t_y_yyy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(178));

    t_y_yyy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(179));

    t_y_xzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(180));

    t_y_xzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(181));

    t_y_xzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(182));

    t_y_xzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(183));

    t_y_xzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(184));

    t_y_xzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(185));

    t_y_xzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(186));

    t_y_xzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(187));

    t_y_xzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(188));

    t_y_xzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(189));

    t_y_xzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(190));

    t_y_xzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(191));

    t_y_xzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(192));

    t_y_xzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(193));

    t_y_xzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(194));

    t_y_xzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(195));

    t_y_xzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(196));

    t_y_xzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(197));

    t_y_xyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(198));

    t_y_xyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(199));

    t_y_xyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(200));

    t_y_xyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(201));

    t_y_xyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(202));

    t_y_xyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(203));

    t_y_xyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(204));

    t_y_xyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(205));

    t_y_xyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(206));

    t_y_xyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(207));

    t_y_xyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(208));

    t_y_xyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(209));

    t_y_xyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(210));

    t_y_xyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(211));

    t_y_xyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(212));

    t_y_xyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(213));

    t_y_xyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(214));

    t_y_xyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(215));

    t_y_xyy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(216));

    t_y_xyy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(217));

    t_y_xyy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(218));

    t_y_xyy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(219));

    t_y_xyy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(220));

    t_y_xyy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(221));

    t_y_xyy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(222));

    t_y_xyy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(223));

    t_y_xyy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(224));

    t_y_xyy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(225));

    t_y_xyy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(226));

    t_y_xyy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(227));

    t_y_xyy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(228));

    t_y_xyy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(229));

    t_y_xyy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(230));

    t_y_xyy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(231));

    t_y_xyy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(232));

    t_y_xyy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(233));

    t_y_xxz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(234));

    t_y_xxz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(235));

    t_y_xxz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(236));

    t_y_xxz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(237));

    t_y_xxz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(238));

    t_y_xxz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(239));

    t_y_xxz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(240));

    t_y_xxz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(241));

    t_y_xxz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(242));

    t_y_xxz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(243));

    t_y_xxz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(244));

    t_y_xxz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(245));

    t_y_xxz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(246));

    t_y_xxz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(247));

    t_y_xxz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(248));

    t_y_xxz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(249));

    t_y_xxz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(250));

    t_y_xxz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(251));

    t_y_xxy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(252));

    t_y_xxy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(253));

    t_y_xxy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(254));

    t_y_xxy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(255));

    t_y_xxy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(256));

    t_y_xxy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(257));

    t_y_xxy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(258));

    t_y_xxy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(259));

    t_y_xxy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(260));

    t_y_xxy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(261));

    t_y_xxy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(262));

    t_y_xxy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(263));

    t_y_xxy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(264));

    t_y_xxy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(265));

    t_y_xxy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(266));

    t_y_xxy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(267));

    t_y_xxy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(268));

    t_y_xxy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(269));

    t_x_zzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(270));

    t_x_zzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(271));

    t_x_zzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(272));

    t_x_zzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(273));

    t_x_zzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(274));

    t_x_zzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(275));

    t_x_zzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(276));

    t_x_zzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(277));

    t_x_zzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(278));

    t_x_zzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(279));

    t_x_zzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(280));

    t_x_zzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(281));

    t_x_zzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(282));

    t_x_zzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(283));

    t_x_zzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(284));

    t_x_zzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(285));

    t_x_zzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(286));

    t_x_zzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(287));

    t_x_yzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(288));

    t_x_yzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(289));

    t_x_yzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(290));

    t_x_yzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(291));

    t_x_yzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(292));

    t_x_yzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(293));

    t_x_yzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(294));

    t_x_yzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(295));

    t_x_yzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(296));

    t_x_yzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(297));

    t_x_yzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(298));

    t_x_yzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(299));

    t_x_yzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(300));

    t_x_yzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(301));

    t_x_yzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(302));

    t_x_yzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(303));

    t_x_yzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(304));

    t_x_yzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(305));

    t_x_yyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(306));

    t_x_yyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(307));

    t_x_yyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(308));

    t_x_yyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(309));

    t_x_yyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(310));

    t_x_yyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(311));

    t_x_yyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(312));

    t_x_yyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(313));

    t_x_yyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(314));

    t_x_yyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(315));

    t_x_yyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(316));

    t_x_yyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(317));

    t_x_yyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(318));

    t_x_yyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(319));

    t_x_yyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(320));

    t_x_yyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(321));

    t_x_yyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(322));

    t_x_yyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(323));

    t_x_yyy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(324));

    t_x_yyy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(325));

    t_x_yyy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(326));

    t_x_yyy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(327));

    t_x_yyy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(328));

    t_x_yyy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(329));

    t_x_yyy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(330));

    t_x_yyy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(331));

    t_x_yyy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(332));

    t_x_yyy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(333));

    t_x_yyy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(334));

    t_x_yyy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(335));

    t_x_yyy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(336));

    t_x_yyy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(337));

    t_x_yyy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(338));

    t_x_yyy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(339));

    t_x_yyy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(340));

    t_x_yyy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(341));

    t_x_xzz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(342));

    t_x_xzz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(343));

    t_x_xzz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(344));

    t_x_xzz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(345));

    t_x_xzz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(346));

    t_x_xzz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(347));

    t_x_xzz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(348));

    t_x_xzz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(349));

    t_x_xzz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(350));

    t_x_xzz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(351));

    t_x_xzz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(352));

    t_x_xzz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(353));

    t_x_xzz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(354));

    t_x_xzz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(355));

    t_x_xzz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(356));

    t_x_xzz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(357));

    t_x_xzz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(358));

    t_x_xzz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(359));

    t_x_xyz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(360));

    t_x_xyz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(361));

    t_x_xyz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(362));

    t_x_xyz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(363));

    t_x_xyz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(364));

    t_x_xyz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(365));

    t_x_xyz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(366));

    t_x_xyz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(367));

    t_x_xyz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(368));

    t_x_xyz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(369));

    t_x_xyz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(370));

    t_x_xyz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(371));

    t_x_xyz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(372));

    t_x_xyz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(373));

    t_x_xyz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(374));

    t_x_xyz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(375));

    t_x_xyz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(376));

    t_x_xyz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(377));

    t_x_xyy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(378));

    t_x_xyy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(379));

    t_x_xyy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(380));

    t_x_xyy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(381));

    t_x_xyy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(382));

    t_x_xyy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(383));

    t_x_xyy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(384));

    t_x_xyy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(385));

    t_x_xyy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(386));

    t_x_xyy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(387));

    t_x_xyy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(388));

    t_x_xyy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(389));

    t_x_xyy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(390));

    t_x_xyy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(391));

    t_x_xyy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(392));

    t_x_xyy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(393));

    t_x_xyy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(394));

    t_x_xyy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(395));

    t_x_xxz_z_zz = intsBufferPFPD.data(intsIndexesPFPD(396));

    t_x_xxz_z_yz = intsBufferPFPD.data(intsIndexesPFPD(397));

    t_x_xxz_z_yy = intsBufferPFPD.data(intsIndexesPFPD(398));

    t_x_xxz_z_xz = intsBufferPFPD.data(intsIndexesPFPD(399));

    t_x_xxz_z_xy = intsBufferPFPD.data(intsIndexesPFPD(400));

    t_x_xxz_z_xx = intsBufferPFPD.data(intsIndexesPFPD(401));

    t_x_xxz_y_zz = intsBufferPFPD.data(intsIndexesPFPD(402));

    t_x_xxz_y_yz = intsBufferPFPD.data(intsIndexesPFPD(403));

    t_x_xxz_y_yy = intsBufferPFPD.data(intsIndexesPFPD(404));

    t_x_xxz_y_xz = intsBufferPFPD.data(intsIndexesPFPD(405));

    t_x_xxz_y_xy = intsBufferPFPD.data(intsIndexesPFPD(406));

    t_x_xxz_y_xx = intsBufferPFPD.data(intsIndexesPFPD(407));

    t_x_xxz_x_zz = intsBufferPFPD.data(intsIndexesPFPD(408));

    t_x_xxz_x_yz = intsBufferPFPD.data(intsIndexesPFPD(409));

    t_x_xxz_x_yy = intsBufferPFPD.data(intsIndexesPFPD(410));

    t_x_xxz_x_xz = intsBufferPFPD.data(intsIndexesPFPD(411));

    t_x_xxz_x_xy = intsBufferPFPD.data(intsIndexesPFPD(412));

    t_x_xxz_x_xx = intsBufferPFPD.data(intsIndexesPFPD(413));

    t_x_xxy_z_zz = intsBufferPFPD.data(intsIndexesPFPD(414));

    t_x_xxy_z_yz = intsBufferPFPD.data(intsIndexesPFPD(415));

    t_x_xxy_z_yy = intsBufferPFPD.data(intsIndexesPFPD(416));

    t_x_xxy_z_xz = intsBufferPFPD.data(intsIndexesPFPD(417));

    t_x_xxy_z_xy = intsBufferPFPD.data(intsIndexesPFPD(418));

    t_x_xxy_z_xx = intsBufferPFPD.data(intsIndexesPFPD(419));

    t_x_xxy_y_zz = intsBufferPFPD.data(intsIndexesPFPD(420));

    t_x_xxy_y_yz = intsBufferPFPD.data(intsIndexesPFPD(421));

    t_x_xxy_y_yy = intsBufferPFPD.data(intsIndexesPFPD(422));

    t_x_xxy_y_xz = intsBufferPFPD.data(intsIndexesPFPD(423));

    t_x_xxy_y_xy = intsBufferPFPD.data(intsIndexesPFPD(424));

    t_x_xxy_y_xx = intsBufferPFPD.data(intsIndexesPFPD(425));

    t_x_xxy_x_zz = intsBufferPFPD.data(intsIndexesPFPD(426));

    t_x_xxy_x_yz = intsBufferPFPD.data(intsIndexesPFPD(427));

    t_x_xxy_x_yy = intsBufferPFPD.data(intsIndexesPFPD(428));

    t_x_xxy_x_xz = intsBufferPFPD.data(intsIndexesPFPD(429));

    t_x_xxy_x_xy = intsBufferPFPD.data(intsIndexesPFPD(430));

    t_x_xxy_x_xx = intsBufferPFPD.data(intsIndexesPFPD(431));

    t_x_xxx_z_zz = intsBufferPFPD.data(intsIndexesPFPD(432));

    t_x_xxx_z_yz = intsBufferPFPD.data(intsIndexesPFPD(433));

    t_x_xxx_z_yy = intsBufferPFPD.data(intsIndexesPFPD(434));

    t_x_xxx_z_xz = intsBufferPFPD.data(intsIndexesPFPD(435));

    t_x_xxx_z_xy = intsBufferPFPD.data(intsIndexesPFPD(436));

    t_x_xxx_z_xx = intsBufferPFPD.data(intsIndexesPFPD(437));

    t_x_xxx_y_zz = intsBufferPFPD.data(intsIndexesPFPD(438));

    t_x_xxx_y_yz = intsBufferPFPD.data(intsIndexesPFPD(439));

    t_x_xxx_y_yy = intsBufferPFPD.data(intsIndexesPFPD(440));

    t_x_xxx_y_xz = intsBufferPFPD.data(intsIndexesPFPD(441));

    t_x_xxx_y_xy = intsBufferPFPD.data(intsIndexesPFPD(442));

    t_x_xxx_y_xx = intsBufferPFPD.data(intsIndexesPFPD(443));

    t_x_xxx_x_zz = intsBufferPFPD.data(intsIndexesPFPD(444));

    t_x_xxx_x_yz = intsBufferPFPD.data(intsIndexesPFPD(445));

    t_x_xxx_x_yy = intsBufferPFPD.data(intsIndexesPFPD(446));

    t_x_xxx_x_xz = intsBufferPFPD.data(intsIndexesPFPD(447));

    t_x_xxx_x_xy = intsBufferPFPD.data(intsIndexesPFPD(448));

    t_x_xxx_x_xx = intsBufferPFPD.data(intsIndexesPFPD(449));

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

    // set up (SGPD) integral components

    t_0_zzzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(0));

    t_0_zzzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(1));

    t_0_zzzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(2));

    t_0_zzzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(3));

    t_0_zzzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(4));

    t_0_zzzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(5));

    t_0_zzzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(6));

    t_0_zzzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(7));

    t_0_zzzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(8));

    t_0_zzzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(9));

    t_0_zzzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(10));

    t_0_zzzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(11));

    t_0_zzzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(12));

    t_0_zzzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(13));

    t_0_zzzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(14));

    t_0_zzzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(15));

    t_0_zzzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(16));

    t_0_zzzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(17));

    t_0_yzzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(18));

    t_0_yzzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(19));

    t_0_yzzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(20));

    t_0_yzzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(21));

    t_0_yzzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(22));

    t_0_yzzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(23));

    t_0_yzzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(24));

    t_0_yzzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(25));

    t_0_yzzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(26));

    t_0_yzzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(27));

    t_0_yzzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(28));

    t_0_yzzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(29));

    t_0_yzzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(30));

    t_0_yzzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(31));

    t_0_yzzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(32));

    t_0_yzzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(33));

    t_0_yzzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(34));

    t_0_yzzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(35));

    t_0_yyzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(36));

    t_0_yyzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(37));

    t_0_yyzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(38));

    t_0_yyzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(39));

    t_0_yyzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(40));

    t_0_yyzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(41));

    t_0_yyzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(42));

    t_0_yyzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(43));

    t_0_yyzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(44));

    t_0_yyzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(45));

    t_0_yyzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(46));

    t_0_yyzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(47));

    t_0_yyzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(48));

    t_0_yyzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(49));

    t_0_yyzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(50));

    t_0_yyzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(51));

    t_0_yyzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(52));

    t_0_yyzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(53));

    t_0_yyyz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(54));

    t_0_yyyz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(55));

    t_0_yyyz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(56));

    t_0_yyyz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(57));

    t_0_yyyz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(58));

    t_0_yyyz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(59));

    t_0_yyyz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(60));

    t_0_yyyz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(61));

    t_0_yyyz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(62));

    t_0_yyyz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(63));

    t_0_yyyz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(64));

    t_0_yyyz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(65));

    t_0_yyyz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(66));

    t_0_yyyz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(67));

    t_0_yyyz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(68));

    t_0_yyyz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(69));

    t_0_yyyz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(70));

    t_0_yyyz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(71));

    t_0_yyyy_z_zz = intsBufferSGPD.data(intsIndexesSGPD(72));

    t_0_yyyy_z_yz = intsBufferSGPD.data(intsIndexesSGPD(73));

    t_0_yyyy_z_yy = intsBufferSGPD.data(intsIndexesSGPD(74));

    t_0_yyyy_z_xz = intsBufferSGPD.data(intsIndexesSGPD(75));

    t_0_yyyy_z_xy = intsBufferSGPD.data(intsIndexesSGPD(76));

    t_0_yyyy_z_xx = intsBufferSGPD.data(intsIndexesSGPD(77));

    t_0_yyyy_y_zz = intsBufferSGPD.data(intsIndexesSGPD(78));

    t_0_yyyy_y_yz = intsBufferSGPD.data(intsIndexesSGPD(79));

    t_0_yyyy_y_yy = intsBufferSGPD.data(intsIndexesSGPD(80));

    t_0_yyyy_y_xz = intsBufferSGPD.data(intsIndexesSGPD(81));

    t_0_yyyy_y_xy = intsBufferSGPD.data(intsIndexesSGPD(82));

    t_0_yyyy_y_xx = intsBufferSGPD.data(intsIndexesSGPD(83));

    t_0_yyyy_x_zz = intsBufferSGPD.data(intsIndexesSGPD(84));

    t_0_yyyy_x_yz = intsBufferSGPD.data(intsIndexesSGPD(85));

    t_0_yyyy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(86));

    t_0_yyyy_x_xz = intsBufferSGPD.data(intsIndexesSGPD(87));

    t_0_yyyy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(88));

    t_0_yyyy_x_xx = intsBufferSGPD.data(intsIndexesSGPD(89));

    t_0_xzzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(90));

    t_0_xzzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(91));

    t_0_xzzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(92));

    t_0_xzzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(93));

    t_0_xzzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(94));

    t_0_xzzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(95));

    t_0_xzzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(96));

    t_0_xzzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(97));

    t_0_xzzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(98));

    t_0_xzzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(99));

    t_0_xzzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(100));

    t_0_xzzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(101));

    t_0_xzzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(102));

    t_0_xzzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(103));

    t_0_xzzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(104));

    t_0_xzzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(105));

    t_0_xzzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(106));

    t_0_xzzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(107));

    t_0_xyzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(108));

    t_0_xyzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(109));

    t_0_xyzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(110));

    t_0_xyzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(111));

    t_0_xyzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(112));

    t_0_xyzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(113));

    t_0_xyzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(114));

    t_0_xyzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(115));

    t_0_xyzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(116));

    t_0_xyzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(117));

    t_0_xyzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(118));

    t_0_xyzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(119));

    t_0_xyzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(120));

    t_0_xyzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(121));

    t_0_xyzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(122));

    t_0_xyzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(123));

    t_0_xyzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(124));

    t_0_xyzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(125));

    t_0_xyyz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(126));

    t_0_xyyz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(127));

    t_0_xyyz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(128));

    t_0_xyyz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(129));

    t_0_xyyz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(130));

    t_0_xyyz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(131));

    t_0_xyyz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(132));

    t_0_xyyz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(133));

    t_0_xyyz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(134));

    t_0_xyyz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(135));

    t_0_xyyz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(136));

    t_0_xyyz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(137));

    t_0_xyyz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(138));

    t_0_xyyz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(139));

    t_0_xyyz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(140));

    t_0_xyyz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(141));

    t_0_xyyz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(142));

    t_0_xyyz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(143));

    t_0_xyyy_z_zz = intsBufferSGPD.data(intsIndexesSGPD(144));

    t_0_xyyy_z_yz = intsBufferSGPD.data(intsIndexesSGPD(145));

    t_0_xyyy_z_yy = intsBufferSGPD.data(intsIndexesSGPD(146));

    t_0_xyyy_z_xz = intsBufferSGPD.data(intsIndexesSGPD(147));

    t_0_xyyy_z_xy = intsBufferSGPD.data(intsIndexesSGPD(148));

    t_0_xyyy_z_xx = intsBufferSGPD.data(intsIndexesSGPD(149));

    t_0_xyyy_y_zz = intsBufferSGPD.data(intsIndexesSGPD(150));

    t_0_xyyy_y_yz = intsBufferSGPD.data(intsIndexesSGPD(151));

    t_0_xyyy_y_yy = intsBufferSGPD.data(intsIndexesSGPD(152));

    t_0_xyyy_y_xz = intsBufferSGPD.data(intsIndexesSGPD(153));

    t_0_xyyy_y_xy = intsBufferSGPD.data(intsIndexesSGPD(154));

    t_0_xyyy_y_xx = intsBufferSGPD.data(intsIndexesSGPD(155));

    t_0_xyyy_x_zz = intsBufferSGPD.data(intsIndexesSGPD(156));

    t_0_xyyy_x_yz = intsBufferSGPD.data(intsIndexesSGPD(157));

    t_0_xyyy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(158));

    t_0_xyyy_x_xz = intsBufferSGPD.data(intsIndexesSGPD(159));

    t_0_xyyy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(160));

    t_0_xyyy_x_xx = intsBufferSGPD.data(intsIndexesSGPD(161));

    t_0_xxzz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(162));

    t_0_xxzz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(163));

    t_0_xxzz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(164));

    t_0_xxzz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(165));

    t_0_xxzz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(166));

    t_0_xxzz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(167));

    t_0_xxzz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(168));

    t_0_xxzz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(169));

    t_0_xxzz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(170));

    t_0_xxzz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(171));

    t_0_xxzz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(172));

    t_0_xxzz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(173));

    t_0_xxzz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(174));

    t_0_xxzz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(175));

    t_0_xxzz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(176));

    t_0_xxzz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(177));

    t_0_xxzz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(178));

    t_0_xxzz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(179));

    t_0_xxyz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(180));

    t_0_xxyz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(181));

    t_0_xxyz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(182));

    t_0_xxyz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(183));

    t_0_xxyz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(184));

    t_0_xxyz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(185));

    t_0_xxyz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(186));

    t_0_xxyz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(187));

    t_0_xxyz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(188));

    t_0_xxyz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(189));

    t_0_xxyz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(190));

    t_0_xxyz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(191));

    t_0_xxyz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(192));

    t_0_xxyz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(193));

    t_0_xxyz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(194));

    t_0_xxyz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(195));

    t_0_xxyz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(196));

    t_0_xxyz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(197));

    t_0_xxyy_z_zz = intsBufferSGPD.data(intsIndexesSGPD(198));

    t_0_xxyy_z_yz = intsBufferSGPD.data(intsIndexesSGPD(199));

    t_0_xxyy_z_yy = intsBufferSGPD.data(intsIndexesSGPD(200));

    t_0_xxyy_z_xz = intsBufferSGPD.data(intsIndexesSGPD(201));

    t_0_xxyy_z_xy = intsBufferSGPD.data(intsIndexesSGPD(202));

    t_0_xxyy_z_xx = intsBufferSGPD.data(intsIndexesSGPD(203));

    t_0_xxyy_y_zz = intsBufferSGPD.data(intsIndexesSGPD(204));

    t_0_xxyy_y_yz = intsBufferSGPD.data(intsIndexesSGPD(205));

    t_0_xxyy_y_yy = intsBufferSGPD.data(intsIndexesSGPD(206));

    t_0_xxyy_y_xz = intsBufferSGPD.data(intsIndexesSGPD(207));

    t_0_xxyy_y_xy = intsBufferSGPD.data(intsIndexesSGPD(208));

    t_0_xxyy_y_xx = intsBufferSGPD.data(intsIndexesSGPD(209));

    t_0_xxyy_x_zz = intsBufferSGPD.data(intsIndexesSGPD(210));

    t_0_xxyy_x_yz = intsBufferSGPD.data(intsIndexesSGPD(211));

    t_0_xxyy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(212));

    t_0_xxyy_x_xz = intsBufferSGPD.data(intsIndexesSGPD(213));

    t_0_xxyy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(214));

    t_0_xxyy_x_xx = intsBufferSGPD.data(intsIndexesSGPD(215));

    t_0_xxxz_z_zz = intsBufferSGPD.data(intsIndexesSGPD(216));

    t_0_xxxz_z_yz = intsBufferSGPD.data(intsIndexesSGPD(217));

    t_0_xxxz_z_yy = intsBufferSGPD.data(intsIndexesSGPD(218));

    t_0_xxxz_z_xz = intsBufferSGPD.data(intsIndexesSGPD(219));

    t_0_xxxz_z_xy = intsBufferSGPD.data(intsIndexesSGPD(220));

    t_0_xxxz_z_xx = intsBufferSGPD.data(intsIndexesSGPD(221));

    t_0_xxxz_y_zz = intsBufferSGPD.data(intsIndexesSGPD(222));

    t_0_xxxz_y_yz = intsBufferSGPD.data(intsIndexesSGPD(223));

    t_0_xxxz_y_yy = intsBufferSGPD.data(intsIndexesSGPD(224));

    t_0_xxxz_y_xz = intsBufferSGPD.data(intsIndexesSGPD(225));

    t_0_xxxz_y_xy = intsBufferSGPD.data(intsIndexesSGPD(226));

    t_0_xxxz_y_xx = intsBufferSGPD.data(intsIndexesSGPD(227));

    t_0_xxxz_x_zz = intsBufferSGPD.data(intsIndexesSGPD(228));

    t_0_xxxz_x_yz = intsBufferSGPD.data(intsIndexesSGPD(229));

    t_0_xxxz_x_yy = intsBufferSGPD.data(intsIndexesSGPD(230));

    t_0_xxxz_x_xz = intsBufferSGPD.data(intsIndexesSGPD(231));

    t_0_xxxz_x_xy = intsBufferSGPD.data(intsIndexesSGPD(232));

    t_0_xxxz_x_xx = intsBufferSGPD.data(intsIndexesSGPD(233));

    t_0_xxxy_z_zz = intsBufferSGPD.data(intsIndexesSGPD(234));

    t_0_xxxy_z_yz = intsBufferSGPD.data(intsIndexesSGPD(235));

    t_0_xxxy_z_yy = intsBufferSGPD.data(intsIndexesSGPD(236));

    t_0_xxxy_z_xz = intsBufferSGPD.data(intsIndexesSGPD(237));

    t_0_xxxy_z_xy = intsBufferSGPD.data(intsIndexesSGPD(238));

    t_0_xxxy_z_xx = intsBufferSGPD.data(intsIndexesSGPD(239));

    t_0_xxxy_y_zz = intsBufferSGPD.data(intsIndexesSGPD(240));

    t_0_xxxy_y_yz = intsBufferSGPD.data(intsIndexesSGPD(241));

    t_0_xxxy_y_yy = intsBufferSGPD.data(intsIndexesSGPD(242));

    t_0_xxxy_y_xz = intsBufferSGPD.data(intsIndexesSGPD(243));

    t_0_xxxy_y_xy = intsBufferSGPD.data(intsIndexesSGPD(244));

    t_0_xxxy_y_xx = intsBufferSGPD.data(intsIndexesSGPD(245));

    t_0_xxxy_x_zz = intsBufferSGPD.data(intsIndexesSGPD(246));

    t_0_xxxy_x_yz = intsBufferSGPD.data(intsIndexesSGPD(247));

    t_0_xxxy_x_yy = intsBufferSGPD.data(intsIndexesSGPD(248));

    t_0_xxxy_x_xz = intsBufferSGPD.data(intsIndexesSGPD(249));

    t_0_xxxy_x_xy = intsBufferSGPD.data(intsIndexesSGPD(250));

    t_0_xxxy_x_xx = intsBufferSGPD.data(intsIndexesSGPD(251));

    t_0_xxxx_z_zz = intsBufferSGPD.data(intsIndexesSGPD(252));

    t_0_xxxx_z_yz = intsBufferSGPD.data(intsIndexesSGPD(253));

    t_0_xxxx_z_yy = intsBufferSGPD.data(intsIndexesSGPD(254));

    t_0_xxxx_z_xz = intsBufferSGPD.data(intsIndexesSGPD(255));

    t_0_xxxx_z_xy = intsBufferSGPD.data(intsIndexesSGPD(256));

    t_0_xxxx_z_xx = intsBufferSGPD.data(intsIndexesSGPD(257));

    t_0_xxxx_y_zz = intsBufferSGPD.data(intsIndexesSGPD(258));

    t_0_xxxx_y_yz = intsBufferSGPD.data(intsIndexesSGPD(259));

    t_0_xxxx_y_yy = intsBufferSGPD.data(intsIndexesSGPD(260));

    t_0_xxxx_y_xz = intsBufferSGPD.data(intsIndexesSGPD(261));

    t_0_xxxx_y_xy = intsBufferSGPD.data(intsIndexesSGPD(262));

    t_0_xxxx_y_xx = intsBufferSGPD.data(intsIndexesSGPD(263));

    t_0_xxxx_x_zz = intsBufferSGPD.data(intsIndexesSGPD(264));

    t_0_xxxx_x_yz = intsBufferSGPD.data(intsIndexesSGPD(265));

    t_0_xxxx_x_yy = intsBufferSGPD.data(intsIndexesSGPD(266));

    t_0_xxxx_x_xz = intsBufferSGPD.data(intsIndexesSGPD(267));

    t_0_xxxx_x_xy = intsBufferSGPD.data(intsIndexesSGPD(268));

    t_0_xxxx_x_xx = intsBufferSGPD.data(intsIndexesSGPD(269));

    #pragma omp simd align(rab_z, t_0_yzz_x_xx, t_0_yzz_x_xy, t_0_yzz_x_xz, t_0_yzz_x_yy,\
                           t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_xx, t_0_yzz_y_xy, t_0_yzz_y_xz,\
                           t_0_yzz_y_yy, t_0_yzz_y_yz, t_0_yzz_y_zz, t_0_yzz_z_xx, t_0_yzz_z_xy,\
                           t_0_yzz_z_xz, t_0_yzz_z_yy, t_0_yzz_z_yz, t_0_yzz_z_zz, t_0_yzzz_x_xx,\
                           t_0_yzzz_x_xy, t_0_yzzz_x_xz, t_0_yzzz_x_yy, t_0_yzzz_x_yz,\
                           t_0_yzzz_x_zz, t_0_yzzz_y_xx, t_0_yzzz_y_xy, t_0_yzzz_y_xz,\
                           t_0_yzzz_y_yy, t_0_yzzz_y_yz, t_0_yzzz_y_zz, t_0_yzzz_z_xx,\
                           t_0_yzzz_z_xy, t_0_yzzz_z_xz, t_0_yzzz_z_yy, t_0_yzzz_z_yz,\
                           t_0_yzzz_z_zz, t_0_zzz_x_xx, t_0_zzz_x_xy, t_0_zzz_x_xz,\
                           t_0_zzz_x_yy, t_0_zzz_x_yz, t_0_zzz_x_zz, t_0_zzz_y_xx, t_0_zzz_y_xy,\
                           t_0_zzz_y_xz, t_0_zzz_y_yy, t_0_zzz_y_yz, t_0_zzz_y_zz, t_0_zzz_z_xx,\
                           t_0_zzz_z_xy, t_0_zzz_z_xz, t_0_zzz_z_yy, t_0_zzz_z_yz, t_0_zzz_z_zz,\
                           t_0_zzzz_x_xx, t_0_zzzz_x_xy, t_0_zzzz_x_xz, t_0_zzzz_x_yy,\
                           t_0_zzzz_x_yz, t_0_zzzz_x_zz, t_0_zzzz_y_xx, t_0_zzzz_y_xy,\
                           t_0_zzzz_y_xz, t_0_zzzz_y_yy, t_0_zzzz_y_yz, t_0_zzzz_y_zz,\
                           t_0_zzzz_z_xx, t_0_zzzz_z_xy, t_0_zzzz_z_xz, t_0_zzzz_z_yy,\
                           t_0_zzzz_z_yz, t_0_zzzz_z_zz, t_z_yzz_x_xx, t_z_yzz_x_xy,\
                           t_z_yzz_x_xz, t_z_yzz_x_yy, t_z_yzz_x_yz, t_z_yzz_x_zz, t_z_yzz_y_xx,\
                           t_z_yzz_y_xy, t_z_yzz_y_xz, t_z_yzz_y_yy, t_z_yzz_y_yz, t_z_yzz_y_zz,\
                           t_z_yzz_z_xx, t_z_yzz_z_xy, t_z_yzz_z_xz, t_z_yzz_z_yy, t_z_yzz_z_yz,\
                           t_z_yzz_z_zz, t_z_zzz_x_xx, t_z_zzz_x_xy, t_z_zzz_x_xz, t_z_zzz_x_yy,\
                           t_z_zzz_x_yz, t_z_zzz_x_zz, t_z_zzz_y_xx, t_z_zzz_y_xy, t_z_zzz_y_xz,\
                           t_z_zzz_y_yy, t_z_zzz_y_yz, t_z_zzz_y_zz, t_z_zzz_z_xx, t_z_zzz_z_xy,\
                           t_z_zzz_z_xz, t_z_zzz_z_yy, t_z_zzz_z_yz, t_z_zzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_z_zz[i] = t_0_zzzz_z_zz[i] - rab_z[i] * t_0_zzz_z_zz[i];

        t_z_zzz_z_yz[i] = t_0_zzzz_z_yz[i] - rab_z[i] * t_0_zzz_z_yz[i];

        t_z_zzz_z_yy[i] = t_0_zzzz_z_yy[i] - rab_z[i] * t_0_zzz_z_yy[i];

        t_z_zzz_z_xz[i] = t_0_zzzz_z_xz[i] - rab_z[i] * t_0_zzz_z_xz[i];

        t_z_zzz_z_xy[i] = t_0_zzzz_z_xy[i] - rab_z[i] * t_0_zzz_z_xy[i];

        t_z_zzz_z_xx[i] = t_0_zzzz_z_xx[i] - rab_z[i] * t_0_zzz_z_xx[i];

        t_z_zzz_y_zz[i] = t_0_zzzz_y_zz[i] - rab_z[i] * t_0_zzz_y_zz[i];

        t_z_zzz_y_yz[i] = t_0_zzzz_y_yz[i] - rab_z[i] * t_0_zzz_y_yz[i];

        t_z_zzz_y_yy[i] = t_0_zzzz_y_yy[i] - rab_z[i] * t_0_zzz_y_yy[i];

        t_z_zzz_y_xz[i] = t_0_zzzz_y_xz[i] - rab_z[i] * t_0_zzz_y_xz[i];

        t_z_zzz_y_xy[i] = t_0_zzzz_y_xy[i] - rab_z[i] * t_0_zzz_y_xy[i];

        t_z_zzz_y_xx[i] = t_0_zzzz_y_xx[i] - rab_z[i] * t_0_zzz_y_xx[i];

        t_z_zzz_x_zz[i] = t_0_zzzz_x_zz[i] - rab_z[i] * t_0_zzz_x_zz[i];

        t_z_zzz_x_yz[i] = t_0_zzzz_x_yz[i] - rab_z[i] * t_0_zzz_x_yz[i];

        t_z_zzz_x_yy[i] = t_0_zzzz_x_yy[i] - rab_z[i] * t_0_zzz_x_yy[i];

        t_z_zzz_x_xz[i] = t_0_zzzz_x_xz[i] - rab_z[i] * t_0_zzz_x_xz[i];

        t_z_zzz_x_xy[i] = t_0_zzzz_x_xy[i] - rab_z[i] * t_0_zzz_x_xy[i];

        t_z_zzz_x_xx[i] = t_0_zzzz_x_xx[i] - rab_z[i] * t_0_zzz_x_xx[i];

        t_z_yzz_z_zz[i] = t_0_yzzz_z_zz[i] - rab_z[i] * t_0_yzz_z_zz[i];

        t_z_yzz_z_yz[i] = t_0_yzzz_z_yz[i] - rab_z[i] * t_0_yzz_z_yz[i];

        t_z_yzz_z_yy[i] = t_0_yzzz_z_yy[i] - rab_z[i] * t_0_yzz_z_yy[i];

        t_z_yzz_z_xz[i] = t_0_yzzz_z_xz[i] - rab_z[i] * t_0_yzz_z_xz[i];

        t_z_yzz_z_xy[i] = t_0_yzzz_z_xy[i] - rab_z[i] * t_0_yzz_z_xy[i];

        t_z_yzz_z_xx[i] = t_0_yzzz_z_xx[i] - rab_z[i] * t_0_yzz_z_xx[i];

        t_z_yzz_y_zz[i] = t_0_yzzz_y_zz[i] - rab_z[i] * t_0_yzz_y_zz[i];

        t_z_yzz_y_yz[i] = t_0_yzzz_y_yz[i] - rab_z[i] * t_0_yzz_y_yz[i];

        t_z_yzz_y_yy[i] = t_0_yzzz_y_yy[i] - rab_z[i] * t_0_yzz_y_yy[i];

        t_z_yzz_y_xz[i] = t_0_yzzz_y_xz[i] - rab_z[i] * t_0_yzz_y_xz[i];

        t_z_yzz_y_xy[i] = t_0_yzzz_y_xy[i] - rab_z[i] * t_0_yzz_y_xy[i];

        t_z_yzz_y_xx[i] = t_0_yzzz_y_xx[i] - rab_z[i] * t_0_yzz_y_xx[i];

        t_z_yzz_x_zz[i] = t_0_yzzz_x_zz[i] - rab_z[i] * t_0_yzz_x_zz[i];

        t_z_yzz_x_yz[i] = t_0_yzzz_x_yz[i] - rab_z[i] * t_0_yzz_x_yz[i];

        t_z_yzz_x_yy[i] = t_0_yzzz_x_yy[i] - rab_z[i] * t_0_yzz_x_yy[i];

        t_z_yzz_x_xz[i] = t_0_yzzz_x_xz[i] - rab_z[i] * t_0_yzz_x_xz[i];

        t_z_yzz_x_xy[i] = t_0_yzzz_x_xy[i] - rab_z[i] * t_0_yzz_x_xy[i];

        t_z_yzz_x_xx[i] = t_0_yzzz_x_xx[i] - rab_z[i] * t_0_yzz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xzz_x_xx, t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy,\
                           t_0_xzz_x_yz, t_0_xzz_x_zz, t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz,\
                           t_0_xzz_y_yy, t_0_xzz_y_yz, t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy,\
                           t_0_xzz_z_xz, t_0_xzz_z_yy, t_0_xzz_z_yz, t_0_xzz_z_zz, t_0_xzzz_x_xx,\
                           t_0_xzzz_x_xy, t_0_xzzz_x_xz, t_0_xzzz_x_yy, t_0_xzzz_x_yz,\
                           t_0_xzzz_x_zz, t_0_xzzz_y_xx, t_0_xzzz_y_xy, t_0_xzzz_y_xz,\
                           t_0_xzzz_y_yy, t_0_xzzz_y_yz, t_0_xzzz_y_zz, t_0_xzzz_z_xx,\
                           t_0_xzzz_z_xy, t_0_xzzz_z_xz, t_0_xzzz_z_yy, t_0_xzzz_z_yz,\
                           t_0_xzzz_z_zz, t_0_yyz_x_xx, t_0_yyz_x_xy, t_0_yyz_x_xz,\
                           t_0_yyz_x_yy, t_0_yyz_x_yz, t_0_yyz_x_zz, t_0_yyz_y_xx, t_0_yyz_y_xy,\
                           t_0_yyz_y_xz, t_0_yyz_y_yy, t_0_yyz_y_yz, t_0_yyz_y_zz, t_0_yyz_z_xx,\
                           t_0_yyz_z_xy, t_0_yyz_z_xz, t_0_yyz_z_yy, t_0_yyz_z_yz, t_0_yyz_z_zz,\
                           t_0_yyzz_x_xx, t_0_yyzz_x_xy, t_0_yyzz_x_xz, t_0_yyzz_x_yy,\
                           t_0_yyzz_x_yz, t_0_yyzz_x_zz, t_0_yyzz_y_xx, t_0_yyzz_y_xy,\
                           t_0_yyzz_y_xz, t_0_yyzz_y_yy, t_0_yyzz_y_yz, t_0_yyzz_y_zz,\
                           t_0_yyzz_z_xx, t_0_yyzz_z_xy, t_0_yyzz_z_xz, t_0_yyzz_z_yy,\
                           t_0_yyzz_z_yz, t_0_yyzz_z_zz, t_z_xzz_x_xx, t_z_xzz_x_xy,\
                           t_z_xzz_x_xz, t_z_xzz_x_yy, t_z_xzz_x_yz, t_z_xzz_x_zz, t_z_xzz_y_xx,\
                           t_z_xzz_y_xy, t_z_xzz_y_xz, t_z_xzz_y_yy, t_z_xzz_y_yz, t_z_xzz_y_zz,\
                           t_z_xzz_z_xx, t_z_xzz_z_xy, t_z_xzz_z_xz, t_z_xzz_z_yy, t_z_xzz_z_yz,\
                           t_z_xzz_z_zz, t_z_yyz_x_xx, t_z_yyz_x_xy, t_z_yyz_x_xz, t_z_yyz_x_yy,\
                           t_z_yyz_x_yz, t_z_yyz_x_zz, t_z_yyz_y_xx, t_z_yyz_y_xy, t_z_yyz_y_xz,\
                           t_z_yyz_y_yy, t_z_yyz_y_yz, t_z_yyz_y_zz, t_z_yyz_z_xx, t_z_yyz_z_xy,\
                           t_z_yyz_z_xz, t_z_yyz_z_yy, t_z_yyz_z_yz, t_z_yyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yyz_z_zz[i] = t_0_yyzz_z_zz[i] - rab_z[i] * t_0_yyz_z_zz[i];

        t_z_yyz_z_yz[i] = t_0_yyzz_z_yz[i] - rab_z[i] * t_0_yyz_z_yz[i];

        t_z_yyz_z_yy[i] = t_0_yyzz_z_yy[i] - rab_z[i] * t_0_yyz_z_yy[i];

        t_z_yyz_z_xz[i] = t_0_yyzz_z_xz[i] - rab_z[i] * t_0_yyz_z_xz[i];

        t_z_yyz_z_xy[i] = t_0_yyzz_z_xy[i] - rab_z[i] * t_0_yyz_z_xy[i];

        t_z_yyz_z_xx[i] = t_0_yyzz_z_xx[i] - rab_z[i] * t_0_yyz_z_xx[i];

        t_z_yyz_y_zz[i] = t_0_yyzz_y_zz[i] - rab_z[i] * t_0_yyz_y_zz[i];

        t_z_yyz_y_yz[i] = t_0_yyzz_y_yz[i] - rab_z[i] * t_0_yyz_y_yz[i];

        t_z_yyz_y_yy[i] = t_0_yyzz_y_yy[i] - rab_z[i] * t_0_yyz_y_yy[i];

        t_z_yyz_y_xz[i] = t_0_yyzz_y_xz[i] - rab_z[i] * t_0_yyz_y_xz[i];

        t_z_yyz_y_xy[i] = t_0_yyzz_y_xy[i] - rab_z[i] * t_0_yyz_y_xy[i];

        t_z_yyz_y_xx[i] = t_0_yyzz_y_xx[i] - rab_z[i] * t_0_yyz_y_xx[i];

        t_z_yyz_x_zz[i] = t_0_yyzz_x_zz[i] - rab_z[i] * t_0_yyz_x_zz[i];

        t_z_yyz_x_yz[i] = t_0_yyzz_x_yz[i] - rab_z[i] * t_0_yyz_x_yz[i];

        t_z_yyz_x_yy[i] = t_0_yyzz_x_yy[i] - rab_z[i] * t_0_yyz_x_yy[i];

        t_z_yyz_x_xz[i] = t_0_yyzz_x_xz[i] - rab_z[i] * t_0_yyz_x_xz[i];

        t_z_yyz_x_xy[i] = t_0_yyzz_x_xy[i] - rab_z[i] * t_0_yyz_x_xy[i];

        t_z_yyz_x_xx[i] = t_0_yyzz_x_xx[i] - rab_z[i] * t_0_yyz_x_xx[i];

        t_z_xzz_z_zz[i] = t_0_xzzz_z_zz[i] - rab_z[i] * t_0_xzz_z_zz[i];

        t_z_xzz_z_yz[i] = t_0_xzzz_z_yz[i] - rab_z[i] * t_0_xzz_z_yz[i];

        t_z_xzz_z_yy[i] = t_0_xzzz_z_yy[i] - rab_z[i] * t_0_xzz_z_yy[i];

        t_z_xzz_z_xz[i] = t_0_xzzz_z_xz[i] - rab_z[i] * t_0_xzz_z_xz[i];

        t_z_xzz_z_xy[i] = t_0_xzzz_z_xy[i] - rab_z[i] * t_0_xzz_z_xy[i];

        t_z_xzz_z_xx[i] = t_0_xzzz_z_xx[i] - rab_z[i] * t_0_xzz_z_xx[i];

        t_z_xzz_y_zz[i] = t_0_xzzz_y_zz[i] - rab_z[i] * t_0_xzz_y_zz[i];

        t_z_xzz_y_yz[i] = t_0_xzzz_y_yz[i] - rab_z[i] * t_0_xzz_y_yz[i];

        t_z_xzz_y_yy[i] = t_0_xzzz_y_yy[i] - rab_z[i] * t_0_xzz_y_yy[i];

        t_z_xzz_y_xz[i] = t_0_xzzz_y_xz[i] - rab_z[i] * t_0_xzz_y_xz[i];

        t_z_xzz_y_xy[i] = t_0_xzzz_y_xy[i] - rab_z[i] * t_0_xzz_y_xy[i];

        t_z_xzz_y_xx[i] = t_0_xzzz_y_xx[i] - rab_z[i] * t_0_xzz_y_xx[i];

        t_z_xzz_x_zz[i] = t_0_xzzz_x_zz[i] - rab_z[i] * t_0_xzz_x_zz[i];

        t_z_xzz_x_yz[i] = t_0_xzzz_x_yz[i] - rab_z[i] * t_0_xzz_x_yz[i];

        t_z_xzz_x_yy[i] = t_0_xzzz_x_yy[i] - rab_z[i] * t_0_xzz_x_yy[i];

        t_z_xzz_x_xz[i] = t_0_xzzz_x_xz[i] - rab_z[i] * t_0_xzz_x_xz[i];

        t_z_xzz_x_xy[i] = t_0_xzzz_x_xy[i] - rab_z[i] * t_0_xzz_x_xy[i];

        t_z_xzz_x_xx[i] = t_0_xzzz_x_xx[i] - rab_z[i] * t_0_xzz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xxz_x_xx, t_0_xxz_x_xy, t_0_xxz_x_xz, t_0_xxz_x_yy,\
                           t_0_xxz_x_yz, t_0_xxz_x_zz, t_0_xxz_y_xx, t_0_xxz_y_xy, t_0_xxz_y_xz,\
                           t_0_xxz_y_yy, t_0_xxz_y_yz, t_0_xxz_y_zz, t_0_xxz_z_xx, t_0_xxz_z_xy,\
                           t_0_xxz_z_xz, t_0_xxz_z_yy, t_0_xxz_z_yz, t_0_xxz_z_zz, t_0_xxzz_x_xx,\
                           t_0_xxzz_x_xy, t_0_xxzz_x_xz, t_0_xxzz_x_yy, t_0_xxzz_x_yz,\
                           t_0_xxzz_x_zz, t_0_xxzz_y_xx, t_0_xxzz_y_xy, t_0_xxzz_y_xz,\
                           t_0_xxzz_y_yy, t_0_xxzz_y_yz, t_0_xxzz_y_zz, t_0_xxzz_z_xx,\
                           t_0_xxzz_z_xy, t_0_xxzz_z_xz, t_0_xxzz_z_yy, t_0_xxzz_z_yz,\
                           t_0_xxzz_z_zz, t_0_xyz_x_xx, t_0_xyz_x_xy, t_0_xyz_x_xz,\
                           t_0_xyz_x_yy, t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx, t_0_xyz_y_xy,\
                           t_0_xyz_y_xz, t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz, t_0_xyz_z_xx,\
                           t_0_xyz_z_xy, t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz, t_0_xyz_z_zz,\
                           t_0_xyzz_x_xx, t_0_xyzz_x_xy, t_0_xyzz_x_xz, t_0_xyzz_x_yy,\
                           t_0_xyzz_x_yz, t_0_xyzz_x_zz, t_0_xyzz_y_xx, t_0_xyzz_y_xy,\
                           t_0_xyzz_y_xz, t_0_xyzz_y_yy, t_0_xyzz_y_yz, t_0_xyzz_y_zz,\
                           t_0_xyzz_z_xx, t_0_xyzz_z_xy, t_0_xyzz_z_xz, t_0_xyzz_z_yy,\
                           t_0_xyzz_z_yz, t_0_xyzz_z_zz, t_z_xxz_x_xx, t_z_xxz_x_xy,\
                           t_z_xxz_x_xz, t_z_xxz_x_yy, t_z_xxz_x_yz, t_z_xxz_x_zz, t_z_xxz_y_xx,\
                           t_z_xxz_y_xy, t_z_xxz_y_xz, t_z_xxz_y_yy, t_z_xxz_y_yz, t_z_xxz_y_zz,\
                           t_z_xxz_z_xx, t_z_xxz_z_xy, t_z_xxz_z_xz, t_z_xxz_z_yy, t_z_xxz_z_yz,\
                           t_z_xxz_z_zz, t_z_xyz_x_xx, t_z_xyz_x_xy, t_z_xyz_x_xz, t_z_xyz_x_yy,\
                           t_z_xyz_x_yz, t_z_xyz_x_zz, t_z_xyz_y_xx, t_z_xyz_y_xy, t_z_xyz_y_xz,\
                           t_z_xyz_y_yy, t_z_xyz_y_yz, t_z_xyz_y_zz, t_z_xyz_z_xx, t_z_xyz_z_xy,\
                           t_z_xyz_z_xz, t_z_xyz_z_yy, t_z_xyz_z_yz, t_z_xyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xyz_z_zz[i] = t_0_xyzz_z_zz[i] - rab_z[i] * t_0_xyz_z_zz[i];

        t_z_xyz_z_yz[i] = t_0_xyzz_z_yz[i] - rab_z[i] * t_0_xyz_z_yz[i];

        t_z_xyz_z_yy[i] = t_0_xyzz_z_yy[i] - rab_z[i] * t_0_xyz_z_yy[i];

        t_z_xyz_z_xz[i] = t_0_xyzz_z_xz[i] - rab_z[i] * t_0_xyz_z_xz[i];

        t_z_xyz_z_xy[i] = t_0_xyzz_z_xy[i] - rab_z[i] * t_0_xyz_z_xy[i];

        t_z_xyz_z_xx[i] = t_0_xyzz_z_xx[i] - rab_z[i] * t_0_xyz_z_xx[i];

        t_z_xyz_y_zz[i] = t_0_xyzz_y_zz[i] - rab_z[i] * t_0_xyz_y_zz[i];

        t_z_xyz_y_yz[i] = t_0_xyzz_y_yz[i] - rab_z[i] * t_0_xyz_y_yz[i];

        t_z_xyz_y_yy[i] = t_0_xyzz_y_yy[i] - rab_z[i] * t_0_xyz_y_yy[i];

        t_z_xyz_y_xz[i] = t_0_xyzz_y_xz[i] - rab_z[i] * t_0_xyz_y_xz[i];

        t_z_xyz_y_xy[i] = t_0_xyzz_y_xy[i] - rab_z[i] * t_0_xyz_y_xy[i];

        t_z_xyz_y_xx[i] = t_0_xyzz_y_xx[i] - rab_z[i] * t_0_xyz_y_xx[i];

        t_z_xyz_x_zz[i] = t_0_xyzz_x_zz[i] - rab_z[i] * t_0_xyz_x_zz[i];

        t_z_xyz_x_yz[i] = t_0_xyzz_x_yz[i] - rab_z[i] * t_0_xyz_x_yz[i];

        t_z_xyz_x_yy[i] = t_0_xyzz_x_yy[i] - rab_z[i] * t_0_xyz_x_yy[i];

        t_z_xyz_x_xz[i] = t_0_xyzz_x_xz[i] - rab_z[i] * t_0_xyz_x_xz[i];

        t_z_xyz_x_xy[i] = t_0_xyzz_x_xy[i] - rab_z[i] * t_0_xyz_x_xy[i];

        t_z_xyz_x_xx[i] = t_0_xyzz_x_xx[i] - rab_z[i] * t_0_xyz_x_xx[i];

        t_z_xxz_z_zz[i] = t_0_xxzz_z_zz[i] - rab_z[i] * t_0_xxz_z_zz[i];

        t_z_xxz_z_yz[i] = t_0_xxzz_z_yz[i] - rab_z[i] * t_0_xxz_z_yz[i];

        t_z_xxz_z_yy[i] = t_0_xxzz_z_yy[i] - rab_z[i] * t_0_xxz_z_yy[i];

        t_z_xxz_z_xz[i] = t_0_xxzz_z_xz[i] - rab_z[i] * t_0_xxz_z_xz[i];

        t_z_xxz_z_xy[i] = t_0_xxzz_z_xy[i] - rab_z[i] * t_0_xxz_z_xy[i];

        t_z_xxz_z_xx[i] = t_0_xxzz_z_xx[i] - rab_z[i] * t_0_xxz_z_xx[i];

        t_z_xxz_y_zz[i] = t_0_xxzz_y_zz[i] - rab_z[i] * t_0_xxz_y_zz[i];

        t_z_xxz_y_yz[i] = t_0_xxzz_y_yz[i] - rab_z[i] * t_0_xxz_y_yz[i];

        t_z_xxz_y_yy[i] = t_0_xxzz_y_yy[i] - rab_z[i] * t_0_xxz_y_yy[i];

        t_z_xxz_y_xz[i] = t_0_xxzz_y_xz[i] - rab_z[i] * t_0_xxz_y_xz[i];

        t_z_xxz_y_xy[i] = t_0_xxzz_y_xy[i] - rab_z[i] * t_0_xxz_y_xy[i];

        t_z_xxz_y_xx[i] = t_0_xxzz_y_xx[i] - rab_z[i] * t_0_xxz_y_xx[i];

        t_z_xxz_x_zz[i] = t_0_xxzz_x_zz[i] - rab_z[i] * t_0_xxz_x_zz[i];

        t_z_xxz_x_yz[i] = t_0_xxzz_x_yz[i] - rab_z[i] * t_0_xxz_x_yz[i];

        t_z_xxz_x_yy[i] = t_0_xxzz_x_yy[i] - rab_z[i] * t_0_xxz_x_yy[i];

        t_z_xxz_x_xz[i] = t_0_xxzz_x_xz[i] - rab_z[i] * t_0_xxz_x_xz[i];

        t_z_xxz_x_xy[i] = t_0_xxzz_x_xy[i] - rab_z[i] * t_0_xxz_x_xy[i];

        t_z_xxz_x_xx[i] = t_0_xxzz_x_xx[i] - rab_z[i] * t_0_xxz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyzz_x_xx, t_0_yyzz_x_xy, t_0_yyzz_x_xz, t_0_yyzz_x_yy,\
                           t_0_yyzz_x_yz, t_0_yyzz_x_zz, t_0_yyzz_y_xx, t_0_yyzz_y_xy,\
                           t_0_yyzz_y_xz, t_0_yyzz_y_yy, t_0_yyzz_y_yz, t_0_yyzz_y_zz,\
                           t_0_yyzz_z_xx, t_0_yyzz_z_xy, t_0_yyzz_z_xz, t_0_yyzz_z_yy,\
                           t_0_yyzz_z_yz, t_0_yyzz_z_zz, t_0_yzz_x_xx, t_0_yzz_x_xy,\
                           t_0_yzz_x_xz, t_0_yzz_x_yy, t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_xx,\
                           t_0_yzz_y_xy, t_0_yzz_y_xz, t_0_yzz_y_yy, t_0_yzz_y_yz, t_0_yzz_y_zz,\
                           t_0_yzz_z_xx, t_0_yzz_z_xy, t_0_yzz_z_xz, t_0_yzz_z_yy, t_0_yzz_z_yz,\
                           t_0_yzz_z_zz, t_0_yzzz_x_xx, t_0_yzzz_x_xy, t_0_yzzz_x_xz,\
                           t_0_yzzz_x_yy, t_0_yzzz_x_yz, t_0_yzzz_x_zz, t_0_yzzz_y_xx,\
                           t_0_yzzz_y_xy, t_0_yzzz_y_xz, t_0_yzzz_y_yy, t_0_yzzz_y_yz,\
                           t_0_yzzz_y_zz, t_0_yzzz_z_xx, t_0_yzzz_z_xy, t_0_yzzz_z_xz,\
                           t_0_yzzz_z_yy, t_0_yzzz_z_yz, t_0_yzzz_z_zz, t_0_zzz_x_xx,\
                           t_0_zzz_x_xy, t_0_zzz_x_xz, t_0_zzz_x_yy, t_0_zzz_x_yz, t_0_zzz_x_zz,\
                           t_0_zzz_y_xx, t_0_zzz_y_xy, t_0_zzz_y_xz, t_0_zzz_y_yy, t_0_zzz_y_yz,\
                           t_0_zzz_y_zz, t_0_zzz_z_xx, t_0_zzz_z_xy, t_0_zzz_z_xz, t_0_zzz_z_yy,\
                           t_0_zzz_z_yz, t_0_zzz_z_zz, t_y_yzz_x_xx, t_y_yzz_x_xy, t_y_yzz_x_xz,\
                           t_y_yzz_x_yy, t_y_yzz_x_yz, t_y_yzz_x_zz, t_y_yzz_y_xx, t_y_yzz_y_xy,\
                           t_y_yzz_y_xz, t_y_yzz_y_yy, t_y_yzz_y_yz, t_y_yzz_y_zz, t_y_yzz_z_xx,\
                           t_y_yzz_z_xy, t_y_yzz_z_xz, t_y_yzz_z_yy, t_y_yzz_z_yz, t_y_yzz_z_zz,\
                           t_y_zzz_x_xx, t_y_zzz_x_xy, t_y_zzz_x_xz, t_y_zzz_x_yy, t_y_zzz_x_yz,\
                           t_y_zzz_x_zz, t_y_zzz_y_xx, t_y_zzz_y_xy, t_y_zzz_y_xz, t_y_zzz_y_yy,\
                           t_y_zzz_y_yz, t_y_zzz_y_zz, t_y_zzz_z_xx, t_y_zzz_z_xy, t_y_zzz_z_xz,\
                           t_y_zzz_z_yy, t_y_zzz_z_yz, t_y_zzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zzz_z_zz[i] = t_0_yzzz_z_zz[i] - rab_y[i] * t_0_zzz_z_zz[i];

        t_y_zzz_z_yz[i] = t_0_yzzz_z_yz[i] - rab_y[i] * t_0_zzz_z_yz[i];

        t_y_zzz_z_yy[i] = t_0_yzzz_z_yy[i] - rab_y[i] * t_0_zzz_z_yy[i];

        t_y_zzz_z_xz[i] = t_0_yzzz_z_xz[i] - rab_y[i] * t_0_zzz_z_xz[i];

        t_y_zzz_z_xy[i] = t_0_yzzz_z_xy[i] - rab_y[i] * t_0_zzz_z_xy[i];

        t_y_zzz_z_xx[i] = t_0_yzzz_z_xx[i] - rab_y[i] * t_0_zzz_z_xx[i];

        t_y_zzz_y_zz[i] = t_0_yzzz_y_zz[i] - rab_y[i] * t_0_zzz_y_zz[i];

        t_y_zzz_y_yz[i] = t_0_yzzz_y_yz[i] - rab_y[i] * t_0_zzz_y_yz[i];

        t_y_zzz_y_yy[i] = t_0_yzzz_y_yy[i] - rab_y[i] * t_0_zzz_y_yy[i];

        t_y_zzz_y_xz[i] = t_0_yzzz_y_xz[i] - rab_y[i] * t_0_zzz_y_xz[i];

        t_y_zzz_y_xy[i] = t_0_yzzz_y_xy[i] - rab_y[i] * t_0_zzz_y_xy[i];

        t_y_zzz_y_xx[i] = t_0_yzzz_y_xx[i] - rab_y[i] * t_0_zzz_y_xx[i];

        t_y_zzz_x_zz[i] = t_0_yzzz_x_zz[i] - rab_y[i] * t_0_zzz_x_zz[i];

        t_y_zzz_x_yz[i] = t_0_yzzz_x_yz[i] - rab_y[i] * t_0_zzz_x_yz[i];

        t_y_zzz_x_yy[i] = t_0_yzzz_x_yy[i] - rab_y[i] * t_0_zzz_x_yy[i];

        t_y_zzz_x_xz[i] = t_0_yzzz_x_xz[i] - rab_y[i] * t_0_zzz_x_xz[i];

        t_y_zzz_x_xy[i] = t_0_yzzz_x_xy[i] - rab_y[i] * t_0_zzz_x_xy[i];

        t_y_zzz_x_xx[i] = t_0_yzzz_x_xx[i] - rab_y[i] * t_0_zzz_x_xx[i];

        t_y_yzz_z_zz[i] = t_0_yyzz_z_zz[i] - rab_y[i] * t_0_yzz_z_zz[i];

        t_y_yzz_z_yz[i] = t_0_yyzz_z_yz[i] - rab_y[i] * t_0_yzz_z_yz[i];

        t_y_yzz_z_yy[i] = t_0_yyzz_z_yy[i] - rab_y[i] * t_0_yzz_z_yy[i];

        t_y_yzz_z_xz[i] = t_0_yyzz_z_xz[i] - rab_y[i] * t_0_yzz_z_xz[i];

        t_y_yzz_z_xy[i] = t_0_yyzz_z_xy[i] - rab_y[i] * t_0_yzz_z_xy[i];

        t_y_yzz_z_xx[i] = t_0_yyzz_z_xx[i] - rab_y[i] * t_0_yzz_z_xx[i];

        t_y_yzz_y_zz[i] = t_0_yyzz_y_zz[i] - rab_y[i] * t_0_yzz_y_zz[i];

        t_y_yzz_y_yz[i] = t_0_yyzz_y_yz[i] - rab_y[i] * t_0_yzz_y_yz[i];

        t_y_yzz_y_yy[i] = t_0_yyzz_y_yy[i] - rab_y[i] * t_0_yzz_y_yy[i];

        t_y_yzz_y_xz[i] = t_0_yyzz_y_xz[i] - rab_y[i] * t_0_yzz_y_xz[i];

        t_y_yzz_y_xy[i] = t_0_yyzz_y_xy[i] - rab_y[i] * t_0_yzz_y_xy[i];

        t_y_yzz_y_xx[i] = t_0_yyzz_y_xx[i] - rab_y[i] * t_0_yzz_y_xx[i];

        t_y_yzz_x_zz[i] = t_0_yyzz_x_zz[i] - rab_y[i] * t_0_yzz_x_zz[i];

        t_y_yzz_x_yz[i] = t_0_yyzz_x_yz[i] - rab_y[i] * t_0_yzz_x_yz[i];

        t_y_yzz_x_yy[i] = t_0_yyzz_x_yy[i] - rab_y[i] * t_0_yzz_x_yy[i];

        t_y_yzz_x_xz[i] = t_0_yyzz_x_xz[i] - rab_y[i] * t_0_yzz_x_xz[i];

        t_y_yzz_x_xy[i] = t_0_yyzz_x_xy[i] - rab_y[i] * t_0_yzz_x_xy[i];

        t_y_yzz_x_xx[i] = t_0_yyzz_x_xx[i] - rab_y[i] * t_0_yzz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyy_x_xx, t_0_yyy_x_xy, t_0_yyy_x_xz, t_0_yyy_x_yy,\
                           t_0_yyy_x_yz, t_0_yyy_x_zz, t_0_yyy_y_xx, t_0_yyy_y_xy, t_0_yyy_y_xz,\
                           t_0_yyy_y_yy, t_0_yyy_y_yz, t_0_yyy_y_zz, t_0_yyy_z_xx, t_0_yyy_z_xy,\
                           t_0_yyy_z_xz, t_0_yyy_z_yy, t_0_yyy_z_yz, t_0_yyy_z_zz, t_0_yyyy_x_xx,\
                           t_0_yyyy_x_xy, t_0_yyyy_x_xz, t_0_yyyy_x_yy, t_0_yyyy_x_yz,\
                           t_0_yyyy_x_zz, t_0_yyyy_y_xx, t_0_yyyy_y_xy, t_0_yyyy_y_xz,\
                           t_0_yyyy_y_yy, t_0_yyyy_y_yz, t_0_yyyy_y_zz, t_0_yyyy_z_xx,\
                           t_0_yyyy_z_xy, t_0_yyyy_z_xz, t_0_yyyy_z_yy, t_0_yyyy_z_yz,\
                           t_0_yyyy_z_zz, t_0_yyyz_x_xx, t_0_yyyz_x_xy, t_0_yyyz_x_xz,\
                           t_0_yyyz_x_yy, t_0_yyyz_x_yz, t_0_yyyz_x_zz, t_0_yyyz_y_xx,\
                           t_0_yyyz_y_xy, t_0_yyyz_y_xz, t_0_yyyz_y_yy, t_0_yyyz_y_yz,\
                           t_0_yyyz_y_zz, t_0_yyyz_z_xx, t_0_yyyz_z_xy, t_0_yyyz_z_xz,\
                           t_0_yyyz_z_yy, t_0_yyyz_z_yz, t_0_yyyz_z_zz, t_0_yyz_x_xx,\
                           t_0_yyz_x_xy, t_0_yyz_x_xz, t_0_yyz_x_yy, t_0_yyz_x_yz, t_0_yyz_x_zz,\
                           t_0_yyz_y_xx, t_0_yyz_y_xy, t_0_yyz_y_xz, t_0_yyz_y_yy, t_0_yyz_y_yz,\
                           t_0_yyz_y_zz, t_0_yyz_z_xx, t_0_yyz_z_xy, t_0_yyz_z_xz, t_0_yyz_z_yy,\
                           t_0_yyz_z_yz, t_0_yyz_z_zz, t_y_yyy_x_xx, t_y_yyy_x_xy, t_y_yyy_x_xz,\
                           t_y_yyy_x_yy, t_y_yyy_x_yz, t_y_yyy_x_zz, t_y_yyy_y_xx, t_y_yyy_y_xy,\
                           t_y_yyy_y_xz, t_y_yyy_y_yy, t_y_yyy_y_yz, t_y_yyy_y_zz, t_y_yyy_z_xx,\
                           t_y_yyy_z_xy, t_y_yyy_z_xz, t_y_yyy_z_yy, t_y_yyy_z_yz, t_y_yyy_z_zz,\
                           t_y_yyz_x_xx, t_y_yyz_x_xy, t_y_yyz_x_xz, t_y_yyz_x_yy, t_y_yyz_x_yz,\
                           t_y_yyz_x_zz, t_y_yyz_y_xx, t_y_yyz_y_xy, t_y_yyz_y_xz, t_y_yyz_y_yy,\
                           t_y_yyz_y_yz, t_y_yyz_y_zz, t_y_yyz_z_xx, t_y_yyz_z_xy, t_y_yyz_z_xz,\
                           t_y_yyz_z_yy, t_y_yyz_z_yz, t_y_yyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yyz_z_zz[i] = t_0_yyyz_z_zz[i] - rab_y[i] * t_0_yyz_z_zz[i];

        t_y_yyz_z_yz[i] = t_0_yyyz_z_yz[i] - rab_y[i] * t_0_yyz_z_yz[i];

        t_y_yyz_z_yy[i] = t_0_yyyz_z_yy[i] - rab_y[i] * t_0_yyz_z_yy[i];

        t_y_yyz_z_xz[i] = t_0_yyyz_z_xz[i] - rab_y[i] * t_0_yyz_z_xz[i];

        t_y_yyz_z_xy[i] = t_0_yyyz_z_xy[i] - rab_y[i] * t_0_yyz_z_xy[i];

        t_y_yyz_z_xx[i] = t_0_yyyz_z_xx[i] - rab_y[i] * t_0_yyz_z_xx[i];

        t_y_yyz_y_zz[i] = t_0_yyyz_y_zz[i] - rab_y[i] * t_0_yyz_y_zz[i];

        t_y_yyz_y_yz[i] = t_0_yyyz_y_yz[i] - rab_y[i] * t_0_yyz_y_yz[i];

        t_y_yyz_y_yy[i] = t_0_yyyz_y_yy[i] - rab_y[i] * t_0_yyz_y_yy[i];

        t_y_yyz_y_xz[i] = t_0_yyyz_y_xz[i] - rab_y[i] * t_0_yyz_y_xz[i];

        t_y_yyz_y_xy[i] = t_0_yyyz_y_xy[i] - rab_y[i] * t_0_yyz_y_xy[i];

        t_y_yyz_y_xx[i] = t_0_yyyz_y_xx[i] - rab_y[i] * t_0_yyz_y_xx[i];

        t_y_yyz_x_zz[i] = t_0_yyyz_x_zz[i] - rab_y[i] * t_0_yyz_x_zz[i];

        t_y_yyz_x_yz[i] = t_0_yyyz_x_yz[i] - rab_y[i] * t_0_yyz_x_yz[i];

        t_y_yyz_x_yy[i] = t_0_yyyz_x_yy[i] - rab_y[i] * t_0_yyz_x_yy[i];

        t_y_yyz_x_xz[i] = t_0_yyyz_x_xz[i] - rab_y[i] * t_0_yyz_x_xz[i];

        t_y_yyz_x_xy[i] = t_0_yyyz_x_xy[i] - rab_y[i] * t_0_yyz_x_xy[i];

        t_y_yyz_x_xx[i] = t_0_yyyz_x_xx[i] - rab_y[i] * t_0_yyz_x_xx[i];

        t_y_yyy_z_zz[i] = t_0_yyyy_z_zz[i] - rab_y[i] * t_0_yyy_z_zz[i];

        t_y_yyy_z_yz[i] = t_0_yyyy_z_yz[i] - rab_y[i] * t_0_yyy_z_yz[i];

        t_y_yyy_z_yy[i] = t_0_yyyy_z_yy[i] - rab_y[i] * t_0_yyy_z_yy[i];

        t_y_yyy_z_xz[i] = t_0_yyyy_z_xz[i] - rab_y[i] * t_0_yyy_z_xz[i];

        t_y_yyy_z_xy[i] = t_0_yyyy_z_xy[i] - rab_y[i] * t_0_yyy_z_xy[i];

        t_y_yyy_z_xx[i] = t_0_yyyy_z_xx[i] - rab_y[i] * t_0_yyy_z_xx[i];

        t_y_yyy_y_zz[i] = t_0_yyyy_y_zz[i] - rab_y[i] * t_0_yyy_y_zz[i];

        t_y_yyy_y_yz[i] = t_0_yyyy_y_yz[i] - rab_y[i] * t_0_yyy_y_yz[i];

        t_y_yyy_y_yy[i] = t_0_yyyy_y_yy[i] - rab_y[i] * t_0_yyy_y_yy[i];

        t_y_yyy_y_xz[i] = t_0_yyyy_y_xz[i] - rab_y[i] * t_0_yyy_y_xz[i];

        t_y_yyy_y_xy[i] = t_0_yyyy_y_xy[i] - rab_y[i] * t_0_yyy_y_xy[i];

        t_y_yyy_y_xx[i] = t_0_yyyy_y_xx[i] - rab_y[i] * t_0_yyy_y_xx[i];

        t_y_yyy_x_zz[i] = t_0_yyyy_x_zz[i] - rab_y[i] * t_0_yyy_x_zz[i];

        t_y_yyy_x_yz[i] = t_0_yyyy_x_yz[i] - rab_y[i] * t_0_yyy_x_yz[i];

        t_y_yyy_x_yy[i] = t_0_yyyy_x_yy[i] - rab_y[i] * t_0_yyy_x_yy[i];

        t_y_yyy_x_xz[i] = t_0_yyyy_x_xz[i] - rab_y[i] * t_0_yyy_x_xz[i];

        t_y_yyy_x_xy[i] = t_0_yyyy_x_xy[i] - rab_y[i] * t_0_yyy_x_xy[i];

        t_y_yyy_x_xx[i] = t_0_yyyy_x_xx[i] - rab_y[i] * t_0_yyy_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyyz_x_xx, t_0_xyyz_x_xy, t_0_xyyz_x_xz, t_0_xyyz_x_yy,\
                           t_0_xyyz_x_yz, t_0_xyyz_x_zz, t_0_xyyz_y_xx, t_0_xyyz_y_xy,\
                           t_0_xyyz_y_xz, t_0_xyyz_y_yy, t_0_xyyz_y_yz, t_0_xyyz_y_zz,\
                           t_0_xyyz_z_xx, t_0_xyyz_z_xy, t_0_xyyz_z_xz, t_0_xyyz_z_yy,\
                           t_0_xyyz_z_yz, t_0_xyyz_z_zz, t_0_xyz_x_xx, t_0_xyz_x_xy,\
                           t_0_xyz_x_xz, t_0_xyz_x_yy, t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx,\
                           t_0_xyz_y_xy, t_0_xyz_y_xz, t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz,\
                           t_0_xyz_z_xx, t_0_xyz_z_xy, t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz,\
                           t_0_xyz_z_zz, t_0_xyzz_x_xx, t_0_xyzz_x_xy, t_0_xyzz_x_xz,\
                           t_0_xyzz_x_yy, t_0_xyzz_x_yz, t_0_xyzz_x_zz, t_0_xyzz_y_xx,\
                           t_0_xyzz_y_xy, t_0_xyzz_y_xz, t_0_xyzz_y_yy, t_0_xyzz_y_yz,\
                           t_0_xyzz_y_zz, t_0_xyzz_z_xx, t_0_xyzz_z_xy, t_0_xyzz_z_xz,\
                           t_0_xyzz_z_yy, t_0_xyzz_z_yz, t_0_xyzz_z_zz, t_0_xzz_x_xx,\
                           t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy, t_0_xzz_x_yz, t_0_xzz_x_zz,\
                           t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz, t_0_xzz_y_yy, t_0_xzz_y_yz,\
                           t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy, t_0_xzz_z_xz, t_0_xzz_z_yy,\
                           t_0_xzz_z_yz, t_0_xzz_z_zz, t_y_xyz_x_xx, t_y_xyz_x_xy, t_y_xyz_x_xz,\
                           t_y_xyz_x_yy, t_y_xyz_x_yz, t_y_xyz_x_zz, t_y_xyz_y_xx, t_y_xyz_y_xy,\
                           t_y_xyz_y_xz, t_y_xyz_y_yy, t_y_xyz_y_yz, t_y_xyz_y_zz, t_y_xyz_z_xx,\
                           t_y_xyz_z_xy, t_y_xyz_z_xz, t_y_xyz_z_yy, t_y_xyz_z_yz, t_y_xyz_z_zz,\
                           t_y_xzz_x_xx, t_y_xzz_x_xy, t_y_xzz_x_xz, t_y_xzz_x_yy, t_y_xzz_x_yz,\
                           t_y_xzz_x_zz, t_y_xzz_y_xx, t_y_xzz_y_xy, t_y_xzz_y_xz, t_y_xzz_y_yy,\
                           t_y_xzz_y_yz, t_y_xzz_y_zz, t_y_xzz_z_xx, t_y_xzz_z_xy, t_y_xzz_z_xz,\
                           t_y_xzz_z_yy, t_y_xzz_z_yz, t_y_xzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xzz_z_zz[i] = t_0_xyzz_z_zz[i] - rab_y[i] * t_0_xzz_z_zz[i];

        t_y_xzz_z_yz[i] = t_0_xyzz_z_yz[i] - rab_y[i] * t_0_xzz_z_yz[i];

        t_y_xzz_z_yy[i] = t_0_xyzz_z_yy[i] - rab_y[i] * t_0_xzz_z_yy[i];

        t_y_xzz_z_xz[i] = t_0_xyzz_z_xz[i] - rab_y[i] * t_0_xzz_z_xz[i];

        t_y_xzz_z_xy[i] = t_0_xyzz_z_xy[i] - rab_y[i] * t_0_xzz_z_xy[i];

        t_y_xzz_z_xx[i] = t_0_xyzz_z_xx[i] - rab_y[i] * t_0_xzz_z_xx[i];

        t_y_xzz_y_zz[i] = t_0_xyzz_y_zz[i] - rab_y[i] * t_0_xzz_y_zz[i];

        t_y_xzz_y_yz[i] = t_0_xyzz_y_yz[i] - rab_y[i] * t_0_xzz_y_yz[i];

        t_y_xzz_y_yy[i] = t_0_xyzz_y_yy[i] - rab_y[i] * t_0_xzz_y_yy[i];

        t_y_xzz_y_xz[i] = t_0_xyzz_y_xz[i] - rab_y[i] * t_0_xzz_y_xz[i];

        t_y_xzz_y_xy[i] = t_0_xyzz_y_xy[i] - rab_y[i] * t_0_xzz_y_xy[i];

        t_y_xzz_y_xx[i] = t_0_xyzz_y_xx[i] - rab_y[i] * t_0_xzz_y_xx[i];

        t_y_xzz_x_zz[i] = t_0_xyzz_x_zz[i] - rab_y[i] * t_0_xzz_x_zz[i];

        t_y_xzz_x_yz[i] = t_0_xyzz_x_yz[i] - rab_y[i] * t_0_xzz_x_yz[i];

        t_y_xzz_x_yy[i] = t_0_xyzz_x_yy[i] - rab_y[i] * t_0_xzz_x_yy[i];

        t_y_xzz_x_xz[i] = t_0_xyzz_x_xz[i] - rab_y[i] * t_0_xzz_x_xz[i];

        t_y_xzz_x_xy[i] = t_0_xyzz_x_xy[i] - rab_y[i] * t_0_xzz_x_xy[i];

        t_y_xzz_x_xx[i] = t_0_xyzz_x_xx[i] - rab_y[i] * t_0_xzz_x_xx[i];

        t_y_xyz_z_zz[i] = t_0_xyyz_z_zz[i] - rab_y[i] * t_0_xyz_z_zz[i];

        t_y_xyz_z_yz[i] = t_0_xyyz_z_yz[i] - rab_y[i] * t_0_xyz_z_yz[i];

        t_y_xyz_z_yy[i] = t_0_xyyz_z_yy[i] - rab_y[i] * t_0_xyz_z_yy[i];

        t_y_xyz_z_xz[i] = t_0_xyyz_z_xz[i] - rab_y[i] * t_0_xyz_z_xz[i];

        t_y_xyz_z_xy[i] = t_0_xyyz_z_xy[i] - rab_y[i] * t_0_xyz_z_xy[i];

        t_y_xyz_z_xx[i] = t_0_xyyz_z_xx[i] - rab_y[i] * t_0_xyz_z_xx[i];

        t_y_xyz_y_zz[i] = t_0_xyyz_y_zz[i] - rab_y[i] * t_0_xyz_y_zz[i];

        t_y_xyz_y_yz[i] = t_0_xyyz_y_yz[i] - rab_y[i] * t_0_xyz_y_yz[i];

        t_y_xyz_y_yy[i] = t_0_xyyz_y_yy[i] - rab_y[i] * t_0_xyz_y_yy[i];

        t_y_xyz_y_xz[i] = t_0_xyyz_y_xz[i] - rab_y[i] * t_0_xyz_y_xz[i];

        t_y_xyz_y_xy[i] = t_0_xyyz_y_xy[i] - rab_y[i] * t_0_xyz_y_xy[i];

        t_y_xyz_y_xx[i] = t_0_xyyz_y_xx[i] - rab_y[i] * t_0_xyz_y_xx[i];

        t_y_xyz_x_zz[i] = t_0_xyyz_x_zz[i] - rab_y[i] * t_0_xyz_x_zz[i];

        t_y_xyz_x_yz[i] = t_0_xyyz_x_yz[i] - rab_y[i] * t_0_xyz_x_yz[i];

        t_y_xyz_x_yy[i] = t_0_xyyz_x_yy[i] - rab_y[i] * t_0_xyz_x_yy[i];

        t_y_xyz_x_xz[i] = t_0_xyyz_x_xz[i] - rab_y[i] * t_0_xyz_x_xz[i];

        t_y_xyz_x_xy[i] = t_0_xyyz_x_xy[i] - rab_y[i] * t_0_xyz_x_xy[i];

        t_y_xyz_x_xx[i] = t_0_xyyz_x_xx[i] - rab_y[i] * t_0_xyz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xxyz_x_xx, t_0_xxyz_x_xy, t_0_xxyz_x_xz, t_0_xxyz_x_yy,\
                           t_0_xxyz_x_yz, t_0_xxyz_x_zz, t_0_xxyz_y_xx, t_0_xxyz_y_xy,\
                           t_0_xxyz_y_xz, t_0_xxyz_y_yy, t_0_xxyz_y_yz, t_0_xxyz_y_zz,\
                           t_0_xxyz_z_xx, t_0_xxyz_z_xy, t_0_xxyz_z_xz, t_0_xxyz_z_yy,\
                           t_0_xxyz_z_yz, t_0_xxyz_z_zz, t_0_xxz_x_xx, t_0_xxz_x_xy,\
                           t_0_xxz_x_xz, t_0_xxz_x_yy, t_0_xxz_x_yz, t_0_xxz_x_zz, t_0_xxz_y_xx,\
                           t_0_xxz_y_xy, t_0_xxz_y_xz, t_0_xxz_y_yy, t_0_xxz_y_yz, t_0_xxz_y_zz,\
                           t_0_xxz_z_xx, t_0_xxz_z_xy, t_0_xxz_z_xz, t_0_xxz_z_yy, t_0_xxz_z_yz,\
                           t_0_xxz_z_zz, t_0_xyy_x_xx, t_0_xyy_x_xy, t_0_xyy_x_xz, t_0_xyy_x_yy,\
                           t_0_xyy_x_yz, t_0_xyy_x_zz, t_0_xyy_y_xx, t_0_xyy_y_xy, t_0_xyy_y_xz,\
                           t_0_xyy_y_yy, t_0_xyy_y_yz, t_0_xyy_y_zz, t_0_xyy_z_xx, t_0_xyy_z_xy,\
                           t_0_xyy_z_xz, t_0_xyy_z_yy, t_0_xyy_z_yz, t_0_xyy_z_zz, t_0_xyyy_x_xx,\
                           t_0_xyyy_x_xy, t_0_xyyy_x_xz, t_0_xyyy_x_yy, t_0_xyyy_x_yz,\
                           t_0_xyyy_x_zz, t_0_xyyy_y_xx, t_0_xyyy_y_xy, t_0_xyyy_y_xz,\
                           t_0_xyyy_y_yy, t_0_xyyy_y_yz, t_0_xyyy_y_zz, t_0_xyyy_z_xx,\
                           t_0_xyyy_z_xy, t_0_xyyy_z_xz, t_0_xyyy_z_yy, t_0_xyyy_z_yz,\
                           t_0_xyyy_z_zz, t_y_xxz_x_xx, t_y_xxz_x_xy, t_y_xxz_x_xz,\
                           t_y_xxz_x_yy, t_y_xxz_x_yz, t_y_xxz_x_zz, t_y_xxz_y_xx, t_y_xxz_y_xy,\
                           t_y_xxz_y_xz, t_y_xxz_y_yy, t_y_xxz_y_yz, t_y_xxz_y_zz, t_y_xxz_z_xx,\
                           t_y_xxz_z_xy, t_y_xxz_z_xz, t_y_xxz_z_yy, t_y_xxz_z_yz, t_y_xxz_z_zz,\
                           t_y_xyy_x_xx, t_y_xyy_x_xy, t_y_xyy_x_xz, t_y_xyy_x_yy, t_y_xyy_x_yz,\
                           t_y_xyy_x_zz, t_y_xyy_y_xx, t_y_xyy_y_xy, t_y_xyy_y_xz, t_y_xyy_y_yy,\
                           t_y_xyy_y_yz, t_y_xyy_y_zz, t_y_xyy_z_xx, t_y_xyy_z_xy, t_y_xyy_z_xz,\
                           t_y_xyy_z_yy, t_y_xyy_z_yz, t_y_xyy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyy_z_zz[i] = t_0_xyyy_z_zz[i] - rab_y[i] * t_0_xyy_z_zz[i];

        t_y_xyy_z_yz[i] = t_0_xyyy_z_yz[i] - rab_y[i] * t_0_xyy_z_yz[i];

        t_y_xyy_z_yy[i] = t_0_xyyy_z_yy[i] - rab_y[i] * t_0_xyy_z_yy[i];

        t_y_xyy_z_xz[i] = t_0_xyyy_z_xz[i] - rab_y[i] * t_0_xyy_z_xz[i];

        t_y_xyy_z_xy[i] = t_0_xyyy_z_xy[i] - rab_y[i] * t_0_xyy_z_xy[i];

        t_y_xyy_z_xx[i] = t_0_xyyy_z_xx[i] - rab_y[i] * t_0_xyy_z_xx[i];

        t_y_xyy_y_zz[i] = t_0_xyyy_y_zz[i] - rab_y[i] * t_0_xyy_y_zz[i];

        t_y_xyy_y_yz[i] = t_0_xyyy_y_yz[i] - rab_y[i] * t_0_xyy_y_yz[i];

        t_y_xyy_y_yy[i] = t_0_xyyy_y_yy[i] - rab_y[i] * t_0_xyy_y_yy[i];

        t_y_xyy_y_xz[i] = t_0_xyyy_y_xz[i] - rab_y[i] * t_0_xyy_y_xz[i];

        t_y_xyy_y_xy[i] = t_0_xyyy_y_xy[i] - rab_y[i] * t_0_xyy_y_xy[i];

        t_y_xyy_y_xx[i] = t_0_xyyy_y_xx[i] - rab_y[i] * t_0_xyy_y_xx[i];

        t_y_xyy_x_zz[i] = t_0_xyyy_x_zz[i] - rab_y[i] * t_0_xyy_x_zz[i];

        t_y_xyy_x_yz[i] = t_0_xyyy_x_yz[i] - rab_y[i] * t_0_xyy_x_yz[i];

        t_y_xyy_x_yy[i] = t_0_xyyy_x_yy[i] - rab_y[i] * t_0_xyy_x_yy[i];

        t_y_xyy_x_xz[i] = t_0_xyyy_x_xz[i] - rab_y[i] * t_0_xyy_x_xz[i];

        t_y_xyy_x_xy[i] = t_0_xyyy_x_xy[i] - rab_y[i] * t_0_xyy_x_xy[i];

        t_y_xyy_x_xx[i] = t_0_xyyy_x_xx[i] - rab_y[i] * t_0_xyy_x_xx[i];

        t_y_xxz_z_zz[i] = t_0_xxyz_z_zz[i] - rab_y[i] * t_0_xxz_z_zz[i];

        t_y_xxz_z_yz[i] = t_0_xxyz_z_yz[i] - rab_y[i] * t_0_xxz_z_yz[i];

        t_y_xxz_z_yy[i] = t_0_xxyz_z_yy[i] - rab_y[i] * t_0_xxz_z_yy[i];

        t_y_xxz_z_xz[i] = t_0_xxyz_z_xz[i] - rab_y[i] * t_0_xxz_z_xz[i];

        t_y_xxz_z_xy[i] = t_0_xxyz_z_xy[i] - rab_y[i] * t_0_xxz_z_xy[i];

        t_y_xxz_z_xx[i] = t_0_xxyz_z_xx[i] - rab_y[i] * t_0_xxz_z_xx[i];

        t_y_xxz_y_zz[i] = t_0_xxyz_y_zz[i] - rab_y[i] * t_0_xxz_y_zz[i];

        t_y_xxz_y_yz[i] = t_0_xxyz_y_yz[i] - rab_y[i] * t_0_xxz_y_yz[i];

        t_y_xxz_y_yy[i] = t_0_xxyz_y_yy[i] - rab_y[i] * t_0_xxz_y_yy[i];

        t_y_xxz_y_xz[i] = t_0_xxyz_y_xz[i] - rab_y[i] * t_0_xxz_y_xz[i];

        t_y_xxz_y_xy[i] = t_0_xxyz_y_xy[i] - rab_y[i] * t_0_xxz_y_xy[i];

        t_y_xxz_y_xx[i] = t_0_xxyz_y_xx[i] - rab_y[i] * t_0_xxz_y_xx[i];

        t_y_xxz_x_zz[i] = t_0_xxyz_x_zz[i] - rab_y[i] * t_0_xxz_x_zz[i];

        t_y_xxz_x_yz[i] = t_0_xxyz_x_yz[i] - rab_y[i] * t_0_xxz_x_yz[i];

        t_y_xxz_x_yy[i] = t_0_xxyz_x_yy[i] - rab_y[i] * t_0_xxz_x_yy[i];

        t_y_xxz_x_xz[i] = t_0_xxyz_x_xz[i] - rab_y[i] * t_0_xxz_x_xz[i];

        t_y_xxz_x_xy[i] = t_0_xxyz_x_xy[i] - rab_y[i] * t_0_xxz_x_xy[i];

        t_y_xxz_x_xx[i] = t_0_xxyz_x_xx[i] - rab_y[i] * t_0_xxz_x_xx[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_0_xxy_x_xx, t_0_xxy_x_xy, t_0_xxy_x_xz, t_0_xxy_x_yy,\
                           t_0_xxy_x_yz, t_0_xxy_x_zz, t_0_xxy_y_xx, t_0_xxy_y_xy, t_0_xxy_y_xz,\
                           t_0_xxy_y_yy, t_0_xxy_y_yz, t_0_xxy_y_zz, t_0_xxy_z_xx, t_0_xxy_z_xy,\
                           t_0_xxy_z_xz, t_0_xxy_z_yy, t_0_xxy_z_yz, t_0_xxy_z_zz, t_0_xxyy_x_xx,\
                           t_0_xxyy_x_xy, t_0_xxyy_x_xz, t_0_xxyy_x_yy, t_0_xxyy_x_yz,\
                           t_0_xxyy_x_zz, t_0_xxyy_y_xx, t_0_xxyy_y_xy, t_0_xxyy_y_xz,\
                           t_0_xxyy_y_yy, t_0_xxyy_y_yz, t_0_xxyy_y_zz, t_0_xxyy_z_xx,\
                           t_0_xxyy_z_xy, t_0_xxyy_z_xz, t_0_xxyy_z_yy, t_0_xxyy_z_yz,\
                           t_0_xxyy_z_zz, t_0_xzzz_x_xx, t_0_xzzz_x_xy, t_0_xzzz_x_xz,\
                           t_0_xzzz_x_yy, t_0_xzzz_x_yz, t_0_xzzz_x_zz, t_0_xzzz_y_xx,\
                           t_0_xzzz_y_xy, t_0_xzzz_y_xz, t_0_xzzz_y_yy, t_0_xzzz_y_yz,\
                           t_0_xzzz_y_zz, t_0_xzzz_z_xx, t_0_xzzz_z_xy, t_0_xzzz_z_xz,\
                           t_0_xzzz_z_yy, t_0_xzzz_z_yz, t_0_xzzz_z_zz, t_0_zzz_x_xx,\
                           t_0_zzz_x_xy, t_0_zzz_x_xz, t_0_zzz_x_yy, t_0_zzz_x_yz, t_0_zzz_x_zz,\
                           t_0_zzz_y_xx, t_0_zzz_y_xy, t_0_zzz_y_xz, t_0_zzz_y_yy, t_0_zzz_y_yz,\
                           t_0_zzz_y_zz, t_0_zzz_z_xx, t_0_zzz_z_xy, t_0_zzz_z_xz, t_0_zzz_z_yy,\
                           t_0_zzz_z_yz, t_0_zzz_z_zz, t_x_zzz_x_xx, t_x_zzz_x_xy, t_x_zzz_x_xz,\
                           t_x_zzz_x_yy, t_x_zzz_x_yz, t_x_zzz_x_zz, t_x_zzz_y_xx, t_x_zzz_y_xy,\
                           t_x_zzz_y_xz, t_x_zzz_y_yy, t_x_zzz_y_yz, t_x_zzz_y_zz, t_x_zzz_z_xx,\
                           t_x_zzz_z_xy, t_x_zzz_z_xz, t_x_zzz_z_yy, t_x_zzz_z_yz, t_x_zzz_z_zz,\
                           t_y_xxy_x_xx, t_y_xxy_x_xy, t_y_xxy_x_xz, t_y_xxy_x_yy, t_y_xxy_x_yz,\
                           t_y_xxy_x_zz, t_y_xxy_y_xx, t_y_xxy_y_xy, t_y_xxy_y_xz, t_y_xxy_y_yy,\
                           t_y_xxy_y_yz, t_y_xxy_y_zz, t_y_xxy_z_xx, t_y_xxy_z_xy, t_y_xxy_z_xz,\
                           t_y_xxy_z_yy, t_y_xxy_z_yz, t_y_xxy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xxy_z_zz[i] = t_0_xxyy_z_zz[i] - rab_y[i] * t_0_xxy_z_zz[i];

        t_y_xxy_z_yz[i] = t_0_xxyy_z_yz[i] - rab_y[i] * t_0_xxy_z_yz[i];

        t_y_xxy_z_yy[i] = t_0_xxyy_z_yy[i] - rab_y[i] * t_0_xxy_z_yy[i];

        t_y_xxy_z_xz[i] = t_0_xxyy_z_xz[i] - rab_y[i] * t_0_xxy_z_xz[i];

        t_y_xxy_z_xy[i] = t_0_xxyy_z_xy[i] - rab_y[i] * t_0_xxy_z_xy[i];

        t_y_xxy_z_xx[i] = t_0_xxyy_z_xx[i] - rab_y[i] * t_0_xxy_z_xx[i];

        t_y_xxy_y_zz[i] = t_0_xxyy_y_zz[i] - rab_y[i] * t_0_xxy_y_zz[i];

        t_y_xxy_y_yz[i] = t_0_xxyy_y_yz[i] - rab_y[i] * t_0_xxy_y_yz[i];

        t_y_xxy_y_yy[i] = t_0_xxyy_y_yy[i] - rab_y[i] * t_0_xxy_y_yy[i];

        t_y_xxy_y_xz[i] = t_0_xxyy_y_xz[i] - rab_y[i] * t_0_xxy_y_xz[i];

        t_y_xxy_y_xy[i] = t_0_xxyy_y_xy[i] - rab_y[i] * t_0_xxy_y_xy[i];

        t_y_xxy_y_xx[i] = t_0_xxyy_y_xx[i] - rab_y[i] * t_0_xxy_y_xx[i];

        t_y_xxy_x_zz[i] = t_0_xxyy_x_zz[i] - rab_y[i] * t_0_xxy_x_zz[i];

        t_y_xxy_x_yz[i] = t_0_xxyy_x_yz[i] - rab_y[i] * t_0_xxy_x_yz[i];

        t_y_xxy_x_yy[i] = t_0_xxyy_x_yy[i] - rab_y[i] * t_0_xxy_x_yy[i];

        t_y_xxy_x_xz[i] = t_0_xxyy_x_xz[i] - rab_y[i] * t_0_xxy_x_xz[i];

        t_y_xxy_x_xy[i] = t_0_xxyy_x_xy[i] - rab_y[i] * t_0_xxy_x_xy[i];

        t_y_xxy_x_xx[i] = t_0_xxyy_x_xx[i] - rab_y[i] * t_0_xxy_x_xx[i];

        t_x_zzz_z_zz[i] = t_0_xzzz_z_zz[i] - rab_x[i] * t_0_zzz_z_zz[i];

        t_x_zzz_z_yz[i] = t_0_xzzz_z_yz[i] - rab_x[i] * t_0_zzz_z_yz[i];

        t_x_zzz_z_yy[i] = t_0_xzzz_z_yy[i] - rab_x[i] * t_0_zzz_z_yy[i];

        t_x_zzz_z_xz[i] = t_0_xzzz_z_xz[i] - rab_x[i] * t_0_zzz_z_xz[i];

        t_x_zzz_z_xy[i] = t_0_xzzz_z_xy[i] - rab_x[i] * t_0_zzz_z_xy[i];

        t_x_zzz_z_xx[i] = t_0_xzzz_z_xx[i] - rab_x[i] * t_0_zzz_z_xx[i];

        t_x_zzz_y_zz[i] = t_0_xzzz_y_zz[i] - rab_x[i] * t_0_zzz_y_zz[i];

        t_x_zzz_y_yz[i] = t_0_xzzz_y_yz[i] - rab_x[i] * t_0_zzz_y_yz[i];

        t_x_zzz_y_yy[i] = t_0_xzzz_y_yy[i] - rab_x[i] * t_0_zzz_y_yy[i];

        t_x_zzz_y_xz[i] = t_0_xzzz_y_xz[i] - rab_x[i] * t_0_zzz_y_xz[i];

        t_x_zzz_y_xy[i] = t_0_xzzz_y_xy[i] - rab_x[i] * t_0_zzz_y_xy[i];

        t_x_zzz_y_xx[i] = t_0_xzzz_y_xx[i] - rab_x[i] * t_0_zzz_y_xx[i];

        t_x_zzz_x_zz[i] = t_0_xzzz_x_zz[i] - rab_x[i] * t_0_zzz_x_zz[i];

        t_x_zzz_x_yz[i] = t_0_xzzz_x_yz[i] - rab_x[i] * t_0_zzz_x_yz[i];

        t_x_zzz_x_yy[i] = t_0_xzzz_x_yy[i] - rab_x[i] * t_0_zzz_x_yy[i];

        t_x_zzz_x_xz[i] = t_0_xzzz_x_xz[i] - rab_x[i] * t_0_zzz_x_xz[i];

        t_x_zzz_x_xy[i] = t_0_xzzz_x_xy[i] - rab_x[i] * t_0_zzz_x_xy[i];

        t_x_zzz_x_xx[i] = t_0_xzzz_x_xx[i] - rab_x[i] * t_0_zzz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyyz_x_xx, t_0_xyyz_x_xy, t_0_xyyz_x_xz, t_0_xyyz_x_yy,\
                           t_0_xyyz_x_yz, t_0_xyyz_x_zz, t_0_xyyz_y_xx, t_0_xyyz_y_xy,\
                           t_0_xyyz_y_xz, t_0_xyyz_y_yy, t_0_xyyz_y_yz, t_0_xyyz_y_zz,\
                           t_0_xyyz_z_xx, t_0_xyyz_z_xy, t_0_xyyz_z_xz, t_0_xyyz_z_yy,\
                           t_0_xyyz_z_yz, t_0_xyyz_z_zz, t_0_xyzz_x_xx, t_0_xyzz_x_xy,\
                           t_0_xyzz_x_xz, t_0_xyzz_x_yy, t_0_xyzz_x_yz, t_0_xyzz_x_zz,\
                           t_0_xyzz_y_xx, t_0_xyzz_y_xy, t_0_xyzz_y_xz, t_0_xyzz_y_yy,\
                           t_0_xyzz_y_yz, t_0_xyzz_y_zz, t_0_xyzz_z_xx, t_0_xyzz_z_xy,\
                           t_0_xyzz_z_xz, t_0_xyzz_z_yy, t_0_xyzz_z_yz, t_0_xyzz_z_zz,\
                           t_0_yyz_x_xx, t_0_yyz_x_xy, t_0_yyz_x_xz, t_0_yyz_x_yy, t_0_yyz_x_yz,\
                           t_0_yyz_x_zz, t_0_yyz_y_xx, t_0_yyz_y_xy, t_0_yyz_y_xz, t_0_yyz_y_yy,\
                           t_0_yyz_y_yz, t_0_yyz_y_zz, t_0_yyz_z_xx, t_0_yyz_z_xy, t_0_yyz_z_xz,\
                           t_0_yyz_z_yy, t_0_yyz_z_yz, t_0_yyz_z_zz, t_0_yzz_x_xx, t_0_yzz_x_xy,\
                           t_0_yzz_x_xz, t_0_yzz_x_yy, t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_xx,\
                           t_0_yzz_y_xy, t_0_yzz_y_xz, t_0_yzz_y_yy, t_0_yzz_y_yz, t_0_yzz_y_zz,\
                           t_0_yzz_z_xx, t_0_yzz_z_xy, t_0_yzz_z_xz, t_0_yzz_z_yy, t_0_yzz_z_yz,\
                           t_0_yzz_z_zz, t_x_yyz_x_xx, t_x_yyz_x_xy, t_x_yyz_x_xz, t_x_yyz_x_yy,\
                           t_x_yyz_x_yz, t_x_yyz_x_zz, t_x_yyz_y_xx, t_x_yyz_y_xy, t_x_yyz_y_xz,\
                           t_x_yyz_y_yy, t_x_yyz_y_yz, t_x_yyz_y_zz, t_x_yyz_z_xx, t_x_yyz_z_xy,\
                           t_x_yyz_z_xz, t_x_yyz_z_yy, t_x_yyz_z_yz, t_x_yyz_z_zz, t_x_yzz_x_xx,\
                           t_x_yzz_x_xy, t_x_yzz_x_xz, t_x_yzz_x_yy, t_x_yzz_x_yz, t_x_yzz_x_zz,\
                           t_x_yzz_y_xx, t_x_yzz_y_xy, t_x_yzz_y_xz, t_x_yzz_y_yy, t_x_yzz_y_yz,\
                           t_x_yzz_y_zz, t_x_yzz_z_xx, t_x_yzz_z_xy, t_x_yzz_z_xz, t_x_yzz_z_yy,\
                           t_x_yzz_z_yz, t_x_yzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yzz_z_zz[i] = t_0_xyzz_z_zz[i] - rab_x[i] * t_0_yzz_z_zz[i];

        t_x_yzz_z_yz[i] = t_0_xyzz_z_yz[i] - rab_x[i] * t_0_yzz_z_yz[i];

        t_x_yzz_z_yy[i] = t_0_xyzz_z_yy[i] - rab_x[i] * t_0_yzz_z_yy[i];

        t_x_yzz_z_xz[i] = t_0_xyzz_z_xz[i] - rab_x[i] * t_0_yzz_z_xz[i];

        t_x_yzz_z_xy[i] = t_0_xyzz_z_xy[i] - rab_x[i] * t_0_yzz_z_xy[i];

        t_x_yzz_z_xx[i] = t_0_xyzz_z_xx[i] - rab_x[i] * t_0_yzz_z_xx[i];

        t_x_yzz_y_zz[i] = t_0_xyzz_y_zz[i] - rab_x[i] * t_0_yzz_y_zz[i];

        t_x_yzz_y_yz[i] = t_0_xyzz_y_yz[i] - rab_x[i] * t_0_yzz_y_yz[i];

        t_x_yzz_y_yy[i] = t_0_xyzz_y_yy[i] - rab_x[i] * t_0_yzz_y_yy[i];

        t_x_yzz_y_xz[i] = t_0_xyzz_y_xz[i] - rab_x[i] * t_0_yzz_y_xz[i];

        t_x_yzz_y_xy[i] = t_0_xyzz_y_xy[i] - rab_x[i] * t_0_yzz_y_xy[i];

        t_x_yzz_y_xx[i] = t_0_xyzz_y_xx[i] - rab_x[i] * t_0_yzz_y_xx[i];

        t_x_yzz_x_zz[i] = t_0_xyzz_x_zz[i] - rab_x[i] * t_0_yzz_x_zz[i];

        t_x_yzz_x_yz[i] = t_0_xyzz_x_yz[i] - rab_x[i] * t_0_yzz_x_yz[i];

        t_x_yzz_x_yy[i] = t_0_xyzz_x_yy[i] - rab_x[i] * t_0_yzz_x_yy[i];

        t_x_yzz_x_xz[i] = t_0_xyzz_x_xz[i] - rab_x[i] * t_0_yzz_x_xz[i];

        t_x_yzz_x_xy[i] = t_0_xyzz_x_xy[i] - rab_x[i] * t_0_yzz_x_xy[i];

        t_x_yzz_x_xx[i] = t_0_xyzz_x_xx[i] - rab_x[i] * t_0_yzz_x_xx[i];

        t_x_yyz_z_zz[i] = t_0_xyyz_z_zz[i] - rab_x[i] * t_0_yyz_z_zz[i];

        t_x_yyz_z_yz[i] = t_0_xyyz_z_yz[i] - rab_x[i] * t_0_yyz_z_yz[i];

        t_x_yyz_z_yy[i] = t_0_xyyz_z_yy[i] - rab_x[i] * t_0_yyz_z_yy[i];

        t_x_yyz_z_xz[i] = t_0_xyyz_z_xz[i] - rab_x[i] * t_0_yyz_z_xz[i];

        t_x_yyz_z_xy[i] = t_0_xyyz_z_xy[i] - rab_x[i] * t_0_yyz_z_xy[i];

        t_x_yyz_z_xx[i] = t_0_xyyz_z_xx[i] - rab_x[i] * t_0_yyz_z_xx[i];

        t_x_yyz_y_zz[i] = t_0_xyyz_y_zz[i] - rab_x[i] * t_0_yyz_y_zz[i];

        t_x_yyz_y_yz[i] = t_0_xyyz_y_yz[i] - rab_x[i] * t_0_yyz_y_yz[i];

        t_x_yyz_y_yy[i] = t_0_xyyz_y_yy[i] - rab_x[i] * t_0_yyz_y_yy[i];

        t_x_yyz_y_xz[i] = t_0_xyyz_y_xz[i] - rab_x[i] * t_0_yyz_y_xz[i];

        t_x_yyz_y_xy[i] = t_0_xyyz_y_xy[i] - rab_x[i] * t_0_yyz_y_xy[i];

        t_x_yyz_y_xx[i] = t_0_xyyz_y_xx[i] - rab_x[i] * t_0_yyz_y_xx[i];

        t_x_yyz_x_zz[i] = t_0_xyyz_x_zz[i] - rab_x[i] * t_0_yyz_x_zz[i];

        t_x_yyz_x_yz[i] = t_0_xyyz_x_yz[i] - rab_x[i] * t_0_yyz_x_yz[i];

        t_x_yyz_x_yy[i] = t_0_xyyz_x_yy[i] - rab_x[i] * t_0_yyz_x_yy[i];

        t_x_yyz_x_xz[i] = t_0_xyyz_x_xz[i] - rab_x[i] * t_0_yyz_x_xz[i];

        t_x_yyz_x_xy[i] = t_0_xyyz_x_xy[i] - rab_x[i] * t_0_yyz_x_xy[i];

        t_x_yyz_x_xx[i] = t_0_xyyz_x_xx[i] - rab_x[i] * t_0_yyz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxzz_x_xx, t_0_xxzz_x_xy, t_0_xxzz_x_xz, t_0_xxzz_x_yy,\
                           t_0_xxzz_x_yz, t_0_xxzz_x_zz, t_0_xxzz_y_xx, t_0_xxzz_y_xy,\
                           t_0_xxzz_y_xz, t_0_xxzz_y_yy, t_0_xxzz_y_yz, t_0_xxzz_y_zz,\
                           t_0_xxzz_z_xx, t_0_xxzz_z_xy, t_0_xxzz_z_xz, t_0_xxzz_z_yy,\
                           t_0_xxzz_z_yz, t_0_xxzz_z_zz, t_0_xyyy_x_xx, t_0_xyyy_x_xy,\
                           t_0_xyyy_x_xz, t_0_xyyy_x_yy, t_0_xyyy_x_yz, t_0_xyyy_x_zz,\
                           t_0_xyyy_y_xx, t_0_xyyy_y_xy, t_0_xyyy_y_xz, t_0_xyyy_y_yy,\
                           t_0_xyyy_y_yz, t_0_xyyy_y_zz, t_0_xyyy_z_xx, t_0_xyyy_z_xy,\
                           t_0_xyyy_z_xz, t_0_xyyy_z_yy, t_0_xyyy_z_yz, t_0_xyyy_z_zz,\
                           t_0_xzz_x_xx, t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy, t_0_xzz_x_yz,\
                           t_0_xzz_x_zz, t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz, t_0_xzz_y_yy,\
                           t_0_xzz_y_yz, t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy, t_0_xzz_z_xz,\
                           t_0_xzz_z_yy, t_0_xzz_z_yz, t_0_xzz_z_zz, t_0_yyy_x_xx, t_0_yyy_x_xy,\
                           t_0_yyy_x_xz, t_0_yyy_x_yy, t_0_yyy_x_yz, t_0_yyy_x_zz, t_0_yyy_y_xx,\
                           t_0_yyy_y_xy, t_0_yyy_y_xz, t_0_yyy_y_yy, t_0_yyy_y_yz, t_0_yyy_y_zz,\
                           t_0_yyy_z_xx, t_0_yyy_z_xy, t_0_yyy_z_xz, t_0_yyy_z_yy, t_0_yyy_z_yz,\
                           t_0_yyy_z_zz, t_x_xzz_x_xx, t_x_xzz_x_xy, t_x_xzz_x_xz, t_x_xzz_x_yy,\
                           t_x_xzz_x_yz, t_x_xzz_x_zz, t_x_xzz_y_xx, t_x_xzz_y_xy, t_x_xzz_y_xz,\
                           t_x_xzz_y_yy, t_x_xzz_y_yz, t_x_xzz_y_zz, t_x_xzz_z_xx, t_x_xzz_z_xy,\
                           t_x_xzz_z_xz, t_x_xzz_z_yy, t_x_xzz_z_yz, t_x_xzz_z_zz, t_x_yyy_x_xx,\
                           t_x_yyy_x_xy, t_x_yyy_x_xz, t_x_yyy_x_yy, t_x_yyy_x_yz, t_x_yyy_x_zz,\
                           t_x_yyy_y_xx, t_x_yyy_y_xy, t_x_yyy_y_xz, t_x_yyy_y_yy, t_x_yyy_y_yz,\
                           t_x_yyy_y_zz, t_x_yyy_z_xx, t_x_yyy_z_xy, t_x_yyy_z_xz, t_x_yyy_z_yy,\
                           t_x_yyy_z_yz, t_x_yyy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yyy_z_zz[i] = t_0_xyyy_z_zz[i] - rab_x[i] * t_0_yyy_z_zz[i];

        t_x_yyy_z_yz[i] = t_0_xyyy_z_yz[i] - rab_x[i] * t_0_yyy_z_yz[i];

        t_x_yyy_z_yy[i] = t_0_xyyy_z_yy[i] - rab_x[i] * t_0_yyy_z_yy[i];

        t_x_yyy_z_xz[i] = t_0_xyyy_z_xz[i] - rab_x[i] * t_0_yyy_z_xz[i];

        t_x_yyy_z_xy[i] = t_0_xyyy_z_xy[i] - rab_x[i] * t_0_yyy_z_xy[i];

        t_x_yyy_z_xx[i] = t_0_xyyy_z_xx[i] - rab_x[i] * t_0_yyy_z_xx[i];

        t_x_yyy_y_zz[i] = t_0_xyyy_y_zz[i] - rab_x[i] * t_0_yyy_y_zz[i];

        t_x_yyy_y_yz[i] = t_0_xyyy_y_yz[i] - rab_x[i] * t_0_yyy_y_yz[i];

        t_x_yyy_y_yy[i] = t_0_xyyy_y_yy[i] - rab_x[i] * t_0_yyy_y_yy[i];

        t_x_yyy_y_xz[i] = t_0_xyyy_y_xz[i] - rab_x[i] * t_0_yyy_y_xz[i];

        t_x_yyy_y_xy[i] = t_0_xyyy_y_xy[i] - rab_x[i] * t_0_yyy_y_xy[i];

        t_x_yyy_y_xx[i] = t_0_xyyy_y_xx[i] - rab_x[i] * t_0_yyy_y_xx[i];

        t_x_yyy_x_zz[i] = t_0_xyyy_x_zz[i] - rab_x[i] * t_0_yyy_x_zz[i];

        t_x_yyy_x_yz[i] = t_0_xyyy_x_yz[i] - rab_x[i] * t_0_yyy_x_yz[i];

        t_x_yyy_x_yy[i] = t_0_xyyy_x_yy[i] - rab_x[i] * t_0_yyy_x_yy[i];

        t_x_yyy_x_xz[i] = t_0_xyyy_x_xz[i] - rab_x[i] * t_0_yyy_x_xz[i];

        t_x_yyy_x_xy[i] = t_0_xyyy_x_xy[i] - rab_x[i] * t_0_yyy_x_xy[i];

        t_x_yyy_x_xx[i] = t_0_xyyy_x_xx[i] - rab_x[i] * t_0_yyy_x_xx[i];

        t_x_xzz_z_zz[i] = t_0_xxzz_z_zz[i] - rab_x[i] * t_0_xzz_z_zz[i];

        t_x_xzz_z_yz[i] = t_0_xxzz_z_yz[i] - rab_x[i] * t_0_xzz_z_yz[i];

        t_x_xzz_z_yy[i] = t_0_xxzz_z_yy[i] - rab_x[i] * t_0_xzz_z_yy[i];

        t_x_xzz_z_xz[i] = t_0_xxzz_z_xz[i] - rab_x[i] * t_0_xzz_z_xz[i];

        t_x_xzz_z_xy[i] = t_0_xxzz_z_xy[i] - rab_x[i] * t_0_xzz_z_xy[i];

        t_x_xzz_z_xx[i] = t_0_xxzz_z_xx[i] - rab_x[i] * t_0_xzz_z_xx[i];

        t_x_xzz_y_zz[i] = t_0_xxzz_y_zz[i] - rab_x[i] * t_0_xzz_y_zz[i];

        t_x_xzz_y_yz[i] = t_0_xxzz_y_yz[i] - rab_x[i] * t_0_xzz_y_yz[i];

        t_x_xzz_y_yy[i] = t_0_xxzz_y_yy[i] - rab_x[i] * t_0_xzz_y_yy[i];

        t_x_xzz_y_xz[i] = t_0_xxzz_y_xz[i] - rab_x[i] * t_0_xzz_y_xz[i];

        t_x_xzz_y_xy[i] = t_0_xxzz_y_xy[i] - rab_x[i] * t_0_xzz_y_xy[i];

        t_x_xzz_y_xx[i] = t_0_xxzz_y_xx[i] - rab_x[i] * t_0_xzz_y_xx[i];

        t_x_xzz_x_zz[i] = t_0_xxzz_x_zz[i] - rab_x[i] * t_0_xzz_x_zz[i];

        t_x_xzz_x_yz[i] = t_0_xxzz_x_yz[i] - rab_x[i] * t_0_xzz_x_yz[i];

        t_x_xzz_x_yy[i] = t_0_xxzz_x_yy[i] - rab_x[i] * t_0_xzz_x_yy[i];

        t_x_xzz_x_xz[i] = t_0_xxzz_x_xz[i] - rab_x[i] * t_0_xzz_x_xz[i];

        t_x_xzz_x_xy[i] = t_0_xxzz_x_xy[i] - rab_x[i] * t_0_xzz_x_xy[i];

        t_x_xzz_x_xx[i] = t_0_xxzz_x_xx[i] - rab_x[i] * t_0_xzz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxyy_x_xx, t_0_xxyy_x_xy, t_0_xxyy_x_xz, t_0_xxyy_x_yy,\
                           t_0_xxyy_x_yz, t_0_xxyy_x_zz, t_0_xxyy_y_xx, t_0_xxyy_y_xy,\
                           t_0_xxyy_y_xz, t_0_xxyy_y_yy, t_0_xxyy_y_yz, t_0_xxyy_y_zz,\
                           t_0_xxyy_z_xx, t_0_xxyy_z_xy, t_0_xxyy_z_xz, t_0_xxyy_z_yy,\
                           t_0_xxyy_z_yz, t_0_xxyy_z_zz, t_0_xxyz_x_xx, t_0_xxyz_x_xy,\
                           t_0_xxyz_x_xz, t_0_xxyz_x_yy, t_0_xxyz_x_yz, t_0_xxyz_x_zz,\
                           t_0_xxyz_y_xx, t_0_xxyz_y_xy, t_0_xxyz_y_xz, t_0_xxyz_y_yy,\
                           t_0_xxyz_y_yz, t_0_xxyz_y_zz, t_0_xxyz_z_xx, t_0_xxyz_z_xy,\
                           t_0_xxyz_z_xz, t_0_xxyz_z_yy, t_0_xxyz_z_yz, t_0_xxyz_z_zz,\
                           t_0_xyy_x_xx, t_0_xyy_x_xy, t_0_xyy_x_xz, t_0_xyy_x_yy, t_0_xyy_x_yz,\
                           t_0_xyy_x_zz, t_0_xyy_y_xx, t_0_xyy_y_xy, t_0_xyy_y_xz, t_0_xyy_y_yy,\
                           t_0_xyy_y_yz, t_0_xyy_y_zz, t_0_xyy_z_xx, t_0_xyy_z_xy, t_0_xyy_z_xz,\
                           t_0_xyy_z_yy, t_0_xyy_z_yz, t_0_xyy_z_zz, t_0_xyz_x_xx, t_0_xyz_x_xy,\
                           t_0_xyz_x_xz, t_0_xyz_x_yy, t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx,\
                           t_0_xyz_y_xy, t_0_xyz_y_xz, t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz,\
                           t_0_xyz_z_xx, t_0_xyz_z_xy, t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz,\
                           t_0_xyz_z_zz, t_x_xyy_x_xx, t_x_xyy_x_xy, t_x_xyy_x_xz, t_x_xyy_x_yy,\
                           t_x_xyy_x_yz, t_x_xyy_x_zz, t_x_xyy_y_xx, t_x_xyy_y_xy, t_x_xyy_y_xz,\
                           t_x_xyy_y_yy, t_x_xyy_y_yz, t_x_xyy_y_zz, t_x_xyy_z_xx, t_x_xyy_z_xy,\
                           t_x_xyy_z_xz, t_x_xyy_z_yy, t_x_xyy_z_yz, t_x_xyy_z_zz, t_x_xyz_x_xx,\
                           t_x_xyz_x_xy, t_x_xyz_x_xz, t_x_xyz_x_yy, t_x_xyz_x_yz, t_x_xyz_x_zz,\
                           t_x_xyz_y_xx, t_x_xyz_y_xy, t_x_xyz_y_xz, t_x_xyz_y_yy, t_x_xyz_y_yz,\
                           t_x_xyz_y_zz, t_x_xyz_z_xx, t_x_xyz_z_xy, t_x_xyz_z_xz, t_x_xyz_z_yy,\
                           t_x_xyz_z_yz, t_x_xyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xyz_z_zz[i] = t_0_xxyz_z_zz[i] - rab_x[i] * t_0_xyz_z_zz[i];

        t_x_xyz_z_yz[i] = t_0_xxyz_z_yz[i] - rab_x[i] * t_0_xyz_z_yz[i];

        t_x_xyz_z_yy[i] = t_0_xxyz_z_yy[i] - rab_x[i] * t_0_xyz_z_yy[i];

        t_x_xyz_z_xz[i] = t_0_xxyz_z_xz[i] - rab_x[i] * t_0_xyz_z_xz[i];

        t_x_xyz_z_xy[i] = t_0_xxyz_z_xy[i] - rab_x[i] * t_0_xyz_z_xy[i];

        t_x_xyz_z_xx[i] = t_0_xxyz_z_xx[i] - rab_x[i] * t_0_xyz_z_xx[i];

        t_x_xyz_y_zz[i] = t_0_xxyz_y_zz[i] - rab_x[i] * t_0_xyz_y_zz[i];

        t_x_xyz_y_yz[i] = t_0_xxyz_y_yz[i] - rab_x[i] * t_0_xyz_y_yz[i];

        t_x_xyz_y_yy[i] = t_0_xxyz_y_yy[i] - rab_x[i] * t_0_xyz_y_yy[i];

        t_x_xyz_y_xz[i] = t_0_xxyz_y_xz[i] - rab_x[i] * t_0_xyz_y_xz[i];

        t_x_xyz_y_xy[i] = t_0_xxyz_y_xy[i] - rab_x[i] * t_0_xyz_y_xy[i];

        t_x_xyz_y_xx[i] = t_0_xxyz_y_xx[i] - rab_x[i] * t_0_xyz_y_xx[i];

        t_x_xyz_x_zz[i] = t_0_xxyz_x_zz[i] - rab_x[i] * t_0_xyz_x_zz[i];

        t_x_xyz_x_yz[i] = t_0_xxyz_x_yz[i] - rab_x[i] * t_0_xyz_x_yz[i];

        t_x_xyz_x_yy[i] = t_0_xxyz_x_yy[i] - rab_x[i] * t_0_xyz_x_yy[i];

        t_x_xyz_x_xz[i] = t_0_xxyz_x_xz[i] - rab_x[i] * t_0_xyz_x_xz[i];

        t_x_xyz_x_xy[i] = t_0_xxyz_x_xy[i] - rab_x[i] * t_0_xyz_x_xy[i];

        t_x_xyz_x_xx[i] = t_0_xxyz_x_xx[i] - rab_x[i] * t_0_xyz_x_xx[i];

        t_x_xyy_z_zz[i] = t_0_xxyy_z_zz[i] - rab_x[i] * t_0_xyy_z_zz[i];

        t_x_xyy_z_yz[i] = t_0_xxyy_z_yz[i] - rab_x[i] * t_0_xyy_z_yz[i];

        t_x_xyy_z_yy[i] = t_0_xxyy_z_yy[i] - rab_x[i] * t_0_xyy_z_yy[i];

        t_x_xyy_z_xz[i] = t_0_xxyy_z_xz[i] - rab_x[i] * t_0_xyy_z_xz[i];

        t_x_xyy_z_xy[i] = t_0_xxyy_z_xy[i] - rab_x[i] * t_0_xyy_z_xy[i];

        t_x_xyy_z_xx[i] = t_0_xxyy_z_xx[i] - rab_x[i] * t_0_xyy_z_xx[i];

        t_x_xyy_y_zz[i] = t_0_xxyy_y_zz[i] - rab_x[i] * t_0_xyy_y_zz[i];

        t_x_xyy_y_yz[i] = t_0_xxyy_y_yz[i] - rab_x[i] * t_0_xyy_y_yz[i];

        t_x_xyy_y_yy[i] = t_0_xxyy_y_yy[i] - rab_x[i] * t_0_xyy_y_yy[i];

        t_x_xyy_y_xz[i] = t_0_xxyy_y_xz[i] - rab_x[i] * t_0_xyy_y_xz[i];

        t_x_xyy_y_xy[i] = t_0_xxyy_y_xy[i] - rab_x[i] * t_0_xyy_y_xy[i];

        t_x_xyy_y_xx[i] = t_0_xxyy_y_xx[i] - rab_x[i] * t_0_xyy_y_xx[i];

        t_x_xyy_x_zz[i] = t_0_xxyy_x_zz[i] - rab_x[i] * t_0_xyy_x_zz[i];

        t_x_xyy_x_yz[i] = t_0_xxyy_x_yz[i] - rab_x[i] * t_0_xyy_x_yz[i];

        t_x_xyy_x_yy[i] = t_0_xxyy_x_yy[i] - rab_x[i] * t_0_xyy_x_yy[i];

        t_x_xyy_x_xz[i] = t_0_xxyy_x_xz[i] - rab_x[i] * t_0_xyy_x_xz[i];

        t_x_xyy_x_xy[i] = t_0_xxyy_x_xy[i] - rab_x[i] * t_0_xyy_x_xy[i];

        t_x_xyy_x_xx[i] = t_0_xxyy_x_xx[i] - rab_x[i] * t_0_xyy_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxxy_x_xx, t_0_xxxy_x_xy, t_0_xxxy_x_xz, t_0_xxxy_x_yy,\
                           t_0_xxxy_x_yz, t_0_xxxy_x_zz, t_0_xxxy_y_xx, t_0_xxxy_y_xy,\
                           t_0_xxxy_y_xz, t_0_xxxy_y_yy, t_0_xxxy_y_yz, t_0_xxxy_y_zz,\
                           t_0_xxxy_z_xx, t_0_xxxy_z_xy, t_0_xxxy_z_xz, t_0_xxxy_z_yy,\
                           t_0_xxxy_z_yz, t_0_xxxy_z_zz, t_0_xxxz_x_xx, t_0_xxxz_x_xy,\
                           t_0_xxxz_x_xz, t_0_xxxz_x_yy, t_0_xxxz_x_yz, t_0_xxxz_x_zz,\
                           t_0_xxxz_y_xx, t_0_xxxz_y_xy, t_0_xxxz_y_xz, t_0_xxxz_y_yy,\
                           t_0_xxxz_y_yz, t_0_xxxz_y_zz, t_0_xxxz_z_xx, t_0_xxxz_z_xy,\
                           t_0_xxxz_z_xz, t_0_xxxz_z_yy, t_0_xxxz_z_yz, t_0_xxxz_z_zz,\
                           t_0_xxy_x_xx, t_0_xxy_x_xy, t_0_xxy_x_xz, t_0_xxy_x_yy, t_0_xxy_x_yz,\
                           t_0_xxy_x_zz, t_0_xxy_y_xx, t_0_xxy_y_xy, t_0_xxy_y_xz, t_0_xxy_y_yy,\
                           t_0_xxy_y_yz, t_0_xxy_y_zz, t_0_xxy_z_xx, t_0_xxy_z_xy, t_0_xxy_z_xz,\
                           t_0_xxy_z_yy, t_0_xxy_z_yz, t_0_xxy_z_zz, t_0_xxz_x_xx, t_0_xxz_x_xy,\
                           t_0_xxz_x_xz, t_0_xxz_x_yy, t_0_xxz_x_yz, t_0_xxz_x_zz, t_0_xxz_y_xx,\
                           t_0_xxz_y_xy, t_0_xxz_y_xz, t_0_xxz_y_yy, t_0_xxz_y_yz, t_0_xxz_y_zz,\
                           t_0_xxz_z_xx, t_0_xxz_z_xy, t_0_xxz_z_xz, t_0_xxz_z_yy, t_0_xxz_z_yz,\
                           t_0_xxz_z_zz, t_x_xxy_x_xx, t_x_xxy_x_xy, t_x_xxy_x_xz, t_x_xxy_x_yy,\
                           t_x_xxy_x_yz, t_x_xxy_x_zz, t_x_xxy_y_xx, t_x_xxy_y_xy, t_x_xxy_y_xz,\
                           t_x_xxy_y_yy, t_x_xxy_y_yz, t_x_xxy_y_zz, t_x_xxy_z_xx, t_x_xxy_z_xy,\
                           t_x_xxy_z_xz, t_x_xxy_z_yy, t_x_xxy_z_yz, t_x_xxy_z_zz, t_x_xxz_x_xx,\
                           t_x_xxz_x_xy, t_x_xxz_x_xz, t_x_xxz_x_yy, t_x_xxz_x_yz, t_x_xxz_x_zz,\
                           t_x_xxz_y_xx, t_x_xxz_y_xy, t_x_xxz_y_xz, t_x_xxz_y_yy, t_x_xxz_y_yz,\
                           t_x_xxz_y_zz, t_x_xxz_z_xx, t_x_xxz_z_xy, t_x_xxz_z_xz, t_x_xxz_z_yy,\
                           t_x_xxz_z_yz, t_x_xxz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxz_z_zz[i] = t_0_xxxz_z_zz[i] - rab_x[i] * t_0_xxz_z_zz[i];

        t_x_xxz_z_yz[i] = t_0_xxxz_z_yz[i] - rab_x[i] * t_0_xxz_z_yz[i];

        t_x_xxz_z_yy[i] = t_0_xxxz_z_yy[i] - rab_x[i] * t_0_xxz_z_yy[i];

        t_x_xxz_z_xz[i] = t_0_xxxz_z_xz[i] - rab_x[i] * t_0_xxz_z_xz[i];

        t_x_xxz_z_xy[i] = t_0_xxxz_z_xy[i] - rab_x[i] * t_0_xxz_z_xy[i];

        t_x_xxz_z_xx[i] = t_0_xxxz_z_xx[i] - rab_x[i] * t_0_xxz_z_xx[i];

        t_x_xxz_y_zz[i] = t_0_xxxz_y_zz[i] - rab_x[i] * t_0_xxz_y_zz[i];

        t_x_xxz_y_yz[i] = t_0_xxxz_y_yz[i] - rab_x[i] * t_0_xxz_y_yz[i];

        t_x_xxz_y_yy[i] = t_0_xxxz_y_yy[i] - rab_x[i] * t_0_xxz_y_yy[i];

        t_x_xxz_y_xz[i] = t_0_xxxz_y_xz[i] - rab_x[i] * t_0_xxz_y_xz[i];

        t_x_xxz_y_xy[i] = t_0_xxxz_y_xy[i] - rab_x[i] * t_0_xxz_y_xy[i];

        t_x_xxz_y_xx[i] = t_0_xxxz_y_xx[i] - rab_x[i] * t_0_xxz_y_xx[i];

        t_x_xxz_x_zz[i] = t_0_xxxz_x_zz[i] - rab_x[i] * t_0_xxz_x_zz[i];

        t_x_xxz_x_yz[i] = t_0_xxxz_x_yz[i] - rab_x[i] * t_0_xxz_x_yz[i];

        t_x_xxz_x_yy[i] = t_0_xxxz_x_yy[i] - rab_x[i] * t_0_xxz_x_yy[i];

        t_x_xxz_x_xz[i] = t_0_xxxz_x_xz[i] - rab_x[i] * t_0_xxz_x_xz[i];

        t_x_xxz_x_xy[i] = t_0_xxxz_x_xy[i] - rab_x[i] * t_0_xxz_x_xy[i];

        t_x_xxz_x_xx[i] = t_0_xxxz_x_xx[i] - rab_x[i] * t_0_xxz_x_xx[i];

        t_x_xxy_z_zz[i] = t_0_xxxy_z_zz[i] - rab_x[i] * t_0_xxy_z_zz[i];

        t_x_xxy_z_yz[i] = t_0_xxxy_z_yz[i] - rab_x[i] * t_0_xxy_z_yz[i];

        t_x_xxy_z_yy[i] = t_0_xxxy_z_yy[i] - rab_x[i] * t_0_xxy_z_yy[i];

        t_x_xxy_z_xz[i] = t_0_xxxy_z_xz[i] - rab_x[i] * t_0_xxy_z_xz[i];

        t_x_xxy_z_xy[i] = t_0_xxxy_z_xy[i] - rab_x[i] * t_0_xxy_z_xy[i];

        t_x_xxy_z_xx[i] = t_0_xxxy_z_xx[i] - rab_x[i] * t_0_xxy_z_xx[i];

        t_x_xxy_y_zz[i] = t_0_xxxy_y_zz[i] - rab_x[i] * t_0_xxy_y_zz[i];

        t_x_xxy_y_yz[i] = t_0_xxxy_y_yz[i] - rab_x[i] * t_0_xxy_y_yz[i];

        t_x_xxy_y_yy[i] = t_0_xxxy_y_yy[i] - rab_x[i] * t_0_xxy_y_yy[i];

        t_x_xxy_y_xz[i] = t_0_xxxy_y_xz[i] - rab_x[i] * t_0_xxy_y_xz[i];

        t_x_xxy_y_xy[i] = t_0_xxxy_y_xy[i] - rab_x[i] * t_0_xxy_y_xy[i];

        t_x_xxy_y_xx[i] = t_0_xxxy_y_xx[i] - rab_x[i] * t_0_xxy_y_xx[i];

        t_x_xxy_x_zz[i] = t_0_xxxy_x_zz[i] - rab_x[i] * t_0_xxy_x_zz[i];

        t_x_xxy_x_yz[i] = t_0_xxxy_x_yz[i] - rab_x[i] * t_0_xxy_x_yz[i];

        t_x_xxy_x_yy[i] = t_0_xxxy_x_yy[i] - rab_x[i] * t_0_xxy_x_yy[i];

        t_x_xxy_x_xz[i] = t_0_xxxy_x_xz[i] - rab_x[i] * t_0_xxy_x_xz[i];

        t_x_xxy_x_xy[i] = t_0_xxxy_x_xy[i] - rab_x[i] * t_0_xxy_x_xy[i];

        t_x_xxy_x_xx[i] = t_0_xxxy_x_xx[i] - rab_x[i] * t_0_xxy_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxx_x_xx, t_0_xxx_x_xy, t_0_xxx_x_xz, t_0_xxx_x_yy,\
                           t_0_xxx_x_yz, t_0_xxx_x_zz, t_0_xxx_y_xx, t_0_xxx_y_xy, t_0_xxx_y_xz,\
                           t_0_xxx_y_yy, t_0_xxx_y_yz, t_0_xxx_y_zz, t_0_xxx_z_xx, t_0_xxx_z_xy,\
                           t_0_xxx_z_xz, t_0_xxx_z_yy, t_0_xxx_z_yz, t_0_xxx_z_zz, t_0_xxxx_x_xx,\
                           t_0_xxxx_x_xy, t_0_xxxx_x_xz, t_0_xxxx_x_yy, t_0_xxxx_x_yz,\
                           t_0_xxxx_x_zz, t_0_xxxx_y_xx, t_0_xxxx_y_xy, t_0_xxxx_y_xz,\
                           t_0_xxxx_y_yy, t_0_xxxx_y_yz, t_0_xxxx_y_zz, t_0_xxxx_z_xx,\
                           t_0_xxxx_z_xy, t_0_xxxx_z_xz, t_0_xxxx_z_yy, t_0_xxxx_z_yz,\
                           t_0_xxxx_z_zz, t_x_xxx_x_xx, t_x_xxx_x_xy, t_x_xxx_x_xz,\
                           t_x_xxx_x_yy, t_x_xxx_x_yz, t_x_xxx_x_zz, t_x_xxx_y_xx, t_x_xxx_y_xy,\
                           t_x_xxx_y_xz, t_x_xxx_y_yy, t_x_xxx_y_yz, t_x_xxx_y_zz, t_x_xxx_z_xx,\
                           t_x_xxx_z_xy, t_x_xxx_z_xz, t_x_xxx_z_yy, t_x_xxx_z_yz, t_x_xxx_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxx_z_zz[i] = t_0_xxxx_z_zz[i] - rab_x[i] * t_0_xxx_z_zz[i];

        t_x_xxx_z_yz[i] = t_0_xxxx_z_yz[i] - rab_x[i] * t_0_xxx_z_yz[i];

        t_x_xxx_z_yy[i] = t_0_xxxx_z_yy[i] - rab_x[i] * t_0_xxx_z_yy[i];

        t_x_xxx_z_xz[i] = t_0_xxxx_z_xz[i] - rab_x[i] * t_0_xxx_z_xz[i];

        t_x_xxx_z_xy[i] = t_0_xxxx_z_xy[i] - rab_x[i] * t_0_xxx_z_xy[i];

        t_x_xxx_z_xx[i] = t_0_xxxx_z_xx[i] - rab_x[i] * t_0_xxx_z_xx[i];

        t_x_xxx_y_zz[i] = t_0_xxxx_y_zz[i] - rab_x[i] * t_0_xxx_y_zz[i];

        t_x_xxx_y_yz[i] = t_0_xxxx_y_yz[i] - rab_x[i] * t_0_xxx_y_yz[i];

        t_x_xxx_y_yy[i] = t_0_xxxx_y_yy[i] - rab_x[i] * t_0_xxx_y_yy[i];

        t_x_xxx_y_xz[i] = t_0_xxxx_y_xz[i] - rab_x[i] * t_0_xxx_y_xz[i];

        t_x_xxx_y_xy[i] = t_0_xxxx_y_xy[i] - rab_x[i] * t_0_xxx_y_xy[i];

        t_x_xxx_y_xx[i] = t_0_xxxx_y_xx[i] - rab_x[i] * t_0_xxx_y_xx[i];

        t_x_xxx_x_zz[i] = t_0_xxxx_x_zz[i] - rab_x[i] * t_0_xxx_x_zz[i];

        t_x_xxx_x_yz[i] = t_0_xxxx_x_yz[i] - rab_x[i] * t_0_xxx_x_yz[i];

        t_x_xxx_x_yy[i] = t_0_xxxx_x_yy[i] - rab_x[i] * t_0_xxx_x_yy[i];

        t_x_xxx_x_xz[i] = t_0_xxxx_x_xz[i] - rab_x[i] * t_0_xxx_x_xz[i];

        t_x_xxx_x_xy[i] = t_0_xxxx_x_xy[i] - rab_x[i] * t_0_xxx_x_xy[i];

        t_x_xxx_x_xx[i] = t_0_xxxx_x_xx[i] - rab_x[i] * t_0_xxx_x_xx[i];
    }
}


} // derirec namespace
