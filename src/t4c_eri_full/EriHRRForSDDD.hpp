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
compHostHRRForSDDD_V0(      BufferHostXY<T>&      intsBufferSDDD,
                      const BufferHostX<int32_t>& intsIndexesSDDD,
                      const BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostXY<T>&      intsBufferSDPF,
                      const BufferHostX<int32_t>& intsIndexesSDPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SDDD) integral components

    t_0_zz_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(0));

    t_0_zz_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(1));

    t_0_zz_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(2));

    t_0_zz_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(3));

    t_0_zz_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(4));

    t_0_zz_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(5));

    t_0_zz_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(6));

    t_0_zz_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(7));

    t_0_zz_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(8));

    t_0_zz_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(9));

    t_0_zz_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(10));

    t_0_zz_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(11));

    t_0_zz_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(12));

    t_0_zz_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(13));

    t_0_zz_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(14));

    t_0_zz_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(15));

    t_0_zz_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(16));

    t_0_zz_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(17));

    t_0_zz_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(18));

    t_0_zz_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(19));

    t_0_zz_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(20));

    t_0_zz_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(21));

    t_0_zz_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(22));

    t_0_zz_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(23));

    t_0_zz_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(24));

    t_0_zz_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(25));

    t_0_zz_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(26));

    t_0_zz_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(27));

    t_0_zz_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(28));

    t_0_zz_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(29));

    t_0_zz_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(30));

    t_0_zz_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(31));

    t_0_zz_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(32));

    t_0_zz_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(33));

    t_0_zz_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(34));

    t_0_zz_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(35));

    t_0_yz_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(36));

    t_0_yz_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(37));

    t_0_yz_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(38));

    t_0_yz_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(39));

    t_0_yz_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(40));

    t_0_yz_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(41));

    t_0_yz_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(42));

    t_0_yz_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(43));

    t_0_yz_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(44));

    t_0_yz_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(45));

    t_0_yz_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(46));

    t_0_yz_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(47));

    t_0_yz_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(48));

    t_0_yz_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(49));

    t_0_yz_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(50));

    t_0_yz_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(51));

    t_0_yz_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(52));

    t_0_yz_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(53));

    t_0_yz_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(54));

    t_0_yz_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(55));

    t_0_yz_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(56));

    t_0_yz_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(57));

    t_0_yz_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(58));

    t_0_yz_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(59));

    t_0_yz_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(60));

    t_0_yz_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(61));

    t_0_yz_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(62));

    t_0_yz_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(63));

    t_0_yz_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(64));

    t_0_yz_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(65));

    t_0_yz_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(66));

    t_0_yz_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(67));

    t_0_yz_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(68));

    t_0_yz_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(69));

    t_0_yz_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(70));

    t_0_yz_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(71));

    t_0_yy_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(72));

    t_0_yy_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(73));

    t_0_yy_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(74));

    t_0_yy_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(75));

    t_0_yy_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(76));

    t_0_yy_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(77));

    t_0_yy_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(78));

    t_0_yy_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(79));

    t_0_yy_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(80));

    t_0_yy_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(81));

    t_0_yy_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(82));

    t_0_yy_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(83));

    t_0_yy_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(84));

    t_0_yy_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(85));

    t_0_yy_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(86));

    t_0_yy_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(87));

    t_0_yy_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(88));

    t_0_yy_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(89));

    t_0_yy_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(90));

    t_0_yy_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(91));

    t_0_yy_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(92));

    t_0_yy_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(93));

    t_0_yy_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(94));

    t_0_yy_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(95));

    t_0_yy_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(96));

    t_0_yy_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(97));

    t_0_yy_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(98));

    t_0_yy_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(99));

    t_0_yy_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(100));

    t_0_yy_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(101));

    t_0_yy_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(102));

    t_0_yy_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(103));

    t_0_yy_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(104));

    t_0_yy_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(105));

    t_0_yy_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(106));

    t_0_yy_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(107));

    t_0_xz_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(108));

    t_0_xz_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(109));

    t_0_xz_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(110));

    t_0_xz_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(111));

    t_0_xz_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(112));

    t_0_xz_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(113));

    t_0_xz_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(114));

    t_0_xz_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(115));

    t_0_xz_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(116));

    t_0_xz_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(117));

    t_0_xz_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(118));

    t_0_xz_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(119));

    t_0_xz_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(120));

    t_0_xz_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(121));

    t_0_xz_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(122));

    t_0_xz_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(123));

    t_0_xz_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(124));

    t_0_xz_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(125));

    t_0_xz_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(126));

    t_0_xz_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(127));

    t_0_xz_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(128));

    t_0_xz_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(129));

    t_0_xz_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(130));

    t_0_xz_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(131));

    t_0_xz_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(132));

    t_0_xz_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(133));

    t_0_xz_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(134));

    t_0_xz_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(135));

    t_0_xz_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(136));

    t_0_xz_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(137));

    t_0_xz_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(138));

    t_0_xz_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(139));

    t_0_xz_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(140));

    t_0_xz_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(141));

    t_0_xz_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(142));

    t_0_xz_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(143));

    t_0_xy_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(144));

    t_0_xy_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(145));

    t_0_xy_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(146));

    t_0_xy_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(147));

    t_0_xy_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(148));

    t_0_xy_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(149));

    t_0_xy_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(150));

    t_0_xy_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(151));

    t_0_xy_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(152));

    t_0_xy_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(153));

    t_0_xy_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(154));

    t_0_xy_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(155));

    t_0_xy_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(156));

    t_0_xy_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(157));

    t_0_xy_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(158));

    t_0_xy_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(159));

    t_0_xy_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(160));

    t_0_xy_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(161));

    t_0_xy_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(162));

    t_0_xy_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(163));

    t_0_xy_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(164));

    t_0_xy_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(165));

    t_0_xy_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(166));

    t_0_xy_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(167));

    t_0_xy_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(168));

    t_0_xy_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(169));

    t_0_xy_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(170));

    t_0_xy_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(171));

    t_0_xy_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(172));

    t_0_xy_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(173));

    t_0_xy_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(174));

    t_0_xy_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(175));

    t_0_xy_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(176));

    t_0_xy_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(177));

    t_0_xy_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(178));

    t_0_xy_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(179));

    t_0_xx_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(180));

    t_0_xx_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(181));

    t_0_xx_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(182));

    t_0_xx_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(183));

    t_0_xx_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(184));

    t_0_xx_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(185));

    t_0_xx_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(186));

    t_0_xx_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(187));

    t_0_xx_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(188));

    t_0_xx_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(189));

    t_0_xx_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(190));

    t_0_xx_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(191));

    t_0_xx_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(192));

    t_0_xx_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(193));

    t_0_xx_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(194));

    t_0_xx_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(195));

    t_0_xx_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(196));

    t_0_xx_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(197));

    t_0_xx_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(198));

    t_0_xx_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(199));

    t_0_xx_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(200));

    t_0_xx_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(201));

    t_0_xx_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(202));

    t_0_xx_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(203));

    t_0_xx_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(204));

    t_0_xx_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(205));

    t_0_xx_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(206));

    t_0_xx_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(207));

    t_0_xx_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(208));

    t_0_xx_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(209));

    t_0_xx_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(210));

    t_0_xx_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(211));

    t_0_xx_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(212));

    t_0_xx_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(213));

    t_0_xx_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(214));

    t_0_xx_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(215));

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

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_zz_x_xx, t_0_zz_x_xxx, t_0_zz_x_xxy,\
                           t_0_zz_x_xxz, t_0_zz_x_xy, t_0_zz_x_xyy, t_0_zz_x_xyz, t_0_zz_x_xz,\
                           t_0_zz_x_xzz, t_0_zz_x_yy, t_0_zz_x_yyy, t_0_zz_x_yyz, t_0_zz_x_yz,\
                           t_0_zz_x_yzz, t_0_zz_x_zz, t_0_zz_x_zzz, t_0_zz_xx_xx, t_0_zz_xx_xy,\
                           t_0_zz_xx_xz, t_0_zz_xx_yy, t_0_zz_xx_yz, t_0_zz_xx_zz, t_0_zz_xy_xx,\
                           t_0_zz_xy_xy, t_0_zz_xy_xz, t_0_zz_xy_yy, t_0_zz_xy_yz, t_0_zz_xy_zz,\
                           t_0_zz_xz_xx, t_0_zz_xz_xy, t_0_zz_xz_xz, t_0_zz_xz_yy, t_0_zz_xz_yz,\
                           t_0_zz_xz_zz, t_0_zz_y_xx, t_0_zz_y_xxy, t_0_zz_y_xxz, t_0_zz_y_xy,\
                           t_0_zz_y_xyy, t_0_zz_y_xyz, t_0_zz_y_xz, t_0_zz_y_xzz, t_0_zz_y_yy,\
                           t_0_zz_y_yyy, t_0_zz_y_yyz, t_0_zz_y_yz, t_0_zz_y_yzz, t_0_zz_y_zz,\
                           t_0_zz_y_zzz, t_0_zz_yy_xx, t_0_zz_yy_xy, t_0_zz_yy_xz, t_0_zz_yy_yy,\
                           t_0_zz_yy_yz, t_0_zz_yy_zz, t_0_zz_yz_xx, t_0_zz_yz_xy, t_0_zz_yz_xz,\
                           t_0_zz_yz_yy, t_0_zz_yz_yz, t_0_zz_yz_zz, t_0_zz_z_xx, t_0_zz_z_xxz,\
                           t_0_zz_z_xy, t_0_zz_z_xyz, t_0_zz_z_xz, t_0_zz_z_xzz, t_0_zz_z_yy,\
                           t_0_zz_z_yyz, t_0_zz_z_yz, t_0_zz_z_yzz, t_0_zz_z_zz, t_0_zz_z_zzz,\
                           t_0_zz_zz_xx, t_0_zz_zz_xy, t_0_zz_zz_xz, t_0_zz_zz_yy, t_0_zz_zz_yz,\
                           t_0_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_zz_zz[i] = t_0_zz_z_zzz[i] - rcd_z[i] * t_0_zz_z_zz[i];

        t_0_zz_zz_yz[i] = t_0_zz_z_yzz[i] - rcd_z[i] * t_0_zz_z_yz[i];

        t_0_zz_zz_yy[i] = t_0_zz_z_yyz[i] - rcd_z[i] * t_0_zz_z_yy[i];

        t_0_zz_zz_xz[i] = t_0_zz_z_xzz[i] - rcd_z[i] * t_0_zz_z_xz[i];

        t_0_zz_zz_xy[i] = t_0_zz_z_xyz[i] - rcd_z[i] * t_0_zz_z_xy[i];

        t_0_zz_zz_xx[i] = t_0_zz_z_xxz[i] - rcd_z[i] * t_0_zz_z_xx[i];

        t_0_zz_yz_zz[i] = t_0_zz_y_zzz[i] - rcd_z[i] * t_0_zz_y_zz[i];

        t_0_zz_yz_yz[i] = t_0_zz_y_yzz[i] - rcd_z[i] * t_0_zz_y_yz[i];

        t_0_zz_yz_yy[i] = t_0_zz_y_yyz[i] - rcd_z[i] * t_0_zz_y_yy[i];

        t_0_zz_yz_xz[i] = t_0_zz_y_xzz[i] - rcd_z[i] * t_0_zz_y_xz[i];

        t_0_zz_yz_xy[i] = t_0_zz_y_xyz[i] - rcd_z[i] * t_0_zz_y_xy[i];

        t_0_zz_yz_xx[i] = t_0_zz_y_xxz[i] - rcd_z[i] * t_0_zz_y_xx[i];

        t_0_zz_yy_zz[i] = t_0_zz_y_yzz[i] - rcd_y[i] * t_0_zz_y_zz[i];

        t_0_zz_yy_yz[i] = t_0_zz_y_yyz[i] - rcd_y[i] * t_0_zz_y_yz[i];

        t_0_zz_yy_yy[i] = t_0_zz_y_yyy[i] - rcd_y[i] * t_0_zz_y_yy[i];

        t_0_zz_yy_xz[i] = t_0_zz_y_xyz[i] - rcd_y[i] * t_0_zz_y_xz[i];

        t_0_zz_yy_xy[i] = t_0_zz_y_xyy[i] - rcd_y[i] * t_0_zz_y_xy[i];

        t_0_zz_yy_xx[i] = t_0_zz_y_xxy[i] - rcd_y[i] * t_0_zz_y_xx[i];

        t_0_zz_xz_zz[i] = t_0_zz_x_zzz[i] - rcd_z[i] * t_0_zz_x_zz[i];

        t_0_zz_xz_yz[i] = t_0_zz_x_yzz[i] - rcd_z[i] * t_0_zz_x_yz[i];

        t_0_zz_xz_yy[i] = t_0_zz_x_yyz[i] - rcd_z[i] * t_0_zz_x_yy[i];

        t_0_zz_xz_xz[i] = t_0_zz_x_xzz[i] - rcd_z[i] * t_0_zz_x_xz[i];

        t_0_zz_xz_xy[i] = t_0_zz_x_xyz[i] - rcd_z[i] * t_0_zz_x_xy[i];

        t_0_zz_xz_xx[i] = t_0_zz_x_xxz[i] - rcd_z[i] * t_0_zz_x_xx[i];

        t_0_zz_xy_zz[i] = t_0_zz_x_yzz[i] - rcd_y[i] * t_0_zz_x_zz[i];

        t_0_zz_xy_yz[i] = t_0_zz_x_yyz[i] - rcd_y[i] * t_0_zz_x_yz[i];

        t_0_zz_xy_yy[i] = t_0_zz_x_yyy[i] - rcd_y[i] * t_0_zz_x_yy[i];

        t_0_zz_xy_xz[i] = t_0_zz_x_xyz[i] - rcd_y[i] * t_0_zz_x_xz[i];

        t_0_zz_xy_xy[i] = t_0_zz_x_xyy[i] - rcd_y[i] * t_0_zz_x_xy[i];

        t_0_zz_xy_xx[i] = t_0_zz_x_xxy[i] - rcd_y[i] * t_0_zz_x_xx[i];

        t_0_zz_xx_zz[i] = t_0_zz_x_xzz[i] - rcd_x[i] * t_0_zz_x_zz[i];

        t_0_zz_xx_yz[i] = t_0_zz_x_xyz[i] - rcd_x[i] * t_0_zz_x_yz[i];

        t_0_zz_xx_yy[i] = t_0_zz_x_xyy[i] - rcd_x[i] * t_0_zz_x_yy[i];

        t_0_zz_xx_xz[i] = t_0_zz_x_xxz[i] - rcd_x[i] * t_0_zz_x_xz[i];

        t_0_zz_xx_xy[i] = t_0_zz_x_xxy[i] - rcd_x[i] * t_0_zz_x_xy[i];

        t_0_zz_xx_xx[i] = t_0_zz_x_xxx[i] - rcd_x[i] * t_0_zz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yz_x_xx, t_0_yz_x_xxx, t_0_yz_x_xxy,\
                           t_0_yz_x_xxz, t_0_yz_x_xy, t_0_yz_x_xyy, t_0_yz_x_xyz, t_0_yz_x_xz,\
                           t_0_yz_x_xzz, t_0_yz_x_yy, t_0_yz_x_yyy, t_0_yz_x_yyz, t_0_yz_x_yz,\
                           t_0_yz_x_yzz, t_0_yz_x_zz, t_0_yz_x_zzz, t_0_yz_xx_xx, t_0_yz_xx_xy,\
                           t_0_yz_xx_xz, t_0_yz_xx_yy, t_0_yz_xx_yz, t_0_yz_xx_zz, t_0_yz_xy_xx,\
                           t_0_yz_xy_xy, t_0_yz_xy_xz, t_0_yz_xy_yy, t_0_yz_xy_yz, t_0_yz_xy_zz,\
                           t_0_yz_xz_xx, t_0_yz_xz_xy, t_0_yz_xz_xz, t_0_yz_xz_yy, t_0_yz_xz_yz,\
                           t_0_yz_xz_zz, t_0_yz_y_xx, t_0_yz_y_xxy, t_0_yz_y_xxz, t_0_yz_y_xy,\
                           t_0_yz_y_xyy, t_0_yz_y_xyz, t_0_yz_y_xz, t_0_yz_y_xzz, t_0_yz_y_yy,\
                           t_0_yz_y_yyy, t_0_yz_y_yyz, t_0_yz_y_yz, t_0_yz_y_yzz, t_0_yz_y_zz,\
                           t_0_yz_y_zzz, t_0_yz_yy_xx, t_0_yz_yy_xy, t_0_yz_yy_xz, t_0_yz_yy_yy,\
                           t_0_yz_yy_yz, t_0_yz_yy_zz, t_0_yz_yz_xx, t_0_yz_yz_xy, t_0_yz_yz_xz,\
                           t_0_yz_yz_yy, t_0_yz_yz_yz, t_0_yz_yz_zz, t_0_yz_z_xx, t_0_yz_z_xxz,\
                           t_0_yz_z_xy, t_0_yz_z_xyz, t_0_yz_z_xz, t_0_yz_z_xzz, t_0_yz_z_yy,\
                           t_0_yz_z_yyz, t_0_yz_z_yz, t_0_yz_z_yzz, t_0_yz_z_zz, t_0_yz_z_zzz,\
                           t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz, t_0_yz_zz_yy, t_0_yz_zz_yz,\
                           t_0_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yz_zz_zz[i] = t_0_yz_z_zzz[i] - rcd_z[i] * t_0_yz_z_zz[i];

        t_0_yz_zz_yz[i] = t_0_yz_z_yzz[i] - rcd_z[i] * t_0_yz_z_yz[i];

        t_0_yz_zz_yy[i] = t_0_yz_z_yyz[i] - rcd_z[i] * t_0_yz_z_yy[i];

        t_0_yz_zz_xz[i] = t_0_yz_z_xzz[i] - rcd_z[i] * t_0_yz_z_xz[i];

        t_0_yz_zz_xy[i] = t_0_yz_z_xyz[i] - rcd_z[i] * t_0_yz_z_xy[i];

        t_0_yz_zz_xx[i] = t_0_yz_z_xxz[i] - rcd_z[i] * t_0_yz_z_xx[i];

        t_0_yz_yz_zz[i] = t_0_yz_y_zzz[i] - rcd_z[i] * t_0_yz_y_zz[i];

        t_0_yz_yz_yz[i] = t_0_yz_y_yzz[i] - rcd_z[i] * t_0_yz_y_yz[i];

        t_0_yz_yz_yy[i] = t_0_yz_y_yyz[i] - rcd_z[i] * t_0_yz_y_yy[i];

        t_0_yz_yz_xz[i] = t_0_yz_y_xzz[i] - rcd_z[i] * t_0_yz_y_xz[i];

        t_0_yz_yz_xy[i] = t_0_yz_y_xyz[i] - rcd_z[i] * t_0_yz_y_xy[i];

        t_0_yz_yz_xx[i] = t_0_yz_y_xxz[i] - rcd_z[i] * t_0_yz_y_xx[i];

        t_0_yz_yy_zz[i] = t_0_yz_y_yzz[i] - rcd_y[i] * t_0_yz_y_zz[i];

        t_0_yz_yy_yz[i] = t_0_yz_y_yyz[i] - rcd_y[i] * t_0_yz_y_yz[i];

        t_0_yz_yy_yy[i] = t_0_yz_y_yyy[i] - rcd_y[i] * t_0_yz_y_yy[i];

        t_0_yz_yy_xz[i] = t_0_yz_y_xyz[i] - rcd_y[i] * t_0_yz_y_xz[i];

        t_0_yz_yy_xy[i] = t_0_yz_y_xyy[i] - rcd_y[i] * t_0_yz_y_xy[i];

        t_0_yz_yy_xx[i] = t_0_yz_y_xxy[i] - rcd_y[i] * t_0_yz_y_xx[i];

        t_0_yz_xz_zz[i] = t_0_yz_x_zzz[i] - rcd_z[i] * t_0_yz_x_zz[i];

        t_0_yz_xz_yz[i] = t_0_yz_x_yzz[i] - rcd_z[i] * t_0_yz_x_yz[i];

        t_0_yz_xz_yy[i] = t_0_yz_x_yyz[i] - rcd_z[i] * t_0_yz_x_yy[i];

        t_0_yz_xz_xz[i] = t_0_yz_x_xzz[i] - rcd_z[i] * t_0_yz_x_xz[i];

        t_0_yz_xz_xy[i] = t_0_yz_x_xyz[i] - rcd_z[i] * t_0_yz_x_xy[i];

        t_0_yz_xz_xx[i] = t_0_yz_x_xxz[i] - rcd_z[i] * t_0_yz_x_xx[i];

        t_0_yz_xy_zz[i] = t_0_yz_x_yzz[i] - rcd_y[i] * t_0_yz_x_zz[i];

        t_0_yz_xy_yz[i] = t_0_yz_x_yyz[i] - rcd_y[i] * t_0_yz_x_yz[i];

        t_0_yz_xy_yy[i] = t_0_yz_x_yyy[i] - rcd_y[i] * t_0_yz_x_yy[i];

        t_0_yz_xy_xz[i] = t_0_yz_x_xyz[i] - rcd_y[i] * t_0_yz_x_xz[i];

        t_0_yz_xy_xy[i] = t_0_yz_x_xyy[i] - rcd_y[i] * t_0_yz_x_xy[i];

        t_0_yz_xy_xx[i] = t_0_yz_x_xxy[i] - rcd_y[i] * t_0_yz_x_xx[i];

        t_0_yz_xx_zz[i] = t_0_yz_x_xzz[i] - rcd_x[i] * t_0_yz_x_zz[i];

        t_0_yz_xx_yz[i] = t_0_yz_x_xyz[i] - rcd_x[i] * t_0_yz_x_yz[i];

        t_0_yz_xx_yy[i] = t_0_yz_x_xyy[i] - rcd_x[i] * t_0_yz_x_yy[i];

        t_0_yz_xx_xz[i] = t_0_yz_x_xxz[i] - rcd_x[i] * t_0_yz_x_xz[i];

        t_0_yz_xx_xy[i] = t_0_yz_x_xxy[i] - rcd_x[i] * t_0_yz_x_xy[i];

        t_0_yz_xx_xx[i] = t_0_yz_x_xxx[i] - rcd_x[i] * t_0_yz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yy_x_xx, t_0_yy_x_xxx, t_0_yy_x_xxy,\
                           t_0_yy_x_xxz, t_0_yy_x_xy, t_0_yy_x_xyy, t_0_yy_x_xyz, t_0_yy_x_xz,\
                           t_0_yy_x_xzz, t_0_yy_x_yy, t_0_yy_x_yyy, t_0_yy_x_yyz, t_0_yy_x_yz,\
                           t_0_yy_x_yzz, t_0_yy_x_zz, t_0_yy_x_zzz, t_0_yy_xx_xx, t_0_yy_xx_xy,\
                           t_0_yy_xx_xz, t_0_yy_xx_yy, t_0_yy_xx_yz, t_0_yy_xx_zz, t_0_yy_xy_xx,\
                           t_0_yy_xy_xy, t_0_yy_xy_xz, t_0_yy_xy_yy, t_0_yy_xy_yz, t_0_yy_xy_zz,\
                           t_0_yy_xz_xx, t_0_yy_xz_xy, t_0_yy_xz_xz, t_0_yy_xz_yy, t_0_yy_xz_yz,\
                           t_0_yy_xz_zz, t_0_yy_y_xx, t_0_yy_y_xxy, t_0_yy_y_xxz, t_0_yy_y_xy,\
                           t_0_yy_y_xyy, t_0_yy_y_xyz, t_0_yy_y_xz, t_0_yy_y_xzz, t_0_yy_y_yy,\
                           t_0_yy_y_yyy, t_0_yy_y_yyz, t_0_yy_y_yz, t_0_yy_y_yzz, t_0_yy_y_zz,\
                           t_0_yy_y_zzz, t_0_yy_yy_xx, t_0_yy_yy_xy, t_0_yy_yy_xz, t_0_yy_yy_yy,\
                           t_0_yy_yy_yz, t_0_yy_yy_zz, t_0_yy_yz_xx, t_0_yy_yz_xy, t_0_yy_yz_xz,\
                           t_0_yy_yz_yy, t_0_yy_yz_yz, t_0_yy_yz_zz, t_0_yy_z_xx, t_0_yy_z_xxz,\
                           t_0_yy_z_xy, t_0_yy_z_xyz, t_0_yy_z_xz, t_0_yy_z_xzz, t_0_yy_z_yy,\
                           t_0_yy_z_yyz, t_0_yy_z_yz, t_0_yy_z_yzz, t_0_yy_z_zz, t_0_yy_z_zzz,\
                           t_0_yy_zz_xx, t_0_yy_zz_xy, t_0_yy_zz_xz, t_0_yy_zz_yy, t_0_yy_zz_yz,\
                           t_0_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yy_zz_zz[i] = t_0_yy_z_zzz[i] - rcd_z[i] * t_0_yy_z_zz[i];

        t_0_yy_zz_yz[i] = t_0_yy_z_yzz[i] - rcd_z[i] * t_0_yy_z_yz[i];

        t_0_yy_zz_yy[i] = t_0_yy_z_yyz[i] - rcd_z[i] * t_0_yy_z_yy[i];

        t_0_yy_zz_xz[i] = t_0_yy_z_xzz[i] - rcd_z[i] * t_0_yy_z_xz[i];

        t_0_yy_zz_xy[i] = t_0_yy_z_xyz[i] - rcd_z[i] * t_0_yy_z_xy[i];

        t_0_yy_zz_xx[i] = t_0_yy_z_xxz[i] - rcd_z[i] * t_0_yy_z_xx[i];

        t_0_yy_yz_zz[i] = t_0_yy_y_zzz[i] - rcd_z[i] * t_0_yy_y_zz[i];

        t_0_yy_yz_yz[i] = t_0_yy_y_yzz[i] - rcd_z[i] * t_0_yy_y_yz[i];

        t_0_yy_yz_yy[i] = t_0_yy_y_yyz[i] - rcd_z[i] * t_0_yy_y_yy[i];

        t_0_yy_yz_xz[i] = t_0_yy_y_xzz[i] - rcd_z[i] * t_0_yy_y_xz[i];

        t_0_yy_yz_xy[i] = t_0_yy_y_xyz[i] - rcd_z[i] * t_0_yy_y_xy[i];

        t_0_yy_yz_xx[i] = t_0_yy_y_xxz[i] - rcd_z[i] * t_0_yy_y_xx[i];

        t_0_yy_yy_zz[i] = t_0_yy_y_yzz[i] - rcd_y[i] * t_0_yy_y_zz[i];

        t_0_yy_yy_yz[i] = t_0_yy_y_yyz[i] - rcd_y[i] * t_0_yy_y_yz[i];

        t_0_yy_yy_yy[i] = t_0_yy_y_yyy[i] - rcd_y[i] * t_0_yy_y_yy[i];

        t_0_yy_yy_xz[i] = t_0_yy_y_xyz[i] - rcd_y[i] * t_0_yy_y_xz[i];

        t_0_yy_yy_xy[i] = t_0_yy_y_xyy[i] - rcd_y[i] * t_0_yy_y_xy[i];

        t_0_yy_yy_xx[i] = t_0_yy_y_xxy[i] - rcd_y[i] * t_0_yy_y_xx[i];

        t_0_yy_xz_zz[i] = t_0_yy_x_zzz[i] - rcd_z[i] * t_0_yy_x_zz[i];

        t_0_yy_xz_yz[i] = t_0_yy_x_yzz[i] - rcd_z[i] * t_0_yy_x_yz[i];

        t_0_yy_xz_yy[i] = t_0_yy_x_yyz[i] - rcd_z[i] * t_0_yy_x_yy[i];

        t_0_yy_xz_xz[i] = t_0_yy_x_xzz[i] - rcd_z[i] * t_0_yy_x_xz[i];

        t_0_yy_xz_xy[i] = t_0_yy_x_xyz[i] - rcd_z[i] * t_0_yy_x_xy[i];

        t_0_yy_xz_xx[i] = t_0_yy_x_xxz[i] - rcd_z[i] * t_0_yy_x_xx[i];

        t_0_yy_xy_zz[i] = t_0_yy_x_yzz[i] - rcd_y[i] * t_0_yy_x_zz[i];

        t_0_yy_xy_yz[i] = t_0_yy_x_yyz[i] - rcd_y[i] * t_0_yy_x_yz[i];

        t_0_yy_xy_yy[i] = t_0_yy_x_yyy[i] - rcd_y[i] * t_0_yy_x_yy[i];

        t_0_yy_xy_xz[i] = t_0_yy_x_xyz[i] - rcd_y[i] * t_0_yy_x_xz[i];

        t_0_yy_xy_xy[i] = t_0_yy_x_xyy[i] - rcd_y[i] * t_0_yy_x_xy[i];

        t_0_yy_xy_xx[i] = t_0_yy_x_xxy[i] - rcd_y[i] * t_0_yy_x_xx[i];

        t_0_yy_xx_zz[i] = t_0_yy_x_xzz[i] - rcd_x[i] * t_0_yy_x_zz[i];

        t_0_yy_xx_yz[i] = t_0_yy_x_xyz[i] - rcd_x[i] * t_0_yy_x_yz[i];

        t_0_yy_xx_yy[i] = t_0_yy_x_xyy[i] - rcd_x[i] * t_0_yy_x_yy[i];

        t_0_yy_xx_xz[i] = t_0_yy_x_xxz[i] - rcd_x[i] * t_0_yy_x_xz[i];

        t_0_yy_xx_xy[i] = t_0_yy_x_xxy[i] - rcd_x[i] * t_0_yy_x_xy[i];

        t_0_yy_xx_xx[i] = t_0_yy_x_xxx[i] - rcd_x[i] * t_0_yy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xz_x_xx, t_0_xz_x_xxx, t_0_xz_x_xxy,\
                           t_0_xz_x_xxz, t_0_xz_x_xy, t_0_xz_x_xyy, t_0_xz_x_xyz, t_0_xz_x_xz,\
                           t_0_xz_x_xzz, t_0_xz_x_yy, t_0_xz_x_yyy, t_0_xz_x_yyz, t_0_xz_x_yz,\
                           t_0_xz_x_yzz, t_0_xz_x_zz, t_0_xz_x_zzz, t_0_xz_xx_xx, t_0_xz_xx_xy,\
                           t_0_xz_xx_xz, t_0_xz_xx_yy, t_0_xz_xx_yz, t_0_xz_xx_zz, t_0_xz_xy_xx,\
                           t_0_xz_xy_xy, t_0_xz_xy_xz, t_0_xz_xy_yy, t_0_xz_xy_yz, t_0_xz_xy_zz,\
                           t_0_xz_xz_xx, t_0_xz_xz_xy, t_0_xz_xz_xz, t_0_xz_xz_yy, t_0_xz_xz_yz,\
                           t_0_xz_xz_zz, t_0_xz_y_xx, t_0_xz_y_xxy, t_0_xz_y_xxz, t_0_xz_y_xy,\
                           t_0_xz_y_xyy, t_0_xz_y_xyz, t_0_xz_y_xz, t_0_xz_y_xzz, t_0_xz_y_yy,\
                           t_0_xz_y_yyy, t_0_xz_y_yyz, t_0_xz_y_yz, t_0_xz_y_yzz, t_0_xz_y_zz,\
                           t_0_xz_y_zzz, t_0_xz_yy_xx, t_0_xz_yy_xy, t_0_xz_yy_xz, t_0_xz_yy_yy,\
                           t_0_xz_yy_yz, t_0_xz_yy_zz, t_0_xz_yz_xx, t_0_xz_yz_xy, t_0_xz_yz_xz,\
                           t_0_xz_yz_yy, t_0_xz_yz_yz, t_0_xz_yz_zz, t_0_xz_z_xx, t_0_xz_z_xxz,\
                           t_0_xz_z_xy, t_0_xz_z_xyz, t_0_xz_z_xz, t_0_xz_z_xzz, t_0_xz_z_yy,\
                           t_0_xz_z_yyz, t_0_xz_z_yz, t_0_xz_z_yzz, t_0_xz_z_zz, t_0_xz_z_zzz,\
                           t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz, t_0_xz_zz_yy, t_0_xz_zz_yz,\
                           t_0_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xz_zz_zz[i] = t_0_xz_z_zzz[i] - rcd_z[i] * t_0_xz_z_zz[i];

        t_0_xz_zz_yz[i] = t_0_xz_z_yzz[i] - rcd_z[i] * t_0_xz_z_yz[i];

        t_0_xz_zz_yy[i] = t_0_xz_z_yyz[i] - rcd_z[i] * t_0_xz_z_yy[i];

        t_0_xz_zz_xz[i] = t_0_xz_z_xzz[i] - rcd_z[i] * t_0_xz_z_xz[i];

        t_0_xz_zz_xy[i] = t_0_xz_z_xyz[i] - rcd_z[i] * t_0_xz_z_xy[i];

        t_0_xz_zz_xx[i] = t_0_xz_z_xxz[i] - rcd_z[i] * t_0_xz_z_xx[i];

        t_0_xz_yz_zz[i] = t_0_xz_y_zzz[i] - rcd_z[i] * t_0_xz_y_zz[i];

        t_0_xz_yz_yz[i] = t_0_xz_y_yzz[i] - rcd_z[i] * t_0_xz_y_yz[i];

        t_0_xz_yz_yy[i] = t_0_xz_y_yyz[i] - rcd_z[i] * t_0_xz_y_yy[i];

        t_0_xz_yz_xz[i] = t_0_xz_y_xzz[i] - rcd_z[i] * t_0_xz_y_xz[i];

        t_0_xz_yz_xy[i] = t_0_xz_y_xyz[i] - rcd_z[i] * t_0_xz_y_xy[i];

        t_0_xz_yz_xx[i] = t_0_xz_y_xxz[i] - rcd_z[i] * t_0_xz_y_xx[i];

        t_0_xz_yy_zz[i] = t_0_xz_y_yzz[i] - rcd_y[i] * t_0_xz_y_zz[i];

        t_0_xz_yy_yz[i] = t_0_xz_y_yyz[i] - rcd_y[i] * t_0_xz_y_yz[i];

        t_0_xz_yy_yy[i] = t_0_xz_y_yyy[i] - rcd_y[i] * t_0_xz_y_yy[i];

        t_0_xz_yy_xz[i] = t_0_xz_y_xyz[i] - rcd_y[i] * t_0_xz_y_xz[i];

        t_0_xz_yy_xy[i] = t_0_xz_y_xyy[i] - rcd_y[i] * t_0_xz_y_xy[i];

        t_0_xz_yy_xx[i] = t_0_xz_y_xxy[i] - rcd_y[i] * t_0_xz_y_xx[i];

        t_0_xz_xz_zz[i] = t_0_xz_x_zzz[i] - rcd_z[i] * t_0_xz_x_zz[i];

        t_0_xz_xz_yz[i] = t_0_xz_x_yzz[i] - rcd_z[i] * t_0_xz_x_yz[i];

        t_0_xz_xz_yy[i] = t_0_xz_x_yyz[i] - rcd_z[i] * t_0_xz_x_yy[i];

        t_0_xz_xz_xz[i] = t_0_xz_x_xzz[i] - rcd_z[i] * t_0_xz_x_xz[i];

        t_0_xz_xz_xy[i] = t_0_xz_x_xyz[i] - rcd_z[i] * t_0_xz_x_xy[i];

        t_0_xz_xz_xx[i] = t_0_xz_x_xxz[i] - rcd_z[i] * t_0_xz_x_xx[i];

        t_0_xz_xy_zz[i] = t_0_xz_x_yzz[i] - rcd_y[i] * t_0_xz_x_zz[i];

        t_0_xz_xy_yz[i] = t_0_xz_x_yyz[i] - rcd_y[i] * t_0_xz_x_yz[i];

        t_0_xz_xy_yy[i] = t_0_xz_x_yyy[i] - rcd_y[i] * t_0_xz_x_yy[i];

        t_0_xz_xy_xz[i] = t_0_xz_x_xyz[i] - rcd_y[i] * t_0_xz_x_xz[i];

        t_0_xz_xy_xy[i] = t_0_xz_x_xyy[i] - rcd_y[i] * t_0_xz_x_xy[i];

        t_0_xz_xy_xx[i] = t_0_xz_x_xxy[i] - rcd_y[i] * t_0_xz_x_xx[i];

        t_0_xz_xx_zz[i] = t_0_xz_x_xzz[i] - rcd_x[i] * t_0_xz_x_zz[i];

        t_0_xz_xx_yz[i] = t_0_xz_x_xyz[i] - rcd_x[i] * t_0_xz_x_yz[i];

        t_0_xz_xx_yy[i] = t_0_xz_x_xyy[i] - rcd_x[i] * t_0_xz_x_yy[i];

        t_0_xz_xx_xz[i] = t_0_xz_x_xxz[i] - rcd_x[i] * t_0_xz_x_xz[i];

        t_0_xz_xx_xy[i] = t_0_xz_x_xxy[i] - rcd_x[i] * t_0_xz_x_xy[i];

        t_0_xz_xx_xx[i] = t_0_xz_x_xxx[i] - rcd_x[i] * t_0_xz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xy_x_xx, t_0_xy_x_xxx, t_0_xy_x_xxy,\
                           t_0_xy_x_xxz, t_0_xy_x_xy, t_0_xy_x_xyy, t_0_xy_x_xyz, t_0_xy_x_xz,\
                           t_0_xy_x_xzz, t_0_xy_x_yy, t_0_xy_x_yyy, t_0_xy_x_yyz, t_0_xy_x_yz,\
                           t_0_xy_x_yzz, t_0_xy_x_zz, t_0_xy_x_zzz, t_0_xy_xx_xx, t_0_xy_xx_xy,\
                           t_0_xy_xx_xz, t_0_xy_xx_yy, t_0_xy_xx_yz, t_0_xy_xx_zz, t_0_xy_xy_xx,\
                           t_0_xy_xy_xy, t_0_xy_xy_xz, t_0_xy_xy_yy, t_0_xy_xy_yz, t_0_xy_xy_zz,\
                           t_0_xy_xz_xx, t_0_xy_xz_xy, t_0_xy_xz_xz, t_0_xy_xz_yy, t_0_xy_xz_yz,\
                           t_0_xy_xz_zz, t_0_xy_y_xx, t_0_xy_y_xxy, t_0_xy_y_xxz, t_0_xy_y_xy,\
                           t_0_xy_y_xyy, t_0_xy_y_xyz, t_0_xy_y_xz, t_0_xy_y_xzz, t_0_xy_y_yy,\
                           t_0_xy_y_yyy, t_0_xy_y_yyz, t_0_xy_y_yz, t_0_xy_y_yzz, t_0_xy_y_zz,\
                           t_0_xy_y_zzz, t_0_xy_yy_xx, t_0_xy_yy_xy, t_0_xy_yy_xz, t_0_xy_yy_yy,\
                           t_0_xy_yy_yz, t_0_xy_yy_zz, t_0_xy_yz_xx, t_0_xy_yz_xy, t_0_xy_yz_xz,\
                           t_0_xy_yz_yy, t_0_xy_yz_yz, t_0_xy_yz_zz, t_0_xy_z_xx, t_0_xy_z_xxz,\
                           t_0_xy_z_xy, t_0_xy_z_xyz, t_0_xy_z_xz, t_0_xy_z_xzz, t_0_xy_z_yy,\
                           t_0_xy_z_yyz, t_0_xy_z_yz, t_0_xy_z_yzz, t_0_xy_z_zz, t_0_xy_z_zzz,\
                           t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz, t_0_xy_zz_yy, t_0_xy_zz_yz,\
                           t_0_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xy_zz_zz[i] = t_0_xy_z_zzz[i] - rcd_z[i] * t_0_xy_z_zz[i];

        t_0_xy_zz_yz[i] = t_0_xy_z_yzz[i] - rcd_z[i] * t_0_xy_z_yz[i];

        t_0_xy_zz_yy[i] = t_0_xy_z_yyz[i] - rcd_z[i] * t_0_xy_z_yy[i];

        t_0_xy_zz_xz[i] = t_0_xy_z_xzz[i] - rcd_z[i] * t_0_xy_z_xz[i];

        t_0_xy_zz_xy[i] = t_0_xy_z_xyz[i] - rcd_z[i] * t_0_xy_z_xy[i];

        t_0_xy_zz_xx[i] = t_0_xy_z_xxz[i] - rcd_z[i] * t_0_xy_z_xx[i];

        t_0_xy_yz_zz[i] = t_0_xy_y_zzz[i] - rcd_z[i] * t_0_xy_y_zz[i];

        t_0_xy_yz_yz[i] = t_0_xy_y_yzz[i] - rcd_z[i] * t_0_xy_y_yz[i];

        t_0_xy_yz_yy[i] = t_0_xy_y_yyz[i] - rcd_z[i] * t_0_xy_y_yy[i];

        t_0_xy_yz_xz[i] = t_0_xy_y_xzz[i] - rcd_z[i] * t_0_xy_y_xz[i];

        t_0_xy_yz_xy[i] = t_0_xy_y_xyz[i] - rcd_z[i] * t_0_xy_y_xy[i];

        t_0_xy_yz_xx[i] = t_0_xy_y_xxz[i] - rcd_z[i] * t_0_xy_y_xx[i];

        t_0_xy_yy_zz[i] = t_0_xy_y_yzz[i] - rcd_y[i] * t_0_xy_y_zz[i];

        t_0_xy_yy_yz[i] = t_0_xy_y_yyz[i] - rcd_y[i] * t_0_xy_y_yz[i];

        t_0_xy_yy_yy[i] = t_0_xy_y_yyy[i] - rcd_y[i] * t_0_xy_y_yy[i];

        t_0_xy_yy_xz[i] = t_0_xy_y_xyz[i] - rcd_y[i] * t_0_xy_y_xz[i];

        t_0_xy_yy_xy[i] = t_0_xy_y_xyy[i] - rcd_y[i] * t_0_xy_y_xy[i];

        t_0_xy_yy_xx[i] = t_0_xy_y_xxy[i] - rcd_y[i] * t_0_xy_y_xx[i];

        t_0_xy_xz_zz[i] = t_0_xy_x_zzz[i] - rcd_z[i] * t_0_xy_x_zz[i];

        t_0_xy_xz_yz[i] = t_0_xy_x_yzz[i] - rcd_z[i] * t_0_xy_x_yz[i];

        t_0_xy_xz_yy[i] = t_0_xy_x_yyz[i] - rcd_z[i] * t_0_xy_x_yy[i];

        t_0_xy_xz_xz[i] = t_0_xy_x_xzz[i] - rcd_z[i] * t_0_xy_x_xz[i];

        t_0_xy_xz_xy[i] = t_0_xy_x_xyz[i] - rcd_z[i] * t_0_xy_x_xy[i];

        t_0_xy_xz_xx[i] = t_0_xy_x_xxz[i] - rcd_z[i] * t_0_xy_x_xx[i];

        t_0_xy_xy_zz[i] = t_0_xy_x_yzz[i] - rcd_y[i] * t_0_xy_x_zz[i];

        t_0_xy_xy_yz[i] = t_0_xy_x_yyz[i] - rcd_y[i] * t_0_xy_x_yz[i];

        t_0_xy_xy_yy[i] = t_0_xy_x_yyy[i] - rcd_y[i] * t_0_xy_x_yy[i];

        t_0_xy_xy_xz[i] = t_0_xy_x_xyz[i] - rcd_y[i] * t_0_xy_x_xz[i];

        t_0_xy_xy_xy[i] = t_0_xy_x_xyy[i] - rcd_y[i] * t_0_xy_x_xy[i];

        t_0_xy_xy_xx[i] = t_0_xy_x_xxy[i] - rcd_y[i] * t_0_xy_x_xx[i];

        t_0_xy_xx_zz[i] = t_0_xy_x_xzz[i] - rcd_x[i] * t_0_xy_x_zz[i];

        t_0_xy_xx_yz[i] = t_0_xy_x_xyz[i] - rcd_x[i] * t_0_xy_x_yz[i];

        t_0_xy_xx_yy[i] = t_0_xy_x_xyy[i] - rcd_x[i] * t_0_xy_x_yy[i];

        t_0_xy_xx_xz[i] = t_0_xy_x_xxz[i] - rcd_x[i] * t_0_xy_x_xz[i];

        t_0_xy_xx_xy[i] = t_0_xy_x_xxy[i] - rcd_x[i] * t_0_xy_x_xy[i];

        t_0_xy_xx_xx[i] = t_0_xy_x_xxx[i] - rcd_x[i] * t_0_xy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_x_xx, t_0_xx_x_xxx, t_0_xx_x_xxy,\
                           t_0_xx_x_xxz, t_0_xx_x_xy, t_0_xx_x_xyy, t_0_xx_x_xyz, t_0_xx_x_xz,\
                           t_0_xx_x_xzz, t_0_xx_x_yy, t_0_xx_x_yyy, t_0_xx_x_yyz, t_0_xx_x_yz,\
                           t_0_xx_x_yzz, t_0_xx_x_zz, t_0_xx_x_zzz, t_0_xx_xx_xx, t_0_xx_xx_xy,\
                           t_0_xx_xx_xz, t_0_xx_xx_yy, t_0_xx_xx_yz, t_0_xx_xx_zz, t_0_xx_xy_xx,\
                           t_0_xx_xy_xy, t_0_xx_xy_xz, t_0_xx_xy_yy, t_0_xx_xy_yz, t_0_xx_xy_zz,\
                           t_0_xx_xz_xx, t_0_xx_xz_xy, t_0_xx_xz_xz, t_0_xx_xz_yy, t_0_xx_xz_yz,\
                           t_0_xx_xz_zz, t_0_xx_y_xx, t_0_xx_y_xxy, t_0_xx_y_xxz, t_0_xx_y_xy,\
                           t_0_xx_y_xyy, t_0_xx_y_xyz, t_0_xx_y_xz, t_0_xx_y_xzz, t_0_xx_y_yy,\
                           t_0_xx_y_yyy, t_0_xx_y_yyz, t_0_xx_y_yz, t_0_xx_y_yzz, t_0_xx_y_zz,\
                           t_0_xx_y_zzz, t_0_xx_yy_xx, t_0_xx_yy_xy, t_0_xx_yy_xz, t_0_xx_yy_yy,\
                           t_0_xx_yy_yz, t_0_xx_yy_zz, t_0_xx_yz_xx, t_0_xx_yz_xy, t_0_xx_yz_xz,\
                           t_0_xx_yz_yy, t_0_xx_yz_yz, t_0_xx_yz_zz, t_0_xx_z_xx, t_0_xx_z_xxz,\
                           t_0_xx_z_xy, t_0_xx_z_xyz, t_0_xx_z_xz, t_0_xx_z_xzz, t_0_xx_z_yy,\
                           t_0_xx_z_yyz, t_0_xx_z_yz, t_0_xx_z_yzz, t_0_xx_z_zz, t_0_xx_z_zzz,\
                           t_0_xx_zz_xx, t_0_xx_zz_xy, t_0_xx_zz_xz, t_0_xx_zz_yy, t_0_xx_zz_yz,\
                           t_0_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xx_zz_zz[i] = t_0_xx_z_zzz[i] - rcd_z[i] * t_0_xx_z_zz[i];

        t_0_xx_zz_yz[i] = t_0_xx_z_yzz[i] - rcd_z[i] * t_0_xx_z_yz[i];

        t_0_xx_zz_yy[i] = t_0_xx_z_yyz[i] - rcd_z[i] * t_0_xx_z_yy[i];

        t_0_xx_zz_xz[i] = t_0_xx_z_xzz[i] - rcd_z[i] * t_0_xx_z_xz[i];

        t_0_xx_zz_xy[i] = t_0_xx_z_xyz[i] - rcd_z[i] * t_0_xx_z_xy[i];

        t_0_xx_zz_xx[i] = t_0_xx_z_xxz[i] - rcd_z[i] * t_0_xx_z_xx[i];

        t_0_xx_yz_zz[i] = t_0_xx_y_zzz[i] - rcd_z[i] * t_0_xx_y_zz[i];

        t_0_xx_yz_yz[i] = t_0_xx_y_yzz[i] - rcd_z[i] * t_0_xx_y_yz[i];

        t_0_xx_yz_yy[i] = t_0_xx_y_yyz[i] - rcd_z[i] * t_0_xx_y_yy[i];

        t_0_xx_yz_xz[i] = t_0_xx_y_xzz[i] - rcd_z[i] * t_0_xx_y_xz[i];

        t_0_xx_yz_xy[i] = t_0_xx_y_xyz[i] - rcd_z[i] * t_0_xx_y_xy[i];

        t_0_xx_yz_xx[i] = t_0_xx_y_xxz[i] - rcd_z[i] * t_0_xx_y_xx[i];

        t_0_xx_yy_zz[i] = t_0_xx_y_yzz[i] - rcd_y[i] * t_0_xx_y_zz[i];

        t_0_xx_yy_yz[i] = t_0_xx_y_yyz[i] - rcd_y[i] * t_0_xx_y_yz[i];

        t_0_xx_yy_yy[i] = t_0_xx_y_yyy[i] - rcd_y[i] * t_0_xx_y_yy[i];

        t_0_xx_yy_xz[i] = t_0_xx_y_xyz[i] - rcd_y[i] * t_0_xx_y_xz[i];

        t_0_xx_yy_xy[i] = t_0_xx_y_xyy[i] - rcd_y[i] * t_0_xx_y_xy[i];

        t_0_xx_yy_xx[i] = t_0_xx_y_xxy[i] - rcd_y[i] * t_0_xx_y_xx[i];

        t_0_xx_xz_zz[i] = t_0_xx_x_zzz[i] - rcd_z[i] * t_0_xx_x_zz[i];

        t_0_xx_xz_yz[i] = t_0_xx_x_yzz[i] - rcd_z[i] * t_0_xx_x_yz[i];

        t_0_xx_xz_yy[i] = t_0_xx_x_yyz[i] - rcd_z[i] * t_0_xx_x_yy[i];

        t_0_xx_xz_xz[i] = t_0_xx_x_xzz[i] - rcd_z[i] * t_0_xx_x_xz[i];

        t_0_xx_xz_xy[i] = t_0_xx_x_xyz[i] - rcd_z[i] * t_0_xx_x_xy[i];

        t_0_xx_xz_xx[i] = t_0_xx_x_xxz[i] - rcd_z[i] * t_0_xx_x_xx[i];

        t_0_xx_xy_zz[i] = t_0_xx_x_yzz[i] - rcd_y[i] * t_0_xx_x_zz[i];

        t_0_xx_xy_yz[i] = t_0_xx_x_yyz[i] - rcd_y[i] * t_0_xx_x_yz[i];

        t_0_xx_xy_yy[i] = t_0_xx_x_yyy[i] - rcd_y[i] * t_0_xx_x_yy[i];

        t_0_xx_xy_xz[i] = t_0_xx_x_xyz[i] - rcd_y[i] * t_0_xx_x_xz[i];

        t_0_xx_xy_xy[i] = t_0_xx_x_xyy[i] - rcd_y[i] * t_0_xx_x_xy[i];

        t_0_xx_xy_xx[i] = t_0_xx_x_xxy[i] - rcd_y[i] * t_0_xx_x_xx[i];

        t_0_xx_xx_zz[i] = t_0_xx_x_xzz[i] - rcd_x[i] * t_0_xx_x_zz[i];

        t_0_xx_xx_yz[i] = t_0_xx_x_xyz[i] - rcd_x[i] * t_0_xx_x_yz[i];

        t_0_xx_xx_yy[i] = t_0_xx_x_xyy[i] - rcd_x[i] * t_0_xx_x_yy[i];

        t_0_xx_xx_xz[i] = t_0_xx_x_xxz[i] - rcd_x[i] * t_0_xx_x_xz[i];

        t_0_xx_xx_xy[i] = t_0_xx_x_xxy[i] - rcd_x[i] * t_0_xx_x_xy[i];

        t_0_xx_xx_xx[i] = t_0_xx_x_xxx[i] - rcd_x[i] * t_0_xx_x_xx[i];
    }
}


} // derirec namespace
