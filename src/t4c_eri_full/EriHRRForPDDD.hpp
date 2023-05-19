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
compHostHRRForPDDD_V0(      BufferHostXY<T>&      intsBufferPDDD,
                      const BufferHostX<int32_t>& intsIndexesPDDD,
                      const BufferHostXY<T>&      intsBufferSDDD,
                      const BufferHostX<int32_t>& intsIndexesSDDD,
                      const BufferHostXY<T>&      intsBufferSFDD,
                      const BufferHostX<int32_t>& intsIndexesSFDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PDDD) integral components

    t_z_zz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(0));

    t_z_zz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(1));

    t_z_zz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(2));

    t_z_zz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(3));

    t_z_zz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(4));

    t_z_zz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(5));

    t_z_zz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(6));

    t_z_zz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(7));

    t_z_zz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(8));

    t_z_zz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(9));

    t_z_zz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(10));

    t_z_zz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(11));

    t_z_zz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(12));

    t_z_zz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(13));

    t_z_zz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(14));

    t_z_zz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(15));

    t_z_zz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(16));

    t_z_zz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(17));

    t_z_zz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(18));

    t_z_zz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(19));

    t_z_zz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(20));

    t_z_zz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(21));

    t_z_zz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(22));

    t_z_zz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(23));

    t_z_zz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(24));

    t_z_zz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(25));

    t_z_zz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(26));

    t_z_zz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(27));

    t_z_zz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(28));

    t_z_zz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(29));

    t_z_zz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(30));

    t_z_zz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(31));

    t_z_zz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(32));

    t_z_zz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(33));

    t_z_zz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(34));

    t_z_zz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(35));

    t_z_yz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(36));

    t_z_yz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(37));

    t_z_yz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(38));

    t_z_yz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(39));

    t_z_yz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(40));

    t_z_yz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(41));

    t_z_yz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(42));

    t_z_yz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(43));

    t_z_yz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(44));

    t_z_yz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(45));

    t_z_yz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(46));

    t_z_yz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(47));

    t_z_yz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(48));

    t_z_yz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(49));

    t_z_yz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(50));

    t_z_yz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(51));

    t_z_yz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(52));

    t_z_yz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(53));

    t_z_yz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(54));

    t_z_yz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(55));

    t_z_yz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(56));

    t_z_yz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(57));

    t_z_yz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(58));

    t_z_yz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(59));

    t_z_yz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(60));

    t_z_yz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(61));

    t_z_yz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(62));

    t_z_yz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(63));

    t_z_yz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(64));

    t_z_yz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(65));

    t_z_yz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(66));

    t_z_yz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(67));

    t_z_yz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(68));

    t_z_yz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(69));

    t_z_yz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(70));

    t_z_yz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(71));

    t_z_yy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(72));

    t_z_yy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(73));

    t_z_yy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(74));

    t_z_yy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(75));

    t_z_yy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(76));

    t_z_yy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(77));

    t_z_yy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(78));

    t_z_yy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(79));

    t_z_yy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(80));

    t_z_yy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(81));

    t_z_yy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(82));

    t_z_yy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(83));

    t_z_yy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(84));

    t_z_yy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(85));

    t_z_yy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(86));

    t_z_yy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(87));

    t_z_yy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(88));

    t_z_yy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(89));

    t_z_yy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(90));

    t_z_yy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(91));

    t_z_yy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(92));

    t_z_yy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(93));

    t_z_yy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(94));

    t_z_yy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(95));

    t_z_yy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(96));

    t_z_yy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(97));

    t_z_yy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(98));

    t_z_yy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(99));

    t_z_yy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(100));

    t_z_yy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(101));

    t_z_yy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(102));

    t_z_yy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(103));

    t_z_yy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(104));

    t_z_yy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(105));

    t_z_yy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(106));

    t_z_yy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(107));

    t_z_xz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(108));

    t_z_xz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(109));

    t_z_xz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(110));

    t_z_xz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(111));

    t_z_xz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(112));

    t_z_xz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(113));

    t_z_xz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(114));

    t_z_xz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(115));

    t_z_xz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(116));

    t_z_xz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(117));

    t_z_xz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(118));

    t_z_xz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(119));

    t_z_xz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(120));

    t_z_xz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(121));

    t_z_xz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(122));

    t_z_xz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(123));

    t_z_xz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(124));

    t_z_xz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(125));

    t_z_xz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(126));

    t_z_xz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(127));

    t_z_xz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(128));

    t_z_xz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(129));

    t_z_xz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(130));

    t_z_xz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(131));

    t_z_xz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(132));

    t_z_xz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(133));

    t_z_xz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(134));

    t_z_xz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(135));

    t_z_xz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(136));

    t_z_xz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(137));

    t_z_xz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(138));

    t_z_xz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(139));

    t_z_xz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(140));

    t_z_xz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(141));

    t_z_xz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(142));

    t_z_xz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(143));

    t_z_xy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(144));

    t_z_xy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(145));

    t_z_xy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(146));

    t_z_xy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(147));

    t_z_xy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(148));

    t_z_xy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(149));

    t_z_xy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(150));

    t_z_xy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(151));

    t_z_xy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(152));

    t_z_xy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(153));

    t_z_xy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(154));

    t_z_xy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(155));

    t_z_xy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(156));

    t_z_xy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(157));

    t_z_xy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(158));

    t_z_xy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(159));

    t_z_xy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(160));

    t_z_xy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(161));

    t_z_xy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(162));

    t_z_xy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(163));

    t_z_xy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(164));

    t_z_xy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(165));

    t_z_xy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(166));

    t_z_xy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(167));

    t_z_xy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(168));

    t_z_xy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(169));

    t_z_xy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(170));

    t_z_xy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(171));

    t_z_xy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(172));

    t_z_xy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(173));

    t_z_xy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(174));

    t_z_xy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(175));

    t_z_xy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(176));

    t_z_xy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(177));

    t_z_xy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(178));

    t_z_xy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(179));

    t_z_xx_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(180));

    t_z_xx_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(181));

    t_z_xx_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(182));

    t_z_xx_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(183));

    t_z_xx_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(184));

    t_z_xx_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(185));

    t_z_xx_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(186));

    t_z_xx_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(187));

    t_z_xx_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(188));

    t_z_xx_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(189));

    t_z_xx_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(190));

    t_z_xx_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(191));

    t_z_xx_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(192));

    t_z_xx_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(193));

    t_z_xx_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(194));

    t_z_xx_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(195));

    t_z_xx_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(196));

    t_z_xx_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(197));

    t_z_xx_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(198));

    t_z_xx_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(199));

    t_z_xx_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(200));

    t_z_xx_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(201));

    t_z_xx_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(202));

    t_z_xx_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(203));

    t_z_xx_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(204));

    t_z_xx_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(205));

    t_z_xx_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(206));

    t_z_xx_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(207));

    t_z_xx_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(208));

    t_z_xx_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(209));

    t_z_xx_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(210));

    t_z_xx_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(211));

    t_z_xx_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(212));

    t_z_xx_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(213));

    t_z_xx_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(214));

    t_z_xx_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(215));

    t_y_zz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(216));

    t_y_zz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(217));

    t_y_zz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(218));

    t_y_zz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(219));

    t_y_zz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(220));

    t_y_zz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(221));

    t_y_zz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(222));

    t_y_zz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(223));

    t_y_zz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(224));

    t_y_zz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(225));

    t_y_zz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(226));

    t_y_zz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(227));

    t_y_zz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(228));

    t_y_zz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(229));

    t_y_zz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(230));

    t_y_zz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(231));

    t_y_zz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(232));

    t_y_zz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(233));

    t_y_zz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(234));

    t_y_zz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(235));

    t_y_zz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(236));

    t_y_zz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(237));

    t_y_zz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(238));

    t_y_zz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(239));

    t_y_zz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(240));

    t_y_zz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(241));

    t_y_zz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(242));

    t_y_zz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(243));

    t_y_zz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(244));

    t_y_zz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(245));

    t_y_zz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(246));

    t_y_zz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(247));

    t_y_zz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(248));

    t_y_zz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(249));

    t_y_zz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(250));

    t_y_zz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(251));

    t_y_yz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(252));

    t_y_yz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(253));

    t_y_yz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(254));

    t_y_yz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(255));

    t_y_yz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(256));

    t_y_yz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(257));

    t_y_yz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(258));

    t_y_yz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(259));

    t_y_yz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(260));

    t_y_yz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(261));

    t_y_yz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(262));

    t_y_yz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(263));

    t_y_yz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(264));

    t_y_yz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(265));

    t_y_yz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(266));

    t_y_yz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(267));

    t_y_yz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(268));

    t_y_yz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(269));

    t_y_yz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(270));

    t_y_yz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(271));

    t_y_yz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(272));

    t_y_yz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(273));

    t_y_yz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(274));

    t_y_yz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(275));

    t_y_yz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(276));

    t_y_yz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(277));

    t_y_yz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(278));

    t_y_yz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(279));

    t_y_yz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(280));

    t_y_yz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(281));

    t_y_yz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(282));

    t_y_yz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(283));

    t_y_yz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(284));

    t_y_yz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(285));

    t_y_yz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(286));

    t_y_yz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(287));

    t_y_yy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(288));

    t_y_yy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(289));

    t_y_yy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(290));

    t_y_yy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(291));

    t_y_yy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(292));

    t_y_yy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(293));

    t_y_yy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(294));

    t_y_yy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(295));

    t_y_yy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(296));

    t_y_yy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(297));

    t_y_yy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(298));

    t_y_yy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(299));

    t_y_yy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(300));

    t_y_yy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(301));

    t_y_yy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(302));

    t_y_yy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(303));

    t_y_yy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(304));

    t_y_yy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(305));

    t_y_yy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(306));

    t_y_yy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(307));

    t_y_yy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(308));

    t_y_yy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(309));

    t_y_yy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(310));

    t_y_yy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(311));

    t_y_yy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(312));

    t_y_yy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(313));

    t_y_yy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(314));

    t_y_yy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(315));

    t_y_yy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(316));

    t_y_yy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(317));

    t_y_yy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(318));

    t_y_yy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(319));

    t_y_yy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(320));

    t_y_yy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(321));

    t_y_yy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(322));

    t_y_yy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(323));

    t_y_xz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(324));

    t_y_xz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(325));

    t_y_xz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(326));

    t_y_xz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(327));

    t_y_xz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(328));

    t_y_xz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(329));

    t_y_xz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(330));

    t_y_xz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(331));

    t_y_xz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(332));

    t_y_xz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(333));

    t_y_xz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(334));

    t_y_xz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(335));

    t_y_xz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(336));

    t_y_xz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(337));

    t_y_xz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(338));

    t_y_xz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(339));

    t_y_xz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(340));

    t_y_xz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(341));

    t_y_xz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(342));

    t_y_xz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(343));

    t_y_xz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(344));

    t_y_xz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(345));

    t_y_xz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(346));

    t_y_xz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(347));

    t_y_xz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(348));

    t_y_xz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(349));

    t_y_xz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(350));

    t_y_xz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(351));

    t_y_xz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(352));

    t_y_xz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(353));

    t_y_xz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(354));

    t_y_xz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(355));

    t_y_xz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(356));

    t_y_xz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(357));

    t_y_xz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(358));

    t_y_xz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(359));

    t_y_xy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(360));

    t_y_xy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(361));

    t_y_xy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(362));

    t_y_xy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(363));

    t_y_xy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(364));

    t_y_xy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(365));

    t_y_xy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(366));

    t_y_xy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(367));

    t_y_xy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(368));

    t_y_xy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(369));

    t_y_xy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(370));

    t_y_xy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(371));

    t_y_xy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(372));

    t_y_xy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(373));

    t_y_xy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(374));

    t_y_xy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(375));

    t_y_xy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(376));

    t_y_xy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(377));

    t_y_xy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(378));

    t_y_xy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(379));

    t_y_xy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(380));

    t_y_xy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(381));

    t_y_xy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(382));

    t_y_xy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(383));

    t_y_xy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(384));

    t_y_xy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(385));

    t_y_xy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(386));

    t_y_xy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(387));

    t_y_xy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(388));

    t_y_xy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(389));

    t_y_xy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(390));

    t_y_xy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(391));

    t_y_xy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(392));

    t_y_xy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(393));

    t_y_xy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(394));

    t_y_xy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(395));

    t_y_xx_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(396));

    t_y_xx_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(397));

    t_y_xx_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(398));

    t_y_xx_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(399));

    t_y_xx_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(400));

    t_y_xx_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(401));

    t_y_xx_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(402));

    t_y_xx_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(403));

    t_y_xx_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(404));

    t_y_xx_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(405));

    t_y_xx_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(406));

    t_y_xx_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(407));

    t_y_xx_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(408));

    t_y_xx_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(409));

    t_y_xx_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(410));

    t_y_xx_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(411));

    t_y_xx_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(412));

    t_y_xx_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(413));

    t_y_xx_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(414));

    t_y_xx_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(415));

    t_y_xx_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(416));

    t_y_xx_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(417));

    t_y_xx_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(418));

    t_y_xx_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(419));

    t_y_xx_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(420));

    t_y_xx_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(421));

    t_y_xx_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(422));

    t_y_xx_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(423));

    t_y_xx_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(424));

    t_y_xx_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(425));

    t_y_xx_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(426));

    t_y_xx_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(427));

    t_y_xx_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(428));

    t_y_xx_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(429));

    t_y_xx_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(430));

    t_y_xx_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(431));

    t_x_zz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(432));

    t_x_zz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(433));

    t_x_zz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(434));

    t_x_zz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(435));

    t_x_zz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(436));

    t_x_zz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(437));

    t_x_zz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(438));

    t_x_zz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(439));

    t_x_zz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(440));

    t_x_zz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(441));

    t_x_zz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(442));

    t_x_zz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(443));

    t_x_zz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(444));

    t_x_zz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(445));

    t_x_zz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(446));

    t_x_zz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(447));

    t_x_zz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(448));

    t_x_zz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(449));

    t_x_zz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(450));

    t_x_zz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(451));

    t_x_zz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(452));

    t_x_zz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(453));

    t_x_zz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(454));

    t_x_zz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(455));

    t_x_zz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(456));

    t_x_zz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(457));

    t_x_zz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(458));

    t_x_zz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(459));

    t_x_zz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(460));

    t_x_zz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(461));

    t_x_zz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(462));

    t_x_zz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(463));

    t_x_zz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(464));

    t_x_zz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(465));

    t_x_zz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(466));

    t_x_zz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(467));

    t_x_yz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(468));

    t_x_yz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(469));

    t_x_yz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(470));

    t_x_yz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(471));

    t_x_yz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(472));

    t_x_yz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(473));

    t_x_yz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(474));

    t_x_yz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(475));

    t_x_yz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(476));

    t_x_yz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(477));

    t_x_yz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(478));

    t_x_yz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(479));

    t_x_yz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(480));

    t_x_yz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(481));

    t_x_yz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(482));

    t_x_yz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(483));

    t_x_yz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(484));

    t_x_yz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(485));

    t_x_yz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(486));

    t_x_yz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(487));

    t_x_yz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(488));

    t_x_yz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(489));

    t_x_yz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(490));

    t_x_yz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(491));

    t_x_yz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(492));

    t_x_yz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(493));

    t_x_yz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(494));

    t_x_yz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(495));

    t_x_yz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(496));

    t_x_yz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(497));

    t_x_yz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(498));

    t_x_yz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(499));

    t_x_yz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(500));

    t_x_yz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(501));

    t_x_yz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(502));

    t_x_yz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(503));

    t_x_yy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(504));

    t_x_yy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(505));

    t_x_yy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(506));

    t_x_yy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(507));

    t_x_yy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(508));

    t_x_yy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(509));

    t_x_yy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(510));

    t_x_yy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(511));

    t_x_yy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(512));

    t_x_yy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(513));

    t_x_yy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(514));

    t_x_yy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(515));

    t_x_yy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(516));

    t_x_yy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(517));

    t_x_yy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(518));

    t_x_yy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(519));

    t_x_yy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(520));

    t_x_yy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(521));

    t_x_yy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(522));

    t_x_yy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(523));

    t_x_yy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(524));

    t_x_yy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(525));

    t_x_yy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(526));

    t_x_yy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(527));

    t_x_yy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(528));

    t_x_yy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(529));

    t_x_yy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(530));

    t_x_yy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(531));

    t_x_yy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(532));

    t_x_yy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(533));

    t_x_yy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(534));

    t_x_yy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(535));

    t_x_yy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(536));

    t_x_yy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(537));

    t_x_yy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(538));

    t_x_yy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(539));

    t_x_xz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(540));

    t_x_xz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(541));

    t_x_xz_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(542));

    t_x_xz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(543));

    t_x_xz_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(544));

    t_x_xz_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(545));

    t_x_xz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(546));

    t_x_xz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(547));

    t_x_xz_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(548));

    t_x_xz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(549));

    t_x_xz_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(550));

    t_x_xz_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(551));

    t_x_xz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(552));

    t_x_xz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(553));

    t_x_xz_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(554));

    t_x_xz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(555));

    t_x_xz_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(556));

    t_x_xz_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(557));

    t_x_xz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(558));

    t_x_xz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(559));

    t_x_xz_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(560));

    t_x_xz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(561));

    t_x_xz_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(562));

    t_x_xz_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(563));

    t_x_xz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(564));

    t_x_xz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(565));

    t_x_xz_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(566));

    t_x_xz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(567));

    t_x_xz_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(568));

    t_x_xz_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(569));

    t_x_xz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(570));

    t_x_xz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(571));

    t_x_xz_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(572));

    t_x_xz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(573));

    t_x_xz_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(574));

    t_x_xz_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(575));

    t_x_xy_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(576));

    t_x_xy_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(577));

    t_x_xy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(578));

    t_x_xy_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(579));

    t_x_xy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(580));

    t_x_xy_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(581));

    t_x_xy_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(582));

    t_x_xy_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(583));

    t_x_xy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(584));

    t_x_xy_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(585));

    t_x_xy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(586));

    t_x_xy_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(587));

    t_x_xy_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(588));

    t_x_xy_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(589));

    t_x_xy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(590));

    t_x_xy_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(591));

    t_x_xy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(592));

    t_x_xy_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(593));

    t_x_xy_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(594));

    t_x_xy_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(595));

    t_x_xy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(596));

    t_x_xy_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(597));

    t_x_xy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(598));

    t_x_xy_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(599));

    t_x_xy_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(600));

    t_x_xy_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(601));

    t_x_xy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(602));

    t_x_xy_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(603));

    t_x_xy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(604));

    t_x_xy_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(605));

    t_x_xy_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(606));

    t_x_xy_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(607));

    t_x_xy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(608));

    t_x_xy_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(609));

    t_x_xy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(610));

    t_x_xy_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(611));

    t_x_xx_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(612));

    t_x_xx_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(613));

    t_x_xx_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(614));

    t_x_xx_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(615));

    t_x_xx_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(616));

    t_x_xx_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(617));

    t_x_xx_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(618));

    t_x_xx_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(619));

    t_x_xx_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(620));

    t_x_xx_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(621));

    t_x_xx_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(622));

    t_x_xx_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(623));

    t_x_xx_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(624));

    t_x_xx_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(625));

    t_x_xx_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(626));

    t_x_xx_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(627));

    t_x_xx_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(628));

    t_x_xx_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(629));

    t_x_xx_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(630));

    t_x_xx_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(631));

    t_x_xx_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(632));

    t_x_xx_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(633));

    t_x_xx_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(634));

    t_x_xx_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(635));

    t_x_xx_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(636));

    t_x_xx_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(637));

    t_x_xx_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(638));

    t_x_xx_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(639));

    t_x_xx_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(640));

    t_x_xx_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(641));

    t_x_xx_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(642));

    t_x_xx_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(643));

    t_x_xx_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(644));

    t_x_xx_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(645));

    t_x_xx_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(646));

    t_x_xx_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(647));

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

    #pragma omp simd align(rab_z, t_0_zz_xx_xx, t_0_zz_xx_xy, t_0_zz_xx_xz, t_0_zz_xx_yy,\
                           t_0_zz_xx_yz, t_0_zz_xx_zz, t_0_zz_xy_xx, t_0_zz_xy_xy, t_0_zz_xy_xz,\
                           t_0_zz_xy_yy, t_0_zz_xy_yz, t_0_zz_xy_zz, t_0_zz_xz_xx, t_0_zz_xz_xy,\
                           t_0_zz_xz_xz, t_0_zz_xz_yy, t_0_zz_xz_yz, t_0_zz_xz_zz, t_0_zz_yy_xx,\
                           t_0_zz_yy_xy, t_0_zz_yy_xz, t_0_zz_yy_yy, t_0_zz_yy_yz, t_0_zz_yy_zz,\
                           t_0_zz_yz_xx, t_0_zz_yz_xy, t_0_zz_yz_xz, t_0_zz_yz_yy, t_0_zz_yz_yz,\
                           t_0_zz_yz_zz, t_0_zz_zz_xx, t_0_zz_zz_xy, t_0_zz_zz_xz, t_0_zz_zz_yy,\
                           t_0_zz_zz_yz, t_0_zz_zz_zz, t_0_zzz_xx_xx, t_0_zzz_xx_xy,\
                           t_0_zzz_xx_xz, t_0_zzz_xx_yy, t_0_zzz_xx_yz, t_0_zzz_xx_zz,\
                           t_0_zzz_xy_xx, t_0_zzz_xy_xy, t_0_zzz_xy_xz, t_0_zzz_xy_yy,\
                           t_0_zzz_xy_yz, t_0_zzz_xy_zz, t_0_zzz_xz_xx, t_0_zzz_xz_xy,\
                           t_0_zzz_xz_xz, t_0_zzz_xz_yy, t_0_zzz_xz_yz, t_0_zzz_xz_zz,\
                           t_0_zzz_yy_xx, t_0_zzz_yy_xy, t_0_zzz_yy_xz, t_0_zzz_yy_yy,\
                           t_0_zzz_yy_yz, t_0_zzz_yy_zz, t_0_zzz_yz_xx, t_0_zzz_yz_xy,\
                           t_0_zzz_yz_xz, t_0_zzz_yz_yy, t_0_zzz_yz_yz, t_0_zzz_yz_zz,\
                           t_0_zzz_zz_xx, t_0_zzz_zz_xy, t_0_zzz_zz_xz, t_0_zzz_zz_yy,\
                           t_0_zzz_zz_yz, t_0_zzz_zz_zz, t_z_zz_xx_xx, t_z_zz_xx_xy,\
                           t_z_zz_xx_xz, t_z_zz_xx_yy, t_z_zz_xx_yz, t_z_zz_xx_zz, t_z_zz_xy_xx,\
                           t_z_zz_xy_xy, t_z_zz_xy_xz, t_z_zz_xy_yy, t_z_zz_xy_yz, t_z_zz_xy_zz,\
                           t_z_zz_xz_xx, t_z_zz_xz_xy, t_z_zz_xz_xz, t_z_zz_xz_yy, t_z_zz_xz_yz,\
                           t_z_zz_xz_zz, t_z_zz_yy_xx, t_z_zz_yy_xy, t_z_zz_yy_xz, t_z_zz_yy_yy,\
                           t_z_zz_yy_yz, t_z_zz_yy_zz, t_z_zz_yz_xx, t_z_zz_yz_xy, t_z_zz_yz_xz,\
                           t_z_zz_yz_yy, t_z_zz_yz_yz, t_z_zz_yz_zz, t_z_zz_zz_xx, t_z_zz_zz_xy,\
                           t_z_zz_zz_xz, t_z_zz_zz_yy, t_z_zz_zz_yz, t_z_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_zz_zz[i] = t_0_zzz_zz_zz[i] - rab_z[i] * t_0_zz_zz_zz[i];

        t_z_zz_zz_yz[i] = t_0_zzz_zz_yz[i] - rab_z[i] * t_0_zz_zz_yz[i];

        t_z_zz_zz_yy[i] = t_0_zzz_zz_yy[i] - rab_z[i] * t_0_zz_zz_yy[i];

        t_z_zz_zz_xz[i] = t_0_zzz_zz_xz[i] - rab_z[i] * t_0_zz_zz_xz[i];

        t_z_zz_zz_xy[i] = t_0_zzz_zz_xy[i] - rab_z[i] * t_0_zz_zz_xy[i];

        t_z_zz_zz_xx[i] = t_0_zzz_zz_xx[i] - rab_z[i] * t_0_zz_zz_xx[i];

        t_z_zz_yz_zz[i] = t_0_zzz_yz_zz[i] - rab_z[i] * t_0_zz_yz_zz[i];

        t_z_zz_yz_yz[i] = t_0_zzz_yz_yz[i] - rab_z[i] * t_0_zz_yz_yz[i];

        t_z_zz_yz_yy[i] = t_0_zzz_yz_yy[i] - rab_z[i] * t_0_zz_yz_yy[i];

        t_z_zz_yz_xz[i] = t_0_zzz_yz_xz[i] - rab_z[i] * t_0_zz_yz_xz[i];

        t_z_zz_yz_xy[i] = t_0_zzz_yz_xy[i] - rab_z[i] * t_0_zz_yz_xy[i];

        t_z_zz_yz_xx[i] = t_0_zzz_yz_xx[i] - rab_z[i] * t_0_zz_yz_xx[i];

        t_z_zz_yy_zz[i] = t_0_zzz_yy_zz[i] - rab_z[i] * t_0_zz_yy_zz[i];

        t_z_zz_yy_yz[i] = t_0_zzz_yy_yz[i] - rab_z[i] * t_0_zz_yy_yz[i];

        t_z_zz_yy_yy[i] = t_0_zzz_yy_yy[i] - rab_z[i] * t_0_zz_yy_yy[i];

        t_z_zz_yy_xz[i] = t_0_zzz_yy_xz[i] - rab_z[i] * t_0_zz_yy_xz[i];

        t_z_zz_yy_xy[i] = t_0_zzz_yy_xy[i] - rab_z[i] * t_0_zz_yy_xy[i];

        t_z_zz_yy_xx[i] = t_0_zzz_yy_xx[i] - rab_z[i] * t_0_zz_yy_xx[i];

        t_z_zz_xz_zz[i] = t_0_zzz_xz_zz[i] - rab_z[i] * t_0_zz_xz_zz[i];

        t_z_zz_xz_yz[i] = t_0_zzz_xz_yz[i] - rab_z[i] * t_0_zz_xz_yz[i];

        t_z_zz_xz_yy[i] = t_0_zzz_xz_yy[i] - rab_z[i] * t_0_zz_xz_yy[i];

        t_z_zz_xz_xz[i] = t_0_zzz_xz_xz[i] - rab_z[i] * t_0_zz_xz_xz[i];

        t_z_zz_xz_xy[i] = t_0_zzz_xz_xy[i] - rab_z[i] * t_0_zz_xz_xy[i];

        t_z_zz_xz_xx[i] = t_0_zzz_xz_xx[i] - rab_z[i] * t_0_zz_xz_xx[i];

        t_z_zz_xy_zz[i] = t_0_zzz_xy_zz[i] - rab_z[i] * t_0_zz_xy_zz[i];

        t_z_zz_xy_yz[i] = t_0_zzz_xy_yz[i] - rab_z[i] * t_0_zz_xy_yz[i];

        t_z_zz_xy_yy[i] = t_0_zzz_xy_yy[i] - rab_z[i] * t_0_zz_xy_yy[i];

        t_z_zz_xy_xz[i] = t_0_zzz_xy_xz[i] - rab_z[i] * t_0_zz_xy_xz[i];

        t_z_zz_xy_xy[i] = t_0_zzz_xy_xy[i] - rab_z[i] * t_0_zz_xy_xy[i];

        t_z_zz_xy_xx[i] = t_0_zzz_xy_xx[i] - rab_z[i] * t_0_zz_xy_xx[i];

        t_z_zz_xx_zz[i] = t_0_zzz_xx_zz[i] - rab_z[i] * t_0_zz_xx_zz[i];

        t_z_zz_xx_yz[i] = t_0_zzz_xx_yz[i] - rab_z[i] * t_0_zz_xx_yz[i];

        t_z_zz_xx_yy[i] = t_0_zzz_xx_yy[i] - rab_z[i] * t_0_zz_xx_yy[i];

        t_z_zz_xx_xz[i] = t_0_zzz_xx_xz[i] - rab_z[i] * t_0_zz_xx_xz[i];

        t_z_zz_xx_xy[i] = t_0_zzz_xx_xy[i] - rab_z[i] * t_0_zz_xx_xy[i];

        t_z_zz_xx_xx[i] = t_0_zzz_xx_xx[i] - rab_z[i] * t_0_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_yz_xx_xx, t_0_yz_xx_xy, t_0_yz_xx_xz, t_0_yz_xx_yy,\
                           t_0_yz_xx_yz, t_0_yz_xx_zz, t_0_yz_xy_xx, t_0_yz_xy_xy, t_0_yz_xy_xz,\
                           t_0_yz_xy_yy, t_0_yz_xy_yz, t_0_yz_xy_zz, t_0_yz_xz_xx, t_0_yz_xz_xy,\
                           t_0_yz_xz_xz, t_0_yz_xz_yy, t_0_yz_xz_yz, t_0_yz_xz_zz, t_0_yz_yy_xx,\
                           t_0_yz_yy_xy, t_0_yz_yy_xz, t_0_yz_yy_yy, t_0_yz_yy_yz, t_0_yz_yy_zz,\
                           t_0_yz_yz_xx, t_0_yz_yz_xy, t_0_yz_yz_xz, t_0_yz_yz_yy, t_0_yz_yz_yz,\
                           t_0_yz_yz_zz, t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz, t_0_yz_zz_yy,\
                           t_0_yz_zz_yz, t_0_yz_zz_zz, t_0_yzz_xx_xx, t_0_yzz_xx_xy,\
                           t_0_yzz_xx_xz, t_0_yzz_xx_yy, t_0_yzz_xx_yz, t_0_yzz_xx_zz,\
                           t_0_yzz_xy_xx, t_0_yzz_xy_xy, t_0_yzz_xy_xz, t_0_yzz_xy_yy,\
                           t_0_yzz_xy_yz, t_0_yzz_xy_zz, t_0_yzz_xz_xx, t_0_yzz_xz_xy,\
                           t_0_yzz_xz_xz, t_0_yzz_xz_yy, t_0_yzz_xz_yz, t_0_yzz_xz_zz,\
                           t_0_yzz_yy_xx, t_0_yzz_yy_xy, t_0_yzz_yy_xz, t_0_yzz_yy_yy,\
                           t_0_yzz_yy_yz, t_0_yzz_yy_zz, t_0_yzz_yz_xx, t_0_yzz_yz_xy,\
                           t_0_yzz_yz_xz, t_0_yzz_yz_yy, t_0_yzz_yz_yz, t_0_yzz_yz_zz,\
                           t_0_yzz_zz_xx, t_0_yzz_zz_xy, t_0_yzz_zz_xz, t_0_yzz_zz_yy,\
                           t_0_yzz_zz_yz, t_0_yzz_zz_zz, t_z_yz_xx_xx, t_z_yz_xx_xy,\
                           t_z_yz_xx_xz, t_z_yz_xx_yy, t_z_yz_xx_yz, t_z_yz_xx_zz, t_z_yz_xy_xx,\
                           t_z_yz_xy_xy, t_z_yz_xy_xz, t_z_yz_xy_yy, t_z_yz_xy_yz, t_z_yz_xy_zz,\
                           t_z_yz_xz_xx, t_z_yz_xz_xy, t_z_yz_xz_xz, t_z_yz_xz_yy, t_z_yz_xz_yz,\
                           t_z_yz_xz_zz, t_z_yz_yy_xx, t_z_yz_yy_xy, t_z_yz_yy_xz, t_z_yz_yy_yy,\
                           t_z_yz_yy_yz, t_z_yz_yy_zz, t_z_yz_yz_xx, t_z_yz_yz_xy, t_z_yz_yz_xz,\
                           t_z_yz_yz_yy, t_z_yz_yz_yz, t_z_yz_yz_zz, t_z_yz_zz_xx, t_z_yz_zz_xy,\
                           t_z_yz_zz_xz, t_z_yz_zz_yy, t_z_yz_zz_yz, t_z_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yz_zz_zz[i] = t_0_yzz_zz_zz[i] - rab_z[i] * t_0_yz_zz_zz[i];

        t_z_yz_zz_yz[i] = t_0_yzz_zz_yz[i] - rab_z[i] * t_0_yz_zz_yz[i];

        t_z_yz_zz_yy[i] = t_0_yzz_zz_yy[i] - rab_z[i] * t_0_yz_zz_yy[i];

        t_z_yz_zz_xz[i] = t_0_yzz_zz_xz[i] - rab_z[i] * t_0_yz_zz_xz[i];

        t_z_yz_zz_xy[i] = t_0_yzz_zz_xy[i] - rab_z[i] * t_0_yz_zz_xy[i];

        t_z_yz_zz_xx[i] = t_0_yzz_zz_xx[i] - rab_z[i] * t_0_yz_zz_xx[i];

        t_z_yz_yz_zz[i] = t_0_yzz_yz_zz[i] - rab_z[i] * t_0_yz_yz_zz[i];

        t_z_yz_yz_yz[i] = t_0_yzz_yz_yz[i] - rab_z[i] * t_0_yz_yz_yz[i];

        t_z_yz_yz_yy[i] = t_0_yzz_yz_yy[i] - rab_z[i] * t_0_yz_yz_yy[i];

        t_z_yz_yz_xz[i] = t_0_yzz_yz_xz[i] - rab_z[i] * t_0_yz_yz_xz[i];

        t_z_yz_yz_xy[i] = t_0_yzz_yz_xy[i] - rab_z[i] * t_0_yz_yz_xy[i];

        t_z_yz_yz_xx[i] = t_0_yzz_yz_xx[i] - rab_z[i] * t_0_yz_yz_xx[i];

        t_z_yz_yy_zz[i] = t_0_yzz_yy_zz[i] - rab_z[i] * t_0_yz_yy_zz[i];

        t_z_yz_yy_yz[i] = t_0_yzz_yy_yz[i] - rab_z[i] * t_0_yz_yy_yz[i];

        t_z_yz_yy_yy[i] = t_0_yzz_yy_yy[i] - rab_z[i] * t_0_yz_yy_yy[i];

        t_z_yz_yy_xz[i] = t_0_yzz_yy_xz[i] - rab_z[i] * t_0_yz_yy_xz[i];

        t_z_yz_yy_xy[i] = t_0_yzz_yy_xy[i] - rab_z[i] * t_0_yz_yy_xy[i];

        t_z_yz_yy_xx[i] = t_0_yzz_yy_xx[i] - rab_z[i] * t_0_yz_yy_xx[i];

        t_z_yz_xz_zz[i] = t_0_yzz_xz_zz[i] - rab_z[i] * t_0_yz_xz_zz[i];

        t_z_yz_xz_yz[i] = t_0_yzz_xz_yz[i] - rab_z[i] * t_0_yz_xz_yz[i];

        t_z_yz_xz_yy[i] = t_0_yzz_xz_yy[i] - rab_z[i] * t_0_yz_xz_yy[i];

        t_z_yz_xz_xz[i] = t_0_yzz_xz_xz[i] - rab_z[i] * t_0_yz_xz_xz[i];

        t_z_yz_xz_xy[i] = t_0_yzz_xz_xy[i] - rab_z[i] * t_0_yz_xz_xy[i];

        t_z_yz_xz_xx[i] = t_0_yzz_xz_xx[i] - rab_z[i] * t_0_yz_xz_xx[i];

        t_z_yz_xy_zz[i] = t_0_yzz_xy_zz[i] - rab_z[i] * t_0_yz_xy_zz[i];

        t_z_yz_xy_yz[i] = t_0_yzz_xy_yz[i] - rab_z[i] * t_0_yz_xy_yz[i];

        t_z_yz_xy_yy[i] = t_0_yzz_xy_yy[i] - rab_z[i] * t_0_yz_xy_yy[i];

        t_z_yz_xy_xz[i] = t_0_yzz_xy_xz[i] - rab_z[i] * t_0_yz_xy_xz[i];

        t_z_yz_xy_xy[i] = t_0_yzz_xy_xy[i] - rab_z[i] * t_0_yz_xy_xy[i];

        t_z_yz_xy_xx[i] = t_0_yzz_xy_xx[i] - rab_z[i] * t_0_yz_xy_xx[i];

        t_z_yz_xx_zz[i] = t_0_yzz_xx_zz[i] - rab_z[i] * t_0_yz_xx_zz[i];

        t_z_yz_xx_yz[i] = t_0_yzz_xx_yz[i] - rab_z[i] * t_0_yz_xx_yz[i];

        t_z_yz_xx_yy[i] = t_0_yzz_xx_yy[i] - rab_z[i] * t_0_yz_xx_yy[i];

        t_z_yz_xx_xz[i] = t_0_yzz_xx_xz[i] - rab_z[i] * t_0_yz_xx_xz[i];

        t_z_yz_xx_xy[i] = t_0_yzz_xx_xy[i] - rab_z[i] * t_0_yz_xx_xy[i];

        t_z_yz_xx_xx[i] = t_0_yzz_xx_xx[i] - rab_z[i] * t_0_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_yy_xx_xx, t_0_yy_xx_xy, t_0_yy_xx_xz, t_0_yy_xx_yy,\
                           t_0_yy_xx_yz, t_0_yy_xx_zz, t_0_yy_xy_xx, t_0_yy_xy_xy, t_0_yy_xy_xz,\
                           t_0_yy_xy_yy, t_0_yy_xy_yz, t_0_yy_xy_zz, t_0_yy_xz_xx, t_0_yy_xz_xy,\
                           t_0_yy_xz_xz, t_0_yy_xz_yy, t_0_yy_xz_yz, t_0_yy_xz_zz, t_0_yy_yy_xx,\
                           t_0_yy_yy_xy, t_0_yy_yy_xz, t_0_yy_yy_yy, t_0_yy_yy_yz, t_0_yy_yy_zz,\
                           t_0_yy_yz_xx, t_0_yy_yz_xy, t_0_yy_yz_xz, t_0_yy_yz_yy, t_0_yy_yz_yz,\
                           t_0_yy_yz_zz, t_0_yy_zz_xx, t_0_yy_zz_xy, t_0_yy_zz_xz, t_0_yy_zz_yy,\
                           t_0_yy_zz_yz, t_0_yy_zz_zz, t_0_yyz_xx_xx, t_0_yyz_xx_xy,\
                           t_0_yyz_xx_xz, t_0_yyz_xx_yy, t_0_yyz_xx_yz, t_0_yyz_xx_zz,\
                           t_0_yyz_xy_xx, t_0_yyz_xy_xy, t_0_yyz_xy_xz, t_0_yyz_xy_yy,\
                           t_0_yyz_xy_yz, t_0_yyz_xy_zz, t_0_yyz_xz_xx, t_0_yyz_xz_xy,\
                           t_0_yyz_xz_xz, t_0_yyz_xz_yy, t_0_yyz_xz_yz, t_0_yyz_xz_zz,\
                           t_0_yyz_yy_xx, t_0_yyz_yy_xy, t_0_yyz_yy_xz, t_0_yyz_yy_yy,\
                           t_0_yyz_yy_yz, t_0_yyz_yy_zz, t_0_yyz_yz_xx, t_0_yyz_yz_xy,\
                           t_0_yyz_yz_xz, t_0_yyz_yz_yy, t_0_yyz_yz_yz, t_0_yyz_yz_zz,\
                           t_0_yyz_zz_xx, t_0_yyz_zz_xy, t_0_yyz_zz_xz, t_0_yyz_zz_yy,\
                           t_0_yyz_zz_yz, t_0_yyz_zz_zz, t_z_yy_xx_xx, t_z_yy_xx_xy,\
                           t_z_yy_xx_xz, t_z_yy_xx_yy, t_z_yy_xx_yz, t_z_yy_xx_zz, t_z_yy_xy_xx,\
                           t_z_yy_xy_xy, t_z_yy_xy_xz, t_z_yy_xy_yy, t_z_yy_xy_yz, t_z_yy_xy_zz,\
                           t_z_yy_xz_xx, t_z_yy_xz_xy, t_z_yy_xz_xz, t_z_yy_xz_yy, t_z_yy_xz_yz,\
                           t_z_yy_xz_zz, t_z_yy_yy_xx, t_z_yy_yy_xy, t_z_yy_yy_xz, t_z_yy_yy_yy,\
                           t_z_yy_yy_yz, t_z_yy_yy_zz, t_z_yy_yz_xx, t_z_yy_yz_xy, t_z_yy_yz_xz,\
                           t_z_yy_yz_yy, t_z_yy_yz_yz, t_z_yy_yz_zz, t_z_yy_zz_xx, t_z_yy_zz_xy,\
                           t_z_yy_zz_xz, t_z_yy_zz_yy, t_z_yy_zz_yz, t_z_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yy_zz_zz[i] = t_0_yyz_zz_zz[i] - rab_z[i] * t_0_yy_zz_zz[i];

        t_z_yy_zz_yz[i] = t_0_yyz_zz_yz[i] - rab_z[i] * t_0_yy_zz_yz[i];

        t_z_yy_zz_yy[i] = t_0_yyz_zz_yy[i] - rab_z[i] * t_0_yy_zz_yy[i];

        t_z_yy_zz_xz[i] = t_0_yyz_zz_xz[i] - rab_z[i] * t_0_yy_zz_xz[i];

        t_z_yy_zz_xy[i] = t_0_yyz_zz_xy[i] - rab_z[i] * t_0_yy_zz_xy[i];

        t_z_yy_zz_xx[i] = t_0_yyz_zz_xx[i] - rab_z[i] * t_0_yy_zz_xx[i];

        t_z_yy_yz_zz[i] = t_0_yyz_yz_zz[i] - rab_z[i] * t_0_yy_yz_zz[i];

        t_z_yy_yz_yz[i] = t_0_yyz_yz_yz[i] - rab_z[i] * t_0_yy_yz_yz[i];

        t_z_yy_yz_yy[i] = t_0_yyz_yz_yy[i] - rab_z[i] * t_0_yy_yz_yy[i];

        t_z_yy_yz_xz[i] = t_0_yyz_yz_xz[i] - rab_z[i] * t_0_yy_yz_xz[i];

        t_z_yy_yz_xy[i] = t_0_yyz_yz_xy[i] - rab_z[i] * t_0_yy_yz_xy[i];

        t_z_yy_yz_xx[i] = t_0_yyz_yz_xx[i] - rab_z[i] * t_0_yy_yz_xx[i];

        t_z_yy_yy_zz[i] = t_0_yyz_yy_zz[i] - rab_z[i] * t_0_yy_yy_zz[i];

        t_z_yy_yy_yz[i] = t_0_yyz_yy_yz[i] - rab_z[i] * t_0_yy_yy_yz[i];

        t_z_yy_yy_yy[i] = t_0_yyz_yy_yy[i] - rab_z[i] * t_0_yy_yy_yy[i];

        t_z_yy_yy_xz[i] = t_0_yyz_yy_xz[i] - rab_z[i] * t_0_yy_yy_xz[i];

        t_z_yy_yy_xy[i] = t_0_yyz_yy_xy[i] - rab_z[i] * t_0_yy_yy_xy[i];

        t_z_yy_yy_xx[i] = t_0_yyz_yy_xx[i] - rab_z[i] * t_0_yy_yy_xx[i];

        t_z_yy_xz_zz[i] = t_0_yyz_xz_zz[i] - rab_z[i] * t_0_yy_xz_zz[i];

        t_z_yy_xz_yz[i] = t_0_yyz_xz_yz[i] - rab_z[i] * t_0_yy_xz_yz[i];

        t_z_yy_xz_yy[i] = t_0_yyz_xz_yy[i] - rab_z[i] * t_0_yy_xz_yy[i];

        t_z_yy_xz_xz[i] = t_0_yyz_xz_xz[i] - rab_z[i] * t_0_yy_xz_xz[i];

        t_z_yy_xz_xy[i] = t_0_yyz_xz_xy[i] - rab_z[i] * t_0_yy_xz_xy[i];

        t_z_yy_xz_xx[i] = t_0_yyz_xz_xx[i] - rab_z[i] * t_0_yy_xz_xx[i];

        t_z_yy_xy_zz[i] = t_0_yyz_xy_zz[i] - rab_z[i] * t_0_yy_xy_zz[i];

        t_z_yy_xy_yz[i] = t_0_yyz_xy_yz[i] - rab_z[i] * t_0_yy_xy_yz[i];

        t_z_yy_xy_yy[i] = t_0_yyz_xy_yy[i] - rab_z[i] * t_0_yy_xy_yy[i];

        t_z_yy_xy_xz[i] = t_0_yyz_xy_xz[i] - rab_z[i] * t_0_yy_xy_xz[i];

        t_z_yy_xy_xy[i] = t_0_yyz_xy_xy[i] - rab_z[i] * t_0_yy_xy_xy[i];

        t_z_yy_xy_xx[i] = t_0_yyz_xy_xx[i] - rab_z[i] * t_0_yy_xy_xx[i];

        t_z_yy_xx_zz[i] = t_0_yyz_xx_zz[i] - rab_z[i] * t_0_yy_xx_zz[i];

        t_z_yy_xx_yz[i] = t_0_yyz_xx_yz[i] - rab_z[i] * t_0_yy_xx_yz[i];

        t_z_yy_xx_yy[i] = t_0_yyz_xx_yy[i] - rab_z[i] * t_0_yy_xx_yy[i];

        t_z_yy_xx_xz[i] = t_0_yyz_xx_xz[i] - rab_z[i] * t_0_yy_xx_xz[i];

        t_z_yy_xx_xy[i] = t_0_yyz_xx_xy[i] - rab_z[i] * t_0_yy_xx_xy[i];

        t_z_yy_xx_xx[i] = t_0_yyz_xx_xx[i] - rab_z[i] * t_0_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xz_xx_xx, t_0_xz_xx_xy, t_0_xz_xx_xz, t_0_xz_xx_yy,\
                           t_0_xz_xx_yz, t_0_xz_xx_zz, t_0_xz_xy_xx, t_0_xz_xy_xy, t_0_xz_xy_xz,\
                           t_0_xz_xy_yy, t_0_xz_xy_yz, t_0_xz_xy_zz, t_0_xz_xz_xx, t_0_xz_xz_xy,\
                           t_0_xz_xz_xz, t_0_xz_xz_yy, t_0_xz_xz_yz, t_0_xz_xz_zz, t_0_xz_yy_xx,\
                           t_0_xz_yy_xy, t_0_xz_yy_xz, t_0_xz_yy_yy, t_0_xz_yy_yz, t_0_xz_yy_zz,\
                           t_0_xz_yz_xx, t_0_xz_yz_xy, t_0_xz_yz_xz, t_0_xz_yz_yy, t_0_xz_yz_yz,\
                           t_0_xz_yz_zz, t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz, t_0_xz_zz_yy,\
                           t_0_xz_zz_yz, t_0_xz_zz_zz, t_0_xzz_xx_xx, t_0_xzz_xx_xy,\
                           t_0_xzz_xx_xz, t_0_xzz_xx_yy, t_0_xzz_xx_yz, t_0_xzz_xx_zz,\
                           t_0_xzz_xy_xx, t_0_xzz_xy_xy, t_0_xzz_xy_xz, t_0_xzz_xy_yy,\
                           t_0_xzz_xy_yz, t_0_xzz_xy_zz, t_0_xzz_xz_xx, t_0_xzz_xz_xy,\
                           t_0_xzz_xz_xz, t_0_xzz_xz_yy, t_0_xzz_xz_yz, t_0_xzz_xz_zz,\
                           t_0_xzz_yy_xx, t_0_xzz_yy_xy, t_0_xzz_yy_xz, t_0_xzz_yy_yy,\
                           t_0_xzz_yy_yz, t_0_xzz_yy_zz, t_0_xzz_yz_xx, t_0_xzz_yz_xy,\
                           t_0_xzz_yz_xz, t_0_xzz_yz_yy, t_0_xzz_yz_yz, t_0_xzz_yz_zz,\
                           t_0_xzz_zz_xx, t_0_xzz_zz_xy, t_0_xzz_zz_xz, t_0_xzz_zz_yy,\
                           t_0_xzz_zz_yz, t_0_xzz_zz_zz, t_z_xz_xx_xx, t_z_xz_xx_xy,\
                           t_z_xz_xx_xz, t_z_xz_xx_yy, t_z_xz_xx_yz, t_z_xz_xx_zz, t_z_xz_xy_xx,\
                           t_z_xz_xy_xy, t_z_xz_xy_xz, t_z_xz_xy_yy, t_z_xz_xy_yz, t_z_xz_xy_zz,\
                           t_z_xz_xz_xx, t_z_xz_xz_xy, t_z_xz_xz_xz, t_z_xz_xz_yy, t_z_xz_xz_yz,\
                           t_z_xz_xz_zz, t_z_xz_yy_xx, t_z_xz_yy_xy, t_z_xz_yy_xz, t_z_xz_yy_yy,\
                           t_z_xz_yy_yz, t_z_xz_yy_zz, t_z_xz_yz_xx, t_z_xz_yz_xy, t_z_xz_yz_xz,\
                           t_z_xz_yz_yy, t_z_xz_yz_yz, t_z_xz_yz_zz, t_z_xz_zz_xx, t_z_xz_zz_xy,\
                           t_z_xz_zz_xz, t_z_xz_zz_yy, t_z_xz_zz_yz, t_z_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xz_zz_zz[i] = t_0_xzz_zz_zz[i] - rab_z[i] * t_0_xz_zz_zz[i];

        t_z_xz_zz_yz[i] = t_0_xzz_zz_yz[i] - rab_z[i] * t_0_xz_zz_yz[i];

        t_z_xz_zz_yy[i] = t_0_xzz_zz_yy[i] - rab_z[i] * t_0_xz_zz_yy[i];

        t_z_xz_zz_xz[i] = t_0_xzz_zz_xz[i] - rab_z[i] * t_0_xz_zz_xz[i];

        t_z_xz_zz_xy[i] = t_0_xzz_zz_xy[i] - rab_z[i] * t_0_xz_zz_xy[i];

        t_z_xz_zz_xx[i] = t_0_xzz_zz_xx[i] - rab_z[i] * t_0_xz_zz_xx[i];

        t_z_xz_yz_zz[i] = t_0_xzz_yz_zz[i] - rab_z[i] * t_0_xz_yz_zz[i];

        t_z_xz_yz_yz[i] = t_0_xzz_yz_yz[i] - rab_z[i] * t_0_xz_yz_yz[i];

        t_z_xz_yz_yy[i] = t_0_xzz_yz_yy[i] - rab_z[i] * t_0_xz_yz_yy[i];

        t_z_xz_yz_xz[i] = t_0_xzz_yz_xz[i] - rab_z[i] * t_0_xz_yz_xz[i];

        t_z_xz_yz_xy[i] = t_0_xzz_yz_xy[i] - rab_z[i] * t_0_xz_yz_xy[i];

        t_z_xz_yz_xx[i] = t_0_xzz_yz_xx[i] - rab_z[i] * t_0_xz_yz_xx[i];

        t_z_xz_yy_zz[i] = t_0_xzz_yy_zz[i] - rab_z[i] * t_0_xz_yy_zz[i];

        t_z_xz_yy_yz[i] = t_0_xzz_yy_yz[i] - rab_z[i] * t_0_xz_yy_yz[i];

        t_z_xz_yy_yy[i] = t_0_xzz_yy_yy[i] - rab_z[i] * t_0_xz_yy_yy[i];

        t_z_xz_yy_xz[i] = t_0_xzz_yy_xz[i] - rab_z[i] * t_0_xz_yy_xz[i];

        t_z_xz_yy_xy[i] = t_0_xzz_yy_xy[i] - rab_z[i] * t_0_xz_yy_xy[i];

        t_z_xz_yy_xx[i] = t_0_xzz_yy_xx[i] - rab_z[i] * t_0_xz_yy_xx[i];

        t_z_xz_xz_zz[i] = t_0_xzz_xz_zz[i] - rab_z[i] * t_0_xz_xz_zz[i];

        t_z_xz_xz_yz[i] = t_0_xzz_xz_yz[i] - rab_z[i] * t_0_xz_xz_yz[i];

        t_z_xz_xz_yy[i] = t_0_xzz_xz_yy[i] - rab_z[i] * t_0_xz_xz_yy[i];

        t_z_xz_xz_xz[i] = t_0_xzz_xz_xz[i] - rab_z[i] * t_0_xz_xz_xz[i];

        t_z_xz_xz_xy[i] = t_0_xzz_xz_xy[i] - rab_z[i] * t_0_xz_xz_xy[i];

        t_z_xz_xz_xx[i] = t_0_xzz_xz_xx[i] - rab_z[i] * t_0_xz_xz_xx[i];

        t_z_xz_xy_zz[i] = t_0_xzz_xy_zz[i] - rab_z[i] * t_0_xz_xy_zz[i];

        t_z_xz_xy_yz[i] = t_0_xzz_xy_yz[i] - rab_z[i] * t_0_xz_xy_yz[i];

        t_z_xz_xy_yy[i] = t_0_xzz_xy_yy[i] - rab_z[i] * t_0_xz_xy_yy[i];

        t_z_xz_xy_xz[i] = t_0_xzz_xy_xz[i] - rab_z[i] * t_0_xz_xy_xz[i];

        t_z_xz_xy_xy[i] = t_0_xzz_xy_xy[i] - rab_z[i] * t_0_xz_xy_xy[i];

        t_z_xz_xy_xx[i] = t_0_xzz_xy_xx[i] - rab_z[i] * t_0_xz_xy_xx[i];

        t_z_xz_xx_zz[i] = t_0_xzz_xx_zz[i] - rab_z[i] * t_0_xz_xx_zz[i];

        t_z_xz_xx_yz[i] = t_0_xzz_xx_yz[i] - rab_z[i] * t_0_xz_xx_yz[i];

        t_z_xz_xx_yy[i] = t_0_xzz_xx_yy[i] - rab_z[i] * t_0_xz_xx_yy[i];

        t_z_xz_xx_xz[i] = t_0_xzz_xx_xz[i] - rab_z[i] * t_0_xz_xx_xz[i];

        t_z_xz_xx_xy[i] = t_0_xzz_xx_xy[i] - rab_z[i] * t_0_xz_xx_xy[i];

        t_z_xz_xx_xx[i] = t_0_xzz_xx_xx[i] - rab_z[i] * t_0_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xy_xx_xx, t_0_xy_xx_xy, t_0_xy_xx_xz, t_0_xy_xx_yy,\
                           t_0_xy_xx_yz, t_0_xy_xx_zz, t_0_xy_xy_xx, t_0_xy_xy_xy, t_0_xy_xy_xz,\
                           t_0_xy_xy_yy, t_0_xy_xy_yz, t_0_xy_xy_zz, t_0_xy_xz_xx, t_0_xy_xz_xy,\
                           t_0_xy_xz_xz, t_0_xy_xz_yy, t_0_xy_xz_yz, t_0_xy_xz_zz, t_0_xy_yy_xx,\
                           t_0_xy_yy_xy, t_0_xy_yy_xz, t_0_xy_yy_yy, t_0_xy_yy_yz, t_0_xy_yy_zz,\
                           t_0_xy_yz_xx, t_0_xy_yz_xy, t_0_xy_yz_xz, t_0_xy_yz_yy, t_0_xy_yz_yz,\
                           t_0_xy_yz_zz, t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz, t_0_xy_zz_yy,\
                           t_0_xy_zz_yz, t_0_xy_zz_zz, t_0_xyz_xx_xx, t_0_xyz_xx_xy,\
                           t_0_xyz_xx_xz, t_0_xyz_xx_yy, t_0_xyz_xx_yz, t_0_xyz_xx_zz,\
                           t_0_xyz_xy_xx, t_0_xyz_xy_xy, t_0_xyz_xy_xz, t_0_xyz_xy_yy,\
                           t_0_xyz_xy_yz, t_0_xyz_xy_zz, t_0_xyz_xz_xx, t_0_xyz_xz_xy,\
                           t_0_xyz_xz_xz, t_0_xyz_xz_yy, t_0_xyz_xz_yz, t_0_xyz_xz_zz,\
                           t_0_xyz_yy_xx, t_0_xyz_yy_xy, t_0_xyz_yy_xz, t_0_xyz_yy_yy,\
                           t_0_xyz_yy_yz, t_0_xyz_yy_zz, t_0_xyz_yz_xx, t_0_xyz_yz_xy,\
                           t_0_xyz_yz_xz, t_0_xyz_yz_yy, t_0_xyz_yz_yz, t_0_xyz_yz_zz,\
                           t_0_xyz_zz_xx, t_0_xyz_zz_xy, t_0_xyz_zz_xz, t_0_xyz_zz_yy,\
                           t_0_xyz_zz_yz, t_0_xyz_zz_zz, t_z_xy_xx_xx, t_z_xy_xx_xy,\
                           t_z_xy_xx_xz, t_z_xy_xx_yy, t_z_xy_xx_yz, t_z_xy_xx_zz, t_z_xy_xy_xx,\
                           t_z_xy_xy_xy, t_z_xy_xy_xz, t_z_xy_xy_yy, t_z_xy_xy_yz, t_z_xy_xy_zz,\
                           t_z_xy_xz_xx, t_z_xy_xz_xy, t_z_xy_xz_xz, t_z_xy_xz_yy, t_z_xy_xz_yz,\
                           t_z_xy_xz_zz, t_z_xy_yy_xx, t_z_xy_yy_xy, t_z_xy_yy_xz, t_z_xy_yy_yy,\
                           t_z_xy_yy_yz, t_z_xy_yy_zz, t_z_xy_yz_xx, t_z_xy_yz_xy, t_z_xy_yz_xz,\
                           t_z_xy_yz_yy, t_z_xy_yz_yz, t_z_xy_yz_zz, t_z_xy_zz_xx, t_z_xy_zz_xy,\
                           t_z_xy_zz_xz, t_z_xy_zz_yy, t_z_xy_zz_yz, t_z_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xy_zz_zz[i] = t_0_xyz_zz_zz[i] - rab_z[i] * t_0_xy_zz_zz[i];

        t_z_xy_zz_yz[i] = t_0_xyz_zz_yz[i] - rab_z[i] * t_0_xy_zz_yz[i];

        t_z_xy_zz_yy[i] = t_0_xyz_zz_yy[i] - rab_z[i] * t_0_xy_zz_yy[i];

        t_z_xy_zz_xz[i] = t_0_xyz_zz_xz[i] - rab_z[i] * t_0_xy_zz_xz[i];

        t_z_xy_zz_xy[i] = t_0_xyz_zz_xy[i] - rab_z[i] * t_0_xy_zz_xy[i];

        t_z_xy_zz_xx[i] = t_0_xyz_zz_xx[i] - rab_z[i] * t_0_xy_zz_xx[i];

        t_z_xy_yz_zz[i] = t_0_xyz_yz_zz[i] - rab_z[i] * t_0_xy_yz_zz[i];

        t_z_xy_yz_yz[i] = t_0_xyz_yz_yz[i] - rab_z[i] * t_0_xy_yz_yz[i];

        t_z_xy_yz_yy[i] = t_0_xyz_yz_yy[i] - rab_z[i] * t_0_xy_yz_yy[i];

        t_z_xy_yz_xz[i] = t_0_xyz_yz_xz[i] - rab_z[i] * t_0_xy_yz_xz[i];

        t_z_xy_yz_xy[i] = t_0_xyz_yz_xy[i] - rab_z[i] * t_0_xy_yz_xy[i];

        t_z_xy_yz_xx[i] = t_0_xyz_yz_xx[i] - rab_z[i] * t_0_xy_yz_xx[i];

        t_z_xy_yy_zz[i] = t_0_xyz_yy_zz[i] - rab_z[i] * t_0_xy_yy_zz[i];

        t_z_xy_yy_yz[i] = t_0_xyz_yy_yz[i] - rab_z[i] * t_0_xy_yy_yz[i];

        t_z_xy_yy_yy[i] = t_0_xyz_yy_yy[i] - rab_z[i] * t_0_xy_yy_yy[i];

        t_z_xy_yy_xz[i] = t_0_xyz_yy_xz[i] - rab_z[i] * t_0_xy_yy_xz[i];

        t_z_xy_yy_xy[i] = t_0_xyz_yy_xy[i] - rab_z[i] * t_0_xy_yy_xy[i];

        t_z_xy_yy_xx[i] = t_0_xyz_yy_xx[i] - rab_z[i] * t_0_xy_yy_xx[i];

        t_z_xy_xz_zz[i] = t_0_xyz_xz_zz[i] - rab_z[i] * t_0_xy_xz_zz[i];

        t_z_xy_xz_yz[i] = t_0_xyz_xz_yz[i] - rab_z[i] * t_0_xy_xz_yz[i];

        t_z_xy_xz_yy[i] = t_0_xyz_xz_yy[i] - rab_z[i] * t_0_xy_xz_yy[i];

        t_z_xy_xz_xz[i] = t_0_xyz_xz_xz[i] - rab_z[i] * t_0_xy_xz_xz[i];

        t_z_xy_xz_xy[i] = t_0_xyz_xz_xy[i] - rab_z[i] * t_0_xy_xz_xy[i];

        t_z_xy_xz_xx[i] = t_0_xyz_xz_xx[i] - rab_z[i] * t_0_xy_xz_xx[i];

        t_z_xy_xy_zz[i] = t_0_xyz_xy_zz[i] - rab_z[i] * t_0_xy_xy_zz[i];

        t_z_xy_xy_yz[i] = t_0_xyz_xy_yz[i] - rab_z[i] * t_0_xy_xy_yz[i];

        t_z_xy_xy_yy[i] = t_0_xyz_xy_yy[i] - rab_z[i] * t_0_xy_xy_yy[i];

        t_z_xy_xy_xz[i] = t_0_xyz_xy_xz[i] - rab_z[i] * t_0_xy_xy_xz[i];

        t_z_xy_xy_xy[i] = t_0_xyz_xy_xy[i] - rab_z[i] * t_0_xy_xy_xy[i];

        t_z_xy_xy_xx[i] = t_0_xyz_xy_xx[i] - rab_z[i] * t_0_xy_xy_xx[i];

        t_z_xy_xx_zz[i] = t_0_xyz_xx_zz[i] - rab_z[i] * t_0_xy_xx_zz[i];

        t_z_xy_xx_yz[i] = t_0_xyz_xx_yz[i] - rab_z[i] * t_0_xy_xx_yz[i];

        t_z_xy_xx_yy[i] = t_0_xyz_xx_yy[i] - rab_z[i] * t_0_xy_xx_yy[i];

        t_z_xy_xx_xz[i] = t_0_xyz_xx_xz[i] - rab_z[i] * t_0_xy_xx_xz[i];

        t_z_xy_xx_xy[i] = t_0_xyz_xx_xy[i] - rab_z[i] * t_0_xy_xx_xy[i];

        t_z_xy_xx_xx[i] = t_0_xyz_xx_xx[i] - rab_z[i] * t_0_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xx_xx_xx, t_0_xx_xx_xy, t_0_xx_xx_xz, t_0_xx_xx_yy,\
                           t_0_xx_xx_yz, t_0_xx_xx_zz, t_0_xx_xy_xx, t_0_xx_xy_xy, t_0_xx_xy_xz,\
                           t_0_xx_xy_yy, t_0_xx_xy_yz, t_0_xx_xy_zz, t_0_xx_xz_xx, t_0_xx_xz_xy,\
                           t_0_xx_xz_xz, t_0_xx_xz_yy, t_0_xx_xz_yz, t_0_xx_xz_zz, t_0_xx_yy_xx,\
                           t_0_xx_yy_xy, t_0_xx_yy_xz, t_0_xx_yy_yy, t_0_xx_yy_yz, t_0_xx_yy_zz,\
                           t_0_xx_yz_xx, t_0_xx_yz_xy, t_0_xx_yz_xz, t_0_xx_yz_yy, t_0_xx_yz_yz,\
                           t_0_xx_yz_zz, t_0_xx_zz_xx, t_0_xx_zz_xy, t_0_xx_zz_xz, t_0_xx_zz_yy,\
                           t_0_xx_zz_yz, t_0_xx_zz_zz, t_0_xxz_xx_xx, t_0_xxz_xx_xy,\
                           t_0_xxz_xx_xz, t_0_xxz_xx_yy, t_0_xxz_xx_yz, t_0_xxz_xx_zz,\
                           t_0_xxz_xy_xx, t_0_xxz_xy_xy, t_0_xxz_xy_xz, t_0_xxz_xy_yy,\
                           t_0_xxz_xy_yz, t_0_xxz_xy_zz, t_0_xxz_xz_xx, t_0_xxz_xz_xy,\
                           t_0_xxz_xz_xz, t_0_xxz_xz_yy, t_0_xxz_xz_yz, t_0_xxz_xz_zz,\
                           t_0_xxz_yy_xx, t_0_xxz_yy_xy, t_0_xxz_yy_xz, t_0_xxz_yy_yy,\
                           t_0_xxz_yy_yz, t_0_xxz_yy_zz, t_0_xxz_yz_xx, t_0_xxz_yz_xy,\
                           t_0_xxz_yz_xz, t_0_xxz_yz_yy, t_0_xxz_yz_yz, t_0_xxz_yz_zz,\
                           t_0_xxz_zz_xx, t_0_xxz_zz_xy, t_0_xxz_zz_xz, t_0_xxz_zz_yy,\
                           t_0_xxz_zz_yz, t_0_xxz_zz_zz, t_z_xx_xx_xx, t_z_xx_xx_xy,\
                           t_z_xx_xx_xz, t_z_xx_xx_yy, t_z_xx_xx_yz, t_z_xx_xx_zz, t_z_xx_xy_xx,\
                           t_z_xx_xy_xy, t_z_xx_xy_xz, t_z_xx_xy_yy, t_z_xx_xy_yz, t_z_xx_xy_zz,\
                           t_z_xx_xz_xx, t_z_xx_xz_xy, t_z_xx_xz_xz, t_z_xx_xz_yy, t_z_xx_xz_yz,\
                           t_z_xx_xz_zz, t_z_xx_yy_xx, t_z_xx_yy_xy, t_z_xx_yy_xz, t_z_xx_yy_yy,\
                           t_z_xx_yy_yz, t_z_xx_yy_zz, t_z_xx_yz_xx, t_z_xx_yz_xy, t_z_xx_yz_xz,\
                           t_z_xx_yz_yy, t_z_xx_yz_yz, t_z_xx_yz_zz, t_z_xx_zz_xx, t_z_xx_zz_xy,\
                           t_z_xx_zz_xz, t_z_xx_zz_yy, t_z_xx_zz_yz, t_z_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xx_zz_zz[i] = t_0_xxz_zz_zz[i] - rab_z[i] * t_0_xx_zz_zz[i];

        t_z_xx_zz_yz[i] = t_0_xxz_zz_yz[i] - rab_z[i] * t_0_xx_zz_yz[i];

        t_z_xx_zz_yy[i] = t_0_xxz_zz_yy[i] - rab_z[i] * t_0_xx_zz_yy[i];

        t_z_xx_zz_xz[i] = t_0_xxz_zz_xz[i] - rab_z[i] * t_0_xx_zz_xz[i];

        t_z_xx_zz_xy[i] = t_0_xxz_zz_xy[i] - rab_z[i] * t_0_xx_zz_xy[i];

        t_z_xx_zz_xx[i] = t_0_xxz_zz_xx[i] - rab_z[i] * t_0_xx_zz_xx[i];

        t_z_xx_yz_zz[i] = t_0_xxz_yz_zz[i] - rab_z[i] * t_0_xx_yz_zz[i];

        t_z_xx_yz_yz[i] = t_0_xxz_yz_yz[i] - rab_z[i] * t_0_xx_yz_yz[i];

        t_z_xx_yz_yy[i] = t_0_xxz_yz_yy[i] - rab_z[i] * t_0_xx_yz_yy[i];

        t_z_xx_yz_xz[i] = t_0_xxz_yz_xz[i] - rab_z[i] * t_0_xx_yz_xz[i];

        t_z_xx_yz_xy[i] = t_0_xxz_yz_xy[i] - rab_z[i] * t_0_xx_yz_xy[i];

        t_z_xx_yz_xx[i] = t_0_xxz_yz_xx[i] - rab_z[i] * t_0_xx_yz_xx[i];

        t_z_xx_yy_zz[i] = t_0_xxz_yy_zz[i] - rab_z[i] * t_0_xx_yy_zz[i];

        t_z_xx_yy_yz[i] = t_0_xxz_yy_yz[i] - rab_z[i] * t_0_xx_yy_yz[i];

        t_z_xx_yy_yy[i] = t_0_xxz_yy_yy[i] - rab_z[i] * t_0_xx_yy_yy[i];

        t_z_xx_yy_xz[i] = t_0_xxz_yy_xz[i] - rab_z[i] * t_0_xx_yy_xz[i];

        t_z_xx_yy_xy[i] = t_0_xxz_yy_xy[i] - rab_z[i] * t_0_xx_yy_xy[i];

        t_z_xx_yy_xx[i] = t_0_xxz_yy_xx[i] - rab_z[i] * t_0_xx_yy_xx[i];

        t_z_xx_xz_zz[i] = t_0_xxz_xz_zz[i] - rab_z[i] * t_0_xx_xz_zz[i];

        t_z_xx_xz_yz[i] = t_0_xxz_xz_yz[i] - rab_z[i] * t_0_xx_xz_yz[i];

        t_z_xx_xz_yy[i] = t_0_xxz_xz_yy[i] - rab_z[i] * t_0_xx_xz_yy[i];

        t_z_xx_xz_xz[i] = t_0_xxz_xz_xz[i] - rab_z[i] * t_0_xx_xz_xz[i];

        t_z_xx_xz_xy[i] = t_0_xxz_xz_xy[i] - rab_z[i] * t_0_xx_xz_xy[i];

        t_z_xx_xz_xx[i] = t_0_xxz_xz_xx[i] - rab_z[i] * t_0_xx_xz_xx[i];

        t_z_xx_xy_zz[i] = t_0_xxz_xy_zz[i] - rab_z[i] * t_0_xx_xy_zz[i];

        t_z_xx_xy_yz[i] = t_0_xxz_xy_yz[i] - rab_z[i] * t_0_xx_xy_yz[i];

        t_z_xx_xy_yy[i] = t_0_xxz_xy_yy[i] - rab_z[i] * t_0_xx_xy_yy[i];

        t_z_xx_xy_xz[i] = t_0_xxz_xy_xz[i] - rab_z[i] * t_0_xx_xy_xz[i];

        t_z_xx_xy_xy[i] = t_0_xxz_xy_xy[i] - rab_z[i] * t_0_xx_xy_xy[i];

        t_z_xx_xy_xx[i] = t_0_xxz_xy_xx[i] - rab_z[i] * t_0_xx_xy_xx[i];

        t_z_xx_xx_zz[i] = t_0_xxz_xx_zz[i] - rab_z[i] * t_0_xx_xx_zz[i];

        t_z_xx_xx_yz[i] = t_0_xxz_xx_yz[i] - rab_z[i] * t_0_xx_xx_yz[i];

        t_z_xx_xx_yy[i] = t_0_xxz_xx_yy[i] - rab_z[i] * t_0_xx_xx_yy[i];

        t_z_xx_xx_xz[i] = t_0_xxz_xx_xz[i] - rab_z[i] * t_0_xx_xx_xz[i];

        t_z_xx_xx_xy[i] = t_0_xxz_xx_xy[i] - rab_z[i] * t_0_xx_xx_xy[i];

        t_z_xx_xx_xx[i] = t_0_xxz_xx_xx[i] - rab_z[i] * t_0_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yzz_xx_xx, t_0_yzz_xx_xy, t_0_yzz_xx_xz, t_0_yzz_xx_yy,\
                           t_0_yzz_xx_yz, t_0_yzz_xx_zz, t_0_yzz_xy_xx, t_0_yzz_xy_xy,\
                           t_0_yzz_xy_xz, t_0_yzz_xy_yy, t_0_yzz_xy_yz, t_0_yzz_xy_zz,\
                           t_0_yzz_xz_xx, t_0_yzz_xz_xy, t_0_yzz_xz_xz, t_0_yzz_xz_yy,\
                           t_0_yzz_xz_yz, t_0_yzz_xz_zz, t_0_yzz_yy_xx, t_0_yzz_yy_xy,\
                           t_0_yzz_yy_xz, t_0_yzz_yy_yy, t_0_yzz_yy_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_xx, t_0_yzz_yz_xy, t_0_yzz_yz_xz, t_0_yzz_yz_yy,\
                           t_0_yzz_yz_yz, t_0_yzz_yz_zz, t_0_yzz_zz_xx, t_0_yzz_zz_xy,\
                           t_0_yzz_zz_xz, t_0_yzz_zz_yy, t_0_yzz_zz_yz, t_0_yzz_zz_zz,\
                           t_0_zz_xx_xx, t_0_zz_xx_xy, t_0_zz_xx_xz, t_0_zz_xx_yy, t_0_zz_xx_yz,\
                           t_0_zz_xx_zz, t_0_zz_xy_xx, t_0_zz_xy_xy, t_0_zz_xy_xz, t_0_zz_xy_yy,\
                           t_0_zz_xy_yz, t_0_zz_xy_zz, t_0_zz_xz_xx, t_0_zz_xz_xy, t_0_zz_xz_xz,\
                           t_0_zz_xz_yy, t_0_zz_xz_yz, t_0_zz_xz_zz, t_0_zz_yy_xx, t_0_zz_yy_xy,\
                           t_0_zz_yy_xz, t_0_zz_yy_yy, t_0_zz_yy_yz, t_0_zz_yy_zz, t_0_zz_yz_xx,\
                           t_0_zz_yz_xy, t_0_zz_yz_xz, t_0_zz_yz_yy, t_0_zz_yz_yz, t_0_zz_yz_zz,\
                           t_0_zz_zz_xx, t_0_zz_zz_xy, t_0_zz_zz_xz, t_0_zz_zz_yy, t_0_zz_zz_yz,\
                           t_0_zz_zz_zz, t_y_zz_xx_xx, t_y_zz_xx_xy, t_y_zz_xx_xz, t_y_zz_xx_yy,\
                           t_y_zz_xx_yz, t_y_zz_xx_zz, t_y_zz_xy_xx, t_y_zz_xy_xy, t_y_zz_xy_xz,\
                           t_y_zz_xy_yy, t_y_zz_xy_yz, t_y_zz_xy_zz, t_y_zz_xz_xx, t_y_zz_xz_xy,\
                           t_y_zz_xz_xz, t_y_zz_xz_yy, t_y_zz_xz_yz, t_y_zz_xz_zz, t_y_zz_yy_xx,\
                           t_y_zz_yy_xy, t_y_zz_yy_xz, t_y_zz_yy_yy, t_y_zz_yy_yz, t_y_zz_yy_zz,\
                           t_y_zz_yz_xx, t_y_zz_yz_xy, t_y_zz_yz_xz, t_y_zz_yz_yy, t_y_zz_yz_yz,\
                           t_y_zz_yz_zz, t_y_zz_zz_xx, t_y_zz_zz_xy, t_y_zz_zz_xz, t_y_zz_zz_yy,\
                           t_y_zz_zz_yz, t_y_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zz_zz_zz[i] = t_0_yzz_zz_zz[i] - rab_y[i] * t_0_zz_zz_zz[i];

        t_y_zz_zz_yz[i] = t_0_yzz_zz_yz[i] - rab_y[i] * t_0_zz_zz_yz[i];

        t_y_zz_zz_yy[i] = t_0_yzz_zz_yy[i] - rab_y[i] * t_0_zz_zz_yy[i];

        t_y_zz_zz_xz[i] = t_0_yzz_zz_xz[i] - rab_y[i] * t_0_zz_zz_xz[i];

        t_y_zz_zz_xy[i] = t_0_yzz_zz_xy[i] - rab_y[i] * t_0_zz_zz_xy[i];

        t_y_zz_zz_xx[i] = t_0_yzz_zz_xx[i] - rab_y[i] * t_0_zz_zz_xx[i];

        t_y_zz_yz_zz[i] = t_0_yzz_yz_zz[i] - rab_y[i] * t_0_zz_yz_zz[i];

        t_y_zz_yz_yz[i] = t_0_yzz_yz_yz[i] - rab_y[i] * t_0_zz_yz_yz[i];

        t_y_zz_yz_yy[i] = t_0_yzz_yz_yy[i] - rab_y[i] * t_0_zz_yz_yy[i];

        t_y_zz_yz_xz[i] = t_0_yzz_yz_xz[i] - rab_y[i] * t_0_zz_yz_xz[i];

        t_y_zz_yz_xy[i] = t_0_yzz_yz_xy[i] - rab_y[i] * t_0_zz_yz_xy[i];

        t_y_zz_yz_xx[i] = t_0_yzz_yz_xx[i] - rab_y[i] * t_0_zz_yz_xx[i];

        t_y_zz_yy_zz[i] = t_0_yzz_yy_zz[i] - rab_y[i] * t_0_zz_yy_zz[i];

        t_y_zz_yy_yz[i] = t_0_yzz_yy_yz[i] - rab_y[i] * t_0_zz_yy_yz[i];

        t_y_zz_yy_yy[i] = t_0_yzz_yy_yy[i] - rab_y[i] * t_0_zz_yy_yy[i];

        t_y_zz_yy_xz[i] = t_0_yzz_yy_xz[i] - rab_y[i] * t_0_zz_yy_xz[i];

        t_y_zz_yy_xy[i] = t_0_yzz_yy_xy[i] - rab_y[i] * t_0_zz_yy_xy[i];

        t_y_zz_yy_xx[i] = t_0_yzz_yy_xx[i] - rab_y[i] * t_0_zz_yy_xx[i];

        t_y_zz_xz_zz[i] = t_0_yzz_xz_zz[i] - rab_y[i] * t_0_zz_xz_zz[i];

        t_y_zz_xz_yz[i] = t_0_yzz_xz_yz[i] - rab_y[i] * t_0_zz_xz_yz[i];

        t_y_zz_xz_yy[i] = t_0_yzz_xz_yy[i] - rab_y[i] * t_0_zz_xz_yy[i];

        t_y_zz_xz_xz[i] = t_0_yzz_xz_xz[i] - rab_y[i] * t_0_zz_xz_xz[i];

        t_y_zz_xz_xy[i] = t_0_yzz_xz_xy[i] - rab_y[i] * t_0_zz_xz_xy[i];

        t_y_zz_xz_xx[i] = t_0_yzz_xz_xx[i] - rab_y[i] * t_0_zz_xz_xx[i];

        t_y_zz_xy_zz[i] = t_0_yzz_xy_zz[i] - rab_y[i] * t_0_zz_xy_zz[i];

        t_y_zz_xy_yz[i] = t_0_yzz_xy_yz[i] - rab_y[i] * t_0_zz_xy_yz[i];

        t_y_zz_xy_yy[i] = t_0_yzz_xy_yy[i] - rab_y[i] * t_0_zz_xy_yy[i];

        t_y_zz_xy_xz[i] = t_0_yzz_xy_xz[i] - rab_y[i] * t_0_zz_xy_xz[i];

        t_y_zz_xy_xy[i] = t_0_yzz_xy_xy[i] - rab_y[i] * t_0_zz_xy_xy[i];

        t_y_zz_xy_xx[i] = t_0_yzz_xy_xx[i] - rab_y[i] * t_0_zz_xy_xx[i];

        t_y_zz_xx_zz[i] = t_0_yzz_xx_zz[i] - rab_y[i] * t_0_zz_xx_zz[i];

        t_y_zz_xx_yz[i] = t_0_yzz_xx_yz[i] - rab_y[i] * t_0_zz_xx_yz[i];

        t_y_zz_xx_yy[i] = t_0_yzz_xx_yy[i] - rab_y[i] * t_0_zz_xx_yy[i];

        t_y_zz_xx_xz[i] = t_0_yzz_xx_xz[i] - rab_y[i] * t_0_zz_xx_xz[i];

        t_y_zz_xx_xy[i] = t_0_yzz_xx_xy[i] - rab_y[i] * t_0_zz_xx_xy[i];

        t_y_zz_xx_xx[i] = t_0_yzz_xx_xx[i] - rab_y[i] * t_0_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyz_xx_xx, t_0_yyz_xx_xy, t_0_yyz_xx_xz, t_0_yyz_xx_yy,\
                           t_0_yyz_xx_yz, t_0_yyz_xx_zz, t_0_yyz_xy_xx, t_0_yyz_xy_xy,\
                           t_0_yyz_xy_xz, t_0_yyz_xy_yy, t_0_yyz_xy_yz, t_0_yyz_xy_zz,\
                           t_0_yyz_xz_xx, t_0_yyz_xz_xy, t_0_yyz_xz_xz, t_0_yyz_xz_yy,\
                           t_0_yyz_xz_yz, t_0_yyz_xz_zz, t_0_yyz_yy_xx, t_0_yyz_yy_xy,\
                           t_0_yyz_yy_xz, t_0_yyz_yy_yy, t_0_yyz_yy_yz, t_0_yyz_yy_zz,\
                           t_0_yyz_yz_xx, t_0_yyz_yz_xy, t_0_yyz_yz_xz, t_0_yyz_yz_yy,\
                           t_0_yyz_yz_yz, t_0_yyz_yz_zz, t_0_yyz_zz_xx, t_0_yyz_zz_xy,\
                           t_0_yyz_zz_xz, t_0_yyz_zz_yy, t_0_yyz_zz_yz, t_0_yyz_zz_zz,\
                           t_0_yz_xx_xx, t_0_yz_xx_xy, t_0_yz_xx_xz, t_0_yz_xx_yy, t_0_yz_xx_yz,\
                           t_0_yz_xx_zz, t_0_yz_xy_xx, t_0_yz_xy_xy, t_0_yz_xy_xz, t_0_yz_xy_yy,\
                           t_0_yz_xy_yz, t_0_yz_xy_zz, t_0_yz_xz_xx, t_0_yz_xz_xy, t_0_yz_xz_xz,\
                           t_0_yz_xz_yy, t_0_yz_xz_yz, t_0_yz_xz_zz, t_0_yz_yy_xx, t_0_yz_yy_xy,\
                           t_0_yz_yy_xz, t_0_yz_yy_yy, t_0_yz_yy_yz, t_0_yz_yy_zz, t_0_yz_yz_xx,\
                           t_0_yz_yz_xy, t_0_yz_yz_xz, t_0_yz_yz_yy, t_0_yz_yz_yz, t_0_yz_yz_zz,\
                           t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz, t_0_yz_zz_yy, t_0_yz_zz_yz,\
                           t_0_yz_zz_zz, t_y_yz_xx_xx, t_y_yz_xx_xy, t_y_yz_xx_xz, t_y_yz_xx_yy,\
                           t_y_yz_xx_yz, t_y_yz_xx_zz, t_y_yz_xy_xx, t_y_yz_xy_xy, t_y_yz_xy_xz,\
                           t_y_yz_xy_yy, t_y_yz_xy_yz, t_y_yz_xy_zz, t_y_yz_xz_xx, t_y_yz_xz_xy,\
                           t_y_yz_xz_xz, t_y_yz_xz_yy, t_y_yz_xz_yz, t_y_yz_xz_zz, t_y_yz_yy_xx,\
                           t_y_yz_yy_xy, t_y_yz_yy_xz, t_y_yz_yy_yy, t_y_yz_yy_yz, t_y_yz_yy_zz,\
                           t_y_yz_yz_xx, t_y_yz_yz_xy, t_y_yz_yz_xz, t_y_yz_yz_yy, t_y_yz_yz_yz,\
                           t_y_yz_yz_zz, t_y_yz_zz_xx, t_y_yz_zz_xy, t_y_yz_zz_xz, t_y_yz_zz_yy,\
                           t_y_yz_zz_yz, t_y_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yz_zz_zz[i] = t_0_yyz_zz_zz[i] - rab_y[i] * t_0_yz_zz_zz[i];

        t_y_yz_zz_yz[i] = t_0_yyz_zz_yz[i] - rab_y[i] * t_0_yz_zz_yz[i];

        t_y_yz_zz_yy[i] = t_0_yyz_zz_yy[i] - rab_y[i] * t_0_yz_zz_yy[i];

        t_y_yz_zz_xz[i] = t_0_yyz_zz_xz[i] - rab_y[i] * t_0_yz_zz_xz[i];

        t_y_yz_zz_xy[i] = t_0_yyz_zz_xy[i] - rab_y[i] * t_0_yz_zz_xy[i];

        t_y_yz_zz_xx[i] = t_0_yyz_zz_xx[i] - rab_y[i] * t_0_yz_zz_xx[i];

        t_y_yz_yz_zz[i] = t_0_yyz_yz_zz[i] - rab_y[i] * t_0_yz_yz_zz[i];

        t_y_yz_yz_yz[i] = t_0_yyz_yz_yz[i] - rab_y[i] * t_0_yz_yz_yz[i];

        t_y_yz_yz_yy[i] = t_0_yyz_yz_yy[i] - rab_y[i] * t_0_yz_yz_yy[i];

        t_y_yz_yz_xz[i] = t_0_yyz_yz_xz[i] - rab_y[i] * t_0_yz_yz_xz[i];

        t_y_yz_yz_xy[i] = t_0_yyz_yz_xy[i] - rab_y[i] * t_0_yz_yz_xy[i];

        t_y_yz_yz_xx[i] = t_0_yyz_yz_xx[i] - rab_y[i] * t_0_yz_yz_xx[i];

        t_y_yz_yy_zz[i] = t_0_yyz_yy_zz[i] - rab_y[i] * t_0_yz_yy_zz[i];

        t_y_yz_yy_yz[i] = t_0_yyz_yy_yz[i] - rab_y[i] * t_0_yz_yy_yz[i];

        t_y_yz_yy_yy[i] = t_0_yyz_yy_yy[i] - rab_y[i] * t_0_yz_yy_yy[i];

        t_y_yz_yy_xz[i] = t_0_yyz_yy_xz[i] - rab_y[i] * t_0_yz_yy_xz[i];

        t_y_yz_yy_xy[i] = t_0_yyz_yy_xy[i] - rab_y[i] * t_0_yz_yy_xy[i];

        t_y_yz_yy_xx[i] = t_0_yyz_yy_xx[i] - rab_y[i] * t_0_yz_yy_xx[i];

        t_y_yz_xz_zz[i] = t_0_yyz_xz_zz[i] - rab_y[i] * t_0_yz_xz_zz[i];

        t_y_yz_xz_yz[i] = t_0_yyz_xz_yz[i] - rab_y[i] * t_0_yz_xz_yz[i];

        t_y_yz_xz_yy[i] = t_0_yyz_xz_yy[i] - rab_y[i] * t_0_yz_xz_yy[i];

        t_y_yz_xz_xz[i] = t_0_yyz_xz_xz[i] - rab_y[i] * t_0_yz_xz_xz[i];

        t_y_yz_xz_xy[i] = t_0_yyz_xz_xy[i] - rab_y[i] * t_0_yz_xz_xy[i];

        t_y_yz_xz_xx[i] = t_0_yyz_xz_xx[i] - rab_y[i] * t_0_yz_xz_xx[i];

        t_y_yz_xy_zz[i] = t_0_yyz_xy_zz[i] - rab_y[i] * t_0_yz_xy_zz[i];

        t_y_yz_xy_yz[i] = t_0_yyz_xy_yz[i] - rab_y[i] * t_0_yz_xy_yz[i];

        t_y_yz_xy_yy[i] = t_0_yyz_xy_yy[i] - rab_y[i] * t_0_yz_xy_yy[i];

        t_y_yz_xy_xz[i] = t_0_yyz_xy_xz[i] - rab_y[i] * t_0_yz_xy_xz[i];

        t_y_yz_xy_xy[i] = t_0_yyz_xy_xy[i] - rab_y[i] * t_0_yz_xy_xy[i];

        t_y_yz_xy_xx[i] = t_0_yyz_xy_xx[i] - rab_y[i] * t_0_yz_xy_xx[i];

        t_y_yz_xx_zz[i] = t_0_yyz_xx_zz[i] - rab_y[i] * t_0_yz_xx_zz[i];

        t_y_yz_xx_yz[i] = t_0_yyz_xx_yz[i] - rab_y[i] * t_0_yz_xx_yz[i];

        t_y_yz_xx_yy[i] = t_0_yyz_xx_yy[i] - rab_y[i] * t_0_yz_xx_yy[i];

        t_y_yz_xx_xz[i] = t_0_yyz_xx_xz[i] - rab_y[i] * t_0_yz_xx_xz[i];

        t_y_yz_xx_xy[i] = t_0_yyz_xx_xy[i] - rab_y[i] * t_0_yz_xx_xy[i];

        t_y_yz_xx_xx[i] = t_0_yyz_xx_xx[i] - rab_y[i] * t_0_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yy_xx_xx, t_0_yy_xx_xy, t_0_yy_xx_xz, t_0_yy_xx_yy,\
                           t_0_yy_xx_yz, t_0_yy_xx_zz, t_0_yy_xy_xx, t_0_yy_xy_xy, t_0_yy_xy_xz,\
                           t_0_yy_xy_yy, t_0_yy_xy_yz, t_0_yy_xy_zz, t_0_yy_xz_xx, t_0_yy_xz_xy,\
                           t_0_yy_xz_xz, t_0_yy_xz_yy, t_0_yy_xz_yz, t_0_yy_xz_zz, t_0_yy_yy_xx,\
                           t_0_yy_yy_xy, t_0_yy_yy_xz, t_0_yy_yy_yy, t_0_yy_yy_yz, t_0_yy_yy_zz,\
                           t_0_yy_yz_xx, t_0_yy_yz_xy, t_0_yy_yz_xz, t_0_yy_yz_yy, t_0_yy_yz_yz,\
                           t_0_yy_yz_zz, t_0_yy_zz_xx, t_0_yy_zz_xy, t_0_yy_zz_xz, t_0_yy_zz_yy,\
                           t_0_yy_zz_yz, t_0_yy_zz_zz, t_0_yyy_xx_xx, t_0_yyy_xx_xy,\
                           t_0_yyy_xx_xz, t_0_yyy_xx_yy, t_0_yyy_xx_yz, t_0_yyy_xx_zz,\
                           t_0_yyy_xy_xx, t_0_yyy_xy_xy, t_0_yyy_xy_xz, t_0_yyy_xy_yy,\
                           t_0_yyy_xy_yz, t_0_yyy_xy_zz, t_0_yyy_xz_xx, t_0_yyy_xz_xy,\
                           t_0_yyy_xz_xz, t_0_yyy_xz_yy, t_0_yyy_xz_yz, t_0_yyy_xz_zz,\
                           t_0_yyy_yy_xx, t_0_yyy_yy_xy, t_0_yyy_yy_xz, t_0_yyy_yy_yy,\
                           t_0_yyy_yy_yz, t_0_yyy_yy_zz, t_0_yyy_yz_xx, t_0_yyy_yz_xy,\
                           t_0_yyy_yz_xz, t_0_yyy_yz_yy, t_0_yyy_yz_yz, t_0_yyy_yz_zz,\
                           t_0_yyy_zz_xx, t_0_yyy_zz_xy, t_0_yyy_zz_xz, t_0_yyy_zz_yy,\
                           t_0_yyy_zz_yz, t_0_yyy_zz_zz, t_y_yy_xx_xx, t_y_yy_xx_xy,\
                           t_y_yy_xx_xz, t_y_yy_xx_yy, t_y_yy_xx_yz, t_y_yy_xx_zz, t_y_yy_xy_xx,\
                           t_y_yy_xy_xy, t_y_yy_xy_xz, t_y_yy_xy_yy, t_y_yy_xy_yz, t_y_yy_xy_zz,\
                           t_y_yy_xz_xx, t_y_yy_xz_xy, t_y_yy_xz_xz, t_y_yy_xz_yy, t_y_yy_xz_yz,\
                           t_y_yy_xz_zz, t_y_yy_yy_xx, t_y_yy_yy_xy, t_y_yy_yy_xz, t_y_yy_yy_yy,\
                           t_y_yy_yy_yz, t_y_yy_yy_zz, t_y_yy_yz_xx, t_y_yy_yz_xy, t_y_yy_yz_xz,\
                           t_y_yy_yz_yy, t_y_yy_yz_yz, t_y_yy_yz_zz, t_y_yy_zz_xx, t_y_yy_zz_xy,\
                           t_y_yy_zz_xz, t_y_yy_zz_yy, t_y_yy_zz_yz, t_y_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yy_zz_zz[i] = t_0_yyy_zz_zz[i] - rab_y[i] * t_0_yy_zz_zz[i];

        t_y_yy_zz_yz[i] = t_0_yyy_zz_yz[i] - rab_y[i] * t_0_yy_zz_yz[i];

        t_y_yy_zz_yy[i] = t_0_yyy_zz_yy[i] - rab_y[i] * t_0_yy_zz_yy[i];

        t_y_yy_zz_xz[i] = t_0_yyy_zz_xz[i] - rab_y[i] * t_0_yy_zz_xz[i];

        t_y_yy_zz_xy[i] = t_0_yyy_zz_xy[i] - rab_y[i] * t_0_yy_zz_xy[i];

        t_y_yy_zz_xx[i] = t_0_yyy_zz_xx[i] - rab_y[i] * t_0_yy_zz_xx[i];

        t_y_yy_yz_zz[i] = t_0_yyy_yz_zz[i] - rab_y[i] * t_0_yy_yz_zz[i];

        t_y_yy_yz_yz[i] = t_0_yyy_yz_yz[i] - rab_y[i] * t_0_yy_yz_yz[i];

        t_y_yy_yz_yy[i] = t_0_yyy_yz_yy[i] - rab_y[i] * t_0_yy_yz_yy[i];

        t_y_yy_yz_xz[i] = t_0_yyy_yz_xz[i] - rab_y[i] * t_0_yy_yz_xz[i];

        t_y_yy_yz_xy[i] = t_0_yyy_yz_xy[i] - rab_y[i] * t_0_yy_yz_xy[i];

        t_y_yy_yz_xx[i] = t_0_yyy_yz_xx[i] - rab_y[i] * t_0_yy_yz_xx[i];

        t_y_yy_yy_zz[i] = t_0_yyy_yy_zz[i] - rab_y[i] * t_0_yy_yy_zz[i];

        t_y_yy_yy_yz[i] = t_0_yyy_yy_yz[i] - rab_y[i] * t_0_yy_yy_yz[i];

        t_y_yy_yy_yy[i] = t_0_yyy_yy_yy[i] - rab_y[i] * t_0_yy_yy_yy[i];

        t_y_yy_yy_xz[i] = t_0_yyy_yy_xz[i] - rab_y[i] * t_0_yy_yy_xz[i];

        t_y_yy_yy_xy[i] = t_0_yyy_yy_xy[i] - rab_y[i] * t_0_yy_yy_xy[i];

        t_y_yy_yy_xx[i] = t_0_yyy_yy_xx[i] - rab_y[i] * t_0_yy_yy_xx[i];

        t_y_yy_xz_zz[i] = t_0_yyy_xz_zz[i] - rab_y[i] * t_0_yy_xz_zz[i];

        t_y_yy_xz_yz[i] = t_0_yyy_xz_yz[i] - rab_y[i] * t_0_yy_xz_yz[i];

        t_y_yy_xz_yy[i] = t_0_yyy_xz_yy[i] - rab_y[i] * t_0_yy_xz_yy[i];

        t_y_yy_xz_xz[i] = t_0_yyy_xz_xz[i] - rab_y[i] * t_0_yy_xz_xz[i];

        t_y_yy_xz_xy[i] = t_0_yyy_xz_xy[i] - rab_y[i] * t_0_yy_xz_xy[i];

        t_y_yy_xz_xx[i] = t_0_yyy_xz_xx[i] - rab_y[i] * t_0_yy_xz_xx[i];

        t_y_yy_xy_zz[i] = t_0_yyy_xy_zz[i] - rab_y[i] * t_0_yy_xy_zz[i];

        t_y_yy_xy_yz[i] = t_0_yyy_xy_yz[i] - rab_y[i] * t_0_yy_xy_yz[i];

        t_y_yy_xy_yy[i] = t_0_yyy_xy_yy[i] - rab_y[i] * t_0_yy_xy_yy[i];

        t_y_yy_xy_xz[i] = t_0_yyy_xy_xz[i] - rab_y[i] * t_0_yy_xy_xz[i];

        t_y_yy_xy_xy[i] = t_0_yyy_xy_xy[i] - rab_y[i] * t_0_yy_xy_xy[i];

        t_y_yy_xy_xx[i] = t_0_yyy_xy_xx[i] - rab_y[i] * t_0_yy_xy_xx[i];

        t_y_yy_xx_zz[i] = t_0_yyy_xx_zz[i] - rab_y[i] * t_0_yy_xx_zz[i];

        t_y_yy_xx_yz[i] = t_0_yyy_xx_yz[i] - rab_y[i] * t_0_yy_xx_yz[i];

        t_y_yy_xx_yy[i] = t_0_yyy_xx_yy[i] - rab_y[i] * t_0_yy_xx_yy[i];

        t_y_yy_xx_xz[i] = t_0_yyy_xx_xz[i] - rab_y[i] * t_0_yy_xx_xz[i];

        t_y_yy_xx_xy[i] = t_0_yyy_xx_xy[i] - rab_y[i] * t_0_yy_xx_xy[i];

        t_y_yy_xx_xx[i] = t_0_yyy_xx_xx[i] - rab_y[i] * t_0_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz, t_0_xyz_xx_yy,\
                           t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx, t_0_xyz_xy_xy,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz, t_0_xyz_xy_zz,\
                           t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz, t_0_xyz_xz_yy,\
                           t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz,\
                           t_0_xz_xx_xx, t_0_xz_xx_xy, t_0_xz_xx_xz, t_0_xz_xx_yy, t_0_xz_xx_yz,\
                           t_0_xz_xx_zz, t_0_xz_xy_xx, t_0_xz_xy_xy, t_0_xz_xy_xz, t_0_xz_xy_yy,\
                           t_0_xz_xy_yz, t_0_xz_xy_zz, t_0_xz_xz_xx, t_0_xz_xz_xy, t_0_xz_xz_xz,\
                           t_0_xz_xz_yy, t_0_xz_xz_yz, t_0_xz_xz_zz, t_0_xz_yy_xx, t_0_xz_yy_xy,\
                           t_0_xz_yy_xz, t_0_xz_yy_yy, t_0_xz_yy_yz, t_0_xz_yy_zz, t_0_xz_yz_xx,\
                           t_0_xz_yz_xy, t_0_xz_yz_xz, t_0_xz_yz_yy, t_0_xz_yz_yz, t_0_xz_yz_zz,\
                           t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz, t_0_xz_zz_yy, t_0_xz_zz_yz,\
                           t_0_xz_zz_zz, t_y_xz_xx_xx, t_y_xz_xx_xy, t_y_xz_xx_xz, t_y_xz_xx_yy,\
                           t_y_xz_xx_yz, t_y_xz_xx_zz, t_y_xz_xy_xx, t_y_xz_xy_xy, t_y_xz_xy_xz,\
                           t_y_xz_xy_yy, t_y_xz_xy_yz, t_y_xz_xy_zz, t_y_xz_xz_xx, t_y_xz_xz_xy,\
                           t_y_xz_xz_xz, t_y_xz_xz_yy, t_y_xz_xz_yz, t_y_xz_xz_zz, t_y_xz_yy_xx,\
                           t_y_xz_yy_xy, t_y_xz_yy_xz, t_y_xz_yy_yy, t_y_xz_yy_yz, t_y_xz_yy_zz,\
                           t_y_xz_yz_xx, t_y_xz_yz_xy, t_y_xz_yz_xz, t_y_xz_yz_yy, t_y_xz_yz_yz,\
                           t_y_xz_yz_zz, t_y_xz_zz_xx, t_y_xz_zz_xy, t_y_xz_zz_xz, t_y_xz_zz_yy,\
                           t_y_xz_zz_yz, t_y_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xz_zz_zz[i] = t_0_xyz_zz_zz[i] - rab_y[i] * t_0_xz_zz_zz[i];

        t_y_xz_zz_yz[i] = t_0_xyz_zz_yz[i] - rab_y[i] * t_0_xz_zz_yz[i];

        t_y_xz_zz_yy[i] = t_0_xyz_zz_yy[i] - rab_y[i] * t_0_xz_zz_yy[i];

        t_y_xz_zz_xz[i] = t_0_xyz_zz_xz[i] - rab_y[i] * t_0_xz_zz_xz[i];

        t_y_xz_zz_xy[i] = t_0_xyz_zz_xy[i] - rab_y[i] * t_0_xz_zz_xy[i];

        t_y_xz_zz_xx[i] = t_0_xyz_zz_xx[i] - rab_y[i] * t_0_xz_zz_xx[i];

        t_y_xz_yz_zz[i] = t_0_xyz_yz_zz[i] - rab_y[i] * t_0_xz_yz_zz[i];

        t_y_xz_yz_yz[i] = t_0_xyz_yz_yz[i] - rab_y[i] * t_0_xz_yz_yz[i];

        t_y_xz_yz_yy[i] = t_0_xyz_yz_yy[i] - rab_y[i] * t_0_xz_yz_yy[i];

        t_y_xz_yz_xz[i] = t_0_xyz_yz_xz[i] - rab_y[i] * t_0_xz_yz_xz[i];

        t_y_xz_yz_xy[i] = t_0_xyz_yz_xy[i] - rab_y[i] * t_0_xz_yz_xy[i];

        t_y_xz_yz_xx[i] = t_0_xyz_yz_xx[i] - rab_y[i] * t_0_xz_yz_xx[i];

        t_y_xz_yy_zz[i] = t_0_xyz_yy_zz[i] - rab_y[i] * t_0_xz_yy_zz[i];

        t_y_xz_yy_yz[i] = t_0_xyz_yy_yz[i] - rab_y[i] * t_0_xz_yy_yz[i];

        t_y_xz_yy_yy[i] = t_0_xyz_yy_yy[i] - rab_y[i] * t_0_xz_yy_yy[i];

        t_y_xz_yy_xz[i] = t_0_xyz_yy_xz[i] - rab_y[i] * t_0_xz_yy_xz[i];

        t_y_xz_yy_xy[i] = t_0_xyz_yy_xy[i] - rab_y[i] * t_0_xz_yy_xy[i];

        t_y_xz_yy_xx[i] = t_0_xyz_yy_xx[i] - rab_y[i] * t_0_xz_yy_xx[i];

        t_y_xz_xz_zz[i] = t_0_xyz_xz_zz[i] - rab_y[i] * t_0_xz_xz_zz[i];

        t_y_xz_xz_yz[i] = t_0_xyz_xz_yz[i] - rab_y[i] * t_0_xz_xz_yz[i];

        t_y_xz_xz_yy[i] = t_0_xyz_xz_yy[i] - rab_y[i] * t_0_xz_xz_yy[i];

        t_y_xz_xz_xz[i] = t_0_xyz_xz_xz[i] - rab_y[i] * t_0_xz_xz_xz[i];

        t_y_xz_xz_xy[i] = t_0_xyz_xz_xy[i] - rab_y[i] * t_0_xz_xz_xy[i];

        t_y_xz_xz_xx[i] = t_0_xyz_xz_xx[i] - rab_y[i] * t_0_xz_xz_xx[i];

        t_y_xz_xy_zz[i] = t_0_xyz_xy_zz[i] - rab_y[i] * t_0_xz_xy_zz[i];

        t_y_xz_xy_yz[i] = t_0_xyz_xy_yz[i] - rab_y[i] * t_0_xz_xy_yz[i];

        t_y_xz_xy_yy[i] = t_0_xyz_xy_yy[i] - rab_y[i] * t_0_xz_xy_yy[i];

        t_y_xz_xy_xz[i] = t_0_xyz_xy_xz[i] - rab_y[i] * t_0_xz_xy_xz[i];

        t_y_xz_xy_xy[i] = t_0_xyz_xy_xy[i] - rab_y[i] * t_0_xz_xy_xy[i];

        t_y_xz_xy_xx[i] = t_0_xyz_xy_xx[i] - rab_y[i] * t_0_xz_xy_xx[i];

        t_y_xz_xx_zz[i] = t_0_xyz_xx_zz[i] - rab_y[i] * t_0_xz_xx_zz[i];

        t_y_xz_xx_yz[i] = t_0_xyz_xx_yz[i] - rab_y[i] * t_0_xz_xx_yz[i];

        t_y_xz_xx_yy[i] = t_0_xyz_xx_yy[i] - rab_y[i] * t_0_xz_xx_yy[i];

        t_y_xz_xx_xz[i] = t_0_xyz_xx_xz[i] - rab_y[i] * t_0_xz_xx_xz[i];

        t_y_xz_xx_xy[i] = t_0_xyz_xx_xy[i] - rab_y[i] * t_0_xz_xx_xy[i];

        t_y_xz_xx_xx[i] = t_0_xyz_xx_xx[i] - rab_y[i] * t_0_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xy_xx_xx, t_0_xy_xx_xy, t_0_xy_xx_xz, t_0_xy_xx_yy,\
                           t_0_xy_xx_yz, t_0_xy_xx_zz, t_0_xy_xy_xx, t_0_xy_xy_xy, t_0_xy_xy_xz,\
                           t_0_xy_xy_yy, t_0_xy_xy_yz, t_0_xy_xy_zz, t_0_xy_xz_xx, t_0_xy_xz_xy,\
                           t_0_xy_xz_xz, t_0_xy_xz_yy, t_0_xy_xz_yz, t_0_xy_xz_zz, t_0_xy_yy_xx,\
                           t_0_xy_yy_xy, t_0_xy_yy_xz, t_0_xy_yy_yy, t_0_xy_yy_yz, t_0_xy_yy_zz,\
                           t_0_xy_yz_xx, t_0_xy_yz_xy, t_0_xy_yz_xz, t_0_xy_yz_yy, t_0_xy_yz_yz,\
                           t_0_xy_yz_zz, t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz, t_0_xy_zz_yy,\
                           t_0_xy_zz_yz, t_0_xy_zz_zz, t_0_xyy_xx_xx, t_0_xyy_xx_xy,\
                           t_0_xyy_xx_xz, t_0_xyy_xx_yy, t_0_xyy_xx_yz, t_0_xyy_xx_zz,\
                           t_0_xyy_xy_xx, t_0_xyy_xy_xy, t_0_xyy_xy_xz, t_0_xyy_xy_yy,\
                           t_0_xyy_xy_yz, t_0_xyy_xy_zz, t_0_xyy_xz_xx, t_0_xyy_xz_xy,\
                           t_0_xyy_xz_xz, t_0_xyy_xz_yy, t_0_xyy_xz_yz, t_0_xyy_xz_zz,\
                           t_0_xyy_yy_xx, t_0_xyy_yy_xy, t_0_xyy_yy_xz, t_0_xyy_yy_yy,\
                           t_0_xyy_yy_yz, t_0_xyy_yy_zz, t_0_xyy_yz_xx, t_0_xyy_yz_xy,\
                           t_0_xyy_yz_xz, t_0_xyy_yz_yy, t_0_xyy_yz_yz, t_0_xyy_yz_zz,\
                           t_0_xyy_zz_xx, t_0_xyy_zz_xy, t_0_xyy_zz_xz, t_0_xyy_zz_yy,\
                           t_0_xyy_zz_yz, t_0_xyy_zz_zz, t_y_xy_xx_xx, t_y_xy_xx_xy,\
                           t_y_xy_xx_xz, t_y_xy_xx_yy, t_y_xy_xx_yz, t_y_xy_xx_zz, t_y_xy_xy_xx,\
                           t_y_xy_xy_xy, t_y_xy_xy_xz, t_y_xy_xy_yy, t_y_xy_xy_yz, t_y_xy_xy_zz,\
                           t_y_xy_xz_xx, t_y_xy_xz_xy, t_y_xy_xz_xz, t_y_xy_xz_yy, t_y_xy_xz_yz,\
                           t_y_xy_xz_zz, t_y_xy_yy_xx, t_y_xy_yy_xy, t_y_xy_yy_xz, t_y_xy_yy_yy,\
                           t_y_xy_yy_yz, t_y_xy_yy_zz, t_y_xy_yz_xx, t_y_xy_yz_xy, t_y_xy_yz_xz,\
                           t_y_xy_yz_yy, t_y_xy_yz_yz, t_y_xy_yz_zz, t_y_xy_zz_xx, t_y_xy_zz_xy,\
                           t_y_xy_zz_xz, t_y_xy_zz_yy, t_y_xy_zz_yz, t_y_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xy_zz_zz[i] = t_0_xyy_zz_zz[i] - rab_y[i] * t_0_xy_zz_zz[i];

        t_y_xy_zz_yz[i] = t_0_xyy_zz_yz[i] - rab_y[i] * t_0_xy_zz_yz[i];

        t_y_xy_zz_yy[i] = t_0_xyy_zz_yy[i] - rab_y[i] * t_0_xy_zz_yy[i];

        t_y_xy_zz_xz[i] = t_0_xyy_zz_xz[i] - rab_y[i] * t_0_xy_zz_xz[i];

        t_y_xy_zz_xy[i] = t_0_xyy_zz_xy[i] - rab_y[i] * t_0_xy_zz_xy[i];

        t_y_xy_zz_xx[i] = t_0_xyy_zz_xx[i] - rab_y[i] * t_0_xy_zz_xx[i];

        t_y_xy_yz_zz[i] = t_0_xyy_yz_zz[i] - rab_y[i] * t_0_xy_yz_zz[i];

        t_y_xy_yz_yz[i] = t_0_xyy_yz_yz[i] - rab_y[i] * t_0_xy_yz_yz[i];

        t_y_xy_yz_yy[i] = t_0_xyy_yz_yy[i] - rab_y[i] * t_0_xy_yz_yy[i];

        t_y_xy_yz_xz[i] = t_0_xyy_yz_xz[i] - rab_y[i] * t_0_xy_yz_xz[i];

        t_y_xy_yz_xy[i] = t_0_xyy_yz_xy[i] - rab_y[i] * t_0_xy_yz_xy[i];

        t_y_xy_yz_xx[i] = t_0_xyy_yz_xx[i] - rab_y[i] * t_0_xy_yz_xx[i];

        t_y_xy_yy_zz[i] = t_0_xyy_yy_zz[i] - rab_y[i] * t_0_xy_yy_zz[i];

        t_y_xy_yy_yz[i] = t_0_xyy_yy_yz[i] - rab_y[i] * t_0_xy_yy_yz[i];

        t_y_xy_yy_yy[i] = t_0_xyy_yy_yy[i] - rab_y[i] * t_0_xy_yy_yy[i];

        t_y_xy_yy_xz[i] = t_0_xyy_yy_xz[i] - rab_y[i] * t_0_xy_yy_xz[i];

        t_y_xy_yy_xy[i] = t_0_xyy_yy_xy[i] - rab_y[i] * t_0_xy_yy_xy[i];

        t_y_xy_yy_xx[i] = t_0_xyy_yy_xx[i] - rab_y[i] * t_0_xy_yy_xx[i];

        t_y_xy_xz_zz[i] = t_0_xyy_xz_zz[i] - rab_y[i] * t_0_xy_xz_zz[i];

        t_y_xy_xz_yz[i] = t_0_xyy_xz_yz[i] - rab_y[i] * t_0_xy_xz_yz[i];

        t_y_xy_xz_yy[i] = t_0_xyy_xz_yy[i] - rab_y[i] * t_0_xy_xz_yy[i];

        t_y_xy_xz_xz[i] = t_0_xyy_xz_xz[i] - rab_y[i] * t_0_xy_xz_xz[i];

        t_y_xy_xz_xy[i] = t_0_xyy_xz_xy[i] - rab_y[i] * t_0_xy_xz_xy[i];

        t_y_xy_xz_xx[i] = t_0_xyy_xz_xx[i] - rab_y[i] * t_0_xy_xz_xx[i];

        t_y_xy_xy_zz[i] = t_0_xyy_xy_zz[i] - rab_y[i] * t_0_xy_xy_zz[i];

        t_y_xy_xy_yz[i] = t_0_xyy_xy_yz[i] - rab_y[i] * t_0_xy_xy_yz[i];

        t_y_xy_xy_yy[i] = t_0_xyy_xy_yy[i] - rab_y[i] * t_0_xy_xy_yy[i];

        t_y_xy_xy_xz[i] = t_0_xyy_xy_xz[i] - rab_y[i] * t_0_xy_xy_xz[i];

        t_y_xy_xy_xy[i] = t_0_xyy_xy_xy[i] - rab_y[i] * t_0_xy_xy_xy[i];

        t_y_xy_xy_xx[i] = t_0_xyy_xy_xx[i] - rab_y[i] * t_0_xy_xy_xx[i];

        t_y_xy_xx_zz[i] = t_0_xyy_xx_zz[i] - rab_y[i] * t_0_xy_xx_zz[i];

        t_y_xy_xx_yz[i] = t_0_xyy_xx_yz[i] - rab_y[i] * t_0_xy_xx_yz[i];

        t_y_xy_xx_yy[i] = t_0_xyy_xx_yy[i] - rab_y[i] * t_0_xy_xx_yy[i];

        t_y_xy_xx_xz[i] = t_0_xyy_xx_xz[i] - rab_y[i] * t_0_xy_xx_xz[i];

        t_y_xy_xx_xy[i] = t_0_xyy_xx_xy[i] - rab_y[i] * t_0_xy_xx_xy[i];

        t_y_xy_xx_xx[i] = t_0_xyy_xx_xx[i] - rab_y[i] * t_0_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xx_xx_xx, t_0_xx_xx_xy, t_0_xx_xx_xz, t_0_xx_xx_yy,\
                           t_0_xx_xx_yz, t_0_xx_xx_zz, t_0_xx_xy_xx, t_0_xx_xy_xy, t_0_xx_xy_xz,\
                           t_0_xx_xy_yy, t_0_xx_xy_yz, t_0_xx_xy_zz, t_0_xx_xz_xx, t_0_xx_xz_xy,\
                           t_0_xx_xz_xz, t_0_xx_xz_yy, t_0_xx_xz_yz, t_0_xx_xz_zz, t_0_xx_yy_xx,\
                           t_0_xx_yy_xy, t_0_xx_yy_xz, t_0_xx_yy_yy, t_0_xx_yy_yz, t_0_xx_yy_zz,\
                           t_0_xx_yz_xx, t_0_xx_yz_xy, t_0_xx_yz_xz, t_0_xx_yz_yy, t_0_xx_yz_yz,\
                           t_0_xx_yz_zz, t_0_xx_zz_xx, t_0_xx_zz_xy, t_0_xx_zz_xz, t_0_xx_zz_yy,\
                           t_0_xx_zz_yz, t_0_xx_zz_zz, t_0_xxy_xx_xx, t_0_xxy_xx_xy,\
                           t_0_xxy_xx_xz, t_0_xxy_xx_yy, t_0_xxy_xx_yz, t_0_xxy_xx_zz,\
                           t_0_xxy_xy_xx, t_0_xxy_xy_xy, t_0_xxy_xy_xz, t_0_xxy_xy_yy,\
                           t_0_xxy_xy_yz, t_0_xxy_xy_zz, t_0_xxy_xz_xx, t_0_xxy_xz_xy,\
                           t_0_xxy_xz_xz, t_0_xxy_xz_yy, t_0_xxy_xz_yz, t_0_xxy_xz_zz,\
                           t_0_xxy_yy_xx, t_0_xxy_yy_xy, t_0_xxy_yy_xz, t_0_xxy_yy_yy,\
                           t_0_xxy_yy_yz, t_0_xxy_yy_zz, t_0_xxy_yz_xx, t_0_xxy_yz_xy,\
                           t_0_xxy_yz_xz, t_0_xxy_yz_yy, t_0_xxy_yz_yz, t_0_xxy_yz_zz,\
                           t_0_xxy_zz_xx, t_0_xxy_zz_xy, t_0_xxy_zz_xz, t_0_xxy_zz_yy,\
                           t_0_xxy_zz_yz, t_0_xxy_zz_zz, t_y_xx_xx_xx, t_y_xx_xx_xy,\
                           t_y_xx_xx_xz, t_y_xx_xx_yy, t_y_xx_xx_yz, t_y_xx_xx_zz, t_y_xx_xy_xx,\
                           t_y_xx_xy_xy, t_y_xx_xy_xz, t_y_xx_xy_yy, t_y_xx_xy_yz, t_y_xx_xy_zz,\
                           t_y_xx_xz_xx, t_y_xx_xz_xy, t_y_xx_xz_xz, t_y_xx_xz_yy, t_y_xx_xz_yz,\
                           t_y_xx_xz_zz, t_y_xx_yy_xx, t_y_xx_yy_xy, t_y_xx_yy_xz, t_y_xx_yy_yy,\
                           t_y_xx_yy_yz, t_y_xx_yy_zz, t_y_xx_yz_xx, t_y_xx_yz_xy, t_y_xx_yz_xz,\
                           t_y_xx_yz_yy, t_y_xx_yz_yz, t_y_xx_yz_zz, t_y_xx_zz_xx, t_y_xx_zz_xy,\
                           t_y_xx_zz_xz, t_y_xx_zz_yy, t_y_xx_zz_yz, t_y_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xx_zz_zz[i] = t_0_xxy_zz_zz[i] - rab_y[i] * t_0_xx_zz_zz[i];

        t_y_xx_zz_yz[i] = t_0_xxy_zz_yz[i] - rab_y[i] * t_0_xx_zz_yz[i];

        t_y_xx_zz_yy[i] = t_0_xxy_zz_yy[i] - rab_y[i] * t_0_xx_zz_yy[i];

        t_y_xx_zz_xz[i] = t_0_xxy_zz_xz[i] - rab_y[i] * t_0_xx_zz_xz[i];

        t_y_xx_zz_xy[i] = t_0_xxy_zz_xy[i] - rab_y[i] * t_0_xx_zz_xy[i];

        t_y_xx_zz_xx[i] = t_0_xxy_zz_xx[i] - rab_y[i] * t_0_xx_zz_xx[i];

        t_y_xx_yz_zz[i] = t_0_xxy_yz_zz[i] - rab_y[i] * t_0_xx_yz_zz[i];

        t_y_xx_yz_yz[i] = t_0_xxy_yz_yz[i] - rab_y[i] * t_0_xx_yz_yz[i];

        t_y_xx_yz_yy[i] = t_0_xxy_yz_yy[i] - rab_y[i] * t_0_xx_yz_yy[i];

        t_y_xx_yz_xz[i] = t_0_xxy_yz_xz[i] - rab_y[i] * t_0_xx_yz_xz[i];

        t_y_xx_yz_xy[i] = t_0_xxy_yz_xy[i] - rab_y[i] * t_0_xx_yz_xy[i];

        t_y_xx_yz_xx[i] = t_0_xxy_yz_xx[i] - rab_y[i] * t_0_xx_yz_xx[i];

        t_y_xx_yy_zz[i] = t_0_xxy_yy_zz[i] - rab_y[i] * t_0_xx_yy_zz[i];

        t_y_xx_yy_yz[i] = t_0_xxy_yy_yz[i] - rab_y[i] * t_0_xx_yy_yz[i];

        t_y_xx_yy_yy[i] = t_0_xxy_yy_yy[i] - rab_y[i] * t_0_xx_yy_yy[i];

        t_y_xx_yy_xz[i] = t_0_xxy_yy_xz[i] - rab_y[i] * t_0_xx_yy_xz[i];

        t_y_xx_yy_xy[i] = t_0_xxy_yy_xy[i] - rab_y[i] * t_0_xx_yy_xy[i];

        t_y_xx_yy_xx[i] = t_0_xxy_yy_xx[i] - rab_y[i] * t_0_xx_yy_xx[i];

        t_y_xx_xz_zz[i] = t_0_xxy_xz_zz[i] - rab_y[i] * t_0_xx_xz_zz[i];

        t_y_xx_xz_yz[i] = t_0_xxy_xz_yz[i] - rab_y[i] * t_0_xx_xz_yz[i];

        t_y_xx_xz_yy[i] = t_0_xxy_xz_yy[i] - rab_y[i] * t_0_xx_xz_yy[i];

        t_y_xx_xz_xz[i] = t_0_xxy_xz_xz[i] - rab_y[i] * t_0_xx_xz_xz[i];

        t_y_xx_xz_xy[i] = t_0_xxy_xz_xy[i] - rab_y[i] * t_0_xx_xz_xy[i];

        t_y_xx_xz_xx[i] = t_0_xxy_xz_xx[i] - rab_y[i] * t_0_xx_xz_xx[i];

        t_y_xx_xy_zz[i] = t_0_xxy_xy_zz[i] - rab_y[i] * t_0_xx_xy_zz[i];

        t_y_xx_xy_yz[i] = t_0_xxy_xy_yz[i] - rab_y[i] * t_0_xx_xy_yz[i];

        t_y_xx_xy_yy[i] = t_0_xxy_xy_yy[i] - rab_y[i] * t_0_xx_xy_yy[i];

        t_y_xx_xy_xz[i] = t_0_xxy_xy_xz[i] - rab_y[i] * t_0_xx_xy_xz[i];

        t_y_xx_xy_xy[i] = t_0_xxy_xy_xy[i] - rab_y[i] * t_0_xx_xy_xy[i];

        t_y_xx_xy_xx[i] = t_0_xxy_xy_xx[i] - rab_y[i] * t_0_xx_xy_xx[i];

        t_y_xx_xx_zz[i] = t_0_xxy_xx_zz[i] - rab_y[i] * t_0_xx_xx_zz[i];

        t_y_xx_xx_yz[i] = t_0_xxy_xx_yz[i] - rab_y[i] * t_0_xx_xx_yz[i];

        t_y_xx_xx_yy[i] = t_0_xxy_xx_yy[i] - rab_y[i] * t_0_xx_xx_yy[i];

        t_y_xx_xx_xz[i] = t_0_xxy_xx_xz[i] - rab_y[i] * t_0_xx_xx_xz[i];

        t_y_xx_xx_xy[i] = t_0_xxy_xx_xy[i] - rab_y[i] * t_0_xx_xx_xy[i];

        t_y_xx_xx_xx[i] = t_0_xxy_xx_xx[i] - rab_y[i] * t_0_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xzz_xx_xx, t_0_xzz_xx_xy, t_0_xzz_xx_xz, t_0_xzz_xx_yy,\
                           t_0_xzz_xx_yz, t_0_xzz_xx_zz, t_0_xzz_xy_xx, t_0_xzz_xy_xy,\
                           t_0_xzz_xy_xz, t_0_xzz_xy_yy, t_0_xzz_xy_yz, t_0_xzz_xy_zz,\
                           t_0_xzz_xz_xx, t_0_xzz_xz_xy, t_0_xzz_xz_xz, t_0_xzz_xz_yy,\
                           t_0_xzz_xz_yz, t_0_xzz_xz_zz, t_0_xzz_yy_xx, t_0_xzz_yy_xy,\
                           t_0_xzz_yy_xz, t_0_xzz_yy_yy, t_0_xzz_yy_yz, t_0_xzz_yy_zz,\
                           t_0_xzz_yz_xx, t_0_xzz_yz_xy, t_0_xzz_yz_xz, t_0_xzz_yz_yy,\
                           t_0_xzz_yz_yz, t_0_xzz_yz_zz, t_0_xzz_zz_xx, t_0_xzz_zz_xy,\
                           t_0_xzz_zz_xz, t_0_xzz_zz_yy, t_0_xzz_zz_yz, t_0_xzz_zz_zz,\
                           t_0_zz_xx_xx, t_0_zz_xx_xy, t_0_zz_xx_xz, t_0_zz_xx_yy, t_0_zz_xx_yz,\
                           t_0_zz_xx_zz, t_0_zz_xy_xx, t_0_zz_xy_xy, t_0_zz_xy_xz, t_0_zz_xy_yy,\
                           t_0_zz_xy_yz, t_0_zz_xy_zz, t_0_zz_xz_xx, t_0_zz_xz_xy, t_0_zz_xz_xz,\
                           t_0_zz_xz_yy, t_0_zz_xz_yz, t_0_zz_xz_zz, t_0_zz_yy_xx, t_0_zz_yy_xy,\
                           t_0_zz_yy_xz, t_0_zz_yy_yy, t_0_zz_yy_yz, t_0_zz_yy_zz, t_0_zz_yz_xx,\
                           t_0_zz_yz_xy, t_0_zz_yz_xz, t_0_zz_yz_yy, t_0_zz_yz_yz, t_0_zz_yz_zz,\
                           t_0_zz_zz_xx, t_0_zz_zz_xy, t_0_zz_zz_xz, t_0_zz_zz_yy, t_0_zz_zz_yz,\
                           t_0_zz_zz_zz, t_x_zz_xx_xx, t_x_zz_xx_xy, t_x_zz_xx_xz, t_x_zz_xx_yy,\
                           t_x_zz_xx_yz, t_x_zz_xx_zz, t_x_zz_xy_xx, t_x_zz_xy_xy, t_x_zz_xy_xz,\
                           t_x_zz_xy_yy, t_x_zz_xy_yz, t_x_zz_xy_zz, t_x_zz_xz_xx, t_x_zz_xz_xy,\
                           t_x_zz_xz_xz, t_x_zz_xz_yy, t_x_zz_xz_yz, t_x_zz_xz_zz, t_x_zz_yy_xx,\
                           t_x_zz_yy_xy, t_x_zz_yy_xz, t_x_zz_yy_yy, t_x_zz_yy_yz, t_x_zz_yy_zz,\
                           t_x_zz_yz_xx, t_x_zz_yz_xy, t_x_zz_yz_xz, t_x_zz_yz_yy, t_x_zz_yz_yz,\
                           t_x_zz_yz_zz, t_x_zz_zz_xx, t_x_zz_zz_xy, t_x_zz_zz_xz, t_x_zz_zz_yy,\
                           t_x_zz_zz_yz, t_x_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zz_zz_zz[i] = t_0_xzz_zz_zz[i] - rab_x[i] * t_0_zz_zz_zz[i];

        t_x_zz_zz_yz[i] = t_0_xzz_zz_yz[i] - rab_x[i] * t_0_zz_zz_yz[i];

        t_x_zz_zz_yy[i] = t_0_xzz_zz_yy[i] - rab_x[i] * t_0_zz_zz_yy[i];

        t_x_zz_zz_xz[i] = t_0_xzz_zz_xz[i] - rab_x[i] * t_0_zz_zz_xz[i];

        t_x_zz_zz_xy[i] = t_0_xzz_zz_xy[i] - rab_x[i] * t_0_zz_zz_xy[i];

        t_x_zz_zz_xx[i] = t_0_xzz_zz_xx[i] - rab_x[i] * t_0_zz_zz_xx[i];

        t_x_zz_yz_zz[i] = t_0_xzz_yz_zz[i] - rab_x[i] * t_0_zz_yz_zz[i];

        t_x_zz_yz_yz[i] = t_0_xzz_yz_yz[i] - rab_x[i] * t_0_zz_yz_yz[i];

        t_x_zz_yz_yy[i] = t_0_xzz_yz_yy[i] - rab_x[i] * t_0_zz_yz_yy[i];

        t_x_zz_yz_xz[i] = t_0_xzz_yz_xz[i] - rab_x[i] * t_0_zz_yz_xz[i];

        t_x_zz_yz_xy[i] = t_0_xzz_yz_xy[i] - rab_x[i] * t_0_zz_yz_xy[i];

        t_x_zz_yz_xx[i] = t_0_xzz_yz_xx[i] - rab_x[i] * t_0_zz_yz_xx[i];

        t_x_zz_yy_zz[i] = t_0_xzz_yy_zz[i] - rab_x[i] * t_0_zz_yy_zz[i];

        t_x_zz_yy_yz[i] = t_0_xzz_yy_yz[i] - rab_x[i] * t_0_zz_yy_yz[i];

        t_x_zz_yy_yy[i] = t_0_xzz_yy_yy[i] - rab_x[i] * t_0_zz_yy_yy[i];

        t_x_zz_yy_xz[i] = t_0_xzz_yy_xz[i] - rab_x[i] * t_0_zz_yy_xz[i];

        t_x_zz_yy_xy[i] = t_0_xzz_yy_xy[i] - rab_x[i] * t_0_zz_yy_xy[i];

        t_x_zz_yy_xx[i] = t_0_xzz_yy_xx[i] - rab_x[i] * t_0_zz_yy_xx[i];

        t_x_zz_xz_zz[i] = t_0_xzz_xz_zz[i] - rab_x[i] * t_0_zz_xz_zz[i];

        t_x_zz_xz_yz[i] = t_0_xzz_xz_yz[i] - rab_x[i] * t_0_zz_xz_yz[i];

        t_x_zz_xz_yy[i] = t_0_xzz_xz_yy[i] - rab_x[i] * t_0_zz_xz_yy[i];

        t_x_zz_xz_xz[i] = t_0_xzz_xz_xz[i] - rab_x[i] * t_0_zz_xz_xz[i];

        t_x_zz_xz_xy[i] = t_0_xzz_xz_xy[i] - rab_x[i] * t_0_zz_xz_xy[i];

        t_x_zz_xz_xx[i] = t_0_xzz_xz_xx[i] - rab_x[i] * t_0_zz_xz_xx[i];

        t_x_zz_xy_zz[i] = t_0_xzz_xy_zz[i] - rab_x[i] * t_0_zz_xy_zz[i];

        t_x_zz_xy_yz[i] = t_0_xzz_xy_yz[i] - rab_x[i] * t_0_zz_xy_yz[i];

        t_x_zz_xy_yy[i] = t_0_xzz_xy_yy[i] - rab_x[i] * t_0_zz_xy_yy[i];

        t_x_zz_xy_xz[i] = t_0_xzz_xy_xz[i] - rab_x[i] * t_0_zz_xy_xz[i];

        t_x_zz_xy_xy[i] = t_0_xzz_xy_xy[i] - rab_x[i] * t_0_zz_xy_xy[i];

        t_x_zz_xy_xx[i] = t_0_xzz_xy_xx[i] - rab_x[i] * t_0_zz_xy_xx[i];

        t_x_zz_xx_zz[i] = t_0_xzz_xx_zz[i] - rab_x[i] * t_0_zz_xx_zz[i];

        t_x_zz_xx_yz[i] = t_0_xzz_xx_yz[i] - rab_x[i] * t_0_zz_xx_yz[i];

        t_x_zz_xx_yy[i] = t_0_xzz_xx_yy[i] - rab_x[i] * t_0_zz_xx_yy[i];

        t_x_zz_xx_xz[i] = t_0_xzz_xx_xz[i] - rab_x[i] * t_0_zz_xx_xz[i];

        t_x_zz_xx_xy[i] = t_0_xzz_xx_xy[i] - rab_x[i] * t_0_zz_xx_xy[i];

        t_x_zz_xx_xx[i] = t_0_xzz_xx_xx[i] - rab_x[i] * t_0_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyz_xx_xx, t_0_xyz_xx_xy, t_0_xyz_xx_xz, t_0_xyz_xx_yy,\
                           t_0_xyz_xx_yz, t_0_xyz_xx_zz, t_0_xyz_xy_xx, t_0_xyz_xy_xy,\
                           t_0_xyz_xy_xz, t_0_xyz_xy_yy, t_0_xyz_xy_yz, t_0_xyz_xy_zz,\
                           t_0_xyz_xz_xx, t_0_xyz_xz_xy, t_0_xyz_xz_xz, t_0_xyz_xz_yy,\
                           t_0_xyz_xz_yz, t_0_xyz_xz_zz, t_0_xyz_yy_xx, t_0_xyz_yy_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yy_yy, t_0_xyz_yy_yz, t_0_xyz_yy_zz,\
                           t_0_xyz_yz_xx, t_0_xyz_yz_xy, t_0_xyz_yz_xz, t_0_xyz_yz_yy,\
                           t_0_xyz_yz_yz, t_0_xyz_yz_zz, t_0_xyz_zz_xx, t_0_xyz_zz_xy,\
                           t_0_xyz_zz_xz, t_0_xyz_zz_yy, t_0_xyz_zz_yz, t_0_xyz_zz_zz,\
                           t_0_yz_xx_xx, t_0_yz_xx_xy, t_0_yz_xx_xz, t_0_yz_xx_yy, t_0_yz_xx_yz,\
                           t_0_yz_xx_zz, t_0_yz_xy_xx, t_0_yz_xy_xy, t_0_yz_xy_xz, t_0_yz_xy_yy,\
                           t_0_yz_xy_yz, t_0_yz_xy_zz, t_0_yz_xz_xx, t_0_yz_xz_xy, t_0_yz_xz_xz,\
                           t_0_yz_xz_yy, t_0_yz_xz_yz, t_0_yz_xz_zz, t_0_yz_yy_xx, t_0_yz_yy_xy,\
                           t_0_yz_yy_xz, t_0_yz_yy_yy, t_0_yz_yy_yz, t_0_yz_yy_zz, t_0_yz_yz_xx,\
                           t_0_yz_yz_xy, t_0_yz_yz_xz, t_0_yz_yz_yy, t_0_yz_yz_yz, t_0_yz_yz_zz,\
                           t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz, t_0_yz_zz_yy, t_0_yz_zz_yz,\
                           t_0_yz_zz_zz, t_x_yz_xx_xx, t_x_yz_xx_xy, t_x_yz_xx_xz, t_x_yz_xx_yy,\
                           t_x_yz_xx_yz, t_x_yz_xx_zz, t_x_yz_xy_xx, t_x_yz_xy_xy, t_x_yz_xy_xz,\
                           t_x_yz_xy_yy, t_x_yz_xy_yz, t_x_yz_xy_zz, t_x_yz_xz_xx, t_x_yz_xz_xy,\
                           t_x_yz_xz_xz, t_x_yz_xz_yy, t_x_yz_xz_yz, t_x_yz_xz_zz, t_x_yz_yy_xx,\
                           t_x_yz_yy_xy, t_x_yz_yy_xz, t_x_yz_yy_yy, t_x_yz_yy_yz, t_x_yz_yy_zz,\
                           t_x_yz_yz_xx, t_x_yz_yz_xy, t_x_yz_yz_xz, t_x_yz_yz_yy, t_x_yz_yz_yz,\
                           t_x_yz_yz_zz, t_x_yz_zz_xx, t_x_yz_zz_xy, t_x_yz_zz_xz, t_x_yz_zz_yy,\
                           t_x_yz_zz_yz, t_x_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yz_zz_zz[i] = t_0_xyz_zz_zz[i] - rab_x[i] * t_0_yz_zz_zz[i];

        t_x_yz_zz_yz[i] = t_0_xyz_zz_yz[i] - rab_x[i] * t_0_yz_zz_yz[i];

        t_x_yz_zz_yy[i] = t_0_xyz_zz_yy[i] - rab_x[i] * t_0_yz_zz_yy[i];

        t_x_yz_zz_xz[i] = t_0_xyz_zz_xz[i] - rab_x[i] * t_0_yz_zz_xz[i];

        t_x_yz_zz_xy[i] = t_0_xyz_zz_xy[i] - rab_x[i] * t_0_yz_zz_xy[i];

        t_x_yz_zz_xx[i] = t_0_xyz_zz_xx[i] - rab_x[i] * t_0_yz_zz_xx[i];

        t_x_yz_yz_zz[i] = t_0_xyz_yz_zz[i] - rab_x[i] * t_0_yz_yz_zz[i];

        t_x_yz_yz_yz[i] = t_0_xyz_yz_yz[i] - rab_x[i] * t_0_yz_yz_yz[i];

        t_x_yz_yz_yy[i] = t_0_xyz_yz_yy[i] - rab_x[i] * t_0_yz_yz_yy[i];

        t_x_yz_yz_xz[i] = t_0_xyz_yz_xz[i] - rab_x[i] * t_0_yz_yz_xz[i];

        t_x_yz_yz_xy[i] = t_0_xyz_yz_xy[i] - rab_x[i] * t_0_yz_yz_xy[i];

        t_x_yz_yz_xx[i] = t_0_xyz_yz_xx[i] - rab_x[i] * t_0_yz_yz_xx[i];

        t_x_yz_yy_zz[i] = t_0_xyz_yy_zz[i] - rab_x[i] * t_0_yz_yy_zz[i];

        t_x_yz_yy_yz[i] = t_0_xyz_yy_yz[i] - rab_x[i] * t_0_yz_yy_yz[i];

        t_x_yz_yy_yy[i] = t_0_xyz_yy_yy[i] - rab_x[i] * t_0_yz_yy_yy[i];

        t_x_yz_yy_xz[i] = t_0_xyz_yy_xz[i] - rab_x[i] * t_0_yz_yy_xz[i];

        t_x_yz_yy_xy[i] = t_0_xyz_yy_xy[i] - rab_x[i] * t_0_yz_yy_xy[i];

        t_x_yz_yy_xx[i] = t_0_xyz_yy_xx[i] - rab_x[i] * t_0_yz_yy_xx[i];

        t_x_yz_xz_zz[i] = t_0_xyz_xz_zz[i] - rab_x[i] * t_0_yz_xz_zz[i];

        t_x_yz_xz_yz[i] = t_0_xyz_xz_yz[i] - rab_x[i] * t_0_yz_xz_yz[i];

        t_x_yz_xz_yy[i] = t_0_xyz_xz_yy[i] - rab_x[i] * t_0_yz_xz_yy[i];

        t_x_yz_xz_xz[i] = t_0_xyz_xz_xz[i] - rab_x[i] * t_0_yz_xz_xz[i];

        t_x_yz_xz_xy[i] = t_0_xyz_xz_xy[i] - rab_x[i] * t_0_yz_xz_xy[i];

        t_x_yz_xz_xx[i] = t_0_xyz_xz_xx[i] - rab_x[i] * t_0_yz_xz_xx[i];

        t_x_yz_xy_zz[i] = t_0_xyz_xy_zz[i] - rab_x[i] * t_0_yz_xy_zz[i];

        t_x_yz_xy_yz[i] = t_0_xyz_xy_yz[i] - rab_x[i] * t_0_yz_xy_yz[i];

        t_x_yz_xy_yy[i] = t_0_xyz_xy_yy[i] - rab_x[i] * t_0_yz_xy_yy[i];

        t_x_yz_xy_xz[i] = t_0_xyz_xy_xz[i] - rab_x[i] * t_0_yz_xy_xz[i];

        t_x_yz_xy_xy[i] = t_0_xyz_xy_xy[i] - rab_x[i] * t_0_yz_xy_xy[i];

        t_x_yz_xy_xx[i] = t_0_xyz_xy_xx[i] - rab_x[i] * t_0_yz_xy_xx[i];

        t_x_yz_xx_zz[i] = t_0_xyz_xx_zz[i] - rab_x[i] * t_0_yz_xx_zz[i];

        t_x_yz_xx_yz[i] = t_0_xyz_xx_yz[i] - rab_x[i] * t_0_yz_xx_yz[i];

        t_x_yz_xx_yy[i] = t_0_xyz_xx_yy[i] - rab_x[i] * t_0_yz_xx_yy[i];

        t_x_yz_xx_xz[i] = t_0_xyz_xx_xz[i] - rab_x[i] * t_0_yz_xx_xz[i];

        t_x_yz_xx_xy[i] = t_0_xyz_xx_xy[i] - rab_x[i] * t_0_yz_xx_xy[i];

        t_x_yz_xx_xx[i] = t_0_xyz_xx_xx[i] - rab_x[i] * t_0_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyy_xx_xx, t_0_xyy_xx_xy, t_0_xyy_xx_xz, t_0_xyy_xx_yy,\
                           t_0_xyy_xx_yz, t_0_xyy_xx_zz, t_0_xyy_xy_xx, t_0_xyy_xy_xy,\
                           t_0_xyy_xy_xz, t_0_xyy_xy_yy, t_0_xyy_xy_yz, t_0_xyy_xy_zz,\
                           t_0_xyy_xz_xx, t_0_xyy_xz_xy, t_0_xyy_xz_xz, t_0_xyy_xz_yy,\
                           t_0_xyy_xz_yz, t_0_xyy_xz_zz, t_0_xyy_yy_xx, t_0_xyy_yy_xy,\
                           t_0_xyy_yy_xz, t_0_xyy_yy_yy, t_0_xyy_yy_yz, t_0_xyy_yy_zz,\
                           t_0_xyy_yz_xx, t_0_xyy_yz_xy, t_0_xyy_yz_xz, t_0_xyy_yz_yy,\
                           t_0_xyy_yz_yz, t_0_xyy_yz_zz, t_0_xyy_zz_xx, t_0_xyy_zz_xy,\
                           t_0_xyy_zz_xz, t_0_xyy_zz_yy, t_0_xyy_zz_yz, t_0_xyy_zz_zz,\
                           t_0_yy_xx_xx, t_0_yy_xx_xy, t_0_yy_xx_xz, t_0_yy_xx_yy, t_0_yy_xx_yz,\
                           t_0_yy_xx_zz, t_0_yy_xy_xx, t_0_yy_xy_xy, t_0_yy_xy_xz, t_0_yy_xy_yy,\
                           t_0_yy_xy_yz, t_0_yy_xy_zz, t_0_yy_xz_xx, t_0_yy_xz_xy, t_0_yy_xz_xz,\
                           t_0_yy_xz_yy, t_0_yy_xz_yz, t_0_yy_xz_zz, t_0_yy_yy_xx, t_0_yy_yy_xy,\
                           t_0_yy_yy_xz, t_0_yy_yy_yy, t_0_yy_yy_yz, t_0_yy_yy_zz, t_0_yy_yz_xx,\
                           t_0_yy_yz_xy, t_0_yy_yz_xz, t_0_yy_yz_yy, t_0_yy_yz_yz, t_0_yy_yz_zz,\
                           t_0_yy_zz_xx, t_0_yy_zz_xy, t_0_yy_zz_xz, t_0_yy_zz_yy, t_0_yy_zz_yz,\
                           t_0_yy_zz_zz, t_x_yy_xx_xx, t_x_yy_xx_xy, t_x_yy_xx_xz, t_x_yy_xx_yy,\
                           t_x_yy_xx_yz, t_x_yy_xx_zz, t_x_yy_xy_xx, t_x_yy_xy_xy, t_x_yy_xy_xz,\
                           t_x_yy_xy_yy, t_x_yy_xy_yz, t_x_yy_xy_zz, t_x_yy_xz_xx, t_x_yy_xz_xy,\
                           t_x_yy_xz_xz, t_x_yy_xz_yy, t_x_yy_xz_yz, t_x_yy_xz_zz, t_x_yy_yy_xx,\
                           t_x_yy_yy_xy, t_x_yy_yy_xz, t_x_yy_yy_yy, t_x_yy_yy_yz, t_x_yy_yy_zz,\
                           t_x_yy_yz_xx, t_x_yy_yz_xy, t_x_yy_yz_xz, t_x_yy_yz_yy, t_x_yy_yz_yz,\
                           t_x_yy_yz_zz, t_x_yy_zz_xx, t_x_yy_zz_xy, t_x_yy_zz_xz, t_x_yy_zz_yy,\
                           t_x_yy_zz_yz, t_x_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yy_zz_zz[i] = t_0_xyy_zz_zz[i] - rab_x[i] * t_0_yy_zz_zz[i];

        t_x_yy_zz_yz[i] = t_0_xyy_zz_yz[i] - rab_x[i] * t_0_yy_zz_yz[i];

        t_x_yy_zz_yy[i] = t_0_xyy_zz_yy[i] - rab_x[i] * t_0_yy_zz_yy[i];

        t_x_yy_zz_xz[i] = t_0_xyy_zz_xz[i] - rab_x[i] * t_0_yy_zz_xz[i];

        t_x_yy_zz_xy[i] = t_0_xyy_zz_xy[i] - rab_x[i] * t_0_yy_zz_xy[i];

        t_x_yy_zz_xx[i] = t_0_xyy_zz_xx[i] - rab_x[i] * t_0_yy_zz_xx[i];

        t_x_yy_yz_zz[i] = t_0_xyy_yz_zz[i] - rab_x[i] * t_0_yy_yz_zz[i];

        t_x_yy_yz_yz[i] = t_0_xyy_yz_yz[i] - rab_x[i] * t_0_yy_yz_yz[i];

        t_x_yy_yz_yy[i] = t_0_xyy_yz_yy[i] - rab_x[i] * t_0_yy_yz_yy[i];

        t_x_yy_yz_xz[i] = t_0_xyy_yz_xz[i] - rab_x[i] * t_0_yy_yz_xz[i];

        t_x_yy_yz_xy[i] = t_0_xyy_yz_xy[i] - rab_x[i] * t_0_yy_yz_xy[i];

        t_x_yy_yz_xx[i] = t_0_xyy_yz_xx[i] - rab_x[i] * t_0_yy_yz_xx[i];

        t_x_yy_yy_zz[i] = t_0_xyy_yy_zz[i] - rab_x[i] * t_0_yy_yy_zz[i];

        t_x_yy_yy_yz[i] = t_0_xyy_yy_yz[i] - rab_x[i] * t_0_yy_yy_yz[i];

        t_x_yy_yy_yy[i] = t_0_xyy_yy_yy[i] - rab_x[i] * t_0_yy_yy_yy[i];

        t_x_yy_yy_xz[i] = t_0_xyy_yy_xz[i] - rab_x[i] * t_0_yy_yy_xz[i];

        t_x_yy_yy_xy[i] = t_0_xyy_yy_xy[i] - rab_x[i] * t_0_yy_yy_xy[i];

        t_x_yy_yy_xx[i] = t_0_xyy_yy_xx[i] - rab_x[i] * t_0_yy_yy_xx[i];

        t_x_yy_xz_zz[i] = t_0_xyy_xz_zz[i] - rab_x[i] * t_0_yy_xz_zz[i];

        t_x_yy_xz_yz[i] = t_0_xyy_xz_yz[i] - rab_x[i] * t_0_yy_xz_yz[i];

        t_x_yy_xz_yy[i] = t_0_xyy_xz_yy[i] - rab_x[i] * t_0_yy_xz_yy[i];

        t_x_yy_xz_xz[i] = t_0_xyy_xz_xz[i] - rab_x[i] * t_0_yy_xz_xz[i];

        t_x_yy_xz_xy[i] = t_0_xyy_xz_xy[i] - rab_x[i] * t_0_yy_xz_xy[i];

        t_x_yy_xz_xx[i] = t_0_xyy_xz_xx[i] - rab_x[i] * t_0_yy_xz_xx[i];

        t_x_yy_xy_zz[i] = t_0_xyy_xy_zz[i] - rab_x[i] * t_0_yy_xy_zz[i];

        t_x_yy_xy_yz[i] = t_0_xyy_xy_yz[i] - rab_x[i] * t_0_yy_xy_yz[i];

        t_x_yy_xy_yy[i] = t_0_xyy_xy_yy[i] - rab_x[i] * t_0_yy_xy_yy[i];

        t_x_yy_xy_xz[i] = t_0_xyy_xy_xz[i] - rab_x[i] * t_0_yy_xy_xz[i];

        t_x_yy_xy_xy[i] = t_0_xyy_xy_xy[i] - rab_x[i] * t_0_yy_xy_xy[i];

        t_x_yy_xy_xx[i] = t_0_xyy_xy_xx[i] - rab_x[i] * t_0_yy_xy_xx[i];

        t_x_yy_xx_zz[i] = t_0_xyy_xx_zz[i] - rab_x[i] * t_0_yy_xx_zz[i];

        t_x_yy_xx_yz[i] = t_0_xyy_xx_yz[i] - rab_x[i] * t_0_yy_xx_yz[i];

        t_x_yy_xx_yy[i] = t_0_xyy_xx_yy[i] - rab_x[i] * t_0_yy_xx_yy[i];

        t_x_yy_xx_xz[i] = t_0_xyy_xx_xz[i] - rab_x[i] * t_0_yy_xx_xz[i];

        t_x_yy_xx_xy[i] = t_0_xyy_xx_xy[i] - rab_x[i] * t_0_yy_xx_xy[i];

        t_x_yy_xx_xx[i] = t_0_xyy_xx_xx[i] - rab_x[i] * t_0_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxz_xx_xx, t_0_xxz_xx_xy, t_0_xxz_xx_xz, t_0_xxz_xx_yy,\
                           t_0_xxz_xx_yz, t_0_xxz_xx_zz, t_0_xxz_xy_xx, t_0_xxz_xy_xy,\
                           t_0_xxz_xy_xz, t_0_xxz_xy_yy, t_0_xxz_xy_yz, t_0_xxz_xy_zz,\
                           t_0_xxz_xz_xx, t_0_xxz_xz_xy, t_0_xxz_xz_xz, t_0_xxz_xz_yy,\
                           t_0_xxz_xz_yz, t_0_xxz_xz_zz, t_0_xxz_yy_xx, t_0_xxz_yy_xy,\
                           t_0_xxz_yy_xz, t_0_xxz_yy_yy, t_0_xxz_yy_yz, t_0_xxz_yy_zz,\
                           t_0_xxz_yz_xx, t_0_xxz_yz_xy, t_0_xxz_yz_xz, t_0_xxz_yz_yy,\
                           t_0_xxz_yz_yz, t_0_xxz_yz_zz, t_0_xxz_zz_xx, t_0_xxz_zz_xy,\
                           t_0_xxz_zz_xz, t_0_xxz_zz_yy, t_0_xxz_zz_yz, t_0_xxz_zz_zz,\
                           t_0_xz_xx_xx, t_0_xz_xx_xy, t_0_xz_xx_xz, t_0_xz_xx_yy, t_0_xz_xx_yz,\
                           t_0_xz_xx_zz, t_0_xz_xy_xx, t_0_xz_xy_xy, t_0_xz_xy_xz, t_0_xz_xy_yy,\
                           t_0_xz_xy_yz, t_0_xz_xy_zz, t_0_xz_xz_xx, t_0_xz_xz_xy, t_0_xz_xz_xz,\
                           t_0_xz_xz_yy, t_0_xz_xz_yz, t_0_xz_xz_zz, t_0_xz_yy_xx, t_0_xz_yy_xy,\
                           t_0_xz_yy_xz, t_0_xz_yy_yy, t_0_xz_yy_yz, t_0_xz_yy_zz, t_0_xz_yz_xx,\
                           t_0_xz_yz_xy, t_0_xz_yz_xz, t_0_xz_yz_yy, t_0_xz_yz_yz, t_0_xz_yz_zz,\
                           t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz, t_0_xz_zz_yy, t_0_xz_zz_yz,\
                           t_0_xz_zz_zz, t_x_xz_xx_xx, t_x_xz_xx_xy, t_x_xz_xx_xz, t_x_xz_xx_yy,\
                           t_x_xz_xx_yz, t_x_xz_xx_zz, t_x_xz_xy_xx, t_x_xz_xy_xy, t_x_xz_xy_xz,\
                           t_x_xz_xy_yy, t_x_xz_xy_yz, t_x_xz_xy_zz, t_x_xz_xz_xx, t_x_xz_xz_xy,\
                           t_x_xz_xz_xz, t_x_xz_xz_yy, t_x_xz_xz_yz, t_x_xz_xz_zz, t_x_xz_yy_xx,\
                           t_x_xz_yy_xy, t_x_xz_yy_xz, t_x_xz_yy_yy, t_x_xz_yy_yz, t_x_xz_yy_zz,\
                           t_x_xz_yz_xx, t_x_xz_yz_xy, t_x_xz_yz_xz, t_x_xz_yz_yy, t_x_xz_yz_yz,\
                           t_x_xz_yz_zz, t_x_xz_zz_xx, t_x_xz_zz_xy, t_x_xz_zz_xz, t_x_xz_zz_yy,\
                           t_x_xz_zz_yz, t_x_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xz_zz_zz[i] = t_0_xxz_zz_zz[i] - rab_x[i] * t_0_xz_zz_zz[i];

        t_x_xz_zz_yz[i] = t_0_xxz_zz_yz[i] - rab_x[i] * t_0_xz_zz_yz[i];

        t_x_xz_zz_yy[i] = t_0_xxz_zz_yy[i] - rab_x[i] * t_0_xz_zz_yy[i];

        t_x_xz_zz_xz[i] = t_0_xxz_zz_xz[i] - rab_x[i] * t_0_xz_zz_xz[i];

        t_x_xz_zz_xy[i] = t_0_xxz_zz_xy[i] - rab_x[i] * t_0_xz_zz_xy[i];

        t_x_xz_zz_xx[i] = t_0_xxz_zz_xx[i] - rab_x[i] * t_0_xz_zz_xx[i];

        t_x_xz_yz_zz[i] = t_0_xxz_yz_zz[i] - rab_x[i] * t_0_xz_yz_zz[i];

        t_x_xz_yz_yz[i] = t_0_xxz_yz_yz[i] - rab_x[i] * t_0_xz_yz_yz[i];

        t_x_xz_yz_yy[i] = t_0_xxz_yz_yy[i] - rab_x[i] * t_0_xz_yz_yy[i];

        t_x_xz_yz_xz[i] = t_0_xxz_yz_xz[i] - rab_x[i] * t_0_xz_yz_xz[i];

        t_x_xz_yz_xy[i] = t_0_xxz_yz_xy[i] - rab_x[i] * t_0_xz_yz_xy[i];

        t_x_xz_yz_xx[i] = t_0_xxz_yz_xx[i] - rab_x[i] * t_0_xz_yz_xx[i];

        t_x_xz_yy_zz[i] = t_0_xxz_yy_zz[i] - rab_x[i] * t_0_xz_yy_zz[i];

        t_x_xz_yy_yz[i] = t_0_xxz_yy_yz[i] - rab_x[i] * t_0_xz_yy_yz[i];

        t_x_xz_yy_yy[i] = t_0_xxz_yy_yy[i] - rab_x[i] * t_0_xz_yy_yy[i];

        t_x_xz_yy_xz[i] = t_0_xxz_yy_xz[i] - rab_x[i] * t_0_xz_yy_xz[i];

        t_x_xz_yy_xy[i] = t_0_xxz_yy_xy[i] - rab_x[i] * t_0_xz_yy_xy[i];

        t_x_xz_yy_xx[i] = t_0_xxz_yy_xx[i] - rab_x[i] * t_0_xz_yy_xx[i];

        t_x_xz_xz_zz[i] = t_0_xxz_xz_zz[i] - rab_x[i] * t_0_xz_xz_zz[i];

        t_x_xz_xz_yz[i] = t_0_xxz_xz_yz[i] - rab_x[i] * t_0_xz_xz_yz[i];

        t_x_xz_xz_yy[i] = t_0_xxz_xz_yy[i] - rab_x[i] * t_0_xz_xz_yy[i];

        t_x_xz_xz_xz[i] = t_0_xxz_xz_xz[i] - rab_x[i] * t_0_xz_xz_xz[i];

        t_x_xz_xz_xy[i] = t_0_xxz_xz_xy[i] - rab_x[i] * t_0_xz_xz_xy[i];

        t_x_xz_xz_xx[i] = t_0_xxz_xz_xx[i] - rab_x[i] * t_0_xz_xz_xx[i];

        t_x_xz_xy_zz[i] = t_0_xxz_xy_zz[i] - rab_x[i] * t_0_xz_xy_zz[i];

        t_x_xz_xy_yz[i] = t_0_xxz_xy_yz[i] - rab_x[i] * t_0_xz_xy_yz[i];

        t_x_xz_xy_yy[i] = t_0_xxz_xy_yy[i] - rab_x[i] * t_0_xz_xy_yy[i];

        t_x_xz_xy_xz[i] = t_0_xxz_xy_xz[i] - rab_x[i] * t_0_xz_xy_xz[i];

        t_x_xz_xy_xy[i] = t_0_xxz_xy_xy[i] - rab_x[i] * t_0_xz_xy_xy[i];

        t_x_xz_xy_xx[i] = t_0_xxz_xy_xx[i] - rab_x[i] * t_0_xz_xy_xx[i];

        t_x_xz_xx_zz[i] = t_0_xxz_xx_zz[i] - rab_x[i] * t_0_xz_xx_zz[i];

        t_x_xz_xx_yz[i] = t_0_xxz_xx_yz[i] - rab_x[i] * t_0_xz_xx_yz[i];

        t_x_xz_xx_yy[i] = t_0_xxz_xx_yy[i] - rab_x[i] * t_0_xz_xx_yy[i];

        t_x_xz_xx_xz[i] = t_0_xxz_xx_xz[i] - rab_x[i] * t_0_xz_xx_xz[i];

        t_x_xz_xx_xy[i] = t_0_xxz_xx_xy[i] - rab_x[i] * t_0_xz_xx_xy[i];

        t_x_xz_xx_xx[i] = t_0_xxz_xx_xx[i] - rab_x[i] * t_0_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxy_xx_xx, t_0_xxy_xx_xy, t_0_xxy_xx_xz, t_0_xxy_xx_yy,\
                           t_0_xxy_xx_yz, t_0_xxy_xx_zz, t_0_xxy_xy_xx, t_0_xxy_xy_xy,\
                           t_0_xxy_xy_xz, t_0_xxy_xy_yy, t_0_xxy_xy_yz, t_0_xxy_xy_zz,\
                           t_0_xxy_xz_xx, t_0_xxy_xz_xy, t_0_xxy_xz_xz, t_0_xxy_xz_yy,\
                           t_0_xxy_xz_yz, t_0_xxy_xz_zz, t_0_xxy_yy_xx, t_0_xxy_yy_xy,\
                           t_0_xxy_yy_xz, t_0_xxy_yy_yy, t_0_xxy_yy_yz, t_0_xxy_yy_zz,\
                           t_0_xxy_yz_xx, t_0_xxy_yz_xy, t_0_xxy_yz_xz, t_0_xxy_yz_yy,\
                           t_0_xxy_yz_yz, t_0_xxy_yz_zz, t_0_xxy_zz_xx, t_0_xxy_zz_xy,\
                           t_0_xxy_zz_xz, t_0_xxy_zz_yy, t_0_xxy_zz_yz, t_0_xxy_zz_zz,\
                           t_0_xy_xx_xx, t_0_xy_xx_xy, t_0_xy_xx_xz, t_0_xy_xx_yy, t_0_xy_xx_yz,\
                           t_0_xy_xx_zz, t_0_xy_xy_xx, t_0_xy_xy_xy, t_0_xy_xy_xz, t_0_xy_xy_yy,\
                           t_0_xy_xy_yz, t_0_xy_xy_zz, t_0_xy_xz_xx, t_0_xy_xz_xy, t_0_xy_xz_xz,\
                           t_0_xy_xz_yy, t_0_xy_xz_yz, t_0_xy_xz_zz, t_0_xy_yy_xx, t_0_xy_yy_xy,\
                           t_0_xy_yy_xz, t_0_xy_yy_yy, t_0_xy_yy_yz, t_0_xy_yy_zz, t_0_xy_yz_xx,\
                           t_0_xy_yz_xy, t_0_xy_yz_xz, t_0_xy_yz_yy, t_0_xy_yz_yz, t_0_xy_yz_zz,\
                           t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz, t_0_xy_zz_yy, t_0_xy_zz_yz,\
                           t_0_xy_zz_zz, t_x_xy_xx_xx, t_x_xy_xx_xy, t_x_xy_xx_xz, t_x_xy_xx_yy,\
                           t_x_xy_xx_yz, t_x_xy_xx_zz, t_x_xy_xy_xx, t_x_xy_xy_xy, t_x_xy_xy_xz,\
                           t_x_xy_xy_yy, t_x_xy_xy_yz, t_x_xy_xy_zz, t_x_xy_xz_xx, t_x_xy_xz_xy,\
                           t_x_xy_xz_xz, t_x_xy_xz_yy, t_x_xy_xz_yz, t_x_xy_xz_zz, t_x_xy_yy_xx,\
                           t_x_xy_yy_xy, t_x_xy_yy_xz, t_x_xy_yy_yy, t_x_xy_yy_yz, t_x_xy_yy_zz,\
                           t_x_xy_yz_xx, t_x_xy_yz_xy, t_x_xy_yz_xz, t_x_xy_yz_yy, t_x_xy_yz_yz,\
                           t_x_xy_yz_zz, t_x_xy_zz_xx, t_x_xy_zz_xy, t_x_xy_zz_xz, t_x_xy_zz_yy,\
                           t_x_xy_zz_yz, t_x_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xy_zz_zz[i] = t_0_xxy_zz_zz[i] - rab_x[i] * t_0_xy_zz_zz[i];

        t_x_xy_zz_yz[i] = t_0_xxy_zz_yz[i] - rab_x[i] * t_0_xy_zz_yz[i];

        t_x_xy_zz_yy[i] = t_0_xxy_zz_yy[i] - rab_x[i] * t_0_xy_zz_yy[i];

        t_x_xy_zz_xz[i] = t_0_xxy_zz_xz[i] - rab_x[i] * t_0_xy_zz_xz[i];

        t_x_xy_zz_xy[i] = t_0_xxy_zz_xy[i] - rab_x[i] * t_0_xy_zz_xy[i];

        t_x_xy_zz_xx[i] = t_0_xxy_zz_xx[i] - rab_x[i] * t_0_xy_zz_xx[i];

        t_x_xy_yz_zz[i] = t_0_xxy_yz_zz[i] - rab_x[i] * t_0_xy_yz_zz[i];

        t_x_xy_yz_yz[i] = t_0_xxy_yz_yz[i] - rab_x[i] * t_0_xy_yz_yz[i];

        t_x_xy_yz_yy[i] = t_0_xxy_yz_yy[i] - rab_x[i] * t_0_xy_yz_yy[i];

        t_x_xy_yz_xz[i] = t_0_xxy_yz_xz[i] - rab_x[i] * t_0_xy_yz_xz[i];

        t_x_xy_yz_xy[i] = t_0_xxy_yz_xy[i] - rab_x[i] * t_0_xy_yz_xy[i];

        t_x_xy_yz_xx[i] = t_0_xxy_yz_xx[i] - rab_x[i] * t_0_xy_yz_xx[i];

        t_x_xy_yy_zz[i] = t_0_xxy_yy_zz[i] - rab_x[i] * t_0_xy_yy_zz[i];

        t_x_xy_yy_yz[i] = t_0_xxy_yy_yz[i] - rab_x[i] * t_0_xy_yy_yz[i];

        t_x_xy_yy_yy[i] = t_0_xxy_yy_yy[i] - rab_x[i] * t_0_xy_yy_yy[i];

        t_x_xy_yy_xz[i] = t_0_xxy_yy_xz[i] - rab_x[i] * t_0_xy_yy_xz[i];

        t_x_xy_yy_xy[i] = t_0_xxy_yy_xy[i] - rab_x[i] * t_0_xy_yy_xy[i];

        t_x_xy_yy_xx[i] = t_0_xxy_yy_xx[i] - rab_x[i] * t_0_xy_yy_xx[i];

        t_x_xy_xz_zz[i] = t_0_xxy_xz_zz[i] - rab_x[i] * t_0_xy_xz_zz[i];

        t_x_xy_xz_yz[i] = t_0_xxy_xz_yz[i] - rab_x[i] * t_0_xy_xz_yz[i];

        t_x_xy_xz_yy[i] = t_0_xxy_xz_yy[i] - rab_x[i] * t_0_xy_xz_yy[i];

        t_x_xy_xz_xz[i] = t_0_xxy_xz_xz[i] - rab_x[i] * t_0_xy_xz_xz[i];

        t_x_xy_xz_xy[i] = t_0_xxy_xz_xy[i] - rab_x[i] * t_0_xy_xz_xy[i];

        t_x_xy_xz_xx[i] = t_0_xxy_xz_xx[i] - rab_x[i] * t_0_xy_xz_xx[i];

        t_x_xy_xy_zz[i] = t_0_xxy_xy_zz[i] - rab_x[i] * t_0_xy_xy_zz[i];

        t_x_xy_xy_yz[i] = t_0_xxy_xy_yz[i] - rab_x[i] * t_0_xy_xy_yz[i];

        t_x_xy_xy_yy[i] = t_0_xxy_xy_yy[i] - rab_x[i] * t_0_xy_xy_yy[i];

        t_x_xy_xy_xz[i] = t_0_xxy_xy_xz[i] - rab_x[i] * t_0_xy_xy_xz[i];

        t_x_xy_xy_xy[i] = t_0_xxy_xy_xy[i] - rab_x[i] * t_0_xy_xy_xy[i];

        t_x_xy_xy_xx[i] = t_0_xxy_xy_xx[i] - rab_x[i] * t_0_xy_xy_xx[i];

        t_x_xy_xx_zz[i] = t_0_xxy_xx_zz[i] - rab_x[i] * t_0_xy_xx_zz[i];

        t_x_xy_xx_yz[i] = t_0_xxy_xx_yz[i] - rab_x[i] * t_0_xy_xx_yz[i];

        t_x_xy_xx_yy[i] = t_0_xxy_xx_yy[i] - rab_x[i] * t_0_xy_xx_yy[i];

        t_x_xy_xx_xz[i] = t_0_xxy_xx_xz[i] - rab_x[i] * t_0_xy_xx_xz[i];

        t_x_xy_xx_xy[i] = t_0_xxy_xx_xy[i] - rab_x[i] * t_0_xy_xx_xy[i];

        t_x_xy_xx_xx[i] = t_0_xxy_xx_xx[i] - rab_x[i] * t_0_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xx_xx_xx, t_0_xx_xx_xy, t_0_xx_xx_xz, t_0_xx_xx_yy,\
                           t_0_xx_xx_yz, t_0_xx_xx_zz, t_0_xx_xy_xx, t_0_xx_xy_xy, t_0_xx_xy_xz,\
                           t_0_xx_xy_yy, t_0_xx_xy_yz, t_0_xx_xy_zz, t_0_xx_xz_xx, t_0_xx_xz_xy,\
                           t_0_xx_xz_xz, t_0_xx_xz_yy, t_0_xx_xz_yz, t_0_xx_xz_zz, t_0_xx_yy_xx,\
                           t_0_xx_yy_xy, t_0_xx_yy_xz, t_0_xx_yy_yy, t_0_xx_yy_yz, t_0_xx_yy_zz,\
                           t_0_xx_yz_xx, t_0_xx_yz_xy, t_0_xx_yz_xz, t_0_xx_yz_yy, t_0_xx_yz_yz,\
                           t_0_xx_yz_zz, t_0_xx_zz_xx, t_0_xx_zz_xy, t_0_xx_zz_xz, t_0_xx_zz_yy,\
                           t_0_xx_zz_yz, t_0_xx_zz_zz, t_0_xxx_xx_xx, t_0_xxx_xx_xy,\
                           t_0_xxx_xx_xz, t_0_xxx_xx_yy, t_0_xxx_xx_yz, t_0_xxx_xx_zz,\
                           t_0_xxx_xy_xx, t_0_xxx_xy_xy, t_0_xxx_xy_xz, t_0_xxx_xy_yy,\
                           t_0_xxx_xy_yz, t_0_xxx_xy_zz, t_0_xxx_xz_xx, t_0_xxx_xz_xy,\
                           t_0_xxx_xz_xz, t_0_xxx_xz_yy, t_0_xxx_xz_yz, t_0_xxx_xz_zz,\
                           t_0_xxx_yy_xx, t_0_xxx_yy_xy, t_0_xxx_yy_xz, t_0_xxx_yy_yy,\
                           t_0_xxx_yy_yz, t_0_xxx_yy_zz, t_0_xxx_yz_xx, t_0_xxx_yz_xy,\
                           t_0_xxx_yz_xz, t_0_xxx_yz_yy, t_0_xxx_yz_yz, t_0_xxx_yz_zz,\
                           t_0_xxx_zz_xx, t_0_xxx_zz_xy, t_0_xxx_zz_xz, t_0_xxx_zz_yy,\
                           t_0_xxx_zz_yz, t_0_xxx_zz_zz, t_x_xx_xx_xx, t_x_xx_xx_xy,\
                           t_x_xx_xx_xz, t_x_xx_xx_yy, t_x_xx_xx_yz, t_x_xx_xx_zz, t_x_xx_xy_xx,\
                           t_x_xx_xy_xy, t_x_xx_xy_xz, t_x_xx_xy_yy, t_x_xx_xy_yz, t_x_xx_xy_zz,\
                           t_x_xx_xz_xx, t_x_xx_xz_xy, t_x_xx_xz_xz, t_x_xx_xz_yy, t_x_xx_xz_yz,\
                           t_x_xx_xz_zz, t_x_xx_yy_xx, t_x_xx_yy_xy, t_x_xx_yy_xz, t_x_xx_yy_yy,\
                           t_x_xx_yy_yz, t_x_xx_yy_zz, t_x_xx_yz_xx, t_x_xx_yz_xy, t_x_xx_yz_xz,\
                           t_x_xx_yz_yy, t_x_xx_yz_yz, t_x_xx_yz_zz, t_x_xx_zz_xx, t_x_xx_zz_xy,\
                           t_x_xx_zz_xz, t_x_xx_zz_yy, t_x_xx_zz_yz, t_x_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xx_zz_zz[i] = t_0_xxx_zz_zz[i] - rab_x[i] * t_0_xx_zz_zz[i];

        t_x_xx_zz_yz[i] = t_0_xxx_zz_yz[i] - rab_x[i] * t_0_xx_zz_yz[i];

        t_x_xx_zz_yy[i] = t_0_xxx_zz_yy[i] - rab_x[i] * t_0_xx_zz_yy[i];

        t_x_xx_zz_xz[i] = t_0_xxx_zz_xz[i] - rab_x[i] * t_0_xx_zz_xz[i];

        t_x_xx_zz_xy[i] = t_0_xxx_zz_xy[i] - rab_x[i] * t_0_xx_zz_xy[i];

        t_x_xx_zz_xx[i] = t_0_xxx_zz_xx[i] - rab_x[i] * t_0_xx_zz_xx[i];

        t_x_xx_yz_zz[i] = t_0_xxx_yz_zz[i] - rab_x[i] * t_0_xx_yz_zz[i];

        t_x_xx_yz_yz[i] = t_0_xxx_yz_yz[i] - rab_x[i] * t_0_xx_yz_yz[i];

        t_x_xx_yz_yy[i] = t_0_xxx_yz_yy[i] - rab_x[i] * t_0_xx_yz_yy[i];

        t_x_xx_yz_xz[i] = t_0_xxx_yz_xz[i] - rab_x[i] * t_0_xx_yz_xz[i];

        t_x_xx_yz_xy[i] = t_0_xxx_yz_xy[i] - rab_x[i] * t_0_xx_yz_xy[i];

        t_x_xx_yz_xx[i] = t_0_xxx_yz_xx[i] - rab_x[i] * t_0_xx_yz_xx[i];

        t_x_xx_yy_zz[i] = t_0_xxx_yy_zz[i] - rab_x[i] * t_0_xx_yy_zz[i];

        t_x_xx_yy_yz[i] = t_0_xxx_yy_yz[i] - rab_x[i] * t_0_xx_yy_yz[i];

        t_x_xx_yy_yy[i] = t_0_xxx_yy_yy[i] - rab_x[i] * t_0_xx_yy_yy[i];

        t_x_xx_yy_xz[i] = t_0_xxx_yy_xz[i] - rab_x[i] * t_0_xx_yy_xz[i];

        t_x_xx_yy_xy[i] = t_0_xxx_yy_xy[i] - rab_x[i] * t_0_xx_yy_xy[i];

        t_x_xx_yy_xx[i] = t_0_xxx_yy_xx[i] - rab_x[i] * t_0_xx_yy_xx[i];

        t_x_xx_xz_zz[i] = t_0_xxx_xz_zz[i] - rab_x[i] * t_0_xx_xz_zz[i];

        t_x_xx_xz_yz[i] = t_0_xxx_xz_yz[i] - rab_x[i] * t_0_xx_xz_yz[i];

        t_x_xx_xz_yy[i] = t_0_xxx_xz_yy[i] - rab_x[i] * t_0_xx_xz_yy[i];

        t_x_xx_xz_xz[i] = t_0_xxx_xz_xz[i] - rab_x[i] * t_0_xx_xz_xz[i];

        t_x_xx_xz_xy[i] = t_0_xxx_xz_xy[i] - rab_x[i] * t_0_xx_xz_xy[i];

        t_x_xx_xz_xx[i] = t_0_xxx_xz_xx[i] - rab_x[i] * t_0_xx_xz_xx[i];

        t_x_xx_xy_zz[i] = t_0_xxx_xy_zz[i] - rab_x[i] * t_0_xx_xy_zz[i];

        t_x_xx_xy_yz[i] = t_0_xxx_xy_yz[i] - rab_x[i] * t_0_xx_xy_yz[i];

        t_x_xx_xy_yy[i] = t_0_xxx_xy_yy[i] - rab_x[i] * t_0_xx_xy_yy[i];

        t_x_xx_xy_xz[i] = t_0_xxx_xy_xz[i] - rab_x[i] * t_0_xx_xy_xz[i];

        t_x_xx_xy_xy[i] = t_0_xxx_xy_xy[i] - rab_x[i] * t_0_xx_xy_xy[i];

        t_x_xx_xy_xx[i] = t_0_xxx_xy_xx[i] - rab_x[i] * t_0_xx_xy_xx[i];

        t_x_xx_xx_zz[i] = t_0_xxx_xx_zz[i] - rab_x[i] * t_0_xx_xx_zz[i];

        t_x_xx_xx_yz[i] = t_0_xxx_xx_yz[i] - rab_x[i] * t_0_xx_xx_yz[i];

        t_x_xx_xx_yy[i] = t_0_xxx_xx_yy[i] - rab_x[i] * t_0_xx_xx_yy[i];

        t_x_xx_xx_xz[i] = t_0_xxx_xx_xz[i] - rab_x[i] * t_0_xx_xx_xz[i];

        t_x_xx_xx_xy[i] = t_0_xxx_xx_xy[i] - rab_x[i] * t_0_xx_xx_xy[i];

        t_x_xx_xx_xx[i] = t_0_xxx_xx_xx[i] - rab_x[i] * t_0_xx_xx_xx[i];
    }
}


} // derirec namespace
