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
compHostHRRForDDDD_V0(      BufferHostXY<T>&      intsBufferDDDD,
                      const BufferHostX<int32_t>& intsIndexesDDDD,
                      const BufferHostXY<T>&      intsBufferPDDD,
                      const BufferHostX<int32_t>& intsIndexesPDDD,
                      const BufferHostXY<T>&      intsBufferPFDD,
                      const BufferHostX<int32_t>& intsIndexesPFDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDDD) integral components

    t_zz_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(0));

    t_zz_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1));

    t_zz_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(2));

    t_zz_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(3));

    t_zz_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(4));

    t_zz_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(5));

    t_zz_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(6));

    t_zz_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(7));

    t_zz_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(8));

    t_zz_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(9));

    t_zz_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(10));

    t_zz_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(11));

    t_zz_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(12));

    t_zz_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(13));

    t_zz_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(14));

    t_zz_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(15));

    t_zz_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(16));

    t_zz_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(17));

    t_zz_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(18));

    t_zz_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(19));

    t_zz_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(20));

    t_zz_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(21));

    t_zz_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(22));

    t_zz_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(23));

    t_zz_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(24));

    t_zz_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(25));

    t_zz_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(26));

    t_zz_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(27));

    t_zz_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(28));

    t_zz_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(29));

    t_zz_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(30));

    t_zz_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(31));

    t_zz_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(32));

    t_zz_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(33));

    t_zz_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(34));

    t_zz_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(35));

    t_zz_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(36));

    t_zz_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(37));

    t_zz_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(38));

    t_zz_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(39));

    t_zz_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(40));

    t_zz_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(41));

    t_zz_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(42));

    t_zz_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(43));

    t_zz_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(44));

    t_zz_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(45));

    t_zz_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(46));

    t_zz_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(47));

    t_zz_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(48));

    t_zz_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(49));

    t_zz_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(50));

    t_zz_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(51));

    t_zz_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(52));

    t_zz_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(53));

    t_zz_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(54));

    t_zz_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(55));

    t_zz_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(56));

    t_zz_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(57));

    t_zz_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(58));

    t_zz_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(59));

    t_zz_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(60));

    t_zz_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(61));

    t_zz_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(62));

    t_zz_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(63));

    t_zz_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(64));

    t_zz_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(65));

    t_zz_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(66));

    t_zz_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(67));

    t_zz_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(68));

    t_zz_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(69));

    t_zz_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(70));

    t_zz_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(71));

    t_zz_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(72));

    t_zz_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(73));

    t_zz_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(74));

    t_zz_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(75));

    t_zz_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(76));

    t_zz_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(77));

    t_zz_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(78));

    t_zz_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(79));

    t_zz_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(80));

    t_zz_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(81));

    t_zz_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(82));

    t_zz_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(83));

    t_zz_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(84));

    t_zz_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(85));

    t_zz_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(86));

    t_zz_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(87));

    t_zz_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(88));

    t_zz_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(89));

    t_zz_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(90));

    t_zz_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(91));

    t_zz_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(92));

    t_zz_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(93));

    t_zz_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(94));

    t_zz_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(95));

    t_zz_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(96));

    t_zz_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(97));

    t_zz_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(98));

    t_zz_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(99));

    t_zz_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(100));

    t_zz_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(101));

    t_zz_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(102));

    t_zz_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(103));

    t_zz_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(104));

    t_zz_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(105));

    t_zz_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(106));

    t_zz_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(107));

    t_zz_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(108));

    t_zz_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(109));

    t_zz_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(110));

    t_zz_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(111));

    t_zz_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(112));

    t_zz_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(113));

    t_zz_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(114));

    t_zz_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(115));

    t_zz_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(116));

    t_zz_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(117));

    t_zz_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(118));

    t_zz_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(119));

    t_zz_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(120));

    t_zz_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(121));

    t_zz_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(122));

    t_zz_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(123));

    t_zz_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(124));

    t_zz_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(125));

    t_zz_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(126));

    t_zz_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(127));

    t_zz_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(128));

    t_zz_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(129));

    t_zz_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(130));

    t_zz_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(131));

    t_zz_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(132));

    t_zz_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(133));

    t_zz_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(134));

    t_zz_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(135));

    t_zz_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(136));

    t_zz_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(137));

    t_zz_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(138));

    t_zz_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(139));

    t_zz_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(140));

    t_zz_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(141));

    t_zz_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(142));

    t_zz_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(143));

    t_zz_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(144));

    t_zz_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(145));

    t_zz_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(146));

    t_zz_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(147));

    t_zz_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(148));

    t_zz_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(149));

    t_zz_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(150));

    t_zz_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(151));

    t_zz_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(152));

    t_zz_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(153));

    t_zz_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(154));

    t_zz_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(155));

    t_zz_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(156));

    t_zz_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(157));

    t_zz_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(158));

    t_zz_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(159));

    t_zz_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(160));

    t_zz_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(161));

    t_zz_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(162));

    t_zz_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(163));

    t_zz_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(164));

    t_zz_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(165));

    t_zz_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(166));

    t_zz_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(167));

    t_zz_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(168));

    t_zz_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(169));

    t_zz_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(170));

    t_zz_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(171));

    t_zz_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(172));

    t_zz_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(173));

    t_zz_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(174));

    t_zz_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(175));

    t_zz_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(176));

    t_zz_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(177));

    t_zz_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(178));

    t_zz_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(179));

    t_zz_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(180));

    t_zz_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(181));

    t_zz_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(182));

    t_zz_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(183));

    t_zz_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(184));

    t_zz_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(185));

    t_zz_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(186));

    t_zz_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(187));

    t_zz_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(188));

    t_zz_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(189));

    t_zz_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(190));

    t_zz_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(191));

    t_zz_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(192));

    t_zz_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(193));

    t_zz_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(194));

    t_zz_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(195));

    t_zz_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(196));

    t_zz_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(197));

    t_zz_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(198));

    t_zz_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(199));

    t_zz_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(200));

    t_zz_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(201));

    t_zz_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(202));

    t_zz_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(203));

    t_zz_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(204));

    t_zz_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(205));

    t_zz_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(206));

    t_zz_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(207));

    t_zz_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(208));

    t_zz_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(209));

    t_zz_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(210));

    t_zz_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(211));

    t_zz_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(212));

    t_zz_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(213));

    t_zz_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(214));

    t_zz_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(215));

    t_yz_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(216));

    t_yz_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(217));

    t_yz_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(218));

    t_yz_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(219));

    t_yz_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(220));

    t_yz_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(221));

    t_yz_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(222));

    t_yz_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(223));

    t_yz_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(224));

    t_yz_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(225));

    t_yz_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(226));

    t_yz_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(227));

    t_yz_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(228));

    t_yz_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(229));

    t_yz_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(230));

    t_yz_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(231));

    t_yz_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(232));

    t_yz_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(233));

    t_yz_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(234));

    t_yz_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(235));

    t_yz_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(236));

    t_yz_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(237));

    t_yz_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(238));

    t_yz_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(239));

    t_yz_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(240));

    t_yz_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(241));

    t_yz_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(242));

    t_yz_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(243));

    t_yz_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(244));

    t_yz_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(245));

    t_yz_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(246));

    t_yz_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(247));

    t_yz_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(248));

    t_yz_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(249));

    t_yz_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(250));

    t_yz_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(251));

    t_yz_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(252));

    t_yz_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(253));

    t_yz_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(254));

    t_yz_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(255));

    t_yz_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(256));

    t_yz_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(257));

    t_yz_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(258));

    t_yz_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(259));

    t_yz_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(260));

    t_yz_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(261));

    t_yz_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(262));

    t_yz_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(263));

    t_yz_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(264));

    t_yz_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(265));

    t_yz_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(266));

    t_yz_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(267));

    t_yz_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(268));

    t_yz_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(269));

    t_yz_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(270));

    t_yz_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(271));

    t_yz_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(272));

    t_yz_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(273));

    t_yz_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(274));

    t_yz_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(275));

    t_yz_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(276));

    t_yz_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(277));

    t_yz_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(278));

    t_yz_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(279));

    t_yz_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(280));

    t_yz_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(281));

    t_yz_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(282));

    t_yz_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(283));

    t_yz_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(284));

    t_yz_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(285));

    t_yz_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(286));

    t_yz_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(287));

    t_yz_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(288));

    t_yz_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(289));

    t_yz_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(290));

    t_yz_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(291));

    t_yz_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(292));

    t_yz_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(293));

    t_yz_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(294));

    t_yz_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(295));

    t_yz_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(296));

    t_yz_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(297));

    t_yz_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(298));

    t_yz_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(299));

    t_yz_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(300));

    t_yz_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(301));

    t_yz_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(302));

    t_yz_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(303));

    t_yz_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(304));

    t_yz_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(305));

    t_yz_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(306));

    t_yz_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(307));

    t_yz_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(308));

    t_yz_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(309));

    t_yz_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(310));

    t_yz_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(311));

    t_yz_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(312));

    t_yz_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(313));

    t_yz_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(314));

    t_yz_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(315));

    t_yz_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(316));

    t_yz_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(317));

    t_yz_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(318));

    t_yz_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(319));

    t_yz_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(320));

    t_yz_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(321));

    t_yz_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(322));

    t_yz_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(323));

    t_yz_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(324));

    t_yz_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(325));

    t_yz_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(326));

    t_yz_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(327));

    t_yz_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(328));

    t_yz_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(329));

    t_yz_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(330));

    t_yz_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(331));

    t_yz_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(332));

    t_yz_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(333));

    t_yz_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(334));

    t_yz_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(335));

    t_yz_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(336));

    t_yz_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(337));

    t_yz_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(338));

    t_yz_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(339));

    t_yz_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(340));

    t_yz_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(341));

    t_yz_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(342));

    t_yz_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(343));

    t_yz_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(344));

    t_yz_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(345));

    t_yz_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(346));

    t_yz_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(347));

    t_yz_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(348));

    t_yz_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(349));

    t_yz_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(350));

    t_yz_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(351));

    t_yz_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(352));

    t_yz_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(353));

    t_yz_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(354));

    t_yz_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(355));

    t_yz_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(356));

    t_yz_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(357));

    t_yz_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(358));

    t_yz_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(359));

    t_yz_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(360));

    t_yz_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(361));

    t_yz_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(362));

    t_yz_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(363));

    t_yz_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(364));

    t_yz_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(365));

    t_yz_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(366));

    t_yz_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(367));

    t_yz_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(368));

    t_yz_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(369));

    t_yz_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(370));

    t_yz_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(371));

    t_yz_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(372));

    t_yz_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(373));

    t_yz_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(374));

    t_yz_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(375));

    t_yz_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(376));

    t_yz_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(377));

    t_yz_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(378));

    t_yz_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(379));

    t_yz_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(380));

    t_yz_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(381));

    t_yz_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(382));

    t_yz_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(383));

    t_yz_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(384));

    t_yz_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(385));

    t_yz_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(386));

    t_yz_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(387));

    t_yz_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(388));

    t_yz_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(389));

    t_yz_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(390));

    t_yz_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(391));

    t_yz_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(392));

    t_yz_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(393));

    t_yz_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(394));

    t_yz_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(395));

    t_yz_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(396));

    t_yz_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(397));

    t_yz_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(398));

    t_yz_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(399));

    t_yz_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(400));

    t_yz_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(401));

    t_yz_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(402));

    t_yz_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(403));

    t_yz_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(404));

    t_yz_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(405));

    t_yz_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(406));

    t_yz_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(407));

    t_yz_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(408));

    t_yz_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(409));

    t_yz_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(410));

    t_yz_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(411));

    t_yz_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(412));

    t_yz_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(413));

    t_yz_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(414));

    t_yz_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(415));

    t_yz_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(416));

    t_yz_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(417));

    t_yz_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(418));

    t_yz_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(419));

    t_yz_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(420));

    t_yz_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(421));

    t_yz_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(422));

    t_yz_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(423));

    t_yz_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(424));

    t_yz_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(425));

    t_yz_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(426));

    t_yz_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(427));

    t_yz_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(428));

    t_yz_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(429));

    t_yz_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(430));

    t_yz_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(431));

    t_yy_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(432));

    t_yy_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(433));

    t_yy_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(434));

    t_yy_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(435));

    t_yy_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(436));

    t_yy_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(437));

    t_yy_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(438));

    t_yy_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(439));

    t_yy_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(440));

    t_yy_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(441));

    t_yy_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(442));

    t_yy_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(443));

    t_yy_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(444));

    t_yy_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(445));

    t_yy_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(446));

    t_yy_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(447));

    t_yy_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(448));

    t_yy_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(449));

    t_yy_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(450));

    t_yy_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(451));

    t_yy_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(452));

    t_yy_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(453));

    t_yy_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(454));

    t_yy_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(455));

    t_yy_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(456));

    t_yy_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(457));

    t_yy_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(458));

    t_yy_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(459));

    t_yy_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(460));

    t_yy_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(461));

    t_yy_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(462));

    t_yy_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(463));

    t_yy_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(464));

    t_yy_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(465));

    t_yy_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(466));

    t_yy_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(467));

    t_yy_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(468));

    t_yy_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(469));

    t_yy_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(470));

    t_yy_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(471));

    t_yy_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(472));

    t_yy_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(473));

    t_yy_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(474));

    t_yy_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(475));

    t_yy_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(476));

    t_yy_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(477));

    t_yy_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(478));

    t_yy_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(479));

    t_yy_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(480));

    t_yy_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(481));

    t_yy_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(482));

    t_yy_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(483));

    t_yy_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(484));

    t_yy_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(485));

    t_yy_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(486));

    t_yy_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(487));

    t_yy_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(488));

    t_yy_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(489));

    t_yy_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(490));

    t_yy_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(491));

    t_yy_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(492));

    t_yy_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(493));

    t_yy_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(494));

    t_yy_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(495));

    t_yy_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(496));

    t_yy_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(497));

    t_yy_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(498));

    t_yy_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(499));

    t_yy_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(500));

    t_yy_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(501));

    t_yy_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(502));

    t_yy_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(503));

    t_yy_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(504));

    t_yy_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(505));

    t_yy_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(506));

    t_yy_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(507));

    t_yy_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(508));

    t_yy_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(509));

    t_yy_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(510));

    t_yy_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(511));

    t_yy_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(512));

    t_yy_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(513));

    t_yy_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(514));

    t_yy_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(515));

    t_yy_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(516));

    t_yy_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(517));

    t_yy_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(518));

    t_yy_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(519));

    t_yy_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(520));

    t_yy_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(521));

    t_yy_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(522));

    t_yy_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(523));

    t_yy_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(524));

    t_yy_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(525));

    t_yy_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(526));

    t_yy_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(527));

    t_yy_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(528));

    t_yy_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(529));

    t_yy_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(530));

    t_yy_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(531));

    t_yy_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(532));

    t_yy_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(533));

    t_yy_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(534));

    t_yy_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(535));

    t_yy_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(536));

    t_yy_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(537));

    t_yy_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(538));

    t_yy_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(539));

    t_yy_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(540));

    t_yy_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(541));

    t_yy_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(542));

    t_yy_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(543));

    t_yy_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(544));

    t_yy_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(545));

    t_yy_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(546));

    t_yy_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(547));

    t_yy_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(548));

    t_yy_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(549));

    t_yy_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(550));

    t_yy_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(551));

    t_yy_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(552));

    t_yy_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(553));

    t_yy_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(554));

    t_yy_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(555));

    t_yy_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(556));

    t_yy_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(557));

    t_yy_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(558));

    t_yy_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(559));

    t_yy_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(560));

    t_yy_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(561));

    t_yy_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(562));

    t_yy_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(563));

    t_yy_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(564));

    t_yy_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(565));

    t_yy_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(566));

    t_yy_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(567));

    t_yy_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(568));

    t_yy_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(569));

    t_yy_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(570));

    t_yy_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(571));

    t_yy_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(572));

    t_yy_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(573));

    t_yy_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(574));

    t_yy_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(575));

    t_yy_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(576));

    t_yy_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(577));

    t_yy_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(578));

    t_yy_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(579));

    t_yy_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(580));

    t_yy_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(581));

    t_yy_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(582));

    t_yy_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(583));

    t_yy_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(584));

    t_yy_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(585));

    t_yy_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(586));

    t_yy_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(587));

    t_yy_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(588));

    t_yy_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(589));

    t_yy_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(590));

    t_yy_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(591));

    t_yy_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(592));

    t_yy_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(593));

    t_yy_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(594));

    t_yy_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(595));

    t_yy_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(596));

    t_yy_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(597));

    t_yy_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(598));

    t_yy_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(599));

    t_yy_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(600));

    t_yy_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(601));

    t_yy_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(602));

    t_yy_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(603));

    t_yy_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(604));

    t_yy_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(605));

    t_yy_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(606));

    t_yy_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(607));

    t_yy_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(608));

    t_yy_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(609));

    t_yy_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(610));

    t_yy_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(611));

    t_yy_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(612));

    t_yy_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(613));

    t_yy_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(614));

    t_yy_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(615));

    t_yy_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(616));

    t_yy_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(617));

    t_yy_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(618));

    t_yy_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(619));

    t_yy_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(620));

    t_yy_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(621));

    t_yy_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(622));

    t_yy_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(623));

    t_yy_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(624));

    t_yy_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(625));

    t_yy_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(626));

    t_yy_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(627));

    t_yy_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(628));

    t_yy_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(629));

    t_yy_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(630));

    t_yy_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(631));

    t_yy_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(632));

    t_yy_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(633));

    t_yy_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(634));

    t_yy_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(635));

    t_yy_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(636));

    t_yy_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(637));

    t_yy_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(638));

    t_yy_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(639));

    t_yy_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(640));

    t_yy_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(641));

    t_yy_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(642));

    t_yy_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(643));

    t_yy_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(644));

    t_yy_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(645));

    t_yy_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(646));

    t_yy_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(647));

    t_xz_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(648));

    t_xz_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(649));

    t_xz_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(650));

    t_xz_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(651));

    t_xz_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(652));

    t_xz_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(653));

    t_xz_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(654));

    t_xz_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(655));

    t_xz_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(656));

    t_xz_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(657));

    t_xz_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(658));

    t_xz_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(659));

    t_xz_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(660));

    t_xz_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(661));

    t_xz_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(662));

    t_xz_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(663));

    t_xz_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(664));

    t_xz_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(665));

    t_xz_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(666));

    t_xz_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(667));

    t_xz_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(668));

    t_xz_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(669));

    t_xz_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(670));

    t_xz_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(671));

    t_xz_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(672));

    t_xz_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(673));

    t_xz_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(674));

    t_xz_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(675));

    t_xz_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(676));

    t_xz_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(677));

    t_xz_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(678));

    t_xz_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(679));

    t_xz_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(680));

    t_xz_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(681));

    t_xz_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(682));

    t_xz_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(683));

    t_xz_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(684));

    t_xz_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(685));

    t_xz_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(686));

    t_xz_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(687));

    t_xz_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(688));

    t_xz_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(689));

    t_xz_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(690));

    t_xz_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(691));

    t_xz_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(692));

    t_xz_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(693));

    t_xz_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(694));

    t_xz_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(695));

    t_xz_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(696));

    t_xz_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(697));

    t_xz_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(698));

    t_xz_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(699));

    t_xz_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(700));

    t_xz_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(701));

    t_xz_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(702));

    t_xz_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(703));

    t_xz_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(704));

    t_xz_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(705));

    t_xz_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(706));

    t_xz_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(707));

    t_xz_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(708));

    t_xz_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(709));

    t_xz_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(710));

    t_xz_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(711));

    t_xz_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(712));

    t_xz_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(713));

    t_xz_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(714));

    t_xz_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(715));

    t_xz_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(716));

    t_xz_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(717));

    t_xz_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(718));

    t_xz_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(719));

    t_xz_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(720));

    t_xz_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(721));

    t_xz_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(722));

    t_xz_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(723));

    t_xz_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(724));

    t_xz_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(725));

    t_xz_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(726));

    t_xz_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(727));

    t_xz_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(728));

    t_xz_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(729));

    t_xz_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(730));

    t_xz_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(731));

    t_xz_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(732));

    t_xz_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(733));

    t_xz_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(734));

    t_xz_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(735));

    t_xz_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(736));

    t_xz_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(737));

    t_xz_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(738));

    t_xz_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(739));

    t_xz_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(740));

    t_xz_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(741));

    t_xz_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(742));

    t_xz_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(743));

    t_xz_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(744));

    t_xz_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(745));

    t_xz_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(746));

    t_xz_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(747));

    t_xz_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(748));

    t_xz_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(749));

    t_xz_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(750));

    t_xz_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(751));

    t_xz_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(752));

    t_xz_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(753));

    t_xz_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(754));

    t_xz_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(755));

    t_xz_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(756));

    t_xz_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(757));

    t_xz_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(758));

    t_xz_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(759));

    t_xz_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(760));

    t_xz_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(761));

    t_xz_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(762));

    t_xz_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(763));

    t_xz_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(764));

    t_xz_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(765));

    t_xz_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(766));

    t_xz_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(767));

    t_xz_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(768));

    t_xz_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(769));

    t_xz_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(770));

    t_xz_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(771));

    t_xz_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(772));

    t_xz_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(773));

    t_xz_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(774));

    t_xz_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(775));

    t_xz_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(776));

    t_xz_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(777));

    t_xz_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(778));

    t_xz_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(779));

    t_xz_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(780));

    t_xz_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(781));

    t_xz_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(782));

    t_xz_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(783));

    t_xz_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(784));

    t_xz_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(785));

    t_xz_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(786));

    t_xz_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(787));

    t_xz_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(788));

    t_xz_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(789));

    t_xz_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(790));

    t_xz_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(791));

    t_xz_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(792));

    t_xz_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(793));

    t_xz_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(794));

    t_xz_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(795));

    t_xz_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(796));

    t_xz_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(797));

    t_xz_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(798));

    t_xz_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(799));

    t_xz_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(800));

    t_xz_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(801));

    t_xz_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(802));

    t_xz_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(803));

    t_xz_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(804));

    t_xz_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(805));

    t_xz_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(806));

    t_xz_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(807));

    t_xz_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(808));

    t_xz_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(809));

    t_xz_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(810));

    t_xz_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(811));

    t_xz_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(812));

    t_xz_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(813));

    t_xz_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(814));

    t_xz_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(815));

    t_xz_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(816));

    t_xz_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(817));

    t_xz_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(818));

    t_xz_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(819));

    t_xz_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(820));

    t_xz_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(821));

    t_xz_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(822));

    t_xz_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(823));

    t_xz_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(824));

    t_xz_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(825));

    t_xz_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(826));

    t_xz_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(827));

    t_xz_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(828));

    t_xz_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(829));

    t_xz_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(830));

    t_xz_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(831));

    t_xz_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(832));

    t_xz_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(833));

    t_xz_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(834));

    t_xz_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(835));

    t_xz_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(836));

    t_xz_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(837));

    t_xz_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(838));

    t_xz_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(839));

    t_xz_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(840));

    t_xz_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(841));

    t_xz_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(842));

    t_xz_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(843));

    t_xz_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(844));

    t_xz_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(845));

    t_xz_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(846));

    t_xz_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(847));

    t_xz_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(848));

    t_xz_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(849));

    t_xz_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(850));

    t_xz_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(851));

    t_xz_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(852));

    t_xz_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(853));

    t_xz_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(854));

    t_xz_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(855));

    t_xz_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(856));

    t_xz_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(857));

    t_xz_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(858));

    t_xz_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(859));

    t_xz_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(860));

    t_xz_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(861));

    t_xz_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(862));

    t_xz_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(863));

    t_xy_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(864));

    t_xy_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(865));

    t_xy_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(866));

    t_xy_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(867));

    t_xy_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(868));

    t_xy_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(869));

    t_xy_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(870));

    t_xy_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(871));

    t_xy_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(872));

    t_xy_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(873));

    t_xy_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(874));

    t_xy_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(875));

    t_xy_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(876));

    t_xy_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(877));

    t_xy_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(878));

    t_xy_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(879));

    t_xy_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(880));

    t_xy_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(881));

    t_xy_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(882));

    t_xy_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(883));

    t_xy_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(884));

    t_xy_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(885));

    t_xy_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(886));

    t_xy_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(887));

    t_xy_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(888));

    t_xy_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(889));

    t_xy_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(890));

    t_xy_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(891));

    t_xy_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(892));

    t_xy_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(893));

    t_xy_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(894));

    t_xy_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(895));

    t_xy_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(896));

    t_xy_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(897));

    t_xy_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(898));

    t_xy_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(899));

    t_xy_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(900));

    t_xy_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(901));

    t_xy_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(902));

    t_xy_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(903));

    t_xy_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(904));

    t_xy_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(905));

    t_xy_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(906));

    t_xy_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(907));

    t_xy_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(908));

    t_xy_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(909));

    t_xy_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(910));

    t_xy_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(911));

    t_xy_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(912));

    t_xy_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(913));

    t_xy_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(914));

    t_xy_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(915));

    t_xy_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(916));

    t_xy_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(917));

    t_xy_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(918));

    t_xy_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(919));

    t_xy_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(920));

    t_xy_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(921));

    t_xy_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(922));

    t_xy_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(923));

    t_xy_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(924));

    t_xy_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(925));

    t_xy_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(926));

    t_xy_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(927));

    t_xy_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(928));

    t_xy_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(929));

    t_xy_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(930));

    t_xy_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(931));

    t_xy_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(932));

    t_xy_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(933));

    t_xy_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(934));

    t_xy_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(935));

    t_xy_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(936));

    t_xy_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(937));

    t_xy_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(938));

    t_xy_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(939));

    t_xy_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(940));

    t_xy_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(941));

    t_xy_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(942));

    t_xy_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(943));

    t_xy_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(944));

    t_xy_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(945));

    t_xy_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(946));

    t_xy_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(947));

    t_xy_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(948));

    t_xy_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(949));

    t_xy_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(950));

    t_xy_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(951));

    t_xy_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(952));

    t_xy_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(953));

    t_xy_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(954));

    t_xy_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(955));

    t_xy_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(956));

    t_xy_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(957));

    t_xy_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(958));

    t_xy_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(959));

    t_xy_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(960));

    t_xy_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(961));

    t_xy_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(962));

    t_xy_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(963));

    t_xy_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(964));

    t_xy_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(965));

    t_xy_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(966));

    t_xy_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(967));

    t_xy_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(968));

    t_xy_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(969));

    t_xy_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(970));

    t_xy_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(971));

    t_xy_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(972));

    t_xy_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(973));

    t_xy_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(974));

    t_xy_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(975));

    t_xy_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(976));

    t_xy_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(977));

    t_xy_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(978));

    t_xy_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(979));

    t_xy_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(980));

    t_xy_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(981));

    t_xy_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(982));

    t_xy_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(983));

    t_xy_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(984));

    t_xy_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(985));

    t_xy_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(986));

    t_xy_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(987));

    t_xy_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(988));

    t_xy_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(989));

    t_xy_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(990));

    t_xy_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(991));

    t_xy_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(992));

    t_xy_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(993));

    t_xy_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(994));

    t_xy_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(995));

    t_xy_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(996));

    t_xy_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(997));

    t_xy_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(998));

    t_xy_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(999));

    t_xy_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1000));

    t_xy_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1001));

    t_xy_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1002));

    t_xy_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1003));

    t_xy_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1004));

    t_xy_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1005));

    t_xy_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1006));

    t_xy_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1007));

    t_xy_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1008));

    t_xy_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1009));

    t_xy_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1010));

    t_xy_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1011));

    t_xy_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1012));

    t_xy_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1013));

    t_xy_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1014));

    t_xy_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1015));

    t_xy_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1016));

    t_xy_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1017));

    t_xy_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1018));

    t_xy_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1019));

    t_xy_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1020));

    t_xy_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1021));

    t_xy_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1022));

    t_xy_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1023));

    t_xy_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1024));

    t_xy_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1025));

    t_xy_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1026));

    t_xy_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1027));

    t_xy_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1028));

    t_xy_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1029));

    t_xy_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1030));

    t_xy_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1031));

    t_xy_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1032));

    t_xy_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1033));

    t_xy_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1034));

    t_xy_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1035));

    t_xy_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1036));

    t_xy_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1037));

    t_xy_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1038));

    t_xy_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1039));

    t_xy_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1040));

    t_xy_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1041));

    t_xy_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1042));

    t_xy_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1043));

    t_xy_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1044));

    t_xy_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1045));

    t_xy_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1046));

    t_xy_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1047));

    t_xy_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1048));

    t_xy_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1049));

    t_xy_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1050));

    t_xy_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1051));

    t_xy_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1052));

    t_xy_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1053));

    t_xy_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1054));

    t_xy_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1055));

    t_xy_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1056));

    t_xy_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1057));

    t_xy_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1058));

    t_xy_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1059));

    t_xy_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1060));

    t_xy_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1061));

    t_xy_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1062));

    t_xy_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1063));

    t_xy_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1064));

    t_xy_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1065));

    t_xy_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1066));

    t_xy_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1067));

    t_xy_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1068));

    t_xy_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1069));

    t_xy_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1070));

    t_xy_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1071));

    t_xy_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1072));

    t_xy_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1073));

    t_xy_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1074));

    t_xy_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1075));

    t_xy_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1076));

    t_xy_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1077));

    t_xy_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1078));

    t_xy_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1079));

    t_xx_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1080));

    t_xx_zz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1081));

    t_xx_zz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1082));

    t_xx_zz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1083));

    t_xx_zz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1084));

    t_xx_zz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1085));

    t_xx_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1086));

    t_xx_zz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1087));

    t_xx_zz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1088));

    t_xx_zz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1089));

    t_xx_zz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1090));

    t_xx_zz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1091));

    t_xx_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1092));

    t_xx_zz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1093));

    t_xx_zz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1094));

    t_xx_zz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1095));

    t_xx_zz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1096));

    t_xx_zz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1097));

    t_xx_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1098));

    t_xx_zz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1099));

    t_xx_zz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1100));

    t_xx_zz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1101));

    t_xx_zz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1102));

    t_xx_zz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1103));

    t_xx_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1104));

    t_xx_zz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1105));

    t_xx_zz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1106));

    t_xx_zz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1107));

    t_xx_zz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1108));

    t_xx_zz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1109));

    t_xx_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1110));

    t_xx_zz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1111));

    t_xx_zz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1112));

    t_xx_zz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1113));

    t_xx_zz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1114));

    t_xx_zz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1115));

    t_xx_yz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1116));

    t_xx_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1117));

    t_xx_yz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1118));

    t_xx_yz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1119));

    t_xx_yz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1120));

    t_xx_yz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1121));

    t_xx_yz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1122));

    t_xx_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1123));

    t_xx_yz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1124));

    t_xx_yz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1125));

    t_xx_yz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1126));

    t_xx_yz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1127));

    t_xx_yz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1128));

    t_xx_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1129));

    t_xx_yz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1130));

    t_xx_yz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1131));

    t_xx_yz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1132));

    t_xx_yz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1133));

    t_xx_yz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1134));

    t_xx_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1135));

    t_xx_yz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1136));

    t_xx_yz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1137));

    t_xx_yz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1138));

    t_xx_yz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1139));

    t_xx_yz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1140));

    t_xx_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1141));

    t_xx_yz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1142));

    t_xx_yz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1143));

    t_xx_yz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1144));

    t_xx_yz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1145));

    t_xx_yz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1146));

    t_xx_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1147));

    t_xx_yz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1148));

    t_xx_yz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1149));

    t_xx_yz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1150));

    t_xx_yz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1151));

    t_xx_yy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1152));

    t_xx_yy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1153));

    t_xx_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1154));

    t_xx_yy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1155));

    t_xx_yy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1156));

    t_xx_yy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1157));

    t_xx_yy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1158));

    t_xx_yy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1159));

    t_xx_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1160));

    t_xx_yy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1161));

    t_xx_yy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1162));

    t_xx_yy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1163));

    t_xx_yy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1164));

    t_xx_yy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1165));

    t_xx_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1166));

    t_xx_yy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1167));

    t_xx_yy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1168));

    t_xx_yy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1169));

    t_xx_yy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1170));

    t_xx_yy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1171));

    t_xx_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1172));

    t_xx_yy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1173));

    t_xx_yy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1174));

    t_xx_yy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1175));

    t_xx_yy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1176));

    t_xx_yy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1177));

    t_xx_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1178));

    t_xx_yy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1179));

    t_xx_yy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1180));

    t_xx_yy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1181));

    t_xx_yy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1182));

    t_xx_yy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1183));

    t_xx_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1184));

    t_xx_yy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1185));

    t_xx_yy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1186));

    t_xx_yy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1187));

    t_xx_xz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1188));

    t_xx_xz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1189));

    t_xx_xz_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1190));

    t_xx_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1191));

    t_xx_xz_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1192));

    t_xx_xz_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1193));

    t_xx_xz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1194));

    t_xx_xz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1195));

    t_xx_xz_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1196));

    t_xx_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1197));

    t_xx_xz_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1198));

    t_xx_xz_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1199));

    t_xx_xz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1200));

    t_xx_xz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1201));

    t_xx_xz_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1202));

    t_xx_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1203));

    t_xx_xz_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1204));

    t_xx_xz_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1205));

    t_xx_xz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1206));

    t_xx_xz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1207));

    t_xx_xz_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1208));

    t_xx_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1209));

    t_xx_xz_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1210));

    t_xx_xz_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1211));

    t_xx_xz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1212));

    t_xx_xz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1213));

    t_xx_xz_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1214));

    t_xx_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1215));

    t_xx_xz_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1216));

    t_xx_xz_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1217));

    t_xx_xz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1218));

    t_xx_xz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1219));

    t_xx_xz_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1220));

    t_xx_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1221));

    t_xx_xz_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1222));

    t_xx_xz_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1223));

    t_xx_xy_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1224));

    t_xx_xy_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1225));

    t_xx_xy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1226));

    t_xx_xy_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1227));

    t_xx_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1228));

    t_xx_xy_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1229));

    t_xx_xy_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1230));

    t_xx_xy_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1231));

    t_xx_xy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1232));

    t_xx_xy_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1233));

    t_xx_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1234));

    t_xx_xy_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1235));

    t_xx_xy_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1236));

    t_xx_xy_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1237));

    t_xx_xy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1238));

    t_xx_xy_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1239));

    t_xx_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1240));

    t_xx_xy_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1241));

    t_xx_xy_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1242));

    t_xx_xy_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1243));

    t_xx_xy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1244));

    t_xx_xy_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1245));

    t_xx_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1246));

    t_xx_xy_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1247));

    t_xx_xy_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1248));

    t_xx_xy_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1249));

    t_xx_xy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1250));

    t_xx_xy_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1251));

    t_xx_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1252));

    t_xx_xy_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1253));

    t_xx_xy_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1254));

    t_xx_xy_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1255));

    t_xx_xy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1256));

    t_xx_xy_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1257));

    t_xx_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1258));

    t_xx_xy_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1259));

    t_xx_xx_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(1260));

    t_xx_xx_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1261));

    t_xx_xx_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(1262));

    t_xx_xx_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(1263));

    t_xx_xx_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(1264));

    t_xx_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(1265));

    t_xx_xx_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(1266));

    t_xx_xx_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(1267));

    t_xx_xx_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(1268));

    t_xx_xx_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(1269));

    t_xx_xx_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(1270));

    t_xx_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(1271));

    t_xx_xx_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(1272));

    t_xx_xx_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(1273));

    t_xx_xx_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(1274));

    t_xx_xx_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(1275));

    t_xx_xx_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(1276));

    t_xx_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(1277));

    t_xx_xx_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(1278));

    t_xx_xx_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(1279));

    t_xx_xx_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(1280));

    t_xx_xx_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(1281));

    t_xx_xx_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(1282));

    t_xx_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(1283));

    t_xx_xx_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(1284));

    t_xx_xx_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(1285));

    t_xx_xx_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(1286));

    t_xx_xx_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(1287));

    t_xx_xx_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(1288));

    t_xx_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(1289));

    t_xx_xx_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(1290));

    t_xx_xx_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(1291));

    t_xx_xx_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(1292));

    t_xx_xx_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(1293));

    t_xx_xx_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(1294));

    t_xx_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(1295));

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

    #pragma omp simd align(rab_z, t_z_zz_xx_xx, t_z_zz_xx_xy, t_z_zz_xx_xz, t_z_zz_xx_yy,\
                           t_z_zz_xx_yz, t_z_zz_xx_zz, t_z_zz_xy_xx, t_z_zz_xy_xy, t_z_zz_xy_xz,\
                           t_z_zz_xy_yy, t_z_zz_xy_yz, t_z_zz_xy_zz, t_z_zz_xz_xx, t_z_zz_xz_xy,\
                           t_z_zz_xz_xz, t_z_zz_xz_yy, t_z_zz_xz_yz, t_z_zz_xz_zz, t_z_zz_yy_xx,\
                           t_z_zz_yy_xy, t_z_zz_yy_xz, t_z_zz_yy_yy, t_z_zz_yy_yz, t_z_zz_yy_zz,\
                           t_z_zz_yz_xx, t_z_zz_yz_xy, t_z_zz_yz_xz, t_z_zz_yz_yy, t_z_zz_yz_yz,\
                           t_z_zz_yz_zz, t_z_zz_zz_xx, t_z_zz_zz_xy, t_z_zz_zz_xz, t_z_zz_zz_yy,\
                           t_z_zz_zz_yz, t_z_zz_zz_zz, t_z_zzz_xx_xx, t_z_zzz_xx_xy,\
                           t_z_zzz_xx_xz, t_z_zzz_xx_yy, t_z_zzz_xx_yz, t_z_zzz_xx_zz,\
                           t_z_zzz_xy_xx, t_z_zzz_xy_xy, t_z_zzz_xy_xz, t_z_zzz_xy_yy,\
                           t_z_zzz_xy_yz, t_z_zzz_xy_zz, t_z_zzz_xz_xx, t_z_zzz_xz_xy,\
                           t_z_zzz_xz_xz, t_z_zzz_xz_yy, t_z_zzz_xz_yz, t_z_zzz_xz_zz,\
                           t_z_zzz_yy_xx, t_z_zzz_yy_xy, t_z_zzz_yy_xz, t_z_zzz_yy_yy,\
                           t_z_zzz_yy_yz, t_z_zzz_yy_zz, t_z_zzz_yz_xx, t_z_zzz_yz_xy,\
                           t_z_zzz_yz_xz, t_z_zzz_yz_yy, t_z_zzz_yz_yz, t_z_zzz_yz_zz,\
                           t_z_zzz_zz_xx, t_z_zzz_zz_xy, t_z_zzz_zz_xz, t_z_zzz_zz_yy,\
                           t_z_zzz_zz_yz, t_z_zzz_zz_zz, t_zz_zz_xx_xx, t_zz_zz_xx_xy,\
                           t_zz_zz_xx_xz, t_zz_zz_xx_yy, t_zz_zz_xx_yz, t_zz_zz_xx_zz,\
                           t_zz_zz_xy_xx, t_zz_zz_xy_xy, t_zz_zz_xy_xz, t_zz_zz_xy_yy,\
                           t_zz_zz_xy_yz, t_zz_zz_xy_zz, t_zz_zz_xz_xx, t_zz_zz_xz_xy,\
                           t_zz_zz_xz_xz, t_zz_zz_xz_yy, t_zz_zz_xz_yz, t_zz_zz_xz_zz,\
                           t_zz_zz_yy_xx, t_zz_zz_yy_xy, t_zz_zz_yy_xz, t_zz_zz_yy_yy,\
                           t_zz_zz_yy_yz, t_zz_zz_yy_zz, t_zz_zz_yz_xx, t_zz_zz_yz_xy,\
                           t_zz_zz_yz_xz, t_zz_zz_yz_yy, t_zz_zz_yz_yz, t_zz_zz_yz_zz,\
                           t_zz_zz_zz_xx, t_zz_zz_zz_xy, t_zz_zz_zz_xz, t_zz_zz_zz_yy,\
                           t_zz_zz_zz_yz, t_zz_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_zz_zz[i] = t_z_zzz_zz_zz[i] - rab_z[i] * t_z_zz_zz_zz[i];

        t_zz_zz_zz_yz[i] = t_z_zzz_zz_yz[i] - rab_z[i] * t_z_zz_zz_yz[i];

        t_zz_zz_zz_yy[i] = t_z_zzz_zz_yy[i] - rab_z[i] * t_z_zz_zz_yy[i];

        t_zz_zz_zz_xz[i] = t_z_zzz_zz_xz[i] - rab_z[i] * t_z_zz_zz_xz[i];

        t_zz_zz_zz_xy[i] = t_z_zzz_zz_xy[i] - rab_z[i] * t_z_zz_zz_xy[i];

        t_zz_zz_zz_xx[i] = t_z_zzz_zz_xx[i] - rab_z[i] * t_z_zz_zz_xx[i];

        t_zz_zz_yz_zz[i] = t_z_zzz_yz_zz[i] - rab_z[i] * t_z_zz_yz_zz[i];

        t_zz_zz_yz_yz[i] = t_z_zzz_yz_yz[i] - rab_z[i] * t_z_zz_yz_yz[i];

        t_zz_zz_yz_yy[i] = t_z_zzz_yz_yy[i] - rab_z[i] * t_z_zz_yz_yy[i];

        t_zz_zz_yz_xz[i] = t_z_zzz_yz_xz[i] - rab_z[i] * t_z_zz_yz_xz[i];

        t_zz_zz_yz_xy[i] = t_z_zzz_yz_xy[i] - rab_z[i] * t_z_zz_yz_xy[i];

        t_zz_zz_yz_xx[i] = t_z_zzz_yz_xx[i] - rab_z[i] * t_z_zz_yz_xx[i];

        t_zz_zz_yy_zz[i] = t_z_zzz_yy_zz[i] - rab_z[i] * t_z_zz_yy_zz[i];

        t_zz_zz_yy_yz[i] = t_z_zzz_yy_yz[i] - rab_z[i] * t_z_zz_yy_yz[i];

        t_zz_zz_yy_yy[i] = t_z_zzz_yy_yy[i] - rab_z[i] * t_z_zz_yy_yy[i];

        t_zz_zz_yy_xz[i] = t_z_zzz_yy_xz[i] - rab_z[i] * t_z_zz_yy_xz[i];

        t_zz_zz_yy_xy[i] = t_z_zzz_yy_xy[i] - rab_z[i] * t_z_zz_yy_xy[i];

        t_zz_zz_yy_xx[i] = t_z_zzz_yy_xx[i] - rab_z[i] * t_z_zz_yy_xx[i];

        t_zz_zz_xz_zz[i] = t_z_zzz_xz_zz[i] - rab_z[i] * t_z_zz_xz_zz[i];

        t_zz_zz_xz_yz[i] = t_z_zzz_xz_yz[i] - rab_z[i] * t_z_zz_xz_yz[i];

        t_zz_zz_xz_yy[i] = t_z_zzz_xz_yy[i] - rab_z[i] * t_z_zz_xz_yy[i];

        t_zz_zz_xz_xz[i] = t_z_zzz_xz_xz[i] - rab_z[i] * t_z_zz_xz_xz[i];

        t_zz_zz_xz_xy[i] = t_z_zzz_xz_xy[i] - rab_z[i] * t_z_zz_xz_xy[i];

        t_zz_zz_xz_xx[i] = t_z_zzz_xz_xx[i] - rab_z[i] * t_z_zz_xz_xx[i];

        t_zz_zz_xy_zz[i] = t_z_zzz_xy_zz[i] - rab_z[i] * t_z_zz_xy_zz[i];

        t_zz_zz_xy_yz[i] = t_z_zzz_xy_yz[i] - rab_z[i] * t_z_zz_xy_yz[i];

        t_zz_zz_xy_yy[i] = t_z_zzz_xy_yy[i] - rab_z[i] * t_z_zz_xy_yy[i];

        t_zz_zz_xy_xz[i] = t_z_zzz_xy_xz[i] - rab_z[i] * t_z_zz_xy_xz[i];

        t_zz_zz_xy_xy[i] = t_z_zzz_xy_xy[i] - rab_z[i] * t_z_zz_xy_xy[i];

        t_zz_zz_xy_xx[i] = t_z_zzz_xy_xx[i] - rab_z[i] * t_z_zz_xy_xx[i];

        t_zz_zz_xx_zz[i] = t_z_zzz_xx_zz[i] - rab_z[i] * t_z_zz_xx_zz[i];

        t_zz_zz_xx_yz[i] = t_z_zzz_xx_yz[i] - rab_z[i] * t_z_zz_xx_yz[i];

        t_zz_zz_xx_yy[i] = t_z_zzz_xx_yy[i] - rab_z[i] * t_z_zz_xx_yy[i];

        t_zz_zz_xx_xz[i] = t_z_zzz_xx_xz[i] - rab_z[i] * t_z_zz_xx_xz[i];

        t_zz_zz_xx_xy[i] = t_z_zzz_xx_xy[i] - rab_z[i] * t_z_zz_xx_xy[i];

        t_zz_zz_xx_xx[i] = t_z_zzz_xx_xx[i] - rab_z[i] * t_z_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_yz_xx_xx, t_z_yz_xx_xy, t_z_yz_xx_xz, t_z_yz_xx_yy,\
                           t_z_yz_xx_yz, t_z_yz_xx_zz, t_z_yz_xy_xx, t_z_yz_xy_xy, t_z_yz_xy_xz,\
                           t_z_yz_xy_yy, t_z_yz_xy_yz, t_z_yz_xy_zz, t_z_yz_xz_xx, t_z_yz_xz_xy,\
                           t_z_yz_xz_xz, t_z_yz_xz_yy, t_z_yz_xz_yz, t_z_yz_xz_zz, t_z_yz_yy_xx,\
                           t_z_yz_yy_xy, t_z_yz_yy_xz, t_z_yz_yy_yy, t_z_yz_yy_yz, t_z_yz_yy_zz,\
                           t_z_yz_yz_xx, t_z_yz_yz_xy, t_z_yz_yz_xz, t_z_yz_yz_yy, t_z_yz_yz_yz,\
                           t_z_yz_yz_zz, t_z_yz_zz_xx, t_z_yz_zz_xy, t_z_yz_zz_xz, t_z_yz_zz_yy,\
                           t_z_yz_zz_yz, t_z_yz_zz_zz, t_z_yzz_xx_xx, t_z_yzz_xx_xy,\
                           t_z_yzz_xx_xz, t_z_yzz_xx_yy, t_z_yzz_xx_yz, t_z_yzz_xx_zz,\
                           t_z_yzz_xy_xx, t_z_yzz_xy_xy, t_z_yzz_xy_xz, t_z_yzz_xy_yy,\
                           t_z_yzz_xy_yz, t_z_yzz_xy_zz, t_z_yzz_xz_xx, t_z_yzz_xz_xy,\
                           t_z_yzz_xz_xz, t_z_yzz_xz_yy, t_z_yzz_xz_yz, t_z_yzz_xz_zz,\
                           t_z_yzz_yy_xx, t_z_yzz_yy_xy, t_z_yzz_yy_xz, t_z_yzz_yy_yy,\
                           t_z_yzz_yy_yz, t_z_yzz_yy_zz, t_z_yzz_yz_xx, t_z_yzz_yz_xy,\
                           t_z_yzz_yz_xz, t_z_yzz_yz_yy, t_z_yzz_yz_yz, t_z_yzz_yz_zz,\
                           t_z_yzz_zz_xx, t_z_yzz_zz_xy, t_z_yzz_zz_xz, t_z_yzz_zz_yy,\
                           t_z_yzz_zz_yz, t_z_yzz_zz_zz, t_zz_yz_xx_xx, t_zz_yz_xx_xy,\
                           t_zz_yz_xx_xz, t_zz_yz_xx_yy, t_zz_yz_xx_yz, t_zz_yz_xx_zz,\
                           t_zz_yz_xy_xx, t_zz_yz_xy_xy, t_zz_yz_xy_xz, t_zz_yz_xy_yy,\
                           t_zz_yz_xy_yz, t_zz_yz_xy_zz, t_zz_yz_xz_xx, t_zz_yz_xz_xy,\
                           t_zz_yz_xz_xz, t_zz_yz_xz_yy, t_zz_yz_xz_yz, t_zz_yz_xz_zz,\
                           t_zz_yz_yy_xx, t_zz_yz_yy_xy, t_zz_yz_yy_xz, t_zz_yz_yy_yy,\
                           t_zz_yz_yy_yz, t_zz_yz_yy_zz, t_zz_yz_yz_xx, t_zz_yz_yz_xy,\
                           t_zz_yz_yz_xz, t_zz_yz_yz_yy, t_zz_yz_yz_yz, t_zz_yz_yz_zz,\
                           t_zz_yz_zz_xx, t_zz_yz_zz_xy, t_zz_yz_zz_xz, t_zz_yz_zz_yy,\
                           t_zz_yz_zz_yz, t_zz_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_yz_zz_zz[i] = t_z_yzz_zz_zz[i] - rab_z[i] * t_z_yz_zz_zz[i];

        t_zz_yz_zz_yz[i] = t_z_yzz_zz_yz[i] - rab_z[i] * t_z_yz_zz_yz[i];

        t_zz_yz_zz_yy[i] = t_z_yzz_zz_yy[i] - rab_z[i] * t_z_yz_zz_yy[i];

        t_zz_yz_zz_xz[i] = t_z_yzz_zz_xz[i] - rab_z[i] * t_z_yz_zz_xz[i];

        t_zz_yz_zz_xy[i] = t_z_yzz_zz_xy[i] - rab_z[i] * t_z_yz_zz_xy[i];

        t_zz_yz_zz_xx[i] = t_z_yzz_zz_xx[i] - rab_z[i] * t_z_yz_zz_xx[i];

        t_zz_yz_yz_zz[i] = t_z_yzz_yz_zz[i] - rab_z[i] * t_z_yz_yz_zz[i];

        t_zz_yz_yz_yz[i] = t_z_yzz_yz_yz[i] - rab_z[i] * t_z_yz_yz_yz[i];

        t_zz_yz_yz_yy[i] = t_z_yzz_yz_yy[i] - rab_z[i] * t_z_yz_yz_yy[i];

        t_zz_yz_yz_xz[i] = t_z_yzz_yz_xz[i] - rab_z[i] * t_z_yz_yz_xz[i];

        t_zz_yz_yz_xy[i] = t_z_yzz_yz_xy[i] - rab_z[i] * t_z_yz_yz_xy[i];

        t_zz_yz_yz_xx[i] = t_z_yzz_yz_xx[i] - rab_z[i] * t_z_yz_yz_xx[i];

        t_zz_yz_yy_zz[i] = t_z_yzz_yy_zz[i] - rab_z[i] * t_z_yz_yy_zz[i];

        t_zz_yz_yy_yz[i] = t_z_yzz_yy_yz[i] - rab_z[i] * t_z_yz_yy_yz[i];

        t_zz_yz_yy_yy[i] = t_z_yzz_yy_yy[i] - rab_z[i] * t_z_yz_yy_yy[i];

        t_zz_yz_yy_xz[i] = t_z_yzz_yy_xz[i] - rab_z[i] * t_z_yz_yy_xz[i];

        t_zz_yz_yy_xy[i] = t_z_yzz_yy_xy[i] - rab_z[i] * t_z_yz_yy_xy[i];

        t_zz_yz_yy_xx[i] = t_z_yzz_yy_xx[i] - rab_z[i] * t_z_yz_yy_xx[i];

        t_zz_yz_xz_zz[i] = t_z_yzz_xz_zz[i] - rab_z[i] * t_z_yz_xz_zz[i];

        t_zz_yz_xz_yz[i] = t_z_yzz_xz_yz[i] - rab_z[i] * t_z_yz_xz_yz[i];

        t_zz_yz_xz_yy[i] = t_z_yzz_xz_yy[i] - rab_z[i] * t_z_yz_xz_yy[i];

        t_zz_yz_xz_xz[i] = t_z_yzz_xz_xz[i] - rab_z[i] * t_z_yz_xz_xz[i];

        t_zz_yz_xz_xy[i] = t_z_yzz_xz_xy[i] - rab_z[i] * t_z_yz_xz_xy[i];

        t_zz_yz_xz_xx[i] = t_z_yzz_xz_xx[i] - rab_z[i] * t_z_yz_xz_xx[i];

        t_zz_yz_xy_zz[i] = t_z_yzz_xy_zz[i] - rab_z[i] * t_z_yz_xy_zz[i];

        t_zz_yz_xy_yz[i] = t_z_yzz_xy_yz[i] - rab_z[i] * t_z_yz_xy_yz[i];

        t_zz_yz_xy_yy[i] = t_z_yzz_xy_yy[i] - rab_z[i] * t_z_yz_xy_yy[i];

        t_zz_yz_xy_xz[i] = t_z_yzz_xy_xz[i] - rab_z[i] * t_z_yz_xy_xz[i];

        t_zz_yz_xy_xy[i] = t_z_yzz_xy_xy[i] - rab_z[i] * t_z_yz_xy_xy[i];

        t_zz_yz_xy_xx[i] = t_z_yzz_xy_xx[i] - rab_z[i] * t_z_yz_xy_xx[i];

        t_zz_yz_xx_zz[i] = t_z_yzz_xx_zz[i] - rab_z[i] * t_z_yz_xx_zz[i];

        t_zz_yz_xx_yz[i] = t_z_yzz_xx_yz[i] - rab_z[i] * t_z_yz_xx_yz[i];

        t_zz_yz_xx_yy[i] = t_z_yzz_xx_yy[i] - rab_z[i] * t_z_yz_xx_yy[i];

        t_zz_yz_xx_xz[i] = t_z_yzz_xx_xz[i] - rab_z[i] * t_z_yz_xx_xz[i];

        t_zz_yz_xx_xy[i] = t_z_yzz_xx_xy[i] - rab_z[i] * t_z_yz_xx_xy[i];

        t_zz_yz_xx_xx[i] = t_z_yzz_xx_xx[i] - rab_z[i] * t_z_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_yy_xx_xx, t_z_yy_xx_xy, t_z_yy_xx_xz, t_z_yy_xx_yy,\
                           t_z_yy_xx_yz, t_z_yy_xx_zz, t_z_yy_xy_xx, t_z_yy_xy_xy, t_z_yy_xy_xz,\
                           t_z_yy_xy_yy, t_z_yy_xy_yz, t_z_yy_xy_zz, t_z_yy_xz_xx, t_z_yy_xz_xy,\
                           t_z_yy_xz_xz, t_z_yy_xz_yy, t_z_yy_xz_yz, t_z_yy_xz_zz, t_z_yy_yy_xx,\
                           t_z_yy_yy_xy, t_z_yy_yy_xz, t_z_yy_yy_yy, t_z_yy_yy_yz, t_z_yy_yy_zz,\
                           t_z_yy_yz_xx, t_z_yy_yz_xy, t_z_yy_yz_xz, t_z_yy_yz_yy, t_z_yy_yz_yz,\
                           t_z_yy_yz_zz, t_z_yy_zz_xx, t_z_yy_zz_xy, t_z_yy_zz_xz, t_z_yy_zz_yy,\
                           t_z_yy_zz_yz, t_z_yy_zz_zz, t_z_yyz_xx_xx, t_z_yyz_xx_xy,\
                           t_z_yyz_xx_xz, t_z_yyz_xx_yy, t_z_yyz_xx_yz, t_z_yyz_xx_zz,\
                           t_z_yyz_xy_xx, t_z_yyz_xy_xy, t_z_yyz_xy_xz, t_z_yyz_xy_yy,\
                           t_z_yyz_xy_yz, t_z_yyz_xy_zz, t_z_yyz_xz_xx, t_z_yyz_xz_xy,\
                           t_z_yyz_xz_xz, t_z_yyz_xz_yy, t_z_yyz_xz_yz, t_z_yyz_xz_zz,\
                           t_z_yyz_yy_xx, t_z_yyz_yy_xy, t_z_yyz_yy_xz, t_z_yyz_yy_yy,\
                           t_z_yyz_yy_yz, t_z_yyz_yy_zz, t_z_yyz_yz_xx, t_z_yyz_yz_xy,\
                           t_z_yyz_yz_xz, t_z_yyz_yz_yy, t_z_yyz_yz_yz, t_z_yyz_yz_zz,\
                           t_z_yyz_zz_xx, t_z_yyz_zz_xy, t_z_yyz_zz_xz, t_z_yyz_zz_yy,\
                           t_z_yyz_zz_yz, t_z_yyz_zz_zz, t_zz_yy_xx_xx, t_zz_yy_xx_xy,\
                           t_zz_yy_xx_xz, t_zz_yy_xx_yy, t_zz_yy_xx_yz, t_zz_yy_xx_zz,\
                           t_zz_yy_xy_xx, t_zz_yy_xy_xy, t_zz_yy_xy_xz, t_zz_yy_xy_yy,\
                           t_zz_yy_xy_yz, t_zz_yy_xy_zz, t_zz_yy_xz_xx, t_zz_yy_xz_xy,\
                           t_zz_yy_xz_xz, t_zz_yy_xz_yy, t_zz_yy_xz_yz, t_zz_yy_xz_zz,\
                           t_zz_yy_yy_xx, t_zz_yy_yy_xy, t_zz_yy_yy_xz, t_zz_yy_yy_yy,\
                           t_zz_yy_yy_yz, t_zz_yy_yy_zz, t_zz_yy_yz_xx, t_zz_yy_yz_xy,\
                           t_zz_yy_yz_xz, t_zz_yy_yz_yy, t_zz_yy_yz_yz, t_zz_yy_yz_zz,\
                           t_zz_yy_zz_xx, t_zz_yy_zz_xy, t_zz_yy_zz_xz, t_zz_yy_zz_yy,\
                           t_zz_yy_zz_yz, t_zz_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_yy_zz_zz[i] = t_z_yyz_zz_zz[i] - rab_z[i] * t_z_yy_zz_zz[i];

        t_zz_yy_zz_yz[i] = t_z_yyz_zz_yz[i] - rab_z[i] * t_z_yy_zz_yz[i];

        t_zz_yy_zz_yy[i] = t_z_yyz_zz_yy[i] - rab_z[i] * t_z_yy_zz_yy[i];

        t_zz_yy_zz_xz[i] = t_z_yyz_zz_xz[i] - rab_z[i] * t_z_yy_zz_xz[i];

        t_zz_yy_zz_xy[i] = t_z_yyz_zz_xy[i] - rab_z[i] * t_z_yy_zz_xy[i];

        t_zz_yy_zz_xx[i] = t_z_yyz_zz_xx[i] - rab_z[i] * t_z_yy_zz_xx[i];

        t_zz_yy_yz_zz[i] = t_z_yyz_yz_zz[i] - rab_z[i] * t_z_yy_yz_zz[i];

        t_zz_yy_yz_yz[i] = t_z_yyz_yz_yz[i] - rab_z[i] * t_z_yy_yz_yz[i];

        t_zz_yy_yz_yy[i] = t_z_yyz_yz_yy[i] - rab_z[i] * t_z_yy_yz_yy[i];

        t_zz_yy_yz_xz[i] = t_z_yyz_yz_xz[i] - rab_z[i] * t_z_yy_yz_xz[i];

        t_zz_yy_yz_xy[i] = t_z_yyz_yz_xy[i] - rab_z[i] * t_z_yy_yz_xy[i];

        t_zz_yy_yz_xx[i] = t_z_yyz_yz_xx[i] - rab_z[i] * t_z_yy_yz_xx[i];

        t_zz_yy_yy_zz[i] = t_z_yyz_yy_zz[i] - rab_z[i] * t_z_yy_yy_zz[i];

        t_zz_yy_yy_yz[i] = t_z_yyz_yy_yz[i] - rab_z[i] * t_z_yy_yy_yz[i];

        t_zz_yy_yy_yy[i] = t_z_yyz_yy_yy[i] - rab_z[i] * t_z_yy_yy_yy[i];

        t_zz_yy_yy_xz[i] = t_z_yyz_yy_xz[i] - rab_z[i] * t_z_yy_yy_xz[i];

        t_zz_yy_yy_xy[i] = t_z_yyz_yy_xy[i] - rab_z[i] * t_z_yy_yy_xy[i];

        t_zz_yy_yy_xx[i] = t_z_yyz_yy_xx[i] - rab_z[i] * t_z_yy_yy_xx[i];

        t_zz_yy_xz_zz[i] = t_z_yyz_xz_zz[i] - rab_z[i] * t_z_yy_xz_zz[i];

        t_zz_yy_xz_yz[i] = t_z_yyz_xz_yz[i] - rab_z[i] * t_z_yy_xz_yz[i];

        t_zz_yy_xz_yy[i] = t_z_yyz_xz_yy[i] - rab_z[i] * t_z_yy_xz_yy[i];

        t_zz_yy_xz_xz[i] = t_z_yyz_xz_xz[i] - rab_z[i] * t_z_yy_xz_xz[i];

        t_zz_yy_xz_xy[i] = t_z_yyz_xz_xy[i] - rab_z[i] * t_z_yy_xz_xy[i];

        t_zz_yy_xz_xx[i] = t_z_yyz_xz_xx[i] - rab_z[i] * t_z_yy_xz_xx[i];

        t_zz_yy_xy_zz[i] = t_z_yyz_xy_zz[i] - rab_z[i] * t_z_yy_xy_zz[i];

        t_zz_yy_xy_yz[i] = t_z_yyz_xy_yz[i] - rab_z[i] * t_z_yy_xy_yz[i];

        t_zz_yy_xy_yy[i] = t_z_yyz_xy_yy[i] - rab_z[i] * t_z_yy_xy_yy[i];

        t_zz_yy_xy_xz[i] = t_z_yyz_xy_xz[i] - rab_z[i] * t_z_yy_xy_xz[i];

        t_zz_yy_xy_xy[i] = t_z_yyz_xy_xy[i] - rab_z[i] * t_z_yy_xy_xy[i];

        t_zz_yy_xy_xx[i] = t_z_yyz_xy_xx[i] - rab_z[i] * t_z_yy_xy_xx[i];

        t_zz_yy_xx_zz[i] = t_z_yyz_xx_zz[i] - rab_z[i] * t_z_yy_xx_zz[i];

        t_zz_yy_xx_yz[i] = t_z_yyz_xx_yz[i] - rab_z[i] * t_z_yy_xx_yz[i];

        t_zz_yy_xx_yy[i] = t_z_yyz_xx_yy[i] - rab_z[i] * t_z_yy_xx_yy[i];

        t_zz_yy_xx_xz[i] = t_z_yyz_xx_xz[i] - rab_z[i] * t_z_yy_xx_xz[i];

        t_zz_yy_xx_xy[i] = t_z_yyz_xx_xy[i] - rab_z[i] * t_z_yy_xx_xy[i];

        t_zz_yy_xx_xx[i] = t_z_yyz_xx_xx[i] - rab_z[i] * t_z_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_xz_xx_xx, t_z_xz_xx_xy, t_z_xz_xx_xz, t_z_xz_xx_yy,\
                           t_z_xz_xx_yz, t_z_xz_xx_zz, t_z_xz_xy_xx, t_z_xz_xy_xy, t_z_xz_xy_xz,\
                           t_z_xz_xy_yy, t_z_xz_xy_yz, t_z_xz_xy_zz, t_z_xz_xz_xx, t_z_xz_xz_xy,\
                           t_z_xz_xz_xz, t_z_xz_xz_yy, t_z_xz_xz_yz, t_z_xz_xz_zz, t_z_xz_yy_xx,\
                           t_z_xz_yy_xy, t_z_xz_yy_xz, t_z_xz_yy_yy, t_z_xz_yy_yz, t_z_xz_yy_zz,\
                           t_z_xz_yz_xx, t_z_xz_yz_xy, t_z_xz_yz_xz, t_z_xz_yz_yy, t_z_xz_yz_yz,\
                           t_z_xz_yz_zz, t_z_xz_zz_xx, t_z_xz_zz_xy, t_z_xz_zz_xz, t_z_xz_zz_yy,\
                           t_z_xz_zz_yz, t_z_xz_zz_zz, t_z_xzz_xx_xx, t_z_xzz_xx_xy,\
                           t_z_xzz_xx_xz, t_z_xzz_xx_yy, t_z_xzz_xx_yz, t_z_xzz_xx_zz,\
                           t_z_xzz_xy_xx, t_z_xzz_xy_xy, t_z_xzz_xy_xz, t_z_xzz_xy_yy,\
                           t_z_xzz_xy_yz, t_z_xzz_xy_zz, t_z_xzz_xz_xx, t_z_xzz_xz_xy,\
                           t_z_xzz_xz_xz, t_z_xzz_xz_yy, t_z_xzz_xz_yz, t_z_xzz_xz_zz,\
                           t_z_xzz_yy_xx, t_z_xzz_yy_xy, t_z_xzz_yy_xz, t_z_xzz_yy_yy,\
                           t_z_xzz_yy_yz, t_z_xzz_yy_zz, t_z_xzz_yz_xx, t_z_xzz_yz_xy,\
                           t_z_xzz_yz_xz, t_z_xzz_yz_yy, t_z_xzz_yz_yz, t_z_xzz_yz_zz,\
                           t_z_xzz_zz_xx, t_z_xzz_zz_xy, t_z_xzz_zz_xz, t_z_xzz_zz_yy,\
                           t_z_xzz_zz_yz, t_z_xzz_zz_zz, t_zz_xz_xx_xx, t_zz_xz_xx_xy,\
                           t_zz_xz_xx_xz, t_zz_xz_xx_yy, t_zz_xz_xx_yz, t_zz_xz_xx_zz,\
                           t_zz_xz_xy_xx, t_zz_xz_xy_xy, t_zz_xz_xy_xz, t_zz_xz_xy_yy,\
                           t_zz_xz_xy_yz, t_zz_xz_xy_zz, t_zz_xz_xz_xx, t_zz_xz_xz_xy,\
                           t_zz_xz_xz_xz, t_zz_xz_xz_yy, t_zz_xz_xz_yz, t_zz_xz_xz_zz,\
                           t_zz_xz_yy_xx, t_zz_xz_yy_xy, t_zz_xz_yy_xz, t_zz_xz_yy_yy,\
                           t_zz_xz_yy_yz, t_zz_xz_yy_zz, t_zz_xz_yz_xx, t_zz_xz_yz_xy,\
                           t_zz_xz_yz_xz, t_zz_xz_yz_yy, t_zz_xz_yz_yz, t_zz_xz_yz_zz,\
                           t_zz_xz_zz_xx, t_zz_xz_zz_xy, t_zz_xz_zz_xz, t_zz_xz_zz_yy,\
                           t_zz_xz_zz_yz, t_zz_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_xz_zz_zz[i] = t_z_xzz_zz_zz[i] - rab_z[i] * t_z_xz_zz_zz[i];

        t_zz_xz_zz_yz[i] = t_z_xzz_zz_yz[i] - rab_z[i] * t_z_xz_zz_yz[i];

        t_zz_xz_zz_yy[i] = t_z_xzz_zz_yy[i] - rab_z[i] * t_z_xz_zz_yy[i];

        t_zz_xz_zz_xz[i] = t_z_xzz_zz_xz[i] - rab_z[i] * t_z_xz_zz_xz[i];

        t_zz_xz_zz_xy[i] = t_z_xzz_zz_xy[i] - rab_z[i] * t_z_xz_zz_xy[i];

        t_zz_xz_zz_xx[i] = t_z_xzz_zz_xx[i] - rab_z[i] * t_z_xz_zz_xx[i];

        t_zz_xz_yz_zz[i] = t_z_xzz_yz_zz[i] - rab_z[i] * t_z_xz_yz_zz[i];

        t_zz_xz_yz_yz[i] = t_z_xzz_yz_yz[i] - rab_z[i] * t_z_xz_yz_yz[i];

        t_zz_xz_yz_yy[i] = t_z_xzz_yz_yy[i] - rab_z[i] * t_z_xz_yz_yy[i];

        t_zz_xz_yz_xz[i] = t_z_xzz_yz_xz[i] - rab_z[i] * t_z_xz_yz_xz[i];

        t_zz_xz_yz_xy[i] = t_z_xzz_yz_xy[i] - rab_z[i] * t_z_xz_yz_xy[i];

        t_zz_xz_yz_xx[i] = t_z_xzz_yz_xx[i] - rab_z[i] * t_z_xz_yz_xx[i];

        t_zz_xz_yy_zz[i] = t_z_xzz_yy_zz[i] - rab_z[i] * t_z_xz_yy_zz[i];

        t_zz_xz_yy_yz[i] = t_z_xzz_yy_yz[i] - rab_z[i] * t_z_xz_yy_yz[i];

        t_zz_xz_yy_yy[i] = t_z_xzz_yy_yy[i] - rab_z[i] * t_z_xz_yy_yy[i];

        t_zz_xz_yy_xz[i] = t_z_xzz_yy_xz[i] - rab_z[i] * t_z_xz_yy_xz[i];

        t_zz_xz_yy_xy[i] = t_z_xzz_yy_xy[i] - rab_z[i] * t_z_xz_yy_xy[i];

        t_zz_xz_yy_xx[i] = t_z_xzz_yy_xx[i] - rab_z[i] * t_z_xz_yy_xx[i];

        t_zz_xz_xz_zz[i] = t_z_xzz_xz_zz[i] - rab_z[i] * t_z_xz_xz_zz[i];

        t_zz_xz_xz_yz[i] = t_z_xzz_xz_yz[i] - rab_z[i] * t_z_xz_xz_yz[i];

        t_zz_xz_xz_yy[i] = t_z_xzz_xz_yy[i] - rab_z[i] * t_z_xz_xz_yy[i];

        t_zz_xz_xz_xz[i] = t_z_xzz_xz_xz[i] - rab_z[i] * t_z_xz_xz_xz[i];

        t_zz_xz_xz_xy[i] = t_z_xzz_xz_xy[i] - rab_z[i] * t_z_xz_xz_xy[i];

        t_zz_xz_xz_xx[i] = t_z_xzz_xz_xx[i] - rab_z[i] * t_z_xz_xz_xx[i];

        t_zz_xz_xy_zz[i] = t_z_xzz_xy_zz[i] - rab_z[i] * t_z_xz_xy_zz[i];

        t_zz_xz_xy_yz[i] = t_z_xzz_xy_yz[i] - rab_z[i] * t_z_xz_xy_yz[i];

        t_zz_xz_xy_yy[i] = t_z_xzz_xy_yy[i] - rab_z[i] * t_z_xz_xy_yy[i];

        t_zz_xz_xy_xz[i] = t_z_xzz_xy_xz[i] - rab_z[i] * t_z_xz_xy_xz[i];

        t_zz_xz_xy_xy[i] = t_z_xzz_xy_xy[i] - rab_z[i] * t_z_xz_xy_xy[i];

        t_zz_xz_xy_xx[i] = t_z_xzz_xy_xx[i] - rab_z[i] * t_z_xz_xy_xx[i];

        t_zz_xz_xx_zz[i] = t_z_xzz_xx_zz[i] - rab_z[i] * t_z_xz_xx_zz[i];

        t_zz_xz_xx_yz[i] = t_z_xzz_xx_yz[i] - rab_z[i] * t_z_xz_xx_yz[i];

        t_zz_xz_xx_yy[i] = t_z_xzz_xx_yy[i] - rab_z[i] * t_z_xz_xx_yy[i];

        t_zz_xz_xx_xz[i] = t_z_xzz_xx_xz[i] - rab_z[i] * t_z_xz_xx_xz[i];

        t_zz_xz_xx_xy[i] = t_z_xzz_xx_xy[i] - rab_z[i] * t_z_xz_xx_xy[i];

        t_zz_xz_xx_xx[i] = t_z_xzz_xx_xx[i] - rab_z[i] * t_z_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_xy_xx_xx, t_z_xy_xx_xy, t_z_xy_xx_xz, t_z_xy_xx_yy,\
                           t_z_xy_xx_yz, t_z_xy_xx_zz, t_z_xy_xy_xx, t_z_xy_xy_xy, t_z_xy_xy_xz,\
                           t_z_xy_xy_yy, t_z_xy_xy_yz, t_z_xy_xy_zz, t_z_xy_xz_xx, t_z_xy_xz_xy,\
                           t_z_xy_xz_xz, t_z_xy_xz_yy, t_z_xy_xz_yz, t_z_xy_xz_zz, t_z_xy_yy_xx,\
                           t_z_xy_yy_xy, t_z_xy_yy_xz, t_z_xy_yy_yy, t_z_xy_yy_yz, t_z_xy_yy_zz,\
                           t_z_xy_yz_xx, t_z_xy_yz_xy, t_z_xy_yz_xz, t_z_xy_yz_yy, t_z_xy_yz_yz,\
                           t_z_xy_yz_zz, t_z_xy_zz_xx, t_z_xy_zz_xy, t_z_xy_zz_xz, t_z_xy_zz_yy,\
                           t_z_xy_zz_yz, t_z_xy_zz_zz, t_z_xyz_xx_xx, t_z_xyz_xx_xy,\
                           t_z_xyz_xx_xz, t_z_xyz_xx_yy, t_z_xyz_xx_yz, t_z_xyz_xx_zz,\
                           t_z_xyz_xy_xx, t_z_xyz_xy_xy, t_z_xyz_xy_xz, t_z_xyz_xy_yy,\
                           t_z_xyz_xy_yz, t_z_xyz_xy_zz, t_z_xyz_xz_xx, t_z_xyz_xz_xy,\
                           t_z_xyz_xz_xz, t_z_xyz_xz_yy, t_z_xyz_xz_yz, t_z_xyz_xz_zz,\
                           t_z_xyz_yy_xx, t_z_xyz_yy_xy, t_z_xyz_yy_xz, t_z_xyz_yy_yy,\
                           t_z_xyz_yy_yz, t_z_xyz_yy_zz, t_z_xyz_yz_xx, t_z_xyz_yz_xy,\
                           t_z_xyz_yz_xz, t_z_xyz_yz_yy, t_z_xyz_yz_yz, t_z_xyz_yz_zz,\
                           t_z_xyz_zz_xx, t_z_xyz_zz_xy, t_z_xyz_zz_xz, t_z_xyz_zz_yy,\
                           t_z_xyz_zz_yz, t_z_xyz_zz_zz, t_zz_xy_xx_xx, t_zz_xy_xx_xy,\
                           t_zz_xy_xx_xz, t_zz_xy_xx_yy, t_zz_xy_xx_yz, t_zz_xy_xx_zz,\
                           t_zz_xy_xy_xx, t_zz_xy_xy_xy, t_zz_xy_xy_xz, t_zz_xy_xy_yy,\
                           t_zz_xy_xy_yz, t_zz_xy_xy_zz, t_zz_xy_xz_xx, t_zz_xy_xz_xy,\
                           t_zz_xy_xz_xz, t_zz_xy_xz_yy, t_zz_xy_xz_yz, t_zz_xy_xz_zz,\
                           t_zz_xy_yy_xx, t_zz_xy_yy_xy, t_zz_xy_yy_xz, t_zz_xy_yy_yy,\
                           t_zz_xy_yy_yz, t_zz_xy_yy_zz, t_zz_xy_yz_xx, t_zz_xy_yz_xy,\
                           t_zz_xy_yz_xz, t_zz_xy_yz_yy, t_zz_xy_yz_yz, t_zz_xy_yz_zz,\
                           t_zz_xy_zz_xx, t_zz_xy_zz_xy, t_zz_xy_zz_xz, t_zz_xy_zz_yy,\
                           t_zz_xy_zz_yz, t_zz_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_xy_zz_zz[i] = t_z_xyz_zz_zz[i] - rab_z[i] * t_z_xy_zz_zz[i];

        t_zz_xy_zz_yz[i] = t_z_xyz_zz_yz[i] - rab_z[i] * t_z_xy_zz_yz[i];

        t_zz_xy_zz_yy[i] = t_z_xyz_zz_yy[i] - rab_z[i] * t_z_xy_zz_yy[i];

        t_zz_xy_zz_xz[i] = t_z_xyz_zz_xz[i] - rab_z[i] * t_z_xy_zz_xz[i];

        t_zz_xy_zz_xy[i] = t_z_xyz_zz_xy[i] - rab_z[i] * t_z_xy_zz_xy[i];

        t_zz_xy_zz_xx[i] = t_z_xyz_zz_xx[i] - rab_z[i] * t_z_xy_zz_xx[i];

        t_zz_xy_yz_zz[i] = t_z_xyz_yz_zz[i] - rab_z[i] * t_z_xy_yz_zz[i];

        t_zz_xy_yz_yz[i] = t_z_xyz_yz_yz[i] - rab_z[i] * t_z_xy_yz_yz[i];

        t_zz_xy_yz_yy[i] = t_z_xyz_yz_yy[i] - rab_z[i] * t_z_xy_yz_yy[i];

        t_zz_xy_yz_xz[i] = t_z_xyz_yz_xz[i] - rab_z[i] * t_z_xy_yz_xz[i];

        t_zz_xy_yz_xy[i] = t_z_xyz_yz_xy[i] - rab_z[i] * t_z_xy_yz_xy[i];

        t_zz_xy_yz_xx[i] = t_z_xyz_yz_xx[i] - rab_z[i] * t_z_xy_yz_xx[i];

        t_zz_xy_yy_zz[i] = t_z_xyz_yy_zz[i] - rab_z[i] * t_z_xy_yy_zz[i];

        t_zz_xy_yy_yz[i] = t_z_xyz_yy_yz[i] - rab_z[i] * t_z_xy_yy_yz[i];

        t_zz_xy_yy_yy[i] = t_z_xyz_yy_yy[i] - rab_z[i] * t_z_xy_yy_yy[i];

        t_zz_xy_yy_xz[i] = t_z_xyz_yy_xz[i] - rab_z[i] * t_z_xy_yy_xz[i];

        t_zz_xy_yy_xy[i] = t_z_xyz_yy_xy[i] - rab_z[i] * t_z_xy_yy_xy[i];

        t_zz_xy_yy_xx[i] = t_z_xyz_yy_xx[i] - rab_z[i] * t_z_xy_yy_xx[i];

        t_zz_xy_xz_zz[i] = t_z_xyz_xz_zz[i] - rab_z[i] * t_z_xy_xz_zz[i];

        t_zz_xy_xz_yz[i] = t_z_xyz_xz_yz[i] - rab_z[i] * t_z_xy_xz_yz[i];

        t_zz_xy_xz_yy[i] = t_z_xyz_xz_yy[i] - rab_z[i] * t_z_xy_xz_yy[i];

        t_zz_xy_xz_xz[i] = t_z_xyz_xz_xz[i] - rab_z[i] * t_z_xy_xz_xz[i];

        t_zz_xy_xz_xy[i] = t_z_xyz_xz_xy[i] - rab_z[i] * t_z_xy_xz_xy[i];

        t_zz_xy_xz_xx[i] = t_z_xyz_xz_xx[i] - rab_z[i] * t_z_xy_xz_xx[i];

        t_zz_xy_xy_zz[i] = t_z_xyz_xy_zz[i] - rab_z[i] * t_z_xy_xy_zz[i];

        t_zz_xy_xy_yz[i] = t_z_xyz_xy_yz[i] - rab_z[i] * t_z_xy_xy_yz[i];

        t_zz_xy_xy_yy[i] = t_z_xyz_xy_yy[i] - rab_z[i] * t_z_xy_xy_yy[i];

        t_zz_xy_xy_xz[i] = t_z_xyz_xy_xz[i] - rab_z[i] * t_z_xy_xy_xz[i];

        t_zz_xy_xy_xy[i] = t_z_xyz_xy_xy[i] - rab_z[i] * t_z_xy_xy_xy[i];

        t_zz_xy_xy_xx[i] = t_z_xyz_xy_xx[i] - rab_z[i] * t_z_xy_xy_xx[i];

        t_zz_xy_xx_zz[i] = t_z_xyz_xx_zz[i] - rab_z[i] * t_z_xy_xx_zz[i];

        t_zz_xy_xx_yz[i] = t_z_xyz_xx_yz[i] - rab_z[i] * t_z_xy_xx_yz[i];

        t_zz_xy_xx_yy[i] = t_z_xyz_xx_yy[i] - rab_z[i] * t_z_xy_xx_yy[i];

        t_zz_xy_xx_xz[i] = t_z_xyz_xx_xz[i] - rab_z[i] * t_z_xy_xx_xz[i];

        t_zz_xy_xx_xy[i] = t_z_xyz_xx_xy[i] - rab_z[i] * t_z_xy_xx_xy[i];

        t_zz_xy_xx_xx[i] = t_z_xyz_xx_xx[i] - rab_z[i] * t_z_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_xx_xx_xx, t_z_xx_xx_xy, t_z_xx_xx_xz, t_z_xx_xx_yy,\
                           t_z_xx_xx_yz, t_z_xx_xx_zz, t_z_xx_xy_xx, t_z_xx_xy_xy, t_z_xx_xy_xz,\
                           t_z_xx_xy_yy, t_z_xx_xy_yz, t_z_xx_xy_zz, t_z_xx_xz_xx, t_z_xx_xz_xy,\
                           t_z_xx_xz_xz, t_z_xx_xz_yy, t_z_xx_xz_yz, t_z_xx_xz_zz, t_z_xx_yy_xx,\
                           t_z_xx_yy_xy, t_z_xx_yy_xz, t_z_xx_yy_yy, t_z_xx_yy_yz, t_z_xx_yy_zz,\
                           t_z_xx_yz_xx, t_z_xx_yz_xy, t_z_xx_yz_xz, t_z_xx_yz_yy, t_z_xx_yz_yz,\
                           t_z_xx_yz_zz, t_z_xx_zz_xx, t_z_xx_zz_xy, t_z_xx_zz_xz, t_z_xx_zz_yy,\
                           t_z_xx_zz_yz, t_z_xx_zz_zz, t_z_xxz_xx_xx, t_z_xxz_xx_xy,\
                           t_z_xxz_xx_xz, t_z_xxz_xx_yy, t_z_xxz_xx_yz, t_z_xxz_xx_zz,\
                           t_z_xxz_xy_xx, t_z_xxz_xy_xy, t_z_xxz_xy_xz, t_z_xxz_xy_yy,\
                           t_z_xxz_xy_yz, t_z_xxz_xy_zz, t_z_xxz_xz_xx, t_z_xxz_xz_xy,\
                           t_z_xxz_xz_xz, t_z_xxz_xz_yy, t_z_xxz_xz_yz, t_z_xxz_xz_zz,\
                           t_z_xxz_yy_xx, t_z_xxz_yy_xy, t_z_xxz_yy_xz, t_z_xxz_yy_yy,\
                           t_z_xxz_yy_yz, t_z_xxz_yy_zz, t_z_xxz_yz_xx, t_z_xxz_yz_xy,\
                           t_z_xxz_yz_xz, t_z_xxz_yz_yy, t_z_xxz_yz_yz, t_z_xxz_yz_zz,\
                           t_z_xxz_zz_xx, t_z_xxz_zz_xy, t_z_xxz_zz_xz, t_z_xxz_zz_yy,\
                           t_z_xxz_zz_yz, t_z_xxz_zz_zz, t_zz_xx_xx_xx, t_zz_xx_xx_xy,\
                           t_zz_xx_xx_xz, t_zz_xx_xx_yy, t_zz_xx_xx_yz, t_zz_xx_xx_zz,\
                           t_zz_xx_xy_xx, t_zz_xx_xy_xy, t_zz_xx_xy_xz, t_zz_xx_xy_yy,\
                           t_zz_xx_xy_yz, t_zz_xx_xy_zz, t_zz_xx_xz_xx, t_zz_xx_xz_xy,\
                           t_zz_xx_xz_xz, t_zz_xx_xz_yy, t_zz_xx_xz_yz, t_zz_xx_xz_zz,\
                           t_zz_xx_yy_xx, t_zz_xx_yy_xy, t_zz_xx_yy_xz, t_zz_xx_yy_yy,\
                           t_zz_xx_yy_yz, t_zz_xx_yy_zz, t_zz_xx_yz_xx, t_zz_xx_yz_xy,\
                           t_zz_xx_yz_xz, t_zz_xx_yz_yy, t_zz_xx_yz_yz, t_zz_xx_yz_zz,\
                           t_zz_xx_zz_xx, t_zz_xx_zz_xy, t_zz_xx_zz_xz, t_zz_xx_zz_yy,\
                           t_zz_xx_zz_yz, t_zz_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_xx_zz_zz[i] = t_z_xxz_zz_zz[i] - rab_z[i] * t_z_xx_zz_zz[i];

        t_zz_xx_zz_yz[i] = t_z_xxz_zz_yz[i] - rab_z[i] * t_z_xx_zz_yz[i];

        t_zz_xx_zz_yy[i] = t_z_xxz_zz_yy[i] - rab_z[i] * t_z_xx_zz_yy[i];

        t_zz_xx_zz_xz[i] = t_z_xxz_zz_xz[i] - rab_z[i] * t_z_xx_zz_xz[i];

        t_zz_xx_zz_xy[i] = t_z_xxz_zz_xy[i] - rab_z[i] * t_z_xx_zz_xy[i];

        t_zz_xx_zz_xx[i] = t_z_xxz_zz_xx[i] - rab_z[i] * t_z_xx_zz_xx[i];

        t_zz_xx_yz_zz[i] = t_z_xxz_yz_zz[i] - rab_z[i] * t_z_xx_yz_zz[i];

        t_zz_xx_yz_yz[i] = t_z_xxz_yz_yz[i] - rab_z[i] * t_z_xx_yz_yz[i];

        t_zz_xx_yz_yy[i] = t_z_xxz_yz_yy[i] - rab_z[i] * t_z_xx_yz_yy[i];

        t_zz_xx_yz_xz[i] = t_z_xxz_yz_xz[i] - rab_z[i] * t_z_xx_yz_xz[i];

        t_zz_xx_yz_xy[i] = t_z_xxz_yz_xy[i] - rab_z[i] * t_z_xx_yz_xy[i];

        t_zz_xx_yz_xx[i] = t_z_xxz_yz_xx[i] - rab_z[i] * t_z_xx_yz_xx[i];

        t_zz_xx_yy_zz[i] = t_z_xxz_yy_zz[i] - rab_z[i] * t_z_xx_yy_zz[i];

        t_zz_xx_yy_yz[i] = t_z_xxz_yy_yz[i] - rab_z[i] * t_z_xx_yy_yz[i];

        t_zz_xx_yy_yy[i] = t_z_xxz_yy_yy[i] - rab_z[i] * t_z_xx_yy_yy[i];

        t_zz_xx_yy_xz[i] = t_z_xxz_yy_xz[i] - rab_z[i] * t_z_xx_yy_xz[i];

        t_zz_xx_yy_xy[i] = t_z_xxz_yy_xy[i] - rab_z[i] * t_z_xx_yy_xy[i];

        t_zz_xx_yy_xx[i] = t_z_xxz_yy_xx[i] - rab_z[i] * t_z_xx_yy_xx[i];

        t_zz_xx_xz_zz[i] = t_z_xxz_xz_zz[i] - rab_z[i] * t_z_xx_xz_zz[i];

        t_zz_xx_xz_yz[i] = t_z_xxz_xz_yz[i] - rab_z[i] * t_z_xx_xz_yz[i];

        t_zz_xx_xz_yy[i] = t_z_xxz_xz_yy[i] - rab_z[i] * t_z_xx_xz_yy[i];

        t_zz_xx_xz_xz[i] = t_z_xxz_xz_xz[i] - rab_z[i] * t_z_xx_xz_xz[i];

        t_zz_xx_xz_xy[i] = t_z_xxz_xz_xy[i] - rab_z[i] * t_z_xx_xz_xy[i];

        t_zz_xx_xz_xx[i] = t_z_xxz_xz_xx[i] - rab_z[i] * t_z_xx_xz_xx[i];

        t_zz_xx_xy_zz[i] = t_z_xxz_xy_zz[i] - rab_z[i] * t_z_xx_xy_zz[i];

        t_zz_xx_xy_yz[i] = t_z_xxz_xy_yz[i] - rab_z[i] * t_z_xx_xy_yz[i];

        t_zz_xx_xy_yy[i] = t_z_xxz_xy_yy[i] - rab_z[i] * t_z_xx_xy_yy[i];

        t_zz_xx_xy_xz[i] = t_z_xxz_xy_xz[i] - rab_z[i] * t_z_xx_xy_xz[i];

        t_zz_xx_xy_xy[i] = t_z_xxz_xy_xy[i] - rab_z[i] * t_z_xx_xy_xy[i];

        t_zz_xx_xy_xx[i] = t_z_xxz_xy_xx[i] - rab_z[i] * t_z_xx_xy_xx[i];

        t_zz_xx_xx_zz[i] = t_z_xxz_xx_zz[i] - rab_z[i] * t_z_xx_xx_zz[i];

        t_zz_xx_xx_yz[i] = t_z_xxz_xx_yz[i] - rab_z[i] * t_z_xx_xx_yz[i];

        t_zz_xx_xx_yy[i] = t_z_xxz_xx_yy[i] - rab_z[i] * t_z_xx_xx_yy[i];

        t_zz_xx_xx_xz[i] = t_z_xxz_xx_xz[i] - rab_z[i] * t_z_xx_xx_xz[i];

        t_zz_xx_xx_xy[i] = t_z_xxz_xx_xy[i] - rab_z[i] * t_z_xx_xx_xy[i];

        t_zz_xx_xx_xx[i] = t_z_xxz_xx_xx[i] - rab_z[i] * t_z_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_zz_xx_xx, t_y_zz_xx_xy, t_y_zz_xx_xz, t_y_zz_xx_yy,\
                           t_y_zz_xx_yz, t_y_zz_xx_zz, t_y_zz_xy_xx, t_y_zz_xy_xy, t_y_zz_xy_xz,\
                           t_y_zz_xy_yy, t_y_zz_xy_yz, t_y_zz_xy_zz, t_y_zz_xz_xx, t_y_zz_xz_xy,\
                           t_y_zz_xz_xz, t_y_zz_xz_yy, t_y_zz_xz_yz, t_y_zz_xz_zz, t_y_zz_yy_xx,\
                           t_y_zz_yy_xy, t_y_zz_yy_xz, t_y_zz_yy_yy, t_y_zz_yy_yz, t_y_zz_yy_zz,\
                           t_y_zz_yz_xx, t_y_zz_yz_xy, t_y_zz_yz_xz, t_y_zz_yz_yy, t_y_zz_yz_yz,\
                           t_y_zz_yz_zz, t_y_zz_zz_xx, t_y_zz_zz_xy, t_y_zz_zz_xz, t_y_zz_zz_yy,\
                           t_y_zz_zz_yz, t_y_zz_zz_zz, t_y_zzz_xx_xx, t_y_zzz_xx_xy,\
                           t_y_zzz_xx_xz, t_y_zzz_xx_yy, t_y_zzz_xx_yz, t_y_zzz_xx_zz,\
                           t_y_zzz_xy_xx, t_y_zzz_xy_xy, t_y_zzz_xy_xz, t_y_zzz_xy_yy,\
                           t_y_zzz_xy_yz, t_y_zzz_xy_zz, t_y_zzz_xz_xx, t_y_zzz_xz_xy,\
                           t_y_zzz_xz_xz, t_y_zzz_xz_yy, t_y_zzz_xz_yz, t_y_zzz_xz_zz,\
                           t_y_zzz_yy_xx, t_y_zzz_yy_xy, t_y_zzz_yy_xz, t_y_zzz_yy_yy,\
                           t_y_zzz_yy_yz, t_y_zzz_yy_zz, t_y_zzz_yz_xx, t_y_zzz_yz_xy,\
                           t_y_zzz_yz_xz, t_y_zzz_yz_yy, t_y_zzz_yz_yz, t_y_zzz_yz_zz,\
                           t_y_zzz_zz_xx, t_y_zzz_zz_xy, t_y_zzz_zz_xz, t_y_zzz_zz_yy,\
                           t_y_zzz_zz_yz, t_y_zzz_zz_zz, t_yz_zz_xx_xx, t_yz_zz_xx_xy,\
                           t_yz_zz_xx_xz, t_yz_zz_xx_yy, t_yz_zz_xx_yz, t_yz_zz_xx_zz,\
                           t_yz_zz_xy_xx, t_yz_zz_xy_xy, t_yz_zz_xy_xz, t_yz_zz_xy_yy,\
                           t_yz_zz_xy_yz, t_yz_zz_xy_zz, t_yz_zz_xz_xx, t_yz_zz_xz_xy,\
                           t_yz_zz_xz_xz, t_yz_zz_xz_yy, t_yz_zz_xz_yz, t_yz_zz_xz_zz,\
                           t_yz_zz_yy_xx, t_yz_zz_yy_xy, t_yz_zz_yy_xz, t_yz_zz_yy_yy,\
                           t_yz_zz_yy_yz, t_yz_zz_yy_zz, t_yz_zz_yz_xx, t_yz_zz_yz_xy,\
                           t_yz_zz_yz_xz, t_yz_zz_yz_yy, t_yz_zz_yz_yz, t_yz_zz_yz_zz,\
                           t_yz_zz_zz_xx, t_yz_zz_zz_xy, t_yz_zz_zz_xz, t_yz_zz_zz_yy,\
                           t_yz_zz_zz_yz, t_yz_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_zz_zz_zz[i] = t_y_zzz_zz_zz[i] - rab_z[i] * t_y_zz_zz_zz[i];

        t_yz_zz_zz_yz[i] = t_y_zzz_zz_yz[i] - rab_z[i] * t_y_zz_zz_yz[i];

        t_yz_zz_zz_yy[i] = t_y_zzz_zz_yy[i] - rab_z[i] * t_y_zz_zz_yy[i];

        t_yz_zz_zz_xz[i] = t_y_zzz_zz_xz[i] - rab_z[i] * t_y_zz_zz_xz[i];

        t_yz_zz_zz_xy[i] = t_y_zzz_zz_xy[i] - rab_z[i] * t_y_zz_zz_xy[i];

        t_yz_zz_zz_xx[i] = t_y_zzz_zz_xx[i] - rab_z[i] * t_y_zz_zz_xx[i];

        t_yz_zz_yz_zz[i] = t_y_zzz_yz_zz[i] - rab_z[i] * t_y_zz_yz_zz[i];

        t_yz_zz_yz_yz[i] = t_y_zzz_yz_yz[i] - rab_z[i] * t_y_zz_yz_yz[i];

        t_yz_zz_yz_yy[i] = t_y_zzz_yz_yy[i] - rab_z[i] * t_y_zz_yz_yy[i];

        t_yz_zz_yz_xz[i] = t_y_zzz_yz_xz[i] - rab_z[i] * t_y_zz_yz_xz[i];

        t_yz_zz_yz_xy[i] = t_y_zzz_yz_xy[i] - rab_z[i] * t_y_zz_yz_xy[i];

        t_yz_zz_yz_xx[i] = t_y_zzz_yz_xx[i] - rab_z[i] * t_y_zz_yz_xx[i];

        t_yz_zz_yy_zz[i] = t_y_zzz_yy_zz[i] - rab_z[i] * t_y_zz_yy_zz[i];

        t_yz_zz_yy_yz[i] = t_y_zzz_yy_yz[i] - rab_z[i] * t_y_zz_yy_yz[i];

        t_yz_zz_yy_yy[i] = t_y_zzz_yy_yy[i] - rab_z[i] * t_y_zz_yy_yy[i];

        t_yz_zz_yy_xz[i] = t_y_zzz_yy_xz[i] - rab_z[i] * t_y_zz_yy_xz[i];

        t_yz_zz_yy_xy[i] = t_y_zzz_yy_xy[i] - rab_z[i] * t_y_zz_yy_xy[i];

        t_yz_zz_yy_xx[i] = t_y_zzz_yy_xx[i] - rab_z[i] * t_y_zz_yy_xx[i];

        t_yz_zz_xz_zz[i] = t_y_zzz_xz_zz[i] - rab_z[i] * t_y_zz_xz_zz[i];

        t_yz_zz_xz_yz[i] = t_y_zzz_xz_yz[i] - rab_z[i] * t_y_zz_xz_yz[i];

        t_yz_zz_xz_yy[i] = t_y_zzz_xz_yy[i] - rab_z[i] * t_y_zz_xz_yy[i];

        t_yz_zz_xz_xz[i] = t_y_zzz_xz_xz[i] - rab_z[i] * t_y_zz_xz_xz[i];

        t_yz_zz_xz_xy[i] = t_y_zzz_xz_xy[i] - rab_z[i] * t_y_zz_xz_xy[i];

        t_yz_zz_xz_xx[i] = t_y_zzz_xz_xx[i] - rab_z[i] * t_y_zz_xz_xx[i];

        t_yz_zz_xy_zz[i] = t_y_zzz_xy_zz[i] - rab_z[i] * t_y_zz_xy_zz[i];

        t_yz_zz_xy_yz[i] = t_y_zzz_xy_yz[i] - rab_z[i] * t_y_zz_xy_yz[i];

        t_yz_zz_xy_yy[i] = t_y_zzz_xy_yy[i] - rab_z[i] * t_y_zz_xy_yy[i];

        t_yz_zz_xy_xz[i] = t_y_zzz_xy_xz[i] - rab_z[i] * t_y_zz_xy_xz[i];

        t_yz_zz_xy_xy[i] = t_y_zzz_xy_xy[i] - rab_z[i] * t_y_zz_xy_xy[i];

        t_yz_zz_xy_xx[i] = t_y_zzz_xy_xx[i] - rab_z[i] * t_y_zz_xy_xx[i];

        t_yz_zz_xx_zz[i] = t_y_zzz_xx_zz[i] - rab_z[i] * t_y_zz_xx_zz[i];

        t_yz_zz_xx_yz[i] = t_y_zzz_xx_yz[i] - rab_z[i] * t_y_zz_xx_yz[i];

        t_yz_zz_xx_yy[i] = t_y_zzz_xx_yy[i] - rab_z[i] * t_y_zz_xx_yy[i];

        t_yz_zz_xx_xz[i] = t_y_zzz_xx_xz[i] - rab_z[i] * t_y_zz_xx_xz[i];

        t_yz_zz_xx_xy[i] = t_y_zzz_xx_xy[i] - rab_z[i] * t_y_zz_xx_xy[i];

        t_yz_zz_xx_xx[i] = t_y_zzz_xx_xx[i] - rab_z[i] * t_y_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_yz_xx_xx, t_y_yz_xx_xy, t_y_yz_xx_xz, t_y_yz_xx_yy,\
                           t_y_yz_xx_yz, t_y_yz_xx_zz, t_y_yz_xy_xx, t_y_yz_xy_xy, t_y_yz_xy_xz,\
                           t_y_yz_xy_yy, t_y_yz_xy_yz, t_y_yz_xy_zz, t_y_yz_xz_xx, t_y_yz_xz_xy,\
                           t_y_yz_xz_xz, t_y_yz_xz_yy, t_y_yz_xz_yz, t_y_yz_xz_zz, t_y_yz_yy_xx,\
                           t_y_yz_yy_xy, t_y_yz_yy_xz, t_y_yz_yy_yy, t_y_yz_yy_yz, t_y_yz_yy_zz,\
                           t_y_yz_yz_xx, t_y_yz_yz_xy, t_y_yz_yz_xz, t_y_yz_yz_yy, t_y_yz_yz_yz,\
                           t_y_yz_yz_zz, t_y_yz_zz_xx, t_y_yz_zz_xy, t_y_yz_zz_xz, t_y_yz_zz_yy,\
                           t_y_yz_zz_yz, t_y_yz_zz_zz, t_y_yzz_xx_xx, t_y_yzz_xx_xy,\
                           t_y_yzz_xx_xz, t_y_yzz_xx_yy, t_y_yzz_xx_yz, t_y_yzz_xx_zz,\
                           t_y_yzz_xy_xx, t_y_yzz_xy_xy, t_y_yzz_xy_xz, t_y_yzz_xy_yy,\
                           t_y_yzz_xy_yz, t_y_yzz_xy_zz, t_y_yzz_xz_xx, t_y_yzz_xz_xy,\
                           t_y_yzz_xz_xz, t_y_yzz_xz_yy, t_y_yzz_xz_yz, t_y_yzz_xz_zz,\
                           t_y_yzz_yy_xx, t_y_yzz_yy_xy, t_y_yzz_yy_xz, t_y_yzz_yy_yy,\
                           t_y_yzz_yy_yz, t_y_yzz_yy_zz, t_y_yzz_yz_xx, t_y_yzz_yz_xy,\
                           t_y_yzz_yz_xz, t_y_yzz_yz_yy, t_y_yzz_yz_yz, t_y_yzz_yz_zz,\
                           t_y_yzz_zz_xx, t_y_yzz_zz_xy, t_y_yzz_zz_xz, t_y_yzz_zz_yy,\
                           t_y_yzz_zz_yz, t_y_yzz_zz_zz, t_yz_yz_xx_xx, t_yz_yz_xx_xy,\
                           t_yz_yz_xx_xz, t_yz_yz_xx_yy, t_yz_yz_xx_yz, t_yz_yz_xx_zz,\
                           t_yz_yz_xy_xx, t_yz_yz_xy_xy, t_yz_yz_xy_xz, t_yz_yz_xy_yy,\
                           t_yz_yz_xy_yz, t_yz_yz_xy_zz, t_yz_yz_xz_xx, t_yz_yz_xz_xy,\
                           t_yz_yz_xz_xz, t_yz_yz_xz_yy, t_yz_yz_xz_yz, t_yz_yz_xz_zz,\
                           t_yz_yz_yy_xx, t_yz_yz_yy_xy, t_yz_yz_yy_xz, t_yz_yz_yy_yy,\
                           t_yz_yz_yy_yz, t_yz_yz_yy_zz, t_yz_yz_yz_xx, t_yz_yz_yz_xy,\
                           t_yz_yz_yz_xz, t_yz_yz_yz_yy, t_yz_yz_yz_yz, t_yz_yz_yz_zz,\
                           t_yz_yz_zz_xx, t_yz_yz_zz_xy, t_yz_yz_zz_xz, t_yz_yz_zz_yy,\
                           t_yz_yz_zz_yz, t_yz_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_yz_zz_zz[i] = t_y_yzz_zz_zz[i] - rab_z[i] * t_y_yz_zz_zz[i];

        t_yz_yz_zz_yz[i] = t_y_yzz_zz_yz[i] - rab_z[i] * t_y_yz_zz_yz[i];

        t_yz_yz_zz_yy[i] = t_y_yzz_zz_yy[i] - rab_z[i] * t_y_yz_zz_yy[i];

        t_yz_yz_zz_xz[i] = t_y_yzz_zz_xz[i] - rab_z[i] * t_y_yz_zz_xz[i];

        t_yz_yz_zz_xy[i] = t_y_yzz_zz_xy[i] - rab_z[i] * t_y_yz_zz_xy[i];

        t_yz_yz_zz_xx[i] = t_y_yzz_zz_xx[i] - rab_z[i] * t_y_yz_zz_xx[i];

        t_yz_yz_yz_zz[i] = t_y_yzz_yz_zz[i] - rab_z[i] * t_y_yz_yz_zz[i];

        t_yz_yz_yz_yz[i] = t_y_yzz_yz_yz[i] - rab_z[i] * t_y_yz_yz_yz[i];

        t_yz_yz_yz_yy[i] = t_y_yzz_yz_yy[i] - rab_z[i] * t_y_yz_yz_yy[i];

        t_yz_yz_yz_xz[i] = t_y_yzz_yz_xz[i] - rab_z[i] * t_y_yz_yz_xz[i];

        t_yz_yz_yz_xy[i] = t_y_yzz_yz_xy[i] - rab_z[i] * t_y_yz_yz_xy[i];

        t_yz_yz_yz_xx[i] = t_y_yzz_yz_xx[i] - rab_z[i] * t_y_yz_yz_xx[i];

        t_yz_yz_yy_zz[i] = t_y_yzz_yy_zz[i] - rab_z[i] * t_y_yz_yy_zz[i];

        t_yz_yz_yy_yz[i] = t_y_yzz_yy_yz[i] - rab_z[i] * t_y_yz_yy_yz[i];

        t_yz_yz_yy_yy[i] = t_y_yzz_yy_yy[i] - rab_z[i] * t_y_yz_yy_yy[i];

        t_yz_yz_yy_xz[i] = t_y_yzz_yy_xz[i] - rab_z[i] * t_y_yz_yy_xz[i];

        t_yz_yz_yy_xy[i] = t_y_yzz_yy_xy[i] - rab_z[i] * t_y_yz_yy_xy[i];

        t_yz_yz_yy_xx[i] = t_y_yzz_yy_xx[i] - rab_z[i] * t_y_yz_yy_xx[i];

        t_yz_yz_xz_zz[i] = t_y_yzz_xz_zz[i] - rab_z[i] * t_y_yz_xz_zz[i];

        t_yz_yz_xz_yz[i] = t_y_yzz_xz_yz[i] - rab_z[i] * t_y_yz_xz_yz[i];

        t_yz_yz_xz_yy[i] = t_y_yzz_xz_yy[i] - rab_z[i] * t_y_yz_xz_yy[i];

        t_yz_yz_xz_xz[i] = t_y_yzz_xz_xz[i] - rab_z[i] * t_y_yz_xz_xz[i];

        t_yz_yz_xz_xy[i] = t_y_yzz_xz_xy[i] - rab_z[i] * t_y_yz_xz_xy[i];

        t_yz_yz_xz_xx[i] = t_y_yzz_xz_xx[i] - rab_z[i] * t_y_yz_xz_xx[i];

        t_yz_yz_xy_zz[i] = t_y_yzz_xy_zz[i] - rab_z[i] * t_y_yz_xy_zz[i];

        t_yz_yz_xy_yz[i] = t_y_yzz_xy_yz[i] - rab_z[i] * t_y_yz_xy_yz[i];

        t_yz_yz_xy_yy[i] = t_y_yzz_xy_yy[i] - rab_z[i] * t_y_yz_xy_yy[i];

        t_yz_yz_xy_xz[i] = t_y_yzz_xy_xz[i] - rab_z[i] * t_y_yz_xy_xz[i];

        t_yz_yz_xy_xy[i] = t_y_yzz_xy_xy[i] - rab_z[i] * t_y_yz_xy_xy[i];

        t_yz_yz_xy_xx[i] = t_y_yzz_xy_xx[i] - rab_z[i] * t_y_yz_xy_xx[i];

        t_yz_yz_xx_zz[i] = t_y_yzz_xx_zz[i] - rab_z[i] * t_y_yz_xx_zz[i];

        t_yz_yz_xx_yz[i] = t_y_yzz_xx_yz[i] - rab_z[i] * t_y_yz_xx_yz[i];

        t_yz_yz_xx_yy[i] = t_y_yzz_xx_yy[i] - rab_z[i] * t_y_yz_xx_yy[i];

        t_yz_yz_xx_xz[i] = t_y_yzz_xx_xz[i] - rab_z[i] * t_y_yz_xx_xz[i];

        t_yz_yz_xx_xy[i] = t_y_yzz_xx_xy[i] - rab_z[i] * t_y_yz_xx_xy[i];

        t_yz_yz_xx_xx[i] = t_y_yzz_xx_xx[i] - rab_z[i] * t_y_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_yy_xx_xx, t_y_yy_xx_xy, t_y_yy_xx_xz, t_y_yy_xx_yy,\
                           t_y_yy_xx_yz, t_y_yy_xx_zz, t_y_yy_xy_xx, t_y_yy_xy_xy, t_y_yy_xy_xz,\
                           t_y_yy_xy_yy, t_y_yy_xy_yz, t_y_yy_xy_zz, t_y_yy_xz_xx, t_y_yy_xz_xy,\
                           t_y_yy_xz_xz, t_y_yy_xz_yy, t_y_yy_xz_yz, t_y_yy_xz_zz, t_y_yy_yy_xx,\
                           t_y_yy_yy_xy, t_y_yy_yy_xz, t_y_yy_yy_yy, t_y_yy_yy_yz, t_y_yy_yy_zz,\
                           t_y_yy_yz_xx, t_y_yy_yz_xy, t_y_yy_yz_xz, t_y_yy_yz_yy, t_y_yy_yz_yz,\
                           t_y_yy_yz_zz, t_y_yy_zz_xx, t_y_yy_zz_xy, t_y_yy_zz_xz, t_y_yy_zz_yy,\
                           t_y_yy_zz_yz, t_y_yy_zz_zz, t_y_yyz_xx_xx, t_y_yyz_xx_xy,\
                           t_y_yyz_xx_xz, t_y_yyz_xx_yy, t_y_yyz_xx_yz, t_y_yyz_xx_zz,\
                           t_y_yyz_xy_xx, t_y_yyz_xy_xy, t_y_yyz_xy_xz, t_y_yyz_xy_yy,\
                           t_y_yyz_xy_yz, t_y_yyz_xy_zz, t_y_yyz_xz_xx, t_y_yyz_xz_xy,\
                           t_y_yyz_xz_xz, t_y_yyz_xz_yy, t_y_yyz_xz_yz, t_y_yyz_xz_zz,\
                           t_y_yyz_yy_xx, t_y_yyz_yy_xy, t_y_yyz_yy_xz, t_y_yyz_yy_yy,\
                           t_y_yyz_yy_yz, t_y_yyz_yy_zz, t_y_yyz_yz_xx, t_y_yyz_yz_xy,\
                           t_y_yyz_yz_xz, t_y_yyz_yz_yy, t_y_yyz_yz_yz, t_y_yyz_yz_zz,\
                           t_y_yyz_zz_xx, t_y_yyz_zz_xy, t_y_yyz_zz_xz, t_y_yyz_zz_yy,\
                           t_y_yyz_zz_yz, t_y_yyz_zz_zz, t_yz_yy_xx_xx, t_yz_yy_xx_xy,\
                           t_yz_yy_xx_xz, t_yz_yy_xx_yy, t_yz_yy_xx_yz, t_yz_yy_xx_zz,\
                           t_yz_yy_xy_xx, t_yz_yy_xy_xy, t_yz_yy_xy_xz, t_yz_yy_xy_yy,\
                           t_yz_yy_xy_yz, t_yz_yy_xy_zz, t_yz_yy_xz_xx, t_yz_yy_xz_xy,\
                           t_yz_yy_xz_xz, t_yz_yy_xz_yy, t_yz_yy_xz_yz, t_yz_yy_xz_zz,\
                           t_yz_yy_yy_xx, t_yz_yy_yy_xy, t_yz_yy_yy_xz, t_yz_yy_yy_yy,\
                           t_yz_yy_yy_yz, t_yz_yy_yy_zz, t_yz_yy_yz_xx, t_yz_yy_yz_xy,\
                           t_yz_yy_yz_xz, t_yz_yy_yz_yy, t_yz_yy_yz_yz, t_yz_yy_yz_zz,\
                           t_yz_yy_zz_xx, t_yz_yy_zz_xy, t_yz_yy_zz_xz, t_yz_yy_zz_yy,\
                           t_yz_yy_zz_yz, t_yz_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_yy_zz_zz[i] = t_y_yyz_zz_zz[i] - rab_z[i] * t_y_yy_zz_zz[i];

        t_yz_yy_zz_yz[i] = t_y_yyz_zz_yz[i] - rab_z[i] * t_y_yy_zz_yz[i];

        t_yz_yy_zz_yy[i] = t_y_yyz_zz_yy[i] - rab_z[i] * t_y_yy_zz_yy[i];

        t_yz_yy_zz_xz[i] = t_y_yyz_zz_xz[i] - rab_z[i] * t_y_yy_zz_xz[i];

        t_yz_yy_zz_xy[i] = t_y_yyz_zz_xy[i] - rab_z[i] * t_y_yy_zz_xy[i];

        t_yz_yy_zz_xx[i] = t_y_yyz_zz_xx[i] - rab_z[i] * t_y_yy_zz_xx[i];

        t_yz_yy_yz_zz[i] = t_y_yyz_yz_zz[i] - rab_z[i] * t_y_yy_yz_zz[i];

        t_yz_yy_yz_yz[i] = t_y_yyz_yz_yz[i] - rab_z[i] * t_y_yy_yz_yz[i];

        t_yz_yy_yz_yy[i] = t_y_yyz_yz_yy[i] - rab_z[i] * t_y_yy_yz_yy[i];

        t_yz_yy_yz_xz[i] = t_y_yyz_yz_xz[i] - rab_z[i] * t_y_yy_yz_xz[i];

        t_yz_yy_yz_xy[i] = t_y_yyz_yz_xy[i] - rab_z[i] * t_y_yy_yz_xy[i];

        t_yz_yy_yz_xx[i] = t_y_yyz_yz_xx[i] - rab_z[i] * t_y_yy_yz_xx[i];

        t_yz_yy_yy_zz[i] = t_y_yyz_yy_zz[i] - rab_z[i] * t_y_yy_yy_zz[i];

        t_yz_yy_yy_yz[i] = t_y_yyz_yy_yz[i] - rab_z[i] * t_y_yy_yy_yz[i];

        t_yz_yy_yy_yy[i] = t_y_yyz_yy_yy[i] - rab_z[i] * t_y_yy_yy_yy[i];

        t_yz_yy_yy_xz[i] = t_y_yyz_yy_xz[i] - rab_z[i] * t_y_yy_yy_xz[i];

        t_yz_yy_yy_xy[i] = t_y_yyz_yy_xy[i] - rab_z[i] * t_y_yy_yy_xy[i];

        t_yz_yy_yy_xx[i] = t_y_yyz_yy_xx[i] - rab_z[i] * t_y_yy_yy_xx[i];

        t_yz_yy_xz_zz[i] = t_y_yyz_xz_zz[i] - rab_z[i] * t_y_yy_xz_zz[i];

        t_yz_yy_xz_yz[i] = t_y_yyz_xz_yz[i] - rab_z[i] * t_y_yy_xz_yz[i];

        t_yz_yy_xz_yy[i] = t_y_yyz_xz_yy[i] - rab_z[i] * t_y_yy_xz_yy[i];

        t_yz_yy_xz_xz[i] = t_y_yyz_xz_xz[i] - rab_z[i] * t_y_yy_xz_xz[i];

        t_yz_yy_xz_xy[i] = t_y_yyz_xz_xy[i] - rab_z[i] * t_y_yy_xz_xy[i];

        t_yz_yy_xz_xx[i] = t_y_yyz_xz_xx[i] - rab_z[i] * t_y_yy_xz_xx[i];

        t_yz_yy_xy_zz[i] = t_y_yyz_xy_zz[i] - rab_z[i] * t_y_yy_xy_zz[i];

        t_yz_yy_xy_yz[i] = t_y_yyz_xy_yz[i] - rab_z[i] * t_y_yy_xy_yz[i];

        t_yz_yy_xy_yy[i] = t_y_yyz_xy_yy[i] - rab_z[i] * t_y_yy_xy_yy[i];

        t_yz_yy_xy_xz[i] = t_y_yyz_xy_xz[i] - rab_z[i] * t_y_yy_xy_xz[i];

        t_yz_yy_xy_xy[i] = t_y_yyz_xy_xy[i] - rab_z[i] * t_y_yy_xy_xy[i];

        t_yz_yy_xy_xx[i] = t_y_yyz_xy_xx[i] - rab_z[i] * t_y_yy_xy_xx[i];

        t_yz_yy_xx_zz[i] = t_y_yyz_xx_zz[i] - rab_z[i] * t_y_yy_xx_zz[i];

        t_yz_yy_xx_yz[i] = t_y_yyz_xx_yz[i] - rab_z[i] * t_y_yy_xx_yz[i];

        t_yz_yy_xx_yy[i] = t_y_yyz_xx_yy[i] - rab_z[i] * t_y_yy_xx_yy[i];

        t_yz_yy_xx_xz[i] = t_y_yyz_xx_xz[i] - rab_z[i] * t_y_yy_xx_xz[i];

        t_yz_yy_xx_xy[i] = t_y_yyz_xx_xy[i] - rab_z[i] * t_y_yy_xx_xy[i];

        t_yz_yy_xx_xx[i] = t_y_yyz_xx_xx[i] - rab_z[i] * t_y_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xz_xx_xx, t_y_xz_xx_xy, t_y_xz_xx_xz, t_y_xz_xx_yy,\
                           t_y_xz_xx_yz, t_y_xz_xx_zz, t_y_xz_xy_xx, t_y_xz_xy_xy, t_y_xz_xy_xz,\
                           t_y_xz_xy_yy, t_y_xz_xy_yz, t_y_xz_xy_zz, t_y_xz_xz_xx, t_y_xz_xz_xy,\
                           t_y_xz_xz_xz, t_y_xz_xz_yy, t_y_xz_xz_yz, t_y_xz_xz_zz, t_y_xz_yy_xx,\
                           t_y_xz_yy_xy, t_y_xz_yy_xz, t_y_xz_yy_yy, t_y_xz_yy_yz, t_y_xz_yy_zz,\
                           t_y_xz_yz_xx, t_y_xz_yz_xy, t_y_xz_yz_xz, t_y_xz_yz_yy, t_y_xz_yz_yz,\
                           t_y_xz_yz_zz, t_y_xz_zz_xx, t_y_xz_zz_xy, t_y_xz_zz_xz, t_y_xz_zz_yy,\
                           t_y_xz_zz_yz, t_y_xz_zz_zz, t_y_xzz_xx_xx, t_y_xzz_xx_xy,\
                           t_y_xzz_xx_xz, t_y_xzz_xx_yy, t_y_xzz_xx_yz, t_y_xzz_xx_zz,\
                           t_y_xzz_xy_xx, t_y_xzz_xy_xy, t_y_xzz_xy_xz, t_y_xzz_xy_yy,\
                           t_y_xzz_xy_yz, t_y_xzz_xy_zz, t_y_xzz_xz_xx, t_y_xzz_xz_xy,\
                           t_y_xzz_xz_xz, t_y_xzz_xz_yy, t_y_xzz_xz_yz, t_y_xzz_xz_zz,\
                           t_y_xzz_yy_xx, t_y_xzz_yy_xy, t_y_xzz_yy_xz, t_y_xzz_yy_yy,\
                           t_y_xzz_yy_yz, t_y_xzz_yy_zz, t_y_xzz_yz_xx, t_y_xzz_yz_xy,\
                           t_y_xzz_yz_xz, t_y_xzz_yz_yy, t_y_xzz_yz_yz, t_y_xzz_yz_zz,\
                           t_y_xzz_zz_xx, t_y_xzz_zz_xy, t_y_xzz_zz_xz, t_y_xzz_zz_yy,\
                           t_y_xzz_zz_yz, t_y_xzz_zz_zz, t_yz_xz_xx_xx, t_yz_xz_xx_xy,\
                           t_yz_xz_xx_xz, t_yz_xz_xx_yy, t_yz_xz_xx_yz, t_yz_xz_xx_zz,\
                           t_yz_xz_xy_xx, t_yz_xz_xy_xy, t_yz_xz_xy_xz, t_yz_xz_xy_yy,\
                           t_yz_xz_xy_yz, t_yz_xz_xy_zz, t_yz_xz_xz_xx, t_yz_xz_xz_xy,\
                           t_yz_xz_xz_xz, t_yz_xz_xz_yy, t_yz_xz_xz_yz, t_yz_xz_xz_zz,\
                           t_yz_xz_yy_xx, t_yz_xz_yy_xy, t_yz_xz_yy_xz, t_yz_xz_yy_yy,\
                           t_yz_xz_yy_yz, t_yz_xz_yy_zz, t_yz_xz_yz_xx, t_yz_xz_yz_xy,\
                           t_yz_xz_yz_xz, t_yz_xz_yz_yy, t_yz_xz_yz_yz, t_yz_xz_yz_zz,\
                           t_yz_xz_zz_xx, t_yz_xz_zz_xy, t_yz_xz_zz_xz, t_yz_xz_zz_yy,\
                           t_yz_xz_zz_yz, t_yz_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_xz_zz_zz[i] = t_y_xzz_zz_zz[i] - rab_z[i] * t_y_xz_zz_zz[i];

        t_yz_xz_zz_yz[i] = t_y_xzz_zz_yz[i] - rab_z[i] * t_y_xz_zz_yz[i];

        t_yz_xz_zz_yy[i] = t_y_xzz_zz_yy[i] - rab_z[i] * t_y_xz_zz_yy[i];

        t_yz_xz_zz_xz[i] = t_y_xzz_zz_xz[i] - rab_z[i] * t_y_xz_zz_xz[i];

        t_yz_xz_zz_xy[i] = t_y_xzz_zz_xy[i] - rab_z[i] * t_y_xz_zz_xy[i];

        t_yz_xz_zz_xx[i] = t_y_xzz_zz_xx[i] - rab_z[i] * t_y_xz_zz_xx[i];

        t_yz_xz_yz_zz[i] = t_y_xzz_yz_zz[i] - rab_z[i] * t_y_xz_yz_zz[i];

        t_yz_xz_yz_yz[i] = t_y_xzz_yz_yz[i] - rab_z[i] * t_y_xz_yz_yz[i];

        t_yz_xz_yz_yy[i] = t_y_xzz_yz_yy[i] - rab_z[i] * t_y_xz_yz_yy[i];

        t_yz_xz_yz_xz[i] = t_y_xzz_yz_xz[i] - rab_z[i] * t_y_xz_yz_xz[i];

        t_yz_xz_yz_xy[i] = t_y_xzz_yz_xy[i] - rab_z[i] * t_y_xz_yz_xy[i];

        t_yz_xz_yz_xx[i] = t_y_xzz_yz_xx[i] - rab_z[i] * t_y_xz_yz_xx[i];

        t_yz_xz_yy_zz[i] = t_y_xzz_yy_zz[i] - rab_z[i] * t_y_xz_yy_zz[i];

        t_yz_xz_yy_yz[i] = t_y_xzz_yy_yz[i] - rab_z[i] * t_y_xz_yy_yz[i];

        t_yz_xz_yy_yy[i] = t_y_xzz_yy_yy[i] - rab_z[i] * t_y_xz_yy_yy[i];

        t_yz_xz_yy_xz[i] = t_y_xzz_yy_xz[i] - rab_z[i] * t_y_xz_yy_xz[i];

        t_yz_xz_yy_xy[i] = t_y_xzz_yy_xy[i] - rab_z[i] * t_y_xz_yy_xy[i];

        t_yz_xz_yy_xx[i] = t_y_xzz_yy_xx[i] - rab_z[i] * t_y_xz_yy_xx[i];

        t_yz_xz_xz_zz[i] = t_y_xzz_xz_zz[i] - rab_z[i] * t_y_xz_xz_zz[i];

        t_yz_xz_xz_yz[i] = t_y_xzz_xz_yz[i] - rab_z[i] * t_y_xz_xz_yz[i];

        t_yz_xz_xz_yy[i] = t_y_xzz_xz_yy[i] - rab_z[i] * t_y_xz_xz_yy[i];

        t_yz_xz_xz_xz[i] = t_y_xzz_xz_xz[i] - rab_z[i] * t_y_xz_xz_xz[i];

        t_yz_xz_xz_xy[i] = t_y_xzz_xz_xy[i] - rab_z[i] * t_y_xz_xz_xy[i];

        t_yz_xz_xz_xx[i] = t_y_xzz_xz_xx[i] - rab_z[i] * t_y_xz_xz_xx[i];

        t_yz_xz_xy_zz[i] = t_y_xzz_xy_zz[i] - rab_z[i] * t_y_xz_xy_zz[i];

        t_yz_xz_xy_yz[i] = t_y_xzz_xy_yz[i] - rab_z[i] * t_y_xz_xy_yz[i];

        t_yz_xz_xy_yy[i] = t_y_xzz_xy_yy[i] - rab_z[i] * t_y_xz_xy_yy[i];

        t_yz_xz_xy_xz[i] = t_y_xzz_xy_xz[i] - rab_z[i] * t_y_xz_xy_xz[i];

        t_yz_xz_xy_xy[i] = t_y_xzz_xy_xy[i] - rab_z[i] * t_y_xz_xy_xy[i];

        t_yz_xz_xy_xx[i] = t_y_xzz_xy_xx[i] - rab_z[i] * t_y_xz_xy_xx[i];

        t_yz_xz_xx_zz[i] = t_y_xzz_xx_zz[i] - rab_z[i] * t_y_xz_xx_zz[i];

        t_yz_xz_xx_yz[i] = t_y_xzz_xx_yz[i] - rab_z[i] * t_y_xz_xx_yz[i];

        t_yz_xz_xx_yy[i] = t_y_xzz_xx_yy[i] - rab_z[i] * t_y_xz_xx_yy[i];

        t_yz_xz_xx_xz[i] = t_y_xzz_xx_xz[i] - rab_z[i] * t_y_xz_xx_xz[i];

        t_yz_xz_xx_xy[i] = t_y_xzz_xx_xy[i] - rab_z[i] * t_y_xz_xx_xy[i];

        t_yz_xz_xx_xx[i] = t_y_xzz_xx_xx[i] - rab_z[i] * t_y_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xy_xx_xx, t_y_xy_xx_xy, t_y_xy_xx_xz, t_y_xy_xx_yy,\
                           t_y_xy_xx_yz, t_y_xy_xx_zz, t_y_xy_xy_xx, t_y_xy_xy_xy, t_y_xy_xy_xz,\
                           t_y_xy_xy_yy, t_y_xy_xy_yz, t_y_xy_xy_zz, t_y_xy_xz_xx, t_y_xy_xz_xy,\
                           t_y_xy_xz_xz, t_y_xy_xz_yy, t_y_xy_xz_yz, t_y_xy_xz_zz, t_y_xy_yy_xx,\
                           t_y_xy_yy_xy, t_y_xy_yy_xz, t_y_xy_yy_yy, t_y_xy_yy_yz, t_y_xy_yy_zz,\
                           t_y_xy_yz_xx, t_y_xy_yz_xy, t_y_xy_yz_xz, t_y_xy_yz_yy, t_y_xy_yz_yz,\
                           t_y_xy_yz_zz, t_y_xy_zz_xx, t_y_xy_zz_xy, t_y_xy_zz_xz, t_y_xy_zz_yy,\
                           t_y_xy_zz_yz, t_y_xy_zz_zz, t_y_xyz_xx_xx, t_y_xyz_xx_xy,\
                           t_y_xyz_xx_xz, t_y_xyz_xx_yy, t_y_xyz_xx_yz, t_y_xyz_xx_zz,\
                           t_y_xyz_xy_xx, t_y_xyz_xy_xy, t_y_xyz_xy_xz, t_y_xyz_xy_yy,\
                           t_y_xyz_xy_yz, t_y_xyz_xy_zz, t_y_xyz_xz_xx, t_y_xyz_xz_xy,\
                           t_y_xyz_xz_xz, t_y_xyz_xz_yy, t_y_xyz_xz_yz, t_y_xyz_xz_zz,\
                           t_y_xyz_yy_xx, t_y_xyz_yy_xy, t_y_xyz_yy_xz, t_y_xyz_yy_yy,\
                           t_y_xyz_yy_yz, t_y_xyz_yy_zz, t_y_xyz_yz_xx, t_y_xyz_yz_xy,\
                           t_y_xyz_yz_xz, t_y_xyz_yz_yy, t_y_xyz_yz_yz, t_y_xyz_yz_zz,\
                           t_y_xyz_zz_xx, t_y_xyz_zz_xy, t_y_xyz_zz_xz, t_y_xyz_zz_yy,\
                           t_y_xyz_zz_yz, t_y_xyz_zz_zz, t_yz_xy_xx_xx, t_yz_xy_xx_xy,\
                           t_yz_xy_xx_xz, t_yz_xy_xx_yy, t_yz_xy_xx_yz, t_yz_xy_xx_zz,\
                           t_yz_xy_xy_xx, t_yz_xy_xy_xy, t_yz_xy_xy_xz, t_yz_xy_xy_yy,\
                           t_yz_xy_xy_yz, t_yz_xy_xy_zz, t_yz_xy_xz_xx, t_yz_xy_xz_xy,\
                           t_yz_xy_xz_xz, t_yz_xy_xz_yy, t_yz_xy_xz_yz, t_yz_xy_xz_zz,\
                           t_yz_xy_yy_xx, t_yz_xy_yy_xy, t_yz_xy_yy_xz, t_yz_xy_yy_yy,\
                           t_yz_xy_yy_yz, t_yz_xy_yy_zz, t_yz_xy_yz_xx, t_yz_xy_yz_xy,\
                           t_yz_xy_yz_xz, t_yz_xy_yz_yy, t_yz_xy_yz_yz, t_yz_xy_yz_zz,\
                           t_yz_xy_zz_xx, t_yz_xy_zz_xy, t_yz_xy_zz_xz, t_yz_xy_zz_yy,\
                           t_yz_xy_zz_yz, t_yz_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_xy_zz_zz[i] = t_y_xyz_zz_zz[i] - rab_z[i] * t_y_xy_zz_zz[i];

        t_yz_xy_zz_yz[i] = t_y_xyz_zz_yz[i] - rab_z[i] * t_y_xy_zz_yz[i];

        t_yz_xy_zz_yy[i] = t_y_xyz_zz_yy[i] - rab_z[i] * t_y_xy_zz_yy[i];

        t_yz_xy_zz_xz[i] = t_y_xyz_zz_xz[i] - rab_z[i] * t_y_xy_zz_xz[i];

        t_yz_xy_zz_xy[i] = t_y_xyz_zz_xy[i] - rab_z[i] * t_y_xy_zz_xy[i];

        t_yz_xy_zz_xx[i] = t_y_xyz_zz_xx[i] - rab_z[i] * t_y_xy_zz_xx[i];

        t_yz_xy_yz_zz[i] = t_y_xyz_yz_zz[i] - rab_z[i] * t_y_xy_yz_zz[i];

        t_yz_xy_yz_yz[i] = t_y_xyz_yz_yz[i] - rab_z[i] * t_y_xy_yz_yz[i];

        t_yz_xy_yz_yy[i] = t_y_xyz_yz_yy[i] - rab_z[i] * t_y_xy_yz_yy[i];

        t_yz_xy_yz_xz[i] = t_y_xyz_yz_xz[i] - rab_z[i] * t_y_xy_yz_xz[i];

        t_yz_xy_yz_xy[i] = t_y_xyz_yz_xy[i] - rab_z[i] * t_y_xy_yz_xy[i];

        t_yz_xy_yz_xx[i] = t_y_xyz_yz_xx[i] - rab_z[i] * t_y_xy_yz_xx[i];

        t_yz_xy_yy_zz[i] = t_y_xyz_yy_zz[i] - rab_z[i] * t_y_xy_yy_zz[i];

        t_yz_xy_yy_yz[i] = t_y_xyz_yy_yz[i] - rab_z[i] * t_y_xy_yy_yz[i];

        t_yz_xy_yy_yy[i] = t_y_xyz_yy_yy[i] - rab_z[i] * t_y_xy_yy_yy[i];

        t_yz_xy_yy_xz[i] = t_y_xyz_yy_xz[i] - rab_z[i] * t_y_xy_yy_xz[i];

        t_yz_xy_yy_xy[i] = t_y_xyz_yy_xy[i] - rab_z[i] * t_y_xy_yy_xy[i];

        t_yz_xy_yy_xx[i] = t_y_xyz_yy_xx[i] - rab_z[i] * t_y_xy_yy_xx[i];

        t_yz_xy_xz_zz[i] = t_y_xyz_xz_zz[i] - rab_z[i] * t_y_xy_xz_zz[i];

        t_yz_xy_xz_yz[i] = t_y_xyz_xz_yz[i] - rab_z[i] * t_y_xy_xz_yz[i];

        t_yz_xy_xz_yy[i] = t_y_xyz_xz_yy[i] - rab_z[i] * t_y_xy_xz_yy[i];

        t_yz_xy_xz_xz[i] = t_y_xyz_xz_xz[i] - rab_z[i] * t_y_xy_xz_xz[i];

        t_yz_xy_xz_xy[i] = t_y_xyz_xz_xy[i] - rab_z[i] * t_y_xy_xz_xy[i];

        t_yz_xy_xz_xx[i] = t_y_xyz_xz_xx[i] - rab_z[i] * t_y_xy_xz_xx[i];

        t_yz_xy_xy_zz[i] = t_y_xyz_xy_zz[i] - rab_z[i] * t_y_xy_xy_zz[i];

        t_yz_xy_xy_yz[i] = t_y_xyz_xy_yz[i] - rab_z[i] * t_y_xy_xy_yz[i];

        t_yz_xy_xy_yy[i] = t_y_xyz_xy_yy[i] - rab_z[i] * t_y_xy_xy_yy[i];

        t_yz_xy_xy_xz[i] = t_y_xyz_xy_xz[i] - rab_z[i] * t_y_xy_xy_xz[i];

        t_yz_xy_xy_xy[i] = t_y_xyz_xy_xy[i] - rab_z[i] * t_y_xy_xy_xy[i];

        t_yz_xy_xy_xx[i] = t_y_xyz_xy_xx[i] - rab_z[i] * t_y_xy_xy_xx[i];

        t_yz_xy_xx_zz[i] = t_y_xyz_xx_zz[i] - rab_z[i] * t_y_xy_xx_zz[i];

        t_yz_xy_xx_yz[i] = t_y_xyz_xx_yz[i] - rab_z[i] * t_y_xy_xx_yz[i];

        t_yz_xy_xx_yy[i] = t_y_xyz_xx_yy[i] - rab_z[i] * t_y_xy_xx_yy[i];

        t_yz_xy_xx_xz[i] = t_y_xyz_xx_xz[i] - rab_z[i] * t_y_xy_xx_xz[i];

        t_yz_xy_xx_xy[i] = t_y_xyz_xx_xy[i] - rab_z[i] * t_y_xy_xx_xy[i];

        t_yz_xy_xx_xx[i] = t_y_xyz_xx_xx[i] - rab_z[i] * t_y_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xx_xx_xx, t_y_xx_xx_xy, t_y_xx_xx_xz, t_y_xx_xx_yy,\
                           t_y_xx_xx_yz, t_y_xx_xx_zz, t_y_xx_xy_xx, t_y_xx_xy_xy, t_y_xx_xy_xz,\
                           t_y_xx_xy_yy, t_y_xx_xy_yz, t_y_xx_xy_zz, t_y_xx_xz_xx, t_y_xx_xz_xy,\
                           t_y_xx_xz_xz, t_y_xx_xz_yy, t_y_xx_xz_yz, t_y_xx_xz_zz, t_y_xx_yy_xx,\
                           t_y_xx_yy_xy, t_y_xx_yy_xz, t_y_xx_yy_yy, t_y_xx_yy_yz, t_y_xx_yy_zz,\
                           t_y_xx_yz_xx, t_y_xx_yz_xy, t_y_xx_yz_xz, t_y_xx_yz_yy, t_y_xx_yz_yz,\
                           t_y_xx_yz_zz, t_y_xx_zz_xx, t_y_xx_zz_xy, t_y_xx_zz_xz, t_y_xx_zz_yy,\
                           t_y_xx_zz_yz, t_y_xx_zz_zz, t_y_xxz_xx_xx, t_y_xxz_xx_xy,\
                           t_y_xxz_xx_xz, t_y_xxz_xx_yy, t_y_xxz_xx_yz, t_y_xxz_xx_zz,\
                           t_y_xxz_xy_xx, t_y_xxz_xy_xy, t_y_xxz_xy_xz, t_y_xxz_xy_yy,\
                           t_y_xxz_xy_yz, t_y_xxz_xy_zz, t_y_xxz_xz_xx, t_y_xxz_xz_xy,\
                           t_y_xxz_xz_xz, t_y_xxz_xz_yy, t_y_xxz_xz_yz, t_y_xxz_xz_zz,\
                           t_y_xxz_yy_xx, t_y_xxz_yy_xy, t_y_xxz_yy_xz, t_y_xxz_yy_yy,\
                           t_y_xxz_yy_yz, t_y_xxz_yy_zz, t_y_xxz_yz_xx, t_y_xxz_yz_xy,\
                           t_y_xxz_yz_xz, t_y_xxz_yz_yy, t_y_xxz_yz_yz, t_y_xxz_yz_zz,\
                           t_y_xxz_zz_xx, t_y_xxz_zz_xy, t_y_xxz_zz_xz, t_y_xxz_zz_yy,\
                           t_y_xxz_zz_yz, t_y_xxz_zz_zz, t_yz_xx_xx_xx, t_yz_xx_xx_xy,\
                           t_yz_xx_xx_xz, t_yz_xx_xx_yy, t_yz_xx_xx_yz, t_yz_xx_xx_zz,\
                           t_yz_xx_xy_xx, t_yz_xx_xy_xy, t_yz_xx_xy_xz, t_yz_xx_xy_yy,\
                           t_yz_xx_xy_yz, t_yz_xx_xy_zz, t_yz_xx_xz_xx, t_yz_xx_xz_xy,\
                           t_yz_xx_xz_xz, t_yz_xx_xz_yy, t_yz_xx_xz_yz, t_yz_xx_xz_zz,\
                           t_yz_xx_yy_xx, t_yz_xx_yy_xy, t_yz_xx_yy_xz, t_yz_xx_yy_yy,\
                           t_yz_xx_yy_yz, t_yz_xx_yy_zz, t_yz_xx_yz_xx, t_yz_xx_yz_xy,\
                           t_yz_xx_yz_xz, t_yz_xx_yz_yy, t_yz_xx_yz_yz, t_yz_xx_yz_zz,\
                           t_yz_xx_zz_xx, t_yz_xx_zz_xy, t_yz_xx_zz_xz, t_yz_xx_zz_yy,\
                           t_yz_xx_zz_yz, t_yz_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_xx_zz_zz[i] = t_y_xxz_zz_zz[i] - rab_z[i] * t_y_xx_zz_zz[i];

        t_yz_xx_zz_yz[i] = t_y_xxz_zz_yz[i] - rab_z[i] * t_y_xx_zz_yz[i];

        t_yz_xx_zz_yy[i] = t_y_xxz_zz_yy[i] - rab_z[i] * t_y_xx_zz_yy[i];

        t_yz_xx_zz_xz[i] = t_y_xxz_zz_xz[i] - rab_z[i] * t_y_xx_zz_xz[i];

        t_yz_xx_zz_xy[i] = t_y_xxz_zz_xy[i] - rab_z[i] * t_y_xx_zz_xy[i];

        t_yz_xx_zz_xx[i] = t_y_xxz_zz_xx[i] - rab_z[i] * t_y_xx_zz_xx[i];

        t_yz_xx_yz_zz[i] = t_y_xxz_yz_zz[i] - rab_z[i] * t_y_xx_yz_zz[i];

        t_yz_xx_yz_yz[i] = t_y_xxz_yz_yz[i] - rab_z[i] * t_y_xx_yz_yz[i];

        t_yz_xx_yz_yy[i] = t_y_xxz_yz_yy[i] - rab_z[i] * t_y_xx_yz_yy[i];

        t_yz_xx_yz_xz[i] = t_y_xxz_yz_xz[i] - rab_z[i] * t_y_xx_yz_xz[i];

        t_yz_xx_yz_xy[i] = t_y_xxz_yz_xy[i] - rab_z[i] * t_y_xx_yz_xy[i];

        t_yz_xx_yz_xx[i] = t_y_xxz_yz_xx[i] - rab_z[i] * t_y_xx_yz_xx[i];

        t_yz_xx_yy_zz[i] = t_y_xxz_yy_zz[i] - rab_z[i] * t_y_xx_yy_zz[i];

        t_yz_xx_yy_yz[i] = t_y_xxz_yy_yz[i] - rab_z[i] * t_y_xx_yy_yz[i];

        t_yz_xx_yy_yy[i] = t_y_xxz_yy_yy[i] - rab_z[i] * t_y_xx_yy_yy[i];

        t_yz_xx_yy_xz[i] = t_y_xxz_yy_xz[i] - rab_z[i] * t_y_xx_yy_xz[i];

        t_yz_xx_yy_xy[i] = t_y_xxz_yy_xy[i] - rab_z[i] * t_y_xx_yy_xy[i];

        t_yz_xx_yy_xx[i] = t_y_xxz_yy_xx[i] - rab_z[i] * t_y_xx_yy_xx[i];

        t_yz_xx_xz_zz[i] = t_y_xxz_xz_zz[i] - rab_z[i] * t_y_xx_xz_zz[i];

        t_yz_xx_xz_yz[i] = t_y_xxz_xz_yz[i] - rab_z[i] * t_y_xx_xz_yz[i];

        t_yz_xx_xz_yy[i] = t_y_xxz_xz_yy[i] - rab_z[i] * t_y_xx_xz_yy[i];

        t_yz_xx_xz_xz[i] = t_y_xxz_xz_xz[i] - rab_z[i] * t_y_xx_xz_xz[i];

        t_yz_xx_xz_xy[i] = t_y_xxz_xz_xy[i] - rab_z[i] * t_y_xx_xz_xy[i];

        t_yz_xx_xz_xx[i] = t_y_xxz_xz_xx[i] - rab_z[i] * t_y_xx_xz_xx[i];

        t_yz_xx_xy_zz[i] = t_y_xxz_xy_zz[i] - rab_z[i] * t_y_xx_xy_zz[i];

        t_yz_xx_xy_yz[i] = t_y_xxz_xy_yz[i] - rab_z[i] * t_y_xx_xy_yz[i];

        t_yz_xx_xy_yy[i] = t_y_xxz_xy_yy[i] - rab_z[i] * t_y_xx_xy_yy[i];

        t_yz_xx_xy_xz[i] = t_y_xxz_xy_xz[i] - rab_z[i] * t_y_xx_xy_xz[i];

        t_yz_xx_xy_xy[i] = t_y_xxz_xy_xy[i] - rab_z[i] * t_y_xx_xy_xy[i];

        t_yz_xx_xy_xx[i] = t_y_xxz_xy_xx[i] - rab_z[i] * t_y_xx_xy_xx[i];

        t_yz_xx_xx_zz[i] = t_y_xxz_xx_zz[i] - rab_z[i] * t_y_xx_xx_zz[i];

        t_yz_xx_xx_yz[i] = t_y_xxz_xx_yz[i] - rab_z[i] * t_y_xx_xx_yz[i];

        t_yz_xx_xx_yy[i] = t_y_xxz_xx_yy[i] - rab_z[i] * t_y_xx_xx_yy[i];

        t_yz_xx_xx_xz[i] = t_y_xxz_xx_xz[i] - rab_z[i] * t_y_xx_xx_xz[i];

        t_yz_xx_xx_xy[i] = t_y_xxz_xx_xy[i] - rab_z[i] * t_y_xx_xx_xy[i];

        t_yz_xx_xx_xx[i] = t_y_xxz_xx_xx[i] - rab_z[i] * t_y_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_yzz_xx_xx, t_y_yzz_xx_xy, t_y_yzz_xx_xz, t_y_yzz_xx_yy,\
                           t_y_yzz_xx_yz, t_y_yzz_xx_zz, t_y_yzz_xy_xx, t_y_yzz_xy_xy,\
                           t_y_yzz_xy_xz, t_y_yzz_xy_yy, t_y_yzz_xy_yz, t_y_yzz_xy_zz,\
                           t_y_yzz_xz_xx, t_y_yzz_xz_xy, t_y_yzz_xz_xz, t_y_yzz_xz_yy,\
                           t_y_yzz_xz_yz, t_y_yzz_xz_zz, t_y_yzz_yy_xx, t_y_yzz_yy_xy,\
                           t_y_yzz_yy_xz, t_y_yzz_yy_yy, t_y_yzz_yy_yz, t_y_yzz_yy_zz,\
                           t_y_yzz_yz_xx, t_y_yzz_yz_xy, t_y_yzz_yz_xz, t_y_yzz_yz_yy,\
                           t_y_yzz_yz_yz, t_y_yzz_yz_zz, t_y_yzz_zz_xx, t_y_yzz_zz_xy,\
                           t_y_yzz_zz_xz, t_y_yzz_zz_yy, t_y_yzz_zz_yz, t_y_yzz_zz_zz,\
                           t_y_zz_xx_xx, t_y_zz_xx_xy, t_y_zz_xx_xz, t_y_zz_xx_yy, t_y_zz_xx_yz,\
                           t_y_zz_xx_zz, t_y_zz_xy_xx, t_y_zz_xy_xy, t_y_zz_xy_xz, t_y_zz_xy_yy,\
                           t_y_zz_xy_yz, t_y_zz_xy_zz, t_y_zz_xz_xx, t_y_zz_xz_xy, t_y_zz_xz_xz,\
                           t_y_zz_xz_yy, t_y_zz_xz_yz, t_y_zz_xz_zz, t_y_zz_yy_xx, t_y_zz_yy_xy,\
                           t_y_zz_yy_xz, t_y_zz_yy_yy, t_y_zz_yy_yz, t_y_zz_yy_zz, t_y_zz_yz_xx,\
                           t_y_zz_yz_xy, t_y_zz_yz_xz, t_y_zz_yz_yy, t_y_zz_yz_yz, t_y_zz_yz_zz,\
                           t_y_zz_zz_xx, t_y_zz_zz_xy, t_y_zz_zz_xz, t_y_zz_zz_yy, t_y_zz_zz_yz,\
                           t_y_zz_zz_zz, t_yy_zz_xx_xx, t_yy_zz_xx_xy, t_yy_zz_xx_xz,\
                           t_yy_zz_xx_yy, t_yy_zz_xx_yz, t_yy_zz_xx_zz, t_yy_zz_xy_xx,\
                           t_yy_zz_xy_xy, t_yy_zz_xy_xz, t_yy_zz_xy_yy, t_yy_zz_xy_yz,\
                           t_yy_zz_xy_zz, t_yy_zz_xz_xx, t_yy_zz_xz_xy, t_yy_zz_xz_xz,\
                           t_yy_zz_xz_yy, t_yy_zz_xz_yz, t_yy_zz_xz_zz, t_yy_zz_yy_xx,\
                           t_yy_zz_yy_xy, t_yy_zz_yy_xz, t_yy_zz_yy_yy, t_yy_zz_yy_yz,\
                           t_yy_zz_yy_zz, t_yy_zz_yz_xx, t_yy_zz_yz_xy, t_yy_zz_yz_xz,\
                           t_yy_zz_yz_yy, t_yy_zz_yz_yz, t_yy_zz_yz_zz, t_yy_zz_zz_xx,\
                           t_yy_zz_zz_xy, t_yy_zz_zz_xz, t_yy_zz_zz_yy, t_yy_zz_zz_yz,\
                           t_yy_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_zz_zz_zz[i] = t_y_yzz_zz_zz[i] - rab_y[i] * t_y_zz_zz_zz[i];

        t_yy_zz_zz_yz[i] = t_y_yzz_zz_yz[i] - rab_y[i] * t_y_zz_zz_yz[i];

        t_yy_zz_zz_yy[i] = t_y_yzz_zz_yy[i] - rab_y[i] * t_y_zz_zz_yy[i];

        t_yy_zz_zz_xz[i] = t_y_yzz_zz_xz[i] - rab_y[i] * t_y_zz_zz_xz[i];

        t_yy_zz_zz_xy[i] = t_y_yzz_zz_xy[i] - rab_y[i] * t_y_zz_zz_xy[i];

        t_yy_zz_zz_xx[i] = t_y_yzz_zz_xx[i] - rab_y[i] * t_y_zz_zz_xx[i];

        t_yy_zz_yz_zz[i] = t_y_yzz_yz_zz[i] - rab_y[i] * t_y_zz_yz_zz[i];

        t_yy_zz_yz_yz[i] = t_y_yzz_yz_yz[i] - rab_y[i] * t_y_zz_yz_yz[i];

        t_yy_zz_yz_yy[i] = t_y_yzz_yz_yy[i] - rab_y[i] * t_y_zz_yz_yy[i];

        t_yy_zz_yz_xz[i] = t_y_yzz_yz_xz[i] - rab_y[i] * t_y_zz_yz_xz[i];

        t_yy_zz_yz_xy[i] = t_y_yzz_yz_xy[i] - rab_y[i] * t_y_zz_yz_xy[i];

        t_yy_zz_yz_xx[i] = t_y_yzz_yz_xx[i] - rab_y[i] * t_y_zz_yz_xx[i];

        t_yy_zz_yy_zz[i] = t_y_yzz_yy_zz[i] - rab_y[i] * t_y_zz_yy_zz[i];

        t_yy_zz_yy_yz[i] = t_y_yzz_yy_yz[i] - rab_y[i] * t_y_zz_yy_yz[i];

        t_yy_zz_yy_yy[i] = t_y_yzz_yy_yy[i] - rab_y[i] * t_y_zz_yy_yy[i];

        t_yy_zz_yy_xz[i] = t_y_yzz_yy_xz[i] - rab_y[i] * t_y_zz_yy_xz[i];

        t_yy_zz_yy_xy[i] = t_y_yzz_yy_xy[i] - rab_y[i] * t_y_zz_yy_xy[i];

        t_yy_zz_yy_xx[i] = t_y_yzz_yy_xx[i] - rab_y[i] * t_y_zz_yy_xx[i];

        t_yy_zz_xz_zz[i] = t_y_yzz_xz_zz[i] - rab_y[i] * t_y_zz_xz_zz[i];

        t_yy_zz_xz_yz[i] = t_y_yzz_xz_yz[i] - rab_y[i] * t_y_zz_xz_yz[i];

        t_yy_zz_xz_yy[i] = t_y_yzz_xz_yy[i] - rab_y[i] * t_y_zz_xz_yy[i];

        t_yy_zz_xz_xz[i] = t_y_yzz_xz_xz[i] - rab_y[i] * t_y_zz_xz_xz[i];

        t_yy_zz_xz_xy[i] = t_y_yzz_xz_xy[i] - rab_y[i] * t_y_zz_xz_xy[i];

        t_yy_zz_xz_xx[i] = t_y_yzz_xz_xx[i] - rab_y[i] * t_y_zz_xz_xx[i];

        t_yy_zz_xy_zz[i] = t_y_yzz_xy_zz[i] - rab_y[i] * t_y_zz_xy_zz[i];

        t_yy_zz_xy_yz[i] = t_y_yzz_xy_yz[i] - rab_y[i] * t_y_zz_xy_yz[i];

        t_yy_zz_xy_yy[i] = t_y_yzz_xy_yy[i] - rab_y[i] * t_y_zz_xy_yy[i];

        t_yy_zz_xy_xz[i] = t_y_yzz_xy_xz[i] - rab_y[i] * t_y_zz_xy_xz[i];

        t_yy_zz_xy_xy[i] = t_y_yzz_xy_xy[i] - rab_y[i] * t_y_zz_xy_xy[i];

        t_yy_zz_xy_xx[i] = t_y_yzz_xy_xx[i] - rab_y[i] * t_y_zz_xy_xx[i];

        t_yy_zz_xx_zz[i] = t_y_yzz_xx_zz[i] - rab_y[i] * t_y_zz_xx_zz[i];

        t_yy_zz_xx_yz[i] = t_y_yzz_xx_yz[i] - rab_y[i] * t_y_zz_xx_yz[i];

        t_yy_zz_xx_yy[i] = t_y_yzz_xx_yy[i] - rab_y[i] * t_y_zz_xx_yy[i];

        t_yy_zz_xx_xz[i] = t_y_yzz_xx_xz[i] - rab_y[i] * t_y_zz_xx_xz[i];

        t_yy_zz_xx_xy[i] = t_y_yzz_xx_xy[i] - rab_y[i] * t_y_zz_xx_xy[i];

        t_yy_zz_xx_xx[i] = t_y_yzz_xx_xx[i] - rab_y[i] * t_y_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_yyz_xx_xx, t_y_yyz_xx_xy, t_y_yyz_xx_xz, t_y_yyz_xx_yy,\
                           t_y_yyz_xx_yz, t_y_yyz_xx_zz, t_y_yyz_xy_xx, t_y_yyz_xy_xy,\
                           t_y_yyz_xy_xz, t_y_yyz_xy_yy, t_y_yyz_xy_yz, t_y_yyz_xy_zz,\
                           t_y_yyz_xz_xx, t_y_yyz_xz_xy, t_y_yyz_xz_xz, t_y_yyz_xz_yy,\
                           t_y_yyz_xz_yz, t_y_yyz_xz_zz, t_y_yyz_yy_xx, t_y_yyz_yy_xy,\
                           t_y_yyz_yy_xz, t_y_yyz_yy_yy, t_y_yyz_yy_yz, t_y_yyz_yy_zz,\
                           t_y_yyz_yz_xx, t_y_yyz_yz_xy, t_y_yyz_yz_xz, t_y_yyz_yz_yy,\
                           t_y_yyz_yz_yz, t_y_yyz_yz_zz, t_y_yyz_zz_xx, t_y_yyz_zz_xy,\
                           t_y_yyz_zz_xz, t_y_yyz_zz_yy, t_y_yyz_zz_yz, t_y_yyz_zz_zz,\
                           t_y_yz_xx_xx, t_y_yz_xx_xy, t_y_yz_xx_xz, t_y_yz_xx_yy, t_y_yz_xx_yz,\
                           t_y_yz_xx_zz, t_y_yz_xy_xx, t_y_yz_xy_xy, t_y_yz_xy_xz, t_y_yz_xy_yy,\
                           t_y_yz_xy_yz, t_y_yz_xy_zz, t_y_yz_xz_xx, t_y_yz_xz_xy, t_y_yz_xz_xz,\
                           t_y_yz_xz_yy, t_y_yz_xz_yz, t_y_yz_xz_zz, t_y_yz_yy_xx, t_y_yz_yy_xy,\
                           t_y_yz_yy_xz, t_y_yz_yy_yy, t_y_yz_yy_yz, t_y_yz_yy_zz, t_y_yz_yz_xx,\
                           t_y_yz_yz_xy, t_y_yz_yz_xz, t_y_yz_yz_yy, t_y_yz_yz_yz, t_y_yz_yz_zz,\
                           t_y_yz_zz_xx, t_y_yz_zz_xy, t_y_yz_zz_xz, t_y_yz_zz_yy, t_y_yz_zz_yz,\
                           t_y_yz_zz_zz, t_yy_yz_xx_xx, t_yy_yz_xx_xy, t_yy_yz_xx_xz,\
                           t_yy_yz_xx_yy, t_yy_yz_xx_yz, t_yy_yz_xx_zz, t_yy_yz_xy_xx,\
                           t_yy_yz_xy_xy, t_yy_yz_xy_xz, t_yy_yz_xy_yy, t_yy_yz_xy_yz,\
                           t_yy_yz_xy_zz, t_yy_yz_xz_xx, t_yy_yz_xz_xy, t_yy_yz_xz_xz,\
                           t_yy_yz_xz_yy, t_yy_yz_xz_yz, t_yy_yz_xz_zz, t_yy_yz_yy_xx,\
                           t_yy_yz_yy_xy, t_yy_yz_yy_xz, t_yy_yz_yy_yy, t_yy_yz_yy_yz,\
                           t_yy_yz_yy_zz, t_yy_yz_yz_xx, t_yy_yz_yz_xy, t_yy_yz_yz_xz,\
                           t_yy_yz_yz_yy, t_yy_yz_yz_yz, t_yy_yz_yz_zz, t_yy_yz_zz_xx,\
                           t_yy_yz_zz_xy, t_yy_yz_zz_xz, t_yy_yz_zz_yy, t_yy_yz_zz_yz,\
                           t_yy_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_yz_zz_zz[i] = t_y_yyz_zz_zz[i] - rab_y[i] * t_y_yz_zz_zz[i];

        t_yy_yz_zz_yz[i] = t_y_yyz_zz_yz[i] - rab_y[i] * t_y_yz_zz_yz[i];

        t_yy_yz_zz_yy[i] = t_y_yyz_zz_yy[i] - rab_y[i] * t_y_yz_zz_yy[i];

        t_yy_yz_zz_xz[i] = t_y_yyz_zz_xz[i] - rab_y[i] * t_y_yz_zz_xz[i];

        t_yy_yz_zz_xy[i] = t_y_yyz_zz_xy[i] - rab_y[i] * t_y_yz_zz_xy[i];

        t_yy_yz_zz_xx[i] = t_y_yyz_zz_xx[i] - rab_y[i] * t_y_yz_zz_xx[i];

        t_yy_yz_yz_zz[i] = t_y_yyz_yz_zz[i] - rab_y[i] * t_y_yz_yz_zz[i];

        t_yy_yz_yz_yz[i] = t_y_yyz_yz_yz[i] - rab_y[i] * t_y_yz_yz_yz[i];

        t_yy_yz_yz_yy[i] = t_y_yyz_yz_yy[i] - rab_y[i] * t_y_yz_yz_yy[i];

        t_yy_yz_yz_xz[i] = t_y_yyz_yz_xz[i] - rab_y[i] * t_y_yz_yz_xz[i];

        t_yy_yz_yz_xy[i] = t_y_yyz_yz_xy[i] - rab_y[i] * t_y_yz_yz_xy[i];

        t_yy_yz_yz_xx[i] = t_y_yyz_yz_xx[i] - rab_y[i] * t_y_yz_yz_xx[i];

        t_yy_yz_yy_zz[i] = t_y_yyz_yy_zz[i] - rab_y[i] * t_y_yz_yy_zz[i];

        t_yy_yz_yy_yz[i] = t_y_yyz_yy_yz[i] - rab_y[i] * t_y_yz_yy_yz[i];

        t_yy_yz_yy_yy[i] = t_y_yyz_yy_yy[i] - rab_y[i] * t_y_yz_yy_yy[i];

        t_yy_yz_yy_xz[i] = t_y_yyz_yy_xz[i] - rab_y[i] * t_y_yz_yy_xz[i];

        t_yy_yz_yy_xy[i] = t_y_yyz_yy_xy[i] - rab_y[i] * t_y_yz_yy_xy[i];

        t_yy_yz_yy_xx[i] = t_y_yyz_yy_xx[i] - rab_y[i] * t_y_yz_yy_xx[i];

        t_yy_yz_xz_zz[i] = t_y_yyz_xz_zz[i] - rab_y[i] * t_y_yz_xz_zz[i];

        t_yy_yz_xz_yz[i] = t_y_yyz_xz_yz[i] - rab_y[i] * t_y_yz_xz_yz[i];

        t_yy_yz_xz_yy[i] = t_y_yyz_xz_yy[i] - rab_y[i] * t_y_yz_xz_yy[i];

        t_yy_yz_xz_xz[i] = t_y_yyz_xz_xz[i] - rab_y[i] * t_y_yz_xz_xz[i];

        t_yy_yz_xz_xy[i] = t_y_yyz_xz_xy[i] - rab_y[i] * t_y_yz_xz_xy[i];

        t_yy_yz_xz_xx[i] = t_y_yyz_xz_xx[i] - rab_y[i] * t_y_yz_xz_xx[i];

        t_yy_yz_xy_zz[i] = t_y_yyz_xy_zz[i] - rab_y[i] * t_y_yz_xy_zz[i];

        t_yy_yz_xy_yz[i] = t_y_yyz_xy_yz[i] - rab_y[i] * t_y_yz_xy_yz[i];

        t_yy_yz_xy_yy[i] = t_y_yyz_xy_yy[i] - rab_y[i] * t_y_yz_xy_yy[i];

        t_yy_yz_xy_xz[i] = t_y_yyz_xy_xz[i] - rab_y[i] * t_y_yz_xy_xz[i];

        t_yy_yz_xy_xy[i] = t_y_yyz_xy_xy[i] - rab_y[i] * t_y_yz_xy_xy[i];

        t_yy_yz_xy_xx[i] = t_y_yyz_xy_xx[i] - rab_y[i] * t_y_yz_xy_xx[i];

        t_yy_yz_xx_zz[i] = t_y_yyz_xx_zz[i] - rab_y[i] * t_y_yz_xx_zz[i];

        t_yy_yz_xx_yz[i] = t_y_yyz_xx_yz[i] - rab_y[i] * t_y_yz_xx_yz[i];

        t_yy_yz_xx_yy[i] = t_y_yyz_xx_yy[i] - rab_y[i] * t_y_yz_xx_yy[i];

        t_yy_yz_xx_xz[i] = t_y_yyz_xx_xz[i] - rab_y[i] * t_y_yz_xx_xz[i];

        t_yy_yz_xx_xy[i] = t_y_yyz_xx_xy[i] - rab_y[i] * t_y_yz_xx_xy[i];

        t_yy_yz_xx_xx[i] = t_y_yyz_xx_xx[i] - rab_y[i] * t_y_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_yy_xx_xx, t_y_yy_xx_xy, t_y_yy_xx_xz, t_y_yy_xx_yy,\
                           t_y_yy_xx_yz, t_y_yy_xx_zz, t_y_yy_xy_xx, t_y_yy_xy_xy, t_y_yy_xy_xz,\
                           t_y_yy_xy_yy, t_y_yy_xy_yz, t_y_yy_xy_zz, t_y_yy_xz_xx, t_y_yy_xz_xy,\
                           t_y_yy_xz_xz, t_y_yy_xz_yy, t_y_yy_xz_yz, t_y_yy_xz_zz, t_y_yy_yy_xx,\
                           t_y_yy_yy_xy, t_y_yy_yy_xz, t_y_yy_yy_yy, t_y_yy_yy_yz, t_y_yy_yy_zz,\
                           t_y_yy_yz_xx, t_y_yy_yz_xy, t_y_yy_yz_xz, t_y_yy_yz_yy, t_y_yy_yz_yz,\
                           t_y_yy_yz_zz, t_y_yy_zz_xx, t_y_yy_zz_xy, t_y_yy_zz_xz, t_y_yy_zz_yy,\
                           t_y_yy_zz_yz, t_y_yy_zz_zz, t_y_yyy_xx_xx, t_y_yyy_xx_xy,\
                           t_y_yyy_xx_xz, t_y_yyy_xx_yy, t_y_yyy_xx_yz, t_y_yyy_xx_zz,\
                           t_y_yyy_xy_xx, t_y_yyy_xy_xy, t_y_yyy_xy_xz, t_y_yyy_xy_yy,\
                           t_y_yyy_xy_yz, t_y_yyy_xy_zz, t_y_yyy_xz_xx, t_y_yyy_xz_xy,\
                           t_y_yyy_xz_xz, t_y_yyy_xz_yy, t_y_yyy_xz_yz, t_y_yyy_xz_zz,\
                           t_y_yyy_yy_xx, t_y_yyy_yy_xy, t_y_yyy_yy_xz, t_y_yyy_yy_yy,\
                           t_y_yyy_yy_yz, t_y_yyy_yy_zz, t_y_yyy_yz_xx, t_y_yyy_yz_xy,\
                           t_y_yyy_yz_xz, t_y_yyy_yz_yy, t_y_yyy_yz_yz, t_y_yyy_yz_zz,\
                           t_y_yyy_zz_xx, t_y_yyy_zz_xy, t_y_yyy_zz_xz, t_y_yyy_zz_yy,\
                           t_y_yyy_zz_yz, t_y_yyy_zz_zz, t_yy_yy_xx_xx, t_yy_yy_xx_xy,\
                           t_yy_yy_xx_xz, t_yy_yy_xx_yy, t_yy_yy_xx_yz, t_yy_yy_xx_zz,\
                           t_yy_yy_xy_xx, t_yy_yy_xy_xy, t_yy_yy_xy_xz, t_yy_yy_xy_yy,\
                           t_yy_yy_xy_yz, t_yy_yy_xy_zz, t_yy_yy_xz_xx, t_yy_yy_xz_xy,\
                           t_yy_yy_xz_xz, t_yy_yy_xz_yy, t_yy_yy_xz_yz, t_yy_yy_xz_zz,\
                           t_yy_yy_yy_xx, t_yy_yy_yy_xy, t_yy_yy_yy_xz, t_yy_yy_yy_yy,\
                           t_yy_yy_yy_yz, t_yy_yy_yy_zz, t_yy_yy_yz_xx, t_yy_yy_yz_xy,\
                           t_yy_yy_yz_xz, t_yy_yy_yz_yy, t_yy_yy_yz_yz, t_yy_yy_yz_zz,\
                           t_yy_yy_zz_xx, t_yy_yy_zz_xy, t_yy_yy_zz_xz, t_yy_yy_zz_yy,\
                           t_yy_yy_zz_yz, t_yy_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_yy_zz_zz[i] = t_y_yyy_zz_zz[i] - rab_y[i] * t_y_yy_zz_zz[i];

        t_yy_yy_zz_yz[i] = t_y_yyy_zz_yz[i] - rab_y[i] * t_y_yy_zz_yz[i];

        t_yy_yy_zz_yy[i] = t_y_yyy_zz_yy[i] - rab_y[i] * t_y_yy_zz_yy[i];

        t_yy_yy_zz_xz[i] = t_y_yyy_zz_xz[i] - rab_y[i] * t_y_yy_zz_xz[i];

        t_yy_yy_zz_xy[i] = t_y_yyy_zz_xy[i] - rab_y[i] * t_y_yy_zz_xy[i];

        t_yy_yy_zz_xx[i] = t_y_yyy_zz_xx[i] - rab_y[i] * t_y_yy_zz_xx[i];

        t_yy_yy_yz_zz[i] = t_y_yyy_yz_zz[i] - rab_y[i] * t_y_yy_yz_zz[i];

        t_yy_yy_yz_yz[i] = t_y_yyy_yz_yz[i] - rab_y[i] * t_y_yy_yz_yz[i];

        t_yy_yy_yz_yy[i] = t_y_yyy_yz_yy[i] - rab_y[i] * t_y_yy_yz_yy[i];

        t_yy_yy_yz_xz[i] = t_y_yyy_yz_xz[i] - rab_y[i] * t_y_yy_yz_xz[i];

        t_yy_yy_yz_xy[i] = t_y_yyy_yz_xy[i] - rab_y[i] * t_y_yy_yz_xy[i];

        t_yy_yy_yz_xx[i] = t_y_yyy_yz_xx[i] - rab_y[i] * t_y_yy_yz_xx[i];

        t_yy_yy_yy_zz[i] = t_y_yyy_yy_zz[i] - rab_y[i] * t_y_yy_yy_zz[i];

        t_yy_yy_yy_yz[i] = t_y_yyy_yy_yz[i] - rab_y[i] * t_y_yy_yy_yz[i];

        t_yy_yy_yy_yy[i] = t_y_yyy_yy_yy[i] - rab_y[i] * t_y_yy_yy_yy[i];

        t_yy_yy_yy_xz[i] = t_y_yyy_yy_xz[i] - rab_y[i] * t_y_yy_yy_xz[i];

        t_yy_yy_yy_xy[i] = t_y_yyy_yy_xy[i] - rab_y[i] * t_y_yy_yy_xy[i];

        t_yy_yy_yy_xx[i] = t_y_yyy_yy_xx[i] - rab_y[i] * t_y_yy_yy_xx[i];

        t_yy_yy_xz_zz[i] = t_y_yyy_xz_zz[i] - rab_y[i] * t_y_yy_xz_zz[i];

        t_yy_yy_xz_yz[i] = t_y_yyy_xz_yz[i] - rab_y[i] * t_y_yy_xz_yz[i];

        t_yy_yy_xz_yy[i] = t_y_yyy_xz_yy[i] - rab_y[i] * t_y_yy_xz_yy[i];

        t_yy_yy_xz_xz[i] = t_y_yyy_xz_xz[i] - rab_y[i] * t_y_yy_xz_xz[i];

        t_yy_yy_xz_xy[i] = t_y_yyy_xz_xy[i] - rab_y[i] * t_y_yy_xz_xy[i];

        t_yy_yy_xz_xx[i] = t_y_yyy_xz_xx[i] - rab_y[i] * t_y_yy_xz_xx[i];

        t_yy_yy_xy_zz[i] = t_y_yyy_xy_zz[i] - rab_y[i] * t_y_yy_xy_zz[i];

        t_yy_yy_xy_yz[i] = t_y_yyy_xy_yz[i] - rab_y[i] * t_y_yy_xy_yz[i];

        t_yy_yy_xy_yy[i] = t_y_yyy_xy_yy[i] - rab_y[i] * t_y_yy_xy_yy[i];

        t_yy_yy_xy_xz[i] = t_y_yyy_xy_xz[i] - rab_y[i] * t_y_yy_xy_xz[i];

        t_yy_yy_xy_xy[i] = t_y_yyy_xy_xy[i] - rab_y[i] * t_y_yy_xy_xy[i];

        t_yy_yy_xy_xx[i] = t_y_yyy_xy_xx[i] - rab_y[i] * t_y_yy_xy_xx[i];

        t_yy_yy_xx_zz[i] = t_y_yyy_xx_zz[i] - rab_y[i] * t_y_yy_xx_zz[i];

        t_yy_yy_xx_yz[i] = t_y_yyy_xx_yz[i] - rab_y[i] * t_y_yy_xx_yz[i];

        t_yy_yy_xx_yy[i] = t_y_yyy_xx_yy[i] - rab_y[i] * t_y_yy_xx_yy[i];

        t_yy_yy_xx_xz[i] = t_y_yyy_xx_xz[i] - rab_y[i] * t_y_yy_xx_xz[i];

        t_yy_yy_xx_xy[i] = t_y_yyy_xx_xy[i] - rab_y[i] * t_y_yy_xx_xy[i];

        t_yy_yy_xx_xx[i] = t_y_yyy_xx_xx[i] - rab_y[i] * t_y_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xyz_xx_xx, t_y_xyz_xx_xy, t_y_xyz_xx_xz, t_y_xyz_xx_yy,\
                           t_y_xyz_xx_yz, t_y_xyz_xx_zz, t_y_xyz_xy_xx, t_y_xyz_xy_xy,\
                           t_y_xyz_xy_xz, t_y_xyz_xy_yy, t_y_xyz_xy_yz, t_y_xyz_xy_zz,\
                           t_y_xyz_xz_xx, t_y_xyz_xz_xy, t_y_xyz_xz_xz, t_y_xyz_xz_yy,\
                           t_y_xyz_xz_yz, t_y_xyz_xz_zz, t_y_xyz_yy_xx, t_y_xyz_yy_xy,\
                           t_y_xyz_yy_xz, t_y_xyz_yy_yy, t_y_xyz_yy_yz, t_y_xyz_yy_zz,\
                           t_y_xyz_yz_xx, t_y_xyz_yz_xy, t_y_xyz_yz_xz, t_y_xyz_yz_yy,\
                           t_y_xyz_yz_yz, t_y_xyz_yz_zz, t_y_xyz_zz_xx, t_y_xyz_zz_xy,\
                           t_y_xyz_zz_xz, t_y_xyz_zz_yy, t_y_xyz_zz_yz, t_y_xyz_zz_zz,\
                           t_y_xz_xx_xx, t_y_xz_xx_xy, t_y_xz_xx_xz, t_y_xz_xx_yy, t_y_xz_xx_yz,\
                           t_y_xz_xx_zz, t_y_xz_xy_xx, t_y_xz_xy_xy, t_y_xz_xy_xz, t_y_xz_xy_yy,\
                           t_y_xz_xy_yz, t_y_xz_xy_zz, t_y_xz_xz_xx, t_y_xz_xz_xy, t_y_xz_xz_xz,\
                           t_y_xz_xz_yy, t_y_xz_xz_yz, t_y_xz_xz_zz, t_y_xz_yy_xx, t_y_xz_yy_xy,\
                           t_y_xz_yy_xz, t_y_xz_yy_yy, t_y_xz_yy_yz, t_y_xz_yy_zz, t_y_xz_yz_xx,\
                           t_y_xz_yz_xy, t_y_xz_yz_xz, t_y_xz_yz_yy, t_y_xz_yz_yz, t_y_xz_yz_zz,\
                           t_y_xz_zz_xx, t_y_xz_zz_xy, t_y_xz_zz_xz, t_y_xz_zz_yy, t_y_xz_zz_yz,\
                           t_y_xz_zz_zz, t_yy_xz_xx_xx, t_yy_xz_xx_xy, t_yy_xz_xx_xz,\
                           t_yy_xz_xx_yy, t_yy_xz_xx_yz, t_yy_xz_xx_zz, t_yy_xz_xy_xx,\
                           t_yy_xz_xy_xy, t_yy_xz_xy_xz, t_yy_xz_xy_yy, t_yy_xz_xy_yz,\
                           t_yy_xz_xy_zz, t_yy_xz_xz_xx, t_yy_xz_xz_xy, t_yy_xz_xz_xz,\
                           t_yy_xz_xz_yy, t_yy_xz_xz_yz, t_yy_xz_xz_zz, t_yy_xz_yy_xx,\
                           t_yy_xz_yy_xy, t_yy_xz_yy_xz, t_yy_xz_yy_yy, t_yy_xz_yy_yz,\
                           t_yy_xz_yy_zz, t_yy_xz_yz_xx, t_yy_xz_yz_xy, t_yy_xz_yz_xz,\
                           t_yy_xz_yz_yy, t_yy_xz_yz_yz, t_yy_xz_yz_zz, t_yy_xz_zz_xx,\
                           t_yy_xz_zz_xy, t_yy_xz_zz_xz, t_yy_xz_zz_yy, t_yy_xz_zz_yz,\
                           t_yy_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_xz_zz_zz[i] = t_y_xyz_zz_zz[i] - rab_y[i] * t_y_xz_zz_zz[i];

        t_yy_xz_zz_yz[i] = t_y_xyz_zz_yz[i] - rab_y[i] * t_y_xz_zz_yz[i];

        t_yy_xz_zz_yy[i] = t_y_xyz_zz_yy[i] - rab_y[i] * t_y_xz_zz_yy[i];

        t_yy_xz_zz_xz[i] = t_y_xyz_zz_xz[i] - rab_y[i] * t_y_xz_zz_xz[i];

        t_yy_xz_zz_xy[i] = t_y_xyz_zz_xy[i] - rab_y[i] * t_y_xz_zz_xy[i];

        t_yy_xz_zz_xx[i] = t_y_xyz_zz_xx[i] - rab_y[i] * t_y_xz_zz_xx[i];

        t_yy_xz_yz_zz[i] = t_y_xyz_yz_zz[i] - rab_y[i] * t_y_xz_yz_zz[i];

        t_yy_xz_yz_yz[i] = t_y_xyz_yz_yz[i] - rab_y[i] * t_y_xz_yz_yz[i];

        t_yy_xz_yz_yy[i] = t_y_xyz_yz_yy[i] - rab_y[i] * t_y_xz_yz_yy[i];

        t_yy_xz_yz_xz[i] = t_y_xyz_yz_xz[i] - rab_y[i] * t_y_xz_yz_xz[i];

        t_yy_xz_yz_xy[i] = t_y_xyz_yz_xy[i] - rab_y[i] * t_y_xz_yz_xy[i];

        t_yy_xz_yz_xx[i] = t_y_xyz_yz_xx[i] - rab_y[i] * t_y_xz_yz_xx[i];

        t_yy_xz_yy_zz[i] = t_y_xyz_yy_zz[i] - rab_y[i] * t_y_xz_yy_zz[i];

        t_yy_xz_yy_yz[i] = t_y_xyz_yy_yz[i] - rab_y[i] * t_y_xz_yy_yz[i];

        t_yy_xz_yy_yy[i] = t_y_xyz_yy_yy[i] - rab_y[i] * t_y_xz_yy_yy[i];

        t_yy_xz_yy_xz[i] = t_y_xyz_yy_xz[i] - rab_y[i] * t_y_xz_yy_xz[i];

        t_yy_xz_yy_xy[i] = t_y_xyz_yy_xy[i] - rab_y[i] * t_y_xz_yy_xy[i];

        t_yy_xz_yy_xx[i] = t_y_xyz_yy_xx[i] - rab_y[i] * t_y_xz_yy_xx[i];

        t_yy_xz_xz_zz[i] = t_y_xyz_xz_zz[i] - rab_y[i] * t_y_xz_xz_zz[i];

        t_yy_xz_xz_yz[i] = t_y_xyz_xz_yz[i] - rab_y[i] * t_y_xz_xz_yz[i];

        t_yy_xz_xz_yy[i] = t_y_xyz_xz_yy[i] - rab_y[i] * t_y_xz_xz_yy[i];

        t_yy_xz_xz_xz[i] = t_y_xyz_xz_xz[i] - rab_y[i] * t_y_xz_xz_xz[i];

        t_yy_xz_xz_xy[i] = t_y_xyz_xz_xy[i] - rab_y[i] * t_y_xz_xz_xy[i];

        t_yy_xz_xz_xx[i] = t_y_xyz_xz_xx[i] - rab_y[i] * t_y_xz_xz_xx[i];

        t_yy_xz_xy_zz[i] = t_y_xyz_xy_zz[i] - rab_y[i] * t_y_xz_xy_zz[i];

        t_yy_xz_xy_yz[i] = t_y_xyz_xy_yz[i] - rab_y[i] * t_y_xz_xy_yz[i];

        t_yy_xz_xy_yy[i] = t_y_xyz_xy_yy[i] - rab_y[i] * t_y_xz_xy_yy[i];

        t_yy_xz_xy_xz[i] = t_y_xyz_xy_xz[i] - rab_y[i] * t_y_xz_xy_xz[i];

        t_yy_xz_xy_xy[i] = t_y_xyz_xy_xy[i] - rab_y[i] * t_y_xz_xy_xy[i];

        t_yy_xz_xy_xx[i] = t_y_xyz_xy_xx[i] - rab_y[i] * t_y_xz_xy_xx[i];

        t_yy_xz_xx_zz[i] = t_y_xyz_xx_zz[i] - rab_y[i] * t_y_xz_xx_zz[i];

        t_yy_xz_xx_yz[i] = t_y_xyz_xx_yz[i] - rab_y[i] * t_y_xz_xx_yz[i];

        t_yy_xz_xx_yy[i] = t_y_xyz_xx_yy[i] - rab_y[i] * t_y_xz_xx_yy[i];

        t_yy_xz_xx_xz[i] = t_y_xyz_xx_xz[i] - rab_y[i] * t_y_xz_xx_xz[i];

        t_yy_xz_xx_xy[i] = t_y_xyz_xx_xy[i] - rab_y[i] * t_y_xz_xx_xy[i];

        t_yy_xz_xx_xx[i] = t_y_xyz_xx_xx[i] - rab_y[i] * t_y_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xy_xx_xx, t_y_xy_xx_xy, t_y_xy_xx_xz, t_y_xy_xx_yy,\
                           t_y_xy_xx_yz, t_y_xy_xx_zz, t_y_xy_xy_xx, t_y_xy_xy_xy, t_y_xy_xy_xz,\
                           t_y_xy_xy_yy, t_y_xy_xy_yz, t_y_xy_xy_zz, t_y_xy_xz_xx, t_y_xy_xz_xy,\
                           t_y_xy_xz_xz, t_y_xy_xz_yy, t_y_xy_xz_yz, t_y_xy_xz_zz, t_y_xy_yy_xx,\
                           t_y_xy_yy_xy, t_y_xy_yy_xz, t_y_xy_yy_yy, t_y_xy_yy_yz, t_y_xy_yy_zz,\
                           t_y_xy_yz_xx, t_y_xy_yz_xy, t_y_xy_yz_xz, t_y_xy_yz_yy, t_y_xy_yz_yz,\
                           t_y_xy_yz_zz, t_y_xy_zz_xx, t_y_xy_zz_xy, t_y_xy_zz_xz, t_y_xy_zz_yy,\
                           t_y_xy_zz_yz, t_y_xy_zz_zz, t_y_xyy_xx_xx, t_y_xyy_xx_xy,\
                           t_y_xyy_xx_xz, t_y_xyy_xx_yy, t_y_xyy_xx_yz, t_y_xyy_xx_zz,\
                           t_y_xyy_xy_xx, t_y_xyy_xy_xy, t_y_xyy_xy_xz, t_y_xyy_xy_yy,\
                           t_y_xyy_xy_yz, t_y_xyy_xy_zz, t_y_xyy_xz_xx, t_y_xyy_xz_xy,\
                           t_y_xyy_xz_xz, t_y_xyy_xz_yy, t_y_xyy_xz_yz, t_y_xyy_xz_zz,\
                           t_y_xyy_yy_xx, t_y_xyy_yy_xy, t_y_xyy_yy_xz, t_y_xyy_yy_yy,\
                           t_y_xyy_yy_yz, t_y_xyy_yy_zz, t_y_xyy_yz_xx, t_y_xyy_yz_xy,\
                           t_y_xyy_yz_xz, t_y_xyy_yz_yy, t_y_xyy_yz_yz, t_y_xyy_yz_zz,\
                           t_y_xyy_zz_xx, t_y_xyy_zz_xy, t_y_xyy_zz_xz, t_y_xyy_zz_yy,\
                           t_y_xyy_zz_yz, t_y_xyy_zz_zz, t_yy_xy_xx_xx, t_yy_xy_xx_xy,\
                           t_yy_xy_xx_xz, t_yy_xy_xx_yy, t_yy_xy_xx_yz, t_yy_xy_xx_zz,\
                           t_yy_xy_xy_xx, t_yy_xy_xy_xy, t_yy_xy_xy_xz, t_yy_xy_xy_yy,\
                           t_yy_xy_xy_yz, t_yy_xy_xy_zz, t_yy_xy_xz_xx, t_yy_xy_xz_xy,\
                           t_yy_xy_xz_xz, t_yy_xy_xz_yy, t_yy_xy_xz_yz, t_yy_xy_xz_zz,\
                           t_yy_xy_yy_xx, t_yy_xy_yy_xy, t_yy_xy_yy_xz, t_yy_xy_yy_yy,\
                           t_yy_xy_yy_yz, t_yy_xy_yy_zz, t_yy_xy_yz_xx, t_yy_xy_yz_xy,\
                           t_yy_xy_yz_xz, t_yy_xy_yz_yy, t_yy_xy_yz_yz, t_yy_xy_yz_zz,\
                           t_yy_xy_zz_xx, t_yy_xy_zz_xy, t_yy_xy_zz_xz, t_yy_xy_zz_yy,\
                           t_yy_xy_zz_yz, t_yy_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_xy_zz_zz[i] = t_y_xyy_zz_zz[i] - rab_y[i] * t_y_xy_zz_zz[i];

        t_yy_xy_zz_yz[i] = t_y_xyy_zz_yz[i] - rab_y[i] * t_y_xy_zz_yz[i];

        t_yy_xy_zz_yy[i] = t_y_xyy_zz_yy[i] - rab_y[i] * t_y_xy_zz_yy[i];

        t_yy_xy_zz_xz[i] = t_y_xyy_zz_xz[i] - rab_y[i] * t_y_xy_zz_xz[i];

        t_yy_xy_zz_xy[i] = t_y_xyy_zz_xy[i] - rab_y[i] * t_y_xy_zz_xy[i];

        t_yy_xy_zz_xx[i] = t_y_xyy_zz_xx[i] - rab_y[i] * t_y_xy_zz_xx[i];

        t_yy_xy_yz_zz[i] = t_y_xyy_yz_zz[i] - rab_y[i] * t_y_xy_yz_zz[i];

        t_yy_xy_yz_yz[i] = t_y_xyy_yz_yz[i] - rab_y[i] * t_y_xy_yz_yz[i];

        t_yy_xy_yz_yy[i] = t_y_xyy_yz_yy[i] - rab_y[i] * t_y_xy_yz_yy[i];

        t_yy_xy_yz_xz[i] = t_y_xyy_yz_xz[i] - rab_y[i] * t_y_xy_yz_xz[i];

        t_yy_xy_yz_xy[i] = t_y_xyy_yz_xy[i] - rab_y[i] * t_y_xy_yz_xy[i];

        t_yy_xy_yz_xx[i] = t_y_xyy_yz_xx[i] - rab_y[i] * t_y_xy_yz_xx[i];

        t_yy_xy_yy_zz[i] = t_y_xyy_yy_zz[i] - rab_y[i] * t_y_xy_yy_zz[i];

        t_yy_xy_yy_yz[i] = t_y_xyy_yy_yz[i] - rab_y[i] * t_y_xy_yy_yz[i];

        t_yy_xy_yy_yy[i] = t_y_xyy_yy_yy[i] - rab_y[i] * t_y_xy_yy_yy[i];

        t_yy_xy_yy_xz[i] = t_y_xyy_yy_xz[i] - rab_y[i] * t_y_xy_yy_xz[i];

        t_yy_xy_yy_xy[i] = t_y_xyy_yy_xy[i] - rab_y[i] * t_y_xy_yy_xy[i];

        t_yy_xy_yy_xx[i] = t_y_xyy_yy_xx[i] - rab_y[i] * t_y_xy_yy_xx[i];

        t_yy_xy_xz_zz[i] = t_y_xyy_xz_zz[i] - rab_y[i] * t_y_xy_xz_zz[i];

        t_yy_xy_xz_yz[i] = t_y_xyy_xz_yz[i] - rab_y[i] * t_y_xy_xz_yz[i];

        t_yy_xy_xz_yy[i] = t_y_xyy_xz_yy[i] - rab_y[i] * t_y_xy_xz_yy[i];

        t_yy_xy_xz_xz[i] = t_y_xyy_xz_xz[i] - rab_y[i] * t_y_xy_xz_xz[i];

        t_yy_xy_xz_xy[i] = t_y_xyy_xz_xy[i] - rab_y[i] * t_y_xy_xz_xy[i];

        t_yy_xy_xz_xx[i] = t_y_xyy_xz_xx[i] - rab_y[i] * t_y_xy_xz_xx[i];

        t_yy_xy_xy_zz[i] = t_y_xyy_xy_zz[i] - rab_y[i] * t_y_xy_xy_zz[i];

        t_yy_xy_xy_yz[i] = t_y_xyy_xy_yz[i] - rab_y[i] * t_y_xy_xy_yz[i];

        t_yy_xy_xy_yy[i] = t_y_xyy_xy_yy[i] - rab_y[i] * t_y_xy_xy_yy[i];

        t_yy_xy_xy_xz[i] = t_y_xyy_xy_xz[i] - rab_y[i] * t_y_xy_xy_xz[i];

        t_yy_xy_xy_xy[i] = t_y_xyy_xy_xy[i] - rab_y[i] * t_y_xy_xy_xy[i];

        t_yy_xy_xy_xx[i] = t_y_xyy_xy_xx[i] - rab_y[i] * t_y_xy_xy_xx[i];

        t_yy_xy_xx_zz[i] = t_y_xyy_xx_zz[i] - rab_y[i] * t_y_xy_xx_zz[i];

        t_yy_xy_xx_yz[i] = t_y_xyy_xx_yz[i] - rab_y[i] * t_y_xy_xx_yz[i];

        t_yy_xy_xx_yy[i] = t_y_xyy_xx_yy[i] - rab_y[i] * t_y_xy_xx_yy[i];

        t_yy_xy_xx_xz[i] = t_y_xyy_xx_xz[i] - rab_y[i] * t_y_xy_xx_xz[i];

        t_yy_xy_xx_xy[i] = t_y_xyy_xx_xy[i] - rab_y[i] * t_y_xy_xx_xy[i];

        t_yy_xy_xx_xx[i] = t_y_xyy_xx_xx[i] - rab_y[i] * t_y_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xx_xx_xx, t_y_xx_xx_xy, t_y_xx_xx_xz, t_y_xx_xx_yy,\
                           t_y_xx_xx_yz, t_y_xx_xx_zz, t_y_xx_xy_xx, t_y_xx_xy_xy, t_y_xx_xy_xz,\
                           t_y_xx_xy_yy, t_y_xx_xy_yz, t_y_xx_xy_zz, t_y_xx_xz_xx, t_y_xx_xz_xy,\
                           t_y_xx_xz_xz, t_y_xx_xz_yy, t_y_xx_xz_yz, t_y_xx_xz_zz, t_y_xx_yy_xx,\
                           t_y_xx_yy_xy, t_y_xx_yy_xz, t_y_xx_yy_yy, t_y_xx_yy_yz, t_y_xx_yy_zz,\
                           t_y_xx_yz_xx, t_y_xx_yz_xy, t_y_xx_yz_xz, t_y_xx_yz_yy, t_y_xx_yz_yz,\
                           t_y_xx_yz_zz, t_y_xx_zz_xx, t_y_xx_zz_xy, t_y_xx_zz_xz, t_y_xx_zz_yy,\
                           t_y_xx_zz_yz, t_y_xx_zz_zz, t_y_xxy_xx_xx, t_y_xxy_xx_xy,\
                           t_y_xxy_xx_xz, t_y_xxy_xx_yy, t_y_xxy_xx_yz, t_y_xxy_xx_zz,\
                           t_y_xxy_xy_xx, t_y_xxy_xy_xy, t_y_xxy_xy_xz, t_y_xxy_xy_yy,\
                           t_y_xxy_xy_yz, t_y_xxy_xy_zz, t_y_xxy_xz_xx, t_y_xxy_xz_xy,\
                           t_y_xxy_xz_xz, t_y_xxy_xz_yy, t_y_xxy_xz_yz, t_y_xxy_xz_zz,\
                           t_y_xxy_yy_xx, t_y_xxy_yy_xy, t_y_xxy_yy_xz, t_y_xxy_yy_yy,\
                           t_y_xxy_yy_yz, t_y_xxy_yy_zz, t_y_xxy_yz_xx, t_y_xxy_yz_xy,\
                           t_y_xxy_yz_xz, t_y_xxy_yz_yy, t_y_xxy_yz_yz, t_y_xxy_yz_zz,\
                           t_y_xxy_zz_xx, t_y_xxy_zz_xy, t_y_xxy_zz_xz, t_y_xxy_zz_yy,\
                           t_y_xxy_zz_yz, t_y_xxy_zz_zz, t_yy_xx_xx_xx, t_yy_xx_xx_xy,\
                           t_yy_xx_xx_xz, t_yy_xx_xx_yy, t_yy_xx_xx_yz, t_yy_xx_xx_zz,\
                           t_yy_xx_xy_xx, t_yy_xx_xy_xy, t_yy_xx_xy_xz, t_yy_xx_xy_yy,\
                           t_yy_xx_xy_yz, t_yy_xx_xy_zz, t_yy_xx_xz_xx, t_yy_xx_xz_xy,\
                           t_yy_xx_xz_xz, t_yy_xx_xz_yy, t_yy_xx_xz_yz, t_yy_xx_xz_zz,\
                           t_yy_xx_yy_xx, t_yy_xx_yy_xy, t_yy_xx_yy_xz, t_yy_xx_yy_yy,\
                           t_yy_xx_yy_yz, t_yy_xx_yy_zz, t_yy_xx_yz_xx, t_yy_xx_yz_xy,\
                           t_yy_xx_yz_xz, t_yy_xx_yz_yy, t_yy_xx_yz_yz, t_yy_xx_yz_zz,\
                           t_yy_xx_zz_xx, t_yy_xx_zz_xy, t_yy_xx_zz_xz, t_yy_xx_zz_yy,\
                           t_yy_xx_zz_yz, t_yy_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_xx_zz_zz[i] = t_y_xxy_zz_zz[i] - rab_y[i] * t_y_xx_zz_zz[i];

        t_yy_xx_zz_yz[i] = t_y_xxy_zz_yz[i] - rab_y[i] * t_y_xx_zz_yz[i];

        t_yy_xx_zz_yy[i] = t_y_xxy_zz_yy[i] - rab_y[i] * t_y_xx_zz_yy[i];

        t_yy_xx_zz_xz[i] = t_y_xxy_zz_xz[i] - rab_y[i] * t_y_xx_zz_xz[i];

        t_yy_xx_zz_xy[i] = t_y_xxy_zz_xy[i] - rab_y[i] * t_y_xx_zz_xy[i];

        t_yy_xx_zz_xx[i] = t_y_xxy_zz_xx[i] - rab_y[i] * t_y_xx_zz_xx[i];

        t_yy_xx_yz_zz[i] = t_y_xxy_yz_zz[i] - rab_y[i] * t_y_xx_yz_zz[i];

        t_yy_xx_yz_yz[i] = t_y_xxy_yz_yz[i] - rab_y[i] * t_y_xx_yz_yz[i];

        t_yy_xx_yz_yy[i] = t_y_xxy_yz_yy[i] - rab_y[i] * t_y_xx_yz_yy[i];

        t_yy_xx_yz_xz[i] = t_y_xxy_yz_xz[i] - rab_y[i] * t_y_xx_yz_xz[i];

        t_yy_xx_yz_xy[i] = t_y_xxy_yz_xy[i] - rab_y[i] * t_y_xx_yz_xy[i];

        t_yy_xx_yz_xx[i] = t_y_xxy_yz_xx[i] - rab_y[i] * t_y_xx_yz_xx[i];

        t_yy_xx_yy_zz[i] = t_y_xxy_yy_zz[i] - rab_y[i] * t_y_xx_yy_zz[i];

        t_yy_xx_yy_yz[i] = t_y_xxy_yy_yz[i] - rab_y[i] * t_y_xx_yy_yz[i];

        t_yy_xx_yy_yy[i] = t_y_xxy_yy_yy[i] - rab_y[i] * t_y_xx_yy_yy[i];

        t_yy_xx_yy_xz[i] = t_y_xxy_yy_xz[i] - rab_y[i] * t_y_xx_yy_xz[i];

        t_yy_xx_yy_xy[i] = t_y_xxy_yy_xy[i] - rab_y[i] * t_y_xx_yy_xy[i];

        t_yy_xx_yy_xx[i] = t_y_xxy_yy_xx[i] - rab_y[i] * t_y_xx_yy_xx[i];

        t_yy_xx_xz_zz[i] = t_y_xxy_xz_zz[i] - rab_y[i] * t_y_xx_xz_zz[i];

        t_yy_xx_xz_yz[i] = t_y_xxy_xz_yz[i] - rab_y[i] * t_y_xx_xz_yz[i];

        t_yy_xx_xz_yy[i] = t_y_xxy_xz_yy[i] - rab_y[i] * t_y_xx_xz_yy[i];

        t_yy_xx_xz_xz[i] = t_y_xxy_xz_xz[i] - rab_y[i] * t_y_xx_xz_xz[i];

        t_yy_xx_xz_xy[i] = t_y_xxy_xz_xy[i] - rab_y[i] * t_y_xx_xz_xy[i];

        t_yy_xx_xz_xx[i] = t_y_xxy_xz_xx[i] - rab_y[i] * t_y_xx_xz_xx[i];

        t_yy_xx_xy_zz[i] = t_y_xxy_xy_zz[i] - rab_y[i] * t_y_xx_xy_zz[i];

        t_yy_xx_xy_yz[i] = t_y_xxy_xy_yz[i] - rab_y[i] * t_y_xx_xy_yz[i];

        t_yy_xx_xy_yy[i] = t_y_xxy_xy_yy[i] - rab_y[i] * t_y_xx_xy_yy[i];

        t_yy_xx_xy_xz[i] = t_y_xxy_xy_xz[i] - rab_y[i] * t_y_xx_xy_xz[i];

        t_yy_xx_xy_xy[i] = t_y_xxy_xy_xy[i] - rab_y[i] * t_y_xx_xy_xy[i];

        t_yy_xx_xy_xx[i] = t_y_xxy_xy_xx[i] - rab_y[i] * t_y_xx_xy_xx[i];

        t_yy_xx_xx_zz[i] = t_y_xxy_xx_zz[i] - rab_y[i] * t_y_xx_xx_zz[i];

        t_yy_xx_xx_yz[i] = t_y_xxy_xx_yz[i] - rab_y[i] * t_y_xx_xx_yz[i];

        t_yy_xx_xx_yy[i] = t_y_xxy_xx_yy[i] - rab_y[i] * t_y_xx_xx_yy[i];

        t_yy_xx_xx_xz[i] = t_y_xxy_xx_xz[i] - rab_y[i] * t_y_xx_xx_xz[i];

        t_yy_xx_xx_xy[i] = t_y_xxy_xx_xy[i] - rab_y[i] * t_y_xx_xx_xy[i];

        t_yy_xx_xx_xx[i] = t_y_xxy_xx_xx[i] - rab_y[i] * t_y_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_zz_xx_xx, t_x_zz_xx_xy, t_x_zz_xx_xz, t_x_zz_xx_yy,\
                           t_x_zz_xx_yz, t_x_zz_xx_zz, t_x_zz_xy_xx, t_x_zz_xy_xy, t_x_zz_xy_xz,\
                           t_x_zz_xy_yy, t_x_zz_xy_yz, t_x_zz_xy_zz, t_x_zz_xz_xx, t_x_zz_xz_xy,\
                           t_x_zz_xz_xz, t_x_zz_xz_yy, t_x_zz_xz_yz, t_x_zz_xz_zz, t_x_zz_yy_xx,\
                           t_x_zz_yy_xy, t_x_zz_yy_xz, t_x_zz_yy_yy, t_x_zz_yy_yz, t_x_zz_yy_zz,\
                           t_x_zz_yz_xx, t_x_zz_yz_xy, t_x_zz_yz_xz, t_x_zz_yz_yy, t_x_zz_yz_yz,\
                           t_x_zz_yz_zz, t_x_zz_zz_xx, t_x_zz_zz_xy, t_x_zz_zz_xz, t_x_zz_zz_yy,\
                           t_x_zz_zz_yz, t_x_zz_zz_zz, t_x_zzz_xx_xx, t_x_zzz_xx_xy,\
                           t_x_zzz_xx_xz, t_x_zzz_xx_yy, t_x_zzz_xx_yz, t_x_zzz_xx_zz,\
                           t_x_zzz_xy_xx, t_x_zzz_xy_xy, t_x_zzz_xy_xz, t_x_zzz_xy_yy,\
                           t_x_zzz_xy_yz, t_x_zzz_xy_zz, t_x_zzz_xz_xx, t_x_zzz_xz_xy,\
                           t_x_zzz_xz_xz, t_x_zzz_xz_yy, t_x_zzz_xz_yz, t_x_zzz_xz_zz,\
                           t_x_zzz_yy_xx, t_x_zzz_yy_xy, t_x_zzz_yy_xz, t_x_zzz_yy_yy,\
                           t_x_zzz_yy_yz, t_x_zzz_yy_zz, t_x_zzz_yz_xx, t_x_zzz_yz_xy,\
                           t_x_zzz_yz_xz, t_x_zzz_yz_yy, t_x_zzz_yz_yz, t_x_zzz_yz_zz,\
                           t_x_zzz_zz_xx, t_x_zzz_zz_xy, t_x_zzz_zz_xz, t_x_zzz_zz_yy,\
                           t_x_zzz_zz_yz, t_x_zzz_zz_zz, t_xz_zz_xx_xx, t_xz_zz_xx_xy,\
                           t_xz_zz_xx_xz, t_xz_zz_xx_yy, t_xz_zz_xx_yz, t_xz_zz_xx_zz,\
                           t_xz_zz_xy_xx, t_xz_zz_xy_xy, t_xz_zz_xy_xz, t_xz_zz_xy_yy,\
                           t_xz_zz_xy_yz, t_xz_zz_xy_zz, t_xz_zz_xz_xx, t_xz_zz_xz_xy,\
                           t_xz_zz_xz_xz, t_xz_zz_xz_yy, t_xz_zz_xz_yz, t_xz_zz_xz_zz,\
                           t_xz_zz_yy_xx, t_xz_zz_yy_xy, t_xz_zz_yy_xz, t_xz_zz_yy_yy,\
                           t_xz_zz_yy_yz, t_xz_zz_yy_zz, t_xz_zz_yz_xx, t_xz_zz_yz_xy,\
                           t_xz_zz_yz_xz, t_xz_zz_yz_yy, t_xz_zz_yz_yz, t_xz_zz_yz_zz,\
                           t_xz_zz_zz_xx, t_xz_zz_zz_xy, t_xz_zz_zz_xz, t_xz_zz_zz_yy,\
                           t_xz_zz_zz_yz, t_xz_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_zz_zz_zz[i] = t_x_zzz_zz_zz[i] - rab_z[i] * t_x_zz_zz_zz[i];

        t_xz_zz_zz_yz[i] = t_x_zzz_zz_yz[i] - rab_z[i] * t_x_zz_zz_yz[i];

        t_xz_zz_zz_yy[i] = t_x_zzz_zz_yy[i] - rab_z[i] * t_x_zz_zz_yy[i];

        t_xz_zz_zz_xz[i] = t_x_zzz_zz_xz[i] - rab_z[i] * t_x_zz_zz_xz[i];

        t_xz_zz_zz_xy[i] = t_x_zzz_zz_xy[i] - rab_z[i] * t_x_zz_zz_xy[i];

        t_xz_zz_zz_xx[i] = t_x_zzz_zz_xx[i] - rab_z[i] * t_x_zz_zz_xx[i];

        t_xz_zz_yz_zz[i] = t_x_zzz_yz_zz[i] - rab_z[i] * t_x_zz_yz_zz[i];

        t_xz_zz_yz_yz[i] = t_x_zzz_yz_yz[i] - rab_z[i] * t_x_zz_yz_yz[i];

        t_xz_zz_yz_yy[i] = t_x_zzz_yz_yy[i] - rab_z[i] * t_x_zz_yz_yy[i];

        t_xz_zz_yz_xz[i] = t_x_zzz_yz_xz[i] - rab_z[i] * t_x_zz_yz_xz[i];

        t_xz_zz_yz_xy[i] = t_x_zzz_yz_xy[i] - rab_z[i] * t_x_zz_yz_xy[i];

        t_xz_zz_yz_xx[i] = t_x_zzz_yz_xx[i] - rab_z[i] * t_x_zz_yz_xx[i];

        t_xz_zz_yy_zz[i] = t_x_zzz_yy_zz[i] - rab_z[i] * t_x_zz_yy_zz[i];

        t_xz_zz_yy_yz[i] = t_x_zzz_yy_yz[i] - rab_z[i] * t_x_zz_yy_yz[i];

        t_xz_zz_yy_yy[i] = t_x_zzz_yy_yy[i] - rab_z[i] * t_x_zz_yy_yy[i];

        t_xz_zz_yy_xz[i] = t_x_zzz_yy_xz[i] - rab_z[i] * t_x_zz_yy_xz[i];

        t_xz_zz_yy_xy[i] = t_x_zzz_yy_xy[i] - rab_z[i] * t_x_zz_yy_xy[i];

        t_xz_zz_yy_xx[i] = t_x_zzz_yy_xx[i] - rab_z[i] * t_x_zz_yy_xx[i];

        t_xz_zz_xz_zz[i] = t_x_zzz_xz_zz[i] - rab_z[i] * t_x_zz_xz_zz[i];

        t_xz_zz_xz_yz[i] = t_x_zzz_xz_yz[i] - rab_z[i] * t_x_zz_xz_yz[i];

        t_xz_zz_xz_yy[i] = t_x_zzz_xz_yy[i] - rab_z[i] * t_x_zz_xz_yy[i];

        t_xz_zz_xz_xz[i] = t_x_zzz_xz_xz[i] - rab_z[i] * t_x_zz_xz_xz[i];

        t_xz_zz_xz_xy[i] = t_x_zzz_xz_xy[i] - rab_z[i] * t_x_zz_xz_xy[i];

        t_xz_zz_xz_xx[i] = t_x_zzz_xz_xx[i] - rab_z[i] * t_x_zz_xz_xx[i];

        t_xz_zz_xy_zz[i] = t_x_zzz_xy_zz[i] - rab_z[i] * t_x_zz_xy_zz[i];

        t_xz_zz_xy_yz[i] = t_x_zzz_xy_yz[i] - rab_z[i] * t_x_zz_xy_yz[i];

        t_xz_zz_xy_yy[i] = t_x_zzz_xy_yy[i] - rab_z[i] * t_x_zz_xy_yy[i];

        t_xz_zz_xy_xz[i] = t_x_zzz_xy_xz[i] - rab_z[i] * t_x_zz_xy_xz[i];

        t_xz_zz_xy_xy[i] = t_x_zzz_xy_xy[i] - rab_z[i] * t_x_zz_xy_xy[i];

        t_xz_zz_xy_xx[i] = t_x_zzz_xy_xx[i] - rab_z[i] * t_x_zz_xy_xx[i];

        t_xz_zz_xx_zz[i] = t_x_zzz_xx_zz[i] - rab_z[i] * t_x_zz_xx_zz[i];

        t_xz_zz_xx_yz[i] = t_x_zzz_xx_yz[i] - rab_z[i] * t_x_zz_xx_yz[i];

        t_xz_zz_xx_yy[i] = t_x_zzz_xx_yy[i] - rab_z[i] * t_x_zz_xx_yy[i];

        t_xz_zz_xx_xz[i] = t_x_zzz_xx_xz[i] - rab_z[i] * t_x_zz_xx_xz[i];

        t_xz_zz_xx_xy[i] = t_x_zzz_xx_xy[i] - rab_z[i] * t_x_zz_xx_xy[i];

        t_xz_zz_xx_xx[i] = t_x_zzz_xx_xx[i] - rab_z[i] * t_x_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_yz_xx_xx, t_x_yz_xx_xy, t_x_yz_xx_xz, t_x_yz_xx_yy,\
                           t_x_yz_xx_yz, t_x_yz_xx_zz, t_x_yz_xy_xx, t_x_yz_xy_xy, t_x_yz_xy_xz,\
                           t_x_yz_xy_yy, t_x_yz_xy_yz, t_x_yz_xy_zz, t_x_yz_xz_xx, t_x_yz_xz_xy,\
                           t_x_yz_xz_xz, t_x_yz_xz_yy, t_x_yz_xz_yz, t_x_yz_xz_zz, t_x_yz_yy_xx,\
                           t_x_yz_yy_xy, t_x_yz_yy_xz, t_x_yz_yy_yy, t_x_yz_yy_yz, t_x_yz_yy_zz,\
                           t_x_yz_yz_xx, t_x_yz_yz_xy, t_x_yz_yz_xz, t_x_yz_yz_yy, t_x_yz_yz_yz,\
                           t_x_yz_yz_zz, t_x_yz_zz_xx, t_x_yz_zz_xy, t_x_yz_zz_xz, t_x_yz_zz_yy,\
                           t_x_yz_zz_yz, t_x_yz_zz_zz, t_x_yzz_xx_xx, t_x_yzz_xx_xy,\
                           t_x_yzz_xx_xz, t_x_yzz_xx_yy, t_x_yzz_xx_yz, t_x_yzz_xx_zz,\
                           t_x_yzz_xy_xx, t_x_yzz_xy_xy, t_x_yzz_xy_xz, t_x_yzz_xy_yy,\
                           t_x_yzz_xy_yz, t_x_yzz_xy_zz, t_x_yzz_xz_xx, t_x_yzz_xz_xy,\
                           t_x_yzz_xz_xz, t_x_yzz_xz_yy, t_x_yzz_xz_yz, t_x_yzz_xz_zz,\
                           t_x_yzz_yy_xx, t_x_yzz_yy_xy, t_x_yzz_yy_xz, t_x_yzz_yy_yy,\
                           t_x_yzz_yy_yz, t_x_yzz_yy_zz, t_x_yzz_yz_xx, t_x_yzz_yz_xy,\
                           t_x_yzz_yz_xz, t_x_yzz_yz_yy, t_x_yzz_yz_yz, t_x_yzz_yz_zz,\
                           t_x_yzz_zz_xx, t_x_yzz_zz_xy, t_x_yzz_zz_xz, t_x_yzz_zz_yy,\
                           t_x_yzz_zz_yz, t_x_yzz_zz_zz, t_xz_yz_xx_xx, t_xz_yz_xx_xy,\
                           t_xz_yz_xx_xz, t_xz_yz_xx_yy, t_xz_yz_xx_yz, t_xz_yz_xx_zz,\
                           t_xz_yz_xy_xx, t_xz_yz_xy_xy, t_xz_yz_xy_xz, t_xz_yz_xy_yy,\
                           t_xz_yz_xy_yz, t_xz_yz_xy_zz, t_xz_yz_xz_xx, t_xz_yz_xz_xy,\
                           t_xz_yz_xz_xz, t_xz_yz_xz_yy, t_xz_yz_xz_yz, t_xz_yz_xz_zz,\
                           t_xz_yz_yy_xx, t_xz_yz_yy_xy, t_xz_yz_yy_xz, t_xz_yz_yy_yy,\
                           t_xz_yz_yy_yz, t_xz_yz_yy_zz, t_xz_yz_yz_xx, t_xz_yz_yz_xy,\
                           t_xz_yz_yz_xz, t_xz_yz_yz_yy, t_xz_yz_yz_yz, t_xz_yz_yz_zz,\
                           t_xz_yz_zz_xx, t_xz_yz_zz_xy, t_xz_yz_zz_xz, t_xz_yz_zz_yy,\
                           t_xz_yz_zz_yz, t_xz_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_yz_zz_zz[i] = t_x_yzz_zz_zz[i] - rab_z[i] * t_x_yz_zz_zz[i];

        t_xz_yz_zz_yz[i] = t_x_yzz_zz_yz[i] - rab_z[i] * t_x_yz_zz_yz[i];

        t_xz_yz_zz_yy[i] = t_x_yzz_zz_yy[i] - rab_z[i] * t_x_yz_zz_yy[i];

        t_xz_yz_zz_xz[i] = t_x_yzz_zz_xz[i] - rab_z[i] * t_x_yz_zz_xz[i];

        t_xz_yz_zz_xy[i] = t_x_yzz_zz_xy[i] - rab_z[i] * t_x_yz_zz_xy[i];

        t_xz_yz_zz_xx[i] = t_x_yzz_zz_xx[i] - rab_z[i] * t_x_yz_zz_xx[i];

        t_xz_yz_yz_zz[i] = t_x_yzz_yz_zz[i] - rab_z[i] * t_x_yz_yz_zz[i];

        t_xz_yz_yz_yz[i] = t_x_yzz_yz_yz[i] - rab_z[i] * t_x_yz_yz_yz[i];

        t_xz_yz_yz_yy[i] = t_x_yzz_yz_yy[i] - rab_z[i] * t_x_yz_yz_yy[i];

        t_xz_yz_yz_xz[i] = t_x_yzz_yz_xz[i] - rab_z[i] * t_x_yz_yz_xz[i];

        t_xz_yz_yz_xy[i] = t_x_yzz_yz_xy[i] - rab_z[i] * t_x_yz_yz_xy[i];

        t_xz_yz_yz_xx[i] = t_x_yzz_yz_xx[i] - rab_z[i] * t_x_yz_yz_xx[i];

        t_xz_yz_yy_zz[i] = t_x_yzz_yy_zz[i] - rab_z[i] * t_x_yz_yy_zz[i];

        t_xz_yz_yy_yz[i] = t_x_yzz_yy_yz[i] - rab_z[i] * t_x_yz_yy_yz[i];

        t_xz_yz_yy_yy[i] = t_x_yzz_yy_yy[i] - rab_z[i] * t_x_yz_yy_yy[i];

        t_xz_yz_yy_xz[i] = t_x_yzz_yy_xz[i] - rab_z[i] * t_x_yz_yy_xz[i];

        t_xz_yz_yy_xy[i] = t_x_yzz_yy_xy[i] - rab_z[i] * t_x_yz_yy_xy[i];

        t_xz_yz_yy_xx[i] = t_x_yzz_yy_xx[i] - rab_z[i] * t_x_yz_yy_xx[i];

        t_xz_yz_xz_zz[i] = t_x_yzz_xz_zz[i] - rab_z[i] * t_x_yz_xz_zz[i];

        t_xz_yz_xz_yz[i] = t_x_yzz_xz_yz[i] - rab_z[i] * t_x_yz_xz_yz[i];

        t_xz_yz_xz_yy[i] = t_x_yzz_xz_yy[i] - rab_z[i] * t_x_yz_xz_yy[i];

        t_xz_yz_xz_xz[i] = t_x_yzz_xz_xz[i] - rab_z[i] * t_x_yz_xz_xz[i];

        t_xz_yz_xz_xy[i] = t_x_yzz_xz_xy[i] - rab_z[i] * t_x_yz_xz_xy[i];

        t_xz_yz_xz_xx[i] = t_x_yzz_xz_xx[i] - rab_z[i] * t_x_yz_xz_xx[i];

        t_xz_yz_xy_zz[i] = t_x_yzz_xy_zz[i] - rab_z[i] * t_x_yz_xy_zz[i];

        t_xz_yz_xy_yz[i] = t_x_yzz_xy_yz[i] - rab_z[i] * t_x_yz_xy_yz[i];

        t_xz_yz_xy_yy[i] = t_x_yzz_xy_yy[i] - rab_z[i] * t_x_yz_xy_yy[i];

        t_xz_yz_xy_xz[i] = t_x_yzz_xy_xz[i] - rab_z[i] * t_x_yz_xy_xz[i];

        t_xz_yz_xy_xy[i] = t_x_yzz_xy_xy[i] - rab_z[i] * t_x_yz_xy_xy[i];

        t_xz_yz_xy_xx[i] = t_x_yzz_xy_xx[i] - rab_z[i] * t_x_yz_xy_xx[i];

        t_xz_yz_xx_zz[i] = t_x_yzz_xx_zz[i] - rab_z[i] * t_x_yz_xx_zz[i];

        t_xz_yz_xx_yz[i] = t_x_yzz_xx_yz[i] - rab_z[i] * t_x_yz_xx_yz[i];

        t_xz_yz_xx_yy[i] = t_x_yzz_xx_yy[i] - rab_z[i] * t_x_yz_xx_yy[i];

        t_xz_yz_xx_xz[i] = t_x_yzz_xx_xz[i] - rab_z[i] * t_x_yz_xx_xz[i];

        t_xz_yz_xx_xy[i] = t_x_yzz_xx_xy[i] - rab_z[i] * t_x_yz_xx_xy[i];

        t_xz_yz_xx_xx[i] = t_x_yzz_xx_xx[i] - rab_z[i] * t_x_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_yy_xx_xx, t_x_yy_xx_xy, t_x_yy_xx_xz, t_x_yy_xx_yy,\
                           t_x_yy_xx_yz, t_x_yy_xx_zz, t_x_yy_xy_xx, t_x_yy_xy_xy, t_x_yy_xy_xz,\
                           t_x_yy_xy_yy, t_x_yy_xy_yz, t_x_yy_xy_zz, t_x_yy_xz_xx, t_x_yy_xz_xy,\
                           t_x_yy_xz_xz, t_x_yy_xz_yy, t_x_yy_xz_yz, t_x_yy_xz_zz, t_x_yy_yy_xx,\
                           t_x_yy_yy_xy, t_x_yy_yy_xz, t_x_yy_yy_yy, t_x_yy_yy_yz, t_x_yy_yy_zz,\
                           t_x_yy_yz_xx, t_x_yy_yz_xy, t_x_yy_yz_xz, t_x_yy_yz_yy, t_x_yy_yz_yz,\
                           t_x_yy_yz_zz, t_x_yy_zz_xx, t_x_yy_zz_xy, t_x_yy_zz_xz, t_x_yy_zz_yy,\
                           t_x_yy_zz_yz, t_x_yy_zz_zz, t_x_yyz_xx_xx, t_x_yyz_xx_xy,\
                           t_x_yyz_xx_xz, t_x_yyz_xx_yy, t_x_yyz_xx_yz, t_x_yyz_xx_zz,\
                           t_x_yyz_xy_xx, t_x_yyz_xy_xy, t_x_yyz_xy_xz, t_x_yyz_xy_yy,\
                           t_x_yyz_xy_yz, t_x_yyz_xy_zz, t_x_yyz_xz_xx, t_x_yyz_xz_xy,\
                           t_x_yyz_xz_xz, t_x_yyz_xz_yy, t_x_yyz_xz_yz, t_x_yyz_xz_zz,\
                           t_x_yyz_yy_xx, t_x_yyz_yy_xy, t_x_yyz_yy_xz, t_x_yyz_yy_yy,\
                           t_x_yyz_yy_yz, t_x_yyz_yy_zz, t_x_yyz_yz_xx, t_x_yyz_yz_xy,\
                           t_x_yyz_yz_xz, t_x_yyz_yz_yy, t_x_yyz_yz_yz, t_x_yyz_yz_zz,\
                           t_x_yyz_zz_xx, t_x_yyz_zz_xy, t_x_yyz_zz_xz, t_x_yyz_zz_yy,\
                           t_x_yyz_zz_yz, t_x_yyz_zz_zz, t_xz_yy_xx_xx, t_xz_yy_xx_xy,\
                           t_xz_yy_xx_xz, t_xz_yy_xx_yy, t_xz_yy_xx_yz, t_xz_yy_xx_zz,\
                           t_xz_yy_xy_xx, t_xz_yy_xy_xy, t_xz_yy_xy_xz, t_xz_yy_xy_yy,\
                           t_xz_yy_xy_yz, t_xz_yy_xy_zz, t_xz_yy_xz_xx, t_xz_yy_xz_xy,\
                           t_xz_yy_xz_xz, t_xz_yy_xz_yy, t_xz_yy_xz_yz, t_xz_yy_xz_zz,\
                           t_xz_yy_yy_xx, t_xz_yy_yy_xy, t_xz_yy_yy_xz, t_xz_yy_yy_yy,\
                           t_xz_yy_yy_yz, t_xz_yy_yy_zz, t_xz_yy_yz_xx, t_xz_yy_yz_xy,\
                           t_xz_yy_yz_xz, t_xz_yy_yz_yy, t_xz_yy_yz_yz, t_xz_yy_yz_zz,\
                           t_xz_yy_zz_xx, t_xz_yy_zz_xy, t_xz_yy_zz_xz, t_xz_yy_zz_yy,\
                           t_xz_yy_zz_yz, t_xz_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_yy_zz_zz[i] = t_x_yyz_zz_zz[i] - rab_z[i] * t_x_yy_zz_zz[i];

        t_xz_yy_zz_yz[i] = t_x_yyz_zz_yz[i] - rab_z[i] * t_x_yy_zz_yz[i];

        t_xz_yy_zz_yy[i] = t_x_yyz_zz_yy[i] - rab_z[i] * t_x_yy_zz_yy[i];

        t_xz_yy_zz_xz[i] = t_x_yyz_zz_xz[i] - rab_z[i] * t_x_yy_zz_xz[i];

        t_xz_yy_zz_xy[i] = t_x_yyz_zz_xy[i] - rab_z[i] * t_x_yy_zz_xy[i];

        t_xz_yy_zz_xx[i] = t_x_yyz_zz_xx[i] - rab_z[i] * t_x_yy_zz_xx[i];

        t_xz_yy_yz_zz[i] = t_x_yyz_yz_zz[i] - rab_z[i] * t_x_yy_yz_zz[i];

        t_xz_yy_yz_yz[i] = t_x_yyz_yz_yz[i] - rab_z[i] * t_x_yy_yz_yz[i];

        t_xz_yy_yz_yy[i] = t_x_yyz_yz_yy[i] - rab_z[i] * t_x_yy_yz_yy[i];

        t_xz_yy_yz_xz[i] = t_x_yyz_yz_xz[i] - rab_z[i] * t_x_yy_yz_xz[i];

        t_xz_yy_yz_xy[i] = t_x_yyz_yz_xy[i] - rab_z[i] * t_x_yy_yz_xy[i];

        t_xz_yy_yz_xx[i] = t_x_yyz_yz_xx[i] - rab_z[i] * t_x_yy_yz_xx[i];

        t_xz_yy_yy_zz[i] = t_x_yyz_yy_zz[i] - rab_z[i] * t_x_yy_yy_zz[i];

        t_xz_yy_yy_yz[i] = t_x_yyz_yy_yz[i] - rab_z[i] * t_x_yy_yy_yz[i];

        t_xz_yy_yy_yy[i] = t_x_yyz_yy_yy[i] - rab_z[i] * t_x_yy_yy_yy[i];

        t_xz_yy_yy_xz[i] = t_x_yyz_yy_xz[i] - rab_z[i] * t_x_yy_yy_xz[i];

        t_xz_yy_yy_xy[i] = t_x_yyz_yy_xy[i] - rab_z[i] * t_x_yy_yy_xy[i];

        t_xz_yy_yy_xx[i] = t_x_yyz_yy_xx[i] - rab_z[i] * t_x_yy_yy_xx[i];

        t_xz_yy_xz_zz[i] = t_x_yyz_xz_zz[i] - rab_z[i] * t_x_yy_xz_zz[i];

        t_xz_yy_xz_yz[i] = t_x_yyz_xz_yz[i] - rab_z[i] * t_x_yy_xz_yz[i];

        t_xz_yy_xz_yy[i] = t_x_yyz_xz_yy[i] - rab_z[i] * t_x_yy_xz_yy[i];

        t_xz_yy_xz_xz[i] = t_x_yyz_xz_xz[i] - rab_z[i] * t_x_yy_xz_xz[i];

        t_xz_yy_xz_xy[i] = t_x_yyz_xz_xy[i] - rab_z[i] * t_x_yy_xz_xy[i];

        t_xz_yy_xz_xx[i] = t_x_yyz_xz_xx[i] - rab_z[i] * t_x_yy_xz_xx[i];

        t_xz_yy_xy_zz[i] = t_x_yyz_xy_zz[i] - rab_z[i] * t_x_yy_xy_zz[i];

        t_xz_yy_xy_yz[i] = t_x_yyz_xy_yz[i] - rab_z[i] * t_x_yy_xy_yz[i];

        t_xz_yy_xy_yy[i] = t_x_yyz_xy_yy[i] - rab_z[i] * t_x_yy_xy_yy[i];

        t_xz_yy_xy_xz[i] = t_x_yyz_xy_xz[i] - rab_z[i] * t_x_yy_xy_xz[i];

        t_xz_yy_xy_xy[i] = t_x_yyz_xy_xy[i] - rab_z[i] * t_x_yy_xy_xy[i];

        t_xz_yy_xy_xx[i] = t_x_yyz_xy_xx[i] - rab_z[i] * t_x_yy_xy_xx[i];

        t_xz_yy_xx_zz[i] = t_x_yyz_xx_zz[i] - rab_z[i] * t_x_yy_xx_zz[i];

        t_xz_yy_xx_yz[i] = t_x_yyz_xx_yz[i] - rab_z[i] * t_x_yy_xx_yz[i];

        t_xz_yy_xx_yy[i] = t_x_yyz_xx_yy[i] - rab_z[i] * t_x_yy_xx_yy[i];

        t_xz_yy_xx_xz[i] = t_x_yyz_xx_xz[i] - rab_z[i] * t_x_yy_xx_xz[i];

        t_xz_yy_xx_xy[i] = t_x_yyz_xx_xy[i] - rab_z[i] * t_x_yy_xx_xy[i];

        t_xz_yy_xx_xx[i] = t_x_yyz_xx_xx[i] - rab_z[i] * t_x_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xz_xx_xx, t_x_xz_xx_xy, t_x_xz_xx_xz, t_x_xz_xx_yy,\
                           t_x_xz_xx_yz, t_x_xz_xx_zz, t_x_xz_xy_xx, t_x_xz_xy_xy, t_x_xz_xy_xz,\
                           t_x_xz_xy_yy, t_x_xz_xy_yz, t_x_xz_xy_zz, t_x_xz_xz_xx, t_x_xz_xz_xy,\
                           t_x_xz_xz_xz, t_x_xz_xz_yy, t_x_xz_xz_yz, t_x_xz_xz_zz, t_x_xz_yy_xx,\
                           t_x_xz_yy_xy, t_x_xz_yy_xz, t_x_xz_yy_yy, t_x_xz_yy_yz, t_x_xz_yy_zz,\
                           t_x_xz_yz_xx, t_x_xz_yz_xy, t_x_xz_yz_xz, t_x_xz_yz_yy, t_x_xz_yz_yz,\
                           t_x_xz_yz_zz, t_x_xz_zz_xx, t_x_xz_zz_xy, t_x_xz_zz_xz, t_x_xz_zz_yy,\
                           t_x_xz_zz_yz, t_x_xz_zz_zz, t_x_xzz_xx_xx, t_x_xzz_xx_xy,\
                           t_x_xzz_xx_xz, t_x_xzz_xx_yy, t_x_xzz_xx_yz, t_x_xzz_xx_zz,\
                           t_x_xzz_xy_xx, t_x_xzz_xy_xy, t_x_xzz_xy_xz, t_x_xzz_xy_yy,\
                           t_x_xzz_xy_yz, t_x_xzz_xy_zz, t_x_xzz_xz_xx, t_x_xzz_xz_xy,\
                           t_x_xzz_xz_xz, t_x_xzz_xz_yy, t_x_xzz_xz_yz, t_x_xzz_xz_zz,\
                           t_x_xzz_yy_xx, t_x_xzz_yy_xy, t_x_xzz_yy_xz, t_x_xzz_yy_yy,\
                           t_x_xzz_yy_yz, t_x_xzz_yy_zz, t_x_xzz_yz_xx, t_x_xzz_yz_xy,\
                           t_x_xzz_yz_xz, t_x_xzz_yz_yy, t_x_xzz_yz_yz, t_x_xzz_yz_zz,\
                           t_x_xzz_zz_xx, t_x_xzz_zz_xy, t_x_xzz_zz_xz, t_x_xzz_zz_yy,\
                           t_x_xzz_zz_yz, t_x_xzz_zz_zz, t_xz_xz_xx_xx, t_xz_xz_xx_xy,\
                           t_xz_xz_xx_xz, t_xz_xz_xx_yy, t_xz_xz_xx_yz, t_xz_xz_xx_zz,\
                           t_xz_xz_xy_xx, t_xz_xz_xy_xy, t_xz_xz_xy_xz, t_xz_xz_xy_yy,\
                           t_xz_xz_xy_yz, t_xz_xz_xy_zz, t_xz_xz_xz_xx, t_xz_xz_xz_xy,\
                           t_xz_xz_xz_xz, t_xz_xz_xz_yy, t_xz_xz_xz_yz, t_xz_xz_xz_zz,\
                           t_xz_xz_yy_xx, t_xz_xz_yy_xy, t_xz_xz_yy_xz, t_xz_xz_yy_yy,\
                           t_xz_xz_yy_yz, t_xz_xz_yy_zz, t_xz_xz_yz_xx, t_xz_xz_yz_xy,\
                           t_xz_xz_yz_xz, t_xz_xz_yz_yy, t_xz_xz_yz_yz, t_xz_xz_yz_zz,\
                           t_xz_xz_zz_xx, t_xz_xz_zz_xy, t_xz_xz_zz_xz, t_xz_xz_zz_yy,\
                           t_xz_xz_zz_yz, t_xz_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_xz_zz_zz[i] = t_x_xzz_zz_zz[i] - rab_z[i] * t_x_xz_zz_zz[i];

        t_xz_xz_zz_yz[i] = t_x_xzz_zz_yz[i] - rab_z[i] * t_x_xz_zz_yz[i];

        t_xz_xz_zz_yy[i] = t_x_xzz_zz_yy[i] - rab_z[i] * t_x_xz_zz_yy[i];

        t_xz_xz_zz_xz[i] = t_x_xzz_zz_xz[i] - rab_z[i] * t_x_xz_zz_xz[i];

        t_xz_xz_zz_xy[i] = t_x_xzz_zz_xy[i] - rab_z[i] * t_x_xz_zz_xy[i];

        t_xz_xz_zz_xx[i] = t_x_xzz_zz_xx[i] - rab_z[i] * t_x_xz_zz_xx[i];

        t_xz_xz_yz_zz[i] = t_x_xzz_yz_zz[i] - rab_z[i] * t_x_xz_yz_zz[i];

        t_xz_xz_yz_yz[i] = t_x_xzz_yz_yz[i] - rab_z[i] * t_x_xz_yz_yz[i];

        t_xz_xz_yz_yy[i] = t_x_xzz_yz_yy[i] - rab_z[i] * t_x_xz_yz_yy[i];

        t_xz_xz_yz_xz[i] = t_x_xzz_yz_xz[i] - rab_z[i] * t_x_xz_yz_xz[i];

        t_xz_xz_yz_xy[i] = t_x_xzz_yz_xy[i] - rab_z[i] * t_x_xz_yz_xy[i];

        t_xz_xz_yz_xx[i] = t_x_xzz_yz_xx[i] - rab_z[i] * t_x_xz_yz_xx[i];

        t_xz_xz_yy_zz[i] = t_x_xzz_yy_zz[i] - rab_z[i] * t_x_xz_yy_zz[i];

        t_xz_xz_yy_yz[i] = t_x_xzz_yy_yz[i] - rab_z[i] * t_x_xz_yy_yz[i];

        t_xz_xz_yy_yy[i] = t_x_xzz_yy_yy[i] - rab_z[i] * t_x_xz_yy_yy[i];

        t_xz_xz_yy_xz[i] = t_x_xzz_yy_xz[i] - rab_z[i] * t_x_xz_yy_xz[i];

        t_xz_xz_yy_xy[i] = t_x_xzz_yy_xy[i] - rab_z[i] * t_x_xz_yy_xy[i];

        t_xz_xz_yy_xx[i] = t_x_xzz_yy_xx[i] - rab_z[i] * t_x_xz_yy_xx[i];

        t_xz_xz_xz_zz[i] = t_x_xzz_xz_zz[i] - rab_z[i] * t_x_xz_xz_zz[i];

        t_xz_xz_xz_yz[i] = t_x_xzz_xz_yz[i] - rab_z[i] * t_x_xz_xz_yz[i];

        t_xz_xz_xz_yy[i] = t_x_xzz_xz_yy[i] - rab_z[i] * t_x_xz_xz_yy[i];

        t_xz_xz_xz_xz[i] = t_x_xzz_xz_xz[i] - rab_z[i] * t_x_xz_xz_xz[i];

        t_xz_xz_xz_xy[i] = t_x_xzz_xz_xy[i] - rab_z[i] * t_x_xz_xz_xy[i];

        t_xz_xz_xz_xx[i] = t_x_xzz_xz_xx[i] - rab_z[i] * t_x_xz_xz_xx[i];

        t_xz_xz_xy_zz[i] = t_x_xzz_xy_zz[i] - rab_z[i] * t_x_xz_xy_zz[i];

        t_xz_xz_xy_yz[i] = t_x_xzz_xy_yz[i] - rab_z[i] * t_x_xz_xy_yz[i];

        t_xz_xz_xy_yy[i] = t_x_xzz_xy_yy[i] - rab_z[i] * t_x_xz_xy_yy[i];

        t_xz_xz_xy_xz[i] = t_x_xzz_xy_xz[i] - rab_z[i] * t_x_xz_xy_xz[i];

        t_xz_xz_xy_xy[i] = t_x_xzz_xy_xy[i] - rab_z[i] * t_x_xz_xy_xy[i];

        t_xz_xz_xy_xx[i] = t_x_xzz_xy_xx[i] - rab_z[i] * t_x_xz_xy_xx[i];

        t_xz_xz_xx_zz[i] = t_x_xzz_xx_zz[i] - rab_z[i] * t_x_xz_xx_zz[i];

        t_xz_xz_xx_yz[i] = t_x_xzz_xx_yz[i] - rab_z[i] * t_x_xz_xx_yz[i];

        t_xz_xz_xx_yy[i] = t_x_xzz_xx_yy[i] - rab_z[i] * t_x_xz_xx_yy[i];

        t_xz_xz_xx_xz[i] = t_x_xzz_xx_xz[i] - rab_z[i] * t_x_xz_xx_xz[i];

        t_xz_xz_xx_xy[i] = t_x_xzz_xx_xy[i] - rab_z[i] * t_x_xz_xx_xy[i];

        t_xz_xz_xx_xx[i] = t_x_xzz_xx_xx[i] - rab_z[i] * t_x_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xy_xx_xx, t_x_xy_xx_xy, t_x_xy_xx_xz, t_x_xy_xx_yy,\
                           t_x_xy_xx_yz, t_x_xy_xx_zz, t_x_xy_xy_xx, t_x_xy_xy_xy, t_x_xy_xy_xz,\
                           t_x_xy_xy_yy, t_x_xy_xy_yz, t_x_xy_xy_zz, t_x_xy_xz_xx, t_x_xy_xz_xy,\
                           t_x_xy_xz_xz, t_x_xy_xz_yy, t_x_xy_xz_yz, t_x_xy_xz_zz, t_x_xy_yy_xx,\
                           t_x_xy_yy_xy, t_x_xy_yy_xz, t_x_xy_yy_yy, t_x_xy_yy_yz, t_x_xy_yy_zz,\
                           t_x_xy_yz_xx, t_x_xy_yz_xy, t_x_xy_yz_xz, t_x_xy_yz_yy, t_x_xy_yz_yz,\
                           t_x_xy_yz_zz, t_x_xy_zz_xx, t_x_xy_zz_xy, t_x_xy_zz_xz, t_x_xy_zz_yy,\
                           t_x_xy_zz_yz, t_x_xy_zz_zz, t_x_xyz_xx_xx, t_x_xyz_xx_xy,\
                           t_x_xyz_xx_xz, t_x_xyz_xx_yy, t_x_xyz_xx_yz, t_x_xyz_xx_zz,\
                           t_x_xyz_xy_xx, t_x_xyz_xy_xy, t_x_xyz_xy_xz, t_x_xyz_xy_yy,\
                           t_x_xyz_xy_yz, t_x_xyz_xy_zz, t_x_xyz_xz_xx, t_x_xyz_xz_xy,\
                           t_x_xyz_xz_xz, t_x_xyz_xz_yy, t_x_xyz_xz_yz, t_x_xyz_xz_zz,\
                           t_x_xyz_yy_xx, t_x_xyz_yy_xy, t_x_xyz_yy_xz, t_x_xyz_yy_yy,\
                           t_x_xyz_yy_yz, t_x_xyz_yy_zz, t_x_xyz_yz_xx, t_x_xyz_yz_xy,\
                           t_x_xyz_yz_xz, t_x_xyz_yz_yy, t_x_xyz_yz_yz, t_x_xyz_yz_zz,\
                           t_x_xyz_zz_xx, t_x_xyz_zz_xy, t_x_xyz_zz_xz, t_x_xyz_zz_yy,\
                           t_x_xyz_zz_yz, t_x_xyz_zz_zz, t_xz_xy_xx_xx, t_xz_xy_xx_xy,\
                           t_xz_xy_xx_xz, t_xz_xy_xx_yy, t_xz_xy_xx_yz, t_xz_xy_xx_zz,\
                           t_xz_xy_xy_xx, t_xz_xy_xy_xy, t_xz_xy_xy_xz, t_xz_xy_xy_yy,\
                           t_xz_xy_xy_yz, t_xz_xy_xy_zz, t_xz_xy_xz_xx, t_xz_xy_xz_xy,\
                           t_xz_xy_xz_xz, t_xz_xy_xz_yy, t_xz_xy_xz_yz, t_xz_xy_xz_zz,\
                           t_xz_xy_yy_xx, t_xz_xy_yy_xy, t_xz_xy_yy_xz, t_xz_xy_yy_yy,\
                           t_xz_xy_yy_yz, t_xz_xy_yy_zz, t_xz_xy_yz_xx, t_xz_xy_yz_xy,\
                           t_xz_xy_yz_xz, t_xz_xy_yz_yy, t_xz_xy_yz_yz, t_xz_xy_yz_zz,\
                           t_xz_xy_zz_xx, t_xz_xy_zz_xy, t_xz_xy_zz_xz, t_xz_xy_zz_yy,\
                           t_xz_xy_zz_yz, t_xz_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_xy_zz_zz[i] = t_x_xyz_zz_zz[i] - rab_z[i] * t_x_xy_zz_zz[i];

        t_xz_xy_zz_yz[i] = t_x_xyz_zz_yz[i] - rab_z[i] * t_x_xy_zz_yz[i];

        t_xz_xy_zz_yy[i] = t_x_xyz_zz_yy[i] - rab_z[i] * t_x_xy_zz_yy[i];

        t_xz_xy_zz_xz[i] = t_x_xyz_zz_xz[i] - rab_z[i] * t_x_xy_zz_xz[i];

        t_xz_xy_zz_xy[i] = t_x_xyz_zz_xy[i] - rab_z[i] * t_x_xy_zz_xy[i];

        t_xz_xy_zz_xx[i] = t_x_xyz_zz_xx[i] - rab_z[i] * t_x_xy_zz_xx[i];

        t_xz_xy_yz_zz[i] = t_x_xyz_yz_zz[i] - rab_z[i] * t_x_xy_yz_zz[i];

        t_xz_xy_yz_yz[i] = t_x_xyz_yz_yz[i] - rab_z[i] * t_x_xy_yz_yz[i];

        t_xz_xy_yz_yy[i] = t_x_xyz_yz_yy[i] - rab_z[i] * t_x_xy_yz_yy[i];

        t_xz_xy_yz_xz[i] = t_x_xyz_yz_xz[i] - rab_z[i] * t_x_xy_yz_xz[i];

        t_xz_xy_yz_xy[i] = t_x_xyz_yz_xy[i] - rab_z[i] * t_x_xy_yz_xy[i];

        t_xz_xy_yz_xx[i] = t_x_xyz_yz_xx[i] - rab_z[i] * t_x_xy_yz_xx[i];

        t_xz_xy_yy_zz[i] = t_x_xyz_yy_zz[i] - rab_z[i] * t_x_xy_yy_zz[i];

        t_xz_xy_yy_yz[i] = t_x_xyz_yy_yz[i] - rab_z[i] * t_x_xy_yy_yz[i];

        t_xz_xy_yy_yy[i] = t_x_xyz_yy_yy[i] - rab_z[i] * t_x_xy_yy_yy[i];

        t_xz_xy_yy_xz[i] = t_x_xyz_yy_xz[i] - rab_z[i] * t_x_xy_yy_xz[i];

        t_xz_xy_yy_xy[i] = t_x_xyz_yy_xy[i] - rab_z[i] * t_x_xy_yy_xy[i];

        t_xz_xy_yy_xx[i] = t_x_xyz_yy_xx[i] - rab_z[i] * t_x_xy_yy_xx[i];

        t_xz_xy_xz_zz[i] = t_x_xyz_xz_zz[i] - rab_z[i] * t_x_xy_xz_zz[i];

        t_xz_xy_xz_yz[i] = t_x_xyz_xz_yz[i] - rab_z[i] * t_x_xy_xz_yz[i];

        t_xz_xy_xz_yy[i] = t_x_xyz_xz_yy[i] - rab_z[i] * t_x_xy_xz_yy[i];

        t_xz_xy_xz_xz[i] = t_x_xyz_xz_xz[i] - rab_z[i] * t_x_xy_xz_xz[i];

        t_xz_xy_xz_xy[i] = t_x_xyz_xz_xy[i] - rab_z[i] * t_x_xy_xz_xy[i];

        t_xz_xy_xz_xx[i] = t_x_xyz_xz_xx[i] - rab_z[i] * t_x_xy_xz_xx[i];

        t_xz_xy_xy_zz[i] = t_x_xyz_xy_zz[i] - rab_z[i] * t_x_xy_xy_zz[i];

        t_xz_xy_xy_yz[i] = t_x_xyz_xy_yz[i] - rab_z[i] * t_x_xy_xy_yz[i];

        t_xz_xy_xy_yy[i] = t_x_xyz_xy_yy[i] - rab_z[i] * t_x_xy_xy_yy[i];

        t_xz_xy_xy_xz[i] = t_x_xyz_xy_xz[i] - rab_z[i] * t_x_xy_xy_xz[i];

        t_xz_xy_xy_xy[i] = t_x_xyz_xy_xy[i] - rab_z[i] * t_x_xy_xy_xy[i];

        t_xz_xy_xy_xx[i] = t_x_xyz_xy_xx[i] - rab_z[i] * t_x_xy_xy_xx[i];

        t_xz_xy_xx_zz[i] = t_x_xyz_xx_zz[i] - rab_z[i] * t_x_xy_xx_zz[i];

        t_xz_xy_xx_yz[i] = t_x_xyz_xx_yz[i] - rab_z[i] * t_x_xy_xx_yz[i];

        t_xz_xy_xx_yy[i] = t_x_xyz_xx_yy[i] - rab_z[i] * t_x_xy_xx_yy[i];

        t_xz_xy_xx_xz[i] = t_x_xyz_xx_xz[i] - rab_z[i] * t_x_xy_xx_xz[i];

        t_xz_xy_xx_xy[i] = t_x_xyz_xx_xy[i] - rab_z[i] * t_x_xy_xx_xy[i];

        t_xz_xy_xx_xx[i] = t_x_xyz_xx_xx[i] - rab_z[i] * t_x_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xx_xx_xx, t_x_xx_xx_xy, t_x_xx_xx_xz, t_x_xx_xx_yy,\
                           t_x_xx_xx_yz, t_x_xx_xx_zz, t_x_xx_xy_xx, t_x_xx_xy_xy, t_x_xx_xy_xz,\
                           t_x_xx_xy_yy, t_x_xx_xy_yz, t_x_xx_xy_zz, t_x_xx_xz_xx, t_x_xx_xz_xy,\
                           t_x_xx_xz_xz, t_x_xx_xz_yy, t_x_xx_xz_yz, t_x_xx_xz_zz, t_x_xx_yy_xx,\
                           t_x_xx_yy_xy, t_x_xx_yy_xz, t_x_xx_yy_yy, t_x_xx_yy_yz, t_x_xx_yy_zz,\
                           t_x_xx_yz_xx, t_x_xx_yz_xy, t_x_xx_yz_xz, t_x_xx_yz_yy, t_x_xx_yz_yz,\
                           t_x_xx_yz_zz, t_x_xx_zz_xx, t_x_xx_zz_xy, t_x_xx_zz_xz, t_x_xx_zz_yy,\
                           t_x_xx_zz_yz, t_x_xx_zz_zz, t_x_xxz_xx_xx, t_x_xxz_xx_xy,\
                           t_x_xxz_xx_xz, t_x_xxz_xx_yy, t_x_xxz_xx_yz, t_x_xxz_xx_zz,\
                           t_x_xxz_xy_xx, t_x_xxz_xy_xy, t_x_xxz_xy_xz, t_x_xxz_xy_yy,\
                           t_x_xxz_xy_yz, t_x_xxz_xy_zz, t_x_xxz_xz_xx, t_x_xxz_xz_xy,\
                           t_x_xxz_xz_xz, t_x_xxz_xz_yy, t_x_xxz_xz_yz, t_x_xxz_xz_zz,\
                           t_x_xxz_yy_xx, t_x_xxz_yy_xy, t_x_xxz_yy_xz, t_x_xxz_yy_yy,\
                           t_x_xxz_yy_yz, t_x_xxz_yy_zz, t_x_xxz_yz_xx, t_x_xxz_yz_xy,\
                           t_x_xxz_yz_xz, t_x_xxz_yz_yy, t_x_xxz_yz_yz, t_x_xxz_yz_zz,\
                           t_x_xxz_zz_xx, t_x_xxz_zz_xy, t_x_xxz_zz_xz, t_x_xxz_zz_yy,\
                           t_x_xxz_zz_yz, t_x_xxz_zz_zz, t_xz_xx_xx_xx, t_xz_xx_xx_xy,\
                           t_xz_xx_xx_xz, t_xz_xx_xx_yy, t_xz_xx_xx_yz, t_xz_xx_xx_zz,\
                           t_xz_xx_xy_xx, t_xz_xx_xy_xy, t_xz_xx_xy_xz, t_xz_xx_xy_yy,\
                           t_xz_xx_xy_yz, t_xz_xx_xy_zz, t_xz_xx_xz_xx, t_xz_xx_xz_xy,\
                           t_xz_xx_xz_xz, t_xz_xx_xz_yy, t_xz_xx_xz_yz, t_xz_xx_xz_zz,\
                           t_xz_xx_yy_xx, t_xz_xx_yy_xy, t_xz_xx_yy_xz, t_xz_xx_yy_yy,\
                           t_xz_xx_yy_yz, t_xz_xx_yy_zz, t_xz_xx_yz_xx, t_xz_xx_yz_xy,\
                           t_xz_xx_yz_xz, t_xz_xx_yz_yy, t_xz_xx_yz_yz, t_xz_xx_yz_zz,\
                           t_xz_xx_zz_xx, t_xz_xx_zz_xy, t_xz_xx_zz_xz, t_xz_xx_zz_yy,\
                           t_xz_xx_zz_yz, t_xz_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_xx_zz_zz[i] = t_x_xxz_zz_zz[i] - rab_z[i] * t_x_xx_zz_zz[i];

        t_xz_xx_zz_yz[i] = t_x_xxz_zz_yz[i] - rab_z[i] * t_x_xx_zz_yz[i];

        t_xz_xx_zz_yy[i] = t_x_xxz_zz_yy[i] - rab_z[i] * t_x_xx_zz_yy[i];

        t_xz_xx_zz_xz[i] = t_x_xxz_zz_xz[i] - rab_z[i] * t_x_xx_zz_xz[i];

        t_xz_xx_zz_xy[i] = t_x_xxz_zz_xy[i] - rab_z[i] * t_x_xx_zz_xy[i];

        t_xz_xx_zz_xx[i] = t_x_xxz_zz_xx[i] - rab_z[i] * t_x_xx_zz_xx[i];

        t_xz_xx_yz_zz[i] = t_x_xxz_yz_zz[i] - rab_z[i] * t_x_xx_yz_zz[i];

        t_xz_xx_yz_yz[i] = t_x_xxz_yz_yz[i] - rab_z[i] * t_x_xx_yz_yz[i];

        t_xz_xx_yz_yy[i] = t_x_xxz_yz_yy[i] - rab_z[i] * t_x_xx_yz_yy[i];

        t_xz_xx_yz_xz[i] = t_x_xxz_yz_xz[i] - rab_z[i] * t_x_xx_yz_xz[i];

        t_xz_xx_yz_xy[i] = t_x_xxz_yz_xy[i] - rab_z[i] * t_x_xx_yz_xy[i];

        t_xz_xx_yz_xx[i] = t_x_xxz_yz_xx[i] - rab_z[i] * t_x_xx_yz_xx[i];

        t_xz_xx_yy_zz[i] = t_x_xxz_yy_zz[i] - rab_z[i] * t_x_xx_yy_zz[i];

        t_xz_xx_yy_yz[i] = t_x_xxz_yy_yz[i] - rab_z[i] * t_x_xx_yy_yz[i];

        t_xz_xx_yy_yy[i] = t_x_xxz_yy_yy[i] - rab_z[i] * t_x_xx_yy_yy[i];

        t_xz_xx_yy_xz[i] = t_x_xxz_yy_xz[i] - rab_z[i] * t_x_xx_yy_xz[i];

        t_xz_xx_yy_xy[i] = t_x_xxz_yy_xy[i] - rab_z[i] * t_x_xx_yy_xy[i];

        t_xz_xx_yy_xx[i] = t_x_xxz_yy_xx[i] - rab_z[i] * t_x_xx_yy_xx[i];

        t_xz_xx_xz_zz[i] = t_x_xxz_xz_zz[i] - rab_z[i] * t_x_xx_xz_zz[i];

        t_xz_xx_xz_yz[i] = t_x_xxz_xz_yz[i] - rab_z[i] * t_x_xx_xz_yz[i];

        t_xz_xx_xz_yy[i] = t_x_xxz_xz_yy[i] - rab_z[i] * t_x_xx_xz_yy[i];

        t_xz_xx_xz_xz[i] = t_x_xxz_xz_xz[i] - rab_z[i] * t_x_xx_xz_xz[i];

        t_xz_xx_xz_xy[i] = t_x_xxz_xz_xy[i] - rab_z[i] * t_x_xx_xz_xy[i];

        t_xz_xx_xz_xx[i] = t_x_xxz_xz_xx[i] - rab_z[i] * t_x_xx_xz_xx[i];

        t_xz_xx_xy_zz[i] = t_x_xxz_xy_zz[i] - rab_z[i] * t_x_xx_xy_zz[i];

        t_xz_xx_xy_yz[i] = t_x_xxz_xy_yz[i] - rab_z[i] * t_x_xx_xy_yz[i];

        t_xz_xx_xy_yy[i] = t_x_xxz_xy_yy[i] - rab_z[i] * t_x_xx_xy_yy[i];

        t_xz_xx_xy_xz[i] = t_x_xxz_xy_xz[i] - rab_z[i] * t_x_xx_xy_xz[i];

        t_xz_xx_xy_xy[i] = t_x_xxz_xy_xy[i] - rab_z[i] * t_x_xx_xy_xy[i];

        t_xz_xx_xy_xx[i] = t_x_xxz_xy_xx[i] - rab_z[i] * t_x_xx_xy_xx[i];

        t_xz_xx_xx_zz[i] = t_x_xxz_xx_zz[i] - rab_z[i] * t_x_xx_xx_zz[i];

        t_xz_xx_xx_yz[i] = t_x_xxz_xx_yz[i] - rab_z[i] * t_x_xx_xx_yz[i];

        t_xz_xx_xx_yy[i] = t_x_xxz_xx_yy[i] - rab_z[i] * t_x_xx_xx_yy[i];

        t_xz_xx_xx_xz[i] = t_x_xxz_xx_xz[i] - rab_z[i] * t_x_xx_xx_xz[i];

        t_xz_xx_xx_xy[i] = t_x_xxz_xx_xy[i] - rab_z[i] * t_x_xx_xx_xy[i];

        t_xz_xx_xx_xx[i] = t_x_xxz_xx_xx[i] - rab_z[i] * t_x_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_yzz_xx_xx, t_x_yzz_xx_xy, t_x_yzz_xx_xz, t_x_yzz_xx_yy,\
                           t_x_yzz_xx_yz, t_x_yzz_xx_zz, t_x_yzz_xy_xx, t_x_yzz_xy_xy,\
                           t_x_yzz_xy_xz, t_x_yzz_xy_yy, t_x_yzz_xy_yz, t_x_yzz_xy_zz,\
                           t_x_yzz_xz_xx, t_x_yzz_xz_xy, t_x_yzz_xz_xz, t_x_yzz_xz_yy,\
                           t_x_yzz_xz_yz, t_x_yzz_xz_zz, t_x_yzz_yy_xx, t_x_yzz_yy_xy,\
                           t_x_yzz_yy_xz, t_x_yzz_yy_yy, t_x_yzz_yy_yz, t_x_yzz_yy_zz,\
                           t_x_yzz_yz_xx, t_x_yzz_yz_xy, t_x_yzz_yz_xz, t_x_yzz_yz_yy,\
                           t_x_yzz_yz_yz, t_x_yzz_yz_zz, t_x_yzz_zz_xx, t_x_yzz_zz_xy,\
                           t_x_yzz_zz_xz, t_x_yzz_zz_yy, t_x_yzz_zz_yz, t_x_yzz_zz_zz,\
                           t_x_zz_xx_xx, t_x_zz_xx_xy, t_x_zz_xx_xz, t_x_zz_xx_yy, t_x_zz_xx_yz,\
                           t_x_zz_xx_zz, t_x_zz_xy_xx, t_x_zz_xy_xy, t_x_zz_xy_xz, t_x_zz_xy_yy,\
                           t_x_zz_xy_yz, t_x_zz_xy_zz, t_x_zz_xz_xx, t_x_zz_xz_xy, t_x_zz_xz_xz,\
                           t_x_zz_xz_yy, t_x_zz_xz_yz, t_x_zz_xz_zz, t_x_zz_yy_xx, t_x_zz_yy_xy,\
                           t_x_zz_yy_xz, t_x_zz_yy_yy, t_x_zz_yy_yz, t_x_zz_yy_zz, t_x_zz_yz_xx,\
                           t_x_zz_yz_xy, t_x_zz_yz_xz, t_x_zz_yz_yy, t_x_zz_yz_yz, t_x_zz_yz_zz,\
                           t_x_zz_zz_xx, t_x_zz_zz_xy, t_x_zz_zz_xz, t_x_zz_zz_yy, t_x_zz_zz_yz,\
                           t_x_zz_zz_zz, t_xy_zz_xx_xx, t_xy_zz_xx_xy, t_xy_zz_xx_xz,\
                           t_xy_zz_xx_yy, t_xy_zz_xx_yz, t_xy_zz_xx_zz, t_xy_zz_xy_xx,\
                           t_xy_zz_xy_xy, t_xy_zz_xy_xz, t_xy_zz_xy_yy, t_xy_zz_xy_yz,\
                           t_xy_zz_xy_zz, t_xy_zz_xz_xx, t_xy_zz_xz_xy, t_xy_zz_xz_xz,\
                           t_xy_zz_xz_yy, t_xy_zz_xz_yz, t_xy_zz_xz_zz, t_xy_zz_yy_xx,\
                           t_xy_zz_yy_xy, t_xy_zz_yy_xz, t_xy_zz_yy_yy, t_xy_zz_yy_yz,\
                           t_xy_zz_yy_zz, t_xy_zz_yz_xx, t_xy_zz_yz_xy, t_xy_zz_yz_xz,\
                           t_xy_zz_yz_yy, t_xy_zz_yz_yz, t_xy_zz_yz_zz, t_xy_zz_zz_xx,\
                           t_xy_zz_zz_xy, t_xy_zz_zz_xz, t_xy_zz_zz_yy, t_xy_zz_zz_yz,\
                           t_xy_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_zz_zz_zz[i] = t_x_yzz_zz_zz[i] - rab_y[i] * t_x_zz_zz_zz[i];

        t_xy_zz_zz_yz[i] = t_x_yzz_zz_yz[i] - rab_y[i] * t_x_zz_zz_yz[i];

        t_xy_zz_zz_yy[i] = t_x_yzz_zz_yy[i] - rab_y[i] * t_x_zz_zz_yy[i];

        t_xy_zz_zz_xz[i] = t_x_yzz_zz_xz[i] - rab_y[i] * t_x_zz_zz_xz[i];

        t_xy_zz_zz_xy[i] = t_x_yzz_zz_xy[i] - rab_y[i] * t_x_zz_zz_xy[i];

        t_xy_zz_zz_xx[i] = t_x_yzz_zz_xx[i] - rab_y[i] * t_x_zz_zz_xx[i];

        t_xy_zz_yz_zz[i] = t_x_yzz_yz_zz[i] - rab_y[i] * t_x_zz_yz_zz[i];

        t_xy_zz_yz_yz[i] = t_x_yzz_yz_yz[i] - rab_y[i] * t_x_zz_yz_yz[i];

        t_xy_zz_yz_yy[i] = t_x_yzz_yz_yy[i] - rab_y[i] * t_x_zz_yz_yy[i];

        t_xy_zz_yz_xz[i] = t_x_yzz_yz_xz[i] - rab_y[i] * t_x_zz_yz_xz[i];

        t_xy_zz_yz_xy[i] = t_x_yzz_yz_xy[i] - rab_y[i] * t_x_zz_yz_xy[i];

        t_xy_zz_yz_xx[i] = t_x_yzz_yz_xx[i] - rab_y[i] * t_x_zz_yz_xx[i];

        t_xy_zz_yy_zz[i] = t_x_yzz_yy_zz[i] - rab_y[i] * t_x_zz_yy_zz[i];

        t_xy_zz_yy_yz[i] = t_x_yzz_yy_yz[i] - rab_y[i] * t_x_zz_yy_yz[i];

        t_xy_zz_yy_yy[i] = t_x_yzz_yy_yy[i] - rab_y[i] * t_x_zz_yy_yy[i];

        t_xy_zz_yy_xz[i] = t_x_yzz_yy_xz[i] - rab_y[i] * t_x_zz_yy_xz[i];

        t_xy_zz_yy_xy[i] = t_x_yzz_yy_xy[i] - rab_y[i] * t_x_zz_yy_xy[i];

        t_xy_zz_yy_xx[i] = t_x_yzz_yy_xx[i] - rab_y[i] * t_x_zz_yy_xx[i];

        t_xy_zz_xz_zz[i] = t_x_yzz_xz_zz[i] - rab_y[i] * t_x_zz_xz_zz[i];

        t_xy_zz_xz_yz[i] = t_x_yzz_xz_yz[i] - rab_y[i] * t_x_zz_xz_yz[i];

        t_xy_zz_xz_yy[i] = t_x_yzz_xz_yy[i] - rab_y[i] * t_x_zz_xz_yy[i];

        t_xy_zz_xz_xz[i] = t_x_yzz_xz_xz[i] - rab_y[i] * t_x_zz_xz_xz[i];

        t_xy_zz_xz_xy[i] = t_x_yzz_xz_xy[i] - rab_y[i] * t_x_zz_xz_xy[i];

        t_xy_zz_xz_xx[i] = t_x_yzz_xz_xx[i] - rab_y[i] * t_x_zz_xz_xx[i];

        t_xy_zz_xy_zz[i] = t_x_yzz_xy_zz[i] - rab_y[i] * t_x_zz_xy_zz[i];

        t_xy_zz_xy_yz[i] = t_x_yzz_xy_yz[i] - rab_y[i] * t_x_zz_xy_yz[i];

        t_xy_zz_xy_yy[i] = t_x_yzz_xy_yy[i] - rab_y[i] * t_x_zz_xy_yy[i];

        t_xy_zz_xy_xz[i] = t_x_yzz_xy_xz[i] - rab_y[i] * t_x_zz_xy_xz[i];

        t_xy_zz_xy_xy[i] = t_x_yzz_xy_xy[i] - rab_y[i] * t_x_zz_xy_xy[i];

        t_xy_zz_xy_xx[i] = t_x_yzz_xy_xx[i] - rab_y[i] * t_x_zz_xy_xx[i];

        t_xy_zz_xx_zz[i] = t_x_yzz_xx_zz[i] - rab_y[i] * t_x_zz_xx_zz[i];

        t_xy_zz_xx_yz[i] = t_x_yzz_xx_yz[i] - rab_y[i] * t_x_zz_xx_yz[i];

        t_xy_zz_xx_yy[i] = t_x_yzz_xx_yy[i] - rab_y[i] * t_x_zz_xx_yy[i];

        t_xy_zz_xx_xz[i] = t_x_yzz_xx_xz[i] - rab_y[i] * t_x_zz_xx_xz[i];

        t_xy_zz_xx_xy[i] = t_x_yzz_xx_xy[i] - rab_y[i] * t_x_zz_xx_xy[i];

        t_xy_zz_xx_xx[i] = t_x_yzz_xx_xx[i] - rab_y[i] * t_x_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_yyz_xx_xx, t_x_yyz_xx_xy, t_x_yyz_xx_xz, t_x_yyz_xx_yy,\
                           t_x_yyz_xx_yz, t_x_yyz_xx_zz, t_x_yyz_xy_xx, t_x_yyz_xy_xy,\
                           t_x_yyz_xy_xz, t_x_yyz_xy_yy, t_x_yyz_xy_yz, t_x_yyz_xy_zz,\
                           t_x_yyz_xz_xx, t_x_yyz_xz_xy, t_x_yyz_xz_xz, t_x_yyz_xz_yy,\
                           t_x_yyz_xz_yz, t_x_yyz_xz_zz, t_x_yyz_yy_xx, t_x_yyz_yy_xy,\
                           t_x_yyz_yy_xz, t_x_yyz_yy_yy, t_x_yyz_yy_yz, t_x_yyz_yy_zz,\
                           t_x_yyz_yz_xx, t_x_yyz_yz_xy, t_x_yyz_yz_xz, t_x_yyz_yz_yy,\
                           t_x_yyz_yz_yz, t_x_yyz_yz_zz, t_x_yyz_zz_xx, t_x_yyz_zz_xy,\
                           t_x_yyz_zz_xz, t_x_yyz_zz_yy, t_x_yyz_zz_yz, t_x_yyz_zz_zz,\
                           t_x_yz_xx_xx, t_x_yz_xx_xy, t_x_yz_xx_xz, t_x_yz_xx_yy, t_x_yz_xx_yz,\
                           t_x_yz_xx_zz, t_x_yz_xy_xx, t_x_yz_xy_xy, t_x_yz_xy_xz, t_x_yz_xy_yy,\
                           t_x_yz_xy_yz, t_x_yz_xy_zz, t_x_yz_xz_xx, t_x_yz_xz_xy, t_x_yz_xz_xz,\
                           t_x_yz_xz_yy, t_x_yz_xz_yz, t_x_yz_xz_zz, t_x_yz_yy_xx, t_x_yz_yy_xy,\
                           t_x_yz_yy_xz, t_x_yz_yy_yy, t_x_yz_yy_yz, t_x_yz_yy_zz, t_x_yz_yz_xx,\
                           t_x_yz_yz_xy, t_x_yz_yz_xz, t_x_yz_yz_yy, t_x_yz_yz_yz, t_x_yz_yz_zz,\
                           t_x_yz_zz_xx, t_x_yz_zz_xy, t_x_yz_zz_xz, t_x_yz_zz_yy, t_x_yz_zz_yz,\
                           t_x_yz_zz_zz, t_xy_yz_xx_xx, t_xy_yz_xx_xy, t_xy_yz_xx_xz,\
                           t_xy_yz_xx_yy, t_xy_yz_xx_yz, t_xy_yz_xx_zz, t_xy_yz_xy_xx,\
                           t_xy_yz_xy_xy, t_xy_yz_xy_xz, t_xy_yz_xy_yy, t_xy_yz_xy_yz,\
                           t_xy_yz_xy_zz, t_xy_yz_xz_xx, t_xy_yz_xz_xy, t_xy_yz_xz_xz,\
                           t_xy_yz_xz_yy, t_xy_yz_xz_yz, t_xy_yz_xz_zz, t_xy_yz_yy_xx,\
                           t_xy_yz_yy_xy, t_xy_yz_yy_xz, t_xy_yz_yy_yy, t_xy_yz_yy_yz,\
                           t_xy_yz_yy_zz, t_xy_yz_yz_xx, t_xy_yz_yz_xy, t_xy_yz_yz_xz,\
                           t_xy_yz_yz_yy, t_xy_yz_yz_yz, t_xy_yz_yz_zz, t_xy_yz_zz_xx,\
                           t_xy_yz_zz_xy, t_xy_yz_zz_xz, t_xy_yz_zz_yy, t_xy_yz_zz_yz,\
                           t_xy_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_yz_zz_zz[i] = t_x_yyz_zz_zz[i] - rab_y[i] * t_x_yz_zz_zz[i];

        t_xy_yz_zz_yz[i] = t_x_yyz_zz_yz[i] - rab_y[i] * t_x_yz_zz_yz[i];

        t_xy_yz_zz_yy[i] = t_x_yyz_zz_yy[i] - rab_y[i] * t_x_yz_zz_yy[i];

        t_xy_yz_zz_xz[i] = t_x_yyz_zz_xz[i] - rab_y[i] * t_x_yz_zz_xz[i];

        t_xy_yz_zz_xy[i] = t_x_yyz_zz_xy[i] - rab_y[i] * t_x_yz_zz_xy[i];

        t_xy_yz_zz_xx[i] = t_x_yyz_zz_xx[i] - rab_y[i] * t_x_yz_zz_xx[i];

        t_xy_yz_yz_zz[i] = t_x_yyz_yz_zz[i] - rab_y[i] * t_x_yz_yz_zz[i];

        t_xy_yz_yz_yz[i] = t_x_yyz_yz_yz[i] - rab_y[i] * t_x_yz_yz_yz[i];

        t_xy_yz_yz_yy[i] = t_x_yyz_yz_yy[i] - rab_y[i] * t_x_yz_yz_yy[i];

        t_xy_yz_yz_xz[i] = t_x_yyz_yz_xz[i] - rab_y[i] * t_x_yz_yz_xz[i];

        t_xy_yz_yz_xy[i] = t_x_yyz_yz_xy[i] - rab_y[i] * t_x_yz_yz_xy[i];

        t_xy_yz_yz_xx[i] = t_x_yyz_yz_xx[i] - rab_y[i] * t_x_yz_yz_xx[i];

        t_xy_yz_yy_zz[i] = t_x_yyz_yy_zz[i] - rab_y[i] * t_x_yz_yy_zz[i];

        t_xy_yz_yy_yz[i] = t_x_yyz_yy_yz[i] - rab_y[i] * t_x_yz_yy_yz[i];

        t_xy_yz_yy_yy[i] = t_x_yyz_yy_yy[i] - rab_y[i] * t_x_yz_yy_yy[i];

        t_xy_yz_yy_xz[i] = t_x_yyz_yy_xz[i] - rab_y[i] * t_x_yz_yy_xz[i];

        t_xy_yz_yy_xy[i] = t_x_yyz_yy_xy[i] - rab_y[i] * t_x_yz_yy_xy[i];

        t_xy_yz_yy_xx[i] = t_x_yyz_yy_xx[i] - rab_y[i] * t_x_yz_yy_xx[i];

        t_xy_yz_xz_zz[i] = t_x_yyz_xz_zz[i] - rab_y[i] * t_x_yz_xz_zz[i];

        t_xy_yz_xz_yz[i] = t_x_yyz_xz_yz[i] - rab_y[i] * t_x_yz_xz_yz[i];

        t_xy_yz_xz_yy[i] = t_x_yyz_xz_yy[i] - rab_y[i] * t_x_yz_xz_yy[i];

        t_xy_yz_xz_xz[i] = t_x_yyz_xz_xz[i] - rab_y[i] * t_x_yz_xz_xz[i];

        t_xy_yz_xz_xy[i] = t_x_yyz_xz_xy[i] - rab_y[i] * t_x_yz_xz_xy[i];

        t_xy_yz_xz_xx[i] = t_x_yyz_xz_xx[i] - rab_y[i] * t_x_yz_xz_xx[i];

        t_xy_yz_xy_zz[i] = t_x_yyz_xy_zz[i] - rab_y[i] * t_x_yz_xy_zz[i];

        t_xy_yz_xy_yz[i] = t_x_yyz_xy_yz[i] - rab_y[i] * t_x_yz_xy_yz[i];

        t_xy_yz_xy_yy[i] = t_x_yyz_xy_yy[i] - rab_y[i] * t_x_yz_xy_yy[i];

        t_xy_yz_xy_xz[i] = t_x_yyz_xy_xz[i] - rab_y[i] * t_x_yz_xy_xz[i];

        t_xy_yz_xy_xy[i] = t_x_yyz_xy_xy[i] - rab_y[i] * t_x_yz_xy_xy[i];

        t_xy_yz_xy_xx[i] = t_x_yyz_xy_xx[i] - rab_y[i] * t_x_yz_xy_xx[i];

        t_xy_yz_xx_zz[i] = t_x_yyz_xx_zz[i] - rab_y[i] * t_x_yz_xx_zz[i];

        t_xy_yz_xx_yz[i] = t_x_yyz_xx_yz[i] - rab_y[i] * t_x_yz_xx_yz[i];

        t_xy_yz_xx_yy[i] = t_x_yyz_xx_yy[i] - rab_y[i] * t_x_yz_xx_yy[i];

        t_xy_yz_xx_xz[i] = t_x_yyz_xx_xz[i] - rab_y[i] * t_x_yz_xx_xz[i];

        t_xy_yz_xx_xy[i] = t_x_yyz_xx_xy[i] - rab_y[i] * t_x_yz_xx_xy[i];

        t_xy_yz_xx_xx[i] = t_x_yyz_xx_xx[i] - rab_y[i] * t_x_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_yy_xx_xx, t_x_yy_xx_xy, t_x_yy_xx_xz, t_x_yy_xx_yy,\
                           t_x_yy_xx_yz, t_x_yy_xx_zz, t_x_yy_xy_xx, t_x_yy_xy_xy, t_x_yy_xy_xz,\
                           t_x_yy_xy_yy, t_x_yy_xy_yz, t_x_yy_xy_zz, t_x_yy_xz_xx, t_x_yy_xz_xy,\
                           t_x_yy_xz_xz, t_x_yy_xz_yy, t_x_yy_xz_yz, t_x_yy_xz_zz, t_x_yy_yy_xx,\
                           t_x_yy_yy_xy, t_x_yy_yy_xz, t_x_yy_yy_yy, t_x_yy_yy_yz, t_x_yy_yy_zz,\
                           t_x_yy_yz_xx, t_x_yy_yz_xy, t_x_yy_yz_xz, t_x_yy_yz_yy, t_x_yy_yz_yz,\
                           t_x_yy_yz_zz, t_x_yy_zz_xx, t_x_yy_zz_xy, t_x_yy_zz_xz, t_x_yy_zz_yy,\
                           t_x_yy_zz_yz, t_x_yy_zz_zz, t_x_yyy_xx_xx, t_x_yyy_xx_xy,\
                           t_x_yyy_xx_xz, t_x_yyy_xx_yy, t_x_yyy_xx_yz, t_x_yyy_xx_zz,\
                           t_x_yyy_xy_xx, t_x_yyy_xy_xy, t_x_yyy_xy_xz, t_x_yyy_xy_yy,\
                           t_x_yyy_xy_yz, t_x_yyy_xy_zz, t_x_yyy_xz_xx, t_x_yyy_xz_xy,\
                           t_x_yyy_xz_xz, t_x_yyy_xz_yy, t_x_yyy_xz_yz, t_x_yyy_xz_zz,\
                           t_x_yyy_yy_xx, t_x_yyy_yy_xy, t_x_yyy_yy_xz, t_x_yyy_yy_yy,\
                           t_x_yyy_yy_yz, t_x_yyy_yy_zz, t_x_yyy_yz_xx, t_x_yyy_yz_xy,\
                           t_x_yyy_yz_xz, t_x_yyy_yz_yy, t_x_yyy_yz_yz, t_x_yyy_yz_zz,\
                           t_x_yyy_zz_xx, t_x_yyy_zz_xy, t_x_yyy_zz_xz, t_x_yyy_zz_yy,\
                           t_x_yyy_zz_yz, t_x_yyy_zz_zz, t_xy_yy_xx_xx, t_xy_yy_xx_xy,\
                           t_xy_yy_xx_xz, t_xy_yy_xx_yy, t_xy_yy_xx_yz, t_xy_yy_xx_zz,\
                           t_xy_yy_xy_xx, t_xy_yy_xy_xy, t_xy_yy_xy_xz, t_xy_yy_xy_yy,\
                           t_xy_yy_xy_yz, t_xy_yy_xy_zz, t_xy_yy_xz_xx, t_xy_yy_xz_xy,\
                           t_xy_yy_xz_xz, t_xy_yy_xz_yy, t_xy_yy_xz_yz, t_xy_yy_xz_zz,\
                           t_xy_yy_yy_xx, t_xy_yy_yy_xy, t_xy_yy_yy_xz, t_xy_yy_yy_yy,\
                           t_xy_yy_yy_yz, t_xy_yy_yy_zz, t_xy_yy_yz_xx, t_xy_yy_yz_xy,\
                           t_xy_yy_yz_xz, t_xy_yy_yz_yy, t_xy_yy_yz_yz, t_xy_yy_yz_zz,\
                           t_xy_yy_zz_xx, t_xy_yy_zz_xy, t_xy_yy_zz_xz, t_xy_yy_zz_yy,\
                           t_xy_yy_zz_yz, t_xy_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_yy_zz_zz[i] = t_x_yyy_zz_zz[i] - rab_y[i] * t_x_yy_zz_zz[i];

        t_xy_yy_zz_yz[i] = t_x_yyy_zz_yz[i] - rab_y[i] * t_x_yy_zz_yz[i];

        t_xy_yy_zz_yy[i] = t_x_yyy_zz_yy[i] - rab_y[i] * t_x_yy_zz_yy[i];

        t_xy_yy_zz_xz[i] = t_x_yyy_zz_xz[i] - rab_y[i] * t_x_yy_zz_xz[i];

        t_xy_yy_zz_xy[i] = t_x_yyy_zz_xy[i] - rab_y[i] * t_x_yy_zz_xy[i];

        t_xy_yy_zz_xx[i] = t_x_yyy_zz_xx[i] - rab_y[i] * t_x_yy_zz_xx[i];

        t_xy_yy_yz_zz[i] = t_x_yyy_yz_zz[i] - rab_y[i] * t_x_yy_yz_zz[i];

        t_xy_yy_yz_yz[i] = t_x_yyy_yz_yz[i] - rab_y[i] * t_x_yy_yz_yz[i];

        t_xy_yy_yz_yy[i] = t_x_yyy_yz_yy[i] - rab_y[i] * t_x_yy_yz_yy[i];

        t_xy_yy_yz_xz[i] = t_x_yyy_yz_xz[i] - rab_y[i] * t_x_yy_yz_xz[i];

        t_xy_yy_yz_xy[i] = t_x_yyy_yz_xy[i] - rab_y[i] * t_x_yy_yz_xy[i];

        t_xy_yy_yz_xx[i] = t_x_yyy_yz_xx[i] - rab_y[i] * t_x_yy_yz_xx[i];

        t_xy_yy_yy_zz[i] = t_x_yyy_yy_zz[i] - rab_y[i] * t_x_yy_yy_zz[i];

        t_xy_yy_yy_yz[i] = t_x_yyy_yy_yz[i] - rab_y[i] * t_x_yy_yy_yz[i];

        t_xy_yy_yy_yy[i] = t_x_yyy_yy_yy[i] - rab_y[i] * t_x_yy_yy_yy[i];

        t_xy_yy_yy_xz[i] = t_x_yyy_yy_xz[i] - rab_y[i] * t_x_yy_yy_xz[i];

        t_xy_yy_yy_xy[i] = t_x_yyy_yy_xy[i] - rab_y[i] * t_x_yy_yy_xy[i];

        t_xy_yy_yy_xx[i] = t_x_yyy_yy_xx[i] - rab_y[i] * t_x_yy_yy_xx[i];

        t_xy_yy_xz_zz[i] = t_x_yyy_xz_zz[i] - rab_y[i] * t_x_yy_xz_zz[i];

        t_xy_yy_xz_yz[i] = t_x_yyy_xz_yz[i] - rab_y[i] * t_x_yy_xz_yz[i];

        t_xy_yy_xz_yy[i] = t_x_yyy_xz_yy[i] - rab_y[i] * t_x_yy_xz_yy[i];

        t_xy_yy_xz_xz[i] = t_x_yyy_xz_xz[i] - rab_y[i] * t_x_yy_xz_xz[i];

        t_xy_yy_xz_xy[i] = t_x_yyy_xz_xy[i] - rab_y[i] * t_x_yy_xz_xy[i];

        t_xy_yy_xz_xx[i] = t_x_yyy_xz_xx[i] - rab_y[i] * t_x_yy_xz_xx[i];

        t_xy_yy_xy_zz[i] = t_x_yyy_xy_zz[i] - rab_y[i] * t_x_yy_xy_zz[i];

        t_xy_yy_xy_yz[i] = t_x_yyy_xy_yz[i] - rab_y[i] * t_x_yy_xy_yz[i];

        t_xy_yy_xy_yy[i] = t_x_yyy_xy_yy[i] - rab_y[i] * t_x_yy_xy_yy[i];

        t_xy_yy_xy_xz[i] = t_x_yyy_xy_xz[i] - rab_y[i] * t_x_yy_xy_xz[i];

        t_xy_yy_xy_xy[i] = t_x_yyy_xy_xy[i] - rab_y[i] * t_x_yy_xy_xy[i];

        t_xy_yy_xy_xx[i] = t_x_yyy_xy_xx[i] - rab_y[i] * t_x_yy_xy_xx[i];

        t_xy_yy_xx_zz[i] = t_x_yyy_xx_zz[i] - rab_y[i] * t_x_yy_xx_zz[i];

        t_xy_yy_xx_yz[i] = t_x_yyy_xx_yz[i] - rab_y[i] * t_x_yy_xx_yz[i];

        t_xy_yy_xx_yy[i] = t_x_yyy_xx_yy[i] - rab_y[i] * t_x_yy_xx_yy[i];

        t_xy_yy_xx_xz[i] = t_x_yyy_xx_xz[i] - rab_y[i] * t_x_yy_xx_xz[i];

        t_xy_yy_xx_xy[i] = t_x_yyy_xx_xy[i] - rab_y[i] * t_x_yy_xx_xy[i];

        t_xy_yy_xx_xx[i] = t_x_yyy_xx_xx[i] - rab_y[i] * t_x_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xyz_xx_xx, t_x_xyz_xx_xy, t_x_xyz_xx_xz, t_x_xyz_xx_yy,\
                           t_x_xyz_xx_yz, t_x_xyz_xx_zz, t_x_xyz_xy_xx, t_x_xyz_xy_xy,\
                           t_x_xyz_xy_xz, t_x_xyz_xy_yy, t_x_xyz_xy_yz, t_x_xyz_xy_zz,\
                           t_x_xyz_xz_xx, t_x_xyz_xz_xy, t_x_xyz_xz_xz, t_x_xyz_xz_yy,\
                           t_x_xyz_xz_yz, t_x_xyz_xz_zz, t_x_xyz_yy_xx, t_x_xyz_yy_xy,\
                           t_x_xyz_yy_xz, t_x_xyz_yy_yy, t_x_xyz_yy_yz, t_x_xyz_yy_zz,\
                           t_x_xyz_yz_xx, t_x_xyz_yz_xy, t_x_xyz_yz_xz, t_x_xyz_yz_yy,\
                           t_x_xyz_yz_yz, t_x_xyz_yz_zz, t_x_xyz_zz_xx, t_x_xyz_zz_xy,\
                           t_x_xyz_zz_xz, t_x_xyz_zz_yy, t_x_xyz_zz_yz, t_x_xyz_zz_zz,\
                           t_x_xz_xx_xx, t_x_xz_xx_xy, t_x_xz_xx_xz, t_x_xz_xx_yy, t_x_xz_xx_yz,\
                           t_x_xz_xx_zz, t_x_xz_xy_xx, t_x_xz_xy_xy, t_x_xz_xy_xz, t_x_xz_xy_yy,\
                           t_x_xz_xy_yz, t_x_xz_xy_zz, t_x_xz_xz_xx, t_x_xz_xz_xy, t_x_xz_xz_xz,\
                           t_x_xz_xz_yy, t_x_xz_xz_yz, t_x_xz_xz_zz, t_x_xz_yy_xx, t_x_xz_yy_xy,\
                           t_x_xz_yy_xz, t_x_xz_yy_yy, t_x_xz_yy_yz, t_x_xz_yy_zz, t_x_xz_yz_xx,\
                           t_x_xz_yz_xy, t_x_xz_yz_xz, t_x_xz_yz_yy, t_x_xz_yz_yz, t_x_xz_yz_zz,\
                           t_x_xz_zz_xx, t_x_xz_zz_xy, t_x_xz_zz_xz, t_x_xz_zz_yy, t_x_xz_zz_yz,\
                           t_x_xz_zz_zz, t_xy_xz_xx_xx, t_xy_xz_xx_xy, t_xy_xz_xx_xz,\
                           t_xy_xz_xx_yy, t_xy_xz_xx_yz, t_xy_xz_xx_zz, t_xy_xz_xy_xx,\
                           t_xy_xz_xy_xy, t_xy_xz_xy_xz, t_xy_xz_xy_yy, t_xy_xz_xy_yz,\
                           t_xy_xz_xy_zz, t_xy_xz_xz_xx, t_xy_xz_xz_xy, t_xy_xz_xz_xz,\
                           t_xy_xz_xz_yy, t_xy_xz_xz_yz, t_xy_xz_xz_zz, t_xy_xz_yy_xx,\
                           t_xy_xz_yy_xy, t_xy_xz_yy_xz, t_xy_xz_yy_yy, t_xy_xz_yy_yz,\
                           t_xy_xz_yy_zz, t_xy_xz_yz_xx, t_xy_xz_yz_xy, t_xy_xz_yz_xz,\
                           t_xy_xz_yz_yy, t_xy_xz_yz_yz, t_xy_xz_yz_zz, t_xy_xz_zz_xx,\
                           t_xy_xz_zz_xy, t_xy_xz_zz_xz, t_xy_xz_zz_yy, t_xy_xz_zz_yz,\
                           t_xy_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_xz_zz_zz[i] = t_x_xyz_zz_zz[i] - rab_y[i] * t_x_xz_zz_zz[i];

        t_xy_xz_zz_yz[i] = t_x_xyz_zz_yz[i] - rab_y[i] * t_x_xz_zz_yz[i];

        t_xy_xz_zz_yy[i] = t_x_xyz_zz_yy[i] - rab_y[i] * t_x_xz_zz_yy[i];

        t_xy_xz_zz_xz[i] = t_x_xyz_zz_xz[i] - rab_y[i] * t_x_xz_zz_xz[i];

        t_xy_xz_zz_xy[i] = t_x_xyz_zz_xy[i] - rab_y[i] * t_x_xz_zz_xy[i];

        t_xy_xz_zz_xx[i] = t_x_xyz_zz_xx[i] - rab_y[i] * t_x_xz_zz_xx[i];

        t_xy_xz_yz_zz[i] = t_x_xyz_yz_zz[i] - rab_y[i] * t_x_xz_yz_zz[i];

        t_xy_xz_yz_yz[i] = t_x_xyz_yz_yz[i] - rab_y[i] * t_x_xz_yz_yz[i];

        t_xy_xz_yz_yy[i] = t_x_xyz_yz_yy[i] - rab_y[i] * t_x_xz_yz_yy[i];

        t_xy_xz_yz_xz[i] = t_x_xyz_yz_xz[i] - rab_y[i] * t_x_xz_yz_xz[i];

        t_xy_xz_yz_xy[i] = t_x_xyz_yz_xy[i] - rab_y[i] * t_x_xz_yz_xy[i];

        t_xy_xz_yz_xx[i] = t_x_xyz_yz_xx[i] - rab_y[i] * t_x_xz_yz_xx[i];

        t_xy_xz_yy_zz[i] = t_x_xyz_yy_zz[i] - rab_y[i] * t_x_xz_yy_zz[i];

        t_xy_xz_yy_yz[i] = t_x_xyz_yy_yz[i] - rab_y[i] * t_x_xz_yy_yz[i];

        t_xy_xz_yy_yy[i] = t_x_xyz_yy_yy[i] - rab_y[i] * t_x_xz_yy_yy[i];

        t_xy_xz_yy_xz[i] = t_x_xyz_yy_xz[i] - rab_y[i] * t_x_xz_yy_xz[i];

        t_xy_xz_yy_xy[i] = t_x_xyz_yy_xy[i] - rab_y[i] * t_x_xz_yy_xy[i];

        t_xy_xz_yy_xx[i] = t_x_xyz_yy_xx[i] - rab_y[i] * t_x_xz_yy_xx[i];

        t_xy_xz_xz_zz[i] = t_x_xyz_xz_zz[i] - rab_y[i] * t_x_xz_xz_zz[i];

        t_xy_xz_xz_yz[i] = t_x_xyz_xz_yz[i] - rab_y[i] * t_x_xz_xz_yz[i];

        t_xy_xz_xz_yy[i] = t_x_xyz_xz_yy[i] - rab_y[i] * t_x_xz_xz_yy[i];

        t_xy_xz_xz_xz[i] = t_x_xyz_xz_xz[i] - rab_y[i] * t_x_xz_xz_xz[i];

        t_xy_xz_xz_xy[i] = t_x_xyz_xz_xy[i] - rab_y[i] * t_x_xz_xz_xy[i];

        t_xy_xz_xz_xx[i] = t_x_xyz_xz_xx[i] - rab_y[i] * t_x_xz_xz_xx[i];

        t_xy_xz_xy_zz[i] = t_x_xyz_xy_zz[i] - rab_y[i] * t_x_xz_xy_zz[i];

        t_xy_xz_xy_yz[i] = t_x_xyz_xy_yz[i] - rab_y[i] * t_x_xz_xy_yz[i];

        t_xy_xz_xy_yy[i] = t_x_xyz_xy_yy[i] - rab_y[i] * t_x_xz_xy_yy[i];

        t_xy_xz_xy_xz[i] = t_x_xyz_xy_xz[i] - rab_y[i] * t_x_xz_xy_xz[i];

        t_xy_xz_xy_xy[i] = t_x_xyz_xy_xy[i] - rab_y[i] * t_x_xz_xy_xy[i];

        t_xy_xz_xy_xx[i] = t_x_xyz_xy_xx[i] - rab_y[i] * t_x_xz_xy_xx[i];

        t_xy_xz_xx_zz[i] = t_x_xyz_xx_zz[i] - rab_y[i] * t_x_xz_xx_zz[i];

        t_xy_xz_xx_yz[i] = t_x_xyz_xx_yz[i] - rab_y[i] * t_x_xz_xx_yz[i];

        t_xy_xz_xx_yy[i] = t_x_xyz_xx_yy[i] - rab_y[i] * t_x_xz_xx_yy[i];

        t_xy_xz_xx_xz[i] = t_x_xyz_xx_xz[i] - rab_y[i] * t_x_xz_xx_xz[i];

        t_xy_xz_xx_xy[i] = t_x_xyz_xx_xy[i] - rab_y[i] * t_x_xz_xx_xy[i];

        t_xy_xz_xx_xx[i] = t_x_xyz_xx_xx[i] - rab_y[i] * t_x_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xy_xx_xx, t_x_xy_xx_xy, t_x_xy_xx_xz, t_x_xy_xx_yy,\
                           t_x_xy_xx_yz, t_x_xy_xx_zz, t_x_xy_xy_xx, t_x_xy_xy_xy, t_x_xy_xy_xz,\
                           t_x_xy_xy_yy, t_x_xy_xy_yz, t_x_xy_xy_zz, t_x_xy_xz_xx, t_x_xy_xz_xy,\
                           t_x_xy_xz_xz, t_x_xy_xz_yy, t_x_xy_xz_yz, t_x_xy_xz_zz, t_x_xy_yy_xx,\
                           t_x_xy_yy_xy, t_x_xy_yy_xz, t_x_xy_yy_yy, t_x_xy_yy_yz, t_x_xy_yy_zz,\
                           t_x_xy_yz_xx, t_x_xy_yz_xy, t_x_xy_yz_xz, t_x_xy_yz_yy, t_x_xy_yz_yz,\
                           t_x_xy_yz_zz, t_x_xy_zz_xx, t_x_xy_zz_xy, t_x_xy_zz_xz, t_x_xy_zz_yy,\
                           t_x_xy_zz_yz, t_x_xy_zz_zz, t_x_xyy_xx_xx, t_x_xyy_xx_xy,\
                           t_x_xyy_xx_xz, t_x_xyy_xx_yy, t_x_xyy_xx_yz, t_x_xyy_xx_zz,\
                           t_x_xyy_xy_xx, t_x_xyy_xy_xy, t_x_xyy_xy_xz, t_x_xyy_xy_yy,\
                           t_x_xyy_xy_yz, t_x_xyy_xy_zz, t_x_xyy_xz_xx, t_x_xyy_xz_xy,\
                           t_x_xyy_xz_xz, t_x_xyy_xz_yy, t_x_xyy_xz_yz, t_x_xyy_xz_zz,\
                           t_x_xyy_yy_xx, t_x_xyy_yy_xy, t_x_xyy_yy_xz, t_x_xyy_yy_yy,\
                           t_x_xyy_yy_yz, t_x_xyy_yy_zz, t_x_xyy_yz_xx, t_x_xyy_yz_xy,\
                           t_x_xyy_yz_xz, t_x_xyy_yz_yy, t_x_xyy_yz_yz, t_x_xyy_yz_zz,\
                           t_x_xyy_zz_xx, t_x_xyy_zz_xy, t_x_xyy_zz_xz, t_x_xyy_zz_yy,\
                           t_x_xyy_zz_yz, t_x_xyy_zz_zz, t_xy_xy_xx_xx, t_xy_xy_xx_xy,\
                           t_xy_xy_xx_xz, t_xy_xy_xx_yy, t_xy_xy_xx_yz, t_xy_xy_xx_zz,\
                           t_xy_xy_xy_xx, t_xy_xy_xy_xy, t_xy_xy_xy_xz, t_xy_xy_xy_yy,\
                           t_xy_xy_xy_yz, t_xy_xy_xy_zz, t_xy_xy_xz_xx, t_xy_xy_xz_xy,\
                           t_xy_xy_xz_xz, t_xy_xy_xz_yy, t_xy_xy_xz_yz, t_xy_xy_xz_zz,\
                           t_xy_xy_yy_xx, t_xy_xy_yy_xy, t_xy_xy_yy_xz, t_xy_xy_yy_yy,\
                           t_xy_xy_yy_yz, t_xy_xy_yy_zz, t_xy_xy_yz_xx, t_xy_xy_yz_xy,\
                           t_xy_xy_yz_xz, t_xy_xy_yz_yy, t_xy_xy_yz_yz, t_xy_xy_yz_zz,\
                           t_xy_xy_zz_xx, t_xy_xy_zz_xy, t_xy_xy_zz_xz, t_xy_xy_zz_yy,\
                           t_xy_xy_zz_yz, t_xy_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_xy_zz_zz[i] = t_x_xyy_zz_zz[i] - rab_y[i] * t_x_xy_zz_zz[i];

        t_xy_xy_zz_yz[i] = t_x_xyy_zz_yz[i] - rab_y[i] * t_x_xy_zz_yz[i];

        t_xy_xy_zz_yy[i] = t_x_xyy_zz_yy[i] - rab_y[i] * t_x_xy_zz_yy[i];

        t_xy_xy_zz_xz[i] = t_x_xyy_zz_xz[i] - rab_y[i] * t_x_xy_zz_xz[i];

        t_xy_xy_zz_xy[i] = t_x_xyy_zz_xy[i] - rab_y[i] * t_x_xy_zz_xy[i];

        t_xy_xy_zz_xx[i] = t_x_xyy_zz_xx[i] - rab_y[i] * t_x_xy_zz_xx[i];

        t_xy_xy_yz_zz[i] = t_x_xyy_yz_zz[i] - rab_y[i] * t_x_xy_yz_zz[i];

        t_xy_xy_yz_yz[i] = t_x_xyy_yz_yz[i] - rab_y[i] * t_x_xy_yz_yz[i];

        t_xy_xy_yz_yy[i] = t_x_xyy_yz_yy[i] - rab_y[i] * t_x_xy_yz_yy[i];

        t_xy_xy_yz_xz[i] = t_x_xyy_yz_xz[i] - rab_y[i] * t_x_xy_yz_xz[i];

        t_xy_xy_yz_xy[i] = t_x_xyy_yz_xy[i] - rab_y[i] * t_x_xy_yz_xy[i];

        t_xy_xy_yz_xx[i] = t_x_xyy_yz_xx[i] - rab_y[i] * t_x_xy_yz_xx[i];

        t_xy_xy_yy_zz[i] = t_x_xyy_yy_zz[i] - rab_y[i] * t_x_xy_yy_zz[i];

        t_xy_xy_yy_yz[i] = t_x_xyy_yy_yz[i] - rab_y[i] * t_x_xy_yy_yz[i];

        t_xy_xy_yy_yy[i] = t_x_xyy_yy_yy[i] - rab_y[i] * t_x_xy_yy_yy[i];

        t_xy_xy_yy_xz[i] = t_x_xyy_yy_xz[i] - rab_y[i] * t_x_xy_yy_xz[i];

        t_xy_xy_yy_xy[i] = t_x_xyy_yy_xy[i] - rab_y[i] * t_x_xy_yy_xy[i];

        t_xy_xy_yy_xx[i] = t_x_xyy_yy_xx[i] - rab_y[i] * t_x_xy_yy_xx[i];

        t_xy_xy_xz_zz[i] = t_x_xyy_xz_zz[i] - rab_y[i] * t_x_xy_xz_zz[i];

        t_xy_xy_xz_yz[i] = t_x_xyy_xz_yz[i] - rab_y[i] * t_x_xy_xz_yz[i];

        t_xy_xy_xz_yy[i] = t_x_xyy_xz_yy[i] - rab_y[i] * t_x_xy_xz_yy[i];

        t_xy_xy_xz_xz[i] = t_x_xyy_xz_xz[i] - rab_y[i] * t_x_xy_xz_xz[i];

        t_xy_xy_xz_xy[i] = t_x_xyy_xz_xy[i] - rab_y[i] * t_x_xy_xz_xy[i];

        t_xy_xy_xz_xx[i] = t_x_xyy_xz_xx[i] - rab_y[i] * t_x_xy_xz_xx[i];

        t_xy_xy_xy_zz[i] = t_x_xyy_xy_zz[i] - rab_y[i] * t_x_xy_xy_zz[i];

        t_xy_xy_xy_yz[i] = t_x_xyy_xy_yz[i] - rab_y[i] * t_x_xy_xy_yz[i];

        t_xy_xy_xy_yy[i] = t_x_xyy_xy_yy[i] - rab_y[i] * t_x_xy_xy_yy[i];

        t_xy_xy_xy_xz[i] = t_x_xyy_xy_xz[i] - rab_y[i] * t_x_xy_xy_xz[i];

        t_xy_xy_xy_xy[i] = t_x_xyy_xy_xy[i] - rab_y[i] * t_x_xy_xy_xy[i];

        t_xy_xy_xy_xx[i] = t_x_xyy_xy_xx[i] - rab_y[i] * t_x_xy_xy_xx[i];

        t_xy_xy_xx_zz[i] = t_x_xyy_xx_zz[i] - rab_y[i] * t_x_xy_xx_zz[i];

        t_xy_xy_xx_yz[i] = t_x_xyy_xx_yz[i] - rab_y[i] * t_x_xy_xx_yz[i];

        t_xy_xy_xx_yy[i] = t_x_xyy_xx_yy[i] - rab_y[i] * t_x_xy_xx_yy[i];

        t_xy_xy_xx_xz[i] = t_x_xyy_xx_xz[i] - rab_y[i] * t_x_xy_xx_xz[i];

        t_xy_xy_xx_xy[i] = t_x_xyy_xx_xy[i] - rab_y[i] * t_x_xy_xx_xy[i];

        t_xy_xy_xx_xx[i] = t_x_xyy_xx_xx[i] - rab_y[i] * t_x_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xx_xx_xx, t_x_xx_xx_xy, t_x_xx_xx_xz, t_x_xx_xx_yy,\
                           t_x_xx_xx_yz, t_x_xx_xx_zz, t_x_xx_xy_xx, t_x_xx_xy_xy, t_x_xx_xy_xz,\
                           t_x_xx_xy_yy, t_x_xx_xy_yz, t_x_xx_xy_zz, t_x_xx_xz_xx, t_x_xx_xz_xy,\
                           t_x_xx_xz_xz, t_x_xx_xz_yy, t_x_xx_xz_yz, t_x_xx_xz_zz, t_x_xx_yy_xx,\
                           t_x_xx_yy_xy, t_x_xx_yy_xz, t_x_xx_yy_yy, t_x_xx_yy_yz, t_x_xx_yy_zz,\
                           t_x_xx_yz_xx, t_x_xx_yz_xy, t_x_xx_yz_xz, t_x_xx_yz_yy, t_x_xx_yz_yz,\
                           t_x_xx_yz_zz, t_x_xx_zz_xx, t_x_xx_zz_xy, t_x_xx_zz_xz, t_x_xx_zz_yy,\
                           t_x_xx_zz_yz, t_x_xx_zz_zz, t_x_xxy_xx_xx, t_x_xxy_xx_xy,\
                           t_x_xxy_xx_xz, t_x_xxy_xx_yy, t_x_xxy_xx_yz, t_x_xxy_xx_zz,\
                           t_x_xxy_xy_xx, t_x_xxy_xy_xy, t_x_xxy_xy_xz, t_x_xxy_xy_yy,\
                           t_x_xxy_xy_yz, t_x_xxy_xy_zz, t_x_xxy_xz_xx, t_x_xxy_xz_xy,\
                           t_x_xxy_xz_xz, t_x_xxy_xz_yy, t_x_xxy_xz_yz, t_x_xxy_xz_zz,\
                           t_x_xxy_yy_xx, t_x_xxy_yy_xy, t_x_xxy_yy_xz, t_x_xxy_yy_yy,\
                           t_x_xxy_yy_yz, t_x_xxy_yy_zz, t_x_xxy_yz_xx, t_x_xxy_yz_xy,\
                           t_x_xxy_yz_xz, t_x_xxy_yz_yy, t_x_xxy_yz_yz, t_x_xxy_yz_zz,\
                           t_x_xxy_zz_xx, t_x_xxy_zz_xy, t_x_xxy_zz_xz, t_x_xxy_zz_yy,\
                           t_x_xxy_zz_yz, t_x_xxy_zz_zz, t_xy_xx_xx_xx, t_xy_xx_xx_xy,\
                           t_xy_xx_xx_xz, t_xy_xx_xx_yy, t_xy_xx_xx_yz, t_xy_xx_xx_zz,\
                           t_xy_xx_xy_xx, t_xy_xx_xy_xy, t_xy_xx_xy_xz, t_xy_xx_xy_yy,\
                           t_xy_xx_xy_yz, t_xy_xx_xy_zz, t_xy_xx_xz_xx, t_xy_xx_xz_xy,\
                           t_xy_xx_xz_xz, t_xy_xx_xz_yy, t_xy_xx_xz_yz, t_xy_xx_xz_zz,\
                           t_xy_xx_yy_xx, t_xy_xx_yy_xy, t_xy_xx_yy_xz, t_xy_xx_yy_yy,\
                           t_xy_xx_yy_yz, t_xy_xx_yy_zz, t_xy_xx_yz_xx, t_xy_xx_yz_xy,\
                           t_xy_xx_yz_xz, t_xy_xx_yz_yy, t_xy_xx_yz_yz, t_xy_xx_yz_zz,\
                           t_xy_xx_zz_xx, t_xy_xx_zz_xy, t_xy_xx_zz_xz, t_xy_xx_zz_yy,\
                           t_xy_xx_zz_yz, t_xy_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_xx_zz_zz[i] = t_x_xxy_zz_zz[i] - rab_y[i] * t_x_xx_zz_zz[i];

        t_xy_xx_zz_yz[i] = t_x_xxy_zz_yz[i] - rab_y[i] * t_x_xx_zz_yz[i];

        t_xy_xx_zz_yy[i] = t_x_xxy_zz_yy[i] - rab_y[i] * t_x_xx_zz_yy[i];

        t_xy_xx_zz_xz[i] = t_x_xxy_zz_xz[i] - rab_y[i] * t_x_xx_zz_xz[i];

        t_xy_xx_zz_xy[i] = t_x_xxy_zz_xy[i] - rab_y[i] * t_x_xx_zz_xy[i];

        t_xy_xx_zz_xx[i] = t_x_xxy_zz_xx[i] - rab_y[i] * t_x_xx_zz_xx[i];

        t_xy_xx_yz_zz[i] = t_x_xxy_yz_zz[i] - rab_y[i] * t_x_xx_yz_zz[i];

        t_xy_xx_yz_yz[i] = t_x_xxy_yz_yz[i] - rab_y[i] * t_x_xx_yz_yz[i];

        t_xy_xx_yz_yy[i] = t_x_xxy_yz_yy[i] - rab_y[i] * t_x_xx_yz_yy[i];

        t_xy_xx_yz_xz[i] = t_x_xxy_yz_xz[i] - rab_y[i] * t_x_xx_yz_xz[i];

        t_xy_xx_yz_xy[i] = t_x_xxy_yz_xy[i] - rab_y[i] * t_x_xx_yz_xy[i];

        t_xy_xx_yz_xx[i] = t_x_xxy_yz_xx[i] - rab_y[i] * t_x_xx_yz_xx[i];

        t_xy_xx_yy_zz[i] = t_x_xxy_yy_zz[i] - rab_y[i] * t_x_xx_yy_zz[i];

        t_xy_xx_yy_yz[i] = t_x_xxy_yy_yz[i] - rab_y[i] * t_x_xx_yy_yz[i];

        t_xy_xx_yy_yy[i] = t_x_xxy_yy_yy[i] - rab_y[i] * t_x_xx_yy_yy[i];

        t_xy_xx_yy_xz[i] = t_x_xxy_yy_xz[i] - rab_y[i] * t_x_xx_yy_xz[i];

        t_xy_xx_yy_xy[i] = t_x_xxy_yy_xy[i] - rab_y[i] * t_x_xx_yy_xy[i];

        t_xy_xx_yy_xx[i] = t_x_xxy_yy_xx[i] - rab_y[i] * t_x_xx_yy_xx[i];

        t_xy_xx_xz_zz[i] = t_x_xxy_xz_zz[i] - rab_y[i] * t_x_xx_xz_zz[i];

        t_xy_xx_xz_yz[i] = t_x_xxy_xz_yz[i] - rab_y[i] * t_x_xx_xz_yz[i];

        t_xy_xx_xz_yy[i] = t_x_xxy_xz_yy[i] - rab_y[i] * t_x_xx_xz_yy[i];

        t_xy_xx_xz_xz[i] = t_x_xxy_xz_xz[i] - rab_y[i] * t_x_xx_xz_xz[i];

        t_xy_xx_xz_xy[i] = t_x_xxy_xz_xy[i] - rab_y[i] * t_x_xx_xz_xy[i];

        t_xy_xx_xz_xx[i] = t_x_xxy_xz_xx[i] - rab_y[i] * t_x_xx_xz_xx[i];

        t_xy_xx_xy_zz[i] = t_x_xxy_xy_zz[i] - rab_y[i] * t_x_xx_xy_zz[i];

        t_xy_xx_xy_yz[i] = t_x_xxy_xy_yz[i] - rab_y[i] * t_x_xx_xy_yz[i];

        t_xy_xx_xy_yy[i] = t_x_xxy_xy_yy[i] - rab_y[i] * t_x_xx_xy_yy[i];

        t_xy_xx_xy_xz[i] = t_x_xxy_xy_xz[i] - rab_y[i] * t_x_xx_xy_xz[i];

        t_xy_xx_xy_xy[i] = t_x_xxy_xy_xy[i] - rab_y[i] * t_x_xx_xy_xy[i];

        t_xy_xx_xy_xx[i] = t_x_xxy_xy_xx[i] - rab_y[i] * t_x_xx_xy_xx[i];

        t_xy_xx_xx_zz[i] = t_x_xxy_xx_zz[i] - rab_y[i] * t_x_xx_xx_zz[i];

        t_xy_xx_xx_yz[i] = t_x_xxy_xx_yz[i] - rab_y[i] * t_x_xx_xx_yz[i];

        t_xy_xx_xx_yy[i] = t_x_xxy_xx_yy[i] - rab_y[i] * t_x_xx_xx_yy[i];

        t_xy_xx_xx_xz[i] = t_x_xxy_xx_xz[i] - rab_y[i] * t_x_xx_xx_xz[i];

        t_xy_xx_xx_xy[i] = t_x_xxy_xx_xy[i] - rab_y[i] * t_x_xx_xx_xy[i];

        t_xy_xx_xx_xx[i] = t_x_xxy_xx_xx[i] - rab_y[i] * t_x_xx_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xzz_xx_xx, t_x_xzz_xx_xy, t_x_xzz_xx_xz, t_x_xzz_xx_yy,\
                           t_x_xzz_xx_yz, t_x_xzz_xx_zz, t_x_xzz_xy_xx, t_x_xzz_xy_xy,\
                           t_x_xzz_xy_xz, t_x_xzz_xy_yy, t_x_xzz_xy_yz, t_x_xzz_xy_zz,\
                           t_x_xzz_xz_xx, t_x_xzz_xz_xy, t_x_xzz_xz_xz, t_x_xzz_xz_yy,\
                           t_x_xzz_xz_yz, t_x_xzz_xz_zz, t_x_xzz_yy_xx, t_x_xzz_yy_xy,\
                           t_x_xzz_yy_xz, t_x_xzz_yy_yy, t_x_xzz_yy_yz, t_x_xzz_yy_zz,\
                           t_x_xzz_yz_xx, t_x_xzz_yz_xy, t_x_xzz_yz_xz, t_x_xzz_yz_yy,\
                           t_x_xzz_yz_yz, t_x_xzz_yz_zz, t_x_xzz_zz_xx, t_x_xzz_zz_xy,\
                           t_x_xzz_zz_xz, t_x_xzz_zz_yy, t_x_xzz_zz_yz, t_x_xzz_zz_zz,\
                           t_x_zz_xx_xx, t_x_zz_xx_xy, t_x_zz_xx_xz, t_x_zz_xx_yy, t_x_zz_xx_yz,\
                           t_x_zz_xx_zz, t_x_zz_xy_xx, t_x_zz_xy_xy, t_x_zz_xy_xz, t_x_zz_xy_yy,\
                           t_x_zz_xy_yz, t_x_zz_xy_zz, t_x_zz_xz_xx, t_x_zz_xz_xy, t_x_zz_xz_xz,\
                           t_x_zz_xz_yy, t_x_zz_xz_yz, t_x_zz_xz_zz, t_x_zz_yy_xx, t_x_zz_yy_xy,\
                           t_x_zz_yy_xz, t_x_zz_yy_yy, t_x_zz_yy_yz, t_x_zz_yy_zz, t_x_zz_yz_xx,\
                           t_x_zz_yz_xy, t_x_zz_yz_xz, t_x_zz_yz_yy, t_x_zz_yz_yz, t_x_zz_yz_zz,\
                           t_x_zz_zz_xx, t_x_zz_zz_xy, t_x_zz_zz_xz, t_x_zz_zz_yy, t_x_zz_zz_yz,\
                           t_x_zz_zz_zz, t_xx_zz_xx_xx, t_xx_zz_xx_xy, t_xx_zz_xx_xz,\
                           t_xx_zz_xx_yy, t_xx_zz_xx_yz, t_xx_zz_xx_zz, t_xx_zz_xy_xx,\
                           t_xx_zz_xy_xy, t_xx_zz_xy_xz, t_xx_zz_xy_yy, t_xx_zz_xy_yz,\
                           t_xx_zz_xy_zz, t_xx_zz_xz_xx, t_xx_zz_xz_xy, t_xx_zz_xz_xz,\
                           t_xx_zz_xz_yy, t_xx_zz_xz_yz, t_xx_zz_xz_zz, t_xx_zz_yy_xx,\
                           t_xx_zz_yy_xy, t_xx_zz_yy_xz, t_xx_zz_yy_yy, t_xx_zz_yy_yz,\
                           t_xx_zz_yy_zz, t_xx_zz_yz_xx, t_xx_zz_yz_xy, t_xx_zz_yz_xz,\
                           t_xx_zz_yz_yy, t_xx_zz_yz_yz, t_xx_zz_yz_zz, t_xx_zz_zz_xx,\
                           t_xx_zz_zz_xy, t_xx_zz_zz_xz, t_xx_zz_zz_yy, t_xx_zz_zz_yz,\
                           t_xx_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_zz_zz_zz[i] = t_x_xzz_zz_zz[i] - rab_x[i] * t_x_zz_zz_zz[i];

        t_xx_zz_zz_yz[i] = t_x_xzz_zz_yz[i] - rab_x[i] * t_x_zz_zz_yz[i];

        t_xx_zz_zz_yy[i] = t_x_xzz_zz_yy[i] - rab_x[i] * t_x_zz_zz_yy[i];

        t_xx_zz_zz_xz[i] = t_x_xzz_zz_xz[i] - rab_x[i] * t_x_zz_zz_xz[i];

        t_xx_zz_zz_xy[i] = t_x_xzz_zz_xy[i] - rab_x[i] * t_x_zz_zz_xy[i];

        t_xx_zz_zz_xx[i] = t_x_xzz_zz_xx[i] - rab_x[i] * t_x_zz_zz_xx[i];

        t_xx_zz_yz_zz[i] = t_x_xzz_yz_zz[i] - rab_x[i] * t_x_zz_yz_zz[i];

        t_xx_zz_yz_yz[i] = t_x_xzz_yz_yz[i] - rab_x[i] * t_x_zz_yz_yz[i];

        t_xx_zz_yz_yy[i] = t_x_xzz_yz_yy[i] - rab_x[i] * t_x_zz_yz_yy[i];

        t_xx_zz_yz_xz[i] = t_x_xzz_yz_xz[i] - rab_x[i] * t_x_zz_yz_xz[i];

        t_xx_zz_yz_xy[i] = t_x_xzz_yz_xy[i] - rab_x[i] * t_x_zz_yz_xy[i];

        t_xx_zz_yz_xx[i] = t_x_xzz_yz_xx[i] - rab_x[i] * t_x_zz_yz_xx[i];

        t_xx_zz_yy_zz[i] = t_x_xzz_yy_zz[i] - rab_x[i] * t_x_zz_yy_zz[i];

        t_xx_zz_yy_yz[i] = t_x_xzz_yy_yz[i] - rab_x[i] * t_x_zz_yy_yz[i];

        t_xx_zz_yy_yy[i] = t_x_xzz_yy_yy[i] - rab_x[i] * t_x_zz_yy_yy[i];

        t_xx_zz_yy_xz[i] = t_x_xzz_yy_xz[i] - rab_x[i] * t_x_zz_yy_xz[i];

        t_xx_zz_yy_xy[i] = t_x_xzz_yy_xy[i] - rab_x[i] * t_x_zz_yy_xy[i];

        t_xx_zz_yy_xx[i] = t_x_xzz_yy_xx[i] - rab_x[i] * t_x_zz_yy_xx[i];

        t_xx_zz_xz_zz[i] = t_x_xzz_xz_zz[i] - rab_x[i] * t_x_zz_xz_zz[i];

        t_xx_zz_xz_yz[i] = t_x_xzz_xz_yz[i] - rab_x[i] * t_x_zz_xz_yz[i];

        t_xx_zz_xz_yy[i] = t_x_xzz_xz_yy[i] - rab_x[i] * t_x_zz_xz_yy[i];

        t_xx_zz_xz_xz[i] = t_x_xzz_xz_xz[i] - rab_x[i] * t_x_zz_xz_xz[i];

        t_xx_zz_xz_xy[i] = t_x_xzz_xz_xy[i] - rab_x[i] * t_x_zz_xz_xy[i];

        t_xx_zz_xz_xx[i] = t_x_xzz_xz_xx[i] - rab_x[i] * t_x_zz_xz_xx[i];

        t_xx_zz_xy_zz[i] = t_x_xzz_xy_zz[i] - rab_x[i] * t_x_zz_xy_zz[i];

        t_xx_zz_xy_yz[i] = t_x_xzz_xy_yz[i] - rab_x[i] * t_x_zz_xy_yz[i];

        t_xx_zz_xy_yy[i] = t_x_xzz_xy_yy[i] - rab_x[i] * t_x_zz_xy_yy[i];

        t_xx_zz_xy_xz[i] = t_x_xzz_xy_xz[i] - rab_x[i] * t_x_zz_xy_xz[i];

        t_xx_zz_xy_xy[i] = t_x_xzz_xy_xy[i] - rab_x[i] * t_x_zz_xy_xy[i];

        t_xx_zz_xy_xx[i] = t_x_xzz_xy_xx[i] - rab_x[i] * t_x_zz_xy_xx[i];

        t_xx_zz_xx_zz[i] = t_x_xzz_xx_zz[i] - rab_x[i] * t_x_zz_xx_zz[i];

        t_xx_zz_xx_yz[i] = t_x_xzz_xx_yz[i] - rab_x[i] * t_x_zz_xx_yz[i];

        t_xx_zz_xx_yy[i] = t_x_xzz_xx_yy[i] - rab_x[i] * t_x_zz_xx_yy[i];

        t_xx_zz_xx_xz[i] = t_x_xzz_xx_xz[i] - rab_x[i] * t_x_zz_xx_xz[i];

        t_xx_zz_xx_xy[i] = t_x_xzz_xx_xy[i] - rab_x[i] * t_x_zz_xx_xy[i];

        t_xx_zz_xx_xx[i] = t_x_xzz_xx_xx[i] - rab_x[i] * t_x_zz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xyz_xx_xx, t_x_xyz_xx_xy, t_x_xyz_xx_xz, t_x_xyz_xx_yy,\
                           t_x_xyz_xx_yz, t_x_xyz_xx_zz, t_x_xyz_xy_xx, t_x_xyz_xy_xy,\
                           t_x_xyz_xy_xz, t_x_xyz_xy_yy, t_x_xyz_xy_yz, t_x_xyz_xy_zz,\
                           t_x_xyz_xz_xx, t_x_xyz_xz_xy, t_x_xyz_xz_xz, t_x_xyz_xz_yy,\
                           t_x_xyz_xz_yz, t_x_xyz_xz_zz, t_x_xyz_yy_xx, t_x_xyz_yy_xy,\
                           t_x_xyz_yy_xz, t_x_xyz_yy_yy, t_x_xyz_yy_yz, t_x_xyz_yy_zz,\
                           t_x_xyz_yz_xx, t_x_xyz_yz_xy, t_x_xyz_yz_xz, t_x_xyz_yz_yy,\
                           t_x_xyz_yz_yz, t_x_xyz_yz_zz, t_x_xyz_zz_xx, t_x_xyz_zz_xy,\
                           t_x_xyz_zz_xz, t_x_xyz_zz_yy, t_x_xyz_zz_yz, t_x_xyz_zz_zz,\
                           t_x_yz_xx_xx, t_x_yz_xx_xy, t_x_yz_xx_xz, t_x_yz_xx_yy, t_x_yz_xx_yz,\
                           t_x_yz_xx_zz, t_x_yz_xy_xx, t_x_yz_xy_xy, t_x_yz_xy_xz, t_x_yz_xy_yy,\
                           t_x_yz_xy_yz, t_x_yz_xy_zz, t_x_yz_xz_xx, t_x_yz_xz_xy, t_x_yz_xz_xz,\
                           t_x_yz_xz_yy, t_x_yz_xz_yz, t_x_yz_xz_zz, t_x_yz_yy_xx, t_x_yz_yy_xy,\
                           t_x_yz_yy_xz, t_x_yz_yy_yy, t_x_yz_yy_yz, t_x_yz_yy_zz, t_x_yz_yz_xx,\
                           t_x_yz_yz_xy, t_x_yz_yz_xz, t_x_yz_yz_yy, t_x_yz_yz_yz, t_x_yz_yz_zz,\
                           t_x_yz_zz_xx, t_x_yz_zz_xy, t_x_yz_zz_xz, t_x_yz_zz_yy, t_x_yz_zz_yz,\
                           t_x_yz_zz_zz, t_xx_yz_xx_xx, t_xx_yz_xx_xy, t_xx_yz_xx_xz,\
                           t_xx_yz_xx_yy, t_xx_yz_xx_yz, t_xx_yz_xx_zz, t_xx_yz_xy_xx,\
                           t_xx_yz_xy_xy, t_xx_yz_xy_xz, t_xx_yz_xy_yy, t_xx_yz_xy_yz,\
                           t_xx_yz_xy_zz, t_xx_yz_xz_xx, t_xx_yz_xz_xy, t_xx_yz_xz_xz,\
                           t_xx_yz_xz_yy, t_xx_yz_xz_yz, t_xx_yz_xz_zz, t_xx_yz_yy_xx,\
                           t_xx_yz_yy_xy, t_xx_yz_yy_xz, t_xx_yz_yy_yy, t_xx_yz_yy_yz,\
                           t_xx_yz_yy_zz, t_xx_yz_yz_xx, t_xx_yz_yz_xy, t_xx_yz_yz_xz,\
                           t_xx_yz_yz_yy, t_xx_yz_yz_yz, t_xx_yz_yz_zz, t_xx_yz_zz_xx,\
                           t_xx_yz_zz_xy, t_xx_yz_zz_xz, t_xx_yz_zz_yy, t_xx_yz_zz_yz,\
                           t_xx_yz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_yz_zz_zz[i] = t_x_xyz_zz_zz[i] - rab_x[i] * t_x_yz_zz_zz[i];

        t_xx_yz_zz_yz[i] = t_x_xyz_zz_yz[i] - rab_x[i] * t_x_yz_zz_yz[i];

        t_xx_yz_zz_yy[i] = t_x_xyz_zz_yy[i] - rab_x[i] * t_x_yz_zz_yy[i];

        t_xx_yz_zz_xz[i] = t_x_xyz_zz_xz[i] - rab_x[i] * t_x_yz_zz_xz[i];

        t_xx_yz_zz_xy[i] = t_x_xyz_zz_xy[i] - rab_x[i] * t_x_yz_zz_xy[i];

        t_xx_yz_zz_xx[i] = t_x_xyz_zz_xx[i] - rab_x[i] * t_x_yz_zz_xx[i];

        t_xx_yz_yz_zz[i] = t_x_xyz_yz_zz[i] - rab_x[i] * t_x_yz_yz_zz[i];

        t_xx_yz_yz_yz[i] = t_x_xyz_yz_yz[i] - rab_x[i] * t_x_yz_yz_yz[i];

        t_xx_yz_yz_yy[i] = t_x_xyz_yz_yy[i] - rab_x[i] * t_x_yz_yz_yy[i];

        t_xx_yz_yz_xz[i] = t_x_xyz_yz_xz[i] - rab_x[i] * t_x_yz_yz_xz[i];

        t_xx_yz_yz_xy[i] = t_x_xyz_yz_xy[i] - rab_x[i] * t_x_yz_yz_xy[i];

        t_xx_yz_yz_xx[i] = t_x_xyz_yz_xx[i] - rab_x[i] * t_x_yz_yz_xx[i];

        t_xx_yz_yy_zz[i] = t_x_xyz_yy_zz[i] - rab_x[i] * t_x_yz_yy_zz[i];

        t_xx_yz_yy_yz[i] = t_x_xyz_yy_yz[i] - rab_x[i] * t_x_yz_yy_yz[i];

        t_xx_yz_yy_yy[i] = t_x_xyz_yy_yy[i] - rab_x[i] * t_x_yz_yy_yy[i];

        t_xx_yz_yy_xz[i] = t_x_xyz_yy_xz[i] - rab_x[i] * t_x_yz_yy_xz[i];

        t_xx_yz_yy_xy[i] = t_x_xyz_yy_xy[i] - rab_x[i] * t_x_yz_yy_xy[i];

        t_xx_yz_yy_xx[i] = t_x_xyz_yy_xx[i] - rab_x[i] * t_x_yz_yy_xx[i];

        t_xx_yz_xz_zz[i] = t_x_xyz_xz_zz[i] - rab_x[i] * t_x_yz_xz_zz[i];

        t_xx_yz_xz_yz[i] = t_x_xyz_xz_yz[i] - rab_x[i] * t_x_yz_xz_yz[i];

        t_xx_yz_xz_yy[i] = t_x_xyz_xz_yy[i] - rab_x[i] * t_x_yz_xz_yy[i];

        t_xx_yz_xz_xz[i] = t_x_xyz_xz_xz[i] - rab_x[i] * t_x_yz_xz_xz[i];

        t_xx_yz_xz_xy[i] = t_x_xyz_xz_xy[i] - rab_x[i] * t_x_yz_xz_xy[i];

        t_xx_yz_xz_xx[i] = t_x_xyz_xz_xx[i] - rab_x[i] * t_x_yz_xz_xx[i];

        t_xx_yz_xy_zz[i] = t_x_xyz_xy_zz[i] - rab_x[i] * t_x_yz_xy_zz[i];

        t_xx_yz_xy_yz[i] = t_x_xyz_xy_yz[i] - rab_x[i] * t_x_yz_xy_yz[i];

        t_xx_yz_xy_yy[i] = t_x_xyz_xy_yy[i] - rab_x[i] * t_x_yz_xy_yy[i];

        t_xx_yz_xy_xz[i] = t_x_xyz_xy_xz[i] - rab_x[i] * t_x_yz_xy_xz[i];

        t_xx_yz_xy_xy[i] = t_x_xyz_xy_xy[i] - rab_x[i] * t_x_yz_xy_xy[i];

        t_xx_yz_xy_xx[i] = t_x_xyz_xy_xx[i] - rab_x[i] * t_x_yz_xy_xx[i];

        t_xx_yz_xx_zz[i] = t_x_xyz_xx_zz[i] - rab_x[i] * t_x_yz_xx_zz[i];

        t_xx_yz_xx_yz[i] = t_x_xyz_xx_yz[i] - rab_x[i] * t_x_yz_xx_yz[i];

        t_xx_yz_xx_yy[i] = t_x_xyz_xx_yy[i] - rab_x[i] * t_x_yz_xx_yy[i];

        t_xx_yz_xx_xz[i] = t_x_xyz_xx_xz[i] - rab_x[i] * t_x_yz_xx_xz[i];

        t_xx_yz_xx_xy[i] = t_x_xyz_xx_xy[i] - rab_x[i] * t_x_yz_xx_xy[i];

        t_xx_yz_xx_xx[i] = t_x_xyz_xx_xx[i] - rab_x[i] * t_x_yz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xyy_xx_xx, t_x_xyy_xx_xy, t_x_xyy_xx_xz, t_x_xyy_xx_yy,\
                           t_x_xyy_xx_yz, t_x_xyy_xx_zz, t_x_xyy_xy_xx, t_x_xyy_xy_xy,\
                           t_x_xyy_xy_xz, t_x_xyy_xy_yy, t_x_xyy_xy_yz, t_x_xyy_xy_zz,\
                           t_x_xyy_xz_xx, t_x_xyy_xz_xy, t_x_xyy_xz_xz, t_x_xyy_xz_yy,\
                           t_x_xyy_xz_yz, t_x_xyy_xz_zz, t_x_xyy_yy_xx, t_x_xyy_yy_xy,\
                           t_x_xyy_yy_xz, t_x_xyy_yy_yy, t_x_xyy_yy_yz, t_x_xyy_yy_zz,\
                           t_x_xyy_yz_xx, t_x_xyy_yz_xy, t_x_xyy_yz_xz, t_x_xyy_yz_yy,\
                           t_x_xyy_yz_yz, t_x_xyy_yz_zz, t_x_xyy_zz_xx, t_x_xyy_zz_xy,\
                           t_x_xyy_zz_xz, t_x_xyy_zz_yy, t_x_xyy_zz_yz, t_x_xyy_zz_zz,\
                           t_x_yy_xx_xx, t_x_yy_xx_xy, t_x_yy_xx_xz, t_x_yy_xx_yy, t_x_yy_xx_yz,\
                           t_x_yy_xx_zz, t_x_yy_xy_xx, t_x_yy_xy_xy, t_x_yy_xy_xz, t_x_yy_xy_yy,\
                           t_x_yy_xy_yz, t_x_yy_xy_zz, t_x_yy_xz_xx, t_x_yy_xz_xy, t_x_yy_xz_xz,\
                           t_x_yy_xz_yy, t_x_yy_xz_yz, t_x_yy_xz_zz, t_x_yy_yy_xx, t_x_yy_yy_xy,\
                           t_x_yy_yy_xz, t_x_yy_yy_yy, t_x_yy_yy_yz, t_x_yy_yy_zz, t_x_yy_yz_xx,\
                           t_x_yy_yz_xy, t_x_yy_yz_xz, t_x_yy_yz_yy, t_x_yy_yz_yz, t_x_yy_yz_zz,\
                           t_x_yy_zz_xx, t_x_yy_zz_xy, t_x_yy_zz_xz, t_x_yy_zz_yy, t_x_yy_zz_yz,\
                           t_x_yy_zz_zz, t_xx_yy_xx_xx, t_xx_yy_xx_xy, t_xx_yy_xx_xz,\
                           t_xx_yy_xx_yy, t_xx_yy_xx_yz, t_xx_yy_xx_zz, t_xx_yy_xy_xx,\
                           t_xx_yy_xy_xy, t_xx_yy_xy_xz, t_xx_yy_xy_yy, t_xx_yy_xy_yz,\
                           t_xx_yy_xy_zz, t_xx_yy_xz_xx, t_xx_yy_xz_xy, t_xx_yy_xz_xz,\
                           t_xx_yy_xz_yy, t_xx_yy_xz_yz, t_xx_yy_xz_zz, t_xx_yy_yy_xx,\
                           t_xx_yy_yy_xy, t_xx_yy_yy_xz, t_xx_yy_yy_yy, t_xx_yy_yy_yz,\
                           t_xx_yy_yy_zz, t_xx_yy_yz_xx, t_xx_yy_yz_xy, t_xx_yy_yz_xz,\
                           t_xx_yy_yz_yy, t_xx_yy_yz_yz, t_xx_yy_yz_zz, t_xx_yy_zz_xx,\
                           t_xx_yy_zz_xy, t_xx_yy_zz_xz, t_xx_yy_zz_yy, t_xx_yy_zz_yz,\
                           t_xx_yy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_yy_zz_zz[i] = t_x_xyy_zz_zz[i] - rab_x[i] * t_x_yy_zz_zz[i];

        t_xx_yy_zz_yz[i] = t_x_xyy_zz_yz[i] - rab_x[i] * t_x_yy_zz_yz[i];

        t_xx_yy_zz_yy[i] = t_x_xyy_zz_yy[i] - rab_x[i] * t_x_yy_zz_yy[i];

        t_xx_yy_zz_xz[i] = t_x_xyy_zz_xz[i] - rab_x[i] * t_x_yy_zz_xz[i];

        t_xx_yy_zz_xy[i] = t_x_xyy_zz_xy[i] - rab_x[i] * t_x_yy_zz_xy[i];

        t_xx_yy_zz_xx[i] = t_x_xyy_zz_xx[i] - rab_x[i] * t_x_yy_zz_xx[i];

        t_xx_yy_yz_zz[i] = t_x_xyy_yz_zz[i] - rab_x[i] * t_x_yy_yz_zz[i];

        t_xx_yy_yz_yz[i] = t_x_xyy_yz_yz[i] - rab_x[i] * t_x_yy_yz_yz[i];

        t_xx_yy_yz_yy[i] = t_x_xyy_yz_yy[i] - rab_x[i] * t_x_yy_yz_yy[i];

        t_xx_yy_yz_xz[i] = t_x_xyy_yz_xz[i] - rab_x[i] * t_x_yy_yz_xz[i];

        t_xx_yy_yz_xy[i] = t_x_xyy_yz_xy[i] - rab_x[i] * t_x_yy_yz_xy[i];

        t_xx_yy_yz_xx[i] = t_x_xyy_yz_xx[i] - rab_x[i] * t_x_yy_yz_xx[i];

        t_xx_yy_yy_zz[i] = t_x_xyy_yy_zz[i] - rab_x[i] * t_x_yy_yy_zz[i];

        t_xx_yy_yy_yz[i] = t_x_xyy_yy_yz[i] - rab_x[i] * t_x_yy_yy_yz[i];

        t_xx_yy_yy_yy[i] = t_x_xyy_yy_yy[i] - rab_x[i] * t_x_yy_yy_yy[i];

        t_xx_yy_yy_xz[i] = t_x_xyy_yy_xz[i] - rab_x[i] * t_x_yy_yy_xz[i];

        t_xx_yy_yy_xy[i] = t_x_xyy_yy_xy[i] - rab_x[i] * t_x_yy_yy_xy[i];

        t_xx_yy_yy_xx[i] = t_x_xyy_yy_xx[i] - rab_x[i] * t_x_yy_yy_xx[i];

        t_xx_yy_xz_zz[i] = t_x_xyy_xz_zz[i] - rab_x[i] * t_x_yy_xz_zz[i];

        t_xx_yy_xz_yz[i] = t_x_xyy_xz_yz[i] - rab_x[i] * t_x_yy_xz_yz[i];

        t_xx_yy_xz_yy[i] = t_x_xyy_xz_yy[i] - rab_x[i] * t_x_yy_xz_yy[i];

        t_xx_yy_xz_xz[i] = t_x_xyy_xz_xz[i] - rab_x[i] * t_x_yy_xz_xz[i];

        t_xx_yy_xz_xy[i] = t_x_xyy_xz_xy[i] - rab_x[i] * t_x_yy_xz_xy[i];

        t_xx_yy_xz_xx[i] = t_x_xyy_xz_xx[i] - rab_x[i] * t_x_yy_xz_xx[i];

        t_xx_yy_xy_zz[i] = t_x_xyy_xy_zz[i] - rab_x[i] * t_x_yy_xy_zz[i];

        t_xx_yy_xy_yz[i] = t_x_xyy_xy_yz[i] - rab_x[i] * t_x_yy_xy_yz[i];

        t_xx_yy_xy_yy[i] = t_x_xyy_xy_yy[i] - rab_x[i] * t_x_yy_xy_yy[i];

        t_xx_yy_xy_xz[i] = t_x_xyy_xy_xz[i] - rab_x[i] * t_x_yy_xy_xz[i];

        t_xx_yy_xy_xy[i] = t_x_xyy_xy_xy[i] - rab_x[i] * t_x_yy_xy_xy[i];

        t_xx_yy_xy_xx[i] = t_x_xyy_xy_xx[i] - rab_x[i] * t_x_yy_xy_xx[i];

        t_xx_yy_xx_zz[i] = t_x_xyy_xx_zz[i] - rab_x[i] * t_x_yy_xx_zz[i];

        t_xx_yy_xx_yz[i] = t_x_xyy_xx_yz[i] - rab_x[i] * t_x_yy_xx_yz[i];

        t_xx_yy_xx_yy[i] = t_x_xyy_xx_yy[i] - rab_x[i] * t_x_yy_xx_yy[i];

        t_xx_yy_xx_xz[i] = t_x_xyy_xx_xz[i] - rab_x[i] * t_x_yy_xx_xz[i];

        t_xx_yy_xx_xy[i] = t_x_xyy_xx_xy[i] - rab_x[i] * t_x_yy_xx_xy[i];

        t_xx_yy_xx_xx[i] = t_x_xyy_xx_xx[i] - rab_x[i] * t_x_yy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xxz_xx_xx, t_x_xxz_xx_xy, t_x_xxz_xx_xz, t_x_xxz_xx_yy,\
                           t_x_xxz_xx_yz, t_x_xxz_xx_zz, t_x_xxz_xy_xx, t_x_xxz_xy_xy,\
                           t_x_xxz_xy_xz, t_x_xxz_xy_yy, t_x_xxz_xy_yz, t_x_xxz_xy_zz,\
                           t_x_xxz_xz_xx, t_x_xxz_xz_xy, t_x_xxz_xz_xz, t_x_xxz_xz_yy,\
                           t_x_xxz_xz_yz, t_x_xxz_xz_zz, t_x_xxz_yy_xx, t_x_xxz_yy_xy,\
                           t_x_xxz_yy_xz, t_x_xxz_yy_yy, t_x_xxz_yy_yz, t_x_xxz_yy_zz,\
                           t_x_xxz_yz_xx, t_x_xxz_yz_xy, t_x_xxz_yz_xz, t_x_xxz_yz_yy,\
                           t_x_xxz_yz_yz, t_x_xxz_yz_zz, t_x_xxz_zz_xx, t_x_xxz_zz_xy,\
                           t_x_xxz_zz_xz, t_x_xxz_zz_yy, t_x_xxz_zz_yz, t_x_xxz_zz_zz,\
                           t_x_xz_xx_xx, t_x_xz_xx_xy, t_x_xz_xx_xz, t_x_xz_xx_yy, t_x_xz_xx_yz,\
                           t_x_xz_xx_zz, t_x_xz_xy_xx, t_x_xz_xy_xy, t_x_xz_xy_xz, t_x_xz_xy_yy,\
                           t_x_xz_xy_yz, t_x_xz_xy_zz, t_x_xz_xz_xx, t_x_xz_xz_xy, t_x_xz_xz_xz,\
                           t_x_xz_xz_yy, t_x_xz_xz_yz, t_x_xz_xz_zz, t_x_xz_yy_xx, t_x_xz_yy_xy,\
                           t_x_xz_yy_xz, t_x_xz_yy_yy, t_x_xz_yy_yz, t_x_xz_yy_zz, t_x_xz_yz_xx,\
                           t_x_xz_yz_xy, t_x_xz_yz_xz, t_x_xz_yz_yy, t_x_xz_yz_yz, t_x_xz_yz_zz,\
                           t_x_xz_zz_xx, t_x_xz_zz_xy, t_x_xz_zz_xz, t_x_xz_zz_yy, t_x_xz_zz_yz,\
                           t_x_xz_zz_zz, t_xx_xz_xx_xx, t_xx_xz_xx_xy, t_xx_xz_xx_xz,\
                           t_xx_xz_xx_yy, t_xx_xz_xx_yz, t_xx_xz_xx_zz, t_xx_xz_xy_xx,\
                           t_xx_xz_xy_xy, t_xx_xz_xy_xz, t_xx_xz_xy_yy, t_xx_xz_xy_yz,\
                           t_xx_xz_xy_zz, t_xx_xz_xz_xx, t_xx_xz_xz_xy, t_xx_xz_xz_xz,\
                           t_xx_xz_xz_yy, t_xx_xz_xz_yz, t_xx_xz_xz_zz, t_xx_xz_yy_xx,\
                           t_xx_xz_yy_xy, t_xx_xz_yy_xz, t_xx_xz_yy_yy, t_xx_xz_yy_yz,\
                           t_xx_xz_yy_zz, t_xx_xz_yz_xx, t_xx_xz_yz_xy, t_xx_xz_yz_xz,\
                           t_xx_xz_yz_yy, t_xx_xz_yz_yz, t_xx_xz_yz_zz, t_xx_xz_zz_xx,\
                           t_xx_xz_zz_xy, t_xx_xz_zz_xz, t_xx_xz_zz_yy, t_xx_xz_zz_yz,\
                           t_xx_xz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_xz_zz_zz[i] = t_x_xxz_zz_zz[i] - rab_x[i] * t_x_xz_zz_zz[i];

        t_xx_xz_zz_yz[i] = t_x_xxz_zz_yz[i] - rab_x[i] * t_x_xz_zz_yz[i];

        t_xx_xz_zz_yy[i] = t_x_xxz_zz_yy[i] - rab_x[i] * t_x_xz_zz_yy[i];

        t_xx_xz_zz_xz[i] = t_x_xxz_zz_xz[i] - rab_x[i] * t_x_xz_zz_xz[i];

        t_xx_xz_zz_xy[i] = t_x_xxz_zz_xy[i] - rab_x[i] * t_x_xz_zz_xy[i];

        t_xx_xz_zz_xx[i] = t_x_xxz_zz_xx[i] - rab_x[i] * t_x_xz_zz_xx[i];

        t_xx_xz_yz_zz[i] = t_x_xxz_yz_zz[i] - rab_x[i] * t_x_xz_yz_zz[i];

        t_xx_xz_yz_yz[i] = t_x_xxz_yz_yz[i] - rab_x[i] * t_x_xz_yz_yz[i];

        t_xx_xz_yz_yy[i] = t_x_xxz_yz_yy[i] - rab_x[i] * t_x_xz_yz_yy[i];

        t_xx_xz_yz_xz[i] = t_x_xxz_yz_xz[i] - rab_x[i] * t_x_xz_yz_xz[i];

        t_xx_xz_yz_xy[i] = t_x_xxz_yz_xy[i] - rab_x[i] * t_x_xz_yz_xy[i];

        t_xx_xz_yz_xx[i] = t_x_xxz_yz_xx[i] - rab_x[i] * t_x_xz_yz_xx[i];

        t_xx_xz_yy_zz[i] = t_x_xxz_yy_zz[i] - rab_x[i] * t_x_xz_yy_zz[i];

        t_xx_xz_yy_yz[i] = t_x_xxz_yy_yz[i] - rab_x[i] * t_x_xz_yy_yz[i];

        t_xx_xz_yy_yy[i] = t_x_xxz_yy_yy[i] - rab_x[i] * t_x_xz_yy_yy[i];

        t_xx_xz_yy_xz[i] = t_x_xxz_yy_xz[i] - rab_x[i] * t_x_xz_yy_xz[i];

        t_xx_xz_yy_xy[i] = t_x_xxz_yy_xy[i] - rab_x[i] * t_x_xz_yy_xy[i];

        t_xx_xz_yy_xx[i] = t_x_xxz_yy_xx[i] - rab_x[i] * t_x_xz_yy_xx[i];

        t_xx_xz_xz_zz[i] = t_x_xxz_xz_zz[i] - rab_x[i] * t_x_xz_xz_zz[i];

        t_xx_xz_xz_yz[i] = t_x_xxz_xz_yz[i] - rab_x[i] * t_x_xz_xz_yz[i];

        t_xx_xz_xz_yy[i] = t_x_xxz_xz_yy[i] - rab_x[i] * t_x_xz_xz_yy[i];

        t_xx_xz_xz_xz[i] = t_x_xxz_xz_xz[i] - rab_x[i] * t_x_xz_xz_xz[i];

        t_xx_xz_xz_xy[i] = t_x_xxz_xz_xy[i] - rab_x[i] * t_x_xz_xz_xy[i];

        t_xx_xz_xz_xx[i] = t_x_xxz_xz_xx[i] - rab_x[i] * t_x_xz_xz_xx[i];

        t_xx_xz_xy_zz[i] = t_x_xxz_xy_zz[i] - rab_x[i] * t_x_xz_xy_zz[i];

        t_xx_xz_xy_yz[i] = t_x_xxz_xy_yz[i] - rab_x[i] * t_x_xz_xy_yz[i];

        t_xx_xz_xy_yy[i] = t_x_xxz_xy_yy[i] - rab_x[i] * t_x_xz_xy_yy[i];

        t_xx_xz_xy_xz[i] = t_x_xxz_xy_xz[i] - rab_x[i] * t_x_xz_xy_xz[i];

        t_xx_xz_xy_xy[i] = t_x_xxz_xy_xy[i] - rab_x[i] * t_x_xz_xy_xy[i];

        t_xx_xz_xy_xx[i] = t_x_xxz_xy_xx[i] - rab_x[i] * t_x_xz_xy_xx[i];

        t_xx_xz_xx_zz[i] = t_x_xxz_xx_zz[i] - rab_x[i] * t_x_xz_xx_zz[i];

        t_xx_xz_xx_yz[i] = t_x_xxz_xx_yz[i] - rab_x[i] * t_x_xz_xx_yz[i];

        t_xx_xz_xx_yy[i] = t_x_xxz_xx_yy[i] - rab_x[i] * t_x_xz_xx_yy[i];

        t_xx_xz_xx_xz[i] = t_x_xxz_xx_xz[i] - rab_x[i] * t_x_xz_xx_xz[i];

        t_xx_xz_xx_xy[i] = t_x_xxz_xx_xy[i] - rab_x[i] * t_x_xz_xx_xy[i];

        t_xx_xz_xx_xx[i] = t_x_xxz_xx_xx[i] - rab_x[i] * t_x_xz_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xxy_xx_xx, t_x_xxy_xx_xy, t_x_xxy_xx_xz, t_x_xxy_xx_yy,\
                           t_x_xxy_xx_yz, t_x_xxy_xx_zz, t_x_xxy_xy_xx, t_x_xxy_xy_xy,\
                           t_x_xxy_xy_xz, t_x_xxy_xy_yy, t_x_xxy_xy_yz, t_x_xxy_xy_zz,\
                           t_x_xxy_xz_xx, t_x_xxy_xz_xy, t_x_xxy_xz_xz, t_x_xxy_xz_yy,\
                           t_x_xxy_xz_yz, t_x_xxy_xz_zz, t_x_xxy_yy_xx, t_x_xxy_yy_xy,\
                           t_x_xxy_yy_xz, t_x_xxy_yy_yy, t_x_xxy_yy_yz, t_x_xxy_yy_zz,\
                           t_x_xxy_yz_xx, t_x_xxy_yz_xy, t_x_xxy_yz_xz, t_x_xxy_yz_yy,\
                           t_x_xxy_yz_yz, t_x_xxy_yz_zz, t_x_xxy_zz_xx, t_x_xxy_zz_xy,\
                           t_x_xxy_zz_xz, t_x_xxy_zz_yy, t_x_xxy_zz_yz, t_x_xxy_zz_zz,\
                           t_x_xy_xx_xx, t_x_xy_xx_xy, t_x_xy_xx_xz, t_x_xy_xx_yy, t_x_xy_xx_yz,\
                           t_x_xy_xx_zz, t_x_xy_xy_xx, t_x_xy_xy_xy, t_x_xy_xy_xz, t_x_xy_xy_yy,\
                           t_x_xy_xy_yz, t_x_xy_xy_zz, t_x_xy_xz_xx, t_x_xy_xz_xy, t_x_xy_xz_xz,\
                           t_x_xy_xz_yy, t_x_xy_xz_yz, t_x_xy_xz_zz, t_x_xy_yy_xx, t_x_xy_yy_xy,\
                           t_x_xy_yy_xz, t_x_xy_yy_yy, t_x_xy_yy_yz, t_x_xy_yy_zz, t_x_xy_yz_xx,\
                           t_x_xy_yz_xy, t_x_xy_yz_xz, t_x_xy_yz_yy, t_x_xy_yz_yz, t_x_xy_yz_zz,\
                           t_x_xy_zz_xx, t_x_xy_zz_xy, t_x_xy_zz_xz, t_x_xy_zz_yy, t_x_xy_zz_yz,\
                           t_x_xy_zz_zz, t_xx_xy_xx_xx, t_xx_xy_xx_xy, t_xx_xy_xx_xz,\
                           t_xx_xy_xx_yy, t_xx_xy_xx_yz, t_xx_xy_xx_zz, t_xx_xy_xy_xx,\
                           t_xx_xy_xy_xy, t_xx_xy_xy_xz, t_xx_xy_xy_yy, t_xx_xy_xy_yz,\
                           t_xx_xy_xy_zz, t_xx_xy_xz_xx, t_xx_xy_xz_xy, t_xx_xy_xz_xz,\
                           t_xx_xy_xz_yy, t_xx_xy_xz_yz, t_xx_xy_xz_zz, t_xx_xy_yy_xx,\
                           t_xx_xy_yy_xy, t_xx_xy_yy_xz, t_xx_xy_yy_yy, t_xx_xy_yy_yz,\
                           t_xx_xy_yy_zz, t_xx_xy_yz_xx, t_xx_xy_yz_xy, t_xx_xy_yz_xz,\
                           t_xx_xy_yz_yy, t_xx_xy_yz_yz, t_xx_xy_yz_zz, t_xx_xy_zz_xx,\
                           t_xx_xy_zz_xy, t_xx_xy_zz_xz, t_xx_xy_zz_yy, t_xx_xy_zz_yz,\
                           t_xx_xy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_xy_zz_zz[i] = t_x_xxy_zz_zz[i] - rab_x[i] * t_x_xy_zz_zz[i];

        t_xx_xy_zz_yz[i] = t_x_xxy_zz_yz[i] - rab_x[i] * t_x_xy_zz_yz[i];

        t_xx_xy_zz_yy[i] = t_x_xxy_zz_yy[i] - rab_x[i] * t_x_xy_zz_yy[i];

        t_xx_xy_zz_xz[i] = t_x_xxy_zz_xz[i] - rab_x[i] * t_x_xy_zz_xz[i];

        t_xx_xy_zz_xy[i] = t_x_xxy_zz_xy[i] - rab_x[i] * t_x_xy_zz_xy[i];

        t_xx_xy_zz_xx[i] = t_x_xxy_zz_xx[i] - rab_x[i] * t_x_xy_zz_xx[i];

        t_xx_xy_yz_zz[i] = t_x_xxy_yz_zz[i] - rab_x[i] * t_x_xy_yz_zz[i];

        t_xx_xy_yz_yz[i] = t_x_xxy_yz_yz[i] - rab_x[i] * t_x_xy_yz_yz[i];

        t_xx_xy_yz_yy[i] = t_x_xxy_yz_yy[i] - rab_x[i] * t_x_xy_yz_yy[i];

        t_xx_xy_yz_xz[i] = t_x_xxy_yz_xz[i] - rab_x[i] * t_x_xy_yz_xz[i];

        t_xx_xy_yz_xy[i] = t_x_xxy_yz_xy[i] - rab_x[i] * t_x_xy_yz_xy[i];

        t_xx_xy_yz_xx[i] = t_x_xxy_yz_xx[i] - rab_x[i] * t_x_xy_yz_xx[i];

        t_xx_xy_yy_zz[i] = t_x_xxy_yy_zz[i] - rab_x[i] * t_x_xy_yy_zz[i];

        t_xx_xy_yy_yz[i] = t_x_xxy_yy_yz[i] - rab_x[i] * t_x_xy_yy_yz[i];

        t_xx_xy_yy_yy[i] = t_x_xxy_yy_yy[i] - rab_x[i] * t_x_xy_yy_yy[i];

        t_xx_xy_yy_xz[i] = t_x_xxy_yy_xz[i] - rab_x[i] * t_x_xy_yy_xz[i];

        t_xx_xy_yy_xy[i] = t_x_xxy_yy_xy[i] - rab_x[i] * t_x_xy_yy_xy[i];

        t_xx_xy_yy_xx[i] = t_x_xxy_yy_xx[i] - rab_x[i] * t_x_xy_yy_xx[i];

        t_xx_xy_xz_zz[i] = t_x_xxy_xz_zz[i] - rab_x[i] * t_x_xy_xz_zz[i];

        t_xx_xy_xz_yz[i] = t_x_xxy_xz_yz[i] - rab_x[i] * t_x_xy_xz_yz[i];

        t_xx_xy_xz_yy[i] = t_x_xxy_xz_yy[i] - rab_x[i] * t_x_xy_xz_yy[i];

        t_xx_xy_xz_xz[i] = t_x_xxy_xz_xz[i] - rab_x[i] * t_x_xy_xz_xz[i];

        t_xx_xy_xz_xy[i] = t_x_xxy_xz_xy[i] - rab_x[i] * t_x_xy_xz_xy[i];

        t_xx_xy_xz_xx[i] = t_x_xxy_xz_xx[i] - rab_x[i] * t_x_xy_xz_xx[i];

        t_xx_xy_xy_zz[i] = t_x_xxy_xy_zz[i] - rab_x[i] * t_x_xy_xy_zz[i];

        t_xx_xy_xy_yz[i] = t_x_xxy_xy_yz[i] - rab_x[i] * t_x_xy_xy_yz[i];

        t_xx_xy_xy_yy[i] = t_x_xxy_xy_yy[i] - rab_x[i] * t_x_xy_xy_yy[i];

        t_xx_xy_xy_xz[i] = t_x_xxy_xy_xz[i] - rab_x[i] * t_x_xy_xy_xz[i];

        t_xx_xy_xy_xy[i] = t_x_xxy_xy_xy[i] - rab_x[i] * t_x_xy_xy_xy[i];

        t_xx_xy_xy_xx[i] = t_x_xxy_xy_xx[i] - rab_x[i] * t_x_xy_xy_xx[i];

        t_xx_xy_xx_zz[i] = t_x_xxy_xx_zz[i] - rab_x[i] * t_x_xy_xx_zz[i];

        t_xx_xy_xx_yz[i] = t_x_xxy_xx_yz[i] - rab_x[i] * t_x_xy_xx_yz[i];

        t_xx_xy_xx_yy[i] = t_x_xxy_xx_yy[i] - rab_x[i] * t_x_xy_xx_yy[i];

        t_xx_xy_xx_xz[i] = t_x_xxy_xx_xz[i] - rab_x[i] * t_x_xy_xx_xz[i];

        t_xx_xy_xx_xy[i] = t_x_xxy_xx_xy[i] - rab_x[i] * t_x_xy_xx_xy[i];

        t_xx_xy_xx_xx[i] = t_x_xxy_xx_xx[i] - rab_x[i] * t_x_xy_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xx_xx_xx, t_x_xx_xx_xy, t_x_xx_xx_xz, t_x_xx_xx_yy,\
                           t_x_xx_xx_yz, t_x_xx_xx_zz, t_x_xx_xy_xx, t_x_xx_xy_xy, t_x_xx_xy_xz,\
                           t_x_xx_xy_yy, t_x_xx_xy_yz, t_x_xx_xy_zz, t_x_xx_xz_xx, t_x_xx_xz_xy,\
                           t_x_xx_xz_xz, t_x_xx_xz_yy, t_x_xx_xz_yz, t_x_xx_xz_zz, t_x_xx_yy_xx,\
                           t_x_xx_yy_xy, t_x_xx_yy_xz, t_x_xx_yy_yy, t_x_xx_yy_yz, t_x_xx_yy_zz,\
                           t_x_xx_yz_xx, t_x_xx_yz_xy, t_x_xx_yz_xz, t_x_xx_yz_yy, t_x_xx_yz_yz,\
                           t_x_xx_yz_zz, t_x_xx_zz_xx, t_x_xx_zz_xy, t_x_xx_zz_xz, t_x_xx_zz_yy,\
                           t_x_xx_zz_yz, t_x_xx_zz_zz, t_x_xxx_xx_xx, t_x_xxx_xx_xy,\
                           t_x_xxx_xx_xz, t_x_xxx_xx_yy, t_x_xxx_xx_yz, t_x_xxx_xx_zz,\
                           t_x_xxx_xy_xx, t_x_xxx_xy_xy, t_x_xxx_xy_xz, t_x_xxx_xy_yy,\
                           t_x_xxx_xy_yz, t_x_xxx_xy_zz, t_x_xxx_xz_xx, t_x_xxx_xz_xy,\
                           t_x_xxx_xz_xz, t_x_xxx_xz_yy, t_x_xxx_xz_yz, t_x_xxx_xz_zz,\
                           t_x_xxx_yy_xx, t_x_xxx_yy_xy, t_x_xxx_yy_xz, t_x_xxx_yy_yy,\
                           t_x_xxx_yy_yz, t_x_xxx_yy_zz, t_x_xxx_yz_xx, t_x_xxx_yz_xy,\
                           t_x_xxx_yz_xz, t_x_xxx_yz_yy, t_x_xxx_yz_yz, t_x_xxx_yz_zz,\
                           t_x_xxx_zz_xx, t_x_xxx_zz_xy, t_x_xxx_zz_xz, t_x_xxx_zz_yy,\
                           t_x_xxx_zz_yz, t_x_xxx_zz_zz, t_xx_xx_xx_xx, t_xx_xx_xx_xy,\
                           t_xx_xx_xx_xz, t_xx_xx_xx_yy, t_xx_xx_xx_yz, t_xx_xx_xx_zz,\
                           t_xx_xx_xy_xx, t_xx_xx_xy_xy, t_xx_xx_xy_xz, t_xx_xx_xy_yy,\
                           t_xx_xx_xy_yz, t_xx_xx_xy_zz, t_xx_xx_xz_xx, t_xx_xx_xz_xy,\
                           t_xx_xx_xz_xz, t_xx_xx_xz_yy, t_xx_xx_xz_yz, t_xx_xx_xz_zz,\
                           t_xx_xx_yy_xx, t_xx_xx_yy_xy, t_xx_xx_yy_xz, t_xx_xx_yy_yy,\
                           t_xx_xx_yy_yz, t_xx_xx_yy_zz, t_xx_xx_yz_xx, t_xx_xx_yz_xy,\
                           t_xx_xx_yz_xz, t_xx_xx_yz_yy, t_xx_xx_yz_yz, t_xx_xx_yz_zz,\
                           t_xx_xx_zz_xx, t_xx_xx_zz_xy, t_xx_xx_zz_xz, t_xx_xx_zz_yy,\
                           t_xx_xx_zz_yz, t_xx_xx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_xx_zz_zz[i] = t_x_xxx_zz_zz[i] - rab_x[i] * t_x_xx_zz_zz[i];

        t_xx_xx_zz_yz[i] = t_x_xxx_zz_yz[i] - rab_x[i] * t_x_xx_zz_yz[i];

        t_xx_xx_zz_yy[i] = t_x_xxx_zz_yy[i] - rab_x[i] * t_x_xx_zz_yy[i];

        t_xx_xx_zz_xz[i] = t_x_xxx_zz_xz[i] - rab_x[i] * t_x_xx_zz_xz[i];

        t_xx_xx_zz_xy[i] = t_x_xxx_zz_xy[i] - rab_x[i] * t_x_xx_zz_xy[i];

        t_xx_xx_zz_xx[i] = t_x_xxx_zz_xx[i] - rab_x[i] * t_x_xx_zz_xx[i];

        t_xx_xx_yz_zz[i] = t_x_xxx_yz_zz[i] - rab_x[i] * t_x_xx_yz_zz[i];

        t_xx_xx_yz_yz[i] = t_x_xxx_yz_yz[i] - rab_x[i] * t_x_xx_yz_yz[i];

        t_xx_xx_yz_yy[i] = t_x_xxx_yz_yy[i] - rab_x[i] * t_x_xx_yz_yy[i];

        t_xx_xx_yz_xz[i] = t_x_xxx_yz_xz[i] - rab_x[i] * t_x_xx_yz_xz[i];

        t_xx_xx_yz_xy[i] = t_x_xxx_yz_xy[i] - rab_x[i] * t_x_xx_yz_xy[i];

        t_xx_xx_yz_xx[i] = t_x_xxx_yz_xx[i] - rab_x[i] * t_x_xx_yz_xx[i];

        t_xx_xx_yy_zz[i] = t_x_xxx_yy_zz[i] - rab_x[i] * t_x_xx_yy_zz[i];

        t_xx_xx_yy_yz[i] = t_x_xxx_yy_yz[i] - rab_x[i] * t_x_xx_yy_yz[i];

        t_xx_xx_yy_yy[i] = t_x_xxx_yy_yy[i] - rab_x[i] * t_x_xx_yy_yy[i];

        t_xx_xx_yy_xz[i] = t_x_xxx_yy_xz[i] - rab_x[i] * t_x_xx_yy_xz[i];

        t_xx_xx_yy_xy[i] = t_x_xxx_yy_xy[i] - rab_x[i] * t_x_xx_yy_xy[i];

        t_xx_xx_yy_xx[i] = t_x_xxx_yy_xx[i] - rab_x[i] * t_x_xx_yy_xx[i];

        t_xx_xx_xz_zz[i] = t_x_xxx_xz_zz[i] - rab_x[i] * t_x_xx_xz_zz[i];

        t_xx_xx_xz_yz[i] = t_x_xxx_xz_yz[i] - rab_x[i] * t_x_xx_xz_yz[i];

        t_xx_xx_xz_yy[i] = t_x_xxx_xz_yy[i] - rab_x[i] * t_x_xx_xz_yy[i];

        t_xx_xx_xz_xz[i] = t_x_xxx_xz_xz[i] - rab_x[i] * t_x_xx_xz_xz[i];

        t_xx_xx_xz_xy[i] = t_x_xxx_xz_xy[i] - rab_x[i] * t_x_xx_xz_xy[i];

        t_xx_xx_xz_xx[i] = t_x_xxx_xz_xx[i] - rab_x[i] * t_x_xx_xz_xx[i];

        t_xx_xx_xy_zz[i] = t_x_xxx_xy_zz[i] - rab_x[i] * t_x_xx_xy_zz[i];

        t_xx_xx_xy_yz[i] = t_x_xxx_xy_yz[i] - rab_x[i] * t_x_xx_xy_yz[i];

        t_xx_xx_xy_yy[i] = t_x_xxx_xy_yy[i] - rab_x[i] * t_x_xx_xy_yy[i];

        t_xx_xx_xy_xz[i] = t_x_xxx_xy_xz[i] - rab_x[i] * t_x_xx_xy_xz[i];

        t_xx_xx_xy_xy[i] = t_x_xxx_xy_xy[i] - rab_x[i] * t_x_xx_xy_xy[i];

        t_xx_xx_xy_xx[i] = t_x_xxx_xy_xx[i] - rab_x[i] * t_x_xx_xy_xx[i];

        t_xx_xx_xx_zz[i] = t_x_xxx_xx_zz[i] - rab_x[i] * t_x_xx_xx_zz[i];

        t_xx_xx_xx_yz[i] = t_x_xxx_xx_yz[i] - rab_x[i] * t_x_xx_xx_yz[i];

        t_xx_xx_xx_yy[i] = t_x_xxx_xx_yy[i] - rab_x[i] * t_x_xx_xx_yy[i];

        t_xx_xx_xx_xz[i] = t_x_xxx_xx_xz[i] - rab_x[i] * t_x_xx_xx_xz[i];

        t_xx_xx_xx_xy[i] = t_x_xxx_xx_xy[i] - rab_x[i] * t_x_xx_xx_xy[i];

        t_xx_xx_xx_xx[i] = t_x_xxx_xx_xx[i] - rab_x[i] * t_x_xx_xx_xx[i];
    }
}


} // derirec namespace
