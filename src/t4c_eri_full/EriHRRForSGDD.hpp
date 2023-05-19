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
compHostHRRForSGDD_V0(      BufferHostXY<T>&      intsBufferSGDD,
                      const BufferHostX<int32_t>& intsIndexesSGDD,
                      const BufferHostXY<T>&      intsBufferSGPD,
                      const BufferHostX<int32_t>& intsIndexesSGPD,
                      const BufferHostXY<T>&      intsBufferSGPF,
                      const BufferHostX<int32_t>& intsIndexesSGPF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SGPF) integral components

    t_0_zzzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(0));

    t_0_zzzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(1));

    t_0_zzzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(2));

    t_0_zzzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(3));

    t_0_zzzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(4));

    t_0_zzzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(5));

    t_0_zzzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(6));

    t_0_zzzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(7));

    t_0_zzzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(8));

    t_0_zzzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(9));

    t_0_zzzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(10));

    t_0_zzzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(11));

    t_0_zzzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(12));

    t_0_zzzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(13));

    t_0_zzzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(14));

    t_0_zzzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(15));

    t_0_zzzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(16));

    t_0_zzzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(17));

    t_0_zzzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(18));

    t_0_zzzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(19));

    t_0_zzzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(20));

    t_0_zzzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(21));

    t_0_zzzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(22));

    t_0_zzzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(23));

    t_0_zzzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(24));

    t_0_yzzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(25));

    t_0_yzzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(26));

    t_0_yzzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(27));

    t_0_yzzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(28));

    t_0_yzzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(29));

    t_0_yzzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(30));

    t_0_yzzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(31));

    t_0_yzzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(32));

    t_0_yzzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(33));

    t_0_yzzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(34));

    t_0_yzzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(35));

    t_0_yzzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(36));

    t_0_yzzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(37));

    t_0_yzzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(38));

    t_0_yzzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(39));

    t_0_yzzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(40));

    t_0_yzzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(41));

    t_0_yzzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(42));

    t_0_yzzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(43));

    t_0_yzzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(44));

    t_0_yzzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(45));

    t_0_yzzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(46));

    t_0_yzzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(47));

    t_0_yzzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(48));

    t_0_yzzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(49));

    t_0_yyzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(50));

    t_0_yyzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(51));

    t_0_yyzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(52));

    t_0_yyzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(53));

    t_0_yyzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(54));

    t_0_yyzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(55));

    t_0_yyzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(56));

    t_0_yyzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(57));

    t_0_yyzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(58));

    t_0_yyzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(59));

    t_0_yyzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(60));

    t_0_yyzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(61));

    t_0_yyzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(62));

    t_0_yyzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(63));

    t_0_yyzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(64));

    t_0_yyzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(65));

    t_0_yyzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(66));

    t_0_yyzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(67));

    t_0_yyzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(68));

    t_0_yyzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(69));

    t_0_yyzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(70));

    t_0_yyzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(71));

    t_0_yyzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(72));

    t_0_yyzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(73));

    t_0_yyzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(74));

    t_0_yyyz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(75));

    t_0_yyyz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(76));

    t_0_yyyz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(77));

    t_0_yyyz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(78));

    t_0_yyyz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(79));

    t_0_yyyz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(80));

    t_0_yyyz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(81));

    t_0_yyyz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(82));

    t_0_yyyz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(83));

    t_0_yyyz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(84));

    t_0_yyyz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(85));

    t_0_yyyz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(86));

    t_0_yyyz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(87));

    t_0_yyyz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(88));

    t_0_yyyz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(89));

    t_0_yyyz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(90));

    t_0_yyyz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(91));

    t_0_yyyz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(92));

    t_0_yyyz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(93));

    t_0_yyyz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(94));

    t_0_yyyz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(95));

    t_0_yyyz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(96));

    t_0_yyyz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(97));

    t_0_yyyz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(98));

    t_0_yyyz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(99));

    t_0_yyyy_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(100));

    t_0_yyyy_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(101));

    t_0_yyyy_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(102));

    t_0_yyyy_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(103));

    t_0_yyyy_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(104));

    t_0_yyyy_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(105));

    t_0_yyyy_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(106));

    t_0_yyyy_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(107));

    t_0_yyyy_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(108));

    t_0_yyyy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(109));

    t_0_yyyy_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(110));

    t_0_yyyy_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(111));

    t_0_yyyy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(112));

    t_0_yyyy_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(113));

    t_0_yyyy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(114));

    t_0_yyyy_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(115));

    t_0_yyyy_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(116));

    t_0_yyyy_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(117));

    t_0_yyyy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(118));

    t_0_yyyy_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(119));

    t_0_yyyy_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(120));

    t_0_yyyy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(121));

    t_0_yyyy_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(122));

    t_0_yyyy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(123));

    t_0_yyyy_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(124));

    t_0_xzzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(125));

    t_0_xzzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(126));

    t_0_xzzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(127));

    t_0_xzzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(128));

    t_0_xzzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(129));

    t_0_xzzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(130));

    t_0_xzzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(131));

    t_0_xzzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(132));

    t_0_xzzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(133));

    t_0_xzzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(134));

    t_0_xzzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(135));

    t_0_xzzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(136));

    t_0_xzzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(137));

    t_0_xzzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(138));

    t_0_xzzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(139));

    t_0_xzzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(140));

    t_0_xzzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(141));

    t_0_xzzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(142));

    t_0_xzzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(143));

    t_0_xzzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(144));

    t_0_xzzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(145));

    t_0_xzzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(146));

    t_0_xzzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(147));

    t_0_xzzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(148));

    t_0_xzzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(149));

    t_0_xyzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(150));

    t_0_xyzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(151));

    t_0_xyzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(152));

    t_0_xyzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(153));

    t_0_xyzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(154));

    t_0_xyzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(155));

    t_0_xyzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(156));

    t_0_xyzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(157));

    t_0_xyzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(158));

    t_0_xyzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(159));

    t_0_xyzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(160));

    t_0_xyzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(161));

    t_0_xyzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(162));

    t_0_xyzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(163));

    t_0_xyzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(164));

    t_0_xyzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(165));

    t_0_xyzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(166));

    t_0_xyzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(167));

    t_0_xyzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(168));

    t_0_xyzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(169));

    t_0_xyzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(170));

    t_0_xyzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(171));

    t_0_xyzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(172));

    t_0_xyzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(173));

    t_0_xyzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(174));

    t_0_xyyz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(175));

    t_0_xyyz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(176));

    t_0_xyyz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(177));

    t_0_xyyz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(178));

    t_0_xyyz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(179));

    t_0_xyyz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(180));

    t_0_xyyz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(181));

    t_0_xyyz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(182));

    t_0_xyyz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(183));

    t_0_xyyz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(184));

    t_0_xyyz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(185));

    t_0_xyyz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(186));

    t_0_xyyz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(187));

    t_0_xyyz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(188));

    t_0_xyyz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(189));

    t_0_xyyz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(190));

    t_0_xyyz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(191));

    t_0_xyyz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(192));

    t_0_xyyz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(193));

    t_0_xyyz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(194));

    t_0_xyyz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(195));

    t_0_xyyz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(196));

    t_0_xyyz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(197));

    t_0_xyyz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(198));

    t_0_xyyz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(199));

    t_0_xyyy_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(200));

    t_0_xyyy_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(201));

    t_0_xyyy_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(202));

    t_0_xyyy_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(203));

    t_0_xyyy_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(204));

    t_0_xyyy_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(205));

    t_0_xyyy_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(206));

    t_0_xyyy_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(207));

    t_0_xyyy_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(208));

    t_0_xyyy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(209));

    t_0_xyyy_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(210));

    t_0_xyyy_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(211));

    t_0_xyyy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(212));

    t_0_xyyy_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(213));

    t_0_xyyy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(214));

    t_0_xyyy_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(215));

    t_0_xyyy_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(216));

    t_0_xyyy_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(217));

    t_0_xyyy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(218));

    t_0_xyyy_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(219));

    t_0_xyyy_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(220));

    t_0_xyyy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(221));

    t_0_xyyy_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(222));

    t_0_xyyy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(223));

    t_0_xyyy_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(224));

    t_0_xxzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(225));

    t_0_xxzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(226));

    t_0_xxzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(227));

    t_0_xxzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(228));

    t_0_xxzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(229));

    t_0_xxzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(230));

    t_0_xxzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(231));

    t_0_xxzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(232));

    t_0_xxzz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(233));

    t_0_xxzz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(234));

    t_0_xxzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(235));

    t_0_xxzz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(236));

    t_0_xxzz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(237));

    t_0_xxzz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(238));

    t_0_xxzz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(239));

    t_0_xxzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(240));

    t_0_xxzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(241));

    t_0_xxzz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(242));

    t_0_xxzz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(243));

    t_0_xxzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(244));

    t_0_xxzz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(245));

    t_0_xxzz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(246));

    t_0_xxzz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(247));

    t_0_xxzz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(248));

    t_0_xxzz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(249));

    t_0_xxyz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(250));

    t_0_xxyz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(251));

    t_0_xxyz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(252));

    t_0_xxyz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(253));

    t_0_xxyz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(254));

    t_0_xxyz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(255));

    t_0_xxyz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(256));

    t_0_xxyz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(257));

    t_0_xxyz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(258));

    t_0_xxyz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(259));

    t_0_xxyz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(260));

    t_0_xxyz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(261));

    t_0_xxyz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(262));

    t_0_xxyz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(263));

    t_0_xxyz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(264));

    t_0_xxyz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(265));

    t_0_xxyz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(266));

    t_0_xxyz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(267));

    t_0_xxyz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(268));

    t_0_xxyz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(269));

    t_0_xxyz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(270));

    t_0_xxyz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(271));

    t_0_xxyz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(272));

    t_0_xxyz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(273));

    t_0_xxyz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(274));

    t_0_xxyy_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(275));

    t_0_xxyy_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(276));

    t_0_xxyy_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(277));

    t_0_xxyy_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(278));

    t_0_xxyy_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(279));

    t_0_xxyy_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(280));

    t_0_xxyy_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(281));

    t_0_xxyy_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(282));

    t_0_xxyy_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(283));

    t_0_xxyy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(284));

    t_0_xxyy_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(285));

    t_0_xxyy_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(286));

    t_0_xxyy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(287));

    t_0_xxyy_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(288));

    t_0_xxyy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(289));

    t_0_xxyy_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(290));

    t_0_xxyy_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(291));

    t_0_xxyy_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(292));

    t_0_xxyy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(293));

    t_0_xxyy_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(294));

    t_0_xxyy_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(295));

    t_0_xxyy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(296));

    t_0_xxyy_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(297));

    t_0_xxyy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(298));

    t_0_xxyy_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(299));

    t_0_xxxz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(300));

    t_0_xxxz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(301));

    t_0_xxxz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(302));

    t_0_xxxz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(303));

    t_0_xxxz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(304));

    t_0_xxxz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(305));

    t_0_xxxz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(306));

    t_0_xxxz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(307));

    t_0_xxxz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(308));

    t_0_xxxz_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(309));

    t_0_xxxz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(310));

    t_0_xxxz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(311));

    t_0_xxxz_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(312));

    t_0_xxxz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(313));

    t_0_xxxz_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(314));

    t_0_xxxz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(315));

    t_0_xxxz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(316));

    t_0_xxxz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(317));

    t_0_xxxz_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(318));

    t_0_xxxz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(319));

    t_0_xxxz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(320));

    t_0_xxxz_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(321));

    t_0_xxxz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(322));

    t_0_xxxz_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(323));

    t_0_xxxz_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(324));

    t_0_xxxy_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(325));

    t_0_xxxy_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(326));

    t_0_xxxy_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(327));

    t_0_xxxy_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(328));

    t_0_xxxy_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(329));

    t_0_xxxy_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(330));

    t_0_xxxy_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(331));

    t_0_xxxy_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(332));

    t_0_xxxy_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(333));

    t_0_xxxy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(334));

    t_0_xxxy_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(335));

    t_0_xxxy_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(336));

    t_0_xxxy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(337));

    t_0_xxxy_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(338));

    t_0_xxxy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(339));

    t_0_xxxy_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(340));

    t_0_xxxy_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(341));

    t_0_xxxy_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(342));

    t_0_xxxy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(343));

    t_0_xxxy_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(344));

    t_0_xxxy_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(345));

    t_0_xxxy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(346));

    t_0_xxxy_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(347));

    t_0_xxxy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(348));

    t_0_xxxy_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(349));

    t_0_xxxx_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(350));

    t_0_xxxx_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(351));

    t_0_xxxx_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(352));

    t_0_xxxx_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(353));

    t_0_xxxx_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(354));

    t_0_xxxx_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(355));

    t_0_xxxx_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(356));

    t_0_xxxx_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(357));

    t_0_xxxx_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(358));

    t_0_xxxx_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(359));

    t_0_xxxx_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(360));

    t_0_xxxx_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(361));

    t_0_xxxx_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(362));

    t_0_xxxx_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(363));

    t_0_xxxx_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(364));

    t_0_xxxx_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(365));

    t_0_xxxx_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(366));

    t_0_xxxx_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(367));

    t_0_xxxx_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(368));

    t_0_xxxx_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(369));

    t_0_xxxx_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(370));

    t_0_xxxx_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(371));

    t_0_xxxx_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(372));

    t_0_xxxx_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(373));

    t_0_xxxx_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(374));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_zzzz_x_xx, t_0_zzzz_x_xxx, t_0_zzzz_x_xxy,\
                           t_0_zzzz_x_xxz, t_0_zzzz_x_xy, t_0_zzzz_x_xyy, t_0_zzzz_x_xyz,\
                           t_0_zzzz_x_xz, t_0_zzzz_x_xzz, t_0_zzzz_x_yy, t_0_zzzz_x_yyy,\
                           t_0_zzzz_x_yyz, t_0_zzzz_x_yz, t_0_zzzz_x_yzz, t_0_zzzz_x_zz,\
                           t_0_zzzz_x_zzz, t_0_zzzz_xx_xx, t_0_zzzz_xx_xy, t_0_zzzz_xx_xz,\
                           t_0_zzzz_xx_yy, t_0_zzzz_xx_yz, t_0_zzzz_xx_zz, t_0_zzzz_xy_xx,\
                           t_0_zzzz_xy_xy, t_0_zzzz_xy_xz, t_0_zzzz_xy_yy, t_0_zzzz_xy_yz,\
                           t_0_zzzz_xy_zz, t_0_zzzz_xz_xx, t_0_zzzz_xz_xy, t_0_zzzz_xz_xz,\
                           t_0_zzzz_xz_yy, t_0_zzzz_xz_yz, t_0_zzzz_xz_zz, t_0_zzzz_y_xx,\
                           t_0_zzzz_y_xxy, t_0_zzzz_y_xxz, t_0_zzzz_y_xy, t_0_zzzz_y_xyy,\
                           t_0_zzzz_y_xyz, t_0_zzzz_y_xz, t_0_zzzz_y_xzz, t_0_zzzz_y_yy,\
                           t_0_zzzz_y_yyy, t_0_zzzz_y_yyz, t_0_zzzz_y_yz, t_0_zzzz_y_yzz,\
                           t_0_zzzz_y_zz, t_0_zzzz_y_zzz, t_0_zzzz_yy_xx, t_0_zzzz_yy_xy,\
                           t_0_zzzz_yy_xz, t_0_zzzz_yy_yy, t_0_zzzz_yy_yz, t_0_zzzz_yy_zz,\
                           t_0_zzzz_yz_xx, t_0_zzzz_yz_xy, t_0_zzzz_yz_xz, t_0_zzzz_yz_yy,\
                           t_0_zzzz_yz_yz, t_0_zzzz_yz_zz, t_0_zzzz_z_xx, t_0_zzzz_z_xxz,\
                           t_0_zzzz_z_xy, t_0_zzzz_z_xyz, t_0_zzzz_z_xz, t_0_zzzz_z_xzz,\
                           t_0_zzzz_z_yy, t_0_zzzz_z_yyz, t_0_zzzz_z_yz, t_0_zzzz_z_yzz,\
                           t_0_zzzz_z_zz, t_0_zzzz_z_zzz, t_0_zzzz_zz_xx, t_0_zzzz_zz_xy,\
                           t_0_zzzz_zz_xz, t_0_zzzz_zz_yy, t_0_zzzz_zz_yz, t_0_zzzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_zz_zz[i] = t_0_zzzz_z_zzz[i] - rcd_z[i] * t_0_zzzz_z_zz[i];

        t_0_zzzz_zz_yz[i] = t_0_zzzz_z_yzz[i] - rcd_z[i] * t_0_zzzz_z_yz[i];

        t_0_zzzz_zz_yy[i] = t_0_zzzz_z_yyz[i] - rcd_z[i] * t_0_zzzz_z_yy[i];

        t_0_zzzz_zz_xz[i] = t_0_zzzz_z_xzz[i] - rcd_z[i] * t_0_zzzz_z_xz[i];

        t_0_zzzz_zz_xy[i] = t_0_zzzz_z_xyz[i] - rcd_z[i] * t_0_zzzz_z_xy[i];

        t_0_zzzz_zz_xx[i] = t_0_zzzz_z_xxz[i] - rcd_z[i] * t_0_zzzz_z_xx[i];

        t_0_zzzz_yz_zz[i] = t_0_zzzz_y_zzz[i] - rcd_z[i] * t_0_zzzz_y_zz[i];

        t_0_zzzz_yz_yz[i] = t_0_zzzz_y_yzz[i] - rcd_z[i] * t_0_zzzz_y_yz[i];

        t_0_zzzz_yz_yy[i] = t_0_zzzz_y_yyz[i] - rcd_z[i] * t_0_zzzz_y_yy[i];

        t_0_zzzz_yz_xz[i] = t_0_zzzz_y_xzz[i] - rcd_z[i] * t_0_zzzz_y_xz[i];

        t_0_zzzz_yz_xy[i] = t_0_zzzz_y_xyz[i] - rcd_z[i] * t_0_zzzz_y_xy[i];

        t_0_zzzz_yz_xx[i] = t_0_zzzz_y_xxz[i] - rcd_z[i] * t_0_zzzz_y_xx[i];

        t_0_zzzz_yy_zz[i] = t_0_zzzz_y_yzz[i] - rcd_y[i] * t_0_zzzz_y_zz[i];

        t_0_zzzz_yy_yz[i] = t_0_zzzz_y_yyz[i] - rcd_y[i] * t_0_zzzz_y_yz[i];

        t_0_zzzz_yy_yy[i] = t_0_zzzz_y_yyy[i] - rcd_y[i] * t_0_zzzz_y_yy[i];

        t_0_zzzz_yy_xz[i] = t_0_zzzz_y_xyz[i] - rcd_y[i] * t_0_zzzz_y_xz[i];

        t_0_zzzz_yy_xy[i] = t_0_zzzz_y_xyy[i] - rcd_y[i] * t_0_zzzz_y_xy[i];

        t_0_zzzz_yy_xx[i] = t_0_zzzz_y_xxy[i] - rcd_y[i] * t_0_zzzz_y_xx[i];

        t_0_zzzz_xz_zz[i] = t_0_zzzz_x_zzz[i] - rcd_z[i] * t_0_zzzz_x_zz[i];

        t_0_zzzz_xz_yz[i] = t_0_zzzz_x_yzz[i] - rcd_z[i] * t_0_zzzz_x_yz[i];

        t_0_zzzz_xz_yy[i] = t_0_zzzz_x_yyz[i] - rcd_z[i] * t_0_zzzz_x_yy[i];

        t_0_zzzz_xz_xz[i] = t_0_zzzz_x_xzz[i] - rcd_z[i] * t_0_zzzz_x_xz[i];

        t_0_zzzz_xz_xy[i] = t_0_zzzz_x_xyz[i] - rcd_z[i] * t_0_zzzz_x_xy[i];

        t_0_zzzz_xz_xx[i] = t_0_zzzz_x_xxz[i] - rcd_z[i] * t_0_zzzz_x_xx[i];

        t_0_zzzz_xy_zz[i] = t_0_zzzz_x_yzz[i] - rcd_y[i] * t_0_zzzz_x_zz[i];

        t_0_zzzz_xy_yz[i] = t_0_zzzz_x_yyz[i] - rcd_y[i] * t_0_zzzz_x_yz[i];

        t_0_zzzz_xy_yy[i] = t_0_zzzz_x_yyy[i] - rcd_y[i] * t_0_zzzz_x_yy[i];

        t_0_zzzz_xy_xz[i] = t_0_zzzz_x_xyz[i] - rcd_y[i] * t_0_zzzz_x_xz[i];

        t_0_zzzz_xy_xy[i] = t_0_zzzz_x_xyy[i] - rcd_y[i] * t_0_zzzz_x_xy[i];

        t_0_zzzz_xy_xx[i] = t_0_zzzz_x_xxy[i] - rcd_y[i] * t_0_zzzz_x_xx[i];

        t_0_zzzz_xx_zz[i] = t_0_zzzz_x_xzz[i] - rcd_x[i] * t_0_zzzz_x_zz[i];

        t_0_zzzz_xx_yz[i] = t_0_zzzz_x_xyz[i] - rcd_x[i] * t_0_zzzz_x_yz[i];

        t_0_zzzz_xx_yy[i] = t_0_zzzz_x_xyy[i] - rcd_x[i] * t_0_zzzz_x_yy[i];

        t_0_zzzz_xx_xz[i] = t_0_zzzz_x_xxz[i] - rcd_x[i] * t_0_zzzz_x_xz[i];

        t_0_zzzz_xx_xy[i] = t_0_zzzz_x_xxy[i] - rcd_x[i] * t_0_zzzz_x_xy[i];

        t_0_zzzz_xx_xx[i] = t_0_zzzz_x_xxx[i] - rcd_x[i] * t_0_zzzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzzz_x_xx, t_0_yzzz_x_xxx, t_0_yzzz_x_xxy,\
                           t_0_yzzz_x_xxz, t_0_yzzz_x_xy, t_0_yzzz_x_xyy, t_0_yzzz_x_xyz,\
                           t_0_yzzz_x_xz, t_0_yzzz_x_xzz, t_0_yzzz_x_yy, t_0_yzzz_x_yyy,\
                           t_0_yzzz_x_yyz, t_0_yzzz_x_yz, t_0_yzzz_x_yzz, t_0_yzzz_x_zz,\
                           t_0_yzzz_x_zzz, t_0_yzzz_xx_xx, t_0_yzzz_xx_xy, t_0_yzzz_xx_xz,\
                           t_0_yzzz_xx_yy, t_0_yzzz_xx_yz, t_0_yzzz_xx_zz, t_0_yzzz_xy_xx,\
                           t_0_yzzz_xy_xy, t_0_yzzz_xy_xz, t_0_yzzz_xy_yy, t_0_yzzz_xy_yz,\
                           t_0_yzzz_xy_zz, t_0_yzzz_xz_xx, t_0_yzzz_xz_xy, t_0_yzzz_xz_xz,\
                           t_0_yzzz_xz_yy, t_0_yzzz_xz_yz, t_0_yzzz_xz_zz, t_0_yzzz_y_xx,\
                           t_0_yzzz_y_xxy, t_0_yzzz_y_xxz, t_0_yzzz_y_xy, t_0_yzzz_y_xyy,\
                           t_0_yzzz_y_xyz, t_0_yzzz_y_xz, t_0_yzzz_y_xzz, t_0_yzzz_y_yy,\
                           t_0_yzzz_y_yyy, t_0_yzzz_y_yyz, t_0_yzzz_y_yz, t_0_yzzz_y_yzz,\
                           t_0_yzzz_y_zz, t_0_yzzz_y_zzz, t_0_yzzz_yy_xx, t_0_yzzz_yy_xy,\
                           t_0_yzzz_yy_xz, t_0_yzzz_yy_yy, t_0_yzzz_yy_yz, t_0_yzzz_yy_zz,\
                           t_0_yzzz_yz_xx, t_0_yzzz_yz_xy, t_0_yzzz_yz_xz, t_0_yzzz_yz_yy,\
                           t_0_yzzz_yz_yz, t_0_yzzz_yz_zz, t_0_yzzz_z_xx, t_0_yzzz_z_xxz,\
                           t_0_yzzz_z_xy, t_0_yzzz_z_xyz, t_0_yzzz_z_xz, t_0_yzzz_z_xzz,\
                           t_0_yzzz_z_yy, t_0_yzzz_z_yyz, t_0_yzzz_z_yz, t_0_yzzz_z_yzz,\
                           t_0_yzzz_z_zz, t_0_yzzz_z_zzz, t_0_yzzz_zz_xx, t_0_yzzz_zz_xy,\
                           t_0_yzzz_zz_xz, t_0_yzzz_zz_yy, t_0_yzzz_zz_yz, t_0_yzzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yzzz_zz_zz[i] = t_0_yzzz_z_zzz[i] - rcd_z[i] * t_0_yzzz_z_zz[i];

        t_0_yzzz_zz_yz[i] = t_0_yzzz_z_yzz[i] - rcd_z[i] * t_0_yzzz_z_yz[i];

        t_0_yzzz_zz_yy[i] = t_0_yzzz_z_yyz[i] - rcd_z[i] * t_0_yzzz_z_yy[i];

        t_0_yzzz_zz_xz[i] = t_0_yzzz_z_xzz[i] - rcd_z[i] * t_0_yzzz_z_xz[i];

        t_0_yzzz_zz_xy[i] = t_0_yzzz_z_xyz[i] - rcd_z[i] * t_0_yzzz_z_xy[i];

        t_0_yzzz_zz_xx[i] = t_0_yzzz_z_xxz[i] - rcd_z[i] * t_0_yzzz_z_xx[i];

        t_0_yzzz_yz_zz[i] = t_0_yzzz_y_zzz[i] - rcd_z[i] * t_0_yzzz_y_zz[i];

        t_0_yzzz_yz_yz[i] = t_0_yzzz_y_yzz[i] - rcd_z[i] * t_0_yzzz_y_yz[i];

        t_0_yzzz_yz_yy[i] = t_0_yzzz_y_yyz[i] - rcd_z[i] * t_0_yzzz_y_yy[i];

        t_0_yzzz_yz_xz[i] = t_0_yzzz_y_xzz[i] - rcd_z[i] * t_0_yzzz_y_xz[i];

        t_0_yzzz_yz_xy[i] = t_0_yzzz_y_xyz[i] - rcd_z[i] * t_0_yzzz_y_xy[i];

        t_0_yzzz_yz_xx[i] = t_0_yzzz_y_xxz[i] - rcd_z[i] * t_0_yzzz_y_xx[i];

        t_0_yzzz_yy_zz[i] = t_0_yzzz_y_yzz[i] - rcd_y[i] * t_0_yzzz_y_zz[i];

        t_0_yzzz_yy_yz[i] = t_0_yzzz_y_yyz[i] - rcd_y[i] * t_0_yzzz_y_yz[i];

        t_0_yzzz_yy_yy[i] = t_0_yzzz_y_yyy[i] - rcd_y[i] * t_0_yzzz_y_yy[i];

        t_0_yzzz_yy_xz[i] = t_0_yzzz_y_xyz[i] - rcd_y[i] * t_0_yzzz_y_xz[i];

        t_0_yzzz_yy_xy[i] = t_0_yzzz_y_xyy[i] - rcd_y[i] * t_0_yzzz_y_xy[i];

        t_0_yzzz_yy_xx[i] = t_0_yzzz_y_xxy[i] - rcd_y[i] * t_0_yzzz_y_xx[i];

        t_0_yzzz_xz_zz[i] = t_0_yzzz_x_zzz[i] - rcd_z[i] * t_0_yzzz_x_zz[i];

        t_0_yzzz_xz_yz[i] = t_0_yzzz_x_yzz[i] - rcd_z[i] * t_0_yzzz_x_yz[i];

        t_0_yzzz_xz_yy[i] = t_0_yzzz_x_yyz[i] - rcd_z[i] * t_0_yzzz_x_yy[i];

        t_0_yzzz_xz_xz[i] = t_0_yzzz_x_xzz[i] - rcd_z[i] * t_0_yzzz_x_xz[i];

        t_0_yzzz_xz_xy[i] = t_0_yzzz_x_xyz[i] - rcd_z[i] * t_0_yzzz_x_xy[i];

        t_0_yzzz_xz_xx[i] = t_0_yzzz_x_xxz[i] - rcd_z[i] * t_0_yzzz_x_xx[i];

        t_0_yzzz_xy_zz[i] = t_0_yzzz_x_yzz[i] - rcd_y[i] * t_0_yzzz_x_zz[i];

        t_0_yzzz_xy_yz[i] = t_0_yzzz_x_yyz[i] - rcd_y[i] * t_0_yzzz_x_yz[i];

        t_0_yzzz_xy_yy[i] = t_0_yzzz_x_yyy[i] - rcd_y[i] * t_0_yzzz_x_yy[i];

        t_0_yzzz_xy_xz[i] = t_0_yzzz_x_xyz[i] - rcd_y[i] * t_0_yzzz_x_xz[i];

        t_0_yzzz_xy_xy[i] = t_0_yzzz_x_xyy[i] - rcd_y[i] * t_0_yzzz_x_xy[i];

        t_0_yzzz_xy_xx[i] = t_0_yzzz_x_xxy[i] - rcd_y[i] * t_0_yzzz_x_xx[i];

        t_0_yzzz_xx_zz[i] = t_0_yzzz_x_xzz[i] - rcd_x[i] * t_0_yzzz_x_zz[i];

        t_0_yzzz_xx_yz[i] = t_0_yzzz_x_xyz[i] - rcd_x[i] * t_0_yzzz_x_yz[i];

        t_0_yzzz_xx_yy[i] = t_0_yzzz_x_xyy[i] - rcd_x[i] * t_0_yzzz_x_yy[i];

        t_0_yzzz_xx_xz[i] = t_0_yzzz_x_xxz[i] - rcd_x[i] * t_0_yzzz_x_xz[i];

        t_0_yzzz_xx_xy[i] = t_0_yzzz_x_xxy[i] - rcd_x[i] * t_0_yzzz_x_xy[i];

        t_0_yzzz_xx_xx[i] = t_0_yzzz_x_xxx[i] - rcd_x[i] * t_0_yzzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyzz_x_xx, t_0_yyzz_x_xxx, t_0_yyzz_x_xxy,\
                           t_0_yyzz_x_xxz, t_0_yyzz_x_xy, t_0_yyzz_x_xyy, t_0_yyzz_x_xyz,\
                           t_0_yyzz_x_xz, t_0_yyzz_x_xzz, t_0_yyzz_x_yy, t_0_yyzz_x_yyy,\
                           t_0_yyzz_x_yyz, t_0_yyzz_x_yz, t_0_yyzz_x_yzz, t_0_yyzz_x_zz,\
                           t_0_yyzz_x_zzz, t_0_yyzz_xx_xx, t_0_yyzz_xx_xy, t_0_yyzz_xx_xz,\
                           t_0_yyzz_xx_yy, t_0_yyzz_xx_yz, t_0_yyzz_xx_zz, t_0_yyzz_xy_xx,\
                           t_0_yyzz_xy_xy, t_0_yyzz_xy_xz, t_0_yyzz_xy_yy, t_0_yyzz_xy_yz,\
                           t_0_yyzz_xy_zz, t_0_yyzz_xz_xx, t_0_yyzz_xz_xy, t_0_yyzz_xz_xz,\
                           t_0_yyzz_xz_yy, t_0_yyzz_xz_yz, t_0_yyzz_xz_zz, t_0_yyzz_y_xx,\
                           t_0_yyzz_y_xxy, t_0_yyzz_y_xxz, t_0_yyzz_y_xy, t_0_yyzz_y_xyy,\
                           t_0_yyzz_y_xyz, t_0_yyzz_y_xz, t_0_yyzz_y_xzz, t_0_yyzz_y_yy,\
                           t_0_yyzz_y_yyy, t_0_yyzz_y_yyz, t_0_yyzz_y_yz, t_0_yyzz_y_yzz,\
                           t_0_yyzz_y_zz, t_0_yyzz_y_zzz, t_0_yyzz_yy_xx, t_0_yyzz_yy_xy,\
                           t_0_yyzz_yy_xz, t_0_yyzz_yy_yy, t_0_yyzz_yy_yz, t_0_yyzz_yy_zz,\
                           t_0_yyzz_yz_xx, t_0_yyzz_yz_xy, t_0_yyzz_yz_xz, t_0_yyzz_yz_yy,\
                           t_0_yyzz_yz_yz, t_0_yyzz_yz_zz, t_0_yyzz_z_xx, t_0_yyzz_z_xxz,\
                           t_0_yyzz_z_xy, t_0_yyzz_z_xyz, t_0_yyzz_z_xz, t_0_yyzz_z_xzz,\
                           t_0_yyzz_z_yy, t_0_yyzz_z_yyz, t_0_yyzz_z_yz, t_0_yyzz_z_yzz,\
                           t_0_yyzz_z_zz, t_0_yyzz_z_zzz, t_0_yyzz_zz_xx, t_0_yyzz_zz_xy,\
                           t_0_yyzz_zz_xz, t_0_yyzz_zz_yy, t_0_yyzz_zz_yz, t_0_yyzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyzz_zz_zz[i] = t_0_yyzz_z_zzz[i] - rcd_z[i] * t_0_yyzz_z_zz[i];

        t_0_yyzz_zz_yz[i] = t_0_yyzz_z_yzz[i] - rcd_z[i] * t_0_yyzz_z_yz[i];

        t_0_yyzz_zz_yy[i] = t_0_yyzz_z_yyz[i] - rcd_z[i] * t_0_yyzz_z_yy[i];

        t_0_yyzz_zz_xz[i] = t_0_yyzz_z_xzz[i] - rcd_z[i] * t_0_yyzz_z_xz[i];

        t_0_yyzz_zz_xy[i] = t_0_yyzz_z_xyz[i] - rcd_z[i] * t_0_yyzz_z_xy[i];

        t_0_yyzz_zz_xx[i] = t_0_yyzz_z_xxz[i] - rcd_z[i] * t_0_yyzz_z_xx[i];

        t_0_yyzz_yz_zz[i] = t_0_yyzz_y_zzz[i] - rcd_z[i] * t_0_yyzz_y_zz[i];

        t_0_yyzz_yz_yz[i] = t_0_yyzz_y_yzz[i] - rcd_z[i] * t_0_yyzz_y_yz[i];

        t_0_yyzz_yz_yy[i] = t_0_yyzz_y_yyz[i] - rcd_z[i] * t_0_yyzz_y_yy[i];

        t_0_yyzz_yz_xz[i] = t_0_yyzz_y_xzz[i] - rcd_z[i] * t_0_yyzz_y_xz[i];

        t_0_yyzz_yz_xy[i] = t_0_yyzz_y_xyz[i] - rcd_z[i] * t_0_yyzz_y_xy[i];

        t_0_yyzz_yz_xx[i] = t_0_yyzz_y_xxz[i] - rcd_z[i] * t_0_yyzz_y_xx[i];

        t_0_yyzz_yy_zz[i] = t_0_yyzz_y_yzz[i] - rcd_y[i] * t_0_yyzz_y_zz[i];

        t_0_yyzz_yy_yz[i] = t_0_yyzz_y_yyz[i] - rcd_y[i] * t_0_yyzz_y_yz[i];

        t_0_yyzz_yy_yy[i] = t_0_yyzz_y_yyy[i] - rcd_y[i] * t_0_yyzz_y_yy[i];

        t_0_yyzz_yy_xz[i] = t_0_yyzz_y_xyz[i] - rcd_y[i] * t_0_yyzz_y_xz[i];

        t_0_yyzz_yy_xy[i] = t_0_yyzz_y_xyy[i] - rcd_y[i] * t_0_yyzz_y_xy[i];

        t_0_yyzz_yy_xx[i] = t_0_yyzz_y_xxy[i] - rcd_y[i] * t_0_yyzz_y_xx[i];

        t_0_yyzz_xz_zz[i] = t_0_yyzz_x_zzz[i] - rcd_z[i] * t_0_yyzz_x_zz[i];

        t_0_yyzz_xz_yz[i] = t_0_yyzz_x_yzz[i] - rcd_z[i] * t_0_yyzz_x_yz[i];

        t_0_yyzz_xz_yy[i] = t_0_yyzz_x_yyz[i] - rcd_z[i] * t_0_yyzz_x_yy[i];

        t_0_yyzz_xz_xz[i] = t_0_yyzz_x_xzz[i] - rcd_z[i] * t_0_yyzz_x_xz[i];

        t_0_yyzz_xz_xy[i] = t_0_yyzz_x_xyz[i] - rcd_z[i] * t_0_yyzz_x_xy[i];

        t_0_yyzz_xz_xx[i] = t_0_yyzz_x_xxz[i] - rcd_z[i] * t_0_yyzz_x_xx[i];

        t_0_yyzz_xy_zz[i] = t_0_yyzz_x_yzz[i] - rcd_y[i] * t_0_yyzz_x_zz[i];

        t_0_yyzz_xy_yz[i] = t_0_yyzz_x_yyz[i] - rcd_y[i] * t_0_yyzz_x_yz[i];

        t_0_yyzz_xy_yy[i] = t_0_yyzz_x_yyy[i] - rcd_y[i] * t_0_yyzz_x_yy[i];

        t_0_yyzz_xy_xz[i] = t_0_yyzz_x_xyz[i] - rcd_y[i] * t_0_yyzz_x_xz[i];

        t_0_yyzz_xy_xy[i] = t_0_yyzz_x_xyy[i] - rcd_y[i] * t_0_yyzz_x_xy[i];

        t_0_yyzz_xy_xx[i] = t_0_yyzz_x_xxy[i] - rcd_y[i] * t_0_yyzz_x_xx[i];

        t_0_yyzz_xx_zz[i] = t_0_yyzz_x_xzz[i] - rcd_x[i] * t_0_yyzz_x_zz[i];

        t_0_yyzz_xx_yz[i] = t_0_yyzz_x_xyz[i] - rcd_x[i] * t_0_yyzz_x_yz[i];

        t_0_yyzz_xx_yy[i] = t_0_yyzz_x_xyy[i] - rcd_x[i] * t_0_yyzz_x_yy[i];

        t_0_yyzz_xx_xz[i] = t_0_yyzz_x_xxz[i] - rcd_x[i] * t_0_yyzz_x_xz[i];

        t_0_yyzz_xx_xy[i] = t_0_yyzz_x_xxy[i] - rcd_x[i] * t_0_yyzz_x_xy[i];

        t_0_yyzz_xx_xx[i] = t_0_yyzz_x_xxx[i] - rcd_x[i] * t_0_yyzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyyz_x_xx, t_0_yyyz_x_xxx, t_0_yyyz_x_xxy,\
                           t_0_yyyz_x_xxz, t_0_yyyz_x_xy, t_0_yyyz_x_xyy, t_0_yyyz_x_xyz,\
                           t_0_yyyz_x_xz, t_0_yyyz_x_xzz, t_0_yyyz_x_yy, t_0_yyyz_x_yyy,\
                           t_0_yyyz_x_yyz, t_0_yyyz_x_yz, t_0_yyyz_x_yzz, t_0_yyyz_x_zz,\
                           t_0_yyyz_x_zzz, t_0_yyyz_xx_xx, t_0_yyyz_xx_xy, t_0_yyyz_xx_xz,\
                           t_0_yyyz_xx_yy, t_0_yyyz_xx_yz, t_0_yyyz_xx_zz, t_0_yyyz_xy_xx,\
                           t_0_yyyz_xy_xy, t_0_yyyz_xy_xz, t_0_yyyz_xy_yy, t_0_yyyz_xy_yz,\
                           t_0_yyyz_xy_zz, t_0_yyyz_xz_xx, t_0_yyyz_xz_xy, t_0_yyyz_xz_xz,\
                           t_0_yyyz_xz_yy, t_0_yyyz_xz_yz, t_0_yyyz_xz_zz, t_0_yyyz_y_xx,\
                           t_0_yyyz_y_xxy, t_0_yyyz_y_xxz, t_0_yyyz_y_xy, t_0_yyyz_y_xyy,\
                           t_0_yyyz_y_xyz, t_0_yyyz_y_xz, t_0_yyyz_y_xzz, t_0_yyyz_y_yy,\
                           t_0_yyyz_y_yyy, t_0_yyyz_y_yyz, t_0_yyyz_y_yz, t_0_yyyz_y_yzz,\
                           t_0_yyyz_y_zz, t_0_yyyz_y_zzz, t_0_yyyz_yy_xx, t_0_yyyz_yy_xy,\
                           t_0_yyyz_yy_xz, t_0_yyyz_yy_yy, t_0_yyyz_yy_yz, t_0_yyyz_yy_zz,\
                           t_0_yyyz_yz_xx, t_0_yyyz_yz_xy, t_0_yyyz_yz_xz, t_0_yyyz_yz_yy,\
                           t_0_yyyz_yz_yz, t_0_yyyz_yz_zz, t_0_yyyz_z_xx, t_0_yyyz_z_xxz,\
                           t_0_yyyz_z_xy, t_0_yyyz_z_xyz, t_0_yyyz_z_xz, t_0_yyyz_z_xzz,\
                           t_0_yyyz_z_yy, t_0_yyyz_z_yyz, t_0_yyyz_z_yz, t_0_yyyz_z_yzz,\
                           t_0_yyyz_z_zz, t_0_yyyz_z_zzz, t_0_yyyz_zz_xx, t_0_yyyz_zz_xy,\
                           t_0_yyyz_zz_xz, t_0_yyyz_zz_yy, t_0_yyyz_zz_yz, t_0_yyyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyyz_zz_zz[i] = t_0_yyyz_z_zzz[i] - rcd_z[i] * t_0_yyyz_z_zz[i];

        t_0_yyyz_zz_yz[i] = t_0_yyyz_z_yzz[i] - rcd_z[i] * t_0_yyyz_z_yz[i];

        t_0_yyyz_zz_yy[i] = t_0_yyyz_z_yyz[i] - rcd_z[i] * t_0_yyyz_z_yy[i];

        t_0_yyyz_zz_xz[i] = t_0_yyyz_z_xzz[i] - rcd_z[i] * t_0_yyyz_z_xz[i];

        t_0_yyyz_zz_xy[i] = t_0_yyyz_z_xyz[i] - rcd_z[i] * t_0_yyyz_z_xy[i];

        t_0_yyyz_zz_xx[i] = t_0_yyyz_z_xxz[i] - rcd_z[i] * t_0_yyyz_z_xx[i];

        t_0_yyyz_yz_zz[i] = t_0_yyyz_y_zzz[i] - rcd_z[i] * t_0_yyyz_y_zz[i];

        t_0_yyyz_yz_yz[i] = t_0_yyyz_y_yzz[i] - rcd_z[i] * t_0_yyyz_y_yz[i];

        t_0_yyyz_yz_yy[i] = t_0_yyyz_y_yyz[i] - rcd_z[i] * t_0_yyyz_y_yy[i];

        t_0_yyyz_yz_xz[i] = t_0_yyyz_y_xzz[i] - rcd_z[i] * t_0_yyyz_y_xz[i];

        t_0_yyyz_yz_xy[i] = t_0_yyyz_y_xyz[i] - rcd_z[i] * t_0_yyyz_y_xy[i];

        t_0_yyyz_yz_xx[i] = t_0_yyyz_y_xxz[i] - rcd_z[i] * t_0_yyyz_y_xx[i];

        t_0_yyyz_yy_zz[i] = t_0_yyyz_y_yzz[i] - rcd_y[i] * t_0_yyyz_y_zz[i];

        t_0_yyyz_yy_yz[i] = t_0_yyyz_y_yyz[i] - rcd_y[i] * t_0_yyyz_y_yz[i];

        t_0_yyyz_yy_yy[i] = t_0_yyyz_y_yyy[i] - rcd_y[i] * t_0_yyyz_y_yy[i];

        t_0_yyyz_yy_xz[i] = t_0_yyyz_y_xyz[i] - rcd_y[i] * t_0_yyyz_y_xz[i];

        t_0_yyyz_yy_xy[i] = t_0_yyyz_y_xyy[i] - rcd_y[i] * t_0_yyyz_y_xy[i];

        t_0_yyyz_yy_xx[i] = t_0_yyyz_y_xxy[i] - rcd_y[i] * t_0_yyyz_y_xx[i];

        t_0_yyyz_xz_zz[i] = t_0_yyyz_x_zzz[i] - rcd_z[i] * t_0_yyyz_x_zz[i];

        t_0_yyyz_xz_yz[i] = t_0_yyyz_x_yzz[i] - rcd_z[i] * t_0_yyyz_x_yz[i];

        t_0_yyyz_xz_yy[i] = t_0_yyyz_x_yyz[i] - rcd_z[i] * t_0_yyyz_x_yy[i];

        t_0_yyyz_xz_xz[i] = t_0_yyyz_x_xzz[i] - rcd_z[i] * t_0_yyyz_x_xz[i];

        t_0_yyyz_xz_xy[i] = t_0_yyyz_x_xyz[i] - rcd_z[i] * t_0_yyyz_x_xy[i];

        t_0_yyyz_xz_xx[i] = t_0_yyyz_x_xxz[i] - rcd_z[i] * t_0_yyyz_x_xx[i];

        t_0_yyyz_xy_zz[i] = t_0_yyyz_x_yzz[i] - rcd_y[i] * t_0_yyyz_x_zz[i];

        t_0_yyyz_xy_yz[i] = t_0_yyyz_x_yyz[i] - rcd_y[i] * t_0_yyyz_x_yz[i];

        t_0_yyyz_xy_yy[i] = t_0_yyyz_x_yyy[i] - rcd_y[i] * t_0_yyyz_x_yy[i];

        t_0_yyyz_xy_xz[i] = t_0_yyyz_x_xyz[i] - rcd_y[i] * t_0_yyyz_x_xz[i];

        t_0_yyyz_xy_xy[i] = t_0_yyyz_x_xyy[i] - rcd_y[i] * t_0_yyyz_x_xy[i];

        t_0_yyyz_xy_xx[i] = t_0_yyyz_x_xxy[i] - rcd_y[i] * t_0_yyyz_x_xx[i];

        t_0_yyyz_xx_zz[i] = t_0_yyyz_x_xzz[i] - rcd_x[i] * t_0_yyyz_x_zz[i];

        t_0_yyyz_xx_yz[i] = t_0_yyyz_x_xyz[i] - rcd_x[i] * t_0_yyyz_x_yz[i];

        t_0_yyyz_xx_yy[i] = t_0_yyyz_x_xyy[i] - rcd_x[i] * t_0_yyyz_x_yy[i];

        t_0_yyyz_xx_xz[i] = t_0_yyyz_x_xxz[i] - rcd_x[i] * t_0_yyyz_x_xz[i];

        t_0_yyyz_xx_xy[i] = t_0_yyyz_x_xxy[i] - rcd_x[i] * t_0_yyyz_x_xy[i];

        t_0_yyyz_xx_xx[i] = t_0_yyyz_x_xxx[i] - rcd_x[i] * t_0_yyyz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyyy_x_xx, t_0_yyyy_x_xxx, t_0_yyyy_x_xxy,\
                           t_0_yyyy_x_xxz, t_0_yyyy_x_xy, t_0_yyyy_x_xyy, t_0_yyyy_x_xyz,\
                           t_0_yyyy_x_xz, t_0_yyyy_x_xzz, t_0_yyyy_x_yy, t_0_yyyy_x_yyy,\
                           t_0_yyyy_x_yyz, t_0_yyyy_x_yz, t_0_yyyy_x_yzz, t_0_yyyy_x_zz,\
                           t_0_yyyy_x_zzz, t_0_yyyy_xx_xx, t_0_yyyy_xx_xy, t_0_yyyy_xx_xz,\
                           t_0_yyyy_xx_yy, t_0_yyyy_xx_yz, t_0_yyyy_xx_zz, t_0_yyyy_xy_xx,\
                           t_0_yyyy_xy_xy, t_0_yyyy_xy_xz, t_0_yyyy_xy_yy, t_0_yyyy_xy_yz,\
                           t_0_yyyy_xy_zz, t_0_yyyy_xz_xx, t_0_yyyy_xz_xy, t_0_yyyy_xz_xz,\
                           t_0_yyyy_xz_yy, t_0_yyyy_xz_yz, t_0_yyyy_xz_zz, t_0_yyyy_y_xx,\
                           t_0_yyyy_y_xxy, t_0_yyyy_y_xxz, t_0_yyyy_y_xy, t_0_yyyy_y_xyy,\
                           t_0_yyyy_y_xyz, t_0_yyyy_y_xz, t_0_yyyy_y_xzz, t_0_yyyy_y_yy,\
                           t_0_yyyy_y_yyy, t_0_yyyy_y_yyz, t_0_yyyy_y_yz, t_0_yyyy_y_yzz,\
                           t_0_yyyy_y_zz, t_0_yyyy_y_zzz, t_0_yyyy_yy_xx, t_0_yyyy_yy_xy,\
                           t_0_yyyy_yy_xz, t_0_yyyy_yy_yy, t_0_yyyy_yy_yz, t_0_yyyy_yy_zz,\
                           t_0_yyyy_yz_xx, t_0_yyyy_yz_xy, t_0_yyyy_yz_xz, t_0_yyyy_yz_yy,\
                           t_0_yyyy_yz_yz, t_0_yyyy_yz_zz, t_0_yyyy_z_xx, t_0_yyyy_z_xxz,\
                           t_0_yyyy_z_xy, t_0_yyyy_z_xyz, t_0_yyyy_z_xz, t_0_yyyy_z_xzz,\
                           t_0_yyyy_z_yy, t_0_yyyy_z_yyz, t_0_yyyy_z_yz, t_0_yyyy_z_yzz,\
                           t_0_yyyy_z_zz, t_0_yyyy_z_zzz, t_0_yyyy_zz_xx, t_0_yyyy_zz_xy,\
                           t_0_yyyy_zz_xz, t_0_yyyy_zz_yy, t_0_yyyy_zz_yz, t_0_yyyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyyy_zz_zz[i] = t_0_yyyy_z_zzz[i] - rcd_z[i] * t_0_yyyy_z_zz[i];

        t_0_yyyy_zz_yz[i] = t_0_yyyy_z_yzz[i] - rcd_z[i] * t_0_yyyy_z_yz[i];

        t_0_yyyy_zz_yy[i] = t_0_yyyy_z_yyz[i] - rcd_z[i] * t_0_yyyy_z_yy[i];

        t_0_yyyy_zz_xz[i] = t_0_yyyy_z_xzz[i] - rcd_z[i] * t_0_yyyy_z_xz[i];

        t_0_yyyy_zz_xy[i] = t_0_yyyy_z_xyz[i] - rcd_z[i] * t_0_yyyy_z_xy[i];

        t_0_yyyy_zz_xx[i] = t_0_yyyy_z_xxz[i] - rcd_z[i] * t_0_yyyy_z_xx[i];

        t_0_yyyy_yz_zz[i] = t_0_yyyy_y_zzz[i] - rcd_z[i] * t_0_yyyy_y_zz[i];

        t_0_yyyy_yz_yz[i] = t_0_yyyy_y_yzz[i] - rcd_z[i] * t_0_yyyy_y_yz[i];

        t_0_yyyy_yz_yy[i] = t_0_yyyy_y_yyz[i] - rcd_z[i] * t_0_yyyy_y_yy[i];

        t_0_yyyy_yz_xz[i] = t_0_yyyy_y_xzz[i] - rcd_z[i] * t_0_yyyy_y_xz[i];

        t_0_yyyy_yz_xy[i] = t_0_yyyy_y_xyz[i] - rcd_z[i] * t_0_yyyy_y_xy[i];

        t_0_yyyy_yz_xx[i] = t_0_yyyy_y_xxz[i] - rcd_z[i] * t_0_yyyy_y_xx[i];

        t_0_yyyy_yy_zz[i] = t_0_yyyy_y_yzz[i] - rcd_y[i] * t_0_yyyy_y_zz[i];

        t_0_yyyy_yy_yz[i] = t_0_yyyy_y_yyz[i] - rcd_y[i] * t_0_yyyy_y_yz[i];

        t_0_yyyy_yy_yy[i] = t_0_yyyy_y_yyy[i] - rcd_y[i] * t_0_yyyy_y_yy[i];

        t_0_yyyy_yy_xz[i] = t_0_yyyy_y_xyz[i] - rcd_y[i] * t_0_yyyy_y_xz[i];

        t_0_yyyy_yy_xy[i] = t_0_yyyy_y_xyy[i] - rcd_y[i] * t_0_yyyy_y_xy[i];

        t_0_yyyy_yy_xx[i] = t_0_yyyy_y_xxy[i] - rcd_y[i] * t_0_yyyy_y_xx[i];

        t_0_yyyy_xz_zz[i] = t_0_yyyy_x_zzz[i] - rcd_z[i] * t_0_yyyy_x_zz[i];

        t_0_yyyy_xz_yz[i] = t_0_yyyy_x_yzz[i] - rcd_z[i] * t_0_yyyy_x_yz[i];

        t_0_yyyy_xz_yy[i] = t_0_yyyy_x_yyz[i] - rcd_z[i] * t_0_yyyy_x_yy[i];

        t_0_yyyy_xz_xz[i] = t_0_yyyy_x_xzz[i] - rcd_z[i] * t_0_yyyy_x_xz[i];

        t_0_yyyy_xz_xy[i] = t_0_yyyy_x_xyz[i] - rcd_z[i] * t_0_yyyy_x_xy[i];

        t_0_yyyy_xz_xx[i] = t_0_yyyy_x_xxz[i] - rcd_z[i] * t_0_yyyy_x_xx[i];

        t_0_yyyy_xy_zz[i] = t_0_yyyy_x_yzz[i] - rcd_y[i] * t_0_yyyy_x_zz[i];

        t_0_yyyy_xy_yz[i] = t_0_yyyy_x_yyz[i] - rcd_y[i] * t_0_yyyy_x_yz[i];

        t_0_yyyy_xy_yy[i] = t_0_yyyy_x_yyy[i] - rcd_y[i] * t_0_yyyy_x_yy[i];

        t_0_yyyy_xy_xz[i] = t_0_yyyy_x_xyz[i] - rcd_y[i] * t_0_yyyy_x_xz[i];

        t_0_yyyy_xy_xy[i] = t_0_yyyy_x_xyy[i] - rcd_y[i] * t_0_yyyy_x_xy[i];

        t_0_yyyy_xy_xx[i] = t_0_yyyy_x_xxy[i] - rcd_y[i] * t_0_yyyy_x_xx[i];

        t_0_yyyy_xx_zz[i] = t_0_yyyy_x_xzz[i] - rcd_x[i] * t_0_yyyy_x_zz[i];

        t_0_yyyy_xx_yz[i] = t_0_yyyy_x_xyz[i] - rcd_x[i] * t_0_yyyy_x_yz[i];

        t_0_yyyy_xx_yy[i] = t_0_yyyy_x_xyy[i] - rcd_x[i] * t_0_yyyy_x_yy[i];

        t_0_yyyy_xx_xz[i] = t_0_yyyy_x_xxz[i] - rcd_x[i] * t_0_yyyy_x_xz[i];

        t_0_yyyy_xx_xy[i] = t_0_yyyy_x_xxy[i] - rcd_x[i] * t_0_yyyy_x_xy[i];

        t_0_yyyy_xx_xx[i] = t_0_yyyy_x_xxx[i] - rcd_x[i] * t_0_yyyy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xzzz_x_xx, t_0_xzzz_x_xxx, t_0_xzzz_x_xxy,\
                           t_0_xzzz_x_xxz, t_0_xzzz_x_xy, t_0_xzzz_x_xyy, t_0_xzzz_x_xyz,\
                           t_0_xzzz_x_xz, t_0_xzzz_x_xzz, t_0_xzzz_x_yy, t_0_xzzz_x_yyy,\
                           t_0_xzzz_x_yyz, t_0_xzzz_x_yz, t_0_xzzz_x_yzz, t_0_xzzz_x_zz,\
                           t_0_xzzz_x_zzz, t_0_xzzz_xx_xx, t_0_xzzz_xx_xy, t_0_xzzz_xx_xz,\
                           t_0_xzzz_xx_yy, t_0_xzzz_xx_yz, t_0_xzzz_xx_zz, t_0_xzzz_xy_xx,\
                           t_0_xzzz_xy_xy, t_0_xzzz_xy_xz, t_0_xzzz_xy_yy, t_0_xzzz_xy_yz,\
                           t_0_xzzz_xy_zz, t_0_xzzz_xz_xx, t_0_xzzz_xz_xy, t_0_xzzz_xz_xz,\
                           t_0_xzzz_xz_yy, t_0_xzzz_xz_yz, t_0_xzzz_xz_zz, t_0_xzzz_y_xx,\
                           t_0_xzzz_y_xxy, t_0_xzzz_y_xxz, t_0_xzzz_y_xy, t_0_xzzz_y_xyy,\
                           t_0_xzzz_y_xyz, t_0_xzzz_y_xz, t_0_xzzz_y_xzz, t_0_xzzz_y_yy,\
                           t_0_xzzz_y_yyy, t_0_xzzz_y_yyz, t_0_xzzz_y_yz, t_0_xzzz_y_yzz,\
                           t_0_xzzz_y_zz, t_0_xzzz_y_zzz, t_0_xzzz_yy_xx, t_0_xzzz_yy_xy,\
                           t_0_xzzz_yy_xz, t_0_xzzz_yy_yy, t_0_xzzz_yy_yz, t_0_xzzz_yy_zz,\
                           t_0_xzzz_yz_xx, t_0_xzzz_yz_xy, t_0_xzzz_yz_xz, t_0_xzzz_yz_yy,\
                           t_0_xzzz_yz_yz, t_0_xzzz_yz_zz, t_0_xzzz_z_xx, t_0_xzzz_z_xxz,\
                           t_0_xzzz_z_xy, t_0_xzzz_z_xyz, t_0_xzzz_z_xz, t_0_xzzz_z_xzz,\
                           t_0_xzzz_z_yy, t_0_xzzz_z_yyz, t_0_xzzz_z_yz, t_0_xzzz_z_yzz,\
                           t_0_xzzz_z_zz, t_0_xzzz_z_zzz, t_0_xzzz_zz_xx, t_0_xzzz_zz_xy,\
                           t_0_xzzz_zz_xz, t_0_xzzz_zz_yy, t_0_xzzz_zz_yz, t_0_xzzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzzz_zz_zz[i] = t_0_xzzz_z_zzz[i] - rcd_z[i] * t_0_xzzz_z_zz[i];

        t_0_xzzz_zz_yz[i] = t_0_xzzz_z_yzz[i] - rcd_z[i] * t_0_xzzz_z_yz[i];

        t_0_xzzz_zz_yy[i] = t_0_xzzz_z_yyz[i] - rcd_z[i] * t_0_xzzz_z_yy[i];

        t_0_xzzz_zz_xz[i] = t_0_xzzz_z_xzz[i] - rcd_z[i] * t_0_xzzz_z_xz[i];

        t_0_xzzz_zz_xy[i] = t_0_xzzz_z_xyz[i] - rcd_z[i] * t_0_xzzz_z_xy[i];

        t_0_xzzz_zz_xx[i] = t_0_xzzz_z_xxz[i] - rcd_z[i] * t_0_xzzz_z_xx[i];

        t_0_xzzz_yz_zz[i] = t_0_xzzz_y_zzz[i] - rcd_z[i] * t_0_xzzz_y_zz[i];

        t_0_xzzz_yz_yz[i] = t_0_xzzz_y_yzz[i] - rcd_z[i] * t_0_xzzz_y_yz[i];

        t_0_xzzz_yz_yy[i] = t_0_xzzz_y_yyz[i] - rcd_z[i] * t_0_xzzz_y_yy[i];

        t_0_xzzz_yz_xz[i] = t_0_xzzz_y_xzz[i] - rcd_z[i] * t_0_xzzz_y_xz[i];

        t_0_xzzz_yz_xy[i] = t_0_xzzz_y_xyz[i] - rcd_z[i] * t_0_xzzz_y_xy[i];

        t_0_xzzz_yz_xx[i] = t_0_xzzz_y_xxz[i] - rcd_z[i] * t_0_xzzz_y_xx[i];

        t_0_xzzz_yy_zz[i] = t_0_xzzz_y_yzz[i] - rcd_y[i] * t_0_xzzz_y_zz[i];

        t_0_xzzz_yy_yz[i] = t_0_xzzz_y_yyz[i] - rcd_y[i] * t_0_xzzz_y_yz[i];

        t_0_xzzz_yy_yy[i] = t_0_xzzz_y_yyy[i] - rcd_y[i] * t_0_xzzz_y_yy[i];

        t_0_xzzz_yy_xz[i] = t_0_xzzz_y_xyz[i] - rcd_y[i] * t_0_xzzz_y_xz[i];

        t_0_xzzz_yy_xy[i] = t_0_xzzz_y_xyy[i] - rcd_y[i] * t_0_xzzz_y_xy[i];

        t_0_xzzz_yy_xx[i] = t_0_xzzz_y_xxy[i] - rcd_y[i] * t_0_xzzz_y_xx[i];

        t_0_xzzz_xz_zz[i] = t_0_xzzz_x_zzz[i] - rcd_z[i] * t_0_xzzz_x_zz[i];

        t_0_xzzz_xz_yz[i] = t_0_xzzz_x_yzz[i] - rcd_z[i] * t_0_xzzz_x_yz[i];

        t_0_xzzz_xz_yy[i] = t_0_xzzz_x_yyz[i] - rcd_z[i] * t_0_xzzz_x_yy[i];

        t_0_xzzz_xz_xz[i] = t_0_xzzz_x_xzz[i] - rcd_z[i] * t_0_xzzz_x_xz[i];

        t_0_xzzz_xz_xy[i] = t_0_xzzz_x_xyz[i] - rcd_z[i] * t_0_xzzz_x_xy[i];

        t_0_xzzz_xz_xx[i] = t_0_xzzz_x_xxz[i] - rcd_z[i] * t_0_xzzz_x_xx[i];

        t_0_xzzz_xy_zz[i] = t_0_xzzz_x_yzz[i] - rcd_y[i] * t_0_xzzz_x_zz[i];

        t_0_xzzz_xy_yz[i] = t_0_xzzz_x_yyz[i] - rcd_y[i] * t_0_xzzz_x_yz[i];

        t_0_xzzz_xy_yy[i] = t_0_xzzz_x_yyy[i] - rcd_y[i] * t_0_xzzz_x_yy[i];

        t_0_xzzz_xy_xz[i] = t_0_xzzz_x_xyz[i] - rcd_y[i] * t_0_xzzz_x_xz[i];

        t_0_xzzz_xy_xy[i] = t_0_xzzz_x_xyy[i] - rcd_y[i] * t_0_xzzz_x_xy[i];

        t_0_xzzz_xy_xx[i] = t_0_xzzz_x_xxy[i] - rcd_y[i] * t_0_xzzz_x_xx[i];

        t_0_xzzz_xx_zz[i] = t_0_xzzz_x_xzz[i] - rcd_x[i] * t_0_xzzz_x_zz[i];

        t_0_xzzz_xx_yz[i] = t_0_xzzz_x_xyz[i] - rcd_x[i] * t_0_xzzz_x_yz[i];

        t_0_xzzz_xx_yy[i] = t_0_xzzz_x_xyy[i] - rcd_x[i] * t_0_xzzz_x_yy[i];

        t_0_xzzz_xx_xz[i] = t_0_xzzz_x_xxz[i] - rcd_x[i] * t_0_xzzz_x_xz[i];

        t_0_xzzz_xx_xy[i] = t_0_xzzz_x_xxy[i] - rcd_x[i] * t_0_xzzz_x_xy[i];

        t_0_xzzz_xx_xx[i] = t_0_xzzz_x_xxx[i] - rcd_x[i] * t_0_xzzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyzz_x_xx, t_0_xyzz_x_xxx, t_0_xyzz_x_xxy,\
                           t_0_xyzz_x_xxz, t_0_xyzz_x_xy, t_0_xyzz_x_xyy, t_0_xyzz_x_xyz,\
                           t_0_xyzz_x_xz, t_0_xyzz_x_xzz, t_0_xyzz_x_yy, t_0_xyzz_x_yyy,\
                           t_0_xyzz_x_yyz, t_0_xyzz_x_yz, t_0_xyzz_x_yzz, t_0_xyzz_x_zz,\
                           t_0_xyzz_x_zzz, t_0_xyzz_xx_xx, t_0_xyzz_xx_xy, t_0_xyzz_xx_xz,\
                           t_0_xyzz_xx_yy, t_0_xyzz_xx_yz, t_0_xyzz_xx_zz, t_0_xyzz_xy_xx,\
                           t_0_xyzz_xy_xy, t_0_xyzz_xy_xz, t_0_xyzz_xy_yy, t_0_xyzz_xy_yz,\
                           t_0_xyzz_xy_zz, t_0_xyzz_xz_xx, t_0_xyzz_xz_xy, t_0_xyzz_xz_xz,\
                           t_0_xyzz_xz_yy, t_0_xyzz_xz_yz, t_0_xyzz_xz_zz, t_0_xyzz_y_xx,\
                           t_0_xyzz_y_xxy, t_0_xyzz_y_xxz, t_0_xyzz_y_xy, t_0_xyzz_y_xyy,\
                           t_0_xyzz_y_xyz, t_0_xyzz_y_xz, t_0_xyzz_y_xzz, t_0_xyzz_y_yy,\
                           t_0_xyzz_y_yyy, t_0_xyzz_y_yyz, t_0_xyzz_y_yz, t_0_xyzz_y_yzz,\
                           t_0_xyzz_y_zz, t_0_xyzz_y_zzz, t_0_xyzz_yy_xx, t_0_xyzz_yy_xy,\
                           t_0_xyzz_yy_xz, t_0_xyzz_yy_yy, t_0_xyzz_yy_yz, t_0_xyzz_yy_zz,\
                           t_0_xyzz_yz_xx, t_0_xyzz_yz_xy, t_0_xyzz_yz_xz, t_0_xyzz_yz_yy,\
                           t_0_xyzz_yz_yz, t_0_xyzz_yz_zz, t_0_xyzz_z_xx, t_0_xyzz_z_xxz,\
                           t_0_xyzz_z_xy, t_0_xyzz_z_xyz, t_0_xyzz_z_xz, t_0_xyzz_z_xzz,\
                           t_0_xyzz_z_yy, t_0_xyzz_z_yyz, t_0_xyzz_z_yz, t_0_xyzz_z_yzz,\
                           t_0_xyzz_z_zz, t_0_xyzz_z_zzz, t_0_xyzz_zz_xx, t_0_xyzz_zz_xy,\
                           t_0_xyzz_zz_xz, t_0_xyzz_zz_yy, t_0_xyzz_zz_yz, t_0_xyzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyzz_zz_zz[i] = t_0_xyzz_z_zzz[i] - rcd_z[i] * t_0_xyzz_z_zz[i];

        t_0_xyzz_zz_yz[i] = t_0_xyzz_z_yzz[i] - rcd_z[i] * t_0_xyzz_z_yz[i];

        t_0_xyzz_zz_yy[i] = t_0_xyzz_z_yyz[i] - rcd_z[i] * t_0_xyzz_z_yy[i];

        t_0_xyzz_zz_xz[i] = t_0_xyzz_z_xzz[i] - rcd_z[i] * t_0_xyzz_z_xz[i];

        t_0_xyzz_zz_xy[i] = t_0_xyzz_z_xyz[i] - rcd_z[i] * t_0_xyzz_z_xy[i];

        t_0_xyzz_zz_xx[i] = t_0_xyzz_z_xxz[i] - rcd_z[i] * t_0_xyzz_z_xx[i];

        t_0_xyzz_yz_zz[i] = t_0_xyzz_y_zzz[i] - rcd_z[i] * t_0_xyzz_y_zz[i];

        t_0_xyzz_yz_yz[i] = t_0_xyzz_y_yzz[i] - rcd_z[i] * t_0_xyzz_y_yz[i];

        t_0_xyzz_yz_yy[i] = t_0_xyzz_y_yyz[i] - rcd_z[i] * t_0_xyzz_y_yy[i];

        t_0_xyzz_yz_xz[i] = t_0_xyzz_y_xzz[i] - rcd_z[i] * t_0_xyzz_y_xz[i];

        t_0_xyzz_yz_xy[i] = t_0_xyzz_y_xyz[i] - rcd_z[i] * t_0_xyzz_y_xy[i];

        t_0_xyzz_yz_xx[i] = t_0_xyzz_y_xxz[i] - rcd_z[i] * t_0_xyzz_y_xx[i];

        t_0_xyzz_yy_zz[i] = t_0_xyzz_y_yzz[i] - rcd_y[i] * t_0_xyzz_y_zz[i];

        t_0_xyzz_yy_yz[i] = t_0_xyzz_y_yyz[i] - rcd_y[i] * t_0_xyzz_y_yz[i];

        t_0_xyzz_yy_yy[i] = t_0_xyzz_y_yyy[i] - rcd_y[i] * t_0_xyzz_y_yy[i];

        t_0_xyzz_yy_xz[i] = t_0_xyzz_y_xyz[i] - rcd_y[i] * t_0_xyzz_y_xz[i];

        t_0_xyzz_yy_xy[i] = t_0_xyzz_y_xyy[i] - rcd_y[i] * t_0_xyzz_y_xy[i];

        t_0_xyzz_yy_xx[i] = t_0_xyzz_y_xxy[i] - rcd_y[i] * t_0_xyzz_y_xx[i];

        t_0_xyzz_xz_zz[i] = t_0_xyzz_x_zzz[i] - rcd_z[i] * t_0_xyzz_x_zz[i];

        t_0_xyzz_xz_yz[i] = t_0_xyzz_x_yzz[i] - rcd_z[i] * t_0_xyzz_x_yz[i];

        t_0_xyzz_xz_yy[i] = t_0_xyzz_x_yyz[i] - rcd_z[i] * t_0_xyzz_x_yy[i];

        t_0_xyzz_xz_xz[i] = t_0_xyzz_x_xzz[i] - rcd_z[i] * t_0_xyzz_x_xz[i];

        t_0_xyzz_xz_xy[i] = t_0_xyzz_x_xyz[i] - rcd_z[i] * t_0_xyzz_x_xy[i];

        t_0_xyzz_xz_xx[i] = t_0_xyzz_x_xxz[i] - rcd_z[i] * t_0_xyzz_x_xx[i];

        t_0_xyzz_xy_zz[i] = t_0_xyzz_x_yzz[i] - rcd_y[i] * t_0_xyzz_x_zz[i];

        t_0_xyzz_xy_yz[i] = t_0_xyzz_x_yyz[i] - rcd_y[i] * t_0_xyzz_x_yz[i];

        t_0_xyzz_xy_yy[i] = t_0_xyzz_x_yyy[i] - rcd_y[i] * t_0_xyzz_x_yy[i];

        t_0_xyzz_xy_xz[i] = t_0_xyzz_x_xyz[i] - rcd_y[i] * t_0_xyzz_x_xz[i];

        t_0_xyzz_xy_xy[i] = t_0_xyzz_x_xyy[i] - rcd_y[i] * t_0_xyzz_x_xy[i];

        t_0_xyzz_xy_xx[i] = t_0_xyzz_x_xxy[i] - rcd_y[i] * t_0_xyzz_x_xx[i];

        t_0_xyzz_xx_zz[i] = t_0_xyzz_x_xzz[i] - rcd_x[i] * t_0_xyzz_x_zz[i];

        t_0_xyzz_xx_yz[i] = t_0_xyzz_x_xyz[i] - rcd_x[i] * t_0_xyzz_x_yz[i];

        t_0_xyzz_xx_yy[i] = t_0_xyzz_x_xyy[i] - rcd_x[i] * t_0_xyzz_x_yy[i];

        t_0_xyzz_xx_xz[i] = t_0_xyzz_x_xxz[i] - rcd_x[i] * t_0_xyzz_x_xz[i];

        t_0_xyzz_xx_xy[i] = t_0_xyzz_x_xxy[i] - rcd_x[i] * t_0_xyzz_x_xy[i];

        t_0_xyzz_xx_xx[i] = t_0_xyzz_x_xxx[i] - rcd_x[i] * t_0_xyzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyz_x_xx, t_0_xyyz_x_xxx, t_0_xyyz_x_xxy,\
                           t_0_xyyz_x_xxz, t_0_xyyz_x_xy, t_0_xyyz_x_xyy, t_0_xyyz_x_xyz,\
                           t_0_xyyz_x_xz, t_0_xyyz_x_xzz, t_0_xyyz_x_yy, t_0_xyyz_x_yyy,\
                           t_0_xyyz_x_yyz, t_0_xyyz_x_yz, t_0_xyyz_x_yzz, t_0_xyyz_x_zz,\
                           t_0_xyyz_x_zzz, t_0_xyyz_xx_xx, t_0_xyyz_xx_xy, t_0_xyyz_xx_xz,\
                           t_0_xyyz_xx_yy, t_0_xyyz_xx_yz, t_0_xyyz_xx_zz, t_0_xyyz_xy_xx,\
                           t_0_xyyz_xy_xy, t_0_xyyz_xy_xz, t_0_xyyz_xy_yy, t_0_xyyz_xy_yz,\
                           t_0_xyyz_xy_zz, t_0_xyyz_xz_xx, t_0_xyyz_xz_xy, t_0_xyyz_xz_xz,\
                           t_0_xyyz_xz_yy, t_0_xyyz_xz_yz, t_0_xyyz_xz_zz, t_0_xyyz_y_xx,\
                           t_0_xyyz_y_xxy, t_0_xyyz_y_xxz, t_0_xyyz_y_xy, t_0_xyyz_y_xyy,\
                           t_0_xyyz_y_xyz, t_0_xyyz_y_xz, t_0_xyyz_y_xzz, t_0_xyyz_y_yy,\
                           t_0_xyyz_y_yyy, t_0_xyyz_y_yyz, t_0_xyyz_y_yz, t_0_xyyz_y_yzz,\
                           t_0_xyyz_y_zz, t_0_xyyz_y_zzz, t_0_xyyz_yy_xx, t_0_xyyz_yy_xy,\
                           t_0_xyyz_yy_xz, t_0_xyyz_yy_yy, t_0_xyyz_yy_yz, t_0_xyyz_yy_zz,\
                           t_0_xyyz_yz_xx, t_0_xyyz_yz_xy, t_0_xyyz_yz_xz, t_0_xyyz_yz_yy,\
                           t_0_xyyz_yz_yz, t_0_xyyz_yz_zz, t_0_xyyz_z_xx, t_0_xyyz_z_xxz,\
                           t_0_xyyz_z_xy, t_0_xyyz_z_xyz, t_0_xyyz_z_xz, t_0_xyyz_z_xzz,\
                           t_0_xyyz_z_yy, t_0_xyyz_z_yyz, t_0_xyyz_z_yz, t_0_xyyz_z_yzz,\
                           t_0_xyyz_z_zz, t_0_xyyz_z_zzz, t_0_xyyz_zz_xx, t_0_xyyz_zz_xy,\
                           t_0_xyyz_zz_xz, t_0_xyyz_zz_yy, t_0_xyyz_zz_yz, t_0_xyyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyz_zz_zz[i] = t_0_xyyz_z_zzz[i] - rcd_z[i] * t_0_xyyz_z_zz[i];

        t_0_xyyz_zz_yz[i] = t_0_xyyz_z_yzz[i] - rcd_z[i] * t_0_xyyz_z_yz[i];

        t_0_xyyz_zz_yy[i] = t_0_xyyz_z_yyz[i] - rcd_z[i] * t_0_xyyz_z_yy[i];

        t_0_xyyz_zz_xz[i] = t_0_xyyz_z_xzz[i] - rcd_z[i] * t_0_xyyz_z_xz[i];

        t_0_xyyz_zz_xy[i] = t_0_xyyz_z_xyz[i] - rcd_z[i] * t_0_xyyz_z_xy[i];

        t_0_xyyz_zz_xx[i] = t_0_xyyz_z_xxz[i] - rcd_z[i] * t_0_xyyz_z_xx[i];

        t_0_xyyz_yz_zz[i] = t_0_xyyz_y_zzz[i] - rcd_z[i] * t_0_xyyz_y_zz[i];

        t_0_xyyz_yz_yz[i] = t_0_xyyz_y_yzz[i] - rcd_z[i] * t_0_xyyz_y_yz[i];

        t_0_xyyz_yz_yy[i] = t_0_xyyz_y_yyz[i] - rcd_z[i] * t_0_xyyz_y_yy[i];

        t_0_xyyz_yz_xz[i] = t_0_xyyz_y_xzz[i] - rcd_z[i] * t_0_xyyz_y_xz[i];

        t_0_xyyz_yz_xy[i] = t_0_xyyz_y_xyz[i] - rcd_z[i] * t_0_xyyz_y_xy[i];

        t_0_xyyz_yz_xx[i] = t_0_xyyz_y_xxz[i] - rcd_z[i] * t_0_xyyz_y_xx[i];

        t_0_xyyz_yy_zz[i] = t_0_xyyz_y_yzz[i] - rcd_y[i] * t_0_xyyz_y_zz[i];

        t_0_xyyz_yy_yz[i] = t_0_xyyz_y_yyz[i] - rcd_y[i] * t_0_xyyz_y_yz[i];

        t_0_xyyz_yy_yy[i] = t_0_xyyz_y_yyy[i] - rcd_y[i] * t_0_xyyz_y_yy[i];

        t_0_xyyz_yy_xz[i] = t_0_xyyz_y_xyz[i] - rcd_y[i] * t_0_xyyz_y_xz[i];

        t_0_xyyz_yy_xy[i] = t_0_xyyz_y_xyy[i] - rcd_y[i] * t_0_xyyz_y_xy[i];

        t_0_xyyz_yy_xx[i] = t_0_xyyz_y_xxy[i] - rcd_y[i] * t_0_xyyz_y_xx[i];

        t_0_xyyz_xz_zz[i] = t_0_xyyz_x_zzz[i] - rcd_z[i] * t_0_xyyz_x_zz[i];

        t_0_xyyz_xz_yz[i] = t_0_xyyz_x_yzz[i] - rcd_z[i] * t_0_xyyz_x_yz[i];

        t_0_xyyz_xz_yy[i] = t_0_xyyz_x_yyz[i] - rcd_z[i] * t_0_xyyz_x_yy[i];

        t_0_xyyz_xz_xz[i] = t_0_xyyz_x_xzz[i] - rcd_z[i] * t_0_xyyz_x_xz[i];

        t_0_xyyz_xz_xy[i] = t_0_xyyz_x_xyz[i] - rcd_z[i] * t_0_xyyz_x_xy[i];

        t_0_xyyz_xz_xx[i] = t_0_xyyz_x_xxz[i] - rcd_z[i] * t_0_xyyz_x_xx[i];

        t_0_xyyz_xy_zz[i] = t_0_xyyz_x_yzz[i] - rcd_y[i] * t_0_xyyz_x_zz[i];

        t_0_xyyz_xy_yz[i] = t_0_xyyz_x_yyz[i] - rcd_y[i] * t_0_xyyz_x_yz[i];

        t_0_xyyz_xy_yy[i] = t_0_xyyz_x_yyy[i] - rcd_y[i] * t_0_xyyz_x_yy[i];

        t_0_xyyz_xy_xz[i] = t_0_xyyz_x_xyz[i] - rcd_y[i] * t_0_xyyz_x_xz[i];

        t_0_xyyz_xy_xy[i] = t_0_xyyz_x_xyy[i] - rcd_y[i] * t_0_xyyz_x_xy[i];

        t_0_xyyz_xy_xx[i] = t_0_xyyz_x_xxy[i] - rcd_y[i] * t_0_xyyz_x_xx[i];

        t_0_xyyz_xx_zz[i] = t_0_xyyz_x_xzz[i] - rcd_x[i] * t_0_xyyz_x_zz[i];

        t_0_xyyz_xx_yz[i] = t_0_xyyz_x_xyz[i] - rcd_x[i] * t_0_xyyz_x_yz[i];

        t_0_xyyz_xx_yy[i] = t_0_xyyz_x_xyy[i] - rcd_x[i] * t_0_xyyz_x_yy[i];

        t_0_xyyz_xx_xz[i] = t_0_xyyz_x_xxz[i] - rcd_x[i] * t_0_xyyz_x_xz[i];

        t_0_xyyz_xx_xy[i] = t_0_xyyz_x_xxy[i] - rcd_x[i] * t_0_xyyz_x_xy[i];

        t_0_xyyz_xx_xx[i] = t_0_xyyz_x_xxx[i] - rcd_x[i] * t_0_xyyz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyy_x_xx, t_0_xyyy_x_xxx, t_0_xyyy_x_xxy,\
                           t_0_xyyy_x_xxz, t_0_xyyy_x_xy, t_0_xyyy_x_xyy, t_0_xyyy_x_xyz,\
                           t_0_xyyy_x_xz, t_0_xyyy_x_xzz, t_0_xyyy_x_yy, t_0_xyyy_x_yyy,\
                           t_0_xyyy_x_yyz, t_0_xyyy_x_yz, t_0_xyyy_x_yzz, t_0_xyyy_x_zz,\
                           t_0_xyyy_x_zzz, t_0_xyyy_xx_xx, t_0_xyyy_xx_xy, t_0_xyyy_xx_xz,\
                           t_0_xyyy_xx_yy, t_0_xyyy_xx_yz, t_0_xyyy_xx_zz, t_0_xyyy_xy_xx,\
                           t_0_xyyy_xy_xy, t_0_xyyy_xy_xz, t_0_xyyy_xy_yy, t_0_xyyy_xy_yz,\
                           t_0_xyyy_xy_zz, t_0_xyyy_xz_xx, t_0_xyyy_xz_xy, t_0_xyyy_xz_xz,\
                           t_0_xyyy_xz_yy, t_0_xyyy_xz_yz, t_0_xyyy_xz_zz, t_0_xyyy_y_xx,\
                           t_0_xyyy_y_xxy, t_0_xyyy_y_xxz, t_0_xyyy_y_xy, t_0_xyyy_y_xyy,\
                           t_0_xyyy_y_xyz, t_0_xyyy_y_xz, t_0_xyyy_y_xzz, t_0_xyyy_y_yy,\
                           t_0_xyyy_y_yyy, t_0_xyyy_y_yyz, t_0_xyyy_y_yz, t_0_xyyy_y_yzz,\
                           t_0_xyyy_y_zz, t_0_xyyy_y_zzz, t_0_xyyy_yy_xx, t_0_xyyy_yy_xy,\
                           t_0_xyyy_yy_xz, t_0_xyyy_yy_yy, t_0_xyyy_yy_yz, t_0_xyyy_yy_zz,\
                           t_0_xyyy_yz_xx, t_0_xyyy_yz_xy, t_0_xyyy_yz_xz, t_0_xyyy_yz_yy,\
                           t_0_xyyy_yz_yz, t_0_xyyy_yz_zz, t_0_xyyy_z_xx, t_0_xyyy_z_xxz,\
                           t_0_xyyy_z_xy, t_0_xyyy_z_xyz, t_0_xyyy_z_xz, t_0_xyyy_z_xzz,\
                           t_0_xyyy_z_yy, t_0_xyyy_z_yyz, t_0_xyyy_z_yz, t_0_xyyy_z_yzz,\
                           t_0_xyyy_z_zz, t_0_xyyy_z_zzz, t_0_xyyy_zz_xx, t_0_xyyy_zz_xy,\
                           t_0_xyyy_zz_xz, t_0_xyyy_zz_yy, t_0_xyyy_zz_yz, t_0_xyyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyy_zz_zz[i] = t_0_xyyy_z_zzz[i] - rcd_z[i] * t_0_xyyy_z_zz[i];

        t_0_xyyy_zz_yz[i] = t_0_xyyy_z_yzz[i] - rcd_z[i] * t_0_xyyy_z_yz[i];

        t_0_xyyy_zz_yy[i] = t_0_xyyy_z_yyz[i] - rcd_z[i] * t_0_xyyy_z_yy[i];

        t_0_xyyy_zz_xz[i] = t_0_xyyy_z_xzz[i] - rcd_z[i] * t_0_xyyy_z_xz[i];

        t_0_xyyy_zz_xy[i] = t_0_xyyy_z_xyz[i] - rcd_z[i] * t_0_xyyy_z_xy[i];

        t_0_xyyy_zz_xx[i] = t_0_xyyy_z_xxz[i] - rcd_z[i] * t_0_xyyy_z_xx[i];

        t_0_xyyy_yz_zz[i] = t_0_xyyy_y_zzz[i] - rcd_z[i] * t_0_xyyy_y_zz[i];

        t_0_xyyy_yz_yz[i] = t_0_xyyy_y_yzz[i] - rcd_z[i] * t_0_xyyy_y_yz[i];

        t_0_xyyy_yz_yy[i] = t_0_xyyy_y_yyz[i] - rcd_z[i] * t_0_xyyy_y_yy[i];

        t_0_xyyy_yz_xz[i] = t_0_xyyy_y_xzz[i] - rcd_z[i] * t_0_xyyy_y_xz[i];

        t_0_xyyy_yz_xy[i] = t_0_xyyy_y_xyz[i] - rcd_z[i] * t_0_xyyy_y_xy[i];

        t_0_xyyy_yz_xx[i] = t_0_xyyy_y_xxz[i] - rcd_z[i] * t_0_xyyy_y_xx[i];

        t_0_xyyy_yy_zz[i] = t_0_xyyy_y_yzz[i] - rcd_y[i] * t_0_xyyy_y_zz[i];

        t_0_xyyy_yy_yz[i] = t_0_xyyy_y_yyz[i] - rcd_y[i] * t_0_xyyy_y_yz[i];

        t_0_xyyy_yy_yy[i] = t_0_xyyy_y_yyy[i] - rcd_y[i] * t_0_xyyy_y_yy[i];

        t_0_xyyy_yy_xz[i] = t_0_xyyy_y_xyz[i] - rcd_y[i] * t_0_xyyy_y_xz[i];

        t_0_xyyy_yy_xy[i] = t_0_xyyy_y_xyy[i] - rcd_y[i] * t_0_xyyy_y_xy[i];

        t_0_xyyy_yy_xx[i] = t_0_xyyy_y_xxy[i] - rcd_y[i] * t_0_xyyy_y_xx[i];

        t_0_xyyy_xz_zz[i] = t_0_xyyy_x_zzz[i] - rcd_z[i] * t_0_xyyy_x_zz[i];

        t_0_xyyy_xz_yz[i] = t_0_xyyy_x_yzz[i] - rcd_z[i] * t_0_xyyy_x_yz[i];

        t_0_xyyy_xz_yy[i] = t_0_xyyy_x_yyz[i] - rcd_z[i] * t_0_xyyy_x_yy[i];

        t_0_xyyy_xz_xz[i] = t_0_xyyy_x_xzz[i] - rcd_z[i] * t_0_xyyy_x_xz[i];

        t_0_xyyy_xz_xy[i] = t_0_xyyy_x_xyz[i] - rcd_z[i] * t_0_xyyy_x_xy[i];

        t_0_xyyy_xz_xx[i] = t_0_xyyy_x_xxz[i] - rcd_z[i] * t_0_xyyy_x_xx[i];

        t_0_xyyy_xy_zz[i] = t_0_xyyy_x_yzz[i] - rcd_y[i] * t_0_xyyy_x_zz[i];

        t_0_xyyy_xy_yz[i] = t_0_xyyy_x_yyz[i] - rcd_y[i] * t_0_xyyy_x_yz[i];

        t_0_xyyy_xy_yy[i] = t_0_xyyy_x_yyy[i] - rcd_y[i] * t_0_xyyy_x_yy[i];

        t_0_xyyy_xy_xz[i] = t_0_xyyy_x_xyz[i] - rcd_y[i] * t_0_xyyy_x_xz[i];

        t_0_xyyy_xy_xy[i] = t_0_xyyy_x_xyy[i] - rcd_y[i] * t_0_xyyy_x_xy[i];

        t_0_xyyy_xy_xx[i] = t_0_xyyy_x_xxy[i] - rcd_y[i] * t_0_xyyy_x_xx[i];

        t_0_xyyy_xx_zz[i] = t_0_xyyy_x_xzz[i] - rcd_x[i] * t_0_xyyy_x_zz[i];

        t_0_xyyy_xx_yz[i] = t_0_xyyy_x_xyz[i] - rcd_x[i] * t_0_xyyy_x_yz[i];

        t_0_xyyy_xx_yy[i] = t_0_xyyy_x_xyy[i] - rcd_x[i] * t_0_xyyy_x_yy[i];

        t_0_xyyy_xx_xz[i] = t_0_xyyy_x_xxz[i] - rcd_x[i] * t_0_xyyy_x_xz[i];

        t_0_xyyy_xx_xy[i] = t_0_xyyy_x_xxy[i] - rcd_x[i] * t_0_xyyy_x_xy[i];

        t_0_xyyy_xx_xx[i] = t_0_xyyy_x_xxx[i] - rcd_x[i] * t_0_xyyy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxzz_x_xx, t_0_xxzz_x_xxx, t_0_xxzz_x_xxy,\
                           t_0_xxzz_x_xxz, t_0_xxzz_x_xy, t_0_xxzz_x_xyy, t_0_xxzz_x_xyz,\
                           t_0_xxzz_x_xz, t_0_xxzz_x_xzz, t_0_xxzz_x_yy, t_0_xxzz_x_yyy,\
                           t_0_xxzz_x_yyz, t_0_xxzz_x_yz, t_0_xxzz_x_yzz, t_0_xxzz_x_zz,\
                           t_0_xxzz_x_zzz, t_0_xxzz_xx_xx, t_0_xxzz_xx_xy, t_0_xxzz_xx_xz,\
                           t_0_xxzz_xx_yy, t_0_xxzz_xx_yz, t_0_xxzz_xx_zz, t_0_xxzz_xy_xx,\
                           t_0_xxzz_xy_xy, t_0_xxzz_xy_xz, t_0_xxzz_xy_yy, t_0_xxzz_xy_yz,\
                           t_0_xxzz_xy_zz, t_0_xxzz_xz_xx, t_0_xxzz_xz_xy, t_0_xxzz_xz_xz,\
                           t_0_xxzz_xz_yy, t_0_xxzz_xz_yz, t_0_xxzz_xz_zz, t_0_xxzz_y_xx,\
                           t_0_xxzz_y_xxy, t_0_xxzz_y_xxz, t_0_xxzz_y_xy, t_0_xxzz_y_xyy,\
                           t_0_xxzz_y_xyz, t_0_xxzz_y_xz, t_0_xxzz_y_xzz, t_0_xxzz_y_yy,\
                           t_0_xxzz_y_yyy, t_0_xxzz_y_yyz, t_0_xxzz_y_yz, t_0_xxzz_y_yzz,\
                           t_0_xxzz_y_zz, t_0_xxzz_y_zzz, t_0_xxzz_yy_xx, t_0_xxzz_yy_xy,\
                           t_0_xxzz_yy_xz, t_0_xxzz_yy_yy, t_0_xxzz_yy_yz, t_0_xxzz_yy_zz,\
                           t_0_xxzz_yz_xx, t_0_xxzz_yz_xy, t_0_xxzz_yz_xz, t_0_xxzz_yz_yy,\
                           t_0_xxzz_yz_yz, t_0_xxzz_yz_zz, t_0_xxzz_z_xx, t_0_xxzz_z_xxz,\
                           t_0_xxzz_z_xy, t_0_xxzz_z_xyz, t_0_xxzz_z_xz, t_0_xxzz_z_xzz,\
                           t_0_xxzz_z_yy, t_0_xxzz_z_yyz, t_0_xxzz_z_yz, t_0_xxzz_z_yzz,\
                           t_0_xxzz_z_zz, t_0_xxzz_z_zzz, t_0_xxzz_zz_xx, t_0_xxzz_zz_xy,\
                           t_0_xxzz_zz_xz, t_0_xxzz_zz_yy, t_0_xxzz_zz_yz, t_0_xxzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxzz_zz_zz[i] = t_0_xxzz_z_zzz[i] - rcd_z[i] * t_0_xxzz_z_zz[i];

        t_0_xxzz_zz_yz[i] = t_0_xxzz_z_yzz[i] - rcd_z[i] * t_0_xxzz_z_yz[i];

        t_0_xxzz_zz_yy[i] = t_0_xxzz_z_yyz[i] - rcd_z[i] * t_0_xxzz_z_yy[i];

        t_0_xxzz_zz_xz[i] = t_0_xxzz_z_xzz[i] - rcd_z[i] * t_0_xxzz_z_xz[i];

        t_0_xxzz_zz_xy[i] = t_0_xxzz_z_xyz[i] - rcd_z[i] * t_0_xxzz_z_xy[i];

        t_0_xxzz_zz_xx[i] = t_0_xxzz_z_xxz[i] - rcd_z[i] * t_0_xxzz_z_xx[i];

        t_0_xxzz_yz_zz[i] = t_0_xxzz_y_zzz[i] - rcd_z[i] * t_0_xxzz_y_zz[i];

        t_0_xxzz_yz_yz[i] = t_0_xxzz_y_yzz[i] - rcd_z[i] * t_0_xxzz_y_yz[i];

        t_0_xxzz_yz_yy[i] = t_0_xxzz_y_yyz[i] - rcd_z[i] * t_0_xxzz_y_yy[i];

        t_0_xxzz_yz_xz[i] = t_0_xxzz_y_xzz[i] - rcd_z[i] * t_0_xxzz_y_xz[i];

        t_0_xxzz_yz_xy[i] = t_0_xxzz_y_xyz[i] - rcd_z[i] * t_0_xxzz_y_xy[i];

        t_0_xxzz_yz_xx[i] = t_0_xxzz_y_xxz[i] - rcd_z[i] * t_0_xxzz_y_xx[i];

        t_0_xxzz_yy_zz[i] = t_0_xxzz_y_yzz[i] - rcd_y[i] * t_0_xxzz_y_zz[i];

        t_0_xxzz_yy_yz[i] = t_0_xxzz_y_yyz[i] - rcd_y[i] * t_0_xxzz_y_yz[i];

        t_0_xxzz_yy_yy[i] = t_0_xxzz_y_yyy[i] - rcd_y[i] * t_0_xxzz_y_yy[i];

        t_0_xxzz_yy_xz[i] = t_0_xxzz_y_xyz[i] - rcd_y[i] * t_0_xxzz_y_xz[i];

        t_0_xxzz_yy_xy[i] = t_0_xxzz_y_xyy[i] - rcd_y[i] * t_0_xxzz_y_xy[i];

        t_0_xxzz_yy_xx[i] = t_0_xxzz_y_xxy[i] - rcd_y[i] * t_0_xxzz_y_xx[i];

        t_0_xxzz_xz_zz[i] = t_0_xxzz_x_zzz[i] - rcd_z[i] * t_0_xxzz_x_zz[i];

        t_0_xxzz_xz_yz[i] = t_0_xxzz_x_yzz[i] - rcd_z[i] * t_0_xxzz_x_yz[i];

        t_0_xxzz_xz_yy[i] = t_0_xxzz_x_yyz[i] - rcd_z[i] * t_0_xxzz_x_yy[i];

        t_0_xxzz_xz_xz[i] = t_0_xxzz_x_xzz[i] - rcd_z[i] * t_0_xxzz_x_xz[i];

        t_0_xxzz_xz_xy[i] = t_0_xxzz_x_xyz[i] - rcd_z[i] * t_0_xxzz_x_xy[i];

        t_0_xxzz_xz_xx[i] = t_0_xxzz_x_xxz[i] - rcd_z[i] * t_0_xxzz_x_xx[i];

        t_0_xxzz_xy_zz[i] = t_0_xxzz_x_yzz[i] - rcd_y[i] * t_0_xxzz_x_zz[i];

        t_0_xxzz_xy_yz[i] = t_0_xxzz_x_yyz[i] - rcd_y[i] * t_0_xxzz_x_yz[i];

        t_0_xxzz_xy_yy[i] = t_0_xxzz_x_yyy[i] - rcd_y[i] * t_0_xxzz_x_yy[i];

        t_0_xxzz_xy_xz[i] = t_0_xxzz_x_xyz[i] - rcd_y[i] * t_0_xxzz_x_xz[i];

        t_0_xxzz_xy_xy[i] = t_0_xxzz_x_xyy[i] - rcd_y[i] * t_0_xxzz_x_xy[i];

        t_0_xxzz_xy_xx[i] = t_0_xxzz_x_xxy[i] - rcd_y[i] * t_0_xxzz_x_xx[i];

        t_0_xxzz_xx_zz[i] = t_0_xxzz_x_xzz[i] - rcd_x[i] * t_0_xxzz_x_zz[i];

        t_0_xxzz_xx_yz[i] = t_0_xxzz_x_xyz[i] - rcd_x[i] * t_0_xxzz_x_yz[i];

        t_0_xxzz_xx_yy[i] = t_0_xxzz_x_xyy[i] - rcd_x[i] * t_0_xxzz_x_yy[i];

        t_0_xxzz_xx_xz[i] = t_0_xxzz_x_xxz[i] - rcd_x[i] * t_0_xxzz_x_xz[i];

        t_0_xxzz_xx_xy[i] = t_0_xxzz_x_xxy[i] - rcd_x[i] * t_0_xxzz_x_xy[i];

        t_0_xxzz_xx_xx[i] = t_0_xxzz_x_xxx[i] - rcd_x[i] * t_0_xxzz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyz_x_xx, t_0_xxyz_x_xxx, t_0_xxyz_x_xxy,\
                           t_0_xxyz_x_xxz, t_0_xxyz_x_xy, t_0_xxyz_x_xyy, t_0_xxyz_x_xyz,\
                           t_0_xxyz_x_xz, t_0_xxyz_x_xzz, t_0_xxyz_x_yy, t_0_xxyz_x_yyy,\
                           t_0_xxyz_x_yyz, t_0_xxyz_x_yz, t_0_xxyz_x_yzz, t_0_xxyz_x_zz,\
                           t_0_xxyz_x_zzz, t_0_xxyz_xx_xx, t_0_xxyz_xx_xy, t_0_xxyz_xx_xz,\
                           t_0_xxyz_xx_yy, t_0_xxyz_xx_yz, t_0_xxyz_xx_zz, t_0_xxyz_xy_xx,\
                           t_0_xxyz_xy_xy, t_0_xxyz_xy_xz, t_0_xxyz_xy_yy, t_0_xxyz_xy_yz,\
                           t_0_xxyz_xy_zz, t_0_xxyz_xz_xx, t_0_xxyz_xz_xy, t_0_xxyz_xz_xz,\
                           t_0_xxyz_xz_yy, t_0_xxyz_xz_yz, t_0_xxyz_xz_zz, t_0_xxyz_y_xx,\
                           t_0_xxyz_y_xxy, t_0_xxyz_y_xxz, t_0_xxyz_y_xy, t_0_xxyz_y_xyy,\
                           t_0_xxyz_y_xyz, t_0_xxyz_y_xz, t_0_xxyz_y_xzz, t_0_xxyz_y_yy,\
                           t_0_xxyz_y_yyy, t_0_xxyz_y_yyz, t_0_xxyz_y_yz, t_0_xxyz_y_yzz,\
                           t_0_xxyz_y_zz, t_0_xxyz_y_zzz, t_0_xxyz_yy_xx, t_0_xxyz_yy_xy,\
                           t_0_xxyz_yy_xz, t_0_xxyz_yy_yy, t_0_xxyz_yy_yz, t_0_xxyz_yy_zz,\
                           t_0_xxyz_yz_xx, t_0_xxyz_yz_xy, t_0_xxyz_yz_xz, t_0_xxyz_yz_yy,\
                           t_0_xxyz_yz_yz, t_0_xxyz_yz_zz, t_0_xxyz_z_xx, t_0_xxyz_z_xxz,\
                           t_0_xxyz_z_xy, t_0_xxyz_z_xyz, t_0_xxyz_z_xz, t_0_xxyz_z_xzz,\
                           t_0_xxyz_z_yy, t_0_xxyz_z_yyz, t_0_xxyz_z_yz, t_0_xxyz_z_yzz,\
                           t_0_xxyz_z_zz, t_0_xxyz_z_zzz, t_0_xxyz_zz_xx, t_0_xxyz_zz_xy,\
                           t_0_xxyz_zz_xz, t_0_xxyz_zz_yy, t_0_xxyz_zz_yz, t_0_xxyz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxyz_zz_zz[i] = t_0_xxyz_z_zzz[i] - rcd_z[i] * t_0_xxyz_z_zz[i];

        t_0_xxyz_zz_yz[i] = t_0_xxyz_z_yzz[i] - rcd_z[i] * t_0_xxyz_z_yz[i];

        t_0_xxyz_zz_yy[i] = t_0_xxyz_z_yyz[i] - rcd_z[i] * t_0_xxyz_z_yy[i];

        t_0_xxyz_zz_xz[i] = t_0_xxyz_z_xzz[i] - rcd_z[i] * t_0_xxyz_z_xz[i];

        t_0_xxyz_zz_xy[i] = t_0_xxyz_z_xyz[i] - rcd_z[i] * t_0_xxyz_z_xy[i];

        t_0_xxyz_zz_xx[i] = t_0_xxyz_z_xxz[i] - rcd_z[i] * t_0_xxyz_z_xx[i];

        t_0_xxyz_yz_zz[i] = t_0_xxyz_y_zzz[i] - rcd_z[i] * t_0_xxyz_y_zz[i];

        t_0_xxyz_yz_yz[i] = t_0_xxyz_y_yzz[i] - rcd_z[i] * t_0_xxyz_y_yz[i];

        t_0_xxyz_yz_yy[i] = t_0_xxyz_y_yyz[i] - rcd_z[i] * t_0_xxyz_y_yy[i];

        t_0_xxyz_yz_xz[i] = t_0_xxyz_y_xzz[i] - rcd_z[i] * t_0_xxyz_y_xz[i];

        t_0_xxyz_yz_xy[i] = t_0_xxyz_y_xyz[i] - rcd_z[i] * t_0_xxyz_y_xy[i];

        t_0_xxyz_yz_xx[i] = t_0_xxyz_y_xxz[i] - rcd_z[i] * t_0_xxyz_y_xx[i];

        t_0_xxyz_yy_zz[i] = t_0_xxyz_y_yzz[i] - rcd_y[i] * t_0_xxyz_y_zz[i];

        t_0_xxyz_yy_yz[i] = t_0_xxyz_y_yyz[i] - rcd_y[i] * t_0_xxyz_y_yz[i];

        t_0_xxyz_yy_yy[i] = t_0_xxyz_y_yyy[i] - rcd_y[i] * t_0_xxyz_y_yy[i];

        t_0_xxyz_yy_xz[i] = t_0_xxyz_y_xyz[i] - rcd_y[i] * t_0_xxyz_y_xz[i];

        t_0_xxyz_yy_xy[i] = t_0_xxyz_y_xyy[i] - rcd_y[i] * t_0_xxyz_y_xy[i];

        t_0_xxyz_yy_xx[i] = t_0_xxyz_y_xxy[i] - rcd_y[i] * t_0_xxyz_y_xx[i];

        t_0_xxyz_xz_zz[i] = t_0_xxyz_x_zzz[i] - rcd_z[i] * t_0_xxyz_x_zz[i];

        t_0_xxyz_xz_yz[i] = t_0_xxyz_x_yzz[i] - rcd_z[i] * t_0_xxyz_x_yz[i];

        t_0_xxyz_xz_yy[i] = t_0_xxyz_x_yyz[i] - rcd_z[i] * t_0_xxyz_x_yy[i];

        t_0_xxyz_xz_xz[i] = t_0_xxyz_x_xzz[i] - rcd_z[i] * t_0_xxyz_x_xz[i];

        t_0_xxyz_xz_xy[i] = t_0_xxyz_x_xyz[i] - rcd_z[i] * t_0_xxyz_x_xy[i];

        t_0_xxyz_xz_xx[i] = t_0_xxyz_x_xxz[i] - rcd_z[i] * t_0_xxyz_x_xx[i];

        t_0_xxyz_xy_zz[i] = t_0_xxyz_x_yzz[i] - rcd_y[i] * t_0_xxyz_x_zz[i];

        t_0_xxyz_xy_yz[i] = t_0_xxyz_x_yyz[i] - rcd_y[i] * t_0_xxyz_x_yz[i];

        t_0_xxyz_xy_yy[i] = t_0_xxyz_x_yyy[i] - rcd_y[i] * t_0_xxyz_x_yy[i];

        t_0_xxyz_xy_xz[i] = t_0_xxyz_x_xyz[i] - rcd_y[i] * t_0_xxyz_x_xz[i];

        t_0_xxyz_xy_xy[i] = t_0_xxyz_x_xyy[i] - rcd_y[i] * t_0_xxyz_x_xy[i];

        t_0_xxyz_xy_xx[i] = t_0_xxyz_x_xxy[i] - rcd_y[i] * t_0_xxyz_x_xx[i];

        t_0_xxyz_xx_zz[i] = t_0_xxyz_x_xzz[i] - rcd_x[i] * t_0_xxyz_x_zz[i];

        t_0_xxyz_xx_yz[i] = t_0_xxyz_x_xyz[i] - rcd_x[i] * t_0_xxyz_x_yz[i];

        t_0_xxyz_xx_yy[i] = t_0_xxyz_x_xyy[i] - rcd_x[i] * t_0_xxyz_x_yy[i];

        t_0_xxyz_xx_xz[i] = t_0_xxyz_x_xxz[i] - rcd_x[i] * t_0_xxyz_x_xz[i];

        t_0_xxyz_xx_xy[i] = t_0_xxyz_x_xxy[i] - rcd_x[i] * t_0_xxyz_x_xy[i];

        t_0_xxyz_xx_xx[i] = t_0_xxyz_x_xxx[i] - rcd_x[i] * t_0_xxyz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyy_x_xx, t_0_xxyy_x_xxx, t_0_xxyy_x_xxy,\
                           t_0_xxyy_x_xxz, t_0_xxyy_x_xy, t_0_xxyy_x_xyy, t_0_xxyy_x_xyz,\
                           t_0_xxyy_x_xz, t_0_xxyy_x_xzz, t_0_xxyy_x_yy, t_0_xxyy_x_yyy,\
                           t_0_xxyy_x_yyz, t_0_xxyy_x_yz, t_0_xxyy_x_yzz, t_0_xxyy_x_zz,\
                           t_0_xxyy_x_zzz, t_0_xxyy_xx_xx, t_0_xxyy_xx_xy, t_0_xxyy_xx_xz,\
                           t_0_xxyy_xx_yy, t_0_xxyy_xx_yz, t_0_xxyy_xx_zz, t_0_xxyy_xy_xx,\
                           t_0_xxyy_xy_xy, t_0_xxyy_xy_xz, t_0_xxyy_xy_yy, t_0_xxyy_xy_yz,\
                           t_0_xxyy_xy_zz, t_0_xxyy_xz_xx, t_0_xxyy_xz_xy, t_0_xxyy_xz_xz,\
                           t_0_xxyy_xz_yy, t_0_xxyy_xz_yz, t_0_xxyy_xz_zz, t_0_xxyy_y_xx,\
                           t_0_xxyy_y_xxy, t_0_xxyy_y_xxz, t_0_xxyy_y_xy, t_0_xxyy_y_xyy,\
                           t_0_xxyy_y_xyz, t_0_xxyy_y_xz, t_0_xxyy_y_xzz, t_0_xxyy_y_yy,\
                           t_0_xxyy_y_yyy, t_0_xxyy_y_yyz, t_0_xxyy_y_yz, t_0_xxyy_y_yzz,\
                           t_0_xxyy_y_zz, t_0_xxyy_y_zzz, t_0_xxyy_yy_xx, t_0_xxyy_yy_xy,\
                           t_0_xxyy_yy_xz, t_0_xxyy_yy_yy, t_0_xxyy_yy_yz, t_0_xxyy_yy_zz,\
                           t_0_xxyy_yz_xx, t_0_xxyy_yz_xy, t_0_xxyy_yz_xz, t_0_xxyy_yz_yy,\
                           t_0_xxyy_yz_yz, t_0_xxyy_yz_zz, t_0_xxyy_z_xx, t_0_xxyy_z_xxz,\
                           t_0_xxyy_z_xy, t_0_xxyy_z_xyz, t_0_xxyy_z_xz, t_0_xxyy_z_xzz,\
                           t_0_xxyy_z_yy, t_0_xxyy_z_yyz, t_0_xxyy_z_yz, t_0_xxyy_z_yzz,\
                           t_0_xxyy_z_zz, t_0_xxyy_z_zzz, t_0_xxyy_zz_xx, t_0_xxyy_zz_xy,\
                           t_0_xxyy_zz_xz, t_0_xxyy_zz_yy, t_0_xxyy_zz_yz, t_0_xxyy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxyy_zz_zz[i] = t_0_xxyy_z_zzz[i] - rcd_z[i] * t_0_xxyy_z_zz[i];

        t_0_xxyy_zz_yz[i] = t_0_xxyy_z_yzz[i] - rcd_z[i] * t_0_xxyy_z_yz[i];

        t_0_xxyy_zz_yy[i] = t_0_xxyy_z_yyz[i] - rcd_z[i] * t_0_xxyy_z_yy[i];

        t_0_xxyy_zz_xz[i] = t_0_xxyy_z_xzz[i] - rcd_z[i] * t_0_xxyy_z_xz[i];

        t_0_xxyy_zz_xy[i] = t_0_xxyy_z_xyz[i] - rcd_z[i] * t_0_xxyy_z_xy[i];

        t_0_xxyy_zz_xx[i] = t_0_xxyy_z_xxz[i] - rcd_z[i] * t_0_xxyy_z_xx[i];

        t_0_xxyy_yz_zz[i] = t_0_xxyy_y_zzz[i] - rcd_z[i] * t_0_xxyy_y_zz[i];

        t_0_xxyy_yz_yz[i] = t_0_xxyy_y_yzz[i] - rcd_z[i] * t_0_xxyy_y_yz[i];

        t_0_xxyy_yz_yy[i] = t_0_xxyy_y_yyz[i] - rcd_z[i] * t_0_xxyy_y_yy[i];

        t_0_xxyy_yz_xz[i] = t_0_xxyy_y_xzz[i] - rcd_z[i] * t_0_xxyy_y_xz[i];

        t_0_xxyy_yz_xy[i] = t_0_xxyy_y_xyz[i] - rcd_z[i] * t_0_xxyy_y_xy[i];

        t_0_xxyy_yz_xx[i] = t_0_xxyy_y_xxz[i] - rcd_z[i] * t_0_xxyy_y_xx[i];

        t_0_xxyy_yy_zz[i] = t_0_xxyy_y_yzz[i] - rcd_y[i] * t_0_xxyy_y_zz[i];

        t_0_xxyy_yy_yz[i] = t_0_xxyy_y_yyz[i] - rcd_y[i] * t_0_xxyy_y_yz[i];

        t_0_xxyy_yy_yy[i] = t_0_xxyy_y_yyy[i] - rcd_y[i] * t_0_xxyy_y_yy[i];

        t_0_xxyy_yy_xz[i] = t_0_xxyy_y_xyz[i] - rcd_y[i] * t_0_xxyy_y_xz[i];

        t_0_xxyy_yy_xy[i] = t_0_xxyy_y_xyy[i] - rcd_y[i] * t_0_xxyy_y_xy[i];

        t_0_xxyy_yy_xx[i] = t_0_xxyy_y_xxy[i] - rcd_y[i] * t_0_xxyy_y_xx[i];

        t_0_xxyy_xz_zz[i] = t_0_xxyy_x_zzz[i] - rcd_z[i] * t_0_xxyy_x_zz[i];

        t_0_xxyy_xz_yz[i] = t_0_xxyy_x_yzz[i] - rcd_z[i] * t_0_xxyy_x_yz[i];

        t_0_xxyy_xz_yy[i] = t_0_xxyy_x_yyz[i] - rcd_z[i] * t_0_xxyy_x_yy[i];

        t_0_xxyy_xz_xz[i] = t_0_xxyy_x_xzz[i] - rcd_z[i] * t_0_xxyy_x_xz[i];

        t_0_xxyy_xz_xy[i] = t_0_xxyy_x_xyz[i] - rcd_z[i] * t_0_xxyy_x_xy[i];

        t_0_xxyy_xz_xx[i] = t_0_xxyy_x_xxz[i] - rcd_z[i] * t_0_xxyy_x_xx[i];

        t_0_xxyy_xy_zz[i] = t_0_xxyy_x_yzz[i] - rcd_y[i] * t_0_xxyy_x_zz[i];

        t_0_xxyy_xy_yz[i] = t_0_xxyy_x_yyz[i] - rcd_y[i] * t_0_xxyy_x_yz[i];

        t_0_xxyy_xy_yy[i] = t_0_xxyy_x_yyy[i] - rcd_y[i] * t_0_xxyy_x_yy[i];

        t_0_xxyy_xy_xz[i] = t_0_xxyy_x_xyz[i] - rcd_y[i] * t_0_xxyy_x_xz[i];

        t_0_xxyy_xy_xy[i] = t_0_xxyy_x_xyy[i] - rcd_y[i] * t_0_xxyy_x_xy[i];

        t_0_xxyy_xy_xx[i] = t_0_xxyy_x_xxy[i] - rcd_y[i] * t_0_xxyy_x_xx[i];

        t_0_xxyy_xx_zz[i] = t_0_xxyy_x_xzz[i] - rcd_x[i] * t_0_xxyy_x_zz[i];

        t_0_xxyy_xx_yz[i] = t_0_xxyy_x_xyz[i] - rcd_x[i] * t_0_xxyy_x_yz[i];

        t_0_xxyy_xx_yy[i] = t_0_xxyy_x_xyy[i] - rcd_x[i] * t_0_xxyy_x_yy[i];

        t_0_xxyy_xx_xz[i] = t_0_xxyy_x_xxz[i] - rcd_x[i] * t_0_xxyy_x_xz[i];

        t_0_xxyy_xx_xy[i] = t_0_xxyy_x_xxy[i] - rcd_x[i] * t_0_xxyy_x_xy[i];

        t_0_xxyy_xx_xx[i] = t_0_xxyy_x_xxx[i] - rcd_x[i] * t_0_xxyy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxz_x_xx, t_0_xxxz_x_xxx, t_0_xxxz_x_xxy,\
                           t_0_xxxz_x_xxz, t_0_xxxz_x_xy, t_0_xxxz_x_xyy, t_0_xxxz_x_xyz,\
                           t_0_xxxz_x_xz, t_0_xxxz_x_xzz, t_0_xxxz_x_yy, t_0_xxxz_x_yyy,\
                           t_0_xxxz_x_yyz, t_0_xxxz_x_yz, t_0_xxxz_x_yzz, t_0_xxxz_x_zz,\
                           t_0_xxxz_x_zzz, t_0_xxxz_xx_xx, t_0_xxxz_xx_xy, t_0_xxxz_xx_xz,\
                           t_0_xxxz_xx_yy, t_0_xxxz_xx_yz, t_0_xxxz_xx_zz, t_0_xxxz_xy_xx,\
                           t_0_xxxz_xy_xy, t_0_xxxz_xy_xz, t_0_xxxz_xy_yy, t_0_xxxz_xy_yz,\
                           t_0_xxxz_xy_zz, t_0_xxxz_xz_xx, t_0_xxxz_xz_xy, t_0_xxxz_xz_xz,\
                           t_0_xxxz_xz_yy, t_0_xxxz_xz_yz, t_0_xxxz_xz_zz, t_0_xxxz_y_xx,\
                           t_0_xxxz_y_xxy, t_0_xxxz_y_xxz, t_0_xxxz_y_xy, t_0_xxxz_y_xyy,\
                           t_0_xxxz_y_xyz, t_0_xxxz_y_xz, t_0_xxxz_y_xzz, t_0_xxxz_y_yy,\
                           t_0_xxxz_y_yyy, t_0_xxxz_y_yyz, t_0_xxxz_y_yz, t_0_xxxz_y_yzz,\
                           t_0_xxxz_y_zz, t_0_xxxz_y_zzz, t_0_xxxz_yy_xx, t_0_xxxz_yy_xy,\
                           t_0_xxxz_yy_xz, t_0_xxxz_yy_yy, t_0_xxxz_yy_yz, t_0_xxxz_yy_zz,\
                           t_0_xxxz_yz_xx, t_0_xxxz_yz_xy, t_0_xxxz_yz_xz, t_0_xxxz_yz_yy,\
                           t_0_xxxz_yz_yz, t_0_xxxz_yz_zz, t_0_xxxz_z_xx, t_0_xxxz_z_xxz,\
                           t_0_xxxz_z_xy, t_0_xxxz_z_xyz, t_0_xxxz_z_xz, t_0_xxxz_z_xzz,\
                           t_0_xxxz_z_yy, t_0_xxxz_z_yyz, t_0_xxxz_z_yz, t_0_xxxz_z_yzz,\
                           t_0_xxxz_z_zz, t_0_xxxz_z_zzz, t_0_xxxz_zz_xx, t_0_xxxz_zz_xy,\
                           t_0_xxxz_zz_xz, t_0_xxxz_zz_yy, t_0_xxxz_zz_yz, t_0_xxxz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxz_zz_zz[i] = t_0_xxxz_z_zzz[i] - rcd_z[i] * t_0_xxxz_z_zz[i];

        t_0_xxxz_zz_yz[i] = t_0_xxxz_z_yzz[i] - rcd_z[i] * t_0_xxxz_z_yz[i];

        t_0_xxxz_zz_yy[i] = t_0_xxxz_z_yyz[i] - rcd_z[i] * t_0_xxxz_z_yy[i];

        t_0_xxxz_zz_xz[i] = t_0_xxxz_z_xzz[i] - rcd_z[i] * t_0_xxxz_z_xz[i];

        t_0_xxxz_zz_xy[i] = t_0_xxxz_z_xyz[i] - rcd_z[i] * t_0_xxxz_z_xy[i];

        t_0_xxxz_zz_xx[i] = t_0_xxxz_z_xxz[i] - rcd_z[i] * t_0_xxxz_z_xx[i];

        t_0_xxxz_yz_zz[i] = t_0_xxxz_y_zzz[i] - rcd_z[i] * t_0_xxxz_y_zz[i];

        t_0_xxxz_yz_yz[i] = t_0_xxxz_y_yzz[i] - rcd_z[i] * t_0_xxxz_y_yz[i];

        t_0_xxxz_yz_yy[i] = t_0_xxxz_y_yyz[i] - rcd_z[i] * t_0_xxxz_y_yy[i];

        t_0_xxxz_yz_xz[i] = t_0_xxxz_y_xzz[i] - rcd_z[i] * t_0_xxxz_y_xz[i];

        t_0_xxxz_yz_xy[i] = t_0_xxxz_y_xyz[i] - rcd_z[i] * t_0_xxxz_y_xy[i];

        t_0_xxxz_yz_xx[i] = t_0_xxxz_y_xxz[i] - rcd_z[i] * t_0_xxxz_y_xx[i];

        t_0_xxxz_yy_zz[i] = t_0_xxxz_y_yzz[i] - rcd_y[i] * t_0_xxxz_y_zz[i];

        t_0_xxxz_yy_yz[i] = t_0_xxxz_y_yyz[i] - rcd_y[i] * t_0_xxxz_y_yz[i];

        t_0_xxxz_yy_yy[i] = t_0_xxxz_y_yyy[i] - rcd_y[i] * t_0_xxxz_y_yy[i];

        t_0_xxxz_yy_xz[i] = t_0_xxxz_y_xyz[i] - rcd_y[i] * t_0_xxxz_y_xz[i];

        t_0_xxxz_yy_xy[i] = t_0_xxxz_y_xyy[i] - rcd_y[i] * t_0_xxxz_y_xy[i];

        t_0_xxxz_yy_xx[i] = t_0_xxxz_y_xxy[i] - rcd_y[i] * t_0_xxxz_y_xx[i];

        t_0_xxxz_xz_zz[i] = t_0_xxxz_x_zzz[i] - rcd_z[i] * t_0_xxxz_x_zz[i];

        t_0_xxxz_xz_yz[i] = t_0_xxxz_x_yzz[i] - rcd_z[i] * t_0_xxxz_x_yz[i];

        t_0_xxxz_xz_yy[i] = t_0_xxxz_x_yyz[i] - rcd_z[i] * t_0_xxxz_x_yy[i];

        t_0_xxxz_xz_xz[i] = t_0_xxxz_x_xzz[i] - rcd_z[i] * t_0_xxxz_x_xz[i];

        t_0_xxxz_xz_xy[i] = t_0_xxxz_x_xyz[i] - rcd_z[i] * t_0_xxxz_x_xy[i];

        t_0_xxxz_xz_xx[i] = t_0_xxxz_x_xxz[i] - rcd_z[i] * t_0_xxxz_x_xx[i];

        t_0_xxxz_xy_zz[i] = t_0_xxxz_x_yzz[i] - rcd_y[i] * t_0_xxxz_x_zz[i];

        t_0_xxxz_xy_yz[i] = t_0_xxxz_x_yyz[i] - rcd_y[i] * t_0_xxxz_x_yz[i];

        t_0_xxxz_xy_yy[i] = t_0_xxxz_x_yyy[i] - rcd_y[i] * t_0_xxxz_x_yy[i];

        t_0_xxxz_xy_xz[i] = t_0_xxxz_x_xyz[i] - rcd_y[i] * t_0_xxxz_x_xz[i];

        t_0_xxxz_xy_xy[i] = t_0_xxxz_x_xyy[i] - rcd_y[i] * t_0_xxxz_x_xy[i];

        t_0_xxxz_xy_xx[i] = t_0_xxxz_x_xxy[i] - rcd_y[i] * t_0_xxxz_x_xx[i];

        t_0_xxxz_xx_zz[i] = t_0_xxxz_x_xzz[i] - rcd_x[i] * t_0_xxxz_x_zz[i];

        t_0_xxxz_xx_yz[i] = t_0_xxxz_x_xyz[i] - rcd_x[i] * t_0_xxxz_x_yz[i];

        t_0_xxxz_xx_yy[i] = t_0_xxxz_x_xyy[i] - rcd_x[i] * t_0_xxxz_x_yy[i];

        t_0_xxxz_xx_xz[i] = t_0_xxxz_x_xxz[i] - rcd_x[i] * t_0_xxxz_x_xz[i];

        t_0_xxxz_xx_xy[i] = t_0_xxxz_x_xxy[i] - rcd_x[i] * t_0_xxxz_x_xy[i];

        t_0_xxxz_xx_xx[i] = t_0_xxxz_x_xxx[i] - rcd_x[i] * t_0_xxxz_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxy_x_xx, t_0_xxxy_x_xxx, t_0_xxxy_x_xxy,\
                           t_0_xxxy_x_xxz, t_0_xxxy_x_xy, t_0_xxxy_x_xyy, t_0_xxxy_x_xyz,\
                           t_0_xxxy_x_xz, t_0_xxxy_x_xzz, t_0_xxxy_x_yy, t_0_xxxy_x_yyy,\
                           t_0_xxxy_x_yyz, t_0_xxxy_x_yz, t_0_xxxy_x_yzz, t_0_xxxy_x_zz,\
                           t_0_xxxy_x_zzz, t_0_xxxy_xx_xx, t_0_xxxy_xx_xy, t_0_xxxy_xx_xz,\
                           t_0_xxxy_xx_yy, t_0_xxxy_xx_yz, t_0_xxxy_xx_zz, t_0_xxxy_xy_xx,\
                           t_0_xxxy_xy_xy, t_0_xxxy_xy_xz, t_0_xxxy_xy_yy, t_0_xxxy_xy_yz,\
                           t_0_xxxy_xy_zz, t_0_xxxy_xz_xx, t_0_xxxy_xz_xy, t_0_xxxy_xz_xz,\
                           t_0_xxxy_xz_yy, t_0_xxxy_xz_yz, t_0_xxxy_xz_zz, t_0_xxxy_y_xx,\
                           t_0_xxxy_y_xxy, t_0_xxxy_y_xxz, t_0_xxxy_y_xy, t_0_xxxy_y_xyy,\
                           t_0_xxxy_y_xyz, t_0_xxxy_y_xz, t_0_xxxy_y_xzz, t_0_xxxy_y_yy,\
                           t_0_xxxy_y_yyy, t_0_xxxy_y_yyz, t_0_xxxy_y_yz, t_0_xxxy_y_yzz,\
                           t_0_xxxy_y_zz, t_0_xxxy_y_zzz, t_0_xxxy_yy_xx, t_0_xxxy_yy_xy,\
                           t_0_xxxy_yy_xz, t_0_xxxy_yy_yy, t_0_xxxy_yy_yz, t_0_xxxy_yy_zz,\
                           t_0_xxxy_yz_xx, t_0_xxxy_yz_xy, t_0_xxxy_yz_xz, t_0_xxxy_yz_yy,\
                           t_0_xxxy_yz_yz, t_0_xxxy_yz_zz, t_0_xxxy_z_xx, t_0_xxxy_z_xxz,\
                           t_0_xxxy_z_xy, t_0_xxxy_z_xyz, t_0_xxxy_z_xz, t_0_xxxy_z_xzz,\
                           t_0_xxxy_z_yy, t_0_xxxy_z_yyz, t_0_xxxy_z_yz, t_0_xxxy_z_yzz,\
                           t_0_xxxy_z_zz, t_0_xxxy_z_zzz, t_0_xxxy_zz_xx, t_0_xxxy_zz_xy,\
                           t_0_xxxy_zz_xz, t_0_xxxy_zz_yy, t_0_xxxy_zz_yz, t_0_xxxy_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxy_zz_zz[i] = t_0_xxxy_z_zzz[i] - rcd_z[i] * t_0_xxxy_z_zz[i];

        t_0_xxxy_zz_yz[i] = t_0_xxxy_z_yzz[i] - rcd_z[i] * t_0_xxxy_z_yz[i];

        t_0_xxxy_zz_yy[i] = t_0_xxxy_z_yyz[i] - rcd_z[i] * t_0_xxxy_z_yy[i];

        t_0_xxxy_zz_xz[i] = t_0_xxxy_z_xzz[i] - rcd_z[i] * t_0_xxxy_z_xz[i];

        t_0_xxxy_zz_xy[i] = t_0_xxxy_z_xyz[i] - rcd_z[i] * t_0_xxxy_z_xy[i];

        t_0_xxxy_zz_xx[i] = t_0_xxxy_z_xxz[i] - rcd_z[i] * t_0_xxxy_z_xx[i];

        t_0_xxxy_yz_zz[i] = t_0_xxxy_y_zzz[i] - rcd_z[i] * t_0_xxxy_y_zz[i];

        t_0_xxxy_yz_yz[i] = t_0_xxxy_y_yzz[i] - rcd_z[i] * t_0_xxxy_y_yz[i];

        t_0_xxxy_yz_yy[i] = t_0_xxxy_y_yyz[i] - rcd_z[i] * t_0_xxxy_y_yy[i];

        t_0_xxxy_yz_xz[i] = t_0_xxxy_y_xzz[i] - rcd_z[i] * t_0_xxxy_y_xz[i];

        t_0_xxxy_yz_xy[i] = t_0_xxxy_y_xyz[i] - rcd_z[i] * t_0_xxxy_y_xy[i];

        t_0_xxxy_yz_xx[i] = t_0_xxxy_y_xxz[i] - rcd_z[i] * t_0_xxxy_y_xx[i];

        t_0_xxxy_yy_zz[i] = t_0_xxxy_y_yzz[i] - rcd_y[i] * t_0_xxxy_y_zz[i];

        t_0_xxxy_yy_yz[i] = t_0_xxxy_y_yyz[i] - rcd_y[i] * t_0_xxxy_y_yz[i];

        t_0_xxxy_yy_yy[i] = t_0_xxxy_y_yyy[i] - rcd_y[i] * t_0_xxxy_y_yy[i];

        t_0_xxxy_yy_xz[i] = t_0_xxxy_y_xyz[i] - rcd_y[i] * t_0_xxxy_y_xz[i];

        t_0_xxxy_yy_xy[i] = t_0_xxxy_y_xyy[i] - rcd_y[i] * t_0_xxxy_y_xy[i];

        t_0_xxxy_yy_xx[i] = t_0_xxxy_y_xxy[i] - rcd_y[i] * t_0_xxxy_y_xx[i];

        t_0_xxxy_xz_zz[i] = t_0_xxxy_x_zzz[i] - rcd_z[i] * t_0_xxxy_x_zz[i];

        t_0_xxxy_xz_yz[i] = t_0_xxxy_x_yzz[i] - rcd_z[i] * t_0_xxxy_x_yz[i];

        t_0_xxxy_xz_yy[i] = t_0_xxxy_x_yyz[i] - rcd_z[i] * t_0_xxxy_x_yy[i];

        t_0_xxxy_xz_xz[i] = t_0_xxxy_x_xzz[i] - rcd_z[i] * t_0_xxxy_x_xz[i];

        t_0_xxxy_xz_xy[i] = t_0_xxxy_x_xyz[i] - rcd_z[i] * t_0_xxxy_x_xy[i];

        t_0_xxxy_xz_xx[i] = t_0_xxxy_x_xxz[i] - rcd_z[i] * t_0_xxxy_x_xx[i];

        t_0_xxxy_xy_zz[i] = t_0_xxxy_x_yzz[i] - rcd_y[i] * t_0_xxxy_x_zz[i];

        t_0_xxxy_xy_yz[i] = t_0_xxxy_x_yyz[i] - rcd_y[i] * t_0_xxxy_x_yz[i];

        t_0_xxxy_xy_yy[i] = t_0_xxxy_x_yyy[i] - rcd_y[i] * t_0_xxxy_x_yy[i];

        t_0_xxxy_xy_xz[i] = t_0_xxxy_x_xyz[i] - rcd_y[i] * t_0_xxxy_x_xz[i];

        t_0_xxxy_xy_xy[i] = t_0_xxxy_x_xyy[i] - rcd_y[i] * t_0_xxxy_x_xy[i];

        t_0_xxxy_xy_xx[i] = t_0_xxxy_x_xxy[i] - rcd_y[i] * t_0_xxxy_x_xx[i];

        t_0_xxxy_xx_zz[i] = t_0_xxxy_x_xzz[i] - rcd_x[i] * t_0_xxxy_x_zz[i];

        t_0_xxxy_xx_yz[i] = t_0_xxxy_x_xyz[i] - rcd_x[i] * t_0_xxxy_x_yz[i];

        t_0_xxxy_xx_yy[i] = t_0_xxxy_x_xyy[i] - rcd_x[i] * t_0_xxxy_x_yy[i];

        t_0_xxxy_xx_xz[i] = t_0_xxxy_x_xxz[i] - rcd_x[i] * t_0_xxxy_x_xz[i];

        t_0_xxxy_xx_xy[i] = t_0_xxxy_x_xxy[i] - rcd_x[i] * t_0_xxxy_x_xy[i];

        t_0_xxxy_xx_xx[i] = t_0_xxxy_x_xxx[i] - rcd_x[i] * t_0_xxxy_x_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_x_xx, t_0_xxxx_x_xxx, t_0_xxxx_x_xxy,\
                           t_0_xxxx_x_xxz, t_0_xxxx_x_xy, t_0_xxxx_x_xyy, t_0_xxxx_x_xyz,\
                           t_0_xxxx_x_xz, t_0_xxxx_x_xzz, t_0_xxxx_x_yy, t_0_xxxx_x_yyy,\
                           t_0_xxxx_x_yyz, t_0_xxxx_x_yz, t_0_xxxx_x_yzz, t_0_xxxx_x_zz,\
                           t_0_xxxx_x_zzz, t_0_xxxx_xx_xx, t_0_xxxx_xx_xy, t_0_xxxx_xx_xz,\
                           t_0_xxxx_xx_yy, t_0_xxxx_xx_yz, t_0_xxxx_xx_zz, t_0_xxxx_xy_xx,\
                           t_0_xxxx_xy_xy, t_0_xxxx_xy_xz, t_0_xxxx_xy_yy, t_0_xxxx_xy_yz,\
                           t_0_xxxx_xy_zz, t_0_xxxx_xz_xx, t_0_xxxx_xz_xy, t_0_xxxx_xz_xz,\
                           t_0_xxxx_xz_yy, t_0_xxxx_xz_yz, t_0_xxxx_xz_zz, t_0_xxxx_y_xx,\
                           t_0_xxxx_y_xxy, t_0_xxxx_y_xxz, t_0_xxxx_y_xy, t_0_xxxx_y_xyy,\
                           t_0_xxxx_y_xyz, t_0_xxxx_y_xz, t_0_xxxx_y_xzz, t_0_xxxx_y_yy,\
                           t_0_xxxx_y_yyy, t_0_xxxx_y_yyz, t_0_xxxx_y_yz, t_0_xxxx_y_yzz,\
                           t_0_xxxx_y_zz, t_0_xxxx_y_zzz, t_0_xxxx_yy_xx, t_0_xxxx_yy_xy,\
                           t_0_xxxx_yy_xz, t_0_xxxx_yy_yy, t_0_xxxx_yy_yz, t_0_xxxx_yy_zz,\
                           t_0_xxxx_yz_xx, t_0_xxxx_yz_xy, t_0_xxxx_yz_xz, t_0_xxxx_yz_yy,\
                           t_0_xxxx_yz_yz, t_0_xxxx_yz_zz, t_0_xxxx_z_xx, t_0_xxxx_z_xxz,\
                           t_0_xxxx_z_xy, t_0_xxxx_z_xyz, t_0_xxxx_z_xz, t_0_xxxx_z_xzz,\
                           t_0_xxxx_z_yy, t_0_xxxx_z_yyz, t_0_xxxx_z_yz, t_0_xxxx_z_yzz,\
                           t_0_xxxx_z_zz, t_0_xxxx_z_zzz, t_0_xxxx_zz_xx, t_0_xxxx_zz_xy,\
                           t_0_xxxx_zz_xz, t_0_xxxx_zz_yy, t_0_xxxx_zz_yz, t_0_xxxx_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxx_zz_zz[i] = t_0_xxxx_z_zzz[i] - rcd_z[i] * t_0_xxxx_z_zz[i];

        t_0_xxxx_zz_yz[i] = t_0_xxxx_z_yzz[i] - rcd_z[i] * t_0_xxxx_z_yz[i];

        t_0_xxxx_zz_yy[i] = t_0_xxxx_z_yyz[i] - rcd_z[i] * t_0_xxxx_z_yy[i];

        t_0_xxxx_zz_xz[i] = t_0_xxxx_z_xzz[i] - rcd_z[i] * t_0_xxxx_z_xz[i];

        t_0_xxxx_zz_xy[i] = t_0_xxxx_z_xyz[i] - rcd_z[i] * t_0_xxxx_z_xy[i];

        t_0_xxxx_zz_xx[i] = t_0_xxxx_z_xxz[i] - rcd_z[i] * t_0_xxxx_z_xx[i];

        t_0_xxxx_yz_zz[i] = t_0_xxxx_y_zzz[i] - rcd_z[i] * t_0_xxxx_y_zz[i];

        t_0_xxxx_yz_yz[i] = t_0_xxxx_y_yzz[i] - rcd_z[i] * t_0_xxxx_y_yz[i];

        t_0_xxxx_yz_yy[i] = t_0_xxxx_y_yyz[i] - rcd_z[i] * t_0_xxxx_y_yy[i];

        t_0_xxxx_yz_xz[i] = t_0_xxxx_y_xzz[i] - rcd_z[i] * t_0_xxxx_y_xz[i];

        t_0_xxxx_yz_xy[i] = t_0_xxxx_y_xyz[i] - rcd_z[i] * t_0_xxxx_y_xy[i];

        t_0_xxxx_yz_xx[i] = t_0_xxxx_y_xxz[i] - rcd_z[i] * t_0_xxxx_y_xx[i];

        t_0_xxxx_yy_zz[i] = t_0_xxxx_y_yzz[i] - rcd_y[i] * t_0_xxxx_y_zz[i];

        t_0_xxxx_yy_yz[i] = t_0_xxxx_y_yyz[i] - rcd_y[i] * t_0_xxxx_y_yz[i];

        t_0_xxxx_yy_yy[i] = t_0_xxxx_y_yyy[i] - rcd_y[i] * t_0_xxxx_y_yy[i];

        t_0_xxxx_yy_xz[i] = t_0_xxxx_y_xyz[i] - rcd_y[i] * t_0_xxxx_y_xz[i];

        t_0_xxxx_yy_xy[i] = t_0_xxxx_y_xyy[i] - rcd_y[i] * t_0_xxxx_y_xy[i];

        t_0_xxxx_yy_xx[i] = t_0_xxxx_y_xxy[i] - rcd_y[i] * t_0_xxxx_y_xx[i];

        t_0_xxxx_xz_zz[i] = t_0_xxxx_x_zzz[i] - rcd_z[i] * t_0_xxxx_x_zz[i];

        t_0_xxxx_xz_yz[i] = t_0_xxxx_x_yzz[i] - rcd_z[i] * t_0_xxxx_x_yz[i];

        t_0_xxxx_xz_yy[i] = t_0_xxxx_x_yyz[i] - rcd_z[i] * t_0_xxxx_x_yy[i];

        t_0_xxxx_xz_xz[i] = t_0_xxxx_x_xzz[i] - rcd_z[i] * t_0_xxxx_x_xz[i];

        t_0_xxxx_xz_xy[i] = t_0_xxxx_x_xyz[i] - rcd_z[i] * t_0_xxxx_x_xy[i];

        t_0_xxxx_xz_xx[i] = t_0_xxxx_x_xxz[i] - rcd_z[i] * t_0_xxxx_x_xx[i];

        t_0_xxxx_xy_zz[i] = t_0_xxxx_x_yzz[i] - rcd_y[i] * t_0_xxxx_x_zz[i];

        t_0_xxxx_xy_yz[i] = t_0_xxxx_x_yyz[i] - rcd_y[i] * t_0_xxxx_x_yz[i];

        t_0_xxxx_xy_yy[i] = t_0_xxxx_x_yyy[i] - rcd_y[i] * t_0_xxxx_x_yy[i];

        t_0_xxxx_xy_xz[i] = t_0_xxxx_x_xyz[i] - rcd_y[i] * t_0_xxxx_x_xz[i];

        t_0_xxxx_xy_xy[i] = t_0_xxxx_x_xyy[i] - rcd_y[i] * t_0_xxxx_x_xy[i];

        t_0_xxxx_xy_xx[i] = t_0_xxxx_x_xxy[i] - rcd_y[i] * t_0_xxxx_x_xx[i];

        t_0_xxxx_xx_zz[i] = t_0_xxxx_x_xzz[i] - rcd_x[i] * t_0_xxxx_x_zz[i];

        t_0_xxxx_xx_yz[i] = t_0_xxxx_x_xyz[i] - rcd_x[i] * t_0_xxxx_x_yz[i];

        t_0_xxxx_xx_yy[i] = t_0_xxxx_x_xyy[i] - rcd_x[i] * t_0_xxxx_x_yy[i];

        t_0_xxxx_xx_xz[i] = t_0_xxxx_x_xxz[i] - rcd_x[i] * t_0_xxxx_x_xz[i];

        t_0_xxxx_xx_xy[i] = t_0_xxxx_x_xxy[i] - rcd_x[i] * t_0_xxxx_x_xy[i];

        t_0_xxxx_xx_xx[i] = t_0_xxxx_x_xxx[i] - rcd_x[i] * t_0_xxxx_x_xx[i];
    }
}


} // derirec namespace
