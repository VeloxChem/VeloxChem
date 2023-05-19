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
compHostHRRForDDSD_V0(      BufferHostXY<T>&      intsBufferDDSD,
                      const BufferHostX<int32_t>& intsIndexesDDSD,
                      const BufferHostXY<T>&      intsBufferPDSD,
                      const BufferHostX<int32_t>& intsIndexesPDSD,
                      const BufferHostXY<T>&      intsBufferPFSD,
                      const BufferHostX<int32_t>& intsIndexesPFSD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDSD) integral components

    t_zz_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(0));

    t_zz_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(1));

    t_zz_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(2));

    t_zz_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(3));

    t_zz_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(4));

    t_zz_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(5));

    t_zz_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(6));

    t_zz_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(7));

    t_zz_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(8));

    t_zz_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(9));

    t_zz_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(10));

    t_zz_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(11));

    t_zz_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(12));

    t_zz_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(13));

    t_zz_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(14));

    t_zz_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(15));

    t_zz_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(16));

    t_zz_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(17));

    t_zz_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(18));

    t_zz_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(19));

    t_zz_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(20));

    t_zz_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(21));

    t_zz_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(22));

    t_zz_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(23));

    t_zz_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(24));

    t_zz_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(25));

    t_zz_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(26));

    t_zz_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(27));

    t_zz_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(28));

    t_zz_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(29));

    t_zz_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(30));

    t_zz_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(31));

    t_zz_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(32));

    t_zz_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(33));

    t_zz_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(34));

    t_zz_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(35));

    t_yz_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(36));

    t_yz_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(37));

    t_yz_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(38));

    t_yz_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(39));

    t_yz_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(40));

    t_yz_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(41));

    t_yz_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(42));

    t_yz_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(43));

    t_yz_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(44));

    t_yz_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(45));

    t_yz_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(46));

    t_yz_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(47));

    t_yz_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(48));

    t_yz_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(49));

    t_yz_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(50));

    t_yz_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(51));

    t_yz_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(52));

    t_yz_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(53));

    t_yz_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(54));

    t_yz_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(55));

    t_yz_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(56));

    t_yz_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(57));

    t_yz_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(58));

    t_yz_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(59));

    t_yz_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(60));

    t_yz_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(61));

    t_yz_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(62));

    t_yz_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(63));

    t_yz_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(64));

    t_yz_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(65));

    t_yz_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(66));

    t_yz_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(67));

    t_yz_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(68));

    t_yz_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(69));

    t_yz_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(70));

    t_yz_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(71));

    t_yy_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(72));

    t_yy_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(73));

    t_yy_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(74));

    t_yy_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(75));

    t_yy_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(76));

    t_yy_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(77));

    t_yy_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(78));

    t_yy_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(79));

    t_yy_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(80));

    t_yy_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(81));

    t_yy_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(82));

    t_yy_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(83));

    t_yy_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(84));

    t_yy_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(85));

    t_yy_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(86));

    t_yy_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(87));

    t_yy_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(88));

    t_yy_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(89));

    t_yy_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(90));

    t_yy_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(91));

    t_yy_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(92));

    t_yy_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(93));

    t_yy_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(94));

    t_yy_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(95));

    t_yy_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(96));

    t_yy_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(97));

    t_yy_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(98));

    t_yy_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(99));

    t_yy_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(100));

    t_yy_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(101));

    t_yy_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(102));

    t_yy_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(103));

    t_yy_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(104));

    t_yy_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(105));

    t_yy_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(106));

    t_yy_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(107));

    t_xz_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(108));

    t_xz_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(109));

    t_xz_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(110));

    t_xz_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(111));

    t_xz_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(112));

    t_xz_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(113));

    t_xz_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(114));

    t_xz_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(115));

    t_xz_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(116));

    t_xz_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(117));

    t_xz_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(118));

    t_xz_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(119));

    t_xz_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(120));

    t_xz_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(121));

    t_xz_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(122));

    t_xz_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(123));

    t_xz_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(124));

    t_xz_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(125));

    t_xz_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(126));

    t_xz_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(127));

    t_xz_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(128));

    t_xz_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(129));

    t_xz_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(130));

    t_xz_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(131));

    t_xz_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(132));

    t_xz_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(133));

    t_xz_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(134));

    t_xz_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(135));

    t_xz_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(136));

    t_xz_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(137));

    t_xz_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(138));

    t_xz_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(139));

    t_xz_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(140));

    t_xz_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(141));

    t_xz_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(142));

    t_xz_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(143));

    t_xy_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(144));

    t_xy_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(145));

    t_xy_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(146));

    t_xy_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(147));

    t_xy_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(148));

    t_xy_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(149));

    t_xy_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(150));

    t_xy_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(151));

    t_xy_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(152));

    t_xy_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(153));

    t_xy_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(154));

    t_xy_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(155));

    t_xy_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(156));

    t_xy_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(157));

    t_xy_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(158));

    t_xy_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(159));

    t_xy_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(160));

    t_xy_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(161));

    t_xy_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(162));

    t_xy_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(163));

    t_xy_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(164));

    t_xy_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(165));

    t_xy_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(166));

    t_xy_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(167));

    t_xy_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(168));

    t_xy_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(169));

    t_xy_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(170));

    t_xy_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(171));

    t_xy_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(172));

    t_xy_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(173));

    t_xy_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(174));

    t_xy_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(175));

    t_xy_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(176));

    t_xy_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(177));

    t_xy_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(178));

    t_xy_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(179));

    t_xx_zz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(180));

    t_xx_zz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(181));

    t_xx_zz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(182));

    t_xx_zz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(183));

    t_xx_zz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(184));

    t_xx_zz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(185));

    t_xx_yz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(186));

    t_xx_yz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(187));

    t_xx_yz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(188));

    t_xx_yz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(189));

    t_xx_yz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(190));

    t_xx_yz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(191));

    t_xx_yy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(192));

    t_xx_yy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(193));

    t_xx_yy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(194));

    t_xx_yy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(195));

    t_xx_yy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(196));

    t_xx_yy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(197));

    t_xx_xz_0_zz = intsBufferDDSD.data(intsIndexesDDSD(198));

    t_xx_xz_0_yz = intsBufferDDSD.data(intsIndexesDDSD(199));

    t_xx_xz_0_yy = intsBufferDDSD.data(intsIndexesDDSD(200));

    t_xx_xz_0_xz = intsBufferDDSD.data(intsIndexesDDSD(201));

    t_xx_xz_0_xy = intsBufferDDSD.data(intsIndexesDDSD(202));

    t_xx_xz_0_xx = intsBufferDDSD.data(intsIndexesDDSD(203));

    t_xx_xy_0_zz = intsBufferDDSD.data(intsIndexesDDSD(204));

    t_xx_xy_0_yz = intsBufferDDSD.data(intsIndexesDDSD(205));

    t_xx_xy_0_yy = intsBufferDDSD.data(intsIndexesDDSD(206));

    t_xx_xy_0_xz = intsBufferDDSD.data(intsIndexesDDSD(207));

    t_xx_xy_0_xy = intsBufferDDSD.data(intsIndexesDDSD(208));

    t_xx_xy_0_xx = intsBufferDDSD.data(intsIndexesDDSD(209));

    t_xx_xx_0_zz = intsBufferDDSD.data(intsIndexesDDSD(210));

    t_xx_xx_0_yz = intsBufferDDSD.data(intsIndexesDDSD(211));

    t_xx_xx_0_yy = intsBufferDDSD.data(intsIndexesDDSD(212));

    t_xx_xx_0_xz = intsBufferDDSD.data(intsIndexesDDSD(213));

    t_xx_xx_0_xy = intsBufferDDSD.data(intsIndexesDDSD(214));

    t_xx_xx_0_xx = intsBufferDDSD.data(intsIndexesDDSD(215));

    // set up (PDSD) integral components

    t_z_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(0));

    t_z_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(1));

    t_z_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(2));

    t_z_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(3));

    t_z_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(4));

    t_z_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(5));

    t_z_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(6));

    t_z_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(7));

    t_z_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(8));

    t_z_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(9));

    t_z_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(10));

    t_z_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(11));

    t_z_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(12));

    t_z_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(13));

    t_z_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(14));

    t_z_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(15));

    t_z_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(16));

    t_z_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(17));

    t_z_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(18));

    t_z_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(19));

    t_z_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(20));

    t_z_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(21));

    t_z_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(22));

    t_z_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(23));

    t_z_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(24));

    t_z_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(25));

    t_z_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(26));

    t_z_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(27));

    t_z_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(28));

    t_z_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(29));

    t_z_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(30));

    t_z_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(31));

    t_z_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(32));

    t_z_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(33));

    t_z_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(34));

    t_z_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(35));

    t_y_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(36));

    t_y_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(37));

    t_y_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(38));

    t_y_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(39));

    t_y_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(40));

    t_y_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(41));

    t_y_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(42));

    t_y_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(43));

    t_y_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(44));

    t_y_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(45));

    t_y_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(46));

    t_y_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(47));

    t_y_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(48));

    t_y_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(49));

    t_y_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(50));

    t_y_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(51));

    t_y_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(52));

    t_y_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(53));

    t_y_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(54));

    t_y_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(55));

    t_y_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(56));

    t_y_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(57));

    t_y_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(58));

    t_y_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(59));

    t_y_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(60));

    t_y_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(61));

    t_y_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(62));

    t_y_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(63));

    t_y_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(64));

    t_y_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(65));

    t_y_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(66));

    t_y_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(67));

    t_y_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(68));

    t_y_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(69));

    t_y_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(70));

    t_y_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(71));

    t_x_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(72));

    t_x_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(73));

    t_x_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(74));

    t_x_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(75));

    t_x_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(76));

    t_x_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(77));

    t_x_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(78));

    t_x_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(79));

    t_x_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(80));

    t_x_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(81));

    t_x_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(82));

    t_x_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(83));

    t_x_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(84));

    t_x_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(85));

    t_x_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(86));

    t_x_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(87));

    t_x_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(88));

    t_x_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(89));

    t_x_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(90));

    t_x_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(91));

    t_x_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(92));

    t_x_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(93));

    t_x_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(94));

    t_x_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(95));

    t_x_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(96));

    t_x_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(97));

    t_x_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(98));

    t_x_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(99));

    t_x_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(100));

    t_x_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(101));

    t_x_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(102));

    t_x_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(103));

    t_x_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(104));

    t_x_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(105));

    t_x_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(106));

    t_x_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(107));

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

    #pragma omp simd align(rab_z, t_z_xx_0_xx, t_z_xx_0_xy, t_z_xx_0_xz, t_z_xx_0_yy,\
                           t_z_xx_0_yz, t_z_xx_0_zz, t_z_xxz_0_xx, t_z_xxz_0_xy, t_z_xxz_0_xz,\
                           t_z_xxz_0_yy, t_z_xxz_0_yz, t_z_xxz_0_zz, t_z_xy_0_xx, t_z_xy_0_xy,\
                           t_z_xy_0_xz, t_z_xy_0_yy, t_z_xy_0_yz, t_z_xy_0_zz, t_z_xyz_0_xx,\
                           t_z_xyz_0_xy, t_z_xyz_0_xz, t_z_xyz_0_yy, t_z_xyz_0_yz, t_z_xyz_0_zz,\
                           t_z_xz_0_xx, t_z_xz_0_xy, t_z_xz_0_xz, t_z_xz_0_yy, t_z_xz_0_yz,\
                           t_z_xz_0_zz, t_z_xzz_0_xx, t_z_xzz_0_xy, t_z_xzz_0_xz, t_z_xzz_0_yy,\
                           t_z_xzz_0_yz, t_z_xzz_0_zz, t_z_yy_0_xx, t_z_yy_0_xy, t_z_yy_0_xz,\
                           t_z_yy_0_yy, t_z_yy_0_yz, t_z_yy_0_zz, t_z_yyz_0_xx, t_z_yyz_0_xy,\
                           t_z_yyz_0_xz, t_z_yyz_0_yy, t_z_yyz_0_yz, t_z_yyz_0_zz, t_z_yz_0_xx,\
                           t_z_yz_0_xy, t_z_yz_0_xz, t_z_yz_0_yy, t_z_yz_0_yz, t_z_yz_0_zz,\
                           t_z_yzz_0_xx, t_z_yzz_0_xy, t_z_yzz_0_xz, t_z_yzz_0_yy, t_z_yzz_0_yz,\
                           t_z_yzz_0_zz, t_z_zz_0_xx, t_z_zz_0_xy, t_z_zz_0_xz, t_z_zz_0_yy,\
                           t_z_zz_0_yz, t_z_zz_0_zz, t_z_zzz_0_xx, t_z_zzz_0_xy, t_z_zzz_0_xz,\
                           t_z_zzz_0_yy, t_z_zzz_0_yz, t_z_zzz_0_zz, t_zz_xx_0_xx, t_zz_xx_0_xy,\
                           t_zz_xx_0_xz, t_zz_xx_0_yy, t_zz_xx_0_yz, t_zz_xx_0_zz, t_zz_xy_0_xx,\
                           t_zz_xy_0_xy, t_zz_xy_0_xz, t_zz_xy_0_yy, t_zz_xy_0_yz, t_zz_xy_0_zz,\
                           t_zz_xz_0_xx, t_zz_xz_0_xy, t_zz_xz_0_xz, t_zz_xz_0_yy, t_zz_xz_0_yz,\
                           t_zz_xz_0_zz, t_zz_yy_0_xx, t_zz_yy_0_xy, t_zz_yy_0_xz, t_zz_yy_0_yy,\
                           t_zz_yy_0_yz, t_zz_yy_0_zz, t_zz_yz_0_xx, t_zz_yz_0_xy, t_zz_yz_0_xz,\
                           t_zz_yz_0_yy, t_zz_yz_0_yz, t_zz_yz_0_zz, t_zz_zz_0_xx, t_zz_zz_0_xy,\
                           t_zz_zz_0_xz, t_zz_zz_0_yy, t_zz_zz_0_yz, t_zz_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_0_zz[i] = t_z_zzz_0_zz[i] - rab_z[i] * t_z_zz_0_zz[i];

        t_zz_zz_0_yz[i] = t_z_zzz_0_yz[i] - rab_z[i] * t_z_zz_0_yz[i];

        t_zz_zz_0_yy[i] = t_z_zzz_0_yy[i] - rab_z[i] * t_z_zz_0_yy[i];

        t_zz_zz_0_xz[i] = t_z_zzz_0_xz[i] - rab_z[i] * t_z_zz_0_xz[i];

        t_zz_zz_0_xy[i] = t_z_zzz_0_xy[i] - rab_z[i] * t_z_zz_0_xy[i];

        t_zz_zz_0_xx[i] = t_z_zzz_0_xx[i] - rab_z[i] * t_z_zz_0_xx[i];

        t_zz_yz_0_zz[i] = t_z_yzz_0_zz[i] - rab_z[i] * t_z_yz_0_zz[i];

        t_zz_yz_0_yz[i] = t_z_yzz_0_yz[i] - rab_z[i] * t_z_yz_0_yz[i];

        t_zz_yz_0_yy[i] = t_z_yzz_0_yy[i] - rab_z[i] * t_z_yz_0_yy[i];

        t_zz_yz_0_xz[i] = t_z_yzz_0_xz[i] - rab_z[i] * t_z_yz_0_xz[i];

        t_zz_yz_0_xy[i] = t_z_yzz_0_xy[i] - rab_z[i] * t_z_yz_0_xy[i];

        t_zz_yz_0_xx[i] = t_z_yzz_0_xx[i] - rab_z[i] * t_z_yz_0_xx[i];

        t_zz_yy_0_zz[i] = t_z_yyz_0_zz[i] - rab_z[i] * t_z_yy_0_zz[i];

        t_zz_yy_0_yz[i] = t_z_yyz_0_yz[i] - rab_z[i] * t_z_yy_0_yz[i];

        t_zz_yy_0_yy[i] = t_z_yyz_0_yy[i] - rab_z[i] * t_z_yy_0_yy[i];

        t_zz_yy_0_xz[i] = t_z_yyz_0_xz[i] - rab_z[i] * t_z_yy_0_xz[i];

        t_zz_yy_0_xy[i] = t_z_yyz_0_xy[i] - rab_z[i] * t_z_yy_0_xy[i];

        t_zz_yy_0_xx[i] = t_z_yyz_0_xx[i] - rab_z[i] * t_z_yy_0_xx[i];

        t_zz_xz_0_zz[i] = t_z_xzz_0_zz[i] - rab_z[i] * t_z_xz_0_zz[i];

        t_zz_xz_0_yz[i] = t_z_xzz_0_yz[i] - rab_z[i] * t_z_xz_0_yz[i];

        t_zz_xz_0_yy[i] = t_z_xzz_0_yy[i] - rab_z[i] * t_z_xz_0_yy[i];

        t_zz_xz_0_xz[i] = t_z_xzz_0_xz[i] - rab_z[i] * t_z_xz_0_xz[i];

        t_zz_xz_0_xy[i] = t_z_xzz_0_xy[i] - rab_z[i] * t_z_xz_0_xy[i];

        t_zz_xz_0_xx[i] = t_z_xzz_0_xx[i] - rab_z[i] * t_z_xz_0_xx[i];

        t_zz_xy_0_zz[i] = t_z_xyz_0_zz[i] - rab_z[i] * t_z_xy_0_zz[i];

        t_zz_xy_0_yz[i] = t_z_xyz_0_yz[i] - rab_z[i] * t_z_xy_0_yz[i];

        t_zz_xy_0_yy[i] = t_z_xyz_0_yy[i] - rab_z[i] * t_z_xy_0_yy[i];

        t_zz_xy_0_xz[i] = t_z_xyz_0_xz[i] - rab_z[i] * t_z_xy_0_xz[i];

        t_zz_xy_0_xy[i] = t_z_xyz_0_xy[i] - rab_z[i] * t_z_xy_0_xy[i];

        t_zz_xy_0_xx[i] = t_z_xyz_0_xx[i] - rab_z[i] * t_z_xy_0_xx[i];

        t_zz_xx_0_zz[i] = t_z_xxz_0_zz[i] - rab_z[i] * t_z_xx_0_zz[i];

        t_zz_xx_0_yz[i] = t_z_xxz_0_yz[i] - rab_z[i] * t_z_xx_0_yz[i];

        t_zz_xx_0_yy[i] = t_z_xxz_0_yy[i] - rab_z[i] * t_z_xx_0_yy[i];

        t_zz_xx_0_xz[i] = t_z_xxz_0_xz[i] - rab_z[i] * t_z_xx_0_xz[i];

        t_zz_xx_0_xy[i] = t_z_xxz_0_xy[i] - rab_z[i] * t_z_xx_0_xy[i];

        t_zz_xx_0_xx[i] = t_z_xxz_0_xx[i] - rab_z[i] * t_z_xx_0_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xx_0_xx, t_y_xx_0_xy, t_y_xx_0_xz, t_y_xx_0_yy,\
                           t_y_xx_0_yz, t_y_xx_0_zz, t_y_xxz_0_xx, t_y_xxz_0_xy, t_y_xxz_0_xz,\
                           t_y_xxz_0_yy, t_y_xxz_0_yz, t_y_xxz_0_zz, t_y_xy_0_xx, t_y_xy_0_xy,\
                           t_y_xy_0_xz, t_y_xy_0_yy, t_y_xy_0_yz, t_y_xy_0_zz, t_y_xyz_0_xx,\
                           t_y_xyz_0_xy, t_y_xyz_0_xz, t_y_xyz_0_yy, t_y_xyz_0_yz, t_y_xyz_0_zz,\
                           t_y_xz_0_xx, t_y_xz_0_xy, t_y_xz_0_xz, t_y_xz_0_yy, t_y_xz_0_yz,\
                           t_y_xz_0_zz, t_y_xzz_0_xx, t_y_xzz_0_xy, t_y_xzz_0_xz, t_y_xzz_0_yy,\
                           t_y_xzz_0_yz, t_y_xzz_0_zz, t_y_yy_0_xx, t_y_yy_0_xy, t_y_yy_0_xz,\
                           t_y_yy_0_yy, t_y_yy_0_yz, t_y_yy_0_zz, t_y_yyz_0_xx, t_y_yyz_0_xy,\
                           t_y_yyz_0_xz, t_y_yyz_0_yy, t_y_yyz_0_yz, t_y_yyz_0_zz, t_y_yz_0_xx,\
                           t_y_yz_0_xy, t_y_yz_0_xz, t_y_yz_0_yy, t_y_yz_0_yz, t_y_yz_0_zz,\
                           t_y_yzz_0_xx, t_y_yzz_0_xy, t_y_yzz_0_xz, t_y_yzz_0_yy, t_y_yzz_0_yz,\
                           t_y_yzz_0_zz, t_y_zz_0_xx, t_y_zz_0_xy, t_y_zz_0_xz, t_y_zz_0_yy,\
                           t_y_zz_0_yz, t_y_zz_0_zz, t_y_zzz_0_xx, t_y_zzz_0_xy, t_y_zzz_0_xz,\
                           t_y_zzz_0_yy, t_y_zzz_0_yz, t_y_zzz_0_zz, t_yz_xx_0_xx, t_yz_xx_0_xy,\
                           t_yz_xx_0_xz, t_yz_xx_0_yy, t_yz_xx_0_yz, t_yz_xx_0_zz, t_yz_xy_0_xx,\
                           t_yz_xy_0_xy, t_yz_xy_0_xz, t_yz_xy_0_yy, t_yz_xy_0_yz, t_yz_xy_0_zz,\
                           t_yz_xz_0_xx, t_yz_xz_0_xy, t_yz_xz_0_xz, t_yz_xz_0_yy, t_yz_xz_0_yz,\
                           t_yz_xz_0_zz, t_yz_yy_0_xx, t_yz_yy_0_xy, t_yz_yy_0_xz, t_yz_yy_0_yy,\
                           t_yz_yy_0_yz, t_yz_yy_0_zz, t_yz_yz_0_xx, t_yz_yz_0_xy, t_yz_yz_0_xz,\
                           t_yz_yz_0_yy, t_yz_yz_0_yz, t_yz_yz_0_zz, t_yz_zz_0_xx, t_yz_zz_0_xy,\
                           t_yz_zz_0_xz, t_yz_zz_0_yy, t_yz_zz_0_yz, t_yz_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_zz_0_zz[i] = t_y_zzz_0_zz[i] - rab_z[i] * t_y_zz_0_zz[i];

        t_yz_zz_0_yz[i] = t_y_zzz_0_yz[i] - rab_z[i] * t_y_zz_0_yz[i];

        t_yz_zz_0_yy[i] = t_y_zzz_0_yy[i] - rab_z[i] * t_y_zz_0_yy[i];

        t_yz_zz_0_xz[i] = t_y_zzz_0_xz[i] - rab_z[i] * t_y_zz_0_xz[i];

        t_yz_zz_0_xy[i] = t_y_zzz_0_xy[i] - rab_z[i] * t_y_zz_0_xy[i];

        t_yz_zz_0_xx[i] = t_y_zzz_0_xx[i] - rab_z[i] * t_y_zz_0_xx[i];

        t_yz_yz_0_zz[i] = t_y_yzz_0_zz[i] - rab_z[i] * t_y_yz_0_zz[i];

        t_yz_yz_0_yz[i] = t_y_yzz_0_yz[i] - rab_z[i] * t_y_yz_0_yz[i];

        t_yz_yz_0_yy[i] = t_y_yzz_0_yy[i] - rab_z[i] * t_y_yz_0_yy[i];

        t_yz_yz_0_xz[i] = t_y_yzz_0_xz[i] - rab_z[i] * t_y_yz_0_xz[i];

        t_yz_yz_0_xy[i] = t_y_yzz_0_xy[i] - rab_z[i] * t_y_yz_0_xy[i];

        t_yz_yz_0_xx[i] = t_y_yzz_0_xx[i] - rab_z[i] * t_y_yz_0_xx[i];

        t_yz_yy_0_zz[i] = t_y_yyz_0_zz[i] - rab_z[i] * t_y_yy_0_zz[i];

        t_yz_yy_0_yz[i] = t_y_yyz_0_yz[i] - rab_z[i] * t_y_yy_0_yz[i];

        t_yz_yy_0_yy[i] = t_y_yyz_0_yy[i] - rab_z[i] * t_y_yy_0_yy[i];

        t_yz_yy_0_xz[i] = t_y_yyz_0_xz[i] - rab_z[i] * t_y_yy_0_xz[i];

        t_yz_yy_0_xy[i] = t_y_yyz_0_xy[i] - rab_z[i] * t_y_yy_0_xy[i];

        t_yz_yy_0_xx[i] = t_y_yyz_0_xx[i] - rab_z[i] * t_y_yy_0_xx[i];

        t_yz_xz_0_zz[i] = t_y_xzz_0_zz[i] - rab_z[i] * t_y_xz_0_zz[i];

        t_yz_xz_0_yz[i] = t_y_xzz_0_yz[i] - rab_z[i] * t_y_xz_0_yz[i];

        t_yz_xz_0_yy[i] = t_y_xzz_0_yy[i] - rab_z[i] * t_y_xz_0_yy[i];

        t_yz_xz_0_xz[i] = t_y_xzz_0_xz[i] - rab_z[i] * t_y_xz_0_xz[i];

        t_yz_xz_0_xy[i] = t_y_xzz_0_xy[i] - rab_z[i] * t_y_xz_0_xy[i];

        t_yz_xz_0_xx[i] = t_y_xzz_0_xx[i] - rab_z[i] * t_y_xz_0_xx[i];

        t_yz_xy_0_zz[i] = t_y_xyz_0_zz[i] - rab_z[i] * t_y_xy_0_zz[i];

        t_yz_xy_0_yz[i] = t_y_xyz_0_yz[i] - rab_z[i] * t_y_xy_0_yz[i];

        t_yz_xy_0_yy[i] = t_y_xyz_0_yy[i] - rab_z[i] * t_y_xy_0_yy[i];

        t_yz_xy_0_xz[i] = t_y_xyz_0_xz[i] - rab_z[i] * t_y_xy_0_xz[i];

        t_yz_xy_0_xy[i] = t_y_xyz_0_xy[i] - rab_z[i] * t_y_xy_0_xy[i];

        t_yz_xy_0_xx[i] = t_y_xyz_0_xx[i] - rab_z[i] * t_y_xy_0_xx[i];

        t_yz_xx_0_zz[i] = t_y_xxz_0_zz[i] - rab_z[i] * t_y_xx_0_zz[i];

        t_yz_xx_0_yz[i] = t_y_xxz_0_yz[i] - rab_z[i] * t_y_xx_0_yz[i];

        t_yz_xx_0_yy[i] = t_y_xxz_0_yy[i] - rab_z[i] * t_y_xx_0_yy[i];

        t_yz_xx_0_xz[i] = t_y_xxz_0_xz[i] - rab_z[i] * t_y_xx_0_xz[i];

        t_yz_xx_0_xy[i] = t_y_xxz_0_xy[i] - rab_z[i] * t_y_xx_0_xy[i];

        t_yz_xx_0_xx[i] = t_y_xxz_0_xx[i] - rab_z[i] * t_y_xx_0_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xx_0_xx, t_y_xx_0_xy, t_y_xx_0_xz, t_y_xx_0_yy,\
                           t_y_xx_0_yz, t_y_xx_0_zz, t_y_xxy_0_xx, t_y_xxy_0_xy, t_y_xxy_0_xz,\
                           t_y_xxy_0_yy, t_y_xxy_0_yz, t_y_xxy_0_zz, t_y_xy_0_xx, t_y_xy_0_xy,\
                           t_y_xy_0_xz, t_y_xy_0_yy, t_y_xy_0_yz, t_y_xy_0_zz, t_y_xyy_0_xx,\
                           t_y_xyy_0_xy, t_y_xyy_0_xz, t_y_xyy_0_yy, t_y_xyy_0_yz, t_y_xyy_0_zz,\
                           t_y_xyz_0_xx, t_y_xyz_0_xy, t_y_xyz_0_xz, t_y_xyz_0_yy, t_y_xyz_0_yz,\
                           t_y_xyz_0_zz, t_y_xz_0_xx, t_y_xz_0_xy, t_y_xz_0_xz, t_y_xz_0_yy,\
                           t_y_xz_0_yz, t_y_xz_0_zz, t_y_yy_0_xx, t_y_yy_0_xy, t_y_yy_0_xz,\
                           t_y_yy_0_yy, t_y_yy_0_yz, t_y_yy_0_zz, t_y_yyy_0_xx, t_y_yyy_0_xy,\
                           t_y_yyy_0_xz, t_y_yyy_0_yy, t_y_yyy_0_yz, t_y_yyy_0_zz, t_y_yyz_0_xx,\
                           t_y_yyz_0_xy, t_y_yyz_0_xz, t_y_yyz_0_yy, t_y_yyz_0_yz, t_y_yyz_0_zz,\
                           t_y_yz_0_xx, t_y_yz_0_xy, t_y_yz_0_xz, t_y_yz_0_yy, t_y_yz_0_yz,\
                           t_y_yz_0_zz, t_y_yzz_0_xx, t_y_yzz_0_xy, t_y_yzz_0_xz, t_y_yzz_0_yy,\
                           t_y_yzz_0_yz, t_y_yzz_0_zz, t_y_zz_0_xx, t_y_zz_0_xy, t_y_zz_0_xz,\
                           t_y_zz_0_yy, t_y_zz_0_yz, t_y_zz_0_zz, t_yy_xx_0_xx, t_yy_xx_0_xy,\
                           t_yy_xx_0_xz, t_yy_xx_0_yy, t_yy_xx_0_yz, t_yy_xx_0_zz, t_yy_xy_0_xx,\
                           t_yy_xy_0_xy, t_yy_xy_0_xz, t_yy_xy_0_yy, t_yy_xy_0_yz, t_yy_xy_0_zz,\
                           t_yy_xz_0_xx, t_yy_xz_0_xy, t_yy_xz_0_xz, t_yy_xz_0_yy, t_yy_xz_0_yz,\
                           t_yy_xz_0_zz, t_yy_yy_0_xx, t_yy_yy_0_xy, t_yy_yy_0_xz, t_yy_yy_0_yy,\
                           t_yy_yy_0_yz, t_yy_yy_0_zz, t_yy_yz_0_xx, t_yy_yz_0_xy, t_yy_yz_0_xz,\
                           t_yy_yz_0_yy, t_yy_yz_0_yz, t_yy_yz_0_zz, t_yy_zz_0_xx, t_yy_zz_0_xy,\
                           t_yy_zz_0_xz, t_yy_zz_0_yy, t_yy_zz_0_yz, t_yy_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_zz_0_zz[i] = t_y_yzz_0_zz[i] - rab_y[i] * t_y_zz_0_zz[i];

        t_yy_zz_0_yz[i] = t_y_yzz_0_yz[i] - rab_y[i] * t_y_zz_0_yz[i];

        t_yy_zz_0_yy[i] = t_y_yzz_0_yy[i] - rab_y[i] * t_y_zz_0_yy[i];

        t_yy_zz_0_xz[i] = t_y_yzz_0_xz[i] - rab_y[i] * t_y_zz_0_xz[i];

        t_yy_zz_0_xy[i] = t_y_yzz_0_xy[i] - rab_y[i] * t_y_zz_0_xy[i];

        t_yy_zz_0_xx[i] = t_y_yzz_0_xx[i] - rab_y[i] * t_y_zz_0_xx[i];

        t_yy_yz_0_zz[i] = t_y_yyz_0_zz[i] - rab_y[i] * t_y_yz_0_zz[i];

        t_yy_yz_0_yz[i] = t_y_yyz_0_yz[i] - rab_y[i] * t_y_yz_0_yz[i];

        t_yy_yz_0_yy[i] = t_y_yyz_0_yy[i] - rab_y[i] * t_y_yz_0_yy[i];

        t_yy_yz_0_xz[i] = t_y_yyz_0_xz[i] - rab_y[i] * t_y_yz_0_xz[i];

        t_yy_yz_0_xy[i] = t_y_yyz_0_xy[i] - rab_y[i] * t_y_yz_0_xy[i];

        t_yy_yz_0_xx[i] = t_y_yyz_0_xx[i] - rab_y[i] * t_y_yz_0_xx[i];

        t_yy_yy_0_zz[i] = t_y_yyy_0_zz[i] - rab_y[i] * t_y_yy_0_zz[i];

        t_yy_yy_0_yz[i] = t_y_yyy_0_yz[i] - rab_y[i] * t_y_yy_0_yz[i];

        t_yy_yy_0_yy[i] = t_y_yyy_0_yy[i] - rab_y[i] * t_y_yy_0_yy[i];

        t_yy_yy_0_xz[i] = t_y_yyy_0_xz[i] - rab_y[i] * t_y_yy_0_xz[i];

        t_yy_yy_0_xy[i] = t_y_yyy_0_xy[i] - rab_y[i] * t_y_yy_0_xy[i];

        t_yy_yy_0_xx[i] = t_y_yyy_0_xx[i] - rab_y[i] * t_y_yy_0_xx[i];

        t_yy_xz_0_zz[i] = t_y_xyz_0_zz[i] - rab_y[i] * t_y_xz_0_zz[i];

        t_yy_xz_0_yz[i] = t_y_xyz_0_yz[i] - rab_y[i] * t_y_xz_0_yz[i];

        t_yy_xz_0_yy[i] = t_y_xyz_0_yy[i] - rab_y[i] * t_y_xz_0_yy[i];

        t_yy_xz_0_xz[i] = t_y_xyz_0_xz[i] - rab_y[i] * t_y_xz_0_xz[i];

        t_yy_xz_0_xy[i] = t_y_xyz_0_xy[i] - rab_y[i] * t_y_xz_0_xy[i];

        t_yy_xz_0_xx[i] = t_y_xyz_0_xx[i] - rab_y[i] * t_y_xz_0_xx[i];

        t_yy_xy_0_zz[i] = t_y_xyy_0_zz[i] - rab_y[i] * t_y_xy_0_zz[i];

        t_yy_xy_0_yz[i] = t_y_xyy_0_yz[i] - rab_y[i] * t_y_xy_0_yz[i];

        t_yy_xy_0_yy[i] = t_y_xyy_0_yy[i] - rab_y[i] * t_y_xy_0_yy[i];

        t_yy_xy_0_xz[i] = t_y_xyy_0_xz[i] - rab_y[i] * t_y_xy_0_xz[i];

        t_yy_xy_0_xy[i] = t_y_xyy_0_xy[i] - rab_y[i] * t_y_xy_0_xy[i];

        t_yy_xy_0_xx[i] = t_y_xyy_0_xx[i] - rab_y[i] * t_y_xy_0_xx[i];

        t_yy_xx_0_zz[i] = t_y_xxy_0_zz[i] - rab_y[i] * t_y_xx_0_zz[i];

        t_yy_xx_0_yz[i] = t_y_xxy_0_yz[i] - rab_y[i] * t_y_xx_0_yz[i];

        t_yy_xx_0_yy[i] = t_y_xxy_0_yy[i] - rab_y[i] * t_y_xx_0_yy[i];

        t_yy_xx_0_xz[i] = t_y_xxy_0_xz[i] - rab_y[i] * t_y_xx_0_xz[i];

        t_yy_xx_0_xy[i] = t_y_xxy_0_xy[i] - rab_y[i] * t_y_xx_0_xy[i];

        t_yy_xx_0_xx[i] = t_y_xxy_0_xx[i] - rab_y[i] * t_y_xx_0_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xx_0_xx, t_x_xx_0_xy, t_x_xx_0_xz, t_x_xx_0_yy,\
                           t_x_xx_0_yz, t_x_xx_0_zz, t_x_xxz_0_xx, t_x_xxz_0_xy, t_x_xxz_0_xz,\
                           t_x_xxz_0_yy, t_x_xxz_0_yz, t_x_xxz_0_zz, t_x_xy_0_xx, t_x_xy_0_xy,\
                           t_x_xy_0_xz, t_x_xy_0_yy, t_x_xy_0_yz, t_x_xy_0_zz, t_x_xyz_0_xx,\
                           t_x_xyz_0_xy, t_x_xyz_0_xz, t_x_xyz_0_yy, t_x_xyz_0_yz, t_x_xyz_0_zz,\
                           t_x_xz_0_xx, t_x_xz_0_xy, t_x_xz_0_xz, t_x_xz_0_yy, t_x_xz_0_yz,\
                           t_x_xz_0_zz, t_x_xzz_0_xx, t_x_xzz_0_xy, t_x_xzz_0_xz, t_x_xzz_0_yy,\
                           t_x_xzz_0_yz, t_x_xzz_0_zz, t_x_yy_0_xx, t_x_yy_0_xy, t_x_yy_0_xz,\
                           t_x_yy_0_yy, t_x_yy_0_yz, t_x_yy_0_zz, t_x_yyz_0_xx, t_x_yyz_0_xy,\
                           t_x_yyz_0_xz, t_x_yyz_0_yy, t_x_yyz_0_yz, t_x_yyz_0_zz, t_x_yz_0_xx,\
                           t_x_yz_0_xy, t_x_yz_0_xz, t_x_yz_0_yy, t_x_yz_0_yz, t_x_yz_0_zz,\
                           t_x_yzz_0_xx, t_x_yzz_0_xy, t_x_yzz_0_xz, t_x_yzz_0_yy, t_x_yzz_0_yz,\
                           t_x_yzz_0_zz, t_x_zz_0_xx, t_x_zz_0_xy, t_x_zz_0_xz, t_x_zz_0_yy,\
                           t_x_zz_0_yz, t_x_zz_0_zz, t_x_zzz_0_xx, t_x_zzz_0_xy, t_x_zzz_0_xz,\
                           t_x_zzz_0_yy, t_x_zzz_0_yz, t_x_zzz_0_zz, t_xz_xx_0_xx, t_xz_xx_0_xy,\
                           t_xz_xx_0_xz, t_xz_xx_0_yy, t_xz_xx_0_yz, t_xz_xx_0_zz, t_xz_xy_0_xx,\
                           t_xz_xy_0_xy, t_xz_xy_0_xz, t_xz_xy_0_yy, t_xz_xy_0_yz, t_xz_xy_0_zz,\
                           t_xz_xz_0_xx, t_xz_xz_0_xy, t_xz_xz_0_xz, t_xz_xz_0_yy, t_xz_xz_0_yz,\
                           t_xz_xz_0_zz, t_xz_yy_0_xx, t_xz_yy_0_xy, t_xz_yy_0_xz, t_xz_yy_0_yy,\
                           t_xz_yy_0_yz, t_xz_yy_0_zz, t_xz_yz_0_xx, t_xz_yz_0_xy, t_xz_yz_0_xz,\
                           t_xz_yz_0_yy, t_xz_yz_0_yz, t_xz_yz_0_zz, t_xz_zz_0_xx, t_xz_zz_0_xy,\
                           t_xz_zz_0_xz, t_xz_zz_0_yy, t_xz_zz_0_yz, t_xz_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_zz_0_zz[i] = t_x_zzz_0_zz[i] - rab_z[i] * t_x_zz_0_zz[i];

        t_xz_zz_0_yz[i] = t_x_zzz_0_yz[i] - rab_z[i] * t_x_zz_0_yz[i];

        t_xz_zz_0_yy[i] = t_x_zzz_0_yy[i] - rab_z[i] * t_x_zz_0_yy[i];

        t_xz_zz_0_xz[i] = t_x_zzz_0_xz[i] - rab_z[i] * t_x_zz_0_xz[i];

        t_xz_zz_0_xy[i] = t_x_zzz_0_xy[i] - rab_z[i] * t_x_zz_0_xy[i];

        t_xz_zz_0_xx[i] = t_x_zzz_0_xx[i] - rab_z[i] * t_x_zz_0_xx[i];

        t_xz_yz_0_zz[i] = t_x_yzz_0_zz[i] - rab_z[i] * t_x_yz_0_zz[i];

        t_xz_yz_0_yz[i] = t_x_yzz_0_yz[i] - rab_z[i] * t_x_yz_0_yz[i];

        t_xz_yz_0_yy[i] = t_x_yzz_0_yy[i] - rab_z[i] * t_x_yz_0_yy[i];

        t_xz_yz_0_xz[i] = t_x_yzz_0_xz[i] - rab_z[i] * t_x_yz_0_xz[i];

        t_xz_yz_0_xy[i] = t_x_yzz_0_xy[i] - rab_z[i] * t_x_yz_0_xy[i];

        t_xz_yz_0_xx[i] = t_x_yzz_0_xx[i] - rab_z[i] * t_x_yz_0_xx[i];

        t_xz_yy_0_zz[i] = t_x_yyz_0_zz[i] - rab_z[i] * t_x_yy_0_zz[i];

        t_xz_yy_0_yz[i] = t_x_yyz_0_yz[i] - rab_z[i] * t_x_yy_0_yz[i];

        t_xz_yy_0_yy[i] = t_x_yyz_0_yy[i] - rab_z[i] * t_x_yy_0_yy[i];

        t_xz_yy_0_xz[i] = t_x_yyz_0_xz[i] - rab_z[i] * t_x_yy_0_xz[i];

        t_xz_yy_0_xy[i] = t_x_yyz_0_xy[i] - rab_z[i] * t_x_yy_0_xy[i];

        t_xz_yy_0_xx[i] = t_x_yyz_0_xx[i] - rab_z[i] * t_x_yy_0_xx[i];

        t_xz_xz_0_zz[i] = t_x_xzz_0_zz[i] - rab_z[i] * t_x_xz_0_zz[i];

        t_xz_xz_0_yz[i] = t_x_xzz_0_yz[i] - rab_z[i] * t_x_xz_0_yz[i];

        t_xz_xz_0_yy[i] = t_x_xzz_0_yy[i] - rab_z[i] * t_x_xz_0_yy[i];

        t_xz_xz_0_xz[i] = t_x_xzz_0_xz[i] - rab_z[i] * t_x_xz_0_xz[i];

        t_xz_xz_0_xy[i] = t_x_xzz_0_xy[i] - rab_z[i] * t_x_xz_0_xy[i];

        t_xz_xz_0_xx[i] = t_x_xzz_0_xx[i] - rab_z[i] * t_x_xz_0_xx[i];

        t_xz_xy_0_zz[i] = t_x_xyz_0_zz[i] - rab_z[i] * t_x_xy_0_zz[i];

        t_xz_xy_0_yz[i] = t_x_xyz_0_yz[i] - rab_z[i] * t_x_xy_0_yz[i];

        t_xz_xy_0_yy[i] = t_x_xyz_0_yy[i] - rab_z[i] * t_x_xy_0_yy[i];

        t_xz_xy_0_xz[i] = t_x_xyz_0_xz[i] - rab_z[i] * t_x_xy_0_xz[i];

        t_xz_xy_0_xy[i] = t_x_xyz_0_xy[i] - rab_z[i] * t_x_xy_0_xy[i];

        t_xz_xy_0_xx[i] = t_x_xyz_0_xx[i] - rab_z[i] * t_x_xy_0_xx[i];

        t_xz_xx_0_zz[i] = t_x_xxz_0_zz[i] - rab_z[i] * t_x_xx_0_zz[i];

        t_xz_xx_0_yz[i] = t_x_xxz_0_yz[i] - rab_z[i] * t_x_xx_0_yz[i];

        t_xz_xx_0_yy[i] = t_x_xxz_0_yy[i] - rab_z[i] * t_x_xx_0_yy[i];

        t_xz_xx_0_xz[i] = t_x_xxz_0_xz[i] - rab_z[i] * t_x_xx_0_xz[i];

        t_xz_xx_0_xy[i] = t_x_xxz_0_xy[i] - rab_z[i] * t_x_xx_0_xy[i];

        t_xz_xx_0_xx[i] = t_x_xxz_0_xx[i] - rab_z[i] * t_x_xx_0_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xx_0_xx, t_x_xx_0_xy, t_x_xx_0_xz, t_x_xx_0_yy,\
                           t_x_xx_0_yz, t_x_xx_0_zz, t_x_xxy_0_xx, t_x_xxy_0_xy, t_x_xxy_0_xz,\
                           t_x_xxy_0_yy, t_x_xxy_0_yz, t_x_xxy_0_zz, t_x_xy_0_xx, t_x_xy_0_xy,\
                           t_x_xy_0_xz, t_x_xy_0_yy, t_x_xy_0_yz, t_x_xy_0_zz, t_x_xyy_0_xx,\
                           t_x_xyy_0_xy, t_x_xyy_0_xz, t_x_xyy_0_yy, t_x_xyy_0_yz, t_x_xyy_0_zz,\
                           t_x_xyz_0_xx, t_x_xyz_0_xy, t_x_xyz_0_xz, t_x_xyz_0_yy, t_x_xyz_0_yz,\
                           t_x_xyz_0_zz, t_x_xz_0_xx, t_x_xz_0_xy, t_x_xz_0_xz, t_x_xz_0_yy,\
                           t_x_xz_0_yz, t_x_xz_0_zz, t_x_yy_0_xx, t_x_yy_0_xy, t_x_yy_0_xz,\
                           t_x_yy_0_yy, t_x_yy_0_yz, t_x_yy_0_zz, t_x_yyy_0_xx, t_x_yyy_0_xy,\
                           t_x_yyy_0_xz, t_x_yyy_0_yy, t_x_yyy_0_yz, t_x_yyy_0_zz, t_x_yyz_0_xx,\
                           t_x_yyz_0_xy, t_x_yyz_0_xz, t_x_yyz_0_yy, t_x_yyz_0_yz, t_x_yyz_0_zz,\
                           t_x_yz_0_xx, t_x_yz_0_xy, t_x_yz_0_xz, t_x_yz_0_yy, t_x_yz_0_yz,\
                           t_x_yz_0_zz, t_x_yzz_0_xx, t_x_yzz_0_xy, t_x_yzz_0_xz, t_x_yzz_0_yy,\
                           t_x_yzz_0_yz, t_x_yzz_0_zz, t_x_zz_0_xx, t_x_zz_0_xy, t_x_zz_0_xz,\
                           t_x_zz_0_yy, t_x_zz_0_yz, t_x_zz_0_zz, t_xy_xx_0_xx, t_xy_xx_0_xy,\
                           t_xy_xx_0_xz, t_xy_xx_0_yy, t_xy_xx_0_yz, t_xy_xx_0_zz, t_xy_xy_0_xx,\
                           t_xy_xy_0_xy, t_xy_xy_0_xz, t_xy_xy_0_yy, t_xy_xy_0_yz, t_xy_xy_0_zz,\
                           t_xy_xz_0_xx, t_xy_xz_0_xy, t_xy_xz_0_xz, t_xy_xz_0_yy, t_xy_xz_0_yz,\
                           t_xy_xz_0_zz, t_xy_yy_0_xx, t_xy_yy_0_xy, t_xy_yy_0_xz, t_xy_yy_0_yy,\
                           t_xy_yy_0_yz, t_xy_yy_0_zz, t_xy_yz_0_xx, t_xy_yz_0_xy, t_xy_yz_0_xz,\
                           t_xy_yz_0_yy, t_xy_yz_0_yz, t_xy_yz_0_zz, t_xy_zz_0_xx, t_xy_zz_0_xy,\
                           t_xy_zz_0_xz, t_xy_zz_0_yy, t_xy_zz_0_yz, t_xy_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_zz_0_zz[i] = t_x_yzz_0_zz[i] - rab_y[i] * t_x_zz_0_zz[i];

        t_xy_zz_0_yz[i] = t_x_yzz_0_yz[i] - rab_y[i] * t_x_zz_0_yz[i];

        t_xy_zz_0_yy[i] = t_x_yzz_0_yy[i] - rab_y[i] * t_x_zz_0_yy[i];

        t_xy_zz_0_xz[i] = t_x_yzz_0_xz[i] - rab_y[i] * t_x_zz_0_xz[i];

        t_xy_zz_0_xy[i] = t_x_yzz_0_xy[i] - rab_y[i] * t_x_zz_0_xy[i];

        t_xy_zz_0_xx[i] = t_x_yzz_0_xx[i] - rab_y[i] * t_x_zz_0_xx[i];

        t_xy_yz_0_zz[i] = t_x_yyz_0_zz[i] - rab_y[i] * t_x_yz_0_zz[i];

        t_xy_yz_0_yz[i] = t_x_yyz_0_yz[i] - rab_y[i] * t_x_yz_0_yz[i];

        t_xy_yz_0_yy[i] = t_x_yyz_0_yy[i] - rab_y[i] * t_x_yz_0_yy[i];

        t_xy_yz_0_xz[i] = t_x_yyz_0_xz[i] - rab_y[i] * t_x_yz_0_xz[i];

        t_xy_yz_0_xy[i] = t_x_yyz_0_xy[i] - rab_y[i] * t_x_yz_0_xy[i];

        t_xy_yz_0_xx[i] = t_x_yyz_0_xx[i] - rab_y[i] * t_x_yz_0_xx[i];

        t_xy_yy_0_zz[i] = t_x_yyy_0_zz[i] - rab_y[i] * t_x_yy_0_zz[i];

        t_xy_yy_0_yz[i] = t_x_yyy_0_yz[i] - rab_y[i] * t_x_yy_0_yz[i];

        t_xy_yy_0_yy[i] = t_x_yyy_0_yy[i] - rab_y[i] * t_x_yy_0_yy[i];

        t_xy_yy_0_xz[i] = t_x_yyy_0_xz[i] - rab_y[i] * t_x_yy_0_xz[i];

        t_xy_yy_0_xy[i] = t_x_yyy_0_xy[i] - rab_y[i] * t_x_yy_0_xy[i];

        t_xy_yy_0_xx[i] = t_x_yyy_0_xx[i] - rab_y[i] * t_x_yy_0_xx[i];

        t_xy_xz_0_zz[i] = t_x_xyz_0_zz[i] - rab_y[i] * t_x_xz_0_zz[i];

        t_xy_xz_0_yz[i] = t_x_xyz_0_yz[i] - rab_y[i] * t_x_xz_0_yz[i];

        t_xy_xz_0_yy[i] = t_x_xyz_0_yy[i] - rab_y[i] * t_x_xz_0_yy[i];

        t_xy_xz_0_xz[i] = t_x_xyz_0_xz[i] - rab_y[i] * t_x_xz_0_xz[i];

        t_xy_xz_0_xy[i] = t_x_xyz_0_xy[i] - rab_y[i] * t_x_xz_0_xy[i];

        t_xy_xz_0_xx[i] = t_x_xyz_0_xx[i] - rab_y[i] * t_x_xz_0_xx[i];

        t_xy_xy_0_zz[i] = t_x_xyy_0_zz[i] - rab_y[i] * t_x_xy_0_zz[i];

        t_xy_xy_0_yz[i] = t_x_xyy_0_yz[i] - rab_y[i] * t_x_xy_0_yz[i];

        t_xy_xy_0_yy[i] = t_x_xyy_0_yy[i] - rab_y[i] * t_x_xy_0_yy[i];

        t_xy_xy_0_xz[i] = t_x_xyy_0_xz[i] - rab_y[i] * t_x_xy_0_xz[i];

        t_xy_xy_0_xy[i] = t_x_xyy_0_xy[i] - rab_y[i] * t_x_xy_0_xy[i];

        t_xy_xy_0_xx[i] = t_x_xyy_0_xx[i] - rab_y[i] * t_x_xy_0_xx[i];

        t_xy_xx_0_zz[i] = t_x_xxy_0_zz[i] - rab_y[i] * t_x_xx_0_zz[i];

        t_xy_xx_0_yz[i] = t_x_xxy_0_yz[i] - rab_y[i] * t_x_xx_0_yz[i];

        t_xy_xx_0_yy[i] = t_x_xxy_0_yy[i] - rab_y[i] * t_x_xx_0_yy[i];

        t_xy_xx_0_xz[i] = t_x_xxy_0_xz[i] - rab_y[i] * t_x_xx_0_xz[i];

        t_xy_xx_0_xy[i] = t_x_xxy_0_xy[i] - rab_y[i] * t_x_xx_0_xy[i];

        t_xy_xx_0_xx[i] = t_x_xxy_0_xx[i] - rab_y[i] * t_x_xx_0_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xx_0_xx, t_x_xx_0_xy, t_x_xx_0_xz, t_x_xx_0_yy,\
                           t_x_xx_0_yz, t_x_xx_0_zz, t_x_xxx_0_xx, t_x_xxx_0_xy, t_x_xxx_0_xz,\
                           t_x_xxx_0_yy, t_x_xxx_0_yz, t_x_xxx_0_zz, t_x_xxy_0_xx, t_x_xxy_0_xy,\
                           t_x_xxy_0_xz, t_x_xxy_0_yy, t_x_xxy_0_yz, t_x_xxy_0_zz, t_x_xxz_0_xx,\
                           t_x_xxz_0_xy, t_x_xxz_0_xz, t_x_xxz_0_yy, t_x_xxz_0_yz, t_x_xxz_0_zz,\
                           t_x_xy_0_xx, t_x_xy_0_xy, t_x_xy_0_xz, t_x_xy_0_yy, t_x_xy_0_yz,\
                           t_x_xy_0_zz, t_x_xyy_0_xx, t_x_xyy_0_xy, t_x_xyy_0_xz, t_x_xyy_0_yy,\
                           t_x_xyy_0_yz, t_x_xyy_0_zz, t_x_xyz_0_xx, t_x_xyz_0_xy, t_x_xyz_0_xz,\
                           t_x_xyz_0_yy, t_x_xyz_0_yz, t_x_xyz_0_zz, t_x_xz_0_xx, t_x_xz_0_xy,\
                           t_x_xz_0_xz, t_x_xz_0_yy, t_x_xz_0_yz, t_x_xz_0_zz, t_x_xzz_0_xx,\
                           t_x_xzz_0_xy, t_x_xzz_0_xz, t_x_xzz_0_yy, t_x_xzz_0_yz, t_x_xzz_0_zz,\
                           t_x_yy_0_xx, t_x_yy_0_xy, t_x_yy_0_xz, t_x_yy_0_yy, t_x_yy_0_yz,\
                           t_x_yy_0_zz, t_x_yz_0_xx, t_x_yz_0_xy, t_x_yz_0_xz, t_x_yz_0_yy,\
                           t_x_yz_0_yz, t_x_yz_0_zz, t_x_zz_0_xx, t_x_zz_0_xy, t_x_zz_0_xz,\
                           t_x_zz_0_yy, t_x_zz_0_yz, t_x_zz_0_zz, t_xx_xx_0_xx, t_xx_xx_0_xy,\
                           t_xx_xx_0_xz, t_xx_xx_0_yy, t_xx_xx_0_yz, t_xx_xx_0_zz, t_xx_xy_0_xx,\
                           t_xx_xy_0_xy, t_xx_xy_0_xz, t_xx_xy_0_yy, t_xx_xy_0_yz, t_xx_xy_0_zz,\
                           t_xx_xz_0_xx, t_xx_xz_0_xy, t_xx_xz_0_xz, t_xx_xz_0_yy, t_xx_xz_0_yz,\
                           t_xx_xz_0_zz, t_xx_yy_0_xx, t_xx_yy_0_xy, t_xx_yy_0_xz, t_xx_yy_0_yy,\
                           t_xx_yy_0_yz, t_xx_yy_0_zz, t_xx_yz_0_xx, t_xx_yz_0_xy, t_xx_yz_0_xz,\
                           t_xx_yz_0_yy, t_xx_yz_0_yz, t_xx_yz_0_zz, t_xx_zz_0_xx, t_xx_zz_0_xy,\
                           t_xx_zz_0_xz, t_xx_zz_0_yy, t_xx_zz_0_yz, t_xx_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_zz_0_zz[i] = t_x_xzz_0_zz[i] - rab_x[i] * t_x_zz_0_zz[i];

        t_xx_zz_0_yz[i] = t_x_xzz_0_yz[i] - rab_x[i] * t_x_zz_0_yz[i];

        t_xx_zz_0_yy[i] = t_x_xzz_0_yy[i] - rab_x[i] * t_x_zz_0_yy[i];

        t_xx_zz_0_xz[i] = t_x_xzz_0_xz[i] - rab_x[i] * t_x_zz_0_xz[i];

        t_xx_zz_0_xy[i] = t_x_xzz_0_xy[i] - rab_x[i] * t_x_zz_0_xy[i];

        t_xx_zz_0_xx[i] = t_x_xzz_0_xx[i] - rab_x[i] * t_x_zz_0_xx[i];

        t_xx_yz_0_zz[i] = t_x_xyz_0_zz[i] - rab_x[i] * t_x_yz_0_zz[i];

        t_xx_yz_0_yz[i] = t_x_xyz_0_yz[i] - rab_x[i] * t_x_yz_0_yz[i];

        t_xx_yz_0_yy[i] = t_x_xyz_0_yy[i] - rab_x[i] * t_x_yz_0_yy[i];

        t_xx_yz_0_xz[i] = t_x_xyz_0_xz[i] - rab_x[i] * t_x_yz_0_xz[i];

        t_xx_yz_0_xy[i] = t_x_xyz_0_xy[i] - rab_x[i] * t_x_yz_0_xy[i];

        t_xx_yz_0_xx[i] = t_x_xyz_0_xx[i] - rab_x[i] * t_x_yz_0_xx[i];

        t_xx_yy_0_zz[i] = t_x_xyy_0_zz[i] - rab_x[i] * t_x_yy_0_zz[i];

        t_xx_yy_0_yz[i] = t_x_xyy_0_yz[i] - rab_x[i] * t_x_yy_0_yz[i];

        t_xx_yy_0_yy[i] = t_x_xyy_0_yy[i] - rab_x[i] * t_x_yy_0_yy[i];

        t_xx_yy_0_xz[i] = t_x_xyy_0_xz[i] - rab_x[i] * t_x_yy_0_xz[i];

        t_xx_yy_0_xy[i] = t_x_xyy_0_xy[i] - rab_x[i] * t_x_yy_0_xy[i];

        t_xx_yy_0_xx[i] = t_x_xyy_0_xx[i] - rab_x[i] * t_x_yy_0_xx[i];

        t_xx_xz_0_zz[i] = t_x_xxz_0_zz[i] - rab_x[i] * t_x_xz_0_zz[i];

        t_xx_xz_0_yz[i] = t_x_xxz_0_yz[i] - rab_x[i] * t_x_xz_0_yz[i];

        t_xx_xz_0_yy[i] = t_x_xxz_0_yy[i] - rab_x[i] * t_x_xz_0_yy[i];

        t_xx_xz_0_xz[i] = t_x_xxz_0_xz[i] - rab_x[i] * t_x_xz_0_xz[i];

        t_xx_xz_0_xy[i] = t_x_xxz_0_xy[i] - rab_x[i] * t_x_xz_0_xy[i];

        t_xx_xz_0_xx[i] = t_x_xxz_0_xx[i] - rab_x[i] * t_x_xz_0_xx[i];

        t_xx_xy_0_zz[i] = t_x_xxy_0_zz[i] - rab_x[i] * t_x_xy_0_zz[i];

        t_xx_xy_0_yz[i] = t_x_xxy_0_yz[i] - rab_x[i] * t_x_xy_0_yz[i];

        t_xx_xy_0_yy[i] = t_x_xxy_0_yy[i] - rab_x[i] * t_x_xy_0_yy[i];

        t_xx_xy_0_xz[i] = t_x_xxy_0_xz[i] - rab_x[i] * t_x_xy_0_xz[i];

        t_xx_xy_0_xy[i] = t_x_xxy_0_xy[i] - rab_x[i] * t_x_xy_0_xy[i];

        t_xx_xy_0_xx[i] = t_x_xxy_0_xx[i] - rab_x[i] * t_x_xy_0_xx[i];

        t_xx_xx_0_zz[i] = t_x_xxx_0_zz[i] - rab_x[i] * t_x_xx_0_zz[i];

        t_xx_xx_0_yz[i] = t_x_xxx_0_yz[i] - rab_x[i] * t_x_xx_0_yz[i];

        t_xx_xx_0_yy[i] = t_x_xxx_0_yy[i] - rab_x[i] * t_x_xx_0_yy[i];

        t_xx_xx_0_xz[i] = t_x_xxx_0_xz[i] - rab_x[i] * t_x_xx_0_xz[i];

        t_xx_xx_0_xy[i] = t_x_xxx_0_xy[i] - rab_x[i] * t_x_xx_0_xy[i];

        t_xx_xx_0_xx[i] = t_x_xxx_0_xx[i] - rab_x[i] * t_x_xx_0_xx[i];
    }
}


} // derirec namespace
