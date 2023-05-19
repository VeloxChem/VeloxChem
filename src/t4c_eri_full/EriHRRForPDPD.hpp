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
compHostHRRForPDPD_V0(      BufferHostXY<T>&      intsBufferPDPD,
                      const BufferHostX<int32_t>& intsIndexesPDPD,
                      const BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PDPD) integral components

    t_z_zz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(0));

    t_z_zz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(1));

    t_z_zz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(2));

    t_z_zz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(3));

    t_z_zz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(4));

    t_z_zz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(5));

    t_z_zz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(6));

    t_z_zz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(7));

    t_z_zz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(8));

    t_z_zz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(9));

    t_z_zz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(10));

    t_z_zz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(11));

    t_z_zz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(12));

    t_z_zz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(13));

    t_z_zz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(14));

    t_z_zz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(15));

    t_z_zz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(16));

    t_z_zz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(17));

    t_z_yz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(18));

    t_z_yz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(19));

    t_z_yz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(20));

    t_z_yz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(21));

    t_z_yz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(22));

    t_z_yz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(23));

    t_z_yz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(24));

    t_z_yz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(25));

    t_z_yz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(26));

    t_z_yz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(27));

    t_z_yz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(28));

    t_z_yz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(29));

    t_z_yz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(30));

    t_z_yz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(31));

    t_z_yz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(32));

    t_z_yz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(33));

    t_z_yz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(34));

    t_z_yz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(35));

    t_z_yy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(36));

    t_z_yy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(37));

    t_z_yy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(38));

    t_z_yy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(39));

    t_z_yy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(40));

    t_z_yy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(41));

    t_z_yy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(42));

    t_z_yy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(43));

    t_z_yy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(44));

    t_z_yy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(45));

    t_z_yy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(46));

    t_z_yy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(47));

    t_z_yy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(48));

    t_z_yy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(49));

    t_z_yy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(50));

    t_z_yy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(51));

    t_z_yy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(52));

    t_z_yy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(53));

    t_z_xz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(54));

    t_z_xz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(55));

    t_z_xz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(56));

    t_z_xz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(57));

    t_z_xz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(58));

    t_z_xz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(59));

    t_z_xz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(60));

    t_z_xz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(61));

    t_z_xz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(62));

    t_z_xz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(63));

    t_z_xz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(64));

    t_z_xz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(65));

    t_z_xz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(66));

    t_z_xz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(67));

    t_z_xz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(68));

    t_z_xz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(69));

    t_z_xz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(70));

    t_z_xz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(71));

    t_z_xy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(72));

    t_z_xy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(73));

    t_z_xy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(74));

    t_z_xy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(75));

    t_z_xy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(76));

    t_z_xy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(77));

    t_z_xy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(78));

    t_z_xy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(79));

    t_z_xy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(80));

    t_z_xy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(81));

    t_z_xy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(82));

    t_z_xy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(83));

    t_z_xy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(84));

    t_z_xy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(85));

    t_z_xy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(86));

    t_z_xy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(87));

    t_z_xy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(88));

    t_z_xy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(89));

    t_z_xx_z_zz = intsBufferPDPD.data(intsIndexesPDPD(90));

    t_z_xx_z_yz = intsBufferPDPD.data(intsIndexesPDPD(91));

    t_z_xx_z_yy = intsBufferPDPD.data(intsIndexesPDPD(92));

    t_z_xx_z_xz = intsBufferPDPD.data(intsIndexesPDPD(93));

    t_z_xx_z_xy = intsBufferPDPD.data(intsIndexesPDPD(94));

    t_z_xx_z_xx = intsBufferPDPD.data(intsIndexesPDPD(95));

    t_z_xx_y_zz = intsBufferPDPD.data(intsIndexesPDPD(96));

    t_z_xx_y_yz = intsBufferPDPD.data(intsIndexesPDPD(97));

    t_z_xx_y_yy = intsBufferPDPD.data(intsIndexesPDPD(98));

    t_z_xx_y_xz = intsBufferPDPD.data(intsIndexesPDPD(99));

    t_z_xx_y_xy = intsBufferPDPD.data(intsIndexesPDPD(100));

    t_z_xx_y_xx = intsBufferPDPD.data(intsIndexesPDPD(101));

    t_z_xx_x_zz = intsBufferPDPD.data(intsIndexesPDPD(102));

    t_z_xx_x_yz = intsBufferPDPD.data(intsIndexesPDPD(103));

    t_z_xx_x_yy = intsBufferPDPD.data(intsIndexesPDPD(104));

    t_z_xx_x_xz = intsBufferPDPD.data(intsIndexesPDPD(105));

    t_z_xx_x_xy = intsBufferPDPD.data(intsIndexesPDPD(106));

    t_z_xx_x_xx = intsBufferPDPD.data(intsIndexesPDPD(107));

    t_y_zz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(108));

    t_y_zz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(109));

    t_y_zz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(110));

    t_y_zz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(111));

    t_y_zz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(112));

    t_y_zz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(113));

    t_y_zz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(114));

    t_y_zz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(115));

    t_y_zz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(116));

    t_y_zz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(117));

    t_y_zz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(118));

    t_y_zz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(119));

    t_y_zz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(120));

    t_y_zz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(121));

    t_y_zz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(122));

    t_y_zz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(123));

    t_y_zz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(124));

    t_y_zz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(125));

    t_y_yz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(126));

    t_y_yz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(127));

    t_y_yz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(128));

    t_y_yz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(129));

    t_y_yz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(130));

    t_y_yz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(131));

    t_y_yz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(132));

    t_y_yz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(133));

    t_y_yz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(134));

    t_y_yz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(135));

    t_y_yz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(136));

    t_y_yz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(137));

    t_y_yz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(138));

    t_y_yz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(139));

    t_y_yz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(140));

    t_y_yz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(141));

    t_y_yz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(142));

    t_y_yz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(143));

    t_y_yy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(144));

    t_y_yy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(145));

    t_y_yy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(146));

    t_y_yy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(147));

    t_y_yy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(148));

    t_y_yy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(149));

    t_y_yy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(150));

    t_y_yy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(151));

    t_y_yy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(152));

    t_y_yy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(153));

    t_y_yy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(154));

    t_y_yy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(155));

    t_y_yy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(156));

    t_y_yy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(157));

    t_y_yy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(158));

    t_y_yy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(159));

    t_y_yy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(160));

    t_y_yy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(161));

    t_y_xz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(162));

    t_y_xz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(163));

    t_y_xz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(164));

    t_y_xz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(165));

    t_y_xz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(166));

    t_y_xz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(167));

    t_y_xz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(168));

    t_y_xz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(169));

    t_y_xz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(170));

    t_y_xz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(171));

    t_y_xz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(172));

    t_y_xz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(173));

    t_y_xz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(174));

    t_y_xz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(175));

    t_y_xz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(176));

    t_y_xz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(177));

    t_y_xz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(178));

    t_y_xz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(179));

    t_y_xy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(180));

    t_y_xy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(181));

    t_y_xy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(182));

    t_y_xy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(183));

    t_y_xy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(184));

    t_y_xy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(185));

    t_y_xy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(186));

    t_y_xy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(187));

    t_y_xy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(188));

    t_y_xy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(189));

    t_y_xy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(190));

    t_y_xy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(191));

    t_y_xy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(192));

    t_y_xy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(193));

    t_y_xy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(194));

    t_y_xy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(195));

    t_y_xy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(196));

    t_y_xy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(197));

    t_y_xx_z_zz = intsBufferPDPD.data(intsIndexesPDPD(198));

    t_y_xx_z_yz = intsBufferPDPD.data(intsIndexesPDPD(199));

    t_y_xx_z_yy = intsBufferPDPD.data(intsIndexesPDPD(200));

    t_y_xx_z_xz = intsBufferPDPD.data(intsIndexesPDPD(201));

    t_y_xx_z_xy = intsBufferPDPD.data(intsIndexesPDPD(202));

    t_y_xx_z_xx = intsBufferPDPD.data(intsIndexesPDPD(203));

    t_y_xx_y_zz = intsBufferPDPD.data(intsIndexesPDPD(204));

    t_y_xx_y_yz = intsBufferPDPD.data(intsIndexesPDPD(205));

    t_y_xx_y_yy = intsBufferPDPD.data(intsIndexesPDPD(206));

    t_y_xx_y_xz = intsBufferPDPD.data(intsIndexesPDPD(207));

    t_y_xx_y_xy = intsBufferPDPD.data(intsIndexesPDPD(208));

    t_y_xx_y_xx = intsBufferPDPD.data(intsIndexesPDPD(209));

    t_y_xx_x_zz = intsBufferPDPD.data(intsIndexesPDPD(210));

    t_y_xx_x_yz = intsBufferPDPD.data(intsIndexesPDPD(211));

    t_y_xx_x_yy = intsBufferPDPD.data(intsIndexesPDPD(212));

    t_y_xx_x_xz = intsBufferPDPD.data(intsIndexesPDPD(213));

    t_y_xx_x_xy = intsBufferPDPD.data(intsIndexesPDPD(214));

    t_y_xx_x_xx = intsBufferPDPD.data(intsIndexesPDPD(215));

    t_x_zz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(216));

    t_x_zz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(217));

    t_x_zz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(218));

    t_x_zz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(219));

    t_x_zz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(220));

    t_x_zz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(221));

    t_x_zz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(222));

    t_x_zz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(223));

    t_x_zz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(224));

    t_x_zz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(225));

    t_x_zz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(226));

    t_x_zz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(227));

    t_x_zz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(228));

    t_x_zz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(229));

    t_x_zz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(230));

    t_x_zz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(231));

    t_x_zz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(232));

    t_x_zz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(233));

    t_x_yz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(234));

    t_x_yz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(235));

    t_x_yz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(236));

    t_x_yz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(237));

    t_x_yz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(238));

    t_x_yz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(239));

    t_x_yz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(240));

    t_x_yz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(241));

    t_x_yz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(242));

    t_x_yz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(243));

    t_x_yz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(244));

    t_x_yz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(245));

    t_x_yz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(246));

    t_x_yz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(247));

    t_x_yz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(248));

    t_x_yz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(249));

    t_x_yz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(250));

    t_x_yz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(251));

    t_x_yy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(252));

    t_x_yy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(253));

    t_x_yy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(254));

    t_x_yy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(255));

    t_x_yy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(256));

    t_x_yy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(257));

    t_x_yy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(258));

    t_x_yy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(259));

    t_x_yy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(260));

    t_x_yy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(261));

    t_x_yy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(262));

    t_x_yy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(263));

    t_x_yy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(264));

    t_x_yy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(265));

    t_x_yy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(266));

    t_x_yy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(267));

    t_x_yy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(268));

    t_x_yy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(269));

    t_x_xz_z_zz = intsBufferPDPD.data(intsIndexesPDPD(270));

    t_x_xz_z_yz = intsBufferPDPD.data(intsIndexesPDPD(271));

    t_x_xz_z_yy = intsBufferPDPD.data(intsIndexesPDPD(272));

    t_x_xz_z_xz = intsBufferPDPD.data(intsIndexesPDPD(273));

    t_x_xz_z_xy = intsBufferPDPD.data(intsIndexesPDPD(274));

    t_x_xz_z_xx = intsBufferPDPD.data(intsIndexesPDPD(275));

    t_x_xz_y_zz = intsBufferPDPD.data(intsIndexesPDPD(276));

    t_x_xz_y_yz = intsBufferPDPD.data(intsIndexesPDPD(277));

    t_x_xz_y_yy = intsBufferPDPD.data(intsIndexesPDPD(278));

    t_x_xz_y_xz = intsBufferPDPD.data(intsIndexesPDPD(279));

    t_x_xz_y_xy = intsBufferPDPD.data(intsIndexesPDPD(280));

    t_x_xz_y_xx = intsBufferPDPD.data(intsIndexesPDPD(281));

    t_x_xz_x_zz = intsBufferPDPD.data(intsIndexesPDPD(282));

    t_x_xz_x_yz = intsBufferPDPD.data(intsIndexesPDPD(283));

    t_x_xz_x_yy = intsBufferPDPD.data(intsIndexesPDPD(284));

    t_x_xz_x_xz = intsBufferPDPD.data(intsIndexesPDPD(285));

    t_x_xz_x_xy = intsBufferPDPD.data(intsIndexesPDPD(286));

    t_x_xz_x_xx = intsBufferPDPD.data(intsIndexesPDPD(287));

    t_x_xy_z_zz = intsBufferPDPD.data(intsIndexesPDPD(288));

    t_x_xy_z_yz = intsBufferPDPD.data(intsIndexesPDPD(289));

    t_x_xy_z_yy = intsBufferPDPD.data(intsIndexesPDPD(290));

    t_x_xy_z_xz = intsBufferPDPD.data(intsIndexesPDPD(291));

    t_x_xy_z_xy = intsBufferPDPD.data(intsIndexesPDPD(292));

    t_x_xy_z_xx = intsBufferPDPD.data(intsIndexesPDPD(293));

    t_x_xy_y_zz = intsBufferPDPD.data(intsIndexesPDPD(294));

    t_x_xy_y_yz = intsBufferPDPD.data(intsIndexesPDPD(295));

    t_x_xy_y_yy = intsBufferPDPD.data(intsIndexesPDPD(296));

    t_x_xy_y_xz = intsBufferPDPD.data(intsIndexesPDPD(297));

    t_x_xy_y_xy = intsBufferPDPD.data(intsIndexesPDPD(298));

    t_x_xy_y_xx = intsBufferPDPD.data(intsIndexesPDPD(299));

    t_x_xy_x_zz = intsBufferPDPD.data(intsIndexesPDPD(300));

    t_x_xy_x_yz = intsBufferPDPD.data(intsIndexesPDPD(301));

    t_x_xy_x_yy = intsBufferPDPD.data(intsIndexesPDPD(302));

    t_x_xy_x_xz = intsBufferPDPD.data(intsIndexesPDPD(303));

    t_x_xy_x_xy = intsBufferPDPD.data(intsIndexesPDPD(304));

    t_x_xy_x_xx = intsBufferPDPD.data(intsIndexesPDPD(305));

    t_x_xx_z_zz = intsBufferPDPD.data(intsIndexesPDPD(306));

    t_x_xx_z_yz = intsBufferPDPD.data(intsIndexesPDPD(307));

    t_x_xx_z_yy = intsBufferPDPD.data(intsIndexesPDPD(308));

    t_x_xx_z_xz = intsBufferPDPD.data(intsIndexesPDPD(309));

    t_x_xx_z_xy = intsBufferPDPD.data(intsIndexesPDPD(310));

    t_x_xx_z_xx = intsBufferPDPD.data(intsIndexesPDPD(311));

    t_x_xx_y_zz = intsBufferPDPD.data(intsIndexesPDPD(312));

    t_x_xx_y_yz = intsBufferPDPD.data(intsIndexesPDPD(313));

    t_x_xx_y_yy = intsBufferPDPD.data(intsIndexesPDPD(314));

    t_x_xx_y_xz = intsBufferPDPD.data(intsIndexesPDPD(315));

    t_x_xx_y_xy = intsBufferPDPD.data(intsIndexesPDPD(316));

    t_x_xx_y_xx = intsBufferPDPD.data(intsIndexesPDPD(317));

    t_x_xx_x_zz = intsBufferPDPD.data(intsIndexesPDPD(318));

    t_x_xx_x_yz = intsBufferPDPD.data(intsIndexesPDPD(319));

    t_x_xx_x_yy = intsBufferPDPD.data(intsIndexesPDPD(320));

    t_x_xx_x_xz = intsBufferPDPD.data(intsIndexesPDPD(321));

    t_x_xx_x_xy = intsBufferPDPD.data(intsIndexesPDPD(322));

    t_x_xx_x_xx = intsBufferPDPD.data(intsIndexesPDPD(323));

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

    #pragma omp simd align(rab_z, t_0_yz_x_xx, t_0_yz_x_xy, t_0_yz_x_xz, t_0_yz_x_yy,\
                           t_0_yz_x_yz, t_0_yz_x_zz, t_0_yz_y_xx, t_0_yz_y_xy, t_0_yz_y_xz,\
                           t_0_yz_y_yy, t_0_yz_y_yz, t_0_yz_y_zz, t_0_yz_z_xx, t_0_yz_z_xy,\
                           t_0_yz_z_xz, t_0_yz_z_yy, t_0_yz_z_yz, t_0_yz_z_zz, t_0_yzz_x_xx,\
                           t_0_yzz_x_xy, t_0_yzz_x_xz, t_0_yzz_x_yy, t_0_yzz_x_yz, t_0_yzz_x_zz,\
                           t_0_yzz_y_xx, t_0_yzz_y_xy, t_0_yzz_y_xz, t_0_yzz_y_yy, t_0_yzz_y_yz,\
                           t_0_yzz_y_zz, t_0_yzz_z_xx, t_0_yzz_z_xy, t_0_yzz_z_xz, t_0_yzz_z_yy,\
                           t_0_yzz_z_yz, t_0_yzz_z_zz, t_0_zz_x_xx, t_0_zz_x_xy, t_0_zz_x_xz,\
                           t_0_zz_x_yy, t_0_zz_x_yz, t_0_zz_x_zz, t_0_zz_y_xx, t_0_zz_y_xy,\
                           t_0_zz_y_xz, t_0_zz_y_yy, t_0_zz_y_yz, t_0_zz_y_zz, t_0_zz_z_xx,\
                           t_0_zz_z_xy, t_0_zz_z_xz, t_0_zz_z_yy, t_0_zz_z_yz, t_0_zz_z_zz,\
                           t_0_zzz_x_xx, t_0_zzz_x_xy, t_0_zzz_x_xz, t_0_zzz_x_yy, t_0_zzz_x_yz,\
                           t_0_zzz_x_zz, t_0_zzz_y_xx, t_0_zzz_y_xy, t_0_zzz_y_xz, t_0_zzz_y_yy,\
                           t_0_zzz_y_yz, t_0_zzz_y_zz, t_0_zzz_z_xx, t_0_zzz_z_xy, t_0_zzz_z_xz,\
                           t_0_zzz_z_yy, t_0_zzz_z_yz, t_0_zzz_z_zz, t_z_yz_x_xx, t_z_yz_x_xy,\
                           t_z_yz_x_xz, t_z_yz_x_yy, t_z_yz_x_yz, t_z_yz_x_zz, t_z_yz_y_xx,\
                           t_z_yz_y_xy, t_z_yz_y_xz, t_z_yz_y_yy, t_z_yz_y_yz, t_z_yz_y_zz,\
                           t_z_yz_z_xx, t_z_yz_z_xy, t_z_yz_z_xz, t_z_yz_z_yy, t_z_yz_z_yz,\
                           t_z_yz_z_zz, t_z_zz_x_xx, t_z_zz_x_xy, t_z_zz_x_xz, t_z_zz_x_yy,\
                           t_z_zz_x_yz, t_z_zz_x_zz, t_z_zz_y_xx, t_z_zz_y_xy, t_z_zz_y_xz,\
                           t_z_zz_y_yy, t_z_zz_y_yz, t_z_zz_y_zz, t_z_zz_z_xx, t_z_zz_z_xy,\
                           t_z_zz_z_xz, t_z_zz_z_yy, t_z_zz_z_yz, t_z_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_z_zz[i] = t_0_zzz_z_zz[i] - rab_z[i] * t_0_zz_z_zz[i];

        t_z_zz_z_yz[i] = t_0_zzz_z_yz[i] - rab_z[i] * t_0_zz_z_yz[i];

        t_z_zz_z_yy[i] = t_0_zzz_z_yy[i] - rab_z[i] * t_0_zz_z_yy[i];

        t_z_zz_z_xz[i] = t_0_zzz_z_xz[i] - rab_z[i] * t_0_zz_z_xz[i];

        t_z_zz_z_xy[i] = t_0_zzz_z_xy[i] - rab_z[i] * t_0_zz_z_xy[i];

        t_z_zz_z_xx[i] = t_0_zzz_z_xx[i] - rab_z[i] * t_0_zz_z_xx[i];

        t_z_zz_y_zz[i] = t_0_zzz_y_zz[i] - rab_z[i] * t_0_zz_y_zz[i];

        t_z_zz_y_yz[i] = t_0_zzz_y_yz[i] - rab_z[i] * t_0_zz_y_yz[i];

        t_z_zz_y_yy[i] = t_0_zzz_y_yy[i] - rab_z[i] * t_0_zz_y_yy[i];

        t_z_zz_y_xz[i] = t_0_zzz_y_xz[i] - rab_z[i] * t_0_zz_y_xz[i];

        t_z_zz_y_xy[i] = t_0_zzz_y_xy[i] - rab_z[i] * t_0_zz_y_xy[i];

        t_z_zz_y_xx[i] = t_0_zzz_y_xx[i] - rab_z[i] * t_0_zz_y_xx[i];

        t_z_zz_x_zz[i] = t_0_zzz_x_zz[i] - rab_z[i] * t_0_zz_x_zz[i];

        t_z_zz_x_yz[i] = t_0_zzz_x_yz[i] - rab_z[i] * t_0_zz_x_yz[i];

        t_z_zz_x_yy[i] = t_0_zzz_x_yy[i] - rab_z[i] * t_0_zz_x_yy[i];

        t_z_zz_x_xz[i] = t_0_zzz_x_xz[i] - rab_z[i] * t_0_zz_x_xz[i];

        t_z_zz_x_xy[i] = t_0_zzz_x_xy[i] - rab_z[i] * t_0_zz_x_xy[i];

        t_z_zz_x_xx[i] = t_0_zzz_x_xx[i] - rab_z[i] * t_0_zz_x_xx[i];

        t_z_yz_z_zz[i] = t_0_yzz_z_zz[i] - rab_z[i] * t_0_yz_z_zz[i];

        t_z_yz_z_yz[i] = t_0_yzz_z_yz[i] - rab_z[i] * t_0_yz_z_yz[i];

        t_z_yz_z_yy[i] = t_0_yzz_z_yy[i] - rab_z[i] * t_0_yz_z_yy[i];

        t_z_yz_z_xz[i] = t_0_yzz_z_xz[i] - rab_z[i] * t_0_yz_z_xz[i];

        t_z_yz_z_xy[i] = t_0_yzz_z_xy[i] - rab_z[i] * t_0_yz_z_xy[i];

        t_z_yz_z_xx[i] = t_0_yzz_z_xx[i] - rab_z[i] * t_0_yz_z_xx[i];

        t_z_yz_y_zz[i] = t_0_yzz_y_zz[i] - rab_z[i] * t_0_yz_y_zz[i];

        t_z_yz_y_yz[i] = t_0_yzz_y_yz[i] - rab_z[i] * t_0_yz_y_yz[i];

        t_z_yz_y_yy[i] = t_0_yzz_y_yy[i] - rab_z[i] * t_0_yz_y_yy[i];

        t_z_yz_y_xz[i] = t_0_yzz_y_xz[i] - rab_z[i] * t_0_yz_y_xz[i];

        t_z_yz_y_xy[i] = t_0_yzz_y_xy[i] - rab_z[i] * t_0_yz_y_xy[i];

        t_z_yz_y_xx[i] = t_0_yzz_y_xx[i] - rab_z[i] * t_0_yz_y_xx[i];

        t_z_yz_x_zz[i] = t_0_yzz_x_zz[i] - rab_z[i] * t_0_yz_x_zz[i];

        t_z_yz_x_yz[i] = t_0_yzz_x_yz[i] - rab_z[i] * t_0_yz_x_yz[i];

        t_z_yz_x_yy[i] = t_0_yzz_x_yy[i] - rab_z[i] * t_0_yz_x_yy[i];

        t_z_yz_x_xz[i] = t_0_yzz_x_xz[i] - rab_z[i] * t_0_yz_x_xz[i];

        t_z_yz_x_xy[i] = t_0_yzz_x_xy[i] - rab_z[i] * t_0_yz_x_xy[i];

        t_z_yz_x_xx[i] = t_0_yzz_x_xx[i] - rab_z[i] * t_0_yz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xz_x_xx, t_0_xz_x_xy, t_0_xz_x_xz, t_0_xz_x_yy,\
                           t_0_xz_x_yz, t_0_xz_x_zz, t_0_xz_y_xx, t_0_xz_y_xy, t_0_xz_y_xz,\
                           t_0_xz_y_yy, t_0_xz_y_yz, t_0_xz_y_zz, t_0_xz_z_xx, t_0_xz_z_xy,\
                           t_0_xz_z_xz, t_0_xz_z_yy, t_0_xz_z_yz, t_0_xz_z_zz, t_0_xzz_x_xx,\
                           t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy, t_0_xzz_x_yz, t_0_xzz_x_zz,\
                           t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz, t_0_xzz_y_yy, t_0_xzz_y_yz,\
                           t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy, t_0_xzz_z_xz, t_0_xzz_z_yy,\
                           t_0_xzz_z_yz, t_0_xzz_z_zz, t_0_yy_x_xx, t_0_yy_x_xy, t_0_yy_x_xz,\
                           t_0_yy_x_yy, t_0_yy_x_yz, t_0_yy_x_zz, t_0_yy_y_xx, t_0_yy_y_xy,\
                           t_0_yy_y_xz, t_0_yy_y_yy, t_0_yy_y_yz, t_0_yy_y_zz, t_0_yy_z_xx,\
                           t_0_yy_z_xy, t_0_yy_z_xz, t_0_yy_z_yy, t_0_yy_z_yz, t_0_yy_z_zz,\
                           t_0_yyz_x_xx, t_0_yyz_x_xy, t_0_yyz_x_xz, t_0_yyz_x_yy, t_0_yyz_x_yz,\
                           t_0_yyz_x_zz, t_0_yyz_y_xx, t_0_yyz_y_xy, t_0_yyz_y_xz, t_0_yyz_y_yy,\
                           t_0_yyz_y_yz, t_0_yyz_y_zz, t_0_yyz_z_xx, t_0_yyz_z_xy, t_0_yyz_z_xz,\
                           t_0_yyz_z_yy, t_0_yyz_z_yz, t_0_yyz_z_zz, t_z_xz_x_xx, t_z_xz_x_xy,\
                           t_z_xz_x_xz, t_z_xz_x_yy, t_z_xz_x_yz, t_z_xz_x_zz, t_z_xz_y_xx,\
                           t_z_xz_y_xy, t_z_xz_y_xz, t_z_xz_y_yy, t_z_xz_y_yz, t_z_xz_y_zz,\
                           t_z_xz_z_xx, t_z_xz_z_xy, t_z_xz_z_xz, t_z_xz_z_yy, t_z_xz_z_yz,\
                           t_z_xz_z_zz, t_z_yy_x_xx, t_z_yy_x_xy, t_z_yy_x_xz, t_z_yy_x_yy,\
                           t_z_yy_x_yz, t_z_yy_x_zz, t_z_yy_y_xx, t_z_yy_y_xy, t_z_yy_y_xz,\
                           t_z_yy_y_yy, t_z_yy_y_yz, t_z_yy_y_zz, t_z_yy_z_xx, t_z_yy_z_xy,\
                           t_z_yy_z_xz, t_z_yy_z_yy, t_z_yy_z_yz, t_z_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_yy_z_zz[i] = t_0_yyz_z_zz[i] - rab_z[i] * t_0_yy_z_zz[i];

        t_z_yy_z_yz[i] = t_0_yyz_z_yz[i] - rab_z[i] * t_0_yy_z_yz[i];

        t_z_yy_z_yy[i] = t_0_yyz_z_yy[i] - rab_z[i] * t_0_yy_z_yy[i];

        t_z_yy_z_xz[i] = t_0_yyz_z_xz[i] - rab_z[i] * t_0_yy_z_xz[i];

        t_z_yy_z_xy[i] = t_0_yyz_z_xy[i] - rab_z[i] * t_0_yy_z_xy[i];

        t_z_yy_z_xx[i] = t_0_yyz_z_xx[i] - rab_z[i] * t_0_yy_z_xx[i];

        t_z_yy_y_zz[i] = t_0_yyz_y_zz[i] - rab_z[i] * t_0_yy_y_zz[i];

        t_z_yy_y_yz[i] = t_0_yyz_y_yz[i] - rab_z[i] * t_0_yy_y_yz[i];

        t_z_yy_y_yy[i] = t_0_yyz_y_yy[i] - rab_z[i] * t_0_yy_y_yy[i];

        t_z_yy_y_xz[i] = t_0_yyz_y_xz[i] - rab_z[i] * t_0_yy_y_xz[i];

        t_z_yy_y_xy[i] = t_0_yyz_y_xy[i] - rab_z[i] * t_0_yy_y_xy[i];

        t_z_yy_y_xx[i] = t_0_yyz_y_xx[i] - rab_z[i] * t_0_yy_y_xx[i];

        t_z_yy_x_zz[i] = t_0_yyz_x_zz[i] - rab_z[i] * t_0_yy_x_zz[i];

        t_z_yy_x_yz[i] = t_0_yyz_x_yz[i] - rab_z[i] * t_0_yy_x_yz[i];

        t_z_yy_x_yy[i] = t_0_yyz_x_yy[i] - rab_z[i] * t_0_yy_x_yy[i];

        t_z_yy_x_xz[i] = t_0_yyz_x_xz[i] - rab_z[i] * t_0_yy_x_xz[i];

        t_z_yy_x_xy[i] = t_0_yyz_x_xy[i] - rab_z[i] * t_0_yy_x_xy[i];

        t_z_yy_x_xx[i] = t_0_yyz_x_xx[i] - rab_z[i] * t_0_yy_x_xx[i];

        t_z_xz_z_zz[i] = t_0_xzz_z_zz[i] - rab_z[i] * t_0_xz_z_zz[i];

        t_z_xz_z_yz[i] = t_0_xzz_z_yz[i] - rab_z[i] * t_0_xz_z_yz[i];

        t_z_xz_z_yy[i] = t_0_xzz_z_yy[i] - rab_z[i] * t_0_xz_z_yy[i];

        t_z_xz_z_xz[i] = t_0_xzz_z_xz[i] - rab_z[i] * t_0_xz_z_xz[i];

        t_z_xz_z_xy[i] = t_0_xzz_z_xy[i] - rab_z[i] * t_0_xz_z_xy[i];

        t_z_xz_z_xx[i] = t_0_xzz_z_xx[i] - rab_z[i] * t_0_xz_z_xx[i];

        t_z_xz_y_zz[i] = t_0_xzz_y_zz[i] - rab_z[i] * t_0_xz_y_zz[i];

        t_z_xz_y_yz[i] = t_0_xzz_y_yz[i] - rab_z[i] * t_0_xz_y_yz[i];

        t_z_xz_y_yy[i] = t_0_xzz_y_yy[i] - rab_z[i] * t_0_xz_y_yy[i];

        t_z_xz_y_xz[i] = t_0_xzz_y_xz[i] - rab_z[i] * t_0_xz_y_xz[i];

        t_z_xz_y_xy[i] = t_0_xzz_y_xy[i] - rab_z[i] * t_0_xz_y_xy[i];

        t_z_xz_y_xx[i] = t_0_xzz_y_xx[i] - rab_z[i] * t_0_xz_y_xx[i];

        t_z_xz_x_zz[i] = t_0_xzz_x_zz[i] - rab_z[i] * t_0_xz_x_zz[i];

        t_z_xz_x_yz[i] = t_0_xzz_x_yz[i] - rab_z[i] * t_0_xz_x_yz[i];

        t_z_xz_x_yy[i] = t_0_xzz_x_yy[i] - rab_z[i] * t_0_xz_x_yy[i];

        t_z_xz_x_xz[i] = t_0_xzz_x_xz[i] - rab_z[i] * t_0_xz_x_xz[i];

        t_z_xz_x_xy[i] = t_0_xzz_x_xy[i] - rab_z[i] * t_0_xz_x_xy[i];

        t_z_xz_x_xx[i] = t_0_xzz_x_xx[i] - rab_z[i] * t_0_xz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_xx_x_xx, t_0_xx_x_xy, t_0_xx_x_xz, t_0_xx_x_yy,\
                           t_0_xx_x_yz, t_0_xx_x_zz, t_0_xx_y_xx, t_0_xx_y_xy, t_0_xx_y_xz,\
                           t_0_xx_y_yy, t_0_xx_y_yz, t_0_xx_y_zz, t_0_xx_z_xx, t_0_xx_z_xy,\
                           t_0_xx_z_xz, t_0_xx_z_yy, t_0_xx_z_yz, t_0_xx_z_zz, t_0_xxz_x_xx,\
                           t_0_xxz_x_xy, t_0_xxz_x_xz, t_0_xxz_x_yy, t_0_xxz_x_yz, t_0_xxz_x_zz,\
                           t_0_xxz_y_xx, t_0_xxz_y_xy, t_0_xxz_y_xz, t_0_xxz_y_yy, t_0_xxz_y_yz,\
                           t_0_xxz_y_zz, t_0_xxz_z_xx, t_0_xxz_z_xy, t_0_xxz_z_xz, t_0_xxz_z_yy,\
                           t_0_xxz_z_yz, t_0_xxz_z_zz, t_0_xy_x_xx, t_0_xy_x_xy, t_0_xy_x_xz,\
                           t_0_xy_x_yy, t_0_xy_x_yz, t_0_xy_x_zz, t_0_xy_y_xx, t_0_xy_y_xy,\
                           t_0_xy_y_xz, t_0_xy_y_yy, t_0_xy_y_yz, t_0_xy_y_zz, t_0_xy_z_xx,\
                           t_0_xy_z_xy, t_0_xy_z_xz, t_0_xy_z_yy, t_0_xy_z_yz, t_0_xy_z_zz,\
                           t_0_xyz_x_xx, t_0_xyz_x_xy, t_0_xyz_x_xz, t_0_xyz_x_yy, t_0_xyz_x_yz,\
                           t_0_xyz_x_zz, t_0_xyz_y_xx, t_0_xyz_y_xy, t_0_xyz_y_xz, t_0_xyz_y_yy,\
                           t_0_xyz_y_yz, t_0_xyz_y_zz, t_0_xyz_z_xx, t_0_xyz_z_xy, t_0_xyz_z_xz,\
                           t_0_xyz_z_yy, t_0_xyz_z_yz, t_0_xyz_z_zz, t_z_xx_x_xx, t_z_xx_x_xy,\
                           t_z_xx_x_xz, t_z_xx_x_yy, t_z_xx_x_yz, t_z_xx_x_zz, t_z_xx_y_xx,\
                           t_z_xx_y_xy, t_z_xx_y_xz, t_z_xx_y_yy, t_z_xx_y_yz, t_z_xx_y_zz,\
                           t_z_xx_z_xx, t_z_xx_z_xy, t_z_xx_z_xz, t_z_xx_z_yy, t_z_xx_z_yz,\
                           t_z_xx_z_zz, t_z_xy_x_xx, t_z_xy_x_xy, t_z_xy_x_xz, t_z_xy_x_yy,\
                           t_z_xy_x_yz, t_z_xy_x_zz, t_z_xy_y_xx, t_z_xy_y_xy, t_z_xy_y_xz,\
                           t_z_xy_y_yy, t_z_xy_y_yz, t_z_xy_y_zz, t_z_xy_z_xx, t_z_xy_z_xy,\
                           t_z_xy_z_xz, t_z_xy_z_yy, t_z_xy_z_yz, t_z_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xy_z_zz[i] = t_0_xyz_z_zz[i] - rab_z[i] * t_0_xy_z_zz[i];

        t_z_xy_z_yz[i] = t_0_xyz_z_yz[i] - rab_z[i] * t_0_xy_z_yz[i];

        t_z_xy_z_yy[i] = t_0_xyz_z_yy[i] - rab_z[i] * t_0_xy_z_yy[i];

        t_z_xy_z_xz[i] = t_0_xyz_z_xz[i] - rab_z[i] * t_0_xy_z_xz[i];

        t_z_xy_z_xy[i] = t_0_xyz_z_xy[i] - rab_z[i] * t_0_xy_z_xy[i];

        t_z_xy_z_xx[i] = t_0_xyz_z_xx[i] - rab_z[i] * t_0_xy_z_xx[i];

        t_z_xy_y_zz[i] = t_0_xyz_y_zz[i] - rab_z[i] * t_0_xy_y_zz[i];

        t_z_xy_y_yz[i] = t_0_xyz_y_yz[i] - rab_z[i] * t_0_xy_y_yz[i];

        t_z_xy_y_yy[i] = t_0_xyz_y_yy[i] - rab_z[i] * t_0_xy_y_yy[i];

        t_z_xy_y_xz[i] = t_0_xyz_y_xz[i] - rab_z[i] * t_0_xy_y_xz[i];

        t_z_xy_y_xy[i] = t_0_xyz_y_xy[i] - rab_z[i] * t_0_xy_y_xy[i];

        t_z_xy_y_xx[i] = t_0_xyz_y_xx[i] - rab_z[i] * t_0_xy_y_xx[i];

        t_z_xy_x_zz[i] = t_0_xyz_x_zz[i] - rab_z[i] * t_0_xy_x_zz[i];

        t_z_xy_x_yz[i] = t_0_xyz_x_yz[i] - rab_z[i] * t_0_xy_x_yz[i];

        t_z_xy_x_yy[i] = t_0_xyz_x_yy[i] - rab_z[i] * t_0_xy_x_yy[i];

        t_z_xy_x_xz[i] = t_0_xyz_x_xz[i] - rab_z[i] * t_0_xy_x_xz[i];

        t_z_xy_x_xy[i] = t_0_xyz_x_xy[i] - rab_z[i] * t_0_xy_x_xy[i];

        t_z_xy_x_xx[i] = t_0_xyz_x_xx[i] - rab_z[i] * t_0_xy_x_xx[i];

        t_z_xx_z_zz[i] = t_0_xxz_z_zz[i] - rab_z[i] * t_0_xx_z_zz[i];

        t_z_xx_z_yz[i] = t_0_xxz_z_yz[i] - rab_z[i] * t_0_xx_z_yz[i];

        t_z_xx_z_yy[i] = t_0_xxz_z_yy[i] - rab_z[i] * t_0_xx_z_yy[i];

        t_z_xx_z_xz[i] = t_0_xxz_z_xz[i] - rab_z[i] * t_0_xx_z_xz[i];

        t_z_xx_z_xy[i] = t_0_xxz_z_xy[i] - rab_z[i] * t_0_xx_z_xy[i];

        t_z_xx_z_xx[i] = t_0_xxz_z_xx[i] - rab_z[i] * t_0_xx_z_xx[i];

        t_z_xx_y_zz[i] = t_0_xxz_y_zz[i] - rab_z[i] * t_0_xx_y_zz[i];

        t_z_xx_y_yz[i] = t_0_xxz_y_yz[i] - rab_z[i] * t_0_xx_y_yz[i];

        t_z_xx_y_yy[i] = t_0_xxz_y_yy[i] - rab_z[i] * t_0_xx_y_yy[i];

        t_z_xx_y_xz[i] = t_0_xxz_y_xz[i] - rab_z[i] * t_0_xx_y_xz[i];

        t_z_xx_y_xy[i] = t_0_xxz_y_xy[i] - rab_z[i] * t_0_xx_y_xy[i];

        t_z_xx_y_xx[i] = t_0_xxz_y_xx[i] - rab_z[i] * t_0_xx_y_xx[i];

        t_z_xx_x_zz[i] = t_0_xxz_x_zz[i] - rab_z[i] * t_0_xx_x_zz[i];

        t_z_xx_x_yz[i] = t_0_xxz_x_yz[i] - rab_z[i] * t_0_xx_x_yz[i];

        t_z_xx_x_yy[i] = t_0_xxz_x_yy[i] - rab_z[i] * t_0_xx_x_yy[i];

        t_z_xx_x_xz[i] = t_0_xxz_x_xz[i] - rab_z[i] * t_0_xx_x_xz[i];

        t_z_xx_x_xy[i] = t_0_xxz_x_xy[i] - rab_z[i] * t_0_xx_x_xy[i];

        t_z_xx_x_xx[i] = t_0_xxz_x_xx[i] - rab_z[i] * t_0_xx_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yyz_x_xx, t_0_yyz_x_xy, t_0_yyz_x_xz, t_0_yyz_x_yy,\
                           t_0_yyz_x_yz, t_0_yyz_x_zz, t_0_yyz_y_xx, t_0_yyz_y_xy, t_0_yyz_y_xz,\
                           t_0_yyz_y_yy, t_0_yyz_y_yz, t_0_yyz_y_zz, t_0_yyz_z_xx, t_0_yyz_z_xy,\
                           t_0_yyz_z_xz, t_0_yyz_z_yy, t_0_yyz_z_yz, t_0_yyz_z_zz, t_0_yz_x_xx,\
                           t_0_yz_x_xy, t_0_yz_x_xz, t_0_yz_x_yy, t_0_yz_x_yz, t_0_yz_x_zz,\
                           t_0_yz_y_xx, t_0_yz_y_xy, t_0_yz_y_xz, t_0_yz_y_yy, t_0_yz_y_yz,\
                           t_0_yz_y_zz, t_0_yz_z_xx, t_0_yz_z_xy, t_0_yz_z_xz, t_0_yz_z_yy,\
                           t_0_yz_z_yz, t_0_yz_z_zz, t_0_yzz_x_xx, t_0_yzz_x_xy, t_0_yzz_x_xz,\
                           t_0_yzz_x_yy, t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_xx, t_0_yzz_y_xy,\
                           t_0_yzz_y_xz, t_0_yzz_y_yy, t_0_yzz_y_yz, t_0_yzz_y_zz, t_0_yzz_z_xx,\
                           t_0_yzz_z_xy, t_0_yzz_z_xz, t_0_yzz_z_yy, t_0_yzz_z_yz, t_0_yzz_z_zz,\
                           t_0_zz_x_xx, t_0_zz_x_xy, t_0_zz_x_xz, t_0_zz_x_yy, t_0_zz_x_yz,\
                           t_0_zz_x_zz, t_0_zz_y_xx, t_0_zz_y_xy, t_0_zz_y_xz, t_0_zz_y_yy,\
                           t_0_zz_y_yz, t_0_zz_y_zz, t_0_zz_z_xx, t_0_zz_z_xy, t_0_zz_z_xz,\
                           t_0_zz_z_yy, t_0_zz_z_yz, t_0_zz_z_zz, t_y_yz_x_xx, t_y_yz_x_xy,\
                           t_y_yz_x_xz, t_y_yz_x_yy, t_y_yz_x_yz, t_y_yz_x_zz, t_y_yz_y_xx,\
                           t_y_yz_y_xy, t_y_yz_y_xz, t_y_yz_y_yy, t_y_yz_y_yz, t_y_yz_y_zz,\
                           t_y_yz_z_xx, t_y_yz_z_xy, t_y_yz_z_xz, t_y_yz_z_yy, t_y_yz_z_yz,\
                           t_y_yz_z_zz, t_y_zz_x_xx, t_y_zz_x_xy, t_y_zz_x_xz, t_y_zz_x_yy,\
                           t_y_zz_x_yz, t_y_zz_x_zz, t_y_zz_y_xx, t_y_zz_y_xy, t_y_zz_y_xz,\
                           t_y_zz_y_yy, t_y_zz_y_yz, t_y_zz_y_zz, t_y_zz_z_xx, t_y_zz_z_xy,\
                           t_y_zz_z_xz, t_y_zz_z_yy, t_y_zz_z_yz, t_y_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zz_z_zz[i] = t_0_yzz_z_zz[i] - rab_y[i] * t_0_zz_z_zz[i];

        t_y_zz_z_yz[i] = t_0_yzz_z_yz[i] - rab_y[i] * t_0_zz_z_yz[i];

        t_y_zz_z_yy[i] = t_0_yzz_z_yy[i] - rab_y[i] * t_0_zz_z_yy[i];

        t_y_zz_z_xz[i] = t_0_yzz_z_xz[i] - rab_y[i] * t_0_zz_z_xz[i];

        t_y_zz_z_xy[i] = t_0_yzz_z_xy[i] - rab_y[i] * t_0_zz_z_xy[i];

        t_y_zz_z_xx[i] = t_0_yzz_z_xx[i] - rab_y[i] * t_0_zz_z_xx[i];

        t_y_zz_y_zz[i] = t_0_yzz_y_zz[i] - rab_y[i] * t_0_zz_y_zz[i];

        t_y_zz_y_yz[i] = t_0_yzz_y_yz[i] - rab_y[i] * t_0_zz_y_yz[i];

        t_y_zz_y_yy[i] = t_0_yzz_y_yy[i] - rab_y[i] * t_0_zz_y_yy[i];

        t_y_zz_y_xz[i] = t_0_yzz_y_xz[i] - rab_y[i] * t_0_zz_y_xz[i];

        t_y_zz_y_xy[i] = t_0_yzz_y_xy[i] - rab_y[i] * t_0_zz_y_xy[i];

        t_y_zz_y_xx[i] = t_0_yzz_y_xx[i] - rab_y[i] * t_0_zz_y_xx[i];

        t_y_zz_x_zz[i] = t_0_yzz_x_zz[i] - rab_y[i] * t_0_zz_x_zz[i];

        t_y_zz_x_yz[i] = t_0_yzz_x_yz[i] - rab_y[i] * t_0_zz_x_yz[i];

        t_y_zz_x_yy[i] = t_0_yzz_x_yy[i] - rab_y[i] * t_0_zz_x_yy[i];

        t_y_zz_x_xz[i] = t_0_yzz_x_xz[i] - rab_y[i] * t_0_zz_x_xz[i];

        t_y_zz_x_xy[i] = t_0_yzz_x_xy[i] - rab_y[i] * t_0_zz_x_xy[i];

        t_y_zz_x_xx[i] = t_0_yzz_x_xx[i] - rab_y[i] * t_0_zz_x_xx[i];

        t_y_yz_z_zz[i] = t_0_yyz_z_zz[i] - rab_y[i] * t_0_yz_z_zz[i];

        t_y_yz_z_yz[i] = t_0_yyz_z_yz[i] - rab_y[i] * t_0_yz_z_yz[i];

        t_y_yz_z_yy[i] = t_0_yyz_z_yy[i] - rab_y[i] * t_0_yz_z_yy[i];

        t_y_yz_z_xz[i] = t_0_yyz_z_xz[i] - rab_y[i] * t_0_yz_z_xz[i];

        t_y_yz_z_xy[i] = t_0_yyz_z_xy[i] - rab_y[i] * t_0_yz_z_xy[i];

        t_y_yz_z_xx[i] = t_0_yyz_z_xx[i] - rab_y[i] * t_0_yz_z_xx[i];

        t_y_yz_y_zz[i] = t_0_yyz_y_zz[i] - rab_y[i] * t_0_yz_y_zz[i];

        t_y_yz_y_yz[i] = t_0_yyz_y_yz[i] - rab_y[i] * t_0_yz_y_yz[i];

        t_y_yz_y_yy[i] = t_0_yyz_y_yy[i] - rab_y[i] * t_0_yz_y_yy[i];

        t_y_yz_y_xz[i] = t_0_yyz_y_xz[i] - rab_y[i] * t_0_yz_y_xz[i];

        t_y_yz_y_xy[i] = t_0_yyz_y_xy[i] - rab_y[i] * t_0_yz_y_xy[i];

        t_y_yz_y_xx[i] = t_0_yyz_y_xx[i] - rab_y[i] * t_0_yz_y_xx[i];

        t_y_yz_x_zz[i] = t_0_yyz_x_zz[i] - rab_y[i] * t_0_yz_x_zz[i];

        t_y_yz_x_yz[i] = t_0_yyz_x_yz[i] - rab_y[i] * t_0_yz_x_yz[i];

        t_y_yz_x_yy[i] = t_0_yyz_x_yy[i] - rab_y[i] * t_0_yz_x_yy[i];

        t_y_yz_x_xz[i] = t_0_yyz_x_xz[i] - rab_y[i] * t_0_yz_x_xz[i];

        t_y_yz_x_xy[i] = t_0_yyz_x_xy[i] - rab_y[i] * t_0_yz_x_xy[i];

        t_y_yz_x_xx[i] = t_0_yyz_x_xx[i] - rab_y[i] * t_0_yz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xyz_x_xx, t_0_xyz_x_xy, t_0_xyz_x_xz, t_0_xyz_x_yy,\
                           t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx, t_0_xyz_y_xy, t_0_xyz_y_xz,\
                           t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz, t_0_xyz_z_xx, t_0_xyz_z_xy,\
                           t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz, t_0_xyz_z_zz, t_0_xz_x_xx,\
                           t_0_xz_x_xy, t_0_xz_x_xz, t_0_xz_x_yy, t_0_xz_x_yz, t_0_xz_x_zz,\
                           t_0_xz_y_xx, t_0_xz_y_xy, t_0_xz_y_xz, t_0_xz_y_yy, t_0_xz_y_yz,\
                           t_0_xz_y_zz, t_0_xz_z_xx, t_0_xz_z_xy, t_0_xz_z_xz, t_0_xz_z_yy,\
                           t_0_xz_z_yz, t_0_xz_z_zz, t_0_yy_x_xx, t_0_yy_x_xy, t_0_yy_x_xz,\
                           t_0_yy_x_yy, t_0_yy_x_yz, t_0_yy_x_zz, t_0_yy_y_xx, t_0_yy_y_xy,\
                           t_0_yy_y_xz, t_0_yy_y_yy, t_0_yy_y_yz, t_0_yy_y_zz, t_0_yy_z_xx,\
                           t_0_yy_z_xy, t_0_yy_z_xz, t_0_yy_z_yy, t_0_yy_z_yz, t_0_yy_z_zz,\
                           t_0_yyy_x_xx, t_0_yyy_x_xy, t_0_yyy_x_xz, t_0_yyy_x_yy, t_0_yyy_x_yz,\
                           t_0_yyy_x_zz, t_0_yyy_y_xx, t_0_yyy_y_xy, t_0_yyy_y_xz, t_0_yyy_y_yy,\
                           t_0_yyy_y_yz, t_0_yyy_y_zz, t_0_yyy_z_xx, t_0_yyy_z_xy, t_0_yyy_z_xz,\
                           t_0_yyy_z_yy, t_0_yyy_z_yz, t_0_yyy_z_zz, t_y_xz_x_xx, t_y_xz_x_xy,\
                           t_y_xz_x_xz, t_y_xz_x_yy, t_y_xz_x_yz, t_y_xz_x_zz, t_y_xz_y_xx,\
                           t_y_xz_y_xy, t_y_xz_y_xz, t_y_xz_y_yy, t_y_xz_y_yz, t_y_xz_y_zz,\
                           t_y_xz_z_xx, t_y_xz_z_xy, t_y_xz_z_xz, t_y_xz_z_yy, t_y_xz_z_yz,\
                           t_y_xz_z_zz, t_y_yy_x_xx, t_y_yy_x_xy, t_y_yy_x_xz, t_y_yy_x_yy,\
                           t_y_yy_x_yz, t_y_yy_x_zz, t_y_yy_y_xx, t_y_yy_y_xy, t_y_yy_y_xz,\
                           t_y_yy_y_yy, t_y_yy_y_yz, t_y_yy_y_zz, t_y_yy_z_xx, t_y_yy_z_xy,\
                           t_y_yy_z_xz, t_y_yy_z_yy, t_y_yy_z_yz, t_y_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yy_z_zz[i] = t_0_yyy_z_zz[i] - rab_y[i] * t_0_yy_z_zz[i];

        t_y_yy_z_yz[i] = t_0_yyy_z_yz[i] - rab_y[i] * t_0_yy_z_yz[i];

        t_y_yy_z_yy[i] = t_0_yyy_z_yy[i] - rab_y[i] * t_0_yy_z_yy[i];

        t_y_yy_z_xz[i] = t_0_yyy_z_xz[i] - rab_y[i] * t_0_yy_z_xz[i];

        t_y_yy_z_xy[i] = t_0_yyy_z_xy[i] - rab_y[i] * t_0_yy_z_xy[i];

        t_y_yy_z_xx[i] = t_0_yyy_z_xx[i] - rab_y[i] * t_0_yy_z_xx[i];

        t_y_yy_y_zz[i] = t_0_yyy_y_zz[i] - rab_y[i] * t_0_yy_y_zz[i];

        t_y_yy_y_yz[i] = t_0_yyy_y_yz[i] - rab_y[i] * t_0_yy_y_yz[i];

        t_y_yy_y_yy[i] = t_0_yyy_y_yy[i] - rab_y[i] * t_0_yy_y_yy[i];

        t_y_yy_y_xz[i] = t_0_yyy_y_xz[i] - rab_y[i] * t_0_yy_y_xz[i];

        t_y_yy_y_xy[i] = t_0_yyy_y_xy[i] - rab_y[i] * t_0_yy_y_xy[i];

        t_y_yy_y_xx[i] = t_0_yyy_y_xx[i] - rab_y[i] * t_0_yy_y_xx[i];

        t_y_yy_x_zz[i] = t_0_yyy_x_zz[i] - rab_y[i] * t_0_yy_x_zz[i];

        t_y_yy_x_yz[i] = t_0_yyy_x_yz[i] - rab_y[i] * t_0_yy_x_yz[i];

        t_y_yy_x_yy[i] = t_0_yyy_x_yy[i] - rab_y[i] * t_0_yy_x_yy[i];

        t_y_yy_x_xz[i] = t_0_yyy_x_xz[i] - rab_y[i] * t_0_yy_x_xz[i];

        t_y_yy_x_xy[i] = t_0_yyy_x_xy[i] - rab_y[i] * t_0_yy_x_xy[i];

        t_y_yy_x_xx[i] = t_0_yyy_x_xx[i] - rab_y[i] * t_0_yy_x_xx[i];

        t_y_xz_z_zz[i] = t_0_xyz_z_zz[i] - rab_y[i] * t_0_xz_z_zz[i];

        t_y_xz_z_yz[i] = t_0_xyz_z_yz[i] - rab_y[i] * t_0_xz_z_yz[i];

        t_y_xz_z_yy[i] = t_0_xyz_z_yy[i] - rab_y[i] * t_0_xz_z_yy[i];

        t_y_xz_z_xz[i] = t_0_xyz_z_xz[i] - rab_y[i] * t_0_xz_z_xz[i];

        t_y_xz_z_xy[i] = t_0_xyz_z_xy[i] - rab_y[i] * t_0_xz_z_xy[i];

        t_y_xz_z_xx[i] = t_0_xyz_z_xx[i] - rab_y[i] * t_0_xz_z_xx[i];

        t_y_xz_y_zz[i] = t_0_xyz_y_zz[i] - rab_y[i] * t_0_xz_y_zz[i];

        t_y_xz_y_yz[i] = t_0_xyz_y_yz[i] - rab_y[i] * t_0_xz_y_yz[i];

        t_y_xz_y_yy[i] = t_0_xyz_y_yy[i] - rab_y[i] * t_0_xz_y_yy[i];

        t_y_xz_y_xz[i] = t_0_xyz_y_xz[i] - rab_y[i] * t_0_xz_y_xz[i];

        t_y_xz_y_xy[i] = t_0_xyz_y_xy[i] - rab_y[i] * t_0_xz_y_xy[i];

        t_y_xz_y_xx[i] = t_0_xyz_y_xx[i] - rab_y[i] * t_0_xz_y_xx[i];

        t_y_xz_x_zz[i] = t_0_xyz_x_zz[i] - rab_y[i] * t_0_xz_x_zz[i];

        t_y_xz_x_yz[i] = t_0_xyz_x_yz[i] - rab_y[i] * t_0_xz_x_yz[i];

        t_y_xz_x_yy[i] = t_0_xyz_x_yy[i] - rab_y[i] * t_0_xz_x_yy[i];

        t_y_xz_x_xz[i] = t_0_xyz_x_xz[i] - rab_y[i] * t_0_xz_x_xz[i];

        t_y_xz_x_xy[i] = t_0_xyz_x_xy[i] - rab_y[i] * t_0_xz_x_xy[i];

        t_y_xz_x_xx[i] = t_0_xyz_x_xx[i] - rab_y[i] * t_0_xz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xx_x_xx, t_0_xx_x_xy, t_0_xx_x_xz, t_0_xx_x_yy,\
                           t_0_xx_x_yz, t_0_xx_x_zz, t_0_xx_y_xx, t_0_xx_y_xy, t_0_xx_y_xz,\
                           t_0_xx_y_yy, t_0_xx_y_yz, t_0_xx_y_zz, t_0_xx_z_xx, t_0_xx_z_xy,\
                           t_0_xx_z_xz, t_0_xx_z_yy, t_0_xx_z_yz, t_0_xx_z_zz, t_0_xxy_x_xx,\
                           t_0_xxy_x_xy, t_0_xxy_x_xz, t_0_xxy_x_yy, t_0_xxy_x_yz, t_0_xxy_x_zz,\
                           t_0_xxy_y_xx, t_0_xxy_y_xy, t_0_xxy_y_xz, t_0_xxy_y_yy, t_0_xxy_y_yz,\
                           t_0_xxy_y_zz, t_0_xxy_z_xx, t_0_xxy_z_xy, t_0_xxy_z_xz, t_0_xxy_z_yy,\
                           t_0_xxy_z_yz, t_0_xxy_z_zz, t_0_xy_x_xx, t_0_xy_x_xy, t_0_xy_x_xz,\
                           t_0_xy_x_yy, t_0_xy_x_yz, t_0_xy_x_zz, t_0_xy_y_xx, t_0_xy_y_xy,\
                           t_0_xy_y_xz, t_0_xy_y_yy, t_0_xy_y_yz, t_0_xy_y_zz, t_0_xy_z_xx,\
                           t_0_xy_z_xy, t_0_xy_z_xz, t_0_xy_z_yy, t_0_xy_z_yz, t_0_xy_z_zz,\
                           t_0_xyy_x_xx, t_0_xyy_x_xy, t_0_xyy_x_xz, t_0_xyy_x_yy, t_0_xyy_x_yz,\
                           t_0_xyy_x_zz, t_0_xyy_y_xx, t_0_xyy_y_xy, t_0_xyy_y_xz, t_0_xyy_y_yy,\
                           t_0_xyy_y_yz, t_0_xyy_y_zz, t_0_xyy_z_xx, t_0_xyy_z_xy, t_0_xyy_z_xz,\
                           t_0_xyy_z_yy, t_0_xyy_z_yz, t_0_xyy_z_zz, t_y_xx_x_xx, t_y_xx_x_xy,\
                           t_y_xx_x_xz, t_y_xx_x_yy, t_y_xx_x_yz, t_y_xx_x_zz, t_y_xx_y_xx,\
                           t_y_xx_y_xy, t_y_xx_y_xz, t_y_xx_y_yy, t_y_xx_y_yz, t_y_xx_y_zz,\
                           t_y_xx_z_xx, t_y_xx_z_xy, t_y_xx_z_xz, t_y_xx_z_yy, t_y_xx_z_yz,\
                           t_y_xx_z_zz, t_y_xy_x_xx, t_y_xy_x_xy, t_y_xy_x_xz, t_y_xy_x_yy,\
                           t_y_xy_x_yz, t_y_xy_x_zz, t_y_xy_y_xx, t_y_xy_y_xy, t_y_xy_y_xz,\
                           t_y_xy_y_yy, t_y_xy_y_yz, t_y_xy_y_zz, t_y_xy_z_xx, t_y_xy_z_xy,\
                           t_y_xy_z_xz, t_y_xy_z_yy, t_y_xy_z_yz, t_y_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xy_z_zz[i] = t_0_xyy_z_zz[i] - rab_y[i] * t_0_xy_z_zz[i];

        t_y_xy_z_yz[i] = t_0_xyy_z_yz[i] - rab_y[i] * t_0_xy_z_yz[i];

        t_y_xy_z_yy[i] = t_0_xyy_z_yy[i] - rab_y[i] * t_0_xy_z_yy[i];

        t_y_xy_z_xz[i] = t_0_xyy_z_xz[i] - rab_y[i] * t_0_xy_z_xz[i];

        t_y_xy_z_xy[i] = t_0_xyy_z_xy[i] - rab_y[i] * t_0_xy_z_xy[i];

        t_y_xy_z_xx[i] = t_0_xyy_z_xx[i] - rab_y[i] * t_0_xy_z_xx[i];

        t_y_xy_y_zz[i] = t_0_xyy_y_zz[i] - rab_y[i] * t_0_xy_y_zz[i];

        t_y_xy_y_yz[i] = t_0_xyy_y_yz[i] - rab_y[i] * t_0_xy_y_yz[i];

        t_y_xy_y_yy[i] = t_0_xyy_y_yy[i] - rab_y[i] * t_0_xy_y_yy[i];

        t_y_xy_y_xz[i] = t_0_xyy_y_xz[i] - rab_y[i] * t_0_xy_y_xz[i];

        t_y_xy_y_xy[i] = t_0_xyy_y_xy[i] - rab_y[i] * t_0_xy_y_xy[i];

        t_y_xy_y_xx[i] = t_0_xyy_y_xx[i] - rab_y[i] * t_0_xy_y_xx[i];

        t_y_xy_x_zz[i] = t_0_xyy_x_zz[i] - rab_y[i] * t_0_xy_x_zz[i];

        t_y_xy_x_yz[i] = t_0_xyy_x_yz[i] - rab_y[i] * t_0_xy_x_yz[i];

        t_y_xy_x_yy[i] = t_0_xyy_x_yy[i] - rab_y[i] * t_0_xy_x_yy[i];

        t_y_xy_x_xz[i] = t_0_xyy_x_xz[i] - rab_y[i] * t_0_xy_x_xz[i];

        t_y_xy_x_xy[i] = t_0_xyy_x_xy[i] - rab_y[i] * t_0_xy_x_xy[i];

        t_y_xy_x_xx[i] = t_0_xyy_x_xx[i] - rab_y[i] * t_0_xy_x_xx[i];

        t_y_xx_z_zz[i] = t_0_xxy_z_zz[i] - rab_y[i] * t_0_xx_z_zz[i];

        t_y_xx_z_yz[i] = t_0_xxy_z_yz[i] - rab_y[i] * t_0_xx_z_yz[i];

        t_y_xx_z_yy[i] = t_0_xxy_z_yy[i] - rab_y[i] * t_0_xx_z_yy[i];

        t_y_xx_z_xz[i] = t_0_xxy_z_xz[i] - rab_y[i] * t_0_xx_z_xz[i];

        t_y_xx_z_xy[i] = t_0_xxy_z_xy[i] - rab_y[i] * t_0_xx_z_xy[i];

        t_y_xx_z_xx[i] = t_0_xxy_z_xx[i] - rab_y[i] * t_0_xx_z_xx[i];

        t_y_xx_y_zz[i] = t_0_xxy_y_zz[i] - rab_y[i] * t_0_xx_y_zz[i];

        t_y_xx_y_yz[i] = t_0_xxy_y_yz[i] - rab_y[i] * t_0_xx_y_yz[i];

        t_y_xx_y_yy[i] = t_0_xxy_y_yy[i] - rab_y[i] * t_0_xx_y_yy[i];

        t_y_xx_y_xz[i] = t_0_xxy_y_xz[i] - rab_y[i] * t_0_xx_y_xz[i];

        t_y_xx_y_xy[i] = t_0_xxy_y_xy[i] - rab_y[i] * t_0_xx_y_xy[i];

        t_y_xx_y_xx[i] = t_0_xxy_y_xx[i] - rab_y[i] * t_0_xx_y_xx[i];

        t_y_xx_x_zz[i] = t_0_xxy_x_zz[i] - rab_y[i] * t_0_xx_x_zz[i];

        t_y_xx_x_yz[i] = t_0_xxy_x_yz[i] - rab_y[i] * t_0_xx_x_yz[i];

        t_y_xx_x_yy[i] = t_0_xxy_x_yy[i] - rab_y[i] * t_0_xx_x_yy[i];

        t_y_xx_x_xz[i] = t_0_xxy_x_xz[i] - rab_y[i] * t_0_xx_x_xz[i];

        t_y_xx_x_xy[i] = t_0_xxy_x_xy[i] - rab_y[i] * t_0_xx_x_xy[i];

        t_y_xx_x_xx[i] = t_0_xxy_x_xx[i] - rab_y[i] * t_0_xx_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xyz_x_xx, t_0_xyz_x_xy, t_0_xyz_x_xz, t_0_xyz_x_yy,\
                           t_0_xyz_x_yz, t_0_xyz_x_zz, t_0_xyz_y_xx, t_0_xyz_y_xy, t_0_xyz_y_xz,\
                           t_0_xyz_y_yy, t_0_xyz_y_yz, t_0_xyz_y_zz, t_0_xyz_z_xx, t_0_xyz_z_xy,\
                           t_0_xyz_z_xz, t_0_xyz_z_yy, t_0_xyz_z_yz, t_0_xyz_z_zz, t_0_xzz_x_xx,\
                           t_0_xzz_x_xy, t_0_xzz_x_xz, t_0_xzz_x_yy, t_0_xzz_x_yz, t_0_xzz_x_zz,\
                           t_0_xzz_y_xx, t_0_xzz_y_xy, t_0_xzz_y_xz, t_0_xzz_y_yy, t_0_xzz_y_yz,\
                           t_0_xzz_y_zz, t_0_xzz_z_xx, t_0_xzz_z_xy, t_0_xzz_z_xz, t_0_xzz_z_yy,\
                           t_0_xzz_z_yz, t_0_xzz_z_zz, t_0_yz_x_xx, t_0_yz_x_xy, t_0_yz_x_xz,\
                           t_0_yz_x_yy, t_0_yz_x_yz, t_0_yz_x_zz, t_0_yz_y_xx, t_0_yz_y_xy,\
                           t_0_yz_y_xz, t_0_yz_y_yy, t_0_yz_y_yz, t_0_yz_y_zz, t_0_yz_z_xx,\
                           t_0_yz_z_xy, t_0_yz_z_xz, t_0_yz_z_yy, t_0_yz_z_yz, t_0_yz_z_zz,\
                           t_0_zz_x_xx, t_0_zz_x_xy, t_0_zz_x_xz, t_0_zz_x_yy, t_0_zz_x_yz,\
                           t_0_zz_x_zz, t_0_zz_y_xx, t_0_zz_y_xy, t_0_zz_y_xz, t_0_zz_y_yy,\
                           t_0_zz_y_yz, t_0_zz_y_zz, t_0_zz_z_xx, t_0_zz_z_xy, t_0_zz_z_xz,\
                           t_0_zz_z_yy, t_0_zz_z_yz, t_0_zz_z_zz, t_x_yz_x_xx, t_x_yz_x_xy,\
                           t_x_yz_x_xz, t_x_yz_x_yy, t_x_yz_x_yz, t_x_yz_x_zz, t_x_yz_y_xx,\
                           t_x_yz_y_xy, t_x_yz_y_xz, t_x_yz_y_yy, t_x_yz_y_yz, t_x_yz_y_zz,\
                           t_x_yz_z_xx, t_x_yz_z_xy, t_x_yz_z_xz, t_x_yz_z_yy, t_x_yz_z_yz,\
                           t_x_yz_z_zz, t_x_zz_x_xx, t_x_zz_x_xy, t_x_zz_x_xz, t_x_zz_x_yy,\
                           t_x_zz_x_yz, t_x_zz_x_zz, t_x_zz_y_xx, t_x_zz_y_xy, t_x_zz_y_xz,\
                           t_x_zz_y_yy, t_x_zz_y_yz, t_x_zz_y_zz, t_x_zz_z_xx, t_x_zz_z_xy,\
                           t_x_zz_z_xz, t_x_zz_z_yy, t_x_zz_z_yz, t_x_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zz_z_zz[i] = t_0_xzz_z_zz[i] - rab_x[i] * t_0_zz_z_zz[i];

        t_x_zz_z_yz[i] = t_0_xzz_z_yz[i] - rab_x[i] * t_0_zz_z_yz[i];

        t_x_zz_z_yy[i] = t_0_xzz_z_yy[i] - rab_x[i] * t_0_zz_z_yy[i];

        t_x_zz_z_xz[i] = t_0_xzz_z_xz[i] - rab_x[i] * t_0_zz_z_xz[i];

        t_x_zz_z_xy[i] = t_0_xzz_z_xy[i] - rab_x[i] * t_0_zz_z_xy[i];

        t_x_zz_z_xx[i] = t_0_xzz_z_xx[i] - rab_x[i] * t_0_zz_z_xx[i];

        t_x_zz_y_zz[i] = t_0_xzz_y_zz[i] - rab_x[i] * t_0_zz_y_zz[i];

        t_x_zz_y_yz[i] = t_0_xzz_y_yz[i] - rab_x[i] * t_0_zz_y_yz[i];

        t_x_zz_y_yy[i] = t_0_xzz_y_yy[i] - rab_x[i] * t_0_zz_y_yy[i];

        t_x_zz_y_xz[i] = t_0_xzz_y_xz[i] - rab_x[i] * t_0_zz_y_xz[i];

        t_x_zz_y_xy[i] = t_0_xzz_y_xy[i] - rab_x[i] * t_0_zz_y_xy[i];

        t_x_zz_y_xx[i] = t_0_xzz_y_xx[i] - rab_x[i] * t_0_zz_y_xx[i];

        t_x_zz_x_zz[i] = t_0_xzz_x_zz[i] - rab_x[i] * t_0_zz_x_zz[i];

        t_x_zz_x_yz[i] = t_0_xzz_x_yz[i] - rab_x[i] * t_0_zz_x_yz[i];

        t_x_zz_x_yy[i] = t_0_xzz_x_yy[i] - rab_x[i] * t_0_zz_x_yy[i];

        t_x_zz_x_xz[i] = t_0_xzz_x_xz[i] - rab_x[i] * t_0_zz_x_xz[i];

        t_x_zz_x_xy[i] = t_0_xzz_x_xy[i] - rab_x[i] * t_0_zz_x_xy[i];

        t_x_zz_x_xx[i] = t_0_xzz_x_xx[i] - rab_x[i] * t_0_zz_x_xx[i];

        t_x_yz_z_zz[i] = t_0_xyz_z_zz[i] - rab_x[i] * t_0_yz_z_zz[i];

        t_x_yz_z_yz[i] = t_0_xyz_z_yz[i] - rab_x[i] * t_0_yz_z_yz[i];

        t_x_yz_z_yy[i] = t_0_xyz_z_yy[i] - rab_x[i] * t_0_yz_z_yy[i];

        t_x_yz_z_xz[i] = t_0_xyz_z_xz[i] - rab_x[i] * t_0_yz_z_xz[i];

        t_x_yz_z_xy[i] = t_0_xyz_z_xy[i] - rab_x[i] * t_0_yz_z_xy[i];

        t_x_yz_z_xx[i] = t_0_xyz_z_xx[i] - rab_x[i] * t_0_yz_z_xx[i];

        t_x_yz_y_zz[i] = t_0_xyz_y_zz[i] - rab_x[i] * t_0_yz_y_zz[i];

        t_x_yz_y_yz[i] = t_0_xyz_y_yz[i] - rab_x[i] * t_0_yz_y_yz[i];

        t_x_yz_y_yy[i] = t_0_xyz_y_yy[i] - rab_x[i] * t_0_yz_y_yy[i];

        t_x_yz_y_xz[i] = t_0_xyz_y_xz[i] - rab_x[i] * t_0_yz_y_xz[i];

        t_x_yz_y_xy[i] = t_0_xyz_y_xy[i] - rab_x[i] * t_0_yz_y_xy[i];

        t_x_yz_y_xx[i] = t_0_xyz_y_xx[i] - rab_x[i] * t_0_yz_y_xx[i];

        t_x_yz_x_zz[i] = t_0_xyz_x_zz[i] - rab_x[i] * t_0_yz_x_zz[i];

        t_x_yz_x_yz[i] = t_0_xyz_x_yz[i] - rab_x[i] * t_0_yz_x_yz[i];

        t_x_yz_x_yy[i] = t_0_xyz_x_yy[i] - rab_x[i] * t_0_yz_x_yy[i];

        t_x_yz_x_xz[i] = t_0_xyz_x_xz[i] - rab_x[i] * t_0_yz_x_xz[i];

        t_x_yz_x_xy[i] = t_0_xyz_x_xy[i] - rab_x[i] * t_0_yz_x_xy[i];

        t_x_yz_x_xx[i] = t_0_xyz_x_xx[i] - rab_x[i] * t_0_yz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xxz_x_xx, t_0_xxz_x_xy, t_0_xxz_x_xz, t_0_xxz_x_yy,\
                           t_0_xxz_x_yz, t_0_xxz_x_zz, t_0_xxz_y_xx, t_0_xxz_y_xy, t_0_xxz_y_xz,\
                           t_0_xxz_y_yy, t_0_xxz_y_yz, t_0_xxz_y_zz, t_0_xxz_z_xx, t_0_xxz_z_xy,\
                           t_0_xxz_z_xz, t_0_xxz_z_yy, t_0_xxz_z_yz, t_0_xxz_z_zz, t_0_xyy_x_xx,\
                           t_0_xyy_x_xy, t_0_xyy_x_xz, t_0_xyy_x_yy, t_0_xyy_x_yz, t_0_xyy_x_zz,\
                           t_0_xyy_y_xx, t_0_xyy_y_xy, t_0_xyy_y_xz, t_0_xyy_y_yy, t_0_xyy_y_yz,\
                           t_0_xyy_y_zz, t_0_xyy_z_xx, t_0_xyy_z_xy, t_0_xyy_z_xz, t_0_xyy_z_yy,\
                           t_0_xyy_z_yz, t_0_xyy_z_zz, t_0_xz_x_xx, t_0_xz_x_xy, t_0_xz_x_xz,\
                           t_0_xz_x_yy, t_0_xz_x_yz, t_0_xz_x_zz, t_0_xz_y_xx, t_0_xz_y_xy,\
                           t_0_xz_y_xz, t_0_xz_y_yy, t_0_xz_y_yz, t_0_xz_y_zz, t_0_xz_z_xx,\
                           t_0_xz_z_xy, t_0_xz_z_xz, t_0_xz_z_yy, t_0_xz_z_yz, t_0_xz_z_zz,\
                           t_0_yy_x_xx, t_0_yy_x_xy, t_0_yy_x_xz, t_0_yy_x_yy, t_0_yy_x_yz,\
                           t_0_yy_x_zz, t_0_yy_y_xx, t_0_yy_y_xy, t_0_yy_y_xz, t_0_yy_y_yy,\
                           t_0_yy_y_yz, t_0_yy_y_zz, t_0_yy_z_xx, t_0_yy_z_xy, t_0_yy_z_xz,\
                           t_0_yy_z_yy, t_0_yy_z_yz, t_0_yy_z_zz, t_x_xz_x_xx, t_x_xz_x_xy,\
                           t_x_xz_x_xz, t_x_xz_x_yy, t_x_xz_x_yz, t_x_xz_x_zz, t_x_xz_y_xx,\
                           t_x_xz_y_xy, t_x_xz_y_xz, t_x_xz_y_yy, t_x_xz_y_yz, t_x_xz_y_zz,\
                           t_x_xz_z_xx, t_x_xz_z_xy, t_x_xz_z_xz, t_x_xz_z_yy, t_x_xz_z_yz,\
                           t_x_xz_z_zz, t_x_yy_x_xx, t_x_yy_x_xy, t_x_yy_x_xz, t_x_yy_x_yy,\
                           t_x_yy_x_yz, t_x_yy_x_zz, t_x_yy_y_xx, t_x_yy_y_xy, t_x_yy_y_xz,\
                           t_x_yy_y_yy, t_x_yy_y_yz, t_x_yy_y_zz, t_x_yy_z_xx, t_x_yy_z_xy,\
                           t_x_yy_z_xz, t_x_yy_z_yy, t_x_yy_z_yz, t_x_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yy_z_zz[i] = t_0_xyy_z_zz[i] - rab_x[i] * t_0_yy_z_zz[i];

        t_x_yy_z_yz[i] = t_0_xyy_z_yz[i] - rab_x[i] * t_0_yy_z_yz[i];

        t_x_yy_z_yy[i] = t_0_xyy_z_yy[i] - rab_x[i] * t_0_yy_z_yy[i];

        t_x_yy_z_xz[i] = t_0_xyy_z_xz[i] - rab_x[i] * t_0_yy_z_xz[i];

        t_x_yy_z_xy[i] = t_0_xyy_z_xy[i] - rab_x[i] * t_0_yy_z_xy[i];

        t_x_yy_z_xx[i] = t_0_xyy_z_xx[i] - rab_x[i] * t_0_yy_z_xx[i];

        t_x_yy_y_zz[i] = t_0_xyy_y_zz[i] - rab_x[i] * t_0_yy_y_zz[i];

        t_x_yy_y_yz[i] = t_0_xyy_y_yz[i] - rab_x[i] * t_0_yy_y_yz[i];

        t_x_yy_y_yy[i] = t_0_xyy_y_yy[i] - rab_x[i] * t_0_yy_y_yy[i];

        t_x_yy_y_xz[i] = t_0_xyy_y_xz[i] - rab_x[i] * t_0_yy_y_xz[i];

        t_x_yy_y_xy[i] = t_0_xyy_y_xy[i] - rab_x[i] * t_0_yy_y_xy[i];

        t_x_yy_y_xx[i] = t_0_xyy_y_xx[i] - rab_x[i] * t_0_yy_y_xx[i];

        t_x_yy_x_zz[i] = t_0_xyy_x_zz[i] - rab_x[i] * t_0_yy_x_zz[i];

        t_x_yy_x_yz[i] = t_0_xyy_x_yz[i] - rab_x[i] * t_0_yy_x_yz[i];

        t_x_yy_x_yy[i] = t_0_xyy_x_yy[i] - rab_x[i] * t_0_yy_x_yy[i];

        t_x_yy_x_xz[i] = t_0_xyy_x_xz[i] - rab_x[i] * t_0_yy_x_xz[i];

        t_x_yy_x_xy[i] = t_0_xyy_x_xy[i] - rab_x[i] * t_0_yy_x_xy[i];

        t_x_yy_x_xx[i] = t_0_xyy_x_xx[i] - rab_x[i] * t_0_yy_x_xx[i];

        t_x_xz_z_zz[i] = t_0_xxz_z_zz[i] - rab_x[i] * t_0_xz_z_zz[i];

        t_x_xz_z_yz[i] = t_0_xxz_z_yz[i] - rab_x[i] * t_0_xz_z_yz[i];

        t_x_xz_z_yy[i] = t_0_xxz_z_yy[i] - rab_x[i] * t_0_xz_z_yy[i];

        t_x_xz_z_xz[i] = t_0_xxz_z_xz[i] - rab_x[i] * t_0_xz_z_xz[i];

        t_x_xz_z_xy[i] = t_0_xxz_z_xy[i] - rab_x[i] * t_0_xz_z_xy[i];

        t_x_xz_z_xx[i] = t_0_xxz_z_xx[i] - rab_x[i] * t_0_xz_z_xx[i];

        t_x_xz_y_zz[i] = t_0_xxz_y_zz[i] - rab_x[i] * t_0_xz_y_zz[i];

        t_x_xz_y_yz[i] = t_0_xxz_y_yz[i] - rab_x[i] * t_0_xz_y_yz[i];

        t_x_xz_y_yy[i] = t_0_xxz_y_yy[i] - rab_x[i] * t_0_xz_y_yy[i];

        t_x_xz_y_xz[i] = t_0_xxz_y_xz[i] - rab_x[i] * t_0_xz_y_xz[i];

        t_x_xz_y_xy[i] = t_0_xxz_y_xy[i] - rab_x[i] * t_0_xz_y_xy[i];

        t_x_xz_y_xx[i] = t_0_xxz_y_xx[i] - rab_x[i] * t_0_xz_y_xx[i];

        t_x_xz_x_zz[i] = t_0_xxz_x_zz[i] - rab_x[i] * t_0_xz_x_zz[i];

        t_x_xz_x_yz[i] = t_0_xxz_x_yz[i] - rab_x[i] * t_0_xz_x_yz[i];

        t_x_xz_x_yy[i] = t_0_xxz_x_yy[i] - rab_x[i] * t_0_xz_x_yy[i];

        t_x_xz_x_xz[i] = t_0_xxz_x_xz[i] - rab_x[i] * t_0_xz_x_xz[i];

        t_x_xz_x_xy[i] = t_0_xxz_x_xy[i] - rab_x[i] * t_0_xz_x_xy[i];

        t_x_xz_x_xx[i] = t_0_xxz_x_xx[i] - rab_x[i] * t_0_xz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xx_x_xx, t_0_xx_x_xy, t_0_xx_x_xz, t_0_xx_x_yy,\
                           t_0_xx_x_yz, t_0_xx_x_zz, t_0_xx_y_xx, t_0_xx_y_xy, t_0_xx_y_xz,\
                           t_0_xx_y_yy, t_0_xx_y_yz, t_0_xx_y_zz, t_0_xx_z_xx, t_0_xx_z_xy,\
                           t_0_xx_z_xz, t_0_xx_z_yy, t_0_xx_z_yz, t_0_xx_z_zz, t_0_xxx_x_xx,\
                           t_0_xxx_x_xy, t_0_xxx_x_xz, t_0_xxx_x_yy, t_0_xxx_x_yz, t_0_xxx_x_zz,\
                           t_0_xxx_y_xx, t_0_xxx_y_xy, t_0_xxx_y_xz, t_0_xxx_y_yy, t_0_xxx_y_yz,\
                           t_0_xxx_y_zz, t_0_xxx_z_xx, t_0_xxx_z_xy, t_0_xxx_z_xz, t_0_xxx_z_yy,\
                           t_0_xxx_z_yz, t_0_xxx_z_zz, t_0_xxy_x_xx, t_0_xxy_x_xy, t_0_xxy_x_xz,\
                           t_0_xxy_x_yy, t_0_xxy_x_yz, t_0_xxy_x_zz, t_0_xxy_y_xx, t_0_xxy_y_xy,\
                           t_0_xxy_y_xz, t_0_xxy_y_yy, t_0_xxy_y_yz, t_0_xxy_y_zz, t_0_xxy_z_xx,\
                           t_0_xxy_z_xy, t_0_xxy_z_xz, t_0_xxy_z_yy, t_0_xxy_z_yz, t_0_xxy_z_zz,\
                           t_0_xy_x_xx, t_0_xy_x_xy, t_0_xy_x_xz, t_0_xy_x_yy, t_0_xy_x_yz,\
                           t_0_xy_x_zz, t_0_xy_y_xx, t_0_xy_y_xy, t_0_xy_y_xz, t_0_xy_y_yy,\
                           t_0_xy_y_yz, t_0_xy_y_zz, t_0_xy_z_xx, t_0_xy_z_xy, t_0_xy_z_xz,\
                           t_0_xy_z_yy, t_0_xy_z_yz, t_0_xy_z_zz, t_x_xx_x_xx, t_x_xx_x_xy,\
                           t_x_xx_x_xz, t_x_xx_x_yy, t_x_xx_x_yz, t_x_xx_x_zz, t_x_xx_y_xx,\
                           t_x_xx_y_xy, t_x_xx_y_xz, t_x_xx_y_yy, t_x_xx_y_yz, t_x_xx_y_zz,\
                           t_x_xx_z_xx, t_x_xx_z_xy, t_x_xx_z_xz, t_x_xx_z_yy, t_x_xx_z_yz,\
                           t_x_xx_z_zz, t_x_xy_x_xx, t_x_xy_x_xy, t_x_xy_x_xz, t_x_xy_x_yy,\
                           t_x_xy_x_yz, t_x_xy_x_zz, t_x_xy_y_xx, t_x_xy_y_xy, t_x_xy_y_xz,\
                           t_x_xy_y_yy, t_x_xy_y_yz, t_x_xy_y_zz, t_x_xy_z_xx, t_x_xy_z_xy,\
                           t_x_xy_z_xz, t_x_xy_z_yy, t_x_xy_z_yz, t_x_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xy_z_zz[i] = t_0_xxy_z_zz[i] - rab_x[i] * t_0_xy_z_zz[i];

        t_x_xy_z_yz[i] = t_0_xxy_z_yz[i] - rab_x[i] * t_0_xy_z_yz[i];

        t_x_xy_z_yy[i] = t_0_xxy_z_yy[i] - rab_x[i] * t_0_xy_z_yy[i];

        t_x_xy_z_xz[i] = t_0_xxy_z_xz[i] - rab_x[i] * t_0_xy_z_xz[i];

        t_x_xy_z_xy[i] = t_0_xxy_z_xy[i] - rab_x[i] * t_0_xy_z_xy[i];

        t_x_xy_z_xx[i] = t_0_xxy_z_xx[i] - rab_x[i] * t_0_xy_z_xx[i];

        t_x_xy_y_zz[i] = t_0_xxy_y_zz[i] - rab_x[i] * t_0_xy_y_zz[i];

        t_x_xy_y_yz[i] = t_0_xxy_y_yz[i] - rab_x[i] * t_0_xy_y_yz[i];

        t_x_xy_y_yy[i] = t_0_xxy_y_yy[i] - rab_x[i] * t_0_xy_y_yy[i];

        t_x_xy_y_xz[i] = t_0_xxy_y_xz[i] - rab_x[i] * t_0_xy_y_xz[i];

        t_x_xy_y_xy[i] = t_0_xxy_y_xy[i] - rab_x[i] * t_0_xy_y_xy[i];

        t_x_xy_y_xx[i] = t_0_xxy_y_xx[i] - rab_x[i] * t_0_xy_y_xx[i];

        t_x_xy_x_zz[i] = t_0_xxy_x_zz[i] - rab_x[i] * t_0_xy_x_zz[i];

        t_x_xy_x_yz[i] = t_0_xxy_x_yz[i] - rab_x[i] * t_0_xy_x_yz[i];

        t_x_xy_x_yy[i] = t_0_xxy_x_yy[i] - rab_x[i] * t_0_xy_x_yy[i];

        t_x_xy_x_xz[i] = t_0_xxy_x_xz[i] - rab_x[i] * t_0_xy_x_xz[i];

        t_x_xy_x_xy[i] = t_0_xxy_x_xy[i] - rab_x[i] * t_0_xy_x_xy[i];

        t_x_xy_x_xx[i] = t_0_xxy_x_xx[i] - rab_x[i] * t_0_xy_x_xx[i];

        t_x_xx_z_zz[i] = t_0_xxx_z_zz[i] - rab_x[i] * t_0_xx_z_zz[i];

        t_x_xx_z_yz[i] = t_0_xxx_z_yz[i] - rab_x[i] * t_0_xx_z_yz[i];

        t_x_xx_z_yy[i] = t_0_xxx_z_yy[i] - rab_x[i] * t_0_xx_z_yy[i];

        t_x_xx_z_xz[i] = t_0_xxx_z_xz[i] - rab_x[i] * t_0_xx_z_xz[i];

        t_x_xx_z_xy[i] = t_0_xxx_z_xy[i] - rab_x[i] * t_0_xx_z_xy[i];

        t_x_xx_z_xx[i] = t_0_xxx_z_xx[i] - rab_x[i] * t_0_xx_z_xx[i];

        t_x_xx_y_zz[i] = t_0_xxx_y_zz[i] - rab_x[i] * t_0_xx_y_zz[i];

        t_x_xx_y_yz[i] = t_0_xxx_y_yz[i] - rab_x[i] * t_0_xx_y_yz[i];

        t_x_xx_y_yy[i] = t_0_xxx_y_yy[i] - rab_x[i] * t_0_xx_y_yy[i];

        t_x_xx_y_xz[i] = t_0_xxx_y_xz[i] - rab_x[i] * t_0_xx_y_xz[i];

        t_x_xx_y_xy[i] = t_0_xxx_y_xy[i] - rab_x[i] * t_0_xx_y_xy[i];

        t_x_xx_y_xx[i] = t_0_xxx_y_xx[i] - rab_x[i] * t_0_xx_y_xx[i];

        t_x_xx_x_zz[i] = t_0_xxx_x_zz[i] - rab_x[i] * t_0_xx_x_zz[i];

        t_x_xx_x_yz[i] = t_0_xxx_x_yz[i] - rab_x[i] * t_0_xx_x_yz[i];

        t_x_xx_x_yy[i] = t_0_xxx_x_yy[i] - rab_x[i] * t_0_xx_x_yy[i];

        t_x_xx_x_xz[i] = t_0_xxx_x_xz[i] - rab_x[i] * t_0_xx_x_xz[i];

        t_x_xx_x_xy[i] = t_0_xxx_x_xy[i] - rab_x[i] * t_0_xx_x_xy[i];

        t_x_xx_x_xx[i] = t_0_xxx_x_xx[i] - rab_x[i] * t_0_xx_x_xx[i];
    }
}


} // derirec namespace
