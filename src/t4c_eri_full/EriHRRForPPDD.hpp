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
compHostHRRForPPDD_V0(      BufferHostXY<T>&      intsBufferPPDD,
                      const BufferHostX<int32_t>& intsIndexesPPDD,
                      const BufferHostXY<T>&      intsBufferSPDD,
                      const BufferHostX<int32_t>& intsIndexesSPDD,
                      const BufferHostXY<T>&      intsBufferSDDD,
                      const BufferHostX<int32_t>& intsIndexesSDDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPDD) integral components

    t_z_z_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(0));

    t_z_z_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(1));

    t_z_z_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(2));

    t_z_z_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(3));

    t_z_z_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(4));

    t_z_z_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(5));

    t_z_z_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(6));

    t_z_z_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(7));

    t_z_z_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(8));

    t_z_z_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(9));

    t_z_z_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(10));

    t_z_z_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(11));

    t_z_z_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(12));

    t_z_z_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(13));

    t_z_z_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(14));

    t_z_z_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(15));

    t_z_z_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(16));

    t_z_z_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(17));

    t_z_z_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(18));

    t_z_z_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(19));

    t_z_z_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(20));

    t_z_z_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(21));

    t_z_z_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(22));

    t_z_z_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(23));

    t_z_z_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(24));

    t_z_z_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(25));

    t_z_z_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(26));

    t_z_z_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(27));

    t_z_z_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(28));

    t_z_z_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(29));

    t_z_z_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(30));

    t_z_z_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(31));

    t_z_z_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(32));

    t_z_z_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(33));

    t_z_z_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(34));

    t_z_z_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(35));

    t_z_y_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(36));

    t_z_y_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(37));

    t_z_y_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(38));

    t_z_y_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(39));

    t_z_y_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(40));

    t_z_y_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(41));

    t_z_y_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(42));

    t_z_y_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(43));

    t_z_y_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(44));

    t_z_y_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(45));

    t_z_y_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(46));

    t_z_y_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(47));

    t_z_y_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(48));

    t_z_y_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(49));

    t_z_y_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(50));

    t_z_y_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(51));

    t_z_y_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(52));

    t_z_y_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(53));

    t_z_y_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(54));

    t_z_y_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(55));

    t_z_y_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(56));

    t_z_y_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(57));

    t_z_y_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(58));

    t_z_y_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(59));

    t_z_y_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(60));

    t_z_y_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(61));

    t_z_y_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(62));

    t_z_y_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(63));

    t_z_y_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(64));

    t_z_y_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(65));

    t_z_y_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(66));

    t_z_y_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(67));

    t_z_y_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(68));

    t_z_y_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(69));

    t_z_y_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(70));

    t_z_y_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(71));

    t_z_x_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(72));

    t_z_x_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(73));

    t_z_x_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(74));

    t_z_x_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(75));

    t_z_x_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(76));

    t_z_x_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(77));

    t_z_x_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(78));

    t_z_x_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(79));

    t_z_x_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(80));

    t_z_x_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(81));

    t_z_x_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(82));

    t_z_x_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(83));

    t_z_x_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(84));

    t_z_x_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(85));

    t_z_x_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(86));

    t_z_x_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(87));

    t_z_x_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(88));

    t_z_x_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(89));

    t_z_x_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(90));

    t_z_x_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(91));

    t_z_x_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(92));

    t_z_x_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(93));

    t_z_x_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(94));

    t_z_x_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(95));

    t_z_x_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(96));

    t_z_x_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(97));

    t_z_x_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(98));

    t_z_x_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(99));

    t_z_x_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(100));

    t_z_x_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(101));

    t_z_x_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(102));

    t_z_x_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(103));

    t_z_x_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(104));

    t_z_x_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(105));

    t_z_x_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(106));

    t_z_x_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(107));

    t_y_z_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(108));

    t_y_z_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(109));

    t_y_z_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(110));

    t_y_z_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(111));

    t_y_z_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(112));

    t_y_z_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(113));

    t_y_z_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(114));

    t_y_z_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(115));

    t_y_z_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(116));

    t_y_z_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(117));

    t_y_z_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(118));

    t_y_z_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(119));

    t_y_z_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(120));

    t_y_z_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(121));

    t_y_z_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(122));

    t_y_z_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(123));

    t_y_z_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(124));

    t_y_z_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(125));

    t_y_z_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(126));

    t_y_z_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(127));

    t_y_z_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(128));

    t_y_z_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(129));

    t_y_z_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(130));

    t_y_z_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(131));

    t_y_z_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(132));

    t_y_z_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(133));

    t_y_z_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(134));

    t_y_z_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(135));

    t_y_z_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(136));

    t_y_z_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(137));

    t_y_z_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(138));

    t_y_z_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(139));

    t_y_z_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(140));

    t_y_z_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(141));

    t_y_z_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(142));

    t_y_z_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(143));

    t_y_y_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(144));

    t_y_y_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(145));

    t_y_y_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(146));

    t_y_y_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(147));

    t_y_y_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(148));

    t_y_y_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(149));

    t_y_y_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(150));

    t_y_y_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(151));

    t_y_y_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(152));

    t_y_y_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(153));

    t_y_y_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(154));

    t_y_y_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(155));

    t_y_y_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(156));

    t_y_y_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(157));

    t_y_y_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(158));

    t_y_y_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(159));

    t_y_y_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(160));

    t_y_y_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(161));

    t_y_y_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(162));

    t_y_y_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(163));

    t_y_y_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(164));

    t_y_y_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(165));

    t_y_y_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(166));

    t_y_y_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(167));

    t_y_y_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(168));

    t_y_y_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(169));

    t_y_y_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(170));

    t_y_y_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(171));

    t_y_y_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(172));

    t_y_y_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(173));

    t_y_y_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(174));

    t_y_y_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(175));

    t_y_y_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(176));

    t_y_y_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(177));

    t_y_y_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(178));

    t_y_y_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(179));

    t_y_x_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(180));

    t_y_x_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(181));

    t_y_x_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(182));

    t_y_x_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(183));

    t_y_x_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(184));

    t_y_x_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(185));

    t_y_x_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(186));

    t_y_x_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(187));

    t_y_x_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(188));

    t_y_x_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(189));

    t_y_x_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(190));

    t_y_x_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(191));

    t_y_x_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(192));

    t_y_x_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(193));

    t_y_x_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(194));

    t_y_x_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(195));

    t_y_x_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(196));

    t_y_x_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(197));

    t_y_x_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(198));

    t_y_x_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(199));

    t_y_x_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(200));

    t_y_x_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(201));

    t_y_x_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(202));

    t_y_x_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(203));

    t_y_x_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(204));

    t_y_x_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(205));

    t_y_x_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(206));

    t_y_x_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(207));

    t_y_x_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(208));

    t_y_x_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(209));

    t_y_x_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(210));

    t_y_x_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(211));

    t_y_x_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(212));

    t_y_x_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(213));

    t_y_x_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(214));

    t_y_x_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(215));

    t_x_z_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(216));

    t_x_z_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(217));

    t_x_z_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(218));

    t_x_z_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(219));

    t_x_z_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(220));

    t_x_z_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(221));

    t_x_z_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(222));

    t_x_z_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(223));

    t_x_z_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(224));

    t_x_z_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(225));

    t_x_z_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(226));

    t_x_z_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(227));

    t_x_z_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(228));

    t_x_z_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(229));

    t_x_z_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(230));

    t_x_z_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(231));

    t_x_z_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(232));

    t_x_z_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(233));

    t_x_z_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(234));

    t_x_z_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(235));

    t_x_z_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(236));

    t_x_z_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(237));

    t_x_z_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(238));

    t_x_z_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(239));

    t_x_z_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(240));

    t_x_z_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(241));

    t_x_z_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(242));

    t_x_z_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(243));

    t_x_z_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(244));

    t_x_z_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(245));

    t_x_z_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(246));

    t_x_z_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(247));

    t_x_z_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(248));

    t_x_z_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(249));

    t_x_z_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(250));

    t_x_z_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(251));

    t_x_y_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(252));

    t_x_y_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(253));

    t_x_y_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(254));

    t_x_y_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(255));

    t_x_y_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(256));

    t_x_y_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(257));

    t_x_y_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(258));

    t_x_y_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(259));

    t_x_y_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(260));

    t_x_y_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(261));

    t_x_y_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(262));

    t_x_y_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(263));

    t_x_y_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(264));

    t_x_y_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(265));

    t_x_y_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(266));

    t_x_y_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(267));

    t_x_y_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(268));

    t_x_y_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(269));

    t_x_y_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(270));

    t_x_y_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(271));

    t_x_y_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(272));

    t_x_y_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(273));

    t_x_y_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(274));

    t_x_y_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(275));

    t_x_y_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(276));

    t_x_y_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(277));

    t_x_y_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(278));

    t_x_y_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(279));

    t_x_y_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(280));

    t_x_y_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(281));

    t_x_y_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(282));

    t_x_y_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(283));

    t_x_y_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(284));

    t_x_y_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(285));

    t_x_y_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(286));

    t_x_y_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(287));

    t_x_x_zz_zz = intsBufferPPDD.data(intsIndexesPPDD(288));

    t_x_x_zz_yz = intsBufferPPDD.data(intsIndexesPPDD(289));

    t_x_x_zz_yy = intsBufferPPDD.data(intsIndexesPPDD(290));

    t_x_x_zz_xz = intsBufferPPDD.data(intsIndexesPPDD(291));

    t_x_x_zz_xy = intsBufferPPDD.data(intsIndexesPPDD(292));

    t_x_x_zz_xx = intsBufferPPDD.data(intsIndexesPPDD(293));

    t_x_x_yz_zz = intsBufferPPDD.data(intsIndexesPPDD(294));

    t_x_x_yz_yz = intsBufferPPDD.data(intsIndexesPPDD(295));

    t_x_x_yz_yy = intsBufferPPDD.data(intsIndexesPPDD(296));

    t_x_x_yz_xz = intsBufferPPDD.data(intsIndexesPPDD(297));

    t_x_x_yz_xy = intsBufferPPDD.data(intsIndexesPPDD(298));

    t_x_x_yz_xx = intsBufferPPDD.data(intsIndexesPPDD(299));

    t_x_x_yy_zz = intsBufferPPDD.data(intsIndexesPPDD(300));

    t_x_x_yy_yz = intsBufferPPDD.data(intsIndexesPPDD(301));

    t_x_x_yy_yy = intsBufferPPDD.data(intsIndexesPPDD(302));

    t_x_x_yy_xz = intsBufferPPDD.data(intsIndexesPPDD(303));

    t_x_x_yy_xy = intsBufferPPDD.data(intsIndexesPPDD(304));

    t_x_x_yy_xx = intsBufferPPDD.data(intsIndexesPPDD(305));

    t_x_x_xz_zz = intsBufferPPDD.data(intsIndexesPPDD(306));

    t_x_x_xz_yz = intsBufferPPDD.data(intsIndexesPPDD(307));

    t_x_x_xz_yy = intsBufferPPDD.data(intsIndexesPPDD(308));

    t_x_x_xz_xz = intsBufferPPDD.data(intsIndexesPPDD(309));

    t_x_x_xz_xy = intsBufferPPDD.data(intsIndexesPPDD(310));

    t_x_x_xz_xx = intsBufferPPDD.data(intsIndexesPPDD(311));

    t_x_x_xy_zz = intsBufferPPDD.data(intsIndexesPPDD(312));

    t_x_x_xy_yz = intsBufferPPDD.data(intsIndexesPPDD(313));

    t_x_x_xy_yy = intsBufferPPDD.data(intsIndexesPPDD(314));

    t_x_x_xy_xz = intsBufferPPDD.data(intsIndexesPPDD(315));

    t_x_x_xy_xy = intsBufferPPDD.data(intsIndexesPPDD(316));

    t_x_x_xy_xx = intsBufferPPDD.data(intsIndexesPPDD(317));

    t_x_x_xx_zz = intsBufferPPDD.data(intsIndexesPPDD(318));

    t_x_x_xx_yz = intsBufferPPDD.data(intsIndexesPPDD(319));

    t_x_x_xx_yy = intsBufferPPDD.data(intsIndexesPPDD(320));

    t_x_x_xx_xz = intsBufferPPDD.data(intsIndexesPPDD(321));

    t_x_x_xx_xy = intsBufferPPDD.data(intsIndexesPPDD(322));

    t_x_x_xx_xx = intsBufferPPDD.data(intsIndexesPPDD(323));

    // set up (SPDD) integral components

    t_0_z_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(0));

    t_0_z_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(1));

    t_0_z_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(2));

    t_0_z_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(3));

    t_0_z_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(4));

    t_0_z_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(5));

    t_0_z_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(6));

    t_0_z_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(7));

    t_0_z_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(8));

    t_0_z_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(9));

    t_0_z_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(10));

    t_0_z_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(11));

    t_0_z_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(12));

    t_0_z_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(13));

    t_0_z_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(14));

    t_0_z_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(15));

    t_0_z_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(16));

    t_0_z_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(17));

    t_0_z_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(18));

    t_0_z_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(19));

    t_0_z_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(20));

    t_0_z_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(21));

    t_0_z_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(22));

    t_0_z_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(23));

    t_0_z_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(24));

    t_0_z_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(25));

    t_0_z_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(26));

    t_0_z_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(27));

    t_0_z_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(28));

    t_0_z_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(29));

    t_0_z_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(30));

    t_0_z_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(31));

    t_0_z_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(32));

    t_0_z_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(33));

    t_0_z_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(34));

    t_0_z_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(35));

    t_0_y_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(36));

    t_0_y_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(37));

    t_0_y_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(38));

    t_0_y_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(39));

    t_0_y_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(40));

    t_0_y_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(41));

    t_0_y_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(42));

    t_0_y_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(43));

    t_0_y_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(44));

    t_0_y_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(45));

    t_0_y_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(46));

    t_0_y_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(47));

    t_0_y_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(48));

    t_0_y_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(49));

    t_0_y_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(50));

    t_0_y_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(51));

    t_0_y_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(52));

    t_0_y_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(53));

    t_0_y_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(54));

    t_0_y_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(55));

    t_0_y_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(56));

    t_0_y_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(57));

    t_0_y_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(58));

    t_0_y_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(59));

    t_0_y_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(60));

    t_0_y_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(61));

    t_0_y_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(62));

    t_0_y_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(63));

    t_0_y_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(64));

    t_0_y_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(65));

    t_0_y_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(66));

    t_0_y_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(67));

    t_0_y_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(68));

    t_0_y_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(69));

    t_0_y_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(70));

    t_0_y_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(71));

    t_0_x_zz_zz = intsBufferSPDD.data(intsIndexesSPDD(72));

    t_0_x_zz_yz = intsBufferSPDD.data(intsIndexesSPDD(73));

    t_0_x_zz_yy = intsBufferSPDD.data(intsIndexesSPDD(74));

    t_0_x_zz_xz = intsBufferSPDD.data(intsIndexesSPDD(75));

    t_0_x_zz_xy = intsBufferSPDD.data(intsIndexesSPDD(76));

    t_0_x_zz_xx = intsBufferSPDD.data(intsIndexesSPDD(77));

    t_0_x_yz_zz = intsBufferSPDD.data(intsIndexesSPDD(78));

    t_0_x_yz_yz = intsBufferSPDD.data(intsIndexesSPDD(79));

    t_0_x_yz_yy = intsBufferSPDD.data(intsIndexesSPDD(80));

    t_0_x_yz_xz = intsBufferSPDD.data(intsIndexesSPDD(81));

    t_0_x_yz_xy = intsBufferSPDD.data(intsIndexesSPDD(82));

    t_0_x_yz_xx = intsBufferSPDD.data(intsIndexesSPDD(83));

    t_0_x_yy_zz = intsBufferSPDD.data(intsIndexesSPDD(84));

    t_0_x_yy_yz = intsBufferSPDD.data(intsIndexesSPDD(85));

    t_0_x_yy_yy = intsBufferSPDD.data(intsIndexesSPDD(86));

    t_0_x_yy_xz = intsBufferSPDD.data(intsIndexesSPDD(87));

    t_0_x_yy_xy = intsBufferSPDD.data(intsIndexesSPDD(88));

    t_0_x_yy_xx = intsBufferSPDD.data(intsIndexesSPDD(89));

    t_0_x_xz_zz = intsBufferSPDD.data(intsIndexesSPDD(90));

    t_0_x_xz_yz = intsBufferSPDD.data(intsIndexesSPDD(91));

    t_0_x_xz_yy = intsBufferSPDD.data(intsIndexesSPDD(92));

    t_0_x_xz_xz = intsBufferSPDD.data(intsIndexesSPDD(93));

    t_0_x_xz_xy = intsBufferSPDD.data(intsIndexesSPDD(94));

    t_0_x_xz_xx = intsBufferSPDD.data(intsIndexesSPDD(95));

    t_0_x_xy_zz = intsBufferSPDD.data(intsIndexesSPDD(96));

    t_0_x_xy_yz = intsBufferSPDD.data(intsIndexesSPDD(97));

    t_0_x_xy_yy = intsBufferSPDD.data(intsIndexesSPDD(98));

    t_0_x_xy_xz = intsBufferSPDD.data(intsIndexesSPDD(99));

    t_0_x_xy_xy = intsBufferSPDD.data(intsIndexesSPDD(100));

    t_0_x_xy_xx = intsBufferSPDD.data(intsIndexesSPDD(101));

    t_0_x_xx_zz = intsBufferSPDD.data(intsIndexesSPDD(102));

    t_0_x_xx_yz = intsBufferSPDD.data(intsIndexesSPDD(103));

    t_0_x_xx_yy = intsBufferSPDD.data(intsIndexesSPDD(104));

    t_0_x_xx_xz = intsBufferSPDD.data(intsIndexesSPDD(105));

    t_0_x_xx_xy = intsBufferSPDD.data(intsIndexesSPDD(106));

    t_0_x_xx_xx = intsBufferSPDD.data(intsIndexesSPDD(107));

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

    #pragma omp simd align(rab_z, t_0_z_xx_xx, t_0_z_xx_xy, t_0_z_xx_xz, t_0_z_xx_yy,\
                           t_0_z_xx_yz, t_0_z_xx_zz, t_0_z_xy_xx, t_0_z_xy_xy, t_0_z_xy_xz,\
                           t_0_z_xy_yy, t_0_z_xy_yz, t_0_z_xy_zz, t_0_z_xz_xx, t_0_z_xz_xy,\
                           t_0_z_xz_xz, t_0_z_xz_yy, t_0_z_xz_yz, t_0_z_xz_zz, t_0_z_yy_xx,\
                           t_0_z_yy_xy, t_0_z_yy_xz, t_0_z_yy_yy, t_0_z_yy_yz, t_0_z_yy_zz,\
                           t_0_z_yz_xx, t_0_z_yz_xy, t_0_z_yz_xz, t_0_z_yz_yy, t_0_z_yz_yz,\
                           t_0_z_yz_zz, t_0_z_zz_xx, t_0_z_zz_xy, t_0_z_zz_xz, t_0_z_zz_yy,\
                           t_0_z_zz_yz, t_0_z_zz_zz, t_0_zz_xx_xx, t_0_zz_xx_xy, t_0_zz_xx_xz,\
                           t_0_zz_xx_yy, t_0_zz_xx_yz, t_0_zz_xx_zz, t_0_zz_xy_xx, t_0_zz_xy_xy,\
                           t_0_zz_xy_xz, t_0_zz_xy_yy, t_0_zz_xy_yz, t_0_zz_xy_zz, t_0_zz_xz_xx,\
                           t_0_zz_xz_xy, t_0_zz_xz_xz, t_0_zz_xz_yy, t_0_zz_xz_yz, t_0_zz_xz_zz,\
                           t_0_zz_yy_xx, t_0_zz_yy_xy, t_0_zz_yy_xz, t_0_zz_yy_yy, t_0_zz_yy_yz,\
                           t_0_zz_yy_zz, t_0_zz_yz_xx, t_0_zz_yz_xy, t_0_zz_yz_xz, t_0_zz_yz_yy,\
                           t_0_zz_yz_yz, t_0_zz_yz_zz, t_0_zz_zz_xx, t_0_zz_zz_xy, t_0_zz_zz_xz,\
                           t_0_zz_zz_yy, t_0_zz_zz_yz, t_0_zz_zz_zz, t_z_z_xx_xx, t_z_z_xx_xy,\
                           t_z_z_xx_xz, t_z_z_xx_yy, t_z_z_xx_yz, t_z_z_xx_zz, t_z_z_xy_xx,\
                           t_z_z_xy_xy, t_z_z_xy_xz, t_z_z_xy_yy, t_z_z_xy_yz, t_z_z_xy_zz,\
                           t_z_z_xz_xx, t_z_z_xz_xy, t_z_z_xz_xz, t_z_z_xz_yy, t_z_z_xz_yz,\
                           t_z_z_xz_zz, t_z_z_yy_xx, t_z_z_yy_xy, t_z_z_yy_xz, t_z_z_yy_yy,\
                           t_z_z_yy_yz, t_z_z_yy_zz, t_z_z_yz_xx, t_z_z_yz_xy, t_z_z_yz_xz,\
                           t_z_z_yz_yy, t_z_z_yz_yz, t_z_z_yz_zz, t_z_z_zz_xx, t_z_z_zz_xy,\
                           t_z_z_zz_xz, t_z_z_zz_yy, t_z_z_zz_yz, t_z_z_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_zz_zz[i] = t_0_zz_zz_zz[i] - rab_z[i] * t_0_z_zz_zz[i];

        t_z_z_zz_yz[i] = t_0_zz_zz_yz[i] - rab_z[i] * t_0_z_zz_yz[i];

        t_z_z_zz_yy[i] = t_0_zz_zz_yy[i] - rab_z[i] * t_0_z_zz_yy[i];

        t_z_z_zz_xz[i] = t_0_zz_zz_xz[i] - rab_z[i] * t_0_z_zz_xz[i];

        t_z_z_zz_xy[i] = t_0_zz_zz_xy[i] - rab_z[i] * t_0_z_zz_xy[i];

        t_z_z_zz_xx[i] = t_0_zz_zz_xx[i] - rab_z[i] * t_0_z_zz_xx[i];

        t_z_z_yz_zz[i] = t_0_zz_yz_zz[i] - rab_z[i] * t_0_z_yz_zz[i];

        t_z_z_yz_yz[i] = t_0_zz_yz_yz[i] - rab_z[i] * t_0_z_yz_yz[i];

        t_z_z_yz_yy[i] = t_0_zz_yz_yy[i] - rab_z[i] * t_0_z_yz_yy[i];

        t_z_z_yz_xz[i] = t_0_zz_yz_xz[i] - rab_z[i] * t_0_z_yz_xz[i];

        t_z_z_yz_xy[i] = t_0_zz_yz_xy[i] - rab_z[i] * t_0_z_yz_xy[i];

        t_z_z_yz_xx[i] = t_0_zz_yz_xx[i] - rab_z[i] * t_0_z_yz_xx[i];

        t_z_z_yy_zz[i] = t_0_zz_yy_zz[i] - rab_z[i] * t_0_z_yy_zz[i];

        t_z_z_yy_yz[i] = t_0_zz_yy_yz[i] - rab_z[i] * t_0_z_yy_yz[i];

        t_z_z_yy_yy[i] = t_0_zz_yy_yy[i] - rab_z[i] * t_0_z_yy_yy[i];

        t_z_z_yy_xz[i] = t_0_zz_yy_xz[i] - rab_z[i] * t_0_z_yy_xz[i];

        t_z_z_yy_xy[i] = t_0_zz_yy_xy[i] - rab_z[i] * t_0_z_yy_xy[i];

        t_z_z_yy_xx[i] = t_0_zz_yy_xx[i] - rab_z[i] * t_0_z_yy_xx[i];

        t_z_z_xz_zz[i] = t_0_zz_xz_zz[i] - rab_z[i] * t_0_z_xz_zz[i];

        t_z_z_xz_yz[i] = t_0_zz_xz_yz[i] - rab_z[i] * t_0_z_xz_yz[i];

        t_z_z_xz_yy[i] = t_0_zz_xz_yy[i] - rab_z[i] * t_0_z_xz_yy[i];

        t_z_z_xz_xz[i] = t_0_zz_xz_xz[i] - rab_z[i] * t_0_z_xz_xz[i];

        t_z_z_xz_xy[i] = t_0_zz_xz_xy[i] - rab_z[i] * t_0_z_xz_xy[i];

        t_z_z_xz_xx[i] = t_0_zz_xz_xx[i] - rab_z[i] * t_0_z_xz_xx[i];

        t_z_z_xy_zz[i] = t_0_zz_xy_zz[i] - rab_z[i] * t_0_z_xy_zz[i];

        t_z_z_xy_yz[i] = t_0_zz_xy_yz[i] - rab_z[i] * t_0_z_xy_yz[i];

        t_z_z_xy_yy[i] = t_0_zz_xy_yy[i] - rab_z[i] * t_0_z_xy_yy[i];

        t_z_z_xy_xz[i] = t_0_zz_xy_xz[i] - rab_z[i] * t_0_z_xy_xz[i];

        t_z_z_xy_xy[i] = t_0_zz_xy_xy[i] - rab_z[i] * t_0_z_xy_xy[i];

        t_z_z_xy_xx[i] = t_0_zz_xy_xx[i] - rab_z[i] * t_0_z_xy_xx[i];

        t_z_z_xx_zz[i] = t_0_zz_xx_zz[i] - rab_z[i] * t_0_z_xx_zz[i];

        t_z_z_xx_yz[i] = t_0_zz_xx_yz[i] - rab_z[i] * t_0_z_xx_yz[i];

        t_z_z_xx_yy[i] = t_0_zz_xx_yy[i] - rab_z[i] * t_0_z_xx_yy[i];

        t_z_z_xx_xz[i] = t_0_zz_xx_xz[i] - rab_z[i] * t_0_z_xx_xz[i];

        t_z_z_xx_xy[i] = t_0_zz_xx_xy[i] - rab_z[i] * t_0_z_xx_xy[i];

        t_z_z_xx_xx[i] = t_0_zz_xx_xx[i] - rab_z[i] * t_0_z_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_y_xx_xx, t_0_y_xx_xy, t_0_y_xx_xz, t_0_y_xx_yy,\
                           t_0_y_xx_yz, t_0_y_xx_zz, t_0_y_xy_xx, t_0_y_xy_xy, t_0_y_xy_xz,\
                           t_0_y_xy_yy, t_0_y_xy_yz, t_0_y_xy_zz, t_0_y_xz_xx, t_0_y_xz_xy,\
                           t_0_y_xz_xz, t_0_y_xz_yy, t_0_y_xz_yz, t_0_y_xz_zz, t_0_y_yy_xx,\
                           t_0_y_yy_xy, t_0_y_yy_xz, t_0_y_yy_yy, t_0_y_yy_yz, t_0_y_yy_zz,\
                           t_0_y_yz_xx, t_0_y_yz_xy, t_0_y_yz_xz, t_0_y_yz_yy, t_0_y_yz_yz,\
                           t_0_y_yz_zz, t_0_y_zz_xx, t_0_y_zz_xy, t_0_y_zz_xz, t_0_y_zz_yy,\
                           t_0_y_zz_yz, t_0_y_zz_zz, t_0_yz_xx_xx, t_0_yz_xx_xy, t_0_yz_xx_xz,\
                           t_0_yz_xx_yy, t_0_yz_xx_yz, t_0_yz_xx_zz, t_0_yz_xy_xx, t_0_yz_xy_xy,\
                           t_0_yz_xy_xz, t_0_yz_xy_yy, t_0_yz_xy_yz, t_0_yz_xy_zz, t_0_yz_xz_xx,\
                           t_0_yz_xz_xy, t_0_yz_xz_xz, t_0_yz_xz_yy, t_0_yz_xz_yz, t_0_yz_xz_zz,\
                           t_0_yz_yy_xx, t_0_yz_yy_xy, t_0_yz_yy_xz, t_0_yz_yy_yy, t_0_yz_yy_yz,\
                           t_0_yz_yy_zz, t_0_yz_yz_xx, t_0_yz_yz_xy, t_0_yz_yz_xz, t_0_yz_yz_yy,\
                           t_0_yz_yz_yz, t_0_yz_yz_zz, t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz,\
                           t_0_yz_zz_yy, t_0_yz_zz_yz, t_0_yz_zz_zz, t_z_y_xx_xx, t_z_y_xx_xy,\
                           t_z_y_xx_xz, t_z_y_xx_yy, t_z_y_xx_yz, t_z_y_xx_zz, t_z_y_xy_xx,\
                           t_z_y_xy_xy, t_z_y_xy_xz, t_z_y_xy_yy, t_z_y_xy_yz, t_z_y_xy_zz,\
                           t_z_y_xz_xx, t_z_y_xz_xy, t_z_y_xz_xz, t_z_y_xz_yy, t_z_y_xz_yz,\
                           t_z_y_xz_zz, t_z_y_yy_xx, t_z_y_yy_xy, t_z_y_yy_xz, t_z_y_yy_yy,\
                           t_z_y_yy_yz, t_z_y_yy_zz, t_z_y_yz_xx, t_z_y_yz_xy, t_z_y_yz_xz,\
                           t_z_y_yz_yy, t_z_y_yz_yz, t_z_y_yz_zz, t_z_y_zz_xx, t_z_y_zz_xy,\
                           t_z_y_zz_xz, t_z_y_zz_yy, t_z_y_zz_yz, t_z_y_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_y_zz_zz[i] = t_0_yz_zz_zz[i] - rab_z[i] * t_0_y_zz_zz[i];

        t_z_y_zz_yz[i] = t_0_yz_zz_yz[i] - rab_z[i] * t_0_y_zz_yz[i];

        t_z_y_zz_yy[i] = t_0_yz_zz_yy[i] - rab_z[i] * t_0_y_zz_yy[i];

        t_z_y_zz_xz[i] = t_0_yz_zz_xz[i] - rab_z[i] * t_0_y_zz_xz[i];

        t_z_y_zz_xy[i] = t_0_yz_zz_xy[i] - rab_z[i] * t_0_y_zz_xy[i];

        t_z_y_zz_xx[i] = t_0_yz_zz_xx[i] - rab_z[i] * t_0_y_zz_xx[i];

        t_z_y_yz_zz[i] = t_0_yz_yz_zz[i] - rab_z[i] * t_0_y_yz_zz[i];

        t_z_y_yz_yz[i] = t_0_yz_yz_yz[i] - rab_z[i] * t_0_y_yz_yz[i];

        t_z_y_yz_yy[i] = t_0_yz_yz_yy[i] - rab_z[i] * t_0_y_yz_yy[i];

        t_z_y_yz_xz[i] = t_0_yz_yz_xz[i] - rab_z[i] * t_0_y_yz_xz[i];

        t_z_y_yz_xy[i] = t_0_yz_yz_xy[i] - rab_z[i] * t_0_y_yz_xy[i];

        t_z_y_yz_xx[i] = t_0_yz_yz_xx[i] - rab_z[i] * t_0_y_yz_xx[i];

        t_z_y_yy_zz[i] = t_0_yz_yy_zz[i] - rab_z[i] * t_0_y_yy_zz[i];

        t_z_y_yy_yz[i] = t_0_yz_yy_yz[i] - rab_z[i] * t_0_y_yy_yz[i];

        t_z_y_yy_yy[i] = t_0_yz_yy_yy[i] - rab_z[i] * t_0_y_yy_yy[i];

        t_z_y_yy_xz[i] = t_0_yz_yy_xz[i] - rab_z[i] * t_0_y_yy_xz[i];

        t_z_y_yy_xy[i] = t_0_yz_yy_xy[i] - rab_z[i] * t_0_y_yy_xy[i];

        t_z_y_yy_xx[i] = t_0_yz_yy_xx[i] - rab_z[i] * t_0_y_yy_xx[i];

        t_z_y_xz_zz[i] = t_0_yz_xz_zz[i] - rab_z[i] * t_0_y_xz_zz[i];

        t_z_y_xz_yz[i] = t_0_yz_xz_yz[i] - rab_z[i] * t_0_y_xz_yz[i];

        t_z_y_xz_yy[i] = t_0_yz_xz_yy[i] - rab_z[i] * t_0_y_xz_yy[i];

        t_z_y_xz_xz[i] = t_0_yz_xz_xz[i] - rab_z[i] * t_0_y_xz_xz[i];

        t_z_y_xz_xy[i] = t_0_yz_xz_xy[i] - rab_z[i] * t_0_y_xz_xy[i];

        t_z_y_xz_xx[i] = t_0_yz_xz_xx[i] - rab_z[i] * t_0_y_xz_xx[i];

        t_z_y_xy_zz[i] = t_0_yz_xy_zz[i] - rab_z[i] * t_0_y_xy_zz[i];

        t_z_y_xy_yz[i] = t_0_yz_xy_yz[i] - rab_z[i] * t_0_y_xy_yz[i];

        t_z_y_xy_yy[i] = t_0_yz_xy_yy[i] - rab_z[i] * t_0_y_xy_yy[i];

        t_z_y_xy_xz[i] = t_0_yz_xy_xz[i] - rab_z[i] * t_0_y_xy_xz[i];

        t_z_y_xy_xy[i] = t_0_yz_xy_xy[i] - rab_z[i] * t_0_y_xy_xy[i];

        t_z_y_xy_xx[i] = t_0_yz_xy_xx[i] - rab_z[i] * t_0_y_xy_xx[i];

        t_z_y_xx_zz[i] = t_0_yz_xx_zz[i] - rab_z[i] * t_0_y_xx_zz[i];

        t_z_y_xx_yz[i] = t_0_yz_xx_yz[i] - rab_z[i] * t_0_y_xx_yz[i];

        t_z_y_xx_yy[i] = t_0_yz_xx_yy[i] - rab_z[i] * t_0_y_xx_yy[i];

        t_z_y_xx_xz[i] = t_0_yz_xx_xz[i] - rab_z[i] * t_0_y_xx_xz[i];

        t_z_y_xx_xy[i] = t_0_yz_xx_xy[i] - rab_z[i] * t_0_y_xx_xy[i];

        t_z_y_xx_xx[i] = t_0_yz_xx_xx[i] - rab_z[i] * t_0_y_xx_xx[i];
    }

    #pragma omp simd align(rab_z, t_0_x_xx_xx, t_0_x_xx_xy, t_0_x_xx_xz, t_0_x_xx_yy,\
                           t_0_x_xx_yz, t_0_x_xx_zz, t_0_x_xy_xx, t_0_x_xy_xy, t_0_x_xy_xz,\
                           t_0_x_xy_yy, t_0_x_xy_yz, t_0_x_xy_zz, t_0_x_xz_xx, t_0_x_xz_xy,\
                           t_0_x_xz_xz, t_0_x_xz_yy, t_0_x_xz_yz, t_0_x_xz_zz, t_0_x_yy_xx,\
                           t_0_x_yy_xy, t_0_x_yy_xz, t_0_x_yy_yy, t_0_x_yy_yz, t_0_x_yy_zz,\
                           t_0_x_yz_xx, t_0_x_yz_xy, t_0_x_yz_xz, t_0_x_yz_yy, t_0_x_yz_yz,\
                           t_0_x_yz_zz, t_0_x_zz_xx, t_0_x_zz_xy, t_0_x_zz_xz, t_0_x_zz_yy,\
                           t_0_x_zz_yz, t_0_x_zz_zz, t_0_xz_xx_xx, t_0_xz_xx_xy, t_0_xz_xx_xz,\
                           t_0_xz_xx_yy, t_0_xz_xx_yz, t_0_xz_xx_zz, t_0_xz_xy_xx, t_0_xz_xy_xy,\
                           t_0_xz_xy_xz, t_0_xz_xy_yy, t_0_xz_xy_yz, t_0_xz_xy_zz, t_0_xz_xz_xx,\
                           t_0_xz_xz_xy, t_0_xz_xz_xz, t_0_xz_xz_yy, t_0_xz_xz_yz, t_0_xz_xz_zz,\
                           t_0_xz_yy_xx, t_0_xz_yy_xy, t_0_xz_yy_xz, t_0_xz_yy_yy, t_0_xz_yy_yz,\
                           t_0_xz_yy_zz, t_0_xz_yz_xx, t_0_xz_yz_xy, t_0_xz_yz_xz, t_0_xz_yz_yy,\
                           t_0_xz_yz_yz, t_0_xz_yz_zz, t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz,\
                           t_0_xz_zz_yy, t_0_xz_zz_yz, t_0_xz_zz_zz, t_z_x_xx_xx, t_z_x_xx_xy,\
                           t_z_x_xx_xz, t_z_x_xx_yy, t_z_x_xx_yz, t_z_x_xx_zz, t_z_x_xy_xx,\
                           t_z_x_xy_xy, t_z_x_xy_xz, t_z_x_xy_yy, t_z_x_xy_yz, t_z_x_xy_zz,\
                           t_z_x_xz_xx, t_z_x_xz_xy, t_z_x_xz_xz, t_z_x_xz_yy, t_z_x_xz_yz,\
                           t_z_x_xz_zz, t_z_x_yy_xx, t_z_x_yy_xy, t_z_x_yy_xz, t_z_x_yy_yy,\
                           t_z_x_yy_yz, t_z_x_yy_zz, t_z_x_yz_xx, t_z_x_yz_xy, t_z_x_yz_xz,\
                           t_z_x_yz_yy, t_z_x_yz_yz, t_z_x_yz_zz, t_z_x_zz_xx, t_z_x_zz_xy,\
                           t_z_x_zz_xz, t_z_x_zz_yy, t_z_x_zz_yz, t_z_x_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_x_zz_zz[i] = t_0_xz_zz_zz[i] - rab_z[i] * t_0_x_zz_zz[i];

        t_z_x_zz_yz[i] = t_0_xz_zz_yz[i] - rab_z[i] * t_0_x_zz_yz[i];

        t_z_x_zz_yy[i] = t_0_xz_zz_yy[i] - rab_z[i] * t_0_x_zz_yy[i];

        t_z_x_zz_xz[i] = t_0_xz_zz_xz[i] - rab_z[i] * t_0_x_zz_xz[i];

        t_z_x_zz_xy[i] = t_0_xz_zz_xy[i] - rab_z[i] * t_0_x_zz_xy[i];

        t_z_x_zz_xx[i] = t_0_xz_zz_xx[i] - rab_z[i] * t_0_x_zz_xx[i];

        t_z_x_yz_zz[i] = t_0_xz_yz_zz[i] - rab_z[i] * t_0_x_yz_zz[i];

        t_z_x_yz_yz[i] = t_0_xz_yz_yz[i] - rab_z[i] * t_0_x_yz_yz[i];

        t_z_x_yz_yy[i] = t_0_xz_yz_yy[i] - rab_z[i] * t_0_x_yz_yy[i];

        t_z_x_yz_xz[i] = t_0_xz_yz_xz[i] - rab_z[i] * t_0_x_yz_xz[i];

        t_z_x_yz_xy[i] = t_0_xz_yz_xy[i] - rab_z[i] * t_0_x_yz_xy[i];

        t_z_x_yz_xx[i] = t_0_xz_yz_xx[i] - rab_z[i] * t_0_x_yz_xx[i];

        t_z_x_yy_zz[i] = t_0_xz_yy_zz[i] - rab_z[i] * t_0_x_yy_zz[i];

        t_z_x_yy_yz[i] = t_0_xz_yy_yz[i] - rab_z[i] * t_0_x_yy_yz[i];

        t_z_x_yy_yy[i] = t_0_xz_yy_yy[i] - rab_z[i] * t_0_x_yy_yy[i];

        t_z_x_yy_xz[i] = t_0_xz_yy_xz[i] - rab_z[i] * t_0_x_yy_xz[i];

        t_z_x_yy_xy[i] = t_0_xz_yy_xy[i] - rab_z[i] * t_0_x_yy_xy[i];

        t_z_x_yy_xx[i] = t_0_xz_yy_xx[i] - rab_z[i] * t_0_x_yy_xx[i];

        t_z_x_xz_zz[i] = t_0_xz_xz_zz[i] - rab_z[i] * t_0_x_xz_zz[i];

        t_z_x_xz_yz[i] = t_0_xz_xz_yz[i] - rab_z[i] * t_0_x_xz_yz[i];

        t_z_x_xz_yy[i] = t_0_xz_xz_yy[i] - rab_z[i] * t_0_x_xz_yy[i];

        t_z_x_xz_xz[i] = t_0_xz_xz_xz[i] - rab_z[i] * t_0_x_xz_xz[i];

        t_z_x_xz_xy[i] = t_0_xz_xz_xy[i] - rab_z[i] * t_0_x_xz_xy[i];

        t_z_x_xz_xx[i] = t_0_xz_xz_xx[i] - rab_z[i] * t_0_x_xz_xx[i];

        t_z_x_xy_zz[i] = t_0_xz_xy_zz[i] - rab_z[i] * t_0_x_xy_zz[i];

        t_z_x_xy_yz[i] = t_0_xz_xy_yz[i] - rab_z[i] * t_0_x_xy_yz[i];

        t_z_x_xy_yy[i] = t_0_xz_xy_yy[i] - rab_z[i] * t_0_x_xy_yy[i];

        t_z_x_xy_xz[i] = t_0_xz_xy_xz[i] - rab_z[i] * t_0_x_xy_xz[i];

        t_z_x_xy_xy[i] = t_0_xz_xy_xy[i] - rab_z[i] * t_0_x_xy_xy[i];

        t_z_x_xy_xx[i] = t_0_xz_xy_xx[i] - rab_z[i] * t_0_x_xy_xx[i];

        t_z_x_xx_zz[i] = t_0_xz_xx_zz[i] - rab_z[i] * t_0_x_xx_zz[i];

        t_z_x_xx_yz[i] = t_0_xz_xx_yz[i] - rab_z[i] * t_0_x_xx_yz[i];

        t_z_x_xx_yy[i] = t_0_xz_xx_yy[i] - rab_z[i] * t_0_x_xx_yy[i];

        t_z_x_xx_xz[i] = t_0_xz_xx_xz[i] - rab_z[i] * t_0_x_xx_xz[i];

        t_z_x_xx_xy[i] = t_0_xz_xx_xy[i] - rab_z[i] * t_0_x_xx_xy[i];

        t_z_x_xx_xx[i] = t_0_xz_xx_xx[i] - rab_z[i] * t_0_x_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_yz_xx_xx, t_0_yz_xx_xy, t_0_yz_xx_xz, t_0_yz_xx_yy,\
                           t_0_yz_xx_yz, t_0_yz_xx_zz, t_0_yz_xy_xx, t_0_yz_xy_xy, t_0_yz_xy_xz,\
                           t_0_yz_xy_yy, t_0_yz_xy_yz, t_0_yz_xy_zz, t_0_yz_xz_xx, t_0_yz_xz_xy,\
                           t_0_yz_xz_xz, t_0_yz_xz_yy, t_0_yz_xz_yz, t_0_yz_xz_zz, t_0_yz_yy_xx,\
                           t_0_yz_yy_xy, t_0_yz_yy_xz, t_0_yz_yy_yy, t_0_yz_yy_yz, t_0_yz_yy_zz,\
                           t_0_yz_yz_xx, t_0_yz_yz_xy, t_0_yz_yz_xz, t_0_yz_yz_yy, t_0_yz_yz_yz,\
                           t_0_yz_yz_zz, t_0_yz_zz_xx, t_0_yz_zz_xy, t_0_yz_zz_xz, t_0_yz_zz_yy,\
                           t_0_yz_zz_yz, t_0_yz_zz_zz, t_0_z_xx_xx, t_0_z_xx_xy, t_0_z_xx_xz,\
                           t_0_z_xx_yy, t_0_z_xx_yz, t_0_z_xx_zz, t_0_z_xy_xx, t_0_z_xy_xy,\
                           t_0_z_xy_xz, t_0_z_xy_yy, t_0_z_xy_yz, t_0_z_xy_zz, t_0_z_xz_xx,\
                           t_0_z_xz_xy, t_0_z_xz_xz, t_0_z_xz_yy, t_0_z_xz_yz, t_0_z_xz_zz,\
                           t_0_z_yy_xx, t_0_z_yy_xy, t_0_z_yy_xz, t_0_z_yy_yy, t_0_z_yy_yz,\
                           t_0_z_yy_zz, t_0_z_yz_xx, t_0_z_yz_xy, t_0_z_yz_xz, t_0_z_yz_yy,\
                           t_0_z_yz_yz, t_0_z_yz_zz, t_0_z_zz_xx, t_0_z_zz_xy, t_0_z_zz_xz,\
                           t_0_z_zz_yy, t_0_z_zz_yz, t_0_z_zz_zz, t_y_z_xx_xx, t_y_z_xx_xy,\
                           t_y_z_xx_xz, t_y_z_xx_yy, t_y_z_xx_yz, t_y_z_xx_zz, t_y_z_xy_xx,\
                           t_y_z_xy_xy, t_y_z_xy_xz, t_y_z_xy_yy, t_y_z_xy_yz, t_y_z_xy_zz,\
                           t_y_z_xz_xx, t_y_z_xz_xy, t_y_z_xz_xz, t_y_z_xz_yy, t_y_z_xz_yz,\
                           t_y_z_xz_zz, t_y_z_yy_xx, t_y_z_yy_xy, t_y_z_yy_xz, t_y_z_yy_yy,\
                           t_y_z_yy_yz, t_y_z_yy_zz, t_y_z_yz_xx, t_y_z_yz_xy, t_y_z_yz_xz,\
                           t_y_z_yz_yy, t_y_z_yz_yz, t_y_z_yz_zz, t_y_z_zz_xx, t_y_z_zz_xy,\
                           t_y_z_zz_xz, t_y_z_zz_yy, t_y_z_zz_yz, t_y_z_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_z_zz_zz[i] = t_0_yz_zz_zz[i] - rab_y[i] * t_0_z_zz_zz[i];

        t_y_z_zz_yz[i] = t_0_yz_zz_yz[i] - rab_y[i] * t_0_z_zz_yz[i];

        t_y_z_zz_yy[i] = t_0_yz_zz_yy[i] - rab_y[i] * t_0_z_zz_yy[i];

        t_y_z_zz_xz[i] = t_0_yz_zz_xz[i] - rab_y[i] * t_0_z_zz_xz[i];

        t_y_z_zz_xy[i] = t_0_yz_zz_xy[i] - rab_y[i] * t_0_z_zz_xy[i];

        t_y_z_zz_xx[i] = t_0_yz_zz_xx[i] - rab_y[i] * t_0_z_zz_xx[i];

        t_y_z_yz_zz[i] = t_0_yz_yz_zz[i] - rab_y[i] * t_0_z_yz_zz[i];

        t_y_z_yz_yz[i] = t_0_yz_yz_yz[i] - rab_y[i] * t_0_z_yz_yz[i];

        t_y_z_yz_yy[i] = t_0_yz_yz_yy[i] - rab_y[i] * t_0_z_yz_yy[i];

        t_y_z_yz_xz[i] = t_0_yz_yz_xz[i] - rab_y[i] * t_0_z_yz_xz[i];

        t_y_z_yz_xy[i] = t_0_yz_yz_xy[i] - rab_y[i] * t_0_z_yz_xy[i];

        t_y_z_yz_xx[i] = t_0_yz_yz_xx[i] - rab_y[i] * t_0_z_yz_xx[i];

        t_y_z_yy_zz[i] = t_0_yz_yy_zz[i] - rab_y[i] * t_0_z_yy_zz[i];

        t_y_z_yy_yz[i] = t_0_yz_yy_yz[i] - rab_y[i] * t_0_z_yy_yz[i];

        t_y_z_yy_yy[i] = t_0_yz_yy_yy[i] - rab_y[i] * t_0_z_yy_yy[i];

        t_y_z_yy_xz[i] = t_0_yz_yy_xz[i] - rab_y[i] * t_0_z_yy_xz[i];

        t_y_z_yy_xy[i] = t_0_yz_yy_xy[i] - rab_y[i] * t_0_z_yy_xy[i];

        t_y_z_yy_xx[i] = t_0_yz_yy_xx[i] - rab_y[i] * t_0_z_yy_xx[i];

        t_y_z_xz_zz[i] = t_0_yz_xz_zz[i] - rab_y[i] * t_0_z_xz_zz[i];

        t_y_z_xz_yz[i] = t_0_yz_xz_yz[i] - rab_y[i] * t_0_z_xz_yz[i];

        t_y_z_xz_yy[i] = t_0_yz_xz_yy[i] - rab_y[i] * t_0_z_xz_yy[i];

        t_y_z_xz_xz[i] = t_0_yz_xz_xz[i] - rab_y[i] * t_0_z_xz_xz[i];

        t_y_z_xz_xy[i] = t_0_yz_xz_xy[i] - rab_y[i] * t_0_z_xz_xy[i];

        t_y_z_xz_xx[i] = t_0_yz_xz_xx[i] - rab_y[i] * t_0_z_xz_xx[i];

        t_y_z_xy_zz[i] = t_0_yz_xy_zz[i] - rab_y[i] * t_0_z_xy_zz[i];

        t_y_z_xy_yz[i] = t_0_yz_xy_yz[i] - rab_y[i] * t_0_z_xy_yz[i];

        t_y_z_xy_yy[i] = t_0_yz_xy_yy[i] - rab_y[i] * t_0_z_xy_yy[i];

        t_y_z_xy_xz[i] = t_0_yz_xy_xz[i] - rab_y[i] * t_0_z_xy_xz[i];

        t_y_z_xy_xy[i] = t_0_yz_xy_xy[i] - rab_y[i] * t_0_z_xy_xy[i];

        t_y_z_xy_xx[i] = t_0_yz_xy_xx[i] - rab_y[i] * t_0_z_xy_xx[i];

        t_y_z_xx_zz[i] = t_0_yz_xx_zz[i] - rab_y[i] * t_0_z_xx_zz[i];

        t_y_z_xx_yz[i] = t_0_yz_xx_yz[i] - rab_y[i] * t_0_z_xx_yz[i];

        t_y_z_xx_yy[i] = t_0_yz_xx_yy[i] - rab_y[i] * t_0_z_xx_yy[i];

        t_y_z_xx_xz[i] = t_0_yz_xx_xz[i] - rab_y[i] * t_0_z_xx_xz[i];

        t_y_z_xx_xy[i] = t_0_yz_xx_xy[i] - rab_y[i] * t_0_z_xx_xy[i];

        t_y_z_xx_xx[i] = t_0_yz_xx_xx[i] - rab_y[i] * t_0_z_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_y_xx_xx, t_0_y_xx_xy, t_0_y_xx_xz, t_0_y_xx_yy,\
                           t_0_y_xx_yz, t_0_y_xx_zz, t_0_y_xy_xx, t_0_y_xy_xy, t_0_y_xy_xz,\
                           t_0_y_xy_yy, t_0_y_xy_yz, t_0_y_xy_zz, t_0_y_xz_xx, t_0_y_xz_xy,\
                           t_0_y_xz_xz, t_0_y_xz_yy, t_0_y_xz_yz, t_0_y_xz_zz, t_0_y_yy_xx,\
                           t_0_y_yy_xy, t_0_y_yy_xz, t_0_y_yy_yy, t_0_y_yy_yz, t_0_y_yy_zz,\
                           t_0_y_yz_xx, t_0_y_yz_xy, t_0_y_yz_xz, t_0_y_yz_yy, t_0_y_yz_yz,\
                           t_0_y_yz_zz, t_0_y_zz_xx, t_0_y_zz_xy, t_0_y_zz_xz, t_0_y_zz_yy,\
                           t_0_y_zz_yz, t_0_y_zz_zz, t_0_yy_xx_xx, t_0_yy_xx_xy, t_0_yy_xx_xz,\
                           t_0_yy_xx_yy, t_0_yy_xx_yz, t_0_yy_xx_zz, t_0_yy_xy_xx, t_0_yy_xy_xy,\
                           t_0_yy_xy_xz, t_0_yy_xy_yy, t_0_yy_xy_yz, t_0_yy_xy_zz, t_0_yy_xz_xx,\
                           t_0_yy_xz_xy, t_0_yy_xz_xz, t_0_yy_xz_yy, t_0_yy_xz_yz, t_0_yy_xz_zz,\
                           t_0_yy_yy_xx, t_0_yy_yy_xy, t_0_yy_yy_xz, t_0_yy_yy_yy, t_0_yy_yy_yz,\
                           t_0_yy_yy_zz, t_0_yy_yz_xx, t_0_yy_yz_xy, t_0_yy_yz_xz, t_0_yy_yz_yy,\
                           t_0_yy_yz_yz, t_0_yy_yz_zz, t_0_yy_zz_xx, t_0_yy_zz_xy, t_0_yy_zz_xz,\
                           t_0_yy_zz_yy, t_0_yy_zz_yz, t_0_yy_zz_zz, t_y_y_xx_xx, t_y_y_xx_xy,\
                           t_y_y_xx_xz, t_y_y_xx_yy, t_y_y_xx_yz, t_y_y_xx_zz, t_y_y_xy_xx,\
                           t_y_y_xy_xy, t_y_y_xy_xz, t_y_y_xy_yy, t_y_y_xy_yz, t_y_y_xy_zz,\
                           t_y_y_xz_xx, t_y_y_xz_xy, t_y_y_xz_xz, t_y_y_xz_yy, t_y_y_xz_yz,\
                           t_y_y_xz_zz, t_y_y_yy_xx, t_y_y_yy_xy, t_y_y_yy_xz, t_y_y_yy_yy,\
                           t_y_y_yy_yz, t_y_y_yy_zz, t_y_y_yz_xx, t_y_y_yz_xy, t_y_y_yz_xz,\
                           t_y_y_yz_yy, t_y_y_yz_yz, t_y_y_yz_zz, t_y_y_zz_xx, t_y_y_zz_xy,\
                           t_y_y_zz_xz, t_y_y_zz_yy, t_y_y_zz_yz, t_y_y_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_y_zz_zz[i] = t_0_yy_zz_zz[i] - rab_y[i] * t_0_y_zz_zz[i];

        t_y_y_zz_yz[i] = t_0_yy_zz_yz[i] - rab_y[i] * t_0_y_zz_yz[i];

        t_y_y_zz_yy[i] = t_0_yy_zz_yy[i] - rab_y[i] * t_0_y_zz_yy[i];

        t_y_y_zz_xz[i] = t_0_yy_zz_xz[i] - rab_y[i] * t_0_y_zz_xz[i];

        t_y_y_zz_xy[i] = t_0_yy_zz_xy[i] - rab_y[i] * t_0_y_zz_xy[i];

        t_y_y_zz_xx[i] = t_0_yy_zz_xx[i] - rab_y[i] * t_0_y_zz_xx[i];

        t_y_y_yz_zz[i] = t_0_yy_yz_zz[i] - rab_y[i] * t_0_y_yz_zz[i];

        t_y_y_yz_yz[i] = t_0_yy_yz_yz[i] - rab_y[i] * t_0_y_yz_yz[i];

        t_y_y_yz_yy[i] = t_0_yy_yz_yy[i] - rab_y[i] * t_0_y_yz_yy[i];

        t_y_y_yz_xz[i] = t_0_yy_yz_xz[i] - rab_y[i] * t_0_y_yz_xz[i];

        t_y_y_yz_xy[i] = t_0_yy_yz_xy[i] - rab_y[i] * t_0_y_yz_xy[i];

        t_y_y_yz_xx[i] = t_0_yy_yz_xx[i] - rab_y[i] * t_0_y_yz_xx[i];

        t_y_y_yy_zz[i] = t_0_yy_yy_zz[i] - rab_y[i] * t_0_y_yy_zz[i];

        t_y_y_yy_yz[i] = t_0_yy_yy_yz[i] - rab_y[i] * t_0_y_yy_yz[i];

        t_y_y_yy_yy[i] = t_0_yy_yy_yy[i] - rab_y[i] * t_0_y_yy_yy[i];

        t_y_y_yy_xz[i] = t_0_yy_yy_xz[i] - rab_y[i] * t_0_y_yy_xz[i];

        t_y_y_yy_xy[i] = t_0_yy_yy_xy[i] - rab_y[i] * t_0_y_yy_xy[i];

        t_y_y_yy_xx[i] = t_0_yy_yy_xx[i] - rab_y[i] * t_0_y_yy_xx[i];

        t_y_y_xz_zz[i] = t_0_yy_xz_zz[i] - rab_y[i] * t_0_y_xz_zz[i];

        t_y_y_xz_yz[i] = t_0_yy_xz_yz[i] - rab_y[i] * t_0_y_xz_yz[i];

        t_y_y_xz_yy[i] = t_0_yy_xz_yy[i] - rab_y[i] * t_0_y_xz_yy[i];

        t_y_y_xz_xz[i] = t_0_yy_xz_xz[i] - rab_y[i] * t_0_y_xz_xz[i];

        t_y_y_xz_xy[i] = t_0_yy_xz_xy[i] - rab_y[i] * t_0_y_xz_xy[i];

        t_y_y_xz_xx[i] = t_0_yy_xz_xx[i] - rab_y[i] * t_0_y_xz_xx[i];

        t_y_y_xy_zz[i] = t_0_yy_xy_zz[i] - rab_y[i] * t_0_y_xy_zz[i];

        t_y_y_xy_yz[i] = t_0_yy_xy_yz[i] - rab_y[i] * t_0_y_xy_yz[i];

        t_y_y_xy_yy[i] = t_0_yy_xy_yy[i] - rab_y[i] * t_0_y_xy_yy[i];

        t_y_y_xy_xz[i] = t_0_yy_xy_xz[i] - rab_y[i] * t_0_y_xy_xz[i];

        t_y_y_xy_xy[i] = t_0_yy_xy_xy[i] - rab_y[i] * t_0_y_xy_xy[i];

        t_y_y_xy_xx[i] = t_0_yy_xy_xx[i] - rab_y[i] * t_0_y_xy_xx[i];

        t_y_y_xx_zz[i] = t_0_yy_xx_zz[i] - rab_y[i] * t_0_y_xx_zz[i];

        t_y_y_xx_yz[i] = t_0_yy_xx_yz[i] - rab_y[i] * t_0_y_xx_yz[i];

        t_y_y_xx_yy[i] = t_0_yy_xx_yy[i] - rab_y[i] * t_0_y_xx_yy[i];

        t_y_y_xx_xz[i] = t_0_yy_xx_xz[i] - rab_y[i] * t_0_y_xx_xz[i];

        t_y_y_xx_xy[i] = t_0_yy_xx_xy[i] - rab_y[i] * t_0_y_xx_xy[i];

        t_y_y_xx_xx[i] = t_0_yy_xx_xx[i] - rab_y[i] * t_0_y_xx_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_x_xx_xx, t_0_x_xx_xy, t_0_x_xx_xz, t_0_x_xx_yy,\
                           t_0_x_xx_yz, t_0_x_xx_zz, t_0_x_xy_xx, t_0_x_xy_xy, t_0_x_xy_xz,\
                           t_0_x_xy_yy, t_0_x_xy_yz, t_0_x_xy_zz, t_0_x_xz_xx, t_0_x_xz_xy,\
                           t_0_x_xz_xz, t_0_x_xz_yy, t_0_x_xz_yz, t_0_x_xz_zz, t_0_x_yy_xx,\
                           t_0_x_yy_xy, t_0_x_yy_xz, t_0_x_yy_yy, t_0_x_yy_yz, t_0_x_yy_zz,\
                           t_0_x_yz_xx, t_0_x_yz_xy, t_0_x_yz_xz, t_0_x_yz_yy, t_0_x_yz_yz,\
                           t_0_x_yz_zz, t_0_x_zz_xx, t_0_x_zz_xy, t_0_x_zz_xz, t_0_x_zz_yy,\
                           t_0_x_zz_yz, t_0_x_zz_zz, t_0_xy_xx_xx, t_0_xy_xx_xy, t_0_xy_xx_xz,\
                           t_0_xy_xx_yy, t_0_xy_xx_yz, t_0_xy_xx_zz, t_0_xy_xy_xx, t_0_xy_xy_xy,\
                           t_0_xy_xy_xz, t_0_xy_xy_yy, t_0_xy_xy_yz, t_0_xy_xy_zz, t_0_xy_xz_xx,\
                           t_0_xy_xz_xy, t_0_xy_xz_xz, t_0_xy_xz_yy, t_0_xy_xz_yz, t_0_xy_xz_zz,\
                           t_0_xy_yy_xx, t_0_xy_yy_xy, t_0_xy_yy_xz, t_0_xy_yy_yy, t_0_xy_yy_yz,\
                           t_0_xy_yy_zz, t_0_xy_yz_xx, t_0_xy_yz_xy, t_0_xy_yz_xz, t_0_xy_yz_yy,\
                           t_0_xy_yz_yz, t_0_xy_yz_zz, t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz,\
                           t_0_xy_zz_yy, t_0_xy_zz_yz, t_0_xy_zz_zz, t_y_x_xx_xx, t_y_x_xx_xy,\
                           t_y_x_xx_xz, t_y_x_xx_yy, t_y_x_xx_yz, t_y_x_xx_zz, t_y_x_xy_xx,\
                           t_y_x_xy_xy, t_y_x_xy_xz, t_y_x_xy_yy, t_y_x_xy_yz, t_y_x_xy_zz,\
                           t_y_x_xz_xx, t_y_x_xz_xy, t_y_x_xz_xz, t_y_x_xz_yy, t_y_x_xz_yz,\
                           t_y_x_xz_zz, t_y_x_yy_xx, t_y_x_yy_xy, t_y_x_yy_xz, t_y_x_yy_yy,\
                           t_y_x_yy_yz, t_y_x_yy_zz, t_y_x_yz_xx, t_y_x_yz_xy, t_y_x_yz_xz,\
                           t_y_x_yz_yy, t_y_x_yz_yz, t_y_x_yz_zz, t_y_x_zz_xx, t_y_x_zz_xy,\
                           t_y_x_zz_xz, t_y_x_zz_yy, t_y_x_zz_yz, t_y_x_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_x_zz_zz[i] = t_0_xy_zz_zz[i] - rab_y[i] * t_0_x_zz_zz[i];

        t_y_x_zz_yz[i] = t_0_xy_zz_yz[i] - rab_y[i] * t_0_x_zz_yz[i];

        t_y_x_zz_yy[i] = t_0_xy_zz_yy[i] - rab_y[i] * t_0_x_zz_yy[i];

        t_y_x_zz_xz[i] = t_0_xy_zz_xz[i] - rab_y[i] * t_0_x_zz_xz[i];

        t_y_x_zz_xy[i] = t_0_xy_zz_xy[i] - rab_y[i] * t_0_x_zz_xy[i];

        t_y_x_zz_xx[i] = t_0_xy_zz_xx[i] - rab_y[i] * t_0_x_zz_xx[i];

        t_y_x_yz_zz[i] = t_0_xy_yz_zz[i] - rab_y[i] * t_0_x_yz_zz[i];

        t_y_x_yz_yz[i] = t_0_xy_yz_yz[i] - rab_y[i] * t_0_x_yz_yz[i];

        t_y_x_yz_yy[i] = t_0_xy_yz_yy[i] - rab_y[i] * t_0_x_yz_yy[i];

        t_y_x_yz_xz[i] = t_0_xy_yz_xz[i] - rab_y[i] * t_0_x_yz_xz[i];

        t_y_x_yz_xy[i] = t_0_xy_yz_xy[i] - rab_y[i] * t_0_x_yz_xy[i];

        t_y_x_yz_xx[i] = t_0_xy_yz_xx[i] - rab_y[i] * t_0_x_yz_xx[i];

        t_y_x_yy_zz[i] = t_0_xy_yy_zz[i] - rab_y[i] * t_0_x_yy_zz[i];

        t_y_x_yy_yz[i] = t_0_xy_yy_yz[i] - rab_y[i] * t_0_x_yy_yz[i];

        t_y_x_yy_yy[i] = t_0_xy_yy_yy[i] - rab_y[i] * t_0_x_yy_yy[i];

        t_y_x_yy_xz[i] = t_0_xy_yy_xz[i] - rab_y[i] * t_0_x_yy_xz[i];

        t_y_x_yy_xy[i] = t_0_xy_yy_xy[i] - rab_y[i] * t_0_x_yy_xy[i];

        t_y_x_yy_xx[i] = t_0_xy_yy_xx[i] - rab_y[i] * t_0_x_yy_xx[i];

        t_y_x_xz_zz[i] = t_0_xy_xz_zz[i] - rab_y[i] * t_0_x_xz_zz[i];

        t_y_x_xz_yz[i] = t_0_xy_xz_yz[i] - rab_y[i] * t_0_x_xz_yz[i];

        t_y_x_xz_yy[i] = t_0_xy_xz_yy[i] - rab_y[i] * t_0_x_xz_yy[i];

        t_y_x_xz_xz[i] = t_0_xy_xz_xz[i] - rab_y[i] * t_0_x_xz_xz[i];

        t_y_x_xz_xy[i] = t_0_xy_xz_xy[i] - rab_y[i] * t_0_x_xz_xy[i];

        t_y_x_xz_xx[i] = t_0_xy_xz_xx[i] - rab_y[i] * t_0_x_xz_xx[i];

        t_y_x_xy_zz[i] = t_0_xy_xy_zz[i] - rab_y[i] * t_0_x_xy_zz[i];

        t_y_x_xy_yz[i] = t_0_xy_xy_yz[i] - rab_y[i] * t_0_x_xy_yz[i];

        t_y_x_xy_yy[i] = t_0_xy_xy_yy[i] - rab_y[i] * t_0_x_xy_yy[i];

        t_y_x_xy_xz[i] = t_0_xy_xy_xz[i] - rab_y[i] * t_0_x_xy_xz[i];

        t_y_x_xy_xy[i] = t_0_xy_xy_xy[i] - rab_y[i] * t_0_x_xy_xy[i];

        t_y_x_xy_xx[i] = t_0_xy_xy_xx[i] - rab_y[i] * t_0_x_xy_xx[i];

        t_y_x_xx_zz[i] = t_0_xy_xx_zz[i] - rab_y[i] * t_0_x_xx_zz[i];

        t_y_x_xx_yz[i] = t_0_xy_xx_yz[i] - rab_y[i] * t_0_x_xx_yz[i];

        t_y_x_xx_yy[i] = t_0_xy_xx_yy[i] - rab_y[i] * t_0_x_xx_yy[i];

        t_y_x_xx_xz[i] = t_0_xy_xx_xz[i] - rab_y[i] * t_0_x_xx_xz[i];

        t_y_x_xx_xy[i] = t_0_xy_xx_xy[i] - rab_y[i] * t_0_x_xx_xy[i];

        t_y_x_xx_xx[i] = t_0_xy_xx_xx[i] - rab_y[i] * t_0_x_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xz_xx_xx, t_0_xz_xx_xy, t_0_xz_xx_xz, t_0_xz_xx_yy,\
                           t_0_xz_xx_yz, t_0_xz_xx_zz, t_0_xz_xy_xx, t_0_xz_xy_xy, t_0_xz_xy_xz,\
                           t_0_xz_xy_yy, t_0_xz_xy_yz, t_0_xz_xy_zz, t_0_xz_xz_xx, t_0_xz_xz_xy,\
                           t_0_xz_xz_xz, t_0_xz_xz_yy, t_0_xz_xz_yz, t_0_xz_xz_zz, t_0_xz_yy_xx,\
                           t_0_xz_yy_xy, t_0_xz_yy_xz, t_0_xz_yy_yy, t_0_xz_yy_yz, t_0_xz_yy_zz,\
                           t_0_xz_yz_xx, t_0_xz_yz_xy, t_0_xz_yz_xz, t_0_xz_yz_yy, t_0_xz_yz_yz,\
                           t_0_xz_yz_zz, t_0_xz_zz_xx, t_0_xz_zz_xy, t_0_xz_zz_xz, t_0_xz_zz_yy,\
                           t_0_xz_zz_yz, t_0_xz_zz_zz, t_0_z_xx_xx, t_0_z_xx_xy, t_0_z_xx_xz,\
                           t_0_z_xx_yy, t_0_z_xx_yz, t_0_z_xx_zz, t_0_z_xy_xx, t_0_z_xy_xy,\
                           t_0_z_xy_xz, t_0_z_xy_yy, t_0_z_xy_yz, t_0_z_xy_zz, t_0_z_xz_xx,\
                           t_0_z_xz_xy, t_0_z_xz_xz, t_0_z_xz_yy, t_0_z_xz_yz, t_0_z_xz_zz,\
                           t_0_z_yy_xx, t_0_z_yy_xy, t_0_z_yy_xz, t_0_z_yy_yy, t_0_z_yy_yz,\
                           t_0_z_yy_zz, t_0_z_yz_xx, t_0_z_yz_xy, t_0_z_yz_xz, t_0_z_yz_yy,\
                           t_0_z_yz_yz, t_0_z_yz_zz, t_0_z_zz_xx, t_0_z_zz_xy, t_0_z_zz_xz,\
                           t_0_z_zz_yy, t_0_z_zz_yz, t_0_z_zz_zz, t_x_z_xx_xx, t_x_z_xx_xy,\
                           t_x_z_xx_xz, t_x_z_xx_yy, t_x_z_xx_yz, t_x_z_xx_zz, t_x_z_xy_xx,\
                           t_x_z_xy_xy, t_x_z_xy_xz, t_x_z_xy_yy, t_x_z_xy_yz, t_x_z_xy_zz,\
                           t_x_z_xz_xx, t_x_z_xz_xy, t_x_z_xz_xz, t_x_z_xz_yy, t_x_z_xz_yz,\
                           t_x_z_xz_zz, t_x_z_yy_xx, t_x_z_yy_xy, t_x_z_yy_xz, t_x_z_yy_yy,\
                           t_x_z_yy_yz, t_x_z_yy_zz, t_x_z_yz_xx, t_x_z_yz_xy, t_x_z_yz_xz,\
                           t_x_z_yz_yy, t_x_z_yz_yz, t_x_z_yz_zz, t_x_z_zz_xx, t_x_z_zz_xy,\
                           t_x_z_zz_xz, t_x_z_zz_yy, t_x_z_zz_yz, t_x_z_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_z_zz_zz[i] = t_0_xz_zz_zz[i] - rab_x[i] * t_0_z_zz_zz[i];

        t_x_z_zz_yz[i] = t_0_xz_zz_yz[i] - rab_x[i] * t_0_z_zz_yz[i];

        t_x_z_zz_yy[i] = t_0_xz_zz_yy[i] - rab_x[i] * t_0_z_zz_yy[i];

        t_x_z_zz_xz[i] = t_0_xz_zz_xz[i] - rab_x[i] * t_0_z_zz_xz[i];

        t_x_z_zz_xy[i] = t_0_xz_zz_xy[i] - rab_x[i] * t_0_z_zz_xy[i];

        t_x_z_zz_xx[i] = t_0_xz_zz_xx[i] - rab_x[i] * t_0_z_zz_xx[i];

        t_x_z_yz_zz[i] = t_0_xz_yz_zz[i] - rab_x[i] * t_0_z_yz_zz[i];

        t_x_z_yz_yz[i] = t_0_xz_yz_yz[i] - rab_x[i] * t_0_z_yz_yz[i];

        t_x_z_yz_yy[i] = t_0_xz_yz_yy[i] - rab_x[i] * t_0_z_yz_yy[i];

        t_x_z_yz_xz[i] = t_0_xz_yz_xz[i] - rab_x[i] * t_0_z_yz_xz[i];

        t_x_z_yz_xy[i] = t_0_xz_yz_xy[i] - rab_x[i] * t_0_z_yz_xy[i];

        t_x_z_yz_xx[i] = t_0_xz_yz_xx[i] - rab_x[i] * t_0_z_yz_xx[i];

        t_x_z_yy_zz[i] = t_0_xz_yy_zz[i] - rab_x[i] * t_0_z_yy_zz[i];

        t_x_z_yy_yz[i] = t_0_xz_yy_yz[i] - rab_x[i] * t_0_z_yy_yz[i];

        t_x_z_yy_yy[i] = t_0_xz_yy_yy[i] - rab_x[i] * t_0_z_yy_yy[i];

        t_x_z_yy_xz[i] = t_0_xz_yy_xz[i] - rab_x[i] * t_0_z_yy_xz[i];

        t_x_z_yy_xy[i] = t_0_xz_yy_xy[i] - rab_x[i] * t_0_z_yy_xy[i];

        t_x_z_yy_xx[i] = t_0_xz_yy_xx[i] - rab_x[i] * t_0_z_yy_xx[i];

        t_x_z_xz_zz[i] = t_0_xz_xz_zz[i] - rab_x[i] * t_0_z_xz_zz[i];

        t_x_z_xz_yz[i] = t_0_xz_xz_yz[i] - rab_x[i] * t_0_z_xz_yz[i];

        t_x_z_xz_yy[i] = t_0_xz_xz_yy[i] - rab_x[i] * t_0_z_xz_yy[i];

        t_x_z_xz_xz[i] = t_0_xz_xz_xz[i] - rab_x[i] * t_0_z_xz_xz[i];

        t_x_z_xz_xy[i] = t_0_xz_xz_xy[i] - rab_x[i] * t_0_z_xz_xy[i];

        t_x_z_xz_xx[i] = t_0_xz_xz_xx[i] - rab_x[i] * t_0_z_xz_xx[i];

        t_x_z_xy_zz[i] = t_0_xz_xy_zz[i] - rab_x[i] * t_0_z_xy_zz[i];

        t_x_z_xy_yz[i] = t_0_xz_xy_yz[i] - rab_x[i] * t_0_z_xy_yz[i];

        t_x_z_xy_yy[i] = t_0_xz_xy_yy[i] - rab_x[i] * t_0_z_xy_yy[i];

        t_x_z_xy_xz[i] = t_0_xz_xy_xz[i] - rab_x[i] * t_0_z_xy_xz[i];

        t_x_z_xy_xy[i] = t_0_xz_xy_xy[i] - rab_x[i] * t_0_z_xy_xy[i];

        t_x_z_xy_xx[i] = t_0_xz_xy_xx[i] - rab_x[i] * t_0_z_xy_xx[i];

        t_x_z_xx_zz[i] = t_0_xz_xx_zz[i] - rab_x[i] * t_0_z_xx_zz[i];

        t_x_z_xx_yz[i] = t_0_xz_xx_yz[i] - rab_x[i] * t_0_z_xx_yz[i];

        t_x_z_xx_yy[i] = t_0_xz_xx_yy[i] - rab_x[i] * t_0_z_xx_yy[i];

        t_x_z_xx_xz[i] = t_0_xz_xx_xz[i] - rab_x[i] * t_0_z_xx_xz[i];

        t_x_z_xx_xy[i] = t_0_xz_xx_xy[i] - rab_x[i] * t_0_z_xx_xy[i];

        t_x_z_xx_xx[i] = t_0_xz_xx_xx[i] - rab_x[i] * t_0_z_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xy_xx_xx, t_0_xy_xx_xy, t_0_xy_xx_xz, t_0_xy_xx_yy,\
                           t_0_xy_xx_yz, t_0_xy_xx_zz, t_0_xy_xy_xx, t_0_xy_xy_xy, t_0_xy_xy_xz,\
                           t_0_xy_xy_yy, t_0_xy_xy_yz, t_0_xy_xy_zz, t_0_xy_xz_xx, t_0_xy_xz_xy,\
                           t_0_xy_xz_xz, t_0_xy_xz_yy, t_0_xy_xz_yz, t_0_xy_xz_zz, t_0_xy_yy_xx,\
                           t_0_xy_yy_xy, t_0_xy_yy_xz, t_0_xy_yy_yy, t_0_xy_yy_yz, t_0_xy_yy_zz,\
                           t_0_xy_yz_xx, t_0_xy_yz_xy, t_0_xy_yz_xz, t_0_xy_yz_yy, t_0_xy_yz_yz,\
                           t_0_xy_yz_zz, t_0_xy_zz_xx, t_0_xy_zz_xy, t_0_xy_zz_xz, t_0_xy_zz_yy,\
                           t_0_xy_zz_yz, t_0_xy_zz_zz, t_0_y_xx_xx, t_0_y_xx_xy, t_0_y_xx_xz,\
                           t_0_y_xx_yy, t_0_y_xx_yz, t_0_y_xx_zz, t_0_y_xy_xx, t_0_y_xy_xy,\
                           t_0_y_xy_xz, t_0_y_xy_yy, t_0_y_xy_yz, t_0_y_xy_zz, t_0_y_xz_xx,\
                           t_0_y_xz_xy, t_0_y_xz_xz, t_0_y_xz_yy, t_0_y_xz_yz, t_0_y_xz_zz,\
                           t_0_y_yy_xx, t_0_y_yy_xy, t_0_y_yy_xz, t_0_y_yy_yy, t_0_y_yy_yz,\
                           t_0_y_yy_zz, t_0_y_yz_xx, t_0_y_yz_xy, t_0_y_yz_xz, t_0_y_yz_yy,\
                           t_0_y_yz_yz, t_0_y_yz_zz, t_0_y_zz_xx, t_0_y_zz_xy, t_0_y_zz_xz,\
                           t_0_y_zz_yy, t_0_y_zz_yz, t_0_y_zz_zz, t_x_y_xx_xx, t_x_y_xx_xy,\
                           t_x_y_xx_xz, t_x_y_xx_yy, t_x_y_xx_yz, t_x_y_xx_zz, t_x_y_xy_xx,\
                           t_x_y_xy_xy, t_x_y_xy_xz, t_x_y_xy_yy, t_x_y_xy_yz, t_x_y_xy_zz,\
                           t_x_y_xz_xx, t_x_y_xz_xy, t_x_y_xz_xz, t_x_y_xz_yy, t_x_y_xz_yz,\
                           t_x_y_xz_zz, t_x_y_yy_xx, t_x_y_yy_xy, t_x_y_yy_xz, t_x_y_yy_yy,\
                           t_x_y_yy_yz, t_x_y_yy_zz, t_x_y_yz_xx, t_x_y_yz_xy, t_x_y_yz_xz,\
                           t_x_y_yz_yy, t_x_y_yz_yz, t_x_y_yz_zz, t_x_y_zz_xx, t_x_y_zz_xy,\
                           t_x_y_zz_xz, t_x_y_zz_yy, t_x_y_zz_yz, t_x_y_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_y_zz_zz[i] = t_0_xy_zz_zz[i] - rab_x[i] * t_0_y_zz_zz[i];

        t_x_y_zz_yz[i] = t_0_xy_zz_yz[i] - rab_x[i] * t_0_y_zz_yz[i];

        t_x_y_zz_yy[i] = t_0_xy_zz_yy[i] - rab_x[i] * t_0_y_zz_yy[i];

        t_x_y_zz_xz[i] = t_0_xy_zz_xz[i] - rab_x[i] * t_0_y_zz_xz[i];

        t_x_y_zz_xy[i] = t_0_xy_zz_xy[i] - rab_x[i] * t_0_y_zz_xy[i];

        t_x_y_zz_xx[i] = t_0_xy_zz_xx[i] - rab_x[i] * t_0_y_zz_xx[i];

        t_x_y_yz_zz[i] = t_0_xy_yz_zz[i] - rab_x[i] * t_0_y_yz_zz[i];

        t_x_y_yz_yz[i] = t_0_xy_yz_yz[i] - rab_x[i] * t_0_y_yz_yz[i];

        t_x_y_yz_yy[i] = t_0_xy_yz_yy[i] - rab_x[i] * t_0_y_yz_yy[i];

        t_x_y_yz_xz[i] = t_0_xy_yz_xz[i] - rab_x[i] * t_0_y_yz_xz[i];

        t_x_y_yz_xy[i] = t_0_xy_yz_xy[i] - rab_x[i] * t_0_y_yz_xy[i];

        t_x_y_yz_xx[i] = t_0_xy_yz_xx[i] - rab_x[i] * t_0_y_yz_xx[i];

        t_x_y_yy_zz[i] = t_0_xy_yy_zz[i] - rab_x[i] * t_0_y_yy_zz[i];

        t_x_y_yy_yz[i] = t_0_xy_yy_yz[i] - rab_x[i] * t_0_y_yy_yz[i];

        t_x_y_yy_yy[i] = t_0_xy_yy_yy[i] - rab_x[i] * t_0_y_yy_yy[i];

        t_x_y_yy_xz[i] = t_0_xy_yy_xz[i] - rab_x[i] * t_0_y_yy_xz[i];

        t_x_y_yy_xy[i] = t_0_xy_yy_xy[i] - rab_x[i] * t_0_y_yy_xy[i];

        t_x_y_yy_xx[i] = t_0_xy_yy_xx[i] - rab_x[i] * t_0_y_yy_xx[i];

        t_x_y_xz_zz[i] = t_0_xy_xz_zz[i] - rab_x[i] * t_0_y_xz_zz[i];

        t_x_y_xz_yz[i] = t_0_xy_xz_yz[i] - rab_x[i] * t_0_y_xz_yz[i];

        t_x_y_xz_yy[i] = t_0_xy_xz_yy[i] - rab_x[i] * t_0_y_xz_yy[i];

        t_x_y_xz_xz[i] = t_0_xy_xz_xz[i] - rab_x[i] * t_0_y_xz_xz[i];

        t_x_y_xz_xy[i] = t_0_xy_xz_xy[i] - rab_x[i] * t_0_y_xz_xy[i];

        t_x_y_xz_xx[i] = t_0_xy_xz_xx[i] - rab_x[i] * t_0_y_xz_xx[i];

        t_x_y_xy_zz[i] = t_0_xy_xy_zz[i] - rab_x[i] * t_0_y_xy_zz[i];

        t_x_y_xy_yz[i] = t_0_xy_xy_yz[i] - rab_x[i] * t_0_y_xy_yz[i];

        t_x_y_xy_yy[i] = t_0_xy_xy_yy[i] - rab_x[i] * t_0_y_xy_yy[i];

        t_x_y_xy_xz[i] = t_0_xy_xy_xz[i] - rab_x[i] * t_0_y_xy_xz[i];

        t_x_y_xy_xy[i] = t_0_xy_xy_xy[i] - rab_x[i] * t_0_y_xy_xy[i];

        t_x_y_xy_xx[i] = t_0_xy_xy_xx[i] - rab_x[i] * t_0_y_xy_xx[i];

        t_x_y_xx_zz[i] = t_0_xy_xx_zz[i] - rab_x[i] * t_0_y_xx_zz[i];

        t_x_y_xx_yz[i] = t_0_xy_xx_yz[i] - rab_x[i] * t_0_y_xx_yz[i];

        t_x_y_xx_yy[i] = t_0_xy_xx_yy[i] - rab_x[i] * t_0_y_xx_yy[i];

        t_x_y_xx_xz[i] = t_0_xy_xx_xz[i] - rab_x[i] * t_0_y_xx_xz[i];

        t_x_y_xx_xy[i] = t_0_xy_xx_xy[i] - rab_x[i] * t_0_y_xx_xy[i];

        t_x_y_xx_xx[i] = t_0_xy_xx_xx[i] - rab_x[i] * t_0_y_xx_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_x_xx_xx, t_0_x_xx_xy, t_0_x_xx_xz, t_0_x_xx_yy,\
                           t_0_x_xx_yz, t_0_x_xx_zz, t_0_x_xy_xx, t_0_x_xy_xy, t_0_x_xy_xz,\
                           t_0_x_xy_yy, t_0_x_xy_yz, t_0_x_xy_zz, t_0_x_xz_xx, t_0_x_xz_xy,\
                           t_0_x_xz_xz, t_0_x_xz_yy, t_0_x_xz_yz, t_0_x_xz_zz, t_0_x_yy_xx,\
                           t_0_x_yy_xy, t_0_x_yy_xz, t_0_x_yy_yy, t_0_x_yy_yz, t_0_x_yy_zz,\
                           t_0_x_yz_xx, t_0_x_yz_xy, t_0_x_yz_xz, t_0_x_yz_yy, t_0_x_yz_yz,\
                           t_0_x_yz_zz, t_0_x_zz_xx, t_0_x_zz_xy, t_0_x_zz_xz, t_0_x_zz_yy,\
                           t_0_x_zz_yz, t_0_x_zz_zz, t_0_xx_xx_xx, t_0_xx_xx_xy, t_0_xx_xx_xz,\
                           t_0_xx_xx_yy, t_0_xx_xx_yz, t_0_xx_xx_zz, t_0_xx_xy_xx, t_0_xx_xy_xy,\
                           t_0_xx_xy_xz, t_0_xx_xy_yy, t_0_xx_xy_yz, t_0_xx_xy_zz, t_0_xx_xz_xx,\
                           t_0_xx_xz_xy, t_0_xx_xz_xz, t_0_xx_xz_yy, t_0_xx_xz_yz, t_0_xx_xz_zz,\
                           t_0_xx_yy_xx, t_0_xx_yy_xy, t_0_xx_yy_xz, t_0_xx_yy_yy, t_0_xx_yy_yz,\
                           t_0_xx_yy_zz, t_0_xx_yz_xx, t_0_xx_yz_xy, t_0_xx_yz_xz, t_0_xx_yz_yy,\
                           t_0_xx_yz_yz, t_0_xx_yz_zz, t_0_xx_zz_xx, t_0_xx_zz_xy, t_0_xx_zz_xz,\
                           t_0_xx_zz_yy, t_0_xx_zz_yz, t_0_xx_zz_zz, t_x_x_xx_xx, t_x_x_xx_xy,\
                           t_x_x_xx_xz, t_x_x_xx_yy, t_x_x_xx_yz, t_x_x_xx_zz, t_x_x_xy_xx,\
                           t_x_x_xy_xy, t_x_x_xy_xz, t_x_x_xy_yy, t_x_x_xy_yz, t_x_x_xy_zz,\
                           t_x_x_xz_xx, t_x_x_xz_xy, t_x_x_xz_xz, t_x_x_xz_yy, t_x_x_xz_yz,\
                           t_x_x_xz_zz, t_x_x_yy_xx, t_x_x_yy_xy, t_x_x_yy_xz, t_x_x_yy_yy,\
                           t_x_x_yy_yz, t_x_x_yy_zz, t_x_x_yz_xx, t_x_x_yz_xy, t_x_x_yz_xz,\
                           t_x_x_yz_yy, t_x_x_yz_yz, t_x_x_yz_zz, t_x_x_zz_xx, t_x_x_zz_xy,\
                           t_x_x_zz_xz, t_x_x_zz_yy, t_x_x_zz_yz, t_x_x_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_x_zz_zz[i] = t_0_xx_zz_zz[i] - rab_x[i] * t_0_x_zz_zz[i];

        t_x_x_zz_yz[i] = t_0_xx_zz_yz[i] - rab_x[i] * t_0_x_zz_yz[i];

        t_x_x_zz_yy[i] = t_0_xx_zz_yy[i] - rab_x[i] * t_0_x_zz_yy[i];

        t_x_x_zz_xz[i] = t_0_xx_zz_xz[i] - rab_x[i] * t_0_x_zz_xz[i];

        t_x_x_zz_xy[i] = t_0_xx_zz_xy[i] - rab_x[i] * t_0_x_zz_xy[i];

        t_x_x_zz_xx[i] = t_0_xx_zz_xx[i] - rab_x[i] * t_0_x_zz_xx[i];

        t_x_x_yz_zz[i] = t_0_xx_yz_zz[i] - rab_x[i] * t_0_x_yz_zz[i];

        t_x_x_yz_yz[i] = t_0_xx_yz_yz[i] - rab_x[i] * t_0_x_yz_yz[i];

        t_x_x_yz_yy[i] = t_0_xx_yz_yy[i] - rab_x[i] * t_0_x_yz_yy[i];

        t_x_x_yz_xz[i] = t_0_xx_yz_xz[i] - rab_x[i] * t_0_x_yz_xz[i];

        t_x_x_yz_xy[i] = t_0_xx_yz_xy[i] - rab_x[i] * t_0_x_yz_xy[i];

        t_x_x_yz_xx[i] = t_0_xx_yz_xx[i] - rab_x[i] * t_0_x_yz_xx[i];

        t_x_x_yy_zz[i] = t_0_xx_yy_zz[i] - rab_x[i] * t_0_x_yy_zz[i];

        t_x_x_yy_yz[i] = t_0_xx_yy_yz[i] - rab_x[i] * t_0_x_yy_yz[i];

        t_x_x_yy_yy[i] = t_0_xx_yy_yy[i] - rab_x[i] * t_0_x_yy_yy[i];

        t_x_x_yy_xz[i] = t_0_xx_yy_xz[i] - rab_x[i] * t_0_x_yy_xz[i];

        t_x_x_yy_xy[i] = t_0_xx_yy_xy[i] - rab_x[i] * t_0_x_yy_xy[i];

        t_x_x_yy_xx[i] = t_0_xx_yy_xx[i] - rab_x[i] * t_0_x_yy_xx[i];

        t_x_x_xz_zz[i] = t_0_xx_xz_zz[i] - rab_x[i] * t_0_x_xz_zz[i];

        t_x_x_xz_yz[i] = t_0_xx_xz_yz[i] - rab_x[i] * t_0_x_xz_yz[i];

        t_x_x_xz_yy[i] = t_0_xx_xz_yy[i] - rab_x[i] * t_0_x_xz_yy[i];

        t_x_x_xz_xz[i] = t_0_xx_xz_xz[i] - rab_x[i] * t_0_x_xz_xz[i];

        t_x_x_xz_xy[i] = t_0_xx_xz_xy[i] - rab_x[i] * t_0_x_xz_xy[i];

        t_x_x_xz_xx[i] = t_0_xx_xz_xx[i] - rab_x[i] * t_0_x_xz_xx[i];

        t_x_x_xy_zz[i] = t_0_xx_xy_zz[i] - rab_x[i] * t_0_x_xy_zz[i];

        t_x_x_xy_yz[i] = t_0_xx_xy_yz[i] - rab_x[i] * t_0_x_xy_yz[i];

        t_x_x_xy_yy[i] = t_0_xx_xy_yy[i] - rab_x[i] * t_0_x_xy_yy[i];

        t_x_x_xy_xz[i] = t_0_xx_xy_xz[i] - rab_x[i] * t_0_x_xy_xz[i];

        t_x_x_xy_xy[i] = t_0_xx_xy_xy[i] - rab_x[i] * t_0_x_xy_xy[i];

        t_x_x_xy_xx[i] = t_0_xx_xy_xx[i] - rab_x[i] * t_0_x_xy_xx[i];

        t_x_x_xx_zz[i] = t_0_xx_xx_zz[i] - rab_x[i] * t_0_x_xx_zz[i];

        t_x_x_xx_yz[i] = t_0_xx_xx_yz[i] - rab_x[i] * t_0_x_xx_yz[i];

        t_x_x_xx_yy[i] = t_0_xx_xx_yy[i] - rab_x[i] * t_0_x_xx_yy[i];

        t_x_x_xx_xz[i] = t_0_xx_xx_xz[i] - rab_x[i] * t_0_x_xx_xz[i];

        t_x_x_xx_xy[i] = t_0_xx_xx_xy[i] - rab_x[i] * t_0_x_xx_xy[i];

        t_x_x_xx_xx[i] = t_0_xx_xx_xx[i] - rab_x[i] * t_0_x_xx_xx[i];
    }
}


} // derirec namespace
