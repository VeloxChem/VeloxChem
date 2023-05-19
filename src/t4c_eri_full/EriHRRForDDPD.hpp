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
compHostHRRForDDPD_V0(      BufferHostXY<T>&      intsBufferDDPD,
                      const BufferHostX<int32_t>& intsIndexesDDPD,
                      const BufferHostXY<T>&      intsBufferPDPD,
                      const BufferHostX<int32_t>& intsIndexesPDPD,
                      const BufferHostXY<T>&      intsBufferPFPD,
                      const BufferHostX<int32_t>& intsIndexesPFPD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDPD) integral components

    t_zz_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(0));

    t_zz_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(1));

    t_zz_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(2));

    t_zz_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(3));

    t_zz_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(4));

    t_zz_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(5));

    t_zz_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(6));

    t_zz_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(7));

    t_zz_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(8));

    t_zz_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(9));

    t_zz_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(10));

    t_zz_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(11));

    t_zz_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(12));

    t_zz_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(13));

    t_zz_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(14));

    t_zz_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(15));

    t_zz_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(16));

    t_zz_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(17));

    t_zz_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(18));

    t_zz_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(19));

    t_zz_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(20));

    t_zz_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(21));

    t_zz_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(22));

    t_zz_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(23));

    t_zz_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(24));

    t_zz_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(25));

    t_zz_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(26));

    t_zz_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(27));

    t_zz_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(28));

    t_zz_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(29));

    t_zz_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(30));

    t_zz_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(31));

    t_zz_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(32));

    t_zz_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(33));

    t_zz_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(34));

    t_zz_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(35));

    t_zz_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(36));

    t_zz_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(37));

    t_zz_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(38));

    t_zz_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(39));

    t_zz_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(40));

    t_zz_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(41));

    t_zz_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(42));

    t_zz_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(43));

    t_zz_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(44));

    t_zz_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(45));

    t_zz_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(46));

    t_zz_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(47));

    t_zz_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(48));

    t_zz_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(49));

    t_zz_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(50));

    t_zz_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(51));

    t_zz_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(52));

    t_zz_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(53));

    t_zz_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(54));

    t_zz_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(55));

    t_zz_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(56));

    t_zz_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(57));

    t_zz_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(58));

    t_zz_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(59));

    t_zz_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(60));

    t_zz_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(61));

    t_zz_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(62));

    t_zz_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(63));

    t_zz_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(64));

    t_zz_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(65));

    t_zz_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(66));

    t_zz_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(67));

    t_zz_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(68));

    t_zz_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(69));

    t_zz_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(70));

    t_zz_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(71));

    t_zz_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(72));

    t_zz_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(73));

    t_zz_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(74));

    t_zz_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(75));

    t_zz_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(76));

    t_zz_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(77));

    t_zz_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(78));

    t_zz_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(79));

    t_zz_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(80));

    t_zz_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(81));

    t_zz_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(82));

    t_zz_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(83));

    t_zz_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(84));

    t_zz_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(85));

    t_zz_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(86));

    t_zz_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(87));

    t_zz_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(88));

    t_zz_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(89));

    t_zz_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(90));

    t_zz_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(91));

    t_zz_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(92));

    t_zz_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(93));

    t_zz_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(94));

    t_zz_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(95));

    t_zz_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(96));

    t_zz_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(97));

    t_zz_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(98));

    t_zz_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(99));

    t_zz_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(100));

    t_zz_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(101));

    t_zz_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(102));

    t_zz_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(103));

    t_zz_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(104));

    t_zz_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(105));

    t_zz_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(106));

    t_zz_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(107));

    t_yz_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(108));

    t_yz_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(109));

    t_yz_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(110));

    t_yz_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(111));

    t_yz_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(112));

    t_yz_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(113));

    t_yz_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(114));

    t_yz_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(115));

    t_yz_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(116));

    t_yz_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(117));

    t_yz_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(118));

    t_yz_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(119));

    t_yz_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(120));

    t_yz_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(121));

    t_yz_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(122));

    t_yz_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(123));

    t_yz_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(124));

    t_yz_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(125));

    t_yz_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(126));

    t_yz_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(127));

    t_yz_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(128));

    t_yz_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(129));

    t_yz_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(130));

    t_yz_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(131));

    t_yz_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(132));

    t_yz_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(133));

    t_yz_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(134));

    t_yz_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(135));

    t_yz_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(136));

    t_yz_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(137));

    t_yz_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(138));

    t_yz_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(139));

    t_yz_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(140));

    t_yz_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(141));

    t_yz_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(142));

    t_yz_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(143));

    t_yz_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(144));

    t_yz_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(145));

    t_yz_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(146));

    t_yz_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(147));

    t_yz_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(148));

    t_yz_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(149));

    t_yz_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(150));

    t_yz_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(151));

    t_yz_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(152));

    t_yz_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(153));

    t_yz_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(154));

    t_yz_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(155));

    t_yz_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(156));

    t_yz_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(157));

    t_yz_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(158));

    t_yz_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(159));

    t_yz_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(160));

    t_yz_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(161));

    t_yz_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(162));

    t_yz_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(163));

    t_yz_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(164));

    t_yz_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(165));

    t_yz_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(166));

    t_yz_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(167));

    t_yz_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(168));

    t_yz_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(169));

    t_yz_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(170));

    t_yz_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(171));

    t_yz_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(172));

    t_yz_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(173));

    t_yz_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(174));

    t_yz_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(175));

    t_yz_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(176));

    t_yz_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(177));

    t_yz_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(178));

    t_yz_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(179));

    t_yz_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(180));

    t_yz_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(181));

    t_yz_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(182));

    t_yz_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(183));

    t_yz_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(184));

    t_yz_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(185));

    t_yz_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(186));

    t_yz_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(187));

    t_yz_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(188));

    t_yz_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(189));

    t_yz_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(190));

    t_yz_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(191));

    t_yz_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(192));

    t_yz_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(193));

    t_yz_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(194));

    t_yz_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(195));

    t_yz_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(196));

    t_yz_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(197));

    t_yz_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(198));

    t_yz_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(199));

    t_yz_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(200));

    t_yz_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(201));

    t_yz_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(202));

    t_yz_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(203));

    t_yz_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(204));

    t_yz_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(205));

    t_yz_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(206));

    t_yz_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(207));

    t_yz_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(208));

    t_yz_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(209));

    t_yz_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(210));

    t_yz_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(211));

    t_yz_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(212));

    t_yz_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(213));

    t_yz_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(214));

    t_yz_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(215));

    t_yy_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(216));

    t_yy_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(217));

    t_yy_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(218));

    t_yy_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(219));

    t_yy_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(220));

    t_yy_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(221));

    t_yy_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(222));

    t_yy_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(223));

    t_yy_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(224));

    t_yy_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(225));

    t_yy_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(226));

    t_yy_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(227));

    t_yy_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(228));

    t_yy_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(229));

    t_yy_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(230));

    t_yy_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(231));

    t_yy_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(232));

    t_yy_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(233));

    t_yy_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(234));

    t_yy_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(235));

    t_yy_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(236));

    t_yy_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(237));

    t_yy_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(238));

    t_yy_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(239));

    t_yy_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(240));

    t_yy_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(241));

    t_yy_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(242));

    t_yy_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(243));

    t_yy_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(244));

    t_yy_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(245));

    t_yy_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(246));

    t_yy_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(247));

    t_yy_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(248));

    t_yy_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(249));

    t_yy_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(250));

    t_yy_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(251));

    t_yy_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(252));

    t_yy_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(253));

    t_yy_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(254));

    t_yy_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(255));

    t_yy_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(256));

    t_yy_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(257));

    t_yy_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(258));

    t_yy_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(259));

    t_yy_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(260));

    t_yy_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(261));

    t_yy_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(262));

    t_yy_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(263));

    t_yy_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(264));

    t_yy_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(265));

    t_yy_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(266));

    t_yy_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(267));

    t_yy_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(268));

    t_yy_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(269));

    t_yy_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(270));

    t_yy_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(271));

    t_yy_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(272));

    t_yy_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(273));

    t_yy_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(274));

    t_yy_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(275));

    t_yy_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(276));

    t_yy_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(277));

    t_yy_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(278));

    t_yy_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(279));

    t_yy_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(280));

    t_yy_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(281));

    t_yy_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(282));

    t_yy_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(283));

    t_yy_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(284));

    t_yy_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(285));

    t_yy_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(286));

    t_yy_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(287));

    t_yy_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(288));

    t_yy_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(289));

    t_yy_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(290));

    t_yy_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(291));

    t_yy_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(292));

    t_yy_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(293));

    t_yy_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(294));

    t_yy_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(295));

    t_yy_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(296));

    t_yy_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(297));

    t_yy_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(298));

    t_yy_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(299));

    t_yy_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(300));

    t_yy_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(301));

    t_yy_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(302));

    t_yy_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(303));

    t_yy_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(304));

    t_yy_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(305));

    t_yy_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(306));

    t_yy_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(307));

    t_yy_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(308));

    t_yy_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(309));

    t_yy_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(310));

    t_yy_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(311));

    t_yy_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(312));

    t_yy_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(313));

    t_yy_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(314));

    t_yy_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(315));

    t_yy_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(316));

    t_yy_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(317));

    t_yy_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(318));

    t_yy_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(319));

    t_yy_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(320));

    t_yy_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(321));

    t_yy_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(322));

    t_yy_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(323));

    t_xz_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(324));

    t_xz_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(325));

    t_xz_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(326));

    t_xz_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(327));

    t_xz_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(328));

    t_xz_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(329));

    t_xz_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(330));

    t_xz_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(331));

    t_xz_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(332));

    t_xz_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(333));

    t_xz_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(334));

    t_xz_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(335));

    t_xz_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(336));

    t_xz_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(337));

    t_xz_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(338));

    t_xz_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(339));

    t_xz_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(340));

    t_xz_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(341));

    t_xz_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(342));

    t_xz_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(343));

    t_xz_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(344));

    t_xz_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(345));

    t_xz_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(346));

    t_xz_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(347));

    t_xz_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(348));

    t_xz_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(349));

    t_xz_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(350));

    t_xz_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(351));

    t_xz_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(352));

    t_xz_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(353));

    t_xz_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(354));

    t_xz_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(355));

    t_xz_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(356));

    t_xz_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(357));

    t_xz_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(358));

    t_xz_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(359));

    t_xz_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(360));

    t_xz_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(361));

    t_xz_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(362));

    t_xz_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(363));

    t_xz_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(364));

    t_xz_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(365));

    t_xz_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(366));

    t_xz_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(367));

    t_xz_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(368));

    t_xz_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(369));

    t_xz_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(370));

    t_xz_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(371));

    t_xz_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(372));

    t_xz_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(373));

    t_xz_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(374));

    t_xz_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(375));

    t_xz_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(376));

    t_xz_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(377));

    t_xz_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(378));

    t_xz_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(379));

    t_xz_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(380));

    t_xz_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(381));

    t_xz_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(382));

    t_xz_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(383));

    t_xz_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(384));

    t_xz_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(385));

    t_xz_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(386));

    t_xz_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(387));

    t_xz_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(388));

    t_xz_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(389));

    t_xz_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(390));

    t_xz_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(391));

    t_xz_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(392));

    t_xz_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(393));

    t_xz_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(394));

    t_xz_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(395));

    t_xz_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(396));

    t_xz_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(397));

    t_xz_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(398));

    t_xz_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(399));

    t_xz_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(400));

    t_xz_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(401));

    t_xz_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(402));

    t_xz_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(403));

    t_xz_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(404));

    t_xz_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(405));

    t_xz_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(406));

    t_xz_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(407));

    t_xz_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(408));

    t_xz_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(409));

    t_xz_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(410));

    t_xz_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(411));

    t_xz_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(412));

    t_xz_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(413));

    t_xz_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(414));

    t_xz_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(415));

    t_xz_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(416));

    t_xz_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(417));

    t_xz_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(418));

    t_xz_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(419));

    t_xz_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(420));

    t_xz_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(421));

    t_xz_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(422));

    t_xz_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(423));

    t_xz_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(424));

    t_xz_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(425));

    t_xz_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(426));

    t_xz_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(427));

    t_xz_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(428));

    t_xz_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(429));

    t_xz_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(430));

    t_xz_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(431));

    t_xy_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(432));

    t_xy_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(433));

    t_xy_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(434));

    t_xy_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(435));

    t_xy_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(436));

    t_xy_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(437));

    t_xy_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(438));

    t_xy_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(439));

    t_xy_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(440));

    t_xy_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(441));

    t_xy_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(442));

    t_xy_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(443));

    t_xy_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(444));

    t_xy_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(445));

    t_xy_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(446));

    t_xy_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(447));

    t_xy_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(448));

    t_xy_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(449));

    t_xy_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(450));

    t_xy_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(451));

    t_xy_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(452));

    t_xy_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(453));

    t_xy_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(454));

    t_xy_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(455));

    t_xy_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(456));

    t_xy_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(457));

    t_xy_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(458));

    t_xy_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(459));

    t_xy_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(460));

    t_xy_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(461));

    t_xy_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(462));

    t_xy_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(463));

    t_xy_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(464));

    t_xy_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(465));

    t_xy_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(466));

    t_xy_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(467));

    t_xy_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(468));

    t_xy_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(469));

    t_xy_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(470));

    t_xy_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(471));

    t_xy_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(472));

    t_xy_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(473));

    t_xy_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(474));

    t_xy_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(475));

    t_xy_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(476));

    t_xy_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(477));

    t_xy_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(478));

    t_xy_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(479));

    t_xy_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(480));

    t_xy_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(481));

    t_xy_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(482));

    t_xy_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(483));

    t_xy_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(484));

    t_xy_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(485));

    t_xy_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(486));

    t_xy_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(487));

    t_xy_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(488));

    t_xy_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(489));

    t_xy_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(490));

    t_xy_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(491));

    t_xy_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(492));

    t_xy_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(493));

    t_xy_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(494));

    t_xy_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(495));

    t_xy_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(496));

    t_xy_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(497));

    t_xy_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(498));

    t_xy_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(499));

    t_xy_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(500));

    t_xy_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(501));

    t_xy_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(502));

    t_xy_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(503));

    t_xy_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(504));

    t_xy_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(505));

    t_xy_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(506));

    t_xy_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(507));

    t_xy_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(508));

    t_xy_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(509));

    t_xy_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(510));

    t_xy_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(511));

    t_xy_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(512));

    t_xy_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(513));

    t_xy_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(514));

    t_xy_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(515));

    t_xy_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(516));

    t_xy_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(517));

    t_xy_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(518));

    t_xy_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(519));

    t_xy_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(520));

    t_xy_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(521));

    t_xy_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(522));

    t_xy_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(523));

    t_xy_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(524));

    t_xy_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(525));

    t_xy_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(526));

    t_xy_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(527));

    t_xy_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(528));

    t_xy_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(529));

    t_xy_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(530));

    t_xy_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(531));

    t_xy_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(532));

    t_xy_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(533));

    t_xy_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(534));

    t_xy_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(535));

    t_xy_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(536));

    t_xy_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(537));

    t_xy_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(538));

    t_xy_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(539));

    t_xx_zz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(540));

    t_xx_zz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(541));

    t_xx_zz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(542));

    t_xx_zz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(543));

    t_xx_zz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(544));

    t_xx_zz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(545));

    t_xx_zz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(546));

    t_xx_zz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(547));

    t_xx_zz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(548));

    t_xx_zz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(549));

    t_xx_zz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(550));

    t_xx_zz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(551));

    t_xx_zz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(552));

    t_xx_zz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(553));

    t_xx_zz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(554));

    t_xx_zz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(555));

    t_xx_zz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(556));

    t_xx_zz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(557));

    t_xx_yz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(558));

    t_xx_yz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(559));

    t_xx_yz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(560));

    t_xx_yz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(561));

    t_xx_yz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(562));

    t_xx_yz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(563));

    t_xx_yz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(564));

    t_xx_yz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(565));

    t_xx_yz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(566));

    t_xx_yz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(567));

    t_xx_yz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(568));

    t_xx_yz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(569));

    t_xx_yz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(570));

    t_xx_yz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(571));

    t_xx_yz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(572));

    t_xx_yz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(573));

    t_xx_yz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(574));

    t_xx_yz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(575));

    t_xx_yy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(576));

    t_xx_yy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(577));

    t_xx_yy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(578));

    t_xx_yy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(579));

    t_xx_yy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(580));

    t_xx_yy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(581));

    t_xx_yy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(582));

    t_xx_yy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(583));

    t_xx_yy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(584));

    t_xx_yy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(585));

    t_xx_yy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(586));

    t_xx_yy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(587));

    t_xx_yy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(588));

    t_xx_yy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(589));

    t_xx_yy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(590));

    t_xx_yy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(591));

    t_xx_yy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(592));

    t_xx_yy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(593));

    t_xx_xz_z_zz = intsBufferDDPD.data(intsIndexesDDPD(594));

    t_xx_xz_z_yz = intsBufferDDPD.data(intsIndexesDDPD(595));

    t_xx_xz_z_yy = intsBufferDDPD.data(intsIndexesDDPD(596));

    t_xx_xz_z_xz = intsBufferDDPD.data(intsIndexesDDPD(597));

    t_xx_xz_z_xy = intsBufferDDPD.data(intsIndexesDDPD(598));

    t_xx_xz_z_xx = intsBufferDDPD.data(intsIndexesDDPD(599));

    t_xx_xz_y_zz = intsBufferDDPD.data(intsIndexesDDPD(600));

    t_xx_xz_y_yz = intsBufferDDPD.data(intsIndexesDDPD(601));

    t_xx_xz_y_yy = intsBufferDDPD.data(intsIndexesDDPD(602));

    t_xx_xz_y_xz = intsBufferDDPD.data(intsIndexesDDPD(603));

    t_xx_xz_y_xy = intsBufferDDPD.data(intsIndexesDDPD(604));

    t_xx_xz_y_xx = intsBufferDDPD.data(intsIndexesDDPD(605));

    t_xx_xz_x_zz = intsBufferDDPD.data(intsIndexesDDPD(606));

    t_xx_xz_x_yz = intsBufferDDPD.data(intsIndexesDDPD(607));

    t_xx_xz_x_yy = intsBufferDDPD.data(intsIndexesDDPD(608));

    t_xx_xz_x_xz = intsBufferDDPD.data(intsIndexesDDPD(609));

    t_xx_xz_x_xy = intsBufferDDPD.data(intsIndexesDDPD(610));

    t_xx_xz_x_xx = intsBufferDDPD.data(intsIndexesDDPD(611));

    t_xx_xy_z_zz = intsBufferDDPD.data(intsIndexesDDPD(612));

    t_xx_xy_z_yz = intsBufferDDPD.data(intsIndexesDDPD(613));

    t_xx_xy_z_yy = intsBufferDDPD.data(intsIndexesDDPD(614));

    t_xx_xy_z_xz = intsBufferDDPD.data(intsIndexesDDPD(615));

    t_xx_xy_z_xy = intsBufferDDPD.data(intsIndexesDDPD(616));

    t_xx_xy_z_xx = intsBufferDDPD.data(intsIndexesDDPD(617));

    t_xx_xy_y_zz = intsBufferDDPD.data(intsIndexesDDPD(618));

    t_xx_xy_y_yz = intsBufferDDPD.data(intsIndexesDDPD(619));

    t_xx_xy_y_yy = intsBufferDDPD.data(intsIndexesDDPD(620));

    t_xx_xy_y_xz = intsBufferDDPD.data(intsIndexesDDPD(621));

    t_xx_xy_y_xy = intsBufferDDPD.data(intsIndexesDDPD(622));

    t_xx_xy_y_xx = intsBufferDDPD.data(intsIndexesDDPD(623));

    t_xx_xy_x_zz = intsBufferDDPD.data(intsIndexesDDPD(624));

    t_xx_xy_x_yz = intsBufferDDPD.data(intsIndexesDDPD(625));

    t_xx_xy_x_yy = intsBufferDDPD.data(intsIndexesDDPD(626));

    t_xx_xy_x_xz = intsBufferDDPD.data(intsIndexesDDPD(627));

    t_xx_xy_x_xy = intsBufferDDPD.data(intsIndexesDDPD(628));

    t_xx_xy_x_xx = intsBufferDDPD.data(intsIndexesDDPD(629));

    t_xx_xx_z_zz = intsBufferDDPD.data(intsIndexesDDPD(630));

    t_xx_xx_z_yz = intsBufferDDPD.data(intsIndexesDDPD(631));

    t_xx_xx_z_yy = intsBufferDDPD.data(intsIndexesDDPD(632));

    t_xx_xx_z_xz = intsBufferDDPD.data(intsIndexesDDPD(633));

    t_xx_xx_z_xy = intsBufferDDPD.data(intsIndexesDDPD(634));

    t_xx_xx_z_xx = intsBufferDDPD.data(intsIndexesDDPD(635));

    t_xx_xx_y_zz = intsBufferDDPD.data(intsIndexesDDPD(636));

    t_xx_xx_y_yz = intsBufferDDPD.data(intsIndexesDDPD(637));

    t_xx_xx_y_yy = intsBufferDDPD.data(intsIndexesDDPD(638));

    t_xx_xx_y_xz = intsBufferDDPD.data(intsIndexesDDPD(639));

    t_xx_xx_y_xy = intsBufferDDPD.data(intsIndexesDDPD(640));

    t_xx_xx_y_xx = intsBufferDDPD.data(intsIndexesDDPD(641));

    t_xx_xx_x_zz = intsBufferDDPD.data(intsIndexesDDPD(642));

    t_xx_xx_x_yz = intsBufferDDPD.data(intsIndexesDDPD(643));

    t_xx_xx_x_yy = intsBufferDDPD.data(intsIndexesDDPD(644));

    t_xx_xx_x_xz = intsBufferDDPD.data(intsIndexesDDPD(645));

    t_xx_xx_x_xy = intsBufferDDPD.data(intsIndexesDDPD(646));

    t_xx_xx_x_xx = intsBufferDDPD.data(intsIndexesDDPD(647));

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

    #pragma omp simd align(rab_z, t_z_yz_x_xx, t_z_yz_x_xy, t_z_yz_x_xz, t_z_yz_x_yy,\
                           t_z_yz_x_yz, t_z_yz_x_zz, t_z_yz_y_xx, t_z_yz_y_xy, t_z_yz_y_xz,\
                           t_z_yz_y_yy, t_z_yz_y_yz, t_z_yz_y_zz, t_z_yz_z_xx, t_z_yz_z_xy,\
                           t_z_yz_z_xz, t_z_yz_z_yy, t_z_yz_z_yz, t_z_yz_z_zz, t_z_yzz_x_xx,\
                           t_z_yzz_x_xy, t_z_yzz_x_xz, t_z_yzz_x_yy, t_z_yzz_x_yz, t_z_yzz_x_zz,\
                           t_z_yzz_y_xx, t_z_yzz_y_xy, t_z_yzz_y_xz, t_z_yzz_y_yy, t_z_yzz_y_yz,\
                           t_z_yzz_y_zz, t_z_yzz_z_xx, t_z_yzz_z_xy, t_z_yzz_z_xz, t_z_yzz_z_yy,\
                           t_z_yzz_z_yz, t_z_yzz_z_zz, t_z_zz_x_xx, t_z_zz_x_xy, t_z_zz_x_xz,\
                           t_z_zz_x_yy, t_z_zz_x_yz, t_z_zz_x_zz, t_z_zz_y_xx, t_z_zz_y_xy,\
                           t_z_zz_y_xz, t_z_zz_y_yy, t_z_zz_y_yz, t_z_zz_y_zz, t_z_zz_z_xx,\
                           t_z_zz_z_xy, t_z_zz_z_xz, t_z_zz_z_yy, t_z_zz_z_yz, t_z_zz_z_zz,\
                           t_z_zzz_x_xx, t_z_zzz_x_xy, t_z_zzz_x_xz, t_z_zzz_x_yy, t_z_zzz_x_yz,\
                           t_z_zzz_x_zz, t_z_zzz_y_xx, t_z_zzz_y_xy, t_z_zzz_y_xz, t_z_zzz_y_yy,\
                           t_z_zzz_y_yz, t_z_zzz_y_zz, t_z_zzz_z_xx, t_z_zzz_z_xy, t_z_zzz_z_xz,\
                           t_z_zzz_z_yy, t_z_zzz_z_yz, t_z_zzz_z_zz, t_zz_yz_x_xx, t_zz_yz_x_xy,\
                           t_zz_yz_x_xz, t_zz_yz_x_yy, t_zz_yz_x_yz, t_zz_yz_x_zz, t_zz_yz_y_xx,\
                           t_zz_yz_y_xy, t_zz_yz_y_xz, t_zz_yz_y_yy, t_zz_yz_y_yz, t_zz_yz_y_zz,\
                           t_zz_yz_z_xx, t_zz_yz_z_xy, t_zz_yz_z_xz, t_zz_yz_z_yy, t_zz_yz_z_yz,\
                           t_zz_yz_z_zz, t_zz_zz_x_xx, t_zz_zz_x_xy, t_zz_zz_x_xz, t_zz_zz_x_yy,\
                           t_zz_zz_x_yz, t_zz_zz_x_zz, t_zz_zz_y_xx, t_zz_zz_y_xy, t_zz_zz_y_xz,\
                           t_zz_zz_y_yy, t_zz_zz_y_yz, t_zz_zz_y_zz, t_zz_zz_z_xx, t_zz_zz_z_xy,\
                           t_zz_zz_z_xz, t_zz_zz_z_yy, t_zz_zz_z_yz, t_zz_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_z_zz[i] = t_z_zzz_z_zz[i] - rab_z[i] * t_z_zz_z_zz[i];

        t_zz_zz_z_yz[i] = t_z_zzz_z_yz[i] - rab_z[i] * t_z_zz_z_yz[i];

        t_zz_zz_z_yy[i] = t_z_zzz_z_yy[i] - rab_z[i] * t_z_zz_z_yy[i];

        t_zz_zz_z_xz[i] = t_z_zzz_z_xz[i] - rab_z[i] * t_z_zz_z_xz[i];

        t_zz_zz_z_xy[i] = t_z_zzz_z_xy[i] - rab_z[i] * t_z_zz_z_xy[i];

        t_zz_zz_z_xx[i] = t_z_zzz_z_xx[i] - rab_z[i] * t_z_zz_z_xx[i];

        t_zz_zz_y_zz[i] = t_z_zzz_y_zz[i] - rab_z[i] * t_z_zz_y_zz[i];

        t_zz_zz_y_yz[i] = t_z_zzz_y_yz[i] - rab_z[i] * t_z_zz_y_yz[i];

        t_zz_zz_y_yy[i] = t_z_zzz_y_yy[i] - rab_z[i] * t_z_zz_y_yy[i];

        t_zz_zz_y_xz[i] = t_z_zzz_y_xz[i] - rab_z[i] * t_z_zz_y_xz[i];

        t_zz_zz_y_xy[i] = t_z_zzz_y_xy[i] - rab_z[i] * t_z_zz_y_xy[i];

        t_zz_zz_y_xx[i] = t_z_zzz_y_xx[i] - rab_z[i] * t_z_zz_y_xx[i];

        t_zz_zz_x_zz[i] = t_z_zzz_x_zz[i] - rab_z[i] * t_z_zz_x_zz[i];

        t_zz_zz_x_yz[i] = t_z_zzz_x_yz[i] - rab_z[i] * t_z_zz_x_yz[i];

        t_zz_zz_x_yy[i] = t_z_zzz_x_yy[i] - rab_z[i] * t_z_zz_x_yy[i];

        t_zz_zz_x_xz[i] = t_z_zzz_x_xz[i] - rab_z[i] * t_z_zz_x_xz[i];

        t_zz_zz_x_xy[i] = t_z_zzz_x_xy[i] - rab_z[i] * t_z_zz_x_xy[i];

        t_zz_zz_x_xx[i] = t_z_zzz_x_xx[i] - rab_z[i] * t_z_zz_x_xx[i];

        t_zz_yz_z_zz[i] = t_z_yzz_z_zz[i] - rab_z[i] * t_z_yz_z_zz[i];

        t_zz_yz_z_yz[i] = t_z_yzz_z_yz[i] - rab_z[i] * t_z_yz_z_yz[i];

        t_zz_yz_z_yy[i] = t_z_yzz_z_yy[i] - rab_z[i] * t_z_yz_z_yy[i];

        t_zz_yz_z_xz[i] = t_z_yzz_z_xz[i] - rab_z[i] * t_z_yz_z_xz[i];

        t_zz_yz_z_xy[i] = t_z_yzz_z_xy[i] - rab_z[i] * t_z_yz_z_xy[i];

        t_zz_yz_z_xx[i] = t_z_yzz_z_xx[i] - rab_z[i] * t_z_yz_z_xx[i];

        t_zz_yz_y_zz[i] = t_z_yzz_y_zz[i] - rab_z[i] * t_z_yz_y_zz[i];

        t_zz_yz_y_yz[i] = t_z_yzz_y_yz[i] - rab_z[i] * t_z_yz_y_yz[i];

        t_zz_yz_y_yy[i] = t_z_yzz_y_yy[i] - rab_z[i] * t_z_yz_y_yy[i];

        t_zz_yz_y_xz[i] = t_z_yzz_y_xz[i] - rab_z[i] * t_z_yz_y_xz[i];

        t_zz_yz_y_xy[i] = t_z_yzz_y_xy[i] - rab_z[i] * t_z_yz_y_xy[i];

        t_zz_yz_y_xx[i] = t_z_yzz_y_xx[i] - rab_z[i] * t_z_yz_y_xx[i];

        t_zz_yz_x_zz[i] = t_z_yzz_x_zz[i] - rab_z[i] * t_z_yz_x_zz[i];

        t_zz_yz_x_yz[i] = t_z_yzz_x_yz[i] - rab_z[i] * t_z_yz_x_yz[i];

        t_zz_yz_x_yy[i] = t_z_yzz_x_yy[i] - rab_z[i] * t_z_yz_x_yy[i];

        t_zz_yz_x_xz[i] = t_z_yzz_x_xz[i] - rab_z[i] * t_z_yz_x_xz[i];

        t_zz_yz_x_xy[i] = t_z_yzz_x_xy[i] - rab_z[i] * t_z_yz_x_xy[i];

        t_zz_yz_x_xx[i] = t_z_yzz_x_xx[i] - rab_z[i] * t_z_yz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_xz_x_xx, t_z_xz_x_xy, t_z_xz_x_xz, t_z_xz_x_yy,\
                           t_z_xz_x_yz, t_z_xz_x_zz, t_z_xz_y_xx, t_z_xz_y_xy, t_z_xz_y_xz,\
                           t_z_xz_y_yy, t_z_xz_y_yz, t_z_xz_y_zz, t_z_xz_z_xx, t_z_xz_z_xy,\
                           t_z_xz_z_xz, t_z_xz_z_yy, t_z_xz_z_yz, t_z_xz_z_zz, t_z_xzz_x_xx,\
                           t_z_xzz_x_xy, t_z_xzz_x_xz, t_z_xzz_x_yy, t_z_xzz_x_yz, t_z_xzz_x_zz,\
                           t_z_xzz_y_xx, t_z_xzz_y_xy, t_z_xzz_y_xz, t_z_xzz_y_yy, t_z_xzz_y_yz,\
                           t_z_xzz_y_zz, t_z_xzz_z_xx, t_z_xzz_z_xy, t_z_xzz_z_xz, t_z_xzz_z_yy,\
                           t_z_xzz_z_yz, t_z_xzz_z_zz, t_z_yy_x_xx, t_z_yy_x_xy, t_z_yy_x_xz,\
                           t_z_yy_x_yy, t_z_yy_x_yz, t_z_yy_x_zz, t_z_yy_y_xx, t_z_yy_y_xy,\
                           t_z_yy_y_xz, t_z_yy_y_yy, t_z_yy_y_yz, t_z_yy_y_zz, t_z_yy_z_xx,\
                           t_z_yy_z_xy, t_z_yy_z_xz, t_z_yy_z_yy, t_z_yy_z_yz, t_z_yy_z_zz,\
                           t_z_yyz_x_xx, t_z_yyz_x_xy, t_z_yyz_x_xz, t_z_yyz_x_yy, t_z_yyz_x_yz,\
                           t_z_yyz_x_zz, t_z_yyz_y_xx, t_z_yyz_y_xy, t_z_yyz_y_xz, t_z_yyz_y_yy,\
                           t_z_yyz_y_yz, t_z_yyz_y_zz, t_z_yyz_z_xx, t_z_yyz_z_xy, t_z_yyz_z_xz,\
                           t_z_yyz_z_yy, t_z_yyz_z_yz, t_z_yyz_z_zz, t_zz_xz_x_xx, t_zz_xz_x_xy,\
                           t_zz_xz_x_xz, t_zz_xz_x_yy, t_zz_xz_x_yz, t_zz_xz_x_zz, t_zz_xz_y_xx,\
                           t_zz_xz_y_xy, t_zz_xz_y_xz, t_zz_xz_y_yy, t_zz_xz_y_yz, t_zz_xz_y_zz,\
                           t_zz_xz_z_xx, t_zz_xz_z_xy, t_zz_xz_z_xz, t_zz_xz_z_yy, t_zz_xz_z_yz,\
                           t_zz_xz_z_zz, t_zz_yy_x_xx, t_zz_yy_x_xy, t_zz_yy_x_xz, t_zz_yy_x_yy,\
                           t_zz_yy_x_yz, t_zz_yy_x_zz, t_zz_yy_y_xx, t_zz_yy_y_xy, t_zz_yy_y_xz,\
                           t_zz_yy_y_yy, t_zz_yy_y_yz, t_zz_yy_y_zz, t_zz_yy_z_xx, t_zz_yy_z_xy,\
                           t_zz_yy_z_xz, t_zz_yy_z_yy, t_zz_yy_z_yz, t_zz_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_yy_z_zz[i] = t_z_yyz_z_zz[i] - rab_z[i] * t_z_yy_z_zz[i];

        t_zz_yy_z_yz[i] = t_z_yyz_z_yz[i] - rab_z[i] * t_z_yy_z_yz[i];

        t_zz_yy_z_yy[i] = t_z_yyz_z_yy[i] - rab_z[i] * t_z_yy_z_yy[i];

        t_zz_yy_z_xz[i] = t_z_yyz_z_xz[i] - rab_z[i] * t_z_yy_z_xz[i];

        t_zz_yy_z_xy[i] = t_z_yyz_z_xy[i] - rab_z[i] * t_z_yy_z_xy[i];

        t_zz_yy_z_xx[i] = t_z_yyz_z_xx[i] - rab_z[i] * t_z_yy_z_xx[i];

        t_zz_yy_y_zz[i] = t_z_yyz_y_zz[i] - rab_z[i] * t_z_yy_y_zz[i];

        t_zz_yy_y_yz[i] = t_z_yyz_y_yz[i] - rab_z[i] * t_z_yy_y_yz[i];

        t_zz_yy_y_yy[i] = t_z_yyz_y_yy[i] - rab_z[i] * t_z_yy_y_yy[i];

        t_zz_yy_y_xz[i] = t_z_yyz_y_xz[i] - rab_z[i] * t_z_yy_y_xz[i];

        t_zz_yy_y_xy[i] = t_z_yyz_y_xy[i] - rab_z[i] * t_z_yy_y_xy[i];

        t_zz_yy_y_xx[i] = t_z_yyz_y_xx[i] - rab_z[i] * t_z_yy_y_xx[i];

        t_zz_yy_x_zz[i] = t_z_yyz_x_zz[i] - rab_z[i] * t_z_yy_x_zz[i];

        t_zz_yy_x_yz[i] = t_z_yyz_x_yz[i] - rab_z[i] * t_z_yy_x_yz[i];

        t_zz_yy_x_yy[i] = t_z_yyz_x_yy[i] - rab_z[i] * t_z_yy_x_yy[i];

        t_zz_yy_x_xz[i] = t_z_yyz_x_xz[i] - rab_z[i] * t_z_yy_x_xz[i];

        t_zz_yy_x_xy[i] = t_z_yyz_x_xy[i] - rab_z[i] * t_z_yy_x_xy[i];

        t_zz_yy_x_xx[i] = t_z_yyz_x_xx[i] - rab_z[i] * t_z_yy_x_xx[i];

        t_zz_xz_z_zz[i] = t_z_xzz_z_zz[i] - rab_z[i] * t_z_xz_z_zz[i];

        t_zz_xz_z_yz[i] = t_z_xzz_z_yz[i] - rab_z[i] * t_z_xz_z_yz[i];

        t_zz_xz_z_yy[i] = t_z_xzz_z_yy[i] - rab_z[i] * t_z_xz_z_yy[i];

        t_zz_xz_z_xz[i] = t_z_xzz_z_xz[i] - rab_z[i] * t_z_xz_z_xz[i];

        t_zz_xz_z_xy[i] = t_z_xzz_z_xy[i] - rab_z[i] * t_z_xz_z_xy[i];

        t_zz_xz_z_xx[i] = t_z_xzz_z_xx[i] - rab_z[i] * t_z_xz_z_xx[i];

        t_zz_xz_y_zz[i] = t_z_xzz_y_zz[i] - rab_z[i] * t_z_xz_y_zz[i];

        t_zz_xz_y_yz[i] = t_z_xzz_y_yz[i] - rab_z[i] * t_z_xz_y_yz[i];

        t_zz_xz_y_yy[i] = t_z_xzz_y_yy[i] - rab_z[i] * t_z_xz_y_yy[i];

        t_zz_xz_y_xz[i] = t_z_xzz_y_xz[i] - rab_z[i] * t_z_xz_y_xz[i];

        t_zz_xz_y_xy[i] = t_z_xzz_y_xy[i] - rab_z[i] * t_z_xz_y_xy[i];

        t_zz_xz_y_xx[i] = t_z_xzz_y_xx[i] - rab_z[i] * t_z_xz_y_xx[i];

        t_zz_xz_x_zz[i] = t_z_xzz_x_zz[i] - rab_z[i] * t_z_xz_x_zz[i];

        t_zz_xz_x_yz[i] = t_z_xzz_x_yz[i] - rab_z[i] * t_z_xz_x_yz[i];

        t_zz_xz_x_yy[i] = t_z_xzz_x_yy[i] - rab_z[i] * t_z_xz_x_yy[i];

        t_zz_xz_x_xz[i] = t_z_xzz_x_xz[i] - rab_z[i] * t_z_xz_x_xz[i];

        t_zz_xz_x_xy[i] = t_z_xzz_x_xy[i] - rab_z[i] * t_z_xz_x_xy[i];

        t_zz_xz_x_xx[i] = t_z_xzz_x_xx[i] - rab_z[i] * t_z_xz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_z_xx_x_xx, t_z_xx_x_xy, t_z_xx_x_xz, t_z_xx_x_yy,\
                           t_z_xx_x_yz, t_z_xx_x_zz, t_z_xx_y_xx, t_z_xx_y_xy, t_z_xx_y_xz,\
                           t_z_xx_y_yy, t_z_xx_y_yz, t_z_xx_y_zz, t_z_xx_z_xx, t_z_xx_z_xy,\
                           t_z_xx_z_xz, t_z_xx_z_yy, t_z_xx_z_yz, t_z_xx_z_zz, t_z_xxz_x_xx,\
                           t_z_xxz_x_xy, t_z_xxz_x_xz, t_z_xxz_x_yy, t_z_xxz_x_yz, t_z_xxz_x_zz,\
                           t_z_xxz_y_xx, t_z_xxz_y_xy, t_z_xxz_y_xz, t_z_xxz_y_yy, t_z_xxz_y_yz,\
                           t_z_xxz_y_zz, t_z_xxz_z_xx, t_z_xxz_z_xy, t_z_xxz_z_xz, t_z_xxz_z_yy,\
                           t_z_xxz_z_yz, t_z_xxz_z_zz, t_z_xy_x_xx, t_z_xy_x_xy, t_z_xy_x_xz,\
                           t_z_xy_x_yy, t_z_xy_x_yz, t_z_xy_x_zz, t_z_xy_y_xx, t_z_xy_y_xy,\
                           t_z_xy_y_xz, t_z_xy_y_yy, t_z_xy_y_yz, t_z_xy_y_zz, t_z_xy_z_xx,\
                           t_z_xy_z_xy, t_z_xy_z_xz, t_z_xy_z_yy, t_z_xy_z_yz, t_z_xy_z_zz,\
                           t_z_xyz_x_xx, t_z_xyz_x_xy, t_z_xyz_x_xz, t_z_xyz_x_yy, t_z_xyz_x_yz,\
                           t_z_xyz_x_zz, t_z_xyz_y_xx, t_z_xyz_y_xy, t_z_xyz_y_xz, t_z_xyz_y_yy,\
                           t_z_xyz_y_yz, t_z_xyz_y_zz, t_z_xyz_z_xx, t_z_xyz_z_xy, t_z_xyz_z_xz,\
                           t_z_xyz_z_yy, t_z_xyz_z_yz, t_z_xyz_z_zz, t_zz_xx_x_xx, t_zz_xx_x_xy,\
                           t_zz_xx_x_xz, t_zz_xx_x_yy, t_zz_xx_x_yz, t_zz_xx_x_zz, t_zz_xx_y_xx,\
                           t_zz_xx_y_xy, t_zz_xx_y_xz, t_zz_xx_y_yy, t_zz_xx_y_yz, t_zz_xx_y_zz,\
                           t_zz_xx_z_xx, t_zz_xx_z_xy, t_zz_xx_z_xz, t_zz_xx_z_yy, t_zz_xx_z_yz,\
                           t_zz_xx_z_zz, t_zz_xy_x_xx, t_zz_xy_x_xy, t_zz_xy_x_xz, t_zz_xy_x_yy,\
                           t_zz_xy_x_yz, t_zz_xy_x_zz, t_zz_xy_y_xx, t_zz_xy_y_xy, t_zz_xy_y_xz,\
                           t_zz_xy_y_yy, t_zz_xy_y_yz, t_zz_xy_y_zz, t_zz_xy_z_xx, t_zz_xy_z_xy,\
                           t_zz_xy_z_xz, t_zz_xy_z_yy, t_zz_xy_z_yz, t_zz_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_xy_z_zz[i] = t_z_xyz_z_zz[i] - rab_z[i] * t_z_xy_z_zz[i];

        t_zz_xy_z_yz[i] = t_z_xyz_z_yz[i] - rab_z[i] * t_z_xy_z_yz[i];

        t_zz_xy_z_yy[i] = t_z_xyz_z_yy[i] - rab_z[i] * t_z_xy_z_yy[i];

        t_zz_xy_z_xz[i] = t_z_xyz_z_xz[i] - rab_z[i] * t_z_xy_z_xz[i];

        t_zz_xy_z_xy[i] = t_z_xyz_z_xy[i] - rab_z[i] * t_z_xy_z_xy[i];

        t_zz_xy_z_xx[i] = t_z_xyz_z_xx[i] - rab_z[i] * t_z_xy_z_xx[i];

        t_zz_xy_y_zz[i] = t_z_xyz_y_zz[i] - rab_z[i] * t_z_xy_y_zz[i];

        t_zz_xy_y_yz[i] = t_z_xyz_y_yz[i] - rab_z[i] * t_z_xy_y_yz[i];

        t_zz_xy_y_yy[i] = t_z_xyz_y_yy[i] - rab_z[i] * t_z_xy_y_yy[i];

        t_zz_xy_y_xz[i] = t_z_xyz_y_xz[i] - rab_z[i] * t_z_xy_y_xz[i];

        t_zz_xy_y_xy[i] = t_z_xyz_y_xy[i] - rab_z[i] * t_z_xy_y_xy[i];

        t_zz_xy_y_xx[i] = t_z_xyz_y_xx[i] - rab_z[i] * t_z_xy_y_xx[i];

        t_zz_xy_x_zz[i] = t_z_xyz_x_zz[i] - rab_z[i] * t_z_xy_x_zz[i];

        t_zz_xy_x_yz[i] = t_z_xyz_x_yz[i] - rab_z[i] * t_z_xy_x_yz[i];

        t_zz_xy_x_yy[i] = t_z_xyz_x_yy[i] - rab_z[i] * t_z_xy_x_yy[i];

        t_zz_xy_x_xz[i] = t_z_xyz_x_xz[i] - rab_z[i] * t_z_xy_x_xz[i];

        t_zz_xy_x_xy[i] = t_z_xyz_x_xy[i] - rab_z[i] * t_z_xy_x_xy[i];

        t_zz_xy_x_xx[i] = t_z_xyz_x_xx[i] - rab_z[i] * t_z_xy_x_xx[i];

        t_zz_xx_z_zz[i] = t_z_xxz_z_zz[i] - rab_z[i] * t_z_xx_z_zz[i];

        t_zz_xx_z_yz[i] = t_z_xxz_z_yz[i] - rab_z[i] * t_z_xx_z_yz[i];

        t_zz_xx_z_yy[i] = t_z_xxz_z_yy[i] - rab_z[i] * t_z_xx_z_yy[i];

        t_zz_xx_z_xz[i] = t_z_xxz_z_xz[i] - rab_z[i] * t_z_xx_z_xz[i];

        t_zz_xx_z_xy[i] = t_z_xxz_z_xy[i] - rab_z[i] * t_z_xx_z_xy[i];

        t_zz_xx_z_xx[i] = t_z_xxz_z_xx[i] - rab_z[i] * t_z_xx_z_xx[i];

        t_zz_xx_y_zz[i] = t_z_xxz_y_zz[i] - rab_z[i] * t_z_xx_y_zz[i];

        t_zz_xx_y_yz[i] = t_z_xxz_y_yz[i] - rab_z[i] * t_z_xx_y_yz[i];

        t_zz_xx_y_yy[i] = t_z_xxz_y_yy[i] - rab_z[i] * t_z_xx_y_yy[i];

        t_zz_xx_y_xz[i] = t_z_xxz_y_xz[i] - rab_z[i] * t_z_xx_y_xz[i];

        t_zz_xx_y_xy[i] = t_z_xxz_y_xy[i] - rab_z[i] * t_z_xx_y_xy[i];

        t_zz_xx_y_xx[i] = t_z_xxz_y_xx[i] - rab_z[i] * t_z_xx_y_xx[i];

        t_zz_xx_x_zz[i] = t_z_xxz_x_zz[i] - rab_z[i] * t_z_xx_x_zz[i];

        t_zz_xx_x_yz[i] = t_z_xxz_x_yz[i] - rab_z[i] * t_z_xx_x_yz[i];

        t_zz_xx_x_yy[i] = t_z_xxz_x_yy[i] - rab_z[i] * t_z_xx_x_yy[i];

        t_zz_xx_x_xz[i] = t_z_xxz_x_xz[i] - rab_z[i] * t_z_xx_x_xz[i];

        t_zz_xx_x_xy[i] = t_z_xxz_x_xy[i] - rab_z[i] * t_z_xx_x_xy[i];

        t_zz_xx_x_xx[i] = t_z_xxz_x_xx[i] - rab_z[i] * t_z_xx_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_yz_x_xx, t_y_yz_x_xy, t_y_yz_x_xz, t_y_yz_x_yy,\
                           t_y_yz_x_yz, t_y_yz_x_zz, t_y_yz_y_xx, t_y_yz_y_xy, t_y_yz_y_xz,\
                           t_y_yz_y_yy, t_y_yz_y_yz, t_y_yz_y_zz, t_y_yz_z_xx, t_y_yz_z_xy,\
                           t_y_yz_z_xz, t_y_yz_z_yy, t_y_yz_z_yz, t_y_yz_z_zz, t_y_yzz_x_xx,\
                           t_y_yzz_x_xy, t_y_yzz_x_xz, t_y_yzz_x_yy, t_y_yzz_x_yz, t_y_yzz_x_zz,\
                           t_y_yzz_y_xx, t_y_yzz_y_xy, t_y_yzz_y_xz, t_y_yzz_y_yy, t_y_yzz_y_yz,\
                           t_y_yzz_y_zz, t_y_yzz_z_xx, t_y_yzz_z_xy, t_y_yzz_z_xz, t_y_yzz_z_yy,\
                           t_y_yzz_z_yz, t_y_yzz_z_zz, t_y_zz_x_xx, t_y_zz_x_xy, t_y_zz_x_xz,\
                           t_y_zz_x_yy, t_y_zz_x_yz, t_y_zz_x_zz, t_y_zz_y_xx, t_y_zz_y_xy,\
                           t_y_zz_y_xz, t_y_zz_y_yy, t_y_zz_y_yz, t_y_zz_y_zz, t_y_zz_z_xx,\
                           t_y_zz_z_xy, t_y_zz_z_xz, t_y_zz_z_yy, t_y_zz_z_yz, t_y_zz_z_zz,\
                           t_y_zzz_x_xx, t_y_zzz_x_xy, t_y_zzz_x_xz, t_y_zzz_x_yy, t_y_zzz_x_yz,\
                           t_y_zzz_x_zz, t_y_zzz_y_xx, t_y_zzz_y_xy, t_y_zzz_y_xz, t_y_zzz_y_yy,\
                           t_y_zzz_y_yz, t_y_zzz_y_zz, t_y_zzz_z_xx, t_y_zzz_z_xy, t_y_zzz_z_xz,\
                           t_y_zzz_z_yy, t_y_zzz_z_yz, t_y_zzz_z_zz, t_yz_yz_x_xx, t_yz_yz_x_xy,\
                           t_yz_yz_x_xz, t_yz_yz_x_yy, t_yz_yz_x_yz, t_yz_yz_x_zz, t_yz_yz_y_xx,\
                           t_yz_yz_y_xy, t_yz_yz_y_xz, t_yz_yz_y_yy, t_yz_yz_y_yz, t_yz_yz_y_zz,\
                           t_yz_yz_z_xx, t_yz_yz_z_xy, t_yz_yz_z_xz, t_yz_yz_z_yy, t_yz_yz_z_yz,\
                           t_yz_yz_z_zz, t_yz_zz_x_xx, t_yz_zz_x_xy, t_yz_zz_x_xz, t_yz_zz_x_yy,\
                           t_yz_zz_x_yz, t_yz_zz_x_zz, t_yz_zz_y_xx, t_yz_zz_y_xy, t_yz_zz_y_xz,\
                           t_yz_zz_y_yy, t_yz_zz_y_yz, t_yz_zz_y_zz, t_yz_zz_z_xx, t_yz_zz_z_xy,\
                           t_yz_zz_z_xz, t_yz_zz_z_yy, t_yz_zz_z_yz, t_yz_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_zz_z_zz[i] = t_y_zzz_z_zz[i] - rab_z[i] * t_y_zz_z_zz[i];

        t_yz_zz_z_yz[i] = t_y_zzz_z_yz[i] - rab_z[i] * t_y_zz_z_yz[i];

        t_yz_zz_z_yy[i] = t_y_zzz_z_yy[i] - rab_z[i] * t_y_zz_z_yy[i];

        t_yz_zz_z_xz[i] = t_y_zzz_z_xz[i] - rab_z[i] * t_y_zz_z_xz[i];

        t_yz_zz_z_xy[i] = t_y_zzz_z_xy[i] - rab_z[i] * t_y_zz_z_xy[i];

        t_yz_zz_z_xx[i] = t_y_zzz_z_xx[i] - rab_z[i] * t_y_zz_z_xx[i];

        t_yz_zz_y_zz[i] = t_y_zzz_y_zz[i] - rab_z[i] * t_y_zz_y_zz[i];

        t_yz_zz_y_yz[i] = t_y_zzz_y_yz[i] - rab_z[i] * t_y_zz_y_yz[i];

        t_yz_zz_y_yy[i] = t_y_zzz_y_yy[i] - rab_z[i] * t_y_zz_y_yy[i];

        t_yz_zz_y_xz[i] = t_y_zzz_y_xz[i] - rab_z[i] * t_y_zz_y_xz[i];

        t_yz_zz_y_xy[i] = t_y_zzz_y_xy[i] - rab_z[i] * t_y_zz_y_xy[i];

        t_yz_zz_y_xx[i] = t_y_zzz_y_xx[i] - rab_z[i] * t_y_zz_y_xx[i];

        t_yz_zz_x_zz[i] = t_y_zzz_x_zz[i] - rab_z[i] * t_y_zz_x_zz[i];

        t_yz_zz_x_yz[i] = t_y_zzz_x_yz[i] - rab_z[i] * t_y_zz_x_yz[i];

        t_yz_zz_x_yy[i] = t_y_zzz_x_yy[i] - rab_z[i] * t_y_zz_x_yy[i];

        t_yz_zz_x_xz[i] = t_y_zzz_x_xz[i] - rab_z[i] * t_y_zz_x_xz[i];

        t_yz_zz_x_xy[i] = t_y_zzz_x_xy[i] - rab_z[i] * t_y_zz_x_xy[i];

        t_yz_zz_x_xx[i] = t_y_zzz_x_xx[i] - rab_z[i] * t_y_zz_x_xx[i];

        t_yz_yz_z_zz[i] = t_y_yzz_z_zz[i] - rab_z[i] * t_y_yz_z_zz[i];

        t_yz_yz_z_yz[i] = t_y_yzz_z_yz[i] - rab_z[i] * t_y_yz_z_yz[i];

        t_yz_yz_z_yy[i] = t_y_yzz_z_yy[i] - rab_z[i] * t_y_yz_z_yy[i];

        t_yz_yz_z_xz[i] = t_y_yzz_z_xz[i] - rab_z[i] * t_y_yz_z_xz[i];

        t_yz_yz_z_xy[i] = t_y_yzz_z_xy[i] - rab_z[i] * t_y_yz_z_xy[i];

        t_yz_yz_z_xx[i] = t_y_yzz_z_xx[i] - rab_z[i] * t_y_yz_z_xx[i];

        t_yz_yz_y_zz[i] = t_y_yzz_y_zz[i] - rab_z[i] * t_y_yz_y_zz[i];

        t_yz_yz_y_yz[i] = t_y_yzz_y_yz[i] - rab_z[i] * t_y_yz_y_yz[i];

        t_yz_yz_y_yy[i] = t_y_yzz_y_yy[i] - rab_z[i] * t_y_yz_y_yy[i];

        t_yz_yz_y_xz[i] = t_y_yzz_y_xz[i] - rab_z[i] * t_y_yz_y_xz[i];

        t_yz_yz_y_xy[i] = t_y_yzz_y_xy[i] - rab_z[i] * t_y_yz_y_xy[i];

        t_yz_yz_y_xx[i] = t_y_yzz_y_xx[i] - rab_z[i] * t_y_yz_y_xx[i];

        t_yz_yz_x_zz[i] = t_y_yzz_x_zz[i] - rab_z[i] * t_y_yz_x_zz[i];

        t_yz_yz_x_yz[i] = t_y_yzz_x_yz[i] - rab_z[i] * t_y_yz_x_yz[i];

        t_yz_yz_x_yy[i] = t_y_yzz_x_yy[i] - rab_z[i] * t_y_yz_x_yy[i];

        t_yz_yz_x_xz[i] = t_y_yzz_x_xz[i] - rab_z[i] * t_y_yz_x_xz[i];

        t_yz_yz_x_xy[i] = t_y_yzz_x_xy[i] - rab_z[i] * t_y_yz_x_xy[i];

        t_yz_yz_x_xx[i] = t_y_yzz_x_xx[i] - rab_z[i] * t_y_yz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xz_x_xx, t_y_xz_x_xy, t_y_xz_x_xz, t_y_xz_x_yy,\
                           t_y_xz_x_yz, t_y_xz_x_zz, t_y_xz_y_xx, t_y_xz_y_xy, t_y_xz_y_xz,\
                           t_y_xz_y_yy, t_y_xz_y_yz, t_y_xz_y_zz, t_y_xz_z_xx, t_y_xz_z_xy,\
                           t_y_xz_z_xz, t_y_xz_z_yy, t_y_xz_z_yz, t_y_xz_z_zz, t_y_xzz_x_xx,\
                           t_y_xzz_x_xy, t_y_xzz_x_xz, t_y_xzz_x_yy, t_y_xzz_x_yz, t_y_xzz_x_zz,\
                           t_y_xzz_y_xx, t_y_xzz_y_xy, t_y_xzz_y_xz, t_y_xzz_y_yy, t_y_xzz_y_yz,\
                           t_y_xzz_y_zz, t_y_xzz_z_xx, t_y_xzz_z_xy, t_y_xzz_z_xz, t_y_xzz_z_yy,\
                           t_y_xzz_z_yz, t_y_xzz_z_zz, t_y_yy_x_xx, t_y_yy_x_xy, t_y_yy_x_xz,\
                           t_y_yy_x_yy, t_y_yy_x_yz, t_y_yy_x_zz, t_y_yy_y_xx, t_y_yy_y_xy,\
                           t_y_yy_y_xz, t_y_yy_y_yy, t_y_yy_y_yz, t_y_yy_y_zz, t_y_yy_z_xx,\
                           t_y_yy_z_xy, t_y_yy_z_xz, t_y_yy_z_yy, t_y_yy_z_yz, t_y_yy_z_zz,\
                           t_y_yyz_x_xx, t_y_yyz_x_xy, t_y_yyz_x_xz, t_y_yyz_x_yy, t_y_yyz_x_yz,\
                           t_y_yyz_x_zz, t_y_yyz_y_xx, t_y_yyz_y_xy, t_y_yyz_y_xz, t_y_yyz_y_yy,\
                           t_y_yyz_y_yz, t_y_yyz_y_zz, t_y_yyz_z_xx, t_y_yyz_z_xy, t_y_yyz_z_xz,\
                           t_y_yyz_z_yy, t_y_yyz_z_yz, t_y_yyz_z_zz, t_yz_xz_x_xx, t_yz_xz_x_xy,\
                           t_yz_xz_x_xz, t_yz_xz_x_yy, t_yz_xz_x_yz, t_yz_xz_x_zz, t_yz_xz_y_xx,\
                           t_yz_xz_y_xy, t_yz_xz_y_xz, t_yz_xz_y_yy, t_yz_xz_y_yz, t_yz_xz_y_zz,\
                           t_yz_xz_z_xx, t_yz_xz_z_xy, t_yz_xz_z_xz, t_yz_xz_z_yy, t_yz_xz_z_yz,\
                           t_yz_xz_z_zz, t_yz_yy_x_xx, t_yz_yy_x_xy, t_yz_yy_x_xz, t_yz_yy_x_yy,\
                           t_yz_yy_x_yz, t_yz_yy_x_zz, t_yz_yy_y_xx, t_yz_yy_y_xy, t_yz_yy_y_xz,\
                           t_yz_yy_y_yy, t_yz_yy_y_yz, t_yz_yy_y_zz, t_yz_yy_z_xx, t_yz_yy_z_xy,\
                           t_yz_yy_z_xz, t_yz_yy_z_yy, t_yz_yy_z_yz, t_yz_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_yy_z_zz[i] = t_y_yyz_z_zz[i] - rab_z[i] * t_y_yy_z_zz[i];

        t_yz_yy_z_yz[i] = t_y_yyz_z_yz[i] - rab_z[i] * t_y_yy_z_yz[i];

        t_yz_yy_z_yy[i] = t_y_yyz_z_yy[i] - rab_z[i] * t_y_yy_z_yy[i];

        t_yz_yy_z_xz[i] = t_y_yyz_z_xz[i] - rab_z[i] * t_y_yy_z_xz[i];

        t_yz_yy_z_xy[i] = t_y_yyz_z_xy[i] - rab_z[i] * t_y_yy_z_xy[i];

        t_yz_yy_z_xx[i] = t_y_yyz_z_xx[i] - rab_z[i] * t_y_yy_z_xx[i];

        t_yz_yy_y_zz[i] = t_y_yyz_y_zz[i] - rab_z[i] * t_y_yy_y_zz[i];

        t_yz_yy_y_yz[i] = t_y_yyz_y_yz[i] - rab_z[i] * t_y_yy_y_yz[i];

        t_yz_yy_y_yy[i] = t_y_yyz_y_yy[i] - rab_z[i] * t_y_yy_y_yy[i];

        t_yz_yy_y_xz[i] = t_y_yyz_y_xz[i] - rab_z[i] * t_y_yy_y_xz[i];

        t_yz_yy_y_xy[i] = t_y_yyz_y_xy[i] - rab_z[i] * t_y_yy_y_xy[i];

        t_yz_yy_y_xx[i] = t_y_yyz_y_xx[i] - rab_z[i] * t_y_yy_y_xx[i];

        t_yz_yy_x_zz[i] = t_y_yyz_x_zz[i] - rab_z[i] * t_y_yy_x_zz[i];

        t_yz_yy_x_yz[i] = t_y_yyz_x_yz[i] - rab_z[i] * t_y_yy_x_yz[i];

        t_yz_yy_x_yy[i] = t_y_yyz_x_yy[i] - rab_z[i] * t_y_yy_x_yy[i];

        t_yz_yy_x_xz[i] = t_y_yyz_x_xz[i] - rab_z[i] * t_y_yy_x_xz[i];

        t_yz_yy_x_xy[i] = t_y_yyz_x_xy[i] - rab_z[i] * t_y_yy_x_xy[i];

        t_yz_yy_x_xx[i] = t_y_yyz_x_xx[i] - rab_z[i] * t_y_yy_x_xx[i];

        t_yz_xz_z_zz[i] = t_y_xzz_z_zz[i] - rab_z[i] * t_y_xz_z_zz[i];

        t_yz_xz_z_yz[i] = t_y_xzz_z_yz[i] - rab_z[i] * t_y_xz_z_yz[i];

        t_yz_xz_z_yy[i] = t_y_xzz_z_yy[i] - rab_z[i] * t_y_xz_z_yy[i];

        t_yz_xz_z_xz[i] = t_y_xzz_z_xz[i] - rab_z[i] * t_y_xz_z_xz[i];

        t_yz_xz_z_xy[i] = t_y_xzz_z_xy[i] - rab_z[i] * t_y_xz_z_xy[i];

        t_yz_xz_z_xx[i] = t_y_xzz_z_xx[i] - rab_z[i] * t_y_xz_z_xx[i];

        t_yz_xz_y_zz[i] = t_y_xzz_y_zz[i] - rab_z[i] * t_y_xz_y_zz[i];

        t_yz_xz_y_yz[i] = t_y_xzz_y_yz[i] - rab_z[i] * t_y_xz_y_yz[i];

        t_yz_xz_y_yy[i] = t_y_xzz_y_yy[i] - rab_z[i] * t_y_xz_y_yy[i];

        t_yz_xz_y_xz[i] = t_y_xzz_y_xz[i] - rab_z[i] * t_y_xz_y_xz[i];

        t_yz_xz_y_xy[i] = t_y_xzz_y_xy[i] - rab_z[i] * t_y_xz_y_xy[i];

        t_yz_xz_y_xx[i] = t_y_xzz_y_xx[i] - rab_z[i] * t_y_xz_y_xx[i];

        t_yz_xz_x_zz[i] = t_y_xzz_x_zz[i] - rab_z[i] * t_y_xz_x_zz[i];

        t_yz_xz_x_yz[i] = t_y_xzz_x_yz[i] - rab_z[i] * t_y_xz_x_yz[i];

        t_yz_xz_x_yy[i] = t_y_xzz_x_yy[i] - rab_z[i] * t_y_xz_x_yy[i];

        t_yz_xz_x_xz[i] = t_y_xzz_x_xz[i] - rab_z[i] * t_y_xz_x_xz[i];

        t_yz_xz_x_xy[i] = t_y_xzz_x_xy[i] - rab_z[i] * t_y_xz_x_xy[i];

        t_yz_xz_x_xx[i] = t_y_xzz_x_xx[i] - rab_z[i] * t_y_xz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_y_xx_x_xx, t_y_xx_x_xy, t_y_xx_x_xz, t_y_xx_x_yy,\
                           t_y_xx_x_yz, t_y_xx_x_zz, t_y_xx_y_xx, t_y_xx_y_xy, t_y_xx_y_xz,\
                           t_y_xx_y_yy, t_y_xx_y_yz, t_y_xx_y_zz, t_y_xx_z_xx, t_y_xx_z_xy,\
                           t_y_xx_z_xz, t_y_xx_z_yy, t_y_xx_z_yz, t_y_xx_z_zz, t_y_xxz_x_xx,\
                           t_y_xxz_x_xy, t_y_xxz_x_xz, t_y_xxz_x_yy, t_y_xxz_x_yz, t_y_xxz_x_zz,\
                           t_y_xxz_y_xx, t_y_xxz_y_xy, t_y_xxz_y_xz, t_y_xxz_y_yy, t_y_xxz_y_yz,\
                           t_y_xxz_y_zz, t_y_xxz_z_xx, t_y_xxz_z_xy, t_y_xxz_z_xz, t_y_xxz_z_yy,\
                           t_y_xxz_z_yz, t_y_xxz_z_zz, t_y_xy_x_xx, t_y_xy_x_xy, t_y_xy_x_xz,\
                           t_y_xy_x_yy, t_y_xy_x_yz, t_y_xy_x_zz, t_y_xy_y_xx, t_y_xy_y_xy,\
                           t_y_xy_y_xz, t_y_xy_y_yy, t_y_xy_y_yz, t_y_xy_y_zz, t_y_xy_z_xx,\
                           t_y_xy_z_xy, t_y_xy_z_xz, t_y_xy_z_yy, t_y_xy_z_yz, t_y_xy_z_zz,\
                           t_y_xyz_x_xx, t_y_xyz_x_xy, t_y_xyz_x_xz, t_y_xyz_x_yy, t_y_xyz_x_yz,\
                           t_y_xyz_x_zz, t_y_xyz_y_xx, t_y_xyz_y_xy, t_y_xyz_y_xz, t_y_xyz_y_yy,\
                           t_y_xyz_y_yz, t_y_xyz_y_zz, t_y_xyz_z_xx, t_y_xyz_z_xy, t_y_xyz_z_xz,\
                           t_y_xyz_z_yy, t_y_xyz_z_yz, t_y_xyz_z_zz, t_yz_xx_x_xx, t_yz_xx_x_xy,\
                           t_yz_xx_x_xz, t_yz_xx_x_yy, t_yz_xx_x_yz, t_yz_xx_x_zz, t_yz_xx_y_xx,\
                           t_yz_xx_y_xy, t_yz_xx_y_xz, t_yz_xx_y_yy, t_yz_xx_y_yz, t_yz_xx_y_zz,\
                           t_yz_xx_z_xx, t_yz_xx_z_xy, t_yz_xx_z_xz, t_yz_xx_z_yy, t_yz_xx_z_yz,\
                           t_yz_xx_z_zz, t_yz_xy_x_xx, t_yz_xy_x_xy, t_yz_xy_x_xz, t_yz_xy_x_yy,\
                           t_yz_xy_x_yz, t_yz_xy_x_zz, t_yz_xy_y_xx, t_yz_xy_y_xy, t_yz_xy_y_xz,\
                           t_yz_xy_y_yy, t_yz_xy_y_yz, t_yz_xy_y_zz, t_yz_xy_z_xx, t_yz_xy_z_xy,\
                           t_yz_xy_z_xz, t_yz_xy_z_yy, t_yz_xy_z_yz, t_yz_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_xy_z_zz[i] = t_y_xyz_z_zz[i] - rab_z[i] * t_y_xy_z_zz[i];

        t_yz_xy_z_yz[i] = t_y_xyz_z_yz[i] - rab_z[i] * t_y_xy_z_yz[i];

        t_yz_xy_z_yy[i] = t_y_xyz_z_yy[i] - rab_z[i] * t_y_xy_z_yy[i];

        t_yz_xy_z_xz[i] = t_y_xyz_z_xz[i] - rab_z[i] * t_y_xy_z_xz[i];

        t_yz_xy_z_xy[i] = t_y_xyz_z_xy[i] - rab_z[i] * t_y_xy_z_xy[i];

        t_yz_xy_z_xx[i] = t_y_xyz_z_xx[i] - rab_z[i] * t_y_xy_z_xx[i];

        t_yz_xy_y_zz[i] = t_y_xyz_y_zz[i] - rab_z[i] * t_y_xy_y_zz[i];

        t_yz_xy_y_yz[i] = t_y_xyz_y_yz[i] - rab_z[i] * t_y_xy_y_yz[i];

        t_yz_xy_y_yy[i] = t_y_xyz_y_yy[i] - rab_z[i] * t_y_xy_y_yy[i];

        t_yz_xy_y_xz[i] = t_y_xyz_y_xz[i] - rab_z[i] * t_y_xy_y_xz[i];

        t_yz_xy_y_xy[i] = t_y_xyz_y_xy[i] - rab_z[i] * t_y_xy_y_xy[i];

        t_yz_xy_y_xx[i] = t_y_xyz_y_xx[i] - rab_z[i] * t_y_xy_y_xx[i];

        t_yz_xy_x_zz[i] = t_y_xyz_x_zz[i] - rab_z[i] * t_y_xy_x_zz[i];

        t_yz_xy_x_yz[i] = t_y_xyz_x_yz[i] - rab_z[i] * t_y_xy_x_yz[i];

        t_yz_xy_x_yy[i] = t_y_xyz_x_yy[i] - rab_z[i] * t_y_xy_x_yy[i];

        t_yz_xy_x_xz[i] = t_y_xyz_x_xz[i] - rab_z[i] * t_y_xy_x_xz[i];

        t_yz_xy_x_xy[i] = t_y_xyz_x_xy[i] - rab_z[i] * t_y_xy_x_xy[i];

        t_yz_xy_x_xx[i] = t_y_xyz_x_xx[i] - rab_z[i] * t_y_xy_x_xx[i];

        t_yz_xx_z_zz[i] = t_y_xxz_z_zz[i] - rab_z[i] * t_y_xx_z_zz[i];

        t_yz_xx_z_yz[i] = t_y_xxz_z_yz[i] - rab_z[i] * t_y_xx_z_yz[i];

        t_yz_xx_z_yy[i] = t_y_xxz_z_yy[i] - rab_z[i] * t_y_xx_z_yy[i];

        t_yz_xx_z_xz[i] = t_y_xxz_z_xz[i] - rab_z[i] * t_y_xx_z_xz[i];

        t_yz_xx_z_xy[i] = t_y_xxz_z_xy[i] - rab_z[i] * t_y_xx_z_xy[i];

        t_yz_xx_z_xx[i] = t_y_xxz_z_xx[i] - rab_z[i] * t_y_xx_z_xx[i];

        t_yz_xx_y_zz[i] = t_y_xxz_y_zz[i] - rab_z[i] * t_y_xx_y_zz[i];

        t_yz_xx_y_yz[i] = t_y_xxz_y_yz[i] - rab_z[i] * t_y_xx_y_yz[i];

        t_yz_xx_y_yy[i] = t_y_xxz_y_yy[i] - rab_z[i] * t_y_xx_y_yy[i];

        t_yz_xx_y_xz[i] = t_y_xxz_y_xz[i] - rab_z[i] * t_y_xx_y_xz[i];

        t_yz_xx_y_xy[i] = t_y_xxz_y_xy[i] - rab_z[i] * t_y_xx_y_xy[i];

        t_yz_xx_y_xx[i] = t_y_xxz_y_xx[i] - rab_z[i] * t_y_xx_y_xx[i];

        t_yz_xx_x_zz[i] = t_y_xxz_x_zz[i] - rab_z[i] * t_y_xx_x_zz[i];

        t_yz_xx_x_yz[i] = t_y_xxz_x_yz[i] - rab_z[i] * t_y_xx_x_yz[i];

        t_yz_xx_x_yy[i] = t_y_xxz_x_yy[i] - rab_z[i] * t_y_xx_x_yy[i];

        t_yz_xx_x_xz[i] = t_y_xxz_x_xz[i] - rab_z[i] * t_y_xx_x_xz[i];

        t_yz_xx_x_xy[i] = t_y_xxz_x_xy[i] - rab_z[i] * t_y_xx_x_xy[i];

        t_yz_xx_x_xx[i] = t_y_xxz_x_xx[i] - rab_z[i] * t_y_xx_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_yyz_x_xx, t_y_yyz_x_xy, t_y_yyz_x_xz, t_y_yyz_x_yy,\
                           t_y_yyz_x_yz, t_y_yyz_x_zz, t_y_yyz_y_xx, t_y_yyz_y_xy, t_y_yyz_y_xz,\
                           t_y_yyz_y_yy, t_y_yyz_y_yz, t_y_yyz_y_zz, t_y_yyz_z_xx, t_y_yyz_z_xy,\
                           t_y_yyz_z_xz, t_y_yyz_z_yy, t_y_yyz_z_yz, t_y_yyz_z_zz, t_y_yz_x_xx,\
                           t_y_yz_x_xy, t_y_yz_x_xz, t_y_yz_x_yy, t_y_yz_x_yz, t_y_yz_x_zz,\
                           t_y_yz_y_xx, t_y_yz_y_xy, t_y_yz_y_xz, t_y_yz_y_yy, t_y_yz_y_yz,\
                           t_y_yz_y_zz, t_y_yz_z_xx, t_y_yz_z_xy, t_y_yz_z_xz, t_y_yz_z_yy,\
                           t_y_yz_z_yz, t_y_yz_z_zz, t_y_yzz_x_xx, t_y_yzz_x_xy, t_y_yzz_x_xz,\
                           t_y_yzz_x_yy, t_y_yzz_x_yz, t_y_yzz_x_zz, t_y_yzz_y_xx, t_y_yzz_y_xy,\
                           t_y_yzz_y_xz, t_y_yzz_y_yy, t_y_yzz_y_yz, t_y_yzz_y_zz, t_y_yzz_z_xx,\
                           t_y_yzz_z_xy, t_y_yzz_z_xz, t_y_yzz_z_yy, t_y_yzz_z_yz, t_y_yzz_z_zz,\
                           t_y_zz_x_xx, t_y_zz_x_xy, t_y_zz_x_xz, t_y_zz_x_yy, t_y_zz_x_yz,\
                           t_y_zz_x_zz, t_y_zz_y_xx, t_y_zz_y_xy, t_y_zz_y_xz, t_y_zz_y_yy,\
                           t_y_zz_y_yz, t_y_zz_y_zz, t_y_zz_z_xx, t_y_zz_z_xy, t_y_zz_z_xz,\
                           t_y_zz_z_yy, t_y_zz_z_yz, t_y_zz_z_zz, t_yy_yz_x_xx, t_yy_yz_x_xy,\
                           t_yy_yz_x_xz, t_yy_yz_x_yy, t_yy_yz_x_yz, t_yy_yz_x_zz, t_yy_yz_y_xx,\
                           t_yy_yz_y_xy, t_yy_yz_y_xz, t_yy_yz_y_yy, t_yy_yz_y_yz, t_yy_yz_y_zz,\
                           t_yy_yz_z_xx, t_yy_yz_z_xy, t_yy_yz_z_xz, t_yy_yz_z_yy, t_yy_yz_z_yz,\
                           t_yy_yz_z_zz, t_yy_zz_x_xx, t_yy_zz_x_xy, t_yy_zz_x_xz, t_yy_zz_x_yy,\
                           t_yy_zz_x_yz, t_yy_zz_x_zz, t_yy_zz_y_xx, t_yy_zz_y_xy, t_yy_zz_y_xz,\
                           t_yy_zz_y_yy, t_yy_zz_y_yz, t_yy_zz_y_zz, t_yy_zz_z_xx, t_yy_zz_z_xy,\
                           t_yy_zz_z_xz, t_yy_zz_z_yy, t_yy_zz_z_yz, t_yy_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_zz_z_zz[i] = t_y_yzz_z_zz[i] - rab_y[i] * t_y_zz_z_zz[i];

        t_yy_zz_z_yz[i] = t_y_yzz_z_yz[i] - rab_y[i] * t_y_zz_z_yz[i];

        t_yy_zz_z_yy[i] = t_y_yzz_z_yy[i] - rab_y[i] * t_y_zz_z_yy[i];

        t_yy_zz_z_xz[i] = t_y_yzz_z_xz[i] - rab_y[i] * t_y_zz_z_xz[i];

        t_yy_zz_z_xy[i] = t_y_yzz_z_xy[i] - rab_y[i] * t_y_zz_z_xy[i];

        t_yy_zz_z_xx[i] = t_y_yzz_z_xx[i] - rab_y[i] * t_y_zz_z_xx[i];

        t_yy_zz_y_zz[i] = t_y_yzz_y_zz[i] - rab_y[i] * t_y_zz_y_zz[i];

        t_yy_zz_y_yz[i] = t_y_yzz_y_yz[i] - rab_y[i] * t_y_zz_y_yz[i];

        t_yy_zz_y_yy[i] = t_y_yzz_y_yy[i] - rab_y[i] * t_y_zz_y_yy[i];

        t_yy_zz_y_xz[i] = t_y_yzz_y_xz[i] - rab_y[i] * t_y_zz_y_xz[i];

        t_yy_zz_y_xy[i] = t_y_yzz_y_xy[i] - rab_y[i] * t_y_zz_y_xy[i];

        t_yy_zz_y_xx[i] = t_y_yzz_y_xx[i] - rab_y[i] * t_y_zz_y_xx[i];

        t_yy_zz_x_zz[i] = t_y_yzz_x_zz[i] - rab_y[i] * t_y_zz_x_zz[i];

        t_yy_zz_x_yz[i] = t_y_yzz_x_yz[i] - rab_y[i] * t_y_zz_x_yz[i];

        t_yy_zz_x_yy[i] = t_y_yzz_x_yy[i] - rab_y[i] * t_y_zz_x_yy[i];

        t_yy_zz_x_xz[i] = t_y_yzz_x_xz[i] - rab_y[i] * t_y_zz_x_xz[i];

        t_yy_zz_x_xy[i] = t_y_yzz_x_xy[i] - rab_y[i] * t_y_zz_x_xy[i];

        t_yy_zz_x_xx[i] = t_y_yzz_x_xx[i] - rab_y[i] * t_y_zz_x_xx[i];

        t_yy_yz_z_zz[i] = t_y_yyz_z_zz[i] - rab_y[i] * t_y_yz_z_zz[i];

        t_yy_yz_z_yz[i] = t_y_yyz_z_yz[i] - rab_y[i] * t_y_yz_z_yz[i];

        t_yy_yz_z_yy[i] = t_y_yyz_z_yy[i] - rab_y[i] * t_y_yz_z_yy[i];

        t_yy_yz_z_xz[i] = t_y_yyz_z_xz[i] - rab_y[i] * t_y_yz_z_xz[i];

        t_yy_yz_z_xy[i] = t_y_yyz_z_xy[i] - rab_y[i] * t_y_yz_z_xy[i];

        t_yy_yz_z_xx[i] = t_y_yyz_z_xx[i] - rab_y[i] * t_y_yz_z_xx[i];

        t_yy_yz_y_zz[i] = t_y_yyz_y_zz[i] - rab_y[i] * t_y_yz_y_zz[i];

        t_yy_yz_y_yz[i] = t_y_yyz_y_yz[i] - rab_y[i] * t_y_yz_y_yz[i];

        t_yy_yz_y_yy[i] = t_y_yyz_y_yy[i] - rab_y[i] * t_y_yz_y_yy[i];

        t_yy_yz_y_xz[i] = t_y_yyz_y_xz[i] - rab_y[i] * t_y_yz_y_xz[i];

        t_yy_yz_y_xy[i] = t_y_yyz_y_xy[i] - rab_y[i] * t_y_yz_y_xy[i];

        t_yy_yz_y_xx[i] = t_y_yyz_y_xx[i] - rab_y[i] * t_y_yz_y_xx[i];

        t_yy_yz_x_zz[i] = t_y_yyz_x_zz[i] - rab_y[i] * t_y_yz_x_zz[i];

        t_yy_yz_x_yz[i] = t_y_yyz_x_yz[i] - rab_y[i] * t_y_yz_x_yz[i];

        t_yy_yz_x_yy[i] = t_y_yyz_x_yy[i] - rab_y[i] * t_y_yz_x_yy[i];

        t_yy_yz_x_xz[i] = t_y_yyz_x_xz[i] - rab_y[i] * t_y_yz_x_xz[i];

        t_yy_yz_x_xy[i] = t_y_yyz_x_xy[i] - rab_y[i] * t_y_yz_x_xy[i];

        t_yy_yz_x_xx[i] = t_y_yyz_x_xx[i] - rab_y[i] * t_y_yz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xyz_x_xx, t_y_xyz_x_xy, t_y_xyz_x_xz, t_y_xyz_x_yy,\
                           t_y_xyz_x_yz, t_y_xyz_x_zz, t_y_xyz_y_xx, t_y_xyz_y_xy, t_y_xyz_y_xz,\
                           t_y_xyz_y_yy, t_y_xyz_y_yz, t_y_xyz_y_zz, t_y_xyz_z_xx, t_y_xyz_z_xy,\
                           t_y_xyz_z_xz, t_y_xyz_z_yy, t_y_xyz_z_yz, t_y_xyz_z_zz, t_y_xz_x_xx,\
                           t_y_xz_x_xy, t_y_xz_x_xz, t_y_xz_x_yy, t_y_xz_x_yz, t_y_xz_x_zz,\
                           t_y_xz_y_xx, t_y_xz_y_xy, t_y_xz_y_xz, t_y_xz_y_yy, t_y_xz_y_yz,\
                           t_y_xz_y_zz, t_y_xz_z_xx, t_y_xz_z_xy, t_y_xz_z_xz, t_y_xz_z_yy,\
                           t_y_xz_z_yz, t_y_xz_z_zz, t_y_yy_x_xx, t_y_yy_x_xy, t_y_yy_x_xz,\
                           t_y_yy_x_yy, t_y_yy_x_yz, t_y_yy_x_zz, t_y_yy_y_xx, t_y_yy_y_xy,\
                           t_y_yy_y_xz, t_y_yy_y_yy, t_y_yy_y_yz, t_y_yy_y_zz, t_y_yy_z_xx,\
                           t_y_yy_z_xy, t_y_yy_z_xz, t_y_yy_z_yy, t_y_yy_z_yz, t_y_yy_z_zz,\
                           t_y_yyy_x_xx, t_y_yyy_x_xy, t_y_yyy_x_xz, t_y_yyy_x_yy, t_y_yyy_x_yz,\
                           t_y_yyy_x_zz, t_y_yyy_y_xx, t_y_yyy_y_xy, t_y_yyy_y_xz, t_y_yyy_y_yy,\
                           t_y_yyy_y_yz, t_y_yyy_y_zz, t_y_yyy_z_xx, t_y_yyy_z_xy, t_y_yyy_z_xz,\
                           t_y_yyy_z_yy, t_y_yyy_z_yz, t_y_yyy_z_zz, t_yy_xz_x_xx, t_yy_xz_x_xy,\
                           t_yy_xz_x_xz, t_yy_xz_x_yy, t_yy_xz_x_yz, t_yy_xz_x_zz, t_yy_xz_y_xx,\
                           t_yy_xz_y_xy, t_yy_xz_y_xz, t_yy_xz_y_yy, t_yy_xz_y_yz, t_yy_xz_y_zz,\
                           t_yy_xz_z_xx, t_yy_xz_z_xy, t_yy_xz_z_xz, t_yy_xz_z_yy, t_yy_xz_z_yz,\
                           t_yy_xz_z_zz, t_yy_yy_x_xx, t_yy_yy_x_xy, t_yy_yy_x_xz, t_yy_yy_x_yy,\
                           t_yy_yy_x_yz, t_yy_yy_x_zz, t_yy_yy_y_xx, t_yy_yy_y_xy, t_yy_yy_y_xz,\
                           t_yy_yy_y_yy, t_yy_yy_y_yz, t_yy_yy_y_zz, t_yy_yy_z_xx, t_yy_yy_z_xy,\
                           t_yy_yy_z_xz, t_yy_yy_z_yy, t_yy_yy_z_yz, t_yy_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_yy_z_zz[i] = t_y_yyy_z_zz[i] - rab_y[i] * t_y_yy_z_zz[i];

        t_yy_yy_z_yz[i] = t_y_yyy_z_yz[i] - rab_y[i] * t_y_yy_z_yz[i];

        t_yy_yy_z_yy[i] = t_y_yyy_z_yy[i] - rab_y[i] * t_y_yy_z_yy[i];

        t_yy_yy_z_xz[i] = t_y_yyy_z_xz[i] - rab_y[i] * t_y_yy_z_xz[i];

        t_yy_yy_z_xy[i] = t_y_yyy_z_xy[i] - rab_y[i] * t_y_yy_z_xy[i];

        t_yy_yy_z_xx[i] = t_y_yyy_z_xx[i] - rab_y[i] * t_y_yy_z_xx[i];

        t_yy_yy_y_zz[i] = t_y_yyy_y_zz[i] - rab_y[i] * t_y_yy_y_zz[i];

        t_yy_yy_y_yz[i] = t_y_yyy_y_yz[i] - rab_y[i] * t_y_yy_y_yz[i];

        t_yy_yy_y_yy[i] = t_y_yyy_y_yy[i] - rab_y[i] * t_y_yy_y_yy[i];

        t_yy_yy_y_xz[i] = t_y_yyy_y_xz[i] - rab_y[i] * t_y_yy_y_xz[i];

        t_yy_yy_y_xy[i] = t_y_yyy_y_xy[i] - rab_y[i] * t_y_yy_y_xy[i];

        t_yy_yy_y_xx[i] = t_y_yyy_y_xx[i] - rab_y[i] * t_y_yy_y_xx[i];

        t_yy_yy_x_zz[i] = t_y_yyy_x_zz[i] - rab_y[i] * t_y_yy_x_zz[i];

        t_yy_yy_x_yz[i] = t_y_yyy_x_yz[i] - rab_y[i] * t_y_yy_x_yz[i];

        t_yy_yy_x_yy[i] = t_y_yyy_x_yy[i] - rab_y[i] * t_y_yy_x_yy[i];

        t_yy_yy_x_xz[i] = t_y_yyy_x_xz[i] - rab_y[i] * t_y_yy_x_xz[i];

        t_yy_yy_x_xy[i] = t_y_yyy_x_xy[i] - rab_y[i] * t_y_yy_x_xy[i];

        t_yy_yy_x_xx[i] = t_y_yyy_x_xx[i] - rab_y[i] * t_y_yy_x_xx[i];

        t_yy_xz_z_zz[i] = t_y_xyz_z_zz[i] - rab_y[i] * t_y_xz_z_zz[i];

        t_yy_xz_z_yz[i] = t_y_xyz_z_yz[i] - rab_y[i] * t_y_xz_z_yz[i];

        t_yy_xz_z_yy[i] = t_y_xyz_z_yy[i] - rab_y[i] * t_y_xz_z_yy[i];

        t_yy_xz_z_xz[i] = t_y_xyz_z_xz[i] - rab_y[i] * t_y_xz_z_xz[i];

        t_yy_xz_z_xy[i] = t_y_xyz_z_xy[i] - rab_y[i] * t_y_xz_z_xy[i];

        t_yy_xz_z_xx[i] = t_y_xyz_z_xx[i] - rab_y[i] * t_y_xz_z_xx[i];

        t_yy_xz_y_zz[i] = t_y_xyz_y_zz[i] - rab_y[i] * t_y_xz_y_zz[i];

        t_yy_xz_y_yz[i] = t_y_xyz_y_yz[i] - rab_y[i] * t_y_xz_y_yz[i];

        t_yy_xz_y_yy[i] = t_y_xyz_y_yy[i] - rab_y[i] * t_y_xz_y_yy[i];

        t_yy_xz_y_xz[i] = t_y_xyz_y_xz[i] - rab_y[i] * t_y_xz_y_xz[i];

        t_yy_xz_y_xy[i] = t_y_xyz_y_xy[i] - rab_y[i] * t_y_xz_y_xy[i];

        t_yy_xz_y_xx[i] = t_y_xyz_y_xx[i] - rab_y[i] * t_y_xz_y_xx[i];

        t_yy_xz_x_zz[i] = t_y_xyz_x_zz[i] - rab_y[i] * t_y_xz_x_zz[i];

        t_yy_xz_x_yz[i] = t_y_xyz_x_yz[i] - rab_y[i] * t_y_xz_x_yz[i];

        t_yy_xz_x_yy[i] = t_y_xyz_x_yy[i] - rab_y[i] * t_y_xz_x_yy[i];

        t_yy_xz_x_xz[i] = t_y_xyz_x_xz[i] - rab_y[i] * t_y_xz_x_xz[i];

        t_yy_xz_x_xy[i] = t_y_xyz_x_xy[i] - rab_y[i] * t_y_xz_x_xy[i];

        t_yy_xz_x_xx[i] = t_y_xyz_x_xx[i] - rab_y[i] * t_y_xz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_y_xx_x_xx, t_y_xx_x_xy, t_y_xx_x_xz, t_y_xx_x_yy,\
                           t_y_xx_x_yz, t_y_xx_x_zz, t_y_xx_y_xx, t_y_xx_y_xy, t_y_xx_y_xz,\
                           t_y_xx_y_yy, t_y_xx_y_yz, t_y_xx_y_zz, t_y_xx_z_xx, t_y_xx_z_xy,\
                           t_y_xx_z_xz, t_y_xx_z_yy, t_y_xx_z_yz, t_y_xx_z_zz, t_y_xxy_x_xx,\
                           t_y_xxy_x_xy, t_y_xxy_x_xz, t_y_xxy_x_yy, t_y_xxy_x_yz, t_y_xxy_x_zz,\
                           t_y_xxy_y_xx, t_y_xxy_y_xy, t_y_xxy_y_xz, t_y_xxy_y_yy, t_y_xxy_y_yz,\
                           t_y_xxy_y_zz, t_y_xxy_z_xx, t_y_xxy_z_xy, t_y_xxy_z_xz, t_y_xxy_z_yy,\
                           t_y_xxy_z_yz, t_y_xxy_z_zz, t_y_xy_x_xx, t_y_xy_x_xy, t_y_xy_x_xz,\
                           t_y_xy_x_yy, t_y_xy_x_yz, t_y_xy_x_zz, t_y_xy_y_xx, t_y_xy_y_xy,\
                           t_y_xy_y_xz, t_y_xy_y_yy, t_y_xy_y_yz, t_y_xy_y_zz, t_y_xy_z_xx,\
                           t_y_xy_z_xy, t_y_xy_z_xz, t_y_xy_z_yy, t_y_xy_z_yz, t_y_xy_z_zz,\
                           t_y_xyy_x_xx, t_y_xyy_x_xy, t_y_xyy_x_xz, t_y_xyy_x_yy, t_y_xyy_x_yz,\
                           t_y_xyy_x_zz, t_y_xyy_y_xx, t_y_xyy_y_xy, t_y_xyy_y_xz, t_y_xyy_y_yy,\
                           t_y_xyy_y_yz, t_y_xyy_y_zz, t_y_xyy_z_xx, t_y_xyy_z_xy, t_y_xyy_z_xz,\
                           t_y_xyy_z_yy, t_y_xyy_z_yz, t_y_xyy_z_zz, t_yy_xx_x_xx, t_yy_xx_x_xy,\
                           t_yy_xx_x_xz, t_yy_xx_x_yy, t_yy_xx_x_yz, t_yy_xx_x_zz, t_yy_xx_y_xx,\
                           t_yy_xx_y_xy, t_yy_xx_y_xz, t_yy_xx_y_yy, t_yy_xx_y_yz, t_yy_xx_y_zz,\
                           t_yy_xx_z_xx, t_yy_xx_z_xy, t_yy_xx_z_xz, t_yy_xx_z_yy, t_yy_xx_z_yz,\
                           t_yy_xx_z_zz, t_yy_xy_x_xx, t_yy_xy_x_xy, t_yy_xy_x_xz, t_yy_xy_x_yy,\
                           t_yy_xy_x_yz, t_yy_xy_x_zz, t_yy_xy_y_xx, t_yy_xy_y_xy, t_yy_xy_y_xz,\
                           t_yy_xy_y_yy, t_yy_xy_y_yz, t_yy_xy_y_zz, t_yy_xy_z_xx, t_yy_xy_z_xy,\
                           t_yy_xy_z_xz, t_yy_xy_z_yy, t_yy_xy_z_yz, t_yy_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_xy_z_zz[i] = t_y_xyy_z_zz[i] - rab_y[i] * t_y_xy_z_zz[i];

        t_yy_xy_z_yz[i] = t_y_xyy_z_yz[i] - rab_y[i] * t_y_xy_z_yz[i];

        t_yy_xy_z_yy[i] = t_y_xyy_z_yy[i] - rab_y[i] * t_y_xy_z_yy[i];

        t_yy_xy_z_xz[i] = t_y_xyy_z_xz[i] - rab_y[i] * t_y_xy_z_xz[i];

        t_yy_xy_z_xy[i] = t_y_xyy_z_xy[i] - rab_y[i] * t_y_xy_z_xy[i];

        t_yy_xy_z_xx[i] = t_y_xyy_z_xx[i] - rab_y[i] * t_y_xy_z_xx[i];

        t_yy_xy_y_zz[i] = t_y_xyy_y_zz[i] - rab_y[i] * t_y_xy_y_zz[i];

        t_yy_xy_y_yz[i] = t_y_xyy_y_yz[i] - rab_y[i] * t_y_xy_y_yz[i];

        t_yy_xy_y_yy[i] = t_y_xyy_y_yy[i] - rab_y[i] * t_y_xy_y_yy[i];

        t_yy_xy_y_xz[i] = t_y_xyy_y_xz[i] - rab_y[i] * t_y_xy_y_xz[i];

        t_yy_xy_y_xy[i] = t_y_xyy_y_xy[i] - rab_y[i] * t_y_xy_y_xy[i];

        t_yy_xy_y_xx[i] = t_y_xyy_y_xx[i] - rab_y[i] * t_y_xy_y_xx[i];

        t_yy_xy_x_zz[i] = t_y_xyy_x_zz[i] - rab_y[i] * t_y_xy_x_zz[i];

        t_yy_xy_x_yz[i] = t_y_xyy_x_yz[i] - rab_y[i] * t_y_xy_x_yz[i];

        t_yy_xy_x_yy[i] = t_y_xyy_x_yy[i] - rab_y[i] * t_y_xy_x_yy[i];

        t_yy_xy_x_xz[i] = t_y_xyy_x_xz[i] - rab_y[i] * t_y_xy_x_xz[i];

        t_yy_xy_x_xy[i] = t_y_xyy_x_xy[i] - rab_y[i] * t_y_xy_x_xy[i];

        t_yy_xy_x_xx[i] = t_y_xyy_x_xx[i] - rab_y[i] * t_y_xy_x_xx[i];

        t_yy_xx_z_zz[i] = t_y_xxy_z_zz[i] - rab_y[i] * t_y_xx_z_zz[i];

        t_yy_xx_z_yz[i] = t_y_xxy_z_yz[i] - rab_y[i] * t_y_xx_z_yz[i];

        t_yy_xx_z_yy[i] = t_y_xxy_z_yy[i] - rab_y[i] * t_y_xx_z_yy[i];

        t_yy_xx_z_xz[i] = t_y_xxy_z_xz[i] - rab_y[i] * t_y_xx_z_xz[i];

        t_yy_xx_z_xy[i] = t_y_xxy_z_xy[i] - rab_y[i] * t_y_xx_z_xy[i];

        t_yy_xx_z_xx[i] = t_y_xxy_z_xx[i] - rab_y[i] * t_y_xx_z_xx[i];

        t_yy_xx_y_zz[i] = t_y_xxy_y_zz[i] - rab_y[i] * t_y_xx_y_zz[i];

        t_yy_xx_y_yz[i] = t_y_xxy_y_yz[i] - rab_y[i] * t_y_xx_y_yz[i];

        t_yy_xx_y_yy[i] = t_y_xxy_y_yy[i] - rab_y[i] * t_y_xx_y_yy[i];

        t_yy_xx_y_xz[i] = t_y_xxy_y_xz[i] - rab_y[i] * t_y_xx_y_xz[i];

        t_yy_xx_y_xy[i] = t_y_xxy_y_xy[i] - rab_y[i] * t_y_xx_y_xy[i];

        t_yy_xx_y_xx[i] = t_y_xxy_y_xx[i] - rab_y[i] * t_y_xx_y_xx[i];

        t_yy_xx_x_zz[i] = t_y_xxy_x_zz[i] - rab_y[i] * t_y_xx_x_zz[i];

        t_yy_xx_x_yz[i] = t_y_xxy_x_yz[i] - rab_y[i] * t_y_xx_x_yz[i];

        t_yy_xx_x_yy[i] = t_y_xxy_x_yy[i] - rab_y[i] * t_y_xx_x_yy[i];

        t_yy_xx_x_xz[i] = t_y_xxy_x_xz[i] - rab_y[i] * t_y_xx_x_xz[i];

        t_yy_xx_x_xy[i] = t_y_xxy_x_xy[i] - rab_y[i] * t_y_xx_x_xy[i];

        t_yy_xx_x_xx[i] = t_y_xxy_x_xx[i] - rab_y[i] * t_y_xx_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_yz_x_xx, t_x_yz_x_xy, t_x_yz_x_xz, t_x_yz_x_yy,\
                           t_x_yz_x_yz, t_x_yz_x_zz, t_x_yz_y_xx, t_x_yz_y_xy, t_x_yz_y_xz,\
                           t_x_yz_y_yy, t_x_yz_y_yz, t_x_yz_y_zz, t_x_yz_z_xx, t_x_yz_z_xy,\
                           t_x_yz_z_xz, t_x_yz_z_yy, t_x_yz_z_yz, t_x_yz_z_zz, t_x_yzz_x_xx,\
                           t_x_yzz_x_xy, t_x_yzz_x_xz, t_x_yzz_x_yy, t_x_yzz_x_yz, t_x_yzz_x_zz,\
                           t_x_yzz_y_xx, t_x_yzz_y_xy, t_x_yzz_y_xz, t_x_yzz_y_yy, t_x_yzz_y_yz,\
                           t_x_yzz_y_zz, t_x_yzz_z_xx, t_x_yzz_z_xy, t_x_yzz_z_xz, t_x_yzz_z_yy,\
                           t_x_yzz_z_yz, t_x_yzz_z_zz, t_x_zz_x_xx, t_x_zz_x_xy, t_x_zz_x_xz,\
                           t_x_zz_x_yy, t_x_zz_x_yz, t_x_zz_x_zz, t_x_zz_y_xx, t_x_zz_y_xy,\
                           t_x_zz_y_xz, t_x_zz_y_yy, t_x_zz_y_yz, t_x_zz_y_zz, t_x_zz_z_xx,\
                           t_x_zz_z_xy, t_x_zz_z_xz, t_x_zz_z_yy, t_x_zz_z_yz, t_x_zz_z_zz,\
                           t_x_zzz_x_xx, t_x_zzz_x_xy, t_x_zzz_x_xz, t_x_zzz_x_yy, t_x_zzz_x_yz,\
                           t_x_zzz_x_zz, t_x_zzz_y_xx, t_x_zzz_y_xy, t_x_zzz_y_xz, t_x_zzz_y_yy,\
                           t_x_zzz_y_yz, t_x_zzz_y_zz, t_x_zzz_z_xx, t_x_zzz_z_xy, t_x_zzz_z_xz,\
                           t_x_zzz_z_yy, t_x_zzz_z_yz, t_x_zzz_z_zz, t_xz_yz_x_xx, t_xz_yz_x_xy,\
                           t_xz_yz_x_xz, t_xz_yz_x_yy, t_xz_yz_x_yz, t_xz_yz_x_zz, t_xz_yz_y_xx,\
                           t_xz_yz_y_xy, t_xz_yz_y_xz, t_xz_yz_y_yy, t_xz_yz_y_yz, t_xz_yz_y_zz,\
                           t_xz_yz_z_xx, t_xz_yz_z_xy, t_xz_yz_z_xz, t_xz_yz_z_yy, t_xz_yz_z_yz,\
                           t_xz_yz_z_zz, t_xz_zz_x_xx, t_xz_zz_x_xy, t_xz_zz_x_xz, t_xz_zz_x_yy,\
                           t_xz_zz_x_yz, t_xz_zz_x_zz, t_xz_zz_y_xx, t_xz_zz_y_xy, t_xz_zz_y_xz,\
                           t_xz_zz_y_yy, t_xz_zz_y_yz, t_xz_zz_y_zz, t_xz_zz_z_xx, t_xz_zz_z_xy,\
                           t_xz_zz_z_xz, t_xz_zz_z_yy, t_xz_zz_z_yz, t_xz_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_zz_z_zz[i] = t_x_zzz_z_zz[i] - rab_z[i] * t_x_zz_z_zz[i];

        t_xz_zz_z_yz[i] = t_x_zzz_z_yz[i] - rab_z[i] * t_x_zz_z_yz[i];

        t_xz_zz_z_yy[i] = t_x_zzz_z_yy[i] - rab_z[i] * t_x_zz_z_yy[i];

        t_xz_zz_z_xz[i] = t_x_zzz_z_xz[i] - rab_z[i] * t_x_zz_z_xz[i];

        t_xz_zz_z_xy[i] = t_x_zzz_z_xy[i] - rab_z[i] * t_x_zz_z_xy[i];

        t_xz_zz_z_xx[i] = t_x_zzz_z_xx[i] - rab_z[i] * t_x_zz_z_xx[i];

        t_xz_zz_y_zz[i] = t_x_zzz_y_zz[i] - rab_z[i] * t_x_zz_y_zz[i];

        t_xz_zz_y_yz[i] = t_x_zzz_y_yz[i] - rab_z[i] * t_x_zz_y_yz[i];

        t_xz_zz_y_yy[i] = t_x_zzz_y_yy[i] - rab_z[i] * t_x_zz_y_yy[i];

        t_xz_zz_y_xz[i] = t_x_zzz_y_xz[i] - rab_z[i] * t_x_zz_y_xz[i];

        t_xz_zz_y_xy[i] = t_x_zzz_y_xy[i] - rab_z[i] * t_x_zz_y_xy[i];

        t_xz_zz_y_xx[i] = t_x_zzz_y_xx[i] - rab_z[i] * t_x_zz_y_xx[i];

        t_xz_zz_x_zz[i] = t_x_zzz_x_zz[i] - rab_z[i] * t_x_zz_x_zz[i];

        t_xz_zz_x_yz[i] = t_x_zzz_x_yz[i] - rab_z[i] * t_x_zz_x_yz[i];

        t_xz_zz_x_yy[i] = t_x_zzz_x_yy[i] - rab_z[i] * t_x_zz_x_yy[i];

        t_xz_zz_x_xz[i] = t_x_zzz_x_xz[i] - rab_z[i] * t_x_zz_x_xz[i];

        t_xz_zz_x_xy[i] = t_x_zzz_x_xy[i] - rab_z[i] * t_x_zz_x_xy[i];

        t_xz_zz_x_xx[i] = t_x_zzz_x_xx[i] - rab_z[i] * t_x_zz_x_xx[i];

        t_xz_yz_z_zz[i] = t_x_yzz_z_zz[i] - rab_z[i] * t_x_yz_z_zz[i];

        t_xz_yz_z_yz[i] = t_x_yzz_z_yz[i] - rab_z[i] * t_x_yz_z_yz[i];

        t_xz_yz_z_yy[i] = t_x_yzz_z_yy[i] - rab_z[i] * t_x_yz_z_yy[i];

        t_xz_yz_z_xz[i] = t_x_yzz_z_xz[i] - rab_z[i] * t_x_yz_z_xz[i];

        t_xz_yz_z_xy[i] = t_x_yzz_z_xy[i] - rab_z[i] * t_x_yz_z_xy[i];

        t_xz_yz_z_xx[i] = t_x_yzz_z_xx[i] - rab_z[i] * t_x_yz_z_xx[i];

        t_xz_yz_y_zz[i] = t_x_yzz_y_zz[i] - rab_z[i] * t_x_yz_y_zz[i];

        t_xz_yz_y_yz[i] = t_x_yzz_y_yz[i] - rab_z[i] * t_x_yz_y_yz[i];

        t_xz_yz_y_yy[i] = t_x_yzz_y_yy[i] - rab_z[i] * t_x_yz_y_yy[i];

        t_xz_yz_y_xz[i] = t_x_yzz_y_xz[i] - rab_z[i] * t_x_yz_y_xz[i];

        t_xz_yz_y_xy[i] = t_x_yzz_y_xy[i] - rab_z[i] * t_x_yz_y_xy[i];

        t_xz_yz_y_xx[i] = t_x_yzz_y_xx[i] - rab_z[i] * t_x_yz_y_xx[i];

        t_xz_yz_x_zz[i] = t_x_yzz_x_zz[i] - rab_z[i] * t_x_yz_x_zz[i];

        t_xz_yz_x_yz[i] = t_x_yzz_x_yz[i] - rab_z[i] * t_x_yz_x_yz[i];

        t_xz_yz_x_yy[i] = t_x_yzz_x_yy[i] - rab_z[i] * t_x_yz_x_yy[i];

        t_xz_yz_x_xz[i] = t_x_yzz_x_xz[i] - rab_z[i] * t_x_yz_x_xz[i];

        t_xz_yz_x_xy[i] = t_x_yzz_x_xy[i] - rab_z[i] * t_x_yz_x_xy[i];

        t_xz_yz_x_xx[i] = t_x_yzz_x_xx[i] - rab_z[i] * t_x_yz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xz_x_xx, t_x_xz_x_xy, t_x_xz_x_xz, t_x_xz_x_yy,\
                           t_x_xz_x_yz, t_x_xz_x_zz, t_x_xz_y_xx, t_x_xz_y_xy, t_x_xz_y_xz,\
                           t_x_xz_y_yy, t_x_xz_y_yz, t_x_xz_y_zz, t_x_xz_z_xx, t_x_xz_z_xy,\
                           t_x_xz_z_xz, t_x_xz_z_yy, t_x_xz_z_yz, t_x_xz_z_zz, t_x_xzz_x_xx,\
                           t_x_xzz_x_xy, t_x_xzz_x_xz, t_x_xzz_x_yy, t_x_xzz_x_yz, t_x_xzz_x_zz,\
                           t_x_xzz_y_xx, t_x_xzz_y_xy, t_x_xzz_y_xz, t_x_xzz_y_yy, t_x_xzz_y_yz,\
                           t_x_xzz_y_zz, t_x_xzz_z_xx, t_x_xzz_z_xy, t_x_xzz_z_xz, t_x_xzz_z_yy,\
                           t_x_xzz_z_yz, t_x_xzz_z_zz, t_x_yy_x_xx, t_x_yy_x_xy, t_x_yy_x_xz,\
                           t_x_yy_x_yy, t_x_yy_x_yz, t_x_yy_x_zz, t_x_yy_y_xx, t_x_yy_y_xy,\
                           t_x_yy_y_xz, t_x_yy_y_yy, t_x_yy_y_yz, t_x_yy_y_zz, t_x_yy_z_xx,\
                           t_x_yy_z_xy, t_x_yy_z_xz, t_x_yy_z_yy, t_x_yy_z_yz, t_x_yy_z_zz,\
                           t_x_yyz_x_xx, t_x_yyz_x_xy, t_x_yyz_x_xz, t_x_yyz_x_yy, t_x_yyz_x_yz,\
                           t_x_yyz_x_zz, t_x_yyz_y_xx, t_x_yyz_y_xy, t_x_yyz_y_xz, t_x_yyz_y_yy,\
                           t_x_yyz_y_yz, t_x_yyz_y_zz, t_x_yyz_z_xx, t_x_yyz_z_xy, t_x_yyz_z_xz,\
                           t_x_yyz_z_yy, t_x_yyz_z_yz, t_x_yyz_z_zz, t_xz_xz_x_xx, t_xz_xz_x_xy,\
                           t_xz_xz_x_xz, t_xz_xz_x_yy, t_xz_xz_x_yz, t_xz_xz_x_zz, t_xz_xz_y_xx,\
                           t_xz_xz_y_xy, t_xz_xz_y_xz, t_xz_xz_y_yy, t_xz_xz_y_yz, t_xz_xz_y_zz,\
                           t_xz_xz_z_xx, t_xz_xz_z_xy, t_xz_xz_z_xz, t_xz_xz_z_yy, t_xz_xz_z_yz,\
                           t_xz_xz_z_zz, t_xz_yy_x_xx, t_xz_yy_x_xy, t_xz_yy_x_xz, t_xz_yy_x_yy,\
                           t_xz_yy_x_yz, t_xz_yy_x_zz, t_xz_yy_y_xx, t_xz_yy_y_xy, t_xz_yy_y_xz,\
                           t_xz_yy_y_yy, t_xz_yy_y_yz, t_xz_yy_y_zz, t_xz_yy_z_xx, t_xz_yy_z_xy,\
                           t_xz_yy_z_xz, t_xz_yy_z_yy, t_xz_yy_z_yz, t_xz_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_yy_z_zz[i] = t_x_yyz_z_zz[i] - rab_z[i] * t_x_yy_z_zz[i];

        t_xz_yy_z_yz[i] = t_x_yyz_z_yz[i] - rab_z[i] * t_x_yy_z_yz[i];

        t_xz_yy_z_yy[i] = t_x_yyz_z_yy[i] - rab_z[i] * t_x_yy_z_yy[i];

        t_xz_yy_z_xz[i] = t_x_yyz_z_xz[i] - rab_z[i] * t_x_yy_z_xz[i];

        t_xz_yy_z_xy[i] = t_x_yyz_z_xy[i] - rab_z[i] * t_x_yy_z_xy[i];

        t_xz_yy_z_xx[i] = t_x_yyz_z_xx[i] - rab_z[i] * t_x_yy_z_xx[i];

        t_xz_yy_y_zz[i] = t_x_yyz_y_zz[i] - rab_z[i] * t_x_yy_y_zz[i];

        t_xz_yy_y_yz[i] = t_x_yyz_y_yz[i] - rab_z[i] * t_x_yy_y_yz[i];

        t_xz_yy_y_yy[i] = t_x_yyz_y_yy[i] - rab_z[i] * t_x_yy_y_yy[i];

        t_xz_yy_y_xz[i] = t_x_yyz_y_xz[i] - rab_z[i] * t_x_yy_y_xz[i];

        t_xz_yy_y_xy[i] = t_x_yyz_y_xy[i] - rab_z[i] * t_x_yy_y_xy[i];

        t_xz_yy_y_xx[i] = t_x_yyz_y_xx[i] - rab_z[i] * t_x_yy_y_xx[i];

        t_xz_yy_x_zz[i] = t_x_yyz_x_zz[i] - rab_z[i] * t_x_yy_x_zz[i];

        t_xz_yy_x_yz[i] = t_x_yyz_x_yz[i] - rab_z[i] * t_x_yy_x_yz[i];

        t_xz_yy_x_yy[i] = t_x_yyz_x_yy[i] - rab_z[i] * t_x_yy_x_yy[i];

        t_xz_yy_x_xz[i] = t_x_yyz_x_xz[i] - rab_z[i] * t_x_yy_x_xz[i];

        t_xz_yy_x_xy[i] = t_x_yyz_x_xy[i] - rab_z[i] * t_x_yy_x_xy[i];

        t_xz_yy_x_xx[i] = t_x_yyz_x_xx[i] - rab_z[i] * t_x_yy_x_xx[i];

        t_xz_xz_z_zz[i] = t_x_xzz_z_zz[i] - rab_z[i] * t_x_xz_z_zz[i];

        t_xz_xz_z_yz[i] = t_x_xzz_z_yz[i] - rab_z[i] * t_x_xz_z_yz[i];

        t_xz_xz_z_yy[i] = t_x_xzz_z_yy[i] - rab_z[i] * t_x_xz_z_yy[i];

        t_xz_xz_z_xz[i] = t_x_xzz_z_xz[i] - rab_z[i] * t_x_xz_z_xz[i];

        t_xz_xz_z_xy[i] = t_x_xzz_z_xy[i] - rab_z[i] * t_x_xz_z_xy[i];

        t_xz_xz_z_xx[i] = t_x_xzz_z_xx[i] - rab_z[i] * t_x_xz_z_xx[i];

        t_xz_xz_y_zz[i] = t_x_xzz_y_zz[i] - rab_z[i] * t_x_xz_y_zz[i];

        t_xz_xz_y_yz[i] = t_x_xzz_y_yz[i] - rab_z[i] * t_x_xz_y_yz[i];

        t_xz_xz_y_yy[i] = t_x_xzz_y_yy[i] - rab_z[i] * t_x_xz_y_yy[i];

        t_xz_xz_y_xz[i] = t_x_xzz_y_xz[i] - rab_z[i] * t_x_xz_y_xz[i];

        t_xz_xz_y_xy[i] = t_x_xzz_y_xy[i] - rab_z[i] * t_x_xz_y_xy[i];

        t_xz_xz_y_xx[i] = t_x_xzz_y_xx[i] - rab_z[i] * t_x_xz_y_xx[i];

        t_xz_xz_x_zz[i] = t_x_xzz_x_zz[i] - rab_z[i] * t_x_xz_x_zz[i];

        t_xz_xz_x_yz[i] = t_x_xzz_x_yz[i] - rab_z[i] * t_x_xz_x_yz[i];

        t_xz_xz_x_yy[i] = t_x_xzz_x_yy[i] - rab_z[i] * t_x_xz_x_yy[i];

        t_xz_xz_x_xz[i] = t_x_xzz_x_xz[i] - rab_z[i] * t_x_xz_x_xz[i];

        t_xz_xz_x_xy[i] = t_x_xzz_x_xy[i] - rab_z[i] * t_x_xz_x_xy[i];

        t_xz_xz_x_xx[i] = t_x_xzz_x_xx[i] - rab_z[i] * t_x_xz_x_xx[i];
    }

    #pragma omp simd align(rab_z, t_x_xx_x_xx, t_x_xx_x_xy, t_x_xx_x_xz, t_x_xx_x_yy,\
                           t_x_xx_x_yz, t_x_xx_x_zz, t_x_xx_y_xx, t_x_xx_y_xy, t_x_xx_y_xz,\
                           t_x_xx_y_yy, t_x_xx_y_yz, t_x_xx_y_zz, t_x_xx_z_xx, t_x_xx_z_xy,\
                           t_x_xx_z_xz, t_x_xx_z_yy, t_x_xx_z_yz, t_x_xx_z_zz, t_x_xxz_x_xx,\
                           t_x_xxz_x_xy, t_x_xxz_x_xz, t_x_xxz_x_yy, t_x_xxz_x_yz, t_x_xxz_x_zz,\
                           t_x_xxz_y_xx, t_x_xxz_y_xy, t_x_xxz_y_xz, t_x_xxz_y_yy, t_x_xxz_y_yz,\
                           t_x_xxz_y_zz, t_x_xxz_z_xx, t_x_xxz_z_xy, t_x_xxz_z_xz, t_x_xxz_z_yy,\
                           t_x_xxz_z_yz, t_x_xxz_z_zz, t_x_xy_x_xx, t_x_xy_x_xy, t_x_xy_x_xz,\
                           t_x_xy_x_yy, t_x_xy_x_yz, t_x_xy_x_zz, t_x_xy_y_xx, t_x_xy_y_xy,\
                           t_x_xy_y_xz, t_x_xy_y_yy, t_x_xy_y_yz, t_x_xy_y_zz, t_x_xy_z_xx,\
                           t_x_xy_z_xy, t_x_xy_z_xz, t_x_xy_z_yy, t_x_xy_z_yz, t_x_xy_z_zz,\
                           t_x_xyz_x_xx, t_x_xyz_x_xy, t_x_xyz_x_xz, t_x_xyz_x_yy, t_x_xyz_x_yz,\
                           t_x_xyz_x_zz, t_x_xyz_y_xx, t_x_xyz_y_xy, t_x_xyz_y_xz, t_x_xyz_y_yy,\
                           t_x_xyz_y_yz, t_x_xyz_y_zz, t_x_xyz_z_xx, t_x_xyz_z_xy, t_x_xyz_z_xz,\
                           t_x_xyz_z_yy, t_x_xyz_z_yz, t_x_xyz_z_zz, t_xz_xx_x_xx, t_xz_xx_x_xy,\
                           t_xz_xx_x_xz, t_xz_xx_x_yy, t_xz_xx_x_yz, t_xz_xx_x_zz, t_xz_xx_y_xx,\
                           t_xz_xx_y_xy, t_xz_xx_y_xz, t_xz_xx_y_yy, t_xz_xx_y_yz, t_xz_xx_y_zz,\
                           t_xz_xx_z_xx, t_xz_xx_z_xy, t_xz_xx_z_xz, t_xz_xx_z_yy, t_xz_xx_z_yz,\
                           t_xz_xx_z_zz, t_xz_xy_x_xx, t_xz_xy_x_xy, t_xz_xy_x_xz, t_xz_xy_x_yy,\
                           t_xz_xy_x_yz, t_xz_xy_x_zz, t_xz_xy_y_xx, t_xz_xy_y_xy, t_xz_xy_y_xz,\
                           t_xz_xy_y_yy, t_xz_xy_y_yz, t_xz_xy_y_zz, t_xz_xy_z_xx, t_xz_xy_z_xy,\
                           t_xz_xy_z_xz, t_xz_xy_z_yy, t_xz_xy_z_yz, t_xz_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_xy_z_zz[i] = t_x_xyz_z_zz[i] - rab_z[i] * t_x_xy_z_zz[i];

        t_xz_xy_z_yz[i] = t_x_xyz_z_yz[i] - rab_z[i] * t_x_xy_z_yz[i];

        t_xz_xy_z_yy[i] = t_x_xyz_z_yy[i] - rab_z[i] * t_x_xy_z_yy[i];

        t_xz_xy_z_xz[i] = t_x_xyz_z_xz[i] - rab_z[i] * t_x_xy_z_xz[i];

        t_xz_xy_z_xy[i] = t_x_xyz_z_xy[i] - rab_z[i] * t_x_xy_z_xy[i];

        t_xz_xy_z_xx[i] = t_x_xyz_z_xx[i] - rab_z[i] * t_x_xy_z_xx[i];

        t_xz_xy_y_zz[i] = t_x_xyz_y_zz[i] - rab_z[i] * t_x_xy_y_zz[i];

        t_xz_xy_y_yz[i] = t_x_xyz_y_yz[i] - rab_z[i] * t_x_xy_y_yz[i];

        t_xz_xy_y_yy[i] = t_x_xyz_y_yy[i] - rab_z[i] * t_x_xy_y_yy[i];

        t_xz_xy_y_xz[i] = t_x_xyz_y_xz[i] - rab_z[i] * t_x_xy_y_xz[i];

        t_xz_xy_y_xy[i] = t_x_xyz_y_xy[i] - rab_z[i] * t_x_xy_y_xy[i];

        t_xz_xy_y_xx[i] = t_x_xyz_y_xx[i] - rab_z[i] * t_x_xy_y_xx[i];

        t_xz_xy_x_zz[i] = t_x_xyz_x_zz[i] - rab_z[i] * t_x_xy_x_zz[i];

        t_xz_xy_x_yz[i] = t_x_xyz_x_yz[i] - rab_z[i] * t_x_xy_x_yz[i];

        t_xz_xy_x_yy[i] = t_x_xyz_x_yy[i] - rab_z[i] * t_x_xy_x_yy[i];

        t_xz_xy_x_xz[i] = t_x_xyz_x_xz[i] - rab_z[i] * t_x_xy_x_xz[i];

        t_xz_xy_x_xy[i] = t_x_xyz_x_xy[i] - rab_z[i] * t_x_xy_x_xy[i];

        t_xz_xy_x_xx[i] = t_x_xyz_x_xx[i] - rab_z[i] * t_x_xy_x_xx[i];

        t_xz_xx_z_zz[i] = t_x_xxz_z_zz[i] - rab_z[i] * t_x_xx_z_zz[i];

        t_xz_xx_z_yz[i] = t_x_xxz_z_yz[i] - rab_z[i] * t_x_xx_z_yz[i];

        t_xz_xx_z_yy[i] = t_x_xxz_z_yy[i] - rab_z[i] * t_x_xx_z_yy[i];

        t_xz_xx_z_xz[i] = t_x_xxz_z_xz[i] - rab_z[i] * t_x_xx_z_xz[i];

        t_xz_xx_z_xy[i] = t_x_xxz_z_xy[i] - rab_z[i] * t_x_xx_z_xy[i];

        t_xz_xx_z_xx[i] = t_x_xxz_z_xx[i] - rab_z[i] * t_x_xx_z_xx[i];

        t_xz_xx_y_zz[i] = t_x_xxz_y_zz[i] - rab_z[i] * t_x_xx_y_zz[i];

        t_xz_xx_y_yz[i] = t_x_xxz_y_yz[i] - rab_z[i] * t_x_xx_y_yz[i];

        t_xz_xx_y_yy[i] = t_x_xxz_y_yy[i] - rab_z[i] * t_x_xx_y_yy[i];

        t_xz_xx_y_xz[i] = t_x_xxz_y_xz[i] - rab_z[i] * t_x_xx_y_xz[i];

        t_xz_xx_y_xy[i] = t_x_xxz_y_xy[i] - rab_z[i] * t_x_xx_y_xy[i];

        t_xz_xx_y_xx[i] = t_x_xxz_y_xx[i] - rab_z[i] * t_x_xx_y_xx[i];

        t_xz_xx_x_zz[i] = t_x_xxz_x_zz[i] - rab_z[i] * t_x_xx_x_zz[i];

        t_xz_xx_x_yz[i] = t_x_xxz_x_yz[i] - rab_z[i] * t_x_xx_x_yz[i];

        t_xz_xx_x_yy[i] = t_x_xxz_x_yy[i] - rab_z[i] * t_x_xx_x_yy[i];

        t_xz_xx_x_xz[i] = t_x_xxz_x_xz[i] - rab_z[i] * t_x_xx_x_xz[i];

        t_xz_xx_x_xy[i] = t_x_xxz_x_xy[i] - rab_z[i] * t_x_xx_x_xy[i];

        t_xz_xx_x_xx[i] = t_x_xxz_x_xx[i] - rab_z[i] * t_x_xx_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_yyz_x_xx, t_x_yyz_x_xy, t_x_yyz_x_xz, t_x_yyz_x_yy,\
                           t_x_yyz_x_yz, t_x_yyz_x_zz, t_x_yyz_y_xx, t_x_yyz_y_xy, t_x_yyz_y_xz,\
                           t_x_yyz_y_yy, t_x_yyz_y_yz, t_x_yyz_y_zz, t_x_yyz_z_xx, t_x_yyz_z_xy,\
                           t_x_yyz_z_xz, t_x_yyz_z_yy, t_x_yyz_z_yz, t_x_yyz_z_zz, t_x_yz_x_xx,\
                           t_x_yz_x_xy, t_x_yz_x_xz, t_x_yz_x_yy, t_x_yz_x_yz, t_x_yz_x_zz,\
                           t_x_yz_y_xx, t_x_yz_y_xy, t_x_yz_y_xz, t_x_yz_y_yy, t_x_yz_y_yz,\
                           t_x_yz_y_zz, t_x_yz_z_xx, t_x_yz_z_xy, t_x_yz_z_xz, t_x_yz_z_yy,\
                           t_x_yz_z_yz, t_x_yz_z_zz, t_x_yzz_x_xx, t_x_yzz_x_xy, t_x_yzz_x_xz,\
                           t_x_yzz_x_yy, t_x_yzz_x_yz, t_x_yzz_x_zz, t_x_yzz_y_xx, t_x_yzz_y_xy,\
                           t_x_yzz_y_xz, t_x_yzz_y_yy, t_x_yzz_y_yz, t_x_yzz_y_zz, t_x_yzz_z_xx,\
                           t_x_yzz_z_xy, t_x_yzz_z_xz, t_x_yzz_z_yy, t_x_yzz_z_yz, t_x_yzz_z_zz,\
                           t_x_zz_x_xx, t_x_zz_x_xy, t_x_zz_x_xz, t_x_zz_x_yy, t_x_zz_x_yz,\
                           t_x_zz_x_zz, t_x_zz_y_xx, t_x_zz_y_xy, t_x_zz_y_xz, t_x_zz_y_yy,\
                           t_x_zz_y_yz, t_x_zz_y_zz, t_x_zz_z_xx, t_x_zz_z_xy, t_x_zz_z_xz,\
                           t_x_zz_z_yy, t_x_zz_z_yz, t_x_zz_z_zz, t_xy_yz_x_xx, t_xy_yz_x_xy,\
                           t_xy_yz_x_xz, t_xy_yz_x_yy, t_xy_yz_x_yz, t_xy_yz_x_zz, t_xy_yz_y_xx,\
                           t_xy_yz_y_xy, t_xy_yz_y_xz, t_xy_yz_y_yy, t_xy_yz_y_yz, t_xy_yz_y_zz,\
                           t_xy_yz_z_xx, t_xy_yz_z_xy, t_xy_yz_z_xz, t_xy_yz_z_yy, t_xy_yz_z_yz,\
                           t_xy_yz_z_zz, t_xy_zz_x_xx, t_xy_zz_x_xy, t_xy_zz_x_xz, t_xy_zz_x_yy,\
                           t_xy_zz_x_yz, t_xy_zz_x_zz, t_xy_zz_y_xx, t_xy_zz_y_xy, t_xy_zz_y_xz,\
                           t_xy_zz_y_yy, t_xy_zz_y_yz, t_xy_zz_y_zz, t_xy_zz_z_xx, t_xy_zz_z_xy,\
                           t_xy_zz_z_xz, t_xy_zz_z_yy, t_xy_zz_z_yz, t_xy_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_zz_z_zz[i] = t_x_yzz_z_zz[i] - rab_y[i] * t_x_zz_z_zz[i];

        t_xy_zz_z_yz[i] = t_x_yzz_z_yz[i] - rab_y[i] * t_x_zz_z_yz[i];

        t_xy_zz_z_yy[i] = t_x_yzz_z_yy[i] - rab_y[i] * t_x_zz_z_yy[i];

        t_xy_zz_z_xz[i] = t_x_yzz_z_xz[i] - rab_y[i] * t_x_zz_z_xz[i];

        t_xy_zz_z_xy[i] = t_x_yzz_z_xy[i] - rab_y[i] * t_x_zz_z_xy[i];

        t_xy_zz_z_xx[i] = t_x_yzz_z_xx[i] - rab_y[i] * t_x_zz_z_xx[i];

        t_xy_zz_y_zz[i] = t_x_yzz_y_zz[i] - rab_y[i] * t_x_zz_y_zz[i];

        t_xy_zz_y_yz[i] = t_x_yzz_y_yz[i] - rab_y[i] * t_x_zz_y_yz[i];

        t_xy_zz_y_yy[i] = t_x_yzz_y_yy[i] - rab_y[i] * t_x_zz_y_yy[i];

        t_xy_zz_y_xz[i] = t_x_yzz_y_xz[i] - rab_y[i] * t_x_zz_y_xz[i];

        t_xy_zz_y_xy[i] = t_x_yzz_y_xy[i] - rab_y[i] * t_x_zz_y_xy[i];

        t_xy_zz_y_xx[i] = t_x_yzz_y_xx[i] - rab_y[i] * t_x_zz_y_xx[i];

        t_xy_zz_x_zz[i] = t_x_yzz_x_zz[i] - rab_y[i] * t_x_zz_x_zz[i];

        t_xy_zz_x_yz[i] = t_x_yzz_x_yz[i] - rab_y[i] * t_x_zz_x_yz[i];

        t_xy_zz_x_yy[i] = t_x_yzz_x_yy[i] - rab_y[i] * t_x_zz_x_yy[i];

        t_xy_zz_x_xz[i] = t_x_yzz_x_xz[i] - rab_y[i] * t_x_zz_x_xz[i];

        t_xy_zz_x_xy[i] = t_x_yzz_x_xy[i] - rab_y[i] * t_x_zz_x_xy[i];

        t_xy_zz_x_xx[i] = t_x_yzz_x_xx[i] - rab_y[i] * t_x_zz_x_xx[i];

        t_xy_yz_z_zz[i] = t_x_yyz_z_zz[i] - rab_y[i] * t_x_yz_z_zz[i];

        t_xy_yz_z_yz[i] = t_x_yyz_z_yz[i] - rab_y[i] * t_x_yz_z_yz[i];

        t_xy_yz_z_yy[i] = t_x_yyz_z_yy[i] - rab_y[i] * t_x_yz_z_yy[i];

        t_xy_yz_z_xz[i] = t_x_yyz_z_xz[i] - rab_y[i] * t_x_yz_z_xz[i];

        t_xy_yz_z_xy[i] = t_x_yyz_z_xy[i] - rab_y[i] * t_x_yz_z_xy[i];

        t_xy_yz_z_xx[i] = t_x_yyz_z_xx[i] - rab_y[i] * t_x_yz_z_xx[i];

        t_xy_yz_y_zz[i] = t_x_yyz_y_zz[i] - rab_y[i] * t_x_yz_y_zz[i];

        t_xy_yz_y_yz[i] = t_x_yyz_y_yz[i] - rab_y[i] * t_x_yz_y_yz[i];

        t_xy_yz_y_yy[i] = t_x_yyz_y_yy[i] - rab_y[i] * t_x_yz_y_yy[i];

        t_xy_yz_y_xz[i] = t_x_yyz_y_xz[i] - rab_y[i] * t_x_yz_y_xz[i];

        t_xy_yz_y_xy[i] = t_x_yyz_y_xy[i] - rab_y[i] * t_x_yz_y_xy[i];

        t_xy_yz_y_xx[i] = t_x_yyz_y_xx[i] - rab_y[i] * t_x_yz_y_xx[i];

        t_xy_yz_x_zz[i] = t_x_yyz_x_zz[i] - rab_y[i] * t_x_yz_x_zz[i];

        t_xy_yz_x_yz[i] = t_x_yyz_x_yz[i] - rab_y[i] * t_x_yz_x_yz[i];

        t_xy_yz_x_yy[i] = t_x_yyz_x_yy[i] - rab_y[i] * t_x_yz_x_yy[i];

        t_xy_yz_x_xz[i] = t_x_yyz_x_xz[i] - rab_y[i] * t_x_yz_x_xz[i];

        t_xy_yz_x_xy[i] = t_x_yyz_x_xy[i] - rab_y[i] * t_x_yz_x_xy[i];

        t_xy_yz_x_xx[i] = t_x_yyz_x_xx[i] - rab_y[i] * t_x_yz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xyz_x_xx, t_x_xyz_x_xy, t_x_xyz_x_xz, t_x_xyz_x_yy,\
                           t_x_xyz_x_yz, t_x_xyz_x_zz, t_x_xyz_y_xx, t_x_xyz_y_xy, t_x_xyz_y_xz,\
                           t_x_xyz_y_yy, t_x_xyz_y_yz, t_x_xyz_y_zz, t_x_xyz_z_xx, t_x_xyz_z_xy,\
                           t_x_xyz_z_xz, t_x_xyz_z_yy, t_x_xyz_z_yz, t_x_xyz_z_zz, t_x_xz_x_xx,\
                           t_x_xz_x_xy, t_x_xz_x_xz, t_x_xz_x_yy, t_x_xz_x_yz, t_x_xz_x_zz,\
                           t_x_xz_y_xx, t_x_xz_y_xy, t_x_xz_y_xz, t_x_xz_y_yy, t_x_xz_y_yz,\
                           t_x_xz_y_zz, t_x_xz_z_xx, t_x_xz_z_xy, t_x_xz_z_xz, t_x_xz_z_yy,\
                           t_x_xz_z_yz, t_x_xz_z_zz, t_x_yy_x_xx, t_x_yy_x_xy, t_x_yy_x_xz,\
                           t_x_yy_x_yy, t_x_yy_x_yz, t_x_yy_x_zz, t_x_yy_y_xx, t_x_yy_y_xy,\
                           t_x_yy_y_xz, t_x_yy_y_yy, t_x_yy_y_yz, t_x_yy_y_zz, t_x_yy_z_xx,\
                           t_x_yy_z_xy, t_x_yy_z_xz, t_x_yy_z_yy, t_x_yy_z_yz, t_x_yy_z_zz,\
                           t_x_yyy_x_xx, t_x_yyy_x_xy, t_x_yyy_x_xz, t_x_yyy_x_yy, t_x_yyy_x_yz,\
                           t_x_yyy_x_zz, t_x_yyy_y_xx, t_x_yyy_y_xy, t_x_yyy_y_xz, t_x_yyy_y_yy,\
                           t_x_yyy_y_yz, t_x_yyy_y_zz, t_x_yyy_z_xx, t_x_yyy_z_xy, t_x_yyy_z_xz,\
                           t_x_yyy_z_yy, t_x_yyy_z_yz, t_x_yyy_z_zz, t_xy_xz_x_xx, t_xy_xz_x_xy,\
                           t_xy_xz_x_xz, t_xy_xz_x_yy, t_xy_xz_x_yz, t_xy_xz_x_zz, t_xy_xz_y_xx,\
                           t_xy_xz_y_xy, t_xy_xz_y_xz, t_xy_xz_y_yy, t_xy_xz_y_yz, t_xy_xz_y_zz,\
                           t_xy_xz_z_xx, t_xy_xz_z_xy, t_xy_xz_z_xz, t_xy_xz_z_yy, t_xy_xz_z_yz,\
                           t_xy_xz_z_zz, t_xy_yy_x_xx, t_xy_yy_x_xy, t_xy_yy_x_xz, t_xy_yy_x_yy,\
                           t_xy_yy_x_yz, t_xy_yy_x_zz, t_xy_yy_y_xx, t_xy_yy_y_xy, t_xy_yy_y_xz,\
                           t_xy_yy_y_yy, t_xy_yy_y_yz, t_xy_yy_y_zz, t_xy_yy_z_xx, t_xy_yy_z_xy,\
                           t_xy_yy_z_xz, t_xy_yy_z_yy, t_xy_yy_z_yz, t_xy_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_yy_z_zz[i] = t_x_yyy_z_zz[i] - rab_y[i] * t_x_yy_z_zz[i];

        t_xy_yy_z_yz[i] = t_x_yyy_z_yz[i] - rab_y[i] * t_x_yy_z_yz[i];

        t_xy_yy_z_yy[i] = t_x_yyy_z_yy[i] - rab_y[i] * t_x_yy_z_yy[i];

        t_xy_yy_z_xz[i] = t_x_yyy_z_xz[i] - rab_y[i] * t_x_yy_z_xz[i];

        t_xy_yy_z_xy[i] = t_x_yyy_z_xy[i] - rab_y[i] * t_x_yy_z_xy[i];

        t_xy_yy_z_xx[i] = t_x_yyy_z_xx[i] - rab_y[i] * t_x_yy_z_xx[i];

        t_xy_yy_y_zz[i] = t_x_yyy_y_zz[i] - rab_y[i] * t_x_yy_y_zz[i];

        t_xy_yy_y_yz[i] = t_x_yyy_y_yz[i] - rab_y[i] * t_x_yy_y_yz[i];

        t_xy_yy_y_yy[i] = t_x_yyy_y_yy[i] - rab_y[i] * t_x_yy_y_yy[i];

        t_xy_yy_y_xz[i] = t_x_yyy_y_xz[i] - rab_y[i] * t_x_yy_y_xz[i];

        t_xy_yy_y_xy[i] = t_x_yyy_y_xy[i] - rab_y[i] * t_x_yy_y_xy[i];

        t_xy_yy_y_xx[i] = t_x_yyy_y_xx[i] - rab_y[i] * t_x_yy_y_xx[i];

        t_xy_yy_x_zz[i] = t_x_yyy_x_zz[i] - rab_y[i] * t_x_yy_x_zz[i];

        t_xy_yy_x_yz[i] = t_x_yyy_x_yz[i] - rab_y[i] * t_x_yy_x_yz[i];

        t_xy_yy_x_yy[i] = t_x_yyy_x_yy[i] - rab_y[i] * t_x_yy_x_yy[i];

        t_xy_yy_x_xz[i] = t_x_yyy_x_xz[i] - rab_y[i] * t_x_yy_x_xz[i];

        t_xy_yy_x_xy[i] = t_x_yyy_x_xy[i] - rab_y[i] * t_x_yy_x_xy[i];

        t_xy_yy_x_xx[i] = t_x_yyy_x_xx[i] - rab_y[i] * t_x_yy_x_xx[i];

        t_xy_xz_z_zz[i] = t_x_xyz_z_zz[i] - rab_y[i] * t_x_xz_z_zz[i];

        t_xy_xz_z_yz[i] = t_x_xyz_z_yz[i] - rab_y[i] * t_x_xz_z_yz[i];

        t_xy_xz_z_yy[i] = t_x_xyz_z_yy[i] - rab_y[i] * t_x_xz_z_yy[i];

        t_xy_xz_z_xz[i] = t_x_xyz_z_xz[i] - rab_y[i] * t_x_xz_z_xz[i];

        t_xy_xz_z_xy[i] = t_x_xyz_z_xy[i] - rab_y[i] * t_x_xz_z_xy[i];

        t_xy_xz_z_xx[i] = t_x_xyz_z_xx[i] - rab_y[i] * t_x_xz_z_xx[i];

        t_xy_xz_y_zz[i] = t_x_xyz_y_zz[i] - rab_y[i] * t_x_xz_y_zz[i];

        t_xy_xz_y_yz[i] = t_x_xyz_y_yz[i] - rab_y[i] * t_x_xz_y_yz[i];

        t_xy_xz_y_yy[i] = t_x_xyz_y_yy[i] - rab_y[i] * t_x_xz_y_yy[i];

        t_xy_xz_y_xz[i] = t_x_xyz_y_xz[i] - rab_y[i] * t_x_xz_y_xz[i];

        t_xy_xz_y_xy[i] = t_x_xyz_y_xy[i] - rab_y[i] * t_x_xz_y_xy[i];

        t_xy_xz_y_xx[i] = t_x_xyz_y_xx[i] - rab_y[i] * t_x_xz_y_xx[i];

        t_xy_xz_x_zz[i] = t_x_xyz_x_zz[i] - rab_y[i] * t_x_xz_x_zz[i];

        t_xy_xz_x_yz[i] = t_x_xyz_x_yz[i] - rab_y[i] * t_x_xz_x_yz[i];

        t_xy_xz_x_yy[i] = t_x_xyz_x_yy[i] - rab_y[i] * t_x_xz_x_yy[i];

        t_xy_xz_x_xz[i] = t_x_xyz_x_xz[i] - rab_y[i] * t_x_xz_x_xz[i];

        t_xy_xz_x_xy[i] = t_x_xyz_x_xy[i] - rab_y[i] * t_x_xz_x_xy[i];

        t_xy_xz_x_xx[i] = t_x_xyz_x_xx[i] - rab_y[i] * t_x_xz_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_x_xx_x_xx, t_x_xx_x_xy, t_x_xx_x_xz, t_x_xx_x_yy,\
                           t_x_xx_x_yz, t_x_xx_x_zz, t_x_xx_y_xx, t_x_xx_y_xy, t_x_xx_y_xz,\
                           t_x_xx_y_yy, t_x_xx_y_yz, t_x_xx_y_zz, t_x_xx_z_xx, t_x_xx_z_xy,\
                           t_x_xx_z_xz, t_x_xx_z_yy, t_x_xx_z_yz, t_x_xx_z_zz, t_x_xxy_x_xx,\
                           t_x_xxy_x_xy, t_x_xxy_x_xz, t_x_xxy_x_yy, t_x_xxy_x_yz, t_x_xxy_x_zz,\
                           t_x_xxy_y_xx, t_x_xxy_y_xy, t_x_xxy_y_xz, t_x_xxy_y_yy, t_x_xxy_y_yz,\
                           t_x_xxy_y_zz, t_x_xxy_z_xx, t_x_xxy_z_xy, t_x_xxy_z_xz, t_x_xxy_z_yy,\
                           t_x_xxy_z_yz, t_x_xxy_z_zz, t_x_xy_x_xx, t_x_xy_x_xy, t_x_xy_x_xz,\
                           t_x_xy_x_yy, t_x_xy_x_yz, t_x_xy_x_zz, t_x_xy_y_xx, t_x_xy_y_xy,\
                           t_x_xy_y_xz, t_x_xy_y_yy, t_x_xy_y_yz, t_x_xy_y_zz, t_x_xy_z_xx,\
                           t_x_xy_z_xy, t_x_xy_z_xz, t_x_xy_z_yy, t_x_xy_z_yz, t_x_xy_z_zz,\
                           t_x_xyy_x_xx, t_x_xyy_x_xy, t_x_xyy_x_xz, t_x_xyy_x_yy, t_x_xyy_x_yz,\
                           t_x_xyy_x_zz, t_x_xyy_y_xx, t_x_xyy_y_xy, t_x_xyy_y_xz, t_x_xyy_y_yy,\
                           t_x_xyy_y_yz, t_x_xyy_y_zz, t_x_xyy_z_xx, t_x_xyy_z_xy, t_x_xyy_z_xz,\
                           t_x_xyy_z_yy, t_x_xyy_z_yz, t_x_xyy_z_zz, t_xy_xx_x_xx, t_xy_xx_x_xy,\
                           t_xy_xx_x_xz, t_xy_xx_x_yy, t_xy_xx_x_yz, t_xy_xx_x_zz, t_xy_xx_y_xx,\
                           t_xy_xx_y_xy, t_xy_xx_y_xz, t_xy_xx_y_yy, t_xy_xx_y_yz, t_xy_xx_y_zz,\
                           t_xy_xx_z_xx, t_xy_xx_z_xy, t_xy_xx_z_xz, t_xy_xx_z_yy, t_xy_xx_z_yz,\
                           t_xy_xx_z_zz, t_xy_xy_x_xx, t_xy_xy_x_xy, t_xy_xy_x_xz, t_xy_xy_x_yy,\
                           t_xy_xy_x_yz, t_xy_xy_x_zz, t_xy_xy_y_xx, t_xy_xy_y_xy, t_xy_xy_y_xz,\
                           t_xy_xy_y_yy, t_xy_xy_y_yz, t_xy_xy_y_zz, t_xy_xy_z_xx, t_xy_xy_z_xy,\
                           t_xy_xy_z_xz, t_xy_xy_z_yy, t_xy_xy_z_yz, t_xy_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_xy_z_zz[i] = t_x_xyy_z_zz[i] - rab_y[i] * t_x_xy_z_zz[i];

        t_xy_xy_z_yz[i] = t_x_xyy_z_yz[i] - rab_y[i] * t_x_xy_z_yz[i];

        t_xy_xy_z_yy[i] = t_x_xyy_z_yy[i] - rab_y[i] * t_x_xy_z_yy[i];

        t_xy_xy_z_xz[i] = t_x_xyy_z_xz[i] - rab_y[i] * t_x_xy_z_xz[i];

        t_xy_xy_z_xy[i] = t_x_xyy_z_xy[i] - rab_y[i] * t_x_xy_z_xy[i];

        t_xy_xy_z_xx[i] = t_x_xyy_z_xx[i] - rab_y[i] * t_x_xy_z_xx[i];

        t_xy_xy_y_zz[i] = t_x_xyy_y_zz[i] - rab_y[i] * t_x_xy_y_zz[i];

        t_xy_xy_y_yz[i] = t_x_xyy_y_yz[i] - rab_y[i] * t_x_xy_y_yz[i];

        t_xy_xy_y_yy[i] = t_x_xyy_y_yy[i] - rab_y[i] * t_x_xy_y_yy[i];

        t_xy_xy_y_xz[i] = t_x_xyy_y_xz[i] - rab_y[i] * t_x_xy_y_xz[i];

        t_xy_xy_y_xy[i] = t_x_xyy_y_xy[i] - rab_y[i] * t_x_xy_y_xy[i];

        t_xy_xy_y_xx[i] = t_x_xyy_y_xx[i] - rab_y[i] * t_x_xy_y_xx[i];

        t_xy_xy_x_zz[i] = t_x_xyy_x_zz[i] - rab_y[i] * t_x_xy_x_zz[i];

        t_xy_xy_x_yz[i] = t_x_xyy_x_yz[i] - rab_y[i] * t_x_xy_x_yz[i];

        t_xy_xy_x_yy[i] = t_x_xyy_x_yy[i] - rab_y[i] * t_x_xy_x_yy[i];

        t_xy_xy_x_xz[i] = t_x_xyy_x_xz[i] - rab_y[i] * t_x_xy_x_xz[i];

        t_xy_xy_x_xy[i] = t_x_xyy_x_xy[i] - rab_y[i] * t_x_xy_x_xy[i];

        t_xy_xy_x_xx[i] = t_x_xyy_x_xx[i] - rab_y[i] * t_x_xy_x_xx[i];

        t_xy_xx_z_zz[i] = t_x_xxy_z_zz[i] - rab_y[i] * t_x_xx_z_zz[i];

        t_xy_xx_z_yz[i] = t_x_xxy_z_yz[i] - rab_y[i] * t_x_xx_z_yz[i];

        t_xy_xx_z_yy[i] = t_x_xxy_z_yy[i] - rab_y[i] * t_x_xx_z_yy[i];

        t_xy_xx_z_xz[i] = t_x_xxy_z_xz[i] - rab_y[i] * t_x_xx_z_xz[i];

        t_xy_xx_z_xy[i] = t_x_xxy_z_xy[i] - rab_y[i] * t_x_xx_z_xy[i];

        t_xy_xx_z_xx[i] = t_x_xxy_z_xx[i] - rab_y[i] * t_x_xx_z_xx[i];

        t_xy_xx_y_zz[i] = t_x_xxy_y_zz[i] - rab_y[i] * t_x_xx_y_zz[i];

        t_xy_xx_y_yz[i] = t_x_xxy_y_yz[i] - rab_y[i] * t_x_xx_y_yz[i];

        t_xy_xx_y_yy[i] = t_x_xxy_y_yy[i] - rab_y[i] * t_x_xx_y_yy[i];

        t_xy_xx_y_xz[i] = t_x_xxy_y_xz[i] - rab_y[i] * t_x_xx_y_xz[i];

        t_xy_xx_y_xy[i] = t_x_xxy_y_xy[i] - rab_y[i] * t_x_xx_y_xy[i];

        t_xy_xx_y_xx[i] = t_x_xxy_y_xx[i] - rab_y[i] * t_x_xx_y_xx[i];

        t_xy_xx_x_zz[i] = t_x_xxy_x_zz[i] - rab_y[i] * t_x_xx_x_zz[i];

        t_xy_xx_x_yz[i] = t_x_xxy_x_yz[i] - rab_y[i] * t_x_xx_x_yz[i];

        t_xy_xx_x_yy[i] = t_x_xxy_x_yy[i] - rab_y[i] * t_x_xx_x_yy[i];

        t_xy_xx_x_xz[i] = t_x_xxy_x_xz[i] - rab_y[i] * t_x_xx_x_xz[i];

        t_xy_xx_x_xy[i] = t_x_xxy_x_xy[i] - rab_y[i] * t_x_xx_x_xy[i];

        t_xy_xx_x_xx[i] = t_x_xxy_x_xx[i] - rab_y[i] * t_x_xx_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xyz_x_xx, t_x_xyz_x_xy, t_x_xyz_x_xz, t_x_xyz_x_yy,\
                           t_x_xyz_x_yz, t_x_xyz_x_zz, t_x_xyz_y_xx, t_x_xyz_y_xy, t_x_xyz_y_xz,\
                           t_x_xyz_y_yy, t_x_xyz_y_yz, t_x_xyz_y_zz, t_x_xyz_z_xx, t_x_xyz_z_xy,\
                           t_x_xyz_z_xz, t_x_xyz_z_yy, t_x_xyz_z_yz, t_x_xyz_z_zz, t_x_xzz_x_xx,\
                           t_x_xzz_x_xy, t_x_xzz_x_xz, t_x_xzz_x_yy, t_x_xzz_x_yz, t_x_xzz_x_zz,\
                           t_x_xzz_y_xx, t_x_xzz_y_xy, t_x_xzz_y_xz, t_x_xzz_y_yy, t_x_xzz_y_yz,\
                           t_x_xzz_y_zz, t_x_xzz_z_xx, t_x_xzz_z_xy, t_x_xzz_z_xz, t_x_xzz_z_yy,\
                           t_x_xzz_z_yz, t_x_xzz_z_zz, t_x_yz_x_xx, t_x_yz_x_xy, t_x_yz_x_xz,\
                           t_x_yz_x_yy, t_x_yz_x_yz, t_x_yz_x_zz, t_x_yz_y_xx, t_x_yz_y_xy,\
                           t_x_yz_y_xz, t_x_yz_y_yy, t_x_yz_y_yz, t_x_yz_y_zz, t_x_yz_z_xx,\
                           t_x_yz_z_xy, t_x_yz_z_xz, t_x_yz_z_yy, t_x_yz_z_yz, t_x_yz_z_zz,\
                           t_x_zz_x_xx, t_x_zz_x_xy, t_x_zz_x_xz, t_x_zz_x_yy, t_x_zz_x_yz,\
                           t_x_zz_x_zz, t_x_zz_y_xx, t_x_zz_y_xy, t_x_zz_y_xz, t_x_zz_y_yy,\
                           t_x_zz_y_yz, t_x_zz_y_zz, t_x_zz_z_xx, t_x_zz_z_xy, t_x_zz_z_xz,\
                           t_x_zz_z_yy, t_x_zz_z_yz, t_x_zz_z_zz, t_xx_yz_x_xx, t_xx_yz_x_xy,\
                           t_xx_yz_x_xz, t_xx_yz_x_yy, t_xx_yz_x_yz, t_xx_yz_x_zz, t_xx_yz_y_xx,\
                           t_xx_yz_y_xy, t_xx_yz_y_xz, t_xx_yz_y_yy, t_xx_yz_y_yz, t_xx_yz_y_zz,\
                           t_xx_yz_z_xx, t_xx_yz_z_xy, t_xx_yz_z_xz, t_xx_yz_z_yy, t_xx_yz_z_yz,\
                           t_xx_yz_z_zz, t_xx_zz_x_xx, t_xx_zz_x_xy, t_xx_zz_x_xz, t_xx_zz_x_yy,\
                           t_xx_zz_x_yz, t_xx_zz_x_zz, t_xx_zz_y_xx, t_xx_zz_y_xy, t_xx_zz_y_xz,\
                           t_xx_zz_y_yy, t_xx_zz_y_yz, t_xx_zz_y_zz, t_xx_zz_z_xx, t_xx_zz_z_xy,\
                           t_xx_zz_z_xz, t_xx_zz_z_yy, t_xx_zz_z_yz, t_xx_zz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_zz_z_zz[i] = t_x_xzz_z_zz[i] - rab_x[i] * t_x_zz_z_zz[i];

        t_xx_zz_z_yz[i] = t_x_xzz_z_yz[i] - rab_x[i] * t_x_zz_z_yz[i];

        t_xx_zz_z_yy[i] = t_x_xzz_z_yy[i] - rab_x[i] * t_x_zz_z_yy[i];

        t_xx_zz_z_xz[i] = t_x_xzz_z_xz[i] - rab_x[i] * t_x_zz_z_xz[i];

        t_xx_zz_z_xy[i] = t_x_xzz_z_xy[i] - rab_x[i] * t_x_zz_z_xy[i];

        t_xx_zz_z_xx[i] = t_x_xzz_z_xx[i] - rab_x[i] * t_x_zz_z_xx[i];

        t_xx_zz_y_zz[i] = t_x_xzz_y_zz[i] - rab_x[i] * t_x_zz_y_zz[i];

        t_xx_zz_y_yz[i] = t_x_xzz_y_yz[i] - rab_x[i] * t_x_zz_y_yz[i];

        t_xx_zz_y_yy[i] = t_x_xzz_y_yy[i] - rab_x[i] * t_x_zz_y_yy[i];

        t_xx_zz_y_xz[i] = t_x_xzz_y_xz[i] - rab_x[i] * t_x_zz_y_xz[i];

        t_xx_zz_y_xy[i] = t_x_xzz_y_xy[i] - rab_x[i] * t_x_zz_y_xy[i];

        t_xx_zz_y_xx[i] = t_x_xzz_y_xx[i] - rab_x[i] * t_x_zz_y_xx[i];

        t_xx_zz_x_zz[i] = t_x_xzz_x_zz[i] - rab_x[i] * t_x_zz_x_zz[i];

        t_xx_zz_x_yz[i] = t_x_xzz_x_yz[i] - rab_x[i] * t_x_zz_x_yz[i];

        t_xx_zz_x_yy[i] = t_x_xzz_x_yy[i] - rab_x[i] * t_x_zz_x_yy[i];

        t_xx_zz_x_xz[i] = t_x_xzz_x_xz[i] - rab_x[i] * t_x_zz_x_xz[i];

        t_xx_zz_x_xy[i] = t_x_xzz_x_xy[i] - rab_x[i] * t_x_zz_x_xy[i];

        t_xx_zz_x_xx[i] = t_x_xzz_x_xx[i] - rab_x[i] * t_x_zz_x_xx[i];

        t_xx_yz_z_zz[i] = t_x_xyz_z_zz[i] - rab_x[i] * t_x_yz_z_zz[i];

        t_xx_yz_z_yz[i] = t_x_xyz_z_yz[i] - rab_x[i] * t_x_yz_z_yz[i];

        t_xx_yz_z_yy[i] = t_x_xyz_z_yy[i] - rab_x[i] * t_x_yz_z_yy[i];

        t_xx_yz_z_xz[i] = t_x_xyz_z_xz[i] - rab_x[i] * t_x_yz_z_xz[i];

        t_xx_yz_z_xy[i] = t_x_xyz_z_xy[i] - rab_x[i] * t_x_yz_z_xy[i];

        t_xx_yz_z_xx[i] = t_x_xyz_z_xx[i] - rab_x[i] * t_x_yz_z_xx[i];

        t_xx_yz_y_zz[i] = t_x_xyz_y_zz[i] - rab_x[i] * t_x_yz_y_zz[i];

        t_xx_yz_y_yz[i] = t_x_xyz_y_yz[i] - rab_x[i] * t_x_yz_y_yz[i];

        t_xx_yz_y_yy[i] = t_x_xyz_y_yy[i] - rab_x[i] * t_x_yz_y_yy[i];

        t_xx_yz_y_xz[i] = t_x_xyz_y_xz[i] - rab_x[i] * t_x_yz_y_xz[i];

        t_xx_yz_y_xy[i] = t_x_xyz_y_xy[i] - rab_x[i] * t_x_yz_y_xy[i];

        t_xx_yz_y_xx[i] = t_x_xyz_y_xx[i] - rab_x[i] * t_x_yz_y_xx[i];

        t_xx_yz_x_zz[i] = t_x_xyz_x_zz[i] - rab_x[i] * t_x_yz_x_zz[i];

        t_xx_yz_x_yz[i] = t_x_xyz_x_yz[i] - rab_x[i] * t_x_yz_x_yz[i];

        t_xx_yz_x_yy[i] = t_x_xyz_x_yy[i] - rab_x[i] * t_x_yz_x_yy[i];

        t_xx_yz_x_xz[i] = t_x_xyz_x_xz[i] - rab_x[i] * t_x_yz_x_xz[i];

        t_xx_yz_x_xy[i] = t_x_xyz_x_xy[i] - rab_x[i] * t_x_yz_x_xy[i];

        t_xx_yz_x_xx[i] = t_x_xyz_x_xx[i] - rab_x[i] * t_x_yz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xxz_x_xx, t_x_xxz_x_xy, t_x_xxz_x_xz, t_x_xxz_x_yy,\
                           t_x_xxz_x_yz, t_x_xxz_x_zz, t_x_xxz_y_xx, t_x_xxz_y_xy, t_x_xxz_y_xz,\
                           t_x_xxz_y_yy, t_x_xxz_y_yz, t_x_xxz_y_zz, t_x_xxz_z_xx, t_x_xxz_z_xy,\
                           t_x_xxz_z_xz, t_x_xxz_z_yy, t_x_xxz_z_yz, t_x_xxz_z_zz, t_x_xyy_x_xx,\
                           t_x_xyy_x_xy, t_x_xyy_x_xz, t_x_xyy_x_yy, t_x_xyy_x_yz, t_x_xyy_x_zz,\
                           t_x_xyy_y_xx, t_x_xyy_y_xy, t_x_xyy_y_xz, t_x_xyy_y_yy, t_x_xyy_y_yz,\
                           t_x_xyy_y_zz, t_x_xyy_z_xx, t_x_xyy_z_xy, t_x_xyy_z_xz, t_x_xyy_z_yy,\
                           t_x_xyy_z_yz, t_x_xyy_z_zz, t_x_xz_x_xx, t_x_xz_x_xy, t_x_xz_x_xz,\
                           t_x_xz_x_yy, t_x_xz_x_yz, t_x_xz_x_zz, t_x_xz_y_xx, t_x_xz_y_xy,\
                           t_x_xz_y_xz, t_x_xz_y_yy, t_x_xz_y_yz, t_x_xz_y_zz, t_x_xz_z_xx,\
                           t_x_xz_z_xy, t_x_xz_z_xz, t_x_xz_z_yy, t_x_xz_z_yz, t_x_xz_z_zz,\
                           t_x_yy_x_xx, t_x_yy_x_xy, t_x_yy_x_xz, t_x_yy_x_yy, t_x_yy_x_yz,\
                           t_x_yy_x_zz, t_x_yy_y_xx, t_x_yy_y_xy, t_x_yy_y_xz, t_x_yy_y_yy,\
                           t_x_yy_y_yz, t_x_yy_y_zz, t_x_yy_z_xx, t_x_yy_z_xy, t_x_yy_z_xz,\
                           t_x_yy_z_yy, t_x_yy_z_yz, t_x_yy_z_zz, t_xx_xz_x_xx, t_xx_xz_x_xy,\
                           t_xx_xz_x_xz, t_xx_xz_x_yy, t_xx_xz_x_yz, t_xx_xz_x_zz, t_xx_xz_y_xx,\
                           t_xx_xz_y_xy, t_xx_xz_y_xz, t_xx_xz_y_yy, t_xx_xz_y_yz, t_xx_xz_y_zz,\
                           t_xx_xz_z_xx, t_xx_xz_z_xy, t_xx_xz_z_xz, t_xx_xz_z_yy, t_xx_xz_z_yz,\
                           t_xx_xz_z_zz, t_xx_yy_x_xx, t_xx_yy_x_xy, t_xx_yy_x_xz, t_xx_yy_x_yy,\
                           t_xx_yy_x_yz, t_xx_yy_x_zz, t_xx_yy_y_xx, t_xx_yy_y_xy, t_xx_yy_y_xz,\
                           t_xx_yy_y_yy, t_xx_yy_y_yz, t_xx_yy_y_zz, t_xx_yy_z_xx, t_xx_yy_z_xy,\
                           t_xx_yy_z_xz, t_xx_yy_z_yy, t_xx_yy_z_yz, t_xx_yy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_yy_z_zz[i] = t_x_xyy_z_zz[i] - rab_x[i] * t_x_yy_z_zz[i];

        t_xx_yy_z_yz[i] = t_x_xyy_z_yz[i] - rab_x[i] * t_x_yy_z_yz[i];

        t_xx_yy_z_yy[i] = t_x_xyy_z_yy[i] - rab_x[i] * t_x_yy_z_yy[i];

        t_xx_yy_z_xz[i] = t_x_xyy_z_xz[i] - rab_x[i] * t_x_yy_z_xz[i];

        t_xx_yy_z_xy[i] = t_x_xyy_z_xy[i] - rab_x[i] * t_x_yy_z_xy[i];

        t_xx_yy_z_xx[i] = t_x_xyy_z_xx[i] - rab_x[i] * t_x_yy_z_xx[i];

        t_xx_yy_y_zz[i] = t_x_xyy_y_zz[i] - rab_x[i] * t_x_yy_y_zz[i];

        t_xx_yy_y_yz[i] = t_x_xyy_y_yz[i] - rab_x[i] * t_x_yy_y_yz[i];

        t_xx_yy_y_yy[i] = t_x_xyy_y_yy[i] - rab_x[i] * t_x_yy_y_yy[i];

        t_xx_yy_y_xz[i] = t_x_xyy_y_xz[i] - rab_x[i] * t_x_yy_y_xz[i];

        t_xx_yy_y_xy[i] = t_x_xyy_y_xy[i] - rab_x[i] * t_x_yy_y_xy[i];

        t_xx_yy_y_xx[i] = t_x_xyy_y_xx[i] - rab_x[i] * t_x_yy_y_xx[i];

        t_xx_yy_x_zz[i] = t_x_xyy_x_zz[i] - rab_x[i] * t_x_yy_x_zz[i];

        t_xx_yy_x_yz[i] = t_x_xyy_x_yz[i] - rab_x[i] * t_x_yy_x_yz[i];

        t_xx_yy_x_yy[i] = t_x_xyy_x_yy[i] - rab_x[i] * t_x_yy_x_yy[i];

        t_xx_yy_x_xz[i] = t_x_xyy_x_xz[i] - rab_x[i] * t_x_yy_x_xz[i];

        t_xx_yy_x_xy[i] = t_x_xyy_x_xy[i] - rab_x[i] * t_x_yy_x_xy[i];

        t_xx_yy_x_xx[i] = t_x_xyy_x_xx[i] - rab_x[i] * t_x_yy_x_xx[i];

        t_xx_xz_z_zz[i] = t_x_xxz_z_zz[i] - rab_x[i] * t_x_xz_z_zz[i];

        t_xx_xz_z_yz[i] = t_x_xxz_z_yz[i] - rab_x[i] * t_x_xz_z_yz[i];

        t_xx_xz_z_yy[i] = t_x_xxz_z_yy[i] - rab_x[i] * t_x_xz_z_yy[i];

        t_xx_xz_z_xz[i] = t_x_xxz_z_xz[i] - rab_x[i] * t_x_xz_z_xz[i];

        t_xx_xz_z_xy[i] = t_x_xxz_z_xy[i] - rab_x[i] * t_x_xz_z_xy[i];

        t_xx_xz_z_xx[i] = t_x_xxz_z_xx[i] - rab_x[i] * t_x_xz_z_xx[i];

        t_xx_xz_y_zz[i] = t_x_xxz_y_zz[i] - rab_x[i] * t_x_xz_y_zz[i];

        t_xx_xz_y_yz[i] = t_x_xxz_y_yz[i] - rab_x[i] * t_x_xz_y_yz[i];

        t_xx_xz_y_yy[i] = t_x_xxz_y_yy[i] - rab_x[i] * t_x_xz_y_yy[i];

        t_xx_xz_y_xz[i] = t_x_xxz_y_xz[i] - rab_x[i] * t_x_xz_y_xz[i];

        t_xx_xz_y_xy[i] = t_x_xxz_y_xy[i] - rab_x[i] * t_x_xz_y_xy[i];

        t_xx_xz_y_xx[i] = t_x_xxz_y_xx[i] - rab_x[i] * t_x_xz_y_xx[i];

        t_xx_xz_x_zz[i] = t_x_xxz_x_zz[i] - rab_x[i] * t_x_xz_x_zz[i];

        t_xx_xz_x_yz[i] = t_x_xxz_x_yz[i] - rab_x[i] * t_x_xz_x_yz[i];

        t_xx_xz_x_yy[i] = t_x_xxz_x_yy[i] - rab_x[i] * t_x_xz_x_yy[i];

        t_xx_xz_x_xz[i] = t_x_xxz_x_xz[i] - rab_x[i] * t_x_xz_x_xz[i];

        t_xx_xz_x_xy[i] = t_x_xxz_x_xy[i] - rab_x[i] * t_x_xz_x_xy[i];

        t_xx_xz_x_xx[i] = t_x_xxz_x_xx[i] - rab_x[i] * t_x_xz_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_x_xx_x_xx, t_x_xx_x_xy, t_x_xx_x_xz, t_x_xx_x_yy,\
                           t_x_xx_x_yz, t_x_xx_x_zz, t_x_xx_y_xx, t_x_xx_y_xy, t_x_xx_y_xz,\
                           t_x_xx_y_yy, t_x_xx_y_yz, t_x_xx_y_zz, t_x_xx_z_xx, t_x_xx_z_xy,\
                           t_x_xx_z_xz, t_x_xx_z_yy, t_x_xx_z_yz, t_x_xx_z_zz, t_x_xxx_x_xx,\
                           t_x_xxx_x_xy, t_x_xxx_x_xz, t_x_xxx_x_yy, t_x_xxx_x_yz, t_x_xxx_x_zz,\
                           t_x_xxx_y_xx, t_x_xxx_y_xy, t_x_xxx_y_xz, t_x_xxx_y_yy, t_x_xxx_y_yz,\
                           t_x_xxx_y_zz, t_x_xxx_z_xx, t_x_xxx_z_xy, t_x_xxx_z_xz, t_x_xxx_z_yy,\
                           t_x_xxx_z_yz, t_x_xxx_z_zz, t_x_xxy_x_xx, t_x_xxy_x_xy, t_x_xxy_x_xz,\
                           t_x_xxy_x_yy, t_x_xxy_x_yz, t_x_xxy_x_zz, t_x_xxy_y_xx, t_x_xxy_y_xy,\
                           t_x_xxy_y_xz, t_x_xxy_y_yy, t_x_xxy_y_yz, t_x_xxy_y_zz, t_x_xxy_z_xx,\
                           t_x_xxy_z_xy, t_x_xxy_z_xz, t_x_xxy_z_yy, t_x_xxy_z_yz, t_x_xxy_z_zz,\
                           t_x_xy_x_xx, t_x_xy_x_xy, t_x_xy_x_xz, t_x_xy_x_yy, t_x_xy_x_yz,\
                           t_x_xy_x_zz, t_x_xy_y_xx, t_x_xy_y_xy, t_x_xy_y_xz, t_x_xy_y_yy,\
                           t_x_xy_y_yz, t_x_xy_y_zz, t_x_xy_z_xx, t_x_xy_z_xy, t_x_xy_z_xz,\
                           t_x_xy_z_yy, t_x_xy_z_yz, t_x_xy_z_zz, t_xx_xx_x_xx, t_xx_xx_x_xy,\
                           t_xx_xx_x_xz, t_xx_xx_x_yy, t_xx_xx_x_yz, t_xx_xx_x_zz, t_xx_xx_y_xx,\
                           t_xx_xx_y_xy, t_xx_xx_y_xz, t_xx_xx_y_yy, t_xx_xx_y_yz, t_xx_xx_y_zz,\
                           t_xx_xx_z_xx, t_xx_xx_z_xy, t_xx_xx_z_xz, t_xx_xx_z_yy, t_xx_xx_z_yz,\
                           t_xx_xx_z_zz, t_xx_xy_x_xx, t_xx_xy_x_xy, t_xx_xy_x_xz, t_xx_xy_x_yy,\
                           t_xx_xy_x_yz, t_xx_xy_x_zz, t_xx_xy_y_xx, t_xx_xy_y_xy, t_xx_xy_y_xz,\
                           t_xx_xy_y_yy, t_xx_xy_y_yz, t_xx_xy_y_zz, t_xx_xy_z_xx, t_xx_xy_z_xy,\
                           t_xx_xy_z_xz, t_xx_xy_z_yy, t_xx_xy_z_yz, t_xx_xy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_xy_z_zz[i] = t_x_xxy_z_zz[i] - rab_x[i] * t_x_xy_z_zz[i];

        t_xx_xy_z_yz[i] = t_x_xxy_z_yz[i] - rab_x[i] * t_x_xy_z_yz[i];

        t_xx_xy_z_yy[i] = t_x_xxy_z_yy[i] - rab_x[i] * t_x_xy_z_yy[i];

        t_xx_xy_z_xz[i] = t_x_xxy_z_xz[i] - rab_x[i] * t_x_xy_z_xz[i];

        t_xx_xy_z_xy[i] = t_x_xxy_z_xy[i] - rab_x[i] * t_x_xy_z_xy[i];

        t_xx_xy_z_xx[i] = t_x_xxy_z_xx[i] - rab_x[i] * t_x_xy_z_xx[i];

        t_xx_xy_y_zz[i] = t_x_xxy_y_zz[i] - rab_x[i] * t_x_xy_y_zz[i];

        t_xx_xy_y_yz[i] = t_x_xxy_y_yz[i] - rab_x[i] * t_x_xy_y_yz[i];

        t_xx_xy_y_yy[i] = t_x_xxy_y_yy[i] - rab_x[i] * t_x_xy_y_yy[i];

        t_xx_xy_y_xz[i] = t_x_xxy_y_xz[i] - rab_x[i] * t_x_xy_y_xz[i];

        t_xx_xy_y_xy[i] = t_x_xxy_y_xy[i] - rab_x[i] * t_x_xy_y_xy[i];

        t_xx_xy_y_xx[i] = t_x_xxy_y_xx[i] - rab_x[i] * t_x_xy_y_xx[i];

        t_xx_xy_x_zz[i] = t_x_xxy_x_zz[i] - rab_x[i] * t_x_xy_x_zz[i];

        t_xx_xy_x_yz[i] = t_x_xxy_x_yz[i] - rab_x[i] * t_x_xy_x_yz[i];

        t_xx_xy_x_yy[i] = t_x_xxy_x_yy[i] - rab_x[i] * t_x_xy_x_yy[i];

        t_xx_xy_x_xz[i] = t_x_xxy_x_xz[i] - rab_x[i] * t_x_xy_x_xz[i];

        t_xx_xy_x_xy[i] = t_x_xxy_x_xy[i] - rab_x[i] * t_x_xy_x_xy[i];

        t_xx_xy_x_xx[i] = t_x_xxy_x_xx[i] - rab_x[i] * t_x_xy_x_xx[i];

        t_xx_xx_z_zz[i] = t_x_xxx_z_zz[i] - rab_x[i] * t_x_xx_z_zz[i];

        t_xx_xx_z_yz[i] = t_x_xxx_z_yz[i] - rab_x[i] * t_x_xx_z_yz[i];

        t_xx_xx_z_yy[i] = t_x_xxx_z_yy[i] - rab_x[i] * t_x_xx_z_yy[i];

        t_xx_xx_z_xz[i] = t_x_xxx_z_xz[i] - rab_x[i] * t_x_xx_z_xz[i];

        t_xx_xx_z_xy[i] = t_x_xxx_z_xy[i] - rab_x[i] * t_x_xx_z_xy[i];

        t_xx_xx_z_xx[i] = t_x_xxx_z_xx[i] - rab_x[i] * t_x_xx_z_xx[i];

        t_xx_xx_y_zz[i] = t_x_xxx_y_zz[i] - rab_x[i] * t_x_xx_y_zz[i];

        t_xx_xx_y_yz[i] = t_x_xxx_y_yz[i] - rab_x[i] * t_x_xx_y_yz[i];

        t_xx_xx_y_yy[i] = t_x_xxx_y_yy[i] - rab_x[i] * t_x_xx_y_yy[i];

        t_xx_xx_y_xz[i] = t_x_xxx_y_xz[i] - rab_x[i] * t_x_xx_y_xz[i];

        t_xx_xx_y_xy[i] = t_x_xxx_y_xy[i] - rab_x[i] * t_x_xx_y_xy[i];

        t_xx_xx_y_xx[i] = t_x_xxx_y_xx[i] - rab_x[i] * t_x_xx_y_xx[i];

        t_xx_xx_x_zz[i] = t_x_xxx_x_zz[i] - rab_x[i] * t_x_xx_x_zz[i];

        t_xx_xx_x_yz[i] = t_x_xxx_x_yz[i] - rab_x[i] * t_x_xx_x_yz[i];

        t_xx_xx_x_yy[i] = t_x_xxx_x_yy[i] - rab_x[i] * t_x_xx_x_yy[i];

        t_xx_xx_x_xz[i] = t_x_xxx_x_xz[i] - rab_x[i] * t_x_xx_x_xz[i];

        t_xx_xx_x_xy[i] = t_x_xxx_x_xy[i] - rab_x[i] * t_x_xx_x_xy[i];

        t_xx_xx_x_xx[i] = t_x_xxx_x_xx[i] - rab_x[i] * t_x_xx_x_xx[i];
    }
}


} // derirec namespace
