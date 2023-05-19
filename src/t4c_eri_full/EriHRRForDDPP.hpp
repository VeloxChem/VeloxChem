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
compHostHRRForDDPP_V0(      BufferHostXY<T>&      intsBufferDDPP,
                      const BufferHostX<int32_t>& intsIndexesDDPP,
                      const BufferHostXY<T>&      intsBufferPDPP,
                      const BufferHostX<int32_t>& intsIndexesPDPP,
                      const BufferHostXY<T>&      intsBufferPFPP,
                      const BufferHostX<int32_t>& intsIndexesPFPP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDPP) integral components

    t_zz_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(0));

    t_zz_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(1));

    t_zz_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(2));

    t_zz_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(3));

    t_zz_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(4));

    t_zz_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(5));

    t_zz_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(6));

    t_zz_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(7));

    t_zz_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(8));

    t_zz_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(9));

    t_zz_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(10));

    t_zz_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(11));

    t_zz_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(12));

    t_zz_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(13));

    t_zz_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(14));

    t_zz_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(15));

    t_zz_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(16));

    t_zz_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(17));

    t_zz_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(18));

    t_zz_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(19));

    t_zz_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(20));

    t_zz_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(21));

    t_zz_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(22));

    t_zz_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(23));

    t_zz_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(24));

    t_zz_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(25));

    t_zz_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(26));

    t_zz_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(27));

    t_zz_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(28));

    t_zz_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(29));

    t_zz_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(30));

    t_zz_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(31));

    t_zz_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(32));

    t_zz_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(33));

    t_zz_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(34));

    t_zz_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(35));

    t_zz_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(36));

    t_zz_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(37));

    t_zz_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(38));

    t_zz_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(39));

    t_zz_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(40));

    t_zz_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(41));

    t_zz_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(42));

    t_zz_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(43));

    t_zz_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(44));

    t_zz_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(45));

    t_zz_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(46));

    t_zz_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(47));

    t_zz_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(48));

    t_zz_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(49));

    t_zz_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(50));

    t_zz_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(51));

    t_zz_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(52));

    t_zz_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(53));

    t_yz_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(54));

    t_yz_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(55));

    t_yz_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(56));

    t_yz_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(57));

    t_yz_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(58));

    t_yz_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(59));

    t_yz_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(60));

    t_yz_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(61));

    t_yz_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(62));

    t_yz_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(63));

    t_yz_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(64));

    t_yz_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(65));

    t_yz_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(66));

    t_yz_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(67));

    t_yz_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(68));

    t_yz_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(69));

    t_yz_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(70));

    t_yz_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(71));

    t_yz_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(72));

    t_yz_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(73));

    t_yz_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(74));

    t_yz_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(75));

    t_yz_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(76));

    t_yz_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(77));

    t_yz_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(78));

    t_yz_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(79));

    t_yz_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(80));

    t_yz_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(81));

    t_yz_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(82));

    t_yz_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(83));

    t_yz_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(84));

    t_yz_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(85));

    t_yz_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(86));

    t_yz_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(87));

    t_yz_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(88));

    t_yz_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(89));

    t_yz_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(90));

    t_yz_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(91));

    t_yz_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(92));

    t_yz_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(93));

    t_yz_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(94));

    t_yz_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(95));

    t_yz_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(96));

    t_yz_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(97));

    t_yz_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(98));

    t_yz_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(99));

    t_yz_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(100));

    t_yz_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(101));

    t_yz_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(102));

    t_yz_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(103));

    t_yz_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(104));

    t_yz_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(105));

    t_yz_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(106));

    t_yz_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(107));

    t_yy_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(108));

    t_yy_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(109));

    t_yy_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(110));

    t_yy_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(111));

    t_yy_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(112));

    t_yy_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(113));

    t_yy_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(114));

    t_yy_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(115));

    t_yy_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(116));

    t_yy_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(117));

    t_yy_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(118));

    t_yy_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(119));

    t_yy_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(120));

    t_yy_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(121));

    t_yy_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(122));

    t_yy_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(123));

    t_yy_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(124));

    t_yy_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(125));

    t_yy_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(126));

    t_yy_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(127));

    t_yy_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(128));

    t_yy_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(129));

    t_yy_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(130));

    t_yy_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(131));

    t_yy_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(132));

    t_yy_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(133));

    t_yy_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(134));

    t_yy_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(135));

    t_yy_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(136));

    t_yy_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(137));

    t_yy_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(138));

    t_yy_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(139));

    t_yy_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(140));

    t_yy_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(141));

    t_yy_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(142));

    t_yy_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(143));

    t_yy_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(144));

    t_yy_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(145));

    t_yy_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(146));

    t_yy_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(147));

    t_yy_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(148));

    t_yy_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(149));

    t_yy_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(150));

    t_yy_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(151));

    t_yy_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(152));

    t_yy_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(153));

    t_yy_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(154));

    t_yy_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(155));

    t_yy_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(156));

    t_yy_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(157));

    t_yy_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(158));

    t_yy_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(159));

    t_yy_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(160));

    t_yy_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(161));

    t_xz_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(162));

    t_xz_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(163));

    t_xz_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(164));

    t_xz_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(165));

    t_xz_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(166));

    t_xz_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(167));

    t_xz_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(168));

    t_xz_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(169));

    t_xz_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(170));

    t_xz_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(171));

    t_xz_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(172));

    t_xz_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(173));

    t_xz_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(174));

    t_xz_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(175));

    t_xz_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(176));

    t_xz_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(177));

    t_xz_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(178));

    t_xz_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(179));

    t_xz_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(180));

    t_xz_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(181));

    t_xz_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(182));

    t_xz_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(183));

    t_xz_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(184));

    t_xz_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(185));

    t_xz_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(186));

    t_xz_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(187));

    t_xz_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(188));

    t_xz_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(189));

    t_xz_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(190));

    t_xz_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(191));

    t_xz_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(192));

    t_xz_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(193));

    t_xz_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(194));

    t_xz_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(195));

    t_xz_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(196));

    t_xz_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(197));

    t_xz_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(198));

    t_xz_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(199));

    t_xz_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(200));

    t_xz_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(201));

    t_xz_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(202));

    t_xz_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(203));

    t_xz_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(204));

    t_xz_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(205));

    t_xz_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(206));

    t_xz_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(207));

    t_xz_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(208));

    t_xz_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(209));

    t_xz_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(210));

    t_xz_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(211));

    t_xz_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(212));

    t_xz_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(213));

    t_xz_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(214));

    t_xz_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(215));

    t_xy_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(216));

    t_xy_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(217));

    t_xy_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(218));

    t_xy_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(219));

    t_xy_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(220));

    t_xy_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(221));

    t_xy_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(222));

    t_xy_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(223));

    t_xy_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(224));

    t_xy_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(225));

    t_xy_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(226));

    t_xy_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(227));

    t_xy_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(228));

    t_xy_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(229));

    t_xy_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(230));

    t_xy_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(231));

    t_xy_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(232));

    t_xy_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(233));

    t_xy_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(234));

    t_xy_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(235));

    t_xy_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(236));

    t_xy_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(237));

    t_xy_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(238));

    t_xy_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(239));

    t_xy_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(240));

    t_xy_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(241));

    t_xy_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(242));

    t_xy_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(243));

    t_xy_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(244));

    t_xy_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(245));

    t_xy_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(246));

    t_xy_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(247));

    t_xy_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(248));

    t_xy_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(249));

    t_xy_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(250));

    t_xy_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(251));

    t_xy_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(252));

    t_xy_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(253));

    t_xy_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(254));

    t_xy_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(255));

    t_xy_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(256));

    t_xy_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(257));

    t_xy_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(258));

    t_xy_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(259));

    t_xy_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(260));

    t_xy_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(261));

    t_xy_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(262));

    t_xy_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(263));

    t_xy_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(264));

    t_xy_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(265));

    t_xy_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(266));

    t_xy_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(267));

    t_xy_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(268));

    t_xy_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(269));

    t_xx_zz_z_z = intsBufferDDPP.data(intsIndexesDDPP(270));

    t_xx_zz_z_y = intsBufferDDPP.data(intsIndexesDDPP(271));

    t_xx_zz_z_x = intsBufferDDPP.data(intsIndexesDDPP(272));

    t_xx_zz_y_z = intsBufferDDPP.data(intsIndexesDDPP(273));

    t_xx_zz_y_y = intsBufferDDPP.data(intsIndexesDDPP(274));

    t_xx_zz_y_x = intsBufferDDPP.data(intsIndexesDDPP(275));

    t_xx_zz_x_z = intsBufferDDPP.data(intsIndexesDDPP(276));

    t_xx_zz_x_y = intsBufferDDPP.data(intsIndexesDDPP(277));

    t_xx_zz_x_x = intsBufferDDPP.data(intsIndexesDDPP(278));

    t_xx_yz_z_z = intsBufferDDPP.data(intsIndexesDDPP(279));

    t_xx_yz_z_y = intsBufferDDPP.data(intsIndexesDDPP(280));

    t_xx_yz_z_x = intsBufferDDPP.data(intsIndexesDDPP(281));

    t_xx_yz_y_z = intsBufferDDPP.data(intsIndexesDDPP(282));

    t_xx_yz_y_y = intsBufferDDPP.data(intsIndexesDDPP(283));

    t_xx_yz_y_x = intsBufferDDPP.data(intsIndexesDDPP(284));

    t_xx_yz_x_z = intsBufferDDPP.data(intsIndexesDDPP(285));

    t_xx_yz_x_y = intsBufferDDPP.data(intsIndexesDDPP(286));

    t_xx_yz_x_x = intsBufferDDPP.data(intsIndexesDDPP(287));

    t_xx_yy_z_z = intsBufferDDPP.data(intsIndexesDDPP(288));

    t_xx_yy_z_y = intsBufferDDPP.data(intsIndexesDDPP(289));

    t_xx_yy_z_x = intsBufferDDPP.data(intsIndexesDDPP(290));

    t_xx_yy_y_z = intsBufferDDPP.data(intsIndexesDDPP(291));

    t_xx_yy_y_y = intsBufferDDPP.data(intsIndexesDDPP(292));

    t_xx_yy_y_x = intsBufferDDPP.data(intsIndexesDDPP(293));

    t_xx_yy_x_z = intsBufferDDPP.data(intsIndexesDDPP(294));

    t_xx_yy_x_y = intsBufferDDPP.data(intsIndexesDDPP(295));

    t_xx_yy_x_x = intsBufferDDPP.data(intsIndexesDDPP(296));

    t_xx_xz_z_z = intsBufferDDPP.data(intsIndexesDDPP(297));

    t_xx_xz_z_y = intsBufferDDPP.data(intsIndexesDDPP(298));

    t_xx_xz_z_x = intsBufferDDPP.data(intsIndexesDDPP(299));

    t_xx_xz_y_z = intsBufferDDPP.data(intsIndexesDDPP(300));

    t_xx_xz_y_y = intsBufferDDPP.data(intsIndexesDDPP(301));

    t_xx_xz_y_x = intsBufferDDPP.data(intsIndexesDDPP(302));

    t_xx_xz_x_z = intsBufferDDPP.data(intsIndexesDDPP(303));

    t_xx_xz_x_y = intsBufferDDPP.data(intsIndexesDDPP(304));

    t_xx_xz_x_x = intsBufferDDPP.data(intsIndexesDDPP(305));

    t_xx_xy_z_z = intsBufferDDPP.data(intsIndexesDDPP(306));

    t_xx_xy_z_y = intsBufferDDPP.data(intsIndexesDDPP(307));

    t_xx_xy_z_x = intsBufferDDPP.data(intsIndexesDDPP(308));

    t_xx_xy_y_z = intsBufferDDPP.data(intsIndexesDDPP(309));

    t_xx_xy_y_y = intsBufferDDPP.data(intsIndexesDDPP(310));

    t_xx_xy_y_x = intsBufferDDPP.data(intsIndexesDDPP(311));

    t_xx_xy_x_z = intsBufferDDPP.data(intsIndexesDDPP(312));

    t_xx_xy_x_y = intsBufferDDPP.data(intsIndexesDDPP(313));

    t_xx_xy_x_x = intsBufferDDPP.data(intsIndexesDDPP(314));

    t_xx_xx_z_z = intsBufferDDPP.data(intsIndexesDDPP(315));

    t_xx_xx_z_y = intsBufferDDPP.data(intsIndexesDDPP(316));

    t_xx_xx_z_x = intsBufferDDPP.data(intsIndexesDDPP(317));

    t_xx_xx_y_z = intsBufferDDPP.data(intsIndexesDDPP(318));

    t_xx_xx_y_y = intsBufferDDPP.data(intsIndexesDDPP(319));

    t_xx_xx_y_x = intsBufferDDPP.data(intsIndexesDDPP(320));

    t_xx_xx_x_z = intsBufferDDPP.data(intsIndexesDDPP(321));

    t_xx_xx_x_y = intsBufferDDPP.data(intsIndexesDDPP(322));

    t_xx_xx_x_x = intsBufferDDPP.data(intsIndexesDDPP(323));

    // set up (PDPP) integral components

    t_z_zz_z_z = intsBufferPDPP.data(intsIndexesPDPP(0));

    t_z_zz_z_y = intsBufferPDPP.data(intsIndexesPDPP(1));

    t_z_zz_z_x = intsBufferPDPP.data(intsIndexesPDPP(2));

    t_z_zz_y_z = intsBufferPDPP.data(intsIndexesPDPP(3));

    t_z_zz_y_y = intsBufferPDPP.data(intsIndexesPDPP(4));

    t_z_zz_y_x = intsBufferPDPP.data(intsIndexesPDPP(5));

    t_z_zz_x_z = intsBufferPDPP.data(intsIndexesPDPP(6));

    t_z_zz_x_y = intsBufferPDPP.data(intsIndexesPDPP(7));

    t_z_zz_x_x = intsBufferPDPP.data(intsIndexesPDPP(8));

    t_z_yz_z_z = intsBufferPDPP.data(intsIndexesPDPP(9));

    t_z_yz_z_y = intsBufferPDPP.data(intsIndexesPDPP(10));

    t_z_yz_z_x = intsBufferPDPP.data(intsIndexesPDPP(11));

    t_z_yz_y_z = intsBufferPDPP.data(intsIndexesPDPP(12));

    t_z_yz_y_y = intsBufferPDPP.data(intsIndexesPDPP(13));

    t_z_yz_y_x = intsBufferPDPP.data(intsIndexesPDPP(14));

    t_z_yz_x_z = intsBufferPDPP.data(intsIndexesPDPP(15));

    t_z_yz_x_y = intsBufferPDPP.data(intsIndexesPDPP(16));

    t_z_yz_x_x = intsBufferPDPP.data(intsIndexesPDPP(17));

    t_z_yy_z_z = intsBufferPDPP.data(intsIndexesPDPP(18));

    t_z_yy_z_y = intsBufferPDPP.data(intsIndexesPDPP(19));

    t_z_yy_z_x = intsBufferPDPP.data(intsIndexesPDPP(20));

    t_z_yy_y_z = intsBufferPDPP.data(intsIndexesPDPP(21));

    t_z_yy_y_y = intsBufferPDPP.data(intsIndexesPDPP(22));

    t_z_yy_y_x = intsBufferPDPP.data(intsIndexesPDPP(23));

    t_z_yy_x_z = intsBufferPDPP.data(intsIndexesPDPP(24));

    t_z_yy_x_y = intsBufferPDPP.data(intsIndexesPDPP(25));

    t_z_yy_x_x = intsBufferPDPP.data(intsIndexesPDPP(26));

    t_z_xz_z_z = intsBufferPDPP.data(intsIndexesPDPP(27));

    t_z_xz_z_y = intsBufferPDPP.data(intsIndexesPDPP(28));

    t_z_xz_z_x = intsBufferPDPP.data(intsIndexesPDPP(29));

    t_z_xz_y_z = intsBufferPDPP.data(intsIndexesPDPP(30));

    t_z_xz_y_y = intsBufferPDPP.data(intsIndexesPDPP(31));

    t_z_xz_y_x = intsBufferPDPP.data(intsIndexesPDPP(32));

    t_z_xz_x_z = intsBufferPDPP.data(intsIndexesPDPP(33));

    t_z_xz_x_y = intsBufferPDPP.data(intsIndexesPDPP(34));

    t_z_xz_x_x = intsBufferPDPP.data(intsIndexesPDPP(35));

    t_z_xy_z_z = intsBufferPDPP.data(intsIndexesPDPP(36));

    t_z_xy_z_y = intsBufferPDPP.data(intsIndexesPDPP(37));

    t_z_xy_z_x = intsBufferPDPP.data(intsIndexesPDPP(38));

    t_z_xy_y_z = intsBufferPDPP.data(intsIndexesPDPP(39));

    t_z_xy_y_y = intsBufferPDPP.data(intsIndexesPDPP(40));

    t_z_xy_y_x = intsBufferPDPP.data(intsIndexesPDPP(41));

    t_z_xy_x_z = intsBufferPDPP.data(intsIndexesPDPP(42));

    t_z_xy_x_y = intsBufferPDPP.data(intsIndexesPDPP(43));

    t_z_xy_x_x = intsBufferPDPP.data(intsIndexesPDPP(44));

    t_z_xx_z_z = intsBufferPDPP.data(intsIndexesPDPP(45));

    t_z_xx_z_y = intsBufferPDPP.data(intsIndexesPDPP(46));

    t_z_xx_z_x = intsBufferPDPP.data(intsIndexesPDPP(47));

    t_z_xx_y_z = intsBufferPDPP.data(intsIndexesPDPP(48));

    t_z_xx_y_y = intsBufferPDPP.data(intsIndexesPDPP(49));

    t_z_xx_y_x = intsBufferPDPP.data(intsIndexesPDPP(50));

    t_z_xx_x_z = intsBufferPDPP.data(intsIndexesPDPP(51));

    t_z_xx_x_y = intsBufferPDPP.data(intsIndexesPDPP(52));

    t_z_xx_x_x = intsBufferPDPP.data(intsIndexesPDPP(53));

    t_y_zz_z_z = intsBufferPDPP.data(intsIndexesPDPP(54));

    t_y_zz_z_y = intsBufferPDPP.data(intsIndexesPDPP(55));

    t_y_zz_z_x = intsBufferPDPP.data(intsIndexesPDPP(56));

    t_y_zz_y_z = intsBufferPDPP.data(intsIndexesPDPP(57));

    t_y_zz_y_y = intsBufferPDPP.data(intsIndexesPDPP(58));

    t_y_zz_y_x = intsBufferPDPP.data(intsIndexesPDPP(59));

    t_y_zz_x_z = intsBufferPDPP.data(intsIndexesPDPP(60));

    t_y_zz_x_y = intsBufferPDPP.data(intsIndexesPDPP(61));

    t_y_zz_x_x = intsBufferPDPP.data(intsIndexesPDPP(62));

    t_y_yz_z_z = intsBufferPDPP.data(intsIndexesPDPP(63));

    t_y_yz_z_y = intsBufferPDPP.data(intsIndexesPDPP(64));

    t_y_yz_z_x = intsBufferPDPP.data(intsIndexesPDPP(65));

    t_y_yz_y_z = intsBufferPDPP.data(intsIndexesPDPP(66));

    t_y_yz_y_y = intsBufferPDPP.data(intsIndexesPDPP(67));

    t_y_yz_y_x = intsBufferPDPP.data(intsIndexesPDPP(68));

    t_y_yz_x_z = intsBufferPDPP.data(intsIndexesPDPP(69));

    t_y_yz_x_y = intsBufferPDPP.data(intsIndexesPDPP(70));

    t_y_yz_x_x = intsBufferPDPP.data(intsIndexesPDPP(71));

    t_y_yy_z_z = intsBufferPDPP.data(intsIndexesPDPP(72));

    t_y_yy_z_y = intsBufferPDPP.data(intsIndexesPDPP(73));

    t_y_yy_z_x = intsBufferPDPP.data(intsIndexesPDPP(74));

    t_y_yy_y_z = intsBufferPDPP.data(intsIndexesPDPP(75));

    t_y_yy_y_y = intsBufferPDPP.data(intsIndexesPDPP(76));

    t_y_yy_y_x = intsBufferPDPP.data(intsIndexesPDPP(77));

    t_y_yy_x_z = intsBufferPDPP.data(intsIndexesPDPP(78));

    t_y_yy_x_y = intsBufferPDPP.data(intsIndexesPDPP(79));

    t_y_yy_x_x = intsBufferPDPP.data(intsIndexesPDPP(80));

    t_y_xz_z_z = intsBufferPDPP.data(intsIndexesPDPP(81));

    t_y_xz_z_y = intsBufferPDPP.data(intsIndexesPDPP(82));

    t_y_xz_z_x = intsBufferPDPP.data(intsIndexesPDPP(83));

    t_y_xz_y_z = intsBufferPDPP.data(intsIndexesPDPP(84));

    t_y_xz_y_y = intsBufferPDPP.data(intsIndexesPDPP(85));

    t_y_xz_y_x = intsBufferPDPP.data(intsIndexesPDPP(86));

    t_y_xz_x_z = intsBufferPDPP.data(intsIndexesPDPP(87));

    t_y_xz_x_y = intsBufferPDPP.data(intsIndexesPDPP(88));

    t_y_xz_x_x = intsBufferPDPP.data(intsIndexesPDPP(89));

    t_y_xy_z_z = intsBufferPDPP.data(intsIndexesPDPP(90));

    t_y_xy_z_y = intsBufferPDPP.data(intsIndexesPDPP(91));

    t_y_xy_z_x = intsBufferPDPP.data(intsIndexesPDPP(92));

    t_y_xy_y_z = intsBufferPDPP.data(intsIndexesPDPP(93));

    t_y_xy_y_y = intsBufferPDPP.data(intsIndexesPDPP(94));

    t_y_xy_y_x = intsBufferPDPP.data(intsIndexesPDPP(95));

    t_y_xy_x_z = intsBufferPDPP.data(intsIndexesPDPP(96));

    t_y_xy_x_y = intsBufferPDPP.data(intsIndexesPDPP(97));

    t_y_xy_x_x = intsBufferPDPP.data(intsIndexesPDPP(98));

    t_y_xx_z_z = intsBufferPDPP.data(intsIndexesPDPP(99));

    t_y_xx_z_y = intsBufferPDPP.data(intsIndexesPDPP(100));

    t_y_xx_z_x = intsBufferPDPP.data(intsIndexesPDPP(101));

    t_y_xx_y_z = intsBufferPDPP.data(intsIndexesPDPP(102));

    t_y_xx_y_y = intsBufferPDPP.data(intsIndexesPDPP(103));

    t_y_xx_y_x = intsBufferPDPP.data(intsIndexesPDPP(104));

    t_y_xx_x_z = intsBufferPDPP.data(intsIndexesPDPP(105));

    t_y_xx_x_y = intsBufferPDPP.data(intsIndexesPDPP(106));

    t_y_xx_x_x = intsBufferPDPP.data(intsIndexesPDPP(107));

    t_x_zz_z_z = intsBufferPDPP.data(intsIndexesPDPP(108));

    t_x_zz_z_y = intsBufferPDPP.data(intsIndexesPDPP(109));

    t_x_zz_z_x = intsBufferPDPP.data(intsIndexesPDPP(110));

    t_x_zz_y_z = intsBufferPDPP.data(intsIndexesPDPP(111));

    t_x_zz_y_y = intsBufferPDPP.data(intsIndexesPDPP(112));

    t_x_zz_y_x = intsBufferPDPP.data(intsIndexesPDPP(113));

    t_x_zz_x_z = intsBufferPDPP.data(intsIndexesPDPP(114));

    t_x_zz_x_y = intsBufferPDPP.data(intsIndexesPDPP(115));

    t_x_zz_x_x = intsBufferPDPP.data(intsIndexesPDPP(116));

    t_x_yz_z_z = intsBufferPDPP.data(intsIndexesPDPP(117));

    t_x_yz_z_y = intsBufferPDPP.data(intsIndexesPDPP(118));

    t_x_yz_z_x = intsBufferPDPP.data(intsIndexesPDPP(119));

    t_x_yz_y_z = intsBufferPDPP.data(intsIndexesPDPP(120));

    t_x_yz_y_y = intsBufferPDPP.data(intsIndexesPDPP(121));

    t_x_yz_y_x = intsBufferPDPP.data(intsIndexesPDPP(122));

    t_x_yz_x_z = intsBufferPDPP.data(intsIndexesPDPP(123));

    t_x_yz_x_y = intsBufferPDPP.data(intsIndexesPDPP(124));

    t_x_yz_x_x = intsBufferPDPP.data(intsIndexesPDPP(125));

    t_x_yy_z_z = intsBufferPDPP.data(intsIndexesPDPP(126));

    t_x_yy_z_y = intsBufferPDPP.data(intsIndexesPDPP(127));

    t_x_yy_z_x = intsBufferPDPP.data(intsIndexesPDPP(128));

    t_x_yy_y_z = intsBufferPDPP.data(intsIndexesPDPP(129));

    t_x_yy_y_y = intsBufferPDPP.data(intsIndexesPDPP(130));

    t_x_yy_y_x = intsBufferPDPP.data(intsIndexesPDPP(131));

    t_x_yy_x_z = intsBufferPDPP.data(intsIndexesPDPP(132));

    t_x_yy_x_y = intsBufferPDPP.data(intsIndexesPDPP(133));

    t_x_yy_x_x = intsBufferPDPP.data(intsIndexesPDPP(134));

    t_x_xz_z_z = intsBufferPDPP.data(intsIndexesPDPP(135));

    t_x_xz_z_y = intsBufferPDPP.data(intsIndexesPDPP(136));

    t_x_xz_z_x = intsBufferPDPP.data(intsIndexesPDPP(137));

    t_x_xz_y_z = intsBufferPDPP.data(intsIndexesPDPP(138));

    t_x_xz_y_y = intsBufferPDPP.data(intsIndexesPDPP(139));

    t_x_xz_y_x = intsBufferPDPP.data(intsIndexesPDPP(140));

    t_x_xz_x_z = intsBufferPDPP.data(intsIndexesPDPP(141));

    t_x_xz_x_y = intsBufferPDPP.data(intsIndexesPDPP(142));

    t_x_xz_x_x = intsBufferPDPP.data(intsIndexesPDPP(143));

    t_x_xy_z_z = intsBufferPDPP.data(intsIndexesPDPP(144));

    t_x_xy_z_y = intsBufferPDPP.data(intsIndexesPDPP(145));

    t_x_xy_z_x = intsBufferPDPP.data(intsIndexesPDPP(146));

    t_x_xy_y_z = intsBufferPDPP.data(intsIndexesPDPP(147));

    t_x_xy_y_y = intsBufferPDPP.data(intsIndexesPDPP(148));

    t_x_xy_y_x = intsBufferPDPP.data(intsIndexesPDPP(149));

    t_x_xy_x_z = intsBufferPDPP.data(intsIndexesPDPP(150));

    t_x_xy_x_y = intsBufferPDPP.data(intsIndexesPDPP(151));

    t_x_xy_x_x = intsBufferPDPP.data(intsIndexesPDPP(152));

    t_x_xx_z_z = intsBufferPDPP.data(intsIndexesPDPP(153));

    t_x_xx_z_y = intsBufferPDPP.data(intsIndexesPDPP(154));

    t_x_xx_z_x = intsBufferPDPP.data(intsIndexesPDPP(155));

    t_x_xx_y_z = intsBufferPDPP.data(intsIndexesPDPP(156));

    t_x_xx_y_y = intsBufferPDPP.data(intsIndexesPDPP(157));

    t_x_xx_y_x = intsBufferPDPP.data(intsIndexesPDPP(158));

    t_x_xx_x_z = intsBufferPDPP.data(intsIndexesPDPP(159));

    t_x_xx_x_y = intsBufferPDPP.data(intsIndexesPDPP(160));

    t_x_xx_x_x = intsBufferPDPP.data(intsIndexesPDPP(161));

    // set up (PFPP) integral components

    t_z_zzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(0));

    t_z_zzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(1));

    t_z_zzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(2));

    t_z_zzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(3));

    t_z_zzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(4));

    t_z_zzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(5));

    t_z_zzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(6));

    t_z_zzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(7));

    t_z_zzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(8));

    t_z_yzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(9));

    t_z_yzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(10));

    t_z_yzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(11));

    t_z_yzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(12));

    t_z_yzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(13));

    t_z_yzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(14));

    t_z_yzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(15));

    t_z_yzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(16));

    t_z_yzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(17));

    t_z_yyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(18));

    t_z_yyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(19));

    t_z_yyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(20));

    t_z_yyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(21));

    t_z_yyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(22));

    t_z_yyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(23));

    t_z_yyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(24));

    t_z_yyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(25));

    t_z_yyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(26));

    t_z_xzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(27));

    t_z_xzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(28));

    t_z_xzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(29));

    t_z_xzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(30));

    t_z_xzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(31));

    t_z_xzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(32));

    t_z_xzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(33));

    t_z_xzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(34));

    t_z_xzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(35));

    t_z_xyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(36));

    t_z_xyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(37));

    t_z_xyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(38));

    t_z_xyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(39));

    t_z_xyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(40));

    t_z_xyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(41));

    t_z_xyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(42));

    t_z_xyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(43));

    t_z_xyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(44));

    t_z_xxz_z_z = intsBufferPFPP.data(intsIndexesPFPP(45));

    t_z_xxz_z_y = intsBufferPFPP.data(intsIndexesPFPP(46));

    t_z_xxz_z_x = intsBufferPFPP.data(intsIndexesPFPP(47));

    t_z_xxz_y_z = intsBufferPFPP.data(intsIndexesPFPP(48));

    t_z_xxz_y_y = intsBufferPFPP.data(intsIndexesPFPP(49));

    t_z_xxz_y_x = intsBufferPFPP.data(intsIndexesPFPP(50));

    t_z_xxz_x_z = intsBufferPFPP.data(intsIndexesPFPP(51));

    t_z_xxz_x_y = intsBufferPFPP.data(intsIndexesPFPP(52));

    t_z_xxz_x_x = intsBufferPFPP.data(intsIndexesPFPP(53));

    t_y_zzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(54));

    t_y_zzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(55));

    t_y_zzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(56));

    t_y_zzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(57));

    t_y_zzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(58));

    t_y_zzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(59));

    t_y_zzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(60));

    t_y_zzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(61));

    t_y_zzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(62));

    t_y_yzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(63));

    t_y_yzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(64));

    t_y_yzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(65));

    t_y_yzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(66));

    t_y_yzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(67));

    t_y_yzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(68));

    t_y_yzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(69));

    t_y_yzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(70));

    t_y_yzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(71));

    t_y_yyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(72));

    t_y_yyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(73));

    t_y_yyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(74));

    t_y_yyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(75));

    t_y_yyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(76));

    t_y_yyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(77));

    t_y_yyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(78));

    t_y_yyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(79));

    t_y_yyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(80));

    t_y_yyy_z_z = intsBufferPFPP.data(intsIndexesPFPP(81));

    t_y_yyy_z_y = intsBufferPFPP.data(intsIndexesPFPP(82));

    t_y_yyy_z_x = intsBufferPFPP.data(intsIndexesPFPP(83));

    t_y_yyy_y_z = intsBufferPFPP.data(intsIndexesPFPP(84));

    t_y_yyy_y_y = intsBufferPFPP.data(intsIndexesPFPP(85));

    t_y_yyy_y_x = intsBufferPFPP.data(intsIndexesPFPP(86));

    t_y_yyy_x_z = intsBufferPFPP.data(intsIndexesPFPP(87));

    t_y_yyy_x_y = intsBufferPFPP.data(intsIndexesPFPP(88));

    t_y_yyy_x_x = intsBufferPFPP.data(intsIndexesPFPP(89));

    t_y_xzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(90));

    t_y_xzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(91));

    t_y_xzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(92));

    t_y_xzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(93));

    t_y_xzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(94));

    t_y_xzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(95));

    t_y_xzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(96));

    t_y_xzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(97));

    t_y_xzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(98));

    t_y_xyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(99));

    t_y_xyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(100));

    t_y_xyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(101));

    t_y_xyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(102));

    t_y_xyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(103));

    t_y_xyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(104));

    t_y_xyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(105));

    t_y_xyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(106));

    t_y_xyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(107));

    t_y_xyy_z_z = intsBufferPFPP.data(intsIndexesPFPP(108));

    t_y_xyy_z_y = intsBufferPFPP.data(intsIndexesPFPP(109));

    t_y_xyy_z_x = intsBufferPFPP.data(intsIndexesPFPP(110));

    t_y_xyy_y_z = intsBufferPFPP.data(intsIndexesPFPP(111));

    t_y_xyy_y_y = intsBufferPFPP.data(intsIndexesPFPP(112));

    t_y_xyy_y_x = intsBufferPFPP.data(intsIndexesPFPP(113));

    t_y_xyy_x_z = intsBufferPFPP.data(intsIndexesPFPP(114));

    t_y_xyy_x_y = intsBufferPFPP.data(intsIndexesPFPP(115));

    t_y_xyy_x_x = intsBufferPFPP.data(intsIndexesPFPP(116));

    t_y_xxz_z_z = intsBufferPFPP.data(intsIndexesPFPP(117));

    t_y_xxz_z_y = intsBufferPFPP.data(intsIndexesPFPP(118));

    t_y_xxz_z_x = intsBufferPFPP.data(intsIndexesPFPP(119));

    t_y_xxz_y_z = intsBufferPFPP.data(intsIndexesPFPP(120));

    t_y_xxz_y_y = intsBufferPFPP.data(intsIndexesPFPP(121));

    t_y_xxz_y_x = intsBufferPFPP.data(intsIndexesPFPP(122));

    t_y_xxz_x_z = intsBufferPFPP.data(intsIndexesPFPP(123));

    t_y_xxz_x_y = intsBufferPFPP.data(intsIndexesPFPP(124));

    t_y_xxz_x_x = intsBufferPFPP.data(intsIndexesPFPP(125));

    t_y_xxy_z_z = intsBufferPFPP.data(intsIndexesPFPP(126));

    t_y_xxy_z_y = intsBufferPFPP.data(intsIndexesPFPP(127));

    t_y_xxy_z_x = intsBufferPFPP.data(intsIndexesPFPP(128));

    t_y_xxy_y_z = intsBufferPFPP.data(intsIndexesPFPP(129));

    t_y_xxy_y_y = intsBufferPFPP.data(intsIndexesPFPP(130));

    t_y_xxy_y_x = intsBufferPFPP.data(intsIndexesPFPP(131));

    t_y_xxy_x_z = intsBufferPFPP.data(intsIndexesPFPP(132));

    t_y_xxy_x_y = intsBufferPFPP.data(intsIndexesPFPP(133));

    t_y_xxy_x_x = intsBufferPFPP.data(intsIndexesPFPP(134));

    t_x_zzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(135));

    t_x_zzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(136));

    t_x_zzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(137));

    t_x_zzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(138));

    t_x_zzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(139));

    t_x_zzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(140));

    t_x_zzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(141));

    t_x_zzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(142));

    t_x_zzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(143));

    t_x_yzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(144));

    t_x_yzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(145));

    t_x_yzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(146));

    t_x_yzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(147));

    t_x_yzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(148));

    t_x_yzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(149));

    t_x_yzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(150));

    t_x_yzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(151));

    t_x_yzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(152));

    t_x_yyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(153));

    t_x_yyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(154));

    t_x_yyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(155));

    t_x_yyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(156));

    t_x_yyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(157));

    t_x_yyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(158));

    t_x_yyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(159));

    t_x_yyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(160));

    t_x_yyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(161));

    t_x_yyy_z_z = intsBufferPFPP.data(intsIndexesPFPP(162));

    t_x_yyy_z_y = intsBufferPFPP.data(intsIndexesPFPP(163));

    t_x_yyy_z_x = intsBufferPFPP.data(intsIndexesPFPP(164));

    t_x_yyy_y_z = intsBufferPFPP.data(intsIndexesPFPP(165));

    t_x_yyy_y_y = intsBufferPFPP.data(intsIndexesPFPP(166));

    t_x_yyy_y_x = intsBufferPFPP.data(intsIndexesPFPP(167));

    t_x_yyy_x_z = intsBufferPFPP.data(intsIndexesPFPP(168));

    t_x_yyy_x_y = intsBufferPFPP.data(intsIndexesPFPP(169));

    t_x_yyy_x_x = intsBufferPFPP.data(intsIndexesPFPP(170));

    t_x_xzz_z_z = intsBufferPFPP.data(intsIndexesPFPP(171));

    t_x_xzz_z_y = intsBufferPFPP.data(intsIndexesPFPP(172));

    t_x_xzz_z_x = intsBufferPFPP.data(intsIndexesPFPP(173));

    t_x_xzz_y_z = intsBufferPFPP.data(intsIndexesPFPP(174));

    t_x_xzz_y_y = intsBufferPFPP.data(intsIndexesPFPP(175));

    t_x_xzz_y_x = intsBufferPFPP.data(intsIndexesPFPP(176));

    t_x_xzz_x_z = intsBufferPFPP.data(intsIndexesPFPP(177));

    t_x_xzz_x_y = intsBufferPFPP.data(intsIndexesPFPP(178));

    t_x_xzz_x_x = intsBufferPFPP.data(intsIndexesPFPP(179));

    t_x_xyz_z_z = intsBufferPFPP.data(intsIndexesPFPP(180));

    t_x_xyz_z_y = intsBufferPFPP.data(intsIndexesPFPP(181));

    t_x_xyz_z_x = intsBufferPFPP.data(intsIndexesPFPP(182));

    t_x_xyz_y_z = intsBufferPFPP.data(intsIndexesPFPP(183));

    t_x_xyz_y_y = intsBufferPFPP.data(intsIndexesPFPP(184));

    t_x_xyz_y_x = intsBufferPFPP.data(intsIndexesPFPP(185));

    t_x_xyz_x_z = intsBufferPFPP.data(intsIndexesPFPP(186));

    t_x_xyz_x_y = intsBufferPFPP.data(intsIndexesPFPP(187));

    t_x_xyz_x_x = intsBufferPFPP.data(intsIndexesPFPP(188));

    t_x_xyy_z_z = intsBufferPFPP.data(intsIndexesPFPP(189));

    t_x_xyy_z_y = intsBufferPFPP.data(intsIndexesPFPP(190));

    t_x_xyy_z_x = intsBufferPFPP.data(intsIndexesPFPP(191));

    t_x_xyy_y_z = intsBufferPFPP.data(intsIndexesPFPP(192));

    t_x_xyy_y_y = intsBufferPFPP.data(intsIndexesPFPP(193));

    t_x_xyy_y_x = intsBufferPFPP.data(intsIndexesPFPP(194));

    t_x_xyy_x_z = intsBufferPFPP.data(intsIndexesPFPP(195));

    t_x_xyy_x_y = intsBufferPFPP.data(intsIndexesPFPP(196));

    t_x_xyy_x_x = intsBufferPFPP.data(intsIndexesPFPP(197));

    t_x_xxz_z_z = intsBufferPFPP.data(intsIndexesPFPP(198));

    t_x_xxz_z_y = intsBufferPFPP.data(intsIndexesPFPP(199));

    t_x_xxz_z_x = intsBufferPFPP.data(intsIndexesPFPP(200));

    t_x_xxz_y_z = intsBufferPFPP.data(intsIndexesPFPP(201));

    t_x_xxz_y_y = intsBufferPFPP.data(intsIndexesPFPP(202));

    t_x_xxz_y_x = intsBufferPFPP.data(intsIndexesPFPP(203));

    t_x_xxz_x_z = intsBufferPFPP.data(intsIndexesPFPP(204));

    t_x_xxz_x_y = intsBufferPFPP.data(intsIndexesPFPP(205));

    t_x_xxz_x_x = intsBufferPFPP.data(intsIndexesPFPP(206));

    t_x_xxy_z_z = intsBufferPFPP.data(intsIndexesPFPP(207));

    t_x_xxy_z_y = intsBufferPFPP.data(intsIndexesPFPP(208));

    t_x_xxy_z_x = intsBufferPFPP.data(intsIndexesPFPP(209));

    t_x_xxy_y_z = intsBufferPFPP.data(intsIndexesPFPP(210));

    t_x_xxy_y_y = intsBufferPFPP.data(intsIndexesPFPP(211));

    t_x_xxy_y_x = intsBufferPFPP.data(intsIndexesPFPP(212));

    t_x_xxy_x_z = intsBufferPFPP.data(intsIndexesPFPP(213));

    t_x_xxy_x_y = intsBufferPFPP.data(intsIndexesPFPP(214));

    t_x_xxy_x_x = intsBufferPFPP.data(intsIndexesPFPP(215));

    t_x_xxx_z_z = intsBufferPFPP.data(intsIndexesPFPP(216));

    t_x_xxx_z_y = intsBufferPFPP.data(intsIndexesPFPP(217));

    t_x_xxx_z_x = intsBufferPFPP.data(intsIndexesPFPP(218));

    t_x_xxx_y_z = intsBufferPFPP.data(intsIndexesPFPP(219));

    t_x_xxx_y_y = intsBufferPFPP.data(intsIndexesPFPP(220));

    t_x_xxx_y_x = intsBufferPFPP.data(intsIndexesPFPP(221));

    t_x_xxx_x_z = intsBufferPFPP.data(intsIndexesPFPP(222));

    t_x_xxx_x_y = intsBufferPFPP.data(intsIndexesPFPP(223));

    t_x_xxx_x_x = intsBufferPFPP.data(intsIndexesPFPP(224));

    #pragma omp simd align(rab_z, t_z_xz_x_x, t_z_xz_x_y, t_z_xz_x_z, t_z_xz_y_x, t_z_xz_y_y,\
                           t_z_xz_y_z, t_z_xz_z_x, t_z_xz_z_y, t_z_xz_z_z, t_z_xzz_x_x,\
                           t_z_xzz_x_y, t_z_xzz_x_z, t_z_xzz_y_x, t_z_xzz_y_y, t_z_xzz_y_z,\
                           t_z_xzz_z_x, t_z_xzz_z_y, t_z_xzz_z_z, t_z_yy_x_x, t_z_yy_x_y,\
                           t_z_yy_x_z, t_z_yy_y_x, t_z_yy_y_y, t_z_yy_y_z, t_z_yy_z_x,\
                           t_z_yy_z_y, t_z_yy_z_z, t_z_yyz_x_x, t_z_yyz_x_y, t_z_yyz_x_z,\
                           t_z_yyz_y_x, t_z_yyz_y_y, t_z_yyz_y_z, t_z_yyz_z_x, t_z_yyz_z_y,\
                           t_z_yyz_z_z, t_z_yz_x_x, t_z_yz_x_y, t_z_yz_x_z, t_z_yz_y_x,\
                           t_z_yz_y_y, t_z_yz_y_z, t_z_yz_z_x, t_z_yz_z_y, t_z_yz_z_z,\
                           t_z_yzz_x_x, t_z_yzz_x_y, t_z_yzz_x_z, t_z_yzz_y_x, t_z_yzz_y_y,\
                           t_z_yzz_y_z, t_z_yzz_z_x, t_z_yzz_z_y, t_z_yzz_z_z, t_z_zz_x_x,\
                           t_z_zz_x_y, t_z_zz_x_z, t_z_zz_y_x, t_z_zz_y_y, t_z_zz_y_z,\
                           t_z_zz_z_x, t_z_zz_z_y, t_z_zz_z_z, t_z_zzz_x_x, t_z_zzz_x_y,\
                           t_z_zzz_x_z, t_z_zzz_y_x, t_z_zzz_y_y, t_z_zzz_y_z, t_z_zzz_z_x,\
                           t_z_zzz_z_y, t_z_zzz_z_z, t_zz_xz_x_x, t_zz_xz_x_y, t_zz_xz_x_z,\
                           t_zz_xz_y_x, t_zz_xz_y_y, t_zz_xz_y_z, t_zz_xz_z_x, t_zz_xz_z_y,\
                           t_zz_xz_z_z, t_zz_yy_x_x, t_zz_yy_x_y, t_zz_yy_x_z, t_zz_yy_y_x,\
                           t_zz_yy_y_y, t_zz_yy_y_z, t_zz_yy_z_x, t_zz_yy_z_y, t_zz_yy_z_z,\
                           t_zz_yz_x_x, t_zz_yz_x_y, t_zz_yz_x_z, t_zz_yz_y_x, t_zz_yz_y_y,\
                           t_zz_yz_y_z, t_zz_yz_z_x, t_zz_yz_z_y, t_zz_yz_z_z, t_zz_zz_x_x,\
                           t_zz_zz_x_y, t_zz_zz_x_z, t_zz_zz_y_x, t_zz_zz_y_y, t_zz_zz_y_z,\
                           t_zz_zz_z_x, t_zz_zz_z_y, t_zz_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_z_z[i] = t_z_zzz_z_z[i] - rab_z[i] * t_z_zz_z_z[i];

        t_zz_zz_z_y[i] = t_z_zzz_z_y[i] - rab_z[i] * t_z_zz_z_y[i];

        t_zz_zz_z_x[i] = t_z_zzz_z_x[i] - rab_z[i] * t_z_zz_z_x[i];

        t_zz_zz_y_z[i] = t_z_zzz_y_z[i] - rab_z[i] * t_z_zz_y_z[i];

        t_zz_zz_y_y[i] = t_z_zzz_y_y[i] - rab_z[i] * t_z_zz_y_y[i];

        t_zz_zz_y_x[i] = t_z_zzz_y_x[i] - rab_z[i] * t_z_zz_y_x[i];

        t_zz_zz_x_z[i] = t_z_zzz_x_z[i] - rab_z[i] * t_z_zz_x_z[i];

        t_zz_zz_x_y[i] = t_z_zzz_x_y[i] - rab_z[i] * t_z_zz_x_y[i];

        t_zz_zz_x_x[i] = t_z_zzz_x_x[i] - rab_z[i] * t_z_zz_x_x[i];

        t_zz_yz_z_z[i] = t_z_yzz_z_z[i] - rab_z[i] * t_z_yz_z_z[i];

        t_zz_yz_z_y[i] = t_z_yzz_z_y[i] - rab_z[i] * t_z_yz_z_y[i];

        t_zz_yz_z_x[i] = t_z_yzz_z_x[i] - rab_z[i] * t_z_yz_z_x[i];

        t_zz_yz_y_z[i] = t_z_yzz_y_z[i] - rab_z[i] * t_z_yz_y_z[i];

        t_zz_yz_y_y[i] = t_z_yzz_y_y[i] - rab_z[i] * t_z_yz_y_y[i];

        t_zz_yz_y_x[i] = t_z_yzz_y_x[i] - rab_z[i] * t_z_yz_y_x[i];

        t_zz_yz_x_z[i] = t_z_yzz_x_z[i] - rab_z[i] * t_z_yz_x_z[i];

        t_zz_yz_x_y[i] = t_z_yzz_x_y[i] - rab_z[i] * t_z_yz_x_y[i];

        t_zz_yz_x_x[i] = t_z_yzz_x_x[i] - rab_z[i] * t_z_yz_x_x[i];

        t_zz_yy_z_z[i] = t_z_yyz_z_z[i] - rab_z[i] * t_z_yy_z_z[i];

        t_zz_yy_z_y[i] = t_z_yyz_z_y[i] - rab_z[i] * t_z_yy_z_y[i];

        t_zz_yy_z_x[i] = t_z_yyz_z_x[i] - rab_z[i] * t_z_yy_z_x[i];

        t_zz_yy_y_z[i] = t_z_yyz_y_z[i] - rab_z[i] * t_z_yy_y_z[i];

        t_zz_yy_y_y[i] = t_z_yyz_y_y[i] - rab_z[i] * t_z_yy_y_y[i];

        t_zz_yy_y_x[i] = t_z_yyz_y_x[i] - rab_z[i] * t_z_yy_y_x[i];

        t_zz_yy_x_z[i] = t_z_yyz_x_z[i] - rab_z[i] * t_z_yy_x_z[i];

        t_zz_yy_x_y[i] = t_z_yyz_x_y[i] - rab_z[i] * t_z_yy_x_y[i];

        t_zz_yy_x_x[i] = t_z_yyz_x_x[i] - rab_z[i] * t_z_yy_x_x[i];

        t_zz_xz_z_z[i] = t_z_xzz_z_z[i] - rab_z[i] * t_z_xz_z_z[i];

        t_zz_xz_z_y[i] = t_z_xzz_z_y[i] - rab_z[i] * t_z_xz_z_y[i];

        t_zz_xz_z_x[i] = t_z_xzz_z_x[i] - rab_z[i] * t_z_xz_z_x[i];

        t_zz_xz_y_z[i] = t_z_xzz_y_z[i] - rab_z[i] * t_z_xz_y_z[i];

        t_zz_xz_y_y[i] = t_z_xzz_y_y[i] - rab_z[i] * t_z_xz_y_y[i];

        t_zz_xz_y_x[i] = t_z_xzz_y_x[i] - rab_z[i] * t_z_xz_y_x[i];

        t_zz_xz_x_z[i] = t_z_xzz_x_z[i] - rab_z[i] * t_z_xz_x_z[i];

        t_zz_xz_x_y[i] = t_z_xzz_x_y[i] - rab_z[i] * t_z_xz_x_y[i];

        t_zz_xz_x_x[i] = t_z_xzz_x_x[i] - rab_z[i] * t_z_xz_x_x[i];
    }

    #pragma omp simd align(rab_z, t_y_yz_x_x, t_y_yz_x_y, t_y_yz_x_z, t_y_yz_y_x, t_y_yz_y_y,\
                           t_y_yz_y_z, t_y_yz_z_x, t_y_yz_z_y, t_y_yz_z_z, t_y_yzz_x_x,\
                           t_y_yzz_x_y, t_y_yzz_x_z, t_y_yzz_y_x, t_y_yzz_y_y, t_y_yzz_y_z,\
                           t_y_yzz_z_x, t_y_yzz_z_y, t_y_yzz_z_z, t_y_zz_x_x, t_y_zz_x_y,\
                           t_y_zz_x_z, t_y_zz_y_x, t_y_zz_y_y, t_y_zz_y_z, t_y_zz_z_x,\
                           t_y_zz_z_y, t_y_zz_z_z, t_y_zzz_x_x, t_y_zzz_x_y, t_y_zzz_x_z,\
                           t_y_zzz_y_x, t_y_zzz_y_y, t_y_zzz_y_z, t_y_zzz_z_x, t_y_zzz_z_y,\
                           t_y_zzz_z_z, t_yz_yz_x_x, t_yz_yz_x_y, t_yz_yz_x_z, t_yz_yz_y_x,\
                           t_yz_yz_y_y, t_yz_yz_y_z, t_yz_yz_z_x, t_yz_yz_z_y, t_yz_yz_z_z,\
                           t_yz_zz_x_x, t_yz_zz_x_y, t_yz_zz_x_z, t_yz_zz_y_x, t_yz_zz_y_y,\
                           t_yz_zz_y_z, t_yz_zz_z_x, t_yz_zz_z_y, t_yz_zz_z_z, t_z_xx_x_x,\
                           t_z_xx_x_y, t_z_xx_x_z, t_z_xx_y_x, t_z_xx_y_y, t_z_xx_y_z,\
                           t_z_xx_z_x, t_z_xx_z_y, t_z_xx_z_z, t_z_xxz_x_x, t_z_xxz_x_y,\
                           t_z_xxz_x_z, t_z_xxz_y_x, t_z_xxz_y_y, t_z_xxz_y_z, t_z_xxz_z_x,\
                           t_z_xxz_z_y, t_z_xxz_z_z, t_z_xy_x_x, t_z_xy_x_y, t_z_xy_x_z,\
                           t_z_xy_y_x, t_z_xy_y_y, t_z_xy_y_z, t_z_xy_z_x, t_z_xy_z_y,\
                           t_z_xy_z_z, t_z_xyz_x_x, t_z_xyz_x_y, t_z_xyz_x_z, t_z_xyz_y_x,\
                           t_z_xyz_y_y, t_z_xyz_y_z, t_z_xyz_z_x, t_z_xyz_z_y, t_z_xyz_z_z,\
                           t_zz_xx_x_x, t_zz_xx_x_y, t_zz_xx_x_z, t_zz_xx_y_x, t_zz_xx_y_y,\
                           t_zz_xx_y_z, t_zz_xx_z_x, t_zz_xx_z_y, t_zz_xx_z_z, t_zz_xy_x_x,\
                           t_zz_xy_x_y, t_zz_xy_x_z, t_zz_xy_y_x, t_zz_xy_y_y, t_zz_xy_y_z,\
                           t_zz_xy_z_x, t_zz_xy_z_y, t_zz_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_xy_z_z[i] = t_z_xyz_z_z[i] - rab_z[i] * t_z_xy_z_z[i];

        t_zz_xy_z_y[i] = t_z_xyz_z_y[i] - rab_z[i] * t_z_xy_z_y[i];

        t_zz_xy_z_x[i] = t_z_xyz_z_x[i] - rab_z[i] * t_z_xy_z_x[i];

        t_zz_xy_y_z[i] = t_z_xyz_y_z[i] - rab_z[i] * t_z_xy_y_z[i];

        t_zz_xy_y_y[i] = t_z_xyz_y_y[i] - rab_z[i] * t_z_xy_y_y[i];

        t_zz_xy_y_x[i] = t_z_xyz_y_x[i] - rab_z[i] * t_z_xy_y_x[i];

        t_zz_xy_x_z[i] = t_z_xyz_x_z[i] - rab_z[i] * t_z_xy_x_z[i];

        t_zz_xy_x_y[i] = t_z_xyz_x_y[i] - rab_z[i] * t_z_xy_x_y[i];

        t_zz_xy_x_x[i] = t_z_xyz_x_x[i] - rab_z[i] * t_z_xy_x_x[i];

        t_zz_xx_z_z[i] = t_z_xxz_z_z[i] - rab_z[i] * t_z_xx_z_z[i];

        t_zz_xx_z_y[i] = t_z_xxz_z_y[i] - rab_z[i] * t_z_xx_z_y[i];

        t_zz_xx_z_x[i] = t_z_xxz_z_x[i] - rab_z[i] * t_z_xx_z_x[i];

        t_zz_xx_y_z[i] = t_z_xxz_y_z[i] - rab_z[i] * t_z_xx_y_z[i];

        t_zz_xx_y_y[i] = t_z_xxz_y_y[i] - rab_z[i] * t_z_xx_y_y[i];

        t_zz_xx_y_x[i] = t_z_xxz_y_x[i] - rab_z[i] * t_z_xx_y_x[i];

        t_zz_xx_x_z[i] = t_z_xxz_x_z[i] - rab_z[i] * t_z_xx_x_z[i];

        t_zz_xx_x_y[i] = t_z_xxz_x_y[i] - rab_z[i] * t_z_xx_x_y[i];

        t_zz_xx_x_x[i] = t_z_xxz_x_x[i] - rab_z[i] * t_z_xx_x_x[i];

        t_yz_zz_z_z[i] = t_y_zzz_z_z[i] - rab_z[i] * t_y_zz_z_z[i];

        t_yz_zz_z_y[i] = t_y_zzz_z_y[i] - rab_z[i] * t_y_zz_z_y[i];

        t_yz_zz_z_x[i] = t_y_zzz_z_x[i] - rab_z[i] * t_y_zz_z_x[i];

        t_yz_zz_y_z[i] = t_y_zzz_y_z[i] - rab_z[i] * t_y_zz_y_z[i];

        t_yz_zz_y_y[i] = t_y_zzz_y_y[i] - rab_z[i] * t_y_zz_y_y[i];

        t_yz_zz_y_x[i] = t_y_zzz_y_x[i] - rab_z[i] * t_y_zz_y_x[i];

        t_yz_zz_x_z[i] = t_y_zzz_x_z[i] - rab_z[i] * t_y_zz_x_z[i];

        t_yz_zz_x_y[i] = t_y_zzz_x_y[i] - rab_z[i] * t_y_zz_x_y[i];

        t_yz_zz_x_x[i] = t_y_zzz_x_x[i] - rab_z[i] * t_y_zz_x_x[i];

        t_yz_yz_z_z[i] = t_y_yzz_z_z[i] - rab_z[i] * t_y_yz_z_z[i];

        t_yz_yz_z_y[i] = t_y_yzz_z_y[i] - rab_z[i] * t_y_yz_z_y[i];

        t_yz_yz_z_x[i] = t_y_yzz_z_x[i] - rab_z[i] * t_y_yz_z_x[i];

        t_yz_yz_y_z[i] = t_y_yzz_y_z[i] - rab_z[i] * t_y_yz_y_z[i];

        t_yz_yz_y_y[i] = t_y_yzz_y_y[i] - rab_z[i] * t_y_yz_y_y[i];

        t_yz_yz_y_x[i] = t_y_yzz_y_x[i] - rab_z[i] * t_y_yz_y_x[i];

        t_yz_yz_x_z[i] = t_y_yzz_x_z[i] - rab_z[i] * t_y_yz_x_z[i];

        t_yz_yz_x_y[i] = t_y_yzz_x_y[i] - rab_z[i] * t_y_yz_x_y[i];

        t_yz_yz_x_x[i] = t_y_yzz_x_x[i] - rab_z[i] * t_y_yz_x_x[i];
    }

    #pragma omp simd align(rab_z, t_y_xx_x_x, t_y_xx_x_y, t_y_xx_x_z, t_y_xx_y_x, t_y_xx_y_y,\
                           t_y_xx_y_z, t_y_xx_z_x, t_y_xx_z_y, t_y_xx_z_z, t_y_xxz_x_x,\
                           t_y_xxz_x_y, t_y_xxz_x_z, t_y_xxz_y_x, t_y_xxz_y_y, t_y_xxz_y_z,\
                           t_y_xxz_z_x, t_y_xxz_z_y, t_y_xxz_z_z, t_y_xy_x_x, t_y_xy_x_y,\
                           t_y_xy_x_z, t_y_xy_y_x, t_y_xy_y_y, t_y_xy_y_z, t_y_xy_z_x,\
                           t_y_xy_z_y, t_y_xy_z_z, t_y_xyz_x_x, t_y_xyz_x_y, t_y_xyz_x_z,\
                           t_y_xyz_y_x, t_y_xyz_y_y, t_y_xyz_y_z, t_y_xyz_z_x, t_y_xyz_z_y,\
                           t_y_xyz_z_z, t_y_xz_x_x, t_y_xz_x_y, t_y_xz_x_z, t_y_xz_y_x,\
                           t_y_xz_y_y, t_y_xz_y_z, t_y_xz_z_x, t_y_xz_z_y, t_y_xz_z_z,\
                           t_y_xzz_x_x, t_y_xzz_x_y, t_y_xzz_x_z, t_y_xzz_y_x, t_y_xzz_y_y,\
                           t_y_xzz_y_z, t_y_xzz_z_x, t_y_xzz_z_y, t_y_xzz_z_z, t_y_yy_x_x,\
                           t_y_yy_x_y, t_y_yy_x_z, t_y_yy_y_x, t_y_yy_y_y, t_y_yy_y_z,\
                           t_y_yy_z_x, t_y_yy_z_y, t_y_yy_z_z, t_y_yyz_x_x, t_y_yyz_x_y,\
                           t_y_yyz_x_z, t_y_yyz_y_x, t_y_yyz_y_y, t_y_yyz_y_z, t_y_yyz_z_x,\
                           t_y_yyz_z_y, t_y_yyz_z_z, t_yz_xx_x_x, t_yz_xx_x_y, t_yz_xx_x_z,\
                           t_yz_xx_y_x, t_yz_xx_y_y, t_yz_xx_y_z, t_yz_xx_z_x, t_yz_xx_z_y,\
                           t_yz_xx_z_z, t_yz_xy_x_x, t_yz_xy_x_y, t_yz_xy_x_z, t_yz_xy_y_x,\
                           t_yz_xy_y_y, t_yz_xy_y_z, t_yz_xy_z_x, t_yz_xy_z_y, t_yz_xy_z_z,\
                           t_yz_xz_x_x, t_yz_xz_x_y, t_yz_xz_x_z, t_yz_xz_y_x, t_yz_xz_y_y,\
                           t_yz_xz_y_z, t_yz_xz_z_x, t_yz_xz_z_y, t_yz_xz_z_z, t_yz_yy_x_x,\
                           t_yz_yy_x_y, t_yz_yy_x_z, t_yz_yy_y_x, t_yz_yy_y_y, t_yz_yy_y_z,\
                           t_yz_yy_z_x, t_yz_yy_z_y, t_yz_yy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yz_yy_z_z[i] = t_y_yyz_z_z[i] - rab_z[i] * t_y_yy_z_z[i];

        t_yz_yy_z_y[i] = t_y_yyz_z_y[i] - rab_z[i] * t_y_yy_z_y[i];

        t_yz_yy_z_x[i] = t_y_yyz_z_x[i] - rab_z[i] * t_y_yy_z_x[i];

        t_yz_yy_y_z[i] = t_y_yyz_y_z[i] - rab_z[i] * t_y_yy_y_z[i];

        t_yz_yy_y_y[i] = t_y_yyz_y_y[i] - rab_z[i] * t_y_yy_y_y[i];

        t_yz_yy_y_x[i] = t_y_yyz_y_x[i] - rab_z[i] * t_y_yy_y_x[i];

        t_yz_yy_x_z[i] = t_y_yyz_x_z[i] - rab_z[i] * t_y_yy_x_z[i];

        t_yz_yy_x_y[i] = t_y_yyz_x_y[i] - rab_z[i] * t_y_yy_x_y[i];

        t_yz_yy_x_x[i] = t_y_yyz_x_x[i] - rab_z[i] * t_y_yy_x_x[i];

        t_yz_xz_z_z[i] = t_y_xzz_z_z[i] - rab_z[i] * t_y_xz_z_z[i];

        t_yz_xz_z_y[i] = t_y_xzz_z_y[i] - rab_z[i] * t_y_xz_z_y[i];

        t_yz_xz_z_x[i] = t_y_xzz_z_x[i] - rab_z[i] * t_y_xz_z_x[i];

        t_yz_xz_y_z[i] = t_y_xzz_y_z[i] - rab_z[i] * t_y_xz_y_z[i];

        t_yz_xz_y_y[i] = t_y_xzz_y_y[i] - rab_z[i] * t_y_xz_y_y[i];

        t_yz_xz_y_x[i] = t_y_xzz_y_x[i] - rab_z[i] * t_y_xz_y_x[i];

        t_yz_xz_x_z[i] = t_y_xzz_x_z[i] - rab_z[i] * t_y_xz_x_z[i];

        t_yz_xz_x_y[i] = t_y_xzz_x_y[i] - rab_z[i] * t_y_xz_x_y[i];

        t_yz_xz_x_x[i] = t_y_xzz_x_x[i] - rab_z[i] * t_y_xz_x_x[i];

        t_yz_xy_z_z[i] = t_y_xyz_z_z[i] - rab_z[i] * t_y_xy_z_z[i];

        t_yz_xy_z_y[i] = t_y_xyz_z_y[i] - rab_z[i] * t_y_xy_z_y[i];

        t_yz_xy_z_x[i] = t_y_xyz_z_x[i] - rab_z[i] * t_y_xy_z_x[i];

        t_yz_xy_y_z[i] = t_y_xyz_y_z[i] - rab_z[i] * t_y_xy_y_z[i];

        t_yz_xy_y_y[i] = t_y_xyz_y_y[i] - rab_z[i] * t_y_xy_y_y[i];

        t_yz_xy_y_x[i] = t_y_xyz_y_x[i] - rab_z[i] * t_y_xy_y_x[i];

        t_yz_xy_x_z[i] = t_y_xyz_x_z[i] - rab_z[i] * t_y_xy_x_z[i];

        t_yz_xy_x_y[i] = t_y_xyz_x_y[i] - rab_z[i] * t_y_xy_x_y[i];

        t_yz_xy_x_x[i] = t_y_xyz_x_x[i] - rab_z[i] * t_y_xy_x_x[i];

        t_yz_xx_z_z[i] = t_y_xxz_z_z[i] - rab_z[i] * t_y_xx_z_z[i];

        t_yz_xx_z_y[i] = t_y_xxz_z_y[i] - rab_z[i] * t_y_xx_z_y[i];

        t_yz_xx_z_x[i] = t_y_xxz_z_x[i] - rab_z[i] * t_y_xx_z_x[i];

        t_yz_xx_y_z[i] = t_y_xxz_y_z[i] - rab_z[i] * t_y_xx_y_z[i];

        t_yz_xx_y_y[i] = t_y_xxz_y_y[i] - rab_z[i] * t_y_xx_y_y[i];

        t_yz_xx_y_x[i] = t_y_xxz_y_x[i] - rab_z[i] * t_y_xx_y_x[i];

        t_yz_xx_x_z[i] = t_y_xxz_x_z[i] - rab_z[i] * t_y_xx_x_z[i];

        t_yz_xx_x_y[i] = t_y_xxz_x_y[i] - rab_z[i] * t_y_xx_x_y[i];

        t_yz_xx_x_x[i] = t_y_xxz_x_x[i] - rab_z[i] * t_y_xx_x_x[i];
    }

    #pragma omp simd align(rab_y, t_y_xyz_x_x, t_y_xyz_x_y, t_y_xyz_x_z, t_y_xyz_y_x,\
                           t_y_xyz_y_y, t_y_xyz_y_z, t_y_xyz_z_x, t_y_xyz_z_y, t_y_xyz_z_z,\
                           t_y_xz_x_x, t_y_xz_x_y, t_y_xz_x_z, t_y_xz_y_x, t_y_xz_y_y,\
                           t_y_xz_y_z, t_y_xz_z_x, t_y_xz_z_y, t_y_xz_z_z, t_y_yy_x_x,\
                           t_y_yy_x_y, t_y_yy_x_z, t_y_yy_y_x, t_y_yy_y_y, t_y_yy_y_z,\
                           t_y_yy_z_x, t_y_yy_z_y, t_y_yy_z_z, t_y_yyy_x_x, t_y_yyy_x_y,\
                           t_y_yyy_x_z, t_y_yyy_y_x, t_y_yyy_y_y, t_y_yyy_y_z, t_y_yyy_z_x,\
                           t_y_yyy_z_y, t_y_yyy_z_z, t_y_yyz_x_x, t_y_yyz_x_y, t_y_yyz_x_z,\
                           t_y_yyz_y_x, t_y_yyz_y_y, t_y_yyz_y_z, t_y_yyz_z_x, t_y_yyz_z_y,\
                           t_y_yyz_z_z, t_y_yz_x_x, t_y_yz_x_y, t_y_yz_x_z, t_y_yz_y_x,\
                           t_y_yz_y_y, t_y_yz_y_z, t_y_yz_z_x, t_y_yz_z_y, t_y_yz_z_z,\
                           t_y_yzz_x_x, t_y_yzz_x_y, t_y_yzz_x_z, t_y_yzz_y_x, t_y_yzz_y_y,\
                           t_y_yzz_y_z, t_y_yzz_z_x, t_y_yzz_z_y, t_y_yzz_z_z, t_y_zz_x_x,\
                           t_y_zz_x_y, t_y_zz_x_z, t_y_zz_y_x, t_y_zz_y_y, t_y_zz_y_z,\
                           t_y_zz_z_x, t_y_zz_z_y, t_y_zz_z_z, t_yy_xz_x_x, t_yy_xz_x_y,\
                           t_yy_xz_x_z, t_yy_xz_y_x, t_yy_xz_y_y, t_yy_xz_y_z, t_yy_xz_z_x,\
                           t_yy_xz_z_y, t_yy_xz_z_z, t_yy_yy_x_x, t_yy_yy_x_y, t_yy_yy_x_z,\
                           t_yy_yy_y_x, t_yy_yy_y_y, t_yy_yy_y_z, t_yy_yy_z_x, t_yy_yy_z_y,\
                           t_yy_yy_z_z, t_yy_yz_x_x, t_yy_yz_x_y, t_yy_yz_x_z, t_yy_yz_y_x,\
                           t_yy_yz_y_y, t_yy_yz_y_z, t_yy_yz_z_x, t_yy_yz_z_y, t_yy_yz_z_z,\
                           t_yy_zz_x_x, t_yy_zz_x_y, t_yy_zz_x_z, t_yy_zz_y_x, t_yy_zz_y_y,\
                           t_yy_zz_y_z, t_yy_zz_z_x, t_yy_zz_z_y, t_yy_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_zz_z_z[i] = t_y_yzz_z_z[i] - rab_y[i] * t_y_zz_z_z[i];

        t_yy_zz_z_y[i] = t_y_yzz_z_y[i] - rab_y[i] * t_y_zz_z_y[i];

        t_yy_zz_z_x[i] = t_y_yzz_z_x[i] - rab_y[i] * t_y_zz_z_x[i];

        t_yy_zz_y_z[i] = t_y_yzz_y_z[i] - rab_y[i] * t_y_zz_y_z[i];

        t_yy_zz_y_y[i] = t_y_yzz_y_y[i] - rab_y[i] * t_y_zz_y_y[i];

        t_yy_zz_y_x[i] = t_y_yzz_y_x[i] - rab_y[i] * t_y_zz_y_x[i];

        t_yy_zz_x_z[i] = t_y_yzz_x_z[i] - rab_y[i] * t_y_zz_x_z[i];

        t_yy_zz_x_y[i] = t_y_yzz_x_y[i] - rab_y[i] * t_y_zz_x_y[i];

        t_yy_zz_x_x[i] = t_y_yzz_x_x[i] - rab_y[i] * t_y_zz_x_x[i];

        t_yy_yz_z_z[i] = t_y_yyz_z_z[i] - rab_y[i] * t_y_yz_z_z[i];

        t_yy_yz_z_y[i] = t_y_yyz_z_y[i] - rab_y[i] * t_y_yz_z_y[i];

        t_yy_yz_z_x[i] = t_y_yyz_z_x[i] - rab_y[i] * t_y_yz_z_x[i];

        t_yy_yz_y_z[i] = t_y_yyz_y_z[i] - rab_y[i] * t_y_yz_y_z[i];

        t_yy_yz_y_y[i] = t_y_yyz_y_y[i] - rab_y[i] * t_y_yz_y_y[i];

        t_yy_yz_y_x[i] = t_y_yyz_y_x[i] - rab_y[i] * t_y_yz_y_x[i];

        t_yy_yz_x_z[i] = t_y_yyz_x_z[i] - rab_y[i] * t_y_yz_x_z[i];

        t_yy_yz_x_y[i] = t_y_yyz_x_y[i] - rab_y[i] * t_y_yz_x_y[i];

        t_yy_yz_x_x[i] = t_y_yyz_x_x[i] - rab_y[i] * t_y_yz_x_x[i];

        t_yy_yy_z_z[i] = t_y_yyy_z_z[i] - rab_y[i] * t_y_yy_z_z[i];

        t_yy_yy_z_y[i] = t_y_yyy_z_y[i] - rab_y[i] * t_y_yy_z_y[i];

        t_yy_yy_z_x[i] = t_y_yyy_z_x[i] - rab_y[i] * t_y_yy_z_x[i];

        t_yy_yy_y_z[i] = t_y_yyy_y_z[i] - rab_y[i] * t_y_yy_y_z[i];

        t_yy_yy_y_y[i] = t_y_yyy_y_y[i] - rab_y[i] * t_y_yy_y_y[i];

        t_yy_yy_y_x[i] = t_y_yyy_y_x[i] - rab_y[i] * t_y_yy_y_x[i];

        t_yy_yy_x_z[i] = t_y_yyy_x_z[i] - rab_y[i] * t_y_yy_x_z[i];

        t_yy_yy_x_y[i] = t_y_yyy_x_y[i] - rab_y[i] * t_y_yy_x_y[i];

        t_yy_yy_x_x[i] = t_y_yyy_x_x[i] - rab_y[i] * t_y_yy_x_x[i];

        t_yy_xz_z_z[i] = t_y_xyz_z_z[i] - rab_y[i] * t_y_xz_z_z[i];

        t_yy_xz_z_y[i] = t_y_xyz_z_y[i] - rab_y[i] * t_y_xz_z_y[i];

        t_yy_xz_z_x[i] = t_y_xyz_z_x[i] - rab_y[i] * t_y_xz_z_x[i];

        t_yy_xz_y_z[i] = t_y_xyz_y_z[i] - rab_y[i] * t_y_xz_y_z[i];

        t_yy_xz_y_y[i] = t_y_xyz_y_y[i] - rab_y[i] * t_y_xz_y_y[i];

        t_yy_xz_y_x[i] = t_y_xyz_y_x[i] - rab_y[i] * t_y_xz_y_x[i];

        t_yy_xz_x_z[i] = t_y_xyz_x_z[i] - rab_y[i] * t_y_xz_x_z[i];

        t_yy_xz_x_y[i] = t_y_xyz_x_y[i] - rab_y[i] * t_y_xz_x_y[i];

        t_yy_xz_x_x[i] = t_y_xyz_x_x[i] - rab_y[i] * t_y_xz_x_x[i];
    }

    #pragma omp simd align(rab_y, rab_z, t_x_yz_x_x, t_x_yz_x_y, t_x_yz_x_z, t_x_yz_y_x,\
                           t_x_yz_y_y, t_x_yz_y_z, t_x_yz_z_x, t_x_yz_z_y, t_x_yz_z_z,\
                           t_x_yzz_x_x, t_x_yzz_x_y, t_x_yzz_x_z, t_x_yzz_y_x, t_x_yzz_y_y,\
                           t_x_yzz_y_z, t_x_yzz_z_x, t_x_yzz_z_y, t_x_yzz_z_z, t_x_zz_x_x,\
                           t_x_zz_x_y, t_x_zz_x_z, t_x_zz_y_x, t_x_zz_y_y, t_x_zz_y_z,\
                           t_x_zz_z_x, t_x_zz_z_y, t_x_zz_z_z, t_x_zzz_x_x, t_x_zzz_x_y,\
                           t_x_zzz_x_z, t_x_zzz_y_x, t_x_zzz_y_y, t_x_zzz_y_z, t_x_zzz_z_x,\
                           t_x_zzz_z_y, t_x_zzz_z_z, t_xz_yz_x_x, t_xz_yz_x_y, t_xz_yz_x_z,\
                           t_xz_yz_y_x, t_xz_yz_y_y, t_xz_yz_y_z, t_xz_yz_z_x, t_xz_yz_z_y,\
                           t_xz_yz_z_z, t_xz_zz_x_x, t_xz_zz_x_y, t_xz_zz_x_z, t_xz_zz_y_x,\
                           t_xz_zz_y_y, t_xz_zz_y_z, t_xz_zz_z_x, t_xz_zz_z_y, t_xz_zz_z_z,\
                           t_y_xx_x_x, t_y_xx_x_y, t_y_xx_x_z, t_y_xx_y_x, t_y_xx_y_y,\
                           t_y_xx_y_z, t_y_xx_z_x, t_y_xx_z_y, t_y_xx_z_z, t_y_xxy_x_x,\
                           t_y_xxy_x_y, t_y_xxy_x_z, t_y_xxy_y_x, t_y_xxy_y_y, t_y_xxy_y_z,\
                           t_y_xxy_z_x, t_y_xxy_z_y, t_y_xxy_z_z, t_y_xy_x_x, t_y_xy_x_y,\
                           t_y_xy_x_z, t_y_xy_y_x, t_y_xy_y_y, t_y_xy_y_z, t_y_xy_z_x,\
                           t_y_xy_z_y, t_y_xy_z_z, t_y_xyy_x_x, t_y_xyy_x_y, t_y_xyy_x_z,\
                           t_y_xyy_y_x, t_y_xyy_y_y, t_y_xyy_y_z, t_y_xyy_z_x, t_y_xyy_z_y,\
                           t_y_xyy_z_z, t_yy_xx_x_x, t_yy_xx_x_y, t_yy_xx_x_z, t_yy_xx_y_x,\
                           t_yy_xx_y_y, t_yy_xx_y_z, t_yy_xx_z_x, t_yy_xx_z_y, t_yy_xx_z_z,\
                           t_yy_xy_x_x, t_yy_xy_x_y, t_yy_xy_x_z, t_yy_xy_y_x, t_yy_xy_y_y,\
                           t_yy_xy_y_z, t_yy_xy_z_x, t_yy_xy_z_y, t_yy_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_xy_z_z[i] = t_y_xyy_z_z[i] - rab_y[i] * t_y_xy_z_z[i];

        t_yy_xy_z_y[i] = t_y_xyy_z_y[i] - rab_y[i] * t_y_xy_z_y[i];

        t_yy_xy_z_x[i] = t_y_xyy_z_x[i] - rab_y[i] * t_y_xy_z_x[i];

        t_yy_xy_y_z[i] = t_y_xyy_y_z[i] - rab_y[i] * t_y_xy_y_z[i];

        t_yy_xy_y_y[i] = t_y_xyy_y_y[i] - rab_y[i] * t_y_xy_y_y[i];

        t_yy_xy_y_x[i] = t_y_xyy_y_x[i] - rab_y[i] * t_y_xy_y_x[i];

        t_yy_xy_x_z[i] = t_y_xyy_x_z[i] - rab_y[i] * t_y_xy_x_z[i];

        t_yy_xy_x_y[i] = t_y_xyy_x_y[i] - rab_y[i] * t_y_xy_x_y[i];

        t_yy_xy_x_x[i] = t_y_xyy_x_x[i] - rab_y[i] * t_y_xy_x_x[i];

        t_yy_xx_z_z[i] = t_y_xxy_z_z[i] - rab_y[i] * t_y_xx_z_z[i];

        t_yy_xx_z_y[i] = t_y_xxy_z_y[i] - rab_y[i] * t_y_xx_z_y[i];

        t_yy_xx_z_x[i] = t_y_xxy_z_x[i] - rab_y[i] * t_y_xx_z_x[i];

        t_yy_xx_y_z[i] = t_y_xxy_y_z[i] - rab_y[i] * t_y_xx_y_z[i];

        t_yy_xx_y_y[i] = t_y_xxy_y_y[i] - rab_y[i] * t_y_xx_y_y[i];

        t_yy_xx_y_x[i] = t_y_xxy_y_x[i] - rab_y[i] * t_y_xx_y_x[i];

        t_yy_xx_x_z[i] = t_y_xxy_x_z[i] - rab_y[i] * t_y_xx_x_z[i];

        t_yy_xx_x_y[i] = t_y_xxy_x_y[i] - rab_y[i] * t_y_xx_x_y[i];

        t_yy_xx_x_x[i] = t_y_xxy_x_x[i] - rab_y[i] * t_y_xx_x_x[i];

        t_xz_zz_z_z[i] = t_x_zzz_z_z[i] - rab_z[i] * t_x_zz_z_z[i];

        t_xz_zz_z_y[i] = t_x_zzz_z_y[i] - rab_z[i] * t_x_zz_z_y[i];

        t_xz_zz_z_x[i] = t_x_zzz_z_x[i] - rab_z[i] * t_x_zz_z_x[i];

        t_xz_zz_y_z[i] = t_x_zzz_y_z[i] - rab_z[i] * t_x_zz_y_z[i];

        t_xz_zz_y_y[i] = t_x_zzz_y_y[i] - rab_z[i] * t_x_zz_y_y[i];

        t_xz_zz_y_x[i] = t_x_zzz_y_x[i] - rab_z[i] * t_x_zz_y_x[i];

        t_xz_zz_x_z[i] = t_x_zzz_x_z[i] - rab_z[i] * t_x_zz_x_z[i];

        t_xz_zz_x_y[i] = t_x_zzz_x_y[i] - rab_z[i] * t_x_zz_x_y[i];

        t_xz_zz_x_x[i] = t_x_zzz_x_x[i] - rab_z[i] * t_x_zz_x_x[i];

        t_xz_yz_z_z[i] = t_x_yzz_z_z[i] - rab_z[i] * t_x_yz_z_z[i];

        t_xz_yz_z_y[i] = t_x_yzz_z_y[i] - rab_z[i] * t_x_yz_z_y[i];

        t_xz_yz_z_x[i] = t_x_yzz_z_x[i] - rab_z[i] * t_x_yz_z_x[i];

        t_xz_yz_y_z[i] = t_x_yzz_y_z[i] - rab_z[i] * t_x_yz_y_z[i];

        t_xz_yz_y_y[i] = t_x_yzz_y_y[i] - rab_z[i] * t_x_yz_y_y[i];

        t_xz_yz_y_x[i] = t_x_yzz_y_x[i] - rab_z[i] * t_x_yz_y_x[i];

        t_xz_yz_x_z[i] = t_x_yzz_x_z[i] - rab_z[i] * t_x_yz_x_z[i];

        t_xz_yz_x_y[i] = t_x_yzz_x_y[i] - rab_z[i] * t_x_yz_x_y[i];

        t_xz_yz_x_x[i] = t_x_yzz_x_x[i] - rab_z[i] * t_x_yz_x_x[i];
    }

    #pragma omp simd align(rab_z, t_x_xx_x_x, t_x_xx_x_y, t_x_xx_x_z, t_x_xx_y_x, t_x_xx_y_y,\
                           t_x_xx_y_z, t_x_xx_z_x, t_x_xx_z_y, t_x_xx_z_z, t_x_xxz_x_x,\
                           t_x_xxz_x_y, t_x_xxz_x_z, t_x_xxz_y_x, t_x_xxz_y_y, t_x_xxz_y_z,\
                           t_x_xxz_z_x, t_x_xxz_z_y, t_x_xxz_z_z, t_x_xy_x_x, t_x_xy_x_y,\
                           t_x_xy_x_z, t_x_xy_y_x, t_x_xy_y_y, t_x_xy_y_z, t_x_xy_z_x,\
                           t_x_xy_z_y, t_x_xy_z_z, t_x_xyz_x_x, t_x_xyz_x_y, t_x_xyz_x_z,\
                           t_x_xyz_y_x, t_x_xyz_y_y, t_x_xyz_y_z, t_x_xyz_z_x, t_x_xyz_z_y,\
                           t_x_xyz_z_z, t_x_xz_x_x, t_x_xz_x_y, t_x_xz_x_z, t_x_xz_y_x,\
                           t_x_xz_y_y, t_x_xz_y_z, t_x_xz_z_x, t_x_xz_z_y, t_x_xz_z_z,\
                           t_x_xzz_x_x, t_x_xzz_x_y, t_x_xzz_x_z, t_x_xzz_y_x, t_x_xzz_y_y,\
                           t_x_xzz_y_z, t_x_xzz_z_x, t_x_xzz_z_y, t_x_xzz_z_z, t_x_yy_x_x,\
                           t_x_yy_x_y, t_x_yy_x_z, t_x_yy_y_x, t_x_yy_y_y, t_x_yy_y_z,\
                           t_x_yy_z_x, t_x_yy_z_y, t_x_yy_z_z, t_x_yyz_x_x, t_x_yyz_x_y,\
                           t_x_yyz_x_z, t_x_yyz_y_x, t_x_yyz_y_y, t_x_yyz_y_z, t_x_yyz_z_x,\
                           t_x_yyz_z_y, t_x_yyz_z_z, t_xz_xx_x_x, t_xz_xx_x_y, t_xz_xx_x_z,\
                           t_xz_xx_y_x, t_xz_xx_y_y, t_xz_xx_y_z, t_xz_xx_z_x, t_xz_xx_z_y,\
                           t_xz_xx_z_z, t_xz_xy_x_x, t_xz_xy_x_y, t_xz_xy_x_z, t_xz_xy_y_x,\
                           t_xz_xy_y_y, t_xz_xy_y_z, t_xz_xy_z_x, t_xz_xy_z_y, t_xz_xy_z_z,\
                           t_xz_xz_x_x, t_xz_xz_x_y, t_xz_xz_x_z, t_xz_xz_y_x, t_xz_xz_y_y,\
                           t_xz_xz_y_z, t_xz_xz_z_x, t_xz_xz_z_y, t_xz_xz_z_z, t_xz_yy_x_x,\
                           t_xz_yy_x_y, t_xz_yy_x_z, t_xz_yy_y_x, t_xz_yy_y_y, t_xz_yy_y_z,\
                           t_xz_yy_z_x, t_xz_yy_z_y, t_xz_yy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xz_yy_z_z[i] = t_x_yyz_z_z[i] - rab_z[i] * t_x_yy_z_z[i];

        t_xz_yy_z_y[i] = t_x_yyz_z_y[i] - rab_z[i] * t_x_yy_z_y[i];

        t_xz_yy_z_x[i] = t_x_yyz_z_x[i] - rab_z[i] * t_x_yy_z_x[i];

        t_xz_yy_y_z[i] = t_x_yyz_y_z[i] - rab_z[i] * t_x_yy_y_z[i];

        t_xz_yy_y_y[i] = t_x_yyz_y_y[i] - rab_z[i] * t_x_yy_y_y[i];

        t_xz_yy_y_x[i] = t_x_yyz_y_x[i] - rab_z[i] * t_x_yy_y_x[i];

        t_xz_yy_x_z[i] = t_x_yyz_x_z[i] - rab_z[i] * t_x_yy_x_z[i];

        t_xz_yy_x_y[i] = t_x_yyz_x_y[i] - rab_z[i] * t_x_yy_x_y[i];

        t_xz_yy_x_x[i] = t_x_yyz_x_x[i] - rab_z[i] * t_x_yy_x_x[i];

        t_xz_xz_z_z[i] = t_x_xzz_z_z[i] - rab_z[i] * t_x_xz_z_z[i];

        t_xz_xz_z_y[i] = t_x_xzz_z_y[i] - rab_z[i] * t_x_xz_z_y[i];

        t_xz_xz_z_x[i] = t_x_xzz_z_x[i] - rab_z[i] * t_x_xz_z_x[i];

        t_xz_xz_y_z[i] = t_x_xzz_y_z[i] - rab_z[i] * t_x_xz_y_z[i];

        t_xz_xz_y_y[i] = t_x_xzz_y_y[i] - rab_z[i] * t_x_xz_y_y[i];

        t_xz_xz_y_x[i] = t_x_xzz_y_x[i] - rab_z[i] * t_x_xz_y_x[i];

        t_xz_xz_x_z[i] = t_x_xzz_x_z[i] - rab_z[i] * t_x_xz_x_z[i];

        t_xz_xz_x_y[i] = t_x_xzz_x_y[i] - rab_z[i] * t_x_xz_x_y[i];

        t_xz_xz_x_x[i] = t_x_xzz_x_x[i] - rab_z[i] * t_x_xz_x_x[i];

        t_xz_xy_z_z[i] = t_x_xyz_z_z[i] - rab_z[i] * t_x_xy_z_z[i];

        t_xz_xy_z_y[i] = t_x_xyz_z_y[i] - rab_z[i] * t_x_xy_z_y[i];

        t_xz_xy_z_x[i] = t_x_xyz_z_x[i] - rab_z[i] * t_x_xy_z_x[i];

        t_xz_xy_y_z[i] = t_x_xyz_y_z[i] - rab_z[i] * t_x_xy_y_z[i];

        t_xz_xy_y_y[i] = t_x_xyz_y_y[i] - rab_z[i] * t_x_xy_y_y[i];

        t_xz_xy_y_x[i] = t_x_xyz_y_x[i] - rab_z[i] * t_x_xy_y_x[i];

        t_xz_xy_x_z[i] = t_x_xyz_x_z[i] - rab_z[i] * t_x_xy_x_z[i];

        t_xz_xy_x_y[i] = t_x_xyz_x_y[i] - rab_z[i] * t_x_xy_x_y[i];

        t_xz_xy_x_x[i] = t_x_xyz_x_x[i] - rab_z[i] * t_x_xy_x_x[i];

        t_xz_xx_z_z[i] = t_x_xxz_z_z[i] - rab_z[i] * t_x_xx_z_z[i];

        t_xz_xx_z_y[i] = t_x_xxz_z_y[i] - rab_z[i] * t_x_xx_z_y[i];

        t_xz_xx_z_x[i] = t_x_xxz_z_x[i] - rab_z[i] * t_x_xx_z_x[i];

        t_xz_xx_y_z[i] = t_x_xxz_y_z[i] - rab_z[i] * t_x_xx_y_z[i];

        t_xz_xx_y_y[i] = t_x_xxz_y_y[i] - rab_z[i] * t_x_xx_y_y[i];

        t_xz_xx_y_x[i] = t_x_xxz_y_x[i] - rab_z[i] * t_x_xx_y_x[i];

        t_xz_xx_x_z[i] = t_x_xxz_x_z[i] - rab_z[i] * t_x_xx_x_z[i];

        t_xz_xx_x_y[i] = t_x_xxz_x_y[i] - rab_z[i] * t_x_xx_x_y[i];

        t_xz_xx_x_x[i] = t_x_xxz_x_x[i] - rab_z[i] * t_x_xx_x_x[i];
    }

    #pragma omp simd align(rab_y, t_x_xyz_x_x, t_x_xyz_x_y, t_x_xyz_x_z, t_x_xyz_y_x,\
                           t_x_xyz_y_y, t_x_xyz_y_z, t_x_xyz_z_x, t_x_xyz_z_y, t_x_xyz_z_z,\
                           t_x_xz_x_x, t_x_xz_x_y, t_x_xz_x_z, t_x_xz_y_x, t_x_xz_y_y,\
                           t_x_xz_y_z, t_x_xz_z_x, t_x_xz_z_y, t_x_xz_z_z, t_x_yy_x_x,\
                           t_x_yy_x_y, t_x_yy_x_z, t_x_yy_y_x, t_x_yy_y_y, t_x_yy_y_z,\
                           t_x_yy_z_x, t_x_yy_z_y, t_x_yy_z_z, t_x_yyy_x_x, t_x_yyy_x_y,\
                           t_x_yyy_x_z, t_x_yyy_y_x, t_x_yyy_y_y, t_x_yyy_y_z, t_x_yyy_z_x,\
                           t_x_yyy_z_y, t_x_yyy_z_z, t_x_yyz_x_x, t_x_yyz_x_y, t_x_yyz_x_z,\
                           t_x_yyz_y_x, t_x_yyz_y_y, t_x_yyz_y_z, t_x_yyz_z_x, t_x_yyz_z_y,\
                           t_x_yyz_z_z, t_x_yz_x_x, t_x_yz_x_y, t_x_yz_x_z, t_x_yz_y_x,\
                           t_x_yz_y_y, t_x_yz_y_z, t_x_yz_z_x, t_x_yz_z_y, t_x_yz_z_z,\
                           t_x_yzz_x_x, t_x_yzz_x_y, t_x_yzz_x_z, t_x_yzz_y_x, t_x_yzz_y_y,\
                           t_x_yzz_y_z, t_x_yzz_z_x, t_x_yzz_z_y, t_x_yzz_z_z, t_x_zz_x_x,\
                           t_x_zz_x_y, t_x_zz_x_z, t_x_zz_y_x, t_x_zz_y_y, t_x_zz_y_z,\
                           t_x_zz_z_x, t_x_zz_z_y, t_x_zz_z_z, t_xy_xz_x_x, t_xy_xz_x_y,\
                           t_xy_xz_x_z, t_xy_xz_y_x, t_xy_xz_y_y, t_xy_xz_y_z, t_xy_xz_z_x,\
                           t_xy_xz_z_y, t_xy_xz_z_z, t_xy_yy_x_x, t_xy_yy_x_y, t_xy_yy_x_z,\
                           t_xy_yy_y_x, t_xy_yy_y_y, t_xy_yy_y_z, t_xy_yy_z_x, t_xy_yy_z_y,\
                           t_xy_yy_z_z, t_xy_yz_x_x, t_xy_yz_x_y, t_xy_yz_x_z, t_xy_yz_y_x,\
                           t_xy_yz_y_y, t_xy_yz_y_z, t_xy_yz_z_x, t_xy_yz_z_y, t_xy_yz_z_z,\
                           t_xy_zz_x_x, t_xy_zz_x_y, t_xy_zz_x_z, t_xy_zz_y_x, t_xy_zz_y_y,\
                           t_xy_zz_y_z, t_xy_zz_z_x, t_xy_zz_z_y, t_xy_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_zz_z_z[i] = t_x_yzz_z_z[i] - rab_y[i] * t_x_zz_z_z[i];

        t_xy_zz_z_y[i] = t_x_yzz_z_y[i] - rab_y[i] * t_x_zz_z_y[i];

        t_xy_zz_z_x[i] = t_x_yzz_z_x[i] - rab_y[i] * t_x_zz_z_x[i];

        t_xy_zz_y_z[i] = t_x_yzz_y_z[i] - rab_y[i] * t_x_zz_y_z[i];

        t_xy_zz_y_y[i] = t_x_yzz_y_y[i] - rab_y[i] * t_x_zz_y_y[i];

        t_xy_zz_y_x[i] = t_x_yzz_y_x[i] - rab_y[i] * t_x_zz_y_x[i];

        t_xy_zz_x_z[i] = t_x_yzz_x_z[i] - rab_y[i] * t_x_zz_x_z[i];

        t_xy_zz_x_y[i] = t_x_yzz_x_y[i] - rab_y[i] * t_x_zz_x_y[i];

        t_xy_zz_x_x[i] = t_x_yzz_x_x[i] - rab_y[i] * t_x_zz_x_x[i];

        t_xy_yz_z_z[i] = t_x_yyz_z_z[i] - rab_y[i] * t_x_yz_z_z[i];

        t_xy_yz_z_y[i] = t_x_yyz_z_y[i] - rab_y[i] * t_x_yz_z_y[i];

        t_xy_yz_z_x[i] = t_x_yyz_z_x[i] - rab_y[i] * t_x_yz_z_x[i];

        t_xy_yz_y_z[i] = t_x_yyz_y_z[i] - rab_y[i] * t_x_yz_y_z[i];

        t_xy_yz_y_y[i] = t_x_yyz_y_y[i] - rab_y[i] * t_x_yz_y_y[i];

        t_xy_yz_y_x[i] = t_x_yyz_y_x[i] - rab_y[i] * t_x_yz_y_x[i];

        t_xy_yz_x_z[i] = t_x_yyz_x_z[i] - rab_y[i] * t_x_yz_x_z[i];

        t_xy_yz_x_y[i] = t_x_yyz_x_y[i] - rab_y[i] * t_x_yz_x_y[i];

        t_xy_yz_x_x[i] = t_x_yyz_x_x[i] - rab_y[i] * t_x_yz_x_x[i];

        t_xy_yy_z_z[i] = t_x_yyy_z_z[i] - rab_y[i] * t_x_yy_z_z[i];

        t_xy_yy_z_y[i] = t_x_yyy_z_y[i] - rab_y[i] * t_x_yy_z_y[i];

        t_xy_yy_z_x[i] = t_x_yyy_z_x[i] - rab_y[i] * t_x_yy_z_x[i];

        t_xy_yy_y_z[i] = t_x_yyy_y_z[i] - rab_y[i] * t_x_yy_y_z[i];

        t_xy_yy_y_y[i] = t_x_yyy_y_y[i] - rab_y[i] * t_x_yy_y_y[i];

        t_xy_yy_y_x[i] = t_x_yyy_y_x[i] - rab_y[i] * t_x_yy_y_x[i];

        t_xy_yy_x_z[i] = t_x_yyy_x_z[i] - rab_y[i] * t_x_yy_x_z[i];

        t_xy_yy_x_y[i] = t_x_yyy_x_y[i] - rab_y[i] * t_x_yy_x_y[i];

        t_xy_yy_x_x[i] = t_x_yyy_x_x[i] - rab_y[i] * t_x_yy_x_x[i];

        t_xy_xz_z_z[i] = t_x_xyz_z_z[i] - rab_y[i] * t_x_xz_z_z[i];

        t_xy_xz_z_y[i] = t_x_xyz_z_y[i] - rab_y[i] * t_x_xz_z_y[i];

        t_xy_xz_z_x[i] = t_x_xyz_z_x[i] - rab_y[i] * t_x_xz_z_x[i];

        t_xy_xz_y_z[i] = t_x_xyz_y_z[i] - rab_y[i] * t_x_xz_y_z[i];

        t_xy_xz_y_y[i] = t_x_xyz_y_y[i] - rab_y[i] * t_x_xz_y_y[i];

        t_xy_xz_y_x[i] = t_x_xyz_y_x[i] - rab_y[i] * t_x_xz_y_x[i];

        t_xy_xz_x_z[i] = t_x_xyz_x_z[i] - rab_y[i] * t_x_xz_x_z[i];

        t_xy_xz_x_y[i] = t_x_xyz_x_y[i] - rab_y[i] * t_x_xz_x_y[i];

        t_xy_xz_x_x[i] = t_x_xyz_x_x[i] - rab_y[i] * t_x_xz_x_x[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_x_xx_x_x, t_x_xx_x_y, t_x_xx_x_z, t_x_xx_y_x,\
                           t_x_xx_y_y, t_x_xx_y_z, t_x_xx_z_x, t_x_xx_z_y, t_x_xx_z_z,\
                           t_x_xxy_x_x, t_x_xxy_x_y, t_x_xxy_x_z, t_x_xxy_y_x, t_x_xxy_y_y,\
                           t_x_xxy_y_z, t_x_xxy_z_x, t_x_xxy_z_y, t_x_xxy_z_z, t_x_xy_x_x,\
                           t_x_xy_x_y, t_x_xy_x_z, t_x_xy_y_x, t_x_xy_y_y, t_x_xy_y_z,\
                           t_x_xy_z_x, t_x_xy_z_y, t_x_xy_z_z, t_x_xyy_x_x, t_x_xyy_x_y,\
                           t_x_xyy_x_z, t_x_xyy_y_x, t_x_xyy_y_y, t_x_xyy_y_z, t_x_xyy_z_x,\
                           t_x_xyy_z_y, t_x_xyy_z_z, t_x_xyz_x_x, t_x_xyz_x_y, t_x_xyz_x_z,\
                           t_x_xyz_y_x, t_x_xyz_y_y, t_x_xyz_y_z, t_x_xyz_z_x, t_x_xyz_z_y,\
                           t_x_xyz_z_z, t_x_xzz_x_x, t_x_xzz_x_y, t_x_xzz_x_z, t_x_xzz_y_x,\
                           t_x_xzz_y_y, t_x_xzz_y_z, t_x_xzz_z_x, t_x_xzz_z_y, t_x_xzz_z_z,\
                           t_x_yz_x_x, t_x_yz_x_y, t_x_yz_x_z, t_x_yz_y_x, t_x_yz_y_y,\
                           t_x_yz_y_z, t_x_yz_z_x, t_x_yz_z_y, t_x_yz_z_z, t_x_zz_x_x,\
                           t_x_zz_x_y, t_x_zz_x_z, t_x_zz_y_x, t_x_zz_y_y, t_x_zz_y_z,\
                           t_x_zz_z_x, t_x_zz_z_y, t_x_zz_z_z, t_xx_yz_x_x, t_xx_yz_x_y,\
                           t_xx_yz_x_z, t_xx_yz_y_x, t_xx_yz_y_y, t_xx_yz_y_z, t_xx_yz_z_x,\
                           t_xx_yz_z_y, t_xx_yz_z_z, t_xx_zz_x_x, t_xx_zz_x_y, t_xx_zz_x_z,\
                           t_xx_zz_y_x, t_xx_zz_y_y, t_xx_zz_y_z, t_xx_zz_z_x, t_xx_zz_z_y,\
                           t_xx_zz_z_z, t_xy_xx_x_x, t_xy_xx_x_y, t_xy_xx_x_z, t_xy_xx_y_x,\
                           t_xy_xx_y_y, t_xy_xx_y_z, t_xy_xx_z_x, t_xy_xx_z_y, t_xy_xx_z_z,\
                           t_xy_xy_x_x, t_xy_xy_x_y, t_xy_xy_x_z, t_xy_xy_y_x, t_xy_xy_y_y,\
                           t_xy_xy_y_z, t_xy_xy_z_x, t_xy_xy_z_y, t_xy_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_xy_z_z[i] = t_x_xyy_z_z[i] - rab_y[i] * t_x_xy_z_z[i];

        t_xy_xy_z_y[i] = t_x_xyy_z_y[i] - rab_y[i] * t_x_xy_z_y[i];

        t_xy_xy_z_x[i] = t_x_xyy_z_x[i] - rab_y[i] * t_x_xy_z_x[i];

        t_xy_xy_y_z[i] = t_x_xyy_y_z[i] - rab_y[i] * t_x_xy_y_z[i];

        t_xy_xy_y_y[i] = t_x_xyy_y_y[i] - rab_y[i] * t_x_xy_y_y[i];

        t_xy_xy_y_x[i] = t_x_xyy_y_x[i] - rab_y[i] * t_x_xy_y_x[i];

        t_xy_xy_x_z[i] = t_x_xyy_x_z[i] - rab_y[i] * t_x_xy_x_z[i];

        t_xy_xy_x_y[i] = t_x_xyy_x_y[i] - rab_y[i] * t_x_xy_x_y[i];

        t_xy_xy_x_x[i] = t_x_xyy_x_x[i] - rab_y[i] * t_x_xy_x_x[i];

        t_xy_xx_z_z[i] = t_x_xxy_z_z[i] - rab_y[i] * t_x_xx_z_z[i];

        t_xy_xx_z_y[i] = t_x_xxy_z_y[i] - rab_y[i] * t_x_xx_z_y[i];

        t_xy_xx_z_x[i] = t_x_xxy_z_x[i] - rab_y[i] * t_x_xx_z_x[i];

        t_xy_xx_y_z[i] = t_x_xxy_y_z[i] - rab_y[i] * t_x_xx_y_z[i];

        t_xy_xx_y_y[i] = t_x_xxy_y_y[i] - rab_y[i] * t_x_xx_y_y[i];

        t_xy_xx_y_x[i] = t_x_xxy_y_x[i] - rab_y[i] * t_x_xx_y_x[i];

        t_xy_xx_x_z[i] = t_x_xxy_x_z[i] - rab_y[i] * t_x_xx_x_z[i];

        t_xy_xx_x_y[i] = t_x_xxy_x_y[i] - rab_y[i] * t_x_xx_x_y[i];

        t_xy_xx_x_x[i] = t_x_xxy_x_x[i] - rab_y[i] * t_x_xx_x_x[i];

        t_xx_zz_z_z[i] = t_x_xzz_z_z[i] - rab_x[i] * t_x_zz_z_z[i];

        t_xx_zz_z_y[i] = t_x_xzz_z_y[i] - rab_x[i] * t_x_zz_z_y[i];

        t_xx_zz_z_x[i] = t_x_xzz_z_x[i] - rab_x[i] * t_x_zz_z_x[i];

        t_xx_zz_y_z[i] = t_x_xzz_y_z[i] - rab_x[i] * t_x_zz_y_z[i];

        t_xx_zz_y_y[i] = t_x_xzz_y_y[i] - rab_x[i] * t_x_zz_y_y[i];

        t_xx_zz_y_x[i] = t_x_xzz_y_x[i] - rab_x[i] * t_x_zz_y_x[i];

        t_xx_zz_x_z[i] = t_x_xzz_x_z[i] - rab_x[i] * t_x_zz_x_z[i];

        t_xx_zz_x_y[i] = t_x_xzz_x_y[i] - rab_x[i] * t_x_zz_x_y[i];

        t_xx_zz_x_x[i] = t_x_xzz_x_x[i] - rab_x[i] * t_x_zz_x_x[i];

        t_xx_yz_z_z[i] = t_x_xyz_z_z[i] - rab_x[i] * t_x_yz_z_z[i];

        t_xx_yz_z_y[i] = t_x_xyz_z_y[i] - rab_x[i] * t_x_yz_z_y[i];

        t_xx_yz_z_x[i] = t_x_xyz_z_x[i] - rab_x[i] * t_x_yz_z_x[i];

        t_xx_yz_y_z[i] = t_x_xyz_y_z[i] - rab_x[i] * t_x_yz_y_z[i];

        t_xx_yz_y_y[i] = t_x_xyz_y_y[i] - rab_x[i] * t_x_yz_y_y[i];

        t_xx_yz_y_x[i] = t_x_xyz_y_x[i] - rab_x[i] * t_x_yz_y_x[i];

        t_xx_yz_x_z[i] = t_x_xyz_x_z[i] - rab_x[i] * t_x_yz_x_z[i];

        t_xx_yz_x_y[i] = t_x_xyz_x_y[i] - rab_x[i] * t_x_yz_x_y[i];

        t_xx_yz_x_x[i] = t_x_xyz_x_x[i] - rab_x[i] * t_x_yz_x_x[i];
    }

    #pragma omp simd align(rab_x, t_x_xx_x_x, t_x_xx_x_y, t_x_xx_x_z, t_x_xx_y_x, t_x_xx_y_y,\
                           t_x_xx_y_z, t_x_xx_z_x, t_x_xx_z_y, t_x_xx_z_z, t_x_xxx_x_x,\
                           t_x_xxx_x_y, t_x_xxx_x_z, t_x_xxx_y_x, t_x_xxx_y_y, t_x_xxx_y_z,\
                           t_x_xxx_z_x, t_x_xxx_z_y, t_x_xxx_z_z, t_x_xxy_x_x, t_x_xxy_x_y,\
                           t_x_xxy_x_z, t_x_xxy_y_x, t_x_xxy_y_y, t_x_xxy_y_z, t_x_xxy_z_x,\
                           t_x_xxy_z_y, t_x_xxy_z_z, t_x_xxz_x_x, t_x_xxz_x_y, t_x_xxz_x_z,\
                           t_x_xxz_y_x, t_x_xxz_y_y, t_x_xxz_y_z, t_x_xxz_z_x, t_x_xxz_z_y,\
                           t_x_xxz_z_z, t_x_xy_x_x, t_x_xy_x_y, t_x_xy_x_z, t_x_xy_y_x,\
                           t_x_xy_y_y, t_x_xy_y_z, t_x_xy_z_x, t_x_xy_z_y, t_x_xy_z_z,\
                           t_x_xyy_x_x, t_x_xyy_x_y, t_x_xyy_x_z, t_x_xyy_y_x, t_x_xyy_y_y,\
                           t_x_xyy_y_z, t_x_xyy_z_x, t_x_xyy_z_y, t_x_xyy_z_z, t_x_xz_x_x,\
                           t_x_xz_x_y, t_x_xz_x_z, t_x_xz_y_x, t_x_xz_y_y, t_x_xz_y_z,\
                           t_x_xz_z_x, t_x_xz_z_y, t_x_xz_z_z, t_x_yy_x_x, t_x_yy_x_y,\
                           t_x_yy_x_z, t_x_yy_y_x, t_x_yy_y_y, t_x_yy_y_z, t_x_yy_z_x,\
                           t_x_yy_z_y, t_x_yy_z_z, t_xx_xx_x_x, t_xx_xx_x_y, t_xx_xx_x_z,\
                           t_xx_xx_y_x, t_xx_xx_y_y, t_xx_xx_y_z, t_xx_xx_z_x, t_xx_xx_z_y,\
                           t_xx_xx_z_z, t_xx_xy_x_x, t_xx_xy_x_y, t_xx_xy_x_z, t_xx_xy_y_x,\
                           t_xx_xy_y_y, t_xx_xy_y_z, t_xx_xy_z_x, t_xx_xy_z_y, t_xx_xy_z_z,\
                           t_xx_xz_x_x, t_xx_xz_x_y, t_xx_xz_x_z, t_xx_xz_y_x, t_xx_xz_y_y,\
                           t_xx_xz_y_z, t_xx_xz_z_x, t_xx_xz_z_y, t_xx_xz_z_z, t_xx_yy_x_x,\
                           t_xx_yy_x_y, t_xx_yy_x_z, t_xx_yy_y_x, t_xx_yy_y_y, t_xx_yy_y_z,\
                           t_xx_yy_z_x, t_xx_yy_z_y, t_xx_yy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xx_yy_z_z[i] = t_x_xyy_z_z[i] - rab_x[i] * t_x_yy_z_z[i];

        t_xx_yy_z_y[i] = t_x_xyy_z_y[i] - rab_x[i] * t_x_yy_z_y[i];

        t_xx_yy_z_x[i] = t_x_xyy_z_x[i] - rab_x[i] * t_x_yy_z_x[i];

        t_xx_yy_y_z[i] = t_x_xyy_y_z[i] - rab_x[i] * t_x_yy_y_z[i];

        t_xx_yy_y_y[i] = t_x_xyy_y_y[i] - rab_x[i] * t_x_yy_y_y[i];

        t_xx_yy_y_x[i] = t_x_xyy_y_x[i] - rab_x[i] * t_x_yy_y_x[i];

        t_xx_yy_x_z[i] = t_x_xyy_x_z[i] - rab_x[i] * t_x_yy_x_z[i];

        t_xx_yy_x_y[i] = t_x_xyy_x_y[i] - rab_x[i] * t_x_yy_x_y[i];

        t_xx_yy_x_x[i] = t_x_xyy_x_x[i] - rab_x[i] * t_x_yy_x_x[i];

        t_xx_xz_z_z[i] = t_x_xxz_z_z[i] - rab_x[i] * t_x_xz_z_z[i];

        t_xx_xz_z_y[i] = t_x_xxz_z_y[i] - rab_x[i] * t_x_xz_z_y[i];

        t_xx_xz_z_x[i] = t_x_xxz_z_x[i] - rab_x[i] * t_x_xz_z_x[i];

        t_xx_xz_y_z[i] = t_x_xxz_y_z[i] - rab_x[i] * t_x_xz_y_z[i];

        t_xx_xz_y_y[i] = t_x_xxz_y_y[i] - rab_x[i] * t_x_xz_y_y[i];

        t_xx_xz_y_x[i] = t_x_xxz_y_x[i] - rab_x[i] * t_x_xz_y_x[i];

        t_xx_xz_x_z[i] = t_x_xxz_x_z[i] - rab_x[i] * t_x_xz_x_z[i];

        t_xx_xz_x_y[i] = t_x_xxz_x_y[i] - rab_x[i] * t_x_xz_x_y[i];

        t_xx_xz_x_x[i] = t_x_xxz_x_x[i] - rab_x[i] * t_x_xz_x_x[i];

        t_xx_xy_z_z[i] = t_x_xxy_z_z[i] - rab_x[i] * t_x_xy_z_z[i];

        t_xx_xy_z_y[i] = t_x_xxy_z_y[i] - rab_x[i] * t_x_xy_z_y[i];

        t_xx_xy_z_x[i] = t_x_xxy_z_x[i] - rab_x[i] * t_x_xy_z_x[i];

        t_xx_xy_y_z[i] = t_x_xxy_y_z[i] - rab_x[i] * t_x_xy_y_z[i];

        t_xx_xy_y_y[i] = t_x_xxy_y_y[i] - rab_x[i] * t_x_xy_y_y[i];

        t_xx_xy_y_x[i] = t_x_xxy_y_x[i] - rab_x[i] * t_x_xy_y_x[i];

        t_xx_xy_x_z[i] = t_x_xxy_x_z[i] - rab_x[i] * t_x_xy_x_z[i];

        t_xx_xy_x_y[i] = t_x_xxy_x_y[i] - rab_x[i] * t_x_xy_x_y[i];

        t_xx_xy_x_x[i] = t_x_xxy_x_x[i] - rab_x[i] * t_x_xy_x_x[i];

        t_xx_xx_z_z[i] = t_x_xxx_z_z[i] - rab_x[i] * t_x_xx_z_z[i];

        t_xx_xx_z_y[i] = t_x_xxx_z_y[i] - rab_x[i] * t_x_xx_z_y[i];

        t_xx_xx_z_x[i] = t_x_xxx_z_x[i] - rab_x[i] * t_x_xx_z_x[i];

        t_xx_xx_y_z[i] = t_x_xxx_y_z[i] - rab_x[i] * t_x_xx_y_z[i];

        t_xx_xx_y_y[i] = t_x_xxx_y_y[i] - rab_x[i] * t_x_xx_y_y[i];

        t_xx_xx_y_x[i] = t_x_xxx_y_x[i] - rab_x[i] * t_x_xx_y_x[i];

        t_xx_xx_x_z[i] = t_x_xxx_x_z[i] - rab_x[i] * t_x_xx_x_z[i];

        t_xx_xx_x_y[i] = t_x_xxx_x_y[i] - rab_x[i] * t_x_xx_x_y[i];

        t_xx_xx_x_x[i] = t_x_xxx_x_x[i] - rab_x[i] * t_x_xx_x_x[i];
    }
}


} // derirec namespace
