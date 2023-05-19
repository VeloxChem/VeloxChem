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
compHostHRRForSGPD_V0(      BufferHostXY<T>&      intsBufferSGPD,
                      const BufferHostX<int32_t>& intsIndexesSGPD,
                      const BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD,
                      const BufferHostXY<T>&      intsBufferSGSF,
                      const BufferHostX<int32_t>& intsIndexesSGSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SGSD) integral components

    t_0_zzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(0));

    t_0_zzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(1));

    t_0_zzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(2));

    t_0_zzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(3));

    t_0_zzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(4));

    t_0_zzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(5));

    t_0_yzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(6));

    t_0_yzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(7));

    t_0_yzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(8));

    t_0_yzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(9));

    t_0_yzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(10));

    t_0_yzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(11));

    t_0_yyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(12));

    t_0_yyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(13));

    t_0_yyzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(14));

    t_0_yyzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(15));

    t_0_yyzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(16));

    t_0_yyzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(17));

    t_0_yyyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(18));

    t_0_yyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(19));

    t_0_yyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(20));

    t_0_yyyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(21));

    t_0_yyyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(22));

    t_0_yyyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(23));

    t_0_yyyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(24));

    t_0_yyyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(25));

    t_0_yyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(26));

    t_0_yyyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(27));

    t_0_yyyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(28));

    t_0_yyyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(29));

    t_0_xzzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(30));

    t_0_xzzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(31));

    t_0_xzzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(32));

    t_0_xzzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(33));

    t_0_xzzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(34));

    t_0_xzzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(35));

    t_0_xyzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(36));

    t_0_xyzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(37));

    t_0_xyzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(38));

    t_0_xyzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(39));

    t_0_xyzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(40));

    t_0_xyzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(41));

    t_0_xyyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(42));

    t_0_xyyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(43));

    t_0_xyyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(44));

    t_0_xyyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(45));

    t_0_xyyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(46));

    t_0_xyyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(47));

    t_0_xyyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(48));

    t_0_xyyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(49));

    t_0_xyyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(50));

    t_0_xyyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(51));

    t_0_xyyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(52));

    t_0_xyyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(53));

    t_0_xxzz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(54));

    t_0_xxzz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(55));

    t_0_xxzz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(56));

    t_0_xxzz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(57));

    t_0_xxzz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(58));

    t_0_xxzz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(59));

    t_0_xxyz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(60));

    t_0_xxyz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(61));

    t_0_xxyz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(62));

    t_0_xxyz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(63));

    t_0_xxyz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(64));

    t_0_xxyz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(65));

    t_0_xxyy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(66));

    t_0_xxyy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(67));

    t_0_xxyy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(68));

    t_0_xxyy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(69));

    t_0_xxyy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(70));

    t_0_xxyy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(71));

    t_0_xxxz_0_zz = intsBufferSGSD.data(intsIndexesSGSD(72));

    t_0_xxxz_0_yz = intsBufferSGSD.data(intsIndexesSGSD(73));

    t_0_xxxz_0_yy = intsBufferSGSD.data(intsIndexesSGSD(74));

    t_0_xxxz_0_xz = intsBufferSGSD.data(intsIndexesSGSD(75));

    t_0_xxxz_0_xy = intsBufferSGSD.data(intsIndexesSGSD(76));

    t_0_xxxz_0_xx = intsBufferSGSD.data(intsIndexesSGSD(77));

    t_0_xxxy_0_zz = intsBufferSGSD.data(intsIndexesSGSD(78));

    t_0_xxxy_0_yz = intsBufferSGSD.data(intsIndexesSGSD(79));

    t_0_xxxy_0_yy = intsBufferSGSD.data(intsIndexesSGSD(80));

    t_0_xxxy_0_xz = intsBufferSGSD.data(intsIndexesSGSD(81));

    t_0_xxxy_0_xy = intsBufferSGSD.data(intsIndexesSGSD(82));

    t_0_xxxy_0_xx = intsBufferSGSD.data(intsIndexesSGSD(83));

    t_0_xxxx_0_zz = intsBufferSGSD.data(intsIndexesSGSD(84));

    t_0_xxxx_0_yz = intsBufferSGSD.data(intsIndexesSGSD(85));

    t_0_xxxx_0_yy = intsBufferSGSD.data(intsIndexesSGSD(86));

    t_0_xxxx_0_xz = intsBufferSGSD.data(intsIndexesSGSD(87));

    t_0_xxxx_0_xy = intsBufferSGSD.data(intsIndexesSGSD(88));

    t_0_xxxx_0_xx = intsBufferSGSD.data(intsIndexesSGSD(89));

    // set up (SGSF) integral components

    t_0_zzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(0));

    t_0_zzzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(1));

    t_0_zzzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(2));

    t_0_zzzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(3));

    t_0_zzzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(4));

    t_0_zzzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(5));

    t_0_zzzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(6));

    t_0_zzzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(7));

    t_0_zzzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(8));

    t_0_zzzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(9));

    t_0_yzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(10));

    t_0_yzzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(11));

    t_0_yzzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(12));

    t_0_yzzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(13));

    t_0_yzzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(14));

    t_0_yzzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(15));

    t_0_yzzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(16));

    t_0_yzzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(17));

    t_0_yzzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(18));

    t_0_yzzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(19));

    t_0_yyzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(20));

    t_0_yyzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(21));

    t_0_yyzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(22));

    t_0_yyzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(23));

    t_0_yyzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(24));

    t_0_yyzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(25));

    t_0_yyzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(26));

    t_0_yyzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(27));

    t_0_yyzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(28));

    t_0_yyzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(29));

    t_0_yyyz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(30));

    t_0_yyyz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(31));

    t_0_yyyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(32));

    t_0_yyyz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(33));

    t_0_yyyz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(34));

    t_0_yyyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(35));

    t_0_yyyz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(36));

    t_0_yyyz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(37));

    t_0_yyyz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(38));

    t_0_yyyz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(39));

    t_0_yyyy_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(40));

    t_0_yyyy_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(41));

    t_0_yyyy_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(42));

    t_0_yyyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(43));

    t_0_yyyy_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(44));

    t_0_yyyy_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(45));

    t_0_yyyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(46));

    t_0_yyyy_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(47));

    t_0_yyyy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(48));

    t_0_yyyy_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(49));

    t_0_xzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(50));

    t_0_xzzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(51));

    t_0_xzzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(52));

    t_0_xzzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(53));

    t_0_xzzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(54));

    t_0_xzzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(55));

    t_0_xzzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(56));

    t_0_xzzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(57));

    t_0_xzzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(58));

    t_0_xzzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(59));

    t_0_xyzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(60));

    t_0_xyzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(61));

    t_0_xyzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(62));

    t_0_xyzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(63));

    t_0_xyzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(64));

    t_0_xyzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(65));

    t_0_xyzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(66));

    t_0_xyzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(67));

    t_0_xyzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(68));

    t_0_xyzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(69));

    t_0_xyyz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(70));

    t_0_xyyz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(71));

    t_0_xyyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(72));

    t_0_xyyz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(73));

    t_0_xyyz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(74));

    t_0_xyyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(75));

    t_0_xyyz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(76));

    t_0_xyyz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(77));

    t_0_xyyz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(78));

    t_0_xyyz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(79));

    t_0_xyyy_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(80));

    t_0_xyyy_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(81));

    t_0_xyyy_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(82));

    t_0_xyyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(83));

    t_0_xyyy_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(84));

    t_0_xyyy_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(85));

    t_0_xyyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(86));

    t_0_xyyy_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(87));

    t_0_xyyy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(88));

    t_0_xyyy_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(89));

    t_0_xxzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(90));

    t_0_xxzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(91));

    t_0_xxzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(92));

    t_0_xxzz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(93));

    t_0_xxzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(94));

    t_0_xxzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(95));

    t_0_xxzz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(96));

    t_0_xxzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(97));

    t_0_xxzz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(98));

    t_0_xxzz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(99));

    t_0_xxyz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(100));

    t_0_xxyz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(101));

    t_0_xxyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(102));

    t_0_xxyz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(103));

    t_0_xxyz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(104));

    t_0_xxyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(105));

    t_0_xxyz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(106));

    t_0_xxyz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(107));

    t_0_xxyz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(108));

    t_0_xxyz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(109));

    t_0_xxyy_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(110));

    t_0_xxyy_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(111));

    t_0_xxyy_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(112));

    t_0_xxyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(113));

    t_0_xxyy_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(114));

    t_0_xxyy_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(115));

    t_0_xxyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(116));

    t_0_xxyy_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(117));

    t_0_xxyy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(118));

    t_0_xxyy_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(119));

    t_0_xxxz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(120));

    t_0_xxxz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(121));

    t_0_xxxz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(122));

    t_0_xxxz_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(123));

    t_0_xxxz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(124));

    t_0_xxxz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(125));

    t_0_xxxz_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(126));

    t_0_xxxz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(127));

    t_0_xxxz_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(128));

    t_0_xxxz_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(129));

    t_0_xxxy_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(130));

    t_0_xxxy_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(131));

    t_0_xxxy_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(132));

    t_0_xxxy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(133));

    t_0_xxxy_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(134));

    t_0_xxxy_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(135));

    t_0_xxxy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(136));

    t_0_xxxy_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(137));

    t_0_xxxy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(138));

    t_0_xxxy_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(139));

    t_0_xxxx_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(140));

    t_0_xxxx_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(141));

    t_0_xxxx_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(142));

    t_0_xxxx_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(143));

    t_0_xxxx_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(144));

    t_0_xxxx_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(145));

    t_0_xxxx_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(146));

    t_0_xxxx_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(147));

    t_0_xxxx_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(148));

    t_0_xxxx_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(149));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzzz_0_xx, t_0_yzzz_0_xxx, t_0_yzzz_0_xxy,\
                           t_0_yzzz_0_xxz, t_0_yzzz_0_xy, t_0_yzzz_0_xyy, t_0_yzzz_0_xyz,\
                           t_0_yzzz_0_xz, t_0_yzzz_0_xzz, t_0_yzzz_0_yy, t_0_yzzz_0_yyy,\
                           t_0_yzzz_0_yyz, t_0_yzzz_0_yz, t_0_yzzz_0_yzz, t_0_yzzz_0_zz,\
                           t_0_yzzz_0_zzz, t_0_yzzz_x_xx, t_0_yzzz_x_xy, t_0_yzzz_x_xz,\
                           t_0_yzzz_x_yy, t_0_yzzz_x_yz, t_0_yzzz_x_zz, t_0_yzzz_y_xx,\
                           t_0_yzzz_y_xy, t_0_yzzz_y_xz, t_0_yzzz_y_yy, t_0_yzzz_y_yz,\
                           t_0_yzzz_y_zz, t_0_yzzz_z_xx, t_0_yzzz_z_xy, t_0_yzzz_z_xz,\
                           t_0_yzzz_z_yy, t_0_yzzz_z_yz, t_0_yzzz_z_zz, t_0_zzzz_0_xx,\
                           t_0_zzzz_0_xxx, t_0_zzzz_0_xxy, t_0_zzzz_0_xxz, t_0_zzzz_0_xy,\
                           t_0_zzzz_0_xyy, t_0_zzzz_0_xyz, t_0_zzzz_0_xz, t_0_zzzz_0_xzz,\
                           t_0_zzzz_0_yy, t_0_zzzz_0_yyy, t_0_zzzz_0_yyz, t_0_zzzz_0_yz,\
                           t_0_zzzz_0_yzz, t_0_zzzz_0_zz, t_0_zzzz_0_zzz, t_0_zzzz_x_xx,\
                           t_0_zzzz_x_xy, t_0_zzzz_x_xz, t_0_zzzz_x_yy, t_0_zzzz_x_yz,\
                           t_0_zzzz_x_zz, t_0_zzzz_y_xx, t_0_zzzz_y_xy, t_0_zzzz_y_xz,\
                           t_0_zzzz_y_yy, t_0_zzzz_y_yz, t_0_zzzz_y_zz, t_0_zzzz_z_xx,\
                           t_0_zzzz_z_xy, t_0_zzzz_z_xz, t_0_zzzz_z_yy, t_0_zzzz_z_yz,\
                           t_0_zzzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_z_zz[i] = t_0_zzzz_0_zzz[i] - rcd_z[i] * t_0_zzzz_0_zz[i];

        t_0_zzzz_z_yz[i] = t_0_zzzz_0_yzz[i] - rcd_z[i] * t_0_zzzz_0_yz[i];

        t_0_zzzz_z_yy[i] = t_0_zzzz_0_yyz[i] - rcd_z[i] * t_0_zzzz_0_yy[i];

        t_0_zzzz_z_xz[i] = t_0_zzzz_0_xzz[i] - rcd_z[i] * t_0_zzzz_0_xz[i];

        t_0_zzzz_z_xy[i] = t_0_zzzz_0_xyz[i] - rcd_z[i] * t_0_zzzz_0_xy[i];

        t_0_zzzz_z_xx[i] = t_0_zzzz_0_xxz[i] - rcd_z[i] * t_0_zzzz_0_xx[i];

        t_0_zzzz_y_zz[i] = t_0_zzzz_0_yzz[i] - rcd_y[i] * t_0_zzzz_0_zz[i];

        t_0_zzzz_y_yz[i] = t_0_zzzz_0_yyz[i] - rcd_y[i] * t_0_zzzz_0_yz[i];

        t_0_zzzz_y_yy[i] = t_0_zzzz_0_yyy[i] - rcd_y[i] * t_0_zzzz_0_yy[i];

        t_0_zzzz_y_xz[i] = t_0_zzzz_0_xyz[i] - rcd_y[i] * t_0_zzzz_0_xz[i];

        t_0_zzzz_y_xy[i] = t_0_zzzz_0_xyy[i] - rcd_y[i] * t_0_zzzz_0_xy[i];

        t_0_zzzz_y_xx[i] = t_0_zzzz_0_xxy[i] - rcd_y[i] * t_0_zzzz_0_xx[i];

        t_0_zzzz_x_zz[i] = t_0_zzzz_0_xzz[i] - rcd_x[i] * t_0_zzzz_0_zz[i];

        t_0_zzzz_x_yz[i] = t_0_zzzz_0_xyz[i] - rcd_x[i] * t_0_zzzz_0_yz[i];

        t_0_zzzz_x_yy[i] = t_0_zzzz_0_xyy[i] - rcd_x[i] * t_0_zzzz_0_yy[i];

        t_0_zzzz_x_xz[i] = t_0_zzzz_0_xxz[i] - rcd_x[i] * t_0_zzzz_0_xz[i];

        t_0_zzzz_x_xy[i] = t_0_zzzz_0_xxy[i] - rcd_x[i] * t_0_zzzz_0_xy[i];

        t_0_zzzz_x_xx[i] = t_0_zzzz_0_xxx[i] - rcd_x[i] * t_0_zzzz_0_xx[i];

        t_0_yzzz_z_zz[i] = t_0_yzzz_0_zzz[i] - rcd_z[i] * t_0_yzzz_0_zz[i];

        t_0_yzzz_z_yz[i] = t_0_yzzz_0_yzz[i] - rcd_z[i] * t_0_yzzz_0_yz[i];

        t_0_yzzz_z_yy[i] = t_0_yzzz_0_yyz[i] - rcd_z[i] * t_0_yzzz_0_yy[i];

        t_0_yzzz_z_xz[i] = t_0_yzzz_0_xzz[i] - rcd_z[i] * t_0_yzzz_0_xz[i];

        t_0_yzzz_z_xy[i] = t_0_yzzz_0_xyz[i] - rcd_z[i] * t_0_yzzz_0_xy[i];

        t_0_yzzz_z_xx[i] = t_0_yzzz_0_xxz[i] - rcd_z[i] * t_0_yzzz_0_xx[i];

        t_0_yzzz_y_zz[i] = t_0_yzzz_0_yzz[i] - rcd_y[i] * t_0_yzzz_0_zz[i];

        t_0_yzzz_y_yz[i] = t_0_yzzz_0_yyz[i] - rcd_y[i] * t_0_yzzz_0_yz[i];

        t_0_yzzz_y_yy[i] = t_0_yzzz_0_yyy[i] - rcd_y[i] * t_0_yzzz_0_yy[i];

        t_0_yzzz_y_xz[i] = t_0_yzzz_0_xyz[i] - rcd_y[i] * t_0_yzzz_0_xz[i];

        t_0_yzzz_y_xy[i] = t_0_yzzz_0_xyy[i] - rcd_y[i] * t_0_yzzz_0_xy[i];

        t_0_yzzz_y_xx[i] = t_0_yzzz_0_xxy[i] - rcd_y[i] * t_0_yzzz_0_xx[i];

        t_0_yzzz_x_zz[i] = t_0_yzzz_0_xzz[i] - rcd_x[i] * t_0_yzzz_0_zz[i];

        t_0_yzzz_x_yz[i] = t_0_yzzz_0_xyz[i] - rcd_x[i] * t_0_yzzz_0_yz[i];

        t_0_yzzz_x_yy[i] = t_0_yzzz_0_xyy[i] - rcd_x[i] * t_0_yzzz_0_yy[i];

        t_0_yzzz_x_xz[i] = t_0_yzzz_0_xxz[i] - rcd_x[i] * t_0_yzzz_0_xz[i];

        t_0_yzzz_x_xy[i] = t_0_yzzz_0_xxy[i] - rcd_x[i] * t_0_yzzz_0_xy[i];

        t_0_yzzz_x_xx[i] = t_0_yzzz_0_xxx[i] - rcd_x[i] * t_0_yzzz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyyz_0_xx, t_0_yyyz_0_xxx, t_0_yyyz_0_xxy,\
                           t_0_yyyz_0_xxz, t_0_yyyz_0_xy, t_0_yyyz_0_xyy, t_0_yyyz_0_xyz,\
                           t_0_yyyz_0_xz, t_0_yyyz_0_xzz, t_0_yyyz_0_yy, t_0_yyyz_0_yyy,\
                           t_0_yyyz_0_yyz, t_0_yyyz_0_yz, t_0_yyyz_0_yzz, t_0_yyyz_0_zz,\
                           t_0_yyyz_0_zzz, t_0_yyyz_x_xx, t_0_yyyz_x_xy, t_0_yyyz_x_xz,\
                           t_0_yyyz_x_yy, t_0_yyyz_x_yz, t_0_yyyz_x_zz, t_0_yyyz_y_xx,\
                           t_0_yyyz_y_xy, t_0_yyyz_y_xz, t_0_yyyz_y_yy, t_0_yyyz_y_yz,\
                           t_0_yyyz_y_zz, t_0_yyyz_z_xx, t_0_yyyz_z_xy, t_0_yyyz_z_xz,\
                           t_0_yyyz_z_yy, t_0_yyyz_z_yz, t_0_yyyz_z_zz, t_0_yyzz_0_xx,\
                           t_0_yyzz_0_xxx, t_0_yyzz_0_xxy, t_0_yyzz_0_xxz, t_0_yyzz_0_xy,\
                           t_0_yyzz_0_xyy, t_0_yyzz_0_xyz, t_0_yyzz_0_xz, t_0_yyzz_0_xzz,\
                           t_0_yyzz_0_yy, t_0_yyzz_0_yyy, t_0_yyzz_0_yyz, t_0_yyzz_0_yz,\
                           t_0_yyzz_0_yzz, t_0_yyzz_0_zz, t_0_yyzz_0_zzz, t_0_yyzz_x_xx,\
                           t_0_yyzz_x_xy, t_0_yyzz_x_xz, t_0_yyzz_x_yy, t_0_yyzz_x_yz,\
                           t_0_yyzz_x_zz, t_0_yyzz_y_xx, t_0_yyzz_y_xy, t_0_yyzz_y_xz,\
                           t_0_yyzz_y_yy, t_0_yyzz_y_yz, t_0_yyzz_y_zz, t_0_yyzz_z_xx,\
                           t_0_yyzz_z_xy, t_0_yyzz_z_xz, t_0_yyzz_z_yy, t_0_yyzz_z_yz,\
                           t_0_yyzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyzz_z_zz[i] = t_0_yyzz_0_zzz[i] - rcd_z[i] * t_0_yyzz_0_zz[i];

        t_0_yyzz_z_yz[i] = t_0_yyzz_0_yzz[i] - rcd_z[i] * t_0_yyzz_0_yz[i];

        t_0_yyzz_z_yy[i] = t_0_yyzz_0_yyz[i] - rcd_z[i] * t_0_yyzz_0_yy[i];

        t_0_yyzz_z_xz[i] = t_0_yyzz_0_xzz[i] - rcd_z[i] * t_0_yyzz_0_xz[i];

        t_0_yyzz_z_xy[i] = t_0_yyzz_0_xyz[i] - rcd_z[i] * t_0_yyzz_0_xy[i];

        t_0_yyzz_z_xx[i] = t_0_yyzz_0_xxz[i] - rcd_z[i] * t_0_yyzz_0_xx[i];

        t_0_yyzz_y_zz[i] = t_0_yyzz_0_yzz[i] - rcd_y[i] * t_0_yyzz_0_zz[i];

        t_0_yyzz_y_yz[i] = t_0_yyzz_0_yyz[i] - rcd_y[i] * t_0_yyzz_0_yz[i];

        t_0_yyzz_y_yy[i] = t_0_yyzz_0_yyy[i] - rcd_y[i] * t_0_yyzz_0_yy[i];

        t_0_yyzz_y_xz[i] = t_0_yyzz_0_xyz[i] - rcd_y[i] * t_0_yyzz_0_xz[i];

        t_0_yyzz_y_xy[i] = t_0_yyzz_0_xyy[i] - rcd_y[i] * t_0_yyzz_0_xy[i];

        t_0_yyzz_y_xx[i] = t_0_yyzz_0_xxy[i] - rcd_y[i] * t_0_yyzz_0_xx[i];

        t_0_yyzz_x_zz[i] = t_0_yyzz_0_xzz[i] - rcd_x[i] * t_0_yyzz_0_zz[i];

        t_0_yyzz_x_yz[i] = t_0_yyzz_0_xyz[i] - rcd_x[i] * t_0_yyzz_0_yz[i];

        t_0_yyzz_x_yy[i] = t_0_yyzz_0_xyy[i] - rcd_x[i] * t_0_yyzz_0_yy[i];

        t_0_yyzz_x_xz[i] = t_0_yyzz_0_xxz[i] - rcd_x[i] * t_0_yyzz_0_xz[i];

        t_0_yyzz_x_xy[i] = t_0_yyzz_0_xxy[i] - rcd_x[i] * t_0_yyzz_0_xy[i];

        t_0_yyzz_x_xx[i] = t_0_yyzz_0_xxx[i] - rcd_x[i] * t_0_yyzz_0_xx[i];

        t_0_yyyz_z_zz[i] = t_0_yyyz_0_zzz[i] - rcd_z[i] * t_0_yyyz_0_zz[i];

        t_0_yyyz_z_yz[i] = t_0_yyyz_0_yzz[i] - rcd_z[i] * t_0_yyyz_0_yz[i];

        t_0_yyyz_z_yy[i] = t_0_yyyz_0_yyz[i] - rcd_z[i] * t_0_yyyz_0_yy[i];

        t_0_yyyz_z_xz[i] = t_0_yyyz_0_xzz[i] - rcd_z[i] * t_0_yyyz_0_xz[i];

        t_0_yyyz_z_xy[i] = t_0_yyyz_0_xyz[i] - rcd_z[i] * t_0_yyyz_0_xy[i];

        t_0_yyyz_z_xx[i] = t_0_yyyz_0_xxz[i] - rcd_z[i] * t_0_yyyz_0_xx[i];

        t_0_yyyz_y_zz[i] = t_0_yyyz_0_yzz[i] - rcd_y[i] * t_0_yyyz_0_zz[i];

        t_0_yyyz_y_yz[i] = t_0_yyyz_0_yyz[i] - rcd_y[i] * t_0_yyyz_0_yz[i];

        t_0_yyyz_y_yy[i] = t_0_yyyz_0_yyy[i] - rcd_y[i] * t_0_yyyz_0_yy[i];

        t_0_yyyz_y_xz[i] = t_0_yyyz_0_xyz[i] - rcd_y[i] * t_0_yyyz_0_xz[i];

        t_0_yyyz_y_xy[i] = t_0_yyyz_0_xyy[i] - rcd_y[i] * t_0_yyyz_0_xy[i];

        t_0_yyyz_y_xx[i] = t_0_yyyz_0_xxy[i] - rcd_y[i] * t_0_yyyz_0_xx[i];

        t_0_yyyz_x_zz[i] = t_0_yyyz_0_xzz[i] - rcd_x[i] * t_0_yyyz_0_zz[i];

        t_0_yyyz_x_yz[i] = t_0_yyyz_0_xyz[i] - rcd_x[i] * t_0_yyyz_0_yz[i];

        t_0_yyyz_x_yy[i] = t_0_yyyz_0_xyy[i] - rcd_x[i] * t_0_yyyz_0_yy[i];

        t_0_yyyz_x_xz[i] = t_0_yyyz_0_xxz[i] - rcd_x[i] * t_0_yyyz_0_xz[i];

        t_0_yyyz_x_xy[i] = t_0_yyyz_0_xxy[i] - rcd_x[i] * t_0_yyyz_0_xy[i];

        t_0_yyyz_x_xx[i] = t_0_yyyz_0_xxx[i] - rcd_x[i] * t_0_yyyz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xzzz_0_xx, t_0_xzzz_0_xxx, t_0_xzzz_0_xxy,\
                           t_0_xzzz_0_xxz, t_0_xzzz_0_xy, t_0_xzzz_0_xyy, t_0_xzzz_0_xyz,\
                           t_0_xzzz_0_xz, t_0_xzzz_0_xzz, t_0_xzzz_0_yy, t_0_xzzz_0_yyy,\
                           t_0_xzzz_0_yyz, t_0_xzzz_0_yz, t_0_xzzz_0_yzz, t_0_xzzz_0_zz,\
                           t_0_xzzz_0_zzz, t_0_xzzz_x_xx, t_0_xzzz_x_xy, t_0_xzzz_x_xz,\
                           t_0_xzzz_x_yy, t_0_xzzz_x_yz, t_0_xzzz_x_zz, t_0_xzzz_y_xx,\
                           t_0_xzzz_y_xy, t_0_xzzz_y_xz, t_0_xzzz_y_yy, t_0_xzzz_y_yz,\
                           t_0_xzzz_y_zz, t_0_xzzz_z_xx, t_0_xzzz_z_xy, t_0_xzzz_z_xz,\
                           t_0_xzzz_z_yy, t_0_xzzz_z_yz, t_0_xzzz_z_zz, t_0_yyyy_0_xx,\
                           t_0_yyyy_0_xxx, t_0_yyyy_0_xxy, t_0_yyyy_0_xxz, t_0_yyyy_0_xy,\
                           t_0_yyyy_0_xyy, t_0_yyyy_0_xyz, t_0_yyyy_0_xz, t_0_yyyy_0_xzz,\
                           t_0_yyyy_0_yy, t_0_yyyy_0_yyy, t_0_yyyy_0_yyz, t_0_yyyy_0_yz,\
                           t_0_yyyy_0_yzz, t_0_yyyy_0_zz, t_0_yyyy_0_zzz, t_0_yyyy_x_xx,\
                           t_0_yyyy_x_xy, t_0_yyyy_x_xz, t_0_yyyy_x_yy, t_0_yyyy_x_yz,\
                           t_0_yyyy_x_zz, t_0_yyyy_y_xx, t_0_yyyy_y_xy, t_0_yyyy_y_xz,\
                           t_0_yyyy_y_yy, t_0_yyyy_y_yz, t_0_yyyy_y_zz, t_0_yyyy_z_xx,\
                           t_0_yyyy_z_xy, t_0_yyyy_z_xz, t_0_yyyy_z_yy, t_0_yyyy_z_yz,\
                           t_0_yyyy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyyy_z_zz[i] = t_0_yyyy_0_zzz[i] - rcd_z[i] * t_0_yyyy_0_zz[i];

        t_0_yyyy_z_yz[i] = t_0_yyyy_0_yzz[i] - rcd_z[i] * t_0_yyyy_0_yz[i];

        t_0_yyyy_z_yy[i] = t_0_yyyy_0_yyz[i] - rcd_z[i] * t_0_yyyy_0_yy[i];

        t_0_yyyy_z_xz[i] = t_0_yyyy_0_xzz[i] - rcd_z[i] * t_0_yyyy_0_xz[i];

        t_0_yyyy_z_xy[i] = t_0_yyyy_0_xyz[i] - rcd_z[i] * t_0_yyyy_0_xy[i];

        t_0_yyyy_z_xx[i] = t_0_yyyy_0_xxz[i] - rcd_z[i] * t_0_yyyy_0_xx[i];

        t_0_yyyy_y_zz[i] = t_0_yyyy_0_yzz[i] - rcd_y[i] * t_0_yyyy_0_zz[i];

        t_0_yyyy_y_yz[i] = t_0_yyyy_0_yyz[i] - rcd_y[i] * t_0_yyyy_0_yz[i];

        t_0_yyyy_y_yy[i] = t_0_yyyy_0_yyy[i] - rcd_y[i] * t_0_yyyy_0_yy[i];

        t_0_yyyy_y_xz[i] = t_0_yyyy_0_xyz[i] - rcd_y[i] * t_0_yyyy_0_xz[i];

        t_0_yyyy_y_xy[i] = t_0_yyyy_0_xyy[i] - rcd_y[i] * t_0_yyyy_0_xy[i];

        t_0_yyyy_y_xx[i] = t_0_yyyy_0_xxy[i] - rcd_y[i] * t_0_yyyy_0_xx[i];

        t_0_yyyy_x_zz[i] = t_0_yyyy_0_xzz[i] - rcd_x[i] * t_0_yyyy_0_zz[i];

        t_0_yyyy_x_yz[i] = t_0_yyyy_0_xyz[i] - rcd_x[i] * t_0_yyyy_0_yz[i];

        t_0_yyyy_x_yy[i] = t_0_yyyy_0_xyy[i] - rcd_x[i] * t_0_yyyy_0_yy[i];

        t_0_yyyy_x_xz[i] = t_0_yyyy_0_xxz[i] - rcd_x[i] * t_0_yyyy_0_xz[i];

        t_0_yyyy_x_xy[i] = t_0_yyyy_0_xxy[i] - rcd_x[i] * t_0_yyyy_0_xy[i];

        t_0_yyyy_x_xx[i] = t_0_yyyy_0_xxx[i] - rcd_x[i] * t_0_yyyy_0_xx[i];

        t_0_xzzz_z_zz[i] = t_0_xzzz_0_zzz[i] - rcd_z[i] * t_0_xzzz_0_zz[i];

        t_0_xzzz_z_yz[i] = t_0_xzzz_0_yzz[i] - rcd_z[i] * t_0_xzzz_0_yz[i];

        t_0_xzzz_z_yy[i] = t_0_xzzz_0_yyz[i] - rcd_z[i] * t_0_xzzz_0_yy[i];

        t_0_xzzz_z_xz[i] = t_0_xzzz_0_xzz[i] - rcd_z[i] * t_0_xzzz_0_xz[i];

        t_0_xzzz_z_xy[i] = t_0_xzzz_0_xyz[i] - rcd_z[i] * t_0_xzzz_0_xy[i];

        t_0_xzzz_z_xx[i] = t_0_xzzz_0_xxz[i] - rcd_z[i] * t_0_xzzz_0_xx[i];

        t_0_xzzz_y_zz[i] = t_0_xzzz_0_yzz[i] - rcd_y[i] * t_0_xzzz_0_zz[i];

        t_0_xzzz_y_yz[i] = t_0_xzzz_0_yyz[i] - rcd_y[i] * t_0_xzzz_0_yz[i];

        t_0_xzzz_y_yy[i] = t_0_xzzz_0_yyy[i] - rcd_y[i] * t_0_xzzz_0_yy[i];

        t_0_xzzz_y_xz[i] = t_0_xzzz_0_xyz[i] - rcd_y[i] * t_0_xzzz_0_xz[i];

        t_0_xzzz_y_xy[i] = t_0_xzzz_0_xyy[i] - rcd_y[i] * t_0_xzzz_0_xy[i];

        t_0_xzzz_y_xx[i] = t_0_xzzz_0_xxy[i] - rcd_y[i] * t_0_xzzz_0_xx[i];

        t_0_xzzz_x_zz[i] = t_0_xzzz_0_xzz[i] - rcd_x[i] * t_0_xzzz_0_zz[i];

        t_0_xzzz_x_yz[i] = t_0_xzzz_0_xyz[i] - rcd_x[i] * t_0_xzzz_0_yz[i];

        t_0_xzzz_x_yy[i] = t_0_xzzz_0_xyy[i] - rcd_x[i] * t_0_xzzz_0_yy[i];

        t_0_xzzz_x_xz[i] = t_0_xzzz_0_xxz[i] - rcd_x[i] * t_0_xzzz_0_xz[i];

        t_0_xzzz_x_xy[i] = t_0_xzzz_0_xxy[i] - rcd_x[i] * t_0_xzzz_0_xy[i];

        t_0_xzzz_x_xx[i] = t_0_xzzz_0_xxx[i] - rcd_x[i] * t_0_xzzz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyz_0_xx, t_0_xyyz_0_xxx, t_0_xyyz_0_xxy,\
                           t_0_xyyz_0_xxz, t_0_xyyz_0_xy, t_0_xyyz_0_xyy, t_0_xyyz_0_xyz,\
                           t_0_xyyz_0_xz, t_0_xyyz_0_xzz, t_0_xyyz_0_yy, t_0_xyyz_0_yyy,\
                           t_0_xyyz_0_yyz, t_0_xyyz_0_yz, t_0_xyyz_0_yzz, t_0_xyyz_0_zz,\
                           t_0_xyyz_0_zzz, t_0_xyyz_x_xx, t_0_xyyz_x_xy, t_0_xyyz_x_xz,\
                           t_0_xyyz_x_yy, t_0_xyyz_x_yz, t_0_xyyz_x_zz, t_0_xyyz_y_xx,\
                           t_0_xyyz_y_xy, t_0_xyyz_y_xz, t_0_xyyz_y_yy, t_0_xyyz_y_yz,\
                           t_0_xyyz_y_zz, t_0_xyyz_z_xx, t_0_xyyz_z_xy, t_0_xyyz_z_xz,\
                           t_0_xyyz_z_yy, t_0_xyyz_z_yz, t_0_xyyz_z_zz, t_0_xyzz_0_xx,\
                           t_0_xyzz_0_xxx, t_0_xyzz_0_xxy, t_0_xyzz_0_xxz, t_0_xyzz_0_xy,\
                           t_0_xyzz_0_xyy, t_0_xyzz_0_xyz, t_0_xyzz_0_xz, t_0_xyzz_0_xzz,\
                           t_0_xyzz_0_yy, t_0_xyzz_0_yyy, t_0_xyzz_0_yyz, t_0_xyzz_0_yz,\
                           t_0_xyzz_0_yzz, t_0_xyzz_0_zz, t_0_xyzz_0_zzz, t_0_xyzz_x_xx,\
                           t_0_xyzz_x_xy, t_0_xyzz_x_xz, t_0_xyzz_x_yy, t_0_xyzz_x_yz,\
                           t_0_xyzz_x_zz, t_0_xyzz_y_xx, t_0_xyzz_y_xy, t_0_xyzz_y_xz,\
                           t_0_xyzz_y_yy, t_0_xyzz_y_yz, t_0_xyzz_y_zz, t_0_xyzz_z_xx,\
                           t_0_xyzz_z_xy, t_0_xyzz_z_xz, t_0_xyzz_z_yy, t_0_xyzz_z_yz,\
                           t_0_xyzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyzz_z_zz[i] = t_0_xyzz_0_zzz[i] - rcd_z[i] * t_0_xyzz_0_zz[i];

        t_0_xyzz_z_yz[i] = t_0_xyzz_0_yzz[i] - rcd_z[i] * t_0_xyzz_0_yz[i];

        t_0_xyzz_z_yy[i] = t_0_xyzz_0_yyz[i] - rcd_z[i] * t_0_xyzz_0_yy[i];

        t_0_xyzz_z_xz[i] = t_0_xyzz_0_xzz[i] - rcd_z[i] * t_0_xyzz_0_xz[i];

        t_0_xyzz_z_xy[i] = t_0_xyzz_0_xyz[i] - rcd_z[i] * t_0_xyzz_0_xy[i];

        t_0_xyzz_z_xx[i] = t_0_xyzz_0_xxz[i] - rcd_z[i] * t_0_xyzz_0_xx[i];

        t_0_xyzz_y_zz[i] = t_0_xyzz_0_yzz[i] - rcd_y[i] * t_0_xyzz_0_zz[i];

        t_0_xyzz_y_yz[i] = t_0_xyzz_0_yyz[i] - rcd_y[i] * t_0_xyzz_0_yz[i];

        t_0_xyzz_y_yy[i] = t_0_xyzz_0_yyy[i] - rcd_y[i] * t_0_xyzz_0_yy[i];

        t_0_xyzz_y_xz[i] = t_0_xyzz_0_xyz[i] - rcd_y[i] * t_0_xyzz_0_xz[i];

        t_0_xyzz_y_xy[i] = t_0_xyzz_0_xyy[i] - rcd_y[i] * t_0_xyzz_0_xy[i];

        t_0_xyzz_y_xx[i] = t_0_xyzz_0_xxy[i] - rcd_y[i] * t_0_xyzz_0_xx[i];

        t_0_xyzz_x_zz[i] = t_0_xyzz_0_xzz[i] - rcd_x[i] * t_0_xyzz_0_zz[i];

        t_0_xyzz_x_yz[i] = t_0_xyzz_0_xyz[i] - rcd_x[i] * t_0_xyzz_0_yz[i];

        t_0_xyzz_x_yy[i] = t_0_xyzz_0_xyy[i] - rcd_x[i] * t_0_xyzz_0_yy[i];

        t_0_xyzz_x_xz[i] = t_0_xyzz_0_xxz[i] - rcd_x[i] * t_0_xyzz_0_xz[i];

        t_0_xyzz_x_xy[i] = t_0_xyzz_0_xxy[i] - rcd_x[i] * t_0_xyzz_0_xy[i];

        t_0_xyzz_x_xx[i] = t_0_xyzz_0_xxx[i] - rcd_x[i] * t_0_xyzz_0_xx[i];

        t_0_xyyz_z_zz[i] = t_0_xyyz_0_zzz[i] - rcd_z[i] * t_0_xyyz_0_zz[i];

        t_0_xyyz_z_yz[i] = t_0_xyyz_0_yzz[i] - rcd_z[i] * t_0_xyyz_0_yz[i];

        t_0_xyyz_z_yy[i] = t_0_xyyz_0_yyz[i] - rcd_z[i] * t_0_xyyz_0_yy[i];

        t_0_xyyz_z_xz[i] = t_0_xyyz_0_xzz[i] - rcd_z[i] * t_0_xyyz_0_xz[i];

        t_0_xyyz_z_xy[i] = t_0_xyyz_0_xyz[i] - rcd_z[i] * t_0_xyyz_0_xy[i];

        t_0_xyyz_z_xx[i] = t_0_xyyz_0_xxz[i] - rcd_z[i] * t_0_xyyz_0_xx[i];

        t_0_xyyz_y_zz[i] = t_0_xyyz_0_yzz[i] - rcd_y[i] * t_0_xyyz_0_zz[i];

        t_0_xyyz_y_yz[i] = t_0_xyyz_0_yyz[i] - rcd_y[i] * t_0_xyyz_0_yz[i];

        t_0_xyyz_y_yy[i] = t_0_xyyz_0_yyy[i] - rcd_y[i] * t_0_xyyz_0_yy[i];

        t_0_xyyz_y_xz[i] = t_0_xyyz_0_xyz[i] - rcd_y[i] * t_0_xyyz_0_xz[i];

        t_0_xyyz_y_xy[i] = t_0_xyyz_0_xyy[i] - rcd_y[i] * t_0_xyyz_0_xy[i];

        t_0_xyyz_y_xx[i] = t_0_xyyz_0_xxy[i] - rcd_y[i] * t_0_xyyz_0_xx[i];

        t_0_xyyz_x_zz[i] = t_0_xyyz_0_xzz[i] - rcd_x[i] * t_0_xyyz_0_zz[i];

        t_0_xyyz_x_yz[i] = t_0_xyyz_0_xyz[i] - rcd_x[i] * t_0_xyyz_0_yz[i];

        t_0_xyyz_x_yy[i] = t_0_xyyz_0_xyy[i] - rcd_x[i] * t_0_xyyz_0_yy[i];

        t_0_xyyz_x_xz[i] = t_0_xyyz_0_xxz[i] - rcd_x[i] * t_0_xyyz_0_xz[i];

        t_0_xyyz_x_xy[i] = t_0_xyyz_0_xxy[i] - rcd_x[i] * t_0_xyyz_0_xy[i];

        t_0_xyyz_x_xx[i] = t_0_xyyz_0_xxx[i] - rcd_x[i] * t_0_xyyz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxzz_0_xx, t_0_xxzz_0_xxx, t_0_xxzz_0_xxy,\
                           t_0_xxzz_0_xxz, t_0_xxzz_0_xy, t_0_xxzz_0_xyy, t_0_xxzz_0_xyz,\
                           t_0_xxzz_0_xz, t_0_xxzz_0_xzz, t_0_xxzz_0_yy, t_0_xxzz_0_yyy,\
                           t_0_xxzz_0_yyz, t_0_xxzz_0_yz, t_0_xxzz_0_yzz, t_0_xxzz_0_zz,\
                           t_0_xxzz_0_zzz, t_0_xxzz_x_xx, t_0_xxzz_x_xy, t_0_xxzz_x_xz,\
                           t_0_xxzz_x_yy, t_0_xxzz_x_yz, t_0_xxzz_x_zz, t_0_xxzz_y_xx,\
                           t_0_xxzz_y_xy, t_0_xxzz_y_xz, t_0_xxzz_y_yy, t_0_xxzz_y_yz,\
                           t_0_xxzz_y_zz, t_0_xxzz_z_xx, t_0_xxzz_z_xy, t_0_xxzz_z_xz,\
                           t_0_xxzz_z_yy, t_0_xxzz_z_yz, t_0_xxzz_z_zz, t_0_xyyy_0_xx,\
                           t_0_xyyy_0_xxx, t_0_xyyy_0_xxy, t_0_xyyy_0_xxz, t_0_xyyy_0_xy,\
                           t_0_xyyy_0_xyy, t_0_xyyy_0_xyz, t_0_xyyy_0_xz, t_0_xyyy_0_xzz,\
                           t_0_xyyy_0_yy, t_0_xyyy_0_yyy, t_0_xyyy_0_yyz, t_0_xyyy_0_yz,\
                           t_0_xyyy_0_yzz, t_0_xyyy_0_zz, t_0_xyyy_0_zzz, t_0_xyyy_x_xx,\
                           t_0_xyyy_x_xy, t_0_xyyy_x_xz, t_0_xyyy_x_yy, t_0_xyyy_x_yz,\
                           t_0_xyyy_x_zz, t_0_xyyy_y_xx, t_0_xyyy_y_xy, t_0_xyyy_y_xz,\
                           t_0_xyyy_y_yy, t_0_xyyy_y_yz, t_0_xyyy_y_zz, t_0_xyyy_z_xx,\
                           t_0_xyyy_z_xy, t_0_xyyy_z_xz, t_0_xyyy_z_yy, t_0_xyyy_z_yz,\
                           t_0_xyyy_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyy_z_zz[i] = t_0_xyyy_0_zzz[i] - rcd_z[i] * t_0_xyyy_0_zz[i];

        t_0_xyyy_z_yz[i] = t_0_xyyy_0_yzz[i] - rcd_z[i] * t_0_xyyy_0_yz[i];

        t_0_xyyy_z_yy[i] = t_0_xyyy_0_yyz[i] - rcd_z[i] * t_0_xyyy_0_yy[i];

        t_0_xyyy_z_xz[i] = t_0_xyyy_0_xzz[i] - rcd_z[i] * t_0_xyyy_0_xz[i];

        t_0_xyyy_z_xy[i] = t_0_xyyy_0_xyz[i] - rcd_z[i] * t_0_xyyy_0_xy[i];

        t_0_xyyy_z_xx[i] = t_0_xyyy_0_xxz[i] - rcd_z[i] * t_0_xyyy_0_xx[i];

        t_0_xyyy_y_zz[i] = t_0_xyyy_0_yzz[i] - rcd_y[i] * t_0_xyyy_0_zz[i];

        t_0_xyyy_y_yz[i] = t_0_xyyy_0_yyz[i] - rcd_y[i] * t_0_xyyy_0_yz[i];

        t_0_xyyy_y_yy[i] = t_0_xyyy_0_yyy[i] - rcd_y[i] * t_0_xyyy_0_yy[i];

        t_0_xyyy_y_xz[i] = t_0_xyyy_0_xyz[i] - rcd_y[i] * t_0_xyyy_0_xz[i];

        t_0_xyyy_y_xy[i] = t_0_xyyy_0_xyy[i] - rcd_y[i] * t_0_xyyy_0_xy[i];

        t_0_xyyy_y_xx[i] = t_0_xyyy_0_xxy[i] - rcd_y[i] * t_0_xyyy_0_xx[i];

        t_0_xyyy_x_zz[i] = t_0_xyyy_0_xzz[i] - rcd_x[i] * t_0_xyyy_0_zz[i];

        t_0_xyyy_x_yz[i] = t_0_xyyy_0_xyz[i] - rcd_x[i] * t_0_xyyy_0_yz[i];

        t_0_xyyy_x_yy[i] = t_0_xyyy_0_xyy[i] - rcd_x[i] * t_0_xyyy_0_yy[i];

        t_0_xyyy_x_xz[i] = t_0_xyyy_0_xxz[i] - rcd_x[i] * t_0_xyyy_0_xz[i];

        t_0_xyyy_x_xy[i] = t_0_xyyy_0_xxy[i] - rcd_x[i] * t_0_xyyy_0_xy[i];

        t_0_xyyy_x_xx[i] = t_0_xyyy_0_xxx[i] - rcd_x[i] * t_0_xyyy_0_xx[i];

        t_0_xxzz_z_zz[i] = t_0_xxzz_0_zzz[i] - rcd_z[i] * t_0_xxzz_0_zz[i];

        t_0_xxzz_z_yz[i] = t_0_xxzz_0_yzz[i] - rcd_z[i] * t_0_xxzz_0_yz[i];

        t_0_xxzz_z_yy[i] = t_0_xxzz_0_yyz[i] - rcd_z[i] * t_0_xxzz_0_yy[i];

        t_0_xxzz_z_xz[i] = t_0_xxzz_0_xzz[i] - rcd_z[i] * t_0_xxzz_0_xz[i];

        t_0_xxzz_z_xy[i] = t_0_xxzz_0_xyz[i] - rcd_z[i] * t_0_xxzz_0_xy[i];

        t_0_xxzz_z_xx[i] = t_0_xxzz_0_xxz[i] - rcd_z[i] * t_0_xxzz_0_xx[i];

        t_0_xxzz_y_zz[i] = t_0_xxzz_0_yzz[i] - rcd_y[i] * t_0_xxzz_0_zz[i];

        t_0_xxzz_y_yz[i] = t_0_xxzz_0_yyz[i] - rcd_y[i] * t_0_xxzz_0_yz[i];

        t_0_xxzz_y_yy[i] = t_0_xxzz_0_yyy[i] - rcd_y[i] * t_0_xxzz_0_yy[i];

        t_0_xxzz_y_xz[i] = t_0_xxzz_0_xyz[i] - rcd_y[i] * t_0_xxzz_0_xz[i];

        t_0_xxzz_y_xy[i] = t_0_xxzz_0_xyy[i] - rcd_y[i] * t_0_xxzz_0_xy[i];

        t_0_xxzz_y_xx[i] = t_0_xxzz_0_xxy[i] - rcd_y[i] * t_0_xxzz_0_xx[i];

        t_0_xxzz_x_zz[i] = t_0_xxzz_0_xzz[i] - rcd_x[i] * t_0_xxzz_0_zz[i];

        t_0_xxzz_x_yz[i] = t_0_xxzz_0_xyz[i] - rcd_x[i] * t_0_xxzz_0_yz[i];

        t_0_xxzz_x_yy[i] = t_0_xxzz_0_xyy[i] - rcd_x[i] * t_0_xxzz_0_yy[i];

        t_0_xxzz_x_xz[i] = t_0_xxzz_0_xxz[i] - rcd_x[i] * t_0_xxzz_0_xz[i];

        t_0_xxzz_x_xy[i] = t_0_xxzz_0_xxy[i] - rcd_x[i] * t_0_xxzz_0_xy[i];

        t_0_xxzz_x_xx[i] = t_0_xxzz_0_xxx[i] - rcd_x[i] * t_0_xxzz_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyy_0_xx, t_0_xxyy_0_xxx, t_0_xxyy_0_xxy,\
                           t_0_xxyy_0_xxz, t_0_xxyy_0_xy, t_0_xxyy_0_xyy, t_0_xxyy_0_xyz,\
                           t_0_xxyy_0_xz, t_0_xxyy_0_xzz, t_0_xxyy_0_yy, t_0_xxyy_0_yyy,\
                           t_0_xxyy_0_yyz, t_0_xxyy_0_yz, t_0_xxyy_0_yzz, t_0_xxyy_0_zz,\
                           t_0_xxyy_0_zzz, t_0_xxyy_x_xx, t_0_xxyy_x_xy, t_0_xxyy_x_xz,\
                           t_0_xxyy_x_yy, t_0_xxyy_x_yz, t_0_xxyy_x_zz, t_0_xxyy_y_xx,\
                           t_0_xxyy_y_xy, t_0_xxyy_y_xz, t_0_xxyy_y_yy, t_0_xxyy_y_yz,\
                           t_0_xxyy_y_zz, t_0_xxyy_z_xx, t_0_xxyy_z_xy, t_0_xxyy_z_xz,\
                           t_0_xxyy_z_yy, t_0_xxyy_z_yz, t_0_xxyy_z_zz, t_0_xxyz_0_xx,\
                           t_0_xxyz_0_xxx, t_0_xxyz_0_xxy, t_0_xxyz_0_xxz, t_0_xxyz_0_xy,\
                           t_0_xxyz_0_xyy, t_0_xxyz_0_xyz, t_0_xxyz_0_xz, t_0_xxyz_0_xzz,\
                           t_0_xxyz_0_yy, t_0_xxyz_0_yyy, t_0_xxyz_0_yyz, t_0_xxyz_0_yz,\
                           t_0_xxyz_0_yzz, t_0_xxyz_0_zz, t_0_xxyz_0_zzz, t_0_xxyz_x_xx,\
                           t_0_xxyz_x_xy, t_0_xxyz_x_xz, t_0_xxyz_x_yy, t_0_xxyz_x_yz,\
                           t_0_xxyz_x_zz, t_0_xxyz_y_xx, t_0_xxyz_y_xy, t_0_xxyz_y_xz,\
                           t_0_xxyz_y_yy, t_0_xxyz_y_yz, t_0_xxyz_y_zz, t_0_xxyz_z_xx,\
                           t_0_xxyz_z_xy, t_0_xxyz_z_xz, t_0_xxyz_z_yy, t_0_xxyz_z_yz,\
                           t_0_xxyz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxyz_z_zz[i] = t_0_xxyz_0_zzz[i] - rcd_z[i] * t_0_xxyz_0_zz[i];

        t_0_xxyz_z_yz[i] = t_0_xxyz_0_yzz[i] - rcd_z[i] * t_0_xxyz_0_yz[i];

        t_0_xxyz_z_yy[i] = t_0_xxyz_0_yyz[i] - rcd_z[i] * t_0_xxyz_0_yy[i];

        t_0_xxyz_z_xz[i] = t_0_xxyz_0_xzz[i] - rcd_z[i] * t_0_xxyz_0_xz[i];

        t_0_xxyz_z_xy[i] = t_0_xxyz_0_xyz[i] - rcd_z[i] * t_0_xxyz_0_xy[i];

        t_0_xxyz_z_xx[i] = t_0_xxyz_0_xxz[i] - rcd_z[i] * t_0_xxyz_0_xx[i];

        t_0_xxyz_y_zz[i] = t_0_xxyz_0_yzz[i] - rcd_y[i] * t_0_xxyz_0_zz[i];

        t_0_xxyz_y_yz[i] = t_0_xxyz_0_yyz[i] - rcd_y[i] * t_0_xxyz_0_yz[i];

        t_0_xxyz_y_yy[i] = t_0_xxyz_0_yyy[i] - rcd_y[i] * t_0_xxyz_0_yy[i];

        t_0_xxyz_y_xz[i] = t_0_xxyz_0_xyz[i] - rcd_y[i] * t_0_xxyz_0_xz[i];

        t_0_xxyz_y_xy[i] = t_0_xxyz_0_xyy[i] - rcd_y[i] * t_0_xxyz_0_xy[i];

        t_0_xxyz_y_xx[i] = t_0_xxyz_0_xxy[i] - rcd_y[i] * t_0_xxyz_0_xx[i];

        t_0_xxyz_x_zz[i] = t_0_xxyz_0_xzz[i] - rcd_x[i] * t_0_xxyz_0_zz[i];

        t_0_xxyz_x_yz[i] = t_0_xxyz_0_xyz[i] - rcd_x[i] * t_0_xxyz_0_yz[i];

        t_0_xxyz_x_yy[i] = t_0_xxyz_0_xyy[i] - rcd_x[i] * t_0_xxyz_0_yy[i];

        t_0_xxyz_x_xz[i] = t_0_xxyz_0_xxz[i] - rcd_x[i] * t_0_xxyz_0_xz[i];

        t_0_xxyz_x_xy[i] = t_0_xxyz_0_xxy[i] - rcd_x[i] * t_0_xxyz_0_xy[i];

        t_0_xxyz_x_xx[i] = t_0_xxyz_0_xxx[i] - rcd_x[i] * t_0_xxyz_0_xx[i];

        t_0_xxyy_z_zz[i] = t_0_xxyy_0_zzz[i] - rcd_z[i] * t_0_xxyy_0_zz[i];

        t_0_xxyy_z_yz[i] = t_0_xxyy_0_yzz[i] - rcd_z[i] * t_0_xxyy_0_yz[i];

        t_0_xxyy_z_yy[i] = t_0_xxyy_0_yyz[i] - rcd_z[i] * t_0_xxyy_0_yy[i];

        t_0_xxyy_z_xz[i] = t_0_xxyy_0_xzz[i] - rcd_z[i] * t_0_xxyy_0_xz[i];

        t_0_xxyy_z_xy[i] = t_0_xxyy_0_xyz[i] - rcd_z[i] * t_0_xxyy_0_xy[i];

        t_0_xxyy_z_xx[i] = t_0_xxyy_0_xxz[i] - rcd_z[i] * t_0_xxyy_0_xx[i];

        t_0_xxyy_y_zz[i] = t_0_xxyy_0_yzz[i] - rcd_y[i] * t_0_xxyy_0_zz[i];

        t_0_xxyy_y_yz[i] = t_0_xxyy_0_yyz[i] - rcd_y[i] * t_0_xxyy_0_yz[i];

        t_0_xxyy_y_yy[i] = t_0_xxyy_0_yyy[i] - rcd_y[i] * t_0_xxyy_0_yy[i];

        t_0_xxyy_y_xz[i] = t_0_xxyy_0_xyz[i] - rcd_y[i] * t_0_xxyy_0_xz[i];

        t_0_xxyy_y_xy[i] = t_0_xxyy_0_xyy[i] - rcd_y[i] * t_0_xxyy_0_xy[i];

        t_0_xxyy_y_xx[i] = t_0_xxyy_0_xxy[i] - rcd_y[i] * t_0_xxyy_0_xx[i];

        t_0_xxyy_x_zz[i] = t_0_xxyy_0_xzz[i] - rcd_x[i] * t_0_xxyy_0_zz[i];

        t_0_xxyy_x_yz[i] = t_0_xxyy_0_xyz[i] - rcd_x[i] * t_0_xxyy_0_yz[i];

        t_0_xxyy_x_yy[i] = t_0_xxyy_0_xyy[i] - rcd_x[i] * t_0_xxyy_0_yy[i];

        t_0_xxyy_x_xz[i] = t_0_xxyy_0_xxz[i] - rcd_x[i] * t_0_xxyy_0_xz[i];

        t_0_xxyy_x_xy[i] = t_0_xxyy_0_xxy[i] - rcd_x[i] * t_0_xxyy_0_xy[i];

        t_0_xxyy_x_xx[i] = t_0_xxyy_0_xxx[i] - rcd_x[i] * t_0_xxyy_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxy_0_xx, t_0_xxxy_0_xxx, t_0_xxxy_0_xxy,\
                           t_0_xxxy_0_xxz, t_0_xxxy_0_xy, t_0_xxxy_0_xyy, t_0_xxxy_0_xyz,\
                           t_0_xxxy_0_xz, t_0_xxxy_0_xzz, t_0_xxxy_0_yy, t_0_xxxy_0_yyy,\
                           t_0_xxxy_0_yyz, t_0_xxxy_0_yz, t_0_xxxy_0_yzz, t_0_xxxy_0_zz,\
                           t_0_xxxy_0_zzz, t_0_xxxy_x_xx, t_0_xxxy_x_xy, t_0_xxxy_x_xz,\
                           t_0_xxxy_x_yy, t_0_xxxy_x_yz, t_0_xxxy_x_zz, t_0_xxxy_y_xx,\
                           t_0_xxxy_y_xy, t_0_xxxy_y_xz, t_0_xxxy_y_yy, t_0_xxxy_y_yz,\
                           t_0_xxxy_y_zz, t_0_xxxy_z_xx, t_0_xxxy_z_xy, t_0_xxxy_z_xz,\
                           t_0_xxxy_z_yy, t_0_xxxy_z_yz, t_0_xxxy_z_zz, t_0_xxxz_0_xx,\
                           t_0_xxxz_0_xxx, t_0_xxxz_0_xxy, t_0_xxxz_0_xxz, t_0_xxxz_0_xy,\
                           t_0_xxxz_0_xyy, t_0_xxxz_0_xyz, t_0_xxxz_0_xz, t_0_xxxz_0_xzz,\
                           t_0_xxxz_0_yy, t_0_xxxz_0_yyy, t_0_xxxz_0_yyz, t_0_xxxz_0_yz,\
                           t_0_xxxz_0_yzz, t_0_xxxz_0_zz, t_0_xxxz_0_zzz, t_0_xxxz_x_xx,\
                           t_0_xxxz_x_xy, t_0_xxxz_x_xz, t_0_xxxz_x_yy, t_0_xxxz_x_yz,\
                           t_0_xxxz_x_zz, t_0_xxxz_y_xx, t_0_xxxz_y_xy, t_0_xxxz_y_xz,\
                           t_0_xxxz_y_yy, t_0_xxxz_y_yz, t_0_xxxz_y_zz, t_0_xxxz_z_xx,\
                           t_0_xxxz_z_xy, t_0_xxxz_z_xz, t_0_xxxz_z_yy, t_0_xxxz_z_yz,\
                           t_0_xxxz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxz_z_zz[i] = t_0_xxxz_0_zzz[i] - rcd_z[i] * t_0_xxxz_0_zz[i];

        t_0_xxxz_z_yz[i] = t_0_xxxz_0_yzz[i] - rcd_z[i] * t_0_xxxz_0_yz[i];

        t_0_xxxz_z_yy[i] = t_0_xxxz_0_yyz[i] - rcd_z[i] * t_0_xxxz_0_yy[i];

        t_0_xxxz_z_xz[i] = t_0_xxxz_0_xzz[i] - rcd_z[i] * t_0_xxxz_0_xz[i];

        t_0_xxxz_z_xy[i] = t_0_xxxz_0_xyz[i] - rcd_z[i] * t_0_xxxz_0_xy[i];

        t_0_xxxz_z_xx[i] = t_0_xxxz_0_xxz[i] - rcd_z[i] * t_0_xxxz_0_xx[i];

        t_0_xxxz_y_zz[i] = t_0_xxxz_0_yzz[i] - rcd_y[i] * t_0_xxxz_0_zz[i];

        t_0_xxxz_y_yz[i] = t_0_xxxz_0_yyz[i] - rcd_y[i] * t_0_xxxz_0_yz[i];

        t_0_xxxz_y_yy[i] = t_0_xxxz_0_yyy[i] - rcd_y[i] * t_0_xxxz_0_yy[i];

        t_0_xxxz_y_xz[i] = t_0_xxxz_0_xyz[i] - rcd_y[i] * t_0_xxxz_0_xz[i];

        t_0_xxxz_y_xy[i] = t_0_xxxz_0_xyy[i] - rcd_y[i] * t_0_xxxz_0_xy[i];

        t_0_xxxz_y_xx[i] = t_0_xxxz_0_xxy[i] - rcd_y[i] * t_0_xxxz_0_xx[i];

        t_0_xxxz_x_zz[i] = t_0_xxxz_0_xzz[i] - rcd_x[i] * t_0_xxxz_0_zz[i];

        t_0_xxxz_x_yz[i] = t_0_xxxz_0_xyz[i] - rcd_x[i] * t_0_xxxz_0_yz[i];

        t_0_xxxz_x_yy[i] = t_0_xxxz_0_xyy[i] - rcd_x[i] * t_0_xxxz_0_yy[i];

        t_0_xxxz_x_xz[i] = t_0_xxxz_0_xxz[i] - rcd_x[i] * t_0_xxxz_0_xz[i];

        t_0_xxxz_x_xy[i] = t_0_xxxz_0_xxy[i] - rcd_x[i] * t_0_xxxz_0_xy[i];

        t_0_xxxz_x_xx[i] = t_0_xxxz_0_xxx[i] - rcd_x[i] * t_0_xxxz_0_xx[i];

        t_0_xxxy_z_zz[i] = t_0_xxxy_0_zzz[i] - rcd_z[i] * t_0_xxxy_0_zz[i];

        t_0_xxxy_z_yz[i] = t_0_xxxy_0_yzz[i] - rcd_z[i] * t_0_xxxy_0_yz[i];

        t_0_xxxy_z_yy[i] = t_0_xxxy_0_yyz[i] - rcd_z[i] * t_0_xxxy_0_yy[i];

        t_0_xxxy_z_xz[i] = t_0_xxxy_0_xzz[i] - rcd_z[i] * t_0_xxxy_0_xz[i];

        t_0_xxxy_z_xy[i] = t_0_xxxy_0_xyz[i] - rcd_z[i] * t_0_xxxy_0_xy[i];

        t_0_xxxy_z_xx[i] = t_0_xxxy_0_xxz[i] - rcd_z[i] * t_0_xxxy_0_xx[i];

        t_0_xxxy_y_zz[i] = t_0_xxxy_0_yzz[i] - rcd_y[i] * t_0_xxxy_0_zz[i];

        t_0_xxxy_y_yz[i] = t_0_xxxy_0_yyz[i] - rcd_y[i] * t_0_xxxy_0_yz[i];

        t_0_xxxy_y_yy[i] = t_0_xxxy_0_yyy[i] - rcd_y[i] * t_0_xxxy_0_yy[i];

        t_0_xxxy_y_xz[i] = t_0_xxxy_0_xyz[i] - rcd_y[i] * t_0_xxxy_0_xz[i];

        t_0_xxxy_y_xy[i] = t_0_xxxy_0_xyy[i] - rcd_y[i] * t_0_xxxy_0_xy[i];

        t_0_xxxy_y_xx[i] = t_0_xxxy_0_xxy[i] - rcd_y[i] * t_0_xxxy_0_xx[i];

        t_0_xxxy_x_zz[i] = t_0_xxxy_0_xzz[i] - rcd_x[i] * t_0_xxxy_0_zz[i];

        t_0_xxxy_x_yz[i] = t_0_xxxy_0_xyz[i] - rcd_x[i] * t_0_xxxy_0_yz[i];

        t_0_xxxy_x_yy[i] = t_0_xxxy_0_xyy[i] - rcd_x[i] * t_0_xxxy_0_yy[i];

        t_0_xxxy_x_xz[i] = t_0_xxxy_0_xxz[i] - rcd_x[i] * t_0_xxxy_0_xz[i];

        t_0_xxxy_x_xy[i] = t_0_xxxy_0_xxy[i] - rcd_x[i] * t_0_xxxy_0_xy[i];

        t_0_xxxy_x_xx[i] = t_0_xxxy_0_xxx[i] - rcd_x[i] * t_0_xxxy_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_0_xx, t_0_xxxx_0_xxx, t_0_xxxx_0_xxy,\
                           t_0_xxxx_0_xxz, t_0_xxxx_0_xy, t_0_xxxx_0_xyy, t_0_xxxx_0_xyz,\
                           t_0_xxxx_0_xz, t_0_xxxx_0_xzz, t_0_xxxx_0_yy, t_0_xxxx_0_yyy,\
                           t_0_xxxx_0_yyz, t_0_xxxx_0_yz, t_0_xxxx_0_yzz, t_0_xxxx_0_zz,\
                           t_0_xxxx_0_zzz, t_0_xxxx_x_xx, t_0_xxxx_x_xy, t_0_xxxx_x_xz,\
                           t_0_xxxx_x_yy, t_0_xxxx_x_yz, t_0_xxxx_x_zz, t_0_xxxx_y_xx,\
                           t_0_xxxx_y_xy, t_0_xxxx_y_xz, t_0_xxxx_y_yy, t_0_xxxx_y_yz,\
                           t_0_xxxx_y_zz, t_0_xxxx_z_xx, t_0_xxxx_z_xy, t_0_xxxx_z_xz,\
                           t_0_xxxx_z_yy, t_0_xxxx_z_yz, t_0_xxxx_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxx_z_zz[i] = t_0_xxxx_0_zzz[i] - rcd_z[i] * t_0_xxxx_0_zz[i];

        t_0_xxxx_z_yz[i] = t_0_xxxx_0_yzz[i] - rcd_z[i] * t_0_xxxx_0_yz[i];

        t_0_xxxx_z_yy[i] = t_0_xxxx_0_yyz[i] - rcd_z[i] * t_0_xxxx_0_yy[i];

        t_0_xxxx_z_xz[i] = t_0_xxxx_0_xzz[i] - rcd_z[i] * t_0_xxxx_0_xz[i];

        t_0_xxxx_z_xy[i] = t_0_xxxx_0_xyz[i] - rcd_z[i] * t_0_xxxx_0_xy[i];

        t_0_xxxx_z_xx[i] = t_0_xxxx_0_xxz[i] - rcd_z[i] * t_0_xxxx_0_xx[i];

        t_0_xxxx_y_zz[i] = t_0_xxxx_0_yzz[i] - rcd_y[i] * t_0_xxxx_0_zz[i];

        t_0_xxxx_y_yz[i] = t_0_xxxx_0_yyz[i] - rcd_y[i] * t_0_xxxx_0_yz[i];

        t_0_xxxx_y_yy[i] = t_0_xxxx_0_yyy[i] - rcd_y[i] * t_0_xxxx_0_yy[i];

        t_0_xxxx_y_xz[i] = t_0_xxxx_0_xyz[i] - rcd_y[i] * t_0_xxxx_0_xz[i];

        t_0_xxxx_y_xy[i] = t_0_xxxx_0_xyy[i] - rcd_y[i] * t_0_xxxx_0_xy[i];

        t_0_xxxx_y_xx[i] = t_0_xxxx_0_xxy[i] - rcd_y[i] * t_0_xxxx_0_xx[i];

        t_0_xxxx_x_zz[i] = t_0_xxxx_0_xzz[i] - rcd_x[i] * t_0_xxxx_0_zz[i];

        t_0_xxxx_x_yz[i] = t_0_xxxx_0_xyz[i] - rcd_x[i] * t_0_xxxx_0_yz[i];

        t_0_xxxx_x_yy[i] = t_0_xxxx_0_xyy[i] - rcd_x[i] * t_0_xxxx_0_yy[i];

        t_0_xxxx_x_xz[i] = t_0_xxxx_0_xxz[i] - rcd_x[i] * t_0_xxxx_0_xz[i];

        t_0_xxxx_x_xy[i] = t_0_xxxx_0_xxy[i] - rcd_x[i] * t_0_xxxx_0_xy[i];

        t_0_xxxx_x_xx[i] = t_0_xxxx_0_xxx[i] - rcd_x[i] * t_0_xxxx_0_xx[i];
    }
}


} // derirec namespace
