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
compHostHRRForSGPF_V0(      BufferHostXY<T>&      intsBufferSGPF,
                      const BufferHostX<int32_t>& intsIndexesSGPF,
                      const BufferHostXY<T>&      intsBufferSGSF,
                      const BufferHostX<int32_t>& intsIndexesSGSF,
                      const BufferHostXY<T>&      intsBufferSGSG,
                      const BufferHostX<int32_t>& intsIndexesSGSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SGSG) integral components

    t_0_zzzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(0));

    t_0_zzzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(1));

    t_0_zzzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(2));

    t_0_zzzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(3));

    t_0_zzzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(4));

    t_0_zzzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(5));

    t_0_zzzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(6));

    t_0_zzzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(7));

    t_0_zzzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(8));

    t_0_zzzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(9));

    t_0_zzzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(10));

    t_0_zzzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(11));

    t_0_zzzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(12));

    t_0_zzzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(13));

    t_0_zzzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(14));

    t_0_yzzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(15));

    t_0_yzzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(16));

    t_0_yzzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(17));

    t_0_yzzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(18));

    t_0_yzzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(19));

    t_0_yzzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(20));

    t_0_yzzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(21));

    t_0_yzzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(22));

    t_0_yzzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(23));

    t_0_yzzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(24));

    t_0_yzzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(25));

    t_0_yzzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(26));

    t_0_yzzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(27));

    t_0_yzzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(28));

    t_0_yzzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(29));

    t_0_yyzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(30));

    t_0_yyzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(31));

    t_0_yyzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(32));

    t_0_yyzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(33));

    t_0_yyzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(34));

    t_0_yyzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(35));

    t_0_yyzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(36));

    t_0_yyzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(37));

    t_0_yyzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(38));

    t_0_yyzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(39));

    t_0_yyzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(40));

    t_0_yyzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(41));

    t_0_yyzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(42));

    t_0_yyzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(43));

    t_0_yyzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(44));

    t_0_yyyz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(45));

    t_0_yyyz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(46));

    t_0_yyyz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(47));

    t_0_yyyz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(48));

    t_0_yyyz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(49));

    t_0_yyyz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(50));

    t_0_yyyz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(51));

    t_0_yyyz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(52));

    t_0_yyyz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(53));

    t_0_yyyz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(54));

    t_0_yyyz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(55));

    t_0_yyyz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(56));

    t_0_yyyz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(57));

    t_0_yyyz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(58));

    t_0_yyyz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(59));

    t_0_yyyy_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(60));

    t_0_yyyy_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(61));

    t_0_yyyy_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(62));

    t_0_yyyy_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(63));

    t_0_yyyy_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(64));

    t_0_yyyy_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(65));

    t_0_yyyy_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(66));

    t_0_yyyy_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(67));

    t_0_yyyy_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(68));

    t_0_yyyy_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(69));

    t_0_yyyy_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(70));

    t_0_yyyy_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(71));

    t_0_yyyy_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(72));

    t_0_yyyy_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(73));

    t_0_yyyy_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(74));

    t_0_xzzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(75));

    t_0_xzzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(76));

    t_0_xzzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(77));

    t_0_xzzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(78));

    t_0_xzzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(79));

    t_0_xzzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(80));

    t_0_xzzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(81));

    t_0_xzzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(82));

    t_0_xzzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(83));

    t_0_xzzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(84));

    t_0_xzzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(85));

    t_0_xzzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(86));

    t_0_xzzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(87));

    t_0_xzzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(88));

    t_0_xzzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(89));

    t_0_xyzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(90));

    t_0_xyzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(91));

    t_0_xyzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(92));

    t_0_xyzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(93));

    t_0_xyzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(94));

    t_0_xyzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(95));

    t_0_xyzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(96));

    t_0_xyzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(97));

    t_0_xyzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(98));

    t_0_xyzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(99));

    t_0_xyzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(100));

    t_0_xyzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(101));

    t_0_xyzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(102));

    t_0_xyzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(103));

    t_0_xyzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(104));

    t_0_xyyz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(105));

    t_0_xyyz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(106));

    t_0_xyyz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(107));

    t_0_xyyz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(108));

    t_0_xyyz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(109));

    t_0_xyyz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(110));

    t_0_xyyz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(111));

    t_0_xyyz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(112));

    t_0_xyyz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(113));

    t_0_xyyz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(114));

    t_0_xyyz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(115));

    t_0_xyyz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(116));

    t_0_xyyz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(117));

    t_0_xyyz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(118));

    t_0_xyyz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(119));

    t_0_xyyy_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(120));

    t_0_xyyy_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(121));

    t_0_xyyy_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(122));

    t_0_xyyy_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(123));

    t_0_xyyy_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(124));

    t_0_xyyy_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(125));

    t_0_xyyy_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(126));

    t_0_xyyy_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(127));

    t_0_xyyy_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(128));

    t_0_xyyy_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(129));

    t_0_xyyy_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(130));

    t_0_xyyy_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(131));

    t_0_xyyy_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(132));

    t_0_xyyy_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(133));

    t_0_xyyy_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(134));

    t_0_xxzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(135));

    t_0_xxzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(136));

    t_0_xxzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(137));

    t_0_xxzz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(138));

    t_0_xxzz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(139));

    t_0_xxzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(140));

    t_0_xxzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(141));

    t_0_xxzz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(142));

    t_0_xxzz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(143));

    t_0_xxzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(144));

    t_0_xxzz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(145));

    t_0_xxzz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(146));

    t_0_xxzz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(147));

    t_0_xxzz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(148));

    t_0_xxzz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(149));

    t_0_xxyz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(150));

    t_0_xxyz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(151));

    t_0_xxyz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(152));

    t_0_xxyz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(153));

    t_0_xxyz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(154));

    t_0_xxyz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(155));

    t_0_xxyz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(156));

    t_0_xxyz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(157));

    t_0_xxyz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(158));

    t_0_xxyz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(159));

    t_0_xxyz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(160));

    t_0_xxyz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(161));

    t_0_xxyz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(162));

    t_0_xxyz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(163));

    t_0_xxyz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(164));

    t_0_xxyy_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(165));

    t_0_xxyy_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(166));

    t_0_xxyy_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(167));

    t_0_xxyy_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(168));

    t_0_xxyy_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(169));

    t_0_xxyy_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(170));

    t_0_xxyy_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(171));

    t_0_xxyy_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(172));

    t_0_xxyy_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(173));

    t_0_xxyy_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(174));

    t_0_xxyy_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(175));

    t_0_xxyy_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(176));

    t_0_xxyy_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(177));

    t_0_xxyy_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(178));

    t_0_xxyy_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(179));

    t_0_xxxz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(180));

    t_0_xxxz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(181));

    t_0_xxxz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(182));

    t_0_xxxz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(183));

    t_0_xxxz_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(184));

    t_0_xxxz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(185));

    t_0_xxxz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(186));

    t_0_xxxz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(187));

    t_0_xxxz_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(188));

    t_0_xxxz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(189));

    t_0_xxxz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(190));

    t_0_xxxz_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(191));

    t_0_xxxz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(192));

    t_0_xxxz_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(193));

    t_0_xxxz_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(194));

    t_0_xxxy_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(195));

    t_0_xxxy_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(196));

    t_0_xxxy_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(197));

    t_0_xxxy_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(198));

    t_0_xxxy_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(199));

    t_0_xxxy_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(200));

    t_0_xxxy_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(201));

    t_0_xxxy_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(202));

    t_0_xxxy_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(203));

    t_0_xxxy_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(204));

    t_0_xxxy_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(205));

    t_0_xxxy_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(206));

    t_0_xxxy_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(207));

    t_0_xxxy_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(208));

    t_0_xxxy_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(209));

    t_0_xxxx_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(210));

    t_0_xxxx_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(211));

    t_0_xxxx_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(212));

    t_0_xxxx_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(213));

    t_0_xxxx_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(214));

    t_0_xxxx_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(215));

    t_0_xxxx_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(216));

    t_0_xxxx_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(217));

    t_0_xxxx_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(218));

    t_0_xxxx_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(219));

    t_0_xxxx_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(220));

    t_0_xxxx_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(221));

    t_0_xxxx_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(222));

    t_0_xxxx_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(223));

    t_0_xxxx_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(224));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzzz_0_xxz, t_0_yzzz_0_xxzz, t_0_yzzz_0_xyz,\
                           t_0_yzzz_0_xyzz, t_0_yzzz_0_xzz, t_0_yzzz_0_xzzz, t_0_yzzz_0_yyy,\
                           t_0_yzzz_0_yyyy, t_0_yzzz_0_yyyz, t_0_yzzz_0_yyz, t_0_yzzz_0_yyzz,\
                           t_0_yzzz_0_yzz, t_0_yzzz_0_yzzz, t_0_yzzz_0_zzz, t_0_yzzz_0_zzzz,\
                           t_0_yzzz_y_xzz, t_0_yzzz_y_yyy, t_0_yzzz_y_yyz, t_0_yzzz_y_yzz,\
                           t_0_yzzz_y_zzz, t_0_yzzz_z_xxz, t_0_yzzz_z_xyz, t_0_yzzz_z_xzz,\
                           t_0_yzzz_z_yyz, t_0_yzzz_z_yzz, t_0_yzzz_z_zzz, t_0_zzzz_0_xxx,\
                           t_0_zzzz_0_xxxx, t_0_zzzz_0_xxxy, t_0_zzzz_0_xxxz, t_0_zzzz_0_xxy,\
                           t_0_zzzz_0_xxyy, t_0_zzzz_0_xxyz, t_0_zzzz_0_xxz, t_0_zzzz_0_xxzz,\
                           t_0_zzzz_0_xyy, t_0_zzzz_0_xyyy, t_0_zzzz_0_xyyz, t_0_zzzz_0_xyz,\
                           t_0_zzzz_0_xyzz, t_0_zzzz_0_xzz, t_0_zzzz_0_xzzz, t_0_zzzz_0_yyy,\
                           t_0_zzzz_0_yyyy, t_0_zzzz_0_yyyz, t_0_zzzz_0_yyz, t_0_zzzz_0_yyzz,\
                           t_0_zzzz_0_yzz, t_0_zzzz_0_yzzz, t_0_zzzz_0_zzz, t_0_zzzz_0_zzzz,\
                           t_0_zzzz_x_xxx, t_0_zzzz_x_xxy, t_0_zzzz_x_xxz, t_0_zzzz_x_xyy,\
                           t_0_zzzz_x_xyz, t_0_zzzz_x_xzz, t_0_zzzz_x_yyy, t_0_zzzz_x_yyz,\
                           t_0_zzzz_x_yzz, t_0_zzzz_x_zzz, t_0_zzzz_y_xxy, t_0_zzzz_y_xxz,\
                           t_0_zzzz_y_xyy, t_0_zzzz_y_xyz, t_0_zzzz_y_xzz, t_0_zzzz_y_yyy,\
                           t_0_zzzz_y_yyz, t_0_zzzz_y_yzz, t_0_zzzz_y_zzz, t_0_zzzz_z_xxz,\
                           t_0_zzzz_z_xyz, t_0_zzzz_z_xzz, t_0_zzzz_z_yyz, t_0_zzzz_z_yzz,\
                           t_0_zzzz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_z_zzz[i] = t_0_zzzz_0_zzzz[i] - rcd_z[i] * t_0_zzzz_0_zzz[i];

        t_0_zzzz_z_yzz[i] = t_0_zzzz_0_yzzz[i] - rcd_z[i] * t_0_zzzz_0_yzz[i];

        t_0_zzzz_z_yyz[i] = t_0_zzzz_0_yyzz[i] - rcd_z[i] * t_0_zzzz_0_yyz[i];

        t_0_zzzz_z_xzz[i] = t_0_zzzz_0_xzzz[i] - rcd_z[i] * t_0_zzzz_0_xzz[i];

        t_0_zzzz_z_xyz[i] = t_0_zzzz_0_xyzz[i] - rcd_z[i] * t_0_zzzz_0_xyz[i];

        t_0_zzzz_z_xxz[i] = t_0_zzzz_0_xxzz[i] - rcd_z[i] * t_0_zzzz_0_xxz[i];

        t_0_zzzz_y_zzz[i] = t_0_zzzz_0_yzzz[i] - rcd_y[i] * t_0_zzzz_0_zzz[i];

        t_0_zzzz_y_yzz[i] = t_0_zzzz_0_yyzz[i] - rcd_y[i] * t_0_zzzz_0_yzz[i];

        t_0_zzzz_y_yyz[i] = t_0_zzzz_0_yyyz[i] - rcd_y[i] * t_0_zzzz_0_yyz[i];

        t_0_zzzz_y_yyy[i] = t_0_zzzz_0_yyyy[i] - rcd_y[i] * t_0_zzzz_0_yyy[i];

        t_0_zzzz_y_xzz[i] = t_0_zzzz_0_xyzz[i] - rcd_y[i] * t_0_zzzz_0_xzz[i];

        t_0_zzzz_y_xyz[i] = t_0_zzzz_0_xyyz[i] - rcd_y[i] * t_0_zzzz_0_xyz[i];

        t_0_zzzz_y_xyy[i] = t_0_zzzz_0_xyyy[i] - rcd_y[i] * t_0_zzzz_0_xyy[i];

        t_0_zzzz_y_xxz[i] = t_0_zzzz_0_xxyz[i] - rcd_y[i] * t_0_zzzz_0_xxz[i];

        t_0_zzzz_y_xxy[i] = t_0_zzzz_0_xxyy[i] - rcd_y[i] * t_0_zzzz_0_xxy[i];

        t_0_zzzz_x_zzz[i] = t_0_zzzz_0_xzzz[i] - rcd_x[i] * t_0_zzzz_0_zzz[i];

        t_0_zzzz_x_yzz[i] = t_0_zzzz_0_xyzz[i] - rcd_x[i] * t_0_zzzz_0_yzz[i];

        t_0_zzzz_x_yyz[i] = t_0_zzzz_0_xyyz[i] - rcd_x[i] * t_0_zzzz_0_yyz[i];

        t_0_zzzz_x_yyy[i] = t_0_zzzz_0_xyyy[i] - rcd_x[i] * t_0_zzzz_0_yyy[i];

        t_0_zzzz_x_xzz[i] = t_0_zzzz_0_xxzz[i] - rcd_x[i] * t_0_zzzz_0_xzz[i];

        t_0_zzzz_x_xyz[i] = t_0_zzzz_0_xxyz[i] - rcd_x[i] * t_0_zzzz_0_xyz[i];

        t_0_zzzz_x_xyy[i] = t_0_zzzz_0_xxyy[i] - rcd_x[i] * t_0_zzzz_0_xyy[i];

        t_0_zzzz_x_xxz[i] = t_0_zzzz_0_xxxz[i] - rcd_x[i] * t_0_zzzz_0_xxz[i];

        t_0_zzzz_x_xxy[i] = t_0_zzzz_0_xxxy[i] - rcd_x[i] * t_0_zzzz_0_xxy[i];

        t_0_zzzz_x_xxx[i] = t_0_zzzz_0_xxxx[i] - rcd_x[i] * t_0_zzzz_0_xxx[i];

        t_0_yzzz_z_zzz[i] = t_0_yzzz_0_zzzz[i] - rcd_z[i] * t_0_yzzz_0_zzz[i];

        t_0_yzzz_z_yzz[i] = t_0_yzzz_0_yzzz[i] - rcd_z[i] * t_0_yzzz_0_yzz[i];

        t_0_yzzz_z_yyz[i] = t_0_yzzz_0_yyzz[i] - rcd_z[i] * t_0_yzzz_0_yyz[i];

        t_0_yzzz_z_xzz[i] = t_0_yzzz_0_xzzz[i] - rcd_z[i] * t_0_yzzz_0_xzz[i];

        t_0_yzzz_z_xyz[i] = t_0_yzzz_0_xyzz[i] - rcd_z[i] * t_0_yzzz_0_xyz[i];

        t_0_yzzz_z_xxz[i] = t_0_yzzz_0_xxzz[i] - rcd_z[i] * t_0_yzzz_0_xxz[i];

        t_0_yzzz_y_zzz[i] = t_0_yzzz_0_yzzz[i] - rcd_y[i] * t_0_yzzz_0_zzz[i];

        t_0_yzzz_y_yzz[i] = t_0_yzzz_0_yyzz[i] - rcd_y[i] * t_0_yzzz_0_yzz[i];

        t_0_yzzz_y_yyz[i] = t_0_yzzz_0_yyyz[i] - rcd_y[i] * t_0_yzzz_0_yyz[i];

        t_0_yzzz_y_yyy[i] = t_0_yzzz_0_yyyy[i] - rcd_y[i] * t_0_yzzz_0_yyy[i];

        t_0_yzzz_y_xzz[i] = t_0_yzzz_0_xyzz[i] - rcd_y[i] * t_0_yzzz_0_xzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyzz_0_xxy, t_0_yyzz_0_xxyy, t_0_yyzz_0_xxyz,\
                           t_0_yyzz_0_xxz, t_0_yyzz_0_xxzz, t_0_yyzz_0_xyy, t_0_yyzz_0_xyyy,\
                           t_0_yyzz_0_xyyz, t_0_yyzz_0_xyz, t_0_yyzz_0_xyzz, t_0_yyzz_0_xzz,\
                           t_0_yyzz_0_xzzz, t_0_yyzz_0_yyy, t_0_yyzz_0_yyyy, t_0_yyzz_0_yyyz,\
                           t_0_yyzz_0_yyz, t_0_yyzz_0_yyzz, t_0_yyzz_0_yzz, t_0_yyzz_0_yzzz,\
                           t_0_yyzz_0_zzz, t_0_yyzz_0_zzzz, t_0_yyzz_x_xyy, t_0_yyzz_x_xyz,\
                           t_0_yyzz_x_xzz, t_0_yyzz_x_yyy, t_0_yyzz_x_yyz, t_0_yyzz_x_yzz,\
                           t_0_yyzz_x_zzz, t_0_yyzz_y_xxy, t_0_yyzz_y_xxz, t_0_yyzz_y_xyy,\
                           t_0_yyzz_y_xyz, t_0_yyzz_y_xzz, t_0_yyzz_y_yyy, t_0_yyzz_y_yyz,\
                           t_0_yyzz_y_yzz, t_0_yyzz_y_zzz, t_0_yyzz_z_xxz, t_0_yyzz_z_xyz,\
                           t_0_yyzz_z_xzz, t_0_yyzz_z_yyz, t_0_yyzz_z_yzz, t_0_yyzz_z_zzz,\
                           t_0_yzzz_0_xxx, t_0_yzzz_0_xxxx, t_0_yzzz_0_xxxy, t_0_yzzz_0_xxxz,\
                           t_0_yzzz_0_xxy, t_0_yzzz_0_xxyy, t_0_yzzz_0_xxyz, t_0_yzzz_0_xxz,\
                           t_0_yzzz_0_xxzz, t_0_yzzz_0_xyy, t_0_yzzz_0_xyyy, t_0_yzzz_0_xyyz,\
                           t_0_yzzz_0_xyz, t_0_yzzz_0_xyzz, t_0_yzzz_0_xzz, t_0_yzzz_0_xzzz,\
                           t_0_yzzz_0_yyy, t_0_yzzz_0_yyz, t_0_yzzz_0_yzz, t_0_yzzz_0_zzz,\
                           t_0_yzzz_x_xxx, t_0_yzzz_x_xxy, t_0_yzzz_x_xxz, t_0_yzzz_x_xyy,\
                           t_0_yzzz_x_xyz, t_0_yzzz_x_xzz, t_0_yzzz_x_yyy, t_0_yzzz_x_yyz,\
                           t_0_yzzz_x_yzz, t_0_yzzz_x_zzz, t_0_yzzz_y_xxy, t_0_yzzz_y_xxz,\
                           t_0_yzzz_y_xyy, t_0_yzzz_y_xyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yzzz_y_xyz[i] = t_0_yzzz_0_xyyz[i] - rcd_y[i] * t_0_yzzz_0_xyz[i];

        t_0_yzzz_y_xyy[i] = t_0_yzzz_0_xyyy[i] - rcd_y[i] * t_0_yzzz_0_xyy[i];

        t_0_yzzz_y_xxz[i] = t_0_yzzz_0_xxyz[i] - rcd_y[i] * t_0_yzzz_0_xxz[i];

        t_0_yzzz_y_xxy[i] = t_0_yzzz_0_xxyy[i] - rcd_y[i] * t_0_yzzz_0_xxy[i];

        t_0_yzzz_x_zzz[i] = t_0_yzzz_0_xzzz[i] - rcd_x[i] * t_0_yzzz_0_zzz[i];

        t_0_yzzz_x_yzz[i] = t_0_yzzz_0_xyzz[i] - rcd_x[i] * t_0_yzzz_0_yzz[i];

        t_0_yzzz_x_yyz[i] = t_0_yzzz_0_xyyz[i] - rcd_x[i] * t_0_yzzz_0_yyz[i];

        t_0_yzzz_x_yyy[i] = t_0_yzzz_0_xyyy[i] - rcd_x[i] * t_0_yzzz_0_yyy[i];

        t_0_yzzz_x_xzz[i] = t_0_yzzz_0_xxzz[i] - rcd_x[i] * t_0_yzzz_0_xzz[i];

        t_0_yzzz_x_xyz[i] = t_0_yzzz_0_xxyz[i] - rcd_x[i] * t_0_yzzz_0_xyz[i];

        t_0_yzzz_x_xyy[i] = t_0_yzzz_0_xxyy[i] - rcd_x[i] * t_0_yzzz_0_xyy[i];

        t_0_yzzz_x_xxz[i] = t_0_yzzz_0_xxxz[i] - rcd_x[i] * t_0_yzzz_0_xxz[i];

        t_0_yzzz_x_xxy[i] = t_0_yzzz_0_xxxy[i] - rcd_x[i] * t_0_yzzz_0_xxy[i];

        t_0_yzzz_x_xxx[i] = t_0_yzzz_0_xxxx[i] - rcd_x[i] * t_0_yzzz_0_xxx[i];

        t_0_yyzz_z_zzz[i] = t_0_yyzz_0_zzzz[i] - rcd_z[i] * t_0_yyzz_0_zzz[i];

        t_0_yyzz_z_yzz[i] = t_0_yyzz_0_yzzz[i] - rcd_z[i] * t_0_yyzz_0_yzz[i];

        t_0_yyzz_z_yyz[i] = t_0_yyzz_0_yyzz[i] - rcd_z[i] * t_0_yyzz_0_yyz[i];

        t_0_yyzz_z_xzz[i] = t_0_yyzz_0_xzzz[i] - rcd_z[i] * t_0_yyzz_0_xzz[i];

        t_0_yyzz_z_xyz[i] = t_0_yyzz_0_xyzz[i] - rcd_z[i] * t_0_yyzz_0_xyz[i];

        t_0_yyzz_z_xxz[i] = t_0_yyzz_0_xxzz[i] - rcd_z[i] * t_0_yyzz_0_xxz[i];

        t_0_yyzz_y_zzz[i] = t_0_yyzz_0_yzzz[i] - rcd_y[i] * t_0_yyzz_0_zzz[i];

        t_0_yyzz_y_yzz[i] = t_0_yyzz_0_yyzz[i] - rcd_y[i] * t_0_yyzz_0_yzz[i];

        t_0_yyzz_y_yyz[i] = t_0_yyzz_0_yyyz[i] - rcd_y[i] * t_0_yyzz_0_yyz[i];

        t_0_yyzz_y_yyy[i] = t_0_yyzz_0_yyyy[i] - rcd_y[i] * t_0_yyzz_0_yyy[i];

        t_0_yyzz_y_xzz[i] = t_0_yyzz_0_xyzz[i] - rcd_y[i] * t_0_yyzz_0_xzz[i];

        t_0_yyzz_y_xyz[i] = t_0_yyzz_0_xyyz[i] - rcd_y[i] * t_0_yyzz_0_xyz[i];

        t_0_yyzz_y_xyy[i] = t_0_yyzz_0_xyyy[i] - rcd_y[i] * t_0_yyzz_0_xyy[i];

        t_0_yyzz_y_xxz[i] = t_0_yyzz_0_xxyz[i] - rcd_y[i] * t_0_yyzz_0_xxz[i];

        t_0_yyzz_y_xxy[i] = t_0_yyzz_0_xxyy[i] - rcd_y[i] * t_0_yyzz_0_xxy[i];

        t_0_yyzz_x_zzz[i] = t_0_yyzz_0_xzzz[i] - rcd_x[i] * t_0_yyzz_0_zzz[i];

        t_0_yyzz_x_yzz[i] = t_0_yyzz_0_xyzz[i] - rcd_x[i] * t_0_yyzz_0_yzz[i];

        t_0_yyzz_x_yyz[i] = t_0_yyzz_0_xyyz[i] - rcd_x[i] * t_0_yyzz_0_yyz[i];

        t_0_yyzz_x_yyy[i] = t_0_yyzz_0_xyyy[i] - rcd_x[i] * t_0_yyzz_0_yyy[i];

        t_0_yyzz_x_xzz[i] = t_0_yyzz_0_xxzz[i] - rcd_x[i] * t_0_yyzz_0_xzz[i];

        t_0_yyzz_x_xyz[i] = t_0_yyzz_0_xxyz[i] - rcd_x[i] * t_0_yyzz_0_xyz[i];

        t_0_yyzz_x_xyy[i] = t_0_yyzz_0_xxyy[i] - rcd_x[i] * t_0_yyzz_0_xyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyyy_0_xxz, t_0_yyyy_0_xxzz, t_0_yyyy_0_xyz,\
                           t_0_yyyy_0_xyzz, t_0_yyyy_0_xzz, t_0_yyyy_0_xzzz, t_0_yyyy_0_yyz,\
                           t_0_yyyy_0_yyzz, t_0_yyyy_0_yzz, t_0_yyyy_0_yzzz, t_0_yyyy_0_zzz,\
                           t_0_yyyy_0_zzzz, t_0_yyyy_y_yzz, t_0_yyyy_y_zzz, t_0_yyyy_z_xxz,\
                           t_0_yyyy_z_xyz, t_0_yyyy_z_xzz, t_0_yyyy_z_yyz, t_0_yyyy_z_yzz,\
                           t_0_yyyy_z_zzz, t_0_yyyz_0_xxx, t_0_yyyz_0_xxxx, t_0_yyyz_0_xxxy,\
                           t_0_yyyz_0_xxxz, t_0_yyyz_0_xxy, t_0_yyyz_0_xxyy, t_0_yyyz_0_xxyz,\
                           t_0_yyyz_0_xxz, t_0_yyyz_0_xxzz, t_0_yyyz_0_xyy, t_0_yyyz_0_xyyy,\
                           t_0_yyyz_0_xyyz, t_0_yyyz_0_xyz, t_0_yyyz_0_xyzz, t_0_yyyz_0_xzz,\
                           t_0_yyyz_0_xzzz, t_0_yyyz_0_yyy, t_0_yyyz_0_yyyy, t_0_yyyz_0_yyyz,\
                           t_0_yyyz_0_yyz, t_0_yyyz_0_yyzz, t_0_yyyz_0_yzz, t_0_yyyz_0_yzzz,\
                           t_0_yyyz_0_zzz, t_0_yyyz_0_zzzz, t_0_yyyz_x_xxx, t_0_yyyz_x_xxy,\
                           t_0_yyyz_x_xxz, t_0_yyyz_x_xyy, t_0_yyyz_x_xyz, t_0_yyyz_x_xzz,\
                           t_0_yyyz_x_yyy, t_0_yyyz_x_yyz, t_0_yyyz_x_yzz, t_0_yyyz_x_zzz,\
                           t_0_yyyz_y_xxy, t_0_yyyz_y_xxz, t_0_yyyz_y_xyy, t_0_yyyz_y_xyz,\
                           t_0_yyyz_y_xzz, t_0_yyyz_y_yyy, t_0_yyyz_y_yyz, t_0_yyyz_y_yzz,\
                           t_0_yyyz_y_zzz, t_0_yyyz_z_xxz, t_0_yyyz_z_xyz, t_0_yyyz_z_xzz,\
                           t_0_yyyz_z_yyz, t_0_yyyz_z_yzz, t_0_yyyz_z_zzz, t_0_yyzz_0_xxx,\
                           t_0_yyzz_0_xxxx, t_0_yyzz_0_xxxy, t_0_yyzz_0_xxxz, t_0_yyzz_0_xxy,\
                           t_0_yyzz_0_xxz, t_0_yyzz_x_xxx, t_0_yyzz_x_xxy, t_0_yyzz_x_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyzz_x_xxz[i] = t_0_yyzz_0_xxxz[i] - rcd_x[i] * t_0_yyzz_0_xxz[i];

        t_0_yyzz_x_xxy[i] = t_0_yyzz_0_xxxy[i] - rcd_x[i] * t_0_yyzz_0_xxy[i];

        t_0_yyzz_x_xxx[i] = t_0_yyzz_0_xxxx[i] - rcd_x[i] * t_0_yyzz_0_xxx[i];

        t_0_yyyz_z_zzz[i] = t_0_yyyz_0_zzzz[i] - rcd_z[i] * t_0_yyyz_0_zzz[i];

        t_0_yyyz_z_yzz[i] = t_0_yyyz_0_yzzz[i] - rcd_z[i] * t_0_yyyz_0_yzz[i];

        t_0_yyyz_z_yyz[i] = t_0_yyyz_0_yyzz[i] - rcd_z[i] * t_0_yyyz_0_yyz[i];

        t_0_yyyz_z_xzz[i] = t_0_yyyz_0_xzzz[i] - rcd_z[i] * t_0_yyyz_0_xzz[i];

        t_0_yyyz_z_xyz[i] = t_0_yyyz_0_xyzz[i] - rcd_z[i] * t_0_yyyz_0_xyz[i];

        t_0_yyyz_z_xxz[i] = t_0_yyyz_0_xxzz[i] - rcd_z[i] * t_0_yyyz_0_xxz[i];

        t_0_yyyz_y_zzz[i] = t_0_yyyz_0_yzzz[i] - rcd_y[i] * t_0_yyyz_0_zzz[i];

        t_0_yyyz_y_yzz[i] = t_0_yyyz_0_yyzz[i] - rcd_y[i] * t_0_yyyz_0_yzz[i];

        t_0_yyyz_y_yyz[i] = t_0_yyyz_0_yyyz[i] - rcd_y[i] * t_0_yyyz_0_yyz[i];

        t_0_yyyz_y_yyy[i] = t_0_yyyz_0_yyyy[i] - rcd_y[i] * t_0_yyyz_0_yyy[i];

        t_0_yyyz_y_xzz[i] = t_0_yyyz_0_xyzz[i] - rcd_y[i] * t_0_yyyz_0_xzz[i];

        t_0_yyyz_y_xyz[i] = t_0_yyyz_0_xyyz[i] - rcd_y[i] * t_0_yyyz_0_xyz[i];

        t_0_yyyz_y_xyy[i] = t_0_yyyz_0_xyyy[i] - rcd_y[i] * t_0_yyyz_0_xyy[i];

        t_0_yyyz_y_xxz[i] = t_0_yyyz_0_xxyz[i] - rcd_y[i] * t_0_yyyz_0_xxz[i];

        t_0_yyyz_y_xxy[i] = t_0_yyyz_0_xxyy[i] - rcd_y[i] * t_0_yyyz_0_xxy[i];

        t_0_yyyz_x_zzz[i] = t_0_yyyz_0_xzzz[i] - rcd_x[i] * t_0_yyyz_0_zzz[i];

        t_0_yyyz_x_yzz[i] = t_0_yyyz_0_xyzz[i] - rcd_x[i] * t_0_yyyz_0_yzz[i];

        t_0_yyyz_x_yyz[i] = t_0_yyyz_0_xyyz[i] - rcd_x[i] * t_0_yyyz_0_yyz[i];

        t_0_yyyz_x_yyy[i] = t_0_yyyz_0_xyyy[i] - rcd_x[i] * t_0_yyyz_0_yyy[i];

        t_0_yyyz_x_xzz[i] = t_0_yyyz_0_xxzz[i] - rcd_x[i] * t_0_yyyz_0_xzz[i];

        t_0_yyyz_x_xyz[i] = t_0_yyyz_0_xxyz[i] - rcd_x[i] * t_0_yyyz_0_xyz[i];

        t_0_yyyz_x_xyy[i] = t_0_yyyz_0_xxyy[i] - rcd_x[i] * t_0_yyyz_0_xyy[i];

        t_0_yyyz_x_xxz[i] = t_0_yyyz_0_xxxz[i] - rcd_x[i] * t_0_yyyz_0_xxz[i];

        t_0_yyyz_x_xxy[i] = t_0_yyyz_0_xxxy[i] - rcd_x[i] * t_0_yyyz_0_xxy[i];

        t_0_yyyz_x_xxx[i] = t_0_yyyz_0_xxxx[i] - rcd_x[i] * t_0_yyyz_0_xxx[i];

        t_0_yyyy_z_zzz[i] = t_0_yyyy_0_zzzz[i] - rcd_z[i] * t_0_yyyy_0_zzz[i];

        t_0_yyyy_z_yzz[i] = t_0_yyyy_0_yzzz[i] - rcd_z[i] * t_0_yyyy_0_yzz[i];

        t_0_yyyy_z_yyz[i] = t_0_yyyy_0_yyzz[i] - rcd_z[i] * t_0_yyyy_0_yyz[i];

        t_0_yyyy_z_xzz[i] = t_0_yyyy_0_xzzz[i] - rcd_z[i] * t_0_yyyy_0_xzz[i];

        t_0_yyyy_z_xyz[i] = t_0_yyyy_0_xyzz[i] - rcd_z[i] * t_0_yyyy_0_xyz[i];

        t_0_yyyy_z_xxz[i] = t_0_yyyy_0_xxzz[i] - rcd_z[i] * t_0_yyyy_0_xxz[i];

        t_0_yyyy_y_zzz[i] = t_0_yyyy_0_yzzz[i] - rcd_y[i] * t_0_yyyy_0_zzz[i];

        t_0_yyyy_y_yzz[i] = t_0_yyyy_0_yyzz[i] - rcd_y[i] * t_0_yyyy_0_yzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xzzz_0_xxy, t_0_xzzz_0_xxyy, t_0_xzzz_0_xxyz,\
                           t_0_xzzz_0_xxz, t_0_xzzz_0_xxzz, t_0_xzzz_0_xyy, t_0_xzzz_0_xyyy,\
                           t_0_xzzz_0_xyyz, t_0_xzzz_0_xyz, t_0_xzzz_0_xyzz, t_0_xzzz_0_xzz,\
                           t_0_xzzz_0_xzzz, t_0_xzzz_0_yyy, t_0_xzzz_0_yyyy, t_0_xzzz_0_yyyz,\
                           t_0_xzzz_0_yyz, t_0_xzzz_0_yyzz, t_0_xzzz_0_yzz, t_0_xzzz_0_yzzz,\
                           t_0_xzzz_0_zzz, t_0_xzzz_0_zzzz, t_0_xzzz_x_yyy, t_0_xzzz_x_yyz,\
                           t_0_xzzz_x_yzz, t_0_xzzz_x_zzz, t_0_xzzz_y_xxy, t_0_xzzz_y_xxz,\
                           t_0_xzzz_y_xyy, t_0_xzzz_y_xyz, t_0_xzzz_y_xzz, t_0_xzzz_y_yyy,\
                           t_0_xzzz_y_yyz, t_0_xzzz_y_yzz, t_0_xzzz_y_zzz, t_0_xzzz_z_xxz,\
                           t_0_xzzz_z_xyz, t_0_xzzz_z_xzz, t_0_xzzz_z_yyz, t_0_xzzz_z_yzz,\
                           t_0_xzzz_z_zzz, t_0_yyyy_0_xxx, t_0_yyyy_0_xxxx, t_0_yyyy_0_xxxy,\
                           t_0_yyyy_0_xxxz, t_0_yyyy_0_xxy, t_0_yyyy_0_xxyy, t_0_yyyy_0_xxyz,\
                           t_0_yyyy_0_xxz, t_0_yyyy_0_xxzz, t_0_yyyy_0_xyy, t_0_yyyy_0_xyyy,\
                           t_0_yyyy_0_xyyz, t_0_yyyy_0_xyz, t_0_yyyy_0_xyzz, t_0_yyyy_0_xzz,\
                           t_0_yyyy_0_xzzz, t_0_yyyy_0_yyy, t_0_yyyy_0_yyyy, t_0_yyyy_0_yyyz,\
                           t_0_yyyy_0_yyz, t_0_yyyy_0_yzz, t_0_yyyy_0_zzz, t_0_yyyy_x_xxx,\
                           t_0_yyyy_x_xxy, t_0_yyyy_x_xxz, t_0_yyyy_x_xyy, t_0_yyyy_x_xyz,\
                           t_0_yyyy_x_xzz, t_0_yyyy_x_yyy, t_0_yyyy_x_yyz, t_0_yyyy_x_yzz,\
                           t_0_yyyy_x_zzz, t_0_yyyy_y_xxy, t_0_yyyy_y_xxz, t_0_yyyy_y_xyy,\
                           t_0_yyyy_y_xyz, t_0_yyyy_y_xzz, t_0_yyyy_y_yyy, t_0_yyyy_y_yyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyyy_y_yyz[i] = t_0_yyyy_0_yyyz[i] - rcd_y[i] * t_0_yyyy_0_yyz[i];

        t_0_yyyy_y_yyy[i] = t_0_yyyy_0_yyyy[i] - rcd_y[i] * t_0_yyyy_0_yyy[i];

        t_0_yyyy_y_xzz[i] = t_0_yyyy_0_xyzz[i] - rcd_y[i] * t_0_yyyy_0_xzz[i];

        t_0_yyyy_y_xyz[i] = t_0_yyyy_0_xyyz[i] - rcd_y[i] * t_0_yyyy_0_xyz[i];

        t_0_yyyy_y_xyy[i] = t_0_yyyy_0_xyyy[i] - rcd_y[i] * t_0_yyyy_0_xyy[i];

        t_0_yyyy_y_xxz[i] = t_0_yyyy_0_xxyz[i] - rcd_y[i] * t_0_yyyy_0_xxz[i];

        t_0_yyyy_y_xxy[i] = t_0_yyyy_0_xxyy[i] - rcd_y[i] * t_0_yyyy_0_xxy[i];

        t_0_yyyy_x_zzz[i] = t_0_yyyy_0_xzzz[i] - rcd_x[i] * t_0_yyyy_0_zzz[i];

        t_0_yyyy_x_yzz[i] = t_0_yyyy_0_xyzz[i] - rcd_x[i] * t_0_yyyy_0_yzz[i];

        t_0_yyyy_x_yyz[i] = t_0_yyyy_0_xyyz[i] - rcd_x[i] * t_0_yyyy_0_yyz[i];

        t_0_yyyy_x_yyy[i] = t_0_yyyy_0_xyyy[i] - rcd_x[i] * t_0_yyyy_0_yyy[i];

        t_0_yyyy_x_xzz[i] = t_0_yyyy_0_xxzz[i] - rcd_x[i] * t_0_yyyy_0_xzz[i];

        t_0_yyyy_x_xyz[i] = t_0_yyyy_0_xxyz[i] - rcd_x[i] * t_0_yyyy_0_xyz[i];

        t_0_yyyy_x_xyy[i] = t_0_yyyy_0_xxyy[i] - rcd_x[i] * t_0_yyyy_0_xyy[i];

        t_0_yyyy_x_xxz[i] = t_0_yyyy_0_xxxz[i] - rcd_x[i] * t_0_yyyy_0_xxz[i];

        t_0_yyyy_x_xxy[i] = t_0_yyyy_0_xxxy[i] - rcd_x[i] * t_0_yyyy_0_xxy[i];

        t_0_yyyy_x_xxx[i] = t_0_yyyy_0_xxxx[i] - rcd_x[i] * t_0_yyyy_0_xxx[i];

        t_0_xzzz_z_zzz[i] = t_0_xzzz_0_zzzz[i] - rcd_z[i] * t_0_xzzz_0_zzz[i];

        t_0_xzzz_z_yzz[i] = t_0_xzzz_0_yzzz[i] - rcd_z[i] * t_0_xzzz_0_yzz[i];

        t_0_xzzz_z_yyz[i] = t_0_xzzz_0_yyzz[i] - rcd_z[i] * t_0_xzzz_0_yyz[i];

        t_0_xzzz_z_xzz[i] = t_0_xzzz_0_xzzz[i] - rcd_z[i] * t_0_xzzz_0_xzz[i];

        t_0_xzzz_z_xyz[i] = t_0_xzzz_0_xyzz[i] - rcd_z[i] * t_0_xzzz_0_xyz[i];

        t_0_xzzz_z_xxz[i] = t_0_xzzz_0_xxzz[i] - rcd_z[i] * t_0_xzzz_0_xxz[i];

        t_0_xzzz_y_zzz[i] = t_0_xzzz_0_yzzz[i] - rcd_y[i] * t_0_xzzz_0_zzz[i];

        t_0_xzzz_y_yzz[i] = t_0_xzzz_0_yyzz[i] - rcd_y[i] * t_0_xzzz_0_yzz[i];

        t_0_xzzz_y_yyz[i] = t_0_xzzz_0_yyyz[i] - rcd_y[i] * t_0_xzzz_0_yyz[i];

        t_0_xzzz_y_yyy[i] = t_0_xzzz_0_yyyy[i] - rcd_y[i] * t_0_xzzz_0_yyy[i];

        t_0_xzzz_y_xzz[i] = t_0_xzzz_0_xyzz[i] - rcd_y[i] * t_0_xzzz_0_xzz[i];

        t_0_xzzz_y_xyz[i] = t_0_xzzz_0_xyyz[i] - rcd_y[i] * t_0_xzzz_0_xyz[i];

        t_0_xzzz_y_xyy[i] = t_0_xzzz_0_xyyy[i] - rcd_y[i] * t_0_xzzz_0_xyy[i];

        t_0_xzzz_y_xxz[i] = t_0_xzzz_0_xxyz[i] - rcd_y[i] * t_0_xzzz_0_xxz[i];

        t_0_xzzz_y_xxy[i] = t_0_xzzz_0_xxyy[i] - rcd_y[i] * t_0_xzzz_0_xxy[i];

        t_0_xzzz_x_zzz[i] = t_0_xzzz_0_xzzz[i] - rcd_x[i] * t_0_xzzz_0_zzz[i];

        t_0_xzzz_x_yzz[i] = t_0_xzzz_0_xyzz[i] - rcd_x[i] * t_0_xzzz_0_yzz[i];

        t_0_xzzz_x_yyz[i] = t_0_xzzz_0_xyyz[i] - rcd_x[i] * t_0_xzzz_0_yyz[i];

        t_0_xzzz_x_yyy[i] = t_0_xzzz_0_xyyy[i] - rcd_x[i] * t_0_xzzz_0_yyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyz_0_xyz, t_0_xyyz_0_xyzz, t_0_xyyz_0_xzz,\
                           t_0_xyyz_0_xzzz, t_0_xyyz_0_yyz, t_0_xyyz_0_yyzz, t_0_xyyz_0_yzz,\
                           t_0_xyyz_0_yzzz, t_0_xyyz_0_zzz, t_0_xyyz_0_zzzz, t_0_xyyz_z_xyz,\
                           t_0_xyyz_z_xzz, t_0_xyyz_z_yyz, t_0_xyyz_z_yzz, t_0_xyyz_z_zzz,\
                           t_0_xyzz_0_xxx, t_0_xyzz_0_xxxx, t_0_xyzz_0_xxxy, t_0_xyzz_0_xxxz,\
                           t_0_xyzz_0_xxy, t_0_xyzz_0_xxyy, t_0_xyzz_0_xxyz, t_0_xyzz_0_xxz,\
                           t_0_xyzz_0_xxzz, t_0_xyzz_0_xyy, t_0_xyzz_0_xyyy, t_0_xyzz_0_xyyz,\
                           t_0_xyzz_0_xyz, t_0_xyzz_0_xyzz, t_0_xyzz_0_xzz, t_0_xyzz_0_xzzz,\
                           t_0_xyzz_0_yyy, t_0_xyzz_0_yyyy, t_0_xyzz_0_yyyz, t_0_xyzz_0_yyz,\
                           t_0_xyzz_0_yyzz, t_0_xyzz_0_yzz, t_0_xyzz_0_yzzz, t_0_xyzz_0_zzz,\
                           t_0_xyzz_0_zzzz, t_0_xyzz_x_xxx, t_0_xyzz_x_xxy, t_0_xyzz_x_xxz,\
                           t_0_xyzz_x_xyy, t_0_xyzz_x_xyz, t_0_xyzz_x_xzz, t_0_xyzz_x_yyy,\
                           t_0_xyzz_x_yyz, t_0_xyzz_x_yzz, t_0_xyzz_x_zzz, t_0_xyzz_y_xxy,\
                           t_0_xyzz_y_xxz, t_0_xyzz_y_xyy, t_0_xyzz_y_xyz, t_0_xyzz_y_xzz,\
                           t_0_xyzz_y_yyy, t_0_xyzz_y_yyz, t_0_xyzz_y_yzz, t_0_xyzz_y_zzz,\
                           t_0_xyzz_z_xxz, t_0_xyzz_z_xyz, t_0_xyzz_z_xzz, t_0_xyzz_z_yyz,\
                           t_0_xyzz_z_yzz, t_0_xyzz_z_zzz, t_0_xzzz_0_xxx, t_0_xzzz_0_xxxx,\
                           t_0_xzzz_0_xxxy, t_0_xzzz_0_xxxz, t_0_xzzz_0_xxy, t_0_xzzz_0_xxyy,\
                           t_0_xzzz_0_xxyz, t_0_xzzz_0_xxz, t_0_xzzz_0_xxzz, t_0_xzzz_0_xyy,\
                           t_0_xzzz_0_xyz, t_0_xzzz_0_xzz, t_0_xzzz_x_xxx, t_0_xzzz_x_xxy,\
                           t_0_xzzz_x_xxz, t_0_xzzz_x_xyy, t_0_xzzz_x_xyz, t_0_xzzz_x_xzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzzz_x_xzz[i] = t_0_xzzz_0_xxzz[i] - rcd_x[i] * t_0_xzzz_0_xzz[i];

        t_0_xzzz_x_xyz[i] = t_0_xzzz_0_xxyz[i] - rcd_x[i] * t_0_xzzz_0_xyz[i];

        t_0_xzzz_x_xyy[i] = t_0_xzzz_0_xxyy[i] - rcd_x[i] * t_0_xzzz_0_xyy[i];

        t_0_xzzz_x_xxz[i] = t_0_xzzz_0_xxxz[i] - rcd_x[i] * t_0_xzzz_0_xxz[i];

        t_0_xzzz_x_xxy[i] = t_0_xzzz_0_xxxy[i] - rcd_x[i] * t_0_xzzz_0_xxy[i];

        t_0_xzzz_x_xxx[i] = t_0_xzzz_0_xxxx[i] - rcd_x[i] * t_0_xzzz_0_xxx[i];

        t_0_xyzz_z_zzz[i] = t_0_xyzz_0_zzzz[i] - rcd_z[i] * t_0_xyzz_0_zzz[i];

        t_0_xyzz_z_yzz[i] = t_0_xyzz_0_yzzz[i] - rcd_z[i] * t_0_xyzz_0_yzz[i];

        t_0_xyzz_z_yyz[i] = t_0_xyzz_0_yyzz[i] - rcd_z[i] * t_0_xyzz_0_yyz[i];

        t_0_xyzz_z_xzz[i] = t_0_xyzz_0_xzzz[i] - rcd_z[i] * t_0_xyzz_0_xzz[i];

        t_0_xyzz_z_xyz[i] = t_0_xyzz_0_xyzz[i] - rcd_z[i] * t_0_xyzz_0_xyz[i];

        t_0_xyzz_z_xxz[i] = t_0_xyzz_0_xxzz[i] - rcd_z[i] * t_0_xyzz_0_xxz[i];

        t_0_xyzz_y_zzz[i] = t_0_xyzz_0_yzzz[i] - rcd_y[i] * t_0_xyzz_0_zzz[i];

        t_0_xyzz_y_yzz[i] = t_0_xyzz_0_yyzz[i] - rcd_y[i] * t_0_xyzz_0_yzz[i];

        t_0_xyzz_y_yyz[i] = t_0_xyzz_0_yyyz[i] - rcd_y[i] * t_0_xyzz_0_yyz[i];

        t_0_xyzz_y_yyy[i] = t_0_xyzz_0_yyyy[i] - rcd_y[i] * t_0_xyzz_0_yyy[i];

        t_0_xyzz_y_xzz[i] = t_0_xyzz_0_xyzz[i] - rcd_y[i] * t_0_xyzz_0_xzz[i];

        t_0_xyzz_y_xyz[i] = t_0_xyzz_0_xyyz[i] - rcd_y[i] * t_0_xyzz_0_xyz[i];

        t_0_xyzz_y_xyy[i] = t_0_xyzz_0_xyyy[i] - rcd_y[i] * t_0_xyzz_0_xyy[i];

        t_0_xyzz_y_xxz[i] = t_0_xyzz_0_xxyz[i] - rcd_y[i] * t_0_xyzz_0_xxz[i];

        t_0_xyzz_y_xxy[i] = t_0_xyzz_0_xxyy[i] - rcd_y[i] * t_0_xyzz_0_xxy[i];

        t_0_xyzz_x_zzz[i] = t_0_xyzz_0_xzzz[i] - rcd_x[i] * t_0_xyzz_0_zzz[i];

        t_0_xyzz_x_yzz[i] = t_0_xyzz_0_xyzz[i] - rcd_x[i] * t_0_xyzz_0_yzz[i];

        t_0_xyzz_x_yyz[i] = t_0_xyzz_0_xyyz[i] - rcd_x[i] * t_0_xyzz_0_yyz[i];

        t_0_xyzz_x_yyy[i] = t_0_xyzz_0_xyyy[i] - rcd_x[i] * t_0_xyzz_0_yyy[i];

        t_0_xyzz_x_xzz[i] = t_0_xyzz_0_xxzz[i] - rcd_x[i] * t_0_xyzz_0_xzz[i];

        t_0_xyzz_x_xyz[i] = t_0_xyzz_0_xxyz[i] - rcd_x[i] * t_0_xyzz_0_xyz[i];

        t_0_xyzz_x_xyy[i] = t_0_xyzz_0_xxyy[i] - rcd_x[i] * t_0_xyzz_0_xyy[i];

        t_0_xyzz_x_xxz[i] = t_0_xyzz_0_xxxz[i] - rcd_x[i] * t_0_xyzz_0_xxz[i];

        t_0_xyzz_x_xxy[i] = t_0_xyzz_0_xxxy[i] - rcd_x[i] * t_0_xyzz_0_xxy[i];

        t_0_xyzz_x_xxx[i] = t_0_xyzz_0_xxxx[i] - rcd_x[i] * t_0_xyzz_0_xxx[i];

        t_0_xyyz_z_zzz[i] = t_0_xyyz_0_zzzz[i] - rcd_z[i] * t_0_xyyz_0_zzz[i];

        t_0_xyyz_z_yzz[i] = t_0_xyyz_0_yzzz[i] - rcd_z[i] * t_0_xyyz_0_yzz[i];

        t_0_xyyz_z_yyz[i] = t_0_xyyz_0_yyzz[i] - rcd_z[i] * t_0_xyyz_0_yyz[i];

        t_0_xyyz_z_xzz[i] = t_0_xyyz_0_xzzz[i] - rcd_z[i] * t_0_xyyz_0_xzz[i];

        t_0_xyyz_z_xyz[i] = t_0_xyyz_0_xyzz[i] - rcd_z[i] * t_0_xyyz_0_xyz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyy_0_xxy, t_0_xyyy_0_xxyy, t_0_xyyy_0_xxyz,\
                           t_0_xyyy_0_xxz, t_0_xyyy_0_xxzz, t_0_xyyy_0_xyy, t_0_xyyy_0_xyyy,\
                           t_0_xyyy_0_xyyz, t_0_xyyy_0_xyz, t_0_xyyy_0_xyzz, t_0_xyyy_0_xzz,\
                           t_0_xyyy_0_xzzz, t_0_xyyy_0_yyy, t_0_xyyy_0_yyyy, t_0_xyyy_0_yyyz,\
                           t_0_xyyy_0_yyz, t_0_xyyy_0_yyzz, t_0_xyyy_0_yzz, t_0_xyyy_0_yzzz,\
                           t_0_xyyy_0_zzz, t_0_xyyy_0_zzzz, t_0_xyyy_x_zzz, t_0_xyyy_y_xxy,\
                           t_0_xyyy_y_xxz, t_0_xyyy_y_xyy, t_0_xyyy_y_xyz, t_0_xyyy_y_xzz,\
                           t_0_xyyy_y_yyy, t_0_xyyy_y_yyz, t_0_xyyy_y_yzz, t_0_xyyy_y_zzz,\
                           t_0_xyyy_z_xxz, t_0_xyyy_z_xyz, t_0_xyyy_z_xzz, t_0_xyyy_z_yyz,\
                           t_0_xyyy_z_yzz, t_0_xyyy_z_zzz, t_0_xyyz_0_xxx, t_0_xyyz_0_xxxx,\
                           t_0_xyyz_0_xxxy, t_0_xyyz_0_xxxz, t_0_xyyz_0_xxy, t_0_xyyz_0_xxyy,\
                           t_0_xyyz_0_xxyz, t_0_xyyz_0_xxz, t_0_xyyz_0_xxzz, t_0_xyyz_0_xyy,\
                           t_0_xyyz_0_xyyy, t_0_xyyz_0_xyyz, t_0_xyyz_0_xyz, t_0_xyyz_0_xyzz,\
                           t_0_xyyz_0_xzz, t_0_xyyz_0_xzzz, t_0_xyyz_0_yyy, t_0_xyyz_0_yyyy,\
                           t_0_xyyz_0_yyyz, t_0_xyyz_0_yyz, t_0_xyyz_0_yyzz, t_0_xyyz_0_yzz,\
                           t_0_xyyz_0_yzzz, t_0_xyyz_0_zzz, t_0_xyyz_x_xxx, t_0_xyyz_x_xxy,\
                           t_0_xyyz_x_xxz, t_0_xyyz_x_xyy, t_0_xyyz_x_xyz, t_0_xyyz_x_xzz,\
                           t_0_xyyz_x_yyy, t_0_xyyz_x_yyz, t_0_xyyz_x_yzz, t_0_xyyz_x_zzz,\
                           t_0_xyyz_y_xxy, t_0_xyyz_y_xxz, t_0_xyyz_y_xyy, t_0_xyyz_y_xyz,\
                           t_0_xyyz_y_xzz, t_0_xyyz_y_yyy, t_0_xyyz_y_yyz, t_0_xyyz_y_yzz,\
                           t_0_xyyz_y_zzz, t_0_xyyz_z_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyz_z_xxz[i] = t_0_xyyz_0_xxzz[i] - rcd_z[i] * t_0_xyyz_0_xxz[i];

        t_0_xyyz_y_zzz[i] = t_0_xyyz_0_yzzz[i] - rcd_y[i] * t_0_xyyz_0_zzz[i];

        t_0_xyyz_y_yzz[i] = t_0_xyyz_0_yyzz[i] - rcd_y[i] * t_0_xyyz_0_yzz[i];

        t_0_xyyz_y_yyz[i] = t_0_xyyz_0_yyyz[i] - rcd_y[i] * t_0_xyyz_0_yyz[i];

        t_0_xyyz_y_yyy[i] = t_0_xyyz_0_yyyy[i] - rcd_y[i] * t_0_xyyz_0_yyy[i];

        t_0_xyyz_y_xzz[i] = t_0_xyyz_0_xyzz[i] - rcd_y[i] * t_0_xyyz_0_xzz[i];

        t_0_xyyz_y_xyz[i] = t_0_xyyz_0_xyyz[i] - rcd_y[i] * t_0_xyyz_0_xyz[i];

        t_0_xyyz_y_xyy[i] = t_0_xyyz_0_xyyy[i] - rcd_y[i] * t_0_xyyz_0_xyy[i];

        t_0_xyyz_y_xxz[i] = t_0_xyyz_0_xxyz[i] - rcd_y[i] * t_0_xyyz_0_xxz[i];

        t_0_xyyz_y_xxy[i] = t_0_xyyz_0_xxyy[i] - rcd_y[i] * t_0_xyyz_0_xxy[i];

        t_0_xyyz_x_zzz[i] = t_0_xyyz_0_xzzz[i] - rcd_x[i] * t_0_xyyz_0_zzz[i];

        t_0_xyyz_x_yzz[i] = t_0_xyyz_0_xyzz[i] - rcd_x[i] * t_0_xyyz_0_yzz[i];

        t_0_xyyz_x_yyz[i] = t_0_xyyz_0_xyyz[i] - rcd_x[i] * t_0_xyyz_0_yyz[i];

        t_0_xyyz_x_yyy[i] = t_0_xyyz_0_xyyy[i] - rcd_x[i] * t_0_xyyz_0_yyy[i];

        t_0_xyyz_x_xzz[i] = t_0_xyyz_0_xxzz[i] - rcd_x[i] * t_0_xyyz_0_xzz[i];

        t_0_xyyz_x_xyz[i] = t_0_xyyz_0_xxyz[i] - rcd_x[i] * t_0_xyyz_0_xyz[i];

        t_0_xyyz_x_xyy[i] = t_0_xyyz_0_xxyy[i] - rcd_x[i] * t_0_xyyz_0_xyy[i];

        t_0_xyyz_x_xxz[i] = t_0_xyyz_0_xxxz[i] - rcd_x[i] * t_0_xyyz_0_xxz[i];

        t_0_xyyz_x_xxy[i] = t_0_xyyz_0_xxxy[i] - rcd_x[i] * t_0_xyyz_0_xxy[i];

        t_0_xyyz_x_xxx[i] = t_0_xyyz_0_xxxx[i] - rcd_x[i] * t_0_xyyz_0_xxx[i];

        t_0_xyyy_z_zzz[i] = t_0_xyyy_0_zzzz[i] - rcd_z[i] * t_0_xyyy_0_zzz[i];

        t_0_xyyy_z_yzz[i] = t_0_xyyy_0_yzzz[i] - rcd_z[i] * t_0_xyyy_0_yzz[i];

        t_0_xyyy_z_yyz[i] = t_0_xyyy_0_yyzz[i] - rcd_z[i] * t_0_xyyy_0_yyz[i];

        t_0_xyyy_z_xzz[i] = t_0_xyyy_0_xzzz[i] - rcd_z[i] * t_0_xyyy_0_xzz[i];

        t_0_xyyy_z_xyz[i] = t_0_xyyy_0_xyzz[i] - rcd_z[i] * t_0_xyyy_0_xyz[i];

        t_0_xyyy_z_xxz[i] = t_0_xyyy_0_xxzz[i] - rcd_z[i] * t_0_xyyy_0_xxz[i];

        t_0_xyyy_y_zzz[i] = t_0_xyyy_0_yzzz[i] - rcd_y[i] * t_0_xyyy_0_zzz[i];

        t_0_xyyy_y_yzz[i] = t_0_xyyy_0_yyzz[i] - rcd_y[i] * t_0_xyyy_0_yzz[i];

        t_0_xyyy_y_yyz[i] = t_0_xyyy_0_yyyz[i] - rcd_y[i] * t_0_xyyy_0_yyz[i];

        t_0_xyyy_y_yyy[i] = t_0_xyyy_0_yyyy[i] - rcd_y[i] * t_0_xyyy_0_yyy[i];

        t_0_xyyy_y_xzz[i] = t_0_xyyy_0_xyzz[i] - rcd_y[i] * t_0_xyyy_0_xzz[i];

        t_0_xyyy_y_xyz[i] = t_0_xyyy_0_xyyz[i] - rcd_y[i] * t_0_xyyy_0_xyz[i];

        t_0_xyyy_y_xyy[i] = t_0_xyyy_0_xyyy[i] - rcd_y[i] * t_0_xyyy_0_xyy[i];

        t_0_xyyy_y_xxz[i] = t_0_xyyy_0_xxyz[i] - rcd_y[i] * t_0_xyyy_0_xxz[i];

        t_0_xyyy_y_xxy[i] = t_0_xyyy_0_xxyy[i] - rcd_y[i] * t_0_xyyy_0_xxy[i];

        t_0_xyyy_x_zzz[i] = t_0_xyyy_0_xzzz[i] - rcd_x[i] * t_0_xyyy_0_zzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyz_0_yzz, t_0_xxyz_0_yzzz, t_0_xxyz_0_zzz,\
                           t_0_xxyz_0_zzzz, t_0_xxyz_z_yzz, t_0_xxyz_z_zzz, t_0_xxzz_0_xxx,\
                           t_0_xxzz_0_xxxx, t_0_xxzz_0_xxxy, t_0_xxzz_0_xxxz, t_0_xxzz_0_xxy,\
                           t_0_xxzz_0_xxyy, t_0_xxzz_0_xxyz, t_0_xxzz_0_xxz, t_0_xxzz_0_xxzz,\
                           t_0_xxzz_0_xyy, t_0_xxzz_0_xyyy, t_0_xxzz_0_xyyz, t_0_xxzz_0_xyz,\
                           t_0_xxzz_0_xyzz, t_0_xxzz_0_xzz, t_0_xxzz_0_xzzz, t_0_xxzz_0_yyy,\
                           t_0_xxzz_0_yyyy, t_0_xxzz_0_yyyz, t_0_xxzz_0_yyz, t_0_xxzz_0_yyzz,\
                           t_0_xxzz_0_yzz, t_0_xxzz_0_yzzz, t_0_xxzz_0_zzz, t_0_xxzz_0_zzzz,\
                           t_0_xxzz_x_xxx, t_0_xxzz_x_xxy, t_0_xxzz_x_xxz, t_0_xxzz_x_xyy,\
                           t_0_xxzz_x_xyz, t_0_xxzz_x_xzz, t_0_xxzz_x_yyy, t_0_xxzz_x_yyz,\
                           t_0_xxzz_x_yzz, t_0_xxzz_x_zzz, t_0_xxzz_y_xxy, t_0_xxzz_y_xxz,\
                           t_0_xxzz_y_xyy, t_0_xxzz_y_xyz, t_0_xxzz_y_xzz, t_0_xxzz_y_yyy,\
                           t_0_xxzz_y_yyz, t_0_xxzz_y_yzz, t_0_xxzz_y_zzz, t_0_xxzz_z_xxz,\
                           t_0_xxzz_z_xyz, t_0_xxzz_z_xzz, t_0_xxzz_z_yyz, t_0_xxzz_z_yzz,\
                           t_0_xxzz_z_zzz, t_0_xyyy_0_xxx, t_0_xyyy_0_xxxx, t_0_xyyy_0_xxxy,\
                           t_0_xyyy_0_xxxz, t_0_xyyy_0_xxy, t_0_xyyy_0_xxyy, t_0_xyyy_0_xxyz,\
                           t_0_xyyy_0_xxz, t_0_xyyy_0_xxzz, t_0_xyyy_0_xyy, t_0_xyyy_0_xyyy,\
                           t_0_xyyy_0_xyyz, t_0_xyyy_0_xyz, t_0_xyyy_0_xyzz, t_0_xyyy_0_xzz,\
                           t_0_xyyy_0_yyy, t_0_xyyy_0_yyz, t_0_xyyy_0_yzz, t_0_xyyy_x_xxx,\
                           t_0_xyyy_x_xxy, t_0_xyyy_x_xxz, t_0_xyyy_x_xyy, t_0_xyyy_x_xyz,\
                           t_0_xyyy_x_xzz, t_0_xyyy_x_yyy, t_0_xyyy_x_yyz, t_0_xyyy_x_yzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyy_x_yzz[i] = t_0_xyyy_0_xyzz[i] - rcd_x[i] * t_0_xyyy_0_yzz[i];

        t_0_xyyy_x_yyz[i] = t_0_xyyy_0_xyyz[i] - rcd_x[i] * t_0_xyyy_0_yyz[i];

        t_0_xyyy_x_yyy[i] = t_0_xyyy_0_xyyy[i] - rcd_x[i] * t_0_xyyy_0_yyy[i];

        t_0_xyyy_x_xzz[i] = t_0_xyyy_0_xxzz[i] - rcd_x[i] * t_0_xyyy_0_xzz[i];

        t_0_xyyy_x_xyz[i] = t_0_xyyy_0_xxyz[i] - rcd_x[i] * t_0_xyyy_0_xyz[i];

        t_0_xyyy_x_xyy[i] = t_0_xyyy_0_xxyy[i] - rcd_x[i] * t_0_xyyy_0_xyy[i];

        t_0_xyyy_x_xxz[i] = t_0_xyyy_0_xxxz[i] - rcd_x[i] * t_0_xyyy_0_xxz[i];

        t_0_xyyy_x_xxy[i] = t_0_xyyy_0_xxxy[i] - rcd_x[i] * t_0_xyyy_0_xxy[i];

        t_0_xyyy_x_xxx[i] = t_0_xyyy_0_xxxx[i] - rcd_x[i] * t_0_xyyy_0_xxx[i];

        t_0_xxzz_z_zzz[i] = t_0_xxzz_0_zzzz[i] - rcd_z[i] * t_0_xxzz_0_zzz[i];

        t_0_xxzz_z_yzz[i] = t_0_xxzz_0_yzzz[i] - rcd_z[i] * t_0_xxzz_0_yzz[i];

        t_0_xxzz_z_yyz[i] = t_0_xxzz_0_yyzz[i] - rcd_z[i] * t_0_xxzz_0_yyz[i];

        t_0_xxzz_z_xzz[i] = t_0_xxzz_0_xzzz[i] - rcd_z[i] * t_0_xxzz_0_xzz[i];

        t_0_xxzz_z_xyz[i] = t_0_xxzz_0_xyzz[i] - rcd_z[i] * t_0_xxzz_0_xyz[i];

        t_0_xxzz_z_xxz[i] = t_0_xxzz_0_xxzz[i] - rcd_z[i] * t_0_xxzz_0_xxz[i];

        t_0_xxzz_y_zzz[i] = t_0_xxzz_0_yzzz[i] - rcd_y[i] * t_0_xxzz_0_zzz[i];

        t_0_xxzz_y_yzz[i] = t_0_xxzz_0_yyzz[i] - rcd_y[i] * t_0_xxzz_0_yzz[i];

        t_0_xxzz_y_yyz[i] = t_0_xxzz_0_yyyz[i] - rcd_y[i] * t_0_xxzz_0_yyz[i];

        t_0_xxzz_y_yyy[i] = t_0_xxzz_0_yyyy[i] - rcd_y[i] * t_0_xxzz_0_yyy[i];

        t_0_xxzz_y_xzz[i] = t_0_xxzz_0_xyzz[i] - rcd_y[i] * t_0_xxzz_0_xzz[i];

        t_0_xxzz_y_xyz[i] = t_0_xxzz_0_xyyz[i] - rcd_y[i] * t_0_xxzz_0_xyz[i];

        t_0_xxzz_y_xyy[i] = t_0_xxzz_0_xyyy[i] - rcd_y[i] * t_0_xxzz_0_xyy[i];

        t_0_xxzz_y_xxz[i] = t_0_xxzz_0_xxyz[i] - rcd_y[i] * t_0_xxzz_0_xxz[i];

        t_0_xxzz_y_xxy[i] = t_0_xxzz_0_xxyy[i] - rcd_y[i] * t_0_xxzz_0_xxy[i];

        t_0_xxzz_x_zzz[i] = t_0_xxzz_0_xzzz[i] - rcd_x[i] * t_0_xxzz_0_zzz[i];

        t_0_xxzz_x_yzz[i] = t_0_xxzz_0_xyzz[i] - rcd_x[i] * t_0_xxzz_0_yzz[i];

        t_0_xxzz_x_yyz[i] = t_0_xxzz_0_xyyz[i] - rcd_x[i] * t_0_xxzz_0_yyz[i];

        t_0_xxzz_x_yyy[i] = t_0_xxzz_0_xyyy[i] - rcd_x[i] * t_0_xxzz_0_yyy[i];

        t_0_xxzz_x_xzz[i] = t_0_xxzz_0_xxzz[i] - rcd_x[i] * t_0_xxzz_0_xzz[i];

        t_0_xxzz_x_xyz[i] = t_0_xxzz_0_xxyz[i] - rcd_x[i] * t_0_xxzz_0_xyz[i];

        t_0_xxzz_x_xyy[i] = t_0_xxzz_0_xxyy[i] - rcd_x[i] * t_0_xxzz_0_xyy[i];

        t_0_xxzz_x_xxz[i] = t_0_xxzz_0_xxxz[i] - rcd_x[i] * t_0_xxzz_0_xxz[i];

        t_0_xxzz_x_xxy[i] = t_0_xxzz_0_xxxy[i] - rcd_x[i] * t_0_xxzz_0_xxy[i];

        t_0_xxzz_x_xxx[i] = t_0_xxzz_0_xxxx[i] - rcd_x[i] * t_0_xxzz_0_xxx[i];

        t_0_xxyz_z_zzz[i] = t_0_xxyz_0_zzzz[i] - rcd_z[i] * t_0_xxyz_0_zzz[i];

        t_0_xxyz_z_yzz[i] = t_0_xxyz_0_yzzz[i] - rcd_z[i] * t_0_xxyz_0_yzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyy_0_xxz, t_0_xxyy_0_xxzz, t_0_xxyy_0_xyy,\
                           t_0_xxyy_0_xyyy, t_0_xxyy_0_xyyz, t_0_xxyy_0_xyz, t_0_xxyy_0_xyzz,\
                           t_0_xxyy_0_xzz, t_0_xxyy_0_xzzz, t_0_xxyy_0_yyy, t_0_xxyy_0_yyyy,\
                           t_0_xxyy_0_yyyz, t_0_xxyy_0_yyz, t_0_xxyy_0_yyzz, t_0_xxyy_0_yzz,\
                           t_0_xxyy_0_yzzz, t_0_xxyy_0_zzz, t_0_xxyy_0_zzzz, t_0_xxyy_y_xyy,\
                           t_0_xxyy_y_xyz, t_0_xxyy_y_xzz, t_0_xxyy_y_yyy, t_0_xxyy_y_yyz,\
                           t_0_xxyy_y_yzz, t_0_xxyy_y_zzz, t_0_xxyy_z_xxz, t_0_xxyy_z_xyz,\
                           t_0_xxyy_z_xzz, t_0_xxyy_z_yyz, t_0_xxyy_z_yzz, t_0_xxyy_z_zzz,\
                           t_0_xxyz_0_xxx, t_0_xxyz_0_xxxx, t_0_xxyz_0_xxxy, t_0_xxyz_0_xxxz,\
                           t_0_xxyz_0_xxy, t_0_xxyz_0_xxyy, t_0_xxyz_0_xxyz, t_0_xxyz_0_xxz,\
                           t_0_xxyz_0_xxzz, t_0_xxyz_0_xyy, t_0_xxyz_0_xyyy, t_0_xxyz_0_xyyz,\
                           t_0_xxyz_0_xyz, t_0_xxyz_0_xyzz, t_0_xxyz_0_xzz, t_0_xxyz_0_xzzz,\
                           t_0_xxyz_0_yyy, t_0_xxyz_0_yyyy, t_0_xxyz_0_yyyz, t_0_xxyz_0_yyz,\
                           t_0_xxyz_0_yyzz, t_0_xxyz_0_yzz, t_0_xxyz_0_yzzz, t_0_xxyz_0_zzz,\
                           t_0_xxyz_x_xxx, t_0_xxyz_x_xxy, t_0_xxyz_x_xxz, t_0_xxyz_x_xyy,\
                           t_0_xxyz_x_xyz, t_0_xxyz_x_xzz, t_0_xxyz_x_yyy, t_0_xxyz_x_yyz,\
                           t_0_xxyz_x_yzz, t_0_xxyz_x_zzz, t_0_xxyz_y_xxy, t_0_xxyz_y_xxz,\
                           t_0_xxyz_y_xyy, t_0_xxyz_y_xyz, t_0_xxyz_y_xzz, t_0_xxyz_y_yyy,\
                           t_0_xxyz_y_yyz, t_0_xxyz_y_yzz, t_0_xxyz_y_zzz, t_0_xxyz_z_xxz,\
                           t_0_xxyz_z_xyz, t_0_xxyz_z_xzz, t_0_xxyz_z_yyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxyz_z_yyz[i] = t_0_xxyz_0_yyzz[i] - rcd_z[i] * t_0_xxyz_0_yyz[i];

        t_0_xxyz_z_xzz[i] = t_0_xxyz_0_xzzz[i] - rcd_z[i] * t_0_xxyz_0_xzz[i];

        t_0_xxyz_z_xyz[i] = t_0_xxyz_0_xyzz[i] - rcd_z[i] * t_0_xxyz_0_xyz[i];

        t_0_xxyz_z_xxz[i] = t_0_xxyz_0_xxzz[i] - rcd_z[i] * t_0_xxyz_0_xxz[i];

        t_0_xxyz_y_zzz[i] = t_0_xxyz_0_yzzz[i] - rcd_y[i] * t_0_xxyz_0_zzz[i];

        t_0_xxyz_y_yzz[i] = t_0_xxyz_0_yyzz[i] - rcd_y[i] * t_0_xxyz_0_yzz[i];

        t_0_xxyz_y_yyz[i] = t_0_xxyz_0_yyyz[i] - rcd_y[i] * t_0_xxyz_0_yyz[i];

        t_0_xxyz_y_yyy[i] = t_0_xxyz_0_yyyy[i] - rcd_y[i] * t_0_xxyz_0_yyy[i];

        t_0_xxyz_y_xzz[i] = t_0_xxyz_0_xyzz[i] - rcd_y[i] * t_0_xxyz_0_xzz[i];

        t_0_xxyz_y_xyz[i] = t_0_xxyz_0_xyyz[i] - rcd_y[i] * t_0_xxyz_0_xyz[i];

        t_0_xxyz_y_xyy[i] = t_0_xxyz_0_xyyy[i] - rcd_y[i] * t_0_xxyz_0_xyy[i];

        t_0_xxyz_y_xxz[i] = t_0_xxyz_0_xxyz[i] - rcd_y[i] * t_0_xxyz_0_xxz[i];

        t_0_xxyz_y_xxy[i] = t_0_xxyz_0_xxyy[i] - rcd_y[i] * t_0_xxyz_0_xxy[i];

        t_0_xxyz_x_zzz[i] = t_0_xxyz_0_xzzz[i] - rcd_x[i] * t_0_xxyz_0_zzz[i];

        t_0_xxyz_x_yzz[i] = t_0_xxyz_0_xyzz[i] - rcd_x[i] * t_0_xxyz_0_yzz[i];

        t_0_xxyz_x_yyz[i] = t_0_xxyz_0_xyyz[i] - rcd_x[i] * t_0_xxyz_0_yyz[i];

        t_0_xxyz_x_yyy[i] = t_0_xxyz_0_xyyy[i] - rcd_x[i] * t_0_xxyz_0_yyy[i];

        t_0_xxyz_x_xzz[i] = t_0_xxyz_0_xxzz[i] - rcd_x[i] * t_0_xxyz_0_xzz[i];

        t_0_xxyz_x_xyz[i] = t_0_xxyz_0_xxyz[i] - rcd_x[i] * t_0_xxyz_0_xyz[i];

        t_0_xxyz_x_xyy[i] = t_0_xxyz_0_xxyy[i] - rcd_x[i] * t_0_xxyz_0_xyy[i];

        t_0_xxyz_x_xxz[i] = t_0_xxyz_0_xxxz[i] - rcd_x[i] * t_0_xxyz_0_xxz[i];

        t_0_xxyz_x_xxy[i] = t_0_xxyz_0_xxxy[i] - rcd_x[i] * t_0_xxyz_0_xxy[i];

        t_0_xxyz_x_xxx[i] = t_0_xxyz_0_xxxx[i] - rcd_x[i] * t_0_xxyz_0_xxx[i];

        t_0_xxyy_z_zzz[i] = t_0_xxyy_0_zzzz[i] - rcd_z[i] * t_0_xxyy_0_zzz[i];

        t_0_xxyy_z_yzz[i] = t_0_xxyy_0_yzzz[i] - rcd_z[i] * t_0_xxyy_0_yzz[i];

        t_0_xxyy_z_yyz[i] = t_0_xxyy_0_yyzz[i] - rcd_z[i] * t_0_xxyy_0_yyz[i];

        t_0_xxyy_z_xzz[i] = t_0_xxyy_0_xzzz[i] - rcd_z[i] * t_0_xxyy_0_xzz[i];

        t_0_xxyy_z_xyz[i] = t_0_xxyy_0_xyzz[i] - rcd_z[i] * t_0_xxyy_0_xyz[i];

        t_0_xxyy_z_xxz[i] = t_0_xxyy_0_xxzz[i] - rcd_z[i] * t_0_xxyy_0_xxz[i];

        t_0_xxyy_y_zzz[i] = t_0_xxyy_0_yzzz[i] - rcd_y[i] * t_0_xxyy_0_zzz[i];

        t_0_xxyy_y_yzz[i] = t_0_xxyy_0_yyzz[i] - rcd_y[i] * t_0_xxyy_0_yzz[i];

        t_0_xxyy_y_yyz[i] = t_0_xxyy_0_yyyz[i] - rcd_y[i] * t_0_xxyy_0_yyz[i];

        t_0_xxyy_y_yyy[i] = t_0_xxyy_0_yyyy[i] - rcd_y[i] * t_0_xxyy_0_yyy[i];

        t_0_xxyy_y_xzz[i] = t_0_xxyy_0_xyzz[i] - rcd_y[i] * t_0_xxyy_0_xzz[i];

        t_0_xxyy_y_xyz[i] = t_0_xxyy_0_xyyz[i] - rcd_y[i] * t_0_xxyy_0_xyz[i];

        t_0_xxyy_y_xyy[i] = t_0_xxyy_0_xyyy[i] - rcd_y[i] * t_0_xxyy_0_xyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxz_0_xxxy, t_0_xxxz_0_xxxz, t_0_xxxz_0_xxy,\
                           t_0_xxxz_0_xxyy, t_0_xxxz_0_xxyz, t_0_xxxz_0_xxz, t_0_xxxz_0_xxzz,\
                           t_0_xxxz_0_xyy, t_0_xxxz_0_xyyy, t_0_xxxz_0_xyyz, t_0_xxxz_0_xyz,\
                           t_0_xxxz_0_xyzz, t_0_xxxz_0_xzz, t_0_xxxz_0_xzzz, t_0_xxxz_0_yyy,\
                           t_0_xxxz_0_yyyy, t_0_xxxz_0_yyyz, t_0_xxxz_0_yyz, t_0_xxxz_0_yyzz,\
                           t_0_xxxz_0_yzz, t_0_xxxz_0_yzzz, t_0_xxxz_0_zzz, t_0_xxxz_0_zzzz,\
                           t_0_xxxz_x_xxy, t_0_xxxz_x_xxz, t_0_xxxz_x_xyy, t_0_xxxz_x_xyz,\
                           t_0_xxxz_x_xzz, t_0_xxxz_x_yyy, t_0_xxxz_x_yyz, t_0_xxxz_x_yzz,\
                           t_0_xxxz_x_zzz, t_0_xxxz_y_xxy, t_0_xxxz_y_xxz, t_0_xxxz_y_xyy,\
                           t_0_xxxz_y_xyz, t_0_xxxz_y_xzz, t_0_xxxz_y_yyy, t_0_xxxz_y_yyz,\
                           t_0_xxxz_y_yzz, t_0_xxxz_y_zzz, t_0_xxxz_z_xxz, t_0_xxxz_z_xyz,\
                           t_0_xxxz_z_xzz, t_0_xxxz_z_yyz, t_0_xxxz_z_yzz, t_0_xxxz_z_zzz,\
                           t_0_xxyy_0_xxx, t_0_xxyy_0_xxxx, t_0_xxyy_0_xxxy, t_0_xxyy_0_xxxz,\
                           t_0_xxyy_0_xxy, t_0_xxyy_0_xxyy, t_0_xxyy_0_xxyz, t_0_xxyy_0_xxz,\
                           t_0_xxyy_0_xxzz, t_0_xxyy_0_xyy, t_0_xxyy_0_xyyy, t_0_xxyy_0_xyyz,\
                           t_0_xxyy_0_xyz, t_0_xxyy_0_xyzz, t_0_xxyy_0_xzz, t_0_xxyy_0_xzzz,\
                           t_0_xxyy_0_yyy, t_0_xxyy_0_yyz, t_0_xxyy_0_yzz, t_0_xxyy_0_zzz,\
                           t_0_xxyy_x_xxx, t_0_xxyy_x_xxy, t_0_xxyy_x_xxz, t_0_xxyy_x_xyy,\
                           t_0_xxyy_x_xyz, t_0_xxyy_x_xzz, t_0_xxyy_x_yyy, t_0_xxyy_x_yyz,\
                           t_0_xxyy_x_yzz, t_0_xxyy_x_zzz, t_0_xxyy_y_xxy, t_0_xxyy_y_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxyy_y_xxz[i] = t_0_xxyy_0_xxyz[i] - rcd_y[i] * t_0_xxyy_0_xxz[i];

        t_0_xxyy_y_xxy[i] = t_0_xxyy_0_xxyy[i] - rcd_y[i] * t_0_xxyy_0_xxy[i];

        t_0_xxyy_x_zzz[i] = t_0_xxyy_0_xzzz[i] - rcd_x[i] * t_0_xxyy_0_zzz[i];

        t_0_xxyy_x_yzz[i] = t_0_xxyy_0_xyzz[i] - rcd_x[i] * t_0_xxyy_0_yzz[i];

        t_0_xxyy_x_yyz[i] = t_0_xxyy_0_xyyz[i] - rcd_x[i] * t_0_xxyy_0_yyz[i];

        t_0_xxyy_x_yyy[i] = t_0_xxyy_0_xyyy[i] - rcd_x[i] * t_0_xxyy_0_yyy[i];

        t_0_xxyy_x_xzz[i] = t_0_xxyy_0_xxzz[i] - rcd_x[i] * t_0_xxyy_0_xzz[i];

        t_0_xxyy_x_xyz[i] = t_0_xxyy_0_xxyz[i] - rcd_x[i] * t_0_xxyy_0_xyz[i];

        t_0_xxyy_x_xyy[i] = t_0_xxyy_0_xxyy[i] - rcd_x[i] * t_0_xxyy_0_xyy[i];

        t_0_xxyy_x_xxz[i] = t_0_xxyy_0_xxxz[i] - rcd_x[i] * t_0_xxyy_0_xxz[i];

        t_0_xxyy_x_xxy[i] = t_0_xxyy_0_xxxy[i] - rcd_x[i] * t_0_xxyy_0_xxy[i];

        t_0_xxyy_x_xxx[i] = t_0_xxyy_0_xxxx[i] - rcd_x[i] * t_0_xxyy_0_xxx[i];

        t_0_xxxz_z_zzz[i] = t_0_xxxz_0_zzzz[i] - rcd_z[i] * t_0_xxxz_0_zzz[i];

        t_0_xxxz_z_yzz[i] = t_0_xxxz_0_yzzz[i] - rcd_z[i] * t_0_xxxz_0_yzz[i];

        t_0_xxxz_z_yyz[i] = t_0_xxxz_0_yyzz[i] - rcd_z[i] * t_0_xxxz_0_yyz[i];

        t_0_xxxz_z_xzz[i] = t_0_xxxz_0_xzzz[i] - rcd_z[i] * t_0_xxxz_0_xzz[i];

        t_0_xxxz_z_xyz[i] = t_0_xxxz_0_xyzz[i] - rcd_z[i] * t_0_xxxz_0_xyz[i];

        t_0_xxxz_z_xxz[i] = t_0_xxxz_0_xxzz[i] - rcd_z[i] * t_0_xxxz_0_xxz[i];

        t_0_xxxz_y_zzz[i] = t_0_xxxz_0_yzzz[i] - rcd_y[i] * t_0_xxxz_0_zzz[i];

        t_0_xxxz_y_yzz[i] = t_0_xxxz_0_yyzz[i] - rcd_y[i] * t_0_xxxz_0_yzz[i];

        t_0_xxxz_y_yyz[i] = t_0_xxxz_0_yyyz[i] - rcd_y[i] * t_0_xxxz_0_yyz[i];

        t_0_xxxz_y_yyy[i] = t_0_xxxz_0_yyyy[i] - rcd_y[i] * t_0_xxxz_0_yyy[i];

        t_0_xxxz_y_xzz[i] = t_0_xxxz_0_xyzz[i] - rcd_y[i] * t_0_xxxz_0_xzz[i];

        t_0_xxxz_y_xyz[i] = t_0_xxxz_0_xyyz[i] - rcd_y[i] * t_0_xxxz_0_xyz[i];

        t_0_xxxz_y_xyy[i] = t_0_xxxz_0_xyyy[i] - rcd_y[i] * t_0_xxxz_0_xyy[i];

        t_0_xxxz_y_xxz[i] = t_0_xxxz_0_xxyz[i] - rcd_y[i] * t_0_xxxz_0_xxz[i];

        t_0_xxxz_y_xxy[i] = t_0_xxxz_0_xxyy[i] - rcd_y[i] * t_0_xxxz_0_xxy[i];

        t_0_xxxz_x_zzz[i] = t_0_xxxz_0_xzzz[i] - rcd_x[i] * t_0_xxxz_0_zzz[i];

        t_0_xxxz_x_yzz[i] = t_0_xxxz_0_xyzz[i] - rcd_x[i] * t_0_xxxz_0_yzz[i];

        t_0_xxxz_x_yyz[i] = t_0_xxxz_0_xyyz[i] - rcd_x[i] * t_0_xxxz_0_yyz[i];

        t_0_xxxz_x_yyy[i] = t_0_xxxz_0_xyyy[i] - rcd_x[i] * t_0_xxxz_0_yyy[i];

        t_0_xxxz_x_xzz[i] = t_0_xxxz_0_xxzz[i] - rcd_x[i] * t_0_xxxz_0_xzz[i];

        t_0_xxxz_x_xyz[i] = t_0_xxxz_0_xxyz[i] - rcd_x[i] * t_0_xxxz_0_xyz[i];

        t_0_xxxz_x_xyy[i] = t_0_xxxz_0_xxyy[i] - rcd_x[i] * t_0_xxxz_0_xyy[i];

        t_0_xxxz_x_xxz[i] = t_0_xxxz_0_xxxz[i] - rcd_x[i] * t_0_xxxz_0_xxz[i];

        t_0_xxxz_x_xxy[i] = t_0_xxxz_0_xxxy[i] - rcd_x[i] * t_0_xxxz_0_xxy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_0_xxz, t_0_xxxx_0_xxzz, t_0_xxxx_0_xyz,\
                           t_0_xxxx_0_xyzz, t_0_xxxx_0_xzz, t_0_xxxx_0_xzzz, t_0_xxxx_0_yyy,\
                           t_0_xxxx_0_yyyy, t_0_xxxx_0_yyyz, t_0_xxxx_0_yyz, t_0_xxxx_0_yyzz,\
                           t_0_xxxx_0_yzz, t_0_xxxx_0_yzzz, t_0_xxxx_0_zzz, t_0_xxxx_0_zzzz,\
                           t_0_xxxx_y_yyy, t_0_xxxx_y_yyz, t_0_xxxx_y_yzz, t_0_xxxx_y_zzz,\
                           t_0_xxxx_z_xxz, t_0_xxxx_z_xyz, t_0_xxxx_z_xzz, t_0_xxxx_z_yyz,\
                           t_0_xxxx_z_yzz, t_0_xxxx_z_zzz, t_0_xxxy_0_xxx, t_0_xxxy_0_xxxx,\
                           t_0_xxxy_0_xxxy, t_0_xxxy_0_xxxz, t_0_xxxy_0_xxy, t_0_xxxy_0_xxyy,\
                           t_0_xxxy_0_xxyz, t_0_xxxy_0_xxz, t_0_xxxy_0_xxzz, t_0_xxxy_0_xyy,\
                           t_0_xxxy_0_xyyy, t_0_xxxy_0_xyyz, t_0_xxxy_0_xyz, t_0_xxxy_0_xyzz,\
                           t_0_xxxy_0_xzz, t_0_xxxy_0_xzzz, t_0_xxxy_0_yyy, t_0_xxxy_0_yyyy,\
                           t_0_xxxy_0_yyyz, t_0_xxxy_0_yyz, t_0_xxxy_0_yyzz, t_0_xxxy_0_yzz,\
                           t_0_xxxy_0_yzzz, t_0_xxxy_0_zzz, t_0_xxxy_0_zzzz, t_0_xxxy_x_xxx,\
                           t_0_xxxy_x_xxy, t_0_xxxy_x_xxz, t_0_xxxy_x_xyy, t_0_xxxy_x_xyz,\
                           t_0_xxxy_x_xzz, t_0_xxxy_x_yyy, t_0_xxxy_x_yyz, t_0_xxxy_x_yzz,\
                           t_0_xxxy_x_zzz, t_0_xxxy_y_xxy, t_0_xxxy_y_xxz, t_0_xxxy_y_xyy,\
                           t_0_xxxy_y_xyz, t_0_xxxy_y_xzz, t_0_xxxy_y_yyy, t_0_xxxy_y_yyz,\
                           t_0_xxxy_y_yzz, t_0_xxxy_y_zzz, t_0_xxxy_z_xxz, t_0_xxxy_z_xyz,\
                           t_0_xxxy_z_xzz, t_0_xxxy_z_yyz, t_0_xxxy_z_yzz, t_0_xxxy_z_zzz,\
                           t_0_xxxz_0_xxx, t_0_xxxz_0_xxxx, t_0_xxxz_x_xxx : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxz_x_xxx[i] = t_0_xxxz_0_xxxx[i] - rcd_x[i] * t_0_xxxz_0_xxx[i];

        t_0_xxxy_z_zzz[i] = t_0_xxxy_0_zzzz[i] - rcd_z[i] * t_0_xxxy_0_zzz[i];

        t_0_xxxy_z_yzz[i] = t_0_xxxy_0_yzzz[i] - rcd_z[i] * t_0_xxxy_0_yzz[i];

        t_0_xxxy_z_yyz[i] = t_0_xxxy_0_yyzz[i] - rcd_z[i] * t_0_xxxy_0_yyz[i];

        t_0_xxxy_z_xzz[i] = t_0_xxxy_0_xzzz[i] - rcd_z[i] * t_0_xxxy_0_xzz[i];

        t_0_xxxy_z_xyz[i] = t_0_xxxy_0_xyzz[i] - rcd_z[i] * t_0_xxxy_0_xyz[i];

        t_0_xxxy_z_xxz[i] = t_0_xxxy_0_xxzz[i] - rcd_z[i] * t_0_xxxy_0_xxz[i];

        t_0_xxxy_y_zzz[i] = t_0_xxxy_0_yzzz[i] - rcd_y[i] * t_0_xxxy_0_zzz[i];

        t_0_xxxy_y_yzz[i] = t_0_xxxy_0_yyzz[i] - rcd_y[i] * t_0_xxxy_0_yzz[i];

        t_0_xxxy_y_yyz[i] = t_0_xxxy_0_yyyz[i] - rcd_y[i] * t_0_xxxy_0_yyz[i];

        t_0_xxxy_y_yyy[i] = t_0_xxxy_0_yyyy[i] - rcd_y[i] * t_0_xxxy_0_yyy[i];

        t_0_xxxy_y_xzz[i] = t_0_xxxy_0_xyzz[i] - rcd_y[i] * t_0_xxxy_0_xzz[i];

        t_0_xxxy_y_xyz[i] = t_0_xxxy_0_xyyz[i] - rcd_y[i] * t_0_xxxy_0_xyz[i];

        t_0_xxxy_y_xyy[i] = t_0_xxxy_0_xyyy[i] - rcd_y[i] * t_0_xxxy_0_xyy[i];

        t_0_xxxy_y_xxz[i] = t_0_xxxy_0_xxyz[i] - rcd_y[i] * t_0_xxxy_0_xxz[i];

        t_0_xxxy_y_xxy[i] = t_0_xxxy_0_xxyy[i] - rcd_y[i] * t_0_xxxy_0_xxy[i];

        t_0_xxxy_x_zzz[i] = t_0_xxxy_0_xzzz[i] - rcd_x[i] * t_0_xxxy_0_zzz[i];

        t_0_xxxy_x_yzz[i] = t_0_xxxy_0_xyzz[i] - rcd_x[i] * t_0_xxxy_0_yzz[i];

        t_0_xxxy_x_yyz[i] = t_0_xxxy_0_xyyz[i] - rcd_x[i] * t_0_xxxy_0_yyz[i];

        t_0_xxxy_x_yyy[i] = t_0_xxxy_0_xyyy[i] - rcd_x[i] * t_0_xxxy_0_yyy[i];

        t_0_xxxy_x_xzz[i] = t_0_xxxy_0_xxzz[i] - rcd_x[i] * t_0_xxxy_0_xzz[i];

        t_0_xxxy_x_xyz[i] = t_0_xxxy_0_xxyz[i] - rcd_x[i] * t_0_xxxy_0_xyz[i];

        t_0_xxxy_x_xyy[i] = t_0_xxxy_0_xxyy[i] - rcd_x[i] * t_0_xxxy_0_xyy[i];

        t_0_xxxy_x_xxz[i] = t_0_xxxy_0_xxxz[i] - rcd_x[i] * t_0_xxxy_0_xxz[i];

        t_0_xxxy_x_xxy[i] = t_0_xxxy_0_xxxy[i] - rcd_x[i] * t_0_xxxy_0_xxy[i];

        t_0_xxxy_x_xxx[i] = t_0_xxxy_0_xxxx[i] - rcd_x[i] * t_0_xxxy_0_xxx[i];

        t_0_xxxx_z_zzz[i] = t_0_xxxx_0_zzzz[i] - rcd_z[i] * t_0_xxxx_0_zzz[i];

        t_0_xxxx_z_yzz[i] = t_0_xxxx_0_yzzz[i] - rcd_z[i] * t_0_xxxx_0_yzz[i];

        t_0_xxxx_z_yyz[i] = t_0_xxxx_0_yyzz[i] - rcd_z[i] * t_0_xxxx_0_yyz[i];

        t_0_xxxx_z_xzz[i] = t_0_xxxx_0_xzzz[i] - rcd_z[i] * t_0_xxxx_0_xzz[i];

        t_0_xxxx_z_xyz[i] = t_0_xxxx_0_xyzz[i] - rcd_z[i] * t_0_xxxx_0_xyz[i];

        t_0_xxxx_z_xxz[i] = t_0_xxxx_0_xxzz[i] - rcd_z[i] * t_0_xxxx_0_xxz[i];

        t_0_xxxx_y_zzz[i] = t_0_xxxx_0_yzzz[i] - rcd_y[i] * t_0_xxxx_0_zzz[i];

        t_0_xxxx_y_yzz[i] = t_0_xxxx_0_yyzz[i] - rcd_y[i] * t_0_xxxx_0_yzz[i];

        t_0_xxxx_y_yyz[i] = t_0_xxxx_0_yyyz[i] - rcd_y[i] * t_0_xxxx_0_yyz[i];

        t_0_xxxx_y_yyy[i] = t_0_xxxx_0_yyyy[i] - rcd_y[i] * t_0_xxxx_0_yyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, t_0_xxxx_0_xxx, t_0_xxxx_0_xxxx, t_0_xxxx_0_xxxy,\
                           t_0_xxxx_0_xxxz, t_0_xxxx_0_xxy, t_0_xxxx_0_xxyy, t_0_xxxx_0_xxyz,\
                           t_0_xxxx_0_xxz, t_0_xxxx_0_xxzz, t_0_xxxx_0_xyy, t_0_xxxx_0_xyyy,\
                           t_0_xxxx_0_xyyz, t_0_xxxx_0_xyz, t_0_xxxx_0_xyzz, t_0_xxxx_0_xzz,\
                           t_0_xxxx_0_xzzz, t_0_xxxx_0_yyy, t_0_xxxx_0_yyz, t_0_xxxx_0_yzz,\
                           t_0_xxxx_0_zzz, t_0_xxxx_x_xxx, t_0_xxxx_x_xxy, t_0_xxxx_x_xxz,\
                           t_0_xxxx_x_xyy, t_0_xxxx_x_xyz, t_0_xxxx_x_xzz, t_0_xxxx_x_yyy,\
                           t_0_xxxx_x_yyz, t_0_xxxx_x_yzz, t_0_xxxx_x_zzz, t_0_xxxx_y_xxy,\
                           t_0_xxxx_y_xxz, t_0_xxxx_y_xyy, t_0_xxxx_y_xyz, t_0_xxxx_y_xzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxx_y_xzz[i] = t_0_xxxx_0_xyzz[i] - rcd_y[i] * t_0_xxxx_0_xzz[i];

        t_0_xxxx_y_xyz[i] = t_0_xxxx_0_xyyz[i] - rcd_y[i] * t_0_xxxx_0_xyz[i];

        t_0_xxxx_y_xyy[i] = t_0_xxxx_0_xyyy[i] - rcd_y[i] * t_0_xxxx_0_xyy[i];

        t_0_xxxx_y_xxz[i] = t_0_xxxx_0_xxyz[i] - rcd_y[i] * t_0_xxxx_0_xxz[i];

        t_0_xxxx_y_xxy[i] = t_0_xxxx_0_xxyy[i] - rcd_y[i] * t_0_xxxx_0_xxy[i];

        t_0_xxxx_x_zzz[i] = t_0_xxxx_0_xzzz[i] - rcd_x[i] * t_0_xxxx_0_zzz[i];

        t_0_xxxx_x_yzz[i] = t_0_xxxx_0_xyzz[i] - rcd_x[i] * t_0_xxxx_0_yzz[i];

        t_0_xxxx_x_yyz[i] = t_0_xxxx_0_xyyz[i] - rcd_x[i] * t_0_xxxx_0_yyz[i];

        t_0_xxxx_x_yyy[i] = t_0_xxxx_0_xyyy[i] - rcd_x[i] * t_0_xxxx_0_yyy[i];

        t_0_xxxx_x_xzz[i] = t_0_xxxx_0_xxzz[i] - rcd_x[i] * t_0_xxxx_0_xzz[i];

        t_0_xxxx_x_xyz[i] = t_0_xxxx_0_xxyz[i] - rcd_x[i] * t_0_xxxx_0_xyz[i];

        t_0_xxxx_x_xyy[i] = t_0_xxxx_0_xxyy[i] - rcd_x[i] * t_0_xxxx_0_xyy[i];

        t_0_xxxx_x_xxz[i] = t_0_xxxx_0_xxxz[i] - rcd_x[i] * t_0_xxxx_0_xxz[i];

        t_0_xxxx_x_xxy[i] = t_0_xxxx_0_xxxy[i] - rcd_x[i] * t_0_xxxx_0_xxy[i];

        t_0_xxxx_x_xxx[i] = t_0_xxxx_0_xxxx[i] - rcd_x[i] * t_0_xxxx_0_xxx[i];
    }
}


} // derirec namespace
