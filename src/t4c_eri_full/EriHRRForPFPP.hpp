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
compHostHRRForPFPP_V0(      BufferHostXY<T>&      intsBufferPFPP,
                      const BufferHostX<int32_t>& intsIndexesPFPP,
                      const BufferHostXY<T>&      intsBufferSFPP,
                      const BufferHostX<int32_t>& intsIndexesSFPP,
                      const BufferHostXY<T>&      intsBufferSGPP,
                      const BufferHostX<int32_t>& intsIndexesSGPP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SFPP) integral components

    t_0_zzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(0));

    t_0_zzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(1));

    t_0_zzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(2));

    t_0_zzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(3));

    t_0_zzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(4));

    t_0_zzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(5));

    t_0_zzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(6));

    t_0_zzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(7));

    t_0_zzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(8));

    t_0_yzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(9));

    t_0_yzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(10));

    t_0_yzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(11));

    t_0_yzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(12));

    t_0_yzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(13));

    t_0_yzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(14));

    t_0_yzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(15));

    t_0_yzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(16));

    t_0_yzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(17));

    t_0_yyz_z_z = intsBufferSFPP.data(intsIndexesSFPP(18));

    t_0_yyz_z_y = intsBufferSFPP.data(intsIndexesSFPP(19));

    t_0_yyz_z_x = intsBufferSFPP.data(intsIndexesSFPP(20));

    t_0_yyz_y_z = intsBufferSFPP.data(intsIndexesSFPP(21));

    t_0_yyz_y_y = intsBufferSFPP.data(intsIndexesSFPP(22));

    t_0_yyz_y_x = intsBufferSFPP.data(intsIndexesSFPP(23));

    t_0_yyz_x_z = intsBufferSFPP.data(intsIndexesSFPP(24));

    t_0_yyz_x_y = intsBufferSFPP.data(intsIndexesSFPP(25));

    t_0_yyz_x_x = intsBufferSFPP.data(intsIndexesSFPP(26));

    t_0_yyy_z_z = intsBufferSFPP.data(intsIndexesSFPP(27));

    t_0_yyy_z_y = intsBufferSFPP.data(intsIndexesSFPP(28));

    t_0_yyy_z_x = intsBufferSFPP.data(intsIndexesSFPP(29));

    t_0_yyy_y_z = intsBufferSFPP.data(intsIndexesSFPP(30));

    t_0_yyy_y_y = intsBufferSFPP.data(intsIndexesSFPP(31));

    t_0_yyy_y_x = intsBufferSFPP.data(intsIndexesSFPP(32));

    t_0_yyy_x_z = intsBufferSFPP.data(intsIndexesSFPP(33));

    t_0_yyy_x_y = intsBufferSFPP.data(intsIndexesSFPP(34));

    t_0_yyy_x_x = intsBufferSFPP.data(intsIndexesSFPP(35));

    t_0_xzz_z_z = intsBufferSFPP.data(intsIndexesSFPP(36));

    t_0_xzz_z_y = intsBufferSFPP.data(intsIndexesSFPP(37));

    t_0_xzz_z_x = intsBufferSFPP.data(intsIndexesSFPP(38));

    t_0_xzz_y_z = intsBufferSFPP.data(intsIndexesSFPP(39));

    t_0_xzz_y_y = intsBufferSFPP.data(intsIndexesSFPP(40));

    t_0_xzz_y_x = intsBufferSFPP.data(intsIndexesSFPP(41));

    t_0_xzz_x_z = intsBufferSFPP.data(intsIndexesSFPP(42));

    t_0_xzz_x_y = intsBufferSFPP.data(intsIndexesSFPP(43));

    t_0_xzz_x_x = intsBufferSFPP.data(intsIndexesSFPP(44));

    t_0_xyz_z_z = intsBufferSFPP.data(intsIndexesSFPP(45));

    t_0_xyz_z_y = intsBufferSFPP.data(intsIndexesSFPP(46));

    t_0_xyz_z_x = intsBufferSFPP.data(intsIndexesSFPP(47));

    t_0_xyz_y_z = intsBufferSFPP.data(intsIndexesSFPP(48));

    t_0_xyz_y_y = intsBufferSFPP.data(intsIndexesSFPP(49));

    t_0_xyz_y_x = intsBufferSFPP.data(intsIndexesSFPP(50));

    t_0_xyz_x_z = intsBufferSFPP.data(intsIndexesSFPP(51));

    t_0_xyz_x_y = intsBufferSFPP.data(intsIndexesSFPP(52));

    t_0_xyz_x_x = intsBufferSFPP.data(intsIndexesSFPP(53));

    t_0_xyy_z_z = intsBufferSFPP.data(intsIndexesSFPP(54));

    t_0_xyy_z_y = intsBufferSFPP.data(intsIndexesSFPP(55));

    t_0_xyy_z_x = intsBufferSFPP.data(intsIndexesSFPP(56));

    t_0_xyy_y_z = intsBufferSFPP.data(intsIndexesSFPP(57));

    t_0_xyy_y_y = intsBufferSFPP.data(intsIndexesSFPP(58));

    t_0_xyy_y_x = intsBufferSFPP.data(intsIndexesSFPP(59));

    t_0_xyy_x_z = intsBufferSFPP.data(intsIndexesSFPP(60));

    t_0_xyy_x_y = intsBufferSFPP.data(intsIndexesSFPP(61));

    t_0_xyy_x_x = intsBufferSFPP.data(intsIndexesSFPP(62));

    t_0_xxz_z_z = intsBufferSFPP.data(intsIndexesSFPP(63));

    t_0_xxz_z_y = intsBufferSFPP.data(intsIndexesSFPP(64));

    t_0_xxz_z_x = intsBufferSFPP.data(intsIndexesSFPP(65));

    t_0_xxz_y_z = intsBufferSFPP.data(intsIndexesSFPP(66));

    t_0_xxz_y_y = intsBufferSFPP.data(intsIndexesSFPP(67));

    t_0_xxz_y_x = intsBufferSFPP.data(intsIndexesSFPP(68));

    t_0_xxz_x_z = intsBufferSFPP.data(intsIndexesSFPP(69));

    t_0_xxz_x_y = intsBufferSFPP.data(intsIndexesSFPP(70));

    t_0_xxz_x_x = intsBufferSFPP.data(intsIndexesSFPP(71));

    t_0_xxy_z_z = intsBufferSFPP.data(intsIndexesSFPP(72));

    t_0_xxy_z_y = intsBufferSFPP.data(intsIndexesSFPP(73));

    t_0_xxy_z_x = intsBufferSFPP.data(intsIndexesSFPP(74));

    t_0_xxy_y_z = intsBufferSFPP.data(intsIndexesSFPP(75));

    t_0_xxy_y_y = intsBufferSFPP.data(intsIndexesSFPP(76));

    t_0_xxy_y_x = intsBufferSFPP.data(intsIndexesSFPP(77));

    t_0_xxy_x_z = intsBufferSFPP.data(intsIndexesSFPP(78));

    t_0_xxy_x_y = intsBufferSFPP.data(intsIndexesSFPP(79));

    t_0_xxy_x_x = intsBufferSFPP.data(intsIndexesSFPP(80));

    t_0_xxx_z_z = intsBufferSFPP.data(intsIndexesSFPP(81));

    t_0_xxx_z_y = intsBufferSFPP.data(intsIndexesSFPP(82));

    t_0_xxx_z_x = intsBufferSFPP.data(intsIndexesSFPP(83));

    t_0_xxx_y_z = intsBufferSFPP.data(intsIndexesSFPP(84));

    t_0_xxx_y_y = intsBufferSFPP.data(intsIndexesSFPP(85));

    t_0_xxx_y_x = intsBufferSFPP.data(intsIndexesSFPP(86));

    t_0_xxx_x_z = intsBufferSFPP.data(intsIndexesSFPP(87));

    t_0_xxx_x_y = intsBufferSFPP.data(intsIndexesSFPP(88));

    t_0_xxx_x_x = intsBufferSFPP.data(intsIndexesSFPP(89));

    // set up (SGPP) integral components

    t_0_zzzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(0));

    t_0_zzzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(1));

    t_0_zzzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(2));

    t_0_zzzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(3));

    t_0_zzzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(4));

    t_0_zzzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(5));

    t_0_zzzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(6));

    t_0_zzzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(7));

    t_0_zzzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(8));

    t_0_yzzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(9));

    t_0_yzzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(10));

    t_0_yzzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(11));

    t_0_yzzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(12));

    t_0_yzzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(13));

    t_0_yzzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(14));

    t_0_yzzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(15));

    t_0_yzzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(16));

    t_0_yzzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(17));

    t_0_yyzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(18));

    t_0_yyzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(19));

    t_0_yyzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(20));

    t_0_yyzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(21));

    t_0_yyzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(22));

    t_0_yyzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(23));

    t_0_yyzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(24));

    t_0_yyzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(25));

    t_0_yyzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(26));

    t_0_yyyz_z_z = intsBufferSGPP.data(intsIndexesSGPP(27));

    t_0_yyyz_z_y = intsBufferSGPP.data(intsIndexesSGPP(28));

    t_0_yyyz_z_x = intsBufferSGPP.data(intsIndexesSGPP(29));

    t_0_yyyz_y_z = intsBufferSGPP.data(intsIndexesSGPP(30));

    t_0_yyyz_y_y = intsBufferSGPP.data(intsIndexesSGPP(31));

    t_0_yyyz_y_x = intsBufferSGPP.data(intsIndexesSGPP(32));

    t_0_yyyz_x_z = intsBufferSGPP.data(intsIndexesSGPP(33));

    t_0_yyyz_x_y = intsBufferSGPP.data(intsIndexesSGPP(34));

    t_0_yyyz_x_x = intsBufferSGPP.data(intsIndexesSGPP(35));

    t_0_yyyy_z_z = intsBufferSGPP.data(intsIndexesSGPP(36));

    t_0_yyyy_z_y = intsBufferSGPP.data(intsIndexesSGPP(37));

    t_0_yyyy_z_x = intsBufferSGPP.data(intsIndexesSGPP(38));

    t_0_yyyy_y_z = intsBufferSGPP.data(intsIndexesSGPP(39));

    t_0_yyyy_y_y = intsBufferSGPP.data(intsIndexesSGPP(40));

    t_0_yyyy_y_x = intsBufferSGPP.data(intsIndexesSGPP(41));

    t_0_yyyy_x_z = intsBufferSGPP.data(intsIndexesSGPP(42));

    t_0_yyyy_x_y = intsBufferSGPP.data(intsIndexesSGPP(43));

    t_0_yyyy_x_x = intsBufferSGPP.data(intsIndexesSGPP(44));

    t_0_xzzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(45));

    t_0_xzzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(46));

    t_0_xzzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(47));

    t_0_xzzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(48));

    t_0_xzzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(49));

    t_0_xzzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(50));

    t_0_xzzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(51));

    t_0_xzzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(52));

    t_0_xzzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(53));

    t_0_xyzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(54));

    t_0_xyzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(55));

    t_0_xyzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(56));

    t_0_xyzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(57));

    t_0_xyzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(58));

    t_0_xyzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(59));

    t_0_xyzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(60));

    t_0_xyzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(61));

    t_0_xyzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(62));

    t_0_xyyz_z_z = intsBufferSGPP.data(intsIndexesSGPP(63));

    t_0_xyyz_z_y = intsBufferSGPP.data(intsIndexesSGPP(64));

    t_0_xyyz_z_x = intsBufferSGPP.data(intsIndexesSGPP(65));

    t_0_xyyz_y_z = intsBufferSGPP.data(intsIndexesSGPP(66));

    t_0_xyyz_y_y = intsBufferSGPP.data(intsIndexesSGPP(67));

    t_0_xyyz_y_x = intsBufferSGPP.data(intsIndexesSGPP(68));

    t_0_xyyz_x_z = intsBufferSGPP.data(intsIndexesSGPP(69));

    t_0_xyyz_x_y = intsBufferSGPP.data(intsIndexesSGPP(70));

    t_0_xyyz_x_x = intsBufferSGPP.data(intsIndexesSGPP(71));

    t_0_xyyy_z_z = intsBufferSGPP.data(intsIndexesSGPP(72));

    t_0_xyyy_z_y = intsBufferSGPP.data(intsIndexesSGPP(73));

    t_0_xyyy_z_x = intsBufferSGPP.data(intsIndexesSGPP(74));

    t_0_xyyy_y_z = intsBufferSGPP.data(intsIndexesSGPP(75));

    t_0_xyyy_y_y = intsBufferSGPP.data(intsIndexesSGPP(76));

    t_0_xyyy_y_x = intsBufferSGPP.data(intsIndexesSGPP(77));

    t_0_xyyy_x_z = intsBufferSGPP.data(intsIndexesSGPP(78));

    t_0_xyyy_x_y = intsBufferSGPP.data(intsIndexesSGPP(79));

    t_0_xyyy_x_x = intsBufferSGPP.data(intsIndexesSGPP(80));

    t_0_xxzz_z_z = intsBufferSGPP.data(intsIndexesSGPP(81));

    t_0_xxzz_z_y = intsBufferSGPP.data(intsIndexesSGPP(82));

    t_0_xxzz_z_x = intsBufferSGPP.data(intsIndexesSGPP(83));

    t_0_xxzz_y_z = intsBufferSGPP.data(intsIndexesSGPP(84));

    t_0_xxzz_y_y = intsBufferSGPP.data(intsIndexesSGPP(85));

    t_0_xxzz_y_x = intsBufferSGPP.data(intsIndexesSGPP(86));

    t_0_xxzz_x_z = intsBufferSGPP.data(intsIndexesSGPP(87));

    t_0_xxzz_x_y = intsBufferSGPP.data(intsIndexesSGPP(88));

    t_0_xxzz_x_x = intsBufferSGPP.data(intsIndexesSGPP(89));

    t_0_xxyz_z_z = intsBufferSGPP.data(intsIndexesSGPP(90));

    t_0_xxyz_z_y = intsBufferSGPP.data(intsIndexesSGPP(91));

    t_0_xxyz_z_x = intsBufferSGPP.data(intsIndexesSGPP(92));

    t_0_xxyz_y_z = intsBufferSGPP.data(intsIndexesSGPP(93));

    t_0_xxyz_y_y = intsBufferSGPP.data(intsIndexesSGPP(94));

    t_0_xxyz_y_x = intsBufferSGPP.data(intsIndexesSGPP(95));

    t_0_xxyz_x_z = intsBufferSGPP.data(intsIndexesSGPP(96));

    t_0_xxyz_x_y = intsBufferSGPP.data(intsIndexesSGPP(97));

    t_0_xxyz_x_x = intsBufferSGPP.data(intsIndexesSGPP(98));

    t_0_xxyy_z_z = intsBufferSGPP.data(intsIndexesSGPP(99));

    t_0_xxyy_z_y = intsBufferSGPP.data(intsIndexesSGPP(100));

    t_0_xxyy_z_x = intsBufferSGPP.data(intsIndexesSGPP(101));

    t_0_xxyy_y_z = intsBufferSGPP.data(intsIndexesSGPP(102));

    t_0_xxyy_y_y = intsBufferSGPP.data(intsIndexesSGPP(103));

    t_0_xxyy_y_x = intsBufferSGPP.data(intsIndexesSGPP(104));

    t_0_xxyy_x_z = intsBufferSGPP.data(intsIndexesSGPP(105));

    t_0_xxyy_x_y = intsBufferSGPP.data(intsIndexesSGPP(106));

    t_0_xxyy_x_x = intsBufferSGPP.data(intsIndexesSGPP(107));

    t_0_xxxz_z_z = intsBufferSGPP.data(intsIndexesSGPP(108));

    t_0_xxxz_z_y = intsBufferSGPP.data(intsIndexesSGPP(109));

    t_0_xxxz_z_x = intsBufferSGPP.data(intsIndexesSGPP(110));

    t_0_xxxz_y_z = intsBufferSGPP.data(intsIndexesSGPP(111));

    t_0_xxxz_y_y = intsBufferSGPP.data(intsIndexesSGPP(112));

    t_0_xxxz_y_x = intsBufferSGPP.data(intsIndexesSGPP(113));

    t_0_xxxz_x_z = intsBufferSGPP.data(intsIndexesSGPP(114));

    t_0_xxxz_x_y = intsBufferSGPP.data(intsIndexesSGPP(115));

    t_0_xxxz_x_x = intsBufferSGPP.data(intsIndexesSGPP(116));

    t_0_xxxy_z_z = intsBufferSGPP.data(intsIndexesSGPP(117));

    t_0_xxxy_z_y = intsBufferSGPP.data(intsIndexesSGPP(118));

    t_0_xxxy_z_x = intsBufferSGPP.data(intsIndexesSGPP(119));

    t_0_xxxy_y_z = intsBufferSGPP.data(intsIndexesSGPP(120));

    t_0_xxxy_y_y = intsBufferSGPP.data(intsIndexesSGPP(121));

    t_0_xxxy_y_x = intsBufferSGPP.data(intsIndexesSGPP(122));

    t_0_xxxy_x_z = intsBufferSGPP.data(intsIndexesSGPP(123));

    t_0_xxxy_x_y = intsBufferSGPP.data(intsIndexesSGPP(124));

    t_0_xxxy_x_x = intsBufferSGPP.data(intsIndexesSGPP(125));

    t_0_xxxx_z_z = intsBufferSGPP.data(intsIndexesSGPP(126));

    t_0_xxxx_z_y = intsBufferSGPP.data(intsIndexesSGPP(127));

    t_0_xxxx_z_x = intsBufferSGPP.data(intsIndexesSGPP(128));

    t_0_xxxx_y_z = intsBufferSGPP.data(intsIndexesSGPP(129));

    t_0_xxxx_y_y = intsBufferSGPP.data(intsIndexesSGPP(130));

    t_0_xxxx_y_x = intsBufferSGPP.data(intsIndexesSGPP(131));

    t_0_xxxx_x_z = intsBufferSGPP.data(intsIndexesSGPP(132));

    t_0_xxxx_x_y = intsBufferSGPP.data(intsIndexesSGPP(133));

    t_0_xxxx_x_x = intsBufferSGPP.data(intsIndexesSGPP(134));

    #pragma omp simd align(rab_z, t_0_xzz_x_x, t_0_xzz_x_y, t_0_xzz_x_z, t_0_xzz_y_x,\
                           t_0_xzz_y_y, t_0_xzz_y_z, t_0_xzz_z_x, t_0_xzz_z_y, t_0_xzz_z_z,\
                           t_0_xzzz_x_x, t_0_xzzz_x_y, t_0_xzzz_x_z, t_0_xzzz_y_x, t_0_xzzz_y_y,\
                           t_0_xzzz_y_z, t_0_xzzz_z_x, t_0_xzzz_z_y, t_0_xzzz_z_z, t_0_yyz_x_x,\
                           t_0_yyz_x_y, t_0_yyz_x_z, t_0_yyz_y_x, t_0_yyz_y_y, t_0_yyz_y_z,\
                           t_0_yyz_z_x, t_0_yyz_z_y, t_0_yyz_z_z, t_0_yyzz_x_x, t_0_yyzz_x_y,\
                           t_0_yyzz_x_z, t_0_yyzz_y_x, t_0_yyzz_y_y, t_0_yyzz_y_z, t_0_yyzz_z_x,\
                           t_0_yyzz_z_y, t_0_yyzz_z_z, t_0_yzz_x_x, t_0_yzz_x_y, t_0_yzz_x_z,\
                           t_0_yzz_y_x, t_0_yzz_y_y, t_0_yzz_y_z, t_0_yzz_z_x, t_0_yzz_z_y,\
                           t_0_yzz_z_z, t_0_yzzz_x_x, t_0_yzzz_x_y, t_0_yzzz_x_z, t_0_yzzz_y_x,\
                           t_0_yzzz_y_y, t_0_yzzz_y_z, t_0_yzzz_z_x, t_0_yzzz_z_y, t_0_yzzz_z_z,\
                           t_0_zzz_x_x, t_0_zzz_x_y, t_0_zzz_x_z, t_0_zzz_y_x, t_0_zzz_y_y,\
                           t_0_zzz_y_z, t_0_zzz_z_x, t_0_zzz_z_y, t_0_zzz_z_z, t_0_zzzz_x_x,\
                           t_0_zzzz_x_y, t_0_zzzz_x_z, t_0_zzzz_y_x, t_0_zzzz_y_y, t_0_zzzz_y_z,\
                           t_0_zzzz_z_x, t_0_zzzz_z_y, t_0_zzzz_z_z, t_z_xzz_x_x, t_z_xzz_x_y,\
                           t_z_xzz_x_z, t_z_xzz_y_x, t_z_xzz_y_y, t_z_xzz_y_z, t_z_xzz_z_x,\
                           t_z_xzz_z_y, t_z_xzz_z_z, t_z_yyz_x_x, t_z_yyz_x_y, t_z_yyz_x_z,\
                           t_z_yyz_y_x, t_z_yyz_y_y, t_z_yyz_y_z, t_z_yyz_z_x, t_z_yyz_z_y,\
                           t_z_yyz_z_z, t_z_yzz_x_x, t_z_yzz_x_y, t_z_yzz_x_z, t_z_yzz_y_x,\
                           t_z_yzz_y_y, t_z_yzz_y_z, t_z_yzz_z_x, t_z_yzz_z_y, t_z_yzz_z_z,\
                           t_z_zzz_x_x, t_z_zzz_x_y, t_z_zzz_x_z, t_z_zzz_y_x, t_z_zzz_y_y,\
                           t_z_zzz_y_z, t_z_zzz_z_x, t_z_zzz_z_y, t_z_zzz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_z_z[i] = t_0_zzzz_z_z[i] - rab_z[i] * t_0_zzz_z_z[i];

        t_z_zzz_z_y[i] = t_0_zzzz_z_y[i] - rab_z[i] * t_0_zzz_z_y[i];

        t_z_zzz_z_x[i] = t_0_zzzz_z_x[i] - rab_z[i] * t_0_zzz_z_x[i];

        t_z_zzz_y_z[i] = t_0_zzzz_y_z[i] - rab_z[i] * t_0_zzz_y_z[i];

        t_z_zzz_y_y[i] = t_0_zzzz_y_y[i] - rab_z[i] * t_0_zzz_y_y[i];

        t_z_zzz_y_x[i] = t_0_zzzz_y_x[i] - rab_z[i] * t_0_zzz_y_x[i];

        t_z_zzz_x_z[i] = t_0_zzzz_x_z[i] - rab_z[i] * t_0_zzz_x_z[i];

        t_z_zzz_x_y[i] = t_0_zzzz_x_y[i] - rab_z[i] * t_0_zzz_x_y[i];

        t_z_zzz_x_x[i] = t_0_zzzz_x_x[i] - rab_z[i] * t_0_zzz_x_x[i];

        t_z_yzz_z_z[i] = t_0_yzzz_z_z[i] - rab_z[i] * t_0_yzz_z_z[i];

        t_z_yzz_z_y[i] = t_0_yzzz_z_y[i] - rab_z[i] * t_0_yzz_z_y[i];

        t_z_yzz_z_x[i] = t_0_yzzz_z_x[i] - rab_z[i] * t_0_yzz_z_x[i];

        t_z_yzz_y_z[i] = t_0_yzzz_y_z[i] - rab_z[i] * t_0_yzz_y_z[i];

        t_z_yzz_y_y[i] = t_0_yzzz_y_y[i] - rab_z[i] * t_0_yzz_y_y[i];

        t_z_yzz_y_x[i] = t_0_yzzz_y_x[i] - rab_z[i] * t_0_yzz_y_x[i];

        t_z_yzz_x_z[i] = t_0_yzzz_x_z[i] - rab_z[i] * t_0_yzz_x_z[i];

        t_z_yzz_x_y[i] = t_0_yzzz_x_y[i] - rab_z[i] * t_0_yzz_x_y[i];

        t_z_yzz_x_x[i] = t_0_yzzz_x_x[i] - rab_z[i] * t_0_yzz_x_x[i];

        t_z_yyz_z_z[i] = t_0_yyzz_z_z[i] - rab_z[i] * t_0_yyz_z_z[i];

        t_z_yyz_z_y[i] = t_0_yyzz_z_y[i] - rab_z[i] * t_0_yyz_z_y[i];

        t_z_yyz_z_x[i] = t_0_yyzz_z_x[i] - rab_z[i] * t_0_yyz_z_x[i];

        t_z_yyz_y_z[i] = t_0_yyzz_y_z[i] - rab_z[i] * t_0_yyz_y_z[i];

        t_z_yyz_y_y[i] = t_0_yyzz_y_y[i] - rab_z[i] * t_0_yyz_y_y[i];

        t_z_yyz_y_x[i] = t_0_yyzz_y_x[i] - rab_z[i] * t_0_yyz_y_x[i];

        t_z_yyz_x_z[i] = t_0_yyzz_x_z[i] - rab_z[i] * t_0_yyz_x_z[i];

        t_z_yyz_x_y[i] = t_0_yyzz_x_y[i] - rab_z[i] * t_0_yyz_x_y[i];

        t_z_yyz_x_x[i] = t_0_yyzz_x_x[i] - rab_z[i] * t_0_yyz_x_x[i];

        t_z_xzz_z_z[i] = t_0_xzzz_z_z[i] - rab_z[i] * t_0_xzz_z_z[i];

        t_z_xzz_z_y[i] = t_0_xzzz_z_y[i] - rab_z[i] * t_0_xzz_z_y[i];

        t_z_xzz_z_x[i] = t_0_xzzz_z_x[i] - rab_z[i] * t_0_xzz_z_x[i];

        t_z_xzz_y_z[i] = t_0_xzzz_y_z[i] - rab_z[i] * t_0_xzz_y_z[i];

        t_z_xzz_y_y[i] = t_0_xzzz_y_y[i] - rab_z[i] * t_0_xzz_y_y[i];

        t_z_xzz_y_x[i] = t_0_xzzz_y_x[i] - rab_z[i] * t_0_xzz_y_x[i];

        t_z_xzz_x_z[i] = t_0_xzzz_x_z[i] - rab_z[i] * t_0_xzz_x_z[i];

        t_z_xzz_x_y[i] = t_0_xzzz_x_y[i] - rab_z[i] * t_0_xzz_x_y[i];

        t_z_xzz_x_x[i] = t_0_xzzz_x_x[i] - rab_z[i] * t_0_xzz_x_x[i];
    }

    #pragma omp simd align(rab_y, rab_z, t_0_xxz_x_x, t_0_xxz_x_y, t_0_xxz_x_z, t_0_xxz_y_x,\
                           t_0_xxz_y_y, t_0_xxz_y_z, t_0_xxz_z_x, t_0_xxz_z_y, t_0_xxz_z_z,\
                           t_0_xxzz_x_x, t_0_xxzz_x_y, t_0_xxzz_x_z, t_0_xxzz_y_x, t_0_xxzz_y_y,\
                           t_0_xxzz_y_z, t_0_xxzz_z_x, t_0_xxzz_z_y, t_0_xxzz_z_z, t_0_xyz_x_x,\
                           t_0_xyz_x_y, t_0_xyz_x_z, t_0_xyz_y_x, t_0_xyz_y_y, t_0_xyz_y_z,\
                           t_0_xyz_z_x, t_0_xyz_z_y, t_0_xyz_z_z, t_0_xyzz_x_x, t_0_xyzz_x_y,\
                           t_0_xyzz_x_z, t_0_xyzz_y_x, t_0_xyzz_y_y, t_0_xyzz_y_z, t_0_xyzz_z_x,\
                           t_0_xyzz_z_y, t_0_xyzz_z_z, t_0_yyzz_x_x, t_0_yyzz_x_y, t_0_yyzz_x_z,\
                           t_0_yyzz_y_x, t_0_yyzz_y_y, t_0_yyzz_y_z, t_0_yyzz_z_x, t_0_yyzz_z_y,\
                           t_0_yyzz_z_z, t_0_yzz_x_x, t_0_yzz_x_y, t_0_yzz_x_z, t_0_yzz_y_x,\
                           t_0_yzz_y_y, t_0_yzz_y_z, t_0_yzz_z_x, t_0_yzz_z_y, t_0_yzz_z_z,\
                           t_0_yzzz_x_x, t_0_yzzz_x_y, t_0_yzzz_x_z, t_0_yzzz_y_x, t_0_yzzz_y_y,\
                           t_0_yzzz_y_z, t_0_yzzz_z_x, t_0_yzzz_z_y, t_0_yzzz_z_z, t_0_zzz_x_x,\
                           t_0_zzz_x_y, t_0_zzz_x_z, t_0_zzz_y_x, t_0_zzz_y_y, t_0_zzz_y_z,\
                           t_0_zzz_z_x, t_0_zzz_z_y, t_0_zzz_z_z, t_y_yzz_x_x, t_y_yzz_x_y,\
                           t_y_yzz_x_z, t_y_yzz_y_x, t_y_yzz_y_y, t_y_yzz_y_z, t_y_yzz_z_x,\
                           t_y_yzz_z_y, t_y_yzz_z_z, t_y_zzz_x_x, t_y_zzz_x_y, t_y_zzz_x_z,\
                           t_y_zzz_y_x, t_y_zzz_y_y, t_y_zzz_y_z, t_y_zzz_z_x, t_y_zzz_z_y,\
                           t_y_zzz_z_z, t_z_xxz_x_x, t_z_xxz_x_y, t_z_xxz_x_z, t_z_xxz_y_x,\
                           t_z_xxz_y_y, t_z_xxz_y_z, t_z_xxz_z_x, t_z_xxz_z_y, t_z_xxz_z_z,\
                           t_z_xyz_x_x, t_z_xyz_x_y, t_z_xyz_x_z, t_z_xyz_y_x, t_z_xyz_y_y,\
                           t_z_xyz_y_z, t_z_xyz_z_x, t_z_xyz_z_y, t_z_xyz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xyz_z_z[i] = t_0_xyzz_z_z[i] - rab_z[i] * t_0_xyz_z_z[i];

        t_z_xyz_z_y[i] = t_0_xyzz_z_y[i] - rab_z[i] * t_0_xyz_z_y[i];

        t_z_xyz_z_x[i] = t_0_xyzz_z_x[i] - rab_z[i] * t_0_xyz_z_x[i];

        t_z_xyz_y_z[i] = t_0_xyzz_y_z[i] - rab_z[i] * t_0_xyz_y_z[i];

        t_z_xyz_y_y[i] = t_0_xyzz_y_y[i] - rab_z[i] * t_0_xyz_y_y[i];

        t_z_xyz_y_x[i] = t_0_xyzz_y_x[i] - rab_z[i] * t_0_xyz_y_x[i];

        t_z_xyz_x_z[i] = t_0_xyzz_x_z[i] - rab_z[i] * t_0_xyz_x_z[i];

        t_z_xyz_x_y[i] = t_0_xyzz_x_y[i] - rab_z[i] * t_0_xyz_x_y[i];

        t_z_xyz_x_x[i] = t_0_xyzz_x_x[i] - rab_z[i] * t_0_xyz_x_x[i];

        t_z_xxz_z_z[i] = t_0_xxzz_z_z[i] - rab_z[i] * t_0_xxz_z_z[i];

        t_z_xxz_z_y[i] = t_0_xxzz_z_y[i] - rab_z[i] * t_0_xxz_z_y[i];

        t_z_xxz_z_x[i] = t_0_xxzz_z_x[i] - rab_z[i] * t_0_xxz_z_x[i];

        t_z_xxz_y_z[i] = t_0_xxzz_y_z[i] - rab_z[i] * t_0_xxz_y_z[i];

        t_z_xxz_y_y[i] = t_0_xxzz_y_y[i] - rab_z[i] * t_0_xxz_y_y[i];

        t_z_xxz_y_x[i] = t_0_xxzz_y_x[i] - rab_z[i] * t_0_xxz_y_x[i];

        t_z_xxz_x_z[i] = t_0_xxzz_x_z[i] - rab_z[i] * t_0_xxz_x_z[i];

        t_z_xxz_x_y[i] = t_0_xxzz_x_y[i] - rab_z[i] * t_0_xxz_x_y[i];

        t_z_xxz_x_x[i] = t_0_xxzz_x_x[i] - rab_z[i] * t_0_xxz_x_x[i];

        t_y_zzz_z_z[i] = t_0_yzzz_z_z[i] - rab_y[i] * t_0_zzz_z_z[i];

        t_y_zzz_z_y[i] = t_0_yzzz_z_y[i] - rab_y[i] * t_0_zzz_z_y[i];

        t_y_zzz_z_x[i] = t_0_yzzz_z_x[i] - rab_y[i] * t_0_zzz_z_x[i];

        t_y_zzz_y_z[i] = t_0_yzzz_y_z[i] - rab_y[i] * t_0_zzz_y_z[i];

        t_y_zzz_y_y[i] = t_0_yzzz_y_y[i] - rab_y[i] * t_0_zzz_y_y[i];

        t_y_zzz_y_x[i] = t_0_yzzz_y_x[i] - rab_y[i] * t_0_zzz_y_x[i];

        t_y_zzz_x_z[i] = t_0_yzzz_x_z[i] - rab_y[i] * t_0_zzz_x_z[i];

        t_y_zzz_x_y[i] = t_0_yzzz_x_y[i] - rab_y[i] * t_0_zzz_x_y[i];

        t_y_zzz_x_x[i] = t_0_yzzz_x_x[i] - rab_y[i] * t_0_zzz_x_x[i];

        t_y_yzz_z_z[i] = t_0_yyzz_z_z[i] - rab_y[i] * t_0_yzz_z_z[i];

        t_y_yzz_z_y[i] = t_0_yyzz_z_y[i] - rab_y[i] * t_0_yzz_z_y[i];

        t_y_yzz_z_x[i] = t_0_yyzz_z_x[i] - rab_y[i] * t_0_yzz_z_x[i];

        t_y_yzz_y_z[i] = t_0_yyzz_y_z[i] - rab_y[i] * t_0_yzz_y_z[i];

        t_y_yzz_y_y[i] = t_0_yyzz_y_y[i] - rab_y[i] * t_0_yzz_y_y[i];

        t_y_yzz_y_x[i] = t_0_yyzz_y_x[i] - rab_y[i] * t_0_yzz_y_x[i];

        t_y_yzz_x_z[i] = t_0_yyzz_x_z[i] - rab_y[i] * t_0_yzz_x_z[i];

        t_y_yzz_x_y[i] = t_0_yyzz_x_y[i] - rab_y[i] * t_0_yzz_x_y[i];

        t_y_yzz_x_x[i] = t_0_yyzz_x_x[i] - rab_y[i] * t_0_yzz_x_x[i];
    }

    #pragma omp simd align(rab_y, t_0_xyyz_x_x, t_0_xyyz_x_y, t_0_xyyz_x_z, t_0_xyyz_y_x,\
                           t_0_xyyz_y_y, t_0_xyyz_y_z, t_0_xyyz_z_x, t_0_xyyz_z_y, t_0_xyyz_z_z,\
                           t_0_xyz_x_x, t_0_xyz_x_y, t_0_xyz_x_z, t_0_xyz_y_x, t_0_xyz_y_y,\
                           t_0_xyz_y_z, t_0_xyz_z_x, t_0_xyz_z_y, t_0_xyz_z_z, t_0_xyzz_x_x,\
                           t_0_xyzz_x_y, t_0_xyzz_x_z, t_0_xyzz_y_x, t_0_xyzz_y_y, t_0_xyzz_y_z,\
                           t_0_xyzz_z_x, t_0_xyzz_z_y, t_0_xyzz_z_z, t_0_xzz_x_x, t_0_xzz_x_y,\
                           t_0_xzz_x_z, t_0_xzz_y_x, t_0_xzz_y_y, t_0_xzz_y_z, t_0_xzz_z_x,\
                           t_0_xzz_z_y, t_0_xzz_z_z, t_0_yyy_x_x, t_0_yyy_x_y, t_0_yyy_x_z,\
                           t_0_yyy_y_x, t_0_yyy_y_y, t_0_yyy_y_z, t_0_yyy_z_x, t_0_yyy_z_y,\
                           t_0_yyy_z_z, t_0_yyyy_x_x, t_0_yyyy_x_y, t_0_yyyy_x_z, t_0_yyyy_y_x,\
                           t_0_yyyy_y_y, t_0_yyyy_y_z, t_0_yyyy_z_x, t_0_yyyy_z_y, t_0_yyyy_z_z,\
                           t_0_yyyz_x_x, t_0_yyyz_x_y, t_0_yyyz_x_z, t_0_yyyz_y_x, t_0_yyyz_y_y,\
                           t_0_yyyz_y_z, t_0_yyyz_z_x, t_0_yyyz_z_y, t_0_yyyz_z_z, t_0_yyz_x_x,\
                           t_0_yyz_x_y, t_0_yyz_x_z, t_0_yyz_y_x, t_0_yyz_y_y, t_0_yyz_y_z,\
                           t_0_yyz_z_x, t_0_yyz_z_y, t_0_yyz_z_z, t_y_xyz_x_x, t_y_xyz_x_y,\
                           t_y_xyz_x_z, t_y_xyz_y_x, t_y_xyz_y_y, t_y_xyz_y_z, t_y_xyz_z_x,\
                           t_y_xyz_z_y, t_y_xyz_z_z, t_y_xzz_x_x, t_y_xzz_x_y, t_y_xzz_x_z,\
                           t_y_xzz_y_x, t_y_xzz_y_y, t_y_xzz_y_z, t_y_xzz_z_x, t_y_xzz_z_y,\
                           t_y_xzz_z_z, t_y_yyy_x_x, t_y_yyy_x_y, t_y_yyy_x_z, t_y_yyy_y_x,\
                           t_y_yyy_y_y, t_y_yyy_y_z, t_y_yyy_z_x, t_y_yyy_z_y, t_y_yyy_z_z,\
                           t_y_yyz_x_x, t_y_yyz_x_y, t_y_yyz_x_z, t_y_yyz_y_x, t_y_yyz_y_y,\
                           t_y_yyz_y_z, t_y_yyz_z_x, t_y_yyz_z_y, t_y_yyz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yyz_z_z[i] = t_0_yyyz_z_z[i] - rab_y[i] * t_0_yyz_z_z[i];

        t_y_yyz_z_y[i] = t_0_yyyz_z_y[i] - rab_y[i] * t_0_yyz_z_y[i];

        t_y_yyz_z_x[i] = t_0_yyyz_z_x[i] - rab_y[i] * t_0_yyz_z_x[i];

        t_y_yyz_y_z[i] = t_0_yyyz_y_z[i] - rab_y[i] * t_0_yyz_y_z[i];

        t_y_yyz_y_y[i] = t_0_yyyz_y_y[i] - rab_y[i] * t_0_yyz_y_y[i];

        t_y_yyz_y_x[i] = t_0_yyyz_y_x[i] - rab_y[i] * t_0_yyz_y_x[i];

        t_y_yyz_x_z[i] = t_0_yyyz_x_z[i] - rab_y[i] * t_0_yyz_x_z[i];

        t_y_yyz_x_y[i] = t_0_yyyz_x_y[i] - rab_y[i] * t_0_yyz_x_y[i];

        t_y_yyz_x_x[i] = t_0_yyyz_x_x[i] - rab_y[i] * t_0_yyz_x_x[i];

        t_y_yyy_z_z[i] = t_0_yyyy_z_z[i] - rab_y[i] * t_0_yyy_z_z[i];

        t_y_yyy_z_y[i] = t_0_yyyy_z_y[i] - rab_y[i] * t_0_yyy_z_y[i];

        t_y_yyy_z_x[i] = t_0_yyyy_z_x[i] - rab_y[i] * t_0_yyy_z_x[i];

        t_y_yyy_y_z[i] = t_0_yyyy_y_z[i] - rab_y[i] * t_0_yyy_y_z[i];

        t_y_yyy_y_y[i] = t_0_yyyy_y_y[i] - rab_y[i] * t_0_yyy_y_y[i];

        t_y_yyy_y_x[i] = t_0_yyyy_y_x[i] - rab_y[i] * t_0_yyy_y_x[i];

        t_y_yyy_x_z[i] = t_0_yyyy_x_z[i] - rab_y[i] * t_0_yyy_x_z[i];

        t_y_yyy_x_y[i] = t_0_yyyy_x_y[i] - rab_y[i] * t_0_yyy_x_y[i];

        t_y_yyy_x_x[i] = t_0_yyyy_x_x[i] - rab_y[i] * t_0_yyy_x_x[i];

        t_y_xzz_z_z[i] = t_0_xyzz_z_z[i] - rab_y[i] * t_0_xzz_z_z[i];

        t_y_xzz_z_y[i] = t_0_xyzz_z_y[i] - rab_y[i] * t_0_xzz_z_y[i];

        t_y_xzz_z_x[i] = t_0_xyzz_z_x[i] - rab_y[i] * t_0_xzz_z_x[i];

        t_y_xzz_y_z[i] = t_0_xyzz_y_z[i] - rab_y[i] * t_0_xzz_y_z[i];

        t_y_xzz_y_y[i] = t_0_xyzz_y_y[i] - rab_y[i] * t_0_xzz_y_y[i];

        t_y_xzz_y_x[i] = t_0_xyzz_y_x[i] - rab_y[i] * t_0_xzz_y_x[i];

        t_y_xzz_x_z[i] = t_0_xyzz_x_z[i] - rab_y[i] * t_0_xzz_x_z[i];

        t_y_xzz_x_y[i] = t_0_xyzz_x_y[i] - rab_y[i] * t_0_xzz_x_y[i];

        t_y_xzz_x_x[i] = t_0_xyzz_x_x[i] - rab_y[i] * t_0_xzz_x_x[i];

        t_y_xyz_z_z[i] = t_0_xyyz_z_z[i] - rab_y[i] * t_0_xyz_z_z[i];

        t_y_xyz_z_y[i] = t_0_xyyz_z_y[i] - rab_y[i] * t_0_xyz_z_y[i];

        t_y_xyz_z_x[i] = t_0_xyyz_z_x[i] - rab_y[i] * t_0_xyz_z_x[i];

        t_y_xyz_y_z[i] = t_0_xyyz_y_z[i] - rab_y[i] * t_0_xyz_y_z[i];

        t_y_xyz_y_y[i] = t_0_xyyz_y_y[i] - rab_y[i] * t_0_xyz_y_y[i];

        t_y_xyz_y_x[i] = t_0_xyyz_y_x[i] - rab_y[i] * t_0_xyz_y_x[i];

        t_y_xyz_x_z[i] = t_0_xyyz_x_z[i] - rab_y[i] * t_0_xyz_x_z[i];

        t_y_xyz_x_y[i] = t_0_xyyz_x_y[i] - rab_y[i] * t_0_xyz_x_y[i];

        t_y_xyz_x_x[i] = t_0_xyyz_x_x[i] - rab_y[i] * t_0_xyz_x_x[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_0_xxy_x_x, t_0_xxy_x_y, t_0_xxy_x_z, t_0_xxy_y_x,\
                           t_0_xxy_y_y, t_0_xxy_y_z, t_0_xxy_z_x, t_0_xxy_z_y, t_0_xxy_z_z,\
                           t_0_xxyy_x_x, t_0_xxyy_x_y, t_0_xxyy_x_z, t_0_xxyy_y_x, t_0_xxyy_y_y,\
                           t_0_xxyy_y_z, t_0_xxyy_z_x, t_0_xxyy_z_y, t_0_xxyy_z_z, t_0_xxyz_x_x,\
                           t_0_xxyz_x_y, t_0_xxyz_x_z, t_0_xxyz_y_x, t_0_xxyz_y_y, t_0_xxyz_y_z,\
                           t_0_xxyz_z_x, t_0_xxyz_z_y, t_0_xxyz_z_z, t_0_xxz_x_x, t_0_xxz_x_y,\
                           t_0_xxz_x_z, t_0_xxz_y_x, t_0_xxz_y_y, t_0_xxz_y_z, t_0_xxz_z_x,\
                           t_0_xxz_z_y, t_0_xxz_z_z, t_0_xyy_x_x, t_0_xyy_x_y, t_0_xyy_x_z,\
                           t_0_xyy_y_x, t_0_xyy_y_y, t_0_xyy_y_z, t_0_xyy_z_x, t_0_xyy_z_y,\
                           t_0_xyy_z_z, t_0_xyyy_x_x, t_0_xyyy_x_y, t_0_xyyy_x_z, t_0_xyyy_y_x,\
                           t_0_xyyy_y_y, t_0_xyyy_y_z, t_0_xyyy_z_x, t_0_xyyy_z_y, t_0_xyyy_z_z,\
                           t_0_xzzz_x_x, t_0_xzzz_x_y, t_0_xzzz_x_z, t_0_xzzz_y_x, t_0_xzzz_y_y,\
                           t_0_xzzz_y_z, t_0_xzzz_z_x, t_0_xzzz_z_y, t_0_xzzz_z_z, t_0_zzz_x_x,\
                           t_0_zzz_x_y, t_0_zzz_x_z, t_0_zzz_y_x, t_0_zzz_y_y, t_0_zzz_y_z,\
                           t_0_zzz_z_x, t_0_zzz_z_y, t_0_zzz_z_z, t_x_zzz_x_x, t_x_zzz_x_y,\
                           t_x_zzz_x_z, t_x_zzz_y_x, t_x_zzz_y_y, t_x_zzz_y_z, t_x_zzz_z_x,\
                           t_x_zzz_z_y, t_x_zzz_z_z, t_y_xxy_x_x, t_y_xxy_x_y, t_y_xxy_x_z,\
                           t_y_xxy_y_x, t_y_xxy_y_y, t_y_xxy_y_z, t_y_xxy_z_x, t_y_xxy_z_y,\
                           t_y_xxy_z_z, t_y_xxz_x_x, t_y_xxz_x_y, t_y_xxz_x_z, t_y_xxz_y_x,\
                           t_y_xxz_y_y, t_y_xxz_y_z, t_y_xxz_z_x, t_y_xxz_z_y, t_y_xxz_z_z,\
                           t_y_xyy_x_x, t_y_xyy_x_y, t_y_xyy_x_z, t_y_xyy_y_x, t_y_xyy_y_y,\
                           t_y_xyy_y_z, t_y_xyy_z_x, t_y_xyy_z_y, t_y_xyy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyy_z_z[i] = t_0_xyyy_z_z[i] - rab_y[i] * t_0_xyy_z_z[i];

        t_y_xyy_z_y[i] = t_0_xyyy_z_y[i] - rab_y[i] * t_0_xyy_z_y[i];

        t_y_xyy_z_x[i] = t_0_xyyy_z_x[i] - rab_y[i] * t_0_xyy_z_x[i];

        t_y_xyy_y_z[i] = t_0_xyyy_y_z[i] - rab_y[i] * t_0_xyy_y_z[i];

        t_y_xyy_y_y[i] = t_0_xyyy_y_y[i] - rab_y[i] * t_0_xyy_y_y[i];

        t_y_xyy_y_x[i] = t_0_xyyy_y_x[i] - rab_y[i] * t_0_xyy_y_x[i];

        t_y_xyy_x_z[i] = t_0_xyyy_x_z[i] - rab_y[i] * t_0_xyy_x_z[i];

        t_y_xyy_x_y[i] = t_0_xyyy_x_y[i] - rab_y[i] * t_0_xyy_x_y[i];

        t_y_xyy_x_x[i] = t_0_xyyy_x_x[i] - rab_y[i] * t_0_xyy_x_x[i];

        t_y_xxz_z_z[i] = t_0_xxyz_z_z[i] - rab_y[i] * t_0_xxz_z_z[i];

        t_y_xxz_z_y[i] = t_0_xxyz_z_y[i] - rab_y[i] * t_0_xxz_z_y[i];

        t_y_xxz_z_x[i] = t_0_xxyz_z_x[i] - rab_y[i] * t_0_xxz_z_x[i];

        t_y_xxz_y_z[i] = t_0_xxyz_y_z[i] - rab_y[i] * t_0_xxz_y_z[i];

        t_y_xxz_y_y[i] = t_0_xxyz_y_y[i] - rab_y[i] * t_0_xxz_y_y[i];

        t_y_xxz_y_x[i] = t_0_xxyz_y_x[i] - rab_y[i] * t_0_xxz_y_x[i];

        t_y_xxz_x_z[i] = t_0_xxyz_x_z[i] - rab_y[i] * t_0_xxz_x_z[i];

        t_y_xxz_x_y[i] = t_0_xxyz_x_y[i] - rab_y[i] * t_0_xxz_x_y[i];

        t_y_xxz_x_x[i] = t_0_xxyz_x_x[i] - rab_y[i] * t_0_xxz_x_x[i];

        t_y_xxy_z_z[i] = t_0_xxyy_z_z[i] - rab_y[i] * t_0_xxy_z_z[i];

        t_y_xxy_z_y[i] = t_0_xxyy_z_y[i] - rab_y[i] * t_0_xxy_z_y[i];

        t_y_xxy_z_x[i] = t_0_xxyy_z_x[i] - rab_y[i] * t_0_xxy_z_x[i];

        t_y_xxy_y_z[i] = t_0_xxyy_y_z[i] - rab_y[i] * t_0_xxy_y_z[i];

        t_y_xxy_y_y[i] = t_0_xxyy_y_y[i] - rab_y[i] * t_0_xxy_y_y[i];

        t_y_xxy_y_x[i] = t_0_xxyy_y_x[i] - rab_y[i] * t_0_xxy_y_x[i];

        t_y_xxy_x_z[i] = t_0_xxyy_x_z[i] - rab_y[i] * t_0_xxy_x_z[i];

        t_y_xxy_x_y[i] = t_0_xxyy_x_y[i] - rab_y[i] * t_0_xxy_x_y[i];

        t_y_xxy_x_x[i] = t_0_xxyy_x_x[i] - rab_y[i] * t_0_xxy_x_x[i];

        t_x_zzz_z_z[i] = t_0_xzzz_z_z[i] - rab_x[i] * t_0_zzz_z_z[i];

        t_x_zzz_z_y[i] = t_0_xzzz_z_y[i] - rab_x[i] * t_0_zzz_z_y[i];

        t_x_zzz_z_x[i] = t_0_xzzz_z_x[i] - rab_x[i] * t_0_zzz_z_x[i];

        t_x_zzz_y_z[i] = t_0_xzzz_y_z[i] - rab_x[i] * t_0_zzz_y_z[i];

        t_x_zzz_y_y[i] = t_0_xzzz_y_y[i] - rab_x[i] * t_0_zzz_y_y[i];

        t_x_zzz_y_x[i] = t_0_xzzz_y_x[i] - rab_x[i] * t_0_zzz_y_x[i];

        t_x_zzz_x_z[i] = t_0_xzzz_x_z[i] - rab_x[i] * t_0_zzz_x_z[i];

        t_x_zzz_x_y[i] = t_0_xzzz_x_y[i] - rab_x[i] * t_0_zzz_x_y[i];

        t_x_zzz_x_x[i] = t_0_xzzz_x_x[i] - rab_x[i] * t_0_zzz_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xxzz_x_x, t_0_xxzz_x_y, t_0_xxzz_x_z, t_0_xxzz_y_x,\
                           t_0_xxzz_y_y, t_0_xxzz_y_z, t_0_xxzz_z_x, t_0_xxzz_z_y, t_0_xxzz_z_z,\
                           t_0_xyyy_x_x, t_0_xyyy_x_y, t_0_xyyy_x_z, t_0_xyyy_y_x, t_0_xyyy_y_y,\
                           t_0_xyyy_y_z, t_0_xyyy_z_x, t_0_xyyy_z_y, t_0_xyyy_z_z, t_0_xyyz_x_x,\
                           t_0_xyyz_x_y, t_0_xyyz_x_z, t_0_xyyz_y_x, t_0_xyyz_y_y, t_0_xyyz_y_z,\
                           t_0_xyyz_z_x, t_0_xyyz_z_y, t_0_xyyz_z_z, t_0_xyzz_x_x, t_0_xyzz_x_y,\
                           t_0_xyzz_x_z, t_0_xyzz_y_x, t_0_xyzz_y_y, t_0_xyzz_y_z, t_0_xyzz_z_x,\
                           t_0_xyzz_z_y, t_0_xyzz_z_z, t_0_xzz_x_x, t_0_xzz_x_y, t_0_xzz_x_z,\
                           t_0_xzz_y_x, t_0_xzz_y_y, t_0_xzz_y_z, t_0_xzz_z_x, t_0_xzz_z_y,\
                           t_0_xzz_z_z, t_0_yyy_x_x, t_0_yyy_x_y, t_0_yyy_x_z, t_0_yyy_y_x,\
                           t_0_yyy_y_y, t_0_yyy_y_z, t_0_yyy_z_x, t_0_yyy_z_y, t_0_yyy_z_z,\
                           t_0_yyz_x_x, t_0_yyz_x_y, t_0_yyz_x_z, t_0_yyz_y_x, t_0_yyz_y_y,\
                           t_0_yyz_y_z, t_0_yyz_z_x, t_0_yyz_z_y, t_0_yyz_z_z, t_0_yzz_x_x,\
                           t_0_yzz_x_y, t_0_yzz_x_z, t_0_yzz_y_x, t_0_yzz_y_y, t_0_yzz_y_z,\
                           t_0_yzz_z_x, t_0_yzz_z_y, t_0_yzz_z_z, t_x_xzz_x_x, t_x_xzz_x_y,\
                           t_x_xzz_x_z, t_x_xzz_y_x, t_x_xzz_y_y, t_x_xzz_y_z, t_x_xzz_z_x,\
                           t_x_xzz_z_y, t_x_xzz_z_z, t_x_yyy_x_x, t_x_yyy_x_y, t_x_yyy_x_z,\
                           t_x_yyy_y_x, t_x_yyy_y_y, t_x_yyy_y_z, t_x_yyy_z_x, t_x_yyy_z_y,\
                           t_x_yyy_z_z, t_x_yyz_x_x, t_x_yyz_x_y, t_x_yyz_x_z, t_x_yyz_y_x,\
                           t_x_yyz_y_y, t_x_yyz_y_z, t_x_yyz_z_x, t_x_yyz_z_y, t_x_yyz_z_z,\
                           t_x_yzz_x_x, t_x_yzz_x_y, t_x_yzz_x_z, t_x_yzz_y_x, t_x_yzz_y_y,\
                           t_x_yzz_y_z, t_x_yzz_z_x, t_x_yzz_z_y, t_x_yzz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_yzz_z_z[i] = t_0_xyzz_z_z[i] - rab_x[i] * t_0_yzz_z_z[i];

        t_x_yzz_z_y[i] = t_0_xyzz_z_y[i] - rab_x[i] * t_0_yzz_z_y[i];

        t_x_yzz_z_x[i] = t_0_xyzz_z_x[i] - rab_x[i] * t_0_yzz_z_x[i];

        t_x_yzz_y_z[i] = t_0_xyzz_y_z[i] - rab_x[i] * t_0_yzz_y_z[i];

        t_x_yzz_y_y[i] = t_0_xyzz_y_y[i] - rab_x[i] * t_0_yzz_y_y[i];

        t_x_yzz_y_x[i] = t_0_xyzz_y_x[i] - rab_x[i] * t_0_yzz_y_x[i];

        t_x_yzz_x_z[i] = t_0_xyzz_x_z[i] - rab_x[i] * t_0_yzz_x_z[i];

        t_x_yzz_x_y[i] = t_0_xyzz_x_y[i] - rab_x[i] * t_0_yzz_x_y[i];

        t_x_yzz_x_x[i] = t_0_xyzz_x_x[i] - rab_x[i] * t_0_yzz_x_x[i];

        t_x_yyz_z_z[i] = t_0_xyyz_z_z[i] - rab_x[i] * t_0_yyz_z_z[i];

        t_x_yyz_z_y[i] = t_0_xyyz_z_y[i] - rab_x[i] * t_0_yyz_z_y[i];

        t_x_yyz_z_x[i] = t_0_xyyz_z_x[i] - rab_x[i] * t_0_yyz_z_x[i];

        t_x_yyz_y_z[i] = t_0_xyyz_y_z[i] - rab_x[i] * t_0_yyz_y_z[i];

        t_x_yyz_y_y[i] = t_0_xyyz_y_y[i] - rab_x[i] * t_0_yyz_y_y[i];

        t_x_yyz_y_x[i] = t_0_xyyz_y_x[i] - rab_x[i] * t_0_yyz_y_x[i];

        t_x_yyz_x_z[i] = t_0_xyyz_x_z[i] - rab_x[i] * t_0_yyz_x_z[i];

        t_x_yyz_x_y[i] = t_0_xyyz_x_y[i] - rab_x[i] * t_0_yyz_x_y[i];

        t_x_yyz_x_x[i] = t_0_xyyz_x_x[i] - rab_x[i] * t_0_yyz_x_x[i];

        t_x_yyy_z_z[i] = t_0_xyyy_z_z[i] - rab_x[i] * t_0_yyy_z_z[i];

        t_x_yyy_z_y[i] = t_0_xyyy_z_y[i] - rab_x[i] * t_0_yyy_z_y[i];

        t_x_yyy_z_x[i] = t_0_xyyy_z_x[i] - rab_x[i] * t_0_yyy_z_x[i];

        t_x_yyy_y_z[i] = t_0_xyyy_y_z[i] - rab_x[i] * t_0_yyy_y_z[i];

        t_x_yyy_y_y[i] = t_0_xyyy_y_y[i] - rab_x[i] * t_0_yyy_y_y[i];

        t_x_yyy_y_x[i] = t_0_xyyy_y_x[i] - rab_x[i] * t_0_yyy_y_x[i];

        t_x_yyy_x_z[i] = t_0_xyyy_x_z[i] - rab_x[i] * t_0_yyy_x_z[i];

        t_x_yyy_x_y[i] = t_0_xyyy_x_y[i] - rab_x[i] * t_0_yyy_x_y[i];

        t_x_yyy_x_x[i] = t_0_xyyy_x_x[i] - rab_x[i] * t_0_yyy_x_x[i];

        t_x_xzz_z_z[i] = t_0_xxzz_z_z[i] - rab_x[i] * t_0_xzz_z_z[i];

        t_x_xzz_z_y[i] = t_0_xxzz_z_y[i] - rab_x[i] * t_0_xzz_z_y[i];

        t_x_xzz_z_x[i] = t_0_xxzz_z_x[i] - rab_x[i] * t_0_xzz_z_x[i];

        t_x_xzz_y_z[i] = t_0_xxzz_y_z[i] - rab_x[i] * t_0_xzz_y_z[i];

        t_x_xzz_y_y[i] = t_0_xxzz_y_y[i] - rab_x[i] * t_0_xzz_y_y[i];

        t_x_xzz_y_x[i] = t_0_xxzz_y_x[i] - rab_x[i] * t_0_xzz_y_x[i];

        t_x_xzz_x_z[i] = t_0_xxzz_x_z[i] - rab_x[i] * t_0_xzz_x_z[i];

        t_x_xzz_x_y[i] = t_0_xxzz_x_y[i] - rab_x[i] * t_0_xzz_x_y[i];

        t_x_xzz_x_x[i] = t_0_xxzz_x_x[i] - rab_x[i] * t_0_xzz_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xxxy_x_x, t_0_xxxy_x_y, t_0_xxxy_x_z, t_0_xxxy_y_x,\
                           t_0_xxxy_y_y, t_0_xxxy_y_z, t_0_xxxy_z_x, t_0_xxxy_z_y, t_0_xxxy_z_z,\
                           t_0_xxxz_x_x, t_0_xxxz_x_y, t_0_xxxz_x_z, t_0_xxxz_y_x, t_0_xxxz_y_y,\
                           t_0_xxxz_y_z, t_0_xxxz_z_x, t_0_xxxz_z_y, t_0_xxxz_z_z, t_0_xxy_x_x,\
                           t_0_xxy_x_y, t_0_xxy_x_z, t_0_xxy_y_x, t_0_xxy_y_y, t_0_xxy_y_z,\
                           t_0_xxy_z_x, t_0_xxy_z_y, t_0_xxy_z_z, t_0_xxyy_x_x, t_0_xxyy_x_y,\
                           t_0_xxyy_x_z, t_0_xxyy_y_x, t_0_xxyy_y_y, t_0_xxyy_y_z, t_0_xxyy_z_x,\
                           t_0_xxyy_z_y, t_0_xxyy_z_z, t_0_xxyz_x_x, t_0_xxyz_x_y, t_0_xxyz_x_z,\
                           t_0_xxyz_y_x, t_0_xxyz_y_y, t_0_xxyz_y_z, t_0_xxyz_z_x, t_0_xxyz_z_y,\
                           t_0_xxyz_z_z, t_0_xxz_x_x, t_0_xxz_x_y, t_0_xxz_x_z, t_0_xxz_y_x,\
                           t_0_xxz_y_y, t_0_xxz_y_z, t_0_xxz_z_x, t_0_xxz_z_y, t_0_xxz_z_z,\
                           t_0_xyy_x_x, t_0_xyy_x_y, t_0_xyy_x_z, t_0_xyy_y_x, t_0_xyy_y_y,\
                           t_0_xyy_y_z, t_0_xyy_z_x, t_0_xyy_z_y, t_0_xyy_z_z, t_0_xyz_x_x,\
                           t_0_xyz_x_y, t_0_xyz_x_z, t_0_xyz_y_x, t_0_xyz_y_y, t_0_xyz_y_z,\
                           t_0_xyz_z_x, t_0_xyz_z_y, t_0_xyz_z_z, t_x_xxy_x_x, t_x_xxy_x_y,\
                           t_x_xxy_x_z, t_x_xxy_y_x, t_x_xxy_y_y, t_x_xxy_y_z, t_x_xxy_z_x,\
                           t_x_xxy_z_y, t_x_xxy_z_z, t_x_xxz_x_x, t_x_xxz_x_y, t_x_xxz_x_z,\
                           t_x_xxz_y_x, t_x_xxz_y_y, t_x_xxz_y_z, t_x_xxz_z_x, t_x_xxz_z_y,\
                           t_x_xxz_z_z, t_x_xyy_x_x, t_x_xyy_x_y, t_x_xyy_x_z, t_x_xyy_y_x,\
                           t_x_xyy_y_y, t_x_xyy_y_z, t_x_xyy_z_x, t_x_xyy_z_y, t_x_xyy_z_z,\
                           t_x_xyz_x_x, t_x_xyz_x_y, t_x_xyz_x_z, t_x_xyz_y_x, t_x_xyz_y_y,\
                           t_x_xyz_y_z, t_x_xyz_z_x, t_x_xyz_z_y, t_x_xyz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xyz_z_z[i] = t_0_xxyz_z_z[i] - rab_x[i] * t_0_xyz_z_z[i];

        t_x_xyz_z_y[i] = t_0_xxyz_z_y[i] - rab_x[i] * t_0_xyz_z_y[i];

        t_x_xyz_z_x[i] = t_0_xxyz_z_x[i] - rab_x[i] * t_0_xyz_z_x[i];

        t_x_xyz_y_z[i] = t_0_xxyz_y_z[i] - rab_x[i] * t_0_xyz_y_z[i];

        t_x_xyz_y_y[i] = t_0_xxyz_y_y[i] - rab_x[i] * t_0_xyz_y_y[i];

        t_x_xyz_y_x[i] = t_0_xxyz_y_x[i] - rab_x[i] * t_0_xyz_y_x[i];

        t_x_xyz_x_z[i] = t_0_xxyz_x_z[i] - rab_x[i] * t_0_xyz_x_z[i];

        t_x_xyz_x_y[i] = t_0_xxyz_x_y[i] - rab_x[i] * t_0_xyz_x_y[i];

        t_x_xyz_x_x[i] = t_0_xxyz_x_x[i] - rab_x[i] * t_0_xyz_x_x[i];

        t_x_xyy_z_z[i] = t_0_xxyy_z_z[i] - rab_x[i] * t_0_xyy_z_z[i];

        t_x_xyy_z_y[i] = t_0_xxyy_z_y[i] - rab_x[i] * t_0_xyy_z_y[i];

        t_x_xyy_z_x[i] = t_0_xxyy_z_x[i] - rab_x[i] * t_0_xyy_z_x[i];

        t_x_xyy_y_z[i] = t_0_xxyy_y_z[i] - rab_x[i] * t_0_xyy_y_z[i];

        t_x_xyy_y_y[i] = t_0_xxyy_y_y[i] - rab_x[i] * t_0_xyy_y_y[i];

        t_x_xyy_y_x[i] = t_0_xxyy_y_x[i] - rab_x[i] * t_0_xyy_y_x[i];

        t_x_xyy_x_z[i] = t_0_xxyy_x_z[i] - rab_x[i] * t_0_xyy_x_z[i];

        t_x_xyy_x_y[i] = t_0_xxyy_x_y[i] - rab_x[i] * t_0_xyy_x_y[i];

        t_x_xyy_x_x[i] = t_0_xxyy_x_x[i] - rab_x[i] * t_0_xyy_x_x[i];

        t_x_xxz_z_z[i] = t_0_xxxz_z_z[i] - rab_x[i] * t_0_xxz_z_z[i];

        t_x_xxz_z_y[i] = t_0_xxxz_z_y[i] - rab_x[i] * t_0_xxz_z_y[i];

        t_x_xxz_z_x[i] = t_0_xxxz_z_x[i] - rab_x[i] * t_0_xxz_z_x[i];

        t_x_xxz_y_z[i] = t_0_xxxz_y_z[i] - rab_x[i] * t_0_xxz_y_z[i];

        t_x_xxz_y_y[i] = t_0_xxxz_y_y[i] - rab_x[i] * t_0_xxz_y_y[i];

        t_x_xxz_y_x[i] = t_0_xxxz_y_x[i] - rab_x[i] * t_0_xxz_y_x[i];

        t_x_xxz_x_z[i] = t_0_xxxz_x_z[i] - rab_x[i] * t_0_xxz_x_z[i];

        t_x_xxz_x_y[i] = t_0_xxxz_x_y[i] - rab_x[i] * t_0_xxz_x_y[i];

        t_x_xxz_x_x[i] = t_0_xxxz_x_x[i] - rab_x[i] * t_0_xxz_x_x[i];

        t_x_xxy_z_z[i] = t_0_xxxy_z_z[i] - rab_x[i] * t_0_xxy_z_z[i];

        t_x_xxy_z_y[i] = t_0_xxxy_z_y[i] - rab_x[i] * t_0_xxy_z_y[i];

        t_x_xxy_z_x[i] = t_0_xxxy_z_x[i] - rab_x[i] * t_0_xxy_z_x[i];

        t_x_xxy_y_z[i] = t_0_xxxy_y_z[i] - rab_x[i] * t_0_xxy_y_z[i];

        t_x_xxy_y_y[i] = t_0_xxxy_y_y[i] - rab_x[i] * t_0_xxy_y_y[i];

        t_x_xxy_y_x[i] = t_0_xxxy_y_x[i] - rab_x[i] * t_0_xxy_y_x[i];

        t_x_xxy_x_z[i] = t_0_xxxy_x_z[i] - rab_x[i] * t_0_xxy_x_z[i];

        t_x_xxy_x_y[i] = t_0_xxxy_x_y[i] - rab_x[i] * t_0_xxy_x_y[i];

        t_x_xxy_x_x[i] = t_0_xxxy_x_x[i] - rab_x[i] * t_0_xxy_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xxx_x_x, t_0_xxx_x_y, t_0_xxx_x_z, t_0_xxx_y_x,\
                           t_0_xxx_y_y, t_0_xxx_y_z, t_0_xxx_z_x, t_0_xxx_z_y, t_0_xxx_z_z,\
                           t_0_xxxx_x_x, t_0_xxxx_x_y, t_0_xxxx_x_z, t_0_xxxx_y_x, t_0_xxxx_y_y,\
                           t_0_xxxx_y_z, t_0_xxxx_z_x, t_0_xxxx_z_y, t_0_xxxx_z_z, t_x_xxx_x_x,\
                           t_x_xxx_x_y, t_x_xxx_x_z, t_x_xxx_y_x, t_x_xxx_y_y, t_x_xxx_y_z,\
                           t_x_xxx_z_x, t_x_xxx_z_y, t_x_xxx_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxx_z_z[i] = t_0_xxxx_z_z[i] - rab_x[i] * t_0_xxx_z_z[i];

        t_x_xxx_z_y[i] = t_0_xxxx_z_y[i] - rab_x[i] * t_0_xxx_z_y[i];

        t_x_xxx_z_x[i] = t_0_xxxx_z_x[i] - rab_x[i] * t_0_xxx_z_x[i];

        t_x_xxx_y_z[i] = t_0_xxxx_y_z[i] - rab_x[i] * t_0_xxx_y_z[i];

        t_x_xxx_y_y[i] = t_0_xxxx_y_y[i] - rab_x[i] * t_0_xxx_y_y[i];

        t_x_xxx_y_x[i] = t_0_xxxx_y_x[i] - rab_x[i] * t_0_xxx_y_x[i];

        t_x_xxx_x_z[i] = t_0_xxxx_x_z[i] - rab_x[i] * t_0_xxx_x_z[i];

        t_x_xxx_x_y[i] = t_0_xxxx_x_y[i] - rab_x[i] * t_0_xxx_x_y[i];

        t_x_xxx_x_x[i] = t_0_xxxx_x_x[i] - rab_x[i] * t_0_xxx_x_x[i];
    }
}


} // derirec namespace
