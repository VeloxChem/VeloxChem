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
compHostHRRForSFPF_V0(      BufferHostXY<T>&      intsBufferSFPF,
                      const BufferHostX<int32_t>& intsIndexesSFPF,
                      const BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF,
                      const BufferHostXY<T>&      intsBufferSFSG,
                      const BufferHostX<int32_t>& intsIndexesSFSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPF) integral components

    t_0_zzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(0));

    t_0_zzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(1));

    t_0_zzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(2));

    t_0_zzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(3));

    t_0_zzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(4));

    t_0_zzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(5));

    t_0_zzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(6));

    t_0_zzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(7));

    t_0_zzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(8));

    t_0_zzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(9));

    t_0_zzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(10));

    t_0_zzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(11));

    t_0_zzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(12));

    t_0_zzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(13));

    t_0_zzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(14));

    t_0_zzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(15));

    t_0_zzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(16));

    t_0_zzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(17));

    t_0_zzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(18));

    t_0_zzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(19));

    t_0_zzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(20));

    t_0_zzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(21));

    t_0_zzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(22));

    t_0_zzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(23));

    t_0_zzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(24));

    t_0_yzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(25));

    t_0_yzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(26));

    t_0_yzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(27));

    t_0_yzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(28));

    t_0_yzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(29));

    t_0_yzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(30));

    t_0_yzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(31));

    t_0_yzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(32));

    t_0_yzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(33));

    t_0_yzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(34));

    t_0_yzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(35));

    t_0_yzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(36));

    t_0_yzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(37));

    t_0_yzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(38));

    t_0_yzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(39));

    t_0_yzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(40));

    t_0_yzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(41));

    t_0_yzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(42));

    t_0_yzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(43));

    t_0_yzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(44));

    t_0_yzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(45));

    t_0_yzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(46));

    t_0_yzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(47));

    t_0_yzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(48));

    t_0_yzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(49));

    t_0_yyz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(50));

    t_0_yyz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(51));

    t_0_yyz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(52));

    t_0_yyz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(53));

    t_0_yyz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(54));

    t_0_yyz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(55));

    t_0_yyz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(56));

    t_0_yyz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(57));

    t_0_yyz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(58));

    t_0_yyz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(59));

    t_0_yyz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(60));

    t_0_yyz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(61));

    t_0_yyz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(62));

    t_0_yyz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(63));

    t_0_yyz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(64));

    t_0_yyz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(65));

    t_0_yyz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(66));

    t_0_yyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(67));

    t_0_yyz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(68));

    t_0_yyz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(69));

    t_0_yyz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(70));

    t_0_yyz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(71));

    t_0_yyz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(72));

    t_0_yyz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(73));

    t_0_yyz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(74));

    t_0_yyy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(75));

    t_0_yyy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(76));

    t_0_yyy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(77));

    t_0_yyy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(78));

    t_0_yyy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(79));

    t_0_yyy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(80));

    t_0_yyy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(81));

    t_0_yyy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(82));

    t_0_yyy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(83));

    t_0_yyy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(84));

    t_0_yyy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(85));

    t_0_yyy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(86));

    t_0_yyy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(87));

    t_0_yyy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(88));

    t_0_yyy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(89));

    t_0_yyy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(90));

    t_0_yyy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(91));

    t_0_yyy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(92));

    t_0_yyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(93));

    t_0_yyy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(94));

    t_0_yyy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(95));

    t_0_yyy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(96));

    t_0_yyy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(97));

    t_0_yyy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(98));

    t_0_yyy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(99));

    t_0_xzz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(100));

    t_0_xzz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(101));

    t_0_xzz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(102));

    t_0_xzz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(103));

    t_0_xzz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(104));

    t_0_xzz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(105));

    t_0_xzz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(106));

    t_0_xzz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(107));

    t_0_xzz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(108));

    t_0_xzz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(109));

    t_0_xzz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(110));

    t_0_xzz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(111));

    t_0_xzz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(112));

    t_0_xzz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(113));

    t_0_xzz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(114));

    t_0_xzz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(115));

    t_0_xzz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(116));

    t_0_xzz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(117));

    t_0_xzz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(118));

    t_0_xzz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(119));

    t_0_xzz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(120));

    t_0_xzz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(121));

    t_0_xzz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(122));

    t_0_xzz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(123));

    t_0_xzz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(124));

    t_0_xyz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(125));

    t_0_xyz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(126));

    t_0_xyz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(127));

    t_0_xyz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(128));

    t_0_xyz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(129));

    t_0_xyz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(130));

    t_0_xyz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(131));

    t_0_xyz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(132));

    t_0_xyz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(133));

    t_0_xyz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(134));

    t_0_xyz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(135));

    t_0_xyz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(136));

    t_0_xyz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(137));

    t_0_xyz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(138));

    t_0_xyz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(139));

    t_0_xyz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(140));

    t_0_xyz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(141));

    t_0_xyz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(142));

    t_0_xyz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(143));

    t_0_xyz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(144));

    t_0_xyz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(145));

    t_0_xyz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(146));

    t_0_xyz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(147));

    t_0_xyz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(148));

    t_0_xyz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(149));

    t_0_xyy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(150));

    t_0_xyy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(151));

    t_0_xyy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(152));

    t_0_xyy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(153));

    t_0_xyy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(154));

    t_0_xyy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(155));

    t_0_xyy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(156));

    t_0_xyy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(157));

    t_0_xyy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(158));

    t_0_xyy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(159));

    t_0_xyy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(160));

    t_0_xyy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(161));

    t_0_xyy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(162));

    t_0_xyy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(163));

    t_0_xyy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(164));

    t_0_xyy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(165));

    t_0_xyy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(166));

    t_0_xyy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(167));

    t_0_xyy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(168));

    t_0_xyy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(169));

    t_0_xyy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(170));

    t_0_xyy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(171));

    t_0_xyy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(172));

    t_0_xyy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(173));

    t_0_xyy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(174));

    t_0_xxz_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(175));

    t_0_xxz_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(176));

    t_0_xxz_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(177));

    t_0_xxz_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(178));

    t_0_xxz_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(179));

    t_0_xxz_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(180));

    t_0_xxz_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(181));

    t_0_xxz_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(182));

    t_0_xxz_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(183));

    t_0_xxz_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(184));

    t_0_xxz_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(185));

    t_0_xxz_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(186));

    t_0_xxz_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(187));

    t_0_xxz_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(188));

    t_0_xxz_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(189));

    t_0_xxz_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(190));

    t_0_xxz_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(191));

    t_0_xxz_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(192));

    t_0_xxz_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(193));

    t_0_xxz_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(194));

    t_0_xxz_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(195));

    t_0_xxz_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(196));

    t_0_xxz_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(197));

    t_0_xxz_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(198));

    t_0_xxz_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(199));

    t_0_xxy_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(200));

    t_0_xxy_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(201));

    t_0_xxy_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(202));

    t_0_xxy_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(203));

    t_0_xxy_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(204));

    t_0_xxy_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(205));

    t_0_xxy_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(206));

    t_0_xxy_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(207));

    t_0_xxy_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(208));

    t_0_xxy_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(209));

    t_0_xxy_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(210));

    t_0_xxy_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(211));

    t_0_xxy_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(212));

    t_0_xxy_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(213));

    t_0_xxy_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(214));

    t_0_xxy_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(215));

    t_0_xxy_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(216));

    t_0_xxy_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(217));

    t_0_xxy_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(218));

    t_0_xxy_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(219));

    t_0_xxy_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(220));

    t_0_xxy_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(221));

    t_0_xxy_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(222));

    t_0_xxy_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(223));

    t_0_xxy_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(224));

    t_0_xxx_z_zzz = intsBufferSFPF.data(intsIndexesSFPF(225));

    t_0_xxx_z_yzz = intsBufferSFPF.data(intsIndexesSFPF(226));

    t_0_xxx_z_yyz = intsBufferSFPF.data(intsIndexesSFPF(227));

    t_0_xxx_z_xzz = intsBufferSFPF.data(intsIndexesSFPF(228));

    t_0_xxx_z_xyz = intsBufferSFPF.data(intsIndexesSFPF(229));

    t_0_xxx_z_xxz = intsBufferSFPF.data(intsIndexesSFPF(230));

    t_0_xxx_y_zzz = intsBufferSFPF.data(intsIndexesSFPF(231));

    t_0_xxx_y_yzz = intsBufferSFPF.data(intsIndexesSFPF(232));

    t_0_xxx_y_yyz = intsBufferSFPF.data(intsIndexesSFPF(233));

    t_0_xxx_y_yyy = intsBufferSFPF.data(intsIndexesSFPF(234));

    t_0_xxx_y_xzz = intsBufferSFPF.data(intsIndexesSFPF(235));

    t_0_xxx_y_xyz = intsBufferSFPF.data(intsIndexesSFPF(236));

    t_0_xxx_y_xyy = intsBufferSFPF.data(intsIndexesSFPF(237));

    t_0_xxx_y_xxz = intsBufferSFPF.data(intsIndexesSFPF(238));

    t_0_xxx_y_xxy = intsBufferSFPF.data(intsIndexesSFPF(239));

    t_0_xxx_x_zzz = intsBufferSFPF.data(intsIndexesSFPF(240));

    t_0_xxx_x_yzz = intsBufferSFPF.data(intsIndexesSFPF(241));

    t_0_xxx_x_yyz = intsBufferSFPF.data(intsIndexesSFPF(242));

    t_0_xxx_x_yyy = intsBufferSFPF.data(intsIndexesSFPF(243));

    t_0_xxx_x_xzz = intsBufferSFPF.data(intsIndexesSFPF(244));

    t_0_xxx_x_xyz = intsBufferSFPF.data(intsIndexesSFPF(245));

    t_0_xxx_x_xyy = intsBufferSFPF.data(intsIndexesSFPF(246));

    t_0_xxx_x_xxz = intsBufferSFPF.data(intsIndexesSFPF(247));

    t_0_xxx_x_xxy = intsBufferSFPF.data(intsIndexesSFPF(248));

    t_0_xxx_x_xxx = intsBufferSFPF.data(intsIndexesSFPF(249));

    // set up (SFSF) integral components

    t_0_zzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(0));

    t_0_zzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(1));

    t_0_zzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(2));

    t_0_zzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(3));

    t_0_zzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(4));

    t_0_zzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(5));

    t_0_zzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(6));

    t_0_zzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(7));

    t_0_zzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(8));

    t_0_zzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(9));

    t_0_yzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(10));

    t_0_yzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(11));

    t_0_yzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(12));

    t_0_yzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(13));

    t_0_yzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(14));

    t_0_yzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(15));

    t_0_yzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(16));

    t_0_yzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(17));

    t_0_yzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(18));

    t_0_yzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(19));

    t_0_yyz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(20));

    t_0_yyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(21));

    t_0_yyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(22));

    t_0_yyz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(23));

    t_0_yyz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(24));

    t_0_yyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(25));

    t_0_yyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(26));

    t_0_yyz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(27));

    t_0_yyz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(28));

    t_0_yyz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(29));

    t_0_yyy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(30));

    t_0_yyy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(31));

    t_0_yyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(32));

    t_0_yyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(33));

    t_0_yyy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(34));

    t_0_yyy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(35));

    t_0_yyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(36));

    t_0_yyy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(37));

    t_0_yyy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(38));

    t_0_yyy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(39));

    t_0_xzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(40));

    t_0_xzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(41));

    t_0_xzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(42));

    t_0_xzz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(43));

    t_0_xzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(44));

    t_0_xzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(45));

    t_0_xzz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(46));

    t_0_xzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(47));

    t_0_xzz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(48));

    t_0_xzz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(49));

    t_0_xyz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(50));

    t_0_xyz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(51));

    t_0_xyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(52));

    t_0_xyz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(53));

    t_0_xyz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(54));

    t_0_xyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(55));

    t_0_xyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(56));

    t_0_xyz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(57));

    t_0_xyz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(58));

    t_0_xyz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(59));

    t_0_xyy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(60));

    t_0_xyy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(61));

    t_0_xyy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(62));

    t_0_xyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(63));

    t_0_xyy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(64));

    t_0_xyy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(65));

    t_0_xyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(66));

    t_0_xyy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(67));

    t_0_xyy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(68));

    t_0_xyy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(69));

    t_0_xxz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(70));

    t_0_xxz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(71));

    t_0_xxz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(72));

    t_0_xxz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(73));

    t_0_xxz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(74));

    t_0_xxz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(75));

    t_0_xxz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(76));

    t_0_xxz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(77));

    t_0_xxz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(78));

    t_0_xxz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(79));

    t_0_xxy_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(80));

    t_0_xxy_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(81));

    t_0_xxy_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(82));

    t_0_xxy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(83));

    t_0_xxy_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(84));

    t_0_xxy_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(85));

    t_0_xxy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(86));

    t_0_xxy_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(87));

    t_0_xxy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(88));

    t_0_xxy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(89));

    t_0_xxx_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(90));

    t_0_xxx_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(91));

    t_0_xxx_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(92));

    t_0_xxx_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(93));

    t_0_xxx_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(94));

    t_0_xxx_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(95));

    t_0_xxx_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(96));

    t_0_xxx_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(97));

    t_0_xxx_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(98));

    t_0_xxx_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(99));

    // set up (SFSG) integral components

    t_0_zzz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(0));

    t_0_zzz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(1));

    t_0_zzz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(2));

    t_0_zzz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(3));

    t_0_zzz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(4));

    t_0_zzz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(5));

    t_0_zzz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(6));

    t_0_zzz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(7));

    t_0_zzz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(8));

    t_0_zzz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(9));

    t_0_zzz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(10));

    t_0_zzz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(11));

    t_0_zzz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(12));

    t_0_zzz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(13));

    t_0_zzz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(14));

    t_0_yzz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(15));

    t_0_yzz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(16));

    t_0_yzz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(17));

    t_0_yzz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(18));

    t_0_yzz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(19));

    t_0_yzz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(20));

    t_0_yzz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(21));

    t_0_yzz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(22));

    t_0_yzz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(23));

    t_0_yzz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(24));

    t_0_yzz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(25));

    t_0_yzz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(26));

    t_0_yzz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(27));

    t_0_yzz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(28));

    t_0_yzz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(29));

    t_0_yyz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(30));

    t_0_yyz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(31));

    t_0_yyz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(32));

    t_0_yyz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(33));

    t_0_yyz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(34));

    t_0_yyz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(35));

    t_0_yyz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(36));

    t_0_yyz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(37));

    t_0_yyz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(38));

    t_0_yyz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(39));

    t_0_yyz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(40));

    t_0_yyz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(41));

    t_0_yyz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(42));

    t_0_yyz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(43));

    t_0_yyz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(44));

    t_0_yyy_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(45));

    t_0_yyy_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(46));

    t_0_yyy_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(47));

    t_0_yyy_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(48));

    t_0_yyy_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(49));

    t_0_yyy_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(50));

    t_0_yyy_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(51));

    t_0_yyy_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(52));

    t_0_yyy_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(53));

    t_0_yyy_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(54));

    t_0_yyy_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(55));

    t_0_yyy_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(56));

    t_0_yyy_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(57));

    t_0_yyy_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(58));

    t_0_yyy_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(59));

    t_0_xzz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(60));

    t_0_xzz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(61));

    t_0_xzz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(62));

    t_0_xzz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(63));

    t_0_xzz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(64));

    t_0_xzz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(65));

    t_0_xzz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(66));

    t_0_xzz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(67));

    t_0_xzz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(68));

    t_0_xzz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(69));

    t_0_xzz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(70));

    t_0_xzz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(71));

    t_0_xzz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(72));

    t_0_xzz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(73));

    t_0_xzz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(74));

    t_0_xyz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(75));

    t_0_xyz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(76));

    t_0_xyz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(77));

    t_0_xyz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(78));

    t_0_xyz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(79));

    t_0_xyz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(80));

    t_0_xyz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(81));

    t_0_xyz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(82));

    t_0_xyz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(83));

    t_0_xyz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(84));

    t_0_xyz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(85));

    t_0_xyz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(86));

    t_0_xyz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(87));

    t_0_xyz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(88));

    t_0_xyz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(89));

    t_0_xyy_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(90));

    t_0_xyy_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(91));

    t_0_xyy_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(92));

    t_0_xyy_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(93));

    t_0_xyy_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(94));

    t_0_xyy_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(95));

    t_0_xyy_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(96));

    t_0_xyy_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(97));

    t_0_xyy_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(98));

    t_0_xyy_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(99));

    t_0_xyy_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(100));

    t_0_xyy_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(101));

    t_0_xyy_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(102));

    t_0_xyy_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(103));

    t_0_xyy_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(104));

    t_0_xxz_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(105));

    t_0_xxz_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(106));

    t_0_xxz_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(107));

    t_0_xxz_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(108));

    t_0_xxz_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(109));

    t_0_xxz_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(110));

    t_0_xxz_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(111));

    t_0_xxz_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(112));

    t_0_xxz_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(113));

    t_0_xxz_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(114));

    t_0_xxz_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(115));

    t_0_xxz_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(116));

    t_0_xxz_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(117));

    t_0_xxz_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(118));

    t_0_xxz_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(119));

    t_0_xxy_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(120));

    t_0_xxy_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(121));

    t_0_xxy_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(122));

    t_0_xxy_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(123));

    t_0_xxy_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(124));

    t_0_xxy_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(125));

    t_0_xxy_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(126));

    t_0_xxy_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(127));

    t_0_xxy_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(128));

    t_0_xxy_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(129));

    t_0_xxy_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(130));

    t_0_xxy_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(131));

    t_0_xxy_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(132));

    t_0_xxy_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(133));

    t_0_xxy_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(134));

    t_0_xxx_0_zzzz = intsBufferSFSG.data(intsIndexesSFSG(135));

    t_0_xxx_0_yzzz = intsBufferSFSG.data(intsIndexesSFSG(136));

    t_0_xxx_0_yyzz = intsBufferSFSG.data(intsIndexesSFSG(137));

    t_0_xxx_0_yyyz = intsBufferSFSG.data(intsIndexesSFSG(138));

    t_0_xxx_0_yyyy = intsBufferSFSG.data(intsIndexesSFSG(139));

    t_0_xxx_0_xzzz = intsBufferSFSG.data(intsIndexesSFSG(140));

    t_0_xxx_0_xyzz = intsBufferSFSG.data(intsIndexesSFSG(141));

    t_0_xxx_0_xyyz = intsBufferSFSG.data(intsIndexesSFSG(142));

    t_0_xxx_0_xyyy = intsBufferSFSG.data(intsIndexesSFSG(143));

    t_0_xxx_0_xxzz = intsBufferSFSG.data(intsIndexesSFSG(144));

    t_0_xxx_0_xxyz = intsBufferSFSG.data(intsIndexesSFSG(145));

    t_0_xxx_0_xxyy = intsBufferSFSG.data(intsIndexesSFSG(146));

    t_0_xxx_0_xxxz = intsBufferSFSG.data(intsIndexesSFSG(147));

    t_0_xxx_0_xxxy = intsBufferSFSG.data(intsIndexesSFSG(148));

    t_0_xxx_0_xxxx = intsBufferSFSG.data(intsIndexesSFSG(149));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yzz_0_xxz, t_0_yzz_0_xxzz, t_0_yzz_0_xyz,\
                           t_0_yzz_0_xyzz, t_0_yzz_0_xzz, t_0_yzz_0_xzzz, t_0_yzz_0_yyy,\
                           t_0_yzz_0_yyyy, t_0_yzz_0_yyyz, t_0_yzz_0_yyz, t_0_yzz_0_yyzz,\
                           t_0_yzz_0_yzz, t_0_yzz_0_yzzz, t_0_yzz_0_zzz, t_0_yzz_0_zzzz,\
                           t_0_yzz_y_xzz, t_0_yzz_y_yyy, t_0_yzz_y_yyz, t_0_yzz_y_yzz,\
                           t_0_yzz_y_zzz, t_0_yzz_z_xxz, t_0_yzz_z_xyz, t_0_yzz_z_xzz,\
                           t_0_yzz_z_yyz, t_0_yzz_z_yzz, t_0_yzz_z_zzz, t_0_zzz_0_xxx,\
                           t_0_zzz_0_xxxx, t_0_zzz_0_xxxy, t_0_zzz_0_xxxz, t_0_zzz_0_xxy,\
                           t_0_zzz_0_xxyy, t_0_zzz_0_xxyz, t_0_zzz_0_xxz, t_0_zzz_0_xxzz,\
                           t_0_zzz_0_xyy, t_0_zzz_0_xyyy, t_0_zzz_0_xyyz, t_0_zzz_0_xyz,\
                           t_0_zzz_0_xyzz, t_0_zzz_0_xzz, t_0_zzz_0_xzzz, t_0_zzz_0_yyy,\
                           t_0_zzz_0_yyyy, t_0_zzz_0_yyyz, t_0_zzz_0_yyz, t_0_zzz_0_yyzz,\
                           t_0_zzz_0_yzz, t_0_zzz_0_yzzz, t_0_zzz_0_zzz, t_0_zzz_0_zzzz,\
                           t_0_zzz_x_xxx, t_0_zzz_x_xxy, t_0_zzz_x_xxz, t_0_zzz_x_xyy,\
                           t_0_zzz_x_xyz, t_0_zzz_x_xzz, t_0_zzz_x_yyy, t_0_zzz_x_yyz,\
                           t_0_zzz_x_yzz, t_0_zzz_x_zzz, t_0_zzz_y_xxy, t_0_zzz_y_xxz,\
                           t_0_zzz_y_xyy, t_0_zzz_y_xyz, t_0_zzz_y_xzz, t_0_zzz_y_yyy,\
                           t_0_zzz_y_yyz, t_0_zzz_y_yzz, t_0_zzz_y_zzz, t_0_zzz_z_xxz,\
                           t_0_zzz_z_xyz, t_0_zzz_z_xzz, t_0_zzz_z_yyz, t_0_zzz_z_yzz,\
                           t_0_zzz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_zzz[i] = t_0_zzz_0_zzzz[i] - rcd_z[i] * t_0_zzz_0_zzz[i];

        t_0_zzz_z_yzz[i] = t_0_zzz_0_yzzz[i] - rcd_z[i] * t_0_zzz_0_yzz[i];

        t_0_zzz_z_yyz[i] = t_0_zzz_0_yyzz[i] - rcd_z[i] * t_0_zzz_0_yyz[i];

        t_0_zzz_z_xzz[i] = t_0_zzz_0_xzzz[i] - rcd_z[i] * t_0_zzz_0_xzz[i];

        t_0_zzz_z_xyz[i] = t_0_zzz_0_xyzz[i] - rcd_z[i] * t_0_zzz_0_xyz[i];

        t_0_zzz_z_xxz[i] = t_0_zzz_0_xxzz[i] - rcd_z[i] * t_0_zzz_0_xxz[i];

        t_0_zzz_y_zzz[i] = t_0_zzz_0_yzzz[i] - rcd_y[i] * t_0_zzz_0_zzz[i];

        t_0_zzz_y_yzz[i] = t_0_zzz_0_yyzz[i] - rcd_y[i] * t_0_zzz_0_yzz[i];

        t_0_zzz_y_yyz[i] = t_0_zzz_0_yyyz[i] - rcd_y[i] * t_0_zzz_0_yyz[i];

        t_0_zzz_y_yyy[i] = t_0_zzz_0_yyyy[i] - rcd_y[i] * t_0_zzz_0_yyy[i];

        t_0_zzz_y_xzz[i] = t_0_zzz_0_xyzz[i] - rcd_y[i] * t_0_zzz_0_xzz[i];

        t_0_zzz_y_xyz[i] = t_0_zzz_0_xyyz[i] - rcd_y[i] * t_0_zzz_0_xyz[i];

        t_0_zzz_y_xyy[i] = t_0_zzz_0_xyyy[i] - rcd_y[i] * t_0_zzz_0_xyy[i];

        t_0_zzz_y_xxz[i] = t_0_zzz_0_xxyz[i] - rcd_y[i] * t_0_zzz_0_xxz[i];

        t_0_zzz_y_xxy[i] = t_0_zzz_0_xxyy[i] - rcd_y[i] * t_0_zzz_0_xxy[i];

        t_0_zzz_x_zzz[i] = t_0_zzz_0_xzzz[i] - rcd_x[i] * t_0_zzz_0_zzz[i];

        t_0_zzz_x_yzz[i] = t_0_zzz_0_xyzz[i] - rcd_x[i] * t_0_zzz_0_yzz[i];

        t_0_zzz_x_yyz[i] = t_0_zzz_0_xyyz[i] - rcd_x[i] * t_0_zzz_0_yyz[i];

        t_0_zzz_x_yyy[i] = t_0_zzz_0_xyyy[i] - rcd_x[i] * t_0_zzz_0_yyy[i];

        t_0_zzz_x_xzz[i] = t_0_zzz_0_xxzz[i] - rcd_x[i] * t_0_zzz_0_xzz[i];

        t_0_zzz_x_xyz[i] = t_0_zzz_0_xxyz[i] - rcd_x[i] * t_0_zzz_0_xyz[i];

        t_0_zzz_x_xyy[i] = t_0_zzz_0_xxyy[i] - rcd_x[i] * t_0_zzz_0_xyy[i];

        t_0_zzz_x_xxz[i] = t_0_zzz_0_xxxz[i] - rcd_x[i] * t_0_zzz_0_xxz[i];

        t_0_zzz_x_xxy[i] = t_0_zzz_0_xxxy[i] - rcd_x[i] * t_0_zzz_0_xxy[i];

        t_0_zzz_x_xxx[i] = t_0_zzz_0_xxxx[i] - rcd_x[i] * t_0_zzz_0_xxx[i];

        t_0_yzz_z_zzz[i] = t_0_yzz_0_zzzz[i] - rcd_z[i] * t_0_yzz_0_zzz[i];

        t_0_yzz_z_yzz[i] = t_0_yzz_0_yzzz[i] - rcd_z[i] * t_0_yzz_0_yzz[i];

        t_0_yzz_z_yyz[i] = t_0_yzz_0_yyzz[i] - rcd_z[i] * t_0_yzz_0_yyz[i];

        t_0_yzz_z_xzz[i] = t_0_yzz_0_xzzz[i] - rcd_z[i] * t_0_yzz_0_xzz[i];

        t_0_yzz_z_xyz[i] = t_0_yzz_0_xyzz[i] - rcd_z[i] * t_0_yzz_0_xyz[i];

        t_0_yzz_z_xxz[i] = t_0_yzz_0_xxzz[i] - rcd_z[i] * t_0_yzz_0_xxz[i];

        t_0_yzz_y_zzz[i] = t_0_yzz_0_yzzz[i] - rcd_y[i] * t_0_yzz_0_zzz[i];

        t_0_yzz_y_yzz[i] = t_0_yzz_0_yyzz[i] - rcd_y[i] * t_0_yzz_0_yzz[i];

        t_0_yzz_y_yyz[i] = t_0_yzz_0_yyyz[i] - rcd_y[i] * t_0_yzz_0_yyz[i];

        t_0_yzz_y_yyy[i] = t_0_yzz_0_yyyy[i] - rcd_y[i] * t_0_yzz_0_yyy[i];

        t_0_yzz_y_xzz[i] = t_0_yzz_0_xyzz[i] - rcd_y[i] * t_0_yzz_0_xzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyz_0_xxy, t_0_yyz_0_xxyy, t_0_yyz_0_xxyz,\
                           t_0_yyz_0_xxz, t_0_yyz_0_xxzz, t_0_yyz_0_xyy, t_0_yyz_0_xyyy,\
                           t_0_yyz_0_xyyz, t_0_yyz_0_xyz, t_0_yyz_0_xyzz, t_0_yyz_0_xzz,\
                           t_0_yyz_0_xzzz, t_0_yyz_0_yyy, t_0_yyz_0_yyyy, t_0_yyz_0_yyyz,\
                           t_0_yyz_0_yyz, t_0_yyz_0_yyzz, t_0_yyz_0_yzz, t_0_yyz_0_yzzz,\
                           t_0_yyz_0_zzz, t_0_yyz_0_zzzz, t_0_yyz_x_xyy, t_0_yyz_x_xyz,\
                           t_0_yyz_x_xzz, t_0_yyz_x_yyy, t_0_yyz_x_yyz, t_0_yyz_x_yzz,\
                           t_0_yyz_x_zzz, t_0_yyz_y_xxy, t_0_yyz_y_xxz, t_0_yyz_y_xyy,\
                           t_0_yyz_y_xyz, t_0_yyz_y_xzz, t_0_yyz_y_yyy, t_0_yyz_y_yyz,\
                           t_0_yyz_y_yzz, t_0_yyz_y_zzz, t_0_yyz_z_xxz, t_0_yyz_z_xyz,\
                           t_0_yyz_z_xzz, t_0_yyz_z_yyz, t_0_yyz_z_yzz, t_0_yyz_z_zzz,\
                           t_0_yzz_0_xxx, t_0_yzz_0_xxxx, t_0_yzz_0_xxxy, t_0_yzz_0_xxxz,\
                           t_0_yzz_0_xxy, t_0_yzz_0_xxyy, t_0_yzz_0_xxyz, t_0_yzz_0_xxz,\
                           t_0_yzz_0_xxzz, t_0_yzz_0_xyy, t_0_yzz_0_xyyy, t_0_yzz_0_xyyz,\
                           t_0_yzz_0_xyz, t_0_yzz_0_xyzz, t_0_yzz_0_xzz, t_0_yzz_0_xzzz,\
                           t_0_yzz_0_yyy, t_0_yzz_0_yyz, t_0_yzz_0_yzz, t_0_yzz_0_zzz,\
                           t_0_yzz_x_xxx, t_0_yzz_x_xxy, t_0_yzz_x_xxz, t_0_yzz_x_xyy,\
                           t_0_yzz_x_xyz, t_0_yzz_x_xzz, t_0_yzz_x_yyy, t_0_yzz_x_yyz,\
                           t_0_yzz_x_yzz, t_0_yzz_x_zzz, t_0_yzz_y_xxy, t_0_yzz_y_xxz,\
                           t_0_yzz_y_xyy, t_0_yzz_y_xyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yzz_y_xyz[i] = t_0_yzz_0_xyyz[i] - rcd_y[i] * t_0_yzz_0_xyz[i];

        t_0_yzz_y_xyy[i] = t_0_yzz_0_xyyy[i] - rcd_y[i] * t_0_yzz_0_xyy[i];

        t_0_yzz_y_xxz[i] = t_0_yzz_0_xxyz[i] - rcd_y[i] * t_0_yzz_0_xxz[i];

        t_0_yzz_y_xxy[i] = t_0_yzz_0_xxyy[i] - rcd_y[i] * t_0_yzz_0_xxy[i];

        t_0_yzz_x_zzz[i] = t_0_yzz_0_xzzz[i] - rcd_x[i] * t_0_yzz_0_zzz[i];

        t_0_yzz_x_yzz[i] = t_0_yzz_0_xyzz[i] - rcd_x[i] * t_0_yzz_0_yzz[i];

        t_0_yzz_x_yyz[i] = t_0_yzz_0_xyyz[i] - rcd_x[i] * t_0_yzz_0_yyz[i];

        t_0_yzz_x_yyy[i] = t_0_yzz_0_xyyy[i] - rcd_x[i] * t_0_yzz_0_yyy[i];

        t_0_yzz_x_xzz[i] = t_0_yzz_0_xxzz[i] - rcd_x[i] * t_0_yzz_0_xzz[i];

        t_0_yzz_x_xyz[i] = t_0_yzz_0_xxyz[i] - rcd_x[i] * t_0_yzz_0_xyz[i];

        t_0_yzz_x_xyy[i] = t_0_yzz_0_xxyy[i] - rcd_x[i] * t_0_yzz_0_xyy[i];

        t_0_yzz_x_xxz[i] = t_0_yzz_0_xxxz[i] - rcd_x[i] * t_0_yzz_0_xxz[i];

        t_0_yzz_x_xxy[i] = t_0_yzz_0_xxxy[i] - rcd_x[i] * t_0_yzz_0_xxy[i];

        t_0_yzz_x_xxx[i] = t_0_yzz_0_xxxx[i] - rcd_x[i] * t_0_yzz_0_xxx[i];

        t_0_yyz_z_zzz[i] = t_0_yyz_0_zzzz[i] - rcd_z[i] * t_0_yyz_0_zzz[i];

        t_0_yyz_z_yzz[i] = t_0_yyz_0_yzzz[i] - rcd_z[i] * t_0_yyz_0_yzz[i];

        t_0_yyz_z_yyz[i] = t_0_yyz_0_yyzz[i] - rcd_z[i] * t_0_yyz_0_yyz[i];

        t_0_yyz_z_xzz[i] = t_0_yyz_0_xzzz[i] - rcd_z[i] * t_0_yyz_0_xzz[i];

        t_0_yyz_z_xyz[i] = t_0_yyz_0_xyzz[i] - rcd_z[i] * t_0_yyz_0_xyz[i];

        t_0_yyz_z_xxz[i] = t_0_yyz_0_xxzz[i] - rcd_z[i] * t_0_yyz_0_xxz[i];

        t_0_yyz_y_zzz[i] = t_0_yyz_0_yzzz[i] - rcd_y[i] * t_0_yyz_0_zzz[i];

        t_0_yyz_y_yzz[i] = t_0_yyz_0_yyzz[i] - rcd_y[i] * t_0_yyz_0_yzz[i];

        t_0_yyz_y_yyz[i] = t_0_yyz_0_yyyz[i] - rcd_y[i] * t_0_yyz_0_yyz[i];

        t_0_yyz_y_yyy[i] = t_0_yyz_0_yyyy[i] - rcd_y[i] * t_0_yyz_0_yyy[i];

        t_0_yyz_y_xzz[i] = t_0_yyz_0_xyzz[i] - rcd_y[i] * t_0_yyz_0_xzz[i];

        t_0_yyz_y_xyz[i] = t_0_yyz_0_xyyz[i] - rcd_y[i] * t_0_yyz_0_xyz[i];

        t_0_yyz_y_xyy[i] = t_0_yyz_0_xyyy[i] - rcd_y[i] * t_0_yyz_0_xyy[i];

        t_0_yyz_y_xxz[i] = t_0_yyz_0_xxyz[i] - rcd_y[i] * t_0_yyz_0_xxz[i];

        t_0_yyz_y_xxy[i] = t_0_yyz_0_xxyy[i] - rcd_y[i] * t_0_yyz_0_xxy[i];

        t_0_yyz_x_zzz[i] = t_0_yyz_0_xzzz[i] - rcd_x[i] * t_0_yyz_0_zzz[i];

        t_0_yyz_x_yzz[i] = t_0_yyz_0_xyzz[i] - rcd_x[i] * t_0_yyz_0_yzz[i];

        t_0_yyz_x_yyz[i] = t_0_yyz_0_xyyz[i] - rcd_x[i] * t_0_yyz_0_yyz[i];

        t_0_yyz_x_yyy[i] = t_0_yyz_0_xyyy[i] - rcd_x[i] * t_0_yyz_0_yyy[i];

        t_0_yyz_x_xzz[i] = t_0_yyz_0_xxzz[i] - rcd_x[i] * t_0_yyz_0_xzz[i];

        t_0_yyz_x_xyz[i] = t_0_yyz_0_xxyz[i] - rcd_x[i] * t_0_yyz_0_xyz[i];

        t_0_yyz_x_xyy[i] = t_0_yyz_0_xxyy[i] - rcd_x[i] * t_0_yyz_0_xyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xzz_0_xxz, t_0_xzz_0_xxzz, t_0_xzz_0_xyz,\
                           t_0_xzz_0_xyzz, t_0_xzz_0_xzz, t_0_xzz_0_xzzz, t_0_xzz_0_yyz,\
                           t_0_xzz_0_yyzz, t_0_xzz_0_yzz, t_0_xzz_0_yzzz, t_0_xzz_0_zzz,\
                           t_0_xzz_0_zzzz, t_0_xzz_y_yzz, t_0_xzz_y_zzz, t_0_xzz_z_xxz,\
                           t_0_xzz_z_xyz, t_0_xzz_z_xzz, t_0_xzz_z_yyz, t_0_xzz_z_yzz,\
                           t_0_xzz_z_zzz, t_0_yyy_0_xxx, t_0_yyy_0_xxxx, t_0_yyy_0_xxxy,\
                           t_0_yyy_0_xxxz, t_0_yyy_0_xxy, t_0_yyy_0_xxyy, t_0_yyy_0_xxyz,\
                           t_0_yyy_0_xxz, t_0_yyy_0_xxzz, t_0_yyy_0_xyy, t_0_yyy_0_xyyy,\
                           t_0_yyy_0_xyyz, t_0_yyy_0_xyz, t_0_yyy_0_xyzz, t_0_yyy_0_xzz,\
                           t_0_yyy_0_xzzz, t_0_yyy_0_yyy, t_0_yyy_0_yyyy, t_0_yyy_0_yyyz,\
                           t_0_yyy_0_yyz, t_0_yyy_0_yyzz, t_0_yyy_0_yzz, t_0_yyy_0_yzzz,\
                           t_0_yyy_0_zzz, t_0_yyy_0_zzzz, t_0_yyy_x_xxx, t_0_yyy_x_xxy,\
                           t_0_yyy_x_xxz, t_0_yyy_x_xyy, t_0_yyy_x_xyz, t_0_yyy_x_xzz,\
                           t_0_yyy_x_yyy, t_0_yyy_x_yyz, t_0_yyy_x_yzz, t_0_yyy_x_zzz,\
                           t_0_yyy_y_xxy, t_0_yyy_y_xxz, t_0_yyy_y_xyy, t_0_yyy_y_xyz,\
                           t_0_yyy_y_xzz, t_0_yyy_y_yyy, t_0_yyy_y_yyz, t_0_yyy_y_yzz,\
                           t_0_yyy_y_zzz, t_0_yyy_z_xxz, t_0_yyy_z_xyz, t_0_yyy_z_xzz,\
                           t_0_yyy_z_yyz, t_0_yyy_z_yzz, t_0_yyy_z_zzz, t_0_yyz_0_xxx,\
                           t_0_yyz_0_xxxx, t_0_yyz_0_xxxy, t_0_yyz_0_xxxz, t_0_yyz_0_xxy,\
                           t_0_yyz_0_xxz, t_0_yyz_x_xxx, t_0_yyz_x_xxy, t_0_yyz_x_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyz_x_xxz[i] = t_0_yyz_0_xxxz[i] - rcd_x[i] * t_0_yyz_0_xxz[i];

        t_0_yyz_x_xxy[i] = t_0_yyz_0_xxxy[i] - rcd_x[i] * t_0_yyz_0_xxy[i];

        t_0_yyz_x_xxx[i] = t_0_yyz_0_xxxx[i] - rcd_x[i] * t_0_yyz_0_xxx[i];

        t_0_yyy_z_zzz[i] = t_0_yyy_0_zzzz[i] - rcd_z[i] * t_0_yyy_0_zzz[i];

        t_0_yyy_z_yzz[i] = t_0_yyy_0_yzzz[i] - rcd_z[i] * t_0_yyy_0_yzz[i];

        t_0_yyy_z_yyz[i] = t_0_yyy_0_yyzz[i] - rcd_z[i] * t_0_yyy_0_yyz[i];

        t_0_yyy_z_xzz[i] = t_0_yyy_0_xzzz[i] - rcd_z[i] * t_0_yyy_0_xzz[i];

        t_0_yyy_z_xyz[i] = t_0_yyy_0_xyzz[i] - rcd_z[i] * t_0_yyy_0_xyz[i];

        t_0_yyy_z_xxz[i] = t_0_yyy_0_xxzz[i] - rcd_z[i] * t_0_yyy_0_xxz[i];

        t_0_yyy_y_zzz[i] = t_0_yyy_0_yzzz[i] - rcd_y[i] * t_0_yyy_0_zzz[i];

        t_0_yyy_y_yzz[i] = t_0_yyy_0_yyzz[i] - rcd_y[i] * t_0_yyy_0_yzz[i];

        t_0_yyy_y_yyz[i] = t_0_yyy_0_yyyz[i] - rcd_y[i] * t_0_yyy_0_yyz[i];

        t_0_yyy_y_yyy[i] = t_0_yyy_0_yyyy[i] - rcd_y[i] * t_0_yyy_0_yyy[i];

        t_0_yyy_y_xzz[i] = t_0_yyy_0_xyzz[i] - rcd_y[i] * t_0_yyy_0_xzz[i];

        t_0_yyy_y_xyz[i] = t_0_yyy_0_xyyz[i] - rcd_y[i] * t_0_yyy_0_xyz[i];

        t_0_yyy_y_xyy[i] = t_0_yyy_0_xyyy[i] - rcd_y[i] * t_0_yyy_0_xyy[i];

        t_0_yyy_y_xxz[i] = t_0_yyy_0_xxyz[i] - rcd_y[i] * t_0_yyy_0_xxz[i];

        t_0_yyy_y_xxy[i] = t_0_yyy_0_xxyy[i] - rcd_y[i] * t_0_yyy_0_xxy[i];

        t_0_yyy_x_zzz[i] = t_0_yyy_0_xzzz[i] - rcd_x[i] * t_0_yyy_0_zzz[i];

        t_0_yyy_x_yzz[i] = t_0_yyy_0_xyzz[i] - rcd_x[i] * t_0_yyy_0_yzz[i];

        t_0_yyy_x_yyz[i] = t_0_yyy_0_xyyz[i] - rcd_x[i] * t_0_yyy_0_yyz[i];

        t_0_yyy_x_yyy[i] = t_0_yyy_0_xyyy[i] - rcd_x[i] * t_0_yyy_0_yyy[i];

        t_0_yyy_x_xzz[i] = t_0_yyy_0_xxzz[i] - rcd_x[i] * t_0_yyy_0_xzz[i];

        t_0_yyy_x_xyz[i] = t_0_yyy_0_xxyz[i] - rcd_x[i] * t_0_yyy_0_xyz[i];

        t_0_yyy_x_xyy[i] = t_0_yyy_0_xxyy[i] - rcd_x[i] * t_0_yyy_0_xyy[i];

        t_0_yyy_x_xxz[i] = t_0_yyy_0_xxxz[i] - rcd_x[i] * t_0_yyy_0_xxz[i];

        t_0_yyy_x_xxy[i] = t_0_yyy_0_xxxy[i] - rcd_x[i] * t_0_yyy_0_xxy[i];

        t_0_yyy_x_xxx[i] = t_0_yyy_0_xxxx[i] - rcd_x[i] * t_0_yyy_0_xxx[i];

        t_0_xzz_z_zzz[i] = t_0_xzz_0_zzzz[i] - rcd_z[i] * t_0_xzz_0_zzz[i];

        t_0_xzz_z_yzz[i] = t_0_xzz_0_yzzz[i] - rcd_z[i] * t_0_xzz_0_yzz[i];

        t_0_xzz_z_yyz[i] = t_0_xzz_0_yyzz[i] - rcd_z[i] * t_0_xzz_0_yyz[i];

        t_0_xzz_z_xzz[i] = t_0_xzz_0_xzzz[i] - rcd_z[i] * t_0_xzz_0_xzz[i];

        t_0_xzz_z_xyz[i] = t_0_xzz_0_xyzz[i] - rcd_z[i] * t_0_xzz_0_xyz[i];

        t_0_xzz_z_xxz[i] = t_0_xzz_0_xxzz[i] - rcd_z[i] * t_0_xzz_0_xxz[i];

        t_0_xzz_y_zzz[i] = t_0_xzz_0_yzzz[i] - rcd_y[i] * t_0_xzz_0_zzz[i];

        t_0_xzz_y_yzz[i] = t_0_xzz_0_yyzz[i] - rcd_y[i] * t_0_xzz_0_yzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyz_0_xxy, t_0_xyz_0_xxyy, t_0_xyz_0_xxyz,\
                           t_0_xyz_0_xxz, t_0_xyz_0_xxzz, t_0_xyz_0_xyy, t_0_xyz_0_xyyy,\
                           t_0_xyz_0_xyyz, t_0_xyz_0_xyz, t_0_xyz_0_xyzz, t_0_xyz_0_xzz,\
                           t_0_xyz_0_xzzz, t_0_xyz_0_yyy, t_0_xyz_0_yyyy, t_0_xyz_0_yyyz,\
                           t_0_xyz_0_yyz, t_0_xyz_0_yyzz, t_0_xyz_0_yzz, t_0_xyz_0_yzzz,\
                           t_0_xyz_0_zzz, t_0_xyz_0_zzzz, t_0_xyz_x_yyy, t_0_xyz_x_yyz,\
                           t_0_xyz_x_yzz, t_0_xyz_x_zzz, t_0_xyz_y_xxy, t_0_xyz_y_xxz,\
                           t_0_xyz_y_xyy, t_0_xyz_y_xyz, t_0_xyz_y_xzz, t_0_xyz_y_yyy,\
                           t_0_xyz_y_yyz, t_0_xyz_y_yzz, t_0_xyz_y_zzz, t_0_xyz_z_xxz,\
                           t_0_xyz_z_xyz, t_0_xyz_z_xzz, t_0_xyz_z_yyz, t_0_xyz_z_yzz,\
                           t_0_xyz_z_zzz, t_0_xzz_0_xxx, t_0_xzz_0_xxxx, t_0_xzz_0_xxxy,\
                           t_0_xzz_0_xxxz, t_0_xzz_0_xxy, t_0_xzz_0_xxyy, t_0_xzz_0_xxyz,\
                           t_0_xzz_0_xxz, t_0_xzz_0_xxzz, t_0_xzz_0_xyy, t_0_xzz_0_xyyy,\
                           t_0_xzz_0_xyyz, t_0_xzz_0_xyz, t_0_xzz_0_xyzz, t_0_xzz_0_xzz,\
                           t_0_xzz_0_xzzz, t_0_xzz_0_yyy, t_0_xzz_0_yyyy, t_0_xzz_0_yyyz,\
                           t_0_xzz_0_yyz, t_0_xzz_0_yzz, t_0_xzz_0_zzz, t_0_xzz_x_xxx,\
                           t_0_xzz_x_xxy, t_0_xzz_x_xxz, t_0_xzz_x_xyy, t_0_xzz_x_xyz,\
                           t_0_xzz_x_xzz, t_0_xzz_x_yyy, t_0_xzz_x_yyz, t_0_xzz_x_yzz,\
                           t_0_xzz_x_zzz, t_0_xzz_y_xxy, t_0_xzz_y_xxz, t_0_xzz_y_xyy,\
                           t_0_xzz_y_xyz, t_0_xzz_y_xzz, t_0_xzz_y_yyy, t_0_xzz_y_yyz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xzz_y_yyz[i] = t_0_xzz_0_yyyz[i] - rcd_y[i] * t_0_xzz_0_yyz[i];

        t_0_xzz_y_yyy[i] = t_0_xzz_0_yyyy[i] - rcd_y[i] * t_0_xzz_0_yyy[i];

        t_0_xzz_y_xzz[i] = t_0_xzz_0_xyzz[i] - rcd_y[i] * t_0_xzz_0_xzz[i];

        t_0_xzz_y_xyz[i] = t_0_xzz_0_xyyz[i] - rcd_y[i] * t_0_xzz_0_xyz[i];

        t_0_xzz_y_xyy[i] = t_0_xzz_0_xyyy[i] - rcd_y[i] * t_0_xzz_0_xyy[i];

        t_0_xzz_y_xxz[i] = t_0_xzz_0_xxyz[i] - rcd_y[i] * t_0_xzz_0_xxz[i];

        t_0_xzz_y_xxy[i] = t_0_xzz_0_xxyy[i] - rcd_y[i] * t_0_xzz_0_xxy[i];

        t_0_xzz_x_zzz[i] = t_0_xzz_0_xzzz[i] - rcd_x[i] * t_0_xzz_0_zzz[i];

        t_0_xzz_x_yzz[i] = t_0_xzz_0_xyzz[i] - rcd_x[i] * t_0_xzz_0_yzz[i];

        t_0_xzz_x_yyz[i] = t_0_xzz_0_xyyz[i] - rcd_x[i] * t_0_xzz_0_yyz[i];

        t_0_xzz_x_yyy[i] = t_0_xzz_0_xyyy[i] - rcd_x[i] * t_0_xzz_0_yyy[i];

        t_0_xzz_x_xzz[i] = t_0_xzz_0_xxzz[i] - rcd_x[i] * t_0_xzz_0_xzz[i];

        t_0_xzz_x_xyz[i] = t_0_xzz_0_xxyz[i] - rcd_x[i] * t_0_xzz_0_xyz[i];

        t_0_xzz_x_xyy[i] = t_0_xzz_0_xxyy[i] - rcd_x[i] * t_0_xzz_0_xyy[i];

        t_0_xzz_x_xxz[i] = t_0_xzz_0_xxxz[i] - rcd_x[i] * t_0_xzz_0_xxz[i];

        t_0_xzz_x_xxy[i] = t_0_xzz_0_xxxy[i] - rcd_x[i] * t_0_xzz_0_xxy[i];

        t_0_xzz_x_xxx[i] = t_0_xzz_0_xxxx[i] - rcd_x[i] * t_0_xzz_0_xxx[i];

        t_0_xyz_z_zzz[i] = t_0_xyz_0_zzzz[i] - rcd_z[i] * t_0_xyz_0_zzz[i];

        t_0_xyz_z_yzz[i] = t_0_xyz_0_yzzz[i] - rcd_z[i] * t_0_xyz_0_yzz[i];

        t_0_xyz_z_yyz[i] = t_0_xyz_0_yyzz[i] - rcd_z[i] * t_0_xyz_0_yyz[i];

        t_0_xyz_z_xzz[i] = t_0_xyz_0_xzzz[i] - rcd_z[i] * t_0_xyz_0_xzz[i];

        t_0_xyz_z_xyz[i] = t_0_xyz_0_xyzz[i] - rcd_z[i] * t_0_xyz_0_xyz[i];

        t_0_xyz_z_xxz[i] = t_0_xyz_0_xxzz[i] - rcd_z[i] * t_0_xyz_0_xxz[i];

        t_0_xyz_y_zzz[i] = t_0_xyz_0_yzzz[i] - rcd_y[i] * t_0_xyz_0_zzz[i];

        t_0_xyz_y_yzz[i] = t_0_xyz_0_yyzz[i] - rcd_y[i] * t_0_xyz_0_yzz[i];

        t_0_xyz_y_yyz[i] = t_0_xyz_0_yyyz[i] - rcd_y[i] * t_0_xyz_0_yyz[i];

        t_0_xyz_y_yyy[i] = t_0_xyz_0_yyyy[i] - rcd_y[i] * t_0_xyz_0_yyy[i];

        t_0_xyz_y_xzz[i] = t_0_xyz_0_xyzz[i] - rcd_y[i] * t_0_xyz_0_xzz[i];

        t_0_xyz_y_xyz[i] = t_0_xyz_0_xyyz[i] - rcd_y[i] * t_0_xyz_0_xyz[i];

        t_0_xyz_y_xyy[i] = t_0_xyz_0_xyyy[i] - rcd_y[i] * t_0_xyz_0_xyy[i];

        t_0_xyz_y_xxz[i] = t_0_xyz_0_xxyz[i] - rcd_y[i] * t_0_xyz_0_xxz[i];

        t_0_xyz_y_xxy[i] = t_0_xyz_0_xxyy[i] - rcd_y[i] * t_0_xyz_0_xxy[i];

        t_0_xyz_x_zzz[i] = t_0_xyz_0_xzzz[i] - rcd_x[i] * t_0_xyz_0_zzz[i];

        t_0_xyz_x_yzz[i] = t_0_xyz_0_xyzz[i] - rcd_x[i] * t_0_xyz_0_yzz[i];

        t_0_xyz_x_yyz[i] = t_0_xyz_0_xyyz[i] - rcd_x[i] * t_0_xyz_0_yyz[i];

        t_0_xyz_x_yyy[i] = t_0_xyz_0_xyyy[i] - rcd_x[i] * t_0_xyz_0_yyy[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxz_0_xyz, t_0_xxz_0_xyzz, t_0_xxz_0_xzz,\
                           t_0_xxz_0_xzzz, t_0_xxz_0_yyz, t_0_xxz_0_yyzz, t_0_xxz_0_yzz,\
                           t_0_xxz_0_yzzz, t_0_xxz_0_zzz, t_0_xxz_0_zzzz, t_0_xxz_z_xyz,\
                           t_0_xxz_z_xzz, t_0_xxz_z_yyz, t_0_xxz_z_yzz, t_0_xxz_z_zzz,\
                           t_0_xyy_0_xxx, t_0_xyy_0_xxxx, t_0_xyy_0_xxxy, t_0_xyy_0_xxxz,\
                           t_0_xyy_0_xxy, t_0_xyy_0_xxyy, t_0_xyy_0_xxyz, t_0_xyy_0_xxz,\
                           t_0_xyy_0_xxzz, t_0_xyy_0_xyy, t_0_xyy_0_xyyy, t_0_xyy_0_xyyz,\
                           t_0_xyy_0_xyz, t_0_xyy_0_xyzz, t_0_xyy_0_xzz, t_0_xyy_0_xzzz,\
                           t_0_xyy_0_yyy, t_0_xyy_0_yyyy, t_0_xyy_0_yyyz, t_0_xyy_0_yyz,\
                           t_0_xyy_0_yyzz, t_0_xyy_0_yzz, t_0_xyy_0_yzzz, t_0_xyy_0_zzz,\
                           t_0_xyy_0_zzzz, t_0_xyy_x_xxx, t_0_xyy_x_xxy, t_0_xyy_x_xxz,\
                           t_0_xyy_x_xyy, t_0_xyy_x_xyz, t_0_xyy_x_xzz, t_0_xyy_x_yyy,\
                           t_0_xyy_x_yyz, t_0_xyy_x_yzz, t_0_xyy_x_zzz, t_0_xyy_y_xxy,\
                           t_0_xyy_y_xxz, t_0_xyy_y_xyy, t_0_xyy_y_xyz, t_0_xyy_y_xzz,\
                           t_0_xyy_y_yyy, t_0_xyy_y_yyz, t_0_xyy_y_yzz, t_0_xyy_y_zzz,\
                           t_0_xyy_z_xxz, t_0_xyy_z_xyz, t_0_xyy_z_xzz, t_0_xyy_z_yyz,\
                           t_0_xyy_z_yzz, t_0_xyy_z_zzz, t_0_xyz_0_xxx, t_0_xyz_0_xxxx,\
                           t_0_xyz_0_xxxy, t_0_xyz_0_xxxz, t_0_xyz_0_xxy, t_0_xyz_0_xxyy,\
                           t_0_xyz_0_xxyz, t_0_xyz_0_xxz, t_0_xyz_0_xxzz, t_0_xyz_0_xyy,\
                           t_0_xyz_0_xyz, t_0_xyz_0_xzz, t_0_xyz_x_xxx, t_0_xyz_x_xxy,\
                           t_0_xyz_x_xxz, t_0_xyz_x_xyy, t_0_xyz_x_xyz, t_0_xyz_x_xzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyz_x_xzz[i] = t_0_xyz_0_xxzz[i] - rcd_x[i] * t_0_xyz_0_xzz[i];

        t_0_xyz_x_xyz[i] = t_0_xyz_0_xxyz[i] - rcd_x[i] * t_0_xyz_0_xyz[i];

        t_0_xyz_x_xyy[i] = t_0_xyz_0_xxyy[i] - rcd_x[i] * t_0_xyz_0_xyy[i];

        t_0_xyz_x_xxz[i] = t_0_xyz_0_xxxz[i] - rcd_x[i] * t_0_xyz_0_xxz[i];

        t_0_xyz_x_xxy[i] = t_0_xyz_0_xxxy[i] - rcd_x[i] * t_0_xyz_0_xxy[i];

        t_0_xyz_x_xxx[i] = t_0_xyz_0_xxxx[i] - rcd_x[i] * t_0_xyz_0_xxx[i];

        t_0_xyy_z_zzz[i] = t_0_xyy_0_zzzz[i] - rcd_z[i] * t_0_xyy_0_zzz[i];

        t_0_xyy_z_yzz[i] = t_0_xyy_0_yzzz[i] - rcd_z[i] * t_0_xyy_0_yzz[i];

        t_0_xyy_z_yyz[i] = t_0_xyy_0_yyzz[i] - rcd_z[i] * t_0_xyy_0_yyz[i];

        t_0_xyy_z_xzz[i] = t_0_xyy_0_xzzz[i] - rcd_z[i] * t_0_xyy_0_xzz[i];

        t_0_xyy_z_xyz[i] = t_0_xyy_0_xyzz[i] - rcd_z[i] * t_0_xyy_0_xyz[i];

        t_0_xyy_z_xxz[i] = t_0_xyy_0_xxzz[i] - rcd_z[i] * t_0_xyy_0_xxz[i];

        t_0_xyy_y_zzz[i] = t_0_xyy_0_yzzz[i] - rcd_y[i] * t_0_xyy_0_zzz[i];

        t_0_xyy_y_yzz[i] = t_0_xyy_0_yyzz[i] - rcd_y[i] * t_0_xyy_0_yzz[i];

        t_0_xyy_y_yyz[i] = t_0_xyy_0_yyyz[i] - rcd_y[i] * t_0_xyy_0_yyz[i];

        t_0_xyy_y_yyy[i] = t_0_xyy_0_yyyy[i] - rcd_y[i] * t_0_xyy_0_yyy[i];

        t_0_xyy_y_xzz[i] = t_0_xyy_0_xyzz[i] - rcd_y[i] * t_0_xyy_0_xzz[i];

        t_0_xyy_y_xyz[i] = t_0_xyy_0_xyyz[i] - rcd_y[i] * t_0_xyy_0_xyz[i];

        t_0_xyy_y_xyy[i] = t_0_xyy_0_xyyy[i] - rcd_y[i] * t_0_xyy_0_xyy[i];

        t_0_xyy_y_xxz[i] = t_0_xyy_0_xxyz[i] - rcd_y[i] * t_0_xyy_0_xxz[i];

        t_0_xyy_y_xxy[i] = t_0_xyy_0_xxyy[i] - rcd_y[i] * t_0_xyy_0_xxy[i];

        t_0_xyy_x_zzz[i] = t_0_xyy_0_xzzz[i] - rcd_x[i] * t_0_xyy_0_zzz[i];

        t_0_xyy_x_yzz[i] = t_0_xyy_0_xyzz[i] - rcd_x[i] * t_0_xyy_0_yzz[i];

        t_0_xyy_x_yyz[i] = t_0_xyy_0_xyyz[i] - rcd_x[i] * t_0_xyy_0_yyz[i];

        t_0_xyy_x_yyy[i] = t_0_xyy_0_xyyy[i] - rcd_x[i] * t_0_xyy_0_yyy[i];

        t_0_xyy_x_xzz[i] = t_0_xyy_0_xxzz[i] - rcd_x[i] * t_0_xyy_0_xzz[i];

        t_0_xyy_x_xyz[i] = t_0_xyy_0_xxyz[i] - rcd_x[i] * t_0_xyy_0_xyz[i];

        t_0_xyy_x_xyy[i] = t_0_xyy_0_xxyy[i] - rcd_x[i] * t_0_xyy_0_xyy[i];

        t_0_xyy_x_xxz[i] = t_0_xyy_0_xxxz[i] - rcd_x[i] * t_0_xyy_0_xxz[i];

        t_0_xyy_x_xxy[i] = t_0_xyy_0_xxxy[i] - rcd_x[i] * t_0_xyy_0_xxy[i];

        t_0_xyy_x_xxx[i] = t_0_xyy_0_xxxx[i] - rcd_x[i] * t_0_xyy_0_xxx[i];

        t_0_xxz_z_zzz[i] = t_0_xxz_0_zzzz[i] - rcd_z[i] * t_0_xxz_0_zzz[i];

        t_0_xxz_z_yzz[i] = t_0_xxz_0_yzzz[i] - rcd_z[i] * t_0_xxz_0_yzz[i];

        t_0_xxz_z_yyz[i] = t_0_xxz_0_yyzz[i] - rcd_z[i] * t_0_xxz_0_yyz[i];

        t_0_xxz_z_xzz[i] = t_0_xxz_0_xzzz[i] - rcd_z[i] * t_0_xxz_0_xzz[i];

        t_0_xxz_z_xyz[i] = t_0_xxz_0_xyzz[i] - rcd_z[i] * t_0_xxz_0_xyz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxy_0_xxy, t_0_xxy_0_xxyy, t_0_xxy_0_xxyz,\
                           t_0_xxy_0_xxz, t_0_xxy_0_xxzz, t_0_xxy_0_xyy, t_0_xxy_0_xyyy,\
                           t_0_xxy_0_xyyz, t_0_xxy_0_xyz, t_0_xxy_0_xyzz, t_0_xxy_0_xzz,\
                           t_0_xxy_0_xzzz, t_0_xxy_0_yyy, t_0_xxy_0_yyyy, t_0_xxy_0_yyyz,\
                           t_0_xxy_0_yyz, t_0_xxy_0_yyzz, t_0_xxy_0_yzz, t_0_xxy_0_yzzz,\
                           t_0_xxy_0_zzz, t_0_xxy_0_zzzz, t_0_xxy_x_zzz, t_0_xxy_y_xxy,\
                           t_0_xxy_y_xxz, t_0_xxy_y_xyy, t_0_xxy_y_xyz, t_0_xxy_y_xzz,\
                           t_0_xxy_y_yyy, t_0_xxy_y_yyz, t_0_xxy_y_yzz, t_0_xxy_y_zzz,\
                           t_0_xxy_z_xxz, t_0_xxy_z_xyz, t_0_xxy_z_xzz, t_0_xxy_z_yyz,\
                           t_0_xxy_z_yzz, t_0_xxy_z_zzz, t_0_xxz_0_xxx, t_0_xxz_0_xxxx,\
                           t_0_xxz_0_xxxy, t_0_xxz_0_xxxz, t_0_xxz_0_xxy, t_0_xxz_0_xxyy,\
                           t_0_xxz_0_xxyz, t_0_xxz_0_xxz, t_0_xxz_0_xxzz, t_0_xxz_0_xyy,\
                           t_0_xxz_0_xyyy, t_0_xxz_0_xyyz, t_0_xxz_0_xyz, t_0_xxz_0_xyzz,\
                           t_0_xxz_0_xzz, t_0_xxz_0_xzzz, t_0_xxz_0_yyy, t_0_xxz_0_yyyy,\
                           t_0_xxz_0_yyyz, t_0_xxz_0_yyz, t_0_xxz_0_yyzz, t_0_xxz_0_yzz,\
                           t_0_xxz_0_yzzz, t_0_xxz_0_zzz, t_0_xxz_x_xxx, t_0_xxz_x_xxy,\
                           t_0_xxz_x_xxz, t_0_xxz_x_xyy, t_0_xxz_x_xyz, t_0_xxz_x_xzz,\
                           t_0_xxz_x_yyy, t_0_xxz_x_yyz, t_0_xxz_x_yzz, t_0_xxz_x_zzz,\
                           t_0_xxz_y_xxy, t_0_xxz_y_xxz, t_0_xxz_y_xyy, t_0_xxz_y_xyz,\
                           t_0_xxz_y_xzz, t_0_xxz_y_yyy, t_0_xxz_y_yyz, t_0_xxz_y_yzz,\
                           t_0_xxz_y_zzz, t_0_xxz_z_xxz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxz_z_xxz[i] = t_0_xxz_0_xxzz[i] - rcd_z[i] * t_0_xxz_0_xxz[i];

        t_0_xxz_y_zzz[i] = t_0_xxz_0_yzzz[i] - rcd_y[i] * t_0_xxz_0_zzz[i];

        t_0_xxz_y_yzz[i] = t_0_xxz_0_yyzz[i] - rcd_y[i] * t_0_xxz_0_yzz[i];

        t_0_xxz_y_yyz[i] = t_0_xxz_0_yyyz[i] - rcd_y[i] * t_0_xxz_0_yyz[i];

        t_0_xxz_y_yyy[i] = t_0_xxz_0_yyyy[i] - rcd_y[i] * t_0_xxz_0_yyy[i];

        t_0_xxz_y_xzz[i] = t_0_xxz_0_xyzz[i] - rcd_y[i] * t_0_xxz_0_xzz[i];

        t_0_xxz_y_xyz[i] = t_0_xxz_0_xyyz[i] - rcd_y[i] * t_0_xxz_0_xyz[i];

        t_0_xxz_y_xyy[i] = t_0_xxz_0_xyyy[i] - rcd_y[i] * t_0_xxz_0_xyy[i];

        t_0_xxz_y_xxz[i] = t_0_xxz_0_xxyz[i] - rcd_y[i] * t_0_xxz_0_xxz[i];

        t_0_xxz_y_xxy[i] = t_0_xxz_0_xxyy[i] - rcd_y[i] * t_0_xxz_0_xxy[i];

        t_0_xxz_x_zzz[i] = t_0_xxz_0_xzzz[i] - rcd_x[i] * t_0_xxz_0_zzz[i];

        t_0_xxz_x_yzz[i] = t_0_xxz_0_xyzz[i] - rcd_x[i] * t_0_xxz_0_yzz[i];

        t_0_xxz_x_yyz[i] = t_0_xxz_0_xyyz[i] - rcd_x[i] * t_0_xxz_0_yyz[i];

        t_0_xxz_x_yyy[i] = t_0_xxz_0_xyyy[i] - rcd_x[i] * t_0_xxz_0_yyy[i];

        t_0_xxz_x_xzz[i] = t_0_xxz_0_xxzz[i] - rcd_x[i] * t_0_xxz_0_xzz[i];

        t_0_xxz_x_xyz[i] = t_0_xxz_0_xxyz[i] - rcd_x[i] * t_0_xxz_0_xyz[i];

        t_0_xxz_x_xyy[i] = t_0_xxz_0_xxyy[i] - rcd_x[i] * t_0_xxz_0_xyy[i];

        t_0_xxz_x_xxz[i] = t_0_xxz_0_xxxz[i] - rcd_x[i] * t_0_xxz_0_xxz[i];

        t_0_xxz_x_xxy[i] = t_0_xxz_0_xxxy[i] - rcd_x[i] * t_0_xxz_0_xxy[i];

        t_0_xxz_x_xxx[i] = t_0_xxz_0_xxxx[i] - rcd_x[i] * t_0_xxz_0_xxx[i];

        t_0_xxy_z_zzz[i] = t_0_xxy_0_zzzz[i] - rcd_z[i] * t_0_xxy_0_zzz[i];

        t_0_xxy_z_yzz[i] = t_0_xxy_0_yzzz[i] - rcd_z[i] * t_0_xxy_0_yzz[i];

        t_0_xxy_z_yyz[i] = t_0_xxy_0_yyzz[i] - rcd_z[i] * t_0_xxy_0_yyz[i];

        t_0_xxy_z_xzz[i] = t_0_xxy_0_xzzz[i] - rcd_z[i] * t_0_xxy_0_xzz[i];

        t_0_xxy_z_xyz[i] = t_0_xxy_0_xyzz[i] - rcd_z[i] * t_0_xxy_0_xyz[i];

        t_0_xxy_z_xxz[i] = t_0_xxy_0_xxzz[i] - rcd_z[i] * t_0_xxy_0_xxz[i];

        t_0_xxy_y_zzz[i] = t_0_xxy_0_yzzz[i] - rcd_y[i] * t_0_xxy_0_zzz[i];

        t_0_xxy_y_yzz[i] = t_0_xxy_0_yyzz[i] - rcd_y[i] * t_0_xxy_0_yzz[i];

        t_0_xxy_y_yyz[i] = t_0_xxy_0_yyyz[i] - rcd_y[i] * t_0_xxy_0_yyz[i];

        t_0_xxy_y_yyy[i] = t_0_xxy_0_yyyy[i] - rcd_y[i] * t_0_xxy_0_yyy[i];

        t_0_xxy_y_xzz[i] = t_0_xxy_0_xyzz[i] - rcd_y[i] * t_0_xxy_0_xzz[i];

        t_0_xxy_y_xyz[i] = t_0_xxy_0_xyyz[i] - rcd_y[i] * t_0_xxy_0_xyz[i];

        t_0_xxy_y_xyy[i] = t_0_xxy_0_xyyy[i] - rcd_y[i] * t_0_xxy_0_xyy[i];

        t_0_xxy_y_xxz[i] = t_0_xxy_0_xxyz[i] - rcd_y[i] * t_0_xxy_0_xxz[i];

        t_0_xxy_y_xxy[i] = t_0_xxy_0_xxyy[i] - rcd_y[i] * t_0_xxy_0_xxy[i];

        t_0_xxy_x_zzz[i] = t_0_xxy_0_xzzz[i] - rcd_x[i] * t_0_xxy_0_zzz[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_0_xxx, t_0_xxx_0_xxxx, t_0_xxx_0_xxxy,\
                           t_0_xxx_0_xxxz, t_0_xxx_0_xxy, t_0_xxx_0_xxyy, t_0_xxx_0_xxyz,\
                           t_0_xxx_0_xxz, t_0_xxx_0_xxzz, t_0_xxx_0_xyy, t_0_xxx_0_xyyy,\
                           t_0_xxx_0_xyyz, t_0_xxx_0_xyz, t_0_xxx_0_xyzz, t_0_xxx_0_xzz,\
                           t_0_xxx_0_xzzz, t_0_xxx_0_yyy, t_0_xxx_0_yyyy, t_0_xxx_0_yyyz,\
                           t_0_xxx_0_yyz, t_0_xxx_0_yyzz, t_0_xxx_0_yzz, t_0_xxx_0_yzzz,\
                           t_0_xxx_0_zzz, t_0_xxx_0_zzzz, t_0_xxx_x_xxx, t_0_xxx_x_xxy,\
                           t_0_xxx_x_xxz, t_0_xxx_x_xyy, t_0_xxx_x_xyz, t_0_xxx_x_xzz,\
                           t_0_xxx_x_yyy, t_0_xxx_x_yyz, t_0_xxx_x_yzz, t_0_xxx_x_zzz,\
                           t_0_xxx_y_xxy, t_0_xxx_y_xxz, t_0_xxx_y_xyy, t_0_xxx_y_xyz,\
                           t_0_xxx_y_xzz, t_0_xxx_y_yyy, t_0_xxx_y_yyz, t_0_xxx_y_yzz,\
                           t_0_xxx_y_zzz, t_0_xxx_z_xxz, t_0_xxx_z_xyz, t_0_xxx_z_xzz,\
                           t_0_xxx_z_yyz, t_0_xxx_z_yzz, t_0_xxx_z_zzz, t_0_xxy_0_xxx,\
                           t_0_xxy_0_xxxx, t_0_xxy_0_xxxy, t_0_xxy_0_xxxz, t_0_xxy_0_xxy,\
                           t_0_xxy_0_xxyy, t_0_xxy_0_xxyz, t_0_xxy_0_xxz, t_0_xxy_0_xxzz,\
                           t_0_xxy_0_xyy, t_0_xxy_0_xyyy, t_0_xxy_0_xyyz, t_0_xxy_0_xyz,\
                           t_0_xxy_0_xyzz, t_0_xxy_0_xzz, t_0_xxy_0_yyy, t_0_xxy_0_yyz,\
                           t_0_xxy_0_yzz, t_0_xxy_x_xxx, t_0_xxy_x_xxy, t_0_xxy_x_xxz,\
                           t_0_xxy_x_xyy, t_0_xxy_x_xyz, t_0_xxy_x_xzz, t_0_xxy_x_yyy,\
                           t_0_xxy_x_yyz, t_0_xxy_x_yzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxy_x_yzz[i] = t_0_xxy_0_xyzz[i] - rcd_x[i] * t_0_xxy_0_yzz[i];

        t_0_xxy_x_yyz[i] = t_0_xxy_0_xyyz[i] - rcd_x[i] * t_0_xxy_0_yyz[i];

        t_0_xxy_x_yyy[i] = t_0_xxy_0_xyyy[i] - rcd_x[i] * t_0_xxy_0_yyy[i];

        t_0_xxy_x_xzz[i] = t_0_xxy_0_xxzz[i] - rcd_x[i] * t_0_xxy_0_xzz[i];

        t_0_xxy_x_xyz[i] = t_0_xxy_0_xxyz[i] - rcd_x[i] * t_0_xxy_0_xyz[i];

        t_0_xxy_x_xyy[i] = t_0_xxy_0_xxyy[i] - rcd_x[i] * t_0_xxy_0_xyy[i];

        t_0_xxy_x_xxz[i] = t_0_xxy_0_xxxz[i] - rcd_x[i] * t_0_xxy_0_xxz[i];

        t_0_xxy_x_xxy[i] = t_0_xxy_0_xxxy[i] - rcd_x[i] * t_0_xxy_0_xxy[i];

        t_0_xxy_x_xxx[i] = t_0_xxy_0_xxxx[i] - rcd_x[i] * t_0_xxy_0_xxx[i];

        t_0_xxx_z_zzz[i] = t_0_xxx_0_zzzz[i] - rcd_z[i] * t_0_xxx_0_zzz[i];

        t_0_xxx_z_yzz[i] = t_0_xxx_0_yzzz[i] - rcd_z[i] * t_0_xxx_0_yzz[i];

        t_0_xxx_z_yyz[i] = t_0_xxx_0_yyzz[i] - rcd_z[i] * t_0_xxx_0_yyz[i];

        t_0_xxx_z_xzz[i] = t_0_xxx_0_xzzz[i] - rcd_z[i] * t_0_xxx_0_xzz[i];

        t_0_xxx_z_xyz[i] = t_0_xxx_0_xyzz[i] - rcd_z[i] * t_0_xxx_0_xyz[i];

        t_0_xxx_z_xxz[i] = t_0_xxx_0_xxzz[i] - rcd_z[i] * t_0_xxx_0_xxz[i];

        t_0_xxx_y_zzz[i] = t_0_xxx_0_yzzz[i] - rcd_y[i] * t_0_xxx_0_zzz[i];

        t_0_xxx_y_yzz[i] = t_0_xxx_0_yyzz[i] - rcd_y[i] * t_0_xxx_0_yzz[i];

        t_0_xxx_y_yyz[i] = t_0_xxx_0_yyyz[i] - rcd_y[i] * t_0_xxx_0_yyz[i];

        t_0_xxx_y_yyy[i] = t_0_xxx_0_yyyy[i] - rcd_y[i] * t_0_xxx_0_yyy[i];

        t_0_xxx_y_xzz[i] = t_0_xxx_0_xyzz[i] - rcd_y[i] * t_0_xxx_0_xzz[i];

        t_0_xxx_y_xyz[i] = t_0_xxx_0_xyyz[i] - rcd_y[i] * t_0_xxx_0_xyz[i];

        t_0_xxx_y_xyy[i] = t_0_xxx_0_xyyy[i] - rcd_y[i] * t_0_xxx_0_xyy[i];

        t_0_xxx_y_xxz[i] = t_0_xxx_0_xxyz[i] - rcd_y[i] * t_0_xxx_0_xxz[i];

        t_0_xxx_y_xxy[i] = t_0_xxx_0_xxyy[i] - rcd_y[i] * t_0_xxx_0_xxy[i];

        t_0_xxx_x_zzz[i] = t_0_xxx_0_xzzz[i] - rcd_x[i] * t_0_xxx_0_zzz[i];

        t_0_xxx_x_yzz[i] = t_0_xxx_0_xyzz[i] - rcd_x[i] * t_0_xxx_0_yzz[i];

        t_0_xxx_x_yyz[i] = t_0_xxx_0_xyyz[i] - rcd_x[i] * t_0_xxx_0_yyz[i];

        t_0_xxx_x_yyy[i] = t_0_xxx_0_xyyy[i] - rcd_x[i] * t_0_xxx_0_yyy[i];

        t_0_xxx_x_xzz[i] = t_0_xxx_0_xxzz[i] - rcd_x[i] * t_0_xxx_0_xzz[i];

        t_0_xxx_x_xyz[i] = t_0_xxx_0_xxyz[i] - rcd_x[i] * t_0_xxx_0_xyz[i];

        t_0_xxx_x_xyy[i] = t_0_xxx_0_xxyy[i] - rcd_x[i] * t_0_xxx_0_xyy[i];

        t_0_xxx_x_xxz[i] = t_0_xxx_0_xxxz[i] - rcd_x[i] * t_0_xxx_0_xxz[i];

        t_0_xxx_x_xxy[i] = t_0_xxx_0_xxxy[i] - rcd_x[i] * t_0_xxx_0_xxy[i];

        t_0_xxx_x_xxx[i] = t_0_xxx_0_xxxx[i] - rcd_x[i] * t_0_xxx_0_xxx[i];
    }
}


} // derirec namespace
