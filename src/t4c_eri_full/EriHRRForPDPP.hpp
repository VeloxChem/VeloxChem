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
compHostHRRForPDPP_V0(      BufferHostXY<T>&      intsBufferPDPP,
                      const BufferHostX<int32_t>& intsIndexesPDPP,
                      const BufferHostXY<T>&      intsBufferSDPP,
                      const BufferHostX<int32_t>& intsIndexesSDPP,
                      const BufferHostXY<T>&      intsBufferSFPP,
                      const BufferHostX<int32_t>& intsIndexesSFPP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SDPP) integral components

    t_0_zz_z_z = intsBufferSDPP.data(intsIndexesSDPP(0));

    t_0_zz_z_y = intsBufferSDPP.data(intsIndexesSDPP(1));

    t_0_zz_z_x = intsBufferSDPP.data(intsIndexesSDPP(2));

    t_0_zz_y_z = intsBufferSDPP.data(intsIndexesSDPP(3));

    t_0_zz_y_y = intsBufferSDPP.data(intsIndexesSDPP(4));

    t_0_zz_y_x = intsBufferSDPP.data(intsIndexesSDPP(5));

    t_0_zz_x_z = intsBufferSDPP.data(intsIndexesSDPP(6));

    t_0_zz_x_y = intsBufferSDPP.data(intsIndexesSDPP(7));

    t_0_zz_x_x = intsBufferSDPP.data(intsIndexesSDPP(8));

    t_0_yz_z_z = intsBufferSDPP.data(intsIndexesSDPP(9));

    t_0_yz_z_y = intsBufferSDPP.data(intsIndexesSDPP(10));

    t_0_yz_z_x = intsBufferSDPP.data(intsIndexesSDPP(11));

    t_0_yz_y_z = intsBufferSDPP.data(intsIndexesSDPP(12));

    t_0_yz_y_y = intsBufferSDPP.data(intsIndexesSDPP(13));

    t_0_yz_y_x = intsBufferSDPP.data(intsIndexesSDPP(14));

    t_0_yz_x_z = intsBufferSDPP.data(intsIndexesSDPP(15));

    t_0_yz_x_y = intsBufferSDPP.data(intsIndexesSDPP(16));

    t_0_yz_x_x = intsBufferSDPP.data(intsIndexesSDPP(17));

    t_0_yy_z_z = intsBufferSDPP.data(intsIndexesSDPP(18));

    t_0_yy_z_y = intsBufferSDPP.data(intsIndexesSDPP(19));

    t_0_yy_z_x = intsBufferSDPP.data(intsIndexesSDPP(20));

    t_0_yy_y_z = intsBufferSDPP.data(intsIndexesSDPP(21));

    t_0_yy_y_y = intsBufferSDPP.data(intsIndexesSDPP(22));

    t_0_yy_y_x = intsBufferSDPP.data(intsIndexesSDPP(23));

    t_0_yy_x_z = intsBufferSDPP.data(intsIndexesSDPP(24));

    t_0_yy_x_y = intsBufferSDPP.data(intsIndexesSDPP(25));

    t_0_yy_x_x = intsBufferSDPP.data(intsIndexesSDPP(26));

    t_0_xz_z_z = intsBufferSDPP.data(intsIndexesSDPP(27));

    t_0_xz_z_y = intsBufferSDPP.data(intsIndexesSDPP(28));

    t_0_xz_z_x = intsBufferSDPP.data(intsIndexesSDPP(29));

    t_0_xz_y_z = intsBufferSDPP.data(intsIndexesSDPP(30));

    t_0_xz_y_y = intsBufferSDPP.data(intsIndexesSDPP(31));

    t_0_xz_y_x = intsBufferSDPP.data(intsIndexesSDPP(32));

    t_0_xz_x_z = intsBufferSDPP.data(intsIndexesSDPP(33));

    t_0_xz_x_y = intsBufferSDPP.data(intsIndexesSDPP(34));

    t_0_xz_x_x = intsBufferSDPP.data(intsIndexesSDPP(35));

    t_0_xy_z_z = intsBufferSDPP.data(intsIndexesSDPP(36));

    t_0_xy_z_y = intsBufferSDPP.data(intsIndexesSDPP(37));

    t_0_xy_z_x = intsBufferSDPP.data(intsIndexesSDPP(38));

    t_0_xy_y_z = intsBufferSDPP.data(intsIndexesSDPP(39));

    t_0_xy_y_y = intsBufferSDPP.data(intsIndexesSDPP(40));

    t_0_xy_y_x = intsBufferSDPP.data(intsIndexesSDPP(41));

    t_0_xy_x_z = intsBufferSDPP.data(intsIndexesSDPP(42));

    t_0_xy_x_y = intsBufferSDPP.data(intsIndexesSDPP(43));

    t_0_xy_x_x = intsBufferSDPP.data(intsIndexesSDPP(44));

    t_0_xx_z_z = intsBufferSDPP.data(intsIndexesSDPP(45));

    t_0_xx_z_y = intsBufferSDPP.data(intsIndexesSDPP(46));

    t_0_xx_z_x = intsBufferSDPP.data(intsIndexesSDPP(47));

    t_0_xx_y_z = intsBufferSDPP.data(intsIndexesSDPP(48));

    t_0_xx_y_y = intsBufferSDPP.data(intsIndexesSDPP(49));

    t_0_xx_y_x = intsBufferSDPP.data(intsIndexesSDPP(50));

    t_0_xx_x_z = intsBufferSDPP.data(intsIndexesSDPP(51));

    t_0_xx_x_y = intsBufferSDPP.data(intsIndexesSDPP(52));

    t_0_xx_x_x = intsBufferSDPP.data(intsIndexesSDPP(53));

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

    #pragma omp simd align(rab_z, t_0_xz_x_x, t_0_xz_x_y, t_0_xz_x_z, t_0_xz_y_x, t_0_xz_y_y,\
                           t_0_xz_y_z, t_0_xz_z_x, t_0_xz_z_y, t_0_xz_z_z, t_0_xzz_x_x,\
                           t_0_xzz_x_y, t_0_xzz_x_z, t_0_xzz_y_x, t_0_xzz_y_y, t_0_xzz_y_z,\
                           t_0_xzz_z_x, t_0_xzz_z_y, t_0_xzz_z_z, t_0_yy_x_x, t_0_yy_x_y,\
                           t_0_yy_x_z, t_0_yy_y_x, t_0_yy_y_y, t_0_yy_y_z, t_0_yy_z_x,\
                           t_0_yy_z_y, t_0_yy_z_z, t_0_yyz_x_x, t_0_yyz_x_y, t_0_yyz_x_z,\
                           t_0_yyz_y_x, t_0_yyz_y_y, t_0_yyz_y_z, t_0_yyz_z_x, t_0_yyz_z_y,\
                           t_0_yyz_z_z, t_0_yz_x_x, t_0_yz_x_y, t_0_yz_x_z, t_0_yz_y_x,\
                           t_0_yz_y_y, t_0_yz_y_z, t_0_yz_z_x, t_0_yz_z_y, t_0_yz_z_z,\
                           t_0_yzz_x_x, t_0_yzz_x_y, t_0_yzz_x_z, t_0_yzz_y_x, t_0_yzz_y_y,\
                           t_0_yzz_y_z, t_0_yzz_z_x, t_0_yzz_z_y, t_0_yzz_z_z, t_0_zz_x_x,\
                           t_0_zz_x_y, t_0_zz_x_z, t_0_zz_y_x, t_0_zz_y_y, t_0_zz_y_z,\
                           t_0_zz_z_x, t_0_zz_z_y, t_0_zz_z_z, t_0_zzz_x_x, t_0_zzz_x_y,\
                           t_0_zzz_x_z, t_0_zzz_y_x, t_0_zzz_y_y, t_0_zzz_y_z, t_0_zzz_z_x,\
                           t_0_zzz_z_y, t_0_zzz_z_z, t_z_xz_x_x, t_z_xz_x_y, t_z_xz_x_z,\
                           t_z_xz_y_x, t_z_xz_y_y, t_z_xz_y_z, t_z_xz_z_x, t_z_xz_z_y,\
                           t_z_xz_z_z, t_z_yy_x_x, t_z_yy_x_y, t_z_yy_x_z, t_z_yy_y_x,\
                           t_z_yy_y_y, t_z_yy_y_z, t_z_yy_z_x, t_z_yy_z_y, t_z_yy_z_z,\
                           t_z_yz_x_x, t_z_yz_x_y, t_z_yz_x_z, t_z_yz_y_x, t_z_yz_y_y,\
                           t_z_yz_y_z, t_z_yz_z_x, t_z_yz_z_y, t_z_yz_z_z, t_z_zz_x_x,\
                           t_z_zz_x_y, t_z_zz_x_z, t_z_zz_y_x, t_z_zz_y_y, t_z_zz_y_z,\
                           t_z_zz_z_x, t_z_zz_z_y, t_z_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_z_z[i] = t_0_zzz_z_z[i] - rab_z[i] * t_0_zz_z_z[i];

        t_z_zz_z_y[i] = t_0_zzz_z_y[i] - rab_z[i] * t_0_zz_z_y[i];

        t_z_zz_z_x[i] = t_0_zzz_z_x[i] - rab_z[i] * t_0_zz_z_x[i];

        t_z_zz_y_z[i] = t_0_zzz_y_z[i] - rab_z[i] * t_0_zz_y_z[i];

        t_z_zz_y_y[i] = t_0_zzz_y_y[i] - rab_z[i] * t_0_zz_y_y[i];

        t_z_zz_y_x[i] = t_0_zzz_y_x[i] - rab_z[i] * t_0_zz_y_x[i];

        t_z_zz_x_z[i] = t_0_zzz_x_z[i] - rab_z[i] * t_0_zz_x_z[i];

        t_z_zz_x_y[i] = t_0_zzz_x_y[i] - rab_z[i] * t_0_zz_x_y[i];

        t_z_zz_x_x[i] = t_0_zzz_x_x[i] - rab_z[i] * t_0_zz_x_x[i];

        t_z_yz_z_z[i] = t_0_yzz_z_z[i] - rab_z[i] * t_0_yz_z_z[i];

        t_z_yz_z_y[i] = t_0_yzz_z_y[i] - rab_z[i] * t_0_yz_z_y[i];

        t_z_yz_z_x[i] = t_0_yzz_z_x[i] - rab_z[i] * t_0_yz_z_x[i];

        t_z_yz_y_z[i] = t_0_yzz_y_z[i] - rab_z[i] * t_0_yz_y_z[i];

        t_z_yz_y_y[i] = t_0_yzz_y_y[i] - rab_z[i] * t_0_yz_y_y[i];

        t_z_yz_y_x[i] = t_0_yzz_y_x[i] - rab_z[i] * t_0_yz_y_x[i];

        t_z_yz_x_z[i] = t_0_yzz_x_z[i] - rab_z[i] * t_0_yz_x_z[i];

        t_z_yz_x_y[i] = t_0_yzz_x_y[i] - rab_z[i] * t_0_yz_x_y[i];

        t_z_yz_x_x[i] = t_0_yzz_x_x[i] - rab_z[i] * t_0_yz_x_x[i];

        t_z_yy_z_z[i] = t_0_yyz_z_z[i] - rab_z[i] * t_0_yy_z_z[i];

        t_z_yy_z_y[i] = t_0_yyz_z_y[i] - rab_z[i] * t_0_yy_z_y[i];

        t_z_yy_z_x[i] = t_0_yyz_z_x[i] - rab_z[i] * t_0_yy_z_x[i];

        t_z_yy_y_z[i] = t_0_yyz_y_z[i] - rab_z[i] * t_0_yy_y_z[i];

        t_z_yy_y_y[i] = t_0_yyz_y_y[i] - rab_z[i] * t_0_yy_y_y[i];

        t_z_yy_y_x[i] = t_0_yyz_y_x[i] - rab_z[i] * t_0_yy_y_x[i];

        t_z_yy_x_z[i] = t_0_yyz_x_z[i] - rab_z[i] * t_0_yy_x_z[i];

        t_z_yy_x_y[i] = t_0_yyz_x_y[i] - rab_z[i] * t_0_yy_x_y[i];

        t_z_yy_x_x[i] = t_0_yyz_x_x[i] - rab_z[i] * t_0_yy_x_x[i];

        t_z_xz_z_z[i] = t_0_xzz_z_z[i] - rab_z[i] * t_0_xz_z_z[i];

        t_z_xz_z_y[i] = t_0_xzz_z_y[i] - rab_z[i] * t_0_xz_z_y[i];

        t_z_xz_z_x[i] = t_0_xzz_z_x[i] - rab_z[i] * t_0_xz_z_x[i];

        t_z_xz_y_z[i] = t_0_xzz_y_z[i] - rab_z[i] * t_0_xz_y_z[i];

        t_z_xz_y_y[i] = t_0_xzz_y_y[i] - rab_z[i] * t_0_xz_y_y[i];

        t_z_xz_y_x[i] = t_0_xzz_y_x[i] - rab_z[i] * t_0_xz_y_x[i];

        t_z_xz_x_z[i] = t_0_xzz_x_z[i] - rab_z[i] * t_0_xz_x_z[i];

        t_z_xz_x_y[i] = t_0_xzz_x_y[i] - rab_z[i] * t_0_xz_x_y[i];

        t_z_xz_x_x[i] = t_0_xzz_x_x[i] - rab_z[i] * t_0_xz_x_x[i];
    }

    #pragma omp simd align(rab_y, rab_z, t_0_xx_x_x, t_0_xx_x_y, t_0_xx_x_z, t_0_xx_y_x,\
                           t_0_xx_y_y, t_0_xx_y_z, t_0_xx_z_x, t_0_xx_z_y, t_0_xx_z_z,\
                           t_0_xxz_x_x, t_0_xxz_x_y, t_0_xxz_x_z, t_0_xxz_y_x, t_0_xxz_y_y,\
                           t_0_xxz_y_z, t_0_xxz_z_x, t_0_xxz_z_y, t_0_xxz_z_z, t_0_xy_x_x,\
                           t_0_xy_x_y, t_0_xy_x_z, t_0_xy_y_x, t_0_xy_y_y, t_0_xy_y_z,\
                           t_0_xy_z_x, t_0_xy_z_y, t_0_xy_z_z, t_0_xyz_x_x, t_0_xyz_x_y,\
                           t_0_xyz_x_z, t_0_xyz_y_x, t_0_xyz_y_y, t_0_xyz_y_z, t_0_xyz_z_x,\
                           t_0_xyz_z_y, t_0_xyz_z_z, t_0_yyz_x_x, t_0_yyz_x_y, t_0_yyz_x_z,\
                           t_0_yyz_y_x, t_0_yyz_y_y, t_0_yyz_y_z, t_0_yyz_z_x, t_0_yyz_z_y,\
                           t_0_yyz_z_z, t_0_yz_x_x, t_0_yz_x_y, t_0_yz_x_z, t_0_yz_y_x,\
                           t_0_yz_y_y, t_0_yz_y_z, t_0_yz_z_x, t_0_yz_z_y, t_0_yz_z_z,\
                           t_0_yzz_x_x, t_0_yzz_x_y, t_0_yzz_x_z, t_0_yzz_y_x, t_0_yzz_y_y,\
                           t_0_yzz_y_z, t_0_yzz_z_x, t_0_yzz_z_y, t_0_yzz_z_z, t_0_zz_x_x,\
                           t_0_zz_x_y, t_0_zz_x_z, t_0_zz_y_x, t_0_zz_y_y, t_0_zz_y_z,\
                           t_0_zz_z_x, t_0_zz_z_y, t_0_zz_z_z, t_y_yz_x_x, t_y_yz_x_y,\
                           t_y_yz_x_z, t_y_yz_y_x, t_y_yz_y_y, t_y_yz_y_z, t_y_yz_z_x,\
                           t_y_yz_z_y, t_y_yz_z_z, t_y_zz_x_x, t_y_zz_x_y, t_y_zz_x_z,\
                           t_y_zz_y_x, t_y_zz_y_y, t_y_zz_y_z, t_y_zz_z_x, t_y_zz_z_y,\
                           t_y_zz_z_z, t_z_xx_x_x, t_z_xx_x_y, t_z_xx_x_z, t_z_xx_y_x,\
                           t_z_xx_y_y, t_z_xx_y_z, t_z_xx_z_x, t_z_xx_z_y, t_z_xx_z_z,\
                           t_z_xy_x_x, t_z_xy_x_y, t_z_xy_x_z, t_z_xy_y_x, t_z_xy_y_y,\
                           t_z_xy_y_z, t_z_xy_z_x, t_z_xy_z_y, t_z_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_xy_z_z[i] = t_0_xyz_z_z[i] - rab_z[i] * t_0_xy_z_z[i];

        t_z_xy_z_y[i] = t_0_xyz_z_y[i] - rab_z[i] * t_0_xy_z_y[i];

        t_z_xy_z_x[i] = t_0_xyz_z_x[i] - rab_z[i] * t_0_xy_z_x[i];

        t_z_xy_y_z[i] = t_0_xyz_y_z[i] - rab_z[i] * t_0_xy_y_z[i];

        t_z_xy_y_y[i] = t_0_xyz_y_y[i] - rab_z[i] * t_0_xy_y_y[i];

        t_z_xy_y_x[i] = t_0_xyz_y_x[i] - rab_z[i] * t_0_xy_y_x[i];

        t_z_xy_x_z[i] = t_0_xyz_x_z[i] - rab_z[i] * t_0_xy_x_z[i];

        t_z_xy_x_y[i] = t_0_xyz_x_y[i] - rab_z[i] * t_0_xy_x_y[i];

        t_z_xy_x_x[i] = t_0_xyz_x_x[i] - rab_z[i] * t_0_xy_x_x[i];

        t_z_xx_z_z[i] = t_0_xxz_z_z[i] - rab_z[i] * t_0_xx_z_z[i];

        t_z_xx_z_y[i] = t_0_xxz_z_y[i] - rab_z[i] * t_0_xx_z_y[i];

        t_z_xx_z_x[i] = t_0_xxz_z_x[i] - rab_z[i] * t_0_xx_z_x[i];

        t_z_xx_y_z[i] = t_0_xxz_y_z[i] - rab_z[i] * t_0_xx_y_z[i];

        t_z_xx_y_y[i] = t_0_xxz_y_y[i] - rab_z[i] * t_0_xx_y_y[i];

        t_z_xx_y_x[i] = t_0_xxz_y_x[i] - rab_z[i] * t_0_xx_y_x[i];

        t_z_xx_x_z[i] = t_0_xxz_x_z[i] - rab_z[i] * t_0_xx_x_z[i];

        t_z_xx_x_y[i] = t_0_xxz_x_y[i] - rab_z[i] * t_0_xx_x_y[i];

        t_z_xx_x_x[i] = t_0_xxz_x_x[i] - rab_z[i] * t_0_xx_x_x[i];

        t_y_zz_z_z[i] = t_0_yzz_z_z[i] - rab_y[i] * t_0_zz_z_z[i];

        t_y_zz_z_y[i] = t_0_yzz_z_y[i] - rab_y[i] * t_0_zz_z_y[i];

        t_y_zz_z_x[i] = t_0_yzz_z_x[i] - rab_y[i] * t_0_zz_z_x[i];

        t_y_zz_y_z[i] = t_0_yzz_y_z[i] - rab_y[i] * t_0_zz_y_z[i];

        t_y_zz_y_y[i] = t_0_yzz_y_y[i] - rab_y[i] * t_0_zz_y_y[i];

        t_y_zz_y_x[i] = t_0_yzz_y_x[i] - rab_y[i] * t_0_zz_y_x[i];

        t_y_zz_x_z[i] = t_0_yzz_x_z[i] - rab_y[i] * t_0_zz_x_z[i];

        t_y_zz_x_y[i] = t_0_yzz_x_y[i] - rab_y[i] * t_0_zz_x_y[i];

        t_y_zz_x_x[i] = t_0_yzz_x_x[i] - rab_y[i] * t_0_zz_x_x[i];

        t_y_yz_z_z[i] = t_0_yyz_z_z[i] - rab_y[i] * t_0_yz_z_z[i];

        t_y_yz_z_y[i] = t_0_yyz_z_y[i] - rab_y[i] * t_0_yz_z_y[i];

        t_y_yz_z_x[i] = t_0_yyz_z_x[i] - rab_y[i] * t_0_yz_z_x[i];

        t_y_yz_y_z[i] = t_0_yyz_y_z[i] - rab_y[i] * t_0_yz_y_z[i];

        t_y_yz_y_y[i] = t_0_yyz_y_y[i] - rab_y[i] * t_0_yz_y_y[i];

        t_y_yz_y_x[i] = t_0_yyz_y_x[i] - rab_y[i] * t_0_yz_y_x[i];

        t_y_yz_x_z[i] = t_0_yyz_x_z[i] - rab_y[i] * t_0_yz_x_z[i];

        t_y_yz_x_y[i] = t_0_yyz_x_y[i] - rab_y[i] * t_0_yz_x_y[i];

        t_y_yz_x_x[i] = t_0_yyz_x_x[i] - rab_y[i] * t_0_yz_x_x[i];
    }

    #pragma omp simd align(rab_y, t_0_xx_x_x, t_0_xx_x_y, t_0_xx_x_z, t_0_xx_y_x, t_0_xx_y_y,\
                           t_0_xx_y_z, t_0_xx_z_x, t_0_xx_z_y, t_0_xx_z_z, t_0_xxy_x_x,\
                           t_0_xxy_x_y, t_0_xxy_x_z, t_0_xxy_y_x, t_0_xxy_y_y, t_0_xxy_y_z,\
                           t_0_xxy_z_x, t_0_xxy_z_y, t_0_xxy_z_z, t_0_xy_x_x, t_0_xy_x_y,\
                           t_0_xy_x_z, t_0_xy_y_x, t_0_xy_y_y, t_0_xy_y_z, t_0_xy_z_x,\
                           t_0_xy_z_y, t_0_xy_z_z, t_0_xyy_x_x, t_0_xyy_x_y, t_0_xyy_x_z,\
                           t_0_xyy_y_x, t_0_xyy_y_y, t_0_xyy_y_z, t_0_xyy_z_x, t_0_xyy_z_y,\
                           t_0_xyy_z_z, t_0_xyz_x_x, t_0_xyz_x_y, t_0_xyz_x_z, t_0_xyz_y_x,\
                           t_0_xyz_y_y, t_0_xyz_y_z, t_0_xyz_z_x, t_0_xyz_z_y, t_0_xyz_z_z,\
                           t_0_xz_x_x, t_0_xz_x_y, t_0_xz_x_z, t_0_xz_y_x, t_0_xz_y_y,\
                           t_0_xz_y_z, t_0_xz_z_x, t_0_xz_z_y, t_0_xz_z_z, t_0_yy_x_x,\
                           t_0_yy_x_y, t_0_yy_x_z, t_0_yy_y_x, t_0_yy_y_y, t_0_yy_y_z,\
                           t_0_yy_z_x, t_0_yy_z_y, t_0_yy_z_z, t_0_yyy_x_x, t_0_yyy_x_y,\
                           t_0_yyy_x_z, t_0_yyy_y_x, t_0_yyy_y_y, t_0_yyy_y_z, t_0_yyy_z_x,\
                           t_0_yyy_z_y, t_0_yyy_z_z, t_y_xx_x_x, t_y_xx_x_y, t_y_xx_x_z,\
                           t_y_xx_y_x, t_y_xx_y_y, t_y_xx_y_z, t_y_xx_z_x, t_y_xx_z_y,\
                           t_y_xx_z_z, t_y_xy_x_x, t_y_xy_x_y, t_y_xy_x_z, t_y_xy_y_x,\
                           t_y_xy_y_y, t_y_xy_y_z, t_y_xy_z_x, t_y_xy_z_y, t_y_xy_z_z,\
                           t_y_xz_x_x, t_y_xz_x_y, t_y_xz_x_z, t_y_xz_y_x, t_y_xz_y_y,\
                           t_y_xz_y_z, t_y_xz_z_x, t_y_xz_z_y, t_y_xz_z_z, t_y_yy_x_x,\
                           t_y_yy_x_y, t_y_yy_x_z, t_y_yy_y_x, t_y_yy_y_y, t_y_yy_y_z,\
                           t_y_yy_z_x, t_y_yy_z_y, t_y_yy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_yy_z_z[i] = t_0_yyy_z_z[i] - rab_y[i] * t_0_yy_z_z[i];

        t_y_yy_z_y[i] = t_0_yyy_z_y[i] - rab_y[i] * t_0_yy_z_y[i];

        t_y_yy_z_x[i] = t_0_yyy_z_x[i] - rab_y[i] * t_0_yy_z_x[i];

        t_y_yy_y_z[i] = t_0_yyy_y_z[i] - rab_y[i] * t_0_yy_y_z[i];

        t_y_yy_y_y[i] = t_0_yyy_y_y[i] - rab_y[i] * t_0_yy_y_y[i];

        t_y_yy_y_x[i] = t_0_yyy_y_x[i] - rab_y[i] * t_0_yy_y_x[i];

        t_y_yy_x_z[i] = t_0_yyy_x_z[i] - rab_y[i] * t_0_yy_x_z[i];

        t_y_yy_x_y[i] = t_0_yyy_x_y[i] - rab_y[i] * t_0_yy_x_y[i];

        t_y_yy_x_x[i] = t_0_yyy_x_x[i] - rab_y[i] * t_0_yy_x_x[i];

        t_y_xz_z_z[i] = t_0_xyz_z_z[i] - rab_y[i] * t_0_xz_z_z[i];

        t_y_xz_z_y[i] = t_0_xyz_z_y[i] - rab_y[i] * t_0_xz_z_y[i];

        t_y_xz_z_x[i] = t_0_xyz_z_x[i] - rab_y[i] * t_0_xz_z_x[i];

        t_y_xz_y_z[i] = t_0_xyz_y_z[i] - rab_y[i] * t_0_xz_y_z[i];

        t_y_xz_y_y[i] = t_0_xyz_y_y[i] - rab_y[i] * t_0_xz_y_y[i];

        t_y_xz_y_x[i] = t_0_xyz_y_x[i] - rab_y[i] * t_0_xz_y_x[i];

        t_y_xz_x_z[i] = t_0_xyz_x_z[i] - rab_y[i] * t_0_xz_x_z[i];

        t_y_xz_x_y[i] = t_0_xyz_x_y[i] - rab_y[i] * t_0_xz_x_y[i];

        t_y_xz_x_x[i] = t_0_xyz_x_x[i] - rab_y[i] * t_0_xz_x_x[i];

        t_y_xy_z_z[i] = t_0_xyy_z_z[i] - rab_y[i] * t_0_xy_z_z[i];

        t_y_xy_z_y[i] = t_0_xyy_z_y[i] - rab_y[i] * t_0_xy_z_y[i];

        t_y_xy_z_x[i] = t_0_xyy_z_x[i] - rab_y[i] * t_0_xy_z_x[i];

        t_y_xy_y_z[i] = t_0_xyy_y_z[i] - rab_y[i] * t_0_xy_y_z[i];

        t_y_xy_y_y[i] = t_0_xyy_y_y[i] - rab_y[i] * t_0_xy_y_y[i];

        t_y_xy_y_x[i] = t_0_xyy_y_x[i] - rab_y[i] * t_0_xy_y_x[i];

        t_y_xy_x_z[i] = t_0_xyy_x_z[i] - rab_y[i] * t_0_xy_x_z[i];

        t_y_xy_x_y[i] = t_0_xyy_x_y[i] - rab_y[i] * t_0_xy_x_y[i];

        t_y_xy_x_x[i] = t_0_xyy_x_x[i] - rab_y[i] * t_0_xy_x_x[i];

        t_y_xx_z_z[i] = t_0_xxy_z_z[i] - rab_y[i] * t_0_xx_z_z[i];

        t_y_xx_z_y[i] = t_0_xxy_z_y[i] - rab_y[i] * t_0_xx_z_y[i];

        t_y_xx_z_x[i] = t_0_xxy_z_x[i] - rab_y[i] * t_0_xx_z_x[i];

        t_y_xx_y_z[i] = t_0_xxy_y_z[i] - rab_y[i] * t_0_xx_y_z[i];

        t_y_xx_y_y[i] = t_0_xxy_y_y[i] - rab_y[i] * t_0_xx_y_y[i];

        t_y_xx_y_x[i] = t_0_xxy_y_x[i] - rab_y[i] * t_0_xx_y_x[i];

        t_y_xx_x_z[i] = t_0_xxy_x_z[i] - rab_y[i] * t_0_xx_x_z[i];

        t_y_xx_x_y[i] = t_0_xxy_x_y[i] - rab_y[i] * t_0_xx_x_y[i];

        t_y_xx_x_x[i] = t_0_xxy_x_x[i] - rab_y[i] * t_0_xx_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xxz_x_x, t_0_xxz_x_y, t_0_xxz_x_z, t_0_xxz_y_x,\
                           t_0_xxz_y_y, t_0_xxz_y_z, t_0_xxz_z_x, t_0_xxz_z_y, t_0_xxz_z_z,\
                           t_0_xyy_x_x, t_0_xyy_x_y, t_0_xyy_x_z, t_0_xyy_y_x, t_0_xyy_y_y,\
                           t_0_xyy_y_z, t_0_xyy_z_x, t_0_xyy_z_y, t_0_xyy_z_z, t_0_xyz_x_x,\
                           t_0_xyz_x_y, t_0_xyz_x_z, t_0_xyz_y_x, t_0_xyz_y_y, t_0_xyz_y_z,\
                           t_0_xyz_z_x, t_0_xyz_z_y, t_0_xyz_z_z, t_0_xz_x_x, t_0_xz_x_y,\
                           t_0_xz_x_z, t_0_xz_y_x, t_0_xz_y_y, t_0_xz_y_z, t_0_xz_z_x,\
                           t_0_xz_z_y, t_0_xz_z_z, t_0_xzz_x_x, t_0_xzz_x_y, t_0_xzz_x_z,\
                           t_0_xzz_y_x, t_0_xzz_y_y, t_0_xzz_y_z, t_0_xzz_z_x, t_0_xzz_z_y,\
                           t_0_xzz_z_z, t_0_yy_x_x, t_0_yy_x_y, t_0_yy_x_z, t_0_yy_y_x,\
                           t_0_yy_y_y, t_0_yy_y_z, t_0_yy_z_x, t_0_yy_z_y, t_0_yy_z_z,\
                           t_0_yz_x_x, t_0_yz_x_y, t_0_yz_x_z, t_0_yz_y_x, t_0_yz_y_y,\
                           t_0_yz_y_z, t_0_yz_z_x, t_0_yz_z_y, t_0_yz_z_z, t_0_zz_x_x,\
                           t_0_zz_x_y, t_0_zz_x_z, t_0_zz_y_x, t_0_zz_y_y, t_0_zz_y_z,\
                           t_0_zz_z_x, t_0_zz_z_y, t_0_zz_z_z, t_x_xz_x_x, t_x_xz_x_y,\
                           t_x_xz_x_z, t_x_xz_y_x, t_x_xz_y_y, t_x_xz_y_z, t_x_xz_z_x,\
                           t_x_xz_z_y, t_x_xz_z_z, t_x_yy_x_x, t_x_yy_x_y, t_x_yy_x_z,\
                           t_x_yy_y_x, t_x_yy_y_y, t_x_yy_y_z, t_x_yy_z_x, t_x_yy_z_y,\
                           t_x_yy_z_z, t_x_yz_x_x, t_x_yz_x_y, t_x_yz_x_z, t_x_yz_y_x,\
                           t_x_yz_y_y, t_x_yz_y_z, t_x_yz_z_x, t_x_yz_z_y, t_x_yz_z_z,\
                           t_x_zz_x_x, t_x_zz_x_y, t_x_zz_x_z, t_x_zz_y_x, t_x_zz_y_y,\
                           t_x_zz_y_z, t_x_zz_z_x, t_x_zz_z_y, t_x_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zz_z_z[i] = t_0_xzz_z_z[i] - rab_x[i] * t_0_zz_z_z[i];

        t_x_zz_z_y[i] = t_0_xzz_z_y[i] - rab_x[i] * t_0_zz_z_y[i];

        t_x_zz_z_x[i] = t_0_xzz_z_x[i] - rab_x[i] * t_0_zz_z_x[i];

        t_x_zz_y_z[i] = t_0_xzz_y_z[i] - rab_x[i] * t_0_zz_y_z[i];

        t_x_zz_y_y[i] = t_0_xzz_y_y[i] - rab_x[i] * t_0_zz_y_y[i];

        t_x_zz_y_x[i] = t_0_xzz_y_x[i] - rab_x[i] * t_0_zz_y_x[i];

        t_x_zz_x_z[i] = t_0_xzz_x_z[i] - rab_x[i] * t_0_zz_x_z[i];

        t_x_zz_x_y[i] = t_0_xzz_x_y[i] - rab_x[i] * t_0_zz_x_y[i];

        t_x_zz_x_x[i] = t_0_xzz_x_x[i] - rab_x[i] * t_0_zz_x_x[i];

        t_x_yz_z_z[i] = t_0_xyz_z_z[i] - rab_x[i] * t_0_yz_z_z[i];

        t_x_yz_z_y[i] = t_0_xyz_z_y[i] - rab_x[i] * t_0_yz_z_y[i];

        t_x_yz_z_x[i] = t_0_xyz_z_x[i] - rab_x[i] * t_0_yz_z_x[i];

        t_x_yz_y_z[i] = t_0_xyz_y_z[i] - rab_x[i] * t_0_yz_y_z[i];

        t_x_yz_y_y[i] = t_0_xyz_y_y[i] - rab_x[i] * t_0_yz_y_y[i];

        t_x_yz_y_x[i] = t_0_xyz_y_x[i] - rab_x[i] * t_0_yz_y_x[i];

        t_x_yz_x_z[i] = t_0_xyz_x_z[i] - rab_x[i] * t_0_yz_x_z[i];

        t_x_yz_x_y[i] = t_0_xyz_x_y[i] - rab_x[i] * t_0_yz_x_y[i];

        t_x_yz_x_x[i] = t_0_xyz_x_x[i] - rab_x[i] * t_0_yz_x_x[i];

        t_x_yy_z_z[i] = t_0_xyy_z_z[i] - rab_x[i] * t_0_yy_z_z[i];

        t_x_yy_z_y[i] = t_0_xyy_z_y[i] - rab_x[i] * t_0_yy_z_y[i];

        t_x_yy_z_x[i] = t_0_xyy_z_x[i] - rab_x[i] * t_0_yy_z_x[i];

        t_x_yy_y_z[i] = t_0_xyy_y_z[i] - rab_x[i] * t_0_yy_y_z[i];

        t_x_yy_y_y[i] = t_0_xyy_y_y[i] - rab_x[i] * t_0_yy_y_y[i];

        t_x_yy_y_x[i] = t_0_xyy_y_x[i] - rab_x[i] * t_0_yy_y_x[i];

        t_x_yy_x_z[i] = t_0_xyy_x_z[i] - rab_x[i] * t_0_yy_x_z[i];

        t_x_yy_x_y[i] = t_0_xyy_x_y[i] - rab_x[i] * t_0_yy_x_y[i];

        t_x_yy_x_x[i] = t_0_xyy_x_x[i] - rab_x[i] * t_0_yy_x_x[i];

        t_x_xz_z_z[i] = t_0_xxz_z_z[i] - rab_x[i] * t_0_xz_z_z[i];

        t_x_xz_z_y[i] = t_0_xxz_z_y[i] - rab_x[i] * t_0_xz_z_y[i];

        t_x_xz_z_x[i] = t_0_xxz_z_x[i] - rab_x[i] * t_0_xz_z_x[i];

        t_x_xz_y_z[i] = t_0_xxz_y_z[i] - rab_x[i] * t_0_xz_y_z[i];

        t_x_xz_y_y[i] = t_0_xxz_y_y[i] - rab_x[i] * t_0_xz_y_y[i];

        t_x_xz_y_x[i] = t_0_xxz_y_x[i] - rab_x[i] * t_0_xz_y_x[i];

        t_x_xz_x_z[i] = t_0_xxz_x_z[i] - rab_x[i] * t_0_xz_x_z[i];

        t_x_xz_x_y[i] = t_0_xxz_x_y[i] - rab_x[i] * t_0_xz_x_y[i];

        t_x_xz_x_x[i] = t_0_xxz_x_x[i] - rab_x[i] * t_0_xz_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xx_x_x, t_0_xx_x_y, t_0_xx_x_z, t_0_xx_y_x, t_0_xx_y_y,\
                           t_0_xx_y_z, t_0_xx_z_x, t_0_xx_z_y, t_0_xx_z_z, t_0_xxx_x_x,\
                           t_0_xxx_x_y, t_0_xxx_x_z, t_0_xxx_y_x, t_0_xxx_y_y, t_0_xxx_y_z,\
                           t_0_xxx_z_x, t_0_xxx_z_y, t_0_xxx_z_z, t_0_xxy_x_x, t_0_xxy_x_y,\
                           t_0_xxy_x_z, t_0_xxy_y_x, t_0_xxy_y_y, t_0_xxy_y_z, t_0_xxy_z_x,\
                           t_0_xxy_z_y, t_0_xxy_z_z, t_0_xy_x_x, t_0_xy_x_y, t_0_xy_x_z,\
                           t_0_xy_y_x, t_0_xy_y_y, t_0_xy_y_z, t_0_xy_z_x, t_0_xy_z_y,\
                           t_0_xy_z_z, t_x_xx_x_x, t_x_xx_x_y, t_x_xx_x_z, t_x_xx_y_x,\
                           t_x_xx_y_y, t_x_xx_y_z, t_x_xx_z_x, t_x_xx_z_y, t_x_xx_z_z,\
                           t_x_xy_x_x, t_x_xy_x_y, t_x_xy_x_z, t_x_xy_y_x, t_x_xy_y_y,\
                           t_x_xy_y_z, t_x_xy_z_x, t_x_xy_z_y, t_x_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xy_z_z[i] = t_0_xxy_z_z[i] - rab_x[i] * t_0_xy_z_z[i];

        t_x_xy_z_y[i] = t_0_xxy_z_y[i] - rab_x[i] * t_0_xy_z_y[i];

        t_x_xy_z_x[i] = t_0_xxy_z_x[i] - rab_x[i] * t_0_xy_z_x[i];

        t_x_xy_y_z[i] = t_0_xxy_y_z[i] - rab_x[i] * t_0_xy_y_z[i];

        t_x_xy_y_y[i] = t_0_xxy_y_y[i] - rab_x[i] * t_0_xy_y_y[i];

        t_x_xy_y_x[i] = t_0_xxy_y_x[i] - rab_x[i] * t_0_xy_y_x[i];

        t_x_xy_x_z[i] = t_0_xxy_x_z[i] - rab_x[i] * t_0_xy_x_z[i];

        t_x_xy_x_y[i] = t_0_xxy_x_y[i] - rab_x[i] * t_0_xy_x_y[i];

        t_x_xy_x_x[i] = t_0_xxy_x_x[i] - rab_x[i] * t_0_xy_x_x[i];

        t_x_xx_z_z[i] = t_0_xxx_z_z[i] - rab_x[i] * t_0_xx_z_z[i];

        t_x_xx_z_y[i] = t_0_xxx_z_y[i] - rab_x[i] * t_0_xx_z_y[i];

        t_x_xx_z_x[i] = t_0_xxx_z_x[i] - rab_x[i] * t_0_xx_z_x[i];

        t_x_xx_y_z[i] = t_0_xxx_y_z[i] - rab_x[i] * t_0_xx_y_z[i];

        t_x_xx_y_y[i] = t_0_xxx_y_y[i] - rab_x[i] * t_0_xx_y_y[i];

        t_x_xx_y_x[i] = t_0_xxx_y_x[i] - rab_x[i] * t_0_xx_y_x[i];

        t_x_xx_x_z[i] = t_0_xxx_x_z[i] - rab_x[i] * t_0_xx_x_z[i];

        t_x_xx_x_y[i] = t_0_xxx_x_y[i] - rab_x[i] * t_0_xx_x_y[i];

        t_x_xx_x_x[i] = t_0_xxx_x_x[i] - rab_x[i] * t_0_xx_x_x[i];
    }
}


} // derirec namespace
