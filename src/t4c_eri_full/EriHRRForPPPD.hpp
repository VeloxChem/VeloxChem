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
compHostHRRForPPPD_V0(      BufferHostXY<T>&      intsBufferPPPD,
                      const BufferHostX<int32_t>& intsIndexesPPPD,
                      const BufferHostXY<T>&      intsBufferSPPD,
                      const BufferHostX<int32_t>& intsIndexesSPPD,
                      const BufferHostXY<T>&      intsBufferSDPD,
                      const BufferHostX<int32_t>& intsIndexesSDPD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPPD) integral components

    t_z_z_z_zz = intsBufferPPPD.data(intsIndexesPPPD(0));

    t_z_z_z_yz = intsBufferPPPD.data(intsIndexesPPPD(1));

    t_z_z_z_yy = intsBufferPPPD.data(intsIndexesPPPD(2));

    t_z_z_z_xz = intsBufferPPPD.data(intsIndexesPPPD(3));

    t_z_z_z_xy = intsBufferPPPD.data(intsIndexesPPPD(4));

    t_z_z_z_xx = intsBufferPPPD.data(intsIndexesPPPD(5));

    t_z_z_y_zz = intsBufferPPPD.data(intsIndexesPPPD(6));

    t_z_z_y_yz = intsBufferPPPD.data(intsIndexesPPPD(7));

    t_z_z_y_yy = intsBufferPPPD.data(intsIndexesPPPD(8));

    t_z_z_y_xz = intsBufferPPPD.data(intsIndexesPPPD(9));

    t_z_z_y_xy = intsBufferPPPD.data(intsIndexesPPPD(10));

    t_z_z_y_xx = intsBufferPPPD.data(intsIndexesPPPD(11));

    t_z_z_x_zz = intsBufferPPPD.data(intsIndexesPPPD(12));

    t_z_z_x_yz = intsBufferPPPD.data(intsIndexesPPPD(13));

    t_z_z_x_yy = intsBufferPPPD.data(intsIndexesPPPD(14));

    t_z_z_x_xz = intsBufferPPPD.data(intsIndexesPPPD(15));

    t_z_z_x_xy = intsBufferPPPD.data(intsIndexesPPPD(16));

    t_z_z_x_xx = intsBufferPPPD.data(intsIndexesPPPD(17));

    t_z_y_z_zz = intsBufferPPPD.data(intsIndexesPPPD(18));

    t_z_y_z_yz = intsBufferPPPD.data(intsIndexesPPPD(19));

    t_z_y_z_yy = intsBufferPPPD.data(intsIndexesPPPD(20));

    t_z_y_z_xz = intsBufferPPPD.data(intsIndexesPPPD(21));

    t_z_y_z_xy = intsBufferPPPD.data(intsIndexesPPPD(22));

    t_z_y_z_xx = intsBufferPPPD.data(intsIndexesPPPD(23));

    t_z_y_y_zz = intsBufferPPPD.data(intsIndexesPPPD(24));

    t_z_y_y_yz = intsBufferPPPD.data(intsIndexesPPPD(25));

    t_z_y_y_yy = intsBufferPPPD.data(intsIndexesPPPD(26));

    t_z_y_y_xz = intsBufferPPPD.data(intsIndexesPPPD(27));

    t_z_y_y_xy = intsBufferPPPD.data(intsIndexesPPPD(28));

    t_z_y_y_xx = intsBufferPPPD.data(intsIndexesPPPD(29));

    t_z_y_x_zz = intsBufferPPPD.data(intsIndexesPPPD(30));

    t_z_y_x_yz = intsBufferPPPD.data(intsIndexesPPPD(31));

    t_z_y_x_yy = intsBufferPPPD.data(intsIndexesPPPD(32));

    t_z_y_x_xz = intsBufferPPPD.data(intsIndexesPPPD(33));

    t_z_y_x_xy = intsBufferPPPD.data(intsIndexesPPPD(34));

    t_z_y_x_xx = intsBufferPPPD.data(intsIndexesPPPD(35));

    t_z_x_z_zz = intsBufferPPPD.data(intsIndexesPPPD(36));

    t_z_x_z_yz = intsBufferPPPD.data(intsIndexesPPPD(37));

    t_z_x_z_yy = intsBufferPPPD.data(intsIndexesPPPD(38));

    t_z_x_z_xz = intsBufferPPPD.data(intsIndexesPPPD(39));

    t_z_x_z_xy = intsBufferPPPD.data(intsIndexesPPPD(40));

    t_z_x_z_xx = intsBufferPPPD.data(intsIndexesPPPD(41));

    t_z_x_y_zz = intsBufferPPPD.data(intsIndexesPPPD(42));

    t_z_x_y_yz = intsBufferPPPD.data(intsIndexesPPPD(43));

    t_z_x_y_yy = intsBufferPPPD.data(intsIndexesPPPD(44));

    t_z_x_y_xz = intsBufferPPPD.data(intsIndexesPPPD(45));

    t_z_x_y_xy = intsBufferPPPD.data(intsIndexesPPPD(46));

    t_z_x_y_xx = intsBufferPPPD.data(intsIndexesPPPD(47));

    t_z_x_x_zz = intsBufferPPPD.data(intsIndexesPPPD(48));

    t_z_x_x_yz = intsBufferPPPD.data(intsIndexesPPPD(49));

    t_z_x_x_yy = intsBufferPPPD.data(intsIndexesPPPD(50));

    t_z_x_x_xz = intsBufferPPPD.data(intsIndexesPPPD(51));

    t_z_x_x_xy = intsBufferPPPD.data(intsIndexesPPPD(52));

    t_z_x_x_xx = intsBufferPPPD.data(intsIndexesPPPD(53));

    t_y_z_z_zz = intsBufferPPPD.data(intsIndexesPPPD(54));

    t_y_z_z_yz = intsBufferPPPD.data(intsIndexesPPPD(55));

    t_y_z_z_yy = intsBufferPPPD.data(intsIndexesPPPD(56));

    t_y_z_z_xz = intsBufferPPPD.data(intsIndexesPPPD(57));

    t_y_z_z_xy = intsBufferPPPD.data(intsIndexesPPPD(58));

    t_y_z_z_xx = intsBufferPPPD.data(intsIndexesPPPD(59));

    t_y_z_y_zz = intsBufferPPPD.data(intsIndexesPPPD(60));

    t_y_z_y_yz = intsBufferPPPD.data(intsIndexesPPPD(61));

    t_y_z_y_yy = intsBufferPPPD.data(intsIndexesPPPD(62));

    t_y_z_y_xz = intsBufferPPPD.data(intsIndexesPPPD(63));

    t_y_z_y_xy = intsBufferPPPD.data(intsIndexesPPPD(64));

    t_y_z_y_xx = intsBufferPPPD.data(intsIndexesPPPD(65));

    t_y_z_x_zz = intsBufferPPPD.data(intsIndexesPPPD(66));

    t_y_z_x_yz = intsBufferPPPD.data(intsIndexesPPPD(67));

    t_y_z_x_yy = intsBufferPPPD.data(intsIndexesPPPD(68));

    t_y_z_x_xz = intsBufferPPPD.data(intsIndexesPPPD(69));

    t_y_z_x_xy = intsBufferPPPD.data(intsIndexesPPPD(70));

    t_y_z_x_xx = intsBufferPPPD.data(intsIndexesPPPD(71));

    t_y_y_z_zz = intsBufferPPPD.data(intsIndexesPPPD(72));

    t_y_y_z_yz = intsBufferPPPD.data(intsIndexesPPPD(73));

    t_y_y_z_yy = intsBufferPPPD.data(intsIndexesPPPD(74));

    t_y_y_z_xz = intsBufferPPPD.data(intsIndexesPPPD(75));

    t_y_y_z_xy = intsBufferPPPD.data(intsIndexesPPPD(76));

    t_y_y_z_xx = intsBufferPPPD.data(intsIndexesPPPD(77));

    t_y_y_y_zz = intsBufferPPPD.data(intsIndexesPPPD(78));

    t_y_y_y_yz = intsBufferPPPD.data(intsIndexesPPPD(79));

    t_y_y_y_yy = intsBufferPPPD.data(intsIndexesPPPD(80));

    t_y_y_y_xz = intsBufferPPPD.data(intsIndexesPPPD(81));

    t_y_y_y_xy = intsBufferPPPD.data(intsIndexesPPPD(82));

    t_y_y_y_xx = intsBufferPPPD.data(intsIndexesPPPD(83));

    t_y_y_x_zz = intsBufferPPPD.data(intsIndexesPPPD(84));

    t_y_y_x_yz = intsBufferPPPD.data(intsIndexesPPPD(85));

    t_y_y_x_yy = intsBufferPPPD.data(intsIndexesPPPD(86));

    t_y_y_x_xz = intsBufferPPPD.data(intsIndexesPPPD(87));

    t_y_y_x_xy = intsBufferPPPD.data(intsIndexesPPPD(88));

    t_y_y_x_xx = intsBufferPPPD.data(intsIndexesPPPD(89));

    t_y_x_z_zz = intsBufferPPPD.data(intsIndexesPPPD(90));

    t_y_x_z_yz = intsBufferPPPD.data(intsIndexesPPPD(91));

    t_y_x_z_yy = intsBufferPPPD.data(intsIndexesPPPD(92));

    t_y_x_z_xz = intsBufferPPPD.data(intsIndexesPPPD(93));

    t_y_x_z_xy = intsBufferPPPD.data(intsIndexesPPPD(94));

    t_y_x_z_xx = intsBufferPPPD.data(intsIndexesPPPD(95));

    t_y_x_y_zz = intsBufferPPPD.data(intsIndexesPPPD(96));

    t_y_x_y_yz = intsBufferPPPD.data(intsIndexesPPPD(97));

    t_y_x_y_yy = intsBufferPPPD.data(intsIndexesPPPD(98));

    t_y_x_y_xz = intsBufferPPPD.data(intsIndexesPPPD(99));

    t_y_x_y_xy = intsBufferPPPD.data(intsIndexesPPPD(100));

    t_y_x_y_xx = intsBufferPPPD.data(intsIndexesPPPD(101));

    t_y_x_x_zz = intsBufferPPPD.data(intsIndexesPPPD(102));

    t_y_x_x_yz = intsBufferPPPD.data(intsIndexesPPPD(103));

    t_y_x_x_yy = intsBufferPPPD.data(intsIndexesPPPD(104));

    t_y_x_x_xz = intsBufferPPPD.data(intsIndexesPPPD(105));

    t_y_x_x_xy = intsBufferPPPD.data(intsIndexesPPPD(106));

    t_y_x_x_xx = intsBufferPPPD.data(intsIndexesPPPD(107));

    t_x_z_z_zz = intsBufferPPPD.data(intsIndexesPPPD(108));

    t_x_z_z_yz = intsBufferPPPD.data(intsIndexesPPPD(109));

    t_x_z_z_yy = intsBufferPPPD.data(intsIndexesPPPD(110));

    t_x_z_z_xz = intsBufferPPPD.data(intsIndexesPPPD(111));

    t_x_z_z_xy = intsBufferPPPD.data(intsIndexesPPPD(112));

    t_x_z_z_xx = intsBufferPPPD.data(intsIndexesPPPD(113));

    t_x_z_y_zz = intsBufferPPPD.data(intsIndexesPPPD(114));

    t_x_z_y_yz = intsBufferPPPD.data(intsIndexesPPPD(115));

    t_x_z_y_yy = intsBufferPPPD.data(intsIndexesPPPD(116));

    t_x_z_y_xz = intsBufferPPPD.data(intsIndexesPPPD(117));

    t_x_z_y_xy = intsBufferPPPD.data(intsIndexesPPPD(118));

    t_x_z_y_xx = intsBufferPPPD.data(intsIndexesPPPD(119));

    t_x_z_x_zz = intsBufferPPPD.data(intsIndexesPPPD(120));

    t_x_z_x_yz = intsBufferPPPD.data(intsIndexesPPPD(121));

    t_x_z_x_yy = intsBufferPPPD.data(intsIndexesPPPD(122));

    t_x_z_x_xz = intsBufferPPPD.data(intsIndexesPPPD(123));

    t_x_z_x_xy = intsBufferPPPD.data(intsIndexesPPPD(124));

    t_x_z_x_xx = intsBufferPPPD.data(intsIndexesPPPD(125));

    t_x_y_z_zz = intsBufferPPPD.data(intsIndexesPPPD(126));

    t_x_y_z_yz = intsBufferPPPD.data(intsIndexesPPPD(127));

    t_x_y_z_yy = intsBufferPPPD.data(intsIndexesPPPD(128));

    t_x_y_z_xz = intsBufferPPPD.data(intsIndexesPPPD(129));

    t_x_y_z_xy = intsBufferPPPD.data(intsIndexesPPPD(130));

    t_x_y_z_xx = intsBufferPPPD.data(intsIndexesPPPD(131));

    t_x_y_y_zz = intsBufferPPPD.data(intsIndexesPPPD(132));

    t_x_y_y_yz = intsBufferPPPD.data(intsIndexesPPPD(133));

    t_x_y_y_yy = intsBufferPPPD.data(intsIndexesPPPD(134));

    t_x_y_y_xz = intsBufferPPPD.data(intsIndexesPPPD(135));

    t_x_y_y_xy = intsBufferPPPD.data(intsIndexesPPPD(136));

    t_x_y_y_xx = intsBufferPPPD.data(intsIndexesPPPD(137));

    t_x_y_x_zz = intsBufferPPPD.data(intsIndexesPPPD(138));

    t_x_y_x_yz = intsBufferPPPD.data(intsIndexesPPPD(139));

    t_x_y_x_yy = intsBufferPPPD.data(intsIndexesPPPD(140));

    t_x_y_x_xz = intsBufferPPPD.data(intsIndexesPPPD(141));

    t_x_y_x_xy = intsBufferPPPD.data(intsIndexesPPPD(142));

    t_x_y_x_xx = intsBufferPPPD.data(intsIndexesPPPD(143));

    t_x_x_z_zz = intsBufferPPPD.data(intsIndexesPPPD(144));

    t_x_x_z_yz = intsBufferPPPD.data(intsIndexesPPPD(145));

    t_x_x_z_yy = intsBufferPPPD.data(intsIndexesPPPD(146));

    t_x_x_z_xz = intsBufferPPPD.data(intsIndexesPPPD(147));

    t_x_x_z_xy = intsBufferPPPD.data(intsIndexesPPPD(148));

    t_x_x_z_xx = intsBufferPPPD.data(intsIndexesPPPD(149));

    t_x_x_y_zz = intsBufferPPPD.data(intsIndexesPPPD(150));

    t_x_x_y_yz = intsBufferPPPD.data(intsIndexesPPPD(151));

    t_x_x_y_yy = intsBufferPPPD.data(intsIndexesPPPD(152));

    t_x_x_y_xz = intsBufferPPPD.data(intsIndexesPPPD(153));

    t_x_x_y_xy = intsBufferPPPD.data(intsIndexesPPPD(154));

    t_x_x_y_xx = intsBufferPPPD.data(intsIndexesPPPD(155));

    t_x_x_x_zz = intsBufferPPPD.data(intsIndexesPPPD(156));

    t_x_x_x_yz = intsBufferPPPD.data(intsIndexesPPPD(157));

    t_x_x_x_yy = intsBufferPPPD.data(intsIndexesPPPD(158));

    t_x_x_x_xz = intsBufferPPPD.data(intsIndexesPPPD(159));

    t_x_x_x_xy = intsBufferPPPD.data(intsIndexesPPPD(160));

    t_x_x_x_xx = intsBufferPPPD.data(intsIndexesPPPD(161));

    // set up (SPPD) integral components

    t_0_z_z_zz = intsBufferSPPD.data(intsIndexesSPPD(0));

    t_0_z_z_yz = intsBufferSPPD.data(intsIndexesSPPD(1));

    t_0_z_z_yy = intsBufferSPPD.data(intsIndexesSPPD(2));

    t_0_z_z_xz = intsBufferSPPD.data(intsIndexesSPPD(3));

    t_0_z_z_xy = intsBufferSPPD.data(intsIndexesSPPD(4));

    t_0_z_z_xx = intsBufferSPPD.data(intsIndexesSPPD(5));

    t_0_z_y_zz = intsBufferSPPD.data(intsIndexesSPPD(6));

    t_0_z_y_yz = intsBufferSPPD.data(intsIndexesSPPD(7));

    t_0_z_y_yy = intsBufferSPPD.data(intsIndexesSPPD(8));

    t_0_z_y_xz = intsBufferSPPD.data(intsIndexesSPPD(9));

    t_0_z_y_xy = intsBufferSPPD.data(intsIndexesSPPD(10));

    t_0_z_y_xx = intsBufferSPPD.data(intsIndexesSPPD(11));

    t_0_z_x_zz = intsBufferSPPD.data(intsIndexesSPPD(12));

    t_0_z_x_yz = intsBufferSPPD.data(intsIndexesSPPD(13));

    t_0_z_x_yy = intsBufferSPPD.data(intsIndexesSPPD(14));

    t_0_z_x_xz = intsBufferSPPD.data(intsIndexesSPPD(15));

    t_0_z_x_xy = intsBufferSPPD.data(intsIndexesSPPD(16));

    t_0_z_x_xx = intsBufferSPPD.data(intsIndexesSPPD(17));

    t_0_y_z_zz = intsBufferSPPD.data(intsIndexesSPPD(18));

    t_0_y_z_yz = intsBufferSPPD.data(intsIndexesSPPD(19));

    t_0_y_z_yy = intsBufferSPPD.data(intsIndexesSPPD(20));

    t_0_y_z_xz = intsBufferSPPD.data(intsIndexesSPPD(21));

    t_0_y_z_xy = intsBufferSPPD.data(intsIndexesSPPD(22));

    t_0_y_z_xx = intsBufferSPPD.data(intsIndexesSPPD(23));

    t_0_y_y_zz = intsBufferSPPD.data(intsIndexesSPPD(24));

    t_0_y_y_yz = intsBufferSPPD.data(intsIndexesSPPD(25));

    t_0_y_y_yy = intsBufferSPPD.data(intsIndexesSPPD(26));

    t_0_y_y_xz = intsBufferSPPD.data(intsIndexesSPPD(27));

    t_0_y_y_xy = intsBufferSPPD.data(intsIndexesSPPD(28));

    t_0_y_y_xx = intsBufferSPPD.data(intsIndexesSPPD(29));

    t_0_y_x_zz = intsBufferSPPD.data(intsIndexesSPPD(30));

    t_0_y_x_yz = intsBufferSPPD.data(intsIndexesSPPD(31));

    t_0_y_x_yy = intsBufferSPPD.data(intsIndexesSPPD(32));

    t_0_y_x_xz = intsBufferSPPD.data(intsIndexesSPPD(33));

    t_0_y_x_xy = intsBufferSPPD.data(intsIndexesSPPD(34));

    t_0_y_x_xx = intsBufferSPPD.data(intsIndexesSPPD(35));

    t_0_x_z_zz = intsBufferSPPD.data(intsIndexesSPPD(36));

    t_0_x_z_yz = intsBufferSPPD.data(intsIndexesSPPD(37));

    t_0_x_z_yy = intsBufferSPPD.data(intsIndexesSPPD(38));

    t_0_x_z_xz = intsBufferSPPD.data(intsIndexesSPPD(39));

    t_0_x_z_xy = intsBufferSPPD.data(intsIndexesSPPD(40));

    t_0_x_z_xx = intsBufferSPPD.data(intsIndexesSPPD(41));

    t_0_x_y_zz = intsBufferSPPD.data(intsIndexesSPPD(42));

    t_0_x_y_yz = intsBufferSPPD.data(intsIndexesSPPD(43));

    t_0_x_y_yy = intsBufferSPPD.data(intsIndexesSPPD(44));

    t_0_x_y_xz = intsBufferSPPD.data(intsIndexesSPPD(45));

    t_0_x_y_xy = intsBufferSPPD.data(intsIndexesSPPD(46));

    t_0_x_y_xx = intsBufferSPPD.data(intsIndexesSPPD(47));

    t_0_x_x_zz = intsBufferSPPD.data(intsIndexesSPPD(48));

    t_0_x_x_yz = intsBufferSPPD.data(intsIndexesSPPD(49));

    t_0_x_x_yy = intsBufferSPPD.data(intsIndexesSPPD(50));

    t_0_x_x_xz = intsBufferSPPD.data(intsIndexesSPPD(51));

    t_0_x_x_xy = intsBufferSPPD.data(intsIndexesSPPD(52));

    t_0_x_x_xx = intsBufferSPPD.data(intsIndexesSPPD(53));

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

    #pragma omp simd align(rab_z, t_0_y_x_xx, t_0_y_x_xy, t_0_y_x_xz, t_0_y_x_yy, t_0_y_x_yz,\
                           t_0_y_x_zz, t_0_y_y_xx, t_0_y_y_xy, t_0_y_y_xz, t_0_y_y_yy,\
                           t_0_y_y_yz, t_0_y_y_zz, t_0_y_z_xx, t_0_y_z_xy, t_0_y_z_xz,\
                           t_0_y_z_yy, t_0_y_z_yz, t_0_y_z_zz, t_0_yz_x_xx, t_0_yz_x_xy,\
                           t_0_yz_x_xz, t_0_yz_x_yy, t_0_yz_x_yz, t_0_yz_x_zz, t_0_yz_y_xx,\
                           t_0_yz_y_xy, t_0_yz_y_xz, t_0_yz_y_yy, t_0_yz_y_yz, t_0_yz_y_zz,\
                           t_0_yz_z_xx, t_0_yz_z_xy, t_0_yz_z_xz, t_0_yz_z_yy, t_0_yz_z_yz,\
                           t_0_yz_z_zz, t_0_z_x_xx, t_0_z_x_xy, t_0_z_x_xz, t_0_z_x_yy,\
                           t_0_z_x_yz, t_0_z_x_zz, t_0_z_y_xx, t_0_z_y_xy, t_0_z_y_xz,\
                           t_0_z_y_yy, t_0_z_y_yz, t_0_z_y_zz, t_0_z_z_xx, t_0_z_z_xy,\
                           t_0_z_z_xz, t_0_z_z_yy, t_0_z_z_yz, t_0_z_z_zz, t_0_zz_x_xx,\
                           t_0_zz_x_xy, t_0_zz_x_xz, t_0_zz_x_yy, t_0_zz_x_yz, t_0_zz_x_zz,\
                           t_0_zz_y_xx, t_0_zz_y_xy, t_0_zz_y_xz, t_0_zz_y_yy, t_0_zz_y_yz,\
                           t_0_zz_y_zz, t_0_zz_z_xx, t_0_zz_z_xy, t_0_zz_z_xz, t_0_zz_z_yy,\
                           t_0_zz_z_yz, t_0_zz_z_zz, t_z_y_x_xx, t_z_y_x_xy, t_z_y_x_xz,\
                           t_z_y_x_yy, t_z_y_x_yz, t_z_y_x_zz, t_z_y_y_xx, t_z_y_y_xy,\
                           t_z_y_y_xz, t_z_y_y_yy, t_z_y_y_yz, t_z_y_y_zz, t_z_y_z_xx,\
                           t_z_y_z_xy, t_z_y_z_xz, t_z_y_z_yy, t_z_y_z_yz, t_z_y_z_zz,\
                           t_z_z_x_xx, t_z_z_x_xy, t_z_z_x_xz, t_z_z_x_yy, t_z_z_x_yz,\
                           t_z_z_x_zz, t_z_z_y_xx, t_z_z_y_xy, t_z_z_y_xz, t_z_z_y_yy,\
                           t_z_z_y_yz, t_z_z_y_zz, t_z_z_z_xx, t_z_z_z_xy, t_z_z_z_xz,\
                           t_z_z_z_yy, t_z_z_z_yz, t_z_z_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_z_zz[i] = t_0_zz_z_zz[i] - rab_z[i] * t_0_z_z_zz[i];

        t_z_z_z_yz[i] = t_0_zz_z_yz[i] - rab_z[i] * t_0_z_z_yz[i];

        t_z_z_z_yy[i] = t_0_zz_z_yy[i] - rab_z[i] * t_0_z_z_yy[i];

        t_z_z_z_xz[i] = t_0_zz_z_xz[i] - rab_z[i] * t_0_z_z_xz[i];

        t_z_z_z_xy[i] = t_0_zz_z_xy[i] - rab_z[i] * t_0_z_z_xy[i];

        t_z_z_z_xx[i] = t_0_zz_z_xx[i] - rab_z[i] * t_0_z_z_xx[i];

        t_z_z_y_zz[i] = t_0_zz_y_zz[i] - rab_z[i] * t_0_z_y_zz[i];

        t_z_z_y_yz[i] = t_0_zz_y_yz[i] - rab_z[i] * t_0_z_y_yz[i];

        t_z_z_y_yy[i] = t_0_zz_y_yy[i] - rab_z[i] * t_0_z_y_yy[i];

        t_z_z_y_xz[i] = t_0_zz_y_xz[i] - rab_z[i] * t_0_z_y_xz[i];

        t_z_z_y_xy[i] = t_0_zz_y_xy[i] - rab_z[i] * t_0_z_y_xy[i];

        t_z_z_y_xx[i] = t_0_zz_y_xx[i] - rab_z[i] * t_0_z_y_xx[i];

        t_z_z_x_zz[i] = t_0_zz_x_zz[i] - rab_z[i] * t_0_z_x_zz[i];

        t_z_z_x_yz[i] = t_0_zz_x_yz[i] - rab_z[i] * t_0_z_x_yz[i];

        t_z_z_x_yy[i] = t_0_zz_x_yy[i] - rab_z[i] * t_0_z_x_yy[i];

        t_z_z_x_xz[i] = t_0_zz_x_xz[i] - rab_z[i] * t_0_z_x_xz[i];

        t_z_z_x_xy[i] = t_0_zz_x_xy[i] - rab_z[i] * t_0_z_x_xy[i];

        t_z_z_x_xx[i] = t_0_zz_x_xx[i] - rab_z[i] * t_0_z_x_xx[i];

        t_z_y_z_zz[i] = t_0_yz_z_zz[i] - rab_z[i] * t_0_y_z_zz[i];

        t_z_y_z_yz[i] = t_0_yz_z_yz[i] - rab_z[i] * t_0_y_z_yz[i];

        t_z_y_z_yy[i] = t_0_yz_z_yy[i] - rab_z[i] * t_0_y_z_yy[i];

        t_z_y_z_xz[i] = t_0_yz_z_xz[i] - rab_z[i] * t_0_y_z_xz[i];

        t_z_y_z_xy[i] = t_0_yz_z_xy[i] - rab_z[i] * t_0_y_z_xy[i];

        t_z_y_z_xx[i] = t_0_yz_z_xx[i] - rab_z[i] * t_0_y_z_xx[i];

        t_z_y_y_zz[i] = t_0_yz_y_zz[i] - rab_z[i] * t_0_y_y_zz[i];

        t_z_y_y_yz[i] = t_0_yz_y_yz[i] - rab_z[i] * t_0_y_y_yz[i];

        t_z_y_y_yy[i] = t_0_yz_y_yy[i] - rab_z[i] * t_0_y_y_yy[i];

        t_z_y_y_xz[i] = t_0_yz_y_xz[i] - rab_z[i] * t_0_y_y_xz[i];

        t_z_y_y_xy[i] = t_0_yz_y_xy[i] - rab_z[i] * t_0_y_y_xy[i];

        t_z_y_y_xx[i] = t_0_yz_y_xx[i] - rab_z[i] * t_0_y_y_xx[i];

        t_z_y_x_zz[i] = t_0_yz_x_zz[i] - rab_z[i] * t_0_y_x_zz[i];

        t_z_y_x_yz[i] = t_0_yz_x_yz[i] - rab_z[i] * t_0_y_x_yz[i];

        t_z_y_x_yy[i] = t_0_yz_x_yy[i] - rab_z[i] * t_0_y_x_yy[i];

        t_z_y_x_xz[i] = t_0_yz_x_xz[i] - rab_z[i] * t_0_y_x_xz[i];

        t_z_y_x_xy[i] = t_0_yz_x_xy[i] - rab_z[i] * t_0_y_x_xy[i];

        t_z_y_x_xx[i] = t_0_yz_x_xx[i] - rab_z[i] * t_0_y_x_xx[i];
    }

    #pragma omp simd align(rab_y, rab_z, t_0_x_x_xx, t_0_x_x_xy, t_0_x_x_xz, t_0_x_x_yy,\
                           t_0_x_x_yz, t_0_x_x_zz, t_0_x_y_xx, t_0_x_y_xy, t_0_x_y_xz,\
                           t_0_x_y_yy, t_0_x_y_yz, t_0_x_y_zz, t_0_x_z_xx, t_0_x_z_xy,\
                           t_0_x_z_xz, t_0_x_z_yy, t_0_x_z_yz, t_0_x_z_zz, t_0_xz_x_xx,\
                           t_0_xz_x_xy, t_0_xz_x_xz, t_0_xz_x_yy, t_0_xz_x_yz, t_0_xz_x_zz,\
                           t_0_xz_y_xx, t_0_xz_y_xy, t_0_xz_y_xz, t_0_xz_y_yy, t_0_xz_y_yz,\
                           t_0_xz_y_zz, t_0_xz_z_xx, t_0_xz_z_xy, t_0_xz_z_xz, t_0_xz_z_yy,\
                           t_0_xz_z_yz, t_0_xz_z_zz, t_0_yz_x_xx, t_0_yz_x_xy, t_0_yz_x_xz,\
                           t_0_yz_x_yy, t_0_yz_x_yz, t_0_yz_x_zz, t_0_yz_y_xx, t_0_yz_y_xy,\
                           t_0_yz_y_xz, t_0_yz_y_yy, t_0_yz_y_yz, t_0_yz_y_zz, t_0_yz_z_xx,\
                           t_0_yz_z_xy, t_0_yz_z_xz, t_0_yz_z_yy, t_0_yz_z_yz, t_0_yz_z_zz,\
                           t_0_z_x_xx, t_0_z_x_xy, t_0_z_x_xz, t_0_z_x_yy, t_0_z_x_yz,\
                           t_0_z_x_zz, t_0_z_y_xx, t_0_z_y_xy, t_0_z_y_xz, t_0_z_y_yy,\
                           t_0_z_y_yz, t_0_z_y_zz, t_0_z_z_xx, t_0_z_z_xy, t_0_z_z_xz,\
                           t_0_z_z_yy, t_0_z_z_yz, t_0_z_z_zz, t_y_z_x_xx, t_y_z_x_xy,\
                           t_y_z_x_xz, t_y_z_x_yy, t_y_z_x_yz, t_y_z_x_zz, t_y_z_y_xx,\
                           t_y_z_y_xy, t_y_z_y_xz, t_y_z_y_yy, t_y_z_y_yz, t_y_z_y_zz,\
                           t_y_z_z_xx, t_y_z_z_xy, t_y_z_z_xz, t_y_z_z_yy, t_y_z_z_yz,\
                           t_y_z_z_zz, t_z_x_x_xx, t_z_x_x_xy, t_z_x_x_xz, t_z_x_x_yy,\
                           t_z_x_x_yz, t_z_x_x_zz, t_z_x_y_xx, t_z_x_y_xy, t_z_x_y_xz,\
                           t_z_x_y_yy, t_z_x_y_yz, t_z_x_y_zz, t_z_x_z_xx, t_z_x_z_xy,\
                           t_z_x_z_xz, t_z_x_z_yy, t_z_x_z_yz, t_z_x_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_x_z_zz[i] = t_0_xz_z_zz[i] - rab_z[i] * t_0_x_z_zz[i];

        t_z_x_z_yz[i] = t_0_xz_z_yz[i] - rab_z[i] * t_0_x_z_yz[i];

        t_z_x_z_yy[i] = t_0_xz_z_yy[i] - rab_z[i] * t_0_x_z_yy[i];

        t_z_x_z_xz[i] = t_0_xz_z_xz[i] - rab_z[i] * t_0_x_z_xz[i];

        t_z_x_z_xy[i] = t_0_xz_z_xy[i] - rab_z[i] * t_0_x_z_xy[i];

        t_z_x_z_xx[i] = t_0_xz_z_xx[i] - rab_z[i] * t_0_x_z_xx[i];

        t_z_x_y_zz[i] = t_0_xz_y_zz[i] - rab_z[i] * t_0_x_y_zz[i];

        t_z_x_y_yz[i] = t_0_xz_y_yz[i] - rab_z[i] * t_0_x_y_yz[i];

        t_z_x_y_yy[i] = t_0_xz_y_yy[i] - rab_z[i] * t_0_x_y_yy[i];

        t_z_x_y_xz[i] = t_0_xz_y_xz[i] - rab_z[i] * t_0_x_y_xz[i];

        t_z_x_y_xy[i] = t_0_xz_y_xy[i] - rab_z[i] * t_0_x_y_xy[i];

        t_z_x_y_xx[i] = t_0_xz_y_xx[i] - rab_z[i] * t_0_x_y_xx[i];

        t_z_x_x_zz[i] = t_0_xz_x_zz[i] - rab_z[i] * t_0_x_x_zz[i];

        t_z_x_x_yz[i] = t_0_xz_x_yz[i] - rab_z[i] * t_0_x_x_yz[i];

        t_z_x_x_yy[i] = t_0_xz_x_yy[i] - rab_z[i] * t_0_x_x_yy[i];

        t_z_x_x_xz[i] = t_0_xz_x_xz[i] - rab_z[i] * t_0_x_x_xz[i];

        t_z_x_x_xy[i] = t_0_xz_x_xy[i] - rab_z[i] * t_0_x_x_xy[i];

        t_z_x_x_xx[i] = t_0_xz_x_xx[i] - rab_z[i] * t_0_x_x_xx[i];

        t_y_z_z_zz[i] = t_0_yz_z_zz[i] - rab_y[i] * t_0_z_z_zz[i];

        t_y_z_z_yz[i] = t_0_yz_z_yz[i] - rab_y[i] * t_0_z_z_yz[i];

        t_y_z_z_yy[i] = t_0_yz_z_yy[i] - rab_y[i] * t_0_z_z_yy[i];

        t_y_z_z_xz[i] = t_0_yz_z_xz[i] - rab_y[i] * t_0_z_z_xz[i];

        t_y_z_z_xy[i] = t_0_yz_z_xy[i] - rab_y[i] * t_0_z_z_xy[i];

        t_y_z_z_xx[i] = t_0_yz_z_xx[i] - rab_y[i] * t_0_z_z_xx[i];

        t_y_z_y_zz[i] = t_0_yz_y_zz[i] - rab_y[i] * t_0_z_y_zz[i];

        t_y_z_y_yz[i] = t_0_yz_y_yz[i] - rab_y[i] * t_0_z_y_yz[i];

        t_y_z_y_yy[i] = t_0_yz_y_yy[i] - rab_y[i] * t_0_z_y_yy[i];

        t_y_z_y_xz[i] = t_0_yz_y_xz[i] - rab_y[i] * t_0_z_y_xz[i];

        t_y_z_y_xy[i] = t_0_yz_y_xy[i] - rab_y[i] * t_0_z_y_xy[i];

        t_y_z_y_xx[i] = t_0_yz_y_xx[i] - rab_y[i] * t_0_z_y_xx[i];

        t_y_z_x_zz[i] = t_0_yz_x_zz[i] - rab_y[i] * t_0_z_x_zz[i];

        t_y_z_x_yz[i] = t_0_yz_x_yz[i] - rab_y[i] * t_0_z_x_yz[i];

        t_y_z_x_yy[i] = t_0_yz_x_yy[i] - rab_y[i] * t_0_z_x_yy[i];

        t_y_z_x_xz[i] = t_0_yz_x_xz[i] - rab_y[i] * t_0_z_x_xz[i];

        t_y_z_x_xy[i] = t_0_yz_x_xy[i] - rab_y[i] * t_0_z_x_xy[i];

        t_y_z_x_xx[i] = t_0_yz_x_xx[i] - rab_y[i] * t_0_z_x_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_x_x_xx, t_0_x_x_xy, t_0_x_x_xz, t_0_x_x_yy, t_0_x_x_yz,\
                           t_0_x_x_zz, t_0_x_y_xx, t_0_x_y_xy, t_0_x_y_xz, t_0_x_y_yy,\
                           t_0_x_y_yz, t_0_x_y_zz, t_0_x_z_xx, t_0_x_z_xy, t_0_x_z_xz,\
                           t_0_x_z_yy, t_0_x_z_yz, t_0_x_z_zz, t_0_xy_x_xx, t_0_xy_x_xy,\
                           t_0_xy_x_xz, t_0_xy_x_yy, t_0_xy_x_yz, t_0_xy_x_zz, t_0_xy_y_xx,\
                           t_0_xy_y_xy, t_0_xy_y_xz, t_0_xy_y_yy, t_0_xy_y_yz, t_0_xy_y_zz,\
                           t_0_xy_z_xx, t_0_xy_z_xy, t_0_xy_z_xz, t_0_xy_z_yy, t_0_xy_z_yz,\
                           t_0_xy_z_zz, t_0_y_x_xx, t_0_y_x_xy, t_0_y_x_xz, t_0_y_x_yy,\
                           t_0_y_x_yz, t_0_y_x_zz, t_0_y_y_xx, t_0_y_y_xy, t_0_y_y_xz,\
                           t_0_y_y_yy, t_0_y_y_yz, t_0_y_y_zz, t_0_y_z_xx, t_0_y_z_xy,\
                           t_0_y_z_xz, t_0_y_z_yy, t_0_y_z_yz, t_0_y_z_zz, t_0_yy_x_xx,\
                           t_0_yy_x_xy, t_0_yy_x_xz, t_0_yy_x_yy, t_0_yy_x_yz, t_0_yy_x_zz,\
                           t_0_yy_y_xx, t_0_yy_y_xy, t_0_yy_y_xz, t_0_yy_y_yy, t_0_yy_y_yz,\
                           t_0_yy_y_zz, t_0_yy_z_xx, t_0_yy_z_xy, t_0_yy_z_xz, t_0_yy_z_yy,\
                           t_0_yy_z_yz, t_0_yy_z_zz, t_y_x_x_xx, t_y_x_x_xy, t_y_x_x_xz,\
                           t_y_x_x_yy, t_y_x_x_yz, t_y_x_x_zz, t_y_x_y_xx, t_y_x_y_xy,\
                           t_y_x_y_xz, t_y_x_y_yy, t_y_x_y_yz, t_y_x_y_zz, t_y_x_z_xx,\
                           t_y_x_z_xy, t_y_x_z_xz, t_y_x_z_yy, t_y_x_z_yz, t_y_x_z_zz,\
                           t_y_y_x_xx, t_y_y_x_xy, t_y_y_x_xz, t_y_y_x_yy, t_y_y_x_yz,\
                           t_y_y_x_zz, t_y_y_y_xx, t_y_y_y_xy, t_y_y_y_xz, t_y_y_y_yy,\
                           t_y_y_y_yz, t_y_y_y_zz, t_y_y_z_xx, t_y_y_z_xy, t_y_y_z_xz,\
                           t_y_y_z_yy, t_y_y_z_yz, t_y_y_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_y_z_zz[i] = t_0_yy_z_zz[i] - rab_y[i] * t_0_y_z_zz[i];

        t_y_y_z_yz[i] = t_0_yy_z_yz[i] - rab_y[i] * t_0_y_z_yz[i];

        t_y_y_z_yy[i] = t_0_yy_z_yy[i] - rab_y[i] * t_0_y_z_yy[i];

        t_y_y_z_xz[i] = t_0_yy_z_xz[i] - rab_y[i] * t_0_y_z_xz[i];

        t_y_y_z_xy[i] = t_0_yy_z_xy[i] - rab_y[i] * t_0_y_z_xy[i];

        t_y_y_z_xx[i] = t_0_yy_z_xx[i] - rab_y[i] * t_0_y_z_xx[i];

        t_y_y_y_zz[i] = t_0_yy_y_zz[i] - rab_y[i] * t_0_y_y_zz[i];

        t_y_y_y_yz[i] = t_0_yy_y_yz[i] - rab_y[i] * t_0_y_y_yz[i];

        t_y_y_y_yy[i] = t_0_yy_y_yy[i] - rab_y[i] * t_0_y_y_yy[i];

        t_y_y_y_xz[i] = t_0_yy_y_xz[i] - rab_y[i] * t_0_y_y_xz[i];

        t_y_y_y_xy[i] = t_0_yy_y_xy[i] - rab_y[i] * t_0_y_y_xy[i];

        t_y_y_y_xx[i] = t_0_yy_y_xx[i] - rab_y[i] * t_0_y_y_xx[i];

        t_y_y_x_zz[i] = t_0_yy_x_zz[i] - rab_y[i] * t_0_y_x_zz[i];

        t_y_y_x_yz[i] = t_0_yy_x_yz[i] - rab_y[i] * t_0_y_x_yz[i];

        t_y_y_x_yy[i] = t_0_yy_x_yy[i] - rab_y[i] * t_0_y_x_yy[i];

        t_y_y_x_xz[i] = t_0_yy_x_xz[i] - rab_y[i] * t_0_y_x_xz[i];

        t_y_y_x_xy[i] = t_0_yy_x_xy[i] - rab_y[i] * t_0_y_x_xy[i];

        t_y_y_x_xx[i] = t_0_yy_x_xx[i] - rab_y[i] * t_0_y_x_xx[i];

        t_y_x_z_zz[i] = t_0_xy_z_zz[i] - rab_y[i] * t_0_x_z_zz[i];

        t_y_x_z_yz[i] = t_0_xy_z_yz[i] - rab_y[i] * t_0_x_z_yz[i];

        t_y_x_z_yy[i] = t_0_xy_z_yy[i] - rab_y[i] * t_0_x_z_yy[i];

        t_y_x_z_xz[i] = t_0_xy_z_xz[i] - rab_y[i] * t_0_x_z_xz[i];

        t_y_x_z_xy[i] = t_0_xy_z_xy[i] - rab_y[i] * t_0_x_z_xy[i];

        t_y_x_z_xx[i] = t_0_xy_z_xx[i] - rab_y[i] * t_0_x_z_xx[i];

        t_y_x_y_zz[i] = t_0_xy_y_zz[i] - rab_y[i] * t_0_x_y_zz[i];

        t_y_x_y_yz[i] = t_0_xy_y_yz[i] - rab_y[i] * t_0_x_y_yz[i];

        t_y_x_y_yy[i] = t_0_xy_y_yy[i] - rab_y[i] * t_0_x_y_yy[i];

        t_y_x_y_xz[i] = t_0_xy_y_xz[i] - rab_y[i] * t_0_x_y_xz[i];

        t_y_x_y_xy[i] = t_0_xy_y_xy[i] - rab_y[i] * t_0_x_y_xy[i];

        t_y_x_y_xx[i] = t_0_xy_y_xx[i] - rab_y[i] * t_0_x_y_xx[i];

        t_y_x_x_zz[i] = t_0_xy_x_zz[i] - rab_y[i] * t_0_x_x_zz[i];

        t_y_x_x_yz[i] = t_0_xy_x_yz[i] - rab_y[i] * t_0_x_x_yz[i];

        t_y_x_x_yy[i] = t_0_xy_x_yy[i] - rab_y[i] * t_0_x_x_yy[i];

        t_y_x_x_xz[i] = t_0_xy_x_xz[i] - rab_y[i] * t_0_x_x_xz[i];

        t_y_x_x_xy[i] = t_0_xy_x_xy[i] - rab_y[i] * t_0_x_x_xy[i];

        t_y_x_x_xx[i] = t_0_xy_x_xx[i] - rab_y[i] * t_0_x_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xy_x_xx, t_0_xy_x_xy, t_0_xy_x_xz, t_0_xy_x_yy,\
                           t_0_xy_x_yz, t_0_xy_x_zz, t_0_xy_y_xx, t_0_xy_y_xy, t_0_xy_y_xz,\
                           t_0_xy_y_yy, t_0_xy_y_yz, t_0_xy_y_zz, t_0_xy_z_xx, t_0_xy_z_xy,\
                           t_0_xy_z_xz, t_0_xy_z_yy, t_0_xy_z_yz, t_0_xy_z_zz, t_0_xz_x_xx,\
                           t_0_xz_x_xy, t_0_xz_x_xz, t_0_xz_x_yy, t_0_xz_x_yz, t_0_xz_x_zz,\
                           t_0_xz_y_xx, t_0_xz_y_xy, t_0_xz_y_xz, t_0_xz_y_yy, t_0_xz_y_yz,\
                           t_0_xz_y_zz, t_0_xz_z_xx, t_0_xz_z_xy, t_0_xz_z_xz, t_0_xz_z_yy,\
                           t_0_xz_z_yz, t_0_xz_z_zz, t_0_y_x_xx, t_0_y_x_xy, t_0_y_x_xz,\
                           t_0_y_x_yy, t_0_y_x_yz, t_0_y_x_zz, t_0_y_y_xx, t_0_y_y_xy,\
                           t_0_y_y_xz, t_0_y_y_yy, t_0_y_y_yz, t_0_y_y_zz, t_0_y_z_xx,\
                           t_0_y_z_xy, t_0_y_z_xz, t_0_y_z_yy, t_0_y_z_yz, t_0_y_z_zz,\
                           t_0_z_x_xx, t_0_z_x_xy, t_0_z_x_xz, t_0_z_x_yy, t_0_z_x_yz,\
                           t_0_z_x_zz, t_0_z_y_xx, t_0_z_y_xy, t_0_z_y_xz, t_0_z_y_yy,\
                           t_0_z_y_yz, t_0_z_y_zz, t_0_z_z_xx, t_0_z_z_xy, t_0_z_z_xz,\
                           t_0_z_z_yy, t_0_z_z_yz, t_0_z_z_zz, t_x_y_x_xx, t_x_y_x_xy,\
                           t_x_y_x_xz, t_x_y_x_yy, t_x_y_x_yz, t_x_y_x_zz, t_x_y_y_xx,\
                           t_x_y_y_xy, t_x_y_y_xz, t_x_y_y_yy, t_x_y_y_yz, t_x_y_y_zz,\
                           t_x_y_z_xx, t_x_y_z_xy, t_x_y_z_xz, t_x_y_z_yy, t_x_y_z_yz,\
                           t_x_y_z_zz, t_x_z_x_xx, t_x_z_x_xy, t_x_z_x_xz, t_x_z_x_yy,\
                           t_x_z_x_yz, t_x_z_x_zz, t_x_z_y_xx, t_x_z_y_xy, t_x_z_y_xz,\
                           t_x_z_y_yy, t_x_z_y_yz, t_x_z_y_zz, t_x_z_z_xx, t_x_z_z_xy,\
                           t_x_z_z_xz, t_x_z_z_yy, t_x_z_z_yz, t_x_z_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_z_z_zz[i] = t_0_xz_z_zz[i] - rab_x[i] * t_0_z_z_zz[i];

        t_x_z_z_yz[i] = t_0_xz_z_yz[i] - rab_x[i] * t_0_z_z_yz[i];

        t_x_z_z_yy[i] = t_0_xz_z_yy[i] - rab_x[i] * t_0_z_z_yy[i];

        t_x_z_z_xz[i] = t_0_xz_z_xz[i] - rab_x[i] * t_0_z_z_xz[i];

        t_x_z_z_xy[i] = t_0_xz_z_xy[i] - rab_x[i] * t_0_z_z_xy[i];

        t_x_z_z_xx[i] = t_0_xz_z_xx[i] - rab_x[i] * t_0_z_z_xx[i];

        t_x_z_y_zz[i] = t_0_xz_y_zz[i] - rab_x[i] * t_0_z_y_zz[i];

        t_x_z_y_yz[i] = t_0_xz_y_yz[i] - rab_x[i] * t_0_z_y_yz[i];

        t_x_z_y_yy[i] = t_0_xz_y_yy[i] - rab_x[i] * t_0_z_y_yy[i];

        t_x_z_y_xz[i] = t_0_xz_y_xz[i] - rab_x[i] * t_0_z_y_xz[i];

        t_x_z_y_xy[i] = t_0_xz_y_xy[i] - rab_x[i] * t_0_z_y_xy[i];

        t_x_z_y_xx[i] = t_0_xz_y_xx[i] - rab_x[i] * t_0_z_y_xx[i];

        t_x_z_x_zz[i] = t_0_xz_x_zz[i] - rab_x[i] * t_0_z_x_zz[i];

        t_x_z_x_yz[i] = t_0_xz_x_yz[i] - rab_x[i] * t_0_z_x_yz[i];

        t_x_z_x_yy[i] = t_0_xz_x_yy[i] - rab_x[i] * t_0_z_x_yy[i];

        t_x_z_x_xz[i] = t_0_xz_x_xz[i] - rab_x[i] * t_0_z_x_xz[i];

        t_x_z_x_xy[i] = t_0_xz_x_xy[i] - rab_x[i] * t_0_z_x_xy[i];

        t_x_z_x_xx[i] = t_0_xz_x_xx[i] - rab_x[i] * t_0_z_x_xx[i];

        t_x_y_z_zz[i] = t_0_xy_z_zz[i] - rab_x[i] * t_0_y_z_zz[i];

        t_x_y_z_yz[i] = t_0_xy_z_yz[i] - rab_x[i] * t_0_y_z_yz[i];

        t_x_y_z_yy[i] = t_0_xy_z_yy[i] - rab_x[i] * t_0_y_z_yy[i];

        t_x_y_z_xz[i] = t_0_xy_z_xz[i] - rab_x[i] * t_0_y_z_xz[i];

        t_x_y_z_xy[i] = t_0_xy_z_xy[i] - rab_x[i] * t_0_y_z_xy[i];

        t_x_y_z_xx[i] = t_0_xy_z_xx[i] - rab_x[i] * t_0_y_z_xx[i];

        t_x_y_y_zz[i] = t_0_xy_y_zz[i] - rab_x[i] * t_0_y_y_zz[i];

        t_x_y_y_yz[i] = t_0_xy_y_yz[i] - rab_x[i] * t_0_y_y_yz[i];

        t_x_y_y_yy[i] = t_0_xy_y_yy[i] - rab_x[i] * t_0_y_y_yy[i];

        t_x_y_y_xz[i] = t_0_xy_y_xz[i] - rab_x[i] * t_0_y_y_xz[i];

        t_x_y_y_xy[i] = t_0_xy_y_xy[i] - rab_x[i] * t_0_y_y_xy[i];

        t_x_y_y_xx[i] = t_0_xy_y_xx[i] - rab_x[i] * t_0_y_y_xx[i];

        t_x_y_x_zz[i] = t_0_xy_x_zz[i] - rab_x[i] * t_0_y_x_zz[i];

        t_x_y_x_yz[i] = t_0_xy_x_yz[i] - rab_x[i] * t_0_y_x_yz[i];

        t_x_y_x_yy[i] = t_0_xy_x_yy[i] - rab_x[i] * t_0_y_x_yy[i];

        t_x_y_x_xz[i] = t_0_xy_x_xz[i] - rab_x[i] * t_0_y_x_xz[i];

        t_x_y_x_xy[i] = t_0_xy_x_xy[i] - rab_x[i] * t_0_y_x_xy[i];

        t_x_y_x_xx[i] = t_0_xy_x_xx[i] - rab_x[i] * t_0_y_x_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_x_x_xx, t_0_x_x_xy, t_0_x_x_xz, t_0_x_x_yy, t_0_x_x_yz,\
                           t_0_x_x_zz, t_0_x_y_xx, t_0_x_y_xy, t_0_x_y_xz, t_0_x_y_yy,\
                           t_0_x_y_yz, t_0_x_y_zz, t_0_x_z_xx, t_0_x_z_xy, t_0_x_z_xz,\
                           t_0_x_z_yy, t_0_x_z_yz, t_0_x_z_zz, t_0_xx_x_xx, t_0_xx_x_xy,\
                           t_0_xx_x_xz, t_0_xx_x_yy, t_0_xx_x_yz, t_0_xx_x_zz, t_0_xx_y_xx,\
                           t_0_xx_y_xy, t_0_xx_y_xz, t_0_xx_y_yy, t_0_xx_y_yz, t_0_xx_y_zz,\
                           t_0_xx_z_xx, t_0_xx_z_xy, t_0_xx_z_xz, t_0_xx_z_yy, t_0_xx_z_yz,\
                           t_0_xx_z_zz, t_x_x_x_xx, t_x_x_x_xy, t_x_x_x_xz, t_x_x_x_yy,\
                           t_x_x_x_yz, t_x_x_x_zz, t_x_x_y_xx, t_x_x_y_xy, t_x_x_y_xz,\
                           t_x_x_y_yy, t_x_x_y_yz, t_x_x_y_zz, t_x_x_z_xx, t_x_x_z_xy,\
                           t_x_x_z_xz, t_x_x_z_yy, t_x_x_z_yz, t_x_x_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_x_z_zz[i] = t_0_xx_z_zz[i] - rab_x[i] * t_0_x_z_zz[i];

        t_x_x_z_yz[i] = t_0_xx_z_yz[i] - rab_x[i] * t_0_x_z_yz[i];

        t_x_x_z_yy[i] = t_0_xx_z_yy[i] - rab_x[i] * t_0_x_z_yy[i];

        t_x_x_z_xz[i] = t_0_xx_z_xz[i] - rab_x[i] * t_0_x_z_xz[i];

        t_x_x_z_xy[i] = t_0_xx_z_xy[i] - rab_x[i] * t_0_x_z_xy[i];

        t_x_x_z_xx[i] = t_0_xx_z_xx[i] - rab_x[i] * t_0_x_z_xx[i];

        t_x_x_y_zz[i] = t_0_xx_y_zz[i] - rab_x[i] * t_0_x_y_zz[i];

        t_x_x_y_yz[i] = t_0_xx_y_yz[i] - rab_x[i] * t_0_x_y_yz[i];

        t_x_x_y_yy[i] = t_0_xx_y_yy[i] - rab_x[i] * t_0_x_y_yy[i];

        t_x_x_y_xz[i] = t_0_xx_y_xz[i] - rab_x[i] * t_0_x_y_xz[i];

        t_x_x_y_xy[i] = t_0_xx_y_xy[i] - rab_x[i] * t_0_x_y_xy[i];

        t_x_x_y_xx[i] = t_0_xx_y_xx[i] - rab_x[i] * t_0_x_y_xx[i];

        t_x_x_x_zz[i] = t_0_xx_x_zz[i] - rab_x[i] * t_0_x_x_zz[i];

        t_x_x_x_yz[i] = t_0_xx_x_yz[i] - rab_x[i] * t_0_x_x_yz[i];

        t_x_x_x_yy[i] = t_0_xx_x_yy[i] - rab_x[i] * t_0_x_x_yy[i];

        t_x_x_x_xz[i] = t_0_xx_x_xz[i] - rab_x[i] * t_0_x_x_xz[i];

        t_x_x_x_xy[i] = t_0_xx_x_xy[i] - rab_x[i] * t_0_x_x_xy[i];

        t_x_x_x_xx[i] = t_0_xx_x_xx[i] - rab_x[i] * t_0_x_x_xx[i];
    }
}


} // derirec namespace
