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
compHostVRRForSGSG_V0(      BufferHostXY<T>&      intsBufferSGSG,
                      const BufferHostX<int32_t>& intsIndexesSGSG0,
                      const BufferHostXY<T>&      intsBufferSDSG0,
                      const BufferHostX<int32_t>& intsIndexesSDSG0,
                      const BufferHostXY<T>&      intsBufferSDSG1,
                      const BufferHostX<int32_t>& intsIndexesSDSG1,
                      const BufferHostXY<T>&      intsBufferSFSF1,
                      const BufferHostX<int32_t>& intsIndexesSFSF1,
                      const BufferHostXY<T>&      intsBufferSFSG0,
                      const BufferHostX<int32_t>& intsIndexesSFSG0,
                      const BufferHostXY<T>&      intsBufferSFSG1,
                      const BufferHostX<int32_t>& intsIndexesSFSG1,
                      const T*                    osFactorsZeta,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

    auto fz_0 = osFactorsBraZeta;

    auto frz2_0 = osFactorsBraRhoZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SGSG]^(0) integral components

    t_0_zzzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(0));

    t_0_zzzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(1));

    t_0_zzzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(2));

    t_0_zzzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(3));

    t_0_zzzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(4));

    t_0_zzzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(5));

    t_0_zzzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(6));

    t_0_zzzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(7));

    t_0_zzzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(8));

    t_0_zzzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(9));

    t_0_zzzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(10));

    t_0_zzzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(11));

    t_0_zzzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(12));

    t_0_zzzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(13));

    t_0_zzzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(14));

    t_0_yzzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(15));

    t_0_yzzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(16));

    t_0_yzzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(17));

    t_0_yzzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(18));

    t_0_yzzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(19));

    t_0_yzzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(20));

    t_0_yzzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(21));

    t_0_yzzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(22));

    t_0_yzzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(23));

    t_0_yzzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(24));

    t_0_yzzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(25));

    t_0_yzzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(26));

    t_0_yzzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(27));

    t_0_yzzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(28));

    t_0_yzzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(29));

    t_0_yyzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(30));

    t_0_yyzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(31));

    t_0_yyzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(32));

    t_0_yyzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(33));

    t_0_yyzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(34));

    t_0_yyzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(35));

    t_0_yyzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(36));

    t_0_yyzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(37));

    t_0_yyzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(38));

    t_0_yyzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(39));

    t_0_yyzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(40));

    t_0_yyzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(41));

    t_0_yyzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(42));

    t_0_yyzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(43));

    t_0_yyzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(44));

    t_0_yyyz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(45));

    t_0_yyyz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(46));

    t_0_yyyz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(47));

    t_0_yyyz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(48));

    t_0_yyyz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(49));

    t_0_yyyz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(50));

    t_0_yyyz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(51));

    t_0_yyyz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(52));

    t_0_yyyz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(53));

    t_0_yyyz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(54));

    t_0_yyyz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(55));

    t_0_yyyz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(56));

    t_0_yyyz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(57));

    t_0_yyyz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(58));

    t_0_yyyz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(59));

    t_0_yyyy_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(60));

    t_0_yyyy_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(61));

    t_0_yyyy_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(62));

    t_0_yyyy_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(63));

    t_0_yyyy_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(64));

    t_0_yyyy_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(65));

    t_0_yyyy_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(66));

    t_0_yyyy_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(67));

    t_0_yyyy_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(68));

    t_0_yyyy_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(69));

    t_0_yyyy_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(70));

    t_0_yyyy_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(71));

    t_0_yyyy_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(72));

    t_0_yyyy_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(73));

    t_0_yyyy_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(74));

    t_0_xzzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(75));

    t_0_xzzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(76));

    t_0_xzzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(77));

    t_0_xzzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(78));

    t_0_xzzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(79));

    t_0_xzzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(80));

    t_0_xzzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(81));

    t_0_xzzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(82));

    t_0_xzzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(83));

    t_0_xzzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(84));

    t_0_xzzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(85));

    t_0_xzzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(86));

    t_0_xzzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(87));

    t_0_xzzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(88));

    t_0_xzzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(89));

    t_0_xyzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(90));

    t_0_xyzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(91));

    t_0_xyzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(92));

    t_0_xyzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(93));

    t_0_xyzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(94));

    t_0_xyzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(95));

    t_0_xyzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(96));

    t_0_xyzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(97));

    t_0_xyzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(98));

    t_0_xyzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(99));

    t_0_xyzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(100));

    t_0_xyzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(101));

    t_0_xyzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(102));

    t_0_xyzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(103));

    t_0_xyzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(104));

    t_0_xyyz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(105));

    t_0_xyyz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(106));

    t_0_xyyz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(107));

    t_0_xyyz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(108));

    t_0_xyyz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(109));

    t_0_xyyz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(110));

    t_0_xyyz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(111));

    t_0_xyyz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(112));

    t_0_xyyz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(113));

    t_0_xyyz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(114));

    t_0_xyyz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(115));

    t_0_xyyz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(116));

    t_0_xyyz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(117));

    t_0_xyyz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(118));

    t_0_xyyz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(119));

    t_0_xyyy_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(120));

    t_0_xyyy_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(121));

    t_0_xyyy_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(122));

    t_0_xyyy_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(123));

    t_0_xyyy_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(124));

    t_0_xyyy_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(125));

    t_0_xyyy_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(126));

    t_0_xyyy_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(127));

    t_0_xyyy_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(128));

    t_0_xyyy_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(129));

    t_0_xyyy_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(130));

    t_0_xyyy_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(131));

    t_0_xyyy_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(132));

    t_0_xyyy_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(133));

    t_0_xyyy_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(134));

    t_0_xxzz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(135));

    t_0_xxzz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(136));

    t_0_xxzz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(137));

    t_0_xxzz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(138));

    t_0_xxzz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(139));

    t_0_xxzz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(140));

    t_0_xxzz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(141));

    t_0_xxzz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(142));

    t_0_xxzz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(143));

    t_0_xxzz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(144));

    t_0_xxzz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(145));

    t_0_xxzz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(146));

    t_0_xxzz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(147));

    t_0_xxzz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(148));

    t_0_xxzz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(149));

    t_0_xxyz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(150));

    t_0_xxyz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(151));

    t_0_xxyz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(152));

    t_0_xxyz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(153));

    t_0_xxyz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(154));

    t_0_xxyz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(155));

    t_0_xxyz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(156));

    t_0_xxyz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(157));

    t_0_xxyz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(158));

    t_0_xxyz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(159));

    t_0_xxyz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(160));

    t_0_xxyz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(161));

    t_0_xxyz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(162));

    t_0_xxyz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(163));

    t_0_xxyz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(164));

    t_0_xxyy_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(165));

    t_0_xxyy_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(166));

    t_0_xxyy_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(167));

    t_0_xxyy_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(168));

    t_0_xxyy_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(169));

    t_0_xxyy_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(170));

    t_0_xxyy_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(171));

    t_0_xxyy_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(172));

    t_0_xxyy_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(173));

    t_0_xxyy_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(174));

    t_0_xxyy_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(175));

    t_0_xxyy_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(176));

    t_0_xxyy_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(177));

    t_0_xxyy_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(178));

    t_0_xxyy_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(179));

    t_0_xxxz_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(180));

    t_0_xxxz_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(181));

    t_0_xxxz_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(182));

    t_0_xxxz_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(183));

    t_0_xxxz_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(184));

    t_0_xxxz_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(185));

    t_0_xxxz_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(186));

    t_0_xxxz_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(187));

    t_0_xxxz_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(188));

    t_0_xxxz_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(189));

    t_0_xxxz_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(190));

    t_0_xxxz_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(191));

    t_0_xxxz_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(192));

    t_0_xxxz_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(193));

    t_0_xxxz_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(194));

    t_0_xxxy_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(195));

    t_0_xxxy_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(196));

    t_0_xxxy_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(197));

    t_0_xxxy_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(198));

    t_0_xxxy_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(199));

    t_0_xxxy_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(200));

    t_0_xxxy_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(201));

    t_0_xxxy_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(202));

    t_0_xxxy_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(203));

    t_0_xxxy_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(204));

    t_0_xxxy_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(205));

    t_0_xxxy_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(206));

    t_0_xxxy_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(207));

    t_0_xxxy_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(208));

    t_0_xxxy_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(209));

    t_0_xxxx_0_zzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(210));

    t_0_xxxx_0_yzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(211));

    t_0_xxxx_0_yyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(212));

    t_0_xxxx_0_yyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(213));

    t_0_xxxx_0_yyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(214));

    t_0_xxxx_0_xzzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(215));

    t_0_xxxx_0_xyzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(216));

    t_0_xxxx_0_xyyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(217));

    t_0_xxxx_0_xyyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(218));

    t_0_xxxx_0_xxzz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(219));

    t_0_xxxx_0_xxyz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(220));

    t_0_xxxx_0_xxyy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(221));

    t_0_xxxx_0_xxxz_0 = intsBufferSGSG0.data(intsIndexesSGSG0(222));

    t_0_xxxx_0_xxxy_0 = intsBufferSGSG0.data(intsIndexesSGSG0(223));

    t_0_xxxx_0_xxxx_0 = intsBufferSGSG0.data(intsIndexesSGSG0(224));

    // set up [SDSG]^(0) integral components

    t_0_zz_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(0));

    t_0_zz_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(1));

    t_0_zz_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(2));

    t_0_zz_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(3));

    t_0_zz_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(4));

    t_0_zz_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(5));

    t_0_zz_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(6));

    t_0_zz_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(7));

    t_0_zz_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(8));

    t_0_zz_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(9));

    t_0_zz_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(10));

    t_0_zz_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(11));

    t_0_zz_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(12));

    t_0_zz_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(13));

    t_0_zz_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(14));

    t_0_yy_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(15));

    t_0_yy_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(16));

    t_0_yy_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(17));

    t_0_yy_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(18));

    t_0_yy_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(19));

    t_0_yy_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(20));

    t_0_yy_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(21));

    t_0_yy_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(22));

    t_0_yy_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(23));

    t_0_yy_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(24));

    t_0_yy_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(25));

    t_0_yy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(26));

    t_0_yy_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(27));

    t_0_yy_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(28));

    t_0_yy_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(29));

    t_0_xx_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(30));

    t_0_xx_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(31));

    t_0_xx_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(32));

    t_0_xx_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(33));

    t_0_xx_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(34));

    t_0_xx_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(35));

    t_0_xx_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(36));

    t_0_xx_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(37));

    t_0_xx_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(38));

    t_0_xx_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(39));

    t_0_xx_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(40));

    t_0_xx_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(41));

    t_0_xx_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(42));

    t_0_xx_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(43));

    t_0_xx_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(44));

    // set up [SDSG]^(1) integral components

    t_0_zz_0_zzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(0));

    t_0_zz_0_yzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(1));

    t_0_zz_0_yyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(2));

    t_0_zz_0_yyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(3));

    t_0_zz_0_yyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(4));

    t_0_zz_0_xzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(5));

    t_0_zz_0_xyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(6));

    t_0_zz_0_xyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(7));

    t_0_zz_0_xyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(8));

    t_0_zz_0_xxzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(9));

    t_0_zz_0_xxyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(10));

    t_0_zz_0_xxyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(11));

    t_0_zz_0_xxxz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(12));

    t_0_zz_0_xxxy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(13));

    t_0_zz_0_xxxx_1 = intsBufferSDSG1.data(intsIndexesSDSG1(14));

    t_0_yy_0_zzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(15));

    t_0_yy_0_yzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(16));

    t_0_yy_0_yyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(17));

    t_0_yy_0_yyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(18));

    t_0_yy_0_yyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(19));

    t_0_yy_0_xzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(20));

    t_0_yy_0_xyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(21));

    t_0_yy_0_xyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(22));

    t_0_yy_0_xyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(23));

    t_0_yy_0_xxzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(24));

    t_0_yy_0_xxyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(25));

    t_0_yy_0_xxyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(26));

    t_0_yy_0_xxxz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(27));

    t_0_yy_0_xxxy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(28));

    t_0_yy_0_xxxx_1 = intsBufferSDSG1.data(intsIndexesSDSG1(29));

    t_0_xx_0_zzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(30));

    t_0_xx_0_yzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(31));

    t_0_xx_0_yyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(32));

    t_0_xx_0_yyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(33));

    t_0_xx_0_yyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(34));

    t_0_xx_0_xzzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(35));

    t_0_xx_0_xyzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(36));

    t_0_xx_0_xyyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(37));

    t_0_xx_0_xyyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(38));

    t_0_xx_0_xxzz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(39));

    t_0_xx_0_xxyz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(40));

    t_0_xx_0_xxyy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(41));

    t_0_xx_0_xxxz_1 = intsBufferSDSG1.data(intsIndexesSDSG1(42));

    t_0_xx_0_xxxy_1 = intsBufferSDSG1.data(intsIndexesSDSG1(43));

    t_0_xx_0_xxxx_1 = intsBufferSDSG1.data(intsIndexesSDSG1(44));

    // set up [SFSF]^(1) integral components

    t_0_zzz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(0));

    t_0_zzz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(1));

    t_0_zzz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(2));

    t_0_zzz_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(3));

    t_0_zzz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(4));

    t_0_zzz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(5));

    t_0_zzz_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(6));

    t_0_zzz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(7));

    t_0_zzz_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(8));

    t_0_zzz_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(9));

    t_0_yzz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(10));

    t_0_yzz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(11));

    t_0_yzz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(12));

    t_0_yzz_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(13));

    t_0_yzz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(14));

    t_0_yzz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(15));

    t_0_yzz_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(16));

    t_0_yzz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(17));

    t_0_yzz_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(18));

    t_0_yyz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(19));

    t_0_yyz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(20));

    t_0_yyz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(21));

    t_0_yyz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(22));

    t_0_yyz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(23));

    t_0_yyz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(24));

    t_0_yyy_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(25));

    t_0_yyy_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(26));

    t_0_yyy_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(27));

    t_0_yyy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(28));

    t_0_yyy_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(29));

    t_0_yyy_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(30));

    t_0_yyy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(31));

    t_0_yyy_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(32));

    t_0_yyy_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(33));

    t_0_yyy_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(34));

    t_0_xzz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(35));

    t_0_xzz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(36));

    t_0_xzz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(37));

    t_0_xzz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(38));

    t_0_xzz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(39));

    t_0_xzz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(40));

    t_0_xyy_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(41));

    t_0_xyy_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(42));

    t_0_xyy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(43));

    t_0_xyy_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(44));

    t_0_xyy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(45));

    t_0_xyy_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(46));

    t_0_xxz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(47));

    t_0_xxz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(48));

    t_0_xxz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(49));

    t_0_xxz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(50));

    t_0_xxz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(51));

    t_0_xxz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(52));

    t_0_xxx_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(53));

    t_0_xxx_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(54));

    t_0_xxx_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(55));

    t_0_xxx_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(56));

    t_0_xxx_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(57));

    t_0_xxx_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(58));

    t_0_xxx_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(59));

    t_0_xxx_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(60));

    t_0_xxx_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(61));

    t_0_xxx_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(62));

    // set up [SFSG]^(0) integral components

    t_0_zzz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(0));

    t_0_zzz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(1));

    t_0_zzz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(2));

    t_0_zzz_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(3));

    t_0_zzz_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(4));

    t_0_zzz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(5));

    t_0_zzz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(6));

    t_0_zzz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(7));

    t_0_zzz_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(8));

    t_0_zzz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(9));

    t_0_zzz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(10));

    t_0_zzz_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(11));

    t_0_zzz_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(12));

    t_0_zzz_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(13));

    t_0_zzz_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(14));

    t_0_yzz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(15));

    t_0_yzz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(16));

    t_0_yzz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(17));

    t_0_yzz_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(18));

    t_0_yzz_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(19));

    t_0_yzz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(20));

    t_0_yzz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(21));

    t_0_yzz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(22));

    t_0_yzz_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(23));

    t_0_yzz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(24));

    t_0_yzz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(25));

    t_0_yzz_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(26));

    t_0_yzz_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(27));

    t_0_yzz_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(28));

    t_0_yzz_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(29));

    t_0_yyz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(30));

    t_0_yyz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(31));

    t_0_yyz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(32));

    t_0_yyz_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(33));

    t_0_yyz_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(34));

    t_0_yyz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(35));

    t_0_yyz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(36));

    t_0_yyz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(37));

    t_0_yyz_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(38));

    t_0_yyz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(39));

    t_0_yyz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(40));

    t_0_yyz_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(41));

    t_0_yyz_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(42));

    t_0_yyz_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(43));

    t_0_yyy_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(44));

    t_0_yyy_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(45));

    t_0_yyy_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(46));

    t_0_yyy_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(47));

    t_0_yyy_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(48));

    t_0_yyy_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(49));

    t_0_yyy_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(50));

    t_0_yyy_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(51));

    t_0_yyy_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(52));

    t_0_yyy_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(53));

    t_0_yyy_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(54));

    t_0_yyy_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(55));

    t_0_yyy_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(56));

    t_0_yyy_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(57));

    t_0_yyy_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(58));

    t_0_xzz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(59));

    t_0_xzz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(60));

    t_0_xzz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(61));

    t_0_xzz_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(62));

    t_0_xzz_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(63));

    t_0_xzz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(64));

    t_0_xzz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(65));

    t_0_xzz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(66));

    t_0_xzz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(67));

    t_0_xzz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(68));

    t_0_xzz_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(69));

    t_0_xzz_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(70));

    t_0_xyy_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(71));

    t_0_xyy_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(72));

    t_0_xyy_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(73));

    t_0_xyy_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(74));

    t_0_xyy_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(75));

    t_0_xyy_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(76));

    t_0_xyy_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(77));

    t_0_xyy_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(78));

    t_0_xyy_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(79));

    t_0_xyy_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(80));

    t_0_xyy_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(81));

    t_0_xyy_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(82));

    t_0_xxz_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(83));

    t_0_xxz_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(84));

    t_0_xxz_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(85));

    t_0_xxz_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(86));

    t_0_xxz_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(87));

    t_0_xxz_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(88));

    t_0_xxz_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(89));

    t_0_xxz_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(90));

    t_0_xxz_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(91));

    t_0_xxz_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(92));

    t_0_xxz_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(93));

    t_0_xxz_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(94));

    t_0_xxz_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(95));

    t_0_xxz_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(96));

    t_0_xxy_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(97));

    t_0_xxy_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(98));

    t_0_xxy_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(99));

    t_0_xxy_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(100));

    t_0_xxy_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(101));

    t_0_xxy_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(102));

    t_0_xxy_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(103));

    t_0_xxy_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(104));

    t_0_xxx_0_zzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(105));

    t_0_xxx_0_yzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(106));

    t_0_xxx_0_yyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(107));

    t_0_xxx_0_yyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(108));

    t_0_xxx_0_yyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(109));

    t_0_xxx_0_xzzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(110));

    t_0_xxx_0_xyzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(111));

    t_0_xxx_0_xyyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(112));

    t_0_xxx_0_xyyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(113));

    t_0_xxx_0_xxzz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(114));

    t_0_xxx_0_xxyz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(115));

    t_0_xxx_0_xxyy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(116));

    t_0_xxx_0_xxxz_0 = intsBufferSFSG0.data(intsIndexesSFSG0(117));

    t_0_xxx_0_xxxy_0 = intsBufferSFSG0.data(intsIndexesSFSG0(118));

    t_0_xxx_0_xxxx_0 = intsBufferSFSG0.data(intsIndexesSFSG0(119));

    // set up [SFSG]^(1) integral components

    t_0_zzz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(0));

    t_0_zzz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(1));

    t_0_zzz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(2));

    t_0_zzz_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(3));

    t_0_zzz_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(4));

    t_0_zzz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(5));

    t_0_zzz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(6));

    t_0_zzz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(7));

    t_0_zzz_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(8));

    t_0_zzz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(9));

    t_0_zzz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(10));

    t_0_zzz_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(11));

    t_0_zzz_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(12));

    t_0_zzz_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(13));

    t_0_zzz_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(14));

    t_0_yzz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(15));

    t_0_yzz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(16));

    t_0_yzz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(17));

    t_0_yzz_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(18));

    t_0_yzz_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(19));

    t_0_yzz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(20));

    t_0_yzz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(21));

    t_0_yzz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(22));

    t_0_yzz_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(23));

    t_0_yzz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(24));

    t_0_yzz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(25));

    t_0_yzz_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(26));

    t_0_yzz_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(27));

    t_0_yzz_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(28));

    t_0_yzz_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(29));

    t_0_yyz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(30));

    t_0_yyz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(31));

    t_0_yyz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(32));

    t_0_yyz_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(33));

    t_0_yyz_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(34));

    t_0_yyz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(35));

    t_0_yyz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(36));

    t_0_yyz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(37));

    t_0_yyz_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(38));

    t_0_yyz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(39));

    t_0_yyz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(40));

    t_0_yyz_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(41));

    t_0_yyz_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(42));

    t_0_yyz_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(43));

    t_0_yyy_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(44));

    t_0_yyy_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(45));

    t_0_yyy_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(46));

    t_0_yyy_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(47));

    t_0_yyy_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(48));

    t_0_yyy_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(49));

    t_0_yyy_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(50));

    t_0_yyy_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(51));

    t_0_yyy_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(52));

    t_0_yyy_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(53));

    t_0_yyy_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(54));

    t_0_yyy_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(55));

    t_0_yyy_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(56));

    t_0_yyy_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(57));

    t_0_yyy_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(58));

    t_0_xzz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(59));

    t_0_xzz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(60));

    t_0_xzz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(61));

    t_0_xzz_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(62));

    t_0_xzz_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(63));

    t_0_xzz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(64));

    t_0_xzz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(65));

    t_0_xzz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(66));

    t_0_xzz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(67));

    t_0_xzz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(68));

    t_0_xzz_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(69));

    t_0_xzz_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(70));

    t_0_xyy_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(71));

    t_0_xyy_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(72));

    t_0_xyy_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(73));

    t_0_xyy_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(74));

    t_0_xyy_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(75));

    t_0_xyy_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(76));

    t_0_xyy_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(77));

    t_0_xyy_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(78));

    t_0_xyy_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(79));

    t_0_xyy_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(80));

    t_0_xyy_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(81));

    t_0_xyy_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(82));

    t_0_xxz_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(83));

    t_0_xxz_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(84));

    t_0_xxz_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(85));

    t_0_xxz_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(86));

    t_0_xxz_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(87));

    t_0_xxz_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(88));

    t_0_xxz_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(89));

    t_0_xxz_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(90));

    t_0_xxz_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(91));

    t_0_xxz_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(92));

    t_0_xxz_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(93));

    t_0_xxz_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(94));

    t_0_xxz_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(95));

    t_0_xxz_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(96));

    t_0_xxy_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(97));

    t_0_xxy_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(98));

    t_0_xxy_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(99));

    t_0_xxy_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(100));

    t_0_xxy_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(101));

    t_0_xxy_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(102));

    t_0_xxy_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(103));

    t_0_xxy_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(104));

    t_0_xxx_0_zzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(105));

    t_0_xxx_0_yzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(106));

    t_0_xxx_0_yyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(107));

    t_0_xxx_0_yyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(108));

    t_0_xxx_0_yyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(109));

    t_0_xxx_0_xzzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(110));

    t_0_xxx_0_xyzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(111));

    t_0_xxx_0_xyyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(112));

    t_0_xxx_0_xyyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(113));

    t_0_xxx_0_xxzz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(114));

    t_0_xxx_0_xxyz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(115));

    t_0_xxx_0_xxyy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(116));

    t_0_xxx_0_xxxz_1 = intsBufferSFSG1.data(intsIndexesSFSG1(117));

    t_0_xxx_0_xxxy_1 = intsBufferSFSG1.data(intsIndexesSFSG1(118));

    t_0_xxx_0_xxxx_1 = intsBufferSFSG1.data(intsIndexesSFSG1(119));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_yyyy_0,\
                               t_0_yy_0_yyyy_1, t_0_yyz_0_yyyy_0, t_0_yyz_0_yyyy_1,\
                               t_0_yyzz_0_xzzz_0, t_0_yyzz_0_yyyy_0, t_0_yyzz_0_yyyz_0,\
                               t_0_yyzz_0_yyzz_0, t_0_yyzz_0_yzzz_0, t_0_yyzz_0_zzzz_0,\
                               t_0_yzz_0_xzzz_0, t_0_yzz_0_xzzz_1, t_0_yzz_0_yyyz_0,\
                               t_0_yzz_0_yyyz_1, t_0_yzz_0_yyz_1, t_0_yzz_0_yyzz_0,\
                               t_0_yzz_0_yyzz_1, t_0_yzz_0_yzz_1, t_0_yzz_0_yzzz_0,\
                               t_0_yzz_0_yzzz_1, t_0_yzz_0_zzz_1, t_0_yzz_0_zzzz_0,\
                               t_0_yzz_0_zzzz_1, t_0_yzzz_0_xxxx_0, t_0_yzzz_0_xxxy_0,\
                               t_0_yzzz_0_xxxz_0, t_0_yzzz_0_xxyy_0, t_0_yzzz_0_xxyz_0,\
                               t_0_yzzz_0_xxzz_0, t_0_yzzz_0_xyyy_0, t_0_yzzz_0_xyyz_0,\
                               t_0_yzzz_0_xyzz_0, t_0_yzzz_0_xzzz_0, t_0_yzzz_0_yyyy_0,\
                               t_0_yzzz_0_yyyz_0, t_0_yzzz_0_yyzz_0, t_0_yzzz_0_yzzz_0,\
                               t_0_yzzz_0_zzzz_0, t_0_zz_0_xxxx_0, t_0_zz_0_xxxx_1,\
                               t_0_zz_0_xxxy_0, t_0_zz_0_xxxy_1, t_0_zz_0_xxxz_0, t_0_zz_0_xxxz_1,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyy_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1,\
                               t_0_zz_0_xxzz_0, t_0_zz_0_xxzz_1, t_0_zz_0_xyyy_0, t_0_zz_0_xyyy_1,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1, t_0_zz_0_xyzz_0, t_0_zz_0_xyzz_1,\
                               t_0_zz_0_xzzz_0, t_0_zz_0_xzzz_1, t_0_zz_0_yyyy_0, t_0_zz_0_yyyy_1,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyyz_1, t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1,\
                               t_0_zz_0_yzzz_0, t_0_zz_0_yzzz_1, t_0_zz_0_zzzz_0, t_0_zz_0_zzzz_1,\
                               t_0_zzz_0_xxx_1, t_0_zzz_0_xxxx_0, t_0_zzz_0_xxxx_1,\
                               t_0_zzz_0_xxxy_0, t_0_zzz_0_xxxy_1, t_0_zzz_0_xxxz_0,\
                               t_0_zzz_0_xxxz_1, t_0_zzz_0_xxy_1, t_0_zzz_0_xxyy_0,\
                               t_0_zzz_0_xxyy_1, t_0_zzz_0_xxyz_0, t_0_zzz_0_xxyz_1,\
                               t_0_zzz_0_xxz_1, t_0_zzz_0_xxzz_0, t_0_zzz_0_xxzz_1,\
                               t_0_zzz_0_xyy_1, t_0_zzz_0_xyyy_0, t_0_zzz_0_xyyy_1,\
                               t_0_zzz_0_xyyz_0, t_0_zzz_0_xyyz_1, t_0_zzz_0_xyz_1,\
                               t_0_zzz_0_xyzz_0, t_0_zzz_0_xyzz_1, t_0_zzz_0_xzz_1,\
                               t_0_zzz_0_xzzz_0, t_0_zzz_0_xzzz_1, t_0_zzz_0_yyy_1,\
                               t_0_zzz_0_yyyy_0, t_0_zzz_0_yyyy_1, t_0_zzz_0_yyyz_0,\
                               t_0_zzz_0_yyyz_1, t_0_zzz_0_yyz_1, t_0_zzz_0_yyzz_0,\
                               t_0_zzz_0_yyzz_1, t_0_zzz_0_yzz_1, t_0_zzz_0_yzzz_0,\
                               t_0_zzz_0_yzzz_1, t_0_zzz_0_zzz_1, t_0_zzz_0_zzzz_0,\
                               t_0_zzz_0_zzzz_1, t_0_zzzz_0_xxxx_0, t_0_zzzz_0_xxxy_0,\
                               t_0_zzzz_0_xxxz_0, t_0_zzzz_0_xxyy_0, t_0_zzzz_0_xxyz_0,\
                               t_0_zzzz_0_xxzz_0, t_0_zzzz_0_xyyy_0, t_0_zzzz_0_xyyz_0,\
                               t_0_zzzz_0_xyzz_0, t_0_zzzz_0_xzzz_0, t_0_zzzz_0_yyyy_0,\
                               t_0_zzzz_0_yyyz_0, t_0_zzzz_0_yyzz_0, t_0_zzzz_0_yzzz_0,\
                               t_0_zzzz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzzz_0[i] += rpb_z[i] * t_0_zzz_0_zzzz_0[i] + rwp_z[i] * t_0_zzz_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_zzzz_0_yzzz_0[i] += rpb_z[i] * t_0_zzz_0_yzzz_0[i] + rwp_z[i] * t_0_zzz_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_zzzz_0_yyzz_0[i] += rpb_z[i] * t_0_zzz_0_yyzz_0[i] + rwp_z[i] * t_0_zzz_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_zzzz_0_yyyz_0[i] += rpb_z[i] * t_0_zzz_0_yyyz_0[i] + rwp_z[i] * t_0_zzz_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_zzzz_0_yyyy_0[i] += rpb_z[i] * t_0_zzz_0_yyyy_0[i] + rwp_z[i] * t_0_zzz_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyyy_1[i];

            t_0_zzzz_0_xzzz_0[i] += rpb_z[i] * t_0_zzz_0_xzzz_0[i] + rwp_z[i] * t_0_zzz_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_zzzz_0_xyzz_0[i] += rpb_z[i] * t_0_zzz_0_xyzz_0[i] + rwp_z[i] * t_0_zzz_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_zzzz_0_xyyz_0[i] += rpb_z[i] * t_0_zzz_0_xyyz_0[i] + rwp_z[i] * t_0_zzz_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_zzzz_0_xyyy_0[i] += rpb_z[i] * t_0_zzz_0_xyyy_0[i] + rwp_z[i] * t_0_zzz_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyyy_1[i];

            t_0_zzzz_0_xxzz_0[i] += rpb_z[i] * t_0_zzz_0_xxzz_0[i] + rwp_z[i] * t_0_zzz_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_zzzz_0_xxyz_0[i] += rpb_z[i] * t_0_zzz_0_xxyz_0[i] + rwp_z[i] * t_0_zzz_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_zzzz_0_xxyy_0[i] += rpb_z[i] * t_0_zzz_0_xxyy_0[i] + rwp_z[i] * t_0_zzz_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxyy_1[i];

            t_0_zzzz_0_xxxz_0[i] += rpb_z[i] * t_0_zzz_0_xxxz_0[i] + rwp_z[i] * t_0_zzz_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_zzzz_0_xxxy_0[i] += rpb_z[i] * t_0_zzz_0_xxxy_0[i] + rwp_z[i] * t_0_zzz_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxy_1[i];

            t_0_zzzz_0_xxxx_0[i] += rpb_z[i] * t_0_zzz_0_xxxx_0[i] + rwp_z[i] * t_0_zzz_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxx_1[i];

            t_0_yzzz_0_zzzz_0[i] += rpb_y[i] * t_0_zzz_0_zzzz_0[i] + rwp_y[i] * t_0_zzz_0_zzzz_1[i];

            t_0_yzzz_0_yzzz_0[i] += rpb_y[i] * t_0_zzz_0_yzzz_0[i] + rwp_y[i] * t_0_zzz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yyzz_0[i] += rpb_y[i] * t_0_zzz_0_yyzz_0[i] + rwp_y[i] * t_0_zzz_0_yyzz_1[i] + fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_yzzz_0_yyyz_0[i] += rpb_y[i] * t_0_zzz_0_yyyz_0[i] + rwp_y[i] * t_0_zzz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_yzzz_0_yyyy_0[i] += rpb_y[i] * t_0_zzz_0_yyyy_0[i] + rwp_y[i] * t_0_zzz_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_yzzz_0_xzzz_0[i] += rpb_y[i] * t_0_zzz_0_xzzz_0[i] + rwp_y[i] * t_0_zzz_0_xzzz_1[i];

            t_0_yzzz_0_xyzz_0[i] += rpb_y[i] * t_0_zzz_0_xyzz_0[i] + rwp_y[i] * t_0_zzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_yzzz_0_xyyz_0[i] += rpb_y[i] * t_0_zzz_0_xyyz_0[i] + rwp_y[i] * t_0_zzz_0_xyyz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_yzzz_0_xyyy_0[i] += rpb_y[i] * t_0_zzz_0_xyyy_0[i] + rwp_y[i] * t_0_zzz_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_yzzz_0_xxzz_0[i] += rpb_y[i] * t_0_zzz_0_xxzz_0[i] + rwp_y[i] * t_0_zzz_0_xxzz_1[i];

            t_0_yzzz_0_xxyz_0[i] += rpb_y[i] * t_0_zzz_0_xxyz_0[i] + rwp_y[i] * t_0_zzz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_yzzz_0_xxyy_0[i] += rpb_y[i] * t_0_zzz_0_xxyy_0[i] + rwp_y[i] * t_0_zzz_0_xxyy_1[i] + fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_yzzz_0_xxxz_0[i] += rpb_y[i] * t_0_zzz_0_xxxz_0[i] + rwp_y[i] * t_0_zzz_0_xxxz_1[i];

            t_0_yzzz_0_xxxy_0[i] += rpb_y[i] * t_0_zzz_0_xxxy_0[i] + rwp_y[i] * t_0_zzz_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_yzzz_0_xxxx_0[i] += rpb_y[i] * t_0_zzz_0_xxxx_0[i] + rwp_y[i] * t_0_zzz_0_xxxx_1[i];

            t_0_yyzz_0_zzzz_0[i] += rpb_y[i] * t_0_yzz_0_zzzz_0[i] + rwp_y[i] * t_0_yzz_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i];

            t_0_yyzz_0_yzzz_0[i] += rpb_y[i] * t_0_yzz_0_yzzz_0[i] + rwp_y[i] * t_0_yzz_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_zzz_1[i];

            t_0_yyzz_0_yyzz_0[i] += rpb_y[i] * t_0_yzz_0_yyzz_0[i] + rwp_y[i] * t_0_yzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_yyzz_0_yyyz_0[i] += rpb_y[i] * t_0_yzz_0_yyyz_0[i] + rwp_y[i] * t_0_yzz_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_yzz_0_yyz_1[i];

            t_0_yyzz_0_yyyy_0[i] += rpb_z[i] * t_0_yyz_0_yyyy_0[i] + rwp_z[i] * t_0_yyz_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i];

            t_0_yyzz_0_xzzz_0[i] += rpb_y[i] * t_0_yzz_0_xzzz_0[i] + rwp_y[i] * t_0_yzz_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_xxxy_0,\
                               t_0_yy_0_xxxy_1, t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1, t_0_yy_0_xxyz_0,\
                               t_0_yy_0_xxyz_1, t_0_yy_0_xxzz_0, t_0_yy_0_xxzz_1, t_0_yy_0_xyyy_0,\
                               t_0_yy_0_xyyy_1, t_0_yy_0_xyyz_0, t_0_yy_0_xyyz_1, t_0_yy_0_xyzz_0,\
                               t_0_yy_0_xyzz_1, t_0_yy_0_xzzz_0, t_0_yy_0_xzzz_1, t_0_yy_0_yyyy_0,\
                               t_0_yy_0_yyyy_1, t_0_yy_0_yyyz_0, t_0_yy_0_yyyz_1, t_0_yy_0_yyzz_0,\
                               t_0_yy_0_yyzz_1, t_0_yy_0_yzzz_0, t_0_yy_0_yzzz_1, t_0_yy_0_zzzz_0,\
                               t_0_yy_0_zzzz_1, t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyy_0_xxy_1, t_0_yyy_0_xxyy_0,\
                               t_0_yyy_0_xxyy_1, t_0_yyy_0_xxyz_0, t_0_yyy_0_xxyz_1,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xxzz_0, t_0_yyy_0_xxzz_1,\
                               t_0_yyy_0_xyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_xyyz_0, t_0_yyy_0_xyyz_1, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xyzz_0, t_0_yyy_0_xyzz_1, t_0_yyy_0_xzz_1,\
                               t_0_yyy_0_xzzz_0, t_0_yyy_0_xzzz_1, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1, t_0_yyy_0_yyyz_0,\
                               t_0_yyy_0_yyyz_1, t_0_yyy_0_yyz_1, t_0_yyy_0_yyzz_0,\
                               t_0_yyy_0_yyzz_1, t_0_yyy_0_yzz_1, t_0_yyy_0_yzzz_0,\
                               t_0_yyy_0_yzzz_1, t_0_yyy_0_zzz_1, t_0_yyy_0_zzzz_0,\
                               t_0_yyy_0_zzzz_1, t_0_yyyy_0_xxyy_0, t_0_yyyy_0_xxyz_0,\
                               t_0_yyyy_0_xxzz_0, t_0_yyyy_0_xyyy_0, t_0_yyyy_0_xyyz_0,\
                               t_0_yyyy_0_xyzz_0, t_0_yyyy_0_xzzz_0, t_0_yyyy_0_yyyy_0,\
                               t_0_yyyy_0_yyyz_0, t_0_yyyy_0_yyzz_0, t_0_yyyy_0_yzzz_0,\
                               t_0_yyyy_0_zzzz_0, t_0_yyyz_0_xxxx_0, t_0_yyyz_0_xxxy_0,\
                               t_0_yyyz_0_xxxz_0, t_0_yyyz_0_xxyy_0, t_0_yyyz_0_xxyz_0,\
                               t_0_yyyz_0_xxzz_0, t_0_yyyz_0_xyyy_0, t_0_yyyz_0_xyyz_0,\
                               t_0_yyyz_0_xyzz_0, t_0_yyyz_0_xzzz_0, t_0_yyyz_0_yyyy_0,\
                               t_0_yyyz_0_yyyz_0, t_0_yyyz_0_yyzz_0, t_0_yyyz_0_yzzz_0,\
                               t_0_yyyz_0_zzzz_0, t_0_yyz_0_xxxy_0, t_0_yyz_0_xxxy_1,\
                               t_0_yyz_0_xxyy_0, t_0_yyz_0_xxyy_1, t_0_yyz_0_xyyy_0,\
                               t_0_yyz_0_xyyy_1, t_0_yyzz_0_xxxx_0, t_0_yyzz_0_xxxy_0,\
                               t_0_yyzz_0_xxxz_0, t_0_yyzz_0_xxyy_0, t_0_yyzz_0_xxyz_0,\
                               t_0_yyzz_0_xxzz_0, t_0_yyzz_0_xyyy_0, t_0_yyzz_0_xyyz_0,\
                               t_0_yyzz_0_xyzz_0, t_0_yzz_0_xxxx_0, t_0_yzz_0_xxxx_1,\
                               t_0_yzz_0_xxxz_0, t_0_yzz_0_xxxz_1, t_0_yzz_0_xxyz_0,\
                               t_0_yzz_0_xxyz_1, t_0_yzz_0_xxz_1, t_0_yzz_0_xxzz_0,\
                               t_0_yzz_0_xxzz_1, t_0_yzz_0_xyyz_0, t_0_yzz_0_xyyz_1,\
                               t_0_yzz_0_xyz_1, t_0_yzz_0_xyzz_0, t_0_yzz_0_xyzz_1,\
                               t_0_yzz_0_xzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxx_1, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxxz_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1, t_0_zz_0_xxzz_0,\
                               t_0_zz_0_xxzz_1, t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1, t_0_zz_0_xyzz_0,\
                               t_0_zz_0_xyzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyzz_0_xyzz_0[i] += rpb_y[i] * t_0_yzz_0_xyzz_0[i] + rwp_y[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xzz_1[i];

            t_0_yyzz_0_xyyz_0[i] += rpb_y[i] * t_0_yzz_0_xyyz_0[i] + rwp_y[i] * t_0_yzz_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fze_0[i] * t_0_yzz_0_xyz_1[i];

            t_0_yyzz_0_xyyy_0[i] += rpb_z[i] * t_0_yyz_0_xyyy_0[i] + rwp_z[i] * t_0_yyz_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i];

            t_0_yyzz_0_xxzz_0[i] += rpb_y[i] * t_0_yzz_0_xxzz_0[i] + rwp_y[i] * t_0_yzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i];

            t_0_yyzz_0_xxyz_0[i] += rpb_y[i] * t_0_yzz_0_xxyz_0[i] + rwp_y[i] * t_0_yzz_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xxz_1[i];

            t_0_yyzz_0_xxyy_0[i] += rpb_z[i] * t_0_yyz_0_xxyy_0[i] + rwp_z[i] * t_0_yyz_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i];

            t_0_yyzz_0_xxxz_0[i] += rpb_y[i] * t_0_yzz_0_xxxz_0[i] + rwp_y[i] * t_0_yzz_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i];

            t_0_yyzz_0_xxxy_0[i] += rpb_z[i] * t_0_yyz_0_xxxy_0[i] + rwp_z[i] * t_0_yyz_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i];

            t_0_yyzz_0_xxxx_0[i] += rpb_y[i] * t_0_yzz_0_xxxx_0[i] + rwp_y[i] * t_0_yzz_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxx_1[i];

            t_0_yyyz_0_zzzz_0[i] += rpb_z[i] * t_0_yyy_0_zzzz_0[i] + rwp_z[i] * t_0_yyy_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_yyyz_0_yzzz_0[i] += rpb_z[i] * t_0_yyy_0_yzzz_0[i] + rwp_z[i] * t_0_yyy_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_yyyz_0_yyzz_0[i] += rpb_z[i] * t_0_yyy_0_yyzz_0[i] + rwp_z[i] * t_0_yyy_0_yyzz_1[i] + fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_yyyz_0_yyyz_0[i] += rpb_z[i] * t_0_yyy_0_yyyz_0[i] + rwp_z[i] * t_0_yyy_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyz_0_yyyy_0[i] += rpb_z[i] * t_0_yyy_0_yyyy_0[i] + rwp_z[i] * t_0_yyy_0_yyyy_1[i];

            t_0_yyyz_0_xzzz_0[i] += rpb_z[i] * t_0_yyy_0_xzzz_0[i] + rwp_z[i] * t_0_yyy_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_yyyz_0_xyzz_0[i] += rpb_z[i] * t_0_yyy_0_xyzz_0[i] + rwp_z[i] * t_0_yyy_0_xyzz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_yyyz_0_xyyz_0[i] += rpb_z[i] * t_0_yyy_0_xyyz_0[i] + rwp_z[i] * t_0_yyy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyz_0_xyyy_0[i] += rpb_z[i] * t_0_yyy_0_xyyy_0[i] + rwp_z[i] * t_0_yyy_0_xyyy_1[i];

            t_0_yyyz_0_xxzz_0[i] += rpb_z[i] * t_0_yyy_0_xxzz_0[i] + rwp_z[i] * t_0_yyy_0_xxzz_1[i] + fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_yyyz_0_xxyz_0[i] += rpb_z[i] * t_0_yyy_0_xxyz_0[i] + rwp_z[i] * t_0_yyy_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxy_1[i];

            t_0_yyyz_0_xxyy_0[i] += rpb_z[i] * t_0_yyy_0_xxyy_0[i] + rwp_z[i] * t_0_yyy_0_xxyy_1[i];

            t_0_yyyz_0_xxxz_0[i] += rpb_z[i] * t_0_yyy_0_xxxz_0[i] + rwp_z[i] * t_0_yyy_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyz_0_xxxy_0[i] += rpb_z[i] * t_0_yyy_0_xxxy_0[i] + rwp_z[i] * t_0_yyy_0_xxxy_1[i];

            t_0_yyyz_0_xxxx_0[i] += rpb_z[i] * t_0_yyy_0_xxxx_0[i] + rwp_z[i] * t_0_yyy_0_xxxx_1[i];

            t_0_yyyy_0_zzzz_0[i] += rpb_y[i] * t_0_yyy_0_zzzz_0[i] + rwp_y[i] * t_0_yyy_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zzzz_1[i];

            t_0_yyyy_0_yzzz_0[i] += rpb_y[i] * t_0_yyy_0_yzzz_0[i] + rwp_y[i] * t_0_yyy_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_yyyy_0_yyzz_0[i] += rpb_y[i] * t_0_yyy_0_yyzz_0[i] + rwp_y[i] * t_0_yyy_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyzz_1[i] + fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_yyyy_0_yyyz_0[i] += rpb_y[i] * t_0_yyy_0_yyyz_0[i] + rwp_y[i] * t_0_yyy_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_yyyy_0_yyyy_0[i] += rpb_y[i] * t_0_yyy_0_yyyy_0[i] + rwp_y[i] * t_0_yyy_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyy_0_xzzz_0[i] += rpb_y[i] * t_0_yyy_0_xzzz_0[i] + rwp_y[i] * t_0_yyy_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xzzz_1[i];

            t_0_yyyy_0_xyzz_0[i] += rpb_y[i] * t_0_yyy_0_xyzz_0[i] + rwp_y[i] * t_0_yyy_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_yyyy_0_xyyz_0[i] += rpb_y[i] * t_0_yyy_0_xyyz_0[i] + rwp_y[i] * t_0_yyy_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyyz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_yyyy_0_xyyy_0[i] += rpb_y[i] * t_0_yyy_0_xyyy_0[i] + rwp_y[i] * t_0_yyy_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyy_0_xxzz_0[i] += rpb_y[i] * t_0_yyy_0_xxzz_0[i] + rwp_y[i] * t_0_yyy_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxzz_1[i];

            t_0_yyyy_0_xxyz_0[i] += rpb_y[i] * t_0_yyy_0_xxyz_0[i] + rwp_y[i] * t_0_yyy_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_yyyy_0_xxyy_0[i] += rpb_y[i] * t_0_yyy_0_xxyy_0[i] + rwp_y[i] * t_0_yyy_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_yyy_0_xxy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rwp_x, rwp_y, t_0_xyyz_0_yyzz_0,\
                               t_0_xyyz_0_yzzz_0, t_0_xyyz_0_zzzz_0, t_0_xyzz_0_xxxx_0,\
                               t_0_xyzz_0_xxxy_0, t_0_xyzz_0_xxxz_0, t_0_xyzz_0_xxyy_0,\
                               t_0_xyzz_0_xxyz_0, t_0_xyzz_0_xxzz_0, t_0_xyzz_0_xyyy_0,\
                               t_0_xyzz_0_xyyz_0, t_0_xyzz_0_xyzz_0, t_0_xyzz_0_xzzz_0,\
                               t_0_xyzz_0_yyyy_0, t_0_xyzz_0_yyyz_0, t_0_xyzz_0_yyzz_0,\
                               t_0_xyzz_0_yzzz_0, t_0_xyzz_0_zzzz_0, t_0_xzz_0_xxxx_0,\
                               t_0_xzz_0_xxxx_1, t_0_xzz_0_xxxz_0, t_0_xzz_0_xxxz_1,\
                               t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1, t_0_xzz_0_xzzz_0,\
                               t_0_xzz_0_xzzz_1, t_0_xzzz_0_xxxx_0, t_0_xzzz_0_xxxy_0,\
                               t_0_xzzz_0_xxxz_0, t_0_xzzz_0_xxyy_0, t_0_xzzz_0_xxyz_0,\
                               t_0_xzzz_0_xxzz_0, t_0_xzzz_0_xyyy_0, t_0_xzzz_0_xyyz_0,\
                               t_0_xzzz_0_xyzz_0, t_0_xzzz_0_xzzz_0, t_0_xzzz_0_yyyy_0,\
                               t_0_xzzz_0_yyyz_0, t_0_xzzz_0_yyzz_0, t_0_xzzz_0_yzzz_0,\
                               t_0_xzzz_0_zzzz_0, t_0_yy_0_xxxx_0, t_0_yy_0_xxxx_1,\
                               t_0_yy_0_xxxy_0, t_0_yy_0_xxxy_1, t_0_yy_0_xxxz_0, t_0_yy_0_xxxz_1,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyyy_0_xxxx_0, t_0_yyyy_0_xxxy_0,\
                               t_0_yyyy_0_xxxz_0, t_0_yyz_0_yyzz_0, t_0_yyz_0_yyzz_1,\
                               t_0_yyz_0_yzzz_0, t_0_yyz_0_yzzz_1, t_0_yyz_0_zzzz_0,\
                               t_0_yyz_0_zzzz_1, t_0_yzz_0_xxxy_0, t_0_yzz_0_xxxy_1,\
                               t_0_yzz_0_xxy_1, t_0_yzz_0_xxyy_0, t_0_yzz_0_xxyy_1,\
                               t_0_yzz_0_xxyz_0, t_0_yzz_0_xxyz_1, t_0_yzz_0_xyy_1,\
                               t_0_yzz_0_xyyy_0, t_0_yzz_0_xyyy_1, t_0_yzz_0_xyyz_0,\
                               t_0_yzz_0_xyyz_1, t_0_yzz_0_xyz_1, t_0_yzz_0_xyzz_0,\
                               t_0_yzz_0_xyzz_1, t_0_yzz_0_yyy_1, t_0_yzz_0_yyyy_0,\
                               t_0_yzz_0_yyyy_1, t_0_yzz_0_yyyz_0, t_0_yzz_0_yyyz_1,\
                               t_0_yzz_0_yyz_1, t_0_yzz_0_yyzz_0, t_0_yzz_0_yyzz_1,\
                               t_0_yzz_0_yzz_1, t_0_yzz_0_yzzz_0, t_0_yzz_0_yzzz_1,\
                               t_0_yzz_0_zzzz_0, t_0_yzz_0_zzzz_1, t_0_zzz_0_xxx_1,\
                               t_0_zzz_0_xxxx_0, t_0_zzz_0_xxxx_1, t_0_zzz_0_xxxy_0,\
                               t_0_zzz_0_xxxy_1, t_0_zzz_0_xxxz_0, t_0_zzz_0_xxxz_1,\
                               t_0_zzz_0_xxy_1, t_0_zzz_0_xxyy_0, t_0_zzz_0_xxyy_1,\
                               t_0_zzz_0_xxyz_0, t_0_zzz_0_xxyz_1, t_0_zzz_0_xxz_1,\
                               t_0_zzz_0_xxzz_0, t_0_zzz_0_xxzz_1, t_0_zzz_0_xyy_1,\
                               t_0_zzz_0_xyyy_0, t_0_zzz_0_xyyy_1, t_0_zzz_0_xyyz_0,\
                               t_0_zzz_0_xyyz_1, t_0_zzz_0_xyz_1, t_0_zzz_0_xyzz_0,\
                               t_0_zzz_0_xyzz_1, t_0_zzz_0_xzz_1, t_0_zzz_0_xzzz_0,\
                               t_0_zzz_0_xzzz_1, t_0_zzz_0_yyy_1, t_0_zzz_0_yyyy_0,\
                               t_0_zzz_0_yyyy_1, t_0_zzz_0_yyyz_0, t_0_zzz_0_yyyz_1,\
                               t_0_zzz_0_yyz_1, t_0_zzz_0_yyzz_0, t_0_zzz_0_yyzz_1,\
                               t_0_zzz_0_yzz_1, t_0_zzz_0_yzzz_0, t_0_zzz_0_yzzz_1,\
                               t_0_zzz_0_zzz_1, t_0_zzz_0_zzzz_0, t_0_zzz_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyyy_0_xxxz_0[i] += rpb_y[i] * t_0_yyy_0_xxxz_0[i] + rwp_y[i] * t_0_yyy_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxz_1[i];

            t_0_yyyy_0_xxxy_0[i] += rpb_y[i] * t_0_yyy_0_xxxy_0[i] + rwp_y[i] * t_0_yyy_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyy_0_xxxx_0[i] += rpb_y[i] * t_0_yyy_0_xxxx_0[i] + rwp_y[i] * t_0_yyy_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxx_1[i];

            t_0_xzzz_0_zzzz_0[i] += rpb_x[i] * t_0_zzz_0_zzzz_0[i] + rwp_x[i] * t_0_zzz_0_zzzz_1[i];

            t_0_xzzz_0_yzzz_0[i] += rpb_x[i] * t_0_zzz_0_yzzz_0[i] + rwp_x[i] * t_0_zzz_0_yzzz_1[i];

            t_0_xzzz_0_yyzz_0[i] += rpb_x[i] * t_0_zzz_0_yyzz_0[i] + rwp_x[i] * t_0_zzz_0_yyzz_1[i];

            t_0_xzzz_0_yyyz_0[i] += rpb_x[i] * t_0_zzz_0_yyyz_0[i] + rwp_x[i] * t_0_zzz_0_yyyz_1[i];

            t_0_xzzz_0_yyyy_0[i] += rpb_x[i] * t_0_zzz_0_yyyy_0[i] + rwp_x[i] * t_0_zzz_0_yyyy_1[i];

            t_0_xzzz_0_xzzz_0[i] += rpb_x[i] * t_0_zzz_0_xzzz_0[i] + rwp_x[i] * t_0_zzz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_xzzz_0_xyzz_0[i] += rpb_x[i] * t_0_zzz_0_xyzz_0[i] + rwp_x[i] * t_0_zzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_xzzz_0_xyyz_0[i] += rpb_x[i] * t_0_zzz_0_xyyz_0[i] + rwp_x[i] * t_0_zzz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_xzzz_0_xyyy_0[i] += rpb_x[i] * t_0_zzz_0_xyyy_0[i] + rwp_x[i] * t_0_zzz_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_xzzz_0_xxzz_0[i] += rpb_x[i] * t_0_zzz_0_xxzz_0[i] + rwp_x[i] * t_0_zzz_0_xxzz_1[i] + fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_xzzz_0_xxyz_0[i] += rpb_x[i] * t_0_zzz_0_xxyz_0[i] + rwp_x[i] * t_0_zzz_0_xxyz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_xzzz_0_xxyy_0[i] += rpb_x[i] * t_0_zzz_0_xxyy_0[i] + rwp_x[i] * t_0_zzz_0_xxyy_1[i] + fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_xzzz_0_xxxz_0[i] += rpb_x[i] * t_0_zzz_0_xxxz_0[i] + rwp_x[i] * t_0_zzz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_xzzz_0_xxxy_0[i] += rpb_x[i] * t_0_zzz_0_xxxy_0[i] + rwp_x[i] * t_0_zzz_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_xzzz_0_xxxx_0[i] += rpb_x[i] * t_0_zzz_0_xxxx_0[i] + rwp_x[i] * t_0_zzz_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_xyzz_0_zzzz_0[i] += rpb_x[i] * t_0_yzz_0_zzzz_0[i] + rwp_x[i] * t_0_yzz_0_zzzz_1[i];

            t_0_xyzz_0_yzzz_0[i] += rpb_x[i] * t_0_yzz_0_yzzz_0[i] + rwp_x[i] * t_0_yzz_0_yzzz_1[i];

            t_0_xyzz_0_yyzz_0[i] += rpb_x[i] * t_0_yzz_0_yyzz_0[i] + rwp_x[i] * t_0_yzz_0_yyzz_1[i];

            t_0_xyzz_0_yyyz_0[i] += rpb_x[i] * t_0_yzz_0_yyyz_0[i] + rwp_x[i] * t_0_yzz_0_yyyz_1[i];

            t_0_xyzz_0_yyyy_0[i] += rpb_x[i] * t_0_yzz_0_yyyy_0[i] + rwp_x[i] * t_0_yzz_0_yyyy_1[i];

            t_0_xyzz_0_xzzz_0[i] += rpb_y[i] * t_0_xzz_0_xzzz_0[i] + rwp_y[i] * t_0_xzz_0_xzzz_1[i];

            t_0_xyzz_0_xyzz_0[i] += rpb_x[i] * t_0_yzz_0_xyzz_0[i] + rwp_x[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyzz_0_xyyz_0[i] += rpb_x[i] * t_0_yzz_0_xyyz_0[i] + rwp_x[i] * t_0_yzz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yyz_1[i];

            t_0_xyzz_0_xyyy_0[i] += rpb_x[i] * t_0_yzz_0_xyyy_0[i] + rwp_x[i] * t_0_yzz_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yyy_1[i];

            t_0_xyzz_0_xxzz_0[i] += rpb_y[i] * t_0_xzz_0_xxzz_0[i] + rwp_y[i] * t_0_xzz_0_xxzz_1[i];

            t_0_xyzz_0_xxyz_0[i] += rpb_x[i] * t_0_yzz_0_xxyz_0[i] + rwp_x[i] * t_0_yzz_0_xxyz_1[i] + fze_0[i] * t_0_yzz_0_xyz_1[i];

            t_0_xyzz_0_xxyy_0[i] += rpb_x[i] * t_0_yzz_0_xxyy_0[i] + rwp_x[i] * t_0_yzz_0_xxyy_1[i] + fze_0[i] * t_0_yzz_0_xyy_1[i];

            t_0_xyzz_0_xxxz_0[i] += rpb_y[i] * t_0_xzz_0_xxxz_0[i] + rwp_y[i] * t_0_xzz_0_xxxz_1[i];

            t_0_xyzz_0_xxxy_0[i] += rpb_x[i] * t_0_yzz_0_xxxy_0[i] + rwp_x[i] * t_0_yzz_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_yzz_0_xxy_1[i];

            t_0_xyzz_0_xxxx_0[i] += rpb_y[i] * t_0_xzz_0_xxxx_0[i] + rwp_y[i] * t_0_xzz_0_xxxx_1[i];

            t_0_xyyz_0_zzzz_0[i] += rpb_x[i] * t_0_yyz_0_zzzz_0[i] + rwp_x[i] * t_0_yyz_0_zzzz_1[i];

            t_0_xyyz_0_yzzz_0[i] += rpb_x[i] * t_0_yyz_0_yzzz_0[i] + rwp_x[i] * t_0_yyz_0_yzzz_1[i];

            t_0_xyyz_0_yyzz_0[i] += rpb_x[i] * t_0_yyz_0_yyzz_0[i] + rwp_x[i] * t_0_yyz_0_yyzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_z, rwp_x, rwp_z, t_0_xx_0_xyyy_0,\
                               t_0_xx_0_xyyy_1, t_0_xxz_0_xyyy_0, t_0_xxz_0_xyyy_1,\
                               t_0_xxzz_0_xyyy_0, t_0_xxzz_0_xyyz_0, t_0_xxzz_0_xyzz_0,\
                               t_0_xxzz_0_xzzz_0, t_0_xxzz_0_yyyy_0, t_0_xxzz_0_yyyz_0,\
                               t_0_xxzz_0_yyzz_0, t_0_xxzz_0_yzzz_0, t_0_xxzz_0_zzzz_0,\
                               t_0_xyy_0_xxxx_0, t_0_xyy_0_xxxx_1, t_0_xyy_0_xxxy_0,\
                               t_0_xyy_0_xxxy_1, t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1,\
                               t_0_xyy_0_xyyy_0, t_0_xyy_0_xyyy_1, t_0_xyyy_0_xxxx_0,\
                               t_0_xyyy_0_xxxy_0, t_0_xyyy_0_xxxz_0, t_0_xyyy_0_xxyy_0,\
                               t_0_xyyy_0_xxyz_0, t_0_xyyy_0_xxzz_0, t_0_xyyy_0_xyyy_0,\
                               t_0_xyyy_0_xyyz_0, t_0_xyyy_0_xyzz_0, t_0_xyyy_0_xzzz_0,\
                               t_0_xyyy_0_yyyy_0, t_0_xyyy_0_yyyz_0, t_0_xyyy_0_yyzz_0,\
                               t_0_xyyy_0_yzzz_0, t_0_xyyy_0_zzzz_0, t_0_xyyz_0_xxxx_0,\
                               t_0_xyyz_0_xxxy_0, t_0_xyyz_0_xxxz_0, t_0_xyyz_0_xxyy_0,\
                               t_0_xyyz_0_xxyz_0, t_0_xyyz_0_xxzz_0, t_0_xyyz_0_xyyy_0,\
                               t_0_xyyz_0_xyyz_0, t_0_xyyz_0_xyzz_0, t_0_xyyz_0_xzzz_0,\
                               t_0_xyyz_0_yyyy_0, t_0_xyyz_0_yyyz_0, t_0_xzz_0_xyyz_0,\
                               t_0_xzz_0_xyyz_1, t_0_xzz_0_xyzz_0, t_0_xzz_0_xyzz_1,\
                               t_0_xzz_0_xzzz_0, t_0_xzz_0_xzzz_1, t_0_xzz_0_yyyy_0,\
                               t_0_xzz_0_yyyy_1, t_0_xzz_0_yyyz_0, t_0_xzz_0_yyyz_1,\
                               t_0_xzz_0_yyz_1, t_0_xzz_0_yyzz_0, t_0_xzz_0_yyzz_1,\
                               t_0_xzz_0_yzz_1, t_0_xzz_0_yzzz_0, t_0_xzz_0_yzzz_1,\
                               t_0_xzz_0_zzz_1, t_0_xzz_0_zzzz_0, t_0_xzz_0_zzzz_1,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyy_0_xxy_1, t_0_yyy_0_xxyy_0,\
                               t_0_yyy_0_xxyy_1, t_0_yyy_0_xxyz_0, t_0_yyy_0_xxyz_1,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xxzz_0, t_0_yyy_0_xxzz_1,\
                               t_0_yyy_0_xyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_xyyz_0, t_0_yyy_0_xyyz_1, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xyzz_0, t_0_yyy_0_xyzz_1, t_0_yyy_0_xzz_1,\
                               t_0_yyy_0_xzzz_0, t_0_yyy_0_xzzz_1, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1, t_0_yyy_0_yyyz_0,\
                               t_0_yyy_0_yyyz_1, t_0_yyy_0_yyz_1, t_0_yyy_0_yyzz_0,\
                               t_0_yyy_0_yyzz_1, t_0_yyy_0_yzz_1, t_0_yyy_0_yzzz_0,\
                               t_0_yyy_0_yzzz_1, t_0_yyy_0_zzz_1, t_0_yyy_0_zzzz_0,\
                               t_0_yyy_0_zzzz_1, t_0_yyz_0_xxxz_0, t_0_yyz_0_xxxz_1,\
                               t_0_yyz_0_xxyz_0, t_0_yyz_0_xxyz_1, t_0_yyz_0_xxz_1,\
                               t_0_yyz_0_xxzz_0, t_0_yyz_0_xxzz_1, t_0_yyz_0_xyyz_0,\
                               t_0_yyz_0_xyyz_1, t_0_yyz_0_xyz_1, t_0_yyz_0_xyzz_0,\
                               t_0_yyz_0_xyzz_1, t_0_yyz_0_xzz_1, t_0_yyz_0_xzzz_0,\
                               t_0_yyz_0_xzzz_1, t_0_yyz_0_yyyy_0, t_0_yyz_0_yyyy_1,\
                               t_0_yyz_0_yyyz_0, t_0_yyz_0_yyyz_1, t_0_yyz_0_yyz_1,\
                               t_0_yyz_0_yzz_1, t_0_yyz_0_zzz_1, t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1,\
                               t_0_zz_0_xyzz_0, t_0_zz_0_xyzz_1, t_0_zz_0_xzzz_0, t_0_zz_0_xzzz_1,\
                               t_0_zz_0_yyyy_0, t_0_zz_0_yyyy_1, t_0_zz_0_yyyz_0, t_0_zz_0_yyyz_1,\
                               t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1, t_0_zz_0_yzzz_0, t_0_zz_0_yzzz_1,\
                               t_0_zz_0_zzzz_0, t_0_zz_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyyz_0_yyyz_0[i] += rpb_x[i] * t_0_yyz_0_yyyz_0[i] + rwp_x[i] * t_0_yyz_0_yyyz_1[i];

            t_0_xyyz_0_yyyy_0[i] += rpb_x[i] * t_0_yyz_0_yyyy_0[i] + rwp_x[i] * t_0_yyz_0_yyyy_1[i];

            t_0_xyyz_0_xzzz_0[i] += rpb_x[i] * t_0_yyz_0_xzzz_0[i] + rwp_x[i] * t_0_yyz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_zzz_1[i];

            t_0_xyyz_0_xyzz_0[i] += rpb_x[i] * t_0_yyz_0_xyzz_0[i] + rwp_x[i] * t_0_yyz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yzz_1[i];

            t_0_xyyz_0_xyyz_0[i] += rpb_x[i] * t_0_yyz_0_xyyz_0[i] + rwp_x[i] * t_0_yyz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyz_0_xyyy_0[i] += rpb_z[i] * t_0_xyy_0_xyyy_0[i] + rwp_z[i] * t_0_xyy_0_xyyy_1[i];

            t_0_xyyz_0_xxzz_0[i] += rpb_x[i] * t_0_yyz_0_xxzz_0[i] + rwp_x[i] * t_0_yyz_0_xxzz_1[i] + fze_0[i] * t_0_yyz_0_xzz_1[i];

            t_0_xyyz_0_xxyz_0[i] += rpb_x[i] * t_0_yyz_0_xxyz_0[i] + rwp_x[i] * t_0_yyz_0_xxyz_1[i] + fze_0[i] * t_0_yyz_0_xyz_1[i];

            t_0_xyyz_0_xxyy_0[i] += rpb_z[i] * t_0_xyy_0_xxyy_0[i] + rwp_z[i] * t_0_xyy_0_xxyy_1[i];

            t_0_xyyz_0_xxxz_0[i] += rpb_x[i] * t_0_yyz_0_xxxz_0[i] + rwp_x[i] * t_0_yyz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_yyz_0_xxz_1[i];

            t_0_xyyz_0_xxxy_0[i] += rpb_z[i] * t_0_xyy_0_xxxy_0[i] + rwp_z[i] * t_0_xyy_0_xxxy_1[i];

            t_0_xyyz_0_xxxx_0[i] += rpb_z[i] * t_0_xyy_0_xxxx_0[i] + rwp_z[i] * t_0_xyy_0_xxxx_1[i];

            t_0_xyyy_0_zzzz_0[i] += rpb_x[i] * t_0_yyy_0_zzzz_0[i] + rwp_x[i] * t_0_yyy_0_zzzz_1[i];

            t_0_xyyy_0_yzzz_0[i] += rpb_x[i] * t_0_yyy_0_yzzz_0[i] + rwp_x[i] * t_0_yyy_0_yzzz_1[i];

            t_0_xyyy_0_yyzz_0[i] += rpb_x[i] * t_0_yyy_0_yyzz_0[i] + rwp_x[i] * t_0_yyy_0_yyzz_1[i];

            t_0_xyyy_0_yyyz_0[i] += rpb_x[i] * t_0_yyy_0_yyyz_0[i] + rwp_x[i] * t_0_yyy_0_yyyz_1[i];

            t_0_xyyy_0_yyyy_0[i] += rpb_x[i] * t_0_yyy_0_yyyy_0[i] + rwp_x[i] * t_0_yyy_0_yyyy_1[i];

            t_0_xyyy_0_xzzz_0[i] += rpb_x[i] * t_0_yyy_0_xzzz_0[i] + rwp_x[i] * t_0_yyy_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_xyyy_0_xyzz_0[i] += rpb_x[i] * t_0_yyy_0_xyzz_0[i] + rwp_x[i] * t_0_yyy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_xyyy_0_xyyz_0[i] += rpb_x[i] * t_0_yyy_0_xyyz_0[i] + rwp_x[i] * t_0_yyy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_xyyy_0_xyyy_0[i] += rpb_x[i] * t_0_yyy_0_xyyy_0[i] + rwp_x[i] * t_0_yyy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xyyy_0_xxzz_0[i] += rpb_x[i] * t_0_yyy_0_xxzz_0[i] + rwp_x[i] * t_0_yyy_0_xxzz_1[i] + fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_xyyy_0_xxyz_0[i] += rpb_x[i] * t_0_yyy_0_xxyz_0[i] + rwp_x[i] * t_0_yyy_0_xxyz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_xyyy_0_xxyy_0[i] += rpb_x[i] * t_0_yyy_0_xxyy_0[i] + rwp_x[i] * t_0_yyy_0_xxyy_1[i] + fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_xyyy_0_xxxz_0[i] += rpb_x[i] * t_0_yyy_0_xxxz_0[i] + rwp_x[i] * t_0_yyy_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_xyyy_0_xxxy_0[i] += rpb_x[i] * t_0_yyy_0_xxxy_0[i] + rwp_x[i] * t_0_yyy_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xxy_1[i];

            t_0_xyyy_0_xxxx_0[i] += rpb_x[i] * t_0_yyy_0_xxxx_0[i] + rwp_x[i] * t_0_yyy_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_xxzz_0_zzzz_0[i] += rpb_x[i] * t_0_xzz_0_zzzz_0[i] + rwp_x[i] * t_0_xzz_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i];

            t_0_xxzz_0_yzzz_0[i] += rpb_x[i] * t_0_xzz_0_yzzz_0[i] + rwp_x[i] * t_0_xzz_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i];

            t_0_xxzz_0_yyzz_0[i] += rpb_x[i] * t_0_xzz_0_yyzz_0[i] + rwp_x[i] * t_0_xzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i];

            t_0_xxzz_0_yyyz_0[i] += rpb_x[i] * t_0_xzz_0_yyyz_0[i] + rwp_x[i] * t_0_xzz_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i];

            t_0_xxzz_0_yyyy_0[i] += rpb_x[i] * t_0_xzz_0_yyyy_0[i] + rwp_x[i] * t_0_xzz_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyy_1[i];

            t_0_xxzz_0_xzzz_0[i] += rpb_x[i] * t_0_xzz_0_xzzz_0[i] + rwp_x[i] * t_0_xzz_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_zzz_1[i];

            t_0_xxzz_0_xyzz_0[i] += rpb_x[i] * t_0_xzz_0_xyzz_0[i] + rwp_x[i] * t_0_xzz_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yzz_1[i];

            t_0_xxzz_0_xyyz_0[i] += rpb_x[i] * t_0_xzz_0_xyyz_0[i] + rwp_x[i] * t_0_xzz_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yyz_1[i];

            t_0_xxzz_0_xyyy_0[i] += rpb_z[i] * t_0_xxz_0_xyyy_0[i] + rwp_z[i] * t_0_xxz_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xyyy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1, t_0_xx_0_xxxy_0,\
                               t_0_xx_0_xxxy_1, t_0_xx_0_xxxz_0, t_0_xx_0_xxxz_1, t_0_xx_0_xxyy_0,\
                               t_0_xx_0_xxyy_1, t_0_xx_0_xxzz_0, t_0_xx_0_xxzz_1, t_0_xx_0_xzzz_0,\
                               t_0_xx_0_xzzz_1, t_0_xxy_0_xxxx_0, t_0_xxy_0_xxxx_1,\
                               t_0_xxy_0_xxxy_0, t_0_xxy_0_xxxy_1, t_0_xxy_0_xxxz_0,\
                               t_0_xxy_0_xxxz_1, t_0_xxy_0_xxyy_0, t_0_xxy_0_xxyy_1,\
                               t_0_xxy_0_xxzz_0, t_0_xxy_0_xxzz_1, t_0_xxy_0_xyyy_0,\
                               t_0_xxy_0_xyyy_1, t_0_xxy_0_xzzz_0, t_0_xxy_0_xzzz_1,\
                               t_0_xxy_0_yyyy_0, t_0_xxy_0_yyyy_1, t_0_xxyy_0_xxxx_0,\
                               t_0_xxyy_0_xxxy_0, t_0_xxyy_0_xxxz_0, t_0_xxyy_0_xxyy_0,\
                               t_0_xxyy_0_xxyz_0, t_0_xxyy_0_xxzz_0, t_0_xxyy_0_xyyy_0,\
                               t_0_xxyy_0_xyyz_0, t_0_xxyy_0_xyzz_0, t_0_xxyy_0_xzzz_0,\
                               t_0_xxyy_0_yyyy_0, t_0_xxyy_0_yyyz_0, t_0_xxyy_0_yyzz_0,\
                               t_0_xxyy_0_yzzz_0, t_0_xxyy_0_zzzz_0, t_0_xxyz_0_xxxx_0,\
                               t_0_xxyz_0_xxxy_0, t_0_xxyz_0_xxxz_0, t_0_xxyz_0_xxyy_0,\
                               t_0_xxyz_0_xxyz_0, t_0_xxyz_0_xxzz_0, t_0_xxyz_0_xyyy_0,\
                               t_0_xxyz_0_xyyz_0, t_0_xxyz_0_xyzz_0, t_0_xxyz_0_xzzz_0,\
                               t_0_xxyz_0_yyyy_0, t_0_xxyz_0_yyyz_0, t_0_xxyz_0_yyzz_0,\
                               t_0_xxyz_0_yzzz_0, t_0_xxyz_0_zzzz_0, t_0_xxz_0_xxxx_0,\
                               t_0_xxz_0_xxxx_1, t_0_xxz_0_xxxy_0, t_0_xxz_0_xxxy_1,\
                               t_0_xxz_0_xxxz_0, t_0_xxz_0_xxxz_1, t_0_xxz_0_xxyy_0,\
                               t_0_xxz_0_xxyy_1, t_0_xxz_0_xxyz_0, t_0_xxz_0_xxyz_1,\
                               t_0_xxz_0_xxz_1, t_0_xxz_0_xxzz_0, t_0_xxz_0_xxzz_1,\
                               t_0_xxz_0_xyyz_0, t_0_xxz_0_xyyz_1, t_0_xxz_0_xyz_1,\
                               t_0_xxz_0_xyzz_0, t_0_xxz_0_xyzz_1, t_0_xxz_0_xzz_1,\
                               t_0_xxz_0_xzzz_0, t_0_xxz_0_xzzz_1, t_0_xxz_0_yyyz_0,\
                               t_0_xxz_0_yyyz_1, t_0_xxz_0_yyz_1, t_0_xxz_0_yyzz_0,\
                               t_0_xxz_0_yyzz_1, t_0_xxz_0_yzz_1, t_0_xxz_0_yzzz_0,\
                               t_0_xxz_0_yzzz_1, t_0_xxz_0_zzz_1, t_0_xxz_0_zzzz_0,\
                               t_0_xxz_0_zzzz_1, t_0_xxzz_0_xxxx_0, t_0_xxzz_0_xxxy_0,\
                               t_0_xxzz_0_xxxz_0, t_0_xxzz_0_xxyy_0, t_0_xxzz_0_xxyz_0,\
                               t_0_xxzz_0_xxzz_0, t_0_xyy_0_xxxy_0, t_0_xyy_0_xxxy_1,\
                               t_0_xyy_0_xxy_1, t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1,\
                               t_0_xyy_0_xxyz_0, t_0_xyy_0_xxyz_1, t_0_xyy_0_xyy_1,\
                               t_0_xyy_0_xyyy_0, t_0_xyy_0_xyyy_1, t_0_xyy_0_xyyz_0,\
                               t_0_xyy_0_xyyz_1, t_0_xyy_0_xyz_1, t_0_xyy_0_xyzz_0,\
                               t_0_xyy_0_xyzz_1, t_0_xyy_0_yyy_1, t_0_xyy_0_yyyy_0,\
                               t_0_xyy_0_yyyy_1, t_0_xyy_0_yyyz_0, t_0_xyy_0_yyyz_1,\
                               t_0_xyy_0_yyz_1, t_0_xyy_0_yyzz_0, t_0_xyy_0_yyzz_1,\
                               t_0_xyy_0_yzz_1, t_0_xyy_0_yzzz_0, t_0_xyy_0_yzzz_1,\
                               t_0_xyy_0_zzzz_0, t_0_xyy_0_zzzz_1, t_0_xzz_0_xxxz_0,\
                               t_0_xzz_0_xxxz_1, t_0_xzz_0_xxyz_0, t_0_xzz_0_xxyz_1,\
                               t_0_xzz_0_xxz_1, t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1,\
                               t_0_xzz_0_xyz_1, t_0_xzz_0_xzz_1, t_0_yy_0_xxxy_0, t_0_yy_0_xxxy_1,\
                               t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1, t_0_yy_0_xxyz_0, t_0_yy_0_xxyz_1,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyy_1, t_0_yy_0_xyyz_0, t_0_yy_0_xyyz_1,\
                               t_0_yy_0_xyzz_0, t_0_yy_0_xyzz_1, t_0_yy_0_yyyy_0, t_0_yy_0_yyyy_1,\
                               t_0_yy_0_yyyz_0, t_0_yy_0_yyyz_1, t_0_yy_0_yyzz_0, t_0_yy_0_yyzz_1,\
                               t_0_yy_0_yzzz_0, t_0_yy_0_yzzz_1, t_0_yy_0_zzzz_0, t_0_yy_0_zzzz_1,\
                               t_0_zz_0_xxxz_0, t_0_zz_0_xxxz_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1,\
                               t_0_zz_0_xxzz_0, t_0_zz_0_xxzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxzz_0_xxzz_0[i] += rpb_x[i] * t_0_xzz_0_xxzz_0[i] + rwp_x[i] * t_0_xzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_xzz_0_xzz_1[i];

            t_0_xxzz_0_xxyz_0[i] += rpb_x[i] * t_0_xzz_0_xxyz_0[i] + rwp_x[i] * t_0_xzz_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fze_0[i] * t_0_xzz_0_xyz_1[i];

            t_0_xxzz_0_xxyy_0[i] += rpb_z[i] * t_0_xxz_0_xxyy_0[i] + rwp_z[i] * t_0_xxz_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxyy_1[i];

            t_0_xxzz_0_xxxz_0[i] += rpb_x[i] * t_0_xzz_0_xxxz_0[i] + rwp_x[i] * t_0_xzz_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_xzz_0_xxz_1[i];

            t_0_xxzz_0_xxxy_0[i] += rpb_z[i] * t_0_xxz_0_xxxy_0[i] + rwp_z[i] * t_0_xxz_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxy_1[i];

            t_0_xxzz_0_xxxx_0[i] += rpb_z[i] * t_0_xxz_0_xxxx_0[i] + rwp_z[i] * t_0_xxz_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i];

            t_0_xxyz_0_zzzz_0[i] += rpb_y[i] * t_0_xxz_0_zzzz_0[i] + rwp_y[i] * t_0_xxz_0_zzzz_1[i];

            t_0_xxyz_0_yzzz_0[i] += rpb_y[i] * t_0_xxz_0_yzzz_0[i] + rwp_y[i] * t_0_xxz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_zzz_1[i];

            t_0_xxyz_0_yyzz_0[i] += rpb_y[i] * t_0_xxz_0_yyzz_0[i] + rwp_y[i] * t_0_xxz_0_yyzz_1[i] + fze_0[i] * t_0_xxz_0_yzz_1[i];

            t_0_xxyz_0_yyyz_0[i] += rpb_y[i] * t_0_xxz_0_yyyz_0[i] + rwp_y[i] * t_0_xxz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_xxz_0_yyz_1[i];

            t_0_xxyz_0_yyyy_0[i] += rpb_z[i] * t_0_xxy_0_yyyy_0[i] + rwp_z[i] * t_0_xxy_0_yyyy_1[i];

            t_0_xxyz_0_xzzz_0[i] += rpb_y[i] * t_0_xxz_0_xzzz_0[i] + rwp_y[i] * t_0_xxz_0_xzzz_1[i];

            t_0_xxyz_0_xyzz_0[i] += rpb_y[i] * t_0_xxz_0_xyzz_0[i] + rwp_y[i] * t_0_xxz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xzz_1[i];

            t_0_xxyz_0_xyyz_0[i] += rpb_y[i] * t_0_xxz_0_xyyz_0[i] + rwp_y[i] * t_0_xxz_0_xyyz_1[i] + fze_0[i] * t_0_xxz_0_xyz_1[i];

            t_0_xxyz_0_xyyy_0[i] += rpb_z[i] * t_0_xxy_0_xyyy_0[i] + rwp_z[i] * t_0_xxy_0_xyyy_1[i];

            t_0_xxyz_0_xxzz_0[i] += rpb_y[i] * t_0_xxz_0_xxzz_0[i] + rwp_y[i] * t_0_xxz_0_xxzz_1[i];

            t_0_xxyz_0_xxyz_0[i] += rpb_y[i] * t_0_xxz_0_xxyz_0[i] + rwp_y[i] * t_0_xxz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xxz_1[i];

            t_0_xxyz_0_xxyy_0[i] += rpb_z[i] * t_0_xxy_0_xxyy_0[i] + rwp_z[i] * t_0_xxy_0_xxyy_1[i];

            t_0_xxyz_0_xxxz_0[i] += rpb_y[i] * t_0_xxz_0_xxxz_0[i] + rwp_y[i] * t_0_xxz_0_xxxz_1[i];

            t_0_xxyz_0_xxxy_0[i] += rpb_z[i] * t_0_xxy_0_xxxy_0[i] + rwp_z[i] * t_0_xxy_0_xxxy_1[i];

            t_0_xxyz_0_xxxx_0[i] += rpb_y[i] * t_0_xxz_0_xxxx_0[i] + rwp_y[i] * t_0_xxz_0_xxxx_1[i];

            t_0_xxyy_0_zzzz_0[i] += rpb_x[i] * t_0_xyy_0_zzzz_0[i] + rwp_x[i] * t_0_xyy_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zzzz_1[i];

            t_0_xxyy_0_yzzz_0[i] += rpb_x[i] * t_0_xyy_0_yzzz_0[i] + rwp_x[i] * t_0_xyy_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yzzz_1[i];

            t_0_xxyy_0_yyzz_0[i] += rpb_x[i] * t_0_xyy_0_yyzz_0[i] + rwp_x[i] * t_0_xyy_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyzz_1[i];

            t_0_xxyy_0_yyyz_0[i] += rpb_x[i] * t_0_xyy_0_yyyz_0[i] + rwp_x[i] * t_0_xyy_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyz_1[i];

            t_0_xxyy_0_yyyy_0[i] += rpb_x[i] * t_0_xyy_0_yyyy_0[i] + rwp_x[i] * t_0_xyy_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i];

            t_0_xxyy_0_xzzz_0[i] += rpb_y[i] * t_0_xxy_0_xzzz_0[i] + rwp_y[i] * t_0_xxy_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xzzz_1[i];

            t_0_xxyy_0_xyzz_0[i] += rpb_x[i] * t_0_xyy_0_xyzz_0[i] + rwp_x[i] * t_0_xyy_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yzz_1[i];

            t_0_xxyy_0_xyyz_0[i] += rpb_x[i] * t_0_xyy_0_xyyz_0[i] + rwp_x[i] * t_0_xyy_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yyz_1[i];

            t_0_xxyy_0_xyyy_0[i] += rpb_x[i] * t_0_xyy_0_xyyy_0[i] + rwp_x[i] * t_0_xyy_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yyy_1[i];

            t_0_xxyy_0_xxzz_0[i] += rpb_y[i] * t_0_xxy_0_xxzz_0[i] + rwp_y[i] * t_0_xxy_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxzz_1[i];

            t_0_xxyy_0_xxyz_0[i] += rpb_x[i] * t_0_xyy_0_xxyz_0[i] + rwp_x[i] * t_0_xyy_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyz_1[i] + fze_0[i] * t_0_xyy_0_xyz_1[i];

            t_0_xxyy_0_xxyy_0[i] += rpb_x[i] * t_0_xyy_0_xxyy_0[i] + rwp_x[i] * t_0_xyy_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_xyy_0_xyy_1[i];

            t_0_xxyy_0_xxxz_0[i] += rpb_y[i] * t_0_xxy_0_xxxz_0[i] + rwp_y[i] * t_0_xxy_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxz_1[i];

            t_0_xxyy_0_xxxy_0[i] += rpb_x[i] * t_0_xyy_0_xxxy_0[i] + rwp_x[i] * t_0_xyy_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_xyy_0_xxy_1[i];

            t_0_xxyy_0_xxxx_0[i] += rpb_y[i] * t_0_xxy_0_xxxx_0[i] + rwp_y[i] * t_0_xxy_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xzzz_0, t_0_xx_0_xzzz_1, t_0_xx_0_yyyy_0,\
                               t_0_xx_0_yyyy_1, t_0_xx_0_yyyz_0, t_0_xx_0_yyyz_1, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yyzz_1, t_0_xx_0_yzzz_0, t_0_xx_0_yzzz_1, t_0_xx_0_zzzz_0,\
                               t_0_xx_0_zzzz_1, t_0_xxx_0_xxx_1, t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1,\
                               t_0_xxx_0_xxxy_0, t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0,\
                               t_0_xxx_0_xxxz_1, t_0_xxx_0_xxy_1, t_0_xxx_0_xxyy_0,\
                               t_0_xxx_0_xxyy_1, t_0_xxx_0_xxyz_0, t_0_xxx_0_xxyz_1,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xxzz_0, t_0_xxx_0_xxzz_1,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyyy_0, t_0_xxx_0_xyyy_1,\
                               t_0_xxx_0_xyyz_0, t_0_xxx_0_xyyz_1, t_0_xxx_0_xyz_1,\
                               t_0_xxx_0_xyzz_0, t_0_xxx_0_xyzz_1, t_0_xxx_0_xzz_1,\
                               t_0_xxx_0_xzzz_0, t_0_xxx_0_xzzz_1, t_0_xxx_0_yyy_1,\
                               t_0_xxx_0_yyyy_0, t_0_xxx_0_yyyy_1, t_0_xxx_0_yyyz_0,\
                               t_0_xxx_0_yyyz_1, t_0_xxx_0_yyz_1, t_0_xxx_0_yyzz_0,\
                               t_0_xxx_0_yyzz_1, t_0_xxx_0_yzz_1, t_0_xxx_0_yzzz_0,\
                               t_0_xxx_0_yzzz_1, t_0_xxx_0_zzz_1, t_0_xxx_0_zzzz_0,\
                               t_0_xxx_0_zzzz_1, t_0_xxxx_0_xzzz_0, t_0_xxxx_0_yyyy_0,\
                               t_0_xxxx_0_yyyz_0, t_0_xxxx_0_yyzz_0, t_0_xxxx_0_yzzz_0,\
                               t_0_xxxx_0_zzzz_0, t_0_xxxy_0_xxxx_0, t_0_xxxy_0_xxxy_0,\
                               t_0_xxxy_0_xxxz_0, t_0_xxxy_0_xxyy_0, t_0_xxxy_0_xxyz_0,\
                               t_0_xxxy_0_xxzz_0, t_0_xxxy_0_xyyy_0, t_0_xxxy_0_xyyz_0,\
                               t_0_xxxy_0_xyzz_0, t_0_xxxy_0_xzzz_0, t_0_xxxy_0_yyyy_0,\
                               t_0_xxxy_0_yyyz_0, t_0_xxxy_0_yyzz_0, t_0_xxxy_0_yzzz_0,\
                               t_0_xxxy_0_zzzz_0, t_0_xxxz_0_xxxx_0, t_0_xxxz_0_xxxy_0,\
                               t_0_xxxz_0_xxxz_0, t_0_xxxz_0_xxyy_0, t_0_xxxz_0_xxyz_0,\
                               t_0_xxxz_0_xxzz_0, t_0_xxxz_0_xyyy_0, t_0_xxxz_0_xyyz_0,\
                               t_0_xxxz_0_xyzz_0, t_0_xxxz_0_xzzz_0, t_0_xxxz_0_yyyy_0,\
                               t_0_xxxz_0_yyyz_0, t_0_xxxz_0_yyzz_0, t_0_xxxz_0_yzzz_0,\
                               t_0_xxxz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxz_0_zzzz_0[i] += rpb_z[i] * t_0_xxx_0_zzzz_0[i] + rwp_z[i] * t_0_xxx_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxz_0_yzzz_0[i] += rpb_z[i] * t_0_xxx_0_yzzz_0[i] + rwp_z[i] * t_0_xxx_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxz_0_yyzz_0[i] += rpb_z[i] * t_0_xxx_0_yyzz_0[i] + rwp_z[i] * t_0_xxx_0_yyzz_1[i] + fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxz_0_yyyz_0[i] += rpb_z[i] * t_0_xxx_0_yyyz_0[i] + rwp_z[i] * t_0_xxx_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxz_0_yyyy_0[i] += rpb_z[i] * t_0_xxx_0_yyyy_0[i] + rwp_z[i] * t_0_xxx_0_yyyy_1[i];

            t_0_xxxz_0_xzzz_0[i] += rpb_z[i] * t_0_xxx_0_xzzz_0[i] + rwp_z[i] * t_0_xxx_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxz_0_xyzz_0[i] += rpb_z[i] * t_0_xxx_0_xyzz_0[i] + rwp_z[i] * t_0_xxx_0_xyzz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxz_0_xyyz_0[i] += rpb_z[i] * t_0_xxx_0_xyyz_0[i] + rwp_z[i] * t_0_xxx_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxz_0_xyyy_0[i] += rpb_z[i] * t_0_xxx_0_xyyy_0[i] + rwp_z[i] * t_0_xxx_0_xyyy_1[i];

            t_0_xxxz_0_xxzz_0[i] += rpb_z[i] * t_0_xxx_0_xxzz_0[i] + rwp_z[i] * t_0_xxx_0_xxzz_1[i] + fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxz_0_xxyz_0[i] += rpb_z[i] * t_0_xxx_0_xxyz_0[i] + rwp_z[i] * t_0_xxx_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxz_0_xxyy_0[i] += rpb_z[i] * t_0_xxx_0_xxyy_0[i] + rwp_z[i] * t_0_xxx_0_xxyy_1[i];

            t_0_xxxz_0_xxxz_0[i] += rpb_z[i] * t_0_xxx_0_xxxz_0[i] + rwp_z[i] * t_0_xxx_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxz_0_xxxy_0[i] += rpb_z[i] * t_0_xxx_0_xxxy_0[i] + rwp_z[i] * t_0_xxx_0_xxxy_1[i];

            t_0_xxxz_0_xxxx_0[i] += rpb_z[i] * t_0_xxx_0_xxxx_0[i] + rwp_z[i] * t_0_xxx_0_xxxx_1[i];

            t_0_xxxy_0_zzzz_0[i] += rpb_y[i] * t_0_xxx_0_zzzz_0[i] + rwp_y[i] * t_0_xxx_0_zzzz_1[i];

            t_0_xxxy_0_yzzz_0[i] += rpb_y[i] * t_0_xxx_0_yzzz_0[i] + rwp_y[i] * t_0_xxx_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxy_0_yyzz_0[i] += rpb_y[i] * t_0_xxx_0_yyzz_0[i] + rwp_y[i] * t_0_xxx_0_yyzz_1[i] + fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxy_0_yyyz_0[i] += rpb_y[i] * t_0_xxx_0_yyyz_0[i] + rwp_y[i] * t_0_xxx_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxy_0_yyyy_0[i] += rpb_y[i] * t_0_xxx_0_yyyy_0[i] + rwp_y[i] * t_0_xxx_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxy_0_xzzz_0[i] += rpb_y[i] * t_0_xxx_0_xzzz_0[i] + rwp_y[i] * t_0_xxx_0_xzzz_1[i];

            t_0_xxxy_0_xyzz_0[i] += rpb_y[i] * t_0_xxx_0_xyzz_0[i] + rwp_y[i] * t_0_xxx_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxy_0_xyyz_0[i] += rpb_y[i] * t_0_xxx_0_xyyz_0[i] + rwp_y[i] * t_0_xxx_0_xyyz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxy_0_xyyy_0[i] += rpb_y[i] * t_0_xxx_0_xyyy_0[i] + rwp_y[i] * t_0_xxx_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxy_0_xxzz_0[i] += rpb_y[i] * t_0_xxx_0_xxzz_0[i] + rwp_y[i] * t_0_xxx_0_xxzz_1[i];

            t_0_xxxy_0_xxyz_0[i] += rpb_y[i] * t_0_xxx_0_xxyz_0[i] + rwp_y[i] * t_0_xxx_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxy_0_xxyy_0[i] += rpb_y[i] * t_0_xxx_0_xxyy_0[i] + rwp_y[i] * t_0_xxx_0_xxyy_1[i] + fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxy_0_xxxz_0[i] += rpb_y[i] * t_0_xxx_0_xxxz_0[i] + rwp_y[i] * t_0_xxx_0_xxxz_1[i];

            t_0_xxxy_0_xxxy_0[i] += rpb_y[i] * t_0_xxx_0_xxxy_0[i] + rwp_y[i] * t_0_xxx_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_xxxx_0[i] += rpb_y[i] * t_0_xxx_0_xxxx_0[i] + rwp_y[i] * t_0_xxx_0_xxxx_1[i];

            t_0_xxxx_0_zzzz_0[i] += rpb_x[i] * t_0_xxx_0_zzzz_0[i] + rwp_x[i] * t_0_xxx_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zzzz_1[i];

            t_0_xxxx_0_yzzz_0[i] += rpb_x[i] * t_0_xxx_0_yzzz_0[i] + rwp_x[i] * t_0_xxx_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yzzz_1[i];

            t_0_xxxx_0_yyzz_0[i] += rpb_x[i] * t_0_xxx_0_yyzz_0[i] + rwp_x[i] * t_0_xxx_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyzz_1[i];

            t_0_xxxx_0_yyyz_0[i] += rpb_x[i] * t_0_xxx_0_yyyz_0[i] + rwp_x[i] * t_0_xxx_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyyz_1[i];

            t_0_xxxx_0_yyyy_0[i] += rpb_x[i] * t_0_xxx_0_yyyy_0[i] + rwp_x[i] * t_0_xxx_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyyy_1[i];

            t_0_xxxx_0_xzzz_0[i] += rpb_x[i] * t_0_xxx_0_xzzz_0[i] + rwp_x[i] * t_0_xxx_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxy_1, t_0_xx_0_xxxz_0, t_0_xx_0_xxxz_1,\
                               t_0_xx_0_xxyy_0, t_0_xx_0_xxyy_1, t_0_xx_0_xxyz_0, t_0_xx_0_xxyz_1,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xxzz_1, t_0_xx_0_xyyy_0, t_0_xx_0_xyyy_1,\
                               t_0_xx_0_xyyz_0, t_0_xx_0_xyyz_1, t_0_xx_0_xyzz_0, t_0_xx_0_xyzz_1,\
                               t_0_xxx_0_xxx_1, t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1,\
                               t_0_xxx_0_xxxy_0, t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0,\
                               t_0_xxx_0_xxxz_1, t_0_xxx_0_xxy_1, t_0_xxx_0_xxyy_0,\
                               t_0_xxx_0_xxyy_1, t_0_xxx_0_xxyz_0, t_0_xxx_0_xxyz_1,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xxzz_0, t_0_xxx_0_xxzz_1,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyyy_0, t_0_xxx_0_xyyy_1,\
                               t_0_xxx_0_xyyz_0, t_0_xxx_0_xyyz_1, t_0_xxx_0_xyz_1,\
                               t_0_xxx_0_xyzz_0, t_0_xxx_0_xyzz_1, t_0_xxx_0_xzz_1,\
                               t_0_xxx_0_yyy_1, t_0_xxx_0_yyz_1, t_0_xxx_0_yzz_1, t_0_xxxx_0_xxxx_0,\
                               t_0_xxxx_0_xxxy_0, t_0_xxxx_0_xxxz_0, t_0_xxxx_0_xxyy_0,\
                               t_0_xxxx_0_xxyz_0, t_0_xxxx_0_xxzz_0, t_0_xxxx_0_xyyy_0,\
                               t_0_xxxx_0_xyyz_0, t_0_xxxx_0_xyzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxx_0_xyzz_0[i] += rpb_x[i] * t_0_xxx_0_xyzz_0[i] + rwp_x[i] * t_0_xxx_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxx_0_xyyz_0[i] += rpb_x[i] * t_0_xxx_0_xyyz_0[i] + rwp_x[i] * t_0_xxx_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxx_0_xyyy_0[i] += rpb_x[i] * t_0_xxx_0_xyyy_0[i] + rwp_x[i] * t_0_xxx_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxx_0_xxzz_0[i] += rpb_x[i] * t_0_xxx_0_xxzz_0[i] + rwp_x[i] * t_0_xxx_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxzz_1[i] + fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxx_0_xxyz_0[i] += rpb_x[i] * t_0_xxx_0_xxyz_0[i] + rwp_x[i] * t_0_xxx_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxyz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxx_0_xxyy_0[i] += rpb_x[i] * t_0_xxx_0_xxyy_0[i] + rwp_x[i] * t_0_xxx_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxyy_1[i] + fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxx_0_xxxz_0[i] += rpb_x[i] * t_0_xxx_0_xxxz_0[i] + rwp_x[i] * t_0_xxx_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxx_0_xxxy_0[i] += rpb_x[i] * t_0_xxx_0_xxxy_0[i] + rwp_x[i] * t_0_xxx_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxx_0_xxxx_0[i] += rpb_x[i] * t_0_xxx_0_xxxx_0[i] + rwp_x[i] * t_0_xxx_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_yyyy_0,\
                               t_0_yy_0_yyyy_1, t_0_yyz_0_yyyy_0, t_0_yyz_0_yyyy_1,\
                               t_0_yyzz_0_xzzz_0, t_0_yyzz_0_yyyy_0, t_0_yyzz_0_yyyz_0,\
                               t_0_yyzz_0_yyzz_0, t_0_yyzz_0_yzzz_0, t_0_yyzz_0_zzzz_0,\
                               t_0_yzz_0_xzzz_0, t_0_yzz_0_xzzz_1, t_0_yzz_0_yyyz_0,\
                               t_0_yzz_0_yyyz_1, t_0_yzz_0_yyz_1, t_0_yzz_0_yyzz_0,\
                               t_0_yzz_0_yyzz_1, t_0_yzz_0_yzz_1, t_0_yzz_0_yzzz_0,\
                               t_0_yzz_0_yzzz_1, t_0_yzz_0_zzz_1, t_0_yzz_0_zzzz_0,\
                               t_0_yzz_0_zzzz_1, t_0_yzzz_0_xxxx_0, t_0_yzzz_0_xxxy_0,\
                               t_0_yzzz_0_xxxz_0, t_0_yzzz_0_xxyy_0, t_0_yzzz_0_xxyz_0,\
                               t_0_yzzz_0_xxzz_0, t_0_yzzz_0_xyyy_0, t_0_yzzz_0_xyyz_0,\
                               t_0_yzzz_0_xyzz_0, t_0_yzzz_0_xzzz_0, t_0_yzzz_0_yyyy_0,\
                               t_0_yzzz_0_yyyz_0, t_0_yzzz_0_yyzz_0, t_0_yzzz_0_yzzz_0,\
                               t_0_yzzz_0_zzzz_0, t_0_zz_0_xxxx_0, t_0_zz_0_xxxx_1,\
                               t_0_zz_0_xxxy_0, t_0_zz_0_xxxy_1, t_0_zz_0_xxxz_0, t_0_zz_0_xxxz_1,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyy_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1,\
                               t_0_zz_0_xxzz_0, t_0_zz_0_xxzz_1, t_0_zz_0_xyyy_0, t_0_zz_0_xyyy_1,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1, t_0_zz_0_xyzz_0, t_0_zz_0_xyzz_1,\
                               t_0_zz_0_xzzz_0, t_0_zz_0_xzzz_1, t_0_zz_0_yyyy_0, t_0_zz_0_yyyy_1,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyyz_1, t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1,\
                               t_0_zz_0_yzzz_0, t_0_zz_0_yzzz_1, t_0_zz_0_zzzz_0, t_0_zz_0_zzzz_1,\
                               t_0_zzz_0_xxx_1, t_0_zzz_0_xxxx_0, t_0_zzz_0_xxxx_1,\
                               t_0_zzz_0_xxxy_0, t_0_zzz_0_xxxy_1, t_0_zzz_0_xxxz_0,\
                               t_0_zzz_0_xxxz_1, t_0_zzz_0_xxy_1, t_0_zzz_0_xxyy_0,\
                               t_0_zzz_0_xxyy_1, t_0_zzz_0_xxyz_0, t_0_zzz_0_xxyz_1,\
                               t_0_zzz_0_xxz_1, t_0_zzz_0_xxzz_0, t_0_zzz_0_xxzz_1,\
                               t_0_zzz_0_xyy_1, t_0_zzz_0_xyyy_0, t_0_zzz_0_xyyy_1,\
                               t_0_zzz_0_xyyz_0, t_0_zzz_0_xyyz_1, t_0_zzz_0_xyz_1,\
                               t_0_zzz_0_xyzz_0, t_0_zzz_0_xyzz_1, t_0_zzz_0_xzz_1,\
                               t_0_zzz_0_xzzz_0, t_0_zzz_0_xzzz_1, t_0_zzz_0_yyy_1,\
                               t_0_zzz_0_yyyy_0, t_0_zzz_0_yyyy_1, t_0_zzz_0_yyyz_0,\
                               t_0_zzz_0_yyyz_1, t_0_zzz_0_yyz_1, t_0_zzz_0_yyzz_0,\
                               t_0_zzz_0_yyzz_1, t_0_zzz_0_yzz_1, t_0_zzz_0_yzzz_0,\
                               t_0_zzz_0_yzzz_1, t_0_zzz_0_zzz_1, t_0_zzz_0_zzzz_0,\
                               t_0_zzz_0_zzzz_1, t_0_zzzz_0_xxxx_0, t_0_zzzz_0_xxxy_0,\
                               t_0_zzzz_0_xxxz_0, t_0_zzzz_0_xxyy_0, t_0_zzzz_0_xxyz_0,\
                               t_0_zzzz_0_xxzz_0, t_0_zzzz_0_xyyy_0, t_0_zzzz_0_xyyz_0,\
                               t_0_zzzz_0_xyzz_0, t_0_zzzz_0_xzzz_0, t_0_zzzz_0_yyyy_0,\
                               t_0_zzzz_0_yyyz_0, t_0_zzzz_0_yyzz_0, t_0_zzzz_0_yzzz_0,\
                               t_0_zzzz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzzz_0[i] = rpb_z[i] * t_0_zzz_0_zzzz_0[i] + rwp_z[i] * t_0_zzz_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_zzzz_0_yzzz_0[i] = rpb_z[i] * t_0_zzz_0_yzzz_0[i] + rwp_z[i] * t_0_zzz_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_zzzz_0_yyzz_0[i] = rpb_z[i] * t_0_zzz_0_yyzz_0[i] + rwp_z[i] * t_0_zzz_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_zzzz_0_yyyz_0[i] = rpb_z[i] * t_0_zzz_0_yyyz_0[i] + rwp_z[i] * t_0_zzz_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_zzzz_0_yyyy_0[i] = rpb_z[i] * t_0_zzz_0_yyyy_0[i] + rwp_z[i] * t_0_zzz_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyyy_1[i];

            t_0_zzzz_0_xzzz_0[i] = rpb_z[i] * t_0_zzz_0_xzzz_0[i] + rwp_z[i] * t_0_zzz_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_zzzz_0_xyzz_0[i] = rpb_z[i] * t_0_zzz_0_xyzz_0[i] + rwp_z[i] * t_0_zzz_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_zzzz_0_xyyz_0[i] = rpb_z[i] * t_0_zzz_0_xyyz_0[i] + rwp_z[i] * t_0_zzz_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_zzzz_0_xyyy_0[i] = rpb_z[i] * t_0_zzz_0_xyyy_0[i] + rwp_z[i] * t_0_zzz_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyyy_1[i];

            t_0_zzzz_0_xxzz_0[i] = rpb_z[i] * t_0_zzz_0_xxzz_0[i] + rwp_z[i] * t_0_zzz_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_zzzz_0_xxyz_0[i] = rpb_z[i] * t_0_zzz_0_xxyz_0[i] + rwp_z[i] * t_0_zzz_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_zzzz_0_xxyy_0[i] = rpb_z[i] * t_0_zzz_0_xxyy_0[i] + rwp_z[i] * t_0_zzz_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxyy_1[i];

            t_0_zzzz_0_xxxz_0[i] = rpb_z[i] * t_0_zzz_0_xxxz_0[i] + rwp_z[i] * t_0_zzz_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_zzzz_0_xxxy_0[i] = rpb_z[i] * t_0_zzz_0_xxxy_0[i] + rwp_z[i] * t_0_zzz_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxy_1[i];

            t_0_zzzz_0_xxxx_0[i] = rpb_z[i] * t_0_zzz_0_xxxx_0[i] + rwp_z[i] * t_0_zzz_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxxx_1[i];

            t_0_yzzz_0_zzzz_0[i] = rpb_y[i] * t_0_zzz_0_zzzz_0[i] + rwp_y[i] * t_0_zzz_0_zzzz_1[i];

            t_0_yzzz_0_yzzz_0[i] = rpb_y[i] * t_0_zzz_0_yzzz_0[i] + rwp_y[i] * t_0_zzz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yyzz_0[i] = rpb_y[i] * t_0_zzz_0_yyzz_0[i] + rwp_y[i] * t_0_zzz_0_yyzz_1[i] + fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_yzzz_0_yyyz_0[i] = rpb_y[i] * t_0_zzz_0_yyyz_0[i] + rwp_y[i] * t_0_zzz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_yzzz_0_yyyy_0[i] = rpb_y[i] * t_0_zzz_0_yyyy_0[i] + rwp_y[i] * t_0_zzz_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_yzzz_0_xzzz_0[i] = rpb_y[i] * t_0_zzz_0_xzzz_0[i] + rwp_y[i] * t_0_zzz_0_xzzz_1[i];

            t_0_yzzz_0_xyzz_0[i] = rpb_y[i] * t_0_zzz_0_xyzz_0[i] + rwp_y[i] * t_0_zzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_yzzz_0_xyyz_0[i] = rpb_y[i] * t_0_zzz_0_xyyz_0[i] + rwp_y[i] * t_0_zzz_0_xyyz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_yzzz_0_xyyy_0[i] = rpb_y[i] * t_0_zzz_0_xyyy_0[i] + rwp_y[i] * t_0_zzz_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_yzzz_0_xxzz_0[i] = rpb_y[i] * t_0_zzz_0_xxzz_0[i] + rwp_y[i] * t_0_zzz_0_xxzz_1[i];

            t_0_yzzz_0_xxyz_0[i] = rpb_y[i] * t_0_zzz_0_xxyz_0[i] + rwp_y[i] * t_0_zzz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_yzzz_0_xxyy_0[i] = rpb_y[i] * t_0_zzz_0_xxyy_0[i] + rwp_y[i] * t_0_zzz_0_xxyy_1[i] + fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_yzzz_0_xxxz_0[i] = rpb_y[i] * t_0_zzz_0_xxxz_0[i] + rwp_y[i] * t_0_zzz_0_xxxz_1[i];

            t_0_yzzz_0_xxxy_0[i] = rpb_y[i] * t_0_zzz_0_xxxy_0[i] + rwp_y[i] * t_0_zzz_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_yzzz_0_xxxx_0[i] = rpb_y[i] * t_0_zzz_0_xxxx_0[i] + rwp_y[i] * t_0_zzz_0_xxxx_1[i];

            t_0_yyzz_0_zzzz_0[i] = rpb_y[i] * t_0_yzz_0_zzzz_0[i] + rwp_y[i] * t_0_yzz_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i];

            t_0_yyzz_0_yzzz_0[i] = rpb_y[i] * t_0_yzz_0_yzzz_0[i] + rwp_y[i] * t_0_yzz_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_zzz_1[i];

            t_0_yyzz_0_yyzz_0[i] = rpb_y[i] * t_0_yzz_0_yyzz_0[i] + rwp_y[i] * t_0_yzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i] + fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_yyzz_0_yyyz_0[i] = rpb_y[i] * t_0_yzz_0_yyyz_0[i] + rwp_y[i] * t_0_yzz_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_yzz_0_yyz_1[i];

            t_0_yyzz_0_yyyy_0[i] = rpb_z[i] * t_0_yyz_0_yyyy_0[i] + rwp_z[i] * t_0_yyz_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i];

            t_0_yyzz_0_xzzz_0[i] = rpb_y[i] * t_0_yzz_0_xzzz_0[i] + rwp_y[i] * t_0_yzz_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_xxxy_0,\
                               t_0_yy_0_xxxy_1, t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1, t_0_yy_0_xxyz_0,\
                               t_0_yy_0_xxyz_1, t_0_yy_0_xxzz_0, t_0_yy_0_xxzz_1, t_0_yy_0_xyyy_0,\
                               t_0_yy_0_xyyy_1, t_0_yy_0_xyyz_0, t_0_yy_0_xyyz_1, t_0_yy_0_xyzz_0,\
                               t_0_yy_0_xyzz_1, t_0_yy_0_xzzz_0, t_0_yy_0_xzzz_1, t_0_yy_0_yyyy_0,\
                               t_0_yy_0_yyyy_1, t_0_yy_0_yyyz_0, t_0_yy_0_yyyz_1, t_0_yy_0_yyzz_0,\
                               t_0_yy_0_yyzz_1, t_0_yy_0_yzzz_0, t_0_yy_0_yzzz_1, t_0_yy_0_zzzz_0,\
                               t_0_yy_0_zzzz_1, t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyy_0_xxy_1, t_0_yyy_0_xxyy_0,\
                               t_0_yyy_0_xxyy_1, t_0_yyy_0_xxyz_0, t_0_yyy_0_xxyz_1,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xxzz_0, t_0_yyy_0_xxzz_1,\
                               t_0_yyy_0_xyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_xyyz_0, t_0_yyy_0_xyyz_1, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xyzz_0, t_0_yyy_0_xyzz_1, t_0_yyy_0_xzz_1,\
                               t_0_yyy_0_xzzz_0, t_0_yyy_0_xzzz_1, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1, t_0_yyy_0_yyyz_0,\
                               t_0_yyy_0_yyyz_1, t_0_yyy_0_yyz_1, t_0_yyy_0_yyzz_0,\
                               t_0_yyy_0_yyzz_1, t_0_yyy_0_yzz_1, t_0_yyy_0_yzzz_0,\
                               t_0_yyy_0_yzzz_1, t_0_yyy_0_zzz_1, t_0_yyy_0_zzzz_0,\
                               t_0_yyy_0_zzzz_1, t_0_yyyy_0_xxyy_0, t_0_yyyy_0_xxyz_0,\
                               t_0_yyyy_0_xxzz_0, t_0_yyyy_0_xyyy_0, t_0_yyyy_0_xyyz_0,\
                               t_0_yyyy_0_xyzz_0, t_0_yyyy_0_xzzz_0, t_0_yyyy_0_yyyy_0,\
                               t_0_yyyy_0_yyyz_0, t_0_yyyy_0_yyzz_0, t_0_yyyy_0_yzzz_0,\
                               t_0_yyyy_0_zzzz_0, t_0_yyyz_0_xxxx_0, t_0_yyyz_0_xxxy_0,\
                               t_0_yyyz_0_xxxz_0, t_0_yyyz_0_xxyy_0, t_0_yyyz_0_xxyz_0,\
                               t_0_yyyz_0_xxzz_0, t_0_yyyz_0_xyyy_0, t_0_yyyz_0_xyyz_0,\
                               t_0_yyyz_0_xyzz_0, t_0_yyyz_0_xzzz_0, t_0_yyyz_0_yyyy_0,\
                               t_0_yyyz_0_yyyz_0, t_0_yyyz_0_yyzz_0, t_0_yyyz_0_yzzz_0,\
                               t_0_yyyz_0_zzzz_0, t_0_yyz_0_xxxy_0, t_0_yyz_0_xxxy_1,\
                               t_0_yyz_0_xxyy_0, t_0_yyz_0_xxyy_1, t_0_yyz_0_xyyy_0,\
                               t_0_yyz_0_xyyy_1, t_0_yyzz_0_xxxx_0, t_0_yyzz_0_xxxy_0,\
                               t_0_yyzz_0_xxxz_0, t_0_yyzz_0_xxyy_0, t_0_yyzz_0_xxyz_0,\
                               t_0_yyzz_0_xxzz_0, t_0_yyzz_0_xyyy_0, t_0_yyzz_0_xyyz_0,\
                               t_0_yyzz_0_xyzz_0, t_0_yzz_0_xxxx_0, t_0_yzz_0_xxxx_1,\
                               t_0_yzz_0_xxxz_0, t_0_yzz_0_xxxz_1, t_0_yzz_0_xxyz_0,\
                               t_0_yzz_0_xxyz_1, t_0_yzz_0_xxz_1, t_0_yzz_0_xxzz_0,\
                               t_0_yzz_0_xxzz_1, t_0_yzz_0_xyyz_0, t_0_yzz_0_xyyz_1,\
                               t_0_yzz_0_xyz_1, t_0_yzz_0_xyzz_0, t_0_yzz_0_xyzz_1,\
                               t_0_yzz_0_xzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxx_1, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxxz_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1, t_0_zz_0_xxzz_0,\
                               t_0_zz_0_xxzz_1, t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1, t_0_zz_0_xyzz_0,\
                               t_0_zz_0_xyzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyzz_0_xyzz_0[i] = rpb_y[i] * t_0_yzz_0_xyzz_0[i] + rwp_y[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xzz_1[i];

            t_0_yyzz_0_xyyz_0[i] = rpb_y[i] * t_0_yzz_0_xyyz_0[i] + rwp_y[i] * t_0_yzz_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fze_0[i] * t_0_yzz_0_xyz_1[i];

            t_0_yyzz_0_xyyy_0[i] = rpb_z[i] * t_0_yyz_0_xyyy_0[i] + rwp_z[i] * t_0_yyz_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i];

            t_0_yyzz_0_xxzz_0[i] = rpb_y[i] * t_0_yzz_0_xxzz_0[i] + rwp_y[i] * t_0_yzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i];

            t_0_yyzz_0_xxyz_0[i] = rpb_y[i] * t_0_yzz_0_xxyz_0[i] + rwp_y[i] * t_0_yzz_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xxz_1[i];

            t_0_yyzz_0_xxyy_0[i] = rpb_z[i] * t_0_yyz_0_xxyy_0[i] + rwp_z[i] * t_0_yyz_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i];

            t_0_yyzz_0_xxxz_0[i] = rpb_y[i] * t_0_yzz_0_xxxz_0[i] + rwp_y[i] * t_0_yzz_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i];

            t_0_yyzz_0_xxxy_0[i] = rpb_z[i] * t_0_yyz_0_xxxy_0[i] + rwp_z[i] * t_0_yyz_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i];

            t_0_yyzz_0_xxxx_0[i] = rpb_y[i] * t_0_yzz_0_xxxx_0[i] + rwp_y[i] * t_0_yzz_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxx_1[i];

            t_0_yyyz_0_zzzz_0[i] = rpb_z[i] * t_0_yyy_0_zzzz_0[i] + rwp_z[i] * t_0_yyy_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_yyyz_0_yzzz_0[i] = rpb_z[i] * t_0_yyy_0_yzzz_0[i] + rwp_z[i] * t_0_yyy_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_yyyz_0_yyzz_0[i] = rpb_z[i] * t_0_yyy_0_yyzz_0[i] + rwp_z[i] * t_0_yyy_0_yyzz_1[i] + fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_yyyz_0_yyyz_0[i] = rpb_z[i] * t_0_yyy_0_yyyz_0[i] + rwp_z[i] * t_0_yyy_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyz_0_yyyy_0[i] = rpb_z[i] * t_0_yyy_0_yyyy_0[i] + rwp_z[i] * t_0_yyy_0_yyyy_1[i];

            t_0_yyyz_0_xzzz_0[i] = rpb_z[i] * t_0_yyy_0_xzzz_0[i] + rwp_z[i] * t_0_yyy_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_yyyz_0_xyzz_0[i] = rpb_z[i] * t_0_yyy_0_xyzz_0[i] + rwp_z[i] * t_0_yyy_0_xyzz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_yyyz_0_xyyz_0[i] = rpb_z[i] * t_0_yyy_0_xyyz_0[i] + rwp_z[i] * t_0_yyy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyz_0_xyyy_0[i] = rpb_z[i] * t_0_yyy_0_xyyy_0[i] + rwp_z[i] * t_0_yyy_0_xyyy_1[i];

            t_0_yyyz_0_xxzz_0[i] = rpb_z[i] * t_0_yyy_0_xxzz_0[i] + rwp_z[i] * t_0_yyy_0_xxzz_1[i] + fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_yyyz_0_xxyz_0[i] = rpb_z[i] * t_0_yyy_0_xxyz_0[i] + rwp_z[i] * t_0_yyy_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxy_1[i];

            t_0_yyyz_0_xxyy_0[i] = rpb_z[i] * t_0_yyy_0_xxyy_0[i] + rwp_z[i] * t_0_yyy_0_xxyy_1[i];

            t_0_yyyz_0_xxxz_0[i] = rpb_z[i] * t_0_yyy_0_xxxz_0[i] + rwp_z[i] * t_0_yyy_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyz_0_xxxy_0[i] = rpb_z[i] * t_0_yyy_0_xxxy_0[i] + rwp_z[i] * t_0_yyy_0_xxxy_1[i];

            t_0_yyyz_0_xxxx_0[i] = rpb_z[i] * t_0_yyy_0_xxxx_0[i] + rwp_z[i] * t_0_yyy_0_xxxx_1[i];

            t_0_yyyy_0_zzzz_0[i] = rpb_y[i] * t_0_yyy_0_zzzz_0[i] + rwp_y[i] * t_0_yyy_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zzzz_1[i];

            t_0_yyyy_0_yzzz_0[i] = rpb_y[i] * t_0_yyy_0_yzzz_0[i] + rwp_y[i] * t_0_yyy_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_yyyy_0_yyzz_0[i] = rpb_y[i] * t_0_yyy_0_yyzz_0[i] + rwp_y[i] * t_0_yyy_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyzz_1[i] + fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_yyyy_0_yyyz_0[i] = rpb_y[i] * t_0_yyy_0_yyyz_0[i] + rwp_y[i] * t_0_yyy_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_yyyy_0_yyyy_0[i] = rpb_y[i] * t_0_yyy_0_yyyy_0[i] + rwp_y[i] * t_0_yyy_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyy_0_xzzz_0[i] = rpb_y[i] * t_0_yyy_0_xzzz_0[i] + rwp_y[i] * t_0_yyy_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xzzz_1[i];

            t_0_yyyy_0_xyzz_0[i] = rpb_y[i] * t_0_yyy_0_xyzz_0[i] + rwp_y[i] * t_0_yyy_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_yyyy_0_xyyz_0[i] = rpb_y[i] * t_0_yyy_0_xyyz_0[i] + rwp_y[i] * t_0_yyy_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyyz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_yyyy_0_xyyy_0[i] = rpb_y[i] * t_0_yyy_0_xyyy_0[i] + rwp_y[i] * t_0_yyy_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyy_0_xxzz_0[i] = rpb_y[i] * t_0_yyy_0_xxzz_0[i] + rwp_y[i] * t_0_yyy_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxzz_1[i];

            t_0_yyyy_0_xxyz_0[i] = rpb_y[i] * t_0_yyy_0_xxyz_0[i] + rwp_y[i] * t_0_yyy_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_yyyy_0_xxyy_0[i] = rpb_y[i] * t_0_yyy_0_xxyy_0[i] + rwp_y[i] * t_0_yyy_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_yyy_0_xxy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rwp_x, rwp_y, t_0_xyyz_0_yyzz_0,\
                               t_0_xyyz_0_yzzz_0, t_0_xyyz_0_zzzz_0, t_0_xyzz_0_xxxx_0,\
                               t_0_xyzz_0_xxxy_0, t_0_xyzz_0_xxxz_0, t_0_xyzz_0_xxyy_0,\
                               t_0_xyzz_0_xxyz_0, t_0_xyzz_0_xxzz_0, t_0_xyzz_0_xyyy_0,\
                               t_0_xyzz_0_xyyz_0, t_0_xyzz_0_xyzz_0, t_0_xyzz_0_xzzz_0,\
                               t_0_xyzz_0_yyyy_0, t_0_xyzz_0_yyyz_0, t_0_xyzz_0_yyzz_0,\
                               t_0_xyzz_0_yzzz_0, t_0_xyzz_0_zzzz_0, t_0_xzz_0_xxxx_0,\
                               t_0_xzz_0_xxxx_1, t_0_xzz_0_xxxz_0, t_0_xzz_0_xxxz_1,\
                               t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1, t_0_xzz_0_xzzz_0,\
                               t_0_xzz_0_xzzz_1, t_0_xzzz_0_xxxx_0, t_0_xzzz_0_xxxy_0,\
                               t_0_xzzz_0_xxxz_0, t_0_xzzz_0_xxyy_0, t_0_xzzz_0_xxyz_0,\
                               t_0_xzzz_0_xxzz_0, t_0_xzzz_0_xyyy_0, t_0_xzzz_0_xyyz_0,\
                               t_0_xzzz_0_xyzz_0, t_0_xzzz_0_xzzz_0, t_0_xzzz_0_yyyy_0,\
                               t_0_xzzz_0_yyyz_0, t_0_xzzz_0_yyzz_0, t_0_xzzz_0_yzzz_0,\
                               t_0_xzzz_0_zzzz_0, t_0_yy_0_xxxx_0, t_0_yy_0_xxxx_1,\
                               t_0_yy_0_xxxy_0, t_0_yy_0_xxxy_1, t_0_yy_0_xxxz_0, t_0_yy_0_xxxz_1,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyyy_0_xxxx_0, t_0_yyyy_0_xxxy_0,\
                               t_0_yyyy_0_xxxz_0, t_0_yyz_0_yyzz_0, t_0_yyz_0_yyzz_1,\
                               t_0_yyz_0_yzzz_0, t_0_yyz_0_yzzz_1, t_0_yyz_0_zzzz_0,\
                               t_0_yyz_0_zzzz_1, t_0_yzz_0_xxxy_0, t_0_yzz_0_xxxy_1,\
                               t_0_yzz_0_xxy_1, t_0_yzz_0_xxyy_0, t_0_yzz_0_xxyy_1,\
                               t_0_yzz_0_xxyz_0, t_0_yzz_0_xxyz_1, t_0_yzz_0_xyy_1,\
                               t_0_yzz_0_xyyy_0, t_0_yzz_0_xyyy_1, t_0_yzz_0_xyyz_0,\
                               t_0_yzz_0_xyyz_1, t_0_yzz_0_xyz_1, t_0_yzz_0_xyzz_0,\
                               t_0_yzz_0_xyzz_1, t_0_yzz_0_yyy_1, t_0_yzz_0_yyyy_0,\
                               t_0_yzz_0_yyyy_1, t_0_yzz_0_yyyz_0, t_0_yzz_0_yyyz_1,\
                               t_0_yzz_0_yyz_1, t_0_yzz_0_yyzz_0, t_0_yzz_0_yyzz_1,\
                               t_0_yzz_0_yzz_1, t_0_yzz_0_yzzz_0, t_0_yzz_0_yzzz_1,\
                               t_0_yzz_0_zzzz_0, t_0_yzz_0_zzzz_1, t_0_zzz_0_xxx_1,\
                               t_0_zzz_0_xxxx_0, t_0_zzz_0_xxxx_1, t_0_zzz_0_xxxy_0,\
                               t_0_zzz_0_xxxy_1, t_0_zzz_0_xxxz_0, t_0_zzz_0_xxxz_1,\
                               t_0_zzz_0_xxy_1, t_0_zzz_0_xxyy_0, t_0_zzz_0_xxyy_1,\
                               t_0_zzz_0_xxyz_0, t_0_zzz_0_xxyz_1, t_0_zzz_0_xxz_1,\
                               t_0_zzz_0_xxzz_0, t_0_zzz_0_xxzz_1, t_0_zzz_0_xyy_1,\
                               t_0_zzz_0_xyyy_0, t_0_zzz_0_xyyy_1, t_0_zzz_0_xyyz_0,\
                               t_0_zzz_0_xyyz_1, t_0_zzz_0_xyz_1, t_0_zzz_0_xyzz_0,\
                               t_0_zzz_0_xyzz_1, t_0_zzz_0_xzz_1, t_0_zzz_0_xzzz_0,\
                               t_0_zzz_0_xzzz_1, t_0_zzz_0_yyy_1, t_0_zzz_0_yyyy_0,\
                               t_0_zzz_0_yyyy_1, t_0_zzz_0_yyyz_0, t_0_zzz_0_yyyz_1,\
                               t_0_zzz_0_yyz_1, t_0_zzz_0_yyzz_0, t_0_zzz_0_yyzz_1,\
                               t_0_zzz_0_yzz_1, t_0_zzz_0_yzzz_0, t_0_zzz_0_yzzz_1,\
                               t_0_zzz_0_zzz_1, t_0_zzz_0_zzzz_0, t_0_zzz_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyyy_0_xxxz_0[i] = rpb_y[i] * t_0_yyy_0_xxxz_0[i] + rwp_y[i] * t_0_yyy_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxz_1[i];

            t_0_yyyy_0_xxxy_0[i] = rpb_y[i] * t_0_yyy_0_xxxy_0[i] + rwp_y[i] * t_0_yyy_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyy_0_xxxx_0[i] = rpb_y[i] * t_0_yyy_0_xxxx_0[i] + rwp_y[i] * t_0_yyy_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxxx_1[i];

            t_0_xzzz_0_zzzz_0[i] = rpb_x[i] * t_0_zzz_0_zzzz_0[i] + rwp_x[i] * t_0_zzz_0_zzzz_1[i];

            t_0_xzzz_0_yzzz_0[i] = rpb_x[i] * t_0_zzz_0_yzzz_0[i] + rwp_x[i] * t_0_zzz_0_yzzz_1[i];

            t_0_xzzz_0_yyzz_0[i] = rpb_x[i] * t_0_zzz_0_yyzz_0[i] + rwp_x[i] * t_0_zzz_0_yyzz_1[i];

            t_0_xzzz_0_yyyz_0[i] = rpb_x[i] * t_0_zzz_0_yyyz_0[i] + rwp_x[i] * t_0_zzz_0_yyyz_1[i];

            t_0_xzzz_0_yyyy_0[i] = rpb_x[i] * t_0_zzz_0_yyyy_0[i] + rwp_x[i] * t_0_zzz_0_yyyy_1[i];

            t_0_xzzz_0_xzzz_0[i] = rpb_x[i] * t_0_zzz_0_xzzz_0[i] + rwp_x[i] * t_0_zzz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zzz_1[i];

            t_0_xzzz_0_xyzz_0[i] = rpb_x[i] * t_0_zzz_0_xyzz_0[i] + rwp_x[i] * t_0_zzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yzz_1[i];

            t_0_xzzz_0_xyyz_0[i] = rpb_x[i] * t_0_zzz_0_xyyz_0[i] + rwp_x[i] * t_0_zzz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyz_1[i];

            t_0_xzzz_0_xyyy_0[i] = rpb_x[i] * t_0_zzz_0_xyyy_0[i] + rwp_x[i] * t_0_zzz_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yyy_1[i];

            t_0_xzzz_0_xxzz_0[i] = rpb_x[i] * t_0_zzz_0_xxzz_0[i] + rwp_x[i] * t_0_zzz_0_xxzz_1[i] + fze_0[i] * t_0_zzz_0_xzz_1[i];

            t_0_xzzz_0_xxyz_0[i] = rpb_x[i] * t_0_zzz_0_xxyz_0[i] + rwp_x[i] * t_0_zzz_0_xxyz_1[i] + fze_0[i] * t_0_zzz_0_xyz_1[i];

            t_0_xzzz_0_xxyy_0[i] = rpb_x[i] * t_0_zzz_0_xxyy_0[i] + rwp_x[i] * t_0_zzz_0_xxyy_1[i] + fze_0[i] * t_0_zzz_0_xyy_1[i];

            t_0_xzzz_0_xxxz_0[i] = rpb_x[i] * t_0_zzz_0_xxxz_0[i] + rwp_x[i] * t_0_zzz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xxz_1[i];

            t_0_xzzz_0_xxxy_0[i] = rpb_x[i] * t_0_zzz_0_xxxy_0[i] + rwp_x[i] * t_0_zzz_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xxy_1[i];

            t_0_xzzz_0_xxxx_0[i] = rpb_x[i] * t_0_zzz_0_xxxx_0[i] + rwp_x[i] * t_0_zzz_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_zzz_0_xxx_1[i];

            t_0_xyzz_0_zzzz_0[i] = rpb_x[i] * t_0_yzz_0_zzzz_0[i] + rwp_x[i] * t_0_yzz_0_zzzz_1[i];

            t_0_xyzz_0_yzzz_0[i] = rpb_x[i] * t_0_yzz_0_yzzz_0[i] + rwp_x[i] * t_0_yzz_0_yzzz_1[i];

            t_0_xyzz_0_yyzz_0[i] = rpb_x[i] * t_0_yzz_0_yyzz_0[i] + rwp_x[i] * t_0_yzz_0_yyzz_1[i];

            t_0_xyzz_0_yyyz_0[i] = rpb_x[i] * t_0_yzz_0_yyyz_0[i] + rwp_x[i] * t_0_yzz_0_yyyz_1[i];

            t_0_xyzz_0_yyyy_0[i] = rpb_x[i] * t_0_yzz_0_yyyy_0[i] + rwp_x[i] * t_0_yzz_0_yyyy_1[i];

            t_0_xyzz_0_xzzz_0[i] = rpb_y[i] * t_0_xzz_0_xzzz_0[i] + rwp_y[i] * t_0_xzz_0_xzzz_1[i];

            t_0_xyzz_0_xyzz_0[i] = rpb_x[i] * t_0_yzz_0_xyzz_0[i] + rwp_x[i] * t_0_yzz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyzz_0_xyyz_0[i] = rpb_x[i] * t_0_yzz_0_xyyz_0[i] + rwp_x[i] * t_0_yzz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yyz_1[i];

            t_0_xyzz_0_xyyy_0[i] = rpb_x[i] * t_0_yzz_0_xyyy_0[i] + rwp_x[i] * t_0_yzz_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yyy_1[i];

            t_0_xyzz_0_xxzz_0[i] = rpb_y[i] * t_0_xzz_0_xxzz_0[i] + rwp_y[i] * t_0_xzz_0_xxzz_1[i];

            t_0_xyzz_0_xxyz_0[i] = rpb_x[i] * t_0_yzz_0_xxyz_0[i] + rwp_x[i] * t_0_yzz_0_xxyz_1[i] + fze_0[i] * t_0_yzz_0_xyz_1[i];

            t_0_xyzz_0_xxyy_0[i] = rpb_x[i] * t_0_yzz_0_xxyy_0[i] + rwp_x[i] * t_0_yzz_0_xxyy_1[i] + fze_0[i] * t_0_yzz_0_xyy_1[i];

            t_0_xyzz_0_xxxz_0[i] = rpb_y[i] * t_0_xzz_0_xxxz_0[i] + rwp_y[i] * t_0_xzz_0_xxxz_1[i];

            t_0_xyzz_0_xxxy_0[i] = rpb_x[i] * t_0_yzz_0_xxxy_0[i] + rwp_x[i] * t_0_yzz_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_yzz_0_xxy_1[i];

            t_0_xyzz_0_xxxx_0[i] = rpb_y[i] * t_0_xzz_0_xxxx_0[i] + rwp_y[i] * t_0_xzz_0_xxxx_1[i];

            t_0_xyyz_0_zzzz_0[i] = rpb_x[i] * t_0_yyz_0_zzzz_0[i] + rwp_x[i] * t_0_yyz_0_zzzz_1[i];

            t_0_xyyz_0_yzzz_0[i] = rpb_x[i] * t_0_yyz_0_yzzz_0[i] + rwp_x[i] * t_0_yyz_0_yzzz_1[i];

            t_0_xyyz_0_yyzz_0[i] = rpb_x[i] * t_0_yyz_0_yyzz_0[i] + rwp_x[i] * t_0_yyz_0_yyzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_z, rwp_x, rwp_z, t_0_xx_0_xyyy_0,\
                               t_0_xx_0_xyyy_1, t_0_xxz_0_xyyy_0, t_0_xxz_0_xyyy_1,\
                               t_0_xxzz_0_xyyy_0, t_0_xxzz_0_xyyz_0, t_0_xxzz_0_xyzz_0,\
                               t_0_xxzz_0_xzzz_0, t_0_xxzz_0_yyyy_0, t_0_xxzz_0_yyyz_0,\
                               t_0_xxzz_0_yyzz_0, t_0_xxzz_0_yzzz_0, t_0_xxzz_0_zzzz_0,\
                               t_0_xyy_0_xxxx_0, t_0_xyy_0_xxxx_1, t_0_xyy_0_xxxy_0,\
                               t_0_xyy_0_xxxy_1, t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1,\
                               t_0_xyy_0_xyyy_0, t_0_xyy_0_xyyy_1, t_0_xyyy_0_xxxx_0,\
                               t_0_xyyy_0_xxxy_0, t_0_xyyy_0_xxxz_0, t_0_xyyy_0_xxyy_0,\
                               t_0_xyyy_0_xxyz_0, t_0_xyyy_0_xxzz_0, t_0_xyyy_0_xyyy_0,\
                               t_0_xyyy_0_xyyz_0, t_0_xyyy_0_xyzz_0, t_0_xyyy_0_xzzz_0,\
                               t_0_xyyy_0_yyyy_0, t_0_xyyy_0_yyyz_0, t_0_xyyy_0_yyzz_0,\
                               t_0_xyyy_0_yzzz_0, t_0_xyyy_0_zzzz_0, t_0_xyyz_0_xxxx_0,\
                               t_0_xyyz_0_xxxy_0, t_0_xyyz_0_xxxz_0, t_0_xyyz_0_xxyy_0,\
                               t_0_xyyz_0_xxyz_0, t_0_xyyz_0_xxzz_0, t_0_xyyz_0_xyyy_0,\
                               t_0_xyyz_0_xyyz_0, t_0_xyyz_0_xyzz_0, t_0_xyyz_0_xzzz_0,\
                               t_0_xyyz_0_yyyy_0, t_0_xyyz_0_yyyz_0, t_0_xzz_0_xyyz_0,\
                               t_0_xzz_0_xyyz_1, t_0_xzz_0_xyzz_0, t_0_xzz_0_xyzz_1,\
                               t_0_xzz_0_xzzz_0, t_0_xzz_0_xzzz_1, t_0_xzz_0_yyyy_0,\
                               t_0_xzz_0_yyyy_1, t_0_xzz_0_yyyz_0, t_0_xzz_0_yyyz_1,\
                               t_0_xzz_0_yyz_1, t_0_xzz_0_yyzz_0, t_0_xzz_0_yyzz_1,\
                               t_0_xzz_0_yzz_1, t_0_xzz_0_yzzz_0, t_0_xzz_0_yzzz_1,\
                               t_0_xzz_0_zzz_1, t_0_xzz_0_zzzz_0, t_0_xzz_0_zzzz_1,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxxx_0, t_0_yyy_0_xxxx_1,\
                               t_0_yyy_0_xxxy_0, t_0_yyy_0_xxxy_1, t_0_yyy_0_xxxz_0,\
                               t_0_yyy_0_xxxz_1, t_0_yyy_0_xxy_1, t_0_yyy_0_xxyy_0,\
                               t_0_yyy_0_xxyy_1, t_0_yyy_0_xxyz_0, t_0_yyy_0_xxyz_1,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xxzz_0, t_0_yyy_0_xxzz_1,\
                               t_0_yyy_0_xyy_1, t_0_yyy_0_xyyy_0, t_0_yyy_0_xyyy_1,\
                               t_0_yyy_0_xyyz_0, t_0_yyy_0_xyyz_1, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xyzz_0, t_0_yyy_0_xyzz_1, t_0_yyy_0_xzz_1,\
                               t_0_yyy_0_xzzz_0, t_0_yyy_0_xzzz_1, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyyy_0, t_0_yyy_0_yyyy_1, t_0_yyy_0_yyyz_0,\
                               t_0_yyy_0_yyyz_1, t_0_yyy_0_yyz_1, t_0_yyy_0_yyzz_0,\
                               t_0_yyy_0_yyzz_1, t_0_yyy_0_yzz_1, t_0_yyy_0_yzzz_0,\
                               t_0_yyy_0_yzzz_1, t_0_yyy_0_zzz_1, t_0_yyy_0_zzzz_0,\
                               t_0_yyy_0_zzzz_1, t_0_yyz_0_xxxz_0, t_0_yyz_0_xxxz_1,\
                               t_0_yyz_0_xxyz_0, t_0_yyz_0_xxyz_1, t_0_yyz_0_xxz_1,\
                               t_0_yyz_0_xxzz_0, t_0_yyz_0_xxzz_1, t_0_yyz_0_xyyz_0,\
                               t_0_yyz_0_xyyz_1, t_0_yyz_0_xyz_1, t_0_yyz_0_xyzz_0,\
                               t_0_yyz_0_xyzz_1, t_0_yyz_0_xzz_1, t_0_yyz_0_xzzz_0,\
                               t_0_yyz_0_xzzz_1, t_0_yyz_0_yyyy_0, t_0_yyz_0_yyyy_1,\
                               t_0_yyz_0_yyyz_0, t_0_yyz_0_yyyz_1, t_0_yyz_0_yyz_1,\
                               t_0_yyz_0_yzz_1, t_0_yyz_0_zzz_1, t_0_zz_0_xyyz_0, t_0_zz_0_xyyz_1,\
                               t_0_zz_0_xyzz_0, t_0_zz_0_xyzz_1, t_0_zz_0_xzzz_0, t_0_zz_0_xzzz_1,\
                               t_0_zz_0_yyyy_0, t_0_zz_0_yyyy_1, t_0_zz_0_yyyz_0, t_0_zz_0_yyyz_1,\
                               t_0_zz_0_yyzz_0, t_0_zz_0_yyzz_1, t_0_zz_0_yzzz_0, t_0_zz_0_yzzz_1,\
                               t_0_zz_0_zzzz_0, t_0_zz_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyyz_0_yyyz_0[i] = rpb_x[i] * t_0_yyz_0_yyyz_0[i] + rwp_x[i] * t_0_yyz_0_yyyz_1[i];

            t_0_xyyz_0_yyyy_0[i] = rpb_x[i] * t_0_yyz_0_yyyy_0[i] + rwp_x[i] * t_0_yyz_0_yyyy_1[i];

            t_0_xyyz_0_xzzz_0[i] = rpb_x[i] * t_0_yyz_0_xzzz_0[i] + rwp_x[i] * t_0_yyz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_zzz_1[i];

            t_0_xyyz_0_xyzz_0[i] = rpb_x[i] * t_0_yyz_0_xyzz_0[i] + rwp_x[i] * t_0_yyz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yzz_1[i];

            t_0_xyyz_0_xyyz_0[i] = rpb_x[i] * t_0_yyz_0_xyyz_0[i] + rwp_x[i] * t_0_yyz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyz_0_xyyy_0[i] = rpb_z[i] * t_0_xyy_0_xyyy_0[i] + rwp_z[i] * t_0_xyy_0_xyyy_1[i];

            t_0_xyyz_0_xxzz_0[i] = rpb_x[i] * t_0_yyz_0_xxzz_0[i] + rwp_x[i] * t_0_yyz_0_xxzz_1[i] + fze_0[i] * t_0_yyz_0_xzz_1[i];

            t_0_xyyz_0_xxyz_0[i] = rpb_x[i] * t_0_yyz_0_xxyz_0[i] + rwp_x[i] * t_0_yyz_0_xxyz_1[i] + fze_0[i] * t_0_yyz_0_xyz_1[i];

            t_0_xyyz_0_xxyy_0[i] = rpb_z[i] * t_0_xyy_0_xxyy_0[i] + rwp_z[i] * t_0_xyy_0_xxyy_1[i];

            t_0_xyyz_0_xxxz_0[i] = rpb_x[i] * t_0_yyz_0_xxxz_0[i] + rwp_x[i] * t_0_yyz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_yyz_0_xxz_1[i];

            t_0_xyyz_0_xxxy_0[i] = rpb_z[i] * t_0_xyy_0_xxxy_0[i] + rwp_z[i] * t_0_xyy_0_xxxy_1[i];

            t_0_xyyz_0_xxxx_0[i] = rpb_z[i] * t_0_xyy_0_xxxx_0[i] + rwp_z[i] * t_0_xyy_0_xxxx_1[i];

            t_0_xyyy_0_zzzz_0[i] = rpb_x[i] * t_0_yyy_0_zzzz_0[i] + rwp_x[i] * t_0_yyy_0_zzzz_1[i];

            t_0_xyyy_0_yzzz_0[i] = rpb_x[i] * t_0_yyy_0_yzzz_0[i] + rwp_x[i] * t_0_yyy_0_yzzz_1[i];

            t_0_xyyy_0_yyzz_0[i] = rpb_x[i] * t_0_yyy_0_yyzz_0[i] + rwp_x[i] * t_0_yyy_0_yyzz_1[i];

            t_0_xyyy_0_yyyz_0[i] = rpb_x[i] * t_0_yyy_0_yyyz_0[i] + rwp_x[i] * t_0_yyy_0_yyyz_1[i];

            t_0_xyyy_0_yyyy_0[i] = rpb_x[i] * t_0_yyy_0_yyyy_0[i] + rwp_x[i] * t_0_yyy_0_yyyy_1[i];

            t_0_xyyy_0_xzzz_0[i] = rpb_x[i] * t_0_yyy_0_xzzz_0[i] + rwp_x[i] * t_0_yyy_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zzz_1[i];

            t_0_xyyy_0_xyzz_0[i] = rpb_x[i] * t_0_yyy_0_xyzz_0[i] + rwp_x[i] * t_0_yyy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yzz_1[i];

            t_0_xyyy_0_xyyz_0[i] = rpb_x[i] * t_0_yyy_0_xyyz_0[i] + rwp_x[i] * t_0_yyy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyz_1[i];

            t_0_xyyy_0_xyyy_0[i] = rpb_x[i] * t_0_yyy_0_xyyy_0[i] + rwp_x[i] * t_0_yyy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yyy_1[i];

            t_0_xyyy_0_xxzz_0[i] = rpb_x[i] * t_0_yyy_0_xxzz_0[i] + rwp_x[i] * t_0_yyy_0_xxzz_1[i] + fze_0[i] * t_0_yyy_0_xzz_1[i];

            t_0_xyyy_0_xxyz_0[i] = rpb_x[i] * t_0_yyy_0_xxyz_0[i] + rwp_x[i] * t_0_yyy_0_xxyz_1[i] + fze_0[i] * t_0_yyy_0_xyz_1[i];

            t_0_xyyy_0_xxyy_0[i] = rpb_x[i] * t_0_yyy_0_xxyy_0[i] + rwp_x[i] * t_0_yyy_0_xxyy_1[i] + fze_0[i] * t_0_yyy_0_xyy_1[i];

            t_0_xyyy_0_xxxz_0[i] = rpb_x[i] * t_0_yyy_0_xxxz_0[i] + rwp_x[i] * t_0_yyy_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xxz_1[i];

            t_0_xyyy_0_xxxy_0[i] = rpb_x[i] * t_0_yyy_0_xxxy_0[i] + rwp_x[i] * t_0_yyy_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xxy_1[i];

            t_0_xyyy_0_xxxx_0[i] = rpb_x[i] * t_0_yyy_0_xxxx_0[i] + rwp_x[i] * t_0_yyy_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_yyy_0_xxx_1[i];

            t_0_xxzz_0_zzzz_0[i] = rpb_x[i] * t_0_xzz_0_zzzz_0[i] + rwp_x[i] * t_0_xzz_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzzz_1[i];

            t_0_xxzz_0_yzzz_0[i] = rpb_x[i] * t_0_xzz_0_yzzz_0[i] + rwp_x[i] * t_0_xzz_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzzz_1[i];

            t_0_xxzz_0_yyzz_0[i] = rpb_x[i] * t_0_xzz_0_yyzz_0[i] + rwp_x[i] * t_0_xzz_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyzz_1[i];

            t_0_xxzz_0_yyyz_0[i] = rpb_x[i] * t_0_xzz_0_yyyz_0[i] + rwp_x[i] * t_0_xzz_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyz_1[i];

            t_0_xxzz_0_yyyy_0[i] = rpb_x[i] * t_0_xzz_0_yyyy_0[i] + rwp_x[i] * t_0_xzz_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyyy_1[i];

            t_0_xxzz_0_xzzz_0[i] = rpb_x[i] * t_0_xzz_0_xzzz_0[i] + rwp_x[i] * t_0_xzz_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_zzz_1[i];

            t_0_xxzz_0_xyzz_0[i] = rpb_x[i] * t_0_xzz_0_xyzz_0[i] + rwp_x[i] * t_0_xzz_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yzz_1[i];

            t_0_xxzz_0_xyyz_0[i] = rpb_x[i] * t_0_xzz_0_xyyz_0[i] + rwp_x[i] * t_0_xzz_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yyz_1[i];

            t_0_xxzz_0_xyyy_0[i] = rpb_z[i] * t_0_xxz_0_xyyy_0[i] + rwp_z[i] * t_0_xxz_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xyyy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1, t_0_xx_0_xxxy_0,\
                               t_0_xx_0_xxxy_1, t_0_xx_0_xxxz_0, t_0_xx_0_xxxz_1, t_0_xx_0_xxyy_0,\
                               t_0_xx_0_xxyy_1, t_0_xx_0_xxzz_0, t_0_xx_0_xxzz_1, t_0_xx_0_xzzz_0,\
                               t_0_xx_0_xzzz_1, t_0_xxy_0_xxxx_0, t_0_xxy_0_xxxx_1,\
                               t_0_xxy_0_xxxy_0, t_0_xxy_0_xxxy_1, t_0_xxy_0_xxxz_0,\
                               t_0_xxy_0_xxxz_1, t_0_xxy_0_xxyy_0, t_0_xxy_0_xxyy_1,\
                               t_0_xxy_0_xxzz_0, t_0_xxy_0_xxzz_1, t_0_xxy_0_xyyy_0,\
                               t_0_xxy_0_xyyy_1, t_0_xxy_0_xzzz_0, t_0_xxy_0_xzzz_1,\
                               t_0_xxy_0_yyyy_0, t_0_xxy_0_yyyy_1, t_0_xxyy_0_xxxx_0,\
                               t_0_xxyy_0_xxxy_0, t_0_xxyy_0_xxxz_0, t_0_xxyy_0_xxyy_0,\
                               t_0_xxyy_0_xxyz_0, t_0_xxyy_0_xxzz_0, t_0_xxyy_0_xyyy_0,\
                               t_0_xxyy_0_xyyz_0, t_0_xxyy_0_xyzz_0, t_0_xxyy_0_xzzz_0,\
                               t_0_xxyy_0_yyyy_0, t_0_xxyy_0_yyyz_0, t_0_xxyy_0_yyzz_0,\
                               t_0_xxyy_0_yzzz_0, t_0_xxyy_0_zzzz_0, t_0_xxyz_0_xxxx_0,\
                               t_0_xxyz_0_xxxy_0, t_0_xxyz_0_xxxz_0, t_0_xxyz_0_xxyy_0,\
                               t_0_xxyz_0_xxyz_0, t_0_xxyz_0_xxzz_0, t_0_xxyz_0_xyyy_0,\
                               t_0_xxyz_0_xyyz_0, t_0_xxyz_0_xyzz_0, t_0_xxyz_0_xzzz_0,\
                               t_0_xxyz_0_yyyy_0, t_0_xxyz_0_yyyz_0, t_0_xxyz_0_yyzz_0,\
                               t_0_xxyz_0_yzzz_0, t_0_xxyz_0_zzzz_0, t_0_xxz_0_xxxx_0,\
                               t_0_xxz_0_xxxx_1, t_0_xxz_0_xxxy_0, t_0_xxz_0_xxxy_1,\
                               t_0_xxz_0_xxxz_0, t_0_xxz_0_xxxz_1, t_0_xxz_0_xxyy_0,\
                               t_0_xxz_0_xxyy_1, t_0_xxz_0_xxyz_0, t_0_xxz_0_xxyz_1,\
                               t_0_xxz_0_xxz_1, t_0_xxz_0_xxzz_0, t_0_xxz_0_xxzz_1,\
                               t_0_xxz_0_xyyz_0, t_0_xxz_0_xyyz_1, t_0_xxz_0_xyz_1,\
                               t_0_xxz_0_xyzz_0, t_0_xxz_0_xyzz_1, t_0_xxz_0_xzz_1,\
                               t_0_xxz_0_xzzz_0, t_0_xxz_0_xzzz_1, t_0_xxz_0_yyyz_0,\
                               t_0_xxz_0_yyyz_1, t_0_xxz_0_yyz_1, t_0_xxz_0_yyzz_0,\
                               t_0_xxz_0_yyzz_1, t_0_xxz_0_yzz_1, t_0_xxz_0_yzzz_0,\
                               t_0_xxz_0_yzzz_1, t_0_xxz_0_zzz_1, t_0_xxz_0_zzzz_0,\
                               t_0_xxz_0_zzzz_1, t_0_xxzz_0_xxxx_0, t_0_xxzz_0_xxxy_0,\
                               t_0_xxzz_0_xxxz_0, t_0_xxzz_0_xxyy_0, t_0_xxzz_0_xxyz_0,\
                               t_0_xxzz_0_xxzz_0, t_0_xyy_0_xxxy_0, t_0_xyy_0_xxxy_1,\
                               t_0_xyy_0_xxy_1, t_0_xyy_0_xxyy_0, t_0_xyy_0_xxyy_1,\
                               t_0_xyy_0_xxyz_0, t_0_xyy_0_xxyz_1, t_0_xyy_0_xyy_1,\
                               t_0_xyy_0_xyyy_0, t_0_xyy_0_xyyy_1, t_0_xyy_0_xyyz_0,\
                               t_0_xyy_0_xyyz_1, t_0_xyy_0_xyz_1, t_0_xyy_0_xyzz_0,\
                               t_0_xyy_0_xyzz_1, t_0_xyy_0_yyy_1, t_0_xyy_0_yyyy_0,\
                               t_0_xyy_0_yyyy_1, t_0_xyy_0_yyyz_0, t_0_xyy_0_yyyz_1,\
                               t_0_xyy_0_yyz_1, t_0_xyy_0_yyzz_0, t_0_xyy_0_yyzz_1,\
                               t_0_xyy_0_yzz_1, t_0_xyy_0_yzzz_0, t_0_xyy_0_yzzz_1,\
                               t_0_xyy_0_zzzz_0, t_0_xyy_0_zzzz_1, t_0_xzz_0_xxxz_0,\
                               t_0_xzz_0_xxxz_1, t_0_xzz_0_xxyz_0, t_0_xzz_0_xxyz_1,\
                               t_0_xzz_0_xxz_1, t_0_xzz_0_xxzz_0, t_0_xzz_0_xxzz_1,\
                               t_0_xzz_0_xyz_1, t_0_xzz_0_xzz_1, t_0_yy_0_xxxy_0, t_0_yy_0_xxxy_1,\
                               t_0_yy_0_xxyy_0, t_0_yy_0_xxyy_1, t_0_yy_0_xxyz_0, t_0_yy_0_xxyz_1,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyy_1, t_0_yy_0_xyyz_0, t_0_yy_0_xyyz_1,\
                               t_0_yy_0_xyzz_0, t_0_yy_0_xyzz_1, t_0_yy_0_yyyy_0, t_0_yy_0_yyyy_1,\
                               t_0_yy_0_yyyz_0, t_0_yy_0_yyyz_1, t_0_yy_0_yyzz_0, t_0_yy_0_yyzz_1,\
                               t_0_yy_0_yzzz_0, t_0_yy_0_yzzz_1, t_0_yy_0_zzzz_0, t_0_yy_0_zzzz_1,\
                               t_0_zz_0_xxxz_0, t_0_zz_0_xxxz_1, t_0_zz_0_xxyz_0, t_0_zz_0_xxyz_1,\
                               t_0_zz_0_xxzz_0, t_0_zz_0_xxzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxzz_0_xxzz_0[i] = rpb_x[i] * t_0_xzz_0_xxzz_0[i] + rwp_x[i] * t_0_xzz_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxzz_1[i] + fze_0[i] * t_0_xzz_0_xzz_1[i];

            t_0_xxzz_0_xxyz_0[i] = rpb_x[i] * t_0_xzz_0_xxyz_0[i] + rwp_x[i] * t_0_xzz_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxyz_1[i] + fze_0[i] * t_0_xzz_0_xyz_1[i];

            t_0_xxzz_0_xxyy_0[i] = rpb_z[i] * t_0_xxz_0_xxyy_0[i] + rwp_z[i] * t_0_xxz_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxyy_1[i];

            t_0_xxzz_0_xxxz_0[i] = rpb_x[i] * t_0_xzz_0_xxxz_0[i] + rwp_x[i] * t_0_xzz_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_xzz_0_xxz_1[i];

            t_0_xxzz_0_xxxy_0[i] = rpb_z[i] * t_0_xxz_0_xxxy_0[i] + rwp_z[i] * t_0_xxz_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxy_1[i];

            t_0_xxzz_0_xxxx_0[i] = rpb_z[i] * t_0_xxz_0_xxxx_0[i] + rwp_z[i] * t_0_xxz_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i];

            t_0_xxyz_0_zzzz_0[i] = rpb_y[i] * t_0_xxz_0_zzzz_0[i] + rwp_y[i] * t_0_xxz_0_zzzz_1[i];

            t_0_xxyz_0_yzzz_0[i] = rpb_y[i] * t_0_xxz_0_yzzz_0[i] + rwp_y[i] * t_0_xxz_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_zzz_1[i];

            t_0_xxyz_0_yyzz_0[i] = rpb_y[i] * t_0_xxz_0_yyzz_0[i] + rwp_y[i] * t_0_xxz_0_yyzz_1[i] + fze_0[i] * t_0_xxz_0_yzz_1[i];

            t_0_xxyz_0_yyyz_0[i] = rpb_y[i] * t_0_xxz_0_yyyz_0[i] + rwp_y[i] * t_0_xxz_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_xxz_0_yyz_1[i];

            t_0_xxyz_0_yyyy_0[i] = rpb_z[i] * t_0_xxy_0_yyyy_0[i] + rwp_z[i] * t_0_xxy_0_yyyy_1[i];

            t_0_xxyz_0_xzzz_0[i] = rpb_y[i] * t_0_xxz_0_xzzz_0[i] + rwp_y[i] * t_0_xxz_0_xzzz_1[i];

            t_0_xxyz_0_xyzz_0[i] = rpb_y[i] * t_0_xxz_0_xyzz_0[i] + rwp_y[i] * t_0_xxz_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xzz_1[i];

            t_0_xxyz_0_xyyz_0[i] = rpb_y[i] * t_0_xxz_0_xyyz_0[i] + rwp_y[i] * t_0_xxz_0_xyyz_1[i] + fze_0[i] * t_0_xxz_0_xyz_1[i];

            t_0_xxyz_0_xyyy_0[i] = rpb_z[i] * t_0_xxy_0_xyyy_0[i] + rwp_z[i] * t_0_xxy_0_xyyy_1[i];

            t_0_xxyz_0_xxzz_0[i] = rpb_y[i] * t_0_xxz_0_xxzz_0[i] + rwp_y[i] * t_0_xxz_0_xxzz_1[i];

            t_0_xxyz_0_xxyz_0[i] = rpb_y[i] * t_0_xxz_0_xxyz_0[i] + rwp_y[i] * t_0_xxz_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xxz_1[i];

            t_0_xxyz_0_xxyy_0[i] = rpb_z[i] * t_0_xxy_0_xxyy_0[i] + rwp_z[i] * t_0_xxy_0_xxyy_1[i];

            t_0_xxyz_0_xxxz_0[i] = rpb_y[i] * t_0_xxz_0_xxxz_0[i] + rwp_y[i] * t_0_xxz_0_xxxz_1[i];

            t_0_xxyz_0_xxxy_0[i] = rpb_z[i] * t_0_xxy_0_xxxy_0[i] + rwp_z[i] * t_0_xxy_0_xxxy_1[i];

            t_0_xxyz_0_xxxx_0[i] = rpb_y[i] * t_0_xxz_0_xxxx_0[i] + rwp_y[i] * t_0_xxz_0_xxxx_1[i];

            t_0_xxyy_0_zzzz_0[i] = rpb_x[i] * t_0_xyy_0_zzzz_0[i] + rwp_x[i] * t_0_xyy_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zzzz_1[i];

            t_0_xxyy_0_yzzz_0[i] = rpb_x[i] * t_0_xyy_0_yzzz_0[i] + rwp_x[i] * t_0_xyy_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yzzz_1[i];

            t_0_xxyy_0_yyzz_0[i] = rpb_x[i] * t_0_xyy_0_yyzz_0[i] + rwp_x[i] * t_0_xyy_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyzz_1[i];

            t_0_xxyy_0_yyyz_0[i] = rpb_x[i] * t_0_xyy_0_yyyz_0[i] + rwp_x[i] * t_0_xyy_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyz_1[i];

            t_0_xxyy_0_yyyy_0[i] = rpb_x[i] * t_0_xyy_0_yyyy_0[i] + rwp_x[i] * t_0_xyy_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyyy_1[i];

            t_0_xxyy_0_xzzz_0[i] = rpb_y[i] * t_0_xxy_0_xzzz_0[i] + rwp_y[i] * t_0_xxy_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xzzz_1[i];

            t_0_xxyy_0_xyzz_0[i] = rpb_x[i] * t_0_xyy_0_xyzz_0[i] + rwp_x[i] * t_0_xyy_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yzz_1[i];

            t_0_xxyy_0_xyyz_0[i] = rpb_x[i] * t_0_xyy_0_xyyz_0[i] + rwp_x[i] * t_0_xyy_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yyz_1[i];

            t_0_xxyy_0_xyyy_0[i] = rpb_x[i] * t_0_xyy_0_xyyy_0[i] + rwp_x[i] * t_0_xyy_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yyy_1[i];

            t_0_xxyy_0_xxzz_0[i] = rpb_y[i] * t_0_xxy_0_xxzz_0[i] + rwp_y[i] * t_0_xxy_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxzz_1[i];

            t_0_xxyy_0_xxyz_0[i] = rpb_x[i] * t_0_xyy_0_xxyz_0[i] + rwp_x[i] * t_0_xyy_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyz_1[i] + fze_0[i] * t_0_xyy_0_xyz_1[i];

            t_0_xxyy_0_xxyy_0[i] = rpb_x[i] * t_0_xyy_0_xxyy_0[i] + rwp_x[i] * t_0_xyy_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxyy_1[i] + fze_0[i] * t_0_xyy_0_xyy_1[i];

            t_0_xxyy_0_xxxz_0[i] = rpb_y[i] * t_0_xxy_0_xxxz_0[i] + rwp_y[i] * t_0_xxy_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxz_1[i];

            t_0_xxyy_0_xxxy_0[i] = rpb_x[i] * t_0_xyy_0_xxxy_0[i] + rwp_x[i] * t_0_xyy_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_xyy_0_xxy_1[i];

            t_0_xxyy_0_xxxx_0[i] = rpb_y[i] * t_0_xxy_0_xxxx_0[i] + rwp_y[i] * t_0_xxy_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xzzz_0, t_0_xx_0_xzzz_1, t_0_xx_0_yyyy_0,\
                               t_0_xx_0_yyyy_1, t_0_xx_0_yyyz_0, t_0_xx_0_yyyz_1, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yyzz_1, t_0_xx_0_yzzz_0, t_0_xx_0_yzzz_1, t_0_xx_0_zzzz_0,\
                               t_0_xx_0_zzzz_1, t_0_xxx_0_xxx_1, t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1,\
                               t_0_xxx_0_xxxy_0, t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0,\
                               t_0_xxx_0_xxxz_1, t_0_xxx_0_xxy_1, t_0_xxx_0_xxyy_0,\
                               t_0_xxx_0_xxyy_1, t_0_xxx_0_xxyz_0, t_0_xxx_0_xxyz_1,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xxzz_0, t_0_xxx_0_xxzz_1,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyyy_0, t_0_xxx_0_xyyy_1,\
                               t_0_xxx_0_xyyz_0, t_0_xxx_0_xyyz_1, t_0_xxx_0_xyz_1,\
                               t_0_xxx_0_xyzz_0, t_0_xxx_0_xyzz_1, t_0_xxx_0_xzz_1,\
                               t_0_xxx_0_xzzz_0, t_0_xxx_0_xzzz_1, t_0_xxx_0_yyy_1,\
                               t_0_xxx_0_yyyy_0, t_0_xxx_0_yyyy_1, t_0_xxx_0_yyyz_0,\
                               t_0_xxx_0_yyyz_1, t_0_xxx_0_yyz_1, t_0_xxx_0_yyzz_0,\
                               t_0_xxx_0_yyzz_1, t_0_xxx_0_yzz_1, t_0_xxx_0_yzzz_0,\
                               t_0_xxx_0_yzzz_1, t_0_xxx_0_zzz_1, t_0_xxx_0_zzzz_0,\
                               t_0_xxx_0_zzzz_1, t_0_xxxx_0_xzzz_0, t_0_xxxx_0_yyyy_0,\
                               t_0_xxxx_0_yyyz_0, t_0_xxxx_0_yyzz_0, t_0_xxxx_0_yzzz_0,\
                               t_0_xxxx_0_zzzz_0, t_0_xxxy_0_xxxx_0, t_0_xxxy_0_xxxy_0,\
                               t_0_xxxy_0_xxxz_0, t_0_xxxy_0_xxyy_0, t_0_xxxy_0_xxyz_0,\
                               t_0_xxxy_0_xxzz_0, t_0_xxxy_0_xyyy_0, t_0_xxxy_0_xyyz_0,\
                               t_0_xxxy_0_xyzz_0, t_0_xxxy_0_xzzz_0, t_0_xxxy_0_yyyy_0,\
                               t_0_xxxy_0_yyyz_0, t_0_xxxy_0_yyzz_0, t_0_xxxy_0_yzzz_0,\
                               t_0_xxxy_0_zzzz_0, t_0_xxxz_0_xxxx_0, t_0_xxxz_0_xxxy_0,\
                               t_0_xxxz_0_xxxz_0, t_0_xxxz_0_xxyy_0, t_0_xxxz_0_xxyz_0,\
                               t_0_xxxz_0_xxzz_0, t_0_xxxz_0_xyyy_0, t_0_xxxz_0_xyyz_0,\
                               t_0_xxxz_0_xyzz_0, t_0_xxxz_0_xzzz_0, t_0_xxxz_0_yyyy_0,\
                               t_0_xxxz_0_yyyz_0, t_0_xxxz_0_yyzz_0, t_0_xxxz_0_yzzz_0,\
                               t_0_xxxz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxz_0_zzzz_0[i] = rpb_z[i] * t_0_xxx_0_zzzz_0[i] + rwp_z[i] * t_0_xxx_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxz_0_yzzz_0[i] = rpb_z[i] * t_0_xxx_0_yzzz_0[i] + rwp_z[i] * t_0_xxx_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxz_0_yyzz_0[i] = rpb_z[i] * t_0_xxx_0_yyzz_0[i] + rwp_z[i] * t_0_xxx_0_yyzz_1[i] + fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxz_0_yyyz_0[i] = rpb_z[i] * t_0_xxx_0_yyyz_0[i] + rwp_z[i] * t_0_xxx_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxz_0_yyyy_0[i] = rpb_z[i] * t_0_xxx_0_yyyy_0[i] + rwp_z[i] * t_0_xxx_0_yyyy_1[i];

            t_0_xxxz_0_xzzz_0[i] = rpb_z[i] * t_0_xxx_0_xzzz_0[i] + rwp_z[i] * t_0_xxx_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxz_0_xyzz_0[i] = rpb_z[i] * t_0_xxx_0_xyzz_0[i] + rwp_z[i] * t_0_xxx_0_xyzz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxz_0_xyyz_0[i] = rpb_z[i] * t_0_xxx_0_xyyz_0[i] + rwp_z[i] * t_0_xxx_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxz_0_xyyy_0[i] = rpb_z[i] * t_0_xxx_0_xyyy_0[i] + rwp_z[i] * t_0_xxx_0_xyyy_1[i];

            t_0_xxxz_0_xxzz_0[i] = rpb_z[i] * t_0_xxx_0_xxzz_0[i] + rwp_z[i] * t_0_xxx_0_xxzz_1[i] + fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxz_0_xxyz_0[i] = rpb_z[i] * t_0_xxx_0_xxyz_0[i] + rwp_z[i] * t_0_xxx_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxz_0_xxyy_0[i] = rpb_z[i] * t_0_xxx_0_xxyy_0[i] + rwp_z[i] * t_0_xxx_0_xxyy_1[i];

            t_0_xxxz_0_xxxz_0[i] = rpb_z[i] * t_0_xxx_0_xxxz_0[i] + rwp_z[i] * t_0_xxx_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxz_0_xxxy_0[i] = rpb_z[i] * t_0_xxx_0_xxxy_0[i] + rwp_z[i] * t_0_xxx_0_xxxy_1[i];

            t_0_xxxz_0_xxxx_0[i] = rpb_z[i] * t_0_xxx_0_xxxx_0[i] + rwp_z[i] * t_0_xxx_0_xxxx_1[i];

            t_0_xxxy_0_zzzz_0[i] = rpb_y[i] * t_0_xxx_0_zzzz_0[i] + rwp_y[i] * t_0_xxx_0_zzzz_1[i];

            t_0_xxxy_0_yzzz_0[i] = rpb_y[i] * t_0_xxx_0_yzzz_0[i] + rwp_y[i] * t_0_xxx_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxy_0_yyzz_0[i] = rpb_y[i] * t_0_xxx_0_yyzz_0[i] + rwp_y[i] * t_0_xxx_0_yyzz_1[i] + fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxy_0_yyyz_0[i] = rpb_y[i] * t_0_xxx_0_yyyz_0[i] + rwp_y[i] * t_0_xxx_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxy_0_yyyy_0[i] = rpb_y[i] * t_0_xxx_0_yyyy_0[i] + rwp_y[i] * t_0_xxx_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxy_0_xzzz_0[i] = rpb_y[i] * t_0_xxx_0_xzzz_0[i] + rwp_y[i] * t_0_xxx_0_xzzz_1[i];

            t_0_xxxy_0_xyzz_0[i] = rpb_y[i] * t_0_xxx_0_xyzz_0[i] + rwp_y[i] * t_0_xxx_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxy_0_xyyz_0[i] = rpb_y[i] * t_0_xxx_0_xyyz_0[i] + rwp_y[i] * t_0_xxx_0_xyyz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxy_0_xyyy_0[i] = rpb_y[i] * t_0_xxx_0_xyyy_0[i] + rwp_y[i] * t_0_xxx_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxy_0_xxzz_0[i] = rpb_y[i] * t_0_xxx_0_xxzz_0[i] + rwp_y[i] * t_0_xxx_0_xxzz_1[i];

            t_0_xxxy_0_xxyz_0[i] = rpb_y[i] * t_0_xxx_0_xxyz_0[i] + rwp_y[i] * t_0_xxx_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxy_0_xxyy_0[i] = rpb_y[i] * t_0_xxx_0_xxyy_0[i] + rwp_y[i] * t_0_xxx_0_xxyy_1[i] + fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxy_0_xxxz_0[i] = rpb_y[i] * t_0_xxx_0_xxxz_0[i] + rwp_y[i] * t_0_xxx_0_xxxz_1[i];

            t_0_xxxy_0_xxxy_0[i] = rpb_y[i] * t_0_xxx_0_xxxy_0[i] + rwp_y[i] * t_0_xxx_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_xxxx_0[i] = rpb_y[i] * t_0_xxx_0_xxxx_0[i] + rwp_y[i] * t_0_xxx_0_xxxx_1[i];

            t_0_xxxx_0_zzzz_0[i] = rpb_x[i] * t_0_xxx_0_zzzz_0[i] + rwp_x[i] * t_0_xxx_0_zzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zzzz_1[i];

            t_0_xxxx_0_yzzz_0[i] = rpb_x[i] * t_0_xxx_0_yzzz_0[i] + rwp_x[i] * t_0_xxx_0_yzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yzzz_1[i];

            t_0_xxxx_0_yyzz_0[i] = rpb_x[i] * t_0_xxx_0_yyzz_0[i] + rwp_x[i] * t_0_xxx_0_yyzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyzz_1[i];

            t_0_xxxx_0_yyyz_0[i] = rpb_x[i] * t_0_xxx_0_yyyz_0[i] + rwp_x[i] * t_0_xxx_0_yyyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyyz_1[i];

            t_0_xxxx_0_yyyy_0[i] = rpb_x[i] * t_0_xxx_0_yyyy_0[i] + rwp_x[i] * t_0_xxx_0_yyyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyyy_1[i];

            t_0_xxxx_0_xzzz_0[i] = rpb_x[i] * t_0_xxx_0_xzzz_0[i] + rwp_x[i] * t_0_xxx_0_xzzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xzzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_xx_0_xxxx_0, t_0_xx_0_xxxx_1,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxy_1, t_0_xx_0_xxxz_0, t_0_xx_0_xxxz_1,\
                               t_0_xx_0_xxyy_0, t_0_xx_0_xxyy_1, t_0_xx_0_xxyz_0, t_0_xx_0_xxyz_1,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xxzz_1, t_0_xx_0_xyyy_0, t_0_xx_0_xyyy_1,\
                               t_0_xx_0_xyyz_0, t_0_xx_0_xyyz_1, t_0_xx_0_xyzz_0, t_0_xx_0_xyzz_1,\
                               t_0_xxx_0_xxx_1, t_0_xxx_0_xxxx_0, t_0_xxx_0_xxxx_1,\
                               t_0_xxx_0_xxxy_0, t_0_xxx_0_xxxy_1, t_0_xxx_0_xxxz_0,\
                               t_0_xxx_0_xxxz_1, t_0_xxx_0_xxy_1, t_0_xxx_0_xxyy_0,\
                               t_0_xxx_0_xxyy_1, t_0_xxx_0_xxyz_0, t_0_xxx_0_xxyz_1,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xxzz_0, t_0_xxx_0_xxzz_1,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyyy_0, t_0_xxx_0_xyyy_1,\
                               t_0_xxx_0_xyyz_0, t_0_xxx_0_xyyz_1, t_0_xxx_0_xyz_1,\
                               t_0_xxx_0_xyzz_0, t_0_xxx_0_xyzz_1, t_0_xxx_0_xzz_1,\
                               t_0_xxx_0_yyy_1, t_0_xxx_0_yyz_1, t_0_xxx_0_yzz_1, t_0_xxxx_0_xxxx_0,\
                               t_0_xxxx_0_xxxy_0, t_0_xxxx_0_xxxz_0, t_0_xxxx_0_xxyy_0,\
                               t_0_xxxx_0_xxyz_0, t_0_xxxx_0_xxzz_0, t_0_xxxx_0_xyyy_0,\
                               t_0_xxxx_0_xyyz_0, t_0_xxxx_0_xyzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxx_0_xyzz_0[i] = rpb_x[i] * t_0_xxx_0_xyzz_0[i] + rwp_x[i] * t_0_xxx_0_xyzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yzz_1[i];

            t_0_xxxx_0_xyyz_0[i] = rpb_x[i] * t_0_xxx_0_xyyz_0[i] + rwp_x[i] * t_0_xxx_0_xyyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyz_1[i];

            t_0_xxxx_0_xyyy_0[i] = rpb_x[i] * t_0_xxx_0_xyyy_0[i] + rwp_x[i] * t_0_xxx_0_xyyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxx_0_xxzz_0[i] = rpb_x[i] * t_0_xxx_0_xxzz_0[i] + rwp_x[i] * t_0_xxx_0_xxzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxzz_1[i] + fze_0[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxx_0_xxyz_0[i] = rpb_x[i] * t_0_xxx_0_xxyz_0[i] + rwp_x[i] * t_0_xxx_0_xxyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxyz_1[i] + fze_0[i] * t_0_xxx_0_xyz_1[i];

            t_0_xxxx_0_xxyy_0[i] = rpb_x[i] * t_0_xxx_0_xxyy_0[i] + rwp_x[i] * t_0_xxx_0_xxyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxyy_1[i] + fze_0[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxx_0_xxxz_0[i] = rpb_x[i] * t_0_xxx_0_xxxz_0[i] + rwp_x[i] * t_0_xxx_0_xxxz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxx_0_xxxy_0[i] = rpb_x[i] * t_0_xxx_0_xxxy_0[i] + rwp_x[i] * t_0_xxx_0_xxxy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxx_0_xxxx_0[i] = rpb_x[i] * t_0_xxx_0_xxxx_0[i] + rwp_x[i] * t_0_xxx_0_xxxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_xxx_0_xxx_1[i];
        }
    }
}


} // derirec namespace
