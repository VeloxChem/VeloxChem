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
compHostVRRForSGSF_V0(      BufferHostXY<T>&      intsBufferSGSF,
                      const BufferHostX<int32_t>& intsIndexesSGSF0,
                      const BufferHostXY<T>&      intsBufferSDSF0,
                      const BufferHostX<int32_t>& intsIndexesSDSF0,
                      const BufferHostXY<T>&      intsBufferSDSF1,
                      const BufferHostX<int32_t>& intsIndexesSDSF1,
                      const BufferHostXY<T>&      intsBufferSFSD1,
                      const BufferHostX<int32_t>& intsIndexesSFSD1,
                      const BufferHostXY<T>&      intsBufferSFSF0,
                      const BufferHostX<int32_t>& intsIndexesSFSF0,
                      const BufferHostXY<T>&      intsBufferSFSF1,
                      const BufferHostX<int32_t>& intsIndexesSFSF1,
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

    // set up [SGSF]^(0) integral components

    t_0_zzzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(0));

    t_0_zzzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(1));

    t_0_zzzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(2));

    t_0_zzzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(3));

    t_0_zzzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(4));

    t_0_zzzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(5));

    t_0_zzzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(6));

    t_0_zzzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(7));

    t_0_zzzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(8));

    t_0_zzzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(9));

    t_0_yzzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(10));

    t_0_yzzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(11));

    t_0_yzzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(12));

    t_0_yzzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(13));

    t_0_yzzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(14));

    t_0_yzzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(15));

    t_0_yzzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(16));

    t_0_yzzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(17));

    t_0_yzzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(18));

    t_0_yzzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(19));

    t_0_yyzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(20));

    t_0_yyzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(21));

    t_0_yyzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(22));

    t_0_yyzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(23));

    t_0_yyzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(24));

    t_0_yyzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(25));

    t_0_yyzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(26));

    t_0_yyzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(27));

    t_0_yyzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(28));

    t_0_yyzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(29));

    t_0_yyyz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(30));

    t_0_yyyz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(31));

    t_0_yyyz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(32));

    t_0_yyyz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(33));

    t_0_yyyz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(34));

    t_0_yyyz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(35));

    t_0_yyyz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(36));

    t_0_yyyz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(37));

    t_0_yyyz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(38));

    t_0_yyyz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(39));

    t_0_yyyy_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(40));

    t_0_yyyy_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(41));

    t_0_yyyy_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(42));

    t_0_yyyy_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(43));

    t_0_yyyy_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(44));

    t_0_yyyy_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(45));

    t_0_yyyy_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(46));

    t_0_yyyy_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(47));

    t_0_yyyy_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(48));

    t_0_yyyy_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(49));

    t_0_xzzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(50));

    t_0_xzzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(51));

    t_0_xzzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(52));

    t_0_xzzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(53));

    t_0_xzzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(54));

    t_0_xzzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(55));

    t_0_xzzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(56));

    t_0_xzzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(57));

    t_0_xzzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(58));

    t_0_xzzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(59));

    t_0_xyzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(60));

    t_0_xyzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(61));

    t_0_xyzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(62));

    t_0_xyzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(63));

    t_0_xyzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(64));

    t_0_xyzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(65));

    t_0_xyzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(66));

    t_0_xyzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(67));

    t_0_xyzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(68));

    t_0_xyzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(69));

    t_0_xyyz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(70));

    t_0_xyyz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(71));

    t_0_xyyz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(72));

    t_0_xyyz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(73));

    t_0_xyyz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(74));

    t_0_xyyz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(75));

    t_0_xyyz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(76));

    t_0_xyyz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(77));

    t_0_xyyz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(78));

    t_0_xyyz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(79));

    t_0_xyyy_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(80));

    t_0_xyyy_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(81));

    t_0_xyyy_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(82));

    t_0_xyyy_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(83));

    t_0_xyyy_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(84));

    t_0_xyyy_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(85));

    t_0_xyyy_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(86));

    t_0_xyyy_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(87));

    t_0_xyyy_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(88));

    t_0_xyyy_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(89));

    t_0_xxzz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(90));

    t_0_xxzz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(91));

    t_0_xxzz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(92));

    t_0_xxzz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(93));

    t_0_xxzz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(94));

    t_0_xxzz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(95));

    t_0_xxzz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(96));

    t_0_xxzz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(97));

    t_0_xxzz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(98));

    t_0_xxzz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(99));

    t_0_xxyz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(100));

    t_0_xxyz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(101));

    t_0_xxyz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(102));

    t_0_xxyz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(103));

    t_0_xxyz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(104));

    t_0_xxyz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(105));

    t_0_xxyz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(106));

    t_0_xxyz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(107));

    t_0_xxyz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(108));

    t_0_xxyz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(109));

    t_0_xxyy_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(110));

    t_0_xxyy_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(111));

    t_0_xxyy_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(112));

    t_0_xxyy_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(113));

    t_0_xxyy_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(114));

    t_0_xxyy_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(115));

    t_0_xxyy_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(116));

    t_0_xxyy_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(117));

    t_0_xxyy_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(118));

    t_0_xxyy_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(119));

    t_0_xxxz_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(120));

    t_0_xxxz_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(121));

    t_0_xxxz_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(122));

    t_0_xxxz_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(123));

    t_0_xxxz_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(124));

    t_0_xxxz_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(125));

    t_0_xxxz_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(126));

    t_0_xxxz_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(127));

    t_0_xxxz_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(128));

    t_0_xxxz_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(129));

    t_0_xxxy_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(130));

    t_0_xxxy_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(131));

    t_0_xxxy_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(132));

    t_0_xxxy_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(133));

    t_0_xxxy_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(134));

    t_0_xxxy_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(135));

    t_0_xxxy_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(136));

    t_0_xxxy_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(137));

    t_0_xxxy_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(138));

    t_0_xxxy_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(139));

    t_0_xxxx_0_zzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(140));

    t_0_xxxx_0_yzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(141));

    t_0_xxxx_0_yyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(142));

    t_0_xxxx_0_yyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(143));

    t_0_xxxx_0_xzz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(144));

    t_0_xxxx_0_xyz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(145));

    t_0_xxxx_0_xyy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(146));

    t_0_xxxx_0_xxz_0 = intsBufferSGSF0.data(intsIndexesSGSF0(147));

    t_0_xxxx_0_xxy_0 = intsBufferSGSF0.data(intsIndexesSGSF0(148));

    t_0_xxxx_0_xxx_0 = intsBufferSGSF0.data(intsIndexesSGSF0(149));

    // set up [SDSF]^(0) integral components

    t_0_zz_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(0));

    t_0_zz_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(1));

    t_0_zz_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(2));

    t_0_zz_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(3));

    t_0_zz_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(4));

    t_0_zz_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(5));

    t_0_zz_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(6));

    t_0_zz_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(7));

    t_0_zz_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(8));

    t_0_zz_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(9));

    t_0_yy_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(10));

    t_0_yy_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(11));

    t_0_yy_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(12));

    t_0_yy_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(13));

    t_0_yy_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(14));

    t_0_yy_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(15));

    t_0_yy_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(16));

    t_0_yy_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(17));

    t_0_yy_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(18));

    t_0_yy_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(19));

    t_0_xx_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(20));

    t_0_xx_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(21));

    t_0_xx_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(22));

    t_0_xx_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(23));

    t_0_xx_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(24));

    t_0_xx_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(25));

    t_0_xx_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(26));

    t_0_xx_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(27));

    t_0_xx_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(28));

    t_0_xx_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(29));

    // set up [SDSF]^(1) integral components

    t_0_zz_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(0));

    t_0_zz_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(1));

    t_0_zz_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(2));

    t_0_zz_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(3));

    t_0_zz_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(4));

    t_0_zz_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(5));

    t_0_zz_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(6));

    t_0_zz_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(7));

    t_0_zz_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(8));

    t_0_zz_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(9));

    t_0_yy_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(10));

    t_0_yy_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(11));

    t_0_yy_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(12));

    t_0_yy_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(13));

    t_0_yy_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(14));

    t_0_yy_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(15));

    t_0_yy_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(16));

    t_0_yy_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(17));

    t_0_yy_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(18));

    t_0_yy_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(19));

    t_0_xx_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(20));

    t_0_xx_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(21));

    t_0_xx_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(22));

    t_0_xx_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(23));

    t_0_xx_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(24));

    t_0_xx_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(25));

    t_0_xx_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(26));

    t_0_xx_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(27));

    t_0_xx_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(28));

    t_0_xx_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(29));

    // set up [SFSD]^(1) integral components

    t_0_zzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(0));

    t_0_zzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(1));

    t_0_zzz_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(2));

    t_0_zzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(3));

    t_0_zzz_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(4));

    t_0_zzz_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(5));

    t_0_yzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(6));

    t_0_yzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(7));

    t_0_yzz_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(8));

    t_0_yzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(9));

    t_0_yzz_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(10));

    t_0_yyz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(11));

    t_0_yyz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(12));

    t_0_yyz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(13));

    t_0_yyy_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(14));

    t_0_yyy_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(15));

    t_0_yyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(16));

    t_0_yyy_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(17));

    t_0_yyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(18));

    t_0_yyy_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(19));

    t_0_xzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(20));

    t_0_xzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(21));

    t_0_xzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(22));

    t_0_xyy_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(23));

    t_0_xyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(24));

    t_0_xyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(25));

    t_0_xxz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(26));

    t_0_xxz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(27));

    t_0_xxz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(28));

    t_0_xxx_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(29));

    t_0_xxx_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(30));

    t_0_xxx_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(31));

    t_0_xxx_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(32));

    t_0_xxx_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(33));

    t_0_xxx_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(34));

    // set up [SFSF]^(0) integral components

    t_0_zzz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(0));

    t_0_zzz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(1));

    t_0_zzz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(2));

    t_0_zzz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(3));

    t_0_zzz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(4));

    t_0_zzz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(5));

    t_0_zzz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(6));

    t_0_zzz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(7));

    t_0_zzz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(8));

    t_0_zzz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(9));

    t_0_yzz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(10));

    t_0_yzz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(11));

    t_0_yzz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(12));

    t_0_yzz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(13));

    t_0_yzz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(14));

    t_0_yzz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(15));

    t_0_yzz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(16));

    t_0_yzz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(17));

    t_0_yzz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(18));

    t_0_yzz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(19));

    t_0_yyz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(20));

    t_0_yyz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(21));

    t_0_yyz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(22));

    t_0_yyz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(23));

    t_0_yyz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(24));

    t_0_yyz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(25));

    t_0_yyz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(26));

    t_0_yyz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(27));

    t_0_yyz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(28));

    t_0_yyy_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(29));

    t_0_yyy_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(30));

    t_0_yyy_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(31));

    t_0_yyy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(32));

    t_0_yyy_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(33));

    t_0_yyy_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(34));

    t_0_yyy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(35));

    t_0_yyy_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(36));

    t_0_yyy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(37));

    t_0_yyy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(38));

    t_0_xzz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(39));

    t_0_xzz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(40));

    t_0_xzz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(41));

    t_0_xzz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(42));

    t_0_xzz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(43));

    t_0_xzz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(44));

    t_0_xzz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(45));

    t_0_xzz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(46));

    t_0_xyy_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(47));

    t_0_xyy_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(48));

    t_0_xyy_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(49));

    t_0_xyy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(50));

    t_0_xyy_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(51));

    t_0_xyy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(52));

    t_0_xyy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(53));

    t_0_xyy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(54));

    t_0_xxz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(55));

    t_0_xxz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(56));

    t_0_xxz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(57));

    t_0_xxz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(58));

    t_0_xxz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(59));

    t_0_xxz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(60));

    t_0_xxz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(61));

    t_0_xxz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(62));

    t_0_xxz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(63));

    t_0_xxy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(64));

    t_0_xxy_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(65));

    t_0_xxy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(66));

    t_0_xxy_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(67));

    t_0_xxy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(68));

    t_0_xxy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(69));

    t_0_xxx_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(70));

    t_0_xxx_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(71));

    t_0_xxx_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(72));

    t_0_xxx_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(73));

    t_0_xxx_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(74));

    t_0_xxx_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(75));

    t_0_xxx_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(76));

    t_0_xxx_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(77));

    t_0_xxx_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(78));

    t_0_xxx_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(79));

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

    t_0_yzz_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(19));

    t_0_yyz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(20));

    t_0_yyz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(21));

    t_0_yyz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(22));

    t_0_yyz_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(23));

    t_0_yyz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(24));

    t_0_yyz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(25));

    t_0_yyz_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(26));

    t_0_yyz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(27));

    t_0_yyz_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(28));

    t_0_yyy_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(29));

    t_0_yyy_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(30));

    t_0_yyy_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(31));

    t_0_yyy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(32));

    t_0_yyy_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(33));

    t_0_yyy_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(34));

    t_0_yyy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(35));

    t_0_yyy_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(36));

    t_0_yyy_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(37));

    t_0_yyy_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(38));

    t_0_xzz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(39));

    t_0_xzz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(40));

    t_0_xzz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(41));

    t_0_xzz_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(42));

    t_0_xzz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(43));

    t_0_xzz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(44));

    t_0_xzz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(45));

    t_0_xzz_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(46));

    t_0_xyy_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(47));

    t_0_xyy_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(48));

    t_0_xyy_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(49));

    t_0_xyy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(50));

    t_0_xyy_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(51));

    t_0_xyy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(52));

    t_0_xyy_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(53));

    t_0_xyy_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(54));

    t_0_xxz_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(55));

    t_0_xxz_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(56));

    t_0_xxz_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(57));

    t_0_xxz_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(58));

    t_0_xxz_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(59));

    t_0_xxz_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(60));

    t_0_xxz_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(61));

    t_0_xxz_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(62));

    t_0_xxz_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(63));

    t_0_xxy_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(64));

    t_0_xxy_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(65));

    t_0_xxy_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(66));

    t_0_xxy_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(67));

    t_0_xxy_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(68));

    t_0_xxy_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(69));

    t_0_xxx_0_zzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(70));

    t_0_xxx_0_yzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(71));

    t_0_xxx_0_yyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(72));

    t_0_xxx_0_yyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(73));

    t_0_xxx_0_xzz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(74));

    t_0_xxx_0_xyz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(75));

    t_0_xxx_0_xyy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(76));

    t_0_xxx_0_xxz_1 = intsBufferSFSF1.data(intsIndexesSFSF1(77));

    t_0_xxx_0_xxy_1 = intsBufferSFSF1.data(intsIndexesSFSF1(78));

    t_0_xxx_0_xxx_1 = intsBufferSFSF1.data(intsIndexesSFSF1(79));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_xxy_0,\
                               t_0_yy_0_xxy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0, t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1,\
                               t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1, t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1,\
                               t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0, t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1,\
                               t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1, t_0_yyyz_0_xyz_0, t_0_yyyz_0_xzz_0,\
                               t_0_yyyz_0_yyy_0, t_0_yyyz_0_yyz_0, t_0_yyyz_0_yzz_0,\
                               t_0_yyyz_0_zzz_0, t_0_yyz_0_xxy_0, t_0_yyz_0_xxy_1, t_0_yyz_0_xyy_0,\
                               t_0_yyz_0_xyy_1, t_0_yyz_0_yyy_0, t_0_yyz_0_yyy_1, t_0_yyzz_0_xxx_0,\
                               t_0_yyzz_0_xxy_0, t_0_yyzz_0_xxz_0, t_0_yyzz_0_xyy_0,\
                               t_0_yyzz_0_xyz_0, t_0_yyzz_0_xzz_0, t_0_yyzz_0_yyy_0,\
                               t_0_yyzz_0_yyz_0, t_0_yyzz_0_yzz_0, t_0_yyzz_0_zzz_0,\
                               t_0_yzz_0_xxx_0, t_0_yzz_0_xxx_1, t_0_yzz_0_xxz_0, t_0_yzz_0_xxz_1,\
                               t_0_yzz_0_xyz_0, t_0_yzz_0_xyz_1, t_0_yzz_0_xz_1, t_0_yzz_0_xzz_0,\
                               t_0_yzz_0_xzz_1, t_0_yzz_0_yyz_0, t_0_yzz_0_yyz_1, t_0_yzz_0_yz_1,\
                               t_0_yzz_0_yzz_0, t_0_yzz_0_yzz_1, t_0_yzz_0_zz_1, t_0_yzz_0_zzz_0,\
                               t_0_yzz_0_zzz_1, t_0_yzzz_0_xxx_0, t_0_yzzz_0_xxy_0,\
                               t_0_yzzz_0_xxz_0, t_0_yzzz_0_xyy_0, t_0_yzzz_0_xyz_0,\
                               t_0_yzzz_0_xzz_0, t_0_yzzz_0_yyy_0, t_0_yzzz_0_yyz_0,\
                               t_0_yzzz_0_yzz_0, t_0_yzzz_0_zzz_0, t_0_zz_0_xxx_0, t_0_zz_0_xxx_1,\
                               t_0_zz_0_xxy_0, t_0_zz_0_xxy_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1,\
                               t_0_zz_0_xyy_0, t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yyy_0, t_0_zz_0_yyy_1,\
                               t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yzz_0, t_0_zz_0_yzz_1,\
                               t_0_zz_0_zzz_0, t_0_zz_0_zzz_1, t_0_zzz_0_xx_1, t_0_zzz_0_xxx_0,\
                               t_0_zzz_0_xxx_1, t_0_zzz_0_xxy_0, t_0_zzz_0_xxy_1, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xxz_1, t_0_zzz_0_xy_1, t_0_zzz_0_xyy_0, t_0_zzz_0_xyy_1,\
                               t_0_zzz_0_xyz_0, t_0_zzz_0_xyz_1, t_0_zzz_0_xz_1, t_0_zzz_0_xzz_0,\
                               t_0_zzz_0_xzz_1, t_0_zzz_0_yy_1, t_0_zzz_0_yyy_0, t_0_zzz_0_yyy_1,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yyz_1, t_0_zzz_0_yz_1, t_0_zzz_0_yzz_0,\
                               t_0_zzz_0_yzz_1, t_0_zzz_0_zz_1, t_0_zzz_0_zzz_0, t_0_zzz_0_zzz_1,\
                               t_0_zzzz_0_xxx_0, t_0_zzzz_0_xxy_0, t_0_zzzz_0_xxz_0,\
                               t_0_zzzz_0_xyy_0, t_0_zzzz_0_xyz_0, t_0_zzzz_0_xzz_0,\
                               t_0_zzzz_0_yyy_0, t_0_zzzz_0_yyz_0, t_0_zzzz_0_yzz_0,\
                               t_0_zzzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzz_0[i] += rpb_z[i] * t_0_zzz_0_zzz_0[i] + rwp_z[i] * t_0_zzz_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_zzzz_0_yzz_0[i] += rpb_z[i] * t_0_zzz_0_yzz_0[i] + rwp_z[i] * t_0_zzz_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yzz_1[i] + fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_zzzz_0_yyz_0[i] += rpb_z[i] * t_0_zzz_0_yyz_0[i] + rwp_z[i] * t_0_zzz_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_zzzz_0_yyy_0[i] += rpb_z[i] * t_0_zzz_0_yyy_0[i] + rwp_z[i] * t_0_zzz_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyy_1[i];

            t_0_zzzz_0_xzz_0[i] += rpb_z[i] * t_0_zzz_0_xzz_0[i] + rwp_z[i] * t_0_zzz_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xzz_1[i] + fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_zzzz_0_xyz_0[i] += rpb_z[i] * t_0_zzz_0_xyz_0[i] + rwp_z[i] * t_0_zzz_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_zzzz_0_xyy_0[i] += rpb_z[i] * t_0_zzz_0_xyy_0[i] + rwp_z[i] * t_0_zzz_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyy_1[i];

            t_0_zzzz_0_xxz_0[i] += rpb_z[i] * t_0_zzz_0_xxz_0[i] + rwp_z[i] * t_0_zzz_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_zzzz_0_xxy_0[i] += rpb_z[i] * t_0_zzz_0_xxy_0[i] + rwp_z[i] * t_0_zzz_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxy_1[i];

            t_0_zzzz_0_xxx_0[i] += rpb_z[i] * t_0_zzz_0_xxx_0[i] + rwp_z[i] * t_0_zzz_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxx_1[i];

            t_0_yzzz_0_zzz_0[i] += rpb_y[i] * t_0_zzz_0_zzz_0[i] + rwp_y[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yzz_0[i] += rpb_y[i] * t_0_zzz_0_yzz_0[i] + rwp_y[i] * t_0_zzz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yyz_0[i] += rpb_y[i] * t_0_zzz_0_yyz_0[i] + rwp_y[i] * t_0_zzz_0_yyz_1[i] + fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_yzzz_0_yyy_0[i] += rpb_y[i] * t_0_zzz_0_yyy_0[i] + rwp_y[i] * t_0_zzz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_yzzz_0_xzz_0[i] += rpb_y[i] * t_0_zzz_0_xzz_0[i] + rwp_y[i] * t_0_zzz_0_xzz_1[i];

            t_0_yzzz_0_xyz_0[i] += rpb_y[i] * t_0_zzz_0_xyz_0[i] + rwp_y[i] * t_0_zzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_yzzz_0_xyy_0[i] += rpb_y[i] * t_0_zzz_0_xyy_0[i] + rwp_y[i] * t_0_zzz_0_xyy_1[i] + fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_yzzz_0_xxz_0[i] += rpb_y[i] * t_0_zzz_0_xxz_0[i] + rwp_y[i] * t_0_zzz_0_xxz_1[i];

            t_0_yzzz_0_xxy_0[i] += rpb_y[i] * t_0_zzz_0_xxy_0[i] + rwp_y[i] * t_0_zzz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_yzzz_0_xxx_0[i] += rpb_y[i] * t_0_zzz_0_xxx_0[i] + rwp_y[i] * t_0_zzz_0_xxx_1[i];

            t_0_yyzz_0_zzz_0[i] += rpb_y[i] * t_0_yzz_0_zzz_0[i] + rwp_y[i] * t_0_yzz_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzz_1[i];

            t_0_yyzz_0_yzz_0[i] += rpb_y[i] * t_0_yzz_0_yzz_0[i] + rwp_y[i] * t_0_yzz_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_zz_1[i];

            t_0_yyzz_0_yyz_0[i] += rpb_y[i] * t_0_yzz_0_yyz_0[i] + rwp_y[i] * t_0_yzz_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_yzz_0_yz_1[i];

            t_0_yyzz_0_yyy_0[i] += rpb_z[i] * t_0_yyz_0_yyy_0[i] + rwp_z[i] * t_0_yyz_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyy_1[i];

            t_0_yyzz_0_xzz_0[i] += rpb_y[i] * t_0_yzz_0_xzz_0[i] + rwp_y[i] * t_0_yzz_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzz_1[i];

            t_0_yyzz_0_xyz_0[i] += rpb_y[i] * t_0_yzz_0_xyz_0[i] + rwp_y[i] * t_0_yzz_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xz_1[i];

            t_0_yyzz_0_xyy_0[i] += rpb_z[i] * t_0_yyz_0_xyy_0[i] + rwp_z[i] * t_0_yyz_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyy_1[i];

            t_0_yyzz_0_xxz_0[i] += rpb_y[i] * t_0_yzz_0_xxz_0[i] + rwp_y[i] * t_0_yzz_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxz_1[i];

            t_0_yyzz_0_xxy_0[i] += rpb_z[i] * t_0_yyz_0_xxy_0[i] + rwp_z[i] * t_0_yyz_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxy_1[i];

            t_0_yyzz_0_xxx_0[i] += rpb_y[i] * t_0_yzz_0_xxx_0[i] + rwp_y[i] * t_0_yzz_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxx_1[i];

            t_0_yyyz_0_zzz_0[i] += rpb_z[i] * t_0_yyy_0_zzz_0[i] + rwp_z[i] * t_0_yyy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_yyyz_0_yzz_0[i] += rpb_z[i] * t_0_yyy_0_yzz_0[i] + rwp_z[i] * t_0_yyy_0_yzz_1[i] + fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_yyyz_0_yyz_0[i] += rpb_z[i] * t_0_yyy_0_yyz_0[i] + rwp_z[i] * t_0_yyy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyz_0_yyy_0[i] += rpb_z[i] * t_0_yyy_0_yyy_0[i] + rwp_z[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyz_0_xzz_0[i] += rpb_z[i] * t_0_yyy_0_xzz_0[i] + rwp_z[i] * t_0_yyy_0_xzz_1[i] + fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_yyyz_0_xyz_0[i] += rpb_z[i] * t_0_yyy_0_xyz_0[i] + rwp_z[i] * t_0_yyy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xyyz_0_yzz_0, t_0_xyyz_0_zzz_0, t_0_xyzz_0_xxx_0,\
                               t_0_xyzz_0_xxy_0, t_0_xyzz_0_xxz_0, t_0_xyzz_0_xyy_0,\
                               t_0_xyzz_0_xyz_0, t_0_xyzz_0_xzz_0, t_0_xyzz_0_yyy_0,\
                               t_0_xyzz_0_yyz_0, t_0_xyzz_0_yzz_0, t_0_xyzz_0_zzz_0,\
                               t_0_xzz_0_xxx_0, t_0_xzz_0_xxx_1, t_0_xzz_0_xxz_0, t_0_xzz_0_xxz_1,\
                               t_0_xzz_0_xzz_0, t_0_xzz_0_xzz_1, t_0_xzzz_0_xxx_0, t_0_xzzz_0_xxy_0,\
                               t_0_xzzz_0_xxz_0, t_0_xzzz_0_xyy_0, t_0_xzzz_0_xyz_0,\
                               t_0_xzzz_0_xzz_0, t_0_xzzz_0_yyy_0, t_0_xzzz_0_yyz_0,\
                               t_0_xzzz_0_yzz_0, t_0_xzzz_0_zzz_0, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xx_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxy_0, t_0_yyy_0_xxy_1, t_0_yyy_0_xxz_0,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyy_1,\
                               t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1, t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1, t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1, t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0,\
                               t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1, t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1,\
                               t_0_yyyy_0_xxx_0, t_0_yyyy_0_xxy_0, t_0_yyyy_0_xxz_0,\
                               t_0_yyyy_0_xyy_0, t_0_yyyy_0_xyz_0, t_0_yyyy_0_xzz_0,\
                               t_0_yyyy_0_yyy_0, t_0_yyyy_0_yyz_0, t_0_yyyy_0_yzz_0,\
                               t_0_yyyy_0_zzz_0, t_0_yyyz_0_xxx_0, t_0_yyyz_0_xxy_0,\
                               t_0_yyyz_0_xxz_0, t_0_yyyz_0_xyy_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_yzz_1, t_0_yyz_0_zzz_0, t_0_yyz_0_zzz_1, t_0_yzz_0_xxy_0,\
                               t_0_yzz_0_xxy_1, t_0_yzz_0_xy_1, t_0_yzz_0_xyy_0, t_0_yzz_0_xyy_1,\
                               t_0_yzz_0_xyz_0, t_0_yzz_0_xyz_1, t_0_yzz_0_yy_1, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyy_1, t_0_yzz_0_yyz_0, t_0_yzz_0_yyz_1, t_0_yzz_0_yz_1,\
                               t_0_yzz_0_yzz_0, t_0_yzz_0_yzz_1, t_0_yzz_0_zzz_0, t_0_yzz_0_zzz_1,\
                               t_0_zzz_0_xx_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxx_1, t_0_zzz_0_xxy_0,\
                               t_0_zzz_0_xxy_1, t_0_zzz_0_xxz_0, t_0_zzz_0_xxz_1, t_0_zzz_0_xy_1,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyy_1, t_0_zzz_0_xyz_0, t_0_zzz_0_xyz_1,\
                               t_0_zzz_0_xz_1, t_0_zzz_0_xzz_0, t_0_zzz_0_xzz_1, t_0_zzz_0_yy_1,\
                               t_0_zzz_0_yyy_0, t_0_zzz_0_yyy_1, t_0_zzz_0_yyz_0, t_0_zzz_0_yyz_1,\
                               t_0_zzz_0_yz_1, t_0_zzz_0_yzz_0, t_0_zzz_0_yzz_1, t_0_zzz_0_zz_1,\
                               t_0_zzz_0_zzz_0, t_0_zzz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyyz_0_xyy_0[i] += rpb_z[i] * t_0_yyy_0_xyy_0[i] + rwp_z[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyz_0_xxz_0[i] += rpb_z[i] * t_0_yyy_0_xxz_0[i] + rwp_z[i] * t_0_yyy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyz_0_xxy_0[i] += rpb_z[i] * t_0_yyy_0_xxy_0[i] + rwp_z[i] * t_0_yyy_0_xxy_1[i];

            t_0_yyyz_0_xxx_0[i] += rpb_z[i] * t_0_yyy_0_xxx_0[i] + rwp_z[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyy_0_zzz_0[i] += rpb_y[i] * t_0_yyy_0_zzz_0[i] + rwp_y[i] * t_0_yyy_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zzz_1[i];

            t_0_yyyy_0_yzz_0[i] += rpb_y[i] * t_0_yyy_0_yzz_0[i] + rwp_y[i] * t_0_yyy_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_yyyy_0_yyz_0[i] += rpb_y[i] * t_0_yyy_0_yyz_0[i] + rwp_y[i] * t_0_yyy_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyz_1[i] + fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_yyyy_0_yyy_0[i] += rpb_y[i] * t_0_yyy_0_yyy_0[i] + rwp_y[i] * t_0_yyy_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyy_0_xzz_0[i] += rpb_y[i] * t_0_yyy_0_xzz_0[i] + rwp_y[i] * t_0_yyy_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xzz_1[i];

            t_0_yyyy_0_xyz_0[i] += rpb_y[i] * t_0_yyy_0_xyz_0[i] + rwp_y[i] * t_0_yyy_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_yyyy_0_xyy_0[i] += rpb_y[i] * t_0_yyy_0_xyy_0[i] + rwp_y[i] * t_0_yyy_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyy_1[i] + fze_0[i] * t_0_yyy_0_xy_1[i];

            t_0_yyyy_0_xxz_0[i] += rpb_y[i] * t_0_yyy_0_xxz_0[i] + rwp_y[i] * t_0_yyy_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxz_1[i];

            t_0_yyyy_0_xxy_0[i] += rpb_y[i] * t_0_yyy_0_xxy_0[i] + rwp_y[i] * t_0_yyy_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyy_0_xxx_0[i] += rpb_y[i] * t_0_yyy_0_xxx_0[i] + rwp_y[i] * t_0_yyy_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxx_1[i];

            t_0_xzzz_0_zzz_0[i] += rpb_x[i] * t_0_zzz_0_zzz_0[i] + rwp_x[i] * t_0_zzz_0_zzz_1[i];

            t_0_xzzz_0_yzz_0[i] += rpb_x[i] * t_0_zzz_0_yzz_0[i] + rwp_x[i] * t_0_zzz_0_yzz_1[i];

            t_0_xzzz_0_yyz_0[i] += rpb_x[i] * t_0_zzz_0_yyz_0[i] + rwp_x[i] * t_0_zzz_0_yyz_1[i];

            t_0_xzzz_0_yyy_0[i] += rpb_x[i] * t_0_zzz_0_yyy_0[i] + rwp_x[i] * t_0_zzz_0_yyy_1[i];

            t_0_xzzz_0_xzz_0[i] += rpb_x[i] * t_0_zzz_0_xzz_0[i] + rwp_x[i] * t_0_zzz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_xyz_0[i] += rpb_x[i] * t_0_zzz_0_xyz_0[i] + rwp_x[i] * t_0_zzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_xzzz_0_xyy_0[i] += rpb_x[i] * t_0_zzz_0_xyy_0[i] + rwp_x[i] * t_0_zzz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_xzzz_0_xxz_0[i] += rpb_x[i] * t_0_zzz_0_xxz_0[i] + rwp_x[i] * t_0_zzz_0_xxz_1[i] + fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_xzzz_0_xxy_0[i] += rpb_x[i] * t_0_zzz_0_xxy_0[i] + rwp_x[i] * t_0_zzz_0_xxy_1[i] + fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_xzzz_0_xxx_0[i] += rpb_x[i] * t_0_zzz_0_xxx_0[i] + rwp_x[i] * t_0_zzz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_xyzz_0_zzz_0[i] += rpb_x[i] * t_0_yzz_0_zzz_0[i] + rwp_x[i] * t_0_yzz_0_zzz_1[i];

            t_0_xyzz_0_yzz_0[i] += rpb_x[i] * t_0_yzz_0_yzz_0[i] + rwp_x[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyzz_0_yyz_0[i] += rpb_x[i] * t_0_yzz_0_yyz_0[i] + rwp_x[i] * t_0_yzz_0_yyz_1[i];

            t_0_xyzz_0_yyy_0[i] += rpb_x[i] * t_0_yzz_0_yyy_0[i] + rwp_x[i] * t_0_yzz_0_yyy_1[i];

            t_0_xyzz_0_xzz_0[i] += rpb_y[i] * t_0_xzz_0_xzz_0[i] + rwp_y[i] * t_0_xzz_0_xzz_1[i];

            t_0_xyzz_0_xyz_0[i] += rpb_x[i] * t_0_yzz_0_xyz_0[i] + rwp_x[i] * t_0_yzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_xyy_0[i] += rpb_x[i] * t_0_yzz_0_xyy_0[i] + rwp_x[i] * t_0_yzz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yy_1[i];

            t_0_xyzz_0_xxz_0[i] += rpb_y[i] * t_0_xzz_0_xxz_0[i] + rwp_y[i] * t_0_xzz_0_xxz_1[i];

            t_0_xyzz_0_xxy_0[i] += rpb_x[i] * t_0_yzz_0_xxy_0[i] + rwp_x[i] * t_0_yzz_0_xxy_1[i] + fze_0[i] * t_0_yzz_0_xy_1[i];

            t_0_xyzz_0_xxx_0[i] += rpb_y[i] * t_0_xzz_0_xxx_0[i] + rwp_y[i] * t_0_xzz_0_xxx_1[i];

            t_0_xyyz_0_zzz_0[i] += rpb_x[i] * t_0_yyz_0_zzz_0[i] + rwp_x[i] * t_0_yyz_0_zzz_1[i];

            t_0_xyyz_0_yzz_0[i] += rpb_x[i] * t_0_yyz_0_yzz_0[i] + rwp_x[i] * t_0_yyz_0_yzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xyy_1, t_0_xxy_0_yyy_0, t_0_xxy_0_yyy_1, t_0_xxyz_0_xxz_0,\
                               t_0_xxyz_0_xyy_0, t_0_xxyz_0_xyz_0, t_0_xxyz_0_xzz_0,\
                               t_0_xxyz_0_yyy_0, t_0_xxyz_0_yyz_0, t_0_xxyz_0_yzz_0,\
                               t_0_xxyz_0_zzz_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxx_1, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxy_1, t_0_xxz_0_xxz_0, t_0_xxz_0_xxz_1, t_0_xxz_0_xyy_0,\
                               t_0_xxz_0_xyy_1, t_0_xxz_0_xyz_0, t_0_xxz_0_xyz_1, t_0_xxz_0_xz_1,\
                               t_0_xxz_0_xzz_0, t_0_xxz_0_xzz_1, t_0_xxz_0_yyz_0, t_0_xxz_0_yyz_1,\
                               t_0_xxz_0_yz_1, t_0_xxz_0_yzz_0, t_0_xxz_0_yzz_1, t_0_xxz_0_zz_1,\
                               t_0_xxz_0_zzz_0, t_0_xxz_0_zzz_1, t_0_xxzz_0_xxx_0, t_0_xxzz_0_xxy_0,\
                               t_0_xxzz_0_xxz_0, t_0_xxzz_0_xyy_0, t_0_xxzz_0_xyz_0,\
                               t_0_xxzz_0_xzz_0, t_0_xxzz_0_yyy_0, t_0_xxzz_0_yyz_0,\
                               t_0_xxzz_0_yzz_0, t_0_xxzz_0_zzz_0, t_0_xyy_0_xxx_0,\
                               t_0_xyy_0_xxx_1, t_0_xyy_0_xxy_0, t_0_xyy_0_xxy_1, t_0_xyy_0_xyy_0,\
                               t_0_xyy_0_xyy_1, t_0_xyyy_0_xxx_0, t_0_xyyy_0_xxy_0,\
                               t_0_xyyy_0_xxz_0, t_0_xyyy_0_xyy_0, t_0_xyyy_0_xyz_0,\
                               t_0_xyyy_0_xzz_0, t_0_xyyy_0_yyy_0, t_0_xyyy_0_yyz_0,\
                               t_0_xyyy_0_yzz_0, t_0_xyyy_0_zzz_0, t_0_xyyz_0_xxx_0,\
                               t_0_xyyz_0_xxy_0, t_0_xyyz_0_xxz_0, t_0_xyyz_0_xyy_0,\
                               t_0_xyyz_0_xyz_0, t_0_xyyz_0_xzz_0, t_0_xyyz_0_yyy_0,\
                               t_0_xyyz_0_yyz_0, t_0_xzz_0_xxz_0, t_0_xzz_0_xxz_1, t_0_xzz_0_xyz_0,\
                               t_0_xzz_0_xyz_1, t_0_xzz_0_xz_1, t_0_xzz_0_xzz_0, t_0_xzz_0_xzz_1,\
                               t_0_xzz_0_yyy_0, t_0_xzz_0_yyy_1, t_0_xzz_0_yyz_0, t_0_xzz_0_yyz_1,\
                               t_0_xzz_0_yz_1, t_0_xzz_0_yzz_0, t_0_xzz_0_yzz_1, t_0_xzz_0_zz_1,\
                               t_0_xzz_0_zzz_0, t_0_xzz_0_zzz_1, t_0_yyy_0_xx_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxy_0, t_0_yyy_0_xxy_1, t_0_yyy_0_xxz_0,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyy_1,\
                               t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1, t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1, t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1, t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0,\
                               t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1, t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1,\
                               t_0_yyz_0_xxz_0, t_0_yyz_0_xxz_1, t_0_yyz_0_xyz_0, t_0_yyz_0_xyz_1,\
                               t_0_yyz_0_xz_1, t_0_yyz_0_xzz_0, t_0_yyz_0_xzz_1, t_0_yyz_0_yyy_0,\
                               t_0_yyz_0_yyy_1, t_0_yyz_0_yyz_0, t_0_yyz_0_yyz_1, t_0_yyz_0_yz_1,\
                               t_0_yyz_0_zz_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xyz_0,\
                               t_0_zz_0_xyz_1, t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yzz_0,\
                               t_0_zz_0_yzz_1, t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyyz_0_yyz_0[i] += rpb_x[i] * t_0_yyz_0_yyz_0[i] + rwp_x[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyz_0_yyy_0[i] += rpb_x[i] * t_0_yyz_0_yyy_0[i] + rwp_x[i] * t_0_yyz_0_yyy_1[i];

            t_0_xyyz_0_xzz_0[i] += rpb_x[i] * t_0_yyz_0_xzz_0[i] + rwp_x[i] * t_0_yyz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_zz_1[i];

            t_0_xyyz_0_xyz_0[i] += rpb_x[i] * t_0_yyz_0_xyz_0[i] + rwp_x[i] * t_0_yyz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_xyy_0[i] += rpb_z[i] * t_0_xyy_0_xyy_0[i] + rwp_z[i] * t_0_xyy_0_xyy_1[i];

            t_0_xyyz_0_xxz_0[i] += rpb_x[i] * t_0_yyz_0_xxz_0[i] + rwp_x[i] * t_0_yyz_0_xxz_1[i] + fze_0[i] * t_0_yyz_0_xz_1[i];

            t_0_xyyz_0_xxy_0[i] += rpb_z[i] * t_0_xyy_0_xxy_0[i] + rwp_z[i] * t_0_xyy_0_xxy_1[i];

            t_0_xyyz_0_xxx_0[i] += rpb_z[i] * t_0_xyy_0_xxx_0[i] + rwp_z[i] * t_0_xyy_0_xxx_1[i];

            t_0_xyyy_0_zzz_0[i] += rpb_x[i] * t_0_yyy_0_zzz_0[i] + rwp_x[i] * t_0_yyy_0_zzz_1[i];

            t_0_xyyy_0_yzz_0[i] += rpb_x[i] * t_0_yyy_0_yzz_0[i] + rwp_x[i] * t_0_yyy_0_yzz_1[i];

            t_0_xyyy_0_yyz_0[i] += rpb_x[i] * t_0_yyy_0_yyz_0[i] + rwp_x[i] * t_0_yyy_0_yyz_1[i];

            t_0_xyyy_0_yyy_0[i] += rpb_x[i] * t_0_yyy_0_yyy_0[i] + rwp_x[i] * t_0_yyy_0_yyy_1[i];

            t_0_xyyy_0_xzz_0[i] += rpb_x[i] * t_0_yyy_0_xzz_0[i] + rwp_x[i] * t_0_yyy_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_xyyy_0_xyz_0[i] += rpb_x[i] * t_0_yyy_0_xyz_0[i] + rwp_x[i] * t_0_yyy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_xyyy_0_xyy_0[i] += rpb_x[i] * t_0_yyy_0_xyy_0[i] + rwp_x[i] * t_0_yyy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xxz_0[i] += rpb_x[i] * t_0_yyy_0_xxz_0[i] + rwp_x[i] * t_0_yyy_0_xxz_1[i] + fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_xyyy_0_xxy_0[i] += rpb_x[i] * t_0_yyy_0_xxy_0[i] + rwp_x[i] * t_0_yyy_0_xxy_1[i] + fze_0[i] * t_0_yyy_0_xy_1[i];

            t_0_xyyy_0_xxx_0[i] += rpb_x[i] * t_0_yyy_0_xxx_0[i] + rwp_x[i] * t_0_yyy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_xxzz_0_zzz_0[i] += rpb_x[i] * t_0_xzz_0_zzz_0[i] + rwp_x[i] * t_0_xzz_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzz_1[i];

            t_0_xxzz_0_yzz_0[i] += rpb_x[i] * t_0_xzz_0_yzz_0[i] + rwp_x[i] * t_0_xzz_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzz_1[i];

            t_0_xxzz_0_yyz_0[i] += rpb_x[i] * t_0_xzz_0_yyz_0[i] + rwp_x[i] * t_0_xzz_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyz_1[i];

            t_0_xxzz_0_yyy_0[i] += rpb_x[i] * t_0_xzz_0_yyy_0[i] + rwp_x[i] * t_0_xzz_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyy_1[i];

            t_0_xxzz_0_xzz_0[i] += rpb_x[i] * t_0_xzz_0_xzz_0[i] + rwp_x[i] * t_0_xzz_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_zz_1[i];

            t_0_xxzz_0_xyz_0[i] += rpb_x[i] * t_0_xzz_0_xyz_0[i] + rwp_x[i] * t_0_xzz_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yz_1[i];

            t_0_xxzz_0_xyy_0[i] += rpb_z[i] * t_0_xxz_0_xyy_0[i] + rwp_z[i] * t_0_xxz_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xyy_1[i];

            t_0_xxzz_0_xxz_0[i] += rpb_x[i] * t_0_xzz_0_xxz_0[i] + rwp_x[i] * t_0_xzz_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_xzz_0_xz_1[i];

            t_0_xxzz_0_xxy_0[i] += rpb_z[i] * t_0_xxz_0_xxy_0[i] + rwp_z[i] * t_0_xxz_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxy_1[i];

            t_0_xxzz_0_xxx_0[i] += rpb_z[i] * t_0_xxz_0_xxx_0[i] + rwp_z[i] * t_0_xxz_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxx_1[i];

            t_0_xxyz_0_zzz_0[i] += rpb_y[i] * t_0_xxz_0_zzz_0[i] + rwp_y[i] * t_0_xxz_0_zzz_1[i];

            t_0_xxyz_0_yzz_0[i] += rpb_y[i] * t_0_xxz_0_yzz_0[i] + rwp_y[i] * t_0_xxz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_zz_1[i];

            t_0_xxyz_0_yyz_0[i] += rpb_y[i] * t_0_xxz_0_yyz_0[i] + rwp_y[i] * t_0_xxz_0_yyz_1[i] + fze_0[i] * t_0_xxz_0_yz_1[i];

            t_0_xxyz_0_yyy_0[i] += rpb_z[i] * t_0_xxy_0_yyy_0[i] + rwp_z[i] * t_0_xxy_0_yyy_1[i];

            t_0_xxyz_0_xzz_0[i] += rpb_y[i] * t_0_xxz_0_xzz_0[i] + rwp_y[i] * t_0_xxz_0_xzz_1[i];

            t_0_xxyz_0_xyz_0[i] += rpb_y[i] * t_0_xxz_0_xyz_0[i] + rwp_y[i] * t_0_xxz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xyy_0[i] += rpb_z[i] * t_0_xxy_0_xyy_0[i] + rwp_z[i] * t_0_xxy_0_xyy_1[i];

            t_0_xxyz_0_xxz_0[i] += rpb_y[i] * t_0_xxz_0_xxz_0[i] + rwp_y[i] * t_0_xxz_0_xxz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxz_0,\
                               t_0_xx_0_xxz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yyy_0,\
                               t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1, t_0_xx_0_yzz_0,\
                               t_0_xx_0_yzz_1, t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_xx_1,\
                               t_0_xxx_0_xxx_0, t_0_xxx_0_xxx_1, t_0_xxx_0_xxy_0, t_0_xxx_0_xxy_1,\
                               t_0_xxx_0_xxz_0, t_0_xxx_0_xxz_1, t_0_xxx_0_xy_1, t_0_xxx_0_xyy_0,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyz_0, t_0_xxx_0_xyz_1, t_0_xxx_0_xz_1,\
                               t_0_xxx_0_xzz_0, t_0_xxx_0_xzz_1, t_0_xxx_0_yy_1, t_0_xxx_0_yyy_0,\
                               t_0_xxx_0_yyy_1, t_0_xxx_0_yyz_0, t_0_xxx_0_yyz_1, t_0_xxx_0_yz_1,\
                               t_0_xxx_0_yzz_0, t_0_xxx_0_yzz_1, t_0_xxx_0_zz_1, t_0_xxx_0_zzz_0,\
                               t_0_xxx_0_zzz_1, t_0_xxxx_0_yyy_0, t_0_xxxx_0_yyz_0,\
                               t_0_xxxx_0_yzz_0, t_0_xxxx_0_zzz_0, t_0_xxxy_0_xxx_0,\
                               t_0_xxxy_0_xxy_0, t_0_xxxy_0_xxz_0, t_0_xxxy_0_xyy_0,\
                               t_0_xxxy_0_xyz_0, t_0_xxxy_0_xzz_0, t_0_xxxy_0_yyy_0,\
                               t_0_xxxy_0_yyz_0, t_0_xxxy_0_yzz_0, t_0_xxxy_0_zzz_0,\
                               t_0_xxxz_0_xxx_0, t_0_xxxz_0_xxy_0, t_0_xxxz_0_xxz_0,\
                               t_0_xxxz_0_xyy_0, t_0_xxxz_0_xyz_0, t_0_xxxz_0_xzz_0,\
                               t_0_xxxz_0_yyy_0, t_0_xxxz_0_yyz_0, t_0_xxxz_0_yzz_0,\
                               t_0_xxxz_0_zzz_0, t_0_xxy_0_xxx_0, t_0_xxy_0_xxx_1, t_0_xxy_0_xxy_0,\
                               t_0_xxy_0_xxy_1, t_0_xxy_0_xxz_0, t_0_xxy_0_xxz_1, t_0_xxy_0_xzz_0,\
                               t_0_xxy_0_xzz_1, t_0_xxyy_0_xxx_0, t_0_xxyy_0_xxy_0,\
                               t_0_xxyy_0_xxz_0, t_0_xxyy_0_xyy_0, t_0_xxyy_0_xyz_0,\
                               t_0_xxyy_0_xzz_0, t_0_xxyy_0_yyy_0, t_0_xxyy_0_yyz_0,\
                               t_0_xxyy_0_yzz_0, t_0_xxyy_0_zzz_0, t_0_xxyz_0_xxx_0,\
                               t_0_xxyz_0_xxy_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxx_1, t_0_xyy_0_xxy_0,\
                               t_0_xyy_0_xxy_1, t_0_xyy_0_xy_1, t_0_xyy_0_xyy_0, t_0_xyy_0_xyy_1,\
                               t_0_xyy_0_xyz_0, t_0_xyy_0_xyz_1, t_0_xyy_0_yy_1, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyy_1, t_0_xyy_0_yyz_0, t_0_xyy_0_yyz_1, t_0_xyy_0_yz_1,\
                               t_0_xyy_0_yzz_0, t_0_xyy_0_yzz_1, t_0_xyy_0_zzz_0, t_0_xyy_0_zzz_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1,\
                               t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zzz_0, t_0_yy_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxyz_0_xxy_0[i] += rpb_z[i] * t_0_xxy_0_xxy_0[i] + rwp_z[i] * t_0_xxy_0_xxy_1[i];

            t_0_xxyz_0_xxx_0[i] += rpb_y[i] * t_0_xxz_0_xxx_0[i] + rwp_y[i] * t_0_xxz_0_xxx_1[i];

            t_0_xxyy_0_zzz_0[i] += rpb_x[i] * t_0_xyy_0_zzz_0[i] + rwp_x[i] * t_0_xyy_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zzz_1[i];

            t_0_xxyy_0_yzz_0[i] += rpb_x[i] * t_0_xyy_0_yzz_0[i] + rwp_x[i] * t_0_xyy_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yzz_1[i];

            t_0_xxyy_0_yyz_0[i] += rpb_x[i] * t_0_xyy_0_yyz_0[i] + rwp_x[i] * t_0_xyy_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyz_1[i];

            t_0_xxyy_0_yyy_0[i] += rpb_x[i] * t_0_xyy_0_yyy_0[i] + rwp_x[i] * t_0_xyy_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyy_1[i];

            t_0_xxyy_0_xzz_0[i] += rpb_y[i] * t_0_xxy_0_xzz_0[i] + rwp_y[i] * t_0_xxy_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xzz_1[i];

            t_0_xxyy_0_xyz_0[i] += rpb_x[i] * t_0_xyy_0_xyz_0[i] + rwp_x[i] * t_0_xyy_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yz_1[i];

            t_0_xxyy_0_xyy_0[i] += rpb_x[i] * t_0_xyy_0_xyy_0[i] + rwp_x[i] * t_0_xyy_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yy_1[i];

            t_0_xxyy_0_xxz_0[i] += rpb_y[i] * t_0_xxy_0_xxz_0[i] + rwp_y[i] * t_0_xxy_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxz_1[i];

            t_0_xxyy_0_xxy_0[i] += rpb_x[i] * t_0_xyy_0_xxy_0[i] + rwp_x[i] * t_0_xyy_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_xyy_0_xy_1[i];

            t_0_xxyy_0_xxx_0[i] += rpb_y[i] * t_0_xxy_0_xxx_0[i] + rwp_y[i] * t_0_xxy_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxx_1[i];

            t_0_xxxz_0_zzz_0[i] += rpb_z[i] * t_0_xxx_0_zzz_0[i] + rwp_z[i] * t_0_xxx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxz_0_yzz_0[i] += rpb_z[i] * t_0_xxx_0_yzz_0[i] + rwp_z[i] * t_0_xxx_0_yzz_1[i] + fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxz_0_yyz_0[i] += rpb_z[i] * t_0_xxx_0_yyz_0[i] + rwp_z[i] * t_0_xxx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxz_0_yyy_0[i] += rpb_z[i] * t_0_xxx_0_yyy_0[i] + rwp_z[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxz_0_xzz_0[i] += rpb_z[i] * t_0_xxx_0_xzz_0[i] + rwp_z[i] * t_0_xxx_0_xzz_1[i] + fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxz_0_xyz_0[i] += rpb_z[i] * t_0_xxx_0_xyz_0[i] + rwp_z[i] * t_0_xxx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxz_0_xyy_0[i] += rpb_z[i] * t_0_xxx_0_xyy_0[i] + rwp_z[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxz_0_xxz_0[i] += rpb_z[i] * t_0_xxx_0_xxz_0[i] + rwp_z[i] * t_0_xxx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxz_0_xxy_0[i] += rpb_z[i] * t_0_xxx_0_xxy_0[i] + rwp_z[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxz_0_xxx_0[i] += rpb_z[i] * t_0_xxx_0_xxx_0[i] + rwp_z[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_zzz_0[i] += rpb_y[i] * t_0_xxx_0_zzz_0[i] + rwp_y[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxy_0_yzz_0[i] += rpb_y[i] * t_0_xxx_0_yzz_0[i] + rwp_y[i] * t_0_xxx_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxy_0_yyz_0[i] += rpb_y[i] * t_0_xxx_0_yyz_0[i] + rwp_y[i] * t_0_xxx_0_yyz_1[i] + fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxy_0_yyy_0[i] += rpb_y[i] * t_0_xxx_0_yyy_0[i] + rwp_y[i] * t_0_xxx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxy_0_xzz_0[i] += rpb_y[i] * t_0_xxx_0_xzz_0[i] + rwp_y[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxy_0_xyz_0[i] += rpb_y[i] * t_0_xxx_0_xyz_0[i] + rwp_y[i] * t_0_xxx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxy_0_xyy_0[i] += rpb_y[i] * t_0_xxx_0_xyy_0[i] + rwp_y[i] * t_0_xxx_0_xyy_1[i] + fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxy_0_xxz_0[i] += rpb_y[i] * t_0_xxx_0_xxz_0[i] + rwp_y[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxy_0_xxy_0[i] += rpb_y[i] * t_0_xxx_0_xxy_0[i] + rwp_y[i] * t_0_xxx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_xxx_0[i] += rpb_y[i] * t_0_xxx_0_xxx_0[i] + rwp_y[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxx_0_zzz_0[i] += rpb_x[i] * t_0_xxx_0_zzz_0[i] + rwp_x[i] * t_0_xxx_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zzz_1[i];

            t_0_xxxx_0_yzz_0[i] += rpb_x[i] * t_0_xxx_0_yzz_0[i] + rwp_x[i] * t_0_xxx_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yzz_1[i];

            t_0_xxxx_0_yyz_0[i] += rpb_x[i] * t_0_xxx_0_yyz_0[i] + rwp_x[i] * t_0_xxx_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyz_1[i];

            t_0_xxxx_0_yyy_0[i] += rpb_x[i] * t_0_xxx_0_yyy_0[i] + rwp_x[i] * t_0_xxx_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1,\
                               t_0_xx_0_xxy_0, t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xxx_0_xx_1, t_0_xxx_0_xxx_0,\
                               t_0_xxx_0_xxx_1, t_0_xxx_0_xxy_0, t_0_xxx_0_xxy_1, t_0_xxx_0_xxz_0,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xy_1, t_0_xxx_0_xyy_0, t_0_xxx_0_xyy_1,\
                               t_0_xxx_0_xyz_0, t_0_xxx_0_xyz_1, t_0_xxx_0_xz_1, t_0_xxx_0_xzz_0,\
                               t_0_xxx_0_xzz_1, t_0_xxx_0_yy_1, t_0_xxx_0_yz_1, t_0_xxx_0_zz_1,\
                               t_0_xxxx_0_xxx_0, t_0_xxxx_0_xxy_0, t_0_xxxx_0_xxz_0,\
                               t_0_xxxx_0_xyy_0, t_0_xxxx_0_xyz_0, t_0_xxxx_0_xzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxx_0_xzz_0[i] += rpb_x[i] * t_0_xxx_0_xzz_0[i] + rwp_x[i] * t_0_xxx_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxx_0_xyz_0[i] += rpb_x[i] * t_0_xxx_0_xyz_0[i] + rwp_x[i] * t_0_xxx_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxx_0_xyy_0[i] += rpb_x[i] * t_0_xxx_0_xyy_0[i] + rwp_x[i] * t_0_xxx_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxx_0_xxz_0[i] += rpb_x[i] * t_0_xxx_0_xxz_0[i] + rwp_x[i] * t_0_xxx_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxz_1[i] + fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxx_0_xxy_0[i] += rpb_x[i] * t_0_xxx_0_xxy_0[i] + rwp_x[i] * t_0_xxx_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxy_1[i] + fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxx_0_xxx_0[i] += rpb_x[i] * t_0_xxx_0_xxx_0[i] + rwp_x[i] * t_0_xxx_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_yy_0_xxy_0,\
                               t_0_yy_0_xxy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1,\
                               t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0, t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1,\
                               t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1, t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1,\
                               t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0, t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1,\
                               t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1, t_0_yyyz_0_xyz_0, t_0_yyyz_0_xzz_0,\
                               t_0_yyyz_0_yyy_0, t_0_yyyz_0_yyz_0, t_0_yyyz_0_yzz_0,\
                               t_0_yyyz_0_zzz_0, t_0_yyz_0_xxy_0, t_0_yyz_0_xxy_1, t_0_yyz_0_xyy_0,\
                               t_0_yyz_0_xyy_1, t_0_yyz_0_yyy_0, t_0_yyz_0_yyy_1, t_0_yyzz_0_xxx_0,\
                               t_0_yyzz_0_xxy_0, t_0_yyzz_0_xxz_0, t_0_yyzz_0_xyy_0,\
                               t_0_yyzz_0_xyz_0, t_0_yyzz_0_xzz_0, t_0_yyzz_0_yyy_0,\
                               t_0_yyzz_0_yyz_0, t_0_yyzz_0_yzz_0, t_0_yyzz_0_zzz_0,\
                               t_0_yzz_0_xxx_0, t_0_yzz_0_xxx_1, t_0_yzz_0_xxz_0, t_0_yzz_0_xxz_1,\
                               t_0_yzz_0_xyz_0, t_0_yzz_0_xyz_1, t_0_yzz_0_xz_1, t_0_yzz_0_xzz_0,\
                               t_0_yzz_0_xzz_1, t_0_yzz_0_yyz_0, t_0_yzz_0_yyz_1, t_0_yzz_0_yz_1,\
                               t_0_yzz_0_yzz_0, t_0_yzz_0_yzz_1, t_0_yzz_0_zz_1, t_0_yzz_0_zzz_0,\
                               t_0_yzz_0_zzz_1, t_0_yzzz_0_xxx_0, t_0_yzzz_0_xxy_0,\
                               t_0_yzzz_0_xxz_0, t_0_yzzz_0_xyy_0, t_0_yzzz_0_xyz_0,\
                               t_0_yzzz_0_xzz_0, t_0_yzzz_0_yyy_0, t_0_yzzz_0_yyz_0,\
                               t_0_yzzz_0_yzz_0, t_0_yzzz_0_zzz_0, t_0_zz_0_xxx_0, t_0_zz_0_xxx_1,\
                               t_0_zz_0_xxy_0, t_0_zz_0_xxy_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1,\
                               t_0_zz_0_xyy_0, t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yyy_0, t_0_zz_0_yyy_1,\
                               t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yzz_0, t_0_zz_0_yzz_1,\
                               t_0_zz_0_zzz_0, t_0_zz_0_zzz_1, t_0_zzz_0_xx_1, t_0_zzz_0_xxx_0,\
                               t_0_zzz_0_xxx_1, t_0_zzz_0_xxy_0, t_0_zzz_0_xxy_1, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xxz_1, t_0_zzz_0_xy_1, t_0_zzz_0_xyy_0, t_0_zzz_0_xyy_1,\
                               t_0_zzz_0_xyz_0, t_0_zzz_0_xyz_1, t_0_zzz_0_xz_1, t_0_zzz_0_xzz_0,\
                               t_0_zzz_0_xzz_1, t_0_zzz_0_yy_1, t_0_zzz_0_yyy_0, t_0_zzz_0_yyy_1,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yyz_1, t_0_zzz_0_yz_1, t_0_zzz_0_yzz_0,\
                               t_0_zzz_0_yzz_1, t_0_zzz_0_zz_1, t_0_zzz_0_zzz_0, t_0_zzz_0_zzz_1,\
                               t_0_zzzz_0_xxx_0, t_0_zzzz_0_xxy_0, t_0_zzzz_0_xxz_0,\
                               t_0_zzzz_0_xyy_0, t_0_zzzz_0_xyz_0, t_0_zzzz_0_xzz_0,\
                               t_0_zzzz_0_yyy_0, t_0_zzzz_0_yyz_0, t_0_zzzz_0_yzz_0,\
                               t_0_zzzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zzz_0[i] = rpb_z[i] * t_0_zzz_0_zzz_0[i] + rwp_z[i] * t_0_zzz_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_zzzz_0_yzz_0[i] = rpb_z[i] * t_0_zzz_0_yzz_0[i] + rwp_z[i] * t_0_zzz_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yzz_1[i] + fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_zzzz_0_yyz_0[i] = rpb_z[i] * t_0_zzz_0_yyz_0[i] + rwp_z[i] * t_0_zzz_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_zzzz_0_yyy_0[i] = rpb_z[i] * t_0_zzz_0_yyy_0[i] + rwp_z[i] * t_0_zzz_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yyy_1[i];

            t_0_zzzz_0_xzz_0[i] = rpb_z[i] * t_0_zzz_0_xzz_0[i] + rwp_z[i] * t_0_zzz_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xzz_1[i] + fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_zzzz_0_xyz_0[i] = rpb_z[i] * t_0_zzz_0_xyz_0[i] + rwp_z[i] * t_0_zzz_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_zzzz_0_xyy_0[i] = rpb_z[i] * t_0_zzz_0_xyy_0[i] + rwp_z[i] * t_0_zzz_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xyy_1[i];

            t_0_zzzz_0_xxz_0[i] = rpb_z[i] * t_0_zzz_0_xxz_0[i] + rwp_z[i] * t_0_zzz_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_zzzz_0_xxy_0[i] = rpb_z[i] * t_0_zzz_0_xxy_0[i] + rwp_z[i] * t_0_zzz_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxy_1[i];

            t_0_zzzz_0_xxx_0[i] = rpb_z[i] * t_0_zzz_0_xxx_0[i] + rwp_z[i] * t_0_zzz_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xxx_1[i];

            t_0_yzzz_0_zzz_0[i] = rpb_y[i] * t_0_zzz_0_zzz_0[i] + rwp_y[i] * t_0_zzz_0_zzz_1[i];

            t_0_yzzz_0_yzz_0[i] = rpb_y[i] * t_0_zzz_0_yzz_0[i] + rwp_y[i] * t_0_zzz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yyz_0[i] = rpb_y[i] * t_0_zzz_0_yyz_0[i] + rwp_y[i] * t_0_zzz_0_yyz_1[i] + fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_yzzz_0_yyy_0[i] = rpb_y[i] * t_0_zzz_0_yyy_0[i] + rwp_y[i] * t_0_zzz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_yzzz_0_xzz_0[i] = rpb_y[i] * t_0_zzz_0_xzz_0[i] + rwp_y[i] * t_0_zzz_0_xzz_1[i];

            t_0_yzzz_0_xyz_0[i] = rpb_y[i] * t_0_zzz_0_xyz_0[i] + rwp_y[i] * t_0_zzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_yzzz_0_xyy_0[i] = rpb_y[i] * t_0_zzz_0_xyy_0[i] + rwp_y[i] * t_0_zzz_0_xyy_1[i] + fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_yzzz_0_xxz_0[i] = rpb_y[i] * t_0_zzz_0_xxz_0[i] + rwp_y[i] * t_0_zzz_0_xxz_1[i];

            t_0_yzzz_0_xxy_0[i] = rpb_y[i] * t_0_zzz_0_xxy_0[i] + rwp_y[i] * t_0_zzz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_yzzz_0_xxx_0[i] = rpb_y[i] * t_0_zzz_0_xxx_0[i] + rwp_y[i] * t_0_zzz_0_xxx_1[i];

            t_0_yyzz_0_zzz_0[i] = rpb_y[i] * t_0_yzz_0_zzz_0[i] + rwp_y[i] * t_0_yzz_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzz_1[i];

            t_0_yyzz_0_yzz_0[i] = rpb_y[i] * t_0_yzz_0_yzz_0[i] + rwp_y[i] * t_0_yzz_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_zz_1[i];

            t_0_yyzz_0_yyz_0[i] = rpb_y[i] * t_0_yzz_0_yyz_0[i] + rwp_y[i] * t_0_yzz_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_yzz_0_yz_1[i];

            t_0_yyzz_0_yyy_0[i] = rpb_z[i] * t_0_yyz_0_yyy_0[i] + rwp_z[i] * t_0_yyz_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyy_1[i];

            t_0_yyzz_0_xzz_0[i] = rpb_y[i] * t_0_yzz_0_xzz_0[i] + rwp_y[i] * t_0_yzz_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzz_1[i];

            t_0_yyzz_0_xyz_0[i] = rpb_y[i] * t_0_yzz_0_xyz_0[i] + rwp_y[i] * t_0_yzz_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_xz_1[i];

            t_0_yyzz_0_xyy_0[i] = rpb_z[i] * t_0_yyz_0_xyy_0[i] + rwp_z[i] * t_0_yyz_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyy_1[i];

            t_0_yyzz_0_xxz_0[i] = rpb_y[i] * t_0_yzz_0_xxz_0[i] + rwp_y[i] * t_0_yzz_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxz_1[i];

            t_0_yyzz_0_xxy_0[i] = rpb_z[i] * t_0_yyz_0_xxy_0[i] + rwp_z[i] * t_0_yyz_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxy_1[i];

            t_0_yyzz_0_xxx_0[i] = rpb_y[i] * t_0_yzz_0_xxx_0[i] + rwp_y[i] * t_0_yzz_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxx_1[i];

            t_0_yyyz_0_zzz_0[i] = rpb_z[i] * t_0_yyy_0_zzz_0[i] + rwp_z[i] * t_0_yyy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_yyyz_0_yzz_0[i] = rpb_z[i] * t_0_yyy_0_yzz_0[i] + rwp_z[i] * t_0_yyy_0_yzz_1[i] + fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_yyyz_0_yyz_0[i] = rpb_z[i] * t_0_yyy_0_yyz_0[i] + rwp_z[i] * t_0_yyy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyz_0_yyy_0[i] = rpb_z[i] * t_0_yyy_0_yyy_0[i] + rwp_z[i] * t_0_yyy_0_yyy_1[i];

            t_0_yyyz_0_xzz_0[i] = rpb_z[i] * t_0_yyy_0_xzz_0[i] + rwp_z[i] * t_0_yyy_0_xzz_1[i] + fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_yyyz_0_xyz_0[i] = rpb_z[i] * t_0_yyy_0_xyz_0[i] + rwp_z[i] * t_0_yyy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xyyz_0_yzz_0, t_0_xyyz_0_zzz_0, t_0_xyzz_0_xxx_0,\
                               t_0_xyzz_0_xxy_0, t_0_xyzz_0_xxz_0, t_0_xyzz_0_xyy_0,\
                               t_0_xyzz_0_xyz_0, t_0_xyzz_0_xzz_0, t_0_xyzz_0_yyy_0,\
                               t_0_xyzz_0_yyz_0, t_0_xyzz_0_yzz_0, t_0_xyzz_0_zzz_0,\
                               t_0_xzz_0_xxx_0, t_0_xzz_0_xxx_1, t_0_xzz_0_xxz_0, t_0_xzz_0_xxz_1,\
                               t_0_xzz_0_xzz_0, t_0_xzz_0_xzz_1, t_0_xzzz_0_xxx_0, t_0_xzzz_0_xxy_0,\
                               t_0_xzzz_0_xxz_0, t_0_xzzz_0_xyy_0, t_0_xzzz_0_xyz_0,\
                               t_0_xzzz_0_xzz_0, t_0_xzzz_0_yyy_0, t_0_xzzz_0_yyz_0,\
                               t_0_xzzz_0_yzz_0, t_0_xzzz_0_zzz_0, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xx_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxy_0, t_0_yyy_0_xxy_1, t_0_yyy_0_xxz_0,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyy_1,\
                               t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1, t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1, t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1, t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0,\
                               t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1, t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1,\
                               t_0_yyyy_0_xxx_0, t_0_yyyy_0_xxy_0, t_0_yyyy_0_xxz_0,\
                               t_0_yyyy_0_xyy_0, t_0_yyyy_0_xyz_0, t_0_yyyy_0_xzz_0,\
                               t_0_yyyy_0_yyy_0, t_0_yyyy_0_yyz_0, t_0_yyyy_0_yzz_0,\
                               t_0_yyyy_0_zzz_0, t_0_yyyz_0_xxx_0, t_0_yyyz_0_xxy_0,\
                               t_0_yyyz_0_xxz_0, t_0_yyyz_0_xyy_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_yzz_1, t_0_yyz_0_zzz_0, t_0_yyz_0_zzz_1, t_0_yzz_0_xxy_0,\
                               t_0_yzz_0_xxy_1, t_0_yzz_0_xy_1, t_0_yzz_0_xyy_0, t_0_yzz_0_xyy_1,\
                               t_0_yzz_0_xyz_0, t_0_yzz_0_xyz_1, t_0_yzz_0_yy_1, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyy_1, t_0_yzz_0_yyz_0, t_0_yzz_0_yyz_1, t_0_yzz_0_yz_1,\
                               t_0_yzz_0_yzz_0, t_0_yzz_0_yzz_1, t_0_yzz_0_zzz_0, t_0_yzz_0_zzz_1,\
                               t_0_zzz_0_xx_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxx_1, t_0_zzz_0_xxy_0,\
                               t_0_zzz_0_xxy_1, t_0_zzz_0_xxz_0, t_0_zzz_0_xxz_1, t_0_zzz_0_xy_1,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyy_1, t_0_zzz_0_xyz_0, t_0_zzz_0_xyz_1,\
                               t_0_zzz_0_xz_1, t_0_zzz_0_xzz_0, t_0_zzz_0_xzz_1, t_0_zzz_0_yy_1,\
                               t_0_zzz_0_yyy_0, t_0_zzz_0_yyy_1, t_0_zzz_0_yyz_0, t_0_zzz_0_yyz_1,\
                               t_0_zzz_0_yz_1, t_0_zzz_0_yzz_0, t_0_zzz_0_yzz_1, t_0_zzz_0_zz_1,\
                               t_0_zzz_0_zzz_0, t_0_zzz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyyz_0_xyy_0[i] = rpb_z[i] * t_0_yyy_0_xyy_0[i] + rwp_z[i] * t_0_yyy_0_xyy_1[i];

            t_0_yyyz_0_xxz_0[i] = rpb_z[i] * t_0_yyy_0_xxz_0[i] + rwp_z[i] * t_0_yyy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyz_0_xxy_0[i] = rpb_z[i] * t_0_yyy_0_xxy_0[i] + rwp_z[i] * t_0_yyy_0_xxy_1[i];

            t_0_yyyz_0_xxx_0[i] = rpb_z[i] * t_0_yyy_0_xxx_0[i] + rwp_z[i] * t_0_yyy_0_xxx_1[i];

            t_0_yyyy_0_zzz_0[i] = rpb_y[i] * t_0_yyy_0_zzz_0[i] + rwp_y[i] * t_0_yyy_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zzz_1[i];

            t_0_yyyy_0_yzz_0[i] = rpb_y[i] * t_0_yyy_0_yzz_0[i] + rwp_y[i] * t_0_yyy_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_yyyy_0_yyz_0[i] = rpb_y[i] * t_0_yyy_0_yyz_0[i] + rwp_y[i] * t_0_yyy_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyz_1[i] + fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_yyyy_0_yyy_0[i] = rpb_y[i] * t_0_yyy_0_yyy_0[i] + rwp_y[i] * t_0_yyy_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyy_0_xzz_0[i] = rpb_y[i] * t_0_yyy_0_xzz_0[i] + rwp_y[i] * t_0_yyy_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xzz_1[i];

            t_0_yyyy_0_xyz_0[i] = rpb_y[i] * t_0_yyy_0_xyz_0[i] + rwp_y[i] * t_0_yyy_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_yyyy_0_xyy_0[i] = rpb_y[i] * t_0_yyy_0_xyy_0[i] + rwp_y[i] * t_0_yyy_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xyy_1[i] + fze_0[i] * t_0_yyy_0_xy_1[i];

            t_0_yyyy_0_xxz_0[i] = rpb_y[i] * t_0_yyy_0_xxz_0[i] + rwp_y[i] * t_0_yyy_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxz_1[i];

            t_0_yyyy_0_xxy_0[i] = rpb_y[i] * t_0_yyy_0_xxy_0[i] + rwp_y[i] * t_0_yyy_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyy_0_xxx_0[i] = rpb_y[i] * t_0_yyy_0_xxx_0[i] + rwp_y[i] * t_0_yyy_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xxx_1[i];

            t_0_xzzz_0_zzz_0[i] = rpb_x[i] * t_0_zzz_0_zzz_0[i] + rwp_x[i] * t_0_zzz_0_zzz_1[i];

            t_0_xzzz_0_yzz_0[i] = rpb_x[i] * t_0_zzz_0_yzz_0[i] + rwp_x[i] * t_0_zzz_0_yzz_1[i];

            t_0_xzzz_0_yyz_0[i] = rpb_x[i] * t_0_zzz_0_yyz_0[i] + rwp_x[i] * t_0_zzz_0_yyz_1[i];

            t_0_xzzz_0_yyy_0[i] = rpb_x[i] * t_0_zzz_0_yyy_0[i] + rwp_x[i] * t_0_zzz_0_yyy_1[i];

            t_0_xzzz_0_xzz_0[i] = rpb_x[i] * t_0_zzz_0_xzz_0[i] + rwp_x[i] * t_0_zzz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_xyz_0[i] = rpb_x[i] * t_0_zzz_0_xyz_0[i] + rwp_x[i] * t_0_zzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yz_1[i];

            t_0_xzzz_0_xyy_0[i] = rpb_x[i] * t_0_zzz_0_xyy_0[i] + rwp_x[i] * t_0_zzz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_yy_1[i];

            t_0_xzzz_0_xxz_0[i] = rpb_x[i] * t_0_zzz_0_xxz_0[i] + rwp_x[i] * t_0_zzz_0_xxz_1[i] + fze_0[i] * t_0_zzz_0_xz_1[i];

            t_0_xzzz_0_xxy_0[i] = rpb_x[i] * t_0_zzz_0_xxy_0[i] + rwp_x[i] * t_0_zzz_0_xxy_1[i] + fze_0[i] * t_0_zzz_0_xy_1[i];

            t_0_xzzz_0_xxx_0[i] = rpb_x[i] * t_0_zzz_0_xxx_0[i] + rwp_x[i] * t_0_zzz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zzz_0_xx_1[i];

            t_0_xyzz_0_zzz_0[i] = rpb_x[i] * t_0_yzz_0_zzz_0[i] + rwp_x[i] * t_0_yzz_0_zzz_1[i];

            t_0_xyzz_0_yzz_0[i] = rpb_x[i] * t_0_yzz_0_yzz_0[i] + rwp_x[i] * t_0_yzz_0_yzz_1[i];

            t_0_xyzz_0_yyz_0[i] = rpb_x[i] * t_0_yzz_0_yyz_0[i] + rwp_x[i] * t_0_yzz_0_yyz_1[i];

            t_0_xyzz_0_yyy_0[i] = rpb_x[i] * t_0_yzz_0_yyy_0[i] + rwp_x[i] * t_0_yzz_0_yyy_1[i];

            t_0_xyzz_0_xzz_0[i] = rpb_y[i] * t_0_xzz_0_xzz_0[i] + rwp_y[i] * t_0_xzz_0_xzz_1[i];

            t_0_xyzz_0_xyz_0[i] = rpb_x[i] * t_0_yzz_0_xyz_0[i] + rwp_x[i] * t_0_yzz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_xyy_0[i] = rpb_x[i] * t_0_yzz_0_xyy_0[i] + rwp_x[i] * t_0_yzz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_yy_1[i];

            t_0_xyzz_0_xxz_0[i] = rpb_y[i] * t_0_xzz_0_xxz_0[i] + rwp_y[i] * t_0_xzz_0_xxz_1[i];

            t_0_xyzz_0_xxy_0[i] = rpb_x[i] * t_0_yzz_0_xxy_0[i] + rwp_x[i] * t_0_yzz_0_xxy_1[i] + fze_0[i] * t_0_yzz_0_xy_1[i];

            t_0_xyzz_0_xxx_0[i] = rpb_y[i] * t_0_xzz_0_xxx_0[i] + rwp_y[i] * t_0_xzz_0_xxx_1[i];

            t_0_xyyz_0_zzz_0[i] = rpb_x[i] * t_0_yyz_0_zzz_0[i] + rwp_x[i] * t_0_yyz_0_zzz_1[i];

            t_0_xyyz_0_yzz_0[i] = rpb_x[i] * t_0_yyz_0_yzz_0[i] + rwp_x[i] * t_0_yyz_0_yzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xyy_1, t_0_xxy_0_yyy_0, t_0_xxy_0_yyy_1, t_0_xxyz_0_xxz_0,\
                               t_0_xxyz_0_xyy_0, t_0_xxyz_0_xyz_0, t_0_xxyz_0_xzz_0,\
                               t_0_xxyz_0_yyy_0, t_0_xxyz_0_yyz_0, t_0_xxyz_0_yzz_0,\
                               t_0_xxyz_0_zzz_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxx_1, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxy_1, t_0_xxz_0_xxz_0, t_0_xxz_0_xxz_1, t_0_xxz_0_xyy_0,\
                               t_0_xxz_0_xyy_1, t_0_xxz_0_xyz_0, t_0_xxz_0_xyz_1, t_0_xxz_0_xz_1,\
                               t_0_xxz_0_xzz_0, t_0_xxz_0_xzz_1, t_0_xxz_0_yyz_0, t_0_xxz_0_yyz_1,\
                               t_0_xxz_0_yz_1, t_0_xxz_0_yzz_0, t_0_xxz_0_yzz_1, t_0_xxz_0_zz_1,\
                               t_0_xxz_0_zzz_0, t_0_xxz_0_zzz_1, t_0_xxzz_0_xxx_0, t_0_xxzz_0_xxy_0,\
                               t_0_xxzz_0_xxz_0, t_0_xxzz_0_xyy_0, t_0_xxzz_0_xyz_0,\
                               t_0_xxzz_0_xzz_0, t_0_xxzz_0_yyy_0, t_0_xxzz_0_yyz_0,\
                               t_0_xxzz_0_yzz_0, t_0_xxzz_0_zzz_0, t_0_xyy_0_xxx_0,\
                               t_0_xyy_0_xxx_1, t_0_xyy_0_xxy_0, t_0_xyy_0_xxy_1, t_0_xyy_0_xyy_0,\
                               t_0_xyy_0_xyy_1, t_0_xyyy_0_xxx_0, t_0_xyyy_0_xxy_0,\
                               t_0_xyyy_0_xxz_0, t_0_xyyy_0_xyy_0, t_0_xyyy_0_xyz_0,\
                               t_0_xyyy_0_xzz_0, t_0_xyyy_0_yyy_0, t_0_xyyy_0_yyz_0,\
                               t_0_xyyy_0_yzz_0, t_0_xyyy_0_zzz_0, t_0_xyyz_0_xxx_0,\
                               t_0_xyyz_0_xxy_0, t_0_xyyz_0_xxz_0, t_0_xyyz_0_xyy_0,\
                               t_0_xyyz_0_xyz_0, t_0_xyyz_0_xzz_0, t_0_xyyz_0_yyy_0,\
                               t_0_xyyz_0_yyz_0, t_0_xzz_0_xxz_0, t_0_xzz_0_xxz_1, t_0_xzz_0_xyz_0,\
                               t_0_xzz_0_xyz_1, t_0_xzz_0_xz_1, t_0_xzz_0_xzz_0, t_0_xzz_0_xzz_1,\
                               t_0_xzz_0_yyy_0, t_0_xzz_0_yyy_1, t_0_xzz_0_yyz_0, t_0_xzz_0_yyz_1,\
                               t_0_xzz_0_yz_1, t_0_xzz_0_yzz_0, t_0_xzz_0_yzz_1, t_0_xzz_0_zz_1,\
                               t_0_xzz_0_zzz_0, t_0_xzz_0_zzz_1, t_0_yyy_0_xx_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxx_1, t_0_yyy_0_xxy_0, t_0_yyy_0_xxy_1, t_0_yyy_0_xxz_0,\
                               t_0_yyy_0_xxz_1, t_0_yyy_0_xy_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyy_1,\
                               t_0_yyy_0_xyz_0, t_0_yyy_0_xyz_1, t_0_yyy_0_xz_1, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_xzz_1, t_0_yyy_0_yy_1, t_0_yyy_0_yyy_0, t_0_yyy_0_yyy_1,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yyz_1, t_0_yyy_0_yz_1, t_0_yyy_0_yzz_0,\
                               t_0_yyy_0_yzz_1, t_0_yyy_0_zz_1, t_0_yyy_0_zzz_0, t_0_yyy_0_zzz_1,\
                               t_0_yyz_0_xxz_0, t_0_yyz_0_xxz_1, t_0_yyz_0_xyz_0, t_0_yyz_0_xyz_1,\
                               t_0_yyz_0_xz_1, t_0_yyz_0_xzz_0, t_0_yyz_0_xzz_1, t_0_yyz_0_yyy_0,\
                               t_0_yyz_0_yyy_1, t_0_yyz_0_yyz_0, t_0_yyz_0_yyz_1, t_0_yyz_0_yz_1,\
                               t_0_yyz_0_zz_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xyz_0,\
                               t_0_zz_0_xyz_1, t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yzz_0,\
                               t_0_zz_0_yzz_1, t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyyz_0_yyz_0[i] = rpb_x[i] * t_0_yyz_0_yyz_0[i] + rwp_x[i] * t_0_yyz_0_yyz_1[i];

            t_0_xyyz_0_yyy_0[i] = rpb_x[i] * t_0_yyz_0_yyy_0[i] + rwp_x[i] * t_0_yyz_0_yyy_1[i];

            t_0_xyyz_0_xzz_0[i] = rpb_x[i] * t_0_yyz_0_xzz_0[i] + rwp_x[i] * t_0_yyz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_zz_1[i];

            t_0_xyyz_0_xyz_0[i] = rpb_x[i] * t_0_yyz_0_xyz_0[i] + rwp_x[i] * t_0_yyz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_xyy_0[i] = rpb_z[i] * t_0_xyy_0_xyy_0[i] + rwp_z[i] * t_0_xyy_0_xyy_1[i];

            t_0_xyyz_0_xxz_0[i] = rpb_x[i] * t_0_yyz_0_xxz_0[i] + rwp_x[i] * t_0_yyz_0_xxz_1[i] + fze_0[i] * t_0_yyz_0_xz_1[i];

            t_0_xyyz_0_xxy_0[i] = rpb_z[i] * t_0_xyy_0_xxy_0[i] + rwp_z[i] * t_0_xyy_0_xxy_1[i];

            t_0_xyyz_0_xxx_0[i] = rpb_z[i] * t_0_xyy_0_xxx_0[i] + rwp_z[i] * t_0_xyy_0_xxx_1[i];

            t_0_xyyy_0_zzz_0[i] = rpb_x[i] * t_0_yyy_0_zzz_0[i] + rwp_x[i] * t_0_yyy_0_zzz_1[i];

            t_0_xyyy_0_yzz_0[i] = rpb_x[i] * t_0_yyy_0_yzz_0[i] + rwp_x[i] * t_0_yyy_0_yzz_1[i];

            t_0_xyyy_0_yyz_0[i] = rpb_x[i] * t_0_yyy_0_yyz_0[i] + rwp_x[i] * t_0_yyy_0_yyz_1[i];

            t_0_xyyy_0_yyy_0[i] = rpb_x[i] * t_0_yyy_0_yyy_0[i] + rwp_x[i] * t_0_yyy_0_yyy_1[i];

            t_0_xyyy_0_xzz_0[i] = rpb_x[i] * t_0_yyy_0_xzz_0[i] + rwp_x[i] * t_0_yyy_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_zz_1[i];

            t_0_xyyy_0_xyz_0[i] = rpb_x[i] * t_0_yyy_0_xyz_0[i] + rwp_x[i] * t_0_yyy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yz_1[i];

            t_0_xyyy_0_xyy_0[i] = rpb_x[i] * t_0_yyy_0_xyy_0[i] + rwp_x[i] * t_0_yyy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xxz_0[i] = rpb_x[i] * t_0_yyy_0_xxz_0[i] + rwp_x[i] * t_0_yyy_0_xxz_1[i] + fze_0[i] * t_0_yyy_0_xz_1[i];

            t_0_xyyy_0_xxy_0[i] = rpb_x[i] * t_0_yyy_0_xxy_0[i] + rwp_x[i] * t_0_yyy_0_xxy_1[i] + fze_0[i] * t_0_yyy_0_xy_1[i];

            t_0_xyyy_0_xxx_0[i] = rpb_x[i] * t_0_yyy_0_xxx_0[i] + rwp_x[i] * t_0_yyy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yyy_0_xx_1[i];

            t_0_xxzz_0_zzz_0[i] = rpb_x[i] * t_0_xzz_0_zzz_0[i] + rwp_x[i] * t_0_xzz_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zzz_1[i];

            t_0_xxzz_0_yzz_0[i] = rpb_x[i] * t_0_xzz_0_yzz_0[i] + rwp_x[i] * t_0_xzz_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yzz_1[i];

            t_0_xxzz_0_yyz_0[i] = rpb_x[i] * t_0_xzz_0_yyz_0[i] + rwp_x[i] * t_0_xzz_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyz_1[i];

            t_0_xxzz_0_yyy_0[i] = rpb_x[i] * t_0_xzz_0_yyy_0[i] + rwp_x[i] * t_0_xzz_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yyy_1[i];

            t_0_xxzz_0_xzz_0[i] = rpb_x[i] * t_0_xzz_0_xzz_0[i] + rwp_x[i] * t_0_xzz_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_zz_1[i];

            t_0_xxzz_0_xyz_0[i] = rpb_x[i] * t_0_xzz_0_xyz_0[i] + rwp_x[i] * t_0_xzz_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_yz_1[i];

            t_0_xxzz_0_xyy_0[i] = rpb_z[i] * t_0_xxz_0_xyy_0[i] + rwp_z[i] * t_0_xxz_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xyy_1[i];

            t_0_xxzz_0_xxz_0[i] = rpb_x[i] * t_0_xzz_0_xxz_0[i] + rwp_x[i] * t_0_xzz_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_xzz_0_xz_1[i];

            t_0_xxzz_0_xxy_0[i] = rpb_z[i] * t_0_xxz_0_xxy_0[i] + rwp_z[i] * t_0_xxz_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxy_1[i];

            t_0_xxzz_0_xxx_0[i] = rpb_z[i] * t_0_xxz_0_xxx_0[i] + rwp_z[i] * t_0_xxz_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxx_1[i];

            t_0_xxyz_0_zzz_0[i] = rpb_y[i] * t_0_xxz_0_zzz_0[i] + rwp_y[i] * t_0_xxz_0_zzz_1[i];

            t_0_xxyz_0_yzz_0[i] = rpb_y[i] * t_0_xxz_0_yzz_0[i] + rwp_y[i] * t_0_xxz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_zz_1[i];

            t_0_xxyz_0_yyz_0[i] = rpb_y[i] * t_0_xxz_0_yyz_0[i] + rwp_y[i] * t_0_xxz_0_yyz_1[i] + fze_0[i] * t_0_xxz_0_yz_1[i];

            t_0_xxyz_0_yyy_0[i] = rpb_z[i] * t_0_xxy_0_yyy_0[i] + rwp_z[i] * t_0_xxy_0_yyy_1[i];

            t_0_xxyz_0_xzz_0[i] = rpb_y[i] * t_0_xxz_0_xzz_0[i] + rwp_y[i] * t_0_xxz_0_xzz_1[i];

            t_0_xxyz_0_xyz_0[i] = rpb_y[i] * t_0_xxz_0_xyz_0[i] + rwp_y[i] * t_0_xxz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xyy_0[i] = rpb_z[i] * t_0_xxy_0_xyy_0[i] + rwp_z[i] * t_0_xxy_0_xyy_1[i];

            t_0_xxyz_0_xxz_0[i] = rpb_y[i] * t_0_xxz_0_xxz_0[i] + rwp_y[i] * t_0_xxz_0_xxz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxz_0,\
                               t_0_xx_0_xxz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yyy_0,\
                               t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1, t_0_xx_0_yzz_0,\
                               t_0_xx_0_yzz_1, t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_xx_1,\
                               t_0_xxx_0_xxx_0, t_0_xxx_0_xxx_1, t_0_xxx_0_xxy_0, t_0_xxx_0_xxy_1,\
                               t_0_xxx_0_xxz_0, t_0_xxx_0_xxz_1, t_0_xxx_0_xy_1, t_0_xxx_0_xyy_0,\
                               t_0_xxx_0_xyy_1, t_0_xxx_0_xyz_0, t_0_xxx_0_xyz_1, t_0_xxx_0_xz_1,\
                               t_0_xxx_0_xzz_0, t_0_xxx_0_xzz_1, t_0_xxx_0_yy_1, t_0_xxx_0_yyy_0,\
                               t_0_xxx_0_yyy_1, t_0_xxx_0_yyz_0, t_0_xxx_0_yyz_1, t_0_xxx_0_yz_1,\
                               t_0_xxx_0_yzz_0, t_0_xxx_0_yzz_1, t_0_xxx_0_zz_1, t_0_xxx_0_zzz_0,\
                               t_0_xxx_0_zzz_1, t_0_xxxx_0_yyy_0, t_0_xxxx_0_yyz_0,\
                               t_0_xxxx_0_yzz_0, t_0_xxxx_0_zzz_0, t_0_xxxy_0_xxx_0,\
                               t_0_xxxy_0_xxy_0, t_0_xxxy_0_xxz_0, t_0_xxxy_0_xyy_0,\
                               t_0_xxxy_0_xyz_0, t_0_xxxy_0_xzz_0, t_0_xxxy_0_yyy_0,\
                               t_0_xxxy_0_yyz_0, t_0_xxxy_0_yzz_0, t_0_xxxy_0_zzz_0,\
                               t_0_xxxz_0_xxx_0, t_0_xxxz_0_xxy_0, t_0_xxxz_0_xxz_0,\
                               t_0_xxxz_0_xyy_0, t_0_xxxz_0_xyz_0, t_0_xxxz_0_xzz_0,\
                               t_0_xxxz_0_yyy_0, t_0_xxxz_0_yyz_0, t_0_xxxz_0_yzz_0,\
                               t_0_xxxz_0_zzz_0, t_0_xxy_0_xxx_0, t_0_xxy_0_xxx_1, t_0_xxy_0_xxy_0,\
                               t_0_xxy_0_xxy_1, t_0_xxy_0_xxz_0, t_0_xxy_0_xxz_1, t_0_xxy_0_xzz_0,\
                               t_0_xxy_0_xzz_1, t_0_xxyy_0_xxx_0, t_0_xxyy_0_xxy_0,\
                               t_0_xxyy_0_xxz_0, t_0_xxyy_0_xyy_0, t_0_xxyy_0_xyz_0,\
                               t_0_xxyy_0_xzz_0, t_0_xxyy_0_yyy_0, t_0_xxyy_0_yyz_0,\
                               t_0_xxyy_0_yzz_0, t_0_xxyy_0_zzz_0, t_0_xxyz_0_xxx_0,\
                               t_0_xxyz_0_xxy_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxx_1, t_0_xyy_0_xxy_0,\
                               t_0_xyy_0_xxy_1, t_0_xyy_0_xy_1, t_0_xyy_0_xyy_0, t_0_xyy_0_xyy_1,\
                               t_0_xyy_0_xyz_0, t_0_xyy_0_xyz_1, t_0_xyy_0_yy_1, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyy_1, t_0_xyy_0_yyz_0, t_0_xyy_0_yyz_1, t_0_xyy_0_yz_1,\
                               t_0_xyy_0_yzz_0, t_0_xyy_0_yzz_1, t_0_xyy_0_zzz_0, t_0_xyy_0_zzz_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1,\
                               t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zzz_0, t_0_yy_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxyz_0_xxy_0[i] = rpb_z[i] * t_0_xxy_0_xxy_0[i] + rwp_z[i] * t_0_xxy_0_xxy_1[i];

            t_0_xxyz_0_xxx_0[i] = rpb_y[i] * t_0_xxz_0_xxx_0[i] + rwp_y[i] * t_0_xxz_0_xxx_1[i];

            t_0_xxyy_0_zzz_0[i] = rpb_x[i] * t_0_xyy_0_zzz_0[i] + rwp_x[i] * t_0_xyy_0_zzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zzz_1[i];

            t_0_xxyy_0_yzz_0[i] = rpb_x[i] * t_0_xyy_0_yzz_0[i] + rwp_x[i] * t_0_xyy_0_yzz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yzz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yzz_1[i];

            t_0_xxyy_0_yyz_0[i] = rpb_x[i] * t_0_xyy_0_yyz_0[i] + rwp_x[i] * t_0_xyy_0_yyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyz_1[i];

            t_0_xxyy_0_yyy_0[i] = rpb_x[i] * t_0_xyy_0_yyy_0[i] + rwp_x[i] * t_0_xyy_0_yyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yyy_1[i];

            t_0_xxyy_0_xzz_0[i] = rpb_y[i] * t_0_xxy_0_xzz_0[i] + rwp_y[i] * t_0_xxy_0_xzz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xzz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xzz_1[i];

            t_0_xxyy_0_xyz_0[i] = rpb_x[i] * t_0_xyy_0_xyz_0[i] + rwp_x[i] * t_0_xyy_0_xyz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yz_1[i];

            t_0_xxyy_0_xyy_0[i] = rpb_x[i] * t_0_xyy_0_xyy_0[i] + rwp_x[i] * t_0_xyy_0_xyy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xyy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_yy_1[i];

            t_0_xxyy_0_xxz_0[i] = rpb_y[i] * t_0_xxy_0_xxz_0[i] + rwp_y[i] * t_0_xxy_0_xxz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxz_1[i];

            t_0_xxyy_0_xxy_0[i] = rpb_x[i] * t_0_xyy_0_xxy_0[i] + rwp_x[i] * t_0_xyy_0_xxy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xxy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_xyy_0_xy_1[i];

            t_0_xxyy_0_xxx_0[i] = rpb_y[i] * t_0_xxy_0_xxx_0[i] + rwp_y[i] * t_0_xxy_0_xxx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xxx_1[i];

            t_0_xxxz_0_zzz_0[i] = rpb_z[i] * t_0_xxx_0_zzz_0[i] + rwp_z[i] * t_0_xxx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxz_0_yzz_0[i] = rpb_z[i] * t_0_xxx_0_yzz_0[i] + rwp_z[i] * t_0_xxx_0_yzz_1[i] + fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxz_0_yyz_0[i] = rpb_z[i] * t_0_xxx_0_yyz_0[i] + rwp_z[i] * t_0_xxx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxz_0_yyy_0[i] = rpb_z[i] * t_0_xxx_0_yyy_0[i] + rwp_z[i] * t_0_xxx_0_yyy_1[i];

            t_0_xxxz_0_xzz_0[i] = rpb_z[i] * t_0_xxx_0_xzz_0[i] + rwp_z[i] * t_0_xxx_0_xzz_1[i] + fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxz_0_xyz_0[i] = rpb_z[i] * t_0_xxx_0_xyz_0[i] + rwp_z[i] * t_0_xxx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxz_0_xyy_0[i] = rpb_z[i] * t_0_xxx_0_xyy_0[i] + rwp_z[i] * t_0_xxx_0_xyy_1[i];

            t_0_xxxz_0_xxz_0[i] = rpb_z[i] * t_0_xxx_0_xxz_0[i] + rwp_z[i] * t_0_xxx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxz_0_xxy_0[i] = rpb_z[i] * t_0_xxx_0_xxy_0[i] + rwp_z[i] * t_0_xxx_0_xxy_1[i];

            t_0_xxxz_0_xxx_0[i] = rpb_z[i] * t_0_xxx_0_xxx_0[i] + rwp_z[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxy_0_zzz_0[i] = rpb_y[i] * t_0_xxx_0_zzz_0[i] + rwp_y[i] * t_0_xxx_0_zzz_1[i];

            t_0_xxxy_0_yzz_0[i] = rpb_y[i] * t_0_xxx_0_yzz_0[i] + rwp_y[i] * t_0_xxx_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxy_0_yyz_0[i] = rpb_y[i] * t_0_xxx_0_yyz_0[i] + rwp_y[i] * t_0_xxx_0_yyz_1[i] + fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxy_0_yyy_0[i] = rpb_y[i] * t_0_xxx_0_yyy_0[i] + rwp_y[i] * t_0_xxx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxy_0_xzz_0[i] = rpb_y[i] * t_0_xxx_0_xzz_0[i] + rwp_y[i] * t_0_xxx_0_xzz_1[i];

            t_0_xxxy_0_xyz_0[i] = rpb_y[i] * t_0_xxx_0_xyz_0[i] + rwp_y[i] * t_0_xxx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxy_0_xyy_0[i] = rpb_y[i] * t_0_xxx_0_xyy_0[i] + rwp_y[i] * t_0_xxx_0_xyy_1[i] + fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxy_0_xxz_0[i] = rpb_y[i] * t_0_xxx_0_xxz_0[i] + rwp_y[i] * t_0_xxx_0_xxz_1[i];

            t_0_xxxy_0_xxy_0[i] = rpb_y[i] * t_0_xxx_0_xxy_0[i] + rwp_y[i] * t_0_xxx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_xxx_0[i] = rpb_y[i] * t_0_xxx_0_xxx_0[i] + rwp_y[i] * t_0_xxx_0_xxx_1[i];

            t_0_xxxx_0_zzz_0[i] = rpb_x[i] * t_0_xxx_0_zzz_0[i] + rwp_x[i] * t_0_xxx_0_zzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zzz_1[i];

            t_0_xxxx_0_yzz_0[i] = rpb_x[i] * t_0_xxx_0_yzz_0[i] + rwp_x[i] * t_0_xxx_0_yzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yzz_1[i];

            t_0_xxxx_0_yyz_0[i] = rpb_x[i] * t_0_xxx_0_yyz_0[i] + rwp_x[i] * t_0_xxx_0_yyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyz_1[i];

            t_0_xxxx_0_yyy_0[i] = rpb_x[i] * t_0_xxx_0_yyy_0[i] + rwp_x[i] * t_0_xxx_0_yyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yyy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1,\
                               t_0_xx_0_xxy_0, t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xxx_0_xx_1, t_0_xxx_0_xxx_0,\
                               t_0_xxx_0_xxx_1, t_0_xxx_0_xxy_0, t_0_xxx_0_xxy_1, t_0_xxx_0_xxz_0,\
                               t_0_xxx_0_xxz_1, t_0_xxx_0_xy_1, t_0_xxx_0_xyy_0, t_0_xxx_0_xyy_1,\
                               t_0_xxx_0_xyz_0, t_0_xxx_0_xyz_1, t_0_xxx_0_xz_1, t_0_xxx_0_xzz_0,\
                               t_0_xxx_0_xzz_1, t_0_xxx_0_yy_1, t_0_xxx_0_yz_1, t_0_xxx_0_zz_1,\
                               t_0_xxxx_0_xxx_0, t_0_xxxx_0_xxy_0, t_0_xxxx_0_xxz_0,\
                               t_0_xxxx_0_xyy_0, t_0_xxxx_0_xyz_0, t_0_xxxx_0_xzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxx_0_xzz_0[i] = rpb_x[i] * t_0_xxx_0_xzz_0[i] + rwp_x[i] * t_0_xxx_0_xzz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xzz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxx_0_xyz_0[i] = rpb_x[i] * t_0_xxx_0_xyz_0[i] + rwp_x[i] * t_0_xxx_0_xyz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yz_1[i];

            t_0_xxxx_0_xyy_0[i] = rpb_x[i] * t_0_xxx_0_xyy_0[i] + rwp_x[i] * t_0_xxx_0_xyy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xyy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxx_0_xxz_0[i] = rpb_x[i] * t_0_xxx_0_xxz_0[i] + rwp_x[i] * t_0_xxx_0_xxz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxz_1[i] + fze_0[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxx_0_xxy_0[i] = rpb_x[i] * t_0_xxx_0_xxy_0[i] + rwp_x[i] * t_0_xxx_0_xxy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxy_1[i] + fze_0[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxx_0_xxx_0[i] = rpb_x[i] * t_0_xxx_0_xxx_0[i] + rwp_x[i] * t_0_xxx_0_xxx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xxx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xxx_0_xx_1[i];
        }
    }
}


} // derirec namespace
