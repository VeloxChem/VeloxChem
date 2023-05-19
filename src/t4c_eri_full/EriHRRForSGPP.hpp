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
compHostHRRForSGPP_V0(      BufferHostXY<T>&      intsBufferSGPP,
                      const BufferHostX<int32_t>& intsIndexesSGPP,
                      const BufferHostXY<T>&      intsBufferSGSP,
                      const BufferHostX<int32_t>& intsIndexesSGSP,
                      const BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SGSP) integral components

    t_0_zzzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(0));

    t_0_zzzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(1));

    t_0_zzzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(2));

    t_0_yzzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(3));

    t_0_yzzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(4));

    t_0_yzzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(5));

    t_0_yyzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(6));

    t_0_yyzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(7));

    t_0_yyzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(8));

    t_0_yyyz_0_z = intsBufferSGSP.data(intsIndexesSGSP(9));

    t_0_yyyz_0_y = intsBufferSGSP.data(intsIndexesSGSP(10));

    t_0_yyyz_0_x = intsBufferSGSP.data(intsIndexesSGSP(11));

    t_0_yyyy_0_z = intsBufferSGSP.data(intsIndexesSGSP(12));

    t_0_yyyy_0_y = intsBufferSGSP.data(intsIndexesSGSP(13));

    t_0_yyyy_0_x = intsBufferSGSP.data(intsIndexesSGSP(14));

    t_0_xzzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(15));

    t_0_xzzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(16));

    t_0_xzzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(17));

    t_0_xyzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(18));

    t_0_xyzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(19));

    t_0_xyzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(20));

    t_0_xyyz_0_z = intsBufferSGSP.data(intsIndexesSGSP(21));

    t_0_xyyz_0_y = intsBufferSGSP.data(intsIndexesSGSP(22));

    t_0_xyyz_0_x = intsBufferSGSP.data(intsIndexesSGSP(23));

    t_0_xyyy_0_z = intsBufferSGSP.data(intsIndexesSGSP(24));

    t_0_xyyy_0_y = intsBufferSGSP.data(intsIndexesSGSP(25));

    t_0_xyyy_0_x = intsBufferSGSP.data(intsIndexesSGSP(26));

    t_0_xxzz_0_z = intsBufferSGSP.data(intsIndexesSGSP(27));

    t_0_xxzz_0_y = intsBufferSGSP.data(intsIndexesSGSP(28));

    t_0_xxzz_0_x = intsBufferSGSP.data(intsIndexesSGSP(29));

    t_0_xxyz_0_z = intsBufferSGSP.data(intsIndexesSGSP(30));

    t_0_xxyz_0_y = intsBufferSGSP.data(intsIndexesSGSP(31));

    t_0_xxyz_0_x = intsBufferSGSP.data(intsIndexesSGSP(32));

    t_0_xxyy_0_z = intsBufferSGSP.data(intsIndexesSGSP(33));

    t_0_xxyy_0_y = intsBufferSGSP.data(intsIndexesSGSP(34));

    t_0_xxyy_0_x = intsBufferSGSP.data(intsIndexesSGSP(35));

    t_0_xxxz_0_z = intsBufferSGSP.data(intsIndexesSGSP(36));

    t_0_xxxz_0_y = intsBufferSGSP.data(intsIndexesSGSP(37));

    t_0_xxxz_0_x = intsBufferSGSP.data(intsIndexesSGSP(38));

    t_0_xxxy_0_z = intsBufferSGSP.data(intsIndexesSGSP(39));

    t_0_xxxy_0_y = intsBufferSGSP.data(intsIndexesSGSP(40));

    t_0_xxxy_0_x = intsBufferSGSP.data(intsIndexesSGSP(41));

    t_0_xxxx_0_z = intsBufferSGSP.data(intsIndexesSGSP(42));

    t_0_xxxx_0_y = intsBufferSGSP.data(intsIndexesSGSP(43));

    t_0_xxxx_0_x = intsBufferSGSP.data(intsIndexesSGSP(44));

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

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_yyyz_0_x, t_0_yyyz_0_xx, t_0_yyyz_0_xy,\
                           t_0_yyyz_0_xz, t_0_yyyz_0_y, t_0_yyyz_0_yy, t_0_yyyz_0_yz,\
                           t_0_yyyz_0_z, t_0_yyyz_0_zz, t_0_yyyz_x_x, t_0_yyyz_x_y,\
                           t_0_yyyz_x_z, t_0_yyyz_y_x, t_0_yyyz_y_y, t_0_yyyz_y_z, t_0_yyyz_z_x,\
                           t_0_yyyz_z_y, t_0_yyyz_z_z, t_0_yyzz_0_x, t_0_yyzz_0_xx,\
                           t_0_yyzz_0_xy, t_0_yyzz_0_xz, t_0_yyzz_0_y, t_0_yyzz_0_yy,\
                           t_0_yyzz_0_yz, t_0_yyzz_0_z, t_0_yyzz_0_zz, t_0_yyzz_x_x,\
                           t_0_yyzz_x_y, t_0_yyzz_x_z, t_0_yyzz_y_x, t_0_yyzz_y_y, t_0_yyzz_y_z,\
                           t_0_yyzz_z_x, t_0_yyzz_z_y, t_0_yyzz_z_z, t_0_yzzz_0_x, t_0_yzzz_0_xx,\
                           t_0_yzzz_0_xy, t_0_yzzz_0_xz, t_0_yzzz_0_y, t_0_yzzz_0_yy,\
                           t_0_yzzz_0_yz, t_0_yzzz_0_z, t_0_yzzz_0_zz, t_0_yzzz_x_x,\
                           t_0_yzzz_x_y, t_0_yzzz_x_z, t_0_yzzz_y_x, t_0_yzzz_y_y, t_0_yzzz_y_z,\
                           t_0_yzzz_z_x, t_0_yzzz_z_y, t_0_yzzz_z_z, t_0_zzzz_0_x, t_0_zzzz_0_xx,\
                           t_0_zzzz_0_xy, t_0_zzzz_0_xz, t_0_zzzz_0_y, t_0_zzzz_0_yy,\
                           t_0_zzzz_0_yz, t_0_zzzz_0_z, t_0_zzzz_0_zz, t_0_zzzz_x_x,\
                           t_0_zzzz_x_y, t_0_zzzz_x_z, t_0_zzzz_y_x, t_0_zzzz_y_y, t_0_zzzz_y_z,\
                           t_0_zzzz_z_x, t_0_zzzz_z_y, t_0_zzzz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_z_z[i] = t_0_zzzz_0_zz[i] - rcd_z[i] * t_0_zzzz_0_z[i];

        t_0_zzzz_z_y[i] = t_0_zzzz_0_yz[i] - rcd_z[i] * t_0_zzzz_0_y[i];

        t_0_zzzz_z_x[i] = t_0_zzzz_0_xz[i] - rcd_z[i] * t_0_zzzz_0_x[i];

        t_0_zzzz_y_z[i] = t_0_zzzz_0_yz[i] - rcd_y[i] * t_0_zzzz_0_z[i];

        t_0_zzzz_y_y[i] = t_0_zzzz_0_yy[i] - rcd_y[i] * t_0_zzzz_0_y[i];

        t_0_zzzz_y_x[i] = t_0_zzzz_0_xy[i] - rcd_y[i] * t_0_zzzz_0_x[i];

        t_0_zzzz_x_z[i] = t_0_zzzz_0_xz[i] - rcd_x[i] * t_0_zzzz_0_z[i];

        t_0_zzzz_x_y[i] = t_0_zzzz_0_xy[i] - rcd_x[i] * t_0_zzzz_0_y[i];

        t_0_zzzz_x_x[i] = t_0_zzzz_0_xx[i] - rcd_x[i] * t_0_zzzz_0_x[i];

        t_0_yzzz_z_z[i] = t_0_yzzz_0_zz[i] - rcd_z[i] * t_0_yzzz_0_z[i];

        t_0_yzzz_z_y[i] = t_0_yzzz_0_yz[i] - rcd_z[i] * t_0_yzzz_0_y[i];

        t_0_yzzz_z_x[i] = t_0_yzzz_0_xz[i] - rcd_z[i] * t_0_yzzz_0_x[i];

        t_0_yzzz_y_z[i] = t_0_yzzz_0_yz[i] - rcd_y[i] * t_0_yzzz_0_z[i];

        t_0_yzzz_y_y[i] = t_0_yzzz_0_yy[i] - rcd_y[i] * t_0_yzzz_0_y[i];

        t_0_yzzz_y_x[i] = t_0_yzzz_0_xy[i] - rcd_y[i] * t_0_yzzz_0_x[i];

        t_0_yzzz_x_z[i] = t_0_yzzz_0_xz[i] - rcd_x[i] * t_0_yzzz_0_z[i];

        t_0_yzzz_x_y[i] = t_0_yzzz_0_xy[i] - rcd_x[i] * t_0_yzzz_0_y[i];

        t_0_yzzz_x_x[i] = t_0_yzzz_0_xx[i] - rcd_x[i] * t_0_yzzz_0_x[i];

        t_0_yyzz_z_z[i] = t_0_yyzz_0_zz[i] - rcd_z[i] * t_0_yyzz_0_z[i];

        t_0_yyzz_z_y[i] = t_0_yyzz_0_yz[i] - rcd_z[i] * t_0_yyzz_0_y[i];

        t_0_yyzz_z_x[i] = t_0_yyzz_0_xz[i] - rcd_z[i] * t_0_yyzz_0_x[i];

        t_0_yyzz_y_z[i] = t_0_yyzz_0_yz[i] - rcd_y[i] * t_0_yyzz_0_z[i];

        t_0_yyzz_y_y[i] = t_0_yyzz_0_yy[i] - rcd_y[i] * t_0_yyzz_0_y[i];

        t_0_yyzz_y_x[i] = t_0_yyzz_0_xy[i] - rcd_y[i] * t_0_yyzz_0_x[i];

        t_0_yyzz_x_z[i] = t_0_yyzz_0_xz[i] - rcd_x[i] * t_0_yyzz_0_z[i];

        t_0_yyzz_x_y[i] = t_0_yyzz_0_xy[i] - rcd_x[i] * t_0_yyzz_0_y[i];

        t_0_yyzz_x_x[i] = t_0_yyzz_0_xx[i] - rcd_x[i] * t_0_yyzz_0_x[i];

        t_0_yyyz_z_z[i] = t_0_yyyz_0_zz[i] - rcd_z[i] * t_0_yyyz_0_z[i];

        t_0_yyyz_z_y[i] = t_0_yyyz_0_yz[i] - rcd_z[i] * t_0_yyyz_0_y[i];

        t_0_yyyz_z_x[i] = t_0_yyyz_0_xz[i] - rcd_z[i] * t_0_yyyz_0_x[i];

        t_0_yyyz_y_z[i] = t_0_yyyz_0_yz[i] - rcd_y[i] * t_0_yyyz_0_z[i];

        t_0_yyyz_y_y[i] = t_0_yyyz_0_yy[i] - rcd_y[i] * t_0_yyyz_0_y[i];

        t_0_yyyz_y_x[i] = t_0_yyyz_0_xy[i] - rcd_y[i] * t_0_yyyz_0_x[i];

        t_0_yyyz_x_z[i] = t_0_yyyz_0_xz[i] - rcd_x[i] * t_0_yyyz_0_z[i];

        t_0_yyyz_x_y[i] = t_0_yyyz_0_xy[i] - rcd_x[i] * t_0_yyyz_0_y[i];

        t_0_yyyz_x_x[i] = t_0_yyyz_0_xx[i] - rcd_x[i] * t_0_yyyz_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xyyz_0_x, t_0_xyyz_0_xx, t_0_xyyz_0_xy,\
                           t_0_xyyz_0_xz, t_0_xyyz_0_y, t_0_xyyz_0_yy, t_0_xyyz_0_yz,\
                           t_0_xyyz_0_z, t_0_xyyz_0_zz, t_0_xyyz_x_x, t_0_xyyz_x_y,\
                           t_0_xyyz_x_z, t_0_xyyz_y_x, t_0_xyyz_y_y, t_0_xyyz_y_z, t_0_xyyz_z_x,\
                           t_0_xyyz_z_y, t_0_xyyz_z_z, t_0_xyzz_0_x, t_0_xyzz_0_xx,\
                           t_0_xyzz_0_xy, t_0_xyzz_0_xz, t_0_xyzz_0_y, t_0_xyzz_0_yy,\
                           t_0_xyzz_0_yz, t_0_xyzz_0_z, t_0_xyzz_0_zz, t_0_xyzz_x_x,\
                           t_0_xyzz_x_y, t_0_xyzz_x_z, t_0_xyzz_y_x, t_0_xyzz_y_y, t_0_xyzz_y_z,\
                           t_0_xyzz_z_x, t_0_xyzz_z_y, t_0_xyzz_z_z, t_0_xzzz_0_x, t_0_xzzz_0_xx,\
                           t_0_xzzz_0_xy, t_0_xzzz_0_xz, t_0_xzzz_0_y, t_0_xzzz_0_yy,\
                           t_0_xzzz_0_yz, t_0_xzzz_0_z, t_0_xzzz_0_zz, t_0_xzzz_x_x,\
                           t_0_xzzz_x_y, t_0_xzzz_x_z, t_0_xzzz_y_x, t_0_xzzz_y_y, t_0_xzzz_y_z,\
                           t_0_xzzz_z_x, t_0_xzzz_z_y, t_0_xzzz_z_z, t_0_yyyy_0_x, t_0_yyyy_0_xx,\
                           t_0_yyyy_0_xy, t_0_yyyy_0_xz, t_0_yyyy_0_y, t_0_yyyy_0_yy,\
                           t_0_yyyy_0_yz, t_0_yyyy_0_z, t_0_yyyy_0_zz, t_0_yyyy_x_x,\
                           t_0_yyyy_x_y, t_0_yyyy_x_z, t_0_yyyy_y_x, t_0_yyyy_y_y, t_0_yyyy_y_z,\
                           t_0_yyyy_z_x, t_0_yyyy_z_y, t_0_yyyy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_yyyy_z_z[i] = t_0_yyyy_0_zz[i] - rcd_z[i] * t_0_yyyy_0_z[i];

        t_0_yyyy_z_y[i] = t_0_yyyy_0_yz[i] - rcd_z[i] * t_0_yyyy_0_y[i];

        t_0_yyyy_z_x[i] = t_0_yyyy_0_xz[i] - rcd_z[i] * t_0_yyyy_0_x[i];

        t_0_yyyy_y_z[i] = t_0_yyyy_0_yz[i] - rcd_y[i] * t_0_yyyy_0_z[i];

        t_0_yyyy_y_y[i] = t_0_yyyy_0_yy[i] - rcd_y[i] * t_0_yyyy_0_y[i];

        t_0_yyyy_y_x[i] = t_0_yyyy_0_xy[i] - rcd_y[i] * t_0_yyyy_0_x[i];

        t_0_yyyy_x_z[i] = t_0_yyyy_0_xz[i] - rcd_x[i] * t_0_yyyy_0_z[i];

        t_0_yyyy_x_y[i] = t_0_yyyy_0_xy[i] - rcd_x[i] * t_0_yyyy_0_y[i];

        t_0_yyyy_x_x[i] = t_0_yyyy_0_xx[i] - rcd_x[i] * t_0_yyyy_0_x[i];

        t_0_xzzz_z_z[i] = t_0_xzzz_0_zz[i] - rcd_z[i] * t_0_xzzz_0_z[i];

        t_0_xzzz_z_y[i] = t_0_xzzz_0_yz[i] - rcd_z[i] * t_0_xzzz_0_y[i];

        t_0_xzzz_z_x[i] = t_0_xzzz_0_xz[i] - rcd_z[i] * t_0_xzzz_0_x[i];

        t_0_xzzz_y_z[i] = t_0_xzzz_0_yz[i] - rcd_y[i] * t_0_xzzz_0_z[i];

        t_0_xzzz_y_y[i] = t_0_xzzz_0_yy[i] - rcd_y[i] * t_0_xzzz_0_y[i];

        t_0_xzzz_y_x[i] = t_0_xzzz_0_xy[i] - rcd_y[i] * t_0_xzzz_0_x[i];

        t_0_xzzz_x_z[i] = t_0_xzzz_0_xz[i] - rcd_x[i] * t_0_xzzz_0_z[i];

        t_0_xzzz_x_y[i] = t_0_xzzz_0_xy[i] - rcd_x[i] * t_0_xzzz_0_y[i];

        t_0_xzzz_x_x[i] = t_0_xzzz_0_xx[i] - rcd_x[i] * t_0_xzzz_0_x[i];

        t_0_xyzz_z_z[i] = t_0_xyzz_0_zz[i] - rcd_z[i] * t_0_xyzz_0_z[i];

        t_0_xyzz_z_y[i] = t_0_xyzz_0_yz[i] - rcd_z[i] * t_0_xyzz_0_y[i];

        t_0_xyzz_z_x[i] = t_0_xyzz_0_xz[i] - rcd_z[i] * t_0_xyzz_0_x[i];

        t_0_xyzz_y_z[i] = t_0_xyzz_0_yz[i] - rcd_y[i] * t_0_xyzz_0_z[i];

        t_0_xyzz_y_y[i] = t_0_xyzz_0_yy[i] - rcd_y[i] * t_0_xyzz_0_y[i];

        t_0_xyzz_y_x[i] = t_0_xyzz_0_xy[i] - rcd_y[i] * t_0_xyzz_0_x[i];

        t_0_xyzz_x_z[i] = t_0_xyzz_0_xz[i] - rcd_x[i] * t_0_xyzz_0_z[i];

        t_0_xyzz_x_y[i] = t_0_xyzz_0_xy[i] - rcd_x[i] * t_0_xyzz_0_y[i];

        t_0_xyzz_x_x[i] = t_0_xyzz_0_xx[i] - rcd_x[i] * t_0_xyzz_0_x[i];

        t_0_xyyz_z_z[i] = t_0_xyyz_0_zz[i] - rcd_z[i] * t_0_xyyz_0_z[i];

        t_0_xyyz_z_y[i] = t_0_xyyz_0_yz[i] - rcd_z[i] * t_0_xyyz_0_y[i];

        t_0_xyyz_z_x[i] = t_0_xyyz_0_xz[i] - rcd_z[i] * t_0_xyyz_0_x[i];

        t_0_xyyz_y_z[i] = t_0_xyyz_0_yz[i] - rcd_y[i] * t_0_xyyz_0_z[i];

        t_0_xyyz_y_y[i] = t_0_xyyz_0_yy[i] - rcd_y[i] * t_0_xyyz_0_y[i];

        t_0_xyyz_y_x[i] = t_0_xyyz_0_xy[i] - rcd_y[i] * t_0_xyyz_0_x[i];

        t_0_xyyz_x_z[i] = t_0_xyyz_0_xz[i] - rcd_x[i] * t_0_xyyz_0_z[i];

        t_0_xyyz_x_y[i] = t_0_xyyz_0_xy[i] - rcd_x[i] * t_0_xyyz_0_y[i];

        t_0_xyyz_x_x[i] = t_0_xyyz_0_xx[i] - rcd_x[i] * t_0_xyyz_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxyy_0_x, t_0_xxyy_0_xx, t_0_xxyy_0_xy,\
                           t_0_xxyy_0_xz, t_0_xxyy_0_y, t_0_xxyy_0_yy, t_0_xxyy_0_yz,\
                           t_0_xxyy_0_z, t_0_xxyy_0_zz, t_0_xxyy_x_x, t_0_xxyy_x_y,\
                           t_0_xxyy_x_z, t_0_xxyy_y_x, t_0_xxyy_y_y, t_0_xxyy_y_z, t_0_xxyy_z_x,\
                           t_0_xxyy_z_y, t_0_xxyy_z_z, t_0_xxyz_0_x, t_0_xxyz_0_xx,\
                           t_0_xxyz_0_xy, t_0_xxyz_0_xz, t_0_xxyz_0_y, t_0_xxyz_0_yy,\
                           t_0_xxyz_0_yz, t_0_xxyz_0_z, t_0_xxyz_0_zz, t_0_xxyz_x_x,\
                           t_0_xxyz_x_y, t_0_xxyz_x_z, t_0_xxyz_y_x, t_0_xxyz_y_y, t_0_xxyz_y_z,\
                           t_0_xxyz_z_x, t_0_xxyz_z_y, t_0_xxyz_z_z, t_0_xxzz_0_x, t_0_xxzz_0_xx,\
                           t_0_xxzz_0_xy, t_0_xxzz_0_xz, t_0_xxzz_0_y, t_0_xxzz_0_yy,\
                           t_0_xxzz_0_yz, t_0_xxzz_0_z, t_0_xxzz_0_zz, t_0_xxzz_x_x,\
                           t_0_xxzz_x_y, t_0_xxzz_x_z, t_0_xxzz_y_x, t_0_xxzz_y_y, t_0_xxzz_y_z,\
                           t_0_xxzz_z_x, t_0_xxzz_z_y, t_0_xxzz_z_z, t_0_xyyy_0_x, t_0_xyyy_0_xx,\
                           t_0_xyyy_0_xy, t_0_xyyy_0_xz, t_0_xyyy_0_y, t_0_xyyy_0_yy,\
                           t_0_xyyy_0_yz, t_0_xyyy_0_z, t_0_xyyy_0_zz, t_0_xyyy_x_x,\
                           t_0_xyyy_x_y, t_0_xyyy_x_z, t_0_xyyy_y_x, t_0_xyyy_y_y, t_0_xyyy_y_z,\
                           t_0_xyyy_z_x, t_0_xyyy_z_y, t_0_xyyy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xyyy_z_z[i] = t_0_xyyy_0_zz[i] - rcd_z[i] * t_0_xyyy_0_z[i];

        t_0_xyyy_z_y[i] = t_0_xyyy_0_yz[i] - rcd_z[i] * t_0_xyyy_0_y[i];

        t_0_xyyy_z_x[i] = t_0_xyyy_0_xz[i] - rcd_z[i] * t_0_xyyy_0_x[i];

        t_0_xyyy_y_z[i] = t_0_xyyy_0_yz[i] - rcd_y[i] * t_0_xyyy_0_z[i];

        t_0_xyyy_y_y[i] = t_0_xyyy_0_yy[i] - rcd_y[i] * t_0_xyyy_0_y[i];

        t_0_xyyy_y_x[i] = t_0_xyyy_0_xy[i] - rcd_y[i] * t_0_xyyy_0_x[i];

        t_0_xyyy_x_z[i] = t_0_xyyy_0_xz[i] - rcd_x[i] * t_0_xyyy_0_z[i];

        t_0_xyyy_x_y[i] = t_0_xyyy_0_xy[i] - rcd_x[i] * t_0_xyyy_0_y[i];

        t_0_xyyy_x_x[i] = t_0_xyyy_0_xx[i] - rcd_x[i] * t_0_xyyy_0_x[i];

        t_0_xxzz_z_z[i] = t_0_xxzz_0_zz[i] - rcd_z[i] * t_0_xxzz_0_z[i];

        t_0_xxzz_z_y[i] = t_0_xxzz_0_yz[i] - rcd_z[i] * t_0_xxzz_0_y[i];

        t_0_xxzz_z_x[i] = t_0_xxzz_0_xz[i] - rcd_z[i] * t_0_xxzz_0_x[i];

        t_0_xxzz_y_z[i] = t_0_xxzz_0_yz[i] - rcd_y[i] * t_0_xxzz_0_z[i];

        t_0_xxzz_y_y[i] = t_0_xxzz_0_yy[i] - rcd_y[i] * t_0_xxzz_0_y[i];

        t_0_xxzz_y_x[i] = t_0_xxzz_0_xy[i] - rcd_y[i] * t_0_xxzz_0_x[i];

        t_0_xxzz_x_z[i] = t_0_xxzz_0_xz[i] - rcd_x[i] * t_0_xxzz_0_z[i];

        t_0_xxzz_x_y[i] = t_0_xxzz_0_xy[i] - rcd_x[i] * t_0_xxzz_0_y[i];

        t_0_xxzz_x_x[i] = t_0_xxzz_0_xx[i] - rcd_x[i] * t_0_xxzz_0_x[i];

        t_0_xxyz_z_z[i] = t_0_xxyz_0_zz[i] - rcd_z[i] * t_0_xxyz_0_z[i];

        t_0_xxyz_z_y[i] = t_0_xxyz_0_yz[i] - rcd_z[i] * t_0_xxyz_0_y[i];

        t_0_xxyz_z_x[i] = t_0_xxyz_0_xz[i] - rcd_z[i] * t_0_xxyz_0_x[i];

        t_0_xxyz_y_z[i] = t_0_xxyz_0_yz[i] - rcd_y[i] * t_0_xxyz_0_z[i];

        t_0_xxyz_y_y[i] = t_0_xxyz_0_yy[i] - rcd_y[i] * t_0_xxyz_0_y[i];

        t_0_xxyz_y_x[i] = t_0_xxyz_0_xy[i] - rcd_y[i] * t_0_xxyz_0_x[i];

        t_0_xxyz_x_z[i] = t_0_xxyz_0_xz[i] - rcd_x[i] * t_0_xxyz_0_z[i];

        t_0_xxyz_x_y[i] = t_0_xxyz_0_xy[i] - rcd_x[i] * t_0_xxyz_0_y[i];

        t_0_xxyz_x_x[i] = t_0_xxyz_0_xx[i] - rcd_x[i] * t_0_xxyz_0_x[i];

        t_0_xxyy_z_z[i] = t_0_xxyy_0_zz[i] - rcd_z[i] * t_0_xxyy_0_z[i];

        t_0_xxyy_z_y[i] = t_0_xxyy_0_yz[i] - rcd_z[i] * t_0_xxyy_0_y[i];

        t_0_xxyy_z_x[i] = t_0_xxyy_0_xz[i] - rcd_z[i] * t_0_xxyy_0_x[i];

        t_0_xxyy_y_z[i] = t_0_xxyy_0_yz[i] - rcd_y[i] * t_0_xxyy_0_z[i];

        t_0_xxyy_y_y[i] = t_0_xxyy_0_yy[i] - rcd_y[i] * t_0_xxyy_0_y[i];

        t_0_xxyy_y_x[i] = t_0_xxyy_0_xy[i] - rcd_y[i] * t_0_xxyy_0_x[i];

        t_0_xxyy_x_z[i] = t_0_xxyy_0_xz[i] - rcd_x[i] * t_0_xxyy_0_z[i];

        t_0_xxyy_x_y[i] = t_0_xxyy_0_xy[i] - rcd_x[i] * t_0_xxyy_0_y[i];

        t_0_xxyy_x_x[i] = t_0_xxyy_0_xx[i] - rcd_x[i] * t_0_xxyy_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_0_x, t_0_xxxx_0_xx, t_0_xxxx_0_xy,\
                           t_0_xxxx_0_xz, t_0_xxxx_0_y, t_0_xxxx_0_yy, t_0_xxxx_0_yz,\
                           t_0_xxxx_0_z, t_0_xxxx_0_zz, t_0_xxxx_x_x, t_0_xxxx_x_y,\
                           t_0_xxxx_x_z, t_0_xxxx_y_x, t_0_xxxx_y_y, t_0_xxxx_y_z, t_0_xxxx_z_x,\
                           t_0_xxxx_z_y, t_0_xxxx_z_z, t_0_xxxy_0_x, t_0_xxxy_0_xx,\
                           t_0_xxxy_0_xy, t_0_xxxy_0_xz, t_0_xxxy_0_y, t_0_xxxy_0_yy,\
                           t_0_xxxy_0_yz, t_0_xxxy_0_z, t_0_xxxy_0_zz, t_0_xxxy_x_x,\
                           t_0_xxxy_x_y, t_0_xxxy_x_z, t_0_xxxy_y_x, t_0_xxxy_y_y, t_0_xxxy_y_z,\
                           t_0_xxxy_z_x, t_0_xxxy_z_y, t_0_xxxy_z_z, t_0_xxxz_0_x, t_0_xxxz_0_xx,\
                           t_0_xxxz_0_xy, t_0_xxxz_0_xz, t_0_xxxz_0_y, t_0_xxxz_0_yy,\
                           t_0_xxxz_0_yz, t_0_xxxz_0_z, t_0_xxxz_0_zz, t_0_xxxz_x_x,\
                           t_0_xxxz_x_y, t_0_xxxz_x_z, t_0_xxxz_y_x, t_0_xxxz_y_y, t_0_xxxz_y_z,\
                           t_0_xxxz_z_x, t_0_xxxz_z_y, t_0_xxxz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xxxz_z_z[i] = t_0_xxxz_0_zz[i] - rcd_z[i] * t_0_xxxz_0_z[i];

        t_0_xxxz_z_y[i] = t_0_xxxz_0_yz[i] - rcd_z[i] * t_0_xxxz_0_y[i];

        t_0_xxxz_z_x[i] = t_0_xxxz_0_xz[i] - rcd_z[i] * t_0_xxxz_0_x[i];

        t_0_xxxz_y_z[i] = t_0_xxxz_0_yz[i] - rcd_y[i] * t_0_xxxz_0_z[i];

        t_0_xxxz_y_y[i] = t_0_xxxz_0_yy[i] - rcd_y[i] * t_0_xxxz_0_y[i];

        t_0_xxxz_y_x[i] = t_0_xxxz_0_xy[i] - rcd_y[i] * t_0_xxxz_0_x[i];

        t_0_xxxz_x_z[i] = t_0_xxxz_0_xz[i] - rcd_x[i] * t_0_xxxz_0_z[i];

        t_0_xxxz_x_y[i] = t_0_xxxz_0_xy[i] - rcd_x[i] * t_0_xxxz_0_y[i];

        t_0_xxxz_x_x[i] = t_0_xxxz_0_xx[i] - rcd_x[i] * t_0_xxxz_0_x[i];

        t_0_xxxy_z_z[i] = t_0_xxxy_0_zz[i] - rcd_z[i] * t_0_xxxy_0_z[i];

        t_0_xxxy_z_y[i] = t_0_xxxy_0_yz[i] - rcd_z[i] * t_0_xxxy_0_y[i];

        t_0_xxxy_z_x[i] = t_0_xxxy_0_xz[i] - rcd_z[i] * t_0_xxxy_0_x[i];

        t_0_xxxy_y_z[i] = t_0_xxxy_0_yz[i] - rcd_y[i] * t_0_xxxy_0_z[i];

        t_0_xxxy_y_y[i] = t_0_xxxy_0_yy[i] - rcd_y[i] * t_0_xxxy_0_y[i];

        t_0_xxxy_y_x[i] = t_0_xxxy_0_xy[i] - rcd_y[i] * t_0_xxxy_0_x[i];

        t_0_xxxy_x_z[i] = t_0_xxxy_0_xz[i] - rcd_x[i] * t_0_xxxy_0_z[i];

        t_0_xxxy_x_y[i] = t_0_xxxy_0_xy[i] - rcd_x[i] * t_0_xxxy_0_y[i];

        t_0_xxxy_x_x[i] = t_0_xxxy_0_xx[i] - rcd_x[i] * t_0_xxxy_0_x[i];

        t_0_xxxx_z_z[i] = t_0_xxxx_0_zz[i] - rcd_z[i] * t_0_xxxx_0_z[i];

        t_0_xxxx_z_y[i] = t_0_xxxx_0_yz[i] - rcd_z[i] * t_0_xxxx_0_y[i];

        t_0_xxxx_z_x[i] = t_0_xxxx_0_xz[i] - rcd_z[i] * t_0_xxxx_0_x[i];

        t_0_xxxx_y_z[i] = t_0_xxxx_0_yz[i] - rcd_y[i] * t_0_xxxx_0_z[i];

        t_0_xxxx_y_y[i] = t_0_xxxx_0_yy[i] - rcd_y[i] * t_0_xxxx_0_y[i];

        t_0_xxxx_y_x[i] = t_0_xxxx_0_xy[i] - rcd_y[i] * t_0_xxxx_0_x[i];

        t_0_xxxx_x_z[i] = t_0_xxxx_0_xz[i] - rcd_x[i] * t_0_xxxx_0_z[i];

        t_0_xxxx_x_y[i] = t_0_xxxx_0_xy[i] - rcd_x[i] * t_0_xxxx_0_y[i];

        t_0_xxxx_x_x[i] = t_0_xxxx_0_xx[i] - rcd_x[i] * t_0_xxxx_0_x[i];
    }
}


} // derirec namespace
