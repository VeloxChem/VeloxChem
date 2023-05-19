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
compHostVRRForSGSD_V0(      BufferHostXY<T>&      intsBufferSGSD,
                      const BufferHostX<int32_t>& intsIndexesSGSD0,
                      const BufferHostXY<T>&      intsBufferSDSD0,
                      const BufferHostX<int32_t>& intsIndexesSDSD0,
                      const BufferHostXY<T>&      intsBufferSDSD1,
                      const BufferHostX<int32_t>& intsIndexesSDSD1,
                      const BufferHostXY<T>&      intsBufferSFSP1,
                      const BufferHostX<int32_t>& intsIndexesSFSP1,
                      const BufferHostXY<T>&      intsBufferSFSD0,
                      const BufferHostX<int32_t>& intsIndexesSFSD0,
                      const BufferHostXY<T>&      intsBufferSFSD1,
                      const BufferHostX<int32_t>& intsIndexesSFSD1,
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

    // set up [SGSD]^(0) integral components

    t_0_zzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(0));

    t_0_zzzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(1));

    t_0_zzzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(2));

    t_0_zzzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(3));

    t_0_zzzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(4));

    t_0_zzzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(5));

    t_0_yzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(6));

    t_0_yzzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(7));

    t_0_yzzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(8));

    t_0_yzzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(9));

    t_0_yzzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(10));

    t_0_yzzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(11));

    t_0_yyzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(12));

    t_0_yyzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(13));

    t_0_yyzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(14));

    t_0_yyzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(15));

    t_0_yyzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(16));

    t_0_yyzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(17));

    t_0_yyyz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(18));

    t_0_yyyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(19));

    t_0_yyyz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(20));

    t_0_yyyz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(21));

    t_0_yyyz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(22));

    t_0_yyyz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(23));

    t_0_yyyy_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(24));

    t_0_yyyy_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(25));

    t_0_yyyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(26));

    t_0_yyyy_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(27));

    t_0_yyyy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(28));

    t_0_yyyy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(29));

    t_0_xzzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(30));

    t_0_xzzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(31));

    t_0_xzzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(32));

    t_0_xzzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(33));

    t_0_xzzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(34));

    t_0_xzzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(35));

    t_0_xyzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(36));

    t_0_xyzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(37));

    t_0_xyzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(38));

    t_0_xyzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(39));

    t_0_xyzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(40));

    t_0_xyzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(41));

    t_0_xyyz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(42));

    t_0_xyyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(43));

    t_0_xyyz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(44));

    t_0_xyyz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(45));

    t_0_xyyz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(46));

    t_0_xyyz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(47));

    t_0_xyyy_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(48));

    t_0_xyyy_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(49));

    t_0_xyyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(50));

    t_0_xyyy_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(51));

    t_0_xyyy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(52));

    t_0_xyyy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(53));

    t_0_xxzz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(54));

    t_0_xxzz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(55));

    t_0_xxzz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(56));

    t_0_xxzz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(57));

    t_0_xxzz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(58));

    t_0_xxzz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(59));

    t_0_xxyz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(60));

    t_0_xxyz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(61));

    t_0_xxyz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(62));

    t_0_xxyz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(63));

    t_0_xxyz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(64));

    t_0_xxyz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(65));

    t_0_xxyy_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(66));

    t_0_xxyy_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(67));

    t_0_xxyy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(68));

    t_0_xxyy_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(69));

    t_0_xxyy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(70));

    t_0_xxyy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(71));

    t_0_xxxz_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(72));

    t_0_xxxz_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(73));

    t_0_xxxz_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(74));

    t_0_xxxz_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(75));

    t_0_xxxz_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(76));

    t_0_xxxz_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(77));

    t_0_xxxy_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(78));

    t_0_xxxy_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(79));

    t_0_xxxy_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(80));

    t_0_xxxy_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(81));

    t_0_xxxy_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(82));

    t_0_xxxy_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(83));

    t_0_xxxx_0_zz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(84));

    t_0_xxxx_0_yz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(85));

    t_0_xxxx_0_yy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(86));

    t_0_xxxx_0_xz_0 = intsBufferSGSD0.data(intsIndexesSGSD0(87));

    t_0_xxxx_0_xy_0 = intsBufferSGSD0.data(intsIndexesSGSD0(88));

    t_0_xxxx_0_xx_0 = intsBufferSGSD0.data(intsIndexesSGSD0(89));

    // set up [SDSD]^(0) integral components

    t_0_zz_0_zz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(0));

    t_0_zz_0_yz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(1));

    t_0_zz_0_yy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(2));

    t_0_zz_0_xz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(3));

    t_0_zz_0_xy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(4));

    t_0_zz_0_xx_0 = intsBufferSDSD0.data(intsIndexesSDSD0(5));

    t_0_yy_0_zz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(6));

    t_0_yy_0_yz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(7));

    t_0_yy_0_yy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(8));

    t_0_yy_0_xz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(9));

    t_0_yy_0_xy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(10));

    t_0_yy_0_xx_0 = intsBufferSDSD0.data(intsIndexesSDSD0(11));

    t_0_xx_0_zz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(12));

    t_0_xx_0_yz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(13));

    t_0_xx_0_yy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(14));

    t_0_xx_0_xz_0 = intsBufferSDSD0.data(intsIndexesSDSD0(15));

    t_0_xx_0_xy_0 = intsBufferSDSD0.data(intsIndexesSDSD0(16));

    t_0_xx_0_xx_0 = intsBufferSDSD0.data(intsIndexesSDSD0(17));

    // set up [SDSD]^(1) integral components

    t_0_zz_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(0));

    t_0_zz_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(1));

    t_0_zz_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(2));

    t_0_zz_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(3));

    t_0_zz_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(4));

    t_0_zz_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(5));

    t_0_yy_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(6));

    t_0_yy_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(7));

    t_0_yy_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(8));

    t_0_yy_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(9));

    t_0_yy_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(10));

    t_0_yy_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(11));

    t_0_xx_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(12));

    t_0_xx_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(13));

    t_0_xx_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(14));

    t_0_xx_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(15));

    t_0_xx_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(16));

    t_0_xx_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(17));

    // set up [SFSP]^(1) integral components

    t_0_zzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(0));

    t_0_zzz_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(1));

    t_0_zzz_0_x_1 = intsBufferSFSP1.data(intsIndexesSFSP1(2));

    t_0_yzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(3));

    t_0_yzz_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(4));

    t_0_yyz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(5));

    t_0_yyy_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(6));

    t_0_yyy_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(7));

    t_0_yyy_0_x_1 = intsBufferSFSP1.data(intsIndexesSFSP1(8));

    t_0_xzz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(9));

    t_0_xyy_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(10));

    t_0_xxz_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(11));

    t_0_xxx_0_z_1 = intsBufferSFSP1.data(intsIndexesSFSP1(12));

    t_0_xxx_0_y_1 = intsBufferSFSP1.data(intsIndexesSFSP1(13));

    t_0_xxx_0_x_1 = intsBufferSFSP1.data(intsIndexesSFSP1(14));

    // set up [SFSD]^(0) integral components

    t_0_zzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(0));

    t_0_zzz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(1));

    t_0_zzz_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(2));

    t_0_zzz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(3));

    t_0_zzz_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(4));

    t_0_zzz_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(5));

    t_0_yzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(6));

    t_0_yzz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(7));

    t_0_yzz_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(8));

    t_0_yzz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(9));

    t_0_yzz_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(10));

    t_0_yzz_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(11));

    t_0_yyz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(12));

    t_0_yyz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(13));

    t_0_yyz_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(14));

    t_0_yyz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(15));

    t_0_yyz_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(16));

    t_0_yyy_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(17));

    t_0_yyy_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(18));

    t_0_yyy_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(19));

    t_0_yyy_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(20));

    t_0_yyy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(21));

    t_0_yyy_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(22));

    t_0_xzz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(23));

    t_0_xzz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(24));

    t_0_xzz_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(25));

    t_0_xzz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(26));

    t_0_xzz_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(27));

    t_0_xyy_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(28));

    t_0_xyy_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(29));

    t_0_xyy_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(30));

    t_0_xyy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(31));

    t_0_xyy_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(32));

    t_0_xxz_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(33));

    t_0_xxz_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(34));

    t_0_xxz_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(35));

    t_0_xxz_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(36));

    t_0_xxz_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(37));

    t_0_xxy_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(38));

    t_0_xxy_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(39));

    t_0_xxy_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(40));

    t_0_xxy_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(41));

    t_0_xxx_0_zz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(42));

    t_0_xxx_0_yz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(43));

    t_0_xxx_0_yy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(44));

    t_0_xxx_0_xz_0 = intsBufferSFSD0.data(intsIndexesSFSD0(45));

    t_0_xxx_0_xy_0 = intsBufferSFSD0.data(intsIndexesSFSD0(46));

    t_0_xxx_0_xx_0 = intsBufferSFSD0.data(intsIndexesSFSD0(47));

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

    t_0_yzz_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(11));

    t_0_yyz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(12));

    t_0_yyz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(13));

    t_0_yyz_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(14));

    t_0_yyz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(15));

    t_0_yyz_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(16));

    t_0_yyy_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(17));

    t_0_yyy_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(18));

    t_0_yyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(19));

    t_0_yyy_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(20));

    t_0_yyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(21));

    t_0_yyy_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(22));

    t_0_xzz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(23));

    t_0_xzz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(24));

    t_0_xzz_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(25));

    t_0_xzz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(26));

    t_0_xzz_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(27));

    t_0_xyy_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(28));

    t_0_xyy_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(29));

    t_0_xyy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(30));

    t_0_xyy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(31));

    t_0_xyy_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(32));

    t_0_xxz_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(33));

    t_0_xxz_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(34));

    t_0_xxz_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(35));

    t_0_xxz_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(36));

    t_0_xxz_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(37));

    t_0_xxy_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(38));

    t_0_xxy_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(39));

    t_0_xxy_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(40));

    t_0_xxy_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(41));

    t_0_xxx_0_zz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(42));

    t_0_xxx_0_yz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(43));

    t_0_xxx_0_yy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(44));

    t_0_xxx_0_xz_1 = intsBufferSFSD1.data(intsIndexesSFSD1(45));

    t_0_xxx_0_xy_1 = intsBufferSFSD1.data(intsIndexesSFSD1(46));

    t_0_xxx_0_xx_1 = intsBufferSFSD1.data(intsIndexesSFSD1(47));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xzzz_0_xx_0, t_0_xzzz_0_xy_0, t_0_xzzz_0_xz_0,\
                               t_0_xzzz_0_yy_0, t_0_xzzz_0_yz_0, t_0_xzzz_0_zz_0, t_0_yy_0_xx_0,\
                               t_0_yy_0_xx_1, t_0_yy_0_xy_0, t_0_yy_0_xy_1, t_0_yy_0_xz_0,\
                               t_0_yy_0_xz_1, t_0_yy_0_yy_0, t_0_yy_0_yy_1, t_0_yy_0_yz_0,\
                               t_0_yy_0_yz_1, t_0_yy_0_zz_0, t_0_yy_0_zz_1, t_0_yyy_0_x_1,\
                               t_0_yyy_0_xx_0, t_0_yyy_0_xx_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1,\
                               t_0_yyy_0_xz_0, t_0_yyy_0_xz_1, t_0_yyy_0_y_1, t_0_yyy_0_yy_0,\
                               t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1, t_0_yyy_0_z_1,\
                               t_0_yyy_0_zz_0, t_0_yyy_0_zz_1, t_0_yyyy_0_xx_0, t_0_yyyy_0_xy_0,\
                               t_0_yyyy_0_xz_0, t_0_yyyy_0_yy_0, t_0_yyyy_0_yz_0, t_0_yyyy_0_zz_0,\
                               t_0_yyyz_0_xx_0, t_0_yyyz_0_xy_0, t_0_yyyz_0_xz_0, t_0_yyyz_0_yy_0,\
                               t_0_yyyz_0_yz_0, t_0_yyyz_0_zz_0, t_0_yyz_0_xy_0, t_0_yyz_0_xy_1,\
                               t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyzz_0_xx_0, t_0_yyzz_0_xy_0,\
                               t_0_yyzz_0_xz_0, t_0_yyzz_0_yy_0, t_0_yyzz_0_yz_0, t_0_yyzz_0_zz_0,\
                               t_0_yzz_0_xx_0, t_0_yzz_0_xx_1, t_0_yzz_0_xz_0, t_0_yzz_0_xz_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_z_1, t_0_yzz_0_zz_0,\
                               t_0_yzz_0_zz_1, t_0_yzzz_0_xx_0, t_0_yzzz_0_xy_0, t_0_yzzz_0_xz_0,\
                               t_0_yzzz_0_yy_0, t_0_yzzz_0_yz_0, t_0_yzzz_0_zz_0, t_0_zz_0_xx_0,\
                               t_0_zz_0_xx_1, t_0_zz_0_xy_0, t_0_zz_0_xy_1, t_0_zz_0_xz_0,\
                               t_0_zz_0_xz_1, t_0_zz_0_yy_0, t_0_zz_0_yy_1, t_0_zz_0_yz_0,\
                               t_0_zz_0_yz_1, t_0_zz_0_zz_0, t_0_zz_0_zz_1, t_0_zzz_0_x_1,\
                               t_0_zzz_0_xx_0, t_0_zzz_0_xx_1, t_0_zzz_0_xy_0, t_0_zzz_0_xy_1,\
                               t_0_zzz_0_xz_0, t_0_zzz_0_xz_1, t_0_zzz_0_y_1, t_0_zzz_0_yy_0,\
                               t_0_zzz_0_yy_1, t_0_zzz_0_yz_0, t_0_zzz_0_yz_1, t_0_zzz_0_z_1,\
                               t_0_zzz_0_zz_0, t_0_zzz_0_zz_1, t_0_zzzz_0_xx_0, t_0_zzzz_0_xy_0,\
                               t_0_zzzz_0_xz_0, t_0_zzzz_0_yy_0, t_0_zzzz_0_yz_0, t_0_zzzz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zz_0[i] += rpb_z[i] * t_0_zzz_0_zz_0[i] + rwp_z[i] * t_0_zzz_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zz_1[i] + fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_zzzz_0_yz_0[i] += rpb_z[i] * t_0_zzz_0_yz_0[i] + rwp_z[i] * t_0_zzz_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_zzzz_0_yy_0[i] += rpb_z[i] * t_0_zzz_0_yy_0[i] + rwp_z[i] * t_0_zzz_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzzz_0_xz_0[i] += rpb_z[i] * t_0_zzz_0_xz_0[i] + rwp_z[i] * t_0_zzz_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_x_1[i];

            t_0_zzzz_0_xy_0[i] += rpb_z[i] * t_0_zzz_0_xy_0[i] + rwp_z[i] * t_0_zzz_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzzz_0_xx_0[i] += rpb_z[i] * t_0_zzz_0_xx_0[i] + rwp_z[i] * t_0_zzz_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzzz_0_zz_0[i] += rpb_y[i] * t_0_zzz_0_zz_0[i] + rwp_y[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yz_0[i] += rpb_y[i] * t_0_zzz_0_yz_0[i] + rwp_y[i] * t_0_zzz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yzzz_0_yy_0[i] += rpb_y[i] * t_0_zzz_0_yy_0[i] + rwp_y[i] * t_0_zzz_0_yy_1[i] + fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_yzzz_0_xz_0[i] += rpb_y[i] * t_0_zzz_0_xz_0[i] + rwp_y[i] * t_0_zzz_0_xz_1[i];

            t_0_yzzz_0_xy_0[i] += rpb_y[i] * t_0_zzz_0_xy_0[i] + rwp_y[i] * t_0_zzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_x_1[i];

            t_0_yzzz_0_xx_0[i] += rpb_y[i] * t_0_zzz_0_xx_0[i] + rwp_y[i] * t_0_zzz_0_xx_1[i];

            t_0_yyzz_0_zz_0[i] += rpb_y[i] * t_0_yzz_0_zz_0[i] + rwp_y[i] * t_0_yzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_yyzz_0_yz_0[i] += rpb_y[i] * t_0_yzz_0_yz_0[i] + rwp_y[i] * t_0_yzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_z_1[i];

            t_0_yyzz_0_yy_0[i] += rpb_z[i] * t_0_yyz_0_yy_0[i] + rwp_z[i] * t_0_yyz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyzz_0_xz_0[i] += rpb_y[i] * t_0_yzz_0_xz_0[i] + rwp_y[i] * t_0_yzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i];

            t_0_yyzz_0_xy_0[i] += rpb_z[i] * t_0_yyz_0_xy_0[i] + rwp_z[i] * t_0_yyz_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyzz_0_xx_0[i] += rpb_y[i] * t_0_yzz_0_xx_0[i] + rwp_y[i] * t_0_yzz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xx_1[i];

            t_0_yyyz_0_zz_0[i] += rpb_z[i] * t_0_yyy_0_zz_0[i] + rwp_z[i] * t_0_yyy_0_zz_1[i] + fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_yyyz_0_yz_0[i] += rpb_z[i] * t_0_yyy_0_yz_0[i] + rwp_z[i] * t_0_yyy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyz_0_yy_0[i] += rpb_z[i] * t_0_yyy_0_yy_0[i] + rwp_z[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyz_0_xz_0[i] += rpb_z[i] * t_0_yyy_0_xz_0[i] + rwp_z[i] * t_0_yyy_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_yyyz_0_xy_0[i] += rpb_z[i] * t_0_yyy_0_xy_0[i] + rwp_z[i] * t_0_yyy_0_xy_1[i];

            t_0_yyyz_0_xx_0[i] += rpb_z[i] * t_0_yyy_0_xx_0[i] + rwp_z[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyy_0_zz_0[i] += rpb_y[i] * t_0_yyy_0_zz_0[i] + rwp_y[i] * t_0_yyy_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyyy_0_yz_0[i] += rpb_y[i] * t_0_yyy_0_yz_0[i] + rwp_y[i] * t_0_yyy_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_yyyy_0_yy_0[i] += rpb_y[i] * t_0_yyy_0_yy_0[i] + rwp_y[i] * t_0_yyy_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yy_1[i] + fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyy_0_xz_0[i] += rpb_y[i] * t_0_yyy_0_xz_0[i] + rwp_y[i] * t_0_yyy_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyyy_0_xy_0[i] += rpb_y[i] * t_0_yyy_0_xy_0[i] + rwp_y[i] * t_0_yyy_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_yyyy_0_xx_0[i] += rpb_y[i] * t_0_yyy_0_xx_0[i] + rwp_y[i] * t_0_yyy_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xx_1[i];

            t_0_xzzz_0_zz_0[i] += rpb_x[i] * t_0_zzz_0_zz_0[i] + rwp_x[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_yz_0[i] += rpb_x[i] * t_0_zzz_0_yz_0[i] + rwp_x[i] * t_0_zzz_0_yz_1[i];

            t_0_xzzz_0_yy_0[i] += rpb_x[i] * t_0_zzz_0_yy_0[i] + rwp_x[i] * t_0_zzz_0_yy_1[i];

            t_0_xzzz_0_xz_0[i] += rpb_x[i] * t_0_zzz_0_xz_0[i] + rwp_x[i] * t_0_zzz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_xzzz_0_xy_0[i] += rpb_x[i] * t_0_zzz_0_xy_0[i] + rwp_x[i] * t_0_zzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_xzzz_0_xx_0[i] += rpb_x[i] * t_0_zzz_0_xx_0[i] + rwp_x[i] * t_0_zzz_0_xx_1[i] + fze_0[i] * t_0_zzz_0_x_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xx_0_xy_0, t_0_xx_0_xy_1,\
                               t_0_xx_0_xz_0, t_0_xx_0_xz_1, t_0_xxy_0_xx_0, t_0_xxy_0_xx_1,\
                               t_0_xxy_0_xy_0, t_0_xxy_0_xy_1, t_0_xxy_0_xz_0, t_0_xxy_0_xz_1,\
                               t_0_xxy_0_yy_0, t_0_xxy_0_yy_1, t_0_xxyy_0_xx_0, t_0_xxyy_0_xy_0,\
                               t_0_xxyy_0_xz_0, t_0_xxyy_0_yy_0, t_0_xxyy_0_yz_0, t_0_xxyy_0_zz_0,\
                               t_0_xxyz_0_xx_0, t_0_xxyz_0_xy_0, t_0_xxyz_0_xz_0, t_0_xxyz_0_yy_0,\
                               t_0_xxyz_0_yz_0, t_0_xxyz_0_zz_0, t_0_xxz_0_xx_0, t_0_xxz_0_xx_1,\
                               t_0_xxz_0_xy_0, t_0_xxz_0_xy_1, t_0_xxz_0_xz_0, t_0_xxz_0_xz_1,\
                               t_0_xxz_0_yz_0, t_0_xxz_0_yz_1, t_0_xxz_0_z_1, t_0_xxz_0_zz_0,\
                               t_0_xxz_0_zz_1, t_0_xxzz_0_xx_0, t_0_xxzz_0_xy_0, t_0_xxzz_0_xz_0,\
                               t_0_xxzz_0_yy_0, t_0_xxzz_0_yz_0, t_0_xxzz_0_zz_0, t_0_xyy_0_xx_0,\
                               t_0_xyy_0_xx_1, t_0_xyy_0_xy_0, t_0_xyy_0_xy_1, t_0_xyy_0_y_1,\
                               t_0_xyy_0_yy_0, t_0_xyy_0_yy_1, t_0_xyy_0_yz_0, t_0_xyy_0_yz_1,\
                               t_0_xyy_0_zz_0, t_0_xyy_0_zz_1, t_0_xyyy_0_xx_0, t_0_xyyy_0_xy_0,\
                               t_0_xyyy_0_xz_0, t_0_xyyy_0_yy_0, t_0_xyyy_0_yz_0, t_0_xyyy_0_zz_0,\
                               t_0_xyyz_0_xx_0, t_0_xyyz_0_xy_0, t_0_xyyz_0_xz_0, t_0_xyyz_0_yy_0,\
                               t_0_xyyz_0_yz_0, t_0_xyyz_0_zz_0, t_0_xyzz_0_xx_0, t_0_xyzz_0_xy_0,\
                               t_0_xyzz_0_xz_0, t_0_xyzz_0_yy_0, t_0_xyzz_0_yz_0, t_0_xyzz_0_zz_0,\
                               t_0_xzz_0_xx_0, t_0_xzz_0_xx_1, t_0_xzz_0_xz_0, t_0_xzz_0_xz_1,\
                               t_0_xzz_0_yy_0, t_0_xzz_0_yy_1, t_0_xzz_0_yz_0, t_0_xzz_0_yz_1,\
                               t_0_xzz_0_z_1, t_0_xzz_0_zz_0, t_0_xzz_0_zz_1, t_0_yy_0_xy_0,\
                               t_0_yy_0_xy_1, t_0_yy_0_yy_0, t_0_yy_0_yy_1, t_0_yy_0_yz_0,\
                               t_0_yy_0_yz_1, t_0_yy_0_zz_0, t_0_yy_0_zz_1, t_0_yyy_0_x_1,\
                               t_0_yyy_0_xx_0, t_0_yyy_0_xx_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1,\
                               t_0_yyy_0_xz_0, t_0_yyy_0_xz_1, t_0_yyy_0_y_1, t_0_yyy_0_yy_0,\
                               t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1, t_0_yyy_0_z_1,\
                               t_0_yyy_0_zz_0, t_0_yyy_0_zz_1, t_0_yyz_0_xz_0, t_0_yyz_0_xz_1,\
                               t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyz_0_yz_0, t_0_yyz_0_yz_1,\
                               t_0_yyz_0_z_1, t_0_yyz_0_zz_0, t_0_yyz_0_zz_1, t_0_yzz_0_xy_0,\
                               t_0_yzz_0_xy_1, t_0_yzz_0_y_1, t_0_yzz_0_yy_0, t_0_yzz_0_yy_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_zz_0, t_0_yzz_0_zz_1,\
                               t_0_zz_0_xz_0, t_0_zz_0_xz_1, t_0_zz_0_yy_0, t_0_zz_0_yy_1,\
                               t_0_zz_0_yz_0, t_0_zz_0_yz_1, t_0_zz_0_zz_0, t_0_zz_0_zz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyzz_0_zz_0[i] += rpb_x[i] * t_0_yzz_0_zz_0[i] + rwp_x[i] * t_0_yzz_0_zz_1[i];

            t_0_xyzz_0_yz_0[i] += rpb_x[i] * t_0_yzz_0_yz_0[i] + rwp_x[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_yy_0[i] += rpb_x[i] * t_0_yzz_0_yy_0[i] + rwp_x[i] * t_0_yzz_0_yy_1[i];

            t_0_xyzz_0_xz_0[i] += rpb_y[i] * t_0_xzz_0_xz_0[i] + rwp_y[i] * t_0_xzz_0_xz_1[i];

            t_0_xyzz_0_xy_0[i] += rpb_x[i] * t_0_yzz_0_xy_0[i] + rwp_x[i] * t_0_yzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_y_1[i];

            t_0_xyzz_0_xx_0[i] += rpb_y[i] * t_0_xzz_0_xx_0[i] + rwp_y[i] * t_0_xzz_0_xx_1[i];

            t_0_xyyz_0_zz_0[i] += rpb_x[i] * t_0_yyz_0_zz_0[i] + rwp_x[i] * t_0_yyz_0_zz_1[i];

            t_0_xyyz_0_yz_0[i] += rpb_x[i] * t_0_yyz_0_yz_0[i] + rwp_x[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_yy_0[i] += rpb_x[i] * t_0_yyz_0_yy_0[i] + rwp_x[i] * t_0_yyz_0_yy_1[i];

            t_0_xyyz_0_xz_0[i] += rpb_x[i] * t_0_yyz_0_xz_0[i] + rwp_x[i] * t_0_yyz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_z_1[i];

            t_0_xyyz_0_xy_0[i] += rpb_z[i] * t_0_xyy_0_xy_0[i] + rwp_z[i] * t_0_xyy_0_xy_1[i];

            t_0_xyyz_0_xx_0[i] += rpb_z[i] * t_0_xyy_0_xx_0[i] + rwp_z[i] * t_0_xyy_0_xx_1[i];

            t_0_xyyy_0_zz_0[i] += rpb_x[i] * t_0_yyy_0_zz_0[i] + rwp_x[i] * t_0_yyy_0_zz_1[i];

            t_0_xyyy_0_yz_0[i] += rpb_x[i] * t_0_yyy_0_yz_0[i] + rwp_x[i] * t_0_yyy_0_yz_1[i];

            t_0_xyyy_0_yy_0[i] += rpb_x[i] * t_0_yyy_0_yy_0[i] + rwp_x[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xz_0[i] += rpb_x[i] * t_0_yyy_0_xz_0[i] + rwp_x[i] * t_0_yyy_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_xyyy_0_xy_0[i] += rpb_x[i] * t_0_yyy_0_xy_0[i] + rwp_x[i] * t_0_yyy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xyyy_0_xx_0[i] += rpb_x[i] * t_0_yyy_0_xx_0[i] + rwp_x[i] * t_0_yyy_0_xx_1[i] + fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_xxzz_0_zz_0[i] += rpb_x[i] * t_0_xzz_0_zz_0[i] + rwp_x[i] * t_0_xzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_xxzz_0_yz_0[i] += rpb_x[i] * t_0_xzz_0_yz_0[i] + rwp_x[i] * t_0_xzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i];

            t_0_xxzz_0_yy_0[i] += rpb_x[i] * t_0_xzz_0_yy_0[i] + rwp_x[i] * t_0_xzz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yy_1[i];

            t_0_xxzz_0_xz_0[i] += rpb_x[i] * t_0_xzz_0_xz_0[i] + rwp_x[i] * t_0_xzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_z_1[i];

            t_0_xxzz_0_xy_0[i] += rpb_z[i] * t_0_xxz_0_xy_0[i] + rwp_z[i] * t_0_xxz_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxzz_0_xx_0[i] += rpb_z[i] * t_0_xxz_0_xx_0[i] + rwp_z[i] * t_0_xxz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxyz_0_zz_0[i] += rpb_y[i] * t_0_xxz_0_zz_0[i] + rwp_y[i] * t_0_xxz_0_zz_1[i];

            t_0_xxyz_0_yz_0[i] += rpb_y[i] * t_0_xxz_0_yz_0[i] + rwp_y[i] * t_0_xxz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_z_1[i];

            t_0_xxyz_0_yy_0[i] += rpb_z[i] * t_0_xxy_0_yy_0[i] + rwp_z[i] * t_0_xxy_0_yy_1[i];

            t_0_xxyz_0_xz_0[i] += rpb_y[i] * t_0_xxz_0_xz_0[i] + rwp_y[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xy_0[i] += rpb_z[i] * t_0_xxy_0_xy_0[i] + rwp_z[i] * t_0_xxy_0_xy_1[i];

            t_0_xxyz_0_xx_0[i] += rpb_y[i] * t_0_xxz_0_xx_0[i] + rwp_y[i] * t_0_xxz_0_xx_1[i];

            t_0_xxyy_0_zz_0[i] += rpb_x[i] * t_0_xyy_0_zz_0[i] + rwp_x[i] * t_0_xyy_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zz_1[i];

            t_0_xxyy_0_yz_0[i] += rpb_x[i] * t_0_xyy_0_yz_0[i] + rwp_x[i] * t_0_xyy_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yz_1[i];

            t_0_xxyy_0_yy_0[i] += rpb_x[i] * t_0_xyy_0_yy_0[i] + rwp_x[i] * t_0_xyy_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_xxyy_0_xz_0[i] += rpb_y[i] * t_0_xxy_0_xz_0[i] + rwp_y[i] * t_0_xxy_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxyy_0_xy_0[i] += rpb_x[i] * t_0_xyy_0_xy_0[i] + rwp_x[i] * t_0_xyy_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_y_1[i];

            t_0_xxyy_0_xx_0[i] += rpb_y[i] * t_0_xxy_0_xx_0[i] + rwp_y[i] * t_0_xxy_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xx_0_xy_0, t_0_xx_0_xy_1,\
                               t_0_xx_0_xz_0, t_0_xx_0_xz_1, t_0_xx_0_yy_0, t_0_xx_0_yy_1,\
                               t_0_xx_0_yz_0, t_0_xx_0_yz_1, t_0_xx_0_zz_0, t_0_xx_0_zz_1,\
                               t_0_xxx_0_x_1, t_0_xxx_0_xx_0, t_0_xxx_0_xx_1, t_0_xxx_0_xy_0,\
                               t_0_xxx_0_xy_1, t_0_xxx_0_xz_0, t_0_xxx_0_xz_1, t_0_xxx_0_y_1,\
                               t_0_xxx_0_yy_0, t_0_xxx_0_yy_1, t_0_xxx_0_yz_0, t_0_xxx_0_yz_1,\
                               t_0_xxx_0_z_1, t_0_xxx_0_zz_0, t_0_xxx_0_zz_1, t_0_xxxx_0_xx_0,\
                               t_0_xxxx_0_xy_0, t_0_xxxx_0_xz_0, t_0_xxxx_0_yy_0, t_0_xxxx_0_yz_0,\
                               t_0_xxxx_0_zz_0, t_0_xxxy_0_xx_0, t_0_xxxy_0_xy_0, t_0_xxxy_0_xz_0,\
                               t_0_xxxy_0_yy_0, t_0_xxxy_0_yz_0, t_0_xxxy_0_zz_0, t_0_xxxz_0_xx_0,\
                               t_0_xxxz_0_xy_0, t_0_xxxz_0_xz_0, t_0_xxxz_0_yy_0, t_0_xxxz_0_yz_0,\
                               t_0_xxxz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxz_0_zz_0[i] += rpb_z[i] * t_0_xxx_0_zz_0[i] + rwp_z[i] * t_0_xxx_0_zz_1[i] + fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxz_0_yz_0[i] += rpb_z[i] * t_0_xxx_0_yz_0[i] + rwp_z[i] * t_0_xxx_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxz_0_yy_0[i] += rpb_z[i] * t_0_xxx_0_yy_0[i] + rwp_z[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxz_0_xz_0[i] += rpb_z[i] * t_0_xxx_0_xz_0[i] + rwp_z[i] * t_0_xxx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxz_0_xy_0[i] += rpb_z[i] * t_0_xxx_0_xy_0[i] + rwp_z[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxz_0_xx_0[i] += rpb_z[i] * t_0_xxx_0_xx_0[i] + rwp_z[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_zz_0[i] += rpb_y[i] * t_0_xxx_0_zz_0[i] + rwp_y[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxy_0_yz_0[i] += rpb_y[i] * t_0_xxx_0_yz_0[i] + rwp_y[i] * t_0_xxx_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxy_0_yy_0[i] += rpb_y[i] * t_0_xxx_0_yy_0[i] + rwp_y[i] * t_0_xxx_0_yy_1[i] + fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxy_0_xz_0[i] += rpb_y[i] * t_0_xxx_0_xz_0[i] + rwp_y[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxy_0_xy_0[i] += rpb_y[i] * t_0_xxx_0_xy_0[i] + rwp_y[i] * t_0_xxx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxy_0_xx_0[i] += rpb_y[i] * t_0_xxx_0_xx_0[i] + rwp_y[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxx_0_zz_0[i] += rpb_x[i] * t_0_xxx_0_zz_0[i] + rwp_x[i] * t_0_xxx_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxxx_0_yz_0[i] += rpb_x[i] * t_0_xxx_0_yz_0[i] + rwp_x[i] * t_0_xxx_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxxx_0_yy_0[i] += rpb_x[i] * t_0_xxx_0_yy_0[i] + rwp_x[i] * t_0_xxx_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxxx_0_xz_0[i] += rpb_x[i] * t_0_xxx_0_xz_0[i] + rwp_x[i] * t_0_xxx_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxx_0_xy_0[i] += rpb_x[i] * t_0_xxx_0_xy_0[i] + rwp_x[i] * t_0_xxx_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxx_0_xx_0[i] += rpb_x[i] * t_0_xxx_0_xx_0[i] + rwp_x[i] * t_0_xxx_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xx_1[i] + fze_0[i] * t_0_xxx_0_x_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xzzz_0_xx_0, t_0_xzzz_0_xy_0, t_0_xzzz_0_xz_0,\
                               t_0_xzzz_0_yy_0, t_0_xzzz_0_yz_0, t_0_xzzz_0_zz_0, t_0_yy_0_xx_0,\
                               t_0_yy_0_xx_1, t_0_yy_0_xy_0, t_0_yy_0_xy_1, t_0_yy_0_xz_0,\
                               t_0_yy_0_xz_1, t_0_yy_0_yy_0, t_0_yy_0_yy_1, t_0_yy_0_yz_0,\
                               t_0_yy_0_yz_1, t_0_yy_0_zz_0, t_0_yy_0_zz_1, t_0_yyy_0_x_1,\
                               t_0_yyy_0_xx_0, t_0_yyy_0_xx_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1,\
                               t_0_yyy_0_xz_0, t_0_yyy_0_xz_1, t_0_yyy_0_y_1, t_0_yyy_0_yy_0,\
                               t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1, t_0_yyy_0_z_1,\
                               t_0_yyy_0_zz_0, t_0_yyy_0_zz_1, t_0_yyyy_0_xx_0, t_0_yyyy_0_xy_0,\
                               t_0_yyyy_0_xz_0, t_0_yyyy_0_yy_0, t_0_yyyy_0_yz_0, t_0_yyyy_0_zz_0,\
                               t_0_yyyz_0_xx_0, t_0_yyyz_0_xy_0, t_0_yyyz_0_xz_0, t_0_yyyz_0_yy_0,\
                               t_0_yyyz_0_yz_0, t_0_yyyz_0_zz_0, t_0_yyz_0_xy_0, t_0_yyz_0_xy_1,\
                               t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyzz_0_xx_0, t_0_yyzz_0_xy_0,\
                               t_0_yyzz_0_xz_0, t_0_yyzz_0_yy_0, t_0_yyzz_0_yz_0, t_0_yyzz_0_zz_0,\
                               t_0_yzz_0_xx_0, t_0_yzz_0_xx_1, t_0_yzz_0_xz_0, t_0_yzz_0_xz_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_z_1, t_0_yzz_0_zz_0,\
                               t_0_yzz_0_zz_1, t_0_yzzz_0_xx_0, t_0_yzzz_0_xy_0, t_0_yzzz_0_xz_0,\
                               t_0_yzzz_0_yy_0, t_0_yzzz_0_yz_0, t_0_yzzz_0_zz_0, t_0_zz_0_xx_0,\
                               t_0_zz_0_xx_1, t_0_zz_0_xy_0, t_0_zz_0_xy_1, t_0_zz_0_xz_0,\
                               t_0_zz_0_xz_1, t_0_zz_0_yy_0, t_0_zz_0_yy_1, t_0_zz_0_yz_0,\
                               t_0_zz_0_yz_1, t_0_zz_0_zz_0, t_0_zz_0_zz_1, t_0_zzz_0_x_1,\
                               t_0_zzz_0_xx_0, t_0_zzz_0_xx_1, t_0_zzz_0_xy_0, t_0_zzz_0_xy_1,\
                               t_0_zzz_0_xz_0, t_0_zzz_0_xz_1, t_0_zzz_0_y_1, t_0_zzz_0_yy_0,\
                               t_0_zzz_0_yy_1, t_0_zzz_0_yz_0, t_0_zzz_0_yz_1, t_0_zzz_0_z_1,\
                               t_0_zzz_0_zz_0, t_0_zzz_0_zz_1, t_0_zzzz_0_xx_0, t_0_zzzz_0_xy_0,\
                               t_0_zzzz_0_xz_0, t_0_zzzz_0_yy_0, t_0_zzzz_0_yz_0, t_0_zzzz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_zz_0[i] = rpb_z[i] * t_0_zzz_0_zz_0[i] + rwp_z[i] * t_0_zzz_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_zz_1[i] + fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_zzzz_0_yz_0[i] = rpb_z[i] * t_0_zzz_0_yz_0[i] + rwp_z[i] * t_0_zzz_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_zzzz_0_yy_0[i] = rpb_z[i] * t_0_zzz_0_yy_0[i] + rwp_z[i] * t_0_zzz_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzzz_0_xz_0[i] = rpb_z[i] * t_0_zzz_0_xz_0[i] + rwp_z[i] * t_0_zzz_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_x_1[i];

            t_0_zzzz_0_xy_0[i] = rpb_z[i] * t_0_zzz_0_xy_0[i] + rwp_z[i] * t_0_zzz_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzzz_0_xx_0[i] = rpb_z[i] * t_0_zzz_0_xx_0[i] + rwp_z[i] * t_0_zzz_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzzz_0_zz_0[i] = rpb_y[i] * t_0_zzz_0_zz_0[i] + rwp_y[i] * t_0_zzz_0_zz_1[i];

            t_0_yzzz_0_yz_0[i] = rpb_y[i] * t_0_zzz_0_yz_0[i] + rwp_y[i] * t_0_zzz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_yzzz_0_yy_0[i] = rpb_y[i] * t_0_zzz_0_yy_0[i] + rwp_y[i] * t_0_zzz_0_yy_1[i] + fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_yzzz_0_xz_0[i] = rpb_y[i] * t_0_zzz_0_xz_0[i] + rwp_y[i] * t_0_zzz_0_xz_1[i];

            t_0_yzzz_0_xy_0[i] = rpb_y[i] * t_0_zzz_0_xy_0[i] + rwp_y[i] * t_0_zzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_x_1[i];

            t_0_yzzz_0_xx_0[i] = rpb_y[i] * t_0_zzz_0_xx_0[i] + rwp_y[i] * t_0_zzz_0_xx_1[i];

            t_0_yyzz_0_zz_0[i] = rpb_y[i] * t_0_yzz_0_zz_0[i] + rwp_y[i] * t_0_yzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_yyzz_0_yz_0[i] = rpb_y[i] * t_0_yzz_0_yz_0[i] + rwp_y[i] * t_0_yzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_z_1[i];

            t_0_yyzz_0_yy_0[i] = rpb_z[i] * t_0_yyz_0_yy_0[i] + rwp_z[i] * t_0_yyz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyzz_0_xz_0[i] = rpb_y[i] * t_0_yzz_0_xz_0[i] + rwp_y[i] * t_0_yzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i];

            t_0_yyzz_0_xy_0[i] = rpb_z[i] * t_0_yyz_0_xy_0[i] + rwp_z[i] * t_0_yyz_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyzz_0_xx_0[i] = rpb_y[i] * t_0_yzz_0_xx_0[i] + rwp_y[i] * t_0_yzz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xx_1[i];

            t_0_yyyz_0_zz_0[i] = rpb_z[i] * t_0_yyy_0_zz_0[i] + rwp_z[i] * t_0_yyy_0_zz_1[i] + fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_yyyz_0_yz_0[i] = rpb_z[i] * t_0_yyy_0_yz_0[i] + rwp_z[i] * t_0_yyy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyz_0_yy_0[i] = rpb_z[i] * t_0_yyy_0_yy_0[i] + rwp_z[i] * t_0_yyy_0_yy_1[i];

            t_0_yyyz_0_xz_0[i] = rpb_z[i] * t_0_yyy_0_xz_0[i] + rwp_z[i] * t_0_yyy_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_yyyz_0_xy_0[i] = rpb_z[i] * t_0_yyy_0_xy_0[i] + rwp_z[i] * t_0_yyy_0_xy_1[i];

            t_0_yyyz_0_xx_0[i] = rpb_z[i] * t_0_yyy_0_xx_0[i] + rwp_z[i] * t_0_yyy_0_xx_1[i];

            t_0_yyyy_0_zz_0[i] = rpb_y[i] * t_0_yyy_0_zz_0[i] + rwp_y[i] * t_0_yyy_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyyy_0_yz_0[i] = rpb_y[i] * t_0_yyy_0_yz_0[i] + rwp_y[i] * t_0_yyy_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_yyyy_0_yy_0[i] = rpb_y[i] * t_0_yyy_0_yy_0[i] + rwp_y[i] * t_0_yyy_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_yy_1[i] + fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_yyyy_0_xz_0[i] = rpb_y[i] * t_0_yyy_0_xz_0[i] + rwp_y[i] * t_0_yyy_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyyy_0_xy_0[i] = rpb_y[i] * t_0_yyy_0_xy_0[i] + rwp_y[i] * t_0_yyy_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_yyyy_0_xx_0[i] = rpb_y[i] * t_0_yyy_0_xx_0[i] + rwp_y[i] * t_0_yyy_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_xx_1[i];

            t_0_xzzz_0_zz_0[i] = rpb_x[i] * t_0_zzz_0_zz_0[i] + rwp_x[i] * t_0_zzz_0_zz_1[i];

            t_0_xzzz_0_yz_0[i] = rpb_x[i] * t_0_zzz_0_yz_0[i] + rwp_x[i] * t_0_zzz_0_yz_1[i];

            t_0_xzzz_0_yy_0[i] = rpb_x[i] * t_0_zzz_0_yy_0[i] + rwp_x[i] * t_0_zzz_0_yy_1[i];

            t_0_xzzz_0_xz_0[i] = rpb_x[i] * t_0_zzz_0_xz_0[i] + rwp_x[i] * t_0_zzz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_z_1[i];

            t_0_xzzz_0_xy_0[i] = rpb_x[i] * t_0_zzz_0_xy_0[i] + rwp_x[i] * t_0_zzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_zzz_0_y_1[i];

            t_0_xzzz_0_xx_0[i] = rpb_x[i] * t_0_zzz_0_xx_0[i] + rwp_x[i] * t_0_zzz_0_xx_1[i] + fze_0[i] * t_0_zzz_0_x_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xx_0_xy_0, t_0_xx_0_xy_1,\
                               t_0_xx_0_xz_0, t_0_xx_0_xz_1, t_0_xxy_0_xx_0, t_0_xxy_0_xx_1,\
                               t_0_xxy_0_xy_0, t_0_xxy_0_xy_1, t_0_xxy_0_xz_0, t_0_xxy_0_xz_1,\
                               t_0_xxy_0_yy_0, t_0_xxy_0_yy_1, t_0_xxyy_0_xx_0, t_0_xxyy_0_xy_0,\
                               t_0_xxyy_0_xz_0, t_0_xxyy_0_yy_0, t_0_xxyy_0_yz_0, t_0_xxyy_0_zz_0,\
                               t_0_xxyz_0_xx_0, t_0_xxyz_0_xy_0, t_0_xxyz_0_xz_0, t_0_xxyz_0_yy_0,\
                               t_0_xxyz_0_yz_0, t_0_xxyz_0_zz_0, t_0_xxz_0_xx_0, t_0_xxz_0_xx_1,\
                               t_0_xxz_0_xy_0, t_0_xxz_0_xy_1, t_0_xxz_0_xz_0, t_0_xxz_0_xz_1,\
                               t_0_xxz_0_yz_0, t_0_xxz_0_yz_1, t_0_xxz_0_z_1, t_0_xxz_0_zz_0,\
                               t_0_xxz_0_zz_1, t_0_xxzz_0_xx_0, t_0_xxzz_0_xy_0, t_0_xxzz_0_xz_0,\
                               t_0_xxzz_0_yy_0, t_0_xxzz_0_yz_0, t_0_xxzz_0_zz_0, t_0_xyy_0_xx_0,\
                               t_0_xyy_0_xx_1, t_0_xyy_0_xy_0, t_0_xyy_0_xy_1, t_0_xyy_0_y_1,\
                               t_0_xyy_0_yy_0, t_0_xyy_0_yy_1, t_0_xyy_0_yz_0, t_0_xyy_0_yz_1,\
                               t_0_xyy_0_zz_0, t_0_xyy_0_zz_1, t_0_xyyy_0_xx_0, t_0_xyyy_0_xy_0,\
                               t_0_xyyy_0_xz_0, t_0_xyyy_0_yy_0, t_0_xyyy_0_yz_0, t_0_xyyy_0_zz_0,\
                               t_0_xyyz_0_xx_0, t_0_xyyz_0_xy_0, t_0_xyyz_0_xz_0, t_0_xyyz_0_yy_0,\
                               t_0_xyyz_0_yz_0, t_0_xyyz_0_zz_0, t_0_xyzz_0_xx_0, t_0_xyzz_0_xy_0,\
                               t_0_xyzz_0_xz_0, t_0_xyzz_0_yy_0, t_0_xyzz_0_yz_0, t_0_xyzz_0_zz_0,\
                               t_0_xzz_0_xx_0, t_0_xzz_0_xx_1, t_0_xzz_0_xz_0, t_0_xzz_0_xz_1,\
                               t_0_xzz_0_yy_0, t_0_xzz_0_yy_1, t_0_xzz_0_yz_0, t_0_xzz_0_yz_1,\
                               t_0_xzz_0_z_1, t_0_xzz_0_zz_0, t_0_xzz_0_zz_1, t_0_yy_0_xy_0,\
                               t_0_yy_0_xy_1, t_0_yy_0_yy_0, t_0_yy_0_yy_1, t_0_yy_0_yz_0,\
                               t_0_yy_0_yz_1, t_0_yy_0_zz_0, t_0_yy_0_zz_1, t_0_yyy_0_x_1,\
                               t_0_yyy_0_xx_0, t_0_yyy_0_xx_1, t_0_yyy_0_xy_0, t_0_yyy_0_xy_1,\
                               t_0_yyy_0_xz_0, t_0_yyy_0_xz_1, t_0_yyy_0_y_1, t_0_yyy_0_yy_0,\
                               t_0_yyy_0_yy_1, t_0_yyy_0_yz_0, t_0_yyy_0_yz_1, t_0_yyy_0_z_1,\
                               t_0_yyy_0_zz_0, t_0_yyy_0_zz_1, t_0_yyz_0_xz_0, t_0_yyz_0_xz_1,\
                               t_0_yyz_0_yy_0, t_0_yyz_0_yy_1, t_0_yyz_0_yz_0, t_0_yyz_0_yz_1,\
                               t_0_yyz_0_z_1, t_0_yyz_0_zz_0, t_0_yyz_0_zz_1, t_0_yzz_0_xy_0,\
                               t_0_yzz_0_xy_1, t_0_yzz_0_y_1, t_0_yzz_0_yy_0, t_0_yzz_0_yy_1,\
                               t_0_yzz_0_yz_0, t_0_yzz_0_yz_1, t_0_yzz_0_zz_0, t_0_yzz_0_zz_1,\
                               t_0_zz_0_xz_0, t_0_zz_0_xz_1, t_0_zz_0_yy_0, t_0_zz_0_yy_1,\
                               t_0_zz_0_yz_0, t_0_zz_0_yz_1, t_0_zz_0_zz_0, t_0_zz_0_zz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xyzz_0_zz_0[i] = rpb_x[i] * t_0_yzz_0_zz_0[i] + rwp_x[i] * t_0_yzz_0_zz_1[i];

            t_0_xyzz_0_yz_0[i] = rpb_x[i] * t_0_yzz_0_yz_0[i] + rwp_x[i] * t_0_yzz_0_yz_1[i];

            t_0_xyzz_0_yy_0[i] = rpb_x[i] * t_0_yzz_0_yy_0[i] + rwp_x[i] * t_0_yzz_0_yy_1[i];

            t_0_xyzz_0_xz_0[i] = rpb_y[i] * t_0_xzz_0_xz_0[i] + rwp_y[i] * t_0_xzz_0_xz_1[i];

            t_0_xyzz_0_xy_0[i] = rpb_x[i] * t_0_yzz_0_xy_0[i] + rwp_x[i] * t_0_yzz_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yzz_0_y_1[i];

            t_0_xyzz_0_xx_0[i] = rpb_y[i] * t_0_xzz_0_xx_0[i] + rwp_y[i] * t_0_xzz_0_xx_1[i];

            t_0_xyyz_0_zz_0[i] = rpb_x[i] * t_0_yyz_0_zz_0[i] + rwp_x[i] * t_0_yyz_0_zz_1[i];

            t_0_xyyz_0_yz_0[i] = rpb_x[i] * t_0_yyz_0_yz_0[i] + rwp_x[i] * t_0_yyz_0_yz_1[i];

            t_0_xyyz_0_yy_0[i] = rpb_x[i] * t_0_yyz_0_yy_0[i] + rwp_x[i] * t_0_yyz_0_yy_1[i];

            t_0_xyyz_0_xz_0[i] = rpb_x[i] * t_0_yyz_0_xz_0[i] + rwp_x[i] * t_0_yyz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyz_0_z_1[i];

            t_0_xyyz_0_xy_0[i] = rpb_z[i] * t_0_xyy_0_xy_0[i] + rwp_z[i] * t_0_xyy_0_xy_1[i];

            t_0_xyyz_0_xx_0[i] = rpb_z[i] * t_0_xyy_0_xx_0[i] + rwp_z[i] * t_0_xyy_0_xx_1[i];

            t_0_xyyy_0_zz_0[i] = rpb_x[i] * t_0_yyy_0_zz_0[i] + rwp_x[i] * t_0_yyy_0_zz_1[i];

            t_0_xyyy_0_yz_0[i] = rpb_x[i] * t_0_yyy_0_yz_0[i] + rwp_x[i] * t_0_yyy_0_yz_1[i];

            t_0_xyyy_0_yy_0[i] = rpb_x[i] * t_0_yyy_0_yy_0[i] + rwp_x[i] * t_0_yyy_0_yy_1[i];

            t_0_xyyy_0_xz_0[i] = rpb_x[i] * t_0_yyy_0_xz_0[i] + rwp_x[i] * t_0_yyy_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_z_1[i];

            t_0_xyyy_0_xy_0[i] = rpb_x[i] * t_0_yyy_0_xy_0[i] + rwp_x[i] * t_0_yyy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_yyy_0_y_1[i];

            t_0_xyyy_0_xx_0[i] = rpb_x[i] * t_0_yyy_0_xx_0[i] + rwp_x[i] * t_0_yyy_0_xx_1[i] + fze_0[i] * t_0_yyy_0_x_1[i];

            t_0_xxzz_0_zz_0[i] = rpb_x[i] * t_0_xzz_0_zz_0[i] + rwp_x[i] * t_0_xzz_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_zz_1[i];

            t_0_xxzz_0_yz_0[i] = rpb_x[i] * t_0_xzz_0_yz_0[i] + rwp_x[i] * t_0_xzz_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yz_1[i];

            t_0_xxzz_0_yy_0[i] = rpb_x[i] * t_0_xzz_0_yy_0[i] + rwp_x[i] * t_0_xzz_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_yy_1[i];

            t_0_xxzz_0_xz_0[i] = rpb_x[i] * t_0_xzz_0_xz_0[i] + rwp_x[i] * t_0_xzz_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xzz_0_z_1[i];

            t_0_xxzz_0_xy_0[i] = rpb_z[i] * t_0_xxz_0_xy_0[i] + rwp_z[i] * t_0_xxz_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxzz_0_xx_0[i] = rpb_z[i] * t_0_xxz_0_xx_0[i] + rwp_z[i] * t_0_xxz_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxyz_0_zz_0[i] = rpb_y[i] * t_0_xxz_0_zz_0[i] + rwp_y[i] * t_0_xxz_0_zz_1[i];

            t_0_xxyz_0_yz_0[i] = rpb_y[i] * t_0_xxz_0_yz_0[i] + rwp_y[i] * t_0_xxz_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxz_0_z_1[i];

            t_0_xxyz_0_yy_0[i] = rpb_z[i] * t_0_xxy_0_yy_0[i] + rwp_z[i] * t_0_xxy_0_yy_1[i];

            t_0_xxyz_0_xz_0[i] = rpb_y[i] * t_0_xxz_0_xz_0[i] + rwp_y[i] * t_0_xxz_0_xz_1[i];

            t_0_xxyz_0_xy_0[i] = rpb_z[i] * t_0_xxy_0_xy_0[i] + rwp_z[i] * t_0_xxy_0_xy_1[i];

            t_0_xxyz_0_xx_0[i] = rpb_y[i] * t_0_xxz_0_xx_0[i] + rwp_y[i] * t_0_xxz_0_xx_1[i];

            t_0_xxyy_0_zz_0[i] = rpb_x[i] * t_0_xyy_0_zz_0[i] + rwp_x[i] * t_0_xyy_0_zz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_zz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_zz_1[i];

            t_0_xxyy_0_yz_0[i] = rpb_x[i] * t_0_xyy_0_yz_0[i] + rwp_x[i] * t_0_xyy_0_yz_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yz_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yz_1[i];

            t_0_xxyy_0_yy_0[i] = rpb_x[i] * t_0_xyy_0_yy_0[i] + rwp_x[i] * t_0_xyy_0_yy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_yy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_yy_1[i];

            t_0_xxyy_0_xz_0[i] = rpb_y[i] * t_0_xxy_0_xz_0[i] + rwp_y[i] * t_0_xxy_0_xz_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xz_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxyy_0_xy_0[i] = rpb_x[i] * t_0_xyy_0_xy_0[i] + rwp_x[i] * t_0_xyy_0_xy_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_xy_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xyy_0_y_1[i];

            t_0_xxyy_0_xx_0[i] = rpb_y[i] * t_0_xxy_0_xx_0[i] + rwp_y[i] * t_0_xxy_0_xx_1[i] + fact_1_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_1_2 * frz2_0[i] * t_0_xx_0_xx_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_xx_0, t_0_xx_0_xx_1, t_0_xx_0_xy_0, t_0_xx_0_xy_1,\
                               t_0_xx_0_xz_0, t_0_xx_0_xz_1, t_0_xx_0_yy_0, t_0_xx_0_yy_1,\
                               t_0_xx_0_yz_0, t_0_xx_0_yz_1, t_0_xx_0_zz_0, t_0_xx_0_zz_1,\
                               t_0_xxx_0_x_1, t_0_xxx_0_xx_0, t_0_xxx_0_xx_1, t_0_xxx_0_xy_0,\
                               t_0_xxx_0_xy_1, t_0_xxx_0_xz_0, t_0_xxx_0_xz_1, t_0_xxx_0_y_1,\
                               t_0_xxx_0_yy_0, t_0_xxx_0_yy_1, t_0_xxx_0_yz_0, t_0_xxx_0_yz_1,\
                               t_0_xxx_0_z_1, t_0_xxx_0_zz_0, t_0_xxx_0_zz_1, t_0_xxxx_0_xx_0,\
                               t_0_xxxx_0_xy_0, t_0_xxxx_0_xz_0, t_0_xxxx_0_yy_0, t_0_xxxx_0_yz_0,\
                               t_0_xxxx_0_zz_0, t_0_xxxy_0_xx_0, t_0_xxxy_0_xy_0, t_0_xxxy_0_xz_0,\
                               t_0_xxxy_0_yy_0, t_0_xxxy_0_yz_0, t_0_xxxy_0_zz_0, t_0_xxxz_0_xx_0,\
                               t_0_xxxz_0_xy_0, t_0_xxxz_0_xz_0, t_0_xxxz_0_yy_0, t_0_xxxz_0_yz_0,\
                               t_0_xxxz_0_zz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxxz_0_zz_0[i] = rpb_z[i] * t_0_xxx_0_zz_0[i] + rwp_z[i] * t_0_xxx_0_zz_1[i] + fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxz_0_yz_0[i] = rpb_z[i] * t_0_xxx_0_yz_0[i] + rwp_z[i] * t_0_xxx_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxz_0_yy_0[i] = rpb_z[i] * t_0_xxx_0_yy_0[i] + rwp_z[i] * t_0_xxx_0_yy_1[i];

            t_0_xxxz_0_xz_0[i] = rpb_z[i] * t_0_xxx_0_xz_0[i] + rwp_z[i] * t_0_xxx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxz_0_xy_0[i] = rpb_z[i] * t_0_xxx_0_xy_0[i] + rwp_z[i] * t_0_xxx_0_xy_1[i];

            t_0_xxxz_0_xx_0[i] = rpb_z[i] * t_0_xxx_0_xx_0[i] + rwp_z[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxy_0_zz_0[i] = rpb_y[i] * t_0_xxx_0_zz_0[i] + rwp_y[i] * t_0_xxx_0_zz_1[i];

            t_0_xxxy_0_yz_0[i] = rpb_y[i] * t_0_xxx_0_yz_0[i] + rwp_y[i] * t_0_xxx_0_yz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxy_0_yy_0[i] = rpb_y[i] * t_0_xxx_0_yy_0[i] + rwp_y[i] * t_0_xxx_0_yy_1[i] + fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxy_0_xz_0[i] = rpb_y[i] * t_0_xxx_0_xz_0[i] + rwp_y[i] * t_0_xxx_0_xz_1[i];

            t_0_xxxy_0_xy_0[i] = rpb_y[i] * t_0_xxx_0_xy_0[i] + rwp_y[i] * t_0_xxx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_x_1[i];

            t_0_xxxy_0_xx_0[i] = rpb_y[i] * t_0_xxx_0_xx_0[i] + rwp_y[i] * t_0_xxx_0_xx_1[i];

            t_0_xxxx_0_zz_0[i] = rpb_x[i] * t_0_xxx_0_zz_0[i] + rwp_x[i] * t_0_xxx_0_zz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_zz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxxx_0_yz_0[i] = rpb_x[i] * t_0_xxx_0_yz_0[i] + rwp_x[i] * t_0_xxx_0_yz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxxx_0_yy_0[i] = rpb_x[i] * t_0_xxx_0_yy_0[i] + rwp_x[i] * t_0_xxx_0_yy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_yy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxxx_0_xz_0[i] = rpb_x[i] * t_0_xxx_0_xz_0[i] + rwp_x[i] * t_0_xxx_0_xz_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xz_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xz_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_z_1[i];

            t_0_xxxx_0_xy_0[i] = rpb_x[i] * t_0_xxx_0_xy_0[i] + rwp_x[i] * t_0_xxx_0_xy_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xy_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xy_1[i] + fact_1_2 * fze_0[i] * t_0_xxx_0_y_1[i];

            t_0_xxxx_0_xx_0[i] = rpb_x[i] * t_0_xxx_0_xx_0[i] + rwp_x[i] * t_0_xxx_0_xx_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_xx_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_xx_1[i] + fze_0[i] * t_0_xxx_0_x_1[i];
        }
    }
}


} // derirec namespace
