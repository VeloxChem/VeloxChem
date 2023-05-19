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
compHostVRRForSDSG_V0(      BufferHostXY<T>&      intsBufferSDSG,
                      const BufferHostX<int32_t>& intsIndexesSDSG0,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
                      const BufferHostXY<T>&      intsBufferSPSF1,
                      const BufferHostX<int32_t>& intsIndexesSPSF1,
                      const BufferHostXY<T>&      intsBufferSPSG0,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSPSG1,
                      const BufferHostX<int32_t>& intsIndexesSPSG1,
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

    t_0_yz_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(15));

    t_0_yz_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(16));

    t_0_yz_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(17));

    t_0_yz_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(18));

    t_0_yz_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(19));

    t_0_yz_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(20));

    t_0_yz_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(21));

    t_0_yz_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(22));

    t_0_yz_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(23));

    t_0_yz_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(24));

    t_0_yz_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(25));

    t_0_yz_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(26));

    t_0_yz_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(27));

    t_0_yz_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(28));

    t_0_yz_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(29));

    t_0_yy_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(30));

    t_0_yy_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(31));

    t_0_yy_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(32));

    t_0_yy_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(33));

    t_0_yy_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(34));

    t_0_yy_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(35));

    t_0_yy_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(36));

    t_0_yy_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(37));

    t_0_yy_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(38));

    t_0_yy_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(39));

    t_0_yy_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(40));

    t_0_yy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(41));

    t_0_yy_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(42));

    t_0_yy_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(43));

    t_0_yy_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(44));

    t_0_xz_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(45));

    t_0_xz_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(46));

    t_0_xz_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(47));

    t_0_xz_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(48));

    t_0_xz_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(49));

    t_0_xz_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(50));

    t_0_xz_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(51));

    t_0_xz_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(52));

    t_0_xz_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(53));

    t_0_xz_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(54));

    t_0_xz_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(55));

    t_0_xz_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(56));

    t_0_xz_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(57));

    t_0_xz_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(58));

    t_0_xz_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(59));

    t_0_xy_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(60));

    t_0_xy_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(61));

    t_0_xy_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(62));

    t_0_xy_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(63));

    t_0_xy_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(64));

    t_0_xy_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(65));

    t_0_xy_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(66));

    t_0_xy_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(67));

    t_0_xy_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(68));

    t_0_xy_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(69));

    t_0_xy_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(70));

    t_0_xy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(71));

    t_0_xy_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(72));

    t_0_xy_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(73));

    t_0_xy_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(74));

    t_0_xx_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(75));

    t_0_xx_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(76));

    t_0_xx_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(77));

    t_0_xx_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(78));

    t_0_xx_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(79));

    t_0_xx_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(80));

    t_0_xx_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(81));

    t_0_xx_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(82));

    t_0_xx_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(83));

    t_0_xx_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(84));

    t_0_xx_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(85));

    t_0_xx_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(86));

    t_0_xx_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(87));

    t_0_xx_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(88));

    t_0_xx_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(89));

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up [SPSF]^(1) integral components

    t_0_z_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(0));

    t_0_z_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(1));

    t_0_z_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(2));

    t_0_z_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(3));

    t_0_z_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(4));

    t_0_z_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(5));

    t_0_z_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(6));

    t_0_z_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(7));

    t_0_z_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(8));

    t_0_z_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(9));

    t_0_y_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(10));

    t_0_y_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(11));

    t_0_y_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(12));

    t_0_y_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(13));

    t_0_y_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(14));

    t_0_y_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(15));

    t_0_y_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(16));

    t_0_y_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(17));

    t_0_y_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(18));

    t_0_y_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(19));

    t_0_x_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(20));

    t_0_x_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(21));

    t_0_x_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(22));

    t_0_x_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(23));

    t_0_x_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(24));

    t_0_x_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(25));

    t_0_x_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(26));

    t_0_x_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(27));

    t_0_x_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(28));

    t_0_x_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(29));

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_z_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_z_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_z_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_z_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_z_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_z_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_z_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

    t_0_y_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(15));

    t_0_y_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(16));

    t_0_y_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(17));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(18));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(19));

    t_0_y_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(20));

    t_0_y_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(21));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(22));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(23));

    t_0_y_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(24));

    t_0_y_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(25));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(26));

    t_0_y_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(27));

    t_0_y_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(28));

    t_0_y_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(29));

    t_0_x_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(30));

    t_0_x_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(31));

    t_0_x_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(32));

    t_0_x_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(33));

    t_0_x_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(34));

    t_0_x_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(35));

    t_0_x_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(36));

    t_0_x_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(37));

    t_0_x_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(38));

    t_0_x_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(39));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(40));

    t_0_x_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(41));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(42));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(43));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(44));

    // set up [SPSG]^(1) integral components

    t_0_z_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(0));

    t_0_z_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(1));

    t_0_z_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(2));

    t_0_z_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(3));

    t_0_z_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(4));

    t_0_z_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(5));

    t_0_z_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(6));

    t_0_z_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(7));

    t_0_z_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(8));

    t_0_z_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(9));

    t_0_z_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(10));

    t_0_z_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(11));

    t_0_z_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(12));

    t_0_z_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(13));

    t_0_z_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(14));

    t_0_y_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(15));

    t_0_y_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(16));

    t_0_y_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(17));

    t_0_y_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(18));

    t_0_y_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(19));

    t_0_y_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(20));

    t_0_y_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(21));

    t_0_y_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(22));

    t_0_y_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(23));

    t_0_y_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(24));

    t_0_y_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(25));

    t_0_y_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(26));

    t_0_y_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(27));

    t_0_y_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(28));

    t_0_y_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(29));

    t_0_x_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(30));

    t_0_x_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(31));

    t_0_x_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(32));

    t_0_x_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(33));

    t_0_x_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(34));

    t_0_x_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(35));

    t_0_x_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(36));

    t_0_x_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(37));

    t_0_x_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(38));

    t_0_x_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(39));

    t_0_x_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(40));

    t_0_x_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(41));

    t_0_x_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(42));

    t_0_x_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(43));

    t_0_x_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(44));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxyy_0,\
                               t_0_y_0_xxyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1, t_0_y_0_xzzz_0,\
                               t_0_y_0_xzzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0, t_0_y_0_yyyy_1,\
                               t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1, t_0_y_0_yyzz_0,\
                               t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0, t_0_y_0_yzzz_1,\
                               t_0_y_0_zzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_yz_0_xxxx_0, t_0_yz_0_xxxy_0, t_0_yz_0_xxxz_0,\
                               t_0_yz_0_xxyy_0, t_0_yz_0_xxyz_0, t_0_yz_0_xxzz_0, t_0_yz_0_xyyy_0,\
                               t_0_yz_0_xyyz_0, t_0_yz_0_xyzz_0, t_0_yz_0_xzzz_0, t_0_yz_0_yyyy_0,\
                               t_0_yz_0_yyyz_0, t_0_yz_0_yyzz_0, t_0_yz_0_yzzz_0, t_0_yz_0_zzzz_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1, t_0_z_0_xxxy_0,\
                               t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1, t_0_z_0_xxy_1,\
                               t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0, t_0_z_0_xxyz_1,\
                               t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1, t_0_z_0_xyy_1,\
                               t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0, t_0_z_0_xyyz_1,\
                               t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1, t_0_z_0_xzz_1,\
                               t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1, t_0_z_0_yyyy_0,\
                               t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1, t_0_z_0_yyz_1,\
                               t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1, t_0_z_0_yzzz_0,\
                               t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0, t_0_z_0_zzzz_1,\
                               t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0, t_0_zz_0_xxyy_0,\
                               t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0, t_0_zz_0_xyyz_0,\
                               t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0, t_0_zz_0_yyyz_0,\
                               t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] += rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] += rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] += rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] += rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] += rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] += rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] += rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] += rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] += rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] += rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] += rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] += rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] += rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] += rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] += rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yz_0_zzzz_0[i] += rpb_y[i] * t_0_z_0_zzzz_0[i] + rwp_y[i] * t_0_z_0_zzzz_1[i];

            t_0_yz_0_yzzz_0[i] += rpb_y[i] * t_0_z_0_yzzz_0[i] + rwp_y[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_yz_0_yyzz_0[i] += rpb_y[i] * t_0_z_0_yyzz_0[i] + rwp_y[i] * t_0_z_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_yz_0_yyyz_0[i] += rpb_y[i] * t_0_z_0_yyyz_0[i] + rwp_y[i] * t_0_z_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_yz_0_yyyy_0[i] += rpb_z[i] * t_0_y_0_yyyy_0[i] + rwp_z[i] * t_0_y_0_yyyy_1[i];

            t_0_yz_0_xzzz_0[i] += rpb_y[i] * t_0_z_0_xzzz_0[i] + rwp_y[i] * t_0_z_0_xzzz_1[i];

            t_0_yz_0_xyzz_0[i] += rpb_y[i] * t_0_z_0_xyzz_0[i] + rwp_y[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_yz_0_xyyz_0[i] += rpb_y[i] * t_0_z_0_xyyz_0[i] + rwp_y[i] * t_0_z_0_xyyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_yz_0_xyyy_0[i] += rpb_z[i] * t_0_y_0_xyyy_0[i] + rwp_z[i] * t_0_y_0_xyyy_1[i];

            t_0_yz_0_xxzz_0[i] += rpb_y[i] * t_0_z_0_xxzz_0[i] + rwp_y[i] * t_0_z_0_xxzz_1[i];

            t_0_yz_0_xxyz_0[i] += rpb_y[i] * t_0_z_0_xxyz_0[i] + rwp_y[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_yz_0_xxyy_0[i] += rpb_z[i] * t_0_y_0_xxyy_0[i] + rwp_z[i] * t_0_y_0_xxyy_1[i];

            t_0_yz_0_xxxz_0[i] += rpb_y[i] * t_0_z_0_xxxz_0[i] + rwp_y[i] * t_0_z_0_xxxz_1[i];

            t_0_yz_0_xxxy_0[i] += rpb_z[i] * t_0_y_0_xxxy_0[i] + rwp_z[i] * t_0_y_0_xxxy_1[i];

            t_0_yz_0_xxxx_0[i] += rpb_y[i] * t_0_z_0_xxxx_0[i] + rwp_y[i] * t_0_z_0_xxxx_1[i];

            t_0_yy_0_zzzz_0[i] += rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] += rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] += rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] += rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] += rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] += rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxx_1, t_0_x_0_xxxy_0, t_0_x_0_xxxy_1, t_0_x_0_xxyy_0,\
                               t_0_x_0_xxyy_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyyy_0,\
                               t_0_x_0_xyyy_1, t_0_x_0_xzzz_0, t_0_x_0_xzzz_1, t_0_xy_0_xxyy_0,\
                               t_0_xy_0_xxyz_0, t_0_xy_0_xxzz_0, t_0_xy_0_xyyy_0, t_0_xy_0_xyyz_0,\
                               t_0_xy_0_xyzz_0, t_0_xy_0_xzzz_0, t_0_xy_0_yyyy_0, t_0_xy_0_yyyz_0,\
                               t_0_xy_0_yyzz_0, t_0_xy_0_yzzz_0, t_0_xy_0_zzzz_0, t_0_xz_0_xxxx_0,\
                               t_0_xz_0_xxxy_0, t_0_xz_0_xxxz_0, t_0_xz_0_xxyy_0, t_0_xz_0_xxyz_0,\
                               t_0_xz_0_xxzz_0, t_0_xz_0_xyyy_0, t_0_xz_0_xyyz_0, t_0_xz_0_xyzz_0,\
                               t_0_xz_0_xzzz_0, t_0_xz_0_yyyy_0, t_0_xz_0_yyyz_0, t_0_xz_0_yyzz_0,\
                               t_0_xz_0_yzzz_0, t_0_xz_0_zzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0,\
                               t_0_y_0_xxxx_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxxz_0,\
                               t_0_y_0_xxxz_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0, t_0_y_0_xxyy_1,\
                               t_0_y_0_xxyz_0, t_0_y_0_xxyz_1, t_0_y_0_xxz_1, t_0_y_0_xxzz_0,\
                               t_0_y_0_xxzz_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_xyyz_0, t_0_y_0_xyyz_1, t_0_y_0_xyz_1, t_0_y_0_xyzz_0,\
                               t_0_y_0_xyzz_1, t_0_y_0_xzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyy_1, t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1,\
                               t_0_y_0_yyzz_0, t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0,\
                               t_0_y_0_yzzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xxxx_0,\
                               t_0_yy_0_xxxy_0, t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0,\
                               t_0_yy_0_xxzz_0, t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0,\
                               t_0_z_0_xxxz_0, t_0_z_0_xxxz_1, t_0_z_0_xxyz_0, t_0_z_0_xxyz_1,\
                               t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyyy_0,\
                               t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1, t_0_z_0_yyz_1,\
                               t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1, t_0_z_0_yzzz_0,\
                               t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0, t_0_z_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yy_0_xyzz_0[i] += rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] += rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] += rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] += rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] += rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] += rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] += rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];

            t_0_yy_0_xxxy_0[i] += rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] += rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xz_0_zzzz_0[i] += rpb_x[i] * t_0_z_0_zzzz_0[i] + rwp_x[i] * t_0_z_0_zzzz_1[i];

            t_0_xz_0_yzzz_0[i] += rpb_x[i] * t_0_z_0_yzzz_0[i] + rwp_x[i] * t_0_z_0_yzzz_1[i];

            t_0_xz_0_yyzz_0[i] += rpb_x[i] * t_0_z_0_yyzz_0[i] + rwp_x[i] * t_0_z_0_yyzz_1[i];

            t_0_xz_0_yyyz_0[i] += rpb_x[i] * t_0_z_0_yyyz_0[i] + rwp_x[i] * t_0_z_0_yyyz_1[i];

            t_0_xz_0_yyyy_0[i] += rpb_x[i] * t_0_z_0_yyyy_0[i] + rwp_x[i] * t_0_z_0_yyyy_1[i];

            t_0_xz_0_xzzz_0[i] += rpb_x[i] * t_0_z_0_xzzz_0[i] + rwp_x[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_xz_0_xyzz_0[i] += rpb_x[i] * t_0_z_0_xyzz_0[i] + rwp_x[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_xz_0_xyyz_0[i] += rpb_x[i] * t_0_z_0_xyyz_0[i] + rwp_x[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_xz_0_xyyy_0[i] += rpb_z[i] * t_0_x_0_xyyy_0[i] + rwp_z[i] * t_0_x_0_xyyy_1[i];

            t_0_xz_0_xxzz_0[i] += rpb_x[i] * t_0_z_0_xxzz_0[i] + rwp_x[i] * t_0_z_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_xz_0_xxyz_0[i] += rpb_x[i] * t_0_z_0_xxyz_0[i] + rwp_x[i] * t_0_z_0_xxyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_xz_0_xxyy_0[i] += rpb_z[i] * t_0_x_0_xxyy_0[i] + rwp_z[i] * t_0_x_0_xxyy_1[i];

            t_0_xz_0_xxxz_0[i] += rpb_x[i] * t_0_z_0_xxxz_0[i] + rwp_x[i] * t_0_z_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_xz_0_xxxy_0[i] += rpb_z[i] * t_0_x_0_xxxy_0[i] + rwp_z[i] * t_0_x_0_xxxy_1[i];

            t_0_xz_0_xxxx_0[i] += rpb_z[i] * t_0_x_0_xxxx_0[i] + rwp_z[i] * t_0_x_0_xxxx_1[i];

            t_0_xy_0_zzzz_0[i] += rpb_x[i] * t_0_y_0_zzzz_0[i] + rwp_x[i] * t_0_y_0_zzzz_1[i];

            t_0_xy_0_yzzz_0[i] += rpb_x[i] * t_0_y_0_yzzz_0[i] + rwp_x[i] * t_0_y_0_yzzz_1[i];

            t_0_xy_0_yyzz_0[i] += rpb_x[i] * t_0_y_0_yyzz_0[i] + rwp_x[i] * t_0_y_0_yyzz_1[i];

            t_0_xy_0_yyyz_0[i] += rpb_x[i] * t_0_y_0_yyyz_0[i] + rwp_x[i] * t_0_y_0_yyyz_1[i];

            t_0_xy_0_yyyy_0[i] += rpb_x[i] * t_0_y_0_yyyy_0[i] + rwp_x[i] * t_0_y_0_yyyy_1[i];

            t_0_xy_0_xzzz_0[i] += rpb_y[i] * t_0_x_0_xzzz_0[i] + rwp_y[i] * t_0_x_0_xzzz_1[i];

            t_0_xy_0_xyzz_0[i] += rpb_x[i] * t_0_y_0_xyzz_0[i] + rwp_x[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_xy_0_xyyz_0[i] += rpb_x[i] * t_0_y_0_xyyz_0[i] + rwp_x[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_xy_0_xyyy_0[i] += rpb_x[i] * t_0_y_0_xyyy_0[i] + rwp_x[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_xy_0_xxzz_0[i] += rpb_y[i] * t_0_x_0_xxzz_0[i] + rwp_y[i] * t_0_x_0_xxzz_1[i];

            t_0_xy_0_xxyz_0[i] += rpb_x[i] * t_0_y_0_xxyz_0[i] + rwp_x[i] * t_0_y_0_xxyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_xy_0_xxyy_0[i] += rpb_x[i] * t_0_y_0_xxyy_0[i] + rwp_x[i] * t_0_y_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xyy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rwp_x, rwp_y, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_x_0_xxx_1, t_0_x_0_xxxx_0, t_0_x_0_xxxx_1,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxy_1, t_0_x_0_xxxz_0, t_0_x_0_xxxz_1,\
                               t_0_x_0_xxy_1, t_0_x_0_xxyy_0, t_0_x_0_xxyy_1, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxyz_1, t_0_x_0_xxz_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1,\
                               t_0_x_0_xyy_1, t_0_x_0_xyyy_0, t_0_x_0_xyyy_1, t_0_x_0_xyyz_0,\
                               t_0_x_0_xyyz_1, t_0_x_0_xyz_1, t_0_x_0_xyzz_0, t_0_x_0_xyzz_1,\
                               t_0_x_0_xzz_1, t_0_x_0_xzzz_0, t_0_x_0_xzzz_1, t_0_x_0_yyy_1,\
                               t_0_x_0_yyyy_0, t_0_x_0_yyyy_1, t_0_x_0_yyyz_0, t_0_x_0_yyyz_1,\
                               t_0_x_0_yyz_1, t_0_x_0_yyzz_0, t_0_x_0_yyzz_1, t_0_x_0_yzz_1,\
                               t_0_x_0_yzzz_0, t_0_x_0_yzzz_1, t_0_x_0_zzz_1, t_0_x_0_zzzz_0,\
                               t_0_x_0_zzzz_1, t_0_xx_0_xxxx_0, t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0,\
                               t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0, t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0,\
                               t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0, t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0,\
                               t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0, t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0,\
                               t_0_xy_0_xxxx_0, t_0_xy_0_xxxy_0, t_0_xy_0_xxxz_0, t_0_y_0_xxxy_0,\
                               t_0_y_0_xxxy_1, t_0_y_0_xxy_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xy_0_xxxz_0[i] += rpb_y[i] * t_0_x_0_xxxz_0[i] + rwp_y[i] * t_0_x_0_xxxz_1[i];

            t_0_xy_0_xxxy_0[i] += rpb_x[i] * t_0_y_0_xxxy_0[i] + rwp_x[i] * t_0_y_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_xy_0_xxxx_0[i] += rpb_y[i] * t_0_x_0_xxxx_0[i] + rwp_y[i] * t_0_x_0_xxxx_1[i];

            t_0_xx_0_zzzz_0[i] += rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] += rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] += rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] += rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] += rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] += rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];

            t_0_xx_0_xyzz_0[i] += rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] += rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] += rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] += rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] += rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] += rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] += rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] += rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] += rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxyy_0,\
                               t_0_y_0_xxyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1, t_0_y_0_xzzz_0,\
                               t_0_y_0_xzzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0, t_0_y_0_yyyy_1,\
                               t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1, t_0_y_0_yyzz_0,\
                               t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0, t_0_y_0_yzzz_1,\
                               t_0_y_0_zzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_yz_0_xxxx_0, t_0_yz_0_xxxy_0, t_0_yz_0_xxxz_0,\
                               t_0_yz_0_xxyy_0, t_0_yz_0_xxyz_0, t_0_yz_0_xxzz_0, t_0_yz_0_xyyy_0,\
                               t_0_yz_0_xyyz_0, t_0_yz_0_xyzz_0, t_0_yz_0_xzzz_0, t_0_yz_0_yyyy_0,\
                               t_0_yz_0_yyyz_0, t_0_yz_0_yyzz_0, t_0_yz_0_yzzz_0, t_0_yz_0_zzzz_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1, t_0_z_0_xxxy_0,\
                               t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1, t_0_z_0_xxy_1,\
                               t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0, t_0_z_0_xxyz_1,\
                               t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1, t_0_z_0_xyy_1,\
                               t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0, t_0_z_0_xyyz_1,\
                               t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1, t_0_z_0_xzz_1,\
                               t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1, t_0_z_0_yyyy_0,\
                               t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1, t_0_z_0_yyz_1,\
                               t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1, t_0_z_0_yzzz_0,\
                               t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0, t_0_z_0_zzzz_1,\
                               t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0, t_0_zz_0_xxyy_0,\
                               t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0, t_0_zz_0_xyyz_0,\
                               t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0, t_0_zz_0_yyyz_0,\
                               t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] = rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] = rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] = rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] = rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] = rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] = rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] = rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] = rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] = rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] = rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] = rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] = rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] = rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] = rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] = rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yz_0_zzzz_0[i] = rpb_y[i] * t_0_z_0_zzzz_0[i] + rwp_y[i] * t_0_z_0_zzzz_1[i];

            t_0_yz_0_yzzz_0[i] = rpb_y[i] * t_0_z_0_yzzz_0[i] + rwp_y[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_yz_0_yyzz_0[i] = rpb_y[i] * t_0_z_0_yyzz_0[i] + rwp_y[i] * t_0_z_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_yz_0_yyyz_0[i] = rpb_y[i] * t_0_z_0_yyyz_0[i] + rwp_y[i] * t_0_z_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_yz_0_yyyy_0[i] = rpb_z[i] * t_0_y_0_yyyy_0[i] + rwp_z[i] * t_0_y_0_yyyy_1[i];

            t_0_yz_0_xzzz_0[i] = rpb_y[i] * t_0_z_0_xzzz_0[i] + rwp_y[i] * t_0_z_0_xzzz_1[i];

            t_0_yz_0_xyzz_0[i] = rpb_y[i] * t_0_z_0_xyzz_0[i] + rwp_y[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_yz_0_xyyz_0[i] = rpb_y[i] * t_0_z_0_xyyz_0[i] + rwp_y[i] * t_0_z_0_xyyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_yz_0_xyyy_0[i] = rpb_z[i] * t_0_y_0_xyyy_0[i] + rwp_z[i] * t_0_y_0_xyyy_1[i];

            t_0_yz_0_xxzz_0[i] = rpb_y[i] * t_0_z_0_xxzz_0[i] + rwp_y[i] * t_0_z_0_xxzz_1[i];

            t_0_yz_0_xxyz_0[i] = rpb_y[i] * t_0_z_0_xxyz_0[i] + rwp_y[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_yz_0_xxyy_0[i] = rpb_z[i] * t_0_y_0_xxyy_0[i] + rwp_z[i] * t_0_y_0_xxyy_1[i];

            t_0_yz_0_xxxz_0[i] = rpb_y[i] * t_0_z_0_xxxz_0[i] + rwp_y[i] * t_0_z_0_xxxz_1[i];

            t_0_yz_0_xxxy_0[i] = rpb_z[i] * t_0_y_0_xxxy_0[i] + rwp_z[i] * t_0_y_0_xxxy_1[i];

            t_0_yz_0_xxxx_0[i] = rpb_y[i] * t_0_z_0_xxxx_0[i] + rwp_y[i] * t_0_z_0_xxxx_1[i];

            t_0_yy_0_zzzz_0[i] = rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] = rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] = rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] = rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] = rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] = rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_x_0_xxxx_0,\
                               t_0_x_0_xxxx_1, t_0_x_0_xxxy_0, t_0_x_0_xxxy_1, t_0_x_0_xxyy_0,\
                               t_0_x_0_xxyy_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyyy_0,\
                               t_0_x_0_xyyy_1, t_0_x_0_xzzz_0, t_0_x_0_xzzz_1, t_0_xy_0_xxyy_0,\
                               t_0_xy_0_xxyz_0, t_0_xy_0_xxzz_0, t_0_xy_0_xyyy_0, t_0_xy_0_xyyz_0,\
                               t_0_xy_0_xyzz_0, t_0_xy_0_xzzz_0, t_0_xy_0_yyyy_0, t_0_xy_0_yyyz_0,\
                               t_0_xy_0_yyzz_0, t_0_xy_0_yzzz_0, t_0_xy_0_zzzz_0, t_0_xz_0_xxxx_0,\
                               t_0_xz_0_xxxy_0, t_0_xz_0_xxxz_0, t_0_xz_0_xxyy_0, t_0_xz_0_xxyz_0,\
                               t_0_xz_0_xxzz_0, t_0_xz_0_xyyy_0, t_0_xz_0_xyyz_0, t_0_xz_0_xyzz_0,\
                               t_0_xz_0_xzzz_0, t_0_xz_0_yyyy_0, t_0_xz_0_yyyz_0, t_0_xz_0_yyzz_0,\
                               t_0_xz_0_yzzz_0, t_0_xz_0_zzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0,\
                               t_0_y_0_xxxx_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxxz_0,\
                               t_0_y_0_xxxz_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0, t_0_y_0_xxyy_1,\
                               t_0_y_0_xxyz_0, t_0_y_0_xxyz_1, t_0_y_0_xxz_1, t_0_y_0_xxzz_0,\
                               t_0_y_0_xxzz_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_xyyz_0, t_0_y_0_xyyz_1, t_0_y_0_xyz_1, t_0_y_0_xyzz_0,\
                               t_0_y_0_xyzz_1, t_0_y_0_xzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyy_1, t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1,\
                               t_0_y_0_yyzz_0, t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0,\
                               t_0_y_0_yzzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xxxx_0,\
                               t_0_yy_0_xxxy_0, t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0,\
                               t_0_yy_0_xxzz_0, t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0,\
                               t_0_z_0_xxxz_0, t_0_z_0_xxxz_1, t_0_z_0_xxyz_0, t_0_z_0_xxyz_1,\
                               t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyyy_0,\
                               t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1, t_0_z_0_yyz_1,\
                               t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1, t_0_z_0_yzzz_0,\
                               t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0, t_0_z_0_zzzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yy_0_xyzz_0[i] = rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] = rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] = rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] = rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] = rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] = rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] = rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];

            t_0_yy_0_xxxy_0[i] = rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] = rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xz_0_zzzz_0[i] = rpb_x[i] * t_0_z_0_zzzz_0[i] + rwp_x[i] * t_0_z_0_zzzz_1[i];

            t_0_xz_0_yzzz_0[i] = rpb_x[i] * t_0_z_0_yzzz_0[i] + rwp_x[i] * t_0_z_0_yzzz_1[i];

            t_0_xz_0_yyzz_0[i] = rpb_x[i] * t_0_z_0_yyzz_0[i] + rwp_x[i] * t_0_z_0_yyzz_1[i];

            t_0_xz_0_yyyz_0[i] = rpb_x[i] * t_0_z_0_yyyz_0[i] + rwp_x[i] * t_0_z_0_yyyz_1[i];

            t_0_xz_0_yyyy_0[i] = rpb_x[i] * t_0_z_0_yyyy_0[i] + rwp_x[i] * t_0_z_0_yyyy_1[i];

            t_0_xz_0_xzzz_0[i] = rpb_x[i] * t_0_z_0_xzzz_0[i] + rwp_x[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_xz_0_xyzz_0[i] = rpb_x[i] * t_0_z_0_xyzz_0[i] + rwp_x[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_xz_0_xyyz_0[i] = rpb_x[i] * t_0_z_0_xyyz_0[i] + rwp_x[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_xz_0_xyyy_0[i] = rpb_z[i] * t_0_x_0_xyyy_0[i] + rwp_z[i] * t_0_x_0_xyyy_1[i];

            t_0_xz_0_xxzz_0[i] = rpb_x[i] * t_0_z_0_xxzz_0[i] + rwp_x[i] * t_0_z_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_xz_0_xxyz_0[i] = rpb_x[i] * t_0_z_0_xxyz_0[i] + rwp_x[i] * t_0_z_0_xxyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_xz_0_xxyy_0[i] = rpb_z[i] * t_0_x_0_xxyy_0[i] + rwp_z[i] * t_0_x_0_xxyy_1[i];

            t_0_xz_0_xxxz_0[i] = rpb_x[i] * t_0_z_0_xxxz_0[i] + rwp_x[i] * t_0_z_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_xz_0_xxxy_0[i] = rpb_z[i] * t_0_x_0_xxxy_0[i] + rwp_z[i] * t_0_x_0_xxxy_1[i];

            t_0_xz_0_xxxx_0[i] = rpb_z[i] * t_0_x_0_xxxx_0[i] + rwp_z[i] * t_0_x_0_xxxx_1[i];

            t_0_xy_0_zzzz_0[i] = rpb_x[i] * t_0_y_0_zzzz_0[i] + rwp_x[i] * t_0_y_0_zzzz_1[i];

            t_0_xy_0_yzzz_0[i] = rpb_x[i] * t_0_y_0_yzzz_0[i] + rwp_x[i] * t_0_y_0_yzzz_1[i];

            t_0_xy_0_yyzz_0[i] = rpb_x[i] * t_0_y_0_yyzz_0[i] + rwp_x[i] * t_0_y_0_yyzz_1[i];

            t_0_xy_0_yyyz_0[i] = rpb_x[i] * t_0_y_0_yyyz_0[i] + rwp_x[i] * t_0_y_0_yyyz_1[i];

            t_0_xy_0_yyyy_0[i] = rpb_x[i] * t_0_y_0_yyyy_0[i] + rwp_x[i] * t_0_y_0_yyyy_1[i];

            t_0_xy_0_xzzz_0[i] = rpb_y[i] * t_0_x_0_xzzz_0[i] + rwp_y[i] * t_0_x_0_xzzz_1[i];

            t_0_xy_0_xyzz_0[i] = rpb_x[i] * t_0_y_0_xyzz_0[i] + rwp_x[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_xy_0_xyyz_0[i] = rpb_x[i] * t_0_y_0_xyyz_0[i] + rwp_x[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_xy_0_xyyy_0[i] = rpb_x[i] * t_0_y_0_xyyy_0[i] + rwp_x[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_xy_0_xxzz_0[i] = rpb_y[i] * t_0_x_0_xxzz_0[i] + rwp_y[i] * t_0_x_0_xxzz_1[i];

            t_0_xy_0_xxyz_0[i] = rpb_x[i] * t_0_y_0_xxyz_0[i] + rwp_x[i] * t_0_y_0_xxyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_xy_0_xxyy_0[i] = rpb_x[i] * t_0_y_0_xxyy_0[i] + rwp_x[i] * t_0_y_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xyy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rwp_x, rwp_y, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_x_0_xxx_1, t_0_x_0_xxxx_0, t_0_x_0_xxxx_1,\
                               t_0_x_0_xxxy_0, t_0_x_0_xxxy_1, t_0_x_0_xxxz_0, t_0_x_0_xxxz_1,\
                               t_0_x_0_xxy_1, t_0_x_0_xxyy_0, t_0_x_0_xxyy_1, t_0_x_0_xxyz_0,\
                               t_0_x_0_xxyz_1, t_0_x_0_xxz_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1,\
                               t_0_x_0_xyy_1, t_0_x_0_xyyy_0, t_0_x_0_xyyy_1, t_0_x_0_xyyz_0,\
                               t_0_x_0_xyyz_1, t_0_x_0_xyz_1, t_0_x_0_xyzz_0, t_0_x_0_xyzz_1,\
                               t_0_x_0_xzz_1, t_0_x_0_xzzz_0, t_0_x_0_xzzz_1, t_0_x_0_yyy_1,\
                               t_0_x_0_yyyy_0, t_0_x_0_yyyy_1, t_0_x_0_yyyz_0, t_0_x_0_yyyz_1,\
                               t_0_x_0_yyz_1, t_0_x_0_yyzz_0, t_0_x_0_yyzz_1, t_0_x_0_yzz_1,\
                               t_0_x_0_yzzz_0, t_0_x_0_yzzz_1, t_0_x_0_zzz_1, t_0_x_0_zzzz_0,\
                               t_0_x_0_zzzz_1, t_0_xx_0_xxxx_0, t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0,\
                               t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0, t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0,\
                               t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0, t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0,\
                               t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0, t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0,\
                               t_0_xy_0_xxxx_0, t_0_xy_0_xxxy_0, t_0_xy_0_xxxz_0, t_0_y_0_xxxy_0,\
                               t_0_y_0_xxxy_1, t_0_y_0_xxy_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xy_0_xxxz_0[i] = rpb_y[i] * t_0_x_0_xxxz_0[i] + rwp_y[i] * t_0_x_0_xxxz_1[i];

            t_0_xy_0_xxxy_0[i] = rpb_x[i] * t_0_y_0_xxxy_0[i] + rwp_x[i] * t_0_y_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_xy_0_xxxx_0[i] = rpb_y[i] * t_0_x_0_xxxx_0[i] + rwp_y[i] * t_0_x_0_xxxx_1[i];

            t_0_xx_0_zzzz_0[i] = rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] = rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] = rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] = rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] = rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] = rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];

            t_0_xx_0_xyzz_0[i] = rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] = rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] = rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] = rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] = rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] = rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] = rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] = rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] = rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSDSG_V1(      BufferHostXY<T>&      intsBufferSDSG,
                      const BufferHostX<int32_t>& intsIndexesSDSG0,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
                      const BufferHostXY<T>&      intsBufferSPSF1,
                      const BufferHostX<int32_t>& intsIndexesSPSF1,
                      const BufferHostXY<T>&      intsBufferSPSG0,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSPSG1,
                      const BufferHostX<int32_t>& intsIndexesSPSG1,
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

    t_0_yz_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(15));

    t_0_yz_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(16));

    t_0_yz_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(17));

    t_0_yz_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(18));

    t_0_yz_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(19));

    t_0_yz_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(20));

    t_0_yz_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(21));

    t_0_yz_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(22));

    t_0_yy_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(23));

    t_0_yy_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(24));

    t_0_yy_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(25));

    t_0_yy_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(26));

    t_0_yy_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(27));

    t_0_yy_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(28));

    t_0_yy_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(29));

    t_0_yy_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(30));

    t_0_yy_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(31));

    t_0_yy_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(32));

    t_0_yy_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(33));

    t_0_yy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(34));

    t_0_yy_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(35));

    t_0_yy_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(36));

    t_0_yy_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(37));

    t_0_xz_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(38));

    t_0_xz_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(39));

    t_0_xz_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(40));

    t_0_xz_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(41));

    t_0_xy_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(42));

    t_0_xy_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(43));

    t_0_xy_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(44));

    t_0_xx_0_zzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(45));

    t_0_xx_0_yzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(46));

    t_0_xx_0_yyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(47));

    t_0_xx_0_yyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(48));

    t_0_xx_0_yyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(49));

    t_0_xx_0_xzzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(50));

    t_0_xx_0_xyzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(51));

    t_0_xx_0_xyyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(52));

    t_0_xx_0_xyyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(53));

    t_0_xx_0_xxzz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(54));

    t_0_xx_0_xxyz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(55));

    t_0_xx_0_xxyy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(56));

    t_0_xx_0_xxxz_0 = intsBufferSDSG0.data(intsIndexesSDSG0(57));

    t_0_xx_0_xxxy_0 = intsBufferSDSG0.data(intsIndexesSDSG0(58));

    t_0_xx_0_xxxx_0 = intsBufferSDSG0.data(intsIndexesSDSG0(59));

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up [SPSF]^(1) integral components

    t_0_z_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(0));

    t_0_z_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(1));

    t_0_z_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(2));

    t_0_z_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(3));

    t_0_z_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(4));

    t_0_z_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(5));

    t_0_z_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(6));

    t_0_z_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(7));

    t_0_z_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(8));

    t_0_z_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(9));

    t_0_y_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(10));

    t_0_y_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(11));

    t_0_y_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(12));

    t_0_y_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(13));

    t_0_y_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(14));

    t_0_y_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(15));

    t_0_y_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(16));

    t_0_y_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(17));

    t_0_y_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(18));

    t_0_y_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(19));

    t_0_x_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(20));

    t_0_x_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(21));

    t_0_x_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(22));

    t_0_x_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(23));

    t_0_x_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(24));

    t_0_x_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(25));

    t_0_x_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(26));

    t_0_x_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(27));

    t_0_x_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(28));

    t_0_x_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(29));

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_z_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_z_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_z_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_z_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_z_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_z_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_z_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

    t_0_y_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(15));

    t_0_y_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(16));

    t_0_y_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(17));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(18));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(19));

    t_0_y_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(20));

    t_0_y_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(21));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(22));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(23));

    t_0_y_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(24));

    t_0_y_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(25));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(26));

    t_0_y_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(27));

    t_0_y_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(28));

    t_0_y_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(29));

    t_0_x_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(30));

    t_0_x_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(31));

    t_0_x_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(32));

    t_0_x_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(33));

    t_0_x_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(34));

    t_0_x_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(35));

    t_0_x_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(36));

    t_0_x_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(37));

    t_0_x_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(38));

    t_0_x_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(39));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(40));

    t_0_x_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(41));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(42));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(43));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(44));

    // set up [SPSG]^(1) integral components

    t_0_z_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(0));

    t_0_z_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(1));

    t_0_z_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(2));

    t_0_z_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(3));

    t_0_z_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(4));

    t_0_z_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(5));

    t_0_z_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(6));

    t_0_z_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(7));

    t_0_z_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(8));

    t_0_z_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(9));

    t_0_z_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(10));

    t_0_z_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(11));

    t_0_z_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(12));

    t_0_z_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(13));

    t_0_z_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(14));

    t_0_y_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(15));

    t_0_y_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(16));

    t_0_y_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(17));

    t_0_y_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(18));

    t_0_y_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(19));

    t_0_y_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(20));

    t_0_y_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(21));

    t_0_y_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(22));

    t_0_y_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(23));

    t_0_y_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(24));

    t_0_y_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(25));

    t_0_y_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(26));

    t_0_y_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(27));

    t_0_y_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(28));

    t_0_y_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(29));

    t_0_x_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(30));

    t_0_x_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(31));

    t_0_x_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(32));

    t_0_x_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(33));

    t_0_x_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(34));

    t_0_x_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(35));

    t_0_x_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(36));

    t_0_x_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(37));

    t_0_x_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(38));

    t_0_x_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(39));

    t_0_x_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(40));

    t_0_x_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(41));

    t_0_x_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(42));

    t_0_x_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(43));

    t_0_x_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(44));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_y_0_xxxz_0, t_0_y_0_xxxz_1, t_0_y_0_xxy_1,\
                               t_0_y_0_xxyy_0, t_0_y_0_xxyy_1, t_0_y_0_xxyz_0, t_0_y_0_xxyz_1,\
                               t_0_y_0_xxz_1, t_0_y_0_xxzz_0, t_0_y_0_xxzz_1, t_0_y_0_xyy_1,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyy_1, t_0_y_0_xyyz_0, t_0_y_0_xyyz_1,\
                               t_0_y_0_xyz_1, t_0_y_0_xyzz_0, t_0_y_0_xyzz_1, t_0_y_0_xzz_1,\
                               t_0_y_0_xzzz_0, t_0_y_0_xzzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyy_1, t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1,\
                               t_0_y_0_yyzz_0, t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0,\
                               t_0_y_0_yzzz_1, t_0_y_0_zzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1,\
                               t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0, t_0_yy_0_xxzz_0,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_yz_0_xxyz_0, t_0_yz_0_xyyz_0, t_0_yz_0_xyzz_0,\
                               t_0_yz_0_yyyy_0, t_0_yz_0_yyyz_0, t_0_yz_0_yyzz_0, t_0_yz_0_yzzz_0,\
                               t_0_yz_0_zzzz_0, t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1,\
                               t_0_z_0_xxxy_0, t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1,\
                               t_0_z_0_xxy_1, t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxyz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xyy_1, t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1,\
                               t_0_z_0_yyyy_0, t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1,\
                               t_0_z_0_yyz_1, t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1,\
                               t_0_z_0_yzzz_0, t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0,\
                               t_0_z_0_zzzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] += rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] += rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] += rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] += rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] += rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] += rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] += rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] += rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] += rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] += rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] += rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] += rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] += rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] += rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] += rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yz_0_zzzz_0[i] += rpb_y[i] * t_0_z_0_zzzz_0[i] + rwp_y[i] * t_0_z_0_zzzz_1[i];

            t_0_yz_0_yzzz_0[i] += rpb_y[i] * t_0_z_0_yzzz_0[i] + rwp_y[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_yz_0_yyzz_0[i] += rpb_y[i] * t_0_z_0_yyzz_0[i] + rwp_y[i] * t_0_z_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_yz_0_yyyz_0[i] += rpb_y[i] * t_0_z_0_yyyz_0[i] + rwp_y[i] * t_0_z_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_yz_0_yyyy_0[i] += rpb_z[i] * t_0_y_0_yyyy_0[i] + rwp_z[i] * t_0_y_0_yyyy_1[i];

            t_0_yz_0_xyzz_0[i] += rpb_y[i] * t_0_z_0_xyzz_0[i] + rwp_y[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_yz_0_xyyz_0[i] += rpb_y[i] * t_0_z_0_xyyz_0[i] + rwp_y[i] * t_0_z_0_xyyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_yz_0_xxyz_0[i] += rpb_y[i] * t_0_z_0_xxyz_0[i] + rwp_y[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_yy_0_zzzz_0[i] += rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] += rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] += rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] += rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] += rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] += rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];

            t_0_yy_0_xyzz_0[i] += rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] += rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] += rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] += rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] += rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] += rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] += rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1, t_0_0_0_yyyz_0,\
                               t_0_0_0_yyyz_1, t_0_0_0_yyzz_0, t_0_0_0_yyzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_yzzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxx_1,\
                               t_0_x_0_xxxx_0, t_0_x_0_xxxx_1, t_0_x_0_xxxy_0, t_0_x_0_xxxy_1,\
                               t_0_x_0_xxxz_0, t_0_x_0_xxxz_1, t_0_x_0_xxy_1, t_0_x_0_xxyy_0,\
                               t_0_x_0_xxyy_1, t_0_x_0_xxyz_0, t_0_x_0_xxyz_1, t_0_x_0_xxz_1,\
                               t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyy_1, t_0_x_0_xyyy_0,\
                               t_0_x_0_xyyy_1, t_0_x_0_xyyz_0, t_0_x_0_xyyz_1, t_0_x_0_xyz_1,\
                               t_0_x_0_xyzz_0, t_0_x_0_xyzz_1, t_0_x_0_xzz_1, t_0_x_0_xzzz_0,\
                               t_0_x_0_xzzz_1, t_0_x_0_yyy_1, t_0_x_0_yyyy_0, t_0_x_0_yyyy_1,\
                               t_0_x_0_yyyz_0, t_0_x_0_yyyz_1, t_0_x_0_yyz_1, t_0_x_0_yyzz_0,\
                               t_0_x_0_yyzz_1, t_0_x_0_yzz_1, t_0_x_0_yzzz_0, t_0_x_0_yzzz_1,\
                               t_0_x_0_zzz_1, t_0_x_0_zzzz_0, t_0_x_0_zzzz_1, t_0_xx_0_xxxx_0,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0, t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0, t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0,\
                               t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0, t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0, t_0_xy_0_xxxy_0, t_0_xy_0_xxyy_0,\
                               t_0_xy_0_xyyy_0, t_0_xz_0_xxxx_0, t_0_xz_0_xxxz_0, t_0_xz_0_xxzz_0,\
                               t_0_xz_0_xzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0, t_0_y_0_xxxx_1,\
                               t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0,\
                               t_0_y_0_xxyy_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_yyy_1, t_0_yy_0_xxxx_0, t_0_yy_0_xxxy_0, t_0_z_0_xxxz_0,\
                               t_0_z_0_xxxz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yy_0_xxxy_0[i] += rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] += rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xz_0_xzzz_0[i] += rpb_x[i] * t_0_z_0_xzzz_0[i] + rwp_x[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_xz_0_xxzz_0[i] += rpb_x[i] * t_0_z_0_xxzz_0[i] + rwp_x[i] * t_0_z_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_xz_0_xxxz_0[i] += rpb_x[i] * t_0_z_0_xxxz_0[i] + rwp_x[i] * t_0_z_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_xz_0_xxxx_0[i] += rpb_z[i] * t_0_x_0_xxxx_0[i] + rwp_z[i] * t_0_x_0_xxxx_1[i];

            t_0_xy_0_xyyy_0[i] += rpb_x[i] * t_0_y_0_xyyy_0[i] + rwp_x[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_xy_0_xxyy_0[i] += rpb_x[i] * t_0_y_0_xxyy_0[i] + rwp_x[i] * t_0_y_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_xy_0_xxxy_0[i] += rpb_x[i] * t_0_y_0_xxxy_0[i] + rwp_x[i] * t_0_y_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_xx_0_zzzz_0[i] += rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] += rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] += rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] += rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] += rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] += rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];

            t_0_xx_0_xyzz_0[i] += rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] += rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] += rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] += rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] += rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] += rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] += rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] += rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] += rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_0_0_xxxx_0,\
                               t_0_0_0_xxxx_1, t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0,\
                               t_0_0_0_xxxz_1, t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0,\
                               t_0_0_0_xxyz_1, t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0,\
                               t_0_0_0_xyyy_1, t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0,\
                               t_0_0_0_xyzz_1, t_0_0_0_xzzz_0, t_0_0_0_xzzz_1, t_0_0_0_yyyy_0,\
                               t_0_0_0_yyyy_1, t_0_0_0_yyyz_0, t_0_0_0_yyyz_1, t_0_0_0_yyzz_0,\
                               t_0_0_0_yyzz_1, t_0_0_0_yzzz_0, t_0_0_0_yzzz_1, t_0_0_0_zzzz_0,\
                               t_0_0_0_zzzz_1, t_0_y_0_xxxz_0, t_0_y_0_xxxz_1, t_0_y_0_xxy_1,\
                               t_0_y_0_xxyy_0, t_0_y_0_xxyy_1, t_0_y_0_xxyz_0, t_0_y_0_xxyz_1,\
                               t_0_y_0_xxz_1, t_0_y_0_xxzz_0, t_0_y_0_xxzz_1, t_0_y_0_xyy_1,\
                               t_0_y_0_xyyy_0, t_0_y_0_xyyy_1, t_0_y_0_xyyz_0, t_0_y_0_xyyz_1,\
                               t_0_y_0_xyz_1, t_0_y_0_xyzz_0, t_0_y_0_xyzz_1, t_0_y_0_xzz_1,\
                               t_0_y_0_xzzz_0, t_0_y_0_xzzz_1, t_0_y_0_yyy_1, t_0_y_0_yyyy_0,\
                               t_0_y_0_yyyy_1, t_0_y_0_yyyz_0, t_0_y_0_yyyz_1, t_0_y_0_yyz_1,\
                               t_0_y_0_yyzz_0, t_0_y_0_yyzz_1, t_0_y_0_yzz_1, t_0_y_0_yzzz_0,\
                               t_0_y_0_yzzz_1, t_0_y_0_zzz_1, t_0_y_0_zzzz_0, t_0_y_0_zzzz_1,\
                               t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0, t_0_yy_0_xxzz_0,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_yz_0_xxyz_0, t_0_yz_0_xyyz_0, t_0_yz_0_xyzz_0,\
                               t_0_yz_0_yyyy_0, t_0_yz_0_yyyz_0, t_0_yz_0_yyzz_0, t_0_yz_0_yzzz_0,\
                               t_0_yz_0_zzzz_0, t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1,\
                               t_0_z_0_xxxy_0, t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1,\
                               t_0_z_0_xxy_1, t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxyz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xyy_1, t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1,\
                               t_0_z_0_yyyy_0, t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1,\
                               t_0_z_0_yyz_1, t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1,\
                               t_0_z_0_yzzz_0, t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0,\
                               t_0_z_0_zzzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] = rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] = rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] = rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] = rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] = rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] = rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] = rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] = rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] = rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] = rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] = rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] = rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] = rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] = rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] = rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yz_0_zzzz_0[i] = rpb_y[i] * t_0_z_0_zzzz_0[i] + rwp_y[i] * t_0_z_0_zzzz_1[i];

            t_0_yz_0_yzzz_0[i] = rpb_y[i] * t_0_z_0_yzzz_0[i] + rwp_y[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_yz_0_yyzz_0[i] = rpb_y[i] * t_0_z_0_yyzz_0[i] + rwp_y[i] * t_0_z_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_yz_0_yyyz_0[i] = rpb_y[i] * t_0_z_0_yyyz_0[i] + rwp_y[i] * t_0_z_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_yz_0_yyyy_0[i] = rpb_z[i] * t_0_y_0_yyyy_0[i] + rwp_z[i] * t_0_y_0_yyyy_1[i];

            t_0_yz_0_xyzz_0[i] = rpb_y[i] * t_0_z_0_xyzz_0[i] + rwp_y[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_yz_0_xyyz_0[i] = rpb_y[i] * t_0_z_0_xyyz_0[i] + rwp_y[i] * t_0_z_0_xyyz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_yz_0_xxyz_0[i] = rpb_y[i] * t_0_z_0_xxyz_0[i] + rwp_y[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_yy_0_zzzz_0[i] = rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] = rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] = rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] = rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] = rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] = rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];

            t_0_yy_0_xyzz_0[i] = rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] = rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] = rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] = rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] = rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] = rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] = rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1, t_0_0_0_yyyz_0,\
                               t_0_0_0_yyyz_1, t_0_0_0_yyzz_0, t_0_0_0_yyzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_yzzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xxx_1,\
                               t_0_x_0_xxxx_0, t_0_x_0_xxxx_1, t_0_x_0_xxxy_0, t_0_x_0_xxxy_1,\
                               t_0_x_0_xxxz_0, t_0_x_0_xxxz_1, t_0_x_0_xxy_1, t_0_x_0_xxyy_0,\
                               t_0_x_0_xxyy_1, t_0_x_0_xxyz_0, t_0_x_0_xxyz_1, t_0_x_0_xxz_1,\
                               t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyy_1, t_0_x_0_xyyy_0,\
                               t_0_x_0_xyyy_1, t_0_x_0_xyyz_0, t_0_x_0_xyyz_1, t_0_x_0_xyz_1,\
                               t_0_x_0_xyzz_0, t_0_x_0_xyzz_1, t_0_x_0_xzz_1, t_0_x_0_xzzz_0,\
                               t_0_x_0_xzzz_1, t_0_x_0_yyy_1, t_0_x_0_yyyy_0, t_0_x_0_yyyy_1,\
                               t_0_x_0_yyyz_0, t_0_x_0_yyyz_1, t_0_x_0_yyz_1, t_0_x_0_yyzz_0,\
                               t_0_x_0_yyzz_1, t_0_x_0_yzz_1, t_0_x_0_yzzz_0, t_0_x_0_yzzz_1,\
                               t_0_x_0_zzz_1, t_0_x_0_zzzz_0, t_0_x_0_zzzz_1, t_0_xx_0_xxxx_0,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0, t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0, t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0,\
                               t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0, t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0, t_0_xy_0_xxxy_0, t_0_xy_0_xxyy_0,\
                               t_0_xy_0_xyyy_0, t_0_xz_0_xxxx_0, t_0_xz_0_xxxz_0, t_0_xz_0_xxzz_0,\
                               t_0_xz_0_xzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0, t_0_y_0_xxxx_1,\
                               t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0,\
                               t_0_y_0_xxyy_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_yyy_1, t_0_yy_0_xxxx_0, t_0_yy_0_xxxy_0, t_0_z_0_xxxz_0,\
                               t_0_z_0_xxxz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yy_0_xxxy_0[i] = rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] = rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xz_0_xzzz_0[i] = rpb_x[i] * t_0_z_0_xzzz_0[i] + rwp_x[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_xz_0_xxzz_0[i] = rpb_x[i] * t_0_z_0_xxzz_0[i] + rwp_x[i] * t_0_z_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_xz_0_xxxz_0[i] = rpb_x[i] * t_0_z_0_xxxz_0[i] + rwp_x[i] * t_0_z_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_xz_0_xxxx_0[i] = rpb_z[i] * t_0_x_0_xxxx_0[i] + rwp_z[i] * t_0_x_0_xxxx_1[i];

            t_0_xy_0_xyyy_0[i] = rpb_x[i] * t_0_y_0_xyyy_0[i] + rwp_x[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_xy_0_xxyy_0[i] = rpb_x[i] * t_0_y_0_xxyy_0[i] + rwp_x[i] * t_0_y_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_xy_0_xxxy_0[i] = rpb_x[i] * t_0_y_0_xxxy_0[i] + rwp_x[i] * t_0_y_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_xx_0_zzzz_0[i] = rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] = rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] = rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] = rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] = rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] = rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];

            t_0_xx_0_xyzz_0[i] = rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] = rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] = rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] = rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] = rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] = rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] = rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] = rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] = rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSDSG_V2(      BufferHostXY<T>&      intsBufferSDSG,
                      const BufferHostX<int32_t>& intsIndexesSDSG0,
                      const BufferHostXY<T>&      intsBufferSSSG0,
                      const BufferHostX<int32_t>& intsIndexesSSSG0,
                      const BufferHostXY<T>&      intsBufferSSSG1,
                      const BufferHostX<int32_t>& intsIndexesSSSG1,
                      const BufferHostXY<T>&      intsBufferSPSF1,
                      const BufferHostX<int32_t>& intsIndexesSPSF1,
                      const BufferHostXY<T>&      intsBufferSPSG0,
                      const BufferHostX<int32_t>& intsIndexesSPSG0,
                      const BufferHostXY<T>&      intsBufferSPSG1,
                      const BufferHostX<int32_t>& intsIndexesSPSG1,
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

    // set up [SSSG]^(0) integral components

    t_0_0_0_zzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(0));

    t_0_0_0_yzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(1));

    t_0_0_0_yyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(2));

    t_0_0_0_yyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(3));

    t_0_0_0_yyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(4));

    t_0_0_0_xzzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(5));

    t_0_0_0_xyzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(6));

    t_0_0_0_xyyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(7));

    t_0_0_0_xyyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(8));

    t_0_0_0_xxzz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(9));

    t_0_0_0_xxyz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(10));

    t_0_0_0_xxyy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(11));

    t_0_0_0_xxxz_0 = intsBufferSSSG0.data(intsIndexesSSSG0(12));

    t_0_0_0_xxxy_0 = intsBufferSSSG0.data(intsIndexesSSSG0(13));

    t_0_0_0_xxxx_0 = intsBufferSSSG0.data(intsIndexesSSSG0(14));

    // set up [SSSG]^(1) integral components

    t_0_0_0_zzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(0));

    t_0_0_0_yzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(1));

    t_0_0_0_yyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(2));

    t_0_0_0_yyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(3));

    t_0_0_0_yyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(4));

    t_0_0_0_xzzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(5));

    t_0_0_0_xyzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(6));

    t_0_0_0_xyyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(7));

    t_0_0_0_xyyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(8));

    t_0_0_0_xxzz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(9));

    t_0_0_0_xxyz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(10));

    t_0_0_0_xxyy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(11));

    t_0_0_0_xxxz_1 = intsBufferSSSG1.data(intsIndexesSSSG1(12));

    t_0_0_0_xxxy_1 = intsBufferSSSG1.data(intsIndexesSSSG1(13));

    t_0_0_0_xxxx_1 = intsBufferSSSG1.data(intsIndexesSSSG1(14));

    // set up [SPSF]^(1) integral components

    t_0_z_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(0));

    t_0_z_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(1));

    t_0_z_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(2));

    t_0_z_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(3));

    t_0_z_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(4));

    t_0_z_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(5));

    t_0_z_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(6));

    t_0_z_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(7));

    t_0_z_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(8));

    t_0_z_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(9));

    t_0_y_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(10));

    t_0_y_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(11));

    t_0_y_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(12));

    t_0_y_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(13));

    t_0_y_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(14));

    t_0_y_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(15));

    t_0_y_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(16));

    t_0_y_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(17));

    t_0_y_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(18));

    t_0_y_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(19));

    t_0_x_0_zzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(20));

    t_0_x_0_yzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(21));

    t_0_x_0_yyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(22));

    t_0_x_0_yyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(23));

    t_0_x_0_xzz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(24));

    t_0_x_0_xyz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(25));

    t_0_x_0_xyy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(26));

    t_0_x_0_xxz_1 = intsBufferSPSF1.data(intsIndexesSPSF1(27));

    t_0_x_0_xxy_1 = intsBufferSPSF1.data(intsIndexesSPSF1(28));

    t_0_x_0_xxx_1 = intsBufferSPSF1.data(intsIndexesSPSF1(29));

    // set up [SPSG]^(0) integral components

    t_0_z_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(0));

    t_0_z_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(1));

    t_0_z_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(2));

    t_0_z_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(3));

    t_0_z_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(4));

    t_0_z_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(5));

    t_0_z_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(6));

    t_0_z_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(7));

    t_0_z_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(8));

    t_0_z_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(9));

    t_0_z_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(10));

    t_0_z_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(11));

    t_0_z_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(12));

    t_0_z_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(13));

    t_0_z_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(14));

    t_0_y_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(15));

    t_0_y_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(16));

    t_0_y_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(17));

    t_0_y_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(18));

    t_0_y_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(19));

    t_0_y_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(20));

    t_0_y_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(21));

    t_0_y_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(22));

    t_0_y_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(23));

    t_0_y_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(24));

    t_0_y_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(25));

    t_0_y_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(26));

    t_0_y_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(27));

    t_0_y_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(28));

    t_0_y_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(29));

    t_0_x_0_zzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(30));

    t_0_x_0_yzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(31));

    t_0_x_0_yyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(32));

    t_0_x_0_yyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(33));

    t_0_x_0_yyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(34));

    t_0_x_0_xzzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(35));

    t_0_x_0_xyzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(36));

    t_0_x_0_xyyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(37));

    t_0_x_0_xyyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(38));

    t_0_x_0_xxzz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(39));

    t_0_x_0_xxyz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(40));

    t_0_x_0_xxyy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(41));

    t_0_x_0_xxxz_0 = intsBufferSPSG0.data(intsIndexesSPSG0(42));

    t_0_x_0_xxxy_0 = intsBufferSPSG0.data(intsIndexesSPSG0(43));

    t_0_x_0_xxxx_0 = intsBufferSPSG0.data(intsIndexesSPSG0(44));

    // set up [SPSG]^(1) integral components

    t_0_z_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(0));

    t_0_z_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(1));

    t_0_z_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(2));

    t_0_z_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(3));

    t_0_z_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(4));

    t_0_z_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(5));

    t_0_z_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(6));

    t_0_z_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(7));

    t_0_z_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(8));

    t_0_z_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(9));

    t_0_z_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(10));

    t_0_z_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(11));

    t_0_z_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(12));

    t_0_z_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(13));

    t_0_z_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(14));

    t_0_y_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(15));

    t_0_y_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(16));

    t_0_y_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(17));

    t_0_y_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(18));

    t_0_y_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(19));

    t_0_y_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(20));

    t_0_y_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(21));

    t_0_y_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(22));

    t_0_y_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(23));

    t_0_y_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(24));

    t_0_y_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(25));

    t_0_y_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(26));

    t_0_y_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(27));

    t_0_y_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(28));

    t_0_y_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(29));

    t_0_x_0_zzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(30));

    t_0_x_0_yzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(31));

    t_0_x_0_yyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(32));

    t_0_x_0_yyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(33));

    t_0_x_0_yyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(34));

    t_0_x_0_xzzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(35));

    t_0_x_0_xyzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(36));

    t_0_x_0_xyyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(37));

    t_0_x_0_xyyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(38));

    t_0_x_0_xxzz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(39));

    t_0_x_0_xxyz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(40));

    t_0_x_0_xxyy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(41));

    t_0_x_0_xxxz_1 = intsBufferSPSG1.data(intsIndexesSPSG1(42));

    t_0_x_0_xxxy_1 = intsBufferSPSG1.data(intsIndexesSPSG1(43));

    t_0_x_0_xxxx_1 = intsBufferSPSG1.data(intsIndexesSPSG1(44));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    const auto fact_2 = static_cast<T>(2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1, t_0_0_0_yyyz_0,\
                               t_0_0_0_yyyz_1, t_0_0_0_yyzz_0, t_0_0_0_yyzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_yzzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xzzz_0,\
                               t_0_x_0_xzzz_1, t_0_x_0_yyyy_0, t_0_x_0_yyyy_1, t_0_x_0_yyyz_0,\
                               t_0_x_0_yyyz_1, t_0_x_0_yyzz_0, t_0_x_0_yyzz_1, t_0_x_0_yzzz_0,\
                               t_0_x_0_yzzz_1, t_0_x_0_zzz_1, t_0_x_0_zzzz_0, t_0_x_0_zzzz_1,\
                               t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0, t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0,\
                               t_0_y_0_xxxx_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxxz_0,\
                               t_0_y_0_xxxz_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0, t_0_y_0_xxyy_1,\
                               t_0_y_0_xxyz_0, t_0_y_0_xxyz_1, t_0_y_0_xxz_1, t_0_y_0_xxzz_0,\
                               t_0_y_0_xxzz_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_xyyz_0, t_0_y_0_xyyz_1, t_0_y_0_xyz_1, t_0_y_0_xyzz_0,\
                               t_0_y_0_xyzz_1, t_0_y_0_xzz_1, t_0_y_0_xzzz_0, t_0_y_0_xzzz_1,\
                               t_0_y_0_yyy_1, t_0_y_0_yyyy_0, t_0_y_0_yyyy_1, t_0_y_0_yyyz_0,\
                               t_0_y_0_yyyz_1, t_0_y_0_yyz_1, t_0_y_0_yyzz_0, t_0_y_0_yyzz_1,\
                               t_0_y_0_yzz_1, t_0_y_0_yzzz_0, t_0_y_0_yzzz_1, t_0_y_0_zzz_1,\
                               t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xxxx_0, t_0_yy_0_xxxy_0,\
                               t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0, t_0_yy_0_xxzz_0,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1,\
                               t_0_z_0_xxxy_0, t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1,\
                               t_0_z_0_xxy_1, t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxyz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xyy_1, t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1,\
                               t_0_z_0_yyyy_0, t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1,\
                               t_0_z_0_yyz_1, t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1,\
                               t_0_z_0_yzzz_0, t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0,\
                               t_0_z_0_zzzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] += rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] += rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] += rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] += rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] += rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] += rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] += rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] += rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] += rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] += rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] += rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] += rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] += rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] += rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] += rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yy_0_zzzz_0[i] += rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] += rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] += rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] += rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] += rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] += rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];

            t_0_yy_0_xyzz_0[i] += rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] += rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] += rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] += rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] += rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] += rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] += rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];

            t_0_yy_0_xxxy_0[i] += rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] += rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xx_0_zzzz_0[i] += rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] += rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] += rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] += rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] += rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] += rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1,\
                               t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1,\
                               t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1,\
                               t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1,\
                               t_0_x_0_xxx_1, t_0_x_0_xxxx_0, t_0_x_0_xxxx_1, t_0_x_0_xxxy_0,\
                               t_0_x_0_xxxy_1, t_0_x_0_xxxz_0, t_0_x_0_xxxz_1, t_0_x_0_xxy_1,\
                               t_0_x_0_xxyy_0, t_0_x_0_xxyy_1, t_0_x_0_xxyz_0, t_0_x_0_xxyz_1,\
                               t_0_x_0_xxz_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyy_1,\
                               t_0_x_0_xyyy_0, t_0_x_0_xyyy_1, t_0_x_0_xyyz_0, t_0_x_0_xyyz_1,\
                               t_0_x_0_xyz_1, t_0_x_0_xyzz_0, t_0_x_0_xyzz_1, t_0_x_0_xzz_1,\
                               t_0_x_0_yyy_1, t_0_x_0_yyz_1, t_0_x_0_yzz_1, t_0_xx_0_xxxx_0,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0, t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0, t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xx_0_xyzz_0[i] += rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] += rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] += rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] += rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] += rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] += rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] += rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] += rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] += rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1, t_0_0_0_xxxy_0,\
                               t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1, t_0_0_0_xxyy_0,\
                               t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1, t_0_0_0_xxzz_0,\
                               t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1, t_0_0_0_xyyz_0,\
                               t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1, t_0_0_0_xzzz_0,\
                               t_0_0_0_xzzz_1, t_0_0_0_yyyy_0, t_0_0_0_yyyy_1, t_0_0_0_yyyz_0,\
                               t_0_0_0_yyyz_1, t_0_0_0_yyzz_0, t_0_0_0_yyzz_1, t_0_0_0_yzzz_0,\
                               t_0_0_0_yzzz_1, t_0_0_0_zzzz_0, t_0_0_0_zzzz_1, t_0_x_0_xzzz_0,\
                               t_0_x_0_xzzz_1, t_0_x_0_yyyy_0, t_0_x_0_yyyy_1, t_0_x_0_yyyz_0,\
                               t_0_x_0_yyyz_1, t_0_x_0_yyzz_0, t_0_x_0_yyzz_1, t_0_x_0_yzzz_0,\
                               t_0_x_0_yzzz_1, t_0_x_0_zzz_1, t_0_x_0_zzzz_0, t_0_x_0_zzzz_1,\
                               t_0_xx_0_xzzz_0, t_0_xx_0_yyyy_0, t_0_xx_0_yyyz_0, t_0_xx_0_yyzz_0,\
                               t_0_xx_0_yzzz_0, t_0_xx_0_zzzz_0, t_0_y_0_xxx_1, t_0_y_0_xxxx_0,\
                               t_0_y_0_xxxx_1, t_0_y_0_xxxy_0, t_0_y_0_xxxy_1, t_0_y_0_xxxz_0,\
                               t_0_y_0_xxxz_1, t_0_y_0_xxy_1, t_0_y_0_xxyy_0, t_0_y_0_xxyy_1,\
                               t_0_y_0_xxyz_0, t_0_y_0_xxyz_1, t_0_y_0_xxz_1, t_0_y_0_xxzz_0,\
                               t_0_y_0_xxzz_1, t_0_y_0_xyy_1, t_0_y_0_xyyy_0, t_0_y_0_xyyy_1,\
                               t_0_y_0_xyyz_0, t_0_y_0_xyyz_1, t_0_y_0_xyz_1, t_0_y_0_xyzz_0,\
                               t_0_y_0_xyzz_1, t_0_y_0_xzz_1, t_0_y_0_xzzz_0, t_0_y_0_xzzz_1,\
                               t_0_y_0_yyy_1, t_0_y_0_yyyy_0, t_0_y_0_yyyy_1, t_0_y_0_yyyz_0,\
                               t_0_y_0_yyyz_1, t_0_y_0_yyz_1, t_0_y_0_yyzz_0, t_0_y_0_yyzz_1,\
                               t_0_y_0_yzz_1, t_0_y_0_yzzz_0, t_0_y_0_yzzz_1, t_0_y_0_zzz_1,\
                               t_0_y_0_zzzz_0, t_0_y_0_zzzz_1, t_0_yy_0_xxxx_0, t_0_yy_0_xxxy_0,\
                               t_0_yy_0_xxxz_0, t_0_yy_0_xxyy_0, t_0_yy_0_xxyz_0, t_0_yy_0_xxzz_0,\
                               t_0_yy_0_xyyy_0, t_0_yy_0_xyyz_0, t_0_yy_0_xyzz_0, t_0_yy_0_xzzz_0,\
                               t_0_yy_0_yyyy_0, t_0_yy_0_yyyz_0, t_0_yy_0_yyzz_0, t_0_yy_0_yzzz_0,\
                               t_0_yy_0_zzzz_0, t_0_z_0_xxx_1, t_0_z_0_xxxx_0, t_0_z_0_xxxx_1,\
                               t_0_z_0_xxxy_0, t_0_z_0_xxxy_1, t_0_z_0_xxxz_0, t_0_z_0_xxxz_1,\
                               t_0_z_0_xxy_1, t_0_z_0_xxyy_0, t_0_z_0_xxyy_1, t_0_z_0_xxyz_0,\
                               t_0_z_0_xxyz_1, t_0_z_0_xxz_1, t_0_z_0_xxzz_0, t_0_z_0_xxzz_1,\
                               t_0_z_0_xyy_1, t_0_z_0_xyyy_0, t_0_z_0_xyyy_1, t_0_z_0_xyyz_0,\
                               t_0_z_0_xyyz_1, t_0_z_0_xyz_1, t_0_z_0_xyzz_0, t_0_z_0_xyzz_1,\
                               t_0_z_0_xzz_1, t_0_z_0_xzzz_0, t_0_z_0_xzzz_1, t_0_z_0_yyy_1,\
                               t_0_z_0_yyyy_0, t_0_z_0_yyyy_1, t_0_z_0_yyyz_0, t_0_z_0_yyyz_1,\
                               t_0_z_0_yyz_1, t_0_z_0_yyzz_0, t_0_z_0_yyzz_1, t_0_z_0_yzz_1,\
                               t_0_z_0_yzzz_0, t_0_z_0_yzzz_1, t_0_z_0_zzz_1, t_0_z_0_zzzz_0,\
                               t_0_z_0_zzzz_1, t_0_zz_0_xxxx_0, t_0_zz_0_xxxy_0, t_0_zz_0_xxxz_0,\
                               t_0_zz_0_xxyy_0, t_0_zz_0_xxyz_0, t_0_zz_0_xxzz_0, t_0_zz_0_xyyy_0,\
                               t_0_zz_0_xyyz_0, t_0_zz_0_xyzz_0, t_0_zz_0_xzzz_0, t_0_zz_0_yyyy_0,\
                               t_0_zz_0_yyyz_0, t_0_zz_0_yyzz_0, t_0_zz_0_yzzz_0, t_0_zz_0_zzzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zz_0_zzzz_0[i] = rpb_z[i] * t_0_z_0_zzzz_0[i] + rwp_z[i] * t_0_z_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i] + fact_2 * fze_0[i] * t_0_z_0_zzz_1[i];

            t_0_zz_0_yzzz_0[i] = rpb_z[i] * t_0_z_0_yzzz_0[i] + rwp_z[i] * t_0_z_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_yzz_1[i];

            t_0_zz_0_yyzz_0[i] = rpb_z[i] * t_0_z_0_yyzz_0[i] + rwp_z[i] * t_0_z_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_z_0_yyz_1[i];

            t_0_zz_0_yyyz_0[i] = rpb_z[i] * t_0_z_0_yyyz_0[i] + rwp_z[i] * t_0_z_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_yyy_1[i];

            t_0_zz_0_yyyy_0[i] = rpb_z[i] * t_0_z_0_yyyy_0[i] + rwp_z[i] * t_0_z_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_zz_0_xzzz_0[i] = rpb_z[i] * t_0_z_0_xzzz_0[i] + rwp_z[i] * t_0_z_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_3_2 * fze_0[i] * t_0_z_0_xzz_1[i];

            t_0_zz_0_xyzz_0[i] = rpb_z[i] * t_0_z_0_xyzz_0[i] + rwp_z[i] * t_0_z_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fze_0[i] * t_0_z_0_xyz_1[i];

            t_0_zz_0_xyyz_0[i] = rpb_z[i] * t_0_z_0_xyyz_0[i] + rwp_z[i] * t_0_z_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xyy_1[i];

            t_0_zz_0_xyyy_0[i] = rpb_z[i] * t_0_z_0_xyyy_0[i] + rwp_z[i] * t_0_z_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i];

            t_0_zz_0_xxzz_0[i] = rpb_z[i] * t_0_z_0_xxzz_0[i] + rwp_z[i] * t_0_z_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_z_0_xxz_1[i];

            t_0_zz_0_xxyz_0[i] = rpb_z[i] * t_0_z_0_xxyz_0[i] + rwp_z[i] * t_0_z_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxy_1[i];

            t_0_zz_0_xxyy_0[i] = rpb_z[i] * t_0_z_0_xxyy_0[i] + rwp_z[i] * t_0_z_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i];

            t_0_zz_0_xxxz_0[i] = rpb_z[i] * t_0_z_0_xxxz_0[i] + rwp_z[i] * t_0_z_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_1_2 * fze_0[i] * t_0_z_0_xxx_1[i];

            t_0_zz_0_xxxy_0[i] = rpb_z[i] * t_0_z_0_xxxy_0[i] + rwp_z[i] * t_0_z_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i];

            t_0_zz_0_xxxx_0[i] = rpb_z[i] * t_0_z_0_xxxx_0[i] + rwp_z[i] * t_0_z_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_yy_0_zzzz_0[i] = rpb_y[i] * t_0_y_0_zzzz_0[i] + rwp_y[i] * t_0_y_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_yy_0_yzzz_0[i] = rpb_y[i] * t_0_y_0_yzzz_0[i] + rwp_y[i] * t_0_y_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_zzz_1[i];

            t_0_yy_0_yyzz_0[i] = rpb_y[i] * t_0_y_0_yyzz_0[i] + rwp_y[i] * t_0_y_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i] + fze_0[i] * t_0_y_0_yzz_1[i];

            t_0_yy_0_yyyz_0[i] = rpb_y[i] * t_0_y_0_yyyz_0[i] + rwp_y[i] * t_0_y_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_yyz_1[i];

            t_0_yy_0_yyyy_0[i] = rpb_y[i] * t_0_y_0_yyyy_0[i] + rwp_y[i] * t_0_y_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i] + fact_2 * fze_0[i] * t_0_y_0_yyy_1[i];

            t_0_yy_0_xzzz_0[i] = rpb_y[i] * t_0_y_0_xzzz_0[i] + rwp_y[i] * t_0_y_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i];

            t_0_yy_0_xyzz_0[i] = rpb_y[i] * t_0_y_0_xyzz_0[i] + rwp_y[i] * t_0_y_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xzz_1[i];

            t_0_yy_0_xyyz_0[i] = rpb_y[i] * t_0_y_0_xyyz_0[i] + rwp_y[i] * t_0_y_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fze_0[i] * t_0_y_0_xyz_1[i];

            t_0_yy_0_xyyy_0[i] = rpb_y[i] * t_0_y_0_xyyy_0[i] + rwp_y[i] * t_0_y_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_3_2 * fze_0[i] * t_0_y_0_xyy_1[i];

            t_0_yy_0_xxzz_0[i] = rpb_y[i] * t_0_y_0_xxzz_0[i] + rwp_y[i] * t_0_y_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i];

            t_0_yy_0_xxyz_0[i] = rpb_y[i] * t_0_y_0_xxyz_0[i] + rwp_y[i] * t_0_y_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxz_1[i];

            t_0_yy_0_xxyy_0[i] = rpb_y[i] * t_0_y_0_xxyy_0[i] + rwp_y[i] * t_0_y_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_y_0_xxy_1[i];

            t_0_yy_0_xxxz_0[i] = rpb_y[i] * t_0_y_0_xxxz_0[i] + rwp_y[i] * t_0_y_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i];

            t_0_yy_0_xxxy_0[i] = rpb_y[i] * t_0_y_0_xxxy_0[i] + rwp_y[i] * t_0_y_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_1_2 * fze_0[i] * t_0_y_0_xxx_1[i];

            t_0_yy_0_xxxx_0[i] = rpb_y[i] * t_0_y_0_xxxx_0[i] + rwp_y[i] * t_0_y_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i];

            t_0_xx_0_zzzz_0[i] = rpb_x[i] * t_0_x_0_zzzz_0[i] + rwp_x[i] * t_0_x_0_zzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_zzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_zzzz_1[i];

            t_0_xx_0_yzzz_0[i] = rpb_x[i] * t_0_x_0_yzzz_0[i] + rwp_x[i] * t_0_x_0_yzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yzzz_1[i];

            t_0_xx_0_yyzz_0[i] = rpb_x[i] * t_0_x_0_yyzz_0[i] + rwp_x[i] * t_0_x_0_yyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyzz_1[i];

            t_0_xx_0_yyyz_0[i] = rpb_x[i] * t_0_x_0_yyyz_0[i] + rwp_x[i] * t_0_x_0_yyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyz_1[i];

            t_0_xx_0_yyyy_0[i] = rpb_x[i] * t_0_x_0_yyyy_0[i] + rwp_x[i] * t_0_x_0_yyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_yyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_yyyy_1[i];

            t_0_xx_0_xzzz_0[i] = rpb_x[i] * t_0_x_0_xzzz_0[i] + rwp_x[i] * t_0_x_0_xzzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xzzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xzzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_zzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_0_0_xxxx_0, t_0_0_0_xxxx_1,\
                               t_0_0_0_xxxy_0, t_0_0_0_xxxy_1, t_0_0_0_xxxz_0, t_0_0_0_xxxz_1,\
                               t_0_0_0_xxyy_0, t_0_0_0_xxyy_1, t_0_0_0_xxyz_0, t_0_0_0_xxyz_1,\
                               t_0_0_0_xxzz_0, t_0_0_0_xxzz_1, t_0_0_0_xyyy_0, t_0_0_0_xyyy_1,\
                               t_0_0_0_xyyz_0, t_0_0_0_xyyz_1, t_0_0_0_xyzz_0, t_0_0_0_xyzz_1,\
                               t_0_x_0_xxx_1, t_0_x_0_xxxx_0, t_0_x_0_xxxx_1, t_0_x_0_xxxy_0,\
                               t_0_x_0_xxxy_1, t_0_x_0_xxxz_0, t_0_x_0_xxxz_1, t_0_x_0_xxy_1,\
                               t_0_x_0_xxyy_0, t_0_x_0_xxyy_1, t_0_x_0_xxyz_0, t_0_x_0_xxyz_1,\
                               t_0_x_0_xxz_1, t_0_x_0_xxzz_0, t_0_x_0_xxzz_1, t_0_x_0_xyy_1,\
                               t_0_x_0_xyyy_0, t_0_x_0_xyyy_1, t_0_x_0_xyyz_0, t_0_x_0_xyyz_1,\
                               t_0_x_0_xyz_1, t_0_x_0_xyzz_0, t_0_x_0_xyzz_1, t_0_x_0_xzz_1,\
                               t_0_x_0_yyy_1, t_0_x_0_yyz_1, t_0_x_0_yzz_1, t_0_xx_0_xxxx_0,\
                               t_0_xx_0_xxxy_0, t_0_xx_0_xxxz_0, t_0_xx_0_xxyy_0, t_0_xx_0_xxyz_0,\
                               t_0_xx_0_xxzz_0, t_0_xx_0_xyyy_0, t_0_xx_0_xyyz_0, t_0_xx_0_xyzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xx_0_xyzz_0[i] = rpb_x[i] * t_0_x_0_xyzz_0[i] + rwp_x[i] * t_0_x_0_xyzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyzz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yzz_1[i];

            t_0_xx_0_xyyz_0[i] = rpb_x[i] * t_0_x_0_xyyz_0[i] + rwp_x[i] * t_0_x_0_xyyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyz_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyz_1[i];

            t_0_xx_0_xyyy_0[i] = rpb_x[i] * t_0_x_0_xyyy_0[i] + rwp_x[i] * t_0_x_0_xyyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xyyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xyyy_1[i] + fact_1_2 * fze_0[i] * t_0_x_0_yyy_1[i];

            t_0_xx_0_xxzz_0[i] = rpb_x[i] * t_0_x_0_xxzz_0[i] + rwp_x[i] * t_0_x_0_xxzz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxzz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxzz_1[i] + fze_0[i] * t_0_x_0_xzz_1[i];

            t_0_xx_0_xxyz_0[i] = rpb_x[i] * t_0_x_0_xxyz_0[i] + rwp_x[i] * t_0_x_0_xxyz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyz_1[i] + fze_0[i] * t_0_x_0_xyz_1[i];

            t_0_xx_0_xxyy_0[i] = rpb_x[i] * t_0_x_0_xxyy_0[i] + rwp_x[i] * t_0_x_0_xxyy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxyy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxyy_1[i] + fze_0[i] * t_0_x_0_xyy_1[i];

            t_0_xx_0_xxxz_0[i] = rpb_x[i] * t_0_x_0_xxxz_0[i] + rwp_x[i] * t_0_x_0_xxxz_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxz_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxz_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxz_1[i];

            t_0_xx_0_xxxy_0[i] = rpb_x[i] * t_0_x_0_xxxy_0[i] + rwp_x[i] * t_0_x_0_xxxy_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxy_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxy_1[i] + fact_3_2 * fze_0[i] * t_0_x_0_xxy_1[i];

            t_0_xx_0_xxxx_0[i] = rpb_x[i] * t_0_x_0_xxxx_0[i] + rwp_x[i] * t_0_x_0_xxxx_1[i] + fact_1_2 * fz_0[i] * t_0_0_0_xxxx_0[i] - fact_1_2 * frz2_0[i] * t_0_0_0_xxxx_1[i] + fact_2 * fze_0[i] * t_0_x_0_xxx_1[i];
        }
    }
}


} // derirec namespace
