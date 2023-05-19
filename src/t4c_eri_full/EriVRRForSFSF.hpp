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
compHostVRRForSFSF_V0(      BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF0,
                      const BufferHostXY<T>&      intsBufferSPSF0,
                      const BufferHostX<int32_t>& intsIndexesSPSF0,
                      const BufferHostXY<T>&      intsBufferSPSF1,
                      const BufferHostX<int32_t>& intsIndexesSPSF1,
                      const BufferHostXY<T>&      intsBufferSDSD1,
                      const BufferHostX<int32_t>& intsIndexesSDSD1,
                      const BufferHostXY<T>&      intsBufferSDSF0,
                      const BufferHostX<int32_t>& intsIndexesSDSF0,
                      const BufferHostXY<T>&      intsBufferSDSF1,
                      const BufferHostX<int32_t>& intsIndexesSDSF1,
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

    t_0_yyz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(29));

    t_0_yyy_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(30));

    t_0_yyy_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(31));

    t_0_yyy_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(32));

    t_0_yyy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(33));

    t_0_yyy_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(34));

    t_0_yyy_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(35));

    t_0_yyy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(36));

    t_0_yyy_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(37));

    t_0_yyy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(38));

    t_0_yyy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(39));

    t_0_xzz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(40));

    t_0_xzz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(41));

    t_0_xzz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(42));

    t_0_xzz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(43));

    t_0_xzz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(44));

    t_0_xzz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(45));

    t_0_xzz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(46));

    t_0_xzz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(47));

    t_0_xzz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(48));

    t_0_xzz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(49));

    t_0_xyz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(50));

    t_0_xyz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(51));

    t_0_xyz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(52));

    t_0_xyz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(53));

    t_0_xyz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(54));

    t_0_xyz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(55));

    t_0_xyz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(56));

    t_0_xyz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(57));

    t_0_xyz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(58));

    t_0_xyz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(59));

    t_0_xyy_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(60));

    t_0_xyy_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(61));

    t_0_xyy_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(62));

    t_0_xyy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(63));

    t_0_xyy_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(64));

    t_0_xyy_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(65));

    t_0_xyy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(66));

    t_0_xyy_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(67));

    t_0_xyy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(68));

    t_0_xyy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(69));

    t_0_xxz_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(70));

    t_0_xxz_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(71));

    t_0_xxz_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(72));

    t_0_xxz_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(73));

    t_0_xxz_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(74));

    t_0_xxz_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(75));

    t_0_xxz_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(76));

    t_0_xxz_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(77));

    t_0_xxz_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(78));

    t_0_xxz_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(79));

    t_0_xxy_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(80));

    t_0_xxy_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(81));

    t_0_xxy_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(82));

    t_0_xxy_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(83));

    t_0_xxy_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(84));

    t_0_xxy_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(85));

    t_0_xxy_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(86));

    t_0_xxy_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(87));

    t_0_xxy_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(88));

    t_0_xxy_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(89));

    t_0_xxx_0_zzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(90));

    t_0_xxx_0_yzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(91));

    t_0_xxx_0_yyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(92));

    t_0_xxx_0_yyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(93));

    t_0_xxx_0_xzz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(94));

    t_0_xxx_0_xyz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(95));

    t_0_xxx_0_xyy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(96));

    t_0_xxx_0_xxz_0 = intsBufferSFSF0.data(intsIndexesSFSF0(97));

    t_0_xxx_0_xxy_0 = intsBufferSFSF0.data(intsIndexesSFSF0(98));

    t_0_xxx_0_xxx_0 = intsBufferSFSF0.data(intsIndexesSFSF0(99));

    // set up [SPSF]^(0) integral components

    t_0_z_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(0));

    t_0_z_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(1));

    t_0_z_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(2));

    t_0_z_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(3));

    t_0_z_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(4));

    t_0_z_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(5));

    t_0_z_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(6));

    t_0_z_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(7));

    t_0_z_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(8));

    t_0_z_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(9));

    t_0_y_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(10));

    t_0_y_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(11));

    t_0_y_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(12));

    t_0_y_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(13));

    t_0_y_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(14));

    t_0_y_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(15));

    t_0_y_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(16));

    t_0_y_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(17));

    t_0_y_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(18));

    t_0_y_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(19));

    t_0_x_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(20));

    t_0_x_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(21));

    t_0_x_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(22));

    t_0_x_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(23));

    t_0_x_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(24));

    t_0_x_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(25));

    t_0_x_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(26));

    t_0_x_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(27));

    t_0_x_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(28));

    t_0_x_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(29));

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

    // set up [SDSD]^(1) integral components

    t_0_zz_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(0));

    t_0_zz_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(1));

    t_0_zz_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(2));

    t_0_zz_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(3));

    t_0_zz_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(4));

    t_0_zz_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(5));

    t_0_yz_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(6));

    t_0_yy_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(7));

    t_0_yy_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(8));

    t_0_yy_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(9));

    t_0_yy_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(10));

    t_0_yy_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(11));

    t_0_yy_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(12));

    t_0_xx_0_zz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(13));

    t_0_xx_0_yz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(14));

    t_0_xx_0_yy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(15));

    t_0_xx_0_xz_1 = intsBufferSDSD1.data(intsIndexesSDSD1(16));

    t_0_xx_0_xy_1 = intsBufferSDSD1.data(intsIndexesSDSD1(17));

    t_0_xx_0_xx_1 = intsBufferSDSD1.data(intsIndexesSDSD1(18));

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

    t_0_yz_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(10));

    t_0_yz_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(11));

    t_0_yz_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(12));

    t_0_yz_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(13));

    t_0_yz_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(14));

    t_0_yy_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(15));

    t_0_yy_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(16));

    t_0_yy_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(17));

    t_0_yy_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(18));

    t_0_yy_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(19));

    t_0_yy_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(20));

    t_0_yy_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(21));

    t_0_yy_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(22));

    t_0_yy_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(23));

    t_0_yy_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(24));

    t_0_xz_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(25));

    t_0_xz_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(26));

    t_0_xz_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(27));

    t_0_xy_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(28));

    t_0_xy_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(29));

    t_0_xx_0_zzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(30));

    t_0_xx_0_yzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(31));

    t_0_xx_0_yyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(32));

    t_0_xx_0_yyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(33));

    t_0_xx_0_xzz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(34));

    t_0_xx_0_xyz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(35));

    t_0_xx_0_xyy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(36));

    t_0_xx_0_xxz_0 = intsBufferSDSF0.data(intsIndexesSDSF0(37));

    t_0_xx_0_xxy_0 = intsBufferSDSF0.data(intsIndexesSDSF0(38));

    t_0_xx_0_xxx_0 = intsBufferSDSF0.data(intsIndexesSDSF0(39));

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

    t_0_yz_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(10));

    t_0_yz_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(11));

    t_0_yz_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(12));

    t_0_yz_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(13));

    t_0_yz_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(14));

    t_0_yy_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(15));

    t_0_yy_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(16));

    t_0_yy_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(17));

    t_0_yy_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(18));

    t_0_yy_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(19));

    t_0_yy_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(20));

    t_0_yy_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(21));

    t_0_yy_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(22));

    t_0_yy_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(23));

    t_0_yy_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(24));

    t_0_xz_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(25));

    t_0_xz_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(26));

    t_0_xz_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(27));

    t_0_xy_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(28));

    t_0_xy_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(29));

    t_0_xx_0_zzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(30));

    t_0_xx_0_yzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(31));

    t_0_xx_0_yyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(32));

    t_0_xx_0_yyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(33));

    t_0_xx_0_xzz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(34));

    t_0_xx_0_xyz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(35));

    t_0_xx_0_xyy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(36));

    t_0_xx_0_xxz_1 = intsBufferSDSF1.data(intsIndexesSDSF1(37));

    t_0_xx_0_xxy_1 = intsBufferSDSF1.data(intsIndexesSDSF1(38));

    t_0_xx_0_xxx_1 = intsBufferSDSF1.data(intsIndexesSDSF1(39));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_y_0_xyz_0,\
                               t_0_y_0_xyz_1, t_0_y_0_xzz_0, t_0_y_0_xzz_1, t_0_y_0_yyy_0,\
                               t_0_y_0_yyy_1, t_0_y_0_yyz_0, t_0_y_0_yyz_1, t_0_y_0_yzz_0,\
                               t_0_y_0_yzz_1, t_0_y_0_zzz_0, t_0_y_0_zzz_1, t_0_yy_0_xx_1,\
                               t_0_yy_0_xxx_0, t_0_yy_0_xxx_1, t_0_yy_0_xxy_0, t_0_yy_0_xxy_1,\
                               t_0_yy_0_xxz_0, t_0_yy_0_xxz_1, t_0_yy_0_xy_1, t_0_yy_0_xyy_0,\
                               t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_xz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1,\
                               t_0_yy_0_yzz_0, t_0_yy_0_yzz_1, t_0_yy_0_zz_1, t_0_yy_0_zzz_0,\
                               t_0_yy_0_zzz_1, t_0_yyy_0_xyz_0, t_0_yyy_0_xzz_0, t_0_yyy_0_yyy_0,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yzz_0, t_0_yyy_0_zzz_0, t_0_yyz_0_xxx_0,\
                               t_0_yyz_0_xxy_0, t_0_yyz_0_xxz_0, t_0_yyz_0_xyy_0, t_0_yyz_0_xyz_0,\
                               t_0_yyz_0_xzz_0, t_0_yyz_0_yyy_0, t_0_yyz_0_yyz_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_zzz_0, t_0_yzz_0_xxx_0, t_0_yzz_0_xxy_0, t_0_yzz_0_xxz_0,\
                               t_0_yzz_0_xyy_0, t_0_yzz_0_xyz_0, t_0_yzz_0_xzz_0, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyz_0, t_0_yzz_0_yzz_0, t_0_yzz_0_zzz_0, t_0_z_0_xxx_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxy_0, t_0_z_0_xxy_1, t_0_z_0_xxz_0,\
                               t_0_z_0_xxz_1, t_0_z_0_xyy_0, t_0_z_0_xyy_1, t_0_z_0_xyz_0,\
                               t_0_z_0_xyz_1, t_0_z_0_xzz_0, t_0_z_0_xzz_1, t_0_z_0_yyy_0,\
                               t_0_z_0_yyy_1, t_0_z_0_yyz_0, t_0_z_0_yyz_1, t_0_z_0_yzz_0,\
                               t_0_z_0_yzz_1, t_0_z_0_zzz_0, t_0_z_0_zzz_1, t_0_zz_0_xx_1,\
                               t_0_zz_0_xxx_0, t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1,\
                               t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0,\
                               t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1,\
                               t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0,\
                               t_0_zz_0_zzz_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxy_0, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyz_0, t_0_zzz_0_xzz_0, t_0_zzz_0_yyy_0,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yzz_0, t_0_zzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_zzz_0[i] += rpb_z[i] * t_0_zz_0_zzz_0[i] + rwp_z[i] * t_0_zz_0_zzz_1[i] + fz_0[i] * t_0_z_0_zzz_0[i] - frz2_0[i] * t_0_z_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_zzz_0_yzz_0[i] += rpb_z[i] * t_0_zz_0_yzz_0[i] + rwp_z[i] * t_0_zz_0_yzz_1[i] + fz_0[i] * t_0_z_0_yzz_0[i] - frz2_0[i] * t_0_z_0_yzz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_zzz_0_yyz_0[i] += rpb_z[i] * t_0_zz_0_yyz_0[i] + rwp_z[i] * t_0_zz_0_yyz_1[i] + fz_0[i] * t_0_z_0_yyz_0[i] - frz2_0[i] * t_0_z_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzz_0_yyy_0[i] += rpb_z[i] * t_0_zz_0_yyy_0[i] + rwp_z[i] * t_0_zz_0_yyy_1[i] + fz_0[i] * t_0_z_0_yyy_0[i] - frz2_0[i] * t_0_z_0_yyy_1[i];

            t_0_zzz_0_xzz_0[i] += rpb_z[i] * t_0_zz_0_xzz_0[i] + rwp_z[i] * t_0_zz_0_xzz_1[i] + fz_0[i] * t_0_z_0_xzz_0[i] - frz2_0[i] * t_0_z_0_xzz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_zzz_0_xyz_0[i] += rpb_z[i] * t_0_zz_0_xyz_0[i] + rwp_z[i] * t_0_zz_0_xyz_1[i] + fz_0[i] * t_0_z_0_xyz_0[i] - frz2_0[i] * t_0_z_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzz_0_xyy_0[i] += rpb_z[i] * t_0_zz_0_xyy_0[i] + rwp_z[i] * t_0_zz_0_xyy_1[i] + fz_0[i] * t_0_z_0_xyy_0[i] - frz2_0[i] * t_0_z_0_xyy_1[i];

            t_0_zzz_0_xxz_0[i] += rpb_z[i] * t_0_zz_0_xxz_0[i] + rwp_z[i] * t_0_zz_0_xxz_1[i] + fz_0[i] * t_0_z_0_xxz_0[i] - frz2_0[i] * t_0_z_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_zzz_0_xxy_0[i] += rpb_z[i] * t_0_zz_0_xxy_0[i] + rwp_z[i] * t_0_zz_0_xxy_1[i] + fz_0[i] * t_0_z_0_xxy_0[i] - frz2_0[i] * t_0_z_0_xxy_1[i];

            t_0_zzz_0_xxx_0[i] += rpb_z[i] * t_0_zz_0_xxx_0[i] + rwp_z[i] * t_0_zz_0_xxx_1[i] + fz_0[i] * t_0_z_0_xxx_0[i] - frz2_0[i] * t_0_z_0_xxx_1[i];

            t_0_yzz_0_zzz_0[i] += rpb_y[i] * t_0_zz_0_zzz_0[i] + rwp_y[i] * t_0_zz_0_zzz_1[i];

            t_0_yzz_0_yzz_0[i] += rpb_y[i] * t_0_zz_0_yzz_0[i] + rwp_y[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_yzz_0_yyz_0[i] += rpb_y[i] * t_0_zz_0_yyz_0[i] + rwp_y[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_yzz_0_yyy_0[i] += rpb_y[i] * t_0_zz_0_yyy_0[i] + rwp_y[i] * t_0_zz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_yzz_0_xzz_0[i] += rpb_y[i] * t_0_zz_0_xzz_0[i] + rwp_y[i] * t_0_zz_0_xzz_1[i];

            t_0_yzz_0_xyz_0[i] += rpb_y[i] * t_0_zz_0_xyz_0[i] + rwp_y[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_yzz_0_xyy_0[i] += rpb_y[i] * t_0_zz_0_xyy_0[i] + rwp_y[i] * t_0_zz_0_xyy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_yzz_0_xxz_0[i] += rpb_y[i] * t_0_zz_0_xxz_0[i] + rwp_y[i] * t_0_zz_0_xxz_1[i];

            t_0_yzz_0_xxy_0[i] += rpb_y[i] * t_0_zz_0_xxy_0[i] + rwp_y[i] * t_0_zz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzz_0_xxx_0[i] += rpb_y[i] * t_0_zz_0_xxx_0[i] + rwp_y[i] * t_0_zz_0_xxx_1[i];

            t_0_yyz_0_zzz_0[i] += rpb_z[i] * t_0_yy_0_zzz_0[i] + rwp_z[i] * t_0_yy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyz_0_yzz_0[i] += rpb_z[i] * t_0_yy_0_yzz_0[i] + rwp_z[i] * t_0_yy_0_yzz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyz_0_yyz_0[i] += rpb_z[i] * t_0_yy_0_yyz_0[i] + rwp_z[i] * t_0_yy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyz_0_yyy_0[i] += rpb_z[i] * t_0_yy_0_yyy_0[i] + rwp_z[i] * t_0_yy_0_yyy_1[i];

            t_0_yyz_0_xzz_0[i] += rpb_z[i] * t_0_yy_0_xzz_0[i] + rwp_z[i] * t_0_yy_0_xzz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyz_0_xyz_0[i] += rpb_z[i] * t_0_yy_0_xyz_0[i] + rwp_z[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyz_0_xyy_0[i] += rpb_z[i] * t_0_yy_0_xyy_0[i] + rwp_z[i] * t_0_yy_0_xyy_1[i];

            t_0_yyz_0_xxz_0[i] += rpb_z[i] * t_0_yy_0_xxz_0[i] + rwp_z[i] * t_0_yy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyz_0_xxy_0[i] += rpb_z[i] * t_0_yy_0_xxy_0[i] + rwp_z[i] * t_0_yy_0_xxy_1[i];

            t_0_yyz_0_xxx_0[i] += rpb_z[i] * t_0_yy_0_xxx_0[i] + rwp_z[i] * t_0_yy_0_xxx_1[i];

            t_0_yyy_0_zzz_0[i] += rpb_y[i] * t_0_yy_0_zzz_0[i] + rwp_y[i] * t_0_yy_0_zzz_1[i] + fz_0[i] * t_0_y_0_zzz_0[i] - frz2_0[i] * t_0_y_0_zzz_1[i];

            t_0_yyy_0_yzz_0[i] += rpb_y[i] * t_0_yy_0_yzz_0[i] + rwp_y[i] * t_0_yy_0_yzz_1[i] + fz_0[i] * t_0_y_0_yzz_0[i] - frz2_0[i] * t_0_y_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyy_0_yyz_0[i] += rpb_y[i] * t_0_yy_0_yyz_0[i] + rwp_y[i] * t_0_yy_0_yyz_1[i] + fz_0[i] * t_0_y_0_yyz_0[i] - frz2_0[i] * t_0_y_0_yyz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyy_0_yyy_0[i] += rpb_y[i] * t_0_yy_0_yyy_0[i] + rwp_y[i] * t_0_yy_0_yyy_1[i] + fz_0[i] * t_0_y_0_yyy_0[i] - frz2_0[i] * t_0_y_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyy_0_xzz_0[i] += rpb_y[i] * t_0_yy_0_xzz_0[i] + rwp_y[i] * t_0_yy_0_xzz_1[i] + fz_0[i] * t_0_y_0_xzz_0[i] - frz2_0[i] * t_0_y_0_xzz_1[i];

            t_0_yyy_0_xyz_0[i] += rpb_y[i] * t_0_yy_0_xyz_0[i] + rwp_y[i] * t_0_yy_0_xyz_1[i] + fz_0[i] * t_0_y_0_xyz_0[i] - frz2_0[i] * t_0_y_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1,\
                               t_0_xx_0_zz_1, t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxz_0_yzz_0,\
                               t_0_xxz_0_zzz_0, t_0_xy_0_xxy_0, t_0_xy_0_xxy_1, t_0_xy_0_xyy_0,\
                               t_0_xy_0_xyy_1, t_0_xyy_0_xxx_0, t_0_xyy_0_xxy_0, t_0_xyy_0_xxz_0,\
                               t_0_xyy_0_xyy_0, t_0_xyy_0_xyz_0, t_0_xyy_0_xzz_0, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyz_0, t_0_xyy_0_yzz_0, t_0_xyy_0_zzz_0, t_0_xyz_0_xxx_0,\
                               t_0_xyz_0_xxy_0, t_0_xyz_0_xxz_0, t_0_xyz_0_xyy_0, t_0_xyz_0_xyz_0,\
                               t_0_xyz_0_xzz_0, t_0_xyz_0_yyy_0, t_0_xyz_0_yyz_0, t_0_xyz_0_yzz_0,\
                               t_0_xyz_0_zzz_0, t_0_xz_0_xxx_0, t_0_xz_0_xxx_1, t_0_xz_0_xxz_0,\
                               t_0_xz_0_xxz_1, t_0_xz_0_xzz_0, t_0_xz_0_xzz_1, t_0_xzz_0_xxx_0,\
                               t_0_xzz_0_xxy_0, t_0_xzz_0_xxz_0, t_0_xzz_0_xyy_0, t_0_xzz_0_xyz_0,\
                               t_0_xzz_0_xzz_0, t_0_xzz_0_yyy_0, t_0_xzz_0_yyz_0, t_0_xzz_0_yzz_0,\
                               t_0_xzz_0_zzz_0, t_0_y_0_xxx_0, t_0_y_0_xxx_1, t_0_y_0_xxy_0,\
                               t_0_y_0_xxy_1, t_0_y_0_xxz_0, t_0_y_0_xxz_1, t_0_y_0_xyy_0,\
                               t_0_y_0_xyy_1, t_0_yy_0_xx_1, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0,\
                               t_0_yy_0_xyz_1, t_0_yy_0_xz_1, t_0_yy_0_xzz_0, t_0_yy_0_xzz_1,\
                               t_0_yy_0_yy_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1, t_0_yy_0_yyz_0,\
                               t_0_yy_0_yyz_1, t_0_yy_0_yz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zz_1, t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxy_0, t_0_yyy_0_xxz_0, t_0_yyy_0_xyy_0, t_0_yz_0_xyz_0,\
                               t_0_yz_0_xyz_1, t_0_yz_0_yyy_0, t_0_yz_0_yyy_1, t_0_yz_0_yyz_0,\
                               t_0_yz_0_yyz_1, t_0_yz_0_yz_1, t_0_yz_0_yzz_0, t_0_yz_0_yzz_1,\
                               t_0_yz_0_zzz_0, t_0_yz_0_zzz_1, t_0_zz_0_xx_1, t_0_zz_0_xxx_0,\
                               t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1, t_0_zz_0_xxz_0,\
                               t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0, t_0_zz_0_xyy_1,\
                               t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1, t_0_zz_0_xzz_0,\
                               t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0, t_0_zz_0_yyy_1,\
                               t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1, t_0_zz_0_yzz_0,\
                               t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyy_0_xyy_0[i] += rpb_y[i] * t_0_yy_0_xyy_0[i] + rwp_y[i] * t_0_yy_0_xyy_1[i] + fz_0[i] * t_0_y_0_xyy_0[i] - frz2_0[i] * t_0_y_0_xyy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyy_0_xxz_0[i] += rpb_y[i] * t_0_yy_0_xxz_0[i] + rwp_y[i] * t_0_yy_0_xxz_1[i] + fz_0[i] * t_0_y_0_xxz_0[i] - frz2_0[i] * t_0_y_0_xxz_1[i];

            t_0_yyy_0_xxy_0[i] += rpb_y[i] * t_0_yy_0_xxy_0[i] + rwp_y[i] * t_0_yy_0_xxy_1[i] + fz_0[i] * t_0_y_0_xxy_0[i] - frz2_0[i] * t_0_y_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyy_0_xxx_0[i] += rpb_y[i] * t_0_yy_0_xxx_0[i] + rwp_y[i] * t_0_yy_0_xxx_1[i] + fz_0[i] * t_0_y_0_xxx_0[i] - frz2_0[i] * t_0_y_0_xxx_1[i];

            t_0_xzz_0_zzz_0[i] += rpb_x[i] * t_0_zz_0_zzz_0[i] + rwp_x[i] * t_0_zz_0_zzz_1[i];

            t_0_xzz_0_yzz_0[i] += rpb_x[i] * t_0_zz_0_yzz_0[i] + rwp_x[i] * t_0_zz_0_yzz_1[i];

            t_0_xzz_0_yyz_0[i] += rpb_x[i] * t_0_zz_0_yyz_0[i] + rwp_x[i] * t_0_zz_0_yyz_1[i];

            t_0_xzz_0_yyy_0[i] += rpb_x[i] * t_0_zz_0_yyy_0[i] + rwp_x[i] * t_0_zz_0_yyy_1[i];

            t_0_xzz_0_xzz_0[i] += rpb_x[i] * t_0_zz_0_xzz_0[i] + rwp_x[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_xzz_0_xyz_0[i] += rpb_x[i] * t_0_zz_0_xyz_0[i] + rwp_x[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_xzz_0_xyy_0[i] += rpb_x[i] * t_0_zz_0_xyy_0[i] + rwp_x[i] * t_0_zz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_xzz_0_xxz_0[i] += rpb_x[i] * t_0_zz_0_xxz_0[i] + rwp_x[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_xzz_0_xxy_0[i] += rpb_x[i] * t_0_zz_0_xxy_0[i] + rwp_x[i] * t_0_zz_0_xxy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_xzz_0_xxx_0[i] += rpb_x[i] * t_0_zz_0_xxx_0[i] + rwp_x[i] * t_0_zz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_xyz_0_zzz_0[i] += rpb_x[i] * t_0_yz_0_zzz_0[i] + rwp_x[i] * t_0_yz_0_zzz_1[i];

            t_0_xyz_0_yzz_0[i] += rpb_x[i] * t_0_yz_0_yzz_0[i] + rwp_x[i] * t_0_yz_0_yzz_1[i];

            t_0_xyz_0_yyz_0[i] += rpb_x[i] * t_0_yz_0_yyz_0[i] + rwp_x[i] * t_0_yz_0_yyz_1[i];

            t_0_xyz_0_yyy_0[i] += rpb_x[i] * t_0_yz_0_yyy_0[i] + rwp_x[i] * t_0_yz_0_yyy_1[i];

            t_0_xyz_0_xzz_0[i] += rpb_y[i] * t_0_xz_0_xzz_0[i] + rwp_y[i] * t_0_xz_0_xzz_1[i];

            t_0_xyz_0_xyz_0[i] += rpb_x[i] * t_0_yz_0_xyz_0[i] + rwp_x[i] * t_0_yz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yz_0_yz_1[i];

            t_0_xyz_0_xyy_0[i] += rpb_z[i] * t_0_xy_0_xyy_0[i] + rwp_z[i] * t_0_xy_0_xyy_1[i];

            t_0_xyz_0_xxz_0[i] += rpb_y[i] * t_0_xz_0_xxz_0[i] + rwp_y[i] * t_0_xz_0_xxz_1[i];

            t_0_xyz_0_xxy_0[i] += rpb_z[i] * t_0_xy_0_xxy_0[i] + rwp_z[i] * t_0_xy_0_xxy_1[i];

            t_0_xyz_0_xxx_0[i] += rpb_y[i] * t_0_xz_0_xxx_0[i] + rwp_y[i] * t_0_xz_0_xxx_1[i];

            t_0_xyy_0_zzz_0[i] += rpb_x[i] * t_0_yy_0_zzz_0[i] + rwp_x[i] * t_0_yy_0_zzz_1[i];

            t_0_xyy_0_yzz_0[i] += rpb_x[i] * t_0_yy_0_yzz_0[i] + rwp_x[i] * t_0_yy_0_yzz_1[i];

            t_0_xyy_0_yyz_0[i] += rpb_x[i] * t_0_yy_0_yyz_0[i] + rwp_x[i] * t_0_yy_0_yyz_1[i];

            t_0_xyy_0_yyy_0[i] += rpb_x[i] * t_0_yy_0_yyy_0[i] + rwp_x[i] * t_0_yy_0_yyy_1[i];

            t_0_xyy_0_xzz_0[i] += rpb_x[i] * t_0_yy_0_xzz_0[i] + rwp_x[i] * t_0_yy_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_xyy_0_xyz_0[i] += rpb_x[i] * t_0_yy_0_xyz_0[i] + rwp_x[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_xyy_0_xyy_0[i] += rpb_x[i] * t_0_yy_0_xyy_0[i] + rwp_x[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_xyy_0_xxz_0[i] += rpb_x[i] * t_0_yy_0_xxz_0[i] + rwp_x[i] * t_0_yy_0_xxz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_xyy_0_xxy_0[i] += rpb_x[i] * t_0_yy_0_xxy_0[i] + rwp_x[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_xyy_0_xxx_0[i] += rpb_x[i] * t_0_yy_0_xxx_0[i] + rwp_x[i] * t_0_yy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_xxz_0_zzz_0[i] += rpb_z[i] * t_0_xx_0_zzz_0[i] + rwp_z[i] * t_0_xx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxz_0_yzz_0[i] += rpb_z[i] * t_0_xx_0_yzz_0[i] + rwp_z[i] * t_0_xx_0_yzz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_xxx_0, t_0_x_0_xxx_1, t_0_x_0_xxy_0, t_0_x_0_xxy_1,\
                               t_0_x_0_xxz_0, t_0_x_0_xxz_1, t_0_x_0_xyy_0, t_0_x_0_xyy_1,\
                               t_0_x_0_xyz_0, t_0_x_0_xyz_1, t_0_x_0_xzz_0, t_0_x_0_xzz_1,\
                               t_0_x_0_yyy_0, t_0_x_0_yyy_1, t_0_x_0_yyz_0, t_0_x_0_yyz_1,\
                               t_0_x_0_yzz_0, t_0_x_0_yzz_1, t_0_x_0_zzz_0, t_0_x_0_zzz_1,\
                               t_0_xx_0_xx_1, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1, t_0_xx_0_xy_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yy_1,\
                               t_0_xx_0_yyy_0, t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1,\
                               t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1, t_0_xx_0_zz_1,\
                               t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_xxx_0, t_0_xxx_0_xxy_0,\
                               t_0_xxx_0_xxz_0, t_0_xxx_0_xyy_0, t_0_xxx_0_xyz_0, t_0_xxx_0_xzz_0,\
                               t_0_xxx_0_yyy_0, t_0_xxx_0_yyz_0, t_0_xxx_0_yzz_0, t_0_xxx_0_zzz_0,\
                               t_0_xxy_0_xxx_0, t_0_xxy_0_xxy_0, t_0_xxy_0_xxz_0, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xyz_0, t_0_xxy_0_xzz_0, t_0_xxy_0_yyy_0, t_0_xxy_0_yyz_0,\
                               t_0_xxy_0_yzz_0, t_0_xxy_0_zzz_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxz_0, t_0_xxz_0_xyy_0, t_0_xxz_0_xyz_0, t_0_xxz_0_xzz_0,\
                               t_0_xxz_0_yyy_0, t_0_xxz_0_yyz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxz_0_yyz_0[i] += rpb_z[i] * t_0_xx_0_yyz_0[i] + rwp_z[i] * t_0_xx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxz_0_yyy_0[i] += rpb_z[i] * t_0_xx_0_yyy_0[i] + rwp_z[i] * t_0_xx_0_yyy_1[i];

            t_0_xxz_0_xzz_0[i] += rpb_z[i] * t_0_xx_0_xzz_0[i] + rwp_z[i] * t_0_xx_0_xzz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxz_0_xyz_0[i] += rpb_z[i] * t_0_xx_0_xyz_0[i] + rwp_z[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxz_0_xyy_0[i] += rpb_z[i] * t_0_xx_0_xyy_0[i] + rwp_z[i] * t_0_xx_0_xyy_1[i];

            t_0_xxz_0_xxz_0[i] += rpb_z[i] * t_0_xx_0_xxz_0[i] + rwp_z[i] * t_0_xx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxz_0_xxy_0[i] += rpb_z[i] * t_0_xx_0_xxy_0[i] + rwp_z[i] * t_0_xx_0_xxy_1[i];

            t_0_xxz_0_xxx_0[i] += rpb_z[i] * t_0_xx_0_xxx_0[i] + rwp_z[i] * t_0_xx_0_xxx_1[i];

            t_0_xxy_0_zzz_0[i] += rpb_y[i] * t_0_xx_0_zzz_0[i] + rwp_y[i] * t_0_xx_0_zzz_1[i];

            t_0_xxy_0_yzz_0[i] += rpb_y[i] * t_0_xx_0_yzz_0[i] + rwp_y[i] * t_0_xx_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxy_0_yyz_0[i] += rpb_y[i] * t_0_xx_0_yyz_0[i] + rwp_y[i] * t_0_xx_0_yyz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxy_0_yyy_0[i] += rpb_y[i] * t_0_xx_0_yyy_0[i] + rwp_y[i] * t_0_xx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxy_0_xzz_0[i] += rpb_y[i] * t_0_xx_0_xzz_0[i] + rwp_y[i] * t_0_xx_0_xzz_1[i];

            t_0_xxy_0_xyz_0[i] += rpb_y[i] * t_0_xx_0_xyz_0[i] + rwp_y[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxy_0_xyy_0[i] += rpb_y[i] * t_0_xx_0_xyy_0[i] + rwp_y[i] * t_0_xx_0_xyy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxy_0_xxz_0[i] += rpb_y[i] * t_0_xx_0_xxz_0[i] + rwp_y[i] * t_0_xx_0_xxz_1[i];

            t_0_xxy_0_xxy_0[i] += rpb_y[i] * t_0_xx_0_xxy_0[i] + rwp_y[i] * t_0_xx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxy_0_xxx_0[i] += rpb_y[i] * t_0_xx_0_xxx_0[i] + rwp_y[i] * t_0_xx_0_xxx_1[i];

            t_0_xxx_0_zzz_0[i] += rpb_x[i] * t_0_xx_0_zzz_0[i] + rwp_x[i] * t_0_xx_0_zzz_1[i] + fz_0[i] * t_0_x_0_zzz_0[i] - frz2_0[i] * t_0_x_0_zzz_1[i];

            t_0_xxx_0_yzz_0[i] += rpb_x[i] * t_0_xx_0_yzz_0[i] + rwp_x[i] * t_0_xx_0_yzz_1[i] + fz_0[i] * t_0_x_0_yzz_0[i] - frz2_0[i] * t_0_x_0_yzz_1[i];

            t_0_xxx_0_yyz_0[i] += rpb_x[i] * t_0_xx_0_yyz_0[i] + rwp_x[i] * t_0_xx_0_yyz_1[i] + fz_0[i] * t_0_x_0_yyz_0[i] - frz2_0[i] * t_0_x_0_yyz_1[i];

            t_0_xxx_0_yyy_0[i] += rpb_x[i] * t_0_xx_0_yyy_0[i] + rwp_x[i] * t_0_xx_0_yyy_1[i] + fz_0[i] * t_0_x_0_yyy_0[i] - frz2_0[i] * t_0_x_0_yyy_1[i];

            t_0_xxx_0_xzz_0[i] += rpb_x[i] * t_0_xx_0_xzz_0[i] + rwp_x[i] * t_0_xx_0_xzz_1[i] + fz_0[i] * t_0_x_0_xzz_0[i] - frz2_0[i] * t_0_x_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxx_0_xyz_0[i] += rpb_x[i] * t_0_xx_0_xyz_0[i] + rwp_x[i] * t_0_xx_0_xyz_1[i] + fz_0[i] * t_0_x_0_xyz_0[i] - frz2_0[i] * t_0_x_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxx_0_xyy_0[i] += rpb_x[i] * t_0_xx_0_xyy_0[i] + rwp_x[i] * t_0_xx_0_xyy_1[i] + fz_0[i] * t_0_x_0_xyy_0[i] - frz2_0[i] * t_0_x_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxx_0_xxz_0[i] += rpb_x[i] * t_0_xx_0_xxz_0[i] + rwp_x[i] * t_0_xx_0_xxz_1[i] + fz_0[i] * t_0_x_0_xxz_0[i] - frz2_0[i] * t_0_x_0_xxz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxx_0_xxy_0[i] += rpb_x[i] * t_0_xx_0_xxy_0[i] + rwp_x[i] * t_0_xx_0_xxy_1[i] + fz_0[i] * t_0_x_0_xxy_0[i] - frz2_0[i] * t_0_x_0_xxy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxx_0_xxx_0[i] += rpb_x[i] * t_0_xx_0_xxx_0[i] + rwp_x[i] * t_0_xx_0_xxx_1[i] + fz_0[i] * t_0_x_0_xxx_0[i] - frz2_0[i] * t_0_x_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_xx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_y_0_xyz_0,\
                               t_0_y_0_xyz_1, t_0_y_0_xzz_0, t_0_y_0_xzz_1, t_0_y_0_yyy_0,\
                               t_0_y_0_yyy_1, t_0_y_0_yyz_0, t_0_y_0_yyz_1, t_0_y_0_yzz_0,\
                               t_0_y_0_yzz_1, t_0_y_0_zzz_0, t_0_y_0_zzz_1, t_0_yy_0_xx_1,\
                               t_0_yy_0_xxx_0, t_0_yy_0_xxx_1, t_0_yy_0_xxy_0, t_0_yy_0_xxy_1,\
                               t_0_yy_0_xxz_0, t_0_yy_0_xxz_1, t_0_yy_0_xy_1, t_0_yy_0_xyy_0,\
                               t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_xz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1,\
                               t_0_yy_0_yzz_0, t_0_yy_0_yzz_1, t_0_yy_0_zz_1, t_0_yy_0_zzz_0,\
                               t_0_yy_0_zzz_1, t_0_yyy_0_xyz_0, t_0_yyy_0_xzz_0, t_0_yyy_0_yyy_0,\
                               t_0_yyy_0_yyz_0, t_0_yyy_0_yzz_0, t_0_yyy_0_zzz_0, t_0_yyz_0_xxx_0,\
                               t_0_yyz_0_xxy_0, t_0_yyz_0_xxz_0, t_0_yyz_0_xyy_0, t_0_yyz_0_xyz_0,\
                               t_0_yyz_0_xzz_0, t_0_yyz_0_yyy_0, t_0_yyz_0_yyz_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_zzz_0, t_0_yzz_0_xxx_0, t_0_yzz_0_xxy_0, t_0_yzz_0_xxz_0,\
                               t_0_yzz_0_xyy_0, t_0_yzz_0_xyz_0, t_0_yzz_0_xzz_0, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyz_0, t_0_yzz_0_yzz_0, t_0_yzz_0_zzz_0, t_0_z_0_xxx_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxy_0, t_0_z_0_xxy_1, t_0_z_0_xxz_0,\
                               t_0_z_0_xxz_1, t_0_z_0_xyy_0, t_0_z_0_xyy_1, t_0_z_0_xyz_0,\
                               t_0_z_0_xyz_1, t_0_z_0_xzz_0, t_0_z_0_xzz_1, t_0_z_0_yyy_0,\
                               t_0_z_0_yyy_1, t_0_z_0_yyz_0, t_0_z_0_yyz_1, t_0_z_0_yzz_0,\
                               t_0_z_0_yzz_1, t_0_z_0_zzz_0, t_0_z_0_zzz_1, t_0_zz_0_xx_1,\
                               t_0_zz_0_xxx_0, t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1,\
                               t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0,\
                               t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1,\
                               t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0,\
                               t_0_zz_0_zzz_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxy_0, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyz_0, t_0_zzz_0_xzz_0, t_0_zzz_0_yyy_0,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yzz_0, t_0_zzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_zzz_0[i] = rpb_z[i] * t_0_zz_0_zzz_0[i] + rwp_z[i] * t_0_zz_0_zzz_1[i] + fz_0[i] * t_0_z_0_zzz_0[i] - frz2_0[i] * t_0_z_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_zzz_0_yzz_0[i] = rpb_z[i] * t_0_zz_0_yzz_0[i] + rwp_z[i] * t_0_zz_0_yzz_1[i] + fz_0[i] * t_0_z_0_yzz_0[i] - frz2_0[i] * t_0_z_0_yzz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_zzz_0_yyz_0[i] = rpb_z[i] * t_0_zz_0_yyz_0[i] + rwp_z[i] * t_0_zz_0_yyz_1[i] + fz_0[i] * t_0_z_0_yyz_0[i] - frz2_0[i] * t_0_z_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzz_0_yyy_0[i] = rpb_z[i] * t_0_zz_0_yyy_0[i] + rwp_z[i] * t_0_zz_0_yyy_1[i] + fz_0[i] * t_0_z_0_yyy_0[i] - frz2_0[i] * t_0_z_0_yyy_1[i];

            t_0_zzz_0_xzz_0[i] = rpb_z[i] * t_0_zz_0_xzz_0[i] + rwp_z[i] * t_0_zz_0_xzz_1[i] + fz_0[i] * t_0_z_0_xzz_0[i] - frz2_0[i] * t_0_z_0_xzz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_zzz_0_xyz_0[i] = rpb_z[i] * t_0_zz_0_xyz_0[i] + rwp_z[i] * t_0_zz_0_xyz_1[i] + fz_0[i] * t_0_z_0_xyz_0[i] - frz2_0[i] * t_0_z_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzz_0_xyy_0[i] = rpb_z[i] * t_0_zz_0_xyy_0[i] + rwp_z[i] * t_0_zz_0_xyy_1[i] + fz_0[i] * t_0_z_0_xyy_0[i] - frz2_0[i] * t_0_z_0_xyy_1[i];

            t_0_zzz_0_xxz_0[i] = rpb_z[i] * t_0_zz_0_xxz_0[i] + rwp_z[i] * t_0_zz_0_xxz_1[i] + fz_0[i] * t_0_z_0_xxz_0[i] - frz2_0[i] * t_0_z_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_zzz_0_xxy_0[i] = rpb_z[i] * t_0_zz_0_xxy_0[i] + rwp_z[i] * t_0_zz_0_xxy_1[i] + fz_0[i] * t_0_z_0_xxy_0[i] - frz2_0[i] * t_0_z_0_xxy_1[i];

            t_0_zzz_0_xxx_0[i] = rpb_z[i] * t_0_zz_0_xxx_0[i] + rwp_z[i] * t_0_zz_0_xxx_1[i] + fz_0[i] * t_0_z_0_xxx_0[i] - frz2_0[i] * t_0_z_0_xxx_1[i];

            t_0_yzz_0_zzz_0[i] = rpb_y[i] * t_0_zz_0_zzz_0[i] + rwp_y[i] * t_0_zz_0_zzz_1[i];

            t_0_yzz_0_yzz_0[i] = rpb_y[i] * t_0_zz_0_yzz_0[i] + rwp_y[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_yzz_0_yyz_0[i] = rpb_y[i] * t_0_zz_0_yyz_0[i] + rwp_y[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_yzz_0_yyy_0[i] = rpb_y[i] * t_0_zz_0_yyy_0[i] + rwp_y[i] * t_0_zz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_yzz_0_xzz_0[i] = rpb_y[i] * t_0_zz_0_xzz_0[i] + rwp_y[i] * t_0_zz_0_xzz_1[i];

            t_0_yzz_0_xyz_0[i] = rpb_y[i] * t_0_zz_0_xyz_0[i] + rwp_y[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_yzz_0_xyy_0[i] = rpb_y[i] * t_0_zz_0_xyy_0[i] + rwp_y[i] * t_0_zz_0_xyy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_yzz_0_xxz_0[i] = rpb_y[i] * t_0_zz_0_xxz_0[i] + rwp_y[i] * t_0_zz_0_xxz_1[i];

            t_0_yzz_0_xxy_0[i] = rpb_y[i] * t_0_zz_0_xxy_0[i] + rwp_y[i] * t_0_zz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzz_0_xxx_0[i] = rpb_y[i] * t_0_zz_0_xxx_0[i] + rwp_y[i] * t_0_zz_0_xxx_1[i];

            t_0_yyz_0_zzz_0[i] = rpb_z[i] * t_0_yy_0_zzz_0[i] + rwp_z[i] * t_0_yy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyz_0_yzz_0[i] = rpb_z[i] * t_0_yy_0_yzz_0[i] + rwp_z[i] * t_0_yy_0_yzz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyz_0_yyz_0[i] = rpb_z[i] * t_0_yy_0_yyz_0[i] + rwp_z[i] * t_0_yy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyz_0_yyy_0[i] = rpb_z[i] * t_0_yy_0_yyy_0[i] + rwp_z[i] * t_0_yy_0_yyy_1[i];

            t_0_yyz_0_xzz_0[i] = rpb_z[i] * t_0_yy_0_xzz_0[i] + rwp_z[i] * t_0_yy_0_xzz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyz_0_xyz_0[i] = rpb_z[i] * t_0_yy_0_xyz_0[i] + rwp_z[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyz_0_xyy_0[i] = rpb_z[i] * t_0_yy_0_xyy_0[i] + rwp_z[i] * t_0_yy_0_xyy_1[i];

            t_0_yyz_0_xxz_0[i] = rpb_z[i] * t_0_yy_0_xxz_0[i] + rwp_z[i] * t_0_yy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyz_0_xxy_0[i] = rpb_z[i] * t_0_yy_0_xxy_0[i] + rwp_z[i] * t_0_yy_0_xxy_1[i];

            t_0_yyz_0_xxx_0[i] = rpb_z[i] * t_0_yy_0_xxx_0[i] + rwp_z[i] * t_0_yy_0_xxx_1[i];

            t_0_yyy_0_zzz_0[i] = rpb_y[i] * t_0_yy_0_zzz_0[i] + rwp_y[i] * t_0_yy_0_zzz_1[i] + fz_0[i] * t_0_y_0_zzz_0[i] - frz2_0[i] * t_0_y_0_zzz_1[i];

            t_0_yyy_0_yzz_0[i] = rpb_y[i] * t_0_yy_0_yzz_0[i] + rwp_y[i] * t_0_yy_0_yzz_1[i] + fz_0[i] * t_0_y_0_yzz_0[i] - frz2_0[i] * t_0_y_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyy_0_yyz_0[i] = rpb_y[i] * t_0_yy_0_yyz_0[i] + rwp_y[i] * t_0_yy_0_yyz_1[i] + fz_0[i] * t_0_y_0_yyz_0[i] - frz2_0[i] * t_0_y_0_yyz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyy_0_yyy_0[i] = rpb_y[i] * t_0_yy_0_yyy_0[i] + rwp_y[i] * t_0_yy_0_yyy_1[i] + fz_0[i] * t_0_y_0_yyy_0[i] - frz2_0[i] * t_0_y_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyy_0_xzz_0[i] = rpb_y[i] * t_0_yy_0_xzz_0[i] + rwp_y[i] * t_0_yy_0_xzz_1[i] + fz_0[i] * t_0_y_0_xzz_0[i] - frz2_0[i] * t_0_y_0_xzz_1[i];

            t_0_yyy_0_xyz_0[i] = rpb_y[i] * t_0_yy_0_xyz_0[i] + rwp_y[i] * t_0_yy_0_xyz_1[i] + fz_0[i] * t_0_y_0_xyz_0[i] - frz2_0[i] * t_0_y_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1,\
                               t_0_xx_0_zz_1, t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxz_0_yzz_0,\
                               t_0_xxz_0_zzz_0, t_0_xy_0_xxy_0, t_0_xy_0_xxy_1, t_0_xy_0_xyy_0,\
                               t_0_xy_0_xyy_1, t_0_xyy_0_xxx_0, t_0_xyy_0_xxy_0, t_0_xyy_0_xxz_0,\
                               t_0_xyy_0_xyy_0, t_0_xyy_0_xyz_0, t_0_xyy_0_xzz_0, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyz_0, t_0_xyy_0_yzz_0, t_0_xyy_0_zzz_0, t_0_xyz_0_xxx_0,\
                               t_0_xyz_0_xxy_0, t_0_xyz_0_xxz_0, t_0_xyz_0_xyy_0, t_0_xyz_0_xyz_0,\
                               t_0_xyz_0_xzz_0, t_0_xyz_0_yyy_0, t_0_xyz_0_yyz_0, t_0_xyz_0_yzz_0,\
                               t_0_xyz_0_zzz_0, t_0_xz_0_xxx_0, t_0_xz_0_xxx_1, t_0_xz_0_xxz_0,\
                               t_0_xz_0_xxz_1, t_0_xz_0_xzz_0, t_0_xz_0_xzz_1, t_0_xzz_0_xxx_0,\
                               t_0_xzz_0_xxy_0, t_0_xzz_0_xxz_0, t_0_xzz_0_xyy_0, t_0_xzz_0_xyz_0,\
                               t_0_xzz_0_xzz_0, t_0_xzz_0_yyy_0, t_0_xzz_0_yyz_0, t_0_xzz_0_yzz_0,\
                               t_0_xzz_0_zzz_0, t_0_y_0_xxx_0, t_0_y_0_xxx_1, t_0_y_0_xxy_0,\
                               t_0_y_0_xxy_1, t_0_y_0_xxz_0, t_0_y_0_xxz_1, t_0_y_0_xyy_0,\
                               t_0_y_0_xyy_1, t_0_yy_0_xx_1, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0,\
                               t_0_yy_0_xyz_1, t_0_yy_0_xz_1, t_0_yy_0_xzz_0, t_0_yy_0_xzz_1,\
                               t_0_yy_0_yy_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1, t_0_yy_0_yyz_0,\
                               t_0_yy_0_yyz_1, t_0_yy_0_yz_1, t_0_yy_0_yzz_0, t_0_yy_0_yzz_1,\
                               t_0_yy_0_zz_1, t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxy_0, t_0_yyy_0_xxz_0, t_0_yyy_0_xyy_0, t_0_yz_0_xyz_0,\
                               t_0_yz_0_xyz_1, t_0_yz_0_yyy_0, t_0_yz_0_yyy_1, t_0_yz_0_yyz_0,\
                               t_0_yz_0_yyz_1, t_0_yz_0_yz_1, t_0_yz_0_yzz_0, t_0_yz_0_yzz_1,\
                               t_0_yz_0_zzz_0, t_0_yz_0_zzz_1, t_0_zz_0_xx_1, t_0_zz_0_xxx_0,\
                               t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1, t_0_zz_0_xxz_0,\
                               t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0, t_0_zz_0_xyy_1,\
                               t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1, t_0_zz_0_xzz_0,\
                               t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0, t_0_zz_0_yyy_1,\
                               t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1, t_0_zz_0_yzz_0,\
                               t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyy_0_xyy_0[i] = rpb_y[i] * t_0_yy_0_xyy_0[i] + rwp_y[i] * t_0_yy_0_xyy_1[i] + fz_0[i] * t_0_y_0_xyy_0[i] - frz2_0[i] * t_0_y_0_xyy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyy_0_xxz_0[i] = rpb_y[i] * t_0_yy_0_xxz_0[i] + rwp_y[i] * t_0_yy_0_xxz_1[i] + fz_0[i] * t_0_y_0_xxz_0[i] - frz2_0[i] * t_0_y_0_xxz_1[i];

            t_0_yyy_0_xxy_0[i] = rpb_y[i] * t_0_yy_0_xxy_0[i] + rwp_y[i] * t_0_yy_0_xxy_1[i] + fz_0[i] * t_0_y_0_xxy_0[i] - frz2_0[i] * t_0_y_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyy_0_xxx_0[i] = rpb_y[i] * t_0_yy_0_xxx_0[i] + rwp_y[i] * t_0_yy_0_xxx_1[i] + fz_0[i] * t_0_y_0_xxx_0[i] - frz2_0[i] * t_0_y_0_xxx_1[i];

            t_0_xzz_0_zzz_0[i] = rpb_x[i] * t_0_zz_0_zzz_0[i] + rwp_x[i] * t_0_zz_0_zzz_1[i];

            t_0_xzz_0_yzz_0[i] = rpb_x[i] * t_0_zz_0_yzz_0[i] + rwp_x[i] * t_0_zz_0_yzz_1[i];

            t_0_xzz_0_yyz_0[i] = rpb_x[i] * t_0_zz_0_yyz_0[i] + rwp_x[i] * t_0_zz_0_yyz_1[i];

            t_0_xzz_0_yyy_0[i] = rpb_x[i] * t_0_zz_0_yyy_0[i] + rwp_x[i] * t_0_zz_0_yyy_1[i];

            t_0_xzz_0_xzz_0[i] = rpb_x[i] * t_0_zz_0_xzz_0[i] + rwp_x[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_xzz_0_xyz_0[i] = rpb_x[i] * t_0_zz_0_xyz_0[i] + rwp_x[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_xzz_0_xyy_0[i] = rpb_x[i] * t_0_zz_0_xyy_0[i] + rwp_x[i] * t_0_zz_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_xzz_0_xxz_0[i] = rpb_x[i] * t_0_zz_0_xxz_0[i] + rwp_x[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_xzz_0_xxy_0[i] = rpb_x[i] * t_0_zz_0_xxy_0[i] + rwp_x[i] * t_0_zz_0_xxy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_xzz_0_xxx_0[i] = rpb_x[i] * t_0_zz_0_xxx_0[i] + rwp_x[i] * t_0_zz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_xyz_0_zzz_0[i] = rpb_x[i] * t_0_yz_0_zzz_0[i] + rwp_x[i] * t_0_yz_0_zzz_1[i];

            t_0_xyz_0_yzz_0[i] = rpb_x[i] * t_0_yz_0_yzz_0[i] + rwp_x[i] * t_0_yz_0_yzz_1[i];

            t_0_xyz_0_yyz_0[i] = rpb_x[i] * t_0_yz_0_yyz_0[i] + rwp_x[i] * t_0_yz_0_yyz_1[i];

            t_0_xyz_0_yyy_0[i] = rpb_x[i] * t_0_yz_0_yyy_0[i] + rwp_x[i] * t_0_yz_0_yyy_1[i];

            t_0_xyz_0_xzz_0[i] = rpb_y[i] * t_0_xz_0_xzz_0[i] + rwp_y[i] * t_0_xz_0_xzz_1[i];

            t_0_xyz_0_xyz_0[i] = rpb_x[i] * t_0_yz_0_xyz_0[i] + rwp_x[i] * t_0_yz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yz_0_yz_1[i];

            t_0_xyz_0_xyy_0[i] = rpb_z[i] * t_0_xy_0_xyy_0[i] + rwp_z[i] * t_0_xy_0_xyy_1[i];

            t_0_xyz_0_xxz_0[i] = rpb_y[i] * t_0_xz_0_xxz_0[i] + rwp_y[i] * t_0_xz_0_xxz_1[i];

            t_0_xyz_0_xxy_0[i] = rpb_z[i] * t_0_xy_0_xxy_0[i] + rwp_z[i] * t_0_xy_0_xxy_1[i];

            t_0_xyz_0_xxx_0[i] = rpb_y[i] * t_0_xz_0_xxx_0[i] + rwp_y[i] * t_0_xz_0_xxx_1[i];

            t_0_xyy_0_zzz_0[i] = rpb_x[i] * t_0_yy_0_zzz_0[i] + rwp_x[i] * t_0_yy_0_zzz_1[i];

            t_0_xyy_0_yzz_0[i] = rpb_x[i] * t_0_yy_0_yzz_0[i] + rwp_x[i] * t_0_yy_0_yzz_1[i];

            t_0_xyy_0_yyz_0[i] = rpb_x[i] * t_0_yy_0_yyz_0[i] + rwp_x[i] * t_0_yy_0_yyz_1[i];

            t_0_xyy_0_yyy_0[i] = rpb_x[i] * t_0_yy_0_yyy_0[i] + rwp_x[i] * t_0_yy_0_yyy_1[i];

            t_0_xyy_0_xzz_0[i] = rpb_x[i] * t_0_yy_0_xzz_0[i] + rwp_x[i] * t_0_yy_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_xyy_0_xyz_0[i] = rpb_x[i] * t_0_yy_0_xyz_0[i] + rwp_x[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_xyy_0_xyy_0[i] = rpb_x[i] * t_0_yy_0_xyy_0[i] + rwp_x[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_xyy_0_xxz_0[i] = rpb_x[i] * t_0_yy_0_xxz_0[i] + rwp_x[i] * t_0_yy_0_xxz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_xyy_0_xxy_0[i] = rpb_x[i] * t_0_yy_0_xxy_0[i] + rwp_x[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_xyy_0_xxx_0[i] = rpb_x[i] * t_0_yy_0_xxx_0[i] + rwp_x[i] * t_0_yy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_xxz_0_zzz_0[i] = rpb_z[i] * t_0_xx_0_zzz_0[i] + rwp_z[i] * t_0_xx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxz_0_yzz_0[i] = rpb_z[i] * t_0_xx_0_yzz_0[i] + rwp_z[i] * t_0_xx_0_yzz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_xxx_0, t_0_x_0_xxx_1, t_0_x_0_xxy_0, t_0_x_0_xxy_1,\
                               t_0_x_0_xxz_0, t_0_x_0_xxz_1, t_0_x_0_xyy_0, t_0_x_0_xyy_1,\
                               t_0_x_0_xyz_0, t_0_x_0_xyz_1, t_0_x_0_xzz_0, t_0_x_0_xzz_1,\
                               t_0_x_0_yyy_0, t_0_x_0_yyy_1, t_0_x_0_yyz_0, t_0_x_0_yyz_1,\
                               t_0_x_0_yzz_0, t_0_x_0_yzz_1, t_0_x_0_zzz_0, t_0_x_0_zzz_1,\
                               t_0_xx_0_xx_1, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1, t_0_xx_0_xy_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yy_1,\
                               t_0_xx_0_yyy_0, t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1,\
                               t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1, t_0_xx_0_zz_1,\
                               t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_xxx_0, t_0_xxx_0_xxy_0,\
                               t_0_xxx_0_xxz_0, t_0_xxx_0_xyy_0, t_0_xxx_0_xyz_0, t_0_xxx_0_xzz_0,\
                               t_0_xxx_0_yyy_0, t_0_xxx_0_yyz_0, t_0_xxx_0_yzz_0, t_0_xxx_0_zzz_0,\
                               t_0_xxy_0_xxx_0, t_0_xxy_0_xxy_0, t_0_xxy_0_xxz_0, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xyz_0, t_0_xxy_0_xzz_0, t_0_xxy_0_yyy_0, t_0_xxy_0_yyz_0,\
                               t_0_xxy_0_yzz_0, t_0_xxy_0_zzz_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxz_0, t_0_xxz_0_xyy_0, t_0_xxz_0_xyz_0, t_0_xxz_0_xzz_0,\
                               t_0_xxz_0_yyy_0, t_0_xxz_0_yyz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxz_0_yyz_0[i] = rpb_z[i] * t_0_xx_0_yyz_0[i] + rwp_z[i] * t_0_xx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxz_0_yyy_0[i] = rpb_z[i] * t_0_xx_0_yyy_0[i] + rwp_z[i] * t_0_xx_0_yyy_1[i];

            t_0_xxz_0_xzz_0[i] = rpb_z[i] * t_0_xx_0_xzz_0[i] + rwp_z[i] * t_0_xx_0_xzz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxz_0_xyz_0[i] = rpb_z[i] * t_0_xx_0_xyz_0[i] + rwp_z[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxz_0_xyy_0[i] = rpb_z[i] * t_0_xx_0_xyy_0[i] + rwp_z[i] * t_0_xx_0_xyy_1[i];

            t_0_xxz_0_xxz_0[i] = rpb_z[i] * t_0_xx_0_xxz_0[i] + rwp_z[i] * t_0_xx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxz_0_xxy_0[i] = rpb_z[i] * t_0_xx_0_xxy_0[i] + rwp_z[i] * t_0_xx_0_xxy_1[i];

            t_0_xxz_0_xxx_0[i] = rpb_z[i] * t_0_xx_0_xxx_0[i] + rwp_z[i] * t_0_xx_0_xxx_1[i];

            t_0_xxy_0_zzz_0[i] = rpb_y[i] * t_0_xx_0_zzz_0[i] + rwp_y[i] * t_0_xx_0_zzz_1[i];

            t_0_xxy_0_yzz_0[i] = rpb_y[i] * t_0_xx_0_yzz_0[i] + rwp_y[i] * t_0_xx_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxy_0_yyz_0[i] = rpb_y[i] * t_0_xx_0_yyz_0[i] + rwp_y[i] * t_0_xx_0_yyz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxy_0_yyy_0[i] = rpb_y[i] * t_0_xx_0_yyy_0[i] + rwp_y[i] * t_0_xx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxy_0_xzz_0[i] = rpb_y[i] * t_0_xx_0_xzz_0[i] + rwp_y[i] * t_0_xx_0_xzz_1[i];

            t_0_xxy_0_xyz_0[i] = rpb_y[i] * t_0_xx_0_xyz_0[i] + rwp_y[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxy_0_xyy_0[i] = rpb_y[i] * t_0_xx_0_xyy_0[i] + rwp_y[i] * t_0_xx_0_xyy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxy_0_xxz_0[i] = rpb_y[i] * t_0_xx_0_xxz_0[i] + rwp_y[i] * t_0_xx_0_xxz_1[i];

            t_0_xxy_0_xxy_0[i] = rpb_y[i] * t_0_xx_0_xxy_0[i] + rwp_y[i] * t_0_xx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxy_0_xxx_0[i] = rpb_y[i] * t_0_xx_0_xxx_0[i] + rwp_y[i] * t_0_xx_0_xxx_1[i];

            t_0_xxx_0_zzz_0[i] = rpb_x[i] * t_0_xx_0_zzz_0[i] + rwp_x[i] * t_0_xx_0_zzz_1[i] + fz_0[i] * t_0_x_0_zzz_0[i] - frz2_0[i] * t_0_x_0_zzz_1[i];

            t_0_xxx_0_yzz_0[i] = rpb_x[i] * t_0_xx_0_yzz_0[i] + rwp_x[i] * t_0_xx_0_yzz_1[i] + fz_0[i] * t_0_x_0_yzz_0[i] - frz2_0[i] * t_0_x_0_yzz_1[i];

            t_0_xxx_0_yyz_0[i] = rpb_x[i] * t_0_xx_0_yyz_0[i] + rwp_x[i] * t_0_xx_0_yyz_1[i] + fz_0[i] * t_0_x_0_yyz_0[i] - frz2_0[i] * t_0_x_0_yyz_1[i];

            t_0_xxx_0_yyy_0[i] = rpb_x[i] * t_0_xx_0_yyy_0[i] + rwp_x[i] * t_0_xx_0_yyy_1[i] + fz_0[i] * t_0_x_0_yyy_0[i] - frz2_0[i] * t_0_x_0_yyy_1[i];

            t_0_xxx_0_xzz_0[i] = rpb_x[i] * t_0_xx_0_xzz_0[i] + rwp_x[i] * t_0_xx_0_xzz_1[i] + fz_0[i] * t_0_x_0_xzz_0[i] - frz2_0[i] * t_0_x_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxx_0_xyz_0[i] = rpb_x[i] * t_0_xx_0_xyz_0[i] + rwp_x[i] * t_0_xx_0_xyz_1[i] + fz_0[i] * t_0_x_0_xyz_0[i] - frz2_0[i] * t_0_x_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxx_0_xyy_0[i] = rpb_x[i] * t_0_xx_0_xyy_0[i] + rwp_x[i] * t_0_xx_0_xyy_1[i] + fz_0[i] * t_0_x_0_xyy_0[i] - frz2_0[i] * t_0_x_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxx_0_xxz_0[i] = rpb_x[i] * t_0_xx_0_xxz_0[i] + rwp_x[i] * t_0_xx_0_xxz_1[i] + fz_0[i] * t_0_x_0_xxz_0[i] - frz2_0[i] * t_0_x_0_xxz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxx_0_xxy_0[i] = rpb_x[i] * t_0_xx_0_xxy_0[i] + rwp_x[i] * t_0_xx_0_xxy_1[i] + fz_0[i] * t_0_x_0_xxy_0[i] - frz2_0[i] * t_0_x_0_xxy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxx_0_xxx_0[i] = rpb_x[i] * t_0_xx_0_xxx_0[i] + rwp_x[i] * t_0_xx_0_xxx_1[i] + fz_0[i] * t_0_x_0_xxx_0[i] - frz2_0[i] * t_0_x_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_xx_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSFSF_V1(      BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF0,
                      const BufferHostXY<T>&      intsBufferSPSF0,
                      const BufferHostX<int32_t>& intsIndexesSPSF0,
                      const BufferHostXY<T>&      intsBufferSPSF1,
                      const BufferHostX<int32_t>& intsIndexesSPSF1,
                      const BufferHostXY<T>&      intsBufferSDSD1,
                      const BufferHostX<int32_t>& intsIndexesSDSD1,
                      const BufferHostXY<T>&      intsBufferSDSF0,
                      const BufferHostX<int32_t>& intsIndexesSDSF0,
                      const BufferHostXY<T>&      intsBufferSDSF1,
                      const BufferHostX<int32_t>& intsIndexesSDSF1,
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

    // set up [SPSF]^(0) integral components

    t_0_z_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(0));

    t_0_z_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(1));

    t_0_z_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(2));

    t_0_z_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(3));

    t_0_z_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(4));

    t_0_z_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(5));

    t_0_z_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(6));

    t_0_z_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(7));

    t_0_z_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(8));

    t_0_z_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(9));

    t_0_y_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(10));

    t_0_y_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(11));

    t_0_y_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(12));

    t_0_y_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(13));

    t_0_y_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(14));

    t_0_y_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(15));

    t_0_y_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(16));

    t_0_y_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(17));

    t_0_y_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(18));

    t_0_y_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(19));

    t_0_x_0_zzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(20));

    t_0_x_0_yzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(21));

    t_0_x_0_yyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(22));

    t_0_x_0_yyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(23));

    t_0_x_0_xzz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(24));

    t_0_x_0_xyz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(25));

    t_0_x_0_xyy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(26));

    t_0_x_0_xxz_0 = intsBufferSPSF0.data(intsIndexesSPSF0(27));

    t_0_x_0_xxy_0 = intsBufferSPSF0.data(intsIndexesSPSF0(28));

    t_0_x_0_xxx_0 = intsBufferSPSF0.data(intsIndexesSPSF0(29));

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

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_y_0_xyy_0,\
                               t_0_y_0_xyy_1, t_0_y_0_xyz_0, t_0_y_0_xyz_1, t_0_y_0_xzz_0,\
                               t_0_y_0_xzz_1, t_0_y_0_yyy_0, t_0_y_0_yyy_1, t_0_y_0_yyz_0,\
                               t_0_y_0_yyz_1, t_0_y_0_yzz_0, t_0_y_0_yzz_1, t_0_y_0_zzz_0,\
                               t_0_y_0_zzz_1, t_0_yy_0_xx_1, t_0_yy_0_xxy_0, t_0_yy_0_xxy_1,\
                               t_0_yy_0_xxz_0, t_0_yy_0_xxz_1, t_0_yy_0_xy_1, t_0_yy_0_xyy_0,\
                               t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_xz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1,\
                               t_0_yy_0_yzz_0, t_0_yy_0_yzz_1, t_0_yy_0_zz_1, t_0_yy_0_zzz_0,\
                               t_0_yy_0_zzz_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyz_0, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_yyy_0, t_0_yyy_0_yyz_0, t_0_yyy_0_yzz_0, t_0_yyy_0_zzz_0,\
                               t_0_yyz_0_xxy_0, t_0_yyz_0_xxz_0, t_0_yyz_0_xyy_0, t_0_yyz_0_xyz_0,\
                               t_0_yyz_0_xzz_0, t_0_yyz_0_yyy_0, t_0_yyz_0_yyz_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_zzz_0, t_0_yzz_0_xxx_0, t_0_yzz_0_xxy_0, t_0_yzz_0_xxz_0,\
                               t_0_yzz_0_xyy_0, t_0_yzz_0_xyz_0, t_0_yzz_0_xzz_0, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyz_0, t_0_yzz_0_yzz_0, t_0_yzz_0_zzz_0, t_0_z_0_xxx_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxy_0, t_0_z_0_xxy_1, t_0_z_0_xxz_0,\
                               t_0_z_0_xxz_1, t_0_z_0_xyy_0, t_0_z_0_xyy_1, t_0_z_0_xyz_0,\
                               t_0_z_0_xyz_1, t_0_z_0_xzz_0, t_0_z_0_xzz_1, t_0_z_0_yyy_0,\
                               t_0_z_0_yyy_1, t_0_z_0_yyz_0, t_0_z_0_yyz_1, t_0_z_0_yzz_0,\
                               t_0_z_0_yzz_1, t_0_z_0_zzz_0, t_0_z_0_zzz_1, t_0_zz_0_xx_1,\
                               t_0_zz_0_xxx_0, t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1,\
                               t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0,\
                               t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1,\
                               t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0,\
                               t_0_zz_0_zzz_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxy_0, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyz_0, t_0_zzz_0_xzz_0, t_0_zzz_0_yyy_0,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yzz_0, t_0_zzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_zzz_0[i] += rpb_z[i] * t_0_zz_0_zzz_0[i] + rwp_z[i] * t_0_zz_0_zzz_1[i] + fz_0[i] * t_0_z_0_zzz_0[i] - frz2_0[i] * t_0_z_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_zzz_0_yzz_0[i] += rpb_z[i] * t_0_zz_0_yzz_0[i] + rwp_z[i] * t_0_zz_0_yzz_1[i] + fz_0[i] * t_0_z_0_yzz_0[i] - frz2_0[i] * t_0_z_0_yzz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_zzz_0_yyz_0[i] += rpb_z[i] * t_0_zz_0_yyz_0[i] + rwp_z[i] * t_0_zz_0_yyz_1[i] + fz_0[i] * t_0_z_0_yyz_0[i] - frz2_0[i] * t_0_z_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzz_0_yyy_0[i] += rpb_z[i] * t_0_zz_0_yyy_0[i] + rwp_z[i] * t_0_zz_0_yyy_1[i] + fz_0[i] * t_0_z_0_yyy_0[i] - frz2_0[i] * t_0_z_0_yyy_1[i];

            t_0_zzz_0_xzz_0[i] += rpb_z[i] * t_0_zz_0_xzz_0[i] + rwp_z[i] * t_0_zz_0_xzz_1[i] + fz_0[i] * t_0_z_0_xzz_0[i] - frz2_0[i] * t_0_z_0_xzz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_zzz_0_xyz_0[i] += rpb_z[i] * t_0_zz_0_xyz_0[i] + rwp_z[i] * t_0_zz_0_xyz_1[i] + fz_0[i] * t_0_z_0_xyz_0[i] - frz2_0[i] * t_0_z_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzz_0_xyy_0[i] += rpb_z[i] * t_0_zz_0_xyy_0[i] + rwp_z[i] * t_0_zz_0_xyy_1[i] + fz_0[i] * t_0_z_0_xyy_0[i] - frz2_0[i] * t_0_z_0_xyy_1[i];

            t_0_zzz_0_xxz_0[i] += rpb_z[i] * t_0_zz_0_xxz_0[i] + rwp_z[i] * t_0_zz_0_xxz_1[i] + fz_0[i] * t_0_z_0_xxz_0[i] - frz2_0[i] * t_0_z_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_zzz_0_xxy_0[i] += rpb_z[i] * t_0_zz_0_xxy_0[i] + rwp_z[i] * t_0_zz_0_xxy_1[i] + fz_0[i] * t_0_z_0_xxy_0[i] - frz2_0[i] * t_0_z_0_xxy_1[i];

            t_0_zzz_0_xxx_0[i] += rpb_z[i] * t_0_zz_0_xxx_0[i] + rwp_z[i] * t_0_zz_0_xxx_1[i] + fz_0[i] * t_0_z_0_xxx_0[i] - frz2_0[i] * t_0_z_0_xxx_1[i];

            t_0_yzz_0_zzz_0[i] += rpb_y[i] * t_0_zz_0_zzz_0[i] + rwp_y[i] * t_0_zz_0_zzz_1[i];

            t_0_yzz_0_yzz_0[i] += rpb_y[i] * t_0_zz_0_yzz_0[i] + rwp_y[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_yzz_0_yyz_0[i] += rpb_y[i] * t_0_zz_0_yyz_0[i] + rwp_y[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_yzz_0_yyy_0[i] += rpb_y[i] * t_0_zz_0_yyy_0[i] + rwp_y[i] * t_0_zz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_yzz_0_xzz_0[i] += rpb_y[i] * t_0_zz_0_xzz_0[i] + rwp_y[i] * t_0_zz_0_xzz_1[i];

            t_0_yzz_0_xyz_0[i] += rpb_y[i] * t_0_zz_0_xyz_0[i] + rwp_y[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_yzz_0_xyy_0[i] += rpb_y[i] * t_0_zz_0_xyy_0[i] + rwp_y[i] * t_0_zz_0_xyy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_yzz_0_xxz_0[i] += rpb_y[i] * t_0_zz_0_xxz_0[i] + rwp_y[i] * t_0_zz_0_xxz_1[i];

            t_0_yzz_0_xxy_0[i] += rpb_y[i] * t_0_zz_0_xxy_0[i] + rwp_y[i] * t_0_zz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzz_0_xxx_0[i] += rpb_y[i] * t_0_zz_0_xxx_0[i] + rwp_y[i] * t_0_zz_0_xxx_1[i];

            t_0_yyz_0_zzz_0[i] += rpb_z[i] * t_0_yy_0_zzz_0[i] + rwp_z[i] * t_0_yy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyz_0_yzz_0[i] += rpb_z[i] * t_0_yy_0_yzz_0[i] + rwp_z[i] * t_0_yy_0_yzz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyz_0_yyz_0[i] += rpb_z[i] * t_0_yy_0_yyz_0[i] + rwp_z[i] * t_0_yy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyz_0_yyy_0[i] += rpb_z[i] * t_0_yy_0_yyy_0[i] + rwp_z[i] * t_0_yy_0_yyy_1[i];

            t_0_yyz_0_xzz_0[i] += rpb_z[i] * t_0_yy_0_xzz_0[i] + rwp_z[i] * t_0_yy_0_xzz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyz_0_xyz_0[i] += rpb_z[i] * t_0_yy_0_xyz_0[i] + rwp_z[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyz_0_xyy_0[i] += rpb_z[i] * t_0_yy_0_xyy_0[i] + rwp_z[i] * t_0_yy_0_xyy_1[i];

            t_0_yyz_0_xxz_0[i] += rpb_z[i] * t_0_yy_0_xxz_0[i] + rwp_z[i] * t_0_yy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyz_0_xxy_0[i] += rpb_z[i] * t_0_yy_0_xxy_0[i] + rwp_z[i] * t_0_yy_0_xxy_1[i];

            t_0_yyy_0_zzz_0[i] += rpb_y[i] * t_0_yy_0_zzz_0[i] + rwp_y[i] * t_0_yy_0_zzz_1[i] + fz_0[i] * t_0_y_0_zzz_0[i] - frz2_0[i] * t_0_y_0_zzz_1[i];

            t_0_yyy_0_yzz_0[i] += rpb_y[i] * t_0_yy_0_yzz_0[i] + rwp_y[i] * t_0_yy_0_yzz_1[i] + fz_0[i] * t_0_y_0_yzz_0[i] - frz2_0[i] * t_0_y_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyy_0_yyz_0[i] += rpb_y[i] * t_0_yy_0_yyz_0[i] + rwp_y[i] * t_0_yy_0_yyz_1[i] + fz_0[i] * t_0_y_0_yyz_0[i] - frz2_0[i] * t_0_y_0_yyz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyy_0_yyy_0[i] += rpb_y[i] * t_0_yy_0_yyy_0[i] + rwp_y[i] * t_0_yy_0_yyy_1[i] + fz_0[i] * t_0_y_0_yyy_0[i] - frz2_0[i] * t_0_y_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyy_0_xzz_0[i] += rpb_y[i] * t_0_yy_0_xzz_0[i] + rwp_y[i] * t_0_yy_0_xzz_1[i] + fz_0[i] * t_0_y_0_xzz_0[i] - frz2_0[i] * t_0_y_0_xzz_1[i];

            t_0_yyy_0_xyz_0[i] += rpb_y[i] * t_0_yy_0_xyz_0[i] + rwp_y[i] * t_0_yy_0_xyz_1[i] + fz_0[i] * t_0_y_0_xyz_0[i] - frz2_0[i] * t_0_y_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyy_0_xyy_0[i] += rpb_y[i] * t_0_yy_0_xyy_0[i] + rwp_y[i] * t_0_yy_0_xyy_1[i] + fz_0[i] * t_0_y_0_xyy_0[i] - frz2_0[i] * t_0_y_0_xyy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_yzz_0, t_0_x_0_yzz_1, t_0_x_0_zzz_0, t_0_x_0_zzz_1,\
                               t_0_xx_0_xx_1, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1, t_0_xx_0_xy_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yy_1,\
                               t_0_xx_0_yyy_0, t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1,\
                               t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1, t_0_xx_0_zz_1,\
                               t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_yzz_0, t_0_xxx_0_zzz_0,\
                               t_0_xxy_0_xxx_0, t_0_xxy_0_xxy_0, t_0_xxy_0_xxz_0, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xzz_0, t_0_xxy_0_yyy_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxz_0, t_0_xxz_0_xyy_0, t_0_xxz_0_xyz_0, t_0_xxz_0_xzz_0,\
                               t_0_xxz_0_yyz_0, t_0_xxz_0_yzz_0, t_0_xxz_0_zzz_0, t_0_xyy_0_xxx_0,\
                               t_0_xyy_0_xxy_0, t_0_xyy_0_xyy_0, t_0_xyy_0_xyz_0, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyz_0, t_0_xyy_0_yzz_0, t_0_xyy_0_zzz_0, t_0_xzz_0_xxx_0,\
                               t_0_xzz_0_xxz_0, t_0_xzz_0_xyz_0, t_0_xzz_0_xzz_0, t_0_xzz_0_yyy_0,\
                               t_0_xzz_0_yyz_0, t_0_xzz_0_yzz_0, t_0_xzz_0_zzz_0, t_0_y_0_xxx_0,\
                               t_0_y_0_xxx_1, t_0_y_0_xxy_0, t_0_y_0_xxy_1, t_0_y_0_xxz_0,\
                               t_0_y_0_xxz_1, t_0_yy_0_xx_1, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0,\
                               t_0_yy_0_xyz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1, t_0_yy_0_yzz_0,\
                               t_0_yy_0_yzz_1, t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxy_0, t_0_yyy_0_xxz_0, t_0_zz_0_xx_1, t_0_zz_0_xxx_0,\
                               t_0_zz_0_xxx_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xyz_0,\
                               t_0_zz_0_xyz_1, t_0_zz_0_xz_1, t_0_zz_0_xzz_0, t_0_zz_0_xzz_1,\
                               t_0_zz_0_yyy_0, t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1,\
                               t_0_zz_0_yz_1, t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1,\
                               t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyy_0_xxz_0[i] += rpb_y[i] * t_0_yy_0_xxz_0[i] + rwp_y[i] * t_0_yy_0_xxz_1[i] + fz_0[i] * t_0_y_0_xxz_0[i] - frz2_0[i] * t_0_y_0_xxz_1[i];

            t_0_yyy_0_xxy_0[i] += rpb_y[i] * t_0_yy_0_xxy_0[i] + rwp_y[i] * t_0_yy_0_xxy_1[i] + fz_0[i] * t_0_y_0_xxy_0[i] - frz2_0[i] * t_0_y_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyy_0_xxx_0[i] += rpb_y[i] * t_0_yy_0_xxx_0[i] + rwp_y[i] * t_0_yy_0_xxx_1[i] + fz_0[i] * t_0_y_0_xxx_0[i] - frz2_0[i] * t_0_y_0_xxx_1[i];

            t_0_xzz_0_zzz_0[i] += rpb_x[i] * t_0_zz_0_zzz_0[i] + rwp_x[i] * t_0_zz_0_zzz_1[i];

            t_0_xzz_0_yzz_0[i] += rpb_x[i] * t_0_zz_0_yzz_0[i] + rwp_x[i] * t_0_zz_0_yzz_1[i];

            t_0_xzz_0_yyz_0[i] += rpb_x[i] * t_0_zz_0_yyz_0[i] + rwp_x[i] * t_0_zz_0_yyz_1[i];

            t_0_xzz_0_yyy_0[i] += rpb_x[i] * t_0_zz_0_yyy_0[i] + rwp_x[i] * t_0_zz_0_yyy_1[i];

            t_0_xzz_0_xzz_0[i] += rpb_x[i] * t_0_zz_0_xzz_0[i] + rwp_x[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_xzz_0_xyz_0[i] += rpb_x[i] * t_0_zz_0_xyz_0[i] + rwp_x[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_xzz_0_xxz_0[i] += rpb_x[i] * t_0_zz_0_xxz_0[i] + rwp_x[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_xzz_0_xxx_0[i] += rpb_x[i] * t_0_zz_0_xxx_0[i] + rwp_x[i] * t_0_zz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_xyy_0_zzz_0[i] += rpb_x[i] * t_0_yy_0_zzz_0[i] + rwp_x[i] * t_0_yy_0_zzz_1[i];

            t_0_xyy_0_yzz_0[i] += rpb_x[i] * t_0_yy_0_yzz_0[i] + rwp_x[i] * t_0_yy_0_yzz_1[i];

            t_0_xyy_0_yyz_0[i] += rpb_x[i] * t_0_yy_0_yyz_0[i] + rwp_x[i] * t_0_yy_0_yyz_1[i];

            t_0_xyy_0_yyy_0[i] += rpb_x[i] * t_0_yy_0_yyy_0[i] + rwp_x[i] * t_0_yy_0_yyy_1[i];

            t_0_xyy_0_xyz_0[i] += rpb_x[i] * t_0_yy_0_xyz_0[i] + rwp_x[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_xyy_0_xyy_0[i] += rpb_x[i] * t_0_yy_0_xyy_0[i] + rwp_x[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_xyy_0_xxy_0[i] += rpb_x[i] * t_0_yy_0_xxy_0[i] + rwp_x[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_xyy_0_xxx_0[i] += rpb_x[i] * t_0_yy_0_xxx_0[i] + rwp_x[i] * t_0_yy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_xxz_0_zzz_0[i] += rpb_z[i] * t_0_xx_0_zzz_0[i] + rwp_z[i] * t_0_xx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxz_0_yzz_0[i] += rpb_z[i] * t_0_xx_0_yzz_0[i] + rwp_z[i] * t_0_xx_0_yzz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxz_0_yyz_0[i] += rpb_z[i] * t_0_xx_0_yyz_0[i] + rwp_z[i] * t_0_xx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxz_0_xzz_0[i] += rpb_z[i] * t_0_xx_0_xzz_0[i] + rwp_z[i] * t_0_xx_0_xzz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxz_0_xyz_0[i] += rpb_z[i] * t_0_xx_0_xyz_0[i] + rwp_z[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxz_0_xyy_0[i] += rpb_z[i] * t_0_xx_0_xyy_0[i] + rwp_z[i] * t_0_xx_0_xyy_1[i];

            t_0_xxz_0_xxz_0[i] += rpb_z[i] * t_0_xx_0_xxz_0[i] + rwp_z[i] * t_0_xx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxz_0_xxy_0[i] += rpb_z[i] * t_0_xx_0_xxy_0[i] + rwp_z[i] * t_0_xx_0_xxy_1[i];

            t_0_xxz_0_xxx_0[i] += rpb_z[i] * t_0_xx_0_xxx_0[i] + rwp_z[i] * t_0_xx_0_xxx_1[i];

            t_0_xxy_0_yyy_0[i] += rpb_y[i] * t_0_xx_0_yyy_0[i] + rwp_y[i] * t_0_xx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxy_0_xzz_0[i] += rpb_y[i] * t_0_xx_0_xzz_0[i] + rwp_y[i] * t_0_xx_0_xzz_1[i];

            t_0_xxy_0_xyy_0[i] += rpb_y[i] * t_0_xx_0_xyy_0[i] + rwp_y[i] * t_0_xx_0_xyy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxy_0_xxz_0[i] += rpb_y[i] * t_0_xx_0_xxz_0[i] + rwp_y[i] * t_0_xx_0_xxz_1[i];

            t_0_xxy_0_xxy_0[i] += rpb_y[i] * t_0_xx_0_xxy_0[i] + rwp_y[i] * t_0_xx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxy_0_xxx_0[i] += rpb_y[i] * t_0_xx_0_xxx_0[i] + rwp_y[i] * t_0_xx_0_xxx_1[i];

            t_0_xxx_0_zzz_0[i] += rpb_x[i] * t_0_xx_0_zzz_0[i] + rwp_x[i] * t_0_xx_0_zzz_1[i] + fz_0[i] * t_0_x_0_zzz_0[i] - frz2_0[i] * t_0_x_0_zzz_1[i];

            t_0_xxx_0_yzz_0[i] += rpb_x[i] * t_0_xx_0_yzz_0[i] + rwp_x[i] * t_0_xx_0_yzz_1[i] + fz_0[i] * t_0_x_0_yzz_0[i] - frz2_0[i] * t_0_x_0_yzz_1[i];
        }
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_x_0_xxx_0, t_0_x_0_xxx_1,\
                               t_0_x_0_xxy_0, t_0_x_0_xxy_1, t_0_x_0_xxz_0, t_0_x_0_xxz_1,\
                               t_0_x_0_xyy_0, t_0_x_0_xyy_1, t_0_x_0_xyz_0, t_0_x_0_xyz_1,\
                               t_0_x_0_xzz_0, t_0_x_0_xzz_1, t_0_x_0_yyy_0, t_0_x_0_yyy_1,\
                               t_0_x_0_yyz_0, t_0_x_0_yyz_1, t_0_xx_0_xx_1, t_0_xx_0_xxx_0,\
                               t_0_xx_0_xxx_1, t_0_xx_0_xxy_0, t_0_xx_0_xxy_1, t_0_xx_0_xxz_0,\
                               t_0_xx_0_xxz_1, t_0_xx_0_xy_1, t_0_xx_0_xyy_0, t_0_xx_0_xyy_1,\
                               t_0_xx_0_xyz_0, t_0_xx_0_xyz_1, t_0_xx_0_xz_1, t_0_xx_0_xzz_0,\
                               t_0_xx_0_xzz_1, t_0_xx_0_yy_1, t_0_xx_0_yyy_0, t_0_xx_0_yyy_1,\
                               t_0_xx_0_yyz_0, t_0_xx_0_yyz_1, t_0_xx_0_yz_1, t_0_xx_0_zz_1,\
                               t_0_xxx_0_xxx_0, t_0_xxx_0_xxy_0, t_0_xxx_0_xxz_0, t_0_xxx_0_xyy_0,\
                               t_0_xxx_0_xyz_0, t_0_xxx_0_xzz_0, t_0_xxx_0_yyy_0, t_0_xxx_0_yyz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxx_0_yyz_0[i] += rpb_x[i] * t_0_xx_0_yyz_0[i] + rwp_x[i] * t_0_xx_0_yyz_1[i] + fz_0[i] * t_0_x_0_yyz_0[i] - frz2_0[i] * t_0_x_0_yyz_1[i];

            t_0_xxx_0_yyy_0[i] += rpb_x[i] * t_0_xx_0_yyy_0[i] + rwp_x[i] * t_0_xx_0_yyy_1[i] + fz_0[i] * t_0_x_0_yyy_0[i] - frz2_0[i] * t_0_x_0_yyy_1[i];

            t_0_xxx_0_xzz_0[i] += rpb_x[i] * t_0_xx_0_xzz_0[i] + rwp_x[i] * t_0_xx_0_xzz_1[i] + fz_0[i] * t_0_x_0_xzz_0[i] - frz2_0[i] * t_0_x_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxx_0_xyz_0[i] += rpb_x[i] * t_0_xx_0_xyz_0[i] + rwp_x[i] * t_0_xx_0_xyz_1[i] + fz_0[i] * t_0_x_0_xyz_0[i] - frz2_0[i] * t_0_x_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxx_0_xyy_0[i] += rpb_x[i] * t_0_xx_0_xyy_0[i] + rwp_x[i] * t_0_xx_0_xyy_1[i] + fz_0[i] * t_0_x_0_xyy_0[i] - frz2_0[i] * t_0_x_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxx_0_xxz_0[i] += rpb_x[i] * t_0_xx_0_xxz_0[i] + rwp_x[i] * t_0_xx_0_xxz_1[i] + fz_0[i] * t_0_x_0_xxz_0[i] - frz2_0[i] * t_0_x_0_xxz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxx_0_xxy_0[i] += rpb_x[i] * t_0_xx_0_xxy_0[i] + rwp_x[i] * t_0_xx_0_xxy_1[i] + fz_0[i] * t_0_x_0_xxy_0[i] - frz2_0[i] * t_0_x_0_xxy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxx_0_xxx_0[i] += rpb_x[i] * t_0_xx_0_xxx_0[i] + rwp_x[i] * t_0_xx_0_xxx_1[i] + fz_0[i] * t_0_x_0_xxx_0[i] - frz2_0[i] * t_0_x_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_xx_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_y, rpb_z, rwp_y, rwp_z, t_0_y_0_xyy_0,\
                               t_0_y_0_xyy_1, t_0_y_0_xyz_0, t_0_y_0_xyz_1, t_0_y_0_xzz_0,\
                               t_0_y_0_xzz_1, t_0_y_0_yyy_0, t_0_y_0_yyy_1, t_0_y_0_yyz_0,\
                               t_0_y_0_yyz_1, t_0_y_0_yzz_0, t_0_y_0_yzz_1, t_0_y_0_zzz_0,\
                               t_0_y_0_zzz_1, t_0_yy_0_xx_1, t_0_yy_0_xxy_0, t_0_yy_0_xxy_1,\
                               t_0_yy_0_xxz_0, t_0_yy_0_xxz_1, t_0_yy_0_xy_1, t_0_yy_0_xyy_0,\
                               t_0_yy_0_xyy_1, t_0_yy_0_xyz_0, t_0_yy_0_xyz_1, t_0_yy_0_xz_1,\
                               t_0_yy_0_xzz_0, t_0_yy_0_xzz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0,\
                               t_0_yy_0_yyy_1, t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1,\
                               t_0_yy_0_yzz_0, t_0_yy_0_yzz_1, t_0_yy_0_zz_1, t_0_yy_0_zzz_0,\
                               t_0_yy_0_zzz_1, t_0_yyy_0_xyy_0, t_0_yyy_0_xyz_0, t_0_yyy_0_xzz_0,\
                               t_0_yyy_0_yyy_0, t_0_yyy_0_yyz_0, t_0_yyy_0_yzz_0, t_0_yyy_0_zzz_0,\
                               t_0_yyz_0_xxy_0, t_0_yyz_0_xxz_0, t_0_yyz_0_xyy_0, t_0_yyz_0_xyz_0,\
                               t_0_yyz_0_xzz_0, t_0_yyz_0_yyy_0, t_0_yyz_0_yyz_0, t_0_yyz_0_yzz_0,\
                               t_0_yyz_0_zzz_0, t_0_yzz_0_xxx_0, t_0_yzz_0_xxy_0, t_0_yzz_0_xxz_0,\
                               t_0_yzz_0_xyy_0, t_0_yzz_0_xyz_0, t_0_yzz_0_xzz_0, t_0_yzz_0_yyy_0,\
                               t_0_yzz_0_yyz_0, t_0_yzz_0_yzz_0, t_0_yzz_0_zzz_0, t_0_z_0_xxx_0,\
                               t_0_z_0_xxx_1, t_0_z_0_xxy_0, t_0_z_0_xxy_1, t_0_z_0_xxz_0,\
                               t_0_z_0_xxz_1, t_0_z_0_xyy_0, t_0_z_0_xyy_1, t_0_z_0_xyz_0,\
                               t_0_z_0_xyz_1, t_0_z_0_xzz_0, t_0_z_0_xzz_1, t_0_z_0_yyy_0,\
                               t_0_z_0_yyy_1, t_0_z_0_yyz_0, t_0_z_0_yyz_1, t_0_z_0_yzz_0,\
                               t_0_z_0_yzz_1, t_0_z_0_zzz_0, t_0_z_0_zzz_1, t_0_zz_0_xx_1,\
                               t_0_zz_0_xxx_0, t_0_zz_0_xxx_1, t_0_zz_0_xxy_0, t_0_zz_0_xxy_1,\
                               t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xy_1, t_0_zz_0_xyy_0,\
                               t_0_zz_0_xyy_1, t_0_zz_0_xyz_0, t_0_zz_0_xyz_1, t_0_zz_0_xz_1,\
                               t_0_zz_0_xzz_0, t_0_zz_0_xzz_1, t_0_zz_0_yy_1, t_0_zz_0_yyy_0,\
                               t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1, t_0_zz_0_yz_1,\
                               t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1, t_0_zz_0_zzz_0,\
                               t_0_zz_0_zzz_1, t_0_zzz_0_xxx_0, t_0_zzz_0_xxy_0, t_0_zzz_0_xxz_0,\
                               t_0_zzz_0_xyy_0, t_0_zzz_0_xyz_0, t_0_zzz_0_xzz_0, t_0_zzz_0_yyy_0,\
                               t_0_zzz_0_yyz_0, t_0_zzz_0_yzz_0, t_0_zzz_0_zzz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_zzz_0[i] = rpb_z[i] * t_0_zz_0_zzz_0[i] + rwp_z[i] * t_0_zz_0_zzz_1[i] + fz_0[i] * t_0_z_0_zzz_0[i] - frz2_0[i] * t_0_z_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_zzz_0_yzz_0[i] = rpb_z[i] * t_0_zz_0_yzz_0[i] + rwp_z[i] * t_0_zz_0_yzz_1[i] + fz_0[i] * t_0_z_0_yzz_0[i] - frz2_0[i] * t_0_z_0_yzz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_zzz_0_yyz_0[i] = rpb_z[i] * t_0_zz_0_yyz_0[i] + rwp_z[i] * t_0_zz_0_yyz_1[i] + fz_0[i] * t_0_z_0_yyz_0[i] - frz2_0[i] * t_0_z_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_zzz_0_yyy_0[i] = rpb_z[i] * t_0_zz_0_yyy_0[i] + rwp_z[i] * t_0_zz_0_yyy_1[i] + fz_0[i] * t_0_z_0_yyy_0[i] - frz2_0[i] * t_0_z_0_yyy_1[i];

            t_0_zzz_0_xzz_0[i] = rpb_z[i] * t_0_zz_0_xzz_0[i] + rwp_z[i] * t_0_zz_0_xzz_1[i] + fz_0[i] * t_0_z_0_xzz_0[i] - frz2_0[i] * t_0_z_0_xzz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_zzz_0_xyz_0[i] = rpb_z[i] * t_0_zz_0_xyz_0[i] + rwp_z[i] * t_0_zz_0_xyz_1[i] + fz_0[i] * t_0_z_0_xyz_0[i] - frz2_0[i] * t_0_z_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_zzz_0_xyy_0[i] = rpb_z[i] * t_0_zz_0_xyy_0[i] + rwp_z[i] * t_0_zz_0_xyy_1[i] + fz_0[i] * t_0_z_0_xyy_0[i] - frz2_0[i] * t_0_z_0_xyy_1[i];

            t_0_zzz_0_xxz_0[i] = rpb_z[i] * t_0_zz_0_xxz_0[i] + rwp_z[i] * t_0_zz_0_xxz_1[i] + fz_0[i] * t_0_z_0_xxz_0[i] - frz2_0[i] * t_0_z_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_zzz_0_xxy_0[i] = rpb_z[i] * t_0_zz_0_xxy_0[i] + rwp_z[i] * t_0_zz_0_xxy_1[i] + fz_0[i] * t_0_z_0_xxy_0[i] - frz2_0[i] * t_0_z_0_xxy_1[i];

            t_0_zzz_0_xxx_0[i] = rpb_z[i] * t_0_zz_0_xxx_0[i] + rwp_z[i] * t_0_zz_0_xxx_1[i] + fz_0[i] * t_0_z_0_xxx_0[i] - frz2_0[i] * t_0_z_0_xxx_1[i];

            t_0_yzz_0_zzz_0[i] = rpb_y[i] * t_0_zz_0_zzz_0[i] + rwp_y[i] * t_0_zz_0_zzz_1[i];

            t_0_yzz_0_yzz_0[i] = rpb_y[i] * t_0_zz_0_yzz_0[i] + rwp_y[i] * t_0_zz_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_yzz_0_yyz_0[i] = rpb_y[i] * t_0_zz_0_yyz_0[i] + rwp_y[i] * t_0_zz_0_yyz_1[i] + fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_yzz_0_yyy_0[i] = rpb_y[i] * t_0_zz_0_yyy_0[i] + rwp_y[i] * t_0_zz_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_yy_1[i];

            t_0_yzz_0_xzz_0[i] = rpb_y[i] * t_0_zz_0_xzz_0[i] + rwp_y[i] * t_0_zz_0_xzz_1[i];

            t_0_yzz_0_xyz_0[i] = rpb_y[i] * t_0_zz_0_xyz_0[i] + rwp_y[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_yzz_0_xyy_0[i] = rpb_y[i] * t_0_zz_0_xyy_0[i] + rwp_y[i] * t_0_zz_0_xyy_1[i] + fze_0[i] * t_0_zz_0_xy_1[i];

            t_0_yzz_0_xxz_0[i] = rpb_y[i] * t_0_zz_0_xxz_0[i] + rwp_y[i] * t_0_zz_0_xxz_1[i];

            t_0_yzz_0_xxy_0[i] = rpb_y[i] * t_0_zz_0_xxy_0[i] + rwp_y[i] * t_0_zz_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_yzz_0_xxx_0[i] = rpb_y[i] * t_0_zz_0_xxx_0[i] + rwp_y[i] * t_0_zz_0_xxx_1[i];

            t_0_yyz_0_zzz_0[i] = rpb_z[i] * t_0_yy_0_zzz_0[i] + rwp_z[i] * t_0_yy_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyz_0_yzz_0[i] = rpb_z[i] * t_0_yy_0_yzz_0[i] + rwp_z[i] * t_0_yy_0_yzz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyz_0_yyz_0[i] = rpb_z[i] * t_0_yy_0_yyz_0[i] + rwp_z[i] * t_0_yy_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyz_0_yyy_0[i] = rpb_z[i] * t_0_yy_0_yyy_0[i] + rwp_z[i] * t_0_yy_0_yyy_1[i];

            t_0_yyz_0_xzz_0[i] = rpb_z[i] * t_0_yy_0_xzz_0[i] + rwp_z[i] * t_0_yy_0_xzz_1[i] + fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyz_0_xyz_0[i] = rpb_z[i] * t_0_yy_0_xyz_0[i] + rwp_z[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_yyz_0_xyy_0[i] = rpb_z[i] * t_0_yy_0_xyy_0[i] + rwp_z[i] * t_0_yy_0_xyy_1[i];

            t_0_yyz_0_xxz_0[i] = rpb_z[i] * t_0_yy_0_xxz_0[i] + rwp_z[i] * t_0_yy_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyz_0_xxy_0[i] = rpb_z[i] * t_0_yy_0_xxy_0[i] + rwp_z[i] * t_0_yy_0_xxy_1[i];

            t_0_yyy_0_zzz_0[i] = rpb_y[i] * t_0_yy_0_zzz_0[i] + rwp_y[i] * t_0_yy_0_zzz_1[i] + fz_0[i] * t_0_y_0_zzz_0[i] - frz2_0[i] * t_0_y_0_zzz_1[i];

            t_0_yyy_0_yzz_0[i] = rpb_y[i] * t_0_yy_0_yzz_0[i] + rwp_y[i] * t_0_yy_0_yzz_1[i] + fz_0[i] * t_0_y_0_yzz_0[i] - frz2_0[i] * t_0_y_0_yzz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_zz_1[i];

            t_0_yyy_0_yyz_0[i] = rpb_y[i] * t_0_yy_0_yyz_0[i] + rwp_y[i] * t_0_yy_0_yyz_1[i] + fz_0[i] * t_0_y_0_yyz_0[i] - frz2_0[i] * t_0_y_0_yyz_1[i] + fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_yyy_0_yyy_0[i] = rpb_y[i] * t_0_yy_0_yyy_0[i] + rwp_y[i] * t_0_yy_0_yyy_1[i] + fz_0[i] * t_0_y_0_yyy_0[i] - frz2_0[i] * t_0_y_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_yyy_0_xzz_0[i] = rpb_y[i] * t_0_yy_0_xzz_0[i] + rwp_y[i] * t_0_yy_0_xzz_1[i] + fz_0[i] * t_0_y_0_xzz_0[i] - frz2_0[i] * t_0_y_0_xzz_1[i];

            t_0_yyy_0_xyz_0[i] = rpb_y[i] * t_0_yy_0_xyz_0[i] + rwp_y[i] * t_0_yy_0_xyz_1[i] + fz_0[i] * t_0_y_0_xyz_0[i] - frz2_0[i] * t_0_y_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xz_1[i];

            t_0_yyy_0_xyy_0[i] = rpb_y[i] * t_0_yy_0_xyy_0[i] + rwp_y[i] * t_0_yy_0_xyy_1[i] + fz_0[i] * t_0_y_0_xyy_0[i] - frz2_0[i] * t_0_y_0_xyy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_yzz_0, t_0_x_0_yzz_1, t_0_x_0_zzz_0, t_0_x_0_zzz_1,\
                               t_0_xx_0_xx_1, t_0_xx_0_xxx_0, t_0_xx_0_xxx_1, t_0_xx_0_xxy_0,\
                               t_0_xx_0_xxy_1, t_0_xx_0_xxz_0, t_0_xx_0_xxz_1, t_0_xx_0_xy_1,\
                               t_0_xx_0_xyy_0, t_0_xx_0_xyy_1, t_0_xx_0_xyz_0, t_0_xx_0_xyz_1,\
                               t_0_xx_0_xz_1, t_0_xx_0_xzz_0, t_0_xx_0_xzz_1, t_0_xx_0_yy_1,\
                               t_0_xx_0_yyy_0, t_0_xx_0_yyy_1, t_0_xx_0_yyz_0, t_0_xx_0_yyz_1,\
                               t_0_xx_0_yz_1, t_0_xx_0_yzz_0, t_0_xx_0_yzz_1, t_0_xx_0_zz_1,\
                               t_0_xx_0_zzz_0, t_0_xx_0_zzz_1, t_0_xxx_0_yzz_0, t_0_xxx_0_zzz_0,\
                               t_0_xxy_0_xxx_0, t_0_xxy_0_xxy_0, t_0_xxy_0_xxz_0, t_0_xxy_0_xyy_0,\
                               t_0_xxy_0_xzz_0, t_0_xxy_0_yyy_0, t_0_xxz_0_xxx_0, t_0_xxz_0_xxy_0,\
                               t_0_xxz_0_xxz_0, t_0_xxz_0_xyy_0, t_0_xxz_0_xyz_0, t_0_xxz_0_xzz_0,\
                               t_0_xxz_0_yyz_0, t_0_xxz_0_yzz_0, t_0_xxz_0_zzz_0, t_0_xyy_0_xxx_0,\
                               t_0_xyy_0_xxy_0, t_0_xyy_0_xyy_0, t_0_xyy_0_xyz_0, t_0_xyy_0_yyy_0,\
                               t_0_xyy_0_yyz_0, t_0_xyy_0_yzz_0, t_0_xyy_0_zzz_0, t_0_xzz_0_xxx_0,\
                               t_0_xzz_0_xxz_0, t_0_xzz_0_xyz_0, t_0_xzz_0_xzz_0, t_0_xzz_0_yyy_0,\
                               t_0_xzz_0_yyz_0, t_0_xzz_0_yzz_0, t_0_xzz_0_zzz_0, t_0_y_0_xxx_0,\
                               t_0_y_0_xxx_1, t_0_y_0_xxy_0, t_0_y_0_xxy_1, t_0_y_0_xxz_0,\
                               t_0_y_0_xxz_1, t_0_yy_0_xx_1, t_0_yy_0_xxx_0, t_0_yy_0_xxx_1,\
                               t_0_yy_0_xxy_0, t_0_yy_0_xxy_1, t_0_yy_0_xxz_0, t_0_yy_0_xxz_1,\
                               t_0_yy_0_xy_1, t_0_yy_0_xyy_0, t_0_yy_0_xyy_1, t_0_yy_0_xyz_0,\
                               t_0_yy_0_xyz_1, t_0_yy_0_yy_1, t_0_yy_0_yyy_0, t_0_yy_0_yyy_1,\
                               t_0_yy_0_yyz_0, t_0_yy_0_yyz_1, t_0_yy_0_yz_1, t_0_yy_0_yzz_0,\
                               t_0_yy_0_yzz_1, t_0_yy_0_zzz_0, t_0_yy_0_zzz_1, t_0_yyy_0_xxx_0,\
                               t_0_yyy_0_xxy_0, t_0_yyy_0_xxz_0, t_0_zz_0_xx_1, t_0_zz_0_xxx_0,\
                               t_0_zz_0_xxx_1, t_0_zz_0_xxz_0, t_0_zz_0_xxz_1, t_0_zz_0_xyz_0,\
                               t_0_zz_0_xyz_1, t_0_zz_0_xz_1, t_0_zz_0_xzz_0, t_0_zz_0_xzz_1,\
                               t_0_zz_0_yyy_0, t_0_zz_0_yyy_1, t_0_zz_0_yyz_0, t_0_zz_0_yyz_1,\
                               t_0_zz_0_yz_1, t_0_zz_0_yzz_0, t_0_zz_0_yzz_1, t_0_zz_0_zz_1,\
                               t_0_zz_0_zzz_0, t_0_zz_0_zzz_1 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_yyy_0_xxz_0[i] = rpb_y[i] * t_0_yy_0_xxz_0[i] + rwp_y[i] * t_0_yy_0_xxz_1[i] + fz_0[i] * t_0_y_0_xxz_0[i] - frz2_0[i] * t_0_y_0_xxz_1[i];

            t_0_yyy_0_xxy_0[i] = rpb_y[i] * t_0_yy_0_xxy_0[i] + rwp_y[i] * t_0_yy_0_xxy_1[i] + fz_0[i] * t_0_y_0_xxy_0[i] - frz2_0[i] * t_0_y_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_yyy_0_xxx_0[i] = rpb_y[i] * t_0_yy_0_xxx_0[i] + rwp_y[i] * t_0_yy_0_xxx_1[i] + fz_0[i] * t_0_y_0_xxx_0[i] - frz2_0[i] * t_0_y_0_xxx_1[i];

            t_0_xzz_0_zzz_0[i] = rpb_x[i] * t_0_zz_0_zzz_0[i] + rwp_x[i] * t_0_zz_0_zzz_1[i];

            t_0_xzz_0_yzz_0[i] = rpb_x[i] * t_0_zz_0_yzz_0[i] + rwp_x[i] * t_0_zz_0_yzz_1[i];

            t_0_xzz_0_yyz_0[i] = rpb_x[i] * t_0_zz_0_yyz_0[i] + rwp_x[i] * t_0_zz_0_yyz_1[i];

            t_0_xzz_0_yyy_0[i] = rpb_x[i] * t_0_zz_0_yyy_0[i] + rwp_x[i] * t_0_zz_0_yyy_1[i];

            t_0_xzz_0_xzz_0[i] = rpb_x[i] * t_0_zz_0_xzz_0[i] + rwp_x[i] * t_0_zz_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_zz_1[i];

            t_0_xzz_0_xyz_0[i] = rpb_x[i] * t_0_zz_0_xyz_0[i] + rwp_x[i] * t_0_zz_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_yz_1[i];

            t_0_xzz_0_xxz_0[i] = rpb_x[i] * t_0_zz_0_xxz_0[i] + rwp_x[i] * t_0_zz_0_xxz_1[i] + fze_0[i] * t_0_zz_0_xz_1[i];

            t_0_xzz_0_xxx_0[i] = rpb_x[i] * t_0_zz_0_xxx_0[i] + rwp_x[i] * t_0_zz_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_zz_0_xx_1[i];

            t_0_xyy_0_zzz_0[i] = rpb_x[i] * t_0_yy_0_zzz_0[i] + rwp_x[i] * t_0_yy_0_zzz_1[i];

            t_0_xyy_0_yzz_0[i] = rpb_x[i] * t_0_yy_0_yzz_0[i] + rwp_x[i] * t_0_yy_0_yzz_1[i];

            t_0_xyy_0_yyz_0[i] = rpb_x[i] * t_0_yy_0_yyz_0[i] + rwp_x[i] * t_0_yy_0_yyz_1[i];

            t_0_xyy_0_yyy_0[i] = rpb_x[i] * t_0_yy_0_yyy_0[i] + rwp_x[i] * t_0_yy_0_yyy_1[i];

            t_0_xyy_0_xyz_0[i] = rpb_x[i] * t_0_yy_0_xyz_0[i] + rwp_x[i] * t_0_yy_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yz_1[i];

            t_0_xyy_0_xyy_0[i] = rpb_x[i] * t_0_yy_0_xyy_0[i] + rwp_x[i] * t_0_yy_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_yy_1[i];

            t_0_xyy_0_xxy_0[i] = rpb_x[i] * t_0_yy_0_xxy_0[i] + rwp_x[i] * t_0_yy_0_xxy_1[i] + fze_0[i] * t_0_yy_0_xy_1[i];

            t_0_xyy_0_xxx_0[i] = rpb_x[i] * t_0_yy_0_xxx_0[i] + rwp_x[i] * t_0_yy_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_yy_0_xx_1[i];

            t_0_xxz_0_zzz_0[i] = rpb_z[i] * t_0_xx_0_zzz_0[i] + rwp_z[i] * t_0_xx_0_zzz_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxz_0_yzz_0[i] = rpb_z[i] * t_0_xx_0_yzz_0[i] + rwp_z[i] * t_0_xx_0_yzz_1[i] + fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxz_0_yyz_0[i] = rpb_z[i] * t_0_xx_0_yyz_0[i] + rwp_z[i] * t_0_xx_0_yyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxz_0_xzz_0[i] = rpb_z[i] * t_0_xx_0_xzz_0[i] + rwp_z[i] * t_0_xx_0_xzz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxz_0_xyz_0[i] = rpb_z[i] * t_0_xx_0_xyz_0[i] + rwp_z[i] * t_0_xx_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxz_0_xyy_0[i] = rpb_z[i] * t_0_xx_0_xyy_0[i] + rwp_z[i] * t_0_xx_0_xyy_1[i];

            t_0_xxz_0_xxz_0[i] = rpb_z[i] * t_0_xx_0_xxz_0[i] + rwp_z[i] * t_0_xx_0_xxz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxz_0_xxy_0[i] = rpb_z[i] * t_0_xx_0_xxy_0[i] + rwp_z[i] * t_0_xx_0_xxy_1[i];

            t_0_xxz_0_xxx_0[i] = rpb_z[i] * t_0_xx_0_xxx_0[i] + rwp_z[i] * t_0_xx_0_xxx_1[i];

            t_0_xxy_0_yyy_0[i] = rpb_y[i] * t_0_xx_0_yyy_0[i] + rwp_y[i] * t_0_xx_0_yyy_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxy_0_xzz_0[i] = rpb_y[i] * t_0_xx_0_xzz_0[i] + rwp_y[i] * t_0_xx_0_xzz_1[i];

            t_0_xxy_0_xyy_0[i] = rpb_y[i] * t_0_xx_0_xyy_0[i] + rwp_y[i] * t_0_xx_0_xyy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxy_0_xxz_0[i] = rpb_y[i] * t_0_xx_0_xxz_0[i] + rwp_y[i] * t_0_xx_0_xxz_1[i];

            t_0_xxy_0_xxy_0[i] = rpb_y[i] * t_0_xx_0_xxy_0[i] + rwp_y[i] * t_0_xx_0_xxy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_xx_1[i];

            t_0_xxy_0_xxx_0[i] = rpb_y[i] * t_0_xx_0_xxx_0[i] + rwp_y[i] * t_0_xx_0_xxx_1[i];

            t_0_xxx_0_zzz_0[i] = rpb_x[i] * t_0_xx_0_zzz_0[i] + rwp_x[i] * t_0_xx_0_zzz_1[i] + fz_0[i] * t_0_x_0_zzz_0[i] - frz2_0[i] * t_0_x_0_zzz_1[i];

            t_0_xxx_0_yzz_0[i] = rpb_x[i] * t_0_xx_0_yzz_0[i] + rwp_x[i] * t_0_xx_0_yzz_1[i] + fz_0[i] * t_0_x_0_yzz_0[i] - frz2_0[i] * t_0_x_0_yzz_1[i];
        }

        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rwp_x, t_0_x_0_xxx_0, t_0_x_0_xxx_1,\
                               t_0_x_0_xxy_0, t_0_x_0_xxy_1, t_0_x_0_xxz_0, t_0_x_0_xxz_1,\
                               t_0_x_0_xyy_0, t_0_x_0_xyy_1, t_0_x_0_xyz_0, t_0_x_0_xyz_1,\
                               t_0_x_0_xzz_0, t_0_x_0_xzz_1, t_0_x_0_yyy_0, t_0_x_0_yyy_1,\
                               t_0_x_0_yyz_0, t_0_x_0_yyz_1, t_0_xx_0_xx_1, t_0_xx_0_xxx_0,\
                               t_0_xx_0_xxx_1, t_0_xx_0_xxy_0, t_0_xx_0_xxy_1, t_0_xx_0_xxz_0,\
                               t_0_xx_0_xxz_1, t_0_xx_0_xy_1, t_0_xx_0_xyy_0, t_0_xx_0_xyy_1,\
                               t_0_xx_0_xyz_0, t_0_xx_0_xyz_1, t_0_xx_0_xz_1, t_0_xx_0_xzz_0,\
                               t_0_xx_0_xzz_1, t_0_xx_0_yy_1, t_0_xx_0_yyy_0, t_0_xx_0_yyy_1,\
                               t_0_xx_0_yyz_0, t_0_xx_0_yyz_1, t_0_xx_0_yz_1, t_0_xx_0_zz_1,\
                               t_0_xxx_0_xxx_0, t_0_xxx_0_xxy_0, t_0_xxx_0_xxz_0, t_0_xxx_0_xyy_0,\
                               t_0_xxx_0_xyz_0, t_0_xxx_0_xzz_0, t_0_xxx_0_yyy_0, t_0_xxx_0_yyz_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_xxx_0_yyz_0[i] = rpb_x[i] * t_0_xx_0_yyz_0[i] + rwp_x[i] * t_0_xx_0_yyz_1[i] + fz_0[i] * t_0_x_0_yyz_0[i] - frz2_0[i] * t_0_x_0_yyz_1[i];

            t_0_xxx_0_yyy_0[i] = rpb_x[i] * t_0_xx_0_yyy_0[i] + rwp_x[i] * t_0_xx_0_yyy_1[i] + fz_0[i] * t_0_x_0_yyy_0[i] - frz2_0[i] * t_0_x_0_yyy_1[i];

            t_0_xxx_0_xzz_0[i] = rpb_x[i] * t_0_xx_0_xzz_0[i] + rwp_x[i] * t_0_xx_0_xzz_1[i] + fz_0[i] * t_0_x_0_xzz_0[i] - frz2_0[i] * t_0_x_0_xzz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_zz_1[i];

            t_0_xxx_0_xyz_0[i] = rpb_x[i] * t_0_xx_0_xyz_0[i] + rwp_x[i] * t_0_xx_0_xyz_1[i] + fz_0[i] * t_0_x_0_xyz_0[i] - frz2_0[i] * t_0_x_0_xyz_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yz_1[i];

            t_0_xxx_0_xyy_0[i] = rpb_x[i] * t_0_xx_0_xyy_0[i] + rwp_x[i] * t_0_xx_0_xyy_1[i] + fz_0[i] * t_0_x_0_xyy_0[i] - frz2_0[i] * t_0_x_0_xyy_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_yy_1[i];

            t_0_xxx_0_xxz_0[i] = rpb_x[i] * t_0_xx_0_xxz_0[i] + rwp_x[i] * t_0_xx_0_xxz_1[i] + fz_0[i] * t_0_x_0_xxz_0[i] - frz2_0[i] * t_0_x_0_xxz_1[i] + fze_0[i] * t_0_xx_0_xz_1[i];

            t_0_xxx_0_xxy_0[i] = rpb_x[i] * t_0_xx_0_xxy_0[i] + rwp_x[i] * t_0_xx_0_xxy_1[i] + fz_0[i] * t_0_x_0_xxy_0[i] - frz2_0[i] * t_0_x_0_xxy_1[i] + fze_0[i] * t_0_xx_0_xy_1[i];

            t_0_xxx_0_xxx_0[i] = rpb_x[i] * t_0_xx_0_xxx_0[i] + rwp_x[i] * t_0_xx_0_xxx_1[i] + fz_0[i] * t_0_x_0_xxx_0[i] - frz2_0[i] * t_0_x_0_xxx_1[i] + fact_3_2 * fze_0[i] * t_0_xx_0_xx_1[i];
        }
    }
}


} // derirec namespace
