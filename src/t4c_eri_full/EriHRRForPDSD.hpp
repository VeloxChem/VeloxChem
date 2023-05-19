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
compHostHRRForPDSD_V0(      BufferHostXY<T>&      intsBufferPDSD,
                      const BufferHostX<int32_t>& intsIndexesPDSD,
                      const BufferHostXY<T>&      intsBufferSDSD,
                      const BufferHostX<int32_t>& intsIndexesSDSD,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PDSD) integral components

    t_z_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(0));

    t_z_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(1));

    t_z_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(2));

    t_z_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(3));

    t_z_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(4));

    t_z_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(5));

    t_z_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(6));

    t_z_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(7));

    t_z_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(8));

    t_z_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(9));

    t_z_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(10));

    t_z_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(11));

    t_z_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(12));

    t_z_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(13));

    t_z_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(14));

    t_z_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(15));

    t_z_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(16));

    t_z_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(17));

    t_z_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(18));

    t_z_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(19));

    t_z_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(20));

    t_z_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(21));

    t_z_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(22));

    t_z_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(23));

    t_z_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(24));

    t_z_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(25));

    t_z_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(26));

    t_z_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(27));

    t_z_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(28));

    t_z_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(29));

    t_z_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(30));

    t_z_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(31));

    t_z_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(32));

    t_z_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(33));

    t_z_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(34));

    t_z_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(35));

    t_y_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(36));

    t_y_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(37));

    t_y_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(38));

    t_y_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(39));

    t_y_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(40));

    t_y_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(41));

    t_y_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(42));

    t_y_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(43));

    t_y_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(44));

    t_y_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(45));

    t_y_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(46));

    t_y_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(47));

    t_y_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(48));

    t_y_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(49));

    t_y_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(50));

    t_y_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(51));

    t_y_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(52));

    t_y_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(53));

    t_y_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(54));

    t_y_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(55));

    t_y_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(56));

    t_y_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(57));

    t_y_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(58));

    t_y_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(59));

    t_y_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(60));

    t_y_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(61));

    t_y_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(62));

    t_y_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(63));

    t_y_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(64));

    t_y_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(65));

    t_y_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(66));

    t_y_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(67));

    t_y_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(68));

    t_y_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(69));

    t_y_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(70));

    t_y_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(71));

    t_x_zz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(72));

    t_x_zz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(73));

    t_x_zz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(74));

    t_x_zz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(75));

    t_x_zz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(76));

    t_x_zz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(77));

    t_x_yz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(78));

    t_x_yz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(79));

    t_x_yz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(80));

    t_x_yz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(81));

    t_x_yz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(82));

    t_x_yz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(83));

    t_x_yy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(84));

    t_x_yy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(85));

    t_x_yy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(86));

    t_x_yy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(87));

    t_x_yy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(88));

    t_x_yy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(89));

    t_x_xz_0_zz = intsBufferPDSD.data(intsIndexesPDSD(90));

    t_x_xz_0_yz = intsBufferPDSD.data(intsIndexesPDSD(91));

    t_x_xz_0_yy = intsBufferPDSD.data(intsIndexesPDSD(92));

    t_x_xz_0_xz = intsBufferPDSD.data(intsIndexesPDSD(93));

    t_x_xz_0_xy = intsBufferPDSD.data(intsIndexesPDSD(94));

    t_x_xz_0_xx = intsBufferPDSD.data(intsIndexesPDSD(95));

    t_x_xy_0_zz = intsBufferPDSD.data(intsIndexesPDSD(96));

    t_x_xy_0_yz = intsBufferPDSD.data(intsIndexesPDSD(97));

    t_x_xy_0_yy = intsBufferPDSD.data(intsIndexesPDSD(98));

    t_x_xy_0_xz = intsBufferPDSD.data(intsIndexesPDSD(99));

    t_x_xy_0_xy = intsBufferPDSD.data(intsIndexesPDSD(100));

    t_x_xy_0_xx = intsBufferPDSD.data(intsIndexesPDSD(101));

    t_x_xx_0_zz = intsBufferPDSD.data(intsIndexesPDSD(102));

    t_x_xx_0_yz = intsBufferPDSD.data(intsIndexesPDSD(103));

    t_x_xx_0_yy = intsBufferPDSD.data(intsIndexesPDSD(104));

    t_x_xx_0_xz = intsBufferPDSD.data(intsIndexesPDSD(105));

    t_x_xx_0_xy = intsBufferPDSD.data(intsIndexesPDSD(106));

    t_x_xx_0_xx = intsBufferPDSD.data(intsIndexesPDSD(107));

    // set up (SDSD) integral components

    t_0_zz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(0));

    t_0_zz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(1));

    t_0_zz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(2));

    t_0_zz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(3));

    t_0_zz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(4));

    t_0_zz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(5));

    t_0_yz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(6));

    t_0_yz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(7));

    t_0_yz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(8));

    t_0_yz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(9));

    t_0_yz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(10));

    t_0_yz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(11));

    t_0_yy_0_zz = intsBufferSDSD.data(intsIndexesSDSD(12));

    t_0_yy_0_yz = intsBufferSDSD.data(intsIndexesSDSD(13));

    t_0_yy_0_yy = intsBufferSDSD.data(intsIndexesSDSD(14));

    t_0_yy_0_xz = intsBufferSDSD.data(intsIndexesSDSD(15));

    t_0_yy_0_xy = intsBufferSDSD.data(intsIndexesSDSD(16));

    t_0_yy_0_xx = intsBufferSDSD.data(intsIndexesSDSD(17));

    t_0_xz_0_zz = intsBufferSDSD.data(intsIndexesSDSD(18));

    t_0_xz_0_yz = intsBufferSDSD.data(intsIndexesSDSD(19));

    t_0_xz_0_yy = intsBufferSDSD.data(intsIndexesSDSD(20));

    t_0_xz_0_xz = intsBufferSDSD.data(intsIndexesSDSD(21));

    t_0_xz_0_xy = intsBufferSDSD.data(intsIndexesSDSD(22));

    t_0_xz_0_xx = intsBufferSDSD.data(intsIndexesSDSD(23));

    t_0_xy_0_zz = intsBufferSDSD.data(intsIndexesSDSD(24));

    t_0_xy_0_yz = intsBufferSDSD.data(intsIndexesSDSD(25));

    t_0_xy_0_yy = intsBufferSDSD.data(intsIndexesSDSD(26));

    t_0_xy_0_xz = intsBufferSDSD.data(intsIndexesSDSD(27));

    t_0_xy_0_xy = intsBufferSDSD.data(intsIndexesSDSD(28));

    t_0_xy_0_xx = intsBufferSDSD.data(intsIndexesSDSD(29));

    t_0_xx_0_zz = intsBufferSDSD.data(intsIndexesSDSD(30));

    t_0_xx_0_yz = intsBufferSDSD.data(intsIndexesSDSD(31));

    t_0_xx_0_yy = intsBufferSDSD.data(intsIndexesSDSD(32));

    t_0_xx_0_xz = intsBufferSDSD.data(intsIndexesSDSD(33));

    t_0_xx_0_xy = intsBufferSDSD.data(intsIndexesSDSD(34));

    t_0_xx_0_xx = intsBufferSDSD.data(intsIndexesSDSD(35));

    // set up (SFSD) integral components

    t_0_zzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(0));

    t_0_zzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(1));

    t_0_zzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(2));

    t_0_zzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(3));

    t_0_zzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(4));

    t_0_zzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(5));

    t_0_yzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(6));

    t_0_yzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(7));

    t_0_yzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(8));

    t_0_yzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(9));

    t_0_yzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(10));

    t_0_yzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(11));

    t_0_yyz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(12));

    t_0_yyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(13));

    t_0_yyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(14));

    t_0_yyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(15));

    t_0_yyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(16));

    t_0_yyz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(17));

    t_0_yyy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(18));

    t_0_yyy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(19));

    t_0_yyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(20));

    t_0_yyy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(21));

    t_0_yyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(22));

    t_0_yyy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(23));

    t_0_xzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(24));

    t_0_xzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(25));

    t_0_xzz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(26));

    t_0_xzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(27));

    t_0_xzz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(28));

    t_0_xzz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(29));

    t_0_xyz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(30));

    t_0_xyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(31));

    t_0_xyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(32));

    t_0_xyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(33));

    t_0_xyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(34));

    t_0_xyz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(35));

    t_0_xyy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(36));

    t_0_xyy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(37));

    t_0_xyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(38));

    t_0_xyy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(39));

    t_0_xyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(40));

    t_0_xyy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(41));

    t_0_xxz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(42));

    t_0_xxz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(43));

    t_0_xxz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(44));

    t_0_xxz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(45));

    t_0_xxz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(46));

    t_0_xxz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(47));

    t_0_xxy_0_zz = intsBufferSFSD.data(intsIndexesSFSD(48));

    t_0_xxy_0_yz = intsBufferSFSD.data(intsIndexesSFSD(49));

    t_0_xxy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(50));

    t_0_xxy_0_xz = intsBufferSFSD.data(intsIndexesSFSD(51));

    t_0_xxy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(52));

    t_0_xxy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(53));

    t_0_xxx_0_zz = intsBufferSFSD.data(intsIndexesSFSD(54));

    t_0_xxx_0_yz = intsBufferSFSD.data(intsIndexesSFSD(55));

    t_0_xxx_0_yy = intsBufferSFSD.data(intsIndexesSFSD(56));

    t_0_xxx_0_xz = intsBufferSFSD.data(intsIndexesSFSD(57));

    t_0_xxx_0_xy = intsBufferSFSD.data(intsIndexesSFSD(58));

    t_0_xxx_0_xx = intsBufferSFSD.data(intsIndexesSFSD(59));

    #pragma omp simd align(rab_z, t_0_xx_0_xx, t_0_xx_0_xy, t_0_xx_0_xz, t_0_xx_0_yy,\
                           t_0_xx_0_yz, t_0_xx_0_zz, t_0_xxz_0_xx, t_0_xxz_0_xy, t_0_xxz_0_xz,\
                           t_0_xxz_0_yy, t_0_xxz_0_yz, t_0_xxz_0_zz, t_0_xy_0_xx, t_0_xy_0_xy,\
                           t_0_xy_0_xz, t_0_xy_0_yy, t_0_xy_0_yz, t_0_xy_0_zz, t_0_xyz_0_xx,\
                           t_0_xyz_0_xy, t_0_xyz_0_xz, t_0_xyz_0_yy, t_0_xyz_0_yz, t_0_xyz_0_zz,\
                           t_0_xz_0_xx, t_0_xz_0_xy, t_0_xz_0_xz, t_0_xz_0_yy, t_0_xz_0_yz,\
                           t_0_xz_0_zz, t_0_xzz_0_xx, t_0_xzz_0_xy, t_0_xzz_0_xz, t_0_xzz_0_yy,\
                           t_0_xzz_0_yz, t_0_xzz_0_zz, t_0_yy_0_xx, t_0_yy_0_xy, t_0_yy_0_xz,\
                           t_0_yy_0_yy, t_0_yy_0_yz, t_0_yy_0_zz, t_0_yyz_0_xx, t_0_yyz_0_xy,\
                           t_0_yyz_0_xz, t_0_yyz_0_yy, t_0_yyz_0_yz, t_0_yyz_0_zz, t_0_yz_0_xx,\
                           t_0_yz_0_xy, t_0_yz_0_xz, t_0_yz_0_yy, t_0_yz_0_yz, t_0_yz_0_zz,\
                           t_0_yzz_0_xx, t_0_yzz_0_xy, t_0_yzz_0_xz, t_0_yzz_0_yy, t_0_yzz_0_yz,\
                           t_0_yzz_0_zz, t_0_zz_0_xx, t_0_zz_0_xy, t_0_zz_0_xz, t_0_zz_0_yy,\
                           t_0_zz_0_yz, t_0_zz_0_zz, t_0_zzz_0_xx, t_0_zzz_0_xy, t_0_zzz_0_xz,\
                           t_0_zzz_0_yy, t_0_zzz_0_yz, t_0_zzz_0_zz, t_z_xx_0_xx, t_z_xx_0_xy,\
                           t_z_xx_0_xz, t_z_xx_0_yy, t_z_xx_0_yz, t_z_xx_0_zz, t_z_xy_0_xx,\
                           t_z_xy_0_xy, t_z_xy_0_xz, t_z_xy_0_yy, t_z_xy_0_yz, t_z_xy_0_zz,\
                           t_z_xz_0_xx, t_z_xz_0_xy, t_z_xz_0_xz, t_z_xz_0_yy, t_z_xz_0_yz,\
                           t_z_xz_0_zz, t_z_yy_0_xx, t_z_yy_0_xy, t_z_yy_0_xz, t_z_yy_0_yy,\
                           t_z_yy_0_yz, t_z_yy_0_zz, t_z_yz_0_xx, t_z_yz_0_xy, t_z_yz_0_xz,\
                           t_z_yz_0_yy, t_z_yz_0_yz, t_z_yz_0_zz, t_z_zz_0_xx, t_z_zz_0_xy,\
                           t_z_zz_0_xz, t_z_zz_0_yy, t_z_zz_0_yz, t_z_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_0_zz[i] = t_0_zzz_0_zz[i] - rab_z[i] * t_0_zz_0_zz[i];

        t_z_zz_0_yz[i] = t_0_zzz_0_yz[i] - rab_z[i] * t_0_zz_0_yz[i];

        t_z_zz_0_yy[i] = t_0_zzz_0_yy[i] - rab_z[i] * t_0_zz_0_yy[i];

        t_z_zz_0_xz[i] = t_0_zzz_0_xz[i] - rab_z[i] * t_0_zz_0_xz[i];

        t_z_zz_0_xy[i] = t_0_zzz_0_xy[i] - rab_z[i] * t_0_zz_0_xy[i];

        t_z_zz_0_xx[i] = t_0_zzz_0_xx[i] - rab_z[i] * t_0_zz_0_xx[i];

        t_z_yz_0_zz[i] = t_0_yzz_0_zz[i] - rab_z[i] * t_0_yz_0_zz[i];

        t_z_yz_0_yz[i] = t_0_yzz_0_yz[i] - rab_z[i] * t_0_yz_0_yz[i];

        t_z_yz_0_yy[i] = t_0_yzz_0_yy[i] - rab_z[i] * t_0_yz_0_yy[i];

        t_z_yz_0_xz[i] = t_0_yzz_0_xz[i] - rab_z[i] * t_0_yz_0_xz[i];

        t_z_yz_0_xy[i] = t_0_yzz_0_xy[i] - rab_z[i] * t_0_yz_0_xy[i];

        t_z_yz_0_xx[i] = t_0_yzz_0_xx[i] - rab_z[i] * t_0_yz_0_xx[i];

        t_z_yy_0_zz[i] = t_0_yyz_0_zz[i] - rab_z[i] * t_0_yy_0_zz[i];

        t_z_yy_0_yz[i] = t_0_yyz_0_yz[i] - rab_z[i] * t_0_yy_0_yz[i];

        t_z_yy_0_yy[i] = t_0_yyz_0_yy[i] - rab_z[i] * t_0_yy_0_yy[i];

        t_z_yy_0_xz[i] = t_0_yyz_0_xz[i] - rab_z[i] * t_0_yy_0_xz[i];

        t_z_yy_0_xy[i] = t_0_yyz_0_xy[i] - rab_z[i] * t_0_yy_0_xy[i];

        t_z_yy_0_xx[i] = t_0_yyz_0_xx[i] - rab_z[i] * t_0_yy_0_xx[i];

        t_z_xz_0_zz[i] = t_0_xzz_0_zz[i] - rab_z[i] * t_0_xz_0_zz[i];

        t_z_xz_0_yz[i] = t_0_xzz_0_yz[i] - rab_z[i] * t_0_xz_0_yz[i];

        t_z_xz_0_yy[i] = t_0_xzz_0_yy[i] - rab_z[i] * t_0_xz_0_yy[i];

        t_z_xz_0_xz[i] = t_0_xzz_0_xz[i] - rab_z[i] * t_0_xz_0_xz[i];

        t_z_xz_0_xy[i] = t_0_xzz_0_xy[i] - rab_z[i] * t_0_xz_0_xy[i];

        t_z_xz_0_xx[i] = t_0_xzz_0_xx[i] - rab_z[i] * t_0_xz_0_xx[i];

        t_z_xy_0_zz[i] = t_0_xyz_0_zz[i] - rab_z[i] * t_0_xy_0_zz[i];

        t_z_xy_0_yz[i] = t_0_xyz_0_yz[i] - rab_z[i] * t_0_xy_0_yz[i];

        t_z_xy_0_yy[i] = t_0_xyz_0_yy[i] - rab_z[i] * t_0_xy_0_yy[i];

        t_z_xy_0_xz[i] = t_0_xyz_0_xz[i] - rab_z[i] * t_0_xy_0_xz[i];

        t_z_xy_0_xy[i] = t_0_xyz_0_xy[i] - rab_z[i] * t_0_xy_0_xy[i];

        t_z_xy_0_xx[i] = t_0_xyz_0_xx[i] - rab_z[i] * t_0_xy_0_xx[i];

        t_z_xx_0_zz[i] = t_0_xxz_0_zz[i] - rab_z[i] * t_0_xx_0_zz[i];

        t_z_xx_0_yz[i] = t_0_xxz_0_yz[i] - rab_z[i] * t_0_xx_0_yz[i];

        t_z_xx_0_yy[i] = t_0_xxz_0_yy[i] - rab_z[i] * t_0_xx_0_yy[i];

        t_z_xx_0_xz[i] = t_0_xxz_0_xz[i] - rab_z[i] * t_0_xx_0_xz[i];

        t_z_xx_0_xy[i] = t_0_xxz_0_xy[i] - rab_z[i] * t_0_xx_0_xy[i];

        t_z_xx_0_xx[i] = t_0_xxz_0_xx[i] - rab_z[i] * t_0_xx_0_xx[i];
    }

    #pragma omp simd align(rab_y, t_0_xx_0_xx, t_0_xx_0_xy, t_0_xx_0_xz, t_0_xx_0_yy,\
                           t_0_xx_0_yz, t_0_xx_0_zz, t_0_xxy_0_xx, t_0_xxy_0_xy, t_0_xxy_0_xz,\
                           t_0_xxy_0_yy, t_0_xxy_0_yz, t_0_xxy_0_zz, t_0_xy_0_xx, t_0_xy_0_xy,\
                           t_0_xy_0_xz, t_0_xy_0_yy, t_0_xy_0_yz, t_0_xy_0_zz, t_0_xyy_0_xx,\
                           t_0_xyy_0_xy, t_0_xyy_0_xz, t_0_xyy_0_yy, t_0_xyy_0_yz, t_0_xyy_0_zz,\
                           t_0_xyz_0_xx, t_0_xyz_0_xy, t_0_xyz_0_xz, t_0_xyz_0_yy, t_0_xyz_0_yz,\
                           t_0_xyz_0_zz, t_0_xz_0_xx, t_0_xz_0_xy, t_0_xz_0_xz, t_0_xz_0_yy,\
                           t_0_xz_0_yz, t_0_xz_0_zz, t_0_yy_0_xx, t_0_yy_0_xy, t_0_yy_0_xz,\
                           t_0_yy_0_yy, t_0_yy_0_yz, t_0_yy_0_zz, t_0_yyy_0_xx, t_0_yyy_0_xy,\
                           t_0_yyy_0_xz, t_0_yyy_0_yy, t_0_yyy_0_yz, t_0_yyy_0_zz, t_0_yyz_0_xx,\
                           t_0_yyz_0_xy, t_0_yyz_0_xz, t_0_yyz_0_yy, t_0_yyz_0_yz, t_0_yyz_0_zz,\
                           t_0_yz_0_xx, t_0_yz_0_xy, t_0_yz_0_xz, t_0_yz_0_yy, t_0_yz_0_yz,\
                           t_0_yz_0_zz, t_0_yzz_0_xx, t_0_yzz_0_xy, t_0_yzz_0_xz, t_0_yzz_0_yy,\
                           t_0_yzz_0_yz, t_0_yzz_0_zz, t_0_zz_0_xx, t_0_zz_0_xy, t_0_zz_0_xz,\
                           t_0_zz_0_yy, t_0_zz_0_yz, t_0_zz_0_zz, t_y_xx_0_xx, t_y_xx_0_xy,\
                           t_y_xx_0_xz, t_y_xx_0_yy, t_y_xx_0_yz, t_y_xx_0_zz, t_y_xy_0_xx,\
                           t_y_xy_0_xy, t_y_xy_0_xz, t_y_xy_0_yy, t_y_xy_0_yz, t_y_xy_0_zz,\
                           t_y_xz_0_xx, t_y_xz_0_xy, t_y_xz_0_xz, t_y_xz_0_yy, t_y_xz_0_yz,\
                           t_y_xz_0_zz, t_y_yy_0_xx, t_y_yy_0_xy, t_y_yy_0_xz, t_y_yy_0_yy,\
                           t_y_yy_0_yz, t_y_yy_0_zz, t_y_yz_0_xx, t_y_yz_0_xy, t_y_yz_0_xz,\
                           t_y_yz_0_yy, t_y_yz_0_yz, t_y_yz_0_zz, t_y_zz_0_xx, t_y_zz_0_xy,\
                           t_y_zz_0_xz, t_y_zz_0_yy, t_y_zz_0_yz, t_y_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_zz_0_zz[i] = t_0_yzz_0_zz[i] - rab_y[i] * t_0_zz_0_zz[i];

        t_y_zz_0_yz[i] = t_0_yzz_0_yz[i] - rab_y[i] * t_0_zz_0_yz[i];

        t_y_zz_0_yy[i] = t_0_yzz_0_yy[i] - rab_y[i] * t_0_zz_0_yy[i];

        t_y_zz_0_xz[i] = t_0_yzz_0_xz[i] - rab_y[i] * t_0_zz_0_xz[i];

        t_y_zz_0_xy[i] = t_0_yzz_0_xy[i] - rab_y[i] * t_0_zz_0_xy[i];

        t_y_zz_0_xx[i] = t_0_yzz_0_xx[i] - rab_y[i] * t_0_zz_0_xx[i];

        t_y_yz_0_zz[i] = t_0_yyz_0_zz[i] - rab_y[i] * t_0_yz_0_zz[i];

        t_y_yz_0_yz[i] = t_0_yyz_0_yz[i] - rab_y[i] * t_0_yz_0_yz[i];

        t_y_yz_0_yy[i] = t_0_yyz_0_yy[i] - rab_y[i] * t_0_yz_0_yy[i];

        t_y_yz_0_xz[i] = t_0_yyz_0_xz[i] - rab_y[i] * t_0_yz_0_xz[i];

        t_y_yz_0_xy[i] = t_0_yyz_0_xy[i] - rab_y[i] * t_0_yz_0_xy[i];

        t_y_yz_0_xx[i] = t_0_yyz_0_xx[i] - rab_y[i] * t_0_yz_0_xx[i];

        t_y_yy_0_zz[i] = t_0_yyy_0_zz[i] - rab_y[i] * t_0_yy_0_zz[i];

        t_y_yy_0_yz[i] = t_0_yyy_0_yz[i] - rab_y[i] * t_0_yy_0_yz[i];

        t_y_yy_0_yy[i] = t_0_yyy_0_yy[i] - rab_y[i] * t_0_yy_0_yy[i];

        t_y_yy_0_xz[i] = t_0_yyy_0_xz[i] - rab_y[i] * t_0_yy_0_xz[i];

        t_y_yy_0_xy[i] = t_0_yyy_0_xy[i] - rab_y[i] * t_0_yy_0_xy[i];

        t_y_yy_0_xx[i] = t_0_yyy_0_xx[i] - rab_y[i] * t_0_yy_0_xx[i];

        t_y_xz_0_zz[i] = t_0_xyz_0_zz[i] - rab_y[i] * t_0_xz_0_zz[i];

        t_y_xz_0_yz[i] = t_0_xyz_0_yz[i] - rab_y[i] * t_0_xz_0_yz[i];

        t_y_xz_0_yy[i] = t_0_xyz_0_yy[i] - rab_y[i] * t_0_xz_0_yy[i];

        t_y_xz_0_xz[i] = t_0_xyz_0_xz[i] - rab_y[i] * t_0_xz_0_xz[i];

        t_y_xz_0_xy[i] = t_0_xyz_0_xy[i] - rab_y[i] * t_0_xz_0_xy[i];

        t_y_xz_0_xx[i] = t_0_xyz_0_xx[i] - rab_y[i] * t_0_xz_0_xx[i];

        t_y_xy_0_zz[i] = t_0_xyy_0_zz[i] - rab_y[i] * t_0_xy_0_zz[i];

        t_y_xy_0_yz[i] = t_0_xyy_0_yz[i] - rab_y[i] * t_0_xy_0_yz[i];

        t_y_xy_0_yy[i] = t_0_xyy_0_yy[i] - rab_y[i] * t_0_xy_0_yy[i];

        t_y_xy_0_xz[i] = t_0_xyy_0_xz[i] - rab_y[i] * t_0_xy_0_xz[i];

        t_y_xy_0_xy[i] = t_0_xyy_0_xy[i] - rab_y[i] * t_0_xy_0_xy[i];

        t_y_xy_0_xx[i] = t_0_xyy_0_xx[i] - rab_y[i] * t_0_xy_0_xx[i];

        t_y_xx_0_zz[i] = t_0_xxy_0_zz[i] - rab_y[i] * t_0_xx_0_zz[i];

        t_y_xx_0_yz[i] = t_0_xxy_0_yz[i] - rab_y[i] * t_0_xx_0_yz[i];

        t_y_xx_0_yy[i] = t_0_xxy_0_yy[i] - rab_y[i] * t_0_xx_0_yy[i];

        t_y_xx_0_xz[i] = t_0_xxy_0_xz[i] - rab_y[i] * t_0_xx_0_xz[i];

        t_y_xx_0_xy[i] = t_0_xxy_0_xy[i] - rab_y[i] * t_0_xx_0_xy[i];

        t_y_xx_0_xx[i] = t_0_xxy_0_xx[i] - rab_y[i] * t_0_xx_0_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_xx_0_xx, t_0_xx_0_xy, t_0_xx_0_xz, t_0_xx_0_yy,\
                           t_0_xx_0_yz, t_0_xx_0_zz, t_0_xxx_0_xx, t_0_xxx_0_xy, t_0_xxx_0_xz,\
                           t_0_xxx_0_yy, t_0_xxx_0_yz, t_0_xxx_0_zz, t_0_xxy_0_xx, t_0_xxy_0_xy,\
                           t_0_xxy_0_xz, t_0_xxy_0_yy, t_0_xxy_0_yz, t_0_xxy_0_zz, t_0_xxz_0_xx,\
                           t_0_xxz_0_xy, t_0_xxz_0_xz, t_0_xxz_0_yy, t_0_xxz_0_yz, t_0_xxz_0_zz,\
                           t_0_xy_0_xx, t_0_xy_0_xy, t_0_xy_0_xz, t_0_xy_0_yy, t_0_xy_0_yz,\
                           t_0_xy_0_zz, t_0_xyy_0_xx, t_0_xyy_0_xy, t_0_xyy_0_xz, t_0_xyy_0_yy,\
                           t_0_xyy_0_yz, t_0_xyy_0_zz, t_0_xyz_0_xx, t_0_xyz_0_xy, t_0_xyz_0_xz,\
                           t_0_xyz_0_yy, t_0_xyz_0_yz, t_0_xyz_0_zz, t_0_xz_0_xx, t_0_xz_0_xy,\
                           t_0_xz_0_xz, t_0_xz_0_yy, t_0_xz_0_yz, t_0_xz_0_zz, t_0_xzz_0_xx,\
                           t_0_xzz_0_xy, t_0_xzz_0_xz, t_0_xzz_0_yy, t_0_xzz_0_yz, t_0_xzz_0_zz,\
                           t_0_yy_0_xx, t_0_yy_0_xy, t_0_yy_0_xz, t_0_yy_0_yy, t_0_yy_0_yz,\
                           t_0_yy_0_zz, t_0_yz_0_xx, t_0_yz_0_xy, t_0_yz_0_xz, t_0_yz_0_yy,\
                           t_0_yz_0_yz, t_0_yz_0_zz, t_0_zz_0_xx, t_0_zz_0_xy, t_0_zz_0_xz,\
                           t_0_zz_0_yy, t_0_zz_0_yz, t_0_zz_0_zz, t_x_xx_0_xx, t_x_xx_0_xy,\
                           t_x_xx_0_xz, t_x_xx_0_yy, t_x_xx_0_yz, t_x_xx_0_zz, t_x_xy_0_xx,\
                           t_x_xy_0_xy, t_x_xy_0_xz, t_x_xy_0_yy, t_x_xy_0_yz, t_x_xy_0_zz,\
                           t_x_xz_0_xx, t_x_xz_0_xy, t_x_xz_0_xz, t_x_xz_0_yy, t_x_xz_0_yz,\
                           t_x_xz_0_zz, t_x_yy_0_xx, t_x_yy_0_xy, t_x_yy_0_xz, t_x_yy_0_yy,\
                           t_x_yy_0_yz, t_x_yy_0_zz, t_x_yz_0_xx, t_x_yz_0_xy, t_x_yz_0_xz,\
                           t_x_yz_0_yy, t_x_yz_0_yz, t_x_yz_0_zz, t_x_zz_0_xx, t_x_zz_0_xy,\
                           t_x_zz_0_xz, t_x_zz_0_yy, t_x_zz_0_yz, t_x_zz_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zz_0_zz[i] = t_0_xzz_0_zz[i] - rab_x[i] * t_0_zz_0_zz[i];

        t_x_zz_0_yz[i] = t_0_xzz_0_yz[i] - rab_x[i] * t_0_zz_0_yz[i];

        t_x_zz_0_yy[i] = t_0_xzz_0_yy[i] - rab_x[i] * t_0_zz_0_yy[i];

        t_x_zz_0_xz[i] = t_0_xzz_0_xz[i] - rab_x[i] * t_0_zz_0_xz[i];

        t_x_zz_0_xy[i] = t_0_xzz_0_xy[i] - rab_x[i] * t_0_zz_0_xy[i];

        t_x_zz_0_xx[i] = t_0_xzz_0_xx[i] - rab_x[i] * t_0_zz_0_xx[i];

        t_x_yz_0_zz[i] = t_0_xyz_0_zz[i] - rab_x[i] * t_0_yz_0_zz[i];

        t_x_yz_0_yz[i] = t_0_xyz_0_yz[i] - rab_x[i] * t_0_yz_0_yz[i];

        t_x_yz_0_yy[i] = t_0_xyz_0_yy[i] - rab_x[i] * t_0_yz_0_yy[i];

        t_x_yz_0_xz[i] = t_0_xyz_0_xz[i] - rab_x[i] * t_0_yz_0_xz[i];

        t_x_yz_0_xy[i] = t_0_xyz_0_xy[i] - rab_x[i] * t_0_yz_0_xy[i];

        t_x_yz_0_xx[i] = t_0_xyz_0_xx[i] - rab_x[i] * t_0_yz_0_xx[i];

        t_x_yy_0_zz[i] = t_0_xyy_0_zz[i] - rab_x[i] * t_0_yy_0_zz[i];

        t_x_yy_0_yz[i] = t_0_xyy_0_yz[i] - rab_x[i] * t_0_yy_0_yz[i];

        t_x_yy_0_yy[i] = t_0_xyy_0_yy[i] - rab_x[i] * t_0_yy_0_yy[i];

        t_x_yy_0_xz[i] = t_0_xyy_0_xz[i] - rab_x[i] * t_0_yy_0_xz[i];

        t_x_yy_0_xy[i] = t_0_xyy_0_xy[i] - rab_x[i] * t_0_yy_0_xy[i];

        t_x_yy_0_xx[i] = t_0_xyy_0_xx[i] - rab_x[i] * t_0_yy_0_xx[i];

        t_x_xz_0_zz[i] = t_0_xxz_0_zz[i] - rab_x[i] * t_0_xz_0_zz[i];

        t_x_xz_0_yz[i] = t_0_xxz_0_yz[i] - rab_x[i] * t_0_xz_0_yz[i];

        t_x_xz_0_yy[i] = t_0_xxz_0_yy[i] - rab_x[i] * t_0_xz_0_yy[i];

        t_x_xz_0_xz[i] = t_0_xxz_0_xz[i] - rab_x[i] * t_0_xz_0_xz[i];

        t_x_xz_0_xy[i] = t_0_xxz_0_xy[i] - rab_x[i] * t_0_xz_0_xy[i];

        t_x_xz_0_xx[i] = t_0_xxz_0_xx[i] - rab_x[i] * t_0_xz_0_xx[i];

        t_x_xy_0_zz[i] = t_0_xxy_0_zz[i] - rab_x[i] * t_0_xy_0_zz[i];

        t_x_xy_0_yz[i] = t_0_xxy_0_yz[i] - rab_x[i] * t_0_xy_0_yz[i];

        t_x_xy_0_yy[i] = t_0_xxy_0_yy[i] - rab_x[i] * t_0_xy_0_yy[i];

        t_x_xy_0_xz[i] = t_0_xxy_0_xz[i] - rab_x[i] * t_0_xy_0_xz[i];

        t_x_xy_0_xy[i] = t_0_xxy_0_xy[i] - rab_x[i] * t_0_xy_0_xy[i];

        t_x_xy_0_xx[i] = t_0_xxy_0_xx[i] - rab_x[i] * t_0_xy_0_xx[i];

        t_x_xx_0_zz[i] = t_0_xxx_0_zz[i] - rab_x[i] * t_0_xx_0_zz[i];

        t_x_xx_0_yz[i] = t_0_xxx_0_yz[i] - rab_x[i] * t_0_xx_0_yz[i];

        t_x_xx_0_yy[i] = t_0_xxx_0_yy[i] - rab_x[i] * t_0_xx_0_yy[i];

        t_x_xx_0_xz[i] = t_0_xxx_0_xz[i] - rab_x[i] * t_0_xx_0_xz[i];

        t_x_xx_0_xy[i] = t_0_xxx_0_xy[i] - rab_x[i] * t_0_xx_0_xy[i];

        t_x_xx_0_xx[i] = t_0_xxx_0_xx[i] - rab_x[i] * t_0_xx_0_xx[i];
    }
}


} // derirec namespace
