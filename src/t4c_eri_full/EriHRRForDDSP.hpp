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
compHostHRRForDDSP_V0(      BufferHostXY<T>&      intsBufferDDSP,
                      const BufferHostX<int32_t>& intsIndexesDDSP,
                      const BufferHostXY<T>&      intsBufferPDSP,
                      const BufferHostX<int32_t>& intsIndexesPDSP,
                      const BufferHostXY<T>&      intsBufferPFSP,
                      const BufferHostX<int32_t>& intsIndexesPFSP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDSP) integral components

    t_zz_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(0));

    t_zz_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(1));

    t_zz_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(2));

    t_zz_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(3));

    t_zz_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(4));

    t_zz_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(5));

    t_zz_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(6));

    t_zz_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(7));

    t_zz_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(8));

    t_zz_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(9));

    t_zz_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(10));

    t_zz_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(11));

    t_zz_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(12));

    t_zz_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(13));

    t_zz_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(14));

    t_zz_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(15));

    t_zz_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(16));

    t_zz_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(17));

    t_yz_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(18));

    t_yz_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(19));

    t_yz_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(20));

    t_yz_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(21));

    t_yz_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(22));

    t_yz_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(23));

    t_yz_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(24));

    t_yz_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(25));

    t_yz_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(26));

    t_yz_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(27));

    t_yz_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(28));

    t_yz_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(29));

    t_yz_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(30));

    t_yz_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(31));

    t_yz_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(32));

    t_yz_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(33));

    t_yz_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(34));

    t_yz_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(35));

    t_yy_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(36));

    t_yy_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(37));

    t_yy_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(38));

    t_yy_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(39));

    t_yy_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(40));

    t_yy_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(41));

    t_yy_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(42));

    t_yy_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(43));

    t_yy_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(44));

    t_yy_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(45));

    t_yy_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(46));

    t_yy_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(47));

    t_yy_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(48));

    t_yy_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(49));

    t_yy_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(50));

    t_yy_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(51));

    t_yy_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(52));

    t_yy_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(53));

    t_xz_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(54));

    t_xz_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(55));

    t_xz_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(56));

    t_xz_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(57));

    t_xz_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(58));

    t_xz_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(59));

    t_xz_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(60));

    t_xz_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(61));

    t_xz_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(62));

    t_xz_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(63));

    t_xz_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(64));

    t_xz_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(65));

    t_xz_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(66));

    t_xz_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(67));

    t_xz_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(68));

    t_xz_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(69));

    t_xz_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(70));

    t_xz_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(71));

    t_xy_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(72));

    t_xy_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(73));

    t_xy_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(74));

    t_xy_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(75));

    t_xy_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(76));

    t_xy_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(77));

    t_xy_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(78));

    t_xy_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(79));

    t_xy_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(80));

    t_xy_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(81));

    t_xy_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(82));

    t_xy_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(83));

    t_xy_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(84));

    t_xy_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(85));

    t_xy_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(86));

    t_xy_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(87));

    t_xy_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(88));

    t_xy_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(89));

    t_xx_zz_0_z = intsBufferDDSP.data(intsIndexesDDSP(90));

    t_xx_zz_0_y = intsBufferDDSP.data(intsIndexesDDSP(91));

    t_xx_zz_0_x = intsBufferDDSP.data(intsIndexesDDSP(92));

    t_xx_yz_0_z = intsBufferDDSP.data(intsIndexesDDSP(93));

    t_xx_yz_0_y = intsBufferDDSP.data(intsIndexesDDSP(94));

    t_xx_yz_0_x = intsBufferDDSP.data(intsIndexesDDSP(95));

    t_xx_yy_0_z = intsBufferDDSP.data(intsIndexesDDSP(96));

    t_xx_yy_0_y = intsBufferDDSP.data(intsIndexesDDSP(97));

    t_xx_yy_0_x = intsBufferDDSP.data(intsIndexesDDSP(98));

    t_xx_xz_0_z = intsBufferDDSP.data(intsIndexesDDSP(99));

    t_xx_xz_0_y = intsBufferDDSP.data(intsIndexesDDSP(100));

    t_xx_xz_0_x = intsBufferDDSP.data(intsIndexesDDSP(101));

    t_xx_xy_0_z = intsBufferDDSP.data(intsIndexesDDSP(102));

    t_xx_xy_0_y = intsBufferDDSP.data(intsIndexesDDSP(103));

    t_xx_xy_0_x = intsBufferDDSP.data(intsIndexesDDSP(104));

    t_xx_xx_0_z = intsBufferDDSP.data(intsIndexesDDSP(105));

    t_xx_xx_0_y = intsBufferDDSP.data(intsIndexesDDSP(106));

    t_xx_xx_0_x = intsBufferDDSP.data(intsIndexesDDSP(107));

    // set up (PDSP) integral components

    t_z_zz_0_z = intsBufferPDSP.data(intsIndexesPDSP(0));

    t_z_zz_0_y = intsBufferPDSP.data(intsIndexesPDSP(1));

    t_z_zz_0_x = intsBufferPDSP.data(intsIndexesPDSP(2));

    t_z_yz_0_z = intsBufferPDSP.data(intsIndexesPDSP(3));

    t_z_yz_0_y = intsBufferPDSP.data(intsIndexesPDSP(4));

    t_z_yz_0_x = intsBufferPDSP.data(intsIndexesPDSP(5));

    t_z_yy_0_z = intsBufferPDSP.data(intsIndexesPDSP(6));

    t_z_yy_0_y = intsBufferPDSP.data(intsIndexesPDSP(7));

    t_z_yy_0_x = intsBufferPDSP.data(intsIndexesPDSP(8));

    t_z_xz_0_z = intsBufferPDSP.data(intsIndexesPDSP(9));

    t_z_xz_0_y = intsBufferPDSP.data(intsIndexesPDSP(10));

    t_z_xz_0_x = intsBufferPDSP.data(intsIndexesPDSP(11));

    t_z_xy_0_z = intsBufferPDSP.data(intsIndexesPDSP(12));

    t_z_xy_0_y = intsBufferPDSP.data(intsIndexesPDSP(13));

    t_z_xy_0_x = intsBufferPDSP.data(intsIndexesPDSP(14));

    t_z_xx_0_z = intsBufferPDSP.data(intsIndexesPDSP(15));

    t_z_xx_0_y = intsBufferPDSP.data(intsIndexesPDSP(16));

    t_z_xx_0_x = intsBufferPDSP.data(intsIndexesPDSP(17));

    t_y_zz_0_z = intsBufferPDSP.data(intsIndexesPDSP(18));

    t_y_zz_0_y = intsBufferPDSP.data(intsIndexesPDSP(19));

    t_y_zz_0_x = intsBufferPDSP.data(intsIndexesPDSP(20));

    t_y_yz_0_z = intsBufferPDSP.data(intsIndexesPDSP(21));

    t_y_yz_0_y = intsBufferPDSP.data(intsIndexesPDSP(22));

    t_y_yz_0_x = intsBufferPDSP.data(intsIndexesPDSP(23));

    t_y_yy_0_z = intsBufferPDSP.data(intsIndexesPDSP(24));

    t_y_yy_0_y = intsBufferPDSP.data(intsIndexesPDSP(25));

    t_y_yy_0_x = intsBufferPDSP.data(intsIndexesPDSP(26));

    t_y_xz_0_z = intsBufferPDSP.data(intsIndexesPDSP(27));

    t_y_xz_0_y = intsBufferPDSP.data(intsIndexesPDSP(28));

    t_y_xz_0_x = intsBufferPDSP.data(intsIndexesPDSP(29));

    t_y_xy_0_z = intsBufferPDSP.data(intsIndexesPDSP(30));

    t_y_xy_0_y = intsBufferPDSP.data(intsIndexesPDSP(31));

    t_y_xy_0_x = intsBufferPDSP.data(intsIndexesPDSP(32));

    t_y_xx_0_z = intsBufferPDSP.data(intsIndexesPDSP(33));

    t_y_xx_0_y = intsBufferPDSP.data(intsIndexesPDSP(34));

    t_y_xx_0_x = intsBufferPDSP.data(intsIndexesPDSP(35));

    t_x_zz_0_z = intsBufferPDSP.data(intsIndexesPDSP(36));

    t_x_zz_0_y = intsBufferPDSP.data(intsIndexesPDSP(37));

    t_x_zz_0_x = intsBufferPDSP.data(intsIndexesPDSP(38));

    t_x_yz_0_z = intsBufferPDSP.data(intsIndexesPDSP(39));

    t_x_yz_0_y = intsBufferPDSP.data(intsIndexesPDSP(40));

    t_x_yz_0_x = intsBufferPDSP.data(intsIndexesPDSP(41));

    t_x_yy_0_z = intsBufferPDSP.data(intsIndexesPDSP(42));

    t_x_yy_0_y = intsBufferPDSP.data(intsIndexesPDSP(43));

    t_x_yy_0_x = intsBufferPDSP.data(intsIndexesPDSP(44));

    t_x_xz_0_z = intsBufferPDSP.data(intsIndexesPDSP(45));

    t_x_xz_0_y = intsBufferPDSP.data(intsIndexesPDSP(46));

    t_x_xz_0_x = intsBufferPDSP.data(intsIndexesPDSP(47));

    t_x_xy_0_z = intsBufferPDSP.data(intsIndexesPDSP(48));

    t_x_xy_0_y = intsBufferPDSP.data(intsIndexesPDSP(49));

    t_x_xy_0_x = intsBufferPDSP.data(intsIndexesPDSP(50));

    t_x_xx_0_z = intsBufferPDSP.data(intsIndexesPDSP(51));

    t_x_xx_0_y = intsBufferPDSP.data(intsIndexesPDSP(52));

    t_x_xx_0_x = intsBufferPDSP.data(intsIndexesPDSP(53));

    // set up (PFSP) integral components

    t_z_zzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(0));

    t_z_zzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(1));

    t_z_zzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(2));

    t_z_yzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(3));

    t_z_yzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(4));

    t_z_yzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(5));

    t_z_yyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(6));

    t_z_yyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(7));

    t_z_yyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(8));

    t_z_xzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(9));

    t_z_xzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(10));

    t_z_xzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(11));

    t_z_xyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(12));

    t_z_xyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(13));

    t_z_xyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(14));

    t_z_xxz_0_z = intsBufferPFSP.data(intsIndexesPFSP(15));

    t_z_xxz_0_y = intsBufferPFSP.data(intsIndexesPFSP(16));

    t_z_xxz_0_x = intsBufferPFSP.data(intsIndexesPFSP(17));

    t_y_zzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(18));

    t_y_zzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(19));

    t_y_zzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(20));

    t_y_yzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(21));

    t_y_yzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(22));

    t_y_yzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(23));

    t_y_yyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(24));

    t_y_yyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(25));

    t_y_yyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(26));

    t_y_yyy_0_z = intsBufferPFSP.data(intsIndexesPFSP(27));

    t_y_yyy_0_y = intsBufferPFSP.data(intsIndexesPFSP(28));

    t_y_yyy_0_x = intsBufferPFSP.data(intsIndexesPFSP(29));

    t_y_xzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(30));

    t_y_xzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(31));

    t_y_xzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(32));

    t_y_xyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(33));

    t_y_xyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(34));

    t_y_xyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(35));

    t_y_xyy_0_z = intsBufferPFSP.data(intsIndexesPFSP(36));

    t_y_xyy_0_y = intsBufferPFSP.data(intsIndexesPFSP(37));

    t_y_xyy_0_x = intsBufferPFSP.data(intsIndexesPFSP(38));

    t_y_xxz_0_z = intsBufferPFSP.data(intsIndexesPFSP(39));

    t_y_xxz_0_y = intsBufferPFSP.data(intsIndexesPFSP(40));

    t_y_xxz_0_x = intsBufferPFSP.data(intsIndexesPFSP(41));

    t_y_xxy_0_z = intsBufferPFSP.data(intsIndexesPFSP(42));

    t_y_xxy_0_y = intsBufferPFSP.data(intsIndexesPFSP(43));

    t_y_xxy_0_x = intsBufferPFSP.data(intsIndexesPFSP(44));

    t_x_zzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(45));

    t_x_zzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(46));

    t_x_zzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(47));

    t_x_yzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(48));

    t_x_yzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(49));

    t_x_yzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(50));

    t_x_yyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(51));

    t_x_yyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(52));

    t_x_yyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(53));

    t_x_yyy_0_z = intsBufferPFSP.data(intsIndexesPFSP(54));

    t_x_yyy_0_y = intsBufferPFSP.data(intsIndexesPFSP(55));

    t_x_yyy_0_x = intsBufferPFSP.data(intsIndexesPFSP(56));

    t_x_xzz_0_z = intsBufferPFSP.data(intsIndexesPFSP(57));

    t_x_xzz_0_y = intsBufferPFSP.data(intsIndexesPFSP(58));

    t_x_xzz_0_x = intsBufferPFSP.data(intsIndexesPFSP(59));

    t_x_xyz_0_z = intsBufferPFSP.data(intsIndexesPFSP(60));

    t_x_xyz_0_y = intsBufferPFSP.data(intsIndexesPFSP(61));

    t_x_xyz_0_x = intsBufferPFSP.data(intsIndexesPFSP(62));

    t_x_xyy_0_z = intsBufferPFSP.data(intsIndexesPFSP(63));

    t_x_xyy_0_y = intsBufferPFSP.data(intsIndexesPFSP(64));

    t_x_xyy_0_x = intsBufferPFSP.data(intsIndexesPFSP(65));

    t_x_xxz_0_z = intsBufferPFSP.data(intsIndexesPFSP(66));

    t_x_xxz_0_y = intsBufferPFSP.data(intsIndexesPFSP(67));

    t_x_xxz_0_x = intsBufferPFSP.data(intsIndexesPFSP(68));

    t_x_xxy_0_z = intsBufferPFSP.data(intsIndexesPFSP(69));

    t_x_xxy_0_y = intsBufferPFSP.data(intsIndexesPFSP(70));

    t_x_xxy_0_x = intsBufferPFSP.data(intsIndexesPFSP(71));

    t_x_xxx_0_z = intsBufferPFSP.data(intsIndexesPFSP(72));

    t_x_xxx_0_y = intsBufferPFSP.data(intsIndexesPFSP(73));

    t_x_xxx_0_x = intsBufferPFSP.data(intsIndexesPFSP(74));

    #pragma omp simd align(rab_z, t_y_xx_0_x, t_y_xx_0_y, t_y_xx_0_z, t_y_xxz_0_x, t_y_xxz_0_y,\
                           t_y_xxz_0_z, t_y_xy_0_x, t_y_xy_0_y, t_y_xy_0_z, t_y_xyz_0_x,\
                           t_y_xyz_0_y, t_y_xyz_0_z, t_y_xz_0_x, t_y_xz_0_y, t_y_xz_0_z,\
                           t_y_xzz_0_x, t_y_xzz_0_y, t_y_xzz_0_z, t_y_yy_0_x, t_y_yy_0_y,\
                           t_y_yy_0_z, t_y_yyz_0_x, t_y_yyz_0_y, t_y_yyz_0_z, t_y_yz_0_x,\
                           t_y_yz_0_y, t_y_yz_0_z, t_y_yzz_0_x, t_y_yzz_0_y, t_y_yzz_0_z,\
                           t_y_zz_0_x, t_y_zz_0_y, t_y_zz_0_z, t_y_zzz_0_x, t_y_zzz_0_y,\
                           t_y_zzz_0_z, t_yz_xx_0_x, t_yz_xx_0_y, t_yz_xx_0_z, t_yz_xy_0_x,\
                           t_yz_xy_0_y, t_yz_xy_0_z, t_yz_xz_0_x, t_yz_xz_0_y, t_yz_xz_0_z,\
                           t_yz_yy_0_x, t_yz_yy_0_y, t_yz_yy_0_z, t_yz_yz_0_x, t_yz_yz_0_y,\
                           t_yz_yz_0_z, t_yz_zz_0_x, t_yz_zz_0_y, t_yz_zz_0_z, t_z_xx_0_x,\
                           t_z_xx_0_y, t_z_xx_0_z, t_z_xxz_0_x, t_z_xxz_0_y, t_z_xxz_0_z,\
                           t_z_xy_0_x, t_z_xy_0_y, t_z_xy_0_z, t_z_xyz_0_x, t_z_xyz_0_y,\
                           t_z_xyz_0_z, t_z_xz_0_x, t_z_xz_0_y, t_z_xz_0_z, t_z_xzz_0_x,\
                           t_z_xzz_0_y, t_z_xzz_0_z, t_z_yy_0_x, t_z_yy_0_y, t_z_yy_0_z,\
                           t_z_yyz_0_x, t_z_yyz_0_y, t_z_yyz_0_z, t_z_yz_0_x, t_z_yz_0_y,\
                           t_z_yz_0_z, t_z_yzz_0_x, t_z_yzz_0_y, t_z_yzz_0_z, t_z_zz_0_x,\
                           t_z_zz_0_y, t_z_zz_0_z, t_z_zzz_0_x, t_z_zzz_0_y, t_z_zzz_0_z,\
                           t_zz_xx_0_x, t_zz_xx_0_y, t_zz_xx_0_z, t_zz_xy_0_x, t_zz_xy_0_y,\
                           t_zz_xy_0_z, t_zz_xz_0_x, t_zz_xz_0_y, t_zz_xz_0_z, t_zz_yy_0_x,\
                           t_zz_yy_0_y, t_zz_yy_0_z, t_zz_yz_0_x, t_zz_yz_0_y, t_zz_yz_0_z,\
                           t_zz_zz_0_x, t_zz_zz_0_y, t_zz_zz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_0_z[i] = t_z_zzz_0_z[i] - rab_z[i] * t_z_zz_0_z[i];

        t_zz_zz_0_y[i] = t_z_zzz_0_y[i] - rab_z[i] * t_z_zz_0_y[i];

        t_zz_zz_0_x[i] = t_z_zzz_0_x[i] - rab_z[i] * t_z_zz_0_x[i];

        t_zz_yz_0_z[i] = t_z_yzz_0_z[i] - rab_z[i] * t_z_yz_0_z[i];

        t_zz_yz_0_y[i] = t_z_yzz_0_y[i] - rab_z[i] * t_z_yz_0_y[i];

        t_zz_yz_0_x[i] = t_z_yzz_0_x[i] - rab_z[i] * t_z_yz_0_x[i];

        t_zz_yy_0_z[i] = t_z_yyz_0_z[i] - rab_z[i] * t_z_yy_0_z[i];

        t_zz_yy_0_y[i] = t_z_yyz_0_y[i] - rab_z[i] * t_z_yy_0_y[i];

        t_zz_yy_0_x[i] = t_z_yyz_0_x[i] - rab_z[i] * t_z_yy_0_x[i];

        t_zz_xz_0_z[i] = t_z_xzz_0_z[i] - rab_z[i] * t_z_xz_0_z[i];

        t_zz_xz_0_y[i] = t_z_xzz_0_y[i] - rab_z[i] * t_z_xz_0_y[i];

        t_zz_xz_0_x[i] = t_z_xzz_0_x[i] - rab_z[i] * t_z_xz_0_x[i];

        t_zz_xy_0_z[i] = t_z_xyz_0_z[i] - rab_z[i] * t_z_xy_0_z[i];

        t_zz_xy_0_y[i] = t_z_xyz_0_y[i] - rab_z[i] * t_z_xy_0_y[i];

        t_zz_xy_0_x[i] = t_z_xyz_0_x[i] - rab_z[i] * t_z_xy_0_x[i];

        t_zz_xx_0_z[i] = t_z_xxz_0_z[i] - rab_z[i] * t_z_xx_0_z[i];

        t_zz_xx_0_y[i] = t_z_xxz_0_y[i] - rab_z[i] * t_z_xx_0_y[i];

        t_zz_xx_0_x[i] = t_z_xxz_0_x[i] - rab_z[i] * t_z_xx_0_x[i];

        t_yz_zz_0_z[i] = t_y_zzz_0_z[i] - rab_z[i] * t_y_zz_0_z[i];

        t_yz_zz_0_y[i] = t_y_zzz_0_y[i] - rab_z[i] * t_y_zz_0_y[i];

        t_yz_zz_0_x[i] = t_y_zzz_0_x[i] - rab_z[i] * t_y_zz_0_x[i];

        t_yz_yz_0_z[i] = t_y_yzz_0_z[i] - rab_z[i] * t_y_yz_0_z[i];

        t_yz_yz_0_y[i] = t_y_yzz_0_y[i] - rab_z[i] * t_y_yz_0_y[i];

        t_yz_yz_0_x[i] = t_y_yzz_0_x[i] - rab_z[i] * t_y_yz_0_x[i];

        t_yz_yy_0_z[i] = t_y_yyz_0_z[i] - rab_z[i] * t_y_yy_0_z[i];

        t_yz_yy_0_y[i] = t_y_yyz_0_y[i] - rab_z[i] * t_y_yy_0_y[i];

        t_yz_yy_0_x[i] = t_y_yyz_0_x[i] - rab_z[i] * t_y_yy_0_x[i];

        t_yz_xz_0_z[i] = t_y_xzz_0_z[i] - rab_z[i] * t_y_xz_0_z[i];

        t_yz_xz_0_y[i] = t_y_xzz_0_y[i] - rab_z[i] * t_y_xz_0_y[i];

        t_yz_xz_0_x[i] = t_y_xzz_0_x[i] - rab_z[i] * t_y_xz_0_x[i];

        t_yz_xy_0_z[i] = t_y_xyz_0_z[i] - rab_z[i] * t_y_xy_0_z[i];

        t_yz_xy_0_y[i] = t_y_xyz_0_y[i] - rab_z[i] * t_y_xy_0_y[i];

        t_yz_xy_0_x[i] = t_y_xyz_0_x[i] - rab_z[i] * t_y_xy_0_x[i];

        t_yz_xx_0_z[i] = t_y_xxz_0_z[i] - rab_z[i] * t_y_xx_0_z[i];

        t_yz_xx_0_y[i] = t_y_xxz_0_y[i] - rab_z[i] * t_y_xx_0_y[i];

        t_yz_xx_0_x[i] = t_y_xxz_0_x[i] - rab_z[i] * t_y_xx_0_x[i];
    }

    #pragma omp simd align(rab_y, rab_z, t_x_xx_0_x, t_x_xx_0_y, t_x_xx_0_z, t_x_xxz_0_x,\
                           t_x_xxz_0_y, t_x_xxz_0_z, t_x_xy_0_x, t_x_xy_0_y, t_x_xy_0_z,\
                           t_x_xyz_0_x, t_x_xyz_0_y, t_x_xyz_0_z, t_x_xz_0_x, t_x_xz_0_y,\
                           t_x_xz_0_z, t_x_xzz_0_x, t_x_xzz_0_y, t_x_xzz_0_z, t_x_yy_0_x,\
                           t_x_yy_0_y, t_x_yy_0_z, t_x_yyz_0_x, t_x_yyz_0_y, t_x_yyz_0_z,\
                           t_x_yz_0_x, t_x_yz_0_y, t_x_yz_0_z, t_x_yzz_0_x, t_x_yzz_0_y,\
                           t_x_yzz_0_z, t_x_zz_0_x, t_x_zz_0_y, t_x_zz_0_z, t_x_zzz_0_x,\
                           t_x_zzz_0_y, t_x_zzz_0_z, t_xz_xx_0_x, t_xz_xx_0_y, t_xz_xx_0_z,\
                           t_xz_xy_0_x, t_xz_xy_0_y, t_xz_xy_0_z, t_xz_xz_0_x, t_xz_xz_0_y,\
                           t_xz_xz_0_z, t_xz_yy_0_x, t_xz_yy_0_y, t_xz_yy_0_z, t_xz_yz_0_x,\
                           t_xz_yz_0_y, t_xz_yz_0_z, t_xz_zz_0_x, t_xz_zz_0_y, t_xz_zz_0_z,\
                           t_y_xx_0_x, t_y_xx_0_y, t_y_xx_0_z, t_y_xxy_0_x, t_y_xxy_0_y,\
                           t_y_xxy_0_z, t_y_xy_0_x, t_y_xy_0_y, t_y_xy_0_z, t_y_xyy_0_x,\
                           t_y_xyy_0_y, t_y_xyy_0_z, t_y_xyz_0_x, t_y_xyz_0_y, t_y_xyz_0_z,\
                           t_y_xz_0_x, t_y_xz_0_y, t_y_xz_0_z, t_y_yy_0_x, t_y_yy_0_y,\
                           t_y_yy_0_z, t_y_yyy_0_x, t_y_yyy_0_y, t_y_yyy_0_z, t_y_yyz_0_x,\
                           t_y_yyz_0_y, t_y_yyz_0_z, t_y_yz_0_x, t_y_yz_0_y, t_y_yz_0_z,\
                           t_y_yzz_0_x, t_y_yzz_0_y, t_y_yzz_0_z, t_y_zz_0_x, t_y_zz_0_y,\
                           t_y_zz_0_z, t_yy_xx_0_x, t_yy_xx_0_y, t_yy_xx_0_z, t_yy_xy_0_x,\
                           t_yy_xy_0_y, t_yy_xy_0_z, t_yy_xz_0_x, t_yy_xz_0_y, t_yy_xz_0_z,\
                           t_yy_yy_0_x, t_yy_yy_0_y, t_yy_yy_0_z, t_yy_yz_0_x, t_yy_yz_0_y,\
                           t_yy_yz_0_z, t_yy_zz_0_x, t_yy_zz_0_y, t_yy_zz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_yy_zz_0_z[i] = t_y_yzz_0_z[i] - rab_y[i] * t_y_zz_0_z[i];

        t_yy_zz_0_y[i] = t_y_yzz_0_y[i] - rab_y[i] * t_y_zz_0_y[i];

        t_yy_zz_0_x[i] = t_y_yzz_0_x[i] - rab_y[i] * t_y_zz_0_x[i];

        t_yy_yz_0_z[i] = t_y_yyz_0_z[i] - rab_y[i] * t_y_yz_0_z[i];

        t_yy_yz_0_y[i] = t_y_yyz_0_y[i] - rab_y[i] * t_y_yz_0_y[i];

        t_yy_yz_0_x[i] = t_y_yyz_0_x[i] - rab_y[i] * t_y_yz_0_x[i];

        t_yy_yy_0_z[i] = t_y_yyy_0_z[i] - rab_y[i] * t_y_yy_0_z[i];

        t_yy_yy_0_y[i] = t_y_yyy_0_y[i] - rab_y[i] * t_y_yy_0_y[i];

        t_yy_yy_0_x[i] = t_y_yyy_0_x[i] - rab_y[i] * t_y_yy_0_x[i];

        t_yy_xz_0_z[i] = t_y_xyz_0_z[i] - rab_y[i] * t_y_xz_0_z[i];

        t_yy_xz_0_y[i] = t_y_xyz_0_y[i] - rab_y[i] * t_y_xz_0_y[i];

        t_yy_xz_0_x[i] = t_y_xyz_0_x[i] - rab_y[i] * t_y_xz_0_x[i];

        t_yy_xy_0_z[i] = t_y_xyy_0_z[i] - rab_y[i] * t_y_xy_0_z[i];

        t_yy_xy_0_y[i] = t_y_xyy_0_y[i] - rab_y[i] * t_y_xy_0_y[i];

        t_yy_xy_0_x[i] = t_y_xyy_0_x[i] - rab_y[i] * t_y_xy_0_x[i];

        t_yy_xx_0_z[i] = t_y_xxy_0_z[i] - rab_y[i] * t_y_xx_0_z[i];

        t_yy_xx_0_y[i] = t_y_xxy_0_y[i] - rab_y[i] * t_y_xx_0_y[i];

        t_yy_xx_0_x[i] = t_y_xxy_0_x[i] - rab_y[i] * t_y_xx_0_x[i];

        t_xz_zz_0_z[i] = t_x_zzz_0_z[i] - rab_z[i] * t_x_zz_0_z[i];

        t_xz_zz_0_y[i] = t_x_zzz_0_y[i] - rab_z[i] * t_x_zz_0_y[i];

        t_xz_zz_0_x[i] = t_x_zzz_0_x[i] - rab_z[i] * t_x_zz_0_x[i];

        t_xz_yz_0_z[i] = t_x_yzz_0_z[i] - rab_z[i] * t_x_yz_0_z[i];

        t_xz_yz_0_y[i] = t_x_yzz_0_y[i] - rab_z[i] * t_x_yz_0_y[i];

        t_xz_yz_0_x[i] = t_x_yzz_0_x[i] - rab_z[i] * t_x_yz_0_x[i];

        t_xz_yy_0_z[i] = t_x_yyz_0_z[i] - rab_z[i] * t_x_yy_0_z[i];

        t_xz_yy_0_y[i] = t_x_yyz_0_y[i] - rab_z[i] * t_x_yy_0_y[i];

        t_xz_yy_0_x[i] = t_x_yyz_0_x[i] - rab_z[i] * t_x_yy_0_x[i];

        t_xz_xz_0_z[i] = t_x_xzz_0_z[i] - rab_z[i] * t_x_xz_0_z[i];

        t_xz_xz_0_y[i] = t_x_xzz_0_y[i] - rab_z[i] * t_x_xz_0_y[i];

        t_xz_xz_0_x[i] = t_x_xzz_0_x[i] - rab_z[i] * t_x_xz_0_x[i];

        t_xz_xy_0_z[i] = t_x_xyz_0_z[i] - rab_z[i] * t_x_xy_0_z[i];

        t_xz_xy_0_y[i] = t_x_xyz_0_y[i] - rab_z[i] * t_x_xy_0_y[i];

        t_xz_xy_0_x[i] = t_x_xyz_0_x[i] - rab_z[i] * t_x_xy_0_x[i];

        t_xz_xx_0_z[i] = t_x_xxz_0_z[i] - rab_z[i] * t_x_xx_0_z[i];

        t_xz_xx_0_y[i] = t_x_xxz_0_y[i] - rab_z[i] * t_x_xx_0_y[i];

        t_xz_xx_0_x[i] = t_x_xxz_0_x[i] - rab_z[i] * t_x_xx_0_x[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_x_xx_0_x, t_x_xx_0_y, t_x_xx_0_z, t_x_xxx_0_x,\
                           t_x_xxx_0_y, t_x_xxx_0_z, t_x_xxy_0_x, t_x_xxy_0_y, t_x_xxy_0_z,\
                           t_x_xxz_0_x, t_x_xxz_0_y, t_x_xxz_0_z, t_x_xy_0_x, t_x_xy_0_y,\
                           t_x_xy_0_z, t_x_xyy_0_x, t_x_xyy_0_y, t_x_xyy_0_z, t_x_xyz_0_x,\
                           t_x_xyz_0_y, t_x_xyz_0_z, t_x_xz_0_x, t_x_xz_0_y, t_x_xz_0_z,\
                           t_x_xzz_0_x, t_x_xzz_0_y, t_x_xzz_0_z, t_x_yy_0_x, t_x_yy_0_y,\
                           t_x_yy_0_z, t_x_yyy_0_x, t_x_yyy_0_y, t_x_yyy_0_z, t_x_yyz_0_x,\
                           t_x_yyz_0_y, t_x_yyz_0_z, t_x_yz_0_x, t_x_yz_0_y, t_x_yz_0_z,\
                           t_x_yzz_0_x, t_x_yzz_0_y, t_x_yzz_0_z, t_x_zz_0_x, t_x_zz_0_y,\
                           t_x_zz_0_z, t_xx_xx_0_x, t_xx_xx_0_y, t_xx_xx_0_z, t_xx_xy_0_x,\
                           t_xx_xy_0_y, t_xx_xy_0_z, t_xx_xz_0_x, t_xx_xz_0_y, t_xx_xz_0_z,\
                           t_xx_yy_0_x, t_xx_yy_0_y, t_xx_yy_0_z, t_xx_yz_0_x, t_xx_yz_0_y,\
                           t_xx_yz_0_z, t_xx_zz_0_x, t_xx_zz_0_y, t_xx_zz_0_z, t_xy_xx_0_x,\
                           t_xy_xx_0_y, t_xy_xx_0_z, t_xy_xy_0_x, t_xy_xy_0_y, t_xy_xy_0_z,\
                           t_xy_xz_0_x, t_xy_xz_0_y, t_xy_xz_0_z, t_xy_yy_0_x, t_xy_yy_0_y,\
                           t_xy_yy_0_z, t_xy_yz_0_x, t_xy_yz_0_y, t_xy_yz_0_z, t_xy_zz_0_x,\
                           t_xy_zz_0_y, t_xy_zz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_xy_zz_0_z[i] = t_x_yzz_0_z[i] - rab_y[i] * t_x_zz_0_z[i];

        t_xy_zz_0_y[i] = t_x_yzz_0_y[i] - rab_y[i] * t_x_zz_0_y[i];

        t_xy_zz_0_x[i] = t_x_yzz_0_x[i] - rab_y[i] * t_x_zz_0_x[i];

        t_xy_yz_0_z[i] = t_x_yyz_0_z[i] - rab_y[i] * t_x_yz_0_z[i];

        t_xy_yz_0_y[i] = t_x_yyz_0_y[i] - rab_y[i] * t_x_yz_0_y[i];

        t_xy_yz_0_x[i] = t_x_yyz_0_x[i] - rab_y[i] * t_x_yz_0_x[i];

        t_xy_yy_0_z[i] = t_x_yyy_0_z[i] - rab_y[i] * t_x_yy_0_z[i];

        t_xy_yy_0_y[i] = t_x_yyy_0_y[i] - rab_y[i] * t_x_yy_0_y[i];

        t_xy_yy_0_x[i] = t_x_yyy_0_x[i] - rab_y[i] * t_x_yy_0_x[i];

        t_xy_xz_0_z[i] = t_x_xyz_0_z[i] - rab_y[i] * t_x_xz_0_z[i];

        t_xy_xz_0_y[i] = t_x_xyz_0_y[i] - rab_y[i] * t_x_xz_0_y[i];

        t_xy_xz_0_x[i] = t_x_xyz_0_x[i] - rab_y[i] * t_x_xz_0_x[i];

        t_xy_xy_0_z[i] = t_x_xyy_0_z[i] - rab_y[i] * t_x_xy_0_z[i];

        t_xy_xy_0_y[i] = t_x_xyy_0_y[i] - rab_y[i] * t_x_xy_0_y[i];

        t_xy_xy_0_x[i] = t_x_xyy_0_x[i] - rab_y[i] * t_x_xy_0_x[i];

        t_xy_xx_0_z[i] = t_x_xxy_0_z[i] - rab_y[i] * t_x_xx_0_z[i];

        t_xy_xx_0_y[i] = t_x_xxy_0_y[i] - rab_y[i] * t_x_xx_0_y[i];

        t_xy_xx_0_x[i] = t_x_xxy_0_x[i] - rab_y[i] * t_x_xx_0_x[i];

        t_xx_zz_0_z[i] = t_x_xzz_0_z[i] - rab_x[i] * t_x_zz_0_z[i];

        t_xx_zz_0_y[i] = t_x_xzz_0_y[i] - rab_x[i] * t_x_zz_0_y[i];

        t_xx_zz_0_x[i] = t_x_xzz_0_x[i] - rab_x[i] * t_x_zz_0_x[i];

        t_xx_yz_0_z[i] = t_x_xyz_0_z[i] - rab_x[i] * t_x_yz_0_z[i];

        t_xx_yz_0_y[i] = t_x_xyz_0_y[i] - rab_x[i] * t_x_yz_0_y[i];

        t_xx_yz_0_x[i] = t_x_xyz_0_x[i] - rab_x[i] * t_x_yz_0_x[i];

        t_xx_yy_0_z[i] = t_x_xyy_0_z[i] - rab_x[i] * t_x_yy_0_z[i];

        t_xx_yy_0_y[i] = t_x_xyy_0_y[i] - rab_x[i] * t_x_yy_0_y[i];

        t_xx_yy_0_x[i] = t_x_xyy_0_x[i] - rab_x[i] * t_x_yy_0_x[i];

        t_xx_xz_0_z[i] = t_x_xxz_0_z[i] - rab_x[i] * t_x_xz_0_z[i];

        t_xx_xz_0_y[i] = t_x_xxz_0_y[i] - rab_x[i] * t_x_xz_0_y[i];

        t_xx_xz_0_x[i] = t_x_xxz_0_x[i] - rab_x[i] * t_x_xz_0_x[i];

        t_xx_xy_0_z[i] = t_x_xxy_0_z[i] - rab_x[i] * t_x_xy_0_z[i];

        t_xx_xy_0_y[i] = t_x_xxy_0_y[i] - rab_x[i] * t_x_xy_0_y[i];

        t_xx_xy_0_x[i] = t_x_xxy_0_x[i] - rab_x[i] * t_x_xy_0_x[i];

        t_xx_xx_0_z[i] = t_x_xxx_0_z[i] - rab_x[i] * t_x_xx_0_z[i];

        t_xx_xx_0_y[i] = t_x_xxx_0_y[i] - rab_x[i] * t_x_xx_0_y[i];

        t_xx_xx_0_x[i] = t_x_xxx_0_x[i] - rab_x[i] * t_x_xx_0_x[i];
    }
}


} // derirec namespace
