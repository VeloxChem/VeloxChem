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
compHostHRRForPFSP_V0(      BufferHostXY<T>&      intsBufferPFSP,
                      const BufferHostX<int32_t>& intsIndexesPFSP,
                      const BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP,
                      const BufferHostXY<T>&      intsBufferSGSP,
                      const BufferHostX<int32_t>& intsIndexesSGSP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SFSP) integral components

    t_0_zzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(0));

    t_0_zzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(1));

    t_0_zzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(2));

    t_0_yzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(3));

    t_0_yzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(4));

    t_0_yzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(5));

    t_0_yyz_0_z = intsBufferSFSP.data(intsIndexesSFSP(6));

    t_0_yyz_0_y = intsBufferSFSP.data(intsIndexesSFSP(7));

    t_0_yyz_0_x = intsBufferSFSP.data(intsIndexesSFSP(8));

    t_0_yyy_0_z = intsBufferSFSP.data(intsIndexesSFSP(9));

    t_0_yyy_0_y = intsBufferSFSP.data(intsIndexesSFSP(10));

    t_0_yyy_0_x = intsBufferSFSP.data(intsIndexesSFSP(11));

    t_0_xzz_0_z = intsBufferSFSP.data(intsIndexesSFSP(12));

    t_0_xzz_0_y = intsBufferSFSP.data(intsIndexesSFSP(13));

    t_0_xzz_0_x = intsBufferSFSP.data(intsIndexesSFSP(14));

    t_0_xyz_0_z = intsBufferSFSP.data(intsIndexesSFSP(15));

    t_0_xyz_0_y = intsBufferSFSP.data(intsIndexesSFSP(16));

    t_0_xyz_0_x = intsBufferSFSP.data(intsIndexesSFSP(17));

    t_0_xyy_0_z = intsBufferSFSP.data(intsIndexesSFSP(18));

    t_0_xyy_0_y = intsBufferSFSP.data(intsIndexesSFSP(19));

    t_0_xyy_0_x = intsBufferSFSP.data(intsIndexesSFSP(20));

    t_0_xxz_0_z = intsBufferSFSP.data(intsIndexesSFSP(21));

    t_0_xxz_0_y = intsBufferSFSP.data(intsIndexesSFSP(22));

    t_0_xxz_0_x = intsBufferSFSP.data(intsIndexesSFSP(23));

    t_0_xxy_0_z = intsBufferSFSP.data(intsIndexesSFSP(24));

    t_0_xxy_0_y = intsBufferSFSP.data(intsIndexesSFSP(25));

    t_0_xxy_0_x = intsBufferSFSP.data(intsIndexesSFSP(26));

    t_0_xxx_0_z = intsBufferSFSP.data(intsIndexesSFSP(27));

    t_0_xxx_0_y = intsBufferSFSP.data(intsIndexesSFSP(28));

    t_0_xxx_0_x = intsBufferSFSP.data(intsIndexesSFSP(29));

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

    #pragma omp simd align(rab_y, rab_z, t_0_xxz_0_x, t_0_xxz_0_y, t_0_xxz_0_z, t_0_xxzz_0_x,\
                           t_0_xxzz_0_y, t_0_xxzz_0_z, t_0_xyyz_0_x, t_0_xyyz_0_y, t_0_xyyz_0_z,\
                           t_0_xyz_0_x, t_0_xyz_0_y, t_0_xyz_0_z, t_0_xyzz_0_x, t_0_xyzz_0_y,\
                           t_0_xyzz_0_z, t_0_xzz_0_x, t_0_xzz_0_y, t_0_xzz_0_z, t_0_xzzz_0_x,\
                           t_0_xzzz_0_y, t_0_xzzz_0_z, t_0_yyy_0_x, t_0_yyy_0_y, t_0_yyy_0_z,\
                           t_0_yyyy_0_x, t_0_yyyy_0_y, t_0_yyyy_0_z, t_0_yyyz_0_x, t_0_yyyz_0_y,\
                           t_0_yyyz_0_z, t_0_yyz_0_x, t_0_yyz_0_y, t_0_yyz_0_z, t_0_yyzz_0_x,\
                           t_0_yyzz_0_y, t_0_yyzz_0_z, t_0_yzz_0_x, t_0_yzz_0_y, t_0_yzz_0_z,\
                           t_0_yzzz_0_x, t_0_yzzz_0_y, t_0_yzzz_0_z, t_0_zzz_0_x, t_0_zzz_0_y,\
                           t_0_zzz_0_z, t_0_zzzz_0_x, t_0_zzzz_0_y, t_0_zzzz_0_z, t_y_xyz_0_x,\
                           t_y_xyz_0_y, t_y_xyz_0_z, t_y_xzz_0_x, t_y_xzz_0_y, t_y_xzz_0_z,\
                           t_y_yyy_0_x, t_y_yyy_0_y, t_y_yyy_0_z, t_y_yyz_0_x, t_y_yyz_0_y,\
                           t_y_yyz_0_z, t_y_yzz_0_x, t_y_yzz_0_y, t_y_yzz_0_z, t_y_zzz_0_x,\
                           t_y_zzz_0_y, t_y_zzz_0_z, t_z_xxz_0_x, t_z_xxz_0_y, t_z_xxz_0_z,\
                           t_z_xyz_0_x, t_z_xyz_0_y, t_z_xyz_0_z, t_z_xzz_0_x, t_z_xzz_0_y,\
                           t_z_xzz_0_z, t_z_yyz_0_x, t_z_yyz_0_y, t_z_yyz_0_z, t_z_yzz_0_x,\
                           t_z_yzz_0_y, t_z_yzz_0_z, t_z_zzz_0_x, t_z_zzz_0_y, t_z_zzz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_0_z[i] = t_0_zzzz_0_z[i] - rab_z[i] * t_0_zzz_0_z[i];

        t_z_zzz_0_y[i] = t_0_zzzz_0_y[i] - rab_z[i] * t_0_zzz_0_y[i];

        t_z_zzz_0_x[i] = t_0_zzzz_0_x[i] - rab_z[i] * t_0_zzz_0_x[i];

        t_z_yzz_0_z[i] = t_0_yzzz_0_z[i] - rab_z[i] * t_0_yzz_0_z[i];

        t_z_yzz_0_y[i] = t_0_yzzz_0_y[i] - rab_z[i] * t_0_yzz_0_y[i];

        t_z_yzz_0_x[i] = t_0_yzzz_0_x[i] - rab_z[i] * t_0_yzz_0_x[i];

        t_z_yyz_0_z[i] = t_0_yyzz_0_z[i] - rab_z[i] * t_0_yyz_0_z[i];

        t_z_yyz_0_y[i] = t_0_yyzz_0_y[i] - rab_z[i] * t_0_yyz_0_y[i];

        t_z_yyz_0_x[i] = t_0_yyzz_0_x[i] - rab_z[i] * t_0_yyz_0_x[i];

        t_z_xzz_0_z[i] = t_0_xzzz_0_z[i] - rab_z[i] * t_0_xzz_0_z[i];

        t_z_xzz_0_y[i] = t_0_xzzz_0_y[i] - rab_z[i] * t_0_xzz_0_y[i];

        t_z_xzz_0_x[i] = t_0_xzzz_0_x[i] - rab_z[i] * t_0_xzz_0_x[i];

        t_z_xyz_0_z[i] = t_0_xyzz_0_z[i] - rab_z[i] * t_0_xyz_0_z[i];

        t_z_xyz_0_y[i] = t_0_xyzz_0_y[i] - rab_z[i] * t_0_xyz_0_y[i];

        t_z_xyz_0_x[i] = t_0_xyzz_0_x[i] - rab_z[i] * t_0_xyz_0_x[i];

        t_z_xxz_0_z[i] = t_0_xxzz_0_z[i] - rab_z[i] * t_0_xxz_0_z[i];

        t_z_xxz_0_y[i] = t_0_xxzz_0_y[i] - rab_z[i] * t_0_xxz_0_y[i];

        t_z_xxz_0_x[i] = t_0_xxzz_0_x[i] - rab_z[i] * t_0_xxz_0_x[i];

        t_y_zzz_0_z[i] = t_0_yzzz_0_z[i] - rab_y[i] * t_0_zzz_0_z[i];

        t_y_zzz_0_y[i] = t_0_yzzz_0_y[i] - rab_y[i] * t_0_zzz_0_y[i];

        t_y_zzz_0_x[i] = t_0_yzzz_0_x[i] - rab_y[i] * t_0_zzz_0_x[i];

        t_y_yzz_0_z[i] = t_0_yyzz_0_z[i] - rab_y[i] * t_0_yzz_0_z[i];

        t_y_yzz_0_y[i] = t_0_yyzz_0_y[i] - rab_y[i] * t_0_yzz_0_y[i];

        t_y_yzz_0_x[i] = t_0_yyzz_0_x[i] - rab_y[i] * t_0_yzz_0_x[i];

        t_y_yyz_0_z[i] = t_0_yyyz_0_z[i] - rab_y[i] * t_0_yyz_0_z[i];

        t_y_yyz_0_y[i] = t_0_yyyz_0_y[i] - rab_y[i] * t_0_yyz_0_y[i];

        t_y_yyz_0_x[i] = t_0_yyyz_0_x[i] - rab_y[i] * t_0_yyz_0_x[i];

        t_y_yyy_0_z[i] = t_0_yyyy_0_z[i] - rab_y[i] * t_0_yyy_0_z[i];

        t_y_yyy_0_y[i] = t_0_yyyy_0_y[i] - rab_y[i] * t_0_yyy_0_y[i];

        t_y_yyy_0_x[i] = t_0_yyyy_0_x[i] - rab_y[i] * t_0_yyy_0_x[i];

        t_y_xzz_0_z[i] = t_0_xyzz_0_z[i] - rab_y[i] * t_0_xzz_0_z[i];

        t_y_xzz_0_y[i] = t_0_xyzz_0_y[i] - rab_y[i] * t_0_xzz_0_y[i];

        t_y_xzz_0_x[i] = t_0_xyzz_0_x[i] - rab_y[i] * t_0_xzz_0_x[i];

        t_y_xyz_0_z[i] = t_0_xyyz_0_z[i] - rab_y[i] * t_0_xyz_0_z[i];

        t_y_xyz_0_y[i] = t_0_xyyz_0_y[i] - rab_y[i] * t_0_xyz_0_y[i];

        t_y_xyz_0_x[i] = t_0_xyyz_0_x[i] - rab_y[i] * t_0_xyz_0_x[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_0_xxxy_0_x, t_0_xxxy_0_y, t_0_xxxy_0_z, t_0_xxxz_0_x,\
                           t_0_xxxz_0_y, t_0_xxxz_0_z, t_0_xxy_0_x, t_0_xxy_0_y, t_0_xxy_0_z,\
                           t_0_xxyy_0_x, t_0_xxyy_0_y, t_0_xxyy_0_z, t_0_xxyz_0_x, t_0_xxyz_0_y,\
                           t_0_xxyz_0_z, t_0_xxz_0_x, t_0_xxz_0_y, t_0_xxz_0_z, t_0_xxzz_0_x,\
                           t_0_xxzz_0_y, t_0_xxzz_0_z, t_0_xyy_0_x, t_0_xyy_0_y, t_0_xyy_0_z,\
                           t_0_xyyy_0_x, t_0_xyyy_0_y, t_0_xyyy_0_z, t_0_xyyz_0_x, t_0_xyyz_0_y,\
                           t_0_xyyz_0_z, t_0_xyz_0_x, t_0_xyz_0_y, t_0_xyz_0_z, t_0_xyzz_0_x,\
                           t_0_xyzz_0_y, t_0_xyzz_0_z, t_0_xzz_0_x, t_0_xzz_0_y, t_0_xzz_0_z,\
                           t_0_xzzz_0_x, t_0_xzzz_0_y, t_0_xzzz_0_z, t_0_yyy_0_x, t_0_yyy_0_y,\
                           t_0_yyy_0_z, t_0_yyz_0_x, t_0_yyz_0_y, t_0_yyz_0_z, t_0_yzz_0_x,\
                           t_0_yzz_0_y, t_0_yzz_0_z, t_0_zzz_0_x, t_0_zzz_0_y, t_0_zzz_0_z,\
                           t_x_xxy_0_x, t_x_xxy_0_y, t_x_xxy_0_z, t_x_xxz_0_x, t_x_xxz_0_y,\
                           t_x_xxz_0_z, t_x_xyy_0_x, t_x_xyy_0_y, t_x_xyy_0_z, t_x_xyz_0_x,\
                           t_x_xyz_0_y, t_x_xyz_0_z, t_x_xzz_0_x, t_x_xzz_0_y, t_x_xzz_0_z,\
                           t_x_yyy_0_x, t_x_yyy_0_y, t_x_yyy_0_z, t_x_yyz_0_x, t_x_yyz_0_y,\
                           t_x_yyz_0_z, t_x_yzz_0_x, t_x_yzz_0_y, t_x_yzz_0_z, t_x_zzz_0_x,\
                           t_x_zzz_0_y, t_x_zzz_0_z, t_y_xxy_0_x, t_y_xxy_0_y, t_y_xxy_0_z,\
                           t_y_xxz_0_x, t_y_xxz_0_y, t_y_xxz_0_z, t_y_xyy_0_x, t_y_xyy_0_y,\
                           t_y_xyy_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_xyy_0_z[i] = t_0_xyyy_0_z[i] - rab_y[i] * t_0_xyy_0_z[i];

        t_y_xyy_0_y[i] = t_0_xyyy_0_y[i] - rab_y[i] * t_0_xyy_0_y[i];

        t_y_xyy_0_x[i] = t_0_xyyy_0_x[i] - rab_y[i] * t_0_xyy_0_x[i];

        t_y_xxz_0_z[i] = t_0_xxyz_0_z[i] - rab_y[i] * t_0_xxz_0_z[i];

        t_y_xxz_0_y[i] = t_0_xxyz_0_y[i] - rab_y[i] * t_0_xxz_0_y[i];

        t_y_xxz_0_x[i] = t_0_xxyz_0_x[i] - rab_y[i] * t_0_xxz_0_x[i];

        t_y_xxy_0_z[i] = t_0_xxyy_0_z[i] - rab_y[i] * t_0_xxy_0_z[i];

        t_y_xxy_0_y[i] = t_0_xxyy_0_y[i] - rab_y[i] * t_0_xxy_0_y[i];

        t_y_xxy_0_x[i] = t_0_xxyy_0_x[i] - rab_y[i] * t_0_xxy_0_x[i];

        t_x_zzz_0_z[i] = t_0_xzzz_0_z[i] - rab_x[i] * t_0_zzz_0_z[i];

        t_x_zzz_0_y[i] = t_0_xzzz_0_y[i] - rab_x[i] * t_0_zzz_0_y[i];

        t_x_zzz_0_x[i] = t_0_xzzz_0_x[i] - rab_x[i] * t_0_zzz_0_x[i];

        t_x_yzz_0_z[i] = t_0_xyzz_0_z[i] - rab_x[i] * t_0_yzz_0_z[i];

        t_x_yzz_0_y[i] = t_0_xyzz_0_y[i] - rab_x[i] * t_0_yzz_0_y[i];

        t_x_yzz_0_x[i] = t_0_xyzz_0_x[i] - rab_x[i] * t_0_yzz_0_x[i];

        t_x_yyz_0_z[i] = t_0_xyyz_0_z[i] - rab_x[i] * t_0_yyz_0_z[i];

        t_x_yyz_0_y[i] = t_0_xyyz_0_y[i] - rab_x[i] * t_0_yyz_0_y[i];

        t_x_yyz_0_x[i] = t_0_xyyz_0_x[i] - rab_x[i] * t_0_yyz_0_x[i];

        t_x_yyy_0_z[i] = t_0_xyyy_0_z[i] - rab_x[i] * t_0_yyy_0_z[i];

        t_x_yyy_0_y[i] = t_0_xyyy_0_y[i] - rab_x[i] * t_0_yyy_0_y[i];

        t_x_yyy_0_x[i] = t_0_xyyy_0_x[i] - rab_x[i] * t_0_yyy_0_x[i];

        t_x_xzz_0_z[i] = t_0_xxzz_0_z[i] - rab_x[i] * t_0_xzz_0_z[i];

        t_x_xzz_0_y[i] = t_0_xxzz_0_y[i] - rab_x[i] * t_0_xzz_0_y[i];

        t_x_xzz_0_x[i] = t_0_xxzz_0_x[i] - rab_x[i] * t_0_xzz_0_x[i];

        t_x_xyz_0_z[i] = t_0_xxyz_0_z[i] - rab_x[i] * t_0_xyz_0_z[i];

        t_x_xyz_0_y[i] = t_0_xxyz_0_y[i] - rab_x[i] * t_0_xyz_0_y[i];

        t_x_xyz_0_x[i] = t_0_xxyz_0_x[i] - rab_x[i] * t_0_xyz_0_x[i];

        t_x_xyy_0_z[i] = t_0_xxyy_0_z[i] - rab_x[i] * t_0_xyy_0_z[i];

        t_x_xyy_0_y[i] = t_0_xxyy_0_y[i] - rab_x[i] * t_0_xyy_0_y[i];

        t_x_xyy_0_x[i] = t_0_xxyy_0_x[i] - rab_x[i] * t_0_xyy_0_x[i];

        t_x_xxz_0_z[i] = t_0_xxxz_0_z[i] - rab_x[i] * t_0_xxz_0_z[i];

        t_x_xxz_0_y[i] = t_0_xxxz_0_y[i] - rab_x[i] * t_0_xxz_0_y[i];

        t_x_xxz_0_x[i] = t_0_xxxz_0_x[i] - rab_x[i] * t_0_xxz_0_x[i];

        t_x_xxy_0_z[i] = t_0_xxxy_0_z[i] - rab_x[i] * t_0_xxy_0_z[i];

        t_x_xxy_0_y[i] = t_0_xxxy_0_y[i] - rab_x[i] * t_0_xxy_0_y[i];

        t_x_xxy_0_x[i] = t_0_xxxy_0_x[i] - rab_x[i] * t_0_xxy_0_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xxx_0_x, t_0_xxx_0_y, t_0_xxx_0_z, t_0_xxxx_0_x,\
                           t_0_xxxx_0_y, t_0_xxxx_0_z, t_x_xxx_0_x, t_x_xxx_0_y, t_x_xxx_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_xxx_0_z[i] = t_0_xxxx_0_z[i] - rab_x[i] * t_0_xxx_0_z[i];

        t_x_xxx_0_y[i] = t_0_xxxx_0_y[i] - rab_x[i] * t_0_xxx_0_y[i];

        t_x_xxx_0_x[i] = t_0_xxxx_0_x[i] - rab_x[i] * t_0_xxx_0_x[i];
    }
}


} // derirec namespace
