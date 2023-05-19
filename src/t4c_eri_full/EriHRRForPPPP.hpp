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
compHostHRRForPPPP_V0(      BufferHostXY<T>&      intsBufferPPPP,
                      const BufferHostX<int32_t>& intsIndexesPPPP,
                      const BufferHostXY<T>&      intsBufferSPPP,
                      const BufferHostX<int32_t>& intsIndexesSPPP,
                      const BufferHostXY<T>&      intsBufferSDPP,
                      const BufferHostX<int32_t>& intsIndexesSDPP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPPP) integral components

    t_z_z_z_z = intsBufferPPPP.data(intsIndexesPPPP(0));

    t_z_z_z_y = intsBufferPPPP.data(intsIndexesPPPP(1));

    t_z_z_z_x = intsBufferPPPP.data(intsIndexesPPPP(2));

    t_z_z_y_z = intsBufferPPPP.data(intsIndexesPPPP(3));

    t_z_z_y_y = intsBufferPPPP.data(intsIndexesPPPP(4));

    t_z_z_y_x = intsBufferPPPP.data(intsIndexesPPPP(5));

    t_z_z_x_z = intsBufferPPPP.data(intsIndexesPPPP(6));

    t_z_z_x_y = intsBufferPPPP.data(intsIndexesPPPP(7));

    t_z_z_x_x = intsBufferPPPP.data(intsIndexesPPPP(8));

    t_z_y_z_z = intsBufferPPPP.data(intsIndexesPPPP(9));

    t_z_y_z_y = intsBufferPPPP.data(intsIndexesPPPP(10));

    t_z_y_z_x = intsBufferPPPP.data(intsIndexesPPPP(11));

    t_z_y_y_z = intsBufferPPPP.data(intsIndexesPPPP(12));

    t_z_y_y_y = intsBufferPPPP.data(intsIndexesPPPP(13));

    t_z_y_y_x = intsBufferPPPP.data(intsIndexesPPPP(14));

    t_z_y_x_z = intsBufferPPPP.data(intsIndexesPPPP(15));

    t_z_y_x_y = intsBufferPPPP.data(intsIndexesPPPP(16));

    t_z_y_x_x = intsBufferPPPP.data(intsIndexesPPPP(17));

    t_z_x_z_z = intsBufferPPPP.data(intsIndexesPPPP(18));

    t_z_x_z_y = intsBufferPPPP.data(intsIndexesPPPP(19));

    t_z_x_z_x = intsBufferPPPP.data(intsIndexesPPPP(20));

    t_z_x_y_z = intsBufferPPPP.data(intsIndexesPPPP(21));

    t_z_x_y_y = intsBufferPPPP.data(intsIndexesPPPP(22));

    t_z_x_y_x = intsBufferPPPP.data(intsIndexesPPPP(23));

    t_z_x_x_z = intsBufferPPPP.data(intsIndexesPPPP(24));

    t_z_x_x_y = intsBufferPPPP.data(intsIndexesPPPP(25));

    t_z_x_x_x = intsBufferPPPP.data(intsIndexesPPPP(26));

    t_y_z_z_z = intsBufferPPPP.data(intsIndexesPPPP(27));

    t_y_z_z_y = intsBufferPPPP.data(intsIndexesPPPP(28));

    t_y_z_z_x = intsBufferPPPP.data(intsIndexesPPPP(29));

    t_y_z_y_z = intsBufferPPPP.data(intsIndexesPPPP(30));

    t_y_z_y_y = intsBufferPPPP.data(intsIndexesPPPP(31));

    t_y_z_y_x = intsBufferPPPP.data(intsIndexesPPPP(32));

    t_y_z_x_z = intsBufferPPPP.data(intsIndexesPPPP(33));

    t_y_z_x_y = intsBufferPPPP.data(intsIndexesPPPP(34));

    t_y_z_x_x = intsBufferPPPP.data(intsIndexesPPPP(35));

    t_y_y_z_z = intsBufferPPPP.data(intsIndexesPPPP(36));

    t_y_y_z_y = intsBufferPPPP.data(intsIndexesPPPP(37));

    t_y_y_z_x = intsBufferPPPP.data(intsIndexesPPPP(38));

    t_y_y_y_z = intsBufferPPPP.data(intsIndexesPPPP(39));

    t_y_y_y_y = intsBufferPPPP.data(intsIndexesPPPP(40));

    t_y_y_y_x = intsBufferPPPP.data(intsIndexesPPPP(41));

    t_y_y_x_z = intsBufferPPPP.data(intsIndexesPPPP(42));

    t_y_y_x_y = intsBufferPPPP.data(intsIndexesPPPP(43));

    t_y_y_x_x = intsBufferPPPP.data(intsIndexesPPPP(44));

    t_y_x_z_z = intsBufferPPPP.data(intsIndexesPPPP(45));

    t_y_x_z_y = intsBufferPPPP.data(intsIndexesPPPP(46));

    t_y_x_z_x = intsBufferPPPP.data(intsIndexesPPPP(47));

    t_y_x_y_z = intsBufferPPPP.data(intsIndexesPPPP(48));

    t_y_x_y_y = intsBufferPPPP.data(intsIndexesPPPP(49));

    t_y_x_y_x = intsBufferPPPP.data(intsIndexesPPPP(50));

    t_y_x_x_z = intsBufferPPPP.data(intsIndexesPPPP(51));

    t_y_x_x_y = intsBufferPPPP.data(intsIndexesPPPP(52));

    t_y_x_x_x = intsBufferPPPP.data(intsIndexesPPPP(53));

    t_x_z_z_z = intsBufferPPPP.data(intsIndexesPPPP(54));

    t_x_z_z_y = intsBufferPPPP.data(intsIndexesPPPP(55));

    t_x_z_z_x = intsBufferPPPP.data(intsIndexesPPPP(56));

    t_x_z_y_z = intsBufferPPPP.data(intsIndexesPPPP(57));

    t_x_z_y_y = intsBufferPPPP.data(intsIndexesPPPP(58));

    t_x_z_y_x = intsBufferPPPP.data(intsIndexesPPPP(59));

    t_x_z_x_z = intsBufferPPPP.data(intsIndexesPPPP(60));

    t_x_z_x_y = intsBufferPPPP.data(intsIndexesPPPP(61));

    t_x_z_x_x = intsBufferPPPP.data(intsIndexesPPPP(62));

    t_x_y_z_z = intsBufferPPPP.data(intsIndexesPPPP(63));

    t_x_y_z_y = intsBufferPPPP.data(intsIndexesPPPP(64));

    t_x_y_z_x = intsBufferPPPP.data(intsIndexesPPPP(65));

    t_x_y_y_z = intsBufferPPPP.data(intsIndexesPPPP(66));

    t_x_y_y_y = intsBufferPPPP.data(intsIndexesPPPP(67));

    t_x_y_y_x = intsBufferPPPP.data(intsIndexesPPPP(68));

    t_x_y_x_z = intsBufferPPPP.data(intsIndexesPPPP(69));

    t_x_y_x_y = intsBufferPPPP.data(intsIndexesPPPP(70));

    t_x_y_x_x = intsBufferPPPP.data(intsIndexesPPPP(71));

    t_x_x_z_z = intsBufferPPPP.data(intsIndexesPPPP(72));

    t_x_x_z_y = intsBufferPPPP.data(intsIndexesPPPP(73));

    t_x_x_z_x = intsBufferPPPP.data(intsIndexesPPPP(74));

    t_x_x_y_z = intsBufferPPPP.data(intsIndexesPPPP(75));

    t_x_x_y_y = intsBufferPPPP.data(intsIndexesPPPP(76));

    t_x_x_y_x = intsBufferPPPP.data(intsIndexesPPPP(77));

    t_x_x_x_z = intsBufferPPPP.data(intsIndexesPPPP(78));

    t_x_x_x_y = intsBufferPPPP.data(intsIndexesPPPP(79));

    t_x_x_x_x = intsBufferPPPP.data(intsIndexesPPPP(80));

    // set up (SPPP) integral components

    t_0_z_z_z = intsBufferSPPP.data(intsIndexesSPPP(0));

    t_0_z_z_y = intsBufferSPPP.data(intsIndexesSPPP(1));

    t_0_z_z_x = intsBufferSPPP.data(intsIndexesSPPP(2));

    t_0_z_y_z = intsBufferSPPP.data(intsIndexesSPPP(3));

    t_0_z_y_y = intsBufferSPPP.data(intsIndexesSPPP(4));

    t_0_z_y_x = intsBufferSPPP.data(intsIndexesSPPP(5));

    t_0_z_x_z = intsBufferSPPP.data(intsIndexesSPPP(6));

    t_0_z_x_y = intsBufferSPPP.data(intsIndexesSPPP(7));

    t_0_z_x_x = intsBufferSPPP.data(intsIndexesSPPP(8));

    t_0_y_z_z = intsBufferSPPP.data(intsIndexesSPPP(9));

    t_0_y_z_y = intsBufferSPPP.data(intsIndexesSPPP(10));

    t_0_y_z_x = intsBufferSPPP.data(intsIndexesSPPP(11));

    t_0_y_y_z = intsBufferSPPP.data(intsIndexesSPPP(12));

    t_0_y_y_y = intsBufferSPPP.data(intsIndexesSPPP(13));

    t_0_y_y_x = intsBufferSPPP.data(intsIndexesSPPP(14));

    t_0_y_x_z = intsBufferSPPP.data(intsIndexesSPPP(15));

    t_0_y_x_y = intsBufferSPPP.data(intsIndexesSPPP(16));

    t_0_y_x_x = intsBufferSPPP.data(intsIndexesSPPP(17));

    t_0_x_z_z = intsBufferSPPP.data(intsIndexesSPPP(18));

    t_0_x_z_y = intsBufferSPPP.data(intsIndexesSPPP(19));

    t_0_x_z_x = intsBufferSPPP.data(intsIndexesSPPP(20));

    t_0_x_y_z = intsBufferSPPP.data(intsIndexesSPPP(21));

    t_0_x_y_y = intsBufferSPPP.data(intsIndexesSPPP(22));

    t_0_x_y_x = intsBufferSPPP.data(intsIndexesSPPP(23));

    t_0_x_x_z = intsBufferSPPP.data(intsIndexesSPPP(24));

    t_0_x_x_y = intsBufferSPPP.data(intsIndexesSPPP(25));

    t_0_x_x_x = intsBufferSPPP.data(intsIndexesSPPP(26));

    // set up (SDPP) integral components

    t_0_zz_z_z = intsBufferSDPP.data(intsIndexesSDPP(0));

    t_0_zz_z_y = intsBufferSDPP.data(intsIndexesSDPP(1));

    t_0_zz_z_x = intsBufferSDPP.data(intsIndexesSDPP(2));

    t_0_zz_y_z = intsBufferSDPP.data(intsIndexesSDPP(3));

    t_0_zz_y_y = intsBufferSDPP.data(intsIndexesSDPP(4));

    t_0_zz_y_x = intsBufferSDPP.data(intsIndexesSDPP(5));

    t_0_zz_x_z = intsBufferSDPP.data(intsIndexesSDPP(6));

    t_0_zz_x_y = intsBufferSDPP.data(intsIndexesSDPP(7));

    t_0_zz_x_x = intsBufferSDPP.data(intsIndexesSDPP(8));

    t_0_yz_z_z = intsBufferSDPP.data(intsIndexesSDPP(9));

    t_0_yz_z_y = intsBufferSDPP.data(intsIndexesSDPP(10));

    t_0_yz_z_x = intsBufferSDPP.data(intsIndexesSDPP(11));

    t_0_yz_y_z = intsBufferSDPP.data(intsIndexesSDPP(12));

    t_0_yz_y_y = intsBufferSDPP.data(intsIndexesSDPP(13));

    t_0_yz_y_x = intsBufferSDPP.data(intsIndexesSDPP(14));

    t_0_yz_x_z = intsBufferSDPP.data(intsIndexesSDPP(15));

    t_0_yz_x_y = intsBufferSDPP.data(intsIndexesSDPP(16));

    t_0_yz_x_x = intsBufferSDPP.data(intsIndexesSDPP(17));

    t_0_yy_z_z = intsBufferSDPP.data(intsIndexesSDPP(18));

    t_0_yy_z_y = intsBufferSDPP.data(intsIndexesSDPP(19));

    t_0_yy_z_x = intsBufferSDPP.data(intsIndexesSDPP(20));

    t_0_yy_y_z = intsBufferSDPP.data(intsIndexesSDPP(21));

    t_0_yy_y_y = intsBufferSDPP.data(intsIndexesSDPP(22));

    t_0_yy_y_x = intsBufferSDPP.data(intsIndexesSDPP(23));

    t_0_yy_x_z = intsBufferSDPP.data(intsIndexesSDPP(24));

    t_0_yy_x_y = intsBufferSDPP.data(intsIndexesSDPP(25));

    t_0_yy_x_x = intsBufferSDPP.data(intsIndexesSDPP(26));

    t_0_xz_z_z = intsBufferSDPP.data(intsIndexesSDPP(27));

    t_0_xz_z_y = intsBufferSDPP.data(intsIndexesSDPP(28));

    t_0_xz_z_x = intsBufferSDPP.data(intsIndexesSDPP(29));

    t_0_xz_y_z = intsBufferSDPP.data(intsIndexesSDPP(30));

    t_0_xz_y_y = intsBufferSDPP.data(intsIndexesSDPP(31));

    t_0_xz_y_x = intsBufferSDPP.data(intsIndexesSDPP(32));

    t_0_xz_x_z = intsBufferSDPP.data(intsIndexesSDPP(33));

    t_0_xz_x_y = intsBufferSDPP.data(intsIndexesSDPP(34));

    t_0_xz_x_x = intsBufferSDPP.data(intsIndexesSDPP(35));

    t_0_xy_z_z = intsBufferSDPP.data(intsIndexesSDPP(36));

    t_0_xy_z_y = intsBufferSDPP.data(intsIndexesSDPP(37));

    t_0_xy_z_x = intsBufferSDPP.data(intsIndexesSDPP(38));

    t_0_xy_y_z = intsBufferSDPP.data(intsIndexesSDPP(39));

    t_0_xy_y_y = intsBufferSDPP.data(intsIndexesSDPP(40));

    t_0_xy_y_x = intsBufferSDPP.data(intsIndexesSDPP(41));

    t_0_xy_x_z = intsBufferSDPP.data(intsIndexesSDPP(42));

    t_0_xy_x_y = intsBufferSDPP.data(intsIndexesSDPP(43));

    t_0_xy_x_x = intsBufferSDPP.data(intsIndexesSDPP(44));

    t_0_xx_z_z = intsBufferSDPP.data(intsIndexesSDPP(45));

    t_0_xx_z_y = intsBufferSDPP.data(intsIndexesSDPP(46));

    t_0_xx_z_x = intsBufferSDPP.data(intsIndexesSDPP(47));

    t_0_xx_y_z = intsBufferSDPP.data(intsIndexesSDPP(48));

    t_0_xx_y_y = intsBufferSDPP.data(intsIndexesSDPP(49));

    t_0_xx_y_x = intsBufferSDPP.data(intsIndexesSDPP(50));

    t_0_xx_x_z = intsBufferSDPP.data(intsIndexesSDPP(51));

    t_0_xx_x_y = intsBufferSDPP.data(intsIndexesSDPP(52));

    t_0_xx_x_x = intsBufferSDPP.data(intsIndexesSDPP(53));

    #pragma omp simd align(rab_y, rab_z, t_0_x_x_x, t_0_x_x_y, t_0_x_x_z, t_0_x_y_x,\
                           t_0_x_y_y, t_0_x_y_z, t_0_x_z_x, t_0_x_z_y, t_0_x_z_z, t_0_xz_x_x,\
                           t_0_xz_x_y, t_0_xz_x_z, t_0_xz_y_x, t_0_xz_y_y, t_0_xz_y_z,\
                           t_0_xz_z_x, t_0_xz_z_y, t_0_xz_z_z, t_0_y_x_x, t_0_y_x_y,\
                           t_0_y_x_z, t_0_y_y_x, t_0_y_y_y, t_0_y_y_z, t_0_y_z_x, t_0_y_z_y,\
                           t_0_y_z_z, t_0_yz_x_x, t_0_yz_x_y, t_0_yz_x_z, t_0_yz_y_x,\
                           t_0_yz_y_y, t_0_yz_y_z, t_0_yz_z_x, t_0_yz_z_y, t_0_yz_z_z,\
                           t_0_z_x_x, t_0_z_x_y, t_0_z_x_z, t_0_z_y_x, t_0_z_y_y, t_0_z_y_z,\
                           t_0_z_z_x, t_0_z_z_y, t_0_z_z_z, t_0_zz_x_x, t_0_zz_x_y,\
                           t_0_zz_x_z, t_0_zz_y_x, t_0_zz_y_y, t_0_zz_y_z, t_0_zz_z_x,\
                           t_0_zz_z_y, t_0_zz_z_z, t_y_z_x_x, t_y_z_x_y, t_y_z_x_z,\
                           t_y_z_y_x, t_y_z_y_y, t_y_z_y_z, t_y_z_z_x, t_y_z_z_y, t_y_z_z_z,\
                           t_z_x_x_x, t_z_x_x_y, t_z_x_x_z, t_z_x_y_x, t_z_x_y_y, t_z_x_y_z,\
                           t_z_x_z_x, t_z_x_z_y, t_z_x_z_z, t_z_y_x_x, t_z_y_x_y, t_z_y_x_z,\
                           t_z_y_y_x, t_z_y_y_y, t_z_y_y_z, t_z_y_z_x, t_z_y_z_y, t_z_y_z_z,\
                           t_z_z_x_x, t_z_z_x_y, t_z_z_x_z, t_z_z_y_x, t_z_z_y_y, t_z_z_y_z,\
                           t_z_z_z_x, t_z_z_z_y, t_z_z_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_z_z[i] = t_0_zz_z_z[i] - rab_z[i] * t_0_z_z_z[i];

        t_z_z_z_y[i] = t_0_zz_z_y[i] - rab_z[i] * t_0_z_z_y[i];

        t_z_z_z_x[i] = t_0_zz_z_x[i] - rab_z[i] * t_0_z_z_x[i];

        t_z_z_y_z[i] = t_0_zz_y_z[i] - rab_z[i] * t_0_z_y_z[i];

        t_z_z_y_y[i] = t_0_zz_y_y[i] - rab_z[i] * t_0_z_y_y[i];

        t_z_z_y_x[i] = t_0_zz_y_x[i] - rab_z[i] * t_0_z_y_x[i];

        t_z_z_x_z[i] = t_0_zz_x_z[i] - rab_z[i] * t_0_z_x_z[i];

        t_z_z_x_y[i] = t_0_zz_x_y[i] - rab_z[i] * t_0_z_x_y[i];

        t_z_z_x_x[i] = t_0_zz_x_x[i] - rab_z[i] * t_0_z_x_x[i];

        t_z_y_z_z[i] = t_0_yz_z_z[i] - rab_z[i] * t_0_y_z_z[i];

        t_z_y_z_y[i] = t_0_yz_z_y[i] - rab_z[i] * t_0_y_z_y[i];

        t_z_y_z_x[i] = t_0_yz_z_x[i] - rab_z[i] * t_0_y_z_x[i];

        t_z_y_y_z[i] = t_0_yz_y_z[i] - rab_z[i] * t_0_y_y_z[i];

        t_z_y_y_y[i] = t_0_yz_y_y[i] - rab_z[i] * t_0_y_y_y[i];

        t_z_y_y_x[i] = t_0_yz_y_x[i] - rab_z[i] * t_0_y_y_x[i];

        t_z_y_x_z[i] = t_0_yz_x_z[i] - rab_z[i] * t_0_y_x_z[i];

        t_z_y_x_y[i] = t_0_yz_x_y[i] - rab_z[i] * t_0_y_x_y[i];

        t_z_y_x_x[i] = t_0_yz_x_x[i] - rab_z[i] * t_0_y_x_x[i];

        t_z_x_z_z[i] = t_0_xz_z_z[i] - rab_z[i] * t_0_x_z_z[i];

        t_z_x_z_y[i] = t_0_xz_z_y[i] - rab_z[i] * t_0_x_z_y[i];

        t_z_x_z_x[i] = t_0_xz_z_x[i] - rab_z[i] * t_0_x_z_x[i];

        t_z_x_y_z[i] = t_0_xz_y_z[i] - rab_z[i] * t_0_x_y_z[i];

        t_z_x_y_y[i] = t_0_xz_y_y[i] - rab_z[i] * t_0_x_y_y[i];

        t_z_x_y_x[i] = t_0_xz_y_x[i] - rab_z[i] * t_0_x_y_x[i];

        t_z_x_x_z[i] = t_0_xz_x_z[i] - rab_z[i] * t_0_x_x_z[i];

        t_z_x_x_y[i] = t_0_xz_x_y[i] - rab_z[i] * t_0_x_x_y[i];

        t_z_x_x_x[i] = t_0_xz_x_x[i] - rab_z[i] * t_0_x_x_x[i];

        t_y_z_z_z[i] = t_0_yz_z_z[i] - rab_y[i] * t_0_z_z_z[i];

        t_y_z_z_y[i] = t_0_yz_z_y[i] - rab_y[i] * t_0_z_z_y[i];

        t_y_z_z_x[i] = t_0_yz_z_x[i] - rab_y[i] * t_0_z_z_x[i];

        t_y_z_y_z[i] = t_0_yz_y_z[i] - rab_y[i] * t_0_z_y_z[i];

        t_y_z_y_y[i] = t_0_yz_y_y[i] - rab_y[i] * t_0_z_y_y[i];

        t_y_z_y_x[i] = t_0_yz_y_x[i] - rab_y[i] * t_0_z_y_x[i];

        t_y_z_x_z[i] = t_0_yz_x_z[i] - rab_y[i] * t_0_z_x_z[i];

        t_y_z_x_y[i] = t_0_yz_x_y[i] - rab_y[i] * t_0_z_x_y[i];

        t_y_z_x_x[i] = t_0_yz_x_x[i] - rab_y[i] * t_0_z_x_x[i];
    }

    #pragma omp simd align(rab_x, rab_y, t_0_x_x_x, t_0_x_x_y, t_0_x_x_z, t_0_x_y_x,\
                           t_0_x_y_y, t_0_x_y_z, t_0_x_z_x, t_0_x_z_y, t_0_x_z_z, t_0_xy_x_x,\
                           t_0_xy_x_y, t_0_xy_x_z, t_0_xy_y_x, t_0_xy_y_y, t_0_xy_y_z,\
                           t_0_xy_z_x, t_0_xy_z_y, t_0_xy_z_z, t_0_xz_x_x, t_0_xz_x_y,\
                           t_0_xz_x_z, t_0_xz_y_x, t_0_xz_y_y, t_0_xz_y_z, t_0_xz_z_x,\
                           t_0_xz_z_y, t_0_xz_z_z, t_0_y_x_x, t_0_y_x_y, t_0_y_x_z,\
                           t_0_y_y_x, t_0_y_y_y, t_0_y_y_z, t_0_y_z_x, t_0_y_z_y, t_0_y_z_z,\
                           t_0_yy_x_x, t_0_yy_x_y, t_0_yy_x_z, t_0_yy_y_x, t_0_yy_y_y,\
                           t_0_yy_y_z, t_0_yy_z_x, t_0_yy_z_y, t_0_yy_z_z, t_0_z_x_x,\
                           t_0_z_x_y, t_0_z_x_z, t_0_z_y_x, t_0_z_y_y, t_0_z_y_z, t_0_z_z_x,\
                           t_0_z_z_y, t_0_z_z_z, t_x_y_x_x, t_x_y_x_y, t_x_y_x_z, t_x_y_y_x,\
                           t_x_y_y_y, t_x_y_y_z, t_x_y_z_x, t_x_y_z_y, t_x_y_z_z, t_x_z_x_x,\
                           t_x_z_x_y, t_x_z_x_z, t_x_z_y_x, t_x_z_y_y, t_x_z_y_z, t_x_z_z_x,\
                           t_x_z_z_y, t_x_z_z_z, t_y_x_x_x, t_y_x_x_y, t_y_x_x_z, t_y_x_y_x,\
                           t_y_x_y_y, t_y_x_y_z, t_y_x_z_x, t_y_x_z_y, t_y_x_z_z, t_y_y_x_x,\
                           t_y_y_x_y, t_y_y_x_z, t_y_y_y_x, t_y_y_y_y, t_y_y_y_z, t_y_y_z_x,\
                           t_y_y_z_y, t_y_y_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_y_y_z_z[i] = t_0_yy_z_z[i] - rab_y[i] * t_0_y_z_z[i];

        t_y_y_z_y[i] = t_0_yy_z_y[i] - rab_y[i] * t_0_y_z_y[i];

        t_y_y_z_x[i] = t_0_yy_z_x[i] - rab_y[i] * t_0_y_z_x[i];

        t_y_y_y_z[i] = t_0_yy_y_z[i] - rab_y[i] * t_0_y_y_z[i];

        t_y_y_y_y[i] = t_0_yy_y_y[i] - rab_y[i] * t_0_y_y_y[i];

        t_y_y_y_x[i] = t_0_yy_y_x[i] - rab_y[i] * t_0_y_y_x[i];

        t_y_y_x_z[i] = t_0_yy_x_z[i] - rab_y[i] * t_0_y_x_z[i];

        t_y_y_x_y[i] = t_0_yy_x_y[i] - rab_y[i] * t_0_y_x_y[i];

        t_y_y_x_x[i] = t_0_yy_x_x[i] - rab_y[i] * t_0_y_x_x[i];

        t_y_x_z_z[i] = t_0_xy_z_z[i] - rab_y[i] * t_0_x_z_z[i];

        t_y_x_z_y[i] = t_0_xy_z_y[i] - rab_y[i] * t_0_x_z_y[i];

        t_y_x_z_x[i] = t_0_xy_z_x[i] - rab_y[i] * t_0_x_z_x[i];

        t_y_x_y_z[i] = t_0_xy_y_z[i] - rab_y[i] * t_0_x_y_z[i];

        t_y_x_y_y[i] = t_0_xy_y_y[i] - rab_y[i] * t_0_x_y_y[i];

        t_y_x_y_x[i] = t_0_xy_y_x[i] - rab_y[i] * t_0_x_y_x[i];

        t_y_x_x_z[i] = t_0_xy_x_z[i] - rab_y[i] * t_0_x_x_z[i];

        t_y_x_x_y[i] = t_0_xy_x_y[i] - rab_y[i] * t_0_x_x_y[i];

        t_y_x_x_x[i] = t_0_xy_x_x[i] - rab_y[i] * t_0_x_x_x[i];

        t_x_z_z_z[i] = t_0_xz_z_z[i] - rab_x[i] * t_0_z_z_z[i];

        t_x_z_z_y[i] = t_0_xz_z_y[i] - rab_x[i] * t_0_z_z_y[i];

        t_x_z_z_x[i] = t_0_xz_z_x[i] - rab_x[i] * t_0_z_z_x[i];

        t_x_z_y_z[i] = t_0_xz_y_z[i] - rab_x[i] * t_0_z_y_z[i];

        t_x_z_y_y[i] = t_0_xz_y_y[i] - rab_x[i] * t_0_z_y_y[i];

        t_x_z_y_x[i] = t_0_xz_y_x[i] - rab_x[i] * t_0_z_y_x[i];

        t_x_z_x_z[i] = t_0_xz_x_z[i] - rab_x[i] * t_0_z_x_z[i];

        t_x_z_x_y[i] = t_0_xz_x_y[i] - rab_x[i] * t_0_z_x_y[i];

        t_x_z_x_x[i] = t_0_xz_x_x[i] - rab_x[i] * t_0_z_x_x[i];

        t_x_y_z_z[i] = t_0_xy_z_z[i] - rab_x[i] * t_0_y_z_z[i];

        t_x_y_z_y[i] = t_0_xy_z_y[i] - rab_x[i] * t_0_y_z_y[i];

        t_x_y_z_x[i] = t_0_xy_z_x[i] - rab_x[i] * t_0_y_z_x[i];

        t_x_y_y_z[i] = t_0_xy_y_z[i] - rab_x[i] * t_0_y_y_z[i];

        t_x_y_y_y[i] = t_0_xy_y_y[i] - rab_x[i] * t_0_y_y_y[i];

        t_x_y_y_x[i] = t_0_xy_y_x[i] - rab_x[i] * t_0_y_y_x[i];

        t_x_y_x_z[i] = t_0_xy_x_z[i] - rab_x[i] * t_0_y_x_z[i];

        t_x_y_x_y[i] = t_0_xy_x_y[i] - rab_x[i] * t_0_y_x_y[i];

        t_x_y_x_x[i] = t_0_xy_x_x[i] - rab_x[i] * t_0_y_x_x[i];
    }

    #pragma omp simd align(rab_x, t_0_x_x_x, t_0_x_x_y, t_0_x_x_z, t_0_x_y_x, t_0_x_y_y,\
                           t_0_x_y_z, t_0_x_z_x, t_0_x_z_y, t_0_x_z_z, t_0_xx_x_x, t_0_xx_x_y,\
                           t_0_xx_x_z, t_0_xx_y_x, t_0_xx_y_y, t_0_xx_y_z, t_0_xx_z_x,\
                           t_0_xx_z_y, t_0_xx_z_z, t_x_x_x_x, t_x_x_x_y, t_x_x_x_z,\
                           t_x_x_y_x, t_x_x_y_y, t_x_x_y_z, t_x_x_z_x, t_x_x_z_y, t_x_x_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_x_z_z[i] = t_0_xx_z_z[i] - rab_x[i] * t_0_x_z_z[i];

        t_x_x_z_y[i] = t_0_xx_z_y[i] - rab_x[i] * t_0_x_z_y[i];

        t_x_x_z_x[i] = t_0_xx_z_x[i] - rab_x[i] * t_0_x_z_x[i];

        t_x_x_y_z[i] = t_0_xx_y_z[i] - rab_x[i] * t_0_x_y_z[i];

        t_x_x_y_y[i] = t_0_xx_y_y[i] - rab_x[i] * t_0_x_y_y[i];

        t_x_x_y_x[i] = t_0_xx_y_x[i] - rab_x[i] * t_0_x_y_x[i];

        t_x_x_x_z[i] = t_0_xx_x_z[i] - rab_x[i] * t_0_x_x_z[i];

        t_x_x_x_y[i] = t_0_xx_x_y[i] - rab_x[i] * t_0_x_x_y[i];

        t_x_x_x_x[i] = t_0_xx_x_x[i] - rab_x[i] * t_0_x_x_x[i];
    }
}


} // derirec namespace
