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
compHostHRRForSDPP_V0(      BufferHostXY<T>&      intsBufferSDPP,
                      const BufferHostX<int32_t>& intsIndexesSDPP,
                      const BufferHostXY<T>&      intsBufferSDSP,
                      const BufferHostX<int32_t>& intsIndexesSDSP,
                      const BufferHostXY<T>&      intsBufferSDSD,
                      const BufferHostX<int32_t>& intsIndexesSDSD,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

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

    // set up (SDSP) integral components

    t_0_zz_0_z = intsBufferSDSP.data(intsIndexesSDSP(0));

    t_0_zz_0_y = intsBufferSDSP.data(intsIndexesSDSP(1));

    t_0_zz_0_x = intsBufferSDSP.data(intsIndexesSDSP(2));

    t_0_yz_0_z = intsBufferSDSP.data(intsIndexesSDSP(3));

    t_0_yz_0_y = intsBufferSDSP.data(intsIndexesSDSP(4));

    t_0_yz_0_x = intsBufferSDSP.data(intsIndexesSDSP(5));

    t_0_yy_0_z = intsBufferSDSP.data(intsIndexesSDSP(6));

    t_0_yy_0_y = intsBufferSDSP.data(intsIndexesSDSP(7));

    t_0_yy_0_x = intsBufferSDSP.data(intsIndexesSDSP(8));

    t_0_xz_0_z = intsBufferSDSP.data(intsIndexesSDSP(9));

    t_0_xz_0_y = intsBufferSDSP.data(intsIndexesSDSP(10));

    t_0_xz_0_x = intsBufferSDSP.data(intsIndexesSDSP(11));

    t_0_xy_0_z = intsBufferSDSP.data(intsIndexesSDSP(12));

    t_0_xy_0_y = intsBufferSDSP.data(intsIndexesSDSP(13));

    t_0_xy_0_x = intsBufferSDSP.data(intsIndexesSDSP(14));

    t_0_xx_0_z = intsBufferSDSP.data(intsIndexesSDSP(15));

    t_0_xx_0_y = intsBufferSDSP.data(intsIndexesSDSP(16));

    t_0_xx_0_x = intsBufferSDSP.data(intsIndexesSDSP(17));

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

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xz_0_x, t_0_xz_0_xx, t_0_xz_0_xy,\
                           t_0_xz_0_xz, t_0_xz_0_y, t_0_xz_0_yy, t_0_xz_0_yz, t_0_xz_0_z,\
                           t_0_xz_0_zz, t_0_xz_x_x, t_0_xz_x_y, t_0_xz_x_z, t_0_xz_y_x,\
                           t_0_xz_y_y, t_0_xz_y_z, t_0_xz_z_x, t_0_xz_z_y, t_0_xz_z_z,\
                           t_0_yy_0_x, t_0_yy_0_xx, t_0_yy_0_xy, t_0_yy_0_xz, t_0_yy_0_y,\
                           t_0_yy_0_yy, t_0_yy_0_yz, t_0_yy_0_z, t_0_yy_0_zz, t_0_yy_x_x,\
                           t_0_yy_x_y, t_0_yy_x_z, t_0_yy_y_x, t_0_yy_y_y, t_0_yy_y_z,\
                           t_0_yy_z_x, t_0_yy_z_y, t_0_yy_z_z, t_0_yz_0_x, t_0_yz_0_xx,\
                           t_0_yz_0_xy, t_0_yz_0_xz, t_0_yz_0_y, t_0_yz_0_yy, t_0_yz_0_yz,\
                           t_0_yz_0_z, t_0_yz_0_zz, t_0_yz_x_x, t_0_yz_x_y, t_0_yz_x_z,\
                           t_0_yz_y_x, t_0_yz_y_y, t_0_yz_y_z, t_0_yz_z_x, t_0_yz_z_y,\
                           t_0_yz_z_z, t_0_zz_0_x, t_0_zz_0_xx, t_0_zz_0_xy, t_0_zz_0_xz,\
                           t_0_zz_0_y, t_0_zz_0_yy, t_0_zz_0_yz, t_0_zz_0_z, t_0_zz_0_zz,\
                           t_0_zz_x_x, t_0_zz_x_y, t_0_zz_x_z, t_0_zz_y_x, t_0_zz_y_y,\
                           t_0_zz_y_z, t_0_zz_z_x, t_0_zz_z_y, t_0_zz_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zz_z_z[i] = t_0_zz_0_zz[i] - rcd_z[i] * t_0_zz_0_z[i];

        t_0_zz_z_y[i] = t_0_zz_0_yz[i] - rcd_z[i] * t_0_zz_0_y[i];

        t_0_zz_z_x[i] = t_0_zz_0_xz[i] - rcd_z[i] * t_0_zz_0_x[i];

        t_0_zz_y_z[i] = t_0_zz_0_yz[i] - rcd_y[i] * t_0_zz_0_z[i];

        t_0_zz_y_y[i] = t_0_zz_0_yy[i] - rcd_y[i] * t_0_zz_0_y[i];

        t_0_zz_y_x[i] = t_0_zz_0_xy[i] - rcd_y[i] * t_0_zz_0_x[i];

        t_0_zz_x_z[i] = t_0_zz_0_xz[i] - rcd_x[i] * t_0_zz_0_z[i];

        t_0_zz_x_y[i] = t_0_zz_0_xy[i] - rcd_x[i] * t_0_zz_0_y[i];

        t_0_zz_x_x[i] = t_0_zz_0_xx[i] - rcd_x[i] * t_0_zz_0_x[i];

        t_0_yz_z_z[i] = t_0_yz_0_zz[i] - rcd_z[i] * t_0_yz_0_z[i];

        t_0_yz_z_y[i] = t_0_yz_0_yz[i] - rcd_z[i] * t_0_yz_0_y[i];

        t_0_yz_z_x[i] = t_0_yz_0_xz[i] - rcd_z[i] * t_0_yz_0_x[i];

        t_0_yz_y_z[i] = t_0_yz_0_yz[i] - rcd_y[i] * t_0_yz_0_z[i];

        t_0_yz_y_y[i] = t_0_yz_0_yy[i] - rcd_y[i] * t_0_yz_0_y[i];

        t_0_yz_y_x[i] = t_0_yz_0_xy[i] - rcd_y[i] * t_0_yz_0_x[i];

        t_0_yz_x_z[i] = t_0_yz_0_xz[i] - rcd_x[i] * t_0_yz_0_z[i];

        t_0_yz_x_y[i] = t_0_yz_0_xy[i] - rcd_x[i] * t_0_yz_0_y[i];

        t_0_yz_x_x[i] = t_0_yz_0_xx[i] - rcd_x[i] * t_0_yz_0_x[i];

        t_0_yy_z_z[i] = t_0_yy_0_zz[i] - rcd_z[i] * t_0_yy_0_z[i];

        t_0_yy_z_y[i] = t_0_yy_0_yz[i] - rcd_z[i] * t_0_yy_0_y[i];

        t_0_yy_z_x[i] = t_0_yy_0_xz[i] - rcd_z[i] * t_0_yy_0_x[i];

        t_0_yy_y_z[i] = t_0_yy_0_yz[i] - rcd_y[i] * t_0_yy_0_z[i];

        t_0_yy_y_y[i] = t_0_yy_0_yy[i] - rcd_y[i] * t_0_yy_0_y[i];

        t_0_yy_y_x[i] = t_0_yy_0_xy[i] - rcd_y[i] * t_0_yy_0_x[i];

        t_0_yy_x_z[i] = t_0_yy_0_xz[i] - rcd_x[i] * t_0_yy_0_z[i];

        t_0_yy_x_y[i] = t_0_yy_0_xy[i] - rcd_x[i] * t_0_yy_0_y[i];

        t_0_yy_x_x[i] = t_0_yy_0_xx[i] - rcd_x[i] * t_0_yy_0_x[i];

        t_0_xz_z_z[i] = t_0_xz_0_zz[i] - rcd_z[i] * t_0_xz_0_z[i];

        t_0_xz_z_y[i] = t_0_xz_0_yz[i] - rcd_z[i] * t_0_xz_0_y[i];

        t_0_xz_z_x[i] = t_0_xz_0_xz[i] - rcd_z[i] * t_0_xz_0_x[i];

        t_0_xz_y_z[i] = t_0_xz_0_yz[i] - rcd_y[i] * t_0_xz_0_z[i];

        t_0_xz_y_y[i] = t_0_xz_0_yy[i] - rcd_y[i] * t_0_xz_0_y[i];

        t_0_xz_y_x[i] = t_0_xz_0_xy[i] - rcd_y[i] * t_0_xz_0_x[i];

        t_0_xz_x_z[i] = t_0_xz_0_xz[i] - rcd_x[i] * t_0_xz_0_z[i];

        t_0_xz_x_y[i] = t_0_xz_0_xy[i] - rcd_x[i] * t_0_xz_0_y[i];

        t_0_xz_x_x[i] = t_0_xz_0_xx[i] - rcd_x[i] * t_0_xz_0_x[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xx_0_x, t_0_xx_0_xx, t_0_xx_0_xy,\
                           t_0_xx_0_xz, t_0_xx_0_y, t_0_xx_0_yy, t_0_xx_0_yz, t_0_xx_0_z,\
                           t_0_xx_0_zz, t_0_xx_x_x, t_0_xx_x_y, t_0_xx_x_z, t_0_xx_y_x,\
                           t_0_xx_y_y, t_0_xx_y_z, t_0_xx_z_x, t_0_xx_z_y, t_0_xx_z_z,\
                           t_0_xy_0_x, t_0_xy_0_xx, t_0_xy_0_xy, t_0_xy_0_xz, t_0_xy_0_y,\
                           t_0_xy_0_yy, t_0_xy_0_yz, t_0_xy_0_z, t_0_xy_0_zz, t_0_xy_x_x,\
                           t_0_xy_x_y, t_0_xy_x_z, t_0_xy_y_x, t_0_xy_y_y, t_0_xy_y_z,\
                           t_0_xy_z_x, t_0_xy_z_y, t_0_xy_z_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_xy_z_z[i] = t_0_xy_0_zz[i] - rcd_z[i] * t_0_xy_0_z[i];

        t_0_xy_z_y[i] = t_0_xy_0_yz[i] - rcd_z[i] * t_0_xy_0_y[i];

        t_0_xy_z_x[i] = t_0_xy_0_xz[i] - rcd_z[i] * t_0_xy_0_x[i];

        t_0_xy_y_z[i] = t_0_xy_0_yz[i] - rcd_y[i] * t_0_xy_0_z[i];

        t_0_xy_y_y[i] = t_0_xy_0_yy[i] - rcd_y[i] * t_0_xy_0_y[i];

        t_0_xy_y_x[i] = t_0_xy_0_xy[i] - rcd_y[i] * t_0_xy_0_x[i];

        t_0_xy_x_z[i] = t_0_xy_0_xz[i] - rcd_x[i] * t_0_xy_0_z[i];

        t_0_xy_x_y[i] = t_0_xy_0_xy[i] - rcd_x[i] * t_0_xy_0_y[i];

        t_0_xy_x_x[i] = t_0_xy_0_xx[i] - rcd_x[i] * t_0_xy_0_x[i];

        t_0_xx_z_z[i] = t_0_xx_0_zz[i] - rcd_z[i] * t_0_xx_0_z[i];

        t_0_xx_z_y[i] = t_0_xx_0_yz[i] - rcd_z[i] * t_0_xx_0_y[i];

        t_0_xx_z_x[i] = t_0_xx_0_xz[i] - rcd_z[i] * t_0_xx_0_x[i];

        t_0_xx_y_z[i] = t_0_xx_0_yz[i] - rcd_y[i] * t_0_xx_0_z[i];

        t_0_xx_y_y[i] = t_0_xx_0_yy[i] - rcd_y[i] * t_0_xx_0_y[i];

        t_0_xx_y_x[i] = t_0_xx_0_xy[i] - rcd_y[i] * t_0_xx_0_x[i];

        t_0_xx_x_z[i] = t_0_xx_0_xz[i] - rcd_x[i] * t_0_xx_0_z[i];

        t_0_xx_x_y[i] = t_0_xx_0_xy[i] - rcd_x[i] * t_0_xx_0_y[i];

        t_0_xx_x_x[i] = t_0_xx_0_xx[i] - rcd_x[i] * t_0_xx_0_x[i];
    }
}


} // derirec namespace
