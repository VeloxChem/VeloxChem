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
compHostHRRForPPSD_V0(      BufferHostXY<T>&      intsBufferPPSD,
                      const BufferHostX<int32_t>& intsIndexesPPSD,
                      const BufferHostXY<T>&      intsBufferSPSD,
                      const BufferHostX<int32_t>& intsIndexesSPSD,
                      const BufferHostXY<T>&      intsBufferSDSD,
                      const BufferHostX<int32_t>& intsIndexesSDSD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PPSD) integral components

    t_z_z_0_zz = intsBufferPPSD.data(intsIndexesPPSD(0));

    t_z_z_0_yz = intsBufferPPSD.data(intsIndexesPPSD(1));

    t_z_z_0_yy = intsBufferPPSD.data(intsIndexesPPSD(2));

    t_z_z_0_xz = intsBufferPPSD.data(intsIndexesPPSD(3));

    t_z_z_0_xy = intsBufferPPSD.data(intsIndexesPPSD(4));

    t_z_z_0_xx = intsBufferPPSD.data(intsIndexesPPSD(5));

    t_z_y_0_zz = intsBufferPPSD.data(intsIndexesPPSD(6));

    t_z_y_0_yz = intsBufferPPSD.data(intsIndexesPPSD(7));

    t_z_y_0_yy = intsBufferPPSD.data(intsIndexesPPSD(8));

    t_z_y_0_xz = intsBufferPPSD.data(intsIndexesPPSD(9));

    t_z_y_0_xy = intsBufferPPSD.data(intsIndexesPPSD(10));

    t_z_y_0_xx = intsBufferPPSD.data(intsIndexesPPSD(11));

    t_z_x_0_zz = intsBufferPPSD.data(intsIndexesPPSD(12));

    t_z_x_0_yz = intsBufferPPSD.data(intsIndexesPPSD(13));

    t_z_x_0_yy = intsBufferPPSD.data(intsIndexesPPSD(14));

    t_z_x_0_xz = intsBufferPPSD.data(intsIndexesPPSD(15));

    t_z_x_0_xy = intsBufferPPSD.data(intsIndexesPPSD(16));

    t_z_x_0_xx = intsBufferPPSD.data(intsIndexesPPSD(17));

    t_y_z_0_zz = intsBufferPPSD.data(intsIndexesPPSD(18));

    t_y_z_0_yz = intsBufferPPSD.data(intsIndexesPPSD(19));

    t_y_z_0_yy = intsBufferPPSD.data(intsIndexesPPSD(20));

    t_y_z_0_xz = intsBufferPPSD.data(intsIndexesPPSD(21));

    t_y_z_0_xy = intsBufferPPSD.data(intsIndexesPPSD(22));

    t_y_z_0_xx = intsBufferPPSD.data(intsIndexesPPSD(23));

    t_y_y_0_zz = intsBufferPPSD.data(intsIndexesPPSD(24));

    t_y_y_0_yz = intsBufferPPSD.data(intsIndexesPPSD(25));

    t_y_y_0_yy = intsBufferPPSD.data(intsIndexesPPSD(26));

    t_y_y_0_xz = intsBufferPPSD.data(intsIndexesPPSD(27));

    t_y_y_0_xy = intsBufferPPSD.data(intsIndexesPPSD(28));

    t_y_y_0_xx = intsBufferPPSD.data(intsIndexesPPSD(29));

    t_y_x_0_zz = intsBufferPPSD.data(intsIndexesPPSD(30));

    t_y_x_0_yz = intsBufferPPSD.data(intsIndexesPPSD(31));

    t_y_x_0_yy = intsBufferPPSD.data(intsIndexesPPSD(32));

    t_y_x_0_xz = intsBufferPPSD.data(intsIndexesPPSD(33));

    t_y_x_0_xy = intsBufferPPSD.data(intsIndexesPPSD(34));

    t_y_x_0_xx = intsBufferPPSD.data(intsIndexesPPSD(35));

    t_x_z_0_zz = intsBufferPPSD.data(intsIndexesPPSD(36));

    t_x_z_0_yz = intsBufferPPSD.data(intsIndexesPPSD(37));

    t_x_z_0_yy = intsBufferPPSD.data(intsIndexesPPSD(38));

    t_x_z_0_xz = intsBufferPPSD.data(intsIndexesPPSD(39));

    t_x_z_0_xy = intsBufferPPSD.data(intsIndexesPPSD(40));

    t_x_z_0_xx = intsBufferPPSD.data(intsIndexesPPSD(41));

    t_x_y_0_zz = intsBufferPPSD.data(intsIndexesPPSD(42));

    t_x_y_0_yz = intsBufferPPSD.data(intsIndexesPPSD(43));

    t_x_y_0_yy = intsBufferPPSD.data(intsIndexesPPSD(44));

    t_x_y_0_xz = intsBufferPPSD.data(intsIndexesPPSD(45));

    t_x_y_0_xy = intsBufferPPSD.data(intsIndexesPPSD(46));

    t_x_y_0_xx = intsBufferPPSD.data(intsIndexesPPSD(47));

    t_x_x_0_zz = intsBufferPPSD.data(intsIndexesPPSD(48));

    t_x_x_0_yz = intsBufferPPSD.data(intsIndexesPPSD(49));

    t_x_x_0_yy = intsBufferPPSD.data(intsIndexesPPSD(50));

    t_x_x_0_xz = intsBufferPPSD.data(intsIndexesPPSD(51));

    t_x_x_0_xy = intsBufferPPSD.data(intsIndexesPPSD(52));

    t_x_x_0_xx = intsBufferPPSD.data(intsIndexesPPSD(53));

    // set up (SPSD) integral components

    t_0_z_0_zz = intsBufferSPSD.data(intsIndexesSPSD(0));

    t_0_z_0_yz = intsBufferSPSD.data(intsIndexesSPSD(1));

    t_0_z_0_yy = intsBufferSPSD.data(intsIndexesSPSD(2));

    t_0_z_0_xz = intsBufferSPSD.data(intsIndexesSPSD(3));

    t_0_z_0_xy = intsBufferSPSD.data(intsIndexesSPSD(4));

    t_0_z_0_xx = intsBufferSPSD.data(intsIndexesSPSD(5));

    t_0_y_0_zz = intsBufferSPSD.data(intsIndexesSPSD(6));

    t_0_y_0_yz = intsBufferSPSD.data(intsIndexesSPSD(7));

    t_0_y_0_yy = intsBufferSPSD.data(intsIndexesSPSD(8));

    t_0_y_0_xz = intsBufferSPSD.data(intsIndexesSPSD(9));

    t_0_y_0_xy = intsBufferSPSD.data(intsIndexesSPSD(10));

    t_0_y_0_xx = intsBufferSPSD.data(intsIndexesSPSD(11));

    t_0_x_0_zz = intsBufferSPSD.data(intsIndexesSPSD(12));

    t_0_x_0_yz = intsBufferSPSD.data(intsIndexesSPSD(13));

    t_0_x_0_yy = intsBufferSPSD.data(intsIndexesSPSD(14));

    t_0_x_0_xz = intsBufferSPSD.data(intsIndexesSPSD(15));

    t_0_x_0_xy = intsBufferSPSD.data(intsIndexesSPSD(16));

    t_0_x_0_xx = intsBufferSPSD.data(intsIndexesSPSD(17));

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

    #pragma omp simd align(rab_y, rab_z, t_0_x_0_xx, t_0_x_0_xy, t_0_x_0_xz, t_0_x_0_yy,\
                           t_0_x_0_yz, t_0_x_0_zz, t_0_xy_0_xx, t_0_xy_0_xy, t_0_xy_0_xz,\
                           t_0_xy_0_yy, t_0_xy_0_yz, t_0_xy_0_zz, t_0_xz_0_xx, t_0_xz_0_xy,\
                           t_0_xz_0_xz, t_0_xz_0_yy, t_0_xz_0_yz, t_0_xz_0_zz, t_0_y_0_xx,\
                           t_0_y_0_xy, t_0_y_0_xz, t_0_y_0_yy, t_0_y_0_yz, t_0_y_0_zz,\
                           t_0_yy_0_xx, t_0_yy_0_xy, t_0_yy_0_xz, t_0_yy_0_yy, t_0_yy_0_yz,\
                           t_0_yy_0_zz, t_0_yz_0_xx, t_0_yz_0_xy, t_0_yz_0_xz, t_0_yz_0_yy,\
                           t_0_yz_0_yz, t_0_yz_0_zz, t_0_z_0_xx, t_0_z_0_xy, t_0_z_0_xz,\
                           t_0_z_0_yy, t_0_z_0_yz, t_0_z_0_zz, t_0_zz_0_xx, t_0_zz_0_xy,\
                           t_0_zz_0_xz, t_0_zz_0_yy, t_0_zz_0_yz, t_0_zz_0_zz, t_y_x_0_xx,\
                           t_y_x_0_xy, t_y_x_0_xz, t_y_x_0_yy, t_y_x_0_yz, t_y_x_0_zz,\
                           t_y_y_0_xx, t_y_y_0_xy, t_y_y_0_xz, t_y_y_0_yy, t_y_y_0_yz,\
                           t_y_y_0_zz, t_y_z_0_xx, t_y_z_0_xy, t_y_z_0_xz, t_y_z_0_yy,\
                           t_y_z_0_yz, t_y_z_0_zz, t_z_x_0_xx, t_z_x_0_xy, t_z_x_0_xz,\
                           t_z_x_0_yy, t_z_x_0_yz, t_z_x_0_zz, t_z_y_0_xx, t_z_y_0_xy,\
                           t_z_y_0_xz, t_z_y_0_yy, t_z_y_0_yz, t_z_y_0_zz, t_z_z_0_xx,\
                           t_z_z_0_xy, t_z_z_0_xz, t_z_z_0_yy, t_z_z_0_yz, t_z_z_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_z_0_zz[i] = t_0_zz_0_zz[i] - rab_z[i] * t_0_z_0_zz[i];

        t_z_z_0_yz[i] = t_0_zz_0_yz[i] - rab_z[i] * t_0_z_0_yz[i];

        t_z_z_0_yy[i] = t_0_zz_0_yy[i] - rab_z[i] * t_0_z_0_yy[i];

        t_z_z_0_xz[i] = t_0_zz_0_xz[i] - rab_z[i] * t_0_z_0_xz[i];

        t_z_z_0_xy[i] = t_0_zz_0_xy[i] - rab_z[i] * t_0_z_0_xy[i];

        t_z_z_0_xx[i] = t_0_zz_0_xx[i] - rab_z[i] * t_0_z_0_xx[i];

        t_z_y_0_zz[i] = t_0_yz_0_zz[i] - rab_z[i] * t_0_y_0_zz[i];

        t_z_y_0_yz[i] = t_0_yz_0_yz[i] - rab_z[i] * t_0_y_0_yz[i];

        t_z_y_0_yy[i] = t_0_yz_0_yy[i] - rab_z[i] * t_0_y_0_yy[i];

        t_z_y_0_xz[i] = t_0_yz_0_xz[i] - rab_z[i] * t_0_y_0_xz[i];

        t_z_y_0_xy[i] = t_0_yz_0_xy[i] - rab_z[i] * t_0_y_0_xy[i];

        t_z_y_0_xx[i] = t_0_yz_0_xx[i] - rab_z[i] * t_0_y_0_xx[i];

        t_z_x_0_zz[i] = t_0_xz_0_zz[i] - rab_z[i] * t_0_x_0_zz[i];

        t_z_x_0_yz[i] = t_0_xz_0_yz[i] - rab_z[i] * t_0_x_0_yz[i];

        t_z_x_0_yy[i] = t_0_xz_0_yy[i] - rab_z[i] * t_0_x_0_yy[i];

        t_z_x_0_xz[i] = t_0_xz_0_xz[i] - rab_z[i] * t_0_x_0_xz[i];

        t_z_x_0_xy[i] = t_0_xz_0_xy[i] - rab_z[i] * t_0_x_0_xy[i];

        t_z_x_0_xx[i] = t_0_xz_0_xx[i] - rab_z[i] * t_0_x_0_xx[i];

        t_y_z_0_zz[i] = t_0_yz_0_zz[i] - rab_y[i] * t_0_z_0_zz[i];

        t_y_z_0_yz[i] = t_0_yz_0_yz[i] - rab_y[i] * t_0_z_0_yz[i];

        t_y_z_0_yy[i] = t_0_yz_0_yy[i] - rab_y[i] * t_0_z_0_yy[i];

        t_y_z_0_xz[i] = t_0_yz_0_xz[i] - rab_y[i] * t_0_z_0_xz[i];

        t_y_z_0_xy[i] = t_0_yz_0_xy[i] - rab_y[i] * t_0_z_0_xy[i];

        t_y_z_0_xx[i] = t_0_yz_0_xx[i] - rab_y[i] * t_0_z_0_xx[i];

        t_y_y_0_zz[i] = t_0_yy_0_zz[i] - rab_y[i] * t_0_y_0_zz[i];

        t_y_y_0_yz[i] = t_0_yy_0_yz[i] - rab_y[i] * t_0_y_0_yz[i];

        t_y_y_0_yy[i] = t_0_yy_0_yy[i] - rab_y[i] * t_0_y_0_yy[i];

        t_y_y_0_xz[i] = t_0_yy_0_xz[i] - rab_y[i] * t_0_y_0_xz[i];

        t_y_y_0_xy[i] = t_0_yy_0_xy[i] - rab_y[i] * t_0_y_0_xy[i];

        t_y_y_0_xx[i] = t_0_yy_0_xx[i] - rab_y[i] * t_0_y_0_xx[i];

        t_y_x_0_zz[i] = t_0_xy_0_zz[i] - rab_y[i] * t_0_x_0_zz[i];

        t_y_x_0_yz[i] = t_0_xy_0_yz[i] - rab_y[i] * t_0_x_0_yz[i];

        t_y_x_0_yy[i] = t_0_xy_0_yy[i] - rab_y[i] * t_0_x_0_yy[i];

        t_y_x_0_xz[i] = t_0_xy_0_xz[i] - rab_y[i] * t_0_x_0_xz[i];

        t_y_x_0_xy[i] = t_0_xy_0_xy[i] - rab_y[i] * t_0_x_0_xy[i];

        t_y_x_0_xx[i] = t_0_xy_0_xx[i] - rab_y[i] * t_0_x_0_xx[i];
    }

    #pragma omp simd align(rab_x, t_0_x_0_xx, t_0_x_0_xy, t_0_x_0_xz, t_0_x_0_yy, t_0_x_0_yz,\
                           t_0_x_0_zz, t_0_xx_0_xx, t_0_xx_0_xy, t_0_xx_0_xz, t_0_xx_0_yy,\
                           t_0_xx_0_yz, t_0_xx_0_zz, t_0_xy_0_xx, t_0_xy_0_xy, t_0_xy_0_xz,\
                           t_0_xy_0_yy, t_0_xy_0_yz, t_0_xy_0_zz, t_0_xz_0_xx, t_0_xz_0_xy,\
                           t_0_xz_0_xz, t_0_xz_0_yy, t_0_xz_0_yz, t_0_xz_0_zz, t_0_y_0_xx,\
                           t_0_y_0_xy, t_0_y_0_xz, t_0_y_0_yy, t_0_y_0_yz, t_0_y_0_zz,\
                           t_0_z_0_xx, t_0_z_0_xy, t_0_z_0_xz, t_0_z_0_yy, t_0_z_0_yz,\
                           t_0_z_0_zz, t_x_x_0_xx, t_x_x_0_xy, t_x_x_0_xz, t_x_x_0_yy,\
                           t_x_x_0_yz, t_x_x_0_zz, t_x_y_0_xx, t_x_y_0_xy, t_x_y_0_xz,\
                           t_x_y_0_yy, t_x_y_0_yz, t_x_y_0_zz, t_x_z_0_xx, t_x_z_0_xy,\
                           t_x_z_0_xz, t_x_z_0_yy, t_x_z_0_yz, t_x_z_0_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_z_0_zz[i] = t_0_xz_0_zz[i] - rab_x[i] * t_0_z_0_zz[i];

        t_x_z_0_yz[i] = t_0_xz_0_yz[i] - rab_x[i] * t_0_z_0_yz[i];

        t_x_z_0_yy[i] = t_0_xz_0_yy[i] - rab_x[i] * t_0_z_0_yy[i];

        t_x_z_0_xz[i] = t_0_xz_0_xz[i] - rab_x[i] * t_0_z_0_xz[i];

        t_x_z_0_xy[i] = t_0_xz_0_xy[i] - rab_x[i] * t_0_z_0_xy[i];

        t_x_z_0_xx[i] = t_0_xz_0_xx[i] - rab_x[i] * t_0_z_0_xx[i];

        t_x_y_0_zz[i] = t_0_xy_0_zz[i] - rab_x[i] * t_0_y_0_zz[i];

        t_x_y_0_yz[i] = t_0_xy_0_yz[i] - rab_x[i] * t_0_y_0_yz[i];

        t_x_y_0_yy[i] = t_0_xy_0_yy[i] - rab_x[i] * t_0_y_0_yy[i];

        t_x_y_0_xz[i] = t_0_xy_0_xz[i] - rab_x[i] * t_0_y_0_xz[i];

        t_x_y_0_xy[i] = t_0_xy_0_xy[i] - rab_x[i] * t_0_y_0_xy[i];

        t_x_y_0_xx[i] = t_0_xy_0_xx[i] - rab_x[i] * t_0_y_0_xx[i];

        t_x_x_0_zz[i] = t_0_xx_0_zz[i] - rab_x[i] * t_0_x_0_zz[i];

        t_x_x_0_yz[i] = t_0_xx_0_yz[i] - rab_x[i] * t_0_x_0_yz[i];

        t_x_x_0_yy[i] = t_0_xx_0_yy[i] - rab_x[i] * t_0_x_0_yy[i];

        t_x_x_0_xz[i] = t_0_xx_0_xz[i] - rab_x[i] * t_0_x_0_xz[i];

        t_x_x_0_xy[i] = t_0_xx_0_xy[i] - rab_x[i] * t_0_x_0_xy[i];

        t_x_x_0_xx[i] = t_0_xx_0_xx[i] - rab_x[i] * t_0_x_0_xx[i];
    }
}


} // derirec namespace
