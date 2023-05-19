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
compHostHRRForSPPD_V0(      BufferHostXY<T>&      intsBufferSPPD,
                      const BufferHostX<int32_t>& intsIndexesSPPD,
                      const BufferHostXY<T>&      intsBufferSPSD,
                      const BufferHostX<int32_t>& intsIndexesSPSD,
                      const BufferHostXY<T>&      intsBufferSPSF,
                      const BufferHostX<int32_t>& intsIndexesSPSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SPPD) integral components

    t_0_z_z_zz = intsBufferSPPD.data(intsIndexesSPPD(0));

    t_0_z_z_yz = intsBufferSPPD.data(intsIndexesSPPD(1));

    t_0_z_z_yy = intsBufferSPPD.data(intsIndexesSPPD(2));

    t_0_z_z_xz = intsBufferSPPD.data(intsIndexesSPPD(3));

    t_0_z_z_xy = intsBufferSPPD.data(intsIndexesSPPD(4));

    t_0_z_z_xx = intsBufferSPPD.data(intsIndexesSPPD(5));

    t_0_z_y_zz = intsBufferSPPD.data(intsIndexesSPPD(6));

    t_0_z_y_yz = intsBufferSPPD.data(intsIndexesSPPD(7));

    t_0_z_y_yy = intsBufferSPPD.data(intsIndexesSPPD(8));

    t_0_z_y_xz = intsBufferSPPD.data(intsIndexesSPPD(9));

    t_0_z_y_xy = intsBufferSPPD.data(intsIndexesSPPD(10));

    t_0_z_y_xx = intsBufferSPPD.data(intsIndexesSPPD(11));

    t_0_z_x_zz = intsBufferSPPD.data(intsIndexesSPPD(12));

    t_0_z_x_yz = intsBufferSPPD.data(intsIndexesSPPD(13));

    t_0_z_x_yy = intsBufferSPPD.data(intsIndexesSPPD(14));

    t_0_z_x_xz = intsBufferSPPD.data(intsIndexesSPPD(15));

    t_0_z_x_xy = intsBufferSPPD.data(intsIndexesSPPD(16));

    t_0_z_x_xx = intsBufferSPPD.data(intsIndexesSPPD(17));

    t_0_y_z_zz = intsBufferSPPD.data(intsIndexesSPPD(18));

    t_0_y_z_yz = intsBufferSPPD.data(intsIndexesSPPD(19));

    t_0_y_z_yy = intsBufferSPPD.data(intsIndexesSPPD(20));

    t_0_y_z_xz = intsBufferSPPD.data(intsIndexesSPPD(21));

    t_0_y_z_xy = intsBufferSPPD.data(intsIndexesSPPD(22));

    t_0_y_z_xx = intsBufferSPPD.data(intsIndexesSPPD(23));

    t_0_y_y_zz = intsBufferSPPD.data(intsIndexesSPPD(24));

    t_0_y_y_yz = intsBufferSPPD.data(intsIndexesSPPD(25));

    t_0_y_y_yy = intsBufferSPPD.data(intsIndexesSPPD(26));

    t_0_y_y_xz = intsBufferSPPD.data(intsIndexesSPPD(27));

    t_0_y_y_xy = intsBufferSPPD.data(intsIndexesSPPD(28));

    t_0_y_y_xx = intsBufferSPPD.data(intsIndexesSPPD(29));

    t_0_y_x_zz = intsBufferSPPD.data(intsIndexesSPPD(30));

    t_0_y_x_yz = intsBufferSPPD.data(intsIndexesSPPD(31));

    t_0_y_x_yy = intsBufferSPPD.data(intsIndexesSPPD(32));

    t_0_y_x_xz = intsBufferSPPD.data(intsIndexesSPPD(33));

    t_0_y_x_xy = intsBufferSPPD.data(intsIndexesSPPD(34));

    t_0_y_x_xx = intsBufferSPPD.data(intsIndexesSPPD(35));

    t_0_x_z_zz = intsBufferSPPD.data(intsIndexesSPPD(36));

    t_0_x_z_yz = intsBufferSPPD.data(intsIndexesSPPD(37));

    t_0_x_z_yy = intsBufferSPPD.data(intsIndexesSPPD(38));

    t_0_x_z_xz = intsBufferSPPD.data(intsIndexesSPPD(39));

    t_0_x_z_xy = intsBufferSPPD.data(intsIndexesSPPD(40));

    t_0_x_z_xx = intsBufferSPPD.data(intsIndexesSPPD(41));

    t_0_x_y_zz = intsBufferSPPD.data(intsIndexesSPPD(42));

    t_0_x_y_yz = intsBufferSPPD.data(intsIndexesSPPD(43));

    t_0_x_y_yy = intsBufferSPPD.data(intsIndexesSPPD(44));

    t_0_x_y_xz = intsBufferSPPD.data(intsIndexesSPPD(45));

    t_0_x_y_xy = intsBufferSPPD.data(intsIndexesSPPD(46));

    t_0_x_y_xx = intsBufferSPPD.data(intsIndexesSPPD(47));

    t_0_x_x_zz = intsBufferSPPD.data(intsIndexesSPPD(48));

    t_0_x_x_yz = intsBufferSPPD.data(intsIndexesSPPD(49));

    t_0_x_x_yy = intsBufferSPPD.data(intsIndexesSPPD(50));

    t_0_x_x_xz = intsBufferSPPD.data(intsIndexesSPPD(51));

    t_0_x_x_xy = intsBufferSPPD.data(intsIndexesSPPD(52));

    t_0_x_x_xx = intsBufferSPPD.data(intsIndexesSPPD(53));

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

    // set up (SPSF) integral components

    t_0_z_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(0));

    t_0_z_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(1));

    t_0_z_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(2));

    t_0_z_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(3));

    t_0_z_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(4));

    t_0_z_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(5));

    t_0_z_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(6));

    t_0_z_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(7));

    t_0_z_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(8));

    t_0_z_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(9));

    t_0_y_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(10));

    t_0_y_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(11));

    t_0_y_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(12));

    t_0_y_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(13));

    t_0_y_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(14));

    t_0_y_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(15));

    t_0_y_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(16));

    t_0_y_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(17));

    t_0_y_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(18));

    t_0_y_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(19));

    t_0_x_0_zzz = intsBufferSPSF.data(intsIndexesSPSF(20));

    t_0_x_0_yzz = intsBufferSPSF.data(intsIndexesSPSF(21));

    t_0_x_0_yyz = intsBufferSPSF.data(intsIndexesSPSF(22));

    t_0_x_0_yyy = intsBufferSPSF.data(intsIndexesSPSF(23));

    t_0_x_0_xzz = intsBufferSPSF.data(intsIndexesSPSF(24));

    t_0_x_0_xyz = intsBufferSPSF.data(intsIndexesSPSF(25));

    t_0_x_0_xyy = intsBufferSPSF.data(intsIndexesSPSF(26));

    t_0_x_0_xxz = intsBufferSPSF.data(intsIndexesSPSF(27));

    t_0_x_0_xxy = intsBufferSPSF.data(intsIndexesSPSF(28));

    t_0_x_0_xxx = intsBufferSPSF.data(intsIndexesSPSF(29));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_y_0_xx, t_0_y_0_xxx, t_0_y_0_xxy,\
                           t_0_y_0_xxz, t_0_y_0_xy, t_0_y_0_xyy, t_0_y_0_xyz, t_0_y_0_xz,\
                           t_0_y_0_xzz, t_0_y_0_yy, t_0_y_0_yyy, t_0_y_0_yyz, t_0_y_0_yz,\
                           t_0_y_0_yzz, t_0_y_0_zz, t_0_y_0_zzz, t_0_y_x_xx, t_0_y_x_xy,\
                           t_0_y_x_xz, t_0_y_x_yy, t_0_y_x_yz, t_0_y_x_zz, t_0_y_y_xx,\
                           t_0_y_y_xy, t_0_y_y_xz, t_0_y_y_yy, t_0_y_y_yz, t_0_y_y_zz,\
                           t_0_y_z_xx, t_0_y_z_xy, t_0_y_z_xz, t_0_y_z_yy, t_0_y_z_yz,\
                           t_0_y_z_zz, t_0_z_0_xx, t_0_z_0_xxx, t_0_z_0_xxy, t_0_z_0_xxz,\
                           t_0_z_0_xy, t_0_z_0_xyy, t_0_z_0_xyz, t_0_z_0_xz, t_0_z_0_xzz,\
                           t_0_z_0_yy, t_0_z_0_yyy, t_0_z_0_yyz, t_0_z_0_yz, t_0_z_0_yzz,\
                           t_0_z_0_zz, t_0_z_0_zzz, t_0_z_x_xx, t_0_z_x_xy, t_0_z_x_xz,\
                           t_0_z_x_yy, t_0_z_x_yz, t_0_z_x_zz, t_0_z_y_xx, t_0_z_y_xy,\
                           t_0_z_y_xz, t_0_z_y_yy, t_0_z_y_yz, t_0_z_y_zz, t_0_z_z_xx,\
                           t_0_z_z_xy, t_0_z_z_xz, t_0_z_z_yy, t_0_z_z_yz, t_0_z_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_z_z_zz[i] = t_0_z_0_zzz[i] - rcd_z[i] * t_0_z_0_zz[i];

        t_0_z_z_yz[i] = t_0_z_0_yzz[i] - rcd_z[i] * t_0_z_0_yz[i];

        t_0_z_z_yy[i] = t_0_z_0_yyz[i] - rcd_z[i] * t_0_z_0_yy[i];

        t_0_z_z_xz[i] = t_0_z_0_xzz[i] - rcd_z[i] * t_0_z_0_xz[i];

        t_0_z_z_xy[i] = t_0_z_0_xyz[i] - rcd_z[i] * t_0_z_0_xy[i];

        t_0_z_z_xx[i] = t_0_z_0_xxz[i] - rcd_z[i] * t_0_z_0_xx[i];

        t_0_z_y_zz[i] = t_0_z_0_yzz[i] - rcd_y[i] * t_0_z_0_zz[i];

        t_0_z_y_yz[i] = t_0_z_0_yyz[i] - rcd_y[i] * t_0_z_0_yz[i];

        t_0_z_y_yy[i] = t_0_z_0_yyy[i] - rcd_y[i] * t_0_z_0_yy[i];

        t_0_z_y_xz[i] = t_0_z_0_xyz[i] - rcd_y[i] * t_0_z_0_xz[i];

        t_0_z_y_xy[i] = t_0_z_0_xyy[i] - rcd_y[i] * t_0_z_0_xy[i];

        t_0_z_y_xx[i] = t_0_z_0_xxy[i] - rcd_y[i] * t_0_z_0_xx[i];

        t_0_z_x_zz[i] = t_0_z_0_xzz[i] - rcd_x[i] * t_0_z_0_zz[i];

        t_0_z_x_yz[i] = t_0_z_0_xyz[i] - rcd_x[i] * t_0_z_0_yz[i];

        t_0_z_x_yy[i] = t_0_z_0_xyy[i] - rcd_x[i] * t_0_z_0_yy[i];

        t_0_z_x_xz[i] = t_0_z_0_xxz[i] - rcd_x[i] * t_0_z_0_xz[i];

        t_0_z_x_xy[i] = t_0_z_0_xxy[i] - rcd_x[i] * t_0_z_0_xy[i];

        t_0_z_x_xx[i] = t_0_z_0_xxx[i] - rcd_x[i] * t_0_z_0_xx[i];

        t_0_y_z_zz[i] = t_0_y_0_zzz[i] - rcd_z[i] * t_0_y_0_zz[i];

        t_0_y_z_yz[i] = t_0_y_0_yzz[i] - rcd_z[i] * t_0_y_0_yz[i];

        t_0_y_z_yy[i] = t_0_y_0_yyz[i] - rcd_z[i] * t_0_y_0_yy[i];

        t_0_y_z_xz[i] = t_0_y_0_xzz[i] - rcd_z[i] * t_0_y_0_xz[i];

        t_0_y_z_xy[i] = t_0_y_0_xyz[i] - rcd_z[i] * t_0_y_0_xy[i];

        t_0_y_z_xx[i] = t_0_y_0_xxz[i] - rcd_z[i] * t_0_y_0_xx[i];

        t_0_y_y_zz[i] = t_0_y_0_yzz[i] - rcd_y[i] * t_0_y_0_zz[i];

        t_0_y_y_yz[i] = t_0_y_0_yyz[i] - rcd_y[i] * t_0_y_0_yz[i];

        t_0_y_y_yy[i] = t_0_y_0_yyy[i] - rcd_y[i] * t_0_y_0_yy[i];

        t_0_y_y_xz[i] = t_0_y_0_xyz[i] - rcd_y[i] * t_0_y_0_xz[i];

        t_0_y_y_xy[i] = t_0_y_0_xyy[i] - rcd_y[i] * t_0_y_0_xy[i];

        t_0_y_y_xx[i] = t_0_y_0_xxy[i] - rcd_y[i] * t_0_y_0_xx[i];

        t_0_y_x_zz[i] = t_0_y_0_xzz[i] - rcd_x[i] * t_0_y_0_zz[i];

        t_0_y_x_yz[i] = t_0_y_0_xyz[i] - rcd_x[i] * t_0_y_0_yz[i];

        t_0_y_x_yy[i] = t_0_y_0_xyy[i] - rcd_x[i] * t_0_y_0_yy[i];

        t_0_y_x_xz[i] = t_0_y_0_xxz[i] - rcd_x[i] * t_0_y_0_xz[i];

        t_0_y_x_xy[i] = t_0_y_0_xxy[i] - rcd_x[i] * t_0_y_0_xy[i];

        t_0_y_x_xx[i] = t_0_y_0_xxx[i] - rcd_x[i] * t_0_y_0_xx[i];
    }

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_x_0_xx, t_0_x_0_xxx, t_0_x_0_xxy,\
                           t_0_x_0_xxz, t_0_x_0_xy, t_0_x_0_xyy, t_0_x_0_xyz, t_0_x_0_xz,\
                           t_0_x_0_xzz, t_0_x_0_yy, t_0_x_0_yyy, t_0_x_0_yyz, t_0_x_0_yz,\
                           t_0_x_0_yzz, t_0_x_0_zz, t_0_x_0_zzz, t_0_x_x_xx, t_0_x_x_xy,\
                           t_0_x_x_xz, t_0_x_x_yy, t_0_x_x_yz, t_0_x_x_zz, t_0_x_y_xx,\
                           t_0_x_y_xy, t_0_x_y_xz, t_0_x_y_yy, t_0_x_y_yz, t_0_x_y_zz,\
                           t_0_x_z_xx, t_0_x_z_xy, t_0_x_z_xz, t_0_x_z_yy, t_0_x_z_yz,\
                           t_0_x_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_x_z_zz[i] = t_0_x_0_zzz[i] - rcd_z[i] * t_0_x_0_zz[i];

        t_0_x_z_yz[i] = t_0_x_0_yzz[i] - rcd_z[i] * t_0_x_0_yz[i];

        t_0_x_z_yy[i] = t_0_x_0_yyz[i] - rcd_z[i] * t_0_x_0_yy[i];

        t_0_x_z_xz[i] = t_0_x_0_xzz[i] - rcd_z[i] * t_0_x_0_xz[i];

        t_0_x_z_xy[i] = t_0_x_0_xyz[i] - rcd_z[i] * t_0_x_0_xy[i];

        t_0_x_z_xx[i] = t_0_x_0_xxz[i] - rcd_z[i] * t_0_x_0_xx[i];

        t_0_x_y_zz[i] = t_0_x_0_yzz[i] - rcd_y[i] * t_0_x_0_zz[i];

        t_0_x_y_yz[i] = t_0_x_0_yyz[i] - rcd_y[i] * t_0_x_0_yz[i];

        t_0_x_y_yy[i] = t_0_x_0_yyy[i] - rcd_y[i] * t_0_x_0_yy[i];

        t_0_x_y_xz[i] = t_0_x_0_xyz[i] - rcd_y[i] * t_0_x_0_xz[i];

        t_0_x_y_xy[i] = t_0_x_0_xyy[i] - rcd_y[i] * t_0_x_0_xy[i];

        t_0_x_y_xx[i] = t_0_x_0_xxy[i] - rcd_y[i] * t_0_x_0_xx[i];

        t_0_x_x_zz[i] = t_0_x_0_xzz[i] - rcd_x[i] * t_0_x_0_zz[i];

        t_0_x_x_yz[i] = t_0_x_0_xyz[i] - rcd_x[i] * t_0_x_0_yz[i];

        t_0_x_x_yy[i] = t_0_x_0_xyy[i] - rcd_x[i] * t_0_x_0_yy[i];

        t_0_x_x_xz[i] = t_0_x_0_xxz[i] - rcd_x[i] * t_0_x_0_xz[i];

        t_0_x_x_xy[i] = t_0_x_0_xxy[i] - rcd_x[i] * t_0_x_0_xy[i];

        t_0_x_x_xx[i] = t_0_x_0_xxx[i] - rcd_x[i] * t_0_x_0_xx[i];
    }
}


} // derirec namespace
