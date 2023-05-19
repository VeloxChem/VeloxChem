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
compHostHRRForSFPD_V0(      BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPD) integral components

    t_0_zzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(0));

    t_0_zzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(1));

    t_0_zzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(2));

    t_0_yzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(3));

    t_0_yzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(4));

    t_0_yzz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(5));

    t_0_yzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(6));

    t_0_yzz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(7));

    t_0_yyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(8));

    t_0_yyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(9));

    t_0_yyz_y_yy = intsBufferSFPD.data(intsIndexesSFPD(10));

    t_0_yyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(11));

    t_0_yyz_x_yy = intsBufferSFPD.data(intsIndexesSFPD(12));

    t_0_yyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(13));

    t_0_yyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(14));

    t_0_xzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(15));

    t_0_xzz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(16));

    t_0_xzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(17));

    t_0_xzz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(18));

    t_0_xyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(19));

    t_0_xyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(20));

    t_0_xyz_y_xy = intsBufferSFPD.data(intsIndexesSFPD(21));

    t_0_xyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(22));

    t_0_xyz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(23));

    t_0_xyz_x_xy = intsBufferSFPD.data(intsIndexesSFPD(24));

    t_0_xyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(25));

    t_0_xyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(26));

    t_0_xyy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(27));

    t_0_xxz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(28));

    t_0_xxz_y_xx = intsBufferSFPD.data(intsIndexesSFPD(29));

    t_0_xxz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(30));

    t_0_xxz_x_xx = intsBufferSFPD.data(intsIndexesSFPD(31));

    t_0_xxy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(32));

    t_0_xxy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(33));

    t_0_xxy_x_xx = intsBufferSFPD.data(intsIndexesSFPD(34));

    t_0_xxx_x_xx = intsBufferSFPD.data(intsIndexesSFPD(35));

    // set up (SFSD) integral components

    t_0_zzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(0));

    t_0_yzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(1));

    t_0_yzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(2));

    t_0_yyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(3));

    t_0_yyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(4));

    t_0_yyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(5));

    t_0_xzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(6));

    t_0_xzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(7));

    t_0_xyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(8));

    t_0_xyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(9));

    t_0_xyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(10));

    t_0_xyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(11));

    t_0_xyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(12));

    t_0_xxz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(13));

    t_0_xxz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(14));

    t_0_xxy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(15));

    t_0_xxy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(16));

    t_0_xxx_0_xx = intsBufferSFSD.data(intsIndexesSFSD(17));

    // set up (SFSF) integral components

    t_0_zzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(0));

    t_0_zzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(1));

    t_0_zzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(2));

    t_0_yzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(3));

    t_0_yzz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(4));

    t_0_yzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(5));

    t_0_yzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(6));

    t_0_yyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(7));

    t_0_yyz_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(8));

    t_0_yyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(9));

    t_0_yyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(10));

    t_0_yyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(11));

    t_0_yyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(12));

    t_0_xzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(13));

    t_0_xzz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(14));

    t_0_xzz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(15));

    t_0_xyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(16));

    t_0_xyz_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(17));

    t_0_xyz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(18));

    t_0_xyz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(19));

    t_0_xyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(20));

    t_0_xyy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(21));

    t_0_xxz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(22));

    t_0_xxz_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(23));

    t_0_xxz_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(24));

    t_0_xxy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(25));

    t_0_xxy_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(26));

    t_0_xxx_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(27));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_0_xx, t_0_xxx_0_xxx, t_0_xxx_x_xx,\
                           t_0_xxy_0_xx, t_0_xxy_0_xxx, t_0_xxy_0_xxy, t_0_xxy_0_xy,\
                           t_0_xxy_x_xx, t_0_xxy_x_xy, t_0_xxy_y_xx, t_0_xxz_0_xx, t_0_xxz_0_xxx,\
                           t_0_xxz_0_xxy, t_0_xxz_0_xxz, t_0_xxz_0_xz, t_0_xxz_x_xx,\
                           t_0_xxz_x_xz, t_0_xxz_y_xx, t_0_xxz_z_xx, t_0_xyy_0_xxy,\
                           t_0_xyy_0_xy, t_0_xyy_0_xyy, t_0_xyy_0_yy, t_0_xyy_x_xy,\
                           t_0_xyy_x_yy, t_0_xyy_y_xy, t_0_xyz_0_xxy, t_0_xyz_0_xxz,\
                           t_0_xyz_0_xy, t_0_xyz_0_xyy, t_0_xyz_0_xyz, t_0_xyz_0_xz,\
                           t_0_xyz_0_yz, t_0_xyz_x_xy, t_0_xyz_x_xz, t_0_xyz_x_yz, t_0_xyz_y_xy,\
                           t_0_xyz_y_xz, t_0_xyz_z_xy, t_0_xzz_0_xxz, t_0_xzz_0_xyz,\
                           t_0_xzz_0_xz, t_0_xzz_0_xzz, t_0_xzz_0_zz, t_0_xzz_x_xz,\
                           t_0_xzz_x_zz, t_0_xzz_y_xz, t_0_xzz_z_xz, t_0_yyy_0_xyy,\
                           t_0_yyy_0_yy, t_0_yyy_0_yyy, t_0_yyy_x_yy, t_0_yyy_y_yy,\
                           t_0_yyz_0_xyy, t_0_yyz_0_xyz, t_0_yyz_0_yy, t_0_yyz_0_yyy,\
                           t_0_yyz_0_yyz, t_0_yyz_0_yz, t_0_yyz_x_yy, t_0_yyz_x_yz,\
                           t_0_yyz_y_yy, t_0_yyz_y_yz, t_0_yyz_z_yy, t_0_yzz_0_xyz,\
                           t_0_yzz_0_xzz, t_0_yzz_0_yyz, t_0_yzz_0_yz, t_0_yzz_0_yzz,\
                           t_0_yzz_0_zz, t_0_yzz_x_yz, t_0_yzz_x_zz, t_0_yzz_y_yz, t_0_yzz_y_zz,\
                           t_0_yzz_z_yz, t_0_zzz_0_xzz, t_0_zzz_0_yzz, t_0_zzz_0_zz,\
                           t_0_zzz_0_zzz, t_0_zzz_x_zz, t_0_zzz_y_zz, t_0_zzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_zz[i] = t_0_zzz_0_zzz[i] - rcd_z[i] * t_0_zzz_0_zz[i];

        t_0_zzz_y_zz[i] = t_0_zzz_0_yzz[i] - rcd_y[i] * t_0_zzz_0_zz[i];

        t_0_zzz_x_zz[i] = t_0_zzz_0_xzz[i] - rcd_x[i] * t_0_zzz_0_zz[i];

        t_0_yzz_z_yz[i] = t_0_yzz_0_yzz[i] - rcd_z[i] * t_0_yzz_0_yz[i];

        t_0_yzz_y_zz[i] = t_0_yzz_0_yzz[i] - rcd_y[i] * t_0_yzz_0_zz[i];

        t_0_yzz_y_yz[i] = t_0_yzz_0_yyz[i] - rcd_y[i] * t_0_yzz_0_yz[i];

        t_0_yzz_x_zz[i] = t_0_yzz_0_xzz[i] - rcd_x[i] * t_0_yzz_0_zz[i];

        t_0_yzz_x_yz[i] = t_0_yzz_0_xyz[i] - rcd_x[i] * t_0_yzz_0_yz[i];

        t_0_yyz_z_yy[i] = t_0_yyz_0_yyz[i] - rcd_z[i] * t_0_yyz_0_yy[i];

        t_0_yyz_y_yz[i] = t_0_yyz_0_yyz[i] - rcd_y[i] * t_0_yyz_0_yz[i];

        t_0_yyz_y_yy[i] = t_0_yyz_0_yyy[i] - rcd_y[i] * t_0_yyz_0_yy[i];

        t_0_yyz_x_yz[i] = t_0_yyz_0_xyz[i] - rcd_x[i] * t_0_yyz_0_yz[i];

        t_0_yyz_x_yy[i] = t_0_yyz_0_xyy[i] - rcd_x[i] * t_0_yyz_0_yy[i];

        t_0_yyy_y_yy[i] = t_0_yyy_0_yyy[i] - rcd_y[i] * t_0_yyy_0_yy[i];

        t_0_yyy_x_yy[i] = t_0_yyy_0_xyy[i] - rcd_x[i] * t_0_yyy_0_yy[i];

        t_0_xzz_z_xz[i] = t_0_xzz_0_xzz[i] - rcd_z[i] * t_0_xzz_0_xz[i];

        t_0_xzz_y_xz[i] = t_0_xzz_0_xyz[i] - rcd_y[i] * t_0_xzz_0_xz[i];

        t_0_xzz_x_zz[i] = t_0_xzz_0_xzz[i] - rcd_x[i] * t_0_xzz_0_zz[i];

        t_0_xzz_x_xz[i] = t_0_xzz_0_xxz[i] - rcd_x[i] * t_0_xzz_0_xz[i];

        t_0_xyz_z_xy[i] = t_0_xyz_0_xyz[i] - rcd_z[i] * t_0_xyz_0_xy[i];

        t_0_xyz_y_xz[i] = t_0_xyz_0_xyz[i] - rcd_y[i] * t_0_xyz_0_xz[i];

        t_0_xyz_y_xy[i] = t_0_xyz_0_xyy[i] - rcd_y[i] * t_0_xyz_0_xy[i];

        t_0_xyz_x_yz[i] = t_0_xyz_0_xyz[i] - rcd_x[i] * t_0_xyz_0_yz[i];

        t_0_xyz_x_xz[i] = t_0_xyz_0_xxz[i] - rcd_x[i] * t_0_xyz_0_xz[i];

        t_0_xyz_x_xy[i] = t_0_xyz_0_xxy[i] - rcd_x[i] * t_0_xyz_0_xy[i];

        t_0_xyy_y_xy[i] = t_0_xyy_0_xyy[i] - rcd_y[i] * t_0_xyy_0_xy[i];

        t_0_xyy_x_yy[i] = t_0_xyy_0_xyy[i] - rcd_x[i] * t_0_xyy_0_yy[i];

        t_0_xyy_x_xy[i] = t_0_xyy_0_xxy[i] - rcd_x[i] * t_0_xyy_0_xy[i];

        t_0_xxz_z_xx[i] = t_0_xxz_0_xxz[i] - rcd_z[i] * t_0_xxz_0_xx[i];

        t_0_xxz_y_xx[i] = t_0_xxz_0_xxy[i] - rcd_y[i] * t_0_xxz_0_xx[i];

        t_0_xxz_x_xz[i] = t_0_xxz_0_xxz[i] - rcd_x[i] * t_0_xxz_0_xz[i];

        t_0_xxz_x_xx[i] = t_0_xxz_0_xxx[i] - rcd_x[i] * t_0_xxz_0_xx[i];

        t_0_xxy_y_xx[i] = t_0_xxy_0_xxy[i] - rcd_y[i] * t_0_xxy_0_xx[i];

        t_0_xxy_x_xy[i] = t_0_xxy_0_xxy[i] - rcd_x[i] * t_0_xxy_0_xy[i];

        t_0_xxy_x_xx[i] = t_0_xxy_0_xxx[i] - rcd_x[i] * t_0_xxy_0_xx[i];

        t_0_xxx_x_xx[i] = t_0_xxx_0_xxx[i] - rcd_x[i] * t_0_xxx_0_xx[i];
    }
}

template <typename T>
auto
compHostHRRForSFPD_V1(      BufferHostXY<T>&      intsBufferSFPD,
                      const BufferHostX<int32_t>& intsIndexesSFPD,
                      const BufferHostXY<T>&      intsBufferSFSD,
                      const BufferHostX<int32_t>& intsIndexesSFSD,
                      const BufferHostXY<T>&      intsBufferSFSF,
                      const BufferHostX<int32_t>& intsIndexesSFSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SFPD) integral components

    t_0_zzz_z_zz = intsBufferSFPD.data(intsIndexesSFPD(0));

    t_0_yzz_z_yz = intsBufferSFPD.data(intsIndexesSFPD(1));

    t_0_yzz_y_zz = intsBufferSFPD.data(intsIndexesSFPD(2));

    t_0_yyz_z_yy = intsBufferSFPD.data(intsIndexesSFPD(3));

    t_0_yyz_y_yz = intsBufferSFPD.data(intsIndexesSFPD(4));

    t_0_yyy_y_yy = intsBufferSFPD.data(intsIndexesSFPD(5));

    t_0_xzz_z_xz = intsBufferSFPD.data(intsIndexesSFPD(6));

    t_0_xzz_x_zz = intsBufferSFPD.data(intsIndexesSFPD(7));

    t_0_xyz_z_xy = intsBufferSFPD.data(intsIndexesSFPD(8));

    t_0_xyz_y_xz = intsBufferSFPD.data(intsIndexesSFPD(9));

    t_0_xyz_x_yz = intsBufferSFPD.data(intsIndexesSFPD(10));

    t_0_xyy_y_xy = intsBufferSFPD.data(intsIndexesSFPD(11));

    t_0_xyy_x_yy = intsBufferSFPD.data(intsIndexesSFPD(12));

    t_0_xxz_z_xx = intsBufferSFPD.data(intsIndexesSFPD(13));

    t_0_xxz_x_xz = intsBufferSFPD.data(intsIndexesSFPD(14));

    t_0_xxy_y_xx = intsBufferSFPD.data(intsIndexesSFPD(15));

    t_0_xxy_x_xy = intsBufferSFPD.data(intsIndexesSFPD(16));

    t_0_xxx_x_xx = intsBufferSFPD.data(intsIndexesSFPD(17));

    // set up (SFSD) integral components

    t_0_zzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(0));

    t_0_yzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(1));

    t_0_yzz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(2));

    t_0_yyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(3));

    t_0_yyz_0_yy = intsBufferSFSD.data(intsIndexesSFSD(4));

    t_0_yyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(5));

    t_0_xzz_0_zz = intsBufferSFSD.data(intsIndexesSFSD(6));

    t_0_xzz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(7));

    t_0_xyz_0_yz = intsBufferSFSD.data(intsIndexesSFSD(8));

    t_0_xyz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(9));

    t_0_xyz_0_xy = intsBufferSFSD.data(intsIndexesSFSD(10));

    t_0_xyy_0_yy = intsBufferSFSD.data(intsIndexesSFSD(11));

    t_0_xyy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(12));

    t_0_xxz_0_xz = intsBufferSFSD.data(intsIndexesSFSD(13));

    t_0_xxz_0_xx = intsBufferSFSD.data(intsIndexesSFSD(14));

    t_0_xxy_0_xy = intsBufferSFSD.data(intsIndexesSFSD(15));

    t_0_xxy_0_xx = intsBufferSFSD.data(intsIndexesSFSD(16));

    t_0_xxx_0_xx = intsBufferSFSD.data(intsIndexesSFSD(17));

    // set up (SFSF) integral components

    t_0_zzz_0_zzz = intsBufferSFSF.data(intsIndexesSFSF(0));

    t_0_yzz_0_yzz = intsBufferSFSF.data(intsIndexesSFSF(1));

    t_0_yyz_0_yyz = intsBufferSFSF.data(intsIndexesSFSF(2));

    t_0_yyy_0_yyy = intsBufferSFSF.data(intsIndexesSFSF(3));

    t_0_xzz_0_xzz = intsBufferSFSF.data(intsIndexesSFSF(4));

    t_0_xyz_0_xyz = intsBufferSFSF.data(intsIndexesSFSF(5));

    t_0_xyy_0_xyy = intsBufferSFSF.data(intsIndexesSFSF(6));

    t_0_xxz_0_xxz = intsBufferSFSF.data(intsIndexesSFSF(7));

    t_0_xxy_0_xxy = intsBufferSFSF.data(intsIndexesSFSF(8));

    t_0_xxx_0_xxx = intsBufferSFSF.data(intsIndexesSFSF(9));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxx_0_xx, t_0_xxx_0_xxx, t_0_xxx_x_xx,\
                           t_0_xxy_0_xx, t_0_xxy_0_xxy, t_0_xxy_0_xy, t_0_xxy_x_xy,\
                           t_0_xxy_y_xx, t_0_xxz_0_xx, t_0_xxz_0_xxz, t_0_xxz_0_xz,\
                           t_0_xxz_x_xz, t_0_xxz_z_xx, t_0_xyy_0_xy, t_0_xyy_0_xyy,\
                           t_0_xyy_0_yy, t_0_xyy_x_yy, t_0_xyy_y_xy, t_0_xyz_0_xy, t_0_xyz_0_xyz,\
                           t_0_xyz_0_xz, t_0_xyz_0_yz, t_0_xyz_x_yz, t_0_xyz_y_xz, t_0_xyz_z_xy,\
                           t_0_xzz_0_xz, t_0_xzz_0_xzz, t_0_xzz_0_zz, t_0_xzz_x_zz,\
                           t_0_xzz_z_xz, t_0_yyy_0_yy, t_0_yyy_0_yyy, t_0_yyy_y_yy,\
                           t_0_yyz_0_yy, t_0_yyz_0_yyz, t_0_yyz_0_yz, t_0_yyz_y_yz,\
                           t_0_yyz_z_yy, t_0_yzz_0_yz, t_0_yzz_0_yzz, t_0_yzz_0_zz,\
                           t_0_yzz_y_zz, t_0_yzz_z_yz, t_0_zzz_0_zz, t_0_zzz_0_zzz,\
                           t_0_zzz_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzz_z_zz[i] = t_0_zzz_0_zzz[i] - rcd_z[i] * t_0_zzz_0_zz[i];

        t_0_yzz_z_yz[i] = t_0_yzz_0_yzz[i] - rcd_z[i] * t_0_yzz_0_yz[i];

        t_0_yzz_y_zz[i] = t_0_yzz_0_yzz[i] - rcd_y[i] * t_0_yzz_0_zz[i];

        t_0_yyz_z_yy[i] = t_0_yyz_0_yyz[i] - rcd_z[i] * t_0_yyz_0_yy[i];

        t_0_yyz_y_yz[i] = t_0_yyz_0_yyz[i] - rcd_y[i] * t_0_yyz_0_yz[i];

        t_0_yyy_y_yy[i] = t_0_yyy_0_yyy[i] - rcd_y[i] * t_0_yyy_0_yy[i];

        t_0_xzz_z_xz[i] = t_0_xzz_0_xzz[i] - rcd_z[i] * t_0_xzz_0_xz[i];

        t_0_xzz_x_zz[i] = t_0_xzz_0_xzz[i] - rcd_x[i] * t_0_xzz_0_zz[i];

        t_0_xyz_z_xy[i] = t_0_xyz_0_xyz[i] - rcd_z[i] * t_0_xyz_0_xy[i];

        t_0_xyz_y_xz[i] = t_0_xyz_0_xyz[i] - rcd_y[i] * t_0_xyz_0_xz[i];

        t_0_xyz_x_yz[i] = t_0_xyz_0_xyz[i] - rcd_x[i] * t_0_xyz_0_yz[i];

        t_0_xyy_y_xy[i] = t_0_xyy_0_xyy[i] - rcd_y[i] * t_0_xyy_0_xy[i];

        t_0_xyy_x_yy[i] = t_0_xyy_0_xyy[i] - rcd_x[i] * t_0_xyy_0_yy[i];

        t_0_xxz_z_xx[i] = t_0_xxz_0_xxz[i] - rcd_z[i] * t_0_xxz_0_xx[i];

        t_0_xxz_x_xz[i] = t_0_xxz_0_xxz[i] - rcd_x[i] * t_0_xxz_0_xz[i];

        t_0_xxy_y_xx[i] = t_0_xxy_0_xxy[i] - rcd_y[i] * t_0_xxy_0_xx[i];

        t_0_xxy_x_xy[i] = t_0_xxy_0_xxy[i] - rcd_x[i] * t_0_xxy_0_xy[i];

        t_0_xxx_x_xx[i] = t_0_xxx_0_xxx[i] - rcd_x[i] * t_0_xxx_0_xx[i];
    }
}


} // derirec namespace
