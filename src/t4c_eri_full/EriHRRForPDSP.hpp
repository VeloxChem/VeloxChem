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
compHostHRRForPDSP_V0(      BufferHostXY<T>&      intsBufferPDSP,
                      const BufferHostX<int32_t>& intsIndexesPDSP,
                      const BufferHostXY<T>&      intsBufferSDSP,
                      const BufferHostX<int32_t>& intsIndexesSDSP,
                      const BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    #pragma omp simd align(rab_y, rab_z, t_0_xx_0_x, t_0_xx_0_y, t_0_xx_0_z, t_0_xxy_0_x,\
                           t_0_xxy_0_y, t_0_xxy_0_z, t_0_xxz_0_x, t_0_xxz_0_y, t_0_xxz_0_z,\
                           t_0_xy_0_x, t_0_xy_0_y, t_0_xy_0_z, t_0_xyy_0_x, t_0_xyy_0_y,\
                           t_0_xyy_0_z, t_0_xyz_0_x, t_0_xyz_0_y, t_0_xyz_0_z, t_0_xz_0_x,\
                           t_0_xz_0_y, t_0_xz_0_z, t_0_xzz_0_x, t_0_xzz_0_y, t_0_xzz_0_z,\
                           t_0_yy_0_x, t_0_yy_0_y, t_0_yy_0_z, t_0_yyy_0_x, t_0_yyy_0_y,\
                           t_0_yyy_0_z, t_0_yyz_0_x, t_0_yyz_0_y, t_0_yyz_0_z, t_0_yz_0_x,\
                           t_0_yz_0_y, t_0_yz_0_z, t_0_yzz_0_x, t_0_yzz_0_y, t_0_yzz_0_z,\
                           t_0_zz_0_x, t_0_zz_0_y, t_0_zz_0_z, t_0_zzz_0_x, t_0_zzz_0_y,\
                           t_0_zzz_0_z, t_y_xx_0_x, t_y_xx_0_y, t_y_xx_0_z, t_y_xy_0_x,\
                           t_y_xy_0_y, t_y_xy_0_z, t_y_xz_0_x, t_y_xz_0_y, t_y_xz_0_z,\
                           t_y_yy_0_x, t_y_yy_0_y, t_y_yy_0_z, t_y_yz_0_x, t_y_yz_0_y,\
                           t_y_yz_0_z, t_y_zz_0_x, t_y_zz_0_y, t_y_zz_0_z, t_z_xx_0_x,\
                           t_z_xx_0_y, t_z_xx_0_z, t_z_xy_0_x, t_z_xy_0_y, t_z_xy_0_z,\
                           t_z_xz_0_x, t_z_xz_0_y, t_z_xz_0_z, t_z_yy_0_x, t_z_yy_0_y,\
                           t_z_yy_0_z, t_z_yz_0_x, t_z_yz_0_y, t_z_yz_0_z, t_z_zz_0_x,\
                           t_z_zz_0_y, t_z_zz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_0_z[i] = t_0_zzz_0_z[i] - rab_z[i] * t_0_zz_0_z[i];

        t_z_zz_0_y[i] = t_0_zzz_0_y[i] - rab_z[i] * t_0_zz_0_y[i];

        t_z_zz_0_x[i] = t_0_zzz_0_x[i] - rab_z[i] * t_0_zz_0_x[i];

        t_z_yz_0_z[i] = t_0_yzz_0_z[i] - rab_z[i] * t_0_yz_0_z[i];

        t_z_yz_0_y[i] = t_0_yzz_0_y[i] - rab_z[i] * t_0_yz_0_y[i];

        t_z_yz_0_x[i] = t_0_yzz_0_x[i] - rab_z[i] * t_0_yz_0_x[i];

        t_z_yy_0_z[i] = t_0_yyz_0_z[i] - rab_z[i] * t_0_yy_0_z[i];

        t_z_yy_0_y[i] = t_0_yyz_0_y[i] - rab_z[i] * t_0_yy_0_y[i];

        t_z_yy_0_x[i] = t_0_yyz_0_x[i] - rab_z[i] * t_0_yy_0_x[i];

        t_z_xz_0_z[i] = t_0_xzz_0_z[i] - rab_z[i] * t_0_xz_0_z[i];

        t_z_xz_0_y[i] = t_0_xzz_0_y[i] - rab_z[i] * t_0_xz_0_y[i];

        t_z_xz_0_x[i] = t_0_xzz_0_x[i] - rab_z[i] * t_0_xz_0_x[i];

        t_z_xy_0_z[i] = t_0_xyz_0_z[i] - rab_z[i] * t_0_xy_0_z[i];

        t_z_xy_0_y[i] = t_0_xyz_0_y[i] - rab_z[i] * t_0_xy_0_y[i];

        t_z_xy_0_x[i] = t_0_xyz_0_x[i] - rab_z[i] * t_0_xy_0_x[i];

        t_z_xx_0_z[i] = t_0_xxz_0_z[i] - rab_z[i] * t_0_xx_0_z[i];

        t_z_xx_0_y[i] = t_0_xxz_0_y[i] - rab_z[i] * t_0_xx_0_y[i];

        t_z_xx_0_x[i] = t_0_xxz_0_x[i] - rab_z[i] * t_0_xx_0_x[i];

        t_y_zz_0_z[i] = t_0_yzz_0_z[i] - rab_y[i] * t_0_zz_0_z[i];

        t_y_zz_0_y[i] = t_0_yzz_0_y[i] - rab_y[i] * t_0_zz_0_y[i];

        t_y_zz_0_x[i] = t_0_yzz_0_x[i] - rab_y[i] * t_0_zz_0_x[i];

        t_y_yz_0_z[i] = t_0_yyz_0_z[i] - rab_y[i] * t_0_yz_0_z[i];

        t_y_yz_0_y[i] = t_0_yyz_0_y[i] - rab_y[i] * t_0_yz_0_y[i];

        t_y_yz_0_x[i] = t_0_yyz_0_x[i] - rab_y[i] * t_0_yz_0_x[i];

        t_y_yy_0_z[i] = t_0_yyy_0_z[i] - rab_y[i] * t_0_yy_0_z[i];

        t_y_yy_0_y[i] = t_0_yyy_0_y[i] - rab_y[i] * t_0_yy_0_y[i];

        t_y_yy_0_x[i] = t_0_yyy_0_x[i] - rab_y[i] * t_0_yy_0_x[i];

        t_y_xz_0_z[i] = t_0_xyz_0_z[i] - rab_y[i] * t_0_xz_0_z[i];

        t_y_xz_0_y[i] = t_0_xyz_0_y[i] - rab_y[i] * t_0_xz_0_y[i];

        t_y_xz_0_x[i] = t_0_xyz_0_x[i] - rab_y[i] * t_0_xz_0_x[i];

        t_y_xy_0_z[i] = t_0_xyy_0_z[i] - rab_y[i] * t_0_xy_0_z[i];

        t_y_xy_0_y[i] = t_0_xyy_0_y[i] - rab_y[i] * t_0_xy_0_y[i];

        t_y_xy_0_x[i] = t_0_xyy_0_x[i] - rab_y[i] * t_0_xy_0_x[i];

        t_y_xx_0_z[i] = t_0_xxy_0_z[i] - rab_y[i] * t_0_xx_0_z[i];

        t_y_xx_0_y[i] = t_0_xxy_0_y[i] - rab_y[i] * t_0_xx_0_y[i];

        t_y_xx_0_x[i] = t_0_xxy_0_x[i] - rab_y[i] * t_0_xx_0_x[i];
    }

    #pragma omp simd align(rab_x, t_0_xx_0_x, t_0_xx_0_y, t_0_xx_0_z, t_0_xxx_0_x, t_0_xxx_0_y,\
                           t_0_xxx_0_z, t_0_xxy_0_x, t_0_xxy_0_y, t_0_xxy_0_z, t_0_xxz_0_x,\
                           t_0_xxz_0_y, t_0_xxz_0_z, t_0_xy_0_x, t_0_xy_0_y, t_0_xy_0_z,\
                           t_0_xyy_0_x, t_0_xyy_0_y, t_0_xyy_0_z, t_0_xyz_0_x, t_0_xyz_0_y,\
                           t_0_xyz_0_z, t_0_xz_0_x, t_0_xz_0_y, t_0_xz_0_z, t_0_xzz_0_x,\
                           t_0_xzz_0_y, t_0_xzz_0_z, t_0_yy_0_x, t_0_yy_0_y, t_0_yy_0_z,\
                           t_0_yz_0_x, t_0_yz_0_y, t_0_yz_0_z, t_0_zz_0_x, t_0_zz_0_y,\
                           t_0_zz_0_z, t_x_xx_0_x, t_x_xx_0_y, t_x_xx_0_z, t_x_xy_0_x,\
                           t_x_xy_0_y, t_x_xy_0_z, t_x_xz_0_x, t_x_xz_0_y, t_x_xz_0_z,\
                           t_x_yy_0_x, t_x_yy_0_y, t_x_yy_0_z, t_x_yz_0_x, t_x_yz_0_y,\
                           t_x_yz_0_z, t_x_zz_0_x, t_x_zz_0_y, t_x_zz_0_z : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_x_zz_0_z[i] = t_0_xzz_0_z[i] - rab_x[i] * t_0_zz_0_z[i];

        t_x_zz_0_y[i] = t_0_xzz_0_y[i] - rab_x[i] * t_0_zz_0_y[i];

        t_x_zz_0_x[i] = t_0_xzz_0_x[i] - rab_x[i] * t_0_zz_0_x[i];

        t_x_yz_0_z[i] = t_0_xyz_0_z[i] - rab_x[i] * t_0_yz_0_z[i];

        t_x_yz_0_y[i] = t_0_xyz_0_y[i] - rab_x[i] * t_0_yz_0_y[i];

        t_x_yz_0_x[i] = t_0_xyz_0_x[i] - rab_x[i] * t_0_yz_0_x[i];

        t_x_yy_0_z[i] = t_0_xyy_0_z[i] - rab_x[i] * t_0_yy_0_z[i];

        t_x_yy_0_y[i] = t_0_xyy_0_y[i] - rab_x[i] * t_0_yy_0_y[i];

        t_x_yy_0_x[i] = t_0_xyy_0_x[i] - rab_x[i] * t_0_yy_0_x[i];

        t_x_xz_0_z[i] = t_0_xxz_0_z[i] - rab_x[i] * t_0_xz_0_z[i];

        t_x_xz_0_y[i] = t_0_xxz_0_y[i] - rab_x[i] * t_0_xz_0_y[i];

        t_x_xz_0_x[i] = t_0_xxz_0_x[i] - rab_x[i] * t_0_xz_0_x[i];

        t_x_xy_0_z[i] = t_0_xxy_0_z[i] - rab_x[i] * t_0_xy_0_z[i];

        t_x_xy_0_y[i] = t_0_xxy_0_y[i] - rab_x[i] * t_0_xy_0_y[i];

        t_x_xy_0_x[i] = t_0_xxy_0_x[i] - rab_x[i] * t_0_xy_0_x[i];

        t_x_xx_0_z[i] = t_0_xxx_0_z[i] - rab_x[i] * t_0_xx_0_z[i];

        t_x_xx_0_y[i] = t_0_xxx_0_y[i] - rab_x[i] * t_0_xx_0_y[i];

        t_x_xx_0_x[i] = t_0_xxx_0_x[i] - rab_x[i] * t_0_xx_0_x[i];
    }
}


} // derirec namespace
