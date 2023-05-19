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
compHostHRRForPDDD_V0(      BufferHostXY<T>&      intsBufferPDDD,
                      const BufferHostX<int32_t>& intsIndexesPDDD,
                      const BufferHostXY<T>&      intsBufferSDDD,
                      const BufferHostX<int32_t>& intsIndexesSDDD,
                      const BufferHostXY<T>&      intsBufferSFDD,
                      const BufferHostX<int32_t>& intsIndexesSFDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (PDDD) integral components

    t_z_zz_zz_zz = intsBufferPDDD.data(intsIndexesPDDD(0));

    t_z_yz_zz_yz = intsBufferPDDD.data(intsIndexesPDDD(1));

    t_z_yy_zz_yy = intsBufferPDDD.data(intsIndexesPDDD(2));

    t_z_xz_zz_xz = intsBufferPDDD.data(intsIndexesPDDD(3));

    t_z_xy_zz_xy = intsBufferPDDD.data(intsIndexesPDDD(4));

    t_z_xx_zz_xx = intsBufferPDDD.data(intsIndexesPDDD(5));

    t_y_zz_yz_zz = intsBufferPDDD.data(intsIndexesPDDD(6));

    t_y_zz_yy_zz = intsBufferPDDD.data(intsIndexesPDDD(7));

    t_y_yz_yz_yz = intsBufferPDDD.data(intsIndexesPDDD(8));

    t_y_yz_yy_yz = intsBufferPDDD.data(intsIndexesPDDD(9));

    t_y_yy_yz_yy = intsBufferPDDD.data(intsIndexesPDDD(10));

    t_y_yy_yy_yy = intsBufferPDDD.data(intsIndexesPDDD(11));

    t_y_xz_yz_xz = intsBufferPDDD.data(intsIndexesPDDD(12));

    t_y_xz_yy_xz = intsBufferPDDD.data(intsIndexesPDDD(13));

    t_y_xy_yz_xy = intsBufferPDDD.data(intsIndexesPDDD(14));

    t_y_xy_yy_xy = intsBufferPDDD.data(intsIndexesPDDD(15));

    t_y_xx_yz_xx = intsBufferPDDD.data(intsIndexesPDDD(16));

    t_y_xx_yy_xx = intsBufferPDDD.data(intsIndexesPDDD(17));

    t_x_zz_xz_zz = intsBufferPDDD.data(intsIndexesPDDD(18));

    t_x_zz_xy_zz = intsBufferPDDD.data(intsIndexesPDDD(19));

    t_x_zz_xx_zz = intsBufferPDDD.data(intsIndexesPDDD(20));

    t_x_yz_xz_yz = intsBufferPDDD.data(intsIndexesPDDD(21));

    t_x_yz_xy_yz = intsBufferPDDD.data(intsIndexesPDDD(22));

    t_x_yz_xx_yz = intsBufferPDDD.data(intsIndexesPDDD(23));

    t_x_yy_xz_yy = intsBufferPDDD.data(intsIndexesPDDD(24));

    t_x_yy_xy_yy = intsBufferPDDD.data(intsIndexesPDDD(25));

    t_x_yy_xx_yy = intsBufferPDDD.data(intsIndexesPDDD(26));

    t_x_xz_xz_xz = intsBufferPDDD.data(intsIndexesPDDD(27));

    t_x_xz_xy_xz = intsBufferPDDD.data(intsIndexesPDDD(28));

    t_x_xz_xx_xz = intsBufferPDDD.data(intsIndexesPDDD(29));

    t_x_xy_xz_xy = intsBufferPDDD.data(intsIndexesPDDD(30));

    t_x_xy_xy_xy = intsBufferPDDD.data(intsIndexesPDDD(31));

    t_x_xy_xx_xy = intsBufferPDDD.data(intsIndexesPDDD(32));

    t_x_xx_xz_xx = intsBufferPDDD.data(intsIndexesPDDD(33));

    t_x_xx_xy_xx = intsBufferPDDD.data(intsIndexesPDDD(34));

    t_x_xx_xx_xx = intsBufferPDDD.data(intsIndexesPDDD(35));

    // set up (SDDD) integral components

    t_0_zz_zz_zz = intsBufferSDDD.data(intsIndexesSDDD(0));

    t_0_zz_yz_zz = intsBufferSDDD.data(intsIndexesSDDD(1));

    t_0_zz_yy_zz = intsBufferSDDD.data(intsIndexesSDDD(2));

    t_0_zz_xz_zz = intsBufferSDDD.data(intsIndexesSDDD(3));

    t_0_zz_xy_zz = intsBufferSDDD.data(intsIndexesSDDD(4));

    t_0_zz_xx_zz = intsBufferSDDD.data(intsIndexesSDDD(5));

    t_0_yz_zz_yz = intsBufferSDDD.data(intsIndexesSDDD(6));

    t_0_yz_yz_yz = intsBufferSDDD.data(intsIndexesSDDD(7));

    t_0_yz_yy_yz = intsBufferSDDD.data(intsIndexesSDDD(8));

    t_0_yz_xz_yz = intsBufferSDDD.data(intsIndexesSDDD(9));

    t_0_yz_xy_yz = intsBufferSDDD.data(intsIndexesSDDD(10));

    t_0_yz_xx_yz = intsBufferSDDD.data(intsIndexesSDDD(11));

    t_0_yy_zz_yy = intsBufferSDDD.data(intsIndexesSDDD(12));

    t_0_yy_yz_yy = intsBufferSDDD.data(intsIndexesSDDD(13));

    t_0_yy_yy_yy = intsBufferSDDD.data(intsIndexesSDDD(14));

    t_0_yy_xz_yy = intsBufferSDDD.data(intsIndexesSDDD(15));

    t_0_yy_xy_yy = intsBufferSDDD.data(intsIndexesSDDD(16));

    t_0_yy_xx_yy = intsBufferSDDD.data(intsIndexesSDDD(17));

    t_0_xz_zz_xz = intsBufferSDDD.data(intsIndexesSDDD(18));

    t_0_xz_yz_xz = intsBufferSDDD.data(intsIndexesSDDD(19));

    t_0_xz_yy_xz = intsBufferSDDD.data(intsIndexesSDDD(20));

    t_0_xz_xz_xz = intsBufferSDDD.data(intsIndexesSDDD(21));

    t_0_xz_xy_xz = intsBufferSDDD.data(intsIndexesSDDD(22));

    t_0_xz_xx_xz = intsBufferSDDD.data(intsIndexesSDDD(23));

    t_0_xy_zz_xy = intsBufferSDDD.data(intsIndexesSDDD(24));

    t_0_xy_yz_xy = intsBufferSDDD.data(intsIndexesSDDD(25));

    t_0_xy_yy_xy = intsBufferSDDD.data(intsIndexesSDDD(26));

    t_0_xy_xz_xy = intsBufferSDDD.data(intsIndexesSDDD(27));

    t_0_xy_xy_xy = intsBufferSDDD.data(intsIndexesSDDD(28));

    t_0_xy_xx_xy = intsBufferSDDD.data(intsIndexesSDDD(29));

    t_0_xx_zz_xx = intsBufferSDDD.data(intsIndexesSDDD(30));

    t_0_xx_yz_xx = intsBufferSDDD.data(intsIndexesSDDD(31));

    t_0_xx_yy_xx = intsBufferSDDD.data(intsIndexesSDDD(32));

    t_0_xx_xz_xx = intsBufferSDDD.data(intsIndexesSDDD(33));

    t_0_xx_xy_xx = intsBufferSDDD.data(intsIndexesSDDD(34));

    t_0_xx_xx_xx = intsBufferSDDD.data(intsIndexesSDDD(35));

    // set up (SFDD) integral components

    t_0_zzz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(0));

    t_0_yzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(1));

    t_0_yzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(2));

    t_0_yzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(3));

    t_0_yyz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(4));

    t_0_yyz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(5));

    t_0_yyz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(6));

    t_0_yyy_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(7));

    t_0_yyy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(8));

    t_0_xzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(9));

    t_0_xzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(10));

    t_0_xzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(11));

    t_0_xzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(12));

    t_0_xyz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(13));

    t_0_xyz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(14));

    t_0_xyz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(15));

    t_0_xyz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(16));

    t_0_xyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(17));

    t_0_xyz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(18));

    t_0_xyy_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(19));

    t_0_xyy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(20));

    t_0_xyy_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(21));

    t_0_xyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(22));

    t_0_xyy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(23));

    t_0_xxz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(24));

    t_0_xxz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(25));

    t_0_xxz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(26));

    t_0_xxz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(27));

    t_0_xxy_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(28));

    t_0_xxy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(29));

    t_0_xxy_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(30));

    t_0_xxy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(31));

    t_0_xxy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(32));

    t_0_xxx_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(33));

    t_0_xxx_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(34));

    t_0_xxx_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(35));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_xx_xx_xx, t_0_xx_xy_xx, t_0_xx_xz_xx,\
                           t_0_xx_yy_xx, t_0_xx_yz_xx, t_0_xx_zz_xx, t_0_xxx_xx_xx,\
                           t_0_xxx_xy_xx, t_0_xxx_xz_xx, t_0_xxy_xx_xy, t_0_xxy_xy_xy,\
                           t_0_xxy_xz_xy, t_0_xxy_yy_xx, t_0_xxy_yz_xx, t_0_xxz_xx_xz,\
                           t_0_xxz_xy_xz, t_0_xxz_xz_xz, t_0_xxz_zz_xx, t_0_xy_xx_xy,\
                           t_0_xy_xy_xy, t_0_xy_xz_xy, t_0_xy_yy_xy, t_0_xy_yz_xy, t_0_xy_zz_xy,\
                           t_0_xyy_xx_yy, t_0_xyy_xy_yy, t_0_xyy_xz_yy, t_0_xyy_yy_xy,\
                           t_0_xyy_yz_xy, t_0_xyz_xx_yz, t_0_xyz_xy_yz, t_0_xyz_xz_yz,\
                           t_0_xyz_yy_xz, t_0_xyz_yz_xz, t_0_xyz_zz_xy, t_0_xz_xx_xz,\
                           t_0_xz_xy_xz, t_0_xz_xz_xz, t_0_xz_yy_xz, t_0_xz_yz_xz, t_0_xz_zz_xz,\
                           t_0_xzz_xx_zz, t_0_xzz_xy_zz, t_0_xzz_xz_zz, t_0_xzz_zz_xz,\
                           t_0_yy_xx_yy, t_0_yy_xy_yy, t_0_yy_xz_yy, t_0_yy_yy_yy, t_0_yy_yz_yy,\
                           t_0_yy_zz_yy, t_0_yyy_yy_yy, t_0_yyy_yz_yy, t_0_yyz_yy_yz,\
                           t_0_yyz_yz_yz, t_0_yyz_zz_yy, t_0_yz_xx_yz, t_0_yz_xy_yz,\
                           t_0_yz_xz_yz, t_0_yz_yy_yz, t_0_yz_yz_yz, t_0_yz_zz_yz, t_0_yzz_yy_zz,\
                           t_0_yzz_yz_zz, t_0_yzz_zz_yz, t_0_zz_xx_zz, t_0_zz_xy_zz,\
                           t_0_zz_xz_zz, t_0_zz_yy_zz, t_0_zz_yz_zz, t_0_zz_zz_zz, t_0_zzz_zz_zz,\
                           t_x_xx_xx_xx, t_x_xx_xy_xx, t_x_xx_xz_xx, t_x_xy_xx_xy, t_x_xy_xy_xy,\
                           t_x_xy_xz_xy, t_x_xz_xx_xz, t_x_xz_xy_xz, t_x_xz_xz_xz, t_x_yy_xx_yy,\
                           t_x_yy_xy_yy, t_x_yy_xz_yy, t_x_yz_xx_yz, t_x_yz_xy_yz, t_x_yz_xz_yz,\
                           t_x_zz_xx_zz, t_x_zz_xy_zz, t_x_zz_xz_zz, t_y_xx_yy_xx, t_y_xx_yz_xx,\
                           t_y_xy_yy_xy, t_y_xy_yz_xy, t_y_xz_yy_xz, t_y_xz_yz_xz, t_y_yy_yy_yy,\
                           t_y_yy_yz_yy, t_y_yz_yy_yz, t_y_yz_yz_yz, t_y_zz_yy_zz, t_y_zz_yz_zz,\
                           t_z_xx_zz_xx, t_z_xy_zz_xy, t_z_xz_zz_xz, t_z_yy_zz_yy, t_z_yz_zz_yz,\
                           t_z_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_zz_zz[i] = t_0_zzz_zz_zz[i] - rab_z[i] * t_0_zz_zz_zz[i];

        t_z_yz_zz_yz[i] = t_0_yzz_zz_yz[i] - rab_z[i] * t_0_yz_zz_yz[i];

        t_z_yy_zz_yy[i] = t_0_yyz_zz_yy[i] - rab_z[i] * t_0_yy_zz_yy[i];

        t_z_xz_zz_xz[i] = t_0_xzz_zz_xz[i] - rab_z[i] * t_0_xz_zz_xz[i];

        t_z_xy_zz_xy[i] = t_0_xyz_zz_xy[i] - rab_z[i] * t_0_xy_zz_xy[i];

        t_z_xx_zz_xx[i] = t_0_xxz_zz_xx[i] - rab_z[i] * t_0_xx_zz_xx[i];

        t_y_zz_yz_zz[i] = t_0_yzz_yz_zz[i] - rab_y[i] * t_0_zz_yz_zz[i];

        t_y_zz_yy_zz[i] = t_0_yzz_yy_zz[i] - rab_y[i] * t_0_zz_yy_zz[i];

        t_y_yz_yz_yz[i] = t_0_yyz_yz_yz[i] - rab_y[i] * t_0_yz_yz_yz[i];

        t_y_yz_yy_yz[i] = t_0_yyz_yy_yz[i] - rab_y[i] * t_0_yz_yy_yz[i];

        t_y_yy_yz_yy[i] = t_0_yyy_yz_yy[i] - rab_y[i] * t_0_yy_yz_yy[i];

        t_y_yy_yy_yy[i] = t_0_yyy_yy_yy[i] - rab_y[i] * t_0_yy_yy_yy[i];

        t_y_xz_yz_xz[i] = t_0_xyz_yz_xz[i] - rab_y[i] * t_0_xz_yz_xz[i];

        t_y_xz_yy_xz[i] = t_0_xyz_yy_xz[i] - rab_y[i] * t_0_xz_yy_xz[i];

        t_y_xy_yz_xy[i] = t_0_xyy_yz_xy[i] - rab_y[i] * t_0_xy_yz_xy[i];

        t_y_xy_yy_xy[i] = t_0_xyy_yy_xy[i] - rab_y[i] * t_0_xy_yy_xy[i];

        t_y_xx_yz_xx[i] = t_0_xxy_yz_xx[i] - rab_y[i] * t_0_xx_yz_xx[i];

        t_y_xx_yy_xx[i] = t_0_xxy_yy_xx[i] - rab_y[i] * t_0_xx_yy_xx[i];

        t_x_zz_xz_zz[i] = t_0_xzz_xz_zz[i] - rab_x[i] * t_0_zz_xz_zz[i];

        t_x_zz_xy_zz[i] = t_0_xzz_xy_zz[i] - rab_x[i] * t_0_zz_xy_zz[i];

        t_x_zz_xx_zz[i] = t_0_xzz_xx_zz[i] - rab_x[i] * t_0_zz_xx_zz[i];

        t_x_yz_xz_yz[i] = t_0_xyz_xz_yz[i] - rab_x[i] * t_0_yz_xz_yz[i];

        t_x_yz_xy_yz[i] = t_0_xyz_xy_yz[i] - rab_x[i] * t_0_yz_xy_yz[i];

        t_x_yz_xx_yz[i] = t_0_xyz_xx_yz[i] - rab_x[i] * t_0_yz_xx_yz[i];

        t_x_yy_xz_yy[i] = t_0_xyy_xz_yy[i] - rab_x[i] * t_0_yy_xz_yy[i];

        t_x_yy_xy_yy[i] = t_0_xyy_xy_yy[i] - rab_x[i] * t_0_yy_xy_yy[i];

        t_x_yy_xx_yy[i] = t_0_xyy_xx_yy[i] - rab_x[i] * t_0_yy_xx_yy[i];

        t_x_xz_xz_xz[i] = t_0_xxz_xz_xz[i] - rab_x[i] * t_0_xz_xz_xz[i];

        t_x_xz_xy_xz[i] = t_0_xxz_xy_xz[i] - rab_x[i] * t_0_xz_xy_xz[i];

        t_x_xz_xx_xz[i] = t_0_xxz_xx_xz[i] - rab_x[i] * t_0_xz_xx_xz[i];

        t_x_xy_xz_xy[i] = t_0_xxy_xz_xy[i] - rab_x[i] * t_0_xy_xz_xy[i];

        t_x_xy_xy_xy[i] = t_0_xxy_xy_xy[i] - rab_x[i] * t_0_xy_xy_xy[i];

        t_x_xy_xx_xy[i] = t_0_xxy_xx_xy[i] - rab_x[i] * t_0_xy_xx_xy[i];

        t_x_xx_xz_xx[i] = t_0_xxx_xz_xx[i] - rab_x[i] * t_0_xx_xz_xx[i];

        t_x_xx_xy_xx[i] = t_0_xxx_xy_xx[i] - rab_x[i] * t_0_xx_xy_xx[i];

        t_x_xx_xx_xx[i] = t_0_xxx_xx_xx[i] - rab_x[i] * t_0_xx_xx_xx[i];
    }
}


} // derirec namespace
