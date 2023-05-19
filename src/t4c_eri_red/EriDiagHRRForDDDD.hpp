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
compHostHRRForDDDD_V0(      BufferHostXY<T>&      intsBufferDDDD,
                      const BufferHostX<int32_t>& intsIndexesDDDD,
                      const BufferHostXY<T>&      intsBufferPDDD,
                      const BufferHostX<int32_t>& intsIndexesPDDD,
                      const BufferHostXY<T>&      intsBufferPFDD,
                      const BufferHostX<int32_t>& intsIndexesPFDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDDD) integral components

    t_zz_zz_zz_zz = intsBufferDDDD.data(intsIndexesDDDD(0));

    t_zz_yz_zz_yz = intsBufferDDDD.data(intsIndexesDDDD(1));

    t_zz_yy_zz_yy = intsBufferDDDD.data(intsIndexesDDDD(2));

    t_zz_xz_zz_xz = intsBufferDDDD.data(intsIndexesDDDD(3));

    t_zz_xy_zz_xy = intsBufferDDDD.data(intsIndexesDDDD(4));

    t_zz_xx_zz_xx = intsBufferDDDD.data(intsIndexesDDDD(5));

    t_yz_zz_yz_zz = intsBufferDDDD.data(intsIndexesDDDD(6));

    t_yz_yz_yz_yz = intsBufferDDDD.data(intsIndexesDDDD(7));

    t_yz_yy_yz_yy = intsBufferDDDD.data(intsIndexesDDDD(8));

    t_yz_xz_yz_xz = intsBufferDDDD.data(intsIndexesDDDD(9));

    t_yz_xy_yz_xy = intsBufferDDDD.data(intsIndexesDDDD(10));

    t_yz_xx_yz_xx = intsBufferDDDD.data(intsIndexesDDDD(11));

    t_yy_zz_yy_zz = intsBufferDDDD.data(intsIndexesDDDD(12));

    t_yy_yz_yy_yz = intsBufferDDDD.data(intsIndexesDDDD(13));

    t_yy_yy_yy_yy = intsBufferDDDD.data(intsIndexesDDDD(14));

    t_yy_xz_yy_xz = intsBufferDDDD.data(intsIndexesDDDD(15));

    t_yy_xy_yy_xy = intsBufferDDDD.data(intsIndexesDDDD(16));

    t_yy_xx_yy_xx = intsBufferDDDD.data(intsIndexesDDDD(17));

    t_xz_zz_xz_zz = intsBufferDDDD.data(intsIndexesDDDD(18));

    t_xz_yz_xz_yz = intsBufferDDDD.data(intsIndexesDDDD(19));

    t_xz_yy_xz_yy = intsBufferDDDD.data(intsIndexesDDDD(20));

    t_xz_xz_xz_xz = intsBufferDDDD.data(intsIndexesDDDD(21));

    t_xz_xy_xz_xy = intsBufferDDDD.data(intsIndexesDDDD(22));

    t_xz_xx_xz_xx = intsBufferDDDD.data(intsIndexesDDDD(23));

    t_xy_zz_xy_zz = intsBufferDDDD.data(intsIndexesDDDD(24));

    t_xy_yz_xy_yz = intsBufferDDDD.data(intsIndexesDDDD(25));

    t_xy_yy_xy_yy = intsBufferDDDD.data(intsIndexesDDDD(26));

    t_xy_xz_xy_xz = intsBufferDDDD.data(intsIndexesDDDD(27));

    t_xy_xy_xy_xy = intsBufferDDDD.data(intsIndexesDDDD(28));

    t_xy_xx_xy_xx = intsBufferDDDD.data(intsIndexesDDDD(29));

    t_xx_zz_xx_zz = intsBufferDDDD.data(intsIndexesDDDD(30));

    t_xx_yz_xx_yz = intsBufferDDDD.data(intsIndexesDDDD(31));

    t_xx_yy_xx_yy = intsBufferDDDD.data(intsIndexesDDDD(32));

    t_xx_xz_xx_xz = intsBufferDDDD.data(intsIndexesDDDD(33));

    t_xx_xy_xx_xy = intsBufferDDDD.data(intsIndexesDDDD(34));

    t_xx_xx_xx_xx = intsBufferDDDD.data(intsIndexesDDDD(35));

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

    // set up (PFDD) integral components

    t_z_zzz_zz_zz = intsBufferPFDD.data(intsIndexesPFDD(0));

    t_z_yzz_zz_yz = intsBufferPFDD.data(intsIndexesPFDD(1));

    t_z_yyz_zz_yy = intsBufferPFDD.data(intsIndexesPFDD(2));

    t_z_xzz_zz_xz = intsBufferPFDD.data(intsIndexesPFDD(3));

    t_z_xyz_zz_xy = intsBufferPFDD.data(intsIndexesPFDD(4));

    t_z_xxz_zz_xx = intsBufferPFDD.data(intsIndexesPFDD(5));

    t_y_zzz_yz_zz = intsBufferPFDD.data(intsIndexesPFDD(6));

    t_y_yzz_yz_yz = intsBufferPFDD.data(intsIndexesPFDD(7));

    t_y_yzz_yy_zz = intsBufferPFDD.data(intsIndexesPFDD(8));

    t_y_yyz_yz_yy = intsBufferPFDD.data(intsIndexesPFDD(9));

    t_y_yyz_yy_yz = intsBufferPFDD.data(intsIndexesPFDD(10));

    t_y_yyy_yy_yy = intsBufferPFDD.data(intsIndexesPFDD(11));

    t_y_xzz_yz_xz = intsBufferPFDD.data(intsIndexesPFDD(12));

    t_y_xyz_yz_xy = intsBufferPFDD.data(intsIndexesPFDD(13));

    t_y_xyz_yy_xz = intsBufferPFDD.data(intsIndexesPFDD(14));

    t_y_xyy_yy_xy = intsBufferPFDD.data(intsIndexesPFDD(15));

    t_y_xxz_yz_xx = intsBufferPFDD.data(intsIndexesPFDD(16));

    t_y_xxy_yy_xx = intsBufferPFDD.data(intsIndexesPFDD(17));

    t_x_zzz_xz_zz = intsBufferPFDD.data(intsIndexesPFDD(18));

    t_x_yzz_xz_yz = intsBufferPFDD.data(intsIndexesPFDD(19));

    t_x_yzz_xy_zz = intsBufferPFDD.data(intsIndexesPFDD(20));

    t_x_yyz_xz_yy = intsBufferPFDD.data(intsIndexesPFDD(21));

    t_x_yyz_xy_yz = intsBufferPFDD.data(intsIndexesPFDD(22));

    t_x_yyy_xy_yy = intsBufferPFDD.data(intsIndexesPFDD(23));

    t_x_xzz_xz_xz = intsBufferPFDD.data(intsIndexesPFDD(24));

    t_x_xzz_xx_zz = intsBufferPFDD.data(intsIndexesPFDD(25));

    t_x_xyz_xz_xy = intsBufferPFDD.data(intsIndexesPFDD(26));

    t_x_xyz_xy_xz = intsBufferPFDD.data(intsIndexesPFDD(27));

    t_x_xyz_xx_yz = intsBufferPFDD.data(intsIndexesPFDD(28));

    t_x_xyy_xy_xy = intsBufferPFDD.data(intsIndexesPFDD(29));

    t_x_xyy_xx_yy = intsBufferPFDD.data(intsIndexesPFDD(30));

    t_x_xxz_xz_xx = intsBufferPFDD.data(intsIndexesPFDD(31));

    t_x_xxz_xx_xz = intsBufferPFDD.data(intsIndexesPFDD(32));

    t_x_xxy_xy_xx = intsBufferPFDD.data(intsIndexesPFDD(33));

    t_x_xxy_xx_xy = intsBufferPFDD.data(intsIndexesPFDD(34));

    t_x_xxx_xx_xx = intsBufferPFDD.data(intsIndexesPFDD(35));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_x_xx_xx_xx, t_x_xx_xy_xx, t_x_xx_xz_xx,\
                           t_x_xxx_xx_xx, t_x_xxy_xx_xy, t_x_xxy_xy_xx, t_x_xxz_xx_xz,\
                           t_x_xxz_xz_xx, t_x_xy_xx_xy, t_x_xy_xy_xy, t_x_xy_xz_xy,\
                           t_x_xyy_xx_yy, t_x_xyy_xy_xy, t_x_xyz_xx_yz, t_x_xyz_xy_xz,\
                           t_x_xyz_xz_xy, t_x_xz_xx_xz, t_x_xz_xy_xz, t_x_xz_xz_xz,\
                           t_x_xzz_xx_zz, t_x_xzz_xz_xz, t_x_yy_xx_yy, t_x_yy_xy_yy,\
                           t_x_yy_xz_yy, t_x_yyy_xy_yy, t_x_yyz_xy_yz, t_x_yyz_xz_yy,\
                           t_x_yz_xx_yz, t_x_yz_xy_yz, t_x_yz_xz_yz, t_x_yzz_xy_zz,\
                           t_x_yzz_xz_yz, t_x_zz_xx_zz, t_x_zz_xy_zz, t_x_zz_xz_zz,\
                           t_x_zzz_xz_zz, t_xx_xx_xx_xx, t_xx_xy_xx_xy, t_xx_xz_xx_xz,\
                           t_xx_yy_xx_yy, t_xx_yz_xx_yz, t_xx_zz_xx_zz, t_xy_xx_xy_xx,\
                           t_xy_xy_xy_xy, t_xy_xz_xy_xz, t_xy_yy_xy_yy, t_xy_yz_xy_yz,\
                           t_xy_zz_xy_zz, t_xz_xx_xz_xx, t_xz_xy_xz_xy, t_xz_xz_xz_xz,\
                           t_xz_yy_xz_yy, t_xz_yz_xz_yz, t_xz_zz_xz_zz, t_y_xx_yy_xx,\
                           t_y_xx_yz_xx, t_y_xxy_yy_xx, t_y_xxz_yz_xx, t_y_xy_yy_xy,\
                           t_y_xy_yz_xy, t_y_xyy_yy_xy, t_y_xyz_yy_xz, t_y_xyz_yz_xy,\
                           t_y_xz_yy_xz, t_y_xz_yz_xz, t_y_xzz_yz_xz, t_y_yy_yy_yy,\
                           t_y_yy_yz_yy, t_y_yyy_yy_yy, t_y_yyz_yy_yz, t_y_yyz_yz_yy,\
                           t_y_yz_yy_yz, t_y_yz_yz_yz, t_y_yzz_yy_zz, t_y_yzz_yz_yz,\
                           t_y_zz_yy_zz, t_y_zz_yz_zz, t_y_zzz_yz_zz, t_yy_xx_yy_xx,\
                           t_yy_xy_yy_xy, t_yy_xz_yy_xz, t_yy_yy_yy_yy, t_yy_yz_yy_yz,\
                           t_yy_zz_yy_zz, t_yz_xx_yz_xx, t_yz_xy_yz_xy, t_yz_xz_yz_xz,\
                           t_yz_yy_yz_yy, t_yz_yz_yz_yz, t_yz_zz_yz_zz, t_z_xx_zz_xx,\
                           t_z_xxz_zz_xx, t_z_xy_zz_xy, t_z_xyz_zz_xy, t_z_xz_zz_xz,\
                           t_z_xzz_zz_xz, t_z_yy_zz_yy, t_z_yyz_zz_yy, t_z_yz_zz_yz,\
                           t_z_yzz_zz_yz, t_z_zz_zz_zz, t_z_zzz_zz_zz, t_zz_xx_zz_xx,\
                           t_zz_xy_zz_xy, t_zz_xz_zz_xz, t_zz_yy_zz_yy, t_zz_yz_zz_yz,\
                           t_zz_zz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_zz_zz[i] = t_z_zzz_zz_zz[i] - rab_z[i] * t_z_zz_zz_zz[i];

        t_zz_yz_zz_yz[i] = t_z_yzz_zz_yz[i] - rab_z[i] * t_z_yz_zz_yz[i];

        t_zz_yy_zz_yy[i] = t_z_yyz_zz_yy[i] - rab_z[i] * t_z_yy_zz_yy[i];

        t_zz_xz_zz_xz[i] = t_z_xzz_zz_xz[i] - rab_z[i] * t_z_xz_zz_xz[i];

        t_zz_xy_zz_xy[i] = t_z_xyz_zz_xy[i] - rab_z[i] * t_z_xy_zz_xy[i];

        t_zz_xx_zz_xx[i] = t_z_xxz_zz_xx[i] - rab_z[i] * t_z_xx_zz_xx[i];

        t_yz_zz_yz_zz[i] = t_y_zzz_yz_zz[i] - rab_z[i] * t_y_zz_yz_zz[i];

        t_yz_yz_yz_yz[i] = t_y_yzz_yz_yz[i] - rab_z[i] * t_y_yz_yz_yz[i];

        t_yz_yy_yz_yy[i] = t_y_yyz_yz_yy[i] - rab_z[i] * t_y_yy_yz_yy[i];

        t_yz_xz_yz_xz[i] = t_y_xzz_yz_xz[i] - rab_z[i] * t_y_xz_yz_xz[i];

        t_yz_xy_yz_xy[i] = t_y_xyz_yz_xy[i] - rab_z[i] * t_y_xy_yz_xy[i];

        t_yz_xx_yz_xx[i] = t_y_xxz_yz_xx[i] - rab_z[i] * t_y_xx_yz_xx[i];

        t_yy_zz_yy_zz[i] = t_y_yzz_yy_zz[i] - rab_y[i] * t_y_zz_yy_zz[i];

        t_yy_yz_yy_yz[i] = t_y_yyz_yy_yz[i] - rab_y[i] * t_y_yz_yy_yz[i];

        t_yy_yy_yy_yy[i] = t_y_yyy_yy_yy[i] - rab_y[i] * t_y_yy_yy_yy[i];

        t_yy_xz_yy_xz[i] = t_y_xyz_yy_xz[i] - rab_y[i] * t_y_xz_yy_xz[i];

        t_yy_xy_yy_xy[i] = t_y_xyy_yy_xy[i] - rab_y[i] * t_y_xy_yy_xy[i];

        t_yy_xx_yy_xx[i] = t_y_xxy_yy_xx[i] - rab_y[i] * t_y_xx_yy_xx[i];

        t_xz_zz_xz_zz[i] = t_x_zzz_xz_zz[i] - rab_z[i] * t_x_zz_xz_zz[i];

        t_xz_yz_xz_yz[i] = t_x_yzz_xz_yz[i] - rab_z[i] * t_x_yz_xz_yz[i];

        t_xz_yy_xz_yy[i] = t_x_yyz_xz_yy[i] - rab_z[i] * t_x_yy_xz_yy[i];

        t_xz_xz_xz_xz[i] = t_x_xzz_xz_xz[i] - rab_z[i] * t_x_xz_xz_xz[i];

        t_xz_xy_xz_xy[i] = t_x_xyz_xz_xy[i] - rab_z[i] * t_x_xy_xz_xy[i];

        t_xz_xx_xz_xx[i] = t_x_xxz_xz_xx[i] - rab_z[i] * t_x_xx_xz_xx[i];

        t_xy_zz_xy_zz[i] = t_x_yzz_xy_zz[i] - rab_y[i] * t_x_zz_xy_zz[i];

        t_xy_yz_xy_yz[i] = t_x_yyz_xy_yz[i] - rab_y[i] * t_x_yz_xy_yz[i];

        t_xy_yy_xy_yy[i] = t_x_yyy_xy_yy[i] - rab_y[i] * t_x_yy_xy_yy[i];

        t_xy_xz_xy_xz[i] = t_x_xyz_xy_xz[i] - rab_y[i] * t_x_xz_xy_xz[i];

        t_xy_xy_xy_xy[i] = t_x_xyy_xy_xy[i] - rab_y[i] * t_x_xy_xy_xy[i];

        t_xy_xx_xy_xx[i] = t_x_xxy_xy_xx[i] - rab_y[i] * t_x_xx_xy_xx[i];

        t_xx_zz_xx_zz[i] = t_x_xzz_xx_zz[i] - rab_x[i] * t_x_zz_xx_zz[i];

        t_xx_yz_xx_yz[i] = t_x_xyz_xx_yz[i] - rab_x[i] * t_x_yz_xx_yz[i];

        t_xx_yy_xx_yy[i] = t_x_xyy_xx_yy[i] - rab_x[i] * t_x_yy_xx_yy[i];

        t_xx_xz_xx_xz[i] = t_x_xxz_xx_xz[i] - rab_x[i] * t_x_xz_xx_xz[i];

        t_xx_xy_xx_xy[i] = t_x_xxy_xx_xy[i] - rab_x[i] * t_x_xy_xx_xy[i];

        t_xx_xx_xx_xx[i] = t_x_xxx_xx_xx[i] - rab_x[i] * t_x_xx_xx_xx[i];
    }
}


} // derirec namespace
