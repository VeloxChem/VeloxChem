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
compHostHRRForPFDD_V0(      BufferHostXY<T>&      intsBufferPFDD,
                      const BufferHostX<int32_t>& intsIndexesPFDD,
                      const BufferHostXY<T>&      intsBufferSFDD,
                      const BufferHostX<int32_t>& intsIndexesSFDD,
                      const BufferHostXY<T>&      intsBufferSGDD,
                      const BufferHostX<int32_t>& intsIndexesSGDD,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SFDD) integral components

    t_0_zzz_zz_zz = intsBufferSFDD.data(intsIndexesSFDD(0));

    t_0_zzz_yz_zz = intsBufferSFDD.data(intsIndexesSFDD(1));

    t_0_zzz_xz_zz = intsBufferSFDD.data(intsIndexesSFDD(2));

    t_0_yzz_zz_yz = intsBufferSFDD.data(intsIndexesSFDD(3));

    t_0_yzz_yz_yz = intsBufferSFDD.data(intsIndexesSFDD(4));

    t_0_yzz_yy_zz = intsBufferSFDD.data(intsIndexesSFDD(5));

    t_0_yzz_xz_yz = intsBufferSFDD.data(intsIndexesSFDD(6));

    t_0_yzz_xy_zz = intsBufferSFDD.data(intsIndexesSFDD(7));

    t_0_yyz_zz_yy = intsBufferSFDD.data(intsIndexesSFDD(8));

    t_0_yyz_yz_yy = intsBufferSFDD.data(intsIndexesSFDD(9));

    t_0_yyz_yy_yz = intsBufferSFDD.data(intsIndexesSFDD(10));

    t_0_yyz_xz_yy = intsBufferSFDD.data(intsIndexesSFDD(11));

    t_0_yyz_xy_yz = intsBufferSFDD.data(intsIndexesSFDD(12));

    t_0_yyy_yy_yy = intsBufferSFDD.data(intsIndexesSFDD(13));

    t_0_yyy_xy_yy = intsBufferSFDD.data(intsIndexesSFDD(14));

    t_0_xzz_zz_xz = intsBufferSFDD.data(intsIndexesSFDD(15));

    t_0_xzz_yz_xz = intsBufferSFDD.data(intsIndexesSFDD(16));

    t_0_xzz_xz_xz = intsBufferSFDD.data(intsIndexesSFDD(17));

    t_0_xzz_xx_zz = intsBufferSFDD.data(intsIndexesSFDD(18));

    t_0_xyz_zz_xy = intsBufferSFDD.data(intsIndexesSFDD(19));

    t_0_xyz_yz_xy = intsBufferSFDD.data(intsIndexesSFDD(20));

    t_0_xyz_yy_xz = intsBufferSFDD.data(intsIndexesSFDD(21));

    t_0_xyz_xz_xy = intsBufferSFDD.data(intsIndexesSFDD(22));

    t_0_xyz_xy_xz = intsBufferSFDD.data(intsIndexesSFDD(23));

    t_0_xyz_xx_yz = intsBufferSFDD.data(intsIndexesSFDD(24));

    t_0_xyy_yy_xy = intsBufferSFDD.data(intsIndexesSFDD(25));

    t_0_xyy_xy_xy = intsBufferSFDD.data(intsIndexesSFDD(26));

    t_0_xyy_xx_yy = intsBufferSFDD.data(intsIndexesSFDD(27));

    t_0_xxz_zz_xx = intsBufferSFDD.data(intsIndexesSFDD(28));

    t_0_xxz_yz_xx = intsBufferSFDD.data(intsIndexesSFDD(29));

    t_0_xxz_xz_xx = intsBufferSFDD.data(intsIndexesSFDD(30));

    t_0_xxz_xx_xz = intsBufferSFDD.data(intsIndexesSFDD(31));

    t_0_xxy_yy_xx = intsBufferSFDD.data(intsIndexesSFDD(32));

    t_0_xxy_xy_xx = intsBufferSFDD.data(intsIndexesSFDD(33));

    t_0_xxy_xx_xy = intsBufferSFDD.data(intsIndexesSFDD(34));

    t_0_xxx_xx_xx = intsBufferSFDD.data(intsIndexesSFDD(35));

    // set up (SGDD) integral components

    t_0_zzzz_zz_zz = intsBufferSGDD.data(intsIndexesSGDD(0));

    t_0_yzzz_zz_yz = intsBufferSGDD.data(intsIndexesSGDD(1));

    t_0_yzzz_yz_zz = intsBufferSGDD.data(intsIndexesSGDD(2));

    t_0_yyzz_zz_yy = intsBufferSGDD.data(intsIndexesSGDD(3));

    t_0_yyzz_yz_yz = intsBufferSGDD.data(intsIndexesSGDD(4));

    t_0_yyzz_yy_zz = intsBufferSGDD.data(intsIndexesSGDD(5));

    t_0_yyyz_yz_yy = intsBufferSGDD.data(intsIndexesSGDD(6));

    t_0_yyyz_yy_yz = intsBufferSGDD.data(intsIndexesSGDD(7));

    t_0_yyyy_yy_yy = intsBufferSGDD.data(intsIndexesSGDD(8));

    t_0_xzzz_zz_xz = intsBufferSGDD.data(intsIndexesSGDD(9));

    t_0_xzzz_xz_zz = intsBufferSGDD.data(intsIndexesSGDD(10));

    t_0_xyzz_zz_xy = intsBufferSGDD.data(intsIndexesSGDD(11));

    t_0_xyzz_yz_xz = intsBufferSGDD.data(intsIndexesSGDD(12));

    t_0_xyzz_xz_yz = intsBufferSGDD.data(intsIndexesSGDD(13));

    t_0_xyzz_xy_zz = intsBufferSGDD.data(intsIndexesSGDD(14));

    t_0_xyyz_yz_xy = intsBufferSGDD.data(intsIndexesSGDD(15));

    t_0_xyyz_yy_xz = intsBufferSGDD.data(intsIndexesSGDD(16));

    t_0_xyyz_xz_yy = intsBufferSGDD.data(intsIndexesSGDD(17));

    t_0_xyyz_xy_yz = intsBufferSGDD.data(intsIndexesSGDD(18));

    t_0_xyyy_yy_xy = intsBufferSGDD.data(intsIndexesSGDD(19));

    t_0_xyyy_xy_yy = intsBufferSGDD.data(intsIndexesSGDD(20));

    t_0_xxzz_zz_xx = intsBufferSGDD.data(intsIndexesSGDD(21));

    t_0_xxzz_xz_xz = intsBufferSGDD.data(intsIndexesSGDD(22));

    t_0_xxzz_xx_zz = intsBufferSGDD.data(intsIndexesSGDD(23));

    t_0_xxyz_yz_xx = intsBufferSGDD.data(intsIndexesSGDD(24));

    t_0_xxyz_xz_xy = intsBufferSGDD.data(intsIndexesSGDD(25));

    t_0_xxyz_xy_xz = intsBufferSGDD.data(intsIndexesSGDD(26));

    t_0_xxyz_xx_yz = intsBufferSGDD.data(intsIndexesSGDD(27));

    t_0_xxyy_yy_xx = intsBufferSGDD.data(intsIndexesSGDD(28));

    t_0_xxyy_xy_xy = intsBufferSGDD.data(intsIndexesSGDD(29));

    t_0_xxyy_xx_yy = intsBufferSGDD.data(intsIndexesSGDD(30));

    t_0_xxxz_xz_xx = intsBufferSGDD.data(intsIndexesSGDD(31));

    t_0_xxxz_xx_xz = intsBufferSGDD.data(intsIndexesSGDD(32));

    t_0_xxxy_xy_xx = intsBufferSGDD.data(intsIndexesSGDD(33));

    t_0_xxxy_xx_xy = intsBufferSGDD.data(intsIndexesSGDD(34));

    t_0_xxxx_xx_xx = intsBufferSGDD.data(intsIndexesSGDD(35));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_xxx_xx_xx, t_0_xxxx_xx_xx, t_0_xxxy_xx_xy,\
                           t_0_xxxy_xy_xx, t_0_xxxz_xx_xz, t_0_xxxz_xz_xx, t_0_xxy_xx_xy,\
                           t_0_xxy_xy_xx, t_0_xxy_yy_xx, t_0_xxyy_xx_yy, t_0_xxyy_xy_xy,\
                           t_0_xxyy_yy_xx, t_0_xxyz_xx_yz, t_0_xxyz_xy_xz, t_0_xxyz_xz_xy,\
                           t_0_xxyz_yz_xx, t_0_xxz_xx_xz, t_0_xxz_xz_xx, t_0_xxz_yz_xx,\
                           t_0_xxz_zz_xx, t_0_xxzz_xx_zz, t_0_xxzz_xz_xz, t_0_xxzz_zz_xx,\
                           t_0_xyy_xx_yy, t_0_xyy_xy_xy, t_0_xyy_yy_xy, t_0_xyyy_xy_yy,\
                           t_0_xyyy_yy_xy, t_0_xyyz_xy_yz, t_0_xyyz_xz_yy, t_0_xyyz_yy_xz,\
                           t_0_xyyz_yz_xy, t_0_xyz_xx_yz, t_0_xyz_xy_xz, t_0_xyz_xz_xy,\
                           t_0_xyz_yy_xz, t_0_xyz_yz_xy, t_0_xyz_zz_xy, t_0_xyzz_xy_zz,\
                           t_0_xyzz_xz_yz, t_0_xyzz_yz_xz, t_0_xyzz_zz_xy, t_0_xzz_xx_zz,\
                           t_0_xzz_xz_xz, t_0_xzz_yz_xz, t_0_xzz_zz_xz, t_0_xzzz_xz_zz,\
                           t_0_xzzz_zz_xz, t_0_yyy_xy_yy, t_0_yyy_yy_yy, t_0_yyyy_yy_yy,\
                           t_0_yyyz_yy_yz, t_0_yyyz_yz_yy, t_0_yyz_xy_yz, t_0_yyz_xz_yy,\
                           t_0_yyz_yy_yz, t_0_yyz_yz_yy, t_0_yyz_zz_yy, t_0_yyzz_yy_zz,\
                           t_0_yyzz_yz_yz, t_0_yyzz_zz_yy, t_0_yzz_xy_zz, t_0_yzz_xz_yz,\
                           t_0_yzz_yy_zz, t_0_yzz_yz_yz, t_0_yzz_zz_yz, t_0_yzzz_yz_zz,\
                           t_0_yzzz_zz_yz, t_0_zzz_xz_zz, t_0_zzz_yz_zz, t_0_zzz_zz_zz,\
                           t_0_zzzz_zz_zz, t_x_xxx_xx_xx, t_x_xxy_xx_xy, t_x_xxy_xy_xx,\
                           t_x_xxz_xx_xz, t_x_xxz_xz_xx, t_x_xyy_xx_yy, t_x_xyy_xy_xy,\
                           t_x_xyz_xx_yz, t_x_xyz_xy_xz, t_x_xyz_xz_xy, t_x_xzz_xx_zz,\
                           t_x_xzz_xz_xz, t_x_yyy_xy_yy, t_x_yyz_xy_yz, t_x_yyz_xz_yy,\
                           t_x_yzz_xy_zz, t_x_yzz_xz_yz, t_x_zzz_xz_zz, t_y_xxy_yy_xx,\
                           t_y_xxz_yz_xx, t_y_xyy_yy_xy, t_y_xyz_yy_xz, t_y_xyz_yz_xy,\
                           t_y_xzz_yz_xz, t_y_yyy_yy_yy, t_y_yyz_yy_yz, t_y_yyz_yz_yy,\
                           t_y_yzz_yy_zz, t_y_yzz_yz_yz, t_y_zzz_yz_zz, t_z_xxz_zz_xx,\
                           t_z_xyz_zz_xy, t_z_xzz_zz_xz, t_z_yyz_zz_yy, t_z_yzz_zz_yz,\
                           t_z_zzz_zz_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_zz_zz[i] = t_0_zzzz_zz_zz[i] - rab_z[i] * t_0_zzz_zz_zz[i];

        t_z_yzz_zz_yz[i] = t_0_yzzz_zz_yz[i] - rab_z[i] * t_0_yzz_zz_yz[i];

        t_z_yyz_zz_yy[i] = t_0_yyzz_zz_yy[i] - rab_z[i] * t_0_yyz_zz_yy[i];

        t_z_xzz_zz_xz[i] = t_0_xzzz_zz_xz[i] - rab_z[i] * t_0_xzz_zz_xz[i];

        t_z_xyz_zz_xy[i] = t_0_xyzz_zz_xy[i] - rab_z[i] * t_0_xyz_zz_xy[i];

        t_z_xxz_zz_xx[i] = t_0_xxzz_zz_xx[i] - rab_z[i] * t_0_xxz_zz_xx[i];

        t_y_zzz_yz_zz[i] = t_0_yzzz_yz_zz[i] - rab_y[i] * t_0_zzz_yz_zz[i];

        t_y_yzz_yz_yz[i] = t_0_yyzz_yz_yz[i] - rab_y[i] * t_0_yzz_yz_yz[i];

        t_y_yzz_yy_zz[i] = t_0_yyzz_yy_zz[i] - rab_y[i] * t_0_yzz_yy_zz[i];

        t_y_yyz_yz_yy[i] = t_0_yyyz_yz_yy[i] - rab_y[i] * t_0_yyz_yz_yy[i];

        t_y_yyz_yy_yz[i] = t_0_yyyz_yy_yz[i] - rab_y[i] * t_0_yyz_yy_yz[i];

        t_y_yyy_yy_yy[i] = t_0_yyyy_yy_yy[i] - rab_y[i] * t_0_yyy_yy_yy[i];

        t_y_xzz_yz_xz[i] = t_0_xyzz_yz_xz[i] - rab_y[i] * t_0_xzz_yz_xz[i];

        t_y_xyz_yz_xy[i] = t_0_xyyz_yz_xy[i] - rab_y[i] * t_0_xyz_yz_xy[i];

        t_y_xyz_yy_xz[i] = t_0_xyyz_yy_xz[i] - rab_y[i] * t_0_xyz_yy_xz[i];

        t_y_xyy_yy_xy[i] = t_0_xyyy_yy_xy[i] - rab_y[i] * t_0_xyy_yy_xy[i];

        t_y_xxz_yz_xx[i] = t_0_xxyz_yz_xx[i] - rab_y[i] * t_0_xxz_yz_xx[i];

        t_y_xxy_yy_xx[i] = t_0_xxyy_yy_xx[i] - rab_y[i] * t_0_xxy_yy_xx[i];

        t_x_zzz_xz_zz[i] = t_0_xzzz_xz_zz[i] - rab_x[i] * t_0_zzz_xz_zz[i];

        t_x_yzz_xz_yz[i] = t_0_xyzz_xz_yz[i] - rab_x[i] * t_0_yzz_xz_yz[i];

        t_x_yzz_xy_zz[i] = t_0_xyzz_xy_zz[i] - rab_x[i] * t_0_yzz_xy_zz[i];

        t_x_yyz_xz_yy[i] = t_0_xyyz_xz_yy[i] - rab_x[i] * t_0_yyz_xz_yy[i];

        t_x_yyz_xy_yz[i] = t_0_xyyz_xy_yz[i] - rab_x[i] * t_0_yyz_xy_yz[i];

        t_x_yyy_xy_yy[i] = t_0_xyyy_xy_yy[i] - rab_x[i] * t_0_yyy_xy_yy[i];

        t_x_xzz_xz_xz[i] = t_0_xxzz_xz_xz[i] - rab_x[i] * t_0_xzz_xz_xz[i];

        t_x_xzz_xx_zz[i] = t_0_xxzz_xx_zz[i] - rab_x[i] * t_0_xzz_xx_zz[i];

        t_x_xyz_xz_xy[i] = t_0_xxyz_xz_xy[i] - rab_x[i] * t_0_xyz_xz_xy[i];

        t_x_xyz_xy_xz[i] = t_0_xxyz_xy_xz[i] - rab_x[i] * t_0_xyz_xy_xz[i];

        t_x_xyz_xx_yz[i] = t_0_xxyz_xx_yz[i] - rab_x[i] * t_0_xyz_xx_yz[i];

        t_x_xyy_xy_xy[i] = t_0_xxyy_xy_xy[i] - rab_x[i] * t_0_xyy_xy_xy[i];

        t_x_xyy_xx_yy[i] = t_0_xxyy_xx_yy[i] - rab_x[i] * t_0_xyy_xx_yy[i];

        t_x_xxz_xz_xx[i] = t_0_xxxz_xz_xx[i] - rab_x[i] * t_0_xxz_xz_xx[i];

        t_x_xxz_xx_xz[i] = t_0_xxxz_xx_xz[i] - rab_x[i] * t_0_xxz_xx_xz[i];

        t_x_xxy_xy_xx[i] = t_0_xxxy_xy_xx[i] - rab_x[i] * t_0_xxy_xy_xx[i];

        t_x_xxy_xx_xy[i] = t_0_xxxy_xx_xy[i] - rab_x[i] * t_0_xxy_xx_xy[i];

        t_x_xxx_xx_xx[i] = t_0_xxxx_xx_xx[i] - rab_x[i] * t_0_xxx_xx_xx[i];
    }
}


} // derirec namespace
