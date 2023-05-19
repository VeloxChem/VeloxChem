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
compHostHRRForDDSS_V0(      BufferHostXY<T>&      intsBufferDDSS,
                      const BufferHostX<int32_t>& intsIndexesDDSS,
                      const BufferHostXY<T>&      intsBufferPDSS,
                      const BufferHostX<int32_t>& intsIndexesPDSS,
                      const BufferHostXY<T>&      intsBufferPFSS,
                      const BufferHostX<int32_t>& intsIndexesPFSS,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

    // set up (DDSS) integral components

    t_zz_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(0));

    t_zz_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(1));

    t_zz_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(2));

    t_zz_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(3));

    t_zz_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(4));

    t_zz_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(5));

    t_yz_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(6));

    t_yz_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(7));

    t_yz_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(8));

    t_yz_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(9));

    t_yz_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(10));

    t_yz_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(11));

    t_yy_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(12));

    t_yy_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(13));

    t_yy_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(14));

    t_yy_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(15));

    t_yy_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(16));

    t_yy_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(17));

    t_xz_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(18));

    t_xz_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(19));

    t_xz_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(20));

    t_xz_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(21));

    t_xz_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(22));

    t_xz_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(23));

    t_xy_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(24));

    t_xy_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(25));

    t_xy_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(26));

    t_xy_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(27));

    t_xy_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(28));

    t_xy_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(29));

    t_xx_zz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(30));

    t_xx_yz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(31));

    t_xx_yy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(32));

    t_xx_xz_0_0 = intsBufferDDSS.data(intsIndexesDDSS(33));

    t_xx_xy_0_0 = intsBufferDDSS.data(intsIndexesDDSS(34));

    t_xx_xx_0_0 = intsBufferDDSS.data(intsIndexesDDSS(35));

    // set up (PDSS) integral components

    t_z_zz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(0));

    t_z_yz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(1));

    t_z_yy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(2));

    t_z_xz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(3));

    t_z_xy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(4));

    t_z_xx_0_0 = intsBufferPDSS.data(intsIndexesPDSS(5));

    t_y_zz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(6));

    t_y_yz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(7));

    t_y_yy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(8));

    t_y_xz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(9));

    t_y_xy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(10));

    t_y_xx_0_0 = intsBufferPDSS.data(intsIndexesPDSS(11));

    t_x_zz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(12));

    t_x_yz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(13));

    t_x_yy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(14));

    t_x_xz_0_0 = intsBufferPDSS.data(intsIndexesPDSS(15));

    t_x_xy_0_0 = intsBufferPDSS.data(intsIndexesPDSS(16));

    t_x_xx_0_0 = intsBufferPDSS.data(intsIndexesPDSS(17));

    // set up (PFSS) integral components

    t_z_zzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(0));

    t_z_yzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(1));

    t_z_yyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(2));

    t_z_xzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(3));

    t_z_xyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(4));

    t_z_xxz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(5));

    t_y_zzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(6));

    t_y_yzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(7));

    t_y_yyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(8));

    t_y_yyy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(9));

    t_y_xzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(10));

    t_y_xyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(11));

    t_y_xyy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(12));

    t_y_xxz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(13));

    t_y_xxy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(14));

    t_x_zzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(15));

    t_x_yzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(16));

    t_x_yyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(17));

    t_x_yyy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(18));

    t_x_xzz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(19));

    t_x_xyz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(20));

    t_x_xyy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(21));

    t_x_xxz_0_0 = intsBufferPFSS.data(intsIndexesPFSS(22));

    t_x_xxy_0_0 = intsBufferPFSS.data(intsIndexesPFSS(23));

    t_x_xxx_0_0 = intsBufferPFSS.data(intsIndexesPFSS(24));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_x_xx_0_0, t_x_xxx_0_0, t_x_xxy_0_0,\
                           t_x_xxz_0_0, t_x_xy_0_0, t_x_xyy_0_0, t_x_xyz_0_0, t_x_xz_0_0,\
                           t_x_xzz_0_0, t_x_yy_0_0, t_x_yyy_0_0, t_x_yyz_0_0, t_x_yz_0_0,\
                           t_x_yzz_0_0, t_x_zz_0_0, t_x_zzz_0_0, t_xx_xx_0_0, t_xx_xy_0_0,\
                           t_xx_xz_0_0, t_xx_yy_0_0, t_xx_yz_0_0, t_xx_zz_0_0, t_xy_xx_0_0,\
                           t_xy_xy_0_0, t_xy_xz_0_0, t_xy_yy_0_0, t_xy_yz_0_0, t_xy_zz_0_0,\
                           t_xz_xx_0_0, t_xz_xy_0_0, t_xz_xz_0_0, t_xz_yy_0_0, t_xz_yz_0_0,\
                           t_xz_zz_0_0, t_y_xx_0_0, t_y_xxy_0_0, t_y_xxz_0_0, t_y_xy_0_0,\
                           t_y_xyy_0_0, t_y_xyz_0_0, t_y_xz_0_0, t_y_xzz_0_0, t_y_yy_0_0,\
                           t_y_yyy_0_0, t_y_yyz_0_0, t_y_yz_0_0, t_y_yzz_0_0, t_y_zz_0_0,\
                           t_y_zzz_0_0, t_yy_xx_0_0, t_yy_xy_0_0, t_yy_xz_0_0, t_yy_yy_0_0,\
                           t_yy_yz_0_0, t_yy_zz_0_0, t_yz_xx_0_0, t_yz_xy_0_0, t_yz_xz_0_0,\
                           t_yz_yy_0_0, t_yz_yz_0_0, t_yz_zz_0_0, t_z_xx_0_0, t_z_xxz_0_0,\
                           t_z_xy_0_0, t_z_xyz_0_0, t_z_xz_0_0, t_z_xzz_0_0, t_z_yy_0_0,\
                           t_z_yyz_0_0, t_z_yz_0_0, t_z_yzz_0_0, t_z_zz_0_0, t_z_zzz_0_0,\
                           t_zz_xx_0_0, t_zz_xy_0_0, t_zz_xz_0_0, t_zz_yy_0_0, t_zz_yz_0_0,\
                           t_zz_zz_0_0 : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_zz_zz_0_0[i] = t_z_zzz_0_0[i] - rab_z[i] * t_z_zz_0_0[i];

        t_zz_yz_0_0[i] = t_z_yzz_0_0[i] - rab_z[i] * t_z_yz_0_0[i];

        t_zz_yy_0_0[i] = t_z_yyz_0_0[i] - rab_z[i] * t_z_yy_0_0[i];

        t_zz_xz_0_0[i] = t_z_xzz_0_0[i] - rab_z[i] * t_z_xz_0_0[i];

        t_zz_xy_0_0[i] = t_z_xyz_0_0[i] - rab_z[i] * t_z_xy_0_0[i];

        t_zz_xx_0_0[i] = t_z_xxz_0_0[i] - rab_z[i] * t_z_xx_0_0[i];

        t_yz_zz_0_0[i] = t_y_zzz_0_0[i] - rab_z[i] * t_y_zz_0_0[i];

        t_yz_yz_0_0[i] = t_y_yzz_0_0[i] - rab_z[i] * t_y_yz_0_0[i];

        t_yz_yy_0_0[i] = t_y_yyz_0_0[i] - rab_z[i] * t_y_yy_0_0[i];

        t_yz_xz_0_0[i] = t_y_xzz_0_0[i] - rab_z[i] * t_y_xz_0_0[i];

        t_yz_xy_0_0[i] = t_y_xyz_0_0[i] - rab_z[i] * t_y_xy_0_0[i];

        t_yz_xx_0_0[i] = t_y_xxz_0_0[i] - rab_z[i] * t_y_xx_0_0[i];

        t_yy_zz_0_0[i] = t_y_yzz_0_0[i] - rab_y[i] * t_y_zz_0_0[i];

        t_yy_yz_0_0[i] = t_y_yyz_0_0[i] - rab_y[i] * t_y_yz_0_0[i];

        t_yy_yy_0_0[i] = t_y_yyy_0_0[i] - rab_y[i] * t_y_yy_0_0[i];

        t_yy_xz_0_0[i] = t_y_xyz_0_0[i] - rab_y[i] * t_y_xz_0_0[i];

        t_yy_xy_0_0[i] = t_y_xyy_0_0[i] - rab_y[i] * t_y_xy_0_0[i];

        t_yy_xx_0_0[i] = t_y_xxy_0_0[i] - rab_y[i] * t_y_xx_0_0[i];

        t_xz_zz_0_0[i] = t_x_zzz_0_0[i] - rab_z[i] * t_x_zz_0_0[i];

        t_xz_yz_0_0[i] = t_x_yzz_0_0[i] - rab_z[i] * t_x_yz_0_0[i];

        t_xz_yy_0_0[i] = t_x_yyz_0_0[i] - rab_z[i] * t_x_yy_0_0[i];

        t_xz_xz_0_0[i] = t_x_xzz_0_0[i] - rab_z[i] * t_x_xz_0_0[i];

        t_xz_xy_0_0[i] = t_x_xyz_0_0[i] - rab_z[i] * t_x_xy_0_0[i];

        t_xz_xx_0_0[i] = t_x_xxz_0_0[i] - rab_z[i] * t_x_xx_0_0[i];

        t_xy_zz_0_0[i] = t_x_yzz_0_0[i] - rab_y[i] * t_x_zz_0_0[i];

        t_xy_yz_0_0[i] = t_x_yyz_0_0[i] - rab_y[i] * t_x_yz_0_0[i];

        t_xy_yy_0_0[i] = t_x_yyy_0_0[i] - rab_y[i] * t_x_yy_0_0[i];

        t_xy_xz_0_0[i] = t_x_xyz_0_0[i] - rab_y[i] * t_x_xz_0_0[i];

        t_xy_xy_0_0[i] = t_x_xyy_0_0[i] - rab_y[i] * t_x_xy_0_0[i];

        t_xy_xx_0_0[i] = t_x_xxy_0_0[i] - rab_y[i] * t_x_xx_0_0[i];

        t_xx_zz_0_0[i] = t_x_xzz_0_0[i] - rab_x[i] * t_x_zz_0_0[i];

        t_xx_yz_0_0[i] = t_x_xyz_0_0[i] - rab_x[i] * t_x_yz_0_0[i];

        t_xx_yy_0_0[i] = t_x_xyy_0_0[i] - rab_x[i] * t_x_yy_0_0[i];

        t_xx_xz_0_0[i] = t_x_xxz_0_0[i] - rab_x[i] * t_x_xz_0_0[i];

        t_xx_xy_0_0[i] = t_x_xxy_0_0[i] - rab_x[i] * t_x_xy_0_0[i];

        t_xx_xx_0_0[i] = t_x_xxx_0_0[i] - rab_x[i] * t_x_xx_0_0[i];
    }
}


} // derirec namespace
