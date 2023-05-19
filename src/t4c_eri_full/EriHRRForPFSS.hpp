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
compHostHRRForPFSS_V0(      BufferHostXY<T>&      intsBufferPFSS,
                      const BufferHostX<int32_t>& intsIndexesPFSS,
                      const BufferHostXY<T>&      intsBufferSFSS,
                      const BufferHostX<int32_t>& intsIndexesSFSS,
                      const BufferHostXY<T>&      intsBufferSGSS,
                      const BufferHostX<int32_t>& intsIndexesSGSS,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SFSS) integral components

    t_0_zzz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(0));

    t_0_yzz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(1));

    t_0_yyz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(2));

    t_0_yyy_0_0 = intsBufferSFSS.data(intsIndexesSFSS(3));

    t_0_xzz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(4));

    t_0_xyz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(5));

    t_0_xyy_0_0 = intsBufferSFSS.data(intsIndexesSFSS(6));

    t_0_xxz_0_0 = intsBufferSFSS.data(intsIndexesSFSS(7));

    t_0_xxy_0_0 = intsBufferSFSS.data(intsIndexesSFSS(8));

    t_0_xxx_0_0 = intsBufferSFSS.data(intsIndexesSFSS(9));

    // set up (SGSS) integral components

    t_0_zzzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(0));

    t_0_yzzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(1));

    t_0_yyzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(2));

    t_0_yyyz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(3));

    t_0_yyyy_0_0 = intsBufferSGSS.data(intsIndexesSGSS(4));

    t_0_xzzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(5));

    t_0_xyzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(6));

    t_0_xyyz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(7));

    t_0_xyyy_0_0 = intsBufferSGSS.data(intsIndexesSGSS(8));

    t_0_xxzz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(9));

    t_0_xxyz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(10));

    t_0_xxyy_0_0 = intsBufferSGSS.data(intsIndexesSGSS(11));

    t_0_xxxz_0_0 = intsBufferSGSS.data(intsIndexesSGSS(12));

    t_0_xxxy_0_0 = intsBufferSGSS.data(intsIndexesSGSS(13));

    t_0_xxxx_0_0 = intsBufferSGSS.data(intsIndexesSGSS(14));

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_xxx_0_0, t_0_xxxx_0_0, t_0_xxxy_0_0,\
                           t_0_xxxz_0_0, t_0_xxy_0_0, t_0_xxyy_0_0, t_0_xxyz_0_0, t_0_xxz_0_0,\
                           t_0_xxzz_0_0, t_0_xyy_0_0, t_0_xyyy_0_0, t_0_xyyz_0_0, t_0_xyz_0_0,\
                           t_0_xyzz_0_0, t_0_xzz_0_0, t_0_xzzz_0_0, t_0_yyy_0_0, t_0_yyyy_0_0,\
                           t_0_yyyz_0_0, t_0_yyz_0_0, t_0_yyzz_0_0, t_0_yzz_0_0, t_0_yzzz_0_0,\
                           t_0_zzz_0_0, t_0_zzzz_0_0, t_x_xxx_0_0, t_x_xxy_0_0, t_x_xxz_0_0,\
                           t_x_xyy_0_0, t_x_xyz_0_0, t_x_xzz_0_0, t_x_yyy_0_0, t_x_yyz_0_0,\
                           t_x_yzz_0_0, t_x_zzz_0_0, t_y_xxy_0_0, t_y_xxz_0_0, t_y_xyy_0_0,\
                           t_y_xyz_0_0, t_y_xzz_0_0, t_y_yyy_0_0, t_y_yyz_0_0, t_y_yzz_0_0,\
                           t_y_zzz_0_0, t_z_xxz_0_0, t_z_xyz_0_0, t_z_xzz_0_0, t_z_yyz_0_0,\
                           t_z_yzz_0_0, t_z_zzz_0_0 : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zzz_0_0[i] = t_0_zzzz_0_0[i] - rab_z[i] * t_0_zzz_0_0[i];

        t_z_yzz_0_0[i] = t_0_yzzz_0_0[i] - rab_z[i] * t_0_yzz_0_0[i];

        t_z_yyz_0_0[i] = t_0_yyzz_0_0[i] - rab_z[i] * t_0_yyz_0_0[i];

        t_z_xzz_0_0[i] = t_0_xzzz_0_0[i] - rab_z[i] * t_0_xzz_0_0[i];

        t_z_xyz_0_0[i] = t_0_xyzz_0_0[i] - rab_z[i] * t_0_xyz_0_0[i];

        t_z_xxz_0_0[i] = t_0_xxzz_0_0[i] - rab_z[i] * t_0_xxz_0_0[i];

        t_y_zzz_0_0[i] = t_0_yzzz_0_0[i] - rab_y[i] * t_0_zzz_0_0[i];

        t_y_yzz_0_0[i] = t_0_yyzz_0_0[i] - rab_y[i] * t_0_yzz_0_0[i];

        t_y_yyz_0_0[i] = t_0_yyyz_0_0[i] - rab_y[i] * t_0_yyz_0_0[i];

        t_y_yyy_0_0[i] = t_0_yyyy_0_0[i] - rab_y[i] * t_0_yyy_0_0[i];

        t_y_xzz_0_0[i] = t_0_xyzz_0_0[i] - rab_y[i] * t_0_xzz_0_0[i];

        t_y_xyz_0_0[i] = t_0_xyyz_0_0[i] - rab_y[i] * t_0_xyz_0_0[i];

        t_y_xyy_0_0[i] = t_0_xyyy_0_0[i] - rab_y[i] * t_0_xyy_0_0[i];

        t_y_xxz_0_0[i] = t_0_xxyz_0_0[i] - rab_y[i] * t_0_xxz_0_0[i];

        t_y_xxy_0_0[i] = t_0_xxyy_0_0[i] - rab_y[i] * t_0_xxy_0_0[i];

        t_x_zzz_0_0[i] = t_0_xzzz_0_0[i] - rab_x[i] * t_0_zzz_0_0[i];

        t_x_yzz_0_0[i] = t_0_xyzz_0_0[i] - rab_x[i] * t_0_yzz_0_0[i];

        t_x_yyz_0_0[i] = t_0_xyyz_0_0[i] - rab_x[i] * t_0_yyz_0_0[i];

        t_x_yyy_0_0[i] = t_0_xyyy_0_0[i] - rab_x[i] * t_0_yyy_0_0[i];

        t_x_xzz_0_0[i] = t_0_xxzz_0_0[i] - rab_x[i] * t_0_xzz_0_0[i];

        t_x_xyz_0_0[i] = t_0_xxyz_0_0[i] - rab_x[i] * t_0_xyz_0_0[i];

        t_x_xyy_0_0[i] = t_0_xxyy_0_0[i] - rab_x[i] * t_0_xyy_0_0[i];

        t_x_xxz_0_0[i] = t_0_xxxz_0_0[i] - rab_x[i] * t_0_xxz_0_0[i];

        t_x_xxy_0_0[i] = t_0_xxxy_0_0[i] - rab_x[i] * t_0_xxy_0_0[i];

        t_x_xxx_0_0[i] = t_0_xxxx_0_0[i] - rab_x[i] * t_0_xxx_0_0[i];
    }
}


} // derirec namespace
