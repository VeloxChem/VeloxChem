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
compHostHRRForPDSS_V0(      BufferHostXY<T>&      intsBufferPDSS,
                      const BufferHostX<int32_t>& intsIndexesPDSS,
                      const BufferHostXY<T>&      intsBufferSDSS,
                      const BufferHostX<int32_t>& intsIndexesSDSS,
                      const BufferHostXY<T>&      intsBufferSFSS,
                      const BufferHostX<int32_t>& intsIndexesSFSS,
                      const BufferHostMY<T, 3>&   rDistancesAB,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(AB) distances

    auto rab_z = rDistancesAB.data(2);

    auto rab_y = rDistancesAB.data(1);

    auto rab_x = rDistancesAB.data(0);

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

    // set up (SDSS) integral components

    t_0_zz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(0));

    t_0_yz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(1));

    t_0_yy_0_0 = intsBufferSDSS.data(intsIndexesSDSS(2));

    t_0_xz_0_0 = intsBufferSDSS.data(intsIndexesSDSS(3));

    t_0_xy_0_0 = intsBufferSDSS.data(intsIndexesSDSS(4));

    t_0_xx_0_0 = intsBufferSDSS.data(intsIndexesSDSS(5));

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

    #pragma omp simd align(rab_x, rab_y, rab_z, t_0_xx_0_0, t_0_xxx_0_0, t_0_xxy_0_0,\
                           t_0_xxz_0_0, t_0_xy_0_0, t_0_xyy_0_0, t_0_xyz_0_0, t_0_xz_0_0,\
                           t_0_xzz_0_0, t_0_yy_0_0, t_0_yyy_0_0, t_0_yyz_0_0, t_0_yz_0_0,\
                           t_0_yzz_0_0, t_0_zz_0_0, t_0_zzz_0_0, t_x_xx_0_0, t_x_xy_0_0,\
                           t_x_xz_0_0, t_x_yy_0_0, t_x_yz_0_0, t_x_zz_0_0, t_y_xx_0_0,\
                           t_y_xy_0_0, t_y_xz_0_0, t_y_yy_0_0, t_y_yz_0_0, t_y_zz_0_0,\
                           t_z_xx_0_0, t_z_xy_0_0, t_z_xz_0_0, t_z_yy_0_0, t_z_yz_0_0,\
                           t_z_zz_0_0 : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_z_zz_0_0[i] = t_0_zzz_0_0[i] - rab_z[i] * t_0_zz_0_0[i];

        t_z_yz_0_0[i] = t_0_yzz_0_0[i] - rab_z[i] * t_0_yz_0_0[i];

        t_z_yy_0_0[i] = t_0_yyz_0_0[i] - rab_z[i] * t_0_yy_0_0[i];

        t_z_xz_0_0[i] = t_0_xzz_0_0[i] - rab_z[i] * t_0_xz_0_0[i];

        t_z_xy_0_0[i] = t_0_xyz_0_0[i] - rab_z[i] * t_0_xy_0_0[i];

        t_z_xx_0_0[i] = t_0_xxz_0_0[i] - rab_z[i] * t_0_xx_0_0[i];

        t_y_zz_0_0[i] = t_0_yzz_0_0[i] - rab_y[i] * t_0_zz_0_0[i];

        t_y_yz_0_0[i] = t_0_yyz_0_0[i] - rab_y[i] * t_0_yz_0_0[i];

        t_y_yy_0_0[i] = t_0_yyy_0_0[i] - rab_y[i] * t_0_yy_0_0[i];

        t_y_xz_0_0[i] = t_0_xyz_0_0[i] - rab_y[i] * t_0_xz_0_0[i];

        t_y_xy_0_0[i] = t_0_xyy_0_0[i] - rab_y[i] * t_0_xy_0_0[i];

        t_y_xx_0_0[i] = t_0_xxy_0_0[i] - rab_y[i] * t_0_xx_0_0[i];

        t_x_zz_0_0[i] = t_0_xzz_0_0[i] - rab_x[i] * t_0_zz_0_0[i];

        t_x_yz_0_0[i] = t_0_xyz_0_0[i] - rab_x[i] * t_0_yz_0_0[i];

        t_x_yy_0_0[i] = t_0_xyy_0_0[i] - rab_x[i] * t_0_yy_0_0[i];

        t_x_xz_0_0[i] = t_0_xxz_0_0[i] - rab_x[i] * t_0_xz_0_0[i];

        t_x_xy_0_0[i] = t_0_xxy_0_0[i] - rab_x[i] * t_0_xy_0_0[i];

        t_x_xx_0_0[i] = t_0_xxx_0_0[i] - rab_x[i] * t_0_xx_0_0[i];
    }
}


} // derirec namespace
