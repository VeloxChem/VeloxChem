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
compHostHRRForSSPD_V0(      BufferHostXY<T>&      intsBufferSSPD,
                      const BufferHostX<int32_t>& intsIndexesSSPD,
                      const BufferHostXY<T>&      intsBufferSSSD,
                      const BufferHostX<int32_t>& intsIndexesSSSD,
                      const BufferHostXY<T>&      intsBufferSSSF,
                      const BufferHostX<int32_t>& intsIndexesSSSF,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SSPD) integral components

    t_0_0_z_zz = intsBufferSSPD.data(intsIndexesSSPD(0));

    t_0_0_z_yz = intsBufferSSPD.data(intsIndexesSSPD(1));

    t_0_0_z_yy = intsBufferSSPD.data(intsIndexesSSPD(2));

    t_0_0_z_xz = intsBufferSSPD.data(intsIndexesSSPD(3));

    t_0_0_z_xy = intsBufferSSPD.data(intsIndexesSSPD(4));

    t_0_0_z_xx = intsBufferSSPD.data(intsIndexesSSPD(5));

    t_0_0_y_zz = intsBufferSSPD.data(intsIndexesSSPD(6));

    t_0_0_y_yz = intsBufferSSPD.data(intsIndexesSSPD(7));

    t_0_0_y_yy = intsBufferSSPD.data(intsIndexesSSPD(8));

    t_0_0_y_xz = intsBufferSSPD.data(intsIndexesSSPD(9));

    t_0_0_y_xy = intsBufferSSPD.data(intsIndexesSSPD(10));

    t_0_0_y_xx = intsBufferSSPD.data(intsIndexesSSPD(11));

    t_0_0_x_zz = intsBufferSSPD.data(intsIndexesSSPD(12));

    t_0_0_x_yz = intsBufferSSPD.data(intsIndexesSSPD(13));

    t_0_0_x_yy = intsBufferSSPD.data(intsIndexesSSPD(14));

    t_0_0_x_xz = intsBufferSSPD.data(intsIndexesSSPD(15));

    t_0_0_x_xy = intsBufferSSPD.data(intsIndexesSSPD(16));

    t_0_0_x_xx = intsBufferSSPD.data(intsIndexesSSPD(17));

    // set up (SSSD) integral components

    t_0_0_0_zz = intsBufferSSSD.data(intsIndexesSSSD(0));

    t_0_0_0_yz = intsBufferSSSD.data(intsIndexesSSSD(1));

    t_0_0_0_yy = intsBufferSSSD.data(intsIndexesSSSD(2));

    t_0_0_0_xz = intsBufferSSSD.data(intsIndexesSSSD(3));

    t_0_0_0_xy = intsBufferSSSD.data(intsIndexesSSSD(4));

    t_0_0_0_xx = intsBufferSSSD.data(intsIndexesSSSD(5));

    // set up (SSSF) integral components

    t_0_0_0_zzz = intsBufferSSSF.data(intsIndexesSSSF(0));

    t_0_0_0_yzz = intsBufferSSSF.data(intsIndexesSSSF(1));

    t_0_0_0_yyz = intsBufferSSSF.data(intsIndexesSSSF(2));

    t_0_0_0_yyy = intsBufferSSSF.data(intsIndexesSSSF(3));

    t_0_0_0_xzz = intsBufferSSSF.data(intsIndexesSSSF(4));

    t_0_0_0_xyz = intsBufferSSSF.data(intsIndexesSSSF(5));

    t_0_0_0_xyy = intsBufferSSSF.data(intsIndexesSSSF(6));

    t_0_0_0_xxz = intsBufferSSSF.data(intsIndexesSSSF(7));

    t_0_0_0_xxy = intsBufferSSSF.data(intsIndexesSSSF(8));

    t_0_0_0_xxx = intsBufferSSSF.data(intsIndexesSSSF(9));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_0_0_xx, t_0_0_0_xxx, t_0_0_0_xxy,\
                           t_0_0_0_xxz, t_0_0_0_xy, t_0_0_0_xyy, t_0_0_0_xyz, t_0_0_0_xz,\
                           t_0_0_0_xzz, t_0_0_0_yy, t_0_0_0_yyy, t_0_0_0_yyz, t_0_0_0_yz,\
                           t_0_0_0_yzz, t_0_0_0_zz, t_0_0_0_zzz, t_0_0_x_xx, t_0_0_x_xy,\
                           t_0_0_x_xz, t_0_0_x_yy, t_0_0_x_yz, t_0_0_x_zz, t_0_0_y_xx,\
                           t_0_0_y_xy, t_0_0_y_xz, t_0_0_y_yy, t_0_0_y_yz, t_0_0_y_zz,\
                           t_0_0_z_xx, t_0_0_z_xy, t_0_0_z_xz, t_0_0_z_yy, t_0_0_z_yz,\
                           t_0_0_z_zz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_0_z_zz[i] = t_0_0_0_zzz[i] - rcd_z[i] * t_0_0_0_zz[i];

        t_0_0_z_yz[i] = t_0_0_0_yzz[i] - rcd_z[i] * t_0_0_0_yz[i];

        t_0_0_z_yy[i] = t_0_0_0_yyz[i] - rcd_z[i] * t_0_0_0_yy[i];

        t_0_0_z_xz[i] = t_0_0_0_xzz[i] - rcd_z[i] * t_0_0_0_xz[i];

        t_0_0_z_xy[i] = t_0_0_0_xyz[i] - rcd_z[i] * t_0_0_0_xy[i];

        t_0_0_z_xx[i] = t_0_0_0_xxz[i] - rcd_z[i] * t_0_0_0_xx[i];

        t_0_0_y_zz[i] = t_0_0_0_yzz[i] - rcd_y[i] * t_0_0_0_zz[i];

        t_0_0_y_yz[i] = t_0_0_0_yyz[i] - rcd_y[i] * t_0_0_0_yz[i];

        t_0_0_y_yy[i] = t_0_0_0_yyy[i] - rcd_y[i] * t_0_0_0_yy[i];

        t_0_0_y_xz[i] = t_0_0_0_xyz[i] - rcd_y[i] * t_0_0_0_xz[i];

        t_0_0_y_xy[i] = t_0_0_0_xyy[i] - rcd_y[i] * t_0_0_0_xy[i];

        t_0_0_y_xx[i] = t_0_0_0_xxy[i] - rcd_y[i] * t_0_0_0_xx[i];

        t_0_0_x_zz[i] = t_0_0_0_xzz[i] - rcd_x[i] * t_0_0_0_zz[i];

        t_0_0_x_yz[i] = t_0_0_0_xyz[i] - rcd_x[i] * t_0_0_0_yz[i];

        t_0_0_x_yy[i] = t_0_0_0_xyy[i] - rcd_x[i] * t_0_0_0_yy[i];

        t_0_0_x_xz[i] = t_0_0_0_xxz[i] - rcd_x[i] * t_0_0_0_xz[i];

        t_0_0_x_xy[i] = t_0_0_0_xxy[i] - rcd_x[i] * t_0_0_0_xy[i];

        t_0_0_x_xx[i] = t_0_0_0_xxx[i] - rcd_x[i] * t_0_0_0_xx[i];
    }
}


} // derirec namespace
