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
compHostHRRForSSPF_V0(      BufferHostXY<T>&      intsBufferSSPF,
                      const BufferHostX<int32_t>& intsIndexesSSPF,
                      const BufferHostXY<T>&      intsBufferSSSF,
                      const BufferHostX<int32_t>& intsIndexesSSSF,
                      const BufferHostXY<T>&      intsBufferSSSG,
                      const BufferHostX<int32_t>& intsIndexesSSSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SSPF) integral components

    t_0_0_z_zzz = intsBufferSSPF.data(intsIndexesSSPF(0));

    t_0_0_z_yzz = intsBufferSSPF.data(intsIndexesSSPF(1));

    t_0_0_z_yyz = intsBufferSSPF.data(intsIndexesSSPF(2));

    t_0_0_z_xzz = intsBufferSSPF.data(intsIndexesSSPF(3));

    t_0_0_z_xyz = intsBufferSSPF.data(intsIndexesSSPF(4));

    t_0_0_z_xxz = intsBufferSSPF.data(intsIndexesSSPF(5));

    t_0_0_y_zzz = intsBufferSSPF.data(intsIndexesSSPF(6));

    t_0_0_y_yzz = intsBufferSSPF.data(intsIndexesSSPF(7));

    t_0_0_y_yyz = intsBufferSSPF.data(intsIndexesSSPF(8));

    t_0_0_y_yyy = intsBufferSSPF.data(intsIndexesSSPF(9));

    t_0_0_y_xzz = intsBufferSSPF.data(intsIndexesSSPF(10));

    t_0_0_y_xyz = intsBufferSSPF.data(intsIndexesSSPF(11));

    t_0_0_y_xyy = intsBufferSSPF.data(intsIndexesSSPF(12));

    t_0_0_y_xxz = intsBufferSSPF.data(intsIndexesSSPF(13));

    t_0_0_y_xxy = intsBufferSSPF.data(intsIndexesSSPF(14));

    t_0_0_x_zzz = intsBufferSSPF.data(intsIndexesSSPF(15));

    t_0_0_x_yzz = intsBufferSSPF.data(intsIndexesSSPF(16));

    t_0_0_x_yyz = intsBufferSSPF.data(intsIndexesSSPF(17));

    t_0_0_x_yyy = intsBufferSSPF.data(intsIndexesSSPF(18));

    t_0_0_x_xzz = intsBufferSSPF.data(intsIndexesSSPF(19));

    t_0_0_x_xyz = intsBufferSSPF.data(intsIndexesSSPF(20));

    t_0_0_x_xyy = intsBufferSSPF.data(intsIndexesSSPF(21));

    t_0_0_x_xxz = intsBufferSSPF.data(intsIndexesSSPF(22));

    t_0_0_x_xxy = intsBufferSSPF.data(intsIndexesSSPF(23));

    t_0_0_x_xxx = intsBufferSSPF.data(intsIndexesSSPF(24));

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

    // set up (SSSG) integral components

    t_0_0_0_zzzz = intsBufferSSSG.data(intsIndexesSSSG(0));

    t_0_0_0_yzzz = intsBufferSSSG.data(intsIndexesSSSG(1));

    t_0_0_0_yyzz = intsBufferSSSG.data(intsIndexesSSSG(2));

    t_0_0_0_yyyz = intsBufferSSSG.data(intsIndexesSSSG(3));

    t_0_0_0_yyyy = intsBufferSSSG.data(intsIndexesSSSG(4));

    t_0_0_0_xzzz = intsBufferSSSG.data(intsIndexesSSSG(5));

    t_0_0_0_xyzz = intsBufferSSSG.data(intsIndexesSSSG(6));

    t_0_0_0_xyyz = intsBufferSSSG.data(intsIndexesSSSG(7));

    t_0_0_0_xyyy = intsBufferSSSG.data(intsIndexesSSSG(8));

    t_0_0_0_xxzz = intsBufferSSSG.data(intsIndexesSSSG(9));

    t_0_0_0_xxyz = intsBufferSSSG.data(intsIndexesSSSG(10));

    t_0_0_0_xxyy = intsBufferSSSG.data(intsIndexesSSSG(11));

    t_0_0_0_xxxz = intsBufferSSSG.data(intsIndexesSSSG(12));

    t_0_0_0_xxxy = intsBufferSSSG.data(intsIndexesSSSG(13));

    t_0_0_0_xxxx = intsBufferSSSG.data(intsIndexesSSSG(14));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_0_0_xxx, t_0_0_0_xxxx, t_0_0_0_xxxy,\
                           t_0_0_0_xxxz, t_0_0_0_xxy, t_0_0_0_xxyy, t_0_0_0_xxyz, t_0_0_0_xxz,\
                           t_0_0_0_xxzz, t_0_0_0_xyy, t_0_0_0_xyyy, t_0_0_0_xyyz, t_0_0_0_xyz,\
                           t_0_0_0_xyzz, t_0_0_0_xzz, t_0_0_0_xzzz, t_0_0_0_yyy, t_0_0_0_yyyy,\
                           t_0_0_0_yyyz, t_0_0_0_yyz, t_0_0_0_yyzz, t_0_0_0_yzz, t_0_0_0_yzzz,\
                           t_0_0_0_zzz, t_0_0_0_zzzz, t_0_0_x_xxx, t_0_0_x_xxy, t_0_0_x_xxz,\
                           t_0_0_x_xyy, t_0_0_x_xyz, t_0_0_x_xzz, t_0_0_x_yyy, t_0_0_x_yyz,\
                           t_0_0_x_yzz, t_0_0_x_zzz, t_0_0_y_xxy, t_0_0_y_xxz, t_0_0_y_xyy,\
                           t_0_0_y_xyz, t_0_0_y_xzz, t_0_0_y_yyy, t_0_0_y_yyz, t_0_0_y_yzz,\
                           t_0_0_y_zzz, t_0_0_z_xxz, t_0_0_z_xyz, t_0_0_z_xzz, t_0_0_z_yyz,\
                           t_0_0_z_yzz, t_0_0_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_0_z_zzz[i] = t_0_0_0_zzzz[i] - rcd_z[i] * t_0_0_0_zzz[i];

        t_0_0_z_yzz[i] = t_0_0_0_yzzz[i] - rcd_z[i] * t_0_0_0_yzz[i];

        t_0_0_z_yyz[i] = t_0_0_0_yyzz[i] - rcd_z[i] * t_0_0_0_yyz[i];

        t_0_0_z_xzz[i] = t_0_0_0_xzzz[i] - rcd_z[i] * t_0_0_0_xzz[i];

        t_0_0_z_xyz[i] = t_0_0_0_xyzz[i] - rcd_z[i] * t_0_0_0_xyz[i];

        t_0_0_z_xxz[i] = t_0_0_0_xxzz[i] - rcd_z[i] * t_0_0_0_xxz[i];

        t_0_0_y_zzz[i] = t_0_0_0_yzzz[i] - rcd_y[i] * t_0_0_0_zzz[i];

        t_0_0_y_yzz[i] = t_0_0_0_yyzz[i] - rcd_y[i] * t_0_0_0_yzz[i];

        t_0_0_y_yyz[i] = t_0_0_0_yyyz[i] - rcd_y[i] * t_0_0_0_yyz[i];

        t_0_0_y_yyy[i] = t_0_0_0_yyyy[i] - rcd_y[i] * t_0_0_0_yyy[i];

        t_0_0_y_xzz[i] = t_0_0_0_xyzz[i] - rcd_y[i] * t_0_0_0_xzz[i];

        t_0_0_y_xyz[i] = t_0_0_0_xyyz[i] - rcd_y[i] * t_0_0_0_xyz[i];

        t_0_0_y_xyy[i] = t_0_0_0_xyyy[i] - rcd_y[i] * t_0_0_0_xyy[i];

        t_0_0_y_xxz[i] = t_0_0_0_xxyz[i] - rcd_y[i] * t_0_0_0_xxz[i];

        t_0_0_y_xxy[i] = t_0_0_0_xxyy[i] - rcd_y[i] * t_0_0_0_xxy[i];

        t_0_0_x_zzz[i] = t_0_0_0_xzzz[i] - rcd_x[i] * t_0_0_0_zzz[i];

        t_0_0_x_yzz[i] = t_0_0_0_xyzz[i] - rcd_x[i] * t_0_0_0_yzz[i];

        t_0_0_x_yyz[i] = t_0_0_0_xyyz[i] - rcd_x[i] * t_0_0_0_yyz[i];

        t_0_0_x_yyy[i] = t_0_0_0_xyyy[i] - rcd_x[i] * t_0_0_0_yyy[i];

        t_0_0_x_xzz[i] = t_0_0_0_xxzz[i] - rcd_x[i] * t_0_0_0_xzz[i];

        t_0_0_x_xyz[i] = t_0_0_0_xxyz[i] - rcd_x[i] * t_0_0_0_xyz[i];

        t_0_0_x_xyy[i] = t_0_0_0_xxyy[i] - rcd_x[i] * t_0_0_0_xyy[i];

        t_0_0_x_xxz[i] = t_0_0_0_xxxz[i] - rcd_x[i] * t_0_0_0_xxz[i];

        t_0_0_x_xxy[i] = t_0_0_0_xxxy[i] - rcd_x[i] * t_0_0_0_xxy[i];

        t_0_0_x_xxx[i] = t_0_0_0_xxxx[i] - rcd_x[i] * t_0_0_0_xxx[i];
    }
}


} // derirec namespace
