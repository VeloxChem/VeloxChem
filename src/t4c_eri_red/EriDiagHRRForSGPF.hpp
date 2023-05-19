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
compHostHRRForSGPF_V0(      BufferHostXY<T>&      intsBufferSGPF,
                      const BufferHostX<int32_t>& intsIndexesSGPF,
                      const BufferHostXY<T>&      intsBufferSGSF,
                      const BufferHostX<int32_t>& intsIndexesSGSF,
                      const BufferHostXY<T>&      intsBufferSGSG,
                      const BufferHostX<int32_t>& intsIndexesSGSG,
                      const BufferHostMY<T, 3>&   rDistancesCD,
                      const int32_t               nBatchPairs) -> void
{
    // set up R(CD) distances

    auto rcd_z = rDistancesCD.data(2);

    auto rcd_y = rDistancesCD.data(1);

    auto rcd_x = rDistancesCD.data(0);

    // set up (SGPF) integral components

    t_0_zzzz_z_zzz = intsBufferSGPF.data(intsIndexesSGPF(0));

    t_0_yzzz_z_yzz = intsBufferSGPF.data(intsIndexesSGPF(1));

    t_0_yzzz_y_zzz = intsBufferSGPF.data(intsIndexesSGPF(2));

    t_0_yyzz_z_yyz = intsBufferSGPF.data(intsIndexesSGPF(3));

    t_0_yyzz_y_yzz = intsBufferSGPF.data(intsIndexesSGPF(4));

    t_0_yyyz_y_yyz = intsBufferSGPF.data(intsIndexesSGPF(5));

    t_0_yyyy_y_yyy = intsBufferSGPF.data(intsIndexesSGPF(6));

    t_0_xzzz_z_xzz = intsBufferSGPF.data(intsIndexesSGPF(7));

    t_0_xzzz_x_zzz = intsBufferSGPF.data(intsIndexesSGPF(8));

    t_0_xyzz_z_xyz = intsBufferSGPF.data(intsIndexesSGPF(9));

    t_0_xyzz_y_xzz = intsBufferSGPF.data(intsIndexesSGPF(10));

    t_0_xyzz_x_yzz = intsBufferSGPF.data(intsIndexesSGPF(11));

    t_0_xyyz_y_xyz = intsBufferSGPF.data(intsIndexesSGPF(12));

    t_0_xyyz_x_yyz = intsBufferSGPF.data(intsIndexesSGPF(13));

    t_0_xyyy_y_xyy = intsBufferSGPF.data(intsIndexesSGPF(14));

    t_0_xyyy_x_yyy = intsBufferSGPF.data(intsIndexesSGPF(15));

    t_0_xxzz_z_xxz = intsBufferSGPF.data(intsIndexesSGPF(16));

    t_0_xxzz_x_xzz = intsBufferSGPF.data(intsIndexesSGPF(17));

    t_0_xxyz_y_xxz = intsBufferSGPF.data(intsIndexesSGPF(18));

    t_0_xxyz_x_xyz = intsBufferSGPF.data(intsIndexesSGPF(19));

    t_0_xxyy_y_xxy = intsBufferSGPF.data(intsIndexesSGPF(20));

    t_0_xxyy_x_xyy = intsBufferSGPF.data(intsIndexesSGPF(21));

    t_0_xxxz_x_xxz = intsBufferSGPF.data(intsIndexesSGPF(22));

    t_0_xxxy_x_xxy = intsBufferSGPF.data(intsIndexesSGPF(23));

    t_0_xxxx_x_xxx = intsBufferSGPF.data(intsIndexesSGPF(24));

    // set up (SGSF) integral components

    t_0_zzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(0));

    t_0_yzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(1));

    t_0_yzzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(2));

    t_0_yyzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(3));

    t_0_yyzz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(4));

    t_0_yyyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(5));

    t_0_yyyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(6));

    t_0_xzzz_0_zzz = intsBufferSGSF.data(intsIndexesSGSF(7));

    t_0_xzzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(8));

    t_0_xyzz_0_yzz = intsBufferSGSF.data(intsIndexesSGSF(9));

    t_0_xyzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(10));

    t_0_xyzz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(11));

    t_0_xyyz_0_yyz = intsBufferSGSF.data(intsIndexesSGSF(12));

    t_0_xyyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(13));

    t_0_xyyy_0_yyy = intsBufferSGSF.data(intsIndexesSGSF(14));

    t_0_xyyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(15));

    t_0_xxzz_0_xzz = intsBufferSGSF.data(intsIndexesSGSF(16));

    t_0_xxzz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(17));

    t_0_xxyz_0_xyz = intsBufferSGSF.data(intsIndexesSGSF(18));

    t_0_xxyz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(19));

    t_0_xxyy_0_xyy = intsBufferSGSF.data(intsIndexesSGSF(20));

    t_0_xxyy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(21));

    t_0_xxxz_0_xxz = intsBufferSGSF.data(intsIndexesSGSF(22));

    t_0_xxxy_0_xxy = intsBufferSGSF.data(intsIndexesSGSF(23));

    t_0_xxxx_0_xxx = intsBufferSGSF.data(intsIndexesSGSF(24));

    // set up (SGSG) integral components

    t_0_zzzz_0_zzzz = intsBufferSGSG.data(intsIndexesSGSG(0));

    t_0_yzzz_0_yzzz = intsBufferSGSG.data(intsIndexesSGSG(1));

    t_0_yyzz_0_yyzz = intsBufferSGSG.data(intsIndexesSGSG(2));

    t_0_yyyz_0_yyyz = intsBufferSGSG.data(intsIndexesSGSG(3));

    t_0_yyyy_0_yyyy = intsBufferSGSG.data(intsIndexesSGSG(4));

    t_0_xzzz_0_xzzz = intsBufferSGSG.data(intsIndexesSGSG(5));

    t_0_xyzz_0_xyzz = intsBufferSGSG.data(intsIndexesSGSG(6));

    t_0_xyyz_0_xyyz = intsBufferSGSG.data(intsIndexesSGSG(7));

    t_0_xyyy_0_xyyy = intsBufferSGSG.data(intsIndexesSGSG(8));

    t_0_xxzz_0_xxzz = intsBufferSGSG.data(intsIndexesSGSG(9));

    t_0_xxyz_0_xxyz = intsBufferSGSG.data(intsIndexesSGSG(10));

    t_0_xxyy_0_xxyy = intsBufferSGSG.data(intsIndexesSGSG(11));

    t_0_xxxz_0_xxxz = intsBufferSGSG.data(intsIndexesSGSG(12));

    t_0_xxxy_0_xxxy = intsBufferSGSG.data(intsIndexesSGSG(13));

    t_0_xxxx_0_xxxx = intsBufferSGSG.data(intsIndexesSGSG(14));

    #pragma omp simd align(rcd_x, rcd_y, rcd_z, t_0_xxxx_0_xxx, t_0_xxxx_0_xxxx, t_0_xxxx_x_xxx,\
                           t_0_xxxy_0_xxxy, t_0_xxxy_0_xxy, t_0_xxxy_x_xxy, t_0_xxxz_0_xxxz,\
                           t_0_xxxz_0_xxz, t_0_xxxz_x_xxz, t_0_xxyy_0_xxy, t_0_xxyy_0_xxyy,\
                           t_0_xxyy_0_xyy, t_0_xxyy_x_xyy, t_0_xxyy_y_xxy, t_0_xxyz_0_xxyz,\
                           t_0_xxyz_0_xxz, t_0_xxyz_0_xyz, t_0_xxyz_x_xyz, t_0_xxyz_y_xxz,\
                           t_0_xxzz_0_xxz, t_0_xxzz_0_xxzz, t_0_xxzz_0_xzz, t_0_xxzz_x_xzz,\
                           t_0_xxzz_z_xxz, t_0_xyyy_0_xyy, t_0_xyyy_0_xyyy, t_0_xyyy_0_yyy,\
                           t_0_xyyy_x_yyy, t_0_xyyy_y_xyy, t_0_xyyz_0_xyyz, t_0_xyyz_0_xyz,\
                           t_0_xyyz_0_yyz, t_0_xyyz_x_yyz, t_0_xyyz_y_xyz, t_0_xyzz_0_xyz,\
                           t_0_xyzz_0_xyzz, t_0_xyzz_0_xzz, t_0_xyzz_0_yzz, t_0_xyzz_x_yzz,\
                           t_0_xyzz_y_xzz, t_0_xyzz_z_xyz, t_0_xzzz_0_xzz, t_0_xzzz_0_xzzz,\
                           t_0_xzzz_0_zzz, t_0_xzzz_x_zzz, t_0_xzzz_z_xzz, t_0_yyyy_0_yyy,\
                           t_0_yyyy_0_yyyy, t_0_yyyy_y_yyy, t_0_yyyz_0_yyyz, t_0_yyyz_0_yyz,\
                           t_0_yyyz_y_yyz, t_0_yyzz_0_yyz, t_0_yyzz_0_yyzz, t_0_yyzz_0_yzz,\
                           t_0_yyzz_y_yzz, t_0_yyzz_z_yyz, t_0_yzzz_0_yzz, t_0_yzzz_0_yzzz,\
                           t_0_yzzz_0_zzz, t_0_yzzz_y_zzz, t_0_yzzz_z_yzz, t_0_zzzz_0_zzz,\
                           t_0_zzzz_0_zzzz, t_0_zzzz_z_zzz : VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        t_0_zzzz_z_zzz[i] = t_0_zzzz_0_zzzz[i] - rcd_z[i] * t_0_zzzz_0_zzz[i];

        t_0_yzzz_z_yzz[i] = t_0_yzzz_0_yzzz[i] - rcd_z[i] * t_0_yzzz_0_yzz[i];

        t_0_yzzz_y_zzz[i] = t_0_yzzz_0_yzzz[i] - rcd_y[i] * t_0_yzzz_0_zzz[i];

        t_0_yyzz_z_yyz[i] = t_0_yyzz_0_yyzz[i] - rcd_z[i] * t_0_yyzz_0_yyz[i];

        t_0_yyzz_y_yzz[i] = t_0_yyzz_0_yyzz[i] - rcd_y[i] * t_0_yyzz_0_yzz[i];

        t_0_yyyz_y_yyz[i] = t_0_yyyz_0_yyyz[i] - rcd_y[i] * t_0_yyyz_0_yyz[i];

        t_0_yyyy_y_yyy[i] = t_0_yyyy_0_yyyy[i] - rcd_y[i] * t_0_yyyy_0_yyy[i];

        t_0_xzzz_z_xzz[i] = t_0_xzzz_0_xzzz[i] - rcd_z[i] * t_0_xzzz_0_xzz[i];

        t_0_xzzz_x_zzz[i] = t_0_xzzz_0_xzzz[i] - rcd_x[i] * t_0_xzzz_0_zzz[i];

        t_0_xyzz_z_xyz[i] = t_0_xyzz_0_xyzz[i] - rcd_z[i] * t_0_xyzz_0_xyz[i];

        t_0_xyzz_y_xzz[i] = t_0_xyzz_0_xyzz[i] - rcd_y[i] * t_0_xyzz_0_xzz[i];

        t_0_xyzz_x_yzz[i] = t_0_xyzz_0_xyzz[i] - rcd_x[i] * t_0_xyzz_0_yzz[i];

        t_0_xyyz_y_xyz[i] = t_0_xyyz_0_xyyz[i] - rcd_y[i] * t_0_xyyz_0_xyz[i];

        t_0_xyyz_x_yyz[i] = t_0_xyyz_0_xyyz[i] - rcd_x[i] * t_0_xyyz_0_yyz[i];

        t_0_xyyy_y_xyy[i] = t_0_xyyy_0_xyyy[i] - rcd_y[i] * t_0_xyyy_0_xyy[i];

        t_0_xyyy_x_yyy[i] = t_0_xyyy_0_xyyy[i] - rcd_x[i] * t_0_xyyy_0_yyy[i];

        t_0_xxzz_z_xxz[i] = t_0_xxzz_0_xxzz[i] - rcd_z[i] * t_0_xxzz_0_xxz[i];

        t_0_xxzz_x_xzz[i] = t_0_xxzz_0_xxzz[i] - rcd_x[i] * t_0_xxzz_0_xzz[i];

        t_0_xxyz_y_xxz[i] = t_0_xxyz_0_xxyz[i] - rcd_y[i] * t_0_xxyz_0_xxz[i];

        t_0_xxyz_x_xyz[i] = t_0_xxyz_0_xxyz[i] - rcd_x[i] * t_0_xxyz_0_xyz[i];

        t_0_xxyy_y_xxy[i] = t_0_xxyy_0_xxyy[i] - rcd_y[i] * t_0_xxyy_0_xxy[i];

        t_0_xxyy_x_xyy[i] = t_0_xxyy_0_xxyy[i] - rcd_x[i] * t_0_xxyy_0_xyy[i];

        t_0_xxxz_x_xxz[i] = t_0_xxxz_0_xxxz[i] - rcd_x[i] * t_0_xxxz_0_xxz[i];

        t_0_xxxy_x_xxy[i] = t_0_xxxy_0_xxxy[i] - rcd_x[i] * t_0_xxxy_0_xxy[i];

        t_0_xxxx_x_xxx[i] = t_0_xxxx_0_xxxx[i] - rcd_x[i] * t_0_xxxx_0_xxx[i];
    }
}


} // derirec namespace
