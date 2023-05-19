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
compHostVRRForSGSS_V0(      BufferHostXY<T>&      intsBufferSGSS,
                      const BufferHostX<int32_t>& intsIndexesSGSS0,
                      const BufferHostXY<T>&      intsBufferSDSS0,
                      const BufferHostX<int32_t>& intsIndexesSDSS0,
                      const BufferHostXY<T>&      intsBufferSDSS1,
                      const BufferHostX<int32_t>& intsIndexesSDSS1,
                      const BufferHostXY<T>&      intsBufferSFSS0,
                      const BufferHostX<int32_t>& intsIndexesSFSS0,
                      const BufferHostXY<T>&      intsBufferSFSS1,
                      const BufferHostX<int32_t>& intsIndexesSFSS1,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fz_0 = osFactorsBraZeta;

    auto frz2_0 = osFactorsBraRhoZeta;

    // set up R(PB) distances

    auto rpb_z = rDistancesPB.data(2);

    auto rpb_y = rDistancesPB.data(1);

    auto rpb_x = rDistancesPB.data(0);

    // set up R(WP) distances

    auto rwp_z = rDistancesWP.data(2);

    auto rwp_y = rDistancesWP.data(1);

    auto rwp_x = rDistancesWP.data(0);

    // set up [SGSS]^(0) integral components

    t_0_zzzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(0));

    t_0_yzzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(1));

    t_0_yyzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(2));

    t_0_yyyz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(3));

    t_0_yyyy_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(4));

    t_0_xzzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(5));

    t_0_xyzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(6));

    t_0_xyyz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(7));

    t_0_xyyy_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(8));

    t_0_xxzz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(9));

    t_0_xxyz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(10));

    t_0_xxyy_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(11));

    t_0_xxxz_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(12));

    t_0_xxxy_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(13));

    t_0_xxxx_0_0_0 = intsBufferSGSS0.data(intsIndexesSGSS0(14));

    // set up [SDSS]^(0) integral components

    t_0_zz_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(0));

    t_0_yy_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(1));

    t_0_xx_0_0_0 = intsBufferSDSS0.data(intsIndexesSDSS0(2));

    // set up [SDSS]^(1) integral components

    t_0_zz_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(0));

    t_0_yy_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(1));

    t_0_xx_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(2));

    // set up [SFSS]^(0) integral components

    t_0_zzz_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(0));

    t_0_yzz_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(1));

    t_0_yyz_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(2));

    t_0_yyy_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(3));

    t_0_xzz_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(4));

    t_0_xyy_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(5));

    t_0_xxz_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(6));

    t_0_xxx_0_0_0 = intsBufferSFSS0.data(intsIndexesSFSS0(7));

    // set up [SFSS]^(1) integral components

    t_0_zzz_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(0));

    t_0_yzz_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(1));

    t_0_yyz_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(2));

    t_0_yyy_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(3));

    t_0_xzz_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(4));

    t_0_xyy_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(5));

    t_0_xxz_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(6));

    t_0_xxx_0_0_1 = intsBufferSFSS1.data(intsIndexesSFSS1(7));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    const auto fact_3_2 = static_cast<T>(3.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_xx_0_0_0, t_0_xx_0_0_1, t_0_xxx_0_0_0, t_0_xxx_0_0_1,\
                               t_0_xxxx_0_0_0, t_0_xxxy_0_0_0, t_0_xxxz_0_0_0, t_0_xxyy_0_0_0,\
                               t_0_xxyz_0_0_0, t_0_xxz_0_0_0, t_0_xxz_0_0_1, t_0_xxzz_0_0_0,\
                               t_0_xyy_0_0_0, t_0_xyy_0_0_1, t_0_xyyy_0_0_0, t_0_xyyz_0_0_0,\
                               t_0_xyzz_0_0_0, t_0_xzz_0_0_0, t_0_xzz_0_0_1, t_0_xzzz_0_0_0,\
                               t_0_yy_0_0_0, t_0_yy_0_0_1, t_0_yyy_0_0_0, t_0_yyy_0_0_1,\
                               t_0_yyyy_0_0_0, t_0_yyyz_0_0_0, t_0_yyz_0_0_0, t_0_yyz_0_0_1,\
                               t_0_yyzz_0_0_0, t_0_yzz_0_0_0, t_0_yzz_0_0_1, t_0_yzzz_0_0_0,\
                               t_0_zz_0_0_0, t_0_zz_0_0_1, t_0_zzz_0_0_0, t_0_zzz_0_0_1,\
                               t_0_zzzz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_0_0[i] += rpb_z[i] * t_0_zzz_0_0_0[i] + rwp_z[i] * t_0_zzz_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_yzzz_0_0_0[i] += rpb_y[i] * t_0_zzz_0_0_0[i] + rwp_y[i] * t_0_zzz_0_0_1[i];

            t_0_yyzz_0_0_0[i] += rpb_y[i] * t_0_yzz_0_0_0[i] + rwp_y[i] * t_0_yzz_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_yyyz_0_0_0[i] += rpb_z[i] * t_0_yyy_0_0_0[i] + rwp_z[i] * t_0_yyy_0_0_1[i];

            t_0_yyyy_0_0_0[i] += rpb_y[i] * t_0_yyy_0_0_0[i] + rwp_y[i] * t_0_yyy_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_0_1[i];

            t_0_xzzz_0_0_0[i] += rpb_x[i] * t_0_zzz_0_0_0[i] + rwp_x[i] * t_0_zzz_0_0_1[i];

            t_0_xyzz_0_0_0[i] += rpb_x[i] * t_0_yzz_0_0_0[i] + rwp_x[i] * t_0_yzz_0_0_1[i];

            t_0_xyyz_0_0_0[i] += rpb_x[i] * t_0_yyz_0_0_0[i] + rwp_x[i] * t_0_yyz_0_0_1[i];

            t_0_xyyy_0_0_0[i] += rpb_x[i] * t_0_yyy_0_0_0[i] + rwp_x[i] * t_0_yyy_0_0_1[i];

            t_0_xxzz_0_0_0[i] += rpb_x[i] * t_0_xzz_0_0_0[i] + rwp_x[i] * t_0_xzz_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_xxyz_0_0_0[i] += rpb_y[i] * t_0_xxz_0_0_0[i] + rwp_y[i] * t_0_xxz_0_0_1[i];

            t_0_xxyy_0_0_0[i] += rpb_x[i] * t_0_xyy_0_0_0[i] + rwp_x[i] * t_0_xyy_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_0_1[i];

            t_0_xxxz_0_0_0[i] += rpb_z[i] * t_0_xxx_0_0_0[i] + rwp_z[i] * t_0_xxx_0_0_1[i];

            t_0_xxxy_0_0_0[i] += rpb_y[i] * t_0_xxx_0_0_0[i] + rwp_y[i] * t_0_xxx_0_0_1[i];

            t_0_xxxx_0_0_0[i] += rpb_x[i] * t_0_xxx_0_0_0[i] + rwp_x[i] * t_0_xxx_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y, rwp_z,\
                               t_0_xx_0_0_0, t_0_xx_0_0_1, t_0_xxx_0_0_0, t_0_xxx_0_0_1,\
                               t_0_xxxx_0_0_0, t_0_xxxy_0_0_0, t_0_xxxz_0_0_0, t_0_xxyy_0_0_0,\
                               t_0_xxyz_0_0_0, t_0_xxz_0_0_0, t_0_xxz_0_0_1, t_0_xxzz_0_0_0,\
                               t_0_xyy_0_0_0, t_0_xyy_0_0_1, t_0_xyyy_0_0_0, t_0_xyyz_0_0_0,\
                               t_0_xyzz_0_0_0, t_0_xzz_0_0_0, t_0_xzz_0_0_1, t_0_xzzz_0_0_0,\
                               t_0_yy_0_0_0, t_0_yy_0_0_1, t_0_yyy_0_0_0, t_0_yyy_0_0_1,\
                               t_0_yyyy_0_0_0, t_0_yyyz_0_0_0, t_0_yyz_0_0_0, t_0_yyz_0_0_1,\
                               t_0_yyzz_0_0_0, t_0_yzz_0_0_0, t_0_yzz_0_0_1, t_0_yzzz_0_0_0,\
                               t_0_zz_0_0_0, t_0_zz_0_0_1, t_0_zzz_0_0_0, t_0_zzz_0_0_1,\
                               t_0_zzzz_0_0_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzzz_0_0_0[i] = rpb_z[i] * t_0_zzz_0_0_0[i] + rwp_z[i] * t_0_zzz_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_yzzz_0_0_0[i] = rpb_y[i] * t_0_zzz_0_0_0[i] + rwp_y[i] * t_0_zzz_0_0_1[i];

            t_0_yyzz_0_0_0[i] = rpb_y[i] * t_0_yzz_0_0_0[i] + rwp_y[i] * t_0_yzz_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_yyyz_0_0_0[i] = rpb_z[i] * t_0_yyy_0_0_0[i] + rwp_z[i] * t_0_yyy_0_0_1[i];

            t_0_yyyy_0_0_0[i] = rpb_y[i] * t_0_yyy_0_0_0[i] + rwp_y[i] * t_0_yyy_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_yy_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_yy_0_0_1[i];

            t_0_xzzz_0_0_0[i] = rpb_x[i] * t_0_zzz_0_0_0[i] + rwp_x[i] * t_0_zzz_0_0_1[i];

            t_0_xyzz_0_0_0[i] = rpb_x[i] * t_0_yzz_0_0_0[i] + rwp_x[i] * t_0_yzz_0_0_1[i];

            t_0_xyyz_0_0_0[i] = rpb_x[i] * t_0_yyz_0_0_0[i] + rwp_x[i] * t_0_yyz_0_0_1[i];

            t_0_xyyy_0_0_0[i] = rpb_x[i] * t_0_yyy_0_0_0[i] + rwp_x[i] * t_0_yyy_0_0_1[i];

            t_0_xxzz_0_0_0[i] = rpb_x[i] * t_0_xzz_0_0_0[i] + rwp_x[i] * t_0_xzz_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_zz_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_zz_0_0_1[i];

            t_0_xxyz_0_0_0[i] = rpb_y[i] * t_0_xxz_0_0_0[i] + rwp_y[i] * t_0_xxz_0_0_1[i];

            t_0_xxyy_0_0_0[i] = rpb_x[i] * t_0_xyy_0_0_0[i] + rwp_x[i] * t_0_xyy_0_0_1[i] + fact_1_2 * fz_0[i] * t_0_yy_0_0_0[i] - fact_1_2 * frz2_0[i] * t_0_yy_0_0_1[i];

            t_0_xxxz_0_0_0[i] = rpb_z[i] * t_0_xxx_0_0_0[i] + rwp_z[i] * t_0_xxx_0_0_1[i];

            t_0_xxxy_0_0_0[i] = rpb_y[i] * t_0_xxx_0_0_0[i] + rwp_y[i] * t_0_xxx_0_0_1[i];

            t_0_xxxx_0_0_0[i] = rpb_x[i] * t_0_xxx_0_0_0[i] + rwp_x[i] * t_0_xxx_0_0_1[i] + fact_3_2 * fz_0[i] * t_0_xx_0_0_0[i] - fact_3_2 * frz2_0[i] * t_0_xx_0_0_1[i];
        }
    }
}


} // derirec namespace
