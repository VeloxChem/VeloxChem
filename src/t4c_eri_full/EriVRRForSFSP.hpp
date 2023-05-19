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
compHostVRRForSFSP_V0(      BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP0,
                      const BufferHostXY<T>&      intsBufferSPSP0,
                      const BufferHostX<int32_t>& intsIndexesSPSP0,
                      const BufferHostXY<T>&      intsBufferSPSP1,
                      const BufferHostX<int32_t>& intsIndexesSPSP1,
                      const BufferHostXY<T>&      intsBufferSDSS1,
                      const BufferHostX<int32_t>& intsIndexesSDSS1,
                      const BufferHostXY<T>&      intsBufferSDSP0,
                      const BufferHostX<int32_t>& intsIndexesSDSP0,
                      const BufferHostXY<T>&      intsBufferSDSP1,
                      const BufferHostX<int32_t>& intsIndexesSDSP1,
                      const T*                    osFactorsZeta,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

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

    // set up [SFSP]^(0) integral components

    t_0_zzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(0));

    t_0_zzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(1));

    t_0_zzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(2));

    t_0_yzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(3));

    t_0_yzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(4));

    t_0_yzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(5));

    t_0_yyz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(6));

    t_0_yyz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(7));

    t_0_yyz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(8));

    t_0_yyy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(9));

    t_0_yyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(10));

    t_0_yyy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(11));

    t_0_xzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(12));

    t_0_xzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(13));

    t_0_xzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(14));

    t_0_xyz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(15));

    t_0_xyz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(16));

    t_0_xyz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(17));

    t_0_xyy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(18));

    t_0_xyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(19));

    t_0_xyy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(20));

    t_0_xxz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(21));

    t_0_xxz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(22));

    t_0_xxz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(23));

    t_0_xxy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(24));

    t_0_xxy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(25));

    t_0_xxy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(26));

    t_0_xxx_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(27));

    t_0_xxx_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(28));

    t_0_xxx_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(29));

    // set up [SPSP]^(0) integral components

    t_0_z_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(0));

    t_0_z_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(1));

    t_0_z_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(2));

    t_0_y_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(3));

    t_0_y_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(4));

    t_0_y_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(5));

    t_0_x_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(6));

    t_0_x_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(7));

    t_0_x_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(8));

    // set up [SPSP]^(1) integral components

    t_0_z_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(0));

    t_0_z_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(1));

    t_0_z_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(2));

    t_0_y_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(3));

    t_0_y_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(4));

    t_0_y_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(5));

    t_0_x_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(6));

    t_0_x_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(7));

    t_0_x_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(8));

    // set up [SDSS]^(1) integral components

    t_0_zz_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(0));

    t_0_yy_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(1));

    t_0_xx_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(2));

    // set up [SDSP]^(0) integral components

    t_0_zz_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(0));

    t_0_zz_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(1));

    t_0_zz_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(2));

    t_0_yz_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(3));

    t_0_yz_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(4));

    t_0_yy_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(5));

    t_0_yy_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(6));

    t_0_yy_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(7));

    t_0_xz_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(8));

    t_0_xx_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(9));

    t_0_xx_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(10));

    t_0_xx_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(11));

    // set up [SDSP]^(1) integral components

    t_0_zz_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(0));

    t_0_zz_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(1));

    t_0_zz_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(2));

    t_0_yz_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(3));

    t_0_yz_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(4));

    t_0_yy_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(5));

    t_0_yy_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(6));

    t_0_yy_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(7));

    t_0_xz_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(8));

    t_0_xx_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(9));

    t_0_xx_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(10));

    t_0_xx_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(11));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxy_0_x_0, t_0_xxy_0_y_0, t_0_xxy_0_z_0, t_0_xxz_0_x_0,\
                               t_0_xxz_0_y_0, t_0_xxz_0_z_0, t_0_xyy_0_x_0, t_0_xyy_0_y_0,\
                               t_0_xyy_0_z_0, t_0_xyz_0_x_0, t_0_xyz_0_y_0, t_0_xyz_0_z_0,\
                               t_0_xz_0_x_0, t_0_xz_0_x_1, t_0_xzz_0_x_0, t_0_xzz_0_y_0,\
                               t_0_xzz_0_z_0, t_0_y_0_x_0, t_0_y_0_x_1, t_0_y_0_y_0,\
                               t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1, t_0_yy_0_0_1,\
                               t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0, t_0_yy_0_y_1,\
                               t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0, t_0_yyy_0_y_0,\
                               t_0_yyy_0_z_0, t_0_yyz_0_x_0, t_0_yyz_0_y_0, t_0_yyz_0_z_0,\
                               t_0_yz_0_y_0, t_0_yz_0_y_1, t_0_yz_0_z_0, t_0_yz_0_z_1,\
                               t_0_yzz_0_x_0, t_0_yzz_0_y_0, t_0_yzz_0_z_0, t_0_z_0_x_0,\
                               t_0_z_0_x_1, t_0_z_0_y_0, t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1,\
                               t_0_zz_0_0_1, t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0,\
                               t_0_zz_0_y_1, t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0,\
                               t_0_zzz_0_y_0, t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] += rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] += rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] += rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] += rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] += rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yzz_0_x_0[i] += rpb_y[i] * t_0_zz_0_x_0[i] + rwp_y[i] * t_0_zz_0_x_1[i];

            t_0_yyz_0_z_0[i] += rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyz_0_y_0[i] += rpb_z[i] * t_0_yy_0_y_0[i] + rwp_z[i] * t_0_yy_0_y_1[i];

            t_0_yyz_0_x_0[i] += rpb_z[i] * t_0_yy_0_x_0[i] + rwp_z[i] * t_0_yy_0_x_1[i];

            t_0_yyy_0_z_0[i] += rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] += rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] += rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] += rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xzz_0_y_0[i] += rpb_x[i] * t_0_zz_0_y_0[i] + rwp_x[i] * t_0_zz_0_y_1[i];

            t_0_xzz_0_x_0[i] += rpb_x[i] * t_0_zz_0_x_0[i] + rwp_x[i] * t_0_zz_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_xyz_0_z_0[i] += rpb_x[i] * t_0_yz_0_z_0[i] + rwp_x[i] * t_0_yz_0_z_1[i];

            t_0_xyz_0_y_0[i] += rpb_x[i] * t_0_yz_0_y_0[i] + rwp_x[i] * t_0_yz_0_y_1[i];

            t_0_xyz_0_x_0[i] += rpb_y[i] * t_0_xz_0_x_0[i] + rwp_y[i] * t_0_xz_0_x_1[i];

            t_0_xyy_0_z_0[i] += rpb_x[i] * t_0_yy_0_z_0[i] + rwp_x[i] * t_0_yy_0_z_1[i];

            t_0_xyy_0_y_0[i] += rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xyy_0_x_0[i] += rpb_x[i] * t_0_yy_0_x_0[i] + rwp_x[i] * t_0_yy_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_xxz_0_z_0[i] += rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxz_0_y_0[i] += rpb_z[i] * t_0_xx_0_y_0[i] + rwp_z[i] * t_0_xx_0_y_1[i];

            t_0_xxz_0_x_0[i] += rpb_z[i] * t_0_xx_0_x_0[i] + rwp_z[i] * t_0_xx_0_x_1[i];

            t_0_xxy_0_z_0[i] += rpb_y[i] * t_0_xx_0_z_0[i] + rwp_y[i] * t_0_xx_0_z_1[i];

            t_0_xxy_0_y_0[i] += rpb_y[i] * t_0_xx_0_y_0[i] + rwp_y[i] * t_0_xx_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxy_0_x_0[i] += rpb_y[i] * t_0_xx_0_x_0[i] + rwp_y[i] * t_0_xx_0_x_1[i];

            t_0_xxx_0_z_0[i] += rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] += rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] += rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxy_0_x_0, t_0_xxy_0_y_0, t_0_xxy_0_z_0, t_0_xxz_0_x_0,\
                               t_0_xxz_0_y_0, t_0_xxz_0_z_0, t_0_xyy_0_x_0, t_0_xyy_0_y_0,\
                               t_0_xyy_0_z_0, t_0_xyz_0_x_0, t_0_xyz_0_y_0, t_0_xyz_0_z_0,\
                               t_0_xz_0_x_0, t_0_xz_0_x_1, t_0_xzz_0_x_0, t_0_xzz_0_y_0,\
                               t_0_xzz_0_z_0, t_0_y_0_x_0, t_0_y_0_x_1, t_0_y_0_y_0,\
                               t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1, t_0_yy_0_0_1,\
                               t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0, t_0_yy_0_y_1,\
                               t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0, t_0_yyy_0_y_0,\
                               t_0_yyy_0_z_0, t_0_yyz_0_x_0, t_0_yyz_0_y_0, t_0_yyz_0_z_0,\
                               t_0_yz_0_y_0, t_0_yz_0_y_1, t_0_yz_0_z_0, t_0_yz_0_z_1,\
                               t_0_yzz_0_x_0, t_0_yzz_0_y_0, t_0_yzz_0_z_0, t_0_z_0_x_0,\
                               t_0_z_0_x_1, t_0_z_0_y_0, t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1,\
                               t_0_zz_0_0_1, t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0,\
                               t_0_zz_0_y_1, t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0,\
                               t_0_zzz_0_y_0, t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] = rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] = rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] = rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] = rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] = rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yzz_0_x_0[i] = rpb_y[i] * t_0_zz_0_x_0[i] + rwp_y[i] * t_0_zz_0_x_1[i];

            t_0_yyz_0_z_0[i] = rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyz_0_y_0[i] = rpb_z[i] * t_0_yy_0_y_0[i] + rwp_z[i] * t_0_yy_0_y_1[i];

            t_0_yyz_0_x_0[i] = rpb_z[i] * t_0_yy_0_x_0[i] + rwp_z[i] * t_0_yy_0_x_1[i];

            t_0_yyy_0_z_0[i] = rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] = rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] = rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] = rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xzz_0_y_0[i] = rpb_x[i] * t_0_zz_0_y_0[i] + rwp_x[i] * t_0_zz_0_y_1[i];

            t_0_xzz_0_x_0[i] = rpb_x[i] * t_0_zz_0_x_0[i] + rwp_x[i] * t_0_zz_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_xyz_0_z_0[i] = rpb_x[i] * t_0_yz_0_z_0[i] + rwp_x[i] * t_0_yz_0_z_1[i];

            t_0_xyz_0_y_0[i] = rpb_x[i] * t_0_yz_0_y_0[i] + rwp_x[i] * t_0_yz_0_y_1[i];

            t_0_xyz_0_x_0[i] = rpb_y[i] * t_0_xz_0_x_0[i] + rwp_y[i] * t_0_xz_0_x_1[i];

            t_0_xyy_0_z_0[i] = rpb_x[i] * t_0_yy_0_z_0[i] + rwp_x[i] * t_0_yy_0_z_1[i];

            t_0_xyy_0_y_0[i] = rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xyy_0_x_0[i] = rpb_x[i] * t_0_yy_0_x_0[i] + rwp_x[i] * t_0_yy_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_xxz_0_z_0[i] = rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxz_0_y_0[i] = rpb_z[i] * t_0_xx_0_y_0[i] + rwp_z[i] * t_0_xx_0_y_1[i];

            t_0_xxz_0_x_0[i] = rpb_z[i] * t_0_xx_0_x_0[i] + rwp_z[i] * t_0_xx_0_x_1[i];

            t_0_xxy_0_z_0[i] = rpb_y[i] * t_0_xx_0_z_0[i] + rwp_y[i] * t_0_xx_0_z_1[i];

            t_0_xxy_0_y_0[i] = rpb_y[i] * t_0_xx_0_y_0[i] + rwp_y[i] * t_0_xx_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxy_0_x_0[i] = rpb_y[i] * t_0_xx_0_x_0[i] + rwp_y[i] * t_0_xx_0_x_1[i];

            t_0_xxx_0_z_0[i] = rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] = rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] = rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSFSP_V1(      BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP0,
                      const BufferHostXY<T>&      intsBufferSPSP0,
                      const BufferHostX<int32_t>& intsIndexesSPSP0,
                      const BufferHostXY<T>&      intsBufferSPSP1,
                      const BufferHostX<int32_t>& intsIndexesSPSP1,
                      const BufferHostXY<T>&      intsBufferSDSS1,
                      const BufferHostX<int32_t>& intsIndexesSDSS1,
                      const BufferHostXY<T>&      intsBufferSDSP0,
                      const BufferHostX<int32_t>& intsIndexesSDSP0,
                      const BufferHostXY<T>&      intsBufferSDSP1,
                      const BufferHostX<int32_t>& intsIndexesSDSP1,
                      const T*                    osFactorsZeta,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

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

    // set up [SFSP]^(0) integral components

    t_0_zzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(0));

    t_0_zzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(1));

    t_0_zzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(2));

    t_0_yzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(3));

    t_0_yzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(4));

    t_0_yzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(5));

    t_0_yyz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(6));

    t_0_yyz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(7));

    t_0_yyy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(8));

    t_0_yyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(9));

    t_0_yyy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(10));

    t_0_xzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(11));

    t_0_xzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(12));

    t_0_xzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(13));

    t_0_xyy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(14));

    t_0_xyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(15));

    t_0_xyy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(16));

    t_0_xxz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(17));

    t_0_xxz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(18));

    t_0_xxy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(19));

    t_0_xxy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(20));

    t_0_xxx_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(21));

    t_0_xxx_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(22));

    t_0_xxx_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(23));

    // set up [SPSP]^(0) integral components

    t_0_z_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(0));

    t_0_z_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(1));

    t_0_z_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(2));

    t_0_y_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(3));

    t_0_y_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(4));

    t_0_y_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(5));

    t_0_x_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(6));

    t_0_x_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(7));

    t_0_x_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(8));

    // set up [SPSP]^(1) integral components

    t_0_z_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(0));

    t_0_z_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(1));

    t_0_z_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(2));

    t_0_y_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(3));

    t_0_y_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(4));

    t_0_y_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(5));

    t_0_x_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(6));

    t_0_x_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(7));

    t_0_x_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(8));

    // set up [SDSS]^(1) integral components

    t_0_zz_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(0));

    t_0_yy_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(1));

    t_0_xx_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(2));

    // set up [SDSP]^(0) integral components

    t_0_zz_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(0));

    t_0_zz_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(1));

    t_0_zz_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(2));

    t_0_yy_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(3));

    t_0_yy_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(4));

    t_0_yy_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(5));

    t_0_xx_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(6));

    t_0_xx_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(7));

    t_0_xx_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(8));

    // set up [SDSP]^(1) integral components

    t_0_zz_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(0));

    t_0_zz_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(1));

    t_0_zz_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(2));

    t_0_yy_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(3));

    t_0_yy_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(4));

    t_0_yy_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(5));

    t_0_xx_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(6));

    t_0_xx_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(7));

    t_0_xx_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(8));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxy_0_x_0, t_0_xxy_0_y_0, t_0_xxz_0_x_0, t_0_xxz_0_z_0,\
                               t_0_xyy_0_x_0, t_0_xyy_0_y_0, t_0_xyy_0_z_0, t_0_xzz_0_x_0,\
                               t_0_xzz_0_y_0, t_0_xzz_0_z_0, t_0_y_0_x_0, t_0_y_0_x_1,\
                               t_0_y_0_y_0, t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1, t_0_yy_0_0_1,\
                               t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0, t_0_yy_0_y_1,\
                               t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0, t_0_yyy_0_y_0,\
                               t_0_yyy_0_z_0, t_0_yyz_0_y_0, t_0_yyz_0_z_0, t_0_yzz_0_x_0,\
                               t_0_yzz_0_y_0, t_0_yzz_0_z_0, t_0_z_0_x_0, t_0_z_0_x_1,\
                               t_0_z_0_y_0, t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1, t_0_zz_0_0_1,\
                               t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0, t_0_zz_0_y_1,\
                               t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0, t_0_zzz_0_y_0,\
                               t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] += rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] += rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] += rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] += rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] += rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yzz_0_x_0[i] += rpb_y[i] * t_0_zz_0_x_0[i] + rwp_y[i] * t_0_zz_0_x_1[i];

            t_0_yyz_0_z_0[i] += rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyz_0_y_0[i] += rpb_z[i] * t_0_yy_0_y_0[i] + rwp_z[i] * t_0_yy_0_y_1[i];

            t_0_yyy_0_z_0[i] += rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] += rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] += rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] += rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xzz_0_y_0[i] += rpb_x[i] * t_0_zz_0_y_0[i] + rwp_x[i] * t_0_zz_0_y_1[i];

            t_0_xzz_0_x_0[i] += rpb_x[i] * t_0_zz_0_x_0[i] + rwp_x[i] * t_0_zz_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_xyy_0_z_0[i] += rpb_x[i] * t_0_yy_0_z_0[i] + rwp_x[i] * t_0_yy_0_z_1[i];

            t_0_xyy_0_y_0[i] += rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xyy_0_x_0[i] += rpb_x[i] * t_0_yy_0_x_0[i] + rwp_x[i] * t_0_yy_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_xxz_0_z_0[i] += rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxz_0_x_0[i] += rpb_z[i] * t_0_xx_0_x_0[i] + rwp_z[i] * t_0_xx_0_x_1[i];

            t_0_xxy_0_y_0[i] += rpb_y[i] * t_0_xx_0_y_0[i] + rwp_y[i] * t_0_xx_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxy_0_x_0[i] += rpb_y[i] * t_0_xx_0_x_0[i] + rwp_y[i] * t_0_xx_0_x_1[i];

            t_0_xxx_0_z_0[i] += rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] += rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] += rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxy_0_x_0, t_0_xxy_0_y_0, t_0_xxz_0_x_0, t_0_xxz_0_z_0,\
                               t_0_xyy_0_x_0, t_0_xyy_0_y_0, t_0_xyy_0_z_0, t_0_xzz_0_x_0,\
                               t_0_xzz_0_y_0, t_0_xzz_0_z_0, t_0_y_0_x_0, t_0_y_0_x_1,\
                               t_0_y_0_y_0, t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1, t_0_yy_0_0_1,\
                               t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0, t_0_yy_0_y_1,\
                               t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0, t_0_yyy_0_y_0,\
                               t_0_yyy_0_z_0, t_0_yyz_0_y_0, t_0_yyz_0_z_0, t_0_yzz_0_x_0,\
                               t_0_yzz_0_y_0, t_0_yzz_0_z_0, t_0_z_0_x_0, t_0_z_0_x_1,\
                               t_0_z_0_y_0, t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1, t_0_zz_0_0_1,\
                               t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0, t_0_zz_0_y_1,\
                               t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0, t_0_zzz_0_y_0,\
                               t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] = rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] = rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] = rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] = rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] = rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yzz_0_x_0[i] = rpb_y[i] * t_0_zz_0_x_0[i] + rwp_y[i] * t_0_zz_0_x_1[i];

            t_0_yyz_0_z_0[i] = rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyz_0_y_0[i] = rpb_z[i] * t_0_yy_0_y_0[i] + rwp_z[i] * t_0_yy_0_y_1[i];

            t_0_yyy_0_z_0[i] = rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] = rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] = rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] = rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xzz_0_y_0[i] = rpb_x[i] * t_0_zz_0_y_0[i] + rwp_x[i] * t_0_zz_0_y_1[i];

            t_0_xzz_0_x_0[i] = rpb_x[i] * t_0_zz_0_x_0[i] + rwp_x[i] * t_0_zz_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_xyy_0_z_0[i] = rpb_x[i] * t_0_yy_0_z_0[i] + rwp_x[i] * t_0_yy_0_z_1[i];

            t_0_xyy_0_y_0[i] = rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xyy_0_x_0[i] = rpb_x[i] * t_0_yy_0_x_0[i] + rwp_x[i] * t_0_yy_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_xxz_0_z_0[i] = rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxz_0_x_0[i] = rpb_z[i] * t_0_xx_0_x_0[i] + rwp_z[i] * t_0_xx_0_x_1[i];

            t_0_xxy_0_y_0[i] = rpb_y[i] * t_0_xx_0_y_0[i] + rwp_y[i] * t_0_xx_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxy_0_x_0[i] = rpb_y[i] * t_0_xx_0_x_0[i] + rwp_y[i] * t_0_xx_0_x_1[i];

            t_0_xxx_0_z_0[i] = rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] = rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] = rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
}

template <typename T>
auto
compHostVRRForSFSP_V2(      BufferHostXY<T>&      intsBufferSFSP,
                      const BufferHostX<int32_t>& intsIndexesSFSP0,
                      const BufferHostXY<T>&      intsBufferSPSP0,
                      const BufferHostX<int32_t>& intsIndexesSPSP0,
                      const BufferHostXY<T>&      intsBufferSPSP1,
                      const BufferHostX<int32_t>& intsIndexesSPSP1,
                      const BufferHostXY<T>&      intsBufferSDSS1,
                      const BufferHostX<int32_t>& intsIndexesSDSS1,
                      const BufferHostXY<T>&      intsBufferSDSP0,
                      const BufferHostX<int32_t>& intsIndexesSDSP0,
                      const BufferHostXY<T>&      intsBufferSDSP1,
                      const BufferHostX<int32_t>& intsIndexesSDSP1,
                      const T*                    osFactorsZeta,
                      const T*                    osFactorsBraZeta,
                      const BufferHostMY<T, 3>&   rDistancesPB,
                      const BufferHostMY<T, 3>&   rDistancesWP,
                      const T*                    osFactorsBraRhoZeta,
                      const bool                  useSummation,
                      const int32_t               nBatchPairs) -> void
{
    // set up Obara-Saika factors

    auto fze_0 = osFactorsZeta;

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

    // set up [SFSP]^(0) integral components

    t_0_zzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(0));

    t_0_zzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(1));

    t_0_zzz_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(2));

    t_0_yzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(3));

    t_0_yzz_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(4));

    t_0_yyz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(5));

    t_0_yyy_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(6));

    t_0_yyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(7));

    t_0_yyy_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(8));

    t_0_xzz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(9));

    t_0_xyy_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(10));

    t_0_xxz_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(11));

    t_0_xxx_0_z_0 = intsBufferSFSP0.data(intsIndexesSFSP0(12));

    t_0_xxx_0_y_0 = intsBufferSFSP0.data(intsIndexesSFSP0(13));

    t_0_xxx_0_x_0 = intsBufferSFSP0.data(intsIndexesSFSP0(14));

    // set up [SPSP]^(0) integral components

    t_0_z_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(0));

    t_0_z_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(1));

    t_0_z_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(2));

    t_0_y_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(3));

    t_0_y_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(4));

    t_0_y_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(5));

    t_0_x_0_z_0 = intsBufferSPSP0.data(intsIndexesSPSP0(6));

    t_0_x_0_y_0 = intsBufferSPSP0.data(intsIndexesSPSP0(7));

    t_0_x_0_x_0 = intsBufferSPSP0.data(intsIndexesSPSP0(8));

    // set up [SPSP]^(1) integral components

    t_0_z_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(0));

    t_0_z_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(1));

    t_0_z_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(2));

    t_0_y_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(3));

    t_0_y_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(4));

    t_0_y_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(5));

    t_0_x_0_z_1 = intsBufferSPSP1.data(intsIndexesSPSP1(6));

    t_0_x_0_y_1 = intsBufferSPSP1.data(intsIndexesSPSP1(7));

    t_0_x_0_x_1 = intsBufferSPSP1.data(intsIndexesSPSP1(8));

    // set up [SDSS]^(1) integral components

    t_0_zz_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(0));

    t_0_yy_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(1));

    t_0_xx_0_0_1 = intsBufferSDSS1.data(intsIndexesSDSS1(2));

    // set up [SDSP]^(0) integral components

    t_0_zz_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(0));

    t_0_zz_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(1));

    t_0_zz_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(2));

    t_0_yy_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(3));

    t_0_yy_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(4));

    t_0_yy_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(5));

    t_0_xx_0_z_0 = intsBufferSDSP0.data(intsIndexesSDSP0(6));

    t_0_xx_0_y_0 = intsBufferSDSP0.data(intsIndexesSDSP0(7));

    t_0_xx_0_x_0 = intsBufferSDSP0.data(intsIndexesSDSP0(8));

    // set up [SDSP]^(1) integral components

    t_0_zz_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(0));

    t_0_zz_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(1));

    t_0_zz_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(2));

    t_0_yy_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(3));

    t_0_yy_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(4));

    t_0_yy_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(5));

    t_0_xx_0_z_1 = intsBufferSDSP1.data(intsIndexesSDSP1(6));

    t_0_xx_0_y_1 = intsBufferSDSP1.data(intsIndexesSDSP1(7));

    t_0_xx_0_x_1 = intsBufferSDSP1.data(intsIndexesSDSP1(8));

    // set up scaling factors

    const auto fact_1_2 = static_cast<T>(1.0 / 2.0);

    if (useSummation)
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxz_0_z_0, t_0_xyy_0_y_0, t_0_xzz_0_z_0, t_0_y_0_x_0,\
                               t_0_y_0_x_1, t_0_y_0_y_0, t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1,\
                               t_0_yy_0_0_1, t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0,\
                               t_0_yy_0_y_1, t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0,\
                               t_0_yyy_0_y_0, t_0_yyy_0_z_0, t_0_yyz_0_z_0, t_0_yzz_0_y_0,\
                               t_0_yzz_0_z_0, t_0_z_0_x_0, t_0_z_0_x_1, t_0_z_0_y_0,\
                               t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1, t_0_zz_0_0_1,\
                               t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0, t_0_zz_0_y_1,\
                               t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0, t_0_zzz_0_y_0,\
                               t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] += rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] += rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] += rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] += rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] += rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yyz_0_z_0[i] += rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_z_0[i] += rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] += rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] += rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] += rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xyy_0_y_0[i] += rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xxz_0_z_0[i] += rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxx_0_z_0[i] += rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] += rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] += rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
    else
    {
        #pragma omp simd align(frz2_0, fz_0, fze_0, rpb_x, rpb_y, rpb_z, rwp_x, rwp_y,\
                               rwp_z, t_0_x_0_x_0, t_0_x_0_x_1, t_0_x_0_y_0, t_0_x_0_y_1,\
                               t_0_x_0_z_0, t_0_x_0_z_1, t_0_xx_0_0_1, t_0_xx_0_x_0,\
                               t_0_xx_0_x_1, t_0_xx_0_y_0, t_0_xx_0_y_1, t_0_xx_0_z_0,\
                               t_0_xx_0_z_1, t_0_xxx_0_x_0, t_0_xxx_0_y_0, t_0_xxx_0_z_0,\
                               t_0_xxz_0_z_0, t_0_xyy_0_y_0, t_0_xzz_0_z_0, t_0_y_0_x_0,\
                               t_0_y_0_x_1, t_0_y_0_y_0, t_0_y_0_y_1, t_0_y_0_z_0, t_0_y_0_z_1,\
                               t_0_yy_0_0_1, t_0_yy_0_x_0, t_0_yy_0_x_1, t_0_yy_0_y_0,\
                               t_0_yy_0_y_1, t_0_yy_0_z_0, t_0_yy_0_z_1, t_0_yyy_0_x_0,\
                               t_0_yyy_0_y_0, t_0_yyy_0_z_0, t_0_yyz_0_z_0, t_0_yzz_0_y_0,\
                               t_0_yzz_0_z_0, t_0_z_0_x_0, t_0_z_0_x_1, t_0_z_0_y_0,\
                               t_0_z_0_y_1, t_0_z_0_z_0, t_0_z_0_z_1, t_0_zz_0_0_1,\
                               t_0_zz_0_x_0, t_0_zz_0_x_1, t_0_zz_0_y_0, t_0_zz_0_y_1,\
                               t_0_zz_0_z_0, t_0_zz_0_z_1, t_0_zzz_0_x_0, t_0_zzz_0_y_0,\
                               t_0_zzz_0_z_0 : VLX_ALIGN)
        for (int32_t i = 0; i < nBatchPairs; i++)
        {
            t_0_zzz_0_z_0[i] = rpb_z[i] * t_0_zz_0_z_0[i] + rwp_z[i] * t_0_zz_0_z_1[i] + fz_0[i] * t_0_z_0_z_0[i] - frz2_0[i] * t_0_z_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_zzz_0_y_0[i] = rpb_z[i] * t_0_zz_0_y_0[i] + rwp_z[i] * t_0_zz_0_y_1[i] + fz_0[i] * t_0_z_0_y_0[i] - frz2_0[i] * t_0_z_0_y_1[i];

            t_0_zzz_0_x_0[i] = rpb_z[i] * t_0_zz_0_x_0[i] + rwp_z[i] * t_0_zz_0_x_1[i] + fz_0[i] * t_0_z_0_x_0[i] - frz2_0[i] * t_0_z_0_x_1[i];

            t_0_yzz_0_z_0[i] = rpb_y[i] * t_0_zz_0_z_0[i] + rwp_y[i] * t_0_zz_0_z_1[i];

            t_0_yzz_0_y_0[i] = rpb_y[i] * t_0_zz_0_y_0[i] + rwp_y[i] * t_0_zz_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_zz_0_0_1[i];

            t_0_yyz_0_z_0[i] = rpb_z[i] * t_0_yy_0_z_0[i] + rwp_z[i] * t_0_yy_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_z_0[i] = rpb_y[i] * t_0_yy_0_z_0[i] + rwp_y[i] * t_0_yy_0_z_1[i] + fz_0[i] * t_0_y_0_z_0[i] - frz2_0[i] * t_0_y_0_z_1[i];

            t_0_yyy_0_y_0[i] = rpb_y[i] * t_0_yy_0_y_0[i] + rwp_y[i] * t_0_yy_0_y_1[i] + fz_0[i] * t_0_y_0_y_0[i] - frz2_0[i] * t_0_y_0_y_1[i] + fact_1_2 * fze_0[i] * t_0_yy_0_0_1[i];

            t_0_yyy_0_x_0[i] = rpb_y[i] * t_0_yy_0_x_0[i] + rwp_y[i] * t_0_yy_0_x_1[i] + fz_0[i] * t_0_y_0_x_0[i] - frz2_0[i] * t_0_y_0_x_1[i];

            t_0_xzz_0_z_0[i] = rpb_x[i] * t_0_zz_0_z_0[i] + rwp_x[i] * t_0_zz_0_z_1[i];

            t_0_xyy_0_y_0[i] = rpb_x[i] * t_0_yy_0_y_0[i] + rwp_x[i] * t_0_yy_0_y_1[i];

            t_0_xxz_0_z_0[i] = rpb_z[i] * t_0_xx_0_z_0[i] + rwp_z[i] * t_0_xx_0_z_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];

            t_0_xxx_0_z_0[i] = rpb_x[i] * t_0_xx_0_z_0[i] + rwp_x[i] * t_0_xx_0_z_1[i] + fz_0[i] * t_0_x_0_z_0[i] - frz2_0[i] * t_0_x_0_z_1[i];

            t_0_xxx_0_y_0[i] = rpb_x[i] * t_0_xx_0_y_0[i] + rwp_x[i] * t_0_xx_0_y_1[i] + fz_0[i] * t_0_x_0_y_0[i] - frz2_0[i] * t_0_x_0_y_1[i];

            t_0_xxx_0_x_0[i] = rpb_x[i] * t_0_xx_0_x_0[i] + rwp_x[i] * t_0_xx_0_x_1[i] + fz_0[i] * t_0_x_0_x_0[i] - frz2_0[i] * t_0_x_0_x_1[i] + fact_1_2 * fze_0[i] * t_0_xx_0_0_1[i];
        }
    }
}


} // derirec namespace
