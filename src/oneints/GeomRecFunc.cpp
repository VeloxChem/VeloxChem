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

#include "GeomRecFunc.hpp"

#include "AngularMomentum.hpp"

namespace geomrecfunc {  // geomrecfunc namespace

void
compGeomForSX(CMemBlock2D<double>&       primBuffer,
              const CMemBlock2D<double>& osFactors,
              const int32_t              nOSFactors,
              const int32_t              iPrimBuffSX,
              const int32_t              iPrimBuffPX,
              const CGtoBlock&           braGtoBlock,
              const CGtoBlock&           ketGtoBlock,
              const int32_t              iContrGto,
              const char                 axis)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();
    
    // set up number of Cartesian components on ket side
    
    const auto kang = ketGtoBlock.getAngularMomentum();
    
    const auto kcomps = angmom::to_CartesianComponents(kang);
    
    // loop over components
    
    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fe = osFactors.data(nOSFactors * idx + nOSFactors - 1);
        
        // loop over ket components
        
        for (int32_t j = 0; j < kcomps; j++)
        {
            if (axis == 'x')
            {
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                #pragma omp simd aligned(ts_0_j, tp_x_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    ts_0_j[k] = tp_x_j[k] * fe[k];
                }
            }
            
            if (axis == 'y')
            {
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                #pragma omp simd aligned(ts_0_j, tp_y_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    ts_0_j[k] = tp_y_j[k] * fe[k];
                }
            }
            
            if (axis == 'z')
            {
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                #pragma omp simd aligned(ts_0_j, tp_z_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    ts_0_j[k] = tp_z_j[k] * fe[k];
                }
            }
        }

        idx++;
    }
}

void
compGeomForPX(CMemBlock2D<double>&       primBuffer,
              const CMemBlock2D<double>& osFactors,
              const int32_t              nOSFactors,
              const int32_t              iPrimBuffPX,
              const int32_t              iPrimBuffDX,
              const int32_t              iPrimBuffSX,
              const CGtoBlock&           braGtoBlock,
              const CGtoBlock&           ketGtoBlock,
              const int32_t              iContrGto,
              const char                 axis)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();
    
    // set up number of Cartesian components on ket side
    
    const auto kang = ketGtoBlock.getAngularMomentum();
    
    const auto kcomps = angmom::to_CartesianComponents(kang);
    
    // loop over components
    
    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fe = osFactors.data(nOSFactors * idx + nOSFactors - 1);
        
        // loop over ket components
        
        for (int32_t j = 0; j < kcomps; j++)
        {
            if (axis == 'x')
            {
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                auto td_xx_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + j);
                
                auto td_xy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + kcomps + j);
                
                auto td_xz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 2 * kcomps + j);
                
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                #pragma omp simd aligned(ts_0_j, tp_x_j, tp_y_j, tp_z_j, td_xx_j, td_xy_j, td_xz_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    tp_x_j[k] = td_xx_j[k] * fe[k] - ts_0_j[k];
                    
                    tp_y_j[k] = td_xy_j[k] * fe[k];
                    
                    tp_z_j[k] = td_xz_j[k] * fe[k];
                }
            }
            
            if (axis == 'y')
            {
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                auto td_xy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + kcomps + j);
                
                auto td_yy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 3 * kcomps + j);
                
                auto td_yz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 4 * kcomps + j);
                
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                #pragma omp simd aligned(ts_0_j, tp_x_j, tp_y_j, tp_z_j, td_yx_j, td_yy_j, td_yz_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    tp_x_j[k] = td_xy_j[k] * fe[k];
                    
                    tp_y_j[k] = td_yy_j[k] * fe[k] - ts_0_j[k];
                    
                    tp_z_j[k] = td_yz_j[k] * fe[k];
                }
            }
            
            if (axis == 'z')
            {
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                auto td_zx_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 2 * kcomps + j);
                
                auto td_zy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 4 * kcomps + j);
                
                auto td_zz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 5 * kcomps + j);
                
                auto ts_0_j = primBuffer.data(iPrimBuffSX + kcomps * idx + j);
                
                #pragma omp simd aligned(ts_0_j, tp_x_j, tp_y_j, tp_z_j, td_zx_j, td_zy_j, td_zz_j, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    tp_x_j[k] = td_zx_j[k] * fe[k];
                    
                    tp_y_j[k] = td_zy_j[k] * fe[k];
                    
                    tp_z_j[k] = td_zz_j[k] * fe[k] - ts_0_j[k];
                }
            }
        }

        idx++;
    }
}

void
compGeomForDX(CMemBlock2D<double>&       primBuffer,
              const CMemBlock2D<double>& osFactors,
              const int32_t              nOSFactors,
              const int32_t              iPrimBuffDX,
              const int32_t              iPrimBuffFX,
              const int32_t              iPrimBuffPX,
              const CGtoBlock&           braGtoBlock,
              const CGtoBlock&           ketGtoBlock,
              const int32_t              iContrGto,
              const char                 axis)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();
    
    // set up number of Cartesian components on ket side
    
    const auto kang = ketGtoBlock.getAngularMomentum();
    
    const auto kcomps = angmom::to_CartesianComponents(kang);
    
    // loop over components
    
    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fe = osFactors.data(nOSFactors * idx + nOSFactors - 1);
        
        // loop over ket components
        
        for (int32_t j = 0; j < kcomps; j++)
        {
            if (axis == 'x')
            {
                auto td_xx_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + j);
                
                auto td_xy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + kcomps + j);
                
                auto td_xz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 2 * kcomps + j);
                
                auto td_yy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 3 * kcomps + j);
                
                auto td_yz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 4 * kcomps + j);
                
                auto td_zz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 5 * kcomps + j);
                
                auto tf_xxx_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + j);
                
                auto tf_xxy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + kcomps + j);
                
                auto tf_xxz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 2 * kcomps + j);
                
                auto tf_xyy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 3 * kcomps + j);
                
                auto tf_xyz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 4 * kcomps + j);
                
                auto tf_xzz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 5 * kcomps + j);
                
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                #pragma omp simd aligned(tp_x_j, tp_y_j, tp_z_j, td_xx_j, td_xy_j, td_xz_j,\
                                         td_yy_j, td_yz_j, td_zz_j, tf_xxx_j, tf_xxy_j,\
                                         tf_xxz_j, tf_xyy_j, tf_xyz_j, tf_xzz, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    td_xx_j[k] = tf_xxx_j[k] * fe[k] - 4.0 * tp_x_j[k];
                    
                    td_xy_j[k] = tf_xxy_j[k] * fe[k] - 2.0 * tp_y_j[k];
                    
                    td_xz_j[k] = tf_xxz_j[k] * fe[k] - 2.0 * tp_z_j[k];
                    
                    td_yy_j[k] = tf_xyy_j[k] * fe[k];
                    
                    td_yz_j[k] = tf_xyz_j[k] * fe[k];
                    
                    td_zz_j[k] = tf_xzz_j[k] * fe[k];
                }
            }
            
            if (axis == 'y')
            {
                auto td_xx_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + j);
                
                auto td_xy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + kcomps + j);
                
                auto td_xz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 2 * kcomps + j);
                
                auto td_yy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 3 * kcomps + j);
                
                auto td_yz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 4 * kcomps + j);
                
                auto td_zz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 5 * kcomps + j);
                
                auto tf_yxx_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + kcomps + j);
                
                auto tf_yxy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 3 * kcomps + j);
                
                auto tf_yxz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 4 * kcomps + j);
                
                auto tf_yyy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 6 * kcomps + j);
                
                auto tf_yyz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 7 * kcomps + j);
                
                auto tf_yzz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 8 * kcomps + j);
                
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                #pragma omp simd aligned(tp_x_j, tp_y_j, tp_z_j, td_xx_j, td_xy_j, td_xz_j,\
                                         td_yy_j, td_yz_j, td_zz_j, tf_yxx_j, tf_yxy_j,\
                                         tf_yxz_j, tf_yyy_j, tf_yyz_j, tf_yzz, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    td_xx_j[k] = tf_yxx_j[k] * fe[k];
                    
                    td_xy_j[k] = tf_yxy_j[k] * fe[k] - 2.0 * tp_x_j[k];
                    
                    td_xz_j[k] = tf_yxz_j[k] * fe[k];
                    
                    td_yy_j[k] = tf_yyy_j[k] * fe[k] - 4.0 * tp_y_j[k];
                    
                    td_yz_j[k] = tf_yyz_j[k] * fe[k] - 2.0 * tp_z_j[k];
                    
                    td_zz_j[k] = tf_yzz_j[k] * fe[k];
                }
            }
            
            if (axis == 'z')
            {
                auto td_xx_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + j);
                
                auto td_xy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + kcomps + j);
                
                auto td_xz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 2 * kcomps + j);
                
                auto td_yy_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 3 * kcomps + j);
                
                auto td_yz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 4 * kcomps + j);
                
                auto td_zz_j = primBuffer.data(iPrimBuffDX + 6 * kcomps * idx + 5 * kcomps + j);
                
                auto tf_zxx_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 2 * kcomps + j);
                
                auto tf_zxy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 4 * kcomps + j);
                
                auto tf_zxz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 5 * kcomps + j);
                
                auto tf_zyy_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 7 * kcomps + j);
                
                auto tf_zyz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 8 * kcomps + j);
                
                auto tf_zzz_j = primBuffer.data(iPrimBuffFX + 10 * kcomps * idx + 9 * kcomps + j);
                
                auto tp_x_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + j);
                
                auto tp_y_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + kcomps + j);
                
                auto tp_z_j = primBuffer.data(iPrimBuffPX + 3 * kcomps * idx + 2 * kcomps + j);
                
                #pragma omp simd aligned(tp_x_j, tp_y_j, tp_z_j, td_xx_j, td_xy_j, td_xz_j,\
                                         td_yy_j, td_yz_j, td_zz_j, tf_zxx_j, tf_zxy_j,\
                                         tf_zxz_j, tf_zyy_j, tf_zyz_j, tf_zzz, fe: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    td_xx_j[k] = tf_zxx_j[k] * fe[k];
                    
                    td_xy_j[k] = tf_zxy_j[k] * fe[k];
                    
                    td_xz_j[k] = tf_zxz_j[k] * fe[k] - 2.0 * tp_x_j[k];
                    
                    td_yy_j[k] = tf_zyy_j[k] * fe[k];
                    
                    td_yz_j[k] = tf_zyz_j[k] * fe[k] - 2.0 * tp_y_j[k];
                    
                    td_zz_j[k] = tf_zzz_j[k] * fe[k] - 4.0 * tp_z_j[k];
                }
            }
        }

        idx++;
    }
}

}  // namespace geomrecfunc
