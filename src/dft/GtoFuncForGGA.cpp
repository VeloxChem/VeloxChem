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

#include "GtoFuncForGGA.hpp"

#include <cmath>

#include "GenFunc.hpp"

namespace ggarec {  // ggarec namespace
    
    void
    compGtoValuesForS(      CMemBlock2D<double>& spherGtoGridBuffer,
                      const double*              gridCoordinatesX,
                      const double*              gridCoordinatesY,
                      const double*              gridCoordinatesZ,
                      const int32_t              gridOffset,
                      const CGtoBlock&           gtoBlock,
                      const int32_t              iContrGto)
    {
        // initialize buffer to zero
        
        spherGtoGridBuffer.zero();
        
        // set up number of grid points
        
        auto ngpnts = spherGtoGridBuffer.size(0);
        
        // set up pointers to primitives data
        
        auto bfnorms = gtoBlock.getNormFactors();
        
        auto bfexps = gtoBlock.getExponents();
        
        // set up pointers to primitives coordinates
        
        auto bfx = gtoBlock.getCoordinatesX();
        
        auto bfy = gtoBlock.getCoordinatesY();
        
        auto bfz = gtoBlock.getCoordinatesZ();
        
        // set up coordinates to primitives positions
        
        auto spos = gtoBlock.getStartPositions();
        
        auto epos = gtoBlock.getEndPositions();
        
        // set up pointer to spherical GTOs values
        
        auto f0_0 = spherGtoGridBuffer.data(0);
        
        auto fx_0 = spherGtoGridBuffer.data(1);
        
        auto fy_0 = spherGtoGridBuffer.data(2);
        
        auto fz_0 = spherGtoGridBuffer.data(3);
        
        // initialize Cartesian buffer to zero
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // primitive data
            
            auto bexp = bfexps[i];
            
            auto bnorm = bfnorms[i];
            
            // primitive position
            
            auto rx = bfx[i];
            
            auto ry = bfy[i];
            
            auto rz = bfz[i];
            
            // loop over grid points
            
            #pragma omp simd
            for (int32_t j = 0; j < ngpnts; j++)
            {
                double dx = gridCoordinatesX[gridOffset +  j] - rx;
                
                double dy = gridCoordinatesY[gridOffset +  j] - ry;
                
                double dz = gridCoordinatesZ[gridOffset +  j] - rz;
                
                double g0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                
                double g1 = -2.0 * bexp * g0;
                
                f0_0[j] += g0;
                
                fx_0[j] += dx * g1;
                
                fy_0[j] += dy * g1;
                
                fz_0[j] += dz * g1;
            }
        }
    }
    
    void
    compGtoValuesForP(      CMemBlock2D<double>& spherGtoGridBuffer,
                            CMemBlock2D<double>& cartGtoGridBuffer,
                      const double*              gridCoordinatesX,
                      const double*              gridCoordinatesY,
                      const double*              gridCoordinatesZ,
                      const int32_t              gridOffset,
                      const CGtoBlock&           gtoBlock,
                      const int32_t              iContrGto)
    {
        // initialize buffers to zero
        
        spherGtoGridBuffer.zero();
        
        cartGtoGridBuffer.zero();
        
        // set up number of grid points
        
        auto ngpnts = spherGtoGridBuffer.size(0);
        
        // set up pointers to primitives data
        
        auto bfnorms = gtoBlock.getNormFactors();
        
        auto bfexps = gtoBlock.getExponents();
        
        // set up pointers to primitives coordinates
        
        auto bfx = gtoBlock.getCoordinatesX();
        
        auto bfy = gtoBlock.getCoordinatesY();
        
        auto bfz = gtoBlock.getCoordinatesZ();
        
        // set up coordinates to primitives positions
        
        auto spos = gtoBlock.getStartPositions();
        
        auto epos = gtoBlock.getEndPositions();
        
        // set up pointer to Cartesian GTOs values
        
        auto f0_x = cartGtoGridBuffer.data(0);
        
        auto fx_x = cartGtoGridBuffer.data(1);
        
        auto fy_x = cartGtoGridBuffer.data(2);
        
        auto fz_x = cartGtoGridBuffer.data(3);
        
        auto f0_y = cartGtoGridBuffer.data(4);
        
        auto fx_y = cartGtoGridBuffer.data(5);
        
        auto fy_y = cartGtoGridBuffer.data(6);
        
        auto fz_y = cartGtoGridBuffer.data(7);
        
        auto f0_z = cartGtoGridBuffer.data(8);
        
        auto fx_z = cartGtoGridBuffer.data(9);
        
        auto fy_z = cartGtoGridBuffer.data(10);
        
        auto fz_z = cartGtoGridBuffer.data(11);
        
        // initialize Cartesian buffer to zero
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // primitive data
            
            auto bexp = bfexps[i];
            
            auto bnorm = bfnorms[i];
            
            // primitive position
            
            auto rx = bfx[i];
            
            auto ry = bfy[i];
            
            auto rz = bfz[i];
            
            // loop over grid points
            
            #pragma omp simd
            for (int32_t j = 0; j < ngpnts; j++)
            {
                double dx = gridCoordinatesX[gridOffset +  j] - rx;
                
                double dy = gridCoordinatesY[gridOffset +  j] - ry;
                
                double dz = gridCoordinatesZ[gridOffset +  j] - rz;
                
                double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                
                double fg_0 = -2.0 * bexp;
                
                // leading p_x
                
                f0_x[j] += f0_0 * dx;
                
                fx_x[j] += f0_0 * (1.0 + dx * dx * fg_0);
                
                fy_x[j] += f0_0 * dx * dy * fg_0;
                
                fz_x[j] += f0_0 * dx * dz * fg_0;
                
                // leading p_y
                
                f0_y[j] += f0_0 * dy;
                
                fx_y[j] += f0_0 * dy * dx * fg_0;
                
                fy_y[j] += f0_0 * (1.0 + dy * dy * fg_0);
                
                fz_y[j] += f0_0 * dy * dz * fg_0;
                
                // leading p_z
                
                f0_z[j] += f0_0 * dz;
                
                fx_z[j] += f0_0 * dz * dx * fg_0;
                
                fy_z[j] += f0_0 * dz * dy * fg_0;
                
                fz_z[j] += f0_0 * (1.0 + dz * dz * fg_0);
            }
        }
        
        genfunc::transform(spherGtoGridBuffer, cartGtoGridBuffer, CSphericalMomentum(1), 0, 0, ngpnts, 4);
    }
    
    void
    compGtoValuesForD(      CMemBlock2D<double>& spherGtoGridBuffer,
                            CMemBlock2D<double>& cartGtoGridBuffer,
                      const double*              gridCoordinatesX,
                      const double*              gridCoordinatesY,
                      const double*              gridCoordinatesZ,
                      const int32_t              gridOffset,
                      const CGtoBlock&           gtoBlock,
                      const int32_t              iContrGto)
    {
        // initialize buffer to zero
        
        spherGtoGridBuffer.zero();
        
        cartGtoGridBuffer.zero();
        
        // set up number of grid points
        
        auto ngpnts = spherGtoGridBuffer.size(0);
        
        // set up pointers to primitives data
        
        auto bfnorms = gtoBlock.getNormFactors();
        
        auto bfexps = gtoBlock.getExponents();
        
        // set up pointers to primitives coordinates
        
        auto bfx = gtoBlock.getCoordinatesX();
        
        auto bfy = gtoBlock.getCoordinatesY();
        
        auto bfz = gtoBlock.getCoordinatesZ();
        
        // set up coordinates to primitives positions
        
        auto spos = gtoBlock.getStartPositions();
        
        auto epos = gtoBlock.getEndPositions();
        
        // set up pointer to spherical GTOs values
        
        auto f0_xx = cartGtoGridBuffer.data(0);
        
        auto fx_xx = cartGtoGridBuffer.data(1);
        
        auto fy_xx = cartGtoGridBuffer.data(2);
        
        auto fz_xx = cartGtoGridBuffer.data(3);
        
        auto f0_xy = cartGtoGridBuffer.data(4);
        
        auto fx_xy = cartGtoGridBuffer.data(5);
        
        auto fy_xy = cartGtoGridBuffer.data(6);
        
        auto fz_xy = cartGtoGridBuffer.data(7);
        
        auto f0_xz = cartGtoGridBuffer.data(8);
        
        auto fx_xz = cartGtoGridBuffer.data(9);
        
        auto fy_xz = cartGtoGridBuffer.data(10);
        
        auto fz_xz = cartGtoGridBuffer.data(11);
        
        auto f0_yy = cartGtoGridBuffer.data(12);
        
        auto fx_yy = cartGtoGridBuffer.data(13);
        
        auto fy_yy = cartGtoGridBuffer.data(14);
        
        auto fz_yy = cartGtoGridBuffer.data(15);
        
        auto f0_yz = cartGtoGridBuffer.data(16);
        
        auto fx_yz = cartGtoGridBuffer.data(17);
        
        auto fy_yz = cartGtoGridBuffer.data(18);
        
        auto fz_yz = cartGtoGridBuffer.data(19);
        
        auto f0_zz = cartGtoGridBuffer.data(20);
        
        auto fx_zz = cartGtoGridBuffer.data(21);
        
        auto fy_zz = cartGtoGridBuffer.data(22);
        
        auto fz_zz = cartGtoGridBuffer.data(23);
        
        // initialize Cartesian buffer to zero
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // primitive data
            
            auto bexp = bfexps[i];
            
            auto bnorm = bfnorms[i];
            
            // primitive position
            
            auto rx = bfx[i];
            
            auto ry = bfy[i];
            
            auto rz = bfz[i];
            
            // loop over grid points
            
            #pragma omp simd
            for (int32_t j = 0; j < ngpnts; j++)
            {
                double dx = gridCoordinatesX[gridOffset +  j] - rx;
                
                double dy = gridCoordinatesY[gridOffset +  j] - ry;
                
                double dz = gridCoordinatesZ[gridOffset +  j] - rz;
                
                double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                
                double fg_0 = -2.0 * bexp;
                
                // leading xx component
                
                double f0_x = dx * f0_0;
                
                f0_xx[j] += f0_x * dx;
                
                fx_xx[j] += f0_x * (2.0 + fg_0 * dx * dx);
                
                fy_xx[j] += f0_x * fg_0 * dx * dy;
                
                fz_xx[j] += f0_x * fg_0 * dx * dz;
                
                // leading xy component
                
                double f0_y = dy * f0_0;
                
                f0_xy[j] += f0_x * dy;
                
                fx_xy[j] += f0_y * (1.0 + fg_0 * dx * dx);
                
                fy_xy[j] += f0_x * (1.0 + fg_0 * dy * dy);
                
                fz_xy[j] += f0_x * fg_0 * dy * dz;
                
                // leading xz component
                
                double f0_z = dz * f0_0;
                
                f0_xz[j] += f0_x * dz;
                
                fx_xz[j] += f0_z * (1.0 + fg_0 * dx * dx);
                
                fy_xz[j] += f0_x * fg_0 * dz * dy;
                
                fz_xz[j] += f0_x * (1.0 + fg_0 * dz * dz);
                
                // leading yy component
                
                f0_yy[j] += f0_y * dy;
                
                fx_yy[j] += f0_y * fg_0 * dy * dx;
                
                fy_yy[j] += f0_y * (2.0 + fg_0 * dy * dy);
                
                fz_yy[j] += f0_y * fg_0 * dy * dz;
                
                // leading yz component
                
                f0_yz[j] += f0_y * dz;
                
                fx_yz[j] += f0_y * fg_0 * dz * dx;
                
                fy_yz[j] += f0_z * (1.0 + fg_0 * dy * dy);
                
                fz_yz[j] += f0_y * (1.0 + fg_0 * dz * dz);
                
                // leading zz component
                
                f0_zz[j] += f0_z * dz;
                
                fx_zz[j] += f0_z * fg_0 * dz * dx;
                
                fy_zz[j] += f0_z * fg_0 * dz * dy;
                
                fz_zz[j] += f0_z * (2.0 + fg_0 * dz * dz);
            }
        }
        
        genfunc::transform(spherGtoGridBuffer, cartGtoGridBuffer, CSphericalMomentum(2), 0, 0, ngpnts, 4);
    }
    
}  // namespace ggarec
