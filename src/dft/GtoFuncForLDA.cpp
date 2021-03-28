//
//                           VELOXCHEM 1.0-RC
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

#include "GtoFuncForLDA.hpp"

#include <cmath>

#include "GenFunc.hpp"

namespace ldarec {  // ldarec namespace
    
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
        
        auto fs_0 = spherGtoGridBuffer.data(0);
        
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
                
                fs_0[j] += bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
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
        
        auto fp_x = cartGtoGridBuffer.data(0);
        
        auto fp_y = cartGtoGridBuffer.data(1);
        
        auto fp_z = cartGtoGridBuffer.data(2);
        
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
                
                double fs_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                
                fp_x[j] += fs_0 * dx;
                
                fp_y[j] += fs_0 * dy;
                
                fp_z[j] += fs_0 * dz;
            }
        }
        
        genfunc::transform(spherGtoGridBuffer, cartGtoGridBuffer, CSphericalMomentum(1), 0, 0, ngpnts, 1);
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
        
        auto fd_xx = cartGtoGridBuffer.data(0);
        
        auto fd_xy = cartGtoGridBuffer.data(1);
        
        auto fd_xz = cartGtoGridBuffer.data(2);
        
        auto fd_yy = cartGtoGridBuffer.data(3);
        
        auto fd_yz = cartGtoGridBuffer.data(4);
        
        auto fd_zz = cartGtoGridBuffer.data(5);
        
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
                
                double fs_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                
                // leading x component
                
                double fr = dx * fs_0;
                
                fd_xx[j] += fr * dx;
                
                fd_xy[j] += fr * dy;
                
                fd_xz[j] += fr * dz;
                
                // leading y component
                
                fr = dy * fs_0;
                
                fd_yy[j] += fr * dy;
                
                fd_yz[j] += fr * dz;
                
                // leading z component
                
                fd_zz[j] += dz * fs_0 * dz;
            }
        }
        
        genfunc::transform(spherGtoGridBuffer, cartGtoGridBuffer, CSphericalMomentum(2), 0, 0, ngpnts, 1);
    }
    
}  // namespace ldarec
