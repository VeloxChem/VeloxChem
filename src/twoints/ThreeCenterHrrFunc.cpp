//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ThreeCenterHrrFunc.hpp"

#include "GenFunc.hpp"
#include "AngularMomentum.hpp"

namespace t3hrrfunc { // t3hrrfunc namespace
    
    void
    compElectronRepulsionForXPP(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern
        
        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 1})) return;
        
        // determine number of components on bra side
        
        auto bcomp = angmom::to_SphericalComponents(braAngularMomentum);
        
        // determine number of contracted pairs on ket side
        
        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
        
        // set up pointers to distances R(CD) = C - D
        
        auto rcdx = cdDistances.data(0);
        
        auto rcdy = cdDistances.data(1);
        
        auto rcdz = cdDistances.data(2);
        
        // get position of integrals in integrals buffer
        
        auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                             {braAngularMomentum, 1, 1});
        
        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 2});
        
        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 1});
        
        // compute primitive integrals up to required order
        
        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SP) integrals
            
            auto g_0_x = contrBuffer.data(g2off + 3 * i);
            
            auto g_0_y = contrBuffer.data(g2off + 3 * i + 1);
            
            auto g_0_z = contrBuffer.data(g2off + 3 * i + 2);
            
            // set up pointers to (X|g(r,r')|SD) integrals
            
            auto g_0_xx = contrBuffer.data(g1off + 6 * i);
            
            auto g_0_xy = contrBuffer.data(g1off + 6 * i + 1);
            
            auto g_0_xz = contrBuffer.data(g1off + 6 * i + 2);
            
            auto g_0_yy = contrBuffer.data(g1off + 6 * i + 3);
            
            auto g_0_yz = contrBuffer.data(g1off + 6 * i + 4);
            
            auto g_0_zz = contrBuffer.data(g1off + 6 * i + 5);
            
            // set up pointers to (X|g(r,r')|PP) integrals
            
            auto g_x_x = contrBuffer.data(goff + 9 * i);
            
            auto g_x_y = contrBuffer.data(goff + 9 * i + 1);
            
            auto g_x_z = contrBuffer.data(goff + 9 * i + 2);
            
            auto g_y_x = contrBuffer.data(goff + 9 * i + 3);
            
            auto g_y_y = contrBuffer.data(goff + 9 * i + 4);
            
            auto g_y_z = contrBuffer.data(goff + 9 * i + 5);
            
            auto g_z_x = contrBuffer.data(goff + 9 * i + 6);
            
            auto g_z_y = contrBuffer.data(goff + 9 * i + 7);
            
            auto g_z_z = contrBuffer.data(goff + 9 * i + 8);
            
            #pragma omp simd aligned(rcdx, rcdy, rcdz, g_0_x, g_0_y, g_0_z,\
                                     g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz,\
                                     g_0_zz, g_x_x, g_x_y, g_x_z, g_y_x, g_y_y,\
                                     g_y_z, g_z_x, g_z_y, g_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component
                
                auto fr = rcdx[j];
                
                g_x_x[j] = g_0_xx[j] - fr * g_0_x[j];
                
                g_x_y[j] = g_0_xy[j] - fr * g_0_y[j];
                
                g_x_z[j] = g_0_xz[j] - fr * g_0_z[j];
                
                // leading y component
                
                fr = rcdy[j];
                
                g_y_x[j] = g_0_xy[j] - fr * g_0_x[j];
                
                g_y_y[j] = g_0_yy[j] - fr * g_0_y[j];
                
                g_y_z[j] = g_0_yz[j] - fr * g_0_z[j];
                
                // leading z component
                
                fr = rcdz[j];
                
                g_z_x[j] = g_0_xz[j] - fr * g_0_x[j];
                
                g_z_y[j] = g_0_yz[j] - fr * g_0_y[j];
                
                g_z_z[j] = g_0_zz[j] - fr * g_0_z[j];
            }
        }
    }
        
} // t3hrrfunc namespace
