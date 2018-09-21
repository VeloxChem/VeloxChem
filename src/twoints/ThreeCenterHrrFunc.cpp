//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
        
        // compute contracted integrals
        
        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SP) integrals
            
            auto g2_0_x = contrBuffer.data(g2off + 3 * i);
            
            auto g2_0_y = contrBuffer.data(g2off + 3 * i + 1);
            
            auto g2_0_z = contrBuffer.data(g2off + 3 * i + 2);
            
            // set up pointers to (X|g(r,r')|SD) integrals
            
            auto g1_0_xx = contrBuffer.data(g1off + 6 * i);
            
            auto g1_0_xy = contrBuffer.data(g1off + 6 * i + 1);
            
            auto g1_0_xz = contrBuffer.data(g1off + 6 * i + 2);
            
            auto g1_0_yy = contrBuffer.data(g1off + 6 * i + 3);
            
            auto g1_0_yz = contrBuffer.data(g1off + 6 * i + 4);
            
            auto g1_0_zz = contrBuffer.data(g1off + 6 * i + 5);
            
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
            
            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_x, g2_0_y, g2_0_z,\
                                     g1_0_xx, g1_0_xy, g1_0_xz, g1_0_yy, g1_0_yz,\
                                     g1_0_zz, g_x_x, g_x_y, g_x_z, g_y_x, g_y_y,\
                                     g_y_z, g_z_x, g_z_y, g_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component
                
                auto fr = rcdx[j];
                
                g_x_x[j] = g1_0_xx[j] - fr * g2_0_x[j];
                
                g_x_y[j] = g1_0_xy[j] - fr * g2_0_y[j];
                
                g_x_z[j] = g1_0_xz[j] - fr * g2_0_z[j];
                
                // leading y component
                
                fr = rcdy[j];
                
                g_y_x[j] = g1_0_xy[j] - fr * g2_0_x[j];
                
                g_y_y[j] = g1_0_yy[j] - fr * g2_0_y[j];
                
                g_y_z[j] = g1_0_yz[j] - fr * g2_0_z[j];
                
                // leading z component
                
                fr = rcdz[j];
                
                g_z_x[j] = g1_0_xz[j] - fr * g2_0_x[j];
                
                g_z_y[j] = g1_0_yz[j] - fr * g2_0_y[j];
                
                g_z_z[j] = g1_0_zz[j] - fr * g2_0_z[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPD(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 2})) return;

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
                                             {braAngularMomentum, 1, 2});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 3});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 2});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SD)^(m) integrals

            auto g2_0_xx = contrBuffer.data(g2off + 6 * i);

            auto g2_0_xy = contrBuffer.data(g2off + 6 * i + 1);

            auto g2_0_xz = contrBuffer.data(g2off + 6 * i + 2);

            auto g2_0_yy = contrBuffer.data(g2off + 6 * i + 3);

            auto g2_0_yz = contrBuffer.data(g2off + 6 * i + 4);

            auto g2_0_zz = contrBuffer.data(g2off + 6 * i + 5);

            // set up pointers to (X|g(r,r')|SF)^(m) integrals

            auto g1_0_xxx = contrBuffer.data(g1off + 10 * i);

            auto g1_0_xxy = contrBuffer.data(g1off + 10 * i + 1);

            auto g1_0_xxz = contrBuffer.data(g1off + 10 * i + 2);

            auto g1_0_xyy = contrBuffer.data(g1off + 10 * i + 3);

            auto g1_0_xyz = contrBuffer.data(g1off + 10 * i + 4);

            auto g1_0_xzz = contrBuffer.data(g1off + 10 * i + 5);

            auto g1_0_yyy = contrBuffer.data(g1off + 10 * i + 6);

            auto g1_0_yyz = contrBuffer.data(g1off + 10 * i + 7);

            auto g1_0_yzz = contrBuffer.data(g1off + 10 * i + 8);

            auto g1_0_zzz = contrBuffer.data(g1off + 10 * i + 9);

            // set up pointers to (X|g(r,r')|PD)^(m) integrals

            auto g_x_xx = contrBuffer.data(goff + 18 * i);

            auto g_x_xy = contrBuffer.data(goff + 18 * i + 1);

            auto g_x_xz = contrBuffer.data(goff + 18 * i + 2);

            auto g_x_yy = contrBuffer.data(goff + 18 * i + 3);

            auto g_x_yz = contrBuffer.data(goff + 18 * i + 4);

            auto g_x_zz = contrBuffer.data(goff + 18 * i + 5);

            auto g_y_xx = contrBuffer.data(goff + 18 * i + 6);

            auto g_y_xy = contrBuffer.data(goff + 18 * i + 7);

            auto g_y_xz = contrBuffer.data(goff + 18 * i + 8);

            auto g_y_yy = contrBuffer.data(goff + 18 * i + 9);

            auto g_y_yz = contrBuffer.data(goff + 18 * i + 10);

            auto g_y_zz = contrBuffer.data(goff + 18 * i + 11);

            auto g_z_xx = contrBuffer.data(goff + 18 * i + 12);

            auto g_z_xy = contrBuffer.data(goff + 18 * i + 13);

            auto g_z_xz = contrBuffer.data(goff + 18 * i + 14);

            auto g_z_yy = contrBuffer.data(goff + 18 * i + 15);

            auto g_z_yz = contrBuffer.data(goff + 18 * i + 16);

            auto g_z_zz = contrBuffer.data(goff + 18 * i + 17);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xx, g2_0_xy, g2_0_xz,\
                                     g2_0_yy, g2_0_yz, g2_0_zz, g1_0_xxx, g1_0_xxy,\
                                     g1_0_xxz, g1_0_xyy, g1_0_xyz, g1_0_xzz, g1_0_yyy,\
                                     g1_0_yyz, g1_0_yzz, g1_0_zzz, g_x_xx, g_x_xy,\
                                     g_x_xz, g_x_yy, g_x_yz, g_x_zz, g_y_xx, g_y_xy,\
                                     g_y_xz, g_y_yy, g_y_yz, g_y_zz, g_z_xx, g_z_xy,\
                                     g_z_xz, g_z_yy, g_z_yz, g_z_zz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xx[j] = g1_0_xxx[j] - fr * g2_0_xx[j];

                g_x_xy[j] = g1_0_xxy[j] - fr * g2_0_xy[j];

                g_x_xz[j] = g1_0_xxz[j] - fr * g2_0_xz[j];

                g_x_yy[j] = g1_0_xyy[j] - fr * g2_0_yy[j];

                g_x_yz[j] = g1_0_xyz[j] - fr * g2_0_yz[j];

                g_x_zz[j] = g1_0_xzz[j] - fr * g2_0_zz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xx[j] = g1_0_xxy[j] - fr * g2_0_xx[j];

                g_y_xy[j] = g1_0_xyy[j] - fr * g2_0_xy[j];

                g_y_xz[j] = g1_0_xyz[j] - fr * g2_0_xz[j];

                g_y_yy[j] = g1_0_yyy[j] - fr * g2_0_yy[j];

                g_y_yz[j] = g1_0_yyz[j] - fr * g2_0_yz[j];

                g_y_zz[j] = g1_0_yzz[j] - fr * g2_0_zz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xx[j] = g1_0_xxz[j] - fr * g2_0_xx[j];

                g_z_xy[j] = g1_0_xyz[j] - fr * g2_0_xy[j];

                g_z_xz[j] = g1_0_xzz[j] - fr * g2_0_xz[j];

                g_z_yy[j] = g1_0_yyz[j] - fr * g2_0_yy[j];

                g_z_yz[j] = g1_0_yzz[j] - fr * g2_0_yz[j];

                g_z_zz[j] = g1_0_zzz[j] - fr * g2_0_zz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPF(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 3})) return;

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
                                             {braAngularMomentum, 1, 3});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 4});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 3});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SF)^(m) integrals

            auto g2_0_xxx = contrBuffer.data(g2off + 10 * i);

            auto g2_0_xxy = contrBuffer.data(g2off + 10 * i + 1);

            auto g2_0_xxz = contrBuffer.data(g2off + 10 * i + 2);

            auto g2_0_xyy = contrBuffer.data(g2off + 10 * i + 3);

            auto g2_0_xyz = contrBuffer.data(g2off + 10 * i + 4);

            auto g2_0_xzz = contrBuffer.data(g2off + 10 * i + 5);

            auto g2_0_yyy = contrBuffer.data(g2off + 10 * i + 6);

            auto g2_0_yyz = contrBuffer.data(g2off + 10 * i + 7);

            auto g2_0_yzz = contrBuffer.data(g2off + 10 * i + 8);

            auto g2_0_zzz = contrBuffer.data(g2off + 10 * i + 9);

            // set up pointers to (X|g(r,r')|SG)^(m) integrals

            auto g1_0_xxxx = contrBuffer.data(g1off + 15 * i);

            auto g1_0_xxxy = contrBuffer.data(g1off + 15 * i + 1);

            auto g1_0_xxxz = contrBuffer.data(g1off + 15 * i + 2);

            auto g1_0_xxyy = contrBuffer.data(g1off + 15 * i + 3);

            auto g1_0_xxyz = contrBuffer.data(g1off + 15 * i + 4);

            auto g1_0_xxzz = contrBuffer.data(g1off + 15 * i + 5);

            auto g1_0_xyyy = contrBuffer.data(g1off + 15 * i + 6);

            auto g1_0_xyyz = contrBuffer.data(g1off + 15 * i + 7);

            auto g1_0_xyzz = contrBuffer.data(g1off + 15 * i + 8);

            auto g1_0_xzzz = contrBuffer.data(g1off + 15 * i + 9);

            auto g1_0_yyyy = contrBuffer.data(g1off + 15 * i + 10);

            auto g1_0_yyyz = contrBuffer.data(g1off + 15 * i + 11);

            auto g1_0_yyzz = contrBuffer.data(g1off + 15 * i + 12);

            auto g1_0_yzzz = contrBuffer.data(g1off + 15 * i + 13);

            auto g1_0_zzzz = contrBuffer.data(g1off + 15 * i + 14);

            // set up pointers to (X|g(r,r')|PF)^(m) integrals

            auto g_x_xxx = contrBuffer.data(goff + 30 * i);

            auto g_x_xxy = contrBuffer.data(goff + 30 * i + 1);

            auto g_x_xxz = contrBuffer.data(goff + 30 * i + 2);

            auto g_x_xyy = contrBuffer.data(goff + 30 * i + 3);

            auto g_x_xyz = contrBuffer.data(goff + 30 * i + 4);

            auto g_x_xzz = contrBuffer.data(goff + 30 * i + 5);

            auto g_x_yyy = contrBuffer.data(goff + 30 * i + 6);

            auto g_x_yyz = contrBuffer.data(goff + 30 * i + 7);

            auto g_x_yzz = contrBuffer.data(goff + 30 * i + 8);

            auto g_x_zzz = contrBuffer.data(goff + 30 * i + 9);

            auto g_y_xxx = contrBuffer.data(goff + 30 * i + 10);

            auto g_y_xxy = contrBuffer.data(goff + 30 * i + 11);

            auto g_y_xxz = contrBuffer.data(goff + 30 * i + 12);

            auto g_y_xyy = contrBuffer.data(goff + 30 * i + 13);

            auto g_y_xyz = contrBuffer.data(goff + 30 * i + 14);

            auto g_y_xzz = contrBuffer.data(goff + 30 * i + 15);

            auto g_y_yyy = contrBuffer.data(goff + 30 * i + 16);

            auto g_y_yyz = contrBuffer.data(goff + 30 * i + 17);

            auto g_y_yzz = contrBuffer.data(goff + 30 * i + 18);

            auto g_y_zzz = contrBuffer.data(goff + 30 * i + 19);

            auto g_z_xxx = contrBuffer.data(goff + 30 * i + 20);

            auto g_z_xxy = contrBuffer.data(goff + 30 * i + 21);

            auto g_z_xxz = contrBuffer.data(goff + 30 * i + 22);

            auto g_z_xyy = contrBuffer.data(goff + 30 * i + 23);

            auto g_z_xyz = contrBuffer.data(goff + 30 * i + 24);

            auto g_z_xzz = contrBuffer.data(goff + 30 * i + 25);

            auto g_z_yyy = contrBuffer.data(goff + 30 * i + 26);

            auto g_z_yyz = contrBuffer.data(goff + 30 * i + 27);

            auto g_z_yzz = contrBuffer.data(goff + 30 * i + 28);

            auto g_z_zzz = contrBuffer.data(goff + 30 * i + 29);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxx, g2_0_xxy, g2_0_xxz,\
                                     g2_0_xyy, g2_0_xyz, g2_0_xzz, g2_0_yyy, g2_0_yyz,\
                                     g2_0_yzz, g2_0_zzz, g1_0_xxxx, g1_0_xxxy,\
                                     g1_0_xxxz, g1_0_xxyy, g1_0_xxyz, g1_0_xxzz,\
                                     g1_0_xyyy, g1_0_xyyz, g1_0_xyzz, g1_0_xzzz,\
                                     g1_0_yyyy, g1_0_yyyz, g1_0_yyzz, g1_0_yzzz,\
                                     g1_0_zzzz, g_x_xxx, g_x_xxy, g_x_xxz, g_x_xyy,\
                                     g_x_xyz, g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz,\
                                     g_x_zzz, g_y_xxx, g_y_xxy, g_y_xxz, g_y_xyy,\
                                     g_y_xyz, g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz,\
                                     g_y_zzz, g_z_xxx, g_z_xxy, g_z_xxz, g_z_xyy,\
                                     g_z_xyz, g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz,\
                                     g_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xxx[j] = g1_0_xxxx[j] - fr * g2_0_xxx[j];

                g_x_xxy[j] = g1_0_xxxy[j] - fr * g2_0_xxy[j];

                g_x_xxz[j] = g1_0_xxxz[j] - fr * g2_0_xxz[j];

                g_x_xyy[j] = g1_0_xxyy[j] - fr * g2_0_xyy[j];

                g_x_xyz[j] = g1_0_xxyz[j] - fr * g2_0_xyz[j];

                g_x_xzz[j] = g1_0_xxzz[j] - fr * g2_0_xzz[j];

                g_x_yyy[j] = g1_0_xyyy[j] - fr * g2_0_yyy[j];

                g_x_yyz[j] = g1_0_xyyz[j] - fr * g2_0_yyz[j];

                g_x_yzz[j] = g1_0_xyzz[j] - fr * g2_0_yzz[j];

                g_x_zzz[j] = g1_0_xzzz[j] - fr * g2_0_zzz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xxx[j] = g1_0_xxxy[j] - fr * g2_0_xxx[j];

                g_y_xxy[j] = g1_0_xxyy[j] - fr * g2_0_xxy[j];

                g_y_xxz[j] = g1_0_xxyz[j] - fr * g2_0_xxz[j];

                g_y_xyy[j] = g1_0_xyyy[j] - fr * g2_0_xyy[j];

                g_y_xyz[j] = g1_0_xyyz[j] - fr * g2_0_xyz[j];

                g_y_xzz[j] = g1_0_xyzz[j] - fr * g2_0_xzz[j];

                g_y_yyy[j] = g1_0_yyyy[j] - fr * g2_0_yyy[j];

                g_y_yyz[j] = g1_0_yyyz[j] - fr * g2_0_yyz[j];

                g_y_yzz[j] = g1_0_yyzz[j] - fr * g2_0_yzz[j];

                g_y_zzz[j] = g1_0_yzzz[j] - fr * g2_0_zzz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xxx[j] = g1_0_xxxz[j] - fr * g2_0_xxx[j];

                g_z_xxy[j] = g1_0_xxyz[j] - fr * g2_0_xxy[j];

                g_z_xxz[j] = g1_0_xxzz[j] - fr * g2_0_xxz[j];

                g_z_xyy[j] = g1_0_xyyz[j] - fr * g2_0_xyy[j];

                g_z_xyz[j] = g1_0_xyzz[j] - fr * g2_0_xyz[j];

                g_z_xzz[j] = g1_0_xzzz[j] - fr * g2_0_xzz[j];

                g_z_yyy[j] = g1_0_yyyz[j] - fr * g2_0_yyy[j];

                g_z_yyz[j] = g1_0_yyzz[j] - fr * g2_0_yyz[j];

                g_z_yzz[j] = g1_0_yzzz[j] - fr * g2_0_yzz[j];

                g_z_zzz[j] = g1_0_zzzz[j] - fr * g2_0_zzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPG(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 4})) return;

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
                                             {braAngularMomentum, 1, 4});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 5});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 4});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SG)^(m) integrals

            auto g2_0_xxxx = contrBuffer.data(g2off + 15 * i);

            auto g2_0_xxxy = contrBuffer.data(g2off + 15 * i + 1);

            auto g2_0_xxxz = contrBuffer.data(g2off + 15 * i + 2);

            auto g2_0_xxyy = contrBuffer.data(g2off + 15 * i + 3);

            auto g2_0_xxyz = contrBuffer.data(g2off + 15 * i + 4);

            auto g2_0_xxzz = contrBuffer.data(g2off + 15 * i + 5);

            auto g2_0_xyyy = contrBuffer.data(g2off + 15 * i + 6);

            auto g2_0_xyyz = contrBuffer.data(g2off + 15 * i + 7);

            auto g2_0_xyzz = contrBuffer.data(g2off + 15 * i + 8);

            auto g2_0_xzzz = contrBuffer.data(g2off + 15 * i + 9);

            auto g2_0_yyyy = contrBuffer.data(g2off + 15 * i + 10);

            auto g2_0_yyyz = contrBuffer.data(g2off + 15 * i + 11);

            auto g2_0_yyzz = contrBuffer.data(g2off + 15 * i + 12);

            auto g2_0_yzzz = contrBuffer.data(g2off + 15 * i + 13);

            auto g2_0_zzzz = contrBuffer.data(g2off + 15 * i + 14);

            // set up pointers to (X|g(r,r')|SH)^(m) integrals

            auto g1_0_xxxxx = contrBuffer.data(g1off + 21 * i);

            auto g1_0_xxxxy = contrBuffer.data(g1off + 21 * i + 1);

            auto g1_0_xxxxz = contrBuffer.data(g1off + 21 * i + 2);

            auto g1_0_xxxyy = contrBuffer.data(g1off + 21 * i + 3);

            auto g1_0_xxxyz = contrBuffer.data(g1off + 21 * i + 4);

            auto g1_0_xxxzz = contrBuffer.data(g1off + 21 * i + 5);

            auto g1_0_xxyyy = contrBuffer.data(g1off + 21 * i + 6);

            auto g1_0_xxyyz = contrBuffer.data(g1off + 21 * i + 7);

            auto g1_0_xxyzz = contrBuffer.data(g1off + 21 * i + 8);

            auto g1_0_xxzzz = contrBuffer.data(g1off + 21 * i + 9);

            auto g1_0_xyyyy = contrBuffer.data(g1off + 21 * i + 10);

            auto g1_0_xyyyz = contrBuffer.data(g1off + 21 * i + 11);

            auto g1_0_xyyzz = contrBuffer.data(g1off + 21 * i + 12);

            auto g1_0_xyzzz = contrBuffer.data(g1off + 21 * i + 13);

            auto g1_0_xzzzz = contrBuffer.data(g1off + 21 * i + 14);

            auto g1_0_yyyyy = contrBuffer.data(g1off + 21 * i + 15);

            auto g1_0_yyyyz = contrBuffer.data(g1off + 21 * i + 16);

            auto g1_0_yyyzz = contrBuffer.data(g1off + 21 * i + 17);

            auto g1_0_yyzzz = contrBuffer.data(g1off + 21 * i + 18);

            auto g1_0_yzzzz = contrBuffer.data(g1off + 21 * i + 19);

            auto g1_0_zzzzz = contrBuffer.data(g1off + 21 * i + 20);

            // set up pointers to (X|g(r,r')|PG)^(m) integrals

            auto g_x_xxxx = contrBuffer.data(goff + 45 * i);

            auto g_x_xxxy = contrBuffer.data(goff + 45 * i + 1);

            auto g_x_xxxz = contrBuffer.data(goff + 45 * i + 2);

            auto g_x_xxyy = contrBuffer.data(goff + 45 * i + 3);

            auto g_x_xxyz = contrBuffer.data(goff + 45 * i + 4);

            auto g_x_xxzz = contrBuffer.data(goff + 45 * i + 5);

            auto g_x_xyyy = contrBuffer.data(goff + 45 * i + 6);

            auto g_x_xyyz = contrBuffer.data(goff + 45 * i + 7);

            auto g_x_xyzz = contrBuffer.data(goff + 45 * i + 8);

            auto g_x_xzzz = contrBuffer.data(goff + 45 * i + 9);

            auto g_x_yyyy = contrBuffer.data(goff + 45 * i + 10);

            auto g_x_yyyz = contrBuffer.data(goff + 45 * i + 11);

            auto g_x_yyzz = contrBuffer.data(goff + 45 * i + 12);

            auto g_x_yzzz = contrBuffer.data(goff + 45 * i + 13);

            auto g_x_zzzz = contrBuffer.data(goff + 45 * i + 14);

            auto g_y_xxxx = contrBuffer.data(goff + 45 * i + 15);

            auto g_y_xxxy = contrBuffer.data(goff + 45 * i + 16);

            auto g_y_xxxz = contrBuffer.data(goff + 45 * i + 17);

            auto g_y_xxyy = contrBuffer.data(goff + 45 * i + 18);

            auto g_y_xxyz = contrBuffer.data(goff + 45 * i + 19);

            auto g_y_xxzz = contrBuffer.data(goff + 45 * i + 20);

            auto g_y_xyyy = contrBuffer.data(goff + 45 * i + 21);

            auto g_y_xyyz = contrBuffer.data(goff + 45 * i + 22);

            auto g_y_xyzz = contrBuffer.data(goff + 45 * i + 23);

            auto g_y_xzzz = contrBuffer.data(goff + 45 * i + 24);

            auto g_y_yyyy = contrBuffer.data(goff + 45 * i + 25);

            auto g_y_yyyz = contrBuffer.data(goff + 45 * i + 26);

            auto g_y_yyzz = contrBuffer.data(goff + 45 * i + 27);

            auto g_y_yzzz = contrBuffer.data(goff + 45 * i + 28);

            auto g_y_zzzz = contrBuffer.data(goff + 45 * i + 29);

            auto g_z_xxxx = contrBuffer.data(goff + 45 * i + 30);

            auto g_z_xxxy = contrBuffer.data(goff + 45 * i + 31);

            auto g_z_xxxz = contrBuffer.data(goff + 45 * i + 32);

            auto g_z_xxyy = contrBuffer.data(goff + 45 * i + 33);

            auto g_z_xxyz = contrBuffer.data(goff + 45 * i + 34);

            auto g_z_xxzz = contrBuffer.data(goff + 45 * i + 35);

            auto g_z_xyyy = contrBuffer.data(goff + 45 * i + 36);

            auto g_z_xyyz = contrBuffer.data(goff + 45 * i + 37);

            auto g_z_xyzz = contrBuffer.data(goff + 45 * i + 38);

            auto g_z_xzzz = contrBuffer.data(goff + 45 * i + 39);

            auto g_z_yyyy = contrBuffer.data(goff + 45 * i + 40);

            auto g_z_yyyz = contrBuffer.data(goff + 45 * i + 41);

            auto g_z_yyzz = contrBuffer.data(goff + 45 * i + 42);

            auto g_z_yzzz = contrBuffer.data(goff + 45 * i + 43);

            auto g_z_zzzz = contrBuffer.data(goff + 45 * i + 44);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxx, g2_0_xxxy, g2_0_xxxz,\
                                     g2_0_xxyy, g2_0_xxyz, g2_0_xxzz, g2_0_xyyy,\
                                     g2_0_xyyz, g2_0_xyzz, g2_0_xzzz, g2_0_yyyy,\
                                     g2_0_yyyz, g2_0_yyzz, g2_0_yzzz, g2_0_zzzz,\
                                     g1_0_xxxxx, g1_0_xxxxy, g1_0_xxxxz, g1_0_xxxyy,\
                                     g1_0_xxxyz, g1_0_xxxzz, g1_0_xxyyy, g1_0_xxyyz,\
                                     g1_0_xxyzz, g1_0_xxzzz, g1_0_xyyyy, g1_0_xyyyz,\
                                     g1_0_xyyzz, g1_0_xyzzz, g1_0_xzzzz, g1_0_yyyyy,\
                                     g1_0_yyyyz, g1_0_yyyzz, g1_0_yyzzz, g1_0_yzzzz,\
                                     g1_0_zzzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz,\
                                     g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy, g_x_xyyz,\
                                     g_x_xyzz, g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz,\
                                     g_x_yzzz, g_x_zzzz, g_y_xxxx, g_y_xxxy, g_y_xxxz,\
                                     g_y_xxyy, g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz,\
                                     g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz,\
                                     g_y_yzzz, g_y_zzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz,\
                                     g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy, g_z_xyyz,\
                                     g_z_xyzz, g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz,\
                                     g_z_yzzz, g_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xxxx[j] = g1_0_xxxxx[j] - fr * g2_0_xxxx[j];

                g_x_xxxy[j] = g1_0_xxxxy[j] - fr * g2_0_xxxy[j];

                g_x_xxxz[j] = g1_0_xxxxz[j] - fr * g2_0_xxxz[j];

                g_x_xxyy[j] = g1_0_xxxyy[j] - fr * g2_0_xxyy[j];

                g_x_xxyz[j] = g1_0_xxxyz[j] - fr * g2_0_xxyz[j];

                g_x_xxzz[j] = g1_0_xxxzz[j] - fr * g2_0_xxzz[j];

                g_x_xyyy[j] = g1_0_xxyyy[j] - fr * g2_0_xyyy[j];

                g_x_xyyz[j] = g1_0_xxyyz[j] - fr * g2_0_xyyz[j];

                g_x_xyzz[j] = g1_0_xxyzz[j] - fr * g2_0_xyzz[j];

                g_x_xzzz[j] = g1_0_xxzzz[j] - fr * g2_0_xzzz[j];

                g_x_yyyy[j] = g1_0_xyyyy[j] - fr * g2_0_yyyy[j];

                g_x_yyyz[j] = g1_0_xyyyz[j] - fr * g2_0_yyyz[j];

                g_x_yyzz[j] = g1_0_xyyzz[j] - fr * g2_0_yyzz[j];

                g_x_yzzz[j] = g1_0_xyzzz[j] - fr * g2_0_yzzz[j];

                g_x_zzzz[j] = g1_0_xzzzz[j] - fr * g2_0_zzzz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xxxx[j] = g1_0_xxxxy[j] - fr * g2_0_xxxx[j];

                g_y_xxxy[j] = g1_0_xxxyy[j] - fr * g2_0_xxxy[j];

                g_y_xxxz[j] = g1_0_xxxyz[j] - fr * g2_0_xxxz[j];

                g_y_xxyy[j] = g1_0_xxyyy[j] - fr * g2_0_xxyy[j];

                g_y_xxyz[j] = g1_0_xxyyz[j] - fr * g2_0_xxyz[j];

                g_y_xxzz[j] = g1_0_xxyzz[j] - fr * g2_0_xxzz[j];

                g_y_xyyy[j] = g1_0_xyyyy[j] - fr * g2_0_xyyy[j];

                g_y_xyyz[j] = g1_0_xyyyz[j] - fr * g2_0_xyyz[j];

                g_y_xyzz[j] = g1_0_xyyzz[j] - fr * g2_0_xyzz[j];

                g_y_xzzz[j] = g1_0_xyzzz[j] - fr * g2_0_xzzz[j];

                g_y_yyyy[j] = g1_0_yyyyy[j] - fr * g2_0_yyyy[j];

                g_y_yyyz[j] = g1_0_yyyyz[j] - fr * g2_0_yyyz[j];

                g_y_yyzz[j] = g1_0_yyyzz[j] - fr * g2_0_yyzz[j];

                g_y_yzzz[j] = g1_0_yyzzz[j] - fr * g2_0_yzzz[j];

                g_y_zzzz[j] = g1_0_yzzzz[j] - fr * g2_0_zzzz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xxxx[j] = g1_0_xxxxz[j] - fr * g2_0_xxxx[j];

                g_z_xxxy[j] = g1_0_xxxyz[j] - fr * g2_0_xxxy[j];

                g_z_xxxz[j] = g1_0_xxxzz[j] - fr * g2_0_xxxz[j];

                g_z_xxyy[j] = g1_0_xxyyz[j] - fr * g2_0_xxyy[j];

                g_z_xxyz[j] = g1_0_xxyzz[j] - fr * g2_0_xxyz[j];

                g_z_xxzz[j] = g1_0_xxzzz[j] - fr * g2_0_xxzz[j];

                g_z_xyyy[j] = g1_0_xyyyz[j] - fr * g2_0_xyyy[j];

                g_z_xyyz[j] = g1_0_xyyzz[j] - fr * g2_0_xyyz[j];

                g_z_xyzz[j] = g1_0_xyzzz[j] - fr * g2_0_xyzz[j];

                g_z_xzzz[j] = g1_0_xzzzz[j] - fr * g2_0_xzzz[j];

                g_z_yyyy[j] = g1_0_yyyyz[j] - fr * g2_0_yyyy[j];

                g_z_yyyz[j] = g1_0_yyyzz[j] - fr * g2_0_yyyz[j];

                g_z_yyzz[j] = g1_0_yyzzz[j] - fr * g2_0_yyzz[j];

                g_z_yzzz[j] = g1_0_yzzzz[j] - fr * g2_0_yzzz[j];

                g_z_zzzz[j] = g1_0_zzzzz[j] - fr * g2_0_zzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPH(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 5})) return;

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
                                             {braAngularMomentum, 1, 5});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 6});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 5});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SH)^(m) integrals

            auto g2_0_xxxxx = contrBuffer.data(g2off + 21 * i);

            auto g2_0_xxxxy = contrBuffer.data(g2off + 21 * i + 1);

            auto g2_0_xxxxz = contrBuffer.data(g2off + 21 * i + 2);

            auto g2_0_xxxyy = contrBuffer.data(g2off + 21 * i + 3);

            auto g2_0_xxxyz = contrBuffer.data(g2off + 21 * i + 4);

            auto g2_0_xxxzz = contrBuffer.data(g2off + 21 * i + 5);

            auto g2_0_xxyyy = contrBuffer.data(g2off + 21 * i + 6);

            auto g2_0_xxyyz = contrBuffer.data(g2off + 21 * i + 7);

            auto g2_0_xxyzz = contrBuffer.data(g2off + 21 * i + 8);

            auto g2_0_xxzzz = contrBuffer.data(g2off + 21 * i + 9);

            auto g2_0_xyyyy = contrBuffer.data(g2off + 21 * i + 10);

            auto g2_0_xyyyz = contrBuffer.data(g2off + 21 * i + 11);

            auto g2_0_xyyzz = contrBuffer.data(g2off + 21 * i + 12);

            auto g2_0_xyzzz = contrBuffer.data(g2off + 21 * i + 13);

            auto g2_0_xzzzz = contrBuffer.data(g2off + 21 * i + 14);

            auto g2_0_yyyyy = contrBuffer.data(g2off + 21 * i + 15);

            auto g2_0_yyyyz = contrBuffer.data(g2off + 21 * i + 16);

            auto g2_0_yyyzz = contrBuffer.data(g2off + 21 * i + 17);

            auto g2_0_yyzzz = contrBuffer.data(g2off + 21 * i + 18);

            auto g2_0_yzzzz = contrBuffer.data(g2off + 21 * i + 19);

            auto g2_0_zzzzz = contrBuffer.data(g2off + 21 * i + 20);

            // set up pointers to (X|g(r,r')|SI)^(m) integrals

            auto g1_0_xxxxxx = contrBuffer.data(g1off + 28 * i);

            auto g1_0_xxxxxy = contrBuffer.data(g1off + 28 * i + 1);

            auto g1_0_xxxxxz = contrBuffer.data(g1off + 28 * i + 2);

            auto g1_0_xxxxyy = contrBuffer.data(g1off + 28 * i + 3);

            auto g1_0_xxxxyz = contrBuffer.data(g1off + 28 * i + 4);

            auto g1_0_xxxxzz = contrBuffer.data(g1off + 28 * i + 5);

            auto g1_0_xxxyyy = contrBuffer.data(g1off + 28 * i + 6);

            auto g1_0_xxxyyz = contrBuffer.data(g1off + 28 * i + 7);

            auto g1_0_xxxyzz = contrBuffer.data(g1off + 28 * i + 8);

            auto g1_0_xxxzzz = contrBuffer.data(g1off + 28 * i + 9);

            auto g1_0_xxyyyy = contrBuffer.data(g1off + 28 * i + 10);

            auto g1_0_xxyyyz = contrBuffer.data(g1off + 28 * i + 11);

            auto g1_0_xxyyzz = contrBuffer.data(g1off + 28 * i + 12);

            auto g1_0_xxyzzz = contrBuffer.data(g1off + 28 * i + 13);

            auto g1_0_xxzzzz = contrBuffer.data(g1off + 28 * i + 14);

            auto g1_0_xyyyyy = contrBuffer.data(g1off + 28 * i + 15);

            auto g1_0_xyyyyz = contrBuffer.data(g1off + 28 * i + 16);

            auto g1_0_xyyyzz = contrBuffer.data(g1off + 28 * i + 17);

            auto g1_0_xyyzzz = contrBuffer.data(g1off + 28 * i + 18);

            auto g1_0_xyzzzz = contrBuffer.data(g1off + 28 * i + 19);

            auto g1_0_xzzzzz = contrBuffer.data(g1off + 28 * i + 20);

            auto g1_0_yyyyyy = contrBuffer.data(g1off + 28 * i + 21);

            auto g1_0_yyyyyz = contrBuffer.data(g1off + 28 * i + 22);

            auto g1_0_yyyyzz = contrBuffer.data(g1off + 28 * i + 23);

            auto g1_0_yyyzzz = contrBuffer.data(g1off + 28 * i + 24);

            auto g1_0_yyzzzz = contrBuffer.data(g1off + 28 * i + 25);

            auto g1_0_yzzzzz = contrBuffer.data(g1off + 28 * i + 26);

            auto g1_0_zzzzzz = contrBuffer.data(g1off + 28 * i + 27);

            // set up pointers to (X|g(r,r')|PH)^(m) integrals

            auto g_x_xxxxx = contrBuffer.data(goff + 63 * i);

            auto g_x_xxxxy = contrBuffer.data(goff + 63 * i + 1);

            auto g_x_xxxxz = contrBuffer.data(goff + 63 * i + 2);

            auto g_x_xxxyy = contrBuffer.data(goff + 63 * i + 3);

            auto g_x_xxxyz = contrBuffer.data(goff + 63 * i + 4);

            auto g_x_xxxzz = contrBuffer.data(goff + 63 * i + 5);

            auto g_x_xxyyy = contrBuffer.data(goff + 63 * i + 6);

            auto g_x_xxyyz = contrBuffer.data(goff + 63 * i + 7);

            auto g_x_xxyzz = contrBuffer.data(goff + 63 * i + 8);

            auto g_x_xxzzz = contrBuffer.data(goff + 63 * i + 9);

            auto g_x_xyyyy = contrBuffer.data(goff + 63 * i + 10);

            auto g_x_xyyyz = contrBuffer.data(goff + 63 * i + 11);

            auto g_x_xyyzz = contrBuffer.data(goff + 63 * i + 12);

            auto g_x_xyzzz = contrBuffer.data(goff + 63 * i + 13);

            auto g_x_xzzzz = contrBuffer.data(goff + 63 * i + 14);

            auto g_x_yyyyy = contrBuffer.data(goff + 63 * i + 15);

            auto g_x_yyyyz = contrBuffer.data(goff + 63 * i + 16);

            auto g_x_yyyzz = contrBuffer.data(goff + 63 * i + 17);

            auto g_x_yyzzz = contrBuffer.data(goff + 63 * i + 18);

            auto g_x_yzzzz = contrBuffer.data(goff + 63 * i + 19);

            auto g_x_zzzzz = contrBuffer.data(goff + 63 * i + 20);

            auto g_y_xxxxx = contrBuffer.data(goff + 63 * i + 21);

            auto g_y_xxxxy = contrBuffer.data(goff + 63 * i + 22);

            auto g_y_xxxxz = contrBuffer.data(goff + 63 * i + 23);

            auto g_y_xxxyy = contrBuffer.data(goff + 63 * i + 24);

            auto g_y_xxxyz = contrBuffer.data(goff + 63 * i + 25);

            auto g_y_xxxzz = contrBuffer.data(goff + 63 * i + 26);

            auto g_y_xxyyy = contrBuffer.data(goff + 63 * i + 27);

            auto g_y_xxyyz = contrBuffer.data(goff + 63 * i + 28);

            auto g_y_xxyzz = contrBuffer.data(goff + 63 * i + 29);

            auto g_y_xxzzz = contrBuffer.data(goff + 63 * i + 30);

            auto g_y_xyyyy = contrBuffer.data(goff + 63 * i + 31);

            auto g_y_xyyyz = contrBuffer.data(goff + 63 * i + 32);

            auto g_y_xyyzz = contrBuffer.data(goff + 63 * i + 33);

            auto g_y_xyzzz = contrBuffer.data(goff + 63 * i + 34);

            auto g_y_xzzzz = contrBuffer.data(goff + 63 * i + 35);

            auto g_y_yyyyy = contrBuffer.data(goff + 63 * i + 36);

            auto g_y_yyyyz = contrBuffer.data(goff + 63 * i + 37);

            auto g_y_yyyzz = contrBuffer.data(goff + 63 * i + 38);

            auto g_y_yyzzz = contrBuffer.data(goff + 63 * i + 39);

            auto g_y_yzzzz = contrBuffer.data(goff + 63 * i + 40);

            auto g_y_zzzzz = contrBuffer.data(goff + 63 * i + 41);

            auto g_z_xxxxx = contrBuffer.data(goff + 63 * i + 42);

            auto g_z_xxxxy = contrBuffer.data(goff + 63 * i + 43);

            auto g_z_xxxxz = contrBuffer.data(goff + 63 * i + 44);

            auto g_z_xxxyy = contrBuffer.data(goff + 63 * i + 45);

            auto g_z_xxxyz = contrBuffer.data(goff + 63 * i + 46);

            auto g_z_xxxzz = contrBuffer.data(goff + 63 * i + 47);

            auto g_z_xxyyy = contrBuffer.data(goff + 63 * i + 48);

            auto g_z_xxyyz = contrBuffer.data(goff + 63 * i + 49);

            auto g_z_xxyzz = contrBuffer.data(goff + 63 * i + 50);

            auto g_z_xxzzz = contrBuffer.data(goff + 63 * i + 51);

            auto g_z_xyyyy = contrBuffer.data(goff + 63 * i + 52);

            auto g_z_xyyyz = contrBuffer.data(goff + 63 * i + 53);

            auto g_z_xyyzz = contrBuffer.data(goff + 63 * i + 54);

            auto g_z_xyzzz = contrBuffer.data(goff + 63 * i + 55);

            auto g_z_xzzzz = contrBuffer.data(goff + 63 * i + 56);

            auto g_z_yyyyy = contrBuffer.data(goff + 63 * i + 57);

            auto g_z_yyyyz = contrBuffer.data(goff + 63 * i + 58);

            auto g_z_yyyzz = contrBuffer.data(goff + 63 * i + 59);

            auto g_z_yyzzz = contrBuffer.data(goff + 63 * i + 60);

            auto g_z_yzzzz = contrBuffer.data(goff + 63 * i + 61);

            auto g_z_zzzzz = contrBuffer.data(goff + 63 * i + 62);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxxx, g2_0_xxxxy,\
                                     g2_0_xxxxz, g2_0_xxxyy, g2_0_xxxyz, g2_0_xxxzz,\
                                     g2_0_xxyyy, g2_0_xxyyz, g2_0_xxyzz, g2_0_xxzzz,\
                                     g2_0_xyyyy, g2_0_xyyyz, g2_0_xyyzz, g2_0_xyzzz,\
                                     g2_0_xzzzz, g2_0_yyyyy, g2_0_yyyyz, g2_0_yyyzz,\
                                     g2_0_yyzzz, g2_0_yzzzz, g2_0_zzzzz, g1_0_xxxxxx,\
                                     g1_0_xxxxxy, g1_0_xxxxxz, g1_0_xxxxyy, g1_0_xxxxyz,\
                                     g1_0_xxxxzz, g1_0_xxxyyy, g1_0_xxxyyz, g1_0_xxxyzz,\
                                     g1_0_xxxzzz, g1_0_xxyyyy, g1_0_xxyyyz, g1_0_xxyyzz,\
                                     g1_0_xxyzzz, g1_0_xxzzzz, g1_0_xyyyyy, g1_0_xyyyyz,\
                                     g1_0_xyyyzz, g1_0_xyyzzz, g1_0_xyzzzz, g1_0_xzzzzz,\
                                     g1_0_yyyyyy, g1_0_yyyyyz, g1_0_yyyyzz, g1_0_yyyzzz,\
                                     g1_0_yyzzzz, g1_0_yzzzzz, g1_0_zzzzzz, g_x_xxxxx,\
                                     g_x_xxxxy, g_x_xxxxz, g_x_xxxyy, g_x_xxxyz,\
                                     g_x_xxxzz, g_x_xxyyy, g_x_xxyyz, g_x_xxyzz,\
                                     g_x_xxzzz, g_x_xyyyy, g_x_xyyyz, g_x_xyyzz,\
                                     g_x_xyzzz, g_x_xzzzz, g_x_yyyyy, g_x_yyyyz,\
                                     g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz,\
                                     g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxyy,\
                                     g_y_xxxyz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyz,\
                                     g_y_xxyzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyz,\
                                     g_y_xyyzz, g_y_xyzzz, g_y_xzzzz, g_y_yyyyy,\
                                     g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz,\
                                     g_y_zzzzz, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz,\
                                     g_z_xxxyy, g_z_xxxyz, g_z_xxxzz, g_z_xxyyy,\
                                     g_z_xxyyz, g_z_xxyzz, g_z_xxzzz, g_z_xyyyy,\
                                     g_z_xyyyz, g_z_xyyzz, g_z_xyzzz, g_z_xzzzz,\
                                     g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz,\
                                     g_z_yzzzz, g_z_zzzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xxxxx[j] = g1_0_xxxxxx[j] - fr * g2_0_xxxxx[j];

                g_x_xxxxy[j] = g1_0_xxxxxy[j] - fr * g2_0_xxxxy[j];

                g_x_xxxxz[j] = g1_0_xxxxxz[j] - fr * g2_0_xxxxz[j];

                g_x_xxxyy[j] = g1_0_xxxxyy[j] - fr * g2_0_xxxyy[j];

                g_x_xxxyz[j] = g1_0_xxxxyz[j] - fr * g2_0_xxxyz[j];

                g_x_xxxzz[j] = g1_0_xxxxzz[j] - fr * g2_0_xxxzz[j];

                g_x_xxyyy[j] = g1_0_xxxyyy[j] - fr * g2_0_xxyyy[j];

                g_x_xxyyz[j] = g1_0_xxxyyz[j] - fr * g2_0_xxyyz[j];

                g_x_xxyzz[j] = g1_0_xxxyzz[j] - fr * g2_0_xxyzz[j];

                g_x_xxzzz[j] = g1_0_xxxzzz[j] - fr * g2_0_xxzzz[j];

                g_x_xyyyy[j] = g1_0_xxyyyy[j] - fr * g2_0_xyyyy[j];

                g_x_xyyyz[j] = g1_0_xxyyyz[j] - fr * g2_0_xyyyz[j];

                g_x_xyyzz[j] = g1_0_xxyyzz[j] - fr * g2_0_xyyzz[j];

                g_x_xyzzz[j] = g1_0_xxyzzz[j] - fr * g2_0_xyzzz[j];

                g_x_xzzzz[j] = g1_0_xxzzzz[j] - fr * g2_0_xzzzz[j];

                g_x_yyyyy[j] = g1_0_xyyyyy[j] - fr * g2_0_yyyyy[j];

                g_x_yyyyz[j] = g1_0_xyyyyz[j] - fr * g2_0_yyyyz[j];

                g_x_yyyzz[j] = g1_0_xyyyzz[j] - fr * g2_0_yyyzz[j];

                g_x_yyzzz[j] = g1_0_xyyzzz[j] - fr * g2_0_yyzzz[j];

                g_x_yzzzz[j] = g1_0_xyzzzz[j] - fr * g2_0_yzzzz[j];

                g_x_zzzzz[j] = g1_0_xzzzzz[j] - fr * g2_0_zzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xxxxx[j] = g1_0_xxxxxy[j] - fr * g2_0_xxxxx[j];

                g_y_xxxxy[j] = g1_0_xxxxyy[j] - fr * g2_0_xxxxy[j];

                g_y_xxxxz[j] = g1_0_xxxxyz[j] - fr * g2_0_xxxxz[j];

                g_y_xxxyy[j] = g1_0_xxxyyy[j] - fr * g2_0_xxxyy[j];

                g_y_xxxyz[j] = g1_0_xxxyyz[j] - fr * g2_0_xxxyz[j];

                g_y_xxxzz[j] = g1_0_xxxyzz[j] - fr * g2_0_xxxzz[j];

                g_y_xxyyy[j] = g1_0_xxyyyy[j] - fr * g2_0_xxyyy[j];

                g_y_xxyyz[j] = g1_0_xxyyyz[j] - fr * g2_0_xxyyz[j];

                g_y_xxyzz[j] = g1_0_xxyyzz[j] - fr * g2_0_xxyzz[j];

                g_y_xxzzz[j] = g1_0_xxyzzz[j] - fr * g2_0_xxzzz[j];

                g_y_xyyyy[j] = g1_0_xyyyyy[j] - fr * g2_0_xyyyy[j];

                g_y_xyyyz[j] = g1_0_xyyyyz[j] - fr * g2_0_xyyyz[j];

                g_y_xyyzz[j] = g1_0_xyyyzz[j] - fr * g2_0_xyyzz[j];

                g_y_xyzzz[j] = g1_0_xyyzzz[j] - fr * g2_0_xyzzz[j];

                g_y_xzzzz[j] = g1_0_xyzzzz[j] - fr * g2_0_xzzzz[j];

                g_y_yyyyy[j] = g1_0_yyyyyy[j] - fr * g2_0_yyyyy[j];

                g_y_yyyyz[j] = g1_0_yyyyyz[j] - fr * g2_0_yyyyz[j];

                g_y_yyyzz[j] = g1_0_yyyyzz[j] - fr * g2_0_yyyzz[j];

                g_y_yyzzz[j] = g1_0_yyyzzz[j] - fr * g2_0_yyzzz[j];

                g_y_yzzzz[j] = g1_0_yyzzzz[j] - fr * g2_0_yzzzz[j];

                g_y_zzzzz[j] = g1_0_yzzzzz[j] - fr * g2_0_zzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xxxxx[j] = g1_0_xxxxxz[j] - fr * g2_0_xxxxx[j];

                g_z_xxxxy[j] = g1_0_xxxxyz[j] - fr * g2_0_xxxxy[j];

                g_z_xxxxz[j] = g1_0_xxxxzz[j] - fr * g2_0_xxxxz[j];

                g_z_xxxyy[j] = g1_0_xxxyyz[j] - fr * g2_0_xxxyy[j];

                g_z_xxxyz[j] = g1_0_xxxyzz[j] - fr * g2_0_xxxyz[j];

                g_z_xxxzz[j] = g1_0_xxxzzz[j] - fr * g2_0_xxxzz[j];

                g_z_xxyyy[j] = g1_0_xxyyyz[j] - fr * g2_0_xxyyy[j];

                g_z_xxyyz[j] = g1_0_xxyyzz[j] - fr * g2_0_xxyyz[j];

                g_z_xxyzz[j] = g1_0_xxyzzz[j] - fr * g2_0_xxyzz[j];

                g_z_xxzzz[j] = g1_0_xxzzzz[j] - fr * g2_0_xxzzz[j];

                g_z_xyyyy[j] = g1_0_xyyyyz[j] - fr * g2_0_xyyyy[j];

                g_z_xyyyz[j] = g1_0_xyyyzz[j] - fr * g2_0_xyyyz[j];

                g_z_xyyzz[j] = g1_0_xyyzzz[j] - fr * g2_0_xyyzz[j];

                g_z_xyzzz[j] = g1_0_xyzzzz[j] - fr * g2_0_xyzzz[j];

                g_z_xzzzz[j] = g1_0_xzzzzz[j] - fr * g2_0_xzzzz[j];

                g_z_yyyyy[j] = g1_0_yyyyyz[j] - fr * g2_0_yyyyy[j];

                g_z_yyyyz[j] = g1_0_yyyyzz[j] - fr * g2_0_yyyyz[j];

                g_z_yyyzz[j] = g1_0_yyyzzz[j] - fr * g2_0_yyyzz[j];

                g_z_yyzzz[j] = g1_0_yyzzzz[j] - fr * g2_0_yyzzz[j];

                g_z_yzzzz[j] = g1_0_yzzzzz[j] - fr * g2_0_yzzzz[j];

                g_z_zzzzz[j] = g1_0_zzzzzz[j] - fr * g2_0_zzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPI(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 6})) return;

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
                                             {braAngularMomentum, 1, 6});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 7});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 6});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SI)^(m) integrals

            auto g2_0_xxxxxx = contrBuffer.data(g2off + 28 * i);

            auto g2_0_xxxxxy = contrBuffer.data(g2off + 28 * i + 1);

            auto g2_0_xxxxxz = contrBuffer.data(g2off + 28 * i + 2);

            auto g2_0_xxxxyy = contrBuffer.data(g2off + 28 * i + 3);

            auto g2_0_xxxxyz = contrBuffer.data(g2off + 28 * i + 4);

            auto g2_0_xxxxzz = contrBuffer.data(g2off + 28 * i + 5);

            auto g2_0_xxxyyy = contrBuffer.data(g2off + 28 * i + 6);

            auto g2_0_xxxyyz = contrBuffer.data(g2off + 28 * i + 7);

            auto g2_0_xxxyzz = contrBuffer.data(g2off + 28 * i + 8);

            auto g2_0_xxxzzz = contrBuffer.data(g2off + 28 * i + 9);

            auto g2_0_xxyyyy = contrBuffer.data(g2off + 28 * i + 10);

            auto g2_0_xxyyyz = contrBuffer.data(g2off + 28 * i + 11);

            auto g2_0_xxyyzz = contrBuffer.data(g2off + 28 * i + 12);

            auto g2_0_xxyzzz = contrBuffer.data(g2off + 28 * i + 13);

            auto g2_0_xxzzzz = contrBuffer.data(g2off + 28 * i + 14);

            auto g2_0_xyyyyy = contrBuffer.data(g2off + 28 * i + 15);

            auto g2_0_xyyyyz = contrBuffer.data(g2off + 28 * i + 16);

            auto g2_0_xyyyzz = contrBuffer.data(g2off + 28 * i + 17);

            auto g2_0_xyyzzz = contrBuffer.data(g2off + 28 * i + 18);

            auto g2_0_xyzzzz = contrBuffer.data(g2off + 28 * i + 19);

            auto g2_0_xzzzzz = contrBuffer.data(g2off + 28 * i + 20);

            auto g2_0_yyyyyy = contrBuffer.data(g2off + 28 * i + 21);

            auto g2_0_yyyyyz = contrBuffer.data(g2off + 28 * i + 22);

            auto g2_0_yyyyzz = contrBuffer.data(g2off + 28 * i + 23);

            auto g2_0_yyyzzz = contrBuffer.data(g2off + 28 * i + 24);

            auto g2_0_yyzzzz = contrBuffer.data(g2off + 28 * i + 25);

            auto g2_0_yzzzzz = contrBuffer.data(g2off + 28 * i + 26);

            auto g2_0_zzzzzz = contrBuffer.data(g2off + 28 * i + 27);

            // set up pointers to (X|g(r,r')|SK)^(m) integrals

            auto g1_0_xxxxxxx = contrBuffer.data(g1off + 36 * i);

            auto g1_0_xxxxxxy = contrBuffer.data(g1off + 36 * i + 1);

            auto g1_0_xxxxxxz = contrBuffer.data(g1off + 36 * i + 2);

            auto g1_0_xxxxxyy = contrBuffer.data(g1off + 36 * i + 3);

            auto g1_0_xxxxxyz = contrBuffer.data(g1off + 36 * i + 4);

            auto g1_0_xxxxxzz = contrBuffer.data(g1off + 36 * i + 5);

            auto g1_0_xxxxyyy = contrBuffer.data(g1off + 36 * i + 6);

            auto g1_0_xxxxyyz = contrBuffer.data(g1off + 36 * i + 7);

            auto g1_0_xxxxyzz = contrBuffer.data(g1off + 36 * i + 8);

            auto g1_0_xxxxzzz = contrBuffer.data(g1off + 36 * i + 9);

            auto g1_0_xxxyyyy = contrBuffer.data(g1off + 36 * i + 10);

            auto g1_0_xxxyyyz = contrBuffer.data(g1off + 36 * i + 11);

            auto g1_0_xxxyyzz = contrBuffer.data(g1off + 36 * i + 12);

            auto g1_0_xxxyzzz = contrBuffer.data(g1off + 36 * i + 13);

            auto g1_0_xxxzzzz = contrBuffer.data(g1off + 36 * i + 14);

            auto g1_0_xxyyyyy = contrBuffer.data(g1off + 36 * i + 15);

            auto g1_0_xxyyyyz = contrBuffer.data(g1off + 36 * i + 16);

            auto g1_0_xxyyyzz = contrBuffer.data(g1off + 36 * i + 17);

            auto g1_0_xxyyzzz = contrBuffer.data(g1off + 36 * i + 18);

            auto g1_0_xxyzzzz = contrBuffer.data(g1off + 36 * i + 19);

            auto g1_0_xxzzzzz = contrBuffer.data(g1off + 36 * i + 20);

            auto g1_0_xyyyyyy = contrBuffer.data(g1off + 36 * i + 21);

            auto g1_0_xyyyyyz = contrBuffer.data(g1off + 36 * i + 22);

            auto g1_0_xyyyyzz = contrBuffer.data(g1off + 36 * i + 23);

            auto g1_0_xyyyzzz = contrBuffer.data(g1off + 36 * i + 24);

            auto g1_0_xyyzzzz = contrBuffer.data(g1off + 36 * i + 25);

            auto g1_0_xyzzzzz = contrBuffer.data(g1off + 36 * i + 26);

            auto g1_0_xzzzzzz = contrBuffer.data(g1off + 36 * i + 27);

            auto g1_0_yyyyyyy = contrBuffer.data(g1off + 36 * i + 28);

            auto g1_0_yyyyyyz = contrBuffer.data(g1off + 36 * i + 29);

            auto g1_0_yyyyyzz = contrBuffer.data(g1off + 36 * i + 30);

            auto g1_0_yyyyzzz = contrBuffer.data(g1off + 36 * i + 31);

            auto g1_0_yyyzzzz = contrBuffer.data(g1off + 36 * i + 32);

            auto g1_0_yyzzzzz = contrBuffer.data(g1off + 36 * i + 33);

            auto g1_0_yzzzzzz = contrBuffer.data(g1off + 36 * i + 34);

            auto g1_0_zzzzzzz = contrBuffer.data(g1off + 36 * i + 35);

            // set up pointers to (X|g(r,r')|PI)^(m) integrals

            auto g_x_xxxxxx = contrBuffer.data(goff + 84 * i);

            auto g_x_xxxxxy = contrBuffer.data(goff + 84 * i + 1);

            auto g_x_xxxxxz = contrBuffer.data(goff + 84 * i + 2);

            auto g_x_xxxxyy = contrBuffer.data(goff + 84 * i + 3);

            auto g_x_xxxxyz = contrBuffer.data(goff + 84 * i + 4);

            auto g_x_xxxxzz = contrBuffer.data(goff + 84 * i + 5);

            auto g_x_xxxyyy = contrBuffer.data(goff + 84 * i + 6);

            auto g_x_xxxyyz = contrBuffer.data(goff + 84 * i + 7);

            auto g_x_xxxyzz = contrBuffer.data(goff + 84 * i + 8);

            auto g_x_xxxzzz = contrBuffer.data(goff + 84 * i + 9);

            auto g_x_xxyyyy = contrBuffer.data(goff + 84 * i + 10);

            auto g_x_xxyyyz = contrBuffer.data(goff + 84 * i + 11);

            auto g_x_xxyyzz = contrBuffer.data(goff + 84 * i + 12);

            auto g_x_xxyzzz = contrBuffer.data(goff + 84 * i + 13);

            auto g_x_xxzzzz = contrBuffer.data(goff + 84 * i + 14);

            auto g_x_xyyyyy = contrBuffer.data(goff + 84 * i + 15);

            auto g_x_xyyyyz = contrBuffer.data(goff + 84 * i + 16);

            auto g_x_xyyyzz = contrBuffer.data(goff + 84 * i + 17);

            auto g_x_xyyzzz = contrBuffer.data(goff + 84 * i + 18);

            auto g_x_xyzzzz = contrBuffer.data(goff + 84 * i + 19);

            auto g_x_xzzzzz = contrBuffer.data(goff + 84 * i + 20);

            auto g_x_yyyyyy = contrBuffer.data(goff + 84 * i + 21);

            auto g_x_yyyyyz = contrBuffer.data(goff + 84 * i + 22);

            auto g_x_yyyyzz = contrBuffer.data(goff + 84 * i + 23);

            auto g_x_yyyzzz = contrBuffer.data(goff + 84 * i + 24);

            auto g_x_yyzzzz = contrBuffer.data(goff + 84 * i + 25);

            auto g_x_yzzzzz = contrBuffer.data(goff + 84 * i + 26);

            auto g_x_zzzzzz = contrBuffer.data(goff + 84 * i + 27);

            auto g_y_xxxxxx = contrBuffer.data(goff + 84 * i + 28);

            auto g_y_xxxxxy = contrBuffer.data(goff + 84 * i + 29);

            auto g_y_xxxxxz = contrBuffer.data(goff + 84 * i + 30);

            auto g_y_xxxxyy = contrBuffer.data(goff + 84 * i + 31);

            auto g_y_xxxxyz = contrBuffer.data(goff + 84 * i + 32);

            auto g_y_xxxxzz = contrBuffer.data(goff + 84 * i + 33);

            auto g_y_xxxyyy = contrBuffer.data(goff + 84 * i + 34);

            auto g_y_xxxyyz = contrBuffer.data(goff + 84 * i + 35);

            auto g_y_xxxyzz = contrBuffer.data(goff + 84 * i + 36);

            auto g_y_xxxzzz = contrBuffer.data(goff + 84 * i + 37);

            auto g_y_xxyyyy = contrBuffer.data(goff + 84 * i + 38);

            auto g_y_xxyyyz = contrBuffer.data(goff + 84 * i + 39);

            auto g_y_xxyyzz = contrBuffer.data(goff + 84 * i + 40);

            auto g_y_xxyzzz = contrBuffer.data(goff + 84 * i + 41);

            auto g_y_xxzzzz = contrBuffer.data(goff + 84 * i + 42);

            auto g_y_xyyyyy = contrBuffer.data(goff + 84 * i + 43);

            auto g_y_xyyyyz = contrBuffer.data(goff + 84 * i + 44);

            auto g_y_xyyyzz = contrBuffer.data(goff + 84 * i + 45);

            auto g_y_xyyzzz = contrBuffer.data(goff + 84 * i + 46);

            auto g_y_xyzzzz = contrBuffer.data(goff + 84 * i + 47);

            auto g_y_xzzzzz = contrBuffer.data(goff + 84 * i + 48);

            auto g_y_yyyyyy = contrBuffer.data(goff + 84 * i + 49);

            auto g_y_yyyyyz = contrBuffer.data(goff + 84 * i + 50);

            auto g_y_yyyyzz = contrBuffer.data(goff + 84 * i + 51);

            auto g_y_yyyzzz = contrBuffer.data(goff + 84 * i + 52);

            auto g_y_yyzzzz = contrBuffer.data(goff + 84 * i + 53);

            auto g_y_yzzzzz = contrBuffer.data(goff + 84 * i + 54);

            auto g_y_zzzzzz = contrBuffer.data(goff + 84 * i + 55);

            auto g_z_xxxxxx = contrBuffer.data(goff + 84 * i + 56);

            auto g_z_xxxxxy = contrBuffer.data(goff + 84 * i + 57);

            auto g_z_xxxxxz = contrBuffer.data(goff + 84 * i + 58);

            auto g_z_xxxxyy = contrBuffer.data(goff + 84 * i + 59);

            auto g_z_xxxxyz = contrBuffer.data(goff + 84 * i + 60);

            auto g_z_xxxxzz = contrBuffer.data(goff + 84 * i + 61);

            auto g_z_xxxyyy = contrBuffer.data(goff + 84 * i + 62);

            auto g_z_xxxyyz = contrBuffer.data(goff + 84 * i + 63);

            auto g_z_xxxyzz = contrBuffer.data(goff + 84 * i + 64);

            auto g_z_xxxzzz = contrBuffer.data(goff + 84 * i + 65);

            auto g_z_xxyyyy = contrBuffer.data(goff + 84 * i + 66);

            auto g_z_xxyyyz = contrBuffer.data(goff + 84 * i + 67);

            auto g_z_xxyyzz = contrBuffer.data(goff + 84 * i + 68);

            auto g_z_xxyzzz = contrBuffer.data(goff + 84 * i + 69);

            auto g_z_xxzzzz = contrBuffer.data(goff + 84 * i + 70);

            auto g_z_xyyyyy = contrBuffer.data(goff + 84 * i + 71);

            auto g_z_xyyyyz = contrBuffer.data(goff + 84 * i + 72);

            auto g_z_xyyyzz = contrBuffer.data(goff + 84 * i + 73);

            auto g_z_xyyzzz = contrBuffer.data(goff + 84 * i + 74);

            auto g_z_xyzzzz = contrBuffer.data(goff + 84 * i + 75);

            auto g_z_xzzzzz = contrBuffer.data(goff + 84 * i + 76);

            auto g_z_yyyyyy = contrBuffer.data(goff + 84 * i + 77);

            auto g_z_yyyyyz = contrBuffer.data(goff + 84 * i + 78);

            auto g_z_yyyyzz = contrBuffer.data(goff + 84 * i + 79);

            auto g_z_yyyzzz = contrBuffer.data(goff + 84 * i + 80);

            auto g_z_yyzzzz = contrBuffer.data(goff + 84 * i + 81);

            auto g_z_yzzzzz = contrBuffer.data(goff + 84 * i + 82);

            auto g_z_zzzzzz = contrBuffer.data(goff + 84 * i + 83);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxxxx, g2_0_xxxxxy,\
                                     g2_0_xxxxxz, g2_0_xxxxyy, g2_0_xxxxyz, g2_0_xxxxzz,\
                                     g2_0_xxxyyy, g2_0_xxxyyz, g2_0_xxxyzz, g2_0_xxxzzz,\
                                     g2_0_xxyyyy, g2_0_xxyyyz, g2_0_xxyyzz, g2_0_xxyzzz,\
                                     g2_0_xxzzzz, g2_0_xyyyyy, g2_0_xyyyyz, g2_0_xyyyzz,\
                                     g2_0_xyyzzz, g2_0_xyzzzz, g2_0_xzzzzz, g2_0_yyyyyy,\
                                     g2_0_yyyyyz, g2_0_yyyyzz, g2_0_yyyzzz, g2_0_yyzzzz,\
                                     g2_0_yzzzzz, g2_0_zzzzzz, g1_0_xxxxxxx, g1_0_xxxxxxy,\
                                     g1_0_xxxxxxz, g1_0_xxxxxyy, g1_0_xxxxxyz,\
                                     g1_0_xxxxxzz, g1_0_xxxxyyy, g1_0_xxxxyyz,\
                                     g1_0_xxxxyzz, g1_0_xxxxzzz, g1_0_xxxyyyy,\
                                     g1_0_xxxyyyz, g1_0_xxxyyzz, g1_0_xxxyzzz,\
                                     g1_0_xxxzzzz, g1_0_xxyyyyy, g1_0_xxyyyyz,\
                                     g1_0_xxyyyzz, g1_0_xxyyzzz, g1_0_xxyzzzz,\
                                     g1_0_xxzzzzz, g1_0_xyyyyyy, g1_0_xyyyyyz,\
                                     g1_0_xyyyyzz, g1_0_xyyyzzz, g1_0_xyyzzzz,\
                                     g1_0_xyzzzzz, g1_0_xzzzzzz, g1_0_yyyyyyy,\
                                     g1_0_yyyyyyz, g1_0_yyyyyzz, g1_0_yyyyzzz,\
                                     g1_0_yyyzzzz, g1_0_yyzzzzz, g1_0_yzzzzzz,\
                                     g1_0_zzzzzzz, g_x_xxxxxx, g_x_xxxxxy, g_x_xxxxxz,\
                                     g_x_xxxxyy, g_x_xxxxyz, g_x_xxxxzz, g_x_xxxyyy,\
                                     g_x_xxxyyz, g_x_xxxyzz, g_x_xxxzzz, g_x_xxyyyy,\
                                     g_x_xxyyyz, g_x_xxyyzz, g_x_xxyzzz, g_x_xxzzzz,\
                                     g_x_xyyyyy, g_x_xyyyyz, g_x_xyyyzz, g_x_xyyzzz,\
                                     g_x_xyzzzz, g_x_xzzzzz, g_x_yyyyyy, g_x_yyyyyz,\
                                     g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz, g_x_yzzzzz,\
                                     g_x_zzzzzz, g_y_xxxxxx, g_y_xxxxxy, g_y_xxxxxz,\
                                     g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxzz, g_y_xxxyyy,\
                                     g_y_xxxyyz, g_y_xxxyzz, g_y_xxxzzz, g_y_xxyyyy,\
                                     g_y_xxyyyz, g_y_xxyyzz, g_y_xxyzzz, g_y_xxzzzz,\
                                     g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyzz, g_y_xyyzzz,\
                                     g_y_xyzzzz, g_y_xzzzzz, g_y_yyyyyy, g_y_yyyyyz,\
                                     g_y_yyyyzz, g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz,\
                                     g_y_zzzzzz, g_z_xxxxxx, g_z_xxxxxy, g_z_xxxxxz,\
                                     g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxzz, g_z_xxxyyy,\
                                     g_z_xxxyyz, g_z_xxxyzz, g_z_xxxzzz, g_z_xxyyyy,\
                                     g_z_xxyyyz, g_z_xxyyzz, g_z_xxyzzz, g_z_xxzzzz,\
                                     g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyzz, g_z_xyyzzz,\
                                     g_z_xyzzzz, g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyz,\
                                     g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz,\
                                     g_z_zzzzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xxxxxx[j] = g1_0_xxxxxxx[j] - fr * g2_0_xxxxxx[j];

                g_x_xxxxxy[j] = g1_0_xxxxxxy[j] - fr * g2_0_xxxxxy[j];

                g_x_xxxxxz[j] = g1_0_xxxxxxz[j] - fr * g2_0_xxxxxz[j];

                g_x_xxxxyy[j] = g1_0_xxxxxyy[j] - fr * g2_0_xxxxyy[j];

                g_x_xxxxyz[j] = g1_0_xxxxxyz[j] - fr * g2_0_xxxxyz[j];

                g_x_xxxxzz[j] = g1_0_xxxxxzz[j] - fr * g2_0_xxxxzz[j];

                g_x_xxxyyy[j] = g1_0_xxxxyyy[j] - fr * g2_0_xxxyyy[j];

                g_x_xxxyyz[j] = g1_0_xxxxyyz[j] - fr * g2_0_xxxyyz[j];

                g_x_xxxyzz[j] = g1_0_xxxxyzz[j] - fr * g2_0_xxxyzz[j];

                g_x_xxxzzz[j] = g1_0_xxxxzzz[j] - fr * g2_0_xxxzzz[j];

                g_x_xxyyyy[j] = g1_0_xxxyyyy[j] - fr * g2_0_xxyyyy[j];

                g_x_xxyyyz[j] = g1_0_xxxyyyz[j] - fr * g2_0_xxyyyz[j];

                g_x_xxyyzz[j] = g1_0_xxxyyzz[j] - fr * g2_0_xxyyzz[j];

                g_x_xxyzzz[j] = g1_0_xxxyzzz[j] - fr * g2_0_xxyzzz[j];

                g_x_xxzzzz[j] = g1_0_xxxzzzz[j] - fr * g2_0_xxzzzz[j];

                g_x_xyyyyy[j] = g1_0_xxyyyyy[j] - fr * g2_0_xyyyyy[j];

                g_x_xyyyyz[j] = g1_0_xxyyyyz[j] - fr * g2_0_xyyyyz[j];

                g_x_xyyyzz[j] = g1_0_xxyyyzz[j] - fr * g2_0_xyyyzz[j];

                g_x_xyyzzz[j] = g1_0_xxyyzzz[j] - fr * g2_0_xyyzzz[j];

                g_x_xyzzzz[j] = g1_0_xxyzzzz[j] - fr * g2_0_xyzzzz[j];

                g_x_xzzzzz[j] = g1_0_xxzzzzz[j] - fr * g2_0_xzzzzz[j];

                g_x_yyyyyy[j] = g1_0_xyyyyyy[j] - fr * g2_0_yyyyyy[j];

                g_x_yyyyyz[j] = g1_0_xyyyyyz[j] - fr * g2_0_yyyyyz[j];

                g_x_yyyyzz[j] = g1_0_xyyyyzz[j] - fr * g2_0_yyyyzz[j];

                g_x_yyyzzz[j] = g1_0_xyyyzzz[j] - fr * g2_0_yyyzzz[j];

                g_x_yyzzzz[j] = g1_0_xyyzzzz[j] - fr * g2_0_yyzzzz[j];

                g_x_yzzzzz[j] = g1_0_xyzzzzz[j] - fr * g2_0_yzzzzz[j];

                g_x_zzzzzz[j] = g1_0_xzzzzzz[j] - fr * g2_0_zzzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xxxxxx[j] = g1_0_xxxxxxy[j] - fr * g2_0_xxxxxx[j];

                g_y_xxxxxy[j] = g1_0_xxxxxyy[j] - fr * g2_0_xxxxxy[j];

                g_y_xxxxxz[j] = g1_0_xxxxxyz[j] - fr * g2_0_xxxxxz[j];

                g_y_xxxxyy[j] = g1_0_xxxxyyy[j] - fr * g2_0_xxxxyy[j];

                g_y_xxxxyz[j] = g1_0_xxxxyyz[j] - fr * g2_0_xxxxyz[j];

                g_y_xxxxzz[j] = g1_0_xxxxyzz[j] - fr * g2_0_xxxxzz[j];

                g_y_xxxyyy[j] = g1_0_xxxyyyy[j] - fr * g2_0_xxxyyy[j];

                g_y_xxxyyz[j] = g1_0_xxxyyyz[j] - fr * g2_0_xxxyyz[j];

                g_y_xxxyzz[j] = g1_0_xxxyyzz[j] - fr * g2_0_xxxyzz[j];

                g_y_xxxzzz[j] = g1_0_xxxyzzz[j] - fr * g2_0_xxxzzz[j];

                g_y_xxyyyy[j] = g1_0_xxyyyyy[j] - fr * g2_0_xxyyyy[j];

                g_y_xxyyyz[j] = g1_0_xxyyyyz[j] - fr * g2_0_xxyyyz[j];

                g_y_xxyyzz[j] = g1_0_xxyyyzz[j] - fr * g2_0_xxyyzz[j];

                g_y_xxyzzz[j] = g1_0_xxyyzzz[j] - fr * g2_0_xxyzzz[j];

                g_y_xxzzzz[j] = g1_0_xxyzzzz[j] - fr * g2_0_xxzzzz[j];

                g_y_xyyyyy[j] = g1_0_xyyyyyy[j] - fr * g2_0_xyyyyy[j];

                g_y_xyyyyz[j] = g1_0_xyyyyyz[j] - fr * g2_0_xyyyyz[j];

                g_y_xyyyzz[j] = g1_0_xyyyyzz[j] - fr * g2_0_xyyyzz[j];

                g_y_xyyzzz[j] = g1_0_xyyyzzz[j] - fr * g2_0_xyyzzz[j];

                g_y_xyzzzz[j] = g1_0_xyyzzzz[j] - fr * g2_0_xyzzzz[j];

                g_y_xzzzzz[j] = g1_0_xyzzzzz[j] - fr * g2_0_xzzzzz[j];

                g_y_yyyyyy[j] = g1_0_yyyyyyy[j] - fr * g2_0_yyyyyy[j];

                g_y_yyyyyz[j] = g1_0_yyyyyyz[j] - fr * g2_0_yyyyyz[j];

                g_y_yyyyzz[j] = g1_0_yyyyyzz[j] - fr * g2_0_yyyyzz[j];

                g_y_yyyzzz[j] = g1_0_yyyyzzz[j] - fr * g2_0_yyyzzz[j];

                g_y_yyzzzz[j] = g1_0_yyyzzzz[j] - fr * g2_0_yyzzzz[j];

                g_y_yzzzzz[j] = g1_0_yyzzzzz[j] - fr * g2_0_yzzzzz[j];

                g_y_zzzzzz[j] = g1_0_yzzzzzz[j] - fr * g2_0_zzzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xxxxxx[j] = g1_0_xxxxxxz[j] - fr * g2_0_xxxxxx[j];

                g_z_xxxxxy[j] = g1_0_xxxxxyz[j] - fr * g2_0_xxxxxy[j];

                g_z_xxxxxz[j] = g1_0_xxxxxzz[j] - fr * g2_0_xxxxxz[j];

                g_z_xxxxyy[j] = g1_0_xxxxyyz[j] - fr * g2_0_xxxxyy[j];

                g_z_xxxxyz[j] = g1_0_xxxxyzz[j] - fr * g2_0_xxxxyz[j];

                g_z_xxxxzz[j] = g1_0_xxxxzzz[j] - fr * g2_0_xxxxzz[j];

                g_z_xxxyyy[j] = g1_0_xxxyyyz[j] - fr * g2_0_xxxyyy[j];

                g_z_xxxyyz[j] = g1_0_xxxyyzz[j] - fr * g2_0_xxxyyz[j];

                g_z_xxxyzz[j] = g1_0_xxxyzzz[j] - fr * g2_0_xxxyzz[j];

                g_z_xxxzzz[j] = g1_0_xxxzzzz[j] - fr * g2_0_xxxzzz[j];

                g_z_xxyyyy[j] = g1_0_xxyyyyz[j] - fr * g2_0_xxyyyy[j];

                g_z_xxyyyz[j] = g1_0_xxyyyzz[j] - fr * g2_0_xxyyyz[j];

                g_z_xxyyzz[j] = g1_0_xxyyzzz[j] - fr * g2_0_xxyyzz[j];

                g_z_xxyzzz[j] = g1_0_xxyzzzz[j] - fr * g2_0_xxyzzz[j];

                g_z_xxzzzz[j] = g1_0_xxzzzzz[j] - fr * g2_0_xxzzzz[j];

                g_z_xyyyyy[j] = g1_0_xyyyyyz[j] - fr * g2_0_xyyyyy[j];

                g_z_xyyyyz[j] = g1_0_xyyyyzz[j] - fr * g2_0_xyyyyz[j];

                g_z_xyyyzz[j] = g1_0_xyyyzzz[j] - fr * g2_0_xyyyzz[j];

                g_z_xyyzzz[j] = g1_0_xyyzzzz[j] - fr * g2_0_xyyzzz[j];

                g_z_xyzzzz[j] = g1_0_xyzzzzz[j] - fr * g2_0_xyzzzz[j];

                g_z_xzzzzz[j] = g1_0_xzzzzzz[j] - fr * g2_0_xzzzzz[j];

                g_z_yyyyyy[j] = g1_0_yyyyyyz[j] - fr * g2_0_yyyyyy[j];

                g_z_yyyyyz[j] = g1_0_yyyyyzz[j] - fr * g2_0_yyyyyz[j];

                g_z_yyyyzz[j] = g1_0_yyyyzzz[j] - fr * g2_0_yyyyzz[j];

                g_z_yyyzzz[j] = g1_0_yyyzzzz[j] - fr * g2_0_yyyzzz[j];

                g_z_yyzzzz[j] = g1_0_yyzzzzz[j] - fr * g2_0_yyzzzz[j];

                g_z_yzzzzz[j] = g1_0_yzzzzzz[j] - fr * g2_0_yzzzzz[j];

                g_z_zzzzzz[j] = g1_0_zzzzzzz[j] - fr * g2_0_zzzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXPK(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 1, 7})) return;

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
                                             {braAngularMomentum, 1, 7});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 8});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 0, 7});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|SK)^(m) integrals

            auto g2_0_xxxxxxx = contrBuffer.data(g2off + 36 * i);

            auto g2_0_xxxxxxy = contrBuffer.data(g2off + 36 * i + 1);

            auto g2_0_xxxxxxz = contrBuffer.data(g2off + 36 * i + 2);

            auto g2_0_xxxxxyy = contrBuffer.data(g2off + 36 * i + 3);

            auto g2_0_xxxxxyz = contrBuffer.data(g2off + 36 * i + 4);

            auto g2_0_xxxxxzz = contrBuffer.data(g2off + 36 * i + 5);

            auto g2_0_xxxxyyy = contrBuffer.data(g2off + 36 * i + 6);

            auto g2_0_xxxxyyz = contrBuffer.data(g2off + 36 * i + 7);

            auto g2_0_xxxxyzz = contrBuffer.data(g2off + 36 * i + 8);

            auto g2_0_xxxxzzz = contrBuffer.data(g2off + 36 * i + 9);

            auto g2_0_xxxyyyy = contrBuffer.data(g2off + 36 * i + 10);

            auto g2_0_xxxyyyz = contrBuffer.data(g2off + 36 * i + 11);

            auto g2_0_xxxyyzz = contrBuffer.data(g2off + 36 * i + 12);

            auto g2_0_xxxyzzz = contrBuffer.data(g2off + 36 * i + 13);

            auto g2_0_xxxzzzz = contrBuffer.data(g2off + 36 * i + 14);

            auto g2_0_xxyyyyy = contrBuffer.data(g2off + 36 * i + 15);

            auto g2_0_xxyyyyz = contrBuffer.data(g2off + 36 * i + 16);

            auto g2_0_xxyyyzz = contrBuffer.data(g2off + 36 * i + 17);

            auto g2_0_xxyyzzz = contrBuffer.data(g2off + 36 * i + 18);

            auto g2_0_xxyzzzz = contrBuffer.data(g2off + 36 * i + 19);

            auto g2_0_xxzzzzz = contrBuffer.data(g2off + 36 * i + 20);

            auto g2_0_xyyyyyy = contrBuffer.data(g2off + 36 * i + 21);

            auto g2_0_xyyyyyz = contrBuffer.data(g2off + 36 * i + 22);

            auto g2_0_xyyyyzz = contrBuffer.data(g2off + 36 * i + 23);

            auto g2_0_xyyyzzz = contrBuffer.data(g2off + 36 * i + 24);

            auto g2_0_xyyzzzz = contrBuffer.data(g2off + 36 * i + 25);

            auto g2_0_xyzzzzz = contrBuffer.data(g2off + 36 * i + 26);

            auto g2_0_xzzzzzz = contrBuffer.data(g2off + 36 * i + 27);

            auto g2_0_yyyyyyy = contrBuffer.data(g2off + 36 * i + 28);

            auto g2_0_yyyyyyz = contrBuffer.data(g2off + 36 * i + 29);

            auto g2_0_yyyyyzz = contrBuffer.data(g2off + 36 * i + 30);

            auto g2_0_yyyyzzz = contrBuffer.data(g2off + 36 * i + 31);

            auto g2_0_yyyzzzz = contrBuffer.data(g2off + 36 * i + 32);

            auto g2_0_yyzzzzz = contrBuffer.data(g2off + 36 * i + 33);

            auto g2_0_yzzzzzz = contrBuffer.data(g2off + 36 * i + 34);

            auto g2_0_zzzzzzz = contrBuffer.data(g2off + 36 * i + 35);

            // set up pointers to (X|g(r,r')|SL)^(m) integrals

            auto g1_0_xxxxxxxx = contrBuffer.data(g1off + 45 * i);

            auto g1_0_xxxxxxxy = contrBuffer.data(g1off + 45 * i + 1);

            auto g1_0_xxxxxxxz = contrBuffer.data(g1off + 45 * i + 2);

            auto g1_0_xxxxxxyy = contrBuffer.data(g1off + 45 * i + 3);

            auto g1_0_xxxxxxyz = contrBuffer.data(g1off + 45 * i + 4);

            auto g1_0_xxxxxxzz = contrBuffer.data(g1off + 45 * i + 5);

            auto g1_0_xxxxxyyy = contrBuffer.data(g1off + 45 * i + 6);

            auto g1_0_xxxxxyyz = contrBuffer.data(g1off + 45 * i + 7);

            auto g1_0_xxxxxyzz = contrBuffer.data(g1off + 45 * i + 8);

            auto g1_0_xxxxxzzz = contrBuffer.data(g1off + 45 * i + 9);

            auto g1_0_xxxxyyyy = contrBuffer.data(g1off + 45 * i + 10);

            auto g1_0_xxxxyyyz = contrBuffer.data(g1off + 45 * i + 11);

            auto g1_0_xxxxyyzz = contrBuffer.data(g1off + 45 * i + 12);

            auto g1_0_xxxxyzzz = contrBuffer.data(g1off + 45 * i + 13);

            auto g1_0_xxxxzzzz = contrBuffer.data(g1off + 45 * i + 14);

            auto g1_0_xxxyyyyy = contrBuffer.data(g1off + 45 * i + 15);

            auto g1_0_xxxyyyyz = contrBuffer.data(g1off + 45 * i + 16);

            auto g1_0_xxxyyyzz = contrBuffer.data(g1off + 45 * i + 17);

            auto g1_0_xxxyyzzz = contrBuffer.data(g1off + 45 * i + 18);

            auto g1_0_xxxyzzzz = contrBuffer.data(g1off + 45 * i + 19);

            auto g1_0_xxxzzzzz = contrBuffer.data(g1off + 45 * i + 20);

            auto g1_0_xxyyyyyy = contrBuffer.data(g1off + 45 * i + 21);

            auto g1_0_xxyyyyyz = contrBuffer.data(g1off + 45 * i + 22);

            auto g1_0_xxyyyyzz = contrBuffer.data(g1off + 45 * i + 23);

            auto g1_0_xxyyyzzz = contrBuffer.data(g1off + 45 * i + 24);

            auto g1_0_xxyyzzzz = contrBuffer.data(g1off + 45 * i + 25);

            auto g1_0_xxyzzzzz = contrBuffer.data(g1off + 45 * i + 26);

            auto g1_0_xxzzzzzz = contrBuffer.data(g1off + 45 * i + 27);

            auto g1_0_xyyyyyyy = contrBuffer.data(g1off + 45 * i + 28);

            auto g1_0_xyyyyyyz = contrBuffer.data(g1off + 45 * i + 29);

            auto g1_0_xyyyyyzz = contrBuffer.data(g1off + 45 * i + 30);

            auto g1_0_xyyyyzzz = contrBuffer.data(g1off + 45 * i + 31);

            auto g1_0_xyyyzzzz = contrBuffer.data(g1off + 45 * i + 32);

            auto g1_0_xyyzzzzz = contrBuffer.data(g1off + 45 * i + 33);

            auto g1_0_xyzzzzzz = contrBuffer.data(g1off + 45 * i + 34);

            auto g1_0_xzzzzzzz = contrBuffer.data(g1off + 45 * i + 35);

            auto g1_0_yyyyyyyy = contrBuffer.data(g1off + 45 * i + 36);

            auto g1_0_yyyyyyyz = contrBuffer.data(g1off + 45 * i + 37);

            auto g1_0_yyyyyyzz = contrBuffer.data(g1off + 45 * i + 38);

            auto g1_0_yyyyyzzz = contrBuffer.data(g1off + 45 * i + 39);

            auto g1_0_yyyyzzzz = contrBuffer.data(g1off + 45 * i + 40);

            auto g1_0_yyyzzzzz = contrBuffer.data(g1off + 45 * i + 41);

            auto g1_0_yyzzzzzz = contrBuffer.data(g1off + 45 * i + 42);

            auto g1_0_yzzzzzzz = contrBuffer.data(g1off + 45 * i + 43);

            auto g1_0_zzzzzzzz = contrBuffer.data(g1off + 45 * i + 44);

            // set up pointers to (X|g(r,r')|PK)^(m) integrals

            auto g_x_xxxxxxx = contrBuffer.data(goff + 108 * i);

            auto g_x_xxxxxxy = contrBuffer.data(goff + 108 * i + 1);

            auto g_x_xxxxxxz = contrBuffer.data(goff + 108 * i + 2);

            auto g_x_xxxxxyy = contrBuffer.data(goff + 108 * i + 3);

            auto g_x_xxxxxyz = contrBuffer.data(goff + 108 * i + 4);

            auto g_x_xxxxxzz = contrBuffer.data(goff + 108 * i + 5);

            auto g_x_xxxxyyy = contrBuffer.data(goff + 108 * i + 6);

            auto g_x_xxxxyyz = contrBuffer.data(goff + 108 * i + 7);

            auto g_x_xxxxyzz = contrBuffer.data(goff + 108 * i + 8);

            auto g_x_xxxxzzz = contrBuffer.data(goff + 108 * i + 9);

            auto g_x_xxxyyyy = contrBuffer.data(goff + 108 * i + 10);

            auto g_x_xxxyyyz = contrBuffer.data(goff + 108 * i + 11);

            auto g_x_xxxyyzz = contrBuffer.data(goff + 108 * i + 12);

            auto g_x_xxxyzzz = contrBuffer.data(goff + 108 * i + 13);

            auto g_x_xxxzzzz = contrBuffer.data(goff + 108 * i + 14);

            auto g_x_xxyyyyy = contrBuffer.data(goff + 108 * i + 15);

            auto g_x_xxyyyyz = contrBuffer.data(goff + 108 * i + 16);

            auto g_x_xxyyyzz = contrBuffer.data(goff + 108 * i + 17);

            auto g_x_xxyyzzz = contrBuffer.data(goff + 108 * i + 18);

            auto g_x_xxyzzzz = contrBuffer.data(goff + 108 * i + 19);

            auto g_x_xxzzzzz = contrBuffer.data(goff + 108 * i + 20);

            auto g_x_xyyyyyy = contrBuffer.data(goff + 108 * i + 21);

            auto g_x_xyyyyyz = contrBuffer.data(goff + 108 * i + 22);

            auto g_x_xyyyyzz = contrBuffer.data(goff + 108 * i + 23);

            auto g_x_xyyyzzz = contrBuffer.data(goff + 108 * i + 24);

            auto g_x_xyyzzzz = contrBuffer.data(goff + 108 * i + 25);

            auto g_x_xyzzzzz = contrBuffer.data(goff + 108 * i + 26);

            auto g_x_xzzzzzz = contrBuffer.data(goff + 108 * i + 27);

            auto g_x_yyyyyyy = contrBuffer.data(goff + 108 * i + 28);

            auto g_x_yyyyyyz = contrBuffer.data(goff + 108 * i + 29);

            auto g_x_yyyyyzz = contrBuffer.data(goff + 108 * i + 30);

            auto g_x_yyyyzzz = contrBuffer.data(goff + 108 * i + 31);

            auto g_x_yyyzzzz = contrBuffer.data(goff + 108 * i + 32);

            auto g_x_yyzzzzz = contrBuffer.data(goff + 108 * i + 33);

            auto g_x_yzzzzzz = contrBuffer.data(goff + 108 * i + 34);

            auto g_x_zzzzzzz = contrBuffer.data(goff + 108 * i + 35);

            auto g_y_xxxxxxx = contrBuffer.data(goff + 108 * i + 36);

            auto g_y_xxxxxxy = contrBuffer.data(goff + 108 * i + 37);

            auto g_y_xxxxxxz = contrBuffer.data(goff + 108 * i + 38);

            auto g_y_xxxxxyy = contrBuffer.data(goff + 108 * i + 39);

            auto g_y_xxxxxyz = contrBuffer.data(goff + 108 * i + 40);

            auto g_y_xxxxxzz = contrBuffer.data(goff + 108 * i + 41);

            auto g_y_xxxxyyy = contrBuffer.data(goff + 108 * i + 42);

            auto g_y_xxxxyyz = contrBuffer.data(goff + 108 * i + 43);

            auto g_y_xxxxyzz = contrBuffer.data(goff + 108 * i + 44);

            auto g_y_xxxxzzz = contrBuffer.data(goff + 108 * i + 45);

            auto g_y_xxxyyyy = contrBuffer.data(goff + 108 * i + 46);

            auto g_y_xxxyyyz = contrBuffer.data(goff + 108 * i + 47);

            auto g_y_xxxyyzz = contrBuffer.data(goff + 108 * i + 48);

            auto g_y_xxxyzzz = contrBuffer.data(goff + 108 * i + 49);

            auto g_y_xxxzzzz = contrBuffer.data(goff + 108 * i + 50);

            auto g_y_xxyyyyy = contrBuffer.data(goff + 108 * i + 51);

            auto g_y_xxyyyyz = contrBuffer.data(goff + 108 * i + 52);

            auto g_y_xxyyyzz = contrBuffer.data(goff + 108 * i + 53);

            auto g_y_xxyyzzz = contrBuffer.data(goff + 108 * i + 54);

            auto g_y_xxyzzzz = contrBuffer.data(goff + 108 * i + 55);

            auto g_y_xxzzzzz = contrBuffer.data(goff + 108 * i + 56);

            auto g_y_xyyyyyy = contrBuffer.data(goff + 108 * i + 57);

            auto g_y_xyyyyyz = contrBuffer.data(goff + 108 * i + 58);

            auto g_y_xyyyyzz = contrBuffer.data(goff + 108 * i + 59);

            auto g_y_xyyyzzz = contrBuffer.data(goff + 108 * i + 60);

            auto g_y_xyyzzzz = contrBuffer.data(goff + 108 * i + 61);

            auto g_y_xyzzzzz = contrBuffer.data(goff + 108 * i + 62);

            auto g_y_xzzzzzz = contrBuffer.data(goff + 108 * i + 63);

            auto g_y_yyyyyyy = contrBuffer.data(goff + 108 * i + 64);

            auto g_y_yyyyyyz = contrBuffer.data(goff + 108 * i + 65);

            auto g_y_yyyyyzz = contrBuffer.data(goff + 108 * i + 66);

            auto g_y_yyyyzzz = contrBuffer.data(goff + 108 * i + 67);

            auto g_y_yyyzzzz = contrBuffer.data(goff + 108 * i + 68);

            auto g_y_yyzzzzz = contrBuffer.data(goff + 108 * i + 69);

            auto g_y_yzzzzzz = contrBuffer.data(goff + 108 * i + 70);

            auto g_y_zzzzzzz = contrBuffer.data(goff + 108 * i + 71);

            auto g_z_xxxxxxx = contrBuffer.data(goff + 108 * i + 72);

            auto g_z_xxxxxxy = contrBuffer.data(goff + 108 * i + 73);

            auto g_z_xxxxxxz = contrBuffer.data(goff + 108 * i + 74);

            auto g_z_xxxxxyy = contrBuffer.data(goff + 108 * i + 75);

            auto g_z_xxxxxyz = contrBuffer.data(goff + 108 * i + 76);

            auto g_z_xxxxxzz = contrBuffer.data(goff + 108 * i + 77);

            auto g_z_xxxxyyy = contrBuffer.data(goff + 108 * i + 78);

            auto g_z_xxxxyyz = contrBuffer.data(goff + 108 * i + 79);

            auto g_z_xxxxyzz = contrBuffer.data(goff + 108 * i + 80);

            auto g_z_xxxxzzz = contrBuffer.data(goff + 108 * i + 81);

            auto g_z_xxxyyyy = contrBuffer.data(goff + 108 * i + 82);

            auto g_z_xxxyyyz = contrBuffer.data(goff + 108 * i + 83);

            auto g_z_xxxyyzz = contrBuffer.data(goff + 108 * i + 84);

            auto g_z_xxxyzzz = contrBuffer.data(goff + 108 * i + 85);

            auto g_z_xxxzzzz = contrBuffer.data(goff + 108 * i + 86);

            auto g_z_xxyyyyy = contrBuffer.data(goff + 108 * i + 87);

            auto g_z_xxyyyyz = contrBuffer.data(goff + 108 * i + 88);

            auto g_z_xxyyyzz = contrBuffer.data(goff + 108 * i + 89);

            auto g_z_xxyyzzz = contrBuffer.data(goff + 108 * i + 90);

            auto g_z_xxyzzzz = contrBuffer.data(goff + 108 * i + 91);

            auto g_z_xxzzzzz = contrBuffer.data(goff + 108 * i + 92);

            auto g_z_xyyyyyy = contrBuffer.data(goff + 108 * i + 93);

            auto g_z_xyyyyyz = contrBuffer.data(goff + 108 * i + 94);

            auto g_z_xyyyyzz = contrBuffer.data(goff + 108 * i + 95);

            auto g_z_xyyyzzz = contrBuffer.data(goff + 108 * i + 96);

            auto g_z_xyyzzzz = contrBuffer.data(goff + 108 * i + 97);

            auto g_z_xyzzzzz = contrBuffer.data(goff + 108 * i + 98);

            auto g_z_xzzzzzz = contrBuffer.data(goff + 108 * i + 99);

            auto g_z_yyyyyyy = contrBuffer.data(goff + 108 * i + 100);

            auto g_z_yyyyyyz = contrBuffer.data(goff + 108 * i + 101);

            auto g_z_yyyyyzz = contrBuffer.data(goff + 108 * i + 102);

            auto g_z_yyyyzzz = contrBuffer.data(goff + 108 * i + 103);

            auto g_z_yyyzzzz = contrBuffer.data(goff + 108 * i + 104);

            auto g_z_yyzzzzz = contrBuffer.data(goff + 108 * i + 105);

            auto g_z_yzzzzzz = contrBuffer.data(goff + 108 * i + 106);

            auto g_z_zzzzzzz = contrBuffer.data(goff + 108 * i + 107);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxxxxx, g2_0_xxxxxxy,\
                                     g2_0_xxxxxxz, g2_0_xxxxxyy, g2_0_xxxxxyz,\
                                     g2_0_xxxxxzz, g2_0_xxxxyyy, g2_0_xxxxyyz,\
                                     g2_0_xxxxyzz, g2_0_xxxxzzz, g2_0_xxxyyyy,\
                                     g2_0_xxxyyyz, g2_0_xxxyyzz, g2_0_xxxyzzz,\
                                     g2_0_xxxzzzz, g2_0_xxyyyyy, g2_0_xxyyyyz,\
                                     g2_0_xxyyyzz, g2_0_xxyyzzz, g2_0_xxyzzzz,\
                                     g2_0_xxzzzzz, g2_0_xyyyyyy, g2_0_xyyyyyz,\
                                     g2_0_xyyyyzz, g2_0_xyyyzzz, g2_0_xyyzzzz,\
                                     g2_0_xyzzzzz, g2_0_xzzzzzz, g2_0_yyyyyyy,\
                                     g2_0_yyyyyyz, g2_0_yyyyyzz, g2_0_yyyyzzz,\
                                     g2_0_yyyzzzz, g2_0_yyzzzzz, g2_0_yzzzzzz,\
                                     g2_0_zzzzzzz, g1_0_xxxxxxxx, g1_0_xxxxxxxy,\
                                     g1_0_xxxxxxxz, g1_0_xxxxxxyy, g1_0_xxxxxxyz,\
                                     g1_0_xxxxxxzz, g1_0_xxxxxyyy, g1_0_xxxxxyyz,\
                                     g1_0_xxxxxyzz, g1_0_xxxxxzzz, g1_0_xxxxyyyy,\
                                     g1_0_xxxxyyyz, g1_0_xxxxyyzz, g1_0_xxxxyzzz,\
                                     g1_0_xxxxzzzz, g1_0_xxxyyyyy, g1_0_xxxyyyyz,\
                                     g1_0_xxxyyyzz, g1_0_xxxyyzzz, g1_0_xxxyzzzz,\
                                     g1_0_xxxzzzzz, g1_0_xxyyyyyy, g1_0_xxyyyyyz,\
                                     g1_0_xxyyyyzz, g1_0_xxyyyzzz, g1_0_xxyyzzzz,\
                                     g1_0_xxyzzzzz, g1_0_xxzzzzzz, g1_0_xyyyyyyy,\
                                     g1_0_xyyyyyyz, g1_0_xyyyyyzz, g1_0_xyyyyzzz,\
                                     g1_0_xyyyzzzz, g1_0_xyyzzzzz, g1_0_xyzzzzzz,\
                                     g1_0_xzzzzzzz, g1_0_yyyyyyyy, g1_0_yyyyyyyz,\
                                     g1_0_yyyyyyzz, g1_0_yyyyyzzz, g1_0_yyyyzzzz,\
                                     g1_0_yyyzzzzz, g1_0_yyzzzzzz, g1_0_yzzzzzzz,\
                                     g1_0_zzzzzzzz, g_x_xxxxxxx, g_x_xxxxxxy, g_x_xxxxxxz,\
                                     g_x_xxxxxyy, g_x_xxxxxyz, g_x_xxxxxzz, g_x_xxxxyyy,\
                                     g_x_xxxxyyz, g_x_xxxxyzz, g_x_xxxxzzz, g_x_xxxyyyy,\
                                     g_x_xxxyyyz, g_x_xxxyyzz, g_x_xxxyzzz, g_x_xxxzzzz,\
                                     g_x_xxyyyyy, g_x_xxyyyyz, g_x_xxyyyzz, g_x_xxyyzzz,\
                                     g_x_xxyzzzz, g_x_xxzzzzz, g_x_xyyyyyy, g_x_xyyyyyz,\
                                     g_x_xyyyyzz, g_x_xyyyzzz, g_x_xyyzzzz, g_x_xyzzzzz,\
                                     g_x_xzzzzzz, g_x_yyyyyyy, g_x_yyyyyyz, g_x_yyyyyzz,\
                                     g_x_yyyyzzz, g_x_yyyzzzz, g_x_yyzzzzz, g_x_yzzzzzz,\
                                     g_x_zzzzzzz, g_y_xxxxxxx, g_y_xxxxxxy, g_y_xxxxxxz,\
                                     g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxzz, g_y_xxxxyyy,\
                                     g_y_xxxxyyz, g_y_xxxxyzz, g_y_xxxxzzz, g_y_xxxyyyy,\
                                     g_y_xxxyyyz, g_y_xxxyyzz, g_y_xxxyzzz, g_y_xxxzzzz,\
                                     g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyzz, g_y_xxyyzzz,\
                                     g_y_xxyzzzz, g_y_xxzzzzz, g_y_xyyyyyy, g_y_xyyyyyz,\
                                     g_y_xyyyyzz, g_y_xyyyzzz, g_y_xyyzzzz, g_y_xyzzzzz,\
                                     g_y_xzzzzzz, g_y_yyyyyyy, g_y_yyyyyyz, g_y_yyyyyzz,\
                                     g_y_yyyyzzz, g_y_yyyzzzz, g_y_yyzzzzz, g_y_yzzzzzz,\
                                     g_y_zzzzzzz, g_z_xxxxxxx, g_z_xxxxxxy, g_z_xxxxxxz,\
                                     g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxzz, g_z_xxxxyyy,\
                                     g_z_xxxxyyz, g_z_xxxxyzz, g_z_xxxxzzz, g_z_xxxyyyy,\
                                     g_z_xxxyyyz, g_z_xxxyyzz, g_z_xxxyzzz, g_z_xxxzzzz,\
                                     g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyzz, g_z_xxyyzzz,\
                                     g_z_xxyzzzz, g_z_xxzzzzz, g_z_xyyyyyy, g_z_xyyyyyz,\
                                     g_z_xyyyyzz, g_z_xyyyzzz, g_z_xyyzzzz, g_z_xyzzzzz,\
                                     g_z_xzzzzzz, g_z_yyyyyyy, g_z_yyyyyyz, g_z_yyyyyzz,\
                                     g_z_yyyyzzz, g_z_yyyzzzz, g_z_yyzzzzz, g_z_yzzzzzz,\
                                     g_z_zzzzzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_x_xxxxxxx[j] = g1_0_xxxxxxxx[j] - fr * g2_0_xxxxxxx[j];

                g_x_xxxxxxy[j] = g1_0_xxxxxxxy[j] - fr * g2_0_xxxxxxy[j];

                g_x_xxxxxxz[j] = g1_0_xxxxxxxz[j] - fr * g2_0_xxxxxxz[j];

                g_x_xxxxxyy[j] = g1_0_xxxxxxyy[j] - fr * g2_0_xxxxxyy[j];

                g_x_xxxxxyz[j] = g1_0_xxxxxxyz[j] - fr * g2_0_xxxxxyz[j];

                g_x_xxxxxzz[j] = g1_0_xxxxxxzz[j] - fr * g2_0_xxxxxzz[j];

                g_x_xxxxyyy[j] = g1_0_xxxxxyyy[j] - fr * g2_0_xxxxyyy[j];

                g_x_xxxxyyz[j] = g1_0_xxxxxyyz[j] - fr * g2_0_xxxxyyz[j];

                g_x_xxxxyzz[j] = g1_0_xxxxxyzz[j] - fr * g2_0_xxxxyzz[j];

                g_x_xxxxzzz[j] = g1_0_xxxxxzzz[j] - fr * g2_0_xxxxzzz[j];

                g_x_xxxyyyy[j] = g1_0_xxxxyyyy[j] - fr * g2_0_xxxyyyy[j];

                g_x_xxxyyyz[j] = g1_0_xxxxyyyz[j] - fr * g2_0_xxxyyyz[j];

                g_x_xxxyyzz[j] = g1_0_xxxxyyzz[j] - fr * g2_0_xxxyyzz[j];

                g_x_xxxyzzz[j] = g1_0_xxxxyzzz[j] - fr * g2_0_xxxyzzz[j];

                g_x_xxxzzzz[j] = g1_0_xxxxzzzz[j] - fr * g2_0_xxxzzzz[j];

                g_x_xxyyyyy[j] = g1_0_xxxyyyyy[j] - fr * g2_0_xxyyyyy[j];

                g_x_xxyyyyz[j] = g1_0_xxxyyyyz[j] - fr * g2_0_xxyyyyz[j];

                g_x_xxyyyzz[j] = g1_0_xxxyyyzz[j] - fr * g2_0_xxyyyzz[j];

                g_x_xxyyzzz[j] = g1_0_xxxyyzzz[j] - fr * g2_0_xxyyzzz[j];

                g_x_xxyzzzz[j] = g1_0_xxxyzzzz[j] - fr * g2_0_xxyzzzz[j];

                g_x_xxzzzzz[j] = g1_0_xxxzzzzz[j] - fr * g2_0_xxzzzzz[j];

                g_x_xyyyyyy[j] = g1_0_xxyyyyyy[j] - fr * g2_0_xyyyyyy[j];

                g_x_xyyyyyz[j] = g1_0_xxyyyyyz[j] - fr * g2_0_xyyyyyz[j];

                g_x_xyyyyzz[j] = g1_0_xxyyyyzz[j] - fr * g2_0_xyyyyzz[j];

                g_x_xyyyzzz[j] = g1_0_xxyyyzzz[j] - fr * g2_0_xyyyzzz[j];

                g_x_xyyzzzz[j] = g1_0_xxyyzzzz[j] - fr * g2_0_xyyzzzz[j];

                g_x_xyzzzzz[j] = g1_0_xxyzzzzz[j] - fr * g2_0_xyzzzzz[j];

                g_x_xzzzzzz[j] = g1_0_xxzzzzzz[j] - fr * g2_0_xzzzzzz[j];

                g_x_yyyyyyy[j] = g1_0_xyyyyyyy[j] - fr * g2_0_yyyyyyy[j];

                g_x_yyyyyyz[j] = g1_0_xyyyyyyz[j] - fr * g2_0_yyyyyyz[j];

                g_x_yyyyyzz[j] = g1_0_xyyyyyzz[j] - fr * g2_0_yyyyyzz[j];

                g_x_yyyyzzz[j] = g1_0_xyyyyzzz[j] - fr * g2_0_yyyyzzz[j];

                g_x_yyyzzzz[j] = g1_0_xyyyzzzz[j] - fr * g2_0_yyyzzzz[j];

                g_x_yyzzzzz[j] = g1_0_xyyzzzzz[j] - fr * g2_0_yyzzzzz[j];

                g_x_yzzzzzz[j] = g1_0_xyzzzzzz[j] - fr * g2_0_yzzzzzz[j];

                g_x_zzzzzzz[j] = g1_0_xzzzzzzz[j] - fr * g2_0_zzzzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_y_xxxxxxx[j] = g1_0_xxxxxxxy[j] - fr * g2_0_xxxxxxx[j];

                g_y_xxxxxxy[j] = g1_0_xxxxxxyy[j] - fr * g2_0_xxxxxxy[j];

                g_y_xxxxxxz[j] = g1_0_xxxxxxyz[j] - fr * g2_0_xxxxxxz[j];

                g_y_xxxxxyy[j] = g1_0_xxxxxyyy[j] - fr * g2_0_xxxxxyy[j];

                g_y_xxxxxyz[j] = g1_0_xxxxxyyz[j] - fr * g2_0_xxxxxyz[j];

                g_y_xxxxxzz[j] = g1_0_xxxxxyzz[j] - fr * g2_0_xxxxxzz[j];

                g_y_xxxxyyy[j] = g1_0_xxxxyyyy[j] - fr * g2_0_xxxxyyy[j];

                g_y_xxxxyyz[j] = g1_0_xxxxyyyz[j] - fr * g2_0_xxxxyyz[j];

                g_y_xxxxyzz[j] = g1_0_xxxxyyzz[j] - fr * g2_0_xxxxyzz[j];

                g_y_xxxxzzz[j] = g1_0_xxxxyzzz[j] - fr * g2_0_xxxxzzz[j];

                g_y_xxxyyyy[j] = g1_0_xxxyyyyy[j] - fr * g2_0_xxxyyyy[j];

                g_y_xxxyyyz[j] = g1_0_xxxyyyyz[j] - fr * g2_0_xxxyyyz[j];

                g_y_xxxyyzz[j] = g1_0_xxxyyyzz[j] - fr * g2_0_xxxyyzz[j];

                g_y_xxxyzzz[j] = g1_0_xxxyyzzz[j] - fr * g2_0_xxxyzzz[j];

                g_y_xxxzzzz[j] = g1_0_xxxyzzzz[j] - fr * g2_0_xxxzzzz[j];

                g_y_xxyyyyy[j] = g1_0_xxyyyyyy[j] - fr * g2_0_xxyyyyy[j];

                g_y_xxyyyyz[j] = g1_0_xxyyyyyz[j] - fr * g2_0_xxyyyyz[j];

                g_y_xxyyyzz[j] = g1_0_xxyyyyzz[j] - fr * g2_0_xxyyyzz[j];

                g_y_xxyyzzz[j] = g1_0_xxyyyzzz[j] - fr * g2_0_xxyyzzz[j];

                g_y_xxyzzzz[j] = g1_0_xxyyzzzz[j] - fr * g2_0_xxyzzzz[j];

                g_y_xxzzzzz[j] = g1_0_xxyzzzzz[j] - fr * g2_0_xxzzzzz[j];

                g_y_xyyyyyy[j] = g1_0_xyyyyyyy[j] - fr * g2_0_xyyyyyy[j];

                g_y_xyyyyyz[j] = g1_0_xyyyyyyz[j] - fr * g2_0_xyyyyyz[j];

                g_y_xyyyyzz[j] = g1_0_xyyyyyzz[j] - fr * g2_0_xyyyyzz[j];

                g_y_xyyyzzz[j] = g1_0_xyyyyzzz[j] - fr * g2_0_xyyyzzz[j];

                g_y_xyyzzzz[j] = g1_0_xyyyzzzz[j] - fr * g2_0_xyyzzzz[j];

                g_y_xyzzzzz[j] = g1_0_xyyzzzzz[j] - fr * g2_0_xyzzzzz[j];

                g_y_xzzzzzz[j] = g1_0_xyzzzzzz[j] - fr * g2_0_xzzzzzz[j];

                g_y_yyyyyyy[j] = g1_0_yyyyyyyy[j] - fr * g2_0_yyyyyyy[j];

                g_y_yyyyyyz[j] = g1_0_yyyyyyyz[j] - fr * g2_0_yyyyyyz[j];

                g_y_yyyyyzz[j] = g1_0_yyyyyyzz[j] - fr * g2_0_yyyyyzz[j];

                g_y_yyyyzzz[j] = g1_0_yyyyyzzz[j] - fr * g2_0_yyyyzzz[j];

                g_y_yyyzzzz[j] = g1_0_yyyyzzzz[j] - fr * g2_0_yyyzzzz[j];

                g_y_yyzzzzz[j] = g1_0_yyyzzzzz[j] - fr * g2_0_yyzzzzz[j];

                g_y_yzzzzzz[j] = g1_0_yyzzzzzz[j] - fr * g2_0_yzzzzzz[j];

                g_y_zzzzzzz[j] = g1_0_yzzzzzzz[j] - fr * g2_0_zzzzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_z_xxxxxxx[j] = g1_0_xxxxxxxz[j] - fr * g2_0_xxxxxxx[j];

                g_z_xxxxxxy[j] = g1_0_xxxxxxyz[j] - fr * g2_0_xxxxxxy[j];

                g_z_xxxxxxz[j] = g1_0_xxxxxxzz[j] - fr * g2_0_xxxxxxz[j];

                g_z_xxxxxyy[j] = g1_0_xxxxxyyz[j] - fr * g2_0_xxxxxyy[j];

                g_z_xxxxxyz[j] = g1_0_xxxxxyzz[j] - fr * g2_0_xxxxxyz[j];

                g_z_xxxxxzz[j] = g1_0_xxxxxzzz[j] - fr * g2_0_xxxxxzz[j];

                g_z_xxxxyyy[j] = g1_0_xxxxyyyz[j] - fr * g2_0_xxxxyyy[j];

                g_z_xxxxyyz[j] = g1_0_xxxxyyzz[j] - fr * g2_0_xxxxyyz[j];

                g_z_xxxxyzz[j] = g1_0_xxxxyzzz[j] - fr * g2_0_xxxxyzz[j];

                g_z_xxxxzzz[j] = g1_0_xxxxzzzz[j] - fr * g2_0_xxxxzzz[j];

                g_z_xxxyyyy[j] = g1_0_xxxyyyyz[j] - fr * g2_0_xxxyyyy[j];

                g_z_xxxyyyz[j] = g1_0_xxxyyyzz[j] - fr * g2_0_xxxyyyz[j];

                g_z_xxxyyzz[j] = g1_0_xxxyyzzz[j] - fr * g2_0_xxxyyzz[j];

                g_z_xxxyzzz[j] = g1_0_xxxyzzzz[j] - fr * g2_0_xxxyzzz[j];

                g_z_xxxzzzz[j] = g1_0_xxxzzzzz[j] - fr * g2_0_xxxzzzz[j];

                g_z_xxyyyyy[j] = g1_0_xxyyyyyz[j] - fr * g2_0_xxyyyyy[j];

                g_z_xxyyyyz[j] = g1_0_xxyyyyzz[j] - fr * g2_0_xxyyyyz[j];

                g_z_xxyyyzz[j] = g1_0_xxyyyzzz[j] - fr * g2_0_xxyyyzz[j];

                g_z_xxyyzzz[j] = g1_0_xxyyzzzz[j] - fr * g2_0_xxyyzzz[j];

                g_z_xxyzzzz[j] = g1_0_xxyzzzzz[j] - fr * g2_0_xxyzzzz[j];

                g_z_xxzzzzz[j] = g1_0_xxzzzzzz[j] - fr * g2_0_xxzzzzz[j];

                g_z_xyyyyyy[j] = g1_0_xyyyyyyz[j] - fr * g2_0_xyyyyyy[j];

                g_z_xyyyyyz[j] = g1_0_xyyyyyzz[j] - fr * g2_0_xyyyyyz[j];

                g_z_xyyyyzz[j] = g1_0_xyyyyzzz[j] - fr * g2_0_xyyyyzz[j];

                g_z_xyyyzzz[j] = g1_0_xyyyzzzz[j] - fr * g2_0_xyyyzzz[j];

                g_z_xyyzzzz[j] = g1_0_xyyzzzzz[j] - fr * g2_0_xyyzzzz[j];

                g_z_xyzzzzz[j] = g1_0_xyzzzzzz[j] - fr * g2_0_xyzzzzz[j];

                g_z_xzzzzzz[j] = g1_0_xzzzzzzz[j] - fr * g2_0_xzzzzzz[j];

                g_z_yyyyyyy[j] = g1_0_yyyyyyyz[j] - fr * g2_0_yyyyyyy[j];

                g_z_yyyyyyz[j] = g1_0_yyyyyyzz[j] - fr * g2_0_yyyyyyz[j];

                g_z_yyyyyzz[j] = g1_0_yyyyyzzz[j] - fr * g2_0_yyyyyzz[j];

                g_z_yyyyzzz[j] = g1_0_yyyyzzzz[j] - fr * g2_0_yyyyzzz[j];

                g_z_yyyzzzz[j] = g1_0_yyyzzzzz[j] - fr * g2_0_yyyzzzz[j];

                g_z_yyzzzzz[j] = g1_0_yyzzzzzz[j] - fr * g2_0_yyzzzzz[j];

                g_z_yzzzzzz[j] = g1_0_yzzzzzzz[j] - fr * g2_0_yzzzzzz[j];

                g_z_zzzzzzz[j] = g1_0_zzzzzzzz[j] - fr * g2_0_zzzzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXDD(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 2, 2})) return;

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
                                             {braAngularMomentum, 2, 2});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 3});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 2});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|PD)^(m) integrals

            auto g2_x_xx = contrBuffer.data(g2off + 18 * i);

            auto g2_x_xy = contrBuffer.data(g2off + 18 * i + 1);

            auto g2_x_xz = contrBuffer.data(g2off + 18 * i + 2);

            auto g2_x_yy = contrBuffer.data(g2off + 18 * i + 3);

            auto g2_x_yz = contrBuffer.data(g2off + 18 * i + 4);

            auto g2_x_zz = contrBuffer.data(g2off + 18 * i + 5);

            auto g2_y_xx = contrBuffer.data(g2off + 18 * i + 6);

            auto g2_y_xy = contrBuffer.data(g2off + 18 * i + 7);

            auto g2_y_xz = contrBuffer.data(g2off + 18 * i + 8);

            auto g2_y_yy = contrBuffer.data(g2off + 18 * i + 9);

            auto g2_y_yz = contrBuffer.data(g2off + 18 * i + 10);

            auto g2_y_zz = contrBuffer.data(g2off + 18 * i + 11);

            auto g2_z_xx = contrBuffer.data(g2off + 18 * i + 12);

            auto g2_z_xy = contrBuffer.data(g2off + 18 * i + 13);

            auto g2_z_xz = contrBuffer.data(g2off + 18 * i + 14);

            auto g2_z_yy = contrBuffer.data(g2off + 18 * i + 15);

            auto g2_z_yz = contrBuffer.data(g2off + 18 * i + 16);

            auto g2_z_zz = contrBuffer.data(g2off + 18 * i + 17);

            // set up pointers to (X|g(r,r')|PF)^(m) integrals

            auto g1_x_xxx = contrBuffer.data(g1off + 30 * i);

            auto g1_x_xxy = contrBuffer.data(g1off + 30 * i + 1);

            auto g1_x_xxz = contrBuffer.data(g1off + 30 * i + 2);

            auto g1_x_xyy = contrBuffer.data(g1off + 30 * i + 3);

            auto g1_x_xyz = contrBuffer.data(g1off + 30 * i + 4);

            auto g1_x_xzz = contrBuffer.data(g1off + 30 * i + 5);

            auto g1_y_xxx = contrBuffer.data(g1off + 30 * i + 10);

            auto g1_y_xxy = contrBuffer.data(g1off + 30 * i + 11);

            auto g1_y_xxz = contrBuffer.data(g1off + 30 * i + 12);

            auto g1_y_xyy = contrBuffer.data(g1off + 30 * i + 13);

            auto g1_y_xyz = contrBuffer.data(g1off + 30 * i + 14);

            auto g1_y_xzz = contrBuffer.data(g1off + 30 * i + 15);

            auto g1_y_yyy = contrBuffer.data(g1off + 30 * i + 16);

            auto g1_y_yyz = contrBuffer.data(g1off + 30 * i + 17);

            auto g1_y_yzz = contrBuffer.data(g1off + 30 * i + 18);

            auto g1_z_xxx = contrBuffer.data(g1off + 30 * i + 20);

            auto g1_z_xxy = contrBuffer.data(g1off + 30 * i + 21);

            auto g1_z_xxz = contrBuffer.data(g1off + 30 * i + 22);

            auto g1_z_xyy = contrBuffer.data(g1off + 30 * i + 23);

            auto g1_z_xyz = contrBuffer.data(g1off + 30 * i + 24);

            auto g1_z_xzz = contrBuffer.data(g1off + 30 * i + 25);

            auto g1_z_yyy = contrBuffer.data(g1off + 30 * i + 26);

            auto g1_z_yyz = contrBuffer.data(g1off + 30 * i + 27);

            auto g1_z_yzz = contrBuffer.data(g1off + 30 * i + 28);

            auto g1_z_zzz = contrBuffer.data(g1off + 30 * i + 29);

            // set up pointers to (X|g(r,r')|DD)^(m) integrals

            auto g_xx_xx = contrBuffer.data(goff + 36 * i);

            auto g_xx_xy = contrBuffer.data(goff + 36 * i + 1);

            auto g_xx_xz = contrBuffer.data(goff + 36 * i + 2);

            auto g_xx_yy = contrBuffer.data(goff + 36 * i + 3);

            auto g_xx_yz = contrBuffer.data(goff + 36 * i + 4);

            auto g_xx_zz = contrBuffer.data(goff + 36 * i + 5);

            auto g_xy_xx = contrBuffer.data(goff + 36 * i + 6);

            auto g_xy_xy = contrBuffer.data(goff + 36 * i + 7);

            auto g_xy_xz = contrBuffer.data(goff + 36 * i + 8);

            auto g_xy_yy = contrBuffer.data(goff + 36 * i + 9);

            auto g_xy_yz = contrBuffer.data(goff + 36 * i + 10);

            auto g_xy_zz = contrBuffer.data(goff + 36 * i + 11);

            auto g_xz_xx = contrBuffer.data(goff + 36 * i + 12);

            auto g_xz_xy = contrBuffer.data(goff + 36 * i + 13);

            auto g_xz_xz = contrBuffer.data(goff + 36 * i + 14);

            auto g_xz_yy = contrBuffer.data(goff + 36 * i + 15);

            auto g_xz_yz = contrBuffer.data(goff + 36 * i + 16);

            auto g_xz_zz = contrBuffer.data(goff + 36 * i + 17);

            auto g_yy_xx = contrBuffer.data(goff + 36 * i + 18);

            auto g_yy_xy = contrBuffer.data(goff + 36 * i + 19);

            auto g_yy_xz = contrBuffer.data(goff + 36 * i + 20);

            auto g_yy_yy = contrBuffer.data(goff + 36 * i + 21);

            auto g_yy_yz = contrBuffer.data(goff + 36 * i + 22);

            auto g_yy_zz = contrBuffer.data(goff + 36 * i + 23);

            auto g_yz_xx = contrBuffer.data(goff + 36 * i + 24);

            auto g_yz_xy = contrBuffer.data(goff + 36 * i + 25);

            auto g_yz_xz = contrBuffer.data(goff + 36 * i + 26);

            auto g_yz_yy = contrBuffer.data(goff + 36 * i + 27);

            auto g_yz_yz = contrBuffer.data(goff + 36 * i + 28);

            auto g_yz_zz = contrBuffer.data(goff + 36 * i + 29);

            auto g_zz_xx = contrBuffer.data(goff + 36 * i + 30);

            auto g_zz_xy = contrBuffer.data(goff + 36 * i + 31);

            auto g_zz_xz = contrBuffer.data(goff + 36 * i + 32);

            auto g_zz_yy = contrBuffer.data(goff + 36 * i + 33);

            auto g_zz_yz = contrBuffer.data(goff + 36 * i + 34);

            auto g_zz_zz = contrBuffer.data(goff + 36 * i + 35);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xx, g2_x_xy, g2_x_xz,\
                                     g2_x_yy, g2_x_yz, g2_x_zz, g2_y_xx, g2_y_xy,\
                                     g2_y_xz, g2_y_yy, g2_y_yz, g2_y_zz, g2_z_xx,\
                                     g2_z_xy, g2_z_xz, g2_z_yy, g2_z_yz, g2_z_zz,\
                                     g1_x_xxx, g1_x_xxy, g1_x_xxz, g1_x_xyy, g1_x_xyz,\
                                     g1_x_xzz, g1_y_xxx, g1_y_xxy, g1_y_xxz, g1_y_xyy,\
                                     g1_y_xyz, g1_y_xzz, g1_y_yyy, g1_y_yyz, g1_y_yzz,\
                                     g1_z_xxx, g1_z_xxy, g1_z_xxz, g1_z_xyy, g1_z_xyz,\
                                     g1_z_xzz, g1_z_yyy, g1_z_yyz, g1_z_yzz, g1_z_zzz,\
                                     g_xx_xx, g_xx_xy, g_xx_xz, g_xx_yy, g_xx_yz,\
                                     g_xx_zz, g_xy_xx, g_xy_xy, g_xy_xz, g_xy_yy,\
                                     g_xy_yz, g_xy_zz, g_xz_xx, g_xz_xy, g_xz_xz,\
                                     g_xz_yy, g_xz_yz, g_xz_zz, g_yy_xx, g_yy_xy,\
                                     g_yy_xz, g_yy_yy, g_yy_yz, g_yy_zz, g_yz_xx,\
                                     g_yz_xy, g_yz_xz, g_yz_yy, g_yz_yz, g_yz_zz,\
                                     g_zz_xx, g_zz_xy, g_zz_xz, g_zz_yy, g_zz_yz,\
                                     g_zz_zz: VLX_ALIGN)
             for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xx_xx[j] = g1_x_xxx[j] - fr * g2_x_xx[j];

                g_xx_xy[j] = g1_x_xxy[j] - fr * g2_x_xy[j];

                g_xx_xz[j] = g1_x_xxz[j] - fr * g2_x_xz[j];

                g_xx_yy[j] = g1_x_xyy[j] - fr * g2_x_yy[j];

                g_xx_yz[j] = g1_x_xyz[j] - fr * g2_x_yz[j];

                g_xx_zz[j] = g1_x_xzz[j] - fr * g2_x_zz[j];

                g_xy_xx[j] = g1_y_xxx[j] - fr * g2_y_xx[j];

                g_xy_xy[j] = g1_y_xxy[j] - fr * g2_y_xy[j];

                g_xy_xz[j] = g1_y_xxz[j] - fr * g2_y_xz[j];

                g_xy_yy[j] = g1_y_xyy[j] - fr * g2_y_yy[j];

                g_xy_yz[j] = g1_y_xyz[j] - fr * g2_y_yz[j];

                g_xy_zz[j] = g1_y_xzz[j] - fr * g2_y_zz[j];

                g_xz_xx[j] = g1_z_xxx[j] - fr * g2_z_xx[j];

                g_xz_xy[j] = g1_z_xxy[j] - fr * g2_z_xy[j];

                g_xz_xz[j] = g1_z_xxz[j] - fr * g2_z_xz[j];

                g_xz_yy[j] = g1_z_xyy[j] - fr * g2_z_yy[j];

                g_xz_yz[j] = g1_z_xyz[j] - fr * g2_z_yz[j];

                g_xz_zz[j] = g1_z_xzz[j] - fr * g2_z_zz[j];

                // leading y component

                fr = rcdy[j];

                g_yy_xx[j] = g1_y_xxy[j] - fr * g2_y_xx[j];

                g_yy_xy[j] = g1_y_xyy[j] - fr * g2_y_xy[j];

                g_yy_xz[j] = g1_y_xyz[j] - fr * g2_y_xz[j];

                g_yy_yy[j] = g1_y_yyy[j] - fr * g2_y_yy[j];

                g_yy_yz[j] = g1_y_yyz[j] - fr * g2_y_yz[j];

                g_yy_zz[j] = g1_y_yzz[j] - fr * g2_y_zz[j];

                g_yz_xx[j] = g1_z_xxy[j] - fr * g2_z_xx[j];

                g_yz_xy[j] = g1_z_xyy[j] - fr * g2_z_xy[j];

                g_yz_xz[j] = g1_z_xyz[j] - fr * g2_z_xz[j];

                g_yz_yy[j] = g1_z_yyy[j] - fr * g2_z_yy[j];

                g_yz_yz[j] = g1_z_yyz[j] - fr * g2_z_yz[j];

                g_yz_zz[j] = g1_z_yzz[j] - fr * g2_z_zz[j];

                // leading z component

                fr = rcdz[j];

                g_zz_xx[j] = g1_z_xxz[j] - fr * g2_z_xx[j];

                g_zz_xy[j] = g1_z_xyz[j] - fr * g2_z_xy[j];

                g_zz_xz[j] = g1_z_xzz[j] - fr * g2_z_xz[j];

                g_zz_yy[j] = g1_z_yyz[j] - fr * g2_z_yy[j];

                g_zz_yz[j] = g1_z_yzz[j] - fr * g2_z_yz[j];

                g_zz_zz[j] = g1_z_zzz[j] - fr * g2_z_zz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXDF(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 2, 3})) return;

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
                                             {braAngularMomentum, 2, 3});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 4});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 3});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|PF)^(m) integrals

            auto g2_x_xxx = contrBuffer.data(g2off + 30 * i);

            auto g2_x_xxy = contrBuffer.data(g2off + 30 * i + 1);

            auto g2_x_xxz = contrBuffer.data(g2off + 30 * i + 2);

            auto g2_x_xyy = contrBuffer.data(g2off + 30 * i + 3);

            auto g2_x_xyz = contrBuffer.data(g2off + 30 * i + 4);

            auto g2_x_xzz = contrBuffer.data(g2off + 30 * i + 5);

            auto g2_x_yyy = contrBuffer.data(g2off + 30 * i + 6);

            auto g2_x_yyz = contrBuffer.data(g2off + 30 * i + 7);

            auto g2_x_yzz = contrBuffer.data(g2off + 30 * i + 8);

            auto g2_x_zzz = contrBuffer.data(g2off + 30 * i + 9);

            auto g2_y_xxx = contrBuffer.data(g2off + 30 * i + 10);

            auto g2_y_xxy = contrBuffer.data(g2off + 30 * i + 11);

            auto g2_y_xxz = contrBuffer.data(g2off + 30 * i + 12);

            auto g2_y_xyy = contrBuffer.data(g2off + 30 * i + 13);

            auto g2_y_xyz = contrBuffer.data(g2off + 30 * i + 14);

            auto g2_y_xzz = contrBuffer.data(g2off + 30 * i + 15);

            auto g2_y_yyy = contrBuffer.data(g2off + 30 * i + 16);

            auto g2_y_yyz = contrBuffer.data(g2off + 30 * i + 17);

            auto g2_y_yzz = contrBuffer.data(g2off + 30 * i + 18);

            auto g2_y_zzz = contrBuffer.data(g2off + 30 * i + 19);

            auto g2_z_xxx = contrBuffer.data(g2off + 30 * i + 20);

            auto g2_z_xxy = contrBuffer.data(g2off + 30 * i + 21);

            auto g2_z_xxz = contrBuffer.data(g2off + 30 * i + 22);

            auto g2_z_xyy = contrBuffer.data(g2off + 30 * i + 23);

            auto g2_z_xyz = contrBuffer.data(g2off + 30 * i + 24);

            auto g2_z_xzz = contrBuffer.data(g2off + 30 * i + 25);

            auto g2_z_yyy = contrBuffer.data(g2off + 30 * i + 26);

            auto g2_z_yyz = contrBuffer.data(g2off + 30 * i + 27);

            auto g2_z_yzz = contrBuffer.data(g2off + 30 * i + 28);

            auto g2_z_zzz = contrBuffer.data(g2off + 30 * i + 29);

            // set up pointers to (X|g(r,r')|PG)^(m) integrals

            auto g1_x_xxxx = contrBuffer.data(g1off + 45 * i);

            auto g1_x_xxxy = contrBuffer.data(g1off + 45 * i + 1);

            auto g1_x_xxxz = contrBuffer.data(g1off + 45 * i + 2);

            auto g1_x_xxyy = contrBuffer.data(g1off + 45 * i + 3);

            auto g1_x_xxyz = contrBuffer.data(g1off + 45 * i + 4);

            auto g1_x_xxzz = contrBuffer.data(g1off + 45 * i + 5);

            auto g1_x_xyyy = contrBuffer.data(g1off + 45 * i + 6);

            auto g1_x_xyyz = contrBuffer.data(g1off + 45 * i + 7);

            auto g1_x_xyzz = contrBuffer.data(g1off + 45 * i + 8);

            auto g1_x_xzzz = contrBuffer.data(g1off + 45 * i + 9);

            auto g1_y_xxxx = contrBuffer.data(g1off + 45 * i + 15);

            auto g1_y_xxxy = contrBuffer.data(g1off + 45 * i + 16);

            auto g1_y_xxxz = contrBuffer.data(g1off + 45 * i + 17);

            auto g1_y_xxyy = contrBuffer.data(g1off + 45 * i + 18);

            auto g1_y_xxyz = contrBuffer.data(g1off + 45 * i + 19);

            auto g1_y_xxzz = contrBuffer.data(g1off + 45 * i + 20);

            auto g1_y_xyyy = contrBuffer.data(g1off + 45 * i + 21);

            auto g1_y_xyyz = contrBuffer.data(g1off + 45 * i + 22);

            auto g1_y_xyzz = contrBuffer.data(g1off + 45 * i + 23);

            auto g1_y_xzzz = contrBuffer.data(g1off + 45 * i + 24);

            auto g1_y_yyyy = contrBuffer.data(g1off + 45 * i + 25);

            auto g1_y_yyyz = contrBuffer.data(g1off + 45 * i + 26);

            auto g1_y_yyzz = contrBuffer.data(g1off + 45 * i + 27);

            auto g1_y_yzzz = contrBuffer.data(g1off + 45 * i + 28);

            auto g1_z_xxxx = contrBuffer.data(g1off + 45 * i + 30);

            auto g1_z_xxxy = contrBuffer.data(g1off + 45 * i + 31);

            auto g1_z_xxxz = contrBuffer.data(g1off + 45 * i + 32);

            auto g1_z_xxyy = contrBuffer.data(g1off + 45 * i + 33);

            auto g1_z_xxyz = contrBuffer.data(g1off + 45 * i + 34);

            auto g1_z_xxzz = contrBuffer.data(g1off + 45 * i + 35);

            auto g1_z_xyyy = contrBuffer.data(g1off + 45 * i + 36);

            auto g1_z_xyyz = contrBuffer.data(g1off + 45 * i + 37);

            auto g1_z_xyzz = contrBuffer.data(g1off + 45 * i + 38);

            auto g1_z_xzzz = contrBuffer.data(g1off + 45 * i + 39);

            auto g1_z_yyyy = contrBuffer.data(g1off + 45 * i + 40);

            auto g1_z_yyyz = contrBuffer.data(g1off + 45 * i + 41);

            auto g1_z_yyzz = contrBuffer.data(g1off + 45 * i + 42);

            auto g1_z_yzzz = contrBuffer.data(g1off + 45 * i + 43);

            auto g1_z_zzzz = contrBuffer.data(g1off + 45 * i + 44);

            // set up pointers to (X|g(r,r')|DF)^(m) integrals

            auto g_xx_xxx = contrBuffer.data(goff + 60 * i);

            auto g_xx_xxy = contrBuffer.data(goff + 60 * i + 1);

            auto g_xx_xxz = contrBuffer.data(goff + 60 * i + 2);

            auto g_xx_xyy = contrBuffer.data(goff + 60 * i + 3);

            auto g_xx_xyz = contrBuffer.data(goff + 60 * i + 4);

            auto g_xx_xzz = contrBuffer.data(goff + 60 * i + 5);

            auto g_xx_yyy = contrBuffer.data(goff + 60 * i + 6);

            auto g_xx_yyz = contrBuffer.data(goff + 60 * i + 7);

            auto g_xx_yzz = contrBuffer.data(goff + 60 * i + 8);

            auto g_xx_zzz = contrBuffer.data(goff + 60 * i + 9);

            auto g_xy_xxx = contrBuffer.data(goff + 60 * i + 10);

            auto g_xy_xxy = contrBuffer.data(goff + 60 * i + 11);

            auto g_xy_xxz = contrBuffer.data(goff + 60 * i + 12);

            auto g_xy_xyy = contrBuffer.data(goff + 60 * i + 13);

            auto g_xy_xyz = contrBuffer.data(goff + 60 * i + 14);

            auto g_xy_xzz = contrBuffer.data(goff + 60 * i + 15);

            auto g_xy_yyy = contrBuffer.data(goff + 60 * i + 16);

            auto g_xy_yyz = contrBuffer.data(goff + 60 * i + 17);

            auto g_xy_yzz = contrBuffer.data(goff + 60 * i + 18);

            auto g_xy_zzz = contrBuffer.data(goff + 60 * i + 19);

            auto g_xz_xxx = contrBuffer.data(goff + 60 * i + 20);

            auto g_xz_xxy = contrBuffer.data(goff + 60 * i + 21);

            auto g_xz_xxz = contrBuffer.data(goff + 60 * i + 22);

            auto g_xz_xyy = contrBuffer.data(goff + 60 * i + 23);

            auto g_xz_xyz = contrBuffer.data(goff + 60 * i + 24);

            auto g_xz_xzz = contrBuffer.data(goff + 60 * i + 25);

            auto g_xz_yyy = contrBuffer.data(goff + 60 * i + 26);

            auto g_xz_yyz = contrBuffer.data(goff + 60 * i + 27);

            auto g_xz_yzz = contrBuffer.data(goff + 60 * i + 28);

            auto g_xz_zzz = contrBuffer.data(goff + 60 * i + 29);

            auto g_yy_xxx = contrBuffer.data(goff + 60 * i + 30);

            auto g_yy_xxy = contrBuffer.data(goff + 60 * i + 31);

            auto g_yy_xxz = contrBuffer.data(goff + 60 * i + 32);

            auto g_yy_xyy = contrBuffer.data(goff + 60 * i + 33);

            auto g_yy_xyz = contrBuffer.data(goff + 60 * i + 34);

            auto g_yy_xzz = contrBuffer.data(goff + 60 * i + 35);

            auto g_yy_yyy = contrBuffer.data(goff + 60 * i + 36);

            auto g_yy_yyz = contrBuffer.data(goff + 60 * i + 37);

            auto g_yy_yzz = contrBuffer.data(goff + 60 * i + 38);

            auto g_yy_zzz = contrBuffer.data(goff + 60 * i + 39);

            auto g_yz_xxx = contrBuffer.data(goff + 60 * i + 40);

            auto g_yz_xxy = contrBuffer.data(goff + 60 * i + 41);

            auto g_yz_xxz = contrBuffer.data(goff + 60 * i + 42);

            auto g_yz_xyy = contrBuffer.data(goff + 60 * i + 43);

            auto g_yz_xyz = contrBuffer.data(goff + 60 * i + 44);

            auto g_yz_xzz = contrBuffer.data(goff + 60 * i + 45);

            auto g_yz_yyy = contrBuffer.data(goff + 60 * i + 46);

            auto g_yz_yyz = contrBuffer.data(goff + 60 * i + 47);

            auto g_yz_yzz = contrBuffer.data(goff + 60 * i + 48);

            auto g_yz_zzz = contrBuffer.data(goff + 60 * i + 49);

            auto g_zz_xxx = contrBuffer.data(goff + 60 * i + 50);

            auto g_zz_xxy = contrBuffer.data(goff + 60 * i + 51);

            auto g_zz_xxz = contrBuffer.data(goff + 60 * i + 52);

            auto g_zz_xyy = contrBuffer.data(goff + 60 * i + 53);

            auto g_zz_xyz = contrBuffer.data(goff + 60 * i + 54);

            auto g_zz_xzz = contrBuffer.data(goff + 60 * i + 55);

            auto g_zz_yyy = contrBuffer.data(goff + 60 * i + 56);

            auto g_zz_yyz = contrBuffer.data(goff + 60 * i + 57);

            auto g_zz_yzz = contrBuffer.data(goff + 60 * i + 58);

            auto g_zz_zzz = contrBuffer.data(goff + 60 * i + 59);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxx, g2_x_xxy, g2_x_xxz,\
                                     g2_x_xyy, g2_x_xyz, g2_x_xzz, g2_x_yyy, g2_x_yyz,\
                                     g2_x_yzz, g2_x_zzz, g2_y_xxx, g2_y_xxy, g2_y_xxz,\
                                     g2_y_xyy, g2_y_xyz, g2_y_xzz, g2_y_yyy, g2_y_yyz,\
                                     g2_y_yzz, g2_y_zzz, g2_z_xxx, g2_z_xxy, g2_z_xxz,\
                                     g2_z_xyy, g2_z_xyz, g2_z_xzz, g2_z_yyy, g2_z_yyz,\
                                     g2_z_yzz, g2_z_zzz, g1_x_xxxx, g1_x_xxxy,\
                                     g1_x_xxxz, g1_x_xxyy, g1_x_xxyz, g1_x_xxzz,\
                                     g1_x_xyyy, g1_x_xyyz, g1_x_xyzz, g1_x_xzzz,\
                                     g1_y_xxxx, g1_y_xxxy, g1_y_xxxz, g1_y_xxyy,\
                                     g1_y_xxyz, g1_y_xxzz, g1_y_xyyy, g1_y_xyyz,\
                                     g1_y_xyzz, g1_y_xzzz, g1_y_yyyy, g1_y_yyyz,\
                                     g1_y_yyzz, g1_y_yzzz, g1_z_xxxx, g1_z_xxxy,\
                                     g1_z_xxxz, g1_z_xxyy, g1_z_xxyz, g1_z_xxzz,\
                                     g1_z_xyyy, g1_z_xyyz, g1_z_xyzz, g1_z_xzzz,\
                                     g1_z_yyyy, g1_z_yyyz, g1_z_yyzz, g1_z_yzzz,\
                                     g1_z_zzzz, g_xx_xxx, g_xx_xxy, g_xx_xxz,\
                                     g_xx_xyy, g_xx_xyz, g_xx_xzz, g_xx_yyy,\
                                     g_xx_yyz, g_xx_yzz, g_xx_zzz, g_xy_xxx,\
                                     g_xy_xxy, g_xy_xxz, g_xy_xyy, g_xy_xyz, g_xy_xzz,\
                                     g_xy_yyy, g_xy_yyz, g_xy_yzz, g_xy_zzz, g_xz_xxx,\
                                     g_xz_xxy, g_xz_xxz, g_xz_xyy, g_xz_xyz, g_xz_xzz,\
                                     g_xz_yyy, g_xz_yyz, g_xz_yzz, g_xz_zzz, g_yy_xxx,\
                                     g_yy_xxy, g_yy_xxz, g_yy_xyy, g_yy_xyz, g_yy_xzz,\
                                     g_yy_yyy, g_yy_yyz, g_yy_yzz, g_yy_zzz, g_yz_xxx,\
                                     g_yz_xxy, g_yz_xxz, g_yz_xyy, g_yz_xyz, g_yz_xzz,\
                                     g_yz_yyy, g_yz_yyz, g_yz_yzz, g_yz_zzz, g_zz_xxx,\
                                     g_zz_xxy, g_zz_xxz, g_zz_xyy, g_zz_xyz, g_zz_xzz,\
                                     g_zz_yyy, g_zz_yyz, g_zz_yzz, g_zz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xx_xxx[j] = g1_x_xxxx[j] - fr * g2_x_xxx[j];

                g_xx_xxy[j] = g1_x_xxxy[j] - fr * g2_x_xxy[j];

                g_xx_xxz[j] = g1_x_xxxz[j] - fr * g2_x_xxz[j];

                g_xx_xyy[j] = g1_x_xxyy[j] - fr * g2_x_xyy[j];

                g_xx_xyz[j] = g1_x_xxyz[j] - fr * g2_x_xyz[j];

                g_xx_xzz[j] = g1_x_xxzz[j] - fr * g2_x_xzz[j];

                g_xx_yyy[j] = g1_x_xyyy[j] - fr * g2_x_yyy[j];

                g_xx_yyz[j] = g1_x_xyyz[j] - fr * g2_x_yyz[j];

                g_xx_yzz[j] = g1_x_xyzz[j] - fr * g2_x_yzz[j];

                g_xx_zzz[j] = g1_x_xzzz[j] - fr * g2_x_zzz[j];

                g_xy_xxx[j] = g1_y_xxxx[j] - fr * g2_y_xxx[j];

                g_xy_xxy[j] = g1_y_xxxy[j] - fr * g2_y_xxy[j];

                g_xy_xxz[j] = g1_y_xxxz[j] - fr * g2_y_xxz[j];

                g_xy_xyy[j] = g1_y_xxyy[j] - fr * g2_y_xyy[j];

                g_xy_xyz[j] = g1_y_xxyz[j] - fr * g2_y_xyz[j];

                g_xy_xzz[j] = g1_y_xxzz[j] - fr * g2_y_xzz[j];

                g_xy_yyy[j] = g1_y_xyyy[j] - fr * g2_y_yyy[j];

                g_xy_yyz[j] = g1_y_xyyz[j] - fr * g2_y_yyz[j];

                g_xy_yzz[j] = g1_y_xyzz[j] - fr * g2_y_yzz[j];

                g_xy_zzz[j] = g1_y_xzzz[j] - fr * g2_y_zzz[j];

                g_xz_xxx[j] = g1_z_xxxx[j] - fr * g2_z_xxx[j];

                g_xz_xxy[j] = g1_z_xxxy[j] - fr * g2_z_xxy[j];

                g_xz_xxz[j] = g1_z_xxxz[j] - fr * g2_z_xxz[j];

                g_xz_xyy[j] = g1_z_xxyy[j] - fr * g2_z_xyy[j];

                g_xz_xyz[j] = g1_z_xxyz[j] - fr * g2_z_xyz[j];

                g_xz_xzz[j] = g1_z_xxzz[j] - fr * g2_z_xzz[j];

                g_xz_yyy[j] = g1_z_xyyy[j] - fr * g2_z_yyy[j];

                g_xz_yyz[j] = g1_z_xyyz[j] - fr * g2_z_yyz[j];

                g_xz_yzz[j] = g1_z_xyzz[j] - fr * g2_z_yzz[j];

                g_xz_zzz[j] = g1_z_xzzz[j] - fr * g2_z_zzz[j];

                // leading y component

                fr = rcdy[j];

                g_yy_xxx[j] = g1_y_xxxy[j] - fr * g2_y_xxx[j];

                g_yy_xxy[j] = g1_y_xxyy[j] - fr * g2_y_xxy[j];

                g_yy_xxz[j] = g1_y_xxyz[j] - fr * g2_y_xxz[j];

                g_yy_xyy[j] = g1_y_xyyy[j] - fr * g2_y_xyy[j];

                g_yy_xyz[j] = g1_y_xyyz[j] - fr * g2_y_xyz[j];

                g_yy_xzz[j] = g1_y_xyzz[j] - fr * g2_y_xzz[j];

                g_yy_yyy[j] = g1_y_yyyy[j] - fr * g2_y_yyy[j];

                g_yy_yyz[j] = g1_y_yyyz[j] - fr * g2_y_yyz[j];

                g_yy_yzz[j] = g1_y_yyzz[j] - fr * g2_y_yzz[j];

                g_yy_zzz[j] = g1_y_yzzz[j] - fr * g2_y_zzz[j];

                g_yz_xxx[j] = g1_z_xxxy[j] - fr * g2_z_xxx[j];

                g_yz_xxy[j] = g1_z_xxyy[j] - fr * g2_z_xxy[j];

                g_yz_xxz[j] = g1_z_xxyz[j] - fr * g2_z_xxz[j];

                g_yz_xyy[j] = g1_z_xyyy[j] - fr * g2_z_xyy[j];

                g_yz_xyz[j] = g1_z_xyyz[j] - fr * g2_z_xyz[j];

                g_yz_xzz[j] = g1_z_xyzz[j] - fr * g2_z_xzz[j];

                g_yz_yyy[j] = g1_z_yyyy[j] - fr * g2_z_yyy[j];

                g_yz_yyz[j] = g1_z_yyyz[j] - fr * g2_z_yyz[j];

                g_yz_yzz[j] = g1_z_yyzz[j] - fr * g2_z_yzz[j];

                g_yz_zzz[j] = g1_z_yzzz[j] - fr * g2_z_zzz[j];

                // leading z component

                fr = rcdz[j];

                g_zz_xxx[j] = g1_z_xxxz[j] - fr * g2_z_xxx[j];

                g_zz_xxy[j] = g1_z_xxyz[j] - fr * g2_z_xxy[j];

                g_zz_xxz[j] = g1_z_xxzz[j] - fr * g2_z_xxz[j];

                g_zz_xyy[j] = g1_z_xyyz[j] - fr * g2_z_xyy[j];

                g_zz_xyz[j] = g1_z_xyzz[j] - fr * g2_z_xyz[j];

                g_zz_xzz[j] = g1_z_xzzz[j] - fr * g2_z_xzz[j];

                g_zz_yyy[j] = g1_z_yyyz[j] - fr * g2_z_yyy[j];

                g_zz_yyz[j] = g1_z_yyzz[j] - fr * g2_z_yyz[j];

                g_zz_yzz[j] = g1_z_yzzz[j] - fr * g2_z_yzz[j];

                g_zz_zzz[j] = g1_z_zzzz[j] - fr * g2_z_zzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXDG(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 2, 4})) return;

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
                                             {braAngularMomentum, 2, 4});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 5});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 4});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|PG)^(m) integrals

            auto g2_x_xxxx = contrBuffer.data(g2off + 45 * i);

            auto g2_x_xxxy = contrBuffer.data(g2off + 45 * i + 1);

            auto g2_x_xxxz = contrBuffer.data(g2off + 45 * i + 2);

            auto g2_x_xxyy = contrBuffer.data(g2off + 45 * i + 3);

            auto g2_x_xxyz = contrBuffer.data(g2off + 45 * i + 4);

            auto g2_x_xxzz = contrBuffer.data(g2off + 45 * i + 5);

            auto g2_x_xyyy = contrBuffer.data(g2off + 45 * i + 6);

            auto g2_x_xyyz = contrBuffer.data(g2off + 45 * i + 7);

            auto g2_x_xyzz = contrBuffer.data(g2off + 45 * i + 8);

            auto g2_x_xzzz = contrBuffer.data(g2off + 45 * i + 9);

            auto g2_x_yyyy = contrBuffer.data(g2off + 45 * i + 10);

            auto g2_x_yyyz = contrBuffer.data(g2off + 45 * i + 11);

            auto g2_x_yyzz = contrBuffer.data(g2off + 45 * i + 12);

            auto g2_x_yzzz = contrBuffer.data(g2off + 45 * i + 13);

            auto g2_x_zzzz = contrBuffer.data(g2off + 45 * i + 14);

            auto g2_y_xxxx = contrBuffer.data(g2off + 45 * i + 15);

            auto g2_y_xxxy = contrBuffer.data(g2off + 45 * i + 16);

            auto g2_y_xxxz = contrBuffer.data(g2off + 45 * i + 17);

            auto g2_y_xxyy = contrBuffer.data(g2off + 45 * i + 18);

            auto g2_y_xxyz = contrBuffer.data(g2off + 45 * i + 19);

            auto g2_y_xxzz = contrBuffer.data(g2off + 45 * i + 20);

            auto g2_y_xyyy = contrBuffer.data(g2off + 45 * i + 21);

            auto g2_y_xyyz = contrBuffer.data(g2off + 45 * i + 22);

            auto g2_y_xyzz = contrBuffer.data(g2off + 45 * i + 23);

            auto g2_y_xzzz = contrBuffer.data(g2off + 45 * i + 24);

            auto g2_y_yyyy = contrBuffer.data(g2off + 45 * i + 25);

            auto g2_y_yyyz = contrBuffer.data(g2off + 45 * i + 26);

            auto g2_y_yyzz = contrBuffer.data(g2off + 45 * i + 27);

            auto g2_y_yzzz = contrBuffer.data(g2off + 45 * i + 28);

            auto g2_y_zzzz = contrBuffer.data(g2off + 45 * i + 29);

            auto g2_z_xxxx = contrBuffer.data(g2off + 45 * i + 30);

            auto g2_z_xxxy = contrBuffer.data(g2off + 45 * i + 31);

            auto g2_z_xxxz = contrBuffer.data(g2off + 45 * i + 32);

            auto g2_z_xxyy = contrBuffer.data(g2off + 45 * i + 33);

            auto g2_z_xxyz = contrBuffer.data(g2off + 45 * i + 34);

            auto g2_z_xxzz = contrBuffer.data(g2off + 45 * i + 35);

            auto g2_z_xyyy = contrBuffer.data(g2off + 45 * i + 36);

            auto g2_z_xyyz = contrBuffer.data(g2off + 45 * i + 37);

            auto g2_z_xyzz = contrBuffer.data(g2off + 45 * i + 38);

            auto g2_z_xzzz = contrBuffer.data(g2off + 45 * i + 39);

            auto g2_z_yyyy = contrBuffer.data(g2off + 45 * i + 40);

            auto g2_z_yyyz = contrBuffer.data(g2off + 45 * i + 41);

            auto g2_z_yyzz = contrBuffer.data(g2off + 45 * i + 42);

            auto g2_z_yzzz = contrBuffer.data(g2off + 45 * i + 43);

            auto g2_z_zzzz = contrBuffer.data(g2off + 45 * i + 44);

            // set up pointers to (X|g(r,r')|PH)^(m) integrals

            auto g1_x_xxxxx = contrBuffer.data(g1off + 63 * i);

            auto g1_x_xxxxy = contrBuffer.data(g1off + 63 * i + 1);

            auto g1_x_xxxxz = contrBuffer.data(g1off + 63 * i + 2);

            auto g1_x_xxxyy = contrBuffer.data(g1off + 63 * i + 3);

            auto g1_x_xxxyz = contrBuffer.data(g1off + 63 * i + 4);

            auto g1_x_xxxzz = contrBuffer.data(g1off + 63 * i + 5);

            auto g1_x_xxyyy = contrBuffer.data(g1off + 63 * i + 6);

            auto g1_x_xxyyz = contrBuffer.data(g1off + 63 * i + 7);

            auto g1_x_xxyzz = contrBuffer.data(g1off + 63 * i + 8);

            auto g1_x_xxzzz = contrBuffer.data(g1off + 63 * i + 9);

            auto g1_x_xyyyy = contrBuffer.data(g1off + 63 * i + 10);

            auto g1_x_xyyyz = contrBuffer.data(g1off + 63 * i + 11);

            auto g1_x_xyyzz = contrBuffer.data(g1off + 63 * i + 12);

            auto g1_x_xyzzz = contrBuffer.data(g1off + 63 * i + 13);

            auto g1_x_xzzzz = contrBuffer.data(g1off + 63 * i + 14);

            auto g1_y_xxxxx = contrBuffer.data(g1off + 63 * i + 21);

            auto g1_y_xxxxy = contrBuffer.data(g1off + 63 * i + 22);

            auto g1_y_xxxxz = contrBuffer.data(g1off + 63 * i + 23);

            auto g1_y_xxxyy = contrBuffer.data(g1off + 63 * i + 24);

            auto g1_y_xxxyz = contrBuffer.data(g1off + 63 * i + 25);

            auto g1_y_xxxzz = contrBuffer.data(g1off + 63 * i + 26);

            auto g1_y_xxyyy = contrBuffer.data(g1off + 63 * i + 27);

            auto g1_y_xxyyz = contrBuffer.data(g1off + 63 * i + 28);

            auto g1_y_xxyzz = contrBuffer.data(g1off + 63 * i + 29);

            auto g1_y_xxzzz = contrBuffer.data(g1off + 63 * i + 30);

            auto g1_y_xyyyy = contrBuffer.data(g1off + 63 * i + 31);

            auto g1_y_xyyyz = contrBuffer.data(g1off + 63 * i + 32);

            auto g1_y_xyyzz = contrBuffer.data(g1off + 63 * i + 33);

            auto g1_y_xyzzz = contrBuffer.data(g1off + 63 * i + 34);

            auto g1_y_xzzzz = contrBuffer.data(g1off + 63 * i + 35);

            auto g1_y_yyyyy = contrBuffer.data(g1off + 63 * i + 36);

            auto g1_y_yyyyz = contrBuffer.data(g1off + 63 * i + 37);

            auto g1_y_yyyzz = contrBuffer.data(g1off + 63 * i + 38);

            auto g1_y_yyzzz = contrBuffer.data(g1off + 63 * i + 39);

            auto g1_y_yzzzz = contrBuffer.data(g1off + 63 * i + 40);

            auto g1_z_xxxxx = contrBuffer.data(g1off + 63 * i + 42);

            auto g1_z_xxxxy = contrBuffer.data(g1off + 63 * i + 43);

            auto g1_z_xxxxz = contrBuffer.data(g1off + 63 * i + 44);

            auto g1_z_xxxyy = contrBuffer.data(g1off + 63 * i + 45);

            auto g1_z_xxxyz = contrBuffer.data(g1off + 63 * i + 46);

            auto g1_z_xxxzz = contrBuffer.data(g1off + 63 * i + 47);

            auto g1_z_xxyyy = contrBuffer.data(g1off + 63 * i + 48);

            auto g1_z_xxyyz = contrBuffer.data(g1off + 63 * i + 49);

            auto g1_z_xxyzz = contrBuffer.data(g1off + 63 * i + 50);

            auto g1_z_xxzzz = contrBuffer.data(g1off + 63 * i + 51);

            auto g1_z_xyyyy = contrBuffer.data(g1off + 63 * i + 52);

            auto g1_z_xyyyz = contrBuffer.data(g1off + 63 * i + 53);

            auto g1_z_xyyzz = contrBuffer.data(g1off + 63 * i + 54);

            auto g1_z_xyzzz = contrBuffer.data(g1off + 63 * i + 55);

            auto g1_z_xzzzz = contrBuffer.data(g1off + 63 * i + 56);

            auto g1_z_yyyyy = contrBuffer.data(g1off + 63 * i + 57);

            auto g1_z_yyyyz = contrBuffer.data(g1off + 63 * i + 58);

            auto g1_z_yyyzz = contrBuffer.data(g1off + 63 * i + 59);

            auto g1_z_yyzzz = contrBuffer.data(g1off + 63 * i + 60);

            auto g1_z_yzzzz = contrBuffer.data(g1off + 63 * i + 61);

            auto g1_z_zzzzz = contrBuffer.data(g1off + 63 * i + 62);

            // set up pointers to (X|g(r,r')|DG)^(m) integrals

            auto g_xx_xxxx = contrBuffer.data(goff + 90 * i);

            auto g_xx_xxxy = contrBuffer.data(goff + 90 * i + 1);

            auto g_xx_xxxz = contrBuffer.data(goff + 90 * i + 2);

            auto g_xx_xxyy = contrBuffer.data(goff + 90 * i + 3);

            auto g_xx_xxyz = contrBuffer.data(goff + 90 * i + 4);

            auto g_xx_xxzz = contrBuffer.data(goff + 90 * i + 5);

            auto g_xx_xyyy = contrBuffer.data(goff + 90 * i + 6);

            auto g_xx_xyyz = contrBuffer.data(goff + 90 * i + 7);

            auto g_xx_xyzz = contrBuffer.data(goff + 90 * i + 8);

            auto g_xx_xzzz = contrBuffer.data(goff + 90 * i + 9);

            auto g_xx_yyyy = contrBuffer.data(goff + 90 * i + 10);

            auto g_xx_yyyz = contrBuffer.data(goff + 90 * i + 11);

            auto g_xx_yyzz = contrBuffer.data(goff + 90 * i + 12);

            auto g_xx_yzzz = contrBuffer.data(goff + 90 * i + 13);

            auto g_xx_zzzz = contrBuffer.data(goff + 90 * i + 14);

            auto g_xy_xxxx = contrBuffer.data(goff + 90 * i + 15);

            auto g_xy_xxxy = contrBuffer.data(goff + 90 * i + 16);

            auto g_xy_xxxz = contrBuffer.data(goff + 90 * i + 17);

            auto g_xy_xxyy = contrBuffer.data(goff + 90 * i + 18);

            auto g_xy_xxyz = contrBuffer.data(goff + 90 * i + 19);

            auto g_xy_xxzz = contrBuffer.data(goff + 90 * i + 20);

            auto g_xy_xyyy = contrBuffer.data(goff + 90 * i + 21);

            auto g_xy_xyyz = contrBuffer.data(goff + 90 * i + 22);

            auto g_xy_xyzz = contrBuffer.data(goff + 90 * i + 23);

            auto g_xy_xzzz = contrBuffer.data(goff + 90 * i + 24);

            auto g_xy_yyyy = contrBuffer.data(goff + 90 * i + 25);

            auto g_xy_yyyz = contrBuffer.data(goff + 90 * i + 26);

            auto g_xy_yyzz = contrBuffer.data(goff + 90 * i + 27);

            auto g_xy_yzzz = contrBuffer.data(goff + 90 * i + 28);

            auto g_xy_zzzz = contrBuffer.data(goff + 90 * i + 29);

            auto g_xz_xxxx = contrBuffer.data(goff + 90 * i + 30);

            auto g_xz_xxxy = contrBuffer.data(goff + 90 * i + 31);

            auto g_xz_xxxz = contrBuffer.data(goff + 90 * i + 32);

            auto g_xz_xxyy = contrBuffer.data(goff + 90 * i + 33);

            auto g_xz_xxyz = contrBuffer.data(goff + 90 * i + 34);

            auto g_xz_xxzz = contrBuffer.data(goff + 90 * i + 35);

            auto g_xz_xyyy = contrBuffer.data(goff + 90 * i + 36);

            auto g_xz_xyyz = contrBuffer.data(goff + 90 * i + 37);

            auto g_xz_xyzz = contrBuffer.data(goff + 90 * i + 38);

            auto g_xz_xzzz = contrBuffer.data(goff + 90 * i + 39);

            auto g_xz_yyyy = contrBuffer.data(goff + 90 * i + 40);

            auto g_xz_yyyz = contrBuffer.data(goff + 90 * i + 41);

            auto g_xz_yyzz = contrBuffer.data(goff + 90 * i + 42);

            auto g_xz_yzzz = contrBuffer.data(goff + 90 * i + 43);

            auto g_xz_zzzz = contrBuffer.data(goff + 90 * i + 44);

            auto g_yy_xxxx = contrBuffer.data(goff + 90 * i + 45);

            auto g_yy_xxxy = contrBuffer.data(goff + 90 * i + 46);

            auto g_yy_xxxz = contrBuffer.data(goff + 90 * i + 47);

            auto g_yy_xxyy = contrBuffer.data(goff + 90 * i + 48);

            auto g_yy_xxyz = contrBuffer.data(goff + 90 * i + 49);

            auto g_yy_xxzz = contrBuffer.data(goff + 90 * i + 50);

            auto g_yy_xyyy = contrBuffer.data(goff + 90 * i + 51);

            auto g_yy_xyyz = contrBuffer.data(goff + 90 * i + 52);

            auto g_yy_xyzz = contrBuffer.data(goff + 90 * i + 53);

            auto g_yy_xzzz = contrBuffer.data(goff + 90 * i + 54);

            auto g_yy_yyyy = contrBuffer.data(goff + 90 * i + 55);

            auto g_yy_yyyz = contrBuffer.data(goff + 90 * i + 56);

            auto g_yy_yyzz = contrBuffer.data(goff + 90 * i + 57);

            auto g_yy_yzzz = contrBuffer.data(goff + 90 * i + 58);

            auto g_yy_zzzz = contrBuffer.data(goff + 90 * i + 59);

            auto g_yz_xxxx = contrBuffer.data(goff + 90 * i + 60);

            auto g_yz_xxxy = contrBuffer.data(goff + 90 * i + 61);

            auto g_yz_xxxz = contrBuffer.data(goff + 90 * i + 62);

            auto g_yz_xxyy = contrBuffer.data(goff + 90 * i + 63);

            auto g_yz_xxyz = contrBuffer.data(goff + 90 * i + 64);

            auto g_yz_xxzz = contrBuffer.data(goff + 90 * i + 65);

            auto g_yz_xyyy = contrBuffer.data(goff + 90 * i + 66);

            auto g_yz_xyyz = contrBuffer.data(goff + 90 * i + 67);

            auto g_yz_xyzz = contrBuffer.data(goff + 90 * i + 68);

            auto g_yz_xzzz = contrBuffer.data(goff + 90 * i + 69);

            auto g_yz_yyyy = contrBuffer.data(goff + 90 * i + 70);

            auto g_yz_yyyz = contrBuffer.data(goff + 90 * i + 71);

            auto g_yz_yyzz = contrBuffer.data(goff + 90 * i + 72);

            auto g_yz_yzzz = contrBuffer.data(goff + 90 * i + 73);

            auto g_yz_zzzz = contrBuffer.data(goff + 90 * i + 74);

            auto g_zz_xxxx = contrBuffer.data(goff + 90 * i + 75);

            auto g_zz_xxxy = contrBuffer.data(goff + 90 * i + 76);

            auto g_zz_xxxz = contrBuffer.data(goff + 90 * i + 77);

            auto g_zz_xxyy = contrBuffer.data(goff + 90 * i + 78);

            auto g_zz_xxyz = contrBuffer.data(goff + 90 * i + 79);

            auto g_zz_xxzz = contrBuffer.data(goff + 90 * i + 80);

            auto g_zz_xyyy = contrBuffer.data(goff + 90 * i + 81);

            auto g_zz_xyyz = contrBuffer.data(goff + 90 * i + 82);

            auto g_zz_xyzz = contrBuffer.data(goff + 90 * i + 83);

            auto g_zz_xzzz = contrBuffer.data(goff + 90 * i + 84);

            auto g_zz_yyyy = contrBuffer.data(goff + 90 * i + 85);

            auto g_zz_yyyz = contrBuffer.data(goff + 90 * i + 86);

            auto g_zz_yyzz = contrBuffer.data(goff + 90 * i + 87);

            auto g_zz_yzzz = contrBuffer.data(goff + 90 * i + 88);

            auto g_zz_zzzz = contrBuffer.data(goff + 90 * i + 89);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxx, g2_x_xxxy, g2_x_xxxz,\
                                     g2_x_xxyy, g2_x_xxyz, g2_x_xxzz, g2_x_xyyy,\
                                     g2_x_xyyz, g2_x_xyzz, g2_x_xzzz, g2_x_yyyy,\
                                     g2_x_yyyz, g2_x_yyzz, g2_x_yzzz, g2_x_zzzz,\
                                     g2_y_xxxx, g2_y_xxxy, g2_y_xxxz, g2_y_xxyy,\
                                     g2_y_xxyz, g2_y_xxzz, g2_y_xyyy, g2_y_xyyz,\
                                     g2_y_xyzz, g2_y_xzzz, g2_y_yyyy, g2_y_yyyz,\
                                     g2_y_yyzz, g2_y_yzzz, g2_y_zzzz, g2_z_xxxx,\
                                     g2_z_xxxy, g2_z_xxxz, g2_z_xxyy, g2_z_xxyz,\
                                     g2_z_xxzz, g2_z_xyyy, g2_z_xyyz, g2_z_xyzz,\
                                     g2_z_xzzz, g2_z_yyyy, g2_z_yyyz, g2_z_yyzz,\
                                     g2_z_yzzz, g2_z_zzzz, g1_x_xxxxx, g1_x_xxxxy,\
                                     g1_x_xxxxz, g1_x_xxxyy, g1_x_xxxyz, g1_x_xxxzz,\
                                     g1_x_xxyyy, g1_x_xxyyz, g1_x_xxyzz, g1_x_xxzzz,\
                                     g1_x_xyyyy, g1_x_xyyyz, g1_x_xyyzz, g1_x_xyzzz,\
                                     g1_x_xzzzz, g1_y_xxxxx, g1_y_xxxxy, g1_y_xxxxz,\
                                     g1_y_xxxyy, g1_y_xxxyz, g1_y_xxxzz, g1_y_xxyyy,\
                                     g1_y_xxyyz, g1_y_xxyzz, g1_y_xxzzz, g1_y_xyyyy,\
                                     g1_y_xyyyz, g1_y_xyyzz, g1_y_xyzzz, g1_y_xzzzz,\
                                     g1_y_yyyyy, g1_y_yyyyz, g1_y_yyyzz, g1_y_yyzzz,\
                                     g1_y_yzzzz, g1_z_xxxxx, g1_z_xxxxy, g1_z_xxxxz,\
                                     g1_z_xxxyy, g1_z_xxxyz, g1_z_xxxzz, g1_z_xxyyy,\
                                     g1_z_xxyyz, g1_z_xxyzz, g1_z_xxzzz, g1_z_xyyyy,\
                                     g1_z_xyyyz, g1_z_xyyzz, g1_z_xyzzz, g1_z_xzzzz,\
                                     g1_z_yyyyy, g1_z_yyyyz, g1_z_yyyzz, g1_z_yyzzz,\
                                     g1_z_yzzzz, g1_z_zzzzz, g_xx_xxxx, g_xx_xxxy,\
                                     g_xx_xxxz, g_xx_xxyy, g_xx_xxyz, g_xx_xxzz,\
                                     g_xx_xyyy, g_xx_xyyz, g_xx_xyzz, g_xx_xzzz,\
                                     g_xx_yyyy, g_xx_yyyz, g_xx_yyzz, g_xx_yzzz,\
                                     g_xx_zzzz, g_xy_xxxx, g_xy_xxxy, g_xy_xxxz,\
                                     g_xy_xxyy, g_xy_xxyz, g_xy_xxzz, g_xy_xyyy,\
                                     g_xy_xyyz, g_xy_xyzz, g_xy_xzzz, g_xy_yyyy,\
                                     g_xy_yyyz, g_xy_yyzz, g_xy_yzzz, g_xy_zzzz,\
                                     g_xz_xxxx, g_xz_xxxy, g_xz_xxxz, g_xz_xxyy,\
                                     g_xz_xxyz, g_xz_xxzz, g_xz_xyyy, g_xz_xyyz,\
                                     g_xz_xyzz, g_xz_xzzz, g_xz_yyyy, g_xz_yyyz,\
                                     g_xz_yyzz, g_xz_yzzz, g_xz_zzzz, g_yy_xxxx,\
                                     g_yy_xxxy, g_yy_xxxz, g_yy_xxyy, g_yy_xxyz,\
                                     g_yy_xxzz, g_yy_xyyy, g_yy_xyyz, g_yy_xyzz,\
                                     g_yy_xzzz, g_yy_yyyy, g_yy_yyyz, g_yy_yyzz,\
                                     g_yy_yzzz, g_yy_zzzz, g_yz_xxxx, g_yz_xxxy,\
                                     g_yz_xxxz, g_yz_xxyy, g_yz_xxyz, g_yz_xxzz,\
                                     g_yz_xyyy, g_yz_xyyz, g_yz_xyzz, g_yz_xzzz,\
                                     g_yz_yyyy, g_yz_yyyz, g_yz_yyzz, g_yz_yzzz,\
                                     g_yz_zzzz, g_zz_xxxx, g_zz_xxxy, g_zz_xxxz,\
                                     g_zz_xxyy, g_zz_xxyz, g_zz_xxzz, g_zz_xyyy,\
                                     g_zz_xyyz, g_zz_xyzz, g_zz_xzzz, g_zz_yyyy,\
                                     g_zz_yyyz, g_zz_yyzz, g_zz_yzzz, g_zz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xx_xxxx[j] = g1_x_xxxxx[j] - fr * g2_x_xxxx[j];

                g_xx_xxxy[j] = g1_x_xxxxy[j] - fr * g2_x_xxxy[j];

                g_xx_xxxz[j] = g1_x_xxxxz[j] - fr * g2_x_xxxz[j];

                g_xx_xxyy[j] = g1_x_xxxyy[j] - fr * g2_x_xxyy[j];

                g_xx_xxyz[j] = g1_x_xxxyz[j] - fr * g2_x_xxyz[j];

                g_xx_xxzz[j] = g1_x_xxxzz[j] - fr * g2_x_xxzz[j];

                g_xx_xyyy[j] = g1_x_xxyyy[j] - fr * g2_x_xyyy[j];

                g_xx_xyyz[j] = g1_x_xxyyz[j] - fr * g2_x_xyyz[j];

                g_xx_xyzz[j] = g1_x_xxyzz[j] - fr * g2_x_xyzz[j];

                g_xx_xzzz[j] = g1_x_xxzzz[j] - fr * g2_x_xzzz[j];

                g_xx_yyyy[j] = g1_x_xyyyy[j] - fr * g2_x_yyyy[j];

                g_xx_yyyz[j] = g1_x_xyyyz[j] - fr * g2_x_yyyz[j];

                g_xx_yyzz[j] = g1_x_xyyzz[j] - fr * g2_x_yyzz[j];

                g_xx_yzzz[j] = g1_x_xyzzz[j] - fr * g2_x_yzzz[j];

                g_xx_zzzz[j] = g1_x_xzzzz[j] - fr * g2_x_zzzz[j];

                g_xy_xxxx[j] = g1_y_xxxxx[j] - fr * g2_y_xxxx[j];

                g_xy_xxxy[j] = g1_y_xxxxy[j] - fr * g2_y_xxxy[j];

                g_xy_xxxz[j] = g1_y_xxxxz[j] - fr * g2_y_xxxz[j];

                g_xy_xxyy[j] = g1_y_xxxyy[j] - fr * g2_y_xxyy[j];

                g_xy_xxyz[j] = g1_y_xxxyz[j] - fr * g2_y_xxyz[j];

                g_xy_xxzz[j] = g1_y_xxxzz[j] - fr * g2_y_xxzz[j];

                g_xy_xyyy[j] = g1_y_xxyyy[j] - fr * g2_y_xyyy[j];

                g_xy_xyyz[j] = g1_y_xxyyz[j] - fr * g2_y_xyyz[j];

                g_xy_xyzz[j] = g1_y_xxyzz[j] - fr * g2_y_xyzz[j];

                g_xy_xzzz[j] = g1_y_xxzzz[j] - fr * g2_y_xzzz[j];

                g_xy_yyyy[j] = g1_y_xyyyy[j] - fr * g2_y_yyyy[j];

                g_xy_yyyz[j] = g1_y_xyyyz[j] - fr * g2_y_yyyz[j];

                g_xy_yyzz[j] = g1_y_xyyzz[j] - fr * g2_y_yyzz[j];

                g_xy_yzzz[j] = g1_y_xyzzz[j] - fr * g2_y_yzzz[j];

                g_xy_zzzz[j] = g1_y_xzzzz[j] - fr * g2_y_zzzz[j];

                g_xz_xxxx[j] = g1_z_xxxxx[j] - fr * g2_z_xxxx[j];

                g_xz_xxxy[j] = g1_z_xxxxy[j] - fr * g2_z_xxxy[j];

                g_xz_xxxz[j] = g1_z_xxxxz[j] - fr * g2_z_xxxz[j];

                g_xz_xxyy[j] = g1_z_xxxyy[j] - fr * g2_z_xxyy[j];

                g_xz_xxyz[j] = g1_z_xxxyz[j] - fr * g2_z_xxyz[j];

                g_xz_xxzz[j] = g1_z_xxxzz[j] - fr * g2_z_xxzz[j];

                g_xz_xyyy[j] = g1_z_xxyyy[j] - fr * g2_z_xyyy[j];

                g_xz_xyyz[j] = g1_z_xxyyz[j] - fr * g2_z_xyyz[j];

                g_xz_xyzz[j] = g1_z_xxyzz[j] - fr * g2_z_xyzz[j];

                g_xz_xzzz[j] = g1_z_xxzzz[j] - fr * g2_z_xzzz[j];

                g_xz_yyyy[j] = g1_z_xyyyy[j] - fr * g2_z_yyyy[j];

                g_xz_yyyz[j] = g1_z_xyyyz[j] - fr * g2_z_yyyz[j];

                g_xz_yyzz[j] = g1_z_xyyzz[j] - fr * g2_z_yyzz[j];

                g_xz_yzzz[j] = g1_z_xyzzz[j] - fr * g2_z_yzzz[j];

                g_xz_zzzz[j] = g1_z_xzzzz[j] - fr * g2_z_zzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yy_xxxx[j] = g1_y_xxxxy[j] - fr * g2_y_xxxx[j];

                g_yy_xxxy[j] = g1_y_xxxyy[j] - fr * g2_y_xxxy[j];

                g_yy_xxxz[j] = g1_y_xxxyz[j] - fr * g2_y_xxxz[j];

                g_yy_xxyy[j] = g1_y_xxyyy[j] - fr * g2_y_xxyy[j];

                g_yy_xxyz[j] = g1_y_xxyyz[j] - fr * g2_y_xxyz[j];

                g_yy_xxzz[j] = g1_y_xxyzz[j] - fr * g2_y_xxzz[j];

                g_yy_xyyy[j] = g1_y_xyyyy[j] - fr * g2_y_xyyy[j];

                g_yy_xyyz[j] = g1_y_xyyyz[j] - fr * g2_y_xyyz[j];

                g_yy_xyzz[j] = g1_y_xyyzz[j] - fr * g2_y_xyzz[j];

                g_yy_xzzz[j] = g1_y_xyzzz[j] - fr * g2_y_xzzz[j];

                g_yy_yyyy[j] = g1_y_yyyyy[j] - fr * g2_y_yyyy[j];

                g_yy_yyyz[j] = g1_y_yyyyz[j] - fr * g2_y_yyyz[j];

                g_yy_yyzz[j] = g1_y_yyyzz[j] - fr * g2_y_yyzz[j];

                g_yy_yzzz[j] = g1_y_yyzzz[j] - fr * g2_y_yzzz[j];

                g_yy_zzzz[j] = g1_y_yzzzz[j] - fr * g2_y_zzzz[j];

                g_yz_xxxx[j] = g1_z_xxxxy[j] - fr * g2_z_xxxx[j];

                g_yz_xxxy[j] = g1_z_xxxyy[j] - fr * g2_z_xxxy[j];

                g_yz_xxxz[j] = g1_z_xxxyz[j] - fr * g2_z_xxxz[j];

                g_yz_xxyy[j] = g1_z_xxyyy[j] - fr * g2_z_xxyy[j];

                g_yz_xxyz[j] = g1_z_xxyyz[j] - fr * g2_z_xxyz[j];

                g_yz_xxzz[j] = g1_z_xxyzz[j] - fr * g2_z_xxzz[j];

                g_yz_xyyy[j] = g1_z_xyyyy[j] - fr * g2_z_xyyy[j];

                g_yz_xyyz[j] = g1_z_xyyyz[j] - fr * g2_z_xyyz[j];

                g_yz_xyzz[j] = g1_z_xyyzz[j] - fr * g2_z_xyzz[j];

                g_yz_xzzz[j] = g1_z_xyzzz[j] - fr * g2_z_xzzz[j];

                g_yz_yyyy[j] = g1_z_yyyyy[j] - fr * g2_z_yyyy[j];

                g_yz_yyyz[j] = g1_z_yyyyz[j] - fr * g2_z_yyyz[j];

                g_yz_yyzz[j] = g1_z_yyyzz[j] - fr * g2_z_yyzz[j];

                g_yz_yzzz[j] = g1_z_yyzzz[j] - fr * g2_z_yzzz[j];

                g_yz_zzzz[j] = g1_z_yzzzz[j] - fr * g2_z_zzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zz_xxxx[j] = g1_z_xxxxz[j] - fr * g2_z_xxxx[j];

                g_zz_xxxy[j] = g1_z_xxxyz[j] - fr * g2_z_xxxy[j];

                g_zz_xxxz[j] = g1_z_xxxzz[j] - fr * g2_z_xxxz[j];

                g_zz_xxyy[j] = g1_z_xxyyz[j] - fr * g2_z_xxyy[j];

                g_zz_xxyz[j] = g1_z_xxyzz[j] - fr * g2_z_xxyz[j];

                g_zz_xxzz[j] = g1_z_xxzzz[j] - fr * g2_z_xxzz[j];

                g_zz_xyyy[j] = g1_z_xyyyz[j] - fr * g2_z_xyyy[j];

                g_zz_xyyz[j] = g1_z_xyyzz[j] - fr * g2_z_xyyz[j];

                g_zz_xyzz[j] = g1_z_xyzzz[j] - fr * g2_z_xyzz[j];

                g_zz_xzzz[j] = g1_z_xzzzz[j] - fr * g2_z_xzzz[j];

                g_zz_yyyy[j] = g1_z_yyyyz[j] - fr * g2_z_yyyy[j];

                g_zz_yyyz[j] = g1_z_yyyzz[j] - fr * g2_z_yyyz[j];

                g_zz_yyzz[j] = g1_z_yyzzz[j] - fr * g2_z_yyzz[j];

                g_zz_yzzz[j] = g1_z_yzzzz[j] - fr * g2_z_yzzz[j];

                g_zz_zzzz[j] = g1_z_zzzzz[j] - fr * g2_z_zzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXDH(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 2, 5})) return;

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
                                             {braAngularMomentum, 2, 5});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 6});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 5});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|PH)^(m) integrals

            auto g2_x_xxxxx = contrBuffer.data(g2off + 63 * i);

            auto g2_x_xxxxy = contrBuffer.data(g2off + 63 * i + 1);

            auto g2_x_xxxxz = contrBuffer.data(g2off + 63 * i + 2);

            auto g2_x_xxxyy = contrBuffer.data(g2off + 63 * i + 3);

            auto g2_x_xxxyz = contrBuffer.data(g2off + 63 * i + 4);

            auto g2_x_xxxzz = contrBuffer.data(g2off + 63 * i + 5);

            auto g2_x_xxyyy = contrBuffer.data(g2off + 63 * i + 6);

            auto g2_x_xxyyz = contrBuffer.data(g2off + 63 * i + 7);

            auto g2_x_xxyzz = contrBuffer.data(g2off + 63 * i + 8);

            auto g2_x_xxzzz = contrBuffer.data(g2off + 63 * i + 9);

            auto g2_x_xyyyy = contrBuffer.data(g2off + 63 * i + 10);

            auto g2_x_xyyyz = contrBuffer.data(g2off + 63 * i + 11);

            auto g2_x_xyyzz = contrBuffer.data(g2off + 63 * i + 12);

            auto g2_x_xyzzz = contrBuffer.data(g2off + 63 * i + 13);

            auto g2_x_xzzzz = contrBuffer.data(g2off + 63 * i + 14);

            auto g2_x_yyyyy = contrBuffer.data(g2off + 63 * i + 15);

            auto g2_x_yyyyz = contrBuffer.data(g2off + 63 * i + 16);

            auto g2_x_yyyzz = contrBuffer.data(g2off + 63 * i + 17);

            auto g2_x_yyzzz = contrBuffer.data(g2off + 63 * i + 18);

            auto g2_x_yzzzz = contrBuffer.data(g2off + 63 * i + 19);

            auto g2_x_zzzzz = contrBuffer.data(g2off + 63 * i + 20);

            auto g2_y_xxxxx = contrBuffer.data(g2off + 63 * i + 21);

            auto g2_y_xxxxy = contrBuffer.data(g2off + 63 * i + 22);

            auto g2_y_xxxxz = contrBuffer.data(g2off + 63 * i + 23);

            auto g2_y_xxxyy = contrBuffer.data(g2off + 63 * i + 24);

            auto g2_y_xxxyz = contrBuffer.data(g2off + 63 * i + 25);

            auto g2_y_xxxzz = contrBuffer.data(g2off + 63 * i + 26);

            auto g2_y_xxyyy = contrBuffer.data(g2off + 63 * i + 27);

            auto g2_y_xxyyz = contrBuffer.data(g2off + 63 * i + 28);

            auto g2_y_xxyzz = contrBuffer.data(g2off + 63 * i + 29);

            auto g2_y_xxzzz = contrBuffer.data(g2off + 63 * i + 30);

            auto g2_y_xyyyy = contrBuffer.data(g2off + 63 * i + 31);

            auto g2_y_xyyyz = contrBuffer.data(g2off + 63 * i + 32);

            auto g2_y_xyyzz = contrBuffer.data(g2off + 63 * i + 33);

            auto g2_y_xyzzz = contrBuffer.data(g2off + 63 * i + 34);

            auto g2_y_xzzzz = contrBuffer.data(g2off + 63 * i + 35);

            auto g2_y_yyyyy = contrBuffer.data(g2off + 63 * i + 36);

            auto g2_y_yyyyz = contrBuffer.data(g2off + 63 * i + 37);

            auto g2_y_yyyzz = contrBuffer.data(g2off + 63 * i + 38);

            auto g2_y_yyzzz = contrBuffer.data(g2off + 63 * i + 39);

            auto g2_y_yzzzz = contrBuffer.data(g2off + 63 * i + 40);

            auto g2_y_zzzzz = contrBuffer.data(g2off + 63 * i + 41);

            auto g2_z_xxxxx = contrBuffer.data(g2off + 63 * i + 42);

            auto g2_z_xxxxy = contrBuffer.data(g2off + 63 * i + 43);

            auto g2_z_xxxxz = contrBuffer.data(g2off + 63 * i + 44);

            auto g2_z_xxxyy = contrBuffer.data(g2off + 63 * i + 45);

            auto g2_z_xxxyz = contrBuffer.data(g2off + 63 * i + 46);

            auto g2_z_xxxzz = contrBuffer.data(g2off + 63 * i + 47);

            auto g2_z_xxyyy = contrBuffer.data(g2off + 63 * i + 48);

            auto g2_z_xxyyz = contrBuffer.data(g2off + 63 * i + 49);

            auto g2_z_xxyzz = contrBuffer.data(g2off + 63 * i + 50);

            auto g2_z_xxzzz = contrBuffer.data(g2off + 63 * i + 51);

            auto g2_z_xyyyy = contrBuffer.data(g2off + 63 * i + 52);

            auto g2_z_xyyyz = contrBuffer.data(g2off + 63 * i + 53);

            auto g2_z_xyyzz = contrBuffer.data(g2off + 63 * i + 54);

            auto g2_z_xyzzz = contrBuffer.data(g2off + 63 * i + 55);

            auto g2_z_xzzzz = contrBuffer.data(g2off + 63 * i + 56);

            auto g2_z_yyyyy = contrBuffer.data(g2off + 63 * i + 57);

            auto g2_z_yyyyz = contrBuffer.data(g2off + 63 * i + 58);

            auto g2_z_yyyzz = contrBuffer.data(g2off + 63 * i + 59);

            auto g2_z_yyzzz = contrBuffer.data(g2off + 63 * i + 60);

            auto g2_z_yzzzz = contrBuffer.data(g2off + 63 * i + 61);

            auto g2_z_zzzzz = contrBuffer.data(g2off + 63 * i + 62);

            // set up pointers to (X|g(r,r')|PI)^(m) integrals

            auto g1_x_xxxxxx = contrBuffer.data(g1off + 84 * i);

            auto g1_x_xxxxxy = contrBuffer.data(g1off + 84 * i + 1);

            auto g1_x_xxxxxz = contrBuffer.data(g1off + 84 * i + 2);

            auto g1_x_xxxxyy = contrBuffer.data(g1off + 84 * i + 3);

            auto g1_x_xxxxyz = contrBuffer.data(g1off + 84 * i + 4);

            auto g1_x_xxxxzz = contrBuffer.data(g1off + 84 * i + 5);

            auto g1_x_xxxyyy = contrBuffer.data(g1off + 84 * i + 6);

            auto g1_x_xxxyyz = contrBuffer.data(g1off + 84 * i + 7);

            auto g1_x_xxxyzz = contrBuffer.data(g1off + 84 * i + 8);

            auto g1_x_xxxzzz = contrBuffer.data(g1off + 84 * i + 9);

            auto g1_x_xxyyyy = contrBuffer.data(g1off + 84 * i + 10);

            auto g1_x_xxyyyz = contrBuffer.data(g1off + 84 * i + 11);

            auto g1_x_xxyyzz = contrBuffer.data(g1off + 84 * i + 12);

            auto g1_x_xxyzzz = contrBuffer.data(g1off + 84 * i + 13);

            auto g1_x_xxzzzz = contrBuffer.data(g1off + 84 * i + 14);

            auto g1_x_xyyyyy = contrBuffer.data(g1off + 84 * i + 15);

            auto g1_x_xyyyyz = contrBuffer.data(g1off + 84 * i + 16);

            auto g1_x_xyyyzz = contrBuffer.data(g1off + 84 * i + 17);

            auto g1_x_xyyzzz = contrBuffer.data(g1off + 84 * i + 18);

            auto g1_x_xyzzzz = contrBuffer.data(g1off + 84 * i + 19);

            auto g1_x_xzzzzz = contrBuffer.data(g1off + 84 * i + 20);

            auto g1_y_xxxxxx = contrBuffer.data(g1off + 84 * i + 28);

            auto g1_y_xxxxxy = contrBuffer.data(g1off + 84 * i + 29);

            auto g1_y_xxxxxz = contrBuffer.data(g1off + 84 * i + 30);

            auto g1_y_xxxxyy = contrBuffer.data(g1off + 84 * i + 31);

            auto g1_y_xxxxyz = contrBuffer.data(g1off + 84 * i + 32);

            auto g1_y_xxxxzz = contrBuffer.data(g1off + 84 * i + 33);

            auto g1_y_xxxyyy = contrBuffer.data(g1off + 84 * i + 34);

            auto g1_y_xxxyyz = contrBuffer.data(g1off + 84 * i + 35);

            auto g1_y_xxxyzz = contrBuffer.data(g1off + 84 * i + 36);

            auto g1_y_xxxzzz = contrBuffer.data(g1off + 84 * i + 37);

            auto g1_y_xxyyyy = contrBuffer.data(g1off + 84 * i + 38);

            auto g1_y_xxyyyz = contrBuffer.data(g1off + 84 * i + 39);

            auto g1_y_xxyyzz = contrBuffer.data(g1off + 84 * i + 40);

            auto g1_y_xxyzzz = contrBuffer.data(g1off + 84 * i + 41);

            auto g1_y_xxzzzz = contrBuffer.data(g1off + 84 * i + 42);

            auto g1_y_xyyyyy = contrBuffer.data(g1off + 84 * i + 43);

            auto g1_y_xyyyyz = contrBuffer.data(g1off + 84 * i + 44);

            auto g1_y_xyyyzz = contrBuffer.data(g1off + 84 * i + 45);

            auto g1_y_xyyzzz = contrBuffer.data(g1off + 84 * i + 46);

            auto g1_y_xyzzzz = contrBuffer.data(g1off + 84 * i + 47);

            auto g1_y_xzzzzz = contrBuffer.data(g1off + 84 * i + 48);

            auto g1_y_yyyyyy = contrBuffer.data(g1off + 84 * i + 49);

            auto g1_y_yyyyyz = contrBuffer.data(g1off + 84 * i + 50);

            auto g1_y_yyyyzz = contrBuffer.data(g1off + 84 * i + 51);

            auto g1_y_yyyzzz = contrBuffer.data(g1off + 84 * i + 52);

            auto g1_y_yyzzzz = contrBuffer.data(g1off + 84 * i + 53);

            auto g1_y_yzzzzz = contrBuffer.data(g1off + 84 * i + 54);

            auto g1_z_xxxxxx = contrBuffer.data(g1off + 84 * i + 56);

            auto g1_z_xxxxxy = contrBuffer.data(g1off + 84 * i + 57);

            auto g1_z_xxxxxz = contrBuffer.data(g1off + 84 * i + 58);

            auto g1_z_xxxxyy = contrBuffer.data(g1off + 84 * i + 59);

            auto g1_z_xxxxyz = contrBuffer.data(g1off + 84 * i + 60);

            auto g1_z_xxxxzz = contrBuffer.data(g1off + 84 * i + 61);

            auto g1_z_xxxyyy = contrBuffer.data(g1off + 84 * i + 62);

            auto g1_z_xxxyyz = contrBuffer.data(g1off + 84 * i + 63);

            auto g1_z_xxxyzz = contrBuffer.data(g1off + 84 * i + 64);

            auto g1_z_xxxzzz = contrBuffer.data(g1off + 84 * i + 65);

            auto g1_z_xxyyyy = contrBuffer.data(g1off + 84 * i + 66);

            auto g1_z_xxyyyz = contrBuffer.data(g1off + 84 * i + 67);

            auto g1_z_xxyyzz = contrBuffer.data(g1off + 84 * i + 68);

            auto g1_z_xxyzzz = contrBuffer.data(g1off + 84 * i + 69);

            auto g1_z_xxzzzz = contrBuffer.data(g1off + 84 * i + 70);

            auto g1_z_xyyyyy = contrBuffer.data(g1off + 84 * i + 71);

            auto g1_z_xyyyyz = contrBuffer.data(g1off + 84 * i + 72);

            auto g1_z_xyyyzz = contrBuffer.data(g1off + 84 * i + 73);

            auto g1_z_xyyzzz = contrBuffer.data(g1off + 84 * i + 74);

            auto g1_z_xyzzzz = contrBuffer.data(g1off + 84 * i + 75);

            auto g1_z_xzzzzz = contrBuffer.data(g1off + 84 * i + 76);

            auto g1_z_yyyyyy = contrBuffer.data(g1off + 84 * i + 77);

            auto g1_z_yyyyyz = contrBuffer.data(g1off + 84 * i + 78);

            auto g1_z_yyyyzz = contrBuffer.data(g1off + 84 * i + 79);

            auto g1_z_yyyzzz = contrBuffer.data(g1off + 84 * i + 80);

            auto g1_z_yyzzzz = contrBuffer.data(g1off + 84 * i + 81);

            auto g1_z_yzzzzz = contrBuffer.data(g1off + 84 * i + 82);

            auto g1_z_zzzzzz = contrBuffer.data(g1off + 84 * i + 83);

            // set up pointers to (X|g(r,r')|DH)^(m) integrals

            auto g_xx_xxxxx = contrBuffer.data(goff + 126 * i);

            auto g_xx_xxxxy = contrBuffer.data(goff + 126 * i + 1);

            auto g_xx_xxxxz = contrBuffer.data(goff + 126 * i + 2);

            auto g_xx_xxxyy = contrBuffer.data(goff + 126 * i + 3);

            auto g_xx_xxxyz = contrBuffer.data(goff + 126 * i + 4);

            auto g_xx_xxxzz = contrBuffer.data(goff + 126 * i + 5);

            auto g_xx_xxyyy = contrBuffer.data(goff + 126 * i + 6);

            auto g_xx_xxyyz = contrBuffer.data(goff + 126 * i + 7);

            auto g_xx_xxyzz = contrBuffer.data(goff + 126 * i + 8);

            auto g_xx_xxzzz = contrBuffer.data(goff + 126 * i + 9);

            auto g_xx_xyyyy = contrBuffer.data(goff + 126 * i + 10);

            auto g_xx_xyyyz = contrBuffer.data(goff + 126 * i + 11);

            auto g_xx_xyyzz = contrBuffer.data(goff + 126 * i + 12);

            auto g_xx_xyzzz = contrBuffer.data(goff + 126 * i + 13);

            auto g_xx_xzzzz = contrBuffer.data(goff + 126 * i + 14);

            auto g_xx_yyyyy = contrBuffer.data(goff + 126 * i + 15);

            auto g_xx_yyyyz = contrBuffer.data(goff + 126 * i + 16);

            auto g_xx_yyyzz = contrBuffer.data(goff + 126 * i + 17);

            auto g_xx_yyzzz = contrBuffer.data(goff + 126 * i + 18);

            auto g_xx_yzzzz = contrBuffer.data(goff + 126 * i + 19);

            auto g_xx_zzzzz = contrBuffer.data(goff + 126 * i + 20);

            auto g_xy_xxxxx = contrBuffer.data(goff + 126 * i + 21);

            auto g_xy_xxxxy = contrBuffer.data(goff + 126 * i + 22);

            auto g_xy_xxxxz = contrBuffer.data(goff + 126 * i + 23);

            auto g_xy_xxxyy = contrBuffer.data(goff + 126 * i + 24);

            auto g_xy_xxxyz = contrBuffer.data(goff + 126 * i + 25);

            auto g_xy_xxxzz = contrBuffer.data(goff + 126 * i + 26);

            auto g_xy_xxyyy = contrBuffer.data(goff + 126 * i + 27);

            auto g_xy_xxyyz = contrBuffer.data(goff + 126 * i + 28);

            auto g_xy_xxyzz = contrBuffer.data(goff + 126 * i + 29);

            auto g_xy_xxzzz = contrBuffer.data(goff + 126 * i + 30);

            auto g_xy_xyyyy = contrBuffer.data(goff + 126 * i + 31);

            auto g_xy_xyyyz = contrBuffer.data(goff + 126 * i + 32);

            auto g_xy_xyyzz = contrBuffer.data(goff + 126 * i + 33);

            auto g_xy_xyzzz = contrBuffer.data(goff + 126 * i + 34);

            auto g_xy_xzzzz = contrBuffer.data(goff + 126 * i + 35);

            auto g_xy_yyyyy = contrBuffer.data(goff + 126 * i + 36);

            auto g_xy_yyyyz = contrBuffer.data(goff + 126 * i + 37);

            auto g_xy_yyyzz = contrBuffer.data(goff + 126 * i + 38);

            auto g_xy_yyzzz = contrBuffer.data(goff + 126 * i + 39);

            auto g_xy_yzzzz = contrBuffer.data(goff + 126 * i + 40);

            auto g_xy_zzzzz = contrBuffer.data(goff + 126 * i + 41);

            auto g_xz_xxxxx = contrBuffer.data(goff + 126 * i + 42);

            auto g_xz_xxxxy = contrBuffer.data(goff + 126 * i + 43);

            auto g_xz_xxxxz = contrBuffer.data(goff + 126 * i + 44);

            auto g_xz_xxxyy = contrBuffer.data(goff + 126 * i + 45);

            auto g_xz_xxxyz = contrBuffer.data(goff + 126 * i + 46);

            auto g_xz_xxxzz = contrBuffer.data(goff + 126 * i + 47);

            auto g_xz_xxyyy = contrBuffer.data(goff + 126 * i + 48);

            auto g_xz_xxyyz = contrBuffer.data(goff + 126 * i + 49);

            auto g_xz_xxyzz = contrBuffer.data(goff + 126 * i + 50);

            auto g_xz_xxzzz = contrBuffer.data(goff + 126 * i + 51);

            auto g_xz_xyyyy = contrBuffer.data(goff + 126 * i + 52);

            auto g_xz_xyyyz = contrBuffer.data(goff + 126 * i + 53);

            auto g_xz_xyyzz = contrBuffer.data(goff + 126 * i + 54);

            auto g_xz_xyzzz = contrBuffer.data(goff + 126 * i + 55);

            auto g_xz_xzzzz = contrBuffer.data(goff + 126 * i + 56);

            auto g_xz_yyyyy = contrBuffer.data(goff + 126 * i + 57);

            auto g_xz_yyyyz = contrBuffer.data(goff + 126 * i + 58);

            auto g_xz_yyyzz = contrBuffer.data(goff + 126 * i + 59);

            auto g_xz_yyzzz = contrBuffer.data(goff + 126 * i + 60);

            auto g_xz_yzzzz = contrBuffer.data(goff + 126 * i + 61);

            auto g_xz_zzzzz = contrBuffer.data(goff + 126 * i + 62);

            auto g_yy_xxxxx = contrBuffer.data(goff + 126 * i + 63);

            auto g_yy_xxxxy = contrBuffer.data(goff + 126 * i + 64);

            auto g_yy_xxxxz = contrBuffer.data(goff + 126 * i + 65);

            auto g_yy_xxxyy = contrBuffer.data(goff + 126 * i + 66);

            auto g_yy_xxxyz = contrBuffer.data(goff + 126 * i + 67);

            auto g_yy_xxxzz = contrBuffer.data(goff + 126 * i + 68);

            auto g_yy_xxyyy = contrBuffer.data(goff + 126 * i + 69);

            auto g_yy_xxyyz = contrBuffer.data(goff + 126 * i + 70);

            auto g_yy_xxyzz = contrBuffer.data(goff + 126 * i + 71);

            auto g_yy_xxzzz = contrBuffer.data(goff + 126 * i + 72);

            auto g_yy_xyyyy = contrBuffer.data(goff + 126 * i + 73);

            auto g_yy_xyyyz = contrBuffer.data(goff + 126 * i + 74);

            auto g_yy_xyyzz = contrBuffer.data(goff + 126 * i + 75);

            auto g_yy_xyzzz = contrBuffer.data(goff + 126 * i + 76);

            auto g_yy_xzzzz = contrBuffer.data(goff + 126 * i + 77);

            auto g_yy_yyyyy = contrBuffer.data(goff + 126 * i + 78);

            auto g_yy_yyyyz = contrBuffer.data(goff + 126 * i + 79);

            auto g_yy_yyyzz = contrBuffer.data(goff + 126 * i + 80);

            auto g_yy_yyzzz = contrBuffer.data(goff + 126 * i + 81);

            auto g_yy_yzzzz = contrBuffer.data(goff + 126 * i + 82);

            auto g_yy_zzzzz = contrBuffer.data(goff + 126 * i + 83);

            auto g_yz_xxxxx = contrBuffer.data(goff + 126 * i + 84);

            auto g_yz_xxxxy = contrBuffer.data(goff + 126 * i + 85);

            auto g_yz_xxxxz = contrBuffer.data(goff + 126 * i + 86);

            auto g_yz_xxxyy = contrBuffer.data(goff + 126 * i + 87);

            auto g_yz_xxxyz = contrBuffer.data(goff + 126 * i + 88);

            auto g_yz_xxxzz = contrBuffer.data(goff + 126 * i + 89);

            auto g_yz_xxyyy = contrBuffer.data(goff + 126 * i + 90);

            auto g_yz_xxyyz = contrBuffer.data(goff + 126 * i + 91);

            auto g_yz_xxyzz = contrBuffer.data(goff + 126 * i + 92);

            auto g_yz_xxzzz = contrBuffer.data(goff + 126 * i + 93);

            auto g_yz_xyyyy = contrBuffer.data(goff + 126 * i + 94);

            auto g_yz_xyyyz = contrBuffer.data(goff + 126 * i + 95);

            auto g_yz_xyyzz = contrBuffer.data(goff + 126 * i + 96);

            auto g_yz_xyzzz = contrBuffer.data(goff + 126 * i + 97);

            auto g_yz_xzzzz = contrBuffer.data(goff + 126 * i + 98);

            auto g_yz_yyyyy = contrBuffer.data(goff + 126 * i + 99);

            auto g_yz_yyyyz = contrBuffer.data(goff + 126 * i + 100);

            auto g_yz_yyyzz = contrBuffer.data(goff + 126 * i + 101);

            auto g_yz_yyzzz = contrBuffer.data(goff + 126 * i + 102);

            auto g_yz_yzzzz = contrBuffer.data(goff + 126 * i + 103);

            auto g_yz_zzzzz = contrBuffer.data(goff + 126 * i + 104);

            auto g_zz_xxxxx = contrBuffer.data(goff + 126 * i + 105);

            auto g_zz_xxxxy = contrBuffer.data(goff + 126 * i + 106);

            auto g_zz_xxxxz = contrBuffer.data(goff + 126 * i + 107);

            auto g_zz_xxxyy = contrBuffer.data(goff + 126 * i + 108);

            auto g_zz_xxxyz = contrBuffer.data(goff + 126 * i + 109);

            auto g_zz_xxxzz = contrBuffer.data(goff + 126 * i + 110);

            auto g_zz_xxyyy = contrBuffer.data(goff + 126 * i + 111);

            auto g_zz_xxyyz = contrBuffer.data(goff + 126 * i + 112);

            auto g_zz_xxyzz = contrBuffer.data(goff + 126 * i + 113);

            auto g_zz_xxzzz = contrBuffer.data(goff + 126 * i + 114);

            auto g_zz_xyyyy = contrBuffer.data(goff + 126 * i + 115);

            auto g_zz_xyyyz = contrBuffer.data(goff + 126 * i + 116);

            auto g_zz_xyyzz = contrBuffer.data(goff + 126 * i + 117);

            auto g_zz_xyzzz = contrBuffer.data(goff + 126 * i + 118);

            auto g_zz_xzzzz = contrBuffer.data(goff + 126 * i + 119);

            auto g_zz_yyyyy = contrBuffer.data(goff + 126 * i + 120);

            auto g_zz_yyyyz = contrBuffer.data(goff + 126 * i + 121);

            auto g_zz_yyyzz = contrBuffer.data(goff + 126 * i + 122);

            auto g_zz_yyzzz = contrBuffer.data(goff + 126 * i + 123);

            auto g_zz_yzzzz = contrBuffer.data(goff + 126 * i + 124);

            auto g_zz_zzzzz = contrBuffer.data(goff + 126 * i + 125);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxxx, g2_x_xxxxy,\
                                     g2_x_xxxxz, g2_x_xxxyy, g2_x_xxxyz, g2_x_xxxzz,\
                                     g2_x_xxyyy, g2_x_xxyyz, g2_x_xxyzz, g2_x_xxzzz,\
                                     g2_x_xyyyy, g2_x_xyyyz, g2_x_xyyzz, g2_x_xyzzz,\
                                     g2_x_xzzzz, g2_x_yyyyy, g2_x_yyyyz, g2_x_yyyzz,\
                                     g2_x_yyzzz, g2_x_yzzzz, g2_x_zzzzz, g2_y_xxxxx,\
                                     g2_y_xxxxy, g2_y_xxxxz, g2_y_xxxyy, g2_y_xxxyz,\
                                     g2_y_xxxzz, g2_y_xxyyy, g2_y_xxyyz, g2_y_xxyzz,\
                                     g2_y_xxzzz, g2_y_xyyyy, g2_y_xyyyz, g2_y_xyyzz,\
                                     g2_y_xyzzz, g2_y_xzzzz, g2_y_yyyyy, g2_y_yyyyz,\
                                     g2_y_yyyzz, g2_y_yyzzz, g2_y_yzzzz, g2_y_zzzzz,\
                                     g2_z_xxxxx, g2_z_xxxxy, g2_z_xxxxz, g2_z_xxxyy,\
                                     g2_z_xxxyz, g2_z_xxxzz, g2_z_xxyyy, g2_z_xxyyz,\
                                     g2_z_xxyzz, g2_z_xxzzz, g2_z_xyyyy, g2_z_xyyyz,\
                                     g2_z_xyyzz, g2_z_xyzzz, g2_z_xzzzz, g2_z_yyyyy,\
                                     g2_z_yyyyz, g2_z_yyyzz, g2_z_yyzzz, g2_z_yzzzz,\
                                     g2_z_zzzzz, g1_x_xxxxxx, g1_x_xxxxxy, g1_x_xxxxxz,\
                                     g1_x_xxxxyy, g1_x_xxxxyz, g1_x_xxxxzz, g1_x_xxxyyy,\
                                     g1_x_xxxyyz, g1_x_xxxyzz, g1_x_xxxzzz, g1_x_xxyyyy,\
                                     g1_x_xxyyyz, g1_x_xxyyzz, g1_x_xxyzzz, g1_x_xxzzzz,\
                                     g1_x_xyyyyy, g1_x_xyyyyz, g1_x_xyyyzz, g1_x_xyyzzz,\
                                     g1_x_xyzzzz, g1_x_xzzzzz, g1_y_xxxxxx, g1_y_xxxxxy,\
                                     g1_y_xxxxxz, g1_y_xxxxyy, g1_y_xxxxyz, g1_y_xxxxzz,\
                                     g1_y_xxxyyy, g1_y_xxxyyz, g1_y_xxxyzz, g1_y_xxxzzz,\
                                     g1_y_xxyyyy, g1_y_xxyyyz, g1_y_xxyyzz, g1_y_xxyzzz,\
                                     g1_y_xxzzzz, g1_y_xyyyyy, g1_y_xyyyyz, g1_y_xyyyzz,\
                                     g1_y_xyyzzz, g1_y_xyzzzz, g1_y_xzzzzz, g1_y_yyyyyy,\
                                     g1_y_yyyyyz, g1_y_yyyyzz, g1_y_yyyzzz, g1_y_yyzzzz,\
                                     g1_y_yzzzzz, g1_z_xxxxxx, g1_z_xxxxxy, g1_z_xxxxxz,\
                                     g1_z_xxxxyy, g1_z_xxxxyz, g1_z_xxxxzz, g1_z_xxxyyy,\
                                     g1_z_xxxyyz, g1_z_xxxyzz, g1_z_xxxzzz, g1_z_xxyyyy,\
                                     g1_z_xxyyyz, g1_z_xxyyzz, g1_z_xxyzzz, g1_z_xxzzzz,\
                                     g1_z_xyyyyy, g1_z_xyyyyz, g1_z_xyyyzz, g1_z_xyyzzz,\
                                     g1_z_xyzzzz, g1_z_xzzzzz, g1_z_yyyyyy, g1_z_yyyyyz,\
                                     g1_z_yyyyzz, g1_z_yyyzzz, g1_z_yyzzzz, g1_z_yzzzzz,\
                                     g1_z_zzzzzz, g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz,\
                                     g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxzz, g_xx_xxyyy,\
                                     g_xx_xxyyz, g_xx_xxyzz, g_xx_xxzzz, g_xx_xyyyy,\
                                     g_xx_xyyyz, g_xx_xyyzz, g_xx_xyzzz, g_xx_xzzzz,\
                                     g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz, g_xx_yyzzz,\
                                     g_xx_yzzzz, g_xx_zzzzz, g_xy_xxxxx, g_xy_xxxxy,\
                                     g_xy_xxxxz, g_xy_xxxyy, g_xy_xxxyz, g_xy_xxxzz,\
                                     g_xy_xxyyy, g_xy_xxyyz, g_xy_xxyzz, g_xy_xxzzz,\
                                     g_xy_xyyyy, g_xy_xyyyz, g_xy_xyyzz, g_xy_xyzzz,\
                                     g_xy_xzzzz, g_xy_yyyyy, g_xy_yyyyz, g_xy_yyyzz,\
                                     g_xy_yyzzz, g_xy_yzzzz, g_xy_zzzzz, g_xz_xxxxx,\
                                     g_xz_xxxxy, g_xz_xxxxz, g_xz_xxxyy, g_xz_xxxyz,\
                                     g_xz_xxxzz, g_xz_xxyyy, g_xz_xxyyz, g_xz_xxyzz,\
                                     g_xz_xxzzz, g_xz_xyyyy, g_xz_xyyyz, g_xz_xyyzz,\
                                     g_xz_xyzzz, g_xz_xzzzz, g_xz_yyyyy, g_xz_yyyyz,\
                                     g_xz_yyyzz, g_xz_yyzzz, g_xz_yzzzz, g_xz_zzzzz,\
                                     g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz, g_yy_xxxyy,\
                                     g_yy_xxxyz, g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyz,\
                                     g_yy_xxyzz, g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyz,\
                                     g_yy_xyyzz, g_yy_xyzzz, g_yy_xzzzz, g_yy_yyyyy,\
                                     g_yy_yyyyz, g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz,\
                                     g_yy_zzzzz, g_yz_xxxxx, g_yz_xxxxy, g_yz_xxxxz,\
                                     g_yz_xxxyy, g_yz_xxxyz, g_yz_xxxzz, g_yz_xxyyy,\
                                     g_yz_xxyyz, g_yz_xxyzz, g_yz_xxzzz, g_yz_xyyyy,\
                                     g_yz_xyyyz, g_yz_xyyzz, g_yz_xyzzz, g_yz_xzzzz,\
                                     g_yz_yyyyy, g_yz_yyyyz, g_yz_yyyzz, g_yz_yyzzz,\
                                     g_yz_yzzzz, g_yz_zzzzz, g_zz_xxxxx, g_zz_xxxxy,\
                                     g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxzz,\
                                     g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyzz, g_zz_xxzzz,\
                                     g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyzz, g_zz_xyzzz,\
                                     g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz,\
                                     g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xx_xxxxx[j] = g1_x_xxxxxx[j] - fr * g2_x_xxxxx[j];

                g_xx_xxxxy[j] = g1_x_xxxxxy[j] - fr * g2_x_xxxxy[j];

                g_xx_xxxxz[j] = g1_x_xxxxxz[j] - fr * g2_x_xxxxz[j];

                g_xx_xxxyy[j] = g1_x_xxxxyy[j] - fr * g2_x_xxxyy[j];

                g_xx_xxxyz[j] = g1_x_xxxxyz[j] - fr * g2_x_xxxyz[j];

                g_xx_xxxzz[j] = g1_x_xxxxzz[j] - fr * g2_x_xxxzz[j];

                g_xx_xxyyy[j] = g1_x_xxxyyy[j] - fr * g2_x_xxyyy[j];

                g_xx_xxyyz[j] = g1_x_xxxyyz[j] - fr * g2_x_xxyyz[j];

                g_xx_xxyzz[j] = g1_x_xxxyzz[j] - fr * g2_x_xxyzz[j];

                g_xx_xxzzz[j] = g1_x_xxxzzz[j] - fr * g2_x_xxzzz[j];

                g_xx_xyyyy[j] = g1_x_xxyyyy[j] - fr * g2_x_xyyyy[j];

                g_xx_xyyyz[j] = g1_x_xxyyyz[j] - fr * g2_x_xyyyz[j];

                g_xx_xyyzz[j] = g1_x_xxyyzz[j] - fr * g2_x_xyyzz[j];

                g_xx_xyzzz[j] = g1_x_xxyzzz[j] - fr * g2_x_xyzzz[j];

                g_xx_xzzzz[j] = g1_x_xxzzzz[j] - fr * g2_x_xzzzz[j];

                g_xx_yyyyy[j] = g1_x_xyyyyy[j] - fr * g2_x_yyyyy[j];

                g_xx_yyyyz[j] = g1_x_xyyyyz[j] - fr * g2_x_yyyyz[j];

                g_xx_yyyzz[j] = g1_x_xyyyzz[j] - fr * g2_x_yyyzz[j];

                g_xx_yyzzz[j] = g1_x_xyyzzz[j] - fr * g2_x_yyzzz[j];

                g_xx_yzzzz[j] = g1_x_xyzzzz[j] - fr * g2_x_yzzzz[j];

                g_xx_zzzzz[j] = g1_x_xzzzzz[j] - fr * g2_x_zzzzz[j];

                g_xy_xxxxx[j] = g1_y_xxxxxx[j] - fr * g2_y_xxxxx[j];

                g_xy_xxxxy[j] = g1_y_xxxxxy[j] - fr * g2_y_xxxxy[j];

                g_xy_xxxxz[j] = g1_y_xxxxxz[j] - fr * g2_y_xxxxz[j];

                g_xy_xxxyy[j] = g1_y_xxxxyy[j] - fr * g2_y_xxxyy[j];

                g_xy_xxxyz[j] = g1_y_xxxxyz[j] - fr * g2_y_xxxyz[j];

                g_xy_xxxzz[j] = g1_y_xxxxzz[j] - fr * g2_y_xxxzz[j];

                g_xy_xxyyy[j] = g1_y_xxxyyy[j] - fr * g2_y_xxyyy[j];

                g_xy_xxyyz[j] = g1_y_xxxyyz[j] - fr * g2_y_xxyyz[j];

                g_xy_xxyzz[j] = g1_y_xxxyzz[j] - fr * g2_y_xxyzz[j];

                g_xy_xxzzz[j] = g1_y_xxxzzz[j] - fr * g2_y_xxzzz[j];

                g_xy_xyyyy[j] = g1_y_xxyyyy[j] - fr * g2_y_xyyyy[j];

                g_xy_xyyyz[j] = g1_y_xxyyyz[j] - fr * g2_y_xyyyz[j];

                g_xy_xyyzz[j] = g1_y_xxyyzz[j] - fr * g2_y_xyyzz[j];

                g_xy_xyzzz[j] = g1_y_xxyzzz[j] - fr * g2_y_xyzzz[j];

                g_xy_xzzzz[j] = g1_y_xxzzzz[j] - fr * g2_y_xzzzz[j];

                g_xy_yyyyy[j] = g1_y_xyyyyy[j] - fr * g2_y_yyyyy[j];

                g_xy_yyyyz[j] = g1_y_xyyyyz[j] - fr * g2_y_yyyyz[j];

                g_xy_yyyzz[j] = g1_y_xyyyzz[j] - fr * g2_y_yyyzz[j];

                g_xy_yyzzz[j] = g1_y_xyyzzz[j] - fr * g2_y_yyzzz[j];

                g_xy_yzzzz[j] = g1_y_xyzzzz[j] - fr * g2_y_yzzzz[j];

                g_xy_zzzzz[j] = g1_y_xzzzzz[j] - fr * g2_y_zzzzz[j];

                g_xz_xxxxx[j] = g1_z_xxxxxx[j] - fr * g2_z_xxxxx[j];

                g_xz_xxxxy[j] = g1_z_xxxxxy[j] - fr * g2_z_xxxxy[j];

                g_xz_xxxxz[j] = g1_z_xxxxxz[j] - fr * g2_z_xxxxz[j];

                g_xz_xxxyy[j] = g1_z_xxxxyy[j] - fr * g2_z_xxxyy[j];

                g_xz_xxxyz[j] = g1_z_xxxxyz[j] - fr * g2_z_xxxyz[j];

                g_xz_xxxzz[j] = g1_z_xxxxzz[j] - fr * g2_z_xxxzz[j];

                g_xz_xxyyy[j] = g1_z_xxxyyy[j] - fr * g2_z_xxyyy[j];

                g_xz_xxyyz[j] = g1_z_xxxyyz[j] - fr * g2_z_xxyyz[j];

                g_xz_xxyzz[j] = g1_z_xxxyzz[j] - fr * g2_z_xxyzz[j];

                g_xz_xxzzz[j] = g1_z_xxxzzz[j] - fr * g2_z_xxzzz[j];

                g_xz_xyyyy[j] = g1_z_xxyyyy[j] - fr * g2_z_xyyyy[j];

                g_xz_xyyyz[j] = g1_z_xxyyyz[j] - fr * g2_z_xyyyz[j];

                g_xz_xyyzz[j] = g1_z_xxyyzz[j] - fr * g2_z_xyyzz[j];

                g_xz_xyzzz[j] = g1_z_xxyzzz[j] - fr * g2_z_xyzzz[j];

                g_xz_xzzzz[j] = g1_z_xxzzzz[j] - fr * g2_z_xzzzz[j];

                g_xz_yyyyy[j] = g1_z_xyyyyy[j] - fr * g2_z_yyyyy[j];

                g_xz_yyyyz[j] = g1_z_xyyyyz[j] - fr * g2_z_yyyyz[j];

                g_xz_yyyzz[j] = g1_z_xyyyzz[j] - fr * g2_z_yyyzz[j];

                g_xz_yyzzz[j] = g1_z_xyyzzz[j] - fr * g2_z_yyzzz[j];

                g_xz_yzzzz[j] = g1_z_xyzzzz[j] - fr * g2_z_yzzzz[j];

                g_xz_zzzzz[j] = g1_z_xzzzzz[j] - fr * g2_z_zzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yy_xxxxx[j] = g1_y_xxxxxy[j] - fr * g2_y_xxxxx[j];

                g_yy_xxxxy[j] = g1_y_xxxxyy[j] - fr * g2_y_xxxxy[j];

                g_yy_xxxxz[j] = g1_y_xxxxyz[j] - fr * g2_y_xxxxz[j];

                g_yy_xxxyy[j] = g1_y_xxxyyy[j] - fr * g2_y_xxxyy[j];

                g_yy_xxxyz[j] = g1_y_xxxyyz[j] - fr * g2_y_xxxyz[j];

                g_yy_xxxzz[j] = g1_y_xxxyzz[j] - fr * g2_y_xxxzz[j];

                g_yy_xxyyy[j] = g1_y_xxyyyy[j] - fr * g2_y_xxyyy[j];

                g_yy_xxyyz[j] = g1_y_xxyyyz[j] - fr * g2_y_xxyyz[j];

                g_yy_xxyzz[j] = g1_y_xxyyzz[j] - fr * g2_y_xxyzz[j];

                g_yy_xxzzz[j] = g1_y_xxyzzz[j] - fr * g2_y_xxzzz[j];

                g_yy_xyyyy[j] = g1_y_xyyyyy[j] - fr * g2_y_xyyyy[j];

                g_yy_xyyyz[j] = g1_y_xyyyyz[j] - fr * g2_y_xyyyz[j];

                g_yy_xyyzz[j] = g1_y_xyyyzz[j] - fr * g2_y_xyyzz[j];

                g_yy_xyzzz[j] = g1_y_xyyzzz[j] - fr * g2_y_xyzzz[j];

                g_yy_xzzzz[j] = g1_y_xyzzzz[j] - fr * g2_y_xzzzz[j];

                g_yy_yyyyy[j] = g1_y_yyyyyy[j] - fr * g2_y_yyyyy[j];

                g_yy_yyyyz[j] = g1_y_yyyyyz[j] - fr * g2_y_yyyyz[j];

                g_yy_yyyzz[j] = g1_y_yyyyzz[j] - fr * g2_y_yyyzz[j];

                g_yy_yyzzz[j] = g1_y_yyyzzz[j] - fr * g2_y_yyzzz[j];

                g_yy_yzzzz[j] = g1_y_yyzzzz[j] - fr * g2_y_yzzzz[j];

                g_yy_zzzzz[j] = g1_y_yzzzzz[j] - fr * g2_y_zzzzz[j];

                g_yz_xxxxx[j] = g1_z_xxxxxy[j] - fr * g2_z_xxxxx[j];

                g_yz_xxxxy[j] = g1_z_xxxxyy[j] - fr * g2_z_xxxxy[j];

                g_yz_xxxxz[j] = g1_z_xxxxyz[j] - fr * g2_z_xxxxz[j];

                g_yz_xxxyy[j] = g1_z_xxxyyy[j] - fr * g2_z_xxxyy[j];

                g_yz_xxxyz[j] = g1_z_xxxyyz[j] - fr * g2_z_xxxyz[j];

                g_yz_xxxzz[j] = g1_z_xxxyzz[j] - fr * g2_z_xxxzz[j];

                g_yz_xxyyy[j] = g1_z_xxyyyy[j] - fr * g2_z_xxyyy[j];

                g_yz_xxyyz[j] = g1_z_xxyyyz[j] - fr * g2_z_xxyyz[j];

                g_yz_xxyzz[j] = g1_z_xxyyzz[j] - fr * g2_z_xxyzz[j];

                g_yz_xxzzz[j] = g1_z_xxyzzz[j] - fr * g2_z_xxzzz[j];

                g_yz_xyyyy[j] = g1_z_xyyyyy[j] - fr * g2_z_xyyyy[j];

                g_yz_xyyyz[j] = g1_z_xyyyyz[j] - fr * g2_z_xyyyz[j];

                g_yz_xyyzz[j] = g1_z_xyyyzz[j] - fr * g2_z_xyyzz[j];

                g_yz_xyzzz[j] = g1_z_xyyzzz[j] - fr * g2_z_xyzzz[j];

                g_yz_xzzzz[j] = g1_z_xyzzzz[j] - fr * g2_z_xzzzz[j];

                g_yz_yyyyy[j] = g1_z_yyyyyy[j] - fr * g2_z_yyyyy[j];

                g_yz_yyyyz[j] = g1_z_yyyyyz[j] - fr * g2_z_yyyyz[j];

                g_yz_yyyzz[j] = g1_z_yyyyzz[j] - fr * g2_z_yyyzz[j];

                g_yz_yyzzz[j] = g1_z_yyyzzz[j] - fr * g2_z_yyzzz[j];

                g_yz_yzzzz[j] = g1_z_yyzzzz[j] - fr * g2_z_yzzzz[j];

                g_yz_zzzzz[j] = g1_z_yzzzzz[j] - fr * g2_z_zzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zz_xxxxx[j] = g1_z_xxxxxz[j] - fr * g2_z_xxxxx[j];

                g_zz_xxxxy[j] = g1_z_xxxxyz[j] - fr * g2_z_xxxxy[j];

                g_zz_xxxxz[j] = g1_z_xxxxzz[j] - fr * g2_z_xxxxz[j];

                g_zz_xxxyy[j] = g1_z_xxxyyz[j] - fr * g2_z_xxxyy[j];

                g_zz_xxxyz[j] = g1_z_xxxyzz[j] - fr * g2_z_xxxyz[j];

                g_zz_xxxzz[j] = g1_z_xxxzzz[j] - fr * g2_z_xxxzz[j];

                g_zz_xxyyy[j] = g1_z_xxyyyz[j] - fr * g2_z_xxyyy[j];

                g_zz_xxyyz[j] = g1_z_xxyyzz[j] - fr * g2_z_xxyyz[j];

                g_zz_xxyzz[j] = g1_z_xxyzzz[j] - fr * g2_z_xxyzz[j];

                g_zz_xxzzz[j] = g1_z_xxzzzz[j] - fr * g2_z_xxzzz[j];

                g_zz_xyyyy[j] = g1_z_xyyyyz[j] - fr * g2_z_xyyyy[j];

                g_zz_xyyyz[j] = g1_z_xyyyzz[j] - fr * g2_z_xyyyz[j];

                g_zz_xyyzz[j] = g1_z_xyyzzz[j] - fr * g2_z_xyyzz[j];

                g_zz_xyzzz[j] = g1_z_xyzzzz[j] - fr * g2_z_xyzzz[j];

                g_zz_xzzzz[j] = g1_z_xzzzzz[j] - fr * g2_z_xzzzz[j];

                g_zz_yyyyy[j] = g1_z_yyyyyz[j] - fr * g2_z_yyyyy[j];

                g_zz_yyyyz[j] = g1_z_yyyyzz[j] - fr * g2_z_yyyyz[j];

                g_zz_yyyzz[j] = g1_z_yyyzzz[j] - fr * g2_z_yyyzz[j];

                g_zz_yyzzz[j] = g1_z_yyzzzz[j] - fr * g2_z_yyzzz[j];

                g_zz_yzzzz[j] = g1_z_yzzzzz[j] - fr * g2_z_yzzzz[j];

                g_zz_zzzzz[j] = g1_z_zzzzzz[j] - fr * g2_z_zzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXDI(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 2, 6})) return;

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
                                             {braAngularMomentum, 2, 6});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 7});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 1, 6});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|PI)^(m) integrals

            auto g2_x_xxxxxx = contrBuffer.data(g2off + 84 * i);

            auto g2_x_xxxxxy = contrBuffer.data(g2off + 84 * i + 1);

            auto g2_x_xxxxxz = contrBuffer.data(g2off + 84 * i + 2);

            auto g2_x_xxxxyy = contrBuffer.data(g2off + 84 * i + 3);

            auto g2_x_xxxxyz = contrBuffer.data(g2off + 84 * i + 4);

            auto g2_x_xxxxzz = contrBuffer.data(g2off + 84 * i + 5);

            auto g2_x_xxxyyy = contrBuffer.data(g2off + 84 * i + 6);

            auto g2_x_xxxyyz = contrBuffer.data(g2off + 84 * i + 7);

            auto g2_x_xxxyzz = contrBuffer.data(g2off + 84 * i + 8);

            auto g2_x_xxxzzz = contrBuffer.data(g2off + 84 * i + 9);

            auto g2_x_xxyyyy = contrBuffer.data(g2off + 84 * i + 10);

            auto g2_x_xxyyyz = contrBuffer.data(g2off + 84 * i + 11);

            auto g2_x_xxyyzz = contrBuffer.data(g2off + 84 * i + 12);

            auto g2_x_xxyzzz = contrBuffer.data(g2off + 84 * i + 13);

            auto g2_x_xxzzzz = contrBuffer.data(g2off + 84 * i + 14);

            auto g2_x_xyyyyy = contrBuffer.data(g2off + 84 * i + 15);

            auto g2_x_xyyyyz = contrBuffer.data(g2off + 84 * i + 16);

            auto g2_x_xyyyzz = contrBuffer.data(g2off + 84 * i + 17);

            auto g2_x_xyyzzz = contrBuffer.data(g2off + 84 * i + 18);

            auto g2_x_xyzzzz = contrBuffer.data(g2off + 84 * i + 19);

            auto g2_x_xzzzzz = contrBuffer.data(g2off + 84 * i + 20);

            auto g2_x_yyyyyy = contrBuffer.data(g2off + 84 * i + 21);

            auto g2_x_yyyyyz = contrBuffer.data(g2off + 84 * i + 22);

            auto g2_x_yyyyzz = contrBuffer.data(g2off + 84 * i + 23);

            auto g2_x_yyyzzz = contrBuffer.data(g2off + 84 * i + 24);

            auto g2_x_yyzzzz = contrBuffer.data(g2off + 84 * i + 25);

            auto g2_x_yzzzzz = contrBuffer.data(g2off + 84 * i + 26);

            auto g2_x_zzzzzz = contrBuffer.data(g2off + 84 * i + 27);

            auto g2_y_xxxxxx = contrBuffer.data(g2off + 84 * i + 28);

            auto g2_y_xxxxxy = contrBuffer.data(g2off + 84 * i + 29);

            auto g2_y_xxxxxz = contrBuffer.data(g2off + 84 * i + 30);

            auto g2_y_xxxxyy = contrBuffer.data(g2off + 84 * i + 31);

            auto g2_y_xxxxyz = contrBuffer.data(g2off + 84 * i + 32);

            auto g2_y_xxxxzz = contrBuffer.data(g2off + 84 * i + 33);

            auto g2_y_xxxyyy = contrBuffer.data(g2off + 84 * i + 34);

            auto g2_y_xxxyyz = contrBuffer.data(g2off + 84 * i + 35);

            auto g2_y_xxxyzz = contrBuffer.data(g2off + 84 * i + 36);

            auto g2_y_xxxzzz = contrBuffer.data(g2off + 84 * i + 37);

            auto g2_y_xxyyyy = contrBuffer.data(g2off + 84 * i + 38);

            auto g2_y_xxyyyz = contrBuffer.data(g2off + 84 * i + 39);

            auto g2_y_xxyyzz = contrBuffer.data(g2off + 84 * i + 40);

            auto g2_y_xxyzzz = contrBuffer.data(g2off + 84 * i + 41);

            auto g2_y_xxzzzz = contrBuffer.data(g2off + 84 * i + 42);

            auto g2_y_xyyyyy = contrBuffer.data(g2off + 84 * i + 43);

            auto g2_y_xyyyyz = contrBuffer.data(g2off + 84 * i + 44);

            auto g2_y_xyyyzz = contrBuffer.data(g2off + 84 * i + 45);

            auto g2_y_xyyzzz = contrBuffer.data(g2off + 84 * i + 46);

            auto g2_y_xyzzzz = contrBuffer.data(g2off + 84 * i + 47);

            auto g2_y_xzzzzz = contrBuffer.data(g2off + 84 * i + 48);

            auto g2_y_yyyyyy = contrBuffer.data(g2off + 84 * i + 49);

            auto g2_y_yyyyyz = contrBuffer.data(g2off + 84 * i + 50);

            auto g2_y_yyyyzz = contrBuffer.data(g2off + 84 * i + 51);

            auto g2_y_yyyzzz = contrBuffer.data(g2off + 84 * i + 52);

            auto g2_y_yyzzzz = contrBuffer.data(g2off + 84 * i + 53);

            auto g2_y_yzzzzz = contrBuffer.data(g2off + 84 * i + 54);

            auto g2_y_zzzzzz = contrBuffer.data(g2off + 84 * i + 55);

            auto g2_z_xxxxxx = contrBuffer.data(g2off + 84 * i + 56);

            auto g2_z_xxxxxy = contrBuffer.data(g2off + 84 * i + 57);

            auto g2_z_xxxxxz = contrBuffer.data(g2off + 84 * i + 58);

            auto g2_z_xxxxyy = contrBuffer.data(g2off + 84 * i + 59);

            auto g2_z_xxxxyz = contrBuffer.data(g2off + 84 * i + 60);

            auto g2_z_xxxxzz = contrBuffer.data(g2off + 84 * i + 61);

            auto g2_z_xxxyyy = contrBuffer.data(g2off + 84 * i + 62);

            auto g2_z_xxxyyz = contrBuffer.data(g2off + 84 * i + 63);

            auto g2_z_xxxyzz = contrBuffer.data(g2off + 84 * i + 64);

            auto g2_z_xxxzzz = contrBuffer.data(g2off + 84 * i + 65);

            auto g2_z_xxyyyy = contrBuffer.data(g2off + 84 * i + 66);

            auto g2_z_xxyyyz = contrBuffer.data(g2off + 84 * i + 67);

            auto g2_z_xxyyzz = contrBuffer.data(g2off + 84 * i + 68);

            auto g2_z_xxyzzz = contrBuffer.data(g2off + 84 * i + 69);

            auto g2_z_xxzzzz = contrBuffer.data(g2off + 84 * i + 70);

            auto g2_z_xyyyyy = contrBuffer.data(g2off + 84 * i + 71);

            auto g2_z_xyyyyz = contrBuffer.data(g2off + 84 * i + 72);

            auto g2_z_xyyyzz = contrBuffer.data(g2off + 84 * i + 73);

            auto g2_z_xyyzzz = contrBuffer.data(g2off + 84 * i + 74);

            auto g2_z_xyzzzz = contrBuffer.data(g2off + 84 * i + 75);

            auto g2_z_xzzzzz = contrBuffer.data(g2off + 84 * i + 76);

            auto g2_z_yyyyyy = contrBuffer.data(g2off + 84 * i + 77);

            auto g2_z_yyyyyz = contrBuffer.data(g2off + 84 * i + 78);

            auto g2_z_yyyyzz = contrBuffer.data(g2off + 84 * i + 79);

            auto g2_z_yyyzzz = contrBuffer.data(g2off + 84 * i + 80);

            auto g2_z_yyzzzz = contrBuffer.data(g2off + 84 * i + 81);

            auto g2_z_yzzzzz = contrBuffer.data(g2off + 84 * i + 82);

            auto g2_z_zzzzzz = contrBuffer.data(g2off + 84 * i + 83);

            // set up pointers to (X|g(r,r')|PK)^(m) integrals

            auto g1_x_xxxxxxx = contrBuffer.data(g1off + 108 * i);

            auto g1_x_xxxxxxy = contrBuffer.data(g1off + 108 * i + 1);

            auto g1_x_xxxxxxz = contrBuffer.data(g1off + 108 * i + 2);

            auto g1_x_xxxxxyy = contrBuffer.data(g1off + 108 * i + 3);

            auto g1_x_xxxxxyz = contrBuffer.data(g1off + 108 * i + 4);

            auto g1_x_xxxxxzz = contrBuffer.data(g1off + 108 * i + 5);

            auto g1_x_xxxxyyy = contrBuffer.data(g1off + 108 * i + 6);

            auto g1_x_xxxxyyz = contrBuffer.data(g1off + 108 * i + 7);

            auto g1_x_xxxxyzz = contrBuffer.data(g1off + 108 * i + 8);

            auto g1_x_xxxxzzz = contrBuffer.data(g1off + 108 * i + 9);

            auto g1_x_xxxyyyy = contrBuffer.data(g1off + 108 * i + 10);

            auto g1_x_xxxyyyz = contrBuffer.data(g1off + 108 * i + 11);

            auto g1_x_xxxyyzz = contrBuffer.data(g1off + 108 * i + 12);

            auto g1_x_xxxyzzz = contrBuffer.data(g1off + 108 * i + 13);

            auto g1_x_xxxzzzz = contrBuffer.data(g1off + 108 * i + 14);

            auto g1_x_xxyyyyy = contrBuffer.data(g1off + 108 * i + 15);

            auto g1_x_xxyyyyz = contrBuffer.data(g1off + 108 * i + 16);

            auto g1_x_xxyyyzz = contrBuffer.data(g1off + 108 * i + 17);

            auto g1_x_xxyyzzz = contrBuffer.data(g1off + 108 * i + 18);

            auto g1_x_xxyzzzz = contrBuffer.data(g1off + 108 * i + 19);

            auto g1_x_xxzzzzz = contrBuffer.data(g1off + 108 * i + 20);

            auto g1_x_xyyyyyy = contrBuffer.data(g1off + 108 * i + 21);

            auto g1_x_xyyyyyz = contrBuffer.data(g1off + 108 * i + 22);

            auto g1_x_xyyyyzz = contrBuffer.data(g1off + 108 * i + 23);

            auto g1_x_xyyyzzz = contrBuffer.data(g1off + 108 * i + 24);

            auto g1_x_xyyzzzz = contrBuffer.data(g1off + 108 * i + 25);

            auto g1_x_xyzzzzz = contrBuffer.data(g1off + 108 * i + 26);

            auto g1_x_xzzzzzz = contrBuffer.data(g1off + 108 * i + 27);

            auto g1_y_xxxxxxx = contrBuffer.data(g1off + 108 * i + 36);

            auto g1_y_xxxxxxy = contrBuffer.data(g1off + 108 * i + 37);

            auto g1_y_xxxxxxz = contrBuffer.data(g1off + 108 * i + 38);

            auto g1_y_xxxxxyy = contrBuffer.data(g1off + 108 * i + 39);

            auto g1_y_xxxxxyz = contrBuffer.data(g1off + 108 * i + 40);

            auto g1_y_xxxxxzz = contrBuffer.data(g1off + 108 * i + 41);

            auto g1_y_xxxxyyy = contrBuffer.data(g1off + 108 * i + 42);

            auto g1_y_xxxxyyz = contrBuffer.data(g1off + 108 * i + 43);

            auto g1_y_xxxxyzz = contrBuffer.data(g1off + 108 * i + 44);

            auto g1_y_xxxxzzz = contrBuffer.data(g1off + 108 * i + 45);

            auto g1_y_xxxyyyy = contrBuffer.data(g1off + 108 * i + 46);

            auto g1_y_xxxyyyz = contrBuffer.data(g1off + 108 * i + 47);

            auto g1_y_xxxyyzz = contrBuffer.data(g1off + 108 * i + 48);

            auto g1_y_xxxyzzz = contrBuffer.data(g1off + 108 * i + 49);

            auto g1_y_xxxzzzz = contrBuffer.data(g1off + 108 * i + 50);

            auto g1_y_xxyyyyy = contrBuffer.data(g1off + 108 * i + 51);

            auto g1_y_xxyyyyz = contrBuffer.data(g1off + 108 * i + 52);

            auto g1_y_xxyyyzz = contrBuffer.data(g1off + 108 * i + 53);

            auto g1_y_xxyyzzz = contrBuffer.data(g1off + 108 * i + 54);

            auto g1_y_xxyzzzz = contrBuffer.data(g1off + 108 * i + 55);

            auto g1_y_xxzzzzz = contrBuffer.data(g1off + 108 * i + 56);

            auto g1_y_xyyyyyy = contrBuffer.data(g1off + 108 * i + 57);

            auto g1_y_xyyyyyz = contrBuffer.data(g1off + 108 * i + 58);

            auto g1_y_xyyyyzz = contrBuffer.data(g1off + 108 * i + 59);

            auto g1_y_xyyyzzz = contrBuffer.data(g1off + 108 * i + 60);

            auto g1_y_xyyzzzz = contrBuffer.data(g1off + 108 * i + 61);

            auto g1_y_xyzzzzz = contrBuffer.data(g1off + 108 * i + 62);

            auto g1_y_xzzzzzz = contrBuffer.data(g1off + 108 * i + 63);

            auto g1_y_yyyyyyy = contrBuffer.data(g1off + 108 * i + 64);

            auto g1_y_yyyyyyz = contrBuffer.data(g1off + 108 * i + 65);

            auto g1_y_yyyyyzz = contrBuffer.data(g1off + 108 * i + 66);

            auto g1_y_yyyyzzz = contrBuffer.data(g1off + 108 * i + 67);

            auto g1_y_yyyzzzz = contrBuffer.data(g1off + 108 * i + 68);

            auto g1_y_yyzzzzz = contrBuffer.data(g1off + 108 * i + 69);

            auto g1_y_yzzzzzz = contrBuffer.data(g1off + 108 * i + 70);

            auto g1_z_xxxxxxx = contrBuffer.data(g1off + 108 * i + 72);

            auto g1_z_xxxxxxy = contrBuffer.data(g1off + 108 * i + 73);

            auto g1_z_xxxxxxz = contrBuffer.data(g1off + 108 * i + 74);

            auto g1_z_xxxxxyy = contrBuffer.data(g1off + 108 * i + 75);

            auto g1_z_xxxxxyz = contrBuffer.data(g1off + 108 * i + 76);

            auto g1_z_xxxxxzz = contrBuffer.data(g1off + 108 * i + 77);

            auto g1_z_xxxxyyy = contrBuffer.data(g1off + 108 * i + 78);

            auto g1_z_xxxxyyz = contrBuffer.data(g1off + 108 * i + 79);

            auto g1_z_xxxxyzz = contrBuffer.data(g1off + 108 * i + 80);

            auto g1_z_xxxxzzz = contrBuffer.data(g1off + 108 * i + 81);

            auto g1_z_xxxyyyy = contrBuffer.data(g1off + 108 * i + 82);

            auto g1_z_xxxyyyz = contrBuffer.data(g1off + 108 * i + 83);

            auto g1_z_xxxyyzz = contrBuffer.data(g1off + 108 * i + 84);

            auto g1_z_xxxyzzz = contrBuffer.data(g1off + 108 * i + 85);

            auto g1_z_xxxzzzz = contrBuffer.data(g1off + 108 * i + 86);

            auto g1_z_xxyyyyy = contrBuffer.data(g1off + 108 * i + 87);

            auto g1_z_xxyyyyz = contrBuffer.data(g1off + 108 * i + 88);

            auto g1_z_xxyyyzz = contrBuffer.data(g1off + 108 * i + 89);

            auto g1_z_xxyyzzz = contrBuffer.data(g1off + 108 * i + 90);

            auto g1_z_xxyzzzz = contrBuffer.data(g1off + 108 * i + 91);

            auto g1_z_xxzzzzz = contrBuffer.data(g1off + 108 * i + 92);

            auto g1_z_xyyyyyy = contrBuffer.data(g1off + 108 * i + 93);

            auto g1_z_xyyyyyz = contrBuffer.data(g1off + 108 * i + 94);

            auto g1_z_xyyyyzz = contrBuffer.data(g1off + 108 * i + 95);

            auto g1_z_xyyyzzz = contrBuffer.data(g1off + 108 * i + 96);

            auto g1_z_xyyzzzz = contrBuffer.data(g1off + 108 * i + 97);

            auto g1_z_xyzzzzz = contrBuffer.data(g1off + 108 * i + 98);

            auto g1_z_xzzzzzz = contrBuffer.data(g1off + 108 * i + 99);

            auto g1_z_yyyyyyy = contrBuffer.data(g1off + 108 * i + 100);

            auto g1_z_yyyyyyz = contrBuffer.data(g1off + 108 * i + 101);

            auto g1_z_yyyyyzz = contrBuffer.data(g1off + 108 * i + 102);

            auto g1_z_yyyyzzz = contrBuffer.data(g1off + 108 * i + 103);

            auto g1_z_yyyzzzz = contrBuffer.data(g1off + 108 * i + 104);

            auto g1_z_yyzzzzz = contrBuffer.data(g1off + 108 * i + 105);

            auto g1_z_yzzzzzz = contrBuffer.data(g1off + 108 * i + 106);

            auto g1_z_zzzzzzz = contrBuffer.data(g1off + 108 * i + 107);

            // set up pointers to (X|g(r,r')|DI)^(m) integrals

            auto g_xx_xxxxxx = contrBuffer.data(goff + 168 * i);

            auto g_xx_xxxxxy = contrBuffer.data(goff + 168 * i + 1);

            auto g_xx_xxxxxz = contrBuffer.data(goff + 168 * i + 2);

            auto g_xx_xxxxyy = contrBuffer.data(goff + 168 * i + 3);

            auto g_xx_xxxxyz = contrBuffer.data(goff + 168 * i + 4);

            auto g_xx_xxxxzz = contrBuffer.data(goff + 168 * i + 5);

            auto g_xx_xxxyyy = contrBuffer.data(goff + 168 * i + 6);

            auto g_xx_xxxyyz = contrBuffer.data(goff + 168 * i + 7);

            auto g_xx_xxxyzz = contrBuffer.data(goff + 168 * i + 8);

            auto g_xx_xxxzzz = contrBuffer.data(goff + 168 * i + 9);

            auto g_xx_xxyyyy = contrBuffer.data(goff + 168 * i + 10);

            auto g_xx_xxyyyz = contrBuffer.data(goff + 168 * i + 11);

            auto g_xx_xxyyzz = contrBuffer.data(goff + 168 * i + 12);

            auto g_xx_xxyzzz = contrBuffer.data(goff + 168 * i + 13);

            auto g_xx_xxzzzz = contrBuffer.data(goff + 168 * i + 14);

            auto g_xx_xyyyyy = contrBuffer.data(goff + 168 * i + 15);

            auto g_xx_xyyyyz = contrBuffer.data(goff + 168 * i + 16);

            auto g_xx_xyyyzz = contrBuffer.data(goff + 168 * i + 17);

            auto g_xx_xyyzzz = contrBuffer.data(goff + 168 * i + 18);

            auto g_xx_xyzzzz = contrBuffer.data(goff + 168 * i + 19);

            auto g_xx_xzzzzz = contrBuffer.data(goff + 168 * i + 20);

            auto g_xx_yyyyyy = contrBuffer.data(goff + 168 * i + 21);

            auto g_xx_yyyyyz = contrBuffer.data(goff + 168 * i + 22);

            auto g_xx_yyyyzz = contrBuffer.data(goff + 168 * i + 23);

            auto g_xx_yyyzzz = contrBuffer.data(goff + 168 * i + 24);

            auto g_xx_yyzzzz = contrBuffer.data(goff + 168 * i + 25);

            auto g_xx_yzzzzz = contrBuffer.data(goff + 168 * i + 26);

            auto g_xx_zzzzzz = contrBuffer.data(goff + 168 * i + 27);

            auto g_xy_xxxxxx = contrBuffer.data(goff + 168 * i + 28);

            auto g_xy_xxxxxy = contrBuffer.data(goff + 168 * i + 29);

            auto g_xy_xxxxxz = contrBuffer.data(goff + 168 * i + 30);

            auto g_xy_xxxxyy = contrBuffer.data(goff + 168 * i + 31);

            auto g_xy_xxxxyz = contrBuffer.data(goff + 168 * i + 32);

            auto g_xy_xxxxzz = contrBuffer.data(goff + 168 * i + 33);

            auto g_xy_xxxyyy = contrBuffer.data(goff + 168 * i + 34);

            auto g_xy_xxxyyz = contrBuffer.data(goff + 168 * i + 35);

            auto g_xy_xxxyzz = contrBuffer.data(goff + 168 * i + 36);

            auto g_xy_xxxzzz = contrBuffer.data(goff + 168 * i + 37);

            auto g_xy_xxyyyy = contrBuffer.data(goff + 168 * i + 38);

            auto g_xy_xxyyyz = contrBuffer.data(goff + 168 * i + 39);

            auto g_xy_xxyyzz = contrBuffer.data(goff + 168 * i + 40);

            auto g_xy_xxyzzz = contrBuffer.data(goff + 168 * i + 41);

            auto g_xy_xxzzzz = contrBuffer.data(goff + 168 * i + 42);

            auto g_xy_xyyyyy = contrBuffer.data(goff + 168 * i + 43);

            auto g_xy_xyyyyz = contrBuffer.data(goff + 168 * i + 44);

            auto g_xy_xyyyzz = contrBuffer.data(goff + 168 * i + 45);

            auto g_xy_xyyzzz = contrBuffer.data(goff + 168 * i + 46);

            auto g_xy_xyzzzz = contrBuffer.data(goff + 168 * i + 47);

            auto g_xy_xzzzzz = contrBuffer.data(goff + 168 * i + 48);

            auto g_xy_yyyyyy = contrBuffer.data(goff + 168 * i + 49);

            auto g_xy_yyyyyz = contrBuffer.data(goff + 168 * i + 50);

            auto g_xy_yyyyzz = contrBuffer.data(goff + 168 * i + 51);

            auto g_xy_yyyzzz = contrBuffer.data(goff + 168 * i + 52);

            auto g_xy_yyzzzz = contrBuffer.data(goff + 168 * i + 53);

            auto g_xy_yzzzzz = contrBuffer.data(goff + 168 * i + 54);

            auto g_xy_zzzzzz = contrBuffer.data(goff + 168 * i + 55);

            auto g_xz_xxxxxx = contrBuffer.data(goff + 168 * i + 56);

            auto g_xz_xxxxxy = contrBuffer.data(goff + 168 * i + 57);

            auto g_xz_xxxxxz = contrBuffer.data(goff + 168 * i + 58);

            auto g_xz_xxxxyy = contrBuffer.data(goff + 168 * i + 59);

            auto g_xz_xxxxyz = contrBuffer.data(goff + 168 * i + 60);

            auto g_xz_xxxxzz = contrBuffer.data(goff + 168 * i + 61);

            auto g_xz_xxxyyy = contrBuffer.data(goff + 168 * i + 62);

            auto g_xz_xxxyyz = contrBuffer.data(goff + 168 * i + 63);

            auto g_xz_xxxyzz = contrBuffer.data(goff + 168 * i + 64);

            auto g_xz_xxxzzz = contrBuffer.data(goff + 168 * i + 65);

            auto g_xz_xxyyyy = contrBuffer.data(goff + 168 * i + 66);

            auto g_xz_xxyyyz = contrBuffer.data(goff + 168 * i + 67);

            auto g_xz_xxyyzz = contrBuffer.data(goff + 168 * i + 68);

            auto g_xz_xxyzzz = contrBuffer.data(goff + 168 * i + 69);

            auto g_xz_xxzzzz = contrBuffer.data(goff + 168 * i + 70);

            auto g_xz_xyyyyy = contrBuffer.data(goff + 168 * i + 71);

            auto g_xz_xyyyyz = contrBuffer.data(goff + 168 * i + 72);

            auto g_xz_xyyyzz = contrBuffer.data(goff + 168 * i + 73);

            auto g_xz_xyyzzz = contrBuffer.data(goff + 168 * i + 74);

            auto g_xz_xyzzzz = contrBuffer.data(goff + 168 * i + 75);

            auto g_xz_xzzzzz = contrBuffer.data(goff + 168 * i + 76);

            auto g_xz_yyyyyy = contrBuffer.data(goff + 168 * i + 77);

            auto g_xz_yyyyyz = contrBuffer.data(goff + 168 * i + 78);

            auto g_xz_yyyyzz = contrBuffer.data(goff + 168 * i + 79);

            auto g_xz_yyyzzz = contrBuffer.data(goff + 168 * i + 80);

            auto g_xz_yyzzzz = contrBuffer.data(goff + 168 * i + 81);

            auto g_xz_yzzzzz = contrBuffer.data(goff + 168 * i + 82);

            auto g_xz_zzzzzz = contrBuffer.data(goff + 168 * i + 83);

            auto g_yy_xxxxxx = contrBuffer.data(goff + 168 * i + 84);

            auto g_yy_xxxxxy = contrBuffer.data(goff + 168 * i + 85);

            auto g_yy_xxxxxz = contrBuffer.data(goff + 168 * i + 86);

            auto g_yy_xxxxyy = contrBuffer.data(goff + 168 * i + 87);

            auto g_yy_xxxxyz = contrBuffer.data(goff + 168 * i + 88);

            auto g_yy_xxxxzz = contrBuffer.data(goff + 168 * i + 89);

            auto g_yy_xxxyyy = contrBuffer.data(goff + 168 * i + 90);

            auto g_yy_xxxyyz = contrBuffer.data(goff + 168 * i + 91);

            auto g_yy_xxxyzz = contrBuffer.data(goff + 168 * i + 92);

            auto g_yy_xxxzzz = contrBuffer.data(goff + 168 * i + 93);

            auto g_yy_xxyyyy = contrBuffer.data(goff + 168 * i + 94);

            auto g_yy_xxyyyz = contrBuffer.data(goff + 168 * i + 95);

            auto g_yy_xxyyzz = contrBuffer.data(goff + 168 * i + 96);

            auto g_yy_xxyzzz = contrBuffer.data(goff + 168 * i + 97);

            auto g_yy_xxzzzz = contrBuffer.data(goff + 168 * i + 98);

            auto g_yy_xyyyyy = contrBuffer.data(goff + 168 * i + 99);

            auto g_yy_xyyyyz = contrBuffer.data(goff + 168 * i + 100);

            auto g_yy_xyyyzz = contrBuffer.data(goff + 168 * i + 101);

            auto g_yy_xyyzzz = contrBuffer.data(goff + 168 * i + 102);

            auto g_yy_xyzzzz = contrBuffer.data(goff + 168 * i + 103);

            auto g_yy_xzzzzz = contrBuffer.data(goff + 168 * i + 104);

            auto g_yy_yyyyyy = contrBuffer.data(goff + 168 * i + 105);

            auto g_yy_yyyyyz = contrBuffer.data(goff + 168 * i + 106);

            auto g_yy_yyyyzz = contrBuffer.data(goff + 168 * i + 107);

            auto g_yy_yyyzzz = contrBuffer.data(goff + 168 * i + 108);

            auto g_yy_yyzzzz = contrBuffer.data(goff + 168 * i + 109);

            auto g_yy_yzzzzz = contrBuffer.data(goff + 168 * i + 110);

            auto g_yy_zzzzzz = contrBuffer.data(goff + 168 * i + 111);

            auto g_yz_xxxxxx = contrBuffer.data(goff + 168 * i + 112);

            auto g_yz_xxxxxy = contrBuffer.data(goff + 168 * i + 113);

            auto g_yz_xxxxxz = contrBuffer.data(goff + 168 * i + 114);

            auto g_yz_xxxxyy = contrBuffer.data(goff + 168 * i + 115);

            auto g_yz_xxxxyz = contrBuffer.data(goff + 168 * i + 116);

            auto g_yz_xxxxzz = contrBuffer.data(goff + 168 * i + 117);

            auto g_yz_xxxyyy = contrBuffer.data(goff + 168 * i + 118);

            auto g_yz_xxxyyz = contrBuffer.data(goff + 168 * i + 119);

            auto g_yz_xxxyzz = contrBuffer.data(goff + 168 * i + 120);

            auto g_yz_xxxzzz = contrBuffer.data(goff + 168 * i + 121);

            auto g_yz_xxyyyy = contrBuffer.data(goff + 168 * i + 122);

            auto g_yz_xxyyyz = contrBuffer.data(goff + 168 * i + 123);

            auto g_yz_xxyyzz = contrBuffer.data(goff + 168 * i + 124);

            auto g_yz_xxyzzz = contrBuffer.data(goff + 168 * i + 125);

            auto g_yz_xxzzzz = contrBuffer.data(goff + 168 * i + 126);

            auto g_yz_xyyyyy = contrBuffer.data(goff + 168 * i + 127);

            auto g_yz_xyyyyz = contrBuffer.data(goff + 168 * i + 128);

            auto g_yz_xyyyzz = contrBuffer.data(goff + 168 * i + 129);

            auto g_yz_xyyzzz = contrBuffer.data(goff + 168 * i + 130);

            auto g_yz_xyzzzz = contrBuffer.data(goff + 168 * i + 131);

            auto g_yz_xzzzzz = contrBuffer.data(goff + 168 * i + 132);

            auto g_yz_yyyyyy = contrBuffer.data(goff + 168 * i + 133);

            auto g_yz_yyyyyz = contrBuffer.data(goff + 168 * i + 134);

            auto g_yz_yyyyzz = contrBuffer.data(goff + 168 * i + 135);

            auto g_yz_yyyzzz = contrBuffer.data(goff + 168 * i + 136);

            auto g_yz_yyzzzz = contrBuffer.data(goff + 168 * i + 137);

            auto g_yz_yzzzzz = contrBuffer.data(goff + 168 * i + 138);

            auto g_yz_zzzzzz = contrBuffer.data(goff + 168 * i + 139);

            auto g_zz_xxxxxx = contrBuffer.data(goff + 168 * i + 140);

            auto g_zz_xxxxxy = contrBuffer.data(goff + 168 * i + 141);

            auto g_zz_xxxxxz = contrBuffer.data(goff + 168 * i + 142);

            auto g_zz_xxxxyy = contrBuffer.data(goff + 168 * i + 143);

            auto g_zz_xxxxyz = contrBuffer.data(goff + 168 * i + 144);

            auto g_zz_xxxxzz = contrBuffer.data(goff + 168 * i + 145);

            auto g_zz_xxxyyy = contrBuffer.data(goff + 168 * i + 146);

            auto g_zz_xxxyyz = contrBuffer.data(goff + 168 * i + 147);

            auto g_zz_xxxyzz = contrBuffer.data(goff + 168 * i + 148);

            auto g_zz_xxxzzz = contrBuffer.data(goff + 168 * i + 149);

            auto g_zz_xxyyyy = contrBuffer.data(goff + 168 * i + 150);

            auto g_zz_xxyyyz = contrBuffer.data(goff + 168 * i + 151);

            auto g_zz_xxyyzz = contrBuffer.data(goff + 168 * i + 152);

            auto g_zz_xxyzzz = contrBuffer.data(goff + 168 * i + 153);

            auto g_zz_xxzzzz = contrBuffer.data(goff + 168 * i + 154);

            auto g_zz_xyyyyy = contrBuffer.data(goff + 168 * i + 155);

            auto g_zz_xyyyyz = contrBuffer.data(goff + 168 * i + 156);

            auto g_zz_xyyyzz = contrBuffer.data(goff + 168 * i + 157);

            auto g_zz_xyyzzz = contrBuffer.data(goff + 168 * i + 158);

            auto g_zz_xyzzzz = contrBuffer.data(goff + 168 * i + 159);

            auto g_zz_xzzzzz = contrBuffer.data(goff + 168 * i + 160);

            auto g_zz_yyyyyy = contrBuffer.data(goff + 168 * i + 161);

            auto g_zz_yyyyyz = contrBuffer.data(goff + 168 * i + 162);

            auto g_zz_yyyyzz = contrBuffer.data(goff + 168 * i + 163);

            auto g_zz_yyyzzz = contrBuffer.data(goff + 168 * i + 164);

            auto g_zz_yyzzzz = contrBuffer.data(goff + 168 * i + 165);

            auto g_zz_yzzzzz = contrBuffer.data(goff + 168 * i + 166);

            auto g_zz_zzzzzz = contrBuffer.data(goff + 168 * i + 167);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxxxx, g2_x_xxxxxy,\
                                     g2_x_xxxxxz, g2_x_xxxxyy, g2_x_xxxxyz, g2_x_xxxxzz,\
                                     g2_x_xxxyyy, g2_x_xxxyyz, g2_x_xxxyzz, g2_x_xxxzzz,\
                                     g2_x_xxyyyy, g2_x_xxyyyz, g2_x_xxyyzz, g2_x_xxyzzz,\
                                     g2_x_xxzzzz, g2_x_xyyyyy, g2_x_xyyyyz, g2_x_xyyyzz,\
                                     g2_x_xyyzzz, g2_x_xyzzzz, g2_x_xzzzzz, g2_x_yyyyyy,\
                                     g2_x_yyyyyz, g2_x_yyyyzz, g2_x_yyyzzz, g2_x_yyzzzz,\
                                     g2_x_yzzzzz, g2_x_zzzzzz, g2_y_xxxxxx, g2_y_xxxxxy,\
                                     g2_y_xxxxxz, g2_y_xxxxyy, g2_y_xxxxyz, g2_y_xxxxzz,\
                                     g2_y_xxxyyy, g2_y_xxxyyz, g2_y_xxxyzz, g2_y_xxxzzz,\
                                     g2_y_xxyyyy, g2_y_xxyyyz, g2_y_xxyyzz, g2_y_xxyzzz,\
                                     g2_y_xxzzzz, g2_y_xyyyyy, g2_y_xyyyyz, g2_y_xyyyzz,\
                                     g2_y_xyyzzz, g2_y_xyzzzz, g2_y_xzzzzz, g2_y_yyyyyy,\
                                     g2_y_yyyyyz, g2_y_yyyyzz, g2_y_yyyzzz, g2_y_yyzzzz,\
                                     g2_y_yzzzzz, g2_y_zzzzzz, g2_z_xxxxxx, g2_z_xxxxxy,\
                                     g2_z_xxxxxz, g2_z_xxxxyy, g2_z_xxxxyz, g2_z_xxxxzz,\
                                     g2_z_xxxyyy, g2_z_xxxyyz, g2_z_xxxyzz, g2_z_xxxzzz,\
                                     g2_z_xxyyyy, g2_z_xxyyyz, g2_z_xxyyzz, g2_z_xxyzzz,\
                                     g2_z_xxzzzz, g2_z_xyyyyy, g2_z_xyyyyz, g2_z_xyyyzz,\
                                     g2_z_xyyzzz, g2_z_xyzzzz, g2_z_xzzzzz, g2_z_yyyyyy,\
                                     g2_z_yyyyyz, g2_z_yyyyzz, g2_z_yyyzzz, g2_z_yyzzzz,\
                                     g2_z_yzzzzz, g2_z_zzzzzz, g1_x_xxxxxxx, g1_x_xxxxxxy,\
                                     g1_x_xxxxxxz, g1_x_xxxxxyy, g1_x_xxxxxyz,\
                                     g1_x_xxxxxzz, g1_x_xxxxyyy, g1_x_xxxxyyz,\
                                     g1_x_xxxxyzz, g1_x_xxxxzzz, g1_x_xxxyyyy,\
                                     g1_x_xxxyyyz, g1_x_xxxyyzz, g1_x_xxxyzzz,\
                                     g1_x_xxxzzzz, g1_x_xxyyyyy, g1_x_xxyyyyz,\
                                     g1_x_xxyyyzz, g1_x_xxyyzzz, g1_x_xxyzzzz,\
                                     g1_x_xxzzzzz, g1_x_xyyyyyy, g1_x_xyyyyyz,\
                                     g1_x_xyyyyzz, g1_x_xyyyzzz, g1_x_xyyzzzz,\
                                     g1_x_xyzzzzz, g1_x_xzzzzzz,  g1_y_xxxxxxx,\
                                     g1_y_xxxxxxy, g1_y_xxxxxxz, g1_y_xxxxxyy,\
                                     g1_y_xxxxxyz, g1_y_xxxxxzz, g1_y_xxxxyyy,\
                                     g1_y_xxxxyyz, g1_y_xxxxyzz, g1_y_xxxxzzz,\
                                     g1_y_xxxyyyy, g1_y_xxxyyyz, g1_y_xxxyyzz,\
                                     g1_y_xxxyzzz, g1_y_xxxzzzz, g1_y_xxyyyyy,\
                                     g1_y_xxyyyyz, g1_y_xxyyyzz, g1_y_xxyyzzz,\
                                     g1_y_xxyzzzz, g1_y_xxzzzzz, g1_y_xyyyyyy,\
                                     g1_y_xyyyyyz, g1_y_xyyyyzz, g1_y_xyyyzzz,\
                                     g1_y_xyyzzzz, g1_y_xyzzzzz, g1_y_xzzzzzz,\
                                     g1_y_yyyyyyy, g1_y_yyyyyyz, g1_y_yyyyyzz,\
                                     g1_y_yyyyzzz, g1_y_yyyzzzz, g1_y_yyzzzzz,\
                                     g1_y_yzzzzzz, g1_z_xxxxxxx, g1_z_xxxxxxy,\
                                     g1_z_xxxxxxz, g1_z_xxxxxyy, g1_z_xxxxxyz,\
                                     g1_z_xxxxxzz, g1_z_xxxxyyy, g1_z_xxxxyyz,\
                                     g1_z_xxxxyzz, g1_z_xxxxzzz, g1_z_xxxyyyy,\
                                     g1_z_xxxyyyz, g1_z_xxxyyzz, g1_z_xxxyzzz,\
                                     g1_z_xxxzzzz, g1_z_xxyyyyy, g1_z_xxyyyyz,\
                                     g1_z_xxyyyzz, g1_z_xxyyzzz, g1_z_xxyzzzz,\
                                     g1_z_xxzzzzz, g1_z_xyyyyyy, g1_z_xyyyyyz,\
                                     g1_z_xyyyyzz, g1_z_xyyyzzz, g1_z_xyyzzzz,\
                                     g1_z_xyzzzzz, g1_z_xzzzzzz, g1_z_yyyyyyy,\
                                     g1_z_yyyyyyz, g1_z_yyyyyzz, g1_z_yyyyzzz,\
                                     g1_z_yyyzzzz, g1_z_yyzzzzz, g1_z_yzzzzzz,\
                                     g1_z_zzzzzzz, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz,\
                                     g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz, g_xx_xxxyyy,\
                                     g_xx_xxxyyz, g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy,\
                                     g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz,\
                                     g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz, g_xx_xyyzzz,\
                                     g_xx_xyzzzz, g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz,\
                                     g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz,\
                                     g_xx_zzzzzz, g_xy_xxxxxx, g_xy_xxxxxy, g_xy_xxxxxz,\
                                     g_xy_xxxxyy, g_xy_xxxxyz, g_xy_xxxxzz, g_xy_xxxyyy,\
                                     g_xy_xxxyyz, g_xy_xxxyzz, g_xy_xxxzzz, g_xy_xxyyyy,\
                                     g_xy_xxyyyz, g_xy_xxyyzz, g_xy_xxyzzz, g_xy_xxzzzz,\
                                     g_xy_xyyyyy, g_xy_xyyyyz, g_xy_xyyyzz, g_xy_xyyzzz,\
                                     g_xy_xyzzzz, g_xy_xzzzzz, g_xy_yyyyyy, g_xy_yyyyyz,\
                                     g_xy_yyyyzz, g_xy_yyyzzz, g_xy_yyzzzz, g_xy_yzzzzz,\
                                     g_xy_zzzzzz, g_xz_xxxxxx, g_xz_xxxxxy, g_xz_xxxxxz,\
                                     g_xz_xxxxyy, g_xz_xxxxyz, g_xz_xxxxzz, g_xz_xxxyyy,\
                                     g_xz_xxxyyz, g_xz_xxxyzz, g_xz_xxxzzz, g_xz_xxyyyy,\
                                     g_xz_xxyyyz, g_xz_xxyyzz, g_xz_xxyzzz, g_xz_xxzzzz,\
                                     g_xz_xyyyyy, g_xz_xyyyyz, g_xz_xyyyzz, g_xz_xyyzzz,\
                                     g_xz_xyzzzz, g_xz_xzzzzz, g_xz_yyyyyy, g_xz_yyyyyz,\
                                     g_xz_yyyyzz, g_xz_yyyzzz, g_xz_yyzzzz, g_xz_yzzzzz,\
                                     g_xz_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz,\
                                     g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz, g_yy_xxxyyy,\
                                     g_yy_xxxyyz, g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy,\
                                     g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz,\
                                     g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz, g_yy_xyyzzz,\
                                     g_yy_xyzzzz, g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz,\
                                     g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz,\
                                     g_yy_zzzzzz, g_yz_xxxxxx, g_yz_xxxxxy, g_yz_xxxxxz,\
                                     g_yz_xxxxyy, g_yz_xxxxyz, g_yz_xxxxzz, g_yz_xxxyyy,\
                                     g_yz_xxxyyz, g_yz_xxxyzz, g_yz_xxxzzz, g_yz_xxyyyy,\
                                     g_yz_xxyyyz, g_yz_xxyyzz, g_yz_xxyzzz, g_yz_xxzzzz,\
                                     g_yz_xyyyyy, g_yz_xyyyyz, g_yz_xyyyzz, g_yz_xyyzzz,\
                                     g_yz_xyzzzz, g_yz_xzzzzz, g_yz_yyyyyy, g_yz_yyyyyz,\
                                     g_yz_yyyyzz, g_yz_yyyzzz, g_yz_yyzzzz, g_yz_yzzzzz,\
                                     g_yz_zzzzzz, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz,\
                                     g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy,\
                                     g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz, g_zz_xxyyyy,\
                                     g_zz_xxyyyz, g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz,\
                                     g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz,\
                                     g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy, g_zz_yyyyyz,\
                                     g_zz_yyyyzz, g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz,\
                                     g_zz_zzzzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xx_xxxxxx[j] = g1_x_xxxxxxx[j] - fr * g2_x_xxxxxx[j];

                g_xx_xxxxxy[j] = g1_x_xxxxxxy[j] - fr * g2_x_xxxxxy[j];

                g_xx_xxxxxz[j] = g1_x_xxxxxxz[j] - fr * g2_x_xxxxxz[j];

                g_xx_xxxxyy[j] = g1_x_xxxxxyy[j] - fr * g2_x_xxxxyy[j];

                g_xx_xxxxyz[j] = g1_x_xxxxxyz[j] - fr * g2_x_xxxxyz[j];

                g_xx_xxxxzz[j] = g1_x_xxxxxzz[j] - fr * g2_x_xxxxzz[j];

                g_xx_xxxyyy[j] = g1_x_xxxxyyy[j] - fr * g2_x_xxxyyy[j];

                g_xx_xxxyyz[j] = g1_x_xxxxyyz[j] - fr * g2_x_xxxyyz[j];

                g_xx_xxxyzz[j] = g1_x_xxxxyzz[j] - fr * g2_x_xxxyzz[j];

                g_xx_xxxzzz[j] = g1_x_xxxxzzz[j] - fr * g2_x_xxxzzz[j];

                g_xx_xxyyyy[j] = g1_x_xxxyyyy[j] - fr * g2_x_xxyyyy[j];

                g_xx_xxyyyz[j] = g1_x_xxxyyyz[j] - fr * g2_x_xxyyyz[j];

                g_xx_xxyyzz[j] = g1_x_xxxyyzz[j] - fr * g2_x_xxyyzz[j];

                g_xx_xxyzzz[j] = g1_x_xxxyzzz[j] - fr * g2_x_xxyzzz[j];

                g_xx_xxzzzz[j] = g1_x_xxxzzzz[j] - fr * g2_x_xxzzzz[j];

                g_xx_xyyyyy[j] = g1_x_xxyyyyy[j] - fr * g2_x_xyyyyy[j];

                g_xx_xyyyyz[j] = g1_x_xxyyyyz[j] - fr * g2_x_xyyyyz[j];

                g_xx_xyyyzz[j] = g1_x_xxyyyzz[j] - fr * g2_x_xyyyzz[j];

                g_xx_xyyzzz[j] = g1_x_xxyyzzz[j] - fr * g2_x_xyyzzz[j];

                g_xx_xyzzzz[j] = g1_x_xxyzzzz[j] - fr * g2_x_xyzzzz[j];

                g_xx_xzzzzz[j] = g1_x_xxzzzzz[j] - fr * g2_x_xzzzzz[j];

                g_xx_yyyyyy[j] = g1_x_xyyyyyy[j] - fr * g2_x_yyyyyy[j];

                g_xx_yyyyyz[j] = g1_x_xyyyyyz[j] - fr * g2_x_yyyyyz[j];

                g_xx_yyyyzz[j] = g1_x_xyyyyzz[j] - fr * g2_x_yyyyzz[j];

                g_xx_yyyzzz[j] = g1_x_xyyyzzz[j] - fr * g2_x_yyyzzz[j];

                g_xx_yyzzzz[j] = g1_x_xyyzzzz[j] - fr * g2_x_yyzzzz[j];

                g_xx_yzzzzz[j] = g1_x_xyzzzzz[j] - fr * g2_x_yzzzzz[j];

                g_xx_zzzzzz[j] = g1_x_xzzzzzz[j] - fr * g2_x_zzzzzz[j];

                g_xy_xxxxxx[j] = g1_y_xxxxxxx[j] - fr * g2_y_xxxxxx[j];

                g_xy_xxxxxy[j] = g1_y_xxxxxxy[j] - fr * g2_y_xxxxxy[j];

                g_xy_xxxxxz[j] = g1_y_xxxxxxz[j] - fr * g2_y_xxxxxz[j];

                g_xy_xxxxyy[j] = g1_y_xxxxxyy[j] - fr * g2_y_xxxxyy[j];

                g_xy_xxxxyz[j] = g1_y_xxxxxyz[j] - fr * g2_y_xxxxyz[j];

                g_xy_xxxxzz[j] = g1_y_xxxxxzz[j] - fr * g2_y_xxxxzz[j];

                g_xy_xxxyyy[j] = g1_y_xxxxyyy[j] - fr * g2_y_xxxyyy[j];

                g_xy_xxxyyz[j] = g1_y_xxxxyyz[j] - fr * g2_y_xxxyyz[j];

                g_xy_xxxyzz[j] = g1_y_xxxxyzz[j] - fr * g2_y_xxxyzz[j];

                g_xy_xxxzzz[j] = g1_y_xxxxzzz[j] - fr * g2_y_xxxzzz[j];

                g_xy_xxyyyy[j] = g1_y_xxxyyyy[j] - fr * g2_y_xxyyyy[j];

                g_xy_xxyyyz[j] = g1_y_xxxyyyz[j] - fr * g2_y_xxyyyz[j];

                g_xy_xxyyzz[j] = g1_y_xxxyyzz[j] - fr * g2_y_xxyyzz[j];

                g_xy_xxyzzz[j] = g1_y_xxxyzzz[j] - fr * g2_y_xxyzzz[j];

                g_xy_xxzzzz[j] = g1_y_xxxzzzz[j] - fr * g2_y_xxzzzz[j];

                g_xy_xyyyyy[j] = g1_y_xxyyyyy[j] - fr * g2_y_xyyyyy[j];

                g_xy_xyyyyz[j] = g1_y_xxyyyyz[j] - fr * g2_y_xyyyyz[j];

                g_xy_xyyyzz[j] = g1_y_xxyyyzz[j] - fr * g2_y_xyyyzz[j];

                g_xy_xyyzzz[j] = g1_y_xxyyzzz[j] - fr * g2_y_xyyzzz[j];

                g_xy_xyzzzz[j] = g1_y_xxyzzzz[j] - fr * g2_y_xyzzzz[j];

                g_xy_xzzzzz[j] = g1_y_xxzzzzz[j] - fr * g2_y_xzzzzz[j];

                g_xy_yyyyyy[j] = g1_y_xyyyyyy[j] - fr * g2_y_yyyyyy[j];

                g_xy_yyyyyz[j] = g1_y_xyyyyyz[j] - fr * g2_y_yyyyyz[j];

                g_xy_yyyyzz[j] = g1_y_xyyyyzz[j] - fr * g2_y_yyyyzz[j];

                g_xy_yyyzzz[j] = g1_y_xyyyzzz[j] - fr * g2_y_yyyzzz[j];

                g_xy_yyzzzz[j] = g1_y_xyyzzzz[j] - fr * g2_y_yyzzzz[j];

                g_xy_yzzzzz[j] = g1_y_xyzzzzz[j] - fr * g2_y_yzzzzz[j];

                g_xy_zzzzzz[j] = g1_y_xzzzzzz[j] - fr * g2_y_zzzzzz[j];

                g_xz_xxxxxx[j] = g1_z_xxxxxxx[j] - fr * g2_z_xxxxxx[j];

                g_xz_xxxxxy[j] = g1_z_xxxxxxy[j] - fr * g2_z_xxxxxy[j];

                g_xz_xxxxxz[j] = g1_z_xxxxxxz[j] - fr * g2_z_xxxxxz[j];

                g_xz_xxxxyy[j] = g1_z_xxxxxyy[j] - fr * g2_z_xxxxyy[j];

                g_xz_xxxxyz[j] = g1_z_xxxxxyz[j] - fr * g2_z_xxxxyz[j];

                g_xz_xxxxzz[j] = g1_z_xxxxxzz[j] - fr * g2_z_xxxxzz[j];

                g_xz_xxxyyy[j] = g1_z_xxxxyyy[j] - fr * g2_z_xxxyyy[j];

                g_xz_xxxyyz[j] = g1_z_xxxxyyz[j] - fr * g2_z_xxxyyz[j];

                g_xz_xxxyzz[j] = g1_z_xxxxyzz[j] - fr * g2_z_xxxyzz[j];

                g_xz_xxxzzz[j] = g1_z_xxxxzzz[j] - fr * g2_z_xxxzzz[j];

                g_xz_xxyyyy[j] = g1_z_xxxyyyy[j] - fr * g2_z_xxyyyy[j];

                g_xz_xxyyyz[j] = g1_z_xxxyyyz[j] - fr * g2_z_xxyyyz[j];

                g_xz_xxyyzz[j] = g1_z_xxxyyzz[j] - fr * g2_z_xxyyzz[j];

                g_xz_xxyzzz[j] = g1_z_xxxyzzz[j] - fr * g2_z_xxyzzz[j];

                g_xz_xxzzzz[j] = g1_z_xxxzzzz[j] - fr * g2_z_xxzzzz[j];

                g_xz_xyyyyy[j] = g1_z_xxyyyyy[j] - fr * g2_z_xyyyyy[j];

                g_xz_xyyyyz[j] = g1_z_xxyyyyz[j] - fr * g2_z_xyyyyz[j];

                g_xz_xyyyzz[j] = g1_z_xxyyyzz[j] - fr * g2_z_xyyyzz[j];

                g_xz_xyyzzz[j] = g1_z_xxyyzzz[j] - fr * g2_z_xyyzzz[j];

                g_xz_xyzzzz[j] = g1_z_xxyzzzz[j] - fr * g2_z_xyzzzz[j];

                g_xz_xzzzzz[j] = g1_z_xxzzzzz[j] - fr * g2_z_xzzzzz[j];

                g_xz_yyyyyy[j] = g1_z_xyyyyyy[j] - fr * g2_z_yyyyyy[j];

                g_xz_yyyyyz[j] = g1_z_xyyyyyz[j] - fr * g2_z_yyyyyz[j];

                g_xz_yyyyzz[j] = g1_z_xyyyyzz[j] - fr * g2_z_yyyyzz[j];

                g_xz_yyyzzz[j] = g1_z_xyyyzzz[j] - fr * g2_z_yyyzzz[j];

                g_xz_yyzzzz[j] = g1_z_xyyzzzz[j] - fr * g2_z_yyzzzz[j];

                g_xz_yzzzzz[j] = g1_z_xyzzzzz[j] - fr * g2_z_yzzzzz[j];

                g_xz_zzzzzz[j] = g1_z_xzzzzzz[j] - fr * g2_z_zzzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yy_xxxxxx[j] = g1_y_xxxxxxy[j] - fr * g2_y_xxxxxx[j];

                g_yy_xxxxxy[j] = g1_y_xxxxxyy[j] - fr * g2_y_xxxxxy[j];

                g_yy_xxxxxz[j] = g1_y_xxxxxyz[j] - fr * g2_y_xxxxxz[j];

                g_yy_xxxxyy[j] = g1_y_xxxxyyy[j] - fr * g2_y_xxxxyy[j];

                g_yy_xxxxyz[j] = g1_y_xxxxyyz[j] - fr * g2_y_xxxxyz[j];

                g_yy_xxxxzz[j] = g1_y_xxxxyzz[j] - fr * g2_y_xxxxzz[j];

                g_yy_xxxyyy[j] = g1_y_xxxyyyy[j] - fr * g2_y_xxxyyy[j];

                g_yy_xxxyyz[j] = g1_y_xxxyyyz[j] - fr * g2_y_xxxyyz[j];

                g_yy_xxxyzz[j] = g1_y_xxxyyzz[j] - fr * g2_y_xxxyzz[j];

                g_yy_xxxzzz[j] = g1_y_xxxyzzz[j] - fr * g2_y_xxxzzz[j];

                g_yy_xxyyyy[j] = g1_y_xxyyyyy[j] - fr * g2_y_xxyyyy[j];

                g_yy_xxyyyz[j] = g1_y_xxyyyyz[j] - fr * g2_y_xxyyyz[j];

                g_yy_xxyyzz[j] = g1_y_xxyyyzz[j] - fr * g2_y_xxyyzz[j];

                g_yy_xxyzzz[j] = g1_y_xxyyzzz[j] - fr * g2_y_xxyzzz[j];

                g_yy_xxzzzz[j] = g1_y_xxyzzzz[j] - fr * g2_y_xxzzzz[j];

                g_yy_xyyyyy[j] = g1_y_xyyyyyy[j] - fr * g2_y_xyyyyy[j];

                g_yy_xyyyyz[j] = g1_y_xyyyyyz[j] - fr * g2_y_xyyyyz[j];

                g_yy_xyyyzz[j] = g1_y_xyyyyzz[j] - fr * g2_y_xyyyzz[j];

                g_yy_xyyzzz[j] = g1_y_xyyyzzz[j] - fr * g2_y_xyyzzz[j];

                g_yy_xyzzzz[j] = g1_y_xyyzzzz[j] - fr * g2_y_xyzzzz[j];

                g_yy_xzzzzz[j] = g1_y_xyzzzzz[j] - fr * g2_y_xzzzzz[j];

                g_yy_yyyyyy[j] = g1_y_yyyyyyy[j] - fr * g2_y_yyyyyy[j];

                g_yy_yyyyyz[j] = g1_y_yyyyyyz[j] - fr * g2_y_yyyyyz[j];

                g_yy_yyyyzz[j] = g1_y_yyyyyzz[j] - fr * g2_y_yyyyzz[j];

                g_yy_yyyzzz[j] = g1_y_yyyyzzz[j] - fr * g2_y_yyyzzz[j];

                g_yy_yyzzzz[j] = g1_y_yyyzzzz[j] - fr * g2_y_yyzzzz[j];

                g_yy_yzzzzz[j] = g1_y_yyzzzzz[j] - fr * g2_y_yzzzzz[j];

                g_yy_zzzzzz[j] = g1_y_yzzzzzz[j] - fr * g2_y_zzzzzz[j];

                g_yz_xxxxxx[j] = g1_z_xxxxxxy[j] - fr * g2_z_xxxxxx[j];

                g_yz_xxxxxy[j] = g1_z_xxxxxyy[j] - fr * g2_z_xxxxxy[j];

                g_yz_xxxxxz[j] = g1_z_xxxxxyz[j] - fr * g2_z_xxxxxz[j];

                g_yz_xxxxyy[j] = g1_z_xxxxyyy[j] - fr * g2_z_xxxxyy[j];

                g_yz_xxxxyz[j] = g1_z_xxxxyyz[j] - fr * g2_z_xxxxyz[j];

                g_yz_xxxxzz[j] = g1_z_xxxxyzz[j] - fr * g2_z_xxxxzz[j];

                g_yz_xxxyyy[j] = g1_z_xxxyyyy[j] - fr * g2_z_xxxyyy[j];

                g_yz_xxxyyz[j] = g1_z_xxxyyyz[j] - fr * g2_z_xxxyyz[j];

                g_yz_xxxyzz[j] = g1_z_xxxyyzz[j] - fr * g2_z_xxxyzz[j];

                g_yz_xxxzzz[j] = g1_z_xxxyzzz[j] - fr * g2_z_xxxzzz[j];

                g_yz_xxyyyy[j] = g1_z_xxyyyyy[j] - fr * g2_z_xxyyyy[j];

                g_yz_xxyyyz[j] = g1_z_xxyyyyz[j] - fr * g2_z_xxyyyz[j];

                g_yz_xxyyzz[j] = g1_z_xxyyyzz[j] - fr * g2_z_xxyyzz[j];

                g_yz_xxyzzz[j] = g1_z_xxyyzzz[j] - fr * g2_z_xxyzzz[j];

                g_yz_xxzzzz[j] = g1_z_xxyzzzz[j] - fr * g2_z_xxzzzz[j];

                g_yz_xyyyyy[j] = g1_z_xyyyyyy[j] - fr * g2_z_xyyyyy[j];

                g_yz_xyyyyz[j] = g1_z_xyyyyyz[j] - fr * g2_z_xyyyyz[j];

                g_yz_xyyyzz[j] = g1_z_xyyyyzz[j] - fr * g2_z_xyyyzz[j];

                g_yz_xyyzzz[j] = g1_z_xyyyzzz[j] - fr * g2_z_xyyzzz[j];

                g_yz_xyzzzz[j] = g1_z_xyyzzzz[j] - fr * g2_z_xyzzzz[j];

                g_yz_xzzzzz[j] = g1_z_xyzzzzz[j] - fr * g2_z_xzzzzz[j];

                g_yz_yyyyyy[j] = g1_z_yyyyyyy[j] - fr * g2_z_yyyyyy[j];

                g_yz_yyyyyz[j] = g1_z_yyyyyyz[j] - fr * g2_z_yyyyyz[j];

                g_yz_yyyyzz[j] = g1_z_yyyyyzz[j] - fr * g2_z_yyyyzz[j];

                g_yz_yyyzzz[j] = g1_z_yyyyzzz[j] - fr * g2_z_yyyzzz[j];

                g_yz_yyzzzz[j] = g1_z_yyyzzzz[j] - fr * g2_z_yyzzzz[j];

                g_yz_yzzzzz[j] = g1_z_yyzzzzz[j] - fr * g2_z_yzzzzz[j];

                g_yz_zzzzzz[j] = g1_z_yzzzzzz[j] - fr * g2_z_zzzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zz_xxxxxx[j] = g1_z_xxxxxxz[j] - fr * g2_z_xxxxxx[j];

                g_zz_xxxxxy[j] = g1_z_xxxxxyz[j] - fr * g2_z_xxxxxy[j];

                g_zz_xxxxxz[j] = g1_z_xxxxxzz[j] - fr * g2_z_xxxxxz[j];

                g_zz_xxxxyy[j] = g1_z_xxxxyyz[j] - fr * g2_z_xxxxyy[j];

                g_zz_xxxxyz[j] = g1_z_xxxxyzz[j] - fr * g2_z_xxxxyz[j];

                g_zz_xxxxzz[j] = g1_z_xxxxzzz[j] - fr * g2_z_xxxxzz[j];

                g_zz_xxxyyy[j] = g1_z_xxxyyyz[j] - fr * g2_z_xxxyyy[j];

                g_zz_xxxyyz[j] = g1_z_xxxyyzz[j] - fr * g2_z_xxxyyz[j];

                g_zz_xxxyzz[j] = g1_z_xxxyzzz[j] - fr * g2_z_xxxyzz[j];

                g_zz_xxxzzz[j] = g1_z_xxxzzzz[j] - fr * g2_z_xxxzzz[j];

                g_zz_xxyyyy[j] = g1_z_xxyyyyz[j] - fr * g2_z_xxyyyy[j];

                g_zz_xxyyyz[j] = g1_z_xxyyyzz[j] - fr * g2_z_xxyyyz[j];

                g_zz_xxyyzz[j] = g1_z_xxyyzzz[j] - fr * g2_z_xxyyzz[j];

                g_zz_xxyzzz[j] = g1_z_xxyzzzz[j] - fr * g2_z_xxyzzz[j];

                g_zz_xxzzzz[j] = g1_z_xxzzzzz[j] - fr * g2_z_xxzzzz[j];

                g_zz_xyyyyy[j] = g1_z_xyyyyyz[j] - fr * g2_z_xyyyyy[j];

                g_zz_xyyyyz[j] = g1_z_xyyyyzz[j] - fr * g2_z_xyyyyz[j];

                g_zz_xyyyzz[j] = g1_z_xyyyzzz[j] - fr * g2_z_xyyyzz[j];

                g_zz_xyyzzz[j] = g1_z_xyyzzzz[j] - fr * g2_z_xyyzzz[j];

                g_zz_xyzzzz[j] = g1_z_xyzzzzz[j] - fr * g2_z_xyzzzz[j];

                g_zz_xzzzzz[j] = g1_z_xzzzzzz[j] - fr * g2_z_xzzzzz[j];

                g_zz_yyyyyy[j] = g1_z_yyyyyyz[j] - fr * g2_z_yyyyyy[j];

                g_zz_yyyyyz[j] = g1_z_yyyyyzz[j] - fr * g2_z_yyyyyz[j];

                g_zz_yyyyzz[j] = g1_z_yyyyzzz[j] - fr * g2_z_yyyyzz[j];

                g_zz_yyyzzz[j] = g1_z_yyyzzzz[j] - fr * g2_z_yyyzzz[j];

                g_zz_yyzzzz[j] = g1_z_yyzzzzz[j] - fr * g2_z_yyzzzz[j];

                g_zz_yzzzzz[j] = g1_z_yzzzzzz[j] - fr * g2_z_yzzzzz[j];

                g_zz_zzzzzz[j] = g1_z_zzzzzzz[j] - fr * g2_z_zzzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXFF(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 3, 3})) return;

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
                                             {braAngularMomentum, 3, 3});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 4});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 3});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|DF)^(m) integrals

            auto g2_xx_xxx = contrBuffer.data(g2off + 60 * i);

            auto g2_xx_xxy = contrBuffer.data(g2off + 60 * i + 1);

            auto g2_xx_xxz = contrBuffer.data(g2off + 60 * i + 2);

            auto g2_xx_xyy = contrBuffer.data(g2off + 60 * i + 3);

            auto g2_xx_xyz = contrBuffer.data(g2off + 60 * i + 4);

            auto g2_xx_xzz = contrBuffer.data(g2off + 60 * i + 5);

            auto g2_xx_yyy = contrBuffer.data(g2off + 60 * i + 6);

            auto g2_xx_yyz = contrBuffer.data(g2off + 60 * i + 7);

            auto g2_xx_yzz = contrBuffer.data(g2off + 60 * i + 8);

            auto g2_xx_zzz = contrBuffer.data(g2off + 60 * i + 9);

            auto g2_xy_xxx = contrBuffer.data(g2off + 60 * i + 10);

            auto g2_xy_xxy = contrBuffer.data(g2off + 60 * i + 11);

            auto g2_xy_xxz = contrBuffer.data(g2off + 60 * i + 12);

            auto g2_xy_xyy = contrBuffer.data(g2off + 60 * i + 13);

            auto g2_xy_xyz = contrBuffer.data(g2off + 60 * i + 14);

            auto g2_xy_xzz = contrBuffer.data(g2off + 60 * i + 15);

            auto g2_xy_yyy = contrBuffer.data(g2off + 60 * i + 16);

            auto g2_xy_yyz = contrBuffer.data(g2off + 60 * i + 17);

            auto g2_xy_yzz = contrBuffer.data(g2off + 60 * i + 18);

            auto g2_xy_zzz = contrBuffer.data(g2off + 60 * i + 19);

            auto g2_xz_xxx = contrBuffer.data(g2off + 60 * i + 20);

            auto g2_xz_xxy = contrBuffer.data(g2off + 60 * i + 21);

            auto g2_xz_xxz = contrBuffer.data(g2off + 60 * i + 22);

            auto g2_xz_xyy = contrBuffer.data(g2off + 60 * i + 23);

            auto g2_xz_xyz = contrBuffer.data(g2off + 60 * i + 24);

            auto g2_xz_xzz = contrBuffer.data(g2off + 60 * i + 25);

            auto g2_xz_yyy = contrBuffer.data(g2off + 60 * i + 26);

            auto g2_xz_yyz = contrBuffer.data(g2off + 60 * i + 27);

            auto g2_xz_yzz = contrBuffer.data(g2off + 60 * i + 28);

            auto g2_xz_zzz = contrBuffer.data(g2off + 60 * i + 29);

            auto g2_yy_xxx = contrBuffer.data(g2off + 60 * i + 30);

            auto g2_yy_xxy = contrBuffer.data(g2off + 60 * i + 31);

            auto g2_yy_xxz = contrBuffer.data(g2off + 60 * i + 32);

            auto g2_yy_xyy = contrBuffer.data(g2off + 60 * i + 33);

            auto g2_yy_xyz = contrBuffer.data(g2off + 60 * i + 34);

            auto g2_yy_xzz = contrBuffer.data(g2off + 60 * i + 35);

            auto g2_yy_yyy = contrBuffer.data(g2off + 60 * i + 36);

            auto g2_yy_yyz = contrBuffer.data(g2off + 60 * i + 37);

            auto g2_yy_yzz = contrBuffer.data(g2off + 60 * i + 38);

            auto g2_yy_zzz = contrBuffer.data(g2off + 60 * i + 39);

            auto g2_yz_xxx = contrBuffer.data(g2off + 60 * i + 40);

            auto g2_yz_xxy = contrBuffer.data(g2off + 60 * i + 41);

            auto g2_yz_xxz = contrBuffer.data(g2off + 60 * i + 42);

            auto g2_yz_xyy = contrBuffer.data(g2off + 60 * i + 43);

            auto g2_yz_xyz = contrBuffer.data(g2off + 60 * i + 44);

            auto g2_yz_xzz = contrBuffer.data(g2off + 60 * i + 45);

            auto g2_yz_yyy = contrBuffer.data(g2off + 60 * i + 46);

            auto g2_yz_yyz = contrBuffer.data(g2off + 60 * i + 47);

            auto g2_yz_yzz = contrBuffer.data(g2off + 60 * i + 48);

            auto g2_yz_zzz = contrBuffer.data(g2off + 60 * i + 49);

            auto g2_zz_xxx = contrBuffer.data(g2off + 60 * i + 50);

            auto g2_zz_xxy = contrBuffer.data(g2off + 60 * i + 51);

            auto g2_zz_xxz = contrBuffer.data(g2off + 60 * i + 52);

            auto g2_zz_xyy = contrBuffer.data(g2off + 60 * i + 53);

            auto g2_zz_xyz = contrBuffer.data(g2off + 60 * i + 54);

            auto g2_zz_xzz = contrBuffer.data(g2off + 60 * i + 55);

            auto g2_zz_yyy = contrBuffer.data(g2off + 60 * i + 56);

            auto g2_zz_yyz = contrBuffer.data(g2off + 60 * i + 57);

            auto g2_zz_yzz = contrBuffer.data(g2off + 60 * i + 58);

            auto g2_zz_zzz = contrBuffer.data(g2off + 60 * i + 59);

            // set up pointers to (X|g(r,r')|DG)^(m) integrals

            auto g1_xx_xxxx = contrBuffer.data(g1off + 90 * i);

            auto g1_xx_xxxy = contrBuffer.data(g1off + 90 * i + 1);

            auto g1_xx_xxxz = contrBuffer.data(g1off + 90 * i + 2);

            auto g1_xx_xxyy = contrBuffer.data(g1off + 90 * i + 3);

            auto g1_xx_xxyz = contrBuffer.data(g1off + 90 * i + 4);

            auto g1_xx_xxzz = contrBuffer.data(g1off + 90 * i + 5);

            auto g1_xx_xyyy = contrBuffer.data(g1off + 90 * i + 6);

            auto g1_xx_xyyz = contrBuffer.data(g1off + 90 * i + 7);

            auto g1_xx_xyzz = contrBuffer.data(g1off + 90 * i + 8);

            auto g1_xx_xzzz = contrBuffer.data(g1off + 90 * i + 9);

            auto g1_xy_xxxx = contrBuffer.data(g1off + 90 * i + 15);

            auto g1_xy_xxxy = contrBuffer.data(g1off + 90 * i + 16);

            auto g1_xy_xxxz = contrBuffer.data(g1off + 90 * i + 17);

            auto g1_xy_xxyy = contrBuffer.data(g1off + 90 * i + 18);

            auto g1_xy_xxyz = contrBuffer.data(g1off + 90 * i + 19);

            auto g1_xy_xxzz = contrBuffer.data(g1off + 90 * i + 20);

            auto g1_xy_xyyy = contrBuffer.data(g1off + 90 * i + 21);

            auto g1_xy_xyyz = contrBuffer.data(g1off + 90 * i + 22);

            auto g1_xy_xyzz = contrBuffer.data(g1off + 90 * i + 23);

            auto g1_xy_xzzz = contrBuffer.data(g1off + 90 * i + 24);

            auto g1_xz_xxxx = contrBuffer.data(g1off + 90 * i + 30);

            auto g1_xz_xxxy = contrBuffer.data(g1off + 90 * i + 31);

            auto g1_xz_xxxz = contrBuffer.data(g1off + 90 * i + 32);

            auto g1_xz_xxyy = contrBuffer.data(g1off + 90 * i + 33);

            auto g1_xz_xxyz = contrBuffer.data(g1off + 90 * i + 34);

            auto g1_xz_xxzz = contrBuffer.data(g1off + 90 * i + 35);

            auto g1_xz_xyyy = contrBuffer.data(g1off + 90 * i + 36);

            auto g1_xz_xyyz = contrBuffer.data(g1off + 90 * i + 37);

            auto g1_xz_xyzz = contrBuffer.data(g1off + 90 * i + 38);

            auto g1_xz_xzzz = contrBuffer.data(g1off + 90 * i + 39);

            auto g1_yy_xxxx = contrBuffer.data(g1off + 90 * i + 45);

            auto g1_yy_xxxy = contrBuffer.data(g1off + 90 * i + 46);

            auto g1_yy_xxxz = contrBuffer.data(g1off + 90 * i + 47);

            auto g1_yy_xxyy = contrBuffer.data(g1off + 90 * i + 48);

            auto g1_yy_xxyz = contrBuffer.data(g1off + 90 * i + 49);

            auto g1_yy_xxzz = contrBuffer.data(g1off + 90 * i + 50);

            auto g1_yy_xyyy = contrBuffer.data(g1off + 90 * i + 51);

            auto g1_yy_xyyz = contrBuffer.data(g1off + 90 * i + 52);

            auto g1_yy_xyzz = contrBuffer.data(g1off + 90 * i + 53);

            auto g1_yy_xzzz = contrBuffer.data(g1off + 90 * i + 54);

            auto g1_yy_yyyy = contrBuffer.data(g1off + 90 * i + 55);

            auto g1_yy_yyyz = contrBuffer.data(g1off + 90 * i + 56);

            auto g1_yy_yyzz = contrBuffer.data(g1off + 90 * i + 57);

            auto g1_yy_yzzz = contrBuffer.data(g1off + 90 * i + 58);

            auto g1_yz_xxxx = contrBuffer.data(g1off + 90 * i + 60);

            auto g1_yz_xxxy = contrBuffer.data(g1off + 90 * i + 61);

            auto g1_yz_xxxz = contrBuffer.data(g1off + 90 * i + 62);

            auto g1_yz_xxyy = contrBuffer.data(g1off + 90 * i + 63);

            auto g1_yz_xxyz = contrBuffer.data(g1off + 90 * i + 64);

            auto g1_yz_xxzz = contrBuffer.data(g1off + 90 * i + 65);

            auto g1_yz_xyyy = contrBuffer.data(g1off + 90 * i + 66);

            auto g1_yz_xyyz = contrBuffer.data(g1off + 90 * i + 67);

            auto g1_yz_xyzz = contrBuffer.data(g1off + 90 * i + 68);

            auto g1_yz_xzzz = contrBuffer.data(g1off + 90 * i + 69);

            auto g1_yz_yyyy = contrBuffer.data(g1off + 90 * i + 70);

            auto g1_yz_yyyz = contrBuffer.data(g1off + 90 * i + 71);

            auto g1_yz_yyzz = contrBuffer.data(g1off + 90 * i + 72);

            auto g1_yz_yzzz = contrBuffer.data(g1off + 90 * i + 73);

            auto g1_zz_xxxx = contrBuffer.data(g1off + 90 * i + 75);

            auto g1_zz_xxxy = contrBuffer.data(g1off + 90 * i + 76);

            auto g1_zz_xxxz = contrBuffer.data(g1off + 90 * i + 77);

            auto g1_zz_xxyy = contrBuffer.data(g1off + 90 * i + 78);

            auto g1_zz_xxyz = contrBuffer.data(g1off + 90 * i + 79);

            auto g1_zz_xxzz = contrBuffer.data(g1off + 90 * i + 80);

            auto g1_zz_xyyy = contrBuffer.data(g1off + 90 * i + 81);

            auto g1_zz_xyyz = contrBuffer.data(g1off + 90 * i + 82);

            auto g1_zz_xyzz = contrBuffer.data(g1off + 90 * i + 83);

            auto g1_zz_xzzz = contrBuffer.data(g1off + 90 * i + 84);

            auto g1_zz_yyyy = contrBuffer.data(g1off + 90 * i + 85);

            auto g1_zz_yyyz = contrBuffer.data(g1off + 90 * i + 86);

            auto g1_zz_yyzz = contrBuffer.data(g1off + 90 * i + 87);

            auto g1_zz_yzzz = contrBuffer.data(g1off + 90 * i + 88);

            auto g1_zz_zzzz = contrBuffer.data(g1off + 90 * i + 89);

            // set up pointers to (X|g(r,r')|FF)^(m) integrals

            auto g_xxx_xxx = contrBuffer.data(goff + 100 * i);

            auto g_xxx_xxy = contrBuffer.data(goff + 100 * i + 1);

            auto g_xxx_xxz = contrBuffer.data(goff + 100 * i + 2);

            auto g_xxx_xyy = contrBuffer.data(goff + 100 * i + 3);

            auto g_xxx_xyz = contrBuffer.data(goff + 100 * i + 4);

            auto g_xxx_xzz = contrBuffer.data(goff + 100 * i + 5);

            auto g_xxx_yyy = contrBuffer.data(goff + 100 * i + 6);

            auto g_xxx_yyz = contrBuffer.data(goff + 100 * i + 7);

            auto g_xxx_yzz = contrBuffer.data(goff + 100 * i + 8);

            auto g_xxx_zzz = contrBuffer.data(goff + 100 * i + 9);

            auto g_xxy_xxx = contrBuffer.data(goff + 100 * i + 10);

            auto g_xxy_xxy = contrBuffer.data(goff + 100 * i + 11);

            auto g_xxy_xxz = contrBuffer.data(goff + 100 * i + 12);

            auto g_xxy_xyy = contrBuffer.data(goff + 100 * i + 13);

            auto g_xxy_xyz = contrBuffer.data(goff + 100 * i + 14);

            auto g_xxy_xzz = contrBuffer.data(goff + 100 * i + 15);

            auto g_xxy_yyy = contrBuffer.data(goff + 100 * i + 16);

            auto g_xxy_yyz = contrBuffer.data(goff + 100 * i + 17);

            auto g_xxy_yzz = contrBuffer.data(goff + 100 * i + 18);

            auto g_xxy_zzz = contrBuffer.data(goff + 100 * i + 19);

            auto g_xxz_xxx = contrBuffer.data(goff + 100 * i + 20);

            auto g_xxz_xxy = contrBuffer.data(goff + 100 * i + 21);

            auto g_xxz_xxz = contrBuffer.data(goff + 100 * i + 22);

            auto g_xxz_xyy = contrBuffer.data(goff + 100 * i + 23);

            auto g_xxz_xyz = contrBuffer.data(goff + 100 * i + 24);

            auto g_xxz_xzz = contrBuffer.data(goff + 100 * i + 25);

            auto g_xxz_yyy = contrBuffer.data(goff + 100 * i + 26);

            auto g_xxz_yyz = contrBuffer.data(goff + 100 * i + 27);

            auto g_xxz_yzz = contrBuffer.data(goff + 100 * i + 28);

            auto g_xxz_zzz = contrBuffer.data(goff + 100 * i + 29);

            auto g_xyy_xxx = contrBuffer.data(goff + 100 * i + 30);

            auto g_xyy_xxy = contrBuffer.data(goff + 100 * i + 31);

            auto g_xyy_xxz = contrBuffer.data(goff + 100 * i + 32);

            auto g_xyy_xyy = contrBuffer.data(goff + 100 * i + 33);

            auto g_xyy_xyz = contrBuffer.data(goff + 100 * i + 34);

            auto g_xyy_xzz = contrBuffer.data(goff + 100 * i + 35);

            auto g_xyy_yyy = contrBuffer.data(goff + 100 * i + 36);

            auto g_xyy_yyz = contrBuffer.data(goff + 100 * i + 37);

            auto g_xyy_yzz = contrBuffer.data(goff + 100 * i + 38);

            auto g_xyy_zzz = contrBuffer.data(goff + 100 * i + 39);

            auto g_xyz_xxx = contrBuffer.data(goff + 100 * i + 40);

            auto g_xyz_xxy = contrBuffer.data(goff + 100 * i + 41);

            auto g_xyz_xxz = contrBuffer.data(goff + 100 * i + 42);

            auto g_xyz_xyy = contrBuffer.data(goff + 100 * i + 43);

            auto g_xyz_xyz = contrBuffer.data(goff + 100 * i + 44);

            auto g_xyz_xzz = contrBuffer.data(goff + 100 * i + 45);

            auto g_xyz_yyy = contrBuffer.data(goff + 100 * i + 46);

            auto g_xyz_yyz = contrBuffer.data(goff + 100 * i + 47);

            auto g_xyz_yzz = contrBuffer.data(goff + 100 * i + 48);

            auto g_xyz_zzz = contrBuffer.data(goff + 100 * i + 49);

            auto g_xzz_xxx = contrBuffer.data(goff + 100 * i + 50);

            auto g_xzz_xxy = contrBuffer.data(goff + 100 * i + 51);

            auto g_xzz_xxz = contrBuffer.data(goff + 100 * i + 52);

            auto g_xzz_xyy = contrBuffer.data(goff + 100 * i + 53);

            auto g_xzz_xyz = contrBuffer.data(goff + 100 * i + 54);

            auto g_xzz_xzz = contrBuffer.data(goff + 100 * i + 55);

            auto g_xzz_yyy = contrBuffer.data(goff + 100 * i + 56);

            auto g_xzz_yyz = contrBuffer.data(goff + 100 * i + 57);

            auto g_xzz_yzz = contrBuffer.data(goff + 100 * i + 58);

            auto g_xzz_zzz = contrBuffer.data(goff + 100 * i + 59);

            auto g_yyy_xxx = contrBuffer.data(goff + 100 * i + 60);

            auto g_yyy_xxy = contrBuffer.data(goff + 100 * i + 61);

            auto g_yyy_xxz = contrBuffer.data(goff + 100 * i + 62);

            auto g_yyy_xyy = contrBuffer.data(goff + 100 * i + 63);

            auto g_yyy_xyz = contrBuffer.data(goff + 100 * i + 64);

            auto g_yyy_xzz = contrBuffer.data(goff + 100 * i + 65);

            auto g_yyy_yyy = contrBuffer.data(goff + 100 * i + 66);

            auto g_yyy_yyz = contrBuffer.data(goff + 100 * i + 67);

            auto g_yyy_yzz = contrBuffer.data(goff + 100 * i + 68);

            auto g_yyy_zzz = contrBuffer.data(goff + 100 * i + 69);

            auto g_yyz_xxx = contrBuffer.data(goff + 100 * i + 70);

            auto g_yyz_xxy = contrBuffer.data(goff + 100 * i + 71);

            auto g_yyz_xxz = contrBuffer.data(goff + 100 * i + 72);

            auto g_yyz_xyy = contrBuffer.data(goff + 100 * i + 73);

            auto g_yyz_xyz = contrBuffer.data(goff + 100 * i + 74);

            auto g_yyz_xzz = contrBuffer.data(goff + 100 * i + 75);

            auto g_yyz_yyy = contrBuffer.data(goff + 100 * i + 76);

            auto g_yyz_yyz = contrBuffer.data(goff + 100 * i + 77);

            auto g_yyz_yzz = contrBuffer.data(goff + 100 * i + 78);

            auto g_yyz_zzz = contrBuffer.data(goff + 100 * i + 79);

            auto g_yzz_xxx = contrBuffer.data(goff + 100 * i + 80);

            auto g_yzz_xxy = contrBuffer.data(goff + 100 * i + 81);

            auto g_yzz_xxz = contrBuffer.data(goff + 100 * i + 82);

            auto g_yzz_xyy = contrBuffer.data(goff + 100 * i + 83);

            auto g_yzz_xyz = contrBuffer.data(goff + 100 * i + 84);

            auto g_yzz_xzz = contrBuffer.data(goff + 100 * i + 85);

            auto g_yzz_yyy = contrBuffer.data(goff + 100 * i + 86);

            auto g_yzz_yyz = contrBuffer.data(goff + 100 * i + 87);

            auto g_yzz_yzz = contrBuffer.data(goff + 100 * i + 88);

            auto g_yzz_zzz = contrBuffer.data(goff + 100 * i + 89);

            auto g_zzz_xxx = contrBuffer.data(goff + 100 * i + 90);

            auto g_zzz_xxy = contrBuffer.data(goff + 100 * i + 91);

            auto g_zzz_xxz = contrBuffer.data(goff + 100 * i + 92);

            auto g_zzz_xyy = contrBuffer.data(goff + 100 * i + 93);

            auto g_zzz_xyz = contrBuffer.data(goff + 100 * i + 94);

            auto g_zzz_xzz = contrBuffer.data(goff + 100 * i + 95);

            auto g_zzz_yyy = contrBuffer.data(goff + 100 * i + 96);

            auto g_zzz_yyz = contrBuffer.data(goff + 100 * i + 97);

            auto g_zzz_yzz = contrBuffer.data(goff + 100 * i + 98);

            auto g_zzz_zzz = contrBuffer.data(goff + 100 * i + 99);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxx, g2_xx_xxy, g2_xx_xxz,\
                                     g2_xx_xyy, g2_xx_xyz, g2_xx_xzz, g2_xx_yyy,\
                                     g2_xx_yyz, g2_xx_yzz, g2_xx_zzz, g2_xy_xxx,\
                                     g2_xy_xxy, g2_xy_xxz, g2_xy_xyy, g2_xy_xyz,\
                                     g2_xy_xzz, g2_xy_yyy, g2_xy_yyz, g2_xy_yzz,\
                                     g2_xy_zzz, g2_xz_xxx, g2_xz_xxy, g2_xz_xxz,\
                                     g2_xz_xyy, g2_xz_xyz, g2_xz_xzz, g2_xz_yyy,\
                                     g2_xz_yyz, g2_xz_yzz, g2_xz_zzz, g2_yy_xxx,\
                                     g2_yy_xxy, g2_yy_xxz, g2_yy_xyy, g2_yy_xyz,\
                                     g2_yy_xzz, g2_yy_yyy, g2_yy_yyz, g2_yy_yzz,\
                                     g2_yy_zzz, g2_yz_xxx, g2_yz_xxy, g2_yz_xxz,\
                                     g2_yz_xyy, g2_yz_xyz, g2_yz_xzz, g2_yz_yyy,\
                                     g2_yz_yyz, g2_yz_yzz, g2_yz_zzz, g2_zz_xxx,\
                                     g2_zz_xxy, g2_zz_xxz, g2_zz_xyy, g2_zz_xyz,\
                                     g2_zz_xzz, g2_zz_yyy, g2_zz_yyz, g2_zz_yzz,\
                                     g2_zz_zzz, g1_xx_xxxx, g1_xx_xxxy, g1_xx_xxxz,\
                                     g1_xx_xxyy, g1_xx_xxyz, g1_xx_xxzz, g1_xx_xyyy,\
                                     g1_xx_xyyz, g1_xx_xyzz, g1_xx_xzzz, g1_xy_xxxx,\
                                     g1_xy_xxxy, g1_xy_xxxz, g1_xy_xxyy, g1_xy_xxyz,\
                                     g1_xy_xxzz, g1_xy_xyyy, g1_xy_xyyz, g1_xy_xyzz,\
                                     g1_xy_xzzz, g1_xz_xxxx, g1_xz_xxxy, g1_xz_xxxz,\
                                     g1_xz_xxyy, g1_xz_xxyz, g1_xz_xxzz, g1_xz_xyyy,\
                                     g1_xz_xyyz, g1_xz_xyzz, g1_xz_xzzz,  g1_yy_xxxx,\
                                     g1_yy_xxxy, g1_yy_xxxz, g1_yy_xxyy, g1_yy_xxyz,\
                                     g1_yy_xxzz, g1_yy_xyyy, g1_yy_xyyz, g1_yy_xyzz,\
                                     g1_yy_xzzz, g1_yy_yyyy, g1_yy_yyyz, g1_yy_yyzz,\
                                     g1_yy_yzzz, g1_yz_xxxx, g1_yz_xxxy, g1_yz_xxxz,\
                                     g1_yz_xxyy, g1_yz_xxyz, g1_yz_xxzz, g1_yz_xyyy,\
                                     g1_yz_xyyz, g1_yz_xyzz, g1_yz_xzzz, g1_yz_yyyy,\
                                     g1_yz_yyyz, g1_yz_yyzz, g1_yz_yzzz,\
                                     g1_zz_xxxx, g1_zz_xxxy, g1_zz_xxxz, g1_zz_xxyy,\
                                     g1_zz_xxyz, g1_zz_xxzz, g1_zz_xyyy, g1_zz_xyyz,\
                                     g1_zz_xyzz, g1_zz_xzzz, g1_zz_yyyy, g1_zz_yyyz,\
                                     g1_zz_yyzz, g1_zz_yzzz, g1_zz_zzzz, g_xxx_xxx,\
                                     g_xxx_xxy, g_xxx_xxz, g_xxx_xyy, g_xxx_xyz,\
                                     g_xxx_xzz, g_xxx_yyy, g_xxx_yyz, g_xxx_yzz,\
                                     g_xxx_zzz, g_xxy_xxx, g_xxy_xxy, g_xxy_xxz,\
                                     g_xxy_xyy, g_xxy_xyz, g_xxy_xzz, g_xxy_yyy,\
                                     g_xxy_yyz, g_xxy_yzz, g_xxy_zzz, g_xxz_xxx,\
                                     g_xxz_xxy, g_xxz_xxz, g_xxz_xyy, g_xxz_xyz,\
                                     g_xxz_xzz, g_xxz_yyy, g_xxz_yyz, g_xxz_yzz,\
                                     g_xxz_zzz, g_xyy_xxx, g_xyy_xxy, g_xyy_xxz,\
                                     g_xyy_xyy, g_xyy_xyz, g_xyy_xzz, g_xyy_yyy,\
                                     g_xyy_yyz, g_xyy_yzz, g_xyy_zzz, g_xyz_xxx,\
                                     g_xyz_xxy, g_xyz_xxz, g_xyz_xyy, g_xyz_xyz,\
                                     g_xyz_xzz, g_xyz_yyy, g_xyz_yyz, g_xyz_yzz,\
                                     g_xyz_zzz, g_xzz_xxx, g_xzz_xxy, g_xzz_xxz,\
                                     g_xzz_xyy, g_xzz_xyz, g_xzz_xzz, g_xzz_yyy,\
                                     g_xzz_yyz, g_xzz_yzz, g_xzz_zzz, g_yyy_xxx,\
                                     g_yyy_xxy, g_yyy_xxz, g_yyy_xyy, g_yyy_xyz,\
                                     g_yyy_xzz, g_yyy_yyy, g_yyy_yyz, g_yyy_yzz,\
                                     g_yyy_zzz, g_yyz_xxx, g_yyz_xxy, g_yyz_xxz,\
                                     g_yyz_xyy, g_yyz_xyz, g_yyz_xzz, g_yyz_yyy,\
                                     g_yyz_yyz, g_yyz_yzz, g_yyz_zzz, g_yzz_xxx,\
                                     g_yzz_xxy, g_yzz_xxz, g_yzz_xyy, g_yzz_xyz,\
                                     g_yzz_xzz, g_yzz_yyy, g_yzz_yyz, g_yzz_yzz,\
                                     g_yzz_zzz, g_zzz_xxx, g_zzz_xxy, g_zzz_xxz,\
                                     g_zzz_xyy, g_zzz_xyz, g_zzz_xzz, g_zzz_yyy,\
                                     g_zzz_yyz, g_zzz_yzz, g_zzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xxx_xxx[j] = g1_xx_xxxx[j] - fr * g2_xx_xxx[j];

                g_xxx_xxy[j] = g1_xx_xxxy[j] - fr * g2_xx_xxy[j];

                g_xxx_xxz[j] = g1_xx_xxxz[j] - fr * g2_xx_xxz[j];

                g_xxx_xyy[j] = g1_xx_xxyy[j] - fr * g2_xx_xyy[j];

                g_xxx_xyz[j] = g1_xx_xxyz[j] - fr * g2_xx_xyz[j];

                g_xxx_xzz[j] = g1_xx_xxzz[j] - fr * g2_xx_xzz[j];

                g_xxx_yyy[j] = g1_xx_xyyy[j] - fr * g2_xx_yyy[j];

                g_xxx_yyz[j] = g1_xx_xyyz[j] - fr * g2_xx_yyz[j];

                g_xxx_yzz[j] = g1_xx_xyzz[j] - fr * g2_xx_yzz[j];

                g_xxx_zzz[j] = g1_xx_xzzz[j] - fr * g2_xx_zzz[j];

                g_xxy_xxx[j] = g1_xy_xxxx[j] - fr * g2_xy_xxx[j];

                g_xxy_xxy[j] = g1_xy_xxxy[j] - fr * g2_xy_xxy[j];

                g_xxy_xxz[j] = g1_xy_xxxz[j] - fr * g2_xy_xxz[j];

                g_xxy_xyy[j] = g1_xy_xxyy[j] - fr * g2_xy_xyy[j];

                g_xxy_xyz[j] = g1_xy_xxyz[j] - fr * g2_xy_xyz[j];

                g_xxy_xzz[j] = g1_xy_xxzz[j] - fr * g2_xy_xzz[j];

                g_xxy_yyy[j] = g1_xy_xyyy[j] - fr * g2_xy_yyy[j];

                g_xxy_yyz[j] = g1_xy_xyyz[j] - fr * g2_xy_yyz[j];

                g_xxy_yzz[j] = g1_xy_xyzz[j] - fr * g2_xy_yzz[j];

                g_xxy_zzz[j] = g1_xy_xzzz[j] - fr * g2_xy_zzz[j];

                g_xxz_xxx[j] = g1_xz_xxxx[j] - fr * g2_xz_xxx[j];

                g_xxz_xxy[j] = g1_xz_xxxy[j] - fr * g2_xz_xxy[j];

                g_xxz_xxz[j] = g1_xz_xxxz[j] - fr * g2_xz_xxz[j];

                g_xxz_xyy[j] = g1_xz_xxyy[j] - fr * g2_xz_xyy[j];

                g_xxz_xyz[j] = g1_xz_xxyz[j] - fr * g2_xz_xyz[j];

                g_xxz_xzz[j] = g1_xz_xxzz[j] - fr * g2_xz_xzz[j];

                g_xxz_yyy[j] = g1_xz_xyyy[j] - fr * g2_xz_yyy[j];

                g_xxz_yyz[j] = g1_xz_xyyz[j] - fr * g2_xz_yyz[j];

                g_xxz_yzz[j] = g1_xz_xyzz[j] - fr * g2_xz_yzz[j];

                g_xxz_zzz[j] = g1_xz_xzzz[j] - fr * g2_xz_zzz[j];

                g_xyy_xxx[j] = g1_yy_xxxx[j] - fr * g2_yy_xxx[j];

                g_xyy_xxy[j] = g1_yy_xxxy[j] - fr * g2_yy_xxy[j];

                g_xyy_xxz[j] = g1_yy_xxxz[j] - fr * g2_yy_xxz[j];

                g_xyy_xyy[j] = g1_yy_xxyy[j] - fr * g2_yy_xyy[j];

                g_xyy_xyz[j] = g1_yy_xxyz[j] - fr * g2_yy_xyz[j];

                g_xyy_xzz[j] = g1_yy_xxzz[j] - fr * g2_yy_xzz[j];

                g_xyy_yyy[j] = g1_yy_xyyy[j] - fr * g2_yy_yyy[j];

                g_xyy_yyz[j] = g1_yy_xyyz[j] - fr * g2_yy_yyz[j];

                g_xyy_yzz[j] = g1_yy_xyzz[j] - fr * g2_yy_yzz[j];

                g_xyy_zzz[j] = g1_yy_xzzz[j] - fr * g2_yy_zzz[j];

                g_xyz_xxx[j] = g1_yz_xxxx[j] - fr * g2_yz_xxx[j];

                g_xyz_xxy[j] = g1_yz_xxxy[j] - fr * g2_yz_xxy[j];

                g_xyz_xxz[j] = g1_yz_xxxz[j] - fr * g2_yz_xxz[j];

                g_xyz_xyy[j] = g1_yz_xxyy[j] - fr * g2_yz_xyy[j];

                g_xyz_xyz[j] = g1_yz_xxyz[j] - fr * g2_yz_xyz[j];

                g_xyz_xzz[j] = g1_yz_xxzz[j] - fr * g2_yz_xzz[j];

                g_xyz_yyy[j] = g1_yz_xyyy[j] - fr * g2_yz_yyy[j];

                g_xyz_yyz[j] = g1_yz_xyyz[j] - fr * g2_yz_yyz[j];

                g_xyz_yzz[j] = g1_yz_xyzz[j] - fr * g2_yz_yzz[j];

                g_xyz_zzz[j] = g1_yz_xzzz[j] - fr * g2_yz_zzz[j];

                g_xzz_xxx[j] = g1_zz_xxxx[j] - fr * g2_zz_xxx[j];

                g_xzz_xxy[j] = g1_zz_xxxy[j] - fr * g2_zz_xxy[j];

                g_xzz_xxz[j] = g1_zz_xxxz[j] - fr * g2_zz_xxz[j];

                g_xzz_xyy[j] = g1_zz_xxyy[j] - fr * g2_zz_xyy[j];

                g_xzz_xyz[j] = g1_zz_xxyz[j] - fr * g2_zz_xyz[j];

                g_xzz_xzz[j] = g1_zz_xxzz[j] - fr * g2_zz_xzz[j];

                g_xzz_yyy[j] = g1_zz_xyyy[j] - fr * g2_zz_yyy[j];

                g_xzz_yyz[j] = g1_zz_xyyz[j] - fr * g2_zz_yyz[j];

                g_xzz_yzz[j] = g1_zz_xyzz[j] - fr * g2_zz_yzz[j];

                g_xzz_zzz[j] = g1_zz_xzzz[j] - fr * g2_zz_zzz[j];

                // leading y component

                fr = rcdy[j];

                g_yyy_xxx[j] = g1_yy_xxxy[j] - fr * g2_yy_xxx[j];

                g_yyy_xxy[j] = g1_yy_xxyy[j] - fr * g2_yy_xxy[j];

                g_yyy_xxz[j] = g1_yy_xxyz[j] - fr * g2_yy_xxz[j];

                g_yyy_xyy[j] = g1_yy_xyyy[j] - fr * g2_yy_xyy[j];

                g_yyy_xyz[j] = g1_yy_xyyz[j] - fr * g2_yy_xyz[j];

                g_yyy_xzz[j] = g1_yy_xyzz[j] - fr * g2_yy_xzz[j];

                g_yyy_yyy[j] = g1_yy_yyyy[j] - fr * g2_yy_yyy[j];

                g_yyy_yyz[j] = g1_yy_yyyz[j] - fr * g2_yy_yyz[j];

                g_yyy_yzz[j] = g1_yy_yyzz[j] - fr * g2_yy_yzz[j];

                g_yyy_zzz[j] = g1_yy_yzzz[j] - fr * g2_yy_zzz[j];

                g_yyz_xxx[j] = g1_yz_xxxy[j] - fr * g2_yz_xxx[j];

                g_yyz_xxy[j] = g1_yz_xxyy[j] - fr * g2_yz_xxy[j];

                g_yyz_xxz[j] = g1_yz_xxyz[j] - fr * g2_yz_xxz[j];

                g_yyz_xyy[j] = g1_yz_xyyy[j] - fr * g2_yz_xyy[j];

                g_yyz_xyz[j] = g1_yz_xyyz[j] - fr * g2_yz_xyz[j];

                g_yyz_xzz[j] = g1_yz_xyzz[j] - fr * g2_yz_xzz[j];

                g_yyz_yyy[j] = g1_yz_yyyy[j] - fr * g2_yz_yyy[j];

                g_yyz_yyz[j] = g1_yz_yyyz[j] - fr * g2_yz_yyz[j];

                g_yyz_yzz[j] = g1_yz_yyzz[j] - fr * g2_yz_yzz[j];

                g_yyz_zzz[j] = g1_yz_yzzz[j] - fr * g2_yz_zzz[j];

                g_yzz_xxx[j] = g1_zz_xxxy[j] - fr * g2_zz_xxx[j];

                g_yzz_xxy[j] = g1_zz_xxyy[j] - fr * g2_zz_xxy[j];

                g_yzz_xxz[j] = g1_zz_xxyz[j] - fr * g2_zz_xxz[j];

                g_yzz_xyy[j] = g1_zz_xyyy[j] - fr * g2_zz_xyy[j];

                g_yzz_xyz[j] = g1_zz_xyyz[j] - fr * g2_zz_xyz[j];

                g_yzz_xzz[j] = g1_zz_xyzz[j] - fr * g2_zz_xzz[j];

                g_yzz_yyy[j] = g1_zz_yyyy[j] - fr * g2_zz_yyy[j];

                g_yzz_yyz[j] = g1_zz_yyyz[j] - fr * g2_zz_yyz[j];

                g_yzz_yzz[j] = g1_zz_yyzz[j] - fr * g2_zz_yzz[j];

                g_yzz_zzz[j] = g1_zz_yzzz[j] - fr * g2_zz_zzz[j];

                // leading z component

                fr = rcdz[j];

                g_zzz_xxx[j] = g1_zz_xxxz[j] - fr * g2_zz_xxx[j];

                g_zzz_xxy[j] = g1_zz_xxyz[j] - fr * g2_zz_xxy[j];

                g_zzz_xxz[j] = g1_zz_xxzz[j] - fr * g2_zz_xxz[j];

                g_zzz_xyy[j] = g1_zz_xyyz[j] - fr * g2_zz_xyy[j];

                g_zzz_xyz[j] = g1_zz_xyzz[j] - fr * g2_zz_xyz[j];

                g_zzz_xzz[j] = g1_zz_xzzz[j] - fr * g2_zz_xzz[j];

                g_zzz_yyy[j] = g1_zz_yyyz[j] - fr * g2_zz_yyy[j];

                g_zzz_yyz[j] = g1_zz_yyzz[j] - fr * g2_zz_yyz[j];

                g_zzz_yzz[j] = g1_zz_yzzz[j] - fr * g2_zz_yzz[j];

                g_zzz_zzz[j] = g1_zz_zzzz[j] - fr * g2_zz_zzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXFG(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 3, 4})) return;

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
                                             {braAngularMomentum, 3, 4});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 5});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 4});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|DG)^(m) integrals

            auto g2_xx_xxxx = contrBuffer.data(g2off + 90 * i);

            auto g2_xx_xxxy = contrBuffer.data(g2off + 90 * i + 1);

            auto g2_xx_xxxz = contrBuffer.data(g2off + 90 * i + 2);

            auto g2_xx_xxyy = contrBuffer.data(g2off + 90 * i + 3);

            auto g2_xx_xxyz = contrBuffer.data(g2off + 90 * i + 4);

            auto g2_xx_xxzz = contrBuffer.data(g2off + 90 * i + 5);

            auto g2_xx_xyyy = contrBuffer.data(g2off + 90 * i + 6);

            auto g2_xx_xyyz = contrBuffer.data(g2off + 90 * i + 7);

            auto g2_xx_xyzz = contrBuffer.data(g2off + 90 * i + 8);

            auto g2_xx_xzzz = contrBuffer.data(g2off + 90 * i + 9);

            auto g2_xx_yyyy = contrBuffer.data(g2off + 90 * i + 10);

            auto g2_xx_yyyz = contrBuffer.data(g2off + 90 * i + 11);

            auto g2_xx_yyzz = contrBuffer.data(g2off + 90 * i + 12);

            auto g2_xx_yzzz = contrBuffer.data(g2off + 90 * i + 13);

            auto g2_xx_zzzz = contrBuffer.data(g2off + 90 * i + 14);

            auto g2_xy_xxxx = contrBuffer.data(g2off + 90 * i + 15);

            auto g2_xy_xxxy = contrBuffer.data(g2off + 90 * i + 16);

            auto g2_xy_xxxz = contrBuffer.data(g2off + 90 * i + 17);

            auto g2_xy_xxyy = contrBuffer.data(g2off + 90 * i + 18);

            auto g2_xy_xxyz = contrBuffer.data(g2off + 90 * i + 19);

            auto g2_xy_xxzz = contrBuffer.data(g2off + 90 * i + 20);

            auto g2_xy_xyyy = contrBuffer.data(g2off + 90 * i + 21);

            auto g2_xy_xyyz = contrBuffer.data(g2off + 90 * i + 22);

            auto g2_xy_xyzz = contrBuffer.data(g2off + 90 * i + 23);

            auto g2_xy_xzzz = contrBuffer.data(g2off + 90 * i + 24);

            auto g2_xy_yyyy = contrBuffer.data(g2off + 90 * i + 25);

            auto g2_xy_yyyz = contrBuffer.data(g2off + 90 * i + 26);

            auto g2_xy_yyzz = contrBuffer.data(g2off + 90 * i + 27);

            auto g2_xy_yzzz = contrBuffer.data(g2off + 90 * i + 28);

            auto g2_xy_zzzz = contrBuffer.data(g2off + 90 * i + 29);

            auto g2_xz_xxxx = contrBuffer.data(g2off + 90 * i + 30);

            auto g2_xz_xxxy = contrBuffer.data(g2off + 90 * i + 31);

            auto g2_xz_xxxz = contrBuffer.data(g2off + 90 * i + 32);

            auto g2_xz_xxyy = contrBuffer.data(g2off + 90 * i + 33);

            auto g2_xz_xxyz = contrBuffer.data(g2off + 90 * i + 34);

            auto g2_xz_xxzz = contrBuffer.data(g2off + 90 * i + 35);

            auto g2_xz_xyyy = contrBuffer.data(g2off + 90 * i + 36);

            auto g2_xz_xyyz = contrBuffer.data(g2off + 90 * i + 37);

            auto g2_xz_xyzz = contrBuffer.data(g2off + 90 * i + 38);

            auto g2_xz_xzzz = contrBuffer.data(g2off + 90 * i + 39);

            auto g2_xz_yyyy = contrBuffer.data(g2off + 90 * i + 40);

            auto g2_xz_yyyz = contrBuffer.data(g2off + 90 * i + 41);

            auto g2_xz_yyzz = contrBuffer.data(g2off + 90 * i + 42);

            auto g2_xz_yzzz = contrBuffer.data(g2off + 90 * i + 43);

            auto g2_xz_zzzz = contrBuffer.data(g2off + 90 * i + 44);

            auto g2_yy_xxxx = contrBuffer.data(g2off + 90 * i + 45);

            auto g2_yy_xxxy = contrBuffer.data(g2off + 90 * i + 46);

            auto g2_yy_xxxz = contrBuffer.data(g2off + 90 * i + 47);

            auto g2_yy_xxyy = contrBuffer.data(g2off + 90 * i + 48);

            auto g2_yy_xxyz = contrBuffer.data(g2off + 90 * i + 49);

            auto g2_yy_xxzz = contrBuffer.data(g2off + 90 * i + 50);

            auto g2_yy_xyyy = contrBuffer.data(g2off + 90 * i + 51);

            auto g2_yy_xyyz = contrBuffer.data(g2off + 90 * i + 52);

            auto g2_yy_xyzz = contrBuffer.data(g2off + 90 * i + 53);

            auto g2_yy_xzzz = contrBuffer.data(g2off + 90 * i + 54);

            auto g2_yy_yyyy = contrBuffer.data(g2off + 90 * i + 55);

            auto g2_yy_yyyz = contrBuffer.data(g2off + 90 * i + 56);

            auto g2_yy_yyzz = contrBuffer.data(g2off + 90 * i + 57);

            auto g2_yy_yzzz = contrBuffer.data(g2off + 90 * i + 58);

            auto g2_yy_zzzz = contrBuffer.data(g2off + 90 * i + 59);

            auto g2_yz_xxxx = contrBuffer.data(g2off + 90 * i + 60);

            auto g2_yz_xxxy = contrBuffer.data(g2off + 90 * i + 61);

            auto g2_yz_xxxz = contrBuffer.data(g2off + 90 * i + 62);

            auto g2_yz_xxyy = contrBuffer.data(g2off + 90 * i + 63);

            auto g2_yz_xxyz = contrBuffer.data(g2off + 90 * i + 64);

            auto g2_yz_xxzz = contrBuffer.data(g2off + 90 * i + 65);

            auto g2_yz_xyyy = contrBuffer.data(g2off + 90 * i + 66);

            auto g2_yz_xyyz = contrBuffer.data(g2off + 90 * i + 67);

            auto g2_yz_xyzz = contrBuffer.data(g2off + 90 * i + 68);

            auto g2_yz_xzzz = contrBuffer.data(g2off + 90 * i + 69);

            auto g2_yz_yyyy = contrBuffer.data(g2off + 90 * i + 70);

            auto g2_yz_yyyz = contrBuffer.data(g2off + 90 * i + 71);

            auto g2_yz_yyzz = contrBuffer.data(g2off + 90 * i + 72);

            auto g2_yz_yzzz = contrBuffer.data(g2off + 90 * i + 73);

            auto g2_yz_zzzz = contrBuffer.data(g2off + 90 * i + 74);

            auto g2_zz_xxxx = contrBuffer.data(g2off + 90 * i + 75);

            auto g2_zz_xxxy = contrBuffer.data(g2off + 90 * i + 76);

            auto g2_zz_xxxz = contrBuffer.data(g2off + 90 * i + 77);

            auto g2_zz_xxyy = contrBuffer.data(g2off + 90 * i + 78);

            auto g2_zz_xxyz = contrBuffer.data(g2off + 90 * i + 79);

            auto g2_zz_xxzz = contrBuffer.data(g2off + 90 * i + 80);

            auto g2_zz_xyyy = contrBuffer.data(g2off + 90 * i + 81);

            auto g2_zz_xyyz = contrBuffer.data(g2off + 90 * i + 82);

            auto g2_zz_xyzz = contrBuffer.data(g2off + 90 * i + 83);

            auto g2_zz_xzzz = contrBuffer.data(g2off + 90 * i + 84);

            auto g2_zz_yyyy = contrBuffer.data(g2off + 90 * i + 85);

            auto g2_zz_yyyz = contrBuffer.data(g2off + 90 * i + 86);

            auto g2_zz_yyzz = contrBuffer.data(g2off + 90 * i + 87);

            auto g2_zz_yzzz = contrBuffer.data(g2off + 90 * i + 88);

            auto g2_zz_zzzz = contrBuffer.data(g2off + 90 * i + 89);

            // set up pointers to (X|g(r,r')|DH)^(m) integrals

            auto g1_xx_xxxxx = contrBuffer.data(g1off + 126 * i);

            auto g1_xx_xxxxy = contrBuffer.data(g1off + 126 * i + 1);

            auto g1_xx_xxxxz = contrBuffer.data(g1off + 126 * i + 2);

            auto g1_xx_xxxyy = contrBuffer.data(g1off + 126 * i + 3);

            auto g1_xx_xxxyz = contrBuffer.data(g1off + 126 * i + 4);

            auto g1_xx_xxxzz = contrBuffer.data(g1off + 126 * i + 5);

            auto g1_xx_xxyyy = contrBuffer.data(g1off + 126 * i + 6);

            auto g1_xx_xxyyz = contrBuffer.data(g1off + 126 * i + 7);

            auto g1_xx_xxyzz = contrBuffer.data(g1off + 126 * i + 8);

            auto g1_xx_xxzzz = contrBuffer.data(g1off + 126 * i + 9);

            auto g1_xx_xyyyy = contrBuffer.data(g1off + 126 * i + 10);

            auto g1_xx_xyyyz = contrBuffer.data(g1off + 126 * i + 11);

            auto g1_xx_xyyzz = contrBuffer.data(g1off + 126 * i + 12);

            auto g1_xx_xyzzz = contrBuffer.data(g1off + 126 * i + 13);

            auto g1_xx_xzzzz = contrBuffer.data(g1off + 126 * i + 14);

            auto g1_xy_xxxxx = contrBuffer.data(g1off + 126 * i + 21);

            auto g1_xy_xxxxy = contrBuffer.data(g1off + 126 * i + 22);

            auto g1_xy_xxxxz = contrBuffer.data(g1off + 126 * i + 23);

            auto g1_xy_xxxyy = contrBuffer.data(g1off + 126 * i + 24);

            auto g1_xy_xxxyz = contrBuffer.data(g1off + 126 * i + 25);

            auto g1_xy_xxxzz = contrBuffer.data(g1off + 126 * i + 26);

            auto g1_xy_xxyyy = contrBuffer.data(g1off + 126 * i + 27);

            auto g1_xy_xxyyz = contrBuffer.data(g1off + 126 * i + 28);

            auto g1_xy_xxyzz = contrBuffer.data(g1off + 126 * i + 29);

            auto g1_xy_xxzzz = contrBuffer.data(g1off + 126 * i + 30);

            auto g1_xy_xyyyy = contrBuffer.data(g1off + 126 * i + 31);

            auto g1_xy_xyyyz = contrBuffer.data(g1off + 126 * i + 32);

            auto g1_xy_xyyzz = contrBuffer.data(g1off + 126 * i + 33);

            auto g1_xy_xyzzz = contrBuffer.data(g1off + 126 * i + 34);

            auto g1_xy_xzzzz = contrBuffer.data(g1off + 126 * i + 35);

            auto g1_xz_xxxxx = contrBuffer.data(g1off + 126 * i + 42);

            auto g1_xz_xxxxy = contrBuffer.data(g1off + 126 * i + 43);

            auto g1_xz_xxxxz = contrBuffer.data(g1off + 126 * i + 44);

            auto g1_xz_xxxyy = contrBuffer.data(g1off + 126 * i + 45);

            auto g1_xz_xxxyz = contrBuffer.data(g1off + 126 * i + 46);

            auto g1_xz_xxxzz = contrBuffer.data(g1off + 126 * i + 47);

            auto g1_xz_xxyyy = contrBuffer.data(g1off + 126 * i + 48);

            auto g1_xz_xxyyz = contrBuffer.data(g1off + 126 * i + 49);

            auto g1_xz_xxyzz = contrBuffer.data(g1off + 126 * i + 50);

            auto g1_xz_xxzzz = contrBuffer.data(g1off + 126 * i + 51);

            auto g1_xz_xyyyy = contrBuffer.data(g1off + 126 * i + 52);

            auto g1_xz_xyyyz = contrBuffer.data(g1off + 126 * i + 53);

            auto g1_xz_xyyzz = contrBuffer.data(g1off + 126 * i + 54);

            auto g1_xz_xyzzz = contrBuffer.data(g1off + 126 * i + 55);

            auto g1_xz_xzzzz = contrBuffer.data(g1off + 126 * i + 56);

            auto g1_yy_xxxxx = contrBuffer.data(g1off + 126 * i + 63);

            auto g1_yy_xxxxy = contrBuffer.data(g1off + 126 * i + 64);

            auto g1_yy_xxxxz = contrBuffer.data(g1off + 126 * i + 65);

            auto g1_yy_xxxyy = contrBuffer.data(g1off + 126 * i + 66);

            auto g1_yy_xxxyz = contrBuffer.data(g1off + 126 * i + 67);

            auto g1_yy_xxxzz = contrBuffer.data(g1off + 126 * i + 68);

            auto g1_yy_xxyyy = contrBuffer.data(g1off + 126 * i + 69);

            auto g1_yy_xxyyz = contrBuffer.data(g1off + 126 * i + 70);

            auto g1_yy_xxyzz = contrBuffer.data(g1off + 126 * i + 71);

            auto g1_yy_xxzzz = contrBuffer.data(g1off + 126 * i + 72);

            auto g1_yy_xyyyy = contrBuffer.data(g1off + 126 * i + 73);

            auto g1_yy_xyyyz = contrBuffer.data(g1off + 126 * i + 74);

            auto g1_yy_xyyzz = contrBuffer.data(g1off + 126 * i + 75);

            auto g1_yy_xyzzz = contrBuffer.data(g1off + 126 * i + 76);

            auto g1_yy_xzzzz = contrBuffer.data(g1off + 126 * i + 77);

            auto g1_yy_yyyyy = contrBuffer.data(g1off + 126 * i + 78);

            auto g1_yy_yyyyz = contrBuffer.data(g1off + 126 * i + 79);

            auto g1_yy_yyyzz = contrBuffer.data(g1off + 126 * i + 80);

            auto g1_yy_yyzzz = contrBuffer.data(g1off + 126 * i + 81);

            auto g1_yy_yzzzz = contrBuffer.data(g1off + 126 * i + 82);

            auto g1_yz_xxxxx = contrBuffer.data(g1off + 126 * i + 84);

            auto g1_yz_xxxxy = contrBuffer.data(g1off + 126 * i + 85);

            auto g1_yz_xxxxz = contrBuffer.data(g1off + 126 * i + 86);

            auto g1_yz_xxxyy = contrBuffer.data(g1off + 126 * i + 87);

            auto g1_yz_xxxyz = contrBuffer.data(g1off + 126 * i + 88);

            auto g1_yz_xxxzz = contrBuffer.data(g1off + 126 * i + 89);

            auto g1_yz_xxyyy = contrBuffer.data(g1off + 126 * i + 90);

            auto g1_yz_xxyyz = contrBuffer.data(g1off + 126 * i + 91);

            auto g1_yz_xxyzz = contrBuffer.data(g1off + 126 * i + 92);

            auto g1_yz_xxzzz = contrBuffer.data(g1off + 126 * i + 93);

            auto g1_yz_xyyyy = contrBuffer.data(g1off + 126 * i + 94);

            auto g1_yz_xyyyz = contrBuffer.data(g1off + 126 * i + 95);

            auto g1_yz_xyyzz = contrBuffer.data(g1off + 126 * i + 96);

            auto g1_yz_xyzzz = contrBuffer.data(g1off + 126 * i + 97);

            auto g1_yz_xzzzz = contrBuffer.data(g1off + 126 * i + 98);

            auto g1_yz_yyyyy = contrBuffer.data(g1off + 126 * i + 99);

            auto g1_yz_yyyyz = contrBuffer.data(g1off + 126 * i + 100);

            auto g1_yz_yyyzz = contrBuffer.data(g1off + 126 * i + 101);

            auto g1_yz_yyzzz = contrBuffer.data(g1off + 126 * i + 102);

            auto g1_yz_yzzzz = contrBuffer.data(g1off + 126 * i + 103);

            auto g1_zz_xxxxx = contrBuffer.data(g1off + 126 * i + 105);

            auto g1_zz_xxxxy = contrBuffer.data(g1off + 126 * i + 106);

            auto g1_zz_xxxxz = contrBuffer.data(g1off + 126 * i + 107);

            auto g1_zz_xxxyy = contrBuffer.data(g1off + 126 * i + 108);

            auto g1_zz_xxxyz = contrBuffer.data(g1off + 126 * i + 109);

            auto g1_zz_xxxzz = contrBuffer.data(g1off + 126 * i + 110);

            auto g1_zz_xxyyy = contrBuffer.data(g1off + 126 * i + 111);

            auto g1_zz_xxyyz = contrBuffer.data(g1off + 126 * i + 112);

            auto g1_zz_xxyzz = contrBuffer.data(g1off + 126 * i + 113);

            auto g1_zz_xxzzz = contrBuffer.data(g1off + 126 * i + 114);

            auto g1_zz_xyyyy = contrBuffer.data(g1off + 126 * i + 115);

            auto g1_zz_xyyyz = contrBuffer.data(g1off + 126 * i + 116);

            auto g1_zz_xyyzz = contrBuffer.data(g1off + 126 * i + 117);

            auto g1_zz_xyzzz = contrBuffer.data(g1off + 126 * i + 118);

            auto g1_zz_xzzzz = contrBuffer.data(g1off + 126 * i + 119);

            auto g1_zz_yyyyy = contrBuffer.data(g1off + 126 * i + 120);

            auto g1_zz_yyyyz = contrBuffer.data(g1off + 126 * i + 121);

            auto g1_zz_yyyzz = contrBuffer.data(g1off + 126 * i + 122);

            auto g1_zz_yyzzz = contrBuffer.data(g1off + 126 * i + 123);

            auto g1_zz_yzzzz = contrBuffer.data(g1off + 126 * i + 124);

            auto g1_zz_zzzzz = contrBuffer.data(g1off + 126 * i + 125);

            // set up pointers to (X|g(r,r')|FG)^(m) integrals

            auto g_xxx_xxxx = contrBuffer.data(goff + 150 * i);

            auto g_xxx_xxxy = contrBuffer.data(goff + 150 * i + 1);

            auto g_xxx_xxxz = contrBuffer.data(goff + 150 * i + 2);

            auto g_xxx_xxyy = contrBuffer.data(goff + 150 * i + 3);

            auto g_xxx_xxyz = contrBuffer.data(goff + 150 * i + 4);

            auto g_xxx_xxzz = contrBuffer.data(goff + 150 * i + 5);

            auto g_xxx_xyyy = contrBuffer.data(goff + 150 * i + 6);

            auto g_xxx_xyyz = contrBuffer.data(goff + 150 * i + 7);

            auto g_xxx_xyzz = contrBuffer.data(goff + 150 * i + 8);

            auto g_xxx_xzzz = contrBuffer.data(goff + 150 * i + 9);

            auto g_xxx_yyyy = contrBuffer.data(goff + 150 * i + 10);

            auto g_xxx_yyyz = contrBuffer.data(goff + 150 * i + 11);

            auto g_xxx_yyzz = contrBuffer.data(goff + 150 * i + 12);

            auto g_xxx_yzzz = contrBuffer.data(goff + 150 * i + 13);

            auto g_xxx_zzzz = contrBuffer.data(goff + 150 * i + 14);

            auto g_xxy_xxxx = contrBuffer.data(goff + 150 * i + 15);

            auto g_xxy_xxxy = contrBuffer.data(goff + 150 * i + 16);

            auto g_xxy_xxxz = contrBuffer.data(goff + 150 * i + 17);

            auto g_xxy_xxyy = contrBuffer.data(goff + 150 * i + 18);

            auto g_xxy_xxyz = contrBuffer.data(goff + 150 * i + 19);

            auto g_xxy_xxzz = contrBuffer.data(goff + 150 * i + 20);

            auto g_xxy_xyyy = contrBuffer.data(goff + 150 * i + 21);

            auto g_xxy_xyyz = contrBuffer.data(goff + 150 * i + 22);

            auto g_xxy_xyzz = contrBuffer.data(goff + 150 * i + 23);

            auto g_xxy_xzzz = contrBuffer.data(goff + 150 * i + 24);

            auto g_xxy_yyyy = contrBuffer.data(goff + 150 * i + 25);

            auto g_xxy_yyyz = contrBuffer.data(goff + 150 * i + 26);

            auto g_xxy_yyzz = contrBuffer.data(goff + 150 * i + 27);

            auto g_xxy_yzzz = contrBuffer.data(goff + 150 * i + 28);

            auto g_xxy_zzzz = contrBuffer.data(goff + 150 * i + 29);

            auto g_xxz_xxxx = contrBuffer.data(goff + 150 * i + 30);

            auto g_xxz_xxxy = contrBuffer.data(goff + 150 * i + 31);

            auto g_xxz_xxxz = contrBuffer.data(goff + 150 * i + 32);

            auto g_xxz_xxyy = contrBuffer.data(goff + 150 * i + 33);

            auto g_xxz_xxyz = contrBuffer.data(goff + 150 * i + 34);

            auto g_xxz_xxzz = contrBuffer.data(goff + 150 * i + 35);

            auto g_xxz_xyyy = contrBuffer.data(goff + 150 * i + 36);

            auto g_xxz_xyyz = contrBuffer.data(goff + 150 * i + 37);

            auto g_xxz_xyzz = contrBuffer.data(goff + 150 * i + 38);

            auto g_xxz_xzzz = contrBuffer.data(goff + 150 * i + 39);

            auto g_xxz_yyyy = contrBuffer.data(goff + 150 * i + 40);

            auto g_xxz_yyyz = contrBuffer.data(goff + 150 * i + 41);

            auto g_xxz_yyzz = contrBuffer.data(goff + 150 * i + 42);

            auto g_xxz_yzzz = contrBuffer.data(goff + 150 * i + 43);

            auto g_xxz_zzzz = contrBuffer.data(goff + 150 * i + 44);

            auto g_xyy_xxxx = contrBuffer.data(goff + 150 * i + 45);

            auto g_xyy_xxxy = contrBuffer.data(goff + 150 * i + 46);

            auto g_xyy_xxxz = contrBuffer.data(goff + 150 * i + 47);

            auto g_xyy_xxyy = contrBuffer.data(goff + 150 * i + 48);

            auto g_xyy_xxyz = contrBuffer.data(goff + 150 * i + 49);

            auto g_xyy_xxzz = contrBuffer.data(goff + 150 * i + 50);

            auto g_xyy_xyyy = contrBuffer.data(goff + 150 * i + 51);

            auto g_xyy_xyyz = contrBuffer.data(goff + 150 * i + 52);

            auto g_xyy_xyzz = contrBuffer.data(goff + 150 * i + 53);

            auto g_xyy_xzzz = contrBuffer.data(goff + 150 * i + 54);

            auto g_xyy_yyyy = contrBuffer.data(goff + 150 * i + 55);

            auto g_xyy_yyyz = contrBuffer.data(goff + 150 * i + 56);

            auto g_xyy_yyzz = contrBuffer.data(goff + 150 * i + 57);

            auto g_xyy_yzzz = contrBuffer.data(goff + 150 * i + 58);

            auto g_xyy_zzzz = contrBuffer.data(goff + 150 * i + 59);

            auto g_xyz_xxxx = contrBuffer.data(goff + 150 * i + 60);

            auto g_xyz_xxxy = contrBuffer.data(goff + 150 * i + 61);

            auto g_xyz_xxxz = contrBuffer.data(goff + 150 * i + 62);

            auto g_xyz_xxyy = contrBuffer.data(goff + 150 * i + 63);

            auto g_xyz_xxyz = contrBuffer.data(goff + 150 * i + 64);

            auto g_xyz_xxzz = contrBuffer.data(goff + 150 * i + 65);

            auto g_xyz_xyyy = contrBuffer.data(goff + 150 * i + 66);

            auto g_xyz_xyyz = contrBuffer.data(goff + 150 * i + 67);

            auto g_xyz_xyzz = contrBuffer.data(goff + 150 * i + 68);

            auto g_xyz_xzzz = contrBuffer.data(goff + 150 * i + 69);

            auto g_xyz_yyyy = contrBuffer.data(goff + 150 * i + 70);

            auto g_xyz_yyyz = contrBuffer.data(goff + 150 * i + 71);

            auto g_xyz_yyzz = contrBuffer.data(goff + 150 * i + 72);

            auto g_xyz_yzzz = contrBuffer.data(goff + 150 * i + 73);

            auto g_xyz_zzzz = contrBuffer.data(goff + 150 * i + 74);

            auto g_xzz_xxxx = contrBuffer.data(goff + 150 * i + 75);

            auto g_xzz_xxxy = contrBuffer.data(goff + 150 * i + 76);

            auto g_xzz_xxxz = contrBuffer.data(goff + 150 * i + 77);

            auto g_xzz_xxyy = contrBuffer.data(goff + 150 * i + 78);

            auto g_xzz_xxyz = contrBuffer.data(goff + 150 * i + 79);

            auto g_xzz_xxzz = contrBuffer.data(goff + 150 * i + 80);

            auto g_xzz_xyyy = contrBuffer.data(goff + 150 * i + 81);

            auto g_xzz_xyyz = contrBuffer.data(goff + 150 * i + 82);

            auto g_xzz_xyzz = contrBuffer.data(goff + 150 * i + 83);

            auto g_xzz_xzzz = contrBuffer.data(goff + 150 * i + 84);

            auto g_xzz_yyyy = contrBuffer.data(goff + 150 * i + 85);

            auto g_xzz_yyyz = contrBuffer.data(goff + 150 * i + 86);

            auto g_xzz_yyzz = contrBuffer.data(goff + 150 * i + 87);

            auto g_xzz_yzzz = contrBuffer.data(goff + 150 * i + 88);

            auto g_xzz_zzzz = contrBuffer.data(goff + 150 * i + 89);

            auto g_yyy_xxxx = contrBuffer.data(goff + 150 * i + 90);

            auto g_yyy_xxxy = contrBuffer.data(goff + 150 * i + 91);

            auto g_yyy_xxxz = contrBuffer.data(goff + 150 * i + 92);

            auto g_yyy_xxyy = contrBuffer.data(goff + 150 * i + 93);

            auto g_yyy_xxyz = contrBuffer.data(goff + 150 * i + 94);

            auto g_yyy_xxzz = contrBuffer.data(goff + 150 * i + 95);

            auto g_yyy_xyyy = contrBuffer.data(goff + 150 * i + 96);

            auto g_yyy_xyyz = contrBuffer.data(goff + 150 * i + 97);

            auto g_yyy_xyzz = contrBuffer.data(goff + 150 * i + 98);

            auto g_yyy_xzzz = contrBuffer.data(goff + 150 * i + 99);

            auto g_yyy_yyyy = contrBuffer.data(goff + 150 * i + 100);

            auto g_yyy_yyyz = contrBuffer.data(goff + 150 * i + 101);

            auto g_yyy_yyzz = contrBuffer.data(goff + 150 * i + 102);

            auto g_yyy_yzzz = contrBuffer.data(goff + 150 * i + 103);

            auto g_yyy_zzzz = contrBuffer.data(goff + 150 * i + 104);

            auto g_yyz_xxxx = contrBuffer.data(goff + 150 * i + 105);

            auto g_yyz_xxxy = contrBuffer.data(goff + 150 * i + 106);

            auto g_yyz_xxxz = contrBuffer.data(goff + 150 * i + 107);

            auto g_yyz_xxyy = contrBuffer.data(goff + 150 * i + 108);

            auto g_yyz_xxyz = contrBuffer.data(goff + 150 * i + 109);

            auto g_yyz_xxzz = contrBuffer.data(goff + 150 * i + 110);

            auto g_yyz_xyyy = contrBuffer.data(goff + 150 * i + 111);

            auto g_yyz_xyyz = contrBuffer.data(goff + 150 * i + 112);

            auto g_yyz_xyzz = contrBuffer.data(goff + 150 * i + 113);

            auto g_yyz_xzzz = contrBuffer.data(goff + 150 * i + 114);

            auto g_yyz_yyyy = contrBuffer.data(goff + 150 * i + 115);

            auto g_yyz_yyyz = contrBuffer.data(goff + 150 * i + 116);

            auto g_yyz_yyzz = contrBuffer.data(goff + 150 * i + 117);

            auto g_yyz_yzzz = contrBuffer.data(goff + 150 * i + 118);

            auto g_yyz_zzzz = contrBuffer.data(goff + 150 * i + 119);

            auto g_yzz_xxxx = contrBuffer.data(goff + 150 * i + 120);

            auto g_yzz_xxxy = contrBuffer.data(goff + 150 * i + 121);

            auto g_yzz_xxxz = contrBuffer.data(goff + 150 * i + 122);

            auto g_yzz_xxyy = contrBuffer.data(goff + 150 * i + 123);

            auto g_yzz_xxyz = contrBuffer.data(goff + 150 * i + 124);

            auto g_yzz_xxzz = contrBuffer.data(goff + 150 * i + 125);

            auto g_yzz_xyyy = contrBuffer.data(goff + 150 * i + 126);

            auto g_yzz_xyyz = contrBuffer.data(goff + 150 * i + 127);

            auto g_yzz_xyzz = contrBuffer.data(goff + 150 * i + 128);

            auto g_yzz_xzzz = contrBuffer.data(goff + 150 * i + 129);

            auto g_yzz_yyyy = contrBuffer.data(goff + 150 * i + 130);

            auto g_yzz_yyyz = contrBuffer.data(goff + 150 * i + 131);

            auto g_yzz_yyzz = contrBuffer.data(goff + 150 * i + 132);

            auto g_yzz_yzzz = contrBuffer.data(goff + 150 * i + 133);

            auto g_yzz_zzzz = contrBuffer.data(goff + 150 * i + 134);

            auto g_zzz_xxxx = contrBuffer.data(goff + 150 * i + 135);

            auto g_zzz_xxxy = contrBuffer.data(goff + 150 * i + 136);

            auto g_zzz_xxxz = contrBuffer.data(goff + 150 * i + 137);

            auto g_zzz_xxyy = contrBuffer.data(goff + 150 * i + 138);

            auto g_zzz_xxyz = contrBuffer.data(goff + 150 * i + 139);

            auto g_zzz_xxzz = contrBuffer.data(goff + 150 * i + 140);

            auto g_zzz_xyyy = contrBuffer.data(goff + 150 * i + 141);

            auto g_zzz_xyyz = contrBuffer.data(goff + 150 * i + 142);

            auto g_zzz_xyzz = contrBuffer.data(goff + 150 * i + 143);

            auto g_zzz_xzzz = contrBuffer.data(goff + 150 * i + 144);

            auto g_zzz_yyyy = contrBuffer.data(goff + 150 * i + 145);

            auto g_zzz_yyyz = contrBuffer.data(goff + 150 * i + 146);

            auto g_zzz_yyzz = contrBuffer.data(goff + 150 * i + 147);

            auto g_zzz_yzzz = contrBuffer.data(goff + 150 * i + 148);

            auto g_zzz_zzzz = contrBuffer.data(goff + 150 * i + 149);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxxx, g2_xx_xxxy,\
                                     g2_xx_xxxz, g2_xx_xxyy, g2_xx_xxyz, g2_xx_xxzz,\
                                     g2_xx_xyyy, g2_xx_xyyz, g2_xx_xyzz, g2_xx_xzzz,\
                                     g2_xx_yyyy, g2_xx_yyyz, g2_xx_yyzz, g2_xx_yzzz,\
                                     g2_xx_zzzz, g2_xy_xxxx, g2_xy_xxxy, g2_xy_xxxz,\
                                     g2_xy_xxyy, g2_xy_xxyz, g2_xy_xxzz, g2_xy_xyyy,\
                                     g2_xy_xyyz, g2_xy_xyzz, g2_xy_xzzz, g2_xy_yyyy,\
                                     g2_xy_yyyz, g2_xy_yyzz, g2_xy_yzzz, g2_xy_zzzz,\
                                     g2_xz_xxxx, g2_xz_xxxy, g2_xz_xxxz, g2_xz_xxyy,\
                                     g2_xz_xxyz, g2_xz_xxzz, g2_xz_xyyy, g2_xz_xyyz,\
                                     g2_xz_xyzz, g2_xz_xzzz, g2_xz_yyyy, g2_xz_yyyz,\
                                     g2_xz_yyzz, g2_xz_yzzz, g2_xz_zzzz, g2_yy_xxxx,\
                                     g2_yy_xxxy, g2_yy_xxxz, g2_yy_xxyy, g2_yy_xxyz,\
                                     g2_yy_xxzz, g2_yy_xyyy, g2_yy_xyyz, g2_yy_xyzz,\
                                     g2_yy_xzzz, g2_yy_yyyy, g2_yy_yyyz, g2_yy_yyzz,\
                                     g2_yy_yzzz, g2_yy_zzzz, g2_yz_xxxx, g2_yz_xxxy,\
                                     g2_yz_xxxz, g2_yz_xxyy, g2_yz_xxyz, g2_yz_xxzz,\
                                     g2_yz_xyyy, g2_yz_xyyz, g2_yz_xyzz, g2_yz_xzzz,\
                                     g2_yz_yyyy, g2_yz_yyyz, g2_yz_yyzz, g2_yz_yzzz,\
                                     g2_yz_zzzz, g2_zz_xxxx, g2_zz_xxxy, g2_zz_xxxz,\
                                     g2_zz_xxyy, g2_zz_xxyz, g2_zz_xxzz, g2_zz_xyyy,\
                                     g2_zz_xyyz, g2_zz_xyzz, g2_zz_xzzz, g2_zz_yyyy,\
                                     g2_zz_yyyz, g2_zz_yyzz, g2_zz_yzzz, g2_zz_zzzz,\
                                     g1_xx_xxxxx, g1_xx_xxxxy, g1_xx_xxxxz, g1_xx_xxxyy,\
                                     g1_xx_xxxyz, g1_xx_xxxzz, g1_xx_xxyyy, g1_xx_xxyyz,\
                                     g1_xx_xxyzz, g1_xx_xxzzz, g1_xx_xyyyy, g1_xx_xyyyz,\
                                     g1_xx_xyyzz, g1_xx_xyzzz, g1_xx_xzzzz, g1_xy_xxxxx,\
                                     g1_xy_xxxxy, g1_xy_xxxxz,\
                                     g1_xy_xxxyy, g1_xy_xxxyz, g1_xy_xxxzz, g1_xy_xxyyy,\
                                     g1_xy_xxyyz, g1_xy_xxyzz, g1_xy_xxzzz, g1_xy_xyyyy,\
                                     g1_xy_xyyyz, g1_xy_xyyzz, g1_xy_xyzzz, g1_xy_xzzzz,\
                                     g1_xz_xxxxx, g1_xz_xxxxy,\
                                     g1_xz_xxxxz, g1_xz_xxxyy, g1_xz_xxxyz, g1_xz_xxxzz,\
                                     g1_xz_xxyyy, g1_xz_xxyyz, g1_xz_xxyzz, g1_xz_xxzzz,\
                                     g1_xz_xyyyy, g1_xz_xyyyz, g1_xz_xyyzz, g1_xz_xyzzz,\
                                     g1_xz_xzzzz, g1_yy_xxxxx,\
                                     g1_yy_xxxxy, g1_yy_xxxxz, g1_yy_xxxyy, g1_yy_xxxyz,\
                                     g1_yy_xxxzz, g1_yy_xxyyy, g1_yy_xxyyz, g1_yy_xxyzz,\
                                     g1_yy_xxzzz, g1_yy_xyyyy, g1_yy_xyyyz, g1_yy_xyyzz,\
                                     g1_yy_xyzzz, g1_yy_xzzzz, g1_yy_yyyyy, g1_yy_yyyyz,\
                                     g1_yy_yyyzz, g1_yy_yyzzz, g1_yy_yzzzz, \
                                     g1_yz_xxxxx, g1_yz_xxxxy, g1_yz_xxxxz, g1_yz_xxxyy,\
                                     g1_yz_xxxyz, g1_yz_xxxzz, g1_yz_xxyyy, g1_yz_xxyyz,\
                                     g1_yz_xxyzz, g1_yz_xxzzz, g1_yz_xyyyy, g1_yz_xyyyz,\
                                     g1_yz_xyyzz, g1_yz_xyzzz, g1_yz_xzzzz, g1_yz_yyyyy,\
                                     g1_yz_yyyyz, g1_yz_yyyzz, g1_yz_yyzzz, g1_yz_yzzzz,\
                                     g1_zz_xxxxx, g1_zz_xxxxy, g1_zz_xxxxz,\
                                     g1_zz_xxxyy, g1_zz_xxxyz, g1_zz_xxxzz, g1_zz_xxyyy,\
                                     g1_zz_xxyyz, g1_zz_xxyzz, g1_zz_xxzzz, g1_zz_xyyyy,\
                                     g1_zz_xyyyz, g1_zz_xyyzz, g1_zz_xyzzz, g1_zz_xzzzz,\
                                     g1_zz_yyyyy, g1_zz_yyyyz, g1_zz_yyyzz, g1_zz_yyzzz,\
                                     g1_zz_yzzzz, g1_zz_zzzzz, g_xxx_xxxx, g_xxx_xxxy,\
                                     g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz,\
                                     g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz, g_xxx_xzzz,\
                                     g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz,\
                                     g_xxx_zzzz, g_xxy_xxxx, g_xxy_xxxy, g_xxy_xxxz,\
                                     g_xxy_xxyy, g_xxy_xxyz, g_xxy_xxzz, g_xxy_xyyy,\
                                     g_xxy_xyyz, g_xxy_xyzz, g_xxy_xzzz, g_xxy_yyyy,\
                                     g_xxy_yyyz, g_xxy_yyzz, g_xxy_yzzz, g_xxy_zzzz,\
                                     g_xxz_xxxx, g_xxz_xxxy, g_xxz_xxxz, g_xxz_xxyy,\
                                     g_xxz_xxyz, g_xxz_xxzz, g_xxz_xyyy, g_xxz_xyyz,\
                                     g_xxz_xyzz, g_xxz_xzzz, g_xxz_yyyy, g_xxz_yyyz,\
                                     g_xxz_yyzz, g_xxz_yzzz, g_xxz_zzzz, g_xyy_xxxx,\
                                     g_xyy_xxxy, g_xyy_xxxz, g_xyy_xxyy, g_xyy_xxyz,\
                                     g_xyy_xxzz, g_xyy_xyyy, g_xyy_xyyz, g_xyy_xyzz,\
                                     g_xyy_xzzz, g_xyy_yyyy, g_xyy_yyyz, g_xyy_yyzz,\
                                     g_xyy_yzzz, g_xyy_zzzz, g_xyz_xxxx, g_xyz_xxxy,\
                                     g_xyz_xxxz, g_xyz_xxyy, g_xyz_xxyz, g_xyz_xxzz,\
                                     g_xyz_xyyy, g_xyz_xyyz, g_xyz_xyzz, g_xyz_xzzz,\
                                     g_xyz_yyyy, g_xyz_yyyz, g_xyz_yyzz, g_xyz_yzzz,\
                                     g_xyz_zzzz, g_xzz_xxxx, g_xzz_xxxy, g_xzz_xxxz,\
                                     g_xzz_xxyy, g_xzz_xxyz, g_xzz_xxzz, g_xzz_xyyy,\
                                     g_xzz_xyyz, g_xzz_xyzz, g_xzz_xzzz, g_xzz_yyyy,\
                                     g_xzz_yyyz, g_xzz_yyzz, g_xzz_yzzz, g_xzz_zzzz,\
                                     g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxyy,\
                                     g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz,\
                                     g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz,\
                                     g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz, g_yyz_xxxx,\
                                     g_yyz_xxxy, g_yyz_xxxz, g_yyz_xxyy, g_yyz_xxyz,\
                                     g_yyz_xxzz, g_yyz_xyyy, g_yyz_xyyz, g_yyz_xyzz,\
                                     g_yyz_xzzz, g_yyz_yyyy, g_yyz_yyyz, g_yyz_yyzz,\
                                     g_yyz_yzzz, g_yyz_zzzz, g_yzz_xxxx, g_yzz_xxxy,\
                                     g_yzz_xxxz, g_yzz_xxyy, g_yzz_xxyz, g_yzz_xxzz,\
                                     g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyzz, g_yzz_xzzz,\
                                     g_yzz_yyyy, g_yzz_yyyz, g_yzz_yyzz, g_yzz_yzzz,\
                                     g_yzz_zzzz, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz,\
                                     g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz, g_zzz_xyyy,\
                                     g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy,\
                                     g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xxx_xxxx[j] = g1_xx_xxxxx[j] - fr * g2_xx_xxxx[j];

                g_xxx_xxxy[j] = g1_xx_xxxxy[j] - fr * g2_xx_xxxy[j];

                g_xxx_xxxz[j] = g1_xx_xxxxz[j] - fr * g2_xx_xxxz[j];

                g_xxx_xxyy[j] = g1_xx_xxxyy[j] - fr * g2_xx_xxyy[j];

                g_xxx_xxyz[j] = g1_xx_xxxyz[j] - fr * g2_xx_xxyz[j];

                g_xxx_xxzz[j] = g1_xx_xxxzz[j] - fr * g2_xx_xxzz[j];

                g_xxx_xyyy[j] = g1_xx_xxyyy[j] - fr * g2_xx_xyyy[j];

                g_xxx_xyyz[j] = g1_xx_xxyyz[j] - fr * g2_xx_xyyz[j];

                g_xxx_xyzz[j] = g1_xx_xxyzz[j] - fr * g2_xx_xyzz[j];

                g_xxx_xzzz[j] = g1_xx_xxzzz[j] - fr * g2_xx_xzzz[j];

                g_xxx_yyyy[j] = g1_xx_xyyyy[j] - fr * g2_xx_yyyy[j];

                g_xxx_yyyz[j] = g1_xx_xyyyz[j] - fr * g2_xx_yyyz[j];

                g_xxx_yyzz[j] = g1_xx_xyyzz[j] - fr * g2_xx_yyzz[j];

                g_xxx_yzzz[j] = g1_xx_xyzzz[j] - fr * g2_xx_yzzz[j];

                g_xxx_zzzz[j] = g1_xx_xzzzz[j] - fr * g2_xx_zzzz[j];

                g_xxy_xxxx[j] = g1_xy_xxxxx[j] - fr * g2_xy_xxxx[j];

                g_xxy_xxxy[j] = g1_xy_xxxxy[j] - fr * g2_xy_xxxy[j];

                g_xxy_xxxz[j] = g1_xy_xxxxz[j] - fr * g2_xy_xxxz[j];

                g_xxy_xxyy[j] = g1_xy_xxxyy[j] - fr * g2_xy_xxyy[j];

                g_xxy_xxyz[j] = g1_xy_xxxyz[j] - fr * g2_xy_xxyz[j];

                g_xxy_xxzz[j] = g1_xy_xxxzz[j] - fr * g2_xy_xxzz[j];

                g_xxy_xyyy[j] = g1_xy_xxyyy[j] - fr * g2_xy_xyyy[j];

                g_xxy_xyyz[j] = g1_xy_xxyyz[j] - fr * g2_xy_xyyz[j];

                g_xxy_xyzz[j] = g1_xy_xxyzz[j] - fr * g2_xy_xyzz[j];

                g_xxy_xzzz[j] = g1_xy_xxzzz[j] - fr * g2_xy_xzzz[j];

                g_xxy_yyyy[j] = g1_xy_xyyyy[j] - fr * g2_xy_yyyy[j];

                g_xxy_yyyz[j] = g1_xy_xyyyz[j] - fr * g2_xy_yyyz[j];

                g_xxy_yyzz[j] = g1_xy_xyyzz[j] - fr * g2_xy_yyzz[j];

                g_xxy_yzzz[j] = g1_xy_xyzzz[j] - fr * g2_xy_yzzz[j];

                g_xxy_zzzz[j] = g1_xy_xzzzz[j] - fr * g2_xy_zzzz[j];

                g_xxz_xxxx[j] = g1_xz_xxxxx[j] - fr * g2_xz_xxxx[j];

                g_xxz_xxxy[j] = g1_xz_xxxxy[j] - fr * g2_xz_xxxy[j];

                g_xxz_xxxz[j] = g1_xz_xxxxz[j] - fr * g2_xz_xxxz[j];

                g_xxz_xxyy[j] = g1_xz_xxxyy[j] - fr * g2_xz_xxyy[j];

                g_xxz_xxyz[j] = g1_xz_xxxyz[j] - fr * g2_xz_xxyz[j];

                g_xxz_xxzz[j] = g1_xz_xxxzz[j] - fr * g2_xz_xxzz[j];

                g_xxz_xyyy[j] = g1_xz_xxyyy[j] - fr * g2_xz_xyyy[j];

                g_xxz_xyyz[j] = g1_xz_xxyyz[j] - fr * g2_xz_xyyz[j];

                g_xxz_xyzz[j] = g1_xz_xxyzz[j] - fr * g2_xz_xyzz[j];

                g_xxz_xzzz[j] = g1_xz_xxzzz[j] - fr * g2_xz_xzzz[j];

                g_xxz_yyyy[j] = g1_xz_xyyyy[j] - fr * g2_xz_yyyy[j];

                g_xxz_yyyz[j] = g1_xz_xyyyz[j] - fr * g2_xz_yyyz[j];

                g_xxz_yyzz[j] = g1_xz_xyyzz[j] - fr * g2_xz_yyzz[j];

                g_xxz_yzzz[j] = g1_xz_xyzzz[j] - fr * g2_xz_yzzz[j];

                g_xxz_zzzz[j] = g1_xz_xzzzz[j] - fr * g2_xz_zzzz[j];

                g_xyy_xxxx[j] = g1_yy_xxxxx[j] - fr * g2_yy_xxxx[j];

                g_xyy_xxxy[j] = g1_yy_xxxxy[j] - fr * g2_yy_xxxy[j];

                g_xyy_xxxz[j] = g1_yy_xxxxz[j] - fr * g2_yy_xxxz[j];

                g_xyy_xxyy[j] = g1_yy_xxxyy[j] - fr * g2_yy_xxyy[j];

                g_xyy_xxyz[j] = g1_yy_xxxyz[j] - fr * g2_yy_xxyz[j];

                g_xyy_xxzz[j] = g1_yy_xxxzz[j] - fr * g2_yy_xxzz[j];

                g_xyy_xyyy[j] = g1_yy_xxyyy[j] - fr * g2_yy_xyyy[j];

                g_xyy_xyyz[j] = g1_yy_xxyyz[j] - fr * g2_yy_xyyz[j];

                g_xyy_xyzz[j] = g1_yy_xxyzz[j] - fr * g2_yy_xyzz[j];

                g_xyy_xzzz[j] = g1_yy_xxzzz[j] - fr * g2_yy_xzzz[j];

                g_xyy_yyyy[j] = g1_yy_xyyyy[j] - fr * g2_yy_yyyy[j];

                g_xyy_yyyz[j] = g1_yy_xyyyz[j] - fr * g2_yy_yyyz[j];

                g_xyy_yyzz[j] = g1_yy_xyyzz[j] - fr * g2_yy_yyzz[j];

                g_xyy_yzzz[j] = g1_yy_xyzzz[j] - fr * g2_yy_yzzz[j];

                g_xyy_zzzz[j] = g1_yy_xzzzz[j] - fr * g2_yy_zzzz[j];

                g_xyz_xxxx[j] = g1_yz_xxxxx[j] - fr * g2_yz_xxxx[j];

                g_xyz_xxxy[j] = g1_yz_xxxxy[j] - fr * g2_yz_xxxy[j];

                g_xyz_xxxz[j] = g1_yz_xxxxz[j] - fr * g2_yz_xxxz[j];

                g_xyz_xxyy[j] = g1_yz_xxxyy[j] - fr * g2_yz_xxyy[j];

                g_xyz_xxyz[j] = g1_yz_xxxyz[j] - fr * g2_yz_xxyz[j];

                g_xyz_xxzz[j] = g1_yz_xxxzz[j] - fr * g2_yz_xxzz[j];

                g_xyz_xyyy[j] = g1_yz_xxyyy[j] - fr * g2_yz_xyyy[j];

                g_xyz_xyyz[j] = g1_yz_xxyyz[j] - fr * g2_yz_xyyz[j];

                g_xyz_xyzz[j] = g1_yz_xxyzz[j] - fr * g2_yz_xyzz[j];

                g_xyz_xzzz[j] = g1_yz_xxzzz[j] - fr * g2_yz_xzzz[j];

                g_xyz_yyyy[j] = g1_yz_xyyyy[j] - fr * g2_yz_yyyy[j];

                g_xyz_yyyz[j] = g1_yz_xyyyz[j] - fr * g2_yz_yyyz[j];

                g_xyz_yyzz[j] = g1_yz_xyyzz[j] - fr * g2_yz_yyzz[j];

                g_xyz_yzzz[j] = g1_yz_xyzzz[j] - fr * g2_yz_yzzz[j];

                g_xyz_zzzz[j] = g1_yz_xzzzz[j] - fr * g2_yz_zzzz[j];

                g_xzz_xxxx[j] = g1_zz_xxxxx[j] - fr * g2_zz_xxxx[j];

                g_xzz_xxxy[j] = g1_zz_xxxxy[j] - fr * g2_zz_xxxy[j];

                g_xzz_xxxz[j] = g1_zz_xxxxz[j] - fr * g2_zz_xxxz[j];

                g_xzz_xxyy[j] = g1_zz_xxxyy[j] - fr * g2_zz_xxyy[j];

                g_xzz_xxyz[j] = g1_zz_xxxyz[j] - fr * g2_zz_xxyz[j];

                g_xzz_xxzz[j] = g1_zz_xxxzz[j] - fr * g2_zz_xxzz[j];

                g_xzz_xyyy[j] = g1_zz_xxyyy[j] - fr * g2_zz_xyyy[j];

                g_xzz_xyyz[j] = g1_zz_xxyyz[j] - fr * g2_zz_xyyz[j];

                g_xzz_xyzz[j] = g1_zz_xxyzz[j] - fr * g2_zz_xyzz[j];

                g_xzz_xzzz[j] = g1_zz_xxzzz[j] - fr * g2_zz_xzzz[j];

                g_xzz_yyyy[j] = g1_zz_xyyyy[j] - fr * g2_zz_yyyy[j];

                g_xzz_yyyz[j] = g1_zz_xyyyz[j] - fr * g2_zz_yyyz[j];

                g_xzz_yyzz[j] = g1_zz_xyyzz[j] - fr * g2_zz_yyzz[j];

                g_xzz_yzzz[j] = g1_zz_xyzzz[j] - fr * g2_zz_yzzz[j];

                g_xzz_zzzz[j] = g1_zz_xzzzz[j] - fr * g2_zz_zzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yyy_xxxx[j] = g1_yy_xxxxy[j] - fr * g2_yy_xxxx[j];

                g_yyy_xxxy[j] = g1_yy_xxxyy[j] - fr * g2_yy_xxxy[j];

                g_yyy_xxxz[j] = g1_yy_xxxyz[j] - fr * g2_yy_xxxz[j];

                g_yyy_xxyy[j] = g1_yy_xxyyy[j] - fr * g2_yy_xxyy[j];

                g_yyy_xxyz[j] = g1_yy_xxyyz[j] - fr * g2_yy_xxyz[j];

                g_yyy_xxzz[j] = g1_yy_xxyzz[j] - fr * g2_yy_xxzz[j];

                g_yyy_xyyy[j] = g1_yy_xyyyy[j] - fr * g2_yy_xyyy[j];

                g_yyy_xyyz[j] = g1_yy_xyyyz[j] - fr * g2_yy_xyyz[j];

                g_yyy_xyzz[j] = g1_yy_xyyzz[j] - fr * g2_yy_xyzz[j];

                g_yyy_xzzz[j] = g1_yy_xyzzz[j] - fr * g2_yy_xzzz[j];

                g_yyy_yyyy[j] = g1_yy_yyyyy[j] - fr * g2_yy_yyyy[j];

                g_yyy_yyyz[j] = g1_yy_yyyyz[j] - fr * g2_yy_yyyz[j];

                g_yyy_yyzz[j] = g1_yy_yyyzz[j] - fr * g2_yy_yyzz[j];

                g_yyy_yzzz[j] = g1_yy_yyzzz[j] - fr * g2_yy_yzzz[j];

                g_yyy_zzzz[j] = g1_yy_yzzzz[j] - fr * g2_yy_zzzz[j];

                g_yyz_xxxx[j] = g1_yz_xxxxy[j] - fr * g2_yz_xxxx[j];

                g_yyz_xxxy[j] = g1_yz_xxxyy[j] - fr * g2_yz_xxxy[j];

                g_yyz_xxxz[j] = g1_yz_xxxyz[j] - fr * g2_yz_xxxz[j];

                g_yyz_xxyy[j] = g1_yz_xxyyy[j] - fr * g2_yz_xxyy[j];

                g_yyz_xxyz[j] = g1_yz_xxyyz[j] - fr * g2_yz_xxyz[j];

                g_yyz_xxzz[j] = g1_yz_xxyzz[j] - fr * g2_yz_xxzz[j];

                g_yyz_xyyy[j] = g1_yz_xyyyy[j] - fr * g2_yz_xyyy[j];

                g_yyz_xyyz[j] = g1_yz_xyyyz[j] - fr * g2_yz_xyyz[j];

                g_yyz_xyzz[j] = g1_yz_xyyzz[j] - fr * g2_yz_xyzz[j];

                g_yyz_xzzz[j] = g1_yz_xyzzz[j] - fr * g2_yz_xzzz[j];

                g_yyz_yyyy[j] = g1_yz_yyyyy[j] - fr * g2_yz_yyyy[j];

                g_yyz_yyyz[j] = g1_yz_yyyyz[j] - fr * g2_yz_yyyz[j];

                g_yyz_yyzz[j] = g1_yz_yyyzz[j] - fr * g2_yz_yyzz[j];

                g_yyz_yzzz[j] = g1_yz_yyzzz[j] - fr * g2_yz_yzzz[j];

                g_yyz_zzzz[j] = g1_yz_yzzzz[j] - fr * g2_yz_zzzz[j];

                g_yzz_xxxx[j] = g1_zz_xxxxy[j] - fr * g2_zz_xxxx[j];

                g_yzz_xxxy[j] = g1_zz_xxxyy[j] - fr * g2_zz_xxxy[j];

                g_yzz_xxxz[j] = g1_zz_xxxyz[j] - fr * g2_zz_xxxz[j];

                g_yzz_xxyy[j] = g1_zz_xxyyy[j] - fr * g2_zz_xxyy[j];

                g_yzz_xxyz[j] = g1_zz_xxyyz[j] - fr * g2_zz_xxyz[j];

                g_yzz_xxzz[j] = g1_zz_xxyzz[j] - fr * g2_zz_xxzz[j];

                g_yzz_xyyy[j] = g1_zz_xyyyy[j] - fr * g2_zz_xyyy[j];

                g_yzz_xyyz[j] = g1_zz_xyyyz[j] - fr * g2_zz_xyyz[j];

                g_yzz_xyzz[j] = g1_zz_xyyzz[j] - fr * g2_zz_xyzz[j];

                g_yzz_xzzz[j] = g1_zz_xyzzz[j] - fr * g2_zz_xzzz[j];

                g_yzz_yyyy[j] = g1_zz_yyyyy[j] - fr * g2_zz_yyyy[j];

                g_yzz_yyyz[j] = g1_zz_yyyyz[j] - fr * g2_zz_yyyz[j];

                g_yzz_yyzz[j] = g1_zz_yyyzz[j] - fr * g2_zz_yyzz[j];

                g_yzz_yzzz[j] = g1_zz_yyzzz[j] - fr * g2_zz_yzzz[j];

                g_yzz_zzzz[j] = g1_zz_yzzzz[j] - fr * g2_zz_zzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zzz_xxxx[j] = g1_zz_xxxxz[j] - fr * g2_zz_xxxx[j];

                g_zzz_xxxy[j] = g1_zz_xxxyz[j] - fr * g2_zz_xxxy[j];

                g_zzz_xxxz[j] = g1_zz_xxxzz[j] - fr * g2_zz_xxxz[j];

                g_zzz_xxyy[j] = g1_zz_xxyyz[j] - fr * g2_zz_xxyy[j];

                g_zzz_xxyz[j] = g1_zz_xxyzz[j] - fr * g2_zz_xxyz[j];

                g_zzz_xxzz[j] = g1_zz_xxzzz[j] - fr * g2_zz_xxzz[j];

                g_zzz_xyyy[j] = g1_zz_xyyyz[j] - fr * g2_zz_xyyy[j];

                g_zzz_xyyz[j] = g1_zz_xyyzz[j] - fr * g2_zz_xyyz[j];

                g_zzz_xyzz[j] = g1_zz_xyzzz[j] - fr * g2_zz_xyzz[j];

                g_zzz_xzzz[j] = g1_zz_xzzzz[j] - fr * g2_zz_xzzz[j];

                g_zzz_yyyy[j] = g1_zz_yyyyz[j] - fr * g2_zz_yyyy[j];

                g_zzz_yyyz[j] = g1_zz_yyyzz[j] - fr * g2_zz_yyyz[j];

                g_zzz_yyzz[j] = g1_zz_yyzzz[j] - fr * g2_zz_yyzz[j];

                g_zzz_yzzz[j] = g1_zz_yzzzz[j] - fr * g2_zz_yzzz[j];

                g_zzz_zzzz[j] = g1_zz_zzzzz[j] - fr * g2_zz_zzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXFH(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 3, 5})) return;

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
                                             {braAngularMomentum, 3, 5});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 6});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 2, 5});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|DH)^(m) integrals

            auto g2_xx_xxxxx = contrBuffer.data(g2off + 126 * i);

            auto g2_xx_xxxxy = contrBuffer.data(g2off + 126 * i + 1);

            auto g2_xx_xxxxz = contrBuffer.data(g2off + 126 * i + 2);

            auto g2_xx_xxxyy = contrBuffer.data(g2off + 126 * i + 3);

            auto g2_xx_xxxyz = contrBuffer.data(g2off + 126 * i + 4);

            auto g2_xx_xxxzz = contrBuffer.data(g2off + 126 * i + 5);

            auto g2_xx_xxyyy = contrBuffer.data(g2off + 126 * i + 6);

            auto g2_xx_xxyyz = contrBuffer.data(g2off + 126 * i + 7);

            auto g2_xx_xxyzz = contrBuffer.data(g2off + 126 * i + 8);

            auto g2_xx_xxzzz = contrBuffer.data(g2off + 126 * i + 9);

            auto g2_xx_xyyyy = contrBuffer.data(g2off + 126 * i + 10);

            auto g2_xx_xyyyz = contrBuffer.data(g2off + 126 * i + 11);

            auto g2_xx_xyyzz = contrBuffer.data(g2off + 126 * i + 12);

            auto g2_xx_xyzzz = contrBuffer.data(g2off + 126 * i + 13);

            auto g2_xx_xzzzz = contrBuffer.data(g2off + 126 * i + 14);

            auto g2_xx_yyyyy = contrBuffer.data(g2off + 126 * i + 15);

            auto g2_xx_yyyyz = contrBuffer.data(g2off + 126 * i + 16);

            auto g2_xx_yyyzz = contrBuffer.data(g2off + 126 * i + 17);

            auto g2_xx_yyzzz = contrBuffer.data(g2off + 126 * i + 18);

            auto g2_xx_yzzzz = contrBuffer.data(g2off + 126 * i + 19);

            auto g2_xx_zzzzz = contrBuffer.data(g2off + 126 * i + 20);

            auto g2_xy_xxxxx = contrBuffer.data(g2off + 126 * i + 21);

            auto g2_xy_xxxxy = contrBuffer.data(g2off + 126 * i + 22);

            auto g2_xy_xxxxz = contrBuffer.data(g2off + 126 * i + 23);

            auto g2_xy_xxxyy = contrBuffer.data(g2off + 126 * i + 24);

            auto g2_xy_xxxyz = contrBuffer.data(g2off + 126 * i + 25);

            auto g2_xy_xxxzz = contrBuffer.data(g2off + 126 * i + 26);

            auto g2_xy_xxyyy = contrBuffer.data(g2off + 126 * i + 27);

            auto g2_xy_xxyyz = contrBuffer.data(g2off + 126 * i + 28);

            auto g2_xy_xxyzz = contrBuffer.data(g2off + 126 * i + 29);

            auto g2_xy_xxzzz = contrBuffer.data(g2off + 126 * i + 30);

            auto g2_xy_xyyyy = contrBuffer.data(g2off + 126 * i + 31);

            auto g2_xy_xyyyz = contrBuffer.data(g2off + 126 * i + 32);

            auto g2_xy_xyyzz = contrBuffer.data(g2off + 126 * i + 33);

            auto g2_xy_xyzzz = contrBuffer.data(g2off + 126 * i + 34);

            auto g2_xy_xzzzz = contrBuffer.data(g2off + 126 * i + 35);

            auto g2_xy_yyyyy = contrBuffer.data(g2off + 126 * i + 36);

            auto g2_xy_yyyyz = contrBuffer.data(g2off + 126 * i + 37);

            auto g2_xy_yyyzz = contrBuffer.data(g2off + 126 * i + 38);

            auto g2_xy_yyzzz = contrBuffer.data(g2off + 126 * i + 39);

            auto g2_xy_yzzzz = contrBuffer.data(g2off + 126 * i + 40);

            auto g2_xy_zzzzz = contrBuffer.data(g2off + 126 * i + 41);

            auto g2_xz_xxxxx = contrBuffer.data(g2off + 126 * i + 42);

            auto g2_xz_xxxxy = contrBuffer.data(g2off + 126 * i + 43);

            auto g2_xz_xxxxz = contrBuffer.data(g2off + 126 * i + 44);

            auto g2_xz_xxxyy = contrBuffer.data(g2off + 126 * i + 45);

            auto g2_xz_xxxyz = contrBuffer.data(g2off + 126 * i + 46);

            auto g2_xz_xxxzz = contrBuffer.data(g2off + 126 * i + 47);

            auto g2_xz_xxyyy = contrBuffer.data(g2off + 126 * i + 48);

            auto g2_xz_xxyyz = contrBuffer.data(g2off + 126 * i + 49);

            auto g2_xz_xxyzz = contrBuffer.data(g2off + 126 * i + 50);

            auto g2_xz_xxzzz = contrBuffer.data(g2off + 126 * i + 51);

            auto g2_xz_xyyyy = contrBuffer.data(g2off + 126 * i + 52);

            auto g2_xz_xyyyz = contrBuffer.data(g2off + 126 * i + 53);

            auto g2_xz_xyyzz = contrBuffer.data(g2off + 126 * i + 54);

            auto g2_xz_xyzzz = contrBuffer.data(g2off + 126 * i + 55);

            auto g2_xz_xzzzz = contrBuffer.data(g2off + 126 * i + 56);

            auto g2_xz_yyyyy = contrBuffer.data(g2off + 126 * i + 57);

            auto g2_xz_yyyyz = contrBuffer.data(g2off + 126 * i + 58);

            auto g2_xz_yyyzz = contrBuffer.data(g2off + 126 * i + 59);

            auto g2_xz_yyzzz = contrBuffer.data(g2off + 126 * i + 60);

            auto g2_xz_yzzzz = contrBuffer.data(g2off + 126 * i + 61);

            auto g2_xz_zzzzz = contrBuffer.data(g2off + 126 * i + 62);

            auto g2_yy_xxxxx = contrBuffer.data(g2off + 126 * i + 63);

            auto g2_yy_xxxxy = contrBuffer.data(g2off + 126 * i + 64);

            auto g2_yy_xxxxz = contrBuffer.data(g2off + 126 * i + 65);

            auto g2_yy_xxxyy = contrBuffer.data(g2off + 126 * i + 66);

            auto g2_yy_xxxyz = contrBuffer.data(g2off + 126 * i + 67);

            auto g2_yy_xxxzz = contrBuffer.data(g2off + 126 * i + 68);

            auto g2_yy_xxyyy = contrBuffer.data(g2off + 126 * i + 69);

            auto g2_yy_xxyyz = contrBuffer.data(g2off + 126 * i + 70);

            auto g2_yy_xxyzz = contrBuffer.data(g2off + 126 * i + 71);

            auto g2_yy_xxzzz = contrBuffer.data(g2off + 126 * i + 72);

            auto g2_yy_xyyyy = contrBuffer.data(g2off + 126 * i + 73);

            auto g2_yy_xyyyz = contrBuffer.data(g2off + 126 * i + 74);

            auto g2_yy_xyyzz = contrBuffer.data(g2off + 126 * i + 75);

            auto g2_yy_xyzzz = contrBuffer.data(g2off + 126 * i + 76);

            auto g2_yy_xzzzz = contrBuffer.data(g2off + 126 * i + 77);

            auto g2_yy_yyyyy = contrBuffer.data(g2off + 126 * i + 78);

            auto g2_yy_yyyyz = contrBuffer.data(g2off + 126 * i + 79);

            auto g2_yy_yyyzz = contrBuffer.data(g2off + 126 * i + 80);

            auto g2_yy_yyzzz = contrBuffer.data(g2off + 126 * i + 81);

            auto g2_yy_yzzzz = contrBuffer.data(g2off + 126 * i + 82);

            auto g2_yy_zzzzz = contrBuffer.data(g2off + 126 * i + 83);

            auto g2_yz_xxxxx = contrBuffer.data(g2off + 126 * i + 84);

            auto g2_yz_xxxxy = contrBuffer.data(g2off + 126 * i + 85);

            auto g2_yz_xxxxz = contrBuffer.data(g2off + 126 * i + 86);

            auto g2_yz_xxxyy = contrBuffer.data(g2off + 126 * i + 87);

            auto g2_yz_xxxyz = contrBuffer.data(g2off + 126 * i + 88);

            auto g2_yz_xxxzz = contrBuffer.data(g2off + 126 * i + 89);

            auto g2_yz_xxyyy = contrBuffer.data(g2off + 126 * i + 90);

            auto g2_yz_xxyyz = contrBuffer.data(g2off + 126 * i + 91);

            auto g2_yz_xxyzz = contrBuffer.data(g2off + 126 * i + 92);

            auto g2_yz_xxzzz = contrBuffer.data(g2off + 126 * i + 93);

            auto g2_yz_xyyyy = contrBuffer.data(g2off + 126 * i + 94);

            auto g2_yz_xyyyz = contrBuffer.data(g2off + 126 * i + 95);

            auto g2_yz_xyyzz = contrBuffer.data(g2off + 126 * i + 96);

            auto g2_yz_xyzzz = contrBuffer.data(g2off + 126 * i + 97);

            auto g2_yz_xzzzz = contrBuffer.data(g2off + 126 * i + 98);

            auto g2_yz_yyyyy = contrBuffer.data(g2off + 126 * i + 99);

            auto g2_yz_yyyyz = contrBuffer.data(g2off + 126 * i + 100);

            auto g2_yz_yyyzz = contrBuffer.data(g2off + 126 * i + 101);

            auto g2_yz_yyzzz = contrBuffer.data(g2off + 126 * i + 102);

            auto g2_yz_yzzzz = contrBuffer.data(g2off + 126 * i + 103);

            auto g2_yz_zzzzz = contrBuffer.data(g2off + 126 * i + 104);

            auto g2_zz_xxxxx = contrBuffer.data(g2off + 126 * i + 105);

            auto g2_zz_xxxxy = contrBuffer.data(g2off + 126 * i + 106);

            auto g2_zz_xxxxz = contrBuffer.data(g2off + 126 * i + 107);

            auto g2_zz_xxxyy = contrBuffer.data(g2off + 126 * i + 108);

            auto g2_zz_xxxyz = contrBuffer.data(g2off + 126 * i + 109);

            auto g2_zz_xxxzz = contrBuffer.data(g2off + 126 * i + 110);

            auto g2_zz_xxyyy = contrBuffer.data(g2off + 126 * i + 111);

            auto g2_zz_xxyyz = contrBuffer.data(g2off + 126 * i + 112);

            auto g2_zz_xxyzz = contrBuffer.data(g2off + 126 * i + 113);

            auto g2_zz_xxzzz = contrBuffer.data(g2off + 126 * i + 114);

            auto g2_zz_xyyyy = contrBuffer.data(g2off + 126 * i + 115);

            auto g2_zz_xyyyz = contrBuffer.data(g2off + 126 * i + 116);

            auto g2_zz_xyyzz = contrBuffer.data(g2off + 126 * i + 117);

            auto g2_zz_xyzzz = contrBuffer.data(g2off + 126 * i + 118);

            auto g2_zz_xzzzz = contrBuffer.data(g2off + 126 * i + 119);

            auto g2_zz_yyyyy = contrBuffer.data(g2off + 126 * i + 120);

            auto g2_zz_yyyyz = contrBuffer.data(g2off + 126 * i + 121);

            auto g2_zz_yyyzz = contrBuffer.data(g2off + 126 * i + 122);

            auto g2_zz_yyzzz = contrBuffer.data(g2off + 126 * i + 123);

            auto g2_zz_yzzzz = contrBuffer.data(g2off + 126 * i + 124);

            auto g2_zz_zzzzz = contrBuffer.data(g2off + 126 * i + 125);

            // set up pointers to (X|g(r,r')|DI)^(m) integrals

            auto g1_xx_xxxxxx = contrBuffer.data(g1off + 168 * i);

            auto g1_xx_xxxxxy = contrBuffer.data(g1off + 168 * i + 1);

            auto g1_xx_xxxxxz = contrBuffer.data(g1off + 168 * i + 2);

            auto g1_xx_xxxxyy = contrBuffer.data(g1off + 168 * i + 3);

            auto g1_xx_xxxxyz = contrBuffer.data(g1off + 168 * i + 4);

            auto g1_xx_xxxxzz = contrBuffer.data(g1off + 168 * i + 5);

            auto g1_xx_xxxyyy = contrBuffer.data(g1off + 168 * i + 6);

            auto g1_xx_xxxyyz = contrBuffer.data(g1off + 168 * i + 7);

            auto g1_xx_xxxyzz = contrBuffer.data(g1off + 168 * i + 8);

            auto g1_xx_xxxzzz = contrBuffer.data(g1off + 168 * i + 9);

            auto g1_xx_xxyyyy = contrBuffer.data(g1off + 168 * i + 10);

            auto g1_xx_xxyyyz = contrBuffer.data(g1off + 168 * i + 11);

            auto g1_xx_xxyyzz = contrBuffer.data(g1off + 168 * i + 12);

            auto g1_xx_xxyzzz = contrBuffer.data(g1off + 168 * i + 13);

            auto g1_xx_xxzzzz = contrBuffer.data(g1off + 168 * i + 14);

            auto g1_xx_xyyyyy = contrBuffer.data(g1off + 168 * i + 15);

            auto g1_xx_xyyyyz = contrBuffer.data(g1off + 168 * i + 16);

            auto g1_xx_xyyyzz = contrBuffer.data(g1off + 168 * i + 17);

            auto g1_xx_xyyzzz = contrBuffer.data(g1off + 168 * i + 18);

            auto g1_xx_xyzzzz = contrBuffer.data(g1off + 168 * i + 19);

            auto g1_xx_xzzzzz = contrBuffer.data(g1off + 168 * i + 20);

            auto g1_xy_xxxxxx = contrBuffer.data(g1off + 168 * i + 28);

            auto g1_xy_xxxxxy = contrBuffer.data(g1off + 168 * i + 29);

            auto g1_xy_xxxxxz = contrBuffer.data(g1off + 168 * i + 30);

            auto g1_xy_xxxxyy = contrBuffer.data(g1off + 168 * i + 31);

            auto g1_xy_xxxxyz = contrBuffer.data(g1off + 168 * i + 32);

            auto g1_xy_xxxxzz = contrBuffer.data(g1off + 168 * i + 33);

            auto g1_xy_xxxyyy = contrBuffer.data(g1off + 168 * i + 34);

            auto g1_xy_xxxyyz = contrBuffer.data(g1off + 168 * i + 35);

            auto g1_xy_xxxyzz = contrBuffer.data(g1off + 168 * i + 36);

            auto g1_xy_xxxzzz = contrBuffer.data(g1off + 168 * i + 37);

            auto g1_xy_xxyyyy = contrBuffer.data(g1off + 168 * i + 38);

            auto g1_xy_xxyyyz = contrBuffer.data(g1off + 168 * i + 39);

            auto g1_xy_xxyyzz = contrBuffer.data(g1off + 168 * i + 40);

            auto g1_xy_xxyzzz = contrBuffer.data(g1off + 168 * i + 41);

            auto g1_xy_xxzzzz = contrBuffer.data(g1off + 168 * i + 42);

            auto g1_xy_xyyyyy = contrBuffer.data(g1off + 168 * i + 43);

            auto g1_xy_xyyyyz = contrBuffer.data(g1off + 168 * i + 44);

            auto g1_xy_xyyyzz = contrBuffer.data(g1off + 168 * i + 45);

            auto g1_xy_xyyzzz = contrBuffer.data(g1off + 168 * i + 46);

            auto g1_xy_xyzzzz = contrBuffer.data(g1off + 168 * i + 47);

            auto g1_xy_xzzzzz = contrBuffer.data(g1off + 168 * i + 48);

            auto g1_xz_xxxxxx = contrBuffer.data(g1off + 168 * i + 56);

            auto g1_xz_xxxxxy = contrBuffer.data(g1off + 168 * i + 57);

            auto g1_xz_xxxxxz = contrBuffer.data(g1off + 168 * i + 58);

            auto g1_xz_xxxxyy = contrBuffer.data(g1off + 168 * i + 59);

            auto g1_xz_xxxxyz = contrBuffer.data(g1off + 168 * i + 60);

            auto g1_xz_xxxxzz = contrBuffer.data(g1off + 168 * i + 61);

            auto g1_xz_xxxyyy = contrBuffer.data(g1off + 168 * i + 62);

            auto g1_xz_xxxyyz = contrBuffer.data(g1off + 168 * i + 63);

            auto g1_xz_xxxyzz = contrBuffer.data(g1off + 168 * i + 64);

            auto g1_xz_xxxzzz = contrBuffer.data(g1off + 168 * i + 65);

            auto g1_xz_xxyyyy = contrBuffer.data(g1off + 168 * i + 66);

            auto g1_xz_xxyyyz = contrBuffer.data(g1off + 168 * i + 67);

            auto g1_xz_xxyyzz = contrBuffer.data(g1off + 168 * i + 68);

            auto g1_xz_xxyzzz = contrBuffer.data(g1off + 168 * i + 69);

            auto g1_xz_xxzzzz = contrBuffer.data(g1off + 168 * i + 70);

            auto g1_xz_xyyyyy = contrBuffer.data(g1off + 168 * i + 71);

            auto g1_xz_xyyyyz = contrBuffer.data(g1off + 168 * i + 72);

            auto g1_xz_xyyyzz = contrBuffer.data(g1off + 168 * i + 73);

            auto g1_xz_xyyzzz = contrBuffer.data(g1off + 168 * i + 74);

            auto g1_xz_xyzzzz = contrBuffer.data(g1off + 168 * i + 75);

            auto g1_xz_xzzzzz = contrBuffer.data(g1off + 168 * i + 76);

            auto g1_yy_xxxxxx = contrBuffer.data(g1off + 168 * i + 84);

            auto g1_yy_xxxxxy = contrBuffer.data(g1off + 168 * i + 85);

            auto g1_yy_xxxxxz = contrBuffer.data(g1off + 168 * i + 86);

            auto g1_yy_xxxxyy = contrBuffer.data(g1off + 168 * i + 87);

            auto g1_yy_xxxxyz = contrBuffer.data(g1off + 168 * i + 88);

            auto g1_yy_xxxxzz = contrBuffer.data(g1off + 168 * i + 89);

            auto g1_yy_xxxyyy = contrBuffer.data(g1off + 168 * i + 90);

            auto g1_yy_xxxyyz = contrBuffer.data(g1off + 168 * i + 91);

            auto g1_yy_xxxyzz = contrBuffer.data(g1off + 168 * i + 92);

            auto g1_yy_xxxzzz = contrBuffer.data(g1off + 168 * i + 93);

            auto g1_yy_xxyyyy = contrBuffer.data(g1off + 168 * i + 94);

            auto g1_yy_xxyyyz = contrBuffer.data(g1off + 168 * i + 95);

            auto g1_yy_xxyyzz = contrBuffer.data(g1off + 168 * i + 96);

            auto g1_yy_xxyzzz = contrBuffer.data(g1off + 168 * i + 97);

            auto g1_yy_xxzzzz = contrBuffer.data(g1off + 168 * i + 98);

            auto g1_yy_xyyyyy = contrBuffer.data(g1off + 168 * i + 99);

            auto g1_yy_xyyyyz = contrBuffer.data(g1off + 168 * i + 100);

            auto g1_yy_xyyyzz = contrBuffer.data(g1off + 168 * i + 101);

            auto g1_yy_xyyzzz = contrBuffer.data(g1off + 168 * i + 102);

            auto g1_yy_xyzzzz = contrBuffer.data(g1off + 168 * i + 103);

            auto g1_yy_xzzzzz = contrBuffer.data(g1off + 168 * i + 104);

            auto g1_yy_yyyyyy = contrBuffer.data(g1off + 168 * i + 105);

            auto g1_yy_yyyyyz = contrBuffer.data(g1off + 168 * i + 106);

            auto g1_yy_yyyyzz = contrBuffer.data(g1off + 168 * i + 107);

            auto g1_yy_yyyzzz = contrBuffer.data(g1off + 168 * i + 108);

            auto g1_yy_yyzzzz = contrBuffer.data(g1off + 168 * i + 109);

            auto g1_yy_yzzzzz = contrBuffer.data(g1off + 168 * i + 110);

            auto g1_yz_xxxxxx = contrBuffer.data(g1off + 168 * i + 112);

            auto g1_yz_xxxxxy = contrBuffer.data(g1off + 168 * i + 113);

            auto g1_yz_xxxxxz = contrBuffer.data(g1off + 168 * i + 114);

            auto g1_yz_xxxxyy = contrBuffer.data(g1off + 168 * i + 115);

            auto g1_yz_xxxxyz = contrBuffer.data(g1off + 168 * i + 116);

            auto g1_yz_xxxxzz = contrBuffer.data(g1off + 168 * i + 117);

            auto g1_yz_xxxyyy = contrBuffer.data(g1off + 168 * i + 118);

            auto g1_yz_xxxyyz = contrBuffer.data(g1off + 168 * i + 119);

            auto g1_yz_xxxyzz = contrBuffer.data(g1off + 168 * i + 120);

            auto g1_yz_xxxzzz = contrBuffer.data(g1off + 168 * i + 121);

            auto g1_yz_xxyyyy = contrBuffer.data(g1off + 168 * i + 122);

            auto g1_yz_xxyyyz = contrBuffer.data(g1off + 168 * i + 123);

            auto g1_yz_xxyyzz = contrBuffer.data(g1off + 168 * i + 124);

            auto g1_yz_xxyzzz = contrBuffer.data(g1off + 168 * i + 125);

            auto g1_yz_xxzzzz = contrBuffer.data(g1off + 168 * i + 126);

            auto g1_yz_xyyyyy = contrBuffer.data(g1off + 168 * i + 127);

            auto g1_yz_xyyyyz = contrBuffer.data(g1off + 168 * i + 128);

            auto g1_yz_xyyyzz = contrBuffer.data(g1off + 168 * i + 129);

            auto g1_yz_xyyzzz = contrBuffer.data(g1off + 168 * i + 130);

            auto g1_yz_xyzzzz = contrBuffer.data(g1off + 168 * i + 131);

            auto g1_yz_xzzzzz = contrBuffer.data(g1off + 168 * i + 132);

            auto g1_yz_yyyyyy = contrBuffer.data(g1off + 168 * i + 133);

            auto g1_yz_yyyyyz = contrBuffer.data(g1off + 168 * i + 134);

            auto g1_yz_yyyyzz = contrBuffer.data(g1off + 168 * i + 135);

            auto g1_yz_yyyzzz = contrBuffer.data(g1off + 168 * i + 136);

            auto g1_yz_yyzzzz = contrBuffer.data(g1off + 168 * i + 137);

            auto g1_yz_yzzzzz = contrBuffer.data(g1off + 168 * i + 138);

            auto g1_zz_xxxxxx = contrBuffer.data(g1off + 168 * i + 140);

            auto g1_zz_xxxxxy = contrBuffer.data(g1off + 168 * i + 141);

            auto g1_zz_xxxxxz = contrBuffer.data(g1off + 168 * i + 142);

            auto g1_zz_xxxxyy = contrBuffer.data(g1off + 168 * i + 143);

            auto g1_zz_xxxxyz = contrBuffer.data(g1off + 168 * i + 144);

            auto g1_zz_xxxxzz = contrBuffer.data(g1off + 168 * i + 145);

            auto g1_zz_xxxyyy = contrBuffer.data(g1off + 168 * i + 146);

            auto g1_zz_xxxyyz = contrBuffer.data(g1off + 168 * i + 147);

            auto g1_zz_xxxyzz = contrBuffer.data(g1off + 168 * i + 148);

            auto g1_zz_xxxzzz = contrBuffer.data(g1off + 168 * i + 149);

            auto g1_zz_xxyyyy = contrBuffer.data(g1off + 168 * i + 150);

            auto g1_zz_xxyyyz = contrBuffer.data(g1off + 168 * i + 151);

            auto g1_zz_xxyyzz = contrBuffer.data(g1off + 168 * i + 152);

            auto g1_zz_xxyzzz = contrBuffer.data(g1off + 168 * i + 153);

            auto g1_zz_xxzzzz = contrBuffer.data(g1off + 168 * i + 154);

            auto g1_zz_xyyyyy = contrBuffer.data(g1off + 168 * i + 155);

            auto g1_zz_xyyyyz = contrBuffer.data(g1off + 168 * i + 156);

            auto g1_zz_xyyyzz = contrBuffer.data(g1off + 168 * i + 157);

            auto g1_zz_xyyzzz = contrBuffer.data(g1off + 168 * i + 158);

            auto g1_zz_xyzzzz = contrBuffer.data(g1off + 168 * i + 159);

            auto g1_zz_xzzzzz = contrBuffer.data(g1off + 168 * i + 160);

            auto g1_zz_yyyyyy = contrBuffer.data(g1off + 168 * i + 161);

            auto g1_zz_yyyyyz = contrBuffer.data(g1off + 168 * i + 162);

            auto g1_zz_yyyyzz = contrBuffer.data(g1off + 168 * i + 163);

            auto g1_zz_yyyzzz = contrBuffer.data(g1off + 168 * i + 164);

            auto g1_zz_yyzzzz = contrBuffer.data(g1off + 168 * i + 165);

            auto g1_zz_yzzzzz = contrBuffer.data(g1off + 168 * i + 166);

            auto g1_zz_zzzzzz = contrBuffer.data(g1off + 168 * i + 167);

            // set up pointers to (X|g(r,r')|FH)^(m) integrals

            auto g_xxx_xxxxx = contrBuffer.data(goff + 210 * i);

            auto g_xxx_xxxxy = contrBuffer.data(goff + 210 * i + 1);

            auto g_xxx_xxxxz = contrBuffer.data(goff + 210 * i + 2);

            auto g_xxx_xxxyy = contrBuffer.data(goff + 210 * i + 3);

            auto g_xxx_xxxyz = contrBuffer.data(goff + 210 * i + 4);

            auto g_xxx_xxxzz = contrBuffer.data(goff + 210 * i + 5);

            auto g_xxx_xxyyy = contrBuffer.data(goff + 210 * i + 6);

            auto g_xxx_xxyyz = contrBuffer.data(goff + 210 * i + 7);

            auto g_xxx_xxyzz = contrBuffer.data(goff + 210 * i + 8);

            auto g_xxx_xxzzz = contrBuffer.data(goff + 210 * i + 9);

            auto g_xxx_xyyyy = contrBuffer.data(goff + 210 * i + 10);

            auto g_xxx_xyyyz = contrBuffer.data(goff + 210 * i + 11);

            auto g_xxx_xyyzz = contrBuffer.data(goff + 210 * i + 12);

            auto g_xxx_xyzzz = contrBuffer.data(goff + 210 * i + 13);

            auto g_xxx_xzzzz = contrBuffer.data(goff + 210 * i + 14);

            auto g_xxx_yyyyy = contrBuffer.data(goff + 210 * i + 15);

            auto g_xxx_yyyyz = contrBuffer.data(goff + 210 * i + 16);

            auto g_xxx_yyyzz = contrBuffer.data(goff + 210 * i + 17);

            auto g_xxx_yyzzz = contrBuffer.data(goff + 210 * i + 18);

            auto g_xxx_yzzzz = contrBuffer.data(goff + 210 * i + 19);

            auto g_xxx_zzzzz = contrBuffer.data(goff + 210 * i + 20);

            auto g_xxy_xxxxx = contrBuffer.data(goff + 210 * i + 21);

            auto g_xxy_xxxxy = contrBuffer.data(goff + 210 * i + 22);

            auto g_xxy_xxxxz = contrBuffer.data(goff + 210 * i + 23);

            auto g_xxy_xxxyy = contrBuffer.data(goff + 210 * i + 24);

            auto g_xxy_xxxyz = contrBuffer.data(goff + 210 * i + 25);

            auto g_xxy_xxxzz = contrBuffer.data(goff + 210 * i + 26);

            auto g_xxy_xxyyy = contrBuffer.data(goff + 210 * i + 27);

            auto g_xxy_xxyyz = contrBuffer.data(goff + 210 * i + 28);

            auto g_xxy_xxyzz = contrBuffer.data(goff + 210 * i + 29);

            auto g_xxy_xxzzz = contrBuffer.data(goff + 210 * i + 30);

            auto g_xxy_xyyyy = contrBuffer.data(goff + 210 * i + 31);

            auto g_xxy_xyyyz = contrBuffer.data(goff + 210 * i + 32);

            auto g_xxy_xyyzz = contrBuffer.data(goff + 210 * i + 33);

            auto g_xxy_xyzzz = contrBuffer.data(goff + 210 * i + 34);

            auto g_xxy_xzzzz = contrBuffer.data(goff + 210 * i + 35);

            auto g_xxy_yyyyy = contrBuffer.data(goff + 210 * i + 36);

            auto g_xxy_yyyyz = contrBuffer.data(goff + 210 * i + 37);

            auto g_xxy_yyyzz = contrBuffer.data(goff + 210 * i + 38);

            auto g_xxy_yyzzz = contrBuffer.data(goff + 210 * i + 39);

            auto g_xxy_yzzzz = contrBuffer.data(goff + 210 * i + 40);

            auto g_xxy_zzzzz = contrBuffer.data(goff + 210 * i + 41);

            auto g_xxz_xxxxx = contrBuffer.data(goff + 210 * i + 42);

            auto g_xxz_xxxxy = contrBuffer.data(goff + 210 * i + 43);

            auto g_xxz_xxxxz = contrBuffer.data(goff + 210 * i + 44);

            auto g_xxz_xxxyy = contrBuffer.data(goff + 210 * i + 45);

            auto g_xxz_xxxyz = contrBuffer.data(goff + 210 * i + 46);

            auto g_xxz_xxxzz = contrBuffer.data(goff + 210 * i + 47);

            auto g_xxz_xxyyy = contrBuffer.data(goff + 210 * i + 48);

            auto g_xxz_xxyyz = contrBuffer.data(goff + 210 * i + 49);

            auto g_xxz_xxyzz = contrBuffer.data(goff + 210 * i + 50);

            auto g_xxz_xxzzz = contrBuffer.data(goff + 210 * i + 51);

            auto g_xxz_xyyyy = contrBuffer.data(goff + 210 * i + 52);

            auto g_xxz_xyyyz = contrBuffer.data(goff + 210 * i + 53);

            auto g_xxz_xyyzz = contrBuffer.data(goff + 210 * i + 54);

            auto g_xxz_xyzzz = contrBuffer.data(goff + 210 * i + 55);

            auto g_xxz_xzzzz = contrBuffer.data(goff + 210 * i + 56);

            auto g_xxz_yyyyy = contrBuffer.data(goff + 210 * i + 57);

            auto g_xxz_yyyyz = contrBuffer.data(goff + 210 * i + 58);

            auto g_xxz_yyyzz = contrBuffer.data(goff + 210 * i + 59);

            auto g_xxz_yyzzz = contrBuffer.data(goff + 210 * i + 60);

            auto g_xxz_yzzzz = contrBuffer.data(goff + 210 * i + 61);

            auto g_xxz_zzzzz = contrBuffer.data(goff + 210 * i + 62);

            auto g_xyy_xxxxx = contrBuffer.data(goff + 210 * i + 63);

            auto g_xyy_xxxxy = contrBuffer.data(goff + 210 * i + 64);

            auto g_xyy_xxxxz = contrBuffer.data(goff + 210 * i + 65);

            auto g_xyy_xxxyy = contrBuffer.data(goff + 210 * i + 66);

            auto g_xyy_xxxyz = contrBuffer.data(goff + 210 * i + 67);

            auto g_xyy_xxxzz = contrBuffer.data(goff + 210 * i + 68);

            auto g_xyy_xxyyy = contrBuffer.data(goff + 210 * i + 69);

            auto g_xyy_xxyyz = contrBuffer.data(goff + 210 * i + 70);

            auto g_xyy_xxyzz = contrBuffer.data(goff + 210 * i + 71);

            auto g_xyy_xxzzz = contrBuffer.data(goff + 210 * i + 72);

            auto g_xyy_xyyyy = contrBuffer.data(goff + 210 * i + 73);

            auto g_xyy_xyyyz = contrBuffer.data(goff + 210 * i + 74);

            auto g_xyy_xyyzz = contrBuffer.data(goff + 210 * i + 75);

            auto g_xyy_xyzzz = contrBuffer.data(goff + 210 * i + 76);

            auto g_xyy_xzzzz = contrBuffer.data(goff + 210 * i + 77);

            auto g_xyy_yyyyy = contrBuffer.data(goff + 210 * i + 78);

            auto g_xyy_yyyyz = contrBuffer.data(goff + 210 * i + 79);

            auto g_xyy_yyyzz = contrBuffer.data(goff + 210 * i + 80);

            auto g_xyy_yyzzz = contrBuffer.data(goff + 210 * i + 81);

            auto g_xyy_yzzzz = contrBuffer.data(goff + 210 * i + 82);

            auto g_xyy_zzzzz = contrBuffer.data(goff + 210 * i + 83);

            auto g_xyz_xxxxx = contrBuffer.data(goff + 210 * i + 84);

            auto g_xyz_xxxxy = contrBuffer.data(goff + 210 * i + 85);

            auto g_xyz_xxxxz = contrBuffer.data(goff + 210 * i + 86);

            auto g_xyz_xxxyy = contrBuffer.data(goff + 210 * i + 87);

            auto g_xyz_xxxyz = contrBuffer.data(goff + 210 * i + 88);

            auto g_xyz_xxxzz = contrBuffer.data(goff + 210 * i + 89);

            auto g_xyz_xxyyy = contrBuffer.data(goff + 210 * i + 90);

            auto g_xyz_xxyyz = contrBuffer.data(goff + 210 * i + 91);

            auto g_xyz_xxyzz = contrBuffer.data(goff + 210 * i + 92);

            auto g_xyz_xxzzz = contrBuffer.data(goff + 210 * i + 93);

            auto g_xyz_xyyyy = contrBuffer.data(goff + 210 * i + 94);

            auto g_xyz_xyyyz = contrBuffer.data(goff + 210 * i + 95);

            auto g_xyz_xyyzz = contrBuffer.data(goff + 210 * i + 96);

            auto g_xyz_xyzzz = contrBuffer.data(goff + 210 * i + 97);

            auto g_xyz_xzzzz = contrBuffer.data(goff + 210 * i + 98);

            auto g_xyz_yyyyy = contrBuffer.data(goff + 210 * i + 99);

            auto g_xyz_yyyyz = contrBuffer.data(goff + 210 * i + 100);

            auto g_xyz_yyyzz = contrBuffer.data(goff + 210 * i + 101);

            auto g_xyz_yyzzz = contrBuffer.data(goff + 210 * i + 102);

            auto g_xyz_yzzzz = contrBuffer.data(goff + 210 * i + 103);

            auto g_xyz_zzzzz = contrBuffer.data(goff + 210 * i + 104);

            auto g_xzz_xxxxx = contrBuffer.data(goff + 210 * i + 105);

            auto g_xzz_xxxxy = contrBuffer.data(goff + 210 * i + 106);

            auto g_xzz_xxxxz = contrBuffer.data(goff + 210 * i + 107);

            auto g_xzz_xxxyy = contrBuffer.data(goff + 210 * i + 108);

            auto g_xzz_xxxyz = contrBuffer.data(goff + 210 * i + 109);

            auto g_xzz_xxxzz = contrBuffer.data(goff + 210 * i + 110);

            auto g_xzz_xxyyy = contrBuffer.data(goff + 210 * i + 111);

            auto g_xzz_xxyyz = contrBuffer.data(goff + 210 * i + 112);

            auto g_xzz_xxyzz = contrBuffer.data(goff + 210 * i + 113);

            auto g_xzz_xxzzz = contrBuffer.data(goff + 210 * i + 114);

            auto g_xzz_xyyyy = contrBuffer.data(goff + 210 * i + 115);

            auto g_xzz_xyyyz = contrBuffer.data(goff + 210 * i + 116);

            auto g_xzz_xyyzz = contrBuffer.data(goff + 210 * i + 117);

            auto g_xzz_xyzzz = contrBuffer.data(goff + 210 * i + 118);

            auto g_xzz_xzzzz = contrBuffer.data(goff + 210 * i + 119);

            auto g_xzz_yyyyy = contrBuffer.data(goff + 210 * i + 120);

            auto g_xzz_yyyyz = contrBuffer.data(goff + 210 * i + 121);

            auto g_xzz_yyyzz = contrBuffer.data(goff + 210 * i + 122);

            auto g_xzz_yyzzz = contrBuffer.data(goff + 210 * i + 123);

            auto g_xzz_yzzzz = contrBuffer.data(goff + 210 * i + 124);

            auto g_xzz_zzzzz = contrBuffer.data(goff + 210 * i + 125);

            auto g_yyy_xxxxx = contrBuffer.data(goff + 210 * i + 126);

            auto g_yyy_xxxxy = contrBuffer.data(goff + 210 * i + 127);

            auto g_yyy_xxxxz = contrBuffer.data(goff + 210 * i + 128);

            auto g_yyy_xxxyy = contrBuffer.data(goff + 210 * i + 129);

            auto g_yyy_xxxyz = contrBuffer.data(goff + 210 * i + 130);

            auto g_yyy_xxxzz = contrBuffer.data(goff + 210 * i + 131);

            auto g_yyy_xxyyy = contrBuffer.data(goff + 210 * i + 132);

            auto g_yyy_xxyyz = contrBuffer.data(goff + 210 * i + 133);

            auto g_yyy_xxyzz = contrBuffer.data(goff + 210 * i + 134);

            auto g_yyy_xxzzz = contrBuffer.data(goff + 210 * i + 135);

            auto g_yyy_xyyyy = contrBuffer.data(goff + 210 * i + 136);

            auto g_yyy_xyyyz = contrBuffer.data(goff + 210 * i + 137);

            auto g_yyy_xyyzz = contrBuffer.data(goff + 210 * i + 138);

            auto g_yyy_xyzzz = contrBuffer.data(goff + 210 * i + 139);

            auto g_yyy_xzzzz = contrBuffer.data(goff + 210 * i + 140);

            auto g_yyy_yyyyy = contrBuffer.data(goff + 210 * i + 141);

            auto g_yyy_yyyyz = contrBuffer.data(goff + 210 * i + 142);

            auto g_yyy_yyyzz = contrBuffer.data(goff + 210 * i + 143);

            auto g_yyy_yyzzz = contrBuffer.data(goff + 210 * i + 144);

            auto g_yyy_yzzzz = contrBuffer.data(goff + 210 * i + 145);

            auto g_yyy_zzzzz = contrBuffer.data(goff + 210 * i + 146);

            auto g_yyz_xxxxx = contrBuffer.data(goff + 210 * i + 147);

            auto g_yyz_xxxxy = contrBuffer.data(goff + 210 * i + 148);

            auto g_yyz_xxxxz = contrBuffer.data(goff + 210 * i + 149);

            auto g_yyz_xxxyy = contrBuffer.data(goff + 210 * i + 150);

            auto g_yyz_xxxyz = contrBuffer.data(goff + 210 * i + 151);

            auto g_yyz_xxxzz = contrBuffer.data(goff + 210 * i + 152);

            auto g_yyz_xxyyy = contrBuffer.data(goff + 210 * i + 153);

            auto g_yyz_xxyyz = contrBuffer.data(goff + 210 * i + 154);

            auto g_yyz_xxyzz = contrBuffer.data(goff + 210 * i + 155);

            auto g_yyz_xxzzz = contrBuffer.data(goff + 210 * i + 156);

            auto g_yyz_xyyyy = contrBuffer.data(goff + 210 * i + 157);

            auto g_yyz_xyyyz = contrBuffer.data(goff + 210 * i + 158);

            auto g_yyz_xyyzz = contrBuffer.data(goff + 210 * i + 159);

            auto g_yyz_xyzzz = contrBuffer.data(goff + 210 * i + 160);

            auto g_yyz_xzzzz = contrBuffer.data(goff + 210 * i + 161);

            auto g_yyz_yyyyy = contrBuffer.data(goff + 210 * i + 162);

            auto g_yyz_yyyyz = contrBuffer.data(goff + 210 * i + 163);

            auto g_yyz_yyyzz = contrBuffer.data(goff + 210 * i + 164);

            auto g_yyz_yyzzz = contrBuffer.data(goff + 210 * i + 165);

            auto g_yyz_yzzzz = contrBuffer.data(goff + 210 * i + 166);

            auto g_yyz_zzzzz = contrBuffer.data(goff + 210 * i + 167);

            auto g_yzz_xxxxx = contrBuffer.data(goff + 210 * i + 168);

            auto g_yzz_xxxxy = contrBuffer.data(goff + 210 * i + 169);

            auto g_yzz_xxxxz = contrBuffer.data(goff + 210 * i + 170);

            auto g_yzz_xxxyy = contrBuffer.data(goff + 210 * i + 171);

            auto g_yzz_xxxyz = contrBuffer.data(goff + 210 * i + 172);

            auto g_yzz_xxxzz = contrBuffer.data(goff + 210 * i + 173);

            auto g_yzz_xxyyy = contrBuffer.data(goff + 210 * i + 174);

            auto g_yzz_xxyyz = contrBuffer.data(goff + 210 * i + 175);

            auto g_yzz_xxyzz = contrBuffer.data(goff + 210 * i + 176);

            auto g_yzz_xxzzz = contrBuffer.data(goff + 210 * i + 177);

            auto g_yzz_xyyyy = contrBuffer.data(goff + 210 * i + 178);

            auto g_yzz_xyyyz = contrBuffer.data(goff + 210 * i + 179);

            auto g_yzz_xyyzz = contrBuffer.data(goff + 210 * i + 180);

            auto g_yzz_xyzzz = contrBuffer.data(goff + 210 * i + 181);

            auto g_yzz_xzzzz = contrBuffer.data(goff + 210 * i + 182);

            auto g_yzz_yyyyy = contrBuffer.data(goff + 210 * i + 183);

            auto g_yzz_yyyyz = contrBuffer.data(goff + 210 * i + 184);

            auto g_yzz_yyyzz = contrBuffer.data(goff + 210 * i + 185);

            auto g_yzz_yyzzz = contrBuffer.data(goff + 210 * i + 186);

            auto g_yzz_yzzzz = contrBuffer.data(goff + 210 * i + 187);

            auto g_yzz_zzzzz = contrBuffer.data(goff + 210 * i + 188);

            auto g_zzz_xxxxx = contrBuffer.data(goff + 210 * i + 189);

            auto g_zzz_xxxxy = contrBuffer.data(goff + 210 * i + 190);

            auto g_zzz_xxxxz = contrBuffer.data(goff + 210 * i + 191);

            auto g_zzz_xxxyy = contrBuffer.data(goff + 210 * i + 192);

            auto g_zzz_xxxyz = contrBuffer.data(goff + 210 * i + 193);

            auto g_zzz_xxxzz = contrBuffer.data(goff + 210 * i + 194);

            auto g_zzz_xxyyy = contrBuffer.data(goff + 210 * i + 195);

            auto g_zzz_xxyyz = contrBuffer.data(goff + 210 * i + 196);

            auto g_zzz_xxyzz = contrBuffer.data(goff + 210 * i + 197);

            auto g_zzz_xxzzz = contrBuffer.data(goff + 210 * i + 198);

            auto g_zzz_xyyyy = contrBuffer.data(goff + 210 * i + 199);

            auto g_zzz_xyyyz = contrBuffer.data(goff + 210 * i + 200);

            auto g_zzz_xyyzz = contrBuffer.data(goff + 210 * i + 201);

            auto g_zzz_xyzzz = contrBuffer.data(goff + 210 * i + 202);

            auto g_zzz_xzzzz = contrBuffer.data(goff + 210 * i + 203);

            auto g_zzz_yyyyy = contrBuffer.data(goff + 210 * i + 204);

            auto g_zzz_yyyyz = contrBuffer.data(goff + 210 * i + 205);

            auto g_zzz_yyyzz = contrBuffer.data(goff + 210 * i + 206);

            auto g_zzz_yyzzz = contrBuffer.data(goff + 210 * i + 207);

            auto g_zzz_yzzzz = contrBuffer.data(goff + 210 * i + 208);

            auto g_zzz_zzzzz = contrBuffer.data(goff + 210 * i + 209);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxxxx, g2_xx_xxxxy,\
                                     g2_xx_xxxxz, g2_xx_xxxyy, g2_xx_xxxyz, g2_xx_xxxzz,\
                                     g2_xx_xxyyy, g2_xx_xxyyz, g2_xx_xxyzz, g2_xx_xxzzz,\
                                     g2_xx_xyyyy, g2_xx_xyyyz, g2_xx_xyyzz, g2_xx_xyzzz,\
                                     g2_xx_xzzzz, g2_xx_yyyyy, g2_xx_yyyyz, g2_xx_yyyzz,\
                                     g2_xx_yyzzz, g2_xx_yzzzz, g2_xx_zzzzz, g2_xy_xxxxx,\
                                     g2_xy_xxxxy, g2_xy_xxxxz, g2_xy_xxxyy, g2_xy_xxxyz,\
                                     g2_xy_xxxzz, g2_xy_xxyyy, g2_xy_xxyyz, g2_xy_xxyzz,\
                                     g2_xy_xxzzz, g2_xy_xyyyy, g2_xy_xyyyz, g2_xy_xyyzz,\
                                     g2_xy_xyzzz, g2_xy_xzzzz, g2_xy_yyyyy, g2_xy_yyyyz,\
                                     g2_xy_yyyzz, g2_xy_yyzzz, g2_xy_yzzzz, g2_xy_zzzzz,\
                                     g2_xz_xxxxx, g2_xz_xxxxy, g2_xz_xxxxz, g2_xz_xxxyy,\
                                     g2_xz_xxxyz, g2_xz_xxxzz, g2_xz_xxyyy, g2_xz_xxyyz,\
                                     g2_xz_xxyzz, g2_xz_xxzzz, g2_xz_xyyyy, g2_xz_xyyyz,\
                                     g2_xz_xyyzz, g2_xz_xyzzz, g2_xz_xzzzz, g2_xz_yyyyy,\
                                     g2_xz_yyyyz, g2_xz_yyyzz, g2_xz_yyzzz, g2_xz_yzzzz,\
                                     g2_xz_zzzzz, g2_yy_xxxxx, g2_yy_xxxxy, g2_yy_xxxxz,\
                                     g2_yy_xxxyy, g2_yy_xxxyz, g2_yy_xxxzz, g2_yy_xxyyy,\
                                     g2_yy_xxyyz, g2_yy_xxyzz, g2_yy_xxzzz, g2_yy_xyyyy,\
                                     g2_yy_xyyyz, g2_yy_xyyzz, g2_yy_xyzzz, g2_yy_xzzzz,\
                                     g2_yy_yyyyy, g2_yy_yyyyz, g2_yy_yyyzz, g2_yy_yyzzz,\
                                     g2_yy_yzzzz, g2_yy_zzzzz, g2_yz_xxxxx, g2_yz_xxxxy,\
                                     g2_yz_xxxxz, g2_yz_xxxyy, g2_yz_xxxyz, g2_yz_xxxzz,\
                                     g2_yz_xxyyy, g2_yz_xxyyz, g2_yz_xxyzz, g2_yz_xxzzz,\
                                     g2_yz_xyyyy, g2_yz_xyyyz, g2_yz_xyyzz, g2_yz_xyzzz,\
                                     g2_yz_xzzzz, g2_yz_yyyyy, g2_yz_yyyyz, g2_yz_yyyzz,\
                                     g2_yz_yyzzz, g2_yz_yzzzz, g2_yz_zzzzz, g2_zz_xxxxx,\
                                     g2_zz_xxxxy, g2_zz_xxxxz, g2_zz_xxxyy, g2_zz_xxxyz,\
                                     g2_zz_xxxzz, g2_zz_xxyyy, g2_zz_xxyyz, g2_zz_xxyzz,\
                                     g2_zz_xxzzz, g2_zz_xyyyy, g2_zz_xyyyz, g2_zz_xyyzz,\
                                     g2_zz_xyzzz, g2_zz_xzzzz, g2_zz_yyyyy, g2_zz_yyyyz,\
                                     g2_zz_yyyzz, g2_zz_yyzzz, g2_zz_yzzzz, g2_zz_zzzzz,\
                                     g1_xx_xxxxxx, g1_xx_xxxxxy, g1_xx_xxxxxz,\
                                     g1_xx_xxxxyy, g1_xx_xxxxyz, g1_xx_xxxxzz,\
                                     g1_xx_xxxyyy, g1_xx_xxxyyz, g1_xx_xxxyzz,\
                                     g1_xx_xxxzzz, g1_xx_xxyyyy, g1_xx_xxyyyz,\
                                     g1_xx_xxyyzz, g1_xx_xxyzzz, g1_xx_xxzzzz,\
                                     g1_xx_xyyyyy, g1_xx_xyyyyz, g1_xx_xyyyzz,\
                                     g1_xx_xyyzzz, g1_xx_xyzzzz, g1_xx_xzzzzz,\
                                     g1_xy_xxxxxx, g1_xy_xxxxxy,\
                                     g1_xy_xxxxxz, g1_xy_xxxxyy, g1_xy_xxxxyz,\
                                     g1_xy_xxxxzz, g1_xy_xxxyyy, g1_xy_xxxyyz,\
                                     g1_xy_xxxyzz, g1_xy_xxxzzz, g1_xy_xxyyyy,\
                                     g1_xy_xxyyyz, g1_xy_xxyyzz, g1_xy_xxyzzz,\
                                     g1_xy_xxzzzz, g1_xy_xyyyyy, g1_xy_xyyyyz,\
                                     g1_xy_xyyyzz, g1_xy_xyyzzz, g1_xy_xyzzzz,\
                                     g1_xy_xzzzzz, g1_xz_xxxxxx,\
                                     g1_xz_xxxxxy, g1_xz_xxxxxz, g1_xz_xxxxyy,\
                                     g1_xz_xxxxyz, g1_xz_xxxxzz, g1_xz_xxxyyy,\
                                     g1_xz_xxxyyz, g1_xz_xxxyzz, g1_xz_xxxzzz,\
                                     g1_xz_xxyyyy, g1_xz_xxyyyz, g1_xz_xxyyzz,\
                                     g1_xz_xxyzzz, g1_xz_xxzzzz, g1_xz_xyyyyy,\
                                     g1_xz_xyyyyz, g1_xz_xyyyzz, g1_xz_xyyzzz,\
                                     g1_xz_xyzzzz, g1_xz_xzzzzz,\
                                     g1_yy_xxxxxx, g1_yy_xxxxxy, g1_yy_xxxxxz,\
                                     g1_yy_xxxxyy, g1_yy_xxxxyz, g1_yy_xxxxzz,\
                                     g1_yy_xxxyyy, g1_yy_xxxyyz, g1_yy_xxxyzz,\
                                     g1_yy_xxxzzz, g1_yy_xxyyyy, g1_yy_xxyyyz,\
                                     g1_yy_xxyyzz, g1_yy_xxyzzz, g1_yy_xxzzzz,\
                                     g1_yy_xyyyyy, g1_yy_xyyyyz, g1_yy_xyyyzz,\
                                     g1_yy_xyyzzz, g1_yy_xyzzzz, g1_yy_xzzzzz,\
                                     g1_yy_yyyyyy, g1_yy_yyyyyz, g1_yy_yyyyzz,\
                                     g1_yy_yyyzzz, g1_yy_yyzzzz, g1_yy_yzzzzz,\
                                     g1_yz_xxxxxx, g1_yz_xxxxxy,\
                                     g1_yz_xxxxxz, g1_yz_xxxxyy, g1_yz_xxxxyz,\
                                     g1_yz_xxxxzz, g1_yz_xxxyyy, g1_yz_xxxyyz,\
                                     g1_yz_xxxyzz, g1_yz_xxxzzz, g1_yz_xxyyyy,\
                                     g1_yz_xxyyyz, g1_yz_xxyyzz, g1_yz_xxyzzz,\
                                     g1_yz_xxzzzz, g1_yz_xyyyyy, g1_yz_xyyyyz,\
                                     g1_yz_xyyyzz, g1_yz_xyyzzz, g1_yz_xyzzzz,\
                                     g1_yz_xzzzzz, g1_yz_yyyyyy, g1_yz_yyyyyz,\
                                     g1_yz_yyyyzz, g1_yz_yyyzzz, g1_yz_yyzzzz,\
                                     g1_yz_yzzzzz, g1_zz_xxxxxx,\
                                     g1_zz_xxxxxy, g1_zz_xxxxxz, g1_zz_xxxxyy,\
                                     g1_zz_xxxxyz, g1_zz_xxxxzz, g1_zz_xxxyyy,\
                                     g1_zz_xxxyyz, g1_zz_xxxyzz, g1_zz_xxxzzz,\
                                     g1_zz_xxyyyy, g1_zz_xxyyyz, g1_zz_xxyyzz,\
                                     g1_zz_xxyzzz, g1_zz_xxzzzz, g1_zz_xyyyyy,\
                                     g1_zz_xyyyyz, g1_zz_xyyyzz, g1_zz_xyyzzz,\
                                     g1_zz_xyzzzz, g1_zz_xzzzzz, g1_zz_yyyyyy,\
                                     g1_zz_yyyyyz, g1_zz_yyyyzz, g1_zz_yyyzzz,\
                                     g1_zz_yyzzzz, g1_zz_yzzzzz, g1_zz_zzzzzz,\
                                     g_xxx_xxxxx, g_xxx_xxxxy, g_xxx_xxxxz, g_xxx_xxxyy,\
                                     g_xxx_xxxyz, g_xxx_xxxzz, g_xxx_xxyyy, g_xxx_xxyyz,\
                                     g_xxx_xxyzz, g_xxx_xxzzz, g_xxx_xyyyy, g_xxx_xyyyz,\
                                     g_xxx_xyyzz, g_xxx_xyzzz, g_xxx_xzzzz, g_xxx_yyyyy,\
                                     g_xxx_yyyyz, g_xxx_yyyzz, g_xxx_yyzzz, g_xxx_yzzzz,\
                                     g_xxx_zzzzz, g_xxy_xxxxx, g_xxy_xxxxy, g_xxy_xxxxz,\
                                     g_xxy_xxxyy, g_xxy_xxxyz, g_xxy_xxxzz, g_xxy_xxyyy,\
                                     g_xxy_xxyyz, g_xxy_xxyzz, g_xxy_xxzzz, g_xxy_xyyyy,\
                                     g_xxy_xyyyz, g_xxy_xyyzz, g_xxy_xyzzz, g_xxy_xzzzz,\
                                     g_xxy_yyyyy, g_xxy_yyyyz, g_xxy_yyyzz, g_xxy_yyzzz,\
                                     g_xxy_yzzzz, g_xxy_zzzzz, g_xxz_xxxxx, g_xxz_xxxxy,\
                                     g_xxz_xxxxz, g_xxz_xxxyy, g_xxz_xxxyz, g_xxz_xxxzz,\
                                     g_xxz_xxyyy, g_xxz_xxyyz, g_xxz_xxyzz, g_xxz_xxzzz,\
                                     g_xxz_xyyyy, g_xxz_xyyyz, g_xxz_xyyzz, g_xxz_xyzzz,\
                                     g_xxz_xzzzz, g_xxz_yyyyy, g_xxz_yyyyz, g_xxz_yyyzz,\
                                     g_xxz_yyzzz, g_xxz_yzzzz, g_xxz_zzzzz, g_xyy_xxxxx,\
                                     g_xyy_xxxxy, g_xyy_xxxxz, g_xyy_xxxyy, g_xyy_xxxyz,\
                                     g_xyy_xxxzz, g_xyy_xxyyy, g_xyy_xxyyz, g_xyy_xxyzz,\
                                     g_xyy_xxzzz, g_xyy_xyyyy, g_xyy_xyyyz, g_xyy_xyyzz,\
                                     g_xyy_xyzzz, g_xyy_xzzzz, g_xyy_yyyyy, g_xyy_yyyyz,\
                                     g_xyy_yyyzz, g_xyy_yyzzz, g_xyy_yzzzz, g_xyy_zzzzz,\
                                     g_xyz_xxxxx, g_xyz_xxxxy, g_xyz_xxxxz, g_xyz_xxxyy,\
                                     g_xyz_xxxyz, g_xyz_xxxzz, g_xyz_xxyyy, g_xyz_xxyyz,\
                                     g_xyz_xxyzz, g_xyz_xxzzz, g_xyz_xyyyy, g_xyz_xyyyz,\
                                     g_xyz_xyyzz, g_xyz_xyzzz, g_xyz_xzzzz, g_xyz_yyyyy,\
                                     g_xyz_yyyyz, g_xyz_yyyzz, g_xyz_yyzzz, g_xyz_yzzzz,\
                                     g_xyz_zzzzz, g_xzz_xxxxx, g_xzz_xxxxy, g_xzz_xxxxz,\
                                     g_xzz_xxxyy, g_xzz_xxxyz, g_xzz_xxxzz, g_xzz_xxyyy,\
                                     g_xzz_xxyyz, g_xzz_xxyzz, g_xzz_xxzzz, g_xzz_xyyyy,\
                                     g_xzz_xyyyz, g_xzz_xyyzz, g_xzz_xyzzz, g_xzz_xzzzz,\
                                     g_xzz_yyyyy, g_xzz_yyyyz, g_xzz_yyyzz, g_xzz_yyzzz,\
                                     g_xzz_yzzzz, g_xzz_zzzzz, g_yyy_xxxxx, g_yyy_xxxxy,\
                                     g_yyy_xxxxz, g_yyy_xxxyy, g_yyy_xxxyz, g_yyy_xxxzz,\
                                     g_yyy_xxyyy, g_yyy_xxyyz, g_yyy_xxyzz, g_yyy_xxzzz,\
                                     g_yyy_xyyyy, g_yyy_xyyyz, g_yyy_xyyzz, g_yyy_xyzzz,\
                                     g_yyy_xzzzz, g_yyy_yyyyy, g_yyy_yyyyz, g_yyy_yyyzz,\
                                     g_yyy_yyzzz, g_yyy_yzzzz, g_yyy_zzzzz, g_yyz_xxxxx,\
                                     g_yyz_xxxxy, g_yyz_xxxxz, g_yyz_xxxyy, g_yyz_xxxyz,\
                                     g_yyz_xxxzz, g_yyz_xxyyy, g_yyz_xxyyz, g_yyz_xxyzz,\
                                     g_yyz_xxzzz, g_yyz_xyyyy, g_yyz_xyyyz, g_yyz_xyyzz,\
                                     g_yyz_xyzzz, g_yyz_xzzzz, g_yyz_yyyyy, g_yyz_yyyyz,\
                                     g_yyz_yyyzz, g_yyz_yyzzz, g_yyz_yzzzz, g_yyz_zzzzz,\
                                     g_yzz_xxxxx, g_yzz_xxxxy, g_yzz_xxxxz, g_yzz_xxxyy,\
                                     g_yzz_xxxyz, g_yzz_xxxzz, g_yzz_xxyyy, g_yzz_xxyyz,\
                                     g_yzz_xxyzz, g_yzz_xxzzz, g_yzz_xyyyy, g_yzz_xyyyz,\
                                     g_yzz_xyyzz, g_yzz_xyzzz, g_yzz_xzzzz, g_yzz_yyyyy,\
                                     g_yzz_yyyyz, g_yzz_yyyzz, g_yzz_yyzzz, g_yzz_yzzzz,\
                                     g_yzz_zzzzz, g_zzz_xxxxx, g_zzz_xxxxy, g_zzz_xxxxz,\
                                     g_zzz_xxxyy, g_zzz_xxxyz, g_zzz_xxxzz, g_zzz_xxyyy,\
                                     g_zzz_xxyyz, g_zzz_xxyzz, g_zzz_xxzzz, g_zzz_xyyyy,\
                                     g_zzz_xyyyz, g_zzz_xyyzz, g_zzz_xyzzz, g_zzz_xzzzz,\
                                     g_zzz_yyyyy, g_zzz_yyyyz, g_zzz_yyyzz, g_zzz_yyzzz,\
                                     g_zzz_yzzzz, g_zzz_zzzzz: VLX_ALIGN)
             for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xxx_xxxxx[j] = g1_xx_xxxxxx[j] - fr * g2_xx_xxxxx[j];

                g_xxx_xxxxy[j] = g1_xx_xxxxxy[j] - fr * g2_xx_xxxxy[j];

                g_xxx_xxxxz[j] = g1_xx_xxxxxz[j] - fr * g2_xx_xxxxz[j];

                g_xxx_xxxyy[j] = g1_xx_xxxxyy[j] - fr * g2_xx_xxxyy[j];

                g_xxx_xxxyz[j] = g1_xx_xxxxyz[j] - fr * g2_xx_xxxyz[j];

                g_xxx_xxxzz[j] = g1_xx_xxxxzz[j] - fr * g2_xx_xxxzz[j];

                g_xxx_xxyyy[j] = g1_xx_xxxyyy[j] - fr * g2_xx_xxyyy[j];

                g_xxx_xxyyz[j] = g1_xx_xxxyyz[j] - fr * g2_xx_xxyyz[j];

                g_xxx_xxyzz[j] = g1_xx_xxxyzz[j] - fr * g2_xx_xxyzz[j];

                g_xxx_xxzzz[j] = g1_xx_xxxzzz[j] - fr * g2_xx_xxzzz[j];

                g_xxx_xyyyy[j] = g1_xx_xxyyyy[j] - fr * g2_xx_xyyyy[j];

                g_xxx_xyyyz[j] = g1_xx_xxyyyz[j] - fr * g2_xx_xyyyz[j];

                g_xxx_xyyzz[j] = g1_xx_xxyyzz[j] - fr * g2_xx_xyyzz[j];

                g_xxx_xyzzz[j] = g1_xx_xxyzzz[j] - fr * g2_xx_xyzzz[j];

                g_xxx_xzzzz[j] = g1_xx_xxzzzz[j] - fr * g2_xx_xzzzz[j];

                g_xxx_yyyyy[j] = g1_xx_xyyyyy[j] - fr * g2_xx_yyyyy[j];

                g_xxx_yyyyz[j] = g1_xx_xyyyyz[j] - fr * g2_xx_yyyyz[j];

                g_xxx_yyyzz[j] = g1_xx_xyyyzz[j] - fr * g2_xx_yyyzz[j];

                g_xxx_yyzzz[j] = g1_xx_xyyzzz[j] - fr * g2_xx_yyzzz[j];

                g_xxx_yzzzz[j] = g1_xx_xyzzzz[j] - fr * g2_xx_yzzzz[j];

                g_xxx_zzzzz[j] = g1_xx_xzzzzz[j] - fr * g2_xx_zzzzz[j];

                g_xxy_xxxxx[j] = g1_xy_xxxxxx[j] - fr * g2_xy_xxxxx[j];

                g_xxy_xxxxy[j] = g1_xy_xxxxxy[j] - fr * g2_xy_xxxxy[j];

                g_xxy_xxxxz[j] = g1_xy_xxxxxz[j] - fr * g2_xy_xxxxz[j];

                g_xxy_xxxyy[j] = g1_xy_xxxxyy[j] - fr * g2_xy_xxxyy[j];

                g_xxy_xxxyz[j] = g1_xy_xxxxyz[j] - fr * g2_xy_xxxyz[j];

                g_xxy_xxxzz[j] = g1_xy_xxxxzz[j] - fr * g2_xy_xxxzz[j];

                g_xxy_xxyyy[j] = g1_xy_xxxyyy[j] - fr * g2_xy_xxyyy[j];

                g_xxy_xxyyz[j] = g1_xy_xxxyyz[j] - fr * g2_xy_xxyyz[j];

                g_xxy_xxyzz[j] = g1_xy_xxxyzz[j] - fr * g2_xy_xxyzz[j];

                g_xxy_xxzzz[j] = g1_xy_xxxzzz[j] - fr * g2_xy_xxzzz[j];

                g_xxy_xyyyy[j] = g1_xy_xxyyyy[j] - fr * g2_xy_xyyyy[j];

                g_xxy_xyyyz[j] = g1_xy_xxyyyz[j] - fr * g2_xy_xyyyz[j];

                g_xxy_xyyzz[j] = g1_xy_xxyyzz[j] - fr * g2_xy_xyyzz[j];

                g_xxy_xyzzz[j] = g1_xy_xxyzzz[j] - fr * g2_xy_xyzzz[j];

                g_xxy_xzzzz[j] = g1_xy_xxzzzz[j] - fr * g2_xy_xzzzz[j];

                g_xxy_yyyyy[j] = g1_xy_xyyyyy[j] - fr * g2_xy_yyyyy[j];

                g_xxy_yyyyz[j] = g1_xy_xyyyyz[j] - fr * g2_xy_yyyyz[j];

                g_xxy_yyyzz[j] = g1_xy_xyyyzz[j] - fr * g2_xy_yyyzz[j];

                g_xxy_yyzzz[j] = g1_xy_xyyzzz[j] - fr * g2_xy_yyzzz[j];

                g_xxy_yzzzz[j] = g1_xy_xyzzzz[j] - fr * g2_xy_yzzzz[j];

                g_xxy_zzzzz[j] = g1_xy_xzzzzz[j] - fr * g2_xy_zzzzz[j];

                g_xxz_xxxxx[j] = g1_xz_xxxxxx[j] - fr * g2_xz_xxxxx[j];

                g_xxz_xxxxy[j] = g1_xz_xxxxxy[j] - fr * g2_xz_xxxxy[j];

                g_xxz_xxxxz[j] = g1_xz_xxxxxz[j] - fr * g2_xz_xxxxz[j];

                g_xxz_xxxyy[j] = g1_xz_xxxxyy[j] - fr * g2_xz_xxxyy[j];

                g_xxz_xxxyz[j] = g1_xz_xxxxyz[j] - fr * g2_xz_xxxyz[j];

                g_xxz_xxxzz[j] = g1_xz_xxxxzz[j] - fr * g2_xz_xxxzz[j];

                g_xxz_xxyyy[j] = g1_xz_xxxyyy[j] - fr * g2_xz_xxyyy[j];

                g_xxz_xxyyz[j] = g1_xz_xxxyyz[j] - fr * g2_xz_xxyyz[j];

                g_xxz_xxyzz[j] = g1_xz_xxxyzz[j] - fr * g2_xz_xxyzz[j];

                g_xxz_xxzzz[j] = g1_xz_xxxzzz[j] - fr * g2_xz_xxzzz[j];

                g_xxz_xyyyy[j] = g1_xz_xxyyyy[j] - fr * g2_xz_xyyyy[j];

                g_xxz_xyyyz[j] = g1_xz_xxyyyz[j] - fr * g2_xz_xyyyz[j];

                g_xxz_xyyzz[j] = g1_xz_xxyyzz[j] - fr * g2_xz_xyyzz[j];

                g_xxz_xyzzz[j] = g1_xz_xxyzzz[j] - fr * g2_xz_xyzzz[j];

                g_xxz_xzzzz[j] = g1_xz_xxzzzz[j] - fr * g2_xz_xzzzz[j];

                g_xxz_yyyyy[j] = g1_xz_xyyyyy[j] - fr * g2_xz_yyyyy[j];

                g_xxz_yyyyz[j] = g1_xz_xyyyyz[j] - fr * g2_xz_yyyyz[j];

                g_xxz_yyyzz[j] = g1_xz_xyyyzz[j] - fr * g2_xz_yyyzz[j];

                g_xxz_yyzzz[j] = g1_xz_xyyzzz[j] - fr * g2_xz_yyzzz[j];

                g_xxz_yzzzz[j] = g1_xz_xyzzzz[j] - fr * g2_xz_yzzzz[j];

                g_xxz_zzzzz[j] = g1_xz_xzzzzz[j] - fr * g2_xz_zzzzz[j];

                g_xyy_xxxxx[j] = g1_yy_xxxxxx[j] - fr * g2_yy_xxxxx[j];

                g_xyy_xxxxy[j] = g1_yy_xxxxxy[j] - fr * g2_yy_xxxxy[j];

                g_xyy_xxxxz[j] = g1_yy_xxxxxz[j] - fr * g2_yy_xxxxz[j];

                g_xyy_xxxyy[j] = g1_yy_xxxxyy[j] - fr * g2_yy_xxxyy[j];

                g_xyy_xxxyz[j] = g1_yy_xxxxyz[j] - fr * g2_yy_xxxyz[j];

                g_xyy_xxxzz[j] = g1_yy_xxxxzz[j] - fr * g2_yy_xxxzz[j];

                g_xyy_xxyyy[j] = g1_yy_xxxyyy[j] - fr * g2_yy_xxyyy[j];

                g_xyy_xxyyz[j] = g1_yy_xxxyyz[j] - fr * g2_yy_xxyyz[j];

                g_xyy_xxyzz[j] = g1_yy_xxxyzz[j] - fr * g2_yy_xxyzz[j];

                g_xyy_xxzzz[j] = g1_yy_xxxzzz[j] - fr * g2_yy_xxzzz[j];

                g_xyy_xyyyy[j] = g1_yy_xxyyyy[j] - fr * g2_yy_xyyyy[j];

                g_xyy_xyyyz[j] = g1_yy_xxyyyz[j] - fr * g2_yy_xyyyz[j];

                g_xyy_xyyzz[j] = g1_yy_xxyyzz[j] - fr * g2_yy_xyyzz[j];

                g_xyy_xyzzz[j] = g1_yy_xxyzzz[j] - fr * g2_yy_xyzzz[j];

                g_xyy_xzzzz[j] = g1_yy_xxzzzz[j] - fr * g2_yy_xzzzz[j];

                g_xyy_yyyyy[j] = g1_yy_xyyyyy[j] - fr * g2_yy_yyyyy[j];

                g_xyy_yyyyz[j] = g1_yy_xyyyyz[j] - fr * g2_yy_yyyyz[j];

                g_xyy_yyyzz[j] = g1_yy_xyyyzz[j] - fr * g2_yy_yyyzz[j];

                g_xyy_yyzzz[j] = g1_yy_xyyzzz[j] - fr * g2_yy_yyzzz[j];

                g_xyy_yzzzz[j] = g1_yy_xyzzzz[j] - fr * g2_yy_yzzzz[j];

                g_xyy_zzzzz[j] = g1_yy_xzzzzz[j] - fr * g2_yy_zzzzz[j];

                g_xyz_xxxxx[j] = g1_yz_xxxxxx[j] - fr * g2_yz_xxxxx[j];

                g_xyz_xxxxy[j] = g1_yz_xxxxxy[j] - fr * g2_yz_xxxxy[j];

                g_xyz_xxxxz[j] = g1_yz_xxxxxz[j] - fr * g2_yz_xxxxz[j];

                g_xyz_xxxyy[j] = g1_yz_xxxxyy[j] - fr * g2_yz_xxxyy[j];

                g_xyz_xxxyz[j] = g1_yz_xxxxyz[j] - fr * g2_yz_xxxyz[j];

                g_xyz_xxxzz[j] = g1_yz_xxxxzz[j] - fr * g2_yz_xxxzz[j];

                g_xyz_xxyyy[j] = g1_yz_xxxyyy[j] - fr * g2_yz_xxyyy[j];

                g_xyz_xxyyz[j] = g1_yz_xxxyyz[j] - fr * g2_yz_xxyyz[j];

                g_xyz_xxyzz[j] = g1_yz_xxxyzz[j] - fr * g2_yz_xxyzz[j];

                g_xyz_xxzzz[j] = g1_yz_xxxzzz[j] - fr * g2_yz_xxzzz[j];

                g_xyz_xyyyy[j] = g1_yz_xxyyyy[j] - fr * g2_yz_xyyyy[j];

                g_xyz_xyyyz[j] = g1_yz_xxyyyz[j] - fr * g2_yz_xyyyz[j];

                g_xyz_xyyzz[j] = g1_yz_xxyyzz[j] - fr * g2_yz_xyyzz[j];

                g_xyz_xyzzz[j] = g1_yz_xxyzzz[j] - fr * g2_yz_xyzzz[j];

                g_xyz_xzzzz[j] = g1_yz_xxzzzz[j] - fr * g2_yz_xzzzz[j];

                g_xyz_yyyyy[j] = g1_yz_xyyyyy[j] - fr * g2_yz_yyyyy[j];

                g_xyz_yyyyz[j] = g1_yz_xyyyyz[j] - fr * g2_yz_yyyyz[j];

                g_xyz_yyyzz[j] = g1_yz_xyyyzz[j] - fr * g2_yz_yyyzz[j];

                g_xyz_yyzzz[j] = g1_yz_xyyzzz[j] - fr * g2_yz_yyzzz[j];

                g_xyz_yzzzz[j] = g1_yz_xyzzzz[j] - fr * g2_yz_yzzzz[j];

                g_xyz_zzzzz[j] = g1_yz_xzzzzz[j] - fr * g2_yz_zzzzz[j];

                g_xzz_xxxxx[j] = g1_zz_xxxxxx[j] - fr * g2_zz_xxxxx[j];

                g_xzz_xxxxy[j] = g1_zz_xxxxxy[j] - fr * g2_zz_xxxxy[j];

                g_xzz_xxxxz[j] = g1_zz_xxxxxz[j] - fr * g2_zz_xxxxz[j];

                g_xzz_xxxyy[j] = g1_zz_xxxxyy[j] - fr * g2_zz_xxxyy[j];

                g_xzz_xxxyz[j] = g1_zz_xxxxyz[j] - fr * g2_zz_xxxyz[j];

                g_xzz_xxxzz[j] = g1_zz_xxxxzz[j] - fr * g2_zz_xxxzz[j];

                g_xzz_xxyyy[j] = g1_zz_xxxyyy[j] - fr * g2_zz_xxyyy[j];

                g_xzz_xxyyz[j] = g1_zz_xxxyyz[j] - fr * g2_zz_xxyyz[j];

                g_xzz_xxyzz[j] = g1_zz_xxxyzz[j] - fr * g2_zz_xxyzz[j];

                g_xzz_xxzzz[j] = g1_zz_xxxzzz[j] - fr * g2_zz_xxzzz[j];

                g_xzz_xyyyy[j] = g1_zz_xxyyyy[j] - fr * g2_zz_xyyyy[j];

                g_xzz_xyyyz[j] = g1_zz_xxyyyz[j] - fr * g2_zz_xyyyz[j];

                g_xzz_xyyzz[j] = g1_zz_xxyyzz[j] - fr * g2_zz_xyyzz[j];

                g_xzz_xyzzz[j] = g1_zz_xxyzzz[j] - fr * g2_zz_xyzzz[j];

                g_xzz_xzzzz[j] = g1_zz_xxzzzz[j] - fr * g2_zz_xzzzz[j];

                g_xzz_yyyyy[j] = g1_zz_xyyyyy[j] - fr * g2_zz_yyyyy[j];

                g_xzz_yyyyz[j] = g1_zz_xyyyyz[j] - fr * g2_zz_yyyyz[j];

                g_xzz_yyyzz[j] = g1_zz_xyyyzz[j] - fr * g2_zz_yyyzz[j];

                g_xzz_yyzzz[j] = g1_zz_xyyzzz[j] - fr * g2_zz_yyzzz[j];

                g_xzz_yzzzz[j] = g1_zz_xyzzzz[j] - fr * g2_zz_yzzzz[j];

                g_xzz_zzzzz[j] = g1_zz_xzzzzz[j] - fr * g2_zz_zzzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yyy_xxxxx[j] = g1_yy_xxxxxy[j] - fr * g2_yy_xxxxx[j];

                g_yyy_xxxxy[j] = g1_yy_xxxxyy[j] - fr * g2_yy_xxxxy[j];

                g_yyy_xxxxz[j] = g1_yy_xxxxyz[j] - fr * g2_yy_xxxxz[j];

                g_yyy_xxxyy[j] = g1_yy_xxxyyy[j] - fr * g2_yy_xxxyy[j];

                g_yyy_xxxyz[j] = g1_yy_xxxyyz[j] - fr * g2_yy_xxxyz[j];

                g_yyy_xxxzz[j] = g1_yy_xxxyzz[j] - fr * g2_yy_xxxzz[j];

                g_yyy_xxyyy[j] = g1_yy_xxyyyy[j] - fr * g2_yy_xxyyy[j];

                g_yyy_xxyyz[j] = g1_yy_xxyyyz[j] - fr * g2_yy_xxyyz[j];

                g_yyy_xxyzz[j] = g1_yy_xxyyzz[j] - fr * g2_yy_xxyzz[j];

                g_yyy_xxzzz[j] = g1_yy_xxyzzz[j] - fr * g2_yy_xxzzz[j];

                g_yyy_xyyyy[j] = g1_yy_xyyyyy[j] - fr * g2_yy_xyyyy[j];

                g_yyy_xyyyz[j] = g1_yy_xyyyyz[j] - fr * g2_yy_xyyyz[j];

                g_yyy_xyyzz[j] = g1_yy_xyyyzz[j] - fr * g2_yy_xyyzz[j];

                g_yyy_xyzzz[j] = g1_yy_xyyzzz[j] - fr * g2_yy_xyzzz[j];

                g_yyy_xzzzz[j] = g1_yy_xyzzzz[j] - fr * g2_yy_xzzzz[j];

                g_yyy_yyyyy[j] = g1_yy_yyyyyy[j] - fr * g2_yy_yyyyy[j];

                g_yyy_yyyyz[j] = g1_yy_yyyyyz[j] - fr * g2_yy_yyyyz[j];

                g_yyy_yyyzz[j] = g1_yy_yyyyzz[j] - fr * g2_yy_yyyzz[j];

                g_yyy_yyzzz[j] = g1_yy_yyyzzz[j] - fr * g2_yy_yyzzz[j];

                g_yyy_yzzzz[j] = g1_yy_yyzzzz[j] - fr * g2_yy_yzzzz[j];

                g_yyy_zzzzz[j] = g1_yy_yzzzzz[j] - fr * g2_yy_zzzzz[j];

                g_yyz_xxxxx[j] = g1_yz_xxxxxy[j] - fr * g2_yz_xxxxx[j];

                g_yyz_xxxxy[j] = g1_yz_xxxxyy[j] - fr * g2_yz_xxxxy[j];

                g_yyz_xxxxz[j] = g1_yz_xxxxyz[j] - fr * g2_yz_xxxxz[j];

                g_yyz_xxxyy[j] = g1_yz_xxxyyy[j] - fr * g2_yz_xxxyy[j];

                g_yyz_xxxyz[j] = g1_yz_xxxyyz[j] - fr * g2_yz_xxxyz[j];

                g_yyz_xxxzz[j] = g1_yz_xxxyzz[j] - fr * g2_yz_xxxzz[j];

                g_yyz_xxyyy[j] = g1_yz_xxyyyy[j] - fr * g2_yz_xxyyy[j];

                g_yyz_xxyyz[j] = g1_yz_xxyyyz[j] - fr * g2_yz_xxyyz[j];

                g_yyz_xxyzz[j] = g1_yz_xxyyzz[j] - fr * g2_yz_xxyzz[j];

                g_yyz_xxzzz[j] = g1_yz_xxyzzz[j] - fr * g2_yz_xxzzz[j];

                g_yyz_xyyyy[j] = g1_yz_xyyyyy[j] - fr * g2_yz_xyyyy[j];

                g_yyz_xyyyz[j] = g1_yz_xyyyyz[j] - fr * g2_yz_xyyyz[j];

                g_yyz_xyyzz[j] = g1_yz_xyyyzz[j] - fr * g2_yz_xyyzz[j];

                g_yyz_xyzzz[j] = g1_yz_xyyzzz[j] - fr * g2_yz_xyzzz[j];

                g_yyz_xzzzz[j] = g1_yz_xyzzzz[j] - fr * g2_yz_xzzzz[j];

                g_yyz_yyyyy[j] = g1_yz_yyyyyy[j] - fr * g2_yz_yyyyy[j];

                g_yyz_yyyyz[j] = g1_yz_yyyyyz[j] - fr * g2_yz_yyyyz[j];

                g_yyz_yyyzz[j] = g1_yz_yyyyzz[j] - fr * g2_yz_yyyzz[j];

                g_yyz_yyzzz[j] = g1_yz_yyyzzz[j] - fr * g2_yz_yyzzz[j];

                g_yyz_yzzzz[j] = g1_yz_yyzzzz[j] - fr * g2_yz_yzzzz[j];

                g_yyz_zzzzz[j] = g1_yz_yzzzzz[j] - fr * g2_yz_zzzzz[j];

                g_yzz_xxxxx[j] = g1_zz_xxxxxy[j] - fr * g2_zz_xxxxx[j];

                g_yzz_xxxxy[j] = g1_zz_xxxxyy[j] - fr * g2_zz_xxxxy[j];

                g_yzz_xxxxz[j] = g1_zz_xxxxyz[j] - fr * g2_zz_xxxxz[j];

                g_yzz_xxxyy[j] = g1_zz_xxxyyy[j] - fr * g2_zz_xxxyy[j];

                g_yzz_xxxyz[j] = g1_zz_xxxyyz[j] - fr * g2_zz_xxxyz[j];

                g_yzz_xxxzz[j] = g1_zz_xxxyzz[j] - fr * g2_zz_xxxzz[j];

                g_yzz_xxyyy[j] = g1_zz_xxyyyy[j] - fr * g2_zz_xxyyy[j];

                g_yzz_xxyyz[j] = g1_zz_xxyyyz[j] - fr * g2_zz_xxyyz[j];

                g_yzz_xxyzz[j] = g1_zz_xxyyzz[j] - fr * g2_zz_xxyzz[j];

                g_yzz_xxzzz[j] = g1_zz_xxyzzz[j] - fr * g2_zz_xxzzz[j];

                g_yzz_xyyyy[j] = g1_zz_xyyyyy[j] - fr * g2_zz_xyyyy[j];

                g_yzz_xyyyz[j] = g1_zz_xyyyyz[j] - fr * g2_zz_xyyyz[j];

                g_yzz_xyyzz[j] = g1_zz_xyyyzz[j] - fr * g2_zz_xyyzz[j];

                g_yzz_xyzzz[j] = g1_zz_xyyzzz[j] - fr * g2_zz_xyzzz[j];

                g_yzz_xzzzz[j] = g1_zz_xyzzzz[j] - fr * g2_zz_xzzzz[j];

                g_yzz_yyyyy[j] = g1_zz_yyyyyy[j] - fr * g2_zz_yyyyy[j];

                g_yzz_yyyyz[j] = g1_zz_yyyyyz[j] - fr * g2_zz_yyyyz[j];

                g_yzz_yyyzz[j] = g1_zz_yyyyzz[j] - fr * g2_zz_yyyzz[j];

                g_yzz_yyzzz[j] = g1_zz_yyyzzz[j] - fr * g2_zz_yyzzz[j];

                g_yzz_yzzzz[j] = g1_zz_yyzzzz[j] - fr * g2_zz_yzzzz[j];

                g_yzz_zzzzz[j] = g1_zz_yzzzzz[j] - fr * g2_zz_zzzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zzz_xxxxx[j] = g1_zz_xxxxxz[j] - fr * g2_zz_xxxxx[j];

                g_zzz_xxxxy[j] = g1_zz_xxxxyz[j] - fr * g2_zz_xxxxy[j];

                g_zzz_xxxxz[j] = g1_zz_xxxxzz[j] - fr * g2_zz_xxxxz[j];

                g_zzz_xxxyy[j] = g1_zz_xxxyyz[j] - fr * g2_zz_xxxyy[j];

                g_zzz_xxxyz[j] = g1_zz_xxxyzz[j] - fr * g2_zz_xxxyz[j];

                g_zzz_xxxzz[j] = g1_zz_xxxzzz[j] - fr * g2_zz_xxxzz[j];

                g_zzz_xxyyy[j] = g1_zz_xxyyyz[j] - fr * g2_zz_xxyyy[j];

                g_zzz_xxyyz[j] = g1_zz_xxyyzz[j] - fr * g2_zz_xxyyz[j];

                g_zzz_xxyzz[j] = g1_zz_xxyzzz[j] - fr * g2_zz_xxyzz[j];

                g_zzz_xxzzz[j] = g1_zz_xxzzzz[j] - fr * g2_zz_xxzzz[j];

                g_zzz_xyyyy[j] = g1_zz_xyyyyz[j] - fr * g2_zz_xyyyy[j];

                g_zzz_xyyyz[j] = g1_zz_xyyyzz[j] - fr * g2_zz_xyyyz[j];

                g_zzz_xyyzz[j] = g1_zz_xyyzzz[j] - fr * g2_zz_xyyzz[j];

                g_zzz_xyzzz[j] = g1_zz_xyzzzz[j] - fr * g2_zz_xyzzz[j];

                g_zzz_xzzzz[j] = g1_zz_xzzzzz[j] - fr * g2_zz_xzzzz[j];

                g_zzz_yyyyy[j] = g1_zz_yyyyyz[j] - fr * g2_zz_yyyyy[j];

                g_zzz_yyyyz[j] = g1_zz_yyyyzz[j] - fr * g2_zz_yyyyz[j];

                g_zzz_yyyzz[j] = g1_zz_yyyzzz[j] - fr * g2_zz_yyyzz[j];

                g_zzz_yyzzz[j] = g1_zz_yyzzzz[j] - fr * g2_zz_yyzzz[j];

                g_zzz_yzzzz[j] = g1_zz_yzzzzz[j] - fr * g2_zz_yzzzz[j];

                g_zzz_zzzzz[j] = g1_zz_zzzzzz[j] - fr * g2_zz_zzzzz[j];
            }
        }
    }
    
    void
    compElectronRepulsionForXGG(      CMemBlock2D<double>&  contrBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  cdDistances,
                                const int32_t               braAngularMomentum,
                                const CGtoPairsBlock&       ketGtoPairsBlock)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {braAngularMomentum, 4, 4})) return;

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
                                             {braAngularMomentum, 4, 4});

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 3, 5});

        auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                              {braAngularMomentum, 3, 4});

        // compute contracted integrals

        for (int32_t i = 0; i < bcomp; i++)
        {
            // set up pointers to (X|g(r,r')|FG)^(m) integrals

            auto g2_xxx_xxxx = contrBuffer.data(g2off + 150 * i);

            auto g2_xxx_xxxy = contrBuffer.data(g2off + 150 * i + 1);

            auto g2_xxx_xxxz = contrBuffer.data(g2off + 150 * i + 2);

            auto g2_xxx_xxyy = contrBuffer.data(g2off + 150 * i + 3);

            auto g2_xxx_xxyz = contrBuffer.data(g2off + 150 * i + 4);

            auto g2_xxx_xxzz = contrBuffer.data(g2off + 150 * i + 5);

            auto g2_xxx_xyyy = contrBuffer.data(g2off + 150 * i + 6);

            auto g2_xxx_xyyz = contrBuffer.data(g2off + 150 * i + 7);

            auto g2_xxx_xyzz = contrBuffer.data(g2off + 150 * i + 8);

            auto g2_xxx_xzzz = contrBuffer.data(g2off + 150 * i + 9);

            auto g2_xxx_yyyy = contrBuffer.data(g2off + 150 * i + 10);

            auto g2_xxx_yyyz = contrBuffer.data(g2off + 150 * i + 11);

            auto g2_xxx_yyzz = contrBuffer.data(g2off + 150 * i + 12);

            auto g2_xxx_yzzz = contrBuffer.data(g2off + 150 * i + 13);

            auto g2_xxx_zzzz = contrBuffer.data(g2off + 150 * i + 14);

            auto g2_xxy_xxxx = contrBuffer.data(g2off + 150 * i + 15);

            auto g2_xxy_xxxy = contrBuffer.data(g2off + 150 * i + 16);

            auto g2_xxy_xxxz = contrBuffer.data(g2off + 150 * i + 17);

            auto g2_xxy_xxyy = contrBuffer.data(g2off + 150 * i + 18);

            auto g2_xxy_xxyz = contrBuffer.data(g2off + 150 * i + 19);

            auto g2_xxy_xxzz = contrBuffer.data(g2off + 150 * i + 20);

            auto g2_xxy_xyyy = contrBuffer.data(g2off + 150 * i + 21);

            auto g2_xxy_xyyz = contrBuffer.data(g2off + 150 * i + 22);

            auto g2_xxy_xyzz = contrBuffer.data(g2off + 150 * i + 23);

            auto g2_xxy_xzzz = contrBuffer.data(g2off + 150 * i + 24);

            auto g2_xxy_yyyy = contrBuffer.data(g2off + 150 * i + 25);

            auto g2_xxy_yyyz = contrBuffer.data(g2off + 150 * i + 26);

            auto g2_xxy_yyzz = contrBuffer.data(g2off + 150 * i + 27);

            auto g2_xxy_yzzz = contrBuffer.data(g2off + 150 * i + 28);

            auto g2_xxy_zzzz = contrBuffer.data(g2off + 150 * i + 29);

            auto g2_xxz_xxxx = contrBuffer.data(g2off + 150 * i + 30);

            auto g2_xxz_xxxy = contrBuffer.data(g2off + 150 * i + 31);

            auto g2_xxz_xxxz = contrBuffer.data(g2off + 150 * i + 32);

            auto g2_xxz_xxyy = contrBuffer.data(g2off + 150 * i + 33);

            auto g2_xxz_xxyz = contrBuffer.data(g2off + 150 * i + 34);

            auto g2_xxz_xxzz = contrBuffer.data(g2off + 150 * i + 35);

            auto g2_xxz_xyyy = contrBuffer.data(g2off + 150 * i + 36);

            auto g2_xxz_xyyz = contrBuffer.data(g2off + 150 * i + 37);

            auto g2_xxz_xyzz = contrBuffer.data(g2off + 150 * i + 38);

            auto g2_xxz_xzzz = contrBuffer.data(g2off + 150 * i + 39);

            auto g2_xxz_yyyy = contrBuffer.data(g2off + 150 * i + 40);

            auto g2_xxz_yyyz = contrBuffer.data(g2off + 150 * i + 41);

            auto g2_xxz_yyzz = contrBuffer.data(g2off + 150 * i + 42);

            auto g2_xxz_yzzz = contrBuffer.data(g2off + 150 * i + 43);

            auto g2_xxz_zzzz = contrBuffer.data(g2off + 150 * i + 44);

            auto g2_xyy_xxxx = contrBuffer.data(g2off + 150 * i + 45);

            auto g2_xyy_xxxy = contrBuffer.data(g2off + 150 * i + 46);

            auto g2_xyy_xxxz = contrBuffer.data(g2off + 150 * i + 47);

            auto g2_xyy_xxyy = contrBuffer.data(g2off + 150 * i + 48);

            auto g2_xyy_xxyz = contrBuffer.data(g2off + 150 * i + 49);

            auto g2_xyy_xxzz = contrBuffer.data(g2off + 150 * i + 50);

            auto g2_xyy_xyyy = contrBuffer.data(g2off + 150 * i + 51);

            auto g2_xyy_xyyz = contrBuffer.data(g2off + 150 * i + 52);

            auto g2_xyy_xyzz = contrBuffer.data(g2off + 150 * i + 53);

            auto g2_xyy_xzzz = contrBuffer.data(g2off + 150 * i + 54);

            auto g2_xyy_yyyy = contrBuffer.data(g2off + 150 * i + 55);

            auto g2_xyy_yyyz = contrBuffer.data(g2off + 150 * i + 56);

            auto g2_xyy_yyzz = contrBuffer.data(g2off + 150 * i + 57);

            auto g2_xyy_yzzz = contrBuffer.data(g2off + 150 * i + 58);

            auto g2_xyy_zzzz = contrBuffer.data(g2off + 150 * i + 59);

            auto g2_xyz_xxxx = contrBuffer.data(g2off + 150 * i + 60);

            auto g2_xyz_xxxy = contrBuffer.data(g2off + 150 * i + 61);

            auto g2_xyz_xxxz = contrBuffer.data(g2off + 150 * i + 62);

            auto g2_xyz_xxyy = contrBuffer.data(g2off + 150 * i + 63);

            auto g2_xyz_xxyz = contrBuffer.data(g2off + 150 * i + 64);

            auto g2_xyz_xxzz = contrBuffer.data(g2off + 150 * i + 65);

            auto g2_xyz_xyyy = contrBuffer.data(g2off + 150 * i + 66);

            auto g2_xyz_xyyz = contrBuffer.data(g2off + 150 * i + 67);

            auto g2_xyz_xyzz = contrBuffer.data(g2off + 150 * i + 68);

            auto g2_xyz_xzzz = contrBuffer.data(g2off + 150 * i + 69);

            auto g2_xyz_yyyy = contrBuffer.data(g2off + 150 * i + 70);

            auto g2_xyz_yyyz = contrBuffer.data(g2off + 150 * i + 71);

            auto g2_xyz_yyzz = contrBuffer.data(g2off + 150 * i + 72);

            auto g2_xyz_yzzz = contrBuffer.data(g2off + 150 * i + 73);

            auto g2_xyz_zzzz = contrBuffer.data(g2off + 150 * i + 74);

            auto g2_xzz_xxxx = contrBuffer.data(g2off + 150 * i + 75);

            auto g2_xzz_xxxy = contrBuffer.data(g2off + 150 * i + 76);

            auto g2_xzz_xxxz = contrBuffer.data(g2off + 150 * i + 77);

            auto g2_xzz_xxyy = contrBuffer.data(g2off + 150 * i + 78);

            auto g2_xzz_xxyz = contrBuffer.data(g2off + 150 * i + 79);

            auto g2_xzz_xxzz = contrBuffer.data(g2off + 150 * i + 80);

            auto g2_xzz_xyyy = contrBuffer.data(g2off + 150 * i + 81);

            auto g2_xzz_xyyz = contrBuffer.data(g2off + 150 * i + 82);

            auto g2_xzz_xyzz = contrBuffer.data(g2off + 150 * i + 83);

            auto g2_xzz_xzzz = contrBuffer.data(g2off + 150 * i + 84);

            auto g2_xzz_yyyy = contrBuffer.data(g2off + 150 * i + 85);

            auto g2_xzz_yyyz = contrBuffer.data(g2off + 150 * i + 86);

            auto g2_xzz_yyzz = contrBuffer.data(g2off + 150 * i + 87);

            auto g2_xzz_yzzz = contrBuffer.data(g2off + 150 * i + 88);

            auto g2_xzz_zzzz = contrBuffer.data(g2off + 150 * i + 89);

            auto g2_yyy_xxxx = contrBuffer.data(g2off + 150 * i + 90);

            auto g2_yyy_xxxy = contrBuffer.data(g2off + 150 * i + 91);

            auto g2_yyy_xxxz = contrBuffer.data(g2off + 150 * i + 92);

            auto g2_yyy_xxyy = contrBuffer.data(g2off + 150 * i + 93);

            auto g2_yyy_xxyz = contrBuffer.data(g2off + 150 * i + 94);

            auto g2_yyy_xxzz = contrBuffer.data(g2off + 150 * i + 95);

            auto g2_yyy_xyyy = contrBuffer.data(g2off + 150 * i + 96);

            auto g2_yyy_xyyz = contrBuffer.data(g2off + 150 * i + 97);

            auto g2_yyy_xyzz = contrBuffer.data(g2off + 150 * i + 98);

            auto g2_yyy_xzzz = contrBuffer.data(g2off + 150 * i + 99);

            auto g2_yyy_yyyy = contrBuffer.data(g2off + 150 * i + 100);

            auto g2_yyy_yyyz = contrBuffer.data(g2off + 150 * i + 101);

            auto g2_yyy_yyzz = contrBuffer.data(g2off + 150 * i + 102);

            auto g2_yyy_yzzz = contrBuffer.data(g2off + 150 * i + 103);

            auto g2_yyy_zzzz = contrBuffer.data(g2off + 150 * i + 104);

            auto g2_yyz_xxxx = contrBuffer.data(g2off + 150 * i + 105);

            auto g2_yyz_xxxy = contrBuffer.data(g2off + 150 * i + 106);

            auto g2_yyz_xxxz = contrBuffer.data(g2off + 150 * i + 107);

            auto g2_yyz_xxyy = contrBuffer.data(g2off + 150 * i + 108);

            auto g2_yyz_xxyz = contrBuffer.data(g2off + 150 * i + 109);

            auto g2_yyz_xxzz = contrBuffer.data(g2off + 150 * i + 110);

            auto g2_yyz_xyyy = contrBuffer.data(g2off + 150 * i + 111);

            auto g2_yyz_xyyz = contrBuffer.data(g2off + 150 * i + 112);

            auto g2_yyz_xyzz = contrBuffer.data(g2off + 150 * i + 113);

            auto g2_yyz_xzzz = contrBuffer.data(g2off + 150 * i + 114);

            auto g2_yyz_yyyy = contrBuffer.data(g2off + 150 * i + 115);

            auto g2_yyz_yyyz = contrBuffer.data(g2off + 150 * i + 116);

            auto g2_yyz_yyzz = contrBuffer.data(g2off + 150 * i + 117);

            auto g2_yyz_yzzz = contrBuffer.data(g2off + 150 * i + 118);

            auto g2_yyz_zzzz = contrBuffer.data(g2off + 150 * i + 119);

            auto g2_yzz_xxxx = contrBuffer.data(g2off + 150 * i + 120);

            auto g2_yzz_xxxy = contrBuffer.data(g2off + 150 * i + 121);

            auto g2_yzz_xxxz = contrBuffer.data(g2off + 150 * i + 122);

            auto g2_yzz_xxyy = contrBuffer.data(g2off + 150 * i + 123);

            auto g2_yzz_xxyz = contrBuffer.data(g2off + 150 * i + 124);

            auto g2_yzz_xxzz = contrBuffer.data(g2off + 150 * i + 125);

            auto g2_yzz_xyyy = contrBuffer.data(g2off + 150 * i + 126);

            auto g2_yzz_xyyz = contrBuffer.data(g2off + 150 * i + 127);

            auto g2_yzz_xyzz = contrBuffer.data(g2off + 150 * i + 128);

            auto g2_yzz_xzzz = contrBuffer.data(g2off + 150 * i + 129);

            auto g2_yzz_yyyy = contrBuffer.data(g2off + 150 * i + 130);

            auto g2_yzz_yyyz = contrBuffer.data(g2off + 150 * i + 131);

            auto g2_yzz_yyzz = contrBuffer.data(g2off + 150 * i + 132);

            auto g2_yzz_yzzz = contrBuffer.data(g2off + 150 * i + 133);

            auto g2_yzz_zzzz = contrBuffer.data(g2off + 150 * i + 134);

            auto g2_zzz_xxxx = contrBuffer.data(g2off + 150 * i + 135);

            auto g2_zzz_xxxy = contrBuffer.data(g2off + 150 * i + 136);

            auto g2_zzz_xxxz = contrBuffer.data(g2off + 150 * i + 137);

            auto g2_zzz_xxyy = contrBuffer.data(g2off + 150 * i + 138);

            auto g2_zzz_xxyz = contrBuffer.data(g2off + 150 * i + 139);

            auto g2_zzz_xxzz = contrBuffer.data(g2off + 150 * i + 140);

            auto g2_zzz_xyyy = contrBuffer.data(g2off + 150 * i + 141);

            auto g2_zzz_xyyz = contrBuffer.data(g2off + 150 * i + 142);

            auto g2_zzz_xyzz = contrBuffer.data(g2off + 150 * i + 143);

            auto g2_zzz_xzzz = contrBuffer.data(g2off + 150 * i + 144);

            auto g2_zzz_yyyy = contrBuffer.data(g2off + 150 * i + 145);

            auto g2_zzz_yyyz = contrBuffer.data(g2off + 150 * i + 146);

            auto g2_zzz_yyzz = contrBuffer.data(g2off + 150 * i + 147);

            auto g2_zzz_yzzz = contrBuffer.data(g2off + 150 * i + 148);

            auto g2_zzz_zzzz = contrBuffer.data(g2off + 150 * i + 149);

            // set up pointers to (X|g(r,r')|FH)^(m) integrals

            auto g1_xxx_xxxxx = contrBuffer.data(g1off + 210 * i);

            auto g1_xxx_xxxxy = contrBuffer.data(g1off + 210 * i + 1);

            auto g1_xxx_xxxxz = contrBuffer.data(g1off + 210 * i + 2);

            auto g1_xxx_xxxyy = contrBuffer.data(g1off + 210 * i + 3);

            auto g1_xxx_xxxyz = contrBuffer.data(g1off + 210 * i + 4);

            auto g1_xxx_xxxzz = contrBuffer.data(g1off + 210 * i + 5);

            auto g1_xxx_xxyyy = contrBuffer.data(g1off + 210 * i + 6);

            auto g1_xxx_xxyyz = contrBuffer.data(g1off + 210 * i + 7);

            auto g1_xxx_xxyzz = contrBuffer.data(g1off + 210 * i + 8);

            auto g1_xxx_xxzzz = contrBuffer.data(g1off + 210 * i + 9);

            auto g1_xxx_xyyyy = contrBuffer.data(g1off + 210 * i + 10);

            auto g1_xxx_xyyyz = contrBuffer.data(g1off + 210 * i + 11);

            auto g1_xxx_xyyzz = contrBuffer.data(g1off + 210 * i + 12);

            auto g1_xxx_xyzzz = contrBuffer.data(g1off + 210 * i + 13);

            auto g1_xxx_xzzzz = contrBuffer.data(g1off + 210 * i + 14);

            auto g1_xxy_xxxxx = contrBuffer.data(g1off + 210 * i + 21);

            auto g1_xxy_xxxxy = contrBuffer.data(g1off + 210 * i + 22);

            auto g1_xxy_xxxxz = contrBuffer.data(g1off + 210 * i + 23);

            auto g1_xxy_xxxyy = contrBuffer.data(g1off + 210 * i + 24);

            auto g1_xxy_xxxyz = contrBuffer.data(g1off + 210 * i + 25);

            auto g1_xxy_xxxzz = contrBuffer.data(g1off + 210 * i + 26);

            auto g1_xxy_xxyyy = contrBuffer.data(g1off + 210 * i + 27);

            auto g1_xxy_xxyyz = contrBuffer.data(g1off + 210 * i + 28);

            auto g1_xxy_xxyzz = contrBuffer.data(g1off + 210 * i + 29);

            auto g1_xxy_xxzzz = contrBuffer.data(g1off + 210 * i + 30);

            auto g1_xxy_xyyyy = contrBuffer.data(g1off + 210 * i + 31);

            auto g1_xxy_xyyyz = contrBuffer.data(g1off + 210 * i + 32);

            auto g1_xxy_xyyzz = contrBuffer.data(g1off + 210 * i + 33);

            auto g1_xxy_xyzzz = contrBuffer.data(g1off + 210 * i + 34);

            auto g1_xxy_xzzzz = contrBuffer.data(g1off + 210 * i + 35);

            auto g1_xxz_xxxxx = contrBuffer.data(g1off + 210 * i + 42);

            auto g1_xxz_xxxxy = contrBuffer.data(g1off + 210 * i + 43);

            auto g1_xxz_xxxxz = contrBuffer.data(g1off + 210 * i + 44);

            auto g1_xxz_xxxyy = contrBuffer.data(g1off + 210 * i + 45);

            auto g1_xxz_xxxyz = contrBuffer.data(g1off + 210 * i + 46);

            auto g1_xxz_xxxzz = contrBuffer.data(g1off + 210 * i + 47);

            auto g1_xxz_xxyyy = contrBuffer.data(g1off + 210 * i + 48);

            auto g1_xxz_xxyyz = contrBuffer.data(g1off + 210 * i + 49);

            auto g1_xxz_xxyzz = contrBuffer.data(g1off + 210 * i + 50);

            auto g1_xxz_xxzzz = contrBuffer.data(g1off + 210 * i + 51);

            auto g1_xxz_xyyyy = contrBuffer.data(g1off + 210 * i + 52);

            auto g1_xxz_xyyyz = contrBuffer.data(g1off + 210 * i + 53);

            auto g1_xxz_xyyzz = contrBuffer.data(g1off + 210 * i + 54);

            auto g1_xxz_xyzzz = contrBuffer.data(g1off + 210 * i + 55);

            auto g1_xxz_xzzzz = contrBuffer.data(g1off + 210 * i + 56);

            auto g1_xyy_xxxxx = contrBuffer.data(g1off + 210 * i + 63);

            auto g1_xyy_xxxxy = contrBuffer.data(g1off + 210 * i + 64);

            auto g1_xyy_xxxxz = contrBuffer.data(g1off + 210 * i + 65);

            auto g1_xyy_xxxyy = contrBuffer.data(g1off + 210 * i + 66);

            auto g1_xyy_xxxyz = contrBuffer.data(g1off + 210 * i + 67);

            auto g1_xyy_xxxzz = contrBuffer.data(g1off + 210 * i + 68);

            auto g1_xyy_xxyyy = contrBuffer.data(g1off + 210 * i + 69);

            auto g1_xyy_xxyyz = contrBuffer.data(g1off + 210 * i + 70);

            auto g1_xyy_xxyzz = contrBuffer.data(g1off + 210 * i + 71);

            auto g1_xyy_xxzzz = contrBuffer.data(g1off + 210 * i + 72);

            auto g1_xyy_xyyyy = contrBuffer.data(g1off + 210 * i + 73);

            auto g1_xyy_xyyyz = contrBuffer.data(g1off + 210 * i + 74);

            auto g1_xyy_xyyzz = contrBuffer.data(g1off + 210 * i + 75);

            auto g1_xyy_xyzzz = contrBuffer.data(g1off + 210 * i + 76);

            auto g1_xyy_xzzzz = contrBuffer.data(g1off + 210 * i + 77);

            auto g1_xyz_xxxxx = contrBuffer.data(g1off + 210 * i + 84);

            auto g1_xyz_xxxxy = contrBuffer.data(g1off + 210 * i + 85);

            auto g1_xyz_xxxxz = contrBuffer.data(g1off + 210 * i + 86);

            auto g1_xyz_xxxyy = contrBuffer.data(g1off + 210 * i + 87);

            auto g1_xyz_xxxyz = contrBuffer.data(g1off + 210 * i + 88);

            auto g1_xyz_xxxzz = contrBuffer.data(g1off + 210 * i + 89);

            auto g1_xyz_xxyyy = contrBuffer.data(g1off + 210 * i + 90);

            auto g1_xyz_xxyyz = contrBuffer.data(g1off + 210 * i + 91);

            auto g1_xyz_xxyzz = contrBuffer.data(g1off + 210 * i + 92);

            auto g1_xyz_xxzzz = contrBuffer.data(g1off + 210 * i + 93);

            auto g1_xyz_xyyyy = contrBuffer.data(g1off + 210 * i + 94);

            auto g1_xyz_xyyyz = contrBuffer.data(g1off + 210 * i + 95);

            auto g1_xyz_xyyzz = contrBuffer.data(g1off + 210 * i + 96);

            auto g1_xyz_xyzzz = contrBuffer.data(g1off + 210 * i + 97);

            auto g1_xyz_xzzzz = contrBuffer.data(g1off + 210 * i + 98);

            auto g1_xzz_xxxxx = contrBuffer.data(g1off + 210 * i + 105);

            auto g1_xzz_xxxxy = contrBuffer.data(g1off + 210 * i + 106);

            auto g1_xzz_xxxxz = contrBuffer.data(g1off + 210 * i + 107);

            auto g1_xzz_xxxyy = contrBuffer.data(g1off + 210 * i + 108);

            auto g1_xzz_xxxyz = contrBuffer.data(g1off + 210 * i + 109);

            auto g1_xzz_xxxzz = contrBuffer.data(g1off + 210 * i + 110);

            auto g1_xzz_xxyyy = contrBuffer.data(g1off + 210 * i + 111);

            auto g1_xzz_xxyyz = contrBuffer.data(g1off + 210 * i + 112);

            auto g1_xzz_xxyzz = contrBuffer.data(g1off + 210 * i + 113);

            auto g1_xzz_xxzzz = contrBuffer.data(g1off + 210 * i + 114);

            auto g1_xzz_xyyyy = contrBuffer.data(g1off + 210 * i + 115);

            auto g1_xzz_xyyyz = contrBuffer.data(g1off + 210 * i + 116);

            auto g1_xzz_xyyzz = contrBuffer.data(g1off + 210 * i + 117);

            auto g1_xzz_xyzzz = contrBuffer.data(g1off + 210 * i + 118);

            auto g1_xzz_xzzzz = contrBuffer.data(g1off + 210 * i + 119);

            auto g1_yyy_xxxxx = contrBuffer.data(g1off + 210 * i + 126);

            auto g1_yyy_xxxxy = contrBuffer.data(g1off + 210 * i + 127);

            auto g1_yyy_xxxxz = contrBuffer.data(g1off + 210 * i + 128);

            auto g1_yyy_xxxyy = contrBuffer.data(g1off + 210 * i + 129);

            auto g1_yyy_xxxyz = contrBuffer.data(g1off + 210 * i + 130);

            auto g1_yyy_xxxzz = contrBuffer.data(g1off + 210 * i + 131);

            auto g1_yyy_xxyyy = contrBuffer.data(g1off + 210 * i + 132);

            auto g1_yyy_xxyyz = contrBuffer.data(g1off + 210 * i + 133);

            auto g1_yyy_xxyzz = contrBuffer.data(g1off + 210 * i + 134);

            auto g1_yyy_xxzzz = contrBuffer.data(g1off + 210 * i + 135);

            auto g1_yyy_xyyyy = contrBuffer.data(g1off + 210 * i + 136);

            auto g1_yyy_xyyyz = contrBuffer.data(g1off + 210 * i + 137);

            auto g1_yyy_xyyzz = contrBuffer.data(g1off + 210 * i + 138);

            auto g1_yyy_xyzzz = contrBuffer.data(g1off + 210 * i + 139);

            auto g1_yyy_xzzzz = contrBuffer.data(g1off + 210 * i + 140);

            auto g1_yyy_yyyyy = contrBuffer.data(g1off + 210 * i + 141);

            auto g1_yyy_yyyyz = contrBuffer.data(g1off + 210 * i + 142);

            auto g1_yyy_yyyzz = contrBuffer.data(g1off + 210 * i + 143);

            auto g1_yyy_yyzzz = contrBuffer.data(g1off + 210 * i + 144);

            auto g1_yyy_yzzzz = contrBuffer.data(g1off + 210 * i + 145);

            auto g1_yyz_xxxxx = contrBuffer.data(g1off + 210 * i + 147);

            auto g1_yyz_xxxxy = contrBuffer.data(g1off + 210 * i + 148);

            auto g1_yyz_xxxxz = contrBuffer.data(g1off + 210 * i + 149);

            auto g1_yyz_xxxyy = contrBuffer.data(g1off + 210 * i + 150);

            auto g1_yyz_xxxyz = contrBuffer.data(g1off + 210 * i + 151);

            auto g1_yyz_xxxzz = contrBuffer.data(g1off + 210 * i + 152);

            auto g1_yyz_xxyyy = contrBuffer.data(g1off + 210 * i + 153);

            auto g1_yyz_xxyyz = contrBuffer.data(g1off + 210 * i + 154);

            auto g1_yyz_xxyzz = contrBuffer.data(g1off + 210 * i + 155);

            auto g1_yyz_xxzzz = contrBuffer.data(g1off + 210 * i + 156);

            auto g1_yyz_xyyyy = contrBuffer.data(g1off + 210 * i + 157);

            auto g1_yyz_xyyyz = contrBuffer.data(g1off + 210 * i + 158);

            auto g1_yyz_xyyzz = contrBuffer.data(g1off + 210 * i + 159);

            auto g1_yyz_xyzzz = contrBuffer.data(g1off + 210 * i + 160);

            auto g1_yyz_xzzzz = contrBuffer.data(g1off + 210 * i + 161);

            auto g1_yyz_yyyyy = contrBuffer.data(g1off + 210 * i + 162);

            auto g1_yyz_yyyyz = contrBuffer.data(g1off + 210 * i + 163);

            auto g1_yyz_yyyzz = contrBuffer.data(g1off + 210 * i + 164);

            auto g1_yyz_yyzzz = contrBuffer.data(g1off + 210 * i + 165);

            auto g1_yyz_yzzzz = contrBuffer.data(g1off + 210 * i + 166);

            auto g1_yzz_xxxxx = contrBuffer.data(g1off + 210 * i + 168);

            auto g1_yzz_xxxxy = contrBuffer.data(g1off + 210 * i + 169);

            auto g1_yzz_xxxxz = contrBuffer.data(g1off + 210 * i + 170);

            auto g1_yzz_xxxyy = contrBuffer.data(g1off + 210 * i + 171);

            auto g1_yzz_xxxyz = contrBuffer.data(g1off + 210 * i + 172);

            auto g1_yzz_xxxzz = contrBuffer.data(g1off + 210 * i + 173);

            auto g1_yzz_xxyyy = contrBuffer.data(g1off + 210 * i + 174);

            auto g1_yzz_xxyyz = contrBuffer.data(g1off + 210 * i + 175);

            auto g1_yzz_xxyzz = contrBuffer.data(g1off + 210 * i + 176);

            auto g1_yzz_xxzzz = contrBuffer.data(g1off + 210 * i + 177);

            auto g1_yzz_xyyyy = contrBuffer.data(g1off + 210 * i + 178);

            auto g1_yzz_xyyyz = contrBuffer.data(g1off + 210 * i + 179);

            auto g1_yzz_xyyzz = contrBuffer.data(g1off + 210 * i + 180);

            auto g1_yzz_xyzzz = contrBuffer.data(g1off + 210 * i + 181);

            auto g1_yzz_xzzzz = contrBuffer.data(g1off + 210 * i + 182);

            auto g1_yzz_yyyyy = contrBuffer.data(g1off + 210 * i + 183);

            auto g1_yzz_yyyyz = contrBuffer.data(g1off + 210 * i + 184);

            auto g1_yzz_yyyzz = contrBuffer.data(g1off + 210 * i + 185);

            auto g1_yzz_yyzzz = contrBuffer.data(g1off + 210 * i + 186);

            auto g1_yzz_yzzzz = contrBuffer.data(g1off + 210 * i + 187);

            auto g1_zzz_xxxxx = contrBuffer.data(g1off + 210 * i + 189);

            auto g1_zzz_xxxxy = contrBuffer.data(g1off + 210 * i + 190);

            auto g1_zzz_xxxxz = contrBuffer.data(g1off + 210 * i + 191);

            auto g1_zzz_xxxyy = contrBuffer.data(g1off + 210 * i + 192);

            auto g1_zzz_xxxyz = contrBuffer.data(g1off + 210 * i + 193);

            auto g1_zzz_xxxzz = contrBuffer.data(g1off + 210 * i + 194);

            auto g1_zzz_xxyyy = contrBuffer.data(g1off + 210 * i + 195);

            auto g1_zzz_xxyyz = contrBuffer.data(g1off + 210 * i + 196);

            auto g1_zzz_xxyzz = contrBuffer.data(g1off + 210 * i + 197);

            auto g1_zzz_xxzzz = contrBuffer.data(g1off + 210 * i + 198);

            auto g1_zzz_xyyyy = contrBuffer.data(g1off + 210 * i + 199);

            auto g1_zzz_xyyyz = contrBuffer.data(g1off + 210 * i + 200);

            auto g1_zzz_xyyzz = contrBuffer.data(g1off + 210 * i + 201);

            auto g1_zzz_xyzzz = contrBuffer.data(g1off + 210 * i + 202);

            auto g1_zzz_xzzzz = contrBuffer.data(g1off + 210 * i + 203);

            auto g1_zzz_yyyyy = contrBuffer.data(g1off + 210 * i + 204);

            auto g1_zzz_yyyyz = contrBuffer.data(g1off + 210 * i + 205);

            auto g1_zzz_yyyzz = contrBuffer.data(g1off + 210 * i + 206);

            auto g1_zzz_yyzzz = contrBuffer.data(g1off + 210 * i + 207);

            auto g1_zzz_yzzzz = contrBuffer.data(g1off + 210 * i + 208);

            auto g1_zzz_zzzzz = contrBuffer.data(g1off + 210 * i + 209);

            // set up pointers to (X|g(r,r')|GG)^(m) integrals

            auto g_xxxx_xxxx = contrBuffer.data(goff + 225 * i);

            auto g_xxxx_xxxy = contrBuffer.data(goff + 225 * i + 1);

            auto g_xxxx_xxxz = contrBuffer.data(goff + 225 * i + 2);

            auto g_xxxx_xxyy = contrBuffer.data(goff + 225 * i + 3);

            auto g_xxxx_xxyz = contrBuffer.data(goff + 225 * i + 4);

            auto g_xxxx_xxzz = contrBuffer.data(goff + 225 * i + 5);

            auto g_xxxx_xyyy = contrBuffer.data(goff + 225 * i + 6);

            auto g_xxxx_xyyz = contrBuffer.data(goff + 225 * i + 7);

            auto g_xxxx_xyzz = contrBuffer.data(goff + 225 * i + 8);

            auto g_xxxx_xzzz = contrBuffer.data(goff + 225 * i + 9);

            auto g_xxxx_yyyy = contrBuffer.data(goff + 225 * i + 10);

            auto g_xxxx_yyyz = contrBuffer.data(goff + 225 * i + 11);

            auto g_xxxx_yyzz = contrBuffer.data(goff + 225 * i + 12);

            auto g_xxxx_yzzz = contrBuffer.data(goff + 225 * i + 13);

            auto g_xxxx_zzzz = contrBuffer.data(goff + 225 * i + 14);

            auto g_xxxy_xxxx = contrBuffer.data(goff + 225 * i + 15);

            auto g_xxxy_xxxy = contrBuffer.data(goff + 225 * i + 16);

            auto g_xxxy_xxxz = contrBuffer.data(goff + 225 * i + 17);

            auto g_xxxy_xxyy = contrBuffer.data(goff + 225 * i + 18);

            auto g_xxxy_xxyz = contrBuffer.data(goff + 225 * i + 19);

            auto g_xxxy_xxzz = contrBuffer.data(goff + 225 * i + 20);

            auto g_xxxy_xyyy = contrBuffer.data(goff + 225 * i + 21);

            auto g_xxxy_xyyz = contrBuffer.data(goff + 225 * i + 22);

            auto g_xxxy_xyzz = contrBuffer.data(goff + 225 * i + 23);

            auto g_xxxy_xzzz = contrBuffer.data(goff + 225 * i + 24);

            auto g_xxxy_yyyy = contrBuffer.data(goff + 225 * i + 25);

            auto g_xxxy_yyyz = contrBuffer.data(goff + 225 * i + 26);

            auto g_xxxy_yyzz = contrBuffer.data(goff + 225 * i + 27);

            auto g_xxxy_yzzz = contrBuffer.data(goff + 225 * i + 28);

            auto g_xxxy_zzzz = contrBuffer.data(goff + 225 * i + 29);

            auto g_xxxz_xxxx = contrBuffer.data(goff + 225 * i + 30);

            auto g_xxxz_xxxy = contrBuffer.data(goff + 225 * i + 31);

            auto g_xxxz_xxxz = contrBuffer.data(goff + 225 * i + 32);

            auto g_xxxz_xxyy = contrBuffer.data(goff + 225 * i + 33);

            auto g_xxxz_xxyz = contrBuffer.data(goff + 225 * i + 34);

            auto g_xxxz_xxzz = contrBuffer.data(goff + 225 * i + 35);

            auto g_xxxz_xyyy = contrBuffer.data(goff + 225 * i + 36);

            auto g_xxxz_xyyz = contrBuffer.data(goff + 225 * i + 37);

            auto g_xxxz_xyzz = contrBuffer.data(goff + 225 * i + 38);

            auto g_xxxz_xzzz = contrBuffer.data(goff + 225 * i + 39);

            auto g_xxxz_yyyy = contrBuffer.data(goff + 225 * i + 40);

            auto g_xxxz_yyyz = contrBuffer.data(goff + 225 * i + 41);

            auto g_xxxz_yyzz = contrBuffer.data(goff + 225 * i + 42);

            auto g_xxxz_yzzz = contrBuffer.data(goff + 225 * i + 43);

            auto g_xxxz_zzzz = contrBuffer.data(goff + 225 * i + 44);

            auto g_xxyy_xxxx = contrBuffer.data(goff + 225 * i + 45);

            auto g_xxyy_xxxy = contrBuffer.data(goff + 225 * i + 46);

            auto g_xxyy_xxxz = contrBuffer.data(goff + 225 * i + 47);

            auto g_xxyy_xxyy = contrBuffer.data(goff + 225 * i + 48);

            auto g_xxyy_xxyz = contrBuffer.data(goff + 225 * i + 49);

            auto g_xxyy_xxzz = contrBuffer.data(goff + 225 * i + 50);

            auto g_xxyy_xyyy = contrBuffer.data(goff + 225 * i + 51);

            auto g_xxyy_xyyz = contrBuffer.data(goff + 225 * i + 52);

            auto g_xxyy_xyzz = contrBuffer.data(goff + 225 * i + 53);

            auto g_xxyy_xzzz = contrBuffer.data(goff + 225 * i + 54);

            auto g_xxyy_yyyy = contrBuffer.data(goff + 225 * i + 55);

            auto g_xxyy_yyyz = contrBuffer.data(goff + 225 * i + 56);

            auto g_xxyy_yyzz = contrBuffer.data(goff + 225 * i + 57);

            auto g_xxyy_yzzz = contrBuffer.data(goff + 225 * i + 58);

            auto g_xxyy_zzzz = contrBuffer.data(goff + 225 * i + 59);

            auto g_xxyz_xxxx = contrBuffer.data(goff + 225 * i + 60);

            auto g_xxyz_xxxy = contrBuffer.data(goff + 225 * i + 61);

            auto g_xxyz_xxxz = contrBuffer.data(goff + 225 * i + 62);

            auto g_xxyz_xxyy = contrBuffer.data(goff + 225 * i + 63);

            auto g_xxyz_xxyz = contrBuffer.data(goff + 225 * i + 64);

            auto g_xxyz_xxzz = contrBuffer.data(goff + 225 * i + 65);

            auto g_xxyz_xyyy = contrBuffer.data(goff + 225 * i + 66);

            auto g_xxyz_xyyz = contrBuffer.data(goff + 225 * i + 67);

            auto g_xxyz_xyzz = contrBuffer.data(goff + 225 * i + 68);

            auto g_xxyz_xzzz = contrBuffer.data(goff + 225 * i + 69);

            auto g_xxyz_yyyy = contrBuffer.data(goff + 225 * i + 70);

            auto g_xxyz_yyyz = contrBuffer.data(goff + 225 * i + 71);

            auto g_xxyz_yyzz = contrBuffer.data(goff + 225 * i + 72);

            auto g_xxyz_yzzz = contrBuffer.data(goff + 225 * i + 73);

            auto g_xxyz_zzzz = contrBuffer.data(goff + 225 * i + 74);

            auto g_xxzz_xxxx = contrBuffer.data(goff + 225 * i + 75);

            auto g_xxzz_xxxy = contrBuffer.data(goff + 225 * i + 76);

            auto g_xxzz_xxxz = contrBuffer.data(goff + 225 * i + 77);

            auto g_xxzz_xxyy = contrBuffer.data(goff + 225 * i + 78);

            auto g_xxzz_xxyz = contrBuffer.data(goff + 225 * i + 79);

            auto g_xxzz_xxzz = contrBuffer.data(goff + 225 * i + 80);

            auto g_xxzz_xyyy = contrBuffer.data(goff + 225 * i + 81);

            auto g_xxzz_xyyz = contrBuffer.data(goff + 225 * i + 82);

            auto g_xxzz_xyzz = contrBuffer.data(goff + 225 * i + 83);

            auto g_xxzz_xzzz = contrBuffer.data(goff + 225 * i + 84);

            auto g_xxzz_yyyy = contrBuffer.data(goff + 225 * i + 85);

            auto g_xxzz_yyyz = contrBuffer.data(goff + 225 * i + 86);

            auto g_xxzz_yyzz = contrBuffer.data(goff + 225 * i + 87);

            auto g_xxzz_yzzz = contrBuffer.data(goff + 225 * i + 88);

            auto g_xxzz_zzzz = contrBuffer.data(goff + 225 * i + 89);

            auto g_xyyy_xxxx = contrBuffer.data(goff + 225 * i + 90);

            auto g_xyyy_xxxy = contrBuffer.data(goff + 225 * i + 91);

            auto g_xyyy_xxxz = contrBuffer.data(goff + 225 * i + 92);

            auto g_xyyy_xxyy = contrBuffer.data(goff + 225 * i + 93);

            auto g_xyyy_xxyz = contrBuffer.data(goff + 225 * i + 94);

            auto g_xyyy_xxzz = contrBuffer.data(goff + 225 * i + 95);

            auto g_xyyy_xyyy = contrBuffer.data(goff + 225 * i + 96);

            auto g_xyyy_xyyz = contrBuffer.data(goff + 225 * i + 97);

            auto g_xyyy_xyzz = contrBuffer.data(goff + 225 * i + 98);

            auto g_xyyy_xzzz = contrBuffer.data(goff + 225 * i + 99);

            auto g_xyyy_yyyy = contrBuffer.data(goff + 225 * i + 100);

            auto g_xyyy_yyyz = contrBuffer.data(goff + 225 * i + 101);

            auto g_xyyy_yyzz = contrBuffer.data(goff + 225 * i + 102);

            auto g_xyyy_yzzz = contrBuffer.data(goff + 225 * i + 103);

            auto g_xyyy_zzzz = contrBuffer.data(goff + 225 * i + 104);

            auto g_xyyz_xxxx = contrBuffer.data(goff + 225 * i + 105);

            auto g_xyyz_xxxy = contrBuffer.data(goff + 225 * i + 106);

            auto g_xyyz_xxxz = contrBuffer.data(goff + 225 * i + 107);

            auto g_xyyz_xxyy = contrBuffer.data(goff + 225 * i + 108);

            auto g_xyyz_xxyz = contrBuffer.data(goff + 225 * i + 109);

            auto g_xyyz_xxzz = contrBuffer.data(goff + 225 * i + 110);

            auto g_xyyz_xyyy = contrBuffer.data(goff + 225 * i + 111);

            auto g_xyyz_xyyz = contrBuffer.data(goff + 225 * i + 112);

            auto g_xyyz_xyzz = contrBuffer.data(goff + 225 * i + 113);

            auto g_xyyz_xzzz = contrBuffer.data(goff + 225 * i + 114);

            auto g_xyyz_yyyy = contrBuffer.data(goff + 225 * i + 115);

            auto g_xyyz_yyyz = contrBuffer.data(goff + 225 * i + 116);

            auto g_xyyz_yyzz = contrBuffer.data(goff + 225 * i + 117);

            auto g_xyyz_yzzz = contrBuffer.data(goff + 225 * i + 118);

            auto g_xyyz_zzzz = contrBuffer.data(goff + 225 * i + 119);

            auto g_xyzz_xxxx = contrBuffer.data(goff + 225 * i + 120);

            auto g_xyzz_xxxy = contrBuffer.data(goff + 225 * i + 121);

            auto g_xyzz_xxxz = contrBuffer.data(goff + 225 * i + 122);

            auto g_xyzz_xxyy = contrBuffer.data(goff + 225 * i + 123);

            auto g_xyzz_xxyz = contrBuffer.data(goff + 225 * i + 124);

            auto g_xyzz_xxzz = contrBuffer.data(goff + 225 * i + 125);

            auto g_xyzz_xyyy = contrBuffer.data(goff + 225 * i + 126);

            auto g_xyzz_xyyz = contrBuffer.data(goff + 225 * i + 127);

            auto g_xyzz_xyzz = contrBuffer.data(goff + 225 * i + 128);

            auto g_xyzz_xzzz = contrBuffer.data(goff + 225 * i + 129);

            auto g_xyzz_yyyy = contrBuffer.data(goff + 225 * i + 130);

            auto g_xyzz_yyyz = contrBuffer.data(goff + 225 * i + 131);

            auto g_xyzz_yyzz = contrBuffer.data(goff + 225 * i + 132);

            auto g_xyzz_yzzz = contrBuffer.data(goff + 225 * i + 133);

            auto g_xyzz_zzzz = contrBuffer.data(goff + 225 * i + 134);

            auto g_xzzz_xxxx = contrBuffer.data(goff + 225 * i + 135);

            auto g_xzzz_xxxy = contrBuffer.data(goff + 225 * i + 136);

            auto g_xzzz_xxxz = contrBuffer.data(goff + 225 * i + 137);

            auto g_xzzz_xxyy = contrBuffer.data(goff + 225 * i + 138);

            auto g_xzzz_xxyz = contrBuffer.data(goff + 225 * i + 139);

            auto g_xzzz_xxzz = contrBuffer.data(goff + 225 * i + 140);

            auto g_xzzz_xyyy = contrBuffer.data(goff + 225 * i + 141);

            auto g_xzzz_xyyz = contrBuffer.data(goff + 225 * i + 142);

            auto g_xzzz_xyzz = contrBuffer.data(goff + 225 * i + 143);

            auto g_xzzz_xzzz = contrBuffer.data(goff + 225 * i + 144);

            auto g_xzzz_yyyy = contrBuffer.data(goff + 225 * i + 145);

            auto g_xzzz_yyyz = contrBuffer.data(goff + 225 * i + 146);

            auto g_xzzz_yyzz = contrBuffer.data(goff + 225 * i + 147);

            auto g_xzzz_yzzz = contrBuffer.data(goff + 225 * i + 148);

            auto g_xzzz_zzzz = contrBuffer.data(goff + 225 * i + 149);

            auto g_yyyy_xxxx = contrBuffer.data(goff + 225 * i + 150);

            auto g_yyyy_xxxy = contrBuffer.data(goff + 225 * i + 151);

            auto g_yyyy_xxxz = contrBuffer.data(goff + 225 * i + 152);

            auto g_yyyy_xxyy = contrBuffer.data(goff + 225 * i + 153);

            auto g_yyyy_xxyz = contrBuffer.data(goff + 225 * i + 154);

            auto g_yyyy_xxzz = contrBuffer.data(goff + 225 * i + 155);

            auto g_yyyy_xyyy = contrBuffer.data(goff + 225 * i + 156);

            auto g_yyyy_xyyz = contrBuffer.data(goff + 225 * i + 157);

            auto g_yyyy_xyzz = contrBuffer.data(goff + 225 * i + 158);

            auto g_yyyy_xzzz = contrBuffer.data(goff + 225 * i + 159);

            auto g_yyyy_yyyy = contrBuffer.data(goff + 225 * i + 160);

            auto g_yyyy_yyyz = contrBuffer.data(goff + 225 * i + 161);

            auto g_yyyy_yyzz = contrBuffer.data(goff + 225 * i + 162);

            auto g_yyyy_yzzz = contrBuffer.data(goff + 225 * i + 163);

            auto g_yyyy_zzzz = contrBuffer.data(goff + 225 * i + 164);

            auto g_yyyz_xxxx = contrBuffer.data(goff + 225 * i + 165);

            auto g_yyyz_xxxy = contrBuffer.data(goff + 225 * i + 166);

            auto g_yyyz_xxxz = contrBuffer.data(goff + 225 * i + 167);

            auto g_yyyz_xxyy = contrBuffer.data(goff + 225 * i + 168);

            auto g_yyyz_xxyz = contrBuffer.data(goff + 225 * i + 169);

            auto g_yyyz_xxzz = contrBuffer.data(goff + 225 * i + 170);

            auto g_yyyz_xyyy = contrBuffer.data(goff + 225 * i + 171);

            auto g_yyyz_xyyz = contrBuffer.data(goff + 225 * i + 172);

            auto g_yyyz_xyzz = contrBuffer.data(goff + 225 * i + 173);

            auto g_yyyz_xzzz = contrBuffer.data(goff + 225 * i + 174);

            auto g_yyyz_yyyy = contrBuffer.data(goff + 225 * i + 175);

            auto g_yyyz_yyyz = contrBuffer.data(goff + 225 * i + 176);

            auto g_yyyz_yyzz = contrBuffer.data(goff + 225 * i + 177);

            auto g_yyyz_yzzz = contrBuffer.data(goff + 225 * i + 178);

            auto g_yyyz_zzzz = contrBuffer.data(goff + 225 * i + 179);

            auto g_yyzz_xxxx = contrBuffer.data(goff + 225 * i + 180);

            auto g_yyzz_xxxy = contrBuffer.data(goff + 225 * i + 181);

            auto g_yyzz_xxxz = contrBuffer.data(goff + 225 * i + 182);

            auto g_yyzz_xxyy = contrBuffer.data(goff + 225 * i + 183);

            auto g_yyzz_xxyz = contrBuffer.data(goff + 225 * i + 184);

            auto g_yyzz_xxzz = contrBuffer.data(goff + 225 * i + 185);

            auto g_yyzz_xyyy = contrBuffer.data(goff + 225 * i + 186);

            auto g_yyzz_xyyz = contrBuffer.data(goff + 225 * i + 187);

            auto g_yyzz_xyzz = contrBuffer.data(goff + 225 * i + 188);

            auto g_yyzz_xzzz = contrBuffer.data(goff + 225 * i + 189);

            auto g_yyzz_yyyy = contrBuffer.data(goff + 225 * i + 190);

            auto g_yyzz_yyyz = contrBuffer.data(goff + 225 * i + 191);

            auto g_yyzz_yyzz = contrBuffer.data(goff + 225 * i + 192);

            auto g_yyzz_yzzz = contrBuffer.data(goff + 225 * i + 193);

            auto g_yyzz_zzzz = contrBuffer.data(goff + 225 * i + 194);

            auto g_yzzz_xxxx = contrBuffer.data(goff + 225 * i + 195);

            auto g_yzzz_xxxy = contrBuffer.data(goff + 225 * i + 196);

            auto g_yzzz_xxxz = contrBuffer.data(goff + 225 * i + 197);

            auto g_yzzz_xxyy = contrBuffer.data(goff + 225 * i + 198);

            auto g_yzzz_xxyz = contrBuffer.data(goff + 225 * i + 199);

            auto g_yzzz_xxzz = contrBuffer.data(goff + 225 * i + 200);

            auto g_yzzz_xyyy = contrBuffer.data(goff + 225 * i + 201);

            auto g_yzzz_xyyz = contrBuffer.data(goff + 225 * i + 202);

            auto g_yzzz_xyzz = contrBuffer.data(goff + 225 * i + 203);

            auto g_yzzz_xzzz = contrBuffer.data(goff + 225 * i + 204);

            auto g_yzzz_yyyy = contrBuffer.data(goff + 225 * i + 205);

            auto g_yzzz_yyyz = contrBuffer.data(goff + 225 * i + 206);

            auto g_yzzz_yyzz = contrBuffer.data(goff + 225 * i + 207);

            auto g_yzzz_yzzz = contrBuffer.data(goff + 225 * i + 208);

            auto g_yzzz_zzzz = contrBuffer.data(goff + 225 * i + 209);

            auto g_zzzz_xxxx = contrBuffer.data(goff + 225 * i + 210);

            auto g_zzzz_xxxy = contrBuffer.data(goff + 225 * i + 211);

            auto g_zzzz_xxxz = contrBuffer.data(goff + 225 * i + 212);

            auto g_zzzz_xxyy = contrBuffer.data(goff + 225 * i + 213);

            auto g_zzzz_xxyz = contrBuffer.data(goff + 225 * i + 214);

            auto g_zzzz_xxzz = contrBuffer.data(goff + 225 * i + 215);

            auto g_zzzz_xyyy = contrBuffer.data(goff + 225 * i + 216);

            auto g_zzzz_xyyz = contrBuffer.data(goff + 225 * i + 217);

            auto g_zzzz_xyzz = contrBuffer.data(goff + 225 * i + 218);

            auto g_zzzz_xzzz = contrBuffer.data(goff + 225 * i + 219);

            auto g_zzzz_yyyy = contrBuffer.data(goff + 225 * i + 220);

            auto g_zzzz_yyyz = contrBuffer.data(goff + 225 * i + 221);

            auto g_zzzz_yyzz = contrBuffer.data(goff + 225 * i + 222);

            auto g_zzzz_yzzz = contrBuffer.data(goff + 225 * i + 223);

            auto g_zzzz_zzzz = contrBuffer.data(goff + 225 * i + 224);

            #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xxx_xxxx, g2_xxx_xxxy,\
                                     g2_xxx_xxxz, g2_xxx_xxyy, g2_xxx_xxyz, g2_xxx_xxzz,\
                                     g2_xxx_xyyy, g2_xxx_xyyz, g2_xxx_xyzz, g2_xxx_xzzz,\
                                     g2_xxx_yyyy, g2_xxx_yyyz, g2_xxx_yyzz, g2_xxx_yzzz,\
                                     g2_xxx_zzzz, g2_xxy_xxxx, g2_xxy_xxxy, g2_xxy_xxxz,\
                                     g2_xxy_xxyy, g2_xxy_xxyz, g2_xxy_xxzz, g2_xxy_xyyy,\
                                     g2_xxy_xyyz, g2_xxy_xyzz, g2_xxy_xzzz, g2_xxy_yyyy,\
                                     g2_xxy_yyyz, g2_xxy_yyzz, g2_xxy_yzzz, g2_xxy_zzzz,\
                                     g2_xxz_xxxx, g2_xxz_xxxy, g2_xxz_xxxz, g2_xxz_xxyy,\
                                     g2_xxz_xxyz, g2_xxz_xxzz, g2_xxz_xyyy, g2_xxz_xyyz,\
                                     g2_xxz_xyzz, g2_xxz_xzzz, g2_xxz_yyyy, g2_xxz_yyyz,\
                                     g2_xxz_yyzz, g2_xxz_yzzz, g2_xxz_zzzz, g2_xyy_xxxx,\
                                     g2_xyy_xxxy, g2_xyy_xxxz, g2_xyy_xxyy, g2_xyy_xxyz,\
                                     g2_xyy_xxzz, g2_xyy_xyyy, g2_xyy_xyyz, g2_xyy_xyzz,\
                                     g2_xyy_xzzz, g2_xyy_yyyy, g2_xyy_yyyz, g2_xyy_yyzz,\
                                     g2_xyy_yzzz, g2_xyy_zzzz, g2_xyz_xxxx, g2_xyz_xxxy,\
                                     g2_xyz_xxxz, g2_xyz_xxyy, g2_xyz_xxyz, g2_xyz_xxzz,\
                                     g2_xyz_xyyy, g2_xyz_xyyz, g2_xyz_xyzz, g2_xyz_xzzz,\
                                     g2_xyz_yyyy, g2_xyz_yyyz, g2_xyz_yyzz, g2_xyz_yzzz,\
                                     g2_xyz_zzzz, g2_xzz_xxxx, g2_xzz_xxxy, g2_xzz_xxxz,\
                                     g2_xzz_xxyy, g2_xzz_xxyz, g2_xzz_xxzz, g2_xzz_xyyy,\
                                     g2_xzz_xyyz, g2_xzz_xyzz, g2_xzz_xzzz, g2_xzz_yyyy,\
                                     g2_xzz_yyyz, g2_xzz_yyzz, g2_xzz_yzzz, g2_xzz_zzzz,\
                                     g2_yyy_xxxx, g2_yyy_xxxy, g2_yyy_xxxz, g2_yyy_xxyy,\
                                     g2_yyy_xxyz, g2_yyy_xxzz, g2_yyy_xyyy, g2_yyy_xyyz,\
                                     g2_yyy_xyzz, g2_yyy_xzzz, g2_yyy_yyyy, g2_yyy_yyyz,\
                                     g2_yyy_yyzz, g2_yyy_yzzz, g2_yyy_zzzz, g2_yyz_xxxx,\
                                     g2_yyz_xxxy, g2_yyz_xxxz, g2_yyz_xxyy, g2_yyz_xxyz,\
                                     g2_yyz_xxzz, g2_yyz_xyyy, g2_yyz_xyyz, g2_yyz_xyzz,\
                                     g2_yyz_xzzz, g2_yyz_yyyy, g2_yyz_yyyz, g2_yyz_yyzz,\
                                     g2_yyz_yzzz, g2_yyz_zzzz, g2_yzz_xxxx, g2_yzz_xxxy,\
                                     g2_yzz_xxxz, g2_yzz_xxyy, g2_yzz_xxyz, g2_yzz_xxzz,\
                                     g2_yzz_xyyy, g2_yzz_xyyz, g2_yzz_xyzz, g2_yzz_xzzz,\
                                     g2_yzz_yyyy, g2_yzz_yyyz, g2_yzz_yyzz, g2_yzz_yzzz,\
                                     g2_yzz_zzzz, g2_zzz_xxxx, g2_zzz_xxxy, g2_zzz_xxxz,\
                                     g2_zzz_xxyy, g2_zzz_xxyz, g2_zzz_xxzz, g2_zzz_xyyy,\
                                     g2_zzz_xyyz, g2_zzz_xyzz, g2_zzz_xzzz, g2_zzz_yyyy,\
                                     g2_zzz_yyyz, g2_zzz_yyzz, g2_zzz_yzzz, g2_zzz_zzzz,\
                                     g1_xxx_xxxxx, g1_xxx_xxxxy, g1_xxx_xxxxz,\
                                     g1_xxx_xxxyy, g1_xxx_xxxyz, g1_xxx_xxxzz,\
                                     g1_xxx_xxyyy, g1_xxx_xxyyz, g1_xxx_xxyzz,\
                                     g1_xxx_xxzzz, g1_xxx_xyyyy, g1_xxx_xyyyz,\
                                     g1_xxx_xyyzz, g1_xxx_xyzzz, g1_xxx_xzzzz,\
                                     g1_xxy_xxxxx, g1_xxy_xxxxy, g1_xxy_xxxxz,\
                                     g1_xxy_xxxyy, g1_xxy_xxxyz, g1_xxy_xxxzz,\
                                     g1_xxy_xxyyy, g1_xxy_xxyyz, g1_xxy_xxyzz,\
                                     g1_xxy_xxzzz, g1_xxy_xyyyy, g1_xxy_xyyyz,\
                                     g1_xxy_xyyzz, g1_xxy_xyzzz, g1_xxy_xzzzz,\
                                     g1_xxz_xxxxx, g1_xxz_xxxxy, g1_xxz_xxxxz,\
                                     g1_xxz_xxxyy, g1_xxz_xxxyz, g1_xxz_xxxzz,\
                                     g1_xxz_xxyyy, g1_xxz_xxyyz, g1_xxz_xxyzz,\
                                     g1_xxz_xxzzz, g1_xxz_xyyyy, g1_xxz_xyyyz,\
                                     g1_xxz_xyyzz, g1_xxz_xyzzz, g1_xxz_xzzzz,\
                                     g1_xyy_xxxxx, g1_xyy_xxxxy, g1_xyy_xxxxz,\
                                     g1_xyy_xxxyy, g1_xyy_xxxyz, g1_xyy_xxxzz,\
                                     g1_xyy_xxyyy, g1_xyy_xxyyz, g1_xyy_xxyzz,\
                                     g1_xyy_xxzzz, g1_xyy_xyyyy, g1_xyy_xyyyz,\
                                     g1_xyy_xyyzz, g1_xyy_xyzzz, g1_xyy_xzzzz,\
                                     g1_xyz_xxxxx, g1_xyz_xxxxy, g1_xyz_xxxxz,\
                                     g1_xyz_xxxyy, g1_xyz_xxxyz, g1_xyz_xxxzz,\
                                     g1_xyz_xxyyy, g1_xyz_xxyyz, g1_xyz_xxyzz,\
                                     g1_xyz_xxzzz, g1_xyz_xyyyy, g1_xyz_xyyyz,\
                                     g1_xyz_xyyzz, g1_xyz_xyzzz, g1_xyz_xzzzz,\
                                     g1_xzz_xxxxx, g1_xzz_xxxxy, g1_xzz_xxxxz,\
                                     g1_xzz_xxxyy, g1_xzz_xxxyz, g1_xzz_xxxzz,\
                                     g1_xzz_xxyyy, g1_xzz_xxyyz, g1_xzz_xxyzz,\
                                     g1_xzz_xxzzz, g1_xzz_xyyyy, g1_xzz_xyyyz,\
                                     g1_xzz_xyyzz, g1_xzz_xyzzz, g1_xzz_xzzzz,\
                                     g1_yyy_xxxxx, g1_yyy_xxxxy, g1_yyy_xxxxz,\
                                     g1_yyy_xxxyy, g1_yyy_xxxyz, g1_yyy_xxxzz,\
                                     g1_yyy_xxyyy, g1_yyy_xxyyz, g1_yyy_xxyzz,\
                                     g1_yyy_xxzzz, g1_yyy_xyyyy, g1_yyy_xyyyz,\
                                     g1_yyy_xyyzz, g1_yyy_xyzzz, g1_yyy_xzzzz,\
                                     g1_yyy_yyyyy, g1_yyy_yyyyz, g1_yyy_yyyzz,\
                                     g1_yyy_yyzzz, g1_yyy_yzzzz,\
                                     g1_yyz_xxxxx, g1_yyz_xxxxy, g1_yyz_xxxxz,\
                                     g1_yyz_xxxyy, g1_yyz_xxxyz, g1_yyz_xxxzz,\
                                     g1_yyz_xxyyy, g1_yyz_xxyyz, g1_yyz_xxyzz,\
                                     g1_yyz_xxzzz, g1_yyz_xyyyy, g1_yyz_xyyyz,\
                                     g1_yyz_xyyzz, g1_yyz_xyzzz, g1_yyz_xzzzz,\
                                     g1_yyz_yyyyy, g1_yyz_yyyyz, g1_yyz_yyyzz,\
                                     g1_yyz_yyzzz, g1_yyz_yzzzz,\
                                     g1_yzz_xxxxx, g1_yzz_xxxxy, g1_yzz_xxxxz,\
                                     g1_yzz_xxxyy, g1_yzz_xxxyz, g1_yzz_xxxzz,\
                                     g1_yzz_xxyyy, g1_yzz_xxyyz, g1_yzz_xxyzz,\
                                     g1_yzz_xxzzz, g1_yzz_xyyyy, g1_yzz_xyyyz,\
                                     g1_yzz_xyyzz, g1_yzz_xyzzz, g1_yzz_xzzzz,\
                                     g1_yzz_yyyyy, g1_yzz_yyyyz, g1_yzz_yyyzz,\
                                     g1_yzz_yyzzz, g1_yzz_yzzzz,\
                                     g1_zzz_xxxxx, g1_zzz_xxxxy, g1_zzz_xxxxz,\
                                     g1_zzz_xxxyy, g1_zzz_xxxyz, g1_zzz_xxxzz,\
                                     g1_zzz_xxyyy, g1_zzz_xxyyz, g1_zzz_xxyzz,\
                                     g1_zzz_xxzzz, g1_zzz_xyyyy, g1_zzz_xyyyz,\
                                     g1_zzz_xyyzz, g1_zzz_xyzzz, g1_zzz_xzzzz,\
                                     g1_zzz_yyyyy, g1_zzz_yyyyz, g1_zzz_yyyzz,\
                                     g1_zzz_yyzzz, g1_zzz_yzzzz, g1_zzz_zzzzz,\
                                     g_xxxx_xxxx, g_xxxx_xxxy, g_xxxx_xxxz, g_xxxx_xxyy,\
                                     g_xxxx_xxyz, g_xxxx_xxzz, g_xxxx_xyyy, g_xxxx_xyyz,\
                                     g_xxxx_xyzz, g_xxxx_xzzz, g_xxxx_yyyy, g_xxxx_yyyz,\
                                     g_xxxx_yyzz, g_xxxx_yzzz, g_xxxx_zzzz, g_xxxy_xxxx,\
                                     g_xxxy_xxxy, g_xxxy_xxxz, g_xxxy_xxyy, g_xxxy_xxyz,\
                                     g_xxxy_xxzz, g_xxxy_xyyy, g_xxxy_xyyz, g_xxxy_xyzz,\
                                     g_xxxy_xzzz, g_xxxy_yyyy, g_xxxy_yyyz, g_xxxy_yyzz,\
                                     g_xxxy_yzzz, g_xxxy_zzzz, g_xxxz_xxxx, g_xxxz_xxxy,\
                                     g_xxxz_xxxz, g_xxxz_xxyy, g_xxxz_xxyz, g_xxxz_xxzz,\
                                     g_xxxz_xyyy, g_xxxz_xyyz, g_xxxz_xyzz, g_xxxz_xzzz,\
                                     g_xxxz_yyyy, g_xxxz_yyyz, g_xxxz_yyzz, g_xxxz_yzzz,\
                                     g_xxxz_zzzz, g_xxyy_xxxx, g_xxyy_xxxy, g_xxyy_xxxz,\
                                     g_xxyy_xxyy, g_xxyy_xxyz, g_xxyy_xxzz, g_xxyy_xyyy,\
                                     g_xxyy_xyyz, g_xxyy_xyzz, g_xxyy_xzzz, g_xxyy_yyyy,\
                                     g_xxyy_yyyz, g_xxyy_yyzz, g_xxyy_yzzz, g_xxyy_zzzz,\
                                     g_xxyz_xxxx, g_xxyz_xxxy, g_xxyz_xxxz, g_xxyz_xxyy,\
                                     g_xxyz_xxyz, g_xxyz_xxzz, g_xxyz_xyyy, g_xxyz_xyyz,\
                                     g_xxyz_xyzz, g_xxyz_xzzz, g_xxyz_yyyy, g_xxyz_yyyz,\
                                     g_xxyz_yyzz, g_xxyz_yzzz, g_xxyz_zzzz, g_xxzz_xxxx,\
                                     g_xxzz_xxxy, g_xxzz_xxxz, g_xxzz_xxyy, g_xxzz_xxyz,\
                                     g_xxzz_xxzz, g_xxzz_xyyy, g_xxzz_xyyz, g_xxzz_xyzz,\
                                     g_xxzz_xzzz, g_xxzz_yyyy, g_xxzz_yyyz, g_xxzz_yyzz,\
                                     g_xxzz_yzzz, g_xxzz_zzzz, g_xyyy_xxxx, g_xyyy_xxxy,\
                                     g_xyyy_xxxz, g_xyyy_xxyy, g_xyyy_xxyz, g_xyyy_xxzz,\
                                     g_xyyy_xyyy, g_xyyy_xyyz, g_xyyy_xyzz, g_xyyy_xzzz,\
                                     g_xyyy_yyyy, g_xyyy_yyyz, g_xyyy_yyzz, g_xyyy_yzzz,\
                                     g_xyyy_zzzz, g_xyyz_xxxx, g_xyyz_xxxy, g_xyyz_xxxz,\
                                     g_xyyz_xxyy, g_xyyz_xxyz, g_xyyz_xxzz, g_xyyz_xyyy,\
                                     g_xyyz_xyyz, g_xyyz_xyzz, g_xyyz_xzzz, g_xyyz_yyyy,\
                                     g_xyyz_yyyz, g_xyyz_yyzz, g_xyyz_yzzz, g_xyyz_zzzz,\
                                     g_xyzz_xxxx, g_xyzz_xxxy, g_xyzz_xxxz, g_xyzz_xxyy,\
                                     g_xyzz_xxyz, g_xyzz_xxzz, g_xyzz_xyyy, g_xyzz_xyyz,\
                                     g_xyzz_xyzz, g_xyzz_xzzz, g_xyzz_yyyy, g_xyzz_yyyz,\
                                     g_xyzz_yyzz, g_xyzz_yzzz, g_xyzz_zzzz, g_xzzz_xxxx,\
                                     g_xzzz_xxxy, g_xzzz_xxxz, g_xzzz_xxyy, g_xzzz_xxyz,\
                                     g_xzzz_xxzz, g_xzzz_xyyy, g_xzzz_xyyz, g_xzzz_xyzz,\
                                     g_xzzz_xzzz, g_xzzz_yyyy, g_xzzz_yyyz, g_xzzz_yyzz,\
                                     g_xzzz_yzzz, g_xzzz_zzzz, g_yyyy_xxxx, g_yyyy_xxxy,\
                                     g_yyyy_xxxz, g_yyyy_xxyy, g_yyyy_xxyz, g_yyyy_xxzz,\
                                     g_yyyy_xyyy, g_yyyy_xyyz, g_yyyy_xyzz, g_yyyy_xzzz,\
                                     g_yyyy_yyyy, g_yyyy_yyyz, g_yyyy_yyzz, g_yyyy_yzzz,\
                                     g_yyyy_zzzz, g_yyyz_xxxx, g_yyyz_xxxy, g_yyyz_xxxz,\
                                     g_yyyz_xxyy, g_yyyz_xxyz, g_yyyz_xxzz, g_yyyz_xyyy,\
                                     g_yyyz_xyyz, g_yyyz_xyzz, g_yyyz_xzzz, g_yyyz_yyyy,\
                                     g_yyyz_yyyz, g_yyyz_yyzz, g_yyyz_yzzz, g_yyyz_zzzz,\
                                     g_yyzz_xxxx, g_yyzz_xxxy, g_yyzz_xxxz, g_yyzz_xxyy,\
                                     g_yyzz_xxyz, g_yyzz_xxzz, g_yyzz_xyyy, g_yyzz_xyyz,\
                                     g_yyzz_xyzz, g_yyzz_xzzz, g_yyzz_yyyy, g_yyzz_yyyz,\
                                     g_yyzz_yyzz, g_yyzz_yzzz, g_yyzz_zzzz, g_yzzz_xxxx,\
                                     g_yzzz_xxxy, g_yzzz_xxxz, g_yzzz_xxyy, g_yzzz_xxyz,\
                                     g_yzzz_xxzz, g_yzzz_xyyy, g_yzzz_xyyz, g_yzzz_xyzz,\
                                     g_yzzz_xzzz, g_yzzz_yyyy, g_yzzz_yyyz, g_yzzz_yyzz,\
                                     g_yzzz_yzzz, g_yzzz_zzzz, g_zzzz_xxxx, g_zzzz_xxxy,\
                                     g_zzzz_xxxz, g_zzzz_xxyy, g_zzzz_xxyz, g_zzzz_xxzz,\
                                     g_zzzz_xyyy, g_zzzz_xyyz, g_zzzz_xyzz, g_zzzz_xzzz,\
                                     g_zzzz_yyyy, g_zzzz_yyyz, g_zzzz_yyzz, g_zzzz_yzzz,\
                                     g_zzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < kdim; j++)
            {
                // leading x component

                double fr = rcdx[j];

                g_xxxx_xxxx[j] = g1_xxx_xxxxx[j] - fr * g2_xxx_xxxx[j];

                g_xxxx_xxxy[j] = g1_xxx_xxxxy[j] - fr * g2_xxx_xxxy[j];

                g_xxxx_xxxz[j] = g1_xxx_xxxxz[j] - fr * g2_xxx_xxxz[j];

                g_xxxx_xxyy[j] = g1_xxx_xxxyy[j] - fr * g2_xxx_xxyy[j];

                g_xxxx_xxyz[j] = g1_xxx_xxxyz[j] - fr * g2_xxx_xxyz[j];

                g_xxxx_xxzz[j] = g1_xxx_xxxzz[j] - fr * g2_xxx_xxzz[j];

                g_xxxx_xyyy[j] = g1_xxx_xxyyy[j] - fr * g2_xxx_xyyy[j];

                g_xxxx_xyyz[j] = g1_xxx_xxyyz[j] - fr * g2_xxx_xyyz[j];

                g_xxxx_xyzz[j] = g1_xxx_xxyzz[j] - fr * g2_xxx_xyzz[j];

                g_xxxx_xzzz[j] = g1_xxx_xxzzz[j] - fr * g2_xxx_xzzz[j];

                g_xxxx_yyyy[j] = g1_xxx_xyyyy[j] - fr * g2_xxx_yyyy[j];

                g_xxxx_yyyz[j] = g1_xxx_xyyyz[j] - fr * g2_xxx_yyyz[j];

                g_xxxx_yyzz[j] = g1_xxx_xyyzz[j] - fr * g2_xxx_yyzz[j];

                g_xxxx_yzzz[j] = g1_xxx_xyzzz[j] - fr * g2_xxx_yzzz[j];

                g_xxxx_zzzz[j] = g1_xxx_xzzzz[j] - fr * g2_xxx_zzzz[j];

                g_xxxy_xxxx[j] = g1_xxy_xxxxx[j] - fr * g2_xxy_xxxx[j];

                g_xxxy_xxxy[j] = g1_xxy_xxxxy[j] - fr * g2_xxy_xxxy[j];

                g_xxxy_xxxz[j] = g1_xxy_xxxxz[j] - fr * g2_xxy_xxxz[j];

                g_xxxy_xxyy[j] = g1_xxy_xxxyy[j] - fr * g2_xxy_xxyy[j];

                g_xxxy_xxyz[j] = g1_xxy_xxxyz[j] - fr * g2_xxy_xxyz[j];

                g_xxxy_xxzz[j] = g1_xxy_xxxzz[j] - fr * g2_xxy_xxzz[j];

                g_xxxy_xyyy[j] = g1_xxy_xxyyy[j] - fr * g2_xxy_xyyy[j];

                g_xxxy_xyyz[j] = g1_xxy_xxyyz[j] - fr * g2_xxy_xyyz[j];

                g_xxxy_xyzz[j] = g1_xxy_xxyzz[j] - fr * g2_xxy_xyzz[j];

                g_xxxy_xzzz[j] = g1_xxy_xxzzz[j] - fr * g2_xxy_xzzz[j];

                g_xxxy_yyyy[j] = g1_xxy_xyyyy[j] - fr * g2_xxy_yyyy[j];

                g_xxxy_yyyz[j] = g1_xxy_xyyyz[j] - fr * g2_xxy_yyyz[j];

                g_xxxy_yyzz[j] = g1_xxy_xyyzz[j] - fr * g2_xxy_yyzz[j];

                g_xxxy_yzzz[j] = g1_xxy_xyzzz[j] - fr * g2_xxy_yzzz[j];

                g_xxxy_zzzz[j] = g1_xxy_xzzzz[j] - fr * g2_xxy_zzzz[j];

                g_xxxz_xxxx[j] = g1_xxz_xxxxx[j] - fr * g2_xxz_xxxx[j];

                g_xxxz_xxxy[j] = g1_xxz_xxxxy[j] - fr * g2_xxz_xxxy[j];

                g_xxxz_xxxz[j] = g1_xxz_xxxxz[j] - fr * g2_xxz_xxxz[j];

                g_xxxz_xxyy[j] = g1_xxz_xxxyy[j] - fr * g2_xxz_xxyy[j];

                g_xxxz_xxyz[j] = g1_xxz_xxxyz[j] - fr * g2_xxz_xxyz[j];

                g_xxxz_xxzz[j] = g1_xxz_xxxzz[j] - fr * g2_xxz_xxzz[j];

                g_xxxz_xyyy[j] = g1_xxz_xxyyy[j] - fr * g2_xxz_xyyy[j];

                g_xxxz_xyyz[j] = g1_xxz_xxyyz[j] - fr * g2_xxz_xyyz[j];

                g_xxxz_xyzz[j] = g1_xxz_xxyzz[j] - fr * g2_xxz_xyzz[j];

                g_xxxz_xzzz[j] = g1_xxz_xxzzz[j] - fr * g2_xxz_xzzz[j];

                g_xxxz_yyyy[j] = g1_xxz_xyyyy[j] - fr * g2_xxz_yyyy[j];

                g_xxxz_yyyz[j] = g1_xxz_xyyyz[j] - fr * g2_xxz_yyyz[j];

                g_xxxz_yyzz[j] = g1_xxz_xyyzz[j] - fr * g2_xxz_yyzz[j];

                g_xxxz_yzzz[j] = g1_xxz_xyzzz[j] - fr * g2_xxz_yzzz[j];

                g_xxxz_zzzz[j] = g1_xxz_xzzzz[j] - fr * g2_xxz_zzzz[j];

                g_xxyy_xxxx[j] = g1_xyy_xxxxx[j] - fr * g2_xyy_xxxx[j];

                g_xxyy_xxxy[j] = g1_xyy_xxxxy[j] - fr * g2_xyy_xxxy[j];

                g_xxyy_xxxz[j] = g1_xyy_xxxxz[j] - fr * g2_xyy_xxxz[j];

                g_xxyy_xxyy[j] = g1_xyy_xxxyy[j] - fr * g2_xyy_xxyy[j];

                g_xxyy_xxyz[j] = g1_xyy_xxxyz[j] - fr * g2_xyy_xxyz[j];

                g_xxyy_xxzz[j] = g1_xyy_xxxzz[j] - fr * g2_xyy_xxzz[j];

                g_xxyy_xyyy[j] = g1_xyy_xxyyy[j] - fr * g2_xyy_xyyy[j];

                g_xxyy_xyyz[j] = g1_xyy_xxyyz[j] - fr * g2_xyy_xyyz[j];

                g_xxyy_xyzz[j] = g1_xyy_xxyzz[j] - fr * g2_xyy_xyzz[j];

                g_xxyy_xzzz[j] = g1_xyy_xxzzz[j] - fr * g2_xyy_xzzz[j];

                g_xxyy_yyyy[j] = g1_xyy_xyyyy[j] - fr * g2_xyy_yyyy[j];

                g_xxyy_yyyz[j] = g1_xyy_xyyyz[j] - fr * g2_xyy_yyyz[j];

                g_xxyy_yyzz[j] = g1_xyy_xyyzz[j] - fr * g2_xyy_yyzz[j];

                g_xxyy_yzzz[j] = g1_xyy_xyzzz[j] - fr * g2_xyy_yzzz[j];

                g_xxyy_zzzz[j] = g1_xyy_xzzzz[j] - fr * g2_xyy_zzzz[j];

                g_xxyz_xxxx[j] = g1_xyz_xxxxx[j] - fr * g2_xyz_xxxx[j];

                g_xxyz_xxxy[j] = g1_xyz_xxxxy[j] - fr * g2_xyz_xxxy[j];

                g_xxyz_xxxz[j] = g1_xyz_xxxxz[j] - fr * g2_xyz_xxxz[j];

                g_xxyz_xxyy[j] = g1_xyz_xxxyy[j] - fr * g2_xyz_xxyy[j];

                g_xxyz_xxyz[j] = g1_xyz_xxxyz[j] - fr * g2_xyz_xxyz[j];

                g_xxyz_xxzz[j] = g1_xyz_xxxzz[j] - fr * g2_xyz_xxzz[j];

                g_xxyz_xyyy[j] = g1_xyz_xxyyy[j] - fr * g2_xyz_xyyy[j];

                g_xxyz_xyyz[j] = g1_xyz_xxyyz[j] - fr * g2_xyz_xyyz[j];

                g_xxyz_xyzz[j] = g1_xyz_xxyzz[j] - fr * g2_xyz_xyzz[j];

                g_xxyz_xzzz[j] = g1_xyz_xxzzz[j] - fr * g2_xyz_xzzz[j];

                g_xxyz_yyyy[j] = g1_xyz_xyyyy[j] - fr * g2_xyz_yyyy[j];

                g_xxyz_yyyz[j] = g1_xyz_xyyyz[j] - fr * g2_xyz_yyyz[j];

                g_xxyz_yyzz[j] = g1_xyz_xyyzz[j] - fr * g2_xyz_yyzz[j];

                g_xxyz_yzzz[j] = g1_xyz_xyzzz[j] - fr * g2_xyz_yzzz[j];

                g_xxyz_zzzz[j] = g1_xyz_xzzzz[j] - fr * g2_xyz_zzzz[j];

                g_xxzz_xxxx[j] = g1_xzz_xxxxx[j] - fr * g2_xzz_xxxx[j];

                g_xxzz_xxxy[j] = g1_xzz_xxxxy[j] - fr * g2_xzz_xxxy[j];

                g_xxzz_xxxz[j] = g1_xzz_xxxxz[j] - fr * g2_xzz_xxxz[j];

                g_xxzz_xxyy[j] = g1_xzz_xxxyy[j] - fr * g2_xzz_xxyy[j];

                g_xxzz_xxyz[j] = g1_xzz_xxxyz[j] - fr * g2_xzz_xxyz[j];

                g_xxzz_xxzz[j] = g1_xzz_xxxzz[j] - fr * g2_xzz_xxzz[j];

                g_xxzz_xyyy[j] = g1_xzz_xxyyy[j] - fr * g2_xzz_xyyy[j];

                g_xxzz_xyyz[j] = g1_xzz_xxyyz[j] - fr * g2_xzz_xyyz[j];

                g_xxzz_xyzz[j] = g1_xzz_xxyzz[j] - fr * g2_xzz_xyzz[j];

                g_xxzz_xzzz[j] = g1_xzz_xxzzz[j] - fr * g2_xzz_xzzz[j];

                g_xxzz_yyyy[j] = g1_xzz_xyyyy[j] - fr * g2_xzz_yyyy[j];

                g_xxzz_yyyz[j] = g1_xzz_xyyyz[j] - fr * g2_xzz_yyyz[j];

                g_xxzz_yyzz[j] = g1_xzz_xyyzz[j] - fr * g2_xzz_yyzz[j];

                g_xxzz_yzzz[j] = g1_xzz_xyzzz[j] - fr * g2_xzz_yzzz[j];

                g_xxzz_zzzz[j] = g1_xzz_xzzzz[j] - fr * g2_xzz_zzzz[j];

                g_xyyy_xxxx[j] = g1_yyy_xxxxx[j] - fr * g2_yyy_xxxx[j];

                g_xyyy_xxxy[j] = g1_yyy_xxxxy[j] - fr * g2_yyy_xxxy[j];

                g_xyyy_xxxz[j] = g1_yyy_xxxxz[j] - fr * g2_yyy_xxxz[j];

                g_xyyy_xxyy[j] = g1_yyy_xxxyy[j] - fr * g2_yyy_xxyy[j];

                g_xyyy_xxyz[j] = g1_yyy_xxxyz[j] - fr * g2_yyy_xxyz[j];

                g_xyyy_xxzz[j] = g1_yyy_xxxzz[j] - fr * g2_yyy_xxzz[j];

                g_xyyy_xyyy[j] = g1_yyy_xxyyy[j] - fr * g2_yyy_xyyy[j];

                g_xyyy_xyyz[j] = g1_yyy_xxyyz[j] - fr * g2_yyy_xyyz[j];

                g_xyyy_xyzz[j] = g1_yyy_xxyzz[j] - fr * g2_yyy_xyzz[j];

                g_xyyy_xzzz[j] = g1_yyy_xxzzz[j] - fr * g2_yyy_xzzz[j];

                g_xyyy_yyyy[j] = g1_yyy_xyyyy[j] - fr * g2_yyy_yyyy[j];

                g_xyyy_yyyz[j] = g1_yyy_xyyyz[j] - fr * g2_yyy_yyyz[j];

                g_xyyy_yyzz[j] = g1_yyy_xyyzz[j] - fr * g2_yyy_yyzz[j];

                g_xyyy_yzzz[j] = g1_yyy_xyzzz[j] - fr * g2_yyy_yzzz[j];

                g_xyyy_zzzz[j] = g1_yyy_xzzzz[j] - fr * g2_yyy_zzzz[j];

                g_xyyz_xxxx[j] = g1_yyz_xxxxx[j] - fr * g2_yyz_xxxx[j];

                g_xyyz_xxxy[j] = g1_yyz_xxxxy[j] - fr * g2_yyz_xxxy[j];

                g_xyyz_xxxz[j] = g1_yyz_xxxxz[j] - fr * g2_yyz_xxxz[j];

                g_xyyz_xxyy[j] = g1_yyz_xxxyy[j] - fr * g2_yyz_xxyy[j];

                g_xyyz_xxyz[j] = g1_yyz_xxxyz[j] - fr * g2_yyz_xxyz[j];

                g_xyyz_xxzz[j] = g1_yyz_xxxzz[j] - fr * g2_yyz_xxzz[j];

                g_xyyz_xyyy[j] = g1_yyz_xxyyy[j] - fr * g2_yyz_xyyy[j];

                g_xyyz_xyyz[j] = g1_yyz_xxyyz[j] - fr * g2_yyz_xyyz[j];

                g_xyyz_xyzz[j] = g1_yyz_xxyzz[j] - fr * g2_yyz_xyzz[j];

                g_xyyz_xzzz[j] = g1_yyz_xxzzz[j] - fr * g2_yyz_xzzz[j];

                g_xyyz_yyyy[j] = g1_yyz_xyyyy[j] - fr * g2_yyz_yyyy[j];

                g_xyyz_yyyz[j] = g1_yyz_xyyyz[j] - fr * g2_yyz_yyyz[j];

                g_xyyz_yyzz[j] = g1_yyz_xyyzz[j] - fr * g2_yyz_yyzz[j];

                g_xyyz_yzzz[j] = g1_yyz_xyzzz[j] - fr * g2_yyz_yzzz[j];

                g_xyyz_zzzz[j] = g1_yyz_xzzzz[j] - fr * g2_yyz_zzzz[j];

                g_xyzz_xxxx[j] = g1_yzz_xxxxx[j] - fr * g2_yzz_xxxx[j];

                g_xyzz_xxxy[j] = g1_yzz_xxxxy[j] - fr * g2_yzz_xxxy[j];

                g_xyzz_xxxz[j] = g1_yzz_xxxxz[j] - fr * g2_yzz_xxxz[j];

                g_xyzz_xxyy[j] = g1_yzz_xxxyy[j] - fr * g2_yzz_xxyy[j];

                g_xyzz_xxyz[j] = g1_yzz_xxxyz[j] - fr * g2_yzz_xxyz[j];

                g_xyzz_xxzz[j] = g1_yzz_xxxzz[j] - fr * g2_yzz_xxzz[j];

                g_xyzz_xyyy[j] = g1_yzz_xxyyy[j] - fr * g2_yzz_xyyy[j];

                g_xyzz_xyyz[j] = g1_yzz_xxyyz[j] - fr * g2_yzz_xyyz[j];

                g_xyzz_xyzz[j] = g1_yzz_xxyzz[j] - fr * g2_yzz_xyzz[j];

                g_xyzz_xzzz[j] = g1_yzz_xxzzz[j] - fr * g2_yzz_xzzz[j];

                g_xyzz_yyyy[j] = g1_yzz_xyyyy[j] - fr * g2_yzz_yyyy[j];

                g_xyzz_yyyz[j] = g1_yzz_xyyyz[j] - fr * g2_yzz_yyyz[j];

                g_xyzz_yyzz[j] = g1_yzz_xyyzz[j] - fr * g2_yzz_yyzz[j];

                g_xyzz_yzzz[j] = g1_yzz_xyzzz[j] - fr * g2_yzz_yzzz[j];

                g_xyzz_zzzz[j] = g1_yzz_xzzzz[j] - fr * g2_yzz_zzzz[j];

                g_xzzz_xxxx[j] = g1_zzz_xxxxx[j] - fr * g2_zzz_xxxx[j];

                g_xzzz_xxxy[j] = g1_zzz_xxxxy[j] - fr * g2_zzz_xxxy[j];

                g_xzzz_xxxz[j] = g1_zzz_xxxxz[j] - fr * g2_zzz_xxxz[j];

                g_xzzz_xxyy[j] = g1_zzz_xxxyy[j] - fr * g2_zzz_xxyy[j];

                g_xzzz_xxyz[j] = g1_zzz_xxxyz[j] - fr * g2_zzz_xxyz[j];

                g_xzzz_xxzz[j] = g1_zzz_xxxzz[j] - fr * g2_zzz_xxzz[j];

                g_xzzz_xyyy[j] = g1_zzz_xxyyy[j] - fr * g2_zzz_xyyy[j];

                g_xzzz_xyyz[j] = g1_zzz_xxyyz[j] - fr * g2_zzz_xyyz[j];

                g_xzzz_xyzz[j] = g1_zzz_xxyzz[j] - fr * g2_zzz_xyzz[j];

                g_xzzz_xzzz[j] = g1_zzz_xxzzz[j] - fr * g2_zzz_xzzz[j];

                g_xzzz_yyyy[j] = g1_zzz_xyyyy[j] - fr * g2_zzz_yyyy[j];

                g_xzzz_yyyz[j] = g1_zzz_xyyyz[j] - fr * g2_zzz_yyyz[j];

                g_xzzz_yyzz[j] = g1_zzz_xyyzz[j] - fr * g2_zzz_yyzz[j];

                g_xzzz_yzzz[j] = g1_zzz_xyzzz[j] - fr * g2_zzz_yzzz[j];

                g_xzzz_zzzz[j] = g1_zzz_xzzzz[j] - fr * g2_zzz_zzzz[j];

                // leading y component

                fr = rcdy[j];

                g_yyyy_xxxx[j] = g1_yyy_xxxxy[j] - fr * g2_yyy_xxxx[j];

                g_yyyy_xxxy[j] = g1_yyy_xxxyy[j] - fr * g2_yyy_xxxy[j];

                g_yyyy_xxxz[j] = g1_yyy_xxxyz[j] - fr * g2_yyy_xxxz[j];

                g_yyyy_xxyy[j] = g1_yyy_xxyyy[j] - fr * g2_yyy_xxyy[j];

                g_yyyy_xxyz[j] = g1_yyy_xxyyz[j] - fr * g2_yyy_xxyz[j];

                g_yyyy_xxzz[j] = g1_yyy_xxyzz[j] - fr * g2_yyy_xxzz[j];

                g_yyyy_xyyy[j] = g1_yyy_xyyyy[j] - fr * g2_yyy_xyyy[j];

                g_yyyy_xyyz[j] = g1_yyy_xyyyz[j] - fr * g2_yyy_xyyz[j];

                g_yyyy_xyzz[j] = g1_yyy_xyyzz[j] - fr * g2_yyy_xyzz[j];

                g_yyyy_xzzz[j] = g1_yyy_xyzzz[j] - fr * g2_yyy_xzzz[j];

                g_yyyy_yyyy[j] = g1_yyy_yyyyy[j] - fr * g2_yyy_yyyy[j];

                g_yyyy_yyyz[j] = g1_yyy_yyyyz[j] - fr * g2_yyy_yyyz[j];

                g_yyyy_yyzz[j] = g1_yyy_yyyzz[j] - fr * g2_yyy_yyzz[j];

                g_yyyy_yzzz[j] = g1_yyy_yyzzz[j] - fr * g2_yyy_yzzz[j];

                g_yyyy_zzzz[j] = g1_yyy_yzzzz[j] - fr * g2_yyy_zzzz[j];

                g_yyyz_xxxx[j] = g1_yyz_xxxxy[j] - fr * g2_yyz_xxxx[j];

                g_yyyz_xxxy[j] = g1_yyz_xxxyy[j] - fr * g2_yyz_xxxy[j];

                g_yyyz_xxxz[j] = g1_yyz_xxxyz[j] - fr * g2_yyz_xxxz[j];

                g_yyyz_xxyy[j] = g1_yyz_xxyyy[j] - fr * g2_yyz_xxyy[j];

                g_yyyz_xxyz[j] = g1_yyz_xxyyz[j] - fr * g2_yyz_xxyz[j];

                g_yyyz_xxzz[j] = g1_yyz_xxyzz[j] - fr * g2_yyz_xxzz[j];

                g_yyyz_xyyy[j] = g1_yyz_xyyyy[j] - fr * g2_yyz_xyyy[j];

                g_yyyz_xyyz[j] = g1_yyz_xyyyz[j] - fr * g2_yyz_xyyz[j];

                g_yyyz_xyzz[j] = g1_yyz_xyyzz[j] - fr * g2_yyz_xyzz[j];

                g_yyyz_xzzz[j] = g1_yyz_xyzzz[j] - fr * g2_yyz_xzzz[j];

                g_yyyz_yyyy[j] = g1_yyz_yyyyy[j] - fr * g2_yyz_yyyy[j];

                g_yyyz_yyyz[j] = g1_yyz_yyyyz[j] - fr * g2_yyz_yyyz[j];

                g_yyyz_yyzz[j] = g1_yyz_yyyzz[j] - fr * g2_yyz_yyzz[j];

                g_yyyz_yzzz[j] = g1_yyz_yyzzz[j] - fr * g2_yyz_yzzz[j];

                g_yyyz_zzzz[j] = g1_yyz_yzzzz[j] - fr * g2_yyz_zzzz[j];

                g_yyzz_xxxx[j] = g1_yzz_xxxxy[j] - fr * g2_yzz_xxxx[j];

                g_yyzz_xxxy[j] = g1_yzz_xxxyy[j] - fr * g2_yzz_xxxy[j];

                g_yyzz_xxxz[j] = g1_yzz_xxxyz[j] - fr * g2_yzz_xxxz[j];

                g_yyzz_xxyy[j] = g1_yzz_xxyyy[j] - fr * g2_yzz_xxyy[j];

                g_yyzz_xxyz[j] = g1_yzz_xxyyz[j] - fr * g2_yzz_xxyz[j];

                g_yyzz_xxzz[j] = g1_yzz_xxyzz[j] - fr * g2_yzz_xxzz[j];

                g_yyzz_xyyy[j] = g1_yzz_xyyyy[j] - fr * g2_yzz_xyyy[j];

                g_yyzz_xyyz[j] = g1_yzz_xyyyz[j] - fr * g2_yzz_xyyz[j];

                g_yyzz_xyzz[j] = g1_yzz_xyyzz[j] - fr * g2_yzz_xyzz[j];

                g_yyzz_xzzz[j] = g1_yzz_xyzzz[j] - fr * g2_yzz_xzzz[j];

                g_yyzz_yyyy[j] = g1_yzz_yyyyy[j] - fr * g2_yzz_yyyy[j];

                g_yyzz_yyyz[j] = g1_yzz_yyyyz[j] - fr * g2_yzz_yyyz[j];

                g_yyzz_yyzz[j] = g1_yzz_yyyzz[j] - fr * g2_yzz_yyzz[j];

                g_yyzz_yzzz[j] = g1_yzz_yyzzz[j] - fr * g2_yzz_yzzz[j];

                g_yyzz_zzzz[j] = g1_yzz_yzzzz[j] - fr * g2_yzz_zzzz[j];

                g_yzzz_xxxx[j] = g1_zzz_xxxxy[j] - fr * g2_zzz_xxxx[j];

                g_yzzz_xxxy[j] = g1_zzz_xxxyy[j] - fr * g2_zzz_xxxy[j];

                g_yzzz_xxxz[j] = g1_zzz_xxxyz[j] - fr * g2_zzz_xxxz[j];

                g_yzzz_xxyy[j] = g1_zzz_xxyyy[j] - fr * g2_zzz_xxyy[j];

                g_yzzz_xxyz[j] = g1_zzz_xxyyz[j] - fr * g2_zzz_xxyz[j];

                g_yzzz_xxzz[j] = g1_zzz_xxyzz[j] - fr * g2_zzz_xxzz[j];

                g_yzzz_xyyy[j] = g1_zzz_xyyyy[j] - fr * g2_zzz_xyyy[j];

                g_yzzz_xyyz[j] = g1_zzz_xyyyz[j] - fr * g2_zzz_xyyz[j];

                g_yzzz_xyzz[j] = g1_zzz_xyyzz[j] - fr * g2_zzz_xyzz[j];

                g_yzzz_xzzz[j] = g1_zzz_xyzzz[j] - fr * g2_zzz_xzzz[j];

                g_yzzz_yyyy[j] = g1_zzz_yyyyy[j] - fr * g2_zzz_yyyy[j];

                g_yzzz_yyyz[j] = g1_zzz_yyyyz[j] - fr * g2_zzz_yyyz[j];

                g_yzzz_yyzz[j] = g1_zzz_yyyzz[j] - fr * g2_zzz_yyzz[j];

                g_yzzz_yzzz[j] = g1_zzz_yyzzz[j] - fr * g2_zzz_yzzz[j];

                g_yzzz_zzzz[j] = g1_zzz_yzzzz[j] - fr * g2_zzz_zzzz[j];

                // leading z component

                fr = rcdz[j];

                g_zzzz_xxxx[j] = g1_zzz_xxxxz[j] - fr * g2_zzz_xxxx[j];

                g_zzzz_xxxy[j] = g1_zzz_xxxyz[j] - fr * g2_zzz_xxxy[j];

                g_zzzz_xxxz[j] = g1_zzz_xxxzz[j] - fr * g2_zzz_xxxz[j];

                g_zzzz_xxyy[j] = g1_zzz_xxyyz[j] - fr * g2_zzz_xxyy[j];

                g_zzzz_xxyz[j] = g1_zzz_xxyzz[j] - fr * g2_zzz_xxyz[j];

                g_zzzz_xxzz[j] = g1_zzz_xxzzz[j] - fr * g2_zzz_xxzz[j];

                g_zzzz_xyyy[j] = g1_zzz_xyyyz[j] - fr * g2_zzz_xyyy[j];

                g_zzzz_xyyz[j] = g1_zzz_xyyzz[j] - fr * g2_zzz_xyyz[j];

                g_zzzz_xyzz[j] = g1_zzz_xyzzz[j] - fr * g2_zzz_xyzz[j];

                g_zzzz_xzzz[j] = g1_zzz_xzzzz[j] - fr * g2_zzz_xzzz[j];

                g_zzzz_yyyy[j] = g1_zzz_yyyyz[j] - fr * g2_zzz_yyyy[j];

                g_zzzz_yyyz[j] = g1_zzz_yyyzz[j] - fr * g2_zzz_yyyz[j];

                g_zzzz_yyzz[j] = g1_zzz_yyzzz[j] - fr * g2_zzz_yyzz[j];

                g_zzzz_yzzz[j] = g1_zzz_yzzzz[j] - fr * g2_zzz_yzzz[j];

                g_zzzz_zzzz[j] = g1_zzz_zzzzz[j] - fr * g2_zzz_zzzz[j];
            }
        }
    }
    
} // t3hrrfunc namespace
