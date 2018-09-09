//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KetHrrFunc.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

namespace kethrrfunc { // kethrrfunc namespace
    
    void
    compElectronRepulsionForSXPP(      CMemBlock2D<double>&  ketBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side
        
        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
        
        if (isBraEqualKet) kdim  = iContrPair + 1;
        
        // set up pointers to distances R(CD) = C - D
        
        auto rcdx = cdDistances.data(0);
        
        auto rcdy = cdDistances.data(1);
        
        auto rcdz = cdDistances.data(2);
        
        // loop over set of data vectors
        
        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 1))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|11)\n");
             
                // determine angular momentum of bra side
                
                auto bang = recPattern[i].first();
                
                auto bcomp = angmom::to_CartesianComponents(bang);
                
                // get position of integrals in integrals buffer
                
                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 1});
                
                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 2});
                
                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 1});
                
                // compute contracted integrals
                
                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SP) integrals
                    
                    auto g2_0_x = ketBuffer.data(g2off + 3 * j);
                    
                    auto g2_0_y = ketBuffer.data(g2off + 3 * j + 1);
                    
                    auto g2_0_z = ketBuffer.data(g2off + 3 * j + 2);
                    
                    // set up pointers to (SX|g(r,r')|SD) integrals
                    
                    auto g1_0_xx = ketBuffer.data(g1off + 6 * j);
                    
                    auto g1_0_xy = ketBuffer.data(g1off + 6 * j + 1);
                    
                    auto g1_0_xz = ketBuffer.data(g1off + 6 * j + 2);
                    
                    auto g1_0_yy = ketBuffer.data(g1off + 6 * j + 3);
                    
                    auto g1_0_yz = ketBuffer.data(g1off + 6 * j + 4);
                    
                    auto g1_0_zz = ketBuffer.data(g1off + 6 * j + 5);
                    
                    // set up pointers to (SX|g(r,r')|PP) integrals
                    
                    auto g_x_x = ketBuffer.data(goff + 9 * j);
                    
                    auto g_x_y = ketBuffer.data(goff + 9 * j + 1);
                    
                    auto g_x_z = ketBuffer.data(goff + 9 * j + 2);
                    
                    auto g_y_x = ketBuffer.data(goff + 9 * j + 3);
                    
                    auto g_y_y = ketBuffer.data(goff + 9 * j + 4);
                    
                    auto g_y_z = ketBuffer.data(goff + 9 * j + 5);
                    
                    auto g_z_x = ketBuffer.data(goff + 9 * j + 6);
                    
                    auto g_z_y = ketBuffer.data(goff + 9 * j + 7);
                    
                    auto g_z_z = ketBuffer.data(goff + 9 * j + 8);
                    
                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_x, g2_0_y, g2_0_z,\
                                             g1_0_xx, g1_0_xy, g1_0_xz, g1_0_yy, g1_0_yz,\
                                             g1_0_zz, g_x_x, g_x_y, g_x_z, g_y_x, g_y_y,\
                                             g_y_z, g_z_x, g_z_y, g_z_z: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component
                        
                        auto fr = rcdx[k];
                        
                        g_x_x[k] = g1_0_xx[k] - fr * g2_0_x[k];
                        
                        g_x_y[k] = g1_0_xy[k] - fr * g2_0_y[k];
                        
                        g_x_z[k] = g1_0_xz[k] - fr * g2_0_z[k];
                        
                        // leading y component
                        
                        fr = rcdy[k];
                        
                        g_y_x[k] = g1_0_xy[k] - fr * g2_0_x[k];
                        
                        g_y_y[k] = g1_0_yy[k] - fr * g2_0_y[k];
                        
                        g_y_z[k] = g1_0_yz[k] - fr * g2_0_z[k];
                        
                        // leading z component
                        
                        fr = rcdz[k];
                        
                        g_z_x[k] = g1_0_xz[k] - fr * g2_0_x[k];
                        
                        g_z_y[k] = g1_0_yz[k] - fr * g2_0_y[k];
                        
                        g_z_z[k] = g1_0_zz[k] - fr * g2_0_z[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXPD(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 2))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|12)\n");
                
                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                         {bang, 1, 2});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                          {bang, 0, 3});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                          {bang, 0, 2});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SD)^(m) integrals

                    auto g2_0_xx = contrBuffer.data(g2off + 6 * j);

                    auto g2_0_xy = contrBuffer.data(g2off + 6 * j + 1);

                    auto g2_0_xz = contrBuffer.data(g2off + 6 * j + 2);

                    auto g2_0_yy = contrBuffer.data(g2off + 6 * j + 3);

                    auto g2_0_yz = contrBuffer.data(g2off + 6 * j + 4);

                    auto g2_0_zz = contrBuffer.data(g2off + 6 * j + 5);

                    // set up pointers to (SX|g(r,r')|SF)^(m) integrals

                    auto g1_0_xxx = contrBuffer.data(g1off + 10 * j);

                    auto g1_0_xxy = contrBuffer.data(g1off + 10 * j + 1);

                    auto g1_0_xxz = contrBuffer.data(g1off + 10 * j + 2);

                    auto g1_0_xyy = contrBuffer.data(g1off + 10 * j + 3);

                    auto g1_0_xyz = contrBuffer.data(g1off + 10 * j + 4);

                    auto g1_0_xzz = contrBuffer.data(g1off + 10 * j + 5);

                    auto g1_0_yyy = contrBuffer.data(g1off + 10 * j + 6);

                    auto g1_0_yyz = contrBuffer.data(g1off + 10 * j + 7);

                    auto g1_0_yzz = contrBuffer.data(g1off + 10 * j + 8);

                    auto g1_0_zzz = contrBuffer.data(g1off + 10 * j + 9);

                    // set up pointers to (SX|g(r,r')|PD)^(m) integrals

                    auto g_x_xx = contrBuffer.data(goff + 18 * j);

                    auto g_x_xy = contrBuffer.data(goff + 18 * j + 1);

                    auto g_x_xz = contrBuffer.data(goff + 18 * j + 2);

                    auto g_x_yy = contrBuffer.data(goff + 18 * j + 3);

                    auto g_x_yz = contrBuffer.data(goff + 18 * j + 4);

                    auto g_x_zz = contrBuffer.data(goff + 18 * j + 5);

                    auto g_y_xx = contrBuffer.data(goff + 18 * j + 6);

                    auto g_y_xy = contrBuffer.data(goff + 18 * j + 7);

                    auto g_y_xz = contrBuffer.data(goff + 18 * j + 8);

                    auto g_y_yy = contrBuffer.data(goff + 18 * j + 9);

                    auto g_y_yz = contrBuffer.data(goff + 18 * j + 10);

                    auto g_y_zz = contrBuffer.data(goff + 18 * j + 11);

                    auto g_z_xx = contrBuffer.data(goff + 18 * j + 12);

                    auto g_z_xy = contrBuffer.data(goff + 18 * j + 13);

                    auto g_z_xz = contrBuffer.data(goff + 18 * j + 14);

                    auto g_z_yy = contrBuffer.data(goff + 18 * j + 15);

                    auto g_z_yz = contrBuffer.data(goff + 18 * j + 16);

                    auto g_z_zz = contrBuffer.data(goff + 18 * j + 17);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xx, g2_0_xy,\
                                             g2_0_xz, g2_0_yy, g2_0_yz, g2_0_zz,\
                                             g1_0_xxx, g1_0_xxy, g1_0_xxz,\
                                             g1_0_xyy, g1_0_xyz, g1_0_xzz,\
                                             g1_0_yyy, g1_0_yyz, g1_0_yzz,\
                                             g1_0_zzz, g_x_xx, g_x_xy, g_x_xz,\
                                             g_x_yy, g_x_yz, g_x_zz, g_y_xx,\
                                             g_y_xy, g_y_xz, g_y_yy, g_y_yz,\
                                             g_y_zz, g_z_xx, g_z_xy, g_z_xz,\
                                             g_z_yy, g_z_yz, g_z_zz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xx[k] = g1_0_xxx[k] - fr * g2_0_xx[k];

                        g_x_xy[k] = g1_0_xxy[k] - fr * g2_0_xy[k];

                        g_x_xz[k] = g1_0_xxz[k] - fr * g2_0_xz[k];

                        g_x_yy[k] = g1_0_xyy[k] - fr * g2_0_yy[k];

                        g_x_yz[k] = g1_0_xyz[k] - fr * g2_0_yz[k];

                        g_x_zz[k] = g1_0_xzz[k] - fr * g2_0_zz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xx[k] = g1_0_xxy[k] - fr * g2_0_xx[k];

                        g_y_xy[k] = g1_0_xyy[k] - fr * g2_0_xy[k];

                        g_y_xz[k] = g1_0_xyz[k] - fr * g2_0_xz[k];

                        g_y_yy[k] = g1_0_yyy[k] - fr * g2_0_yy[k];

                        g_y_yz[k] = g1_0_yyz[k] - fr * g2_0_yz[k];

                        g_y_zz[k] = g1_0_yzz[k] - fr * g2_0_zz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xx[k] = g1_0_xxz[k] - fr * g2_0_xx[k];

                        g_z_xy[k] = g1_0_xyz[k] - fr * g2_0_xy[k];

                        g_z_xz[k] = g1_0_xzz[k] - fr * g2_0_xz[k];

                        g_z_yy[k] = g1_0_yyz[k] - fr * g2_0_yy[k];

                        g_z_yz[k] = g1_0_yzz[k] - fr * g2_0_yz[k];

                        g_z_zz[k] = g1_0_zzz[k] - fr * g2_0_zz[k];
                    }
                }
            }
        }
    }

    void
    compElectronRepulsionForSXPF(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 3))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|13)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 3});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 4});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 3});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SF)^(m) integrals

                    auto g2_0_xxx = contrBuffer.data(g2off + 10 * j);

                    auto g2_0_xxy = contrBuffer.data(g2off + 10 * j + 1);

                    auto g2_0_xxz = contrBuffer.data(g2off + 10 * j + 2);

                    auto g2_0_xyy = contrBuffer.data(g2off + 10 * j + 3);

                    auto g2_0_xyz = contrBuffer.data(g2off + 10 * j + 4);

                    auto g2_0_xzz = contrBuffer.data(g2off + 10 * j + 5);

                    auto g2_0_yyy = contrBuffer.data(g2off + 10 * j + 6);

                    auto g2_0_yyz = contrBuffer.data(g2off + 10 * j + 7);

                    auto g2_0_yzz = contrBuffer.data(g2off + 10 * j + 8);

                    auto g2_0_zzz = contrBuffer.data(g2off + 10 * j + 9);

                    // set up pointers to (SX|g(r,r')|SG)^(m) integrals

                    auto g1_0_xxxx = contrBuffer.data(g1off + 15 * j);

                    auto g1_0_xxxy = contrBuffer.data(g1off + 15 * j + 1);

                    auto g1_0_xxxz = contrBuffer.data(g1off + 15 * j + 2);

                    auto g1_0_xxyy = contrBuffer.data(g1off + 15 * j + 3);

                    auto g1_0_xxyz = contrBuffer.data(g1off + 15 * j + 4);

                    auto g1_0_xxzz = contrBuffer.data(g1off + 15 * j + 5);

                    auto g1_0_xyyy = contrBuffer.data(g1off + 15 * j + 6);

                    auto g1_0_xyyz = contrBuffer.data(g1off + 15 * j + 7);

                    auto g1_0_xyzz = contrBuffer.data(g1off + 15 * j + 8);

                    auto g1_0_xzzz = contrBuffer.data(g1off + 15 * j + 9);

                    auto g1_0_yyyy = contrBuffer.data(g1off + 15 * j + 10);

                    auto g1_0_yyyz = contrBuffer.data(g1off + 15 * j + 11);

                    auto g1_0_yyzz = contrBuffer.data(g1off + 15 * j + 12);

                    auto g1_0_yzzz = contrBuffer.data(g1off + 15 * j + 13);

                    auto g1_0_zzzz = contrBuffer.data(g1off + 15 * j + 14);

                    // set up pointers to (SX|g(r,r')|PF)^(m) integrals

                    auto g_x_xxx = contrBuffer.data(goff + 30 * j);

                    auto g_x_xxy = contrBuffer.data(goff + 30 * j + 1);

                    auto g_x_xxz = contrBuffer.data(goff + 30 * j + 2);

                    auto g_x_xyy = contrBuffer.data(goff + 30 * j + 3);

                    auto g_x_xyz = contrBuffer.data(goff + 30 * j + 4);

                    auto g_x_xzz = contrBuffer.data(goff + 30 * j + 5);

                    auto g_x_yyy = contrBuffer.data(goff + 30 * j + 6);

                    auto g_x_yyz = contrBuffer.data(goff + 30 * j + 7);

                    auto g_x_yzz = contrBuffer.data(goff + 30 * j + 8);

                    auto g_x_zzz = contrBuffer.data(goff + 30 * j + 9);

                    auto g_y_xxx = contrBuffer.data(goff + 30 * j + 10);

                    auto g_y_xxy = contrBuffer.data(goff + 30 * j + 11);

                    auto g_y_xxz = contrBuffer.data(goff + 30 * j + 12);

                    auto g_y_xyy = contrBuffer.data(goff + 30 * j + 13);

                    auto g_y_xyz = contrBuffer.data(goff + 30 * j + 14);

                    auto g_y_xzz = contrBuffer.data(goff + 30 * j + 15);

                    auto g_y_yyy = contrBuffer.data(goff + 30 * j + 16);

                    auto g_y_yyz = contrBuffer.data(goff + 30 * j + 17);

                    auto g_y_yzz = contrBuffer.data(goff + 30 * j + 18);

                    auto g_y_zzz = contrBuffer.data(goff + 30 * j + 19);

                    auto g_z_xxx = contrBuffer.data(goff + 30 * j + 20);

                    auto g_z_xxy = contrBuffer.data(goff + 30 * j + 21);

                    auto g_z_xxz = contrBuffer.data(goff + 30 * j + 22);

                    auto g_z_xyy = contrBuffer.data(goff + 30 * j + 23);

                    auto g_z_xyz = contrBuffer.data(goff + 30 * j + 24);

                    auto g_z_xzz = contrBuffer.data(goff + 30 * j + 25);

                    auto g_z_yyy = contrBuffer.data(goff + 30 * j + 26);

                    auto g_z_yyz = contrBuffer.data(goff + 30 * j + 27);

                    auto g_z_yzz = contrBuffer.data(goff + 30 * j + 28);

                    auto g_z_zzz = contrBuffer.data(goff + 30 * j + 29);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxx, g2_0_xxy,\
                                             g2_0_xxz, g2_0_xyy, g2_0_xyz, g2_0_xzz,\
                                             g2_0_yyy, g2_0_yyz, g2_0_yzz, g2_0_zzz,\
                                             g1_0_xxxx, g1_0_xxxy, g1_0_xxxz, g1_0_xxyy,\
                                             g1_0_xxyz, g1_0_xxzz, g1_0_xyyy, g1_0_xyyz,\
                                             g1_0_xyzz, g1_0_xzzz, g1_0_yyyy, g1_0_yyyz,\
                                             g1_0_yyzz, g1_0_yzzz, g1_0_zzzz, g_x_xxx,\
                                             g_x_xxy, g_x_xxz, g_x_xyy, g_x_xyz,\
                                             g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz,\
                                             g_x_zzz, g_y_xxx, g_y_xxy, g_y_xxz,\
                                             g_y_xyy, g_y_xyz, g_y_xzz, g_y_yyy,\
                                             g_y_yyz, g_y_yzz, g_y_zzz, g_z_xxx,\
                                             g_z_xxy, g_z_xxz, g_z_xyy, g_z_xyz,\
                                             g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz,\
                                             g_z_zzz: VLX_ALIGN)
                     for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xxx[k] = g1_0_xxxx[k] - fr * g2_0_xxx[k];

                        g_x_xxy[k] = g1_0_xxxy[k] - fr * g2_0_xxy[k];

                        g_x_xxz[k] = g1_0_xxxz[k] - fr * g2_0_xxz[k];

                        g_x_xyy[k] = g1_0_xxyy[k] - fr * g2_0_xyy[k];

                        g_x_xyz[k] = g1_0_xxyz[k] - fr * g2_0_xyz[k];

                        g_x_xzz[k] = g1_0_xxzz[k] - fr * g2_0_xzz[k];

                        g_x_yyy[k] = g1_0_xyyy[k] - fr * g2_0_yyy[k];

                        g_x_yyz[k] = g1_0_xyyz[k] - fr * g2_0_yyz[k];

                        g_x_yzz[k] = g1_0_xyzz[k] - fr * g2_0_yzz[k];

                        g_x_zzz[k] = g1_0_xzzz[k] - fr * g2_0_zzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xxx[k] = g1_0_xxxy[k] - fr * g2_0_xxx[k];

                        g_y_xxy[k] = g1_0_xxyy[k] - fr * g2_0_xxy[k];

                        g_y_xxz[k] = g1_0_xxyz[k] - fr * g2_0_xxz[k];

                        g_y_xyy[k] = g1_0_xyyy[k] - fr * g2_0_xyy[k];

                        g_y_xyz[k] = g1_0_xyyz[k] - fr * g2_0_xyz[k];

                        g_y_xzz[k] = g1_0_xyzz[k] - fr * g2_0_xzz[k];

                        g_y_yyy[k] = g1_0_yyyy[k] - fr * g2_0_yyy[k];

                        g_y_yyz[k] = g1_0_yyyz[k] - fr * g2_0_yyz[k];

                        g_y_yzz[k] = g1_0_yyzz[k] - fr * g2_0_yzz[k];

                        g_y_zzz[k] = g1_0_yzzz[k] - fr * g2_0_zzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xxx[k] = g1_0_xxxz[k] - fr * g2_0_xxx[k];

                        g_z_xxy[k] = g1_0_xxyz[k] - fr * g2_0_xxy[k];

                        g_z_xxz[k] = g1_0_xxzz[k] - fr * g2_0_xxz[k];

                        g_z_xyy[k] = g1_0_xyyz[k] - fr * g2_0_xyy[k];

                        g_z_xyz[k] = g1_0_xyzz[k] - fr * g2_0_xyz[k];

                        g_z_xzz[k] = g1_0_xzzz[k] - fr * g2_0_xzz[k];

                        g_z_yyy[k] = g1_0_yyyz[k] - fr * g2_0_yyy[k];

                        g_z_yyz[k] = g1_0_yyzz[k] - fr * g2_0_yyz[k];

                        g_z_yzz[k] = g1_0_yzzz[k] - fr * g2_0_yzz[k];

                        g_z_zzz[k] = g1_0_zzzz[k] - fr * g2_0_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXPG(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 4))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|14)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 4});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 5});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 4});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SG)^(m) integrals

                    auto g2_0_xxxx = contrBuffer.data(g2off + 15 * j);

                    auto g2_0_xxxy = contrBuffer.data(g2off + 15 * j + 1);

                    auto g2_0_xxxz = contrBuffer.data(g2off + 15 * j + 2);

                    auto g2_0_xxyy = contrBuffer.data(g2off + 15 * j + 3);

                    auto g2_0_xxyz = contrBuffer.data(g2off + 15 * j + 4);

                    auto g2_0_xxzz = contrBuffer.data(g2off + 15 * j + 5);

                    auto g2_0_xyyy = contrBuffer.data(g2off + 15 * j + 6);

                    auto g2_0_xyyz = contrBuffer.data(g2off + 15 * j + 7);

                    auto g2_0_xyzz = contrBuffer.data(g2off + 15 * j + 8);

                    auto g2_0_xzzz = contrBuffer.data(g2off + 15 * j + 9);

                    auto g2_0_yyyy = contrBuffer.data(g2off + 15 * j + 10);

                    auto g2_0_yyyz = contrBuffer.data(g2off + 15 * j + 11);

                    auto g2_0_yyzz = contrBuffer.data(g2off + 15 * j + 12);

                    auto g2_0_yzzz = contrBuffer.data(g2off + 15 * j + 13);

                    auto g2_0_zzzz = contrBuffer.data(g2off + 15 * j + 14);

                    // set up pointers to (SX|g(r,r')|SH)^(m) integrals

                    auto g1_0_xxxxx = contrBuffer.data(g1off + 21 * j);

                    auto g1_0_xxxxy = contrBuffer.data(g1off + 21 * j + 1);

                    auto g1_0_xxxxz = contrBuffer.data(g1off + 21 * j + 2);

                    auto g1_0_xxxyy = contrBuffer.data(g1off + 21 * j + 3);

                    auto g1_0_xxxyz = contrBuffer.data(g1off + 21 * j + 4);

                    auto g1_0_xxxzz = contrBuffer.data(g1off + 21 * j + 5);

                    auto g1_0_xxyyy = contrBuffer.data(g1off + 21 * j + 6);

                    auto g1_0_xxyyz = contrBuffer.data(g1off + 21 * j + 7);

                    auto g1_0_xxyzz = contrBuffer.data(g1off + 21 * j + 8);

                    auto g1_0_xxzzz = contrBuffer.data(g1off + 21 * j + 9);

                    auto g1_0_xyyyy = contrBuffer.data(g1off + 21 * j + 10);

                    auto g1_0_xyyyz = contrBuffer.data(g1off + 21 * j + 11);

                    auto g1_0_xyyzz = contrBuffer.data(g1off + 21 * j + 12);

                    auto g1_0_xyzzz = contrBuffer.data(g1off + 21 * j + 13);

                    auto g1_0_xzzzz = contrBuffer.data(g1off + 21 * j + 14);

                    auto g1_0_yyyyy = contrBuffer.data(g1off + 21 * j + 15);

                    auto g1_0_yyyyz = contrBuffer.data(g1off + 21 * j + 16);

                    auto g1_0_yyyzz = contrBuffer.data(g1off + 21 * j + 17);

                    auto g1_0_yyzzz = contrBuffer.data(g1off + 21 * j + 18);

                    auto g1_0_yzzzz = contrBuffer.data(g1off + 21 * j + 19);

                    auto g1_0_zzzzz = contrBuffer.data(g1off + 21 * j + 20);

                    // set up pointers to (SX|g(r,r')|PG)^(m) integrals

                    auto g_x_xxxx = contrBuffer.data(goff + 45 * j);

                    auto g_x_xxxy = contrBuffer.data(goff + 45 * j + 1);

                    auto g_x_xxxz = contrBuffer.data(goff + 45 * j + 2);

                    auto g_x_xxyy = contrBuffer.data(goff + 45 * j + 3);

                    auto g_x_xxyz = contrBuffer.data(goff + 45 * j + 4);

                    auto g_x_xxzz = contrBuffer.data(goff + 45 * j + 5);

                    auto g_x_xyyy = contrBuffer.data(goff + 45 * j + 6);

                    auto g_x_xyyz = contrBuffer.data(goff + 45 * j + 7);

                    auto g_x_xyzz = contrBuffer.data(goff + 45 * j + 8);

                    auto g_x_xzzz = contrBuffer.data(goff + 45 * j + 9);

                    auto g_x_yyyy = contrBuffer.data(goff + 45 * j + 10);

                    auto g_x_yyyz = contrBuffer.data(goff + 45 * j + 11);

                    auto g_x_yyzz = contrBuffer.data(goff + 45 * j + 12);

                    auto g_x_yzzz = contrBuffer.data(goff + 45 * j + 13);

                    auto g_x_zzzz = contrBuffer.data(goff + 45 * j + 14);

                    auto g_y_xxxx = contrBuffer.data(goff + 45 * j + 15);

                    auto g_y_xxxy = contrBuffer.data(goff + 45 * j + 16);

                    auto g_y_xxxz = contrBuffer.data(goff + 45 * j + 17);

                    auto g_y_xxyy = contrBuffer.data(goff + 45 * j + 18);

                    auto g_y_xxyz = contrBuffer.data(goff + 45 * j + 19);

                    auto g_y_xxzz = contrBuffer.data(goff + 45 * j + 20);

                    auto g_y_xyyy = contrBuffer.data(goff + 45 * j + 21);

                    auto g_y_xyyz = contrBuffer.data(goff + 45 * j + 22);

                    auto g_y_xyzz = contrBuffer.data(goff + 45 * j + 23);

                    auto g_y_xzzz = contrBuffer.data(goff + 45 * j + 24);

                    auto g_y_yyyy = contrBuffer.data(goff + 45 * j + 25);

                    auto g_y_yyyz = contrBuffer.data(goff + 45 * j + 26);

                    auto g_y_yyzz = contrBuffer.data(goff + 45 * j + 27);

                    auto g_y_yzzz = contrBuffer.data(goff + 45 * j + 28);

                    auto g_y_zzzz = contrBuffer.data(goff + 45 * j + 29);

                    auto g_z_xxxx = contrBuffer.data(goff + 45 * j + 30);

                    auto g_z_xxxy = contrBuffer.data(goff + 45 * j + 31);

                    auto g_z_xxxz = contrBuffer.data(goff + 45 * j + 32);

                    auto g_z_xxyy = contrBuffer.data(goff + 45 * j + 33);

                    auto g_z_xxyz = contrBuffer.data(goff + 45 * j + 34);

                    auto g_z_xxzz = contrBuffer.data(goff + 45 * j + 35);

                    auto g_z_xyyy = contrBuffer.data(goff + 45 * j + 36);

                    auto g_z_xyyz = contrBuffer.data(goff + 45 * j + 37);

                    auto g_z_xyzz = contrBuffer.data(goff + 45 * j + 38);

                    auto g_z_xzzz = contrBuffer.data(goff + 45 * j + 39);

                    auto g_z_yyyy = contrBuffer.data(goff + 45 * j + 40);

                    auto g_z_yyyz = contrBuffer.data(goff + 45 * j + 41);

                    auto g_z_yyzz = contrBuffer.data(goff + 45 * j + 42);

                    auto g_z_yzzz = contrBuffer.data(goff + 45 * j + 43);

                    auto g_z_zzzz = contrBuffer.data(goff + 45 * j + 44);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxx, g2_0_xxxy,\
                                             g2_0_xxxz, g2_0_xxyy, g2_0_xxyz, g2_0_xxzz,\
                                             g2_0_xyyy, g2_0_xyyz, g2_0_xyzz, g2_0_xzzz,\
                                             g2_0_yyyy, g2_0_yyyz, g2_0_yyzz, g2_0_yzzz,\
                                             g2_0_zzzz, g1_0_xxxxx, g1_0_xxxxy,\
                                             g1_0_xxxxz, g1_0_xxxyy, g1_0_xxxyz,\
                                             g1_0_xxxzz, g1_0_xxyyy, g1_0_xxyyz,\
                                             g1_0_xxyzz, g1_0_xxzzz, g1_0_xyyyy,\
                                             g1_0_xyyyz, g1_0_xyyzz, g1_0_xyzzz,\
                                             g1_0_xzzzz, g1_0_yyyyy, g1_0_yyyyz,\
                                             g1_0_yyyzz, g1_0_yyzzz, g1_0_yzzzz,\
                                             g1_0_zzzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz,\
                                             g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy,\
                                             g_x_xyyz, g_x_xyzz, g_x_xzzz, g_x_yyyy,\
                                             g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz,\
                                             g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxyy,\
                                             g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz,\
                                             g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz,\
                                             g_y_yyzz, g_y_yzzz, g_y_zzzz, g_z_xxxx,\
                                             g_z_xxxy, g_z_xxxz, g_z_xxyy, g_z_xxyz,\
                                             g_z_xxzz, g_z_xyyy, g_z_xyyz, g_z_xyzz,\
                                             g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz,\
                                             g_z_yzzz, g_z_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xxxx[k] = g1_0_xxxxx[k] - fr * g2_0_xxxx[k];

                        g_x_xxxy[k] = g1_0_xxxxy[k] - fr * g2_0_xxxy[k];

                        g_x_xxxz[k] = g1_0_xxxxz[k] - fr * g2_0_xxxz[k];

                        g_x_xxyy[k] = g1_0_xxxyy[k] - fr * g2_0_xxyy[k];

                        g_x_xxyz[k] = g1_0_xxxyz[k] - fr * g2_0_xxyz[k];

                        g_x_xxzz[k] = g1_0_xxxzz[k] - fr * g2_0_xxzz[k];

                        g_x_xyyy[k] = g1_0_xxyyy[k] - fr * g2_0_xyyy[k];

                        g_x_xyyz[k] = g1_0_xxyyz[k] - fr * g2_0_xyyz[k];

                        g_x_xyzz[k] = g1_0_xxyzz[k] - fr * g2_0_xyzz[k];

                        g_x_xzzz[k] = g1_0_xxzzz[k] - fr * g2_0_xzzz[k];

                        g_x_yyyy[k] = g1_0_xyyyy[k] - fr * g2_0_yyyy[k];

                        g_x_yyyz[k] = g1_0_xyyyz[k] - fr * g2_0_yyyz[k];

                        g_x_yyzz[k] = g1_0_xyyzz[k] - fr * g2_0_yyzz[k];

                        g_x_yzzz[k] = g1_0_xyzzz[k] - fr * g2_0_yzzz[k];

                        g_x_zzzz[k] = g1_0_xzzzz[k] - fr * g2_0_zzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xxxx[k] = g1_0_xxxxy[k] - fr * g2_0_xxxx[k];

                        g_y_xxxy[k] = g1_0_xxxyy[k] - fr * g2_0_xxxy[k];

                        g_y_xxxz[k] = g1_0_xxxyz[k] - fr * g2_0_xxxz[k];

                        g_y_xxyy[k] = g1_0_xxyyy[k] - fr * g2_0_xxyy[k];

                        g_y_xxyz[k] = g1_0_xxyyz[k] - fr * g2_0_xxyz[k];

                        g_y_xxzz[k] = g1_0_xxyzz[k] - fr * g2_0_xxzz[k];

                        g_y_xyyy[k] = g1_0_xyyyy[k] - fr * g2_0_xyyy[k];

                        g_y_xyyz[k] = g1_0_xyyyz[k] - fr * g2_0_xyyz[k];

                        g_y_xyzz[k] = g1_0_xyyzz[k] - fr * g2_0_xyzz[k];

                        g_y_xzzz[k] = g1_0_xyzzz[k] - fr * g2_0_xzzz[k];

                        g_y_yyyy[k] = g1_0_yyyyy[k] - fr * g2_0_yyyy[k];

                        g_y_yyyz[k] = g1_0_yyyyz[k] - fr * g2_0_yyyz[k];

                        g_y_yyzz[k] = g1_0_yyyzz[k] - fr * g2_0_yyzz[k];

                        g_y_yzzz[k] = g1_0_yyzzz[k] - fr * g2_0_yzzz[k];

                        g_y_zzzz[k] = g1_0_yzzzz[k] - fr * g2_0_zzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xxxx[k] = g1_0_xxxxz[k] - fr * g2_0_xxxx[k];

                        g_z_xxxy[k] = g1_0_xxxyz[k] - fr * g2_0_xxxy[k];

                        g_z_xxxz[k] = g1_0_xxxzz[k] - fr * g2_0_xxxz[k];

                        g_z_xxyy[k] = g1_0_xxyyz[k] - fr * g2_0_xxyy[k];

                        g_z_xxyz[k] = g1_0_xxyzz[k] - fr * g2_0_xxyz[k];

                        g_z_xxzz[k] = g1_0_xxzzz[k] - fr * g2_0_xxzz[k];

                        g_z_xyyy[k] = g1_0_xyyyz[k] - fr * g2_0_xyyy[k];

                        g_z_xyyz[k] = g1_0_xyyzz[k] - fr * g2_0_xyyz[k];

                        g_z_xyzz[k] = g1_0_xyzzz[k] - fr * g2_0_xyzz[k];

                        g_z_xzzz[k] = g1_0_xzzzz[k] - fr * g2_0_xzzz[k];

                        g_z_yyyy[k] = g1_0_yyyyz[k] - fr * g2_0_yyyy[k];

                        g_z_yyyz[k] = g1_0_yyyzz[k] - fr * g2_0_yyyz[k];

                        g_z_yyzz[k] = g1_0_yyzzz[k] - fr * g2_0_yyzz[k];

                        g_z_yzzz[k] = g1_0_yzzzz[k] - fr * g2_0_yzzz[k];

                        g_z_zzzz[k] = g1_0_zzzzz[k] - fr * g2_0_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXPH(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 5))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|15)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 5});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 6});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 5});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SH)^(m) integrals

                    auto g2_0_xxxxx = contrBuffer.data(g2off + 21 * j);

                    auto g2_0_xxxxy = contrBuffer.data(g2off + 21 * j + 1);

                    auto g2_0_xxxxz = contrBuffer.data(g2off + 21 * j + 2);

                    auto g2_0_xxxyy = contrBuffer.data(g2off + 21 * j + 3);

                    auto g2_0_xxxyz = contrBuffer.data(g2off + 21 * j + 4);

                    auto g2_0_xxxzz = contrBuffer.data(g2off + 21 * j + 5);

                    auto g2_0_xxyyy = contrBuffer.data(g2off + 21 * j + 6);

                    auto g2_0_xxyyz = contrBuffer.data(g2off + 21 * j + 7);

                    auto g2_0_xxyzz = contrBuffer.data(g2off + 21 * j + 8);

                    auto g2_0_xxzzz = contrBuffer.data(g2off + 21 * j + 9);

                    auto g2_0_xyyyy = contrBuffer.data(g2off + 21 * j + 10);

                    auto g2_0_xyyyz = contrBuffer.data(g2off + 21 * j + 11);

                    auto g2_0_xyyzz = contrBuffer.data(g2off + 21 * j + 12);

                    auto g2_0_xyzzz = contrBuffer.data(g2off + 21 * j + 13);

                    auto g2_0_xzzzz = contrBuffer.data(g2off + 21 * j + 14);

                    auto g2_0_yyyyy = contrBuffer.data(g2off + 21 * j + 15);

                    auto g2_0_yyyyz = contrBuffer.data(g2off + 21 * j + 16);

                    auto g2_0_yyyzz = contrBuffer.data(g2off + 21 * j + 17);

                    auto g2_0_yyzzz = contrBuffer.data(g2off + 21 * j + 18);

                    auto g2_0_yzzzz = contrBuffer.data(g2off + 21 * j + 19);

                    auto g2_0_zzzzz = contrBuffer.data(g2off + 21 * j + 20);

                    // set up pointers to (SX|g(r,r')|SI)^(m) integrals

                    auto g1_0_xxxxxx = contrBuffer.data(g1off + 28 * j);

                    auto g1_0_xxxxxy = contrBuffer.data(g1off + 28 * j + 1);

                    auto g1_0_xxxxxz = contrBuffer.data(g1off + 28 * j + 2);

                    auto g1_0_xxxxyy = contrBuffer.data(g1off + 28 * j + 3);

                    auto g1_0_xxxxyz = contrBuffer.data(g1off + 28 * j + 4);

                    auto g1_0_xxxxzz = contrBuffer.data(g1off + 28 * j + 5);

                    auto g1_0_xxxyyy = contrBuffer.data(g1off + 28 * j + 6);

                    auto g1_0_xxxyyz = contrBuffer.data(g1off + 28 * j + 7);

                    auto g1_0_xxxyzz = contrBuffer.data(g1off + 28 * j + 8);

                    auto g1_0_xxxzzz = contrBuffer.data(g1off + 28 * j + 9);

                    auto g1_0_xxyyyy = contrBuffer.data(g1off + 28 * j + 10);

                    auto g1_0_xxyyyz = contrBuffer.data(g1off + 28 * j + 11);

                    auto g1_0_xxyyzz = contrBuffer.data(g1off + 28 * j + 12);

                    auto g1_0_xxyzzz = contrBuffer.data(g1off + 28 * j + 13);

                    auto g1_0_xxzzzz = contrBuffer.data(g1off + 28 * j + 14);

                    auto g1_0_xyyyyy = contrBuffer.data(g1off + 28 * j + 15);

                    auto g1_0_xyyyyz = contrBuffer.data(g1off + 28 * j + 16);

                    auto g1_0_xyyyzz = contrBuffer.data(g1off + 28 * j + 17);

                    auto g1_0_xyyzzz = contrBuffer.data(g1off + 28 * j + 18);

                    auto g1_0_xyzzzz = contrBuffer.data(g1off + 28 * j + 19);

                    auto g1_0_xzzzzz = contrBuffer.data(g1off + 28 * j + 20);

                    auto g1_0_yyyyyy = contrBuffer.data(g1off + 28 * j + 21);

                    auto g1_0_yyyyyz = contrBuffer.data(g1off + 28 * j + 22);

                    auto g1_0_yyyyzz = contrBuffer.data(g1off + 28 * j + 23);

                    auto g1_0_yyyzzz = contrBuffer.data(g1off + 28 * j + 24);

                    auto g1_0_yyzzzz = contrBuffer.data(g1off + 28 * j + 25);

                    auto g1_0_yzzzzz = contrBuffer.data(g1off + 28 * j + 26);

                    auto g1_0_zzzzzz = contrBuffer.data(g1off + 28 * j + 27);

                    // set up pointers to (SX|g(r,r')|PH)^(m) integrals

                    auto g_x_xxxxx = contrBuffer.data(goff + 63 * j);

                    auto g_x_xxxxy = contrBuffer.data(goff + 63 * j + 1);

                    auto g_x_xxxxz = contrBuffer.data(goff + 63 * j + 2);

                    auto g_x_xxxyy = contrBuffer.data(goff + 63 * j + 3);

                    auto g_x_xxxyz = contrBuffer.data(goff + 63 * j + 4);

                    auto g_x_xxxzz = contrBuffer.data(goff + 63 * j + 5);

                    auto g_x_xxyyy = contrBuffer.data(goff + 63 * j + 6);

                    auto g_x_xxyyz = contrBuffer.data(goff + 63 * j + 7);

                    auto g_x_xxyzz = contrBuffer.data(goff + 63 * j + 8);

                    auto g_x_xxzzz = contrBuffer.data(goff + 63 * j + 9);

                    auto g_x_xyyyy = contrBuffer.data(goff + 63 * j + 10);

                    auto g_x_xyyyz = contrBuffer.data(goff + 63 * j + 11);

                    auto g_x_xyyzz = contrBuffer.data(goff + 63 * j + 12);

                    auto g_x_xyzzz = contrBuffer.data(goff + 63 * j + 13);

                    auto g_x_xzzzz = contrBuffer.data(goff + 63 * j + 14);

                    auto g_x_yyyyy = contrBuffer.data(goff + 63 * j + 15);

                    auto g_x_yyyyz = contrBuffer.data(goff + 63 * j + 16);

                    auto g_x_yyyzz = contrBuffer.data(goff + 63 * j + 17);

                    auto g_x_yyzzz = contrBuffer.data(goff + 63 * j + 18);

                    auto g_x_yzzzz = contrBuffer.data(goff + 63 * j + 19);

                    auto g_x_zzzzz = contrBuffer.data(goff + 63 * j + 20);

                    auto g_y_xxxxx = contrBuffer.data(goff + 63 * j + 21);

                    auto g_y_xxxxy = contrBuffer.data(goff + 63 * j + 22);

                    auto g_y_xxxxz = contrBuffer.data(goff + 63 * j + 23);

                    auto g_y_xxxyy = contrBuffer.data(goff + 63 * j + 24);

                    auto g_y_xxxyz = contrBuffer.data(goff + 63 * j + 25);

                    auto g_y_xxxzz = contrBuffer.data(goff + 63 * j + 26);

                    auto g_y_xxyyy = contrBuffer.data(goff + 63 * j + 27);

                    auto g_y_xxyyz = contrBuffer.data(goff + 63 * j + 28);

                    auto g_y_xxyzz = contrBuffer.data(goff + 63 * j + 29);

                    auto g_y_xxzzz = contrBuffer.data(goff + 63 * j + 30);

                    auto g_y_xyyyy = contrBuffer.data(goff + 63 * j + 31);

                    auto g_y_xyyyz = contrBuffer.data(goff + 63 * j + 32);

                    auto g_y_xyyzz = contrBuffer.data(goff + 63 * j + 33);

                    auto g_y_xyzzz = contrBuffer.data(goff + 63 * j + 34);

                    auto g_y_xzzzz = contrBuffer.data(goff + 63 * j + 35);

                    auto g_y_yyyyy = contrBuffer.data(goff + 63 * j + 36);

                    auto g_y_yyyyz = contrBuffer.data(goff + 63 * j + 37);

                    auto g_y_yyyzz = contrBuffer.data(goff + 63 * j + 38);

                    auto g_y_yyzzz = contrBuffer.data(goff + 63 * j + 39);

                    auto g_y_yzzzz = contrBuffer.data(goff + 63 * j + 40);

                    auto g_y_zzzzz = contrBuffer.data(goff + 63 * j + 41);

                    auto g_z_xxxxx = contrBuffer.data(goff + 63 * j + 42);

                    auto g_z_xxxxy = contrBuffer.data(goff + 63 * j + 43);

                    auto g_z_xxxxz = contrBuffer.data(goff + 63 * j + 44);

                    auto g_z_xxxyy = contrBuffer.data(goff + 63 * j + 45);

                    auto g_z_xxxyz = contrBuffer.data(goff + 63 * j + 46);

                    auto g_z_xxxzz = contrBuffer.data(goff + 63 * j + 47);

                    auto g_z_xxyyy = contrBuffer.data(goff + 63 * j + 48);

                    auto g_z_xxyyz = contrBuffer.data(goff + 63 * j + 49);

                    auto g_z_xxyzz = contrBuffer.data(goff + 63 * j + 50);

                    auto g_z_xxzzz = contrBuffer.data(goff + 63 * j + 51);

                    auto g_z_xyyyy = contrBuffer.data(goff + 63 * j + 52);

                    auto g_z_xyyyz = contrBuffer.data(goff + 63 * j + 53);

                    auto g_z_xyyzz = contrBuffer.data(goff + 63 * j + 54);

                    auto g_z_xyzzz = contrBuffer.data(goff + 63 * j + 55);

                    auto g_z_xzzzz = contrBuffer.data(goff + 63 * j + 56);

                    auto g_z_yyyyy = contrBuffer.data(goff + 63 * j + 57);

                    auto g_z_yyyyz = contrBuffer.data(goff + 63 * j + 58);

                    auto g_z_yyyzz = contrBuffer.data(goff + 63 * j + 59);

                    auto g_z_yyzzz = contrBuffer.data(goff + 63 * j + 60);

                    auto g_z_yzzzz = contrBuffer.data(goff + 63 * j + 61);

                    auto g_z_zzzzz = contrBuffer.data(goff + 63 * j + 62);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxxx, g2_0_xxxxy,\
                                             g2_0_xxxxz, g2_0_xxxyy, g2_0_xxxyz,\
                                             g2_0_xxxzz, g2_0_xxyyy, g2_0_xxyyz,\
                                             g2_0_xxyzz, g2_0_xxzzz, g2_0_xyyyy,\
                                             g2_0_xyyyz, g2_0_xyyzz, g2_0_xyzzz,\
                                             g2_0_xzzzz, g2_0_yyyyy, g2_0_yyyyz,\
                                             g2_0_yyyzz, g2_0_yyzzz, g2_0_yzzzz,\
                                             g2_0_zzzzz, g1_0_xxxxxx, g1_0_xxxxxy,\
                                             g1_0_xxxxxz, g1_0_xxxxyy, g1_0_xxxxyz,\
                                             g1_0_xxxxzz, g1_0_xxxyyy, g1_0_xxxyyz,\
                                             g1_0_xxxyzz, g1_0_xxxzzz, g1_0_xxyyyy,\
                                             g1_0_xxyyyz, g1_0_xxyyzz, g1_0_xxyzzz,\
                                             g1_0_xxzzzz, g1_0_xyyyyy, g1_0_xyyyyz,\
                                             g1_0_xyyyzz, g1_0_xyyzzz, g1_0_xyzzzz,\
                                             g1_0_xzzzzz, g1_0_yyyyyy, g1_0_yyyyyz,\
                                             g1_0_yyyyzz, g1_0_yyyzzz, g1_0_yyzzzz,\
                                             g1_0_yzzzzz, g1_0_zzzzzz, g_x_xxxxx,\
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
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xxxxx[k] = g1_0_xxxxxx[k] - fr * g2_0_xxxxx[k];

                        g_x_xxxxy[k] = g1_0_xxxxxy[k] - fr * g2_0_xxxxy[k];

                        g_x_xxxxz[k] = g1_0_xxxxxz[k] - fr * g2_0_xxxxz[k];

                        g_x_xxxyy[k] = g1_0_xxxxyy[k] - fr * g2_0_xxxyy[k];

                        g_x_xxxyz[k] = g1_0_xxxxyz[k] - fr * g2_0_xxxyz[k];

                        g_x_xxxzz[k] = g1_0_xxxxzz[k] - fr * g2_0_xxxzz[k];

                        g_x_xxyyy[k] = g1_0_xxxyyy[k] - fr * g2_0_xxyyy[k];

                        g_x_xxyyz[k] = g1_0_xxxyyz[k] - fr * g2_0_xxyyz[k];

                        g_x_xxyzz[k] = g1_0_xxxyzz[k] - fr * g2_0_xxyzz[k];

                        g_x_xxzzz[k] = g1_0_xxxzzz[k] - fr * g2_0_xxzzz[k];

                        g_x_xyyyy[k] = g1_0_xxyyyy[k] - fr * g2_0_xyyyy[k];

                        g_x_xyyyz[k] = g1_0_xxyyyz[k] - fr * g2_0_xyyyz[k];

                        g_x_xyyzz[k] = g1_0_xxyyzz[k] - fr * g2_0_xyyzz[k];

                        g_x_xyzzz[k] = g1_0_xxyzzz[k] - fr * g2_0_xyzzz[k];

                        g_x_xzzzz[k] = g1_0_xxzzzz[k] - fr * g2_0_xzzzz[k];

                        g_x_yyyyy[k] = g1_0_xyyyyy[k] - fr * g2_0_yyyyy[k];

                        g_x_yyyyz[k] = g1_0_xyyyyz[k] - fr * g2_0_yyyyz[k];

                        g_x_yyyzz[k] = g1_0_xyyyzz[k] - fr * g2_0_yyyzz[k];

                        g_x_yyzzz[k] = g1_0_xyyzzz[k] - fr * g2_0_yyzzz[k];

                        g_x_yzzzz[k] = g1_0_xyzzzz[k] - fr * g2_0_yzzzz[k];

                        g_x_zzzzz[k] = g1_0_xzzzzz[k] - fr * g2_0_zzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xxxxx[k] = g1_0_xxxxxy[k] - fr * g2_0_xxxxx[k];

                        g_y_xxxxy[k] = g1_0_xxxxyy[k] - fr * g2_0_xxxxy[k];

                        g_y_xxxxz[k] = g1_0_xxxxyz[k] - fr * g2_0_xxxxz[k];

                        g_y_xxxyy[k] = g1_0_xxxyyy[k] - fr * g2_0_xxxyy[k];

                        g_y_xxxyz[k] = g1_0_xxxyyz[k] - fr * g2_0_xxxyz[k];

                        g_y_xxxzz[k] = g1_0_xxxyzz[k] - fr * g2_0_xxxzz[k];

                        g_y_xxyyy[k] = g1_0_xxyyyy[k] - fr * g2_0_xxyyy[k];

                        g_y_xxyyz[k] = g1_0_xxyyyz[k] - fr * g2_0_xxyyz[k];

                        g_y_xxyzz[k] = g1_0_xxyyzz[k] - fr * g2_0_xxyzz[k];

                        g_y_xxzzz[k] = g1_0_xxyzzz[k] - fr * g2_0_xxzzz[k];

                        g_y_xyyyy[k] = g1_0_xyyyyy[k] - fr * g2_0_xyyyy[k];

                        g_y_xyyyz[k] = g1_0_xyyyyz[k] - fr * g2_0_xyyyz[k];

                        g_y_xyyzz[k] = g1_0_xyyyzz[k] - fr * g2_0_xyyzz[k];

                        g_y_xyzzz[k] = g1_0_xyyzzz[k] - fr * g2_0_xyzzz[k];

                        g_y_xzzzz[k] = g1_0_xyzzzz[k] - fr * g2_0_xzzzz[k];

                        g_y_yyyyy[k] = g1_0_yyyyyy[k] - fr * g2_0_yyyyy[k];

                        g_y_yyyyz[k] = g1_0_yyyyyz[k] - fr * g2_0_yyyyz[k];

                        g_y_yyyzz[k] = g1_0_yyyyzz[k] - fr * g2_0_yyyzz[k];

                        g_y_yyzzz[k] = g1_0_yyyzzz[k] - fr * g2_0_yyzzz[k];

                        g_y_yzzzz[k] = g1_0_yyzzzz[k] - fr * g2_0_yzzzz[k];

                        g_y_zzzzz[k] = g1_0_yzzzzz[k] - fr * g2_0_zzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xxxxx[k] = g1_0_xxxxxz[k] - fr * g2_0_xxxxx[k];

                        g_z_xxxxy[k] = g1_0_xxxxyz[k] - fr * g2_0_xxxxy[k];

                        g_z_xxxxz[k] = g1_0_xxxxzz[k] - fr * g2_0_xxxxz[k];

                        g_z_xxxyy[k] = g1_0_xxxyyz[k] - fr * g2_0_xxxyy[k];

                        g_z_xxxyz[k] = g1_0_xxxyzz[k] - fr * g2_0_xxxyz[k];

                        g_z_xxxzz[k] = g1_0_xxxzzz[k] - fr * g2_0_xxxzz[k];

                        g_z_xxyyy[k] = g1_0_xxyyyz[k] - fr * g2_0_xxyyy[k];

                        g_z_xxyyz[k] = g1_0_xxyyzz[k] - fr * g2_0_xxyyz[k];

                        g_z_xxyzz[k] = g1_0_xxyzzz[k] - fr * g2_0_xxyzz[k];

                        g_z_xxzzz[k] = g1_0_xxzzzz[k] - fr * g2_0_xxzzz[k];

                        g_z_xyyyy[k] = g1_0_xyyyyz[k] - fr * g2_0_xyyyy[k];

                        g_z_xyyyz[k] = g1_0_xyyyzz[k] - fr * g2_0_xyyyz[k];

                        g_z_xyyzz[k] = g1_0_xyyzzz[k] - fr * g2_0_xyyzz[k];

                        g_z_xyzzz[k] = g1_0_xyzzzz[k] - fr * g2_0_xyzzz[k];

                        g_z_xzzzz[k] = g1_0_xzzzzz[k] - fr * g2_0_xzzzz[k];

                        g_z_yyyyy[k] = g1_0_yyyyyz[k] - fr * g2_0_yyyyy[k];

                        g_z_yyyyz[k] = g1_0_yyyyzz[k] - fr * g2_0_yyyyz[k];

                        g_z_yyyzz[k] = g1_0_yyyzzz[k] - fr * g2_0_yyyzz[k];

                        g_z_yyzzz[k] = g1_0_yyzzzz[k] - fr * g2_0_yyzzz[k];

                        g_z_yzzzz[k] = g1_0_yzzzzz[k] - fr * g2_0_yzzzz[k];

                        g_z_zzzzz[k] = g1_0_zzzzzz[k] - fr * g2_0_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXPI(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 6))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|16)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 6});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 7});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 6});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SI)^(m) integrals

                    auto g2_0_xxxxxx = contrBuffer.data(g2off + 28 * j);

                    auto g2_0_xxxxxy = contrBuffer.data(g2off + 28 * j + 1);

                    auto g2_0_xxxxxz = contrBuffer.data(g2off + 28 * j + 2);

                    auto g2_0_xxxxyy = contrBuffer.data(g2off + 28 * j + 3);

                    auto g2_0_xxxxyz = contrBuffer.data(g2off + 28 * j + 4);

                    auto g2_0_xxxxzz = contrBuffer.data(g2off + 28 * j + 5);

                    auto g2_0_xxxyyy = contrBuffer.data(g2off + 28 * j + 6);

                    auto g2_0_xxxyyz = contrBuffer.data(g2off + 28 * j + 7);

                    auto g2_0_xxxyzz = contrBuffer.data(g2off + 28 * j + 8);

                    auto g2_0_xxxzzz = contrBuffer.data(g2off + 28 * j + 9);

                    auto g2_0_xxyyyy = contrBuffer.data(g2off + 28 * j + 10);

                    auto g2_0_xxyyyz = contrBuffer.data(g2off + 28 * j + 11);

                    auto g2_0_xxyyzz = contrBuffer.data(g2off + 28 * j + 12);

                    auto g2_0_xxyzzz = contrBuffer.data(g2off + 28 * j + 13);

                    auto g2_0_xxzzzz = contrBuffer.data(g2off + 28 * j + 14);

                    auto g2_0_xyyyyy = contrBuffer.data(g2off + 28 * j + 15);

                    auto g2_0_xyyyyz = contrBuffer.data(g2off + 28 * j + 16);

                    auto g2_0_xyyyzz = contrBuffer.data(g2off + 28 * j + 17);

                    auto g2_0_xyyzzz = contrBuffer.data(g2off + 28 * j + 18);

                    auto g2_0_xyzzzz = contrBuffer.data(g2off + 28 * j + 19);

                    auto g2_0_xzzzzz = contrBuffer.data(g2off + 28 * j + 20);

                    auto g2_0_yyyyyy = contrBuffer.data(g2off + 28 * j + 21);

                    auto g2_0_yyyyyz = contrBuffer.data(g2off + 28 * j + 22);

                    auto g2_0_yyyyzz = contrBuffer.data(g2off + 28 * j + 23);

                    auto g2_0_yyyzzz = contrBuffer.data(g2off + 28 * j + 24);

                    auto g2_0_yyzzzz = contrBuffer.data(g2off + 28 * j + 25);

                    auto g2_0_yzzzzz = contrBuffer.data(g2off + 28 * j + 26);

                    auto g2_0_zzzzzz = contrBuffer.data(g2off + 28 * j + 27);

                    // set up pointers to (SX|g(r,r')|SK)^(m) integrals

                    auto g1_0_xxxxxxx = contrBuffer.data(g1off + 36 * j);

                    auto g1_0_xxxxxxy = contrBuffer.data(g1off + 36 * j + 1);

                    auto g1_0_xxxxxxz = contrBuffer.data(g1off + 36 * j + 2);

                    auto g1_0_xxxxxyy = contrBuffer.data(g1off + 36 * j + 3);

                    auto g1_0_xxxxxyz = contrBuffer.data(g1off + 36 * j + 4);

                    auto g1_0_xxxxxzz = contrBuffer.data(g1off + 36 * j + 5);

                    auto g1_0_xxxxyyy = contrBuffer.data(g1off + 36 * j + 6);

                    auto g1_0_xxxxyyz = contrBuffer.data(g1off + 36 * j + 7);

                    auto g1_0_xxxxyzz = contrBuffer.data(g1off + 36 * j + 8);

                    auto g1_0_xxxxzzz = contrBuffer.data(g1off + 36 * j + 9);

                    auto g1_0_xxxyyyy = contrBuffer.data(g1off + 36 * j + 10);

                    auto g1_0_xxxyyyz = contrBuffer.data(g1off + 36 * j + 11);

                    auto g1_0_xxxyyzz = contrBuffer.data(g1off + 36 * j + 12);

                    auto g1_0_xxxyzzz = contrBuffer.data(g1off + 36 * j + 13);

                    auto g1_0_xxxzzzz = contrBuffer.data(g1off + 36 * j + 14);

                    auto g1_0_xxyyyyy = contrBuffer.data(g1off + 36 * j + 15);

                    auto g1_0_xxyyyyz = contrBuffer.data(g1off + 36 * j + 16);

                    auto g1_0_xxyyyzz = contrBuffer.data(g1off + 36 * j + 17);

                    auto g1_0_xxyyzzz = contrBuffer.data(g1off + 36 * j + 18);

                    auto g1_0_xxyzzzz = contrBuffer.data(g1off + 36 * j + 19);

                    auto g1_0_xxzzzzz = contrBuffer.data(g1off + 36 * j + 20);

                    auto g1_0_xyyyyyy = contrBuffer.data(g1off + 36 * j + 21);

                    auto g1_0_xyyyyyz = contrBuffer.data(g1off + 36 * j + 22);

                    auto g1_0_xyyyyzz = contrBuffer.data(g1off + 36 * j + 23);

                    auto g1_0_xyyyzzz = contrBuffer.data(g1off + 36 * j + 24);

                    auto g1_0_xyyzzzz = contrBuffer.data(g1off + 36 * j + 25);

                    auto g1_0_xyzzzzz = contrBuffer.data(g1off + 36 * j + 26);

                    auto g1_0_xzzzzzz = contrBuffer.data(g1off + 36 * j + 27);

                    auto g1_0_yyyyyyy = contrBuffer.data(g1off + 36 * j + 28);

                    auto g1_0_yyyyyyz = contrBuffer.data(g1off + 36 * j + 29);

                    auto g1_0_yyyyyzz = contrBuffer.data(g1off + 36 * j + 30);

                    auto g1_0_yyyyzzz = contrBuffer.data(g1off + 36 * j + 31);

                    auto g1_0_yyyzzzz = contrBuffer.data(g1off + 36 * j + 32);

                    auto g1_0_yyzzzzz = contrBuffer.data(g1off + 36 * j + 33);

                    auto g1_0_yzzzzzz = contrBuffer.data(g1off + 36 * j + 34);

                    auto g1_0_zzzzzzz = contrBuffer.data(g1off + 36 * j + 35);

                    // set up pointers to (SX|g(r,r')|PI)^(m) integrals

                    auto g_x_xxxxxx = contrBuffer.data(goff + 84 * j);

                    auto g_x_xxxxxy = contrBuffer.data(goff + 84 * j + 1);

                    auto g_x_xxxxxz = contrBuffer.data(goff + 84 * j + 2);

                    auto g_x_xxxxyy = contrBuffer.data(goff + 84 * j + 3);

                    auto g_x_xxxxyz = contrBuffer.data(goff + 84 * j + 4);

                    auto g_x_xxxxzz = contrBuffer.data(goff + 84 * j + 5);

                    auto g_x_xxxyyy = contrBuffer.data(goff + 84 * j + 6);

                    auto g_x_xxxyyz = contrBuffer.data(goff + 84 * j + 7);

                    auto g_x_xxxyzz = contrBuffer.data(goff + 84 * j + 8);

                    auto g_x_xxxzzz = contrBuffer.data(goff + 84 * j + 9);

                    auto g_x_xxyyyy = contrBuffer.data(goff + 84 * j + 10);

                    auto g_x_xxyyyz = contrBuffer.data(goff + 84 * j + 11);

                    auto g_x_xxyyzz = contrBuffer.data(goff + 84 * j + 12);

                    auto g_x_xxyzzz = contrBuffer.data(goff + 84 * j + 13);

                    auto g_x_xxzzzz = contrBuffer.data(goff + 84 * j + 14);

                    auto g_x_xyyyyy = contrBuffer.data(goff + 84 * j + 15);

                    auto g_x_xyyyyz = contrBuffer.data(goff + 84 * j + 16);

                    auto g_x_xyyyzz = contrBuffer.data(goff + 84 * j + 17);

                    auto g_x_xyyzzz = contrBuffer.data(goff + 84 * j + 18);

                    auto g_x_xyzzzz = contrBuffer.data(goff + 84 * j + 19);

                    auto g_x_xzzzzz = contrBuffer.data(goff + 84 * j + 20);

                    auto g_x_yyyyyy = contrBuffer.data(goff + 84 * j + 21);

                    auto g_x_yyyyyz = contrBuffer.data(goff + 84 * j + 22);

                    auto g_x_yyyyzz = contrBuffer.data(goff + 84 * j + 23);

                    auto g_x_yyyzzz = contrBuffer.data(goff + 84 * j + 24);

                    auto g_x_yyzzzz = contrBuffer.data(goff + 84 * j + 25);

                    auto g_x_yzzzzz = contrBuffer.data(goff + 84 * j + 26);

                    auto g_x_zzzzzz = contrBuffer.data(goff + 84 * j + 27);

                    auto g_y_xxxxxx = contrBuffer.data(goff + 84 * j + 28);

                    auto g_y_xxxxxy = contrBuffer.data(goff + 84 * j + 29);

                    auto g_y_xxxxxz = contrBuffer.data(goff + 84 * j + 30);

                    auto g_y_xxxxyy = contrBuffer.data(goff + 84 * j + 31);

                    auto g_y_xxxxyz = contrBuffer.data(goff + 84 * j + 32);

                    auto g_y_xxxxzz = contrBuffer.data(goff + 84 * j + 33);

                    auto g_y_xxxyyy = contrBuffer.data(goff + 84 * j + 34);

                    auto g_y_xxxyyz = contrBuffer.data(goff + 84 * j + 35);

                    auto g_y_xxxyzz = contrBuffer.data(goff + 84 * j + 36);

                    auto g_y_xxxzzz = contrBuffer.data(goff + 84 * j + 37);

                    auto g_y_xxyyyy = contrBuffer.data(goff + 84 * j + 38);

                    auto g_y_xxyyyz = contrBuffer.data(goff + 84 * j + 39);

                    auto g_y_xxyyzz = contrBuffer.data(goff + 84 * j + 40);

                    auto g_y_xxyzzz = contrBuffer.data(goff + 84 * j + 41);

                    auto g_y_xxzzzz = contrBuffer.data(goff + 84 * j + 42);

                    auto g_y_xyyyyy = contrBuffer.data(goff + 84 * j + 43);

                    auto g_y_xyyyyz = contrBuffer.data(goff + 84 * j + 44);

                    auto g_y_xyyyzz = contrBuffer.data(goff + 84 * j + 45);

                    auto g_y_xyyzzz = contrBuffer.data(goff + 84 * j + 46);

                    auto g_y_xyzzzz = contrBuffer.data(goff + 84 * j + 47);

                    auto g_y_xzzzzz = contrBuffer.data(goff + 84 * j + 48);

                    auto g_y_yyyyyy = contrBuffer.data(goff + 84 * j + 49);

                    auto g_y_yyyyyz = contrBuffer.data(goff + 84 * j + 50);

                    auto g_y_yyyyzz = contrBuffer.data(goff + 84 * j + 51);

                    auto g_y_yyyzzz = contrBuffer.data(goff + 84 * j + 52);

                    auto g_y_yyzzzz = contrBuffer.data(goff + 84 * j + 53);

                    auto g_y_yzzzzz = contrBuffer.data(goff + 84 * j + 54);

                    auto g_y_zzzzzz = contrBuffer.data(goff + 84 * j + 55);

                    auto g_z_xxxxxx = contrBuffer.data(goff + 84 * j + 56);

                    auto g_z_xxxxxy = contrBuffer.data(goff + 84 * j + 57);

                    auto g_z_xxxxxz = contrBuffer.data(goff + 84 * j + 58);

                    auto g_z_xxxxyy = contrBuffer.data(goff + 84 * j + 59);

                    auto g_z_xxxxyz = contrBuffer.data(goff + 84 * j + 60);

                    auto g_z_xxxxzz = contrBuffer.data(goff + 84 * j + 61);

                    auto g_z_xxxyyy = contrBuffer.data(goff + 84 * j + 62);

                    auto g_z_xxxyyz = contrBuffer.data(goff + 84 * j + 63);

                    auto g_z_xxxyzz = contrBuffer.data(goff + 84 * j + 64);

                    auto g_z_xxxzzz = contrBuffer.data(goff + 84 * j + 65);

                    auto g_z_xxyyyy = contrBuffer.data(goff + 84 * j + 66);

                    auto g_z_xxyyyz = contrBuffer.data(goff + 84 * j + 67);

                    auto g_z_xxyyzz = contrBuffer.data(goff + 84 * j + 68);

                    auto g_z_xxyzzz = contrBuffer.data(goff + 84 * j + 69);

                    auto g_z_xxzzzz = contrBuffer.data(goff + 84 * j + 70);

                    auto g_z_xyyyyy = contrBuffer.data(goff + 84 * j + 71);

                    auto g_z_xyyyyz = contrBuffer.data(goff + 84 * j + 72);

                    auto g_z_xyyyzz = contrBuffer.data(goff + 84 * j + 73);

                    auto g_z_xyyzzz = contrBuffer.data(goff + 84 * j + 74);

                    auto g_z_xyzzzz = contrBuffer.data(goff + 84 * j + 75);

                    auto g_z_xzzzzz = contrBuffer.data(goff + 84 * j + 76);

                    auto g_z_yyyyyy = contrBuffer.data(goff + 84 * j + 77);

                    auto g_z_yyyyyz = contrBuffer.data(goff + 84 * j + 78);

                    auto g_z_yyyyzz = contrBuffer.data(goff + 84 * j + 79);

                    auto g_z_yyyzzz = contrBuffer.data(goff + 84 * j + 80);

                    auto g_z_yyzzzz = contrBuffer.data(goff + 84 * j + 81);

                    auto g_z_yzzzzz = contrBuffer.data(goff + 84 * j + 82);

                    auto g_z_zzzzzz = contrBuffer.data(goff + 84 * j + 83);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_0_xxxxxx, g2_0_xxxxxy,\
                                             g2_0_xxxxxz, g2_0_xxxxyy, g2_0_xxxxyz,\
                                             g2_0_xxxxzz, g2_0_xxxyyy, g2_0_xxxyyz,\
                                             g2_0_xxxyzz, g2_0_xxxzzz, g2_0_xxyyyy,\
                                             g2_0_xxyyyz, g2_0_xxyyzz, g2_0_xxyzzz,\
                                             g2_0_xxzzzz, g2_0_xyyyyy, g2_0_xyyyyz,\
                                             g2_0_xyyyzz, g2_0_xyyzzz, g2_0_xyzzzz,\
                                             g2_0_xzzzzz, g2_0_yyyyyy, g2_0_yyyyyz,\
                                             g2_0_yyyyzz, g2_0_yyyzzz, g2_0_yyzzzz,\
                                             g2_0_yzzzzz, g2_0_zzzzzz, g1_0_xxxxxxx,\
                                             g1_0_xxxxxxy, g1_0_xxxxxxz, g1_0_xxxxxyy,\
                                             g1_0_xxxxxyz, g1_0_xxxxxzz, g1_0_xxxxyyy,\
                                             g1_0_xxxxyyz, g1_0_xxxxyzz, g1_0_xxxxzzz,\
                                             g1_0_xxxyyyy, g1_0_xxxyyyz, g1_0_xxxyyzz,\
                                             g1_0_xxxyzzz, g1_0_xxxzzzz, g1_0_xxyyyyy,\
                                             g1_0_xxyyyyz, g1_0_xxyyyzz, g1_0_xxyyzzz,\
                                             g1_0_xxyzzzz, g1_0_xxzzzzz, g1_0_xyyyyyy,\
                                             g1_0_xyyyyyz, g1_0_xyyyyzz, g1_0_xyyyzzz,\
                                             g1_0_xyyzzzz, g1_0_xyzzzzz, g1_0_xzzzzzz,\
                                             g1_0_yyyyyyy, g1_0_yyyyyyz, g1_0_yyyyyzz,\
                                             g1_0_yyyyzzz, g1_0_yyyzzzz, g1_0_yyzzzzz,\
                                             g1_0_yzzzzzz, g1_0_zzzzzzz, g_x_xxxxxx,\
                                             g_x_xxxxxy, g_x_xxxxxz, g_x_xxxxyy,\
                                             g_x_xxxxyz, g_x_xxxxzz, g_x_xxxyyy,\
                                             g_x_xxxyyz, g_x_xxxyzz, g_x_xxxzzz,\
                                             g_x_xxyyyy, g_x_xxyyyz, g_x_xxyyzz,\
                                             g_x_xxyzzz, g_x_xxzzzz, g_x_xyyyyy,\
                                             g_x_xyyyyz, g_x_xyyyzz, g_x_xyyzzz,\
                                             g_x_xyzzzz, g_x_xzzzzz, g_x_yyyyyy,\
                                             g_x_yyyyyz, g_x_yyyyzz, g_x_yyyzzz,\
                                             g_x_yyzzzz, g_x_yzzzzz, g_x_zzzzzz,\
                                             g_y_xxxxxx, g_y_xxxxxy, g_y_xxxxxz,\
                                             g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxzz,\
                                             g_y_xxxyyy, g_y_xxxyyz, g_y_xxxyzz,\
                                             g_y_xxxzzz, g_y_xxyyyy, g_y_xxyyyz,\
                                             g_y_xxyyzz, g_y_xxyzzz, g_y_xxzzzz,\
                                             g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyzz,\
                                             g_y_xyyzzz, g_y_xyzzzz, g_y_xzzzzz,\
                                             g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyzz,\
                                             g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz,\
                                             g_y_zzzzzz, g_z_xxxxxx, g_z_xxxxxy,\
                                             g_z_xxxxxz, g_z_xxxxyy, g_z_xxxxyz,\
                                             g_z_xxxxzz, g_z_xxxyyy, g_z_xxxyyz,\
                                             g_z_xxxyzz, g_z_xxxzzz, g_z_xxyyyy,\
                                             g_z_xxyyyz, g_z_xxyyzz, g_z_xxyzzz,\
                                             g_z_xxzzzz, g_z_xyyyyy, g_z_xyyyyz,\
                                             g_z_xyyyzz, g_z_xyyzzz, g_z_xyzzzz,\
                                             g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyz,\
                                             g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz,\
                                             g_z_yzzzzz, g_z_zzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xxxxxx[k] = g1_0_xxxxxxx[k] - fr * g2_0_xxxxxx[k];

                        g_x_xxxxxy[k] = g1_0_xxxxxxy[k] - fr * g2_0_xxxxxy[k];

                        g_x_xxxxxz[k] = g1_0_xxxxxxz[k] - fr * g2_0_xxxxxz[k];

                        g_x_xxxxyy[k] = g1_0_xxxxxyy[k] - fr * g2_0_xxxxyy[k];

                        g_x_xxxxyz[k] = g1_0_xxxxxyz[k] - fr * g2_0_xxxxyz[k];

                        g_x_xxxxzz[k] = g1_0_xxxxxzz[k] - fr * g2_0_xxxxzz[k];

                        g_x_xxxyyy[k] = g1_0_xxxxyyy[k] - fr * g2_0_xxxyyy[k];

                        g_x_xxxyyz[k] = g1_0_xxxxyyz[k] - fr * g2_0_xxxyyz[k];

                        g_x_xxxyzz[k] = g1_0_xxxxyzz[k] - fr * g2_0_xxxyzz[k];

                        g_x_xxxzzz[k] = g1_0_xxxxzzz[k] - fr * g2_0_xxxzzz[k];

                        g_x_xxyyyy[k] = g1_0_xxxyyyy[k] - fr * g2_0_xxyyyy[k];

                        g_x_xxyyyz[k] = g1_0_xxxyyyz[k] - fr * g2_0_xxyyyz[k];

                        g_x_xxyyzz[k] = g1_0_xxxyyzz[k] - fr * g2_0_xxyyzz[k];

                        g_x_xxyzzz[k] = g1_0_xxxyzzz[k] - fr * g2_0_xxyzzz[k];

                        g_x_xxzzzz[k] = g1_0_xxxzzzz[k] - fr * g2_0_xxzzzz[k];

                        g_x_xyyyyy[k] = g1_0_xxyyyyy[k] - fr * g2_0_xyyyyy[k];

                        g_x_xyyyyz[k] = g1_0_xxyyyyz[k] - fr * g2_0_xyyyyz[k];

                        g_x_xyyyzz[k] = g1_0_xxyyyzz[k] - fr * g2_0_xyyyzz[k];

                        g_x_xyyzzz[k] = g1_0_xxyyzzz[k] - fr * g2_0_xyyzzz[k];

                        g_x_xyzzzz[k] = g1_0_xxyzzzz[k] - fr * g2_0_xyzzzz[k];

                        g_x_xzzzzz[k] = g1_0_xxzzzzz[k] - fr * g2_0_xzzzzz[k];

                        g_x_yyyyyy[k] = g1_0_xyyyyyy[k] - fr * g2_0_yyyyyy[k];

                        g_x_yyyyyz[k] = g1_0_xyyyyyz[k] - fr * g2_0_yyyyyz[k];

                        g_x_yyyyzz[k] = g1_0_xyyyyzz[k] - fr * g2_0_yyyyzz[k];

                        g_x_yyyzzz[k] = g1_0_xyyyzzz[k] - fr * g2_0_yyyzzz[k];

                        g_x_yyzzzz[k] = g1_0_xyyzzzz[k] - fr * g2_0_yyzzzz[k];

                        g_x_yzzzzz[k] = g1_0_xyzzzzz[k] - fr * g2_0_yzzzzz[k];

                        g_x_zzzzzz[k] = g1_0_xzzzzzz[k] - fr * g2_0_zzzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xxxxxx[k] = g1_0_xxxxxxy[k] - fr * g2_0_xxxxxx[k];

                        g_y_xxxxxy[k] = g1_0_xxxxxyy[k] - fr * g2_0_xxxxxy[k];

                        g_y_xxxxxz[k] = g1_0_xxxxxyz[k] - fr * g2_0_xxxxxz[k];

                        g_y_xxxxyy[k] = g1_0_xxxxyyy[k] - fr * g2_0_xxxxyy[k];

                        g_y_xxxxyz[k] = g1_0_xxxxyyz[k] - fr * g2_0_xxxxyz[k];

                        g_y_xxxxzz[k] = g1_0_xxxxyzz[k] - fr * g2_0_xxxxzz[k];

                        g_y_xxxyyy[k] = g1_0_xxxyyyy[k] - fr * g2_0_xxxyyy[k];

                        g_y_xxxyyz[k] = g1_0_xxxyyyz[k] - fr * g2_0_xxxyyz[k];

                        g_y_xxxyzz[k] = g1_0_xxxyyzz[k] - fr * g2_0_xxxyzz[k];

                        g_y_xxxzzz[k] = g1_0_xxxyzzz[k] - fr * g2_0_xxxzzz[k];

                        g_y_xxyyyy[k] = g1_0_xxyyyyy[k] - fr * g2_0_xxyyyy[k];

                        g_y_xxyyyz[k] = g1_0_xxyyyyz[k] - fr * g2_0_xxyyyz[k];

                        g_y_xxyyzz[k] = g1_0_xxyyyzz[k] - fr * g2_0_xxyyzz[k];

                        g_y_xxyzzz[k] = g1_0_xxyyzzz[k] - fr * g2_0_xxyzzz[k];

                        g_y_xxzzzz[k] = g1_0_xxyzzzz[k] - fr * g2_0_xxzzzz[k];

                        g_y_xyyyyy[k] = g1_0_xyyyyyy[k] - fr * g2_0_xyyyyy[k];

                        g_y_xyyyyz[k] = g1_0_xyyyyyz[k] - fr * g2_0_xyyyyz[k];

                        g_y_xyyyzz[k] = g1_0_xyyyyzz[k] - fr * g2_0_xyyyzz[k];

                        g_y_xyyzzz[k] = g1_0_xyyyzzz[k] - fr * g2_0_xyyzzz[k];

                        g_y_xyzzzz[k] = g1_0_xyyzzzz[k] - fr * g2_0_xyzzzz[k];

                        g_y_xzzzzz[k] = g1_0_xyzzzzz[k] - fr * g2_0_xzzzzz[k];

                        g_y_yyyyyy[k] = g1_0_yyyyyyy[k] - fr * g2_0_yyyyyy[k];

                        g_y_yyyyyz[k] = g1_0_yyyyyyz[k] - fr * g2_0_yyyyyz[k];

                        g_y_yyyyzz[k] = g1_0_yyyyyzz[k] - fr * g2_0_yyyyzz[k];

                        g_y_yyyzzz[k] = g1_0_yyyyzzz[k] - fr * g2_0_yyyzzz[k];

                        g_y_yyzzzz[k] = g1_0_yyyzzzz[k] - fr * g2_0_yyzzzz[k];

                        g_y_yzzzzz[k] = g1_0_yyzzzzz[k] - fr * g2_0_yzzzzz[k];

                        g_y_zzzzzz[k] = g1_0_yzzzzzz[k] - fr * g2_0_zzzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xxxxxx[k] = g1_0_xxxxxxz[k] - fr * g2_0_xxxxxx[k];

                        g_z_xxxxxy[k] = g1_0_xxxxxyz[k] - fr * g2_0_xxxxxy[k];

                        g_z_xxxxxz[k] = g1_0_xxxxxzz[k] - fr * g2_0_xxxxxz[k];

                        g_z_xxxxyy[k] = g1_0_xxxxyyz[k] - fr * g2_0_xxxxyy[k];

                        g_z_xxxxyz[k] = g1_0_xxxxyzz[k] - fr * g2_0_xxxxyz[k];

                        g_z_xxxxzz[k] = g1_0_xxxxzzz[k] - fr * g2_0_xxxxzz[k];

                        g_z_xxxyyy[k] = g1_0_xxxyyyz[k] - fr * g2_0_xxxyyy[k];

                        g_z_xxxyyz[k] = g1_0_xxxyyzz[k] - fr * g2_0_xxxyyz[k];

                        g_z_xxxyzz[k] = g1_0_xxxyzzz[k] - fr * g2_0_xxxyzz[k];

                        g_z_xxxzzz[k] = g1_0_xxxzzzz[k] - fr * g2_0_xxxzzz[k];

                        g_z_xxyyyy[k] = g1_0_xxyyyyz[k] - fr * g2_0_xxyyyy[k];

                        g_z_xxyyyz[k] = g1_0_xxyyyzz[k] - fr * g2_0_xxyyyz[k];

                        g_z_xxyyzz[k] = g1_0_xxyyzzz[k] - fr * g2_0_xxyyzz[k];

                        g_z_xxyzzz[k] = g1_0_xxyzzzz[k] - fr * g2_0_xxyzzz[k];

                        g_z_xxzzzz[k] = g1_0_xxzzzzz[k] - fr * g2_0_xxzzzz[k];

                        g_z_xyyyyy[k] = g1_0_xyyyyyz[k] - fr * g2_0_xyyyyy[k];

                        g_z_xyyyyz[k] = g1_0_xyyyyzz[k] - fr * g2_0_xyyyyz[k];

                        g_z_xyyyzz[k] = g1_0_xyyyzzz[k] - fr * g2_0_xyyyzz[k];

                        g_z_xyyzzz[k] = g1_0_xyyzzzz[k] - fr * g2_0_xyyzzz[k];

                        g_z_xyzzzz[k] = g1_0_xyzzzzz[k] - fr * g2_0_xyzzzz[k];

                        g_z_xzzzzz[k] = g1_0_xzzzzzz[k] - fr * g2_0_xzzzzz[k];

                        g_z_yyyyyy[k] = g1_0_yyyyyyz[k] - fr * g2_0_yyyyyy[k];

                        g_z_yyyyyz[k] = g1_0_yyyyyzz[k] - fr * g2_0_yyyyyz[k];

                        g_z_yyyyzz[k] = g1_0_yyyyzzz[k] - fr * g2_0_yyyyzz[k];

                        g_z_yyyzzz[k] = g1_0_yyyzzzz[k] - fr * g2_0_yyyzzz[k];

                        g_z_yyzzzz[k] = g1_0_yyzzzzz[k] - fr * g2_0_yyzzzz[k];

                        g_z_yzzzzz[k] = g1_0_yzzzzzz[k] - fr * g2_0_yzzzzz[k];

                        g_z_zzzzzz[k] = g1_0_zzzzzzz[k] - fr * g2_0_zzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXPK(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 1) && (recPattern[i].third() == 7))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|17)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 1, 7});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 8});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 0, 7});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|SK)^(m) integrals

                    auto g2_0_xxxxxxx = contrBuffer.data(g2off + 36 * j);

                    auto g2_0_xxxxxxy = contrBuffer.data(g2off + 36 * j + 1);

                    auto g2_0_xxxxxxz = contrBuffer.data(g2off + 36 * j + 2);

                    auto g2_0_xxxxxyy = contrBuffer.data(g2off + 36 * j + 3);

                    auto g2_0_xxxxxyz = contrBuffer.data(g2off + 36 * j + 4);

                    auto g2_0_xxxxxzz = contrBuffer.data(g2off + 36 * j + 5);

                    auto g2_0_xxxxyyy = contrBuffer.data(g2off + 36 * j + 6);

                    auto g2_0_xxxxyyz = contrBuffer.data(g2off + 36 * j + 7);

                    auto g2_0_xxxxyzz = contrBuffer.data(g2off + 36 * j + 8);

                    auto g2_0_xxxxzzz = contrBuffer.data(g2off + 36 * j + 9);

                    auto g2_0_xxxyyyy = contrBuffer.data(g2off + 36 * j + 10);

                    auto g2_0_xxxyyyz = contrBuffer.data(g2off + 36 * j + 11);

                    auto g2_0_xxxyyzz = contrBuffer.data(g2off + 36 * j + 12);

                    auto g2_0_xxxyzzz = contrBuffer.data(g2off + 36 * j + 13);

                    auto g2_0_xxxzzzz = contrBuffer.data(g2off + 36 * j + 14);

                    auto g2_0_xxyyyyy = contrBuffer.data(g2off + 36 * j + 15);

                    auto g2_0_xxyyyyz = contrBuffer.data(g2off + 36 * j + 16);

                    auto g2_0_xxyyyzz = contrBuffer.data(g2off + 36 * j + 17);

                    auto g2_0_xxyyzzz = contrBuffer.data(g2off + 36 * j + 18);

                    auto g2_0_xxyzzzz = contrBuffer.data(g2off + 36 * j + 19);

                    auto g2_0_xxzzzzz = contrBuffer.data(g2off + 36 * j + 20);

                    auto g2_0_xyyyyyy = contrBuffer.data(g2off + 36 * j + 21);

                    auto g2_0_xyyyyyz = contrBuffer.data(g2off + 36 * j + 22);

                    auto g2_0_xyyyyzz = contrBuffer.data(g2off + 36 * j + 23);

                    auto g2_0_xyyyzzz = contrBuffer.data(g2off + 36 * j + 24);

                    auto g2_0_xyyzzzz = contrBuffer.data(g2off + 36 * j + 25);

                    auto g2_0_xyzzzzz = contrBuffer.data(g2off + 36 * j + 26);

                    auto g2_0_xzzzzzz = contrBuffer.data(g2off + 36 * j + 27);

                    auto g2_0_yyyyyyy = contrBuffer.data(g2off + 36 * j + 28);

                    auto g2_0_yyyyyyz = contrBuffer.data(g2off + 36 * j + 29);

                    auto g2_0_yyyyyzz = contrBuffer.data(g2off + 36 * j + 30);

                    auto g2_0_yyyyzzz = contrBuffer.data(g2off + 36 * j + 31);

                    auto g2_0_yyyzzzz = contrBuffer.data(g2off + 36 * j + 32);

                    auto g2_0_yyzzzzz = contrBuffer.data(g2off + 36 * j + 33);

                    auto g2_0_yzzzzzz = contrBuffer.data(g2off + 36 * j + 34);

                    auto g2_0_zzzzzzz = contrBuffer.data(g2off + 36 * j + 35);

                    // set up pointers to (SX|g(r,r')|SL)^(m) integrals

                    auto g1_0_xxxxxxxx = contrBuffer.data(g1off + 45 * j);

                    auto g1_0_xxxxxxxy = contrBuffer.data(g1off + 45 * j + 1);

                    auto g1_0_xxxxxxxz = contrBuffer.data(g1off + 45 * j + 2);

                    auto g1_0_xxxxxxyy = contrBuffer.data(g1off + 45 * j + 3);

                    auto g1_0_xxxxxxyz = contrBuffer.data(g1off + 45 * j + 4);

                    auto g1_0_xxxxxxzz = contrBuffer.data(g1off + 45 * j + 5);

                    auto g1_0_xxxxxyyy = contrBuffer.data(g1off + 45 * j + 6);

                    auto g1_0_xxxxxyyz = contrBuffer.data(g1off + 45 * j + 7);

                    auto g1_0_xxxxxyzz = contrBuffer.data(g1off + 45 * j + 8);

                    auto g1_0_xxxxxzzz = contrBuffer.data(g1off + 45 * j + 9);

                    auto g1_0_xxxxyyyy = contrBuffer.data(g1off + 45 * j + 10);

                    auto g1_0_xxxxyyyz = contrBuffer.data(g1off + 45 * j + 11);

                    auto g1_0_xxxxyyzz = contrBuffer.data(g1off + 45 * j + 12);

                    auto g1_0_xxxxyzzz = contrBuffer.data(g1off + 45 * j + 13);

                    auto g1_0_xxxxzzzz = contrBuffer.data(g1off + 45 * j + 14);

                    auto g1_0_xxxyyyyy = contrBuffer.data(g1off + 45 * j + 15);

                    auto g1_0_xxxyyyyz = contrBuffer.data(g1off + 45 * j + 16);

                    auto g1_0_xxxyyyzz = contrBuffer.data(g1off + 45 * j + 17);

                    auto g1_0_xxxyyzzz = contrBuffer.data(g1off + 45 * j + 18);

                    auto g1_0_xxxyzzzz = contrBuffer.data(g1off + 45 * j + 19);

                    auto g1_0_xxxzzzzz = contrBuffer.data(g1off + 45 * j + 20);

                    auto g1_0_xxyyyyyy = contrBuffer.data(g1off + 45 * j + 21);

                    auto g1_0_xxyyyyyz = contrBuffer.data(g1off + 45 * j + 22);

                    auto g1_0_xxyyyyzz = contrBuffer.data(g1off + 45 * j + 23);

                    auto g1_0_xxyyyzzz = contrBuffer.data(g1off + 45 * j + 24);

                    auto g1_0_xxyyzzzz = contrBuffer.data(g1off + 45 * j + 25);

                    auto g1_0_xxyzzzzz = contrBuffer.data(g1off + 45 * j + 26);

                    auto g1_0_xxzzzzzz = contrBuffer.data(g1off + 45 * j + 27);

                    auto g1_0_xyyyyyyy = contrBuffer.data(g1off + 45 * j + 28);

                    auto g1_0_xyyyyyyz = contrBuffer.data(g1off + 45 * j + 29);

                    auto g1_0_xyyyyyzz = contrBuffer.data(g1off + 45 * j + 30);

                    auto g1_0_xyyyyzzz = contrBuffer.data(g1off + 45 * j + 31);

                    auto g1_0_xyyyzzzz = contrBuffer.data(g1off + 45 * j + 32);

                    auto g1_0_xyyzzzzz = contrBuffer.data(g1off + 45 * j + 33);

                    auto g1_0_xyzzzzzz = contrBuffer.data(g1off + 45 * j + 34);

                    auto g1_0_xzzzzzzz = contrBuffer.data(g1off + 45 * j + 35);

                    auto g1_0_yyyyyyyy = contrBuffer.data(g1off + 45 * j + 36);

                    auto g1_0_yyyyyyyz = contrBuffer.data(g1off + 45 * j + 37);

                    auto g1_0_yyyyyyzz = contrBuffer.data(g1off + 45 * j + 38);

                    auto g1_0_yyyyyzzz = contrBuffer.data(g1off + 45 * j + 39);

                    auto g1_0_yyyyzzzz = contrBuffer.data(g1off + 45 * j + 40);

                    auto g1_0_yyyzzzzz = contrBuffer.data(g1off + 45 * j + 41);

                    auto g1_0_yyzzzzzz = contrBuffer.data(g1off + 45 * j + 42);

                    auto g1_0_yzzzzzzz = contrBuffer.data(g1off + 45 * j + 43);

                    auto g1_0_zzzzzzzz = contrBuffer.data(g1off + 45 * j + 44);

                    // set up pointers to (SX|g(r,r')|PK)^(m) integrals

                    auto g_x_xxxxxxx = contrBuffer.data(goff + 108 * j);

                    auto g_x_xxxxxxy = contrBuffer.data(goff + 108 * j + 1);

                    auto g_x_xxxxxxz = contrBuffer.data(goff + 108 * j + 2);

                    auto g_x_xxxxxyy = contrBuffer.data(goff + 108 * j + 3);

                    auto g_x_xxxxxyz = contrBuffer.data(goff + 108 * j + 4);

                    auto g_x_xxxxxzz = contrBuffer.data(goff + 108 * j + 5);

                    auto g_x_xxxxyyy = contrBuffer.data(goff + 108 * j + 6);

                    auto g_x_xxxxyyz = contrBuffer.data(goff + 108 * j + 7);

                    auto g_x_xxxxyzz = contrBuffer.data(goff + 108 * j + 8);

                    auto g_x_xxxxzzz = contrBuffer.data(goff + 108 * j + 9);

                    auto g_x_xxxyyyy = contrBuffer.data(goff + 108 * j + 10);

                    auto g_x_xxxyyyz = contrBuffer.data(goff + 108 * j + 11);

                    auto g_x_xxxyyzz = contrBuffer.data(goff + 108 * j + 12);

                    auto g_x_xxxyzzz = contrBuffer.data(goff + 108 * j + 13);

                    auto g_x_xxxzzzz = contrBuffer.data(goff + 108 * j + 14);

                    auto g_x_xxyyyyy = contrBuffer.data(goff + 108 * j + 15);

                    auto g_x_xxyyyyz = contrBuffer.data(goff + 108 * j + 16);

                    auto g_x_xxyyyzz = contrBuffer.data(goff + 108 * j + 17);

                    auto g_x_xxyyzzz = contrBuffer.data(goff + 108 * j + 18);

                    auto g_x_xxyzzzz = contrBuffer.data(goff + 108 * j + 19);

                    auto g_x_xxzzzzz = contrBuffer.data(goff + 108 * j + 20);

                    auto g_x_xyyyyyy = contrBuffer.data(goff + 108 * j + 21);

                    auto g_x_xyyyyyz = contrBuffer.data(goff + 108 * j + 22);

                    auto g_x_xyyyyzz = contrBuffer.data(goff + 108 * j + 23);

                    auto g_x_xyyyzzz = contrBuffer.data(goff + 108 * j + 24);

                    auto g_x_xyyzzzz = contrBuffer.data(goff + 108 * j + 25);

                    auto g_x_xyzzzzz = contrBuffer.data(goff + 108 * j + 26);

                    auto g_x_xzzzzzz = contrBuffer.data(goff + 108 * j + 27);

                    auto g_x_yyyyyyy = contrBuffer.data(goff + 108 * j + 28);

                    auto g_x_yyyyyyz = contrBuffer.data(goff + 108 * j + 29);

                    auto g_x_yyyyyzz = contrBuffer.data(goff + 108 * j + 30);

                    auto g_x_yyyyzzz = contrBuffer.data(goff + 108 * j + 31);

                    auto g_x_yyyzzzz = contrBuffer.data(goff + 108 * j + 32);

                    auto g_x_yyzzzzz = contrBuffer.data(goff + 108 * j + 33);

                    auto g_x_yzzzzzz = contrBuffer.data(goff + 108 * j + 34);

                    auto g_x_zzzzzzz = contrBuffer.data(goff + 108 * j + 35);

                    auto g_y_xxxxxxx = contrBuffer.data(goff + 108 * j + 36);

                    auto g_y_xxxxxxy = contrBuffer.data(goff + 108 * j + 37);

                    auto g_y_xxxxxxz = contrBuffer.data(goff + 108 * j + 38);

                    auto g_y_xxxxxyy = contrBuffer.data(goff + 108 * j + 39);

                    auto g_y_xxxxxyz = contrBuffer.data(goff + 108 * j + 40);

                    auto g_y_xxxxxzz = contrBuffer.data(goff + 108 * j + 41);

                    auto g_y_xxxxyyy = contrBuffer.data(goff + 108 * j + 42);

                    auto g_y_xxxxyyz = contrBuffer.data(goff + 108 * j + 43);

                    auto g_y_xxxxyzz = contrBuffer.data(goff + 108 * j + 44);

                    auto g_y_xxxxzzz = contrBuffer.data(goff + 108 * j + 45);

                    auto g_y_xxxyyyy = contrBuffer.data(goff + 108 * j + 46);

                    auto g_y_xxxyyyz = contrBuffer.data(goff + 108 * j + 47);

                    auto g_y_xxxyyzz = contrBuffer.data(goff + 108 * j + 48);

                    auto g_y_xxxyzzz = contrBuffer.data(goff + 108 * j + 49);

                    auto g_y_xxxzzzz = contrBuffer.data(goff + 108 * j + 50);

                    auto g_y_xxyyyyy = contrBuffer.data(goff + 108 * j + 51);

                    auto g_y_xxyyyyz = contrBuffer.data(goff + 108 * j + 52);

                    auto g_y_xxyyyzz = contrBuffer.data(goff + 108 * j + 53);

                    auto g_y_xxyyzzz = contrBuffer.data(goff + 108 * j + 54);

                    auto g_y_xxyzzzz = contrBuffer.data(goff + 108 * j + 55);

                    auto g_y_xxzzzzz = contrBuffer.data(goff + 108 * j + 56);

                    auto g_y_xyyyyyy = contrBuffer.data(goff + 108 * j + 57);

                    auto g_y_xyyyyyz = contrBuffer.data(goff + 108 * j + 58);

                    auto g_y_xyyyyzz = contrBuffer.data(goff + 108 * j + 59);

                    auto g_y_xyyyzzz = contrBuffer.data(goff + 108 * j + 60);

                    auto g_y_xyyzzzz = contrBuffer.data(goff + 108 * j + 61);

                    auto g_y_xyzzzzz = contrBuffer.data(goff + 108 * j + 62);

                    auto g_y_xzzzzzz = contrBuffer.data(goff + 108 * j + 63);

                    auto g_y_yyyyyyy = contrBuffer.data(goff + 108 * j + 64);

                    auto g_y_yyyyyyz = contrBuffer.data(goff + 108 * j + 65);

                    auto g_y_yyyyyzz = contrBuffer.data(goff + 108 * j + 66);

                    auto g_y_yyyyzzz = contrBuffer.data(goff + 108 * j + 67);

                    auto g_y_yyyzzzz = contrBuffer.data(goff + 108 * j + 68);

                    auto g_y_yyzzzzz = contrBuffer.data(goff + 108 * j + 69);

                    auto g_y_yzzzzzz = contrBuffer.data(goff + 108 * j + 70);

                    auto g_y_zzzzzzz = contrBuffer.data(goff + 108 * j + 71);

                    auto g_z_xxxxxxx = contrBuffer.data(goff + 108 * j + 72);

                    auto g_z_xxxxxxy = contrBuffer.data(goff + 108 * j + 73);

                    auto g_z_xxxxxxz = contrBuffer.data(goff + 108 * j + 74);

                    auto g_z_xxxxxyy = contrBuffer.data(goff + 108 * j + 75);

                    auto g_z_xxxxxyz = contrBuffer.data(goff + 108 * j + 76);

                    auto g_z_xxxxxzz = contrBuffer.data(goff + 108 * j + 77);

                    auto g_z_xxxxyyy = contrBuffer.data(goff + 108 * j + 78);

                    auto g_z_xxxxyyz = contrBuffer.data(goff + 108 * j + 79);

                    auto g_z_xxxxyzz = contrBuffer.data(goff + 108 * j + 80);

                    auto g_z_xxxxzzz = contrBuffer.data(goff + 108 * j + 81);

                    auto g_z_xxxyyyy = contrBuffer.data(goff + 108 * j + 82);

                    auto g_z_xxxyyyz = contrBuffer.data(goff + 108 * j + 83);

                    auto g_z_xxxyyzz = contrBuffer.data(goff + 108 * j + 84);

                    auto g_z_xxxyzzz = contrBuffer.data(goff + 108 * j + 85);

                    auto g_z_xxxzzzz = contrBuffer.data(goff + 108 * j + 86);

                    auto g_z_xxyyyyy = contrBuffer.data(goff + 108 * j + 87);

                    auto g_z_xxyyyyz = contrBuffer.data(goff + 108 * j + 88);

                    auto g_z_xxyyyzz = contrBuffer.data(goff + 108 * j + 89);

                    auto g_z_xxyyzzz = contrBuffer.data(goff + 108 * j + 90);

                    auto g_z_xxyzzzz = contrBuffer.data(goff + 108 * j + 91);

                    auto g_z_xxzzzzz = contrBuffer.data(goff + 108 * j + 92);

                    auto g_z_xyyyyyy = contrBuffer.data(goff + 108 * j + 93);

                    auto g_z_xyyyyyz = contrBuffer.data(goff + 108 * j + 94);

                    auto g_z_xyyyyzz = contrBuffer.data(goff + 108 * j + 95);

                    auto g_z_xyyyzzz = contrBuffer.data(goff + 108 * j + 96);

                    auto g_z_xyyzzzz = contrBuffer.data(goff + 108 * j + 97);

                    auto g_z_xyzzzzz = contrBuffer.data(goff + 108 * j + 98);

                    auto g_z_xzzzzzz = contrBuffer.data(goff + 108 * j + 99);

                    auto g_z_yyyyyyy = contrBuffer.data(goff + 108 * j + 100);

                    auto g_z_yyyyyyz = contrBuffer.data(goff + 108 * j + 101);

                    auto g_z_yyyyyzz = contrBuffer.data(goff + 108 * j + 102);

                    auto g_z_yyyyzzz = contrBuffer.data(goff + 108 * j + 103);

                    auto g_z_yyyzzzz = contrBuffer.data(goff + 108 * j + 104);

                    auto g_z_yyzzzzz = contrBuffer.data(goff + 108 * j + 105);

                    auto g_z_yzzzzzz = contrBuffer.data(goff + 108 * j + 106);

                    auto g_z_zzzzzzz = contrBuffer.data(goff + 108 * j + 107);

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
                                             g1_0_zzzzzzzz, g_x_xxxxxxx, g_x_xxxxxxy,\
                                             g_x_xxxxxxz, g_x_xxxxxyy, g_x_xxxxxyz,\
                                             g_x_xxxxxzz, g_x_xxxxyyy, g_x_xxxxyyz,\
                                             g_x_xxxxyzz, g_x_xxxxzzz, g_x_xxxyyyy,\
                                             g_x_xxxyyyz, g_x_xxxyyzz, g_x_xxxyzzz,\
                                             g_x_xxxzzzz, g_x_xxyyyyy, g_x_xxyyyyz,\
                                             g_x_xxyyyzz, g_x_xxyyzzz, g_x_xxyzzzz,\
                                             g_x_xxzzzzz, g_x_xyyyyyy, g_x_xyyyyyz,\
                                             g_x_xyyyyzz, g_x_xyyyzzz, g_x_xyyzzzz,\
                                             g_x_xyzzzzz, g_x_xzzzzzz, g_x_yyyyyyy,\
                                             g_x_yyyyyyz, g_x_yyyyyzz, g_x_yyyyzzz,\
                                             g_x_yyyzzzz, g_x_yyzzzzz, g_x_yzzzzzz,\
                                             g_x_zzzzzzz, g_y_xxxxxxx, g_y_xxxxxxy,\
                                             g_y_xxxxxxz, g_y_xxxxxyy, g_y_xxxxxyz,\
                                             g_y_xxxxxzz, g_y_xxxxyyy, g_y_xxxxyyz,\
                                             g_y_xxxxyzz, g_y_xxxxzzz, g_y_xxxyyyy,\
                                             g_y_xxxyyyz, g_y_xxxyyzz, g_y_xxxyzzz,\
                                             g_y_xxxzzzz, g_y_xxyyyyy, g_y_xxyyyyz,\
                                             g_y_xxyyyzz, g_y_xxyyzzz, g_y_xxyzzzz,\
                                             g_y_xxzzzzz, g_y_xyyyyyy, g_y_xyyyyyz,\
                                             g_y_xyyyyzz, g_y_xyyyzzz, g_y_xyyzzzz,\
                                             g_y_xyzzzzz, g_y_xzzzzzz, g_y_yyyyyyy,\
                                             g_y_yyyyyyz, g_y_yyyyyzz, g_y_yyyyzzz,\
                                             g_y_yyyzzzz, g_y_yyzzzzz, g_y_yzzzzzz,\
                                             g_y_zzzzzzz, g_z_xxxxxxx, g_z_xxxxxxy,\
                                             g_z_xxxxxxz, g_z_xxxxxyy, g_z_xxxxxyz,\
                                             g_z_xxxxxzz, g_z_xxxxyyy, g_z_xxxxyyz,\
                                             g_z_xxxxyzz, g_z_xxxxzzz, g_z_xxxyyyy,\
                                             g_z_xxxyyyz, g_z_xxxyyzz, g_z_xxxyzzz,\
                                             g_z_xxxzzzz, g_z_xxyyyyy, g_z_xxyyyyz,\
                                             g_z_xxyyyzz, g_z_xxyyzzz, g_z_xxyzzzz,\
                                             g_z_xxzzzzz, g_z_xyyyyyy, g_z_xyyyyyz,\
                                             g_z_xyyyyzz, g_z_xyyyzzz, g_z_xyyzzzz,\
                                             g_z_xyzzzzz, g_z_xzzzzzz, g_z_yyyyyyy,\
                                             g_z_yyyyyyz, g_z_yyyyyzz, g_z_yyyyzzz,\
                                             g_z_yyyzzzz, g_z_yyzzzzz, g_z_yzzzzzz,\
                                             g_z_zzzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_x_xxxxxxx[k] = g1_0_xxxxxxxx[k] - fr * g2_0_xxxxxxx[k];

                        g_x_xxxxxxy[k] = g1_0_xxxxxxxy[k] - fr * g2_0_xxxxxxy[k];

                        g_x_xxxxxxz[k] = g1_0_xxxxxxxz[k] - fr * g2_0_xxxxxxz[k];

                        g_x_xxxxxyy[k] = g1_0_xxxxxxyy[k] - fr * g2_0_xxxxxyy[k];

                        g_x_xxxxxyz[k] = g1_0_xxxxxxyz[k] - fr * g2_0_xxxxxyz[k];

                        g_x_xxxxxzz[k] = g1_0_xxxxxxzz[k] - fr * g2_0_xxxxxzz[k];

                        g_x_xxxxyyy[k] = g1_0_xxxxxyyy[k] - fr * g2_0_xxxxyyy[k];

                        g_x_xxxxyyz[k] = g1_0_xxxxxyyz[k] - fr * g2_0_xxxxyyz[k];

                        g_x_xxxxyzz[k] = g1_0_xxxxxyzz[k] - fr * g2_0_xxxxyzz[k];

                        g_x_xxxxzzz[k] = g1_0_xxxxxzzz[k] - fr * g2_0_xxxxzzz[k];

                        g_x_xxxyyyy[k] = g1_0_xxxxyyyy[k] - fr * g2_0_xxxyyyy[k];

                        g_x_xxxyyyz[k] = g1_0_xxxxyyyz[k] - fr * g2_0_xxxyyyz[k];

                        g_x_xxxyyzz[k] = g1_0_xxxxyyzz[k] - fr * g2_0_xxxyyzz[k];

                        g_x_xxxyzzz[k] = g1_0_xxxxyzzz[k] - fr * g2_0_xxxyzzz[k];

                        g_x_xxxzzzz[k] = g1_0_xxxxzzzz[k] - fr * g2_0_xxxzzzz[k];

                        g_x_xxyyyyy[k] = g1_0_xxxyyyyy[k] - fr * g2_0_xxyyyyy[k];

                        g_x_xxyyyyz[k] = g1_0_xxxyyyyz[k] - fr * g2_0_xxyyyyz[k];

                        g_x_xxyyyzz[k] = g1_0_xxxyyyzz[k] - fr * g2_0_xxyyyzz[k];

                        g_x_xxyyzzz[k] = g1_0_xxxyyzzz[k] - fr * g2_0_xxyyzzz[k];

                        g_x_xxyzzzz[k] = g1_0_xxxyzzzz[k] - fr * g2_0_xxyzzzz[k];

                        g_x_xxzzzzz[k] = g1_0_xxxzzzzz[k] - fr * g2_0_xxzzzzz[k];

                        g_x_xyyyyyy[k] = g1_0_xxyyyyyy[k] - fr * g2_0_xyyyyyy[k];

                        g_x_xyyyyyz[k] = g1_0_xxyyyyyz[k] - fr * g2_0_xyyyyyz[k];

                        g_x_xyyyyzz[k] = g1_0_xxyyyyzz[k] - fr * g2_0_xyyyyzz[k];

                        g_x_xyyyzzz[k] = g1_0_xxyyyzzz[k] - fr * g2_0_xyyyzzz[k];

                        g_x_xyyzzzz[k] = g1_0_xxyyzzzz[k] - fr * g2_0_xyyzzzz[k];

                        g_x_xyzzzzz[k] = g1_0_xxyzzzzz[k] - fr * g2_0_xyzzzzz[k];

                        g_x_xzzzzzz[k] = g1_0_xxzzzzzz[k] - fr * g2_0_xzzzzzz[k];

                        g_x_yyyyyyy[k] = g1_0_xyyyyyyy[k] - fr * g2_0_yyyyyyy[k];

                        g_x_yyyyyyz[k] = g1_0_xyyyyyyz[k] - fr * g2_0_yyyyyyz[k];

                        g_x_yyyyyzz[k] = g1_0_xyyyyyzz[k] - fr * g2_0_yyyyyzz[k];

                        g_x_yyyyzzz[k] = g1_0_xyyyyzzz[k] - fr * g2_0_yyyyzzz[k];

                        g_x_yyyzzzz[k] = g1_0_xyyyzzzz[k] - fr * g2_0_yyyzzzz[k];

                        g_x_yyzzzzz[k] = g1_0_xyyzzzzz[k] - fr * g2_0_yyzzzzz[k];

                        g_x_yzzzzzz[k] = g1_0_xyzzzzzz[k] - fr * g2_0_yzzzzzz[k];

                        g_x_zzzzzzz[k] = g1_0_xzzzzzzz[k] - fr * g2_0_zzzzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_y_xxxxxxx[k] = g1_0_xxxxxxxy[k] - fr * g2_0_xxxxxxx[k];

                        g_y_xxxxxxy[k] = g1_0_xxxxxxyy[k] - fr * g2_0_xxxxxxy[k];

                        g_y_xxxxxxz[k] = g1_0_xxxxxxyz[k] - fr * g2_0_xxxxxxz[k];

                        g_y_xxxxxyy[k] = g1_0_xxxxxyyy[k] - fr * g2_0_xxxxxyy[k];

                        g_y_xxxxxyz[k] = g1_0_xxxxxyyz[k] - fr * g2_0_xxxxxyz[k];

                        g_y_xxxxxzz[k] = g1_0_xxxxxyzz[k] - fr * g2_0_xxxxxzz[k];

                        g_y_xxxxyyy[k] = g1_0_xxxxyyyy[k] - fr * g2_0_xxxxyyy[k];

                        g_y_xxxxyyz[k] = g1_0_xxxxyyyz[k] - fr * g2_0_xxxxyyz[k];

                        g_y_xxxxyzz[k] = g1_0_xxxxyyzz[k] - fr * g2_0_xxxxyzz[k];

                        g_y_xxxxzzz[k] = g1_0_xxxxyzzz[k] - fr * g2_0_xxxxzzz[k];

                        g_y_xxxyyyy[k] = g1_0_xxxyyyyy[k] - fr * g2_0_xxxyyyy[k];

                        g_y_xxxyyyz[k] = g1_0_xxxyyyyz[k] - fr * g2_0_xxxyyyz[k];

                        g_y_xxxyyzz[k] = g1_0_xxxyyyzz[k] - fr * g2_0_xxxyyzz[k];

                        g_y_xxxyzzz[k] = g1_0_xxxyyzzz[k] - fr * g2_0_xxxyzzz[k];

                        g_y_xxxzzzz[k] = g1_0_xxxyzzzz[k] - fr * g2_0_xxxzzzz[k];

                        g_y_xxyyyyy[k] = g1_0_xxyyyyyy[k] - fr * g2_0_xxyyyyy[k];

                        g_y_xxyyyyz[k] = g1_0_xxyyyyyz[k] - fr * g2_0_xxyyyyz[k];

                        g_y_xxyyyzz[k] = g1_0_xxyyyyzz[k] - fr * g2_0_xxyyyzz[k];

                        g_y_xxyyzzz[k] = g1_0_xxyyyzzz[k] - fr * g2_0_xxyyzzz[k];

                        g_y_xxyzzzz[k] = g1_0_xxyyzzzz[k] - fr * g2_0_xxyzzzz[k];

                        g_y_xxzzzzz[k] = g1_0_xxyzzzzz[k] - fr * g2_0_xxzzzzz[k];

                        g_y_xyyyyyy[k] = g1_0_xyyyyyyy[k] - fr * g2_0_xyyyyyy[k];

                        g_y_xyyyyyz[k] = g1_0_xyyyyyyz[k] - fr * g2_0_xyyyyyz[k];

                        g_y_xyyyyzz[k] = g1_0_xyyyyyzz[k] - fr * g2_0_xyyyyzz[k];

                        g_y_xyyyzzz[k] = g1_0_xyyyyzzz[k] - fr * g2_0_xyyyzzz[k];

                        g_y_xyyzzzz[k] = g1_0_xyyyzzzz[k] - fr * g2_0_xyyzzzz[k];

                        g_y_xyzzzzz[k] = g1_0_xyyzzzzz[k] - fr * g2_0_xyzzzzz[k];

                        g_y_xzzzzzz[k] = g1_0_xyzzzzzz[k] - fr * g2_0_xzzzzzz[k];

                        g_y_yyyyyyy[k] = g1_0_yyyyyyyy[k] - fr * g2_0_yyyyyyy[k];

                        g_y_yyyyyyz[k] = g1_0_yyyyyyyz[k] - fr * g2_0_yyyyyyz[k];

                        g_y_yyyyyzz[k] = g1_0_yyyyyyzz[k] - fr * g2_0_yyyyyzz[k];

                        g_y_yyyyzzz[k] = g1_0_yyyyyzzz[k] - fr * g2_0_yyyyzzz[k];

                        g_y_yyyzzzz[k] = g1_0_yyyyzzzz[k] - fr * g2_0_yyyzzzz[k];

                        g_y_yyzzzzz[k] = g1_0_yyyzzzzz[k] - fr * g2_0_yyzzzzz[k];

                        g_y_yzzzzzz[k] = g1_0_yyzzzzzz[k] - fr * g2_0_yzzzzzz[k];

                        g_y_zzzzzzz[k] = g1_0_yzzzzzzz[k] - fr * g2_0_zzzzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_z_xxxxxxx[k] = g1_0_xxxxxxxz[k] - fr * g2_0_xxxxxxx[k];

                        g_z_xxxxxxy[k] = g1_0_xxxxxxyz[k] - fr * g2_0_xxxxxxy[k];

                        g_z_xxxxxxz[k] = g1_0_xxxxxxzz[k] - fr * g2_0_xxxxxxz[k];

                        g_z_xxxxxyy[k] = g1_0_xxxxxyyz[k] - fr * g2_0_xxxxxyy[k];

                        g_z_xxxxxyz[k] = g1_0_xxxxxyzz[k] - fr * g2_0_xxxxxyz[k];

                        g_z_xxxxxzz[k] = g1_0_xxxxxzzz[k] - fr * g2_0_xxxxxzz[k];

                        g_z_xxxxyyy[k] = g1_0_xxxxyyyz[k] - fr * g2_0_xxxxyyy[k];

                        g_z_xxxxyyz[k] = g1_0_xxxxyyzz[k] - fr * g2_0_xxxxyyz[k];

                        g_z_xxxxyzz[k] = g1_0_xxxxyzzz[k] - fr * g2_0_xxxxyzz[k];

                        g_z_xxxxzzz[k] = g1_0_xxxxzzzz[k] - fr * g2_0_xxxxzzz[k];

                        g_z_xxxyyyy[k] = g1_0_xxxyyyyz[k] - fr * g2_0_xxxyyyy[k];

                        g_z_xxxyyyz[k] = g1_0_xxxyyyzz[k] - fr * g2_0_xxxyyyz[k];

                        g_z_xxxyyzz[k] = g1_0_xxxyyzzz[k] - fr * g2_0_xxxyyzz[k];

                        g_z_xxxyzzz[k] = g1_0_xxxyzzzz[k] - fr * g2_0_xxxyzzz[k];

                        g_z_xxxzzzz[k] = g1_0_xxxzzzzz[k] - fr * g2_0_xxxzzzz[k];

                        g_z_xxyyyyy[k] = g1_0_xxyyyyyz[k] - fr * g2_0_xxyyyyy[k];

                        g_z_xxyyyyz[k] = g1_0_xxyyyyzz[k] - fr * g2_0_xxyyyyz[k];

                        g_z_xxyyyzz[k] = g1_0_xxyyyzzz[k] - fr * g2_0_xxyyyzz[k];

                        g_z_xxyyzzz[k] = g1_0_xxyyzzzz[k] - fr * g2_0_xxyyzzz[k];

                        g_z_xxyzzzz[k] = g1_0_xxyzzzzz[k] - fr * g2_0_xxyzzzz[k];

                        g_z_xxzzzzz[k] = g1_0_xxzzzzzz[k] - fr * g2_0_xxzzzzz[k];

                        g_z_xyyyyyy[k] = g1_0_xyyyyyyz[k] - fr * g2_0_xyyyyyy[k];

                        g_z_xyyyyyz[k] = g1_0_xyyyyyzz[k] - fr * g2_0_xyyyyyz[k];

                        g_z_xyyyyzz[k] = g1_0_xyyyyzzz[k] - fr * g2_0_xyyyyzz[k];

                        g_z_xyyyzzz[k] = g1_0_xyyyzzzz[k] - fr * g2_0_xyyyzzz[k];

                        g_z_xyyzzzz[k] = g1_0_xyyzzzzz[k] - fr * g2_0_xyyzzzz[k];

                        g_z_xyzzzzz[k] = g1_0_xyzzzzzz[k] - fr * g2_0_xyzzzzz[k];

                        g_z_xzzzzzz[k] = g1_0_xzzzzzzz[k] - fr * g2_0_xzzzzzz[k];

                        g_z_yyyyyyy[k] = g1_0_yyyyyyyz[k] - fr * g2_0_yyyyyyy[k];

                        g_z_yyyyyyz[k] = g1_0_yyyyyyzz[k] - fr * g2_0_yyyyyyz[k];

                        g_z_yyyyyzz[k] = g1_0_yyyyyzzz[k] - fr * g2_0_yyyyyzz[k];

                        g_z_yyyyzzz[k] = g1_0_yyyyzzzz[k] - fr * g2_0_yyyyzzz[k];

                        g_z_yyyzzzz[k] = g1_0_yyyzzzzz[k] - fr * g2_0_yyyzzzz[k];

                        g_z_yyzzzzz[k] = g1_0_yyzzzzzz[k] - fr * g2_0_yyzzzzz[k];

                        g_z_yzzzzzz[k] = g1_0_yzzzzzzz[k] - fr * g2_0_yzzzzzz[k];

                        g_z_zzzzzzz[k] = g1_0_zzzzzzzz[k] - fr * g2_0_zzzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXDD(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 2) && (recPattern[i].third() == 2))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|22)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 2, 2});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 3});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 2});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|PD)^(m) integrals

                    auto g2_x_xx = contrBuffer.data(g2off + 18 * j);

                    auto g2_x_xy = contrBuffer.data(g2off + 18 * j + 1);

                    auto g2_x_xz = contrBuffer.data(g2off + 18 * j + 2);

                    auto g2_x_yy = contrBuffer.data(g2off + 18 * j + 3);

                    auto g2_x_yz = contrBuffer.data(g2off + 18 * j + 4);

                    auto g2_x_zz = contrBuffer.data(g2off + 18 * j + 5);

                    auto g2_y_xx = contrBuffer.data(g2off + 18 * j + 6);

                    auto g2_y_xy = contrBuffer.data(g2off + 18 * j + 7);

                    auto g2_y_xz = contrBuffer.data(g2off + 18 * j + 8);

                    auto g2_y_yy = contrBuffer.data(g2off + 18 * j + 9);

                    auto g2_y_yz = contrBuffer.data(g2off + 18 * j + 10);

                    auto g2_y_zz = contrBuffer.data(g2off + 18 * j + 11);

                    auto g2_z_xx = contrBuffer.data(g2off + 18 * j + 12);

                    auto g2_z_xy = contrBuffer.data(g2off + 18 * j + 13);

                    auto g2_z_xz = contrBuffer.data(g2off + 18 * j + 14);

                    auto g2_z_yy = contrBuffer.data(g2off + 18 * j + 15);

                    auto g2_z_yz = contrBuffer.data(g2off + 18 * j + 16);

                    auto g2_z_zz = contrBuffer.data(g2off + 18 * j + 17);

                    // set up pointers to (SX|g(r,r')|PF)^(m) integrals

                    auto g1_x_xxx = contrBuffer.data(g1off + 30 * j);

                    auto g1_x_xxy = contrBuffer.data(g1off + 30 * j + 1);

                    auto g1_x_xxz = contrBuffer.data(g1off + 30 * j + 2);

                    auto g1_x_xyy = contrBuffer.data(g1off + 30 * j + 3);

                    auto g1_x_xyz = contrBuffer.data(g1off + 30 * j + 4);

                    auto g1_x_xzz = contrBuffer.data(g1off + 30 * j + 5);

                    auto g1_y_xxx = contrBuffer.data(g1off + 30 * j + 10);

                    auto g1_y_xxy = contrBuffer.data(g1off + 30 * j + 11);

                    auto g1_y_xxz = contrBuffer.data(g1off + 30 * j + 12);

                    auto g1_y_xyy = contrBuffer.data(g1off + 30 * j + 13);

                    auto g1_y_xyz = contrBuffer.data(g1off + 30 * j + 14);

                    auto g1_y_xzz = contrBuffer.data(g1off + 30 * j + 15);

                    auto g1_y_yyy = contrBuffer.data(g1off + 30 * j + 16);

                    auto g1_y_yyz = contrBuffer.data(g1off + 30 * j + 17);

                    auto g1_y_yzz = contrBuffer.data(g1off + 30 * j + 18);

                    auto g1_z_xxx = contrBuffer.data(g1off + 30 * j + 20);

                    auto g1_z_xxy = contrBuffer.data(g1off + 30 * j + 21);

                    auto g1_z_xxz = contrBuffer.data(g1off + 30 * j + 22);

                    auto g1_z_xyy = contrBuffer.data(g1off + 30 * j + 23);

                    auto g1_z_xyz = contrBuffer.data(g1off + 30 * j + 24);

                    auto g1_z_xzz = contrBuffer.data(g1off + 30 * j + 25);

                    auto g1_z_yyy = contrBuffer.data(g1off + 30 * j + 26);

                    auto g1_z_yyz = contrBuffer.data(g1off + 30 * j + 27);

                    auto g1_z_yzz = contrBuffer.data(g1off + 30 * j + 28);

                    auto g1_z_zzz = contrBuffer.data(g1off + 30 * j + 29);

                    // set up pointers to (SX|g(r,r')|DD)^(m) integrals

                    auto g_xx_xx = contrBuffer.data(goff + 36 * j);

                    auto g_xx_xy = contrBuffer.data(goff + 36 * j + 1);

                    auto g_xx_xz = contrBuffer.data(goff + 36 * j + 2);

                    auto g_xx_yy = contrBuffer.data(goff + 36 * j + 3);

                    auto g_xx_yz = contrBuffer.data(goff + 36 * j + 4);

                    auto g_xx_zz = contrBuffer.data(goff + 36 * j + 5);

                    auto g_xy_xx = contrBuffer.data(goff + 36 * j + 6);

                    auto g_xy_xy = contrBuffer.data(goff + 36 * j + 7);

                    auto g_xy_xz = contrBuffer.data(goff + 36 * j + 8);

                    auto g_xy_yy = contrBuffer.data(goff + 36 * j + 9);

                    auto g_xy_yz = contrBuffer.data(goff + 36 * j + 10);

                    auto g_xy_zz = contrBuffer.data(goff + 36 * j + 11);

                    auto g_xz_xx = contrBuffer.data(goff + 36 * j + 12);

                    auto g_xz_xy = contrBuffer.data(goff + 36 * j + 13);

                    auto g_xz_xz = contrBuffer.data(goff + 36 * j + 14);

                    auto g_xz_yy = contrBuffer.data(goff + 36 * j + 15);

                    auto g_xz_yz = contrBuffer.data(goff + 36 * j + 16);

                    auto g_xz_zz = contrBuffer.data(goff + 36 * j + 17);

                    auto g_yy_xx = contrBuffer.data(goff + 36 * j + 18);

                    auto g_yy_xy = contrBuffer.data(goff + 36 * j + 19);

                    auto g_yy_xz = contrBuffer.data(goff + 36 * j + 20);

                    auto g_yy_yy = contrBuffer.data(goff + 36 * j + 21);

                    auto g_yy_yz = contrBuffer.data(goff + 36 * j + 22);

                    auto g_yy_zz = contrBuffer.data(goff + 36 * j + 23);

                    auto g_yz_xx = contrBuffer.data(goff + 36 * j + 24);

                    auto g_yz_xy = contrBuffer.data(goff + 36 * j + 25);

                    auto g_yz_xz = contrBuffer.data(goff + 36 * j + 26);

                    auto g_yz_yy = contrBuffer.data(goff + 36 * j + 27);

                    auto g_yz_yz = contrBuffer.data(goff + 36 * j + 28);

                    auto g_yz_zz = contrBuffer.data(goff + 36 * j + 29);

                    auto g_zz_xx = contrBuffer.data(goff + 36 * j + 30);

                    auto g_zz_xy = contrBuffer.data(goff + 36 * j + 31);

                    auto g_zz_xz = contrBuffer.data(goff + 36 * j + 32);

                    auto g_zz_yy = contrBuffer.data(goff + 36 * j + 33);

                    auto g_zz_yz = contrBuffer.data(goff + 36 * j + 34);

                    auto g_zz_zz = contrBuffer.data(goff + 36 * j + 35);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xx, g2_x_xy,\
                                             g2_x_xz, g2_x_yy, g2_x_yz, g2_x_zz,\
                                             g2_y_xx, g2_y_xy, g2_y_xz, g2_y_yy,\
                                             g2_y_yz, g2_y_zz, g2_z_xx, g2_z_xy,\
                                             g2_z_xz, g2_z_yy, g2_z_yz, g2_z_zz,\
                                             g1_x_xxx, g1_x_xxy, g1_x_xxz, g1_x_xyy,\
                                             g1_x_xyz, g1_x_xzz, g1_y_xxx, g1_y_xxy,\
                                             g1_y_xxz, g1_y_xyy, g1_y_xyz, g1_y_xzz,\
                                             g1_y_yyy, g1_y_yyz, g1_y_yzz,\
                                             g1_z_xxx, g1_z_xxy, g1_z_xxz, g1_z_xyy,\
                                             g1_z_xyz, g1_z_xzz, g1_z_yyy, g1_z_yyz,\
                                             g1_z_yzz, g1_z_zzz, g_xx_xx, g_xx_xy,\
                                             g_xx_xz, g_xx_yy, g_xx_yz, g_xx_zz,\
                                             g_xy_xx, g_xy_xy, g_xy_xz, g_xy_yy,\
                                             g_xy_yz, g_xy_zz, g_xz_xx, g_xz_xy,\
                                             g_xz_xz, g_xz_yy, g_xz_yz, g_xz_zz,\
                                             g_yy_xx, g_yy_xy, g_yy_xz, g_yy_yy,\
                                             g_yy_yz, g_yy_zz, g_yz_xx, g_yz_xy,\
                                             g_yz_xz, g_yz_yy, g_yz_yz, g_yz_zz,\
                                             g_zz_xx, g_zz_xy, g_zz_xz, g_zz_yy,\
                                             g_zz_yz, g_zz_zz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xx_xx[k] = g1_x_xxx[k] - fr * g2_x_xx[k];

                        g_xx_xy[k] = g1_x_xxy[k] - fr * g2_x_xy[k];

                        g_xx_xz[k] = g1_x_xxz[k] - fr * g2_x_xz[k];

                        g_xx_yy[k] = g1_x_xyy[k] - fr * g2_x_yy[k];

                        g_xx_yz[k] = g1_x_xyz[k] - fr * g2_x_yz[k];

                        g_xx_zz[k] = g1_x_xzz[k] - fr * g2_x_zz[k];

                        g_xy_xx[k] = g1_y_xxx[k] - fr * g2_y_xx[k];

                        g_xy_xy[k] = g1_y_xxy[k] - fr * g2_y_xy[k];

                        g_xy_xz[k] = g1_y_xxz[k] - fr * g2_y_xz[k];

                        g_xy_yy[k] = g1_y_xyy[k] - fr * g2_y_yy[k];

                        g_xy_yz[k] = g1_y_xyz[k] - fr * g2_y_yz[k];

                        g_xy_zz[k] = g1_y_xzz[k] - fr * g2_y_zz[k];

                        g_xz_xx[k] = g1_z_xxx[k] - fr * g2_z_xx[k];

                        g_xz_xy[k] = g1_z_xxy[k] - fr * g2_z_xy[k];

                        g_xz_xz[k] = g1_z_xxz[k] - fr * g2_z_xz[k];

                        g_xz_yy[k] = g1_z_xyy[k] - fr * g2_z_yy[k];

                        g_xz_yz[k] = g1_z_xyz[k] - fr * g2_z_yz[k];

                        g_xz_zz[k] = g1_z_xzz[k] - fr * g2_z_zz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yy_xx[k] = g1_y_xxy[k] - fr * g2_y_xx[k];

                        g_yy_xy[k] = g1_y_xyy[k] - fr * g2_y_xy[k];

                        g_yy_xz[k] = g1_y_xyz[k] - fr * g2_y_xz[k];

                        g_yy_yy[k] = g1_y_yyy[k] - fr * g2_y_yy[k];

                        g_yy_yz[k] = g1_y_yyz[k] - fr * g2_y_yz[k];

                        g_yy_zz[k] = g1_y_yzz[k] - fr * g2_y_zz[k];

                        g_yz_xx[k] = g1_z_xxy[k] - fr * g2_z_xx[k];

                        g_yz_xy[k] = g1_z_xyy[k] - fr * g2_z_xy[k];

                        g_yz_xz[k] = g1_z_xyz[k] - fr * g2_z_xz[k];

                        g_yz_yy[k] = g1_z_yyy[k] - fr * g2_z_yy[k];

                        g_yz_yz[k] = g1_z_yyz[k] - fr * g2_z_yz[k];

                        g_yz_zz[k] = g1_z_yzz[k] - fr * g2_z_zz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zz_xx[k] = g1_z_xxz[k] - fr * g2_z_xx[k];

                        g_zz_xy[k] = g1_z_xyz[k] - fr * g2_z_xy[k];

                        g_zz_xz[k] = g1_z_xzz[k] - fr * g2_z_xz[k];

                        g_zz_yy[k] = g1_z_yyz[k] - fr * g2_z_yy[k];

                        g_zz_yz[k] = g1_z_yzz[k] - fr * g2_z_yz[k];

                        g_zz_zz[k] = g1_z_zzz[k] - fr * g2_z_zz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXDF(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 2) && (recPattern[i].third() == 3))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|23)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 2, 3});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 4});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 3});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|PF)^(m) integrals

                    auto g2_x_xxx = contrBuffer.data(g2off + 30 * j);

                    auto g2_x_xxy = contrBuffer.data(g2off + 30 * j + 1);

                    auto g2_x_xxz = contrBuffer.data(g2off + 30 * j + 2);

                    auto g2_x_xyy = contrBuffer.data(g2off + 30 * j + 3);

                    auto g2_x_xyz = contrBuffer.data(g2off + 30 * j + 4);

                    auto g2_x_xzz = contrBuffer.data(g2off + 30 * j + 5);

                    auto g2_x_yyy = contrBuffer.data(g2off + 30 * j + 6);

                    auto g2_x_yyz = contrBuffer.data(g2off + 30 * j + 7);

                    auto g2_x_yzz = contrBuffer.data(g2off + 30 * j + 8);

                    auto g2_x_zzz = contrBuffer.data(g2off + 30 * j + 9);

                    auto g2_y_xxx = contrBuffer.data(g2off + 30 * j + 10);

                    auto g2_y_xxy = contrBuffer.data(g2off + 30 * j + 11);

                    auto g2_y_xxz = contrBuffer.data(g2off + 30 * j + 12);

                    auto g2_y_xyy = contrBuffer.data(g2off + 30 * j + 13);

                    auto g2_y_xyz = contrBuffer.data(g2off + 30 * j + 14);

                    auto g2_y_xzz = contrBuffer.data(g2off + 30 * j + 15);

                    auto g2_y_yyy = contrBuffer.data(g2off + 30 * j + 16);

                    auto g2_y_yyz = contrBuffer.data(g2off + 30 * j + 17);

                    auto g2_y_yzz = contrBuffer.data(g2off + 30 * j + 18);

                    auto g2_y_zzz = contrBuffer.data(g2off + 30 * j + 19);

                    auto g2_z_xxx = contrBuffer.data(g2off + 30 * j + 20);

                    auto g2_z_xxy = contrBuffer.data(g2off + 30 * j + 21);

                    auto g2_z_xxz = contrBuffer.data(g2off + 30 * j + 22);

                    auto g2_z_xyy = contrBuffer.data(g2off + 30 * j + 23);

                    auto g2_z_xyz = contrBuffer.data(g2off + 30 * j + 24);

                    auto g2_z_xzz = contrBuffer.data(g2off + 30 * j + 25);

                    auto g2_z_yyy = contrBuffer.data(g2off + 30 * j + 26);

                    auto g2_z_yyz = contrBuffer.data(g2off + 30 * j + 27);

                    auto g2_z_yzz = contrBuffer.data(g2off + 30 * j + 28);

                    auto g2_z_zzz = contrBuffer.data(g2off + 30 * j + 29);

                    // set up pointers to (SX|g(r,r')|PG)^(m) integrals

                    auto g1_x_xxxx = contrBuffer.data(g1off + 45 * j);

                    auto g1_x_xxxy = contrBuffer.data(g1off + 45 * j + 1);

                    auto g1_x_xxxz = contrBuffer.data(g1off + 45 * j + 2);

                    auto g1_x_xxyy = contrBuffer.data(g1off + 45 * j + 3);

                    auto g1_x_xxyz = contrBuffer.data(g1off + 45 * j + 4);

                    auto g1_x_xxzz = contrBuffer.data(g1off + 45 * j + 5);

                    auto g1_x_xyyy = contrBuffer.data(g1off + 45 * j + 6);

                    auto g1_x_xyyz = contrBuffer.data(g1off + 45 * j + 7);

                    auto g1_x_xyzz = contrBuffer.data(g1off + 45 * j + 8);

                    auto g1_x_xzzz = contrBuffer.data(g1off + 45 * j + 9);

                    auto g1_y_xxxx = contrBuffer.data(g1off + 45 * j + 15);

                    auto g1_y_xxxy = contrBuffer.data(g1off + 45 * j + 16);

                    auto g1_y_xxxz = contrBuffer.data(g1off + 45 * j + 17);

                    auto g1_y_xxyy = contrBuffer.data(g1off + 45 * j + 18);

                    auto g1_y_xxyz = contrBuffer.data(g1off + 45 * j + 19);

                    auto g1_y_xxzz = contrBuffer.data(g1off + 45 * j + 20);

                    auto g1_y_xyyy = contrBuffer.data(g1off + 45 * j + 21);

                    auto g1_y_xyyz = contrBuffer.data(g1off + 45 * j + 22);

                    auto g1_y_xyzz = contrBuffer.data(g1off + 45 * j + 23);

                    auto g1_y_xzzz = contrBuffer.data(g1off + 45 * j + 24);

                    auto g1_y_yyyy = contrBuffer.data(g1off + 45 * j + 25);

                    auto g1_y_yyyz = contrBuffer.data(g1off + 45 * j + 26);

                    auto g1_y_yyzz = contrBuffer.data(g1off + 45 * j + 27);

                    auto g1_y_yzzz = contrBuffer.data(g1off + 45 * j + 28);

                    auto g1_z_xxxx = contrBuffer.data(g1off + 45 * j + 30);

                    auto g1_z_xxxy = contrBuffer.data(g1off + 45 * j + 31);

                    auto g1_z_xxxz = contrBuffer.data(g1off + 45 * j + 32);

                    auto g1_z_xxyy = contrBuffer.data(g1off + 45 * j + 33);

                    auto g1_z_xxyz = contrBuffer.data(g1off + 45 * j + 34);

                    auto g1_z_xxzz = contrBuffer.data(g1off + 45 * j + 35);

                    auto g1_z_xyyy = contrBuffer.data(g1off + 45 * j + 36);

                    auto g1_z_xyyz = contrBuffer.data(g1off + 45 * j + 37);

                    auto g1_z_xyzz = contrBuffer.data(g1off + 45 * j + 38);

                    auto g1_z_xzzz = contrBuffer.data(g1off + 45 * j + 39);

                    auto g1_z_yyyy = contrBuffer.data(g1off + 45 * j + 40);

                    auto g1_z_yyyz = contrBuffer.data(g1off + 45 * j + 41);

                    auto g1_z_yyzz = contrBuffer.data(g1off + 45 * j + 42);

                    auto g1_z_yzzz = contrBuffer.data(g1off + 45 * j + 43);

                    auto g1_z_zzzz = contrBuffer.data(g1off + 45 * j + 44);

                    // set up pointers to (SX|g(r,r')|DF)^(m) integrals

                    auto g_xx_xxx = contrBuffer.data(goff + 60 * j);

                    auto g_xx_xxy = contrBuffer.data(goff + 60 * j + 1);

                    auto g_xx_xxz = contrBuffer.data(goff + 60 * j + 2);

                    auto g_xx_xyy = contrBuffer.data(goff + 60 * j + 3);

                    auto g_xx_xyz = contrBuffer.data(goff + 60 * j + 4);

                    auto g_xx_xzz = contrBuffer.data(goff + 60 * j + 5);

                    auto g_xx_yyy = contrBuffer.data(goff + 60 * j + 6);

                    auto g_xx_yyz = contrBuffer.data(goff + 60 * j + 7);

                    auto g_xx_yzz = contrBuffer.data(goff + 60 * j + 8);

                    auto g_xx_zzz = contrBuffer.data(goff + 60 * j + 9);

                    auto g_xy_xxx = contrBuffer.data(goff + 60 * j + 10);

                    auto g_xy_xxy = contrBuffer.data(goff + 60 * j + 11);

                    auto g_xy_xxz = contrBuffer.data(goff + 60 * j + 12);

                    auto g_xy_xyy = contrBuffer.data(goff + 60 * j + 13);

                    auto g_xy_xyz = contrBuffer.data(goff + 60 * j + 14);

                    auto g_xy_xzz = contrBuffer.data(goff + 60 * j + 15);

                    auto g_xy_yyy = contrBuffer.data(goff + 60 * j + 16);

                    auto g_xy_yyz = contrBuffer.data(goff + 60 * j + 17);

                    auto g_xy_yzz = contrBuffer.data(goff + 60 * j + 18);

                    auto g_xy_zzz = contrBuffer.data(goff + 60 * j + 19);

                    auto g_xz_xxx = contrBuffer.data(goff + 60 * j + 20);

                    auto g_xz_xxy = contrBuffer.data(goff + 60 * j + 21);

                    auto g_xz_xxz = contrBuffer.data(goff + 60 * j + 22);

                    auto g_xz_xyy = contrBuffer.data(goff + 60 * j + 23);

                    auto g_xz_xyz = contrBuffer.data(goff + 60 * j + 24);

                    auto g_xz_xzz = contrBuffer.data(goff + 60 * j + 25);

                    auto g_xz_yyy = contrBuffer.data(goff + 60 * j + 26);

                    auto g_xz_yyz = contrBuffer.data(goff + 60 * j + 27);

                    auto g_xz_yzz = contrBuffer.data(goff + 60 * j + 28);

                    auto g_xz_zzz = contrBuffer.data(goff + 60 * j + 29);

                    auto g_yy_xxx = contrBuffer.data(goff + 60 * j + 30);

                    auto g_yy_xxy = contrBuffer.data(goff + 60 * j + 31);

                    auto g_yy_xxz = contrBuffer.data(goff + 60 * j + 32);

                    auto g_yy_xyy = contrBuffer.data(goff + 60 * j + 33);

                    auto g_yy_xyz = contrBuffer.data(goff + 60 * j + 34);

                    auto g_yy_xzz = contrBuffer.data(goff + 60 * j + 35);

                    auto g_yy_yyy = contrBuffer.data(goff + 60 * j + 36);

                    auto g_yy_yyz = contrBuffer.data(goff + 60 * j + 37);

                    auto g_yy_yzz = contrBuffer.data(goff + 60 * j + 38);

                    auto g_yy_zzz = contrBuffer.data(goff + 60 * j + 39);

                    auto g_yz_xxx = contrBuffer.data(goff + 60 * j + 40);

                    auto g_yz_xxy = contrBuffer.data(goff + 60 * j + 41);

                    auto g_yz_xxz = contrBuffer.data(goff + 60 * j + 42);

                    auto g_yz_xyy = contrBuffer.data(goff + 60 * j + 43);

                    auto g_yz_xyz = contrBuffer.data(goff + 60 * j + 44);

                    auto g_yz_xzz = contrBuffer.data(goff + 60 * j + 45);

                    auto g_yz_yyy = contrBuffer.data(goff + 60 * j + 46);

                    auto g_yz_yyz = contrBuffer.data(goff + 60 * j + 47);

                    auto g_yz_yzz = contrBuffer.data(goff + 60 * j + 48);

                    auto g_yz_zzz = contrBuffer.data(goff + 60 * j + 49);

                    auto g_zz_xxx = contrBuffer.data(goff + 60 * j + 50);

                    auto g_zz_xxy = contrBuffer.data(goff + 60 * j + 51);

                    auto g_zz_xxz = contrBuffer.data(goff + 60 * j + 52);

                    auto g_zz_xyy = contrBuffer.data(goff + 60 * j + 53);

                    auto g_zz_xyz = contrBuffer.data(goff + 60 * j + 54);

                    auto g_zz_xzz = contrBuffer.data(goff + 60 * j + 55);

                    auto g_zz_yyy = contrBuffer.data(goff + 60 * j + 56);

                    auto g_zz_yyz = contrBuffer.data(goff + 60 * j + 57);

                    auto g_zz_yzz = contrBuffer.data(goff + 60 * j + 58);

                    auto g_zz_zzz = contrBuffer.data(goff + 60 * j + 59);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxx, g2_x_xxy,\
                                             g2_x_xxz, g2_x_xyy, g2_x_xyz, g2_x_xzz,\
                                             g2_x_yyy, g2_x_yyz, g2_x_yzz, g2_x_zzz,\
                                             g2_y_xxx, g2_y_xxy, g2_y_xxz, g2_y_xyy,\
                                             g2_y_xyz, g2_y_xzz, g2_y_yyy, g2_y_yyz,\
                                             g2_y_yzz, g2_y_zzz, g2_z_xxx, g2_z_xxy,\
                                             g2_z_xxz, g2_z_xyy, g2_z_xyz, g2_z_xzz,\
                                             g2_z_yyy, g2_z_yyz, g2_z_yzz, g2_z_zzz,\
                                             g1_x_xxxx, g1_x_xxxy, g1_x_xxxz, g1_x_xxyy,\
                                             g1_x_xxyz, g1_x_xxzz, g1_x_xyyy, g1_x_xyyz,\
                                             g1_x_xyzz, g1_x_xzzz, g1_y_xxxx,\
                                             g1_y_xxxy, g1_y_xxxz, g1_y_xxyy, g1_y_xxyz,\
                                             g1_y_xxzz, g1_y_xyyy, g1_y_xyyz, g1_y_xyzz,\
                                             g1_y_xzzz, g1_y_yyyy, g1_y_yyyz, g1_y_yyzz,\
                                             g1_y_yzzz, g1_z_xxxx, g1_z_xxxy,\
                                             g1_z_xxxz, g1_z_xxyy, g1_z_xxyz, g1_z_xxzz,\
                                             g1_z_xyyy, g1_z_xyyz, g1_z_xyzz, g1_z_xzzz,\
                                             g1_z_yyyy, g1_z_yyyz, g1_z_yyzz, g1_z_yzzz,\
                                             g1_z_zzzz, g_xx_xxx, g_xx_xxy, g_xx_xxz,\
                                             g_xx_xyy, g_xx_xyz, g_xx_xzz, g_xx_yyy,\
                                             g_xx_yyz, g_xx_yzz, g_xx_zzz, g_xy_xxx,\
                                             g_xy_xxy, g_xy_xxz, g_xy_xyy, g_xy_xyz,\
                                             g_xy_xzz, g_xy_yyy, g_xy_yyz, g_xy_yzz,\
                                             g_xy_zzz, g_xz_xxx, g_xz_xxy, g_xz_xxz,\
                                             g_xz_xyy, g_xz_xyz, g_xz_xzz, g_xz_yyy,\
                                             g_xz_yyz, g_xz_yzz, g_xz_zzz, g_yy_xxx,\
                                             g_yy_xxy, g_yy_xxz, g_yy_xyy, g_yy_xyz,\
                                             g_yy_xzz, g_yy_yyy, g_yy_yyz, g_yy_yzz,\
                                             g_yy_zzz, g_yz_xxx, g_yz_xxy, g_yz_xxz,\
                                             g_yz_xyy, g_yz_xyz, g_yz_xzz, g_yz_yyy,\
                                             g_yz_yyz, g_yz_yzz, g_yz_zzz, g_zz_xxx,\
                                             g_zz_xxy, g_zz_xxz, g_zz_xyy, g_zz_xyz,\
                                             g_zz_xzz, g_zz_yyy, g_zz_yyz, g_zz_yzz,\
                                             g_zz_zzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xx_xxx[k] = g1_x_xxxx[k] - fr * g2_x_xxx[k];

                        g_xx_xxy[k] = g1_x_xxxy[k] - fr * g2_x_xxy[k];

                        g_xx_xxz[k] = g1_x_xxxz[k] - fr * g2_x_xxz[k];

                        g_xx_xyy[k] = g1_x_xxyy[k] - fr * g2_x_xyy[k];

                        g_xx_xyz[k] = g1_x_xxyz[k] - fr * g2_x_xyz[k];

                        g_xx_xzz[k] = g1_x_xxzz[k] - fr * g2_x_xzz[k];

                        g_xx_yyy[k] = g1_x_xyyy[k] - fr * g2_x_yyy[k];

                        g_xx_yyz[k] = g1_x_xyyz[k] - fr * g2_x_yyz[k];

                        g_xx_yzz[k] = g1_x_xyzz[k] - fr * g2_x_yzz[k];

                        g_xx_zzz[k] = g1_x_xzzz[k] - fr * g2_x_zzz[k];

                        g_xy_xxx[k] = g1_y_xxxx[k] - fr * g2_y_xxx[k];

                        g_xy_xxy[k] = g1_y_xxxy[k] - fr * g2_y_xxy[k];

                        g_xy_xxz[k] = g1_y_xxxz[k] - fr * g2_y_xxz[k];

                        g_xy_xyy[k] = g1_y_xxyy[k] - fr * g2_y_xyy[k];

                        g_xy_xyz[k] = g1_y_xxyz[k] - fr * g2_y_xyz[k];

                        g_xy_xzz[k] = g1_y_xxzz[k] - fr * g2_y_xzz[k];

                        g_xy_yyy[k] = g1_y_xyyy[k] - fr * g2_y_yyy[k];

                        g_xy_yyz[k] = g1_y_xyyz[k] - fr * g2_y_yyz[k];

                        g_xy_yzz[k] = g1_y_xyzz[k] - fr * g2_y_yzz[k];

                        g_xy_zzz[k] = g1_y_xzzz[k] - fr * g2_y_zzz[k];

                        g_xz_xxx[k] = g1_z_xxxx[k] - fr * g2_z_xxx[k];

                        g_xz_xxy[k] = g1_z_xxxy[k] - fr * g2_z_xxy[k];

                        g_xz_xxz[k] = g1_z_xxxz[k] - fr * g2_z_xxz[k];

                        g_xz_xyy[k] = g1_z_xxyy[k] - fr * g2_z_xyy[k];

                        g_xz_xyz[k] = g1_z_xxyz[k] - fr * g2_z_xyz[k];

                        g_xz_xzz[k] = g1_z_xxzz[k] - fr * g2_z_xzz[k];

                        g_xz_yyy[k] = g1_z_xyyy[k] - fr * g2_z_yyy[k];

                        g_xz_yyz[k] = g1_z_xyyz[k] - fr * g2_z_yyz[k];

                        g_xz_yzz[k] = g1_z_xyzz[k] - fr * g2_z_yzz[k];

                        g_xz_zzz[k] = g1_z_xzzz[k] - fr * g2_z_zzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yy_xxx[k] = g1_y_xxxy[k] - fr * g2_y_xxx[k];

                        g_yy_xxy[k] = g1_y_xxyy[k] - fr * g2_y_xxy[k];

                        g_yy_xxz[k] = g1_y_xxyz[k] - fr * g2_y_xxz[k];

                        g_yy_xyy[k] = g1_y_xyyy[k] - fr * g2_y_xyy[k];

                        g_yy_xyz[k] = g1_y_xyyz[k] - fr * g2_y_xyz[k];

                        g_yy_xzz[k] = g1_y_xyzz[k] - fr * g2_y_xzz[k];

                        g_yy_yyy[k] = g1_y_yyyy[k] - fr * g2_y_yyy[k];

                        g_yy_yyz[k] = g1_y_yyyz[k] - fr * g2_y_yyz[k];

                        g_yy_yzz[k] = g1_y_yyzz[k] - fr * g2_y_yzz[k];

                        g_yy_zzz[k] = g1_y_yzzz[k] - fr * g2_y_zzz[k];

                        g_yz_xxx[k] = g1_z_xxxy[k] - fr * g2_z_xxx[k];

                        g_yz_xxy[k] = g1_z_xxyy[k] - fr * g2_z_xxy[k];

                        g_yz_xxz[k] = g1_z_xxyz[k] - fr * g2_z_xxz[k];

                        g_yz_xyy[k] = g1_z_xyyy[k] - fr * g2_z_xyy[k];

                        g_yz_xyz[k] = g1_z_xyyz[k] - fr * g2_z_xyz[k];

                        g_yz_xzz[k] = g1_z_xyzz[k] - fr * g2_z_xzz[k];

                        g_yz_yyy[k] = g1_z_yyyy[k] - fr * g2_z_yyy[k];

                        g_yz_yyz[k] = g1_z_yyyz[k] - fr * g2_z_yyz[k];

                        g_yz_yzz[k] = g1_z_yyzz[k] - fr * g2_z_yzz[k];

                        g_yz_zzz[k] = g1_z_yzzz[k] - fr * g2_z_zzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zz_xxx[k] = g1_z_xxxz[k] - fr * g2_z_xxx[k];

                        g_zz_xxy[k] = g1_z_xxyz[k] - fr * g2_z_xxy[k];

                        g_zz_xxz[k] = g1_z_xxzz[k] - fr * g2_z_xxz[k];

                        g_zz_xyy[k] = g1_z_xyyz[k] - fr * g2_z_xyy[k];

                        g_zz_xyz[k] = g1_z_xyzz[k] - fr * g2_z_xyz[k];

                        g_zz_xzz[k] = g1_z_xzzz[k] - fr * g2_z_xzz[k];

                        g_zz_yyy[k] = g1_z_yyyz[k] - fr * g2_z_yyy[k];

                        g_zz_yyz[k] = g1_z_yyzz[k] - fr * g2_z_yyz[k];

                        g_zz_yzz[k] = g1_z_yzzz[k] - fr * g2_z_yzz[k];

                        g_zz_zzz[k] = g1_z_zzzz[k] - fr * g2_z_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXDG(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 2) && (recPattern[i].third() == 4))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|24)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 2, 4});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 5});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 4});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|PG)^(m) integrals

                    auto g2_x_xxxx = contrBuffer.data(g2off + 45 * j);

                    auto g2_x_xxxy = contrBuffer.data(g2off + 45 * j + 1);

                    auto g2_x_xxxz = contrBuffer.data(g2off + 45 * j + 2);

                    auto g2_x_xxyy = contrBuffer.data(g2off + 45 * j + 3);

                    auto g2_x_xxyz = contrBuffer.data(g2off + 45 * j + 4);

                    auto g2_x_xxzz = contrBuffer.data(g2off + 45 * j + 5);

                    auto g2_x_xyyy = contrBuffer.data(g2off + 45 * j + 6);

                    auto g2_x_xyyz = contrBuffer.data(g2off + 45 * j + 7);

                    auto g2_x_xyzz = contrBuffer.data(g2off + 45 * j + 8);

                    auto g2_x_xzzz = contrBuffer.data(g2off + 45 * j + 9);

                    auto g2_x_yyyy = contrBuffer.data(g2off + 45 * j + 10);

                    auto g2_x_yyyz = contrBuffer.data(g2off + 45 * j + 11);

                    auto g2_x_yyzz = contrBuffer.data(g2off + 45 * j + 12);

                    auto g2_x_yzzz = contrBuffer.data(g2off + 45 * j + 13);

                    auto g2_x_zzzz = contrBuffer.data(g2off + 45 * j + 14);

                    auto g2_y_xxxx = contrBuffer.data(g2off + 45 * j + 15);

                    auto g2_y_xxxy = contrBuffer.data(g2off + 45 * j + 16);

                    auto g2_y_xxxz = contrBuffer.data(g2off + 45 * j + 17);

                    auto g2_y_xxyy = contrBuffer.data(g2off + 45 * j + 18);

                    auto g2_y_xxyz = contrBuffer.data(g2off + 45 * j + 19);

                    auto g2_y_xxzz = contrBuffer.data(g2off + 45 * j + 20);

                    auto g2_y_xyyy = contrBuffer.data(g2off + 45 * j + 21);

                    auto g2_y_xyyz = contrBuffer.data(g2off + 45 * j + 22);

                    auto g2_y_xyzz = contrBuffer.data(g2off + 45 * j + 23);

                    auto g2_y_xzzz = contrBuffer.data(g2off + 45 * j + 24);

                    auto g2_y_yyyy = contrBuffer.data(g2off + 45 * j + 25);

                    auto g2_y_yyyz = contrBuffer.data(g2off + 45 * j + 26);

                    auto g2_y_yyzz = contrBuffer.data(g2off + 45 * j + 27);

                    auto g2_y_yzzz = contrBuffer.data(g2off + 45 * j + 28);

                    auto g2_y_zzzz = contrBuffer.data(g2off + 45 * j + 29);

                    auto g2_z_xxxx = contrBuffer.data(g2off + 45 * j + 30);

                    auto g2_z_xxxy = contrBuffer.data(g2off + 45 * j + 31);

                    auto g2_z_xxxz = contrBuffer.data(g2off + 45 * j + 32);

                    auto g2_z_xxyy = contrBuffer.data(g2off + 45 * j + 33);

                    auto g2_z_xxyz = contrBuffer.data(g2off + 45 * j + 34);

                    auto g2_z_xxzz = contrBuffer.data(g2off + 45 * j + 35);

                    auto g2_z_xyyy = contrBuffer.data(g2off + 45 * j + 36);

                    auto g2_z_xyyz = contrBuffer.data(g2off + 45 * j + 37);

                    auto g2_z_xyzz = contrBuffer.data(g2off + 45 * j + 38);

                    auto g2_z_xzzz = contrBuffer.data(g2off + 45 * j + 39);

                    auto g2_z_yyyy = contrBuffer.data(g2off + 45 * j + 40);

                    auto g2_z_yyyz = contrBuffer.data(g2off + 45 * j + 41);

                    auto g2_z_yyzz = contrBuffer.data(g2off + 45 * j + 42);

                    auto g2_z_yzzz = contrBuffer.data(g2off + 45 * j + 43);

                    auto g2_z_zzzz = contrBuffer.data(g2off + 45 * j + 44);

                    // set up pointers to (SX|g(r,r')|PH)^(m) integrals

                    auto g1_x_xxxxx = contrBuffer.data(g1off + 63 * j);

                    auto g1_x_xxxxy = contrBuffer.data(g1off + 63 * j + 1);

                    auto g1_x_xxxxz = contrBuffer.data(g1off + 63 * j + 2);

                    auto g1_x_xxxyy = contrBuffer.data(g1off + 63 * j + 3);

                    auto g1_x_xxxyz = contrBuffer.data(g1off + 63 * j + 4);

                    auto g1_x_xxxzz = contrBuffer.data(g1off + 63 * j + 5);

                    auto g1_x_xxyyy = contrBuffer.data(g1off + 63 * j + 6);

                    auto g1_x_xxyyz = contrBuffer.data(g1off + 63 * j + 7);

                    auto g1_x_xxyzz = contrBuffer.data(g1off + 63 * j + 8);

                    auto g1_x_xxzzz = contrBuffer.data(g1off + 63 * j + 9);

                    auto g1_x_xyyyy = contrBuffer.data(g1off + 63 * j + 10);

                    auto g1_x_xyyyz = contrBuffer.data(g1off + 63 * j + 11);

                    auto g1_x_xyyzz = contrBuffer.data(g1off + 63 * j + 12);

                    auto g1_x_xyzzz = contrBuffer.data(g1off + 63 * j + 13);

                    auto g1_x_xzzzz = contrBuffer.data(g1off + 63 * j + 14);

                    auto g1_y_xxxxx = contrBuffer.data(g1off + 63 * j + 21);

                    auto g1_y_xxxxy = contrBuffer.data(g1off + 63 * j + 22);

                    auto g1_y_xxxxz = contrBuffer.data(g1off + 63 * j + 23);

                    auto g1_y_xxxyy = contrBuffer.data(g1off + 63 * j + 24);

                    auto g1_y_xxxyz = contrBuffer.data(g1off + 63 * j + 25);

                    auto g1_y_xxxzz = contrBuffer.data(g1off + 63 * j + 26);

                    auto g1_y_xxyyy = contrBuffer.data(g1off + 63 * j + 27);

                    auto g1_y_xxyyz = contrBuffer.data(g1off + 63 * j + 28);

                    auto g1_y_xxyzz = contrBuffer.data(g1off + 63 * j + 29);

                    auto g1_y_xxzzz = contrBuffer.data(g1off + 63 * j + 30);

                    auto g1_y_xyyyy = contrBuffer.data(g1off + 63 * j + 31);

                    auto g1_y_xyyyz = contrBuffer.data(g1off + 63 * j + 32);

                    auto g1_y_xyyzz = contrBuffer.data(g1off + 63 * j + 33);

                    auto g1_y_xyzzz = contrBuffer.data(g1off + 63 * j + 34);

                    auto g1_y_xzzzz = contrBuffer.data(g1off + 63 * j + 35);

                    auto g1_y_yyyyy = contrBuffer.data(g1off + 63 * j + 36);

                    auto g1_y_yyyyz = contrBuffer.data(g1off + 63 * j + 37);

                    auto g1_y_yyyzz = contrBuffer.data(g1off + 63 * j + 38);

                    auto g1_y_yyzzz = contrBuffer.data(g1off + 63 * j + 39);

                    auto g1_y_yzzzz = contrBuffer.data(g1off + 63 * j + 40);

                    auto g1_z_xxxxx = contrBuffer.data(g1off + 63 * j + 42);

                    auto g1_z_xxxxy = contrBuffer.data(g1off + 63 * j + 43);

                    auto g1_z_xxxxz = contrBuffer.data(g1off + 63 * j + 44);

                    auto g1_z_xxxyy = contrBuffer.data(g1off + 63 * j + 45);

                    auto g1_z_xxxyz = contrBuffer.data(g1off + 63 * j + 46);

                    auto g1_z_xxxzz = contrBuffer.data(g1off + 63 * j + 47);

                    auto g1_z_xxyyy = contrBuffer.data(g1off + 63 * j + 48);

                    auto g1_z_xxyyz = contrBuffer.data(g1off + 63 * j + 49);

                    auto g1_z_xxyzz = contrBuffer.data(g1off + 63 * j + 50);

                    auto g1_z_xxzzz = contrBuffer.data(g1off + 63 * j + 51);

                    auto g1_z_xyyyy = contrBuffer.data(g1off + 63 * j + 52);

                    auto g1_z_xyyyz = contrBuffer.data(g1off + 63 * j + 53);

                    auto g1_z_xyyzz = contrBuffer.data(g1off + 63 * j + 54);

                    auto g1_z_xyzzz = contrBuffer.data(g1off + 63 * j + 55);

                    auto g1_z_xzzzz = contrBuffer.data(g1off + 63 * j + 56);

                    auto g1_z_yyyyy = contrBuffer.data(g1off + 63 * j + 57);

                    auto g1_z_yyyyz = contrBuffer.data(g1off + 63 * j + 58);

                    auto g1_z_yyyzz = contrBuffer.data(g1off + 63 * j + 59);

                    auto g1_z_yyzzz = contrBuffer.data(g1off + 63 * j + 60);

                    auto g1_z_yzzzz = contrBuffer.data(g1off + 63 * j + 61);

                    auto g1_z_zzzzz = contrBuffer.data(g1off + 63 * j + 62);

                    // set up pointers to (SX|g(r,r')|DG)^(m) integrals

                    auto g_xx_xxxx = contrBuffer.data(goff + 90 * j);

                    auto g_xx_xxxy = contrBuffer.data(goff + 90 * j + 1);

                    auto g_xx_xxxz = contrBuffer.data(goff + 90 * j + 2);

                    auto g_xx_xxyy = contrBuffer.data(goff + 90 * j + 3);

                    auto g_xx_xxyz = contrBuffer.data(goff + 90 * j + 4);

                    auto g_xx_xxzz = contrBuffer.data(goff + 90 * j + 5);

                    auto g_xx_xyyy = contrBuffer.data(goff + 90 * j + 6);

                    auto g_xx_xyyz = contrBuffer.data(goff + 90 * j + 7);

                    auto g_xx_xyzz = contrBuffer.data(goff + 90 * j + 8);

                    auto g_xx_xzzz = contrBuffer.data(goff + 90 * j + 9);

                    auto g_xx_yyyy = contrBuffer.data(goff + 90 * j + 10);

                    auto g_xx_yyyz = contrBuffer.data(goff + 90 * j + 11);

                    auto g_xx_yyzz = contrBuffer.data(goff + 90 * j + 12);

                    auto g_xx_yzzz = contrBuffer.data(goff + 90 * j + 13);

                    auto g_xx_zzzz = contrBuffer.data(goff + 90 * j + 14);

                    auto g_xy_xxxx = contrBuffer.data(goff + 90 * j + 15);

                    auto g_xy_xxxy = contrBuffer.data(goff + 90 * j + 16);

                    auto g_xy_xxxz = contrBuffer.data(goff + 90 * j + 17);

                    auto g_xy_xxyy = contrBuffer.data(goff + 90 * j + 18);

                    auto g_xy_xxyz = contrBuffer.data(goff + 90 * j + 19);

                    auto g_xy_xxzz = contrBuffer.data(goff + 90 * j + 20);

                    auto g_xy_xyyy = contrBuffer.data(goff + 90 * j + 21);

                    auto g_xy_xyyz = contrBuffer.data(goff + 90 * j + 22);

                    auto g_xy_xyzz = contrBuffer.data(goff + 90 * j + 23);

                    auto g_xy_xzzz = contrBuffer.data(goff + 90 * j + 24);

                    auto g_xy_yyyy = contrBuffer.data(goff + 90 * j + 25);

                    auto g_xy_yyyz = contrBuffer.data(goff + 90 * j + 26);

                    auto g_xy_yyzz = contrBuffer.data(goff + 90 * j + 27);

                    auto g_xy_yzzz = contrBuffer.data(goff + 90 * j + 28);

                    auto g_xy_zzzz = contrBuffer.data(goff + 90 * j + 29);

                    auto g_xz_xxxx = contrBuffer.data(goff + 90 * j + 30);

                    auto g_xz_xxxy = contrBuffer.data(goff + 90 * j + 31);

                    auto g_xz_xxxz = contrBuffer.data(goff + 90 * j + 32);

                    auto g_xz_xxyy = contrBuffer.data(goff + 90 * j + 33);

                    auto g_xz_xxyz = contrBuffer.data(goff + 90 * j + 34);

                    auto g_xz_xxzz = contrBuffer.data(goff + 90 * j + 35);

                    auto g_xz_xyyy = contrBuffer.data(goff + 90 * j + 36);

                    auto g_xz_xyyz = contrBuffer.data(goff + 90 * j + 37);

                    auto g_xz_xyzz = contrBuffer.data(goff + 90 * j + 38);

                    auto g_xz_xzzz = contrBuffer.data(goff + 90 * j + 39);

                    auto g_xz_yyyy = contrBuffer.data(goff + 90 * j + 40);

                    auto g_xz_yyyz = contrBuffer.data(goff + 90 * j + 41);

                    auto g_xz_yyzz = contrBuffer.data(goff + 90 * j + 42);

                    auto g_xz_yzzz = contrBuffer.data(goff + 90 * j + 43);

                    auto g_xz_zzzz = contrBuffer.data(goff + 90 * j + 44);

                    auto g_yy_xxxx = contrBuffer.data(goff + 90 * j + 45);

                    auto g_yy_xxxy = contrBuffer.data(goff + 90 * j + 46);

                    auto g_yy_xxxz = contrBuffer.data(goff + 90 * j + 47);

                    auto g_yy_xxyy = contrBuffer.data(goff + 90 * j + 48);

                    auto g_yy_xxyz = contrBuffer.data(goff + 90 * j + 49);

                    auto g_yy_xxzz = contrBuffer.data(goff + 90 * j + 50);

                    auto g_yy_xyyy = contrBuffer.data(goff + 90 * j + 51);

                    auto g_yy_xyyz = contrBuffer.data(goff + 90 * j + 52);

                    auto g_yy_xyzz = contrBuffer.data(goff + 90 * j + 53);

                    auto g_yy_xzzz = contrBuffer.data(goff + 90 * j + 54);

                    auto g_yy_yyyy = contrBuffer.data(goff + 90 * j + 55);

                    auto g_yy_yyyz = contrBuffer.data(goff + 90 * j + 56);

                    auto g_yy_yyzz = contrBuffer.data(goff + 90 * j + 57);

                    auto g_yy_yzzz = contrBuffer.data(goff + 90 * j + 58);

                    auto g_yy_zzzz = contrBuffer.data(goff + 90 * j + 59);

                    auto g_yz_xxxx = contrBuffer.data(goff + 90 * j + 60);

                    auto g_yz_xxxy = contrBuffer.data(goff + 90 * j + 61);

                    auto g_yz_xxxz = contrBuffer.data(goff + 90 * j + 62);

                    auto g_yz_xxyy = contrBuffer.data(goff + 90 * j + 63);

                    auto g_yz_xxyz = contrBuffer.data(goff + 90 * j + 64);

                    auto g_yz_xxzz = contrBuffer.data(goff + 90 * j + 65);

                    auto g_yz_xyyy = contrBuffer.data(goff + 90 * j + 66);

                    auto g_yz_xyyz = contrBuffer.data(goff + 90 * j + 67);

                    auto g_yz_xyzz = contrBuffer.data(goff + 90 * j + 68);

                    auto g_yz_xzzz = contrBuffer.data(goff + 90 * j + 69);

                    auto g_yz_yyyy = contrBuffer.data(goff + 90 * j + 70);

                    auto g_yz_yyyz = contrBuffer.data(goff + 90 * j + 71);

                    auto g_yz_yyzz = contrBuffer.data(goff + 90 * j + 72);

                    auto g_yz_yzzz = contrBuffer.data(goff + 90 * j + 73);

                    auto g_yz_zzzz = contrBuffer.data(goff + 90 * j + 74);

                    auto g_zz_xxxx = contrBuffer.data(goff + 90 * j + 75);

                    auto g_zz_xxxy = contrBuffer.data(goff + 90 * j + 76);

                    auto g_zz_xxxz = contrBuffer.data(goff + 90 * j + 77);

                    auto g_zz_xxyy = contrBuffer.data(goff + 90 * j + 78);

                    auto g_zz_xxyz = contrBuffer.data(goff + 90 * j + 79);

                    auto g_zz_xxzz = contrBuffer.data(goff + 90 * j + 80);

                    auto g_zz_xyyy = contrBuffer.data(goff + 90 * j + 81);

                    auto g_zz_xyyz = contrBuffer.data(goff + 90 * j + 82);

                    auto g_zz_xyzz = contrBuffer.data(goff + 90 * j + 83);

                    auto g_zz_xzzz = contrBuffer.data(goff + 90 * j + 84);

                    auto g_zz_yyyy = contrBuffer.data(goff + 90 * j + 85);

                    auto g_zz_yyyz = contrBuffer.data(goff + 90 * j + 86);

                    auto g_zz_yyzz = contrBuffer.data(goff + 90 * j + 87);

                    auto g_zz_yzzz = contrBuffer.data(goff + 90 * j + 88);

                    auto g_zz_zzzz = contrBuffer.data(goff + 90 * j + 89);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxx, g2_x_xxxy,\
                                             g2_x_xxxz, g2_x_xxyy, g2_x_xxyz, g2_x_xxzz,\
                                             g2_x_xyyy, g2_x_xyyz, g2_x_xyzz, g2_x_xzzz,\
                                             g2_x_yyyy, g2_x_yyyz, g2_x_yyzz, g2_x_yzzz,\
                                             g2_x_zzzz, g2_y_xxxx, g2_y_xxxy, g2_y_xxxz,\
                                             g2_y_xxyy, g2_y_xxyz, g2_y_xxzz, g2_y_xyyy,\
                                             g2_y_xyyz, g2_y_xyzz, g2_y_xzzz, g2_y_yyyy,\
                                             g2_y_yyyz, g2_y_yyzz, g2_y_yzzz, g2_y_zzzz,\
                                             g2_z_xxxx, g2_z_xxxy, g2_z_xxxz, g2_z_xxyy,\
                                             g2_z_xxyz, g2_z_xxzz, g2_z_xyyy, g2_z_xyyz,\
                                             g2_z_xyzz, g2_z_xzzz, g2_z_yyyy, g2_z_yyyz,\
                                             g2_z_yyzz, g2_z_yzzz, g2_z_zzzz, g1_x_xxxxx,\
                                             g1_x_xxxxy, g1_x_xxxxz, g1_x_xxxyy,\
                                             g1_x_xxxyz, g1_x_xxxzz, g1_x_xxyyy,\
                                             g1_x_xxyyz, g1_x_xxyzz, g1_x_xxzzz,\
                                             g1_x_xyyyy, g1_x_xyyyz, g1_x_xyyzz,\
                                             g1_x_xyzzz, g1_x_xzzzz, g1_y_xxxxx,\
                                             g1_y_xxxxy, g1_y_xxxxz, g1_y_xxxyy,\
                                             g1_y_xxxyz, g1_y_xxxzz, g1_y_xxyyy,\
                                             g1_y_xxyyz, g1_y_xxyzz, g1_y_xxzzz,\
                                             g1_y_xyyyy, g1_y_xyyyz, g1_y_xyyzz,\
                                             g1_y_xyzzz, g1_y_xzzzz, g1_y_yyyyy,\
                                             g1_y_yyyyz, g1_y_yyyzz, g1_y_yyzzz,\
                                             g1_y_yzzzz, g1_z_xxxxx,\
                                             g1_z_xxxxy, g1_z_xxxxz, g1_z_xxxyy,\
                                             g1_z_xxxyz, g1_z_xxxzz, g1_z_xxyyy,\
                                             g1_z_xxyyz, g1_z_xxyzz, g1_z_xxzzz,\
                                             g1_z_xyyyy, g1_z_xyyyz, g1_z_xyyzz,\
                                             g1_z_xyzzz, g1_z_xzzzz, g1_z_yyyyy,\
                                             g1_z_yyyyz, g1_z_yyyzz, g1_z_yyzzz,\
                                             g1_z_yzzzz, g1_z_zzzzz, g_xx_xxxx,\
                                             g_xx_xxxy, g_xx_xxxz, g_xx_xxyy, g_xx_xxyz,\
                                             g_xx_xxzz, g_xx_xyyy, g_xx_xyyz, g_xx_xyzz,\
                                             g_xx_xzzz, g_xx_yyyy, g_xx_yyyz, g_xx_yyzz,\
                                             g_xx_yzzz, g_xx_zzzz, g_xy_xxxx, g_xy_xxxy,\
                                             g_xy_xxxz, g_xy_xxyy, g_xy_xxyz, g_xy_xxzz,\
                                             g_xy_xyyy, g_xy_xyyz, g_xy_xyzz, g_xy_xzzz,\
                                             g_xy_yyyy, g_xy_yyyz, g_xy_yyzz, g_xy_yzzz,\
                                             g_xy_zzzz, g_xz_xxxx, g_xz_xxxy, g_xz_xxxz,\
                                             g_xz_xxyy, g_xz_xxyz, g_xz_xxzz, g_xz_xyyy,\
                                             g_xz_xyyz, g_xz_xyzz, g_xz_xzzz, g_xz_yyyy,\
                                             g_xz_yyyz, g_xz_yyzz, g_xz_yzzz, g_xz_zzzz,\
                                             g_yy_xxxx, g_yy_xxxy, g_yy_xxxz, g_yy_xxyy,\
                                             g_yy_xxyz, g_yy_xxzz, g_yy_xyyy, g_yy_xyyz,\
                                             g_yy_xyzz, g_yy_xzzz, g_yy_yyyy, g_yy_yyyz,\
                                             g_yy_yyzz, g_yy_yzzz, g_yy_zzzz, g_yz_xxxx,\
                                             g_yz_xxxy, g_yz_xxxz, g_yz_xxyy, g_yz_xxyz,\
                                             g_yz_xxzz, g_yz_xyyy, g_yz_xyyz, g_yz_xyzz,\
                                             g_yz_xzzz, g_yz_yyyy, g_yz_yyyz, g_yz_yyzz,\
                                             g_yz_yzzz, g_yz_zzzz, g_zz_xxxx, g_zz_xxxy,\
                                             g_zz_xxxz, g_zz_xxyy, g_zz_xxyz, g_zz_xxzz,\
                                             g_zz_xyyy, g_zz_xyyz, g_zz_xyzz, g_zz_xzzz,\
                                             g_zz_yyyy, g_zz_yyyz, g_zz_yyzz, g_zz_yzzz,\
                                             g_zz_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xx_xxxx[k] = g1_x_xxxxx[k] - fr * g2_x_xxxx[k];

                        g_xx_xxxy[k] = g1_x_xxxxy[k] - fr * g2_x_xxxy[k];

                        g_xx_xxxz[k] = g1_x_xxxxz[k] - fr * g2_x_xxxz[k];

                        g_xx_xxyy[k] = g1_x_xxxyy[k] - fr * g2_x_xxyy[k];

                        g_xx_xxyz[k] = g1_x_xxxyz[k] - fr * g2_x_xxyz[k];

                        g_xx_xxzz[k] = g1_x_xxxzz[k] - fr * g2_x_xxzz[k];

                        g_xx_xyyy[k] = g1_x_xxyyy[k] - fr * g2_x_xyyy[k];

                        g_xx_xyyz[k] = g1_x_xxyyz[k] - fr * g2_x_xyyz[k];

                        g_xx_xyzz[k] = g1_x_xxyzz[k] - fr * g2_x_xyzz[k];

                        g_xx_xzzz[k] = g1_x_xxzzz[k] - fr * g2_x_xzzz[k];

                        g_xx_yyyy[k] = g1_x_xyyyy[k] - fr * g2_x_yyyy[k];

                        g_xx_yyyz[k] = g1_x_xyyyz[k] - fr * g2_x_yyyz[k];

                        g_xx_yyzz[k] = g1_x_xyyzz[k] - fr * g2_x_yyzz[k];

                        g_xx_yzzz[k] = g1_x_xyzzz[k] - fr * g2_x_yzzz[k];

                        g_xx_zzzz[k] = g1_x_xzzzz[k] - fr * g2_x_zzzz[k];

                        g_xy_xxxx[k] = g1_y_xxxxx[k] - fr * g2_y_xxxx[k];

                        g_xy_xxxy[k] = g1_y_xxxxy[k] - fr * g2_y_xxxy[k];

                        g_xy_xxxz[k] = g1_y_xxxxz[k] - fr * g2_y_xxxz[k];

                        g_xy_xxyy[k] = g1_y_xxxyy[k] - fr * g2_y_xxyy[k];

                        g_xy_xxyz[k] = g1_y_xxxyz[k] - fr * g2_y_xxyz[k];

                        g_xy_xxzz[k] = g1_y_xxxzz[k] - fr * g2_y_xxzz[k];

                        g_xy_xyyy[k] = g1_y_xxyyy[k] - fr * g2_y_xyyy[k];

                        g_xy_xyyz[k] = g1_y_xxyyz[k] - fr * g2_y_xyyz[k];

                        g_xy_xyzz[k] = g1_y_xxyzz[k] - fr * g2_y_xyzz[k];

                        g_xy_xzzz[k] = g1_y_xxzzz[k] - fr * g2_y_xzzz[k];

                        g_xy_yyyy[k] = g1_y_xyyyy[k] - fr * g2_y_yyyy[k];

                        g_xy_yyyz[k] = g1_y_xyyyz[k] - fr * g2_y_yyyz[k];

                        g_xy_yyzz[k] = g1_y_xyyzz[k] - fr * g2_y_yyzz[k];

                        g_xy_yzzz[k] = g1_y_xyzzz[k] - fr * g2_y_yzzz[k];

                        g_xy_zzzz[k] = g1_y_xzzzz[k] - fr * g2_y_zzzz[k];

                        g_xz_xxxx[k] = g1_z_xxxxx[k] - fr * g2_z_xxxx[k];

                        g_xz_xxxy[k] = g1_z_xxxxy[k] - fr * g2_z_xxxy[k];

                        g_xz_xxxz[k] = g1_z_xxxxz[k] - fr * g2_z_xxxz[k];

                        g_xz_xxyy[k] = g1_z_xxxyy[k] - fr * g2_z_xxyy[k];

                        g_xz_xxyz[k] = g1_z_xxxyz[k] - fr * g2_z_xxyz[k];

                        g_xz_xxzz[k] = g1_z_xxxzz[k] - fr * g2_z_xxzz[k];

                        g_xz_xyyy[k] = g1_z_xxyyy[k] - fr * g2_z_xyyy[k];

                        g_xz_xyyz[k] = g1_z_xxyyz[k] - fr * g2_z_xyyz[k];

                        g_xz_xyzz[k] = g1_z_xxyzz[k] - fr * g2_z_xyzz[k];

                        g_xz_xzzz[k] = g1_z_xxzzz[k] - fr * g2_z_xzzz[k];

                        g_xz_yyyy[k] = g1_z_xyyyy[k] - fr * g2_z_yyyy[k];

                        g_xz_yyyz[k] = g1_z_xyyyz[k] - fr * g2_z_yyyz[k];

                        g_xz_yyzz[k] = g1_z_xyyzz[k] - fr * g2_z_yyzz[k];

                        g_xz_yzzz[k] = g1_z_xyzzz[k] - fr * g2_z_yzzz[k];

                        g_xz_zzzz[k] = g1_z_xzzzz[k] - fr * g2_z_zzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yy_xxxx[k] = g1_y_xxxxy[k] - fr * g2_y_xxxx[k];

                        g_yy_xxxy[k] = g1_y_xxxyy[k] - fr * g2_y_xxxy[k];

                        g_yy_xxxz[k] = g1_y_xxxyz[k] - fr * g2_y_xxxz[k];

                        g_yy_xxyy[k] = g1_y_xxyyy[k] - fr * g2_y_xxyy[k];

                        g_yy_xxyz[k] = g1_y_xxyyz[k] - fr * g2_y_xxyz[k];

                        g_yy_xxzz[k] = g1_y_xxyzz[k] - fr * g2_y_xxzz[k];

                        g_yy_xyyy[k] = g1_y_xyyyy[k] - fr * g2_y_xyyy[k];

                        g_yy_xyyz[k] = g1_y_xyyyz[k] - fr * g2_y_xyyz[k];

                        g_yy_xyzz[k] = g1_y_xyyzz[k] - fr * g2_y_xyzz[k];

                        g_yy_xzzz[k] = g1_y_xyzzz[k] - fr * g2_y_xzzz[k];

                        g_yy_yyyy[k] = g1_y_yyyyy[k] - fr * g2_y_yyyy[k];

                        g_yy_yyyz[k] = g1_y_yyyyz[k] - fr * g2_y_yyyz[k];

                        g_yy_yyzz[k] = g1_y_yyyzz[k] - fr * g2_y_yyzz[k];

                        g_yy_yzzz[k] = g1_y_yyzzz[k] - fr * g2_y_yzzz[k];

                        g_yy_zzzz[k] = g1_y_yzzzz[k] - fr * g2_y_zzzz[k];

                        g_yz_xxxx[k] = g1_z_xxxxy[k] - fr * g2_z_xxxx[k];

                        g_yz_xxxy[k] = g1_z_xxxyy[k] - fr * g2_z_xxxy[k];

                        g_yz_xxxz[k] = g1_z_xxxyz[k] - fr * g2_z_xxxz[k];

                        g_yz_xxyy[k] = g1_z_xxyyy[k] - fr * g2_z_xxyy[k];

                        g_yz_xxyz[k] = g1_z_xxyyz[k] - fr * g2_z_xxyz[k];

                        g_yz_xxzz[k] = g1_z_xxyzz[k] - fr * g2_z_xxzz[k];

                        g_yz_xyyy[k] = g1_z_xyyyy[k] - fr * g2_z_xyyy[k];

                        g_yz_xyyz[k] = g1_z_xyyyz[k] - fr * g2_z_xyyz[k];

                        g_yz_xyzz[k] = g1_z_xyyzz[k] - fr * g2_z_xyzz[k];

                        g_yz_xzzz[k] = g1_z_xyzzz[k] - fr * g2_z_xzzz[k];

                        g_yz_yyyy[k] = g1_z_yyyyy[k] - fr * g2_z_yyyy[k];

                        g_yz_yyyz[k] = g1_z_yyyyz[k] - fr * g2_z_yyyz[k];

                        g_yz_yyzz[k] = g1_z_yyyzz[k] - fr * g2_z_yyzz[k];

                        g_yz_yzzz[k] = g1_z_yyzzz[k] - fr * g2_z_yzzz[k];

                        g_yz_zzzz[k] = g1_z_yzzzz[k] - fr * g2_z_zzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zz_xxxx[k] = g1_z_xxxxz[k] - fr * g2_z_xxxx[k];

                        g_zz_xxxy[k] = g1_z_xxxyz[k] - fr * g2_z_xxxy[k];

                        g_zz_xxxz[k] = g1_z_xxxzz[k] - fr * g2_z_xxxz[k];

                        g_zz_xxyy[k] = g1_z_xxyyz[k] - fr * g2_z_xxyy[k];

                        g_zz_xxyz[k] = g1_z_xxyzz[k] - fr * g2_z_xxyz[k];

                        g_zz_xxzz[k] = g1_z_xxzzz[k] - fr * g2_z_xxzz[k];

                        g_zz_xyyy[k] = g1_z_xyyyz[k] - fr * g2_z_xyyy[k];

                        g_zz_xyyz[k] = g1_z_xyyzz[k] - fr * g2_z_xyyz[k];

                        g_zz_xyzz[k] = g1_z_xyzzz[k] - fr * g2_z_xyzz[k];

                        g_zz_xzzz[k] = g1_z_xzzzz[k] - fr * g2_z_xzzz[k];

                        g_zz_yyyy[k] = g1_z_yyyyz[k] - fr * g2_z_yyyy[k];

                        g_zz_yyyz[k] = g1_z_yyyzz[k] - fr * g2_z_yyyz[k];

                        g_zz_yyzz[k] = g1_z_yyzzz[k] - fr * g2_z_yyzz[k];

                        g_zz_yzzz[k] = g1_z_yzzzz[k] - fr * g2_z_yzzz[k];

                        g_zz_zzzz[k] = g1_z_zzzzz[k] - fr * g2_z_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXDH(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 2) && (recPattern[i].third() == 5))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|25)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 2, 5});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 6});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 5});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|PH)^(m) integrals

                    auto g2_x_xxxxx = contrBuffer.data(g2off + 63 * j);

                    auto g2_x_xxxxy = contrBuffer.data(g2off + 63 * j + 1);

                    auto g2_x_xxxxz = contrBuffer.data(g2off + 63 * j + 2);

                    auto g2_x_xxxyy = contrBuffer.data(g2off + 63 * j + 3);

                    auto g2_x_xxxyz = contrBuffer.data(g2off + 63 * j + 4);

                    auto g2_x_xxxzz = contrBuffer.data(g2off + 63 * j + 5);

                    auto g2_x_xxyyy = contrBuffer.data(g2off + 63 * j + 6);

                    auto g2_x_xxyyz = contrBuffer.data(g2off + 63 * j + 7);

                    auto g2_x_xxyzz = contrBuffer.data(g2off + 63 * j + 8);

                    auto g2_x_xxzzz = contrBuffer.data(g2off + 63 * j + 9);

                    auto g2_x_xyyyy = contrBuffer.data(g2off + 63 * j + 10);

                    auto g2_x_xyyyz = contrBuffer.data(g2off + 63 * j + 11);

                    auto g2_x_xyyzz = contrBuffer.data(g2off + 63 * j + 12);

                    auto g2_x_xyzzz = contrBuffer.data(g2off + 63 * j + 13);

                    auto g2_x_xzzzz = contrBuffer.data(g2off + 63 * j + 14);

                    auto g2_x_yyyyy = contrBuffer.data(g2off + 63 * j + 15);

                    auto g2_x_yyyyz = contrBuffer.data(g2off + 63 * j + 16);

                    auto g2_x_yyyzz = contrBuffer.data(g2off + 63 * j + 17);

                    auto g2_x_yyzzz = contrBuffer.data(g2off + 63 * j + 18);

                    auto g2_x_yzzzz = contrBuffer.data(g2off + 63 * j + 19);

                    auto g2_x_zzzzz = contrBuffer.data(g2off + 63 * j + 20);

                    auto g2_y_xxxxx = contrBuffer.data(g2off + 63 * j + 21);

                    auto g2_y_xxxxy = contrBuffer.data(g2off + 63 * j + 22);

                    auto g2_y_xxxxz = contrBuffer.data(g2off + 63 * j + 23);

                    auto g2_y_xxxyy = contrBuffer.data(g2off + 63 * j + 24);

                    auto g2_y_xxxyz = contrBuffer.data(g2off + 63 * j + 25);

                    auto g2_y_xxxzz = contrBuffer.data(g2off + 63 * j + 26);

                    auto g2_y_xxyyy = contrBuffer.data(g2off + 63 * j + 27);

                    auto g2_y_xxyyz = contrBuffer.data(g2off + 63 * j + 28);

                    auto g2_y_xxyzz = contrBuffer.data(g2off + 63 * j + 29);

                    auto g2_y_xxzzz = contrBuffer.data(g2off + 63 * j + 30);

                    auto g2_y_xyyyy = contrBuffer.data(g2off + 63 * j + 31);

                    auto g2_y_xyyyz = contrBuffer.data(g2off + 63 * j + 32);

                    auto g2_y_xyyzz = contrBuffer.data(g2off + 63 * j + 33);

                    auto g2_y_xyzzz = contrBuffer.data(g2off + 63 * j + 34);

                    auto g2_y_xzzzz = contrBuffer.data(g2off + 63 * j + 35);

                    auto g2_y_yyyyy = contrBuffer.data(g2off + 63 * j + 36);

                    auto g2_y_yyyyz = contrBuffer.data(g2off + 63 * j + 37);

                    auto g2_y_yyyzz = contrBuffer.data(g2off + 63 * j + 38);

                    auto g2_y_yyzzz = contrBuffer.data(g2off + 63 * j + 39);

                    auto g2_y_yzzzz = contrBuffer.data(g2off + 63 * j + 40);

                    auto g2_y_zzzzz = contrBuffer.data(g2off + 63 * j + 41);

                    auto g2_z_xxxxx = contrBuffer.data(g2off + 63 * j + 42);

                    auto g2_z_xxxxy = contrBuffer.data(g2off + 63 * j + 43);

                    auto g2_z_xxxxz = contrBuffer.data(g2off + 63 * j + 44);

                    auto g2_z_xxxyy = contrBuffer.data(g2off + 63 * j + 45);

                    auto g2_z_xxxyz = contrBuffer.data(g2off + 63 * j + 46);

                    auto g2_z_xxxzz = contrBuffer.data(g2off + 63 * j + 47);

                    auto g2_z_xxyyy = contrBuffer.data(g2off + 63 * j + 48);

                    auto g2_z_xxyyz = contrBuffer.data(g2off + 63 * j + 49);

                    auto g2_z_xxyzz = contrBuffer.data(g2off + 63 * j + 50);

                    auto g2_z_xxzzz = contrBuffer.data(g2off + 63 * j + 51);

                    auto g2_z_xyyyy = contrBuffer.data(g2off + 63 * j + 52);

                    auto g2_z_xyyyz = contrBuffer.data(g2off + 63 * j + 53);

                    auto g2_z_xyyzz = contrBuffer.data(g2off + 63 * j + 54);

                    auto g2_z_xyzzz = contrBuffer.data(g2off + 63 * j + 55);

                    auto g2_z_xzzzz = contrBuffer.data(g2off + 63 * j + 56);

                    auto g2_z_yyyyy = contrBuffer.data(g2off + 63 * j + 57);

                    auto g2_z_yyyyz = contrBuffer.data(g2off + 63 * j + 58);

                    auto g2_z_yyyzz = contrBuffer.data(g2off + 63 * j + 59);

                    auto g2_z_yyzzz = contrBuffer.data(g2off + 63 * j + 60);

                    auto g2_z_yzzzz = contrBuffer.data(g2off + 63 * j + 61);

                    auto g2_z_zzzzz = contrBuffer.data(g2off + 63 * j + 62);

                    // set up pointers to (SX|g(r,r')|PI)^(m) integrals

                    auto g1_x_xxxxxx = contrBuffer.data(g1off + 84 * j);

                    auto g1_x_xxxxxy = contrBuffer.data(g1off + 84 * j + 1);

                    auto g1_x_xxxxxz = contrBuffer.data(g1off + 84 * j + 2);

                    auto g1_x_xxxxyy = contrBuffer.data(g1off + 84 * j + 3);

                    auto g1_x_xxxxyz = contrBuffer.data(g1off + 84 * j + 4);

                    auto g1_x_xxxxzz = contrBuffer.data(g1off + 84 * j + 5);

                    auto g1_x_xxxyyy = contrBuffer.data(g1off + 84 * j + 6);

                    auto g1_x_xxxyyz = contrBuffer.data(g1off + 84 * j + 7);

                    auto g1_x_xxxyzz = contrBuffer.data(g1off + 84 * j + 8);

                    auto g1_x_xxxzzz = contrBuffer.data(g1off + 84 * j + 9);

                    auto g1_x_xxyyyy = contrBuffer.data(g1off + 84 * j + 10);

                    auto g1_x_xxyyyz = contrBuffer.data(g1off + 84 * j + 11);

                    auto g1_x_xxyyzz = contrBuffer.data(g1off + 84 * j + 12);

                    auto g1_x_xxyzzz = contrBuffer.data(g1off + 84 * j + 13);

                    auto g1_x_xxzzzz = contrBuffer.data(g1off + 84 * j + 14);

                    auto g1_x_xyyyyy = contrBuffer.data(g1off + 84 * j + 15);

                    auto g1_x_xyyyyz = contrBuffer.data(g1off + 84 * j + 16);

                    auto g1_x_xyyyzz = contrBuffer.data(g1off + 84 * j + 17);

                    auto g1_x_xyyzzz = contrBuffer.data(g1off + 84 * j + 18);

                    auto g1_x_xyzzzz = contrBuffer.data(g1off + 84 * j + 19);

                    auto g1_x_xzzzzz = contrBuffer.data(g1off + 84 * j + 20);

                    auto g1_y_xxxxxx = contrBuffer.data(g1off + 84 * j + 28);

                    auto g1_y_xxxxxy = contrBuffer.data(g1off + 84 * j + 29);

                    auto g1_y_xxxxxz = contrBuffer.data(g1off + 84 * j + 30);

                    auto g1_y_xxxxyy = contrBuffer.data(g1off + 84 * j + 31);

                    auto g1_y_xxxxyz = contrBuffer.data(g1off + 84 * j + 32);

                    auto g1_y_xxxxzz = contrBuffer.data(g1off + 84 * j + 33);

                    auto g1_y_xxxyyy = contrBuffer.data(g1off + 84 * j + 34);

                    auto g1_y_xxxyyz = contrBuffer.data(g1off + 84 * j + 35);

                    auto g1_y_xxxyzz = contrBuffer.data(g1off + 84 * j + 36);

                    auto g1_y_xxxzzz = contrBuffer.data(g1off + 84 * j + 37);

                    auto g1_y_xxyyyy = contrBuffer.data(g1off + 84 * j + 38);

                    auto g1_y_xxyyyz = contrBuffer.data(g1off + 84 * j + 39);

                    auto g1_y_xxyyzz = contrBuffer.data(g1off + 84 * j + 40);

                    auto g1_y_xxyzzz = contrBuffer.data(g1off + 84 * j + 41);

                    auto g1_y_xxzzzz = contrBuffer.data(g1off + 84 * j + 42);

                    auto g1_y_xyyyyy = contrBuffer.data(g1off + 84 * j + 43);

                    auto g1_y_xyyyyz = contrBuffer.data(g1off + 84 * j + 44);

                    auto g1_y_xyyyzz = contrBuffer.data(g1off + 84 * j + 45);

                    auto g1_y_xyyzzz = contrBuffer.data(g1off + 84 * j + 46);

                    auto g1_y_xyzzzz = contrBuffer.data(g1off + 84 * j + 47);

                    auto g1_y_xzzzzz = contrBuffer.data(g1off + 84 * j + 48);

                    auto g1_y_yyyyyy = contrBuffer.data(g1off + 84 * j + 49);

                    auto g1_y_yyyyyz = contrBuffer.data(g1off + 84 * j + 50);

                    auto g1_y_yyyyzz = contrBuffer.data(g1off + 84 * j + 51);

                    auto g1_y_yyyzzz = contrBuffer.data(g1off + 84 * j + 52);

                    auto g1_y_yyzzzz = contrBuffer.data(g1off + 84 * j + 53);

                    auto g1_y_yzzzzz = contrBuffer.data(g1off + 84 * j + 54);

                    auto g1_z_xxxxxx = contrBuffer.data(g1off + 84 * j + 56);

                    auto g1_z_xxxxxy = contrBuffer.data(g1off + 84 * j + 57);

                    auto g1_z_xxxxxz = contrBuffer.data(g1off + 84 * j + 58);

                    auto g1_z_xxxxyy = contrBuffer.data(g1off + 84 * j + 59);

                    auto g1_z_xxxxyz = contrBuffer.data(g1off + 84 * j + 60);

                    auto g1_z_xxxxzz = contrBuffer.data(g1off + 84 * j + 61);

                    auto g1_z_xxxyyy = contrBuffer.data(g1off + 84 * j + 62);

                    auto g1_z_xxxyyz = contrBuffer.data(g1off + 84 * j + 63);

                    auto g1_z_xxxyzz = contrBuffer.data(g1off + 84 * j + 64);

                    auto g1_z_xxxzzz = contrBuffer.data(g1off + 84 * j + 65);

                    auto g1_z_xxyyyy = contrBuffer.data(g1off + 84 * j + 66);

                    auto g1_z_xxyyyz = contrBuffer.data(g1off + 84 * j + 67);

                    auto g1_z_xxyyzz = contrBuffer.data(g1off + 84 * j + 68);

                    auto g1_z_xxyzzz = contrBuffer.data(g1off + 84 * j + 69);

                    auto g1_z_xxzzzz = contrBuffer.data(g1off + 84 * j + 70);

                    auto g1_z_xyyyyy = contrBuffer.data(g1off + 84 * j + 71);

                    auto g1_z_xyyyyz = contrBuffer.data(g1off + 84 * j + 72);

                    auto g1_z_xyyyzz = contrBuffer.data(g1off + 84 * j + 73);

                    auto g1_z_xyyzzz = contrBuffer.data(g1off + 84 * j + 74);

                    auto g1_z_xyzzzz = contrBuffer.data(g1off + 84 * j + 75);

                    auto g1_z_xzzzzz = contrBuffer.data(g1off + 84 * j + 76);

                    auto g1_z_yyyyyy = contrBuffer.data(g1off + 84 * j + 77);

                    auto g1_z_yyyyyz = contrBuffer.data(g1off + 84 * j + 78);

                    auto g1_z_yyyyzz = contrBuffer.data(g1off + 84 * j + 79);

                    auto g1_z_yyyzzz = contrBuffer.data(g1off + 84 * j + 80);

                    auto g1_z_yyzzzz = contrBuffer.data(g1off + 84 * j + 81);

                    auto g1_z_yzzzzz = contrBuffer.data(g1off + 84 * j + 82);

                    auto g1_z_zzzzzz = contrBuffer.data(g1off + 84 * j + 83);

                    // set up pointers to (SX|g(r,r')|DH)^(m) integrals

                    auto g_xx_xxxxx = contrBuffer.data(goff + 126 * j);

                    auto g_xx_xxxxy = contrBuffer.data(goff + 126 * j + 1);

                    auto g_xx_xxxxz = contrBuffer.data(goff + 126 * j + 2);

                    auto g_xx_xxxyy = contrBuffer.data(goff + 126 * j + 3);

                    auto g_xx_xxxyz = contrBuffer.data(goff + 126 * j + 4);

                    auto g_xx_xxxzz = contrBuffer.data(goff + 126 * j + 5);

                    auto g_xx_xxyyy = contrBuffer.data(goff + 126 * j + 6);

                    auto g_xx_xxyyz = contrBuffer.data(goff + 126 * j + 7);

                    auto g_xx_xxyzz = contrBuffer.data(goff + 126 * j + 8);

                    auto g_xx_xxzzz = contrBuffer.data(goff + 126 * j + 9);

                    auto g_xx_xyyyy = contrBuffer.data(goff + 126 * j + 10);

                    auto g_xx_xyyyz = contrBuffer.data(goff + 126 * j + 11);

                    auto g_xx_xyyzz = contrBuffer.data(goff + 126 * j + 12);

                    auto g_xx_xyzzz = contrBuffer.data(goff + 126 * j + 13);

                    auto g_xx_xzzzz = contrBuffer.data(goff + 126 * j + 14);

                    auto g_xx_yyyyy = contrBuffer.data(goff + 126 * j + 15);

                    auto g_xx_yyyyz = contrBuffer.data(goff + 126 * j + 16);

                    auto g_xx_yyyzz = contrBuffer.data(goff + 126 * j + 17);

                    auto g_xx_yyzzz = contrBuffer.data(goff + 126 * j + 18);

                    auto g_xx_yzzzz = contrBuffer.data(goff + 126 * j + 19);

                    auto g_xx_zzzzz = contrBuffer.data(goff + 126 * j + 20);

                    auto g_xy_xxxxx = contrBuffer.data(goff + 126 * j + 21);

                    auto g_xy_xxxxy = contrBuffer.data(goff + 126 * j + 22);

                    auto g_xy_xxxxz = contrBuffer.data(goff + 126 * j + 23);

                    auto g_xy_xxxyy = contrBuffer.data(goff + 126 * j + 24);

                    auto g_xy_xxxyz = contrBuffer.data(goff + 126 * j + 25);

                    auto g_xy_xxxzz = contrBuffer.data(goff + 126 * j + 26);

                    auto g_xy_xxyyy = contrBuffer.data(goff + 126 * j + 27);

                    auto g_xy_xxyyz = contrBuffer.data(goff + 126 * j + 28);

                    auto g_xy_xxyzz = contrBuffer.data(goff + 126 * j + 29);

                    auto g_xy_xxzzz = contrBuffer.data(goff + 126 * j + 30);

                    auto g_xy_xyyyy = contrBuffer.data(goff + 126 * j + 31);

                    auto g_xy_xyyyz = contrBuffer.data(goff + 126 * j + 32);

                    auto g_xy_xyyzz = contrBuffer.data(goff + 126 * j + 33);

                    auto g_xy_xyzzz = contrBuffer.data(goff + 126 * j + 34);

                    auto g_xy_xzzzz = contrBuffer.data(goff + 126 * j + 35);

                    auto g_xy_yyyyy = contrBuffer.data(goff + 126 * j + 36);

                    auto g_xy_yyyyz = contrBuffer.data(goff + 126 * j + 37);

                    auto g_xy_yyyzz = contrBuffer.data(goff + 126 * j + 38);

                    auto g_xy_yyzzz = contrBuffer.data(goff + 126 * j + 39);

                    auto g_xy_yzzzz = contrBuffer.data(goff + 126 * j + 40);

                    auto g_xy_zzzzz = contrBuffer.data(goff + 126 * j + 41);

                    auto g_xz_xxxxx = contrBuffer.data(goff + 126 * j + 42);

                    auto g_xz_xxxxy = contrBuffer.data(goff + 126 * j + 43);

                    auto g_xz_xxxxz = contrBuffer.data(goff + 126 * j + 44);

                    auto g_xz_xxxyy = contrBuffer.data(goff + 126 * j + 45);

                    auto g_xz_xxxyz = contrBuffer.data(goff + 126 * j + 46);

                    auto g_xz_xxxzz = contrBuffer.data(goff + 126 * j + 47);

                    auto g_xz_xxyyy = contrBuffer.data(goff + 126 * j + 48);

                    auto g_xz_xxyyz = contrBuffer.data(goff + 126 * j + 49);

                    auto g_xz_xxyzz = contrBuffer.data(goff + 126 * j + 50);

                    auto g_xz_xxzzz = contrBuffer.data(goff + 126 * j + 51);

                    auto g_xz_xyyyy = contrBuffer.data(goff + 126 * j + 52);

                    auto g_xz_xyyyz = contrBuffer.data(goff + 126 * j + 53);

                    auto g_xz_xyyzz = contrBuffer.data(goff + 126 * j + 54);

                    auto g_xz_xyzzz = contrBuffer.data(goff + 126 * j + 55);

                    auto g_xz_xzzzz = contrBuffer.data(goff + 126 * j + 56);

                    auto g_xz_yyyyy = contrBuffer.data(goff + 126 * j + 57);

                    auto g_xz_yyyyz = contrBuffer.data(goff + 126 * j + 58);

                    auto g_xz_yyyzz = contrBuffer.data(goff + 126 * j + 59);

                    auto g_xz_yyzzz = contrBuffer.data(goff + 126 * j + 60);

                    auto g_xz_yzzzz = contrBuffer.data(goff + 126 * j + 61);

                    auto g_xz_zzzzz = contrBuffer.data(goff + 126 * j + 62);

                    auto g_yy_xxxxx = contrBuffer.data(goff + 126 * j + 63);

                    auto g_yy_xxxxy = contrBuffer.data(goff + 126 * j + 64);

                    auto g_yy_xxxxz = contrBuffer.data(goff + 126 * j + 65);

                    auto g_yy_xxxyy = contrBuffer.data(goff + 126 * j + 66);

                    auto g_yy_xxxyz = contrBuffer.data(goff + 126 * j + 67);

                    auto g_yy_xxxzz = contrBuffer.data(goff + 126 * j + 68);

                    auto g_yy_xxyyy = contrBuffer.data(goff + 126 * j + 69);

                    auto g_yy_xxyyz = contrBuffer.data(goff + 126 * j + 70);

                    auto g_yy_xxyzz = contrBuffer.data(goff + 126 * j + 71);

                    auto g_yy_xxzzz = contrBuffer.data(goff + 126 * j + 72);

                    auto g_yy_xyyyy = contrBuffer.data(goff + 126 * j + 73);

                    auto g_yy_xyyyz = contrBuffer.data(goff + 126 * j + 74);

                    auto g_yy_xyyzz = contrBuffer.data(goff + 126 * j + 75);

                    auto g_yy_xyzzz = contrBuffer.data(goff + 126 * j + 76);

                    auto g_yy_xzzzz = contrBuffer.data(goff + 126 * j + 77);

                    auto g_yy_yyyyy = contrBuffer.data(goff + 126 * j + 78);

                    auto g_yy_yyyyz = contrBuffer.data(goff + 126 * j + 79);

                    auto g_yy_yyyzz = contrBuffer.data(goff + 126 * j + 80);

                    auto g_yy_yyzzz = contrBuffer.data(goff + 126 * j + 81);

                    auto g_yy_yzzzz = contrBuffer.data(goff + 126 * j + 82);

                    auto g_yy_zzzzz = contrBuffer.data(goff + 126 * j + 83);

                    auto g_yz_xxxxx = contrBuffer.data(goff + 126 * j + 84);

                    auto g_yz_xxxxy = contrBuffer.data(goff + 126 * j + 85);

                    auto g_yz_xxxxz = contrBuffer.data(goff + 126 * j + 86);

                    auto g_yz_xxxyy = contrBuffer.data(goff + 126 * j + 87);

                    auto g_yz_xxxyz = contrBuffer.data(goff + 126 * j + 88);

                    auto g_yz_xxxzz = contrBuffer.data(goff + 126 * j + 89);

                    auto g_yz_xxyyy = contrBuffer.data(goff + 126 * j + 90);

                    auto g_yz_xxyyz = contrBuffer.data(goff + 126 * j + 91);

                    auto g_yz_xxyzz = contrBuffer.data(goff + 126 * j + 92);

                    auto g_yz_xxzzz = contrBuffer.data(goff + 126 * j + 93);

                    auto g_yz_xyyyy = contrBuffer.data(goff + 126 * j + 94);

                    auto g_yz_xyyyz = contrBuffer.data(goff + 126 * j + 95);

                    auto g_yz_xyyzz = contrBuffer.data(goff + 126 * j + 96);

                    auto g_yz_xyzzz = contrBuffer.data(goff + 126 * j + 97);

                    auto g_yz_xzzzz = contrBuffer.data(goff + 126 * j + 98);

                    auto g_yz_yyyyy = contrBuffer.data(goff + 126 * j + 99);

                    auto g_yz_yyyyz = contrBuffer.data(goff + 126 * j + 100);

                    auto g_yz_yyyzz = contrBuffer.data(goff + 126 * j + 101);

                    auto g_yz_yyzzz = contrBuffer.data(goff + 126 * j + 102);

                    auto g_yz_yzzzz = contrBuffer.data(goff + 126 * j + 103);

                    auto g_yz_zzzzz = contrBuffer.data(goff + 126 * j + 104);

                    auto g_zz_xxxxx = contrBuffer.data(goff + 126 * j + 105);

                    auto g_zz_xxxxy = contrBuffer.data(goff + 126 * j + 106);

                    auto g_zz_xxxxz = contrBuffer.data(goff + 126 * j + 107);

                    auto g_zz_xxxyy = contrBuffer.data(goff + 126 * j + 108);

                    auto g_zz_xxxyz = contrBuffer.data(goff + 126 * j + 109);

                    auto g_zz_xxxzz = contrBuffer.data(goff + 126 * j + 110);

                    auto g_zz_xxyyy = contrBuffer.data(goff + 126 * j + 111);

                    auto g_zz_xxyyz = contrBuffer.data(goff + 126 * j + 112);

                    auto g_zz_xxyzz = contrBuffer.data(goff + 126 * j + 113);

                    auto g_zz_xxzzz = contrBuffer.data(goff + 126 * j + 114);

                    auto g_zz_xyyyy = contrBuffer.data(goff + 126 * j + 115);

                    auto g_zz_xyyyz = contrBuffer.data(goff + 126 * j + 116);

                    auto g_zz_xyyzz = contrBuffer.data(goff + 126 * j + 117);

                    auto g_zz_xyzzz = contrBuffer.data(goff + 126 * j + 118);

                    auto g_zz_xzzzz = contrBuffer.data(goff + 126 * j + 119);

                    auto g_zz_yyyyy = contrBuffer.data(goff + 126 * j + 120);

                    auto g_zz_yyyyz = contrBuffer.data(goff + 126 * j + 121);

                    auto g_zz_yyyzz = contrBuffer.data(goff + 126 * j + 122);

                    auto g_zz_yyzzz = contrBuffer.data(goff + 126 * j + 123);

                    auto g_zz_yzzzz = contrBuffer.data(goff + 126 * j + 124);

                    auto g_zz_zzzzz = contrBuffer.data(goff + 126 * j + 125);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxxx, g2_x_xxxxy,\
                                             g2_x_xxxxz, g2_x_xxxyy, g2_x_xxxyz,\
                                             g2_x_xxxzz, g2_x_xxyyy, g2_x_xxyyz,\
                                             g2_x_xxyzz, g2_x_xxzzz, g2_x_xyyyy,\
                                             g2_x_xyyyz, g2_x_xyyzz, g2_x_xyzzz,\
                                             g2_x_xzzzz, g2_x_yyyyy, g2_x_yyyyz,\
                                             g2_x_yyyzz, g2_x_yyzzz, g2_x_yzzzz,\
                                             g2_x_zzzzz, g2_y_xxxxx, g2_y_xxxxy,\
                                             g2_y_xxxxz, g2_y_xxxyy, g2_y_xxxyz,\
                                             g2_y_xxxzz, g2_y_xxyyy, g2_y_xxyyz,\
                                             g2_y_xxyzz, g2_y_xxzzz, g2_y_xyyyy,\
                                             g2_y_xyyyz, g2_y_xyyzz, g2_y_xyzzz,\
                                             g2_y_xzzzz, g2_y_yyyyy, g2_y_yyyyz,\
                                             g2_y_yyyzz, g2_y_yyzzz, g2_y_yzzzz,\
                                             g2_y_zzzzz, g2_z_xxxxx, g2_z_xxxxy,\
                                             g2_z_xxxxz, g2_z_xxxyy, g2_z_xxxyz,\
                                             g2_z_xxxzz, g2_z_xxyyy, g2_z_xxyyz,\
                                             g2_z_xxyzz, g2_z_xxzzz, g2_z_xyyyy,\
                                             g2_z_xyyyz, g2_z_xyyzz, g2_z_xyzzz,\
                                             g2_z_xzzzz, g2_z_yyyyy, g2_z_yyyyz,\
                                             g2_z_yyyzz, g2_z_yyzzz, g2_z_yzzzz,\
                                             g2_z_zzzzz, g1_x_xxxxxx, g1_x_xxxxxy,\
                                             g1_x_xxxxxz, g1_x_xxxxyy, g1_x_xxxxyz,\
                                             g1_x_xxxxzz, g1_x_xxxyyy, g1_x_xxxyyz,\
                                             g1_x_xxxyzz, g1_x_xxxzzz, g1_x_xxyyyy,\
                                             g1_x_xxyyyz, g1_x_xxyyzz, g1_x_xxyzzz,\
                                             g1_x_xxzzzz, g1_x_xyyyyy, g1_x_xyyyyz,\
                                             g1_x_xyyyzz, g1_x_xyyzzz, g1_x_xyzzzz,\
                                             g1_x_xzzzzz, g1_y_xxxxxx,\
                                             g1_y_xxxxxy, g1_y_xxxxxz, g1_y_xxxxyy,\
                                             g1_y_xxxxyz, g1_y_xxxxzz, g1_y_xxxyyy,\
                                             g1_y_xxxyyz, g1_y_xxxyzz, g1_y_xxxzzz,\
                                             g1_y_xxyyyy, g1_y_xxyyyz, g1_y_xxyyzz,\
                                             g1_y_xxyzzz, g1_y_xxzzzz, g1_y_xyyyyy,\
                                             g1_y_xyyyyz, g1_y_xyyyzz, g1_y_xyyzzz,\
                                             g1_y_xyzzzz, g1_y_xzzzzz, g1_y_yyyyyy,\
                                             g1_y_yyyyyz, g1_y_yyyyzz, g1_y_yyyzzz,\
                                             g1_y_yyzzzz, g1_y_yzzzzz,\
                                             g1_z_xxxxxx, g1_z_xxxxxy, g1_z_xxxxxz,\
                                             g1_z_xxxxyy, g1_z_xxxxyz, g1_z_xxxxzz,\
                                             g1_z_xxxyyy, g1_z_xxxyyz, g1_z_xxxyzz,\
                                             g1_z_xxxzzz, g1_z_xxyyyy, g1_z_xxyyyz,\
                                             g1_z_xxyyzz, g1_z_xxyzzz, g1_z_xxzzzz,\
                                             g1_z_xyyyyy, g1_z_xyyyyz, g1_z_xyyyzz,\
                                             g1_z_xyyzzz, g1_z_xyzzzz, g1_z_xzzzzz,\
                                             g1_z_yyyyyy, g1_z_yyyyyz, g1_z_yyyyzz,\
                                             g1_z_yyyzzz, g1_z_yyzzzz, g1_z_yzzzzz,\
                                             g1_z_zzzzzz, g_xx_xxxxx, g_xx_xxxxy,\
                                             g_xx_xxxxz, g_xx_xxxyy, g_xx_xxxyz,\
                                             g_xx_xxxzz, g_xx_xxyyy, g_xx_xxyyz,\
                                             g_xx_xxyzz, g_xx_xxzzz, g_xx_xyyyy,\
                                             g_xx_xyyyz, g_xx_xyyzz, g_xx_xyzzz,\
                                             g_xx_xzzzz, g_xx_yyyyy, g_xx_yyyyz,\
                                             g_xx_yyyzz, g_xx_yyzzz, g_xx_yzzzz,\
                                             g_xx_zzzzz, g_xy_xxxxx, g_xy_xxxxy,\
                                             g_xy_xxxxz, g_xy_xxxyy, g_xy_xxxyz,\
                                             g_xy_xxxzz, g_xy_xxyyy, g_xy_xxyyz,\
                                             g_xy_xxyzz, g_xy_xxzzz, g_xy_xyyyy,\
                                             g_xy_xyyyz, g_xy_xyyzz, g_xy_xyzzz,\
                                             g_xy_xzzzz, g_xy_yyyyy, g_xy_yyyyz,\
                                             g_xy_yyyzz, g_xy_yyzzz, g_xy_yzzzz,\
                                             g_xy_zzzzz, g_xz_xxxxx, g_xz_xxxxy,\
                                             g_xz_xxxxz, g_xz_xxxyy, g_xz_xxxyz,\
                                             g_xz_xxxzz, g_xz_xxyyy, g_xz_xxyyz,\
                                             g_xz_xxyzz, g_xz_xxzzz, g_xz_xyyyy,\
                                             g_xz_xyyyz, g_xz_xyyzz, g_xz_xyzzz,\
                                             g_xz_xzzzz, g_xz_yyyyy, g_xz_yyyyz,\
                                             g_xz_yyyzz, g_xz_yyzzz, g_xz_yzzzz,\
                                             g_xz_zzzzz, g_yy_xxxxx, g_yy_xxxxy,\
                                             g_yy_xxxxz, g_yy_xxxyy, g_yy_xxxyz,\
                                             g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyz,\
                                             g_yy_xxyzz, g_yy_xxzzz, g_yy_xyyyy,\
                                             g_yy_xyyyz, g_yy_xyyzz, g_yy_xyzzz,\
                                             g_yy_xzzzz, g_yy_yyyyy, g_yy_yyyyz,\
                                             g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz,\
                                             g_yy_zzzzz, g_yz_xxxxx, g_yz_xxxxy,\
                                             g_yz_xxxxz, g_yz_xxxyy, g_yz_xxxyz,\
                                             g_yz_xxxzz, g_yz_xxyyy, g_yz_xxyyz,\
                                             g_yz_xxyzz, g_yz_xxzzz, g_yz_xyyyy,\
                                             g_yz_xyyyz, g_yz_xyyzz, g_yz_xyzzz,\
                                             g_yz_xzzzz, g_yz_yyyyy, g_yz_yyyyz,\
                                             g_yz_yyyzz, g_yz_yyzzz, g_yz_yzzzz,\
                                             g_yz_zzzzz, g_zz_xxxxx, g_zz_xxxxy,\
                                             g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyz,\
                                             g_zz_xxxzz, g_zz_xxyyy, g_zz_xxyyz,\
                                             g_zz_xxyzz, g_zz_xxzzz, g_zz_xyyyy,\
                                             g_zz_xyyyz, g_zz_xyyzz, g_zz_xyzzz,\
                                             g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyz,\
                                             g_zz_yyyzz, g_zz_yyzzz, g_zz_yzzzz,\
                                             g_zz_zzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xx_xxxxx[k] = g1_x_xxxxxx[k] - fr * g2_x_xxxxx[k];

                        g_xx_xxxxy[k] = g1_x_xxxxxy[k] - fr * g2_x_xxxxy[k];

                        g_xx_xxxxz[k] = g1_x_xxxxxz[k] - fr * g2_x_xxxxz[k];

                        g_xx_xxxyy[k] = g1_x_xxxxyy[k] - fr * g2_x_xxxyy[k];

                        g_xx_xxxyz[k] = g1_x_xxxxyz[k] - fr * g2_x_xxxyz[k];

                        g_xx_xxxzz[k] = g1_x_xxxxzz[k] - fr * g2_x_xxxzz[k];

                        g_xx_xxyyy[k] = g1_x_xxxyyy[k] - fr * g2_x_xxyyy[k];

                        g_xx_xxyyz[k] = g1_x_xxxyyz[k] - fr * g2_x_xxyyz[k];

                        g_xx_xxyzz[k] = g1_x_xxxyzz[k] - fr * g2_x_xxyzz[k];

                        g_xx_xxzzz[k] = g1_x_xxxzzz[k] - fr * g2_x_xxzzz[k];

                        g_xx_xyyyy[k] = g1_x_xxyyyy[k] - fr * g2_x_xyyyy[k];

                        g_xx_xyyyz[k] = g1_x_xxyyyz[k] - fr * g2_x_xyyyz[k];

                        g_xx_xyyzz[k] = g1_x_xxyyzz[k] - fr * g2_x_xyyzz[k];

                        g_xx_xyzzz[k] = g1_x_xxyzzz[k] - fr * g2_x_xyzzz[k];

                        g_xx_xzzzz[k] = g1_x_xxzzzz[k] - fr * g2_x_xzzzz[k];

                        g_xx_yyyyy[k] = g1_x_xyyyyy[k] - fr * g2_x_yyyyy[k];

                        g_xx_yyyyz[k] = g1_x_xyyyyz[k] - fr * g2_x_yyyyz[k];

                        g_xx_yyyzz[k] = g1_x_xyyyzz[k] - fr * g2_x_yyyzz[k];

                        g_xx_yyzzz[k] = g1_x_xyyzzz[k] - fr * g2_x_yyzzz[k];

                        g_xx_yzzzz[k] = g1_x_xyzzzz[k] - fr * g2_x_yzzzz[k];

                        g_xx_zzzzz[k] = g1_x_xzzzzz[k] - fr * g2_x_zzzzz[k];

                        g_xy_xxxxx[k] = g1_y_xxxxxx[k] - fr * g2_y_xxxxx[k];

                        g_xy_xxxxy[k] = g1_y_xxxxxy[k] - fr * g2_y_xxxxy[k];

                        g_xy_xxxxz[k] = g1_y_xxxxxz[k] - fr * g2_y_xxxxz[k];

                        g_xy_xxxyy[k] = g1_y_xxxxyy[k] - fr * g2_y_xxxyy[k];

                        g_xy_xxxyz[k] = g1_y_xxxxyz[k] - fr * g2_y_xxxyz[k];

                        g_xy_xxxzz[k] = g1_y_xxxxzz[k] - fr * g2_y_xxxzz[k];

                        g_xy_xxyyy[k] = g1_y_xxxyyy[k] - fr * g2_y_xxyyy[k];

                        g_xy_xxyyz[k] = g1_y_xxxyyz[k] - fr * g2_y_xxyyz[k];

                        g_xy_xxyzz[k] = g1_y_xxxyzz[k] - fr * g2_y_xxyzz[k];

                        g_xy_xxzzz[k] = g1_y_xxxzzz[k] - fr * g2_y_xxzzz[k];

                        g_xy_xyyyy[k] = g1_y_xxyyyy[k] - fr * g2_y_xyyyy[k];

                        g_xy_xyyyz[k] = g1_y_xxyyyz[k] - fr * g2_y_xyyyz[k];

                        g_xy_xyyzz[k] = g1_y_xxyyzz[k] - fr * g2_y_xyyzz[k];

                        g_xy_xyzzz[k] = g1_y_xxyzzz[k] - fr * g2_y_xyzzz[k];

                        g_xy_xzzzz[k] = g1_y_xxzzzz[k] - fr * g2_y_xzzzz[k];

                        g_xy_yyyyy[k] = g1_y_xyyyyy[k] - fr * g2_y_yyyyy[k];

                        g_xy_yyyyz[k] = g1_y_xyyyyz[k] - fr * g2_y_yyyyz[k];

                        g_xy_yyyzz[k] = g1_y_xyyyzz[k] - fr * g2_y_yyyzz[k];

                        g_xy_yyzzz[k] = g1_y_xyyzzz[k] - fr * g2_y_yyzzz[k];

                        g_xy_yzzzz[k] = g1_y_xyzzzz[k] - fr * g2_y_yzzzz[k];

                        g_xy_zzzzz[k] = g1_y_xzzzzz[k] - fr * g2_y_zzzzz[k];

                        g_xz_xxxxx[k] = g1_z_xxxxxx[k] - fr * g2_z_xxxxx[k];

                        g_xz_xxxxy[k] = g1_z_xxxxxy[k] - fr * g2_z_xxxxy[k];

                        g_xz_xxxxz[k] = g1_z_xxxxxz[k] - fr * g2_z_xxxxz[k];

                        g_xz_xxxyy[k] = g1_z_xxxxyy[k] - fr * g2_z_xxxyy[k];

                        g_xz_xxxyz[k] = g1_z_xxxxyz[k] - fr * g2_z_xxxyz[k];

                        g_xz_xxxzz[k] = g1_z_xxxxzz[k] - fr * g2_z_xxxzz[k];

                        g_xz_xxyyy[k] = g1_z_xxxyyy[k] - fr * g2_z_xxyyy[k];

                        g_xz_xxyyz[k] = g1_z_xxxyyz[k] - fr * g2_z_xxyyz[k];

                        g_xz_xxyzz[k] = g1_z_xxxyzz[k] - fr * g2_z_xxyzz[k];

                        g_xz_xxzzz[k] = g1_z_xxxzzz[k] - fr * g2_z_xxzzz[k];

                        g_xz_xyyyy[k] = g1_z_xxyyyy[k] - fr * g2_z_xyyyy[k];

                        g_xz_xyyyz[k] = g1_z_xxyyyz[k] - fr * g2_z_xyyyz[k];

                        g_xz_xyyzz[k] = g1_z_xxyyzz[k] - fr * g2_z_xyyzz[k];

                        g_xz_xyzzz[k] = g1_z_xxyzzz[k] - fr * g2_z_xyzzz[k];

                        g_xz_xzzzz[k] = g1_z_xxzzzz[k] - fr * g2_z_xzzzz[k];

                        g_xz_yyyyy[k] = g1_z_xyyyyy[k] - fr * g2_z_yyyyy[k];

                        g_xz_yyyyz[k] = g1_z_xyyyyz[k] - fr * g2_z_yyyyz[k];

                        g_xz_yyyzz[k] = g1_z_xyyyzz[k] - fr * g2_z_yyyzz[k];

                        g_xz_yyzzz[k] = g1_z_xyyzzz[k] - fr * g2_z_yyzzz[k];

                        g_xz_yzzzz[k] = g1_z_xyzzzz[k] - fr * g2_z_yzzzz[k];

                        g_xz_zzzzz[k] = g1_z_xzzzzz[k] - fr * g2_z_zzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yy_xxxxx[k] = g1_y_xxxxxy[k] - fr * g2_y_xxxxx[k];

                        g_yy_xxxxy[k] = g1_y_xxxxyy[k] - fr * g2_y_xxxxy[k];

                        g_yy_xxxxz[k] = g1_y_xxxxyz[k] - fr * g2_y_xxxxz[k];

                        g_yy_xxxyy[k] = g1_y_xxxyyy[k] - fr * g2_y_xxxyy[k];

                        g_yy_xxxyz[k] = g1_y_xxxyyz[k] - fr * g2_y_xxxyz[k];

                        g_yy_xxxzz[k] = g1_y_xxxyzz[k] - fr * g2_y_xxxzz[k];

                        g_yy_xxyyy[k] = g1_y_xxyyyy[k] - fr * g2_y_xxyyy[k];

                        g_yy_xxyyz[k] = g1_y_xxyyyz[k] - fr * g2_y_xxyyz[k];

                        g_yy_xxyzz[k] = g1_y_xxyyzz[k] - fr * g2_y_xxyzz[k];

                        g_yy_xxzzz[k] = g1_y_xxyzzz[k] - fr * g2_y_xxzzz[k];

                        g_yy_xyyyy[k] = g1_y_xyyyyy[k] - fr * g2_y_xyyyy[k];

                        g_yy_xyyyz[k] = g1_y_xyyyyz[k] - fr * g2_y_xyyyz[k];

                        g_yy_xyyzz[k] = g1_y_xyyyzz[k] - fr * g2_y_xyyzz[k];

                        g_yy_xyzzz[k] = g1_y_xyyzzz[k] - fr * g2_y_xyzzz[k];

                        g_yy_xzzzz[k] = g1_y_xyzzzz[k] - fr * g2_y_xzzzz[k];

                        g_yy_yyyyy[k] = g1_y_yyyyyy[k] - fr * g2_y_yyyyy[k];

                        g_yy_yyyyz[k] = g1_y_yyyyyz[k] - fr * g2_y_yyyyz[k];

                        g_yy_yyyzz[k] = g1_y_yyyyzz[k] - fr * g2_y_yyyzz[k];

                        g_yy_yyzzz[k] = g1_y_yyyzzz[k] - fr * g2_y_yyzzz[k];

                        g_yy_yzzzz[k] = g1_y_yyzzzz[k] - fr * g2_y_yzzzz[k];

                        g_yy_zzzzz[k] = g1_y_yzzzzz[k] - fr * g2_y_zzzzz[k];

                        g_yz_xxxxx[k] = g1_z_xxxxxy[k] - fr * g2_z_xxxxx[k];

                        g_yz_xxxxy[k] = g1_z_xxxxyy[k] - fr * g2_z_xxxxy[k];

                        g_yz_xxxxz[k] = g1_z_xxxxyz[k] - fr * g2_z_xxxxz[k];

                        g_yz_xxxyy[k] = g1_z_xxxyyy[k] - fr * g2_z_xxxyy[k];

                        g_yz_xxxyz[k] = g1_z_xxxyyz[k] - fr * g2_z_xxxyz[k];

                        g_yz_xxxzz[k] = g1_z_xxxyzz[k] - fr * g2_z_xxxzz[k];

                        g_yz_xxyyy[k] = g1_z_xxyyyy[k] - fr * g2_z_xxyyy[k];

                        g_yz_xxyyz[k] = g1_z_xxyyyz[k] - fr * g2_z_xxyyz[k];

                        g_yz_xxyzz[k] = g1_z_xxyyzz[k] - fr * g2_z_xxyzz[k];

                        g_yz_xxzzz[k] = g1_z_xxyzzz[k] - fr * g2_z_xxzzz[k];

                        g_yz_xyyyy[k] = g1_z_xyyyyy[k] - fr * g2_z_xyyyy[k];

                        g_yz_xyyyz[k] = g1_z_xyyyyz[k] - fr * g2_z_xyyyz[k];

                        g_yz_xyyzz[k] = g1_z_xyyyzz[k] - fr * g2_z_xyyzz[k];

                        g_yz_xyzzz[k] = g1_z_xyyzzz[k] - fr * g2_z_xyzzz[k];

                        g_yz_xzzzz[k] = g1_z_xyzzzz[k] - fr * g2_z_xzzzz[k];

                        g_yz_yyyyy[k] = g1_z_yyyyyy[k] - fr * g2_z_yyyyy[k];

                        g_yz_yyyyz[k] = g1_z_yyyyyz[k] - fr * g2_z_yyyyz[k];

                        g_yz_yyyzz[k] = g1_z_yyyyzz[k] - fr * g2_z_yyyzz[k];

                        g_yz_yyzzz[k] = g1_z_yyyzzz[k] - fr * g2_z_yyzzz[k];

                        g_yz_yzzzz[k] = g1_z_yyzzzz[k] - fr * g2_z_yzzzz[k];

                        g_yz_zzzzz[k] = g1_z_yzzzzz[k] - fr * g2_z_zzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zz_xxxxx[k] = g1_z_xxxxxz[k] - fr * g2_z_xxxxx[k];

                        g_zz_xxxxy[k] = g1_z_xxxxyz[k] - fr * g2_z_xxxxy[k];

                        g_zz_xxxxz[k] = g1_z_xxxxzz[k] - fr * g2_z_xxxxz[k];

                        g_zz_xxxyy[k] = g1_z_xxxyyz[k] - fr * g2_z_xxxyy[k];

                        g_zz_xxxyz[k] = g1_z_xxxyzz[k] - fr * g2_z_xxxyz[k];

                        g_zz_xxxzz[k] = g1_z_xxxzzz[k] - fr * g2_z_xxxzz[k];

                        g_zz_xxyyy[k] = g1_z_xxyyyz[k] - fr * g2_z_xxyyy[k];

                        g_zz_xxyyz[k] = g1_z_xxyyzz[k] - fr * g2_z_xxyyz[k];

                        g_zz_xxyzz[k] = g1_z_xxyzzz[k] - fr * g2_z_xxyzz[k];

                        g_zz_xxzzz[k] = g1_z_xxzzzz[k] - fr * g2_z_xxzzz[k];

                        g_zz_xyyyy[k] = g1_z_xyyyyz[k] - fr * g2_z_xyyyy[k];

                        g_zz_xyyyz[k] = g1_z_xyyyzz[k] - fr * g2_z_xyyyz[k];

                        g_zz_xyyzz[k] = g1_z_xyyzzz[k] - fr * g2_z_xyyzz[k];

                        g_zz_xyzzz[k] = g1_z_xyzzzz[k] - fr * g2_z_xyzzz[k];

                        g_zz_xzzzz[k] = g1_z_xzzzzz[k] - fr * g2_z_xzzzz[k];

                        g_zz_yyyyy[k] = g1_z_yyyyyz[k] - fr * g2_z_yyyyy[k];

                        g_zz_yyyyz[k] = g1_z_yyyyzz[k] - fr * g2_z_yyyyz[k];

                        g_zz_yyyzz[k] = g1_z_yyyzzz[k] - fr * g2_z_yyyzz[k];

                        g_zz_yyzzz[k] = g1_z_yyzzzz[k] - fr * g2_z_yyzzz[k];

                        g_zz_yzzzz[k] = g1_z_yzzzzz[k] - fr * g2_z_yzzzz[k];

                        g_zz_zzzzz[k] = g1_z_zzzzzz[k] - fr * g2_z_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXDI(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 2) && (recPattern[i].third() == 6))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|26)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 2, 6});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 7});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 1, 6});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|PI)^(m) integrals

                    auto g2_x_xxxxxx = contrBuffer.data(g2off + 84 * j);

                    auto g2_x_xxxxxy = contrBuffer.data(g2off + 84 * j + 1);

                    auto g2_x_xxxxxz = contrBuffer.data(g2off + 84 * j + 2);

                    auto g2_x_xxxxyy = contrBuffer.data(g2off + 84 * j + 3);

                    auto g2_x_xxxxyz = contrBuffer.data(g2off + 84 * j + 4);

                    auto g2_x_xxxxzz = contrBuffer.data(g2off + 84 * j + 5);

                    auto g2_x_xxxyyy = contrBuffer.data(g2off + 84 * j + 6);

                    auto g2_x_xxxyyz = contrBuffer.data(g2off + 84 * j + 7);

                    auto g2_x_xxxyzz = contrBuffer.data(g2off + 84 * j + 8);

                    auto g2_x_xxxzzz = contrBuffer.data(g2off + 84 * j + 9);

                    auto g2_x_xxyyyy = contrBuffer.data(g2off + 84 * j + 10);

                    auto g2_x_xxyyyz = contrBuffer.data(g2off + 84 * j + 11);

                    auto g2_x_xxyyzz = contrBuffer.data(g2off + 84 * j + 12);

                    auto g2_x_xxyzzz = contrBuffer.data(g2off + 84 * j + 13);

                    auto g2_x_xxzzzz = contrBuffer.data(g2off + 84 * j + 14);

                    auto g2_x_xyyyyy = contrBuffer.data(g2off + 84 * j + 15);

                    auto g2_x_xyyyyz = contrBuffer.data(g2off + 84 * j + 16);

                    auto g2_x_xyyyzz = contrBuffer.data(g2off + 84 * j + 17);

                    auto g2_x_xyyzzz = contrBuffer.data(g2off + 84 * j + 18);

                    auto g2_x_xyzzzz = contrBuffer.data(g2off + 84 * j + 19);

                    auto g2_x_xzzzzz = contrBuffer.data(g2off + 84 * j + 20);

                    auto g2_x_yyyyyy = contrBuffer.data(g2off + 84 * j + 21);

                    auto g2_x_yyyyyz = contrBuffer.data(g2off + 84 * j + 22);

                    auto g2_x_yyyyzz = contrBuffer.data(g2off + 84 * j + 23);

                    auto g2_x_yyyzzz = contrBuffer.data(g2off + 84 * j + 24);

                    auto g2_x_yyzzzz = contrBuffer.data(g2off + 84 * j + 25);

                    auto g2_x_yzzzzz = contrBuffer.data(g2off + 84 * j + 26);

                    auto g2_x_zzzzzz = contrBuffer.data(g2off + 84 * j + 27);

                    auto g2_y_xxxxxx = contrBuffer.data(g2off + 84 * j + 28);

                    auto g2_y_xxxxxy = contrBuffer.data(g2off + 84 * j + 29);

                    auto g2_y_xxxxxz = contrBuffer.data(g2off + 84 * j + 30);

                    auto g2_y_xxxxyy = contrBuffer.data(g2off + 84 * j + 31);

                    auto g2_y_xxxxyz = contrBuffer.data(g2off + 84 * j + 32);

                    auto g2_y_xxxxzz = contrBuffer.data(g2off + 84 * j + 33);

                    auto g2_y_xxxyyy = contrBuffer.data(g2off + 84 * j + 34);

                    auto g2_y_xxxyyz = contrBuffer.data(g2off + 84 * j + 35);

                    auto g2_y_xxxyzz = contrBuffer.data(g2off + 84 * j + 36);

                    auto g2_y_xxxzzz = contrBuffer.data(g2off + 84 * j + 37);

                    auto g2_y_xxyyyy = contrBuffer.data(g2off + 84 * j + 38);

                    auto g2_y_xxyyyz = contrBuffer.data(g2off + 84 * j + 39);

                    auto g2_y_xxyyzz = contrBuffer.data(g2off + 84 * j + 40);

                    auto g2_y_xxyzzz = contrBuffer.data(g2off + 84 * j + 41);

                    auto g2_y_xxzzzz = contrBuffer.data(g2off + 84 * j + 42);

                    auto g2_y_xyyyyy = contrBuffer.data(g2off + 84 * j + 43);

                    auto g2_y_xyyyyz = contrBuffer.data(g2off + 84 * j + 44);

                    auto g2_y_xyyyzz = contrBuffer.data(g2off + 84 * j + 45);

                    auto g2_y_xyyzzz = contrBuffer.data(g2off + 84 * j + 46);

                    auto g2_y_xyzzzz = contrBuffer.data(g2off + 84 * j + 47);

                    auto g2_y_xzzzzz = contrBuffer.data(g2off + 84 * j + 48);

                    auto g2_y_yyyyyy = contrBuffer.data(g2off + 84 * j + 49);

                    auto g2_y_yyyyyz = contrBuffer.data(g2off + 84 * j + 50);

                    auto g2_y_yyyyzz = contrBuffer.data(g2off + 84 * j + 51);

                    auto g2_y_yyyzzz = contrBuffer.data(g2off + 84 * j + 52);

                    auto g2_y_yyzzzz = contrBuffer.data(g2off + 84 * j + 53);

                    auto g2_y_yzzzzz = contrBuffer.data(g2off + 84 * j + 54);

                    auto g2_y_zzzzzz = contrBuffer.data(g2off + 84 * j + 55);

                    auto g2_z_xxxxxx = contrBuffer.data(g2off + 84 * j + 56);

                    auto g2_z_xxxxxy = contrBuffer.data(g2off + 84 * j + 57);

                    auto g2_z_xxxxxz = contrBuffer.data(g2off + 84 * j + 58);

                    auto g2_z_xxxxyy = contrBuffer.data(g2off + 84 * j + 59);

                    auto g2_z_xxxxyz = contrBuffer.data(g2off + 84 * j + 60);

                    auto g2_z_xxxxzz = contrBuffer.data(g2off + 84 * j + 61);

                    auto g2_z_xxxyyy = contrBuffer.data(g2off + 84 * j + 62);

                    auto g2_z_xxxyyz = contrBuffer.data(g2off + 84 * j + 63);

                    auto g2_z_xxxyzz = contrBuffer.data(g2off + 84 * j + 64);

                    auto g2_z_xxxzzz = contrBuffer.data(g2off + 84 * j + 65);

                    auto g2_z_xxyyyy = contrBuffer.data(g2off + 84 * j + 66);

                    auto g2_z_xxyyyz = contrBuffer.data(g2off + 84 * j + 67);

                    auto g2_z_xxyyzz = contrBuffer.data(g2off + 84 * j + 68);

                    auto g2_z_xxyzzz = contrBuffer.data(g2off + 84 * j + 69);

                    auto g2_z_xxzzzz = contrBuffer.data(g2off + 84 * j + 70);

                    auto g2_z_xyyyyy = contrBuffer.data(g2off + 84 * j + 71);

                    auto g2_z_xyyyyz = contrBuffer.data(g2off + 84 * j + 72);

                    auto g2_z_xyyyzz = contrBuffer.data(g2off + 84 * j + 73);

                    auto g2_z_xyyzzz = contrBuffer.data(g2off + 84 * j + 74);

                    auto g2_z_xyzzzz = contrBuffer.data(g2off + 84 * j + 75);

                    auto g2_z_xzzzzz = contrBuffer.data(g2off + 84 * j + 76);

                    auto g2_z_yyyyyy = contrBuffer.data(g2off + 84 * j + 77);

                    auto g2_z_yyyyyz = contrBuffer.data(g2off + 84 * j + 78);

                    auto g2_z_yyyyzz = contrBuffer.data(g2off + 84 * j + 79);

                    auto g2_z_yyyzzz = contrBuffer.data(g2off + 84 * j + 80);

                    auto g2_z_yyzzzz = contrBuffer.data(g2off + 84 * j + 81);

                    auto g2_z_yzzzzz = contrBuffer.data(g2off + 84 * j + 82);

                    auto g2_z_zzzzzz = contrBuffer.data(g2off + 84 * j + 83);

                    // set up pointers to (SX|g(r,r')|PK)^(m) integrals

                    auto g1_x_xxxxxxx = contrBuffer.data(g1off + 108 * j);

                    auto g1_x_xxxxxxy = contrBuffer.data(g1off + 108 * j + 1);

                    auto g1_x_xxxxxxz = contrBuffer.data(g1off + 108 * j + 2);

                    auto g1_x_xxxxxyy = contrBuffer.data(g1off + 108 * j + 3);

                    auto g1_x_xxxxxyz = contrBuffer.data(g1off + 108 * j + 4);

                    auto g1_x_xxxxxzz = contrBuffer.data(g1off + 108 * j + 5);

                    auto g1_x_xxxxyyy = contrBuffer.data(g1off + 108 * j + 6);

                    auto g1_x_xxxxyyz = contrBuffer.data(g1off + 108 * j + 7);

                    auto g1_x_xxxxyzz = contrBuffer.data(g1off + 108 * j + 8);

                    auto g1_x_xxxxzzz = contrBuffer.data(g1off + 108 * j + 9);

                    auto g1_x_xxxyyyy = contrBuffer.data(g1off + 108 * j + 10);

                    auto g1_x_xxxyyyz = contrBuffer.data(g1off + 108 * j + 11);

                    auto g1_x_xxxyyzz = contrBuffer.data(g1off + 108 * j + 12);

                    auto g1_x_xxxyzzz = contrBuffer.data(g1off + 108 * j + 13);

                    auto g1_x_xxxzzzz = contrBuffer.data(g1off + 108 * j + 14);

                    auto g1_x_xxyyyyy = contrBuffer.data(g1off + 108 * j + 15);

                    auto g1_x_xxyyyyz = contrBuffer.data(g1off + 108 * j + 16);

                    auto g1_x_xxyyyzz = contrBuffer.data(g1off + 108 * j + 17);

                    auto g1_x_xxyyzzz = contrBuffer.data(g1off + 108 * j + 18);

                    auto g1_x_xxyzzzz = contrBuffer.data(g1off + 108 * j + 19);

                    auto g1_x_xxzzzzz = contrBuffer.data(g1off + 108 * j + 20);

                    auto g1_x_xyyyyyy = contrBuffer.data(g1off + 108 * j + 21);

                    auto g1_x_xyyyyyz = contrBuffer.data(g1off + 108 * j + 22);

                    auto g1_x_xyyyyzz = contrBuffer.data(g1off + 108 * j + 23);

                    auto g1_x_xyyyzzz = contrBuffer.data(g1off + 108 * j + 24);

                    auto g1_x_xyyzzzz = contrBuffer.data(g1off + 108 * j + 25);

                    auto g1_x_xyzzzzz = contrBuffer.data(g1off + 108 * j + 26);

                    auto g1_x_xzzzzzz = contrBuffer.data(g1off + 108 * j + 27);

                    auto g1_y_xxxxxxx = contrBuffer.data(g1off + 108 * j + 36);

                    auto g1_y_xxxxxxy = contrBuffer.data(g1off + 108 * j + 37);

                    auto g1_y_xxxxxxz = contrBuffer.data(g1off + 108 * j + 38);

                    auto g1_y_xxxxxyy = contrBuffer.data(g1off + 108 * j + 39);

                    auto g1_y_xxxxxyz = contrBuffer.data(g1off + 108 * j + 40);

                    auto g1_y_xxxxxzz = contrBuffer.data(g1off + 108 * j + 41);

                    auto g1_y_xxxxyyy = contrBuffer.data(g1off + 108 * j + 42);

                    auto g1_y_xxxxyyz = contrBuffer.data(g1off + 108 * j + 43);

                    auto g1_y_xxxxyzz = contrBuffer.data(g1off + 108 * j + 44);

                    auto g1_y_xxxxzzz = contrBuffer.data(g1off + 108 * j + 45);

                    auto g1_y_xxxyyyy = contrBuffer.data(g1off + 108 * j + 46);

                    auto g1_y_xxxyyyz = contrBuffer.data(g1off + 108 * j + 47);

                    auto g1_y_xxxyyzz = contrBuffer.data(g1off + 108 * j + 48);

                    auto g1_y_xxxyzzz = contrBuffer.data(g1off + 108 * j + 49);

                    auto g1_y_xxxzzzz = contrBuffer.data(g1off + 108 * j + 50);

                    auto g1_y_xxyyyyy = contrBuffer.data(g1off + 108 * j + 51);

                    auto g1_y_xxyyyyz = contrBuffer.data(g1off + 108 * j + 52);

                    auto g1_y_xxyyyzz = contrBuffer.data(g1off + 108 * j + 53);

                    auto g1_y_xxyyzzz = contrBuffer.data(g1off + 108 * j + 54);

                    auto g1_y_xxyzzzz = contrBuffer.data(g1off + 108 * j + 55);

                    auto g1_y_xxzzzzz = contrBuffer.data(g1off + 108 * j + 56);

                    auto g1_y_xyyyyyy = contrBuffer.data(g1off + 108 * j + 57);

                    auto g1_y_xyyyyyz = contrBuffer.data(g1off + 108 * j + 58);

                    auto g1_y_xyyyyzz = contrBuffer.data(g1off + 108 * j + 59);

                    auto g1_y_xyyyzzz = contrBuffer.data(g1off + 108 * j + 60);

                    auto g1_y_xyyzzzz = contrBuffer.data(g1off + 108 * j + 61);

                    auto g1_y_xyzzzzz = contrBuffer.data(g1off + 108 * j + 62);

                    auto g1_y_xzzzzzz = contrBuffer.data(g1off + 108 * j + 63);

                    auto g1_y_yyyyyyy = contrBuffer.data(g1off + 108 * j + 64);

                    auto g1_y_yyyyyyz = contrBuffer.data(g1off + 108 * j + 65);

                    auto g1_y_yyyyyzz = contrBuffer.data(g1off + 108 * j + 66);

                    auto g1_y_yyyyzzz = contrBuffer.data(g1off + 108 * j + 67);

                    auto g1_y_yyyzzzz = contrBuffer.data(g1off + 108 * j + 68);

                    auto g1_y_yyzzzzz = contrBuffer.data(g1off + 108 * j + 69);

                    auto g1_y_yzzzzzz = contrBuffer.data(g1off + 108 * j + 70);

                    auto g1_z_xxxxxxx = contrBuffer.data(g1off + 108 * j + 72);

                    auto g1_z_xxxxxxy = contrBuffer.data(g1off + 108 * j + 73);

                    auto g1_z_xxxxxxz = contrBuffer.data(g1off + 108 * j + 74);

                    auto g1_z_xxxxxyy = contrBuffer.data(g1off + 108 * j + 75);

                    auto g1_z_xxxxxyz = contrBuffer.data(g1off + 108 * j + 76);

                    auto g1_z_xxxxxzz = contrBuffer.data(g1off + 108 * j + 77);

                    auto g1_z_xxxxyyy = contrBuffer.data(g1off + 108 * j + 78);

                    auto g1_z_xxxxyyz = contrBuffer.data(g1off + 108 * j + 79);

                    auto g1_z_xxxxyzz = contrBuffer.data(g1off + 108 * j + 80);

                    auto g1_z_xxxxzzz = contrBuffer.data(g1off + 108 * j + 81);

                    auto g1_z_xxxyyyy = contrBuffer.data(g1off + 108 * j + 82);

                    auto g1_z_xxxyyyz = contrBuffer.data(g1off + 108 * j + 83);

                    auto g1_z_xxxyyzz = contrBuffer.data(g1off + 108 * j + 84);

                    auto g1_z_xxxyzzz = contrBuffer.data(g1off + 108 * j + 85);

                    auto g1_z_xxxzzzz = contrBuffer.data(g1off + 108 * j + 86);

                    auto g1_z_xxyyyyy = contrBuffer.data(g1off + 108 * j + 87);

                    auto g1_z_xxyyyyz = contrBuffer.data(g1off + 108 * j + 88);

                    auto g1_z_xxyyyzz = contrBuffer.data(g1off + 108 * j + 89);

                    auto g1_z_xxyyzzz = contrBuffer.data(g1off + 108 * j + 90);

                    auto g1_z_xxyzzzz = contrBuffer.data(g1off + 108 * j + 91);

                    auto g1_z_xxzzzzz = contrBuffer.data(g1off + 108 * j + 92);

                    auto g1_z_xyyyyyy = contrBuffer.data(g1off + 108 * j + 93);

                    auto g1_z_xyyyyyz = contrBuffer.data(g1off + 108 * j + 94);

                    auto g1_z_xyyyyzz = contrBuffer.data(g1off + 108 * j + 95);

                    auto g1_z_xyyyzzz = contrBuffer.data(g1off + 108 * j + 96);

                    auto g1_z_xyyzzzz = contrBuffer.data(g1off + 108 * j + 97);

                    auto g1_z_xyzzzzz = contrBuffer.data(g1off + 108 * j + 98);

                    auto g1_z_xzzzzzz = contrBuffer.data(g1off + 108 * j + 99);

                    auto g1_z_yyyyyyy = contrBuffer.data(g1off + 108 * j + 100);

                    auto g1_z_yyyyyyz = contrBuffer.data(g1off + 108 * j + 101);

                    auto g1_z_yyyyyzz = contrBuffer.data(g1off + 108 * j + 102);

                    auto g1_z_yyyyzzz = contrBuffer.data(g1off + 108 * j + 103);

                    auto g1_z_yyyzzzz = contrBuffer.data(g1off + 108 * j + 104);

                    auto g1_z_yyzzzzz = contrBuffer.data(g1off + 108 * j + 105);

                    auto g1_z_yzzzzzz = contrBuffer.data(g1off + 108 * j + 106);

                    auto g1_z_zzzzzzz = contrBuffer.data(g1off + 108 * j + 107);

                    // set up pointers to (SX|g(r,r')|DI)^(m) integrals

                    auto g_xx_xxxxxx = contrBuffer.data(goff + 168 * j);

                    auto g_xx_xxxxxy = contrBuffer.data(goff + 168 * j + 1);

                    auto g_xx_xxxxxz = contrBuffer.data(goff + 168 * j + 2);

                    auto g_xx_xxxxyy = contrBuffer.data(goff + 168 * j + 3);

                    auto g_xx_xxxxyz = contrBuffer.data(goff + 168 * j + 4);

                    auto g_xx_xxxxzz = contrBuffer.data(goff + 168 * j + 5);

                    auto g_xx_xxxyyy = contrBuffer.data(goff + 168 * j + 6);

                    auto g_xx_xxxyyz = contrBuffer.data(goff + 168 * j + 7);

                    auto g_xx_xxxyzz = contrBuffer.data(goff + 168 * j + 8);

                    auto g_xx_xxxzzz = contrBuffer.data(goff + 168 * j + 9);

                    auto g_xx_xxyyyy = contrBuffer.data(goff + 168 * j + 10);

                    auto g_xx_xxyyyz = contrBuffer.data(goff + 168 * j + 11);

                    auto g_xx_xxyyzz = contrBuffer.data(goff + 168 * j + 12);

                    auto g_xx_xxyzzz = contrBuffer.data(goff + 168 * j + 13);

                    auto g_xx_xxzzzz = contrBuffer.data(goff + 168 * j + 14);

                    auto g_xx_xyyyyy = contrBuffer.data(goff + 168 * j + 15);

                    auto g_xx_xyyyyz = contrBuffer.data(goff + 168 * j + 16);

                    auto g_xx_xyyyzz = contrBuffer.data(goff + 168 * j + 17);

                    auto g_xx_xyyzzz = contrBuffer.data(goff + 168 * j + 18);

                    auto g_xx_xyzzzz = contrBuffer.data(goff + 168 * j + 19);

                    auto g_xx_xzzzzz = contrBuffer.data(goff + 168 * j + 20);

                    auto g_xx_yyyyyy = contrBuffer.data(goff + 168 * j + 21);

                    auto g_xx_yyyyyz = contrBuffer.data(goff + 168 * j + 22);

                    auto g_xx_yyyyzz = contrBuffer.data(goff + 168 * j + 23);

                    auto g_xx_yyyzzz = contrBuffer.data(goff + 168 * j + 24);

                    auto g_xx_yyzzzz = contrBuffer.data(goff + 168 * j + 25);

                    auto g_xx_yzzzzz = contrBuffer.data(goff + 168 * j + 26);

                    auto g_xx_zzzzzz = contrBuffer.data(goff + 168 * j + 27);

                    auto g_xy_xxxxxx = contrBuffer.data(goff + 168 * j + 28);

                    auto g_xy_xxxxxy = contrBuffer.data(goff + 168 * j + 29);

                    auto g_xy_xxxxxz = contrBuffer.data(goff + 168 * j + 30);

                    auto g_xy_xxxxyy = contrBuffer.data(goff + 168 * j + 31);

                    auto g_xy_xxxxyz = contrBuffer.data(goff + 168 * j + 32);

                    auto g_xy_xxxxzz = contrBuffer.data(goff + 168 * j + 33);

                    auto g_xy_xxxyyy = contrBuffer.data(goff + 168 * j + 34);

                    auto g_xy_xxxyyz = contrBuffer.data(goff + 168 * j + 35);

                    auto g_xy_xxxyzz = contrBuffer.data(goff + 168 * j + 36);

                    auto g_xy_xxxzzz = contrBuffer.data(goff + 168 * j + 37);

                    auto g_xy_xxyyyy = contrBuffer.data(goff + 168 * j + 38);

                    auto g_xy_xxyyyz = contrBuffer.data(goff + 168 * j + 39);

                    auto g_xy_xxyyzz = contrBuffer.data(goff + 168 * j + 40);

                    auto g_xy_xxyzzz = contrBuffer.data(goff + 168 * j + 41);

                    auto g_xy_xxzzzz = contrBuffer.data(goff + 168 * j + 42);

                    auto g_xy_xyyyyy = contrBuffer.data(goff + 168 * j + 43);

                    auto g_xy_xyyyyz = contrBuffer.data(goff + 168 * j + 44);

                    auto g_xy_xyyyzz = contrBuffer.data(goff + 168 * j + 45);

                    auto g_xy_xyyzzz = contrBuffer.data(goff + 168 * j + 46);

                    auto g_xy_xyzzzz = contrBuffer.data(goff + 168 * j + 47);

                    auto g_xy_xzzzzz = contrBuffer.data(goff + 168 * j + 48);

                    auto g_xy_yyyyyy = contrBuffer.data(goff + 168 * j + 49);

                    auto g_xy_yyyyyz = contrBuffer.data(goff + 168 * j + 50);

                    auto g_xy_yyyyzz = contrBuffer.data(goff + 168 * j + 51);

                    auto g_xy_yyyzzz = contrBuffer.data(goff + 168 * j + 52);

                    auto g_xy_yyzzzz = contrBuffer.data(goff + 168 * j + 53);

                    auto g_xy_yzzzzz = contrBuffer.data(goff + 168 * j + 54);

                    auto g_xy_zzzzzz = contrBuffer.data(goff + 168 * j + 55);

                    auto g_xz_xxxxxx = contrBuffer.data(goff + 168 * j + 56);

                    auto g_xz_xxxxxy = contrBuffer.data(goff + 168 * j + 57);

                    auto g_xz_xxxxxz = contrBuffer.data(goff + 168 * j + 58);

                    auto g_xz_xxxxyy = contrBuffer.data(goff + 168 * j + 59);

                    auto g_xz_xxxxyz = contrBuffer.data(goff + 168 * j + 60);

                    auto g_xz_xxxxzz = contrBuffer.data(goff + 168 * j + 61);

                    auto g_xz_xxxyyy = contrBuffer.data(goff + 168 * j + 62);

                    auto g_xz_xxxyyz = contrBuffer.data(goff + 168 * j + 63);

                    auto g_xz_xxxyzz = contrBuffer.data(goff + 168 * j + 64);

                    auto g_xz_xxxzzz = contrBuffer.data(goff + 168 * j + 65);

                    auto g_xz_xxyyyy = contrBuffer.data(goff + 168 * j + 66);

                    auto g_xz_xxyyyz = contrBuffer.data(goff + 168 * j + 67);

                    auto g_xz_xxyyzz = contrBuffer.data(goff + 168 * j + 68);

                    auto g_xz_xxyzzz = contrBuffer.data(goff + 168 * j + 69);

                    auto g_xz_xxzzzz = contrBuffer.data(goff + 168 * j + 70);

                    auto g_xz_xyyyyy = contrBuffer.data(goff + 168 * j + 71);

                    auto g_xz_xyyyyz = contrBuffer.data(goff + 168 * j + 72);

                    auto g_xz_xyyyzz = contrBuffer.data(goff + 168 * j + 73);

                    auto g_xz_xyyzzz = contrBuffer.data(goff + 168 * j + 74);

                    auto g_xz_xyzzzz = contrBuffer.data(goff + 168 * j + 75);

                    auto g_xz_xzzzzz = contrBuffer.data(goff + 168 * j + 76);

                    auto g_xz_yyyyyy = contrBuffer.data(goff + 168 * j + 77);

                    auto g_xz_yyyyyz = contrBuffer.data(goff + 168 * j + 78);

                    auto g_xz_yyyyzz = contrBuffer.data(goff + 168 * j + 79);

                    auto g_xz_yyyzzz = contrBuffer.data(goff + 168 * j + 80);

                    auto g_xz_yyzzzz = contrBuffer.data(goff + 168 * j + 81);

                    auto g_xz_yzzzzz = contrBuffer.data(goff + 168 * j + 82);

                    auto g_xz_zzzzzz = contrBuffer.data(goff + 168 * j + 83);

                    auto g_yy_xxxxxx = contrBuffer.data(goff + 168 * j + 84);

                    auto g_yy_xxxxxy = contrBuffer.data(goff + 168 * j + 85);

                    auto g_yy_xxxxxz = contrBuffer.data(goff + 168 * j + 86);

                    auto g_yy_xxxxyy = contrBuffer.data(goff + 168 * j + 87);

                    auto g_yy_xxxxyz = contrBuffer.data(goff + 168 * j + 88);

                    auto g_yy_xxxxzz = contrBuffer.data(goff + 168 * j + 89);

                    auto g_yy_xxxyyy = contrBuffer.data(goff + 168 * j + 90);

                    auto g_yy_xxxyyz = contrBuffer.data(goff + 168 * j + 91);

                    auto g_yy_xxxyzz = contrBuffer.data(goff + 168 * j + 92);

                    auto g_yy_xxxzzz = contrBuffer.data(goff + 168 * j + 93);

                    auto g_yy_xxyyyy = contrBuffer.data(goff + 168 * j + 94);

                    auto g_yy_xxyyyz = contrBuffer.data(goff + 168 * j + 95);

                    auto g_yy_xxyyzz = contrBuffer.data(goff + 168 * j + 96);

                    auto g_yy_xxyzzz = contrBuffer.data(goff + 168 * j + 97);

                    auto g_yy_xxzzzz = contrBuffer.data(goff + 168 * j + 98);

                    auto g_yy_xyyyyy = contrBuffer.data(goff + 168 * j + 99);

                    auto g_yy_xyyyyz = contrBuffer.data(goff + 168 * j + 100);

                    auto g_yy_xyyyzz = contrBuffer.data(goff + 168 * j + 101);

                    auto g_yy_xyyzzz = contrBuffer.data(goff + 168 * j + 102);

                    auto g_yy_xyzzzz = contrBuffer.data(goff + 168 * j + 103);

                    auto g_yy_xzzzzz = contrBuffer.data(goff + 168 * j + 104);

                    auto g_yy_yyyyyy = contrBuffer.data(goff + 168 * j + 105);

                    auto g_yy_yyyyyz = contrBuffer.data(goff + 168 * j + 106);

                    auto g_yy_yyyyzz = contrBuffer.data(goff + 168 * j + 107);

                    auto g_yy_yyyzzz = contrBuffer.data(goff + 168 * j + 108);

                    auto g_yy_yyzzzz = contrBuffer.data(goff + 168 * j + 109);

                    auto g_yy_yzzzzz = contrBuffer.data(goff + 168 * j + 110);

                    auto g_yy_zzzzzz = contrBuffer.data(goff + 168 * j + 111);

                    auto g_yz_xxxxxx = contrBuffer.data(goff + 168 * j + 112);

                    auto g_yz_xxxxxy = contrBuffer.data(goff + 168 * j + 113);

                    auto g_yz_xxxxxz = contrBuffer.data(goff + 168 * j + 114);

                    auto g_yz_xxxxyy = contrBuffer.data(goff + 168 * j + 115);

                    auto g_yz_xxxxyz = contrBuffer.data(goff + 168 * j + 116);

                    auto g_yz_xxxxzz = contrBuffer.data(goff + 168 * j + 117);

                    auto g_yz_xxxyyy = contrBuffer.data(goff + 168 * j + 118);

                    auto g_yz_xxxyyz = contrBuffer.data(goff + 168 * j + 119);

                    auto g_yz_xxxyzz = contrBuffer.data(goff + 168 * j + 120);

                    auto g_yz_xxxzzz = contrBuffer.data(goff + 168 * j + 121);

                    auto g_yz_xxyyyy = contrBuffer.data(goff + 168 * j + 122);

                    auto g_yz_xxyyyz = contrBuffer.data(goff + 168 * j + 123);

                    auto g_yz_xxyyzz = contrBuffer.data(goff + 168 * j + 124);

                    auto g_yz_xxyzzz = contrBuffer.data(goff + 168 * j + 125);

                    auto g_yz_xxzzzz = contrBuffer.data(goff + 168 * j + 126);

                    auto g_yz_xyyyyy = contrBuffer.data(goff + 168 * j + 127);

                    auto g_yz_xyyyyz = contrBuffer.data(goff + 168 * j + 128);

                    auto g_yz_xyyyzz = contrBuffer.data(goff + 168 * j + 129);

                    auto g_yz_xyyzzz = contrBuffer.data(goff + 168 * j + 130);

                    auto g_yz_xyzzzz = contrBuffer.data(goff + 168 * j + 131);

                    auto g_yz_xzzzzz = contrBuffer.data(goff + 168 * j + 132);

                    auto g_yz_yyyyyy = contrBuffer.data(goff + 168 * j + 133);

                    auto g_yz_yyyyyz = contrBuffer.data(goff + 168 * j + 134);

                    auto g_yz_yyyyzz = contrBuffer.data(goff + 168 * j + 135);

                    auto g_yz_yyyzzz = contrBuffer.data(goff + 168 * j + 136);

                    auto g_yz_yyzzzz = contrBuffer.data(goff + 168 * j + 137);

                    auto g_yz_yzzzzz = contrBuffer.data(goff + 168 * j + 138);

                    auto g_yz_zzzzzz = contrBuffer.data(goff + 168 * j + 139);

                    auto g_zz_xxxxxx = contrBuffer.data(goff + 168 * j + 140);

                    auto g_zz_xxxxxy = contrBuffer.data(goff + 168 * j + 141);

                    auto g_zz_xxxxxz = contrBuffer.data(goff + 168 * j + 142);

                    auto g_zz_xxxxyy = contrBuffer.data(goff + 168 * j + 143);

                    auto g_zz_xxxxyz = contrBuffer.data(goff + 168 * j + 144);

                    auto g_zz_xxxxzz = contrBuffer.data(goff + 168 * j + 145);

                    auto g_zz_xxxyyy = contrBuffer.data(goff + 168 * j + 146);

                    auto g_zz_xxxyyz = contrBuffer.data(goff + 168 * j + 147);

                    auto g_zz_xxxyzz = contrBuffer.data(goff + 168 * j + 148);

                    auto g_zz_xxxzzz = contrBuffer.data(goff + 168 * j + 149);

                    auto g_zz_xxyyyy = contrBuffer.data(goff + 168 * j + 150);

                    auto g_zz_xxyyyz = contrBuffer.data(goff + 168 * j + 151);

                    auto g_zz_xxyyzz = contrBuffer.data(goff + 168 * j + 152);

                    auto g_zz_xxyzzz = contrBuffer.data(goff + 168 * j + 153);

                    auto g_zz_xxzzzz = contrBuffer.data(goff + 168 * j + 154);

                    auto g_zz_xyyyyy = contrBuffer.data(goff + 168 * j + 155);

                    auto g_zz_xyyyyz = contrBuffer.data(goff + 168 * j + 156);

                    auto g_zz_xyyyzz = contrBuffer.data(goff + 168 * j + 157);

                    auto g_zz_xyyzzz = contrBuffer.data(goff + 168 * j + 158);

                    auto g_zz_xyzzzz = contrBuffer.data(goff + 168 * j + 159);

                    auto g_zz_xzzzzz = contrBuffer.data(goff + 168 * j + 160);

                    auto g_zz_yyyyyy = contrBuffer.data(goff + 168 * j + 161);

                    auto g_zz_yyyyyz = contrBuffer.data(goff + 168 * j + 162);

                    auto g_zz_yyyyzz = contrBuffer.data(goff + 168 * j + 163);

                    auto g_zz_yyyzzz = contrBuffer.data(goff + 168 * j + 164);

                    auto g_zz_yyzzzz = contrBuffer.data(goff + 168 * j + 165);

                    auto g_zz_yzzzzz = contrBuffer.data(goff + 168 * j + 166);

                    auto g_zz_zzzzzz = contrBuffer.data(goff + 168 * j + 167);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_x_xxxxxx, g2_x_xxxxxy,\
                                             g2_x_xxxxxz, g2_x_xxxxyy, g2_x_xxxxyz,\
                                             g2_x_xxxxzz, g2_x_xxxyyy, g2_x_xxxyyz,\
                                             g2_x_xxxyzz, g2_x_xxxzzz, g2_x_xxyyyy,\
                                             g2_x_xxyyyz, g2_x_xxyyzz, g2_x_xxyzzz,\
                                             g2_x_xxzzzz, g2_x_xyyyyy, g2_x_xyyyyz,\
                                             g2_x_xyyyzz, g2_x_xyyzzz, g2_x_xyzzzz,\
                                             g2_x_xzzzzz, g2_x_yyyyyy, g2_x_yyyyyz,\
                                             g2_x_yyyyzz, g2_x_yyyzzz, g2_x_yyzzzz,\
                                             g2_x_yzzzzz, g2_x_zzzzzz, g2_y_xxxxxx,\
                                             g2_y_xxxxxy, g2_y_xxxxxz, g2_y_xxxxyy,\
                                             g2_y_xxxxyz, g2_y_xxxxzz, g2_y_xxxyyy,\
                                             g2_y_xxxyyz, g2_y_xxxyzz, g2_y_xxxzzz,\
                                             g2_y_xxyyyy, g2_y_xxyyyz, g2_y_xxyyzz,\
                                             g2_y_xxyzzz, g2_y_xxzzzz, g2_y_xyyyyy,\
                                             g2_y_xyyyyz, g2_y_xyyyzz, g2_y_xyyzzz,\
                                             g2_y_xyzzzz, g2_y_xzzzzz, g2_y_yyyyyy,\
                                             g2_y_yyyyyz, g2_y_yyyyzz, g2_y_yyyzzz,\
                                             g2_y_yyzzzz, g2_y_yzzzzz, g2_y_zzzzzz,\
                                             g2_z_xxxxxx, g2_z_xxxxxy, g2_z_xxxxxz,\
                                             g2_z_xxxxyy, g2_z_xxxxyz, g2_z_xxxxzz,\
                                             g2_z_xxxyyy, g2_z_xxxyyz, g2_z_xxxyzz,\
                                             g2_z_xxxzzz, g2_z_xxyyyy, g2_z_xxyyyz,\
                                             g2_z_xxyyzz, g2_z_xxyzzz, g2_z_xxzzzz,\
                                             g2_z_xyyyyy, g2_z_xyyyyz, g2_z_xyyyzz,\
                                             g2_z_xyyzzz, g2_z_xyzzzz, g2_z_xzzzzz,\
                                             g2_z_yyyyyy, g2_z_yyyyyz, g2_z_yyyyzz,\
                                             g2_z_yyyzzz, g2_z_yyzzzz, g2_z_yzzzzz,\
                                             g2_z_zzzzzz, g1_x_xxxxxxx, g1_x_xxxxxxy,\
                                             g1_x_xxxxxxz, g1_x_xxxxxyy, g1_x_xxxxxyz,\
                                             g1_x_xxxxxzz, g1_x_xxxxyyy, g1_x_xxxxyyz,\
                                             g1_x_xxxxyzz, g1_x_xxxxzzz, g1_x_xxxyyyy,\
                                             g1_x_xxxyyyz, g1_x_xxxyyzz, g1_x_xxxyzzz,\
                                             g1_x_xxxzzzz, g1_x_xxyyyyy, g1_x_xxyyyyz,\
                                             g1_x_xxyyyzz, g1_x_xxyyzzz, g1_x_xxyzzzz,\
                                             g1_x_xxzzzzz, g1_x_xyyyyyy, g1_x_xyyyyyz,\
                                             g1_x_xyyyyzz, g1_x_xyyyzzz, g1_x_xyyzzzz,\
                                             g1_x_xyzzzzz, g1_x_xzzzzzz, g1_y_xxxxxxx,\
                                             g1_y_xxxxxxy,\
                                             g1_y_xxxxxxz, g1_y_xxxxxyy, g1_y_xxxxxyz,\
                                             g1_y_xxxxxzz, g1_y_xxxxyyy, g1_y_xxxxyyz,\
                                             g1_y_xxxxyzz, g1_y_xxxxzzz, g1_y_xxxyyyy,\
                                             g1_y_xxxyyyz, g1_y_xxxyyzz, g1_y_xxxyzzz,\
                                             g1_y_xxxzzzz, g1_y_xxyyyyy, g1_y_xxyyyyz,\
                                             g1_y_xxyyyzz, g1_y_xxyyzzz, g1_y_xxyzzzz,\
                                             g1_y_xxzzzzz, g1_y_xyyyyyy, g1_y_xyyyyyz,\
                                             g1_y_xyyyyzz, g1_y_xyyyzzz, g1_y_xyyzzzz,\
                                             g1_y_xyzzzzz, g1_y_xzzzzzz, g1_y_yyyyyyy,\
                                             g1_y_yyyyyyz, g1_y_yyyyyzz, g1_y_yyyyzzz,\
                                             g1_y_yyyzzzz, g1_y_yyzzzzz, g1_y_yzzzzzz,\
                                             g1_z_xxxxxxx, g1_z_xxxxxxy,\
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
                                             g1_z_zzzzzzz, g_xx_xxxxxx, g_xx_xxxxxy,\
                                             g_xx_xxxxxz, g_xx_xxxxyy, g_xx_xxxxyz,\
                                             g_xx_xxxxzz, g_xx_xxxyyy, g_xx_xxxyyz,\
                                             g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy,\
                                             g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz,\
                                             g_xx_xxzzzz, g_xx_xyyyyy, g_xx_xyyyyz,\
                                             g_xx_xyyyzz, g_xx_xyyzzz, g_xx_xyzzzz,\
                                             g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz,\
                                             g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz,\
                                             g_xx_yzzzzz, g_xx_zzzzzz, g_xy_xxxxxx,\
                                             g_xy_xxxxxy, g_xy_xxxxxz, g_xy_xxxxyy,\
                                             g_xy_xxxxyz, g_xy_xxxxzz, g_xy_xxxyyy,\
                                             g_xy_xxxyyz, g_xy_xxxyzz, g_xy_xxxzzz,\
                                             g_xy_xxyyyy, g_xy_xxyyyz, g_xy_xxyyzz,\
                                             g_xy_xxyzzz, g_xy_xxzzzz, g_xy_xyyyyy,\
                                             g_xy_xyyyyz, g_xy_xyyyzz, g_xy_xyyzzz,\
                                             g_xy_xyzzzz, g_xy_xzzzzz, g_xy_yyyyyy,\
                                             g_xy_yyyyyz, g_xy_yyyyzz, g_xy_yyyzzz,\
                                             g_xy_yyzzzz, g_xy_yzzzzz, g_xy_zzzzzz,\
                                             g_xz_xxxxxx, g_xz_xxxxxy, g_xz_xxxxxz,\
                                             g_xz_xxxxyy, g_xz_xxxxyz, g_xz_xxxxzz,\
                                             g_xz_xxxyyy, g_xz_xxxyyz, g_xz_xxxyzz,\
                                             g_xz_xxxzzz, g_xz_xxyyyy, g_xz_xxyyyz,\
                                             g_xz_xxyyzz, g_xz_xxyzzz, g_xz_xxzzzz,\
                                             g_xz_xyyyyy, g_xz_xyyyyz, g_xz_xyyyzz,\
                                             g_xz_xyyzzz, g_xz_xyzzzz, g_xz_xzzzzz,\
                                             g_xz_yyyyyy, g_xz_yyyyyz, g_xz_yyyyzz,\
                                             g_xz_yyyzzz, g_xz_yyzzzz, g_xz_yzzzzz,\
                                             g_xz_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy,\
                                             g_yy_xxxxxz, g_yy_xxxxyy, g_yy_xxxxyz,\
                                             g_yy_xxxxzz, g_yy_xxxyyy, g_yy_xxxyyz,\
                                             g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy,\
                                             g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz,\
                                             g_yy_xxzzzz, g_yy_xyyyyy, g_yy_xyyyyz,\
                                             g_yy_xyyyzz, g_yy_xyyzzz, g_yy_xyzzzz,\
                                             g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz,\
                                             g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz,\
                                             g_yy_yzzzzz, g_yy_zzzzzz, g_yz_xxxxxx,\
                                             g_yz_xxxxxy, g_yz_xxxxxz, g_yz_xxxxyy,\
                                             g_yz_xxxxyz, g_yz_xxxxzz, g_yz_xxxyyy,\
                                             g_yz_xxxyyz, g_yz_xxxyzz, g_yz_xxxzzz,\
                                             g_yz_xxyyyy, g_yz_xxyyyz, g_yz_xxyyzz,\
                                             g_yz_xxyzzz, g_yz_xxzzzz, g_yz_xyyyyy,\
                                             g_yz_xyyyyz, g_yz_xyyyzz, g_yz_xyyzzz,\
                                             g_yz_xyzzzz, g_yz_xzzzzz, g_yz_yyyyyy,\
                                             g_yz_yyyyyz, g_yz_yyyyzz, g_yz_yyyzzz,\
                                             g_yz_yyzzzz, g_yz_yzzzzz, g_yz_zzzzzz,\
                                             g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz,\
                                             g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz,\
                                             g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyzz,\
                                             g_zz_xxxzzz, g_zz_xxyyyy, g_zz_xxyyyz,\
                                             g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz,\
                                             g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz,\
                                             g_zz_xyyzzz, g_zz_xyzzzz, g_zz_xzzzzz,\
                                             g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyzz,\
                                             g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz,\
                                             g_zz_zzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xx_xxxxxx[k] = g1_x_xxxxxxx[k] - fr * g2_x_xxxxxx[k];

                        g_xx_xxxxxy[k] = g1_x_xxxxxxy[k] - fr * g2_x_xxxxxy[k];

                        g_xx_xxxxxz[k] = g1_x_xxxxxxz[k] - fr * g2_x_xxxxxz[k];

                        g_xx_xxxxyy[k] = g1_x_xxxxxyy[k] - fr * g2_x_xxxxyy[k];

                        g_xx_xxxxyz[k] = g1_x_xxxxxyz[k] - fr * g2_x_xxxxyz[k];

                        g_xx_xxxxzz[k] = g1_x_xxxxxzz[k] - fr * g2_x_xxxxzz[k];

                        g_xx_xxxyyy[k] = g1_x_xxxxyyy[k] - fr * g2_x_xxxyyy[k];

                        g_xx_xxxyyz[k] = g1_x_xxxxyyz[k] - fr * g2_x_xxxyyz[k];

                        g_xx_xxxyzz[k] = g1_x_xxxxyzz[k] - fr * g2_x_xxxyzz[k];

                        g_xx_xxxzzz[k] = g1_x_xxxxzzz[k] - fr * g2_x_xxxzzz[k];

                        g_xx_xxyyyy[k] = g1_x_xxxyyyy[k] - fr * g2_x_xxyyyy[k];

                        g_xx_xxyyyz[k] = g1_x_xxxyyyz[k] - fr * g2_x_xxyyyz[k];

                        g_xx_xxyyzz[k] = g1_x_xxxyyzz[k] - fr * g2_x_xxyyzz[k];

                        g_xx_xxyzzz[k] = g1_x_xxxyzzz[k] - fr * g2_x_xxyzzz[k];

                        g_xx_xxzzzz[k] = g1_x_xxxzzzz[k] - fr * g2_x_xxzzzz[k];

                        g_xx_xyyyyy[k] = g1_x_xxyyyyy[k] - fr * g2_x_xyyyyy[k];

                        g_xx_xyyyyz[k] = g1_x_xxyyyyz[k] - fr * g2_x_xyyyyz[k];

                        g_xx_xyyyzz[k] = g1_x_xxyyyzz[k] - fr * g2_x_xyyyzz[k];

                        g_xx_xyyzzz[k] = g1_x_xxyyzzz[k] - fr * g2_x_xyyzzz[k];

                        g_xx_xyzzzz[k] = g1_x_xxyzzzz[k] - fr * g2_x_xyzzzz[k];

                        g_xx_xzzzzz[k] = g1_x_xxzzzzz[k] - fr * g2_x_xzzzzz[k];

                        g_xx_yyyyyy[k] = g1_x_xyyyyyy[k] - fr * g2_x_yyyyyy[k];

                        g_xx_yyyyyz[k] = g1_x_xyyyyyz[k] - fr * g2_x_yyyyyz[k];

                        g_xx_yyyyzz[k] = g1_x_xyyyyzz[k] - fr * g2_x_yyyyzz[k];

                        g_xx_yyyzzz[k] = g1_x_xyyyzzz[k] - fr * g2_x_yyyzzz[k];

                        g_xx_yyzzzz[k] = g1_x_xyyzzzz[k] - fr * g2_x_yyzzzz[k];

                        g_xx_yzzzzz[k] = g1_x_xyzzzzz[k] - fr * g2_x_yzzzzz[k];

                        g_xx_zzzzzz[k] = g1_x_xzzzzzz[k] - fr * g2_x_zzzzzz[k];

                        g_xy_xxxxxx[k] = g1_y_xxxxxxx[k] - fr * g2_y_xxxxxx[k];

                        g_xy_xxxxxy[k] = g1_y_xxxxxxy[k] - fr * g2_y_xxxxxy[k];

                        g_xy_xxxxxz[k] = g1_y_xxxxxxz[k] - fr * g2_y_xxxxxz[k];

                        g_xy_xxxxyy[k] = g1_y_xxxxxyy[k] - fr * g2_y_xxxxyy[k];

                        g_xy_xxxxyz[k] = g1_y_xxxxxyz[k] - fr * g2_y_xxxxyz[k];

                        g_xy_xxxxzz[k] = g1_y_xxxxxzz[k] - fr * g2_y_xxxxzz[k];

                        g_xy_xxxyyy[k] = g1_y_xxxxyyy[k] - fr * g2_y_xxxyyy[k];

                        g_xy_xxxyyz[k] = g1_y_xxxxyyz[k] - fr * g2_y_xxxyyz[k];

                        g_xy_xxxyzz[k] = g1_y_xxxxyzz[k] - fr * g2_y_xxxyzz[k];

                        g_xy_xxxzzz[k] = g1_y_xxxxzzz[k] - fr * g2_y_xxxzzz[k];

                        g_xy_xxyyyy[k] = g1_y_xxxyyyy[k] - fr * g2_y_xxyyyy[k];

                        g_xy_xxyyyz[k] = g1_y_xxxyyyz[k] - fr * g2_y_xxyyyz[k];

                        g_xy_xxyyzz[k] = g1_y_xxxyyzz[k] - fr * g2_y_xxyyzz[k];

                        g_xy_xxyzzz[k] = g1_y_xxxyzzz[k] - fr * g2_y_xxyzzz[k];

                        g_xy_xxzzzz[k] = g1_y_xxxzzzz[k] - fr * g2_y_xxzzzz[k];

                        g_xy_xyyyyy[k] = g1_y_xxyyyyy[k] - fr * g2_y_xyyyyy[k];

                        g_xy_xyyyyz[k] = g1_y_xxyyyyz[k] - fr * g2_y_xyyyyz[k];

                        g_xy_xyyyzz[k] = g1_y_xxyyyzz[k] - fr * g2_y_xyyyzz[k];

                        g_xy_xyyzzz[k] = g1_y_xxyyzzz[k] - fr * g2_y_xyyzzz[k];

                        g_xy_xyzzzz[k] = g1_y_xxyzzzz[k] - fr * g2_y_xyzzzz[k];

                        g_xy_xzzzzz[k] = g1_y_xxzzzzz[k] - fr * g2_y_xzzzzz[k];

                        g_xy_yyyyyy[k] = g1_y_xyyyyyy[k] - fr * g2_y_yyyyyy[k];

                        g_xy_yyyyyz[k] = g1_y_xyyyyyz[k] - fr * g2_y_yyyyyz[k];

                        g_xy_yyyyzz[k] = g1_y_xyyyyzz[k] - fr * g2_y_yyyyzz[k];

                        g_xy_yyyzzz[k] = g1_y_xyyyzzz[k] - fr * g2_y_yyyzzz[k];

                        g_xy_yyzzzz[k] = g1_y_xyyzzzz[k] - fr * g2_y_yyzzzz[k];

                        g_xy_yzzzzz[k] = g1_y_xyzzzzz[k] - fr * g2_y_yzzzzz[k];

                        g_xy_zzzzzz[k] = g1_y_xzzzzzz[k] - fr * g2_y_zzzzzz[k];

                        g_xz_xxxxxx[k] = g1_z_xxxxxxx[k] - fr * g2_z_xxxxxx[k];

                        g_xz_xxxxxy[k] = g1_z_xxxxxxy[k] - fr * g2_z_xxxxxy[k];

                        g_xz_xxxxxz[k] = g1_z_xxxxxxz[k] - fr * g2_z_xxxxxz[k];

                        g_xz_xxxxyy[k] = g1_z_xxxxxyy[k] - fr * g2_z_xxxxyy[k];

                        g_xz_xxxxyz[k] = g1_z_xxxxxyz[k] - fr * g2_z_xxxxyz[k];

                        g_xz_xxxxzz[k] = g1_z_xxxxxzz[k] - fr * g2_z_xxxxzz[k];

                        g_xz_xxxyyy[k] = g1_z_xxxxyyy[k] - fr * g2_z_xxxyyy[k];

                        g_xz_xxxyyz[k] = g1_z_xxxxyyz[k] - fr * g2_z_xxxyyz[k];

                        g_xz_xxxyzz[k] = g1_z_xxxxyzz[k] - fr * g2_z_xxxyzz[k];

                        g_xz_xxxzzz[k] = g1_z_xxxxzzz[k] - fr * g2_z_xxxzzz[k];

                        g_xz_xxyyyy[k] = g1_z_xxxyyyy[k] - fr * g2_z_xxyyyy[k];

                        g_xz_xxyyyz[k] = g1_z_xxxyyyz[k] - fr * g2_z_xxyyyz[k];

                        g_xz_xxyyzz[k] = g1_z_xxxyyzz[k] - fr * g2_z_xxyyzz[k];

                        g_xz_xxyzzz[k] = g1_z_xxxyzzz[k] - fr * g2_z_xxyzzz[k];

                        g_xz_xxzzzz[k] = g1_z_xxxzzzz[k] - fr * g2_z_xxzzzz[k];

                        g_xz_xyyyyy[k] = g1_z_xxyyyyy[k] - fr * g2_z_xyyyyy[k];

                        g_xz_xyyyyz[k] = g1_z_xxyyyyz[k] - fr * g2_z_xyyyyz[k];

                        g_xz_xyyyzz[k] = g1_z_xxyyyzz[k] - fr * g2_z_xyyyzz[k];

                        g_xz_xyyzzz[k] = g1_z_xxyyzzz[k] - fr * g2_z_xyyzzz[k];

                        g_xz_xyzzzz[k] = g1_z_xxyzzzz[k] - fr * g2_z_xyzzzz[k];

                        g_xz_xzzzzz[k] = g1_z_xxzzzzz[k] - fr * g2_z_xzzzzz[k];

                        g_xz_yyyyyy[k] = g1_z_xyyyyyy[k] - fr * g2_z_yyyyyy[k];

                        g_xz_yyyyyz[k] = g1_z_xyyyyyz[k] - fr * g2_z_yyyyyz[k];

                        g_xz_yyyyzz[k] = g1_z_xyyyyzz[k] - fr * g2_z_yyyyzz[k];

                        g_xz_yyyzzz[k] = g1_z_xyyyzzz[k] - fr * g2_z_yyyzzz[k];

                        g_xz_yyzzzz[k] = g1_z_xyyzzzz[k] - fr * g2_z_yyzzzz[k];

                        g_xz_yzzzzz[k] = g1_z_xyzzzzz[k] - fr * g2_z_yzzzzz[k];

                        g_xz_zzzzzz[k] = g1_z_xzzzzzz[k] - fr * g2_z_zzzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yy_xxxxxx[k] = g1_y_xxxxxxy[k] - fr * g2_y_xxxxxx[k];

                        g_yy_xxxxxy[k] = g1_y_xxxxxyy[k] - fr * g2_y_xxxxxy[k];

                        g_yy_xxxxxz[k] = g1_y_xxxxxyz[k] - fr * g2_y_xxxxxz[k];

                        g_yy_xxxxyy[k] = g1_y_xxxxyyy[k] - fr * g2_y_xxxxyy[k];

                        g_yy_xxxxyz[k] = g1_y_xxxxyyz[k] - fr * g2_y_xxxxyz[k];

                        g_yy_xxxxzz[k] = g1_y_xxxxyzz[k] - fr * g2_y_xxxxzz[k];

                        g_yy_xxxyyy[k] = g1_y_xxxyyyy[k] - fr * g2_y_xxxyyy[k];

                        g_yy_xxxyyz[k] = g1_y_xxxyyyz[k] - fr * g2_y_xxxyyz[k];

                        g_yy_xxxyzz[k] = g1_y_xxxyyzz[k] - fr * g2_y_xxxyzz[k];

                        g_yy_xxxzzz[k] = g1_y_xxxyzzz[k] - fr * g2_y_xxxzzz[k];

                        g_yy_xxyyyy[k] = g1_y_xxyyyyy[k] - fr * g2_y_xxyyyy[k];

                        g_yy_xxyyyz[k] = g1_y_xxyyyyz[k] - fr * g2_y_xxyyyz[k];

                        g_yy_xxyyzz[k] = g1_y_xxyyyzz[k] - fr * g2_y_xxyyzz[k];

                        g_yy_xxyzzz[k] = g1_y_xxyyzzz[k] - fr * g2_y_xxyzzz[k];

                        g_yy_xxzzzz[k] = g1_y_xxyzzzz[k] - fr * g2_y_xxzzzz[k];

                        g_yy_xyyyyy[k] = g1_y_xyyyyyy[k] - fr * g2_y_xyyyyy[k];

                        g_yy_xyyyyz[k] = g1_y_xyyyyyz[k] - fr * g2_y_xyyyyz[k];

                        g_yy_xyyyzz[k] = g1_y_xyyyyzz[k] - fr * g2_y_xyyyzz[k];

                        g_yy_xyyzzz[k] = g1_y_xyyyzzz[k] - fr * g2_y_xyyzzz[k];

                        g_yy_xyzzzz[k] = g1_y_xyyzzzz[k] - fr * g2_y_xyzzzz[k];

                        g_yy_xzzzzz[k] = g1_y_xyzzzzz[k] - fr * g2_y_xzzzzz[k];

                        g_yy_yyyyyy[k] = g1_y_yyyyyyy[k] - fr * g2_y_yyyyyy[k];

                        g_yy_yyyyyz[k] = g1_y_yyyyyyz[k] - fr * g2_y_yyyyyz[k];

                        g_yy_yyyyzz[k] = g1_y_yyyyyzz[k] - fr * g2_y_yyyyzz[k];

                        g_yy_yyyzzz[k] = g1_y_yyyyzzz[k] - fr * g2_y_yyyzzz[k];

                        g_yy_yyzzzz[k] = g1_y_yyyzzzz[k] - fr * g2_y_yyzzzz[k];

                        g_yy_yzzzzz[k] = g1_y_yyzzzzz[k] - fr * g2_y_yzzzzz[k];

                        g_yy_zzzzzz[k] = g1_y_yzzzzzz[k] - fr * g2_y_zzzzzz[k];

                        g_yz_xxxxxx[k] = g1_z_xxxxxxy[k] - fr * g2_z_xxxxxx[k];

                        g_yz_xxxxxy[k] = g1_z_xxxxxyy[k] - fr * g2_z_xxxxxy[k];

                        g_yz_xxxxxz[k] = g1_z_xxxxxyz[k] - fr * g2_z_xxxxxz[k];

                        g_yz_xxxxyy[k] = g1_z_xxxxyyy[k] - fr * g2_z_xxxxyy[k];

                        g_yz_xxxxyz[k] = g1_z_xxxxyyz[k] - fr * g2_z_xxxxyz[k];

                        g_yz_xxxxzz[k] = g1_z_xxxxyzz[k] - fr * g2_z_xxxxzz[k];

                        g_yz_xxxyyy[k] = g1_z_xxxyyyy[k] - fr * g2_z_xxxyyy[k];

                        g_yz_xxxyyz[k] = g1_z_xxxyyyz[k] - fr * g2_z_xxxyyz[k];

                        g_yz_xxxyzz[k] = g1_z_xxxyyzz[k] - fr * g2_z_xxxyzz[k];

                        g_yz_xxxzzz[k] = g1_z_xxxyzzz[k] - fr * g2_z_xxxzzz[k];

                        g_yz_xxyyyy[k] = g1_z_xxyyyyy[k] - fr * g2_z_xxyyyy[k];

                        g_yz_xxyyyz[k] = g1_z_xxyyyyz[k] - fr * g2_z_xxyyyz[k];

                        g_yz_xxyyzz[k] = g1_z_xxyyyzz[k] - fr * g2_z_xxyyzz[k];

                        g_yz_xxyzzz[k] = g1_z_xxyyzzz[k] - fr * g2_z_xxyzzz[k];

                        g_yz_xxzzzz[k] = g1_z_xxyzzzz[k] - fr * g2_z_xxzzzz[k];

                        g_yz_xyyyyy[k] = g1_z_xyyyyyy[k] - fr * g2_z_xyyyyy[k];

                        g_yz_xyyyyz[k] = g1_z_xyyyyyz[k] - fr * g2_z_xyyyyz[k];

                        g_yz_xyyyzz[k] = g1_z_xyyyyzz[k] - fr * g2_z_xyyyzz[k];

                        g_yz_xyyzzz[k] = g1_z_xyyyzzz[k] - fr * g2_z_xyyzzz[k];

                        g_yz_xyzzzz[k] = g1_z_xyyzzzz[k] - fr * g2_z_xyzzzz[k];

                        g_yz_xzzzzz[k] = g1_z_xyzzzzz[k] - fr * g2_z_xzzzzz[k];

                        g_yz_yyyyyy[k] = g1_z_yyyyyyy[k] - fr * g2_z_yyyyyy[k];

                        g_yz_yyyyyz[k] = g1_z_yyyyyyz[k] - fr * g2_z_yyyyyz[k];

                        g_yz_yyyyzz[k] = g1_z_yyyyyzz[k] - fr * g2_z_yyyyzz[k];

                        g_yz_yyyzzz[k] = g1_z_yyyyzzz[k] - fr * g2_z_yyyzzz[k];

                        g_yz_yyzzzz[k] = g1_z_yyyzzzz[k] - fr * g2_z_yyzzzz[k];

                        g_yz_yzzzzz[k] = g1_z_yyzzzzz[k] - fr * g2_z_yzzzzz[k];

                        g_yz_zzzzzz[k] = g1_z_yzzzzzz[k] - fr * g2_z_zzzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zz_xxxxxx[k] = g1_z_xxxxxxz[k] - fr * g2_z_xxxxxx[k];

                        g_zz_xxxxxy[k] = g1_z_xxxxxyz[k] - fr * g2_z_xxxxxy[k];

                        g_zz_xxxxxz[k] = g1_z_xxxxxzz[k] - fr * g2_z_xxxxxz[k];

                        g_zz_xxxxyy[k] = g1_z_xxxxyyz[k] - fr * g2_z_xxxxyy[k];

                        g_zz_xxxxyz[k] = g1_z_xxxxyzz[k] - fr * g2_z_xxxxyz[k];

                        g_zz_xxxxzz[k] = g1_z_xxxxzzz[k] - fr * g2_z_xxxxzz[k];

                        g_zz_xxxyyy[k] = g1_z_xxxyyyz[k] - fr * g2_z_xxxyyy[k];

                        g_zz_xxxyyz[k] = g1_z_xxxyyzz[k] - fr * g2_z_xxxyyz[k];

                        g_zz_xxxyzz[k] = g1_z_xxxyzzz[k] - fr * g2_z_xxxyzz[k];

                        g_zz_xxxzzz[k] = g1_z_xxxzzzz[k] - fr * g2_z_xxxzzz[k];

                        g_zz_xxyyyy[k] = g1_z_xxyyyyz[k] - fr * g2_z_xxyyyy[k];

                        g_zz_xxyyyz[k] = g1_z_xxyyyzz[k] - fr * g2_z_xxyyyz[k];

                        g_zz_xxyyzz[k] = g1_z_xxyyzzz[k] - fr * g2_z_xxyyzz[k];

                        g_zz_xxyzzz[k] = g1_z_xxyzzzz[k] - fr * g2_z_xxyzzz[k];

                        g_zz_xxzzzz[k] = g1_z_xxzzzzz[k] - fr * g2_z_xxzzzz[k];

                        g_zz_xyyyyy[k] = g1_z_xyyyyyz[k] - fr * g2_z_xyyyyy[k];

                        g_zz_xyyyyz[k] = g1_z_xyyyyzz[k] - fr * g2_z_xyyyyz[k];

                        g_zz_xyyyzz[k] = g1_z_xyyyzzz[k] - fr * g2_z_xyyyzz[k];

                        g_zz_xyyzzz[k] = g1_z_xyyzzzz[k] - fr * g2_z_xyyzzz[k];

                        g_zz_xyzzzz[k] = g1_z_xyzzzzz[k] - fr * g2_z_xyzzzz[k];

                        g_zz_xzzzzz[k] = g1_z_xzzzzzz[k] - fr * g2_z_xzzzzz[k];

                        g_zz_yyyyyy[k] = g1_z_yyyyyyz[k] - fr * g2_z_yyyyyy[k];

                        g_zz_yyyyyz[k] = g1_z_yyyyyzz[k] - fr * g2_z_yyyyyz[k];

                        g_zz_yyyyzz[k] = g1_z_yyyyzzz[k] - fr * g2_z_yyyyzz[k];

                        g_zz_yyyzzz[k] = g1_z_yyyzzzz[k] - fr * g2_z_yyyzzz[k];

                        g_zz_yyzzzz[k] = g1_z_yyzzzzz[k] - fr * g2_z_yyzzzz[k];

                        g_zz_yzzzzz[k] = g1_z_yzzzzzz[k] - fr * g2_z_yzzzzz[k];

                        g_zz_zzzzzz[k] = g1_z_zzzzzzz[k] - fr * g2_z_zzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXFF(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 3) && (recPattern[i].third() == 3))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|33)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 3, 3});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 4});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 3});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|DF)^(m) integrals

                    auto g2_xx_xxx = contrBuffer.data(g2off + 60 * j);

                    auto g2_xx_xxy = contrBuffer.data(g2off + 60 * j + 1);

                    auto g2_xx_xxz = contrBuffer.data(g2off + 60 * j + 2);

                    auto g2_xx_xyy = contrBuffer.data(g2off + 60 * j + 3);

                    auto g2_xx_xyz = contrBuffer.data(g2off + 60 * j + 4);

                    auto g2_xx_xzz = contrBuffer.data(g2off + 60 * j + 5);

                    auto g2_xx_yyy = contrBuffer.data(g2off + 60 * j + 6);

                    auto g2_xx_yyz = contrBuffer.data(g2off + 60 * j + 7);

                    auto g2_xx_yzz = contrBuffer.data(g2off + 60 * j + 8);

                    auto g2_xx_zzz = contrBuffer.data(g2off + 60 * j + 9);

                    auto g2_xy_xxx = contrBuffer.data(g2off + 60 * j + 10);

                    auto g2_xy_xxy = contrBuffer.data(g2off + 60 * j + 11);

                    auto g2_xy_xxz = contrBuffer.data(g2off + 60 * j + 12);

                    auto g2_xy_xyy = contrBuffer.data(g2off + 60 * j + 13);

                    auto g2_xy_xyz = contrBuffer.data(g2off + 60 * j + 14);

                    auto g2_xy_xzz = contrBuffer.data(g2off + 60 * j + 15);

                    auto g2_xy_yyy = contrBuffer.data(g2off + 60 * j + 16);

                    auto g2_xy_yyz = contrBuffer.data(g2off + 60 * j + 17);

                    auto g2_xy_yzz = contrBuffer.data(g2off + 60 * j + 18);

                    auto g2_xy_zzz = contrBuffer.data(g2off + 60 * j + 19);

                    auto g2_xz_xxx = contrBuffer.data(g2off + 60 * j + 20);

                    auto g2_xz_xxy = contrBuffer.data(g2off + 60 * j + 21);

                    auto g2_xz_xxz = contrBuffer.data(g2off + 60 * j + 22);

                    auto g2_xz_xyy = contrBuffer.data(g2off + 60 * j + 23);

                    auto g2_xz_xyz = contrBuffer.data(g2off + 60 * j + 24);

                    auto g2_xz_xzz = contrBuffer.data(g2off + 60 * j + 25);

                    auto g2_xz_yyy = contrBuffer.data(g2off + 60 * j + 26);

                    auto g2_xz_yyz = contrBuffer.data(g2off + 60 * j + 27);

                    auto g2_xz_yzz = contrBuffer.data(g2off + 60 * j + 28);

                    auto g2_xz_zzz = contrBuffer.data(g2off + 60 * j + 29);

                    auto g2_yy_xxx = contrBuffer.data(g2off + 60 * j + 30);

                    auto g2_yy_xxy = contrBuffer.data(g2off + 60 * j + 31);

                    auto g2_yy_xxz = contrBuffer.data(g2off + 60 * j + 32);

                    auto g2_yy_xyy = contrBuffer.data(g2off + 60 * j + 33);

                    auto g2_yy_xyz = contrBuffer.data(g2off + 60 * j + 34);

                    auto g2_yy_xzz = contrBuffer.data(g2off + 60 * j + 35);

                    auto g2_yy_yyy = contrBuffer.data(g2off + 60 * j + 36);

                    auto g2_yy_yyz = contrBuffer.data(g2off + 60 * j + 37);

                    auto g2_yy_yzz = contrBuffer.data(g2off + 60 * j + 38);

                    auto g2_yy_zzz = contrBuffer.data(g2off + 60 * j + 39);

                    auto g2_yz_xxx = contrBuffer.data(g2off + 60 * j + 40);

                    auto g2_yz_xxy = contrBuffer.data(g2off + 60 * j + 41);

                    auto g2_yz_xxz = contrBuffer.data(g2off + 60 * j + 42);

                    auto g2_yz_xyy = contrBuffer.data(g2off + 60 * j + 43);

                    auto g2_yz_xyz = contrBuffer.data(g2off + 60 * j + 44);

                    auto g2_yz_xzz = contrBuffer.data(g2off + 60 * j + 45);

                    auto g2_yz_yyy = contrBuffer.data(g2off + 60 * j + 46);

                    auto g2_yz_yyz = contrBuffer.data(g2off + 60 * j + 47);

                    auto g2_yz_yzz = contrBuffer.data(g2off + 60 * j + 48);

                    auto g2_yz_zzz = contrBuffer.data(g2off + 60 * j + 49);

                    auto g2_zz_xxx = contrBuffer.data(g2off + 60 * j + 50);

                    auto g2_zz_xxy = contrBuffer.data(g2off + 60 * j + 51);

                    auto g2_zz_xxz = contrBuffer.data(g2off + 60 * j + 52);

                    auto g2_zz_xyy = contrBuffer.data(g2off + 60 * j + 53);

                    auto g2_zz_xyz = contrBuffer.data(g2off + 60 * j + 54);

                    auto g2_zz_xzz = contrBuffer.data(g2off + 60 * j + 55);

                    auto g2_zz_yyy = contrBuffer.data(g2off + 60 * j + 56);

                    auto g2_zz_yyz = contrBuffer.data(g2off + 60 * j + 57);

                    auto g2_zz_yzz = contrBuffer.data(g2off + 60 * j + 58);

                    auto g2_zz_zzz = contrBuffer.data(g2off + 60 * j + 59);

                    // set up pointers to (SX|g(r,r')|DG)^(m) integrals

                    auto g1_xx_xxxx = contrBuffer.data(g1off + 90 * j);

                    auto g1_xx_xxxy = contrBuffer.data(g1off + 90 * j + 1);

                    auto g1_xx_xxxz = contrBuffer.data(g1off + 90 * j + 2);

                    auto g1_xx_xxyy = contrBuffer.data(g1off + 90 * j + 3);

                    auto g1_xx_xxyz = contrBuffer.data(g1off + 90 * j + 4);

                    auto g1_xx_xxzz = contrBuffer.data(g1off + 90 * j + 5);

                    auto g1_xx_xyyy = contrBuffer.data(g1off + 90 * j + 6);

                    auto g1_xx_xyyz = contrBuffer.data(g1off + 90 * j + 7);

                    auto g1_xx_xyzz = contrBuffer.data(g1off + 90 * j + 8);

                    auto g1_xx_xzzz = contrBuffer.data(g1off + 90 * j + 9);

                    auto g1_xy_xxxx = contrBuffer.data(g1off + 90 * j + 15);

                    auto g1_xy_xxxy = contrBuffer.data(g1off + 90 * j + 16);

                    auto g1_xy_xxxz = contrBuffer.data(g1off + 90 * j + 17);

                    auto g1_xy_xxyy = contrBuffer.data(g1off + 90 * j + 18);

                    auto g1_xy_xxyz = contrBuffer.data(g1off + 90 * j + 19);

                    auto g1_xy_xxzz = contrBuffer.data(g1off + 90 * j + 20);

                    auto g1_xy_xyyy = contrBuffer.data(g1off + 90 * j + 21);

                    auto g1_xy_xyyz = contrBuffer.data(g1off + 90 * j + 22);

                    auto g1_xy_xyzz = contrBuffer.data(g1off + 90 * j + 23);

                    auto g1_xy_xzzz = contrBuffer.data(g1off + 90 * j + 24);

                    auto g1_xz_xxxx = contrBuffer.data(g1off + 90 * j + 30);

                    auto g1_xz_xxxy = contrBuffer.data(g1off + 90 * j + 31);

                    auto g1_xz_xxxz = contrBuffer.data(g1off + 90 * j + 32);

                    auto g1_xz_xxyy = contrBuffer.data(g1off + 90 * j + 33);

                    auto g1_xz_xxyz = contrBuffer.data(g1off + 90 * j + 34);

                    auto g1_xz_xxzz = contrBuffer.data(g1off + 90 * j + 35);

                    auto g1_xz_xyyy = contrBuffer.data(g1off + 90 * j + 36);

                    auto g1_xz_xyyz = contrBuffer.data(g1off + 90 * j + 37);

                    auto g1_xz_xyzz = contrBuffer.data(g1off + 90 * j + 38);

                    auto g1_xz_xzzz = contrBuffer.data(g1off + 90 * j + 39);

                    auto g1_yy_xxxx = contrBuffer.data(g1off + 90 * j + 45);

                    auto g1_yy_xxxy = contrBuffer.data(g1off + 90 * j + 46);

                    auto g1_yy_xxxz = contrBuffer.data(g1off + 90 * j + 47);

                    auto g1_yy_xxyy = contrBuffer.data(g1off + 90 * j + 48);

                    auto g1_yy_xxyz = contrBuffer.data(g1off + 90 * j + 49);

                    auto g1_yy_xxzz = contrBuffer.data(g1off + 90 * j + 50);

                    auto g1_yy_xyyy = contrBuffer.data(g1off + 90 * j + 51);

                    auto g1_yy_xyyz = contrBuffer.data(g1off + 90 * j + 52);

                    auto g1_yy_xyzz = contrBuffer.data(g1off + 90 * j + 53);

                    auto g1_yy_xzzz = contrBuffer.data(g1off + 90 * j + 54);

                    auto g1_yy_yyyy = contrBuffer.data(g1off + 90 * j + 55);

                    auto g1_yy_yyyz = contrBuffer.data(g1off + 90 * j + 56);

                    auto g1_yy_yyzz = contrBuffer.data(g1off + 90 * j + 57);

                    auto g1_yy_yzzz = contrBuffer.data(g1off + 90 * j + 58);

                    auto g1_yz_xxxx = contrBuffer.data(g1off + 90 * j + 60);

                    auto g1_yz_xxxy = contrBuffer.data(g1off + 90 * j + 61);

                    auto g1_yz_xxxz = contrBuffer.data(g1off + 90 * j + 62);

                    auto g1_yz_xxyy = contrBuffer.data(g1off + 90 * j + 63);

                    auto g1_yz_xxyz = contrBuffer.data(g1off + 90 * j + 64);

                    auto g1_yz_xxzz = contrBuffer.data(g1off + 90 * j + 65);

                    auto g1_yz_xyyy = contrBuffer.data(g1off + 90 * j + 66);

                    auto g1_yz_xyyz = contrBuffer.data(g1off + 90 * j + 67);

                    auto g1_yz_xyzz = contrBuffer.data(g1off + 90 * j + 68);

                    auto g1_yz_xzzz = contrBuffer.data(g1off + 90 * j + 69);

                    auto g1_yz_yyyy = contrBuffer.data(g1off + 90 * j + 70);

                    auto g1_yz_yyyz = contrBuffer.data(g1off + 90 * j + 71);

                    auto g1_yz_yyzz = contrBuffer.data(g1off + 90 * j + 72);

                    auto g1_yz_yzzz = contrBuffer.data(g1off + 90 * j + 73);

                    auto g1_zz_xxxx = contrBuffer.data(g1off + 90 * j + 75);

                    auto g1_zz_xxxy = contrBuffer.data(g1off + 90 * j + 76);

                    auto g1_zz_xxxz = contrBuffer.data(g1off + 90 * j + 77);

                    auto g1_zz_xxyy = contrBuffer.data(g1off + 90 * j + 78);

                    auto g1_zz_xxyz = contrBuffer.data(g1off + 90 * j + 79);

                    auto g1_zz_xxzz = contrBuffer.data(g1off + 90 * j + 80);

                    auto g1_zz_xyyy = contrBuffer.data(g1off + 90 * j + 81);

                    auto g1_zz_xyyz = contrBuffer.data(g1off + 90 * j + 82);

                    auto g1_zz_xyzz = contrBuffer.data(g1off + 90 * j + 83);

                    auto g1_zz_xzzz = contrBuffer.data(g1off + 90 * j + 84);

                    auto g1_zz_yyyy = contrBuffer.data(g1off + 90 * j + 85);

                    auto g1_zz_yyyz = contrBuffer.data(g1off + 90 * j + 86);

                    auto g1_zz_yyzz = contrBuffer.data(g1off + 90 * j + 87);

                    auto g1_zz_yzzz = contrBuffer.data(g1off + 90 * j + 88);

                    auto g1_zz_zzzz = contrBuffer.data(g1off + 90 * j + 89);

                    // set up pointers to (SX|g(r,r')|FF)^(m) integrals

                    auto g_xxx_xxx = contrBuffer.data(goff + 100 * j);

                    auto g_xxx_xxy = contrBuffer.data(goff + 100 * j + 1);

                    auto g_xxx_xxz = contrBuffer.data(goff + 100 * j + 2);

                    auto g_xxx_xyy = contrBuffer.data(goff + 100 * j + 3);

                    auto g_xxx_xyz = contrBuffer.data(goff + 100 * j + 4);

                    auto g_xxx_xzz = contrBuffer.data(goff + 100 * j + 5);

                    auto g_xxx_yyy = contrBuffer.data(goff + 100 * j + 6);

                    auto g_xxx_yyz = contrBuffer.data(goff + 100 * j + 7);

                    auto g_xxx_yzz = contrBuffer.data(goff + 100 * j + 8);

                    auto g_xxx_zzz = contrBuffer.data(goff + 100 * j + 9);

                    auto g_xxy_xxx = contrBuffer.data(goff + 100 * j + 10);

                    auto g_xxy_xxy = contrBuffer.data(goff + 100 * j + 11);

                    auto g_xxy_xxz = contrBuffer.data(goff + 100 * j + 12);

                    auto g_xxy_xyy = contrBuffer.data(goff + 100 * j + 13);

                    auto g_xxy_xyz = contrBuffer.data(goff + 100 * j + 14);

                    auto g_xxy_xzz = contrBuffer.data(goff + 100 * j + 15);

                    auto g_xxy_yyy = contrBuffer.data(goff + 100 * j + 16);

                    auto g_xxy_yyz = contrBuffer.data(goff + 100 * j + 17);

                    auto g_xxy_yzz = contrBuffer.data(goff + 100 * j + 18);

                    auto g_xxy_zzz = contrBuffer.data(goff + 100 * j + 19);

                    auto g_xxz_xxx = contrBuffer.data(goff + 100 * j + 20);

                    auto g_xxz_xxy = contrBuffer.data(goff + 100 * j + 21);

                    auto g_xxz_xxz = contrBuffer.data(goff + 100 * j + 22);

                    auto g_xxz_xyy = contrBuffer.data(goff + 100 * j + 23);

                    auto g_xxz_xyz = contrBuffer.data(goff + 100 * j + 24);

                    auto g_xxz_xzz = contrBuffer.data(goff + 100 * j + 25);

                    auto g_xxz_yyy = contrBuffer.data(goff + 100 * j + 26);

                    auto g_xxz_yyz = contrBuffer.data(goff + 100 * j + 27);

                    auto g_xxz_yzz = contrBuffer.data(goff + 100 * j + 28);

                    auto g_xxz_zzz = contrBuffer.data(goff + 100 * j + 29);

                    auto g_xyy_xxx = contrBuffer.data(goff + 100 * j + 30);

                    auto g_xyy_xxy = contrBuffer.data(goff + 100 * j + 31);

                    auto g_xyy_xxz = contrBuffer.data(goff + 100 * j + 32);

                    auto g_xyy_xyy = contrBuffer.data(goff + 100 * j + 33);

                    auto g_xyy_xyz = contrBuffer.data(goff + 100 * j + 34);

                    auto g_xyy_xzz = contrBuffer.data(goff + 100 * j + 35);

                    auto g_xyy_yyy = contrBuffer.data(goff + 100 * j + 36);

                    auto g_xyy_yyz = contrBuffer.data(goff + 100 * j + 37);

                    auto g_xyy_yzz = contrBuffer.data(goff + 100 * j + 38);

                    auto g_xyy_zzz = contrBuffer.data(goff + 100 * j + 39);

                    auto g_xyz_xxx = contrBuffer.data(goff + 100 * j + 40);

                    auto g_xyz_xxy = contrBuffer.data(goff + 100 * j + 41);

                    auto g_xyz_xxz = contrBuffer.data(goff + 100 * j + 42);

                    auto g_xyz_xyy = contrBuffer.data(goff + 100 * j + 43);

                    auto g_xyz_xyz = contrBuffer.data(goff + 100 * j + 44);

                    auto g_xyz_xzz = contrBuffer.data(goff + 100 * j + 45);

                    auto g_xyz_yyy = contrBuffer.data(goff + 100 * j + 46);

                    auto g_xyz_yyz = contrBuffer.data(goff + 100 * j + 47);

                    auto g_xyz_yzz = contrBuffer.data(goff + 100 * j + 48);

                    auto g_xyz_zzz = contrBuffer.data(goff + 100 * j + 49);

                    auto g_xzz_xxx = contrBuffer.data(goff + 100 * j + 50);

                    auto g_xzz_xxy = contrBuffer.data(goff + 100 * j + 51);

                    auto g_xzz_xxz = contrBuffer.data(goff + 100 * j + 52);

                    auto g_xzz_xyy = contrBuffer.data(goff + 100 * j + 53);

                    auto g_xzz_xyz = contrBuffer.data(goff + 100 * j + 54);

                    auto g_xzz_xzz = contrBuffer.data(goff + 100 * j + 55);

                    auto g_xzz_yyy = contrBuffer.data(goff + 100 * j + 56);

                    auto g_xzz_yyz = contrBuffer.data(goff + 100 * j + 57);

                    auto g_xzz_yzz = contrBuffer.data(goff + 100 * j + 58);

                    auto g_xzz_zzz = contrBuffer.data(goff + 100 * j + 59);

                    auto g_yyy_xxx = contrBuffer.data(goff + 100 * j + 60);

                    auto g_yyy_xxy = contrBuffer.data(goff + 100 * j + 61);

                    auto g_yyy_xxz = contrBuffer.data(goff + 100 * j + 62);

                    auto g_yyy_xyy = contrBuffer.data(goff + 100 * j + 63);

                    auto g_yyy_xyz = contrBuffer.data(goff + 100 * j + 64);

                    auto g_yyy_xzz = contrBuffer.data(goff + 100 * j + 65);

                    auto g_yyy_yyy = contrBuffer.data(goff + 100 * j + 66);

                    auto g_yyy_yyz = contrBuffer.data(goff + 100 * j + 67);

                    auto g_yyy_yzz = contrBuffer.data(goff + 100 * j + 68);

                    auto g_yyy_zzz = contrBuffer.data(goff + 100 * j + 69);

                    auto g_yyz_xxx = contrBuffer.data(goff + 100 * j + 70);

                    auto g_yyz_xxy = contrBuffer.data(goff + 100 * j + 71);

                    auto g_yyz_xxz = contrBuffer.data(goff + 100 * j + 72);

                    auto g_yyz_xyy = contrBuffer.data(goff + 100 * j + 73);

                    auto g_yyz_xyz = contrBuffer.data(goff + 100 * j + 74);

                    auto g_yyz_xzz = contrBuffer.data(goff + 100 * j + 75);

                    auto g_yyz_yyy = contrBuffer.data(goff + 100 * j + 76);

                    auto g_yyz_yyz = contrBuffer.data(goff + 100 * j + 77);

                    auto g_yyz_yzz = contrBuffer.data(goff + 100 * j + 78);

                    auto g_yyz_zzz = contrBuffer.data(goff + 100 * j + 79);

                    auto g_yzz_xxx = contrBuffer.data(goff + 100 * j + 80);

                    auto g_yzz_xxy = contrBuffer.data(goff + 100 * j + 81);

                    auto g_yzz_xxz = contrBuffer.data(goff + 100 * j + 82);

                    auto g_yzz_xyy = contrBuffer.data(goff + 100 * j + 83);

                    auto g_yzz_xyz = contrBuffer.data(goff + 100 * j + 84);

                    auto g_yzz_xzz = contrBuffer.data(goff + 100 * j + 85);

                    auto g_yzz_yyy = contrBuffer.data(goff + 100 * j + 86);

                    auto g_yzz_yyz = contrBuffer.data(goff + 100 * j + 87);

                    auto g_yzz_yzz = contrBuffer.data(goff + 100 * j + 88);

                    auto g_yzz_zzz = contrBuffer.data(goff + 100 * j + 89);

                    auto g_zzz_xxx = contrBuffer.data(goff + 100 * j + 90);

                    auto g_zzz_xxy = contrBuffer.data(goff + 100 * j + 91);

                    auto g_zzz_xxz = contrBuffer.data(goff + 100 * j + 92);

                    auto g_zzz_xyy = contrBuffer.data(goff + 100 * j + 93);

                    auto g_zzz_xyz = contrBuffer.data(goff + 100 * j + 94);

                    auto g_zzz_xzz = contrBuffer.data(goff + 100 * j + 95);

                    auto g_zzz_yyy = contrBuffer.data(goff + 100 * j + 96);

                    auto g_zzz_yyz = contrBuffer.data(goff + 100 * j + 97);

                    auto g_zzz_yzz = contrBuffer.data(goff + 100 * j + 98);

                    auto g_zzz_zzz = contrBuffer.data(goff + 100 * j + 99);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxx, g2_xx_xxy,\
                                             g2_xx_xxz, g2_xx_xyy, g2_xx_xyz, g2_xx_xzz,\
                                             g2_xx_yyy, g2_xx_yyz, g2_xx_yzz, g2_xx_zzz,\
                                             g2_xy_xxx, g2_xy_xxy, g2_xy_xxz, g2_xy_xyy,\
                                             g2_xy_xyz, g2_xy_xzz, g2_xy_yyy, g2_xy_yyz,\
                                             g2_xy_yzz, g2_xy_zzz, g2_xz_xxx, g2_xz_xxy,\
                                             g2_xz_xxz, g2_xz_xyy, g2_xz_xyz, g2_xz_xzz,\
                                             g2_xz_yyy, g2_xz_yyz, g2_xz_yzz, g2_xz_zzz,\
                                             g2_yy_xxx, g2_yy_xxy, g2_yy_xxz, g2_yy_xyy,\
                                             g2_yy_xyz, g2_yy_xzz, g2_yy_yyy, g2_yy_yyz,\
                                             g2_yy_yzz, g2_yy_zzz, g2_yz_xxx, g2_yz_xxy,\
                                             g2_yz_xxz, g2_yz_xyy, g2_yz_xyz, g2_yz_xzz,\
                                             g2_yz_yyy, g2_yz_yyz, g2_yz_yzz, g2_yz_zzz,\
                                             g2_zz_xxx, g2_zz_xxy, g2_zz_xxz, g2_zz_xyy,\
                                             g2_zz_xyz, g2_zz_xzz, g2_zz_yyy, g2_zz_yyz,\
                                             g2_zz_yzz, g2_zz_zzz, g1_xx_xxxx,\
                                             g1_xx_xxxy, g1_xx_xxxz, g1_xx_xxyy,\
                                             g1_xx_xxyz, g1_xx_xxzz, g1_xx_xyyy,\
                                             g1_xx_xyyz, g1_xx_xyzz, g1_xx_xzzz,\
                                             g1_xy_xxxx,\
                                             g1_xy_xxxy, g1_xy_xxxz, g1_xy_xxyy,\
                                             g1_xy_xxyz, g1_xy_xxzz, g1_xy_xyyy,\
                                             g1_xy_xyyz, g1_xy_xyzz, g1_xy_xzzz,\
                                             g1_xz_xxxx,\
                                             g1_xz_xxxy, g1_xz_xxxz, g1_xz_xxyy,\
                                             g1_xz_xxyz, g1_xz_xxzz, g1_xz_xyyy,\
                                             g1_xz_xyyz, g1_xz_xyzz, g1_xz_xzzz,\
                                             g1_yy_xxxx,\
                                             g1_yy_xxxy, g1_yy_xxxz, g1_yy_xxyy,\
                                             g1_yy_xxyz, g1_yy_xxzz, g1_yy_xyyy,\
                                             g1_yy_xyyz, g1_yy_xyzz, g1_yy_xzzz,\
                                             g1_yy_yyyy, g1_yy_yyyz, g1_yy_yyzz,\
                                             g1_yy_yzzz, g1_yz_xxxx,\
                                             g1_yz_xxxy, g1_yz_xxxz, g1_yz_xxyy,\
                                             g1_yz_xxyz, g1_yz_xxzz, g1_yz_xyyy,\
                                             g1_yz_xyyz, g1_yz_xyzz, g1_yz_xzzz,\
                                             g1_yz_yyyy, g1_yz_yyyz, g1_yz_yyzz,\
                                             g1_yz_yzzz, g1_zz_xxxx,\
                                             g1_zz_xxxy, g1_zz_xxxz, g1_zz_xxyy,\
                                             g1_zz_xxyz, g1_zz_xxzz, g1_zz_xyyy,\
                                             g1_zz_xyyz, g1_zz_xyzz, g1_zz_xzzz,\
                                             g1_zz_yyyy, g1_zz_yyyz, g1_zz_yyzz,\
                                             g1_zz_yzzz, g1_zz_zzzz, g_xxx_xxx,\
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
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xxx_xxx[k] = g1_xx_xxxx[k] - fr * g2_xx_xxx[k];

                        g_xxx_xxy[k] = g1_xx_xxxy[k] - fr * g2_xx_xxy[k];

                        g_xxx_xxz[k] = g1_xx_xxxz[k] - fr * g2_xx_xxz[k];

                        g_xxx_xyy[k] = g1_xx_xxyy[k] - fr * g2_xx_xyy[k];

                        g_xxx_xyz[k] = g1_xx_xxyz[k] - fr * g2_xx_xyz[k];

                        g_xxx_xzz[k] = g1_xx_xxzz[k] - fr * g2_xx_xzz[k];

                        g_xxx_yyy[k] = g1_xx_xyyy[k] - fr * g2_xx_yyy[k];

                        g_xxx_yyz[k] = g1_xx_xyyz[k] - fr * g2_xx_yyz[k];

                        g_xxx_yzz[k] = g1_xx_xyzz[k] - fr * g2_xx_yzz[k];

                        g_xxx_zzz[k] = g1_xx_xzzz[k] - fr * g2_xx_zzz[k];

                        g_xxy_xxx[k] = g1_xy_xxxx[k] - fr * g2_xy_xxx[k];

                        g_xxy_xxy[k] = g1_xy_xxxy[k] - fr * g2_xy_xxy[k];

                        g_xxy_xxz[k] = g1_xy_xxxz[k] - fr * g2_xy_xxz[k];

                        g_xxy_xyy[k] = g1_xy_xxyy[k] - fr * g2_xy_xyy[k];

                        g_xxy_xyz[k] = g1_xy_xxyz[k] - fr * g2_xy_xyz[k];

                        g_xxy_xzz[k] = g1_xy_xxzz[k] - fr * g2_xy_xzz[k];

                        g_xxy_yyy[k] = g1_xy_xyyy[k] - fr * g2_xy_yyy[k];

                        g_xxy_yyz[k] = g1_xy_xyyz[k] - fr * g2_xy_yyz[k];

                        g_xxy_yzz[k] = g1_xy_xyzz[k] - fr * g2_xy_yzz[k];

                        g_xxy_zzz[k] = g1_xy_xzzz[k] - fr * g2_xy_zzz[k];

                        g_xxz_xxx[k] = g1_xz_xxxx[k] - fr * g2_xz_xxx[k];

                        g_xxz_xxy[k] = g1_xz_xxxy[k] - fr * g2_xz_xxy[k];

                        g_xxz_xxz[k] = g1_xz_xxxz[k] - fr * g2_xz_xxz[k];

                        g_xxz_xyy[k] = g1_xz_xxyy[k] - fr * g2_xz_xyy[k];

                        g_xxz_xyz[k] = g1_xz_xxyz[k] - fr * g2_xz_xyz[k];

                        g_xxz_xzz[k] = g1_xz_xxzz[k] - fr * g2_xz_xzz[k];

                        g_xxz_yyy[k] = g1_xz_xyyy[k] - fr * g2_xz_yyy[k];

                        g_xxz_yyz[k] = g1_xz_xyyz[k] - fr * g2_xz_yyz[k];

                        g_xxz_yzz[k] = g1_xz_xyzz[k] - fr * g2_xz_yzz[k];

                        g_xxz_zzz[k] = g1_xz_xzzz[k] - fr * g2_xz_zzz[k];

                        g_xyy_xxx[k] = g1_yy_xxxx[k] - fr * g2_yy_xxx[k];

                        g_xyy_xxy[k] = g1_yy_xxxy[k] - fr * g2_yy_xxy[k];

                        g_xyy_xxz[k] = g1_yy_xxxz[k] - fr * g2_yy_xxz[k];

                        g_xyy_xyy[k] = g1_yy_xxyy[k] - fr * g2_yy_xyy[k];

                        g_xyy_xyz[k] = g1_yy_xxyz[k] - fr * g2_yy_xyz[k];

                        g_xyy_xzz[k] = g1_yy_xxzz[k] - fr * g2_yy_xzz[k];

                        g_xyy_yyy[k] = g1_yy_xyyy[k] - fr * g2_yy_yyy[k];

                        g_xyy_yyz[k] = g1_yy_xyyz[k] - fr * g2_yy_yyz[k];

                        g_xyy_yzz[k] = g1_yy_xyzz[k] - fr * g2_yy_yzz[k];

                        g_xyy_zzz[k] = g1_yy_xzzz[k] - fr * g2_yy_zzz[k];

                        g_xyz_xxx[k] = g1_yz_xxxx[k] - fr * g2_yz_xxx[k];

                        g_xyz_xxy[k] = g1_yz_xxxy[k] - fr * g2_yz_xxy[k];

                        g_xyz_xxz[k] = g1_yz_xxxz[k] - fr * g2_yz_xxz[k];

                        g_xyz_xyy[k] = g1_yz_xxyy[k] - fr * g2_yz_xyy[k];

                        g_xyz_xyz[k] = g1_yz_xxyz[k] - fr * g2_yz_xyz[k];

                        g_xyz_xzz[k] = g1_yz_xxzz[k] - fr * g2_yz_xzz[k];

                        g_xyz_yyy[k] = g1_yz_xyyy[k] - fr * g2_yz_yyy[k];

                        g_xyz_yyz[k] = g1_yz_xyyz[k] - fr * g2_yz_yyz[k];

                        g_xyz_yzz[k] = g1_yz_xyzz[k] - fr * g2_yz_yzz[k];

                        g_xyz_zzz[k] = g1_yz_xzzz[k] - fr * g2_yz_zzz[k];

                        g_xzz_xxx[k] = g1_zz_xxxx[k] - fr * g2_zz_xxx[k];

                        g_xzz_xxy[k] = g1_zz_xxxy[k] - fr * g2_zz_xxy[k];

                        g_xzz_xxz[k] = g1_zz_xxxz[k] - fr * g2_zz_xxz[k];

                        g_xzz_xyy[k] = g1_zz_xxyy[k] - fr * g2_zz_xyy[k];

                        g_xzz_xyz[k] = g1_zz_xxyz[k] - fr * g2_zz_xyz[k];

                        g_xzz_xzz[k] = g1_zz_xxzz[k] - fr * g2_zz_xzz[k];

                        g_xzz_yyy[k] = g1_zz_xyyy[k] - fr * g2_zz_yyy[k];

                        g_xzz_yyz[k] = g1_zz_xyyz[k] - fr * g2_zz_yyz[k];

                        g_xzz_yzz[k] = g1_zz_xyzz[k] - fr * g2_zz_yzz[k];

                        g_xzz_zzz[k] = g1_zz_xzzz[k] - fr * g2_zz_zzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yyy_xxx[k] = g1_yy_xxxy[k] - fr * g2_yy_xxx[k];

                        g_yyy_xxy[k] = g1_yy_xxyy[k] - fr * g2_yy_xxy[k];

                        g_yyy_xxz[k] = g1_yy_xxyz[k] - fr * g2_yy_xxz[k];

                        g_yyy_xyy[k] = g1_yy_xyyy[k] - fr * g2_yy_xyy[k];

                        g_yyy_xyz[k] = g1_yy_xyyz[k] - fr * g2_yy_xyz[k];

                        g_yyy_xzz[k] = g1_yy_xyzz[k] - fr * g2_yy_xzz[k];

                        g_yyy_yyy[k] = g1_yy_yyyy[k] - fr * g2_yy_yyy[k];

                        g_yyy_yyz[k] = g1_yy_yyyz[k] - fr * g2_yy_yyz[k];

                        g_yyy_yzz[k] = g1_yy_yyzz[k] - fr * g2_yy_yzz[k];

                        g_yyy_zzz[k] = g1_yy_yzzz[k] - fr * g2_yy_zzz[k];

                        g_yyz_xxx[k] = g1_yz_xxxy[k] - fr * g2_yz_xxx[k];

                        g_yyz_xxy[k] = g1_yz_xxyy[k] - fr * g2_yz_xxy[k];

                        g_yyz_xxz[k] = g1_yz_xxyz[k] - fr * g2_yz_xxz[k];

                        g_yyz_xyy[k] = g1_yz_xyyy[k] - fr * g2_yz_xyy[k];

                        g_yyz_xyz[k] = g1_yz_xyyz[k] - fr * g2_yz_xyz[k];

                        g_yyz_xzz[k] = g1_yz_xyzz[k] - fr * g2_yz_xzz[k];

                        g_yyz_yyy[k] = g1_yz_yyyy[k] - fr * g2_yz_yyy[k];

                        g_yyz_yyz[k] = g1_yz_yyyz[k] - fr * g2_yz_yyz[k];

                        g_yyz_yzz[k] = g1_yz_yyzz[k] - fr * g2_yz_yzz[k];

                        g_yyz_zzz[k] = g1_yz_yzzz[k] - fr * g2_yz_zzz[k];

                        g_yzz_xxx[k] = g1_zz_xxxy[k] - fr * g2_zz_xxx[k];

                        g_yzz_xxy[k] = g1_zz_xxyy[k] - fr * g2_zz_xxy[k];

                        g_yzz_xxz[k] = g1_zz_xxyz[k] - fr * g2_zz_xxz[k];

                        g_yzz_xyy[k] = g1_zz_xyyy[k] - fr * g2_zz_xyy[k];

                        g_yzz_xyz[k] = g1_zz_xyyz[k] - fr * g2_zz_xyz[k];

                        g_yzz_xzz[k] = g1_zz_xyzz[k] - fr * g2_zz_xzz[k];

                        g_yzz_yyy[k] = g1_zz_yyyy[k] - fr * g2_zz_yyy[k];

                        g_yzz_yyz[k] = g1_zz_yyyz[k] - fr * g2_zz_yyz[k];

                        g_yzz_yzz[k] = g1_zz_yyzz[k] - fr * g2_zz_yzz[k];

                        g_yzz_zzz[k] = g1_zz_yzzz[k] - fr * g2_zz_zzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zzz_xxx[k] = g1_zz_xxxz[k] - fr * g2_zz_xxx[k];

                        g_zzz_xxy[k] = g1_zz_xxyz[k] - fr * g2_zz_xxy[k];

                        g_zzz_xxz[k] = g1_zz_xxzz[k] - fr * g2_zz_xxz[k];

                        g_zzz_xyy[k] = g1_zz_xyyz[k] - fr * g2_zz_xyy[k];

                        g_zzz_xyz[k] = g1_zz_xyzz[k] - fr * g2_zz_xyz[k];

                        g_zzz_xzz[k] = g1_zz_xzzz[k] - fr * g2_zz_xzz[k];

                        g_zzz_yyy[k] = g1_zz_yyyz[k] - fr * g2_zz_yyy[k];

                        g_zzz_yyz[k] = g1_zz_yyzz[k] - fr * g2_zz_yyz[k];

                        g_zzz_yzz[k] = g1_zz_yzzz[k] - fr * g2_zz_yzz[k];

                        g_zzz_zzz[k] = g1_zz_zzzz[k] - fr * g2_zz_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXFG(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 3) && (recPattern[i].third() == 4))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|34)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 3, 4});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 5});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 4});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|DG)^(m) integrals

                    auto g2_xx_xxxx = contrBuffer.data(g2off + 90 * j);

                    auto g2_xx_xxxy = contrBuffer.data(g2off + 90 * j + 1);

                    auto g2_xx_xxxz = contrBuffer.data(g2off + 90 * j + 2);

                    auto g2_xx_xxyy = contrBuffer.data(g2off + 90 * j + 3);

                    auto g2_xx_xxyz = contrBuffer.data(g2off + 90 * j + 4);

                    auto g2_xx_xxzz = contrBuffer.data(g2off + 90 * j + 5);

                    auto g2_xx_xyyy = contrBuffer.data(g2off + 90 * j + 6);

                    auto g2_xx_xyyz = contrBuffer.data(g2off + 90 * j + 7);

                    auto g2_xx_xyzz = contrBuffer.data(g2off + 90 * j + 8);

                    auto g2_xx_xzzz = contrBuffer.data(g2off + 90 * j + 9);

                    auto g2_xx_yyyy = contrBuffer.data(g2off + 90 * j + 10);

                    auto g2_xx_yyyz = contrBuffer.data(g2off + 90 * j + 11);

                    auto g2_xx_yyzz = contrBuffer.data(g2off + 90 * j + 12);

                    auto g2_xx_yzzz = contrBuffer.data(g2off + 90 * j + 13);

                    auto g2_xx_zzzz = contrBuffer.data(g2off + 90 * j + 14);

                    auto g2_xy_xxxx = contrBuffer.data(g2off + 90 * j + 15);

                    auto g2_xy_xxxy = contrBuffer.data(g2off + 90 * j + 16);

                    auto g2_xy_xxxz = contrBuffer.data(g2off + 90 * j + 17);

                    auto g2_xy_xxyy = contrBuffer.data(g2off + 90 * j + 18);

                    auto g2_xy_xxyz = contrBuffer.data(g2off + 90 * j + 19);

                    auto g2_xy_xxzz = contrBuffer.data(g2off + 90 * j + 20);

                    auto g2_xy_xyyy = contrBuffer.data(g2off + 90 * j + 21);

                    auto g2_xy_xyyz = contrBuffer.data(g2off + 90 * j + 22);

                    auto g2_xy_xyzz = contrBuffer.data(g2off + 90 * j + 23);

                    auto g2_xy_xzzz = contrBuffer.data(g2off + 90 * j + 24);

                    auto g2_xy_yyyy = contrBuffer.data(g2off + 90 * j + 25);

                    auto g2_xy_yyyz = contrBuffer.data(g2off + 90 * j + 26);

                    auto g2_xy_yyzz = contrBuffer.data(g2off + 90 * j + 27);

                    auto g2_xy_yzzz = contrBuffer.data(g2off + 90 * j + 28);

                    auto g2_xy_zzzz = contrBuffer.data(g2off + 90 * j + 29);

                    auto g2_xz_xxxx = contrBuffer.data(g2off + 90 * j + 30);

                    auto g2_xz_xxxy = contrBuffer.data(g2off + 90 * j + 31);

                    auto g2_xz_xxxz = contrBuffer.data(g2off + 90 * j + 32);

                    auto g2_xz_xxyy = contrBuffer.data(g2off + 90 * j + 33);

                    auto g2_xz_xxyz = contrBuffer.data(g2off + 90 * j + 34);

                    auto g2_xz_xxzz = contrBuffer.data(g2off + 90 * j + 35);

                    auto g2_xz_xyyy = contrBuffer.data(g2off + 90 * j + 36);

                    auto g2_xz_xyyz = contrBuffer.data(g2off + 90 * j + 37);

                    auto g2_xz_xyzz = contrBuffer.data(g2off + 90 * j + 38);

                    auto g2_xz_xzzz = contrBuffer.data(g2off + 90 * j + 39);

                    auto g2_xz_yyyy = contrBuffer.data(g2off + 90 * j + 40);

                    auto g2_xz_yyyz = contrBuffer.data(g2off + 90 * j + 41);

                    auto g2_xz_yyzz = contrBuffer.data(g2off + 90 * j + 42);

                    auto g2_xz_yzzz = contrBuffer.data(g2off + 90 * j + 43);

                    auto g2_xz_zzzz = contrBuffer.data(g2off + 90 * j + 44);

                    auto g2_yy_xxxx = contrBuffer.data(g2off + 90 * j + 45);

                    auto g2_yy_xxxy = contrBuffer.data(g2off + 90 * j + 46);

                    auto g2_yy_xxxz = contrBuffer.data(g2off + 90 * j + 47);

                    auto g2_yy_xxyy = contrBuffer.data(g2off + 90 * j + 48);

                    auto g2_yy_xxyz = contrBuffer.data(g2off + 90 * j + 49);

                    auto g2_yy_xxzz = contrBuffer.data(g2off + 90 * j + 50);

                    auto g2_yy_xyyy = contrBuffer.data(g2off + 90 * j + 51);

                    auto g2_yy_xyyz = contrBuffer.data(g2off + 90 * j + 52);

                    auto g2_yy_xyzz = contrBuffer.data(g2off + 90 * j + 53);

                    auto g2_yy_xzzz = contrBuffer.data(g2off + 90 * j + 54);

                    auto g2_yy_yyyy = contrBuffer.data(g2off + 90 * j + 55);

                    auto g2_yy_yyyz = contrBuffer.data(g2off + 90 * j + 56);

                    auto g2_yy_yyzz = contrBuffer.data(g2off + 90 * j + 57);

                    auto g2_yy_yzzz = contrBuffer.data(g2off + 90 * j + 58);

                    auto g2_yy_zzzz = contrBuffer.data(g2off + 90 * j + 59);

                    auto g2_yz_xxxx = contrBuffer.data(g2off + 90 * j + 60);

                    auto g2_yz_xxxy = contrBuffer.data(g2off + 90 * j + 61);

                    auto g2_yz_xxxz = contrBuffer.data(g2off + 90 * j + 62);

                    auto g2_yz_xxyy = contrBuffer.data(g2off + 90 * j + 63);

                    auto g2_yz_xxyz = contrBuffer.data(g2off + 90 * j + 64);

                    auto g2_yz_xxzz = contrBuffer.data(g2off + 90 * j + 65);

                    auto g2_yz_xyyy = contrBuffer.data(g2off + 90 * j + 66);

                    auto g2_yz_xyyz = contrBuffer.data(g2off + 90 * j + 67);

                    auto g2_yz_xyzz = contrBuffer.data(g2off + 90 * j + 68);

                    auto g2_yz_xzzz = contrBuffer.data(g2off + 90 * j + 69);

                    auto g2_yz_yyyy = contrBuffer.data(g2off + 90 * j + 70);

                    auto g2_yz_yyyz = contrBuffer.data(g2off + 90 * j + 71);

                    auto g2_yz_yyzz = contrBuffer.data(g2off + 90 * j + 72);

                    auto g2_yz_yzzz = contrBuffer.data(g2off + 90 * j + 73);

                    auto g2_yz_zzzz = contrBuffer.data(g2off + 90 * j + 74);

                    auto g2_zz_xxxx = contrBuffer.data(g2off + 90 * j + 75);

                    auto g2_zz_xxxy = contrBuffer.data(g2off + 90 * j + 76);

                    auto g2_zz_xxxz = contrBuffer.data(g2off + 90 * j + 77);

                    auto g2_zz_xxyy = contrBuffer.data(g2off + 90 * j + 78);

                    auto g2_zz_xxyz = contrBuffer.data(g2off + 90 * j + 79);

                    auto g2_zz_xxzz = contrBuffer.data(g2off + 90 * j + 80);

                    auto g2_zz_xyyy = contrBuffer.data(g2off + 90 * j + 81);

                    auto g2_zz_xyyz = contrBuffer.data(g2off + 90 * j + 82);

                    auto g2_zz_xyzz = contrBuffer.data(g2off + 90 * j + 83);

                    auto g2_zz_xzzz = contrBuffer.data(g2off + 90 * j + 84);

                    auto g2_zz_yyyy = contrBuffer.data(g2off + 90 * j + 85);

                    auto g2_zz_yyyz = contrBuffer.data(g2off + 90 * j + 86);

                    auto g2_zz_yyzz = contrBuffer.data(g2off + 90 * j + 87);

                    auto g2_zz_yzzz = contrBuffer.data(g2off + 90 * j + 88);

                    auto g2_zz_zzzz = contrBuffer.data(g2off + 90 * j + 89);

                    // set up pointers to (SX|g(r,r')|DH)^(m) integrals

                    auto g1_xx_xxxxx = contrBuffer.data(g1off + 126 * j);

                    auto g1_xx_xxxxy = contrBuffer.data(g1off + 126 * j + 1);

                    auto g1_xx_xxxxz = contrBuffer.data(g1off + 126 * j + 2);

                    auto g1_xx_xxxyy = contrBuffer.data(g1off + 126 * j + 3);

                    auto g1_xx_xxxyz = contrBuffer.data(g1off + 126 * j + 4);

                    auto g1_xx_xxxzz = contrBuffer.data(g1off + 126 * j + 5);

                    auto g1_xx_xxyyy = contrBuffer.data(g1off + 126 * j + 6);

                    auto g1_xx_xxyyz = contrBuffer.data(g1off + 126 * j + 7);

                    auto g1_xx_xxyzz = contrBuffer.data(g1off + 126 * j + 8);

                    auto g1_xx_xxzzz = contrBuffer.data(g1off + 126 * j + 9);

                    auto g1_xx_xyyyy = contrBuffer.data(g1off + 126 * j + 10);

                    auto g1_xx_xyyyz = contrBuffer.data(g1off + 126 * j + 11);

                    auto g1_xx_xyyzz = contrBuffer.data(g1off + 126 * j + 12);

                    auto g1_xx_xyzzz = contrBuffer.data(g1off + 126 * j + 13);

                    auto g1_xx_xzzzz = contrBuffer.data(g1off + 126 * j + 14);

                    auto g1_xy_xxxxx = contrBuffer.data(g1off + 126 * j + 21);

                    auto g1_xy_xxxxy = contrBuffer.data(g1off + 126 * j + 22);

                    auto g1_xy_xxxxz = contrBuffer.data(g1off + 126 * j + 23);

                    auto g1_xy_xxxyy = contrBuffer.data(g1off + 126 * j + 24);

                    auto g1_xy_xxxyz = contrBuffer.data(g1off + 126 * j + 25);

                    auto g1_xy_xxxzz = contrBuffer.data(g1off + 126 * j + 26);

                    auto g1_xy_xxyyy = contrBuffer.data(g1off + 126 * j + 27);

                    auto g1_xy_xxyyz = contrBuffer.data(g1off + 126 * j + 28);

                    auto g1_xy_xxyzz = contrBuffer.data(g1off + 126 * j + 29);

                    auto g1_xy_xxzzz = contrBuffer.data(g1off + 126 * j + 30);

                    auto g1_xy_xyyyy = contrBuffer.data(g1off + 126 * j + 31);

                    auto g1_xy_xyyyz = contrBuffer.data(g1off + 126 * j + 32);

                    auto g1_xy_xyyzz = contrBuffer.data(g1off + 126 * j + 33);

                    auto g1_xy_xyzzz = contrBuffer.data(g1off + 126 * j + 34);

                    auto g1_xy_xzzzz = contrBuffer.data(g1off + 126 * j + 35);

                    auto g1_xz_xxxxx = contrBuffer.data(g1off + 126 * j + 42);

                    auto g1_xz_xxxxy = contrBuffer.data(g1off + 126 * j + 43);

                    auto g1_xz_xxxxz = contrBuffer.data(g1off + 126 * j + 44);

                    auto g1_xz_xxxyy = contrBuffer.data(g1off + 126 * j + 45);

                    auto g1_xz_xxxyz = contrBuffer.data(g1off + 126 * j + 46);

                    auto g1_xz_xxxzz = contrBuffer.data(g1off + 126 * j + 47);

                    auto g1_xz_xxyyy = contrBuffer.data(g1off + 126 * j + 48);

                    auto g1_xz_xxyyz = contrBuffer.data(g1off + 126 * j + 49);

                    auto g1_xz_xxyzz = contrBuffer.data(g1off + 126 * j + 50);

                    auto g1_xz_xxzzz = contrBuffer.data(g1off + 126 * j + 51);

                    auto g1_xz_xyyyy = contrBuffer.data(g1off + 126 * j + 52);

                    auto g1_xz_xyyyz = contrBuffer.data(g1off + 126 * j + 53);

                    auto g1_xz_xyyzz = contrBuffer.data(g1off + 126 * j + 54);

                    auto g1_xz_xyzzz = contrBuffer.data(g1off + 126 * j + 55);

                    auto g1_xz_xzzzz = contrBuffer.data(g1off + 126 * j + 56);

                    auto g1_yy_xxxxx = contrBuffer.data(g1off + 126 * j + 63);

                    auto g1_yy_xxxxy = contrBuffer.data(g1off + 126 * j + 64);

                    auto g1_yy_xxxxz = contrBuffer.data(g1off + 126 * j + 65);

                    auto g1_yy_xxxyy = contrBuffer.data(g1off + 126 * j + 66);

                    auto g1_yy_xxxyz = contrBuffer.data(g1off + 126 * j + 67);

                    auto g1_yy_xxxzz = contrBuffer.data(g1off + 126 * j + 68);

                    auto g1_yy_xxyyy = contrBuffer.data(g1off + 126 * j + 69);

                    auto g1_yy_xxyyz = contrBuffer.data(g1off + 126 * j + 70);

                    auto g1_yy_xxyzz = contrBuffer.data(g1off + 126 * j + 71);

                    auto g1_yy_xxzzz = contrBuffer.data(g1off + 126 * j + 72);

                    auto g1_yy_xyyyy = contrBuffer.data(g1off + 126 * j + 73);

                    auto g1_yy_xyyyz = contrBuffer.data(g1off + 126 * j + 74);

                    auto g1_yy_xyyzz = contrBuffer.data(g1off + 126 * j + 75);

                    auto g1_yy_xyzzz = contrBuffer.data(g1off + 126 * j + 76);

                    auto g1_yy_xzzzz = contrBuffer.data(g1off + 126 * j + 77);

                    auto g1_yy_yyyyy = contrBuffer.data(g1off + 126 * j + 78);

                    auto g1_yy_yyyyz = contrBuffer.data(g1off + 126 * j + 79);

                    auto g1_yy_yyyzz = contrBuffer.data(g1off + 126 * j + 80);

                    auto g1_yy_yyzzz = contrBuffer.data(g1off + 126 * j + 81);

                    auto g1_yy_yzzzz = contrBuffer.data(g1off + 126 * j + 82);

                    auto g1_yz_xxxxx = contrBuffer.data(g1off + 126 * j + 84);

                    auto g1_yz_xxxxy = contrBuffer.data(g1off + 126 * j + 85);

                    auto g1_yz_xxxxz = contrBuffer.data(g1off + 126 * j + 86);

                    auto g1_yz_xxxyy = contrBuffer.data(g1off + 126 * j + 87);

                    auto g1_yz_xxxyz = contrBuffer.data(g1off + 126 * j + 88);

                    auto g1_yz_xxxzz = contrBuffer.data(g1off + 126 * j + 89);

                    auto g1_yz_xxyyy = contrBuffer.data(g1off + 126 * j + 90);

                    auto g1_yz_xxyyz = contrBuffer.data(g1off + 126 * j + 91);

                    auto g1_yz_xxyzz = contrBuffer.data(g1off + 126 * j + 92);

                    auto g1_yz_xxzzz = contrBuffer.data(g1off + 126 * j + 93);

                    auto g1_yz_xyyyy = contrBuffer.data(g1off + 126 * j + 94);

                    auto g1_yz_xyyyz = contrBuffer.data(g1off + 126 * j + 95);

                    auto g1_yz_xyyzz = contrBuffer.data(g1off + 126 * j + 96);

                    auto g1_yz_xyzzz = contrBuffer.data(g1off + 126 * j + 97);

                    auto g1_yz_xzzzz = contrBuffer.data(g1off + 126 * j + 98);

                    auto g1_yz_yyyyy = contrBuffer.data(g1off + 126 * j + 99);

                    auto g1_yz_yyyyz = contrBuffer.data(g1off + 126 * j + 100);

                    auto g1_yz_yyyzz = contrBuffer.data(g1off + 126 * j + 101);

                    auto g1_yz_yyzzz = contrBuffer.data(g1off + 126 * j + 102);

                    auto g1_yz_yzzzz = contrBuffer.data(g1off + 126 * j + 103);

                    auto g1_zz_xxxxx = contrBuffer.data(g1off + 126 * j + 105);

                    auto g1_zz_xxxxy = contrBuffer.data(g1off + 126 * j + 106);

                    auto g1_zz_xxxxz = contrBuffer.data(g1off + 126 * j + 107);

                    auto g1_zz_xxxyy = contrBuffer.data(g1off + 126 * j + 108);

                    auto g1_zz_xxxyz = contrBuffer.data(g1off + 126 * j + 109);

                    auto g1_zz_xxxzz = contrBuffer.data(g1off + 126 * j + 110);

                    auto g1_zz_xxyyy = contrBuffer.data(g1off + 126 * j + 111);

                    auto g1_zz_xxyyz = contrBuffer.data(g1off + 126 * j + 112);

                    auto g1_zz_xxyzz = contrBuffer.data(g1off + 126 * j + 113);

                    auto g1_zz_xxzzz = contrBuffer.data(g1off + 126 * j + 114);

                    auto g1_zz_xyyyy = contrBuffer.data(g1off + 126 * j + 115);

                    auto g1_zz_xyyyz = contrBuffer.data(g1off + 126 * j + 116);

                    auto g1_zz_xyyzz = contrBuffer.data(g1off + 126 * j + 117);

                    auto g1_zz_xyzzz = contrBuffer.data(g1off + 126 * j + 118);

                    auto g1_zz_xzzzz = contrBuffer.data(g1off + 126 * j + 119);

                    auto g1_zz_yyyyy = contrBuffer.data(g1off + 126 * j + 120);

                    auto g1_zz_yyyyz = contrBuffer.data(g1off + 126 * j + 121);

                    auto g1_zz_yyyzz = contrBuffer.data(g1off + 126 * j + 122);

                    auto g1_zz_yyzzz = contrBuffer.data(g1off + 126 * j + 123);

                    auto g1_zz_yzzzz = contrBuffer.data(g1off + 126 * j + 124);

                    auto g1_zz_zzzzz = contrBuffer.data(g1off + 126 * j + 125);

                    // set up pointers to (SX|g(r,r')|FG)^(m) integrals

                    auto g_xxx_xxxx = contrBuffer.data(goff + 150 * j);

                    auto g_xxx_xxxy = contrBuffer.data(goff + 150 * j + 1);

                    auto g_xxx_xxxz = contrBuffer.data(goff + 150 * j + 2);

                    auto g_xxx_xxyy = contrBuffer.data(goff + 150 * j + 3);

                    auto g_xxx_xxyz = contrBuffer.data(goff + 150 * j + 4);

                    auto g_xxx_xxzz = contrBuffer.data(goff + 150 * j + 5);

                    auto g_xxx_xyyy = contrBuffer.data(goff + 150 * j + 6);

                    auto g_xxx_xyyz = contrBuffer.data(goff + 150 * j + 7);

                    auto g_xxx_xyzz = contrBuffer.data(goff + 150 * j + 8);

                    auto g_xxx_xzzz = contrBuffer.data(goff + 150 * j + 9);

                    auto g_xxx_yyyy = contrBuffer.data(goff + 150 * j + 10);

                    auto g_xxx_yyyz = contrBuffer.data(goff + 150 * j + 11);

                    auto g_xxx_yyzz = contrBuffer.data(goff + 150 * j + 12);

                    auto g_xxx_yzzz = contrBuffer.data(goff + 150 * j + 13);

                    auto g_xxx_zzzz = contrBuffer.data(goff + 150 * j + 14);

                    auto g_xxy_xxxx = contrBuffer.data(goff + 150 * j + 15);

                    auto g_xxy_xxxy = contrBuffer.data(goff + 150 * j + 16);

                    auto g_xxy_xxxz = contrBuffer.data(goff + 150 * j + 17);

                    auto g_xxy_xxyy = contrBuffer.data(goff + 150 * j + 18);

                    auto g_xxy_xxyz = contrBuffer.data(goff + 150 * j + 19);

                    auto g_xxy_xxzz = contrBuffer.data(goff + 150 * j + 20);

                    auto g_xxy_xyyy = contrBuffer.data(goff + 150 * j + 21);

                    auto g_xxy_xyyz = contrBuffer.data(goff + 150 * j + 22);

                    auto g_xxy_xyzz = contrBuffer.data(goff + 150 * j + 23);

                    auto g_xxy_xzzz = contrBuffer.data(goff + 150 * j + 24);

                    auto g_xxy_yyyy = contrBuffer.data(goff + 150 * j + 25);

                    auto g_xxy_yyyz = contrBuffer.data(goff + 150 * j + 26);

                    auto g_xxy_yyzz = contrBuffer.data(goff + 150 * j + 27);

                    auto g_xxy_yzzz = contrBuffer.data(goff + 150 * j + 28);

                    auto g_xxy_zzzz = contrBuffer.data(goff + 150 * j + 29);

                    auto g_xxz_xxxx = contrBuffer.data(goff + 150 * j + 30);

                    auto g_xxz_xxxy = contrBuffer.data(goff + 150 * j + 31);

                    auto g_xxz_xxxz = contrBuffer.data(goff + 150 * j + 32);

                    auto g_xxz_xxyy = contrBuffer.data(goff + 150 * j + 33);

                    auto g_xxz_xxyz = contrBuffer.data(goff + 150 * j + 34);

                    auto g_xxz_xxzz = contrBuffer.data(goff + 150 * j + 35);

                    auto g_xxz_xyyy = contrBuffer.data(goff + 150 * j + 36);

                    auto g_xxz_xyyz = contrBuffer.data(goff + 150 * j + 37);

                    auto g_xxz_xyzz = contrBuffer.data(goff + 150 * j + 38);

                    auto g_xxz_xzzz = contrBuffer.data(goff + 150 * j + 39);

                    auto g_xxz_yyyy = contrBuffer.data(goff + 150 * j + 40);

                    auto g_xxz_yyyz = contrBuffer.data(goff + 150 * j + 41);

                    auto g_xxz_yyzz = contrBuffer.data(goff + 150 * j + 42);

                    auto g_xxz_yzzz = contrBuffer.data(goff + 150 * j + 43);

                    auto g_xxz_zzzz = contrBuffer.data(goff + 150 * j + 44);

                    auto g_xyy_xxxx = contrBuffer.data(goff + 150 * j + 45);

                    auto g_xyy_xxxy = contrBuffer.data(goff + 150 * j + 46);

                    auto g_xyy_xxxz = contrBuffer.data(goff + 150 * j + 47);

                    auto g_xyy_xxyy = contrBuffer.data(goff + 150 * j + 48);

                    auto g_xyy_xxyz = contrBuffer.data(goff + 150 * j + 49);

                    auto g_xyy_xxzz = contrBuffer.data(goff + 150 * j + 50);

                    auto g_xyy_xyyy = contrBuffer.data(goff + 150 * j + 51);

                    auto g_xyy_xyyz = contrBuffer.data(goff + 150 * j + 52);

                    auto g_xyy_xyzz = contrBuffer.data(goff + 150 * j + 53);

                    auto g_xyy_xzzz = contrBuffer.data(goff + 150 * j + 54);

                    auto g_xyy_yyyy = contrBuffer.data(goff + 150 * j + 55);

                    auto g_xyy_yyyz = contrBuffer.data(goff + 150 * j + 56);

                    auto g_xyy_yyzz = contrBuffer.data(goff + 150 * j + 57);

                    auto g_xyy_yzzz = contrBuffer.data(goff + 150 * j + 58);

                    auto g_xyy_zzzz = contrBuffer.data(goff + 150 * j + 59);

                    auto g_xyz_xxxx = contrBuffer.data(goff + 150 * j + 60);

                    auto g_xyz_xxxy = contrBuffer.data(goff + 150 * j + 61);

                    auto g_xyz_xxxz = contrBuffer.data(goff + 150 * j + 62);

                    auto g_xyz_xxyy = contrBuffer.data(goff + 150 * j + 63);

                    auto g_xyz_xxyz = contrBuffer.data(goff + 150 * j + 64);

                    auto g_xyz_xxzz = contrBuffer.data(goff + 150 * j + 65);

                    auto g_xyz_xyyy = contrBuffer.data(goff + 150 * j + 66);

                    auto g_xyz_xyyz = contrBuffer.data(goff + 150 * j + 67);

                    auto g_xyz_xyzz = contrBuffer.data(goff + 150 * j + 68);

                    auto g_xyz_xzzz = contrBuffer.data(goff + 150 * j + 69);

                    auto g_xyz_yyyy = contrBuffer.data(goff + 150 * j + 70);

                    auto g_xyz_yyyz = contrBuffer.data(goff + 150 * j + 71);

                    auto g_xyz_yyzz = contrBuffer.data(goff + 150 * j + 72);

                    auto g_xyz_yzzz = contrBuffer.data(goff + 150 * j + 73);

                    auto g_xyz_zzzz = contrBuffer.data(goff + 150 * j + 74);

                    auto g_xzz_xxxx = contrBuffer.data(goff + 150 * j + 75);

                    auto g_xzz_xxxy = contrBuffer.data(goff + 150 * j + 76);

                    auto g_xzz_xxxz = contrBuffer.data(goff + 150 * j + 77);

                    auto g_xzz_xxyy = contrBuffer.data(goff + 150 * j + 78);

                    auto g_xzz_xxyz = contrBuffer.data(goff + 150 * j + 79);

                    auto g_xzz_xxzz = contrBuffer.data(goff + 150 * j + 80);

                    auto g_xzz_xyyy = contrBuffer.data(goff + 150 * j + 81);

                    auto g_xzz_xyyz = contrBuffer.data(goff + 150 * j + 82);

                    auto g_xzz_xyzz = contrBuffer.data(goff + 150 * j + 83);

                    auto g_xzz_xzzz = contrBuffer.data(goff + 150 * j + 84);

                    auto g_xzz_yyyy = contrBuffer.data(goff + 150 * j + 85);

                    auto g_xzz_yyyz = contrBuffer.data(goff + 150 * j + 86);

                    auto g_xzz_yyzz = contrBuffer.data(goff + 150 * j + 87);

                    auto g_xzz_yzzz = contrBuffer.data(goff + 150 * j + 88);

                    auto g_xzz_zzzz = contrBuffer.data(goff + 150 * j + 89);

                    auto g_yyy_xxxx = contrBuffer.data(goff + 150 * j + 90);

                    auto g_yyy_xxxy = contrBuffer.data(goff + 150 * j + 91);

                    auto g_yyy_xxxz = contrBuffer.data(goff + 150 * j + 92);

                    auto g_yyy_xxyy = contrBuffer.data(goff + 150 * j + 93);

                    auto g_yyy_xxyz = contrBuffer.data(goff + 150 * j + 94);

                    auto g_yyy_xxzz = contrBuffer.data(goff + 150 * j + 95);

                    auto g_yyy_xyyy = contrBuffer.data(goff + 150 * j + 96);

                    auto g_yyy_xyyz = contrBuffer.data(goff + 150 * j + 97);

                    auto g_yyy_xyzz = contrBuffer.data(goff + 150 * j + 98);

                    auto g_yyy_xzzz = contrBuffer.data(goff + 150 * j + 99);

                    auto g_yyy_yyyy = contrBuffer.data(goff + 150 * j + 100);

                    auto g_yyy_yyyz = contrBuffer.data(goff + 150 * j + 101);

                    auto g_yyy_yyzz = contrBuffer.data(goff + 150 * j + 102);

                    auto g_yyy_yzzz = contrBuffer.data(goff + 150 * j + 103);

                    auto g_yyy_zzzz = contrBuffer.data(goff + 150 * j + 104);

                    auto g_yyz_xxxx = contrBuffer.data(goff + 150 * j + 105);

                    auto g_yyz_xxxy = contrBuffer.data(goff + 150 * j + 106);

                    auto g_yyz_xxxz = contrBuffer.data(goff + 150 * j + 107);

                    auto g_yyz_xxyy = contrBuffer.data(goff + 150 * j + 108);

                    auto g_yyz_xxyz = contrBuffer.data(goff + 150 * j + 109);

                    auto g_yyz_xxzz = contrBuffer.data(goff + 150 * j + 110);

                    auto g_yyz_xyyy = contrBuffer.data(goff + 150 * j + 111);

                    auto g_yyz_xyyz = contrBuffer.data(goff + 150 * j + 112);

                    auto g_yyz_xyzz = contrBuffer.data(goff + 150 * j + 113);

                    auto g_yyz_xzzz = contrBuffer.data(goff + 150 * j + 114);

                    auto g_yyz_yyyy = contrBuffer.data(goff + 150 * j + 115);

                    auto g_yyz_yyyz = contrBuffer.data(goff + 150 * j + 116);

                    auto g_yyz_yyzz = contrBuffer.data(goff + 150 * j + 117);

                    auto g_yyz_yzzz = contrBuffer.data(goff + 150 * j + 118);

                    auto g_yyz_zzzz = contrBuffer.data(goff + 150 * j + 119);

                    auto g_yzz_xxxx = contrBuffer.data(goff + 150 * j + 120);

                    auto g_yzz_xxxy = contrBuffer.data(goff + 150 * j + 121);

                    auto g_yzz_xxxz = contrBuffer.data(goff + 150 * j + 122);

                    auto g_yzz_xxyy = contrBuffer.data(goff + 150 * j + 123);

                    auto g_yzz_xxyz = contrBuffer.data(goff + 150 * j + 124);

                    auto g_yzz_xxzz = contrBuffer.data(goff + 150 * j + 125);

                    auto g_yzz_xyyy = contrBuffer.data(goff + 150 * j + 126);

                    auto g_yzz_xyyz = contrBuffer.data(goff + 150 * j + 127);

                    auto g_yzz_xyzz = contrBuffer.data(goff + 150 * j + 128);

                    auto g_yzz_xzzz = contrBuffer.data(goff + 150 * j + 129);

                    auto g_yzz_yyyy = contrBuffer.data(goff + 150 * j + 130);

                    auto g_yzz_yyyz = contrBuffer.data(goff + 150 * j + 131);

                    auto g_yzz_yyzz = contrBuffer.data(goff + 150 * j + 132);

                    auto g_yzz_yzzz = contrBuffer.data(goff + 150 * j + 133);

                    auto g_yzz_zzzz = contrBuffer.data(goff + 150 * j + 134);

                    auto g_zzz_xxxx = contrBuffer.data(goff + 150 * j + 135);

                    auto g_zzz_xxxy = contrBuffer.data(goff + 150 * j + 136);

                    auto g_zzz_xxxz = contrBuffer.data(goff + 150 * j + 137);

                    auto g_zzz_xxyy = contrBuffer.data(goff + 150 * j + 138);

                    auto g_zzz_xxyz = contrBuffer.data(goff + 150 * j + 139);

                    auto g_zzz_xxzz = contrBuffer.data(goff + 150 * j + 140);

                    auto g_zzz_xyyy = contrBuffer.data(goff + 150 * j + 141);

                    auto g_zzz_xyyz = contrBuffer.data(goff + 150 * j + 142);

                    auto g_zzz_xyzz = contrBuffer.data(goff + 150 * j + 143);

                    auto g_zzz_xzzz = contrBuffer.data(goff + 150 * j + 144);

                    auto g_zzz_yyyy = contrBuffer.data(goff + 150 * j + 145);

                    auto g_zzz_yyyz = contrBuffer.data(goff + 150 * j + 146);

                    auto g_zzz_yyzz = contrBuffer.data(goff + 150 * j + 147);

                    auto g_zzz_yzzz = contrBuffer.data(goff + 150 * j + 148);

                    auto g_zzz_zzzz = contrBuffer.data(goff + 150 * j + 149);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxxx, g2_xx_xxxy,\
                                             g2_xx_xxxz, g2_xx_xxyy, g2_xx_xxyz,\
                                             g2_xx_xxzz, g2_xx_xyyy, g2_xx_xyyz,\
                                             g2_xx_xyzz, g2_xx_xzzz, g2_xx_yyyy,\
                                             g2_xx_yyyz, g2_xx_yyzz, g2_xx_yzzz,\
                                             g2_xx_zzzz, g2_xy_xxxx, g2_xy_xxxy,\
                                             g2_xy_xxxz, g2_xy_xxyy, g2_xy_xxyz,\
                                             g2_xy_xxzz, g2_xy_xyyy, g2_xy_xyyz,\
                                             g2_xy_xyzz, g2_xy_xzzz, g2_xy_yyyy,\
                                             g2_xy_yyyz, g2_xy_yyzz, g2_xy_yzzz,\
                                             g2_xy_zzzz, g2_xz_xxxx, g2_xz_xxxy,\
                                             g2_xz_xxxz, g2_xz_xxyy, g2_xz_xxyz,\
                                             g2_xz_xxzz, g2_xz_xyyy, g2_xz_xyyz,\
                                             g2_xz_xyzz, g2_xz_xzzz, g2_xz_yyyy,\
                                             g2_xz_yyyz, g2_xz_yyzz, g2_xz_yzzz,\
                                             g2_xz_zzzz, g2_yy_xxxx, g2_yy_xxxy,\
                                             g2_yy_xxxz, g2_yy_xxyy, g2_yy_xxyz,\
                                             g2_yy_xxzz, g2_yy_xyyy, g2_yy_xyyz,\
                                             g2_yy_xyzz, g2_yy_xzzz, g2_yy_yyyy,\
                                             g2_yy_yyyz, g2_yy_yyzz, g2_yy_yzzz,\
                                             g2_yy_zzzz, g2_yz_xxxx, g2_yz_xxxy,\
                                             g2_yz_xxxz, g2_yz_xxyy, g2_yz_xxyz,\
                                             g2_yz_xxzz, g2_yz_xyyy, g2_yz_xyyz,\
                                             g2_yz_xyzz, g2_yz_xzzz, g2_yz_yyyy,\
                                             g2_yz_yyyz, g2_yz_yyzz, g2_yz_yzzz,\
                                             g2_yz_zzzz, g2_zz_xxxx, g2_zz_xxxy,\
                                             g2_zz_xxxz, g2_zz_xxyy, g2_zz_xxyz,\
                                             g2_zz_xxzz, g2_zz_xyyy, g2_zz_xyyz,\
                                             g2_zz_xyzz, g2_zz_xzzz, g2_zz_yyyy,\
                                             g2_zz_yyyz, g2_zz_yyzz, g2_zz_yzzz,\
                                             g2_zz_zzzz, g1_xx_xxxxx, g1_xx_xxxxy,\
                                             g1_xx_xxxxz, g1_xx_xxxyy, g1_xx_xxxyz,\
                                             g1_xx_xxxzz, g1_xx_xxyyy, g1_xx_xxyyz,\
                                             g1_xx_xxyzz, g1_xx_xxzzz, g1_xx_xyyyy,\
                                             g1_xx_xyyyz, g1_xx_xyyzz, g1_xx_xyzzz,\
                                             g1_xx_xzzzz, g1_xy_xxxxx, g1_xy_xxxxy,\
                                             g1_xy_xxxxz, g1_xy_xxxyy, g1_xy_xxxyz,\
                                             g1_xy_xxxzz, g1_xy_xxyyy, g1_xy_xxyyz,\
                                             g1_xy_xxyzz, g1_xy_xxzzz, g1_xy_xyyyy,\
                                             g1_xy_xyyyz, g1_xy_xyyzz, g1_xy_xyzzz,\
                                             g1_xy_xzzzz, g1_xz_xxxxx, g1_xz_xxxxy,\
                                             g1_xz_xxxxz, g1_xz_xxxyy, g1_xz_xxxyz,\
                                             g1_xz_xxxzz, g1_xz_xxyyy, g1_xz_xxyyz,\
                                             g1_xz_xxyzz, g1_xz_xxzzz, g1_xz_xyyyy,\
                                             g1_xz_xyyyz, g1_xz_xyyzz, g1_xz_xyzzz,\
                                             g1_xz_xzzzz, g1_yy_xxxxx, g1_yy_xxxxy,\
                                             g1_yy_xxxxz, g1_yy_xxxyy, g1_yy_xxxyz,\
                                             g1_yy_xxxzz, g1_yy_xxyyy, g1_yy_xxyyz,\
                                             g1_yy_xxyzz, g1_yy_xxzzz, g1_yy_xyyyy,\
                                             g1_yy_xyyyz, g1_yy_xyyzz, g1_yy_xyzzz,\
                                             g1_yy_xzzzz, g1_yy_yyyyy, g1_yy_yyyyz,\
                                             g1_yy_yyyzz, g1_yy_yyzzz, g1_yy_yzzzz,\
                                             g1_yz_xxxxx, g1_yz_xxxxy,\
                                             g1_yz_xxxxz, g1_yz_xxxyy, g1_yz_xxxyz,\
                                             g1_yz_xxxzz, g1_yz_xxyyy, g1_yz_xxyyz,\
                                             g1_yz_xxyzz, g1_yz_xxzzz, g1_yz_xyyyy,\
                                             g1_yz_xyyyz, g1_yz_xyyzz, g1_yz_xyzzz,\
                                             g1_yz_xzzzz, g1_yz_yyyyy, g1_yz_yyyyz,\
                                             g1_yz_yyyzz, g1_yz_yyzzz, g1_yz_yzzzz,\
                                             g1_zz_xxxxx, g1_zz_xxxxy,\
                                             g1_zz_xxxxz, g1_zz_xxxyy, g1_zz_xxxyz,\
                                             g1_zz_xxxzz, g1_zz_xxyyy, g1_zz_xxyyz,\
                                             g1_zz_xxyzz, g1_zz_xxzzz, g1_zz_xyyyy,\
                                             g1_zz_xyyyz, g1_zz_xyyzz, g1_zz_xyzzz,\
                                             g1_zz_xzzzz, g1_zz_yyyyy, g1_zz_yyyyz,\
                                             g1_zz_yyyzz, g1_zz_yyzzz, g1_zz_yzzzz,\
                                             g1_zz_zzzzz, g_xxx_xxxx, g_xxx_xxxy,\
                                             g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz,\
                                             g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz,\
                                             g_xxx_xyzz, g_xxx_xzzz, g_xxx_yyyy,\
                                             g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz,\
                                             g_xxx_zzzz, g_xxy_xxxx, g_xxy_xxxy,\
                                             g_xxy_xxxz, g_xxy_xxyy, g_xxy_xxyz,\
                                             g_xxy_xxzz, g_xxy_xyyy, g_xxy_xyyz,\
                                             g_xxy_xyzz, g_xxy_xzzz, g_xxy_yyyy,\
                                             g_xxy_yyyz, g_xxy_yyzz, g_xxy_yzzz,\
                                             g_xxy_zzzz, g_xxz_xxxx, g_xxz_xxxy,\
                                             g_xxz_xxxz, g_xxz_xxyy, g_xxz_xxyz,\
                                             g_xxz_xxzz, g_xxz_xyyy, g_xxz_xyyz,\
                                             g_xxz_xyzz, g_xxz_xzzz, g_xxz_yyyy,\
                                             g_xxz_yyyz, g_xxz_yyzz, g_xxz_yzzz,\
                                             g_xxz_zzzz, g_xyy_xxxx, g_xyy_xxxy,\
                                             g_xyy_xxxz, g_xyy_xxyy, g_xyy_xxyz,\
                                             g_xyy_xxzz, g_xyy_xyyy, g_xyy_xyyz,\
                                             g_xyy_xyzz, g_xyy_xzzz, g_xyy_yyyy,\
                                             g_xyy_yyyz, g_xyy_yyzz, g_xyy_yzzz,\
                                             g_xyy_zzzz, g_xyz_xxxx, g_xyz_xxxy,\
                                             g_xyz_xxxz, g_xyz_xxyy, g_xyz_xxyz,\
                                             g_xyz_xxzz, g_xyz_xyyy, g_xyz_xyyz,\
                                             g_xyz_xyzz, g_xyz_xzzz, g_xyz_yyyy,\
                                             g_xyz_yyyz, g_xyz_yyzz, g_xyz_yzzz,\
                                             g_xyz_zzzz, g_xzz_xxxx, g_xzz_xxxy,\
                                             g_xzz_xxxz, g_xzz_xxyy, g_xzz_xxyz,\
                                             g_xzz_xxzz, g_xzz_xyyy, g_xzz_xyyz,\
                                             g_xzz_xyzz, g_xzz_xzzz, g_xzz_yyyy,\
                                             g_xzz_yyyz, g_xzz_yyzz, g_xzz_yzzz,\
                                             g_xzz_zzzz, g_yyy_xxxx, g_yyy_xxxy,\
                                             g_yyy_xxxz, g_yyy_xxyy, g_yyy_xxyz,\
                                             g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz,\
                                             g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy,\
                                             g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz,\
                                             g_yyy_zzzz, g_yyz_xxxx, g_yyz_xxxy,\
                                             g_yyz_xxxz, g_yyz_xxyy, g_yyz_xxyz,\
                                             g_yyz_xxzz, g_yyz_xyyy, g_yyz_xyyz,\
                                             g_yyz_xyzz, g_yyz_xzzz, g_yyz_yyyy,\
                                             g_yyz_yyyz, g_yyz_yyzz, g_yyz_yzzz,\
                                             g_yyz_zzzz, g_yzz_xxxx, g_yzz_xxxy,\
                                             g_yzz_xxxz, g_yzz_xxyy, g_yzz_xxyz,\
                                             g_yzz_xxzz, g_yzz_xyyy, g_yzz_xyyz,\
                                             g_yzz_xyzz, g_yzz_xzzz, g_yzz_yyyy,\
                                             g_yzz_yyyz, g_yzz_yyzz, g_yzz_yzzz,\
                                             g_yzz_zzzz, g_zzz_xxxx, g_zzz_xxxy,\
                                             g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz,\
                                             g_zzz_xxzz, g_zzz_xyyy, g_zzz_xyyz,\
                                             g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy,\
                                             g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz,\
                                             g_zzz_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xxx_xxxx[k] = g1_xx_xxxxx[k] - fr * g2_xx_xxxx[k];

                        g_xxx_xxxy[k] = g1_xx_xxxxy[k] - fr * g2_xx_xxxy[k];

                        g_xxx_xxxz[k] = g1_xx_xxxxz[k] - fr * g2_xx_xxxz[k];

                        g_xxx_xxyy[k] = g1_xx_xxxyy[k] - fr * g2_xx_xxyy[k];

                        g_xxx_xxyz[k] = g1_xx_xxxyz[k] - fr * g2_xx_xxyz[k];

                        g_xxx_xxzz[k] = g1_xx_xxxzz[k] - fr * g2_xx_xxzz[k];

                        g_xxx_xyyy[k] = g1_xx_xxyyy[k] - fr * g2_xx_xyyy[k];

                        g_xxx_xyyz[k] = g1_xx_xxyyz[k] - fr * g2_xx_xyyz[k];

                        g_xxx_xyzz[k] = g1_xx_xxyzz[k] - fr * g2_xx_xyzz[k];

                        g_xxx_xzzz[k] = g1_xx_xxzzz[k] - fr * g2_xx_xzzz[k];

                        g_xxx_yyyy[k] = g1_xx_xyyyy[k] - fr * g2_xx_yyyy[k];

                        g_xxx_yyyz[k] = g1_xx_xyyyz[k] - fr * g2_xx_yyyz[k];

                        g_xxx_yyzz[k] = g1_xx_xyyzz[k] - fr * g2_xx_yyzz[k];

                        g_xxx_yzzz[k] = g1_xx_xyzzz[k] - fr * g2_xx_yzzz[k];

                        g_xxx_zzzz[k] = g1_xx_xzzzz[k] - fr * g2_xx_zzzz[k];

                        g_xxy_xxxx[k] = g1_xy_xxxxx[k] - fr * g2_xy_xxxx[k];

                        g_xxy_xxxy[k] = g1_xy_xxxxy[k] - fr * g2_xy_xxxy[k];

                        g_xxy_xxxz[k] = g1_xy_xxxxz[k] - fr * g2_xy_xxxz[k];

                        g_xxy_xxyy[k] = g1_xy_xxxyy[k] - fr * g2_xy_xxyy[k];

                        g_xxy_xxyz[k] = g1_xy_xxxyz[k] - fr * g2_xy_xxyz[k];

                        g_xxy_xxzz[k] = g1_xy_xxxzz[k] - fr * g2_xy_xxzz[k];

                        g_xxy_xyyy[k] = g1_xy_xxyyy[k] - fr * g2_xy_xyyy[k];

                        g_xxy_xyyz[k] = g1_xy_xxyyz[k] - fr * g2_xy_xyyz[k];

                        g_xxy_xyzz[k] = g1_xy_xxyzz[k] - fr * g2_xy_xyzz[k];

                        g_xxy_xzzz[k] = g1_xy_xxzzz[k] - fr * g2_xy_xzzz[k];

                        g_xxy_yyyy[k] = g1_xy_xyyyy[k] - fr * g2_xy_yyyy[k];

                        g_xxy_yyyz[k] = g1_xy_xyyyz[k] - fr * g2_xy_yyyz[k];

                        g_xxy_yyzz[k] = g1_xy_xyyzz[k] - fr * g2_xy_yyzz[k];

                        g_xxy_yzzz[k] = g1_xy_xyzzz[k] - fr * g2_xy_yzzz[k];

                        g_xxy_zzzz[k] = g1_xy_xzzzz[k] - fr * g2_xy_zzzz[k];

                        g_xxz_xxxx[k] = g1_xz_xxxxx[k] - fr * g2_xz_xxxx[k];

                        g_xxz_xxxy[k] = g1_xz_xxxxy[k] - fr * g2_xz_xxxy[k];

                        g_xxz_xxxz[k] = g1_xz_xxxxz[k] - fr * g2_xz_xxxz[k];

                        g_xxz_xxyy[k] = g1_xz_xxxyy[k] - fr * g2_xz_xxyy[k];

                        g_xxz_xxyz[k] = g1_xz_xxxyz[k] - fr * g2_xz_xxyz[k];

                        g_xxz_xxzz[k] = g1_xz_xxxzz[k] - fr * g2_xz_xxzz[k];

                        g_xxz_xyyy[k] = g1_xz_xxyyy[k] - fr * g2_xz_xyyy[k];

                        g_xxz_xyyz[k] = g1_xz_xxyyz[k] - fr * g2_xz_xyyz[k];

                        g_xxz_xyzz[k] = g1_xz_xxyzz[k] - fr * g2_xz_xyzz[k];

                        g_xxz_xzzz[k] = g1_xz_xxzzz[k] - fr * g2_xz_xzzz[k];

                        g_xxz_yyyy[k] = g1_xz_xyyyy[k] - fr * g2_xz_yyyy[k];

                        g_xxz_yyyz[k] = g1_xz_xyyyz[k] - fr * g2_xz_yyyz[k];

                        g_xxz_yyzz[k] = g1_xz_xyyzz[k] - fr * g2_xz_yyzz[k];

                        g_xxz_yzzz[k] = g1_xz_xyzzz[k] - fr * g2_xz_yzzz[k];

                        g_xxz_zzzz[k] = g1_xz_xzzzz[k] - fr * g2_xz_zzzz[k];

                        g_xyy_xxxx[k] = g1_yy_xxxxx[k] - fr * g2_yy_xxxx[k];

                        g_xyy_xxxy[k] = g1_yy_xxxxy[k] - fr * g2_yy_xxxy[k];

                        g_xyy_xxxz[k] = g1_yy_xxxxz[k] - fr * g2_yy_xxxz[k];

                        g_xyy_xxyy[k] = g1_yy_xxxyy[k] - fr * g2_yy_xxyy[k];

                        g_xyy_xxyz[k] = g1_yy_xxxyz[k] - fr * g2_yy_xxyz[k];

                        g_xyy_xxzz[k] = g1_yy_xxxzz[k] - fr * g2_yy_xxzz[k];

                        g_xyy_xyyy[k] = g1_yy_xxyyy[k] - fr * g2_yy_xyyy[k];

                        g_xyy_xyyz[k] = g1_yy_xxyyz[k] - fr * g2_yy_xyyz[k];

                        g_xyy_xyzz[k] = g1_yy_xxyzz[k] - fr * g2_yy_xyzz[k];

                        g_xyy_xzzz[k] = g1_yy_xxzzz[k] - fr * g2_yy_xzzz[k];

                        g_xyy_yyyy[k] = g1_yy_xyyyy[k] - fr * g2_yy_yyyy[k];

                        g_xyy_yyyz[k] = g1_yy_xyyyz[k] - fr * g2_yy_yyyz[k];

                        g_xyy_yyzz[k] = g1_yy_xyyzz[k] - fr * g2_yy_yyzz[k];

                        g_xyy_yzzz[k] = g1_yy_xyzzz[k] - fr * g2_yy_yzzz[k];

                        g_xyy_zzzz[k] = g1_yy_xzzzz[k] - fr * g2_yy_zzzz[k];

                        g_xyz_xxxx[k] = g1_yz_xxxxx[k] - fr * g2_yz_xxxx[k];

                        g_xyz_xxxy[k] = g1_yz_xxxxy[k] - fr * g2_yz_xxxy[k];

                        g_xyz_xxxz[k] = g1_yz_xxxxz[k] - fr * g2_yz_xxxz[k];

                        g_xyz_xxyy[k] = g1_yz_xxxyy[k] - fr * g2_yz_xxyy[k];

                        g_xyz_xxyz[k] = g1_yz_xxxyz[k] - fr * g2_yz_xxyz[k];

                        g_xyz_xxzz[k] = g1_yz_xxxzz[k] - fr * g2_yz_xxzz[k];

                        g_xyz_xyyy[k] = g1_yz_xxyyy[k] - fr * g2_yz_xyyy[k];

                        g_xyz_xyyz[k] = g1_yz_xxyyz[k] - fr * g2_yz_xyyz[k];

                        g_xyz_xyzz[k] = g1_yz_xxyzz[k] - fr * g2_yz_xyzz[k];

                        g_xyz_xzzz[k] = g1_yz_xxzzz[k] - fr * g2_yz_xzzz[k];

                        g_xyz_yyyy[k] = g1_yz_xyyyy[k] - fr * g2_yz_yyyy[k];

                        g_xyz_yyyz[k] = g1_yz_xyyyz[k] - fr * g2_yz_yyyz[k];

                        g_xyz_yyzz[k] = g1_yz_xyyzz[k] - fr * g2_yz_yyzz[k];

                        g_xyz_yzzz[k] = g1_yz_xyzzz[k] - fr * g2_yz_yzzz[k];

                        g_xyz_zzzz[k] = g1_yz_xzzzz[k] - fr * g2_yz_zzzz[k];

                        g_xzz_xxxx[k] = g1_zz_xxxxx[k] - fr * g2_zz_xxxx[k];

                        g_xzz_xxxy[k] = g1_zz_xxxxy[k] - fr * g2_zz_xxxy[k];

                        g_xzz_xxxz[k] = g1_zz_xxxxz[k] - fr * g2_zz_xxxz[k];

                        g_xzz_xxyy[k] = g1_zz_xxxyy[k] - fr * g2_zz_xxyy[k];

                        g_xzz_xxyz[k] = g1_zz_xxxyz[k] - fr * g2_zz_xxyz[k];

                        g_xzz_xxzz[k] = g1_zz_xxxzz[k] - fr * g2_zz_xxzz[k];

                        g_xzz_xyyy[k] = g1_zz_xxyyy[k] - fr * g2_zz_xyyy[k];

                        g_xzz_xyyz[k] = g1_zz_xxyyz[k] - fr * g2_zz_xyyz[k];

                        g_xzz_xyzz[k] = g1_zz_xxyzz[k] - fr * g2_zz_xyzz[k];

                        g_xzz_xzzz[k] = g1_zz_xxzzz[k] - fr * g2_zz_xzzz[k];

                        g_xzz_yyyy[k] = g1_zz_xyyyy[k] - fr * g2_zz_yyyy[k];

                        g_xzz_yyyz[k] = g1_zz_xyyyz[k] - fr * g2_zz_yyyz[k];

                        g_xzz_yyzz[k] = g1_zz_xyyzz[k] - fr * g2_zz_yyzz[k];

                        g_xzz_yzzz[k] = g1_zz_xyzzz[k] - fr * g2_zz_yzzz[k];

                        g_xzz_zzzz[k] = g1_zz_xzzzz[k] - fr * g2_zz_zzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yyy_xxxx[k] = g1_yy_xxxxy[k] - fr * g2_yy_xxxx[k];

                        g_yyy_xxxy[k] = g1_yy_xxxyy[k] - fr * g2_yy_xxxy[k];

                        g_yyy_xxxz[k] = g1_yy_xxxyz[k] - fr * g2_yy_xxxz[k];

                        g_yyy_xxyy[k] = g1_yy_xxyyy[k] - fr * g2_yy_xxyy[k];

                        g_yyy_xxyz[k] = g1_yy_xxyyz[k] - fr * g2_yy_xxyz[k];

                        g_yyy_xxzz[k] = g1_yy_xxyzz[k] - fr * g2_yy_xxzz[k];

                        g_yyy_xyyy[k] = g1_yy_xyyyy[k] - fr * g2_yy_xyyy[k];

                        g_yyy_xyyz[k] = g1_yy_xyyyz[k] - fr * g2_yy_xyyz[k];

                        g_yyy_xyzz[k] = g1_yy_xyyzz[k] - fr * g2_yy_xyzz[k];

                        g_yyy_xzzz[k] = g1_yy_xyzzz[k] - fr * g2_yy_xzzz[k];

                        g_yyy_yyyy[k] = g1_yy_yyyyy[k] - fr * g2_yy_yyyy[k];

                        g_yyy_yyyz[k] = g1_yy_yyyyz[k] - fr * g2_yy_yyyz[k];

                        g_yyy_yyzz[k] = g1_yy_yyyzz[k] - fr * g2_yy_yyzz[k];

                        g_yyy_yzzz[k] = g1_yy_yyzzz[k] - fr * g2_yy_yzzz[k];

                        g_yyy_zzzz[k] = g1_yy_yzzzz[k] - fr * g2_yy_zzzz[k];

                        g_yyz_xxxx[k] = g1_yz_xxxxy[k] - fr * g2_yz_xxxx[k];

                        g_yyz_xxxy[k] = g1_yz_xxxyy[k] - fr * g2_yz_xxxy[k];

                        g_yyz_xxxz[k] = g1_yz_xxxyz[k] - fr * g2_yz_xxxz[k];

                        g_yyz_xxyy[k] = g1_yz_xxyyy[k] - fr * g2_yz_xxyy[k];

                        g_yyz_xxyz[k] = g1_yz_xxyyz[k] - fr * g2_yz_xxyz[k];

                        g_yyz_xxzz[k] = g1_yz_xxyzz[k] - fr * g2_yz_xxzz[k];

                        g_yyz_xyyy[k] = g1_yz_xyyyy[k] - fr * g2_yz_xyyy[k];

                        g_yyz_xyyz[k] = g1_yz_xyyyz[k] - fr * g2_yz_xyyz[k];

                        g_yyz_xyzz[k] = g1_yz_xyyzz[k] - fr * g2_yz_xyzz[k];

                        g_yyz_xzzz[k] = g1_yz_xyzzz[k] - fr * g2_yz_xzzz[k];

                        g_yyz_yyyy[k] = g1_yz_yyyyy[k] - fr * g2_yz_yyyy[k];

                        g_yyz_yyyz[k] = g1_yz_yyyyz[k] - fr * g2_yz_yyyz[k];

                        g_yyz_yyzz[k] = g1_yz_yyyzz[k] - fr * g2_yz_yyzz[k];

                        g_yyz_yzzz[k] = g1_yz_yyzzz[k] - fr * g2_yz_yzzz[k];

                        g_yyz_zzzz[k] = g1_yz_yzzzz[k] - fr * g2_yz_zzzz[k];

                        g_yzz_xxxx[k] = g1_zz_xxxxy[k] - fr * g2_zz_xxxx[k];

                        g_yzz_xxxy[k] = g1_zz_xxxyy[k] - fr * g2_zz_xxxy[k];

                        g_yzz_xxxz[k] = g1_zz_xxxyz[k] - fr * g2_zz_xxxz[k];

                        g_yzz_xxyy[k] = g1_zz_xxyyy[k] - fr * g2_zz_xxyy[k];

                        g_yzz_xxyz[k] = g1_zz_xxyyz[k] - fr * g2_zz_xxyz[k];

                        g_yzz_xxzz[k] = g1_zz_xxyzz[k] - fr * g2_zz_xxzz[k];

                        g_yzz_xyyy[k] = g1_zz_xyyyy[k] - fr * g2_zz_xyyy[k];

                        g_yzz_xyyz[k] = g1_zz_xyyyz[k] - fr * g2_zz_xyyz[k];

                        g_yzz_xyzz[k] = g1_zz_xyyzz[k] - fr * g2_zz_xyzz[k];

                        g_yzz_xzzz[k] = g1_zz_xyzzz[k] - fr * g2_zz_xzzz[k];

                        g_yzz_yyyy[k] = g1_zz_yyyyy[k] - fr * g2_zz_yyyy[k];

                        g_yzz_yyyz[k] = g1_zz_yyyyz[k] - fr * g2_zz_yyyz[k];

                        g_yzz_yyzz[k] = g1_zz_yyyzz[k] - fr * g2_zz_yyzz[k];

                        g_yzz_yzzz[k] = g1_zz_yyzzz[k] - fr * g2_zz_yzzz[k];

                        g_yzz_zzzz[k] = g1_zz_yzzzz[k] - fr * g2_zz_zzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zzz_xxxx[k] = g1_zz_xxxxz[k] - fr * g2_zz_xxxx[k];

                        g_zzz_xxxy[k] = g1_zz_xxxyz[k] - fr * g2_zz_xxxy[k];

                        g_zzz_xxxz[k] = g1_zz_xxxzz[k] - fr * g2_zz_xxxz[k];

                        g_zzz_xxyy[k] = g1_zz_xxyyz[k] - fr * g2_zz_xxyy[k];

                        g_zzz_xxyz[k] = g1_zz_xxyzz[k] - fr * g2_zz_xxyz[k];

                        g_zzz_xxzz[k] = g1_zz_xxzzz[k] - fr * g2_zz_xxzz[k];

                        g_zzz_xyyy[k] = g1_zz_xyyyz[k] - fr * g2_zz_xyyy[k];

                        g_zzz_xyyz[k] = g1_zz_xyyzz[k] - fr * g2_zz_xyyz[k];

                        g_zzz_xyzz[k] = g1_zz_xyzzz[k] - fr * g2_zz_xyzz[k];

                        g_zzz_xzzz[k] = g1_zz_xzzzz[k] - fr * g2_zz_xzzz[k];

                        g_zzz_yyyy[k] = g1_zz_yyyyz[k] - fr * g2_zz_yyyy[k];

                        g_zzz_yyyz[k] = g1_zz_yyyzz[k] - fr * g2_zz_yyyz[k];

                        g_zzz_yyzz[k] = g1_zz_yyzzz[k] - fr * g2_zz_yyzz[k];

                        g_zzz_yzzz[k] = g1_zz_yzzzz[k] - fr * g2_zz_yzzz[k];

                        g_zzz_zzzz[k] = g1_zz_zzzzz[k] - fr * g2_zz_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXFH(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 3) && (recPattern[i].third() == 5))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|35)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 3, 5});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 6});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 2, 5});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|DH)^(m) integrals

                    auto g2_xx_xxxxx = contrBuffer.data(g2off + 126 * j);

                    auto g2_xx_xxxxy = contrBuffer.data(g2off + 126 * j + 1);

                    auto g2_xx_xxxxz = contrBuffer.data(g2off + 126 * j + 2);

                    auto g2_xx_xxxyy = contrBuffer.data(g2off + 126 * j + 3);

                    auto g2_xx_xxxyz = contrBuffer.data(g2off + 126 * j + 4);

                    auto g2_xx_xxxzz = contrBuffer.data(g2off + 126 * j + 5);

                    auto g2_xx_xxyyy = contrBuffer.data(g2off + 126 * j + 6);

                    auto g2_xx_xxyyz = contrBuffer.data(g2off + 126 * j + 7);

                    auto g2_xx_xxyzz = contrBuffer.data(g2off + 126 * j + 8);

                    auto g2_xx_xxzzz = contrBuffer.data(g2off + 126 * j + 9);

                    auto g2_xx_xyyyy = contrBuffer.data(g2off + 126 * j + 10);

                    auto g2_xx_xyyyz = contrBuffer.data(g2off + 126 * j + 11);

                    auto g2_xx_xyyzz = contrBuffer.data(g2off + 126 * j + 12);

                    auto g2_xx_xyzzz = contrBuffer.data(g2off + 126 * j + 13);

                    auto g2_xx_xzzzz = contrBuffer.data(g2off + 126 * j + 14);

                    auto g2_xx_yyyyy = contrBuffer.data(g2off + 126 * j + 15);

                    auto g2_xx_yyyyz = contrBuffer.data(g2off + 126 * j + 16);

                    auto g2_xx_yyyzz = contrBuffer.data(g2off + 126 * j + 17);

                    auto g2_xx_yyzzz = contrBuffer.data(g2off + 126 * j + 18);

                    auto g2_xx_yzzzz = contrBuffer.data(g2off + 126 * j + 19);

                    auto g2_xx_zzzzz = contrBuffer.data(g2off + 126 * j + 20);

                    auto g2_xy_xxxxx = contrBuffer.data(g2off + 126 * j + 21);

                    auto g2_xy_xxxxy = contrBuffer.data(g2off + 126 * j + 22);

                    auto g2_xy_xxxxz = contrBuffer.data(g2off + 126 * j + 23);

                    auto g2_xy_xxxyy = contrBuffer.data(g2off + 126 * j + 24);

                    auto g2_xy_xxxyz = contrBuffer.data(g2off + 126 * j + 25);

                    auto g2_xy_xxxzz = contrBuffer.data(g2off + 126 * j + 26);

                    auto g2_xy_xxyyy = contrBuffer.data(g2off + 126 * j + 27);

                    auto g2_xy_xxyyz = contrBuffer.data(g2off + 126 * j + 28);

                    auto g2_xy_xxyzz = contrBuffer.data(g2off + 126 * j + 29);

                    auto g2_xy_xxzzz = contrBuffer.data(g2off + 126 * j + 30);

                    auto g2_xy_xyyyy = contrBuffer.data(g2off + 126 * j + 31);

                    auto g2_xy_xyyyz = contrBuffer.data(g2off + 126 * j + 32);

                    auto g2_xy_xyyzz = contrBuffer.data(g2off + 126 * j + 33);

                    auto g2_xy_xyzzz = contrBuffer.data(g2off + 126 * j + 34);

                    auto g2_xy_xzzzz = contrBuffer.data(g2off + 126 * j + 35);

                    auto g2_xy_yyyyy = contrBuffer.data(g2off + 126 * j + 36);

                    auto g2_xy_yyyyz = contrBuffer.data(g2off + 126 * j + 37);

                    auto g2_xy_yyyzz = contrBuffer.data(g2off + 126 * j + 38);

                    auto g2_xy_yyzzz = contrBuffer.data(g2off + 126 * j + 39);

                    auto g2_xy_yzzzz = contrBuffer.data(g2off + 126 * j + 40);

                    auto g2_xy_zzzzz = contrBuffer.data(g2off + 126 * j + 41);

                    auto g2_xz_xxxxx = contrBuffer.data(g2off + 126 * j + 42);

                    auto g2_xz_xxxxy = contrBuffer.data(g2off + 126 * j + 43);

                    auto g2_xz_xxxxz = contrBuffer.data(g2off + 126 * j + 44);

                    auto g2_xz_xxxyy = contrBuffer.data(g2off + 126 * j + 45);

                    auto g2_xz_xxxyz = contrBuffer.data(g2off + 126 * j + 46);

                    auto g2_xz_xxxzz = contrBuffer.data(g2off + 126 * j + 47);

                    auto g2_xz_xxyyy = contrBuffer.data(g2off + 126 * j + 48);

                    auto g2_xz_xxyyz = contrBuffer.data(g2off + 126 * j + 49);

                    auto g2_xz_xxyzz = contrBuffer.data(g2off + 126 * j + 50);

                    auto g2_xz_xxzzz = contrBuffer.data(g2off + 126 * j + 51);

                    auto g2_xz_xyyyy = contrBuffer.data(g2off + 126 * j + 52);

                    auto g2_xz_xyyyz = contrBuffer.data(g2off + 126 * j + 53);

                    auto g2_xz_xyyzz = contrBuffer.data(g2off + 126 * j + 54);

                    auto g2_xz_xyzzz = contrBuffer.data(g2off + 126 * j + 55);

                    auto g2_xz_xzzzz = contrBuffer.data(g2off + 126 * j + 56);

                    auto g2_xz_yyyyy = contrBuffer.data(g2off + 126 * j + 57);

                    auto g2_xz_yyyyz = contrBuffer.data(g2off + 126 * j + 58);

                    auto g2_xz_yyyzz = contrBuffer.data(g2off + 126 * j + 59);

                    auto g2_xz_yyzzz = contrBuffer.data(g2off + 126 * j + 60);

                    auto g2_xz_yzzzz = contrBuffer.data(g2off + 126 * j + 61);

                    auto g2_xz_zzzzz = contrBuffer.data(g2off + 126 * j + 62);

                    auto g2_yy_xxxxx = contrBuffer.data(g2off + 126 * j + 63);

                    auto g2_yy_xxxxy = contrBuffer.data(g2off + 126 * j + 64);

                    auto g2_yy_xxxxz = contrBuffer.data(g2off + 126 * j + 65);

                    auto g2_yy_xxxyy = contrBuffer.data(g2off + 126 * j + 66);

                    auto g2_yy_xxxyz = contrBuffer.data(g2off + 126 * j + 67);

                    auto g2_yy_xxxzz = contrBuffer.data(g2off + 126 * j + 68);

                    auto g2_yy_xxyyy = contrBuffer.data(g2off + 126 * j + 69);

                    auto g2_yy_xxyyz = contrBuffer.data(g2off + 126 * j + 70);

                    auto g2_yy_xxyzz = contrBuffer.data(g2off + 126 * j + 71);

                    auto g2_yy_xxzzz = contrBuffer.data(g2off + 126 * j + 72);

                    auto g2_yy_xyyyy = contrBuffer.data(g2off + 126 * j + 73);

                    auto g2_yy_xyyyz = contrBuffer.data(g2off + 126 * j + 74);

                    auto g2_yy_xyyzz = contrBuffer.data(g2off + 126 * j + 75);

                    auto g2_yy_xyzzz = contrBuffer.data(g2off + 126 * j + 76);

                    auto g2_yy_xzzzz = contrBuffer.data(g2off + 126 * j + 77);

                    auto g2_yy_yyyyy = contrBuffer.data(g2off + 126 * j + 78);

                    auto g2_yy_yyyyz = contrBuffer.data(g2off + 126 * j + 79);

                    auto g2_yy_yyyzz = contrBuffer.data(g2off + 126 * j + 80);

                    auto g2_yy_yyzzz = contrBuffer.data(g2off + 126 * j + 81);

                    auto g2_yy_yzzzz = contrBuffer.data(g2off + 126 * j + 82);

                    auto g2_yy_zzzzz = contrBuffer.data(g2off + 126 * j + 83);

                    auto g2_yz_xxxxx = contrBuffer.data(g2off + 126 * j + 84);

                    auto g2_yz_xxxxy = contrBuffer.data(g2off + 126 * j + 85);

                    auto g2_yz_xxxxz = contrBuffer.data(g2off + 126 * j + 86);

                    auto g2_yz_xxxyy = contrBuffer.data(g2off + 126 * j + 87);

                    auto g2_yz_xxxyz = contrBuffer.data(g2off + 126 * j + 88);

                    auto g2_yz_xxxzz = contrBuffer.data(g2off + 126 * j + 89);

                    auto g2_yz_xxyyy = contrBuffer.data(g2off + 126 * j + 90);

                    auto g2_yz_xxyyz = contrBuffer.data(g2off + 126 * j + 91);

                    auto g2_yz_xxyzz = contrBuffer.data(g2off + 126 * j + 92);

                    auto g2_yz_xxzzz = contrBuffer.data(g2off + 126 * j + 93);

                    auto g2_yz_xyyyy = contrBuffer.data(g2off + 126 * j + 94);

                    auto g2_yz_xyyyz = contrBuffer.data(g2off + 126 * j + 95);

                    auto g2_yz_xyyzz = contrBuffer.data(g2off + 126 * j + 96);

                    auto g2_yz_xyzzz = contrBuffer.data(g2off + 126 * j + 97);

                    auto g2_yz_xzzzz = contrBuffer.data(g2off + 126 * j + 98);

                    auto g2_yz_yyyyy = contrBuffer.data(g2off + 126 * j + 99);

                    auto g2_yz_yyyyz = contrBuffer.data(g2off + 126 * j + 100);

                    auto g2_yz_yyyzz = contrBuffer.data(g2off + 126 * j + 101);

                    auto g2_yz_yyzzz = contrBuffer.data(g2off + 126 * j + 102);

                    auto g2_yz_yzzzz = contrBuffer.data(g2off + 126 * j + 103);

                    auto g2_yz_zzzzz = contrBuffer.data(g2off + 126 * j + 104);

                    auto g2_zz_xxxxx = contrBuffer.data(g2off + 126 * j + 105);

                    auto g2_zz_xxxxy = contrBuffer.data(g2off + 126 * j + 106);

                    auto g2_zz_xxxxz = contrBuffer.data(g2off + 126 * j + 107);

                    auto g2_zz_xxxyy = contrBuffer.data(g2off + 126 * j + 108);

                    auto g2_zz_xxxyz = contrBuffer.data(g2off + 126 * j + 109);

                    auto g2_zz_xxxzz = contrBuffer.data(g2off + 126 * j + 110);

                    auto g2_zz_xxyyy = contrBuffer.data(g2off + 126 * j + 111);

                    auto g2_zz_xxyyz = contrBuffer.data(g2off + 126 * j + 112);

                    auto g2_zz_xxyzz = contrBuffer.data(g2off + 126 * j + 113);

                    auto g2_zz_xxzzz = contrBuffer.data(g2off + 126 * j + 114);

                    auto g2_zz_xyyyy = contrBuffer.data(g2off + 126 * j + 115);

                    auto g2_zz_xyyyz = contrBuffer.data(g2off + 126 * j + 116);

                    auto g2_zz_xyyzz = contrBuffer.data(g2off + 126 * j + 117);

                    auto g2_zz_xyzzz = contrBuffer.data(g2off + 126 * j + 118);

                    auto g2_zz_xzzzz = contrBuffer.data(g2off + 126 * j + 119);

                    auto g2_zz_yyyyy = contrBuffer.data(g2off + 126 * j + 120);

                    auto g2_zz_yyyyz = contrBuffer.data(g2off + 126 * j + 121);

                    auto g2_zz_yyyzz = contrBuffer.data(g2off + 126 * j + 122);

                    auto g2_zz_yyzzz = contrBuffer.data(g2off + 126 * j + 123);

                    auto g2_zz_yzzzz = contrBuffer.data(g2off + 126 * j + 124);

                    auto g2_zz_zzzzz = contrBuffer.data(g2off + 126 * j + 125);

                    // set up pointers to (SX|g(r,r')|DI)^(m) integrals

                    auto g1_xx_xxxxxx = contrBuffer.data(g1off + 168 * j);

                    auto g1_xx_xxxxxy = contrBuffer.data(g1off + 168 * j + 1);

                    auto g1_xx_xxxxxz = contrBuffer.data(g1off + 168 * j + 2);

                    auto g1_xx_xxxxyy = contrBuffer.data(g1off + 168 * j + 3);

                    auto g1_xx_xxxxyz = contrBuffer.data(g1off + 168 * j + 4);

                    auto g1_xx_xxxxzz = contrBuffer.data(g1off + 168 * j + 5);

                    auto g1_xx_xxxyyy = contrBuffer.data(g1off + 168 * j + 6);

                    auto g1_xx_xxxyyz = contrBuffer.data(g1off + 168 * j + 7);

                    auto g1_xx_xxxyzz = contrBuffer.data(g1off + 168 * j + 8);

                    auto g1_xx_xxxzzz = contrBuffer.data(g1off + 168 * j + 9);

                    auto g1_xx_xxyyyy = contrBuffer.data(g1off + 168 * j + 10);

                    auto g1_xx_xxyyyz = contrBuffer.data(g1off + 168 * j + 11);

                    auto g1_xx_xxyyzz = contrBuffer.data(g1off + 168 * j + 12);

                    auto g1_xx_xxyzzz = contrBuffer.data(g1off + 168 * j + 13);

                    auto g1_xx_xxzzzz = contrBuffer.data(g1off + 168 * j + 14);

                    auto g1_xx_xyyyyy = contrBuffer.data(g1off + 168 * j + 15);

                    auto g1_xx_xyyyyz = contrBuffer.data(g1off + 168 * j + 16);

                    auto g1_xx_xyyyzz = contrBuffer.data(g1off + 168 * j + 17);

                    auto g1_xx_xyyzzz = contrBuffer.data(g1off + 168 * j + 18);

                    auto g1_xx_xyzzzz = contrBuffer.data(g1off + 168 * j + 19);

                    auto g1_xx_xzzzzz = contrBuffer.data(g1off + 168 * j + 20);

                    auto g1_xy_xxxxxx = contrBuffer.data(g1off + 168 * j + 28);

                    auto g1_xy_xxxxxy = contrBuffer.data(g1off + 168 * j + 29);

                    auto g1_xy_xxxxxz = contrBuffer.data(g1off + 168 * j + 30);

                    auto g1_xy_xxxxyy = contrBuffer.data(g1off + 168 * j + 31);

                    auto g1_xy_xxxxyz = contrBuffer.data(g1off + 168 * j + 32);

                    auto g1_xy_xxxxzz = contrBuffer.data(g1off + 168 * j + 33);

                    auto g1_xy_xxxyyy = contrBuffer.data(g1off + 168 * j + 34);

                    auto g1_xy_xxxyyz = contrBuffer.data(g1off + 168 * j + 35);

                    auto g1_xy_xxxyzz = contrBuffer.data(g1off + 168 * j + 36);

                    auto g1_xy_xxxzzz = contrBuffer.data(g1off + 168 * j + 37);

                    auto g1_xy_xxyyyy = contrBuffer.data(g1off + 168 * j + 38);

                    auto g1_xy_xxyyyz = contrBuffer.data(g1off + 168 * j + 39);

                    auto g1_xy_xxyyzz = contrBuffer.data(g1off + 168 * j + 40);

                    auto g1_xy_xxyzzz = contrBuffer.data(g1off + 168 * j + 41);

                    auto g1_xy_xxzzzz = contrBuffer.data(g1off + 168 * j + 42);

                    auto g1_xy_xyyyyy = contrBuffer.data(g1off + 168 * j + 43);

                    auto g1_xy_xyyyyz = contrBuffer.data(g1off + 168 * j + 44);

                    auto g1_xy_xyyyzz = contrBuffer.data(g1off + 168 * j + 45);

                    auto g1_xy_xyyzzz = contrBuffer.data(g1off + 168 * j + 46);

                    auto g1_xy_xyzzzz = contrBuffer.data(g1off + 168 * j + 47);

                    auto g1_xy_xzzzzz = contrBuffer.data(g1off + 168 * j + 48);

                    auto g1_xz_xxxxxx = contrBuffer.data(g1off + 168 * j + 56);

                    auto g1_xz_xxxxxy = contrBuffer.data(g1off + 168 * j + 57);

                    auto g1_xz_xxxxxz = contrBuffer.data(g1off + 168 * j + 58);

                    auto g1_xz_xxxxyy = contrBuffer.data(g1off + 168 * j + 59);

                    auto g1_xz_xxxxyz = contrBuffer.data(g1off + 168 * j + 60);

                    auto g1_xz_xxxxzz = contrBuffer.data(g1off + 168 * j + 61);

                    auto g1_xz_xxxyyy = contrBuffer.data(g1off + 168 * j + 62);

                    auto g1_xz_xxxyyz = contrBuffer.data(g1off + 168 * j + 63);

                    auto g1_xz_xxxyzz = contrBuffer.data(g1off + 168 * j + 64);

                    auto g1_xz_xxxzzz = contrBuffer.data(g1off + 168 * j + 65);

                    auto g1_xz_xxyyyy = contrBuffer.data(g1off + 168 * j + 66);

                    auto g1_xz_xxyyyz = contrBuffer.data(g1off + 168 * j + 67);

                    auto g1_xz_xxyyzz = contrBuffer.data(g1off + 168 * j + 68);

                    auto g1_xz_xxyzzz = contrBuffer.data(g1off + 168 * j + 69);

                    auto g1_xz_xxzzzz = contrBuffer.data(g1off + 168 * j + 70);

                    auto g1_xz_xyyyyy = contrBuffer.data(g1off + 168 * j + 71);

                    auto g1_xz_xyyyyz = contrBuffer.data(g1off + 168 * j + 72);

                    auto g1_xz_xyyyzz = contrBuffer.data(g1off + 168 * j + 73);

                    auto g1_xz_xyyzzz = contrBuffer.data(g1off + 168 * j + 74);

                    auto g1_xz_xyzzzz = contrBuffer.data(g1off + 168 * j + 75);

                    auto g1_xz_xzzzzz = contrBuffer.data(g1off + 168 * j + 76);

                    auto g1_yy_xxxxxx = contrBuffer.data(g1off + 168 * j + 84);

                    auto g1_yy_xxxxxy = contrBuffer.data(g1off + 168 * j + 85);

                    auto g1_yy_xxxxxz = contrBuffer.data(g1off + 168 * j + 86);

                    auto g1_yy_xxxxyy = contrBuffer.data(g1off + 168 * j + 87);

                    auto g1_yy_xxxxyz = contrBuffer.data(g1off + 168 * j + 88);

                    auto g1_yy_xxxxzz = contrBuffer.data(g1off + 168 * j + 89);

                    auto g1_yy_xxxyyy = contrBuffer.data(g1off + 168 * j + 90);

                    auto g1_yy_xxxyyz = contrBuffer.data(g1off + 168 * j + 91);

                    auto g1_yy_xxxyzz = contrBuffer.data(g1off + 168 * j + 92);

                    auto g1_yy_xxxzzz = contrBuffer.data(g1off + 168 * j + 93);

                    auto g1_yy_xxyyyy = contrBuffer.data(g1off + 168 * j + 94);

                    auto g1_yy_xxyyyz = contrBuffer.data(g1off + 168 * j + 95);

                    auto g1_yy_xxyyzz = contrBuffer.data(g1off + 168 * j + 96);

                    auto g1_yy_xxyzzz = contrBuffer.data(g1off + 168 * j + 97);

                    auto g1_yy_xxzzzz = contrBuffer.data(g1off + 168 * j + 98);

                    auto g1_yy_xyyyyy = contrBuffer.data(g1off + 168 * j + 99);

                    auto g1_yy_xyyyyz = contrBuffer.data(g1off + 168 * j + 100);

                    auto g1_yy_xyyyzz = contrBuffer.data(g1off + 168 * j + 101);

                    auto g1_yy_xyyzzz = contrBuffer.data(g1off + 168 * j + 102);

                    auto g1_yy_xyzzzz = contrBuffer.data(g1off + 168 * j + 103);

                    auto g1_yy_xzzzzz = contrBuffer.data(g1off + 168 * j + 104);

                    auto g1_yy_yyyyyy = contrBuffer.data(g1off + 168 * j + 105);

                    auto g1_yy_yyyyyz = contrBuffer.data(g1off + 168 * j + 106);

                    auto g1_yy_yyyyzz = contrBuffer.data(g1off + 168 * j + 107);

                    auto g1_yy_yyyzzz = contrBuffer.data(g1off + 168 * j + 108);

                    auto g1_yy_yyzzzz = contrBuffer.data(g1off + 168 * j + 109);

                    auto g1_yy_yzzzzz = contrBuffer.data(g1off + 168 * j + 110);

                    auto g1_yz_xxxxxx = contrBuffer.data(g1off + 168 * j + 112);

                    auto g1_yz_xxxxxy = contrBuffer.data(g1off + 168 * j + 113);

                    auto g1_yz_xxxxxz = contrBuffer.data(g1off + 168 * j + 114);

                    auto g1_yz_xxxxyy = contrBuffer.data(g1off + 168 * j + 115);

                    auto g1_yz_xxxxyz = contrBuffer.data(g1off + 168 * j + 116);

                    auto g1_yz_xxxxzz = contrBuffer.data(g1off + 168 * j + 117);

                    auto g1_yz_xxxyyy = contrBuffer.data(g1off + 168 * j + 118);

                    auto g1_yz_xxxyyz = contrBuffer.data(g1off + 168 * j + 119);

                    auto g1_yz_xxxyzz = contrBuffer.data(g1off + 168 * j + 120);

                    auto g1_yz_xxxzzz = contrBuffer.data(g1off + 168 * j + 121);

                    auto g1_yz_xxyyyy = contrBuffer.data(g1off + 168 * j + 122);

                    auto g1_yz_xxyyyz = contrBuffer.data(g1off + 168 * j + 123);

                    auto g1_yz_xxyyzz = contrBuffer.data(g1off + 168 * j + 124);

                    auto g1_yz_xxyzzz = contrBuffer.data(g1off + 168 * j + 125);

                    auto g1_yz_xxzzzz = contrBuffer.data(g1off + 168 * j + 126);

                    auto g1_yz_xyyyyy = contrBuffer.data(g1off + 168 * j + 127);

                    auto g1_yz_xyyyyz = contrBuffer.data(g1off + 168 * j + 128);

                    auto g1_yz_xyyyzz = contrBuffer.data(g1off + 168 * j + 129);

                    auto g1_yz_xyyzzz = contrBuffer.data(g1off + 168 * j + 130);

                    auto g1_yz_xyzzzz = contrBuffer.data(g1off + 168 * j + 131);

                    auto g1_yz_xzzzzz = contrBuffer.data(g1off + 168 * j + 132);

                    auto g1_yz_yyyyyy = contrBuffer.data(g1off + 168 * j + 133);

                    auto g1_yz_yyyyyz = contrBuffer.data(g1off + 168 * j + 134);

                    auto g1_yz_yyyyzz = contrBuffer.data(g1off + 168 * j + 135);

                    auto g1_yz_yyyzzz = contrBuffer.data(g1off + 168 * j + 136);

                    auto g1_yz_yyzzzz = contrBuffer.data(g1off + 168 * j + 137);

                    auto g1_yz_yzzzzz = contrBuffer.data(g1off + 168 * j + 138);

                    auto g1_zz_xxxxxx = contrBuffer.data(g1off + 168 * j + 140);

                    auto g1_zz_xxxxxy = contrBuffer.data(g1off + 168 * j + 141);

                    auto g1_zz_xxxxxz = contrBuffer.data(g1off + 168 * j + 142);

                    auto g1_zz_xxxxyy = contrBuffer.data(g1off + 168 * j + 143);

                    auto g1_zz_xxxxyz = contrBuffer.data(g1off + 168 * j + 144);

                    auto g1_zz_xxxxzz = contrBuffer.data(g1off + 168 * j + 145);

                    auto g1_zz_xxxyyy = contrBuffer.data(g1off + 168 * j + 146);

                    auto g1_zz_xxxyyz = contrBuffer.data(g1off + 168 * j + 147);

                    auto g1_zz_xxxyzz = contrBuffer.data(g1off + 168 * j + 148);

                    auto g1_zz_xxxzzz = contrBuffer.data(g1off + 168 * j + 149);

                    auto g1_zz_xxyyyy = contrBuffer.data(g1off + 168 * j + 150);

                    auto g1_zz_xxyyyz = contrBuffer.data(g1off + 168 * j + 151);

                    auto g1_zz_xxyyzz = contrBuffer.data(g1off + 168 * j + 152);

                    auto g1_zz_xxyzzz = contrBuffer.data(g1off + 168 * j + 153);

                    auto g1_zz_xxzzzz = contrBuffer.data(g1off + 168 * j + 154);

                    auto g1_zz_xyyyyy = contrBuffer.data(g1off + 168 * j + 155);

                    auto g1_zz_xyyyyz = contrBuffer.data(g1off + 168 * j + 156);

                    auto g1_zz_xyyyzz = contrBuffer.data(g1off + 168 * j + 157);

                    auto g1_zz_xyyzzz = contrBuffer.data(g1off + 168 * j + 158);

                    auto g1_zz_xyzzzz = contrBuffer.data(g1off + 168 * j + 159);

                    auto g1_zz_xzzzzz = contrBuffer.data(g1off + 168 * j + 160);

                    auto g1_zz_yyyyyy = contrBuffer.data(g1off + 168 * j + 161);

                    auto g1_zz_yyyyyz = contrBuffer.data(g1off + 168 * j + 162);

                    auto g1_zz_yyyyzz = contrBuffer.data(g1off + 168 * j + 163);

                    auto g1_zz_yyyzzz = contrBuffer.data(g1off + 168 * j + 164);

                    auto g1_zz_yyzzzz = contrBuffer.data(g1off + 168 * j + 165);

                    auto g1_zz_yzzzzz = contrBuffer.data(g1off + 168 * j + 166);

                    auto g1_zz_zzzzzz = contrBuffer.data(g1off + 168 * j + 167);

                    // set up pointers to (SX|g(r,r')|FH)^(m) integrals

                    auto g_xxx_xxxxx = contrBuffer.data(goff + 210 * j);

                    auto g_xxx_xxxxy = contrBuffer.data(goff + 210 * j + 1);

                    auto g_xxx_xxxxz = contrBuffer.data(goff + 210 * j + 2);

                    auto g_xxx_xxxyy = contrBuffer.data(goff + 210 * j + 3);

                    auto g_xxx_xxxyz = contrBuffer.data(goff + 210 * j + 4);

                    auto g_xxx_xxxzz = contrBuffer.data(goff + 210 * j + 5);

                    auto g_xxx_xxyyy = contrBuffer.data(goff + 210 * j + 6);

                    auto g_xxx_xxyyz = contrBuffer.data(goff + 210 * j + 7);

                    auto g_xxx_xxyzz = contrBuffer.data(goff + 210 * j + 8);

                    auto g_xxx_xxzzz = contrBuffer.data(goff + 210 * j + 9);

                    auto g_xxx_xyyyy = contrBuffer.data(goff + 210 * j + 10);

                    auto g_xxx_xyyyz = contrBuffer.data(goff + 210 * j + 11);

                    auto g_xxx_xyyzz = contrBuffer.data(goff + 210 * j + 12);

                    auto g_xxx_xyzzz = contrBuffer.data(goff + 210 * j + 13);

                    auto g_xxx_xzzzz = contrBuffer.data(goff + 210 * j + 14);

                    auto g_xxx_yyyyy = contrBuffer.data(goff + 210 * j + 15);

                    auto g_xxx_yyyyz = contrBuffer.data(goff + 210 * j + 16);

                    auto g_xxx_yyyzz = contrBuffer.data(goff + 210 * j + 17);

                    auto g_xxx_yyzzz = contrBuffer.data(goff + 210 * j + 18);

                    auto g_xxx_yzzzz = contrBuffer.data(goff + 210 * j + 19);

                    auto g_xxx_zzzzz = contrBuffer.data(goff + 210 * j + 20);

                    auto g_xxy_xxxxx = contrBuffer.data(goff + 210 * j + 21);

                    auto g_xxy_xxxxy = contrBuffer.data(goff + 210 * j + 22);

                    auto g_xxy_xxxxz = contrBuffer.data(goff + 210 * j + 23);

                    auto g_xxy_xxxyy = contrBuffer.data(goff + 210 * j + 24);

                    auto g_xxy_xxxyz = contrBuffer.data(goff + 210 * j + 25);

                    auto g_xxy_xxxzz = contrBuffer.data(goff + 210 * j + 26);

                    auto g_xxy_xxyyy = contrBuffer.data(goff + 210 * j + 27);

                    auto g_xxy_xxyyz = contrBuffer.data(goff + 210 * j + 28);

                    auto g_xxy_xxyzz = contrBuffer.data(goff + 210 * j + 29);

                    auto g_xxy_xxzzz = contrBuffer.data(goff + 210 * j + 30);

                    auto g_xxy_xyyyy = contrBuffer.data(goff + 210 * j + 31);

                    auto g_xxy_xyyyz = contrBuffer.data(goff + 210 * j + 32);

                    auto g_xxy_xyyzz = contrBuffer.data(goff + 210 * j + 33);

                    auto g_xxy_xyzzz = contrBuffer.data(goff + 210 * j + 34);

                    auto g_xxy_xzzzz = contrBuffer.data(goff + 210 * j + 35);

                    auto g_xxy_yyyyy = contrBuffer.data(goff + 210 * j + 36);

                    auto g_xxy_yyyyz = contrBuffer.data(goff + 210 * j + 37);

                    auto g_xxy_yyyzz = contrBuffer.data(goff + 210 * j + 38);

                    auto g_xxy_yyzzz = contrBuffer.data(goff + 210 * j + 39);

                    auto g_xxy_yzzzz = contrBuffer.data(goff + 210 * j + 40);

                    auto g_xxy_zzzzz = contrBuffer.data(goff + 210 * j + 41);

                    auto g_xxz_xxxxx = contrBuffer.data(goff + 210 * j + 42);

                    auto g_xxz_xxxxy = contrBuffer.data(goff + 210 * j + 43);

                    auto g_xxz_xxxxz = contrBuffer.data(goff + 210 * j + 44);

                    auto g_xxz_xxxyy = contrBuffer.data(goff + 210 * j + 45);

                    auto g_xxz_xxxyz = contrBuffer.data(goff + 210 * j + 46);

                    auto g_xxz_xxxzz = contrBuffer.data(goff + 210 * j + 47);

                    auto g_xxz_xxyyy = contrBuffer.data(goff + 210 * j + 48);

                    auto g_xxz_xxyyz = contrBuffer.data(goff + 210 * j + 49);

                    auto g_xxz_xxyzz = contrBuffer.data(goff + 210 * j + 50);

                    auto g_xxz_xxzzz = contrBuffer.data(goff + 210 * j + 51);

                    auto g_xxz_xyyyy = contrBuffer.data(goff + 210 * j + 52);

                    auto g_xxz_xyyyz = contrBuffer.data(goff + 210 * j + 53);

                    auto g_xxz_xyyzz = contrBuffer.data(goff + 210 * j + 54);

                    auto g_xxz_xyzzz = contrBuffer.data(goff + 210 * j + 55);

                    auto g_xxz_xzzzz = contrBuffer.data(goff + 210 * j + 56);

                    auto g_xxz_yyyyy = contrBuffer.data(goff + 210 * j + 57);

                    auto g_xxz_yyyyz = contrBuffer.data(goff + 210 * j + 58);

                    auto g_xxz_yyyzz = contrBuffer.data(goff + 210 * j + 59);

                    auto g_xxz_yyzzz = contrBuffer.data(goff + 210 * j + 60);

                    auto g_xxz_yzzzz = contrBuffer.data(goff + 210 * j + 61);

                    auto g_xxz_zzzzz = contrBuffer.data(goff + 210 * j + 62);

                    auto g_xyy_xxxxx = contrBuffer.data(goff + 210 * j + 63);

                    auto g_xyy_xxxxy = contrBuffer.data(goff + 210 * j + 64);

                    auto g_xyy_xxxxz = contrBuffer.data(goff + 210 * j + 65);

                    auto g_xyy_xxxyy = contrBuffer.data(goff + 210 * j + 66);

                    auto g_xyy_xxxyz = contrBuffer.data(goff + 210 * j + 67);

                    auto g_xyy_xxxzz = contrBuffer.data(goff + 210 * j + 68);

                    auto g_xyy_xxyyy = contrBuffer.data(goff + 210 * j + 69);

                    auto g_xyy_xxyyz = contrBuffer.data(goff + 210 * j + 70);

                    auto g_xyy_xxyzz = contrBuffer.data(goff + 210 * j + 71);

                    auto g_xyy_xxzzz = contrBuffer.data(goff + 210 * j + 72);

                    auto g_xyy_xyyyy = contrBuffer.data(goff + 210 * j + 73);

                    auto g_xyy_xyyyz = contrBuffer.data(goff + 210 * j + 74);

                    auto g_xyy_xyyzz = contrBuffer.data(goff + 210 * j + 75);

                    auto g_xyy_xyzzz = contrBuffer.data(goff + 210 * j + 76);

                    auto g_xyy_xzzzz = contrBuffer.data(goff + 210 * j + 77);

                    auto g_xyy_yyyyy = contrBuffer.data(goff + 210 * j + 78);

                    auto g_xyy_yyyyz = contrBuffer.data(goff + 210 * j + 79);

                    auto g_xyy_yyyzz = contrBuffer.data(goff + 210 * j + 80);

                    auto g_xyy_yyzzz = contrBuffer.data(goff + 210 * j + 81);

                    auto g_xyy_yzzzz = contrBuffer.data(goff + 210 * j + 82);

                    auto g_xyy_zzzzz = contrBuffer.data(goff + 210 * j + 83);

                    auto g_xyz_xxxxx = contrBuffer.data(goff + 210 * j + 84);

                    auto g_xyz_xxxxy = contrBuffer.data(goff + 210 * j + 85);

                    auto g_xyz_xxxxz = contrBuffer.data(goff + 210 * j + 86);

                    auto g_xyz_xxxyy = contrBuffer.data(goff + 210 * j + 87);

                    auto g_xyz_xxxyz = contrBuffer.data(goff + 210 * j + 88);

                    auto g_xyz_xxxzz = contrBuffer.data(goff + 210 * j + 89);

                    auto g_xyz_xxyyy = contrBuffer.data(goff + 210 * j + 90);

                    auto g_xyz_xxyyz = contrBuffer.data(goff + 210 * j + 91);

                    auto g_xyz_xxyzz = contrBuffer.data(goff + 210 * j + 92);

                    auto g_xyz_xxzzz = contrBuffer.data(goff + 210 * j + 93);

                    auto g_xyz_xyyyy = contrBuffer.data(goff + 210 * j + 94);

                    auto g_xyz_xyyyz = contrBuffer.data(goff + 210 * j + 95);

                    auto g_xyz_xyyzz = contrBuffer.data(goff + 210 * j + 96);

                    auto g_xyz_xyzzz = contrBuffer.data(goff + 210 * j + 97);

                    auto g_xyz_xzzzz = contrBuffer.data(goff + 210 * j + 98);

                    auto g_xyz_yyyyy = contrBuffer.data(goff + 210 * j + 99);

                    auto g_xyz_yyyyz = contrBuffer.data(goff + 210 * j + 100);

                    auto g_xyz_yyyzz = contrBuffer.data(goff + 210 * j + 101);

                    auto g_xyz_yyzzz = contrBuffer.data(goff + 210 * j + 102);

                    auto g_xyz_yzzzz = contrBuffer.data(goff + 210 * j + 103);

                    auto g_xyz_zzzzz = contrBuffer.data(goff + 210 * j + 104);

                    auto g_xzz_xxxxx = contrBuffer.data(goff + 210 * j + 105);

                    auto g_xzz_xxxxy = contrBuffer.data(goff + 210 * j + 106);

                    auto g_xzz_xxxxz = contrBuffer.data(goff + 210 * j + 107);

                    auto g_xzz_xxxyy = contrBuffer.data(goff + 210 * j + 108);

                    auto g_xzz_xxxyz = contrBuffer.data(goff + 210 * j + 109);

                    auto g_xzz_xxxzz = contrBuffer.data(goff + 210 * j + 110);

                    auto g_xzz_xxyyy = contrBuffer.data(goff + 210 * j + 111);

                    auto g_xzz_xxyyz = contrBuffer.data(goff + 210 * j + 112);

                    auto g_xzz_xxyzz = contrBuffer.data(goff + 210 * j + 113);

                    auto g_xzz_xxzzz = contrBuffer.data(goff + 210 * j + 114);

                    auto g_xzz_xyyyy = contrBuffer.data(goff + 210 * j + 115);

                    auto g_xzz_xyyyz = contrBuffer.data(goff + 210 * j + 116);

                    auto g_xzz_xyyzz = contrBuffer.data(goff + 210 * j + 117);

                    auto g_xzz_xyzzz = contrBuffer.data(goff + 210 * j + 118);

                    auto g_xzz_xzzzz = contrBuffer.data(goff + 210 * j + 119);

                    auto g_xzz_yyyyy = contrBuffer.data(goff + 210 * j + 120);

                    auto g_xzz_yyyyz = contrBuffer.data(goff + 210 * j + 121);

                    auto g_xzz_yyyzz = contrBuffer.data(goff + 210 * j + 122);

                    auto g_xzz_yyzzz = contrBuffer.data(goff + 210 * j + 123);

                    auto g_xzz_yzzzz = contrBuffer.data(goff + 210 * j + 124);

                    auto g_xzz_zzzzz = contrBuffer.data(goff + 210 * j + 125);

                    auto g_yyy_xxxxx = contrBuffer.data(goff + 210 * j + 126);

                    auto g_yyy_xxxxy = contrBuffer.data(goff + 210 * j + 127);

                    auto g_yyy_xxxxz = contrBuffer.data(goff + 210 * j + 128);

                    auto g_yyy_xxxyy = contrBuffer.data(goff + 210 * j + 129);

                    auto g_yyy_xxxyz = contrBuffer.data(goff + 210 * j + 130);

                    auto g_yyy_xxxzz = contrBuffer.data(goff + 210 * j + 131);

                    auto g_yyy_xxyyy = contrBuffer.data(goff + 210 * j + 132);

                    auto g_yyy_xxyyz = contrBuffer.data(goff + 210 * j + 133);

                    auto g_yyy_xxyzz = contrBuffer.data(goff + 210 * j + 134);

                    auto g_yyy_xxzzz = contrBuffer.data(goff + 210 * j + 135);

                    auto g_yyy_xyyyy = contrBuffer.data(goff + 210 * j + 136);

                    auto g_yyy_xyyyz = contrBuffer.data(goff + 210 * j + 137);

                    auto g_yyy_xyyzz = contrBuffer.data(goff + 210 * j + 138);

                    auto g_yyy_xyzzz = contrBuffer.data(goff + 210 * j + 139);

                    auto g_yyy_xzzzz = contrBuffer.data(goff + 210 * j + 140);

                    auto g_yyy_yyyyy = contrBuffer.data(goff + 210 * j + 141);

                    auto g_yyy_yyyyz = contrBuffer.data(goff + 210 * j + 142);

                    auto g_yyy_yyyzz = contrBuffer.data(goff + 210 * j + 143);

                    auto g_yyy_yyzzz = contrBuffer.data(goff + 210 * j + 144);

                    auto g_yyy_yzzzz = contrBuffer.data(goff + 210 * j + 145);

                    auto g_yyy_zzzzz = contrBuffer.data(goff + 210 * j + 146);

                    auto g_yyz_xxxxx = contrBuffer.data(goff + 210 * j + 147);

                    auto g_yyz_xxxxy = contrBuffer.data(goff + 210 * j + 148);

                    auto g_yyz_xxxxz = contrBuffer.data(goff + 210 * j + 149);

                    auto g_yyz_xxxyy = contrBuffer.data(goff + 210 * j + 150);

                    auto g_yyz_xxxyz = contrBuffer.data(goff + 210 * j + 151);

                    auto g_yyz_xxxzz = contrBuffer.data(goff + 210 * j + 152);

                    auto g_yyz_xxyyy = contrBuffer.data(goff + 210 * j + 153);

                    auto g_yyz_xxyyz = contrBuffer.data(goff + 210 * j + 154);

                    auto g_yyz_xxyzz = contrBuffer.data(goff + 210 * j + 155);

                    auto g_yyz_xxzzz = contrBuffer.data(goff + 210 * j + 156);

                    auto g_yyz_xyyyy = contrBuffer.data(goff + 210 * j + 157);

                    auto g_yyz_xyyyz = contrBuffer.data(goff + 210 * j + 158);

                    auto g_yyz_xyyzz = contrBuffer.data(goff + 210 * j + 159);

                    auto g_yyz_xyzzz = contrBuffer.data(goff + 210 * j + 160);

                    auto g_yyz_xzzzz = contrBuffer.data(goff + 210 * j + 161);

                    auto g_yyz_yyyyy = contrBuffer.data(goff + 210 * j + 162);

                    auto g_yyz_yyyyz = contrBuffer.data(goff + 210 * j + 163);

                    auto g_yyz_yyyzz = contrBuffer.data(goff + 210 * j + 164);

                    auto g_yyz_yyzzz = contrBuffer.data(goff + 210 * j + 165);

                    auto g_yyz_yzzzz = contrBuffer.data(goff + 210 * j + 166);

                    auto g_yyz_zzzzz = contrBuffer.data(goff + 210 * j + 167);

                    auto g_yzz_xxxxx = contrBuffer.data(goff + 210 * j + 168);

                    auto g_yzz_xxxxy = contrBuffer.data(goff + 210 * j + 169);

                    auto g_yzz_xxxxz = contrBuffer.data(goff + 210 * j + 170);

                    auto g_yzz_xxxyy = contrBuffer.data(goff + 210 * j + 171);

                    auto g_yzz_xxxyz = contrBuffer.data(goff + 210 * j + 172);

                    auto g_yzz_xxxzz = contrBuffer.data(goff + 210 * j + 173);

                    auto g_yzz_xxyyy = contrBuffer.data(goff + 210 * j + 174);

                    auto g_yzz_xxyyz = contrBuffer.data(goff + 210 * j + 175);

                    auto g_yzz_xxyzz = contrBuffer.data(goff + 210 * j + 176);

                    auto g_yzz_xxzzz = contrBuffer.data(goff + 210 * j + 177);

                    auto g_yzz_xyyyy = contrBuffer.data(goff + 210 * j + 178);

                    auto g_yzz_xyyyz = contrBuffer.data(goff + 210 * j + 179);

                    auto g_yzz_xyyzz = contrBuffer.data(goff + 210 * j + 180);

                    auto g_yzz_xyzzz = contrBuffer.data(goff + 210 * j + 181);

                    auto g_yzz_xzzzz = contrBuffer.data(goff + 210 * j + 182);

                    auto g_yzz_yyyyy = contrBuffer.data(goff + 210 * j + 183);

                    auto g_yzz_yyyyz = contrBuffer.data(goff + 210 * j + 184);

                    auto g_yzz_yyyzz = contrBuffer.data(goff + 210 * j + 185);

                    auto g_yzz_yyzzz = contrBuffer.data(goff + 210 * j + 186);

                    auto g_yzz_yzzzz = contrBuffer.data(goff + 210 * j + 187);

                    auto g_yzz_zzzzz = contrBuffer.data(goff + 210 * j + 188);

                    auto g_zzz_xxxxx = contrBuffer.data(goff + 210 * j + 189);

                    auto g_zzz_xxxxy = contrBuffer.data(goff + 210 * j + 190);

                    auto g_zzz_xxxxz = contrBuffer.data(goff + 210 * j + 191);

                    auto g_zzz_xxxyy = contrBuffer.data(goff + 210 * j + 192);

                    auto g_zzz_xxxyz = contrBuffer.data(goff + 210 * j + 193);

                    auto g_zzz_xxxzz = contrBuffer.data(goff + 210 * j + 194);

                    auto g_zzz_xxyyy = contrBuffer.data(goff + 210 * j + 195);

                    auto g_zzz_xxyyz = contrBuffer.data(goff + 210 * j + 196);

                    auto g_zzz_xxyzz = contrBuffer.data(goff + 210 * j + 197);

                    auto g_zzz_xxzzz = contrBuffer.data(goff + 210 * j + 198);

                    auto g_zzz_xyyyy = contrBuffer.data(goff + 210 * j + 199);

                    auto g_zzz_xyyyz = contrBuffer.data(goff + 210 * j + 200);

                    auto g_zzz_xyyzz = contrBuffer.data(goff + 210 * j + 201);

                    auto g_zzz_xyzzz = contrBuffer.data(goff + 210 * j + 202);

                    auto g_zzz_xzzzz = contrBuffer.data(goff + 210 * j + 203);

                    auto g_zzz_yyyyy = contrBuffer.data(goff + 210 * j + 204);

                    auto g_zzz_yyyyz = contrBuffer.data(goff + 210 * j + 205);

                    auto g_zzz_yyyzz = contrBuffer.data(goff + 210 * j + 206);

                    auto g_zzz_yyzzz = contrBuffer.data(goff + 210 * j + 207);

                    auto g_zzz_yzzzz = contrBuffer.data(goff + 210 * j + 208);

                    auto g_zzz_zzzzz = contrBuffer.data(goff + 210 * j + 209);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xx_xxxxx, g2_xx_xxxxy,\
                                             g2_xx_xxxxz, g2_xx_xxxyy, g2_xx_xxxyz,\
                                             g2_xx_xxxzz, g2_xx_xxyyy, g2_xx_xxyyz,\
                                             g2_xx_xxyzz, g2_xx_xxzzz, g2_xx_xyyyy,\
                                             g2_xx_xyyyz, g2_xx_xyyzz, g2_xx_xyzzz,\
                                             g2_xx_xzzzz, g2_xx_yyyyy, g2_xx_yyyyz,\
                                             g2_xx_yyyzz, g2_xx_yyzzz, g2_xx_yzzzz,\
                                             g2_xx_zzzzz, g2_xy_xxxxx, g2_xy_xxxxy,\
                                             g2_xy_xxxxz, g2_xy_xxxyy, g2_xy_xxxyz,\
                                             g2_xy_xxxzz, g2_xy_xxyyy, g2_xy_xxyyz,\
                                             g2_xy_xxyzz, g2_xy_xxzzz, g2_xy_xyyyy,\
                                             g2_xy_xyyyz, g2_xy_xyyzz, g2_xy_xyzzz,\
                                             g2_xy_xzzzz, g2_xy_yyyyy, g2_xy_yyyyz,\
                                             g2_xy_yyyzz, g2_xy_yyzzz, g2_xy_yzzzz,\
                                             g2_xy_zzzzz, g2_xz_xxxxx, g2_xz_xxxxy,\
                                             g2_xz_xxxxz, g2_xz_xxxyy, g2_xz_xxxyz,\
                                             g2_xz_xxxzz, g2_xz_xxyyy, g2_xz_xxyyz,\
                                             g2_xz_xxyzz, g2_xz_xxzzz, g2_xz_xyyyy,\
                                             g2_xz_xyyyz, g2_xz_xyyzz, g2_xz_xyzzz,\
                                             g2_xz_xzzzz, g2_xz_yyyyy, g2_xz_yyyyz,\
                                             g2_xz_yyyzz, g2_xz_yyzzz, g2_xz_yzzzz,\
                                             g2_xz_zzzzz, g2_yy_xxxxx, g2_yy_xxxxy,\
                                             g2_yy_xxxxz, g2_yy_xxxyy, g2_yy_xxxyz,\
                                             g2_yy_xxxzz, g2_yy_xxyyy, g2_yy_xxyyz,\
                                             g2_yy_xxyzz, g2_yy_xxzzz, g2_yy_xyyyy,\
                                             g2_yy_xyyyz, g2_yy_xyyzz, g2_yy_xyzzz,\
                                             g2_yy_xzzzz, g2_yy_yyyyy, g2_yy_yyyyz,\
                                             g2_yy_yyyzz, g2_yy_yyzzz, g2_yy_yzzzz,\
                                             g2_yy_zzzzz, g2_yz_xxxxx, g2_yz_xxxxy,\
                                             g2_yz_xxxxz, g2_yz_xxxyy, g2_yz_xxxyz,\
                                             g2_yz_xxxzz, g2_yz_xxyyy, g2_yz_xxyyz,\
                                             g2_yz_xxyzz, g2_yz_xxzzz, g2_yz_xyyyy,\
                                             g2_yz_xyyyz, g2_yz_xyyzz, g2_yz_xyzzz,\
                                             g2_yz_xzzzz, g2_yz_yyyyy, g2_yz_yyyyz,\
                                             g2_yz_yyyzz, g2_yz_yyzzz, g2_yz_yzzzz,\
                                             g2_yz_zzzzz, g2_zz_xxxxx, g2_zz_xxxxy,\
                                             g2_zz_xxxxz, g2_zz_xxxyy, g2_zz_xxxyz,\
                                             g2_zz_xxxzz, g2_zz_xxyyy, g2_zz_xxyyz,\
                                             g2_zz_xxyzz, g2_zz_xxzzz, g2_zz_xyyyy,\
                                             g2_zz_xyyyz, g2_zz_xyyzz, g2_zz_xyzzz,\
                                             g2_zz_xzzzz, g2_zz_yyyyy, g2_zz_yyyyz,\
                                             g2_zz_yyyzz, g2_zz_yyzzz, g2_zz_yzzzz,\
                                             g2_zz_zzzzz, g1_xx_xxxxxx, g1_xx_xxxxxy,\
                                             g1_xx_xxxxxz, g1_xx_xxxxyy, g1_xx_xxxxyz,\
                                             g1_xx_xxxxzz, g1_xx_xxxyyy, g1_xx_xxxyyz,\
                                             g1_xx_xxxyzz, g1_xx_xxxzzz, g1_xx_xxyyyy,\
                                             g1_xx_xxyyyz, g1_xx_xxyyzz, g1_xx_xxyzzz,\
                                             g1_xx_xxzzzz, g1_xx_xyyyyy, g1_xx_xyyyyz,\
                                             g1_xx_xyyyzz, g1_xx_xyyzzz, g1_xx_xyzzzz,\
                                             g1_xx_xzzzzz, g1_xy_xxxxxx,\
                                             g1_xy_xxxxxy, g1_xy_xxxxxz, g1_xy_xxxxyy,\
                                             g1_xy_xxxxyz, g1_xy_xxxxzz, g1_xy_xxxyyy,\
                                             g1_xy_xxxyyz, g1_xy_xxxyzz, g1_xy_xxxzzz,\
                                             g1_xy_xxyyyy, g1_xy_xxyyyz, g1_xy_xxyyzz,\
                                             g1_xy_xxyzzz, g1_xy_xxzzzz, g1_xy_xyyyyy,\
                                             g1_xy_xyyyyz, g1_xy_xyyyzz, g1_xy_xyyzzz,\
                                             g1_xy_xyzzzz, g1_xy_xzzzzz,\
                                             g1_xz_xxxxxx, g1_xz_xxxxxy, g1_xz_xxxxxz,\
                                             g1_xz_xxxxyy, g1_xz_xxxxyz, g1_xz_xxxxzz,\
                                             g1_xz_xxxyyy, g1_xz_xxxyyz, g1_xz_xxxyzz,\
                                             g1_xz_xxxzzz, g1_xz_xxyyyy, g1_xz_xxyyyz,\
                                             g1_xz_xxyyzz, g1_xz_xxyzzz, g1_xz_xxzzzz,\
                                             g1_xz_xyyyyy, g1_xz_xyyyyz, g1_xz_xyyyzz,\
                                             g1_xz_xyyzzz, g1_xz_xyzzzz, g1_xz_xzzzzz,\
                                             g1_yy_xxxxxx, g1_yy_xxxxxy,\
                                             g1_yy_xxxxxz, g1_yy_xxxxyy, g1_yy_xxxxyz,\
                                             g1_yy_xxxxzz, g1_yy_xxxyyy, g1_yy_xxxyyz,\
                                             g1_yy_xxxyzz, g1_yy_xxxzzz, g1_yy_xxyyyy,\
                                             g1_yy_xxyyyz, g1_yy_xxyyzz, g1_yy_xxyzzz,\
                                             g1_yy_xxzzzz, g1_yy_xyyyyy, g1_yy_xyyyyz,\
                                             g1_yy_xyyyzz, g1_yy_xyyzzz, g1_yy_xyzzzz,\
                                             g1_yy_xzzzzz, g1_yy_yyyyyy, g1_yy_yyyyyz,\
                                             g1_yy_yyyyzz, g1_yy_yyyzzz, g1_yy_yyzzzz,\
                                             g1_yy_yzzzzz, g1_yz_xxxxxx,\
                                             g1_yz_xxxxxy, g1_yz_xxxxxz, g1_yz_xxxxyy,\
                                             g1_yz_xxxxyz, g1_yz_xxxxzz, g1_yz_xxxyyy,\
                                             g1_yz_xxxyyz, g1_yz_xxxyzz, g1_yz_xxxzzz,\
                                             g1_yz_xxyyyy, g1_yz_xxyyyz, g1_yz_xxyyzz,\
                                             g1_yz_xxyzzz, g1_yz_xxzzzz, g1_yz_xyyyyy,\
                                             g1_yz_xyyyyz, g1_yz_xyyyzz, g1_yz_xyyzzz,\
                                             g1_yz_xyzzzz, g1_yz_xzzzzz, g1_yz_yyyyyy,\
                                             g1_yz_yyyyyz, g1_yz_yyyyzz, g1_yz_yyyzzz,\
                                             g1_yz_yyzzzz, g1_yz_yzzzzz,\
                                             g1_zz_xxxxxx, g1_zz_xxxxxy, g1_zz_xxxxxz,\
                                             g1_zz_xxxxyy, g1_zz_xxxxyz, g1_zz_xxxxzz,\
                                             g1_zz_xxxyyy, g1_zz_xxxyyz, g1_zz_xxxyzz,\
                                             g1_zz_xxxzzz, g1_zz_xxyyyy, g1_zz_xxyyyz,\
                                             g1_zz_xxyyzz, g1_zz_xxyzzz, g1_zz_xxzzzz,\
                                             g1_zz_xyyyyy, g1_zz_xyyyyz, g1_zz_xyyyzz,\
                                             g1_zz_xyyzzz, g1_zz_xyzzzz, g1_zz_xzzzzz,\
                                             g1_zz_yyyyyy, g1_zz_yyyyyz, g1_zz_yyyyzz,\
                                             g1_zz_yyyzzz, g1_zz_yyzzzz, g1_zz_yzzzzz,\
                                             g1_zz_zzzzzz, g_xxx_xxxxx, g_xxx_xxxxy,\
                                             g_xxx_xxxxz, g_xxx_xxxyy, g_xxx_xxxyz,\
                                             g_xxx_xxxzz, g_xxx_xxyyy, g_xxx_xxyyz,\
                                             g_xxx_xxyzz, g_xxx_xxzzz, g_xxx_xyyyy,\
                                             g_xxx_xyyyz, g_xxx_xyyzz, g_xxx_xyzzz,\
                                             g_xxx_xzzzz, g_xxx_yyyyy, g_xxx_yyyyz,\
                                             g_xxx_yyyzz, g_xxx_yyzzz, g_xxx_yzzzz,\
                                             g_xxx_zzzzz, g_xxy_xxxxx, g_xxy_xxxxy,\
                                             g_xxy_xxxxz, g_xxy_xxxyy, g_xxy_xxxyz,\
                                             g_xxy_xxxzz, g_xxy_xxyyy, g_xxy_xxyyz,\
                                             g_xxy_xxyzz, g_xxy_xxzzz, g_xxy_xyyyy,\
                                             g_xxy_xyyyz, g_xxy_xyyzz, g_xxy_xyzzz,\
                                             g_xxy_xzzzz, g_xxy_yyyyy, g_xxy_yyyyz,\
                                             g_xxy_yyyzz, g_xxy_yyzzz, g_xxy_yzzzz,\
                                             g_xxy_zzzzz, g_xxz_xxxxx, g_xxz_xxxxy,\
                                             g_xxz_xxxxz, g_xxz_xxxyy, g_xxz_xxxyz,\
                                             g_xxz_xxxzz, g_xxz_xxyyy, g_xxz_xxyyz,\
                                             g_xxz_xxyzz, g_xxz_xxzzz, g_xxz_xyyyy,\
                                             g_xxz_xyyyz, g_xxz_xyyzz, g_xxz_xyzzz,\
                                             g_xxz_xzzzz, g_xxz_yyyyy, g_xxz_yyyyz,\
                                             g_xxz_yyyzz, g_xxz_yyzzz, g_xxz_yzzzz,\
                                             g_xxz_zzzzz, g_xyy_xxxxx, g_xyy_xxxxy,\
                                             g_xyy_xxxxz, g_xyy_xxxyy, g_xyy_xxxyz,\
                                             g_xyy_xxxzz, g_xyy_xxyyy, g_xyy_xxyyz,\
                                             g_xyy_xxyzz, g_xyy_xxzzz, g_xyy_xyyyy,\
                                             g_xyy_xyyyz, g_xyy_xyyzz, g_xyy_xyzzz,\
                                             g_xyy_xzzzz, g_xyy_yyyyy, g_xyy_yyyyz,\
                                             g_xyy_yyyzz, g_xyy_yyzzz, g_xyy_yzzzz,\
                                             g_xyy_zzzzz, g_xyz_xxxxx, g_xyz_xxxxy,\
                                             g_xyz_xxxxz, g_xyz_xxxyy, g_xyz_xxxyz,\
                                             g_xyz_xxxzz, g_xyz_xxyyy, g_xyz_xxyyz,\
                                             g_xyz_xxyzz, g_xyz_xxzzz, g_xyz_xyyyy,\
                                             g_xyz_xyyyz, g_xyz_xyyzz, g_xyz_xyzzz,\
                                             g_xyz_xzzzz, g_xyz_yyyyy, g_xyz_yyyyz,\
                                             g_xyz_yyyzz, g_xyz_yyzzz, g_xyz_yzzzz,\
                                             g_xyz_zzzzz, g_xzz_xxxxx, g_xzz_xxxxy,\
                                             g_xzz_xxxxz, g_xzz_xxxyy, g_xzz_xxxyz,\
                                             g_xzz_xxxzz, g_xzz_xxyyy, g_xzz_xxyyz,\
                                             g_xzz_xxyzz, g_xzz_xxzzz, g_xzz_xyyyy,\
                                             g_xzz_xyyyz, g_xzz_xyyzz, g_xzz_xyzzz,\
                                             g_xzz_xzzzz, g_xzz_yyyyy, g_xzz_yyyyz,\
                                             g_xzz_yyyzz, g_xzz_yyzzz, g_xzz_yzzzz,\
                                             g_xzz_zzzzz, g_yyy_xxxxx, g_yyy_xxxxy,\
                                             g_yyy_xxxxz, g_yyy_xxxyy, g_yyy_xxxyz,\
                                             g_yyy_xxxzz, g_yyy_xxyyy, g_yyy_xxyyz,\
                                             g_yyy_xxyzz, g_yyy_xxzzz, g_yyy_xyyyy,\
                                             g_yyy_xyyyz, g_yyy_xyyzz, g_yyy_xyzzz,\
                                             g_yyy_xzzzz, g_yyy_yyyyy, g_yyy_yyyyz,\
                                             g_yyy_yyyzz, g_yyy_yyzzz, g_yyy_yzzzz,\
                                             g_yyy_zzzzz, g_yyz_xxxxx, g_yyz_xxxxy,\
                                             g_yyz_xxxxz, g_yyz_xxxyy, g_yyz_xxxyz,\
                                             g_yyz_xxxzz, g_yyz_xxyyy, g_yyz_xxyyz,\
                                             g_yyz_xxyzz, g_yyz_xxzzz, g_yyz_xyyyy,\
                                             g_yyz_xyyyz, g_yyz_xyyzz, g_yyz_xyzzz,\
                                             g_yyz_xzzzz, g_yyz_yyyyy, g_yyz_yyyyz,\
                                             g_yyz_yyyzz, g_yyz_yyzzz, g_yyz_yzzzz,\
                                             g_yyz_zzzzz, g_yzz_xxxxx, g_yzz_xxxxy,\
                                             g_yzz_xxxxz, g_yzz_xxxyy, g_yzz_xxxyz,\
                                             g_yzz_xxxzz, g_yzz_xxyyy, g_yzz_xxyyz,\
                                             g_yzz_xxyzz, g_yzz_xxzzz, g_yzz_xyyyy,\
                                             g_yzz_xyyyz, g_yzz_xyyzz, g_yzz_xyzzz,\
                                             g_yzz_xzzzz, g_yzz_yyyyy, g_yzz_yyyyz,\
                                             g_yzz_yyyzz, g_yzz_yyzzz, g_yzz_yzzzz,\
                                             g_yzz_zzzzz, g_zzz_xxxxx, g_zzz_xxxxy,\
                                             g_zzz_xxxxz, g_zzz_xxxyy, g_zzz_xxxyz,\
                                             g_zzz_xxxzz, g_zzz_xxyyy, g_zzz_xxyyz,\
                                             g_zzz_xxyzz, g_zzz_xxzzz, g_zzz_xyyyy,\
                                             g_zzz_xyyyz, g_zzz_xyyzz, g_zzz_xyzzz,\
                                             g_zzz_xzzzz, g_zzz_yyyyy, g_zzz_yyyyz,\
                                             g_zzz_yyyzz, g_zzz_yyzzz, g_zzz_yzzzz,\
                                             g_zzz_zzzzz: VLX_ALIGN)
                     for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xxx_xxxxx[k] = g1_xx_xxxxxx[k] - fr * g2_xx_xxxxx[k];

                        g_xxx_xxxxy[k] = g1_xx_xxxxxy[k] - fr * g2_xx_xxxxy[k];

                        g_xxx_xxxxz[k] = g1_xx_xxxxxz[k] - fr * g2_xx_xxxxz[k];

                        g_xxx_xxxyy[k] = g1_xx_xxxxyy[k] - fr * g2_xx_xxxyy[k];

                        g_xxx_xxxyz[k] = g1_xx_xxxxyz[k] - fr * g2_xx_xxxyz[k];

                        g_xxx_xxxzz[k] = g1_xx_xxxxzz[k] - fr * g2_xx_xxxzz[k];

                        g_xxx_xxyyy[k] = g1_xx_xxxyyy[k] - fr * g2_xx_xxyyy[k];

                        g_xxx_xxyyz[k] = g1_xx_xxxyyz[k] - fr * g2_xx_xxyyz[k];

                        g_xxx_xxyzz[k] = g1_xx_xxxyzz[k] - fr * g2_xx_xxyzz[k];

                        g_xxx_xxzzz[k] = g1_xx_xxxzzz[k] - fr * g2_xx_xxzzz[k];

                        g_xxx_xyyyy[k] = g1_xx_xxyyyy[k] - fr * g2_xx_xyyyy[k];

                        g_xxx_xyyyz[k] = g1_xx_xxyyyz[k] - fr * g2_xx_xyyyz[k];

                        g_xxx_xyyzz[k] = g1_xx_xxyyzz[k] - fr * g2_xx_xyyzz[k];

                        g_xxx_xyzzz[k] = g1_xx_xxyzzz[k] - fr * g2_xx_xyzzz[k];

                        g_xxx_xzzzz[k] = g1_xx_xxzzzz[k] - fr * g2_xx_xzzzz[k];

                        g_xxx_yyyyy[k] = g1_xx_xyyyyy[k] - fr * g2_xx_yyyyy[k];

                        g_xxx_yyyyz[k] = g1_xx_xyyyyz[k] - fr * g2_xx_yyyyz[k];

                        g_xxx_yyyzz[k] = g1_xx_xyyyzz[k] - fr * g2_xx_yyyzz[k];

                        g_xxx_yyzzz[k] = g1_xx_xyyzzz[k] - fr * g2_xx_yyzzz[k];

                        g_xxx_yzzzz[k] = g1_xx_xyzzzz[k] - fr * g2_xx_yzzzz[k];

                        g_xxx_zzzzz[k] = g1_xx_xzzzzz[k] - fr * g2_xx_zzzzz[k];

                        g_xxy_xxxxx[k] = g1_xy_xxxxxx[k] - fr * g2_xy_xxxxx[k];

                        g_xxy_xxxxy[k] = g1_xy_xxxxxy[k] - fr * g2_xy_xxxxy[k];

                        g_xxy_xxxxz[k] = g1_xy_xxxxxz[k] - fr * g2_xy_xxxxz[k];

                        g_xxy_xxxyy[k] = g1_xy_xxxxyy[k] - fr * g2_xy_xxxyy[k];

                        g_xxy_xxxyz[k] = g1_xy_xxxxyz[k] - fr * g2_xy_xxxyz[k];

                        g_xxy_xxxzz[k] = g1_xy_xxxxzz[k] - fr * g2_xy_xxxzz[k];

                        g_xxy_xxyyy[k] = g1_xy_xxxyyy[k] - fr * g2_xy_xxyyy[k];

                        g_xxy_xxyyz[k] = g1_xy_xxxyyz[k] - fr * g2_xy_xxyyz[k];

                        g_xxy_xxyzz[k] = g1_xy_xxxyzz[k] - fr * g2_xy_xxyzz[k];

                        g_xxy_xxzzz[k] = g1_xy_xxxzzz[k] - fr * g2_xy_xxzzz[k];

                        g_xxy_xyyyy[k] = g1_xy_xxyyyy[k] - fr * g2_xy_xyyyy[k];

                        g_xxy_xyyyz[k] = g1_xy_xxyyyz[k] - fr * g2_xy_xyyyz[k];

                        g_xxy_xyyzz[k] = g1_xy_xxyyzz[k] - fr * g2_xy_xyyzz[k];

                        g_xxy_xyzzz[k] = g1_xy_xxyzzz[k] - fr * g2_xy_xyzzz[k];

                        g_xxy_xzzzz[k] = g1_xy_xxzzzz[k] - fr * g2_xy_xzzzz[k];

                        g_xxy_yyyyy[k] = g1_xy_xyyyyy[k] - fr * g2_xy_yyyyy[k];

                        g_xxy_yyyyz[k] = g1_xy_xyyyyz[k] - fr * g2_xy_yyyyz[k];

                        g_xxy_yyyzz[k] = g1_xy_xyyyzz[k] - fr * g2_xy_yyyzz[k];

                        g_xxy_yyzzz[k] = g1_xy_xyyzzz[k] - fr * g2_xy_yyzzz[k];

                        g_xxy_yzzzz[k] = g1_xy_xyzzzz[k] - fr * g2_xy_yzzzz[k];

                        g_xxy_zzzzz[k] = g1_xy_xzzzzz[k] - fr * g2_xy_zzzzz[k];

                        g_xxz_xxxxx[k] = g1_xz_xxxxxx[k] - fr * g2_xz_xxxxx[k];

                        g_xxz_xxxxy[k] = g1_xz_xxxxxy[k] - fr * g2_xz_xxxxy[k];

                        g_xxz_xxxxz[k] = g1_xz_xxxxxz[k] - fr * g2_xz_xxxxz[k];

                        g_xxz_xxxyy[k] = g1_xz_xxxxyy[k] - fr * g2_xz_xxxyy[k];

                        g_xxz_xxxyz[k] = g1_xz_xxxxyz[k] - fr * g2_xz_xxxyz[k];

                        g_xxz_xxxzz[k] = g1_xz_xxxxzz[k] - fr * g2_xz_xxxzz[k];

                        g_xxz_xxyyy[k] = g1_xz_xxxyyy[k] - fr * g2_xz_xxyyy[k];

                        g_xxz_xxyyz[k] = g1_xz_xxxyyz[k] - fr * g2_xz_xxyyz[k];

                        g_xxz_xxyzz[k] = g1_xz_xxxyzz[k] - fr * g2_xz_xxyzz[k];

                        g_xxz_xxzzz[k] = g1_xz_xxxzzz[k] - fr * g2_xz_xxzzz[k];

                        g_xxz_xyyyy[k] = g1_xz_xxyyyy[k] - fr * g2_xz_xyyyy[k];

                        g_xxz_xyyyz[k] = g1_xz_xxyyyz[k] - fr * g2_xz_xyyyz[k];

                        g_xxz_xyyzz[k] = g1_xz_xxyyzz[k] - fr * g2_xz_xyyzz[k];

                        g_xxz_xyzzz[k] = g1_xz_xxyzzz[k] - fr * g2_xz_xyzzz[k];

                        g_xxz_xzzzz[k] = g1_xz_xxzzzz[k] - fr * g2_xz_xzzzz[k];

                        g_xxz_yyyyy[k] = g1_xz_xyyyyy[k] - fr * g2_xz_yyyyy[k];

                        g_xxz_yyyyz[k] = g1_xz_xyyyyz[k] - fr * g2_xz_yyyyz[k];

                        g_xxz_yyyzz[k] = g1_xz_xyyyzz[k] - fr * g2_xz_yyyzz[k];

                        g_xxz_yyzzz[k] = g1_xz_xyyzzz[k] - fr * g2_xz_yyzzz[k];

                        g_xxz_yzzzz[k] = g1_xz_xyzzzz[k] - fr * g2_xz_yzzzz[k];

                        g_xxz_zzzzz[k] = g1_xz_xzzzzz[k] - fr * g2_xz_zzzzz[k];

                        g_xyy_xxxxx[k] = g1_yy_xxxxxx[k] - fr * g2_yy_xxxxx[k];

                        g_xyy_xxxxy[k] = g1_yy_xxxxxy[k] - fr * g2_yy_xxxxy[k];

                        g_xyy_xxxxz[k] = g1_yy_xxxxxz[k] - fr * g2_yy_xxxxz[k];

                        g_xyy_xxxyy[k] = g1_yy_xxxxyy[k] - fr * g2_yy_xxxyy[k];

                        g_xyy_xxxyz[k] = g1_yy_xxxxyz[k] - fr * g2_yy_xxxyz[k];

                        g_xyy_xxxzz[k] = g1_yy_xxxxzz[k] - fr * g2_yy_xxxzz[k];

                        g_xyy_xxyyy[k] = g1_yy_xxxyyy[k] - fr * g2_yy_xxyyy[k];

                        g_xyy_xxyyz[k] = g1_yy_xxxyyz[k] - fr * g2_yy_xxyyz[k];

                        g_xyy_xxyzz[k] = g1_yy_xxxyzz[k] - fr * g2_yy_xxyzz[k];

                        g_xyy_xxzzz[k] = g1_yy_xxxzzz[k] - fr * g2_yy_xxzzz[k];

                        g_xyy_xyyyy[k] = g1_yy_xxyyyy[k] - fr * g2_yy_xyyyy[k];

                        g_xyy_xyyyz[k] = g1_yy_xxyyyz[k] - fr * g2_yy_xyyyz[k];

                        g_xyy_xyyzz[k] = g1_yy_xxyyzz[k] - fr * g2_yy_xyyzz[k];

                        g_xyy_xyzzz[k] = g1_yy_xxyzzz[k] - fr * g2_yy_xyzzz[k];

                        g_xyy_xzzzz[k] = g1_yy_xxzzzz[k] - fr * g2_yy_xzzzz[k];

                        g_xyy_yyyyy[k] = g1_yy_xyyyyy[k] - fr * g2_yy_yyyyy[k];

                        g_xyy_yyyyz[k] = g1_yy_xyyyyz[k] - fr * g2_yy_yyyyz[k];

                        g_xyy_yyyzz[k] = g1_yy_xyyyzz[k] - fr * g2_yy_yyyzz[k];

                        g_xyy_yyzzz[k] = g1_yy_xyyzzz[k] - fr * g2_yy_yyzzz[k];

                        g_xyy_yzzzz[k] = g1_yy_xyzzzz[k] - fr * g2_yy_yzzzz[k];

                        g_xyy_zzzzz[k] = g1_yy_xzzzzz[k] - fr * g2_yy_zzzzz[k];

                        g_xyz_xxxxx[k] = g1_yz_xxxxxx[k] - fr * g2_yz_xxxxx[k];

                        g_xyz_xxxxy[k] = g1_yz_xxxxxy[k] - fr * g2_yz_xxxxy[k];

                        g_xyz_xxxxz[k] = g1_yz_xxxxxz[k] - fr * g2_yz_xxxxz[k];

                        g_xyz_xxxyy[k] = g1_yz_xxxxyy[k] - fr * g2_yz_xxxyy[k];

                        g_xyz_xxxyz[k] = g1_yz_xxxxyz[k] - fr * g2_yz_xxxyz[k];

                        g_xyz_xxxzz[k] = g1_yz_xxxxzz[k] - fr * g2_yz_xxxzz[k];

                        g_xyz_xxyyy[k] = g1_yz_xxxyyy[k] - fr * g2_yz_xxyyy[k];

                        g_xyz_xxyyz[k] = g1_yz_xxxyyz[k] - fr * g2_yz_xxyyz[k];

                        g_xyz_xxyzz[k] = g1_yz_xxxyzz[k] - fr * g2_yz_xxyzz[k];

                        g_xyz_xxzzz[k] = g1_yz_xxxzzz[k] - fr * g2_yz_xxzzz[k];

                        g_xyz_xyyyy[k] = g1_yz_xxyyyy[k] - fr * g2_yz_xyyyy[k];

                        g_xyz_xyyyz[k] = g1_yz_xxyyyz[k] - fr * g2_yz_xyyyz[k];

                        g_xyz_xyyzz[k] = g1_yz_xxyyzz[k] - fr * g2_yz_xyyzz[k];

                        g_xyz_xyzzz[k] = g1_yz_xxyzzz[k] - fr * g2_yz_xyzzz[k];

                        g_xyz_xzzzz[k] = g1_yz_xxzzzz[k] - fr * g2_yz_xzzzz[k];

                        g_xyz_yyyyy[k] = g1_yz_xyyyyy[k] - fr * g2_yz_yyyyy[k];

                        g_xyz_yyyyz[k] = g1_yz_xyyyyz[k] - fr * g2_yz_yyyyz[k];

                        g_xyz_yyyzz[k] = g1_yz_xyyyzz[k] - fr * g2_yz_yyyzz[k];

                        g_xyz_yyzzz[k] = g1_yz_xyyzzz[k] - fr * g2_yz_yyzzz[k];

                        g_xyz_yzzzz[k] = g1_yz_xyzzzz[k] - fr * g2_yz_yzzzz[k];

                        g_xyz_zzzzz[k] = g1_yz_xzzzzz[k] - fr * g2_yz_zzzzz[k];

                        g_xzz_xxxxx[k] = g1_zz_xxxxxx[k] - fr * g2_zz_xxxxx[k];

                        g_xzz_xxxxy[k] = g1_zz_xxxxxy[k] - fr * g2_zz_xxxxy[k];

                        g_xzz_xxxxz[k] = g1_zz_xxxxxz[k] - fr * g2_zz_xxxxz[k];

                        g_xzz_xxxyy[k] = g1_zz_xxxxyy[k] - fr * g2_zz_xxxyy[k];

                        g_xzz_xxxyz[k] = g1_zz_xxxxyz[k] - fr * g2_zz_xxxyz[k];

                        g_xzz_xxxzz[k] = g1_zz_xxxxzz[k] - fr * g2_zz_xxxzz[k];

                        g_xzz_xxyyy[k] = g1_zz_xxxyyy[k] - fr * g2_zz_xxyyy[k];

                        g_xzz_xxyyz[k] = g1_zz_xxxyyz[k] - fr * g2_zz_xxyyz[k];

                        g_xzz_xxyzz[k] = g1_zz_xxxyzz[k] - fr * g2_zz_xxyzz[k];

                        g_xzz_xxzzz[k] = g1_zz_xxxzzz[k] - fr * g2_zz_xxzzz[k];

                        g_xzz_xyyyy[k] = g1_zz_xxyyyy[k] - fr * g2_zz_xyyyy[k];

                        g_xzz_xyyyz[k] = g1_zz_xxyyyz[k] - fr * g2_zz_xyyyz[k];

                        g_xzz_xyyzz[k] = g1_zz_xxyyzz[k] - fr * g2_zz_xyyzz[k];

                        g_xzz_xyzzz[k] = g1_zz_xxyzzz[k] - fr * g2_zz_xyzzz[k];

                        g_xzz_xzzzz[k] = g1_zz_xxzzzz[k] - fr * g2_zz_xzzzz[k];

                        g_xzz_yyyyy[k] = g1_zz_xyyyyy[k] - fr * g2_zz_yyyyy[k];

                        g_xzz_yyyyz[k] = g1_zz_xyyyyz[k] - fr * g2_zz_yyyyz[k];

                        g_xzz_yyyzz[k] = g1_zz_xyyyzz[k] - fr * g2_zz_yyyzz[k];

                        g_xzz_yyzzz[k] = g1_zz_xyyzzz[k] - fr * g2_zz_yyzzz[k];

                        g_xzz_yzzzz[k] = g1_zz_xyzzzz[k] - fr * g2_zz_yzzzz[k];

                        g_xzz_zzzzz[k] = g1_zz_xzzzzz[k] - fr * g2_zz_zzzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yyy_xxxxx[k] = g1_yy_xxxxxy[k] - fr * g2_yy_xxxxx[k];

                        g_yyy_xxxxy[k] = g1_yy_xxxxyy[k] - fr * g2_yy_xxxxy[k];

                        g_yyy_xxxxz[k] = g1_yy_xxxxyz[k] - fr * g2_yy_xxxxz[k];

                        g_yyy_xxxyy[k] = g1_yy_xxxyyy[k] - fr * g2_yy_xxxyy[k];

                        g_yyy_xxxyz[k] = g1_yy_xxxyyz[k] - fr * g2_yy_xxxyz[k];

                        g_yyy_xxxzz[k] = g1_yy_xxxyzz[k] - fr * g2_yy_xxxzz[k];

                        g_yyy_xxyyy[k] = g1_yy_xxyyyy[k] - fr * g2_yy_xxyyy[k];

                        g_yyy_xxyyz[k] = g1_yy_xxyyyz[k] - fr * g2_yy_xxyyz[k];

                        g_yyy_xxyzz[k] = g1_yy_xxyyzz[k] - fr * g2_yy_xxyzz[k];

                        g_yyy_xxzzz[k] = g1_yy_xxyzzz[k] - fr * g2_yy_xxzzz[k];

                        g_yyy_xyyyy[k] = g1_yy_xyyyyy[k] - fr * g2_yy_xyyyy[k];

                        g_yyy_xyyyz[k] = g1_yy_xyyyyz[k] - fr * g2_yy_xyyyz[k];

                        g_yyy_xyyzz[k] = g1_yy_xyyyzz[k] - fr * g2_yy_xyyzz[k];

                        g_yyy_xyzzz[k] = g1_yy_xyyzzz[k] - fr * g2_yy_xyzzz[k];

                        g_yyy_xzzzz[k] = g1_yy_xyzzzz[k] - fr * g2_yy_xzzzz[k];

                        g_yyy_yyyyy[k] = g1_yy_yyyyyy[k] - fr * g2_yy_yyyyy[k];

                        g_yyy_yyyyz[k] = g1_yy_yyyyyz[k] - fr * g2_yy_yyyyz[k];

                        g_yyy_yyyzz[k] = g1_yy_yyyyzz[k] - fr * g2_yy_yyyzz[k];

                        g_yyy_yyzzz[k] = g1_yy_yyyzzz[k] - fr * g2_yy_yyzzz[k];

                        g_yyy_yzzzz[k] = g1_yy_yyzzzz[k] - fr * g2_yy_yzzzz[k];

                        g_yyy_zzzzz[k] = g1_yy_yzzzzz[k] - fr * g2_yy_zzzzz[k];

                        g_yyz_xxxxx[k] = g1_yz_xxxxxy[k] - fr * g2_yz_xxxxx[k];

                        g_yyz_xxxxy[k] = g1_yz_xxxxyy[k] - fr * g2_yz_xxxxy[k];

                        g_yyz_xxxxz[k] = g1_yz_xxxxyz[k] - fr * g2_yz_xxxxz[k];

                        g_yyz_xxxyy[k] = g1_yz_xxxyyy[k] - fr * g2_yz_xxxyy[k];

                        g_yyz_xxxyz[k] = g1_yz_xxxyyz[k] - fr * g2_yz_xxxyz[k];

                        g_yyz_xxxzz[k] = g1_yz_xxxyzz[k] - fr * g2_yz_xxxzz[k];

                        g_yyz_xxyyy[k] = g1_yz_xxyyyy[k] - fr * g2_yz_xxyyy[k];

                        g_yyz_xxyyz[k] = g1_yz_xxyyyz[k] - fr * g2_yz_xxyyz[k];

                        g_yyz_xxyzz[k] = g1_yz_xxyyzz[k] - fr * g2_yz_xxyzz[k];

                        g_yyz_xxzzz[k] = g1_yz_xxyzzz[k] - fr * g2_yz_xxzzz[k];

                        g_yyz_xyyyy[k] = g1_yz_xyyyyy[k] - fr * g2_yz_xyyyy[k];

                        g_yyz_xyyyz[k] = g1_yz_xyyyyz[k] - fr * g2_yz_xyyyz[k];

                        g_yyz_xyyzz[k] = g1_yz_xyyyzz[k] - fr * g2_yz_xyyzz[k];

                        g_yyz_xyzzz[k] = g1_yz_xyyzzz[k] - fr * g2_yz_xyzzz[k];

                        g_yyz_xzzzz[k] = g1_yz_xyzzzz[k] - fr * g2_yz_xzzzz[k];

                        g_yyz_yyyyy[k] = g1_yz_yyyyyy[k] - fr * g2_yz_yyyyy[k];

                        g_yyz_yyyyz[k] = g1_yz_yyyyyz[k] - fr * g2_yz_yyyyz[k];

                        g_yyz_yyyzz[k] = g1_yz_yyyyzz[k] - fr * g2_yz_yyyzz[k];

                        g_yyz_yyzzz[k] = g1_yz_yyyzzz[k] - fr * g2_yz_yyzzz[k];

                        g_yyz_yzzzz[k] = g1_yz_yyzzzz[k] - fr * g2_yz_yzzzz[k];

                        g_yyz_zzzzz[k] = g1_yz_yzzzzz[k] - fr * g2_yz_zzzzz[k];

                        g_yzz_xxxxx[k] = g1_zz_xxxxxy[k] - fr * g2_zz_xxxxx[k];

                        g_yzz_xxxxy[k] = g1_zz_xxxxyy[k] - fr * g2_zz_xxxxy[k];

                        g_yzz_xxxxz[k] = g1_zz_xxxxyz[k] - fr * g2_zz_xxxxz[k];

                        g_yzz_xxxyy[k] = g1_zz_xxxyyy[k] - fr * g2_zz_xxxyy[k];

                        g_yzz_xxxyz[k] = g1_zz_xxxyyz[k] - fr * g2_zz_xxxyz[k];

                        g_yzz_xxxzz[k] = g1_zz_xxxyzz[k] - fr * g2_zz_xxxzz[k];

                        g_yzz_xxyyy[k] = g1_zz_xxyyyy[k] - fr * g2_zz_xxyyy[k];

                        g_yzz_xxyyz[k] = g1_zz_xxyyyz[k] - fr * g2_zz_xxyyz[k];

                        g_yzz_xxyzz[k] = g1_zz_xxyyzz[k] - fr * g2_zz_xxyzz[k];

                        g_yzz_xxzzz[k] = g1_zz_xxyzzz[k] - fr * g2_zz_xxzzz[k];

                        g_yzz_xyyyy[k] = g1_zz_xyyyyy[k] - fr * g2_zz_xyyyy[k];

                        g_yzz_xyyyz[k] = g1_zz_xyyyyz[k] - fr * g2_zz_xyyyz[k];

                        g_yzz_xyyzz[k] = g1_zz_xyyyzz[k] - fr * g2_zz_xyyzz[k];

                        g_yzz_xyzzz[k] = g1_zz_xyyzzz[k] - fr * g2_zz_xyzzz[k];

                        g_yzz_xzzzz[k] = g1_zz_xyzzzz[k] - fr * g2_zz_xzzzz[k];

                        g_yzz_yyyyy[k] = g1_zz_yyyyyy[k] - fr * g2_zz_yyyyy[k];

                        g_yzz_yyyyz[k] = g1_zz_yyyyyz[k] - fr * g2_zz_yyyyz[k];

                        g_yzz_yyyzz[k] = g1_zz_yyyyzz[k] - fr * g2_zz_yyyzz[k];

                        g_yzz_yyzzz[k] = g1_zz_yyyzzz[k] - fr * g2_zz_yyzzz[k];

                        g_yzz_yzzzz[k] = g1_zz_yyzzzz[k] - fr * g2_zz_yzzzz[k];

                        g_yzz_zzzzz[k] = g1_zz_yzzzzz[k] - fr * g2_zz_zzzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zzz_xxxxx[k] = g1_zz_xxxxxz[k] - fr * g2_zz_xxxxx[k];

                        g_zzz_xxxxy[k] = g1_zz_xxxxyz[k] - fr * g2_zz_xxxxy[k];

                        g_zzz_xxxxz[k] = g1_zz_xxxxzz[k] - fr * g2_zz_xxxxz[k];

                        g_zzz_xxxyy[k] = g1_zz_xxxyyz[k] - fr * g2_zz_xxxyy[k];

                        g_zzz_xxxyz[k] = g1_zz_xxxyzz[k] - fr * g2_zz_xxxyz[k];

                        g_zzz_xxxzz[k] = g1_zz_xxxzzz[k] - fr * g2_zz_xxxzz[k];

                        g_zzz_xxyyy[k] = g1_zz_xxyyyz[k] - fr * g2_zz_xxyyy[k];

                        g_zzz_xxyyz[k] = g1_zz_xxyyzz[k] - fr * g2_zz_xxyyz[k];

                        g_zzz_xxyzz[k] = g1_zz_xxyzzz[k] - fr * g2_zz_xxyzz[k];

                        g_zzz_xxzzz[k] = g1_zz_xxzzzz[k] - fr * g2_zz_xxzzz[k];

                        g_zzz_xyyyy[k] = g1_zz_xyyyyz[k] - fr * g2_zz_xyyyy[k];

                        g_zzz_xyyyz[k] = g1_zz_xyyyzz[k] - fr * g2_zz_xyyyz[k];

                        g_zzz_xyyzz[k] = g1_zz_xyyzzz[k] - fr * g2_zz_xyyzz[k];

                        g_zzz_xyzzz[k] = g1_zz_xyzzzz[k] - fr * g2_zz_xyzzz[k];

                        g_zzz_xzzzz[k] = g1_zz_xzzzzz[k] - fr * g2_zz_xzzzz[k];

                        g_zzz_yyyyy[k] = g1_zz_yyyyyz[k] - fr * g2_zz_yyyyy[k];

                        g_zzz_yyyyz[k] = g1_zz_yyyyzz[k] - fr * g2_zz_yyyyz[k];

                        g_zzz_yyyzz[k] = g1_zz_yyyzzz[k] - fr * g2_zz_yyyzz[k];

                        g_zzz_yyzzz[k] = g1_zz_yyzzzz[k] - fr * g2_zz_yyzzz[k];

                        g_zzz_yzzzz[k] = g1_zz_yzzzzz[k] - fr * g2_zz_yzzzz[k];

                        g_zzz_zzzzz[k] = g1_zz_zzzzzz[k] - fr * g2_zz_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForSXGG(      CMemBlock2D<double>&  contrBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  cdDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set up pointers to distances R(CD) = C - D

        auto rcdx = cdDistances.data(0);

        auto rcdy = cdDistances.data(1);

        auto rcdz = cdDistances.data(2);

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].second() == 4) && (recPattern[i].third() == 4))
            {
                if (iContrPair == 0) printf("-> applying ket HRR for (0X|44)\n");

                // determine angular momentum of bra side

                auto bang  = recPattern[i].first();

                auto bcomp = angmom::to_CartesianComponents(bang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {bang, 4, 4});

                auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 3, 5});

                auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                      {bang, 3, 4});

                // compute contracted integrals

                for (int32_t j = 0; j < bcomp; j++)
                {
                    // set up pointers to (SX|g(r,r')|FG)^(m) integrals

                    auto g2_xxx_xxxx = contrBuffer.data(g2off + 150 * j);

                    auto g2_xxx_xxxy = contrBuffer.data(g2off + 150 * j + 1);

                    auto g2_xxx_xxxz = contrBuffer.data(g2off + 150 * j + 2);

                    auto g2_xxx_xxyy = contrBuffer.data(g2off + 150 * j + 3);

                    auto g2_xxx_xxyz = contrBuffer.data(g2off + 150 * j + 4);

                    auto g2_xxx_xxzz = contrBuffer.data(g2off + 150 * j + 5);

                    auto g2_xxx_xyyy = contrBuffer.data(g2off + 150 * j + 6);

                    auto g2_xxx_xyyz = contrBuffer.data(g2off + 150 * j + 7);

                    auto g2_xxx_xyzz = contrBuffer.data(g2off + 150 * j + 8);

                    auto g2_xxx_xzzz = contrBuffer.data(g2off + 150 * j + 9);

                    auto g2_xxx_yyyy = contrBuffer.data(g2off + 150 * j + 10);

                    auto g2_xxx_yyyz = contrBuffer.data(g2off + 150 * j + 11);

                    auto g2_xxx_yyzz = contrBuffer.data(g2off + 150 * j + 12);

                    auto g2_xxx_yzzz = contrBuffer.data(g2off + 150 * j + 13);

                    auto g2_xxx_zzzz = contrBuffer.data(g2off + 150 * j + 14);

                    auto g2_xxy_xxxx = contrBuffer.data(g2off + 150 * j + 15);

                    auto g2_xxy_xxxy = contrBuffer.data(g2off + 150 * j + 16);

                    auto g2_xxy_xxxz = contrBuffer.data(g2off + 150 * j + 17);

                    auto g2_xxy_xxyy = contrBuffer.data(g2off + 150 * j + 18);

                    auto g2_xxy_xxyz = contrBuffer.data(g2off + 150 * j + 19);

                    auto g2_xxy_xxzz = contrBuffer.data(g2off + 150 * j + 20);

                    auto g2_xxy_xyyy = contrBuffer.data(g2off + 150 * j + 21);

                    auto g2_xxy_xyyz = contrBuffer.data(g2off + 150 * j + 22);

                    auto g2_xxy_xyzz = contrBuffer.data(g2off + 150 * j + 23);

                    auto g2_xxy_xzzz = contrBuffer.data(g2off + 150 * j + 24);

                    auto g2_xxy_yyyy = contrBuffer.data(g2off + 150 * j + 25);

                    auto g2_xxy_yyyz = contrBuffer.data(g2off + 150 * j + 26);

                    auto g2_xxy_yyzz = contrBuffer.data(g2off + 150 * j + 27);

                    auto g2_xxy_yzzz = contrBuffer.data(g2off + 150 * j + 28);

                    auto g2_xxy_zzzz = contrBuffer.data(g2off + 150 * j + 29);

                    auto g2_xxz_xxxx = contrBuffer.data(g2off + 150 * j + 30);

                    auto g2_xxz_xxxy = contrBuffer.data(g2off + 150 * j + 31);

                    auto g2_xxz_xxxz = contrBuffer.data(g2off + 150 * j + 32);

                    auto g2_xxz_xxyy = contrBuffer.data(g2off + 150 * j + 33);

                    auto g2_xxz_xxyz = contrBuffer.data(g2off + 150 * j + 34);

                    auto g2_xxz_xxzz = contrBuffer.data(g2off + 150 * j + 35);

                    auto g2_xxz_xyyy = contrBuffer.data(g2off + 150 * j + 36);

                    auto g2_xxz_xyyz = contrBuffer.data(g2off + 150 * j + 37);

                    auto g2_xxz_xyzz = contrBuffer.data(g2off + 150 * j + 38);

                    auto g2_xxz_xzzz = contrBuffer.data(g2off + 150 * j + 39);

                    auto g2_xxz_yyyy = contrBuffer.data(g2off + 150 * j + 40);

                    auto g2_xxz_yyyz = contrBuffer.data(g2off + 150 * j + 41);

                    auto g2_xxz_yyzz = contrBuffer.data(g2off + 150 * j + 42);

                    auto g2_xxz_yzzz = contrBuffer.data(g2off + 150 * j + 43);

                    auto g2_xxz_zzzz = contrBuffer.data(g2off + 150 * j + 44);

                    auto g2_xyy_xxxx = contrBuffer.data(g2off + 150 * j + 45);

                    auto g2_xyy_xxxy = contrBuffer.data(g2off + 150 * j + 46);

                    auto g2_xyy_xxxz = contrBuffer.data(g2off + 150 * j + 47);

                    auto g2_xyy_xxyy = contrBuffer.data(g2off + 150 * j + 48);

                    auto g2_xyy_xxyz = contrBuffer.data(g2off + 150 * j + 49);

                    auto g2_xyy_xxzz = contrBuffer.data(g2off + 150 * j + 50);

                    auto g2_xyy_xyyy = contrBuffer.data(g2off + 150 * j + 51);

                    auto g2_xyy_xyyz = contrBuffer.data(g2off + 150 * j + 52);

                    auto g2_xyy_xyzz = contrBuffer.data(g2off + 150 * j + 53);

                    auto g2_xyy_xzzz = contrBuffer.data(g2off + 150 * j + 54);

                    auto g2_xyy_yyyy = contrBuffer.data(g2off + 150 * j + 55);

                    auto g2_xyy_yyyz = contrBuffer.data(g2off + 150 * j + 56);

                    auto g2_xyy_yyzz = contrBuffer.data(g2off + 150 * j + 57);

                    auto g2_xyy_yzzz = contrBuffer.data(g2off + 150 * j + 58);

                    auto g2_xyy_zzzz = contrBuffer.data(g2off + 150 * j + 59);

                    auto g2_xyz_xxxx = contrBuffer.data(g2off + 150 * j + 60);

                    auto g2_xyz_xxxy = contrBuffer.data(g2off + 150 * j + 61);

                    auto g2_xyz_xxxz = contrBuffer.data(g2off + 150 * j + 62);

                    auto g2_xyz_xxyy = contrBuffer.data(g2off + 150 * j + 63);

                    auto g2_xyz_xxyz = contrBuffer.data(g2off + 150 * j + 64);

                    auto g2_xyz_xxzz = contrBuffer.data(g2off + 150 * j + 65);

                    auto g2_xyz_xyyy = contrBuffer.data(g2off + 150 * j + 66);

                    auto g2_xyz_xyyz = contrBuffer.data(g2off + 150 * j + 67);

                    auto g2_xyz_xyzz = contrBuffer.data(g2off + 150 * j + 68);

                    auto g2_xyz_xzzz = contrBuffer.data(g2off + 150 * j + 69);

                    auto g2_xyz_yyyy = contrBuffer.data(g2off + 150 * j + 70);

                    auto g2_xyz_yyyz = contrBuffer.data(g2off + 150 * j + 71);

                    auto g2_xyz_yyzz = contrBuffer.data(g2off + 150 * j + 72);

                    auto g2_xyz_yzzz = contrBuffer.data(g2off + 150 * j + 73);

                    auto g2_xyz_zzzz = contrBuffer.data(g2off + 150 * j + 74);

                    auto g2_xzz_xxxx = contrBuffer.data(g2off + 150 * j + 75);

                    auto g2_xzz_xxxy = contrBuffer.data(g2off + 150 * j + 76);

                    auto g2_xzz_xxxz = contrBuffer.data(g2off + 150 * j + 77);

                    auto g2_xzz_xxyy = contrBuffer.data(g2off + 150 * j + 78);

                    auto g2_xzz_xxyz = contrBuffer.data(g2off + 150 * j + 79);

                    auto g2_xzz_xxzz = contrBuffer.data(g2off + 150 * j + 80);

                    auto g2_xzz_xyyy = contrBuffer.data(g2off + 150 * j + 81);

                    auto g2_xzz_xyyz = contrBuffer.data(g2off + 150 * j + 82);

                    auto g2_xzz_xyzz = contrBuffer.data(g2off + 150 * j + 83);

                    auto g2_xzz_xzzz = contrBuffer.data(g2off + 150 * j + 84);

                    auto g2_xzz_yyyy = contrBuffer.data(g2off + 150 * j + 85);

                    auto g2_xzz_yyyz = contrBuffer.data(g2off + 150 * j + 86);

                    auto g2_xzz_yyzz = contrBuffer.data(g2off + 150 * j + 87);

                    auto g2_xzz_yzzz = contrBuffer.data(g2off + 150 * j + 88);

                    auto g2_xzz_zzzz = contrBuffer.data(g2off + 150 * j + 89);

                    auto g2_yyy_xxxx = contrBuffer.data(g2off + 150 * j + 90);

                    auto g2_yyy_xxxy = contrBuffer.data(g2off + 150 * j + 91);

                    auto g2_yyy_xxxz = contrBuffer.data(g2off + 150 * j + 92);

                    auto g2_yyy_xxyy = contrBuffer.data(g2off + 150 * j + 93);

                    auto g2_yyy_xxyz = contrBuffer.data(g2off + 150 * j + 94);

                    auto g2_yyy_xxzz = contrBuffer.data(g2off + 150 * j + 95);

                    auto g2_yyy_xyyy = contrBuffer.data(g2off + 150 * j + 96);

                    auto g2_yyy_xyyz = contrBuffer.data(g2off + 150 * j + 97);

                    auto g2_yyy_xyzz = contrBuffer.data(g2off + 150 * j + 98);

                    auto g2_yyy_xzzz = contrBuffer.data(g2off + 150 * j + 99);

                    auto g2_yyy_yyyy = contrBuffer.data(g2off + 150 * j + 100);

                    auto g2_yyy_yyyz = contrBuffer.data(g2off + 150 * j + 101);

                    auto g2_yyy_yyzz = contrBuffer.data(g2off + 150 * j + 102);

                    auto g2_yyy_yzzz = contrBuffer.data(g2off + 150 * j + 103);

                    auto g2_yyy_zzzz = contrBuffer.data(g2off + 150 * j + 104);

                    auto g2_yyz_xxxx = contrBuffer.data(g2off + 150 * j + 105);

                    auto g2_yyz_xxxy = contrBuffer.data(g2off + 150 * j + 106);

                    auto g2_yyz_xxxz = contrBuffer.data(g2off + 150 * j + 107);

                    auto g2_yyz_xxyy = contrBuffer.data(g2off + 150 * j + 108);

                    auto g2_yyz_xxyz = contrBuffer.data(g2off + 150 * j + 109);

                    auto g2_yyz_xxzz = contrBuffer.data(g2off + 150 * j + 110);

                    auto g2_yyz_xyyy = contrBuffer.data(g2off + 150 * j + 111);

                    auto g2_yyz_xyyz = contrBuffer.data(g2off + 150 * j + 112);

                    auto g2_yyz_xyzz = contrBuffer.data(g2off + 150 * j + 113);

                    auto g2_yyz_xzzz = contrBuffer.data(g2off + 150 * j + 114);

                    auto g2_yyz_yyyy = contrBuffer.data(g2off + 150 * j + 115);

                    auto g2_yyz_yyyz = contrBuffer.data(g2off + 150 * j + 116);

                    auto g2_yyz_yyzz = contrBuffer.data(g2off + 150 * j + 117);

                    auto g2_yyz_yzzz = contrBuffer.data(g2off + 150 * j + 118);

                    auto g2_yyz_zzzz = contrBuffer.data(g2off + 150 * j + 119);

                    auto g2_yzz_xxxx = contrBuffer.data(g2off + 150 * j + 120);

                    auto g2_yzz_xxxy = contrBuffer.data(g2off + 150 * j + 121);

                    auto g2_yzz_xxxz = contrBuffer.data(g2off + 150 * j + 122);

                    auto g2_yzz_xxyy = contrBuffer.data(g2off + 150 * j + 123);

                    auto g2_yzz_xxyz = contrBuffer.data(g2off + 150 * j + 124);

                    auto g2_yzz_xxzz = contrBuffer.data(g2off + 150 * j + 125);

                    auto g2_yzz_xyyy = contrBuffer.data(g2off + 150 * j + 126);

                    auto g2_yzz_xyyz = contrBuffer.data(g2off + 150 * j + 127);

                    auto g2_yzz_xyzz = contrBuffer.data(g2off + 150 * j + 128);

                    auto g2_yzz_xzzz = contrBuffer.data(g2off + 150 * j + 129);

                    auto g2_yzz_yyyy = contrBuffer.data(g2off + 150 * j + 130);

                    auto g2_yzz_yyyz = contrBuffer.data(g2off + 150 * j + 131);

                    auto g2_yzz_yyzz = contrBuffer.data(g2off + 150 * j + 132);

                    auto g2_yzz_yzzz = contrBuffer.data(g2off + 150 * j + 133);

                    auto g2_yzz_zzzz = contrBuffer.data(g2off + 150 * j + 134);

                    auto g2_zzz_xxxx = contrBuffer.data(g2off + 150 * j + 135);

                    auto g2_zzz_xxxy = contrBuffer.data(g2off + 150 * j + 136);

                    auto g2_zzz_xxxz = contrBuffer.data(g2off + 150 * j + 137);

                    auto g2_zzz_xxyy = contrBuffer.data(g2off + 150 * j + 138);

                    auto g2_zzz_xxyz = contrBuffer.data(g2off + 150 * j + 139);

                    auto g2_zzz_xxzz = contrBuffer.data(g2off + 150 * j + 140);

                    auto g2_zzz_xyyy = contrBuffer.data(g2off + 150 * j + 141);

                    auto g2_zzz_xyyz = contrBuffer.data(g2off + 150 * j + 142);

                    auto g2_zzz_xyzz = contrBuffer.data(g2off + 150 * j + 143);

                    auto g2_zzz_xzzz = contrBuffer.data(g2off + 150 * j + 144);

                    auto g2_zzz_yyyy = contrBuffer.data(g2off + 150 * j + 145);

                    auto g2_zzz_yyyz = contrBuffer.data(g2off + 150 * j + 146);

                    auto g2_zzz_yyzz = contrBuffer.data(g2off + 150 * j + 147);

                    auto g2_zzz_yzzz = contrBuffer.data(g2off + 150 * j + 148);

                    auto g2_zzz_zzzz = contrBuffer.data(g2off + 150 * j + 149);

                    // set up pointers to (SX|g(r,r')|FH)^(m) integrals

                    auto g1_xxx_xxxxx = contrBuffer.data(g1off + 210 * j);

                    auto g1_xxx_xxxxy = contrBuffer.data(g1off + 210 * j + 1);

                    auto g1_xxx_xxxxz = contrBuffer.data(g1off + 210 * j + 2);

                    auto g1_xxx_xxxyy = contrBuffer.data(g1off + 210 * j + 3);

                    auto g1_xxx_xxxyz = contrBuffer.data(g1off + 210 * j + 4);

                    auto g1_xxx_xxxzz = contrBuffer.data(g1off + 210 * j + 5);

                    auto g1_xxx_xxyyy = contrBuffer.data(g1off + 210 * j + 6);

                    auto g1_xxx_xxyyz = contrBuffer.data(g1off + 210 * j + 7);

                    auto g1_xxx_xxyzz = contrBuffer.data(g1off + 210 * j + 8);

                    auto g1_xxx_xxzzz = contrBuffer.data(g1off + 210 * j + 9);

                    auto g1_xxx_xyyyy = contrBuffer.data(g1off + 210 * j + 10);

                    auto g1_xxx_xyyyz = contrBuffer.data(g1off + 210 * j + 11);

                    auto g1_xxx_xyyzz = contrBuffer.data(g1off + 210 * j + 12);

                    auto g1_xxx_xyzzz = contrBuffer.data(g1off + 210 * j + 13);

                    auto g1_xxx_xzzzz = contrBuffer.data(g1off + 210 * j + 14);

                    auto g1_xxy_xxxxx = contrBuffer.data(g1off + 210 * j + 21);

                    auto g1_xxy_xxxxy = contrBuffer.data(g1off + 210 * j + 22);

                    auto g1_xxy_xxxxz = contrBuffer.data(g1off + 210 * j + 23);

                    auto g1_xxy_xxxyy = contrBuffer.data(g1off + 210 * j + 24);

                    auto g1_xxy_xxxyz = contrBuffer.data(g1off + 210 * j + 25);

                    auto g1_xxy_xxxzz = contrBuffer.data(g1off + 210 * j + 26);

                    auto g1_xxy_xxyyy = contrBuffer.data(g1off + 210 * j + 27);

                    auto g1_xxy_xxyyz = contrBuffer.data(g1off + 210 * j + 28);

                    auto g1_xxy_xxyzz = contrBuffer.data(g1off + 210 * j + 29);

                    auto g1_xxy_xxzzz = contrBuffer.data(g1off + 210 * j + 30);

                    auto g1_xxy_xyyyy = contrBuffer.data(g1off + 210 * j + 31);

                    auto g1_xxy_xyyyz = contrBuffer.data(g1off + 210 * j + 32);

                    auto g1_xxy_xyyzz = contrBuffer.data(g1off + 210 * j + 33);

                    auto g1_xxy_xyzzz = contrBuffer.data(g1off + 210 * j + 34);

                    auto g1_xxy_xzzzz = contrBuffer.data(g1off + 210 * j + 35);

                    auto g1_xxz_xxxxx = contrBuffer.data(g1off + 210 * j + 42);

                    auto g1_xxz_xxxxy = contrBuffer.data(g1off + 210 * j + 43);

                    auto g1_xxz_xxxxz = contrBuffer.data(g1off + 210 * j + 44);

                    auto g1_xxz_xxxyy = contrBuffer.data(g1off + 210 * j + 45);

                    auto g1_xxz_xxxyz = contrBuffer.data(g1off + 210 * j + 46);

                    auto g1_xxz_xxxzz = contrBuffer.data(g1off + 210 * j + 47);

                    auto g1_xxz_xxyyy = contrBuffer.data(g1off + 210 * j + 48);

                    auto g1_xxz_xxyyz = contrBuffer.data(g1off + 210 * j + 49);

                    auto g1_xxz_xxyzz = contrBuffer.data(g1off + 210 * j + 50);

                    auto g1_xxz_xxzzz = contrBuffer.data(g1off + 210 * j + 51);

                    auto g1_xxz_xyyyy = contrBuffer.data(g1off + 210 * j + 52);

                    auto g1_xxz_xyyyz = contrBuffer.data(g1off + 210 * j + 53);

                    auto g1_xxz_xyyzz = contrBuffer.data(g1off + 210 * j + 54);

                    auto g1_xxz_xyzzz = contrBuffer.data(g1off + 210 * j + 55);

                    auto g1_xxz_xzzzz = contrBuffer.data(g1off + 210 * j + 56);

                    auto g1_xyy_xxxxx = contrBuffer.data(g1off + 210 * j + 63);

                    auto g1_xyy_xxxxy = contrBuffer.data(g1off + 210 * j + 64);

                    auto g1_xyy_xxxxz = contrBuffer.data(g1off + 210 * j + 65);

                    auto g1_xyy_xxxyy = contrBuffer.data(g1off + 210 * j + 66);

                    auto g1_xyy_xxxyz = contrBuffer.data(g1off + 210 * j + 67);

                    auto g1_xyy_xxxzz = contrBuffer.data(g1off + 210 * j + 68);

                    auto g1_xyy_xxyyy = contrBuffer.data(g1off + 210 * j + 69);

                    auto g1_xyy_xxyyz = contrBuffer.data(g1off + 210 * j + 70);

                    auto g1_xyy_xxyzz = contrBuffer.data(g1off + 210 * j + 71);

                    auto g1_xyy_xxzzz = contrBuffer.data(g1off + 210 * j + 72);

                    auto g1_xyy_xyyyy = contrBuffer.data(g1off + 210 * j + 73);

                    auto g1_xyy_xyyyz = contrBuffer.data(g1off + 210 * j + 74);

                    auto g1_xyy_xyyzz = contrBuffer.data(g1off + 210 * j + 75);

                    auto g1_xyy_xyzzz = contrBuffer.data(g1off + 210 * j + 76);

                    auto g1_xyy_xzzzz = contrBuffer.data(g1off + 210 * j + 77);

                    auto g1_xyz_xxxxx = contrBuffer.data(g1off + 210 * j + 84);

                    auto g1_xyz_xxxxy = contrBuffer.data(g1off + 210 * j + 85);

                    auto g1_xyz_xxxxz = contrBuffer.data(g1off + 210 * j + 86);

                    auto g1_xyz_xxxyy = contrBuffer.data(g1off + 210 * j + 87);

                    auto g1_xyz_xxxyz = contrBuffer.data(g1off + 210 * j + 88);

                    auto g1_xyz_xxxzz = contrBuffer.data(g1off + 210 * j + 89);

                    auto g1_xyz_xxyyy = contrBuffer.data(g1off + 210 * j + 90);

                    auto g1_xyz_xxyyz = contrBuffer.data(g1off + 210 * j + 91);

                    auto g1_xyz_xxyzz = contrBuffer.data(g1off + 210 * j + 92);

                    auto g1_xyz_xxzzz = contrBuffer.data(g1off + 210 * j + 93);

                    auto g1_xyz_xyyyy = contrBuffer.data(g1off + 210 * j + 94);

                    auto g1_xyz_xyyyz = contrBuffer.data(g1off + 210 * j + 95);

                    auto g1_xyz_xyyzz = contrBuffer.data(g1off + 210 * j + 96);

                    auto g1_xyz_xyzzz = contrBuffer.data(g1off + 210 * j + 97);

                    auto g1_xyz_xzzzz = contrBuffer.data(g1off + 210 * j + 98);

                    auto g1_xzz_xxxxx = contrBuffer.data(g1off + 210 * j + 105);

                    auto g1_xzz_xxxxy = contrBuffer.data(g1off + 210 * j + 106);

                    auto g1_xzz_xxxxz = contrBuffer.data(g1off + 210 * j + 107);

                    auto g1_xzz_xxxyy = contrBuffer.data(g1off + 210 * j + 108);

                    auto g1_xzz_xxxyz = contrBuffer.data(g1off + 210 * j + 109);

                    auto g1_xzz_xxxzz = contrBuffer.data(g1off + 210 * j + 110);

                    auto g1_xzz_xxyyy = contrBuffer.data(g1off + 210 * j + 111);

                    auto g1_xzz_xxyyz = contrBuffer.data(g1off + 210 * j + 112);

                    auto g1_xzz_xxyzz = contrBuffer.data(g1off + 210 * j + 113);

                    auto g1_xzz_xxzzz = contrBuffer.data(g1off + 210 * j + 114);

                    auto g1_xzz_xyyyy = contrBuffer.data(g1off + 210 * j + 115);

                    auto g1_xzz_xyyyz = contrBuffer.data(g1off + 210 * j + 116);

                    auto g1_xzz_xyyzz = contrBuffer.data(g1off + 210 * j + 117);

                    auto g1_xzz_xyzzz = contrBuffer.data(g1off + 210 * j + 118);

                    auto g1_xzz_xzzzz = contrBuffer.data(g1off + 210 * j + 119);

                    auto g1_yyy_xxxxx = contrBuffer.data(g1off + 210 * j + 126);

                    auto g1_yyy_xxxxy = contrBuffer.data(g1off + 210 * j + 127);

                    auto g1_yyy_xxxxz = contrBuffer.data(g1off + 210 * j + 128);

                    auto g1_yyy_xxxyy = contrBuffer.data(g1off + 210 * j + 129);

                    auto g1_yyy_xxxyz = contrBuffer.data(g1off + 210 * j + 130);

                    auto g1_yyy_xxxzz = contrBuffer.data(g1off + 210 * j + 131);

                    auto g1_yyy_xxyyy = contrBuffer.data(g1off + 210 * j + 132);

                    auto g1_yyy_xxyyz = contrBuffer.data(g1off + 210 * j + 133);

                    auto g1_yyy_xxyzz = contrBuffer.data(g1off + 210 * j + 134);

                    auto g1_yyy_xxzzz = contrBuffer.data(g1off + 210 * j + 135);

                    auto g1_yyy_xyyyy = contrBuffer.data(g1off + 210 * j + 136);

                    auto g1_yyy_xyyyz = contrBuffer.data(g1off + 210 * j + 137);

                    auto g1_yyy_xyyzz = contrBuffer.data(g1off + 210 * j + 138);

                    auto g1_yyy_xyzzz = contrBuffer.data(g1off + 210 * j + 139);

                    auto g1_yyy_xzzzz = contrBuffer.data(g1off + 210 * j + 140);

                    auto g1_yyy_yyyyy = contrBuffer.data(g1off + 210 * j + 141);

                    auto g1_yyy_yyyyz = contrBuffer.data(g1off + 210 * j + 142);

                    auto g1_yyy_yyyzz = contrBuffer.data(g1off + 210 * j + 143);

                    auto g1_yyy_yyzzz = contrBuffer.data(g1off + 210 * j + 144);

                    auto g1_yyy_yzzzz = contrBuffer.data(g1off + 210 * j + 145);

                    auto g1_yyz_xxxxx = contrBuffer.data(g1off + 210 * j + 147);

                    auto g1_yyz_xxxxy = contrBuffer.data(g1off + 210 * j + 148);

                    auto g1_yyz_xxxxz = contrBuffer.data(g1off + 210 * j + 149);

                    auto g1_yyz_xxxyy = contrBuffer.data(g1off + 210 * j + 150);

                    auto g1_yyz_xxxyz = contrBuffer.data(g1off + 210 * j + 151);

                    auto g1_yyz_xxxzz = contrBuffer.data(g1off + 210 * j + 152);

                    auto g1_yyz_xxyyy = contrBuffer.data(g1off + 210 * j + 153);

                    auto g1_yyz_xxyyz = contrBuffer.data(g1off + 210 * j + 154);

                    auto g1_yyz_xxyzz = contrBuffer.data(g1off + 210 * j + 155);

                    auto g1_yyz_xxzzz = contrBuffer.data(g1off + 210 * j + 156);

                    auto g1_yyz_xyyyy = contrBuffer.data(g1off + 210 * j + 157);

                    auto g1_yyz_xyyyz = contrBuffer.data(g1off + 210 * j + 158);

                    auto g1_yyz_xyyzz = contrBuffer.data(g1off + 210 * j + 159);

                    auto g1_yyz_xyzzz = contrBuffer.data(g1off + 210 * j + 160);

                    auto g1_yyz_xzzzz = contrBuffer.data(g1off + 210 * j + 161);

                    auto g1_yyz_yyyyy = contrBuffer.data(g1off + 210 * j + 162);

                    auto g1_yyz_yyyyz = contrBuffer.data(g1off + 210 * j + 163);

                    auto g1_yyz_yyyzz = contrBuffer.data(g1off + 210 * j + 164);

                    auto g1_yyz_yyzzz = contrBuffer.data(g1off + 210 * j + 165);

                    auto g1_yyz_yzzzz = contrBuffer.data(g1off + 210 * j + 166);

                    auto g1_yzz_xxxxx = contrBuffer.data(g1off + 210 * j + 168);

                    auto g1_yzz_xxxxy = contrBuffer.data(g1off + 210 * j + 169);

                    auto g1_yzz_xxxxz = contrBuffer.data(g1off + 210 * j + 170);

                    auto g1_yzz_xxxyy = contrBuffer.data(g1off + 210 * j + 171);

                    auto g1_yzz_xxxyz = contrBuffer.data(g1off + 210 * j + 172);

                    auto g1_yzz_xxxzz = contrBuffer.data(g1off + 210 * j + 173);

                    auto g1_yzz_xxyyy = contrBuffer.data(g1off + 210 * j + 174);

                    auto g1_yzz_xxyyz = contrBuffer.data(g1off + 210 * j + 175);

                    auto g1_yzz_xxyzz = contrBuffer.data(g1off + 210 * j + 176);

                    auto g1_yzz_xxzzz = contrBuffer.data(g1off + 210 * j + 177);

                    auto g1_yzz_xyyyy = contrBuffer.data(g1off + 210 * j + 178);

                    auto g1_yzz_xyyyz = contrBuffer.data(g1off + 210 * j + 179);

                    auto g1_yzz_xyyzz = contrBuffer.data(g1off + 210 * j + 180);

                    auto g1_yzz_xyzzz = contrBuffer.data(g1off + 210 * j + 181);

                    auto g1_yzz_xzzzz = contrBuffer.data(g1off + 210 * j + 182);

                    auto g1_yzz_yyyyy = contrBuffer.data(g1off + 210 * j + 183);

                    auto g1_yzz_yyyyz = contrBuffer.data(g1off + 210 * j + 184);

                    auto g1_yzz_yyyzz = contrBuffer.data(g1off + 210 * j + 185);

                    auto g1_yzz_yyzzz = contrBuffer.data(g1off + 210 * j + 186);

                    auto g1_yzz_yzzzz = contrBuffer.data(g1off + 210 * j + 187);

                    auto g1_zzz_xxxxx = contrBuffer.data(g1off + 210 * j + 189);

                    auto g1_zzz_xxxxy = contrBuffer.data(g1off + 210 * j + 190);

                    auto g1_zzz_xxxxz = contrBuffer.data(g1off + 210 * j + 191);

                    auto g1_zzz_xxxyy = contrBuffer.data(g1off + 210 * j + 192);

                    auto g1_zzz_xxxyz = contrBuffer.data(g1off + 210 * j + 193);

                    auto g1_zzz_xxxzz = contrBuffer.data(g1off + 210 * j + 194);

                    auto g1_zzz_xxyyy = contrBuffer.data(g1off + 210 * j + 195);

                    auto g1_zzz_xxyyz = contrBuffer.data(g1off + 210 * j + 196);

                    auto g1_zzz_xxyzz = contrBuffer.data(g1off + 210 * j + 197);

                    auto g1_zzz_xxzzz = contrBuffer.data(g1off + 210 * j + 198);

                    auto g1_zzz_xyyyy = contrBuffer.data(g1off + 210 * j + 199);

                    auto g1_zzz_xyyyz = contrBuffer.data(g1off + 210 * j + 200);

                    auto g1_zzz_xyyzz = contrBuffer.data(g1off + 210 * j + 201);

                    auto g1_zzz_xyzzz = contrBuffer.data(g1off + 210 * j + 202);

                    auto g1_zzz_xzzzz = contrBuffer.data(g1off + 210 * j + 203);

                    auto g1_zzz_yyyyy = contrBuffer.data(g1off + 210 * j + 204);

                    auto g1_zzz_yyyyz = contrBuffer.data(g1off + 210 * j + 205);

                    auto g1_zzz_yyyzz = contrBuffer.data(g1off + 210 * j + 206);

                    auto g1_zzz_yyzzz = contrBuffer.data(g1off + 210 * j + 207);

                    auto g1_zzz_yzzzz = contrBuffer.data(g1off + 210 * j + 208);

                    auto g1_zzz_zzzzz = contrBuffer.data(g1off + 210 * j + 209);

                    // set up pointers to (SX|g(r,r')|GG)^(m) integrals

                    auto g_xxxx_xxxx = contrBuffer.data(goff + 225 * j);

                    auto g_xxxx_xxxy = contrBuffer.data(goff + 225 * j + 1);

                    auto g_xxxx_xxxz = contrBuffer.data(goff + 225 * j + 2);

                    auto g_xxxx_xxyy = contrBuffer.data(goff + 225 * j + 3);

                    auto g_xxxx_xxyz = contrBuffer.data(goff + 225 * j + 4);

                    auto g_xxxx_xxzz = contrBuffer.data(goff + 225 * j + 5);

                    auto g_xxxx_xyyy = contrBuffer.data(goff + 225 * j + 6);

                    auto g_xxxx_xyyz = contrBuffer.data(goff + 225 * j + 7);

                    auto g_xxxx_xyzz = contrBuffer.data(goff + 225 * j + 8);

                    auto g_xxxx_xzzz = contrBuffer.data(goff + 225 * j + 9);

                    auto g_xxxx_yyyy = contrBuffer.data(goff + 225 * j + 10);

                    auto g_xxxx_yyyz = contrBuffer.data(goff + 225 * j + 11);

                    auto g_xxxx_yyzz = contrBuffer.data(goff + 225 * j + 12);

                    auto g_xxxx_yzzz = contrBuffer.data(goff + 225 * j + 13);

                    auto g_xxxx_zzzz = contrBuffer.data(goff + 225 * j + 14);

                    auto g_xxxy_xxxx = contrBuffer.data(goff + 225 * j + 15);

                    auto g_xxxy_xxxy = contrBuffer.data(goff + 225 * j + 16);

                    auto g_xxxy_xxxz = contrBuffer.data(goff + 225 * j + 17);

                    auto g_xxxy_xxyy = contrBuffer.data(goff + 225 * j + 18);

                    auto g_xxxy_xxyz = contrBuffer.data(goff + 225 * j + 19);

                    auto g_xxxy_xxzz = contrBuffer.data(goff + 225 * j + 20);

                    auto g_xxxy_xyyy = contrBuffer.data(goff + 225 * j + 21);

                    auto g_xxxy_xyyz = contrBuffer.data(goff + 225 * j + 22);

                    auto g_xxxy_xyzz = contrBuffer.data(goff + 225 * j + 23);

                    auto g_xxxy_xzzz = contrBuffer.data(goff + 225 * j + 24);

                    auto g_xxxy_yyyy = contrBuffer.data(goff + 225 * j + 25);

                    auto g_xxxy_yyyz = contrBuffer.data(goff + 225 * j + 26);

                    auto g_xxxy_yyzz = contrBuffer.data(goff + 225 * j + 27);

                    auto g_xxxy_yzzz = contrBuffer.data(goff + 225 * j + 28);

                    auto g_xxxy_zzzz = contrBuffer.data(goff + 225 * j + 29);

                    auto g_xxxz_xxxx = contrBuffer.data(goff + 225 * j + 30);

                    auto g_xxxz_xxxy = contrBuffer.data(goff + 225 * j + 31);

                    auto g_xxxz_xxxz = contrBuffer.data(goff + 225 * j + 32);

                    auto g_xxxz_xxyy = contrBuffer.data(goff + 225 * j + 33);

                    auto g_xxxz_xxyz = contrBuffer.data(goff + 225 * j + 34);

                    auto g_xxxz_xxzz = contrBuffer.data(goff + 225 * j + 35);

                    auto g_xxxz_xyyy = contrBuffer.data(goff + 225 * j + 36);

                    auto g_xxxz_xyyz = contrBuffer.data(goff + 225 * j + 37);

                    auto g_xxxz_xyzz = contrBuffer.data(goff + 225 * j + 38);

                    auto g_xxxz_xzzz = contrBuffer.data(goff + 225 * j + 39);

                    auto g_xxxz_yyyy = contrBuffer.data(goff + 225 * j + 40);

                    auto g_xxxz_yyyz = contrBuffer.data(goff + 225 * j + 41);

                    auto g_xxxz_yyzz = contrBuffer.data(goff + 225 * j + 42);

                    auto g_xxxz_yzzz = contrBuffer.data(goff + 225 * j + 43);

                    auto g_xxxz_zzzz = contrBuffer.data(goff + 225 * j + 44);

                    auto g_xxyy_xxxx = contrBuffer.data(goff + 225 * j + 45);

                    auto g_xxyy_xxxy = contrBuffer.data(goff + 225 * j + 46);

                    auto g_xxyy_xxxz = contrBuffer.data(goff + 225 * j + 47);

                    auto g_xxyy_xxyy = contrBuffer.data(goff + 225 * j + 48);

                    auto g_xxyy_xxyz = contrBuffer.data(goff + 225 * j + 49);

                    auto g_xxyy_xxzz = contrBuffer.data(goff + 225 * j + 50);

                    auto g_xxyy_xyyy = contrBuffer.data(goff + 225 * j + 51);

                    auto g_xxyy_xyyz = contrBuffer.data(goff + 225 * j + 52);

                    auto g_xxyy_xyzz = contrBuffer.data(goff + 225 * j + 53);

                    auto g_xxyy_xzzz = contrBuffer.data(goff + 225 * j + 54);

                    auto g_xxyy_yyyy = contrBuffer.data(goff + 225 * j + 55);

                    auto g_xxyy_yyyz = contrBuffer.data(goff + 225 * j + 56);

                    auto g_xxyy_yyzz = contrBuffer.data(goff + 225 * j + 57);

                    auto g_xxyy_yzzz = contrBuffer.data(goff + 225 * j + 58);

                    auto g_xxyy_zzzz = contrBuffer.data(goff + 225 * j + 59);

                    auto g_xxyz_xxxx = contrBuffer.data(goff + 225 * j + 60);

                    auto g_xxyz_xxxy = contrBuffer.data(goff + 225 * j + 61);

                    auto g_xxyz_xxxz = contrBuffer.data(goff + 225 * j + 62);

                    auto g_xxyz_xxyy = contrBuffer.data(goff + 225 * j + 63);

                    auto g_xxyz_xxyz = contrBuffer.data(goff + 225 * j + 64);

                    auto g_xxyz_xxzz = contrBuffer.data(goff + 225 * j + 65);

                    auto g_xxyz_xyyy = contrBuffer.data(goff + 225 * j + 66);

                    auto g_xxyz_xyyz = contrBuffer.data(goff + 225 * j + 67);

                    auto g_xxyz_xyzz = contrBuffer.data(goff + 225 * j + 68);

                    auto g_xxyz_xzzz = contrBuffer.data(goff + 225 * j + 69);

                    auto g_xxyz_yyyy = contrBuffer.data(goff + 225 * j + 70);

                    auto g_xxyz_yyyz = contrBuffer.data(goff + 225 * j + 71);

                    auto g_xxyz_yyzz = contrBuffer.data(goff + 225 * j + 72);

                    auto g_xxyz_yzzz = contrBuffer.data(goff + 225 * j + 73);

                    auto g_xxyz_zzzz = contrBuffer.data(goff + 225 * j + 74);

                    auto g_xxzz_xxxx = contrBuffer.data(goff + 225 * j + 75);

                    auto g_xxzz_xxxy = contrBuffer.data(goff + 225 * j + 76);

                    auto g_xxzz_xxxz = contrBuffer.data(goff + 225 * j + 77);

                    auto g_xxzz_xxyy = contrBuffer.data(goff + 225 * j + 78);

                    auto g_xxzz_xxyz = contrBuffer.data(goff + 225 * j + 79);

                    auto g_xxzz_xxzz = contrBuffer.data(goff + 225 * j + 80);

                    auto g_xxzz_xyyy = contrBuffer.data(goff + 225 * j + 81);

                    auto g_xxzz_xyyz = contrBuffer.data(goff + 225 * j + 82);

                    auto g_xxzz_xyzz = contrBuffer.data(goff + 225 * j + 83);

                    auto g_xxzz_xzzz = contrBuffer.data(goff + 225 * j + 84);

                    auto g_xxzz_yyyy = contrBuffer.data(goff + 225 * j + 85);

                    auto g_xxzz_yyyz = contrBuffer.data(goff + 225 * j + 86);

                    auto g_xxzz_yyzz = contrBuffer.data(goff + 225 * j + 87);

                    auto g_xxzz_yzzz = contrBuffer.data(goff + 225 * j + 88);

                    auto g_xxzz_zzzz = contrBuffer.data(goff + 225 * j + 89);

                    auto g_xyyy_xxxx = contrBuffer.data(goff + 225 * j + 90);

                    auto g_xyyy_xxxy = contrBuffer.data(goff + 225 * j + 91);

                    auto g_xyyy_xxxz = contrBuffer.data(goff + 225 * j + 92);

                    auto g_xyyy_xxyy = contrBuffer.data(goff + 225 * j + 93);

                    auto g_xyyy_xxyz = contrBuffer.data(goff + 225 * j + 94);

                    auto g_xyyy_xxzz = contrBuffer.data(goff + 225 * j + 95);

                    auto g_xyyy_xyyy = contrBuffer.data(goff + 225 * j + 96);

                    auto g_xyyy_xyyz = contrBuffer.data(goff + 225 * j + 97);

                    auto g_xyyy_xyzz = contrBuffer.data(goff + 225 * j + 98);

                    auto g_xyyy_xzzz = contrBuffer.data(goff + 225 * j + 99);

                    auto g_xyyy_yyyy = contrBuffer.data(goff + 225 * j + 100);

                    auto g_xyyy_yyyz = contrBuffer.data(goff + 225 * j + 101);

                    auto g_xyyy_yyzz = contrBuffer.data(goff + 225 * j + 102);

                    auto g_xyyy_yzzz = contrBuffer.data(goff + 225 * j + 103);

                    auto g_xyyy_zzzz = contrBuffer.data(goff + 225 * j + 104);

                    auto g_xyyz_xxxx = contrBuffer.data(goff + 225 * j + 105);

                    auto g_xyyz_xxxy = contrBuffer.data(goff + 225 * j + 106);

                    auto g_xyyz_xxxz = contrBuffer.data(goff + 225 * j + 107);

                    auto g_xyyz_xxyy = contrBuffer.data(goff + 225 * j + 108);

                    auto g_xyyz_xxyz = contrBuffer.data(goff + 225 * j + 109);

                    auto g_xyyz_xxzz = contrBuffer.data(goff + 225 * j + 110);

                    auto g_xyyz_xyyy = contrBuffer.data(goff + 225 * j + 111);

                    auto g_xyyz_xyyz = contrBuffer.data(goff + 225 * j + 112);

                    auto g_xyyz_xyzz = contrBuffer.data(goff + 225 * j + 113);

                    auto g_xyyz_xzzz = contrBuffer.data(goff + 225 * j + 114);

                    auto g_xyyz_yyyy = contrBuffer.data(goff + 225 * j + 115);

                    auto g_xyyz_yyyz = contrBuffer.data(goff + 225 * j + 116);

                    auto g_xyyz_yyzz = contrBuffer.data(goff + 225 * j + 117);

                    auto g_xyyz_yzzz = contrBuffer.data(goff + 225 * j + 118);

                    auto g_xyyz_zzzz = contrBuffer.data(goff + 225 * j + 119);

                    auto g_xyzz_xxxx = contrBuffer.data(goff + 225 * j + 120);

                    auto g_xyzz_xxxy = contrBuffer.data(goff + 225 * j + 121);

                    auto g_xyzz_xxxz = contrBuffer.data(goff + 225 * j + 122);

                    auto g_xyzz_xxyy = contrBuffer.data(goff + 225 * j + 123);

                    auto g_xyzz_xxyz = contrBuffer.data(goff + 225 * j + 124);

                    auto g_xyzz_xxzz = contrBuffer.data(goff + 225 * j + 125);

                    auto g_xyzz_xyyy = contrBuffer.data(goff + 225 * j + 126);

                    auto g_xyzz_xyyz = contrBuffer.data(goff + 225 * j + 127);

                    auto g_xyzz_xyzz = contrBuffer.data(goff + 225 * j + 128);

                    auto g_xyzz_xzzz = contrBuffer.data(goff + 225 * j + 129);

                    auto g_xyzz_yyyy = contrBuffer.data(goff + 225 * j + 130);

                    auto g_xyzz_yyyz = contrBuffer.data(goff + 225 * j + 131);

                    auto g_xyzz_yyzz = contrBuffer.data(goff + 225 * j + 132);

                    auto g_xyzz_yzzz = contrBuffer.data(goff + 225 * j + 133);

                    auto g_xyzz_zzzz = contrBuffer.data(goff + 225 * j + 134);

                    auto g_xzzz_xxxx = contrBuffer.data(goff + 225 * j + 135);

                    auto g_xzzz_xxxy = contrBuffer.data(goff + 225 * j + 136);

                    auto g_xzzz_xxxz = contrBuffer.data(goff + 225 * j + 137);

                    auto g_xzzz_xxyy = contrBuffer.data(goff + 225 * j + 138);

                    auto g_xzzz_xxyz = contrBuffer.data(goff + 225 * j + 139);

                    auto g_xzzz_xxzz = contrBuffer.data(goff + 225 * j + 140);

                    auto g_xzzz_xyyy = contrBuffer.data(goff + 225 * j + 141);

                    auto g_xzzz_xyyz = contrBuffer.data(goff + 225 * j + 142);

                    auto g_xzzz_xyzz = contrBuffer.data(goff + 225 * j + 143);

                    auto g_xzzz_xzzz = contrBuffer.data(goff + 225 * j + 144);

                    auto g_xzzz_yyyy = contrBuffer.data(goff + 225 * j + 145);

                    auto g_xzzz_yyyz = contrBuffer.data(goff + 225 * j + 146);

                    auto g_xzzz_yyzz = contrBuffer.data(goff + 225 * j + 147);

                    auto g_xzzz_yzzz = contrBuffer.data(goff + 225 * j + 148);

                    auto g_xzzz_zzzz = contrBuffer.data(goff + 225 * j + 149);

                    auto g_yyyy_xxxx = contrBuffer.data(goff + 225 * j + 150);

                    auto g_yyyy_xxxy = contrBuffer.data(goff + 225 * j + 151);

                    auto g_yyyy_xxxz = contrBuffer.data(goff + 225 * j + 152);

                    auto g_yyyy_xxyy = contrBuffer.data(goff + 225 * j + 153);

                    auto g_yyyy_xxyz = contrBuffer.data(goff + 225 * j + 154);

                    auto g_yyyy_xxzz = contrBuffer.data(goff + 225 * j + 155);

                    auto g_yyyy_xyyy = contrBuffer.data(goff + 225 * j + 156);

                    auto g_yyyy_xyyz = contrBuffer.data(goff + 225 * j + 157);

                    auto g_yyyy_xyzz = contrBuffer.data(goff + 225 * j + 158);

                    auto g_yyyy_xzzz = contrBuffer.data(goff + 225 * j + 159);

                    auto g_yyyy_yyyy = contrBuffer.data(goff + 225 * j + 160);

                    auto g_yyyy_yyyz = contrBuffer.data(goff + 225 * j + 161);

                    auto g_yyyy_yyzz = contrBuffer.data(goff + 225 * j + 162);

                    auto g_yyyy_yzzz = contrBuffer.data(goff + 225 * j + 163);

                    auto g_yyyy_zzzz = contrBuffer.data(goff + 225 * j + 164);

                    auto g_yyyz_xxxx = contrBuffer.data(goff + 225 * j + 165);

                    auto g_yyyz_xxxy = contrBuffer.data(goff + 225 * j + 166);

                    auto g_yyyz_xxxz = contrBuffer.data(goff + 225 * j + 167);

                    auto g_yyyz_xxyy = contrBuffer.data(goff + 225 * j + 168);

                    auto g_yyyz_xxyz = contrBuffer.data(goff + 225 * j + 169);

                    auto g_yyyz_xxzz = contrBuffer.data(goff + 225 * j + 170);

                    auto g_yyyz_xyyy = contrBuffer.data(goff + 225 * j + 171);

                    auto g_yyyz_xyyz = contrBuffer.data(goff + 225 * j + 172);

                    auto g_yyyz_xyzz = contrBuffer.data(goff + 225 * j + 173);

                    auto g_yyyz_xzzz = contrBuffer.data(goff + 225 * j + 174);

                    auto g_yyyz_yyyy = contrBuffer.data(goff + 225 * j + 175);

                    auto g_yyyz_yyyz = contrBuffer.data(goff + 225 * j + 176);

                    auto g_yyyz_yyzz = contrBuffer.data(goff + 225 * j + 177);

                    auto g_yyyz_yzzz = contrBuffer.data(goff + 225 * j + 178);

                    auto g_yyyz_zzzz = contrBuffer.data(goff + 225 * j + 179);

                    auto g_yyzz_xxxx = contrBuffer.data(goff + 225 * j + 180);

                    auto g_yyzz_xxxy = contrBuffer.data(goff + 225 * j + 181);

                    auto g_yyzz_xxxz = contrBuffer.data(goff + 225 * j + 182);

                    auto g_yyzz_xxyy = contrBuffer.data(goff + 225 * j + 183);

                    auto g_yyzz_xxyz = contrBuffer.data(goff + 225 * j + 184);

                    auto g_yyzz_xxzz = contrBuffer.data(goff + 225 * j + 185);

                    auto g_yyzz_xyyy = contrBuffer.data(goff + 225 * j + 186);

                    auto g_yyzz_xyyz = contrBuffer.data(goff + 225 * j + 187);

                    auto g_yyzz_xyzz = contrBuffer.data(goff + 225 * j + 188);

                    auto g_yyzz_xzzz = contrBuffer.data(goff + 225 * j + 189);

                    auto g_yyzz_yyyy = contrBuffer.data(goff + 225 * j + 190);

                    auto g_yyzz_yyyz = contrBuffer.data(goff + 225 * j + 191);

                    auto g_yyzz_yyzz = contrBuffer.data(goff + 225 * j + 192);

                    auto g_yyzz_yzzz = contrBuffer.data(goff + 225 * j + 193);

                    auto g_yyzz_zzzz = contrBuffer.data(goff + 225 * j + 194);

                    auto g_yzzz_xxxx = contrBuffer.data(goff + 225 * j + 195);

                    auto g_yzzz_xxxy = contrBuffer.data(goff + 225 * j + 196);

                    auto g_yzzz_xxxz = contrBuffer.data(goff + 225 * j + 197);

                    auto g_yzzz_xxyy = contrBuffer.data(goff + 225 * j + 198);

                    auto g_yzzz_xxyz = contrBuffer.data(goff + 225 * j + 199);

                    auto g_yzzz_xxzz = contrBuffer.data(goff + 225 * j + 200);

                    auto g_yzzz_xyyy = contrBuffer.data(goff + 225 * j + 201);

                    auto g_yzzz_xyyz = contrBuffer.data(goff + 225 * j + 202);

                    auto g_yzzz_xyzz = contrBuffer.data(goff + 225 * j + 203);

                    auto g_yzzz_xzzz = contrBuffer.data(goff + 225 * j + 204);

                    auto g_yzzz_yyyy = contrBuffer.data(goff + 225 * j + 205);

                    auto g_yzzz_yyyz = contrBuffer.data(goff + 225 * j + 206);

                    auto g_yzzz_yyzz = contrBuffer.data(goff + 225 * j + 207);

                    auto g_yzzz_yzzz = contrBuffer.data(goff + 225 * j + 208);

                    auto g_yzzz_zzzz = contrBuffer.data(goff + 225 * j + 209);

                    auto g_zzzz_xxxx = contrBuffer.data(goff + 225 * j + 210);

                    auto g_zzzz_xxxy = contrBuffer.data(goff + 225 * j + 211);

                    auto g_zzzz_xxxz = contrBuffer.data(goff + 225 * j + 212);

                    auto g_zzzz_xxyy = contrBuffer.data(goff + 225 * j + 213);

                    auto g_zzzz_xxyz = contrBuffer.data(goff + 225 * j + 214);

                    auto g_zzzz_xxzz = contrBuffer.data(goff + 225 * j + 215);

                    auto g_zzzz_xyyy = contrBuffer.data(goff + 225 * j + 216);

                    auto g_zzzz_xyyz = contrBuffer.data(goff + 225 * j + 217);

                    auto g_zzzz_xyzz = contrBuffer.data(goff + 225 * j + 218);

                    auto g_zzzz_xzzz = contrBuffer.data(goff + 225 * j + 219);

                    auto g_zzzz_yyyy = contrBuffer.data(goff + 225 * j + 220);

                    auto g_zzzz_yyyz = contrBuffer.data(goff + 225 * j + 221);

                    auto g_zzzz_yyzz = contrBuffer.data(goff + 225 * j + 222);

                    auto g_zzzz_yzzz = contrBuffer.data(goff + 225 * j + 223);

                    auto g_zzzz_zzzz = contrBuffer.data(goff + 225 * j + 224);

                    #pragma omp simd aligned(rcdx, rcdy, rcdz, g2_xxx_xxxx, g2_xxx_xxxy,\
                                             g2_xxx_xxxz, g2_xxx_xxyy, g2_xxx_xxyz,\
                                             g2_xxx_xxzz, g2_xxx_xyyy, g2_xxx_xyyz,\
                                             g2_xxx_xyzz, g2_xxx_xzzz, g2_xxx_yyyy,\
                                             g2_xxx_yyyz, g2_xxx_yyzz, g2_xxx_yzzz,\
                                             g2_xxx_zzzz, g2_xxy_xxxx, g2_xxy_xxxy,\
                                             g2_xxy_xxxz, g2_xxy_xxyy, g2_xxy_xxyz,\
                                             g2_xxy_xxzz, g2_xxy_xyyy, g2_xxy_xyyz,\
                                             g2_xxy_xyzz, g2_xxy_xzzz, g2_xxy_yyyy,\
                                             g2_xxy_yyyz, g2_xxy_yyzz, g2_xxy_yzzz,\
                                             g2_xxy_zzzz, g2_xxz_xxxx, g2_xxz_xxxy,\
                                             g2_xxz_xxxz, g2_xxz_xxyy, g2_xxz_xxyz,\
                                             g2_xxz_xxzz, g2_xxz_xyyy, g2_xxz_xyyz,\
                                             g2_xxz_xyzz, g2_xxz_xzzz, g2_xxz_yyyy,\
                                             g2_xxz_yyyz, g2_xxz_yyzz, g2_xxz_yzzz,\
                                             g2_xxz_zzzz, g2_xyy_xxxx, g2_xyy_xxxy,\
                                             g2_xyy_xxxz, g2_xyy_xxyy, g2_xyy_xxyz,\
                                             g2_xyy_xxzz, g2_xyy_xyyy, g2_xyy_xyyz,\
                                             g2_xyy_xyzz, g2_xyy_xzzz, g2_xyy_yyyy,\
                                             g2_xyy_yyyz, g2_xyy_yyzz, g2_xyy_yzzz,\
                                             g2_xyy_zzzz, g2_xyz_xxxx, g2_xyz_xxxy,\
                                             g2_xyz_xxxz, g2_xyz_xxyy, g2_xyz_xxyz,\
                                             g2_xyz_xxzz, g2_xyz_xyyy, g2_xyz_xyyz,\
                                             g2_xyz_xyzz, g2_xyz_xzzz, g2_xyz_yyyy,\
                                             g2_xyz_yyyz, g2_xyz_yyzz, g2_xyz_yzzz,\
                                             g2_xyz_zzzz, g2_xzz_xxxx, g2_xzz_xxxy,\
                                             g2_xzz_xxxz, g2_xzz_xxyy, g2_xzz_xxyz,\
                                             g2_xzz_xxzz, g2_xzz_xyyy, g2_xzz_xyyz,\
                                             g2_xzz_xyzz, g2_xzz_xzzz, g2_xzz_yyyy,\
                                             g2_xzz_yyyz, g2_xzz_yyzz, g2_xzz_yzzz,\
                                             g2_xzz_zzzz, g2_yyy_xxxx, g2_yyy_xxxy,\
                                             g2_yyy_xxxz, g2_yyy_xxyy, g2_yyy_xxyz,\
                                             g2_yyy_xxzz, g2_yyy_xyyy, g2_yyy_xyyz,\
                                             g2_yyy_xyzz, g2_yyy_xzzz, g2_yyy_yyyy,\
                                             g2_yyy_yyyz, g2_yyy_yyzz, g2_yyy_yzzz,\
                                             g2_yyy_zzzz, g2_yyz_xxxx, g2_yyz_xxxy,\
                                             g2_yyz_xxxz, g2_yyz_xxyy, g2_yyz_xxyz,\
                                             g2_yyz_xxzz, g2_yyz_xyyy, g2_yyz_xyyz,\
                                             g2_yyz_xyzz, g2_yyz_xzzz, g2_yyz_yyyy,\
                                             g2_yyz_yyyz, g2_yyz_yyzz, g2_yyz_yzzz,\
                                             g2_yyz_zzzz, g2_yzz_xxxx, g2_yzz_xxxy,\
                                             g2_yzz_xxxz, g2_yzz_xxyy, g2_yzz_xxyz,\
                                             g2_yzz_xxzz, g2_yzz_xyyy, g2_yzz_xyyz,\
                                             g2_yzz_xyzz, g2_yzz_xzzz, g2_yzz_yyyy,\
                                             g2_yzz_yyyz, g2_yzz_yyzz, g2_yzz_yzzz,\
                                             g2_yzz_zzzz, g2_zzz_xxxx, g2_zzz_xxxy,\
                                             g2_zzz_xxxz, g2_zzz_xxyy, g2_zzz_xxyz,\
                                             g2_zzz_xxzz, g2_zzz_xyyy, g2_zzz_xyyz,\
                                             g2_zzz_xyzz, g2_zzz_xzzz, g2_zzz_yyyy,\
                                             g2_zzz_yyyz, g2_zzz_yyzz, g2_zzz_yzzz,\
                                             g2_zzz_zzzz, g1_xxx_xxxxx, g1_xxx_xxxxy,\
                                             g1_xxx_xxxxz, g1_xxx_xxxyy, g1_xxx_xxxyz,\
                                             g1_xxx_xxxzz, g1_xxx_xxyyy, g1_xxx_xxyyz,\
                                             g1_xxx_xxyzz, g1_xxx_xxzzz, g1_xxx_xyyyy,\
                                             g1_xxx_xyyyz, g1_xxx_xyyzz, g1_xxx_xyzzz,\
                                             g1_xxx_xzzzz, g1_xxy_xxxxx, g1_xxy_xxxxy,\
                                             g1_xxy_xxxxz, g1_xxy_xxxyy, g1_xxy_xxxyz,\
                                             g1_xxy_xxxzz, g1_xxy_xxyyy, g1_xxy_xxyyz,\
                                             g1_xxy_xxyzz, g1_xxy_xxzzz, g1_xxy_xyyyy,\
                                             g1_xxy_xyyyz, g1_xxy_xyyzz, g1_xxy_xyzzz,\
                                             g1_xxy_xzzzz, g1_xxz_xxxxx, g1_xxz_xxxxy,\
                                             g1_xxz_xxxxz, g1_xxz_xxxyy, g1_xxz_xxxyz,\
                                             g1_xxz_xxxzz, g1_xxz_xxyyy, g1_xxz_xxyyz,\
                                             g1_xxz_xxyzz, g1_xxz_xxzzz, g1_xxz_xyyyy,\
                                             g1_xxz_xyyyz, g1_xxz_xyyzz, g1_xxz_xyzzz,\
                                             g1_xxz_xzzzz, g1_xyy_xxxxx, g1_xyy_xxxxy,\
                                             g1_xyy_xxxxz, g1_xyy_xxxyy, g1_xyy_xxxyz,\
                                             g1_xyy_xxxzz, g1_xyy_xxyyy, g1_xyy_xxyyz,\
                                             g1_xyy_xxyzz, g1_xyy_xxzzz, g1_xyy_xyyyy,\
                                             g1_xyy_xyyyz, g1_xyy_xyyzz, g1_xyy_xyzzz,\
                                             g1_xyy_xzzzz, g1_xyz_xxxxx, g1_xyz_xxxxy,\
                                             g1_xyz_xxxxz, g1_xyz_xxxyy, g1_xyz_xxxyz,\
                                             g1_xyz_xxxzz, g1_xyz_xxyyy, g1_xyz_xxyyz,\
                                             g1_xyz_xxyzz, g1_xyz_xxzzz, g1_xyz_xyyyy,\
                                             g1_xyz_xyyyz, g1_xyz_xyyzz, g1_xyz_xyzzz,\
                                             g1_xyz_xzzzz, g1_xzz_xxxxx, g1_xzz_xxxxy,\
                                             g1_xzz_xxxxz, g1_xzz_xxxyy, g1_xzz_xxxyz,\
                                             g1_xzz_xxxzz, g1_xzz_xxyyy, g1_xzz_xxyyz,\
                                             g1_xzz_xxyzz, g1_xzz_xxzzz, g1_xzz_xyyyy,\
                                             g1_xzz_xyyyz, g1_xzz_xyyzz, g1_xzz_xyzzz,\
                                             g1_xzz_xzzzz, g1_yyy_xxxxx, g1_yyy_xxxxy,\
                                             g1_yyy_xxxxz, g1_yyy_xxxyy, g1_yyy_xxxyz,\
                                             g1_yyy_xxxzz, g1_yyy_xxyyy, g1_yyy_xxyyz,\
                                             g1_yyy_xxyzz, g1_yyy_xxzzz, g1_yyy_xyyyy,\
                                             g1_yyy_xyyyz, g1_yyy_xyyzz, g1_yyy_xyzzz,\
                                             g1_yyy_xzzzz, g1_yyy_yyyyy, g1_yyy_yyyyz,\
                                             g1_yyy_yyyzz, g1_yyy_yyzzz, g1_yyy_yzzzz,\
                                             g1_yyz_xxxxx, g1_yyz_xxxxy,\
                                             g1_yyz_xxxxz, g1_yyz_xxxyy, g1_yyz_xxxyz,\
                                             g1_yyz_xxxzz, g1_yyz_xxyyy, g1_yyz_xxyyz,\
                                             g1_yyz_xxyzz, g1_yyz_xxzzz, g1_yyz_xyyyy,\
                                             g1_yyz_xyyyz, g1_yyz_xyyzz, g1_yyz_xyzzz,\
                                             g1_yyz_xzzzz, g1_yyz_yyyyy, g1_yyz_yyyyz,\
                                             g1_yyz_yyyzz, g1_yyz_yyzzz, g1_yyz_yzzzz,\
                                             g1_yzz_xxxxx, g1_yzz_xxxxy,\
                                             g1_yzz_xxxxz, g1_yzz_xxxyy, g1_yzz_xxxyz,\
                                             g1_yzz_xxxzz, g1_yzz_xxyyy, g1_yzz_xxyyz,\
                                             g1_yzz_xxyzz, g1_yzz_xxzzz, g1_yzz_xyyyy,\
                                             g1_yzz_xyyyz, g1_yzz_xyyzz, g1_yzz_xyzzz,\
                                             g1_yzz_xzzzz, g1_yzz_yyyyy, g1_yzz_yyyyz,\
                                             g1_yzz_yyyzz, g1_yzz_yyzzz, g1_yzz_yzzzz,\
                                             g1_zzz_xxxxx, g1_zzz_xxxxy,\
                                             g1_zzz_xxxxz, g1_zzz_xxxyy, g1_zzz_xxxyz,\
                                             g1_zzz_xxxzz, g1_zzz_xxyyy, g1_zzz_xxyyz,\
                                             g1_zzz_xxyzz, g1_zzz_xxzzz, g1_zzz_xyyyy,\
                                             g1_zzz_xyyyz, g1_zzz_xyyzz, g1_zzz_xyzzz,\
                                             g1_zzz_xzzzz, g1_zzz_yyyyy, g1_zzz_yyyyz,\
                                             g1_zzz_yyyzz, g1_zzz_yyzzz, g1_zzz_yzzzz,\
                                             g1_zzz_zzzzz, g_xxxx_xxxx, g_xxxx_xxxy,\
                                             g_xxxx_xxxz, g_xxxx_xxyy, g_xxxx_xxyz,\
                                             g_xxxx_xxzz, g_xxxx_xyyy, g_xxxx_xyyz,\
                                             g_xxxx_xyzz, g_xxxx_xzzz, g_xxxx_yyyy,\
                                             g_xxxx_yyyz, g_xxxx_yyzz, g_xxxx_yzzz,\
                                             g_xxxx_zzzz, g_xxxy_xxxx, g_xxxy_xxxy,\
                                             g_xxxy_xxxz, g_xxxy_xxyy, g_xxxy_xxyz,\
                                             g_xxxy_xxzz, g_xxxy_xyyy, g_xxxy_xyyz,\
                                             g_xxxy_xyzz, g_xxxy_xzzz, g_xxxy_yyyy,\
                                             g_xxxy_yyyz, g_xxxy_yyzz, g_xxxy_yzzz,\
                                             g_xxxy_zzzz, g_xxxz_xxxx, g_xxxz_xxxy,\
                                             g_xxxz_xxxz, g_xxxz_xxyy, g_xxxz_xxyz,\
                                             g_xxxz_xxzz, g_xxxz_xyyy, g_xxxz_xyyz,\
                                             g_xxxz_xyzz, g_xxxz_xzzz, g_xxxz_yyyy,\
                                             g_xxxz_yyyz, g_xxxz_yyzz, g_xxxz_yzzz,\
                                             g_xxxz_zzzz, g_xxyy_xxxx, g_xxyy_xxxy,\
                                             g_xxyy_xxxz, g_xxyy_xxyy, g_xxyy_xxyz,\
                                             g_xxyy_xxzz, g_xxyy_xyyy, g_xxyy_xyyz,\
                                             g_xxyy_xyzz, g_xxyy_xzzz, g_xxyy_yyyy,\
                                             g_xxyy_yyyz, g_xxyy_yyzz, g_xxyy_yzzz,\
                                             g_xxyy_zzzz, g_xxyz_xxxx, g_xxyz_xxxy,\
                                             g_xxyz_xxxz, g_xxyz_xxyy, g_xxyz_xxyz,\
                                             g_xxyz_xxzz, g_xxyz_xyyy, g_xxyz_xyyz,\
                                             g_xxyz_xyzz, g_xxyz_xzzz, g_xxyz_yyyy,\
                                             g_xxyz_yyyz, g_xxyz_yyzz, g_xxyz_yzzz,\
                                             g_xxyz_zzzz, g_xxzz_xxxx, g_xxzz_xxxy,\
                                             g_xxzz_xxxz, g_xxzz_xxyy, g_xxzz_xxyz,\
                                             g_xxzz_xxzz, g_xxzz_xyyy, g_xxzz_xyyz,\
                                             g_xxzz_xyzz, g_xxzz_xzzz, g_xxzz_yyyy,\
                                             g_xxzz_yyyz, g_xxzz_yyzz, g_xxzz_yzzz,\
                                             g_xxzz_zzzz, g_xyyy_xxxx, g_xyyy_xxxy,\
                                             g_xyyy_xxxz, g_xyyy_xxyy, g_xyyy_xxyz,\
                                             g_xyyy_xxzz, g_xyyy_xyyy, g_xyyy_xyyz,\
                                             g_xyyy_xyzz, g_xyyy_xzzz, g_xyyy_yyyy,\
                                             g_xyyy_yyyz, g_xyyy_yyzz, g_xyyy_yzzz,\
                                             g_xyyy_zzzz, g_xyyz_xxxx, g_xyyz_xxxy,\
                                             g_xyyz_xxxz, g_xyyz_xxyy, g_xyyz_xxyz,\
                                             g_xyyz_xxzz, g_xyyz_xyyy, g_xyyz_xyyz,\
                                             g_xyyz_xyzz, g_xyyz_xzzz, g_xyyz_yyyy,\
                                             g_xyyz_yyyz, g_xyyz_yyzz, g_xyyz_yzzz,\
                                             g_xyyz_zzzz, g_xyzz_xxxx, g_xyzz_xxxy,\
                                             g_xyzz_xxxz, g_xyzz_xxyy, g_xyzz_xxyz,\
                                             g_xyzz_xxzz, g_xyzz_xyyy, g_xyzz_xyyz,\
                                             g_xyzz_xyzz, g_xyzz_xzzz, g_xyzz_yyyy,\
                                             g_xyzz_yyyz, g_xyzz_yyzz, g_xyzz_yzzz,\
                                             g_xyzz_zzzz, g_xzzz_xxxx, g_xzzz_xxxy,\
                                             g_xzzz_xxxz, g_xzzz_xxyy, g_xzzz_xxyz,\
                                             g_xzzz_xxzz, g_xzzz_xyyy, g_xzzz_xyyz,\
                                             g_xzzz_xyzz, g_xzzz_xzzz, g_xzzz_yyyy,\
                                             g_xzzz_yyyz, g_xzzz_yyzz, g_xzzz_yzzz,\
                                             g_xzzz_zzzz, g_yyyy_xxxx, g_yyyy_xxxy,\
                                             g_yyyy_xxxz, g_yyyy_xxyy, g_yyyy_xxyz,\
                                             g_yyyy_xxzz, g_yyyy_xyyy, g_yyyy_xyyz,\
                                             g_yyyy_xyzz, g_yyyy_xzzz, g_yyyy_yyyy,\
                                             g_yyyy_yyyz, g_yyyy_yyzz, g_yyyy_yzzz,\
                                             g_yyyy_zzzz, g_yyyz_xxxx, g_yyyz_xxxy,\
                                             g_yyyz_xxxz, g_yyyz_xxyy, g_yyyz_xxyz,\
                                             g_yyyz_xxzz, g_yyyz_xyyy, g_yyyz_xyyz,\
                                             g_yyyz_xyzz, g_yyyz_xzzz, g_yyyz_yyyy,\
                                             g_yyyz_yyyz, g_yyyz_yyzz, g_yyyz_yzzz,\
                                             g_yyyz_zzzz, g_yyzz_xxxx, g_yyzz_xxxy,\
                                             g_yyzz_xxxz, g_yyzz_xxyy, g_yyzz_xxyz,\
                                             g_yyzz_xxzz, g_yyzz_xyyy, g_yyzz_xyyz,\
                                             g_yyzz_xyzz, g_yyzz_xzzz, g_yyzz_yyyy,\
                                             g_yyzz_yyyz, g_yyzz_yyzz, g_yyzz_yzzz,\
                                             g_yyzz_zzzz, g_yzzz_xxxx, g_yzzz_xxxy,\
                                             g_yzzz_xxxz, g_yzzz_xxyy, g_yzzz_xxyz,\
                                             g_yzzz_xxzz, g_yzzz_xyyy, g_yzzz_xyyz,\
                                             g_yzzz_xyzz, g_yzzz_xzzz, g_yzzz_yyyy,\
                                             g_yzzz_yyyz, g_yzzz_yyzz, g_yzzz_yzzz,\
                                             g_yzzz_zzzz, g_zzzz_xxxx, g_zzzz_xxxy,\
                                             g_zzzz_xxxz, g_zzzz_xxyy, g_zzzz_xxyz,\
                                             g_zzzz_xxzz, g_zzzz_xyyy, g_zzzz_xyyz,\
                                             g_zzzz_xyzz, g_zzzz_xzzz, g_zzzz_yyyy,\
                                             g_zzzz_yyyz, g_zzzz_yyzz, g_zzzz_yzzz,\
                                             g_zzzz_zzzz: VLX_ALIGN)
                     for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        double fr = rcdx[k];

                        g_xxxx_xxxx[k] = g1_xxx_xxxxx[k] - fr * g2_xxx_xxxx[k];

                        g_xxxx_xxxy[k] = g1_xxx_xxxxy[k] - fr * g2_xxx_xxxy[k];

                        g_xxxx_xxxz[k] = g1_xxx_xxxxz[k] - fr * g2_xxx_xxxz[k];

                        g_xxxx_xxyy[k] = g1_xxx_xxxyy[k] - fr * g2_xxx_xxyy[k];

                        g_xxxx_xxyz[k] = g1_xxx_xxxyz[k] - fr * g2_xxx_xxyz[k];

                        g_xxxx_xxzz[k] = g1_xxx_xxxzz[k] - fr * g2_xxx_xxzz[k];

                        g_xxxx_xyyy[k] = g1_xxx_xxyyy[k] - fr * g2_xxx_xyyy[k];

                        g_xxxx_xyyz[k] = g1_xxx_xxyyz[k] - fr * g2_xxx_xyyz[k];

                        g_xxxx_xyzz[k] = g1_xxx_xxyzz[k] - fr * g2_xxx_xyzz[k];

                        g_xxxx_xzzz[k] = g1_xxx_xxzzz[k] - fr * g2_xxx_xzzz[k];

                        g_xxxx_yyyy[k] = g1_xxx_xyyyy[k] - fr * g2_xxx_yyyy[k];

                        g_xxxx_yyyz[k] = g1_xxx_xyyyz[k] - fr * g2_xxx_yyyz[k];

                        g_xxxx_yyzz[k] = g1_xxx_xyyzz[k] - fr * g2_xxx_yyzz[k];

                        g_xxxx_yzzz[k] = g1_xxx_xyzzz[k] - fr * g2_xxx_yzzz[k];

                        g_xxxx_zzzz[k] = g1_xxx_xzzzz[k] - fr * g2_xxx_zzzz[k];

                        g_xxxy_xxxx[k] = g1_xxy_xxxxx[k] - fr * g2_xxy_xxxx[k];

                        g_xxxy_xxxy[k] = g1_xxy_xxxxy[k] - fr * g2_xxy_xxxy[k];

                        g_xxxy_xxxz[k] = g1_xxy_xxxxz[k] - fr * g2_xxy_xxxz[k];

                        g_xxxy_xxyy[k] = g1_xxy_xxxyy[k] - fr * g2_xxy_xxyy[k];

                        g_xxxy_xxyz[k] = g1_xxy_xxxyz[k] - fr * g2_xxy_xxyz[k];

                        g_xxxy_xxzz[k] = g1_xxy_xxxzz[k] - fr * g2_xxy_xxzz[k];

                        g_xxxy_xyyy[k] = g1_xxy_xxyyy[k] - fr * g2_xxy_xyyy[k];

                        g_xxxy_xyyz[k] = g1_xxy_xxyyz[k] - fr * g2_xxy_xyyz[k];

                        g_xxxy_xyzz[k] = g1_xxy_xxyzz[k] - fr * g2_xxy_xyzz[k];

                        g_xxxy_xzzz[k] = g1_xxy_xxzzz[k] - fr * g2_xxy_xzzz[k];

                        g_xxxy_yyyy[k] = g1_xxy_xyyyy[k] - fr * g2_xxy_yyyy[k];

                        g_xxxy_yyyz[k] = g1_xxy_xyyyz[k] - fr * g2_xxy_yyyz[k];

                        g_xxxy_yyzz[k] = g1_xxy_xyyzz[k] - fr * g2_xxy_yyzz[k];

                        g_xxxy_yzzz[k] = g1_xxy_xyzzz[k] - fr * g2_xxy_yzzz[k];

                        g_xxxy_zzzz[k] = g1_xxy_xzzzz[k] - fr * g2_xxy_zzzz[k];

                        g_xxxz_xxxx[k] = g1_xxz_xxxxx[k] - fr * g2_xxz_xxxx[k];

                        g_xxxz_xxxy[k] = g1_xxz_xxxxy[k] - fr * g2_xxz_xxxy[k];

                        g_xxxz_xxxz[k] = g1_xxz_xxxxz[k] - fr * g2_xxz_xxxz[k];

                        g_xxxz_xxyy[k] = g1_xxz_xxxyy[k] - fr * g2_xxz_xxyy[k];

                        g_xxxz_xxyz[k] = g1_xxz_xxxyz[k] - fr * g2_xxz_xxyz[k];

                        g_xxxz_xxzz[k] = g1_xxz_xxxzz[k] - fr * g2_xxz_xxzz[k];

                        g_xxxz_xyyy[k] = g1_xxz_xxyyy[k] - fr * g2_xxz_xyyy[k];

                        g_xxxz_xyyz[k] = g1_xxz_xxyyz[k] - fr * g2_xxz_xyyz[k];

                        g_xxxz_xyzz[k] = g1_xxz_xxyzz[k] - fr * g2_xxz_xyzz[k];

                        g_xxxz_xzzz[k] = g1_xxz_xxzzz[k] - fr * g2_xxz_xzzz[k];

                        g_xxxz_yyyy[k] = g1_xxz_xyyyy[k] - fr * g2_xxz_yyyy[k];

                        g_xxxz_yyyz[k] = g1_xxz_xyyyz[k] - fr * g2_xxz_yyyz[k];

                        g_xxxz_yyzz[k] = g1_xxz_xyyzz[k] - fr * g2_xxz_yyzz[k];

                        g_xxxz_yzzz[k] = g1_xxz_xyzzz[k] - fr * g2_xxz_yzzz[k];

                        g_xxxz_zzzz[k] = g1_xxz_xzzzz[k] - fr * g2_xxz_zzzz[k];

                        g_xxyy_xxxx[k] = g1_xyy_xxxxx[k] - fr * g2_xyy_xxxx[k];

                        g_xxyy_xxxy[k] = g1_xyy_xxxxy[k] - fr * g2_xyy_xxxy[k];

                        g_xxyy_xxxz[k] = g1_xyy_xxxxz[k] - fr * g2_xyy_xxxz[k];

                        g_xxyy_xxyy[k] = g1_xyy_xxxyy[k] - fr * g2_xyy_xxyy[k];

                        g_xxyy_xxyz[k] = g1_xyy_xxxyz[k] - fr * g2_xyy_xxyz[k];

                        g_xxyy_xxzz[k] = g1_xyy_xxxzz[k] - fr * g2_xyy_xxzz[k];

                        g_xxyy_xyyy[k] = g1_xyy_xxyyy[k] - fr * g2_xyy_xyyy[k];

                        g_xxyy_xyyz[k] = g1_xyy_xxyyz[k] - fr * g2_xyy_xyyz[k];

                        g_xxyy_xyzz[k] = g1_xyy_xxyzz[k] - fr * g2_xyy_xyzz[k];

                        g_xxyy_xzzz[k] = g1_xyy_xxzzz[k] - fr * g2_xyy_xzzz[k];

                        g_xxyy_yyyy[k] = g1_xyy_xyyyy[k] - fr * g2_xyy_yyyy[k];

                        g_xxyy_yyyz[k] = g1_xyy_xyyyz[k] - fr * g2_xyy_yyyz[k];

                        g_xxyy_yyzz[k] = g1_xyy_xyyzz[k] - fr * g2_xyy_yyzz[k];

                        g_xxyy_yzzz[k] = g1_xyy_xyzzz[k] - fr * g2_xyy_yzzz[k];

                        g_xxyy_zzzz[k] = g1_xyy_xzzzz[k] - fr * g2_xyy_zzzz[k];

                        g_xxyz_xxxx[k] = g1_xyz_xxxxx[k] - fr * g2_xyz_xxxx[k];

                        g_xxyz_xxxy[k] = g1_xyz_xxxxy[k] - fr * g2_xyz_xxxy[k];

                        g_xxyz_xxxz[k] = g1_xyz_xxxxz[k] - fr * g2_xyz_xxxz[k];

                        g_xxyz_xxyy[k] = g1_xyz_xxxyy[k] - fr * g2_xyz_xxyy[k];

                        g_xxyz_xxyz[k] = g1_xyz_xxxyz[k] - fr * g2_xyz_xxyz[k];

                        g_xxyz_xxzz[k] = g1_xyz_xxxzz[k] - fr * g2_xyz_xxzz[k];

                        g_xxyz_xyyy[k] = g1_xyz_xxyyy[k] - fr * g2_xyz_xyyy[k];

                        g_xxyz_xyyz[k] = g1_xyz_xxyyz[k] - fr * g2_xyz_xyyz[k];

                        g_xxyz_xyzz[k] = g1_xyz_xxyzz[k] - fr * g2_xyz_xyzz[k];

                        g_xxyz_xzzz[k] = g1_xyz_xxzzz[k] - fr * g2_xyz_xzzz[k];

                        g_xxyz_yyyy[k] = g1_xyz_xyyyy[k] - fr * g2_xyz_yyyy[k];

                        g_xxyz_yyyz[k] = g1_xyz_xyyyz[k] - fr * g2_xyz_yyyz[k];

                        g_xxyz_yyzz[k] = g1_xyz_xyyzz[k] - fr * g2_xyz_yyzz[k];

                        g_xxyz_yzzz[k] = g1_xyz_xyzzz[k] - fr * g2_xyz_yzzz[k];

                        g_xxyz_zzzz[k] = g1_xyz_xzzzz[k] - fr * g2_xyz_zzzz[k];

                        g_xxzz_xxxx[k] = g1_xzz_xxxxx[k] - fr * g2_xzz_xxxx[k];

                        g_xxzz_xxxy[k] = g1_xzz_xxxxy[k] - fr * g2_xzz_xxxy[k];

                        g_xxzz_xxxz[k] = g1_xzz_xxxxz[k] - fr * g2_xzz_xxxz[k];

                        g_xxzz_xxyy[k] = g1_xzz_xxxyy[k] - fr * g2_xzz_xxyy[k];

                        g_xxzz_xxyz[k] = g1_xzz_xxxyz[k] - fr * g2_xzz_xxyz[k];

                        g_xxzz_xxzz[k] = g1_xzz_xxxzz[k] - fr * g2_xzz_xxzz[k];

                        g_xxzz_xyyy[k] = g1_xzz_xxyyy[k] - fr * g2_xzz_xyyy[k];

                        g_xxzz_xyyz[k] = g1_xzz_xxyyz[k] - fr * g2_xzz_xyyz[k];

                        g_xxzz_xyzz[k] = g1_xzz_xxyzz[k] - fr * g2_xzz_xyzz[k];

                        g_xxzz_xzzz[k] = g1_xzz_xxzzz[k] - fr * g2_xzz_xzzz[k];

                        g_xxzz_yyyy[k] = g1_xzz_xyyyy[k] - fr * g2_xzz_yyyy[k];

                        g_xxzz_yyyz[k] = g1_xzz_xyyyz[k] - fr * g2_xzz_yyyz[k];

                        g_xxzz_yyzz[k] = g1_xzz_xyyzz[k] - fr * g2_xzz_yyzz[k];

                        g_xxzz_yzzz[k] = g1_xzz_xyzzz[k] - fr * g2_xzz_yzzz[k];

                        g_xxzz_zzzz[k] = g1_xzz_xzzzz[k] - fr * g2_xzz_zzzz[k];

                        g_xyyy_xxxx[k] = g1_yyy_xxxxx[k] - fr * g2_yyy_xxxx[k];

                        g_xyyy_xxxy[k] = g1_yyy_xxxxy[k] - fr * g2_yyy_xxxy[k];

                        g_xyyy_xxxz[k] = g1_yyy_xxxxz[k] - fr * g2_yyy_xxxz[k];

                        g_xyyy_xxyy[k] = g1_yyy_xxxyy[k] - fr * g2_yyy_xxyy[k];

                        g_xyyy_xxyz[k] = g1_yyy_xxxyz[k] - fr * g2_yyy_xxyz[k];

                        g_xyyy_xxzz[k] = g1_yyy_xxxzz[k] - fr * g2_yyy_xxzz[k];

                        g_xyyy_xyyy[k] = g1_yyy_xxyyy[k] - fr * g2_yyy_xyyy[k];

                        g_xyyy_xyyz[k] = g1_yyy_xxyyz[k] - fr * g2_yyy_xyyz[k];

                        g_xyyy_xyzz[k] = g1_yyy_xxyzz[k] - fr * g2_yyy_xyzz[k];

                        g_xyyy_xzzz[k] = g1_yyy_xxzzz[k] - fr * g2_yyy_xzzz[k];

                        g_xyyy_yyyy[k] = g1_yyy_xyyyy[k] - fr * g2_yyy_yyyy[k];

                        g_xyyy_yyyz[k] = g1_yyy_xyyyz[k] - fr * g2_yyy_yyyz[k];

                        g_xyyy_yyzz[k] = g1_yyy_xyyzz[k] - fr * g2_yyy_yyzz[k];

                        g_xyyy_yzzz[k] = g1_yyy_xyzzz[k] - fr * g2_yyy_yzzz[k];

                        g_xyyy_zzzz[k] = g1_yyy_xzzzz[k] - fr * g2_yyy_zzzz[k];

                        g_xyyz_xxxx[k] = g1_yyz_xxxxx[k] - fr * g2_yyz_xxxx[k];

                        g_xyyz_xxxy[k] = g1_yyz_xxxxy[k] - fr * g2_yyz_xxxy[k];

                        g_xyyz_xxxz[k] = g1_yyz_xxxxz[k] - fr * g2_yyz_xxxz[k];

                        g_xyyz_xxyy[k] = g1_yyz_xxxyy[k] - fr * g2_yyz_xxyy[k];

                        g_xyyz_xxyz[k] = g1_yyz_xxxyz[k] - fr * g2_yyz_xxyz[k];

                        g_xyyz_xxzz[k] = g1_yyz_xxxzz[k] - fr * g2_yyz_xxzz[k];

                        g_xyyz_xyyy[k] = g1_yyz_xxyyy[k] - fr * g2_yyz_xyyy[k];

                        g_xyyz_xyyz[k] = g1_yyz_xxyyz[k] - fr * g2_yyz_xyyz[k];

                        g_xyyz_xyzz[k] = g1_yyz_xxyzz[k] - fr * g2_yyz_xyzz[k];

                        g_xyyz_xzzz[k] = g1_yyz_xxzzz[k] - fr * g2_yyz_xzzz[k];

                        g_xyyz_yyyy[k] = g1_yyz_xyyyy[k] - fr * g2_yyz_yyyy[k];

                        g_xyyz_yyyz[k] = g1_yyz_xyyyz[k] - fr * g2_yyz_yyyz[k];

                        g_xyyz_yyzz[k] = g1_yyz_xyyzz[k] - fr * g2_yyz_yyzz[k];

                        g_xyyz_yzzz[k] = g1_yyz_xyzzz[k] - fr * g2_yyz_yzzz[k];

                        g_xyyz_zzzz[k] = g1_yyz_xzzzz[k] - fr * g2_yyz_zzzz[k];

                        g_xyzz_xxxx[k] = g1_yzz_xxxxx[k] - fr * g2_yzz_xxxx[k];

                        g_xyzz_xxxy[k] = g1_yzz_xxxxy[k] - fr * g2_yzz_xxxy[k];

                        g_xyzz_xxxz[k] = g1_yzz_xxxxz[k] - fr * g2_yzz_xxxz[k];

                        g_xyzz_xxyy[k] = g1_yzz_xxxyy[k] - fr * g2_yzz_xxyy[k];

                        g_xyzz_xxyz[k] = g1_yzz_xxxyz[k] - fr * g2_yzz_xxyz[k];

                        g_xyzz_xxzz[k] = g1_yzz_xxxzz[k] - fr * g2_yzz_xxzz[k];

                        g_xyzz_xyyy[k] = g1_yzz_xxyyy[k] - fr * g2_yzz_xyyy[k];

                        g_xyzz_xyyz[k] = g1_yzz_xxyyz[k] - fr * g2_yzz_xyyz[k];

                        g_xyzz_xyzz[k] = g1_yzz_xxyzz[k] - fr * g2_yzz_xyzz[k];

                        g_xyzz_xzzz[k] = g1_yzz_xxzzz[k] - fr * g2_yzz_xzzz[k];

                        g_xyzz_yyyy[k] = g1_yzz_xyyyy[k] - fr * g2_yzz_yyyy[k];

                        g_xyzz_yyyz[k] = g1_yzz_xyyyz[k] - fr * g2_yzz_yyyz[k];

                        g_xyzz_yyzz[k] = g1_yzz_xyyzz[k] - fr * g2_yzz_yyzz[k];

                        g_xyzz_yzzz[k] = g1_yzz_xyzzz[k] - fr * g2_yzz_yzzz[k];

                        g_xyzz_zzzz[k] = g1_yzz_xzzzz[k] - fr * g2_yzz_zzzz[k];

                        g_xzzz_xxxx[k] = g1_zzz_xxxxx[k] - fr * g2_zzz_xxxx[k];

                        g_xzzz_xxxy[k] = g1_zzz_xxxxy[k] - fr * g2_zzz_xxxy[k];

                        g_xzzz_xxxz[k] = g1_zzz_xxxxz[k] - fr * g2_zzz_xxxz[k];

                        g_xzzz_xxyy[k] = g1_zzz_xxxyy[k] - fr * g2_zzz_xxyy[k];

                        g_xzzz_xxyz[k] = g1_zzz_xxxyz[k] - fr * g2_zzz_xxyz[k];

                        g_xzzz_xxzz[k] = g1_zzz_xxxzz[k] - fr * g2_zzz_xxzz[k];

                        g_xzzz_xyyy[k] = g1_zzz_xxyyy[k] - fr * g2_zzz_xyyy[k];

                        g_xzzz_xyyz[k] = g1_zzz_xxyyz[k] - fr * g2_zzz_xyyz[k];

                        g_xzzz_xyzz[k] = g1_zzz_xxyzz[k] - fr * g2_zzz_xyzz[k];

                        g_xzzz_xzzz[k] = g1_zzz_xxzzz[k] - fr * g2_zzz_xzzz[k];

                        g_xzzz_yyyy[k] = g1_zzz_xyyyy[k] - fr * g2_zzz_yyyy[k];

                        g_xzzz_yyyz[k] = g1_zzz_xyyyz[k] - fr * g2_zzz_yyyz[k];

                        g_xzzz_yyzz[k] = g1_zzz_xyyzz[k] - fr * g2_zzz_yyzz[k];

                        g_xzzz_yzzz[k] = g1_zzz_xyzzz[k] - fr * g2_zzz_yzzz[k];

                        g_xzzz_zzzz[k] = g1_zzz_xzzzz[k] - fr * g2_zzz_zzzz[k];

                        // leading y component

                        fr = rcdy[k];

                        g_yyyy_xxxx[k] = g1_yyy_xxxxy[k] - fr * g2_yyy_xxxx[k];

                        g_yyyy_xxxy[k] = g1_yyy_xxxyy[k] - fr * g2_yyy_xxxy[k];

                        g_yyyy_xxxz[k] = g1_yyy_xxxyz[k] - fr * g2_yyy_xxxz[k];

                        g_yyyy_xxyy[k] = g1_yyy_xxyyy[k] - fr * g2_yyy_xxyy[k];

                        g_yyyy_xxyz[k] = g1_yyy_xxyyz[k] - fr * g2_yyy_xxyz[k];

                        g_yyyy_xxzz[k] = g1_yyy_xxyzz[k] - fr * g2_yyy_xxzz[k];

                        g_yyyy_xyyy[k] = g1_yyy_xyyyy[k] - fr * g2_yyy_xyyy[k];

                        g_yyyy_xyyz[k] = g1_yyy_xyyyz[k] - fr * g2_yyy_xyyz[k];

                        g_yyyy_xyzz[k] = g1_yyy_xyyzz[k] - fr * g2_yyy_xyzz[k];

                        g_yyyy_xzzz[k] = g1_yyy_xyzzz[k] - fr * g2_yyy_xzzz[k];

                        g_yyyy_yyyy[k] = g1_yyy_yyyyy[k] - fr * g2_yyy_yyyy[k];

                        g_yyyy_yyyz[k] = g1_yyy_yyyyz[k] - fr * g2_yyy_yyyz[k];

                        g_yyyy_yyzz[k] = g1_yyy_yyyzz[k] - fr * g2_yyy_yyzz[k];

                        g_yyyy_yzzz[k] = g1_yyy_yyzzz[k] - fr * g2_yyy_yzzz[k];

                        g_yyyy_zzzz[k] = g1_yyy_yzzzz[k] - fr * g2_yyy_zzzz[k];

                        g_yyyz_xxxx[k] = g1_yyz_xxxxy[k] - fr * g2_yyz_xxxx[k];

                        g_yyyz_xxxy[k] = g1_yyz_xxxyy[k] - fr * g2_yyz_xxxy[k];

                        g_yyyz_xxxz[k] = g1_yyz_xxxyz[k] - fr * g2_yyz_xxxz[k];

                        g_yyyz_xxyy[k] = g1_yyz_xxyyy[k] - fr * g2_yyz_xxyy[k];

                        g_yyyz_xxyz[k] = g1_yyz_xxyyz[k] - fr * g2_yyz_xxyz[k];

                        g_yyyz_xxzz[k] = g1_yyz_xxyzz[k] - fr * g2_yyz_xxzz[k];

                        g_yyyz_xyyy[k] = g1_yyz_xyyyy[k] - fr * g2_yyz_xyyy[k];

                        g_yyyz_xyyz[k] = g1_yyz_xyyyz[k] - fr * g2_yyz_xyyz[k];

                        g_yyyz_xyzz[k] = g1_yyz_xyyzz[k] - fr * g2_yyz_xyzz[k];

                        g_yyyz_xzzz[k] = g1_yyz_xyzzz[k] - fr * g2_yyz_xzzz[k];

                        g_yyyz_yyyy[k] = g1_yyz_yyyyy[k] - fr * g2_yyz_yyyy[k];

                        g_yyyz_yyyz[k] = g1_yyz_yyyyz[k] - fr * g2_yyz_yyyz[k];

                        g_yyyz_yyzz[k] = g1_yyz_yyyzz[k] - fr * g2_yyz_yyzz[k];

                        g_yyyz_yzzz[k] = g1_yyz_yyzzz[k] - fr * g2_yyz_yzzz[k];

                        g_yyyz_zzzz[k] = g1_yyz_yzzzz[k] - fr * g2_yyz_zzzz[k];

                        g_yyzz_xxxx[k] = g1_yzz_xxxxy[k] - fr * g2_yzz_xxxx[k];

                        g_yyzz_xxxy[k] = g1_yzz_xxxyy[k] - fr * g2_yzz_xxxy[k];

                        g_yyzz_xxxz[k] = g1_yzz_xxxyz[k] - fr * g2_yzz_xxxz[k];

                        g_yyzz_xxyy[k] = g1_yzz_xxyyy[k] - fr * g2_yzz_xxyy[k];

                        g_yyzz_xxyz[k] = g1_yzz_xxyyz[k] - fr * g2_yzz_xxyz[k];

                        g_yyzz_xxzz[k] = g1_yzz_xxyzz[k] - fr * g2_yzz_xxzz[k];

                        g_yyzz_xyyy[k] = g1_yzz_xyyyy[k] - fr * g2_yzz_xyyy[k];

                        g_yyzz_xyyz[k] = g1_yzz_xyyyz[k] - fr * g2_yzz_xyyz[k];

                        g_yyzz_xyzz[k] = g1_yzz_xyyzz[k] - fr * g2_yzz_xyzz[k];

                        g_yyzz_xzzz[k] = g1_yzz_xyzzz[k] - fr * g2_yzz_xzzz[k];

                        g_yyzz_yyyy[k] = g1_yzz_yyyyy[k] - fr * g2_yzz_yyyy[k];

                        g_yyzz_yyyz[k] = g1_yzz_yyyyz[k] - fr * g2_yzz_yyyz[k];

                        g_yyzz_yyzz[k] = g1_yzz_yyyzz[k] - fr * g2_yzz_yyzz[k];

                        g_yyzz_yzzz[k] = g1_yzz_yyzzz[k] - fr * g2_yzz_yzzz[k];

                        g_yyzz_zzzz[k] = g1_yzz_yzzzz[k] - fr * g2_yzz_zzzz[k];

                        g_yzzz_xxxx[k] = g1_zzz_xxxxy[k] - fr * g2_zzz_xxxx[k];

                        g_yzzz_xxxy[k] = g1_zzz_xxxyy[k] - fr * g2_zzz_xxxy[k];

                        g_yzzz_xxxz[k] = g1_zzz_xxxyz[k] - fr * g2_zzz_xxxz[k];

                        g_yzzz_xxyy[k] = g1_zzz_xxyyy[k] - fr * g2_zzz_xxyy[k];

                        g_yzzz_xxyz[k] = g1_zzz_xxyyz[k] - fr * g2_zzz_xxyz[k];

                        g_yzzz_xxzz[k] = g1_zzz_xxyzz[k] - fr * g2_zzz_xxzz[k];

                        g_yzzz_xyyy[k] = g1_zzz_xyyyy[k] - fr * g2_zzz_xyyy[k];

                        g_yzzz_xyyz[k] = g1_zzz_xyyyz[k] - fr * g2_zzz_xyyz[k];

                        g_yzzz_xyzz[k] = g1_zzz_xyyzz[k] - fr * g2_zzz_xyzz[k];

                        g_yzzz_xzzz[k] = g1_zzz_xyzzz[k] - fr * g2_zzz_xzzz[k];

                        g_yzzz_yyyy[k] = g1_zzz_yyyyy[k] - fr * g2_zzz_yyyy[k];

                        g_yzzz_yyyz[k] = g1_zzz_yyyyz[k] - fr * g2_zzz_yyyz[k];

                        g_yzzz_yyzz[k] = g1_zzz_yyyzz[k] - fr * g2_zzz_yyzz[k];

                        g_yzzz_yzzz[k] = g1_zzz_yyzzz[k] - fr * g2_zzz_yzzz[k];

                        g_yzzz_zzzz[k] = g1_zzz_yzzzz[k] - fr * g2_zzz_zzzz[k];

                        // leading z component

                        fr = rcdz[k];

                        g_zzzz_xxxx[k] = g1_zzz_xxxxz[k] - fr * g2_zzz_xxxx[k];

                        g_zzzz_xxxy[k] = g1_zzz_xxxyz[k] - fr * g2_zzz_xxxy[k];

                        g_zzzz_xxxz[k] = g1_zzz_xxxzz[k] - fr * g2_zzz_xxxz[k];

                        g_zzzz_xxyy[k] = g1_zzz_xxyyz[k] - fr * g2_zzz_xxyy[k];

                        g_zzzz_xxyz[k] = g1_zzz_xxyzz[k] - fr * g2_zzz_xxyz[k];

                        g_zzzz_xxzz[k] = g1_zzz_xxzzz[k] - fr * g2_zzz_xxzz[k];

                        g_zzzz_xyyy[k] = g1_zzz_xyyyz[k] - fr * g2_zzz_xyyy[k];

                        g_zzzz_xyyz[k] = g1_zzz_xyyzz[k] - fr * g2_zzz_xyyz[k];

                        g_zzzz_xyzz[k] = g1_zzz_xyzzz[k] - fr * g2_zzz_xyzz[k];

                        g_zzzz_xzzz[k] = g1_zzz_xzzzz[k] - fr * g2_zzz_xzzz[k];

                        g_zzzz_yyyy[k] = g1_zzz_yyyyz[k] - fr * g2_zzz_yyyy[k];

                        g_zzzz_yyyz[k] = g1_zzz_yyyzz[k] - fr * g2_zzz_yyyz[k];

                        g_zzzz_yyzz[k] = g1_zzz_yyzzz[k] - fr * g2_zzz_yyzz[k];

                        g_zzzz_yzzz[k] = g1_zzz_yzzzz[k] - fr * g2_zzz_yzzz[k];

                        g_zzzz_zzzz[k] = g1_zzz_zzzzz[k] - fr * g2_zzz_zzzz[k];
                    }
                }
            }
        }
    }
    
} // kethrrfunc namespace
