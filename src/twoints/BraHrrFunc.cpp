//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "BraHrrFunc.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

namespace brahrrfunc { // brahrrfunc namespace

    void
    compElectronRepulsionForPPXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side
        
        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
        
        if (isBraEqualKet) kdim  = iContrPair + 1;
        
        // set up distances R(AB) = A - B
        
        auto abx = (abDistances.data(0))[iContrPair];
        
        auto aby = (abDistances.data(1))[iContrPair];
        
        auto abz = (abDistances.data(2))[iContrPair];
        
        // loop over set of data vectors
        
        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 1))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (11|XX)\n");
                
                // determine angular momentum of ket side
                
                auto cang = recPattern[i].third();
                
                auto dang = recPattern[i].fourth();
                
                auto kcomp = angmom::to_SphericalComponents(cang, dang);
                
                // get position of integrals in integrals buffer
                
                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 1, cang, dang});
                
                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 2, cang, dang});
                
                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 1, cang, dang});
                
                // compute contracted integrals
                
                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SP|g(r,r')|XX) integrals
                    
                    auto g2_0_x = braBuffer.data(g2off + j);
                    
                    auto g2_0_y = braBuffer.data(g2off + kcomp + j);
                    
                    auto g2_0_z = braBuffer.data(g2off + 2 * kcomp + j);
                    
                    // set up pointers to (SD|g(r,r')|XX) integrals
                    
                    auto g1_0_xx = braBuffer.data(g1off + j);
                    
                    auto g1_0_xy = braBuffer.data(g1off + kcomp + j);
                    
                    auto g1_0_xz = braBuffer.data(g1off + 2 * kcomp + j);
                    
                    auto g1_0_yy = braBuffer.data(g1off + 3 * kcomp + j);
                    
                    auto g1_0_yz = braBuffer.data(g1off + 4 * kcomp + j);
                    
                    auto g1_0_zz = braBuffer.data(g1off + 5 * kcomp + j);
                    
                    // set up pointers to (PP|g(r,r')|XX) integrals
                    
                    auto g_x_x = braBuffer.data(goff + j);
                    
                    auto g_x_y = braBuffer.data(goff + kcomp + j);
                    
                    auto g_x_z = braBuffer.data(goff + 2 * kcomp + j);
                    
                    auto g_y_x = braBuffer.data(goff + 3 * kcomp + j);
                    
                    auto g_y_y = braBuffer.data(goff + 4 * kcomp + j);
                    
                    auto g_y_z = braBuffer.data(goff + 5 * kcomp + j);
                    
                    auto g_z_x = braBuffer.data(goff + 6 * kcomp + j);
                    
                    auto g_z_y = braBuffer.data(goff + 7 * kcomp + j);
                    
                    auto g_z_z = braBuffer.data(goff + 8 * kcomp + j);
                    
                    #pragma omp simd aligned(g2_0_x, g2_0_y, g2_0_z, g1_0_xx,\
                                             g1_0_xy, g1_0_xz, g1_0_yy, g1_0_yz,\
                                             g1_0_zz, g_x_x, g_x_y, g_x_z, g_y_x,\
                                             g_y_y, g_y_z, g_z_x, g_z_y,\
                                             g_z_z: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component
                        
                        g_x_x[k] = g1_0_xx[k] - abx * g2_0_x[k];
                        
                        g_x_y[k] = g1_0_xy[k] - abx * g2_0_y[k];
                        
                        g_x_z[k] = g1_0_xz[k] - abx * g2_0_z[k];
                        
                        // leading y component
                        
                        g_y_x[k] = g1_0_xy[k] - aby * g2_0_x[k];
                        
                        g_y_y[k] = g1_0_yy[k] - aby * g2_0_y[k];
                        
                        g_y_z[k] = g1_0_yz[k] - aby * g2_0_z[k];
                        
                        // leading z component
                        
                        g_z_x[k] = g1_0_xz[k] - abz * g2_0_x[k];
                        
                        g_z_y[k] = g1_0_yz[k] - abz * g2_0_y[k];
                        
                        g_z_z[k] = g1_0_zz[k] - abz * g2_0_z[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPDXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 2))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (12|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 2, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 3, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 2, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SD|g(r,r')|XX) integrals

                    auto g2_0_xx = braBuffer.data(g2off + j);

                    auto g2_0_xy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_yy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_yz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_zz = braBuffer.data(g2off + 5 * kcomp + j);

                    // set up pointers to (SF|g(r,r')|XX) integrals

                    auto g1_0_xxx = braBuffer.data(g1off + j);

                    auto g1_0_xxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_yyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_yyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_yzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_zzz = braBuffer.data(g1off + 9 * kcomp + j);

                    // set up pointers to (PD|g(r,r')|XX) integrals

                    auto g_x_xx = braBuffer.data(goff + j);

                    auto g_x_xy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_yy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_yz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_zz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_y_xx = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_y_xy = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_y_xz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_y_yy = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_y_yz = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_y_zz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_z_xx = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_z_xy = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_z_xz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_z_yy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_z_yz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_z_zz = braBuffer.data(goff + 17 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xx, g2_0_xy, g2_0_xz, g2_0_yy,\
                                             g2_0_yz, g2_0_zz, g1_0_xxx, g1_0_xxy,\
                                             g1_0_xxz, g1_0_xyy, g1_0_xyz, g1_0_xzz,\
                                             g1_0_yyy, g1_0_yyz, g1_0_yzz, g1_0_zzz,\
                                             g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz,\
                                             g_x_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy,\
                                             g_y_yz, g_y_zz, g_z_xx, g_z_xy, g_z_xz,\
                                             g_z_yy, g_z_yz, g_z_zz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xx[k] = g1_0_xxx[k] - abx * g2_0_xx[k];

                        g_x_xy[k] = g1_0_xxy[k] - abx * g2_0_xy[k];

                        g_x_xz[k] = g1_0_xxz[k] - abx * g2_0_xz[k];

                        g_x_yy[k] = g1_0_xyy[k] - abx * g2_0_yy[k];

                        g_x_yz[k] = g1_0_xyz[k] - abx * g2_0_yz[k];

                        g_x_zz[k] = g1_0_xzz[k] - abx * g2_0_zz[k];

                        // leading y component

                        g_y_xx[k] = g1_0_xxy[k] - aby * g2_0_xx[k];

                        g_y_xy[k] = g1_0_xyy[k] - aby * g2_0_xy[k];

                        g_y_xz[k] = g1_0_xyz[k] - aby * g2_0_xz[k];

                        g_y_yy[k] = g1_0_yyy[k] - aby * g2_0_yy[k];

                        g_y_yz[k] = g1_0_yyz[k] - aby * g2_0_yz[k];

                        g_y_zz[k] = g1_0_yzz[k] - aby * g2_0_zz[k];

                        // leading z component

                        g_z_xx[k] = g1_0_xxz[k] - abz * g2_0_xx[k];

                        g_z_xy[k] = g1_0_xyz[k] - abz * g2_0_xy[k];

                        g_z_xz[k] = g1_0_xzz[k] - abz * g2_0_xz[k];

                        g_z_yy[k] = g1_0_yyz[k] - abz * g2_0_yy[k];

                        g_z_yz[k] = g1_0_yzz[k] - abz * g2_0_yz[k];

                        g_z_zz[k] = g1_0_zzz[k] - abz * g2_0_zz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPFXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 3))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (13|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 3, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 4, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 3, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SF|g(r,r')|XX) integrals

                    auto g2_0_xxx = braBuffer.data(g2off + j);

                    auto g2_0_xxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_xyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_xyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_xzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_0_yyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_0_yyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_0_yzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_0_zzz = braBuffer.data(g2off + 9 * kcomp + j);

                    // set up pointers to (SG|g(r,r')|XX) integrals

                    auto g1_0_xxxx = braBuffer.data(g1off + j);

                    auto g1_0_xxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_xyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_xyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_xyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_xzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_0_yyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_0_yyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_0_yyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_0_yzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_0_zzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    // set up pointers to (PF|g(r,r')|XX) integrals

                    auto g_x_xxx = braBuffer.data(goff + j);

                    auto g_x_xxy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_xyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_xyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_xzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_x_yyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_x_yyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_x_yzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_x_zzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_y_xxx = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_y_xxy = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_y_xxz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_y_xyy = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_y_xyz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_y_xzz = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_y_yyy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_y_yyz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_y_yzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_y_zzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_z_xxx = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_z_xxy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_z_xxz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_z_xyy = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_z_xyz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_z_xzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_z_yyy = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_z_yyz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_z_yzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_z_zzz = braBuffer.data(goff + 29 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xxx, g2_0_xxy, g2_0_xxz, g2_0_xyy,\
                                             g2_0_xyz, g2_0_xzz, g2_0_yyy, g2_0_yyz,\
                                             g2_0_yzz, g2_0_zzz, g1_0_xxxx, g1_0_xxxy,\
                                             g1_0_xxxz, g1_0_xxyy, g1_0_xxyz, g1_0_xxzz,\
                                             g1_0_xyyy, g1_0_xyyz, g1_0_xyzz, g1_0_xzzz,\
                                             g1_0_yyyy, g1_0_yyyz, g1_0_yyzz, g1_0_yzzz,\
                                             g1_0_zzzz, g_x_xxx, g_x_xxy, g_x_xxz,\
                                             g_x_xyy, g_x_xyz, g_x_xzz, g_x_yyy,\
                                             g_x_yyz, g_x_yzz, g_x_zzz, g_y_xxx,\
                                             g_y_xxy, g_y_xxz, g_y_xyy, g_y_xyz,\
                                             g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz,\
                                             g_y_zzz, g_z_xxx, g_z_xxy, g_z_xxz,\
                                             g_z_xyy, g_z_xyz, g_z_xzz, g_z_yyy,\
                                             g_z_yyz, g_z_yzz, g_z_zzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xxx[k] = g1_0_xxxx[k] - abx * g2_0_xxx[k];

                        g_x_xxy[k] = g1_0_xxxy[k] - abx * g2_0_xxy[k];

                        g_x_xxz[k] = g1_0_xxxz[k] - abx * g2_0_xxz[k];

                        g_x_xyy[k] = g1_0_xxyy[k] - abx * g2_0_xyy[k];

                        g_x_xyz[k] = g1_0_xxyz[k] - abx * g2_0_xyz[k];

                        g_x_xzz[k] = g1_0_xxzz[k] - abx * g2_0_xzz[k];

                        g_x_yyy[k] = g1_0_xyyy[k] - abx * g2_0_yyy[k];

                        g_x_yyz[k] = g1_0_xyyz[k] - abx * g2_0_yyz[k];

                        g_x_yzz[k] = g1_0_xyzz[k] - abx * g2_0_yzz[k];

                        g_x_zzz[k] = g1_0_xzzz[k] - abx * g2_0_zzz[k];

                        // leading y component

                        g_y_xxx[k] = g1_0_xxxy[k] - aby * g2_0_xxx[k];

                        g_y_xxy[k] = g1_0_xxyy[k] - aby * g2_0_xxy[k];

                        g_y_xxz[k] = g1_0_xxyz[k] - aby * g2_0_xxz[k];

                        g_y_xyy[k] = g1_0_xyyy[k] - aby * g2_0_xyy[k];

                        g_y_xyz[k] = g1_0_xyyz[k] - aby * g2_0_xyz[k];

                        g_y_xzz[k] = g1_0_xyzz[k] - aby * g2_0_xzz[k];

                        g_y_yyy[k] = g1_0_yyyy[k] - aby * g2_0_yyy[k];

                        g_y_yyz[k] = g1_0_yyyz[k] - aby * g2_0_yyz[k];

                        g_y_yzz[k] = g1_0_yyzz[k] - aby * g2_0_yzz[k];

                        g_y_zzz[k] = g1_0_yzzz[k] - aby * g2_0_zzz[k];

                        // leading z component

                        g_z_xxx[k] = g1_0_xxxz[k] - abz * g2_0_xxx[k];

                        g_z_xxy[k] = g1_0_xxyz[k] - abz * g2_0_xxy[k];

                        g_z_xxz[k] = g1_0_xxzz[k] - abz * g2_0_xxz[k];

                        g_z_xyy[k] = g1_0_xyyz[k] - abz * g2_0_xyy[k];

                        g_z_xyz[k] = g1_0_xyzz[k] - abz * g2_0_xyz[k];

                        g_z_xzz[k] = g1_0_xzzz[k] - abz * g2_0_xzz[k];

                        g_z_yyy[k] = g1_0_yyyz[k] - abz * g2_0_yyy[k];

                        g_z_yyz[k] = g1_0_yyzz[k] - abz * g2_0_yyz[k];

                        g_z_yzz[k] = g1_0_yzzz[k] - abz * g2_0_yzz[k];

                        g_z_zzz[k] = g1_0_zzzz[k] - abz * g2_0_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPGXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 4))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (14|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 4, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 5, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 4, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SG|g(r,r')|XX) integrals

                    auto g2_0_xxxx = braBuffer.data(g2off + j);

                    auto g2_0_xxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_xxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_xxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_xxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_0_xyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_0_xyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_0_xyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_0_xzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_0_yyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_0_yyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_0_yyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_0_yzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_0_zzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    // set up pointers to (SH|g(r,r')|XX) integrals

                    auto g1_0_xxxxx = braBuffer.data(g1off + j);

                    auto g1_0_xxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_xxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_xxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_xxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_xxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_0_xyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_0_xyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_0_xyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_0_xyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_0_xzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_0_yyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_0_yyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_0_yyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_0_yyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_0_yzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_0_zzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    // set up pointers to (PG|g(r,r')|XX) integrals

                    auto g_x_xxxx = braBuffer.data(goff + j);

                    auto g_x_xxxy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_xxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_xxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_xxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_x_xyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_x_xyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_x_xyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_x_xzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_x_yyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_x_yyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_x_yyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_x_yzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_x_zzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_y_xxxx = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_y_xxxy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_y_xxxz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_y_xxyy = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_y_xxyz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_y_xxzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_y_xyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_y_xyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_y_xyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_y_xzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_y_yyyy = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_y_yyyz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_y_yyzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_y_yzzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_y_zzzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_z_xxxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_z_xxxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_z_xxxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_z_xxyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_z_xxyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_z_xxzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_z_xyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_z_xyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_z_xyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_z_xzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_z_yyyy = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_z_yyyz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_z_yyzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_z_yzzz = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_z_zzzz = braBuffer.data(goff + 44 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xxxx, g2_0_xxxy, g2_0_xxxz, g2_0_xxyy,\
                                             g2_0_xxyz, g2_0_xxzz, g2_0_xyyy, g2_0_xyyz,\
                                             g2_0_xyzz, g2_0_xzzz, g2_0_yyyy, g2_0_yyyz,\
                                             g2_0_yyzz, g2_0_yzzz, g2_0_zzzz, g1_0_xxxxx,\
                                             g1_0_xxxxy, g1_0_xxxxz, g1_0_xxxyy,\
                                             g1_0_xxxyz, g1_0_xxxzz, g1_0_xxyyy,\
                                             g1_0_xxyyz, g1_0_xxyzz, g1_0_xxzzz,\
                                             g1_0_xyyyy, g1_0_xyyyz, g1_0_xyyzz,\
                                             g1_0_xyzzz, g1_0_xzzzz, g1_0_yyyyy,\
                                             g1_0_yyyyz, g1_0_yyyzz, g1_0_yyzzz,\
                                             g1_0_yzzzz, g1_0_zzzzz, g_x_xxxx,\
                                             g_x_xxxy, g_x_xxxz, g_x_xxyy, g_x_xxyz,\
                                             g_x_xxzz, g_x_xyyy, g_x_xyyz, g_x_xyzz,\
                                             g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz,\
                                             g_x_yzzz, g_x_zzzz, g_y_xxxx, g_y_xxxy,\
                                             g_y_xxxz, g_y_xxyy, g_y_xxyz, g_y_xxzz,\
                                             g_y_xyyy, g_y_xyyz, g_y_xyzz, g_y_xzzz,\
                                             g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz,\
                                             g_y_zzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz,\
                                             g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy,\
                                             g_z_xyyz, g_z_xyzz, g_z_xzzz, g_z_yyyy,\
                                             g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xxxx[k] = g1_0_xxxxx[k] - abx * g2_0_xxxx[k];

                        g_x_xxxy[k] = g1_0_xxxxy[k] - abx * g2_0_xxxy[k];

                        g_x_xxxz[k] = g1_0_xxxxz[k] - abx * g2_0_xxxz[k];

                        g_x_xxyy[k] = g1_0_xxxyy[k] - abx * g2_0_xxyy[k];

                        g_x_xxyz[k] = g1_0_xxxyz[k] - abx * g2_0_xxyz[k];

                        g_x_xxzz[k] = g1_0_xxxzz[k] - abx * g2_0_xxzz[k];

                        g_x_xyyy[k] = g1_0_xxyyy[k] - abx * g2_0_xyyy[k];

                        g_x_xyyz[k] = g1_0_xxyyz[k] - abx * g2_0_xyyz[k];

                        g_x_xyzz[k] = g1_0_xxyzz[k] - abx * g2_0_xyzz[k];

                        g_x_xzzz[k] = g1_0_xxzzz[k] - abx * g2_0_xzzz[k];

                        g_x_yyyy[k] = g1_0_xyyyy[k] - abx * g2_0_yyyy[k];

                        g_x_yyyz[k] = g1_0_xyyyz[k] - abx * g2_0_yyyz[k];

                        g_x_yyzz[k] = g1_0_xyyzz[k] - abx * g2_0_yyzz[k];

                        g_x_yzzz[k] = g1_0_xyzzz[k] - abx * g2_0_yzzz[k];

                        g_x_zzzz[k] = g1_0_xzzzz[k] - abx * g2_0_zzzz[k];

                        // leading y component

                        g_y_xxxx[k] = g1_0_xxxxy[k] - aby * g2_0_xxxx[k];

                        g_y_xxxy[k] = g1_0_xxxyy[k] - aby * g2_0_xxxy[k];

                        g_y_xxxz[k] = g1_0_xxxyz[k] - aby * g2_0_xxxz[k];

                        g_y_xxyy[k] = g1_0_xxyyy[k] - aby * g2_0_xxyy[k];

                        g_y_xxyz[k] = g1_0_xxyyz[k] - aby * g2_0_xxyz[k];

                        g_y_xxzz[k] = g1_0_xxyzz[k] - aby * g2_0_xxzz[k];

                        g_y_xyyy[k] = g1_0_xyyyy[k] - aby * g2_0_xyyy[k];

                        g_y_xyyz[k] = g1_0_xyyyz[k] - aby * g2_0_xyyz[k];

                        g_y_xyzz[k] = g1_0_xyyzz[k] - aby * g2_0_xyzz[k];

                        g_y_xzzz[k] = g1_0_xyzzz[k] - aby * g2_0_xzzz[k];

                        g_y_yyyy[k] = g1_0_yyyyy[k] - aby * g2_0_yyyy[k];

                        g_y_yyyz[k] = g1_0_yyyyz[k] - aby * g2_0_yyyz[k];

                        g_y_yyzz[k] = g1_0_yyyzz[k] - aby * g2_0_yyzz[k];

                        g_y_yzzz[k] = g1_0_yyzzz[k] - aby * g2_0_yzzz[k];

                        g_y_zzzz[k] = g1_0_yzzzz[k] - aby * g2_0_zzzz[k];

                        // leading z component

                        g_z_xxxx[k] = g1_0_xxxxz[k] - abz * g2_0_xxxx[k];

                        g_z_xxxy[k] = g1_0_xxxyz[k] - abz * g2_0_xxxy[k];

                        g_z_xxxz[k] = g1_0_xxxzz[k] - abz * g2_0_xxxz[k];

                        g_z_xxyy[k] = g1_0_xxyyz[k] - abz * g2_0_xxyy[k];

                        g_z_xxyz[k] = g1_0_xxyzz[k] - abz * g2_0_xxyz[k];

                        g_z_xxzz[k] = g1_0_xxzzz[k] - abz * g2_0_xxzz[k];

                        g_z_xyyy[k] = g1_0_xyyyz[k] - abz * g2_0_xyyy[k];

                        g_z_xyyz[k] = g1_0_xyyzz[k] - abz * g2_0_xyyz[k];

                        g_z_xyzz[k] = g1_0_xyzzz[k] - abz * g2_0_xyzz[k];

                        g_z_xzzz[k] = g1_0_xzzzz[k] - abz * g2_0_xzzz[k];

                        g_z_yyyy[k] = g1_0_yyyyz[k] - abz * g2_0_yyyy[k];

                        g_z_yyyz[k] = g1_0_yyyzz[k] - abz * g2_0_yyyz[k];

                        g_z_yyzz[k] = g1_0_yyzzz[k] - abz * g2_0_yyzz[k];

                        g_z_yzzz[k] = g1_0_yzzzz[k] - abz * g2_0_yzzz[k];

                        g_z_zzzz[k] = g1_0_zzzzz[k] - abz * g2_0_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPHXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 5))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (15|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 5, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 6, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 5, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SH|g(r,r')|XX) integrals

                    auto g2_0_xxxxx = braBuffer.data(g2off + j);

                    auto g2_0_xxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_xxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_xxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_xxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_0_xxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_0_xxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_0_xxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_0_xxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_0_xyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_0_xyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_0_xyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_0_xyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_0_xzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_0_yyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_0_yyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_0_yyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_0_yyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_0_yzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_0_zzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    // set up pointers to (SI|g(r,r')|XX) integrals

                    auto g1_0_xxxxxx = braBuffer.data(g1off + j);

                    auto g1_0_xxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_xxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_xxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_xxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_xxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_0_xxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_0_xxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_0_xxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_0_xxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_0_xxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_0_xyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_0_xyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_0_xyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_0_xyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_0_xyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_0_xzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_0_yyyyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_0_yyyyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_0_yyyyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_0_yyyzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_0_yyzzzz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_0_yzzzzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_0_zzzzzz = braBuffer.data(g1off + 27 * kcomp + j);

                    // set up pointers to (PH|g(r,r')|XX) integrals

                    auto g_x_xxxxx = braBuffer.data(goff + j);

                    auto g_x_xxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_xxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_xxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_xxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_x_xxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_x_xxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_x_xxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_x_xxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_x_xyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_x_xyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_x_xyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_x_xyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_x_xzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_x_yyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_x_yyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_x_yyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_x_yyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_x_yzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_x_zzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_y_xxxxx = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_y_xxxxy = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_y_xxxxz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_y_xxxyy = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_y_xxxyz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_y_xxxzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_y_xxyyy = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_y_xxyyz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_y_xxyzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_y_xxzzz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_y_xyyyy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_y_xyyyz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_y_xyyzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_y_xyzzz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_y_xzzzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_y_yyyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_y_yyyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_y_yyyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_y_yyzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_y_yzzzz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_y_zzzzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_z_xxxxx = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_z_xxxxy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_z_xxxxz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_z_xxxyy = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_z_xxxyz = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_z_xxxzz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_z_xxyyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_z_xxyyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_z_xxyzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_z_xxzzz = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_z_xyyyy = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_z_xyyyz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_z_xyyzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_z_xyzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_z_xzzzz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_z_yyyyy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_z_yyyyz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_z_yyyzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_z_yyzzz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_z_yzzzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_z_zzzzz = braBuffer.data(goff + 62 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xxxxx, g2_0_xxxxy, g2_0_xxxxz,\
                                             g2_0_xxxyy, g2_0_xxxyz, g2_0_xxxzz,\
                                             g2_0_xxyyy, g2_0_xxyyz, g2_0_xxyzz,\
                                             g2_0_xxzzz, g2_0_xyyyy, g2_0_xyyyz,\
                                             g2_0_xyyzz, g2_0_xyzzz, g2_0_xzzzz,\
                                             g2_0_yyyyy, g2_0_yyyyz, g2_0_yyyzz,\
                                             g2_0_yyzzz, g2_0_yzzzz, g2_0_zzzzz,\
                                             g1_0_xxxxxx, g1_0_xxxxxy, g1_0_xxxxxz,\
                                             g1_0_xxxxyy, g1_0_xxxxyz, g1_0_xxxxzz,\
                                             g1_0_xxxyyy, g1_0_xxxyyz, g1_0_xxxyzz,\
                                             g1_0_xxxzzz, g1_0_xxyyyy, g1_0_xxyyyz,\
                                             g1_0_xxyyzz, g1_0_xxyzzz, g1_0_xxzzzz,\
                                             g1_0_xyyyyy, g1_0_xyyyyz, g1_0_xyyyzz,\
                                             g1_0_xyyzzz, g1_0_xyzzzz, g1_0_xzzzzz,\
                                             g1_0_yyyyyy, g1_0_yyyyyz, g1_0_yyyyzz,\
                                             g1_0_yyyzzz, g1_0_yyzzzz, g1_0_yzzzzz,\
                                             g1_0_zzzzzz, g_x_xxxxx, g_x_xxxxy,\
                                             g_x_xxxxz, g_x_xxxyy, g_x_xxxyz, g_x_xxxzz,\
                                             g_x_xxyyy, g_x_xxyyz, g_x_xxyzz, g_x_xxzzz,\
                                             g_x_xyyyy, g_x_xyyyz, g_x_xyyzz, g_x_xyzzz,\
                                             g_x_xzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz,\
                                             g_x_yyzzz, g_x_yzzzz, g_x_zzzzz, g_y_xxxxx,\
                                             g_y_xxxxy, g_y_xxxxz, g_y_xxxyy, g_y_xxxyz,\
                                             g_y_xxxzz, g_y_xxyyy, g_y_xxyyz, g_y_xxyzz,\
                                             g_y_xxzzz, g_y_xyyyy, g_y_xyyyz, g_y_xyyzz,\
                                             g_y_xyzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyz,\
                                             g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz,\
                                             g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxyy,\
                                             g_z_xxxyz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyz,\
                                             g_z_xxyzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyz,\
                                             g_z_xyyzz, g_z_xyzzz, g_z_xzzzz, g_z_yyyyy,\
                                             g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz,\
                                             g_z_zzzzz: VLX_ALIGN)
                     for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xxxxx[k] = g1_0_xxxxxx[k] - abx * g2_0_xxxxx[k];

                        g_x_xxxxy[k] = g1_0_xxxxxy[k] - abx * g2_0_xxxxy[k];

                        g_x_xxxxz[k] = g1_0_xxxxxz[k] - abx * g2_0_xxxxz[k];

                        g_x_xxxyy[k] = g1_0_xxxxyy[k] - abx * g2_0_xxxyy[k];

                        g_x_xxxyz[k] = g1_0_xxxxyz[k] - abx * g2_0_xxxyz[k];

                        g_x_xxxzz[k] = g1_0_xxxxzz[k] - abx * g2_0_xxxzz[k];

                        g_x_xxyyy[k] = g1_0_xxxyyy[k] - abx * g2_0_xxyyy[k];

                        g_x_xxyyz[k] = g1_0_xxxyyz[k] - abx * g2_0_xxyyz[k];

                        g_x_xxyzz[k] = g1_0_xxxyzz[k] - abx * g2_0_xxyzz[k];

                        g_x_xxzzz[k] = g1_0_xxxzzz[k] - abx * g2_0_xxzzz[k];

                        g_x_xyyyy[k] = g1_0_xxyyyy[k] - abx * g2_0_xyyyy[k];

                        g_x_xyyyz[k] = g1_0_xxyyyz[k] - abx * g2_0_xyyyz[k];

                        g_x_xyyzz[k] = g1_0_xxyyzz[k] - abx * g2_0_xyyzz[k];

                        g_x_xyzzz[k] = g1_0_xxyzzz[k] - abx * g2_0_xyzzz[k];

                        g_x_xzzzz[k] = g1_0_xxzzzz[k] - abx * g2_0_xzzzz[k];

                        g_x_yyyyy[k] = g1_0_xyyyyy[k] - abx * g2_0_yyyyy[k];

                        g_x_yyyyz[k] = g1_0_xyyyyz[k] - abx * g2_0_yyyyz[k];

                        g_x_yyyzz[k] = g1_0_xyyyzz[k] - abx * g2_0_yyyzz[k];

                        g_x_yyzzz[k] = g1_0_xyyzzz[k] - abx * g2_0_yyzzz[k];

                        g_x_yzzzz[k] = g1_0_xyzzzz[k] - abx * g2_0_yzzzz[k];

                        g_x_zzzzz[k] = g1_0_xzzzzz[k] - abx * g2_0_zzzzz[k];

                        // leading y component

                        g_y_xxxxx[k] = g1_0_xxxxxy[k] - aby * g2_0_xxxxx[k];

                        g_y_xxxxy[k] = g1_0_xxxxyy[k] - aby * g2_0_xxxxy[k];

                        g_y_xxxxz[k] = g1_0_xxxxyz[k] - aby * g2_0_xxxxz[k];

                        g_y_xxxyy[k] = g1_0_xxxyyy[k] - aby * g2_0_xxxyy[k];

                        g_y_xxxyz[k] = g1_0_xxxyyz[k] - aby * g2_0_xxxyz[k];

                        g_y_xxxzz[k] = g1_0_xxxyzz[k] - aby * g2_0_xxxzz[k];

                        g_y_xxyyy[k] = g1_0_xxyyyy[k] - aby * g2_0_xxyyy[k];

                        g_y_xxyyz[k] = g1_0_xxyyyz[k] - aby * g2_0_xxyyz[k];

                        g_y_xxyzz[k] = g1_0_xxyyzz[k] - aby * g2_0_xxyzz[k];

                        g_y_xxzzz[k] = g1_0_xxyzzz[k] - aby * g2_0_xxzzz[k];

                        g_y_xyyyy[k] = g1_0_xyyyyy[k] - aby * g2_0_xyyyy[k];

                        g_y_xyyyz[k] = g1_0_xyyyyz[k] - aby * g2_0_xyyyz[k];

                        g_y_xyyzz[k] = g1_0_xyyyzz[k] - aby * g2_0_xyyzz[k];

                        g_y_xyzzz[k] = g1_0_xyyzzz[k] - aby * g2_0_xyzzz[k];

                        g_y_xzzzz[k] = g1_0_xyzzzz[k] - aby * g2_0_xzzzz[k];

                        g_y_yyyyy[k] = g1_0_yyyyyy[k] - aby * g2_0_yyyyy[k];

                        g_y_yyyyz[k] = g1_0_yyyyyz[k] - aby * g2_0_yyyyz[k];

                        g_y_yyyzz[k] = g1_0_yyyyzz[k] - aby * g2_0_yyyzz[k];

                        g_y_yyzzz[k] = g1_0_yyyzzz[k] - aby * g2_0_yyzzz[k];

                        g_y_yzzzz[k] = g1_0_yyzzzz[k] - aby * g2_0_yzzzz[k];

                        g_y_zzzzz[k] = g1_0_yzzzzz[k] - aby * g2_0_zzzzz[k];

                        // leading z component

                        g_z_xxxxx[k] = g1_0_xxxxxz[k] - abz * g2_0_xxxxx[k];

                        g_z_xxxxy[k] = g1_0_xxxxyz[k] - abz * g2_0_xxxxy[k];

                        g_z_xxxxz[k] = g1_0_xxxxzz[k] - abz * g2_0_xxxxz[k];

                        g_z_xxxyy[k] = g1_0_xxxyyz[k] - abz * g2_0_xxxyy[k];

                        g_z_xxxyz[k] = g1_0_xxxyzz[k] - abz * g2_0_xxxyz[k];

                        g_z_xxxzz[k] = g1_0_xxxzzz[k] - abz * g2_0_xxxzz[k];

                        g_z_xxyyy[k] = g1_0_xxyyyz[k] - abz * g2_0_xxyyy[k];

                        g_z_xxyyz[k] = g1_0_xxyyzz[k] - abz * g2_0_xxyyz[k];

                        g_z_xxyzz[k] = g1_0_xxyzzz[k] - abz * g2_0_xxyzz[k];

                        g_z_xxzzz[k] = g1_0_xxzzzz[k] - abz * g2_0_xxzzz[k];

                        g_z_xyyyy[k] = g1_0_xyyyyz[k] - abz * g2_0_xyyyy[k];

                        g_z_xyyyz[k] = g1_0_xyyyzz[k] - abz * g2_0_xyyyz[k];

                        g_z_xyyzz[k] = g1_0_xyyzzz[k] - abz * g2_0_xyyzz[k];

                        g_z_xyzzz[k] = g1_0_xyzzzz[k] - abz * g2_0_xyzzz[k];

                        g_z_xzzzz[k] = g1_0_xzzzzz[k] - abz * g2_0_xzzzz[k];

                        g_z_yyyyy[k] = g1_0_yyyyyz[k] - abz * g2_0_yyyyy[k];

                        g_z_yyyyz[k] = g1_0_yyyyzz[k] - abz * g2_0_yyyyz[k];

                        g_z_yyyzz[k] = g1_0_yyyzzz[k] - abz * g2_0_yyyzz[k];

                        g_z_yyzzz[k] = g1_0_yyzzzz[k] - abz * g2_0_yyzzz[k];

                        g_z_yzzzz[k] = g1_0_yzzzzz[k] - abz * g2_0_yzzzz[k];

                        g_z_zzzzz[k] = g1_0_zzzzzz[k] - abz * g2_0_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPIXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 6))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (16|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 6, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 7, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 6, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SI|g(r,r')|XX) integrals

                    auto g2_0_xxxxxx = braBuffer.data(g2off + j);

                    auto g2_0_xxxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xxxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_xxxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_xxxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_xxxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_0_xxxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_0_xxxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_0_xxxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_0_xxxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_0_xxyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_0_xxyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_0_xxyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_0_xxyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_0_xxzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_0_xyyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_0_xyyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_0_xyyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_0_xyyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_0_xyzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_0_xzzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_0_yyyyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_0_yyyyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_0_yyyyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_0_yyyzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_0_yyzzzz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_0_yzzzzz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_0_zzzzzz = braBuffer.data(g2off + 27 * kcomp + j);

                    // set up pointers to (SK|g(r,r')|XX) integrals

                    auto g1_0_xxxxxxx = braBuffer.data(g1off + j);

                    auto g1_0_xxxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xxxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xxxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xxxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_xxxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_xxxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_xxxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_xxxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_0_xxxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_0_xxxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_0_xxxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_0_xxxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_0_xxxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_0_xxyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_0_xxyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_0_xxyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_0_xxyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_0_xxyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_0_xxzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_0_xyyyyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_0_xyyyyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_0_xyyyyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_0_xyyyzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_0_xyyzzzz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_0_xyzzzzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_0_xzzzzzz = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_0_yyyyyyy = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_0_yyyyyyz = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_0_yyyyyzz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_0_yyyyzzz = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_0_yyyzzzz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_0_yyzzzzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_0_yzzzzzz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_0_zzzzzzz = braBuffer.data(g1off + 35 * kcomp + j);

                    // set up pointers to (PI|g(r,r')|XX) integrals

                    auto g_x_xxxxxx = braBuffer.data(goff + j);

                    auto g_x_xxxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xxxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_xxxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_xxxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_xxxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_x_xxxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_x_xxxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_x_xxxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_x_xxxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_x_xxyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_x_xxyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_x_xxyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_x_xxyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_x_xxzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_x_xyyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_x_xyyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_x_xyyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_x_xyyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_x_xyzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_x_xzzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_x_yyyyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_x_yyyyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_x_yyyyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_x_yyyzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_x_yyzzzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_x_yzzzzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_x_zzzzzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_y_xxxxxx = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_y_xxxxxy = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_y_xxxxxz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_y_xxxxyy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_y_xxxxyz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_y_xxxxzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_y_xxxyyy = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_y_xxxyyz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_y_xxxyzz = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_y_xxxzzz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_y_xxyyyy = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_y_xxyyyz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_y_xxyyzz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_y_xxyzzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_y_xxzzzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_y_xyyyyy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_y_xyyyyz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_y_xyyyzz = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_y_xyyzzz = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_y_xyzzzz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_y_xzzzzz = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_y_yyyyyy = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_y_yyyyyz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_y_yyyyzz = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_y_yyyzzz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_y_yyzzzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_y_yzzzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_y_zzzzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_z_xxxxxx = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_z_xxxxxy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_z_xxxxxz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_z_xxxxyy = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_z_xxxxyz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_z_xxxxzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_z_xxxyyy = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_z_xxxyyz = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_z_xxxyzz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_z_xxxzzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_z_xxyyyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_z_xxyyyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_z_xxyyzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_z_xxyzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_z_xxzzzz = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_z_xyyyyy = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_z_xyyyyz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_z_xyyyzz = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_z_xyyzzz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_z_xyzzzz = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_z_xzzzzz = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_z_yyyyyy = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_z_yyyyyz = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_z_yyyyzz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_z_yyyzzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_z_yyzzzz = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_z_yzzzzz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_z_zzzzzz = braBuffer.data(goff + 83 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xxxxxx, g2_0_xxxxxy, g2_0_xxxxxz,\
                                             g2_0_xxxxyy, g2_0_xxxxyz, g2_0_xxxxzz,\
                                             g2_0_xxxyyy, g2_0_xxxyyz, g2_0_xxxyzz,\
                                             g2_0_xxxzzz, g2_0_xxyyyy, g2_0_xxyyyz,\
                                             g2_0_xxyyzz, g2_0_xxyzzz, g2_0_xxzzzz,\
                                             g2_0_xyyyyy, g2_0_xyyyyz, g2_0_xyyyzz,\
                                             g2_0_xyyzzz, g2_0_xyzzzz, g2_0_xzzzzz,\
                                             g2_0_yyyyyy, g2_0_yyyyyz, g2_0_yyyyzz,\
                                             g2_0_yyyzzz, g2_0_yyzzzz, g2_0_yzzzzz,\
                                             g2_0_zzzzzz, g1_0_xxxxxxx, g1_0_xxxxxxy,\
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
                                             g1_0_zzzzzzz, g_x_xxxxxx, g_x_xxxxxy,\
                                             g_x_xxxxxz, g_x_xxxxyy, g_x_xxxxyz,\
                                             g_x_xxxxzz, g_x_xxxyyy, g_x_xxxyyz,\
                                             g_x_xxxyzz, g_x_xxxzzz, g_x_xxyyyy,\
                                             g_x_xxyyyz, g_x_xxyyzz, g_x_xxyzzz,\
                                             g_x_xxzzzz, g_x_xyyyyy, g_x_xyyyyz,\
                                             g_x_xyyyzz, g_x_xyyzzz, g_x_xyzzzz,\
                                             g_x_xzzzzz, g_x_yyyyyy, g_x_yyyyyz,\
                                             g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz,\
                                             g_x_yzzzzz, g_x_zzzzzz, g_y_xxxxxx,\
                                             g_y_xxxxxy, g_y_xxxxxz, g_y_xxxxyy,\
                                             g_y_xxxxyz, g_y_xxxxzz, g_y_xxxyyy,\
                                             g_y_xxxyyz, g_y_xxxyzz, g_y_xxxzzz,\
                                             g_y_xxyyyy, g_y_xxyyyz, g_y_xxyyzz,\
                                             g_y_xxyzzz, g_y_xxzzzz, g_y_xyyyyy,\
                                             g_y_xyyyyz, g_y_xyyyzz, g_y_xyyzzz,\
                                             g_y_xyzzzz, g_y_xzzzzz, g_y_yyyyyy,\
                                             g_y_yyyyyz, g_y_yyyyzz, g_y_yyyzzz,\
                                             g_y_yyzzzz, g_y_yzzzzz, g_y_zzzzzz,\
                                             g_z_xxxxxx, g_z_xxxxxy, g_z_xxxxxz,\
                                             g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxzz,\
                                             g_z_xxxyyy, g_z_xxxyyz, g_z_xxxyzz,\
                                             g_z_xxxzzz, g_z_xxyyyy, g_z_xxyyyz,\
                                             g_z_xxyyzz, g_z_xxyzzz, g_z_xxzzzz,\
                                             g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyzz,\
                                             g_z_xyyzzz, g_z_xyzzzz, g_z_xzzzzz,\
                                             g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyzz,\
                                             g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz,\
                                             g_z_zzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xxxxxx[k] = g1_0_xxxxxxx[k] - abx * g2_0_xxxxxx[k];

                        g_x_xxxxxy[k] = g1_0_xxxxxxy[k] - abx * g2_0_xxxxxy[k];

                        g_x_xxxxxz[k] = g1_0_xxxxxxz[k] - abx * g2_0_xxxxxz[k];

                        g_x_xxxxyy[k] = g1_0_xxxxxyy[k] - abx * g2_0_xxxxyy[k];

                        g_x_xxxxyz[k] = g1_0_xxxxxyz[k] - abx * g2_0_xxxxyz[k];

                        g_x_xxxxzz[k] = g1_0_xxxxxzz[k] - abx * g2_0_xxxxzz[k];

                        g_x_xxxyyy[k] = g1_0_xxxxyyy[k] - abx * g2_0_xxxyyy[k];

                        g_x_xxxyyz[k] = g1_0_xxxxyyz[k] - abx * g2_0_xxxyyz[k];

                        g_x_xxxyzz[k] = g1_0_xxxxyzz[k] - abx * g2_0_xxxyzz[k];

                        g_x_xxxzzz[k] = g1_0_xxxxzzz[k] - abx * g2_0_xxxzzz[k];

                        g_x_xxyyyy[k] = g1_0_xxxyyyy[k] - abx * g2_0_xxyyyy[k];

                        g_x_xxyyyz[k] = g1_0_xxxyyyz[k] - abx * g2_0_xxyyyz[k];

                        g_x_xxyyzz[k] = g1_0_xxxyyzz[k] - abx * g2_0_xxyyzz[k];

                        g_x_xxyzzz[k] = g1_0_xxxyzzz[k] - abx * g2_0_xxyzzz[k];

                        g_x_xxzzzz[k] = g1_0_xxxzzzz[k] - abx * g2_0_xxzzzz[k];

                        g_x_xyyyyy[k] = g1_0_xxyyyyy[k] - abx * g2_0_xyyyyy[k];

                        g_x_xyyyyz[k] = g1_0_xxyyyyz[k] - abx * g2_0_xyyyyz[k];

                        g_x_xyyyzz[k] = g1_0_xxyyyzz[k] - abx * g2_0_xyyyzz[k];

                        g_x_xyyzzz[k] = g1_0_xxyyzzz[k] - abx * g2_0_xyyzzz[k];

                        g_x_xyzzzz[k] = g1_0_xxyzzzz[k] - abx * g2_0_xyzzzz[k];

                        g_x_xzzzzz[k] = g1_0_xxzzzzz[k] - abx * g2_0_xzzzzz[k];

                        g_x_yyyyyy[k] = g1_0_xyyyyyy[k] - abx * g2_0_yyyyyy[k];

                        g_x_yyyyyz[k] = g1_0_xyyyyyz[k] - abx * g2_0_yyyyyz[k];

                        g_x_yyyyzz[k] = g1_0_xyyyyzz[k] - abx * g2_0_yyyyzz[k];

                        g_x_yyyzzz[k] = g1_0_xyyyzzz[k] - abx * g2_0_yyyzzz[k];

                        g_x_yyzzzz[k] = g1_0_xyyzzzz[k] - abx * g2_0_yyzzzz[k];

                        g_x_yzzzzz[k] = g1_0_xyzzzzz[k] - abx * g2_0_yzzzzz[k];

                        g_x_zzzzzz[k] = g1_0_xzzzzzz[k] - abx * g2_0_zzzzzz[k];

                        // leading y component

                        g_y_xxxxxx[k] = g1_0_xxxxxxy[k] - aby * g2_0_xxxxxx[k];

                        g_y_xxxxxy[k] = g1_0_xxxxxyy[k] - aby * g2_0_xxxxxy[k];

                        g_y_xxxxxz[k] = g1_0_xxxxxyz[k] - aby * g2_0_xxxxxz[k];

                        g_y_xxxxyy[k] = g1_0_xxxxyyy[k] - aby * g2_0_xxxxyy[k];

                        g_y_xxxxyz[k] = g1_0_xxxxyyz[k] - aby * g2_0_xxxxyz[k];

                        g_y_xxxxzz[k] = g1_0_xxxxyzz[k] - aby * g2_0_xxxxzz[k];

                        g_y_xxxyyy[k] = g1_0_xxxyyyy[k] - aby * g2_0_xxxyyy[k];

                        g_y_xxxyyz[k] = g1_0_xxxyyyz[k] - aby * g2_0_xxxyyz[k];

                        g_y_xxxyzz[k] = g1_0_xxxyyzz[k] - aby * g2_0_xxxyzz[k];

                        g_y_xxxzzz[k] = g1_0_xxxyzzz[k] - aby * g2_0_xxxzzz[k];

                        g_y_xxyyyy[k] = g1_0_xxyyyyy[k] - aby * g2_0_xxyyyy[k];

                        g_y_xxyyyz[k] = g1_0_xxyyyyz[k] - aby * g2_0_xxyyyz[k];

                        g_y_xxyyzz[k] = g1_0_xxyyyzz[k] - aby * g2_0_xxyyzz[k];

                        g_y_xxyzzz[k] = g1_0_xxyyzzz[k] - aby * g2_0_xxyzzz[k];

                        g_y_xxzzzz[k] = g1_0_xxyzzzz[k] - aby * g2_0_xxzzzz[k];

                        g_y_xyyyyy[k] = g1_0_xyyyyyy[k] - aby * g2_0_xyyyyy[k];

                        g_y_xyyyyz[k] = g1_0_xyyyyyz[k] - aby * g2_0_xyyyyz[k];

                        g_y_xyyyzz[k] = g1_0_xyyyyzz[k] - aby * g2_0_xyyyzz[k];

                        g_y_xyyzzz[k] = g1_0_xyyyzzz[k] - aby * g2_0_xyyzzz[k];

                        g_y_xyzzzz[k] = g1_0_xyyzzzz[k] - aby * g2_0_xyzzzz[k];

                        g_y_xzzzzz[k] = g1_0_xyzzzzz[k] - aby * g2_0_xzzzzz[k];

                        g_y_yyyyyy[k] = g1_0_yyyyyyy[k] - aby * g2_0_yyyyyy[k];

                        g_y_yyyyyz[k] = g1_0_yyyyyyz[k] - aby * g2_0_yyyyyz[k];

                        g_y_yyyyzz[k] = g1_0_yyyyyzz[k] - aby * g2_0_yyyyzz[k];

                        g_y_yyyzzz[k] = g1_0_yyyyzzz[k] - aby * g2_0_yyyzzz[k];

                        g_y_yyzzzz[k] = g1_0_yyyzzzz[k] - aby * g2_0_yyzzzz[k];

                        g_y_yzzzzz[k] = g1_0_yyzzzzz[k] - aby * g2_0_yzzzzz[k];

                        g_y_zzzzzz[k] = g1_0_yzzzzzz[k] - aby * g2_0_zzzzzz[k];

                        // leading z component

                        g_z_xxxxxx[k] = g1_0_xxxxxxz[k] - abz * g2_0_xxxxxx[k];

                        g_z_xxxxxy[k] = g1_0_xxxxxyz[k] - abz * g2_0_xxxxxy[k];

                        g_z_xxxxxz[k] = g1_0_xxxxxzz[k] - abz * g2_0_xxxxxz[k];

                        g_z_xxxxyy[k] = g1_0_xxxxyyz[k] - abz * g2_0_xxxxyy[k];

                        g_z_xxxxyz[k] = g1_0_xxxxyzz[k] - abz * g2_0_xxxxyz[k];

                        g_z_xxxxzz[k] = g1_0_xxxxzzz[k] - abz * g2_0_xxxxzz[k];

                        g_z_xxxyyy[k] = g1_0_xxxyyyz[k] - abz * g2_0_xxxyyy[k];

                        g_z_xxxyyz[k] = g1_0_xxxyyzz[k] - abz * g2_0_xxxyyz[k];

                        g_z_xxxyzz[k] = g1_0_xxxyzzz[k] - abz * g2_0_xxxyzz[k];

                        g_z_xxxzzz[k] = g1_0_xxxzzzz[k] - abz * g2_0_xxxzzz[k];

                        g_z_xxyyyy[k] = g1_0_xxyyyyz[k] - abz * g2_0_xxyyyy[k];

                        g_z_xxyyyz[k] = g1_0_xxyyyzz[k] - abz * g2_0_xxyyyz[k];

                        g_z_xxyyzz[k] = g1_0_xxyyzzz[k] - abz * g2_0_xxyyzz[k];

                        g_z_xxyzzz[k] = g1_0_xxyzzzz[k] - abz * g2_0_xxyzzz[k];

                        g_z_xxzzzz[k] = g1_0_xxzzzzz[k] - abz * g2_0_xxzzzz[k];

                        g_z_xyyyyy[k] = g1_0_xyyyyyz[k] - abz * g2_0_xyyyyy[k];

                        g_z_xyyyyz[k] = g1_0_xyyyyzz[k] - abz * g2_0_xyyyyz[k];

                        g_z_xyyyzz[k] = g1_0_xyyyzzz[k] - abz * g2_0_xyyyzz[k];

                        g_z_xyyzzz[k] = g1_0_xyyzzzz[k] - abz * g2_0_xyyzzz[k];

                        g_z_xyzzzz[k] = g1_0_xyzzzzz[k] - abz * g2_0_xyzzzz[k];

                        g_z_xzzzzz[k] = g1_0_xzzzzzz[k] - abz * g2_0_xzzzzz[k];

                        g_z_yyyyyy[k] = g1_0_yyyyyyz[k] - abz * g2_0_yyyyyy[k];

                        g_z_yyyyyz[k] = g1_0_yyyyyzz[k] - abz * g2_0_yyyyyz[k];

                        g_z_yyyyzz[k] = g1_0_yyyyzzz[k] - abz * g2_0_yyyyzz[k];

                        g_z_yyyzzz[k] = g1_0_yyyzzzz[k] - abz * g2_0_yyyzzz[k];

                        g_z_yyzzzz[k] = g1_0_yyzzzzz[k] - abz * g2_0_yyzzzz[k];

                        g_z_yzzzzz[k] = g1_0_yzzzzzz[k] - abz * g2_0_yzzzzz[k];

                        g_z_zzzzzz[k] = g1_0_zzzzzzz[k] - abz * g2_0_zzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForPKXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 1) && (recPattern[i].second() == 7))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (17|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {1, 7, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 8, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {0, 7, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (SK|g(r,r')|XX) integrals

                    auto g2_0_xxxxxxx = braBuffer.data(g2off + j);

                    auto g2_0_xxxxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_0_xxxxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_0_xxxxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_0_xxxxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_0_xxxxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_0_xxxxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_0_xxxxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_0_xxxxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_0_xxxxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_0_xxxyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_0_xxxyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_0_xxxyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_0_xxxyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_0_xxxzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_0_xxyyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_0_xxyyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_0_xxyyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_0_xxyyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_0_xxyzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_0_xxzzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_0_xyyyyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_0_xyyyyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_0_xyyyyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_0_xyyyzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_0_xyyzzzz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_0_xyzzzzz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_0_xzzzzzz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_0_yyyyyyy = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_0_yyyyyyz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_0_yyyyyzz = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_0_yyyyzzz = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_0_yyyzzzz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_0_yyzzzzz = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_0_yzzzzzz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_0_zzzzzzz = braBuffer.data(g2off + 35 * kcomp + j);

                    // set up pointers to (SL|g(r,r')|XX) integrals

                    auto g1_0_xxxxxxxx = braBuffer.data(g1off + j);

                    auto g1_0_xxxxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_0_xxxxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_0_xxxxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_0_xxxxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_0_xxxxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_0_xxxxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_0_xxxxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_0_xxxxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_0_xxxxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_0_xxxxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_0_xxxxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_0_xxxxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_0_xxxxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_0_xxxxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_0_xxxyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_0_xxxyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_0_xxxyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_0_xxxyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_0_xxxyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_0_xxxzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_0_xxyyyyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_0_xxyyyyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_0_xxyyyyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_0_xxyyyzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_0_xxyyzzzz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_0_xxyzzzzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_0_xxzzzzzz = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_0_xyyyyyyy = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_0_xyyyyyyz = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_0_xyyyyyzz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_0_xyyyyzzz = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_0_xyyyzzzz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_0_xyyzzzzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_0_xyzzzzzz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_0_xzzzzzzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_0_yyyyyyyy = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_0_yyyyyyyz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_0_yyyyyyzz = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_0_yyyyyzzz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_0_yyyyzzzz = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_0_yyyzzzzz = braBuffer.data(g1off + 41 * kcomp + j);

                    auto g1_0_yyzzzzzz = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_0_yzzzzzzz = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_0_zzzzzzzz = braBuffer.data(g1off + 44 * kcomp + j);

                    // set up pointers to (PK|g(r,r')|XX) integrals

                    auto g_x_xxxxxxx = braBuffer.data(goff + j);

                    auto g_x_xxxxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_x_xxxxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_x_xxxxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_x_xxxxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_x_xxxxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_x_xxxxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_x_xxxxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_x_xxxxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_x_xxxxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_x_xxxyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_x_xxxyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_x_xxxyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_x_xxxyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_x_xxxzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_x_xxyyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_x_xxyyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_x_xxyyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_x_xxyyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_x_xxyzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_x_xxzzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_x_xyyyyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_x_xyyyyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_x_xyyyyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_x_xyyyzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_x_xyyzzzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_x_xyzzzzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_x_xzzzzzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_x_yyyyyyy = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_x_yyyyyyz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_x_yyyyyzz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_x_yyyyzzz = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_x_yyyzzzz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_x_yyzzzzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_x_yzzzzzz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_x_zzzzzzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_y_xxxxxxx = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_y_xxxxxxy = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_y_xxxxxxz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_y_xxxxxyy = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_y_xxxxxyz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_y_xxxxxzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_y_xxxxyyy = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_y_xxxxyyz = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_y_xxxxyzz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_y_xxxxzzz = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_y_xxxyyyy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_y_xxxyyyz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_y_xxxyyzz = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_y_xxxyzzz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_y_xxxzzzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_y_xxyyyyy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_y_xxyyyyz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_y_xxyyyzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_y_xxyyzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_y_xxyzzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_y_xxzzzzz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_y_xyyyyyy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_y_xyyyyyz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_y_xyyyyzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_y_xyyyzzz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_y_xyyzzzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_y_xyzzzzz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_y_xzzzzzz = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_y_yyyyyyy = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_y_yyyyyyz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_y_yyyyyzz = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_y_yyyyzzz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_y_yyyzzzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_y_yyzzzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_y_yzzzzzz = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_y_zzzzzzz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_z_xxxxxxx = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_z_xxxxxxy = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_z_xxxxxxz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_z_xxxxxyy = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_z_xxxxxyz = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_z_xxxxxzz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_z_xxxxyyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_z_xxxxyyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_z_xxxxyzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_z_xxxxzzz = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_z_xxxyyyy = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_z_xxxyyyz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_z_xxxyyzz = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_z_xxxyzzz = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_z_xxxzzzz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_z_xxyyyyy = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_z_xxyyyyz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_z_xxyyyzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_z_xxyyzzz = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_z_xxyzzzz = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_z_xxzzzzz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_z_xyyyyyy = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_z_xyyyyyz = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_z_xyyyyzz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_z_xyyyzzz = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_z_xyyzzzz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_z_xyzzzzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_z_xzzzzzz = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_z_yyyyyyy = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_z_yyyyyyz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_z_yyyyyzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_z_yyyyzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_z_yyyzzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_z_yyzzzzz = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_z_yzzzzzz = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_z_zzzzzzz = braBuffer.data(goff + 107 * kcomp + j);

                    #pragma omp simd aligned(g2_0_xxxxxxx, g2_0_xxxxxxy, g2_0_xxxxxxz,\
                                             g2_0_xxxxxyy, g2_0_xxxxxyz, g2_0_xxxxxzz,\
                                             g2_0_xxxxyyy, g2_0_xxxxyyz, g2_0_xxxxyzz,\
                                             g2_0_xxxxzzz, g2_0_xxxyyyy, g2_0_xxxyyyz,\
                                             g2_0_xxxyyzz, g2_0_xxxyzzz, g2_0_xxxzzzz,\
                                             g2_0_xxyyyyy, g2_0_xxyyyyz, g2_0_xxyyyzz,\
                                             g2_0_xxyyzzz, g2_0_xxyzzzz, g2_0_xxzzzzz,\
                                             g2_0_xyyyyyy, g2_0_xyyyyyz, g2_0_xyyyyzz,\
                                             g2_0_xyyyzzz, g2_0_xyyzzzz, g2_0_xyzzzzz,\
                                             g2_0_xzzzzzz, g2_0_yyyyyyy, g2_0_yyyyyyz,\
                                             g2_0_yyyyyzz, g2_0_yyyyzzz, g2_0_yyyzzzz,\
                                             g2_0_yyzzzzz, g2_0_yzzzzzz, g2_0_zzzzzzz,\
                                             g1_0_xxxxxxxx, g1_0_xxxxxxxy, g1_0_xxxxxxxz,\
                                             g1_0_xxxxxxyy, g1_0_xxxxxxyz, g1_0_xxxxxxzz,\
                                             g1_0_xxxxxyyy, g1_0_xxxxxyyz, g1_0_xxxxxyzz,\
                                             g1_0_xxxxxzzz, g1_0_xxxxyyyy, g1_0_xxxxyyyz,\
                                             g1_0_xxxxyyzz, g1_0_xxxxyzzz, g1_0_xxxxzzzz,\
                                             g1_0_xxxyyyyy, g1_0_xxxyyyyz, g1_0_xxxyyyzz,\
                                             g1_0_xxxyyzzz, g1_0_xxxyzzzz, g1_0_xxxzzzzz,\
                                             g1_0_xxyyyyyy, g1_0_xxyyyyyz, g1_0_xxyyyyzz,\
                                             g1_0_xxyyyzzz, g1_0_xxyyzzzz, g1_0_xxyzzzzz,\
                                             g1_0_xxzzzzzz, g1_0_xyyyyyyy, g1_0_xyyyyyyz,\
                                             g1_0_xyyyyyzz, g1_0_xyyyyzzz, g1_0_xyyyzzzz,\
                                             g1_0_xyyzzzzz, g1_0_xyzzzzzz, g1_0_xzzzzzzz,\
                                             g1_0_yyyyyyyy, g1_0_yyyyyyyz, g1_0_yyyyyyzz,\
                                             g1_0_yyyyyzzz, g1_0_yyyyzzzz, g1_0_yyyzzzzz,\
                                             g1_0_yyzzzzzz, g1_0_yzzzzzzz, g1_0_zzzzzzzz,\
                                             g_x_xxxxxxx, g_x_xxxxxxy, g_x_xxxxxxz,\
                                             g_x_xxxxxyy, g_x_xxxxxyz, g_x_xxxxxzz,\
                                             g_x_xxxxyyy, g_x_xxxxyyz, g_x_xxxxyzz,\
                                             g_x_xxxxzzz, g_x_xxxyyyy, g_x_xxxyyyz,\
                                             g_x_xxxyyzz, g_x_xxxyzzz, g_x_xxxzzzz,\
                                             g_x_xxyyyyy, g_x_xxyyyyz, g_x_xxyyyzz,\
                                             g_x_xxyyzzz, g_x_xxyzzzz, g_x_xxzzzzz,\
                                             g_x_xyyyyyy, g_x_xyyyyyz, g_x_xyyyyzz,\
                                             g_x_xyyyzzz, g_x_xyyzzzz, g_x_xyzzzzz,\
                                             g_x_xzzzzzz, g_x_yyyyyyy, g_x_yyyyyyz,\
                                             g_x_yyyyyzz, g_x_yyyyzzz, g_x_yyyzzzz,\
                                             g_x_yyzzzzz, g_x_yzzzzzz, g_x_zzzzzzz,\
                                             g_y_xxxxxxx, g_y_xxxxxxy, g_y_xxxxxxz,\
                                             g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxzz,\
                                             g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyzz,\
                                             g_y_xxxxzzz, g_y_xxxyyyy, g_y_xxxyyyz,\
                                             g_y_xxxyyzz, g_y_xxxyzzz, g_y_xxxzzzz,\
                                             g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyzz,\
                                             g_y_xxyyzzz, g_y_xxyzzzz, g_y_xxzzzzz,\
                                             g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyzz,\
                                             g_y_xyyyzzz, g_y_xyyzzzz, g_y_xyzzzzz,\
                                             g_y_xzzzzzz, g_y_yyyyyyy, g_y_yyyyyyz,\
                                             g_y_yyyyyzz, g_y_yyyyzzz, g_y_yyyzzzz,\
                                             g_y_yyzzzzz, g_y_yzzzzzz, g_y_zzzzzzz,\
                                             g_z_xxxxxxx, g_z_xxxxxxy, g_z_xxxxxxz,\
                                             g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxzz,\
                                             g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyzz,\
                                             g_z_xxxxzzz, g_z_xxxyyyy, g_z_xxxyyyz,\
                                             g_z_xxxyyzz, g_z_xxxyzzz, g_z_xxxzzzz,\
                                             g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyzz,\
                                             g_z_xxyyzzz, g_z_xxyzzzz, g_z_xxzzzzz,\
                                             g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyzz,\
                                             g_z_xyyyzzz, g_z_xyyzzzz, g_z_xyzzzzz,\
                                             g_z_xzzzzzz, g_z_yyyyyyy, g_z_yyyyyyz,\
                                             g_z_yyyyyzz, g_z_yyyyzzz, g_z_yyyzzzz,\
                                             g_z_yyzzzzz, g_z_yzzzzzz, g_z_zzzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_x_xxxxxxx[k] = g1_0_xxxxxxxx[k] - abx * g2_0_xxxxxxx[k];

                        g_x_xxxxxxy[k] = g1_0_xxxxxxxy[k] - abx * g2_0_xxxxxxy[k];

                        g_x_xxxxxxz[k] = g1_0_xxxxxxxz[k] - abx * g2_0_xxxxxxz[k];

                        g_x_xxxxxyy[k] = g1_0_xxxxxxyy[k] - abx * g2_0_xxxxxyy[k];

                        g_x_xxxxxyz[k] = g1_0_xxxxxxyz[k] - abx * g2_0_xxxxxyz[k];

                        g_x_xxxxxzz[k] = g1_0_xxxxxxzz[k] - abx * g2_0_xxxxxzz[k];

                        g_x_xxxxyyy[k] = g1_0_xxxxxyyy[k] - abx * g2_0_xxxxyyy[k];

                        g_x_xxxxyyz[k] = g1_0_xxxxxyyz[k] - abx * g2_0_xxxxyyz[k];

                        g_x_xxxxyzz[k] = g1_0_xxxxxyzz[k] - abx * g2_0_xxxxyzz[k];

                        g_x_xxxxzzz[k] = g1_0_xxxxxzzz[k] - abx * g2_0_xxxxzzz[k];

                        g_x_xxxyyyy[k] = g1_0_xxxxyyyy[k] - abx * g2_0_xxxyyyy[k];

                        g_x_xxxyyyz[k] = g1_0_xxxxyyyz[k] - abx * g2_0_xxxyyyz[k];

                        g_x_xxxyyzz[k] = g1_0_xxxxyyzz[k] - abx * g2_0_xxxyyzz[k];

                        g_x_xxxyzzz[k] = g1_0_xxxxyzzz[k] - abx * g2_0_xxxyzzz[k];

                        g_x_xxxzzzz[k] = g1_0_xxxxzzzz[k] - abx * g2_0_xxxzzzz[k];

                        g_x_xxyyyyy[k] = g1_0_xxxyyyyy[k] - abx * g2_0_xxyyyyy[k];

                        g_x_xxyyyyz[k] = g1_0_xxxyyyyz[k] - abx * g2_0_xxyyyyz[k];

                        g_x_xxyyyzz[k] = g1_0_xxxyyyzz[k] - abx * g2_0_xxyyyzz[k];

                        g_x_xxyyzzz[k] = g1_0_xxxyyzzz[k] - abx * g2_0_xxyyzzz[k];

                        g_x_xxyzzzz[k] = g1_0_xxxyzzzz[k] - abx * g2_0_xxyzzzz[k];

                        g_x_xxzzzzz[k] = g1_0_xxxzzzzz[k] - abx * g2_0_xxzzzzz[k];

                        g_x_xyyyyyy[k] = g1_0_xxyyyyyy[k] - abx * g2_0_xyyyyyy[k];

                        g_x_xyyyyyz[k] = g1_0_xxyyyyyz[k] - abx * g2_0_xyyyyyz[k];

                        g_x_xyyyyzz[k] = g1_0_xxyyyyzz[k] - abx * g2_0_xyyyyzz[k];

                        g_x_xyyyzzz[k] = g1_0_xxyyyzzz[k] - abx * g2_0_xyyyzzz[k];

                        g_x_xyyzzzz[k] = g1_0_xxyyzzzz[k] - abx * g2_0_xyyzzzz[k];

                        g_x_xyzzzzz[k] = g1_0_xxyzzzzz[k] - abx * g2_0_xyzzzzz[k];

                        g_x_xzzzzzz[k] = g1_0_xxzzzzzz[k] - abx * g2_0_xzzzzzz[k];

                        g_x_yyyyyyy[k] = g1_0_xyyyyyyy[k] - abx * g2_0_yyyyyyy[k];

                        g_x_yyyyyyz[k] = g1_0_xyyyyyyz[k] - abx * g2_0_yyyyyyz[k];

                        g_x_yyyyyzz[k] = g1_0_xyyyyyzz[k] - abx * g2_0_yyyyyzz[k];

                        g_x_yyyyzzz[k] = g1_0_xyyyyzzz[k] - abx * g2_0_yyyyzzz[k];

                        g_x_yyyzzzz[k] = g1_0_xyyyzzzz[k] - abx * g2_0_yyyzzzz[k];

                        g_x_yyzzzzz[k] = g1_0_xyyzzzzz[k] - abx * g2_0_yyzzzzz[k];

                        g_x_yzzzzzz[k] = g1_0_xyzzzzzz[k] - abx * g2_0_yzzzzzz[k];

                        g_x_zzzzzzz[k] = g1_0_xzzzzzzz[k] - abx * g2_0_zzzzzzz[k];

                        // leading y component

                        g_y_xxxxxxx[k] = g1_0_xxxxxxxy[k] - aby * g2_0_xxxxxxx[k];

                        g_y_xxxxxxy[k] = g1_0_xxxxxxyy[k] - aby * g2_0_xxxxxxy[k];

                        g_y_xxxxxxz[k] = g1_0_xxxxxxyz[k] - aby * g2_0_xxxxxxz[k];

                        g_y_xxxxxyy[k] = g1_0_xxxxxyyy[k] - aby * g2_0_xxxxxyy[k];

                        g_y_xxxxxyz[k] = g1_0_xxxxxyyz[k] - aby * g2_0_xxxxxyz[k];

                        g_y_xxxxxzz[k] = g1_0_xxxxxyzz[k] - aby * g2_0_xxxxxzz[k];

                        g_y_xxxxyyy[k] = g1_0_xxxxyyyy[k] - aby * g2_0_xxxxyyy[k];

                        g_y_xxxxyyz[k] = g1_0_xxxxyyyz[k] - aby * g2_0_xxxxyyz[k];

                        g_y_xxxxyzz[k] = g1_0_xxxxyyzz[k] - aby * g2_0_xxxxyzz[k];

                        g_y_xxxxzzz[k] = g1_0_xxxxyzzz[k] - aby * g2_0_xxxxzzz[k];

                        g_y_xxxyyyy[k] = g1_0_xxxyyyyy[k] - aby * g2_0_xxxyyyy[k];

                        g_y_xxxyyyz[k] = g1_0_xxxyyyyz[k] - aby * g2_0_xxxyyyz[k];

                        g_y_xxxyyzz[k] = g1_0_xxxyyyzz[k] - aby * g2_0_xxxyyzz[k];

                        g_y_xxxyzzz[k] = g1_0_xxxyyzzz[k] - aby * g2_0_xxxyzzz[k];

                        g_y_xxxzzzz[k] = g1_0_xxxyzzzz[k] - aby * g2_0_xxxzzzz[k];

                        g_y_xxyyyyy[k] = g1_0_xxyyyyyy[k] - aby * g2_0_xxyyyyy[k];

                        g_y_xxyyyyz[k] = g1_0_xxyyyyyz[k] - aby * g2_0_xxyyyyz[k];

                        g_y_xxyyyzz[k] = g1_0_xxyyyyzz[k] - aby * g2_0_xxyyyzz[k];

                        g_y_xxyyzzz[k] = g1_0_xxyyyzzz[k] - aby * g2_0_xxyyzzz[k];

                        g_y_xxyzzzz[k] = g1_0_xxyyzzzz[k] - aby * g2_0_xxyzzzz[k];

                        g_y_xxzzzzz[k] = g1_0_xxyzzzzz[k] - aby * g2_0_xxzzzzz[k];

                        g_y_xyyyyyy[k] = g1_0_xyyyyyyy[k] - aby * g2_0_xyyyyyy[k];

                        g_y_xyyyyyz[k] = g1_0_xyyyyyyz[k] - aby * g2_0_xyyyyyz[k];

                        g_y_xyyyyzz[k] = g1_0_xyyyyyzz[k] - aby * g2_0_xyyyyzz[k];

                        g_y_xyyyzzz[k] = g1_0_xyyyyzzz[k] - aby * g2_0_xyyyzzz[k];

                        g_y_xyyzzzz[k] = g1_0_xyyyzzzz[k] - aby * g2_0_xyyzzzz[k];

                        g_y_xyzzzzz[k] = g1_0_xyyzzzzz[k] - aby * g2_0_xyzzzzz[k];

                        g_y_xzzzzzz[k] = g1_0_xyzzzzzz[k] - aby * g2_0_xzzzzzz[k];

                        g_y_yyyyyyy[k] = g1_0_yyyyyyyy[k] - aby * g2_0_yyyyyyy[k];

                        g_y_yyyyyyz[k] = g1_0_yyyyyyyz[k] - aby * g2_0_yyyyyyz[k];

                        g_y_yyyyyzz[k] = g1_0_yyyyyyzz[k] - aby * g2_0_yyyyyzz[k];

                        g_y_yyyyzzz[k] = g1_0_yyyyyzzz[k] - aby * g2_0_yyyyzzz[k];

                        g_y_yyyzzzz[k] = g1_0_yyyyzzzz[k] - aby * g2_0_yyyzzzz[k];

                        g_y_yyzzzzz[k] = g1_0_yyyzzzzz[k] - aby * g2_0_yyzzzzz[k];

                        g_y_yzzzzzz[k] = g1_0_yyzzzzzz[k] - aby * g2_0_yzzzzzz[k];

                        g_y_zzzzzzz[k] = g1_0_yzzzzzzz[k] - aby * g2_0_zzzzzzz[k];

                        // leading z component

                        g_z_xxxxxxx[k] = g1_0_xxxxxxxz[k] - abz * g2_0_xxxxxxx[k];

                        g_z_xxxxxxy[k] = g1_0_xxxxxxyz[k] - abz * g2_0_xxxxxxy[k];

                        g_z_xxxxxxz[k] = g1_0_xxxxxxzz[k] - abz * g2_0_xxxxxxz[k];

                        g_z_xxxxxyy[k] = g1_0_xxxxxyyz[k] - abz * g2_0_xxxxxyy[k];

                        g_z_xxxxxyz[k] = g1_0_xxxxxyzz[k] - abz * g2_0_xxxxxyz[k];

                        g_z_xxxxxzz[k] = g1_0_xxxxxzzz[k] - abz * g2_0_xxxxxzz[k];

                        g_z_xxxxyyy[k] = g1_0_xxxxyyyz[k] - abz * g2_0_xxxxyyy[k];

                        g_z_xxxxyyz[k] = g1_0_xxxxyyzz[k] - abz * g2_0_xxxxyyz[k];

                        g_z_xxxxyzz[k] = g1_0_xxxxyzzz[k] - abz * g2_0_xxxxyzz[k];

                        g_z_xxxxzzz[k] = g1_0_xxxxzzzz[k] - abz * g2_0_xxxxzzz[k];

                        g_z_xxxyyyy[k] = g1_0_xxxyyyyz[k] - abz * g2_0_xxxyyyy[k];

                        g_z_xxxyyyz[k] = g1_0_xxxyyyzz[k] - abz * g2_0_xxxyyyz[k];

                        g_z_xxxyyzz[k] = g1_0_xxxyyzzz[k] - abz * g2_0_xxxyyzz[k];

                        g_z_xxxyzzz[k] = g1_0_xxxyzzzz[k] - abz * g2_0_xxxyzzz[k];

                        g_z_xxxzzzz[k] = g1_0_xxxzzzzz[k] - abz * g2_0_xxxzzzz[k];

                        g_z_xxyyyyy[k] = g1_0_xxyyyyyz[k] - abz * g2_0_xxyyyyy[k];

                        g_z_xxyyyyz[k] = g1_0_xxyyyyzz[k] - abz * g2_0_xxyyyyz[k];

                        g_z_xxyyyzz[k] = g1_0_xxyyyzzz[k] - abz * g2_0_xxyyyzz[k];

                        g_z_xxyyzzz[k] = g1_0_xxyyzzzz[k] - abz * g2_0_xxyyzzz[k];

                        g_z_xxyzzzz[k] = g1_0_xxyzzzzz[k] - abz * g2_0_xxyzzzz[k];

                        g_z_xxzzzzz[k] = g1_0_xxzzzzzz[k] - abz * g2_0_xxzzzzz[k];

                        g_z_xyyyyyy[k] = g1_0_xyyyyyyz[k] - abz * g2_0_xyyyyyy[k];

                        g_z_xyyyyyz[k] = g1_0_xyyyyyzz[k] - abz * g2_0_xyyyyyz[k];

                        g_z_xyyyyzz[k] = g1_0_xyyyyzzz[k] - abz * g2_0_xyyyyzz[k];

                        g_z_xyyyzzz[k] = g1_0_xyyyzzzz[k] - abz * g2_0_xyyyzzz[k];

                        g_z_xyyzzzz[k] = g1_0_xyyzzzzz[k] - abz * g2_0_xyyzzzz[k];

                        g_z_xyzzzzz[k] = g1_0_xyzzzzzz[k] - abz * g2_0_xyzzzzz[k];

                        g_z_xzzzzzz[k] = g1_0_xzzzzzzz[k] - abz * g2_0_xzzzzzz[k];

                        g_z_yyyyyyy[k] = g1_0_yyyyyyyz[k] - abz * g2_0_yyyyyyy[k];

                        g_z_yyyyyyz[k] = g1_0_yyyyyyzz[k] - abz * g2_0_yyyyyyz[k];

                        g_z_yyyyyzz[k] = g1_0_yyyyyzzz[k] - abz * g2_0_yyyyyzz[k];

                        g_z_yyyyzzz[k] = g1_0_yyyyzzzz[k] - abz * g2_0_yyyyzzz[k];

                        g_z_yyyzzzz[k] = g1_0_yyyzzzzz[k] - abz * g2_0_yyyzzzz[k];

                        g_z_yyzzzzz[k] = g1_0_yyzzzzzz[k] - abz * g2_0_yyzzzzz[k];

                        g_z_yzzzzzz[k] = g1_0_yzzzzzzz[k] - abz * g2_0_yzzzzzz[k];

                        g_z_zzzzzzz[k] = g1_0_zzzzzzzz[k] - abz * g2_0_zzzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForDDXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 2) && (recPattern[i].second() == 2))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (22|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {2, 2, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 3, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 2, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (PD|g(r,r')|XX) integrals

                    auto g2_x_xx = braBuffer.data(g2off + j);

                    auto g2_x_xy = braBuffer.data(g2off + kcomp + j);

                    auto g2_x_xz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_x_yy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_x_yz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_x_zz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_y_xx = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_y_xy = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_y_xz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_y_yy = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_y_yz = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_y_zz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_z_xx = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_z_xy = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_z_xz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_z_yy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_z_yz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_z_zz = braBuffer.data(g2off + 17 * kcomp + j);

                    // set up pointers to (PF|g(r,r')|XX) integrals

                    auto g1_x_xxx = braBuffer.data(g1off + j);

                    auto g1_x_xxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_x_xxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_x_xyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_x_xyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_x_xzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_y_xxx = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_y_xxy = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_y_xxz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_y_xyy = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_y_xyz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_y_xzz = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_y_yyy = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_y_yyz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_y_yzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_z_xxx = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_z_xxy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_z_xxz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_z_xyy = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_z_xyz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_z_xzz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_z_yyy = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_z_yyz = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_z_yzz = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_z_zzz = braBuffer.data(g1off + 29 * kcomp + j);

                    // set up pointers to (DD|g(r,r')|XX) integrals

                    auto g_xx_xx = braBuffer.data(goff + j);

                    auto g_xx_xy = braBuffer.data(goff + kcomp + j);

                    auto g_xx_xz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xx_yy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xx_yz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xx_zz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xy_xx = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xy_xy = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xy_xz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xy_yy = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xy_yz = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xy_zz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xz_xx = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xz_xy = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xz_xz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xz_yy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xz_yz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xz_zz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_yy_xx = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_yy_xy = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_yy_xz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_yy_yy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_yy_yz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_yy_zz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_yz_xx = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_yz_xy = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_yz_xz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_yz_yy = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_yz_yz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_yz_zz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_zz_xx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_zz_xy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_zz_xz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_zz_yy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_zz_yz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_zz_zz = braBuffer.data(goff + 35 * kcomp + j);

                    #pragma omp simd aligned(g2_x_xx, g2_x_xy, g2_x_xz, g2_x_yy,\
                                             g2_x_yz, g2_x_zz, g2_y_xx, g2_y_xy,\
                                             g2_y_xz, g2_y_yy, g2_y_yz, g2_y_zz,\
                                             g2_z_xx, g2_z_xy, g2_z_xz, g2_z_yy,\
                                             g2_z_yz, g2_z_zz, g1_x_xxx, g1_x_xxy,\
                                             g1_x_xxz, g1_x_xyy, g1_x_xyz, g1_x_xzz,\
                                             g1_y_xxx, g1_y_xxy, g1_y_xxz, g1_y_xyy,\
                                             g1_y_xyz, g1_y_xzz, g1_y_yyy, g1_y_yyz,\
                                             g1_y_yzz, g1_z_xxx, g1_z_xxy,\
                                             g1_z_xxz, g1_z_xyy, g1_z_xyz, g1_z_xzz,\
                                             g1_z_yyy, g1_z_yyz, g1_z_yzz, g1_z_zzz,\
                                             g_xx_xx, g_xx_xy, g_xx_xz, g_xx_yy,\
                                             g_xx_yz, g_xx_zz, g_xy_xx, g_xy_xy,\
                                             g_xy_xz, g_xy_yy, g_xy_yz, g_xy_zz,\
                                             g_xz_xx, g_xz_xy, g_xz_xz, g_xz_yy,\
                                             g_xz_yz, g_xz_zz, g_yy_xx, g_yy_xy,\
                                             g_yy_xz, g_yy_yy, g_yy_yz, g_yy_zz,\
                                             g_yz_xx, g_yz_xy, g_yz_xz, g_yz_yy,\
                                             g_yz_yz, g_yz_zz, g_zz_xx, g_zz_xy,\
                                             g_zz_xz, g_zz_yy, g_zz_yz, g_zz_zz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xx_xx[k] = g1_x_xxx[k] - abx * g2_x_xx[k];

                        g_xx_xy[k] = g1_x_xxy[k] - abx * g2_x_xy[k];

                        g_xx_xz[k] = g1_x_xxz[k] - abx * g2_x_xz[k];

                        g_xx_yy[k] = g1_x_xyy[k] - abx * g2_x_yy[k];

                        g_xx_yz[k] = g1_x_xyz[k] - abx * g2_x_yz[k];

                        g_xx_zz[k] = g1_x_xzz[k] - abx * g2_x_zz[k];

                        g_xy_xx[k] = g1_y_xxx[k] - abx * g2_y_xx[k];

                        g_xy_xy[k] = g1_y_xxy[k] - abx * g2_y_xy[k];

                        g_xy_xz[k] = g1_y_xxz[k] - abx * g2_y_xz[k];

                        g_xy_yy[k] = g1_y_xyy[k] - abx * g2_y_yy[k];

                        g_xy_yz[k] = g1_y_xyz[k] - abx * g2_y_yz[k];

                        g_xy_zz[k] = g1_y_xzz[k] - abx * g2_y_zz[k];

                        g_xz_xx[k] = g1_z_xxx[k] - abx * g2_z_xx[k];

                        g_xz_xy[k] = g1_z_xxy[k] - abx * g2_z_xy[k];

                        g_xz_xz[k] = g1_z_xxz[k] - abx * g2_z_xz[k];

                        g_xz_yy[k] = g1_z_xyy[k] - abx * g2_z_yy[k];

                        g_xz_yz[k] = g1_z_xyz[k] - abx * g2_z_yz[k];

                        g_xz_zz[k] = g1_z_xzz[k] - abx * g2_z_zz[k];

                        // leading y component

                        g_yy_xx[k] = g1_y_xxy[k] - aby * g2_y_xx[k];

                        g_yy_xy[k] = g1_y_xyy[k] - aby * g2_y_xy[k];

                        g_yy_xz[k] = g1_y_xyz[k] - aby * g2_y_xz[k];

                        g_yy_yy[k] = g1_y_yyy[k] - aby * g2_y_yy[k];

                        g_yy_yz[k] = g1_y_yyz[k] - aby * g2_y_yz[k];

                        g_yy_zz[k] = g1_y_yzz[k] - aby * g2_y_zz[k];

                        g_yz_xx[k] = g1_z_xxy[k] - aby * g2_z_xx[k];

                        g_yz_xy[k] = g1_z_xyy[k] - aby * g2_z_xy[k];

                        g_yz_xz[k] = g1_z_xyz[k] - aby * g2_z_xz[k];

                        g_yz_yy[k] = g1_z_yyy[k] - aby * g2_z_yy[k];

                        g_yz_yz[k] = g1_z_yyz[k] - aby * g2_z_yz[k];

                        g_yz_zz[k] = g1_z_yzz[k] - aby * g2_z_zz[k];

                        // leading z component

                        g_zz_xx[k] = g1_z_xxz[k] - abz * g2_z_xx[k];

                        g_zz_xy[k] = g1_z_xyz[k] - abz * g2_z_xy[k];

                        g_zz_xz[k] = g1_z_xzz[k] - abz * g2_z_xz[k];

                        g_zz_yy[k] = g1_z_yyz[k] - abz * g2_z_yy[k];

                        g_zz_yz[k] = g1_z_yzz[k] - abz * g2_z_yz[k];

                        g_zz_zz[k] = g1_z_zzz[k] - abz * g2_z_zz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForDFXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 2) && (recPattern[i].second() == 3))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (23|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {2, 3, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 4, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 3, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (PF|g(r,r')|XX) integrals

                    auto g2_x_xxx = braBuffer.data(g2off + j);

                    auto g2_x_xxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_x_xxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_x_xyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_x_xyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_x_xzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_x_yyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_x_yyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_x_yzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_x_zzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_y_xxx = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_y_xxy = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_y_xxz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_y_xyy = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_y_xyz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_y_xzz = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_y_yyy = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_y_yyz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_y_yzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_y_zzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_z_xxx = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_z_xxy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_z_xxz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_z_xyy = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_z_xyz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_z_xzz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_z_yyy = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_z_yyz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_z_yzz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_z_zzz = braBuffer.data(g2off + 29 * kcomp + j);

                    // set up pointers to (PG|g(r,r')|XX) integrals

                    auto g1_x_xxxx = braBuffer.data(g1off + j);

                    auto g1_x_xxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_x_xxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_x_xxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_x_xxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_x_xxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_x_xyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_x_xyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_x_xyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_x_xzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_y_xxxx = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_y_xxxy = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_y_xxxz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_y_xxyy = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_y_xxyz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_y_xxzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_y_xyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_y_xyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_y_xyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_y_xzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_y_yyyy = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_y_yyyz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_y_yyzz = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_y_yzzz = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_z_xxxx = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_z_xxxy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_z_xxxz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_z_xxyy = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_z_xxyz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_z_xxzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_z_xyyy = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_z_xyyz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_z_xyzz = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_z_xzzz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_z_yyyy = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_z_yyyz = braBuffer.data(g1off + 41 * kcomp + j);

                    auto g1_z_yyzz = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_z_yzzz = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_z_zzzz = braBuffer.data(g1off + 44 * kcomp + j);

                    // set up pointers to (DF|g(r,r')|XX) integrals

                    auto g_xx_xxx = braBuffer.data(goff + j);

                    auto g_xx_xxy = braBuffer.data(goff + kcomp + j);

                    auto g_xx_xxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xx_xyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xx_xyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xx_xzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xx_yyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xx_yyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xx_yzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xx_zzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xy_xxx = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xy_xxy = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xy_xxz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xy_xyy = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xy_xyz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xy_xzz = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xy_yyy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xy_yyz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xy_yzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xy_zzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xz_xxx = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xz_xxy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xz_xxz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xz_xyy = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xz_xyz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xz_xzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xz_yyy = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xz_yyz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xz_yzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xz_zzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_yy_xxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_yy_xxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_yy_xxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_yy_xyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_yy_xyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_yy_xzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_yy_yyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_yy_yyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_yy_yzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_yy_zzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_yz_xxx = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_yz_xxy = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_yz_xxz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_yz_xyy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_yz_xyz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_yz_xzz = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_yz_yyy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_yz_yyz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_yz_yzz = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_yz_zzz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_zz_xxx = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_zz_xxy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_zz_xxz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_zz_xyy = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_zz_xyz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_zz_xzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_zz_yyy = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_zz_yyz = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_zz_yzz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_zz_zzz = braBuffer.data(goff + 59 * kcomp + j);

                    #pragma omp simd aligned(g2_x_xxx, g2_x_xxy, g2_x_xxz, g2_x_xyy,\
                                             g2_x_xyz, g2_x_xzz, g2_x_yyy, g2_x_yyz,\
                                             g2_x_yzz, g2_x_zzz, g2_y_xxx, g2_y_xxy,\
                                             g2_y_xxz, g2_y_xyy, g2_y_xyz, g2_y_xzz,\
                                             g2_y_yyy, g2_y_yyz, g2_y_yzz, g2_y_zzz,\
                                             g2_z_xxx, g2_z_xxy, g2_z_xxz, g2_z_xyy,\
                                             g2_z_xyz, g2_z_xzz, g2_z_yyy, g2_z_yyz,\
                                             g2_z_yzz, g2_z_zzz, g1_x_xxxx, g1_x_xxxy,\
                                             g1_x_xxxz, g1_x_xxyy, g1_x_xxyz, g1_x_xxzz,\
                                             g1_x_xyyy, g1_x_xyyz, g1_x_xyzz, g1_x_xzzz,\
                                             g1_y_xxxx, g1_y_xxxy, g1_y_xxxz,\
                                             g1_y_xxyy, g1_y_xxyz, g1_y_xxzz, g1_y_xyyy,\
                                             g1_y_xyyz, g1_y_xyzz, g1_y_xzzz, g1_y_yyyy,\
                                             g1_y_yyyz, g1_y_yyzz, g1_y_yzzz,\
                                             g1_z_xxxx, g1_z_xxxy, g1_z_xxxz, g1_z_xxyy,\
                                             g1_z_xxyz, g1_z_xxzz, g1_z_xyyy, g1_z_xyyz,\
                                             g1_z_xyzz, g1_z_xzzz, g1_z_yyyy, g1_z_yyyz,\
                                             g1_z_yyzz, g1_z_yzzz, g1_z_zzzz, g_xx_xxx,\
                                             g_xx_xxy, g_xx_xxz, g_xx_xyy, g_xx_xyz,\
                                             g_xx_xzz, g_xx_yyy, g_xx_yyz, g_xx_yzz,\
                                             g_xx_zzz, g_xy_xxx, g_xy_xxy, g_xy_xxz,\
                                             g_xy_xyy, g_xy_xyz, g_xy_xzz, g_xy_yyy,\
                                             g_xy_yyz, g_xy_yzz, g_xy_zzz, g_xz_xxx,\
                                             g_xz_xxy, g_xz_xxz, g_xz_xyy, g_xz_xyz,\
                                             g_xz_xzz, g_xz_yyy, g_xz_yyz, g_xz_yzz,\
                                             g_xz_zzz, g_yy_xxx, g_yy_xxy, g_yy_xxz,\
                                             g_yy_xyy, g_yy_xyz, g_yy_xzz, g_yy_yyy,\
                                             g_yy_yyz, g_yy_yzz, g_yy_zzz, g_yz_xxx,\
                                             g_yz_xxy, g_yz_xxz, g_yz_xyy, g_yz_xyz,\
                                             g_yz_xzz, g_yz_yyy, g_yz_yyz, g_yz_yzz,\
                                             g_yz_zzz, g_zz_xxx, g_zz_xxy, g_zz_xxz,\
                                             g_zz_xyy, g_zz_xyz, g_zz_xzz, g_zz_yyy,\
                                             g_zz_yyz, g_zz_yzz, g_zz_zzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xx_xxx[k] = g1_x_xxxx[k] - abx * g2_x_xxx[k];

                        g_xx_xxy[k] = g1_x_xxxy[k] - abx * g2_x_xxy[k];

                        g_xx_xxz[k] = g1_x_xxxz[k] - abx * g2_x_xxz[k];

                        g_xx_xyy[k] = g1_x_xxyy[k] - abx * g2_x_xyy[k];

                        g_xx_xyz[k] = g1_x_xxyz[k] - abx * g2_x_xyz[k];

                        g_xx_xzz[k] = g1_x_xxzz[k] - abx * g2_x_xzz[k];

                        g_xx_yyy[k] = g1_x_xyyy[k] - abx * g2_x_yyy[k];

                        g_xx_yyz[k] = g1_x_xyyz[k] - abx * g2_x_yyz[k];

                        g_xx_yzz[k] = g1_x_xyzz[k] - abx * g2_x_yzz[k];

                        g_xx_zzz[k] = g1_x_xzzz[k] - abx * g2_x_zzz[k];

                        g_xy_xxx[k] = g1_y_xxxx[k] - abx * g2_y_xxx[k];

                        g_xy_xxy[k] = g1_y_xxxy[k] - abx * g2_y_xxy[k];

                        g_xy_xxz[k] = g1_y_xxxz[k] - abx * g2_y_xxz[k];

                        g_xy_xyy[k] = g1_y_xxyy[k] - abx * g2_y_xyy[k];

                        g_xy_xyz[k] = g1_y_xxyz[k] - abx * g2_y_xyz[k];

                        g_xy_xzz[k] = g1_y_xxzz[k] - abx * g2_y_xzz[k];

                        g_xy_yyy[k] = g1_y_xyyy[k] - abx * g2_y_yyy[k];

                        g_xy_yyz[k] = g1_y_xyyz[k] - abx * g2_y_yyz[k];

                        g_xy_yzz[k] = g1_y_xyzz[k] - abx * g2_y_yzz[k];

                        g_xy_zzz[k] = g1_y_xzzz[k] - abx * g2_y_zzz[k];

                        g_xz_xxx[k] = g1_z_xxxx[k] - abx * g2_z_xxx[k];

                        g_xz_xxy[k] = g1_z_xxxy[k] - abx * g2_z_xxy[k];

                        g_xz_xxz[k] = g1_z_xxxz[k] - abx * g2_z_xxz[k];

                        g_xz_xyy[k] = g1_z_xxyy[k] - abx * g2_z_xyy[k];

                        g_xz_xyz[k] = g1_z_xxyz[k] - abx * g2_z_xyz[k];

                        g_xz_xzz[k] = g1_z_xxzz[k] - abx * g2_z_xzz[k];

                        g_xz_yyy[k] = g1_z_xyyy[k] - abx * g2_z_yyy[k];

                        g_xz_yyz[k] = g1_z_xyyz[k] - abx * g2_z_yyz[k];

                        g_xz_yzz[k] = g1_z_xyzz[k] - abx * g2_z_yzz[k];

                        g_xz_zzz[k] = g1_z_xzzz[k] - abx * g2_z_zzz[k];

                        // leading y component

                        g_yy_xxx[k] = g1_y_xxxy[k] - aby * g2_y_xxx[k];

                        g_yy_xxy[k] = g1_y_xxyy[k] - aby * g2_y_xxy[k];

                        g_yy_xxz[k] = g1_y_xxyz[k] - aby * g2_y_xxz[k];

                        g_yy_xyy[k] = g1_y_xyyy[k] - aby * g2_y_xyy[k];

                        g_yy_xyz[k] = g1_y_xyyz[k] - aby * g2_y_xyz[k];

                        g_yy_xzz[k] = g1_y_xyzz[k] - aby * g2_y_xzz[k];

                        g_yy_yyy[k] = g1_y_yyyy[k] - aby * g2_y_yyy[k];

                        g_yy_yyz[k] = g1_y_yyyz[k] - aby * g2_y_yyz[k];

                        g_yy_yzz[k] = g1_y_yyzz[k] - aby * g2_y_yzz[k];

                        g_yy_zzz[k] = g1_y_yzzz[k] - aby * g2_y_zzz[k];

                        g_yz_xxx[k] = g1_z_xxxy[k] - aby * g2_z_xxx[k];

                        g_yz_xxy[k] = g1_z_xxyy[k] - aby * g2_z_xxy[k];

                        g_yz_xxz[k] = g1_z_xxyz[k] - aby * g2_z_xxz[k];

                        g_yz_xyy[k] = g1_z_xyyy[k] - aby * g2_z_xyy[k];

                        g_yz_xyz[k] = g1_z_xyyz[k] - aby * g2_z_xyz[k];

                        g_yz_xzz[k] = g1_z_xyzz[k] - aby * g2_z_xzz[k];

                        g_yz_yyy[k] = g1_z_yyyy[k] - aby * g2_z_yyy[k];

                        g_yz_yyz[k] = g1_z_yyyz[k] - aby * g2_z_yyz[k];

                        g_yz_yzz[k] = g1_z_yyzz[k] - aby * g2_z_yzz[k];

                        g_yz_zzz[k] = g1_z_yzzz[k] - aby * g2_z_zzz[k];

                        // leading z component

                        g_zz_xxx[k] = g1_z_xxxz[k] - abz * g2_z_xxx[k];

                        g_zz_xxy[k] = g1_z_xxyz[k] - abz * g2_z_xxy[k];

                        g_zz_xxz[k] = g1_z_xxzz[k] - abz * g2_z_xxz[k];

                        g_zz_xyy[k] = g1_z_xyyz[k] - abz * g2_z_xyy[k];

                        g_zz_xyz[k] = g1_z_xyzz[k] - abz * g2_z_xyz[k];

                        g_zz_xzz[k] = g1_z_xzzz[k] - abz * g2_z_xzz[k];

                        g_zz_yyy[k] = g1_z_yyyz[k] - abz * g2_z_yyy[k];

                        g_zz_yyz[k] = g1_z_yyzz[k] - abz * g2_z_yyz[k];

                        g_zz_yzz[k] = g1_z_yzzz[k] - abz * g2_z_yzz[k];

                        g_zz_zzz[k] = g1_z_zzzz[k] - abz * g2_z_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForDGXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 2) && (recPattern[i].second() == 4))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (24|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {2, 4, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 5, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 4, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (PG|g(r,r')|XX) integrals

                    auto g2_x_xxxx = braBuffer.data(g2off + j);

                    auto g2_x_xxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_x_xxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_x_xxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_x_xxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_x_xxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_x_xyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_x_xyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_x_xyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_x_xzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_x_yyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_x_yyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_x_yyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_x_yzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_x_zzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_y_xxxx = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_y_xxxy = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_y_xxxz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_y_xxyy = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_y_xxyz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_y_xxzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_y_xyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_y_xyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_y_xyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_y_xzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_y_yyyy = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_y_yyyz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_y_yyzz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_y_yzzz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_y_zzzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_z_xxxx = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_z_xxxy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_z_xxxz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_z_xxyy = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_z_xxyz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_z_xxzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_z_xyyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_z_xyyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_z_xyzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_z_xzzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_z_yyyy = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_z_yyyz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_z_yyzz = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_z_yzzz = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_z_zzzz = braBuffer.data(g2off + 44 * kcomp + j);

                    // set up pointers to (PH|g(r,r')|XX) integrals

                    auto g1_x_xxxxx = braBuffer.data(g1off + j);

                    auto g1_x_xxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_x_xxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_x_xxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_x_xxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_x_xxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_x_xxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_x_xxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_x_xxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_x_xxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_x_xyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_x_xyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_x_xyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_x_xyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_x_xzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_y_xxxxx = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_y_xxxxy = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_y_xxxxz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_y_xxxyy = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_y_xxxyz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_y_xxxzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_y_xxyyy = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_y_xxyyz = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_y_xxyzz = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_y_xxzzz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_y_xyyyy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_y_xyyyz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_y_xyyzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_y_xyzzz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_y_xzzzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_y_yyyyy = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_y_yyyyz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_y_yyyzz = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_y_yyzzz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_y_yzzzz = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_z_xxxxx = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_z_xxxxy = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_z_xxxxz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_z_xxxyy = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_z_xxxyz = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_z_xxxzz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_z_xxyyy = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_z_xxyyz = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_z_xxyzz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_z_xxzzz = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_z_xyyyy = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_z_xyyyz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_z_xyyzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_z_xyzzz = braBuffer.data(g1off + 55 * kcomp + j);

                    auto g1_z_xzzzz = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_z_yyyyy = braBuffer.data(g1off + 57 * kcomp + j);

                    auto g1_z_yyyyz = braBuffer.data(g1off + 58 * kcomp + j);

                    auto g1_z_yyyzz = braBuffer.data(g1off + 59 * kcomp + j);

                    auto g1_z_yyzzz = braBuffer.data(g1off + 60 * kcomp + j);

                    auto g1_z_yzzzz = braBuffer.data(g1off + 61 * kcomp + j);

                    auto g1_z_zzzzz = braBuffer.data(g1off + 62 * kcomp + j);

                    // set up pointers to (DG|g(r,r')|XX) integrals

                    auto g_xx_xxxx = braBuffer.data(goff + j);

                    auto g_xx_xxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xx_xxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xx_xxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xx_xxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xx_xxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xx_xyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xx_xyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xx_xyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xx_xzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xx_yyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xx_yyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xx_yyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xx_yzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xx_zzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xy_xxxx = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xy_xxxy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xy_xxxz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xy_xxyy = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xy_xxyz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xy_xxzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xy_xyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xy_xyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xy_xyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xy_xzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xy_yyyy = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xy_yyyz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xy_yyzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xy_yzzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xy_zzzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xz_xxxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xz_xxxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xz_xxxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xz_xxyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xz_xxyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xz_xxzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xz_xyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xz_xyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xz_xyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xz_xzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xz_yyyy = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xz_yyyz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xz_yyzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xz_yzzz = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xz_zzzz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_yy_xxxx = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_yy_xxxy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_yy_xxxz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_yy_xxyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_yy_xxyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_yy_xxzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_yy_xyyy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_yy_xyyz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_yy_xyzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_yy_xzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_yy_yyyy = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_yy_yyyz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_yy_yyzz = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_yy_yzzz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_yy_zzzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_yz_xxxx = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_yz_xxxy = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_yz_xxxz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_yz_xxyy = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_yz_xxyz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_yz_xxzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_yz_xyyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_yz_xyyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_yz_xyzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_yz_xzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_yz_yyyy = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_yz_yyyz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_yz_yyzz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_yz_yzzz = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_yz_zzzz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_zz_xxxx = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_zz_xxxy = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_zz_xxxz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_zz_xxyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_zz_xxyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_zz_xxzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_zz_xyyy = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_zz_xyyz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_zz_xyzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_zz_xzzz = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_zz_yyyy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_zz_yyyz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_zz_yyzz = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_zz_yzzz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_zz_zzzz = braBuffer.data(goff + 89 * kcomp + j);

                    #pragma omp simd aligned(g2_x_xxxx, g2_x_xxxy, g2_x_xxxz, g2_x_xxyy,\
                                             g2_x_xxyz, g2_x_xxzz, g2_x_xyyy, g2_x_xyyz,\
                                             g2_x_xyzz, g2_x_xzzz, g2_x_yyyy, g2_x_yyyz,\
                                             g2_x_yyzz, g2_x_yzzz, g2_x_zzzz, g2_y_xxxx,\
                                             g2_y_xxxy, g2_y_xxxz, g2_y_xxyy, g2_y_xxyz,\
                                             g2_y_xxzz, g2_y_xyyy, g2_y_xyyz, g2_y_xyzz,\
                                             g2_y_xzzz, g2_y_yyyy, g2_y_yyyz, g2_y_yyzz,\
                                             g2_y_yzzz, g2_y_zzzz, g2_z_xxxx, g2_z_xxxy,\
                                             g2_z_xxxz, g2_z_xxyy, g2_z_xxyz, g2_z_xxzz,\
                                             g2_z_xyyy, g2_z_xyyz, g2_z_xyzz, g2_z_xzzz,\
                                             g2_z_yyyy, g2_z_yyyz, g2_z_yyzz, g2_z_yzzz,\
                                             g2_z_zzzz, g1_x_xxxxx, g1_x_xxxxy,\
                                             g1_x_xxxxz, g1_x_xxxyy, g1_x_xxxyz,\
                                             g1_x_xxxzz, g1_x_xxyyy, g1_x_xxyyz,\
                                             g1_x_xxyzz, g1_x_xxzzz, g1_x_xyyyy,\
                                             g1_x_xyyyz, g1_x_xyyzz, g1_x_xyzzz,\
                                             g1_x_xzzzz, g1_y_xxxxx, g1_y_xxxxy,\
                                             g1_y_xxxxz, g1_y_xxxyy, g1_y_xxxyz,\
                                             g1_y_xxxzz, g1_y_xxyyy, g1_y_xxyyz,\
                                             g1_y_xxyzz, g1_y_xxzzz, g1_y_xyyyy,\
                                             g1_y_xyyyz, g1_y_xyyzz, g1_y_xyzzz,\
                                             g1_y_xzzzz, g1_y_yyyyy, g1_y_yyyyz,\
                                             g1_y_yyyzz, g1_y_yyzzz, g1_y_yzzzz,\
                                             g1_z_xxxxx, g1_z_xxxxy,\
                                             g1_z_xxxxz, g1_z_xxxyy, g1_z_xxxyz,\
                                             g1_z_xxxzz, g1_z_xxyyy, g1_z_xxyyz,\
                                             g1_z_xxyzz, g1_z_xxzzz, g1_z_xyyyy,\
                                             g1_z_xyyyz, g1_z_xyyzz, g1_z_xyzzz,\
                                             g1_z_xzzzz, g1_z_yyyyy, g1_z_yyyyz,\
                                             g1_z_yyyzz, g1_z_yyzzz, g1_z_yzzzz,\
                                             g1_z_zzzzz, g_xx_xxxx, g_xx_xxxy,\
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
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xx_xxxx[k] = g1_x_xxxxx[k] - abx * g2_x_xxxx[k];

                        g_xx_xxxy[k] = g1_x_xxxxy[k] - abx * g2_x_xxxy[k];

                        g_xx_xxxz[k] = g1_x_xxxxz[k] - abx * g2_x_xxxz[k];

                        g_xx_xxyy[k] = g1_x_xxxyy[k] - abx * g2_x_xxyy[k];

                        g_xx_xxyz[k] = g1_x_xxxyz[k] - abx * g2_x_xxyz[k];

                        g_xx_xxzz[k] = g1_x_xxxzz[k] - abx * g2_x_xxzz[k];

                        g_xx_xyyy[k] = g1_x_xxyyy[k] - abx * g2_x_xyyy[k];

                        g_xx_xyyz[k] = g1_x_xxyyz[k] - abx * g2_x_xyyz[k];

                        g_xx_xyzz[k] = g1_x_xxyzz[k] - abx * g2_x_xyzz[k];

                        g_xx_xzzz[k] = g1_x_xxzzz[k] - abx * g2_x_xzzz[k];

                        g_xx_yyyy[k] = g1_x_xyyyy[k] - abx * g2_x_yyyy[k];

                        g_xx_yyyz[k] = g1_x_xyyyz[k] - abx * g2_x_yyyz[k];

                        g_xx_yyzz[k] = g1_x_xyyzz[k] - abx * g2_x_yyzz[k];

                        g_xx_yzzz[k] = g1_x_xyzzz[k] - abx * g2_x_yzzz[k];

                        g_xx_zzzz[k] = g1_x_xzzzz[k] - abx * g2_x_zzzz[k];

                        g_xy_xxxx[k] = g1_y_xxxxx[k] - abx * g2_y_xxxx[k];

                        g_xy_xxxy[k] = g1_y_xxxxy[k] - abx * g2_y_xxxy[k];

                        g_xy_xxxz[k] = g1_y_xxxxz[k] - abx * g2_y_xxxz[k];

                        g_xy_xxyy[k] = g1_y_xxxyy[k] - abx * g2_y_xxyy[k];

                        g_xy_xxyz[k] = g1_y_xxxyz[k] - abx * g2_y_xxyz[k];

                        g_xy_xxzz[k] = g1_y_xxxzz[k] - abx * g2_y_xxzz[k];

                        g_xy_xyyy[k] = g1_y_xxyyy[k] - abx * g2_y_xyyy[k];

                        g_xy_xyyz[k] = g1_y_xxyyz[k] - abx * g2_y_xyyz[k];

                        g_xy_xyzz[k] = g1_y_xxyzz[k] - abx * g2_y_xyzz[k];

                        g_xy_xzzz[k] = g1_y_xxzzz[k] - abx * g2_y_xzzz[k];

                        g_xy_yyyy[k] = g1_y_xyyyy[k] - abx * g2_y_yyyy[k];

                        g_xy_yyyz[k] = g1_y_xyyyz[k] - abx * g2_y_yyyz[k];

                        g_xy_yyzz[k] = g1_y_xyyzz[k] - abx * g2_y_yyzz[k];

                        g_xy_yzzz[k] = g1_y_xyzzz[k] - abx * g2_y_yzzz[k];

                        g_xy_zzzz[k] = g1_y_xzzzz[k] - abx * g2_y_zzzz[k];

                        g_xz_xxxx[k] = g1_z_xxxxx[k] - abx * g2_z_xxxx[k];

                        g_xz_xxxy[k] = g1_z_xxxxy[k] - abx * g2_z_xxxy[k];

                        g_xz_xxxz[k] = g1_z_xxxxz[k] - abx * g2_z_xxxz[k];

                        g_xz_xxyy[k] = g1_z_xxxyy[k] - abx * g2_z_xxyy[k];

                        g_xz_xxyz[k] = g1_z_xxxyz[k] - abx * g2_z_xxyz[k];

                        g_xz_xxzz[k] = g1_z_xxxzz[k] - abx * g2_z_xxzz[k];

                        g_xz_xyyy[k] = g1_z_xxyyy[k] - abx * g2_z_xyyy[k];

                        g_xz_xyyz[k] = g1_z_xxyyz[k] - abx * g2_z_xyyz[k];

                        g_xz_xyzz[k] = g1_z_xxyzz[k] - abx * g2_z_xyzz[k];

                        g_xz_xzzz[k] = g1_z_xxzzz[k] - abx * g2_z_xzzz[k];

                        g_xz_yyyy[k] = g1_z_xyyyy[k] - abx * g2_z_yyyy[k];

                        g_xz_yyyz[k] = g1_z_xyyyz[k] - abx * g2_z_yyyz[k];

                        g_xz_yyzz[k] = g1_z_xyyzz[k] - abx * g2_z_yyzz[k];

                        g_xz_yzzz[k] = g1_z_xyzzz[k] - abx * g2_z_yzzz[k];

                        g_xz_zzzz[k] = g1_z_xzzzz[k] - abx * g2_z_zzzz[k];

                        // leading y component

                        g_yy_xxxx[k] = g1_y_xxxxy[k] - aby * g2_y_xxxx[k];

                        g_yy_xxxy[k] = g1_y_xxxyy[k] - aby * g2_y_xxxy[k];

                        g_yy_xxxz[k] = g1_y_xxxyz[k] - aby * g2_y_xxxz[k];

                        g_yy_xxyy[k] = g1_y_xxyyy[k] - aby * g2_y_xxyy[k];

                        g_yy_xxyz[k] = g1_y_xxyyz[k] - aby * g2_y_xxyz[k];

                        g_yy_xxzz[k] = g1_y_xxyzz[k] - aby * g2_y_xxzz[k];

                        g_yy_xyyy[k] = g1_y_xyyyy[k] - aby * g2_y_xyyy[k];

                        g_yy_xyyz[k] = g1_y_xyyyz[k] - aby * g2_y_xyyz[k];

                        g_yy_xyzz[k] = g1_y_xyyzz[k] - aby * g2_y_xyzz[k];

                        g_yy_xzzz[k] = g1_y_xyzzz[k] - aby * g2_y_xzzz[k];

                        g_yy_yyyy[k] = g1_y_yyyyy[k] - aby * g2_y_yyyy[k];

                        g_yy_yyyz[k] = g1_y_yyyyz[k] - aby * g2_y_yyyz[k];

                        g_yy_yyzz[k] = g1_y_yyyzz[k] - aby * g2_y_yyzz[k];

                        g_yy_yzzz[k] = g1_y_yyzzz[k] - aby * g2_y_yzzz[k];

                        g_yy_zzzz[k] = g1_y_yzzzz[k] - aby * g2_y_zzzz[k];

                        g_yz_xxxx[k] = g1_z_xxxxy[k] - aby * g2_z_xxxx[k];

                        g_yz_xxxy[k] = g1_z_xxxyy[k] - aby * g2_z_xxxy[k];

                        g_yz_xxxz[k] = g1_z_xxxyz[k] - aby * g2_z_xxxz[k];

                        g_yz_xxyy[k] = g1_z_xxyyy[k] - aby * g2_z_xxyy[k];

                        g_yz_xxyz[k] = g1_z_xxyyz[k] - aby * g2_z_xxyz[k];

                        g_yz_xxzz[k] = g1_z_xxyzz[k] - aby * g2_z_xxzz[k];

                        g_yz_xyyy[k] = g1_z_xyyyy[k] - aby * g2_z_xyyy[k];

                        g_yz_xyyz[k] = g1_z_xyyyz[k] - aby * g2_z_xyyz[k];

                        g_yz_xyzz[k] = g1_z_xyyzz[k] - aby * g2_z_xyzz[k];

                        g_yz_xzzz[k] = g1_z_xyzzz[k] - aby * g2_z_xzzz[k];

                        g_yz_yyyy[k] = g1_z_yyyyy[k] - aby * g2_z_yyyy[k];

                        g_yz_yyyz[k] = g1_z_yyyyz[k] - aby * g2_z_yyyz[k];

                        g_yz_yyzz[k] = g1_z_yyyzz[k] - aby * g2_z_yyzz[k];

                        g_yz_yzzz[k] = g1_z_yyzzz[k] - aby * g2_z_yzzz[k];

                        g_yz_zzzz[k] = g1_z_yzzzz[k] - aby * g2_z_zzzz[k];

                        // leading z component

                        g_zz_xxxx[k] = g1_z_xxxxz[k] - abz * g2_z_xxxx[k];

                        g_zz_xxxy[k] = g1_z_xxxyz[k] - abz * g2_z_xxxy[k];

                        g_zz_xxxz[k] = g1_z_xxxzz[k] - abz * g2_z_xxxz[k];

                        g_zz_xxyy[k] = g1_z_xxyyz[k] - abz * g2_z_xxyy[k];

                        g_zz_xxyz[k] = g1_z_xxyzz[k] - abz * g2_z_xxyz[k];

                        g_zz_xxzz[k] = g1_z_xxzzz[k] - abz * g2_z_xxzz[k];

                        g_zz_xyyy[k] = g1_z_xyyyz[k] - abz * g2_z_xyyy[k];

                        g_zz_xyyz[k] = g1_z_xyyzz[k] - abz * g2_z_xyyz[k];

                        g_zz_xyzz[k] = g1_z_xyzzz[k] - abz * g2_z_xyzz[k];

                        g_zz_xzzz[k] = g1_z_xzzzz[k] - abz * g2_z_xzzz[k];

                        g_zz_yyyy[k] = g1_z_yyyyz[k] - abz * g2_z_yyyy[k];

                        g_zz_yyyz[k] = g1_z_yyyzz[k] - abz * g2_z_yyyz[k];

                        g_zz_yyzz[k] = g1_z_yyzzz[k] - abz * g2_z_yyzz[k];

                        g_zz_yzzz[k] = g1_z_yzzzz[k] - abz * g2_z_yzzz[k];

                        g_zz_zzzz[k] = g1_z_zzzzz[k] - abz * g2_z_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForDHXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 2) && (recPattern[i].second() == 5))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (25|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {2, 5, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 6, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 5, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (PH|g(r,r')|XX) integrals

                    auto g2_x_xxxxx = braBuffer.data(g2off + j);

                    auto g2_x_xxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_x_xxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_x_xxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_x_xxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_x_xxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_x_xxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_x_xxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_x_xxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_x_xxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_x_xyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_x_xyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_x_xyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_x_xyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_x_xzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_x_yyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_x_yyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_x_yyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_x_yyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_x_yzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_x_zzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_y_xxxxx = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_y_xxxxy = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_y_xxxxz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_y_xxxyy = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_y_xxxyz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_y_xxxzz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_y_xxyyy = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_y_xxyyz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_y_xxyzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_y_xxzzz = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_y_xyyyy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_y_xyyyz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_y_xyyzz = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_y_xyzzz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_y_xzzzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_y_yyyyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_y_yyyyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_y_yyyzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_y_yyzzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_y_yzzzz = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_y_zzzzz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_z_xxxxx = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_z_xxxxy = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_z_xxxxz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_z_xxxyy = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_z_xxxyz = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_z_xxxzz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_z_xxyyy = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_z_xxyyz = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_z_xxyzz = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_z_xxzzz = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_z_xyyyy = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_z_xyyyz = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_z_xyyzz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_z_xyzzz = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_z_xzzzz = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_z_yyyyy = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_z_yyyyz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_z_yyyzz = braBuffer.data(g2off + 59 * kcomp + j);

                    auto g2_z_yyzzz = braBuffer.data(g2off + 60 * kcomp + j);

                    auto g2_z_yzzzz = braBuffer.data(g2off + 61 * kcomp + j);

                    auto g2_z_zzzzz = braBuffer.data(g2off + 62 * kcomp + j);

                    // set up pointers to (PI|g(r,r')|XX) integrals

                    auto g1_x_xxxxxx = braBuffer.data(g1off + j);

                    auto g1_x_xxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_x_xxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_x_xxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_x_xxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_x_xxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_x_xxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_x_xxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_x_xxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_x_xxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_x_xxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_x_xxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_x_xxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_x_xxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_x_xxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_x_xyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_x_xyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_x_xyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_x_xyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_x_xyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_x_xzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_y_xxxxxx = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_y_xxxxxy = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_y_xxxxxz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_y_xxxxyy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_y_xxxxyz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_y_xxxxzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_y_xxxyyy = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_y_xxxyyz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_y_xxxyzz = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_y_xxxzzz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_y_xxyyyy = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_y_xxyyyz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_y_xxyyzz = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_y_xxyzzz = braBuffer.data(g1off + 41 * kcomp + j);

                    auto g1_y_xxzzzz = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_y_xyyyyy = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_y_xyyyyz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_y_xyyyzz = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_y_xyyzzz = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_y_xyzzzz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_y_xzzzzz = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_y_yyyyyy = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_y_yyyyyz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_y_yyyyzz = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_y_yyyzzz = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_y_yyzzzz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_y_yzzzzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_z_xxxxxx = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_z_xxxxxy = braBuffer.data(g1off + 57 * kcomp + j);

                    auto g1_z_xxxxxz = braBuffer.data(g1off + 58 * kcomp + j);

                    auto g1_z_xxxxyy = braBuffer.data(g1off + 59 * kcomp + j);

                    auto g1_z_xxxxyz = braBuffer.data(g1off + 60 * kcomp + j);

                    auto g1_z_xxxxzz = braBuffer.data(g1off + 61 * kcomp + j);

                    auto g1_z_xxxyyy = braBuffer.data(g1off + 62 * kcomp + j);

                    auto g1_z_xxxyyz = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_z_xxxyzz = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_z_xxxzzz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_z_xxyyyy = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_z_xxyyyz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_z_xxyyzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_z_xxyzzz = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_z_xxzzzz = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_z_xyyyyy = braBuffer.data(g1off + 71 * kcomp + j);

                    auto g1_z_xyyyyz = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_z_xyyyzz = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_z_xyyzzz = braBuffer.data(g1off + 74 * kcomp + j);

                    auto g1_z_xyzzzz = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_z_xzzzzz = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_z_yyyyyy = braBuffer.data(g1off + 77 * kcomp + j);

                    auto g1_z_yyyyyz = braBuffer.data(g1off + 78 * kcomp + j);

                    auto g1_z_yyyyzz = braBuffer.data(g1off + 79 * kcomp + j);

                    auto g1_z_yyyzzz = braBuffer.data(g1off + 80 * kcomp + j);

                    auto g1_z_yyzzzz = braBuffer.data(g1off + 81 * kcomp + j);

                    auto g1_z_yzzzzz = braBuffer.data(g1off + 82 * kcomp + j);

                    auto g1_z_zzzzzz = braBuffer.data(g1off + 83 * kcomp + j);

                    // set up pointers to (DH|g(r,r')|XX) integrals

                    auto g_xx_xxxxx = braBuffer.data(goff + j);

                    auto g_xx_xxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xx_xxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xx_xxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xx_xxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xx_xxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xx_xxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xx_xxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xx_xxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xx_xxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xx_xyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xx_xyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xx_xyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xx_xyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xx_xzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xx_yyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xx_yyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xx_yyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xx_yyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xx_yzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xx_zzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xy_xxxxx = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xy_xxxxy = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xy_xxxxz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xy_xxxyy = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xy_xxxyz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xy_xxxzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xy_xxyyy = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xy_xxyyz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xy_xxyzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xy_xxzzz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xy_xyyyy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xy_xyyyz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xy_xyyzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xy_xyzzz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xy_xzzzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xy_yyyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xy_yyyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xy_yyyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xy_yyzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xy_yzzzz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xy_zzzzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xz_xxxxx = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xz_xxxxy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xz_xxxxz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xz_xxxyy = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xz_xxxyz = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xz_xxxzz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xz_xxyyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xz_xxyyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xz_xxyzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xz_xxzzz = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xz_xyyyy = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xz_xyyyz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xz_xyyzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xz_xyzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xz_xzzzz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xz_yyyyy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xz_yyyyz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xz_yyyzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_xz_yyzzz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_xz_yzzzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_xz_zzzzz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_yy_xxxxx = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_yy_xxxxy = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_yy_xxxxz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_yy_xxxyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_yy_xxxyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_yy_xxxzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_yy_xxyyy = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_yy_xxyyz = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_yy_xxyzz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_yy_xxzzz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_yy_xyyyy = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_yy_xyyyz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_yy_xyyzz = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_yy_xyzzz = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_yy_xzzzz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_yy_yyyyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_yy_yyyyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_yy_yyyzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_yy_yyzzz = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_yy_yzzzz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_yy_zzzzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_yz_xxxxx = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_yz_xxxxy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_yz_xxxxz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_yz_xxxyy = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_yz_xxxyz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_yz_xxxzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_yz_xxyyy = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_yz_xxyyz = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_yz_xxyzz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_yz_xxzzz = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_yz_xyyyy = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_yz_xyyyz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_yz_xyyzz = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_yz_xyzzz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_yz_xzzzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_yz_yyyyy = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_yz_yyyyz = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_yz_yyyzz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_yz_yyzzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_yz_yzzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_yz_zzzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_zz_xxxxx = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_zz_xxxxy = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_zz_xxxxz = braBuffer.data(goff + 107 * kcomp + j);

                    auto g_zz_xxxyy = braBuffer.data(goff + 108 * kcomp + j);

                    auto g_zz_xxxyz = braBuffer.data(goff + 109 * kcomp + j);

                    auto g_zz_xxxzz = braBuffer.data(goff + 110 * kcomp + j);

                    auto g_zz_xxyyy = braBuffer.data(goff + 111 * kcomp + j);

                    auto g_zz_xxyyz = braBuffer.data(goff + 112 * kcomp + j);

                    auto g_zz_xxyzz = braBuffer.data(goff + 113 * kcomp + j);

                    auto g_zz_xxzzz = braBuffer.data(goff + 114 * kcomp + j);

                    auto g_zz_xyyyy = braBuffer.data(goff + 115 * kcomp + j);

                    auto g_zz_xyyyz = braBuffer.data(goff + 116 * kcomp + j);

                    auto g_zz_xyyzz = braBuffer.data(goff + 117 * kcomp + j);

                    auto g_zz_xyzzz = braBuffer.data(goff + 118 * kcomp + j);

                    auto g_zz_xzzzz = braBuffer.data(goff + 119 * kcomp + j);

                    auto g_zz_yyyyy = braBuffer.data(goff + 120 * kcomp + j);

                    auto g_zz_yyyyz = braBuffer.data(goff + 121 * kcomp + j);

                    auto g_zz_yyyzz = braBuffer.data(goff + 122 * kcomp + j);

                    auto g_zz_yyzzz = braBuffer.data(goff + 123 * kcomp + j);

                    auto g_zz_yzzzz = braBuffer.data(goff + 124 * kcomp + j);

                    auto g_zz_zzzzz = braBuffer.data(goff + 125 * kcomp + j);

                    #pragma omp simd aligned(g2_x_xxxxx, g2_x_xxxxy, g2_x_xxxxz,\
                                             g2_x_xxxyy, g2_x_xxxyz, g2_x_xxxzz,\
                                             g2_x_xxyyy, g2_x_xxyyz, g2_x_xxyzz,\
                                             g2_x_xxzzz, g2_x_xyyyy, g2_x_xyyyz,\
                                             g2_x_xyyzz, g2_x_xyzzz, g2_x_xzzzz,\
                                             g2_x_yyyyy, g2_x_yyyyz, g2_x_yyyzz,\
                                             g2_x_yyzzz, g2_x_yzzzz, g2_x_zzzzz,\
                                             g2_y_xxxxx, g2_y_xxxxy, g2_y_xxxxz,\
                                             g2_y_xxxyy, g2_y_xxxyz, g2_y_xxxzz,\
                                             g2_y_xxyyy, g2_y_xxyyz, g2_y_xxyzz,\
                                             g2_y_xxzzz, g2_y_xyyyy, g2_y_xyyyz,\
                                             g2_y_xyyzz, g2_y_xyzzz, g2_y_xzzzz,\
                                             g2_y_yyyyy, g2_y_yyyyz, g2_y_yyyzz,\
                                             g2_y_yyzzz, g2_y_yzzzz, g2_y_zzzzz,\
                                             g2_z_xxxxx, g2_z_xxxxy, g2_z_xxxxz,\
                                             g2_z_xxxyy, g2_z_xxxyz, g2_z_xxxzz,\
                                             g2_z_xxyyy, g2_z_xxyyz, g2_z_xxyzz,\
                                             g2_z_xxzzz, g2_z_xyyyy, g2_z_xyyyz,\
                                             g2_z_xyyzz, g2_z_xyzzz, g2_z_xzzzz,\
                                             g2_z_yyyyy, g2_z_yyyyz, g2_z_yyyzz,\
                                             g2_z_yyzzz, g2_z_yzzzz, g2_z_zzzzz,\
                                             g1_x_xxxxxx, g1_x_xxxxxy, g1_x_xxxxxz,\
                                             g1_x_xxxxyy, g1_x_xxxxyz, g1_x_xxxxzz,\
                                             g1_x_xxxyyy, g1_x_xxxyyz, g1_x_xxxyzz,\
                                             g1_x_xxxzzz, g1_x_xxyyyy, g1_x_xxyyyz,\
                                             g1_x_xxyyzz, g1_x_xxyzzz, g1_x_xxzzzz,\
                                             g1_x_xyyyyy, g1_x_xyyyyz, g1_x_xyyyzz,\
                                             g1_x_xyyzzz, g1_x_xyzzzz, g1_x_xzzzzz,\
                                             g1_y_xxxxxx, g1_y_xxxxxy,\
                                             g1_y_xxxxxz, g1_y_xxxxyy, g1_y_xxxxyz,\
                                             g1_y_xxxxzz, g1_y_xxxyyy, g1_y_xxxyyz,\
                                             g1_y_xxxyzz, g1_y_xxxzzz, g1_y_xxyyyy,\
                                             g1_y_xxyyyz, g1_y_xxyyzz, g1_y_xxyzzz,\
                                             g1_y_xxzzzz, g1_y_xyyyyy, g1_y_xyyyyz,\
                                             g1_y_xyyyzz, g1_y_xyyzzz, g1_y_xyzzzz,\
                                             g1_y_xzzzzz, g1_y_yyyyyy, g1_y_yyyyyz,\
                                             g1_y_yyyyzz, g1_y_yyyzzz, g1_y_yyzzzz,\
                                             g1_y_yzzzzz, g1_z_xxxxxx,\
                                             g1_z_xxxxxy, g1_z_xxxxxz, g1_z_xxxxyy,\
                                             g1_z_xxxxyz, g1_z_xxxxzz, g1_z_xxxyyy,\
                                             g1_z_xxxyyz, g1_z_xxxyzz, g1_z_xxxzzz,\
                                             g1_z_xxyyyy, g1_z_xxyyyz, g1_z_xxyyzz,\
                                             g1_z_xxyzzz, g1_z_xxzzzz, g1_z_xyyyyy,\
                                             g1_z_xyyyyz, g1_z_xyyyzz, g1_z_xyyzzz,\
                                             g1_z_xyzzzz, g1_z_xzzzzz, g1_z_yyyyyy,\
                                             g1_z_yyyyyz, g1_z_yyyyzz, g1_z_yyyzzz,\
                                             g1_z_yyzzzz, g1_z_yzzzzz, g1_z_zzzzzz,\
                                             g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz,\
                                             g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxzz,\
                                             g_xx_xxyyy, g_xx_xxyyz, g_xx_xxyzz,\
                                             g_xx_xxzzz, g_xx_xyyyy, g_xx_xyyyz,\
                                             g_xx_xyyzz, g_xx_xyzzz, g_xx_xzzzz,\
                                             g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz,\
                                             g_xx_yyzzz, g_xx_yzzzz, g_xx_zzzzz,\
                                             g_xy_xxxxx, g_xy_xxxxy, g_xy_xxxxz,\
                                             g_xy_xxxyy, g_xy_xxxyz, g_xy_xxxzz,\
                                             g_xy_xxyyy, g_xy_xxyyz, g_xy_xxyzz,\
                                             g_xy_xxzzz, g_xy_xyyyy, g_xy_xyyyz,\
                                             g_xy_xyyzz, g_xy_xyzzz, g_xy_xzzzz,\
                                             g_xy_yyyyy, g_xy_yyyyz, g_xy_yyyzz,\
                                             g_xy_yyzzz, g_xy_yzzzz, g_xy_zzzzz,\
                                             g_xz_xxxxx, g_xz_xxxxy, g_xz_xxxxz,\
                                             g_xz_xxxyy, g_xz_xxxyz, g_xz_xxxzz,\
                                             g_xz_xxyyy, g_xz_xxyyz, g_xz_xxyzz,\
                                             g_xz_xxzzz, g_xz_xyyyy, g_xz_xyyyz,\
                                             g_xz_xyyzz, g_xz_xyzzz, g_xz_xzzzz,\
                                             g_xz_yyyyy, g_xz_yyyyz, g_xz_yyyzz,\
                                             g_xz_yyzzz, g_xz_yzzzz, g_xz_zzzzz,\
                                             g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz,\
                                             g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxzz,\
                                             g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyzz,\
                                             g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyz,\
                                             g_yy_xyyzz, g_yy_xyzzz, g_yy_xzzzz,\
                                             g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyzz,\
                                             g_yy_yyzzz, g_yy_yzzzz, g_yy_zzzzz,\
                                             g_yz_xxxxx, g_yz_xxxxy, g_yz_xxxxz,\
                                             g_yz_xxxyy, g_yz_xxxyz, g_yz_xxxzz,\
                                             g_yz_xxyyy, g_yz_xxyyz, g_yz_xxyzz,\
                                             g_yz_xxzzz, g_yz_xyyyy, g_yz_xyyyz,\
                                             g_yz_xyyzz, g_yz_xyzzz, g_yz_xzzzz,\
                                             g_yz_yyyyy, g_yz_yyyyz, g_yz_yyyzz,\
                                             g_yz_yyzzz, g_yz_yzzzz, g_yz_zzzzz,\
                                             g_zz_xxxxx, g_zz_xxxxy, g_zz_xxxxz,\
                                             g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxzz,\
                                             g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyzz,\
                                             g_zz_xxzzz, g_zz_xyyyy, g_zz_xyyyz,\
                                             g_zz_xyyzz, g_zz_xyzzz, g_zz_xzzzz,\
                                             g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz,\
                                             g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xx_xxxxx[k] = g1_x_xxxxxx[k] - abx * g2_x_xxxxx[k];

                        g_xx_xxxxy[k] = g1_x_xxxxxy[k] - abx * g2_x_xxxxy[k];

                        g_xx_xxxxz[k] = g1_x_xxxxxz[k] - abx * g2_x_xxxxz[k];

                        g_xx_xxxyy[k] = g1_x_xxxxyy[k] - abx * g2_x_xxxyy[k];

                        g_xx_xxxyz[k] = g1_x_xxxxyz[k] - abx * g2_x_xxxyz[k];

                        g_xx_xxxzz[k] = g1_x_xxxxzz[k] - abx * g2_x_xxxzz[k];

                        g_xx_xxyyy[k] = g1_x_xxxyyy[k] - abx * g2_x_xxyyy[k];

                        g_xx_xxyyz[k] = g1_x_xxxyyz[k] - abx * g2_x_xxyyz[k];

                        g_xx_xxyzz[k] = g1_x_xxxyzz[k] - abx * g2_x_xxyzz[k];

                        g_xx_xxzzz[k] = g1_x_xxxzzz[k] - abx * g2_x_xxzzz[k];

                        g_xx_xyyyy[k] = g1_x_xxyyyy[k] - abx * g2_x_xyyyy[k];

                        g_xx_xyyyz[k] = g1_x_xxyyyz[k] - abx * g2_x_xyyyz[k];

                        g_xx_xyyzz[k] = g1_x_xxyyzz[k] - abx * g2_x_xyyzz[k];

                        g_xx_xyzzz[k] = g1_x_xxyzzz[k] - abx * g2_x_xyzzz[k];

                        g_xx_xzzzz[k] = g1_x_xxzzzz[k] - abx * g2_x_xzzzz[k];

                        g_xx_yyyyy[k] = g1_x_xyyyyy[k] - abx * g2_x_yyyyy[k];

                        g_xx_yyyyz[k] = g1_x_xyyyyz[k] - abx * g2_x_yyyyz[k];

                        g_xx_yyyzz[k] = g1_x_xyyyzz[k] - abx * g2_x_yyyzz[k];

                        g_xx_yyzzz[k] = g1_x_xyyzzz[k] - abx * g2_x_yyzzz[k];

                        g_xx_yzzzz[k] = g1_x_xyzzzz[k] - abx * g2_x_yzzzz[k];

                        g_xx_zzzzz[k] = g1_x_xzzzzz[k] - abx * g2_x_zzzzz[k];

                        g_xy_xxxxx[k] = g1_y_xxxxxx[k] - abx * g2_y_xxxxx[k];

                        g_xy_xxxxy[k] = g1_y_xxxxxy[k] - abx * g2_y_xxxxy[k];

                        g_xy_xxxxz[k] = g1_y_xxxxxz[k] - abx * g2_y_xxxxz[k];

                        g_xy_xxxyy[k] = g1_y_xxxxyy[k] - abx * g2_y_xxxyy[k];

                        g_xy_xxxyz[k] = g1_y_xxxxyz[k] - abx * g2_y_xxxyz[k];

                        g_xy_xxxzz[k] = g1_y_xxxxzz[k] - abx * g2_y_xxxzz[k];

                        g_xy_xxyyy[k] = g1_y_xxxyyy[k] - abx * g2_y_xxyyy[k];

                        g_xy_xxyyz[k] = g1_y_xxxyyz[k] - abx * g2_y_xxyyz[k];

                        g_xy_xxyzz[k] = g1_y_xxxyzz[k] - abx * g2_y_xxyzz[k];

                        g_xy_xxzzz[k] = g1_y_xxxzzz[k] - abx * g2_y_xxzzz[k];

                        g_xy_xyyyy[k] = g1_y_xxyyyy[k] - abx * g2_y_xyyyy[k];

                        g_xy_xyyyz[k] = g1_y_xxyyyz[k] - abx * g2_y_xyyyz[k];

                        g_xy_xyyzz[k] = g1_y_xxyyzz[k] - abx * g2_y_xyyzz[k];

                        g_xy_xyzzz[k] = g1_y_xxyzzz[k] - abx * g2_y_xyzzz[k];

                        g_xy_xzzzz[k] = g1_y_xxzzzz[k] - abx * g2_y_xzzzz[k];

                        g_xy_yyyyy[k] = g1_y_xyyyyy[k] - abx * g2_y_yyyyy[k];

                        g_xy_yyyyz[k] = g1_y_xyyyyz[k] - abx * g2_y_yyyyz[k];

                        g_xy_yyyzz[k] = g1_y_xyyyzz[k] - abx * g2_y_yyyzz[k];

                        g_xy_yyzzz[k] = g1_y_xyyzzz[k] - abx * g2_y_yyzzz[k];

                        g_xy_yzzzz[k] = g1_y_xyzzzz[k] - abx * g2_y_yzzzz[k];

                        g_xy_zzzzz[k] = g1_y_xzzzzz[k] - abx * g2_y_zzzzz[k];

                        g_xz_xxxxx[k] = g1_z_xxxxxx[k] - abx * g2_z_xxxxx[k];

                        g_xz_xxxxy[k] = g1_z_xxxxxy[k] - abx * g2_z_xxxxy[k];

                        g_xz_xxxxz[k] = g1_z_xxxxxz[k] - abx * g2_z_xxxxz[k];

                        g_xz_xxxyy[k] = g1_z_xxxxyy[k] - abx * g2_z_xxxyy[k];

                        g_xz_xxxyz[k] = g1_z_xxxxyz[k] - abx * g2_z_xxxyz[k];

                        g_xz_xxxzz[k] = g1_z_xxxxzz[k] - abx * g2_z_xxxzz[k];

                        g_xz_xxyyy[k] = g1_z_xxxyyy[k] - abx * g2_z_xxyyy[k];

                        g_xz_xxyyz[k] = g1_z_xxxyyz[k] - abx * g2_z_xxyyz[k];

                        g_xz_xxyzz[k] = g1_z_xxxyzz[k] - abx * g2_z_xxyzz[k];

                        g_xz_xxzzz[k] = g1_z_xxxzzz[k] - abx * g2_z_xxzzz[k];

                        g_xz_xyyyy[k] = g1_z_xxyyyy[k] - abx * g2_z_xyyyy[k];

                        g_xz_xyyyz[k] = g1_z_xxyyyz[k] - abx * g2_z_xyyyz[k];

                        g_xz_xyyzz[k] = g1_z_xxyyzz[k] - abx * g2_z_xyyzz[k];

                        g_xz_xyzzz[k] = g1_z_xxyzzz[k] - abx * g2_z_xyzzz[k];

                        g_xz_xzzzz[k] = g1_z_xxzzzz[k] - abx * g2_z_xzzzz[k];

                        g_xz_yyyyy[k] = g1_z_xyyyyy[k] - abx * g2_z_yyyyy[k];

                        g_xz_yyyyz[k] = g1_z_xyyyyz[k] - abx * g2_z_yyyyz[k];

                        g_xz_yyyzz[k] = g1_z_xyyyzz[k] - abx * g2_z_yyyzz[k];

                        g_xz_yyzzz[k] = g1_z_xyyzzz[k] - abx * g2_z_yyzzz[k];

                        g_xz_yzzzz[k] = g1_z_xyzzzz[k] - abx * g2_z_yzzzz[k];

                        g_xz_zzzzz[k] = g1_z_xzzzzz[k] - abx * g2_z_zzzzz[k];

                        // leading y component

                        g_yy_xxxxx[k] = g1_y_xxxxxy[k] - aby * g2_y_xxxxx[k];

                        g_yy_xxxxy[k] = g1_y_xxxxyy[k] - aby * g2_y_xxxxy[k];

                        g_yy_xxxxz[k] = g1_y_xxxxyz[k] - aby * g2_y_xxxxz[k];

                        g_yy_xxxyy[k] = g1_y_xxxyyy[k] - aby * g2_y_xxxyy[k];

                        g_yy_xxxyz[k] = g1_y_xxxyyz[k] - aby * g2_y_xxxyz[k];

                        g_yy_xxxzz[k] = g1_y_xxxyzz[k] - aby * g2_y_xxxzz[k];

                        g_yy_xxyyy[k] = g1_y_xxyyyy[k] - aby * g2_y_xxyyy[k];

                        g_yy_xxyyz[k] = g1_y_xxyyyz[k] - aby * g2_y_xxyyz[k];

                        g_yy_xxyzz[k] = g1_y_xxyyzz[k] - aby * g2_y_xxyzz[k];

                        g_yy_xxzzz[k] = g1_y_xxyzzz[k] - aby * g2_y_xxzzz[k];

                        g_yy_xyyyy[k] = g1_y_xyyyyy[k] - aby * g2_y_xyyyy[k];

                        g_yy_xyyyz[k] = g1_y_xyyyyz[k] - aby * g2_y_xyyyz[k];

                        g_yy_xyyzz[k] = g1_y_xyyyzz[k] - aby * g2_y_xyyzz[k];

                        g_yy_xyzzz[k] = g1_y_xyyzzz[k] - aby * g2_y_xyzzz[k];

                        g_yy_xzzzz[k] = g1_y_xyzzzz[k] - aby * g2_y_xzzzz[k];

                        g_yy_yyyyy[k] = g1_y_yyyyyy[k] - aby * g2_y_yyyyy[k];

                        g_yy_yyyyz[k] = g1_y_yyyyyz[k] - aby * g2_y_yyyyz[k];

                        g_yy_yyyzz[k] = g1_y_yyyyzz[k] - aby * g2_y_yyyzz[k];

                        g_yy_yyzzz[k] = g1_y_yyyzzz[k] - aby * g2_y_yyzzz[k];

                        g_yy_yzzzz[k] = g1_y_yyzzzz[k] - aby * g2_y_yzzzz[k];

                        g_yy_zzzzz[k] = g1_y_yzzzzz[k] - aby * g2_y_zzzzz[k];

                        g_yz_xxxxx[k] = g1_z_xxxxxy[k] - aby * g2_z_xxxxx[k];

                        g_yz_xxxxy[k] = g1_z_xxxxyy[k] - aby * g2_z_xxxxy[k];

                        g_yz_xxxxz[k] = g1_z_xxxxyz[k] - aby * g2_z_xxxxz[k];

                        g_yz_xxxyy[k] = g1_z_xxxyyy[k] - aby * g2_z_xxxyy[k];

                        g_yz_xxxyz[k] = g1_z_xxxyyz[k] - aby * g2_z_xxxyz[k];

                        g_yz_xxxzz[k] = g1_z_xxxyzz[k] - aby * g2_z_xxxzz[k];

                        g_yz_xxyyy[k] = g1_z_xxyyyy[k] - aby * g2_z_xxyyy[k];

                        g_yz_xxyyz[k] = g1_z_xxyyyz[k] - aby * g2_z_xxyyz[k];

                        g_yz_xxyzz[k] = g1_z_xxyyzz[k] - aby * g2_z_xxyzz[k];

                        g_yz_xxzzz[k] = g1_z_xxyzzz[k] - aby * g2_z_xxzzz[k];

                        g_yz_xyyyy[k] = g1_z_xyyyyy[k] - aby * g2_z_xyyyy[k];

                        g_yz_xyyyz[k] = g1_z_xyyyyz[k] - aby * g2_z_xyyyz[k];

                        g_yz_xyyzz[k] = g1_z_xyyyzz[k] - aby * g2_z_xyyzz[k];

                        g_yz_xyzzz[k] = g1_z_xyyzzz[k] - aby * g2_z_xyzzz[k];

                        g_yz_xzzzz[k] = g1_z_xyzzzz[k] - aby * g2_z_xzzzz[k];

                        g_yz_yyyyy[k] = g1_z_yyyyyy[k] - aby * g2_z_yyyyy[k];

                        g_yz_yyyyz[k] = g1_z_yyyyyz[k] - aby * g2_z_yyyyz[k];

                        g_yz_yyyzz[k] = g1_z_yyyyzz[k] - aby * g2_z_yyyzz[k];

                        g_yz_yyzzz[k] = g1_z_yyyzzz[k] - aby * g2_z_yyzzz[k];

                        g_yz_yzzzz[k] = g1_z_yyzzzz[k] - aby * g2_z_yzzzz[k];

                        g_yz_zzzzz[k] = g1_z_yzzzzz[k] - aby * g2_z_zzzzz[k];

                        // leading z component

                        g_zz_xxxxx[k] = g1_z_xxxxxz[k] - abz * g2_z_xxxxx[k];

                        g_zz_xxxxy[k] = g1_z_xxxxyz[k] - abz * g2_z_xxxxy[k];

                        g_zz_xxxxz[k] = g1_z_xxxxzz[k] - abz * g2_z_xxxxz[k];

                        g_zz_xxxyy[k] = g1_z_xxxyyz[k] - abz * g2_z_xxxyy[k];

                        g_zz_xxxyz[k] = g1_z_xxxyzz[k] - abz * g2_z_xxxyz[k];

                        g_zz_xxxzz[k] = g1_z_xxxzzz[k] - abz * g2_z_xxxzz[k];

                        g_zz_xxyyy[k] = g1_z_xxyyyz[k] - abz * g2_z_xxyyy[k];

                        g_zz_xxyyz[k] = g1_z_xxyyzz[k] - abz * g2_z_xxyyz[k];

                        g_zz_xxyzz[k] = g1_z_xxyzzz[k] - abz * g2_z_xxyzz[k];

                        g_zz_xxzzz[k] = g1_z_xxzzzz[k] - abz * g2_z_xxzzz[k];

                        g_zz_xyyyy[k] = g1_z_xyyyyz[k] - abz * g2_z_xyyyy[k];

                        g_zz_xyyyz[k] = g1_z_xyyyzz[k] - abz * g2_z_xyyyz[k];

                        g_zz_xyyzz[k] = g1_z_xyyzzz[k] - abz * g2_z_xyyzz[k];

                        g_zz_xyzzz[k] = g1_z_xyzzzz[k] - abz * g2_z_xyzzz[k];

                        g_zz_xzzzz[k] = g1_z_xzzzzz[k] - abz * g2_z_xzzzz[k];

                        g_zz_yyyyy[k] = g1_z_yyyyyz[k] - abz * g2_z_yyyyy[k];

                        g_zz_yyyyz[k] = g1_z_yyyyzz[k] - abz * g2_z_yyyyz[k];

                        g_zz_yyyzz[k] = g1_z_yyyzzz[k] - abz * g2_z_yyyzz[k];

                        g_zz_yyzzz[k] = g1_z_yyzzzz[k] - abz * g2_z_yyzzz[k];

                        g_zz_yzzzz[k] = g1_z_yzzzzz[k] - abz * g2_z_yzzzz[k];

                        g_zz_zzzzz[k] = g1_z_zzzzzz[k] - abz * g2_z_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForDIXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 2) && (recPattern[i].second() == 6))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (26|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {2, 6, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 7, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {1, 6, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (PI|g(r,r')|XX) integrals

                    auto g2_x_xxxxxx = braBuffer.data(g2off + j);

                    auto g2_x_xxxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_x_xxxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_x_xxxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_x_xxxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_x_xxxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_x_xxxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_x_xxxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_x_xxxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_x_xxxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_x_xxyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_x_xxyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_x_xxyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_x_xxyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_x_xxzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_x_xyyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_x_xyyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_x_xyyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_x_xyyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_x_xyzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_x_xzzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_x_yyyyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_x_yyyyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_x_yyyyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_x_yyyzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_x_yyzzzz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_x_yzzzzz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_x_zzzzzz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_y_xxxxxx = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_y_xxxxxy = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_y_xxxxxz = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_y_xxxxyy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_y_xxxxyz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_y_xxxxzz = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_y_xxxyyy = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_y_xxxyyz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_y_xxxyzz = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_y_xxxzzz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_y_xxyyyy = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_y_xxyyyz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_y_xxyyzz = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_y_xxyzzz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_y_xxzzzz = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_y_xyyyyy = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_y_xyyyyz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_y_xyyyzz = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_y_xyyzzz = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_y_xyzzzz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_y_xzzzzz = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_y_yyyyyy = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_y_yyyyyz = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_y_yyyyzz = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_y_yyyzzz = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_y_yyzzzz = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_y_yzzzzz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_y_zzzzzz = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_z_xxxxxx = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_z_xxxxxy = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_z_xxxxxz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_z_xxxxyy = braBuffer.data(g2off + 59 * kcomp + j);

                    auto g2_z_xxxxyz = braBuffer.data(g2off + 60 * kcomp + j);

                    auto g2_z_xxxxzz = braBuffer.data(g2off + 61 * kcomp + j);

                    auto g2_z_xxxyyy = braBuffer.data(g2off + 62 * kcomp + j);

                    auto g2_z_xxxyyz = braBuffer.data(g2off + 63 * kcomp + j);

                    auto g2_z_xxxyzz = braBuffer.data(g2off + 64 * kcomp + j);

                    auto g2_z_xxxzzz = braBuffer.data(g2off + 65 * kcomp + j);

                    auto g2_z_xxyyyy = braBuffer.data(g2off + 66 * kcomp + j);

                    auto g2_z_xxyyyz = braBuffer.data(g2off + 67 * kcomp + j);

                    auto g2_z_xxyyzz = braBuffer.data(g2off + 68 * kcomp + j);

                    auto g2_z_xxyzzz = braBuffer.data(g2off + 69 * kcomp + j);

                    auto g2_z_xxzzzz = braBuffer.data(g2off + 70 * kcomp + j);

                    auto g2_z_xyyyyy = braBuffer.data(g2off + 71 * kcomp + j);

                    auto g2_z_xyyyyz = braBuffer.data(g2off + 72 * kcomp + j);

                    auto g2_z_xyyyzz = braBuffer.data(g2off + 73 * kcomp + j);

                    auto g2_z_xyyzzz = braBuffer.data(g2off + 74 * kcomp + j);

                    auto g2_z_xyzzzz = braBuffer.data(g2off + 75 * kcomp + j);

                    auto g2_z_xzzzzz = braBuffer.data(g2off + 76 * kcomp + j);

                    auto g2_z_yyyyyy = braBuffer.data(g2off + 77 * kcomp + j);

                    auto g2_z_yyyyyz = braBuffer.data(g2off + 78 * kcomp + j);

                    auto g2_z_yyyyzz = braBuffer.data(g2off + 79 * kcomp + j);

                    auto g2_z_yyyzzz = braBuffer.data(g2off + 80 * kcomp + j);

                    auto g2_z_yyzzzz = braBuffer.data(g2off + 81 * kcomp + j);

                    auto g2_z_yzzzzz = braBuffer.data(g2off + 82 * kcomp + j);

                    auto g2_z_zzzzzz = braBuffer.data(g2off + 83 * kcomp + j);

                    // set up pointers to (PK|g(r,r')|XX) integrals

                    auto g1_x_xxxxxxx = braBuffer.data(g1off + j);

                    auto g1_x_xxxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_x_xxxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_x_xxxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_x_xxxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_x_xxxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_x_xxxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_x_xxxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_x_xxxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_x_xxxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_x_xxxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_x_xxxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_x_xxxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_x_xxxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_x_xxxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_x_xxyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_x_xxyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_x_xxyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_x_xxyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_x_xxyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_x_xxzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_x_xyyyyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_x_xyyyyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_x_xyyyyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_x_xyyyzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_x_xyyzzzz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_x_xyzzzzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_x_xzzzzzz = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_y_xxxxxxx = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_y_xxxxxxy = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_y_xxxxxxz = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_y_xxxxxyy = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_y_xxxxxyz = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_y_xxxxxzz = braBuffer.data(g1off + 41 * kcomp + j);

                    auto g1_y_xxxxyyy = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_y_xxxxyyz = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_y_xxxxyzz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_y_xxxxzzz = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_y_xxxyyyy = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_y_xxxyyyz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_y_xxxyyzz = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_y_xxxyzzz = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_y_xxxzzzz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_y_xxyyyyy = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_y_xxyyyyz = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_y_xxyyyzz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_y_xxyyzzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_y_xxyzzzz = braBuffer.data(g1off + 55 * kcomp + j);

                    auto g1_y_xxzzzzz = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_y_xyyyyyy = braBuffer.data(g1off + 57 * kcomp + j);

                    auto g1_y_xyyyyyz = braBuffer.data(g1off + 58 * kcomp + j);

                    auto g1_y_xyyyyzz = braBuffer.data(g1off + 59 * kcomp + j);

                    auto g1_y_xyyyzzz = braBuffer.data(g1off + 60 * kcomp + j);

                    auto g1_y_xyyzzzz = braBuffer.data(g1off + 61 * kcomp + j);

                    auto g1_y_xyzzzzz = braBuffer.data(g1off + 62 * kcomp + j);

                    auto g1_y_xzzzzzz = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_y_yyyyyyy = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_y_yyyyyyz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_y_yyyyyzz = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_y_yyyyzzz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_y_yyyzzzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_y_yyzzzzz = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_y_yzzzzzz = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_z_xxxxxxx = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_z_xxxxxxy = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_z_xxxxxxz = braBuffer.data(g1off + 74 * kcomp + j);

                    auto g1_z_xxxxxyy = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_z_xxxxxyz = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_z_xxxxxzz = braBuffer.data(g1off + 77 * kcomp + j);

                    auto g1_z_xxxxyyy = braBuffer.data(g1off + 78 * kcomp + j);

                    auto g1_z_xxxxyyz = braBuffer.data(g1off + 79 * kcomp + j);

                    auto g1_z_xxxxyzz = braBuffer.data(g1off + 80 * kcomp + j);

                    auto g1_z_xxxxzzz = braBuffer.data(g1off + 81 * kcomp + j);

                    auto g1_z_xxxyyyy = braBuffer.data(g1off + 82 * kcomp + j);

                    auto g1_z_xxxyyyz = braBuffer.data(g1off + 83 * kcomp + j);

                    auto g1_z_xxxyyzz = braBuffer.data(g1off + 84 * kcomp + j);

                    auto g1_z_xxxyzzz = braBuffer.data(g1off + 85 * kcomp + j);

                    auto g1_z_xxxzzzz = braBuffer.data(g1off + 86 * kcomp + j);

                    auto g1_z_xxyyyyy = braBuffer.data(g1off + 87 * kcomp + j);

                    auto g1_z_xxyyyyz = braBuffer.data(g1off + 88 * kcomp + j);

                    auto g1_z_xxyyyzz = braBuffer.data(g1off + 89 * kcomp + j);

                    auto g1_z_xxyyzzz = braBuffer.data(g1off + 90 * kcomp + j);

                    auto g1_z_xxyzzzz = braBuffer.data(g1off + 91 * kcomp + j);

                    auto g1_z_xxzzzzz = braBuffer.data(g1off + 92 * kcomp + j);

                    auto g1_z_xyyyyyy = braBuffer.data(g1off + 93 * kcomp + j);

                    auto g1_z_xyyyyyz = braBuffer.data(g1off + 94 * kcomp + j);

                    auto g1_z_xyyyyzz = braBuffer.data(g1off + 95 * kcomp + j);

                    auto g1_z_xyyyzzz = braBuffer.data(g1off + 96 * kcomp + j);

                    auto g1_z_xyyzzzz = braBuffer.data(g1off + 97 * kcomp + j);

                    auto g1_z_xyzzzzz = braBuffer.data(g1off + 98 * kcomp + j);

                    auto g1_z_xzzzzzz = braBuffer.data(g1off + 99 * kcomp + j);

                    auto g1_z_yyyyyyy = braBuffer.data(g1off + 100 * kcomp + j);

                    auto g1_z_yyyyyyz = braBuffer.data(g1off + 101 * kcomp + j);

                    auto g1_z_yyyyyzz = braBuffer.data(g1off + 102 * kcomp + j);

                    auto g1_z_yyyyzzz = braBuffer.data(g1off + 103 * kcomp + j);

                    auto g1_z_yyyzzzz = braBuffer.data(g1off + 104 * kcomp + j);

                    auto g1_z_yyzzzzz = braBuffer.data(g1off + 105 * kcomp + j);

                    auto g1_z_yzzzzzz = braBuffer.data(g1off + 106 * kcomp + j);

                    auto g1_z_zzzzzzz = braBuffer.data(g1off + 107 * kcomp + j);

                    // set up pointers to (DI|g(r,r')|XX) integrals

                    auto g_xx_xxxxxx = braBuffer.data(goff + j);

                    auto g_xx_xxxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xx_xxxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xx_xxxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xx_xxxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xx_xxxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xx_xxxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xx_xxxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xx_xxxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xx_xxxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xx_xxyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xx_xxyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xx_xxyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xx_xxyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xx_xxzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xx_xyyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xx_xyyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xx_xyyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xx_xyyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xx_xyzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xx_xzzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xx_yyyyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xx_yyyyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xx_yyyyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xx_yyyzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xx_yyzzzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xx_yzzzzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xx_zzzzzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xy_xxxxxx = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xy_xxxxxy = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xy_xxxxxz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xy_xxxxyy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xy_xxxxyz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xy_xxxxzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xy_xxxyyy = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xy_xxxyyz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xy_xxxyzz = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xy_xxxzzz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xy_xxyyyy = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xy_xxyyyz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xy_xxyyzz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xy_xxyzzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xy_xxzzzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xy_xyyyyy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xy_xyyyyz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xy_xyyyzz = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xy_xyyzzz = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xy_xyzzzz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xy_xzzzzz = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xy_yyyyyy = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xy_yyyyyz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xy_yyyyzz = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xy_yyyzzz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xy_yyzzzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xy_yzzzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xy_zzzzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xz_xxxxxx = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xz_xxxxxy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xz_xxxxxz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xz_xxxxyy = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_xz_xxxxyz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_xz_xxxxzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_xz_xxxyyy = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_xz_xxxyyz = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_xz_xxxyzz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_xz_xxxzzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_xz_xxyyyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_xz_xxyyyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_xz_xxyyzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_xz_xxyzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_xz_xxzzzz = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_xz_xyyyyy = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_xz_xyyyyz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_xz_xyyyzz = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_xz_xyyzzz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_xz_xyzzzz = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_xz_xzzzzz = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_xz_yyyyyy = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_xz_yyyyyz = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_xz_yyyyzz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_xz_yyyzzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_xz_yyzzzz = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_xz_yzzzzz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_xz_zzzzzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_yy_xxxxxx = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_yy_xxxxxy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_yy_xxxxxz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_yy_xxxxyy = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_yy_xxxxyz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_yy_xxxxzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_yy_xxxyyy = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_yy_xxxyyz = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_yy_xxxyzz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_yy_xxxzzz = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_yy_xxyyyy = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_yy_xxyyyz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_yy_xxyyzz = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_yy_xxyzzz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_yy_xxzzzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_yy_xyyyyy = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_yy_xyyyyz = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_yy_xyyyzz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_yy_xyyzzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_yy_xyzzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_yy_xzzzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_yy_yyyyyy = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_yy_yyyyyz = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_yy_yyyyzz = braBuffer.data(goff + 107 * kcomp + j);

                    auto g_yy_yyyzzz = braBuffer.data(goff + 108 * kcomp + j);

                    auto g_yy_yyzzzz = braBuffer.data(goff + 109 * kcomp + j);

                    auto g_yy_yzzzzz = braBuffer.data(goff + 110 * kcomp + j);

                    auto g_yy_zzzzzz = braBuffer.data(goff + 111 * kcomp + j);

                    auto g_yz_xxxxxx = braBuffer.data(goff + 112 * kcomp + j);

                    auto g_yz_xxxxxy = braBuffer.data(goff + 113 * kcomp + j);

                    auto g_yz_xxxxxz = braBuffer.data(goff + 114 * kcomp + j);

                    auto g_yz_xxxxyy = braBuffer.data(goff + 115 * kcomp + j);

                    auto g_yz_xxxxyz = braBuffer.data(goff + 116 * kcomp + j);

                    auto g_yz_xxxxzz = braBuffer.data(goff + 117 * kcomp + j);

                    auto g_yz_xxxyyy = braBuffer.data(goff + 118 * kcomp + j);

                    auto g_yz_xxxyyz = braBuffer.data(goff + 119 * kcomp + j);

                    auto g_yz_xxxyzz = braBuffer.data(goff + 120 * kcomp + j);

                    auto g_yz_xxxzzz = braBuffer.data(goff + 121 * kcomp + j);

                    auto g_yz_xxyyyy = braBuffer.data(goff + 122 * kcomp + j);

                    auto g_yz_xxyyyz = braBuffer.data(goff + 123 * kcomp + j);

                    auto g_yz_xxyyzz = braBuffer.data(goff + 124 * kcomp + j);

                    auto g_yz_xxyzzz = braBuffer.data(goff + 125 * kcomp + j);

                    auto g_yz_xxzzzz = braBuffer.data(goff + 126 * kcomp + j);

                    auto g_yz_xyyyyy = braBuffer.data(goff + 127 * kcomp + j);

                    auto g_yz_xyyyyz = braBuffer.data(goff + 128 * kcomp + j);

                    auto g_yz_xyyyzz = braBuffer.data(goff + 129 * kcomp + j);

                    auto g_yz_xyyzzz = braBuffer.data(goff + 130 * kcomp + j);

                    auto g_yz_xyzzzz = braBuffer.data(goff + 131 * kcomp + j);

                    auto g_yz_xzzzzz = braBuffer.data(goff + 132 * kcomp + j);

                    auto g_yz_yyyyyy = braBuffer.data(goff + 133 * kcomp + j);

                    auto g_yz_yyyyyz = braBuffer.data(goff + 134 * kcomp + j);

                    auto g_yz_yyyyzz = braBuffer.data(goff + 135 * kcomp + j);

                    auto g_yz_yyyzzz = braBuffer.data(goff + 136 * kcomp + j);

                    auto g_yz_yyzzzz = braBuffer.data(goff + 137 * kcomp + j);

                    auto g_yz_yzzzzz = braBuffer.data(goff + 138 * kcomp + j);

                    auto g_yz_zzzzzz = braBuffer.data(goff + 139 * kcomp + j);

                    auto g_zz_xxxxxx = braBuffer.data(goff + 140 * kcomp + j);

                    auto g_zz_xxxxxy = braBuffer.data(goff + 141 * kcomp + j);

                    auto g_zz_xxxxxz = braBuffer.data(goff + 142 * kcomp + j);

                    auto g_zz_xxxxyy = braBuffer.data(goff + 143 * kcomp + j);

                    auto g_zz_xxxxyz = braBuffer.data(goff + 144 * kcomp + j);

                    auto g_zz_xxxxzz = braBuffer.data(goff + 145 * kcomp + j);

                    auto g_zz_xxxyyy = braBuffer.data(goff + 146 * kcomp + j);

                    auto g_zz_xxxyyz = braBuffer.data(goff + 147 * kcomp + j);

                    auto g_zz_xxxyzz = braBuffer.data(goff + 148 * kcomp + j);

                    auto g_zz_xxxzzz = braBuffer.data(goff + 149 * kcomp + j);

                    auto g_zz_xxyyyy = braBuffer.data(goff + 150 * kcomp + j);

                    auto g_zz_xxyyyz = braBuffer.data(goff + 151 * kcomp + j);

                    auto g_zz_xxyyzz = braBuffer.data(goff + 152 * kcomp + j);

                    auto g_zz_xxyzzz = braBuffer.data(goff + 153 * kcomp + j);

                    auto g_zz_xxzzzz = braBuffer.data(goff + 154 * kcomp + j);

                    auto g_zz_xyyyyy = braBuffer.data(goff + 155 * kcomp + j);

                    auto g_zz_xyyyyz = braBuffer.data(goff + 156 * kcomp + j);

                    auto g_zz_xyyyzz = braBuffer.data(goff + 157 * kcomp + j);

                    auto g_zz_xyyzzz = braBuffer.data(goff + 158 * kcomp + j);

                    auto g_zz_xyzzzz = braBuffer.data(goff + 159 * kcomp + j);

                    auto g_zz_xzzzzz = braBuffer.data(goff + 160 * kcomp + j);

                    auto g_zz_yyyyyy = braBuffer.data(goff + 161 * kcomp + j);

                    auto g_zz_yyyyyz = braBuffer.data(goff + 162 * kcomp + j);

                    auto g_zz_yyyyzz = braBuffer.data(goff + 163 * kcomp + j);

                    auto g_zz_yyyzzz = braBuffer.data(goff + 164 * kcomp + j);

                    auto g_zz_yyzzzz = braBuffer.data(goff + 165 * kcomp + j);

                    auto g_zz_yzzzzz = braBuffer.data(goff + 166 * kcomp + j);

                    auto g_zz_zzzzzz = braBuffer.data(goff + 167 * kcomp + j);

                    #pragma omp simd aligned(g2_x_xxxxxx, g2_x_xxxxxy, g2_x_xxxxxz,\
                                             g2_x_xxxxyy, g2_x_xxxxyz, g2_x_xxxxzz,\
                                             g2_x_xxxyyy, g2_x_xxxyyz, g2_x_xxxyzz,\
                                             g2_x_xxxzzz, g2_x_xxyyyy, g2_x_xxyyyz,\
                                             g2_x_xxyyzz, g2_x_xxyzzz, g2_x_xxzzzz,\
                                             g2_x_xyyyyy, g2_x_xyyyyz, g2_x_xyyyzz,\
                                             g2_x_xyyzzz, g2_x_xyzzzz, g2_x_xzzzzz,\
                                             g2_x_yyyyyy, g2_x_yyyyyz, g2_x_yyyyzz,\
                                             g2_x_yyyzzz, g2_x_yyzzzz, g2_x_yzzzzz,\
                                             g2_x_zzzzzz, g2_y_xxxxxx, g2_y_xxxxxy,\
                                             g2_y_xxxxxz, g2_y_xxxxyy, g2_y_xxxxyz,\
                                             g2_y_xxxxzz, g2_y_xxxyyy, g2_y_xxxyyz,\
                                             g2_y_xxxyzz, g2_y_xxxzzz, g2_y_xxyyyy,\
                                             g2_y_xxyyyz, g2_y_xxyyzz, g2_y_xxyzzz,\
                                             g2_y_xxzzzz, g2_y_xyyyyy, g2_y_xyyyyz,\
                                             g2_y_xyyyzz, g2_y_xyyzzz, g2_y_xyzzzz,\
                                             g2_y_xzzzzz, g2_y_yyyyyy, g2_y_yyyyyz,\
                                             g2_y_yyyyzz, g2_y_yyyzzz, g2_y_yyzzzz,\
                                             g2_y_yzzzzz, g2_y_zzzzzz, g2_z_xxxxxx,\
                                             g2_z_xxxxxy, g2_z_xxxxxz, g2_z_xxxxyy,\
                                             g2_z_xxxxyz, g2_z_xxxxzz, g2_z_xxxyyy,\
                                             g2_z_xxxyyz, g2_z_xxxyzz, g2_z_xxxzzz,\
                                             g2_z_xxyyyy, g2_z_xxyyyz, g2_z_xxyyzz,\
                                             g2_z_xxyzzz, g2_z_xxzzzz, g2_z_xyyyyy,\
                                             g2_z_xyyyyz, g2_z_xyyyzz, g2_z_xyyzzz,\
                                             g2_z_xyzzzz, g2_z_xzzzzz, g2_z_yyyyyy,\
                                             g2_z_yyyyyz, g2_z_yyyyzz, g2_z_yyyzzz,\
                                             g2_z_yyzzzz, g2_z_yzzzzz, g2_z_zzzzzz,\
                                             g1_x_xxxxxxx, g1_x_xxxxxxy, g1_x_xxxxxxz,\
                                             g1_x_xxxxxyy, g1_x_xxxxxyz, g1_x_xxxxxzz,\
                                             g1_x_xxxxyyy, g1_x_xxxxyyz, g1_x_xxxxyzz,\
                                             g1_x_xxxxzzz, g1_x_xxxyyyy, g1_x_xxxyyyz,\
                                             g1_x_xxxyyzz, g1_x_xxxyzzz, g1_x_xxxzzzz,\
                                             g1_x_xxyyyyy, g1_x_xxyyyyz, g1_x_xxyyyzz,\
                                             g1_x_xxyyzzz, g1_x_xxyzzzz, g1_x_xxzzzzz,\
                                             g1_x_xyyyyyy, g1_x_xyyyyyz, g1_x_xyyyyzz,\
                                             g1_x_xyyyzzz, g1_x_xyyzzzz, g1_x_xyzzzzz,\
                                             g1_x_xzzzzzz, \
                                             g1_y_xxxxxxx, g1_y_xxxxxxy, g1_y_xxxxxxz,\
                                             g1_y_xxxxxyy, g1_y_xxxxxyz, g1_y_xxxxxzz,\
                                             g1_y_xxxxyyy, g1_y_xxxxyyz, g1_y_xxxxyzz,\
                                             g1_y_xxxxzzz, g1_y_xxxyyyy, g1_y_xxxyyyz,\
                                             g1_y_xxxyyzz, g1_y_xxxyzzz, g1_y_xxxzzzz,\
                                             g1_y_xxyyyyy, g1_y_xxyyyyz, g1_y_xxyyyzz,\
                                             g1_y_xxyyzzz, g1_y_xxyzzzz, g1_y_xxzzzzz,\
                                             g1_y_xyyyyyy, g1_y_xyyyyyz, g1_y_xyyyyzz,\
                                             g1_y_xyyyzzz, g1_y_xyyzzzz, g1_y_xyzzzzz,\
                                             g1_y_xzzzzzz, g1_y_yyyyyyy, g1_y_yyyyyyz,\
                                             g1_y_yyyyyzz, g1_y_yyyyzzz, g1_y_yyyzzzz,\
                                             g1_y_yyzzzzz, g1_y_yzzzzzz, \
                                             g1_z_xxxxxxx, g1_z_xxxxxxy, g1_z_xxxxxxz,\
                                             g1_z_xxxxxyy, g1_z_xxxxxyz, g1_z_xxxxxzz,\
                                             g1_z_xxxxyyy, g1_z_xxxxyyz, g1_z_xxxxyzz,\
                                             g1_z_xxxxzzz, g1_z_xxxyyyy, g1_z_xxxyyyz,\
                                             g1_z_xxxyyzz, g1_z_xxxyzzz, g1_z_xxxzzzz,\
                                             g1_z_xxyyyyy, g1_z_xxyyyyz, g1_z_xxyyyzz,\
                                             g1_z_xxyyzzz, g1_z_xxyzzzz, g1_z_xxzzzzz,\
                                             g1_z_xyyyyyy, g1_z_xyyyyyz, g1_z_xyyyyzz,\
                                             g1_z_xyyyzzz, g1_z_xyyzzzz, g1_z_xyzzzzz,\
                                             g1_z_xzzzzzz, g1_z_yyyyyyy, g1_z_yyyyyyz,\
                                             g1_z_yyyyyzz, g1_z_yyyyzzz, g1_z_yyyzzzz,\
                                             g1_z_yyzzzzz, g1_z_yzzzzzz, g1_z_zzzzzzz,\
                                             g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz,\
                                             g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz,\
                                             g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyzz,\
                                             g_xx_xxxzzz, g_xx_xxyyyy, g_xx_xxyyyz,\
                                             g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz,\
                                             g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz,\
                                             g_xx_xyyzzz, g_xx_xyzzzz, g_xx_xzzzzz,\
                                             g_xx_yyyyyy, g_xx_yyyyyz, g_xx_yyyyzz,\
                                             g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz,\
                                             g_xx_zzzzzz, g_xy_xxxxxx, g_xy_xxxxxy,\
                                             g_xy_xxxxxz, g_xy_xxxxyy, g_xy_xxxxyz,\
                                             g_xy_xxxxzz, g_xy_xxxyyy, g_xy_xxxyyz,\
                                             g_xy_xxxyzz, g_xy_xxxzzz, g_xy_xxyyyy,\
                                             g_xy_xxyyyz, g_xy_xxyyzz, g_xy_xxyzzz,\
                                             g_xy_xxzzzz, g_xy_xyyyyy, g_xy_xyyyyz,\
                                             g_xy_xyyyzz, g_xy_xyyzzz, g_xy_xyzzzz,\
                                             g_xy_xzzzzz, g_xy_yyyyyy, g_xy_yyyyyz,\
                                             g_xy_yyyyzz, g_xy_yyyzzz, g_xy_yyzzzz,\
                                             g_xy_yzzzzz, g_xy_zzzzzz, g_xz_xxxxxx,\
                                             g_xz_xxxxxy, g_xz_xxxxxz, g_xz_xxxxyy,\
                                             g_xz_xxxxyz, g_xz_xxxxzz, g_xz_xxxyyy,\
                                             g_xz_xxxyyz, g_xz_xxxyzz, g_xz_xxxzzz,\
                                             g_xz_xxyyyy, g_xz_xxyyyz, g_xz_xxyyzz,\
                                             g_xz_xxyzzz, g_xz_xxzzzz, g_xz_xyyyyy,\
                                             g_xz_xyyyyz, g_xz_xyyyzz, g_xz_xyyzzz,\
                                             g_xz_xyzzzz, g_xz_xzzzzz, g_xz_yyyyyy,\
                                             g_xz_yyyyyz, g_xz_yyyyzz, g_xz_yyyzzz,\
                                             g_xz_yyzzzz, g_xz_yzzzzz, g_xz_zzzzzz,\
                                             g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz,\
                                             g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz,\
                                             g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyzz,\
                                             g_yy_xxxzzz, g_yy_xxyyyy, g_yy_xxyyyz,\
                                             g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz,\
                                             g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz,\
                                             g_yy_xyyzzz, g_yy_xyzzzz, g_yy_xzzzzz,\
                                             g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyzz,\
                                             g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz,\
                                             g_yy_zzzzzz, g_yz_xxxxxx, g_yz_xxxxxy,\
                                             g_yz_xxxxxz, g_yz_xxxxyy, g_yz_xxxxyz,\
                                             g_yz_xxxxzz, g_yz_xxxyyy, g_yz_xxxyyz,\
                                             g_yz_xxxyzz, g_yz_xxxzzz, g_yz_xxyyyy,\
                                             g_yz_xxyyyz, g_yz_xxyyzz, g_yz_xxyzzz,\
                                             g_yz_xxzzzz, g_yz_xyyyyy, g_yz_xyyyyz,\
                                             g_yz_xyyyzz, g_yz_xyyzzz, g_yz_xyzzzz,\
                                             g_yz_xzzzzz, g_yz_yyyyyy, g_yz_yyyyyz,\
                                             g_yz_yyyyzz, g_yz_yyyzzz, g_yz_yyzzzz,\
                                             g_yz_yzzzzz, g_yz_zzzzzz, g_zz_xxxxxx,\
                                             g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxyy,\
                                             g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy,\
                                             g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz,\
                                             g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyzz,\
                                             g_zz_xxyzzz, g_zz_xxzzzz, g_zz_xyyyyy,\
                                             g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz,\
                                             g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy,\
                                             g_zz_yyyyyz, g_zz_yyyyzz, g_zz_yyyzzz,\
                                             g_zz_yyzzzz, g_zz_yzzzzz, g_zz_zzzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xx_xxxxxx[k] = g1_x_xxxxxxx[k] - abx * g2_x_xxxxxx[k];

                        g_xx_xxxxxy[k] = g1_x_xxxxxxy[k] - abx * g2_x_xxxxxy[k];

                        g_xx_xxxxxz[k] = g1_x_xxxxxxz[k] - abx * g2_x_xxxxxz[k];

                        g_xx_xxxxyy[k] = g1_x_xxxxxyy[k] - abx * g2_x_xxxxyy[k];

                        g_xx_xxxxyz[k] = g1_x_xxxxxyz[k] - abx * g2_x_xxxxyz[k];

                        g_xx_xxxxzz[k] = g1_x_xxxxxzz[k] - abx * g2_x_xxxxzz[k];

                        g_xx_xxxyyy[k] = g1_x_xxxxyyy[k] - abx * g2_x_xxxyyy[k];

                        g_xx_xxxyyz[k] = g1_x_xxxxyyz[k] - abx * g2_x_xxxyyz[k];

                        g_xx_xxxyzz[k] = g1_x_xxxxyzz[k] - abx * g2_x_xxxyzz[k];

                        g_xx_xxxzzz[k] = g1_x_xxxxzzz[k] - abx * g2_x_xxxzzz[k];

                        g_xx_xxyyyy[k] = g1_x_xxxyyyy[k] - abx * g2_x_xxyyyy[k];

                        g_xx_xxyyyz[k] = g1_x_xxxyyyz[k] - abx * g2_x_xxyyyz[k];

                        g_xx_xxyyzz[k] = g1_x_xxxyyzz[k] - abx * g2_x_xxyyzz[k];

                        g_xx_xxyzzz[k] = g1_x_xxxyzzz[k] - abx * g2_x_xxyzzz[k];

                        g_xx_xxzzzz[k] = g1_x_xxxzzzz[k] - abx * g2_x_xxzzzz[k];

                        g_xx_xyyyyy[k] = g1_x_xxyyyyy[k] - abx * g2_x_xyyyyy[k];

                        g_xx_xyyyyz[k] = g1_x_xxyyyyz[k] - abx * g2_x_xyyyyz[k];

                        g_xx_xyyyzz[k] = g1_x_xxyyyzz[k] - abx * g2_x_xyyyzz[k];

                        g_xx_xyyzzz[k] = g1_x_xxyyzzz[k] - abx * g2_x_xyyzzz[k];

                        g_xx_xyzzzz[k] = g1_x_xxyzzzz[k] - abx * g2_x_xyzzzz[k];

                        g_xx_xzzzzz[k] = g1_x_xxzzzzz[k] - abx * g2_x_xzzzzz[k];

                        g_xx_yyyyyy[k] = g1_x_xyyyyyy[k] - abx * g2_x_yyyyyy[k];

                        g_xx_yyyyyz[k] = g1_x_xyyyyyz[k] - abx * g2_x_yyyyyz[k];

                        g_xx_yyyyzz[k] = g1_x_xyyyyzz[k] - abx * g2_x_yyyyzz[k];

                        g_xx_yyyzzz[k] = g1_x_xyyyzzz[k] - abx * g2_x_yyyzzz[k];

                        g_xx_yyzzzz[k] = g1_x_xyyzzzz[k] - abx * g2_x_yyzzzz[k];

                        g_xx_yzzzzz[k] = g1_x_xyzzzzz[k] - abx * g2_x_yzzzzz[k];

                        g_xx_zzzzzz[k] = g1_x_xzzzzzz[k] - abx * g2_x_zzzzzz[k];

                        g_xy_xxxxxx[k] = g1_y_xxxxxxx[k] - abx * g2_y_xxxxxx[k];

                        g_xy_xxxxxy[k] = g1_y_xxxxxxy[k] - abx * g2_y_xxxxxy[k];

                        g_xy_xxxxxz[k] = g1_y_xxxxxxz[k] - abx * g2_y_xxxxxz[k];

                        g_xy_xxxxyy[k] = g1_y_xxxxxyy[k] - abx * g2_y_xxxxyy[k];

                        g_xy_xxxxyz[k] = g1_y_xxxxxyz[k] - abx * g2_y_xxxxyz[k];

                        g_xy_xxxxzz[k] = g1_y_xxxxxzz[k] - abx * g2_y_xxxxzz[k];

                        g_xy_xxxyyy[k] = g1_y_xxxxyyy[k] - abx * g2_y_xxxyyy[k];

                        g_xy_xxxyyz[k] = g1_y_xxxxyyz[k] - abx * g2_y_xxxyyz[k];

                        g_xy_xxxyzz[k] = g1_y_xxxxyzz[k] - abx * g2_y_xxxyzz[k];

                        g_xy_xxxzzz[k] = g1_y_xxxxzzz[k] - abx * g2_y_xxxzzz[k];

                        g_xy_xxyyyy[k] = g1_y_xxxyyyy[k] - abx * g2_y_xxyyyy[k];

                        g_xy_xxyyyz[k] = g1_y_xxxyyyz[k] - abx * g2_y_xxyyyz[k];

                        g_xy_xxyyzz[k] = g1_y_xxxyyzz[k] - abx * g2_y_xxyyzz[k];

                        g_xy_xxyzzz[k] = g1_y_xxxyzzz[k] - abx * g2_y_xxyzzz[k];

                        g_xy_xxzzzz[k] = g1_y_xxxzzzz[k] - abx * g2_y_xxzzzz[k];

                        g_xy_xyyyyy[k] = g1_y_xxyyyyy[k] - abx * g2_y_xyyyyy[k];

                        g_xy_xyyyyz[k] = g1_y_xxyyyyz[k] - abx * g2_y_xyyyyz[k];

                        g_xy_xyyyzz[k] = g1_y_xxyyyzz[k] - abx * g2_y_xyyyzz[k];

                        g_xy_xyyzzz[k] = g1_y_xxyyzzz[k] - abx * g2_y_xyyzzz[k];

                        g_xy_xyzzzz[k] = g1_y_xxyzzzz[k] - abx * g2_y_xyzzzz[k];

                        g_xy_xzzzzz[k] = g1_y_xxzzzzz[k] - abx * g2_y_xzzzzz[k];

                        g_xy_yyyyyy[k] = g1_y_xyyyyyy[k] - abx * g2_y_yyyyyy[k];

                        g_xy_yyyyyz[k] = g1_y_xyyyyyz[k] - abx * g2_y_yyyyyz[k];

                        g_xy_yyyyzz[k] = g1_y_xyyyyzz[k] - abx * g2_y_yyyyzz[k];

                        g_xy_yyyzzz[k] = g1_y_xyyyzzz[k] - abx * g2_y_yyyzzz[k];

                        g_xy_yyzzzz[k] = g1_y_xyyzzzz[k] - abx * g2_y_yyzzzz[k];

                        g_xy_yzzzzz[k] = g1_y_xyzzzzz[k] - abx * g2_y_yzzzzz[k];

                        g_xy_zzzzzz[k] = g1_y_xzzzzzz[k] - abx * g2_y_zzzzzz[k];

                        g_xz_xxxxxx[k] = g1_z_xxxxxxx[k] - abx * g2_z_xxxxxx[k];

                        g_xz_xxxxxy[k] = g1_z_xxxxxxy[k] - abx * g2_z_xxxxxy[k];

                        g_xz_xxxxxz[k] = g1_z_xxxxxxz[k] - abx * g2_z_xxxxxz[k];

                        g_xz_xxxxyy[k] = g1_z_xxxxxyy[k] - abx * g2_z_xxxxyy[k];

                        g_xz_xxxxyz[k] = g1_z_xxxxxyz[k] - abx * g2_z_xxxxyz[k];

                        g_xz_xxxxzz[k] = g1_z_xxxxxzz[k] - abx * g2_z_xxxxzz[k];

                        g_xz_xxxyyy[k] = g1_z_xxxxyyy[k] - abx * g2_z_xxxyyy[k];

                        g_xz_xxxyyz[k] = g1_z_xxxxyyz[k] - abx * g2_z_xxxyyz[k];

                        g_xz_xxxyzz[k] = g1_z_xxxxyzz[k] - abx * g2_z_xxxyzz[k];

                        g_xz_xxxzzz[k] = g1_z_xxxxzzz[k] - abx * g2_z_xxxzzz[k];

                        g_xz_xxyyyy[k] = g1_z_xxxyyyy[k] - abx * g2_z_xxyyyy[k];

                        g_xz_xxyyyz[k] = g1_z_xxxyyyz[k] - abx * g2_z_xxyyyz[k];

                        g_xz_xxyyzz[k] = g1_z_xxxyyzz[k] - abx * g2_z_xxyyzz[k];

                        g_xz_xxyzzz[k] = g1_z_xxxyzzz[k] - abx * g2_z_xxyzzz[k];

                        g_xz_xxzzzz[k] = g1_z_xxxzzzz[k] - abx * g2_z_xxzzzz[k];

                        g_xz_xyyyyy[k] = g1_z_xxyyyyy[k] - abx * g2_z_xyyyyy[k];

                        g_xz_xyyyyz[k] = g1_z_xxyyyyz[k] - abx * g2_z_xyyyyz[k];

                        g_xz_xyyyzz[k] = g1_z_xxyyyzz[k] - abx * g2_z_xyyyzz[k];

                        g_xz_xyyzzz[k] = g1_z_xxyyzzz[k] - abx * g2_z_xyyzzz[k];

                        g_xz_xyzzzz[k] = g1_z_xxyzzzz[k] - abx * g2_z_xyzzzz[k];

                        g_xz_xzzzzz[k] = g1_z_xxzzzzz[k] - abx * g2_z_xzzzzz[k];

                        g_xz_yyyyyy[k] = g1_z_xyyyyyy[k] - abx * g2_z_yyyyyy[k];

                        g_xz_yyyyyz[k] = g1_z_xyyyyyz[k] - abx * g2_z_yyyyyz[k];

                        g_xz_yyyyzz[k] = g1_z_xyyyyzz[k] - abx * g2_z_yyyyzz[k];

                        g_xz_yyyzzz[k] = g1_z_xyyyzzz[k] - abx * g2_z_yyyzzz[k];

                        g_xz_yyzzzz[k] = g1_z_xyyzzzz[k] - abx * g2_z_yyzzzz[k];

                        g_xz_yzzzzz[k] = g1_z_xyzzzzz[k] - abx * g2_z_yzzzzz[k];

                        g_xz_zzzzzz[k] = g1_z_xzzzzzz[k] - abx * g2_z_zzzzzz[k];

                        // leading y component

                        g_yy_xxxxxx[k] = g1_y_xxxxxxy[k] - aby * g2_y_xxxxxx[k];

                        g_yy_xxxxxy[k] = g1_y_xxxxxyy[k] - aby * g2_y_xxxxxy[k];

                        g_yy_xxxxxz[k] = g1_y_xxxxxyz[k] - aby * g2_y_xxxxxz[k];

                        g_yy_xxxxyy[k] = g1_y_xxxxyyy[k] - aby * g2_y_xxxxyy[k];

                        g_yy_xxxxyz[k] = g1_y_xxxxyyz[k] - aby * g2_y_xxxxyz[k];

                        g_yy_xxxxzz[k] = g1_y_xxxxyzz[k] - aby * g2_y_xxxxzz[k];

                        g_yy_xxxyyy[k] = g1_y_xxxyyyy[k] - aby * g2_y_xxxyyy[k];

                        g_yy_xxxyyz[k] = g1_y_xxxyyyz[k] - aby * g2_y_xxxyyz[k];

                        g_yy_xxxyzz[k] = g1_y_xxxyyzz[k] - aby * g2_y_xxxyzz[k];

                        g_yy_xxxzzz[k] = g1_y_xxxyzzz[k] - aby * g2_y_xxxzzz[k];

                        g_yy_xxyyyy[k] = g1_y_xxyyyyy[k] - aby * g2_y_xxyyyy[k];

                        g_yy_xxyyyz[k] = g1_y_xxyyyyz[k] - aby * g2_y_xxyyyz[k];

                        g_yy_xxyyzz[k] = g1_y_xxyyyzz[k] - aby * g2_y_xxyyzz[k];

                        g_yy_xxyzzz[k] = g1_y_xxyyzzz[k] - aby * g2_y_xxyzzz[k];

                        g_yy_xxzzzz[k] = g1_y_xxyzzzz[k] - aby * g2_y_xxzzzz[k];

                        g_yy_xyyyyy[k] = g1_y_xyyyyyy[k] - aby * g2_y_xyyyyy[k];

                        g_yy_xyyyyz[k] = g1_y_xyyyyyz[k] - aby * g2_y_xyyyyz[k];

                        g_yy_xyyyzz[k] = g1_y_xyyyyzz[k] - aby * g2_y_xyyyzz[k];

                        g_yy_xyyzzz[k] = g1_y_xyyyzzz[k] - aby * g2_y_xyyzzz[k];

                        g_yy_xyzzzz[k] = g1_y_xyyzzzz[k] - aby * g2_y_xyzzzz[k];

                        g_yy_xzzzzz[k] = g1_y_xyzzzzz[k] - aby * g2_y_xzzzzz[k];

                        g_yy_yyyyyy[k] = g1_y_yyyyyyy[k] - aby * g2_y_yyyyyy[k];

                        g_yy_yyyyyz[k] = g1_y_yyyyyyz[k] - aby * g2_y_yyyyyz[k];

                        g_yy_yyyyzz[k] = g1_y_yyyyyzz[k] - aby * g2_y_yyyyzz[k];

                        g_yy_yyyzzz[k] = g1_y_yyyyzzz[k] - aby * g2_y_yyyzzz[k];

                        g_yy_yyzzzz[k] = g1_y_yyyzzzz[k] - aby * g2_y_yyzzzz[k];

                        g_yy_yzzzzz[k] = g1_y_yyzzzzz[k] - aby * g2_y_yzzzzz[k];

                        g_yy_zzzzzz[k] = g1_y_yzzzzzz[k] - aby * g2_y_zzzzzz[k];

                        g_yz_xxxxxx[k] = g1_z_xxxxxxy[k] - aby * g2_z_xxxxxx[k];

                        g_yz_xxxxxy[k] = g1_z_xxxxxyy[k] - aby * g2_z_xxxxxy[k];

                        g_yz_xxxxxz[k] = g1_z_xxxxxyz[k] - aby * g2_z_xxxxxz[k];

                        g_yz_xxxxyy[k] = g1_z_xxxxyyy[k] - aby * g2_z_xxxxyy[k];

                        g_yz_xxxxyz[k] = g1_z_xxxxyyz[k] - aby * g2_z_xxxxyz[k];

                        g_yz_xxxxzz[k] = g1_z_xxxxyzz[k] - aby * g2_z_xxxxzz[k];

                        g_yz_xxxyyy[k] = g1_z_xxxyyyy[k] - aby * g2_z_xxxyyy[k];

                        g_yz_xxxyyz[k] = g1_z_xxxyyyz[k] - aby * g2_z_xxxyyz[k];

                        g_yz_xxxyzz[k] = g1_z_xxxyyzz[k] - aby * g2_z_xxxyzz[k];

                        g_yz_xxxzzz[k] = g1_z_xxxyzzz[k] - aby * g2_z_xxxzzz[k];

                        g_yz_xxyyyy[k] = g1_z_xxyyyyy[k] - aby * g2_z_xxyyyy[k];

                        g_yz_xxyyyz[k] = g1_z_xxyyyyz[k] - aby * g2_z_xxyyyz[k];

                        g_yz_xxyyzz[k] = g1_z_xxyyyzz[k] - aby * g2_z_xxyyzz[k];

                        g_yz_xxyzzz[k] = g1_z_xxyyzzz[k] - aby * g2_z_xxyzzz[k];

                        g_yz_xxzzzz[k] = g1_z_xxyzzzz[k] - aby * g2_z_xxzzzz[k];

                        g_yz_xyyyyy[k] = g1_z_xyyyyyy[k] - aby * g2_z_xyyyyy[k];

                        g_yz_xyyyyz[k] = g1_z_xyyyyyz[k] - aby * g2_z_xyyyyz[k];

                        g_yz_xyyyzz[k] = g1_z_xyyyyzz[k] - aby * g2_z_xyyyzz[k];

                        g_yz_xyyzzz[k] = g1_z_xyyyzzz[k] - aby * g2_z_xyyzzz[k];

                        g_yz_xyzzzz[k] = g1_z_xyyzzzz[k] - aby * g2_z_xyzzzz[k];

                        g_yz_xzzzzz[k] = g1_z_xyzzzzz[k] - aby * g2_z_xzzzzz[k];

                        g_yz_yyyyyy[k] = g1_z_yyyyyyy[k] - aby * g2_z_yyyyyy[k];

                        g_yz_yyyyyz[k] = g1_z_yyyyyyz[k] - aby * g2_z_yyyyyz[k];

                        g_yz_yyyyzz[k] = g1_z_yyyyyzz[k] - aby * g2_z_yyyyzz[k];

                        g_yz_yyyzzz[k] = g1_z_yyyyzzz[k] - aby * g2_z_yyyzzz[k];

                        g_yz_yyzzzz[k] = g1_z_yyyzzzz[k] - aby * g2_z_yyzzzz[k];

                        g_yz_yzzzzz[k] = g1_z_yyzzzzz[k] - aby * g2_z_yzzzzz[k];

                        g_yz_zzzzzz[k] = g1_z_yzzzzzz[k] - aby * g2_z_zzzzzz[k];

                        // leading z component

                        g_zz_xxxxxx[k] = g1_z_xxxxxxz[k] - abz * g2_z_xxxxxx[k];

                        g_zz_xxxxxy[k] = g1_z_xxxxxyz[k] - abz * g2_z_xxxxxy[k];

                        g_zz_xxxxxz[k] = g1_z_xxxxxzz[k] - abz * g2_z_xxxxxz[k];

                        g_zz_xxxxyy[k] = g1_z_xxxxyyz[k] - abz * g2_z_xxxxyy[k];

                        g_zz_xxxxyz[k] = g1_z_xxxxyzz[k] - abz * g2_z_xxxxyz[k];

                        g_zz_xxxxzz[k] = g1_z_xxxxzzz[k] - abz * g2_z_xxxxzz[k];

                        g_zz_xxxyyy[k] = g1_z_xxxyyyz[k] - abz * g2_z_xxxyyy[k];

                        g_zz_xxxyyz[k] = g1_z_xxxyyzz[k] - abz * g2_z_xxxyyz[k];

                        g_zz_xxxyzz[k] = g1_z_xxxyzzz[k] - abz * g2_z_xxxyzz[k];

                        g_zz_xxxzzz[k] = g1_z_xxxzzzz[k] - abz * g2_z_xxxzzz[k];

                        g_zz_xxyyyy[k] = g1_z_xxyyyyz[k] - abz * g2_z_xxyyyy[k];

                        g_zz_xxyyyz[k] = g1_z_xxyyyzz[k] - abz * g2_z_xxyyyz[k];

                        g_zz_xxyyzz[k] = g1_z_xxyyzzz[k] - abz * g2_z_xxyyzz[k];

                        g_zz_xxyzzz[k] = g1_z_xxyzzzz[k] - abz * g2_z_xxyzzz[k];

                        g_zz_xxzzzz[k] = g1_z_xxzzzzz[k] - abz * g2_z_xxzzzz[k];

                        g_zz_xyyyyy[k] = g1_z_xyyyyyz[k] - abz * g2_z_xyyyyy[k];

                        g_zz_xyyyyz[k] = g1_z_xyyyyzz[k] - abz * g2_z_xyyyyz[k];

                        g_zz_xyyyzz[k] = g1_z_xyyyzzz[k] - abz * g2_z_xyyyzz[k];

                        g_zz_xyyzzz[k] = g1_z_xyyzzzz[k] - abz * g2_z_xyyzzz[k];

                        g_zz_xyzzzz[k] = g1_z_xyzzzzz[k] - abz * g2_z_xyzzzz[k];

                        g_zz_xzzzzz[k] = g1_z_xzzzzzz[k] - abz * g2_z_xzzzzz[k];

                        g_zz_yyyyyy[k] = g1_z_yyyyyyz[k] - abz * g2_z_yyyyyy[k];

                        g_zz_yyyyyz[k] = g1_z_yyyyyzz[k] - abz * g2_z_yyyyyz[k];

                        g_zz_yyyyzz[k] = g1_z_yyyyzzz[k] - abz * g2_z_yyyyzz[k];

                        g_zz_yyyzzz[k] = g1_z_yyyzzzz[k] - abz * g2_z_yyyzzz[k];

                        g_zz_yyzzzz[k] = g1_z_yyzzzzz[k] - abz * g2_z_yyzzzz[k];

                        g_zz_yzzzzz[k] = g1_z_yzzzzzz[k] - abz * g2_z_yzzzzz[k];

                        g_zz_zzzzzz[k] = g1_z_zzzzzzz[k] - abz * g2_z_zzzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForFFXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 3) && (recPattern[i].second() == 3))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (33|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {3, 3, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 4, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 3, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (DF|g(r,r')|XX) integrals

                    auto g2_xx_xxx = braBuffer.data(g2off + j);

                    auto g2_xx_xxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_xx_xxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_xx_xyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_xx_xyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_xx_xzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_xx_yyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_xx_yyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_xx_yzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_xx_zzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_xy_xxx = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_xy_xxy = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_xy_xxz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_xy_xyy = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_xy_xyz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_xy_xzz = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_xy_yyy = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_xy_yyz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_xy_yzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_xy_zzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_xz_xxx = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_xz_xxy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_xz_xxz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_xz_xyy = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_xz_xyz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_xz_xzz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_xz_yyy = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_xz_yyz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_xz_yzz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_xz_zzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_yy_xxx = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_yy_xxy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_yy_xxz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_yy_xyy = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_yy_xyz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_yy_xzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_yy_yyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_yy_yyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_yy_yzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_yy_zzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_yz_xxx = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_yz_xxy = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_yz_xxz = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_yz_xyy = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_yz_xyz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_yz_xzz = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_yz_yyy = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_yz_yyz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_yz_yzz = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_yz_zzz = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_zz_xxx = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_zz_xxy = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_zz_xxz = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_zz_xyy = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_zz_xyz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_zz_xzz = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_zz_yyy = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_zz_yyz = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_zz_yzz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_zz_zzz = braBuffer.data(g2off + 59 * kcomp + j);

                    // set up pointers to (DG|g(r,r')|XX) integrals

                    auto g1_xx_xxxx = braBuffer.data(g1off + j);

                    auto g1_xx_xxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_xx_xxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_xx_xxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_xx_xxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_xx_xxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_xx_xyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_xx_xyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_xx_xyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_xx_xzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_xy_xxxx = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_xy_xxxy = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_xy_xxxz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_xy_xxyy = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_xy_xxyz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_xy_xxzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_xy_xyyy = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_xy_xyyz = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_xy_xyzz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_xy_xzzz = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_xz_xxxx = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_xz_xxxy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_xz_xxxz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_xz_xxyy = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_xz_xxyz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_xz_xxzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_xz_xyyy = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_xz_xyyz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_xz_xyzz = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_xz_xzzz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_yy_xxxx = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_yy_xxxy = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_yy_xxxz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_yy_xxyy = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_yy_xxyz = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_yy_xxzz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_yy_xyyy = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_yy_xyyz = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_yy_xyzz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_yy_xzzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_yy_yyyy = braBuffer.data(g1off + 55 * kcomp + j);

                    auto g1_yy_yyyz = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_yy_yyzz = braBuffer.data(g1off + 57 * kcomp + j);

                    auto g1_yy_yzzz = braBuffer.data(g1off + 58 * kcomp + j);

                    auto g1_yz_xxxx = braBuffer.data(g1off + 60 * kcomp + j);

                    auto g1_yz_xxxy = braBuffer.data(g1off + 61 * kcomp + j);

                    auto g1_yz_xxxz = braBuffer.data(g1off + 62 * kcomp + j);

                    auto g1_yz_xxyy = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_yz_xxyz = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_yz_xxzz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_yz_xyyy = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_yz_xyyz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_yz_xyzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_yz_xzzz = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_yz_yyyy = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_yz_yyyz = braBuffer.data(g1off + 71 * kcomp + j);

                    auto g1_yz_yyzz = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_yz_yzzz = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_zz_xxxx = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_zz_xxxy = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_zz_xxxz = braBuffer.data(g1off + 77 * kcomp + j);

                    auto g1_zz_xxyy = braBuffer.data(g1off + 78 * kcomp + j);

                    auto g1_zz_xxyz = braBuffer.data(g1off + 79 * kcomp + j);

                    auto g1_zz_xxzz = braBuffer.data(g1off + 80 * kcomp + j);

                    auto g1_zz_xyyy = braBuffer.data(g1off + 81 * kcomp + j);

                    auto g1_zz_xyyz = braBuffer.data(g1off + 82 * kcomp + j);

                    auto g1_zz_xyzz = braBuffer.data(g1off + 83 * kcomp + j);

                    auto g1_zz_xzzz = braBuffer.data(g1off + 84 * kcomp + j);

                    auto g1_zz_yyyy = braBuffer.data(g1off + 85 * kcomp + j);

                    auto g1_zz_yyyz = braBuffer.data(g1off + 86 * kcomp + j);

                    auto g1_zz_yyzz = braBuffer.data(g1off + 87 * kcomp + j);

                    auto g1_zz_yzzz = braBuffer.data(g1off + 88 * kcomp + j);

                    auto g1_zz_zzzz = braBuffer.data(g1off + 89 * kcomp + j);

                    // set up pointers to (FF|g(r,r')|XX) integrals

                    auto g_xxx_xxx = braBuffer.data(goff + j);

                    auto g_xxx_xxy = braBuffer.data(goff + kcomp + j);

                    auto g_xxx_xxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xxx_xyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xxx_xyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xxx_xzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xxx_yyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xxx_yyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xxx_yzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xxx_zzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xxy_xxx = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xxy_xxy = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xxy_xxz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xxy_xyy = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xxy_xyz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xxy_xzz = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xxy_yyy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xxy_yyz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xxy_yzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xxy_zzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xxz_xxx = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xxz_xxy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xxz_xxz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xxz_xyy = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xxz_xyz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xxz_xzz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xxz_yyy = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xxz_yyz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xxz_yzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xxz_zzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xyy_xxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xyy_xxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xyy_xxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xyy_xyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xyy_xyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xyy_xzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xyy_yyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xyy_yyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xyy_yzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xyy_zzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xyz_xxx = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xyz_xxy = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xyz_xxz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xyz_xyy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xyz_xyz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xyz_xzz = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xyz_yyy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xyz_yyz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xyz_yzz = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xyz_zzz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xzz_xxx = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xzz_xxy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xzz_xxz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xzz_xyy = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xzz_xyz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xzz_xzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xzz_yyy = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xzz_yyz = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xzz_yzz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xzz_zzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_yyy_xxx = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_yyy_xxy = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_yyy_xxz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_yyy_xyy = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_yyy_xyz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_yyy_xzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_yyy_yyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_yyy_yyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_yyy_yzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_yyy_zzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_yyz_xxx = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_yyz_xxy = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_yyz_xxz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_yyz_xyy = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_yyz_xyz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_yyz_xzz = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_yyz_yyy = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_yyz_yyz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_yyz_yzz = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_yyz_zzz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_yzz_xxx = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_yzz_xxy = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_yzz_xxz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_yzz_xyy = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_yzz_xyz = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_yzz_xzz = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_yzz_yyy = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_yzz_yyz = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_yzz_yzz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_yzz_zzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_zzz_xxx = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_zzz_xxy = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_zzz_xxz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_zzz_xyy = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_zzz_xyz = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_zzz_xzz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_zzz_yyy = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_zzz_yyz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_zzz_yzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_zzz_zzz = braBuffer.data(goff + 99 * kcomp + j);

                    #pragma omp simd aligned(g2_xx_xxx, g2_xx_xxy, g2_xx_xxz, g2_xx_xyy,\
                                             g2_xx_xyz, g2_xx_xzz, g2_xx_yyy, g2_xx_yyz,\
                                             g2_xx_yzz, g2_xx_zzz, g2_xy_xxx, g2_xy_xxy,\
                                             g2_xy_xxz, g2_xy_xyy, g2_xy_xyz, g2_xy_xzz,\
                                             g2_xy_yyy, g2_xy_yyz, g2_xy_yzz, g2_xy_zzz,\
                                             g2_xz_xxx, g2_xz_xxy, g2_xz_xxz, g2_xz_xyy,\
                                             g2_xz_xyz, g2_xz_xzz, g2_xz_yyy, g2_xz_yyz,\
                                             g2_xz_yzz, g2_xz_zzz, g2_yy_xxx, g2_yy_xxy,\
                                             g2_yy_xxz, g2_yy_xyy, g2_yy_xyz, g2_yy_xzz,\
                                             g2_yy_yyy, g2_yy_yyz, g2_yy_yzz, g2_yy_zzz,\
                                             g2_yz_xxx, g2_yz_xxy, g2_yz_xxz, g2_yz_xyy,\
                                             g2_yz_xyz, g2_yz_xzz, g2_yz_yyy, g2_yz_yyz,\
                                             g2_yz_yzz, g2_yz_zzz, g2_zz_xxx, g2_zz_xxy,\
                                             g2_zz_xxz, g2_zz_xyy, g2_zz_xyz, g2_zz_xzz,\
                                             g2_zz_yyy, g2_zz_yyz, g2_zz_yzz, g2_zz_zzz,\
                                             g1_xx_xxxx, g1_xx_xxxy, g1_xx_xxxz,\
                                             g1_xx_xxyy, g1_xx_xxyz, g1_xx_xxzz,\
                                             g1_xx_xyyy, g1_xx_xyyz, g1_xx_xyzz,\
                                             g1_xx_xzzz, \
                                             g1_xy_xxxx, g1_xy_xxxy, g1_xy_xxxz,\
                                             g1_xy_xxyy, g1_xy_xxyz, g1_xy_xxzz,\
                                             g1_xy_xyyy, g1_xy_xyyz, g1_xy_xyzz,\
                                             g1_xy_xzzz, \
                                             g1_xz_xxxx, g1_xz_xxxy, g1_xz_xxxz,\
                                             g1_xz_xxyy, g1_xz_xxyz, g1_xz_xxzz,\
                                             g1_xz_xyyy, g1_xz_xyyz, g1_xz_xyzz,\
                                             g1_xz_xzzz, \
                                             g1_yy_xxxx, g1_yy_xxxy, g1_yy_xxxz,\
                                             g1_yy_xxyy, g1_yy_xxyz, g1_yy_xxzz,\
                                             g1_yy_xyyy, g1_yy_xyyz, g1_yy_xyzz,\
                                             g1_yy_xzzz, g1_yy_yyyy, g1_yy_yyyz,\
                                             g1_yy_yyzz, g1_yy_yzzz, \
                                             g1_yz_xxxx, g1_yz_xxxy, g1_yz_xxxz,\
                                             g1_yz_xxyy, g1_yz_xxyz, g1_yz_xxzz,\
                                             g1_yz_xyyy, g1_yz_xyyz, g1_yz_xyzz,\
                                             g1_yz_xzzz, g1_yz_yyyy, g1_yz_yyyz,\
                                             g1_yz_yyzz, g1_yz_yzzz, \
                                             g1_zz_xxxx, g1_zz_xxxy, g1_zz_xxxz,\
                                             g1_zz_xxyy, g1_zz_xxyz, g1_zz_xxzz,\
                                             g1_zz_xyyy, g1_zz_xyyz, g1_zz_xyzz,\
                                             g1_zz_xzzz, g1_zz_yyyy, g1_zz_yyyz,\
                                             g1_zz_yyzz, g1_zz_yzzz, g1_zz_zzzz,\
                                             g_xxx_xxx, g_xxx_xxy, g_xxx_xxz, g_xxx_xyy,\
                                             g_xxx_xyz, g_xxx_xzz, g_xxx_yyy, g_xxx_yyz,\
                                             g_xxx_yzz, g_xxx_zzz, g_xxy_xxx, g_xxy_xxy,\
                                             g_xxy_xxz, g_xxy_xyy, g_xxy_xyz, g_xxy_xzz,\
                                             g_xxy_yyy, g_xxy_yyz, g_xxy_yzz, g_xxy_zzz,\
                                             g_xxz_xxx, g_xxz_xxy, g_xxz_xxz, g_xxz_xyy,\
                                             g_xxz_xyz, g_xxz_xzz, g_xxz_yyy, g_xxz_yyz,\
                                             g_xxz_yzz, g_xxz_zzz, g_xyy_xxx, g_xyy_xxy,\
                                             g_xyy_xxz, g_xyy_xyy, g_xyy_xyz, g_xyy_xzz,\
                                             g_xyy_yyy, g_xyy_yyz, g_xyy_yzz, g_xyy_zzz,\
                                             g_xyz_xxx, g_xyz_xxy, g_xyz_xxz, g_xyz_xyy,\
                                             g_xyz_xyz, g_xyz_xzz, g_xyz_yyy, g_xyz_yyz,\
                                             g_xyz_yzz, g_xyz_zzz, g_xzz_xxx, g_xzz_xxy,\
                                             g_xzz_xxz, g_xzz_xyy, g_xzz_xyz, g_xzz_xzz,\
                                             g_xzz_yyy, g_xzz_yyz, g_xzz_yzz, g_xzz_zzz,\
                                             g_yyy_xxx, g_yyy_xxy, g_yyy_xxz, g_yyy_xyy,\
                                             g_yyy_xyz, g_yyy_xzz, g_yyy_yyy, g_yyy_yyz,\
                                             g_yyy_yzz, g_yyy_zzz, g_yyz_xxx, g_yyz_xxy,\
                                             g_yyz_xxz, g_yyz_xyy, g_yyz_xyz, g_yyz_xzz,\
                                             g_yyz_yyy, g_yyz_yyz, g_yyz_yzz, g_yyz_zzz,\
                                             g_yzz_xxx, g_yzz_xxy, g_yzz_xxz, g_yzz_xyy,\
                                             g_yzz_xyz, g_yzz_xzz, g_yzz_yyy, g_yzz_yyz,\
                                             g_yzz_yzz, g_yzz_zzz, g_zzz_xxx, g_zzz_xxy,\
                                             g_zzz_xxz, g_zzz_xyy, g_zzz_xyz, g_zzz_xzz,\
                                             g_zzz_yyy, g_zzz_yyz, g_zzz_yzz, g_zzz_zzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xxx_xxx[k] = g1_xx_xxxx[k] - abx * g2_xx_xxx[k];

                        g_xxx_xxy[k] = g1_xx_xxxy[k] - abx * g2_xx_xxy[k];

                        g_xxx_xxz[k] = g1_xx_xxxz[k] - abx * g2_xx_xxz[k];

                        g_xxx_xyy[k] = g1_xx_xxyy[k] - abx * g2_xx_xyy[k];

                        g_xxx_xyz[k] = g1_xx_xxyz[k] - abx * g2_xx_xyz[k];

                        g_xxx_xzz[k] = g1_xx_xxzz[k] - abx * g2_xx_xzz[k];

                        g_xxx_yyy[k] = g1_xx_xyyy[k] - abx * g2_xx_yyy[k];

                        g_xxx_yyz[k] = g1_xx_xyyz[k] - abx * g2_xx_yyz[k];

                        g_xxx_yzz[k] = g1_xx_xyzz[k] - abx * g2_xx_yzz[k];

                        g_xxx_zzz[k] = g1_xx_xzzz[k] - abx * g2_xx_zzz[k];

                        g_xxy_xxx[k] = g1_xy_xxxx[k] - abx * g2_xy_xxx[k];

                        g_xxy_xxy[k] = g1_xy_xxxy[k] - abx * g2_xy_xxy[k];

                        g_xxy_xxz[k] = g1_xy_xxxz[k] - abx * g2_xy_xxz[k];

                        g_xxy_xyy[k] = g1_xy_xxyy[k] - abx * g2_xy_xyy[k];

                        g_xxy_xyz[k] = g1_xy_xxyz[k] - abx * g2_xy_xyz[k];

                        g_xxy_xzz[k] = g1_xy_xxzz[k] - abx * g2_xy_xzz[k];

                        g_xxy_yyy[k] = g1_xy_xyyy[k] - abx * g2_xy_yyy[k];

                        g_xxy_yyz[k] = g1_xy_xyyz[k] - abx * g2_xy_yyz[k];

                        g_xxy_yzz[k] = g1_xy_xyzz[k] - abx * g2_xy_yzz[k];

                        g_xxy_zzz[k] = g1_xy_xzzz[k] - abx * g2_xy_zzz[k];

                        g_xxz_xxx[k] = g1_xz_xxxx[k] - abx * g2_xz_xxx[k];

                        g_xxz_xxy[k] = g1_xz_xxxy[k] - abx * g2_xz_xxy[k];

                        g_xxz_xxz[k] = g1_xz_xxxz[k] - abx * g2_xz_xxz[k];

                        g_xxz_xyy[k] = g1_xz_xxyy[k] - abx * g2_xz_xyy[k];

                        g_xxz_xyz[k] = g1_xz_xxyz[k] - abx * g2_xz_xyz[k];

                        g_xxz_xzz[k] = g1_xz_xxzz[k] - abx * g2_xz_xzz[k];

                        g_xxz_yyy[k] = g1_xz_xyyy[k] - abx * g2_xz_yyy[k];

                        g_xxz_yyz[k] = g1_xz_xyyz[k] - abx * g2_xz_yyz[k];

                        g_xxz_yzz[k] = g1_xz_xyzz[k] - abx * g2_xz_yzz[k];

                        g_xxz_zzz[k] = g1_xz_xzzz[k] - abx * g2_xz_zzz[k];

                        g_xyy_xxx[k] = g1_yy_xxxx[k] - abx * g2_yy_xxx[k];

                        g_xyy_xxy[k] = g1_yy_xxxy[k] - abx * g2_yy_xxy[k];

                        g_xyy_xxz[k] = g1_yy_xxxz[k] - abx * g2_yy_xxz[k];

                        g_xyy_xyy[k] = g1_yy_xxyy[k] - abx * g2_yy_xyy[k];

                        g_xyy_xyz[k] = g1_yy_xxyz[k] - abx * g2_yy_xyz[k];

                        g_xyy_xzz[k] = g1_yy_xxzz[k] - abx * g2_yy_xzz[k];

                        g_xyy_yyy[k] = g1_yy_xyyy[k] - abx * g2_yy_yyy[k];

                        g_xyy_yyz[k] = g1_yy_xyyz[k] - abx * g2_yy_yyz[k];

                        g_xyy_yzz[k] = g1_yy_xyzz[k] - abx * g2_yy_yzz[k];

                        g_xyy_zzz[k] = g1_yy_xzzz[k] - abx * g2_yy_zzz[k];

                        g_xyz_xxx[k] = g1_yz_xxxx[k] - abx * g2_yz_xxx[k];

                        g_xyz_xxy[k] = g1_yz_xxxy[k] - abx * g2_yz_xxy[k];

                        g_xyz_xxz[k] = g1_yz_xxxz[k] - abx * g2_yz_xxz[k];

                        g_xyz_xyy[k] = g1_yz_xxyy[k] - abx * g2_yz_xyy[k];

                        g_xyz_xyz[k] = g1_yz_xxyz[k] - abx * g2_yz_xyz[k];

                        g_xyz_xzz[k] = g1_yz_xxzz[k] - abx * g2_yz_xzz[k];

                        g_xyz_yyy[k] = g1_yz_xyyy[k] - abx * g2_yz_yyy[k];

                        g_xyz_yyz[k] = g1_yz_xyyz[k] - abx * g2_yz_yyz[k];

                        g_xyz_yzz[k] = g1_yz_xyzz[k] - abx * g2_yz_yzz[k];

                        g_xyz_zzz[k] = g1_yz_xzzz[k] - abx * g2_yz_zzz[k];

                        g_xzz_xxx[k] = g1_zz_xxxx[k] - abx * g2_zz_xxx[k];

                        g_xzz_xxy[k] = g1_zz_xxxy[k] - abx * g2_zz_xxy[k];

                        g_xzz_xxz[k] = g1_zz_xxxz[k] - abx * g2_zz_xxz[k];

                        g_xzz_xyy[k] = g1_zz_xxyy[k] - abx * g2_zz_xyy[k];

                        g_xzz_xyz[k] = g1_zz_xxyz[k] - abx * g2_zz_xyz[k];

                        g_xzz_xzz[k] = g1_zz_xxzz[k] - abx * g2_zz_xzz[k];

                        g_xzz_yyy[k] = g1_zz_xyyy[k] - abx * g2_zz_yyy[k];

                        g_xzz_yyz[k] = g1_zz_xyyz[k] - abx * g2_zz_yyz[k];

                        g_xzz_yzz[k] = g1_zz_xyzz[k] - abx * g2_zz_yzz[k];

                        g_xzz_zzz[k] = g1_zz_xzzz[k] - abx * g2_zz_zzz[k];

                        // leading y component

                        g_yyy_xxx[k] = g1_yy_xxxy[k] - aby * g2_yy_xxx[k];

                        g_yyy_xxy[k] = g1_yy_xxyy[k] - aby * g2_yy_xxy[k];

                        g_yyy_xxz[k] = g1_yy_xxyz[k] - aby * g2_yy_xxz[k];

                        g_yyy_xyy[k] = g1_yy_xyyy[k] - aby * g2_yy_xyy[k];

                        g_yyy_xyz[k] = g1_yy_xyyz[k] - aby * g2_yy_xyz[k];

                        g_yyy_xzz[k] = g1_yy_xyzz[k] - aby * g2_yy_xzz[k];

                        g_yyy_yyy[k] = g1_yy_yyyy[k] - aby * g2_yy_yyy[k];

                        g_yyy_yyz[k] = g1_yy_yyyz[k] - aby * g2_yy_yyz[k];

                        g_yyy_yzz[k] = g1_yy_yyzz[k] - aby * g2_yy_yzz[k];

                        g_yyy_zzz[k] = g1_yy_yzzz[k] - aby * g2_yy_zzz[k];

                        g_yyz_xxx[k] = g1_yz_xxxy[k] - aby * g2_yz_xxx[k];

                        g_yyz_xxy[k] = g1_yz_xxyy[k] - aby * g2_yz_xxy[k];

                        g_yyz_xxz[k] = g1_yz_xxyz[k] - aby * g2_yz_xxz[k];

                        g_yyz_xyy[k] = g1_yz_xyyy[k] - aby * g2_yz_xyy[k];

                        g_yyz_xyz[k] = g1_yz_xyyz[k] - aby * g2_yz_xyz[k];

                        g_yyz_xzz[k] = g1_yz_xyzz[k] - aby * g2_yz_xzz[k];

                        g_yyz_yyy[k] = g1_yz_yyyy[k] - aby * g2_yz_yyy[k];

                        g_yyz_yyz[k] = g1_yz_yyyz[k] - aby * g2_yz_yyz[k];

                        g_yyz_yzz[k] = g1_yz_yyzz[k] - aby * g2_yz_yzz[k];

                        g_yyz_zzz[k] = g1_yz_yzzz[k] - aby * g2_yz_zzz[k];

                        g_yzz_xxx[k] = g1_zz_xxxy[k] - aby * g2_zz_xxx[k];

                        g_yzz_xxy[k] = g1_zz_xxyy[k] - aby * g2_zz_xxy[k];

                        g_yzz_xxz[k] = g1_zz_xxyz[k] - aby * g2_zz_xxz[k];

                        g_yzz_xyy[k] = g1_zz_xyyy[k] - aby * g2_zz_xyy[k];

                        g_yzz_xyz[k] = g1_zz_xyyz[k] - aby * g2_zz_xyz[k];

                        g_yzz_xzz[k] = g1_zz_xyzz[k] - aby * g2_zz_xzz[k];

                        g_yzz_yyy[k] = g1_zz_yyyy[k] - aby * g2_zz_yyy[k];

                        g_yzz_yyz[k] = g1_zz_yyyz[k] - aby * g2_zz_yyz[k];

                        g_yzz_yzz[k] = g1_zz_yyzz[k] - aby * g2_zz_yzz[k];

                        g_yzz_zzz[k] = g1_zz_yzzz[k] - aby * g2_zz_zzz[k];

                        // leading z component

                        g_zzz_xxx[k] = g1_zz_xxxz[k] - abz * g2_zz_xxx[k];

                        g_zzz_xxy[k] = g1_zz_xxyz[k] - abz * g2_zz_xxy[k];

                        g_zzz_xxz[k] = g1_zz_xxzz[k] - abz * g2_zz_xxz[k];

                        g_zzz_xyy[k] = g1_zz_xyyz[k] - abz * g2_zz_xyy[k];

                        g_zzz_xyz[k] = g1_zz_xyzz[k] - abz * g2_zz_xyz[k];

                        g_zzz_xzz[k] = g1_zz_xzzz[k] - abz * g2_zz_xzz[k];

                        g_zzz_yyy[k] = g1_zz_yyyz[k] - abz * g2_zz_yyy[k];

                        g_zzz_yyz[k] = g1_zz_yyzz[k] - abz * g2_zz_yyz[k];

                        g_zzz_yzz[k] = g1_zz_yzzz[k] - abz * g2_zz_yzz[k];

                        g_zzz_zzz[k] = g1_zz_zzzz[k] - abz * g2_zz_zzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForFGXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 3) && (recPattern[i].second() == 4))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (34|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {3, 4, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 5, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 4, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (DG|g(r,r')|XX) integrals

                    auto g2_xx_xxxx = braBuffer.data(g2off + j);

                    auto g2_xx_xxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_xx_xxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_xx_xxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_xx_xxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_xx_xxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_xx_xyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_xx_xyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_xx_xyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_xx_xzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_xx_yyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_xx_yyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_xx_yyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_xx_yzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_xx_zzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_xy_xxxx = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_xy_xxxy = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_xy_xxxz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_xy_xxyy = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_xy_xxyz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_xy_xxzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_xy_xyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_xy_xyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_xy_xyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_xy_xzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_xy_yyyy = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_xy_yyyz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_xy_yyzz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_xy_yzzz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_xy_zzzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_xz_xxxx = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_xz_xxxy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_xz_xxxz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_xz_xxyy = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_xz_xxyz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_xz_xxzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_xz_xyyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_xz_xyyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_xz_xyzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_xz_xzzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_xz_yyyy = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_xz_yyyz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_xz_yyzz = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_xz_yzzz = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_xz_zzzz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_yy_xxxx = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_yy_xxxy = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_yy_xxxz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_yy_xxyy = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_yy_xxyz = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_yy_xxzz = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_yy_xyyy = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_yy_xyyz = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_yy_xyzz = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_yy_xzzz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_yy_yyyy = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_yy_yyyz = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_yy_yyzz = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_yy_yzzz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_yy_zzzz = braBuffer.data(g2off + 59 * kcomp + j);

                    auto g2_yz_xxxx = braBuffer.data(g2off + 60 * kcomp + j);

                    auto g2_yz_xxxy = braBuffer.data(g2off + 61 * kcomp + j);

                    auto g2_yz_xxxz = braBuffer.data(g2off + 62 * kcomp + j);

                    auto g2_yz_xxyy = braBuffer.data(g2off + 63 * kcomp + j);

                    auto g2_yz_xxyz = braBuffer.data(g2off + 64 * kcomp + j);

                    auto g2_yz_xxzz = braBuffer.data(g2off + 65 * kcomp + j);

                    auto g2_yz_xyyy = braBuffer.data(g2off + 66 * kcomp + j);

                    auto g2_yz_xyyz = braBuffer.data(g2off + 67 * kcomp + j);

                    auto g2_yz_xyzz = braBuffer.data(g2off + 68 * kcomp + j);

                    auto g2_yz_xzzz = braBuffer.data(g2off + 69 * kcomp + j);

                    auto g2_yz_yyyy = braBuffer.data(g2off + 70 * kcomp + j);

                    auto g2_yz_yyyz = braBuffer.data(g2off + 71 * kcomp + j);

                    auto g2_yz_yyzz = braBuffer.data(g2off + 72 * kcomp + j);

                    auto g2_yz_yzzz = braBuffer.data(g2off + 73 * kcomp + j);

                    auto g2_yz_zzzz = braBuffer.data(g2off + 74 * kcomp + j);

                    auto g2_zz_xxxx = braBuffer.data(g2off + 75 * kcomp + j);

                    auto g2_zz_xxxy = braBuffer.data(g2off + 76 * kcomp + j);

                    auto g2_zz_xxxz = braBuffer.data(g2off + 77 * kcomp + j);

                    auto g2_zz_xxyy = braBuffer.data(g2off + 78 * kcomp + j);

                    auto g2_zz_xxyz = braBuffer.data(g2off + 79 * kcomp + j);

                    auto g2_zz_xxzz = braBuffer.data(g2off + 80 * kcomp + j);

                    auto g2_zz_xyyy = braBuffer.data(g2off + 81 * kcomp + j);

                    auto g2_zz_xyyz = braBuffer.data(g2off + 82 * kcomp + j);

                    auto g2_zz_xyzz = braBuffer.data(g2off + 83 * kcomp + j);

                    auto g2_zz_xzzz = braBuffer.data(g2off + 84 * kcomp + j);

                    auto g2_zz_yyyy = braBuffer.data(g2off + 85 * kcomp + j);

                    auto g2_zz_yyyz = braBuffer.data(g2off + 86 * kcomp + j);

                    auto g2_zz_yyzz = braBuffer.data(g2off + 87 * kcomp + j);

                    auto g2_zz_yzzz = braBuffer.data(g2off + 88 * kcomp + j);

                    auto g2_zz_zzzz = braBuffer.data(g2off + 89 * kcomp + j);

                    // set up pointers to (DH|g(r,r')|XX) integrals

                    auto g1_xx_xxxxx = braBuffer.data(g1off + j);

                    auto g1_xx_xxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_xx_xxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_xx_xxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_xx_xxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_xx_xxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_xx_xxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_xx_xxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_xx_xxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_xx_xxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_xx_xyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_xx_xyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_xx_xyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_xx_xyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_xx_xzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_xy_xxxxx = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_xy_xxxxy = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_xy_xxxxz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_xy_xxxyy = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_xy_xxxyz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_xy_xxxzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_xy_xxyyy = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_xy_xxyyz = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_xy_xxyzz = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_xy_xxzzz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_xy_xyyyy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_xy_xyyyz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_xy_xyyzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_xy_xyzzz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_xy_xzzzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_xz_xxxxx = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_xz_xxxxy = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_xz_xxxxz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_xz_xxxyy = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_xz_xxxyz = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_xz_xxxzz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_xz_xxyyy = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_xz_xxyyz = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_xz_xxyzz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_xz_xxzzz = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_xz_xyyyy = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_xz_xyyyz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_xz_xyyzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_xz_xyzzz = braBuffer.data(g1off + 55 * kcomp + j);

                    auto g1_xz_xzzzz = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_yy_xxxxx = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_yy_xxxxy = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_yy_xxxxz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_yy_xxxyy = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_yy_xxxyz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_yy_xxxzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_yy_xxyyy = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_yy_xxyyz = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_yy_xxyzz = braBuffer.data(g1off + 71 * kcomp + j);

                    auto g1_yy_xxzzz = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_yy_xyyyy = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_yy_xyyyz = braBuffer.data(g1off + 74 * kcomp + j);

                    auto g1_yy_xyyzz = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_yy_xyzzz = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_yy_xzzzz = braBuffer.data(g1off + 77 * kcomp + j);

                    auto g1_yy_yyyyy = braBuffer.data(g1off + 78 * kcomp + j);

                    auto g1_yy_yyyyz = braBuffer.data(g1off + 79 * kcomp + j);

                    auto g1_yy_yyyzz = braBuffer.data(g1off + 80 * kcomp + j);

                    auto g1_yy_yyzzz = braBuffer.data(g1off + 81 * kcomp + j);

                    auto g1_yy_yzzzz = braBuffer.data(g1off + 82 * kcomp + j);

                    auto g1_yz_xxxxx = braBuffer.data(g1off + 84 * kcomp + j);

                    auto g1_yz_xxxxy = braBuffer.data(g1off + 85 * kcomp + j);

                    auto g1_yz_xxxxz = braBuffer.data(g1off + 86 * kcomp + j);

                    auto g1_yz_xxxyy = braBuffer.data(g1off + 87 * kcomp + j);

                    auto g1_yz_xxxyz = braBuffer.data(g1off + 88 * kcomp + j);

                    auto g1_yz_xxxzz = braBuffer.data(g1off + 89 * kcomp + j);

                    auto g1_yz_xxyyy = braBuffer.data(g1off + 90 * kcomp + j);

                    auto g1_yz_xxyyz = braBuffer.data(g1off + 91 * kcomp + j);

                    auto g1_yz_xxyzz = braBuffer.data(g1off + 92 * kcomp + j);

                    auto g1_yz_xxzzz = braBuffer.data(g1off + 93 * kcomp + j);

                    auto g1_yz_xyyyy = braBuffer.data(g1off + 94 * kcomp + j);

                    auto g1_yz_xyyyz = braBuffer.data(g1off + 95 * kcomp + j);

                    auto g1_yz_xyyzz = braBuffer.data(g1off + 96 * kcomp + j);

                    auto g1_yz_xyzzz = braBuffer.data(g1off + 97 * kcomp + j);

                    auto g1_yz_xzzzz = braBuffer.data(g1off + 98 * kcomp + j);

                    auto g1_yz_yyyyy = braBuffer.data(g1off + 99 * kcomp + j);

                    auto g1_yz_yyyyz = braBuffer.data(g1off + 100 * kcomp + j);

                    auto g1_yz_yyyzz = braBuffer.data(g1off + 101 * kcomp + j);

                    auto g1_yz_yyzzz = braBuffer.data(g1off + 102 * kcomp + j);

                    auto g1_yz_yzzzz = braBuffer.data(g1off + 103 * kcomp + j);

                    auto g1_zz_xxxxx = braBuffer.data(g1off + 105 * kcomp + j);

                    auto g1_zz_xxxxy = braBuffer.data(g1off + 106 * kcomp + j);

                    auto g1_zz_xxxxz = braBuffer.data(g1off + 107 * kcomp + j);

                    auto g1_zz_xxxyy = braBuffer.data(g1off + 108 * kcomp + j);

                    auto g1_zz_xxxyz = braBuffer.data(g1off + 109 * kcomp + j);

                    auto g1_zz_xxxzz = braBuffer.data(g1off + 110 * kcomp + j);

                    auto g1_zz_xxyyy = braBuffer.data(g1off + 111 * kcomp + j);

                    auto g1_zz_xxyyz = braBuffer.data(g1off + 112 * kcomp + j);

                    auto g1_zz_xxyzz = braBuffer.data(g1off + 113 * kcomp + j);

                    auto g1_zz_xxzzz = braBuffer.data(g1off + 114 * kcomp + j);

                    auto g1_zz_xyyyy = braBuffer.data(g1off + 115 * kcomp + j);

                    auto g1_zz_xyyyz = braBuffer.data(g1off + 116 * kcomp + j);

                    auto g1_zz_xyyzz = braBuffer.data(g1off + 117 * kcomp + j);

                    auto g1_zz_xyzzz = braBuffer.data(g1off + 118 * kcomp + j);

                    auto g1_zz_xzzzz = braBuffer.data(g1off + 119 * kcomp + j);

                    auto g1_zz_yyyyy = braBuffer.data(g1off + 120 * kcomp + j);

                    auto g1_zz_yyyyz = braBuffer.data(g1off + 121 * kcomp + j);

                    auto g1_zz_yyyzz = braBuffer.data(g1off + 122 * kcomp + j);

                    auto g1_zz_yyzzz = braBuffer.data(g1off + 123 * kcomp + j);

                    auto g1_zz_yzzzz = braBuffer.data(g1off + 124 * kcomp + j);

                    auto g1_zz_zzzzz = braBuffer.data(g1off + 125 * kcomp + j);

                    // set up pointers to (FG|g(r,r')|XX) integrals

                    auto g_xxx_xxxx = braBuffer.data(goff + j);

                    auto g_xxx_xxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xxx_xxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xxx_xxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xxx_xxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xxx_xxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xxx_xyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xxx_xyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xxx_xyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xxx_xzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xxx_yyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xxx_yyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xxx_yyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xxx_yzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xxx_zzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xxy_xxxx = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xxy_xxxy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xxy_xxxz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xxy_xxyy = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xxy_xxyz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xxy_xxzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xxy_xyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xxy_xyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xxy_xyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xxy_xzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xxy_yyyy = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xxy_yyyz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xxy_yyzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xxy_yzzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xxy_zzzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xxz_xxxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xxz_xxxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xxz_xxxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xxz_xxyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xxz_xxyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xxz_xxzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xxz_xyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xxz_xyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xxz_xyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xxz_xzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xxz_yyyy = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xxz_yyyz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xxz_yyzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xxz_yzzz = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xxz_zzzz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xyy_xxxx = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xyy_xxxy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xyy_xxxz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xyy_xxyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xyy_xxyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xyy_xxzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xyy_xyyy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xyy_xyyz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xyy_xyzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xyy_xzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xyy_yyyy = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xyy_yyyz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xyy_yyzz = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xyy_yzzz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xyy_zzzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_xyz_xxxx = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_xyz_xxxy = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_xyz_xxxz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_xyz_xxyy = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_xyz_xxyz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_xyz_xxzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_xyz_xyyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_xyz_xyyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_xyz_xyzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_xyz_xzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_xyz_yyyy = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_xyz_yyyz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_xyz_yyzz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_xyz_yzzz = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_xyz_zzzz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_xzz_xxxx = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_xzz_xxxy = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_xzz_xxxz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_xzz_xxyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_xzz_xxyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_xzz_xxzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_xzz_xyyy = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_xzz_xyyz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_xzz_xyzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_xzz_xzzz = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_xzz_yyyy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_xzz_yyyz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_xzz_yyzz = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_xzz_yzzz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_xzz_zzzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_yyy_xxxx = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_yyy_xxxy = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_yyy_xxxz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_yyy_xxyy = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_yyy_xxyz = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_yyy_xxzz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_yyy_xyyy = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_yyy_xyyz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_yyy_xyzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_yyy_xzzz = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_yyy_yyyy = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_yyy_yyyz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_yyy_yyzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_yyy_yzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_yyy_zzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_yyz_xxxx = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_yyz_xxxy = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_yyz_xxxz = braBuffer.data(goff + 107 * kcomp + j);

                    auto g_yyz_xxyy = braBuffer.data(goff + 108 * kcomp + j);

                    auto g_yyz_xxyz = braBuffer.data(goff + 109 * kcomp + j);

                    auto g_yyz_xxzz = braBuffer.data(goff + 110 * kcomp + j);

                    auto g_yyz_xyyy = braBuffer.data(goff + 111 * kcomp + j);

                    auto g_yyz_xyyz = braBuffer.data(goff + 112 * kcomp + j);

                    auto g_yyz_xyzz = braBuffer.data(goff + 113 * kcomp + j);

                    auto g_yyz_xzzz = braBuffer.data(goff + 114 * kcomp + j);

                    auto g_yyz_yyyy = braBuffer.data(goff + 115 * kcomp + j);

                    auto g_yyz_yyyz = braBuffer.data(goff + 116 * kcomp + j);

                    auto g_yyz_yyzz = braBuffer.data(goff + 117 * kcomp + j);

                    auto g_yyz_yzzz = braBuffer.data(goff + 118 * kcomp + j);

                    auto g_yyz_zzzz = braBuffer.data(goff + 119 * kcomp + j);

                    auto g_yzz_xxxx = braBuffer.data(goff + 120 * kcomp + j);

                    auto g_yzz_xxxy = braBuffer.data(goff + 121 * kcomp + j);

                    auto g_yzz_xxxz = braBuffer.data(goff + 122 * kcomp + j);

                    auto g_yzz_xxyy = braBuffer.data(goff + 123 * kcomp + j);

                    auto g_yzz_xxyz = braBuffer.data(goff + 124 * kcomp + j);

                    auto g_yzz_xxzz = braBuffer.data(goff + 125 * kcomp + j);

                    auto g_yzz_xyyy = braBuffer.data(goff + 126 * kcomp + j);

                    auto g_yzz_xyyz = braBuffer.data(goff + 127 * kcomp + j);

                    auto g_yzz_xyzz = braBuffer.data(goff + 128 * kcomp + j);

                    auto g_yzz_xzzz = braBuffer.data(goff + 129 * kcomp + j);

                    auto g_yzz_yyyy = braBuffer.data(goff + 130 * kcomp + j);

                    auto g_yzz_yyyz = braBuffer.data(goff + 131 * kcomp + j);

                    auto g_yzz_yyzz = braBuffer.data(goff + 132 * kcomp + j);

                    auto g_yzz_yzzz = braBuffer.data(goff + 133 * kcomp + j);

                    auto g_yzz_zzzz = braBuffer.data(goff + 134 * kcomp + j);

                    auto g_zzz_xxxx = braBuffer.data(goff + 135 * kcomp + j);

                    auto g_zzz_xxxy = braBuffer.data(goff + 136 * kcomp + j);

                    auto g_zzz_xxxz = braBuffer.data(goff + 137 * kcomp + j);

                    auto g_zzz_xxyy = braBuffer.data(goff + 138 * kcomp + j);

                    auto g_zzz_xxyz = braBuffer.data(goff + 139 * kcomp + j);

                    auto g_zzz_xxzz = braBuffer.data(goff + 140 * kcomp + j);

                    auto g_zzz_xyyy = braBuffer.data(goff + 141 * kcomp + j);

                    auto g_zzz_xyyz = braBuffer.data(goff + 142 * kcomp + j);

                    auto g_zzz_xyzz = braBuffer.data(goff + 143 * kcomp + j);

                    auto g_zzz_xzzz = braBuffer.data(goff + 144 * kcomp + j);

                    auto g_zzz_yyyy = braBuffer.data(goff + 145 * kcomp + j);

                    auto g_zzz_yyyz = braBuffer.data(goff + 146 * kcomp + j);

                    auto g_zzz_yyzz = braBuffer.data(goff + 147 * kcomp + j);

                    auto g_zzz_yzzz = braBuffer.data(goff + 148 * kcomp + j);

                    auto g_zzz_zzzz = braBuffer.data(goff + 149 * kcomp + j);

                    #pragma omp simd aligned(g2_xx_xxxx, g2_xx_xxxy, g2_xx_xxxz,\
                                             g2_xx_xxyy, g2_xx_xxyz, g2_xx_xxzz,\
                                             g2_xx_xyyy, g2_xx_xyyz, g2_xx_xyzz,\
                                             g2_xx_xzzz, g2_xx_yyyy, g2_xx_yyyz,\
                                             g2_xx_yyzz, g2_xx_yzzz, g2_xx_zzzz,\
                                             g2_xy_xxxx, g2_xy_xxxy, g2_xy_xxxz,\
                                             g2_xy_xxyy, g2_xy_xxyz, g2_xy_xxzz,\
                                             g2_xy_xyyy, g2_xy_xyyz, g2_xy_xyzz,\
                                             g2_xy_xzzz, g2_xy_yyyy, g2_xy_yyyz,\
                                             g2_xy_yyzz, g2_xy_yzzz, g2_xy_zzzz,\
                                             g2_xz_xxxx, g2_xz_xxxy, g2_xz_xxxz,\
                                             g2_xz_xxyy, g2_xz_xxyz, g2_xz_xxzz,\
                                             g2_xz_xyyy, g2_xz_xyyz, g2_xz_xyzz,\
                                             g2_xz_xzzz, g2_xz_yyyy, g2_xz_yyyz,\
                                             g2_xz_yyzz, g2_xz_yzzz, g2_xz_zzzz,\
                                             g2_yy_xxxx, g2_yy_xxxy, g2_yy_xxxz,\
                                             g2_yy_xxyy, g2_yy_xxyz, g2_yy_xxzz,\
                                             g2_yy_xyyy, g2_yy_xyyz, g2_yy_xyzz,\
                                             g2_yy_xzzz, g2_yy_yyyy, g2_yy_yyyz,\
                                             g2_yy_yyzz, g2_yy_yzzz, g2_yy_zzzz,\
                                             g2_yz_xxxx, g2_yz_xxxy, g2_yz_xxxz,\
                                             g2_yz_xxyy, g2_yz_xxyz, g2_yz_xxzz,\
                                             g2_yz_xyyy, g2_yz_xyyz, g2_yz_xyzz,\
                                             g2_yz_xzzz, g2_yz_yyyy, g2_yz_yyyz,\
                                             g2_yz_yyzz, g2_yz_yzzz, g2_yz_zzzz,\
                                             g2_zz_xxxx, g2_zz_xxxy, g2_zz_xxxz,\
                                             g2_zz_xxyy, g2_zz_xxyz, g2_zz_xxzz,\
                                             g2_zz_xyyy, g2_zz_xyyz, g2_zz_xyzz,\
                                             g2_zz_xzzz, g2_zz_yyyy, g2_zz_yyyz,\
                                             g2_zz_yyzz, g2_zz_yzzz, g2_zz_zzzz,\
                                             g1_xx_xxxxx, g1_xx_xxxxy, g1_xx_xxxxz,\
                                             g1_xx_xxxyy, g1_xx_xxxyz, g1_xx_xxxzz,\
                                             g1_xx_xxyyy, g1_xx_xxyyz, g1_xx_xxyzz,\
                                             g1_xx_xxzzz, g1_xx_xyyyy, g1_xx_xyyyz,\
                                             g1_xx_xyyzz, g1_xx_xyzzz, g1_xx_xzzzz,\
                                             g1_xy_xxxxx, g1_xy_xxxxy, g1_xy_xxxxz,\
                                             g1_xy_xxxyy, g1_xy_xxxyz, g1_xy_xxxzz,\
                                             g1_xy_xxyyy, g1_xy_xxyyz, g1_xy_xxyzz,\
                                             g1_xy_xxzzz, g1_xy_xyyyy, g1_xy_xyyyz,\
                                             g1_xy_xyyzz, g1_xy_xyzzz, g1_xy_xzzzz,\
                                             g1_xz_xxxxx, g1_xz_xxxxy, g1_xz_xxxxz,\
                                             g1_xz_xxxyy, g1_xz_xxxyz, g1_xz_xxxzz,\
                                             g1_xz_xxyyy, g1_xz_xxyyz, g1_xz_xxyzz,\
                                             g1_xz_xxzzz, g1_xz_xyyyy, g1_xz_xyyyz,\
                                             g1_xz_xyyzz, g1_xz_xyzzz, g1_xz_xzzzz,\
                                             g1_yy_xxxxx, g1_yy_xxxxy, g1_yy_xxxxz,\
                                             g1_yy_xxxyy, g1_yy_xxxyz, g1_yy_xxxzz,\
                                             g1_yy_xxyyy, g1_yy_xxyyz, g1_yy_xxyzz,\
                                             g1_yy_xxzzz, g1_yy_xyyyy, g1_yy_xyyyz,\
                                             g1_yy_xyyzz, g1_yy_xyzzz, g1_yy_xzzzz,\
                                             g1_yy_yyyyy, g1_yy_yyyyz, g1_yy_yyyzz,\
                                             g1_yy_yyzzz, g1_yy_yzzzz, \
                                             g1_yz_xxxxx, g1_yz_xxxxy, g1_yz_xxxxz,\
                                             g1_yz_xxxyy, g1_yz_xxxyz, g1_yz_xxxzz,\
                                             g1_yz_xxyyy, g1_yz_xxyyz, g1_yz_xxyzz,\
                                             g1_yz_xxzzz, g1_yz_xyyyy, g1_yz_xyyyz,\
                                             g1_yz_xyyzz, g1_yz_xyzzz, g1_yz_xzzzz,\
                                             g1_yz_yyyyy, g1_yz_yyyyz, g1_yz_yyyzz,\
                                             g1_yz_yyzzz, g1_yz_yzzzz, \
                                             g1_zz_xxxxx, g1_zz_xxxxy, g1_zz_xxxxz,\
                                             g1_zz_xxxyy, g1_zz_xxxyz, g1_zz_xxxzz,\
                                             g1_zz_xxyyy, g1_zz_xxyyz, g1_zz_xxyzz,\
                                             g1_zz_xxzzz, g1_zz_xyyyy, g1_zz_xyyyz,\
                                             g1_zz_xyyzz, g1_zz_xyzzz, g1_zz_xzzzz,\
                                             g1_zz_yyyyy, g1_zz_yyyyz, g1_zz_yyyzz,\
                                             g1_zz_yyzzz, g1_zz_yzzzz, g1_zz_zzzzz,\
                                             g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz,\
                                             g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz,\
                                             g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz,\
                                             g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz,\
                                             g_xxx_yyzz, g_xxx_yzzz, g_xxx_zzzz,\
                                             g_xxy_xxxx, g_xxy_xxxy, g_xxy_xxxz,\
                                             g_xxy_xxyy, g_xxy_xxyz, g_xxy_xxzz,\
                                             g_xxy_xyyy, g_xxy_xyyz, g_xxy_xyzz,\
                                             g_xxy_xzzz, g_xxy_yyyy, g_xxy_yyyz,\
                                             g_xxy_yyzz, g_xxy_yzzz, g_xxy_zzzz,\
                                             g_xxz_xxxx, g_xxz_xxxy, g_xxz_xxxz,\
                                             g_xxz_xxyy, g_xxz_xxyz, g_xxz_xxzz,\
                                             g_xxz_xyyy, g_xxz_xyyz, g_xxz_xyzz,\
                                             g_xxz_xzzz, g_xxz_yyyy, g_xxz_yyyz,\
                                             g_xxz_yyzz, g_xxz_yzzz, g_xxz_zzzz,\
                                             g_xyy_xxxx, g_xyy_xxxy, g_xyy_xxxz,\
                                             g_xyy_xxyy, g_xyy_xxyz, g_xyy_xxzz,\
                                             g_xyy_xyyy, g_xyy_xyyz, g_xyy_xyzz,\
                                             g_xyy_xzzz, g_xyy_yyyy, g_xyy_yyyz,\
                                             g_xyy_yyzz, g_xyy_yzzz, g_xyy_zzzz,\
                                             g_xyz_xxxx, g_xyz_xxxy, g_xyz_xxxz,\
                                             g_xyz_xxyy, g_xyz_xxyz, g_xyz_xxzz,\
                                             g_xyz_xyyy, g_xyz_xyyz, g_xyz_xyzz,\
                                             g_xyz_xzzz, g_xyz_yyyy, g_xyz_yyyz,\
                                             g_xyz_yyzz, g_xyz_yzzz, g_xyz_zzzz,\
                                             g_xzz_xxxx, g_xzz_xxxy, g_xzz_xxxz,\
                                             g_xzz_xxyy, g_xzz_xxyz, g_xzz_xxzz,\
                                             g_xzz_xyyy, g_xzz_xyyz, g_xzz_xyzz,\
                                             g_xzz_xzzz, g_xzz_yyyy, g_xzz_yyyz,\
                                             g_xzz_yyzz, g_xzz_yzzz, g_xzz_zzzz,\
                                             g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz,\
                                             g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz,\
                                             g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyzz,\
                                             g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz,\
                                             g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz,\
                                             g_yyz_xxxx, g_yyz_xxxy, g_yyz_xxxz,\
                                             g_yyz_xxyy, g_yyz_xxyz, g_yyz_xxzz,\
                                             g_yyz_xyyy, g_yyz_xyyz, g_yyz_xyzz,\
                                             g_yyz_xzzz, g_yyz_yyyy, g_yyz_yyyz,\
                                             g_yyz_yyzz, g_yyz_yzzz, g_yyz_zzzz,\
                                             g_yzz_xxxx, g_yzz_xxxy, g_yzz_xxxz,\
                                             g_yzz_xxyy, g_yzz_xxyz, g_yzz_xxzz,\
                                             g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyzz,\
                                             g_yzz_xzzz, g_yzz_yyyy, g_yzz_yyyz,\
                                             g_yzz_yyzz, g_yzz_yzzz, g_yzz_zzzz,\
                                             g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz,\
                                             g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz,\
                                             g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz,\
                                             g_zzz_xzzz, g_zzz_yyyy, g_zzz_yyyz,\
                                             g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xxx_xxxx[k] = g1_xx_xxxxx[k] - abx * g2_xx_xxxx[k];

                        g_xxx_xxxy[k] = g1_xx_xxxxy[k] - abx * g2_xx_xxxy[k];

                        g_xxx_xxxz[k] = g1_xx_xxxxz[k] - abx * g2_xx_xxxz[k];

                        g_xxx_xxyy[k] = g1_xx_xxxyy[k] - abx * g2_xx_xxyy[k];

                        g_xxx_xxyz[k] = g1_xx_xxxyz[k] - abx * g2_xx_xxyz[k];

                        g_xxx_xxzz[k] = g1_xx_xxxzz[k] - abx * g2_xx_xxzz[k];

                        g_xxx_xyyy[k] = g1_xx_xxyyy[k] - abx * g2_xx_xyyy[k];

                        g_xxx_xyyz[k] = g1_xx_xxyyz[k] - abx * g2_xx_xyyz[k];

                        g_xxx_xyzz[k] = g1_xx_xxyzz[k] - abx * g2_xx_xyzz[k];

                        g_xxx_xzzz[k] = g1_xx_xxzzz[k] - abx * g2_xx_xzzz[k];

                        g_xxx_yyyy[k] = g1_xx_xyyyy[k] - abx * g2_xx_yyyy[k];

                        g_xxx_yyyz[k] = g1_xx_xyyyz[k] - abx * g2_xx_yyyz[k];

                        g_xxx_yyzz[k] = g1_xx_xyyzz[k] - abx * g2_xx_yyzz[k];

                        g_xxx_yzzz[k] = g1_xx_xyzzz[k] - abx * g2_xx_yzzz[k];

                        g_xxx_zzzz[k] = g1_xx_xzzzz[k] - abx * g2_xx_zzzz[k];

                        g_xxy_xxxx[k] = g1_xy_xxxxx[k] - abx * g2_xy_xxxx[k];

                        g_xxy_xxxy[k] = g1_xy_xxxxy[k] - abx * g2_xy_xxxy[k];

                        g_xxy_xxxz[k] = g1_xy_xxxxz[k] - abx * g2_xy_xxxz[k];

                        g_xxy_xxyy[k] = g1_xy_xxxyy[k] - abx * g2_xy_xxyy[k];

                        g_xxy_xxyz[k] = g1_xy_xxxyz[k] - abx * g2_xy_xxyz[k];

                        g_xxy_xxzz[k] = g1_xy_xxxzz[k] - abx * g2_xy_xxzz[k];

                        g_xxy_xyyy[k] = g1_xy_xxyyy[k] - abx * g2_xy_xyyy[k];

                        g_xxy_xyyz[k] = g1_xy_xxyyz[k] - abx * g2_xy_xyyz[k];

                        g_xxy_xyzz[k] = g1_xy_xxyzz[k] - abx * g2_xy_xyzz[k];

                        g_xxy_xzzz[k] = g1_xy_xxzzz[k] - abx * g2_xy_xzzz[k];

                        g_xxy_yyyy[k] = g1_xy_xyyyy[k] - abx * g2_xy_yyyy[k];

                        g_xxy_yyyz[k] = g1_xy_xyyyz[k] - abx * g2_xy_yyyz[k];

                        g_xxy_yyzz[k] = g1_xy_xyyzz[k] - abx * g2_xy_yyzz[k];

                        g_xxy_yzzz[k] = g1_xy_xyzzz[k] - abx * g2_xy_yzzz[k];

                        g_xxy_zzzz[k] = g1_xy_xzzzz[k] - abx * g2_xy_zzzz[k];

                        g_xxz_xxxx[k] = g1_xz_xxxxx[k] - abx * g2_xz_xxxx[k];

                        g_xxz_xxxy[k] = g1_xz_xxxxy[k] - abx * g2_xz_xxxy[k];

                        g_xxz_xxxz[k] = g1_xz_xxxxz[k] - abx * g2_xz_xxxz[k];

                        g_xxz_xxyy[k] = g1_xz_xxxyy[k] - abx * g2_xz_xxyy[k];

                        g_xxz_xxyz[k] = g1_xz_xxxyz[k] - abx * g2_xz_xxyz[k];

                        g_xxz_xxzz[k] = g1_xz_xxxzz[k] - abx * g2_xz_xxzz[k];

                        g_xxz_xyyy[k] = g1_xz_xxyyy[k] - abx * g2_xz_xyyy[k];

                        g_xxz_xyyz[k] = g1_xz_xxyyz[k] - abx * g2_xz_xyyz[k];

                        g_xxz_xyzz[k] = g1_xz_xxyzz[k] - abx * g2_xz_xyzz[k];

                        g_xxz_xzzz[k] = g1_xz_xxzzz[k] - abx * g2_xz_xzzz[k];

                        g_xxz_yyyy[k] = g1_xz_xyyyy[k] - abx * g2_xz_yyyy[k];

                        g_xxz_yyyz[k] = g1_xz_xyyyz[k] - abx * g2_xz_yyyz[k];

                        g_xxz_yyzz[k] = g1_xz_xyyzz[k] - abx * g2_xz_yyzz[k];

                        g_xxz_yzzz[k] = g1_xz_xyzzz[k] - abx * g2_xz_yzzz[k];

                        g_xxz_zzzz[k] = g1_xz_xzzzz[k] - abx * g2_xz_zzzz[k];

                        g_xyy_xxxx[k] = g1_yy_xxxxx[k] - abx * g2_yy_xxxx[k];

                        g_xyy_xxxy[k] = g1_yy_xxxxy[k] - abx * g2_yy_xxxy[k];

                        g_xyy_xxxz[k] = g1_yy_xxxxz[k] - abx * g2_yy_xxxz[k];

                        g_xyy_xxyy[k] = g1_yy_xxxyy[k] - abx * g2_yy_xxyy[k];

                        g_xyy_xxyz[k] = g1_yy_xxxyz[k] - abx * g2_yy_xxyz[k];

                        g_xyy_xxzz[k] = g1_yy_xxxzz[k] - abx * g2_yy_xxzz[k];

                        g_xyy_xyyy[k] = g1_yy_xxyyy[k] - abx * g2_yy_xyyy[k];

                        g_xyy_xyyz[k] = g1_yy_xxyyz[k] - abx * g2_yy_xyyz[k];

                        g_xyy_xyzz[k] = g1_yy_xxyzz[k] - abx * g2_yy_xyzz[k];

                        g_xyy_xzzz[k] = g1_yy_xxzzz[k] - abx * g2_yy_xzzz[k];

                        g_xyy_yyyy[k] = g1_yy_xyyyy[k] - abx * g2_yy_yyyy[k];

                        g_xyy_yyyz[k] = g1_yy_xyyyz[k] - abx * g2_yy_yyyz[k];

                        g_xyy_yyzz[k] = g1_yy_xyyzz[k] - abx * g2_yy_yyzz[k];

                        g_xyy_yzzz[k] = g1_yy_xyzzz[k] - abx * g2_yy_yzzz[k];

                        g_xyy_zzzz[k] = g1_yy_xzzzz[k] - abx * g2_yy_zzzz[k];

                        g_xyz_xxxx[k] = g1_yz_xxxxx[k] - abx * g2_yz_xxxx[k];

                        g_xyz_xxxy[k] = g1_yz_xxxxy[k] - abx * g2_yz_xxxy[k];

                        g_xyz_xxxz[k] = g1_yz_xxxxz[k] - abx * g2_yz_xxxz[k];

                        g_xyz_xxyy[k] = g1_yz_xxxyy[k] - abx * g2_yz_xxyy[k];

                        g_xyz_xxyz[k] = g1_yz_xxxyz[k] - abx * g2_yz_xxyz[k];

                        g_xyz_xxzz[k] = g1_yz_xxxzz[k] - abx * g2_yz_xxzz[k];

                        g_xyz_xyyy[k] = g1_yz_xxyyy[k] - abx * g2_yz_xyyy[k];

                        g_xyz_xyyz[k] = g1_yz_xxyyz[k] - abx * g2_yz_xyyz[k];

                        g_xyz_xyzz[k] = g1_yz_xxyzz[k] - abx * g2_yz_xyzz[k];

                        g_xyz_xzzz[k] = g1_yz_xxzzz[k] - abx * g2_yz_xzzz[k];

                        g_xyz_yyyy[k] = g1_yz_xyyyy[k] - abx * g2_yz_yyyy[k];

                        g_xyz_yyyz[k] = g1_yz_xyyyz[k] - abx * g2_yz_yyyz[k];

                        g_xyz_yyzz[k] = g1_yz_xyyzz[k] - abx * g2_yz_yyzz[k];

                        g_xyz_yzzz[k] = g1_yz_xyzzz[k] - abx * g2_yz_yzzz[k];

                        g_xyz_zzzz[k] = g1_yz_xzzzz[k] - abx * g2_yz_zzzz[k];

                        g_xzz_xxxx[k] = g1_zz_xxxxx[k] - abx * g2_zz_xxxx[k];

                        g_xzz_xxxy[k] = g1_zz_xxxxy[k] - abx * g2_zz_xxxy[k];

                        g_xzz_xxxz[k] = g1_zz_xxxxz[k] - abx * g2_zz_xxxz[k];

                        g_xzz_xxyy[k] = g1_zz_xxxyy[k] - abx * g2_zz_xxyy[k];

                        g_xzz_xxyz[k] = g1_zz_xxxyz[k] - abx * g2_zz_xxyz[k];

                        g_xzz_xxzz[k] = g1_zz_xxxzz[k] - abx * g2_zz_xxzz[k];

                        g_xzz_xyyy[k] = g1_zz_xxyyy[k] - abx * g2_zz_xyyy[k];

                        g_xzz_xyyz[k] = g1_zz_xxyyz[k] - abx * g2_zz_xyyz[k];

                        g_xzz_xyzz[k] = g1_zz_xxyzz[k] - abx * g2_zz_xyzz[k];

                        g_xzz_xzzz[k] = g1_zz_xxzzz[k] - abx * g2_zz_xzzz[k];

                        g_xzz_yyyy[k] = g1_zz_xyyyy[k] - abx * g2_zz_yyyy[k];

                        g_xzz_yyyz[k] = g1_zz_xyyyz[k] - abx * g2_zz_yyyz[k];

                        g_xzz_yyzz[k] = g1_zz_xyyzz[k] - abx * g2_zz_yyzz[k];

                        g_xzz_yzzz[k] = g1_zz_xyzzz[k] - abx * g2_zz_yzzz[k];

                        g_xzz_zzzz[k] = g1_zz_xzzzz[k] - abx * g2_zz_zzzz[k];

                        // leading y component

                        g_yyy_xxxx[k] = g1_yy_xxxxy[k] - aby * g2_yy_xxxx[k];

                        g_yyy_xxxy[k] = g1_yy_xxxyy[k] - aby * g2_yy_xxxy[k];

                        g_yyy_xxxz[k] = g1_yy_xxxyz[k] - aby * g2_yy_xxxz[k];

                        g_yyy_xxyy[k] = g1_yy_xxyyy[k] - aby * g2_yy_xxyy[k];

                        g_yyy_xxyz[k] = g1_yy_xxyyz[k] - aby * g2_yy_xxyz[k];

                        g_yyy_xxzz[k] = g1_yy_xxyzz[k] - aby * g2_yy_xxzz[k];

                        g_yyy_xyyy[k] = g1_yy_xyyyy[k] - aby * g2_yy_xyyy[k];

                        g_yyy_xyyz[k] = g1_yy_xyyyz[k] - aby * g2_yy_xyyz[k];

                        g_yyy_xyzz[k] = g1_yy_xyyzz[k] - aby * g2_yy_xyzz[k];

                        g_yyy_xzzz[k] = g1_yy_xyzzz[k] - aby * g2_yy_xzzz[k];

                        g_yyy_yyyy[k] = g1_yy_yyyyy[k] - aby * g2_yy_yyyy[k];

                        g_yyy_yyyz[k] = g1_yy_yyyyz[k] - aby * g2_yy_yyyz[k];

                        g_yyy_yyzz[k] = g1_yy_yyyzz[k] - aby * g2_yy_yyzz[k];

                        g_yyy_yzzz[k] = g1_yy_yyzzz[k] - aby * g2_yy_yzzz[k];

                        g_yyy_zzzz[k] = g1_yy_yzzzz[k] - aby * g2_yy_zzzz[k];

                        g_yyz_xxxx[k] = g1_yz_xxxxy[k] - aby * g2_yz_xxxx[k];

                        g_yyz_xxxy[k] = g1_yz_xxxyy[k] - aby * g2_yz_xxxy[k];

                        g_yyz_xxxz[k] = g1_yz_xxxyz[k] - aby * g2_yz_xxxz[k];

                        g_yyz_xxyy[k] = g1_yz_xxyyy[k] - aby * g2_yz_xxyy[k];

                        g_yyz_xxyz[k] = g1_yz_xxyyz[k] - aby * g2_yz_xxyz[k];

                        g_yyz_xxzz[k] = g1_yz_xxyzz[k] - aby * g2_yz_xxzz[k];

                        g_yyz_xyyy[k] = g1_yz_xyyyy[k] - aby * g2_yz_xyyy[k];

                        g_yyz_xyyz[k] = g1_yz_xyyyz[k] - aby * g2_yz_xyyz[k];

                        g_yyz_xyzz[k] = g1_yz_xyyzz[k] - aby * g2_yz_xyzz[k];

                        g_yyz_xzzz[k] = g1_yz_xyzzz[k] - aby * g2_yz_xzzz[k];

                        g_yyz_yyyy[k] = g1_yz_yyyyy[k] - aby * g2_yz_yyyy[k];

                        g_yyz_yyyz[k] = g1_yz_yyyyz[k] - aby * g2_yz_yyyz[k];

                        g_yyz_yyzz[k] = g1_yz_yyyzz[k] - aby * g2_yz_yyzz[k];

                        g_yyz_yzzz[k] = g1_yz_yyzzz[k] - aby * g2_yz_yzzz[k];

                        g_yyz_zzzz[k] = g1_yz_yzzzz[k] - aby * g2_yz_zzzz[k];

                        g_yzz_xxxx[k] = g1_zz_xxxxy[k] - aby * g2_zz_xxxx[k];

                        g_yzz_xxxy[k] = g1_zz_xxxyy[k] - aby * g2_zz_xxxy[k];

                        g_yzz_xxxz[k] = g1_zz_xxxyz[k] - aby * g2_zz_xxxz[k];

                        g_yzz_xxyy[k] = g1_zz_xxyyy[k] - aby * g2_zz_xxyy[k];

                        g_yzz_xxyz[k] = g1_zz_xxyyz[k] - aby * g2_zz_xxyz[k];

                        g_yzz_xxzz[k] = g1_zz_xxyzz[k] - aby * g2_zz_xxzz[k];

                        g_yzz_xyyy[k] = g1_zz_xyyyy[k] - aby * g2_zz_xyyy[k];

                        g_yzz_xyyz[k] = g1_zz_xyyyz[k] - aby * g2_zz_xyyz[k];

                        g_yzz_xyzz[k] = g1_zz_xyyzz[k] - aby * g2_zz_xyzz[k];

                        g_yzz_xzzz[k] = g1_zz_xyzzz[k] - aby * g2_zz_xzzz[k];

                        g_yzz_yyyy[k] = g1_zz_yyyyy[k] - aby * g2_zz_yyyy[k];

                        g_yzz_yyyz[k] = g1_zz_yyyyz[k] - aby * g2_zz_yyyz[k];

                        g_yzz_yyzz[k] = g1_zz_yyyzz[k] - aby * g2_zz_yyzz[k];

                        g_yzz_yzzz[k] = g1_zz_yyzzz[k] - aby * g2_zz_yzzz[k];

                        g_yzz_zzzz[k] = g1_zz_yzzzz[k] - aby * g2_zz_zzzz[k];

                        // leading z component

                        g_zzz_xxxx[k] = g1_zz_xxxxz[k] - abz * g2_zz_xxxx[k];

                        g_zzz_xxxy[k] = g1_zz_xxxyz[k] - abz * g2_zz_xxxy[k];

                        g_zzz_xxxz[k] = g1_zz_xxxzz[k] - abz * g2_zz_xxxz[k];

                        g_zzz_xxyy[k] = g1_zz_xxyyz[k] - abz * g2_zz_xxyy[k];

                        g_zzz_xxyz[k] = g1_zz_xxyzz[k] - abz * g2_zz_xxyz[k];

                        g_zzz_xxzz[k] = g1_zz_xxzzz[k] - abz * g2_zz_xxzz[k];

                        g_zzz_xyyy[k] = g1_zz_xyyyz[k] - abz * g2_zz_xyyy[k];

                        g_zzz_xyyz[k] = g1_zz_xyyzz[k] - abz * g2_zz_xyyz[k];

                        g_zzz_xyzz[k] = g1_zz_xyzzz[k] - abz * g2_zz_xyzz[k];

                        g_zzz_xzzz[k] = g1_zz_xzzzz[k] - abz * g2_zz_xzzz[k];

                        g_zzz_yyyy[k] = g1_zz_yyyyz[k] - abz * g2_zz_yyyy[k];

                        g_zzz_yyyz[k] = g1_zz_yyyzz[k] - abz * g2_zz_yyyz[k];

                        g_zzz_yyzz[k] = g1_zz_yyzzz[k] - abz * g2_zz_yyzz[k];

                        g_zzz_yzzz[k] = g1_zz_yzzzz[k] - abz * g2_zz_yzzz[k];

                        g_zzz_zzzz[k] = g1_zz_zzzzz[k] - abz * g2_zz_zzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForFHXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 3) && (recPattern[i].second() == 5))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (35|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {3, 5, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 6, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {2, 5, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (DH|g(r,r')|XX) integrals

                    auto g2_xx_xxxxx = braBuffer.data(g2off + j);

                    auto g2_xx_xxxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_xx_xxxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_xx_xxxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_xx_xxxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_xx_xxxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_xx_xxyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_xx_xxyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_xx_xxyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_xx_xxzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_xx_xyyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_xx_xyyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_xx_xyyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_xx_xyzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_xx_xzzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_xx_yyyyy = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_xx_yyyyz = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_xx_yyyzz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_xx_yyzzz = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_xx_yzzzz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_xx_zzzzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_xy_xxxxx = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_xy_xxxxy = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_xy_xxxxz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_xy_xxxyy = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_xy_xxxyz = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_xy_xxxzz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_xy_xxyyy = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_xy_xxyyz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_xy_xxyzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_xy_xxzzz = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_xy_xyyyy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_xy_xyyyz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_xy_xyyzz = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_xy_xyzzz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_xy_xzzzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_xy_yyyyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_xy_yyyyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_xy_yyyzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_xy_yyzzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_xy_yzzzz = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_xy_zzzzz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_xz_xxxxx = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_xz_xxxxy = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_xz_xxxxz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_xz_xxxyy = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_xz_xxxyz = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_xz_xxxzz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_xz_xxyyy = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_xz_xxyyz = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_xz_xxyzz = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_xz_xxzzz = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_xz_xyyyy = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_xz_xyyyz = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_xz_xyyzz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_xz_xyzzz = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_xz_xzzzz = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_xz_yyyyy = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_xz_yyyyz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_xz_yyyzz = braBuffer.data(g2off + 59 * kcomp + j);

                    auto g2_xz_yyzzz = braBuffer.data(g2off + 60 * kcomp + j);

                    auto g2_xz_yzzzz = braBuffer.data(g2off + 61 * kcomp + j);

                    auto g2_xz_zzzzz = braBuffer.data(g2off + 62 * kcomp + j);

                    auto g2_yy_xxxxx = braBuffer.data(g2off + 63 * kcomp + j);

                    auto g2_yy_xxxxy = braBuffer.data(g2off + 64 * kcomp + j);

                    auto g2_yy_xxxxz = braBuffer.data(g2off + 65 * kcomp + j);

                    auto g2_yy_xxxyy = braBuffer.data(g2off + 66 * kcomp + j);

                    auto g2_yy_xxxyz = braBuffer.data(g2off + 67 * kcomp + j);

                    auto g2_yy_xxxzz = braBuffer.data(g2off + 68 * kcomp + j);

                    auto g2_yy_xxyyy = braBuffer.data(g2off + 69 * kcomp + j);

                    auto g2_yy_xxyyz = braBuffer.data(g2off + 70 * kcomp + j);

                    auto g2_yy_xxyzz = braBuffer.data(g2off + 71 * kcomp + j);

                    auto g2_yy_xxzzz = braBuffer.data(g2off + 72 * kcomp + j);

                    auto g2_yy_xyyyy = braBuffer.data(g2off + 73 * kcomp + j);

                    auto g2_yy_xyyyz = braBuffer.data(g2off + 74 * kcomp + j);

                    auto g2_yy_xyyzz = braBuffer.data(g2off + 75 * kcomp + j);

                    auto g2_yy_xyzzz = braBuffer.data(g2off + 76 * kcomp + j);

                    auto g2_yy_xzzzz = braBuffer.data(g2off + 77 * kcomp + j);

                    auto g2_yy_yyyyy = braBuffer.data(g2off + 78 * kcomp + j);

                    auto g2_yy_yyyyz = braBuffer.data(g2off + 79 * kcomp + j);

                    auto g2_yy_yyyzz = braBuffer.data(g2off + 80 * kcomp + j);

                    auto g2_yy_yyzzz = braBuffer.data(g2off + 81 * kcomp + j);

                    auto g2_yy_yzzzz = braBuffer.data(g2off + 82 * kcomp + j);

                    auto g2_yy_zzzzz = braBuffer.data(g2off + 83 * kcomp + j);

                    auto g2_yz_xxxxx = braBuffer.data(g2off + 84 * kcomp + j);

                    auto g2_yz_xxxxy = braBuffer.data(g2off + 85 * kcomp + j);

                    auto g2_yz_xxxxz = braBuffer.data(g2off + 86 * kcomp + j);

                    auto g2_yz_xxxyy = braBuffer.data(g2off + 87 * kcomp + j);

                    auto g2_yz_xxxyz = braBuffer.data(g2off + 88 * kcomp + j);

                    auto g2_yz_xxxzz = braBuffer.data(g2off + 89 * kcomp + j);

                    auto g2_yz_xxyyy = braBuffer.data(g2off + 90 * kcomp + j);

                    auto g2_yz_xxyyz = braBuffer.data(g2off + 91 * kcomp + j);

                    auto g2_yz_xxyzz = braBuffer.data(g2off + 92 * kcomp + j);

                    auto g2_yz_xxzzz = braBuffer.data(g2off + 93 * kcomp + j);

                    auto g2_yz_xyyyy = braBuffer.data(g2off + 94 * kcomp + j);

                    auto g2_yz_xyyyz = braBuffer.data(g2off + 95 * kcomp + j);

                    auto g2_yz_xyyzz = braBuffer.data(g2off + 96 * kcomp + j);

                    auto g2_yz_xyzzz = braBuffer.data(g2off + 97 * kcomp + j);

                    auto g2_yz_xzzzz = braBuffer.data(g2off + 98 * kcomp + j);

                    auto g2_yz_yyyyy = braBuffer.data(g2off + 99 * kcomp + j);

                    auto g2_yz_yyyyz = braBuffer.data(g2off + 100 * kcomp + j);

                    auto g2_yz_yyyzz = braBuffer.data(g2off + 101 * kcomp + j);

                    auto g2_yz_yyzzz = braBuffer.data(g2off + 102 * kcomp + j);

                    auto g2_yz_yzzzz = braBuffer.data(g2off + 103 * kcomp + j);

                    auto g2_yz_zzzzz = braBuffer.data(g2off + 104 * kcomp + j);

                    auto g2_zz_xxxxx = braBuffer.data(g2off + 105 * kcomp + j);

                    auto g2_zz_xxxxy = braBuffer.data(g2off + 106 * kcomp + j);

                    auto g2_zz_xxxxz = braBuffer.data(g2off + 107 * kcomp + j);

                    auto g2_zz_xxxyy = braBuffer.data(g2off + 108 * kcomp + j);

                    auto g2_zz_xxxyz = braBuffer.data(g2off + 109 * kcomp + j);

                    auto g2_zz_xxxzz = braBuffer.data(g2off + 110 * kcomp + j);

                    auto g2_zz_xxyyy = braBuffer.data(g2off + 111 * kcomp + j);

                    auto g2_zz_xxyyz = braBuffer.data(g2off + 112 * kcomp + j);

                    auto g2_zz_xxyzz = braBuffer.data(g2off + 113 * kcomp + j);

                    auto g2_zz_xxzzz = braBuffer.data(g2off + 114 * kcomp + j);

                    auto g2_zz_xyyyy = braBuffer.data(g2off + 115 * kcomp + j);

                    auto g2_zz_xyyyz = braBuffer.data(g2off + 116 * kcomp + j);

                    auto g2_zz_xyyzz = braBuffer.data(g2off + 117 * kcomp + j);

                    auto g2_zz_xyzzz = braBuffer.data(g2off + 118 * kcomp + j);

                    auto g2_zz_xzzzz = braBuffer.data(g2off + 119 * kcomp + j);

                    auto g2_zz_yyyyy = braBuffer.data(g2off + 120 * kcomp + j);

                    auto g2_zz_yyyyz = braBuffer.data(g2off + 121 * kcomp + j);

                    auto g2_zz_yyyzz = braBuffer.data(g2off + 122 * kcomp + j);

                    auto g2_zz_yyzzz = braBuffer.data(g2off + 123 * kcomp + j);

                    auto g2_zz_yzzzz = braBuffer.data(g2off + 124 * kcomp + j);

                    auto g2_zz_zzzzz = braBuffer.data(g2off + 125 * kcomp + j);

                    // set up pointers to (DI|g(r,r')|XX) integrals

                    auto g1_xx_xxxxxx = braBuffer.data(g1off + j);

                    auto g1_xx_xxxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_xx_xxxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_xx_xxxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_xx_xxxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_xx_xxxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_xx_xxxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_xx_xxxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_xx_xxxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_xx_xxxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_xx_xxyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_xx_xxyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_xx_xxyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_xx_xxyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_xx_xxzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_xx_xyyyyy = braBuffer.data(g1off + 15 * kcomp + j);

                    auto g1_xx_xyyyyz = braBuffer.data(g1off + 16 * kcomp + j);

                    auto g1_xx_xyyyzz = braBuffer.data(g1off + 17 * kcomp + j);

                    auto g1_xx_xyyzzz = braBuffer.data(g1off + 18 * kcomp + j);

                    auto g1_xx_xyzzzz = braBuffer.data(g1off + 19 * kcomp + j);

                    auto g1_xx_xzzzzz = braBuffer.data(g1off + 20 * kcomp + j);

                    auto g1_xy_xxxxxx = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_xy_xxxxxy = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_xy_xxxxxz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_xy_xxxxyy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_xy_xxxxyz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_xy_xxxxzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_xy_xxxyyy = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_xy_xxxyyz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_xy_xxxyzz = braBuffer.data(g1off + 36 * kcomp + j);

                    auto g1_xy_xxxzzz = braBuffer.data(g1off + 37 * kcomp + j);

                    auto g1_xy_xxyyyy = braBuffer.data(g1off + 38 * kcomp + j);

                    auto g1_xy_xxyyyz = braBuffer.data(g1off + 39 * kcomp + j);

                    auto g1_xy_xxyyzz = braBuffer.data(g1off + 40 * kcomp + j);

                    auto g1_xy_xxyzzz = braBuffer.data(g1off + 41 * kcomp + j);

                    auto g1_xy_xxzzzz = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_xy_xyyyyy = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_xy_xyyyyz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_xy_xyyyzz = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_xy_xyyzzz = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_xy_xyzzzz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_xy_xzzzzz = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_xz_xxxxxx = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_xz_xxxxxy = braBuffer.data(g1off + 57 * kcomp + j);

                    auto g1_xz_xxxxxz = braBuffer.data(g1off + 58 * kcomp + j);

                    auto g1_xz_xxxxyy = braBuffer.data(g1off + 59 * kcomp + j);

                    auto g1_xz_xxxxyz = braBuffer.data(g1off + 60 * kcomp + j);

                    auto g1_xz_xxxxzz = braBuffer.data(g1off + 61 * kcomp + j);

                    auto g1_xz_xxxyyy = braBuffer.data(g1off + 62 * kcomp + j);

                    auto g1_xz_xxxyyz = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_xz_xxxyzz = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_xz_xxxzzz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_xz_xxyyyy = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_xz_xxyyyz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_xz_xxyyzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_xz_xxyzzz = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_xz_xxzzzz = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_xz_xyyyyy = braBuffer.data(g1off + 71 * kcomp + j);

                    auto g1_xz_xyyyyz = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_xz_xyyyzz = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_xz_xyyzzz = braBuffer.data(g1off + 74 * kcomp + j);

                    auto g1_xz_xyzzzz = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_xz_xzzzzz = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_yy_xxxxxx = braBuffer.data(g1off + 84 * kcomp + j);

                    auto g1_yy_xxxxxy = braBuffer.data(g1off + 85 * kcomp + j);

                    auto g1_yy_xxxxxz = braBuffer.data(g1off + 86 * kcomp + j);

                    auto g1_yy_xxxxyy = braBuffer.data(g1off + 87 * kcomp + j);

                    auto g1_yy_xxxxyz = braBuffer.data(g1off + 88 * kcomp + j);

                    auto g1_yy_xxxxzz = braBuffer.data(g1off + 89 * kcomp + j);

                    auto g1_yy_xxxyyy = braBuffer.data(g1off + 90 * kcomp + j);

                    auto g1_yy_xxxyyz = braBuffer.data(g1off + 91 * kcomp + j);

                    auto g1_yy_xxxyzz = braBuffer.data(g1off + 92 * kcomp + j);

                    auto g1_yy_xxxzzz = braBuffer.data(g1off + 93 * kcomp + j);

                    auto g1_yy_xxyyyy = braBuffer.data(g1off + 94 * kcomp + j);

                    auto g1_yy_xxyyyz = braBuffer.data(g1off + 95 * kcomp + j);

                    auto g1_yy_xxyyzz = braBuffer.data(g1off + 96 * kcomp + j);

                    auto g1_yy_xxyzzz = braBuffer.data(g1off + 97 * kcomp + j);

                    auto g1_yy_xxzzzz = braBuffer.data(g1off + 98 * kcomp + j);

                    auto g1_yy_xyyyyy = braBuffer.data(g1off + 99 * kcomp + j);

                    auto g1_yy_xyyyyz = braBuffer.data(g1off + 100 * kcomp + j);

                    auto g1_yy_xyyyzz = braBuffer.data(g1off + 101 * kcomp + j);

                    auto g1_yy_xyyzzz = braBuffer.data(g1off + 102 * kcomp + j);

                    auto g1_yy_xyzzzz = braBuffer.data(g1off + 103 * kcomp + j);

                    auto g1_yy_xzzzzz = braBuffer.data(g1off + 104 * kcomp + j);

                    auto g1_yy_yyyyyy = braBuffer.data(g1off + 105 * kcomp + j);

                    auto g1_yy_yyyyyz = braBuffer.data(g1off + 106 * kcomp + j);

                    auto g1_yy_yyyyzz = braBuffer.data(g1off + 107 * kcomp + j);

                    auto g1_yy_yyyzzz = braBuffer.data(g1off + 108 * kcomp + j);

                    auto g1_yy_yyzzzz = braBuffer.data(g1off + 109 * kcomp + j);

                    auto g1_yy_yzzzzz = braBuffer.data(g1off + 110 * kcomp + j);

                    auto g1_yz_xxxxxx = braBuffer.data(g1off + 112 * kcomp + j);

                    auto g1_yz_xxxxxy = braBuffer.data(g1off + 113 * kcomp + j);

                    auto g1_yz_xxxxxz = braBuffer.data(g1off + 114 * kcomp + j);

                    auto g1_yz_xxxxyy = braBuffer.data(g1off + 115 * kcomp + j);

                    auto g1_yz_xxxxyz = braBuffer.data(g1off + 116 * kcomp + j);

                    auto g1_yz_xxxxzz = braBuffer.data(g1off + 117 * kcomp + j);

                    auto g1_yz_xxxyyy = braBuffer.data(g1off + 118 * kcomp + j);

                    auto g1_yz_xxxyyz = braBuffer.data(g1off + 119 * kcomp + j);

                    auto g1_yz_xxxyzz = braBuffer.data(g1off + 120 * kcomp + j);

                    auto g1_yz_xxxzzz = braBuffer.data(g1off + 121 * kcomp + j);

                    auto g1_yz_xxyyyy = braBuffer.data(g1off + 122 * kcomp + j);

                    auto g1_yz_xxyyyz = braBuffer.data(g1off + 123 * kcomp + j);

                    auto g1_yz_xxyyzz = braBuffer.data(g1off + 124 * kcomp + j);

                    auto g1_yz_xxyzzz = braBuffer.data(g1off + 125 * kcomp + j);

                    auto g1_yz_xxzzzz = braBuffer.data(g1off + 126 * kcomp + j);

                    auto g1_yz_xyyyyy = braBuffer.data(g1off + 127 * kcomp + j);

                    auto g1_yz_xyyyyz = braBuffer.data(g1off + 128 * kcomp + j);

                    auto g1_yz_xyyyzz = braBuffer.data(g1off + 129 * kcomp + j);

                    auto g1_yz_xyyzzz = braBuffer.data(g1off + 130 * kcomp + j);

                    auto g1_yz_xyzzzz = braBuffer.data(g1off + 131 * kcomp + j);

                    auto g1_yz_xzzzzz = braBuffer.data(g1off + 132 * kcomp + j);

                    auto g1_yz_yyyyyy = braBuffer.data(g1off + 133 * kcomp + j);

                    auto g1_yz_yyyyyz = braBuffer.data(g1off + 134 * kcomp + j);

                    auto g1_yz_yyyyzz = braBuffer.data(g1off + 135 * kcomp + j);

                    auto g1_yz_yyyzzz = braBuffer.data(g1off + 136 * kcomp + j);

                    auto g1_yz_yyzzzz = braBuffer.data(g1off + 137 * kcomp + j);

                    auto g1_yz_yzzzzz = braBuffer.data(g1off + 138 * kcomp + j);

                    auto g1_zz_xxxxxx = braBuffer.data(g1off + 140 * kcomp + j);

                    auto g1_zz_xxxxxy = braBuffer.data(g1off + 141 * kcomp + j);

                    auto g1_zz_xxxxxz = braBuffer.data(g1off + 142 * kcomp + j);

                    auto g1_zz_xxxxyy = braBuffer.data(g1off + 143 * kcomp + j);

                    auto g1_zz_xxxxyz = braBuffer.data(g1off + 144 * kcomp + j);

                    auto g1_zz_xxxxzz = braBuffer.data(g1off + 145 * kcomp + j);

                    auto g1_zz_xxxyyy = braBuffer.data(g1off + 146 * kcomp + j);

                    auto g1_zz_xxxyyz = braBuffer.data(g1off + 147 * kcomp + j);

                    auto g1_zz_xxxyzz = braBuffer.data(g1off + 148 * kcomp + j);

                    auto g1_zz_xxxzzz = braBuffer.data(g1off + 149 * kcomp + j);

                    auto g1_zz_xxyyyy = braBuffer.data(g1off + 150 * kcomp + j);

                    auto g1_zz_xxyyyz = braBuffer.data(g1off + 151 * kcomp + j);

                    auto g1_zz_xxyyzz = braBuffer.data(g1off + 152 * kcomp + j);

                    auto g1_zz_xxyzzz = braBuffer.data(g1off + 153 * kcomp + j);

                    auto g1_zz_xxzzzz = braBuffer.data(g1off + 154 * kcomp + j);

                    auto g1_zz_xyyyyy = braBuffer.data(g1off + 155 * kcomp + j);

                    auto g1_zz_xyyyyz = braBuffer.data(g1off + 156 * kcomp + j);

                    auto g1_zz_xyyyzz = braBuffer.data(g1off + 157 * kcomp + j);

                    auto g1_zz_xyyzzz = braBuffer.data(g1off + 158 * kcomp + j);

                    auto g1_zz_xyzzzz = braBuffer.data(g1off + 159 * kcomp + j);

                    auto g1_zz_xzzzzz = braBuffer.data(g1off + 160 * kcomp + j);

                    auto g1_zz_yyyyyy = braBuffer.data(g1off + 161 * kcomp + j);

                    auto g1_zz_yyyyyz = braBuffer.data(g1off + 162 * kcomp + j);

                    auto g1_zz_yyyyzz = braBuffer.data(g1off + 163 * kcomp + j);

                    auto g1_zz_yyyzzz = braBuffer.data(g1off + 164 * kcomp + j);

                    auto g1_zz_yyzzzz = braBuffer.data(g1off + 165 * kcomp + j);

                    auto g1_zz_yzzzzz = braBuffer.data(g1off + 166 * kcomp + j);

                    auto g1_zz_zzzzzz = braBuffer.data(g1off + 167 * kcomp + j);

                    // set up pointers to (FH|g(r,r')|XX) integrals

                    auto g_xxx_xxxxx = braBuffer.data(goff + j);

                    auto g_xxx_xxxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xxx_xxxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xxx_xxxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xxx_xxxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xxx_xxxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xxx_xxyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xxx_xxyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xxx_xxyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xxx_xxzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xxx_xyyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xxx_xyyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xxx_xyyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xxx_xyzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xxx_xzzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xxx_yyyyy = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xxx_yyyyz = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xxx_yyyzz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xxx_yyzzz = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xxx_yzzzz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xxx_zzzzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xxy_xxxxx = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xxy_xxxxy = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xxy_xxxxz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xxy_xxxyy = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xxy_xxxyz = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xxy_xxxzz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xxy_xxyyy = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xxy_xxyyz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xxy_xxyzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xxy_xxzzz = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xxy_xyyyy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xxy_xyyyz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xxy_xyyzz = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xxy_xyzzz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xxy_xzzzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xxy_yyyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xxy_yyyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xxy_yyyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xxy_yyzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xxy_yzzzz = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xxy_zzzzz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xxz_xxxxx = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xxz_xxxxy = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xxz_xxxxz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xxz_xxxyy = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xxz_xxxyz = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xxz_xxxzz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xxz_xxyyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xxz_xxyyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xxz_xxyzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xxz_xxzzz = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xxz_xyyyy = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xxz_xyyyz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xxz_xyyzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xxz_xyzzz = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xxz_xzzzz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xxz_yyyyy = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xxz_yyyyz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xxz_yyyzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_xxz_yyzzz = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_xxz_yzzzz = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_xxz_zzzzz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_xyy_xxxxx = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_xyy_xxxxy = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_xyy_xxxxz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_xyy_xxxyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_xyy_xxxyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_xyy_xxxzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_xyy_xxyyy = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_xyy_xxyyz = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_xyy_xxyzz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_xyy_xxzzz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_xyy_xyyyy = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_xyy_xyyyz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_xyy_xyyzz = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_xyy_xyzzz = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_xyy_xzzzz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_xyy_yyyyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_xyy_yyyyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_xyy_yyyzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_xyy_yyzzz = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_xyy_yzzzz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_xyy_zzzzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_xyz_xxxxx = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_xyz_xxxxy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_xyz_xxxxz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_xyz_xxxyy = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_xyz_xxxyz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_xyz_xxxzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_xyz_xxyyy = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_xyz_xxyyz = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_xyz_xxyzz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_xyz_xxzzz = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_xyz_xyyyy = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_xyz_xyyyz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_xyz_xyyzz = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_xyz_xyzzz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_xyz_xzzzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_xyz_yyyyy = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_xyz_yyyyz = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_xyz_yyyzz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_xyz_yyzzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_xyz_yzzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_xyz_zzzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_xzz_xxxxx = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_xzz_xxxxy = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_xzz_xxxxz = braBuffer.data(goff + 107 * kcomp + j);

                    auto g_xzz_xxxyy = braBuffer.data(goff + 108 * kcomp + j);

                    auto g_xzz_xxxyz = braBuffer.data(goff + 109 * kcomp + j);

                    auto g_xzz_xxxzz = braBuffer.data(goff + 110 * kcomp + j);

                    auto g_xzz_xxyyy = braBuffer.data(goff + 111 * kcomp + j);

                    auto g_xzz_xxyyz = braBuffer.data(goff + 112 * kcomp + j);

                    auto g_xzz_xxyzz = braBuffer.data(goff + 113 * kcomp + j);

                    auto g_xzz_xxzzz = braBuffer.data(goff + 114 * kcomp + j);

                    auto g_xzz_xyyyy = braBuffer.data(goff + 115 * kcomp + j);

                    auto g_xzz_xyyyz = braBuffer.data(goff + 116 * kcomp + j);

                    auto g_xzz_xyyzz = braBuffer.data(goff + 117 * kcomp + j);

                    auto g_xzz_xyzzz = braBuffer.data(goff + 118 * kcomp + j);

                    auto g_xzz_xzzzz = braBuffer.data(goff + 119 * kcomp + j);

                    auto g_xzz_yyyyy = braBuffer.data(goff + 120 * kcomp + j);

                    auto g_xzz_yyyyz = braBuffer.data(goff + 121 * kcomp + j);

                    auto g_xzz_yyyzz = braBuffer.data(goff + 122 * kcomp + j);

                    auto g_xzz_yyzzz = braBuffer.data(goff + 123 * kcomp + j);

                    auto g_xzz_yzzzz = braBuffer.data(goff + 124 * kcomp + j);

                    auto g_xzz_zzzzz = braBuffer.data(goff + 125 * kcomp + j);

                    auto g_yyy_xxxxx = braBuffer.data(goff + 126 * kcomp + j);

                    auto g_yyy_xxxxy = braBuffer.data(goff + 127 * kcomp + j);

                    auto g_yyy_xxxxz = braBuffer.data(goff + 128 * kcomp + j);

                    auto g_yyy_xxxyy = braBuffer.data(goff + 129 * kcomp + j);

                    auto g_yyy_xxxyz = braBuffer.data(goff + 130 * kcomp + j);

                    auto g_yyy_xxxzz = braBuffer.data(goff + 131 * kcomp + j);

                    auto g_yyy_xxyyy = braBuffer.data(goff + 132 * kcomp + j);

                    auto g_yyy_xxyyz = braBuffer.data(goff + 133 * kcomp + j);

                    auto g_yyy_xxyzz = braBuffer.data(goff + 134 * kcomp + j);

                    auto g_yyy_xxzzz = braBuffer.data(goff + 135 * kcomp + j);

                    auto g_yyy_xyyyy = braBuffer.data(goff + 136 * kcomp + j);

                    auto g_yyy_xyyyz = braBuffer.data(goff + 137 * kcomp + j);

                    auto g_yyy_xyyzz = braBuffer.data(goff + 138 * kcomp + j);

                    auto g_yyy_xyzzz = braBuffer.data(goff + 139 * kcomp + j);

                    auto g_yyy_xzzzz = braBuffer.data(goff + 140 * kcomp + j);

                    auto g_yyy_yyyyy = braBuffer.data(goff + 141 * kcomp + j);

                    auto g_yyy_yyyyz = braBuffer.data(goff + 142 * kcomp + j);

                    auto g_yyy_yyyzz = braBuffer.data(goff + 143 * kcomp + j);

                    auto g_yyy_yyzzz = braBuffer.data(goff + 144 * kcomp + j);

                    auto g_yyy_yzzzz = braBuffer.data(goff + 145 * kcomp + j);

                    auto g_yyy_zzzzz = braBuffer.data(goff + 146 * kcomp + j);

                    auto g_yyz_xxxxx = braBuffer.data(goff + 147 * kcomp + j);

                    auto g_yyz_xxxxy = braBuffer.data(goff + 148 * kcomp + j);

                    auto g_yyz_xxxxz = braBuffer.data(goff + 149 * kcomp + j);

                    auto g_yyz_xxxyy = braBuffer.data(goff + 150 * kcomp + j);

                    auto g_yyz_xxxyz = braBuffer.data(goff + 151 * kcomp + j);

                    auto g_yyz_xxxzz = braBuffer.data(goff + 152 * kcomp + j);

                    auto g_yyz_xxyyy = braBuffer.data(goff + 153 * kcomp + j);

                    auto g_yyz_xxyyz = braBuffer.data(goff + 154 * kcomp + j);

                    auto g_yyz_xxyzz = braBuffer.data(goff + 155 * kcomp + j);

                    auto g_yyz_xxzzz = braBuffer.data(goff + 156 * kcomp + j);

                    auto g_yyz_xyyyy = braBuffer.data(goff + 157 * kcomp + j);

                    auto g_yyz_xyyyz = braBuffer.data(goff + 158 * kcomp + j);

                    auto g_yyz_xyyzz = braBuffer.data(goff + 159 * kcomp + j);

                    auto g_yyz_xyzzz = braBuffer.data(goff + 160 * kcomp + j);

                    auto g_yyz_xzzzz = braBuffer.data(goff + 161 * kcomp + j);

                    auto g_yyz_yyyyy = braBuffer.data(goff + 162 * kcomp + j);

                    auto g_yyz_yyyyz = braBuffer.data(goff + 163 * kcomp + j);

                    auto g_yyz_yyyzz = braBuffer.data(goff + 164 * kcomp + j);

                    auto g_yyz_yyzzz = braBuffer.data(goff + 165 * kcomp + j);

                    auto g_yyz_yzzzz = braBuffer.data(goff + 166 * kcomp + j);

                    auto g_yyz_zzzzz = braBuffer.data(goff + 167 * kcomp + j);

                    auto g_yzz_xxxxx = braBuffer.data(goff + 168 * kcomp + j);

                    auto g_yzz_xxxxy = braBuffer.data(goff + 169 * kcomp + j);

                    auto g_yzz_xxxxz = braBuffer.data(goff + 170 * kcomp + j);

                    auto g_yzz_xxxyy = braBuffer.data(goff + 171 * kcomp + j);

                    auto g_yzz_xxxyz = braBuffer.data(goff + 172 * kcomp + j);

                    auto g_yzz_xxxzz = braBuffer.data(goff + 173 * kcomp + j);

                    auto g_yzz_xxyyy = braBuffer.data(goff + 174 * kcomp + j);

                    auto g_yzz_xxyyz = braBuffer.data(goff + 175 * kcomp + j);

                    auto g_yzz_xxyzz = braBuffer.data(goff + 176 * kcomp + j);

                    auto g_yzz_xxzzz = braBuffer.data(goff + 177 * kcomp + j);

                    auto g_yzz_xyyyy = braBuffer.data(goff + 178 * kcomp + j);

                    auto g_yzz_xyyyz = braBuffer.data(goff + 179 * kcomp + j);

                    auto g_yzz_xyyzz = braBuffer.data(goff + 180 * kcomp + j);

                    auto g_yzz_xyzzz = braBuffer.data(goff + 181 * kcomp + j);

                    auto g_yzz_xzzzz = braBuffer.data(goff + 182 * kcomp + j);

                    auto g_yzz_yyyyy = braBuffer.data(goff + 183 * kcomp + j);

                    auto g_yzz_yyyyz = braBuffer.data(goff + 184 * kcomp + j);

                    auto g_yzz_yyyzz = braBuffer.data(goff + 185 * kcomp + j);

                    auto g_yzz_yyzzz = braBuffer.data(goff + 186 * kcomp + j);

                    auto g_yzz_yzzzz = braBuffer.data(goff + 187 * kcomp + j);

                    auto g_yzz_zzzzz = braBuffer.data(goff + 188 * kcomp + j);

                    auto g_zzz_xxxxx = braBuffer.data(goff + 189 * kcomp + j);

                    auto g_zzz_xxxxy = braBuffer.data(goff + 190 * kcomp + j);

                    auto g_zzz_xxxxz = braBuffer.data(goff + 191 * kcomp + j);

                    auto g_zzz_xxxyy = braBuffer.data(goff + 192 * kcomp + j);

                    auto g_zzz_xxxyz = braBuffer.data(goff + 193 * kcomp + j);

                    auto g_zzz_xxxzz = braBuffer.data(goff + 194 * kcomp + j);

                    auto g_zzz_xxyyy = braBuffer.data(goff + 195 * kcomp + j);

                    auto g_zzz_xxyyz = braBuffer.data(goff + 196 * kcomp + j);

                    auto g_zzz_xxyzz = braBuffer.data(goff + 197 * kcomp + j);

                    auto g_zzz_xxzzz = braBuffer.data(goff + 198 * kcomp + j);

                    auto g_zzz_xyyyy = braBuffer.data(goff + 199 * kcomp + j);

                    auto g_zzz_xyyyz = braBuffer.data(goff + 200 * kcomp + j);

                    auto g_zzz_xyyzz = braBuffer.data(goff + 201 * kcomp + j);

                    auto g_zzz_xyzzz = braBuffer.data(goff + 202 * kcomp + j);

                    auto g_zzz_xzzzz = braBuffer.data(goff + 203 * kcomp + j);

                    auto g_zzz_yyyyy = braBuffer.data(goff + 204 * kcomp + j);

                    auto g_zzz_yyyyz = braBuffer.data(goff + 205 * kcomp + j);

                    auto g_zzz_yyyzz = braBuffer.data(goff + 206 * kcomp + j);

                    auto g_zzz_yyzzz = braBuffer.data(goff + 207 * kcomp + j);

                    auto g_zzz_yzzzz = braBuffer.data(goff + 208 * kcomp + j);

                    auto g_zzz_zzzzz = braBuffer.data(goff + 209 * kcomp + j);

                    #pragma omp simd aligned(g2_xx_xxxxx, g2_xx_xxxxy, g2_xx_xxxxz,\
                                             g2_xx_xxxyy, g2_xx_xxxyz, g2_xx_xxxzz,\
                                             g2_xx_xxyyy, g2_xx_xxyyz, g2_xx_xxyzz,\
                                             g2_xx_xxzzz, g2_xx_xyyyy, g2_xx_xyyyz,\
                                             g2_xx_xyyzz, g2_xx_xyzzz, g2_xx_xzzzz,\
                                             g2_xx_yyyyy, g2_xx_yyyyz, g2_xx_yyyzz,\
                                             g2_xx_yyzzz, g2_xx_yzzzz, g2_xx_zzzzz,\
                                             g2_xy_xxxxx, g2_xy_xxxxy, g2_xy_xxxxz,\
                                             g2_xy_xxxyy, g2_xy_xxxyz, g2_xy_xxxzz,\
                                             g2_xy_xxyyy, g2_xy_xxyyz, g2_xy_xxyzz,\
                                             g2_xy_xxzzz, g2_xy_xyyyy, g2_xy_xyyyz,\
                                             g2_xy_xyyzz, g2_xy_xyzzz, g2_xy_xzzzz,\
                                             g2_xy_yyyyy, g2_xy_yyyyz, g2_xy_yyyzz,\
                                             g2_xy_yyzzz, g2_xy_yzzzz, g2_xy_zzzzz,\
                                             g2_xz_xxxxx, g2_xz_xxxxy, g2_xz_xxxxz,\
                                             g2_xz_xxxyy, g2_xz_xxxyz, g2_xz_xxxzz,\
                                             g2_xz_xxyyy, g2_xz_xxyyz, g2_xz_xxyzz,\
                                             g2_xz_xxzzz, g2_xz_xyyyy, g2_xz_xyyyz,\
                                             g2_xz_xyyzz, g2_xz_xyzzz, g2_xz_xzzzz,\
                                             g2_xz_yyyyy, g2_xz_yyyyz, g2_xz_yyyzz,\
                                             g2_xz_yyzzz, g2_xz_yzzzz, g2_xz_zzzzz,\
                                             g2_yy_xxxxx, g2_yy_xxxxy, g2_yy_xxxxz,\
                                             g2_yy_xxxyy, g2_yy_xxxyz, g2_yy_xxxzz,\
                                             g2_yy_xxyyy, g2_yy_xxyyz, g2_yy_xxyzz,\
                                             g2_yy_xxzzz, g2_yy_xyyyy, g2_yy_xyyyz,\
                                             g2_yy_xyyzz, g2_yy_xyzzz, g2_yy_xzzzz,\
                                             g2_yy_yyyyy, g2_yy_yyyyz, g2_yy_yyyzz,\
                                             g2_yy_yyzzz, g2_yy_yzzzz, g2_yy_zzzzz,\
                                             g2_yz_xxxxx, g2_yz_xxxxy, g2_yz_xxxxz,\
                                             g2_yz_xxxyy, g2_yz_xxxyz, g2_yz_xxxzz,\
                                             g2_yz_xxyyy, g2_yz_xxyyz, g2_yz_xxyzz,\
                                             g2_yz_xxzzz, g2_yz_xyyyy, g2_yz_xyyyz,\
                                             g2_yz_xyyzz, g2_yz_xyzzz, g2_yz_xzzzz,\
                                             g2_yz_yyyyy, g2_yz_yyyyz, g2_yz_yyyzz,\
                                             g2_yz_yyzzz, g2_yz_yzzzz, g2_yz_zzzzz,\
                                             g2_zz_xxxxx, g2_zz_xxxxy, g2_zz_xxxxz,\
                                             g2_zz_xxxyy, g2_zz_xxxyz, g2_zz_xxxzz,\
                                             g2_zz_xxyyy, g2_zz_xxyyz, g2_zz_xxyzz,\
                                             g2_zz_xxzzz, g2_zz_xyyyy, g2_zz_xyyyz,\
                                             g2_zz_xyyzz, g2_zz_xyzzz, g2_zz_xzzzz,\
                                             g2_zz_yyyyy, g2_zz_yyyyz, g2_zz_yyyzz,\
                                             g2_zz_yyzzz, g2_zz_yzzzz, g2_zz_zzzzz,\
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
                                             g_xxx_xxxxx, g_xxx_xxxxy, g_xxx_xxxxz,\
                                             g_xxx_xxxyy, g_xxx_xxxyz, g_xxx_xxxzz,\
                                             g_xxx_xxyyy, g_xxx_xxyyz, g_xxx_xxyzz,\
                                             g_xxx_xxzzz, g_xxx_xyyyy, g_xxx_xyyyz,\
                                             g_xxx_xyyzz, g_xxx_xyzzz, g_xxx_xzzzz,\
                                             g_xxx_yyyyy, g_xxx_yyyyz, g_xxx_yyyzz,\
                                             g_xxx_yyzzz, g_xxx_yzzzz, g_xxx_zzzzz,\
                                             g_xxy_xxxxx, g_xxy_xxxxy, g_xxy_xxxxz,\
                                             g_xxy_xxxyy, g_xxy_xxxyz, g_xxy_xxxzz,\
                                             g_xxy_xxyyy, g_xxy_xxyyz, g_xxy_xxyzz,\
                                             g_xxy_xxzzz, g_xxy_xyyyy, g_xxy_xyyyz,\
                                             g_xxy_xyyzz, g_xxy_xyzzz, g_xxy_xzzzz,\
                                             g_xxy_yyyyy, g_xxy_yyyyz, g_xxy_yyyzz,\
                                             g_xxy_yyzzz, g_xxy_yzzzz, g_xxy_zzzzz,\
                                             g_xxz_xxxxx, g_xxz_xxxxy, g_xxz_xxxxz,\
                                             g_xxz_xxxyy, g_xxz_xxxyz, g_xxz_xxxzz,\
                                             g_xxz_xxyyy, g_xxz_xxyyz, g_xxz_xxyzz,\
                                             g_xxz_xxzzz, g_xxz_xyyyy, g_xxz_xyyyz,\
                                             g_xxz_xyyzz, g_xxz_xyzzz, g_xxz_xzzzz,\
                                             g_xxz_yyyyy, g_xxz_yyyyz, g_xxz_yyyzz,\
                                             g_xxz_yyzzz, g_xxz_yzzzz, g_xxz_zzzzz,\
                                             g_xyy_xxxxx, g_xyy_xxxxy, g_xyy_xxxxz,\
                                             g_xyy_xxxyy, g_xyy_xxxyz, g_xyy_xxxzz,\
                                             g_xyy_xxyyy, g_xyy_xxyyz, g_xyy_xxyzz,\
                                             g_xyy_xxzzz, g_xyy_xyyyy, g_xyy_xyyyz,\
                                             g_xyy_xyyzz, g_xyy_xyzzz, g_xyy_xzzzz,\
                                             g_xyy_yyyyy, g_xyy_yyyyz, g_xyy_yyyzz,\
                                             g_xyy_yyzzz, g_xyy_yzzzz, g_xyy_zzzzz,\
                                             g_xyz_xxxxx, g_xyz_xxxxy, g_xyz_xxxxz,\
                                             g_xyz_xxxyy, g_xyz_xxxyz, g_xyz_xxxzz,\
                                             g_xyz_xxyyy, g_xyz_xxyyz, g_xyz_xxyzz,\
                                             g_xyz_xxzzz, g_xyz_xyyyy, g_xyz_xyyyz,\
                                             g_xyz_xyyzz, g_xyz_xyzzz, g_xyz_xzzzz,\
                                             g_xyz_yyyyy, g_xyz_yyyyz, g_xyz_yyyzz,\
                                             g_xyz_yyzzz, g_xyz_yzzzz, g_xyz_zzzzz,\
                                             g_xzz_xxxxx, g_xzz_xxxxy, g_xzz_xxxxz,\
                                             g_xzz_xxxyy, g_xzz_xxxyz, g_xzz_xxxzz,\
                                             g_xzz_xxyyy, g_xzz_xxyyz, g_xzz_xxyzz,\
                                             g_xzz_xxzzz, g_xzz_xyyyy, g_xzz_xyyyz,\
                                             g_xzz_xyyzz, g_xzz_xyzzz, g_xzz_xzzzz,\
                                             g_xzz_yyyyy, g_xzz_yyyyz, g_xzz_yyyzz,\
                                             g_xzz_yyzzz, g_xzz_yzzzz, g_xzz_zzzzz,\
                                             g_yyy_xxxxx, g_yyy_xxxxy, g_yyy_xxxxz,\
                                             g_yyy_xxxyy, g_yyy_xxxyz, g_yyy_xxxzz,\
                                             g_yyy_xxyyy, g_yyy_xxyyz, g_yyy_xxyzz,\
                                             g_yyy_xxzzz, g_yyy_xyyyy, g_yyy_xyyyz,\
                                             g_yyy_xyyzz, g_yyy_xyzzz, g_yyy_xzzzz,\
                                             g_yyy_yyyyy, g_yyy_yyyyz, g_yyy_yyyzz,\
                                             g_yyy_yyzzz, g_yyy_yzzzz, g_yyy_zzzzz,\
                                             g_yyz_xxxxx, g_yyz_xxxxy, g_yyz_xxxxz,\
                                             g_yyz_xxxyy, g_yyz_xxxyz, g_yyz_xxxzz,\
                                             g_yyz_xxyyy, g_yyz_xxyyz, g_yyz_xxyzz,\
                                             g_yyz_xxzzz, g_yyz_xyyyy, g_yyz_xyyyz,\
                                             g_yyz_xyyzz, g_yyz_xyzzz, g_yyz_xzzzz,\
                                             g_yyz_yyyyy, g_yyz_yyyyz, g_yyz_yyyzz,\
                                             g_yyz_yyzzz, g_yyz_yzzzz, g_yyz_zzzzz,\
                                             g_yzz_xxxxx, g_yzz_xxxxy, g_yzz_xxxxz,\
                                             g_yzz_xxxyy, g_yzz_xxxyz, g_yzz_xxxzz,\
                                             g_yzz_xxyyy, g_yzz_xxyyz, g_yzz_xxyzz,\
                                             g_yzz_xxzzz, g_yzz_xyyyy, g_yzz_xyyyz,\
                                             g_yzz_xyyzz, g_yzz_xyzzz, g_yzz_xzzzz,\
                                             g_yzz_yyyyy, g_yzz_yyyyz, g_yzz_yyyzz,\
                                             g_yzz_yyzzz, g_yzz_yzzzz, g_yzz_zzzzz,\
                                             g_zzz_xxxxx, g_zzz_xxxxy, g_zzz_xxxxz,\
                                             g_zzz_xxxyy, g_zzz_xxxyz, g_zzz_xxxzz,\
                                             g_zzz_xxyyy, g_zzz_xxyyz, g_zzz_xxyzz,\
                                             g_zzz_xxzzz, g_zzz_xyyyy, g_zzz_xyyyz,\
                                             g_zzz_xyyzz, g_zzz_xyzzz, g_zzz_xzzzz,\
                                             g_zzz_yyyyy, g_zzz_yyyyz, g_zzz_yyyzz,\
                                             g_zzz_yyzzz, g_zzz_yzzzz, g_zzz_zzzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xxx_xxxxx[k] = g1_xx_xxxxxx[k] - abx * g2_xx_xxxxx[k];

                        g_xxx_xxxxy[k] = g1_xx_xxxxxy[k] - abx * g2_xx_xxxxy[k];

                        g_xxx_xxxxz[k] = g1_xx_xxxxxz[k] - abx * g2_xx_xxxxz[k];

                        g_xxx_xxxyy[k] = g1_xx_xxxxyy[k] - abx * g2_xx_xxxyy[k];

                        g_xxx_xxxyz[k] = g1_xx_xxxxyz[k] - abx * g2_xx_xxxyz[k];

                        g_xxx_xxxzz[k] = g1_xx_xxxxzz[k] - abx * g2_xx_xxxzz[k];

                        g_xxx_xxyyy[k] = g1_xx_xxxyyy[k] - abx * g2_xx_xxyyy[k];

                        g_xxx_xxyyz[k] = g1_xx_xxxyyz[k] - abx * g2_xx_xxyyz[k];

                        g_xxx_xxyzz[k] = g1_xx_xxxyzz[k] - abx * g2_xx_xxyzz[k];

                        g_xxx_xxzzz[k] = g1_xx_xxxzzz[k] - abx * g2_xx_xxzzz[k];

                        g_xxx_xyyyy[k] = g1_xx_xxyyyy[k] - abx * g2_xx_xyyyy[k];

                        g_xxx_xyyyz[k] = g1_xx_xxyyyz[k] - abx * g2_xx_xyyyz[k];

                        g_xxx_xyyzz[k] = g1_xx_xxyyzz[k] - abx * g2_xx_xyyzz[k];

                        g_xxx_xyzzz[k] = g1_xx_xxyzzz[k] - abx * g2_xx_xyzzz[k];

                        g_xxx_xzzzz[k] = g1_xx_xxzzzz[k] - abx * g2_xx_xzzzz[k];

                        g_xxx_yyyyy[k] = g1_xx_xyyyyy[k] - abx * g2_xx_yyyyy[k];

                        g_xxx_yyyyz[k] = g1_xx_xyyyyz[k] - abx * g2_xx_yyyyz[k];

                        g_xxx_yyyzz[k] = g1_xx_xyyyzz[k] - abx * g2_xx_yyyzz[k];

                        g_xxx_yyzzz[k] = g1_xx_xyyzzz[k] - abx * g2_xx_yyzzz[k];

                        g_xxx_yzzzz[k] = g1_xx_xyzzzz[k] - abx * g2_xx_yzzzz[k];

                        g_xxx_zzzzz[k] = g1_xx_xzzzzz[k] - abx * g2_xx_zzzzz[k];

                        g_xxy_xxxxx[k] = g1_xy_xxxxxx[k] - abx * g2_xy_xxxxx[k];

                        g_xxy_xxxxy[k] = g1_xy_xxxxxy[k] - abx * g2_xy_xxxxy[k];

                        g_xxy_xxxxz[k] = g1_xy_xxxxxz[k] - abx * g2_xy_xxxxz[k];

                        g_xxy_xxxyy[k] = g1_xy_xxxxyy[k] - abx * g2_xy_xxxyy[k];

                        g_xxy_xxxyz[k] = g1_xy_xxxxyz[k] - abx * g2_xy_xxxyz[k];

                        g_xxy_xxxzz[k] = g1_xy_xxxxzz[k] - abx * g2_xy_xxxzz[k];

                        g_xxy_xxyyy[k] = g1_xy_xxxyyy[k] - abx * g2_xy_xxyyy[k];

                        g_xxy_xxyyz[k] = g1_xy_xxxyyz[k] - abx * g2_xy_xxyyz[k];

                        g_xxy_xxyzz[k] = g1_xy_xxxyzz[k] - abx * g2_xy_xxyzz[k];

                        g_xxy_xxzzz[k] = g1_xy_xxxzzz[k] - abx * g2_xy_xxzzz[k];

                        g_xxy_xyyyy[k] = g1_xy_xxyyyy[k] - abx * g2_xy_xyyyy[k];

                        g_xxy_xyyyz[k] = g1_xy_xxyyyz[k] - abx * g2_xy_xyyyz[k];

                        g_xxy_xyyzz[k] = g1_xy_xxyyzz[k] - abx * g2_xy_xyyzz[k];

                        g_xxy_xyzzz[k] = g1_xy_xxyzzz[k] - abx * g2_xy_xyzzz[k];

                        g_xxy_xzzzz[k] = g1_xy_xxzzzz[k] - abx * g2_xy_xzzzz[k];

                        g_xxy_yyyyy[k] = g1_xy_xyyyyy[k] - abx * g2_xy_yyyyy[k];

                        g_xxy_yyyyz[k] = g1_xy_xyyyyz[k] - abx * g2_xy_yyyyz[k];

                        g_xxy_yyyzz[k] = g1_xy_xyyyzz[k] - abx * g2_xy_yyyzz[k];

                        g_xxy_yyzzz[k] = g1_xy_xyyzzz[k] - abx * g2_xy_yyzzz[k];

                        g_xxy_yzzzz[k] = g1_xy_xyzzzz[k] - abx * g2_xy_yzzzz[k];

                        g_xxy_zzzzz[k] = g1_xy_xzzzzz[k] - abx * g2_xy_zzzzz[k];

                        g_xxz_xxxxx[k] = g1_xz_xxxxxx[k] - abx * g2_xz_xxxxx[k];

                        g_xxz_xxxxy[k] = g1_xz_xxxxxy[k] - abx * g2_xz_xxxxy[k];

                        g_xxz_xxxxz[k] = g1_xz_xxxxxz[k] - abx * g2_xz_xxxxz[k];

                        g_xxz_xxxyy[k] = g1_xz_xxxxyy[k] - abx * g2_xz_xxxyy[k];

                        g_xxz_xxxyz[k] = g1_xz_xxxxyz[k] - abx * g2_xz_xxxyz[k];

                        g_xxz_xxxzz[k] = g1_xz_xxxxzz[k] - abx * g2_xz_xxxzz[k];

                        g_xxz_xxyyy[k] = g1_xz_xxxyyy[k] - abx * g2_xz_xxyyy[k];

                        g_xxz_xxyyz[k] = g1_xz_xxxyyz[k] - abx * g2_xz_xxyyz[k];

                        g_xxz_xxyzz[k] = g1_xz_xxxyzz[k] - abx * g2_xz_xxyzz[k];

                        g_xxz_xxzzz[k] = g1_xz_xxxzzz[k] - abx * g2_xz_xxzzz[k];

                        g_xxz_xyyyy[k] = g1_xz_xxyyyy[k] - abx * g2_xz_xyyyy[k];

                        g_xxz_xyyyz[k] = g1_xz_xxyyyz[k] - abx * g2_xz_xyyyz[k];

                        g_xxz_xyyzz[k] = g1_xz_xxyyzz[k] - abx * g2_xz_xyyzz[k];

                        g_xxz_xyzzz[k] = g1_xz_xxyzzz[k] - abx * g2_xz_xyzzz[k];

                        g_xxz_xzzzz[k] = g1_xz_xxzzzz[k] - abx * g2_xz_xzzzz[k];

                        g_xxz_yyyyy[k] = g1_xz_xyyyyy[k] - abx * g2_xz_yyyyy[k];

                        g_xxz_yyyyz[k] = g1_xz_xyyyyz[k] - abx * g2_xz_yyyyz[k];

                        g_xxz_yyyzz[k] = g1_xz_xyyyzz[k] - abx * g2_xz_yyyzz[k];

                        g_xxz_yyzzz[k] = g1_xz_xyyzzz[k] - abx * g2_xz_yyzzz[k];

                        g_xxz_yzzzz[k] = g1_xz_xyzzzz[k] - abx * g2_xz_yzzzz[k];

                        g_xxz_zzzzz[k] = g1_xz_xzzzzz[k] - abx * g2_xz_zzzzz[k];

                        g_xyy_xxxxx[k] = g1_yy_xxxxxx[k] - abx * g2_yy_xxxxx[k];

                        g_xyy_xxxxy[k] = g1_yy_xxxxxy[k] - abx * g2_yy_xxxxy[k];

                        g_xyy_xxxxz[k] = g1_yy_xxxxxz[k] - abx * g2_yy_xxxxz[k];

                        g_xyy_xxxyy[k] = g1_yy_xxxxyy[k] - abx * g2_yy_xxxyy[k];

                        g_xyy_xxxyz[k] = g1_yy_xxxxyz[k] - abx * g2_yy_xxxyz[k];

                        g_xyy_xxxzz[k] = g1_yy_xxxxzz[k] - abx * g2_yy_xxxzz[k];

                        g_xyy_xxyyy[k] = g1_yy_xxxyyy[k] - abx * g2_yy_xxyyy[k];

                        g_xyy_xxyyz[k] = g1_yy_xxxyyz[k] - abx * g2_yy_xxyyz[k];

                        g_xyy_xxyzz[k] = g1_yy_xxxyzz[k] - abx * g2_yy_xxyzz[k];

                        g_xyy_xxzzz[k] = g1_yy_xxxzzz[k] - abx * g2_yy_xxzzz[k];

                        g_xyy_xyyyy[k] = g1_yy_xxyyyy[k] - abx * g2_yy_xyyyy[k];

                        g_xyy_xyyyz[k] = g1_yy_xxyyyz[k] - abx * g2_yy_xyyyz[k];

                        g_xyy_xyyzz[k] = g1_yy_xxyyzz[k] - abx * g2_yy_xyyzz[k];

                        g_xyy_xyzzz[k] = g1_yy_xxyzzz[k] - abx * g2_yy_xyzzz[k];

                        g_xyy_xzzzz[k] = g1_yy_xxzzzz[k] - abx * g2_yy_xzzzz[k];

                        g_xyy_yyyyy[k] = g1_yy_xyyyyy[k] - abx * g2_yy_yyyyy[k];

                        g_xyy_yyyyz[k] = g1_yy_xyyyyz[k] - abx * g2_yy_yyyyz[k];

                        g_xyy_yyyzz[k] = g1_yy_xyyyzz[k] - abx * g2_yy_yyyzz[k];

                        g_xyy_yyzzz[k] = g1_yy_xyyzzz[k] - abx * g2_yy_yyzzz[k];

                        g_xyy_yzzzz[k] = g1_yy_xyzzzz[k] - abx * g2_yy_yzzzz[k];

                        g_xyy_zzzzz[k] = g1_yy_xzzzzz[k] - abx * g2_yy_zzzzz[k];

                        g_xyz_xxxxx[k] = g1_yz_xxxxxx[k] - abx * g2_yz_xxxxx[k];

                        g_xyz_xxxxy[k] = g1_yz_xxxxxy[k] - abx * g2_yz_xxxxy[k];

                        g_xyz_xxxxz[k] = g1_yz_xxxxxz[k] - abx * g2_yz_xxxxz[k];

                        g_xyz_xxxyy[k] = g1_yz_xxxxyy[k] - abx * g2_yz_xxxyy[k];

                        g_xyz_xxxyz[k] = g1_yz_xxxxyz[k] - abx * g2_yz_xxxyz[k];

                        g_xyz_xxxzz[k] = g1_yz_xxxxzz[k] - abx * g2_yz_xxxzz[k];

                        g_xyz_xxyyy[k] = g1_yz_xxxyyy[k] - abx * g2_yz_xxyyy[k];

                        g_xyz_xxyyz[k] = g1_yz_xxxyyz[k] - abx * g2_yz_xxyyz[k];

                        g_xyz_xxyzz[k] = g1_yz_xxxyzz[k] - abx * g2_yz_xxyzz[k];

                        g_xyz_xxzzz[k] = g1_yz_xxxzzz[k] - abx * g2_yz_xxzzz[k];

                        g_xyz_xyyyy[k] = g1_yz_xxyyyy[k] - abx * g2_yz_xyyyy[k];

                        g_xyz_xyyyz[k] = g1_yz_xxyyyz[k] - abx * g2_yz_xyyyz[k];

                        g_xyz_xyyzz[k] = g1_yz_xxyyzz[k] - abx * g2_yz_xyyzz[k];

                        g_xyz_xyzzz[k] = g1_yz_xxyzzz[k] - abx * g2_yz_xyzzz[k];

                        g_xyz_xzzzz[k] = g1_yz_xxzzzz[k] - abx * g2_yz_xzzzz[k];

                        g_xyz_yyyyy[k] = g1_yz_xyyyyy[k] - abx * g2_yz_yyyyy[k];

                        g_xyz_yyyyz[k] = g1_yz_xyyyyz[k] - abx * g2_yz_yyyyz[k];

                        g_xyz_yyyzz[k] = g1_yz_xyyyzz[k] - abx * g2_yz_yyyzz[k];

                        g_xyz_yyzzz[k] = g1_yz_xyyzzz[k] - abx * g2_yz_yyzzz[k];

                        g_xyz_yzzzz[k] = g1_yz_xyzzzz[k] - abx * g2_yz_yzzzz[k];

                        g_xyz_zzzzz[k] = g1_yz_xzzzzz[k] - abx * g2_yz_zzzzz[k];

                        g_xzz_xxxxx[k] = g1_zz_xxxxxx[k] - abx * g2_zz_xxxxx[k];

                        g_xzz_xxxxy[k] = g1_zz_xxxxxy[k] - abx * g2_zz_xxxxy[k];

                        g_xzz_xxxxz[k] = g1_zz_xxxxxz[k] - abx * g2_zz_xxxxz[k];

                        g_xzz_xxxyy[k] = g1_zz_xxxxyy[k] - abx * g2_zz_xxxyy[k];

                        g_xzz_xxxyz[k] = g1_zz_xxxxyz[k] - abx * g2_zz_xxxyz[k];

                        g_xzz_xxxzz[k] = g1_zz_xxxxzz[k] - abx * g2_zz_xxxzz[k];

                        g_xzz_xxyyy[k] = g1_zz_xxxyyy[k] - abx * g2_zz_xxyyy[k];

                        g_xzz_xxyyz[k] = g1_zz_xxxyyz[k] - abx * g2_zz_xxyyz[k];

                        g_xzz_xxyzz[k] = g1_zz_xxxyzz[k] - abx * g2_zz_xxyzz[k];

                        g_xzz_xxzzz[k] = g1_zz_xxxzzz[k] - abx * g2_zz_xxzzz[k];

                        g_xzz_xyyyy[k] = g1_zz_xxyyyy[k] - abx * g2_zz_xyyyy[k];

                        g_xzz_xyyyz[k] = g1_zz_xxyyyz[k] - abx * g2_zz_xyyyz[k];

                        g_xzz_xyyzz[k] = g1_zz_xxyyzz[k] - abx * g2_zz_xyyzz[k];

                        g_xzz_xyzzz[k] = g1_zz_xxyzzz[k] - abx * g2_zz_xyzzz[k];

                        g_xzz_xzzzz[k] = g1_zz_xxzzzz[k] - abx * g2_zz_xzzzz[k];

                        g_xzz_yyyyy[k] = g1_zz_xyyyyy[k] - abx * g2_zz_yyyyy[k];

                        g_xzz_yyyyz[k] = g1_zz_xyyyyz[k] - abx * g2_zz_yyyyz[k];

                        g_xzz_yyyzz[k] = g1_zz_xyyyzz[k] - abx * g2_zz_yyyzz[k];

                        g_xzz_yyzzz[k] = g1_zz_xyyzzz[k] - abx * g2_zz_yyzzz[k];

                        g_xzz_yzzzz[k] = g1_zz_xyzzzz[k] - abx * g2_zz_yzzzz[k];

                        g_xzz_zzzzz[k] = g1_zz_xzzzzz[k] - abx * g2_zz_zzzzz[k];

                        // leading y component

                        g_yyy_xxxxx[k] = g1_yy_xxxxxy[k] - aby * g2_yy_xxxxx[k];

                        g_yyy_xxxxy[k] = g1_yy_xxxxyy[k] - aby * g2_yy_xxxxy[k];

                        g_yyy_xxxxz[k] = g1_yy_xxxxyz[k] - aby * g2_yy_xxxxz[k];

                        g_yyy_xxxyy[k] = g1_yy_xxxyyy[k] - aby * g2_yy_xxxyy[k];

                        g_yyy_xxxyz[k] = g1_yy_xxxyyz[k] - aby * g2_yy_xxxyz[k];

                        g_yyy_xxxzz[k] = g1_yy_xxxyzz[k] - aby * g2_yy_xxxzz[k];

                        g_yyy_xxyyy[k] = g1_yy_xxyyyy[k] - aby * g2_yy_xxyyy[k];

                        g_yyy_xxyyz[k] = g1_yy_xxyyyz[k] - aby * g2_yy_xxyyz[k];

                        g_yyy_xxyzz[k] = g1_yy_xxyyzz[k] - aby * g2_yy_xxyzz[k];

                        g_yyy_xxzzz[k] = g1_yy_xxyzzz[k] - aby * g2_yy_xxzzz[k];

                        g_yyy_xyyyy[k] = g1_yy_xyyyyy[k] - aby * g2_yy_xyyyy[k];

                        g_yyy_xyyyz[k] = g1_yy_xyyyyz[k] - aby * g2_yy_xyyyz[k];

                        g_yyy_xyyzz[k] = g1_yy_xyyyzz[k] - aby * g2_yy_xyyzz[k];

                        g_yyy_xyzzz[k] = g1_yy_xyyzzz[k] - aby * g2_yy_xyzzz[k];

                        g_yyy_xzzzz[k] = g1_yy_xyzzzz[k] - aby * g2_yy_xzzzz[k];

                        g_yyy_yyyyy[k] = g1_yy_yyyyyy[k] - aby * g2_yy_yyyyy[k];

                        g_yyy_yyyyz[k] = g1_yy_yyyyyz[k] - aby * g2_yy_yyyyz[k];

                        g_yyy_yyyzz[k] = g1_yy_yyyyzz[k] - aby * g2_yy_yyyzz[k];

                        g_yyy_yyzzz[k] = g1_yy_yyyzzz[k] - aby * g2_yy_yyzzz[k];

                        g_yyy_yzzzz[k] = g1_yy_yyzzzz[k] - aby * g2_yy_yzzzz[k];

                        g_yyy_zzzzz[k] = g1_yy_yzzzzz[k] - aby * g2_yy_zzzzz[k];

                        g_yyz_xxxxx[k] = g1_yz_xxxxxy[k] - aby * g2_yz_xxxxx[k];

                        g_yyz_xxxxy[k] = g1_yz_xxxxyy[k] - aby * g2_yz_xxxxy[k];

                        g_yyz_xxxxz[k] = g1_yz_xxxxyz[k] - aby * g2_yz_xxxxz[k];

                        g_yyz_xxxyy[k] = g1_yz_xxxyyy[k] - aby * g2_yz_xxxyy[k];

                        g_yyz_xxxyz[k] = g1_yz_xxxyyz[k] - aby * g2_yz_xxxyz[k];

                        g_yyz_xxxzz[k] = g1_yz_xxxyzz[k] - aby * g2_yz_xxxzz[k];

                        g_yyz_xxyyy[k] = g1_yz_xxyyyy[k] - aby * g2_yz_xxyyy[k];

                        g_yyz_xxyyz[k] = g1_yz_xxyyyz[k] - aby * g2_yz_xxyyz[k];

                        g_yyz_xxyzz[k] = g1_yz_xxyyzz[k] - aby * g2_yz_xxyzz[k];

                        g_yyz_xxzzz[k] = g1_yz_xxyzzz[k] - aby * g2_yz_xxzzz[k];

                        g_yyz_xyyyy[k] = g1_yz_xyyyyy[k] - aby * g2_yz_xyyyy[k];

                        g_yyz_xyyyz[k] = g1_yz_xyyyyz[k] - aby * g2_yz_xyyyz[k];

                        g_yyz_xyyzz[k] = g1_yz_xyyyzz[k] - aby * g2_yz_xyyzz[k];

                        g_yyz_xyzzz[k] = g1_yz_xyyzzz[k] - aby * g2_yz_xyzzz[k];

                        g_yyz_xzzzz[k] = g1_yz_xyzzzz[k] - aby * g2_yz_xzzzz[k];

                        g_yyz_yyyyy[k] = g1_yz_yyyyyy[k] - aby * g2_yz_yyyyy[k];

                        g_yyz_yyyyz[k] = g1_yz_yyyyyz[k] - aby * g2_yz_yyyyz[k];

                        g_yyz_yyyzz[k] = g1_yz_yyyyzz[k] - aby * g2_yz_yyyzz[k];

                        g_yyz_yyzzz[k] = g1_yz_yyyzzz[k] - aby * g2_yz_yyzzz[k];

                        g_yyz_yzzzz[k] = g1_yz_yyzzzz[k] - aby * g2_yz_yzzzz[k];

                        g_yyz_zzzzz[k] = g1_yz_yzzzzz[k] - aby * g2_yz_zzzzz[k];

                        g_yzz_xxxxx[k] = g1_zz_xxxxxy[k] - aby * g2_zz_xxxxx[k];

                        g_yzz_xxxxy[k] = g1_zz_xxxxyy[k] - aby * g2_zz_xxxxy[k];

                        g_yzz_xxxxz[k] = g1_zz_xxxxyz[k] - aby * g2_zz_xxxxz[k];

                        g_yzz_xxxyy[k] = g1_zz_xxxyyy[k] - aby * g2_zz_xxxyy[k];

                        g_yzz_xxxyz[k] = g1_zz_xxxyyz[k] - aby * g2_zz_xxxyz[k];

                        g_yzz_xxxzz[k] = g1_zz_xxxyzz[k] - aby * g2_zz_xxxzz[k];

                        g_yzz_xxyyy[k] = g1_zz_xxyyyy[k] - aby * g2_zz_xxyyy[k];

                        g_yzz_xxyyz[k] = g1_zz_xxyyyz[k] - aby * g2_zz_xxyyz[k];

                        g_yzz_xxyzz[k] = g1_zz_xxyyzz[k] - aby * g2_zz_xxyzz[k];

                        g_yzz_xxzzz[k] = g1_zz_xxyzzz[k] - aby * g2_zz_xxzzz[k];

                        g_yzz_xyyyy[k] = g1_zz_xyyyyy[k] - aby * g2_zz_xyyyy[k];

                        g_yzz_xyyyz[k] = g1_zz_xyyyyz[k] - aby * g2_zz_xyyyz[k];

                        g_yzz_xyyzz[k] = g1_zz_xyyyzz[k] - aby * g2_zz_xyyzz[k];

                        g_yzz_xyzzz[k] = g1_zz_xyyzzz[k] - aby * g2_zz_xyzzz[k];

                        g_yzz_xzzzz[k] = g1_zz_xyzzzz[k] - aby * g2_zz_xzzzz[k];

                        g_yzz_yyyyy[k] = g1_zz_yyyyyy[k] - aby * g2_zz_yyyyy[k];

                        g_yzz_yyyyz[k] = g1_zz_yyyyyz[k] - aby * g2_zz_yyyyz[k];

                        g_yzz_yyyzz[k] = g1_zz_yyyyzz[k] - aby * g2_zz_yyyzz[k];

                        g_yzz_yyzzz[k] = g1_zz_yyyzzz[k] - aby * g2_zz_yyzzz[k];

                        g_yzz_yzzzz[k] = g1_zz_yyzzzz[k] - aby * g2_zz_yzzzz[k];

                        g_yzz_zzzzz[k] = g1_zz_yzzzzz[k] - aby * g2_zz_zzzzz[k];

                        // leading z component

                        g_zzz_xxxxx[k] = g1_zz_xxxxxz[k] - abz * g2_zz_xxxxx[k];

                        g_zzz_xxxxy[k] = g1_zz_xxxxyz[k] - abz * g2_zz_xxxxy[k];

                        g_zzz_xxxxz[k] = g1_zz_xxxxzz[k] - abz * g2_zz_xxxxz[k];

                        g_zzz_xxxyy[k] = g1_zz_xxxyyz[k] - abz * g2_zz_xxxyy[k];

                        g_zzz_xxxyz[k] = g1_zz_xxxyzz[k] - abz * g2_zz_xxxyz[k];

                        g_zzz_xxxzz[k] = g1_zz_xxxzzz[k] - abz * g2_zz_xxxzz[k];

                        g_zzz_xxyyy[k] = g1_zz_xxyyyz[k] - abz * g2_zz_xxyyy[k];

                        g_zzz_xxyyz[k] = g1_zz_xxyyzz[k] - abz * g2_zz_xxyyz[k];

                        g_zzz_xxyzz[k] = g1_zz_xxyzzz[k] - abz * g2_zz_xxyzz[k];

                        g_zzz_xxzzz[k] = g1_zz_xxzzzz[k] - abz * g2_zz_xxzzz[k];

                        g_zzz_xyyyy[k] = g1_zz_xyyyyz[k] - abz * g2_zz_xyyyy[k];

                        g_zzz_xyyyz[k] = g1_zz_xyyyzz[k] - abz * g2_zz_xyyyz[k];

                        g_zzz_xyyzz[k] = g1_zz_xyyzzz[k] - abz * g2_zz_xyyzz[k];

                        g_zzz_xyzzz[k] = g1_zz_xyzzzz[k] - abz * g2_zz_xyzzz[k];

                        g_zzz_xzzzz[k] = g1_zz_xzzzzz[k] - abz * g2_zz_xzzzz[k];

                        g_zzz_yyyyy[k] = g1_zz_yyyyyz[k] - abz * g2_zz_yyyyy[k];

                        g_zzz_yyyyz[k] = g1_zz_yyyyzz[k] - abz * g2_zz_yyyyz[k];

                        g_zzz_yyyzz[k] = g1_zz_yyyzzz[k] - abz * g2_zz_yyyzz[k];

                        g_zzz_yyzzz[k] = g1_zz_yyzzzz[k] - abz * g2_zz_yyzzz[k];

                        g_zzz_yzzzz[k] = g1_zz_yzzzzz[k] - abz * g2_zz_yzzzz[k];

                        g_zzz_zzzzz[k] = g1_zz_zzzzzz[k] - abz * g2_zz_zzzzz[k];
                    }
                }
            }
        }
    }
    
    void
    compElectronRepulsionForGGXX(      CMemBlock2D<double>&  braBuffer,
                                 const CVecFourIndexes&      recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const bool                  isBraEqualKet,
                                 const int32_t               iContrPair)
    {
        // determine dimensions of ket side

        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

        if (isBraEqualKet) kdim  = iContrPair + 1;

        // set distances R(AB) = A - B

        auto abx = (abDistances.data(0))[iContrPair];

        auto aby = (abDistances.data(1))[iContrPair];

        auto abz = (abDistances.data(2))[iContrPair];

        // loop over set of data vectors

        for (size_t i = 0; i < recPattern.size(); i++)
        {
            if ((recPattern[i].first() == 4) && (recPattern[i].second() == 4))
            {
                if (iContrPair == 0) printf("-> applying bra HRR for (44|XX)\n");

                // determine angular momentum of ket side

                auto cang  = recPattern[i].third();

                auto dang  = recPattern[i].fourth();

                auto kcomp = angmom::to_SphericalComponents(cang, dang);

                // get position of integrals in integrals buffer

                auto goff = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {4, 4, cang, dang});

                auto g1off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {3, 5, cang, dang});

                auto g2off = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                         {3, 4, cang, dang});

                // compute contracted integrals

                for (int32_t j = 0; j < kcomp; j++)
                {
                    // set up pointers to (FG|g(r,r')|XX) integrals

                    auto g2_xxx_xxxx = braBuffer.data(g2off + j);

                    auto g2_xxx_xxxy = braBuffer.data(g2off + kcomp + j);

                    auto g2_xxx_xxxz = braBuffer.data(g2off + 2 * kcomp + j);

                    auto g2_xxx_xxyy = braBuffer.data(g2off + 3 * kcomp + j);

                    auto g2_xxx_xxyz = braBuffer.data(g2off + 4 * kcomp + j);

                    auto g2_xxx_xxzz = braBuffer.data(g2off + 5 * kcomp + j);

                    auto g2_xxx_xyyy = braBuffer.data(g2off + 6 * kcomp + j);

                    auto g2_xxx_xyyz = braBuffer.data(g2off + 7 * kcomp + j);

                    auto g2_xxx_xyzz = braBuffer.data(g2off + 8 * kcomp + j);

                    auto g2_xxx_xzzz = braBuffer.data(g2off + 9 * kcomp + j);

                    auto g2_xxx_yyyy = braBuffer.data(g2off + 10 * kcomp + j);

                    auto g2_xxx_yyyz = braBuffer.data(g2off + 11 * kcomp + j);

                    auto g2_xxx_yyzz = braBuffer.data(g2off + 12 * kcomp + j);

                    auto g2_xxx_yzzz = braBuffer.data(g2off + 13 * kcomp + j);

                    auto g2_xxx_zzzz = braBuffer.data(g2off + 14 * kcomp + j);

                    auto g2_xxy_xxxx = braBuffer.data(g2off + 15 * kcomp + j);

                    auto g2_xxy_xxxy = braBuffer.data(g2off + 16 * kcomp + j);

                    auto g2_xxy_xxxz = braBuffer.data(g2off + 17 * kcomp + j);

                    auto g2_xxy_xxyy = braBuffer.data(g2off + 18 * kcomp + j);

                    auto g2_xxy_xxyz = braBuffer.data(g2off + 19 * kcomp + j);

                    auto g2_xxy_xxzz = braBuffer.data(g2off + 20 * kcomp + j);

                    auto g2_xxy_xyyy = braBuffer.data(g2off + 21 * kcomp + j);

                    auto g2_xxy_xyyz = braBuffer.data(g2off + 22 * kcomp + j);

                    auto g2_xxy_xyzz = braBuffer.data(g2off + 23 * kcomp + j);

                    auto g2_xxy_xzzz = braBuffer.data(g2off + 24 * kcomp + j);

                    auto g2_xxy_yyyy = braBuffer.data(g2off + 25 * kcomp + j);

                    auto g2_xxy_yyyz = braBuffer.data(g2off + 26 * kcomp + j);

                    auto g2_xxy_yyzz = braBuffer.data(g2off + 27 * kcomp + j);

                    auto g2_xxy_yzzz = braBuffer.data(g2off + 28 * kcomp + j);

                    auto g2_xxy_zzzz = braBuffer.data(g2off + 29 * kcomp + j);

                    auto g2_xxz_xxxx = braBuffer.data(g2off + 30 * kcomp + j);

                    auto g2_xxz_xxxy = braBuffer.data(g2off + 31 * kcomp + j);

                    auto g2_xxz_xxxz = braBuffer.data(g2off + 32 * kcomp + j);

                    auto g2_xxz_xxyy = braBuffer.data(g2off + 33 * kcomp + j);

                    auto g2_xxz_xxyz = braBuffer.data(g2off + 34 * kcomp + j);

                    auto g2_xxz_xxzz = braBuffer.data(g2off + 35 * kcomp + j);

                    auto g2_xxz_xyyy = braBuffer.data(g2off + 36 * kcomp + j);

                    auto g2_xxz_xyyz = braBuffer.data(g2off + 37 * kcomp + j);

                    auto g2_xxz_xyzz = braBuffer.data(g2off + 38 * kcomp + j);

                    auto g2_xxz_xzzz = braBuffer.data(g2off + 39 * kcomp + j);

                    auto g2_xxz_yyyy = braBuffer.data(g2off + 40 * kcomp + j);

                    auto g2_xxz_yyyz = braBuffer.data(g2off + 41 * kcomp + j);

                    auto g2_xxz_yyzz = braBuffer.data(g2off + 42 * kcomp + j);

                    auto g2_xxz_yzzz = braBuffer.data(g2off + 43 * kcomp + j);

                    auto g2_xxz_zzzz = braBuffer.data(g2off + 44 * kcomp + j);

                    auto g2_xyy_xxxx = braBuffer.data(g2off + 45 * kcomp + j);

                    auto g2_xyy_xxxy = braBuffer.data(g2off + 46 * kcomp + j);

                    auto g2_xyy_xxxz = braBuffer.data(g2off + 47 * kcomp + j);

                    auto g2_xyy_xxyy = braBuffer.data(g2off + 48 * kcomp + j);

                    auto g2_xyy_xxyz = braBuffer.data(g2off + 49 * kcomp + j);

                    auto g2_xyy_xxzz = braBuffer.data(g2off + 50 * kcomp + j);

                    auto g2_xyy_xyyy = braBuffer.data(g2off + 51 * kcomp + j);

                    auto g2_xyy_xyyz = braBuffer.data(g2off + 52 * kcomp + j);

                    auto g2_xyy_xyzz = braBuffer.data(g2off + 53 * kcomp + j);

                    auto g2_xyy_xzzz = braBuffer.data(g2off + 54 * kcomp + j);

                    auto g2_xyy_yyyy = braBuffer.data(g2off + 55 * kcomp + j);

                    auto g2_xyy_yyyz = braBuffer.data(g2off + 56 * kcomp + j);

                    auto g2_xyy_yyzz = braBuffer.data(g2off + 57 * kcomp + j);

                    auto g2_xyy_yzzz = braBuffer.data(g2off + 58 * kcomp + j);

                    auto g2_xyy_zzzz = braBuffer.data(g2off + 59 * kcomp + j);

                    auto g2_xyz_xxxx = braBuffer.data(g2off + 60 * kcomp + j);

                    auto g2_xyz_xxxy = braBuffer.data(g2off + 61 * kcomp + j);

                    auto g2_xyz_xxxz = braBuffer.data(g2off + 62 * kcomp + j);

                    auto g2_xyz_xxyy = braBuffer.data(g2off + 63 * kcomp + j);

                    auto g2_xyz_xxyz = braBuffer.data(g2off + 64 * kcomp + j);

                    auto g2_xyz_xxzz = braBuffer.data(g2off + 65 * kcomp + j);

                    auto g2_xyz_xyyy = braBuffer.data(g2off + 66 * kcomp + j);

                    auto g2_xyz_xyyz = braBuffer.data(g2off + 67 * kcomp + j);

                    auto g2_xyz_xyzz = braBuffer.data(g2off + 68 * kcomp + j);

                    auto g2_xyz_xzzz = braBuffer.data(g2off + 69 * kcomp + j);

                    auto g2_xyz_yyyy = braBuffer.data(g2off + 70 * kcomp + j);

                    auto g2_xyz_yyyz = braBuffer.data(g2off + 71 * kcomp + j);

                    auto g2_xyz_yyzz = braBuffer.data(g2off + 72 * kcomp + j);

                    auto g2_xyz_yzzz = braBuffer.data(g2off + 73 * kcomp + j);

                    auto g2_xyz_zzzz = braBuffer.data(g2off + 74 * kcomp + j);

                    auto g2_xzz_xxxx = braBuffer.data(g2off + 75 * kcomp + j);

                    auto g2_xzz_xxxy = braBuffer.data(g2off + 76 * kcomp + j);

                    auto g2_xzz_xxxz = braBuffer.data(g2off + 77 * kcomp + j);

                    auto g2_xzz_xxyy = braBuffer.data(g2off + 78 * kcomp + j);

                    auto g2_xzz_xxyz = braBuffer.data(g2off + 79 * kcomp + j);

                    auto g2_xzz_xxzz = braBuffer.data(g2off + 80 * kcomp + j);

                    auto g2_xzz_xyyy = braBuffer.data(g2off + 81 * kcomp + j);

                    auto g2_xzz_xyyz = braBuffer.data(g2off + 82 * kcomp + j);

                    auto g2_xzz_xyzz = braBuffer.data(g2off + 83 * kcomp + j);

                    auto g2_xzz_xzzz = braBuffer.data(g2off + 84 * kcomp + j);

                    auto g2_xzz_yyyy = braBuffer.data(g2off + 85 * kcomp + j);

                    auto g2_xzz_yyyz = braBuffer.data(g2off + 86 * kcomp + j);

                    auto g2_xzz_yyzz = braBuffer.data(g2off + 87 * kcomp + j);

                    auto g2_xzz_yzzz = braBuffer.data(g2off + 88 * kcomp + j);

                    auto g2_xzz_zzzz = braBuffer.data(g2off + 89 * kcomp + j);

                    auto g2_yyy_xxxx = braBuffer.data(g2off + 90 * kcomp + j);

                    auto g2_yyy_xxxy = braBuffer.data(g2off + 91 * kcomp + j);

                    auto g2_yyy_xxxz = braBuffer.data(g2off + 92 * kcomp + j);

                    auto g2_yyy_xxyy = braBuffer.data(g2off + 93 * kcomp + j);

                    auto g2_yyy_xxyz = braBuffer.data(g2off + 94 * kcomp + j);

                    auto g2_yyy_xxzz = braBuffer.data(g2off + 95 * kcomp + j);

                    auto g2_yyy_xyyy = braBuffer.data(g2off + 96 * kcomp + j);

                    auto g2_yyy_xyyz = braBuffer.data(g2off + 97 * kcomp + j);

                    auto g2_yyy_xyzz = braBuffer.data(g2off + 98 * kcomp + j);

                    auto g2_yyy_xzzz = braBuffer.data(g2off + 99 * kcomp + j);

                    auto g2_yyy_yyyy = braBuffer.data(g2off + 100 * kcomp + j);

                    auto g2_yyy_yyyz = braBuffer.data(g2off + 101 * kcomp + j);

                    auto g2_yyy_yyzz = braBuffer.data(g2off + 102 * kcomp + j);

                    auto g2_yyy_yzzz = braBuffer.data(g2off + 103 * kcomp + j);

                    auto g2_yyy_zzzz = braBuffer.data(g2off + 104 * kcomp + j);

                    auto g2_yyz_xxxx = braBuffer.data(g2off + 105 * kcomp + j);

                    auto g2_yyz_xxxy = braBuffer.data(g2off + 106 * kcomp + j);

                    auto g2_yyz_xxxz = braBuffer.data(g2off + 107 * kcomp + j);

                    auto g2_yyz_xxyy = braBuffer.data(g2off + 108 * kcomp + j);

                    auto g2_yyz_xxyz = braBuffer.data(g2off + 109 * kcomp + j);

                    auto g2_yyz_xxzz = braBuffer.data(g2off + 110 * kcomp + j);

                    auto g2_yyz_xyyy = braBuffer.data(g2off + 111 * kcomp + j);

                    auto g2_yyz_xyyz = braBuffer.data(g2off + 112 * kcomp + j);

                    auto g2_yyz_xyzz = braBuffer.data(g2off + 113 * kcomp + j);

                    auto g2_yyz_xzzz = braBuffer.data(g2off + 114 * kcomp + j);

                    auto g2_yyz_yyyy = braBuffer.data(g2off + 115 * kcomp + j);

                    auto g2_yyz_yyyz = braBuffer.data(g2off + 116 * kcomp + j);

                    auto g2_yyz_yyzz = braBuffer.data(g2off + 117 * kcomp + j);

                    auto g2_yyz_yzzz = braBuffer.data(g2off + 118 * kcomp + j);

                    auto g2_yyz_zzzz = braBuffer.data(g2off + 119 * kcomp + j);

                    auto g2_yzz_xxxx = braBuffer.data(g2off + 120 * kcomp + j);

                    auto g2_yzz_xxxy = braBuffer.data(g2off + 121 * kcomp + j);

                    auto g2_yzz_xxxz = braBuffer.data(g2off + 122 * kcomp + j);

                    auto g2_yzz_xxyy = braBuffer.data(g2off + 123 * kcomp + j);

                    auto g2_yzz_xxyz = braBuffer.data(g2off + 124 * kcomp + j);

                    auto g2_yzz_xxzz = braBuffer.data(g2off + 125 * kcomp + j);

                    auto g2_yzz_xyyy = braBuffer.data(g2off + 126 * kcomp + j);

                    auto g2_yzz_xyyz = braBuffer.data(g2off + 127 * kcomp + j);

                    auto g2_yzz_xyzz = braBuffer.data(g2off + 128 * kcomp + j);

                    auto g2_yzz_xzzz = braBuffer.data(g2off + 129 * kcomp + j);

                    auto g2_yzz_yyyy = braBuffer.data(g2off + 130 * kcomp + j);

                    auto g2_yzz_yyyz = braBuffer.data(g2off + 131 * kcomp + j);

                    auto g2_yzz_yyzz = braBuffer.data(g2off + 132 * kcomp + j);

                    auto g2_yzz_yzzz = braBuffer.data(g2off + 133 * kcomp + j);

                    auto g2_yzz_zzzz = braBuffer.data(g2off + 134 * kcomp + j);

                    auto g2_zzz_xxxx = braBuffer.data(g2off + 135 * kcomp + j);

                    auto g2_zzz_xxxy = braBuffer.data(g2off + 136 * kcomp + j);

                    auto g2_zzz_xxxz = braBuffer.data(g2off + 137 * kcomp + j);

                    auto g2_zzz_xxyy = braBuffer.data(g2off + 138 * kcomp + j);

                    auto g2_zzz_xxyz = braBuffer.data(g2off + 139 * kcomp + j);

                    auto g2_zzz_xxzz = braBuffer.data(g2off + 140 * kcomp + j);

                    auto g2_zzz_xyyy = braBuffer.data(g2off + 141 * kcomp + j);

                    auto g2_zzz_xyyz = braBuffer.data(g2off + 142 * kcomp + j);

                    auto g2_zzz_xyzz = braBuffer.data(g2off + 143 * kcomp + j);

                    auto g2_zzz_xzzz = braBuffer.data(g2off + 144 * kcomp + j);

                    auto g2_zzz_yyyy = braBuffer.data(g2off + 145 * kcomp + j);

                    auto g2_zzz_yyyz = braBuffer.data(g2off + 146 * kcomp + j);

                    auto g2_zzz_yyzz = braBuffer.data(g2off + 147 * kcomp + j);

                    auto g2_zzz_yzzz = braBuffer.data(g2off + 148 * kcomp + j);

                    auto g2_zzz_zzzz = braBuffer.data(g2off + 149 * kcomp + j);

                    // set up pointers to (FH|g(r,r')|XX) integrals

                    auto g1_xxx_xxxxx = braBuffer.data(g1off + j);

                    auto g1_xxx_xxxxy = braBuffer.data(g1off + kcomp + j);

                    auto g1_xxx_xxxxz = braBuffer.data(g1off + 2 * kcomp + j);

                    auto g1_xxx_xxxyy = braBuffer.data(g1off + 3 * kcomp + j);

                    auto g1_xxx_xxxyz = braBuffer.data(g1off + 4 * kcomp + j);

                    auto g1_xxx_xxxzz = braBuffer.data(g1off + 5 * kcomp + j);

                    auto g1_xxx_xxyyy = braBuffer.data(g1off + 6 * kcomp + j);

                    auto g1_xxx_xxyyz = braBuffer.data(g1off + 7 * kcomp + j);

                    auto g1_xxx_xxyzz = braBuffer.data(g1off + 8 * kcomp + j);

                    auto g1_xxx_xxzzz = braBuffer.data(g1off + 9 * kcomp + j);

                    auto g1_xxx_xyyyy = braBuffer.data(g1off + 10 * kcomp + j);

                    auto g1_xxx_xyyyz = braBuffer.data(g1off + 11 * kcomp + j);

                    auto g1_xxx_xyyzz = braBuffer.data(g1off + 12 * kcomp + j);

                    auto g1_xxx_xyzzz = braBuffer.data(g1off + 13 * kcomp + j);

                    auto g1_xxx_xzzzz = braBuffer.data(g1off + 14 * kcomp + j);

                    auto g1_xxy_xxxxx = braBuffer.data(g1off + 21 * kcomp + j);

                    auto g1_xxy_xxxxy = braBuffer.data(g1off + 22 * kcomp + j);

                    auto g1_xxy_xxxxz = braBuffer.data(g1off + 23 * kcomp + j);

                    auto g1_xxy_xxxyy = braBuffer.data(g1off + 24 * kcomp + j);

                    auto g1_xxy_xxxyz = braBuffer.data(g1off + 25 * kcomp + j);

                    auto g1_xxy_xxxzz = braBuffer.data(g1off + 26 * kcomp + j);

                    auto g1_xxy_xxyyy = braBuffer.data(g1off + 27 * kcomp + j);

                    auto g1_xxy_xxyyz = braBuffer.data(g1off + 28 * kcomp + j);

                    auto g1_xxy_xxyzz = braBuffer.data(g1off + 29 * kcomp + j);

                    auto g1_xxy_xxzzz = braBuffer.data(g1off + 30 * kcomp + j);

                    auto g1_xxy_xyyyy = braBuffer.data(g1off + 31 * kcomp + j);

                    auto g1_xxy_xyyyz = braBuffer.data(g1off + 32 * kcomp + j);

                    auto g1_xxy_xyyzz = braBuffer.data(g1off + 33 * kcomp + j);

                    auto g1_xxy_xyzzz = braBuffer.data(g1off + 34 * kcomp + j);

                    auto g1_xxy_xzzzz = braBuffer.data(g1off + 35 * kcomp + j);

                    auto g1_xxz_xxxxx = braBuffer.data(g1off + 42 * kcomp + j);

                    auto g1_xxz_xxxxy = braBuffer.data(g1off + 43 * kcomp + j);

                    auto g1_xxz_xxxxz = braBuffer.data(g1off + 44 * kcomp + j);

                    auto g1_xxz_xxxyy = braBuffer.data(g1off + 45 * kcomp + j);

                    auto g1_xxz_xxxyz = braBuffer.data(g1off + 46 * kcomp + j);

                    auto g1_xxz_xxxzz = braBuffer.data(g1off + 47 * kcomp + j);

                    auto g1_xxz_xxyyy = braBuffer.data(g1off + 48 * kcomp + j);

                    auto g1_xxz_xxyyz = braBuffer.data(g1off + 49 * kcomp + j);

                    auto g1_xxz_xxyzz = braBuffer.data(g1off + 50 * kcomp + j);

                    auto g1_xxz_xxzzz = braBuffer.data(g1off + 51 * kcomp + j);

                    auto g1_xxz_xyyyy = braBuffer.data(g1off + 52 * kcomp + j);

                    auto g1_xxz_xyyyz = braBuffer.data(g1off + 53 * kcomp + j);

                    auto g1_xxz_xyyzz = braBuffer.data(g1off + 54 * kcomp + j);

                    auto g1_xxz_xyzzz = braBuffer.data(g1off + 55 * kcomp + j);

                    auto g1_xxz_xzzzz = braBuffer.data(g1off + 56 * kcomp + j);

                    auto g1_xyy_xxxxx = braBuffer.data(g1off + 63 * kcomp + j);

                    auto g1_xyy_xxxxy = braBuffer.data(g1off + 64 * kcomp + j);

                    auto g1_xyy_xxxxz = braBuffer.data(g1off + 65 * kcomp + j);

                    auto g1_xyy_xxxyy = braBuffer.data(g1off + 66 * kcomp + j);

                    auto g1_xyy_xxxyz = braBuffer.data(g1off + 67 * kcomp + j);

                    auto g1_xyy_xxxzz = braBuffer.data(g1off + 68 * kcomp + j);

                    auto g1_xyy_xxyyy = braBuffer.data(g1off + 69 * kcomp + j);

                    auto g1_xyy_xxyyz = braBuffer.data(g1off + 70 * kcomp + j);

                    auto g1_xyy_xxyzz = braBuffer.data(g1off + 71 * kcomp + j);

                    auto g1_xyy_xxzzz = braBuffer.data(g1off + 72 * kcomp + j);

                    auto g1_xyy_xyyyy = braBuffer.data(g1off + 73 * kcomp + j);

                    auto g1_xyy_xyyyz = braBuffer.data(g1off + 74 * kcomp + j);

                    auto g1_xyy_xyyzz = braBuffer.data(g1off + 75 * kcomp + j);

                    auto g1_xyy_xyzzz = braBuffer.data(g1off + 76 * kcomp + j);

                    auto g1_xyy_xzzzz = braBuffer.data(g1off + 77 * kcomp + j);

                    auto g1_xyz_xxxxx = braBuffer.data(g1off + 84 * kcomp + j);

                    auto g1_xyz_xxxxy = braBuffer.data(g1off + 85 * kcomp + j);

                    auto g1_xyz_xxxxz = braBuffer.data(g1off + 86 * kcomp + j);

                    auto g1_xyz_xxxyy = braBuffer.data(g1off + 87 * kcomp + j);

                    auto g1_xyz_xxxyz = braBuffer.data(g1off + 88 * kcomp + j);

                    auto g1_xyz_xxxzz = braBuffer.data(g1off + 89 * kcomp + j);

                    auto g1_xyz_xxyyy = braBuffer.data(g1off + 90 * kcomp + j);

                    auto g1_xyz_xxyyz = braBuffer.data(g1off + 91 * kcomp + j);

                    auto g1_xyz_xxyzz = braBuffer.data(g1off + 92 * kcomp + j);

                    auto g1_xyz_xxzzz = braBuffer.data(g1off + 93 * kcomp + j);

                    auto g1_xyz_xyyyy = braBuffer.data(g1off + 94 * kcomp + j);

                    auto g1_xyz_xyyyz = braBuffer.data(g1off + 95 * kcomp + j);

                    auto g1_xyz_xyyzz = braBuffer.data(g1off + 96 * kcomp + j);

                    auto g1_xyz_xyzzz = braBuffer.data(g1off + 97 * kcomp + j);

                    auto g1_xyz_xzzzz = braBuffer.data(g1off + 98 * kcomp + j);

                    auto g1_xzz_xxxxx = braBuffer.data(g1off + 105 * kcomp + j);

                    auto g1_xzz_xxxxy = braBuffer.data(g1off + 106 * kcomp + j);

                    auto g1_xzz_xxxxz = braBuffer.data(g1off + 107 * kcomp + j);

                    auto g1_xzz_xxxyy = braBuffer.data(g1off + 108 * kcomp + j);

                    auto g1_xzz_xxxyz = braBuffer.data(g1off + 109 * kcomp + j);

                    auto g1_xzz_xxxzz = braBuffer.data(g1off + 110 * kcomp + j);

                    auto g1_xzz_xxyyy = braBuffer.data(g1off + 111 * kcomp + j);

                    auto g1_xzz_xxyyz = braBuffer.data(g1off + 112 * kcomp + j);

                    auto g1_xzz_xxyzz = braBuffer.data(g1off + 113 * kcomp + j);

                    auto g1_xzz_xxzzz = braBuffer.data(g1off + 114 * kcomp + j);

                    auto g1_xzz_xyyyy = braBuffer.data(g1off + 115 * kcomp + j);

                    auto g1_xzz_xyyyz = braBuffer.data(g1off + 116 * kcomp + j);

                    auto g1_xzz_xyyzz = braBuffer.data(g1off + 117 * kcomp + j);

                    auto g1_xzz_xyzzz = braBuffer.data(g1off + 118 * kcomp + j);

                    auto g1_xzz_xzzzz = braBuffer.data(g1off + 119 * kcomp + j);

                    auto g1_yyy_xxxxx = braBuffer.data(g1off + 126 * kcomp + j);

                    auto g1_yyy_xxxxy = braBuffer.data(g1off + 127 * kcomp + j);

                    auto g1_yyy_xxxxz = braBuffer.data(g1off + 128 * kcomp + j);

                    auto g1_yyy_xxxyy = braBuffer.data(g1off + 129 * kcomp + j);

                    auto g1_yyy_xxxyz = braBuffer.data(g1off + 130 * kcomp + j);

                    auto g1_yyy_xxxzz = braBuffer.data(g1off + 131 * kcomp + j);

                    auto g1_yyy_xxyyy = braBuffer.data(g1off + 132 * kcomp + j);

                    auto g1_yyy_xxyyz = braBuffer.data(g1off + 133 * kcomp + j);

                    auto g1_yyy_xxyzz = braBuffer.data(g1off + 134 * kcomp + j);

                    auto g1_yyy_xxzzz = braBuffer.data(g1off + 135 * kcomp + j);

                    auto g1_yyy_xyyyy = braBuffer.data(g1off + 136 * kcomp + j);

                    auto g1_yyy_xyyyz = braBuffer.data(g1off + 137 * kcomp + j);

                    auto g1_yyy_xyyzz = braBuffer.data(g1off + 138 * kcomp + j);

                    auto g1_yyy_xyzzz = braBuffer.data(g1off + 139 * kcomp + j);

                    auto g1_yyy_xzzzz = braBuffer.data(g1off + 140 * kcomp + j);

                    auto g1_yyy_yyyyy = braBuffer.data(g1off + 141 * kcomp + j);

                    auto g1_yyy_yyyyz = braBuffer.data(g1off + 142 * kcomp + j);

                    auto g1_yyy_yyyzz = braBuffer.data(g1off + 143 * kcomp + j);

                    auto g1_yyy_yyzzz = braBuffer.data(g1off + 144 * kcomp + j);

                    auto g1_yyy_yzzzz = braBuffer.data(g1off + 145 * kcomp + j);

                    auto g1_yyz_xxxxx = braBuffer.data(g1off + 147 * kcomp + j);

                    auto g1_yyz_xxxxy = braBuffer.data(g1off + 148 * kcomp + j);

                    auto g1_yyz_xxxxz = braBuffer.data(g1off + 149 * kcomp + j);

                    auto g1_yyz_xxxyy = braBuffer.data(g1off + 150 * kcomp + j);

                    auto g1_yyz_xxxyz = braBuffer.data(g1off + 151 * kcomp + j);

                    auto g1_yyz_xxxzz = braBuffer.data(g1off + 152 * kcomp + j);

                    auto g1_yyz_xxyyy = braBuffer.data(g1off + 153 * kcomp + j);

                    auto g1_yyz_xxyyz = braBuffer.data(g1off + 154 * kcomp + j);

                    auto g1_yyz_xxyzz = braBuffer.data(g1off + 155 * kcomp + j);

                    auto g1_yyz_xxzzz = braBuffer.data(g1off + 156 * kcomp + j);

                    auto g1_yyz_xyyyy = braBuffer.data(g1off + 157 * kcomp + j);

                    auto g1_yyz_xyyyz = braBuffer.data(g1off + 158 * kcomp + j);

                    auto g1_yyz_xyyzz = braBuffer.data(g1off + 159 * kcomp + j);

                    auto g1_yyz_xyzzz = braBuffer.data(g1off + 160 * kcomp + j);

                    auto g1_yyz_xzzzz = braBuffer.data(g1off + 161 * kcomp + j);

                    auto g1_yyz_yyyyy = braBuffer.data(g1off + 162 * kcomp + j);

                    auto g1_yyz_yyyyz = braBuffer.data(g1off + 163 * kcomp + j);

                    auto g1_yyz_yyyzz = braBuffer.data(g1off + 164 * kcomp + j);

                    auto g1_yyz_yyzzz = braBuffer.data(g1off + 165 * kcomp + j);

                    auto g1_yyz_yzzzz = braBuffer.data(g1off + 166 * kcomp + j);

                    auto g1_yzz_xxxxx = braBuffer.data(g1off + 168 * kcomp + j);

                    auto g1_yzz_xxxxy = braBuffer.data(g1off + 169 * kcomp + j);

                    auto g1_yzz_xxxxz = braBuffer.data(g1off + 170 * kcomp + j);

                    auto g1_yzz_xxxyy = braBuffer.data(g1off + 171 * kcomp + j);

                    auto g1_yzz_xxxyz = braBuffer.data(g1off + 172 * kcomp + j);

                    auto g1_yzz_xxxzz = braBuffer.data(g1off + 173 * kcomp + j);

                    auto g1_yzz_xxyyy = braBuffer.data(g1off + 174 * kcomp + j);

                    auto g1_yzz_xxyyz = braBuffer.data(g1off + 175 * kcomp + j);

                    auto g1_yzz_xxyzz = braBuffer.data(g1off + 176 * kcomp + j);

                    auto g1_yzz_xxzzz = braBuffer.data(g1off + 177 * kcomp + j);

                    auto g1_yzz_xyyyy = braBuffer.data(g1off + 178 * kcomp + j);

                    auto g1_yzz_xyyyz = braBuffer.data(g1off + 179 * kcomp + j);

                    auto g1_yzz_xyyzz = braBuffer.data(g1off + 180 * kcomp + j);

                    auto g1_yzz_xyzzz = braBuffer.data(g1off + 181 * kcomp + j);

                    auto g1_yzz_xzzzz = braBuffer.data(g1off + 182 * kcomp + j);

                    auto g1_yzz_yyyyy = braBuffer.data(g1off + 183 * kcomp + j);

                    auto g1_yzz_yyyyz = braBuffer.data(g1off + 184 * kcomp + j);

                    auto g1_yzz_yyyzz = braBuffer.data(g1off + 185 * kcomp + j);

                    auto g1_yzz_yyzzz = braBuffer.data(g1off + 186 * kcomp + j);

                    auto g1_yzz_yzzzz = braBuffer.data(g1off + 187 * kcomp + j);

                    auto g1_zzz_xxxxx = braBuffer.data(g1off + 189 * kcomp + j);

                    auto g1_zzz_xxxxy = braBuffer.data(g1off + 190 * kcomp + j);

                    auto g1_zzz_xxxxz = braBuffer.data(g1off + 191 * kcomp + j);

                    auto g1_zzz_xxxyy = braBuffer.data(g1off + 192 * kcomp + j);

                    auto g1_zzz_xxxyz = braBuffer.data(g1off + 193 * kcomp + j);

                    auto g1_zzz_xxxzz = braBuffer.data(g1off + 194 * kcomp + j);

                    auto g1_zzz_xxyyy = braBuffer.data(g1off + 195 * kcomp + j);

                    auto g1_zzz_xxyyz = braBuffer.data(g1off + 196 * kcomp + j);

                    auto g1_zzz_xxyzz = braBuffer.data(g1off + 197 * kcomp + j);

                    auto g1_zzz_xxzzz = braBuffer.data(g1off + 198 * kcomp + j);

                    auto g1_zzz_xyyyy = braBuffer.data(g1off + 199 * kcomp + j);

                    auto g1_zzz_xyyyz = braBuffer.data(g1off + 200 * kcomp + j);

                    auto g1_zzz_xyyzz = braBuffer.data(g1off + 201 * kcomp + j);

                    auto g1_zzz_xyzzz = braBuffer.data(g1off + 202 * kcomp + j);

                    auto g1_zzz_xzzzz = braBuffer.data(g1off + 203 * kcomp + j);

                    auto g1_zzz_yyyyy = braBuffer.data(g1off + 204 * kcomp + j);

                    auto g1_zzz_yyyyz = braBuffer.data(g1off + 205 * kcomp + j);

                    auto g1_zzz_yyyzz = braBuffer.data(g1off + 206 * kcomp + j);

                    auto g1_zzz_yyzzz = braBuffer.data(g1off + 207 * kcomp + j);

                    auto g1_zzz_yzzzz = braBuffer.data(g1off + 208 * kcomp + j);

                    auto g1_zzz_zzzzz = braBuffer.data(g1off + 209 * kcomp + j);

                    // set up pointers to (GG|g(r,r')|XX) integrals

                    auto g_xxxx_xxxx = braBuffer.data(goff + j);

                    auto g_xxxx_xxxy = braBuffer.data(goff + kcomp + j);

                    auto g_xxxx_xxxz = braBuffer.data(goff + 2 * kcomp + j);

                    auto g_xxxx_xxyy = braBuffer.data(goff + 3 * kcomp + j);

                    auto g_xxxx_xxyz = braBuffer.data(goff + 4 * kcomp + j);

                    auto g_xxxx_xxzz = braBuffer.data(goff + 5 * kcomp + j);

                    auto g_xxxx_xyyy = braBuffer.data(goff + 6 * kcomp + j);

                    auto g_xxxx_xyyz = braBuffer.data(goff + 7 * kcomp + j);

                    auto g_xxxx_xyzz = braBuffer.data(goff + 8 * kcomp + j);

                    auto g_xxxx_xzzz = braBuffer.data(goff + 9 * kcomp + j);

                    auto g_xxxx_yyyy = braBuffer.data(goff + 10 * kcomp + j);

                    auto g_xxxx_yyyz = braBuffer.data(goff + 11 * kcomp + j);

                    auto g_xxxx_yyzz = braBuffer.data(goff + 12 * kcomp + j);

                    auto g_xxxx_yzzz = braBuffer.data(goff + 13 * kcomp + j);

                    auto g_xxxx_zzzz = braBuffer.data(goff + 14 * kcomp + j);

                    auto g_xxxy_xxxx = braBuffer.data(goff + 15 * kcomp + j);

                    auto g_xxxy_xxxy = braBuffer.data(goff + 16 * kcomp + j);

                    auto g_xxxy_xxxz = braBuffer.data(goff + 17 * kcomp + j);

                    auto g_xxxy_xxyy = braBuffer.data(goff + 18 * kcomp + j);

                    auto g_xxxy_xxyz = braBuffer.data(goff + 19 * kcomp + j);

                    auto g_xxxy_xxzz = braBuffer.data(goff + 20 * kcomp + j);

                    auto g_xxxy_xyyy = braBuffer.data(goff + 21 * kcomp + j);

                    auto g_xxxy_xyyz = braBuffer.data(goff + 22 * kcomp + j);

                    auto g_xxxy_xyzz = braBuffer.data(goff + 23 * kcomp + j);

                    auto g_xxxy_xzzz = braBuffer.data(goff + 24 * kcomp + j);

                    auto g_xxxy_yyyy = braBuffer.data(goff + 25 * kcomp + j);

                    auto g_xxxy_yyyz = braBuffer.data(goff + 26 * kcomp + j);

                    auto g_xxxy_yyzz = braBuffer.data(goff + 27 * kcomp + j);

                    auto g_xxxy_yzzz = braBuffer.data(goff + 28 * kcomp + j);

                    auto g_xxxy_zzzz = braBuffer.data(goff + 29 * kcomp + j);

                    auto g_xxxz_xxxx = braBuffer.data(goff + 30 * kcomp + j);

                    auto g_xxxz_xxxy = braBuffer.data(goff + 31 * kcomp + j);

                    auto g_xxxz_xxxz = braBuffer.data(goff + 32 * kcomp + j);

                    auto g_xxxz_xxyy = braBuffer.data(goff + 33 * kcomp + j);

                    auto g_xxxz_xxyz = braBuffer.data(goff + 34 * kcomp + j);

                    auto g_xxxz_xxzz = braBuffer.data(goff + 35 * kcomp + j);

                    auto g_xxxz_xyyy = braBuffer.data(goff + 36 * kcomp + j);

                    auto g_xxxz_xyyz = braBuffer.data(goff + 37 * kcomp + j);

                    auto g_xxxz_xyzz = braBuffer.data(goff + 38 * kcomp + j);

                    auto g_xxxz_xzzz = braBuffer.data(goff + 39 * kcomp + j);

                    auto g_xxxz_yyyy = braBuffer.data(goff + 40 * kcomp + j);

                    auto g_xxxz_yyyz = braBuffer.data(goff + 41 * kcomp + j);

                    auto g_xxxz_yyzz = braBuffer.data(goff + 42 * kcomp + j);

                    auto g_xxxz_yzzz = braBuffer.data(goff + 43 * kcomp + j);

                    auto g_xxxz_zzzz = braBuffer.data(goff + 44 * kcomp + j);

                    auto g_xxyy_xxxx = braBuffer.data(goff + 45 * kcomp + j);

                    auto g_xxyy_xxxy = braBuffer.data(goff + 46 * kcomp + j);

                    auto g_xxyy_xxxz = braBuffer.data(goff + 47 * kcomp + j);

                    auto g_xxyy_xxyy = braBuffer.data(goff + 48 * kcomp + j);

                    auto g_xxyy_xxyz = braBuffer.data(goff + 49 * kcomp + j);

                    auto g_xxyy_xxzz = braBuffer.data(goff + 50 * kcomp + j);

                    auto g_xxyy_xyyy = braBuffer.data(goff + 51 * kcomp + j);

                    auto g_xxyy_xyyz = braBuffer.data(goff + 52 * kcomp + j);

                    auto g_xxyy_xyzz = braBuffer.data(goff + 53 * kcomp + j);

                    auto g_xxyy_xzzz = braBuffer.data(goff + 54 * kcomp + j);

                    auto g_xxyy_yyyy = braBuffer.data(goff + 55 * kcomp + j);

                    auto g_xxyy_yyyz = braBuffer.data(goff + 56 * kcomp + j);

                    auto g_xxyy_yyzz = braBuffer.data(goff + 57 * kcomp + j);

                    auto g_xxyy_yzzz = braBuffer.data(goff + 58 * kcomp + j);

                    auto g_xxyy_zzzz = braBuffer.data(goff + 59 * kcomp + j);

                    auto g_xxyz_xxxx = braBuffer.data(goff + 60 * kcomp + j);

                    auto g_xxyz_xxxy = braBuffer.data(goff + 61 * kcomp + j);

                    auto g_xxyz_xxxz = braBuffer.data(goff + 62 * kcomp + j);

                    auto g_xxyz_xxyy = braBuffer.data(goff + 63 * kcomp + j);

                    auto g_xxyz_xxyz = braBuffer.data(goff + 64 * kcomp + j);

                    auto g_xxyz_xxzz = braBuffer.data(goff + 65 * kcomp + j);

                    auto g_xxyz_xyyy = braBuffer.data(goff + 66 * kcomp + j);

                    auto g_xxyz_xyyz = braBuffer.data(goff + 67 * kcomp + j);

                    auto g_xxyz_xyzz = braBuffer.data(goff + 68 * kcomp + j);

                    auto g_xxyz_xzzz = braBuffer.data(goff + 69 * kcomp + j);

                    auto g_xxyz_yyyy = braBuffer.data(goff + 70 * kcomp + j);

                    auto g_xxyz_yyyz = braBuffer.data(goff + 71 * kcomp + j);

                    auto g_xxyz_yyzz = braBuffer.data(goff + 72 * kcomp + j);

                    auto g_xxyz_yzzz = braBuffer.data(goff + 73 * kcomp + j);

                    auto g_xxyz_zzzz = braBuffer.data(goff + 74 * kcomp + j);

                    auto g_xxzz_xxxx = braBuffer.data(goff + 75 * kcomp + j);

                    auto g_xxzz_xxxy = braBuffer.data(goff + 76 * kcomp + j);

                    auto g_xxzz_xxxz = braBuffer.data(goff + 77 * kcomp + j);

                    auto g_xxzz_xxyy = braBuffer.data(goff + 78 * kcomp + j);

                    auto g_xxzz_xxyz = braBuffer.data(goff + 79 * kcomp + j);

                    auto g_xxzz_xxzz = braBuffer.data(goff + 80 * kcomp + j);

                    auto g_xxzz_xyyy = braBuffer.data(goff + 81 * kcomp + j);

                    auto g_xxzz_xyyz = braBuffer.data(goff + 82 * kcomp + j);

                    auto g_xxzz_xyzz = braBuffer.data(goff + 83 * kcomp + j);

                    auto g_xxzz_xzzz = braBuffer.data(goff + 84 * kcomp + j);

                    auto g_xxzz_yyyy = braBuffer.data(goff + 85 * kcomp + j);

                    auto g_xxzz_yyyz = braBuffer.data(goff + 86 * kcomp + j);

                    auto g_xxzz_yyzz = braBuffer.data(goff + 87 * kcomp + j);

                    auto g_xxzz_yzzz = braBuffer.data(goff + 88 * kcomp + j);

                    auto g_xxzz_zzzz = braBuffer.data(goff + 89 * kcomp + j);

                    auto g_xyyy_xxxx = braBuffer.data(goff + 90 * kcomp + j);

                    auto g_xyyy_xxxy = braBuffer.data(goff + 91 * kcomp + j);

                    auto g_xyyy_xxxz = braBuffer.data(goff + 92 * kcomp + j);

                    auto g_xyyy_xxyy = braBuffer.data(goff + 93 * kcomp + j);

                    auto g_xyyy_xxyz = braBuffer.data(goff + 94 * kcomp + j);

                    auto g_xyyy_xxzz = braBuffer.data(goff + 95 * kcomp + j);

                    auto g_xyyy_xyyy = braBuffer.data(goff + 96 * kcomp + j);

                    auto g_xyyy_xyyz = braBuffer.data(goff + 97 * kcomp + j);

                    auto g_xyyy_xyzz = braBuffer.data(goff + 98 * kcomp + j);

                    auto g_xyyy_xzzz = braBuffer.data(goff + 99 * kcomp + j);

                    auto g_xyyy_yyyy = braBuffer.data(goff + 100 * kcomp + j);

                    auto g_xyyy_yyyz = braBuffer.data(goff + 101 * kcomp + j);

                    auto g_xyyy_yyzz = braBuffer.data(goff + 102 * kcomp + j);

                    auto g_xyyy_yzzz = braBuffer.data(goff + 103 * kcomp + j);

                    auto g_xyyy_zzzz = braBuffer.data(goff + 104 * kcomp + j);

                    auto g_xyyz_xxxx = braBuffer.data(goff + 105 * kcomp + j);

                    auto g_xyyz_xxxy = braBuffer.data(goff + 106 * kcomp + j);

                    auto g_xyyz_xxxz = braBuffer.data(goff + 107 * kcomp + j);

                    auto g_xyyz_xxyy = braBuffer.data(goff + 108 * kcomp + j);

                    auto g_xyyz_xxyz = braBuffer.data(goff + 109 * kcomp + j);

                    auto g_xyyz_xxzz = braBuffer.data(goff + 110 * kcomp + j);

                    auto g_xyyz_xyyy = braBuffer.data(goff + 111 * kcomp + j);

                    auto g_xyyz_xyyz = braBuffer.data(goff + 112 * kcomp + j);

                    auto g_xyyz_xyzz = braBuffer.data(goff + 113 * kcomp + j);

                    auto g_xyyz_xzzz = braBuffer.data(goff + 114 * kcomp + j);

                    auto g_xyyz_yyyy = braBuffer.data(goff + 115 * kcomp + j);

                    auto g_xyyz_yyyz = braBuffer.data(goff + 116 * kcomp + j);

                    auto g_xyyz_yyzz = braBuffer.data(goff + 117 * kcomp + j);

                    auto g_xyyz_yzzz = braBuffer.data(goff + 118 * kcomp + j);

                    auto g_xyyz_zzzz = braBuffer.data(goff + 119 * kcomp + j);

                    auto g_xyzz_xxxx = braBuffer.data(goff + 120 * kcomp + j);

                    auto g_xyzz_xxxy = braBuffer.data(goff + 121 * kcomp + j);

                    auto g_xyzz_xxxz = braBuffer.data(goff + 122 * kcomp + j);

                    auto g_xyzz_xxyy = braBuffer.data(goff + 123 * kcomp + j);

                    auto g_xyzz_xxyz = braBuffer.data(goff + 124 * kcomp + j);

                    auto g_xyzz_xxzz = braBuffer.data(goff + 125 * kcomp + j);

                    auto g_xyzz_xyyy = braBuffer.data(goff + 126 * kcomp + j);

                    auto g_xyzz_xyyz = braBuffer.data(goff + 127 * kcomp + j);

                    auto g_xyzz_xyzz = braBuffer.data(goff + 128 * kcomp + j);

                    auto g_xyzz_xzzz = braBuffer.data(goff + 129 * kcomp + j);

                    auto g_xyzz_yyyy = braBuffer.data(goff + 130 * kcomp + j);

                    auto g_xyzz_yyyz = braBuffer.data(goff + 131 * kcomp + j);

                    auto g_xyzz_yyzz = braBuffer.data(goff + 132 * kcomp + j);

                    auto g_xyzz_yzzz = braBuffer.data(goff + 133 * kcomp + j);

                    auto g_xyzz_zzzz = braBuffer.data(goff + 134 * kcomp + j);

                    auto g_xzzz_xxxx = braBuffer.data(goff + 135 * kcomp + j);

                    auto g_xzzz_xxxy = braBuffer.data(goff + 136 * kcomp + j);

                    auto g_xzzz_xxxz = braBuffer.data(goff + 137 * kcomp + j);

                    auto g_xzzz_xxyy = braBuffer.data(goff + 138 * kcomp + j);

                    auto g_xzzz_xxyz = braBuffer.data(goff + 139 * kcomp + j);

                    auto g_xzzz_xxzz = braBuffer.data(goff + 140 * kcomp + j);

                    auto g_xzzz_xyyy = braBuffer.data(goff + 141 * kcomp + j);

                    auto g_xzzz_xyyz = braBuffer.data(goff + 142 * kcomp + j);

                    auto g_xzzz_xyzz = braBuffer.data(goff + 143 * kcomp + j);

                    auto g_xzzz_xzzz = braBuffer.data(goff + 144 * kcomp + j);

                    auto g_xzzz_yyyy = braBuffer.data(goff + 145 * kcomp + j);

                    auto g_xzzz_yyyz = braBuffer.data(goff + 146 * kcomp + j);

                    auto g_xzzz_yyzz = braBuffer.data(goff + 147 * kcomp + j);

                    auto g_xzzz_yzzz = braBuffer.data(goff + 148 * kcomp + j);

                    auto g_xzzz_zzzz = braBuffer.data(goff + 149 * kcomp + j);

                    auto g_yyyy_xxxx = braBuffer.data(goff + 150 * kcomp + j);

                    auto g_yyyy_xxxy = braBuffer.data(goff + 151 * kcomp + j);

                    auto g_yyyy_xxxz = braBuffer.data(goff + 152 * kcomp + j);

                    auto g_yyyy_xxyy = braBuffer.data(goff + 153 * kcomp + j);

                    auto g_yyyy_xxyz = braBuffer.data(goff + 154 * kcomp + j);

                    auto g_yyyy_xxzz = braBuffer.data(goff + 155 * kcomp + j);

                    auto g_yyyy_xyyy = braBuffer.data(goff + 156 * kcomp + j);

                    auto g_yyyy_xyyz = braBuffer.data(goff + 157 * kcomp + j);

                    auto g_yyyy_xyzz = braBuffer.data(goff + 158 * kcomp + j);

                    auto g_yyyy_xzzz = braBuffer.data(goff + 159 * kcomp + j);

                    auto g_yyyy_yyyy = braBuffer.data(goff + 160 * kcomp + j);

                    auto g_yyyy_yyyz = braBuffer.data(goff + 161 * kcomp + j);

                    auto g_yyyy_yyzz = braBuffer.data(goff + 162 * kcomp + j);

                    auto g_yyyy_yzzz = braBuffer.data(goff + 163 * kcomp + j);

                    auto g_yyyy_zzzz = braBuffer.data(goff + 164 * kcomp + j);

                    auto g_yyyz_xxxx = braBuffer.data(goff + 165 * kcomp + j);

                    auto g_yyyz_xxxy = braBuffer.data(goff + 166 * kcomp + j);

                    auto g_yyyz_xxxz = braBuffer.data(goff + 167 * kcomp + j);

                    auto g_yyyz_xxyy = braBuffer.data(goff + 168 * kcomp + j);

                    auto g_yyyz_xxyz = braBuffer.data(goff + 169 * kcomp + j);

                    auto g_yyyz_xxzz = braBuffer.data(goff + 170 * kcomp + j);

                    auto g_yyyz_xyyy = braBuffer.data(goff + 171 * kcomp + j);

                    auto g_yyyz_xyyz = braBuffer.data(goff + 172 * kcomp + j);

                    auto g_yyyz_xyzz = braBuffer.data(goff + 173 * kcomp + j);

                    auto g_yyyz_xzzz = braBuffer.data(goff + 174 * kcomp + j);

                    auto g_yyyz_yyyy = braBuffer.data(goff + 175 * kcomp + j);

                    auto g_yyyz_yyyz = braBuffer.data(goff + 176 * kcomp + j);

                    auto g_yyyz_yyzz = braBuffer.data(goff + 177 * kcomp + j);

                    auto g_yyyz_yzzz = braBuffer.data(goff + 178 * kcomp + j);

                    auto g_yyyz_zzzz = braBuffer.data(goff + 179 * kcomp + j);

                    auto g_yyzz_xxxx = braBuffer.data(goff + 180 * kcomp + j);

                    auto g_yyzz_xxxy = braBuffer.data(goff + 181 * kcomp + j);

                    auto g_yyzz_xxxz = braBuffer.data(goff + 182 * kcomp + j);

                    auto g_yyzz_xxyy = braBuffer.data(goff + 183 * kcomp + j);

                    auto g_yyzz_xxyz = braBuffer.data(goff + 184 * kcomp + j);

                    auto g_yyzz_xxzz = braBuffer.data(goff + 185 * kcomp + j);

                    auto g_yyzz_xyyy = braBuffer.data(goff + 186 * kcomp + j);

                    auto g_yyzz_xyyz = braBuffer.data(goff + 187 * kcomp + j);

                    auto g_yyzz_xyzz = braBuffer.data(goff + 188 * kcomp + j);

                    auto g_yyzz_xzzz = braBuffer.data(goff + 189 * kcomp + j);

                    auto g_yyzz_yyyy = braBuffer.data(goff + 190 * kcomp + j);

                    auto g_yyzz_yyyz = braBuffer.data(goff + 191 * kcomp + j);

                    auto g_yyzz_yyzz = braBuffer.data(goff + 192 * kcomp + j);

                    auto g_yyzz_yzzz = braBuffer.data(goff + 193 * kcomp + j);

                    auto g_yyzz_zzzz = braBuffer.data(goff + 194 * kcomp + j);

                    auto g_yzzz_xxxx = braBuffer.data(goff + 195 * kcomp + j);

                    auto g_yzzz_xxxy = braBuffer.data(goff + 196 * kcomp + j);

                    auto g_yzzz_xxxz = braBuffer.data(goff + 197 * kcomp + j);

                    auto g_yzzz_xxyy = braBuffer.data(goff + 198 * kcomp + j);

                    auto g_yzzz_xxyz = braBuffer.data(goff + 199 * kcomp + j);

                    auto g_yzzz_xxzz = braBuffer.data(goff + 200 * kcomp + j);

                    auto g_yzzz_xyyy = braBuffer.data(goff + 201 * kcomp + j);

                    auto g_yzzz_xyyz = braBuffer.data(goff + 202 * kcomp + j);

                    auto g_yzzz_xyzz = braBuffer.data(goff + 203 * kcomp + j);

                    auto g_yzzz_xzzz = braBuffer.data(goff + 204 * kcomp + j);

                    auto g_yzzz_yyyy = braBuffer.data(goff + 205 * kcomp + j);

                    auto g_yzzz_yyyz = braBuffer.data(goff + 206 * kcomp + j);

                    auto g_yzzz_yyzz = braBuffer.data(goff + 207 * kcomp + j);

                    auto g_yzzz_yzzz = braBuffer.data(goff + 208 * kcomp + j);

                    auto g_yzzz_zzzz = braBuffer.data(goff + 209 * kcomp + j);

                    auto g_zzzz_xxxx = braBuffer.data(goff + 210 * kcomp + j);

                    auto g_zzzz_xxxy = braBuffer.data(goff + 211 * kcomp + j);

                    auto g_zzzz_xxxz = braBuffer.data(goff + 212 * kcomp + j);

                    auto g_zzzz_xxyy = braBuffer.data(goff + 213 * kcomp + j);

                    auto g_zzzz_xxyz = braBuffer.data(goff + 214 * kcomp + j);

                    auto g_zzzz_xxzz = braBuffer.data(goff + 215 * kcomp + j);

                    auto g_zzzz_xyyy = braBuffer.data(goff + 216 * kcomp + j);

                    auto g_zzzz_xyyz = braBuffer.data(goff + 217 * kcomp + j);

                    auto g_zzzz_xyzz = braBuffer.data(goff + 218 * kcomp + j);

                    auto g_zzzz_xzzz = braBuffer.data(goff + 219 * kcomp + j);

                    auto g_zzzz_yyyy = braBuffer.data(goff + 220 * kcomp + j);

                    auto g_zzzz_yyyz = braBuffer.data(goff + 221 * kcomp + j);

                    auto g_zzzz_yyzz = braBuffer.data(goff + 222 * kcomp + j);

                    auto g_zzzz_yzzz = braBuffer.data(goff + 223 * kcomp + j);

                    auto g_zzzz_zzzz = braBuffer.data(goff + 224 * kcomp + j);

                    #pragma omp simd aligned(g2_xxx_xxxx, g2_xxx_xxxy, g2_xxx_xxxz,\
                                             g2_xxx_xxyy, g2_xxx_xxyz, g2_xxx_xxzz,\
                                             g2_xxx_xyyy, g2_xxx_xyyz, g2_xxx_xyzz,\
                                             g2_xxx_xzzz, g2_xxx_yyyy, g2_xxx_yyyz,\
                                             g2_xxx_yyzz, g2_xxx_yzzz, g2_xxx_zzzz,\
                                             g2_xxy_xxxx, g2_xxy_xxxy, g2_xxy_xxxz,\
                                             g2_xxy_xxyy, g2_xxy_xxyz, g2_xxy_xxzz,\
                                             g2_xxy_xyyy, g2_xxy_xyyz, g2_xxy_xyzz,\
                                             g2_xxy_xzzz, g2_xxy_yyyy, g2_xxy_yyyz,\
                                             g2_xxy_yyzz, g2_xxy_yzzz, g2_xxy_zzzz,\
                                             g2_xxz_xxxx, g2_xxz_xxxy, g2_xxz_xxxz,\
                                             g2_xxz_xxyy, g2_xxz_xxyz, g2_xxz_xxzz,\
                                             g2_xxz_xyyy, g2_xxz_xyyz, g2_xxz_xyzz,\
                                             g2_xxz_xzzz, g2_xxz_yyyy, g2_xxz_yyyz,\
                                             g2_xxz_yyzz, g2_xxz_yzzz, g2_xxz_zzzz,\
                                             g2_xyy_xxxx, g2_xyy_xxxy, g2_xyy_xxxz,\
                                             g2_xyy_xxyy, g2_xyy_xxyz, g2_xyy_xxzz,\
                                             g2_xyy_xyyy, g2_xyy_xyyz, g2_xyy_xyzz,\
                                             g2_xyy_xzzz, g2_xyy_yyyy, g2_xyy_yyyz,\
                                             g2_xyy_yyzz, g2_xyy_yzzz, g2_xyy_zzzz,\
                                             g2_xyz_xxxx, g2_xyz_xxxy, g2_xyz_xxxz,\
                                             g2_xyz_xxyy, g2_xyz_xxyz, g2_xyz_xxzz,\
                                             g2_xyz_xyyy, g2_xyz_xyyz, g2_xyz_xyzz,\
                                             g2_xyz_xzzz, g2_xyz_yyyy, g2_xyz_yyyz,\
                                             g2_xyz_yyzz, g2_xyz_yzzz, g2_xyz_zzzz,\
                                             g2_xzz_xxxx, g2_xzz_xxxy, g2_xzz_xxxz,\
                                             g2_xzz_xxyy, g2_xzz_xxyz, g2_xzz_xxzz,\
                                             g2_xzz_xyyy, g2_xzz_xyyz, g2_xzz_xyzz,\
                                             g2_xzz_xzzz, g2_xzz_yyyy, g2_xzz_yyyz,\
                                             g2_xzz_yyzz, g2_xzz_yzzz, g2_xzz_zzzz,\
                                             g2_yyy_xxxx, g2_yyy_xxxy, g2_yyy_xxxz,\
                                             g2_yyy_xxyy, g2_yyy_xxyz, g2_yyy_xxzz,\
                                             g2_yyy_xyyy, g2_yyy_xyyz, g2_yyy_xyzz,\
                                             g2_yyy_xzzz, g2_yyy_yyyy, g2_yyy_yyyz,\
                                             g2_yyy_yyzz, g2_yyy_yzzz, g2_yyy_zzzz,\
                                             g2_yyz_xxxx, g2_yyz_xxxy, g2_yyz_xxxz,\
                                             g2_yyz_xxyy, g2_yyz_xxyz, g2_yyz_xxzz,\
                                             g2_yyz_xyyy, g2_yyz_xyyz, g2_yyz_xyzz,\
                                             g2_yyz_xzzz, g2_yyz_yyyy, g2_yyz_yyyz,\
                                             g2_yyz_yyzz, g2_yyz_yzzz, g2_yyz_zzzz,\
                                             g2_yzz_xxxx, g2_yzz_xxxy, g2_yzz_xxxz,\
                                             g2_yzz_xxyy, g2_yzz_xxyz, g2_yzz_xxzz,\
                                             g2_yzz_xyyy, g2_yzz_xyyz, g2_yzz_xyzz,\
                                             g2_yzz_xzzz, g2_yzz_yyyy, g2_yzz_yyyz,\
                                             g2_yzz_yyzz, g2_yzz_yzzz, g2_yzz_zzzz,\
                                             g2_zzz_xxxx, g2_zzz_xxxy, g2_zzz_xxxz,\
                                             g2_zzz_xxyy, g2_zzz_xxyz, g2_zzz_xxzz,\
                                             g2_zzz_xyyy, g2_zzz_xyyz, g2_zzz_xyzz,\
                                             g2_zzz_xzzz, g2_zzz_yyyy, g2_zzz_yyyz,\
                                             g2_zzz_yyzz, g2_zzz_yzzz, g2_zzz_zzzz,\
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
                                             g1_yyy_yyzzz, g1_yyy_yzzzz, \
                                             g1_yyz_xxxxx, g1_yyz_xxxxy, g1_yyz_xxxxz,\
                                             g1_yyz_xxxyy, g1_yyz_xxxyz, g1_yyz_xxxzz,\
                                             g1_yyz_xxyyy, g1_yyz_xxyyz, g1_yyz_xxyzz,\
                                             g1_yyz_xxzzz, g1_yyz_xyyyy, g1_yyz_xyyyz,\
                                             g1_yyz_xyyzz, g1_yyz_xyzzz, g1_yyz_xzzzz,\
                                             g1_yyz_yyyyy, g1_yyz_yyyyz, g1_yyz_yyyzz,\
                                             g1_yyz_yyzzz, g1_yyz_yzzzz, \
                                             g1_yzz_xxxxx, g1_yzz_xxxxy, g1_yzz_xxxxz,\
                                             g1_yzz_xxxyy, g1_yzz_xxxyz, g1_yzz_xxxzz,\
                                             g1_yzz_xxyyy, g1_yzz_xxyyz, g1_yzz_xxyzz,\
                                             g1_yzz_xxzzz, g1_yzz_xyyyy, g1_yzz_xyyyz,\
                                             g1_yzz_xyyzz, g1_yzz_xyzzz, g1_yzz_xzzzz,\
                                             g1_yzz_yyyyy, g1_yzz_yyyyz, g1_yzz_yyyzz,\
                                             g1_yzz_yyzzz, g1_yzz_yzzzz, \
                                             g1_zzz_xxxxx, g1_zzz_xxxxy, g1_zzz_xxxxz,\
                                             g1_zzz_xxxyy, g1_zzz_xxxyz, g1_zzz_xxxzz,\
                                             g1_zzz_xxyyy, g1_zzz_xxyyz, g1_zzz_xxyzz,\
                                             g1_zzz_xxzzz, g1_zzz_xyyyy, g1_zzz_xyyyz,\
                                             g1_zzz_xyyzz, g1_zzz_xyzzz, g1_zzz_xzzzz,\
                                             g1_zzz_yyyyy, g1_zzz_yyyyz, g1_zzz_yyyzz,\
                                             g1_zzz_yyzzz, g1_zzz_yzzzz, g1_zzz_zzzzz,\
                                             g_xxxx_xxxx, g_xxxx_xxxy, g_xxxx_xxxz,\
                                             g_xxxx_xxyy, g_xxxx_xxyz, g_xxxx_xxzz,\
                                             g_xxxx_xyyy, g_xxxx_xyyz, g_xxxx_xyzz,\
                                             g_xxxx_xzzz, g_xxxx_yyyy, g_xxxx_yyyz,\
                                             g_xxxx_yyzz, g_xxxx_yzzz, g_xxxx_zzzz,\
                                             g_xxxy_xxxx, g_xxxy_xxxy, g_xxxy_xxxz,\
                                             g_xxxy_xxyy, g_xxxy_xxyz, g_xxxy_xxzz,\
                                             g_xxxy_xyyy, g_xxxy_xyyz, g_xxxy_xyzz,\
                                             g_xxxy_xzzz, g_xxxy_yyyy, g_xxxy_yyyz,\
                                             g_xxxy_yyzz, g_xxxy_yzzz, g_xxxy_zzzz,\
                                             g_xxxz_xxxx, g_xxxz_xxxy, g_xxxz_xxxz,\
                                             g_xxxz_xxyy, g_xxxz_xxyz, g_xxxz_xxzz,\
                                             g_xxxz_xyyy, g_xxxz_xyyz, g_xxxz_xyzz,\
                                             g_xxxz_xzzz, g_xxxz_yyyy, g_xxxz_yyyz,\
                                             g_xxxz_yyzz, g_xxxz_yzzz, g_xxxz_zzzz,\
                                             g_xxyy_xxxx, g_xxyy_xxxy, g_xxyy_xxxz,\
                                             g_xxyy_xxyy, g_xxyy_xxyz, g_xxyy_xxzz,\
                                             g_xxyy_xyyy, g_xxyy_xyyz, g_xxyy_xyzz,\
                                             g_xxyy_xzzz, g_xxyy_yyyy, g_xxyy_yyyz,\
                                             g_xxyy_yyzz, g_xxyy_yzzz, g_xxyy_zzzz,\
                                             g_xxyz_xxxx, g_xxyz_xxxy, g_xxyz_xxxz,\
                                             g_xxyz_xxyy, g_xxyz_xxyz, g_xxyz_xxzz,\
                                             g_xxyz_xyyy, g_xxyz_xyyz, g_xxyz_xyzz,\
                                             g_xxyz_xzzz, g_xxyz_yyyy, g_xxyz_yyyz,\
                                             g_xxyz_yyzz, g_xxyz_yzzz, g_xxyz_zzzz,\
                                             g_xxzz_xxxx, g_xxzz_xxxy, g_xxzz_xxxz,\
                                             g_xxzz_xxyy, g_xxzz_xxyz, g_xxzz_xxzz,\
                                             g_xxzz_xyyy, g_xxzz_xyyz, g_xxzz_xyzz,\
                                             g_xxzz_xzzz, g_xxzz_yyyy, g_xxzz_yyyz,\
                                             g_xxzz_yyzz, g_xxzz_yzzz, g_xxzz_zzzz,\
                                             g_xyyy_xxxx, g_xyyy_xxxy, g_xyyy_xxxz,\
                                             g_xyyy_xxyy, g_xyyy_xxyz, g_xyyy_xxzz,\
                                             g_xyyy_xyyy, g_xyyy_xyyz, g_xyyy_xyzz,\
                                             g_xyyy_xzzz, g_xyyy_yyyy, g_xyyy_yyyz,\
                                             g_xyyy_yyzz, g_xyyy_yzzz, g_xyyy_zzzz,\
                                             g_xyyz_xxxx, g_xyyz_xxxy, g_xyyz_xxxz,\
                                             g_xyyz_xxyy, g_xyyz_xxyz, g_xyyz_xxzz,\
                                             g_xyyz_xyyy, g_xyyz_xyyz, g_xyyz_xyzz,\
                                             g_xyyz_xzzz, g_xyyz_yyyy, g_xyyz_yyyz,\
                                             g_xyyz_yyzz, g_xyyz_yzzz, g_xyyz_zzzz,\
                                             g_xyzz_xxxx, g_xyzz_xxxy, g_xyzz_xxxz,\
                                             g_xyzz_xxyy, g_xyzz_xxyz, g_xyzz_xxzz,\
                                             g_xyzz_xyyy, g_xyzz_xyyz, g_xyzz_xyzz,\
                                             g_xyzz_xzzz, g_xyzz_yyyy, g_xyzz_yyyz,\
                                             g_xyzz_yyzz, g_xyzz_yzzz, g_xyzz_zzzz,\
                                             g_xzzz_xxxx, g_xzzz_xxxy, g_xzzz_xxxz,\
                                             g_xzzz_xxyy, g_xzzz_xxyz, g_xzzz_xxzz,\
                                             g_xzzz_xyyy, g_xzzz_xyyz, g_xzzz_xyzz,\
                                             g_xzzz_xzzz, g_xzzz_yyyy, g_xzzz_yyyz,\
                                             g_xzzz_yyzz, g_xzzz_yzzz, g_xzzz_zzzz,\
                                             g_yyyy_xxxx, g_yyyy_xxxy, g_yyyy_xxxz,\
                                             g_yyyy_xxyy, g_yyyy_xxyz, g_yyyy_xxzz,\
                                             g_yyyy_xyyy, g_yyyy_xyyz, g_yyyy_xyzz,\
                                             g_yyyy_xzzz, g_yyyy_yyyy, g_yyyy_yyyz,\
                                             g_yyyy_yyzz, g_yyyy_yzzz, g_yyyy_zzzz,\
                                             g_yyyz_xxxx, g_yyyz_xxxy, g_yyyz_xxxz,\
                                             g_yyyz_xxyy, g_yyyz_xxyz, g_yyyz_xxzz,\
                                             g_yyyz_xyyy, g_yyyz_xyyz, g_yyyz_xyzz,\
                                             g_yyyz_xzzz, g_yyyz_yyyy, g_yyyz_yyyz,\
                                             g_yyyz_yyzz, g_yyyz_yzzz, g_yyyz_zzzz,\
                                             g_yyzz_xxxx, g_yyzz_xxxy, g_yyzz_xxxz,\
                                             g_yyzz_xxyy, g_yyzz_xxyz, g_yyzz_xxzz,\
                                             g_yyzz_xyyy, g_yyzz_xyyz, g_yyzz_xyzz,\
                                             g_yyzz_xzzz, g_yyzz_yyyy, g_yyzz_yyyz,\
                                             g_yyzz_yyzz, g_yyzz_yzzz, g_yyzz_zzzz,\
                                             g_yzzz_xxxx, g_yzzz_xxxy, g_yzzz_xxxz,\
                                             g_yzzz_xxyy, g_yzzz_xxyz, g_yzzz_xxzz,\
                                             g_yzzz_xyyy, g_yzzz_xyyz, g_yzzz_xyzz,\
                                             g_yzzz_xzzz, g_yzzz_yyyy, g_yzzz_yyyz,\
                                             g_yzzz_yyzz, g_yzzz_yzzz, g_yzzz_zzzz,\
                                             g_zzzz_xxxx, g_zzzz_xxxy, g_zzzz_xxxz,\
                                             g_zzzz_xxyy, g_zzzz_xxyz, g_zzzz_xxzz,\
                                             g_zzzz_xyyy, g_zzzz_xyyz, g_zzzz_xyzz,\
                                             g_zzzz_xzzz, g_zzzz_yyyy, g_zzzz_yyyz,\
                                             g_zzzz_yyzz, g_zzzz_yzzz, g_zzzz_zzzz: VLX_ALIGN)
                    for (int32_t k = 0; k < kdim; k++)
                    {
                        // leading x component

                        g_xxxx_xxxx[k] = g1_xxx_xxxxx[k] - abx * g2_xxx_xxxx[k];

                        g_xxxx_xxxy[k] = g1_xxx_xxxxy[k] - abx * g2_xxx_xxxy[k];

                        g_xxxx_xxxz[k] = g1_xxx_xxxxz[k] - abx * g2_xxx_xxxz[k];

                        g_xxxx_xxyy[k] = g1_xxx_xxxyy[k] - abx * g2_xxx_xxyy[k];

                        g_xxxx_xxyz[k] = g1_xxx_xxxyz[k] - abx * g2_xxx_xxyz[k];

                        g_xxxx_xxzz[k] = g1_xxx_xxxzz[k] - abx * g2_xxx_xxzz[k];

                        g_xxxx_xyyy[k] = g1_xxx_xxyyy[k] - abx * g2_xxx_xyyy[k];

                        g_xxxx_xyyz[k] = g1_xxx_xxyyz[k] - abx * g2_xxx_xyyz[k];

                        g_xxxx_xyzz[k] = g1_xxx_xxyzz[k] - abx * g2_xxx_xyzz[k];

                        g_xxxx_xzzz[k] = g1_xxx_xxzzz[k] - abx * g2_xxx_xzzz[k];

                        g_xxxx_yyyy[k] = g1_xxx_xyyyy[k] - abx * g2_xxx_yyyy[k];

                        g_xxxx_yyyz[k] = g1_xxx_xyyyz[k] - abx * g2_xxx_yyyz[k];

                        g_xxxx_yyzz[k] = g1_xxx_xyyzz[k] - abx * g2_xxx_yyzz[k];

                        g_xxxx_yzzz[k] = g1_xxx_xyzzz[k] - abx * g2_xxx_yzzz[k];

                        g_xxxx_zzzz[k] = g1_xxx_xzzzz[k] - abx * g2_xxx_zzzz[k];

                        g_xxxy_xxxx[k] = g1_xxy_xxxxx[k] - abx * g2_xxy_xxxx[k];

                        g_xxxy_xxxy[k] = g1_xxy_xxxxy[k] - abx * g2_xxy_xxxy[k];

                        g_xxxy_xxxz[k] = g1_xxy_xxxxz[k] - abx * g2_xxy_xxxz[k];

                        g_xxxy_xxyy[k] = g1_xxy_xxxyy[k] - abx * g2_xxy_xxyy[k];

                        g_xxxy_xxyz[k] = g1_xxy_xxxyz[k] - abx * g2_xxy_xxyz[k];

                        g_xxxy_xxzz[k] = g1_xxy_xxxzz[k] - abx * g2_xxy_xxzz[k];

                        g_xxxy_xyyy[k] = g1_xxy_xxyyy[k] - abx * g2_xxy_xyyy[k];

                        g_xxxy_xyyz[k] = g1_xxy_xxyyz[k] - abx * g2_xxy_xyyz[k];

                        g_xxxy_xyzz[k] = g1_xxy_xxyzz[k] - abx * g2_xxy_xyzz[k];

                        g_xxxy_xzzz[k] = g1_xxy_xxzzz[k] - abx * g2_xxy_xzzz[k];

                        g_xxxy_yyyy[k] = g1_xxy_xyyyy[k] - abx * g2_xxy_yyyy[k];

                        g_xxxy_yyyz[k] = g1_xxy_xyyyz[k] - abx * g2_xxy_yyyz[k];

                        g_xxxy_yyzz[k] = g1_xxy_xyyzz[k] - abx * g2_xxy_yyzz[k];

                        g_xxxy_yzzz[k] = g1_xxy_xyzzz[k] - abx * g2_xxy_yzzz[k];

                        g_xxxy_zzzz[k] = g1_xxy_xzzzz[k] - abx * g2_xxy_zzzz[k];

                        g_xxxz_xxxx[k] = g1_xxz_xxxxx[k] - abx * g2_xxz_xxxx[k];

                        g_xxxz_xxxy[k] = g1_xxz_xxxxy[k] - abx * g2_xxz_xxxy[k];

                        g_xxxz_xxxz[k] = g1_xxz_xxxxz[k] - abx * g2_xxz_xxxz[k];

                        g_xxxz_xxyy[k] = g1_xxz_xxxyy[k] - abx * g2_xxz_xxyy[k];

                        g_xxxz_xxyz[k] = g1_xxz_xxxyz[k] - abx * g2_xxz_xxyz[k];

                        g_xxxz_xxzz[k] = g1_xxz_xxxzz[k] - abx * g2_xxz_xxzz[k];

                        g_xxxz_xyyy[k] = g1_xxz_xxyyy[k] - abx * g2_xxz_xyyy[k];

                        g_xxxz_xyyz[k] = g1_xxz_xxyyz[k] - abx * g2_xxz_xyyz[k];

                        g_xxxz_xyzz[k] = g1_xxz_xxyzz[k] - abx * g2_xxz_xyzz[k];

                        g_xxxz_xzzz[k] = g1_xxz_xxzzz[k] - abx * g2_xxz_xzzz[k];

                        g_xxxz_yyyy[k] = g1_xxz_xyyyy[k] - abx * g2_xxz_yyyy[k];

                        g_xxxz_yyyz[k] = g1_xxz_xyyyz[k] - abx * g2_xxz_yyyz[k];

                        g_xxxz_yyzz[k] = g1_xxz_xyyzz[k] - abx * g2_xxz_yyzz[k];

                        g_xxxz_yzzz[k] = g1_xxz_xyzzz[k] - abx * g2_xxz_yzzz[k];

                        g_xxxz_zzzz[k] = g1_xxz_xzzzz[k] - abx * g2_xxz_zzzz[k];

                        g_xxyy_xxxx[k] = g1_xyy_xxxxx[k] - abx * g2_xyy_xxxx[k];

                        g_xxyy_xxxy[k] = g1_xyy_xxxxy[k] - abx * g2_xyy_xxxy[k];

                        g_xxyy_xxxz[k] = g1_xyy_xxxxz[k] - abx * g2_xyy_xxxz[k];

                        g_xxyy_xxyy[k] = g1_xyy_xxxyy[k] - abx * g2_xyy_xxyy[k];

                        g_xxyy_xxyz[k] = g1_xyy_xxxyz[k] - abx * g2_xyy_xxyz[k];

                        g_xxyy_xxzz[k] = g1_xyy_xxxzz[k] - abx * g2_xyy_xxzz[k];

                        g_xxyy_xyyy[k] = g1_xyy_xxyyy[k] - abx * g2_xyy_xyyy[k];

                        g_xxyy_xyyz[k] = g1_xyy_xxyyz[k] - abx * g2_xyy_xyyz[k];

                        g_xxyy_xyzz[k] = g1_xyy_xxyzz[k] - abx * g2_xyy_xyzz[k];

                        g_xxyy_xzzz[k] = g1_xyy_xxzzz[k] - abx * g2_xyy_xzzz[k];

                        g_xxyy_yyyy[k] = g1_xyy_xyyyy[k] - abx * g2_xyy_yyyy[k];

                        g_xxyy_yyyz[k] = g1_xyy_xyyyz[k] - abx * g2_xyy_yyyz[k];

                        g_xxyy_yyzz[k] = g1_xyy_xyyzz[k] - abx * g2_xyy_yyzz[k];

                        g_xxyy_yzzz[k] = g1_xyy_xyzzz[k] - abx * g2_xyy_yzzz[k];

                        g_xxyy_zzzz[k] = g1_xyy_xzzzz[k] - abx * g2_xyy_zzzz[k];

                        g_xxyz_xxxx[k] = g1_xyz_xxxxx[k] - abx * g2_xyz_xxxx[k];

                        g_xxyz_xxxy[k] = g1_xyz_xxxxy[k] - abx * g2_xyz_xxxy[k];

                        g_xxyz_xxxz[k] = g1_xyz_xxxxz[k] - abx * g2_xyz_xxxz[k];

                        g_xxyz_xxyy[k] = g1_xyz_xxxyy[k] - abx * g2_xyz_xxyy[k];

                        g_xxyz_xxyz[k] = g1_xyz_xxxyz[k] - abx * g2_xyz_xxyz[k];

                        g_xxyz_xxzz[k] = g1_xyz_xxxzz[k] - abx * g2_xyz_xxzz[k];

                        g_xxyz_xyyy[k] = g1_xyz_xxyyy[k] - abx * g2_xyz_xyyy[k];

                        g_xxyz_xyyz[k] = g1_xyz_xxyyz[k] - abx * g2_xyz_xyyz[k];

                        g_xxyz_xyzz[k] = g1_xyz_xxyzz[k] - abx * g2_xyz_xyzz[k];

                        g_xxyz_xzzz[k] = g1_xyz_xxzzz[k] - abx * g2_xyz_xzzz[k];

                        g_xxyz_yyyy[k] = g1_xyz_xyyyy[k] - abx * g2_xyz_yyyy[k];

                        g_xxyz_yyyz[k] = g1_xyz_xyyyz[k] - abx * g2_xyz_yyyz[k];

                        g_xxyz_yyzz[k] = g1_xyz_xyyzz[k] - abx * g2_xyz_yyzz[k];

                        g_xxyz_yzzz[k] = g1_xyz_xyzzz[k] - abx * g2_xyz_yzzz[k];

                        g_xxyz_zzzz[k] = g1_xyz_xzzzz[k] - abx * g2_xyz_zzzz[k];

                        g_xxzz_xxxx[k] = g1_xzz_xxxxx[k] - abx * g2_xzz_xxxx[k];

                        g_xxzz_xxxy[k] = g1_xzz_xxxxy[k] - abx * g2_xzz_xxxy[k];

                        g_xxzz_xxxz[k] = g1_xzz_xxxxz[k] - abx * g2_xzz_xxxz[k];

                        g_xxzz_xxyy[k] = g1_xzz_xxxyy[k] - abx * g2_xzz_xxyy[k];

                        g_xxzz_xxyz[k] = g1_xzz_xxxyz[k] - abx * g2_xzz_xxyz[k];

                        g_xxzz_xxzz[k] = g1_xzz_xxxzz[k] - abx * g2_xzz_xxzz[k];

                        g_xxzz_xyyy[k] = g1_xzz_xxyyy[k] - abx * g2_xzz_xyyy[k];

                        g_xxzz_xyyz[k] = g1_xzz_xxyyz[k] - abx * g2_xzz_xyyz[k];

                        g_xxzz_xyzz[k] = g1_xzz_xxyzz[k] - abx * g2_xzz_xyzz[k];

                        g_xxzz_xzzz[k] = g1_xzz_xxzzz[k] - abx * g2_xzz_xzzz[k];

                        g_xxzz_yyyy[k] = g1_xzz_xyyyy[k] - abx * g2_xzz_yyyy[k];

                        g_xxzz_yyyz[k] = g1_xzz_xyyyz[k] - abx * g2_xzz_yyyz[k];

                        g_xxzz_yyzz[k] = g1_xzz_xyyzz[k] - abx * g2_xzz_yyzz[k];

                        g_xxzz_yzzz[k] = g1_xzz_xyzzz[k] - abx * g2_xzz_yzzz[k];

                        g_xxzz_zzzz[k] = g1_xzz_xzzzz[k] - abx * g2_xzz_zzzz[k];

                        g_xyyy_xxxx[k] = g1_yyy_xxxxx[k] - abx * g2_yyy_xxxx[k];

                        g_xyyy_xxxy[k] = g1_yyy_xxxxy[k] - abx * g2_yyy_xxxy[k];

                        g_xyyy_xxxz[k] = g1_yyy_xxxxz[k] - abx * g2_yyy_xxxz[k];

                        g_xyyy_xxyy[k] = g1_yyy_xxxyy[k] - abx * g2_yyy_xxyy[k];

                        g_xyyy_xxyz[k] = g1_yyy_xxxyz[k] - abx * g2_yyy_xxyz[k];

                        g_xyyy_xxzz[k] = g1_yyy_xxxzz[k] - abx * g2_yyy_xxzz[k];

                        g_xyyy_xyyy[k] = g1_yyy_xxyyy[k] - abx * g2_yyy_xyyy[k];

                        g_xyyy_xyyz[k] = g1_yyy_xxyyz[k] - abx * g2_yyy_xyyz[k];

                        g_xyyy_xyzz[k] = g1_yyy_xxyzz[k] - abx * g2_yyy_xyzz[k];

                        g_xyyy_xzzz[k] = g1_yyy_xxzzz[k] - abx * g2_yyy_xzzz[k];

                        g_xyyy_yyyy[k] = g1_yyy_xyyyy[k] - abx * g2_yyy_yyyy[k];

                        g_xyyy_yyyz[k] = g1_yyy_xyyyz[k] - abx * g2_yyy_yyyz[k];

                        g_xyyy_yyzz[k] = g1_yyy_xyyzz[k] - abx * g2_yyy_yyzz[k];

                        g_xyyy_yzzz[k] = g1_yyy_xyzzz[k] - abx * g2_yyy_yzzz[k];

                        g_xyyy_zzzz[k] = g1_yyy_xzzzz[k] - abx * g2_yyy_zzzz[k];

                        g_xyyz_xxxx[k] = g1_yyz_xxxxx[k] - abx * g2_yyz_xxxx[k];

                        g_xyyz_xxxy[k] = g1_yyz_xxxxy[k] - abx * g2_yyz_xxxy[k];

                        g_xyyz_xxxz[k] = g1_yyz_xxxxz[k] - abx * g2_yyz_xxxz[k];

                        g_xyyz_xxyy[k] = g1_yyz_xxxyy[k] - abx * g2_yyz_xxyy[k];

                        g_xyyz_xxyz[k] = g1_yyz_xxxyz[k] - abx * g2_yyz_xxyz[k];

                        g_xyyz_xxzz[k] = g1_yyz_xxxzz[k] - abx * g2_yyz_xxzz[k];

                        g_xyyz_xyyy[k] = g1_yyz_xxyyy[k] - abx * g2_yyz_xyyy[k];

                        g_xyyz_xyyz[k] = g1_yyz_xxyyz[k] - abx * g2_yyz_xyyz[k];

                        g_xyyz_xyzz[k] = g1_yyz_xxyzz[k] - abx * g2_yyz_xyzz[k];

                        g_xyyz_xzzz[k] = g1_yyz_xxzzz[k] - abx * g2_yyz_xzzz[k];

                        g_xyyz_yyyy[k] = g1_yyz_xyyyy[k] - abx * g2_yyz_yyyy[k];

                        g_xyyz_yyyz[k] = g1_yyz_xyyyz[k] - abx * g2_yyz_yyyz[k];

                        g_xyyz_yyzz[k] = g1_yyz_xyyzz[k] - abx * g2_yyz_yyzz[k];

                        g_xyyz_yzzz[k] = g1_yyz_xyzzz[k] - abx * g2_yyz_yzzz[k];

                        g_xyyz_zzzz[k] = g1_yyz_xzzzz[k] - abx * g2_yyz_zzzz[k];

                        g_xyzz_xxxx[k] = g1_yzz_xxxxx[k] - abx * g2_yzz_xxxx[k];

                        g_xyzz_xxxy[k] = g1_yzz_xxxxy[k] - abx * g2_yzz_xxxy[k];

                        g_xyzz_xxxz[k] = g1_yzz_xxxxz[k] - abx * g2_yzz_xxxz[k];

                        g_xyzz_xxyy[k] = g1_yzz_xxxyy[k] - abx * g2_yzz_xxyy[k];

                        g_xyzz_xxyz[k] = g1_yzz_xxxyz[k] - abx * g2_yzz_xxyz[k];

                        g_xyzz_xxzz[k] = g1_yzz_xxxzz[k] - abx * g2_yzz_xxzz[k];

                        g_xyzz_xyyy[k] = g1_yzz_xxyyy[k] - abx * g2_yzz_xyyy[k];

                        g_xyzz_xyyz[k] = g1_yzz_xxyyz[k] - abx * g2_yzz_xyyz[k];

                        g_xyzz_xyzz[k] = g1_yzz_xxyzz[k] - abx * g2_yzz_xyzz[k];

                        g_xyzz_xzzz[k] = g1_yzz_xxzzz[k] - abx * g2_yzz_xzzz[k];

                        g_xyzz_yyyy[k] = g1_yzz_xyyyy[k] - abx * g2_yzz_yyyy[k];

                        g_xyzz_yyyz[k] = g1_yzz_xyyyz[k] - abx * g2_yzz_yyyz[k];

                        g_xyzz_yyzz[k] = g1_yzz_xyyzz[k] - abx * g2_yzz_yyzz[k];

                        g_xyzz_yzzz[k] = g1_yzz_xyzzz[k] - abx * g2_yzz_yzzz[k];

                        g_xyzz_zzzz[k] = g1_yzz_xzzzz[k] - abx * g2_yzz_zzzz[k];

                        g_xzzz_xxxx[k] = g1_zzz_xxxxx[k] - abx * g2_zzz_xxxx[k];

                        g_xzzz_xxxy[k] = g1_zzz_xxxxy[k] - abx * g2_zzz_xxxy[k];

                        g_xzzz_xxxz[k] = g1_zzz_xxxxz[k] - abx * g2_zzz_xxxz[k];

                        g_xzzz_xxyy[k] = g1_zzz_xxxyy[k] - abx * g2_zzz_xxyy[k];

                        g_xzzz_xxyz[k] = g1_zzz_xxxyz[k] - abx * g2_zzz_xxyz[k];

                        g_xzzz_xxzz[k] = g1_zzz_xxxzz[k] - abx * g2_zzz_xxzz[k];

                        g_xzzz_xyyy[k] = g1_zzz_xxyyy[k] - abx * g2_zzz_xyyy[k];

                        g_xzzz_xyyz[k] = g1_zzz_xxyyz[k] - abx * g2_zzz_xyyz[k];

                        g_xzzz_xyzz[k] = g1_zzz_xxyzz[k] - abx * g2_zzz_xyzz[k];

                        g_xzzz_xzzz[k] = g1_zzz_xxzzz[k] - abx * g2_zzz_xzzz[k];

                        g_xzzz_yyyy[k] = g1_zzz_xyyyy[k] - abx * g2_zzz_yyyy[k];

                        g_xzzz_yyyz[k] = g1_zzz_xyyyz[k] - abx * g2_zzz_yyyz[k];

                        g_xzzz_yyzz[k] = g1_zzz_xyyzz[k] - abx * g2_zzz_yyzz[k];

                        g_xzzz_yzzz[k] = g1_zzz_xyzzz[k] - abx * g2_zzz_yzzz[k];

                        g_xzzz_zzzz[k] = g1_zzz_xzzzz[k] - abx * g2_zzz_zzzz[k];

                        // leading y component

                        g_yyyy_xxxx[k] = g1_yyy_xxxxy[k] - aby * g2_yyy_xxxx[k];

                        g_yyyy_xxxy[k] = g1_yyy_xxxyy[k] - aby * g2_yyy_xxxy[k];

                        g_yyyy_xxxz[k] = g1_yyy_xxxyz[k] - aby * g2_yyy_xxxz[k];

                        g_yyyy_xxyy[k] = g1_yyy_xxyyy[k] - aby * g2_yyy_xxyy[k];

                        g_yyyy_xxyz[k] = g1_yyy_xxyyz[k] - aby * g2_yyy_xxyz[k];

                        g_yyyy_xxzz[k] = g1_yyy_xxyzz[k] - aby * g2_yyy_xxzz[k];

                        g_yyyy_xyyy[k] = g1_yyy_xyyyy[k] - aby * g2_yyy_xyyy[k];

                        g_yyyy_xyyz[k] = g1_yyy_xyyyz[k] - aby * g2_yyy_xyyz[k];

                        g_yyyy_xyzz[k] = g1_yyy_xyyzz[k] - aby * g2_yyy_xyzz[k];

                        g_yyyy_xzzz[k] = g1_yyy_xyzzz[k] - aby * g2_yyy_xzzz[k];

                        g_yyyy_yyyy[k] = g1_yyy_yyyyy[k] - aby * g2_yyy_yyyy[k];

                        g_yyyy_yyyz[k] = g1_yyy_yyyyz[k] - aby * g2_yyy_yyyz[k];

                        g_yyyy_yyzz[k] = g1_yyy_yyyzz[k] - aby * g2_yyy_yyzz[k];

                        g_yyyy_yzzz[k] = g1_yyy_yyzzz[k] - aby * g2_yyy_yzzz[k];

                        g_yyyy_zzzz[k] = g1_yyy_yzzzz[k] - aby * g2_yyy_zzzz[k];

                        g_yyyz_xxxx[k] = g1_yyz_xxxxy[k] - aby * g2_yyz_xxxx[k];

                        g_yyyz_xxxy[k] = g1_yyz_xxxyy[k] - aby * g2_yyz_xxxy[k];

                        g_yyyz_xxxz[k] = g1_yyz_xxxyz[k] - aby * g2_yyz_xxxz[k];

                        g_yyyz_xxyy[k] = g1_yyz_xxyyy[k] - aby * g2_yyz_xxyy[k];

                        g_yyyz_xxyz[k] = g1_yyz_xxyyz[k] - aby * g2_yyz_xxyz[k];

                        g_yyyz_xxzz[k] = g1_yyz_xxyzz[k] - aby * g2_yyz_xxzz[k];

                        g_yyyz_xyyy[k] = g1_yyz_xyyyy[k] - aby * g2_yyz_xyyy[k];

                        g_yyyz_xyyz[k] = g1_yyz_xyyyz[k] - aby * g2_yyz_xyyz[k];

                        g_yyyz_xyzz[k] = g1_yyz_xyyzz[k] - aby * g2_yyz_xyzz[k];

                        g_yyyz_xzzz[k] = g1_yyz_xyzzz[k] - aby * g2_yyz_xzzz[k];

                        g_yyyz_yyyy[k] = g1_yyz_yyyyy[k] - aby * g2_yyz_yyyy[k];

                        g_yyyz_yyyz[k] = g1_yyz_yyyyz[k] - aby * g2_yyz_yyyz[k];

                        g_yyyz_yyzz[k] = g1_yyz_yyyzz[k] - aby * g2_yyz_yyzz[k];

                        g_yyyz_yzzz[k] = g1_yyz_yyzzz[k] - aby * g2_yyz_yzzz[k];

                        g_yyyz_zzzz[k] = g1_yyz_yzzzz[k] - aby * g2_yyz_zzzz[k];

                        g_yyzz_xxxx[k] = g1_yzz_xxxxy[k] - aby * g2_yzz_xxxx[k];

                        g_yyzz_xxxy[k] = g1_yzz_xxxyy[k] - aby * g2_yzz_xxxy[k];

                        g_yyzz_xxxz[k] = g1_yzz_xxxyz[k] - aby * g2_yzz_xxxz[k];

                        g_yyzz_xxyy[k] = g1_yzz_xxyyy[k] - aby * g2_yzz_xxyy[k];

                        g_yyzz_xxyz[k] = g1_yzz_xxyyz[k] - aby * g2_yzz_xxyz[k];

                        g_yyzz_xxzz[k] = g1_yzz_xxyzz[k] - aby * g2_yzz_xxzz[k];

                        g_yyzz_xyyy[k] = g1_yzz_xyyyy[k] - aby * g2_yzz_xyyy[k];

                        g_yyzz_xyyz[k] = g1_yzz_xyyyz[k] - aby * g2_yzz_xyyz[k];

                        g_yyzz_xyzz[k] = g1_yzz_xyyzz[k] - aby * g2_yzz_xyzz[k];

                        g_yyzz_xzzz[k] = g1_yzz_xyzzz[k] - aby * g2_yzz_xzzz[k];

                        g_yyzz_yyyy[k] = g1_yzz_yyyyy[k] - aby * g2_yzz_yyyy[k];

                        g_yyzz_yyyz[k] = g1_yzz_yyyyz[k] - aby * g2_yzz_yyyz[k];

                        g_yyzz_yyzz[k] = g1_yzz_yyyzz[k] - aby * g2_yzz_yyzz[k];

                        g_yyzz_yzzz[k] = g1_yzz_yyzzz[k] - aby * g2_yzz_yzzz[k];

                        g_yyzz_zzzz[k] = g1_yzz_yzzzz[k] - aby * g2_yzz_zzzz[k];

                        g_yzzz_xxxx[k] = g1_zzz_xxxxy[k] - aby * g2_zzz_xxxx[k];

                        g_yzzz_xxxy[k] = g1_zzz_xxxyy[k] - aby * g2_zzz_xxxy[k];

                        g_yzzz_xxxz[k] = g1_zzz_xxxyz[k] - aby * g2_zzz_xxxz[k];

                        g_yzzz_xxyy[k] = g1_zzz_xxyyy[k] - aby * g2_zzz_xxyy[k];

                        g_yzzz_xxyz[k] = g1_zzz_xxyyz[k] - aby * g2_zzz_xxyz[k];

                        g_yzzz_xxzz[k] = g1_zzz_xxyzz[k] - aby * g2_zzz_xxzz[k];

                        g_yzzz_xyyy[k] = g1_zzz_xyyyy[k] - aby * g2_zzz_xyyy[k];

                        g_yzzz_xyyz[k] = g1_zzz_xyyyz[k] - aby * g2_zzz_xyyz[k];

                        g_yzzz_xyzz[k] = g1_zzz_xyyzz[k] - aby * g2_zzz_xyzz[k];

                        g_yzzz_xzzz[k] = g1_zzz_xyzzz[k] - aby * g2_zzz_xzzz[k];

                        g_yzzz_yyyy[k] = g1_zzz_yyyyy[k] - aby * g2_zzz_yyyy[k];

                        g_yzzz_yyyz[k] = g1_zzz_yyyyz[k] - aby * g2_zzz_yyyz[k];

                        g_yzzz_yyzz[k] = g1_zzz_yyyzz[k] - aby * g2_zzz_yyzz[k];

                        g_yzzz_yzzz[k] = g1_zzz_yyzzz[k] - aby * g2_zzz_yzzz[k];

                        g_yzzz_zzzz[k] = g1_zzz_yzzzz[k] - aby * g2_zzz_zzzz[k];

                        // leading z component

                        g_zzzz_xxxx[k] = g1_zzz_xxxxz[k] - abz * g2_zzz_xxxx[k];

                        g_zzzz_xxxy[k] = g1_zzz_xxxyz[k] - abz * g2_zzz_xxxy[k];

                        g_zzzz_xxxz[k] = g1_zzz_xxxzz[k] - abz * g2_zzz_xxxz[k];

                        g_zzzz_xxyy[k] = g1_zzz_xxyyz[k] - abz * g2_zzz_xxyy[k];

                        g_zzzz_xxyz[k] = g1_zzz_xxyzz[k] - abz * g2_zzz_xxyz[k];

                        g_zzzz_xxzz[k] = g1_zzz_xxzzz[k] - abz * g2_zzz_xxzz[k];

                        g_zzzz_xyyy[k] = g1_zzz_xyyyz[k] - abz * g2_zzz_xyyy[k];

                        g_zzzz_xyyz[k] = g1_zzz_xyyzz[k] - abz * g2_zzz_xyyz[k];

                        g_zzzz_xyzz[k] = g1_zzz_xyzzz[k] - abz * g2_zzz_xyzz[k];

                        g_zzzz_xzzz[k] = g1_zzz_xzzzz[k] - abz * g2_zzz_xzzz[k];

                        g_zzzz_yyyy[k] = g1_zzz_yyyyz[k] - abz * g2_zzz_yyyy[k];

                        g_zzzz_yyyz[k] = g1_zzz_yyyzz[k] - abz * g2_zzz_yyyz[k];

                        g_zzzz_yyzz[k] = g1_zzz_yyzzz[k] - abz * g2_zzz_yyzz[k];

                        g_zzzz_yzzz[k] = g1_zzz_yzzzz[k] - abz * g2_zzz_yzzz[k];

                        g_zzzz_zzzz[k] = g1_zzz_zzzzz[k] - abz * g2_zzz_zzzz[k];
                    }
                }
            }
        }
    }
    
    
} // brahrrfunc namespace
