//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFunc.hpp"

#include "MathConst.hpp"
#include "GenFunc.hpp"

namespace ediprecfunc { // ediprecfunc namespace
    
    void
    compElectricDipoleForSS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  abDistances,
                            const CMemBlock2D<double>&  pcDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // set up pointers to primitives data on ket side
        
        auto knorm = ketGtoBlock.getNormFactors();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // fetch up pi values
        
        auto fpi = mathconst::getPiValue();
        
        // get position of integrals in primitves buffer
        
        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});
        
        auto doff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bnorm[i];
            
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // set up primitives buffer data
            
            auto fovl = primBuffer.data(soff + idx);
            
            auto fdipx = primBuffer.data(doff + idx);
            
            auto fdipy = primBuffer.data(doff + bdim + idx);
            
            auto fdipz = primBuffer.data(doff + 2 * bdim + idx);
            
            #pragma omp simd aligned(fovl, fx, fz, knorm, abx, aby,\
                                     abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fovl[j] = fb * knorm[j] * std::pow(fpi * fx[j], 1.5)
                
                        * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                             abz[j] * abz[j]));
                
                fdipx[j] = pcx[j] * fovl[j];
                
                fdipy[j] = pcy[j] * fovl[j];
                
                fdipz[j] = pcz[j] * fovl[j];
            }
            
            idx++;
        }
    }
    
    void
    compElectricDipoleForSP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {0, 1, 0})) return;

        if (iContrGto == 0) printf(" * VRR: (0|M|1)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto doff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto d1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|M|P) integrals

            auto dx_0_x = primBuffer.data(doff + 3 * idx);

            auto dx_0_y = primBuffer.data(doff + 3 * idx + 1);

            auto dx_0_z = primBuffer.data(doff + 3 * idx + 2);

            auto dy_0_x = primBuffer.data(doff + 3 * bdim + 3 * idx);

            auto dy_0_y = primBuffer.data(doff + 3 * bdim + 3 * idx + 1);

            auto dy_0_z = primBuffer.data(doff + 3 * bdim + 3 * idx + 2);

            auto dz_0_x = primBuffer.data(doff + 6 * bdim + 3 * idx);

            auto dz_0_y = primBuffer.data(doff + 6 * bdim + 3 * idx + 1);

            auto dz_0_z = primBuffer.data(doff + 6 * bdim + 3 * idx + 2);

            // set up pointers to (S|M|S) integrals

            auto dx_0_0 = primBuffer.data(d1off + idx);

            auto dy_0_0 = primBuffer.data(d1off + bdim + idx);

            auto dz_0_0 = primBuffer.data(d1off + 2 * bdim + idx);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(soff + idx);

            #pragma omp simd aligned(pbx, pby, pbz, dx_0_x, dx_0_y, dx_0_z, dy_0_x,\
                                     dy_0_y, dy_0_z, dz_0_x, dz_0_y, dz_0_z, dx_0_0,\
                                     dy_0_0, dz_0_0, s_0_0: VLX_ALIGN)
             for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                // electric dipole integrals

                // leading x component

                double fr = pbx[j];

                dx_0_x[j] = fr * dx_0_0[j] + f2t * s_0_0[j];

                dy_0_x[j] = fr * dy_0_0[j];

                dz_0_x[j] = fr * dz_0_0[j];

                // leading y component

                fr = pby[j];

                dx_0_y[j] = fr * dx_0_0[j];

                dy_0_y[j] = fr * dy_0_0[j] + f2t * s_0_0[j];

                dz_0_y[j] = fr * dz_0_0[j];

                // leading z component

                fr = pbz[j];

                dx_0_z[j] = fr * dx_0_0[j];

                dy_0_z[j] = fr * dy_0_0[j];

                dz_0_z[j] = fr * dz_0_0[j] + f2t * s_0_0[j];
            }
            
            idx++; 
        }
    }
    
    void
    compElectricDipoleForSD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern
        
        if (!genfunc::isInVector(recPattern, {0, 2, 0})) return;
        
        if (iContrGto == 0) printf(" * VRR: (0|M|2)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // get position of integrals in primitves buffer
        
        auto doff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});
        
        auto d1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});
        
        auto d2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});
        
        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});
        
        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            // set up pointers to distances R(PB)
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            // set up pointers to (S|M|D) integrals
            
            auto dx_0_xx = primBuffer.data(doff + 6 * idx);
            
            auto dx_0_xy = primBuffer.data(doff + 6 * idx + 1);
            
            auto dx_0_xz = primBuffer.data(doff + 6 * idx + 2);
            
            auto dx_0_yy = primBuffer.data(doff + 6 * idx + 3);
            
            auto dx_0_yz = primBuffer.data(doff + 6 * idx + 4);
            
            auto dx_0_zz = primBuffer.data(doff + 6 * idx + 5);
            
            auto dy_0_xx = primBuffer.data(doff + 6 * bdim + 6 * idx);
            
            auto dy_0_xy = primBuffer.data(doff + 6 * bdim + 6 * idx + 1);
            
            auto dy_0_xz = primBuffer.data(doff + 6 * bdim + 6 * idx + 2);
            
            auto dy_0_yy = primBuffer.data(doff + 6 * bdim + 6 * idx + 3);
            
            auto dy_0_yz = primBuffer.data(doff + 6 * bdim + 6 * idx + 4);
            
            auto dy_0_zz = primBuffer.data(doff + 6 * bdim + 6 * idx + 5);
            
            auto dz_0_xx = primBuffer.data(doff + 12 * bdim + 6 * idx);
            
            auto dz_0_xy = primBuffer.data(doff + 12 * bdim + 6 * idx + 1);
            
            auto dz_0_xz = primBuffer.data(doff + 12 * bdim + 6 * idx + 2);
            
            auto dz_0_yy = primBuffer.data(doff + 12 * bdim + 6 * idx + 3);
            
            auto dz_0_yz = primBuffer.data(doff + 12 * bdim + 6 * idx + 4);
            
            auto dz_0_zz = primBuffer.data(doff + 12 * bdim + 6 * idx + 5);
            
            // set up pointers to (S|M|P) integrals
            
            auto dx_0_x = primBuffer.data(d1off + 3 * idx);
            
            auto dx_0_y = primBuffer.data(d1off + 3 * idx + 1);
            
            auto dx_0_z = primBuffer.data(d1off + 3 * idx + 2);
            
            auto dy_0_x = primBuffer.data(d1off + 3 * bdim + 3 * idx);
            
            auto dy_0_y = primBuffer.data(d1off + 3 * bdim + 3 * idx + 1);
            
            auto dy_0_z = primBuffer.data(d1off + 3 * bdim + 3 * idx + 2);
            
            auto dz_0_x = primBuffer.data(d1off + 6 * bdim + 3 * idx);
            
            auto dz_0_y = primBuffer.data(d1off + 6 * bdim + 3 * idx + 1);
            
            auto dz_0_z = primBuffer.data(d1off + 6 * bdim + 3 * idx + 2);
            
            // set up pointers to (S|M|S) integrals
            
            auto dx_0_0 = primBuffer.data(d2off + idx);
            
            auto dy_0_0 = primBuffer.data(d2off + bdim + idx);
            
            auto dz_0_0 = primBuffer.data(d2off + 2 * bdim + idx);
            
            // set up pointers to (S|P) integrals
            
            auto s_0_x = primBuffer.data(soff + 3 * idx);
            
            auto s_0_y = primBuffer.data(soff + 3 * idx + 1);
            
            auto s_0_z = primBuffer.data(soff + 3 * idx + 2);
            
            // set up pointers to (S|S) integrals
            
            auto s_0_0 = primBuffer.data(s1off + idx);
            
            #pragma omp simd aligned(fx, pbx, pby, pbz, dx_0_xx, dx_0_xy, dx_0_xz,\
                                     dx_0_yy, dx_0_yz, dx_0_zz, dy_0_xx, dy_0_xy,\
                                     dy_0_xz, dy_0_yy, dy_0_yz, dy_0_zz, dz_0_xx,\
                                     dz_0_xy, dz_0_xz, dz_0_yy, dz_0_yz, dz_0_zz,\
                                     dx_0_x, dx_0_y, dx_0_z, dy_0_x, dy_0_y, dy_0_z,\
                                     dz_0_x, dz_0_y, dz_0_z, dx_0_0, dy_0_0, dz_0_0,\
                                     s_0_x, s_0_y, s_0_z, s_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor
                
                double f2t = 0.50 * fx[j];
                
                // overlap integrals
                
                // leading x component
                
                s_0_x[j] = pbx[j] * s_0_0[j];
                
                // leading y component
                
                s_0_y[j] = pby[j] * s_0_0[j];
                
                // leading z component
                
                s_0_z[j] = pbz[j] * s_0_0[j];
                
                // electric dipole integrals
                
                // leading x component
                
                double fr = pbx[j];
                
                dx_0_xx[j] = fr * dx_0_x[j] + f2t * dx_0_0[j] + f2t * s_0_x[j];
                
                dy_0_xx[j] = fr * dy_0_x[j] + f2t * dy_0_0[j];
                
                dz_0_xx[j] = fr * dz_0_x[j] + f2t * dz_0_0[j];
                
                dx_0_xy[j] = fr * dx_0_y[j] + f2t * s_0_y[j];
                
                dy_0_xy[j] = fr * dy_0_y[j];
                
                dz_0_xy[j] = fr * dz_0_y[j];
                
                dx_0_xz[j] = fr * dx_0_z[j] + f2t * s_0_z[j];
                
                dy_0_xz[j] = fr * dy_0_z[j];
                
                dz_0_xz[j] = fr * dz_0_z[j];
                
                // leading y component
                
                fr = pby[j];
                
                dx_0_yy[j] = fr * dx_0_y[j] + f2t * dx_0_0[j];
                
                dy_0_yy[j] = fr * dy_0_y[j] + f2t * dy_0_0[j] + f2t * s_0_y[j];
                
                dz_0_yy[j] = fr * dz_0_y[j] + f2t * dz_0_0[j];
                
                dx_0_yz[j] = fr * dx_0_z[j];
                
                dy_0_yz[j] = fr * dy_0_z[j] + f2t * s_0_z[j];
                
                dz_0_yz[j] = fr * dz_0_z[j];
                
                // leading z component
                
                fr = pbz[j];
                
                dx_0_zz[j] = fr * dx_0_z[j] + f2t * dx_0_0[j];
                
                dy_0_zz[j] = fr * dy_0_z[j] + f2t * dy_0_0[j];
                
                dz_0_zz[j] = fr * dz_0_z[j] + f2t * dz_0_0[j] + f2t * s_0_z[j];
                
            }
            
            idx++;
        }
    }
    
    void
    compElectricDipoleForSF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {0, 3, 0})) return;

        if (iContrGto == 0) printf(" * VRR: (0|M|3)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto doff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto d1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto d2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|M|F) integrals

            auto dx_0_xxx = primBuffer.data(doff + 10 * idx);

            auto dx_0_xxy = primBuffer.data(doff + 10 * idx + 1);

            auto dx_0_xxz = primBuffer.data(doff + 10 * idx + 2);

            auto dx_0_xyy = primBuffer.data(doff + 10 * idx + 3);

            auto dx_0_xyz = primBuffer.data(doff + 10 * idx + 4);

            auto dx_0_xzz = primBuffer.data(doff + 10 * idx + 5);

            auto dx_0_yyy = primBuffer.data(doff + 10 * idx + 6);

            auto dx_0_yyz = primBuffer.data(doff + 10 * idx + 7);

            auto dx_0_yzz = primBuffer.data(doff + 10 * idx + 8);

            auto dx_0_zzz = primBuffer.data(doff + 10 * idx + 9);

            auto dy_0_xxx = primBuffer.data(doff + 10 * bdim + 10 * idx);

            auto dy_0_xxy = primBuffer.data(doff + 10 * bdim + 10 * idx + 1);

            auto dy_0_xxz = primBuffer.data(doff + 10 * bdim + 10 * idx + 2);

            auto dy_0_xyy = primBuffer.data(doff + 10 * bdim + 10 * idx + 3);

            auto dy_0_xyz = primBuffer.data(doff + 10 * bdim + 10 * idx + 4);

            auto dy_0_xzz = primBuffer.data(doff + 10 * bdim + 10 * idx + 5);

            auto dy_0_yyy = primBuffer.data(doff + 10 * bdim + 10 * idx + 6);

            auto dy_0_yyz = primBuffer.data(doff + 10 * bdim + 10 * idx + 7);

            auto dy_0_yzz = primBuffer.data(doff + 10 * bdim + 10 * idx + 8);

            auto dy_0_zzz = primBuffer.data(doff + 10 * bdim + 10 * idx + 9);

            auto dz_0_xxx = primBuffer.data(doff + 20 * bdim + 10 * idx);

            auto dz_0_xxy = primBuffer.data(doff + 20 * bdim + 10 * idx + 1);

            auto dz_0_xxz = primBuffer.data(doff + 20 * bdim + 10 * idx + 2);

            auto dz_0_xyy = primBuffer.data(doff + 20 * bdim + 10 * idx + 3);

            auto dz_0_xyz = primBuffer.data(doff + 20 * bdim + 10 * idx + 4);

            auto dz_0_xzz = primBuffer.data(doff + 20 * bdim + 10 * idx + 5);

            auto dz_0_yyy = primBuffer.data(doff + 20 * bdim + 10 * idx + 6);

            auto dz_0_yyz = primBuffer.data(doff + 20 * bdim + 10 * idx + 7);

            auto dz_0_yzz = primBuffer.data(doff + 20 * bdim + 10 * idx + 8);

            auto dz_0_zzz = primBuffer.data(doff + 20 * bdim + 10 * idx + 9);

            // set up pointers to (S|M|D) integrals

            auto dx_0_xx = primBuffer.data(d1off + 6 * idx);

            auto dx_0_xy = primBuffer.data(d1off + 6 * idx + 1);

            auto dx_0_xz = primBuffer.data(d1off + 6 * idx + 2);

            auto dx_0_yy = primBuffer.data(d1off + 6 * idx + 3);

            auto dx_0_yz = primBuffer.data(d1off + 6 * idx + 4);

            auto dx_0_zz = primBuffer.data(d1off + 6 * idx + 5);

            auto dy_0_xx = primBuffer.data(d1off + 6 * bdim + 6 * idx);

            auto dy_0_xy = primBuffer.data(d1off + 6 * bdim + 6 * idx + 1);

            auto dy_0_xz = primBuffer.data(d1off + 6 * bdim + 6 * idx + 2);

            auto dy_0_yy = primBuffer.data(d1off + 6 * bdim + 6 * idx + 3);

            auto dy_0_yz = primBuffer.data(d1off + 6 * bdim + 6 * idx + 4);

            auto dy_0_zz = primBuffer.data(d1off + 6 * bdim + 6 * idx + 5);

            auto dz_0_xx = primBuffer.data(d1off + 12 * bdim + 6 * idx);

            auto dz_0_xy = primBuffer.data(d1off + 12 * bdim + 6 * idx + 1);

            auto dz_0_xz = primBuffer.data(d1off + 12 * bdim + 6 * idx + 2);

            auto dz_0_yy = primBuffer.data(d1off + 12 * bdim + 6 * idx + 3);

            auto dz_0_yz = primBuffer.data(d1off + 12 * bdim + 6 * idx + 4);

            auto dz_0_zz = primBuffer.data(d1off + 12 * bdim + 6 * idx + 5);

            // set up pointers to (S|M|P) integrals

            auto dx_0_x = primBuffer.data(d2off + 3 * idx);

            auto dx_0_y = primBuffer.data(d2off + 3 * idx + 1);

            auto dx_0_z = primBuffer.data(d2off + 3 * idx + 2);

            auto dy_0_x = primBuffer.data(d2off + 3 * bdim + 3 * idx);

            auto dy_0_y = primBuffer.data(d2off + 3 * bdim + 3 * idx + 1);

            auto dy_0_z = primBuffer.data(d2off + 3 * bdim + 3 * idx + 2);

            auto dz_0_x = primBuffer.data(d2off + 6 * bdim + 3 * idx);

            auto dz_0_y = primBuffer.data(d2off + 6 * bdim + 3 * idx + 1);

            auto dz_0_z = primBuffer.data(d2off + 6 * bdim + 3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(soff + 6 * idx);

            auto s_0_xy = primBuffer.data(soff + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(soff + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(soff + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(soff + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(soff + 6 * idx + 5);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s1off + 3 * idx);

            auto s_0_y = primBuffer.data(s1off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s1off + 3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(s2off + idx);

            #pragma omp simd aligned(fx, pbx, pby, pbz, dx_0_xxx, dx_0_xxy, dx_0_xxz,\
                                     dx_0_xyy, dx_0_xyz, dx_0_xzz, dx_0_yyy, dx_0_yyz,\
                                     dx_0_yzz, dx_0_zzz, dy_0_xxx, dy_0_xxy, dy_0_xxz,\
                                     dy_0_xyy, dy_0_xyz, dy_0_xzz, dy_0_yyy, dy_0_yyz,\
                                     dy_0_yzz, dy_0_zzz, dz_0_xxx, dz_0_xxy, dz_0_xxz,\
                                     dz_0_xyy, dz_0_xyz, dz_0_xzz, dz_0_yyy, dz_0_yyz,\
                                     dz_0_yzz, dz_0_zzz, dx_0_xx, dx_0_xy, dx_0_xz,\
                                     dx_0_yy, dx_0_yz, dx_0_zz, dy_0_xx, dy_0_xy,\
                                     dy_0_xz, dy_0_yy, dy_0_yz, dy_0_zz, dz_0_xx,\
                                     dz_0_xy, dz_0_xz, dz_0_yy, dz_0_yz, dz_0_zz,\
                                     dx_0_x, dx_0_y, dx_0_z, dy_0_x, dy_0_y, dy_0_z,\
                                     dz_0_x, dz_0_y, dz_0_z, s_0_xx, s_0_xy, s_0_xz,\
                                     s_0_yy, s_0_yz, s_0_zz, s_0_x, s_0_y, s_0_z,\
                                     s_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                // overlap integrals

                // leading x component

                double fr = pbx[j];

                s_0_xx[j] = fr * s_0_x[j] + f2t * s_0_0[j];

                s_0_xy[j] = fr * s_0_y[j];

                s_0_xz[j] = fr * s_0_z[j];

                // leading y component

                fr = pby[j];

                s_0_yy[j] = fr * s_0_y[j] + f2t * s_0_0[j];

                s_0_yz[j] = fr * s_0_z[j];

                // leading z component

                s_0_zz[j] = pbz[j] * s_0_z[j] + f2t * s_0_0[j];

                // electric dipole integrals

                // leading x component

                fr = pbx[j];

                dx_0_xxx[j] = fr * dx_0_xx[j] + 2.0 * f2t * dx_0_x[j] + f2t * s_0_xx[j];

                dy_0_xxx[j] = fr * dy_0_xx[j] + 2.0 * f2t * dy_0_x[j];

                dz_0_xxx[j] = fr * dz_0_xx[j] + 2.0 * f2t * dz_0_x[j];

                dx_0_xxy[j] = fr * dx_0_xy[j] + f2t * dx_0_y[j] + f2t * s_0_xy[j];

                dy_0_xxy[j] = fr * dy_0_xy[j] + f2t * dy_0_y[j];

                dz_0_xxy[j] = fr * dz_0_xy[j] + f2t * dz_0_y[j];

                dx_0_xxz[j] = fr * dx_0_xz[j] + f2t * dx_0_z[j] + f2t * s_0_xz[j];

                dy_0_xxz[j] = fr * dy_0_xz[j] + f2t * dy_0_z[j];

                dz_0_xxz[j] = fr * dz_0_xz[j] + f2t * dz_0_z[j];

                dx_0_xyy[j] = fr * dx_0_yy[j] + f2t * s_0_yy[j];

                dy_0_xyy[j] = fr * dy_0_yy[j];

                dz_0_xyy[j] = fr * dz_0_yy[j];

                dx_0_xyz[j] = fr * dx_0_yz[j] + f2t * s_0_yz[j];

                dy_0_xyz[j] = fr * dy_0_yz[j];

                dz_0_xyz[j] = fr * dz_0_yz[j];

                dx_0_xzz[j] = fr * dx_0_zz[j] + f2t * s_0_zz[j];

                dy_0_xzz[j] = fr * dy_0_zz[j];

                dz_0_xzz[j] = fr * dz_0_zz[j];

                // leading y component

                fr = pby[j];

                dx_0_yyy[j] = fr * dx_0_yy[j] + 2.0 * f2t * dx_0_y[j];

                dy_0_yyy[j] = fr * dy_0_yy[j] + 2.0 * f2t * dy_0_y[j] + f2t * s_0_yy[j];

                dz_0_yyy[j] = fr * dz_0_yy[j] + 2.0 * f2t * dz_0_y[j];

                dx_0_yyz[j] = fr * dx_0_yz[j] + f2t * dx_0_z[j];

                dy_0_yyz[j] = fr * dy_0_yz[j] + f2t * dy_0_z[j] + f2t * s_0_yz[j];

                dz_0_yyz[j] = fr * dz_0_yz[j] + f2t * dz_0_z[j];

                dx_0_yzz[j] = fr * dx_0_zz[j];

                dy_0_yzz[j] = fr * dy_0_zz[j] + f2t * s_0_zz[j];

                dz_0_yzz[j] = fr * dz_0_zz[j];

                // leading z component

                fr = pbz[j];

                dx_0_zzz[j] = fr * dx_0_zz[j] + 2.0 * f2t * dx_0_z[j];

                dy_0_zzz[j] = fr * dy_0_zz[j] + 2.0 * f2t * dy_0_z[j];

                dz_0_zzz[j] = fr * dz_0_zz[j] + 2.0 * f2t * dz_0_z[j] + f2t * s_0_zz[j];

            }

            idx++;
        }
    }
    
    void
    compElectricDipoleForSG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {0, 4, 0})) return;

        if (iContrGto == 0) printf(" * VRR: (0|M|4)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto doff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 0});

        auto d1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto d2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|M|G) integrals

            auto dx_0_xxxx = primBuffer.data(doff + 15 * idx);

            auto dx_0_xxxy = primBuffer.data(doff + 15 * idx + 1);

            auto dx_0_xxxz = primBuffer.data(doff + 15 * idx + 2);

            auto dx_0_xxyy = primBuffer.data(doff + 15 * idx + 3);

            auto dx_0_xxyz = primBuffer.data(doff + 15 * idx + 4);

            auto dx_0_xxzz = primBuffer.data(doff + 15 * idx + 5);

            auto dx_0_xyyy = primBuffer.data(doff + 15 * idx + 6);

            auto dx_0_xyyz = primBuffer.data(doff + 15 * idx + 7);

            auto dx_0_xyzz = primBuffer.data(doff + 15 * idx + 8);

            auto dx_0_xzzz = primBuffer.data(doff + 15 * idx + 9);

            auto dx_0_yyyy = primBuffer.data(doff + 15 * idx + 10);

            auto dx_0_yyyz = primBuffer.data(doff + 15 * idx + 11);

            auto dx_0_yyzz = primBuffer.data(doff + 15 * idx + 12);

            auto dx_0_yzzz = primBuffer.data(doff + 15 * idx + 13);

            auto dx_0_zzzz = primBuffer.data(doff + 15 * idx + 14);

            auto dy_0_xxxx = primBuffer.data(doff + 15 * bdim + 15 * idx);

            auto dy_0_xxxy = primBuffer.data(doff + 15 * bdim + 15 * idx + 1);

            auto dy_0_xxxz = primBuffer.data(doff + 15 * bdim + 15 * idx + 2);

            auto dy_0_xxyy = primBuffer.data(doff + 15 * bdim + 15 * idx + 3);

            auto dy_0_xxyz = primBuffer.data(doff + 15 * bdim + 15 * idx + 4);

            auto dy_0_xxzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 5);

            auto dy_0_xyyy = primBuffer.data(doff + 15 * bdim + 15 * idx + 6);

            auto dy_0_xyyz = primBuffer.data(doff + 15 * bdim + 15 * idx + 7);

            auto dy_0_xyzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 8);

            auto dy_0_xzzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 9);

            auto dy_0_yyyy = primBuffer.data(doff + 15 * bdim + 15 * idx + 10);

            auto dy_0_yyyz = primBuffer.data(doff + 15 * bdim + 15 * idx + 11);

            auto dy_0_yyzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 12);

            auto dy_0_yzzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 13);

            auto dy_0_zzzz = primBuffer.data(doff + 15 * bdim + 15 * idx + 14);

            auto dz_0_xxxx = primBuffer.data(doff + 30 * bdim + 15 * idx);

            auto dz_0_xxxy = primBuffer.data(doff + 30 * bdim + 15 * idx + 1);

            auto dz_0_xxxz = primBuffer.data(doff + 30 * bdim + 15 * idx + 2);

            auto dz_0_xxyy = primBuffer.data(doff + 30 * bdim + 15 * idx + 3);

            auto dz_0_xxyz = primBuffer.data(doff + 30 * bdim + 15 * idx + 4);

            auto dz_0_xxzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 5);

            auto dz_0_xyyy = primBuffer.data(doff + 30 * bdim + 15 * idx + 6);

            auto dz_0_xyyz = primBuffer.data(doff + 30 * bdim + 15 * idx + 7);

            auto dz_0_xyzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 8);

            auto dz_0_xzzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 9);

            auto dz_0_yyyy = primBuffer.data(doff + 30 * bdim + 15 * idx + 10);

            auto dz_0_yyyz = primBuffer.data(doff + 30 * bdim + 15 * idx + 11);

            auto dz_0_yyzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 12);

            auto dz_0_yzzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 13);

            auto dz_0_zzzz = primBuffer.data(doff + 30 * bdim + 15 * idx + 14);

            // set up pointers to (S|M|F) integrals

            auto dx_0_xxx = primBuffer.data(d1off + 10 * idx);

            auto dx_0_xxy = primBuffer.data(d1off + 10 * idx + 1);

            auto dx_0_xxz = primBuffer.data(d1off + 10 * idx + 2);

            auto dx_0_xyy = primBuffer.data(d1off + 10 * idx + 3);

            auto dx_0_xyz = primBuffer.data(d1off + 10 * idx + 4);

            auto dx_0_xzz = primBuffer.data(d1off + 10 * idx + 5);

            auto dx_0_yyy = primBuffer.data(d1off + 10 * idx + 6);

            auto dx_0_yyz = primBuffer.data(d1off + 10 * idx + 7);

            auto dx_0_yzz = primBuffer.data(d1off + 10 * idx + 8);

            auto dx_0_zzz = primBuffer.data(d1off + 10 * idx + 9);

            auto dy_0_xxx = primBuffer.data(d1off + 10 * bdim + 10 * idx);

            auto dy_0_xxy = primBuffer.data(d1off + 10 * bdim + 10 * idx + 1);

            auto dy_0_xxz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 2);

            auto dy_0_xyy = primBuffer.data(d1off + 10 * bdim + 10 * idx + 3);

            auto dy_0_xyz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 4);

            auto dy_0_xzz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 5);

            auto dy_0_yyy = primBuffer.data(d1off + 10 * bdim + 10 * idx + 6);

            auto dy_0_yyz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 7);

            auto dy_0_yzz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 8);

            auto dy_0_zzz = primBuffer.data(d1off + 10 * bdim + 10 * idx + 9);

            auto dz_0_xxx = primBuffer.data(d1off + 20 * bdim + 10 * idx);

            auto dz_0_xxy = primBuffer.data(d1off + 20 * bdim + 10 * idx + 1);

            auto dz_0_xxz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 2);

            auto dz_0_xyy = primBuffer.data(d1off + 20 * bdim + 10 * idx + 3);

            auto dz_0_xyz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 4);

            auto dz_0_xzz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 5);

            auto dz_0_yyy = primBuffer.data(d1off + 20 * bdim + 10 * idx + 6);

            auto dz_0_yyz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 7);

            auto dz_0_yzz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 8);

            auto dz_0_zzz = primBuffer.data(d1off + 20 * bdim + 10 * idx + 9);

            // set up pointers to (S|M|D) integrals

            auto dx_0_xx = primBuffer.data(d2off + 6 * idx);

            auto dx_0_xy = primBuffer.data(d2off + 6 * idx + 1);

            auto dx_0_xz = primBuffer.data(d2off + 6 * idx + 2);

            auto dx_0_yy = primBuffer.data(d2off + 6 * idx + 3);

            auto dx_0_yz = primBuffer.data(d2off + 6 * idx + 4);

            auto dx_0_zz = primBuffer.data(d2off + 6 * idx + 5);

            auto dy_0_xx = primBuffer.data(d2off + 6 * bdim + 6 * idx);

            auto dy_0_xy = primBuffer.data(d2off + 6 * bdim + 6 * idx + 1);

            auto dy_0_xz = primBuffer.data(d2off + 6 * bdim + 6 * idx + 2);

            auto dy_0_yy = primBuffer.data(d2off + 6 * bdim + 6 * idx + 3);

            auto dy_0_yz = primBuffer.data(d2off + 6 * bdim + 6 * idx + 4);

            auto dy_0_zz = primBuffer.data(d2off + 6 * bdim + 6 * idx + 5);

            auto dz_0_xx = primBuffer.data(d2off + 12 * bdim + 6 * idx);

            auto dz_0_xy = primBuffer.data(d2off + 12 * bdim + 6 * idx + 1);

            auto dz_0_xz = primBuffer.data(d2off + 12 * bdim + 6 * idx + 2);

            auto dz_0_yy = primBuffer.data(d2off + 12 * bdim + 6 * idx + 3);

            auto dz_0_yz = primBuffer.data(d2off + 12 * bdim + 6 * idx + 4);

            auto dz_0_zz = primBuffer.data(d2off + 12 * bdim + 6 * idx + 5);

            // set up pointers to (S|F) integrals

            auto s_0_xxx = primBuffer.data(soff + 10 * idx);

            auto s_0_xxy = primBuffer.data(soff + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(soff + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(soff + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(soff + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(soff + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(soff + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(soff + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(soff + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(soff + 10 * idx + 9);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(s1off + 6 * idx);

            auto s_0_xy = primBuffer.data(s1off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(s1off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(s1off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(s1off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(s1off + 6 * idx + 5);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s2off + 3 * idx);

            auto s_0_y = primBuffer.data(s2off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s2off + 3 * idx + 2);

            #pragma omp simd aligned(fx, pbx, pby, pbz, dx_0_xxxx, dx_0_xxxy, dx_0_xxxz,\
                                     dx_0_xxyy, dx_0_xxyz, dx_0_xxzz, dx_0_xyyy,\
                                     dx_0_xyyz, dx_0_xyzz, dx_0_xzzz, dx_0_yyyy,\
                                     dx_0_yyyz, dx_0_yyzz, dx_0_yzzz, dx_0_zzzz,\
                                     dy_0_xxxx, dy_0_xxxy, dy_0_xxxz, dy_0_xxyy,\
                                     dy_0_xxyz, dy_0_xxzz, dy_0_xyyy, dy_0_xyyz,\
                                     dy_0_xyzz, dy_0_xzzz, dy_0_yyyy, dy_0_yyyz,\
                                     dy_0_yyzz, dy_0_yzzz, dy_0_zzzz, dz_0_xxxx,\
                                     dz_0_xxxy, dz_0_xxxz, dz_0_xxyy, dz_0_xxyz,\
                                     dz_0_xxzz, dz_0_xyyy, dz_0_xyyz, dz_0_xyzz,\
                                     dz_0_xzzz, dz_0_yyyy, dz_0_yyyz, dz_0_yyzz,\
                                     dz_0_yzzz, dz_0_zzzz, dx_0_xxx, dx_0_xxy,\
                                     dx_0_xxz, dx_0_xyy, dx_0_xyz, dx_0_xzz, dx_0_yyy,\
                                     dx_0_yyz, dx_0_yzz, dx_0_zzz, dy_0_xxx, dy_0_xxy,\
                                     dy_0_xxz, dy_0_xyy, dy_0_xyz, dy_0_xzz, dy_0_yyy,\
                                     dy_0_yyz, dy_0_yzz, dy_0_zzz, dz_0_xxx, dz_0_xxy,\
                                     dz_0_xxz, dz_0_xyy, dz_0_xyz, dz_0_xzz, dz_0_yyy,\
                                     dz_0_yyz, dz_0_yzz, dz_0_zzz, dx_0_xx, dx_0_xy,\
                                     dx_0_xz, dx_0_yy, dx_0_yz, dx_0_zz, dy_0_xx,\
                                     dy_0_xy, dy_0_xz, dy_0_yy, dy_0_yz, dy_0_zz,\
                                     dz_0_xx, dz_0_xy, dz_0_xz, dz_0_yy, dz_0_yz,\
                                     dz_0_zz, s_0_xxx, s_0_xxy, s_0_xxz, s_0_xyy,\
                                     s_0_xyz, s_0_xzz, s_0_yyy, s_0_yyz, s_0_yzz,\
                                     s_0_zzz, s_0_xx, s_0_xy, s_0_xz, s_0_yy, s_0_yz,\
                                     s_0_zz, s_0_x, s_0_y, s_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                // overlap integrals

                // leading x component

                double fr = pbx[j];

                s_0_xxx[j] = fr * s_0_xx[j] + 2.0 * f2t * s_0_x[j];

                s_0_xxy[j] = fr * s_0_xy[j] + f2t * s_0_y[j];

                s_0_xxz[j] = fr * s_0_xz[j] + f2t * s_0_z[j];

                s_0_xyy[j] = fr * s_0_yy[j];

                s_0_xyz[j] = fr * s_0_yz[j];

                s_0_xzz[j] = fr * s_0_zz[j];

                // leading y component

                fr = pby[j];

                s_0_yyy[j] = fr * s_0_yy[j] + 2.0 * f2t * s_0_y[j];

                s_0_yyz[j] = fr * s_0_yz[j] + f2t * s_0_z[j];

                s_0_yzz[j] = fr * s_0_zz[j];

                // leading z component

                s_0_zzz[j] = pbz[j] * s_0_zz[j] + 2.0 * f2t * s_0_z[j];

                // electric dipole integrals

                // leading x component

                fr = pbx[j];

                dx_0_xxxx[j] = fr * dx_0_xxx[j] + 3.0 * f2t * dx_0_xx[j] + f2t * s_0_xxx[j];

                dy_0_xxxx[j] = fr * dy_0_xxx[j] + 3.0 * f2t * dy_0_xx[j];

                dz_0_xxxx[j] = fr * dz_0_xxx[j] + 3.0 * f2t * dz_0_xx[j];

                dx_0_xxxy[j] = fr * dx_0_xxy[j] + 2.0 * f2t * dx_0_xy[j] + f2t * s_0_xxy[j];

                dy_0_xxxy[j] = fr * dy_0_xxy[j] + 2.0 * f2t * dy_0_xy[j];

                dz_0_xxxy[j] = fr * dz_0_xxy[j] + 2.0 * f2t * dz_0_xy[j];

                dx_0_xxxz[j] = fr * dx_0_xxz[j] + 2.0 * f2t * dx_0_xz[j] + f2t * s_0_xxz[j];

                dy_0_xxxz[j] = fr * dy_0_xxz[j] + 2.0 * f2t * dy_0_xz[j];

                dz_0_xxxz[j] = fr * dz_0_xxz[j] + 2.0 * f2t * dz_0_xz[j];

                dx_0_xxyy[j] = fr * dx_0_xyy[j] + f2t * dx_0_yy[j] + f2t * s_0_xyy[j];

                dy_0_xxyy[j] = fr * dy_0_xyy[j] + f2t * dy_0_yy[j];

                dz_0_xxyy[j] = fr * dz_0_xyy[j] + f2t * dz_0_yy[j];

                dx_0_xxyz[j] = fr * dx_0_xyz[j] + f2t * dx_0_yz[j] + f2t * s_0_xyz[j];

                dy_0_xxyz[j] = fr * dy_0_xyz[j] + f2t * dy_0_yz[j];

                dz_0_xxyz[j] = fr * dz_0_xyz[j] + f2t * dz_0_yz[j];

                dx_0_xxzz[j] = fr * dx_0_xzz[j] + f2t * dx_0_zz[j] + f2t * s_0_xzz[j];

                dy_0_xxzz[j] = fr * dy_0_xzz[j] + f2t * dy_0_zz[j];

                dz_0_xxzz[j] = fr * dz_0_xzz[j] + f2t * dz_0_zz[j];

                dx_0_xyyy[j] = fr * dx_0_yyy[j] + f2t * s_0_yyy[j];

                dy_0_xyyy[j] = fr * dy_0_yyy[j];

                dz_0_xyyy[j] = fr * dz_0_yyy[j];

                dx_0_xyyz[j] = fr * dx_0_yyz[j] + f2t * s_0_yyz[j];

                dy_0_xyyz[j] = fr * dy_0_yyz[j];

                dz_0_xyyz[j] = fr * dz_0_yyz[j];

                dx_0_xyzz[j] = fr * dx_0_yzz[j] + f2t * s_0_yzz[j];

                dy_0_xyzz[j] = fr * dy_0_yzz[j];

                dz_0_xyzz[j] = fr * dz_0_yzz[j];

                dx_0_xzzz[j] = fr * dx_0_zzz[j] + f2t * s_0_zzz[j];

                dy_0_xzzz[j] = fr * dy_0_zzz[j];

                dz_0_xzzz[j] = fr * dz_0_zzz[j];

                // leading y component

                fr = pby[j];

                dx_0_yyyy[j] = fr * dx_0_yyy[j] + 3.0 * f2t * dx_0_yy[j];

                dy_0_yyyy[j] = fr * dy_0_yyy[j] + 3.0 * f2t * dy_0_yy[j] + f2t * s_0_yyy[j];

                dz_0_yyyy[j] = fr * dz_0_yyy[j] + 3.0 * f2t * dz_0_yy[j];

                dx_0_yyyz[j] = fr * dx_0_yyz[j] + 2.0 * f2t * dx_0_yz[j];

                dy_0_yyyz[j] = fr * dy_0_yyz[j] + 2.0 * f2t * dy_0_yz[j] + f2t * s_0_yyz[j];

                dz_0_yyyz[j] = fr * dz_0_yyz[j] + 2.0 * f2t * dz_0_yz[j];

                dx_0_yyzz[j] = fr * dx_0_yzz[j] + f2t * dx_0_zz[j];

                dy_0_yyzz[j] = fr * dy_0_yzz[j] + f2t * dy_0_zz[j] + f2t * s_0_yzz[j];

                dz_0_yyzz[j] = fr * dz_0_yzz[j] + f2t * dz_0_zz[j];

                dx_0_yzzz[j] = fr * dx_0_zzz[j];

                dy_0_yzzz[j] = fr * dy_0_zzz[j] + f2t * s_0_zzz[j];

                dz_0_yzzz[j] = fr * dz_0_zzz[j];

                // leading z component

                fr = pbz[j];

                dx_0_zzzz[j] = fr * dx_0_zzz[j] + 3.0 * f2t * dx_0_zz[j];

                dy_0_zzzz[j] = fr * dy_0_zzz[j] + 3.0 * f2t * dy_0_zz[j];

                dz_0_zzzz[j] = fr * dz_0_zzz[j] + 3.0 * f2t * dz_0_zz[j] + f2t * s_0_zzz[j];

            }

            idx++;
        }
    }
    
} // ediprecfunc namespace
