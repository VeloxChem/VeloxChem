//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "NuclearPotentialRecFunc.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "GenFunc.hpp"

namespace npotrecfunc { // npotrecfunc namespace

    void
    compNuclearPotentialForSS(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CBoysFunction&        bfTable,
                                    CMemBlock<double>&    bfArguments,
                                    CMemBlock2D<double>&  bfValues,
                              const int32_t               bfOrder,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  abDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        if (iContrGto == 0) printf(" ==> VRR(0, 0)\n");
        
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto knorm = ketGtoBlock.getNormFactors();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up Obara-Saika prefactors
            
            auto fg = osFactors.data(3 * idx + 2);
            
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // compute Boys function argument
            
            auto fargs = bfArguments.data();
            
            #pragma omp simd aligned(fg, pcx, pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fg[j] * (pcx[j] * pcx[j] + pcy[j] * pcy[j] +
                                    
                                    pcz[j] * pcz[j]);
            }
            
            // evaluate Boys function values
            
            bfTable.compute(bfValues, bfArguments, bfOrder);
            
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(3 * idx);
            
            auto fz = osFactors.data(3 * idx + 1);
            
            auto fb = bnorm[i];
            
            // fetch up pi values
            
            auto fpi = 2.0 * mathconst::getPiValue();
            
            // compute overlap scaling factor

            #pragma omp simd aligned(fx, fz, knorm, abx, aby, abz,\
                                     fargs: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fb * knorm[j] * fpi * fx[j]
                
                         * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                              abz[j] * abz[j]));
            }
            
            // distribute (s|A(0)|s) integrals
            
            for (int32_t j = 0; j <= bfOrder; j++)
            {
                auto pidx = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {0, 0, j});
                
                auto t_0_0 = primBuffer.data(pidx + idx);
                
                auto bvals = bfValues.data(j);
                
                #pragma omp simd aligned(t_0_0, bvals, fargs: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    t_0_0[k] = bvals[k] * fargs[k];
                }
            }
            
            idx++;
        }
    }
    
    void
    compNuclearPotentialForSP(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  pbDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern
        
        if (!genfunc::isInVector(recPattern, {0, 1, 0})) return;
        
        if (iContrGto == 0) printf(" ==> VRR(0, 1)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute primitive integrals up to required order 
        
        auto bord = genfunc::maxOrderOfPair(recPattern, 0, 1);
        
        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer
            
            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 1, i});
            
            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i});
            
            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i + 1});
            
            // loop over contracted GTO on bra side
            
            int32_t idx = 0;
            
            //std::cout << "### Prim buffer: " << primBuffer;
            
            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to distances R(PB)
                
                auto pbx = pbDistances.data(3 * idx);
                
                auto pby = pbDistances.data(3 * idx + 1);
                
                auto pbz = pbDistances.data(3 * idx + 2);
                
                // set up pointers to distances R(PC)
                
                auto pcx = pcDistances.data(3 * idx);
                
                auto pcy = pcDistances.data(3 * idx + 1);
                
                auto pcz = pcDistances.data(3 * idx + 2);
                
                // set up pointers to (S|A(0)|S)^(m) integrals
                
                auto t00_0_0 = primBuffer.data(t10off + idx);
                
                // set up pointers to (S|A(0)|S)^(m+1) integrals
                
                auto t01_0_0 = primBuffer.data(t11off + idx);
                
                // set up pointers to (S|A(0)|P)^(m) integrals
                
                auto t_0_x = primBuffer.data(toff + 3 * idx);
                
                auto t_0_y = primBuffer.data(toff + 3 * idx + 1);
                
                auto t_0_z = primBuffer.data(toff + 3 * idx + 2);
                
                #pragma omp simd aligned(pbx, pby, pbz, pcx, pcy, pcz, t_0_x, \
                                         t_0_y, t_0_z, t00_0_0, t01_0_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    double fpb = t00_0_0[k];
                    
                    double fpc = t01_0_0[k];
                    
                    t_0_x[k] =  fpb * pbx[k] - fpc * pcx[k];
                    
                    t_0_y[k] =  fpb * pby[k] - fpc * pcy[k];
                    
                    t_0_z[k] =  fpb * pbz[k] - fpc * pcz[k];
                }
                
                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForPS(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern
        
        if (!genfunc::isInVector(recPattern, {1, 0, 0})) return;
        
        if (iContrGto == 0) printf(" ==> VRR(1, 0)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute primitive integrals up to required order
        
        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 0);
        
        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer
            
            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {1, 0, i});
            
            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i});
            
            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i + 1});
            
            // loop over contracted GTO on bra side
            
            int32_t idx = 0;
            
            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to distances R(PA)
                
                auto pax = paDistances.data(3 * idx);
                
                auto pay = paDistances.data(3 * idx + 1);
                
                auto paz = paDistances.data(3 * idx + 2);
                
                // set up pointers to distances R(PC)
                
                auto pcx = pcDistances.data(3 * idx);
                
                auto pcy = pcDistances.data(3 * idx + 1);
                
                auto pcz = pcDistances.data(3 * idx + 2);
                
                // set up pointers to (S|A(0)|S)^(m) integrals
                
                auto t00_0_0 = primBuffer.data(t10off + idx);
                
                // set up pointers to (S|A(0)|S)^(m+1) integrals
                
                auto t01_0_0 = primBuffer.data(t11off + idx);
                
                // set up pointers to (S|A(0)|P)^(m) integrals
                
                auto t_x_0 = primBuffer.data(toff + 3 * idx);
                
                auto t_y_0 = primBuffer.data(toff + 3 * idx + 1);
                
                auto t_z_0 = primBuffer.data(toff + 3 * idx + 2);
                
                #pragma omp simd aligned(pax, pay, paz, pcx, pcy, pcz, t_x_0, \
                                         t_y_0, t_z_0, t00_0_0, t01_0_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    double fpa = t00_0_0[k];
                    
                    double fpc = t01_0_0[k];
                    
                    t_x_0[k] = fpa * pax[k] - fpc * pcx[k];
                    
                    t_y_0[k] = fpa * pay[k] - fpc * pcy[k];
                    
                    t_z_0[k] = fpa * paz[k] - fpc * pcz[k];
                }
                
                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForPP(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern
        
        if (!genfunc::isInVector(recPattern, {1, 1, 0})) return;
        
        if (iContrGto == 0) printf(" ==> VRR(1, 1)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute primitive integrals up to required order
        
        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 1);
        
        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer
            
            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {1, 1, i});
            
            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 1, i});
            
            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 1, i + 1});
            
            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i});
            
            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                   {0, 0, i + 1});
            
            // loop over contracted GTO on bra side
            
            int32_t idx = 0;
            
            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors
                
                auto fx = osFactors.data(3 * idx);
                
                // set up pointers to distances R(PA)
                
                auto pax = paDistances.data(3 * idx);
                
                auto pay = paDistances.data(3 * idx + 1);
                
                auto paz = paDistances.data(3 * idx + 2);
                
                // set up pointers to distances R(PC)
                
                auto pcx = pcDistances.data(3 * idx);
                
                auto pcy = pcDistances.data(3 * idx + 1);
                
                auto pcz = pcDistances.data(3 * idx + 2);
                
                // set up pointers to (S|A(0)|S)^(m) integrals
                
                auto tk0_0_0 = primBuffer.data(tk0off + idx);
                
                // set up pointers to (S|A(0)|S)^(m+1) integrals
                
                auto tk1_0_0 = primBuffer.data(tk1off + idx);
                
                // set up pointers to (S|A(0)|P)^(m) integrals
                
                auto t10_0_x = primBuffer.data(t10off + 3 * idx);
                
                auto t10_0_y = primBuffer.data(t10off + 3 * idx + 1);
                
                auto t10_0_z = primBuffer.data(t10off + 3 * idx + 2);
                
                // set up pointers to (S|A(0)|P)^(m+1) integrals
                
                auto t11_0_x = primBuffer.data(t11off + 3 * idx);
                
                auto t11_0_y = primBuffer.data(t11off + 3 * idx + 1);
                
                auto t11_0_z = primBuffer.data(t11off + 3 * idx + 2);
                
                // set up pointers to (P|A(0)|P)^(m) integrals
                
                auto t_x_x = primBuffer.data(toff + 9 * idx);
                
                auto t_x_y = primBuffer.data(toff + 9 * idx + 1);
                
                auto t_x_z = primBuffer.data(toff + 9 * idx + 2);
                
                auto t_y_x = primBuffer.data(toff + 9 * idx + 3);
                
                auto t_y_y = primBuffer.data(toff + 9 * idx + 4);
                
                auto t_y_z = primBuffer.data(toff + 9 * idx + 5);
                
                auto t_z_x = primBuffer.data(toff + 9 * idx + 6);
                
                auto t_z_y = primBuffer.data(toff + 9 * idx + 7);
                
                auto t_z_z = primBuffer.data(toff + 9 * idx + 8);
                
                #pragma omp simd aligned(fx, pax, pay, paz, pcx, pcy, pcz,\
                                         tk0_0_0, tk1_0_0, t10_0_x, t10_0_y,\
                                         t10_0_z, t11_0_x, t11_0_y, t11_0_z,\
                                         t_x_x, t_x_y, t_x_z, t_y_x, t_y_y,\
                                         t_y_z, t_z_x, t_z_y, t_z_z: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor
                    
                    double f2t = 0.50 * fx[k] * (tk0_0_0[k] - tk1_0_0[k]);
                    
                    // leading x component
                    
                    double fra = pax[k];
                    
                    double frc = pcx[k];
                    
                    t_x_x[k] = fra * t10_0_x[k] - frc * t11_0_x[k] + f2t;
                    
                    t_x_y[k] = fra * t10_0_y[k] - frc * t11_0_y[k];
                    
                    t_x_z[k] = fra * t10_0_z[k] - frc * t11_0_z[k];
                    
                    // leading y component
                    
                    fra = pay[k];
                    
                    frc = pcy[k];
                    
                    t_y_x[k] = fra * t10_0_x[k] - frc * t11_0_x[k];
                    
                    t_y_y[k] = fra * t10_0_y[k] - frc * t11_0_y[k] + f2t;
                    
                    t_y_z[k] = fra * t10_0_z[k] - frc * t11_0_z[k];
                    
                    // leading z component
                    
                    fra = paz[k];
                    
                    frc = pcz[k];
                    
                    t_z_x[k] = fra * t10_0_x[k] - frc * t11_0_x[k];
                    
                    t_z_y[k] = fra * t10_0_y[k] - frc * t11_0_y[k];
                    
                    t_z_z[k] = fra * t10_0_z[k] - frc * t11_0_z[k] + f2t;
                }
                
                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForSD(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  pbDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {0, 2, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (0|2)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 0, 2);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PB)

                auto pbx = pbDistances.data(3 * idx);

                auto pby = pbDistances.data(3 * idx + 1);

                auto pbz = pbDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (S|A(0)|S)^(m) integrals

                auto t20_0_0 = primBuffer.data(t20off + idx);

                // set up pointers to (S|A(0)|S)^(m+1) integrals

                auto t21_0_0 = primBuffer.data(t21off + idx);

                // set up pointers to (S|A(0)|P)^(m) integrals

                auto t10_0_x = primBuffer.data(t10off + 3 * idx);

                auto t10_0_y = primBuffer.data(t10off + 3 * idx + 1);

                auto t10_0_z = primBuffer.data(t10off + 3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m+1) integrals

                auto t11_0_x = primBuffer.data(t11off + 3 * idx);

                auto t11_0_y = primBuffer.data(t11off + 3 * idx + 1);

                auto t11_0_z = primBuffer.data(t11off + 3 * idx + 2);

                // set up pointers to (D|A(0)|S)^(m) integrals

                auto t_0_xx = primBuffer.data(toff + 6 * idx);

                auto t_0_xy = primBuffer.data(toff + 6 * idx + 1);

                auto t_0_xz = primBuffer.data(toff + 6 * idx + 2);

                auto t_0_yy = primBuffer.data(toff + 6 * idx + 3);

                auto t_0_yz = primBuffer.data(toff + 6 * idx + 4);

                auto t_0_zz = primBuffer.data(toff + 6 * idx + 5);

                #pragma omp simd aligned(fx, pbx, pby, pbz, t20_0_0, t21_0_0, t10_0_x,\
                                         t10_0_y, t10_0_z, t11_0_x, t11_0_y, t11_0_z,\
                                         t_0_xx, t_0_xy, t_0_xz, t_0_yy, t_0_yz,\
                                         t_0_zz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k] * (t20_0_0[k] - t21_0_0[k]);

                    // leading x component

                    double frb = pbx[k];

                    double frc = pcx[k];

                    t_0_xx[k] = frb * t10_0_x[k] - frc * t11_0_x[k] + f2t;

                    t_0_xy[k] = frb * t10_0_y[k] - frc * t11_0_y[k];

                    t_0_xz[k] = frb * t10_0_z[k] - frc * t11_0_z[k];

                    // leading y component

                    frb = pby[k];

                    frc = pcy[k];

                    t_0_yy[k] = frb * t10_0_y[k] - frc * t11_0_y[k] + f2t;

                    t_0_yz[k] = frb * t10_0_z[k] - frc * t11_0_z[k];

                    // leading z component

                    t_0_zz[k] = pbz[k] * t10_0_z[k] - pcz[k] * t11_0_z[k] + f2t;
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForDS(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {2, 0, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (2|0)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 2, 0);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PA)

                auto pax = paDistances.data(3 * idx);

                auto pay = paDistances.data(3 * idx + 1);

                auto paz = paDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (S|A(0)|S)^(m) integrals

                auto t20_0_0 = primBuffer.data(t20off + idx);

                // set up pointers to (S|A(0)|S)^(m+1) integrals

                auto t21_0_0 = primBuffer.data(t21off + idx);

                // set up pointers to (P|A(0)|S)^(m) integrals

                auto t10_x_0 = primBuffer.data(t10off + 3 * idx);

                auto t10_y_0 = primBuffer.data(t10off + 3 * idx + 1);

                auto t10_z_0 = primBuffer.data(t10off + 3 * idx + 2);

                // set up pointers to (P|A(0)|S)^(m+1) integrals

                auto t11_x_0 = primBuffer.data(t11off + 3 * idx);

                auto t11_y_0 = primBuffer.data(t11off + 3 * idx + 1);

                auto t11_z_0 = primBuffer.data(t11off + 3 * idx + 2);

                // set up pointers to (D|A(0)|S)^(m) integrals

                auto t_xx_0 = primBuffer.data(toff + 6 * idx);

                auto t_xy_0 = primBuffer.data(toff + 6 * idx + 1);

                auto t_xz_0 = primBuffer.data(toff + 6 * idx + 2);

                auto t_yy_0 = primBuffer.data(toff + 6 * idx + 3);

                auto t_yz_0 = primBuffer.data(toff + 6 * idx + 4);

                auto t_zz_0 = primBuffer.data(toff + 6 * idx + 5);

                #pragma omp simd aligned(fx, pax, pay, paz, t20_0_0, t21_0_0, t10_x_0,\
                                         t10_y_0, t10_z_0, t11_x_0, t11_y_0, t11_z_0,\
                                         t_xx_0, t_xy_0, t_xz_0, t_yy_0, t_yz_0,\
                                         t_zz_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k] * (t20_0_0[k] - t21_0_0[k]);

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xx_0[k] = fra * t10_x_0[k] - frc * t11_x_0[k] + f2t;

                    t_xy_0[k] = fra * t10_y_0[k] - frc * t11_y_0[k];

                    t_xz_0[k] = fra * t10_z_0[k] - frc * t11_z_0[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yy_0[k] = fra * t10_y_0[k] - frc * t11_y_0[k] + f2t;

                    t_yz_0[k] = fra * t10_z_0[k] - frc * t11_z_0[k];

                    // leading z component

                    t_zz_0[k] = paz[k] * t10_z_0[k] - pcz[k] * t11_z_0[k] + f2t;
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForPD(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {1, 2, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (1|2)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 2);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PA)

                auto pax = paDistances.data(3 * idx);

                auto pay = paDistances.data(3 * idx + 1);

                auto paz = paDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m) integrals

                auto tk0_0_x = primBuffer.data(tk0off + 3 * idx);

                auto tk0_0_y = primBuffer.data(tk0off + 3 * idx + 1);

                auto tk0_0_z = primBuffer.data(tk0off + 3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m+1) integrals

                auto tk1_0_x = primBuffer.data(tk1off + 3 * idx);

                auto tk1_0_y = primBuffer.data(tk1off + 3 * idx + 1);

                auto tk1_0_z = primBuffer.data(tk1off + 3 * idx + 2);

                // set up pointers to (S|A(0)|D)^(m) integrals

                auto t10_0_xx = primBuffer.data(t10off + 6 * idx);

                auto t10_0_xy = primBuffer.data(t10off + 6 * idx + 1);

                auto t10_0_xz = primBuffer.data(t10off + 6 * idx + 2);

                auto t10_0_yy = primBuffer.data(t10off + 6 * idx + 3);

                auto t10_0_yz = primBuffer.data(t10off + 6 * idx + 4);

                auto t10_0_zz = primBuffer.data(t10off + 6 * idx + 5);

                // set up pointers to (S|A(0)|D)^(m+1) integrals

                auto t11_0_xx = primBuffer.data(t11off + 6 * idx);

                auto t11_0_xy = primBuffer.data(t11off + 6 * idx + 1);

                auto t11_0_xz = primBuffer.data(t11off + 6 * idx + 2);

                auto t11_0_yy = primBuffer.data(t11off + 6 * idx + 3);

                auto t11_0_yz = primBuffer.data(t11off + 6 * idx + 4);

                auto t11_0_zz = primBuffer.data(t11off + 6 * idx + 5);

                // set up pointers to (P|A(0)|D)^(m) integrals

                auto t_x_xx = primBuffer.data(toff + 18 * idx);

                auto t_x_xy = primBuffer.data(toff + 18 * idx + 1);

                auto t_x_xz = primBuffer.data(toff + 18 * idx + 2);

                auto t_x_yy = primBuffer.data(toff + 18 * idx + 3);

                auto t_x_yz = primBuffer.data(toff + 18 * idx + 4);

                auto t_x_zz = primBuffer.data(toff + 18 * idx + 5);

                auto t_y_xx = primBuffer.data(toff + 18 * idx + 6);

                auto t_y_xy = primBuffer.data(toff + 18 * idx + 7);

                auto t_y_xz = primBuffer.data(toff + 18 * idx + 8);

                auto t_y_yy = primBuffer.data(toff + 18 * idx + 9);

                auto t_y_yz = primBuffer.data(toff + 18 * idx + 10);

                auto t_y_zz = primBuffer.data(toff + 18 * idx + 11);

                auto t_z_xx = primBuffer.data(toff + 18 * idx + 12);

                auto t_z_xy = primBuffer.data(toff + 18 * idx + 13);

                auto t_z_xz = primBuffer.data(toff + 18 * idx + 14);

                auto t_z_yy = primBuffer.data(toff + 18 * idx + 15);

                auto t_z_yz = primBuffer.data(toff + 18 * idx + 16);

                auto t_z_zz = primBuffer.data(toff + 18 * idx + 17);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_0_x, tk0_0_y, tk0_0_z,\
                                         tk1_0_x, tk1_0_y, tk1_0_z, t10_0_xx, t10_0_xy,\
                                         t10_0_xz, t10_0_yy, t10_0_yz, t10_0_zz,\
                                         t11_0_xx, t11_0_xy, t11_0_xz, t11_0_yy,\
                                         t11_0_yz, t11_0_zz, t_x_xx, t_x_xy, t_x_xz,\
                                         t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy,\
                                         t_y_xz, t_y_yy, t_y_yz, t_y_zz, t_z_xx,\
                                         t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_x_xx[k] = fra * t10_0_xx[k] - frc * t11_0_xx[k] + f2t * (2.0 * tk0_0_x[k] - 2.0 * tk1_0_x[k]);

                    t_x_xy[k] = fra * t10_0_xy[k] - frc * t11_0_xy[k] + f2t * (tk0_0_y[k] - tk1_0_y[k]);

                    t_x_xz[k] = fra * t10_0_xz[k] - frc * t11_0_xz[k] + f2t * (tk0_0_z[k] - tk1_0_z[k]);

                    t_x_yy[k] = fra * t10_0_yy[k] - frc * t11_0_yy[k];

                    t_x_yz[k] = fra * t10_0_yz[k] - frc * t11_0_yz[k];

                    t_x_zz[k] = fra * t10_0_zz[k] - frc * t11_0_zz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_y_xx[k] = fra * t10_0_xx[k]- frc * t11_0_xx[k];

                    t_y_xy[k] = fra * t10_0_xy[k]- frc * t11_0_xy[k] + f2t * (tk0_0_x[k] - tk1_0_x[k]);

                    t_y_xz[k] = fra * t10_0_xz[k]- frc * t11_0_xz[k];

                    t_y_yy[k] = fra * t10_0_yy[k]- frc * t11_0_yy[k] + f2t * (2.0 * tk0_0_y[k] - 2.0 * tk1_0_y[k]);

                    t_y_yz[k] = fra * t10_0_yz[k]- frc * t11_0_yz[k] + f2t * (tk0_0_z[k] - tk1_0_z[k]);

                    t_y_zz[k] = fra * t10_0_zz[k]- frc * t11_0_zz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_z_xx[k] = fra * t10_0_xx[k] - frc * t11_0_xx[k];

                    t_z_xy[k] = fra * t10_0_xy[k] - frc * t11_0_xy[k];

                    t_z_xz[k] = fra * t10_0_xz[k] - frc * t11_0_xz[k] + f2t * (tk0_0_x[k] - tk1_0_x[k]);

                    t_z_yy[k] = fra * t10_0_yy[k] - frc * t11_0_yy[k];

                    t_z_yz[k] = fra * t10_0_yz[k] - frc * t11_0_yz[k] + f2t * (tk0_0_y[k] - tk1_0_y[k]);

                    t_z_zz[k] = fra * t10_0_zz[k] - frc * t11_0_zz[k] + f2t * (2.0 * tk0_0_z[k] - 2.0 * tk1_0_z[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForDP(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {2, 1, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (2|1)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 2, 1);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PA)

                auto pax = paDistances.data(3 * idx);

                auto pay = paDistances.data(3 * idx + 1);

                auto paz = paDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (P|A(0)|S)^(m) integrals

                auto tk0_x_0 = primBuffer.data(tk0off + 3 * idx);

                auto tk0_y_0 = primBuffer.data(tk0off + 3 * idx + 1);

                auto tk0_z_0 = primBuffer.data(tk0off + 3 * idx + 2);

                // set up pointers to (P|A(0)|S)^(m+1) integrals

                auto tk1_x_0 = primBuffer.data(tk1off + 3 * idx);

                auto tk1_y_0 = primBuffer.data(tk1off + 3 * idx + 1);

                auto tk1_z_0 = primBuffer.data(tk1off + 3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m) integrals

                auto t20_0_x = primBuffer.data(t20off + 3 * idx);

                auto t20_0_y = primBuffer.data(t20off + 3 * idx + 1);

                auto t20_0_z = primBuffer.data(t20off + 3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m+1) integrals

                auto t21_0_x = primBuffer.data(t21off + 3 * idx);

                auto t21_0_y = primBuffer.data(t21off + 3 * idx + 1);

                auto t21_0_z = primBuffer.data(t21off + 3 * idx + 2);

                // set up pointers to (P|A(0)|P)^(m) integrals

                auto t10_x_x = primBuffer.data(t10off + 9 * idx);

                auto t10_x_y = primBuffer.data(t10off + 9 * idx + 1);

                auto t10_x_z = primBuffer.data(t10off + 9 * idx + 2);

                auto t10_y_x = primBuffer.data(t10off + 9 * idx + 3);

                auto t10_y_y = primBuffer.data(t10off + 9 * idx + 4);

                auto t10_y_z = primBuffer.data(t10off + 9 * idx + 5);

                auto t10_z_x = primBuffer.data(t10off + 9 * idx + 6);

                auto t10_z_y = primBuffer.data(t10off + 9 * idx + 7);

                auto t10_z_z = primBuffer.data(t10off + 9 * idx + 8);

                // set up pointers to (P|A(0)|P)^(m+1) integrals

                auto t11_x_x = primBuffer.data(t11off + 9 * idx);

                auto t11_x_y = primBuffer.data(t11off + 9 * idx + 1);

                auto t11_x_z = primBuffer.data(t11off + 9 * idx + 2);

                auto t11_y_x = primBuffer.data(t11off + 9 * idx + 3);

                auto t11_y_y = primBuffer.data(t11off + 9 * idx + 4);

                auto t11_y_z = primBuffer.data(t11off + 9 * idx + 5);

                auto t11_z_x = primBuffer.data(t11off + 9 * idx + 6);

                auto t11_z_y = primBuffer.data(t11off + 9 * idx + 7);

                auto t11_z_z = primBuffer.data(t11off + 9 * idx + 8);

                // set up pointers to (D|A(0)|P)^(m) integrals

                auto t_xx_x = primBuffer.data(toff + 18 * idx);

                auto t_xx_y = primBuffer.data(toff + 18 * idx + 1);

                auto t_xx_z = primBuffer.data(toff + 18 * idx + 2);

                auto t_xy_x = primBuffer.data(toff + 18 * idx + 3);

                auto t_xy_y = primBuffer.data(toff + 18 * idx + 4);

                auto t_xy_z = primBuffer.data(toff + 18 * idx + 5);

                auto t_xz_x = primBuffer.data(toff + 18 * idx + 6);

                auto t_xz_y = primBuffer.data(toff + 18 * idx + 7);

                auto t_xz_z = primBuffer.data(toff + 18 * idx + 8);

                auto t_yy_x = primBuffer.data(toff + 18 * idx + 9);

                auto t_yy_y = primBuffer.data(toff + 18 * idx + 10);

                auto t_yy_z = primBuffer.data(toff + 18 * idx + 11);

                auto t_yz_x = primBuffer.data(toff + 18 * idx + 12);

                auto t_yz_y = primBuffer.data(toff + 18 * idx + 13);

                auto t_yz_z = primBuffer.data(toff + 18 * idx + 14);

                auto t_zz_x = primBuffer.data(toff + 18 * idx + 15);

                auto t_zz_y = primBuffer.data(toff + 18 * idx + 16);

                auto t_zz_z = primBuffer.data(toff + 18 * idx + 17);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_x_0, tk0_y_0, tk0_z_0,\
                                         tk1_x_0, tk1_y_0, tk1_z_0, t20_0_x, t20_0_y,\
                                         t20_0_z, t21_0_x, t21_0_y, t21_0_z, t10_x_x,\
                                         t10_x_y, t10_x_z, t10_y_x, t10_y_y, t10_y_z,\
                                         t10_z_x, t10_z_y, t10_z_z, t11_x_x, t11_x_y,\
                                         t11_x_z, t11_y_x, t11_y_y, t11_y_z, t11_z_x,\
                                         t11_z_y, t11_z_z, t_xx_x, t_xx_y, t_xx_z,\
                                         t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y,\
                                         t_xz_z, t_yy_x, t_yy_y, t_yy_z, t_yz_x,\
                                         t_yz_y, t_yz_z, t_zz_x, t_zz_y, t_zz_z: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xx_x[k] = fra * t10_x_x[k] - frc * t11_x_x[k] + f2t * (t20_0_x[k] - t21_0_x[k] + tk0_x_0[k] - tk1_x_0[k]);

                    t_xx_y[k] = fra * t10_x_y[k] - frc * t11_x_y[k] + f2t * (t20_0_y[k] - t21_0_y[k]);

                    t_xx_z[k] = fra * t10_x_z[k] - frc * t11_x_z[k] + f2t * (t20_0_z[k] - t21_0_z[k]);

                    t_xy_x[k] = fra * t10_y_x[k] - frc * t11_y_x[k] + f2t * (tk0_y_0[k] - tk1_y_0[k]);

                    t_xy_y[k] = fra * t10_y_y[k] - frc * t11_y_y[k];

                    t_xy_z[k] = fra * t10_y_z[k] - frc * t11_y_z[k];

                    t_xz_x[k] = fra * t10_z_x[k] - frc * t11_z_x[k] + f2t * (tk0_z_0[k] - tk1_z_0[k]);

                    t_xz_y[k] = fra * t10_z_y[k] - frc * t11_z_y[k];

                    t_xz_z[k] = fra * t10_z_z[k] - frc * t11_z_z[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yy_x[k] = fra * t10_y_x[k]- frc * t11_y_x[k] + f2t * (t20_0_x[k] - t21_0_x[k]);

                    t_yy_y[k] = fra * t10_y_y[k]- frc * t11_y_y[k] + f2t * (t20_0_y[k] - t21_0_y[k] + tk0_y_0[k] - tk1_y_0[k]);

                    t_yy_z[k] = fra * t10_y_z[k]- frc * t11_y_z[k] + f2t * (t20_0_z[k] - t21_0_z[k]);

                    t_yz_x[k] = fra * t10_z_x[k]- frc * t11_z_x[k];

                    t_yz_y[k] = fra * t10_z_y[k]- frc * t11_z_y[k] + f2t * (tk0_z_0[k] - tk1_z_0[k]);

                    t_yz_z[k] = fra * t10_z_z[k]- frc * t11_z_z[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zz_x[k] = fra * t10_z_x[k] - frc * t11_z_x[k] + f2t * (t20_0_x[k] - t21_0_x[k]);

                    t_zz_y[k] = fra * t10_z_y[k] - frc * t11_z_y[k] + f2t * (t20_0_y[k] - t21_0_y[k]);

                    t_zz_z[k] = fra * t10_z_z[k] - frc * t11_z_z[k] + f2t * (t20_0_z[k] - t21_0_z[k] + tk0_z_0[k] - tk1_z_0[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForDD(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {2, 2, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (2|2)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 2, 2);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PA)

                auto pax = paDistances.data(3 * idx);

                auto pay = paDistances.data(3 * idx + 1);

                auto paz = paDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (P|A(0)|P)^(m) integrals

                auto tk0_x_x = primBuffer.data(tk0off + 9 * idx);

                auto tk0_x_y = primBuffer.data(tk0off + 9 * idx + 1);

                auto tk0_x_z = primBuffer.data(tk0off + 9 * idx + 2);

                auto tk0_y_x = primBuffer.data(tk0off + 9 * idx + 3);

                auto tk0_y_y = primBuffer.data(tk0off + 9 * idx + 4);

                auto tk0_y_z = primBuffer.data(tk0off + 9 * idx + 5);

                auto tk0_z_x = primBuffer.data(tk0off + 9 * idx + 6);

                auto tk0_z_y = primBuffer.data(tk0off + 9 * idx + 7);

                auto tk0_z_z = primBuffer.data(tk0off + 9 * idx + 8);

                // set up pointers to (P|A(0)|P)^(m+1) integrals

                auto tk1_x_x = primBuffer.data(tk1off + 9 * idx);

                auto tk1_x_y = primBuffer.data(tk1off + 9 * idx + 1);

                auto tk1_x_z = primBuffer.data(tk1off + 9 * idx + 2);

                auto tk1_y_x = primBuffer.data(tk1off + 9 * idx + 3);

                auto tk1_y_y = primBuffer.data(tk1off + 9 * idx + 4);

                auto tk1_y_z = primBuffer.data(tk1off + 9 * idx + 5);

                auto tk1_z_x = primBuffer.data(tk1off + 9 * idx + 6);

                auto tk1_z_y = primBuffer.data(tk1off + 9 * idx + 7);

                auto tk1_z_z = primBuffer.data(tk1off + 9 * idx + 8);

                // set up pointers to (S|A(0)|D)^(m) integrals

                auto t20_0_xx = primBuffer.data(t20off + 6 * idx);

                auto t20_0_xy = primBuffer.data(t20off + 6 * idx + 1);

                auto t20_0_xz = primBuffer.data(t20off + 6 * idx + 2);

                auto t20_0_yy = primBuffer.data(t20off + 6 * idx + 3);

                auto t20_0_yz = primBuffer.data(t20off + 6 * idx + 4);

                auto t20_0_zz = primBuffer.data(t20off + 6 * idx + 5);

                // set up pointers to (S|A(0)|D)^(m+1) integrals

                auto t21_0_xx = primBuffer.data(t21off + 6 * idx);

                auto t21_0_xy = primBuffer.data(t21off + 6 * idx + 1);

                auto t21_0_xz = primBuffer.data(t21off + 6 * idx + 2);

                auto t21_0_yy = primBuffer.data(t21off + 6 * idx + 3);

                auto t21_0_yz = primBuffer.data(t21off + 6 * idx + 4);

                auto t21_0_zz = primBuffer.data(t21off + 6 * idx + 5);

                // set up pointers to (P|A(0)|D)^(m) integrals

                auto t10_x_xx = primBuffer.data(t10off + 18 * idx);

                auto t10_x_xy = primBuffer.data(t10off + 18 * idx + 1);

                auto t10_x_xz = primBuffer.data(t10off + 18 * idx + 2);

                auto t10_x_yy = primBuffer.data(t10off + 18 * idx + 3);

                auto t10_x_yz = primBuffer.data(t10off + 18 * idx + 4);

                auto t10_x_zz = primBuffer.data(t10off + 18 * idx + 5);

                auto t10_y_xx = primBuffer.data(t10off + 18 * idx + 6);

                auto t10_y_xy = primBuffer.data(t10off + 18 * idx + 7);

                auto t10_y_xz = primBuffer.data(t10off + 18 * idx + 8);

                auto t10_y_yy = primBuffer.data(t10off + 18 * idx + 9);

                auto t10_y_yz = primBuffer.data(t10off + 18 * idx + 10);

                auto t10_y_zz = primBuffer.data(t10off + 18 * idx + 11);

                auto t10_z_xx = primBuffer.data(t10off + 18 * idx + 12);

                auto t10_z_xy = primBuffer.data(t10off + 18 * idx + 13);

                auto t10_z_xz = primBuffer.data(t10off + 18 * idx + 14);

                auto t10_z_yy = primBuffer.data(t10off + 18 * idx + 15);

                auto t10_z_yz = primBuffer.data(t10off + 18 * idx + 16);

                auto t10_z_zz = primBuffer.data(t10off + 18 * idx + 17);

                // set up pointers to (P|A(0)|D)^(m+1) integrals

                auto t11_x_xx = primBuffer.data(t11off + 18 * idx);

                auto t11_x_xy = primBuffer.data(t11off + 18 * idx + 1);

                auto t11_x_xz = primBuffer.data(t11off + 18 * idx + 2);

                auto t11_x_yy = primBuffer.data(t11off + 18 * idx + 3);

                auto t11_x_yz = primBuffer.data(t11off + 18 * idx + 4);

                auto t11_x_zz = primBuffer.data(t11off + 18 * idx + 5);

                auto t11_y_xx = primBuffer.data(t11off + 18 * idx + 6);

                auto t11_y_xy = primBuffer.data(t11off + 18 * idx + 7);

                auto t11_y_xz = primBuffer.data(t11off + 18 * idx + 8);

                auto t11_y_yy = primBuffer.data(t11off + 18 * idx + 9);

                auto t11_y_yz = primBuffer.data(t11off + 18 * idx + 10);

                auto t11_y_zz = primBuffer.data(t11off + 18 * idx + 11);

                auto t11_z_xx = primBuffer.data(t11off + 18 * idx + 12);

                auto t11_z_xy = primBuffer.data(t11off + 18 * idx + 13);

                auto t11_z_xz = primBuffer.data(t11off + 18 * idx + 14);

                auto t11_z_yy = primBuffer.data(t11off + 18 * idx + 15);

                auto t11_z_yz = primBuffer.data(t11off + 18 * idx + 16);

                auto t11_z_zz = primBuffer.data(t11off + 18 * idx + 17);

                // set up pointers to (D|A(0)|D)^(m) integrals

                auto t_xx_xx = primBuffer.data(toff + 36 * idx);

                auto t_xx_xy = primBuffer.data(toff + 36 * idx + 1);

                auto t_xx_xz = primBuffer.data(toff + 36 * idx + 2);

                auto t_xx_yy = primBuffer.data(toff + 36 * idx + 3);

                auto t_xx_yz = primBuffer.data(toff + 36 * idx + 4);

                auto t_xx_zz = primBuffer.data(toff + 36 * idx + 5);

                auto t_xy_xx = primBuffer.data(toff + 36 * idx + 6);

                auto t_xy_xy = primBuffer.data(toff + 36 * idx + 7);

                auto t_xy_xz = primBuffer.data(toff + 36 * idx + 8);

                auto t_xy_yy = primBuffer.data(toff + 36 * idx + 9);

                auto t_xy_yz = primBuffer.data(toff + 36 * idx + 10);

                auto t_xy_zz = primBuffer.data(toff + 36 * idx + 11);

                auto t_xz_xx = primBuffer.data(toff + 36 * idx + 12);

                auto t_xz_xy = primBuffer.data(toff + 36 * idx + 13);

                auto t_xz_xz = primBuffer.data(toff + 36 * idx + 14);

                auto t_xz_yy = primBuffer.data(toff + 36 * idx + 15);

                auto t_xz_yz = primBuffer.data(toff + 36 * idx + 16);

                auto t_xz_zz = primBuffer.data(toff + 36 * idx + 17);

                auto t_yy_xx = primBuffer.data(toff + 36 * idx + 18);

                auto t_yy_xy = primBuffer.data(toff + 36 * idx + 19);

                auto t_yy_xz = primBuffer.data(toff + 36 * idx + 20);

                auto t_yy_yy = primBuffer.data(toff + 36 * idx + 21);

                auto t_yy_yz = primBuffer.data(toff + 36 * idx + 22);

                auto t_yy_zz = primBuffer.data(toff + 36 * idx + 23);

                auto t_yz_xx = primBuffer.data(toff + 36 * idx + 24);

                auto t_yz_xy = primBuffer.data(toff + 36 * idx + 25);

                auto t_yz_xz = primBuffer.data(toff + 36 * idx + 26);

                auto t_yz_yy = primBuffer.data(toff + 36 * idx + 27);

                auto t_yz_yz = primBuffer.data(toff + 36 * idx + 28);

                auto t_yz_zz = primBuffer.data(toff + 36 * idx + 29);

                auto t_zz_xx = primBuffer.data(toff + 36 * idx + 30);

                auto t_zz_xy = primBuffer.data(toff + 36 * idx + 31);

                auto t_zz_xz = primBuffer.data(toff + 36 * idx + 32);

                auto t_zz_yy = primBuffer.data(toff + 36 * idx + 33);

                auto t_zz_yz = primBuffer.data(toff + 36 * idx + 34);

                auto t_zz_zz = primBuffer.data(toff + 36 * idx + 35);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_x_x, tk0_x_y, tk0_x_z,\
                                         tk0_y_x, tk0_y_y, tk0_y_z, tk0_z_x, tk0_z_y,\
                                         tk0_z_z, tk1_x_x, tk1_x_y, tk1_x_z, tk1_y_x,\
                                         tk1_y_y, tk1_y_z, tk1_z_x, tk1_z_y, tk1_z_z,\
                                         t20_0_xx, t20_0_xy, t20_0_xz, t20_0_yy,\
                                         t20_0_yz, t20_0_zz, t21_0_xx, t21_0_xy,\
                                         t21_0_xz, t21_0_yy, t21_0_yz, t21_0_zz,\
                                         t10_x_xx, t10_x_xy, t10_x_xz, t10_x_yy,\
                                         t10_x_yz, t10_x_zz, t10_y_xx, t10_y_xy,\
                                         t10_y_xz, t10_y_yy, t10_y_yz, t10_y_zz,\
                                         t10_z_xx, t10_z_xy, t10_z_xz, t10_z_yy,\
                                         t10_z_yz, t10_z_zz, t11_x_xx, t11_x_xy,\
                                         t11_x_xz, t11_x_yy, t11_x_yz, t11_x_zz,\
                                         t11_y_xx, t11_y_xy, t11_y_xz, t11_y_yy,\
                                         t11_y_yz, t11_y_zz, t11_z_xx, t11_z_xy,\
                                         t11_z_xz, t11_z_yy, t11_z_yz, t11_z_zz,\
                                         t_xx_xx, t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz,\
                                         t_xx_zz, t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy,\
                                         t_xy_yz, t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz,\
                                         t_xz_yy, t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy,\
                                         t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx,\
                                         t_yz_xy, t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz,\
                                         t_zz_xx, t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz,\
                                         t_zz_zz: VLX_ALIGN)
                 for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xx_xx[k] = fra * t10_x_xx[k] - frc * t11_x_xx[k] + f2t * (t20_0_xx[k] - t21_0_xx[k] + 2.0 * tk0_x_x[k] - 2.0 * tk1_x_x[k]);

                    t_xx_xy[k] = fra * t10_x_xy[k] - frc * t11_x_xy[k] + f2t * (t20_0_xy[k] - t21_0_xy[k] + tk0_x_y[k] - tk1_x_y[k]);

                    t_xx_xz[k] = fra * t10_x_xz[k] - frc * t11_x_xz[k] + f2t * (t20_0_xz[k] - t21_0_xz[k] + tk0_x_z[k] - tk1_x_z[k]);

                    t_xx_yy[k] = fra * t10_x_yy[k] - frc * t11_x_yy[k] + f2t * (t20_0_yy[k] - t21_0_yy[k]);

                    t_xx_yz[k] = fra * t10_x_yz[k] - frc * t11_x_yz[k] + f2t * (t20_0_yz[k] - t21_0_yz[k]);

                    t_xx_zz[k] = fra * t10_x_zz[k] - frc * t11_x_zz[k] + f2t * (t20_0_zz[k] - t21_0_zz[k]);

                    t_xy_xx[k] = fra * t10_y_xx[k] - frc * t11_y_xx[k] + f2t * (2.0 * tk0_y_x[k] - 2.0 * tk1_y_x[k]);

                    t_xy_xy[k] = fra * t10_y_xy[k] - frc * t11_y_xy[k] + f2t * (tk0_y_y[k] - tk1_y_y[k]);

                    t_xy_xz[k] = fra * t10_y_xz[k] - frc * t11_y_xz[k] + f2t * (tk0_y_z[k] - tk1_y_z[k]);

                    t_xy_yy[k] = fra * t10_y_yy[k] - frc * t11_y_yy[k];

                    t_xy_yz[k] = fra * t10_y_yz[k] - frc * t11_y_yz[k];

                    t_xy_zz[k] = fra * t10_y_zz[k] - frc * t11_y_zz[k];

                    t_xz_xx[k] = fra * t10_z_xx[k] - frc * t11_z_xx[k] + f2t * (2.0 * tk0_z_x[k] - 2.0 * tk1_z_x[k]);

                    t_xz_xy[k] = fra * t10_z_xy[k] - frc * t11_z_xy[k] + f2t * (tk0_z_y[k] - tk1_z_y[k]);

                    t_xz_xz[k] = fra * t10_z_xz[k] - frc * t11_z_xz[k] + f2t * (tk0_z_z[k] - tk1_z_z[k]);

                    t_xz_yy[k] = fra * t10_z_yy[k] - frc * t11_z_yy[k];

                    t_xz_yz[k] = fra * t10_z_yz[k] - frc * t11_z_yz[k];

                    t_xz_zz[k] = fra * t10_z_zz[k] - frc * t11_z_zz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yy_xx[k] = fra * t10_y_xx[k] - frc * t11_y_xx[k] + f2t * (t20_0_xx[k] - t21_0_xx[k]);

                    t_yy_xy[k] = fra * t10_y_xy[k] - frc * t11_y_xy[k] + f2t * (t20_0_xy[k] - t21_0_xy[k] + tk0_y_x[k] - tk1_y_x[k]);

                    t_yy_xz[k] = fra * t10_y_xz[k] - frc * t11_y_xz[k] + f2t * (t20_0_xz[k] - t21_0_xz[k]);

                    t_yy_yy[k] = fra * t10_y_yy[k] - frc * t11_y_yy[k] + f2t * (t20_0_yy[k] - t21_0_yy[k] + 2.0 * tk0_y_y[k] - 2.0 * tk1_y_y[k]);

                    t_yy_yz[k] = fra * t10_y_yz[k] - frc * t11_y_yz[k] + f2t * (t20_0_yz[k] - t21_0_yz[k] + tk0_y_z[k] - tk1_y_z[k]);

                    t_yy_zz[k] = fra * t10_y_zz[k] - frc * t11_y_zz[k] + f2t * (t20_0_zz[k] - t21_0_zz[k]);

                    t_yz_xx[k] = fra * t10_z_xx[k] - frc * t11_z_xx[k];

                    t_yz_xy[k] = fra * t10_z_xy[k] - frc * t11_z_xy[k] + f2t * (tk0_z_x[k] - tk1_z_x[k]);

                    t_yz_xz[k] = fra * t10_z_xz[k] - frc * t11_z_xz[k];

                    t_yz_yy[k] = fra * t10_z_yy[k] - frc * t11_z_yy[k] + f2t * (2.0 * tk0_z_y[k] - 2.0 * tk1_z_y[k]);

                    t_yz_yz[k] = fra * t10_z_yz[k] - frc * t11_z_yz[k] + f2t * (tk0_z_z[k] - tk1_z_z[k]);

                    t_yz_zz[k] = fra * t10_z_zz[k] - frc * t11_z_zz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zz_xx[k] = fra * t10_z_xx[k] - frc * t11_z_xx[k] + f2t * (t20_0_xx[k] - t21_0_xx[k]);

                    t_zz_xy[k] = fra * t10_z_xy[k] - frc * t11_z_xy[k] + f2t * (t20_0_xy[k] - t21_0_xy[k]);

                    t_zz_xz[k] = fra * t10_z_xz[k] - frc * t11_z_xz[k] + f2t * (t20_0_xz[k] - t21_0_xz[k] + tk0_z_x[k] - tk1_z_x[k]);

                    t_zz_yy[k] = fra * t10_z_yy[k] - frc * t11_z_yy[k] + f2t * (t20_0_yy[k] - t21_0_yy[k]);

                    t_zz_yz[k] = fra * t10_z_yz[k] - frc * t11_z_yz[k] + f2t * (t20_0_yz[k] - t21_0_yz[k] + tk0_z_y[k] - tk1_z_y[k]);

                    t_zz_zz[k] = fra * t10_z_zz[k] - frc * t11_z_zz[k] + f2t * (t20_0_zz[k] - t21_0_zz[k] + 2.0 * tk0_z_z[k] - 2.0 * tk1_z_z[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForSF(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  pbDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {0, 3, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (0|3)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 0, 3);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 1, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PB)

                auto pbx = pbDistances.data(3 * idx);

                auto pby = pbDistances.data(3 * idx + 1);

                auto pbz = pbDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m) integrals

                auto t20_0_x = primBuffer.data(t20off + 3 * idx);

                auto t20_0_y = primBuffer.data(t20off + 3 * idx + 1);

                auto t20_0_z = primBuffer.data(t20off + 3 * idx + 2);

                // set up pointers to (S|A(0)|P)^(m+1) integrals

                auto t21_0_x = primBuffer.data(t21off + 3 * idx);

                auto t21_0_y = primBuffer.data(t21off + 3 * idx + 1);

                auto t21_0_z = primBuffer.data(t21off + 3 * idx + 2);

                // set up pointers to (S|A(0)|D)^(m) integrals

                auto t10_0_xx = primBuffer.data(t10off + 6 * idx);

                auto t10_0_xy = primBuffer.data(t10off + 6 * idx + 1);

                auto t10_0_xz = primBuffer.data(t10off + 6 * idx + 2);

                auto t10_0_yy = primBuffer.data(t10off + 6 * idx + 3);

                auto t10_0_yz = primBuffer.data(t10off + 6 * idx + 4);

                auto t10_0_zz = primBuffer.data(t10off + 6 * idx + 5);

                // set up pointers to (S|A(0)|D)^(m+1) integrals

                auto t11_0_xx = primBuffer.data(t11off + 6 * idx);

                auto t11_0_xy = primBuffer.data(t11off + 6 * idx + 1);

                auto t11_0_xz = primBuffer.data(t11off + 6 * idx + 2);

                auto t11_0_yy = primBuffer.data(t11off + 6 * idx + 3);

                auto t11_0_yz = primBuffer.data(t11off + 6 * idx + 4);

                auto t11_0_zz = primBuffer.data(t11off + 6 * idx + 5);

                // set up pointers to (S|A(0)|F)^(m) integrals

                auto t_0_xxx = primBuffer.data(toff + 10 * idx);

                auto t_0_xxy = primBuffer.data(toff + 10 * idx + 1);

                auto t_0_xxz = primBuffer.data(toff + 10 * idx + 2);

                auto t_0_xyy = primBuffer.data(toff + 10 * idx + 3);

                auto t_0_xyz = primBuffer.data(toff + 10 * idx + 4);

                auto t_0_xzz = primBuffer.data(toff + 10 * idx + 5);

                auto t_0_yyy = primBuffer.data(toff + 10 * idx + 6);

                auto t_0_yyz = primBuffer.data(toff + 10 * idx + 7);

                auto t_0_yzz = primBuffer.data(toff + 10 * idx + 8);

                auto t_0_zzz = primBuffer.data(toff + 10 * idx + 9);

                #pragma omp simd aligned(fx, pbx, pby, pbz, t20_0_x, t20_0_y, t20_0_z,\
                                         t21_0_x, t21_0_y, t21_0_z, t10_0_xx, t10_0_xy,\
                                         t10_0_xz, t10_0_yy, t10_0_yz, t10_0_zz,\
                                         t11_0_xx, t11_0_xy, t11_0_xz, t11_0_yy,\
                                         t11_0_yz, t11_0_zz, t_0_xxx, t_0_xxy,\
                                         t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy,\
                                         t_0_yyz, t_0_yzz, t_0_zzz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pbx[k];

                    double frc = pcx[k];

                    t_0_xxx[k] = fra * t10_0_xx[k] - frc * t11_0_xx[k] + f2t * (2.0 * t20_0_x[k] - 2.0 * t21_0_x[k]);

                    t_0_xxy[k] = fra * t10_0_xy[k] - frc * t11_0_xy[k] + f2t * (t20_0_y[k] - t21_0_y[k]);

                    t_0_xxz[k] = fra * t10_0_xz[k] - frc * t11_0_xz[k] + f2t * (t20_0_z[k] - t21_0_z[k]);

                    t_0_xyy[k] = fra * t10_0_yy[k] - frc * t11_0_yy[k];

                    t_0_xyz[k] = fra * t10_0_yz[k] - frc * t11_0_yz[k];

                    t_0_xzz[k] = fra * t10_0_zz[k] - frc * t11_0_zz[k];

                    // leading y component

                    fra = pby[k];

                    frc = pcy[k];

                    t_0_yyy[k] = fra * t10_0_yy[k] - frc * t10_0_yy[k] + f2t * (2.0 * t20_0_y[k] - 2.0 * t21_0_y[k]);

                    t_0_yyz[k] = fra * t10_0_yz[k] - frc * t10_0_yz[k] + f2t * (t20_0_z[k] - t21_0_z[k]);

                    t_0_yzz[k] = fra * t10_0_zz[k] - frc * t10_0_zz[k];

                    // leading z component

                    t_0_zzz[k] = pbz[k] * t10_0_zz[k] - pcz[k] * t10_0_zz[k] + f2t * (2.0 * t20_0_z[k] - 2.0 * t21_0_z[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForFS(      CMemBlock2D<double>&  primBuffer,
                              const CVecThreeIndexes&     recPattern,
                              const std::vector<int32_t>& recIndexes,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
    {
        // skip integrals if not included in recursion pattern

        if (!genfunc::isInVector(recPattern, {3, 0, 0})) return;

        if (iContrGto  == 0) printf(" * VRR: (3|0)\n");

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 3, 0);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 0, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 0, i + 1});

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to distances R(PA)

                auto pax = paDistances.data(3 * idx);

                auto pay = paDistances.data(3 * idx + 1);

                auto paz = paDistances.data(3 * idx + 2);

                // set up pointers to distances R(PC)

                auto pcx = pcDistances.data(3 * idx);

                auto pcy = pcDistances.data(3 * idx + 1);

                auto pcz = pcDistances.data(3 * idx + 2);

                // set up pointers to (P|A(0)|S)^(m) integrals

                auto t20_x_0 = primBuffer.data(t20off + 3 * idx);

                auto t20_y_0 = primBuffer.data(t20off + 3 * idx + 1);

                auto t20_z_0 = primBuffer.data(t20off + 3 * idx + 2);

                // set up pointers to (P|A(0)|S)^(m+1) integrals

                auto t21_x_0 = primBuffer.data(t21off + 3 * idx);

                auto t21_y_0 = primBuffer.data(t21off + 3 * idx + 1);

                auto t21_z_0 = primBuffer.data(t21off + 3 * idx + 2);

                // set up pointers to (D|A(0)|S)^(m) integrals

                auto t10_xx_0 = primBuffer.data(t10off + 6 * idx);

                auto t10_xy_0 = primBuffer.data(t10off + 6 * idx + 1);

                auto t10_xz_0 = primBuffer.data(t10off + 6 * idx + 2);

                auto t10_yy_0 = primBuffer.data(t10off + 6 * idx + 3);

                auto t10_yz_0 = primBuffer.data(t10off + 6 * idx + 4);

                auto t10_zz_0 = primBuffer.data(t10off + 6 * idx + 5);

                // set up pointers to (D|A(0)|S)^(m+1) integrals

                auto t11_xx_0 = primBuffer.data(t11off + 6 * idx);

                auto t11_xy_0 = primBuffer.data(t11off + 6 * idx + 1);

                auto t11_xz_0 = primBuffer.data(t11off + 6 * idx + 2);

                auto t11_yy_0 = primBuffer.data(t11off + 6 * idx + 3);

                auto t11_yz_0 = primBuffer.data(t11off + 6 * idx + 4);

                auto t11_zz_0 = primBuffer.data(t11off + 6 * idx + 5);

                // set up pointers to (F|A(0)|S)^(m) integrals

                auto t_xxx_0 = primBuffer.data(toff + 10 * idx);

                auto t_xxy_0 = primBuffer.data(toff + 10 * idx + 1);

                auto t_xxz_0 = primBuffer.data(toff + 10 * idx + 2);

                auto t_xyy_0 = primBuffer.data(toff + 10 * idx + 3);

                auto t_xyz_0 = primBuffer.data(toff + 10 * idx + 4);

                auto t_xzz_0 = primBuffer.data(toff + 10 * idx + 5);

                auto t_yyy_0 = primBuffer.data(toff + 10 * idx + 6);

                auto t_yyz_0 = primBuffer.data(toff + 10 * idx + 7);

                auto t_yzz_0 = primBuffer.data(toff + 10 * idx + 8);

                auto t_zzz_0 = primBuffer.data(toff + 10 * idx + 9);

                #pragma omp simd aligned(fx, pax, pay, paz, t20_x_0, t20_y_0, t20_z_0,\
                                         t21_x_0, t21_y_0, t21_z_0, t10_xx_0, t10_xy_0,\
                                         t10_xz_0, t10_yy_0, t10_yz_0, t10_zz_0,\
                                         t11_xx_0, t11_xy_0, t11_xz_0, t11_yy_0,\
                                         t11_yz_0, t11_zz_0, t_xxx_0, t_xxy_0,\
                                         t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, t_yyy_0,\
                                         t_yyz_0, t_yzz_0, t_zzz_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxx_0[k] = fra * t10_xx_0[k] - frc * t11_xx_0[k] + f2t * (2.0 * t20_x_0[k] - 2.0 * t21_x_0[k]);

                    t_xxy_0[k] = fra * t10_xy_0[k] - frc * t11_xy_0[k] + f2t * (t20_y_0[k] - t21_y_0[k]);

                    t_xxz_0[k] = fra * t10_xz_0[k] - frc * t11_xz_0[k] + f2t * (t20_z_0[k] - t21_z_0[k]);

                    t_xyy_0[k] = fra * t10_yy_0[k] - frc * t11_yy_0[k];

                    t_xyz_0[k] = fra * t10_yz_0[k] - frc * t11_yz_0[k];

                    t_xzz_0[k] = fra * t10_zz_0[k] - frc * t11_zz_0[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyy_0[k] = fra * t10_yy_0[k] - frc * t10_yy_0[k] + f2t * (2.0 * t20_y_0[k] - 2.0 * t21_y_0[k]);

                    t_yyz_0[k] = fra * t10_yz_0[k] - frc * t10_yz_0[k] + f2t * (t20_z_0[k] - t21_z_0[k]);

                    t_yzz_0[k] = fra * t10_zz_0[k] - frc * t10_zz_0[k];

                    // leading z component

                    t_zzz_0[k] = paz[k] * t10_zz_0[k] - pcz[k] * t10_zz_0[k] + f2t * (2.0 * t20_z_0[k] - 2.0 * t21_z_0[k]);
                }

                idx++;
            }
        }
    }
    
} // npotrecfunc namespace
