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
            
            #pragma omp simd aligned(fargs, fg, pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_0_zz, pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_zz_0, pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz,\
                                         pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_yz_y, t_yz_z, t_zz_x, t_zz_y, t_zz_z,\
                                         pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_zz_zz, pcx, pcy, pcz: VLX_ALIGN)
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
                                         t_0_yyz, t_0_yzz, t_0_zzz, pcx, pcy, pcz: VLX_ALIGN)
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

                    t_0_yyy[k] = fra * t10_0_yy[k] - frc * t11_0_yy[k] + f2t * (2.0 * t20_0_y[k] - 2.0 * t21_0_y[k]);

                    t_0_yyz[k] = fra * t10_0_yz[k] - frc * t11_0_yz[k] + f2t * (t20_0_z[k] - t21_0_z[k]);

                    t_0_yzz[k] = fra * t10_0_zz[k] - frc * t11_0_zz[k];

                    // leading z component

                    t_0_zzz[k] = pbz[k] * t10_0_zz[k] - pcz[k] * t11_0_zz[k] + f2t * (2.0 * t20_0_z[k] - 2.0 * t21_0_z[k]);
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
                                         t_yyz_0, t_yzz_0, t_zzz_0, pcx, pcy, pcz: VLX_ALIGN)
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

                    t_yyy_0[k] = fra * t10_yy_0[k] - frc * t11_yy_0[k] + f2t * (2.0 * t20_y_0[k] - 2.0 * t21_y_0[k]);

                    t_yyz_0[k] = fra * t10_yz_0[k] - frc * t11_yz_0[k] + f2t * (t20_z_0[k] - t21_z_0[k]);

                    t_yzz_0[k] = fra * t10_zz_0[k] - frc * t11_zz_0[k];

                    // leading z component

                    t_zzz_0[k] = paz[k] * t10_zz_0[k] - pcz[k] * t11_zz_0[k] + f2t * (2.0 * t20_z_0[k] - 2.0 * t21_z_0[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForPF(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {1, 3, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 3);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i + 1});

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

                // set up pointers to (S|A(0)|D)^(m) integrals

                auto tk0_0_xx = primBuffer.data(tk0off + 6 * idx);

                auto tk0_0_xy = primBuffer.data(tk0off + 6 * idx + 1);

                auto tk0_0_xz = primBuffer.data(tk0off + 6 * idx + 2);

                auto tk0_0_yy = primBuffer.data(tk0off + 6 * idx + 3);

                auto tk0_0_yz = primBuffer.data(tk0off + 6 * idx + 4);

                auto tk0_0_zz = primBuffer.data(tk0off + 6 * idx + 5);

                // set up pointers to (S|A(0)|D)^(m+1) integrals

                auto tk1_0_xx = primBuffer.data(tk1off + 6 * idx);

                auto tk1_0_xy = primBuffer.data(tk1off + 6 * idx + 1);

                auto tk1_0_xz = primBuffer.data(tk1off + 6 * idx + 2);

                auto tk1_0_yy = primBuffer.data(tk1off + 6 * idx + 3);

                auto tk1_0_yz = primBuffer.data(tk1off + 6 * idx + 4);

                auto tk1_0_zz = primBuffer.data(tk1off + 6 * idx + 5);

                // set up pointers to (S|A(0)|F)^(m) integrals

                auto t10_0_xxx = primBuffer.data(t10off + 10 * idx);

                auto t10_0_xxy = primBuffer.data(t10off + 10 * idx + 1);

                auto t10_0_xxz = primBuffer.data(t10off + 10 * idx + 2);

                auto t10_0_xyy = primBuffer.data(t10off + 10 * idx + 3);

                auto t10_0_xyz = primBuffer.data(t10off + 10 * idx + 4);

                auto t10_0_xzz = primBuffer.data(t10off + 10 * idx + 5);

                auto t10_0_yyy = primBuffer.data(t10off + 10 * idx + 6);

                auto t10_0_yyz = primBuffer.data(t10off + 10 * idx + 7);

                auto t10_0_yzz = primBuffer.data(t10off + 10 * idx + 8);

                auto t10_0_zzz = primBuffer.data(t10off + 10 * idx + 9);

                // set up pointers to (S|A(0)|F)^(m+1) integrals

                auto t11_0_xxx = primBuffer.data(t11off + 10 * idx);

                auto t11_0_xxy = primBuffer.data(t11off + 10 * idx + 1);

                auto t11_0_xxz = primBuffer.data(t11off + 10 * idx + 2);

                auto t11_0_xyy = primBuffer.data(t11off + 10 * idx + 3);

                auto t11_0_xyz = primBuffer.data(t11off + 10 * idx + 4);

                auto t11_0_xzz = primBuffer.data(t11off + 10 * idx + 5);

                auto t11_0_yyy = primBuffer.data(t11off + 10 * idx + 6);

                auto t11_0_yyz = primBuffer.data(t11off + 10 * idx + 7);

                auto t11_0_yzz = primBuffer.data(t11off + 10 * idx + 8);

                auto t11_0_zzz = primBuffer.data(t11off + 10 * idx + 9);

                // set up pointers to (P|A(0)|F)^(m) integrals

                auto t_x_xxx = primBuffer.data(toff + 30 * idx);

                auto t_x_xxy = primBuffer.data(toff + 30 * idx + 1);

                auto t_x_xxz = primBuffer.data(toff + 30 * idx + 2);

                auto t_x_xyy = primBuffer.data(toff + 30 * idx + 3);

                auto t_x_xyz = primBuffer.data(toff + 30 * idx + 4);

                auto t_x_xzz = primBuffer.data(toff + 30 * idx + 5);

                auto t_x_yyy = primBuffer.data(toff + 30 * idx + 6);

                auto t_x_yyz = primBuffer.data(toff + 30 * idx + 7);

                auto t_x_yzz = primBuffer.data(toff + 30 * idx + 8);

                auto t_x_zzz = primBuffer.data(toff + 30 * idx + 9);

                auto t_y_xxx = primBuffer.data(toff + 30 * idx + 10);

                auto t_y_xxy = primBuffer.data(toff + 30 * idx + 11);

                auto t_y_xxz = primBuffer.data(toff + 30 * idx + 12);

                auto t_y_xyy = primBuffer.data(toff + 30 * idx + 13);

                auto t_y_xyz = primBuffer.data(toff + 30 * idx + 14);

                auto t_y_xzz = primBuffer.data(toff + 30 * idx + 15);

                auto t_y_yyy = primBuffer.data(toff + 30 * idx + 16);

                auto t_y_yyz = primBuffer.data(toff + 30 * idx + 17);

                auto t_y_yzz = primBuffer.data(toff + 30 * idx + 18);

                auto t_y_zzz = primBuffer.data(toff + 30 * idx + 19);

                auto t_z_xxx = primBuffer.data(toff + 30 * idx + 20);

                auto t_z_xxy = primBuffer.data(toff + 30 * idx + 21);

                auto t_z_xxz = primBuffer.data(toff + 30 * idx + 22);

                auto t_z_xyy = primBuffer.data(toff + 30 * idx + 23);

                auto t_z_xyz = primBuffer.data(toff + 30 * idx + 24);

                auto t_z_xzz = primBuffer.data(toff + 30 * idx + 25);

                auto t_z_yyy = primBuffer.data(toff + 30 * idx + 26);

                auto t_z_yyz = primBuffer.data(toff + 30 * idx + 27);

                auto t_z_yzz = primBuffer.data(toff + 30 * idx + 28);

                auto t_z_zzz = primBuffer.data(toff + 30 * idx + 29);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_0_xx, tk0_0_xy,\
                                         tk0_0_xz, tk0_0_yy, tk0_0_yz, tk0_0_zz,\
                                         tk1_0_xx, tk1_0_xy, tk1_0_xz, tk1_0_yy,\
                                         tk1_0_yz, tk1_0_zz, t10_0_xxx, t10_0_xxy,\
                                         t10_0_xxz, t10_0_xyy, t10_0_xyz, t10_0_xzz,\
                                         t10_0_yyy, t10_0_yyz, t10_0_yzz, t10_0_zzz,\
                                         t11_0_xxx, t11_0_xxy, t11_0_xxz, t11_0_xyy,\
                                         t11_0_xyz, t11_0_xzz, t11_0_yyy, t11_0_yyz,\
                                         t11_0_yzz, t11_0_zzz, t_x_xxx, t_x_xxy,\
                                         t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy,\
                                         t_x_yyz, t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy,\
                                         t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy,\
                                         t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy,\
                                         t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy,\
                                         t_z_yyz, t_z_yzz, t_z_zzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_x_xxx[k] = fra * t10_0_xxx[k] - frc * t11_0_xxx[k] + f2t * (3.0 * tk0_0_xx[k] - 3.0 * tk1_0_xx[k]);

                    t_x_xxy[k] = fra * t10_0_xxy[k] - frc * t11_0_xxy[k] + f2t * (2.0 * tk0_0_xy[k] - 2.0 * tk1_0_xy[k]);

                    t_x_xxz[k] = fra * t10_0_xxz[k] - frc * t11_0_xxz[k] + f2t * (2.0 * tk0_0_xz[k] - 2.0 * tk1_0_xz[k]);

                    t_x_xyy[k] = fra * t10_0_xyy[k] - frc * t11_0_xyy[k] + f2t * (tk0_0_yy[k] - tk1_0_yy[k]);

                    t_x_xyz[k] = fra * t10_0_xyz[k] - frc * t11_0_xyz[k] + f2t * (tk0_0_yz[k] - tk1_0_yz[k]);

                    t_x_xzz[k] = fra * t10_0_xzz[k] - frc * t11_0_xzz[k] + f2t * (tk0_0_zz[k] - tk1_0_zz[k]);

                    t_x_yyy[k] = fra * t10_0_yyy[k] - frc * t11_0_yyy[k];

                    t_x_yyz[k] = fra * t10_0_yyz[k] - frc * t11_0_yyz[k];

                    t_x_yzz[k] = fra * t10_0_yzz[k] - frc * t11_0_yzz[k];

                    t_x_zzz[k] = fra * t10_0_zzz[k] - frc * t11_0_zzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_y_xxx[k] = fra * t10_0_xxx[k] - frc * t11_0_xxx[k];

                    t_y_xxy[k] = fra * t10_0_xxy[k] - frc * t11_0_xxy[k] + f2t * (tk0_0_xx[k] - tk1_0_xx[k]);

                    t_y_xxz[k] = fra * t10_0_xxz[k] - frc * t11_0_xxz[k];

                    t_y_xyy[k] = fra * t10_0_xyy[k] - frc * t11_0_xyy[k] + f2t * (2.0 * tk0_0_xy[k] - 2.0 * tk1_0_xy[k]);

                    t_y_xyz[k] = fra * t10_0_xyz[k] - frc * t11_0_xyz[k] + f2t * (tk0_0_xz[k] - tk1_0_xz[k]);

                    t_y_xzz[k] = fra * t10_0_xzz[k] - frc * t11_0_xzz[k];

                    t_y_yyy[k] = fra * t10_0_yyy[k] - frc * t11_0_yyy[k] + f2t * (3.0 * tk0_0_yy[k] - 3.0 * tk1_0_yy[k]);

                    t_y_yyz[k] = fra * t10_0_yyz[k] - frc * t11_0_yyz[k] + f2t * (2.0 * tk0_0_yz[k] - 2.0 * tk1_0_yz[k]);

                    t_y_yzz[k] = fra * t10_0_yzz[k] - frc * t11_0_yzz[k] + f2t * (tk0_0_zz[k] - tk1_0_zz[k]);

                    t_y_zzz[k] = fra * t10_0_zzz[k] - frc * t11_0_zzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_z_xxx[k] = fra * t10_0_xxx[k] - frc * t11_0_xxx[k];

                    t_z_xxy[k] = fra * t10_0_xxy[k] - frc * t11_0_xxy[k];

                    t_z_xxz[k] = fra * t10_0_xxz[k] - frc * t11_0_xxz[k] + f2t * (tk0_0_xx[k] - tk1_0_xx[k]);

                    t_z_xyy[k] = fra * t10_0_xyy[k] - frc * t11_0_xyy[k];

                    t_z_xyz[k] = fra * t10_0_xyz[k] - frc * t11_0_xyz[k] + f2t * (tk0_0_xy[k] - tk1_0_xy[k]);

                    t_z_xzz[k] = fra * t10_0_xzz[k] - frc * t11_0_xzz[k] + f2t * (2.0 * tk0_0_xz[k] - 2.0 * tk1_0_xz[k]);

                    t_z_yyy[k] = fra * t10_0_yyy[k] - frc * t11_0_yyy[k];

                    t_z_yyz[k] = fra * t10_0_yyz[k] - frc * t11_0_yyz[k] + f2t * (tk0_0_yy[k] - tk1_0_yy[k]);

                    t_z_yzz[k] = fra * t10_0_yzz[k] - frc * t11_0_yzz[k] + f2t * (2.0 * tk0_0_yz[k] - 2.0 * tk1_0_yz[k]);

                    t_z_zzz[k] = fra * t10_0_zzz[k] - frc * t11_0_zzz[k] + f2t * (3.0 * tk0_0_zz[k] - 3.0 * tk1_0_zz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForFP(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {3, 1, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 3, 1);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 1, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 1, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i + 1});

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

                // set up pointers to (D|A(0)|S)^(m) integrals

                auto tk0_xx_0 = primBuffer.data(tk0off + 6 * idx);

                auto tk0_xy_0 = primBuffer.data(tk0off + 6 * idx + 1);

                auto tk0_xz_0 = primBuffer.data(tk0off + 6 * idx + 2);

                auto tk0_yy_0 = primBuffer.data(tk0off + 6 * idx + 3);

                auto tk0_yz_0 = primBuffer.data(tk0off + 6 * idx + 4);

                auto tk0_zz_0 = primBuffer.data(tk0off + 6 * idx + 5);

                // set up pointers to (D|A(0)|S)^(m+1) integrals

                auto tk1_xx_0 = primBuffer.data(tk1off + 6 * idx);

                auto tk1_xy_0 = primBuffer.data(tk1off + 6 * idx + 1);

                auto tk1_xz_0 = primBuffer.data(tk1off + 6 * idx + 2);

                auto tk1_yy_0 = primBuffer.data(tk1off + 6 * idx + 3);

                auto tk1_yz_0 = primBuffer.data(tk1off + 6 * idx + 4);

                auto tk1_zz_0 = primBuffer.data(tk1off + 6 * idx + 5);

                // set up pointers to (P|A(0)|P)^(m) integrals

                auto t20_x_x = primBuffer.data(t20off + 9 * idx);

                auto t20_x_y = primBuffer.data(t20off + 9 * idx + 1);

                auto t20_x_z = primBuffer.data(t20off + 9 * idx + 2);

                auto t20_y_x = primBuffer.data(t20off + 9 * idx + 3);

                auto t20_y_y = primBuffer.data(t20off + 9 * idx + 4);

                auto t20_y_z = primBuffer.data(t20off + 9 * idx + 5);

                auto t20_z_x = primBuffer.data(t20off + 9 * idx + 6);

                auto t20_z_y = primBuffer.data(t20off + 9 * idx + 7);

                auto t20_z_z = primBuffer.data(t20off + 9 * idx + 8);

                // set up pointers to (P|A(0)|P)^(m+1) integrals

                auto t21_x_x = primBuffer.data(t21off + 9 * idx);

                auto t21_x_y = primBuffer.data(t21off + 9 * idx + 1);

                auto t21_x_z = primBuffer.data(t21off + 9 * idx + 2);

                auto t21_y_x = primBuffer.data(t21off + 9 * idx + 3);

                auto t21_y_y = primBuffer.data(t21off + 9 * idx + 4);

                auto t21_y_z = primBuffer.data(t21off + 9 * idx + 5);

                auto t21_z_x = primBuffer.data(t21off + 9 * idx + 6);

                auto t21_z_y = primBuffer.data(t21off + 9 * idx + 7);

                auto t21_z_z = primBuffer.data(t21off + 9 * idx + 8);

                // set up pointers to (D|A(0)|P)^(m) integrals

                auto t10_xx_x = primBuffer.data(t10off + 18 * idx);

                auto t10_xx_y = primBuffer.data(t10off + 18 * idx + 1);

                auto t10_xx_z = primBuffer.data(t10off + 18 * idx + 2);

                auto t10_xy_x = primBuffer.data(t10off + 18 * idx + 3);

                auto t10_xy_y = primBuffer.data(t10off + 18 * idx + 4);

                auto t10_xy_z = primBuffer.data(t10off + 18 * idx + 5);

                auto t10_xz_x = primBuffer.data(t10off + 18 * idx + 6);

                auto t10_xz_y = primBuffer.data(t10off + 18 * idx + 7);

                auto t10_xz_z = primBuffer.data(t10off + 18 * idx + 8);

                auto t10_yy_x = primBuffer.data(t10off + 18 * idx + 9);

                auto t10_yy_y = primBuffer.data(t10off + 18 * idx + 10);

                auto t10_yy_z = primBuffer.data(t10off + 18 * idx + 11);

                auto t10_yz_x = primBuffer.data(t10off + 18 * idx + 12);

                auto t10_yz_y = primBuffer.data(t10off + 18 * idx + 13);

                auto t10_yz_z = primBuffer.data(t10off + 18 * idx + 14);

                auto t10_zz_x = primBuffer.data(t10off + 18 * idx + 15);

                auto t10_zz_y = primBuffer.data(t10off + 18 * idx + 16);

                auto t10_zz_z = primBuffer.data(t10off + 18 * idx + 17);

                // set up pointers to (D|A(0)|P)^(m+1) integrals

                auto t11_xx_x = primBuffer.data(t11off + 18 * idx);

                auto t11_xx_y = primBuffer.data(t11off + 18 * idx + 1);

                auto t11_xx_z = primBuffer.data(t11off + 18 * idx + 2);

                auto t11_xy_x = primBuffer.data(t11off + 18 * idx + 3);

                auto t11_xy_y = primBuffer.data(t11off + 18 * idx + 4);

                auto t11_xy_z = primBuffer.data(t11off + 18 * idx + 5);

                auto t11_xz_x = primBuffer.data(t11off + 18 * idx + 6);

                auto t11_xz_y = primBuffer.data(t11off + 18 * idx + 7);

                auto t11_xz_z = primBuffer.data(t11off + 18 * idx + 8);

                auto t11_yy_x = primBuffer.data(t11off + 18 * idx + 9);

                auto t11_yy_y = primBuffer.data(t11off + 18 * idx + 10);

                auto t11_yy_z = primBuffer.data(t11off + 18 * idx + 11);

                auto t11_yz_x = primBuffer.data(t11off + 18 * idx + 12);

                auto t11_yz_y = primBuffer.data(t11off + 18 * idx + 13);

                auto t11_yz_z = primBuffer.data(t11off + 18 * idx + 14);

                auto t11_zz_x = primBuffer.data(t11off + 18 * idx + 15);

                auto t11_zz_y = primBuffer.data(t11off + 18 * idx + 16);

                auto t11_zz_z = primBuffer.data(t11off + 18 * idx + 17);

                // set up pointers to (F|A(0)|P)^(m) integrals

                auto t_xxx_x = primBuffer.data(toff + 30 * idx);

                auto t_xxx_y = primBuffer.data(toff + 30 * idx + 1);

                auto t_xxx_z = primBuffer.data(toff + 30 * idx + 2);

                auto t_xxy_x = primBuffer.data(toff + 30 * idx + 3);

                auto t_xxy_y = primBuffer.data(toff + 30 * idx + 4);

                auto t_xxy_z = primBuffer.data(toff + 30 * idx + 5);

                auto t_xxz_x = primBuffer.data(toff + 30 * idx + 6);

                auto t_xxz_y = primBuffer.data(toff + 30 * idx + 7);

                auto t_xxz_z = primBuffer.data(toff + 30 * idx + 8);

                auto t_xyy_x = primBuffer.data(toff + 30 * idx + 9);

                auto t_xyy_y = primBuffer.data(toff + 30 * idx + 10);

                auto t_xyy_z = primBuffer.data(toff + 30 * idx + 11);

                auto t_xyz_x = primBuffer.data(toff + 30 * idx + 12);

                auto t_xyz_y = primBuffer.data(toff + 30 * idx + 13);

                auto t_xyz_z = primBuffer.data(toff + 30 * idx + 14);

                auto t_xzz_x = primBuffer.data(toff + 30 * idx + 15);

                auto t_xzz_y = primBuffer.data(toff + 30 * idx + 16);

                auto t_xzz_z = primBuffer.data(toff + 30 * idx + 17);

                auto t_yyy_x = primBuffer.data(toff + 30 * idx + 18);

                auto t_yyy_y = primBuffer.data(toff + 30 * idx + 19);

                auto t_yyy_z = primBuffer.data(toff + 30 * idx + 20);

                auto t_yyz_x = primBuffer.data(toff + 30 * idx + 21);

                auto t_yyz_y = primBuffer.data(toff + 30 * idx + 22);

                auto t_yyz_z = primBuffer.data(toff + 30 * idx + 23);

                auto t_yzz_x = primBuffer.data(toff + 30 * idx + 24);

                auto t_yzz_y = primBuffer.data(toff + 30 * idx + 25);

                auto t_yzz_z = primBuffer.data(toff + 30 * idx + 26);

                auto t_zzz_x = primBuffer.data(toff + 30 * idx + 27);

                auto t_zzz_y = primBuffer.data(toff + 30 * idx + 28);

                auto t_zzz_z = primBuffer.data(toff + 30 * idx + 29);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xx_0, tk0_xy_0,\
                                         tk0_xz_0, tk0_yy_0, tk0_yz_0, tk0_zz_0,\
                                         tk1_xx_0, tk1_xy_0, tk1_xz_0, tk1_yy_0,\
                                         tk1_yz_0, tk1_zz_0, t20_x_x, t20_x_y,\
                                         t20_x_z, t20_y_x, t20_y_y, t20_y_z, t20_z_x,\
                                         t20_z_y, t20_z_z, t21_x_x, t21_x_y, t21_x_z,\
                                         t21_y_x, t21_y_y, t21_y_z, t21_z_x, t21_z_y,\
                                         t21_z_z, t10_xx_x, t10_xx_y, t10_xx_z,\
                                         t10_xy_x, t10_xy_y, t10_xy_z, t10_xz_x,\
                                         t10_xz_y, t10_xz_z, t10_yy_x, t10_yy_y,\
                                         t10_yy_z, t10_yz_x, t10_yz_y, t10_yz_z,\
                                         t10_zz_x, t10_zz_y, t10_zz_z, t11_xx_x,\
                                         t11_xx_y, t11_xx_z, t11_xy_x, t11_xy_y,\
                                         t11_xy_z, t11_xz_x, t11_xz_y, t11_xz_z,\
                                         t11_yy_x, t11_yy_y, t11_yy_z, t11_yz_x,\
                                         t11_yz_y, t11_yz_z, t11_zz_x, t11_zz_y,\
                                         t11_zz_z, t_xxx_x, t_xxx_y, t_xxx_z, t_xxy_x,\
                                         t_xxy_y, t_xxy_z, t_xxz_x, t_xxz_y, t_xxz_z,\
                                         t_xyy_x, t_xyy_y, t_xyy_z, t_xyz_x, t_xyz_y,\
                                         t_xyz_z, t_xzz_x, t_xzz_y, t_xzz_z, t_yyy_x,\
                                         t_yyy_y, t_yyy_z, t_yyz_x, t_yyz_y, t_yyz_z,\
                                         t_yzz_x, t_yzz_y, t_yzz_z, t_zzz_x, t_zzz_y,\
                                         t_zzz_z, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxx_x[k] = fra * t10_xx_x[k] - frc * t11_xx_x[k] + f2t * (2.0 * t20_x_x[k] - 2.0 * t21_x_x[k] + tk0_xx_0[k] - tk1_xx_0[k]);

                    t_xxx_y[k] = fra * t10_xx_y[k] - frc * t11_xx_y[k] + f2t * (2.0 * t20_x_y[k] - 2.0 * t21_x_y[k]);

                    t_xxx_z[k] = fra * t10_xx_z[k] - frc * t11_xx_z[k] + f2t * (2.0 * t20_x_z[k] - 2.0 * t21_x_z[k]);

                    t_xxy_x[k] = fra * t10_xy_x[k] - frc * t11_xy_x[k] + f2t * (t20_y_x[k] - t21_y_x[k] + tk0_xy_0[k] - tk1_xy_0[k]);

                    t_xxy_y[k] = fra * t10_xy_y[k] - frc * t11_xy_y[k] + f2t * (t20_y_y[k] - t21_y_y[k]);

                    t_xxy_z[k] = fra * t10_xy_z[k] - frc * t11_xy_z[k] + f2t * (t20_y_z[k] - t21_y_z[k]);

                    t_xxz_x[k] = fra * t10_xz_x[k] - frc * t11_xz_x[k] + f2t * (t20_z_x[k] - t21_z_x[k] + tk0_xz_0[k] - tk1_xz_0[k]);

                    t_xxz_y[k] = fra * t10_xz_y[k] - frc * t11_xz_y[k] + f2t * (t20_z_y[k] - t21_z_y[k]);

                    t_xxz_z[k] = fra * t10_xz_z[k] - frc * t11_xz_z[k] + f2t * (t20_z_z[k] - t21_z_z[k]);

                    t_xyy_x[k] = fra * t10_yy_x[k] - frc * t11_yy_x[k] + f2t * (tk0_yy_0[k] - tk1_yy_0[k]);

                    t_xyy_y[k] = fra * t10_yy_y[k] - frc * t11_yy_y[k];

                    t_xyy_z[k] = fra * t10_yy_z[k] - frc * t11_yy_z[k];

                    t_xyz_x[k] = fra * t10_yz_x[k] - frc * t11_yz_x[k] + f2t * (tk0_yz_0[k] - tk1_yz_0[k]);

                    t_xyz_y[k] = fra * t10_yz_y[k] - frc * t11_yz_y[k];

                    t_xyz_z[k] = fra * t10_yz_z[k] - frc * t11_yz_z[k];

                    t_xzz_x[k] = fra * t10_zz_x[k] - frc * t11_zz_x[k] + f2t * (tk0_zz_0[k] - tk1_zz_0[k]);

                    t_xzz_y[k] = fra * t10_zz_y[k] - frc * t11_zz_y[k];

                    t_xzz_z[k] = fra * t10_zz_z[k] - frc * t11_zz_z[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyy_x[k] = fra * t10_yy_x[k] - frc * t11_yy_x[k] + f2t * (2.0 * t20_y_x[k] - 2.0 * t21_y_x[k]);

                    t_yyy_y[k] = fra * t10_yy_y[k] - frc * t11_yy_y[k] + f2t * (2.0 * t20_y_y[k] - 2.0 * t21_y_y[k] + tk0_yy_0[k] - tk1_yy_0[k]);

                    t_yyy_z[k] = fra * t10_yy_z[k] - frc * t11_yy_z[k] + f2t * (2.0 * t20_y_z[k] - 2.0 * t21_y_z[k]);

                    t_yyz_x[k] = fra * t10_yz_x[k] - frc * t11_yz_x[k] + f2t * (t20_z_x[k] - t21_z_x[k]);

                    t_yyz_y[k] = fra * t10_yz_y[k] - frc * t11_yz_y[k] + f2t * (t20_z_y[k] - t21_z_y[k] + tk0_yz_0[k] - tk1_yz_0[k]);

                    t_yyz_z[k] = fra * t10_yz_z[k] - frc * t11_yz_z[k] + f2t * (t20_z_z[k] - t21_z_z[k]);

                    t_yzz_x[k] = fra * t10_zz_x[k] - frc * t11_zz_x[k];

                    t_yzz_y[k] = fra * t10_zz_y[k] - frc * t11_zz_y[k] + f2t * (tk0_zz_0[k] - tk1_zz_0[k]);

                    t_yzz_z[k] = fra * t10_zz_z[k] - frc * t11_zz_z[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzz_x[k] = fra * t10_zz_x[k] - frc * t11_zz_x[k] + f2t * (2.0 * t20_z_x[k] - 2.0 * t21_z_x[k]);

                    t_zzz_y[k] = fra * t10_zz_y[k] - frc * t11_zz_y[k] + f2t * (2.0 * t20_z_y[k] - 2.0 * t21_z_y[k]);

                    t_zzz_z[k] = fra * t10_zz_z[k] - frc * t11_zz_z[k] + f2t * (2.0 * t20_z_z[k] - 2.0 * t21_z_z[k] + tk0_zz_0[k] - tk1_zz_0[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForDF(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {2, 3, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 2, 3);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i + 1});

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

                // set up pointers to (P|A(0)|D)^(m) integrals

                auto tk0_x_xx = primBuffer.data(tk0off + 18 * idx);

                auto tk0_x_xy = primBuffer.data(tk0off + 18 * idx + 1);

                auto tk0_x_xz = primBuffer.data(tk0off + 18 * idx + 2);

                auto tk0_x_yy = primBuffer.data(tk0off + 18 * idx + 3);

                auto tk0_x_yz = primBuffer.data(tk0off + 18 * idx + 4);

                auto tk0_x_zz = primBuffer.data(tk0off + 18 * idx + 5);

                auto tk0_y_xx = primBuffer.data(tk0off + 18 * idx + 6);

                auto tk0_y_xy = primBuffer.data(tk0off + 18 * idx + 7);

                auto tk0_y_xz = primBuffer.data(tk0off + 18 * idx + 8);

                auto tk0_y_yy = primBuffer.data(tk0off + 18 * idx + 9);

                auto tk0_y_yz = primBuffer.data(tk0off + 18 * idx + 10);

                auto tk0_y_zz = primBuffer.data(tk0off + 18 * idx + 11);

                auto tk0_z_xx = primBuffer.data(tk0off + 18 * idx + 12);

                auto tk0_z_xy = primBuffer.data(tk0off + 18 * idx + 13);

                auto tk0_z_xz = primBuffer.data(tk0off + 18 * idx + 14);

                auto tk0_z_yy = primBuffer.data(tk0off + 18 * idx + 15);

                auto tk0_z_yz = primBuffer.data(tk0off + 18 * idx + 16);

                auto tk0_z_zz = primBuffer.data(tk0off + 18 * idx + 17);

                // set up pointers to (P|A(0)|D)^(m+1) integrals

                auto tk1_x_xx = primBuffer.data(tk1off + 18 * idx);

                auto tk1_x_xy = primBuffer.data(tk1off + 18 * idx + 1);

                auto tk1_x_xz = primBuffer.data(tk1off + 18 * idx + 2);

                auto tk1_x_yy = primBuffer.data(tk1off + 18 * idx + 3);

                auto tk1_x_yz = primBuffer.data(tk1off + 18 * idx + 4);

                auto tk1_x_zz = primBuffer.data(tk1off + 18 * idx + 5);

                auto tk1_y_xx = primBuffer.data(tk1off + 18 * idx + 6);

                auto tk1_y_xy = primBuffer.data(tk1off + 18 * idx + 7);

                auto tk1_y_xz = primBuffer.data(tk1off + 18 * idx + 8);

                auto tk1_y_yy = primBuffer.data(tk1off + 18 * idx + 9);

                auto tk1_y_yz = primBuffer.data(tk1off + 18 * idx + 10);

                auto tk1_y_zz = primBuffer.data(tk1off + 18 * idx + 11);

                auto tk1_z_xx = primBuffer.data(tk1off + 18 * idx + 12);

                auto tk1_z_xy = primBuffer.data(tk1off + 18 * idx + 13);

                auto tk1_z_xz = primBuffer.data(tk1off + 18 * idx + 14);

                auto tk1_z_yy = primBuffer.data(tk1off + 18 * idx + 15);

                auto tk1_z_yz = primBuffer.data(tk1off + 18 * idx + 16);

                auto tk1_z_zz = primBuffer.data(tk1off + 18 * idx + 17);

                // set up pointers to (S|A(0)|F)^(m) integrals

                auto t20_0_xxx = primBuffer.data(t20off + 10 * idx);

                auto t20_0_xxy = primBuffer.data(t20off + 10 * idx + 1);

                auto t20_0_xxz = primBuffer.data(t20off + 10 * idx + 2);

                auto t20_0_xyy = primBuffer.data(t20off + 10 * idx + 3);

                auto t20_0_xyz = primBuffer.data(t20off + 10 * idx + 4);

                auto t20_0_xzz = primBuffer.data(t20off + 10 * idx + 5);

                auto t20_0_yyy = primBuffer.data(t20off + 10 * idx + 6);

                auto t20_0_yyz = primBuffer.data(t20off + 10 * idx + 7);

                auto t20_0_yzz = primBuffer.data(t20off + 10 * idx + 8);

                auto t20_0_zzz = primBuffer.data(t20off + 10 * idx + 9);

                // set up pointers to (S|A(0)|F)^(m+1) integrals

                auto t21_0_xxx = primBuffer.data(t21off + 10 * idx);

                auto t21_0_xxy = primBuffer.data(t21off + 10 * idx + 1);

                auto t21_0_xxz = primBuffer.data(t21off + 10 * idx + 2);

                auto t21_0_xyy = primBuffer.data(t21off + 10 * idx + 3);

                auto t21_0_xyz = primBuffer.data(t21off + 10 * idx + 4);

                auto t21_0_xzz = primBuffer.data(t21off + 10 * idx + 5);

                auto t21_0_yyy = primBuffer.data(t21off + 10 * idx + 6);

                auto t21_0_yyz = primBuffer.data(t21off + 10 * idx + 7);

                auto t21_0_yzz = primBuffer.data(t21off + 10 * idx + 8);

                auto t21_0_zzz = primBuffer.data(t21off + 10 * idx + 9);

                // set up pointers to (P|A(0)|F)^(m) integrals

                auto t10_x_xxx = primBuffer.data(t10off + 30 * idx);

                auto t10_x_xxy = primBuffer.data(t10off + 30 * idx + 1);

                auto t10_x_xxz = primBuffer.data(t10off + 30 * idx + 2);

                auto t10_x_xyy = primBuffer.data(t10off + 30 * idx + 3);

                auto t10_x_xyz = primBuffer.data(t10off + 30 * idx + 4);

                auto t10_x_xzz = primBuffer.data(t10off + 30 * idx + 5);

                auto t10_x_yyy = primBuffer.data(t10off + 30 * idx + 6);

                auto t10_x_yyz = primBuffer.data(t10off + 30 * idx + 7);

                auto t10_x_yzz = primBuffer.data(t10off + 30 * idx + 8);

                auto t10_x_zzz = primBuffer.data(t10off + 30 * idx + 9);

                auto t10_y_xxx = primBuffer.data(t10off + 30 * idx + 10);

                auto t10_y_xxy = primBuffer.data(t10off + 30 * idx + 11);

                auto t10_y_xxz = primBuffer.data(t10off + 30 * idx + 12);

                auto t10_y_xyy = primBuffer.data(t10off + 30 * idx + 13);

                auto t10_y_xyz = primBuffer.data(t10off + 30 * idx + 14);

                auto t10_y_xzz = primBuffer.data(t10off + 30 * idx + 15);

                auto t10_y_yyy = primBuffer.data(t10off + 30 * idx + 16);

                auto t10_y_yyz = primBuffer.data(t10off + 30 * idx + 17);

                auto t10_y_yzz = primBuffer.data(t10off + 30 * idx + 18);

                auto t10_y_zzz = primBuffer.data(t10off + 30 * idx + 19);

                auto t10_z_xxx = primBuffer.data(t10off + 30 * idx + 20);

                auto t10_z_xxy = primBuffer.data(t10off + 30 * idx + 21);

                auto t10_z_xxz = primBuffer.data(t10off + 30 * idx + 22);

                auto t10_z_xyy = primBuffer.data(t10off + 30 * idx + 23);

                auto t10_z_xyz = primBuffer.data(t10off + 30 * idx + 24);

                auto t10_z_xzz = primBuffer.data(t10off + 30 * idx + 25);

                auto t10_z_yyy = primBuffer.data(t10off + 30 * idx + 26);

                auto t10_z_yyz = primBuffer.data(t10off + 30 * idx + 27);

                auto t10_z_yzz = primBuffer.data(t10off + 30 * idx + 28);

                auto t10_z_zzz = primBuffer.data(t10off + 30 * idx + 29);

                // set up pointers to (P|A(0)|F)^(m+1) integrals

                auto t11_x_xxx = primBuffer.data(t11off + 30 * idx);

                auto t11_x_xxy = primBuffer.data(t11off + 30 * idx + 1);

                auto t11_x_xxz = primBuffer.data(t11off + 30 * idx + 2);

                auto t11_x_xyy = primBuffer.data(t11off + 30 * idx + 3);

                auto t11_x_xyz = primBuffer.data(t11off + 30 * idx + 4);

                auto t11_x_xzz = primBuffer.data(t11off + 30 * idx + 5);

                auto t11_x_yyy = primBuffer.data(t11off + 30 * idx + 6);

                auto t11_x_yyz = primBuffer.data(t11off + 30 * idx + 7);

                auto t11_x_yzz = primBuffer.data(t11off + 30 * idx + 8);

                auto t11_x_zzz = primBuffer.data(t11off + 30 * idx + 9);

                auto t11_y_xxx = primBuffer.data(t11off + 30 * idx + 10);

                auto t11_y_xxy = primBuffer.data(t11off + 30 * idx + 11);

                auto t11_y_xxz = primBuffer.data(t11off + 30 * idx + 12);

                auto t11_y_xyy = primBuffer.data(t11off + 30 * idx + 13);

                auto t11_y_xyz = primBuffer.data(t11off + 30 * idx + 14);

                auto t11_y_xzz = primBuffer.data(t11off + 30 * idx + 15);

                auto t11_y_yyy = primBuffer.data(t11off + 30 * idx + 16);

                auto t11_y_yyz = primBuffer.data(t11off + 30 * idx + 17);

                auto t11_y_yzz = primBuffer.data(t11off + 30 * idx + 18);

                auto t11_y_zzz = primBuffer.data(t11off + 30 * idx + 19);

                auto t11_z_xxx = primBuffer.data(t11off + 30 * idx + 20);

                auto t11_z_xxy = primBuffer.data(t11off + 30 * idx + 21);

                auto t11_z_xxz = primBuffer.data(t11off + 30 * idx + 22);

                auto t11_z_xyy = primBuffer.data(t11off + 30 * idx + 23);

                auto t11_z_xyz = primBuffer.data(t11off + 30 * idx + 24);

                auto t11_z_xzz = primBuffer.data(t11off + 30 * idx + 25);

                auto t11_z_yyy = primBuffer.data(t11off + 30 * idx + 26);

                auto t11_z_yyz = primBuffer.data(t11off + 30 * idx + 27);

                auto t11_z_yzz = primBuffer.data(t11off + 30 * idx + 28);

                auto t11_z_zzz = primBuffer.data(t11off + 30 * idx + 29);

                // set up pointers to (D|A(0)|F)^(m) integrals

                auto t_xx_xxx = primBuffer.data(toff + 60 * idx);

                auto t_xx_xxy = primBuffer.data(toff + 60 * idx + 1);

                auto t_xx_xxz = primBuffer.data(toff + 60 * idx + 2);

                auto t_xx_xyy = primBuffer.data(toff + 60 * idx + 3);

                auto t_xx_xyz = primBuffer.data(toff + 60 * idx + 4);

                auto t_xx_xzz = primBuffer.data(toff + 60 * idx + 5);

                auto t_xx_yyy = primBuffer.data(toff + 60 * idx + 6);

                auto t_xx_yyz = primBuffer.data(toff + 60 * idx + 7);

                auto t_xx_yzz = primBuffer.data(toff + 60 * idx + 8);

                auto t_xx_zzz = primBuffer.data(toff + 60 * idx + 9);

                auto t_xy_xxx = primBuffer.data(toff + 60 * idx + 10);

                auto t_xy_xxy = primBuffer.data(toff + 60 * idx + 11);

                auto t_xy_xxz = primBuffer.data(toff + 60 * idx + 12);

                auto t_xy_xyy = primBuffer.data(toff + 60 * idx + 13);

                auto t_xy_xyz = primBuffer.data(toff + 60 * idx + 14);

                auto t_xy_xzz = primBuffer.data(toff + 60 * idx + 15);

                auto t_xy_yyy = primBuffer.data(toff + 60 * idx + 16);

                auto t_xy_yyz = primBuffer.data(toff + 60 * idx + 17);

                auto t_xy_yzz = primBuffer.data(toff + 60 * idx + 18);

                auto t_xy_zzz = primBuffer.data(toff + 60 * idx + 19);

                auto t_xz_xxx = primBuffer.data(toff + 60 * idx + 20);

                auto t_xz_xxy = primBuffer.data(toff + 60 * idx + 21);

                auto t_xz_xxz = primBuffer.data(toff + 60 * idx + 22);

                auto t_xz_xyy = primBuffer.data(toff + 60 * idx + 23);

                auto t_xz_xyz = primBuffer.data(toff + 60 * idx + 24);

                auto t_xz_xzz = primBuffer.data(toff + 60 * idx + 25);

                auto t_xz_yyy = primBuffer.data(toff + 60 * idx + 26);

                auto t_xz_yyz = primBuffer.data(toff + 60 * idx + 27);

                auto t_xz_yzz = primBuffer.data(toff + 60 * idx + 28);

                auto t_xz_zzz = primBuffer.data(toff + 60 * idx + 29);

                auto t_yy_xxx = primBuffer.data(toff + 60 * idx + 30);

                auto t_yy_xxy = primBuffer.data(toff + 60 * idx + 31);

                auto t_yy_xxz = primBuffer.data(toff + 60 * idx + 32);

                auto t_yy_xyy = primBuffer.data(toff + 60 * idx + 33);

                auto t_yy_xyz = primBuffer.data(toff + 60 * idx + 34);

                auto t_yy_xzz = primBuffer.data(toff + 60 * idx + 35);

                auto t_yy_yyy = primBuffer.data(toff + 60 * idx + 36);

                auto t_yy_yyz = primBuffer.data(toff + 60 * idx + 37);

                auto t_yy_yzz = primBuffer.data(toff + 60 * idx + 38);

                auto t_yy_zzz = primBuffer.data(toff + 60 * idx + 39);

                auto t_yz_xxx = primBuffer.data(toff + 60 * idx + 40);

                auto t_yz_xxy = primBuffer.data(toff + 60 * idx + 41);

                auto t_yz_xxz = primBuffer.data(toff + 60 * idx + 42);

                auto t_yz_xyy = primBuffer.data(toff + 60 * idx + 43);

                auto t_yz_xyz = primBuffer.data(toff + 60 * idx + 44);

                auto t_yz_xzz = primBuffer.data(toff + 60 * idx + 45);

                auto t_yz_yyy = primBuffer.data(toff + 60 * idx + 46);

                auto t_yz_yyz = primBuffer.data(toff + 60 * idx + 47);

                auto t_yz_yzz = primBuffer.data(toff + 60 * idx + 48);

                auto t_yz_zzz = primBuffer.data(toff + 60 * idx + 49);

                auto t_zz_xxx = primBuffer.data(toff + 60 * idx + 50);

                auto t_zz_xxy = primBuffer.data(toff + 60 * idx + 51);

                auto t_zz_xxz = primBuffer.data(toff + 60 * idx + 52);

                auto t_zz_xyy = primBuffer.data(toff + 60 * idx + 53);

                auto t_zz_xyz = primBuffer.data(toff + 60 * idx + 54);

                auto t_zz_xzz = primBuffer.data(toff + 60 * idx + 55);

                auto t_zz_yyy = primBuffer.data(toff + 60 * idx + 56);

                auto t_zz_yyz = primBuffer.data(toff + 60 * idx + 57);

                auto t_zz_yzz = primBuffer.data(toff + 60 * idx + 58);

                auto t_zz_zzz = primBuffer.data(toff + 60 * idx + 59);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_x_xx, tk0_x_xy,\
                                         tk0_x_xz, tk0_x_yy, tk0_x_yz, tk0_x_zz,\
                                         tk0_y_xx, tk0_y_xy, tk0_y_xz, tk0_y_yy,\
                                         tk0_y_yz, tk0_y_zz, tk0_z_xx, tk0_z_xy,\
                                         tk0_z_xz, tk0_z_yy, tk0_z_yz, tk0_z_zz,\
                                         tk1_x_xx, tk1_x_xy, tk1_x_xz, tk1_x_yy,\
                                         tk1_x_yz, tk1_x_zz, tk1_y_xx, tk1_y_xy,\
                                         tk1_y_xz, tk1_y_yy, tk1_y_yz, tk1_y_zz,\
                                         tk1_z_xx, tk1_z_xy, tk1_z_xz, tk1_z_yy,\
                                         tk1_z_yz, tk1_z_zz, t20_0_xxx, t20_0_xxy,\
                                         t20_0_xxz, t20_0_xyy, t20_0_xyz, t20_0_xzz,\
                                         t20_0_yyy, t20_0_yyz, t20_0_yzz, t20_0_zzz,\
                                         t21_0_xxx, t21_0_xxy, t21_0_xxz, t21_0_xyy,\
                                         t21_0_xyz, t21_0_xzz, t21_0_yyy, t21_0_yyz,\
                                         t21_0_yzz, t21_0_zzz, t10_x_xxx, t10_x_xxy,\
                                         t10_x_xxz, t10_x_xyy, t10_x_xyz, t10_x_xzz,\
                                         t10_x_yyy, t10_x_yyz, t10_x_yzz, t10_x_zzz,\
                                         t10_y_xxx, t10_y_xxy, t10_y_xxz, t10_y_xyy,\
                                         t10_y_xyz, t10_y_xzz, t10_y_yyy, t10_y_yyz,\
                                         t10_y_yzz, t10_y_zzz, t10_z_xxx, t10_z_xxy,\
                                         t10_z_xxz, t10_z_xyy, t10_z_xyz, t10_z_xzz,\
                                         t10_z_yyy, t10_z_yyz, t10_z_yzz, t10_z_zzz,\
                                         t11_x_xxx, t11_x_xxy, t11_x_xxz, t11_x_xyy,\
                                         t11_x_xyz, t11_x_xzz, t11_x_yyy, t11_x_yyz,\
                                         t11_x_yzz, t11_x_zzz, t11_y_xxx, t11_y_xxy,\
                                         t11_y_xxz, t11_y_xyy, t11_y_xyz, t11_y_xzz,\
                                         t11_y_yyy, t11_y_yyz, t11_y_yzz, t11_y_zzz,\
                                         t11_z_xxx, t11_z_xxy, t11_z_xxz, t11_z_xyy,\
                                         t11_z_xyz, t11_z_xzz, t11_z_yyy, t11_z_yyz,\
                                         t11_z_yzz, t11_z_zzz, t_xx_xxx, t_xx_xxy,\
                                         t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz,\
                                         t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz,\
                                         t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy,\
                                         t_xy_xyz, t_xy_xzz, t_xy_yyy, t_xy_yyz,\
                                         t_xy_yzz, t_xy_zzz, t_xz_xxx, t_xz_xxy,\
                                         t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz,\
                                         t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz,\
                                         t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy,\
                                         t_yy_xyz, t_yy_xzz, t_yy_yyy, t_yy_yyz,\
                                         t_yy_yzz, t_yy_zzz, t_yz_xxx, t_yz_xxy,\
                                         t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz,\
                                         t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz,\
                                         t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy,\
                                         t_zz_xyz, t_zz_xzz, t_zz_yyy, t_zz_yyz,\
                                         t_zz_yzz, t_zz_zzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xx_xxx[k] = fra * t10_x_xxx[k] - frc * t11_x_xxx[k] + f2t * (t20_0_xxx[k] - t21_0_xxx[k] + 3.0 * tk0_x_xx[k] - 3.0 * tk1_x_xx[k]);

                    t_xx_xxy[k] = fra * t10_x_xxy[k] - frc * t11_x_xxy[k] + f2t * (t20_0_xxy[k] - t21_0_xxy[k] + 2.0 * tk0_x_xy[k] - 2.0 * tk1_x_xy[k]);

                    t_xx_xxz[k] = fra * t10_x_xxz[k] - frc * t11_x_xxz[k] + f2t * (t20_0_xxz[k] - t21_0_xxz[k] + 2.0 * tk0_x_xz[k] - 2.0 * tk1_x_xz[k]);

                    t_xx_xyy[k] = fra * t10_x_xyy[k] - frc * t11_x_xyy[k] + f2t * (t20_0_xyy[k] - t21_0_xyy[k] + tk0_x_yy[k] - tk1_x_yy[k]);

                    t_xx_xyz[k] = fra * t10_x_xyz[k] - frc * t11_x_xyz[k] + f2t * (t20_0_xyz[k] - t21_0_xyz[k] + tk0_x_yz[k] - tk1_x_yz[k]);

                    t_xx_xzz[k] = fra * t10_x_xzz[k] - frc * t11_x_xzz[k] + f2t * (t20_0_xzz[k] - t21_0_xzz[k] + tk0_x_zz[k] - tk1_x_zz[k]);

                    t_xx_yyy[k] = fra * t10_x_yyy[k] - frc * t11_x_yyy[k] + f2t * (t20_0_yyy[k] - t21_0_yyy[k]);

                    t_xx_yyz[k] = fra * t10_x_yyz[k] - frc * t11_x_yyz[k] + f2t * (t20_0_yyz[k] - t21_0_yyz[k]);

                    t_xx_yzz[k] = fra * t10_x_yzz[k] - frc * t11_x_yzz[k] + f2t * (t20_0_yzz[k] - t21_0_yzz[k]);

                    t_xx_zzz[k] = fra * t10_x_zzz[k] - frc * t11_x_zzz[k] + f2t * (t20_0_zzz[k] - t21_0_zzz[k]);

                    t_xy_xxx[k] = fra * t10_y_xxx[k] - frc * t11_y_xxx[k] + f2t * (3.0 * tk0_y_xx[k] - 3.0 * tk1_y_xx[k]);

                    t_xy_xxy[k] = fra * t10_y_xxy[k] - frc * t11_y_xxy[k] + f2t * (2.0 * tk0_y_xy[k] - 2.0 * tk1_y_xy[k]);

                    t_xy_xxz[k] = fra * t10_y_xxz[k] - frc * t11_y_xxz[k] + f2t * (2.0 * tk0_y_xz[k] - 2.0 * tk1_y_xz[k]);

                    t_xy_xyy[k] = fra * t10_y_xyy[k] - frc * t11_y_xyy[k] + f2t * (tk0_y_yy[k] - tk1_y_yy[k]);

                    t_xy_xyz[k] = fra * t10_y_xyz[k] - frc * t11_y_xyz[k] + f2t * (tk0_y_yz[k] - tk1_y_yz[k]);

                    t_xy_xzz[k] = fra * t10_y_xzz[k] - frc * t11_y_xzz[k] + f2t * (tk0_y_zz[k] - tk1_y_zz[k]);

                    t_xy_yyy[k] = fra * t10_y_yyy[k] - frc * t11_y_yyy[k];

                    t_xy_yyz[k] = fra * t10_y_yyz[k] - frc * t11_y_yyz[k];

                    t_xy_yzz[k] = fra * t10_y_yzz[k] - frc * t11_y_yzz[k];

                    t_xy_zzz[k] = fra * t10_y_zzz[k] - frc * t11_y_zzz[k];

                    t_xz_xxx[k] = fra * t10_z_xxx[k] - frc * t11_z_xxx[k] + f2t * (3.0 * tk0_z_xx[k] - 3.0 * tk1_z_xx[k]);

                    t_xz_xxy[k] = fra * t10_z_xxy[k] - frc * t11_z_xxy[k] + f2t * (2.0 * tk0_z_xy[k] - 2.0 * tk1_z_xy[k]);

                    t_xz_xxz[k] = fra * t10_z_xxz[k] - frc * t11_z_xxz[k] + f2t * (2.0 * tk0_z_xz[k] - 2.0 * tk1_z_xz[k]);

                    t_xz_xyy[k] = fra * t10_z_xyy[k] - frc * t11_z_xyy[k] + f2t * (tk0_z_yy[k] - tk1_z_yy[k]);

                    t_xz_xyz[k] = fra * t10_z_xyz[k] - frc * t11_z_xyz[k] + f2t * (tk0_z_yz[k] - tk1_z_yz[k]);

                    t_xz_xzz[k] = fra * t10_z_xzz[k] - frc * t11_z_xzz[k] + f2t * (tk0_z_zz[k] - tk1_z_zz[k]);

                    t_xz_yyy[k] = fra * t10_z_yyy[k] - frc * t11_z_yyy[k];

                    t_xz_yyz[k] = fra * t10_z_yyz[k] - frc * t11_z_yyz[k];

                    t_xz_yzz[k] = fra * t10_z_yzz[k] - frc * t11_z_yzz[k];

                    t_xz_zzz[k] = fra * t10_z_zzz[k] - frc * t11_z_zzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yy_xxx[k] = fra * t10_y_xxx[k] - frc * t11_y_xxx[k] + f2t * (t20_0_xxx[k] - t21_0_xxx[k]);

                    t_yy_xxy[k] = fra * t10_y_xxy[k] - frc * t11_y_xxy[k] + f2t * (t20_0_xxy[k] - t21_0_xxy[k] + tk0_y_xx[k] - tk1_y_xx[k]);

                    t_yy_xxz[k] = fra * t10_y_xxz[k] - frc * t11_y_xxz[k] + f2t * (t20_0_xxz[k] - t21_0_xxz[k]);

                    t_yy_xyy[k] = fra * t10_y_xyy[k] - frc * t11_y_xyy[k] + f2t * (t20_0_xyy[k] - t21_0_xyy[k] + 2.0 * tk0_y_xy[k] - 2.0 * tk1_y_xy[k]);

                    t_yy_xyz[k] = fra * t10_y_xyz[k] - frc * t11_y_xyz[k] + f2t * (t20_0_xyz[k] - t21_0_xyz[k] + tk0_y_xz[k] - tk1_y_xz[k]);

                    t_yy_xzz[k] = fra * t10_y_xzz[k] - frc * t11_y_xzz[k] + f2t * (t20_0_xzz[k] - t21_0_xzz[k]);

                    t_yy_yyy[k] = fra * t10_y_yyy[k] - frc * t11_y_yyy[k] + f2t * (t20_0_yyy[k] - t21_0_yyy[k] + 3.0 * tk0_y_yy[k] - 3.0 * tk1_y_yy[k]);

                    t_yy_yyz[k] = fra * t10_y_yyz[k] - frc * t11_y_yyz[k] + f2t * (t20_0_yyz[k] - t21_0_yyz[k] + 2.0 * tk0_y_yz[k] - 2.0 * tk1_y_yz[k]);

                    t_yy_yzz[k] = fra * t10_y_yzz[k] - frc * t11_y_yzz[k] + f2t * (t20_0_yzz[k] - t21_0_yzz[k] + tk0_y_zz[k] - tk1_y_zz[k]);

                    t_yy_zzz[k] = fra * t10_y_zzz[k] - frc * t11_y_zzz[k] + f2t * (t20_0_zzz[k] - t21_0_zzz[k]);

                    t_yz_xxx[k] = fra * t10_z_xxx[k] - frc * t11_z_xxx[k];

                    t_yz_xxy[k] = fra * t10_z_xxy[k] - frc * t11_z_xxy[k] + f2t * (tk0_z_xx[k] - tk1_z_xx[k]);

                    t_yz_xxz[k] = fra * t10_z_xxz[k] - frc * t11_z_xxz[k];

                    t_yz_xyy[k] = fra * t10_z_xyy[k] - frc * t11_z_xyy[k] + f2t * (2.0 * tk0_z_xy[k] - 2.0 * tk1_z_xy[k]);

                    t_yz_xyz[k] = fra * t10_z_xyz[k] - frc * t11_z_xyz[k] + f2t * (tk0_z_xz[k] - tk1_z_xz[k]);

                    t_yz_xzz[k] = fra * t10_z_xzz[k] - frc * t11_z_xzz[k];

                    t_yz_yyy[k] = fra * t10_z_yyy[k] - frc * t11_z_yyy[k] + f2t * (3.0 * tk0_z_yy[k] - 3.0 * tk1_z_yy[k]);

                    t_yz_yyz[k] = fra * t10_z_yyz[k] - frc * t11_z_yyz[k] + f2t * (2.0 * tk0_z_yz[k] - 2.0 * tk1_z_yz[k]);

                    t_yz_yzz[k] = fra * t10_z_yzz[k] - frc * t11_z_yzz[k] + f2t * (tk0_z_zz[k] - tk1_z_zz[k]);

                    t_yz_zzz[k] = fra * t10_z_zzz[k] - frc * t11_z_zzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zz_xxx[k] = fra * t10_z_xxx[k] - frc * t11_z_xxx[k] + f2t * (t20_0_xxx[k] - t21_0_xxx[k]);

                    t_zz_xxy[k] = fra * t10_z_xxy[k] - frc * t11_z_xxy[k] + f2t * (t20_0_xxy[k] - t21_0_xxy[k]);

                    t_zz_xxz[k] = fra * t10_z_xxz[k] - frc * t11_z_xxz[k] + f2t * (t20_0_xxz[k] - t21_0_xxz[k] + tk0_z_xx[k] - tk1_z_xx[k]);

                    t_zz_xyy[k] = fra * t10_z_xyy[k] - frc * t11_z_xyy[k] + f2t * (t20_0_xyy[k] - t21_0_xyy[k]);

                    t_zz_xyz[k] = fra * t10_z_xyz[k] - frc * t11_z_xyz[k] + f2t * (t20_0_xyz[k] - t21_0_xyz[k] + tk0_z_xy[k] - tk1_z_xy[k]);

                    t_zz_xzz[k] = fra * t10_z_xzz[k] - frc * t11_z_xzz[k] + f2t * (t20_0_xzz[k] - t21_0_xzz[k] + 2.0 * tk0_z_xz[k] - 2.0 * tk1_z_xz[k]);

                    t_zz_yyy[k] = fra * t10_z_yyy[k] - frc * t11_z_yyy[k] + f2t * (t20_0_yyy[k] - t21_0_yyy[k]);

                    t_zz_yyz[k] = fra * t10_z_yyz[k] - frc * t11_z_yyz[k] + f2t * (t20_0_yyz[k] - t21_0_yyz[k] + tk0_z_yy[k] - tk1_z_yy[k]);

                    t_zz_yzz[k] = fra * t10_z_yzz[k] - frc * t11_z_yzz[k] + f2t * (t20_0_yzz[k] - t21_0_yzz[k] + 2.0 * tk0_z_yz[k] - 2.0 * tk1_z_yz[k]);

                    t_zz_zzz[k] = fra * t10_z_zzz[k] - frc * t11_z_zzz[k] + f2t * (t20_0_zzz[k] - t21_0_zzz[k] + 3.0 * tk0_z_zz[k] - 3.0 * tk1_z_zz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForFD(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {3, 2, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 3, 2);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 2, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 2, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i + 1});

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

                // set up pointers to (D|A(0)|P)^(m) integrals

                auto tk0_xx_x = primBuffer.data(tk0off + 18 * idx);

                auto tk0_xx_y = primBuffer.data(tk0off + 18 * idx + 1);

                auto tk0_xx_z = primBuffer.data(tk0off + 18 * idx + 2);

                auto tk0_xy_x = primBuffer.data(tk0off + 18 * idx + 3);

                auto tk0_xy_y = primBuffer.data(tk0off + 18 * idx + 4);

                auto tk0_xy_z = primBuffer.data(tk0off + 18 * idx + 5);

                auto tk0_xz_x = primBuffer.data(tk0off + 18 * idx + 6);

                auto tk0_xz_y = primBuffer.data(tk0off + 18 * idx + 7);

                auto tk0_xz_z = primBuffer.data(tk0off + 18 * idx + 8);

                auto tk0_yy_x = primBuffer.data(tk0off + 18 * idx + 9);

                auto tk0_yy_y = primBuffer.data(tk0off + 18 * idx + 10);

                auto tk0_yy_z = primBuffer.data(tk0off + 18 * idx + 11);

                auto tk0_yz_x = primBuffer.data(tk0off + 18 * idx + 12);

                auto tk0_yz_y = primBuffer.data(tk0off + 18 * idx + 13);

                auto tk0_yz_z = primBuffer.data(tk0off + 18 * idx + 14);

                auto tk0_zz_x = primBuffer.data(tk0off + 18 * idx + 15);

                auto tk0_zz_y = primBuffer.data(tk0off + 18 * idx + 16);

                auto tk0_zz_z = primBuffer.data(tk0off + 18 * idx + 17);

                // set up pointers to (D|A(0)|P)^(m+1) integrals

                auto tk1_xx_x = primBuffer.data(tk1off + 18 * idx);

                auto tk1_xx_y = primBuffer.data(tk1off + 18 * idx + 1);

                auto tk1_xx_z = primBuffer.data(tk1off + 18 * idx + 2);

                auto tk1_xy_x = primBuffer.data(tk1off + 18 * idx + 3);

                auto tk1_xy_y = primBuffer.data(tk1off + 18 * idx + 4);

                auto tk1_xy_z = primBuffer.data(tk1off + 18 * idx + 5);

                auto tk1_xz_x = primBuffer.data(tk1off + 18 * idx + 6);

                auto tk1_xz_y = primBuffer.data(tk1off + 18 * idx + 7);

                auto tk1_xz_z = primBuffer.data(tk1off + 18 * idx + 8);

                auto tk1_yy_x = primBuffer.data(tk1off + 18 * idx + 9);

                auto tk1_yy_y = primBuffer.data(tk1off + 18 * idx + 10);

                auto tk1_yy_z = primBuffer.data(tk1off + 18 * idx + 11);

                auto tk1_yz_x = primBuffer.data(tk1off + 18 * idx + 12);

                auto tk1_yz_y = primBuffer.data(tk1off + 18 * idx + 13);

                auto tk1_yz_z = primBuffer.data(tk1off + 18 * idx + 14);

                auto tk1_zz_x = primBuffer.data(tk1off + 18 * idx + 15);

                auto tk1_zz_y = primBuffer.data(tk1off + 18 * idx + 16);

                auto tk1_zz_z = primBuffer.data(tk1off + 18 * idx + 17);

                // set up pointers to (P|A(0)|D)^(m) integrals

                auto t20_x_xx = primBuffer.data(t20off + 18 * idx);

                auto t20_x_xy = primBuffer.data(t20off + 18 * idx + 1);

                auto t20_x_xz = primBuffer.data(t20off + 18 * idx + 2);

                auto t20_x_yy = primBuffer.data(t20off + 18 * idx + 3);

                auto t20_x_yz = primBuffer.data(t20off + 18 * idx + 4);

                auto t20_x_zz = primBuffer.data(t20off + 18 * idx + 5);

                auto t20_y_xx = primBuffer.data(t20off + 18 * idx + 6);

                auto t20_y_xy = primBuffer.data(t20off + 18 * idx + 7);

                auto t20_y_xz = primBuffer.data(t20off + 18 * idx + 8);

                auto t20_y_yy = primBuffer.data(t20off + 18 * idx + 9);

                auto t20_y_yz = primBuffer.data(t20off + 18 * idx + 10);

                auto t20_y_zz = primBuffer.data(t20off + 18 * idx + 11);

                auto t20_z_xx = primBuffer.data(t20off + 18 * idx + 12);

                auto t20_z_xy = primBuffer.data(t20off + 18 * idx + 13);

                auto t20_z_xz = primBuffer.data(t20off + 18 * idx + 14);

                auto t20_z_yy = primBuffer.data(t20off + 18 * idx + 15);

                auto t20_z_yz = primBuffer.data(t20off + 18 * idx + 16);

                auto t20_z_zz = primBuffer.data(t20off + 18 * idx + 17);

                // set up pointers to (P|A(0)|D)^(m+1) integrals

                auto t21_x_xx = primBuffer.data(t21off + 18 * idx);

                auto t21_x_xy = primBuffer.data(t21off + 18 * idx + 1);

                auto t21_x_xz = primBuffer.data(t21off + 18 * idx + 2);

                auto t21_x_yy = primBuffer.data(t21off + 18 * idx + 3);

                auto t21_x_yz = primBuffer.data(t21off + 18 * idx + 4);

                auto t21_x_zz = primBuffer.data(t21off + 18 * idx + 5);

                auto t21_y_xx = primBuffer.data(t21off + 18 * idx + 6);

                auto t21_y_xy = primBuffer.data(t21off + 18 * idx + 7);

                auto t21_y_xz = primBuffer.data(t21off + 18 * idx + 8);

                auto t21_y_yy = primBuffer.data(t21off + 18 * idx + 9);

                auto t21_y_yz = primBuffer.data(t21off + 18 * idx + 10);

                auto t21_y_zz = primBuffer.data(t21off + 18 * idx + 11);

                auto t21_z_xx = primBuffer.data(t21off + 18 * idx + 12);

                auto t21_z_xy = primBuffer.data(t21off + 18 * idx + 13);

                auto t21_z_xz = primBuffer.data(t21off + 18 * idx + 14);

                auto t21_z_yy = primBuffer.data(t21off + 18 * idx + 15);

                auto t21_z_yz = primBuffer.data(t21off + 18 * idx + 16);

                auto t21_z_zz = primBuffer.data(t21off + 18 * idx + 17);

                // set up pointers to (D|A(0)|D)^(m) integrals

                auto t10_xx_xx = primBuffer.data(t10off + 36 * idx);

                auto t10_xx_xy = primBuffer.data(t10off + 36 * idx + 1);

                auto t10_xx_xz = primBuffer.data(t10off + 36 * idx + 2);

                auto t10_xx_yy = primBuffer.data(t10off + 36 * idx + 3);

                auto t10_xx_yz = primBuffer.data(t10off + 36 * idx + 4);

                auto t10_xx_zz = primBuffer.data(t10off + 36 * idx + 5);

                auto t10_xy_xx = primBuffer.data(t10off + 36 * idx + 6);

                auto t10_xy_xy = primBuffer.data(t10off + 36 * idx + 7);

                auto t10_xy_xz = primBuffer.data(t10off + 36 * idx + 8);

                auto t10_xy_yy = primBuffer.data(t10off + 36 * idx + 9);

                auto t10_xy_yz = primBuffer.data(t10off + 36 * idx + 10);

                auto t10_xy_zz = primBuffer.data(t10off + 36 * idx + 11);

                auto t10_xz_xx = primBuffer.data(t10off + 36 * idx + 12);

                auto t10_xz_xy = primBuffer.data(t10off + 36 * idx + 13);

                auto t10_xz_xz = primBuffer.data(t10off + 36 * idx + 14);

                auto t10_xz_yy = primBuffer.data(t10off + 36 * idx + 15);

                auto t10_xz_yz = primBuffer.data(t10off + 36 * idx + 16);

                auto t10_xz_zz = primBuffer.data(t10off + 36 * idx + 17);

                auto t10_yy_xx = primBuffer.data(t10off + 36 * idx + 18);

                auto t10_yy_xy = primBuffer.data(t10off + 36 * idx + 19);

                auto t10_yy_xz = primBuffer.data(t10off + 36 * idx + 20);

                auto t10_yy_yy = primBuffer.data(t10off + 36 * idx + 21);

                auto t10_yy_yz = primBuffer.data(t10off + 36 * idx + 22);

                auto t10_yy_zz = primBuffer.data(t10off + 36 * idx + 23);

                auto t10_yz_xx = primBuffer.data(t10off + 36 * idx + 24);

                auto t10_yz_xy = primBuffer.data(t10off + 36 * idx + 25);

                auto t10_yz_xz = primBuffer.data(t10off + 36 * idx + 26);

                auto t10_yz_yy = primBuffer.data(t10off + 36 * idx + 27);

                auto t10_yz_yz = primBuffer.data(t10off + 36 * idx + 28);

                auto t10_yz_zz = primBuffer.data(t10off + 36 * idx + 29);

                auto t10_zz_xx = primBuffer.data(t10off + 36 * idx + 30);

                auto t10_zz_xy = primBuffer.data(t10off + 36 * idx + 31);

                auto t10_zz_xz = primBuffer.data(t10off + 36 * idx + 32);

                auto t10_zz_yy = primBuffer.data(t10off + 36 * idx + 33);

                auto t10_zz_yz = primBuffer.data(t10off + 36 * idx + 34);

                auto t10_zz_zz = primBuffer.data(t10off + 36 * idx + 35);

                // set up pointers to (D|A(0)|D)^(m+1) integrals

                auto t11_xx_xx = primBuffer.data(t11off + 36 * idx);

                auto t11_xx_xy = primBuffer.data(t11off + 36 * idx + 1);

                auto t11_xx_xz = primBuffer.data(t11off + 36 * idx + 2);

                auto t11_xx_yy = primBuffer.data(t11off + 36 * idx + 3);

                auto t11_xx_yz = primBuffer.data(t11off + 36 * idx + 4);

                auto t11_xx_zz = primBuffer.data(t11off + 36 * idx + 5);

                auto t11_xy_xx = primBuffer.data(t11off + 36 * idx + 6);

                auto t11_xy_xy = primBuffer.data(t11off + 36 * idx + 7);

                auto t11_xy_xz = primBuffer.data(t11off + 36 * idx + 8);

                auto t11_xy_yy = primBuffer.data(t11off + 36 * idx + 9);

                auto t11_xy_yz = primBuffer.data(t11off + 36 * idx + 10);

                auto t11_xy_zz = primBuffer.data(t11off + 36 * idx + 11);

                auto t11_xz_xx = primBuffer.data(t11off + 36 * idx + 12);

                auto t11_xz_xy = primBuffer.data(t11off + 36 * idx + 13);

                auto t11_xz_xz = primBuffer.data(t11off + 36 * idx + 14);

                auto t11_xz_yy = primBuffer.data(t11off + 36 * idx + 15);

                auto t11_xz_yz = primBuffer.data(t11off + 36 * idx + 16);

                auto t11_xz_zz = primBuffer.data(t11off + 36 * idx + 17);

                auto t11_yy_xx = primBuffer.data(t11off + 36 * idx + 18);

                auto t11_yy_xy = primBuffer.data(t11off + 36 * idx + 19);

                auto t11_yy_xz = primBuffer.data(t11off + 36 * idx + 20);

                auto t11_yy_yy = primBuffer.data(t11off + 36 * idx + 21);

                auto t11_yy_yz = primBuffer.data(t11off + 36 * idx + 22);

                auto t11_yy_zz = primBuffer.data(t11off + 36 * idx + 23);

                auto t11_yz_xx = primBuffer.data(t11off + 36 * idx + 24);

                auto t11_yz_xy = primBuffer.data(t11off + 36 * idx + 25);

                auto t11_yz_xz = primBuffer.data(t11off + 36 * idx + 26);

                auto t11_yz_yy = primBuffer.data(t11off + 36 * idx + 27);

                auto t11_yz_yz = primBuffer.data(t11off + 36 * idx + 28);

                auto t11_yz_zz = primBuffer.data(t11off + 36 * idx + 29);

                auto t11_zz_xx = primBuffer.data(t11off + 36 * idx + 30);

                auto t11_zz_xy = primBuffer.data(t11off + 36 * idx + 31);

                auto t11_zz_xz = primBuffer.data(t11off + 36 * idx + 32);

                auto t11_zz_yy = primBuffer.data(t11off + 36 * idx + 33);

                auto t11_zz_yz = primBuffer.data(t11off + 36 * idx + 34);

                auto t11_zz_zz = primBuffer.data(t11off + 36 * idx + 35);

                // set up pointers to (F|A(0)|D)^(m) integrals

                auto t_xxx_xx = primBuffer.data(toff + 60 * idx);

                auto t_xxx_xy = primBuffer.data(toff + 60 * idx + 1);

                auto t_xxx_xz = primBuffer.data(toff + 60 * idx + 2);

                auto t_xxx_yy = primBuffer.data(toff + 60 * idx + 3);

                auto t_xxx_yz = primBuffer.data(toff + 60 * idx + 4);

                auto t_xxx_zz = primBuffer.data(toff + 60 * idx + 5);

                auto t_xxy_xx = primBuffer.data(toff + 60 * idx + 6);

                auto t_xxy_xy = primBuffer.data(toff + 60 * idx + 7);

                auto t_xxy_xz = primBuffer.data(toff + 60 * idx + 8);

                auto t_xxy_yy = primBuffer.data(toff + 60 * idx + 9);

                auto t_xxy_yz = primBuffer.data(toff + 60 * idx + 10);

                auto t_xxy_zz = primBuffer.data(toff + 60 * idx + 11);

                auto t_xxz_xx = primBuffer.data(toff + 60 * idx + 12);

                auto t_xxz_xy = primBuffer.data(toff + 60 * idx + 13);

                auto t_xxz_xz = primBuffer.data(toff + 60 * idx + 14);

                auto t_xxz_yy = primBuffer.data(toff + 60 * idx + 15);

                auto t_xxz_yz = primBuffer.data(toff + 60 * idx + 16);

                auto t_xxz_zz = primBuffer.data(toff + 60 * idx + 17);

                auto t_xyy_xx = primBuffer.data(toff + 60 * idx + 18);

                auto t_xyy_xy = primBuffer.data(toff + 60 * idx + 19);

                auto t_xyy_xz = primBuffer.data(toff + 60 * idx + 20);

                auto t_xyy_yy = primBuffer.data(toff + 60 * idx + 21);

                auto t_xyy_yz = primBuffer.data(toff + 60 * idx + 22);

                auto t_xyy_zz = primBuffer.data(toff + 60 * idx + 23);

                auto t_xyz_xx = primBuffer.data(toff + 60 * idx + 24);

                auto t_xyz_xy = primBuffer.data(toff + 60 * idx + 25);

                auto t_xyz_xz = primBuffer.data(toff + 60 * idx + 26);

                auto t_xyz_yy = primBuffer.data(toff + 60 * idx + 27);

                auto t_xyz_yz = primBuffer.data(toff + 60 * idx + 28);

                auto t_xyz_zz = primBuffer.data(toff + 60 * idx + 29);

                auto t_xzz_xx = primBuffer.data(toff + 60 * idx + 30);

                auto t_xzz_xy = primBuffer.data(toff + 60 * idx + 31);

                auto t_xzz_xz = primBuffer.data(toff + 60 * idx + 32);

                auto t_xzz_yy = primBuffer.data(toff + 60 * idx + 33);

                auto t_xzz_yz = primBuffer.data(toff + 60 * idx + 34);

                auto t_xzz_zz = primBuffer.data(toff + 60 * idx + 35);

                auto t_yyy_xx = primBuffer.data(toff + 60 * idx + 36);

                auto t_yyy_xy = primBuffer.data(toff + 60 * idx + 37);

                auto t_yyy_xz = primBuffer.data(toff + 60 * idx + 38);

                auto t_yyy_yy = primBuffer.data(toff + 60 * idx + 39);

                auto t_yyy_yz = primBuffer.data(toff + 60 * idx + 40);

                auto t_yyy_zz = primBuffer.data(toff + 60 * idx + 41);

                auto t_yyz_xx = primBuffer.data(toff + 60 * idx + 42);

                auto t_yyz_xy = primBuffer.data(toff + 60 * idx + 43);

                auto t_yyz_xz = primBuffer.data(toff + 60 * idx + 44);

                auto t_yyz_yy = primBuffer.data(toff + 60 * idx + 45);

                auto t_yyz_yz = primBuffer.data(toff + 60 * idx + 46);

                auto t_yyz_zz = primBuffer.data(toff + 60 * idx + 47);

                auto t_yzz_xx = primBuffer.data(toff + 60 * idx + 48);

                auto t_yzz_xy = primBuffer.data(toff + 60 * idx + 49);

                auto t_yzz_xz = primBuffer.data(toff + 60 * idx + 50);

                auto t_yzz_yy = primBuffer.data(toff + 60 * idx + 51);

                auto t_yzz_yz = primBuffer.data(toff + 60 * idx + 52);

                auto t_yzz_zz = primBuffer.data(toff + 60 * idx + 53);

                auto t_zzz_xx = primBuffer.data(toff + 60 * idx + 54);

                auto t_zzz_xy = primBuffer.data(toff + 60 * idx + 55);

                auto t_zzz_xz = primBuffer.data(toff + 60 * idx + 56);

                auto t_zzz_yy = primBuffer.data(toff + 60 * idx + 57);

                auto t_zzz_yz = primBuffer.data(toff + 60 * idx + 58);

                auto t_zzz_zz = primBuffer.data(toff + 60 * idx + 59);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xx_x, tk0_xx_y,\
                                         tk0_xx_z, tk0_xy_x, tk0_xy_y, tk0_xy_z,\
                                         tk0_xz_x, tk0_xz_y, tk0_xz_z, tk0_yy_x,\
                                         tk0_yy_y, tk0_yy_z, tk0_yz_x, tk0_yz_y,\
                                         tk0_yz_z, tk0_zz_x, tk0_zz_y, tk0_zz_z,\
                                         tk1_xx_x, tk1_xx_y, tk1_xx_z, tk1_xy_x,\
                                         tk1_xy_y, tk1_xy_z, tk1_xz_x, tk1_xz_y,\
                                         tk1_xz_z, tk1_yy_x, tk1_yy_y, tk1_yy_z,\
                                         tk1_yz_x, tk1_yz_y, tk1_yz_z, tk1_zz_x,\
                                         tk1_zz_y, tk1_zz_z, t20_x_xx, t20_x_xy,\
                                         t20_x_xz, t20_x_yy, t20_x_yz, t20_x_zz,\
                                         t20_y_xx, t20_y_xy, t20_y_xz, t20_y_yy,\
                                         t20_y_yz, t20_y_zz, t20_z_xx, t20_z_xy,\
                                         t20_z_xz, t20_z_yy, t20_z_yz, t20_z_zz,\
                                         t21_x_xx, t21_x_xy, t21_x_xz, t21_x_yy,\
                                         t21_x_yz, t21_x_zz, t21_y_xx, t21_y_xy,\
                                         t21_y_xz, t21_y_yy, t21_y_yz, t21_y_zz,\
                                         t21_z_xx, t21_z_xy, t21_z_xz, t21_z_yy,\
                                         t21_z_yz, t21_z_zz, t10_xx_xx, t10_xx_xy,\
                                         t10_xx_xz, t10_xx_yy, t10_xx_yz, t10_xx_zz,\
                                         t10_xy_xx, t10_xy_xy, t10_xy_xz, t10_xy_yy,\
                                         t10_xy_yz, t10_xy_zz, t10_xz_xx, t10_xz_xy,\
                                         t10_xz_xz, t10_xz_yy, t10_xz_yz, t10_xz_zz,\
                                         t10_yy_xx, t10_yy_xy, t10_yy_xz, t10_yy_yy,\
                                         t10_yy_yz, t10_yy_zz, t10_yz_xx, t10_yz_xy,\
                                         t10_yz_xz, t10_yz_yy, t10_yz_yz, t10_yz_zz,\
                                         t10_zz_xx, t10_zz_xy, t10_zz_xz, t10_zz_yy,\
                                         t10_zz_yz, t10_zz_zz, t11_xx_xx, t11_xx_xy,\
                                         t11_xx_xz, t11_xx_yy, t11_xx_yz, t11_xx_zz,\
                                         t11_xy_xx, t11_xy_xy, t11_xy_xz, t11_xy_yy,\
                                         t11_xy_yz, t11_xy_zz, t11_xz_xx, t11_xz_xy,\
                                         t11_xz_xz, t11_xz_yy, t11_xz_yz, t11_xz_zz,\
                                         t11_yy_xx, t11_yy_xy, t11_yy_xz, t11_yy_yy,\
                                         t11_yy_yz, t11_yy_zz, t11_yz_xx, t11_yz_xy,\
                                         t11_yz_xz, t11_yz_yy, t11_yz_yz, t11_yz_zz,\
                                         t11_zz_xx, t11_zz_xy, t11_zz_xz, t11_zz_yy,\
                                         t11_zz_yz, t11_zz_zz, t_xxx_xx, t_xxx_xy,\
                                         t_xxx_xz, t_xxx_yy, t_xxx_yz, t_xxx_zz,\
                                         t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy,\
                                         t_xxy_yz, t_xxy_zz, t_xxz_xx, t_xxz_xy,\
                                         t_xxz_xz, t_xxz_yy, t_xxz_yz, t_xxz_zz,\
                                         t_xyy_xx, t_xyy_xy, t_xyy_xz, t_xyy_yy,\
                                         t_xyy_yz, t_xyy_zz, t_xyz_xx, t_xyz_xy,\
                                         t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz,\
                                         t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy,\
                                         t_xzz_yz, t_xzz_zz, t_yyy_xx, t_yyy_xy,\
                                         t_yyy_xz, t_yyy_yy, t_yyy_yz, t_yyy_zz,\
                                         t_yyz_xx, t_yyz_xy, t_yyz_xz, t_yyz_yy,\
                                         t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy,\
                                         t_yzz_xz, t_yzz_yy, t_yzz_yz, t_yzz_zz,\
                                         t_zzz_xx, t_zzz_xy, t_zzz_xz, t_zzz_yy,\
                                         t_zzz_yz, t_zzz_zz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxx_xx[k] = fra * t10_xx_xx[k] - frc * t11_xx_xx[k] + f2t * (2.0 * t20_x_xx[k] - 2.0 * t21_x_xx[k] + 2.0 * tk0_xx_x[k] - 2.0 * tk1_xx_x[k]);

                    t_xxx_xy[k] = fra * t10_xx_xy[k] - frc * t11_xx_xy[k] + f2t * (2.0 * t20_x_xy[k] - 2.0 * t21_x_xy[k] + tk0_xx_y[k] - tk1_xx_y[k]);

                    t_xxx_xz[k] = fra * t10_xx_xz[k] - frc * t11_xx_xz[k] + f2t * (2.0 * t20_x_xz[k] - 2.0 * t21_x_xz[k] + tk0_xx_z[k] - tk1_xx_z[k]);

                    t_xxx_yy[k] = fra * t10_xx_yy[k] - frc * t11_xx_yy[k] + f2t * (2.0 * t20_x_yy[k] - 2.0 * t21_x_yy[k]);

                    t_xxx_yz[k] = fra * t10_xx_yz[k] - frc * t11_xx_yz[k] + f2t * (2.0 * t20_x_yz[k] - 2.0 * t21_x_yz[k]);

                    t_xxx_zz[k] = fra * t10_xx_zz[k] - frc * t11_xx_zz[k] + f2t * (2.0 * t20_x_zz[k] - 2.0 * t21_x_zz[k]);

                    t_xxy_xx[k] = fra * t10_xy_xx[k] - frc * t11_xy_xx[k] + f2t * (t20_y_xx[k] - t21_y_xx[k] + 2.0 * tk0_xy_x[k] - 2.0 * tk1_xy_x[k]);

                    t_xxy_xy[k] = fra * t10_xy_xy[k] - frc * t11_xy_xy[k] + f2t * (t20_y_xy[k] - t21_y_xy[k] + tk0_xy_y[k] - tk1_xy_y[k]);

                    t_xxy_xz[k] = fra * t10_xy_xz[k] - frc * t11_xy_xz[k] + f2t * (t20_y_xz[k] - t21_y_xz[k] + tk0_xy_z[k] - tk1_xy_z[k]);

                    t_xxy_yy[k] = fra * t10_xy_yy[k] - frc * t11_xy_yy[k] + f2t * (t20_y_yy[k] - t21_y_yy[k]);

                    t_xxy_yz[k] = fra * t10_xy_yz[k] - frc * t11_xy_yz[k] + f2t * (t20_y_yz[k] - t21_y_yz[k]);

                    t_xxy_zz[k] = fra * t10_xy_zz[k] - frc * t11_xy_zz[k] + f2t * (t20_y_zz[k] - t21_y_zz[k]);

                    t_xxz_xx[k] = fra * t10_xz_xx[k] - frc * t11_xz_xx[k] + f2t * (t20_z_xx[k] - t21_z_xx[k] + 2.0 * tk0_xz_x[k] - 2.0 * tk1_xz_x[k]);

                    t_xxz_xy[k] = fra * t10_xz_xy[k] - frc * t11_xz_xy[k] + f2t * (t20_z_xy[k] - t21_z_xy[k] + tk0_xz_y[k] - tk1_xz_y[k]);

                    t_xxz_xz[k] = fra * t10_xz_xz[k] - frc * t11_xz_xz[k] + f2t * (t20_z_xz[k] - t21_z_xz[k] + tk0_xz_z[k] - tk1_xz_z[k]);

                    t_xxz_yy[k] = fra * t10_xz_yy[k] - frc * t11_xz_yy[k] + f2t * (t20_z_yy[k] - t21_z_yy[k]);

                    t_xxz_yz[k] = fra * t10_xz_yz[k] - frc * t11_xz_yz[k] + f2t * (t20_z_yz[k] - t21_z_yz[k]);

                    t_xxz_zz[k] = fra * t10_xz_zz[k] - frc * t11_xz_zz[k] + f2t * (t20_z_zz[k] - t21_z_zz[k]);

                    t_xyy_xx[k] = fra * t10_yy_xx[k] - frc * t11_yy_xx[k] + f2t * (2.0 * tk0_yy_x[k] - 2.0 * tk1_yy_x[k]);

                    t_xyy_xy[k] = fra * t10_yy_xy[k] - frc * t11_yy_xy[k] + f2t * (tk0_yy_y[k] - tk1_yy_y[k]);

                    t_xyy_xz[k] = fra * t10_yy_xz[k] - frc * t11_yy_xz[k] + f2t * (tk0_yy_z[k] - tk1_yy_z[k]);

                    t_xyy_yy[k] = fra * t10_yy_yy[k] - frc * t11_yy_yy[k];

                    t_xyy_yz[k] = fra * t10_yy_yz[k] - frc * t11_yy_yz[k];

                    t_xyy_zz[k] = fra * t10_yy_zz[k] - frc * t11_yy_zz[k];

                    t_xyz_xx[k] = fra * t10_yz_xx[k] - frc * t11_yz_xx[k] + f2t * (2.0 * tk0_yz_x[k] - 2.0 * tk1_yz_x[k]);

                    t_xyz_xy[k] = fra * t10_yz_xy[k] - frc * t11_yz_xy[k] + f2t * (tk0_yz_y[k] - tk1_yz_y[k]);

                    t_xyz_xz[k] = fra * t10_yz_xz[k] - frc * t11_yz_xz[k] + f2t * (tk0_yz_z[k] - tk1_yz_z[k]);

                    t_xyz_yy[k] = fra * t10_yz_yy[k] - frc * t11_yz_yy[k];

                    t_xyz_yz[k] = fra * t10_yz_yz[k] - frc * t11_yz_yz[k];

                    t_xyz_zz[k] = fra * t10_yz_zz[k] - frc * t11_yz_zz[k];

                    t_xzz_xx[k] = fra * t10_zz_xx[k] - frc * t11_zz_xx[k] + f2t * (2.0 * tk0_zz_x[k] - 2.0 * tk1_zz_x[k]);

                    t_xzz_xy[k] = fra * t10_zz_xy[k] - frc * t11_zz_xy[k] + f2t * (tk0_zz_y[k] - tk1_zz_y[k]);

                    t_xzz_xz[k] = fra * t10_zz_xz[k] - frc * t11_zz_xz[k] + f2t * (tk0_zz_z[k] - tk1_zz_z[k]);

                    t_xzz_yy[k] = fra * t10_zz_yy[k] - frc * t11_zz_yy[k];

                    t_xzz_yz[k] = fra * t10_zz_yz[k] - frc * t11_zz_yz[k];

                    t_xzz_zz[k] = fra * t10_zz_zz[k] - frc * t11_zz_zz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyy_xx[k] = fra * t10_yy_xx[k] - frc * t11_yy_xx[k] + f2t * (2.0 * t20_y_xx[k] - 2.0 * t21_y_xx[k]);

                    t_yyy_xy[k] = fra * t10_yy_xy[k] - frc * t11_yy_xy[k] + f2t * (2.0 * t20_y_xy[k] - 2.0 * t21_y_xy[k] + tk0_yy_x[k] - tk1_yy_x[k]);

                    t_yyy_xz[k] = fra * t10_yy_xz[k] - frc * t11_yy_xz[k] + f2t * (2.0 * t20_y_xz[k] - 2.0 * t21_y_xz[k]);

                    t_yyy_yy[k] = fra * t10_yy_yy[k] - frc * t11_yy_yy[k] + f2t * (2.0 * t20_y_yy[k] - 2.0 * t21_y_yy[k] + 2.0 * tk0_yy_y[k] - 2.0 * tk1_yy_y[k]);

                    t_yyy_yz[k] = fra * t10_yy_yz[k] - frc * t11_yy_yz[k] + f2t * (2.0 * t20_y_yz[k] - 2.0 * t21_y_yz[k] + tk0_yy_z[k] - tk1_yy_z[k]);

                    t_yyy_zz[k] = fra * t10_yy_zz[k] - frc * t11_yy_zz[k] + f2t * (2.0 * t20_y_zz[k] - 2.0 * t21_y_zz[k]);

                    t_yyz_xx[k] = fra * t10_yz_xx[k] - frc * t11_yz_xx[k] + f2t * (t20_z_xx[k] - t21_z_xx[k]);

                    t_yyz_xy[k] = fra * t10_yz_xy[k] - frc * t11_yz_xy[k] + f2t * (t20_z_xy[k] - t21_z_xy[k] + tk0_yz_x[k] - tk1_yz_x[k]);

                    t_yyz_xz[k] = fra * t10_yz_xz[k] - frc * t11_yz_xz[k] + f2t * (t20_z_xz[k] - t21_z_xz[k]);

                    t_yyz_yy[k] = fra * t10_yz_yy[k] - frc * t11_yz_yy[k] + f2t * (t20_z_yy[k] - t21_z_yy[k] + 2.0 * tk0_yz_y[k] - 2.0 * tk1_yz_y[k]);

                    t_yyz_yz[k] = fra * t10_yz_yz[k] - frc * t11_yz_yz[k] + f2t * (t20_z_yz[k] - t21_z_yz[k] + tk0_yz_z[k] - tk1_yz_z[k]);

                    t_yyz_zz[k] = fra * t10_yz_zz[k] - frc * t11_yz_zz[k] + f2t * (t20_z_zz[k] - t21_z_zz[k]);

                    t_yzz_xx[k] = fra * t10_zz_xx[k] - frc * t11_zz_xx[k];

                    t_yzz_xy[k] = fra * t10_zz_xy[k] - frc * t11_zz_xy[k] + f2t * (tk0_zz_x[k] - tk1_zz_x[k]);

                    t_yzz_xz[k] = fra * t10_zz_xz[k] - frc * t11_zz_xz[k];

                    t_yzz_yy[k] = fra * t10_zz_yy[k] - frc * t11_zz_yy[k] + f2t * (2.0 * tk0_zz_y[k] - 2.0 * tk1_zz_y[k]);

                    t_yzz_yz[k] = fra * t10_zz_yz[k] - frc * t11_zz_yz[k] + f2t * (tk0_zz_z[k] - tk1_zz_z[k]);

                    t_yzz_zz[k] = fra * t10_zz_zz[k] - frc * t11_zz_zz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzz_xx[k] = fra * t10_zz_xx[k] - frc * t11_zz_xx[k] + f2t * (2.0 * t20_z_xx[k] - 2.0 * t21_z_xx[k]);

                    t_zzz_xy[k] = fra * t10_zz_xy[k] - frc * t11_zz_xy[k] + f2t * (2.0 * t20_z_xy[k] - 2.0 * t21_z_xy[k]);

                    t_zzz_xz[k] = fra * t10_zz_xz[k] - frc * t11_zz_xz[k] + f2t * (2.0 * t20_z_xz[k] - 2.0 * t21_z_xz[k] + tk0_zz_x[k] - tk1_zz_x[k]);

                    t_zzz_yy[k] = fra * t10_zz_yy[k] - frc * t11_zz_yy[k] + f2t * (2.0 * t20_z_yy[k] - 2.0 * t21_z_yy[k]);

                    t_zzz_yz[k] = fra * t10_zz_yz[k] - frc * t11_zz_yz[k] + f2t * (2.0 * t20_z_yz[k] - 2.0 * t21_z_yz[k] + tk0_zz_y[k] - tk1_zz_y[k]);

                    t_zzz_zz[k] = fra * t10_zz_zz[k] - frc * t11_zz_zz[k] + f2t * (2.0 * t20_z_zz[k] - 2.0 * t21_z_zz[k] + 2.0 * tk0_zz_z[k] - 2.0 * tk1_zz_z[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForFF(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {3, 3, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 3, 3);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 3, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i + 1});

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

                // set up pointers to (D|A(0)|D)^(m) integrals

                auto tk0_xx_xx = primBuffer.data(tk0off + 36 * idx);

                auto tk0_xx_xy = primBuffer.data(tk0off + 36 * idx + 1);

                auto tk0_xx_xz = primBuffer.data(tk0off + 36 * idx + 2);

                auto tk0_xx_yy = primBuffer.data(tk0off + 36 * idx + 3);

                auto tk0_xx_yz = primBuffer.data(tk0off + 36 * idx + 4);

                auto tk0_xx_zz = primBuffer.data(tk0off + 36 * idx + 5);

                auto tk0_xy_xx = primBuffer.data(tk0off + 36 * idx + 6);

                auto tk0_xy_xy = primBuffer.data(tk0off + 36 * idx + 7);

                auto tk0_xy_xz = primBuffer.data(tk0off + 36 * idx + 8);

                auto tk0_xy_yy = primBuffer.data(tk0off + 36 * idx + 9);

                auto tk0_xy_yz = primBuffer.data(tk0off + 36 * idx + 10);

                auto tk0_xy_zz = primBuffer.data(tk0off + 36 * idx + 11);

                auto tk0_xz_xx = primBuffer.data(tk0off + 36 * idx + 12);

                auto tk0_xz_xy = primBuffer.data(tk0off + 36 * idx + 13);

                auto tk0_xz_xz = primBuffer.data(tk0off + 36 * idx + 14);

                auto tk0_xz_yy = primBuffer.data(tk0off + 36 * idx + 15);

                auto tk0_xz_yz = primBuffer.data(tk0off + 36 * idx + 16);

                auto tk0_xz_zz = primBuffer.data(tk0off + 36 * idx + 17);

                auto tk0_yy_xx = primBuffer.data(tk0off + 36 * idx + 18);

                auto tk0_yy_xy = primBuffer.data(tk0off + 36 * idx + 19);

                auto tk0_yy_xz = primBuffer.data(tk0off + 36 * idx + 20);

                auto tk0_yy_yy = primBuffer.data(tk0off + 36 * idx + 21);

                auto tk0_yy_yz = primBuffer.data(tk0off + 36 * idx + 22);

                auto tk0_yy_zz = primBuffer.data(tk0off + 36 * idx + 23);

                auto tk0_yz_xx = primBuffer.data(tk0off + 36 * idx + 24);

                auto tk0_yz_xy = primBuffer.data(tk0off + 36 * idx + 25);

                auto tk0_yz_xz = primBuffer.data(tk0off + 36 * idx + 26);

                auto tk0_yz_yy = primBuffer.data(tk0off + 36 * idx + 27);

                auto tk0_yz_yz = primBuffer.data(tk0off + 36 * idx + 28);

                auto tk0_yz_zz = primBuffer.data(tk0off + 36 * idx + 29);

                auto tk0_zz_xx = primBuffer.data(tk0off + 36 * idx + 30);

                auto tk0_zz_xy = primBuffer.data(tk0off + 36 * idx + 31);

                auto tk0_zz_xz = primBuffer.data(tk0off + 36 * idx + 32);

                auto tk0_zz_yy = primBuffer.data(tk0off + 36 * idx + 33);

                auto tk0_zz_yz = primBuffer.data(tk0off + 36 * idx + 34);

                auto tk0_zz_zz = primBuffer.data(tk0off + 36 * idx + 35);

                // set up pointers to (D|A(0)|D)^(m+1) integrals

                auto tk1_xx_xx = primBuffer.data(tk1off + 36 * idx);

                auto tk1_xx_xy = primBuffer.data(tk1off + 36 * idx + 1);

                auto tk1_xx_xz = primBuffer.data(tk1off + 36 * idx + 2);

                auto tk1_xx_yy = primBuffer.data(tk1off + 36 * idx + 3);

                auto tk1_xx_yz = primBuffer.data(tk1off + 36 * idx + 4);

                auto tk1_xx_zz = primBuffer.data(tk1off + 36 * idx + 5);

                auto tk1_xy_xx = primBuffer.data(tk1off + 36 * idx + 6);

                auto tk1_xy_xy = primBuffer.data(tk1off + 36 * idx + 7);

                auto tk1_xy_xz = primBuffer.data(tk1off + 36 * idx + 8);

                auto tk1_xy_yy = primBuffer.data(tk1off + 36 * idx + 9);

                auto tk1_xy_yz = primBuffer.data(tk1off + 36 * idx + 10);

                auto tk1_xy_zz = primBuffer.data(tk1off + 36 * idx + 11);

                auto tk1_xz_xx = primBuffer.data(tk1off + 36 * idx + 12);

                auto tk1_xz_xy = primBuffer.data(tk1off + 36 * idx + 13);

                auto tk1_xz_xz = primBuffer.data(tk1off + 36 * idx + 14);

                auto tk1_xz_yy = primBuffer.data(tk1off + 36 * idx + 15);

                auto tk1_xz_yz = primBuffer.data(tk1off + 36 * idx + 16);

                auto tk1_xz_zz = primBuffer.data(tk1off + 36 * idx + 17);

                auto tk1_yy_xx = primBuffer.data(tk1off + 36 * idx + 18);

                auto tk1_yy_xy = primBuffer.data(tk1off + 36 * idx + 19);

                auto tk1_yy_xz = primBuffer.data(tk1off + 36 * idx + 20);

                auto tk1_yy_yy = primBuffer.data(tk1off + 36 * idx + 21);

                auto tk1_yy_yz = primBuffer.data(tk1off + 36 * idx + 22);

                auto tk1_yy_zz = primBuffer.data(tk1off + 36 * idx + 23);

                auto tk1_yz_xx = primBuffer.data(tk1off + 36 * idx + 24);

                auto tk1_yz_xy = primBuffer.data(tk1off + 36 * idx + 25);

                auto tk1_yz_xz = primBuffer.data(tk1off + 36 * idx + 26);

                auto tk1_yz_yy = primBuffer.data(tk1off + 36 * idx + 27);

                auto tk1_yz_yz = primBuffer.data(tk1off + 36 * idx + 28);

                auto tk1_yz_zz = primBuffer.data(tk1off + 36 * idx + 29);

                auto tk1_zz_xx = primBuffer.data(tk1off + 36 * idx + 30);

                auto tk1_zz_xy = primBuffer.data(tk1off + 36 * idx + 31);

                auto tk1_zz_xz = primBuffer.data(tk1off + 36 * idx + 32);

                auto tk1_zz_yy = primBuffer.data(tk1off + 36 * idx + 33);

                auto tk1_zz_yz = primBuffer.data(tk1off + 36 * idx + 34);

                auto tk1_zz_zz = primBuffer.data(tk1off + 36 * idx + 35);

                // set up pointers to (P|A(0)|F)^(m) integrals

                auto t20_x_xxx = primBuffer.data(t20off + 30 * idx);

                auto t20_x_xxy = primBuffer.data(t20off + 30 * idx + 1);

                auto t20_x_xxz = primBuffer.data(t20off + 30 * idx + 2);

                auto t20_x_xyy = primBuffer.data(t20off + 30 * idx + 3);

                auto t20_x_xyz = primBuffer.data(t20off + 30 * idx + 4);

                auto t20_x_xzz = primBuffer.data(t20off + 30 * idx + 5);

                auto t20_x_yyy = primBuffer.data(t20off + 30 * idx + 6);

                auto t20_x_yyz = primBuffer.data(t20off + 30 * idx + 7);

                auto t20_x_yzz = primBuffer.data(t20off + 30 * idx + 8);

                auto t20_x_zzz = primBuffer.data(t20off + 30 * idx + 9);

                auto t20_y_xxx = primBuffer.data(t20off + 30 * idx + 10);

                auto t20_y_xxy = primBuffer.data(t20off + 30 * idx + 11);

                auto t20_y_xxz = primBuffer.data(t20off + 30 * idx + 12);

                auto t20_y_xyy = primBuffer.data(t20off + 30 * idx + 13);

                auto t20_y_xyz = primBuffer.data(t20off + 30 * idx + 14);

                auto t20_y_xzz = primBuffer.data(t20off + 30 * idx + 15);

                auto t20_y_yyy = primBuffer.data(t20off + 30 * idx + 16);

                auto t20_y_yyz = primBuffer.data(t20off + 30 * idx + 17);

                auto t20_y_yzz = primBuffer.data(t20off + 30 * idx + 18);

                auto t20_y_zzz = primBuffer.data(t20off + 30 * idx + 19);

                auto t20_z_xxx = primBuffer.data(t20off + 30 * idx + 20);

                auto t20_z_xxy = primBuffer.data(t20off + 30 * idx + 21);

                auto t20_z_xxz = primBuffer.data(t20off + 30 * idx + 22);

                auto t20_z_xyy = primBuffer.data(t20off + 30 * idx + 23);

                auto t20_z_xyz = primBuffer.data(t20off + 30 * idx + 24);

                auto t20_z_xzz = primBuffer.data(t20off + 30 * idx + 25);

                auto t20_z_yyy = primBuffer.data(t20off + 30 * idx + 26);

                auto t20_z_yyz = primBuffer.data(t20off + 30 * idx + 27);

                auto t20_z_yzz = primBuffer.data(t20off + 30 * idx + 28);

                auto t20_z_zzz = primBuffer.data(t20off + 30 * idx + 29);

                // set up pointers to (P|A(0)|F)^(m+1) integrals

                auto t21_x_xxx = primBuffer.data(t21off + 30 * idx);

                auto t21_x_xxy = primBuffer.data(t21off + 30 * idx + 1);

                auto t21_x_xxz = primBuffer.data(t21off + 30 * idx + 2);

                auto t21_x_xyy = primBuffer.data(t21off + 30 * idx + 3);

                auto t21_x_xyz = primBuffer.data(t21off + 30 * idx + 4);

                auto t21_x_xzz = primBuffer.data(t21off + 30 * idx + 5);

                auto t21_x_yyy = primBuffer.data(t21off + 30 * idx + 6);

                auto t21_x_yyz = primBuffer.data(t21off + 30 * idx + 7);

                auto t21_x_yzz = primBuffer.data(t21off + 30 * idx + 8);

                auto t21_x_zzz = primBuffer.data(t21off + 30 * idx + 9);

                auto t21_y_xxx = primBuffer.data(t21off + 30 * idx + 10);

                auto t21_y_xxy = primBuffer.data(t21off + 30 * idx + 11);

                auto t21_y_xxz = primBuffer.data(t21off + 30 * idx + 12);

                auto t21_y_xyy = primBuffer.data(t21off + 30 * idx + 13);

                auto t21_y_xyz = primBuffer.data(t21off + 30 * idx + 14);

                auto t21_y_xzz = primBuffer.data(t21off + 30 * idx + 15);

                auto t21_y_yyy = primBuffer.data(t21off + 30 * idx + 16);

                auto t21_y_yyz = primBuffer.data(t21off + 30 * idx + 17);

                auto t21_y_yzz = primBuffer.data(t21off + 30 * idx + 18);

                auto t21_y_zzz = primBuffer.data(t21off + 30 * idx + 19);

                auto t21_z_xxx = primBuffer.data(t21off + 30 * idx + 20);

                auto t21_z_xxy = primBuffer.data(t21off + 30 * idx + 21);

                auto t21_z_xxz = primBuffer.data(t21off + 30 * idx + 22);

                auto t21_z_xyy = primBuffer.data(t21off + 30 * idx + 23);

                auto t21_z_xyz = primBuffer.data(t21off + 30 * idx + 24);

                auto t21_z_xzz = primBuffer.data(t21off + 30 * idx + 25);

                auto t21_z_yyy = primBuffer.data(t21off + 30 * idx + 26);

                auto t21_z_yyz = primBuffer.data(t21off + 30 * idx + 27);

                auto t21_z_yzz = primBuffer.data(t21off + 30 * idx + 28);

                auto t21_z_zzz = primBuffer.data(t21off + 30 * idx + 29);

                // set up pointers to (D|A(0)|F)^(m) integrals

                auto t10_xx_xxx = primBuffer.data(t10off + 60 * idx);

                auto t10_xx_xxy = primBuffer.data(t10off + 60 * idx + 1);

                auto t10_xx_xxz = primBuffer.data(t10off + 60 * idx + 2);

                auto t10_xx_xyy = primBuffer.data(t10off + 60 * idx + 3);

                auto t10_xx_xyz = primBuffer.data(t10off + 60 * idx + 4);

                auto t10_xx_xzz = primBuffer.data(t10off + 60 * idx + 5);

                auto t10_xx_yyy = primBuffer.data(t10off + 60 * idx + 6);

                auto t10_xx_yyz = primBuffer.data(t10off + 60 * idx + 7);

                auto t10_xx_yzz = primBuffer.data(t10off + 60 * idx + 8);

                auto t10_xx_zzz = primBuffer.data(t10off + 60 * idx + 9);

                auto t10_xy_xxx = primBuffer.data(t10off + 60 * idx + 10);

                auto t10_xy_xxy = primBuffer.data(t10off + 60 * idx + 11);

                auto t10_xy_xxz = primBuffer.data(t10off + 60 * idx + 12);

                auto t10_xy_xyy = primBuffer.data(t10off + 60 * idx + 13);

                auto t10_xy_xyz = primBuffer.data(t10off + 60 * idx + 14);

                auto t10_xy_xzz = primBuffer.data(t10off + 60 * idx + 15);

                auto t10_xy_yyy = primBuffer.data(t10off + 60 * idx + 16);

                auto t10_xy_yyz = primBuffer.data(t10off + 60 * idx + 17);

                auto t10_xy_yzz = primBuffer.data(t10off + 60 * idx + 18);

                auto t10_xy_zzz = primBuffer.data(t10off + 60 * idx + 19);

                auto t10_xz_xxx = primBuffer.data(t10off + 60 * idx + 20);

                auto t10_xz_xxy = primBuffer.data(t10off + 60 * idx + 21);

                auto t10_xz_xxz = primBuffer.data(t10off + 60 * idx + 22);

                auto t10_xz_xyy = primBuffer.data(t10off + 60 * idx + 23);

                auto t10_xz_xyz = primBuffer.data(t10off + 60 * idx + 24);

                auto t10_xz_xzz = primBuffer.data(t10off + 60 * idx + 25);

                auto t10_xz_yyy = primBuffer.data(t10off + 60 * idx + 26);

                auto t10_xz_yyz = primBuffer.data(t10off + 60 * idx + 27);

                auto t10_xz_yzz = primBuffer.data(t10off + 60 * idx + 28);

                auto t10_xz_zzz = primBuffer.data(t10off + 60 * idx + 29);

                auto t10_yy_xxx = primBuffer.data(t10off + 60 * idx + 30);

                auto t10_yy_xxy = primBuffer.data(t10off + 60 * idx + 31);

                auto t10_yy_xxz = primBuffer.data(t10off + 60 * idx + 32);

                auto t10_yy_xyy = primBuffer.data(t10off + 60 * idx + 33);

                auto t10_yy_xyz = primBuffer.data(t10off + 60 * idx + 34);

                auto t10_yy_xzz = primBuffer.data(t10off + 60 * idx + 35);

                auto t10_yy_yyy = primBuffer.data(t10off + 60 * idx + 36);

                auto t10_yy_yyz = primBuffer.data(t10off + 60 * idx + 37);

                auto t10_yy_yzz = primBuffer.data(t10off + 60 * idx + 38);

                auto t10_yy_zzz = primBuffer.data(t10off + 60 * idx + 39);

                auto t10_yz_xxx = primBuffer.data(t10off + 60 * idx + 40);

                auto t10_yz_xxy = primBuffer.data(t10off + 60 * idx + 41);

                auto t10_yz_xxz = primBuffer.data(t10off + 60 * idx + 42);

                auto t10_yz_xyy = primBuffer.data(t10off + 60 * idx + 43);

                auto t10_yz_xyz = primBuffer.data(t10off + 60 * idx + 44);

                auto t10_yz_xzz = primBuffer.data(t10off + 60 * idx + 45);

                auto t10_yz_yyy = primBuffer.data(t10off + 60 * idx + 46);

                auto t10_yz_yyz = primBuffer.data(t10off + 60 * idx + 47);

                auto t10_yz_yzz = primBuffer.data(t10off + 60 * idx + 48);

                auto t10_yz_zzz = primBuffer.data(t10off + 60 * idx + 49);

                auto t10_zz_xxx = primBuffer.data(t10off + 60 * idx + 50);

                auto t10_zz_xxy = primBuffer.data(t10off + 60 * idx + 51);

                auto t10_zz_xxz = primBuffer.data(t10off + 60 * idx + 52);

                auto t10_zz_xyy = primBuffer.data(t10off + 60 * idx + 53);

                auto t10_zz_xyz = primBuffer.data(t10off + 60 * idx + 54);

                auto t10_zz_xzz = primBuffer.data(t10off + 60 * idx + 55);

                auto t10_zz_yyy = primBuffer.data(t10off + 60 * idx + 56);

                auto t10_zz_yyz = primBuffer.data(t10off + 60 * idx + 57);

                auto t10_zz_yzz = primBuffer.data(t10off + 60 * idx + 58);

                auto t10_zz_zzz = primBuffer.data(t10off + 60 * idx + 59);

                // set up pointers to (D|A(0)|F)^(m+1) integrals

                auto t11_xx_xxx = primBuffer.data(t11off + 60 * idx);

                auto t11_xx_xxy = primBuffer.data(t11off + 60 * idx + 1);

                auto t11_xx_xxz = primBuffer.data(t11off + 60 * idx + 2);

                auto t11_xx_xyy = primBuffer.data(t11off + 60 * idx + 3);

                auto t11_xx_xyz = primBuffer.data(t11off + 60 * idx + 4);

                auto t11_xx_xzz = primBuffer.data(t11off + 60 * idx + 5);

                auto t11_xx_yyy = primBuffer.data(t11off + 60 * idx + 6);

                auto t11_xx_yyz = primBuffer.data(t11off + 60 * idx + 7);

                auto t11_xx_yzz = primBuffer.data(t11off + 60 * idx + 8);

                auto t11_xx_zzz = primBuffer.data(t11off + 60 * idx + 9);

                auto t11_xy_xxx = primBuffer.data(t11off + 60 * idx + 10);

                auto t11_xy_xxy = primBuffer.data(t11off + 60 * idx + 11);

                auto t11_xy_xxz = primBuffer.data(t11off + 60 * idx + 12);

                auto t11_xy_xyy = primBuffer.data(t11off + 60 * idx + 13);

                auto t11_xy_xyz = primBuffer.data(t11off + 60 * idx + 14);

                auto t11_xy_xzz = primBuffer.data(t11off + 60 * idx + 15);

                auto t11_xy_yyy = primBuffer.data(t11off + 60 * idx + 16);

                auto t11_xy_yyz = primBuffer.data(t11off + 60 * idx + 17);

                auto t11_xy_yzz = primBuffer.data(t11off + 60 * idx + 18);

                auto t11_xy_zzz = primBuffer.data(t11off + 60 * idx + 19);

                auto t11_xz_xxx = primBuffer.data(t11off + 60 * idx + 20);

                auto t11_xz_xxy = primBuffer.data(t11off + 60 * idx + 21);

                auto t11_xz_xxz = primBuffer.data(t11off + 60 * idx + 22);

                auto t11_xz_xyy = primBuffer.data(t11off + 60 * idx + 23);

                auto t11_xz_xyz = primBuffer.data(t11off + 60 * idx + 24);

                auto t11_xz_xzz = primBuffer.data(t11off + 60 * idx + 25);

                auto t11_xz_yyy = primBuffer.data(t11off + 60 * idx + 26);

                auto t11_xz_yyz = primBuffer.data(t11off + 60 * idx + 27);

                auto t11_xz_yzz = primBuffer.data(t11off + 60 * idx + 28);

                auto t11_xz_zzz = primBuffer.data(t11off + 60 * idx + 29);

                auto t11_yy_xxx = primBuffer.data(t11off + 60 * idx + 30);

                auto t11_yy_xxy = primBuffer.data(t11off + 60 * idx + 31);

                auto t11_yy_xxz = primBuffer.data(t11off + 60 * idx + 32);

                auto t11_yy_xyy = primBuffer.data(t11off + 60 * idx + 33);

                auto t11_yy_xyz = primBuffer.data(t11off + 60 * idx + 34);

                auto t11_yy_xzz = primBuffer.data(t11off + 60 * idx + 35);

                auto t11_yy_yyy = primBuffer.data(t11off + 60 * idx + 36);

                auto t11_yy_yyz = primBuffer.data(t11off + 60 * idx + 37);

                auto t11_yy_yzz = primBuffer.data(t11off + 60 * idx + 38);

                auto t11_yy_zzz = primBuffer.data(t11off + 60 * idx + 39);

                auto t11_yz_xxx = primBuffer.data(t11off + 60 * idx + 40);

                auto t11_yz_xxy = primBuffer.data(t11off + 60 * idx + 41);

                auto t11_yz_xxz = primBuffer.data(t11off + 60 * idx + 42);

                auto t11_yz_xyy = primBuffer.data(t11off + 60 * idx + 43);

                auto t11_yz_xyz = primBuffer.data(t11off + 60 * idx + 44);

                auto t11_yz_xzz = primBuffer.data(t11off + 60 * idx + 45);

                auto t11_yz_yyy = primBuffer.data(t11off + 60 * idx + 46);

                auto t11_yz_yyz = primBuffer.data(t11off + 60 * idx + 47);

                auto t11_yz_yzz = primBuffer.data(t11off + 60 * idx + 48);

                auto t11_yz_zzz = primBuffer.data(t11off + 60 * idx + 49);

                auto t11_zz_xxx = primBuffer.data(t11off + 60 * idx + 50);

                auto t11_zz_xxy = primBuffer.data(t11off + 60 * idx + 51);

                auto t11_zz_xxz = primBuffer.data(t11off + 60 * idx + 52);

                auto t11_zz_xyy = primBuffer.data(t11off + 60 * idx + 53);

                auto t11_zz_xyz = primBuffer.data(t11off + 60 * idx + 54);

                auto t11_zz_xzz = primBuffer.data(t11off + 60 * idx + 55);

                auto t11_zz_yyy = primBuffer.data(t11off + 60 * idx + 56);

                auto t11_zz_yyz = primBuffer.data(t11off + 60 * idx + 57);

                auto t11_zz_yzz = primBuffer.data(t11off + 60 * idx + 58);

                auto t11_zz_zzz = primBuffer.data(t11off + 60 * idx + 59);

                // set up pointers to (F|A(0)|F)^(m) integrals

                auto t_xxx_xxx = primBuffer.data(toff + 100 * idx);

                auto t_xxx_xxy = primBuffer.data(toff + 100 * idx + 1);

                auto t_xxx_xxz = primBuffer.data(toff + 100 * idx + 2);

                auto t_xxx_xyy = primBuffer.data(toff + 100 * idx + 3);

                auto t_xxx_xyz = primBuffer.data(toff + 100 * idx + 4);

                auto t_xxx_xzz = primBuffer.data(toff + 100 * idx + 5);

                auto t_xxx_yyy = primBuffer.data(toff + 100 * idx + 6);

                auto t_xxx_yyz = primBuffer.data(toff + 100 * idx + 7);

                auto t_xxx_yzz = primBuffer.data(toff + 100 * idx + 8);

                auto t_xxx_zzz = primBuffer.data(toff + 100 * idx + 9);

                auto t_xxy_xxx = primBuffer.data(toff + 100 * idx + 10);

                auto t_xxy_xxy = primBuffer.data(toff + 100 * idx + 11);

                auto t_xxy_xxz = primBuffer.data(toff + 100 * idx + 12);

                auto t_xxy_xyy = primBuffer.data(toff + 100 * idx + 13);

                auto t_xxy_xyz = primBuffer.data(toff + 100 * idx + 14);

                auto t_xxy_xzz = primBuffer.data(toff + 100 * idx + 15);

                auto t_xxy_yyy = primBuffer.data(toff + 100 * idx + 16);

                auto t_xxy_yyz = primBuffer.data(toff + 100 * idx + 17);

                auto t_xxy_yzz = primBuffer.data(toff + 100 * idx + 18);

                auto t_xxy_zzz = primBuffer.data(toff + 100 * idx + 19);

                auto t_xxz_xxx = primBuffer.data(toff + 100 * idx + 20);

                auto t_xxz_xxy = primBuffer.data(toff + 100 * idx + 21);

                auto t_xxz_xxz = primBuffer.data(toff + 100 * idx + 22);

                auto t_xxz_xyy = primBuffer.data(toff + 100 * idx + 23);

                auto t_xxz_xyz = primBuffer.data(toff + 100 * idx + 24);

                auto t_xxz_xzz = primBuffer.data(toff + 100 * idx + 25);

                auto t_xxz_yyy = primBuffer.data(toff + 100 * idx + 26);

                auto t_xxz_yyz = primBuffer.data(toff + 100 * idx + 27);

                auto t_xxz_yzz = primBuffer.data(toff + 100 * idx + 28);

                auto t_xxz_zzz = primBuffer.data(toff + 100 * idx + 29);

                auto t_xyy_xxx = primBuffer.data(toff + 100 * idx + 30);

                auto t_xyy_xxy = primBuffer.data(toff + 100 * idx + 31);

                auto t_xyy_xxz = primBuffer.data(toff + 100 * idx + 32);

                auto t_xyy_xyy = primBuffer.data(toff + 100 * idx + 33);

                auto t_xyy_xyz = primBuffer.data(toff + 100 * idx + 34);

                auto t_xyy_xzz = primBuffer.data(toff + 100 * idx + 35);

                auto t_xyy_yyy = primBuffer.data(toff + 100 * idx + 36);

                auto t_xyy_yyz = primBuffer.data(toff + 100 * idx + 37);

                auto t_xyy_yzz = primBuffer.data(toff + 100 * idx + 38);

                auto t_xyy_zzz = primBuffer.data(toff + 100 * idx + 39);

                auto t_xyz_xxx = primBuffer.data(toff + 100 * idx + 40);

                auto t_xyz_xxy = primBuffer.data(toff + 100 * idx + 41);

                auto t_xyz_xxz = primBuffer.data(toff + 100 * idx + 42);

                auto t_xyz_xyy = primBuffer.data(toff + 100 * idx + 43);

                auto t_xyz_xyz = primBuffer.data(toff + 100 * idx + 44);

                auto t_xyz_xzz = primBuffer.data(toff + 100 * idx + 45);

                auto t_xyz_yyy = primBuffer.data(toff + 100 * idx + 46);

                auto t_xyz_yyz = primBuffer.data(toff + 100 * idx + 47);

                auto t_xyz_yzz = primBuffer.data(toff + 100 * idx + 48);

                auto t_xyz_zzz = primBuffer.data(toff + 100 * idx + 49);

                auto t_xzz_xxx = primBuffer.data(toff + 100 * idx + 50);

                auto t_xzz_xxy = primBuffer.data(toff + 100 * idx + 51);

                auto t_xzz_xxz = primBuffer.data(toff + 100 * idx + 52);

                auto t_xzz_xyy = primBuffer.data(toff + 100 * idx + 53);

                auto t_xzz_xyz = primBuffer.data(toff + 100 * idx + 54);

                auto t_xzz_xzz = primBuffer.data(toff + 100 * idx + 55);

                auto t_xzz_yyy = primBuffer.data(toff + 100 * idx + 56);

                auto t_xzz_yyz = primBuffer.data(toff + 100 * idx + 57);

                auto t_xzz_yzz = primBuffer.data(toff + 100 * idx + 58);

                auto t_xzz_zzz = primBuffer.data(toff + 100 * idx + 59);

                auto t_yyy_xxx = primBuffer.data(toff + 100 * idx + 60);

                auto t_yyy_xxy = primBuffer.data(toff + 100 * idx + 61);

                auto t_yyy_xxz = primBuffer.data(toff + 100 * idx + 62);

                auto t_yyy_xyy = primBuffer.data(toff + 100 * idx + 63);

                auto t_yyy_xyz = primBuffer.data(toff + 100 * idx + 64);

                auto t_yyy_xzz = primBuffer.data(toff + 100 * idx + 65);

                auto t_yyy_yyy = primBuffer.data(toff + 100 * idx + 66);

                auto t_yyy_yyz = primBuffer.data(toff + 100 * idx + 67);

                auto t_yyy_yzz = primBuffer.data(toff + 100 * idx + 68);

                auto t_yyy_zzz = primBuffer.data(toff + 100 * idx + 69);

                auto t_yyz_xxx = primBuffer.data(toff + 100 * idx + 70);

                auto t_yyz_xxy = primBuffer.data(toff + 100 * idx + 71);

                auto t_yyz_xxz = primBuffer.data(toff + 100 * idx + 72);

                auto t_yyz_xyy = primBuffer.data(toff + 100 * idx + 73);

                auto t_yyz_xyz = primBuffer.data(toff + 100 * idx + 74);

                auto t_yyz_xzz = primBuffer.data(toff + 100 * idx + 75);

                auto t_yyz_yyy = primBuffer.data(toff + 100 * idx + 76);

                auto t_yyz_yyz = primBuffer.data(toff + 100 * idx + 77);

                auto t_yyz_yzz = primBuffer.data(toff + 100 * idx + 78);

                auto t_yyz_zzz = primBuffer.data(toff + 100 * idx + 79);

                auto t_yzz_xxx = primBuffer.data(toff + 100 * idx + 80);

                auto t_yzz_xxy = primBuffer.data(toff + 100 * idx + 81);

                auto t_yzz_xxz = primBuffer.data(toff + 100 * idx + 82);

                auto t_yzz_xyy = primBuffer.data(toff + 100 * idx + 83);

                auto t_yzz_xyz = primBuffer.data(toff + 100 * idx + 84);

                auto t_yzz_xzz = primBuffer.data(toff + 100 * idx + 85);

                auto t_yzz_yyy = primBuffer.data(toff + 100 * idx + 86);

                auto t_yzz_yyz = primBuffer.data(toff + 100 * idx + 87);

                auto t_yzz_yzz = primBuffer.data(toff + 100 * idx + 88);

                auto t_yzz_zzz = primBuffer.data(toff + 100 * idx + 89);

                auto t_zzz_xxx = primBuffer.data(toff + 100 * idx + 90);

                auto t_zzz_xxy = primBuffer.data(toff + 100 * idx + 91);

                auto t_zzz_xxz = primBuffer.data(toff + 100 * idx + 92);

                auto t_zzz_xyy = primBuffer.data(toff + 100 * idx + 93);

                auto t_zzz_xyz = primBuffer.data(toff + 100 * idx + 94);

                auto t_zzz_xzz = primBuffer.data(toff + 100 * idx + 95);

                auto t_zzz_yyy = primBuffer.data(toff + 100 * idx + 96);

                auto t_zzz_yyz = primBuffer.data(toff + 100 * idx + 97);

                auto t_zzz_yzz = primBuffer.data(toff + 100 * idx + 98);

                auto t_zzz_zzz = primBuffer.data(toff + 100 * idx + 99);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xx_xx, tk0_xx_xy,\
                                         tk0_xx_xz, tk0_xx_yy, tk0_xx_yz, tk0_xx_zz,\
                                         tk0_xy_xx, tk0_xy_xy, tk0_xy_xz, tk0_xy_yy,\
                                         tk0_xy_yz, tk0_xy_zz, tk0_xz_xx, tk0_xz_xy,\
                                         tk0_xz_xz, tk0_xz_yy, tk0_xz_yz, tk0_xz_zz,\
                                         tk0_yy_xx, tk0_yy_xy, tk0_yy_xz, tk0_yy_yy,\
                                         tk0_yy_yz, tk0_yy_zz, tk0_yz_xx, tk0_yz_xy,\
                                         tk0_yz_xz, tk0_yz_yy, tk0_yz_yz, tk0_yz_zz,\
                                         tk0_zz_xx, tk0_zz_xy, tk0_zz_xz, tk0_zz_yy,\
                                         tk0_zz_yz, tk0_zz_zz, tk1_xx_xx, tk1_xx_xy,\
                                         tk1_xx_xz, tk1_xx_yy, tk1_xx_yz, tk1_xx_zz,\
                                         tk1_xy_xx, tk1_xy_xy, tk1_xy_xz, tk1_xy_yy,\
                                         tk1_xy_yz, tk1_xy_zz, tk1_xz_xx, tk1_xz_xy,\
                                         tk1_xz_xz, tk1_xz_yy, tk1_xz_yz, tk1_xz_zz,\
                                         tk1_yy_xx, tk1_yy_xy, tk1_yy_xz, tk1_yy_yy,\
                                         tk1_yy_yz, tk1_yy_zz, tk1_yz_xx, tk1_yz_xy,\
                                         tk1_yz_xz, tk1_yz_yy, tk1_yz_yz, tk1_yz_zz,\
                                         tk1_zz_xx, tk1_zz_xy, tk1_zz_xz, tk1_zz_yy,\
                                         tk1_zz_yz, tk1_zz_zz, t20_x_xxx, t20_x_xxy,\
                                         t20_x_xxz, t20_x_xyy, t20_x_xyz, t20_x_xzz,\
                                         t20_x_yyy, t20_x_yyz, t20_x_yzz, t20_x_zzz,\
                                         t20_y_xxx, t20_y_xxy, t20_y_xxz, t20_y_xyy,\
                                         t20_y_xyz, t20_y_xzz, t20_y_yyy, t20_y_yyz,\
                                         t20_y_yzz, t20_y_zzz, t20_z_xxx, t20_z_xxy,\
                                         t20_z_xxz, t20_z_xyy, t20_z_xyz, t20_z_xzz,\
                                         t20_z_yyy, t20_z_yyz, t20_z_yzz, t20_z_zzz,\
                                         t21_x_xxx, t21_x_xxy, t21_x_xxz, t21_x_xyy,\
                                         t21_x_xyz, t21_x_xzz, t21_x_yyy, t21_x_yyz,\
                                         t21_x_yzz, t21_x_zzz, t21_y_xxx, t21_y_xxy,\
                                         t21_y_xxz, t21_y_xyy, t21_y_xyz, t21_y_xzz,\
                                         t21_y_yyy, t21_y_yyz, t21_y_yzz, t21_y_zzz,\
                                         t21_z_xxx, t21_z_xxy, t21_z_xxz, t21_z_xyy,\
                                         t21_z_xyz, t21_z_xzz, t21_z_yyy, t21_z_yyz,\
                                         t21_z_yzz, t21_z_zzz, t10_xx_xxx, t10_xx_xxy,\
                                         t10_xx_xxz, t10_xx_xyy, t10_xx_xyz, t10_xx_xzz,\
                                         t10_xx_yyy, t10_xx_yyz, t10_xx_yzz, t10_xx_zzz,\
                                         t10_xy_xxx, t10_xy_xxy, t10_xy_xxz, t10_xy_xyy,\
                                         t10_xy_xyz, t10_xy_xzz, t10_xy_yyy, t10_xy_yyz,\
                                         t10_xy_yzz, t10_xy_zzz, t10_xz_xxx, t10_xz_xxy,\
                                         t10_xz_xxz, t10_xz_xyy, t10_xz_xyz, t10_xz_xzz,\
                                         t10_xz_yyy, t10_xz_yyz, t10_xz_yzz, t10_xz_zzz,\
                                         t10_yy_xxx, t10_yy_xxy, t10_yy_xxz, t10_yy_xyy,\
                                         t10_yy_xyz, t10_yy_xzz, t10_yy_yyy, t10_yy_yyz,\
                                         t10_yy_yzz, t10_yy_zzz, t10_yz_xxx, t10_yz_xxy,\
                                         t10_yz_xxz, t10_yz_xyy, t10_yz_xyz, t10_yz_xzz,\
                                         t10_yz_yyy, t10_yz_yyz, t10_yz_yzz, t10_yz_zzz,\
                                         t10_zz_xxx, t10_zz_xxy, t10_zz_xxz, t10_zz_xyy,\
                                         t10_zz_xyz, t10_zz_xzz, t10_zz_yyy, t10_zz_yyz,\
                                         t10_zz_yzz, t10_zz_zzz, t11_xx_xxx, t11_xx_xxy,\
                                         t11_xx_xxz, t11_xx_xyy, t11_xx_xyz, t11_xx_xzz,\
                                         t11_xx_yyy, t11_xx_yyz, t11_xx_yzz, t11_xx_zzz,\
                                         t11_xy_xxx, t11_xy_xxy, t11_xy_xxz, t11_xy_xyy,\
                                         t11_xy_xyz, t11_xy_xzz, t11_xy_yyy, t11_xy_yyz,\
                                         t11_xy_yzz, t11_xy_zzz, t11_xz_xxx, t11_xz_xxy,\
                                         t11_xz_xxz, t11_xz_xyy, t11_xz_xyz, t11_xz_xzz,\
                                         t11_xz_yyy, t11_xz_yyz, t11_xz_yzz, t11_xz_zzz,\
                                         t11_yy_xxx, t11_yy_xxy, t11_yy_xxz, t11_yy_xyy,\
                                         t11_yy_xyz, t11_yy_xzz, t11_yy_yyy, t11_yy_yyz,\
                                         t11_yy_yzz, t11_yy_zzz, t11_yz_xxx, t11_yz_xxy,\
                                         t11_yz_xxz, t11_yz_xyy, t11_yz_xyz, t11_yz_xzz,\
                                         t11_yz_yyy, t11_yz_yyz, t11_yz_yzz, t11_yz_zzz,\
                                         t11_zz_xxx, t11_zz_xxy, t11_zz_xxz, t11_zz_xyy,\
                                         t11_zz_xyz, t11_zz_xzz, t11_zz_yyy, t11_zz_yyz,\
                                         t11_zz_yzz, t11_zz_zzz, t_xxx_xxx, t_xxx_xxy,\
                                         t_xxx_xxz, t_xxx_xyy, t_xxx_xyz, t_xxx_xzz,\
                                         t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz,\
                                         t_xxy_xxx, t_xxy_xxy, t_xxy_xxz, t_xxy_xyy,\
                                         t_xxy_xyz, t_xxy_xzz, t_xxy_yyy, t_xxy_yyz,\
                                         t_xxy_yzz, t_xxy_zzz, t_xxz_xxx, t_xxz_xxy,\
                                         t_xxz_xxz, t_xxz_xyy, t_xxz_xyz, t_xxz_xzz,\
                                         t_xxz_yyy, t_xxz_yyz, t_xxz_yzz, t_xxz_zzz,\
                                         t_xyy_xxx, t_xyy_xxy, t_xyy_xxz, t_xyy_xyy,\
                                         t_xyy_xyz, t_xyy_xzz, t_xyy_yyy, t_xyy_yyz,\
                                         t_xyy_yzz, t_xyy_zzz, t_xyz_xxx, t_xyz_xxy,\
                                         t_xyz_xxz, t_xyz_xyy, t_xyz_xyz, t_xyz_xzz,\
                                         t_xyz_yyy, t_xyz_yyz, t_xyz_yzz, t_xyz_zzz,\
                                         t_xzz_xxx, t_xzz_xxy, t_xzz_xxz, t_xzz_xyy,\
                                         t_xzz_xyz, t_xzz_xzz, t_xzz_yyy, t_xzz_yyz,\
                                         t_xzz_yzz, t_xzz_zzz, t_yyy_xxx, t_yyy_xxy,\
                                         t_yyy_xxz, t_yyy_xyy, t_yyy_xyz, t_yyy_xzz,\
                                         t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, t_yyy_zzz,\
                                         t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy,\
                                         t_yyz_xyz, t_yyz_xzz, t_yyz_yyy, t_yyz_yyz,\
                                         t_yyz_yzz, t_yyz_zzz, t_yzz_xxx, t_yzz_xxy,\
                                         t_yzz_xxz, t_yzz_xyy, t_yzz_xyz, t_yzz_xzz,\
                                         t_yzz_yyy, t_yzz_yyz, t_yzz_yzz, t_yzz_zzz,\
                                         t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, t_zzz_xyy,\
                                         t_zzz_xyz, t_zzz_xzz, t_zzz_yyy, t_zzz_yyz,\
                                         t_zzz_yzz, t_zzz_zzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxx_xxx[k] = fra * t10_xx_xxx[k] - frc * t11_xx_xxx[k] + f2t * (2.0 * t20_x_xxx[k] - 2.0 * t21_x_xxx[k] + 3.0 * tk0_xx_xx[k] - 3.0 * tk1_xx_xx[k]);

                    t_xxx_xxy[k] = fra * t10_xx_xxy[k] - frc * t11_xx_xxy[k] + f2t * (2.0 * t20_x_xxy[k] - 2.0 * t21_x_xxy[k] + 2.0 * tk0_xx_xy[k] - 2.0 * tk1_xx_xy[k]);

                    t_xxx_xxz[k] = fra * t10_xx_xxz[k] - frc * t11_xx_xxz[k] + f2t * (2.0 * t20_x_xxz[k] - 2.0 * t21_x_xxz[k] + 2.0 * tk0_xx_xz[k] - 2.0 * tk1_xx_xz[k]);

                    t_xxx_xyy[k] = fra * t10_xx_xyy[k] - frc * t11_xx_xyy[k] + f2t * (2.0 * t20_x_xyy[k] - 2.0 * t21_x_xyy[k] + tk0_xx_yy[k] - tk1_xx_yy[k]);

                    t_xxx_xyz[k] = fra * t10_xx_xyz[k] - frc * t11_xx_xyz[k] + f2t * (2.0 * t20_x_xyz[k] - 2.0 * t21_x_xyz[k] + tk0_xx_yz[k] - tk1_xx_yz[k]);

                    t_xxx_xzz[k] = fra * t10_xx_xzz[k] - frc * t11_xx_xzz[k] + f2t * (2.0 * t20_x_xzz[k] - 2.0 * t21_x_xzz[k] + tk0_xx_zz[k] - tk1_xx_zz[k]);

                    t_xxx_yyy[k] = fra * t10_xx_yyy[k] - frc * t11_xx_yyy[k] + f2t * (2.0 * t20_x_yyy[k] - 2.0 * t21_x_yyy[k]);

                    t_xxx_yyz[k] = fra * t10_xx_yyz[k] - frc * t11_xx_yyz[k] + f2t * (2.0 * t20_x_yyz[k] - 2.0 * t21_x_yyz[k]);

                    t_xxx_yzz[k] = fra * t10_xx_yzz[k] - frc * t11_xx_yzz[k] + f2t * (2.0 * t20_x_yzz[k] - 2.0 * t21_x_yzz[k]);

                    t_xxx_zzz[k] = fra * t10_xx_zzz[k] - frc * t11_xx_zzz[k] + f2t * (2.0 * t20_x_zzz[k] - 2.0 * t21_x_zzz[k]);

                    t_xxy_xxx[k] = fra * t10_xy_xxx[k] - frc * t11_xy_xxx[k] + f2t * (t20_y_xxx[k] - t21_y_xxx[k] + 3.0 * tk0_xy_xx[k] - 3.0 * tk1_xy_xx[k]);

                    t_xxy_xxy[k] = fra * t10_xy_xxy[k] - frc * t11_xy_xxy[k] + f2t * (t20_y_xxy[k] - t21_y_xxy[k] + 2.0 * tk0_xy_xy[k] - 2.0 * tk1_xy_xy[k]);

                    t_xxy_xxz[k] = fra * t10_xy_xxz[k] - frc * t11_xy_xxz[k] + f2t * (t20_y_xxz[k] - t21_y_xxz[k] + 2.0 * tk0_xy_xz[k] - 2.0 * tk1_xy_xz[k]);

                    t_xxy_xyy[k] = fra * t10_xy_xyy[k] - frc * t11_xy_xyy[k] + f2t * (t20_y_xyy[k] - t21_y_xyy[k] + tk0_xy_yy[k] - tk1_xy_yy[k]);

                    t_xxy_xyz[k] = fra * t10_xy_xyz[k] - frc * t11_xy_xyz[k] + f2t * (t20_y_xyz[k] - t21_y_xyz[k] + tk0_xy_yz[k] - tk1_xy_yz[k]);

                    t_xxy_xzz[k] = fra * t10_xy_xzz[k] - frc * t11_xy_xzz[k] + f2t * (t20_y_xzz[k] - t21_y_xzz[k] + tk0_xy_zz[k] - tk1_xy_zz[k]);

                    t_xxy_yyy[k] = fra * t10_xy_yyy[k] - frc * t11_xy_yyy[k] + f2t * (t20_y_yyy[k] - t21_y_yyy[k]);

                    t_xxy_yyz[k] = fra * t10_xy_yyz[k] - frc * t11_xy_yyz[k] + f2t * (t20_y_yyz[k] - t21_y_yyz[k]);

                    t_xxy_yzz[k] = fra * t10_xy_yzz[k] - frc * t11_xy_yzz[k] + f2t * (t20_y_yzz[k] - t21_y_yzz[k]);

                    t_xxy_zzz[k] = fra * t10_xy_zzz[k] - frc * t11_xy_zzz[k] + f2t * (t20_y_zzz[k] - t21_y_zzz[k]);

                    t_xxz_xxx[k] = fra * t10_xz_xxx[k] - frc * t11_xz_xxx[k] + f2t * (t20_z_xxx[k] - t21_z_xxx[k] + 3.0 * tk0_xz_xx[k] - 3.0 * tk1_xz_xx[k]);

                    t_xxz_xxy[k] = fra * t10_xz_xxy[k] - frc * t11_xz_xxy[k] + f2t * (t20_z_xxy[k] - t21_z_xxy[k] + 2.0 * tk0_xz_xy[k] - 2.0 * tk1_xz_xy[k]);

                    t_xxz_xxz[k] = fra * t10_xz_xxz[k] - frc * t11_xz_xxz[k] + f2t * (t20_z_xxz[k] - t21_z_xxz[k] + 2.0 * tk0_xz_xz[k] - 2.0 * tk1_xz_xz[k]);

                    t_xxz_xyy[k] = fra * t10_xz_xyy[k] - frc * t11_xz_xyy[k] + f2t * (t20_z_xyy[k] - t21_z_xyy[k] + tk0_xz_yy[k] - tk1_xz_yy[k]);

                    t_xxz_xyz[k] = fra * t10_xz_xyz[k] - frc * t11_xz_xyz[k] + f2t * (t20_z_xyz[k] - t21_z_xyz[k] + tk0_xz_yz[k] - tk1_xz_yz[k]);

                    t_xxz_xzz[k] = fra * t10_xz_xzz[k] - frc * t11_xz_xzz[k] + f2t * (t20_z_xzz[k] - t21_z_xzz[k] + tk0_xz_zz[k] - tk1_xz_zz[k]);

                    t_xxz_yyy[k] = fra * t10_xz_yyy[k] - frc * t11_xz_yyy[k] + f2t * (t20_z_yyy[k] - t21_z_yyy[k]);

                    t_xxz_yyz[k] = fra * t10_xz_yyz[k] - frc * t11_xz_yyz[k] + f2t * (t20_z_yyz[k] - t21_z_yyz[k]);

                    t_xxz_yzz[k] = fra * t10_xz_yzz[k] - frc * t11_xz_yzz[k] + f2t * (t20_z_yzz[k] - t21_z_yzz[k]);

                    t_xxz_zzz[k] = fra * t10_xz_zzz[k] - frc * t11_xz_zzz[k] + f2t * (t20_z_zzz[k] - t21_z_zzz[k]);

                    t_xyy_xxx[k] = fra * t10_yy_xxx[k] - frc * t11_yy_xxx[k] + f2t * (3.0 * tk0_yy_xx[k] - 3.0 * tk1_yy_xx[k]);

                    t_xyy_xxy[k] = fra * t10_yy_xxy[k] - frc * t11_yy_xxy[k] + f2t * (2.0 * tk0_yy_xy[k] - 2.0 * tk1_yy_xy[k]);

                    t_xyy_xxz[k] = fra * t10_yy_xxz[k] - frc * t11_yy_xxz[k] + f2t * (2.0 * tk0_yy_xz[k] - 2.0 * tk1_yy_xz[k]);

                    t_xyy_xyy[k] = fra * t10_yy_xyy[k] - frc * t11_yy_xyy[k] + f2t * (tk0_yy_yy[k] - tk1_yy_yy[k]);

                    t_xyy_xyz[k] = fra * t10_yy_xyz[k] - frc * t11_yy_xyz[k] + f2t * (tk0_yy_yz[k] - tk1_yy_yz[k]);

                    t_xyy_xzz[k] = fra * t10_yy_xzz[k] - frc * t11_yy_xzz[k] + f2t * (tk0_yy_zz[k] - tk1_yy_zz[k]);

                    t_xyy_yyy[k] = fra * t10_yy_yyy[k] - frc * t11_yy_yyy[k];

                    t_xyy_yyz[k] = fra * t10_yy_yyz[k] - frc * t11_yy_yyz[k];

                    t_xyy_yzz[k] = fra * t10_yy_yzz[k] - frc * t11_yy_yzz[k];

                    t_xyy_zzz[k] = fra * t10_yy_zzz[k] - frc * t11_yy_zzz[k];

                    t_xyz_xxx[k] = fra * t10_yz_xxx[k] - frc * t11_yz_xxx[k] + f2t * (3.0 * tk0_yz_xx[k] - 3.0 * tk1_yz_xx[k]);

                    t_xyz_xxy[k] = fra * t10_yz_xxy[k] - frc * t11_yz_xxy[k] + f2t * (2.0 * tk0_yz_xy[k] - 2.0 * tk1_yz_xy[k]);

                    t_xyz_xxz[k] = fra * t10_yz_xxz[k] - frc * t11_yz_xxz[k] + f2t * (2.0 * tk0_yz_xz[k] - 2.0 * tk1_yz_xz[k]);

                    t_xyz_xyy[k] = fra * t10_yz_xyy[k] - frc * t11_yz_xyy[k] + f2t * (tk0_yz_yy[k] - tk1_yz_yy[k]);

                    t_xyz_xyz[k] = fra * t10_yz_xyz[k] - frc * t11_yz_xyz[k] + f2t * (tk0_yz_yz[k] - tk1_yz_yz[k]);

                    t_xyz_xzz[k] = fra * t10_yz_xzz[k] - frc * t11_yz_xzz[k] + f2t * (tk0_yz_zz[k] - tk1_yz_zz[k]);

                    t_xyz_yyy[k] = fra * t10_yz_yyy[k] - frc * t11_yz_yyy[k];

                    t_xyz_yyz[k] = fra * t10_yz_yyz[k] - frc * t11_yz_yyz[k];

                    t_xyz_yzz[k] = fra * t10_yz_yzz[k] - frc * t11_yz_yzz[k];

                    t_xyz_zzz[k] = fra * t10_yz_zzz[k] - frc * t11_yz_zzz[k];

                    t_xzz_xxx[k] = fra * t10_zz_xxx[k] - frc * t11_zz_xxx[k] + f2t * (3.0 * tk0_zz_xx[k] - 3.0 * tk1_zz_xx[k]);

                    t_xzz_xxy[k] = fra * t10_zz_xxy[k] - frc * t11_zz_xxy[k] + f2t * (2.0 * tk0_zz_xy[k] - 2.0 * tk1_zz_xy[k]);

                    t_xzz_xxz[k] = fra * t10_zz_xxz[k] - frc * t11_zz_xxz[k] + f2t * (2.0 * tk0_zz_xz[k] - 2.0 * tk1_zz_xz[k]);

                    t_xzz_xyy[k] = fra * t10_zz_xyy[k] - frc * t11_zz_xyy[k] + f2t * (tk0_zz_yy[k] - tk1_zz_yy[k]);

                    t_xzz_xyz[k] = fra * t10_zz_xyz[k] - frc * t11_zz_xyz[k] + f2t * (tk0_zz_yz[k] - tk1_zz_yz[k]);

                    t_xzz_xzz[k] = fra * t10_zz_xzz[k] - frc * t11_zz_xzz[k] + f2t * (tk0_zz_zz[k] - tk1_zz_zz[k]);

                    t_xzz_yyy[k] = fra * t10_zz_yyy[k] - frc * t11_zz_yyy[k];

                    t_xzz_yyz[k] = fra * t10_zz_yyz[k] - frc * t11_zz_yyz[k];

                    t_xzz_yzz[k] = fra * t10_zz_yzz[k] - frc * t11_zz_yzz[k];

                    t_xzz_zzz[k] = fra * t10_zz_zzz[k] - frc * t11_zz_zzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyy_xxx[k] = fra * t10_yy_xxx[k] - frc * t11_yy_xxx[k] + f2t * (2.0 * t20_y_xxx[k] - 2.0 * t21_y_xxx[k]);

                    t_yyy_xxy[k] = fra * t10_yy_xxy[k] - frc * t11_yy_xxy[k] + f2t * (2.0 * t20_y_xxy[k] - 2.0 * t21_y_xxy[k] + tk0_yy_xx[k] - tk1_yy_xx[k]);

                    t_yyy_xxz[k] = fra * t10_yy_xxz[k] - frc * t11_yy_xxz[k] + f2t * (2.0 * t20_y_xxz[k] - 2.0 * t21_y_xxz[k]);

                    t_yyy_xyy[k] = fra * t10_yy_xyy[k] - frc * t11_yy_xyy[k] + f2t * (2.0 * t20_y_xyy[k] - 2.0 * t21_y_xyy[k] + 2.0 * tk0_yy_xy[k] - 2.0 * tk1_yy_xy[k]);

                    t_yyy_xyz[k] = fra * t10_yy_xyz[k] - frc * t11_yy_xyz[k] + f2t * (2.0 * t20_y_xyz[k] - 2.0 * t21_y_xyz[k] + tk0_yy_xz[k] - tk1_yy_xz[k]);

                    t_yyy_xzz[k] = fra * t10_yy_xzz[k] - frc * t11_yy_xzz[k] + f2t * (2.0 * t20_y_xzz[k] - 2.0 * t21_y_xzz[k]);

                    t_yyy_yyy[k] = fra * t10_yy_yyy[k] - frc * t11_yy_yyy[k] + f2t * (2.0 * t20_y_yyy[k] - 2.0 * t21_y_yyy[k] + 3.0 * tk0_yy_yy[k] - 3.0 * tk1_yy_yy[k]);

                    t_yyy_yyz[k] = fra * t10_yy_yyz[k] - frc * t11_yy_yyz[k] + f2t * (2.0 * t20_y_yyz[k] - 2.0 * t21_y_yyz[k] + 2.0 * tk0_yy_yz[k] - 2.0 * tk1_yy_yz[k]);

                    t_yyy_yzz[k] = fra * t10_yy_yzz[k] - frc * t11_yy_yzz[k] + f2t * (2.0 * t20_y_yzz[k] - 2.0 * t21_y_yzz[k] + tk0_yy_zz[k] - tk1_yy_zz[k]);

                    t_yyy_zzz[k] = fra * t10_yy_zzz[k] - frc * t11_yy_zzz[k] + f2t * (2.0 * t20_y_zzz[k] - 2.0 * t21_y_zzz[k]);

                    t_yyz_xxx[k] = fra * t10_yz_xxx[k] - frc * t11_yz_xxx[k] + f2t * (t20_z_xxx[k] - t21_z_xxx[k]);

                    t_yyz_xxy[k] = fra * t10_yz_xxy[k] - frc * t11_yz_xxy[k] + f2t * (t20_z_xxy[k] - t21_z_xxy[k] + tk0_yz_xx[k] - tk1_yz_xx[k]);

                    t_yyz_xxz[k] = fra * t10_yz_xxz[k] - frc * t11_yz_xxz[k] + f2t * (t20_z_xxz[k] - t21_z_xxz[k]);

                    t_yyz_xyy[k] = fra * t10_yz_xyy[k] - frc * t11_yz_xyy[k] + f2t * (t20_z_xyy[k] - t21_z_xyy[k] + 2.0 * tk0_yz_xy[k] - 2.0 * tk1_yz_xy[k]);

                    t_yyz_xyz[k] = fra * t10_yz_xyz[k] - frc * t11_yz_xyz[k] + f2t * (t20_z_xyz[k] - t21_z_xyz[k] + tk0_yz_xz[k] - tk1_yz_xz[k]);

                    t_yyz_xzz[k] = fra * t10_yz_xzz[k] - frc * t11_yz_xzz[k] + f2t * (t20_z_xzz[k] - t21_z_xzz[k]);

                    t_yyz_yyy[k] = fra * t10_yz_yyy[k] - frc * t11_yz_yyy[k] + f2t * (t20_z_yyy[k] - t21_z_yyy[k] + 3.0 * tk0_yz_yy[k] - 3.0 * tk1_yz_yy[k]);

                    t_yyz_yyz[k] = fra * t10_yz_yyz[k] - frc * t11_yz_yyz[k] + f2t * (t20_z_yyz[k] - t21_z_yyz[k] + 2.0 * tk0_yz_yz[k] - 2.0 * tk1_yz_yz[k]);

                    t_yyz_yzz[k] = fra * t10_yz_yzz[k] - frc * t11_yz_yzz[k] + f2t * (t20_z_yzz[k] - t21_z_yzz[k] + tk0_yz_zz[k] - tk1_yz_zz[k]);

                    t_yyz_zzz[k] = fra * t10_yz_zzz[k] - frc * t11_yz_zzz[k] + f2t * (t20_z_zzz[k] - t21_z_zzz[k]);

                    t_yzz_xxx[k] = fra * t10_zz_xxx[k] - frc * t11_zz_xxx[k];

                    t_yzz_xxy[k] = fra * t10_zz_xxy[k] - frc * t11_zz_xxy[k] + f2t * (tk0_zz_xx[k] - tk1_zz_xx[k]);

                    t_yzz_xxz[k] = fra * t10_zz_xxz[k] - frc * t11_zz_xxz[k];

                    t_yzz_xyy[k] = fra * t10_zz_xyy[k] - frc * t11_zz_xyy[k] + f2t * (2.0 * tk0_zz_xy[k] - 2.0 * tk1_zz_xy[k]);

                    t_yzz_xyz[k] = fra * t10_zz_xyz[k] - frc * t11_zz_xyz[k] + f2t * (tk0_zz_xz[k] - tk1_zz_xz[k]);

                    t_yzz_xzz[k] = fra * t10_zz_xzz[k] - frc * t11_zz_xzz[k];

                    t_yzz_yyy[k] = fra * t10_zz_yyy[k] - frc * t11_zz_yyy[k] + f2t * (3.0 * tk0_zz_yy[k] - 3.0 * tk1_zz_yy[k]);

                    t_yzz_yyz[k] = fra * t10_zz_yyz[k] - frc * t11_zz_yyz[k] + f2t * (2.0 * tk0_zz_yz[k] - 2.0 * tk1_zz_yz[k]);

                    t_yzz_yzz[k] = fra * t10_zz_yzz[k] - frc * t11_zz_yzz[k] + f2t * (tk0_zz_zz[k] - tk1_zz_zz[k]);

                    t_yzz_zzz[k] = fra * t10_zz_zzz[k] - frc * t11_zz_zzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzz_xxx[k] = fra * t10_zz_xxx[k] - frc * t11_zz_xxx[k] + f2t * (2.0 * t20_z_xxx[k] - 2.0 * t21_z_xxx[k]);

                    t_zzz_xxy[k] = fra * t10_zz_xxy[k] - frc * t11_zz_xxy[k] + f2t * (2.0 * t20_z_xxy[k] - 2.0 * t21_z_xxy[k]);

                    t_zzz_xxz[k] = fra * t10_zz_xxz[k] - frc * t11_zz_xxz[k] + f2t * (2.0 * t20_z_xxz[k] - 2.0 * t21_z_xxz[k] + tk0_zz_xx[k] - tk1_zz_xx[k]);

                    t_zzz_xyy[k] = fra * t10_zz_xyy[k] - frc * t11_zz_xyy[k] + f2t * (2.0 * t20_z_xyy[k] - 2.0 * t21_z_xyy[k]);

                    t_zzz_xyz[k] = fra * t10_zz_xyz[k] - frc * t11_zz_xyz[k] + f2t * (2.0 * t20_z_xyz[k] - 2.0 * t21_z_xyz[k] + tk0_zz_xy[k] - tk1_zz_xy[k]);

                    t_zzz_xzz[k] = fra * t10_zz_xzz[k] - frc * t11_zz_xzz[k] + f2t * (2.0 * t20_z_xzz[k] - 2.0 * t21_z_xzz[k] + 2.0 * tk0_zz_xz[k] - 2.0 * tk1_zz_xz[k]);

                    t_zzz_yyy[k] = fra * t10_zz_yyy[k] - frc * t11_zz_yyy[k] + f2t * (2.0 * t20_z_yyy[k] - 2.0 * t21_z_yyy[k]);

                    t_zzz_yyz[k] = fra * t10_zz_yyz[k] - frc * t11_zz_yyz[k] + f2t * (2.0 * t20_z_yyz[k] - 2.0 * t21_z_yyz[k] + tk0_zz_yy[k] - tk1_zz_yy[k]);

                    t_zzz_yzz[k] = fra * t10_zz_yzz[k] - frc * t11_zz_yzz[k] + f2t * (2.0 * t20_z_yzz[k] - 2.0 * t21_z_yzz[k] + 2.0 * tk0_zz_yz[k] - 2.0 * tk1_zz_yz[k]);

                    t_zzz_zzz[k] = fra * t10_zz_zzz[k] - frc * t11_zz_zzz[k] + f2t * (2.0 * t20_z_zzz[k] - 2.0 * t21_z_zzz[k] + 3.0 * tk0_zz_zz[k] - 3.0 * tk1_zz_zz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForSG(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {0, 4, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 0, 4);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 4, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 2, i + 1});

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

                // set up pointers to (S|A(0)|F)^(m) integrals

                auto t10_0_xxx = primBuffer.data(t10off + 10 * idx);

                auto t10_0_xxy = primBuffer.data(t10off + 10 * idx + 1);

                auto t10_0_xxz = primBuffer.data(t10off + 10 * idx + 2);

                auto t10_0_xyy = primBuffer.data(t10off + 10 * idx + 3);

                auto t10_0_xyz = primBuffer.data(t10off + 10 * idx + 4);

                auto t10_0_xzz = primBuffer.data(t10off + 10 * idx + 5);

                auto t10_0_yyy = primBuffer.data(t10off + 10 * idx + 6);

                auto t10_0_yyz = primBuffer.data(t10off + 10 * idx + 7);

                auto t10_0_yzz = primBuffer.data(t10off + 10 * idx + 8);

                auto t10_0_zzz = primBuffer.data(t10off + 10 * idx + 9);

                // set up pointers to (S|A(0)|F)^(m+1) integrals

                auto t11_0_xxx = primBuffer.data(t11off + 10 * idx);

                auto t11_0_xxy = primBuffer.data(t11off + 10 * idx + 1);

                auto t11_0_xxz = primBuffer.data(t11off + 10 * idx + 2);

                auto t11_0_xyy = primBuffer.data(t11off + 10 * idx + 3);

                auto t11_0_xyz = primBuffer.data(t11off + 10 * idx + 4);

                auto t11_0_xzz = primBuffer.data(t11off + 10 * idx + 5);

                auto t11_0_yyy = primBuffer.data(t11off + 10 * idx + 6);

                auto t11_0_yyz = primBuffer.data(t11off + 10 * idx + 7);

                auto t11_0_yzz = primBuffer.data(t11off + 10 * idx + 8);

                auto t11_0_zzz = primBuffer.data(t11off + 10 * idx + 9);

                // set up pointers to (S|A(0)|G)^(m) integrals

                auto t_0_xxxx = primBuffer.data(toff + 15 * idx);

                auto t_0_xxxy = primBuffer.data(toff + 15 * idx + 1);

                auto t_0_xxxz = primBuffer.data(toff + 15 * idx + 2);

                auto t_0_xxyy = primBuffer.data(toff + 15 * idx + 3);

                auto t_0_xxyz = primBuffer.data(toff + 15 * idx + 4);

                auto t_0_xxzz = primBuffer.data(toff + 15 * idx + 5);

                auto t_0_xyyy = primBuffer.data(toff + 15 * idx + 6);

                auto t_0_xyyz = primBuffer.data(toff + 15 * idx + 7);

                auto t_0_xyzz = primBuffer.data(toff + 15 * idx + 8);

                auto t_0_xzzz = primBuffer.data(toff + 15 * idx + 9);

                auto t_0_yyyy = primBuffer.data(toff + 15 * idx + 10);

                auto t_0_yyyz = primBuffer.data(toff + 15 * idx + 11);

                auto t_0_yyzz = primBuffer.data(toff + 15 * idx + 12);

                auto t_0_yzzz = primBuffer.data(toff + 15 * idx + 13);

                auto t_0_zzzz = primBuffer.data(toff + 15 * idx + 14);

                #pragma omp simd aligned(fx, pbx, pby, pbz, t20_0_xx, t20_0_xy,\
                                         t20_0_xz, t20_0_yy, t20_0_yz, t20_0_zz,\
                                         t21_0_xx, t21_0_xy, t21_0_xz, t21_0_yy,\
                                         t21_0_yz, t21_0_zz, t10_0_xxx, t10_0_xxy,\
                                         t10_0_xxz, t10_0_xyy, t10_0_xyz, t10_0_xzz,\
                                         t10_0_yyy, t10_0_yyz, t10_0_yzz, t10_0_zzz,\
                                         t11_0_xxx, t11_0_xxy, t11_0_xxz, t11_0_xyy,\
                                         t11_0_xyz, t11_0_xzz, t11_0_yyy, t11_0_yyz,\
                                         t11_0_yzz, t11_0_zzz, t_0_xxxx, t_0_xxxy,\
                                         t_0_xxxz, t_0_xxyy, t_0_xxyz, t_0_xxzz,\
                                         t_0_xyyy, t_0_xyyz, t_0_xyzz, t_0_xzzz,\
                                         t_0_yyyy, t_0_yyyz, t_0_yyzz, t_0_yzzz,\
                                         t_0_zzzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pbx[k];

                    double frc = pcx[k];

                    t_0_xxxx[k] = fra * t10_0_xxx[k] - frc * t11_0_xxx[k] + f2t * (3.0 * t20_0_xx[k] - 3.0 * t21_0_xx[k]);

                    t_0_xxxy[k] = fra * t10_0_xxy[k] - frc * t11_0_xxy[k] + f2t * (2.0 * t20_0_xy[k] - 2.0 * t21_0_xy[k]);

                    t_0_xxxz[k] = fra * t10_0_xxz[k] - frc * t11_0_xxz[k] + f2t * (2.0 * t20_0_xz[k] - 2.0 * t21_0_xz[k]);

                    t_0_xxyy[k] = fra * t10_0_xyy[k] - frc * t11_0_xyy[k] + f2t * (t20_0_yy[k] - t21_0_yy[k]);

                    t_0_xxyz[k] = fra * t10_0_xyz[k] - frc * t11_0_xyz[k] + f2t * (t20_0_yz[k] - t21_0_yz[k]);

                    t_0_xxzz[k] = fra * t10_0_xzz[k] - frc * t11_0_xzz[k] + f2t * (t20_0_zz[k] - t21_0_zz[k]);

                    t_0_xyyy[k] = fra * t10_0_yyy[k] - frc * t11_0_yyy[k];

                    t_0_xyyz[k] = fra * t10_0_yyz[k] - frc * t11_0_yyz[k];

                    t_0_xyzz[k] = fra * t10_0_yzz[k] - frc * t11_0_yzz[k];

                    t_0_xzzz[k] = fra * t10_0_zzz[k] - frc * t11_0_zzz[k];

                    // leading y component

                    fra = pby[k];

                    frc = pcy[k];

                    t_0_yyyy[k] = fra * t10_0_yyy[k] - frc * t11_0_yyy[k] + f2t * (3.0 * t20_0_yy[k] - 3.0 * t21_0_yy[k]);

                    t_0_yyyz[k] = fra * t10_0_yyz[k] - frc * t11_0_yyz[k] + f2t * (2.0 * t20_0_yz[k] - 2.0 * t21_0_yz[k]);

                    t_0_yyzz[k] = fra * t10_0_yzz[k] - frc * t11_0_yzz[k] + f2t * (t20_0_zz[k] - t21_0_zz[k]);

                    t_0_yzzz[k] = fra * t10_0_zzz[k] - frc * t11_0_zzz[k];

                    // leading z component

                    t_0_zzzz[k] = pbz[k] * t10_0_zzz[k] - pcz[k] * t11_0_zzz[k] + f2t * (3.0 * t20_0_zz[k] - 3.0 * t21_0_zz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForGS(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {4, 0, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 4, 0);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {4, 0, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 0, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 0, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 0, i + 1});

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

                // set up pointers to (D|A(0)|S)^(m) integrals

                auto t20_xx_0 = primBuffer.data(t20off + 6 * idx);

                auto t20_xy_0 = primBuffer.data(t20off + 6 * idx + 1);

                auto t20_xz_0 = primBuffer.data(t20off + 6 * idx + 2);

                auto t20_yy_0 = primBuffer.data(t20off + 6 * idx + 3);

                auto t20_yz_0 = primBuffer.data(t20off + 6 * idx + 4);

                auto t20_zz_0 = primBuffer.data(t20off + 6 * idx + 5);

                // set up pointers to (D|A(0)|S)^(m+1) integrals

                auto t21_xx_0 = primBuffer.data(t21off + 6 * idx);

                auto t21_xy_0 = primBuffer.data(t21off + 6 * idx + 1);

                auto t21_xz_0 = primBuffer.data(t21off + 6 * idx + 2);

                auto t21_yy_0 = primBuffer.data(t21off + 6 * idx + 3);

                auto t21_yz_0 = primBuffer.data(t21off + 6 * idx + 4);

                auto t21_zz_0 = primBuffer.data(t21off + 6 * idx + 5);

                // set up pointers to (F|A(0)|S)^(m) integrals

                auto t10_xxx_0 = primBuffer.data(t10off + 10 * idx);

                auto t10_xxy_0 = primBuffer.data(t10off + 10 * idx + 1);

                auto t10_xxz_0 = primBuffer.data(t10off + 10 * idx + 2);

                auto t10_xyy_0 = primBuffer.data(t10off + 10 * idx + 3);

                auto t10_xyz_0 = primBuffer.data(t10off + 10 * idx + 4);

                auto t10_xzz_0 = primBuffer.data(t10off + 10 * idx + 5);

                auto t10_yyy_0 = primBuffer.data(t10off + 10 * idx + 6);

                auto t10_yyz_0 = primBuffer.data(t10off + 10 * idx + 7);

                auto t10_yzz_0 = primBuffer.data(t10off + 10 * idx + 8);

                auto t10_zzz_0 = primBuffer.data(t10off + 10 * idx + 9);

                // set up pointers to (F|A(0)|S)^(m+1) integrals

                auto t11_xxx_0 = primBuffer.data(t11off + 10 * idx);

                auto t11_xxy_0 = primBuffer.data(t11off + 10 * idx + 1);

                auto t11_xxz_0 = primBuffer.data(t11off + 10 * idx + 2);

                auto t11_xyy_0 = primBuffer.data(t11off + 10 * idx + 3);

                auto t11_xyz_0 = primBuffer.data(t11off + 10 * idx + 4);

                auto t11_xzz_0 = primBuffer.data(t11off + 10 * idx + 5);

                auto t11_yyy_0 = primBuffer.data(t11off + 10 * idx + 6);

                auto t11_yyz_0 = primBuffer.data(t11off + 10 * idx + 7);

                auto t11_yzz_0 = primBuffer.data(t11off + 10 * idx + 8);

                auto t11_zzz_0 = primBuffer.data(t11off + 10 * idx + 9);

                // set up pointers to (G|A(0)|S)^(m) integrals

                auto t_xxxx_0 = primBuffer.data(toff + 15 * idx);

                auto t_xxxy_0 = primBuffer.data(toff + 15 * idx + 1);

                auto t_xxxz_0 = primBuffer.data(toff + 15 * idx + 2);

                auto t_xxyy_0 = primBuffer.data(toff + 15 * idx + 3);

                auto t_xxyz_0 = primBuffer.data(toff + 15 * idx + 4);

                auto t_xxzz_0 = primBuffer.data(toff + 15 * idx + 5);

                auto t_xyyy_0 = primBuffer.data(toff + 15 * idx + 6);

                auto t_xyyz_0 = primBuffer.data(toff + 15 * idx + 7);

                auto t_xyzz_0 = primBuffer.data(toff + 15 * idx + 8);

                auto t_xzzz_0 = primBuffer.data(toff + 15 * idx + 9);

                auto t_yyyy_0 = primBuffer.data(toff + 15 * idx + 10);

                auto t_yyyz_0 = primBuffer.data(toff + 15 * idx + 11);

                auto t_yyzz_0 = primBuffer.data(toff + 15 * idx + 12);

                auto t_yzzz_0 = primBuffer.data(toff + 15 * idx + 13);

                auto t_zzzz_0 = primBuffer.data(toff + 15 * idx + 14);

                #pragma omp simd aligned(fx, pax, pay, paz, t20_xx_0, t20_xy_0,\
                                         t20_xz_0, t20_yy_0, t20_yz_0, t20_zz_0,\
                                         t21_xx_0, t21_xy_0, t21_xz_0, t21_yy_0,\
                                         t21_yz_0, t21_zz_0, t10_xxx_0, t10_xxy_0,\
                                         t10_xxz_0, t10_xyy_0, t10_xyz_0, t10_xzz_0,\
                                         t10_yyy_0, t10_yyz_0, t10_yzz_0, t10_zzz_0,\
                                         t11_xxx_0, t11_xxy_0, t11_xxz_0, t11_xyy_0,\
                                         t11_xyz_0, t11_xzz_0, t11_yyy_0, t11_yyz_0,\
                                         t11_yzz_0, t11_zzz_0, t_xxxx_0, t_xxxy_0,\
                                         t_xxxz_0, t_xxyy_0, t_xxyz_0, t_xxzz_0,\
                                         t_xyyy_0, t_xyyz_0, t_xyzz_0, t_xzzz_0,\
                                         t_yyyy_0, t_yyyz_0, t_yyzz_0, t_yzzz_0,\
                                         t_zzzz_0, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxxx_0[k] = fra * t10_xxx_0[k] - frc * t11_xxx_0[k] + f2t * (3.0 * t20_xx_0[k] - 3.0 * t21_xx_0[k]);

                    t_xxxy_0[k] = fra * t10_xxy_0[k] - frc * t11_xxy_0[k] + f2t * (2.0 * t20_xy_0[k] - 2.0 * t21_xy_0[k]);

                    t_xxxz_0[k] = fra * t10_xxz_0[k] - frc * t11_xxz_0[k] + f2t * (2.0 * t20_xz_0[k] - 2.0 * t21_xz_0[k]);

                    t_xxyy_0[k] = fra * t10_xyy_0[k] - frc * t11_xyy_0[k] + f2t * (t20_yy_0[k] - t21_yy_0[k]);

                    t_xxyz_0[k] = fra * t10_xyz_0[k] - frc * t11_xyz_0[k] + f2t * (t20_yz_0[k] - t21_yz_0[k]);

                    t_xxzz_0[k] = fra * t10_xzz_0[k] - frc * t11_xzz_0[k] + f2t * (t20_zz_0[k] - t21_zz_0[k]);

                    t_xyyy_0[k] = fra * t10_yyy_0[k] - frc * t11_yyy_0[k];

                    t_xyyz_0[k] = fra * t10_yyz_0[k] - frc * t11_yyz_0[k];

                    t_xyzz_0[k] = fra * t10_yzz_0[k] - frc * t11_yzz_0[k];

                    t_xzzz_0[k] = fra * t10_zzz_0[k] - frc * t11_zzz_0[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyyy_0[k] = fra * t10_yyy_0[k] - frc * t11_yyy_0[k] + f2t * (3.0 * t20_yy_0[k] - 3.0 * t21_yy_0[k]);

                    t_yyyz_0[k] = fra * t10_yyz_0[k] - frc * t11_yyz_0[k] + f2t * (2.0 * t20_yz_0[k] - 2.0 * t21_yz_0[k]);

                    t_yyzz_0[k] = fra * t10_yzz_0[k] - frc * t11_yzz_0[k] + f2t * (t20_zz_0[k] - t21_zz_0[k]);

                    t_yzzz_0[k] = fra * t10_zzz_0[k] - frc * t11_zzz_0[k];

                    // leading z component

                    t_zzzz_0[k] = paz[k] * t10_zzz_0[k] - pcz[k] * t11_zzz_0[k] + f2t * (3.0 * t20_zz_0[k] - 3.0 * t21_zz_0[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForPG(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {1, 4, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 4);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 4, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 4, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 4, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 3, i + 1});

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

                // set up pointers to (S|A(0)|F)^(m) integrals

                auto tk0_0_xxx = primBuffer.data(tk0off + 10 * idx);

                auto tk0_0_xxy = primBuffer.data(tk0off + 10 * idx + 1);

                auto tk0_0_xxz = primBuffer.data(tk0off + 10 * idx + 2);

                auto tk0_0_xyy = primBuffer.data(tk0off + 10 * idx + 3);

                auto tk0_0_xyz = primBuffer.data(tk0off + 10 * idx + 4);

                auto tk0_0_xzz = primBuffer.data(tk0off + 10 * idx + 5);

                auto tk0_0_yyy = primBuffer.data(tk0off + 10 * idx + 6);

                auto tk0_0_yyz = primBuffer.data(tk0off + 10 * idx + 7);

                auto tk0_0_yzz = primBuffer.data(tk0off + 10 * idx + 8);

                auto tk0_0_zzz = primBuffer.data(tk0off + 10 * idx + 9);

                // set up pointers to (S|A(0)|F)^(m+1) integrals

                auto tk1_0_xxx = primBuffer.data(tk1off + 10 * idx);

                auto tk1_0_xxy = primBuffer.data(tk1off + 10 * idx + 1);

                auto tk1_0_xxz = primBuffer.data(tk1off + 10 * idx + 2);

                auto tk1_0_xyy = primBuffer.data(tk1off + 10 * idx + 3);

                auto tk1_0_xyz = primBuffer.data(tk1off + 10 * idx + 4);

                auto tk1_0_xzz = primBuffer.data(tk1off + 10 * idx + 5);

                auto tk1_0_yyy = primBuffer.data(tk1off + 10 * idx + 6);

                auto tk1_0_yyz = primBuffer.data(tk1off + 10 * idx + 7);

                auto tk1_0_yzz = primBuffer.data(tk1off + 10 * idx + 8);

                auto tk1_0_zzz = primBuffer.data(tk1off + 10 * idx + 9);

                // set up pointers to (S|A(0)|G)^(m) integrals

                auto t10_0_xxxx = primBuffer.data(t10off + 15 * idx);

                auto t10_0_xxxy = primBuffer.data(t10off + 15 * idx + 1);

                auto t10_0_xxxz = primBuffer.data(t10off + 15 * idx + 2);

                auto t10_0_xxyy = primBuffer.data(t10off + 15 * idx + 3);

                auto t10_0_xxyz = primBuffer.data(t10off + 15 * idx + 4);

                auto t10_0_xxzz = primBuffer.data(t10off + 15 * idx + 5);

                auto t10_0_xyyy = primBuffer.data(t10off + 15 * idx + 6);

                auto t10_0_xyyz = primBuffer.data(t10off + 15 * idx + 7);

                auto t10_0_xyzz = primBuffer.data(t10off + 15 * idx + 8);

                auto t10_0_xzzz = primBuffer.data(t10off + 15 * idx + 9);

                auto t10_0_yyyy = primBuffer.data(t10off + 15 * idx + 10);

                auto t10_0_yyyz = primBuffer.data(t10off + 15 * idx + 11);

                auto t10_0_yyzz = primBuffer.data(t10off + 15 * idx + 12);

                auto t10_0_yzzz = primBuffer.data(t10off + 15 * idx + 13);

                auto t10_0_zzzz = primBuffer.data(t10off + 15 * idx + 14);

                // set up pointers to (S|A(0)|G)^(m+1) integrals

                auto t11_0_xxxx = primBuffer.data(t11off + 15 * idx);

                auto t11_0_xxxy = primBuffer.data(t11off + 15 * idx + 1);

                auto t11_0_xxxz = primBuffer.data(t11off + 15 * idx + 2);

                auto t11_0_xxyy = primBuffer.data(t11off + 15 * idx + 3);

                auto t11_0_xxyz = primBuffer.data(t11off + 15 * idx + 4);

                auto t11_0_xxzz = primBuffer.data(t11off + 15 * idx + 5);

                auto t11_0_xyyy = primBuffer.data(t11off + 15 * idx + 6);

                auto t11_0_xyyz = primBuffer.data(t11off + 15 * idx + 7);

                auto t11_0_xyzz = primBuffer.data(t11off + 15 * idx + 8);

                auto t11_0_xzzz = primBuffer.data(t11off + 15 * idx + 9);

                auto t11_0_yyyy = primBuffer.data(t11off + 15 * idx + 10);

                auto t11_0_yyyz = primBuffer.data(t11off + 15 * idx + 11);

                auto t11_0_yyzz = primBuffer.data(t11off + 15 * idx + 12);

                auto t11_0_yzzz = primBuffer.data(t11off + 15 * idx + 13);

                auto t11_0_zzzz = primBuffer.data(t11off + 15 * idx + 14);

                // set up pointers to (P|A(0)|G)^(m) integrals

                auto t_x_xxxx = primBuffer.data(toff + 45 * idx);

                auto t_x_xxxy = primBuffer.data(toff + 45 * idx + 1);

                auto t_x_xxxz = primBuffer.data(toff + 45 * idx + 2);

                auto t_x_xxyy = primBuffer.data(toff + 45 * idx + 3);

                auto t_x_xxyz = primBuffer.data(toff + 45 * idx + 4);

                auto t_x_xxzz = primBuffer.data(toff + 45 * idx + 5);

                auto t_x_xyyy = primBuffer.data(toff + 45 * idx + 6);

                auto t_x_xyyz = primBuffer.data(toff + 45 * idx + 7);

                auto t_x_xyzz = primBuffer.data(toff + 45 * idx + 8);

                auto t_x_xzzz = primBuffer.data(toff + 45 * idx + 9);

                auto t_x_yyyy = primBuffer.data(toff + 45 * idx + 10);

                auto t_x_yyyz = primBuffer.data(toff + 45 * idx + 11);

                auto t_x_yyzz = primBuffer.data(toff + 45 * idx + 12);

                auto t_x_yzzz = primBuffer.data(toff + 45 * idx + 13);

                auto t_x_zzzz = primBuffer.data(toff + 45 * idx + 14);

                auto t_y_xxxx = primBuffer.data(toff + 45 * idx + 15);

                auto t_y_xxxy = primBuffer.data(toff + 45 * idx + 16);

                auto t_y_xxxz = primBuffer.data(toff + 45 * idx + 17);

                auto t_y_xxyy = primBuffer.data(toff + 45 * idx + 18);

                auto t_y_xxyz = primBuffer.data(toff + 45 * idx + 19);

                auto t_y_xxzz = primBuffer.data(toff + 45 * idx + 20);

                auto t_y_xyyy = primBuffer.data(toff + 45 * idx + 21);

                auto t_y_xyyz = primBuffer.data(toff + 45 * idx + 22);

                auto t_y_xyzz = primBuffer.data(toff + 45 * idx + 23);

                auto t_y_xzzz = primBuffer.data(toff + 45 * idx + 24);

                auto t_y_yyyy = primBuffer.data(toff + 45 * idx + 25);

                auto t_y_yyyz = primBuffer.data(toff + 45 * idx + 26);

                auto t_y_yyzz = primBuffer.data(toff + 45 * idx + 27);

                auto t_y_yzzz = primBuffer.data(toff + 45 * idx + 28);

                auto t_y_zzzz = primBuffer.data(toff + 45 * idx + 29);

                auto t_z_xxxx = primBuffer.data(toff + 45 * idx + 30);

                auto t_z_xxxy = primBuffer.data(toff + 45 * idx + 31);

                auto t_z_xxxz = primBuffer.data(toff + 45 * idx + 32);

                auto t_z_xxyy = primBuffer.data(toff + 45 * idx + 33);

                auto t_z_xxyz = primBuffer.data(toff + 45 * idx + 34);

                auto t_z_xxzz = primBuffer.data(toff + 45 * idx + 35);

                auto t_z_xyyy = primBuffer.data(toff + 45 * idx + 36);

                auto t_z_xyyz = primBuffer.data(toff + 45 * idx + 37);

                auto t_z_xyzz = primBuffer.data(toff + 45 * idx + 38);

                auto t_z_xzzz = primBuffer.data(toff + 45 * idx + 39);

                auto t_z_yyyy = primBuffer.data(toff + 45 * idx + 40);

                auto t_z_yyyz = primBuffer.data(toff + 45 * idx + 41);

                auto t_z_yyzz = primBuffer.data(toff + 45 * idx + 42);

                auto t_z_yzzz = primBuffer.data(toff + 45 * idx + 43);

                auto t_z_zzzz = primBuffer.data(toff + 45 * idx + 44);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_0_xxx, tk0_0_xxy,\
                                         tk0_0_xxz, tk0_0_xyy, tk0_0_xyz, tk0_0_xzz,\
                                         tk0_0_yyy, tk0_0_yyz, tk0_0_yzz, tk0_0_zzz,\
                                         tk1_0_xxx, tk1_0_xxy, tk1_0_xxz, tk1_0_xyy,\
                                         tk1_0_xyz, tk1_0_xzz, tk1_0_yyy, tk1_0_yyz,\
                                         tk1_0_yzz, tk1_0_zzz, t10_0_xxxx, t10_0_xxxy,\
                                         t10_0_xxxz, t10_0_xxyy, t10_0_xxyz, t10_0_xxzz,\
                                         t10_0_xyyy, t10_0_xyyz, t10_0_xyzz, t10_0_xzzz,\
                                         t10_0_yyyy, t10_0_yyyz, t10_0_yyzz, t10_0_yzzz,\
                                         t10_0_zzzz, t11_0_xxxx, t11_0_xxxy, t11_0_xxxz,\
                                         t11_0_xxyy, t11_0_xxyz, t11_0_xxzz, t11_0_xyyy,\
                                         t11_0_xyyz, t11_0_xyzz, t11_0_xzzz, t11_0_yyyy,\
                                         t11_0_yyyz, t11_0_yyzz, t11_0_yzzz, t11_0_zzzz,\
                                         t_x_xxxx, t_x_xxxy, t_x_xxxz, t_x_xxyy,\
                                         t_x_xxyz, t_x_xxzz, t_x_xyyy, t_x_xyyz,\
                                         t_x_xyzz, t_x_xzzz, t_x_yyyy, t_x_yyyz,\
                                         t_x_yyzz, t_x_yzzz, t_x_zzzz, t_y_xxxx,\
                                         t_y_xxxy, t_y_xxxz, t_y_xxyy, t_y_xxyz,\
                                         t_y_xxzz, t_y_xyyy, t_y_xyyz, t_y_xyzz,\
                                         t_y_xzzz, t_y_yyyy, t_y_yyyz, t_y_yyzz,\
                                         t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy,\
                                         t_z_xxxz, t_z_xxyy, t_z_xxyz, t_z_xxzz,\
                                         t_z_xyyy, t_z_xyyz, t_z_xyzz, t_z_xzzz,\
                                         t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz,\
                                         t_z_zzzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_x_xxxx[k] = fra * t10_0_xxxx[k] - frc * t11_0_xxxx[k] + f2t * (4.0 * tk0_0_xxx[k] - 4.0 * tk1_0_xxx[k]);

                    t_x_xxxy[k] = fra * t10_0_xxxy[k] - frc * t11_0_xxxy[k] + f2t * (3.0 * tk0_0_xxy[k] - 3.0 * tk1_0_xxy[k]);

                    t_x_xxxz[k] = fra * t10_0_xxxz[k] - frc * t11_0_xxxz[k] + f2t * (3.0 * tk0_0_xxz[k] - 3.0 * tk1_0_xxz[k]);

                    t_x_xxyy[k] = fra * t10_0_xxyy[k] - frc * t11_0_xxyy[k] + f2t * (2.0 * tk0_0_xyy[k] - 2.0 * tk1_0_xyy[k]);

                    t_x_xxyz[k] = fra * t10_0_xxyz[k] - frc * t11_0_xxyz[k] + f2t * (2.0 * tk0_0_xyz[k] - 2.0 * tk1_0_xyz[k]);

                    t_x_xxzz[k] = fra * t10_0_xxzz[k] - frc * t11_0_xxzz[k] + f2t * (2.0 * tk0_0_xzz[k] - 2.0 * tk1_0_xzz[k]);

                    t_x_xyyy[k] = fra * t10_0_xyyy[k] - frc * t11_0_xyyy[k] + f2t * (tk0_0_yyy[k] - tk1_0_yyy[k]);

                    t_x_xyyz[k] = fra * t10_0_xyyz[k] - frc * t11_0_xyyz[k] + f2t * (tk0_0_yyz[k] - tk1_0_yyz[k]);

                    t_x_xyzz[k] = fra * t10_0_xyzz[k] - frc * t11_0_xyzz[k] + f2t * (tk0_0_yzz[k] - tk1_0_yzz[k]);

                    t_x_xzzz[k] = fra * t10_0_xzzz[k] - frc * t11_0_xzzz[k] + f2t * (tk0_0_zzz[k] - tk1_0_zzz[k]);

                    t_x_yyyy[k] = fra * t10_0_yyyy[k] - frc * t11_0_yyyy[k];

                    t_x_yyyz[k] = fra * t10_0_yyyz[k] - frc * t11_0_yyyz[k];

                    t_x_yyzz[k] = fra * t10_0_yyzz[k] - frc * t11_0_yyzz[k];

                    t_x_yzzz[k] = fra * t10_0_yzzz[k] - frc * t11_0_yzzz[k];

                    t_x_zzzz[k] = fra * t10_0_zzzz[k] - frc * t11_0_zzzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_y_xxxx[k] = fra * t10_0_xxxx[k] - frc * t11_0_xxxx[k];

                    t_y_xxxy[k] = fra * t10_0_xxxy[k] - frc * t11_0_xxxy[k] + f2t * (tk0_0_xxx[k] - tk1_0_xxx[k]);

                    t_y_xxxz[k] = fra * t10_0_xxxz[k] - frc * t11_0_xxxz[k];

                    t_y_xxyy[k] = fra * t10_0_xxyy[k] - frc * t11_0_xxyy[k] + f2t * (2.0 * tk0_0_xxy[k] - 2.0 * tk1_0_xxy[k]);

                    t_y_xxyz[k] = fra * t10_0_xxyz[k] - frc * t11_0_xxyz[k] + f2t * (tk0_0_xxz[k] - tk1_0_xxz[k]);

                    t_y_xxzz[k] = fra * t10_0_xxzz[k] - frc * t11_0_xxzz[k];

                    t_y_xyyy[k] = fra * t10_0_xyyy[k] - frc * t11_0_xyyy[k] + f2t * (3.0 * tk0_0_xyy[k] - 3.0 * tk1_0_xyy[k]);

                    t_y_xyyz[k] = fra * t10_0_xyyz[k] - frc * t11_0_xyyz[k] + f2t * (2.0 * tk0_0_xyz[k] - 2.0 * tk1_0_xyz[k]);

                    t_y_xyzz[k] = fra * t10_0_xyzz[k] - frc * t11_0_xyzz[k] + f2t * (tk0_0_xzz[k] - tk1_0_xzz[k]);

                    t_y_xzzz[k] = fra * t10_0_xzzz[k] - frc * t11_0_xzzz[k];

                    t_y_yyyy[k] = fra * t10_0_yyyy[k] - frc * t11_0_yyyy[k] + f2t * (4.0 * tk0_0_yyy[k] - 4.0 * tk1_0_yyy[k]);

                    t_y_yyyz[k] = fra * t10_0_yyyz[k] - frc * t11_0_yyyz[k] + f2t * (3.0 * tk0_0_yyz[k] - 3.0 * tk1_0_yyz[k]);

                    t_y_yyzz[k] = fra * t10_0_yyzz[k] - frc * t11_0_yyzz[k] + f2t * (2.0 * tk0_0_yzz[k] - 2.0 * tk1_0_yzz[k]);

                    t_y_yzzz[k] = fra * t10_0_yzzz[k] - frc * t11_0_yzzz[k] + f2t * (tk0_0_zzz[k] - tk1_0_zzz[k]);

                    t_y_zzzz[k] = fra * t10_0_zzzz[k] - frc * t11_0_zzzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_z_xxxx[k] = fra * t10_0_xxxx[k] - frc * t11_0_xxxx[k];

                    t_z_xxxy[k] = fra * t10_0_xxxy[k] - frc * t11_0_xxxy[k];

                    t_z_xxxz[k] = fra * t10_0_xxxz[k] - frc * t11_0_xxxz[k] + f2t * (tk0_0_xxx[k] - tk1_0_xxx[k]);

                    t_z_xxyy[k] = fra * t10_0_xxyy[k] - frc * t11_0_xxyy[k];

                    t_z_xxyz[k] = fra * t10_0_xxyz[k] - frc * t11_0_xxyz[k] + f2t * (tk0_0_xxy[k] - tk1_0_xxy[k]);

                    t_z_xxzz[k] = fra * t10_0_xxzz[k] - frc * t11_0_xxzz[k] + f2t * (2.0 * tk0_0_xxz[k] - 2.0 * tk1_0_xxz[k]);

                    t_z_xyyy[k] = fra * t10_0_xyyy[k] - frc * t11_0_xyyy[k];

                    t_z_xyyz[k] = fra * t10_0_xyyz[k] - frc * t11_0_xyyz[k] + f2t * (tk0_0_xyy[k] - tk1_0_xyy[k]);

                    t_z_xyzz[k] = fra * t10_0_xyzz[k] - frc * t11_0_xyzz[k] + f2t * (2.0 * tk0_0_xyz[k] - 2.0 * tk1_0_xyz[k]);

                    t_z_xzzz[k] = fra * t10_0_xzzz[k] - frc * t11_0_xzzz[k] + f2t * (3.0 * tk0_0_xzz[k] - 3.0 * tk1_0_xzz[k]);

                    t_z_yyyy[k] = fra * t10_0_yyyy[k] - frc * t11_0_yyyy[k];

                    t_z_yyyz[k] = fra * t10_0_yyyz[k] - frc * t11_0_yyyz[k] + f2t * (tk0_0_yyy[k] - tk1_0_yyy[k]);

                    t_z_yyzz[k] = fra * t10_0_yyzz[k] - frc * t11_0_yyzz[k] + f2t * (2.0 * tk0_0_yyz[k] - 2.0 * tk1_0_yyz[k]);

                    t_z_yzzz[k] = fra * t10_0_yzzz[k] - frc * t11_0_yzzz[k] + f2t * (3.0 * tk0_0_yzz[k] - 3.0 * tk1_0_yzz[k]);

                    t_z_zzzz[k] = fra * t10_0_zzzz[k] - frc * t11_0_zzzz[k] + f2t * (4.0 * tk0_0_zzz[k] - 4.0 * tk1_0_zzz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForGP(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {4, 1, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 4, 1);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {4, 1, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 1, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 1, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 1, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 0, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 0, i + 1});

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

                // set up pointers to (F|A(0)|S)^(m) integrals

                auto tk0_xxx_0 = primBuffer.data(tk0off + 10 * idx);

                auto tk0_xxy_0 = primBuffer.data(tk0off + 10 * idx + 1);

                auto tk0_xxz_0 = primBuffer.data(tk0off + 10 * idx + 2);

                auto tk0_xyy_0 = primBuffer.data(tk0off + 10 * idx + 3);

                auto tk0_xyz_0 = primBuffer.data(tk0off + 10 * idx + 4);

                auto tk0_xzz_0 = primBuffer.data(tk0off + 10 * idx + 5);

                auto tk0_yyy_0 = primBuffer.data(tk0off + 10 * idx + 6);

                auto tk0_yyz_0 = primBuffer.data(tk0off + 10 * idx + 7);

                auto tk0_yzz_0 = primBuffer.data(tk0off + 10 * idx + 8);

                auto tk0_zzz_0 = primBuffer.data(tk0off + 10 * idx + 9);

                // set up pointers to (F|A(0)|S)^(m+1) integrals

                auto tk1_xxx_0 = primBuffer.data(tk1off + 10 * idx);

                auto tk1_xxy_0 = primBuffer.data(tk1off + 10 * idx + 1);

                auto tk1_xxz_0 = primBuffer.data(tk1off + 10 * idx + 2);

                auto tk1_xyy_0 = primBuffer.data(tk1off + 10 * idx + 3);

                auto tk1_xyz_0 = primBuffer.data(tk1off + 10 * idx + 4);

                auto tk1_xzz_0 = primBuffer.data(tk1off + 10 * idx + 5);

                auto tk1_yyy_0 = primBuffer.data(tk1off + 10 * idx + 6);

                auto tk1_yyz_0 = primBuffer.data(tk1off + 10 * idx + 7);

                auto tk1_yzz_0 = primBuffer.data(tk1off + 10 * idx + 8);

                auto tk1_zzz_0 = primBuffer.data(tk1off + 10 * idx + 9);

                // set up pointers to (D|A(0)|P)^(m) integrals

                auto t20_xx_x = primBuffer.data(t20off + 18 * idx);

                auto t20_xx_y = primBuffer.data(t20off + 18 * idx + 1);

                auto t20_xx_z = primBuffer.data(t20off + 18 * idx + 2);

                auto t20_xy_x = primBuffer.data(t20off + 18 * idx + 3);

                auto t20_xy_y = primBuffer.data(t20off + 18 * idx + 4);

                auto t20_xy_z = primBuffer.data(t20off + 18 * idx + 5);

                auto t20_xz_x = primBuffer.data(t20off + 18 * idx + 6);

                auto t20_xz_y = primBuffer.data(t20off + 18 * idx + 7);

                auto t20_xz_z = primBuffer.data(t20off + 18 * idx + 8);

                auto t20_yy_x = primBuffer.data(t20off + 18 * idx + 9);

                auto t20_yy_y = primBuffer.data(t20off + 18 * idx + 10);

                auto t20_yy_z = primBuffer.data(t20off + 18 * idx + 11);

                auto t20_yz_x = primBuffer.data(t20off + 18 * idx + 12);

                auto t20_yz_y = primBuffer.data(t20off + 18 * idx + 13);

                auto t20_yz_z = primBuffer.data(t20off + 18 * idx + 14);

                auto t20_zz_x = primBuffer.data(t20off + 18 * idx + 15);

                auto t20_zz_y = primBuffer.data(t20off + 18 * idx + 16);

                auto t20_zz_z = primBuffer.data(t20off + 18 * idx + 17);

                // set up pointers to (D|A(0)|P)^(m+1) integrals

                auto t21_xx_x = primBuffer.data(t21off + 18 * idx);

                auto t21_xx_y = primBuffer.data(t21off + 18 * idx + 1);

                auto t21_xx_z = primBuffer.data(t21off + 18 * idx + 2);

                auto t21_xy_x = primBuffer.data(t21off + 18 * idx + 3);

                auto t21_xy_y = primBuffer.data(t21off + 18 * idx + 4);

                auto t21_xy_z = primBuffer.data(t21off + 18 * idx + 5);

                auto t21_xz_x = primBuffer.data(t21off + 18 * idx + 6);

                auto t21_xz_y = primBuffer.data(t21off + 18 * idx + 7);

                auto t21_xz_z = primBuffer.data(t21off + 18 * idx + 8);

                auto t21_yy_x = primBuffer.data(t21off + 18 * idx + 9);

                auto t21_yy_y = primBuffer.data(t21off + 18 * idx + 10);

                auto t21_yy_z = primBuffer.data(t21off + 18 * idx + 11);

                auto t21_yz_x = primBuffer.data(t21off + 18 * idx + 12);

                auto t21_yz_y = primBuffer.data(t21off + 18 * idx + 13);

                auto t21_yz_z = primBuffer.data(t21off + 18 * idx + 14);

                auto t21_zz_x = primBuffer.data(t21off + 18 * idx + 15);

                auto t21_zz_y = primBuffer.data(t21off + 18 * idx + 16);

                auto t21_zz_z = primBuffer.data(t21off + 18 * idx + 17);

                // set up pointers to (F|A(0)|P)^(m) integrals

                auto t10_xxx_x = primBuffer.data(t10off + 30 * idx);

                auto t10_xxx_y = primBuffer.data(t10off + 30 * idx + 1);

                auto t10_xxx_z = primBuffer.data(t10off + 30 * idx + 2);

                auto t10_xxy_x = primBuffer.data(t10off + 30 * idx + 3);

                auto t10_xxy_y = primBuffer.data(t10off + 30 * idx + 4);

                auto t10_xxy_z = primBuffer.data(t10off + 30 * idx + 5);

                auto t10_xxz_x = primBuffer.data(t10off + 30 * idx + 6);

                auto t10_xxz_y = primBuffer.data(t10off + 30 * idx + 7);

                auto t10_xxz_z = primBuffer.data(t10off + 30 * idx + 8);

                auto t10_xyy_x = primBuffer.data(t10off + 30 * idx + 9);

                auto t10_xyy_y = primBuffer.data(t10off + 30 * idx + 10);

                auto t10_xyy_z = primBuffer.data(t10off + 30 * idx + 11);

                auto t10_xyz_x = primBuffer.data(t10off + 30 * idx + 12);

                auto t10_xyz_y = primBuffer.data(t10off + 30 * idx + 13);

                auto t10_xyz_z = primBuffer.data(t10off + 30 * idx + 14);

                auto t10_xzz_x = primBuffer.data(t10off + 30 * idx + 15);

                auto t10_xzz_y = primBuffer.data(t10off + 30 * idx + 16);

                auto t10_xzz_z = primBuffer.data(t10off + 30 * idx + 17);

                auto t10_yyy_x = primBuffer.data(t10off + 30 * idx + 18);

                auto t10_yyy_y = primBuffer.data(t10off + 30 * idx + 19);

                auto t10_yyy_z = primBuffer.data(t10off + 30 * idx + 20);

                auto t10_yyz_x = primBuffer.data(t10off + 30 * idx + 21);

                auto t10_yyz_y = primBuffer.data(t10off + 30 * idx + 22);

                auto t10_yyz_z = primBuffer.data(t10off + 30 * idx + 23);

                auto t10_yzz_x = primBuffer.data(t10off + 30 * idx + 24);

                auto t10_yzz_y = primBuffer.data(t10off + 30 * idx + 25);

                auto t10_yzz_z = primBuffer.data(t10off + 30 * idx + 26);

                auto t10_zzz_x = primBuffer.data(t10off + 30 * idx + 27);

                auto t10_zzz_y = primBuffer.data(t10off + 30 * idx + 28);

                auto t10_zzz_z = primBuffer.data(t10off + 30 * idx + 29);

                // set up pointers to (F|A(0)|P)^(m+1) integrals

                auto t11_xxx_x = primBuffer.data(t11off + 30 * idx);

                auto t11_xxx_y = primBuffer.data(t11off + 30 * idx + 1);

                auto t11_xxx_z = primBuffer.data(t11off + 30 * idx + 2);

                auto t11_xxy_x = primBuffer.data(t11off + 30 * idx + 3);

                auto t11_xxy_y = primBuffer.data(t11off + 30 * idx + 4);

                auto t11_xxy_z = primBuffer.data(t11off + 30 * idx + 5);

                auto t11_xxz_x = primBuffer.data(t11off + 30 * idx + 6);

                auto t11_xxz_y = primBuffer.data(t11off + 30 * idx + 7);

                auto t11_xxz_z = primBuffer.data(t11off + 30 * idx + 8);

                auto t11_xyy_x = primBuffer.data(t11off + 30 * idx + 9);

                auto t11_xyy_y = primBuffer.data(t11off + 30 * idx + 10);

                auto t11_xyy_z = primBuffer.data(t11off + 30 * idx + 11);

                auto t11_xyz_x = primBuffer.data(t11off + 30 * idx + 12);

                auto t11_xyz_y = primBuffer.data(t11off + 30 * idx + 13);

                auto t11_xyz_z = primBuffer.data(t11off + 30 * idx + 14);

                auto t11_xzz_x = primBuffer.data(t11off + 30 * idx + 15);

                auto t11_xzz_y = primBuffer.data(t11off + 30 * idx + 16);

                auto t11_xzz_z = primBuffer.data(t11off + 30 * idx + 17);

                auto t11_yyy_x = primBuffer.data(t11off + 30 * idx + 18);

                auto t11_yyy_y = primBuffer.data(t11off + 30 * idx + 19);

                auto t11_yyy_z = primBuffer.data(t11off + 30 * idx + 20);

                auto t11_yyz_x = primBuffer.data(t11off + 30 * idx + 21);

                auto t11_yyz_y = primBuffer.data(t11off + 30 * idx + 22);

                auto t11_yyz_z = primBuffer.data(t11off + 30 * idx + 23);

                auto t11_yzz_x = primBuffer.data(t11off + 30 * idx + 24);

                auto t11_yzz_y = primBuffer.data(t11off + 30 * idx + 25);

                auto t11_yzz_z = primBuffer.data(t11off + 30 * idx + 26);

                auto t11_zzz_x = primBuffer.data(t11off + 30 * idx + 27);

                auto t11_zzz_y = primBuffer.data(t11off + 30 * idx + 28);

                auto t11_zzz_z = primBuffer.data(t11off + 30 * idx + 29);

                // set up pointers to (G|A(0)|P)^(m) integrals

                auto t_xxxx_x = primBuffer.data(toff + 45 * idx);

                auto t_xxxx_y = primBuffer.data(toff + 45 * idx + 1);

                auto t_xxxx_z = primBuffer.data(toff + 45 * idx + 2);

                auto t_xxxy_x = primBuffer.data(toff + 45 * idx + 3);

                auto t_xxxy_y = primBuffer.data(toff + 45 * idx + 4);

                auto t_xxxy_z = primBuffer.data(toff + 45 * idx + 5);

                auto t_xxxz_x = primBuffer.data(toff + 45 * idx + 6);

                auto t_xxxz_y = primBuffer.data(toff + 45 * idx + 7);

                auto t_xxxz_z = primBuffer.data(toff + 45 * idx + 8);

                auto t_xxyy_x = primBuffer.data(toff + 45 * idx + 9);

                auto t_xxyy_y = primBuffer.data(toff + 45 * idx + 10);

                auto t_xxyy_z = primBuffer.data(toff + 45 * idx + 11);

                auto t_xxyz_x = primBuffer.data(toff + 45 * idx + 12);

                auto t_xxyz_y = primBuffer.data(toff + 45 * idx + 13);

                auto t_xxyz_z = primBuffer.data(toff + 45 * idx + 14);

                auto t_xxzz_x = primBuffer.data(toff + 45 * idx + 15);

                auto t_xxzz_y = primBuffer.data(toff + 45 * idx + 16);

                auto t_xxzz_z = primBuffer.data(toff + 45 * idx + 17);

                auto t_xyyy_x = primBuffer.data(toff + 45 * idx + 18);

                auto t_xyyy_y = primBuffer.data(toff + 45 * idx + 19);

                auto t_xyyy_z = primBuffer.data(toff + 45 * idx + 20);

                auto t_xyyz_x = primBuffer.data(toff + 45 * idx + 21);

                auto t_xyyz_y = primBuffer.data(toff + 45 * idx + 22);

                auto t_xyyz_z = primBuffer.data(toff + 45 * idx + 23);

                auto t_xyzz_x = primBuffer.data(toff + 45 * idx + 24);

                auto t_xyzz_y = primBuffer.data(toff + 45 * idx + 25);

                auto t_xyzz_z = primBuffer.data(toff + 45 * idx + 26);

                auto t_xzzz_x = primBuffer.data(toff + 45 * idx + 27);

                auto t_xzzz_y = primBuffer.data(toff + 45 * idx + 28);

                auto t_xzzz_z = primBuffer.data(toff + 45 * idx + 29);

                auto t_yyyy_x = primBuffer.data(toff + 45 * idx + 30);

                auto t_yyyy_y = primBuffer.data(toff + 45 * idx + 31);

                auto t_yyyy_z = primBuffer.data(toff + 45 * idx + 32);

                auto t_yyyz_x = primBuffer.data(toff + 45 * idx + 33);

                auto t_yyyz_y = primBuffer.data(toff + 45 * idx + 34);

                auto t_yyyz_z = primBuffer.data(toff + 45 * idx + 35);

                auto t_yyzz_x = primBuffer.data(toff + 45 * idx + 36);

                auto t_yyzz_y = primBuffer.data(toff + 45 * idx + 37);

                auto t_yyzz_z = primBuffer.data(toff + 45 * idx + 38);

                auto t_yzzz_x = primBuffer.data(toff + 45 * idx + 39);

                auto t_yzzz_y = primBuffer.data(toff + 45 * idx + 40);

                auto t_yzzz_z = primBuffer.data(toff + 45 * idx + 41);

                auto t_zzzz_x = primBuffer.data(toff + 45 * idx + 42);

                auto t_zzzz_y = primBuffer.data(toff + 45 * idx + 43);

                auto t_zzzz_z = primBuffer.data(toff + 45 * idx + 44);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xxx_0, tk0_xxy_0,\
                                         tk0_xxz_0, tk0_xyy_0, tk0_xyz_0, tk0_xzz_0,\
                                         tk0_yyy_0, tk0_yyz_0, tk0_yzz_0, tk0_zzz_0,\
                                         tk1_xxx_0, tk1_xxy_0, tk1_xxz_0, tk1_xyy_0,\
                                         tk1_xyz_0, tk1_xzz_0, tk1_yyy_0, tk1_yyz_0,\
                                         tk1_yzz_0, tk1_zzz_0, t20_xx_x, t20_xx_y,\
                                         t20_xx_z, t20_xy_x, t20_xy_y, t20_xy_z,\
                                         t20_xz_x, t20_xz_y, t20_xz_z, t20_yy_x,\
                                         t20_yy_y, t20_yy_z, t20_yz_x, t20_yz_y,\
                                         t20_yz_z, t20_zz_x, t20_zz_y, t20_zz_z,\
                                         t21_xx_x, t21_xx_y, t21_xx_z, t21_xy_x,\
                                         t21_xy_y, t21_xy_z, t21_xz_x, t21_xz_y,\
                                         t21_xz_z, t21_yy_x, t21_yy_y, t21_yy_z,\
                                         t21_yz_x, t21_yz_y, t21_yz_z, t21_zz_x,\
                                         t21_zz_y, t21_zz_z, t10_xxx_x, t10_xxx_y,\
                                         t10_xxx_z, t10_xxy_x, t10_xxy_y, t10_xxy_z,\
                                         t10_xxz_x, t10_xxz_y, t10_xxz_z, t10_xyy_x,\
                                         t10_xyy_y, t10_xyy_z, t10_xyz_x, t10_xyz_y,\
                                         t10_xyz_z, t10_xzz_x, t10_xzz_y, t10_xzz_z,\
                                         t10_yyy_x, t10_yyy_y, t10_yyy_z, t10_yyz_x,\
                                         t10_yyz_y, t10_yyz_z, t10_yzz_x, t10_yzz_y,\
                                         t10_yzz_z, t10_zzz_x, t10_zzz_y, t10_zzz_z,\
                                         t11_xxx_x, t11_xxx_y, t11_xxx_z, t11_xxy_x,\
                                         t11_xxy_y, t11_xxy_z, t11_xxz_x, t11_xxz_y,\
                                         t11_xxz_z, t11_xyy_x, t11_xyy_y, t11_xyy_z,\
                                         t11_xyz_x, t11_xyz_y, t11_xyz_z, t11_xzz_x,\
                                         t11_xzz_y, t11_xzz_z, t11_yyy_x, t11_yyy_y,\
                                         t11_yyy_z, t11_yyz_x, t11_yyz_y, t11_yyz_z,\
                                         t11_yzz_x, t11_yzz_y, t11_yzz_z, t11_zzz_x,\
                                         t11_zzz_y, t11_zzz_z, t_xxxx_x, t_xxxx_y,\
                                         t_xxxx_z, t_xxxy_x, t_xxxy_y, t_xxxy_z,\
                                         t_xxxz_x, t_xxxz_y, t_xxxz_z, t_xxyy_x,\
                                         t_xxyy_y, t_xxyy_z, t_xxyz_x, t_xxyz_y,\
                                         t_xxyz_z, t_xxzz_x, t_xxzz_y, t_xxzz_z,\
                                         t_xyyy_x, t_xyyy_y, t_xyyy_z, t_xyyz_x,\
                                         t_xyyz_y, t_xyyz_z, t_xyzz_x, t_xyzz_y,\
                                         t_xyzz_z, t_xzzz_x, t_xzzz_y, t_xzzz_z,\
                                         t_yyyy_x, t_yyyy_y, t_yyyy_z, t_yyyz_x,\
                                         t_yyyz_y, t_yyyz_z, t_yyzz_x, t_yyzz_y,\
                                         t_yyzz_z, t_yzzz_x, t_yzzz_y, t_yzzz_z,\
                                         t_zzzz_x, t_zzzz_y, t_zzzz_z, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxxx_x[k] = fra * t10_xxx_x[k] - frc * t11_xxx_x[k] + f2t * (3.0 * t20_xx_x[k] - 3.0 * t21_xx_x[k] + tk0_xxx_0[k] - tk1_xxx_0[k]);

                    t_xxxx_y[k] = fra * t10_xxx_y[k] - frc * t11_xxx_y[k] + f2t * (3.0 * t20_xx_y[k] - 3.0 * t21_xx_y[k]);

                    t_xxxx_z[k] = fra * t10_xxx_z[k] - frc * t11_xxx_z[k] + f2t * (3.0 * t20_xx_z[k] - 3.0 * t21_xx_z[k]);

                    t_xxxy_x[k] = fra * t10_xxy_x[k] - frc * t11_xxy_x[k] + f2t * (2.0 * t20_xy_x[k] - 2.0 * t21_xy_x[k] + tk0_xxy_0[k] - tk1_xxy_0[k]);

                    t_xxxy_y[k] = fra * t10_xxy_y[k] - frc * t11_xxy_y[k] + f2t * (2.0 * t20_xy_y[k] - 2.0 * t21_xy_y[k]);

                    t_xxxy_z[k] = fra * t10_xxy_z[k] - frc * t11_xxy_z[k] + f2t * (2.0 * t20_xy_z[k] - 2.0 * t21_xy_z[k]);

                    t_xxxz_x[k] = fra * t10_xxz_x[k] - frc * t11_xxz_x[k] + f2t * (2.0 * t20_xz_x[k] - 2.0 * t21_xz_x[k] + tk0_xxz_0[k] - tk1_xxz_0[k]);

                    t_xxxz_y[k] = fra * t10_xxz_y[k] - frc * t11_xxz_y[k] + f2t * (2.0 * t20_xz_y[k] - 2.0 * t21_xz_y[k]);

                    t_xxxz_z[k] = fra * t10_xxz_z[k] - frc * t11_xxz_z[k] + f2t * (2.0 * t20_xz_z[k] - 2.0 * t21_xz_z[k]);

                    t_xxyy_x[k] = fra * t10_xyy_x[k] - frc * t11_xyy_x[k] + f2t * (t20_yy_x[k] - t21_yy_x[k] + tk0_xyy_0[k] - tk1_xyy_0[k]);

                    t_xxyy_y[k] = fra * t10_xyy_y[k] - frc * t11_xyy_y[k] + f2t * (t20_yy_y[k] - t21_yy_y[k]);

                    t_xxyy_z[k] = fra * t10_xyy_z[k] - frc * t11_xyy_z[k] + f2t * (t20_yy_z[k] - t21_yy_z[k]);

                    t_xxyz_x[k] = fra * t10_xyz_x[k] - frc * t11_xyz_x[k] + f2t * (t20_yz_x[k] - t21_yz_x[k] + tk0_xyz_0[k] - tk1_xyz_0[k]);

                    t_xxyz_y[k] = fra * t10_xyz_y[k] - frc * t11_xyz_y[k] + f2t * (t20_yz_y[k] - t21_yz_y[k]);

                    t_xxyz_z[k] = fra * t10_xyz_z[k] - frc * t11_xyz_z[k] + f2t * (t20_yz_z[k] - t21_yz_z[k]);

                    t_xxzz_x[k] = fra * t10_xzz_x[k] - frc * t11_xzz_x[k] + f2t * (t20_zz_x[k] - t21_zz_x[k] + tk0_xzz_0[k] - tk1_xzz_0[k]);

                    t_xxzz_y[k] = fra * t10_xzz_y[k] - frc * t11_xzz_y[k] + f2t * (t20_zz_y[k] - t21_zz_y[k]);

                    t_xxzz_z[k] = fra * t10_xzz_z[k] - frc * t11_xzz_z[k] + f2t * (t20_zz_z[k] - t21_zz_z[k]);

                    t_xyyy_x[k] = fra * t10_yyy_x[k] - frc * t11_yyy_x[k] + f2t * (tk0_yyy_0[k] - tk1_yyy_0[k]);

                    t_xyyy_y[k] = fra * t10_yyy_y[k] - frc * t11_yyy_y[k];

                    t_xyyy_z[k] = fra * t10_yyy_z[k] - frc * t11_yyy_z[k];

                    t_xyyz_x[k] = fra * t10_yyz_x[k] - frc * t11_yyz_x[k] + f2t * (tk0_yyz_0[k] - tk1_yyz_0[k]);

                    t_xyyz_y[k] = fra * t10_yyz_y[k] - frc * t11_yyz_y[k];

                    t_xyyz_z[k] = fra * t10_yyz_z[k] - frc * t11_yyz_z[k];

                    t_xyzz_x[k] = fra * t10_yzz_x[k] - frc * t11_yzz_x[k] + f2t * (tk0_yzz_0[k] - tk1_yzz_0[k]);

                    t_xyzz_y[k] = fra * t10_yzz_y[k] - frc * t11_yzz_y[k];

                    t_xyzz_z[k] = fra * t10_yzz_z[k] - frc * t11_yzz_z[k];

                    t_xzzz_x[k] = fra * t10_zzz_x[k] - frc * t11_zzz_x[k] + f2t * (tk0_zzz_0[k] - tk1_zzz_0[k]);

                    t_xzzz_y[k] = fra * t10_zzz_y[k] - frc * t11_zzz_y[k];

                    t_xzzz_z[k] = fra * t10_zzz_z[k] - frc * t11_zzz_z[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyyy_x[k] = fra * t10_yyy_x[k] - frc * t11_yyy_x[k] + f2t * (3.0 * t20_yy_x[k] - 3.0 * t21_yy_x[k]);

                    t_yyyy_y[k] = fra * t10_yyy_y[k] - frc * t11_yyy_y[k] + f2t * (3.0 * t20_yy_y[k] - 3.0 * t21_yy_y[k] + tk0_yyy_0[k] - tk1_yyy_0[k]);

                    t_yyyy_z[k] = fra * t10_yyy_z[k] - frc * t11_yyy_z[k] + f2t * (3.0 * t20_yy_z[k] - 3.0 * t21_yy_z[k]);

                    t_yyyz_x[k] = fra * t10_yyz_x[k] - frc * t11_yyz_x[k] + f2t * (2.0 * t20_yz_x[k] - 2.0 * t21_yz_x[k]);

                    t_yyyz_y[k] = fra * t10_yyz_y[k] - frc * t11_yyz_y[k] + f2t * (2.0 * t20_yz_y[k] - 2.0 * t21_yz_y[k] + tk0_yyz_0[k] - tk1_yyz_0[k]);

                    t_yyyz_z[k] = fra * t10_yyz_z[k] - frc * t11_yyz_z[k] + f2t * (2.0 * t20_yz_z[k] - 2.0 * t21_yz_z[k]);

                    t_yyzz_x[k] = fra * t10_yzz_x[k] - frc * t11_yzz_x[k] + f2t * (t20_zz_x[k] - t21_zz_x[k]);

                    t_yyzz_y[k] = fra * t10_yzz_y[k] - frc * t11_yzz_y[k] + f2t * (t20_zz_y[k] - t21_zz_y[k] + tk0_yzz_0[k] - tk1_yzz_0[k]);

                    t_yyzz_z[k] = fra * t10_yzz_z[k] - frc * t11_yzz_z[k] + f2t * (t20_zz_z[k] - t21_zz_z[k]);

                    t_yzzz_x[k] = fra * t10_zzz_x[k] - frc * t11_zzz_x[k];

                    t_yzzz_y[k] = fra * t10_zzz_y[k] - frc * t11_zzz_y[k] + f2t * (tk0_zzz_0[k] - tk1_zzz_0[k]);

                    t_yzzz_z[k] = fra * t10_zzz_z[k] - frc * t11_zzz_z[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzzz_x[k] = fra * t10_zzz_x[k] - frc * t11_zzz_x[k] + f2t * (3.0 * t20_zz_x[k] - 3.0 * t21_zz_x[k]);

                    t_zzzz_y[k] = fra * t10_zzz_y[k] - frc * t11_zzz_y[k] + f2t * (3.0 * t20_zz_y[k] - 3.0 * t21_zz_y[k]);

                    t_zzzz_z[k] = fra * t10_zzz_z[k] - frc * t11_zzz_z[k] + f2t * (3.0 * t20_zz_z[k] - 3.0 * t21_zz_z[k] + tk0_zzz_0[k] - tk1_zzz_0[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForDG(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {2, 4, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 2, 4);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 4, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 4, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 4, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 4, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 4, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 3, i + 1});

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

                // set up pointers to (P|A(0)|F)^(m) integrals

                auto tk0_x_xxx = primBuffer.data(tk0off + 30 * idx);

                auto tk0_x_xxy = primBuffer.data(tk0off + 30 * idx + 1);

                auto tk0_x_xxz = primBuffer.data(tk0off + 30 * idx + 2);

                auto tk0_x_xyy = primBuffer.data(tk0off + 30 * idx + 3);

                auto tk0_x_xyz = primBuffer.data(tk0off + 30 * idx + 4);

                auto tk0_x_xzz = primBuffer.data(tk0off + 30 * idx + 5);

                auto tk0_x_yyy = primBuffer.data(tk0off + 30 * idx + 6);

                auto tk0_x_yyz = primBuffer.data(tk0off + 30 * idx + 7);

                auto tk0_x_yzz = primBuffer.data(tk0off + 30 * idx + 8);

                auto tk0_x_zzz = primBuffer.data(tk0off + 30 * idx + 9);

                auto tk0_y_xxx = primBuffer.data(tk0off + 30 * idx + 10);

                auto tk0_y_xxy = primBuffer.data(tk0off + 30 * idx + 11);

                auto tk0_y_xxz = primBuffer.data(tk0off + 30 * idx + 12);

                auto tk0_y_xyy = primBuffer.data(tk0off + 30 * idx + 13);

                auto tk0_y_xyz = primBuffer.data(tk0off + 30 * idx + 14);

                auto tk0_y_xzz = primBuffer.data(tk0off + 30 * idx + 15);

                auto tk0_y_yyy = primBuffer.data(tk0off + 30 * idx + 16);

                auto tk0_y_yyz = primBuffer.data(tk0off + 30 * idx + 17);

                auto tk0_y_yzz = primBuffer.data(tk0off + 30 * idx + 18);

                auto tk0_y_zzz = primBuffer.data(tk0off + 30 * idx + 19);

                auto tk0_z_xxx = primBuffer.data(tk0off + 30 * idx + 20);

                auto tk0_z_xxy = primBuffer.data(tk0off + 30 * idx + 21);

                auto tk0_z_xxz = primBuffer.data(tk0off + 30 * idx + 22);

                auto tk0_z_xyy = primBuffer.data(tk0off + 30 * idx + 23);

                auto tk0_z_xyz = primBuffer.data(tk0off + 30 * idx + 24);

                auto tk0_z_xzz = primBuffer.data(tk0off + 30 * idx + 25);

                auto tk0_z_yyy = primBuffer.data(tk0off + 30 * idx + 26);

                auto tk0_z_yyz = primBuffer.data(tk0off + 30 * idx + 27);

                auto tk0_z_yzz = primBuffer.data(tk0off + 30 * idx + 28);

                auto tk0_z_zzz = primBuffer.data(tk0off + 30 * idx + 29);

                // set up pointers to (P|A(0)|F)^(m+1) integrals

                auto tk1_x_xxx = primBuffer.data(tk1off + 30 * idx);

                auto tk1_x_xxy = primBuffer.data(tk1off + 30 * idx + 1);

                auto tk1_x_xxz = primBuffer.data(tk1off + 30 * idx + 2);

                auto tk1_x_xyy = primBuffer.data(tk1off + 30 * idx + 3);

                auto tk1_x_xyz = primBuffer.data(tk1off + 30 * idx + 4);

                auto tk1_x_xzz = primBuffer.data(tk1off + 30 * idx + 5);

                auto tk1_x_yyy = primBuffer.data(tk1off + 30 * idx + 6);

                auto tk1_x_yyz = primBuffer.data(tk1off + 30 * idx + 7);

                auto tk1_x_yzz = primBuffer.data(tk1off + 30 * idx + 8);

                auto tk1_x_zzz = primBuffer.data(tk1off + 30 * idx + 9);

                auto tk1_y_xxx = primBuffer.data(tk1off + 30 * idx + 10);

                auto tk1_y_xxy = primBuffer.data(tk1off + 30 * idx + 11);

                auto tk1_y_xxz = primBuffer.data(tk1off + 30 * idx + 12);

                auto tk1_y_xyy = primBuffer.data(tk1off + 30 * idx + 13);

                auto tk1_y_xyz = primBuffer.data(tk1off + 30 * idx + 14);

                auto tk1_y_xzz = primBuffer.data(tk1off + 30 * idx + 15);

                auto tk1_y_yyy = primBuffer.data(tk1off + 30 * idx + 16);

                auto tk1_y_yyz = primBuffer.data(tk1off + 30 * idx + 17);

                auto tk1_y_yzz = primBuffer.data(tk1off + 30 * idx + 18);

                auto tk1_y_zzz = primBuffer.data(tk1off + 30 * idx + 19);

                auto tk1_z_xxx = primBuffer.data(tk1off + 30 * idx + 20);

                auto tk1_z_xxy = primBuffer.data(tk1off + 30 * idx + 21);

                auto tk1_z_xxz = primBuffer.data(tk1off + 30 * idx + 22);

                auto tk1_z_xyy = primBuffer.data(tk1off + 30 * idx + 23);

                auto tk1_z_xyz = primBuffer.data(tk1off + 30 * idx + 24);

                auto tk1_z_xzz = primBuffer.data(tk1off + 30 * idx + 25);

                auto tk1_z_yyy = primBuffer.data(tk1off + 30 * idx + 26);

                auto tk1_z_yyz = primBuffer.data(tk1off + 30 * idx + 27);

                auto tk1_z_yzz = primBuffer.data(tk1off + 30 * idx + 28);

                auto tk1_z_zzz = primBuffer.data(tk1off + 30 * idx + 29);

                // set up pointers to (S|A(0)|G)^(m) integrals

                auto t20_0_xxxx = primBuffer.data(t20off + 15 * idx);

                auto t20_0_xxxy = primBuffer.data(t20off + 15 * idx + 1);

                auto t20_0_xxxz = primBuffer.data(t20off + 15 * idx + 2);

                auto t20_0_xxyy = primBuffer.data(t20off + 15 * idx + 3);

                auto t20_0_xxyz = primBuffer.data(t20off + 15 * idx + 4);

                auto t20_0_xxzz = primBuffer.data(t20off + 15 * idx + 5);

                auto t20_0_xyyy = primBuffer.data(t20off + 15 * idx + 6);

                auto t20_0_xyyz = primBuffer.data(t20off + 15 * idx + 7);

                auto t20_0_xyzz = primBuffer.data(t20off + 15 * idx + 8);

                auto t20_0_xzzz = primBuffer.data(t20off + 15 * idx + 9);

                auto t20_0_yyyy = primBuffer.data(t20off + 15 * idx + 10);

                auto t20_0_yyyz = primBuffer.data(t20off + 15 * idx + 11);

                auto t20_0_yyzz = primBuffer.data(t20off + 15 * idx + 12);

                auto t20_0_yzzz = primBuffer.data(t20off + 15 * idx + 13);

                auto t20_0_zzzz = primBuffer.data(t20off + 15 * idx + 14);

                // set up pointers to (S|A(0)|G)^(m+1) integrals

                auto t21_0_xxxx = primBuffer.data(t21off + 15 * idx);

                auto t21_0_xxxy = primBuffer.data(t21off + 15 * idx + 1);

                auto t21_0_xxxz = primBuffer.data(t21off + 15 * idx + 2);

                auto t21_0_xxyy = primBuffer.data(t21off + 15 * idx + 3);

                auto t21_0_xxyz = primBuffer.data(t21off + 15 * idx + 4);

                auto t21_0_xxzz = primBuffer.data(t21off + 15 * idx + 5);

                auto t21_0_xyyy = primBuffer.data(t21off + 15 * idx + 6);

                auto t21_0_xyyz = primBuffer.data(t21off + 15 * idx + 7);

                auto t21_0_xyzz = primBuffer.data(t21off + 15 * idx + 8);

                auto t21_0_xzzz = primBuffer.data(t21off + 15 * idx + 9);

                auto t21_0_yyyy = primBuffer.data(t21off + 15 * idx + 10);

                auto t21_0_yyyz = primBuffer.data(t21off + 15 * idx + 11);

                auto t21_0_yyzz = primBuffer.data(t21off + 15 * idx + 12);

                auto t21_0_yzzz = primBuffer.data(t21off + 15 * idx + 13);

                auto t21_0_zzzz = primBuffer.data(t21off + 15 * idx + 14);

                // set up pointers to (P|A(0)|G)^(m) integrals

                auto t10_x_xxxx = primBuffer.data(t10off + 45 * idx);

                auto t10_x_xxxy = primBuffer.data(t10off + 45 * idx + 1);

                auto t10_x_xxxz = primBuffer.data(t10off + 45 * idx + 2);

                auto t10_x_xxyy = primBuffer.data(t10off + 45 * idx + 3);

                auto t10_x_xxyz = primBuffer.data(t10off + 45 * idx + 4);

                auto t10_x_xxzz = primBuffer.data(t10off + 45 * idx + 5);

                auto t10_x_xyyy = primBuffer.data(t10off + 45 * idx + 6);

                auto t10_x_xyyz = primBuffer.data(t10off + 45 * idx + 7);

                auto t10_x_xyzz = primBuffer.data(t10off + 45 * idx + 8);

                auto t10_x_xzzz = primBuffer.data(t10off + 45 * idx + 9);

                auto t10_x_yyyy = primBuffer.data(t10off + 45 * idx + 10);

                auto t10_x_yyyz = primBuffer.data(t10off + 45 * idx + 11);

                auto t10_x_yyzz = primBuffer.data(t10off + 45 * idx + 12);

                auto t10_x_yzzz = primBuffer.data(t10off + 45 * idx + 13);

                auto t10_x_zzzz = primBuffer.data(t10off + 45 * idx + 14);

                auto t10_y_xxxx = primBuffer.data(t10off + 45 * idx + 15);

                auto t10_y_xxxy = primBuffer.data(t10off + 45 * idx + 16);

                auto t10_y_xxxz = primBuffer.data(t10off + 45 * idx + 17);

                auto t10_y_xxyy = primBuffer.data(t10off + 45 * idx + 18);

                auto t10_y_xxyz = primBuffer.data(t10off + 45 * idx + 19);

                auto t10_y_xxzz = primBuffer.data(t10off + 45 * idx + 20);

                auto t10_y_xyyy = primBuffer.data(t10off + 45 * idx + 21);

                auto t10_y_xyyz = primBuffer.data(t10off + 45 * idx + 22);

                auto t10_y_xyzz = primBuffer.data(t10off + 45 * idx + 23);

                auto t10_y_xzzz = primBuffer.data(t10off + 45 * idx + 24);

                auto t10_y_yyyy = primBuffer.data(t10off + 45 * idx + 25);

                auto t10_y_yyyz = primBuffer.data(t10off + 45 * idx + 26);

                auto t10_y_yyzz = primBuffer.data(t10off + 45 * idx + 27);

                auto t10_y_yzzz = primBuffer.data(t10off + 45 * idx + 28);

                auto t10_y_zzzz = primBuffer.data(t10off + 45 * idx + 29);

                auto t10_z_xxxx = primBuffer.data(t10off + 45 * idx + 30);

                auto t10_z_xxxy = primBuffer.data(t10off + 45 * idx + 31);

                auto t10_z_xxxz = primBuffer.data(t10off + 45 * idx + 32);

                auto t10_z_xxyy = primBuffer.data(t10off + 45 * idx + 33);

                auto t10_z_xxyz = primBuffer.data(t10off + 45 * idx + 34);

                auto t10_z_xxzz = primBuffer.data(t10off + 45 * idx + 35);

                auto t10_z_xyyy = primBuffer.data(t10off + 45 * idx + 36);

                auto t10_z_xyyz = primBuffer.data(t10off + 45 * idx + 37);

                auto t10_z_xyzz = primBuffer.data(t10off + 45 * idx + 38);

                auto t10_z_xzzz = primBuffer.data(t10off + 45 * idx + 39);

                auto t10_z_yyyy = primBuffer.data(t10off + 45 * idx + 40);

                auto t10_z_yyyz = primBuffer.data(t10off + 45 * idx + 41);

                auto t10_z_yyzz = primBuffer.data(t10off + 45 * idx + 42);

                auto t10_z_yzzz = primBuffer.data(t10off + 45 * idx + 43);

                auto t10_z_zzzz = primBuffer.data(t10off + 45 * idx + 44);

                // set up pointers to (P|A(0)|G)^(m+1) integrals

                auto t11_x_xxxx = primBuffer.data(t11off + 45 * idx);

                auto t11_x_xxxy = primBuffer.data(t11off + 45 * idx + 1);

                auto t11_x_xxxz = primBuffer.data(t11off + 45 * idx + 2);

                auto t11_x_xxyy = primBuffer.data(t11off + 45 * idx + 3);

                auto t11_x_xxyz = primBuffer.data(t11off + 45 * idx + 4);

                auto t11_x_xxzz = primBuffer.data(t11off + 45 * idx + 5);

                auto t11_x_xyyy = primBuffer.data(t11off + 45 * idx + 6);

                auto t11_x_xyyz = primBuffer.data(t11off + 45 * idx + 7);

                auto t11_x_xyzz = primBuffer.data(t11off + 45 * idx + 8);

                auto t11_x_xzzz = primBuffer.data(t11off + 45 * idx + 9);

                auto t11_x_yyyy = primBuffer.data(t11off + 45 * idx + 10);

                auto t11_x_yyyz = primBuffer.data(t11off + 45 * idx + 11);

                auto t11_x_yyzz = primBuffer.data(t11off + 45 * idx + 12);

                auto t11_x_yzzz = primBuffer.data(t11off + 45 * idx + 13);

                auto t11_x_zzzz = primBuffer.data(t11off + 45 * idx + 14);

                auto t11_y_xxxx = primBuffer.data(t11off + 45 * idx + 15);

                auto t11_y_xxxy = primBuffer.data(t11off + 45 * idx + 16);

                auto t11_y_xxxz = primBuffer.data(t11off + 45 * idx + 17);

                auto t11_y_xxyy = primBuffer.data(t11off + 45 * idx + 18);

                auto t11_y_xxyz = primBuffer.data(t11off + 45 * idx + 19);

                auto t11_y_xxzz = primBuffer.data(t11off + 45 * idx + 20);

                auto t11_y_xyyy = primBuffer.data(t11off + 45 * idx + 21);

                auto t11_y_xyyz = primBuffer.data(t11off + 45 * idx + 22);

                auto t11_y_xyzz = primBuffer.data(t11off + 45 * idx + 23);

                auto t11_y_xzzz = primBuffer.data(t11off + 45 * idx + 24);

                auto t11_y_yyyy = primBuffer.data(t11off + 45 * idx + 25);

                auto t11_y_yyyz = primBuffer.data(t11off + 45 * idx + 26);

                auto t11_y_yyzz = primBuffer.data(t11off + 45 * idx + 27);

                auto t11_y_yzzz = primBuffer.data(t11off + 45 * idx + 28);

                auto t11_y_zzzz = primBuffer.data(t11off + 45 * idx + 29);

                auto t11_z_xxxx = primBuffer.data(t11off + 45 * idx + 30);

                auto t11_z_xxxy = primBuffer.data(t11off + 45 * idx + 31);

                auto t11_z_xxxz = primBuffer.data(t11off + 45 * idx + 32);

                auto t11_z_xxyy = primBuffer.data(t11off + 45 * idx + 33);

                auto t11_z_xxyz = primBuffer.data(t11off + 45 * idx + 34);

                auto t11_z_xxzz = primBuffer.data(t11off + 45 * idx + 35);

                auto t11_z_xyyy = primBuffer.data(t11off + 45 * idx + 36);

                auto t11_z_xyyz = primBuffer.data(t11off + 45 * idx + 37);

                auto t11_z_xyzz = primBuffer.data(t11off + 45 * idx + 38);

                auto t11_z_xzzz = primBuffer.data(t11off + 45 * idx + 39);

                auto t11_z_yyyy = primBuffer.data(t11off + 45 * idx + 40);

                auto t11_z_yyyz = primBuffer.data(t11off + 45 * idx + 41);

                auto t11_z_yyzz = primBuffer.data(t11off + 45 * idx + 42);

                auto t11_z_yzzz = primBuffer.data(t11off + 45 * idx + 43);

                auto t11_z_zzzz = primBuffer.data(t11off + 45 * idx + 44);

                // set up pointers to (D|A(0)|G)^(m) integrals

                auto t_xx_xxxx = primBuffer.data(toff + 90 * idx);

                auto t_xx_xxxy = primBuffer.data(toff + 90 * idx + 1);

                auto t_xx_xxxz = primBuffer.data(toff + 90 * idx + 2);

                auto t_xx_xxyy = primBuffer.data(toff + 90 * idx + 3);

                auto t_xx_xxyz = primBuffer.data(toff + 90 * idx + 4);

                auto t_xx_xxzz = primBuffer.data(toff + 90 * idx + 5);

                auto t_xx_xyyy = primBuffer.data(toff + 90 * idx + 6);

                auto t_xx_xyyz = primBuffer.data(toff + 90 * idx + 7);

                auto t_xx_xyzz = primBuffer.data(toff + 90 * idx + 8);

                auto t_xx_xzzz = primBuffer.data(toff + 90 * idx + 9);

                auto t_xx_yyyy = primBuffer.data(toff + 90 * idx + 10);

                auto t_xx_yyyz = primBuffer.data(toff + 90 * idx + 11);

                auto t_xx_yyzz = primBuffer.data(toff + 90 * idx + 12);

                auto t_xx_yzzz = primBuffer.data(toff + 90 * idx + 13);

                auto t_xx_zzzz = primBuffer.data(toff + 90 * idx + 14);

                auto t_xy_xxxx = primBuffer.data(toff + 90 * idx + 15);

                auto t_xy_xxxy = primBuffer.data(toff + 90 * idx + 16);

                auto t_xy_xxxz = primBuffer.data(toff + 90 * idx + 17);

                auto t_xy_xxyy = primBuffer.data(toff + 90 * idx + 18);

                auto t_xy_xxyz = primBuffer.data(toff + 90 * idx + 19);

                auto t_xy_xxzz = primBuffer.data(toff + 90 * idx + 20);

                auto t_xy_xyyy = primBuffer.data(toff + 90 * idx + 21);

                auto t_xy_xyyz = primBuffer.data(toff + 90 * idx + 22);

                auto t_xy_xyzz = primBuffer.data(toff + 90 * idx + 23);

                auto t_xy_xzzz = primBuffer.data(toff + 90 * idx + 24);

                auto t_xy_yyyy = primBuffer.data(toff + 90 * idx + 25);

                auto t_xy_yyyz = primBuffer.data(toff + 90 * idx + 26);

                auto t_xy_yyzz = primBuffer.data(toff + 90 * idx + 27);

                auto t_xy_yzzz = primBuffer.data(toff + 90 * idx + 28);

                auto t_xy_zzzz = primBuffer.data(toff + 90 * idx + 29);

                auto t_xz_xxxx = primBuffer.data(toff + 90 * idx + 30);

                auto t_xz_xxxy = primBuffer.data(toff + 90 * idx + 31);

                auto t_xz_xxxz = primBuffer.data(toff + 90 * idx + 32);

                auto t_xz_xxyy = primBuffer.data(toff + 90 * idx + 33);

                auto t_xz_xxyz = primBuffer.data(toff + 90 * idx + 34);

                auto t_xz_xxzz = primBuffer.data(toff + 90 * idx + 35);

                auto t_xz_xyyy = primBuffer.data(toff + 90 * idx + 36);

                auto t_xz_xyyz = primBuffer.data(toff + 90 * idx + 37);

                auto t_xz_xyzz = primBuffer.data(toff + 90 * idx + 38);

                auto t_xz_xzzz = primBuffer.data(toff + 90 * idx + 39);

                auto t_xz_yyyy = primBuffer.data(toff + 90 * idx + 40);

                auto t_xz_yyyz = primBuffer.data(toff + 90 * idx + 41);

                auto t_xz_yyzz = primBuffer.data(toff + 90 * idx + 42);

                auto t_xz_yzzz = primBuffer.data(toff + 90 * idx + 43);

                auto t_xz_zzzz = primBuffer.data(toff + 90 * idx + 44);

                auto t_yy_xxxx = primBuffer.data(toff + 90 * idx + 45);

                auto t_yy_xxxy = primBuffer.data(toff + 90 * idx + 46);

                auto t_yy_xxxz = primBuffer.data(toff + 90 * idx + 47);

                auto t_yy_xxyy = primBuffer.data(toff + 90 * idx + 48);

                auto t_yy_xxyz = primBuffer.data(toff + 90 * idx + 49);

                auto t_yy_xxzz = primBuffer.data(toff + 90 * idx + 50);

                auto t_yy_xyyy = primBuffer.data(toff + 90 * idx + 51);

                auto t_yy_xyyz = primBuffer.data(toff + 90 * idx + 52);

                auto t_yy_xyzz = primBuffer.data(toff + 90 * idx + 53);

                auto t_yy_xzzz = primBuffer.data(toff + 90 * idx + 54);

                auto t_yy_yyyy = primBuffer.data(toff + 90 * idx + 55);

                auto t_yy_yyyz = primBuffer.data(toff + 90 * idx + 56);

                auto t_yy_yyzz = primBuffer.data(toff + 90 * idx + 57);

                auto t_yy_yzzz = primBuffer.data(toff + 90 * idx + 58);

                auto t_yy_zzzz = primBuffer.data(toff + 90 * idx + 59);

                auto t_yz_xxxx = primBuffer.data(toff + 90 * idx + 60);

                auto t_yz_xxxy = primBuffer.data(toff + 90 * idx + 61);

                auto t_yz_xxxz = primBuffer.data(toff + 90 * idx + 62);

                auto t_yz_xxyy = primBuffer.data(toff + 90 * idx + 63);

                auto t_yz_xxyz = primBuffer.data(toff + 90 * idx + 64);

                auto t_yz_xxzz = primBuffer.data(toff + 90 * idx + 65);

                auto t_yz_xyyy = primBuffer.data(toff + 90 * idx + 66);

                auto t_yz_xyyz = primBuffer.data(toff + 90 * idx + 67);

                auto t_yz_xyzz = primBuffer.data(toff + 90 * idx + 68);

                auto t_yz_xzzz = primBuffer.data(toff + 90 * idx + 69);

                auto t_yz_yyyy = primBuffer.data(toff + 90 * idx + 70);

                auto t_yz_yyyz = primBuffer.data(toff + 90 * idx + 71);

                auto t_yz_yyzz = primBuffer.data(toff + 90 * idx + 72);

                auto t_yz_yzzz = primBuffer.data(toff + 90 * idx + 73);

                auto t_yz_zzzz = primBuffer.data(toff + 90 * idx + 74);

                auto t_zz_xxxx = primBuffer.data(toff + 90 * idx + 75);

                auto t_zz_xxxy = primBuffer.data(toff + 90 * idx + 76);

                auto t_zz_xxxz = primBuffer.data(toff + 90 * idx + 77);

                auto t_zz_xxyy = primBuffer.data(toff + 90 * idx + 78);

                auto t_zz_xxyz = primBuffer.data(toff + 90 * idx + 79);

                auto t_zz_xxzz = primBuffer.data(toff + 90 * idx + 80);

                auto t_zz_xyyy = primBuffer.data(toff + 90 * idx + 81);

                auto t_zz_xyyz = primBuffer.data(toff + 90 * idx + 82);

                auto t_zz_xyzz = primBuffer.data(toff + 90 * idx + 83);

                auto t_zz_xzzz = primBuffer.data(toff + 90 * idx + 84);

                auto t_zz_yyyy = primBuffer.data(toff + 90 * idx + 85);

                auto t_zz_yyyz = primBuffer.data(toff + 90 * idx + 86);

                auto t_zz_yyzz = primBuffer.data(toff + 90 * idx + 87);

                auto t_zz_yzzz = primBuffer.data(toff + 90 * idx + 88);

                auto t_zz_zzzz = primBuffer.data(toff + 90 * idx + 89);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_x_xxx, tk0_x_xxy,\
                                         tk0_x_xxz, tk0_x_xyy, tk0_x_xyz, tk0_x_xzz,\
                                         tk0_x_yyy, tk0_x_yyz, tk0_x_yzz, tk0_x_zzz,\
                                         tk0_y_xxx, tk0_y_xxy, tk0_y_xxz, tk0_y_xyy,\
                                         tk0_y_xyz, tk0_y_xzz, tk0_y_yyy, tk0_y_yyz,\
                                         tk0_y_yzz, tk0_y_zzz, tk0_z_xxx, tk0_z_xxy,\
                                         tk0_z_xxz, tk0_z_xyy, tk0_z_xyz, tk0_z_xzz,\
                                         tk0_z_yyy, tk0_z_yyz, tk0_z_yzz, tk0_z_zzz,\
                                         tk1_x_xxx, tk1_x_xxy, tk1_x_xxz, tk1_x_xyy,\
                                         tk1_x_xyz, tk1_x_xzz, tk1_x_yyy, tk1_x_yyz,\
                                         tk1_x_yzz, tk1_x_zzz, tk1_y_xxx, tk1_y_xxy,\
                                         tk1_y_xxz, tk1_y_xyy, tk1_y_xyz, tk1_y_xzz,\
                                         tk1_y_yyy, tk1_y_yyz, tk1_y_yzz, tk1_y_zzz,\
                                         tk1_z_xxx, tk1_z_xxy, tk1_z_xxz, tk1_z_xyy,\
                                         tk1_z_xyz, tk1_z_xzz, tk1_z_yyy, tk1_z_yyz,\
                                         tk1_z_yzz, tk1_z_zzz, t20_0_xxxx, t20_0_xxxy,\
                                         t20_0_xxxz, t20_0_xxyy, t20_0_xxyz, t20_0_xxzz,\
                                         t20_0_xyyy, t20_0_xyyz, t20_0_xyzz, t20_0_xzzz,\
                                         t20_0_yyyy, t20_0_yyyz, t20_0_yyzz, t20_0_yzzz,\
                                         t20_0_zzzz, t21_0_xxxx, t21_0_xxxy, t21_0_xxxz,\
                                         t21_0_xxyy, t21_0_xxyz, t21_0_xxzz, t21_0_xyyy,\
                                         t21_0_xyyz, t21_0_xyzz, t21_0_xzzz, t21_0_yyyy,\
                                         t21_0_yyyz, t21_0_yyzz, t21_0_yzzz, t21_0_zzzz,\
                                         t10_x_xxxx, t10_x_xxxy, t10_x_xxxz, t10_x_xxyy,\
                                         t10_x_xxyz, t10_x_xxzz, t10_x_xyyy, t10_x_xyyz,\
                                         t10_x_xyzz, t10_x_xzzz, t10_x_yyyy, t10_x_yyyz,\
                                         t10_x_yyzz, t10_x_yzzz, t10_x_zzzz, t10_y_xxxx,\
                                         t10_y_xxxy, t10_y_xxxz, t10_y_xxyy, t10_y_xxyz,\
                                         t10_y_xxzz, t10_y_xyyy, t10_y_xyyz, t10_y_xyzz,\
                                         t10_y_xzzz, t10_y_yyyy, t10_y_yyyz, t10_y_yyzz,\
                                         t10_y_yzzz, t10_y_zzzz, t10_z_xxxx, t10_z_xxxy,\
                                         t10_z_xxxz, t10_z_xxyy, t10_z_xxyz, t10_z_xxzz,\
                                         t10_z_xyyy, t10_z_xyyz, t10_z_xyzz, t10_z_xzzz,\
                                         t10_z_yyyy, t10_z_yyyz, t10_z_yyzz, t10_z_yzzz,\
                                         t10_z_zzzz, t11_x_xxxx, t11_x_xxxy, t11_x_xxxz,\
                                         t11_x_xxyy, t11_x_xxyz, t11_x_xxzz, t11_x_xyyy,\
                                         t11_x_xyyz, t11_x_xyzz, t11_x_xzzz, t11_x_yyyy,\
                                         t11_x_yyyz, t11_x_yyzz, t11_x_yzzz, t11_x_zzzz,\
                                         t11_y_xxxx, t11_y_xxxy, t11_y_xxxz, t11_y_xxyy,\
                                         t11_y_xxyz, t11_y_xxzz, t11_y_xyyy, t11_y_xyyz,\
                                         t11_y_xyzz, t11_y_xzzz, t11_y_yyyy, t11_y_yyyz,\
                                         t11_y_yyzz, t11_y_yzzz, t11_y_zzzz, t11_z_xxxx,\
                                         t11_z_xxxy, t11_z_xxxz, t11_z_xxyy, t11_z_xxyz,\
                                         t11_z_xxzz, t11_z_xyyy, t11_z_xyyz, t11_z_xyzz,\
                                         t11_z_xzzz, t11_z_yyyy, t11_z_yyyz, t11_z_yyzz,\
                                         t11_z_yzzz, t11_z_zzzz, t_xx_xxxx, t_xx_xxxy,\
                                         t_xx_xxxz, t_xx_xxyy, t_xx_xxyz, t_xx_xxzz,\
                                         t_xx_xyyy, t_xx_xyyz, t_xx_xyzz, t_xx_xzzz,\
                                         t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz,\
                                         t_xx_zzzz, t_xy_xxxx, t_xy_xxxy, t_xy_xxxz,\
                                         t_xy_xxyy, t_xy_xxyz, t_xy_xxzz, t_xy_xyyy,\
                                         t_xy_xyyz, t_xy_xyzz, t_xy_xzzz, t_xy_yyyy,\
                                         t_xy_yyyz, t_xy_yyzz, t_xy_yzzz, t_xy_zzzz,\
                                         t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy,\
                                         t_xz_xxyz, t_xz_xxzz, t_xz_xyyy, t_xz_xyyz,\
                                         t_xz_xyzz, t_xz_xzzz, t_xz_yyyy, t_xz_yyyz,\
                                         t_xz_yyzz, t_xz_yzzz, t_xz_zzzz, t_yy_xxxx,\
                                         t_yy_xxxy, t_yy_xxxz, t_yy_xxyy, t_yy_xxyz,\
                                         t_yy_xxzz, t_yy_xyyy, t_yy_xyyz, t_yy_xyzz,\
                                         t_yy_xzzz, t_yy_yyyy, t_yy_yyyz, t_yy_yyzz,\
                                         t_yy_yzzz, t_yy_zzzz, t_yz_xxxx, t_yz_xxxy,\
                                         t_yz_xxxz, t_yz_xxyy, t_yz_xxyz, t_yz_xxzz,\
                                         t_yz_xyyy, t_yz_xyyz, t_yz_xyzz, t_yz_xzzz,\
                                         t_yz_yyyy, t_yz_yyyz, t_yz_yyzz, t_yz_yzzz,\
                                         t_yz_zzzz, t_zz_xxxx, t_zz_xxxy, t_zz_xxxz,\
                                         t_zz_xxyy, t_zz_xxyz, t_zz_xxzz, t_zz_xyyy,\
                                         t_zz_xyyz, t_zz_xyzz, t_zz_xzzz, t_zz_yyyy,\
                                         t_zz_yyyz, t_zz_yyzz, t_zz_yzzz, t_zz_zzzz,\
                                         pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xx_xxxx[k] = fra * t10_x_xxxx[k] - frc * t11_x_xxxx[k] + f2t * (t20_0_xxxx[k] - t21_0_xxxx[k] + 4.0 * tk0_x_xxx[k] - 4.0 * tk1_x_xxx[k]);

                    t_xx_xxxy[k] = fra * t10_x_xxxy[k] - frc * t11_x_xxxy[k] + f2t * (t20_0_xxxy[k] - t21_0_xxxy[k] + 3.0 * tk0_x_xxy[k] - 3.0 * tk1_x_xxy[k]);

                    t_xx_xxxz[k] = fra * t10_x_xxxz[k] - frc * t11_x_xxxz[k] + f2t * (t20_0_xxxz[k] - t21_0_xxxz[k] + 3.0 * tk0_x_xxz[k] - 3.0 * tk1_x_xxz[k]);

                    t_xx_xxyy[k] = fra * t10_x_xxyy[k] - frc * t11_x_xxyy[k] + f2t * (t20_0_xxyy[k] - t21_0_xxyy[k] + 2.0 * tk0_x_xyy[k] - 2.0 * tk1_x_xyy[k]);

                    t_xx_xxyz[k] = fra * t10_x_xxyz[k] - frc * t11_x_xxyz[k] + f2t * (t20_0_xxyz[k] - t21_0_xxyz[k] + 2.0 * tk0_x_xyz[k] - 2.0 * tk1_x_xyz[k]);

                    t_xx_xxzz[k] = fra * t10_x_xxzz[k] - frc * t11_x_xxzz[k] + f2t * (t20_0_xxzz[k] - t21_0_xxzz[k] + 2.0 * tk0_x_xzz[k] - 2.0 * tk1_x_xzz[k]);

                    t_xx_xyyy[k] = fra * t10_x_xyyy[k] - frc * t11_x_xyyy[k] + f2t * (t20_0_xyyy[k] - t21_0_xyyy[k] + tk0_x_yyy[k] - tk1_x_yyy[k]);

                    t_xx_xyyz[k] = fra * t10_x_xyyz[k] - frc * t11_x_xyyz[k] + f2t * (t20_0_xyyz[k] - t21_0_xyyz[k] + tk0_x_yyz[k] - tk1_x_yyz[k]);

                    t_xx_xyzz[k] = fra * t10_x_xyzz[k] - frc * t11_x_xyzz[k] + f2t * (t20_0_xyzz[k] - t21_0_xyzz[k] + tk0_x_yzz[k] - tk1_x_yzz[k]);

                    t_xx_xzzz[k] = fra * t10_x_xzzz[k] - frc * t11_x_xzzz[k] + f2t * (t20_0_xzzz[k] - t21_0_xzzz[k] + tk0_x_zzz[k] - tk1_x_zzz[k]);

                    t_xx_yyyy[k] = fra * t10_x_yyyy[k] - frc * t11_x_yyyy[k] + f2t * (t20_0_yyyy[k] - t21_0_yyyy[k]);

                    t_xx_yyyz[k] = fra * t10_x_yyyz[k] - frc * t11_x_yyyz[k] + f2t * (t20_0_yyyz[k] - t21_0_yyyz[k]);

                    t_xx_yyzz[k] = fra * t10_x_yyzz[k] - frc * t11_x_yyzz[k] + f2t * (t20_0_yyzz[k] - t21_0_yyzz[k]);

                    t_xx_yzzz[k] = fra * t10_x_yzzz[k] - frc * t11_x_yzzz[k] + f2t * (t20_0_yzzz[k] - t21_0_yzzz[k]);

                    t_xx_zzzz[k] = fra * t10_x_zzzz[k] - frc * t11_x_zzzz[k] + f2t * (t20_0_zzzz[k] - t21_0_zzzz[k]);

                    t_xy_xxxx[k] = fra * t10_y_xxxx[k] - frc * t11_y_xxxx[k] + f2t * (4.0 * tk0_y_xxx[k] - 4.0 * tk1_y_xxx[k]);

                    t_xy_xxxy[k] = fra * t10_y_xxxy[k] - frc * t11_y_xxxy[k] + f2t * (3.0 * tk0_y_xxy[k] - 3.0 * tk1_y_xxy[k]);

                    t_xy_xxxz[k] = fra * t10_y_xxxz[k] - frc * t11_y_xxxz[k] + f2t * (3.0 * tk0_y_xxz[k] - 3.0 * tk1_y_xxz[k]);

                    t_xy_xxyy[k] = fra * t10_y_xxyy[k] - frc * t11_y_xxyy[k] + f2t * (2.0 * tk0_y_xyy[k] - 2.0 * tk1_y_xyy[k]);

                    t_xy_xxyz[k] = fra * t10_y_xxyz[k] - frc * t11_y_xxyz[k] + f2t * (2.0 * tk0_y_xyz[k] - 2.0 * tk1_y_xyz[k]);

                    t_xy_xxzz[k] = fra * t10_y_xxzz[k] - frc * t11_y_xxzz[k] + f2t * (2.0 * tk0_y_xzz[k] - 2.0 * tk1_y_xzz[k]);

                    t_xy_xyyy[k] = fra * t10_y_xyyy[k] - frc * t11_y_xyyy[k] + f2t * (tk0_y_yyy[k] - tk1_y_yyy[k]);

                    t_xy_xyyz[k] = fra * t10_y_xyyz[k] - frc * t11_y_xyyz[k] + f2t * (tk0_y_yyz[k] - tk1_y_yyz[k]);

                    t_xy_xyzz[k] = fra * t10_y_xyzz[k] - frc * t11_y_xyzz[k] + f2t * (tk0_y_yzz[k] - tk1_y_yzz[k]);

                    t_xy_xzzz[k] = fra * t10_y_xzzz[k] - frc * t11_y_xzzz[k] + f2t * (tk0_y_zzz[k] - tk1_y_zzz[k]);

                    t_xy_yyyy[k] = fra * t10_y_yyyy[k] - frc * t11_y_yyyy[k];

                    t_xy_yyyz[k] = fra * t10_y_yyyz[k] - frc * t11_y_yyyz[k];

                    t_xy_yyzz[k] = fra * t10_y_yyzz[k] - frc * t11_y_yyzz[k];

                    t_xy_yzzz[k] = fra * t10_y_yzzz[k] - frc * t11_y_yzzz[k];

                    t_xy_zzzz[k] = fra * t10_y_zzzz[k] - frc * t11_y_zzzz[k];

                    t_xz_xxxx[k] = fra * t10_z_xxxx[k] - frc * t11_z_xxxx[k] + f2t * (4.0 * tk0_z_xxx[k] - 4.0 * tk1_z_xxx[k]);

                    t_xz_xxxy[k] = fra * t10_z_xxxy[k] - frc * t11_z_xxxy[k] + f2t * (3.0 * tk0_z_xxy[k] - 3.0 * tk1_z_xxy[k]);

                    t_xz_xxxz[k] = fra * t10_z_xxxz[k] - frc * t11_z_xxxz[k] + f2t * (3.0 * tk0_z_xxz[k] - 3.0 * tk1_z_xxz[k]);

                    t_xz_xxyy[k] = fra * t10_z_xxyy[k] - frc * t11_z_xxyy[k] + f2t * (2.0 * tk0_z_xyy[k] - 2.0 * tk1_z_xyy[k]);

                    t_xz_xxyz[k] = fra * t10_z_xxyz[k] - frc * t11_z_xxyz[k] + f2t * (2.0 * tk0_z_xyz[k] - 2.0 * tk1_z_xyz[k]);

                    t_xz_xxzz[k] = fra * t10_z_xxzz[k] - frc * t11_z_xxzz[k] + f2t * (2.0 * tk0_z_xzz[k] - 2.0 * tk1_z_xzz[k]);

                    t_xz_xyyy[k] = fra * t10_z_xyyy[k] - frc * t11_z_xyyy[k] + f2t * (tk0_z_yyy[k] - tk1_z_yyy[k]);

                    t_xz_xyyz[k] = fra * t10_z_xyyz[k] - frc * t11_z_xyyz[k] + f2t * (tk0_z_yyz[k] - tk1_z_yyz[k]);

                    t_xz_xyzz[k] = fra * t10_z_xyzz[k] - frc * t11_z_xyzz[k] + f2t * (tk0_z_yzz[k] - tk1_z_yzz[k]);

                    t_xz_xzzz[k] = fra * t10_z_xzzz[k] - frc * t11_z_xzzz[k] + f2t * (tk0_z_zzz[k] - tk1_z_zzz[k]);

                    t_xz_yyyy[k] = fra * t10_z_yyyy[k] - frc * t11_z_yyyy[k];

                    t_xz_yyyz[k] = fra * t10_z_yyyz[k] - frc * t11_z_yyyz[k];

                    t_xz_yyzz[k] = fra * t10_z_yyzz[k] - frc * t11_z_yyzz[k];

                    t_xz_yzzz[k] = fra * t10_z_yzzz[k] - frc * t11_z_yzzz[k];

                    t_xz_zzzz[k] = fra * t10_z_zzzz[k] - frc * t11_z_zzzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yy_xxxx[k] = fra * t10_y_xxxx[k] - frc * t11_y_xxxx[k] + f2t * (t20_0_xxxx[k] - t21_0_xxxx[k]);

                    t_yy_xxxy[k] = fra * t10_y_xxxy[k] - frc * t11_y_xxxy[k] + f2t * (t20_0_xxxy[k] - t21_0_xxxy[k] + tk0_y_xxx[k] - tk1_y_xxx[k]);

                    t_yy_xxxz[k] = fra * t10_y_xxxz[k] - frc * t11_y_xxxz[k] + f2t * (t20_0_xxxz[k] - t21_0_xxxz[k]);

                    t_yy_xxyy[k] = fra * t10_y_xxyy[k] - frc * t11_y_xxyy[k] + f2t * (t20_0_xxyy[k] - t21_0_xxyy[k] + 2.0 * tk0_y_xxy[k] - 2.0 * tk1_y_xxy[k]);

                    t_yy_xxyz[k] = fra * t10_y_xxyz[k] - frc * t11_y_xxyz[k] + f2t * (t20_0_xxyz[k] - t21_0_xxyz[k] + tk0_y_xxz[k] - tk1_y_xxz[k]);

                    t_yy_xxzz[k] = fra * t10_y_xxzz[k] - frc * t11_y_xxzz[k] + f2t * (t20_0_xxzz[k] - t21_0_xxzz[k]);

                    t_yy_xyyy[k] = fra * t10_y_xyyy[k] - frc * t11_y_xyyy[k] + f2t * (t20_0_xyyy[k] - t21_0_xyyy[k] + 3.0 * tk0_y_xyy[k] - 3.0 * tk1_y_xyy[k]);

                    t_yy_xyyz[k] = fra * t10_y_xyyz[k] - frc * t11_y_xyyz[k] + f2t * (t20_0_xyyz[k] - t21_0_xyyz[k] + 2.0 * tk0_y_xyz[k] - 2.0 * tk1_y_xyz[k]);

                    t_yy_xyzz[k] = fra * t10_y_xyzz[k] - frc * t11_y_xyzz[k] + f2t * (t20_0_xyzz[k] - t21_0_xyzz[k] + tk0_y_xzz[k] - tk1_y_xzz[k]);

                    t_yy_xzzz[k] = fra * t10_y_xzzz[k] - frc * t11_y_xzzz[k] + f2t * (t20_0_xzzz[k] - t21_0_xzzz[k]);

                    t_yy_yyyy[k] = fra * t10_y_yyyy[k] - frc * t11_y_yyyy[k] + f2t * (t20_0_yyyy[k] - t21_0_yyyy[k] + 4.0 * tk0_y_yyy[k] - 4.0 * tk1_y_yyy[k]);

                    t_yy_yyyz[k] = fra * t10_y_yyyz[k] - frc * t11_y_yyyz[k] + f2t * (t20_0_yyyz[k] - t21_0_yyyz[k] + 3.0 * tk0_y_yyz[k] - 3.0 * tk1_y_yyz[k]);

                    t_yy_yyzz[k] = fra * t10_y_yyzz[k] - frc * t11_y_yyzz[k] + f2t * (t20_0_yyzz[k] - t21_0_yyzz[k] + 2.0 * tk0_y_yzz[k] - 2.0 * tk1_y_yzz[k]);

                    t_yy_yzzz[k] = fra * t10_y_yzzz[k] - frc * t11_y_yzzz[k] + f2t * (t20_0_yzzz[k] - t21_0_yzzz[k] + tk0_y_zzz[k] - tk1_y_zzz[k]);

                    t_yy_zzzz[k] = fra * t10_y_zzzz[k] - frc * t11_y_zzzz[k] + f2t * (t20_0_zzzz[k] - t21_0_zzzz[k]);

                    t_yz_xxxx[k] = fra * t10_z_xxxx[k] - frc * t11_z_xxxx[k];

                    t_yz_xxxy[k] = fra * t10_z_xxxy[k] - frc * t11_z_xxxy[k] + f2t * (tk0_z_xxx[k] - tk1_z_xxx[k]);

                    t_yz_xxxz[k] = fra * t10_z_xxxz[k] - frc * t11_z_xxxz[k];

                    t_yz_xxyy[k] = fra * t10_z_xxyy[k] - frc * t11_z_xxyy[k] + f2t * (2.0 * tk0_z_xxy[k] - 2.0 * tk1_z_xxy[k]);

                    t_yz_xxyz[k] = fra * t10_z_xxyz[k] - frc * t11_z_xxyz[k] + f2t * (tk0_z_xxz[k] - tk1_z_xxz[k]);

                    t_yz_xxzz[k] = fra * t10_z_xxzz[k] - frc * t11_z_xxzz[k];

                    t_yz_xyyy[k] = fra * t10_z_xyyy[k] - frc * t11_z_xyyy[k] + f2t * (3.0 * tk0_z_xyy[k] - 3.0 * tk1_z_xyy[k]);

                    t_yz_xyyz[k] = fra * t10_z_xyyz[k] - frc * t11_z_xyyz[k] + f2t * (2.0 * tk0_z_xyz[k] - 2.0 * tk1_z_xyz[k]);

                    t_yz_xyzz[k] = fra * t10_z_xyzz[k] - frc * t11_z_xyzz[k] + f2t * (tk0_z_xzz[k] - tk1_z_xzz[k]);

                    t_yz_xzzz[k] = fra * t10_z_xzzz[k] - frc * t11_z_xzzz[k];

                    t_yz_yyyy[k] = fra * t10_z_yyyy[k] - frc * t11_z_yyyy[k] + f2t * (4.0 * tk0_z_yyy[k] - 4.0 * tk1_z_yyy[k]);

                    t_yz_yyyz[k] = fra * t10_z_yyyz[k] - frc * t11_z_yyyz[k] + f2t * (3.0 * tk0_z_yyz[k] - 3.0 * tk1_z_yyz[k]);

                    t_yz_yyzz[k] = fra * t10_z_yyzz[k] - frc * t11_z_yyzz[k] + f2t * (2.0 * tk0_z_yzz[k] - 2.0 * tk1_z_yzz[k]);

                    t_yz_yzzz[k] = fra * t10_z_yzzz[k] - frc * t11_z_yzzz[k] + f2t * (tk0_z_zzz[k] - tk1_z_zzz[k]);

                    t_yz_zzzz[k] = fra * t10_z_zzzz[k] - frc * t11_z_zzzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zz_xxxx[k] = fra * t10_z_xxxx[k] - frc * t11_z_xxxx[k] + f2t * (t20_0_xxxx[k] - t21_0_xxxx[k]);

                    t_zz_xxxy[k] = fra * t10_z_xxxy[k] - frc * t11_z_xxxy[k] + f2t * (t20_0_xxxy[k] - t21_0_xxxy[k]);

                    t_zz_xxxz[k] = fra * t10_z_xxxz[k] - frc * t11_z_xxxz[k] + f2t * (t20_0_xxxz[k] - t21_0_xxxz[k] + tk0_z_xxx[k] - tk1_z_xxx[k]);

                    t_zz_xxyy[k] = fra * t10_z_xxyy[k] - frc * t11_z_xxyy[k] + f2t * (t20_0_xxyy[k] - t21_0_xxyy[k]);

                    t_zz_xxyz[k] = fra * t10_z_xxyz[k] - frc * t11_z_xxyz[k] + f2t * (t20_0_xxyz[k] - t21_0_xxyz[k] + tk0_z_xxy[k] - tk1_z_xxy[k]);

                    t_zz_xxzz[k] = fra * t10_z_xxzz[k] - frc * t11_z_xxzz[k] + f2t * (t20_0_xxzz[k] - t21_0_xxzz[k] + 2.0 * tk0_z_xxz[k] - 2.0 * tk1_z_xxz[k]);

                    t_zz_xyyy[k] = fra * t10_z_xyyy[k] - frc * t11_z_xyyy[k] + f2t * (t20_0_xyyy[k] - t21_0_xyyy[k]);

                    t_zz_xyyz[k] = fra * t10_z_xyyz[k] - frc * t11_z_xyyz[k] + f2t * (t20_0_xyyz[k] - t21_0_xyyz[k] + tk0_z_xyy[k] - tk1_z_xyy[k]);

                    t_zz_xyzz[k] = fra * t10_z_xyzz[k] - frc * t11_z_xyzz[k] + f2t * (t20_0_xyzz[k] - t21_0_xyzz[k] + 2.0 * tk0_z_xyz[k] - 2.0 * tk1_z_xyz[k]);

                    t_zz_xzzz[k] = fra * t10_z_xzzz[k] - frc * t11_z_xzzz[k] + f2t * (t20_0_xzzz[k] - t21_0_xzzz[k] + 3.0 * tk0_z_xzz[k] - 3.0 * tk1_z_xzz[k]);

                    t_zz_yyyy[k] = fra * t10_z_yyyy[k] - frc * t11_z_yyyy[k] + f2t * (t20_0_yyyy[k] - t21_0_yyyy[k]);

                    t_zz_yyyz[k] = fra * t10_z_yyyz[k] - frc * t11_z_yyyz[k] + f2t * (t20_0_yyyz[k] - t21_0_yyyz[k] + tk0_z_yyy[k] - tk1_z_yyy[k]);

                    t_zz_yyzz[k] = fra * t10_z_yyzz[k] - frc * t11_z_yyzz[k] + f2t * (t20_0_yyzz[k] - t21_0_yyzz[k] + 2.0 * tk0_z_yyz[k] - 2.0 * tk1_z_yyz[k]);

                    t_zz_yzzz[k] = fra * t10_z_yzzz[k] - frc * t11_z_yzzz[k] + f2t * (t20_0_yzzz[k] - t21_0_yzzz[k] + 3.0 * tk0_z_yzz[k] - 3.0 * tk1_z_yzz[k]);

                    t_zz_zzzz[k] = fra * t10_z_zzzz[k] - frc * t11_z_zzzz[k] + f2t * (t20_0_zzzz[k] - t21_0_zzzz[k] + 4.0 * tk0_z_zzz[k] - 4.0 * tk1_z_zzz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForGD(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {4, 2, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 4, 2);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {4, 2, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 2, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 2, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 2, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 1, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 1, i + 1});

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

                // set up pointers to (F|A(0)|P)^(m) integrals

                auto tk0_xxx_x = primBuffer.data(tk0off + 30 * idx);

                auto tk0_xxx_y = primBuffer.data(tk0off + 30 * idx + 1);

                auto tk0_xxx_z = primBuffer.data(tk0off + 30 * idx + 2);

                auto tk0_xxy_x = primBuffer.data(tk0off + 30 * idx + 3);

                auto tk0_xxy_y = primBuffer.data(tk0off + 30 * idx + 4);

                auto tk0_xxy_z = primBuffer.data(tk0off + 30 * idx + 5);

                auto tk0_xxz_x = primBuffer.data(tk0off + 30 * idx + 6);

                auto tk0_xxz_y = primBuffer.data(tk0off + 30 * idx + 7);

                auto tk0_xxz_z = primBuffer.data(tk0off + 30 * idx + 8);

                auto tk0_xyy_x = primBuffer.data(tk0off + 30 * idx + 9);

                auto tk0_xyy_y = primBuffer.data(tk0off + 30 * idx + 10);

                auto tk0_xyy_z = primBuffer.data(tk0off + 30 * idx + 11);

                auto tk0_xyz_x = primBuffer.data(tk0off + 30 * idx + 12);

                auto tk0_xyz_y = primBuffer.data(tk0off + 30 * idx + 13);

                auto tk0_xyz_z = primBuffer.data(tk0off + 30 * idx + 14);

                auto tk0_xzz_x = primBuffer.data(tk0off + 30 * idx + 15);

                auto tk0_xzz_y = primBuffer.data(tk0off + 30 * idx + 16);

                auto tk0_xzz_z = primBuffer.data(tk0off + 30 * idx + 17);

                auto tk0_yyy_x = primBuffer.data(tk0off + 30 * idx + 18);

                auto tk0_yyy_y = primBuffer.data(tk0off + 30 * idx + 19);

                auto tk0_yyy_z = primBuffer.data(tk0off + 30 * idx + 20);

                auto tk0_yyz_x = primBuffer.data(tk0off + 30 * idx + 21);

                auto tk0_yyz_y = primBuffer.data(tk0off + 30 * idx + 22);

                auto tk0_yyz_z = primBuffer.data(tk0off + 30 * idx + 23);

                auto tk0_yzz_x = primBuffer.data(tk0off + 30 * idx + 24);

                auto tk0_yzz_y = primBuffer.data(tk0off + 30 * idx + 25);

                auto tk0_yzz_z = primBuffer.data(tk0off + 30 * idx + 26);

                auto tk0_zzz_x = primBuffer.data(tk0off + 30 * idx + 27);

                auto tk0_zzz_y = primBuffer.data(tk0off + 30 * idx + 28);

                auto tk0_zzz_z = primBuffer.data(tk0off + 30 * idx + 29);

                // set up pointers to (F|A(0)|P)^(m+1) integrals

                auto tk1_xxx_x = primBuffer.data(tk1off + 30 * idx);

                auto tk1_xxx_y = primBuffer.data(tk1off + 30 * idx + 1);

                auto tk1_xxx_z = primBuffer.data(tk1off + 30 * idx + 2);

                auto tk1_xxy_x = primBuffer.data(tk1off + 30 * idx + 3);

                auto tk1_xxy_y = primBuffer.data(tk1off + 30 * idx + 4);

                auto tk1_xxy_z = primBuffer.data(tk1off + 30 * idx + 5);

                auto tk1_xxz_x = primBuffer.data(tk1off + 30 * idx + 6);

                auto tk1_xxz_y = primBuffer.data(tk1off + 30 * idx + 7);

                auto tk1_xxz_z = primBuffer.data(tk1off + 30 * idx + 8);

                auto tk1_xyy_x = primBuffer.data(tk1off + 30 * idx + 9);

                auto tk1_xyy_y = primBuffer.data(tk1off + 30 * idx + 10);

                auto tk1_xyy_z = primBuffer.data(tk1off + 30 * idx + 11);

                auto tk1_xyz_x = primBuffer.data(tk1off + 30 * idx + 12);

                auto tk1_xyz_y = primBuffer.data(tk1off + 30 * idx + 13);

                auto tk1_xyz_z = primBuffer.data(tk1off + 30 * idx + 14);

                auto tk1_xzz_x = primBuffer.data(tk1off + 30 * idx + 15);

                auto tk1_xzz_y = primBuffer.data(tk1off + 30 * idx + 16);

                auto tk1_xzz_z = primBuffer.data(tk1off + 30 * idx + 17);

                auto tk1_yyy_x = primBuffer.data(tk1off + 30 * idx + 18);

                auto tk1_yyy_y = primBuffer.data(tk1off + 30 * idx + 19);

                auto tk1_yyy_z = primBuffer.data(tk1off + 30 * idx + 20);

                auto tk1_yyz_x = primBuffer.data(tk1off + 30 * idx + 21);

                auto tk1_yyz_y = primBuffer.data(tk1off + 30 * idx + 22);

                auto tk1_yyz_z = primBuffer.data(tk1off + 30 * idx + 23);

                auto tk1_yzz_x = primBuffer.data(tk1off + 30 * idx + 24);

                auto tk1_yzz_y = primBuffer.data(tk1off + 30 * idx + 25);

                auto tk1_yzz_z = primBuffer.data(tk1off + 30 * idx + 26);

                auto tk1_zzz_x = primBuffer.data(tk1off + 30 * idx + 27);

                auto tk1_zzz_y = primBuffer.data(tk1off + 30 * idx + 28);

                auto tk1_zzz_z = primBuffer.data(tk1off + 30 * idx + 29);

                // set up pointers to (D|A(0)|D)^(m) integrals

                auto t20_xx_xx = primBuffer.data(t20off + 36 * idx);

                auto t20_xx_xy = primBuffer.data(t20off + 36 * idx + 1);

                auto t20_xx_xz = primBuffer.data(t20off + 36 * idx + 2);

                auto t20_xx_yy = primBuffer.data(t20off + 36 * idx + 3);

                auto t20_xx_yz = primBuffer.data(t20off + 36 * idx + 4);

                auto t20_xx_zz = primBuffer.data(t20off + 36 * idx + 5);

                auto t20_xy_xx = primBuffer.data(t20off + 36 * idx + 6);

                auto t20_xy_xy = primBuffer.data(t20off + 36 * idx + 7);

                auto t20_xy_xz = primBuffer.data(t20off + 36 * idx + 8);

                auto t20_xy_yy = primBuffer.data(t20off + 36 * idx + 9);

                auto t20_xy_yz = primBuffer.data(t20off + 36 * idx + 10);

                auto t20_xy_zz = primBuffer.data(t20off + 36 * idx + 11);

                auto t20_xz_xx = primBuffer.data(t20off + 36 * idx + 12);

                auto t20_xz_xy = primBuffer.data(t20off + 36 * idx + 13);

                auto t20_xz_xz = primBuffer.data(t20off + 36 * idx + 14);

                auto t20_xz_yy = primBuffer.data(t20off + 36 * idx + 15);

                auto t20_xz_yz = primBuffer.data(t20off + 36 * idx + 16);

                auto t20_xz_zz = primBuffer.data(t20off + 36 * idx + 17);

                auto t20_yy_xx = primBuffer.data(t20off + 36 * idx + 18);

                auto t20_yy_xy = primBuffer.data(t20off + 36 * idx + 19);

                auto t20_yy_xz = primBuffer.data(t20off + 36 * idx + 20);

                auto t20_yy_yy = primBuffer.data(t20off + 36 * idx + 21);

                auto t20_yy_yz = primBuffer.data(t20off + 36 * idx + 22);

                auto t20_yy_zz = primBuffer.data(t20off + 36 * idx + 23);

                auto t20_yz_xx = primBuffer.data(t20off + 36 * idx + 24);

                auto t20_yz_xy = primBuffer.data(t20off + 36 * idx + 25);

                auto t20_yz_xz = primBuffer.data(t20off + 36 * idx + 26);

                auto t20_yz_yy = primBuffer.data(t20off + 36 * idx + 27);

                auto t20_yz_yz = primBuffer.data(t20off + 36 * idx + 28);

                auto t20_yz_zz = primBuffer.data(t20off + 36 * idx + 29);

                auto t20_zz_xx = primBuffer.data(t20off + 36 * idx + 30);

                auto t20_zz_xy = primBuffer.data(t20off + 36 * idx + 31);

                auto t20_zz_xz = primBuffer.data(t20off + 36 * idx + 32);

                auto t20_zz_yy = primBuffer.data(t20off + 36 * idx + 33);

                auto t20_zz_yz = primBuffer.data(t20off + 36 * idx + 34);

                auto t20_zz_zz = primBuffer.data(t20off + 36 * idx + 35);

                // set up pointers to (D|A(0)|D)^(m+1) integrals

                auto t21_xx_xx = primBuffer.data(t21off + 36 * idx);

                auto t21_xx_xy = primBuffer.data(t21off + 36 * idx + 1);

                auto t21_xx_xz = primBuffer.data(t21off + 36 * idx + 2);

                auto t21_xx_yy = primBuffer.data(t21off + 36 * idx + 3);

                auto t21_xx_yz = primBuffer.data(t21off + 36 * idx + 4);

                auto t21_xx_zz = primBuffer.data(t21off + 36 * idx + 5);

                auto t21_xy_xx = primBuffer.data(t21off + 36 * idx + 6);

                auto t21_xy_xy = primBuffer.data(t21off + 36 * idx + 7);

                auto t21_xy_xz = primBuffer.data(t21off + 36 * idx + 8);

                auto t21_xy_yy = primBuffer.data(t21off + 36 * idx + 9);

                auto t21_xy_yz = primBuffer.data(t21off + 36 * idx + 10);

                auto t21_xy_zz = primBuffer.data(t21off + 36 * idx + 11);

                auto t21_xz_xx = primBuffer.data(t21off + 36 * idx + 12);

                auto t21_xz_xy = primBuffer.data(t21off + 36 * idx + 13);

                auto t21_xz_xz = primBuffer.data(t21off + 36 * idx + 14);

                auto t21_xz_yy = primBuffer.data(t21off + 36 * idx + 15);

                auto t21_xz_yz = primBuffer.data(t21off + 36 * idx + 16);

                auto t21_xz_zz = primBuffer.data(t21off + 36 * idx + 17);

                auto t21_yy_xx = primBuffer.data(t21off + 36 * idx + 18);

                auto t21_yy_xy = primBuffer.data(t21off + 36 * idx + 19);

                auto t21_yy_xz = primBuffer.data(t21off + 36 * idx + 20);

                auto t21_yy_yy = primBuffer.data(t21off + 36 * idx + 21);

                auto t21_yy_yz = primBuffer.data(t21off + 36 * idx + 22);

                auto t21_yy_zz = primBuffer.data(t21off + 36 * idx + 23);

                auto t21_yz_xx = primBuffer.data(t21off + 36 * idx + 24);

                auto t21_yz_xy = primBuffer.data(t21off + 36 * idx + 25);

                auto t21_yz_xz = primBuffer.data(t21off + 36 * idx + 26);

                auto t21_yz_yy = primBuffer.data(t21off + 36 * idx + 27);

                auto t21_yz_yz = primBuffer.data(t21off + 36 * idx + 28);

                auto t21_yz_zz = primBuffer.data(t21off + 36 * idx + 29);

                auto t21_zz_xx = primBuffer.data(t21off + 36 * idx + 30);

                auto t21_zz_xy = primBuffer.data(t21off + 36 * idx + 31);

                auto t21_zz_xz = primBuffer.data(t21off + 36 * idx + 32);

                auto t21_zz_yy = primBuffer.data(t21off + 36 * idx + 33);

                auto t21_zz_yz = primBuffer.data(t21off + 36 * idx + 34);

                auto t21_zz_zz = primBuffer.data(t21off + 36 * idx + 35);

                // set up pointers to (F|A(0)|D)^(m) integrals

                auto t10_xxx_xx = primBuffer.data(t10off + 60 * idx);

                auto t10_xxx_xy = primBuffer.data(t10off + 60 * idx + 1);

                auto t10_xxx_xz = primBuffer.data(t10off + 60 * idx + 2);

                auto t10_xxx_yy = primBuffer.data(t10off + 60 * idx + 3);

                auto t10_xxx_yz = primBuffer.data(t10off + 60 * idx + 4);

                auto t10_xxx_zz = primBuffer.data(t10off + 60 * idx + 5);

                auto t10_xxy_xx = primBuffer.data(t10off + 60 * idx + 6);

                auto t10_xxy_xy = primBuffer.data(t10off + 60 * idx + 7);

                auto t10_xxy_xz = primBuffer.data(t10off + 60 * idx + 8);

                auto t10_xxy_yy = primBuffer.data(t10off + 60 * idx + 9);

                auto t10_xxy_yz = primBuffer.data(t10off + 60 * idx + 10);

                auto t10_xxy_zz = primBuffer.data(t10off + 60 * idx + 11);

                auto t10_xxz_xx = primBuffer.data(t10off + 60 * idx + 12);

                auto t10_xxz_xy = primBuffer.data(t10off + 60 * idx + 13);

                auto t10_xxz_xz = primBuffer.data(t10off + 60 * idx + 14);

                auto t10_xxz_yy = primBuffer.data(t10off + 60 * idx + 15);

                auto t10_xxz_yz = primBuffer.data(t10off + 60 * idx + 16);

                auto t10_xxz_zz = primBuffer.data(t10off + 60 * idx + 17);

                auto t10_xyy_xx = primBuffer.data(t10off + 60 * idx + 18);

                auto t10_xyy_xy = primBuffer.data(t10off + 60 * idx + 19);

                auto t10_xyy_xz = primBuffer.data(t10off + 60 * idx + 20);

                auto t10_xyy_yy = primBuffer.data(t10off + 60 * idx + 21);

                auto t10_xyy_yz = primBuffer.data(t10off + 60 * idx + 22);

                auto t10_xyy_zz = primBuffer.data(t10off + 60 * idx + 23);

                auto t10_xyz_xx = primBuffer.data(t10off + 60 * idx + 24);

                auto t10_xyz_xy = primBuffer.data(t10off + 60 * idx + 25);

                auto t10_xyz_xz = primBuffer.data(t10off + 60 * idx + 26);

                auto t10_xyz_yy = primBuffer.data(t10off + 60 * idx + 27);

                auto t10_xyz_yz = primBuffer.data(t10off + 60 * idx + 28);

                auto t10_xyz_zz = primBuffer.data(t10off + 60 * idx + 29);

                auto t10_xzz_xx = primBuffer.data(t10off + 60 * idx + 30);

                auto t10_xzz_xy = primBuffer.data(t10off + 60 * idx + 31);

                auto t10_xzz_xz = primBuffer.data(t10off + 60 * idx + 32);

                auto t10_xzz_yy = primBuffer.data(t10off + 60 * idx + 33);

                auto t10_xzz_yz = primBuffer.data(t10off + 60 * idx + 34);

                auto t10_xzz_zz = primBuffer.data(t10off + 60 * idx + 35);

                auto t10_yyy_xx = primBuffer.data(t10off + 60 * idx + 36);

                auto t10_yyy_xy = primBuffer.data(t10off + 60 * idx + 37);

                auto t10_yyy_xz = primBuffer.data(t10off + 60 * idx + 38);

                auto t10_yyy_yy = primBuffer.data(t10off + 60 * idx + 39);

                auto t10_yyy_yz = primBuffer.data(t10off + 60 * idx + 40);

                auto t10_yyy_zz = primBuffer.data(t10off + 60 * idx + 41);

                auto t10_yyz_xx = primBuffer.data(t10off + 60 * idx + 42);

                auto t10_yyz_xy = primBuffer.data(t10off + 60 * idx + 43);

                auto t10_yyz_xz = primBuffer.data(t10off + 60 * idx + 44);

                auto t10_yyz_yy = primBuffer.data(t10off + 60 * idx + 45);

                auto t10_yyz_yz = primBuffer.data(t10off + 60 * idx + 46);

                auto t10_yyz_zz = primBuffer.data(t10off + 60 * idx + 47);

                auto t10_yzz_xx = primBuffer.data(t10off + 60 * idx + 48);

                auto t10_yzz_xy = primBuffer.data(t10off + 60 * idx + 49);

                auto t10_yzz_xz = primBuffer.data(t10off + 60 * idx + 50);

                auto t10_yzz_yy = primBuffer.data(t10off + 60 * idx + 51);

                auto t10_yzz_yz = primBuffer.data(t10off + 60 * idx + 52);

                auto t10_yzz_zz = primBuffer.data(t10off + 60 * idx + 53);

                auto t10_zzz_xx = primBuffer.data(t10off + 60 * idx + 54);

                auto t10_zzz_xy = primBuffer.data(t10off + 60 * idx + 55);

                auto t10_zzz_xz = primBuffer.data(t10off + 60 * idx + 56);

                auto t10_zzz_yy = primBuffer.data(t10off + 60 * idx + 57);

                auto t10_zzz_yz = primBuffer.data(t10off + 60 * idx + 58);

                auto t10_zzz_zz = primBuffer.data(t10off + 60 * idx + 59);

                // set up pointers to (F|A(0)|D)^(m+1) integrals

                auto t11_xxx_xx = primBuffer.data(t11off + 60 * idx);

                auto t11_xxx_xy = primBuffer.data(t11off + 60 * idx + 1);

                auto t11_xxx_xz = primBuffer.data(t11off + 60 * idx + 2);

                auto t11_xxx_yy = primBuffer.data(t11off + 60 * idx + 3);

                auto t11_xxx_yz = primBuffer.data(t11off + 60 * idx + 4);

                auto t11_xxx_zz = primBuffer.data(t11off + 60 * idx + 5);

                auto t11_xxy_xx = primBuffer.data(t11off + 60 * idx + 6);

                auto t11_xxy_xy = primBuffer.data(t11off + 60 * idx + 7);

                auto t11_xxy_xz = primBuffer.data(t11off + 60 * idx + 8);

                auto t11_xxy_yy = primBuffer.data(t11off + 60 * idx + 9);

                auto t11_xxy_yz = primBuffer.data(t11off + 60 * idx + 10);

                auto t11_xxy_zz = primBuffer.data(t11off + 60 * idx + 11);

                auto t11_xxz_xx = primBuffer.data(t11off + 60 * idx + 12);

                auto t11_xxz_xy = primBuffer.data(t11off + 60 * idx + 13);

                auto t11_xxz_xz = primBuffer.data(t11off + 60 * idx + 14);

                auto t11_xxz_yy = primBuffer.data(t11off + 60 * idx + 15);

                auto t11_xxz_yz = primBuffer.data(t11off + 60 * idx + 16);

                auto t11_xxz_zz = primBuffer.data(t11off + 60 * idx + 17);

                auto t11_xyy_xx = primBuffer.data(t11off + 60 * idx + 18);

                auto t11_xyy_xy = primBuffer.data(t11off + 60 * idx + 19);

                auto t11_xyy_xz = primBuffer.data(t11off + 60 * idx + 20);

                auto t11_xyy_yy = primBuffer.data(t11off + 60 * idx + 21);

                auto t11_xyy_yz = primBuffer.data(t11off + 60 * idx + 22);

                auto t11_xyy_zz = primBuffer.data(t11off + 60 * idx + 23);

                auto t11_xyz_xx = primBuffer.data(t11off + 60 * idx + 24);

                auto t11_xyz_xy = primBuffer.data(t11off + 60 * idx + 25);

                auto t11_xyz_xz = primBuffer.data(t11off + 60 * idx + 26);

                auto t11_xyz_yy = primBuffer.data(t11off + 60 * idx + 27);

                auto t11_xyz_yz = primBuffer.data(t11off + 60 * idx + 28);

                auto t11_xyz_zz = primBuffer.data(t11off + 60 * idx + 29);

                auto t11_xzz_xx = primBuffer.data(t11off + 60 * idx + 30);

                auto t11_xzz_xy = primBuffer.data(t11off + 60 * idx + 31);

                auto t11_xzz_xz = primBuffer.data(t11off + 60 * idx + 32);

                auto t11_xzz_yy = primBuffer.data(t11off + 60 * idx + 33);

                auto t11_xzz_yz = primBuffer.data(t11off + 60 * idx + 34);

                auto t11_xzz_zz = primBuffer.data(t11off + 60 * idx + 35);

                auto t11_yyy_xx = primBuffer.data(t11off + 60 * idx + 36);

                auto t11_yyy_xy = primBuffer.data(t11off + 60 * idx + 37);

                auto t11_yyy_xz = primBuffer.data(t11off + 60 * idx + 38);

                auto t11_yyy_yy = primBuffer.data(t11off + 60 * idx + 39);

                auto t11_yyy_yz = primBuffer.data(t11off + 60 * idx + 40);

                auto t11_yyy_zz = primBuffer.data(t11off + 60 * idx + 41);

                auto t11_yyz_xx = primBuffer.data(t11off + 60 * idx + 42);

                auto t11_yyz_xy = primBuffer.data(t11off + 60 * idx + 43);

                auto t11_yyz_xz = primBuffer.data(t11off + 60 * idx + 44);

                auto t11_yyz_yy = primBuffer.data(t11off + 60 * idx + 45);

                auto t11_yyz_yz = primBuffer.data(t11off + 60 * idx + 46);

                auto t11_yyz_zz = primBuffer.data(t11off + 60 * idx + 47);

                auto t11_yzz_xx = primBuffer.data(t11off + 60 * idx + 48);

                auto t11_yzz_xy = primBuffer.data(t11off + 60 * idx + 49);

                auto t11_yzz_xz = primBuffer.data(t11off + 60 * idx + 50);

                auto t11_yzz_yy = primBuffer.data(t11off + 60 * idx + 51);

                auto t11_yzz_yz = primBuffer.data(t11off + 60 * idx + 52);

                auto t11_yzz_zz = primBuffer.data(t11off + 60 * idx + 53);

                auto t11_zzz_xx = primBuffer.data(t11off + 60 * idx + 54);

                auto t11_zzz_xy = primBuffer.data(t11off + 60 * idx + 55);

                auto t11_zzz_xz = primBuffer.data(t11off + 60 * idx + 56);

                auto t11_zzz_yy = primBuffer.data(t11off + 60 * idx + 57);

                auto t11_zzz_yz = primBuffer.data(t11off + 60 * idx + 58);

                auto t11_zzz_zz = primBuffer.data(t11off + 60 * idx + 59);

                // set up pointers to (G|A(0)|D)^(m) integrals

                auto t_xxxx_xx = primBuffer.data(toff + 90 * idx);

                auto t_xxxx_xy = primBuffer.data(toff + 90 * idx + 1);

                auto t_xxxx_xz = primBuffer.data(toff + 90 * idx + 2);

                auto t_xxxx_yy = primBuffer.data(toff + 90 * idx + 3);

                auto t_xxxx_yz = primBuffer.data(toff + 90 * idx + 4);

                auto t_xxxx_zz = primBuffer.data(toff + 90 * idx + 5);

                auto t_xxxy_xx = primBuffer.data(toff + 90 * idx + 6);

                auto t_xxxy_xy = primBuffer.data(toff + 90 * idx + 7);

                auto t_xxxy_xz = primBuffer.data(toff + 90 * idx + 8);

                auto t_xxxy_yy = primBuffer.data(toff + 90 * idx + 9);

                auto t_xxxy_yz = primBuffer.data(toff + 90 * idx + 10);

                auto t_xxxy_zz = primBuffer.data(toff + 90 * idx + 11);

                auto t_xxxz_xx = primBuffer.data(toff + 90 * idx + 12);

                auto t_xxxz_xy = primBuffer.data(toff + 90 * idx + 13);

                auto t_xxxz_xz = primBuffer.data(toff + 90 * idx + 14);

                auto t_xxxz_yy = primBuffer.data(toff + 90 * idx + 15);

                auto t_xxxz_yz = primBuffer.data(toff + 90 * idx + 16);

                auto t_xxxz_zz = primBuffer.data(toff + 90 * idx + 17);

                auto t_xxyy_xx = primBuffer.data(toff + 90 * idx + 18);

                auto t_xxyy_xy = primBuffer.data(toff + 90 * idx + 19);

                auto t_xxyy_xz = primBuffer.data(toff + 90 * idx + 20);

                auto t_xxyy_yy = primBuffer.data(toff + 90 * idx + 21);

                auto t_xxyy_yz = primBuffer.data(toff + 90 * idx + 22);

                auto t_xxyy_zz = primBuffer.data(toff + 90 * idx + 23);

                auto t_xxyz_xx = primBuffer.data(toff + 90 * idx + 24);

                auto t_xxyz_xy = primBuffer.data(toff + 90 * idx + 25);

                auto t_xxyz_xz = primBuffer.data(toff + 90 * idx + 26);

                auto t_xxyz_yy = primBuffer.data(toff + 90 * idx + 27);

                auto t_xxyz_yz = primBuffer.data(toff + 90 * idx + 28);

                auto t_xxyz_zz = primBuffer.data(toff + 90 * idx + 29);

                auto t_xxzz_xx = primBuffer.data(toff + 90 * idx + 30);

                auto t_xxzz_xy = primBuffer.data(toff + 90 * idx + 31);

                auto t_xxzz_xz = primBuffer.data(toff + 90 * idx + 32);

                auto t_xxzz_yy = primBuffer.data(toff + 90 * idx + 33);

                auto t_xxzz_yz = primBuffer.data(toff + 90 * idx + 34);

                auto t_xxzz_zz = primBuffer.data(toff + 90 * idx + 35);

                auto t_xyyy_xx = primBuffer.data(toff + 90 * idx + 36);

                auto t_xyyy_xy = primBuffer.data(toff + 90 * idx + 37);

                auto t_xyyy_xz = primBuffer.data(toff + 90 * idx + 38);

                auto t_xyyy_yy = primBuffer.data(toff + 90 * idx + 39);

                auto t_xyyy_yz = primBuffer.data(toff + 90 * idx + 40);

                auto t_xyyy_zz = primBuffer.data(toff + 90 * idx + 41);

                auto t_xyyz_xx = primBuffer.data(toff + 90 * idx + 42);

                auto t_xyyz_xy = primBuffer.data(toff + 90 * idx + 43);

                auto t_xyyz_xz = primBuffer.data(toff + 90 * idx + 44);

                auto t_xyyz_yy = primBuffer.data(toff + 90 * idx + 45);

                auto t_xyyz_yz = primBuffer.data(toff + 90 * idx + 46);

                auto t_xyyz_zz = primBuffer.data(toff + 90 * idx + 47);

                auto t_xyzz_xx = primBuffer.data(toff + 90 * idx + 48);

                auto t_xyzz_xy = primBuffer.data(toff + 90 * idx + 49);

                auto t_xyzz_xz = primBuffer.data(toff + 90 * idx + 50);

                auto t_xyzz_yy = primBuffer.data(toff + 90 * idx + 51);

                auto t_xyzz_yz = primBuffer.data(toff + 90 * idx + 52);

                auto t_xyzz_zz = primBuffer.data(toff + 90 * idx + 53);

                auto t_xzzz_xx = primBuffer.data(toff + 90 * idx + 54);

                auto t_xzzz_xy = primBuffer.data(toff + 90 * idx + 55);

                auto t_xzzz_xz = primBuffer.data(toff + 90 * idx + 56);

                auto t_xzzz_yy = primBuffer.data(toff + 90 * idx + 57);

                auto t_xzzz_yz = primBuffer.data(toff + 90 * idx + 58);

                auto t_xzzz_zz = primBuffer.data(toff + 90 * idx + 59);

                auto t_yyyy_xx = primBuffer.data(toff + 90 * idx + 60);

                auto t_yyyy_xy = primBuffer.data(toff + 90 * idx + 61);

                auto t_yyyy_xz = primBuffer.data(toff + 90 * idx + 62);

                auto t_yyyy_yy = primBuffer.data(toff + 90 * idx + 63);

                auto t_yyyy_yz = primBuffer.data(toff + 90 * idx + 64);

                auto t_yyyy_zz = primBuffer.data(toff + 90 * idx + 65);

                auto t_yyyz_xx = primBuffer.data(toff + 90 * idx + 66);

                auto t_yyyz_xy = primBuffer.data(toff + 90 * idx + 67);

                auto t_yyyz_xz = primBuffer.data(toff + 90 * idx + 68);

                auto t_yyyz_yy = primBuffer.data(toff + 90 * idx + 69);

                auto t_yyyz_yz = primBuffer.data(toff + 90 * idx + 70);

                auto t_yyyz_zz = primBuffer.data(toff + 90 * idx + 71);

                auto t_yyzz_xx = primBuffer.data(toff + 90 * idx + 72);

                auto t_yyzz_xy = primBuffer.data(toff + 90 * idx + 73);

                auto t_yyzz_xz = primBuffer.data(toff + 90 * idx + 74);

                auto t_yyzz_yy = primBuffer.data(toff + 90 * idx + 75);

                auto t_yyzz_yz = primBuffer.data(toff + 90 * idx + 76);

                auto t_yyzz_zz = primBuffer.data(toff + 90 * idx + 77);

                auto t_yzzz_xx = primBuffer.data(toff + 90 * idx + 78);

                auto t_yzzz_xy = primBuffer.data(toff + 90 * idx + 79);

                auto t_yzzz_xz = primBuffer.data(toff + 90 * idx + 80);

                auto t_yzzz_yy = primBuffer.data(toff + 90 * idx + 81);

                auto t_yzzz_yz = primBuffer.data(toff + 90 * idx + 82);

                auto t_yzzz_zz = primBuffer.data(toff + 90 * idx + 83);

                auto t_zzzz_xx = primBuffer.data(toff + 90 * idx + 84);

                auto t_zzzz_xy = primBuffer.data(toff + 90 * idx + 85);

                auto t_zzzz_xz = primBuffer.data(toff + 90 * idx + 86);

                auto t_zzzz_yy = primBuffer.data(toff + 90 * idx + 87);

                auto t_zzzz_yz = primBuffer.data(toff + 90 * idx + 88);

                auto t_zzzz_zz = primBuffer.data(toff + 90 * idx + 89);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xxx_x, tk0_xxx_y,\
                                         tk0_xxx_z, tk0_xxy_x, tk0_xxy_y, tk0_xxy_z,\
                                         tk0_xxz_x, tk0_xxz_y, tk0_xxz_z, tk0_xyy_x,\
                                         tk0_xyy_y, tk0_xyy_z, tk0_xyz_x, tk0_xyz_y,\
                                         tk0_xyz_z, tk0_xzz_x, tk0_xzz_y, tk0_xzz_z,\
                                         tk0_yyy_x, tk0_yyy_y, tk0_yyy_z, tk0_yyz_x,\
                                         tk0_yyz_y, tk0_yyz_z, tk0_yzz_x, tk0_yzz_y,\
                                         tk0_yzz_z, tk0_zzz_x, tk0_zzz_y, tk0_zzz_z,\
                                         tk1_xxx_x, tk1_xxx_y, tk1_xxx_z, tk1_xxy_x,\
                                         tk1_xxy_y, tk1_xxy_z, tk1_xxz_x, tk1_xxz_y,\
                                         tk1_xxz_z, tk1_xyy_x, tk1_xyy_y, tk1_xyy_z,\
                                         tk1_xyz_x, tk1_xyz_y, tk1_xyz_z, tk1_xzz_x,\
                                         tk1_xzz_y, tk1_xzz_z, tk1_yyy_x, tk1_yyy_y,\
                                         tk1_yyy_z, tk1_yyz_x, tk1_yyz_y, tk1_yyz_z,\
                                         tk1_yzz_x, tk1_yzz_y, tk1_yzz_z, tk1_zzz_x,\
                                         tk1_zzz_y, tk1_zzz_z, t20_xx_xx, t20_xx_xy,\
                                         t20_xx_xz, t20_xx_yy, t20_xx_yz, t20_xx_zz,\
                                         t20_xy_xx, t20_xy_xy, t20_xy_xz, t20_xy_yy,\
                                         t20_xy_yz, t20_xy_zz, t20_xz_xx, t20_xz_xy,\
                                         t20_xz_xz, t20_xz_yy, t20_xz_yz, t20_xz_zz,\
                                         t20_yy_xx, t20_yy_xy, t20_yy_xz, t20_yy_yy,\
                                         t20_yy_yz, t20_yy_zz, t20_yz_xx, t20_yz_xy,\
                                         t20_yz_xz, t20_yz_yy, t20_yz_yz, t20_yz_zz,\
                                         t20_zz_xx, t20_zz_xy, t20_zz_xz, t20_zz_yy,\
                                         t20_zz_yz, t20_zz_zz, t21_xx_xx, t21_xx_xy,\
                                         t21_xx_xz, t21_xx_yy, t21_xx_yz, t21_xx_zz,\
                                         t21_xy_xx, t21_xy_xy, t21_xy_xz, t21_xy_yy,\
                                         t21_xy_yz, t21_xy_zz, t21_xz_xx, t21_xz_xy,\
                                         t21_xz_xz, t21_xz_yy, t21_xz_yz, t21_xz_zz,\
                                         t21_yy_xx, t21_yy_xy, t21_yy_xz, t21_yy_yy,\
                                         t21_yy_yz, t21_yy_zz, t21_yz_xx, t21_yz_xy,\
                                         t21_yz_xz, t21_yz_yy, t21_yz_yz, t21_yz_zz,\
                                         t21_zz_xx, t21_zz_xy, t21_zz_xz, t21_zz_yy,\
                                         t21_zz_yz, t21_zz_zz, t10_xxx_xx, t10_xxx_xy,\
                                         t10_xxx_xz, t10_xxx_yy, t10_xxx_yz, t10_xxx_zz,\
                                         t10_xxy_xx, t10_xxy_xy, t10_xxy_xz, t10_xxy_yy,\
                                         t10_xxy_yz, t10_xxy_zz, t10_xxz_xx, t10_xxz_xy,\
                                         t10_xxz_xz, t10_xxz_yy, t10_xxz_yz, t10_xxz_zz,\
                                         t10_xyy_xx, t10_xyy_xy, t10_xyy_xz, t10_xyy_yy,\
                                         t10_xyy_yz, t10_xyy_zz, t10_xyz_xx, t10_xyz_xy,\
                                         t10_xyz_xz, t10_xyz_yy, t10_xyz_yz, t10_xyz_zz,\
                                         t10_xzz_xx, t10_xzz_xy, t10_xzz_xz, t10_xzz_yy,\
                                         t10_xzz_yz, t10_xzz_zz, t10_yyy_xx, t10_yyy_xy,\
                                         t10_yyy_xz, t10_yyy_yy, t10_yyy_yz, t10_yyy_zz,\
                                         t10_yyz_xx, t10_yyz_xy, t10_yyz_xz, t10_yyz_yy,\
                                         t10_yyz_yz, t10_yyz_zz, t10_yzz_xx, t10_yzz_xy,\
                                         t10_yzz_xz, t10_yzz_yy, t10_yzz_yz, t10_yzz_zz,\
                                         t10_zzz_xx, t10_zzz_xy, t10_zzz_xz, t10_zzz_yy,\
                                         t10_zzz_yz, t10_zzz_zz, t11_xxx_xx, t11_xxx_xy,\
                                         t11_xxx_xz, t11_xxx_yy, t11_xxx_yz, t11_xxx_zz,\
                                         t11_xxy_xx, t11_xxy_xy, t11_xxy_xz, t11_xxy_yy,\
                                         t11_xxy_yz, t11_xxy_zz, t11_xxz_xx, t11_xxz_xy,\
                                         t11_xxz_xz, t11_xxz_yy, t11_xxz_yz, t11_xxz_zz,\
                                         t11_xyy_xx, t11_xyy_xy, t11_xyy_xz, t11_xyy_yy,\
                                         t11_xyy_yz, t11_xyy_zz, t11_xyz_xx, t11_xyz_xy,\
                                         t11_xyz_xz, t11_xyz_yy, t11_xyz_yz, t11_xyz_zz,\
                                         t11_xzz_xx, t11_xzz_xy, t11_xzz_xz, t11_xzz_yy,\
                                         t11_xzz_yz, t11_xzz_zz, t11_yyy_xx, t11_yyy_xy,\
                                         t11_yyy_xz, t11_yyy_yy, t11_yyy_yz, t11_yyy_zz,\
                                         t11_yyz_xx, t11_yyz_xy, t11_yyz_xz, t11_yyz_yy,\
                                         t11_yyz_yz, t11_yyz_zz, t11_yzz_xx, t11_yzz_xy,\
                                         t11_yzz_xz, t11_yzz_yy, t11_yzz_yz, t11_yzz_zz,\
                                         t11_zzz_xx, t11_zzz_xy, t11_zzz_xz, t11_zzz_yy,\
                                         t11_zzz_yz, t11_zzz_zz, t_xxxx_xx, t_xxxx_xy,\
                                         t_xxxx_xz, t_xxxx_yy, t_xxxx_yz, t_xxxx_zz,\
                                         t_xxxy_xx, t_xxxy_xy, t_xxxy_xz, t_xxxy_yy,\
                                         t_xxxy_yz, t_xxxy_zz, t_xxxz_xx, t_xxxz_xy,\
                                         t_xxxz_xz, t_xxxz_yy, t_xxxz_yz, t_xxxz_zz,\
                                         t_xxyy_xx, t_xxyy_xy, t_xxyy_xz, t_xxyy_yy,\
                                         t_xxyy_yz, t_xxyy_zz, t_xxyz_xx, t_xxyz_xy,\
                                         t_xxyz_xz, t_xxyz_yy, t_xxyz_yz, t_xxyz_zz,\
                                         t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_yy,\
                                         t_xxzz_yz, t_xxzz_zz, t_xyyy_xx, t_xyyy_xy,\
                                         t_xyyy_xz, t_xyyy_yy, t_xyyy_yz, t_xyyy_zz,\
                                         t_xyyz_xx, t_xyyz_xy, t_xyyz_xz, t_xyyz_yy,\
                                         t_xyyz_yz, t_xyyz_zz, t_xyzz_xx, t_xyzz_xy,\
                                         t_xyzz_xz, t_xyzz_yy, t_xyzz_yz, t_xyzz_zz,\
                                         t_xzzz_xx, t_xzzz_xy, t_xzzz_xz, t_xzzz_yy,\
                                         t_xzzz_yz, t_xzzz_zz, t_yyyy_xx, t_yyyy_xy,\
                                         t_yyyy_xz, t_yyyy_yy, t_yyyy_yz, t_yyyy_zz,\
                                         t_yyyz_xx, t_yyyz_xy, t_yyyz_xz, t_yyyz_yy,\
                                         t_yyyz_yz, t_yyyz_zz, t_yyzz_xx, t_yyzz_xy,\
                                         t_yyzz_xz, t_yyzz_yy, t_yyzz_yz, t_yyzz_zz,\
                                         t_yzzz_xx, t_yzzz_xy, t_yzzz_xz, t_yzzz_yy,\
                                         t_yzzz_yz, t_yzzz_zz, t_zzzz_xx, t_zzzz_xy,\
                                         t_zzzz_xz, t_zzzz_yy, t_zzzz_yz, t_zzzz_zz,\
                                         pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxxx_xx[k] = fra * t10_xxx_xx[k] - frc * t11_xxx_xx[k] + f2t * (3.0 * t20_xx_xx[k] - 3.0 * t21_xx_xx[k] + 2.0 * tk0_xxx_x[k] - 2.0 * tk1_xxx_x[k]);

                    t_xxxx_xy[k] = fra * t10_xxx_xy[k] - frc * t11_xxx_xy[k] + f2t * (3.0 * t20_xx_xy[k] - 3.0 * t21_xx_xy[k] + tk0_xxx_y[k] - tk1_xxx_y[k]);

                    t_xxxx_xz[k] = fra * t10_xxx_xz[k] - frc * t11_xxx_xz[k] + f2t * (3.0 * t20_xx_xz[k] - 3.0 * t21_xx_xz[k] + tk0_xxx_z[k] - tk1_xxx_z[k]);

                    t_xxxx_yy[k] = fra * t10_xxx_yy[k] - frc * t11_xxx_yy[k] + f2t * (3.0 * t20_xx_yy[k] - 3.0 * t21_xx_yy[k]);

                    t_xxxx_yz[k] = fra * t10_xxx_yz[k] - frc * t11_xxx_yz[k] + f2t * (3.0 * t20_xx_yz[k] - 3.0 * t21_xx_yz[k]);

                    t_xxxx_zz[k] = fra * t10_xxx_zz[k] - frc * t11_xxx_zz[k] + f2t * (3.0 * t20_xx_zz[k] - 3.0 * t21_xx_zz[k]);

                    t_xxxy_xx[k] = fra * t10_xxy_xx[k] - frc * t11_xxy_xx[k] + f2t * (2.0 * t20_xy_xx[k] - 2.0 * t21_xy_xx[k] + 2.0 * tk0_xxy_x[k] - 2.0 * tk1_xxy_x[k]);

                    t_xxxy_xy[k] = fra * t10_xxy_xy[k] - frc * t11_xxy_xy[k] + f2t * (2.0 * t20_xy_xy[k] - 2.0 * t21_xy_xy[k] + tk0_xxy_y[k] - tk1_xxy_y[k]);

                    t_xxxy_xz[k] = fra * t10_xxy_xz[k] - frc * t11_xxy_xz[k] + f2t * (2.0 * t20_xy_xz[k] - 2.0 * t21_xy_xz[k] + tk0_xxy_z[k] - tk1_xxy_z[k]);

                    t_xxxy_yy[k] = fra * t10_xxy_yy[k] - frc * t11_xxy_yy[k] + f2t * (2.0 * t20_xy_yy[k] - 2.0 * t21_xy_yy[k]);

                    t_xxxy_yz[k] = fra * t10_xxy_yz[k] - frc * t11_xxy_yz[k] + f2t * (2.0 * t20_xy_yz[k] - 2.0 * t21_xy_yz[k]);

                    t_xxxy_zz[k] = fra * t10_xxy_zz[k] - frc * t11_xxy_zz[k] + f2t * (2.0 * t20_xy_zz[k] - 2.0 * t21_xy_zz[k]);

                    t_xxxz_xx[k] = fra * t10_xxz_xx[k] - frc * t11_xxz_xx[k] + f2t * (2.0 * t20_xz_xx[k] - 2.0 * t21_xz_xx[k] + 2.0 * tk0_xxz_x[k] - 2.0 * tk1_xxz_x[k]);

                    t_xxxz_xy[k] = fra * t10_xxz_xy[k] - frc * t11_xxz_xy[k] + f2t * (2.0 * t20_xz_xy[k] - 2.0 * t21_xz_xy[k] + tk0_xxz_y[k] - tk1_xxz_y[k]);

                    t_xxxz_xz[k] = fra * t10_xxz_xz[k] - frc * t11_xxz_xz[k] + f2t * (2.0 * t20_xz_xz[k] - 2.0 * t21_xz_xz[k] + tk0_xxz_z[k] - tk1_xxz_z[k]);

                    t_xxxz_yy[k] = fra * t10_xxz_yy[k] - frc * t11_xxz_yy[k] + f2t * (2.0 * t20_xz_yy[k] - 2.0 * t21_xz_yy[k]);

                    t_xxxz_yz[k] = fra * t10_xxz_yz[k] - frc * t11_xxz_yz[k] + f2t * (2.0 * t20_xz_yz[k] - 2.0 * t21_xz_yz[k]);

                    t_xxxz_zz[k] = fra * t10_xxz_zz[k] - frc * t11_xxz_zz[k] + f2t * (2.0 * t20_xz_zz[k] - 2.0 * t21_xz_zz[k]);

                    t_xxyy_xx[k] = fra * t10_xyy_xx[k] - frc * t11_xyy_xx[k] + f2t * (t20_yy_xx[k] - t21_yy_xx[k] + 2.0 * tk0_xyy_x[k] - 2.0 * tk1_xyy_x[k]);

                    t_xxyy_xy[k] = fra * t10_xyy_xy[k] - frc * t11_xyy_xy[k] + f2t * (t20_yy_xy[k] - t21_yy_xy[k] + tk0_xyy_y[k] - tk1_xyy_y[k]);

                    t_xxyy_xz[k] = fra * t10_xyy_xz[k] - frc * t11_xyy_xz[k] + f2t * (t20_yy_xz[k] - t21_yy_xz[k] + tk0_xyy_z[k] - tk1_xyy_z[k]);

                    t_xxyy_yy[k] = fra * t10_xyy_yy[k] - frc * t11_xyy_yy[k] + f2t * (t20_yy_yy[k] - t21_yy_yy[k]);

                    t_xxyy_yz[k] = fra * t10_xyy_yz[k] - frc * t11_xyy_yz[k] + f2t * (t20_yy_yz[k] - t21_yy_yz[k]);

                    t_xxyy_zz[k] = fra * t10_xyy_zz[k] - frc * t11_xyy_zz[k] + f2t * (t20_yy_zz[k] - t21_yy_zz[k]);

                    t_xxyz_xx[k] = fra * t10_xyz_xx[k] - frc * t11_xyz_xx[k] + f2t * (t20_yz_xx[k] - t21_yz_xx[k] + 2.0 * tk0_xyz_x[k] - 2.0 * tk1_xyz_x[k]);

                    t_xxyz_xy[k] = fra * t10_xyz_xy[k] - frc * t11_xyz_xy[k] + f2t * (t20_yz_xy[k] - t21_yz_xy[k] + tk0_xyz_y[k] - tk1_xyz_y[k]);

                    t_xxyz_xz[k] = fra * t10_xyz_xz[k] - frc * t11_xyz_xz[k] + f2t * (t20_yz_xz[k] - t21_yz_xz[k] + tk0_xyz_z[k] - tk1_xyz_z[k]);

                    t_xxyz_yy[k] = fra * t10_xyz_yy[k] - frc * t11_xyz_yy[k] + f2t * (t20_yz_yy[k] - t21_yz_yy[k]);

                    t_xxyz_yz[k] = fra * t10_xyz_yz[k] - frc * t11_xyz_yz[k] + f2t * (t20_yz_yz[k] - t21_yz_yz[k]);

                    t_xxyz_zz[k] = fra * t10_xyz_zz[k] - frc * t11_xyz_zz[k] + f2t * (t20_yz_zz[k] - t21_yz_zz[k]);

                    t_xxzz_xx[k] = fra * t10_xzz_xx[k] - frc * t11_xzz_xx[k] + f2t * (t20_zz_xx[k] - t21_zz_xx[k] + 2.0 * tk0_xzz_x[k] - 2.0 * tk1_xzz_x[k]);

                    t_xxzz_xy[k] = fra * t10_xzz_xy[k] - frc * t11_xzz_xy[k] + f2t * (t20_zz_xy[k] - t21_zz_xy[k] + tk0_xzz_y[k] - tk1_xzz_y[k]);

                    t_xxzz_xz[k] = fra * t10_xzz_xz[k] - frc * t11_xzz_xz[k] + f2t * (t20_zz_xz[k] - t21_zz_xz[k] + tk0_xzz_z[k] - tk1_xzz_z[k]);

                    t_xxzz_yy[k] = fra * t10_xzz_yy[k] - frc * t11_xzz_yy[k] + f2t * (t20_zz_yy[k] - t21_zz_yy[k]);

                    t_xxzz_yz[k] = fra * t10_xzz_yz[k] - frc * t11_xzz_yz[k] + f2t * (t20_zz_yz[k] - t21_zz_yz[k]);

                    t_xxzz_zz[k] = fra * t10_xzz_zz[k] - frc * t11_xzz_zz[k] + f2t * (t20_zz_zz[k] - t21_zz_zz[k]);

                    t_xyyy_xx[k] = fra * t10_yyy_xx[k] - frc * t11_yyy_xx[k] + f2t * (2.0 * tk0_yyy_x[k] - 2.0 * tk1_yyy_x[k]);

                    t_xyyy_xy[k] = fra * t10_yyy_xy[k] - frc * t11_yyy_xy[k] + f2t * (tk0_yyy_y[k] - tk1_yyy_y[k]);

                    t_xyyy_xz[k] = fra * t10_yyy_xz[k] - frc * t11_yyy_xz[k] + f2t * (tk0_yyy_z[k] - tk1_yyy_z[k]);

                    t_xyyy_yy[k] = fra * t10_yyy_yy[k] - frc * t11_yyy_yy[k];

                    t_xyyy_yz[k] = fra * t10_yyy_yz[k] - frc * t11_yyy_yz[k];

                    t_xyyy_zz[k] = fra * t10_yyy_zz[k] - frc * t11_yyy_zz[k];

                    t_xyyz_xx[k] = fra * t10_yyz_xx[k] - frc * t11_yyz_xx[k] + f2t * (2.0 * tk0_yyz_x[k] - 2.0 * tk1_yyz_x[k]);

                    t_xyyz_xy[k] = fra * t10_yyz_xy[k] - frc * t11_yyz_xy[k] + f2t * (tk0_yyz_y[k] - tk1_yyz_y[k]);

                    t_xyyz_xz[k] = fra * t10_yyz_xz[k] - frc * t11_yyz_xz[k] + f2t * (tk0_yyz_z[k] - tk1_yyz_z[k]);

                    t_xyyz_yy[k] = fra * t10_yyz_yy[k] - frc * t11_yyz_yy[k];

                    t_xyyz_yz[k] = fra * t10_yyz_yz[k] - frc * t11_yyz_yz[k];

                    t_xyyz_zz[k] = fra * t10_yyz_zz[k] - frc * t11_yyz_zz[k];

                    t_xyzz_xx[k] = fra * t10_yzz_xx[k] - frc * t11_yzz_xx[k] + f2t * (2.0 * tk0_yzz_x[k] - 2.0 * tk1_yzz_x[k]);

                    t_xyzz_xy[k] = fra * t10_yzz_xy[k] - frc * t11_yzz_xy[k] + f2t * (tk0_yzz_y[k] - tk1_yzz_y[k]);

                    t_xyzz_xz[k] = fra * t10_yzz_xz[k] - frc * t11_yzz_xz[k] + f2t * (tk0_yzz_z[k] - tk1_yzz_z[k]);

                    t_xyzz_yy[k] = fra * t10_yzz_yy[k] - frc * t11_yzz_yy[k];

                    t_xyzz_yz[k] = fra * t10_yzz_yz[k] - frc * t11_yzz_yz[k];

                    t_xyzz_zz[k] = fra * t10_yzz_zz[k] - frc * t11_yzz_zz[k];

                    t_xzzz_xx[k] = fra * t10_zzz_xx[k] - frc * t11_zzz_xx[k] + f2t * (2.0 * tk0_zzz_x[k] - 2.0 * tk1_zzz_x[k]);

                    t_xzzz_xy[k] = fra * t10_zzz_xy[k] - frc * t11_zzz_xy[k] + f2t * (tk0_zzz_y[k] - tk1_zzz_y[k]);

                    t_xzzz_xz[k] = fra * t10_zzz_xz[k] - frc * t11_zzz_xz[k] + f2t * (tk0_zzz_z[k] - tk1_zzz_z[k]);

                    t_xzzz_yy[k] = fra * t10_zzz_yy[k] - frc * t11_zzz_yy[k];

                    t_xzzz_yz[k] = fra * t10_zzz_yz[k] - frc * t11_zzz_yz[k];

                    t_xzzz_zz[k] = fra * t10_zzz_zz[k] - frc * t11_zzz_zz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyyy_xx[k] = fra * t10_yyy_xx[k] - frc * t11_yyy_xx[k] + f2t * (3.0 * t20_yy_xx[k] - 3.0 * t21_yy_xx[k]);

                    t_yyyy_xy[k] = fra * t10_yyy_xy[k] - frc * t11_yyy_xy[k] + f2t * (3.0 * t20_yy_xy[k] - 3.0 * t21_yy_xy[k] + tk0_yyy_x[k] - tk1_yyy_x[k]);

                    t_yyyy_xz[k] = fra * t10_yyy_xz[k] - frc * t11_yyy_xz[k] + f2t * (3.0 * t20_yy_xz[k] - 3.0 * t21_yy_xz[k]);

                    t_yyyy_yy[k] = fra * t10_yyy_yy[k] - frc * t11_yyy_yy[k] + f2t * (3.0 * t20_yy_yy[k] - 3.0 * t21_yy_yy[k] + 2.0 * tk0_yyy_y[k] - 2.0 * tk1_yyy_y[k]);

                    t_yyyy_yz[k] = fra * t10_yyy_yz[k] - frc * t11_yyy_yz[k] + f2t * (3.0 * t20_yy_yz[k] - 3.0 * t21_yy_yz[k] + tk0_yyy_z[k] - tk1_yyy_z[k]);

                    t_yyyy_zz[k] = fra * t10_yyy_zz[k] - frc * t11_yyy_zz[k] + f2t * (3.0 * t20_yy_zz[k] - 3.0 * t21_yy_zz[k]);

                    t_yyyz_xx[k] = fra * t10_yyz_xx[k] - frc * t11_yyz_xx[k] + f2t * (2.0 * t20_yz_xx[k] - 2.0 * t21_yz_xx[k]);

                    t_yyyz_xy[k] = fra * t10_yyz_xy[k] - frc * t11_yyz_xy[k] + f2t * (2.0 * t20_yz_xy[k] - 2.0 * t21_yz_xy[k] + tk0_yyz_x[k] - tk1_yyz_x[k]);

                    t_yyyz_xz[k] = fra * t10_yyz_xz[k] - frc * t11_yyz_xz[k] + f2t * (2.0 * t20_yz_xz[k] - 2.0 * t21_yz_xz[k]);

                    t_yyyz_yy[k] = fra * t10_yyz_yy[k] - frc * t11_yyz_yy[k] + f2t * (2.0 * t20_yz_yy[k] - 2.0 * t21_yz_yy[k] + 2.0 * tk0_yyz_y[k] - 2.0 * tk1_yyz_y[k]);

                    t_yyyz_yz[k] = fra * t10_yyz_yz[k] - frc * t11_yyz_yz[k] + f2t * (2.0 * t20_yz_yz[k] - 2.0 * t21_yz_yz[k] + tk0_yyz_z[k] - tk1_yyz_z[k]);

                    t_yyyz_zz[k] = fra * t10_yyz_zz[k] - frc * t11_yyz_zz[k] + f2t * (2.0 * t20_yz_zz[k] - 2.0 * t21_yz_zz[k]);

                    t_yyzz_xx[k] = fra * t10_yzz_xx[k] - frc * t11_yzz_xx[k] + f2t * (t20_zz_xx[k] - t21_zz_xx[k]);

                    t_yyzz_xy[k] = fra * t10_yzz_xy[k] - frc * t11_yzz_xy[k] + f2t * (t20_zz_xy[k] - t21_zz_xy[k] + tk0_yzz_x[k] - tk1_yzz_x[k]);

                    t_yyzz_xz[k] = fra * t10_yzz_xz[k] - frc * t11_yzz_xz[k] + f2t * (t20_zz_xz[k] - t21_zz_xz[k]);

                    t_yyzz_yy[k] = fra * t10_yzz_yy[k] - frc * t11_yzz_yy[k] + f2t * (t20_zz_yy[k] - t21_zz_yy[k] + 2.0 * tk0_yzz_y[k] - 2.0 * tk1_yzz_y[k]);

                    t_yyzz_yz[k] = fra * t10_yzz_yz[k] - frc * t11_yzz_yz[k] + f2t * (t20_zz_yz[k] - t21_zz_yz[k] + tk0_yzz_z[k] - tk1_yzz_z[k]);

                    t_yyzz_zz[k] = fra * t10_yzz_zz[k] - frc * t11_yzz_zz[k] + f2t * (t20_zz_zz[k] - t21_zz_zz[k]);

                    t_yzzz_xx[k] = fra * t10_zzz_xx[k] - frc * t11_zzz_xx[k];

                    t_yzzz_xy[k] = fra * t10_zzz_xy[k] - frc * t11_zzz_xy[k] + f2t * (tk0_zzz_x[k] - tk1_zzz_x[k]);

                    t_yzzz_xz[k] = fra * t10_zzz_xz[k] - frc * t11_zzz_xz[k];

                    t_yzzz_yy[k] = fra * t10_zzz_yy[k] - frc * t11_zzz_yy[k] + f2t * (2.0 * tk0_zzz_y[k] - 2.0 * tk1_zzz_y[k]);

                    t_yzzz_yz[k] = fra * t10_zzz_yz[k] - frc * t11_zzz_yz[k] + f2t * (tk0_zzz_z[k] - tk1_zzz_z[k]);

                    t_yzzz_zz[k] = fra * t10_zzz_zz[k] - frc * t11_zzz_zz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzzz_xx[k] = fra * t10_zzz_xx[k] - frc * t11_zzz_xx[k] + f2t * (3.0 * t20_zz_xx[k] - 3.0 * t21_zz_xx[k]);

                    t_zzzz_xy[k] = fra * t10_zzz_xy[k] - frc * t11_zzz_xy[k] + f2t * (3.0 * t20_zz_xy[k] - 3.0 * t21_zz_xy[k]);

                    t_zzzz_xz[k] = fra * t10_zzz_xz[k] - frc * t11_zzz_xz[k] + f2t * (3.0 * t20_zz_xz[k] - 3.0 * t21_zz_xz[k] + tk0_zzz_x[k] - tk1_zzz_x[k]);

                    t_zzzz_yy[k] = fra * t10_zzz_yy[k] - frc * t11_zzz_yy[k] + f2t * (3.0 * t20_zz_yy[k] - 3.0 * t21_zz_yy[k]);

                    t_zzzz_yz[k] = fra * t10_zzz_yz[k] - frc * t11_zzz_yz[k] + f2t * (3.0 * t20_zz_yz[k] - 3.0 * t21_zz_yz[k] + tk0_zzz_y[k] - tk1_zzz_y[k]);

                    t_zzzz_zz[k] = fra * t10_zzz_zz[k] - frc * t11_zzz_zz[k] + f2t * (3.0 * t20_zz_zz[k] - 3.0 * t21_zz_zz[k] + 2.0 * tk0_zzz_z[k] - 2.0 * tk1_zzz_z[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForFG(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {3, 4, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 3, 4);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 4, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 4, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 4, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 4, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {1, 4, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i + 1});

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

                // set up pointers to (D|A(0)|F)^(m) integrals

                auto tk0_xx_xxx = primBuffer.data(tk0off + 60 * idx);

                auto tk0_xx_xxy = primBuffer.data(tk0off + 60 * idx + 1);

                auto tk0_xx_xxz = primBuffer.data(tk0off + 60 * idx + 2);

                auto tk0_xx_xyy = primBuffer.data(tk0off + 60 * idx + 3);

                auto tk0_xx_xyz = primBuffer.data(tk0off + 60 * idx + 4);

                auto tk0_xx_xzz = primBuffer.data(tk0off + 60 * idx + 5);

                auto tk0_xx_yyy = primBuffer.data(tk0off + 60 * idx + 6);

                auto tk0_xx_yyz = primBuffer.data(tk0off + 60 * idx + 7);

                auto tk0_xx_yzz = primBuffer.data(tk0off + 60 * idx + 8);

                auto tk0_xx_zzz = primBuffer.data(tk0off + 60 * idx + 9);

                auto tk0_xy_xxx = primBuffer.data(tk0off + 60 * idx + 10);

                auto tk0_xy_xxy = primBuffer.data(tk0off + 60 * idx + 11);

                auto tk0_xy_xxz = primBuffer.data(tk0off + 60 * idx + 12);

                auto tk0_xy_xyy = primBuffer.data(tk0off + 60 * idx + 13);

                auto tk0_xy_xyz = primBuffer.data(tk0off + 60 * idx + 14);

                auto tk0_xy_xzz = primBuffer.data(tk0off + 60 * idx + 15);

                auto tk0_xy_yyy = primBuffer.data(tk0off + 60 * idx + 16);

                auto tk0_xy_yyz = primBuffer.data(tk0off + 60 * idx + 17);

                auto tk0_xy_yzz = primBuffer.data(tk0off + 60 * idx + 18);

                auto tk0_xy_zzz = primBuffer.data(tk0off + 60 * idx + 19);

                auto tk0_xz_xxx = primBuffer.data(tk0off + 60 * idx + 20);

                auto tk0_xz_xxy = primBuffer.data(tk0off + 60 * idx + 21);

                auto tk0_xz_xxz = primBuffer.data(tk0off + 60 * idx + 22);

                auto tk0_xz_xyy = primBuffer.data(tk0off + 60 * idx + 23);

                auto tk0_xz_xyz = primBuffer.data(tk0off + 60 * idx + 24);

                auto tk0_xz_xzz = primBuffer.data(tk0off + 60 * idx + 25);

                auto tk0_xz_yyy = primBuffer.data(tk0off + 60 * idx + 26);

                auto tk0_xz_yyz = primBuffer.data(tk0off + 60 * idx + 27);

                auto tk0_xz_yzz = primBuffer.data(tk0off + 60 * idx + 28);

                auto tk0_xz_zzz = primBuffer.data(tk0off + 60 * idx + 29);

                auto tk0_yy_xxx = primBuffer.data(tk0off + 60 * idx + 30);

                auto tk0_yy_xxy = primBuffer.data(tk0off + 60 * idx + 31);

                auto tk0_yy_xxz = primBuffer.data(tk0off + 60 * idx + 32);

                auto tk0_yy_xyy = primBuffer.data(tk0off + 60 * idx + 33);

                auto tk0_yy_xyz = primBuffer.data(tk0off + 60 * idx + 34);

                auto tk0_yy_xzz = primBuffer.data(tk0off + 60 * idx + 35);

                auto tk0_yy_yyy = primBuffer.data(tk0off + 60 * idx + 36);

                auto tk0_yy_yyz = primBuffer.data(tk0off + 60 * idx + 37);

                auto tk0_yy_yzz = primBuffer.data(tk0off + 60 * idx + 38);

                auto tk0_yy_zzz = primBuffer.data(tk0off + 60 * idx + 39);

                auto tk0_yz_xxx = primBuffer.data(tk0off + 60 * idx + 40);

                auto tk0_yz_xxy = primBuffer.data(tk0off + 60 * idx + 41);

                auto tk0_yz_xxz = primBuffer.data(tk0off + 60 * idx + 42);

                auto tk0_yz_xyy = primBuffer.data(tk0off + 60 * idx + 43);

                auto tk0_yz_xyz = primBuffer.data(tk0off + 60 * idx + 44);

                auto tk0_yz_xzz = primBuffer.data(tk0off + 60 * idx + 45);

                auto tk0_yz_yyy = primBuffer.data(tk0off + 60 * idx + 46);

                auto tk0_yz_yyz = primBuffer.data(tk0off + 60 * idx + 47);

                auto tk0_yz_yzz = primBuffer.data(tk0off + 60 * idx + 48);

                auto tk0_yz_zzz = primBuffer.data(tk0off + 60 * idx + 49);

                auto tk0_zz_xxx = primBuffer.data(tk0off + 60 * idx + 50);

                auto tk0_zz_xxy = primBuffer.data(tk0off + 60 * idx + 51);

                auto tk0_zz_xxz = primBuffer.data(tk0off + 60 * idx + 52);

                auto tk0_zz_xyy = primBuffer.data(tk0off + 60 * idx + 53);

                auto tk0_zz_xyz = primBuffer.data(tk0off + 60 * idx + 54);

                auto tk0_zz_xzz = primBuffer.data(tk0off + 60 * idx + 55);

                auto tk0_zz_yyy = primBuffer.data(tk0off + 60 * idx + 56);

                auto tk0_zz_yyz = primBuffer.data(tk0off + 60 * idx + 57);

                auto tk0_zz_yzz = primBuffer.data(tk0off + 60 * idx + 58);

                auto tk0_zz_zzz = primBuffer.data(tk0off + 60 * idx + 59);

                // set up pointers to (D|A(0)|F)^(m+1) integrals

                auto tk1_xx_xxx = primBuffer.data(tk1off + 60 * idx);

                auto tk1_xx_xxy = primBuffer.data(tk1off + 60 * idx + 1);

                auto tk1_xx_xxz = primBuffer.data(tk1off + 60 * idx + 2);

                auto tk1_xx_xyy = primBuffer.data(tk1off + 60 * idx + 3);

                auto tk1_xx_xyz = primBuffer.data(tk1off + 60 * idx + 4);

                auto tk1_xx_xzz = primBuffer.data(tk1off + 60 * idx + 5);

                auto tk1_xx_yyy = primBuffer.data(tk1off + 60 * idx + 6);

                auto tk1_xx_yyz = primBuffer.data(tk1off + 60 * idx + 7);

                auto tk1_xx_yzz = primBuffer.data(tk1off + 60 * idx + 8);

                auto tk1_xx_zzz = primBuffer.data(tk1off + 60 * idx + 9);

                auto tk1_xy_xxx = primBuffer.data(tk1off + 60 * idx + 10);

                auto tk1_xy_xxy = primBuffer.data(tk1off + 60 * idx + 11);

                auto tk1_xy_xxz = primBuffer.data(tk1off + 60 * idx + 12);

                auto tk1_xy_xyy = primBuffer.data(tk1off + 60 * idx + 13);

                auto tk1_xy_xyz = primBuffer.data(tk1off + 60 * idx + 14);

                auto tk1_xy_xzz = primBuffer.data(tk1off + 60 * idx + 15);

                auto tk1_xy_yyy = primBuffer.data(tk1off + 60 * idx + 16);

                auto tk1_xy_yyz = primBuffer.data(tk1off + 60 * idx + 17);

                auto tk1_xy_yzz = primBuffer.data(tk1off + 60 * idx + 18);

                auto tk1_xy_zzz = primBuffer.data(tk1off + 60 * idx + 19);

                auto tk1_xz_xxx = primBuffer.data(tk1off + 60 * idx + 20);

                auto tk1_xz_xxy = primBuffer.data(tk1off + 60 * idx + 21);

                auto tk1_xz_xxz = primBuffer.data(tk1off + 60 * idx + 22);

                auto tk1_xz_xyy = primBuffer.data(tk1off + 60 * idx + 23);

                auto tk1_xz_xyz = primBuffer.data(tk1off + 60 * idx + 24);

                auto tk1_xz_xzz = primBuffer.data(tk1off + 60 * idx + 25);

                auto tk1_xz_yyy = primBuffer.data(tk1off + 60 * idx + 26);

                auto tk1_xz_yyz = primBuffer.data(tk1off + 60 * idx + 27);

                auto tk1_xz_yzz = primBuffer.data(tk1off + 60 * idx + 28);

                auto tk1_xz_zzz = primBuffer.data(tk1off + 60 * idx + 29);

                auto tk1_yy_xxx = primBuffer.data(tk1off + 60 * idx + 30);

                auto tk1_yy_xxy = primBuffer.data(tk1off + 60 * idx + 31);

                auto tk1_yy_xxz = primBuffer.data(tk1off + 60 * idx + 32);

                auto tk1_yy_xyy = primBuffer.data(tk1off + 60 * idx + 33);

                auto tk1_yy_xyz = primBuffer.data(tk1off + 60 * idx + 34);

                auto tk1_yy_xzz = primBuffer.data(tk1off + 60 * idx + 35);

                auto tk1_yy_yyy = primBuffer.data(tk1off + 60 * idx + 36);

                auto tk1_yy_yyz = primBuffer.data(tk1off + 60 * idx + 37);

                auto tk1_yy_yzz = primBuffer.data(tk1off + 60 * idx + 38);

                auto tk1_yy_zzz = primBuffer.data(tk1off + 60 * idx + 39);

                auto tk1_yz_xxx = primBuffer.data(tk1off + 60 * idx + 40);

                auto tk1_yz_xxy = primBuffer.data(tk1off + 60 * idx + 41);

                auto tk1_yz_xxz = primBuffer.data(tk1off + 60 * idx + 42);

                auto tk1_yz_xyy = primBuffer.data(tk1off + 60 * idx + 43);

                auto tk1_yz_xyz = primBuffer.data(tk1off + 60 * idx + 44);

                auto tk1_yz_xzz = primBuffer.data(tk1off + 60 * idx + 45);

                auto tk1_yz_yyy = primBuffer.data(tk1off + 60 * idx + 46);

                auto tk1_yz_yyz = primBuffer.data(tk1off + 60 * idx + 47);

                auto tk1_yz_yzz = primBuffer.data(tk1off + 60 * idx + 48);

                auto tk1_yz_zzz = primBuffer.data(tk1off + 60 * idx + 49);

                auto tk1_zz_xxx = primBuffer.data(tk1off + 60 * idx + 50);

                auto tk1_zz_xxy = primBuffer.data(tk1off + 60 * idx + 51);

                auto tk1_zz_xxz = primBuffer.data(tk1off + 60 * idx + 52);

                auto tk1_zz_xyy = primBuffer.data(tk1off + 60 * idx + 53);

                auto tk1_zz_xyz = primBuffer.data(tk1off + 60 * idx + 54);

                auto tk1_zz_xzz = primBuffer.data(tk1off + 60 * idx + 55);

                auto tk1_zz_yyy = primBuffer.data(tk1off + 60 * idx + 56);

                auto tk1_zz_yyz = primBuffer.data(tk1off + 60 * idx + 57);

                auto tk1_zz_yzz = primBuffer.data(tk1off + 60 * idx + 58);

                auto tk1_zz_zzz = primBuffer.data(tk1off + 60 * idx + 59);

                // set up pointers to (P|A(0)|G)^(m) integrals

                auto t20_x_xxxx = primBuffer.data(t20off + 45 * idx);

                auto t20_x_xxxy = primBuffer.data(t20off + 45 * idx + 1);

                auto t20_x_xxxz = primBuffer.data(t20off + 45 * idx + 2);

                auto t20_x_xxyy = primBuffer.data(t20off + 45 * idx + 3);

                auto t20_x_xxyz = primBuffer.data(t20off + 45 * idx + 4);

                auto t20_x_xxzz = primBuffer.data(t20off + 45 * idx + 5);

                auto t20_x_xyyy = primBuffer.data(t20off + 45 * idx + 6);

                auto t20_x_xyyz = primBuffer.data(t20off + 45 * idx + 7);

                auto t20_x_xyzz = primBuffer.data(t20off + 45 * idx + 8);

                auto t20_x_xzzz = primBuffer.data(t20off + 45 * idx + 9);

                auto t20_x_yyyy = primBuffer.data(t20off + 45 * idx + 10);

                auto t20_x_yyyz = primBuffer.data(t20off + 45 * idx + 11);

                auto t20_x_yyzz = primBuffer.data(t20off + 45 * idx + 12);

                auto t20_x_yzzz = primBuffer.data(t20off + 45 * idx + 13);

                auto t20_x_zzzz = primBuffer.data(t20off + 45 * idx + 14);

                auto t20_y_xxxx = primBuffer.data(t20off + 45 * idx + 15);

                auto t20_y_xxxy = primBuffer.data(t20off + 45 * idx + 16);

                auto t20_y_xxxz = primBuffer.data(t20off + 45 * idx + 17);

                auto t20_y_xxyy = primBuffer.data(t20off + 45 * idx + 18);

                auto t20_y_xxyz = primBuffer.data(t20off + 45 * idx + 19);

                auto t20_y_xxzz = primBuffer.data(t20off + 45 * idx + 20);

                auto t20_y_xyyy = primBuffer.data(t20off + 45 * idx + 21);

                auto t20_y_xyyz = primBuffer.data(t20off + 45 * idx + 22);

                auto t20_y_xyzz = primBuffer.data(t20off + 45 * idx + 23);

                auto t20_y_xzzz = primBuffer.data(t20off + 45 * idx + 24);

                auto t20_y_yyyy = primBuffer.data(t20off + 45 * idx + 25);

                auto t20_y_yyyz = primBuffer.data(t20off + 45 * idx + 26);

                auto t20_y_yyzz = primBuffer.data(t20off + 45 * idx + 27);

                auto t20_y_yzzz = primBuffer.data(t20off + 45 * idx + 28);

                auto t20_y_zzzz = primBuffer.data(t20off + 45 * idx + 29);

                auto t20_z_xxxx = primBuffer.data(t20off + 45 * idx + 30);

                auto t20_z_xxxy = primBuffer.data(t20off + 45 * idx + 31);

                auto t20_z_xxxz = primBuffer.data(t20off + 45 * idx + 32);

                auto t20_z_xxyy = primBuffer.data(t20off + 45 * idx + 33);

                auto t20_z_xxyz = primBuffer.data(t20off + 45 * idx + 34);

                auto t20_z_xxzz = primBuffer.data(t20off + 45 * idx + 35);

                auto t20_z_xyyy = primBuffer.data(t20off + 45 * idx + 36);

                auto t20_z_xyyz = primBuffer.data(t20off + 45 * idx + 37);

                auto t20_z_xyzz = primBuffer.data(t20off + 45 * idx + 38);

                auto t20_z_xzzz = primBuffer.data(t20off + 45 * idx + 39);

                auto t20_z_yyyy = primBuffer.data(t20off + 45 * idx + 40);

                auto t20_z_yyyz = primBuffer.data(t20off + 45 * idx + 41);

                auto t20_z_yyzz = primBuffer.data(t20off + 45 * idx + 42);

                auto t20_z_yzzz = primBuffer.data(t20off + 45 * idx + 43);

                auto t20_z_zzzz = primBuffer.data(t20off + 45 * idx + 44);

                // set up pointers to (P|A(0)|G)^(m+1) integrals

                auto t21_x_xxxx = primBuffer.data(t21off + 45 * idx);

                auto t21_x_xxxy = primBuffer.data(t21off + 45 * idx + 1);

                auto t21_x_xxxz = primBuffer.data(t21off + 45 * idx + 2);

                auto t21_x_xxyy = primBuffer.data(t21off + 45 * idx + 3);

                auto t21_x_xxyz = primBuffer.data(t21off + 45 * idx + 4);

                auto t21_x_xxzz = primBuffer.data(t21off + 45 * idx + 5);

                auto t21_x_xyyy = primBuffer.data(t21off + 45 * idx + 6);

                auto t21_x_xyyz = primBuffer.data(t21off + 45 * idx + 7);

                auto t21_x_xyzz = primBuffer.data(t21off + 45 * idx + 8);

                auto t21_x_xzzz = primBuffer.data(t21off + 45 * idx + 9);

                auto t21_x_yyyy = primBuffer.data(t21off + 45 * idx + 10);

                auto t21_x_yyyz = primBuffer.data(t21off + 45 * idx + 11);

                auto t21_x_yyzz = primBuffer.data(t21off + 45 * idx + 12);

                auto t21_x_yzzz = primBuffer.data(t21off + 45 * idx + 13);

                auto t21_x_zzzz = primBuffer.data(t21off + 45 * idx + 14);

                auto t21_y_xxxx = primBuffer.data(t21off + 45 * idx + 15);

                auto t21_y_xxxy = primBuffer.data(t21off + 45 * idx + 16);

                auto t21_y_xxxz = primBuffer.data(t21off + 45 * idx + 17);

                auto t21_y_xxyy = primBuffer.data(t21off + 45 * idx + 18);

                auto t21_y_xxyz = primBuffer.data(t21off + 45 * idx + 19);

                auto t21_y_xxzz = primBuffer.data(t21off + 45 * idx + 20);

                auto t21_y_xyyy = primBuffer.data(t21off + 45 * idx + 21);

                auto t21_y_xyyz = primBuffer.data(t21off + 45 * idx + 22);

                auto t21_y_xyzz = primBuffer.data(t21off + 45 * idx + 23);

                auto t21_y_xzzz = primBuffer.data(t21off + 45 * idx + 24);

                auto t21_y_yyyy = primBuffer.data(t21off + 45 * idx + 25);

                auto t21_y_yyyz = primBuffer.data(t21off + 45 * idx + 26);

                auto t21_y_yyzz = primBuffer.data(t21off + 45 * idx + 27);

                auto t21_y_yzzz = primBuffer.data(t21off + 45 * idx + 28);

                auto t21_y_zzzz = primBuffer.data(t21off + 45 * idx + 29);

                auto t21_z_xxxx = primBuffer.data(t21off + 45 * idx + 30);

                auto t21_z_xxxy = primBuffer.data(t21off + 45 * idx + 31);

                auto t21_z_xxxz = primBuffer.data(t21off + 45 * idx + 32);

                auto t21_z_xxyy = primBuffer.data(t21off + 45 * idx + 33);

                auto t21_z_xxyz = primBuffer.data(t21off + 45 * idx + 34);

                auto t21_z_xxzz = primBuffer.data(t21off + 45 * idx + 35);

                auto t21_z_xyyy = primBuffer.data(t21off + 45 * idx + 36);

                auto t21_z_xyyz = primBuffer.data(t21off + 45 * idx + 37);

                auto t21_z_xyzz = primBuffer.data(t21off + 45 * idx + 38);

                auto t21_z_xzzz = primBuffer.data(t21off + 45 * idx + 39);

                auto t21_z_yyyy = primBuffer.data(t21off + 45 * idx + 40);

                auto t21_z_yyyz = primBuffer.data(t21off + 45 * idx + 41);

                auto t21_z_yyzz = primBuffer.data(t21off + 45 * idx + 42);

                auto t21_z_yzzz = primBuffer.data(t21off + 45 * idx + 43);

                auto t21_z_zzzz = primBuffer.data(t21off + 45 * idx + 44);

                // set up pointers to (D|A(0)|G)^(m) integrals

                auto t10_xx_xxxx = primBuffer.data(t10off + 90 * idx);

                auto t10_xx_xxxy = primBuffer.data(t10off + 90 * idx + 1);

                auto t10_xx_xxxz = primBuffer.data(t10off + 90 * idx + 2);

                auto t10_xx_xxyy = primBuffer.data(t10off + 90 * idx + 3);

                auto t10_xx_xxyz = primBuffer.data(t10off + 90 * idx + 4);

                auto t10_xx_xxzz = primBuffer.data(t10off + 90 * idx + 5);

                auto t10_xx_xyyy = primBuffer.data(t10off + 90 * idx + 6);

                auto t10_xx_xyyz = primBuffer.data(t10off + 90 * idx + 7);

                auto t10_xx_xyzz = primBuffer.data(t10off + 90 * idx + 8);

                auto t10_xx_xzzz = primBuffer.data(t10off + 90 * idx + 9);

                auto t10_xx_yyyy = primBuffer.data(t10off + 90 * idx + 10);

                auto t10_xx_yyyz = primBuffer.data(t10off + 90 * idx + 11);

                auto t10_xx_yyzz = primBuffer.data(t10off + 90 * idx + 12);

                auto t10_xx_yzzz = primBuffer.data(t10off + 90 * idx + 13);

                auto t10_xx_zzzz = primBuffer.data(t10off + 90 * idx + 14);

                auto t10_xy_xxxx = primBuffer.data(t10off + 90 * idx + 15);

                auto t10_xy_xxxy = primBuffer.data(t10off + 90 * idx + 16);

                auto t10_xy_xxxz = primBuffer.data(t10off + 90 * idx + 17);

                auto t10_xy_xxyy = primBuffer.data(t10off + 90 * idx + 18);

                auto t10_xy_xxyz = primBuffer.data(t10off + 90 * idx + 19);

                auto t10_xy_xxzz = primBuffer.data(t10off + 90 * idx + 20);

                auto t10_xy_xyyy = primBuffer.data(t10off + 90 * idx + 21);

                auto t10_xy_xyyz = primBuffer.data(t10off + 90 * idx + 22);

                auto t10_xy_xyzz = primBuffer.data(t10off + 90 * idx + 23);

                auto t10_xy_xzzz = primBuffer.data(t10off + 90 * idx + 24);

                auto t10_xy_yyyy = primBuffer.data(t10off + 90 * idx + 25);

                auto t10_xy_yyyz = primBuffer.data(t10off + 90 * idx + 26);

                auto t10_xy_yyzz = primBuffer.data(t10off + 90 * idx + 27);

                auto t10_xy_yzzz = primBuffer.data(t10off + 90 * idx + 28);

                auto t10_xy_zzzz = primBuffer.data(t10off + 90 * idx + 29);

                auto t10_xz_xxxx = primBuffer.data(t10off + 90 * idx + 30);

                auto t10_xz_xxxy = primBuffer.data(t10off + 90 * idx + 31);

                auto t10_xz_xxxz = primBuffer.data(t10off + 90 * idx + 32);

                auto t10_xz_xxyy = primBuffer.data(t10off + 90 * idx + 33);

                auto t10_xz_xxyz = primBuffer.data(t10off + 90 * idx + 34);

                auto t10_xz_xxzz = primBuffer.data(t10off + 90 * idx + 35);

                auto t10_xz_xyyy = primBuffer.data(t10off + 90 * idx + 36);

                auto t10_xz_xyyz = primBuffer.data(t10off + 90 * idx + 37);

                auto t10_xz_xyzz = primBuffer.data(t10off + 90 * idx + 38);

                auto t10_xz_xzzz = primBuffer.data(t10off + 90 * idx + 39);

                auto t10_xz_yyyy = primBuffer.data(t10off + 90 * idx + 40);

                auto t10_xz_yyyz = primBuffer.data(t10off + 90 * idx + 41);

                auto t10_xz_yyzz = primBuffer.data(t10off + 90 * idx + 42);

                auto t10_xz_yzzz = primBuffer.data(t10off + 90 * idx + 43);

                auto t10_xz_zzzz = primBuffer.data(t10off + 90 * idx + 44);

                auto t10_yy_xxxx = primBuffer.data(t10off + 90 * idx + 45);

                auto t10_yy_xxxy = primBuffer.data(t10off + 90 * idx + 46);

                auto t10_yy_xxxz = primBuffer.data(t10off + 90 * idx + 47);

                auto t10_yy_xxyy = primBuffer.data(t10off + 90 * idx + 48);

                auto t10_yy_xxyz = primBuffer.data(t10off + 90 * idx + 49);

                auto t10_yy_xxzz = primBuffer.data(t10off + 90 * idx + 50);

                auto t10_yy_xyyy = primBuffer.data(t10off + 90 * idx + 51);

                auto t10_yy_xyyz = primBuffer.data(t10off + 90 * idx + 52);

                auto t10_yy_xyzz = primBuffer.data(t10off + 90 * idx + 53);

                auto t10_yy_xzzz = primBuffer.data(t10off + 90 * idx + 54);

                auto t10_yy_yyyy = primBuffer.data(t10off + 90 * idx + 55);

                auto t10_yy_yyyz = primBuffer.data(t10off + 90 * idx + 56);

                auto t10_yy_yyzz = primBuffer.data(t10off + 90 * idx + 57);

                auto t10_yy_yzzz = primBuffer.data(t10off + 90 * idx + 58);

                auto t10_yy_zzzz = primBuffer.data(t10off + 90 * idx + 59);

                auto t10_yz_xxxx = primBuffer.data(t10off + 90 * idx + 60);

                auto t10_yz_xxxy = primBuffer.data(t10off + 90 * idx + 61);

                auto t10_yz_xxxz = primBuffer.data(t10off + 90 * idx + 62);

                auto t10_yz_xxyy = primBuffer.data(t10off + 90 * idx + 63);

                auto t10_yz_xxyz = primBuffer.data(t10off + 90 * idx + 64);

                auto t10_yz_xxzz = primBuffer.data(t10off + 90 * idx + 65);

                auto t10_yz_xyyy = primBuffer.data(t10off + 90 * idx + 66);

                auto t10_yz_xyyz = primBuffer.data(t10off + 90 * idx + 67);

                auto t10_yz_xyzz = primBuffer.data(t10off + 90 * idx + 68);

                auto t10_yz_xzzz = primBuffer.data(t10off + 90 * idx + 69);

                auto t10_yz_yyyy = primBuffer.data(t10off + 90 * idx + 70);

                auto t10_yz_yyyz = primBuffer.data(t10off + 90 * idx + 71);

                auto t10_yz_yyzz = primBuffer.data(t10off + 90 * idx + 72);

                auto t10_yz_yzzz = primBuffer.data(t10off + 90 * idx + 73);

                auto t10_yz_zzzz = primBuffer.data(t10off + 90 * idx + 74);

                auto t10_zz_xxxx = primBuffer.data(t10off + 90 * idx + 75);

                auto t10_zz_xxxy = primBuffer.data(t10off + 90 * idx + 76);

                auto t10_zz_xxxz = primBuffer.data(t10off + 90 * idx + 77);

                auto t10_zz_xxyy = primBuffer.data(t10off + 90 * idx + 78);

                auto t10_zz_xxyz = primBuffer.data(t10off + 90 * idx + 79);

                auto t10_zz_xxzz = primBuffer.data(t10off + 90 * idx + 80);

                auto t10_zz_xyyy = primBuffer.data(t10off + 90 * idx + 81);

                auto t10_zz_xyyz = primBuffer.data(t10off + 90 * idx + 82);

                auto t10_zz_xyzz = primBuffer.data(t10off + 90 * idx + 83);

                auto t10_zz_xzzz = primBuffer.data(t10off + 90 * idx + 84);

                auto t10_zz_yyyy = primBuffer.data(t10off + 90 * idx + 85);

                auto t10_zz_yyyz = primBuffer.data(t10off + 90 * idx + 86);

                auto t10_zz_yyzz = primBuffer.data(t10off + 90 * idx + 87);

                auto t10_zz_yzzz = primBuffer.data(t10off + 90 * idx + 88);

                auto t10_zz_zzzz = primBuffer.data(t10off + 90 * idx + 89);

                // set up pointers to (D|A(0)|G)^(m+1) integrals

                auto t11_xx_xxxx = primBuffer.data(t11off + 90 * idx);

                auto t11_xx_xxxy = primBuffer.data(t11off + 90 * idx + 1);

                auto t11_xx_xxxz = primBuffer.data(t11off + 90 * idx + 2);

                auto t11_xx_xxyy = primBuffer.data(t11off + 90 * idx + 3);

                auto t11_xx_xxyz = primBuffer.data(t11off + 90 * idx + 4);

                auto t11_xx_xxzz = primBuffer.data(t11off + 90 * idx + 5);

                auto t11_xx_xyyy = primBuffer.data(t11off + 90 * idx + 6);

                auto t11_xx_xyyz = primBuffer.data(t11off + 90 * idx + 7);

                auto t11_xx_xyzz = primBuffer.data(t11off + 90 * idx + 8);

                auto t11_xx_xzzz = primBuffer.data(t11off + 90 * idx + 9);

                auto t11_xx_yyyy = primBuffer.data(t11off + 90 * idx + 10);

                auto t11_xx_yyyz = primBuffer.data(t11off + 90 * idx + 11);

                auto t11_xx_yyzz = primBuffer.data(t11off + 90 * idx + 12);

                auto t11_xx_yzzz = primBuffer.data(t11off + 90 * idx + 13);

                auto t11_xx_zzzz = primBuffer.data(t11off + 90 * idx + 14);

                auto t11_xy_xxxx = primBuffer.data(t11off + 90 * idx + 15);

                auto t11_xy_xxxy = primBuffer.data(t11off + 90 * idx + 16);

                auto t11_xy_xxxz = primBuffer.data(t11off + 90 * idx + 17);

                auto t11_xy_xxyy = primBuffer.data(t11off + 90 * idx + 18);

                auto t11_xy_xxyz = primBuffer.data(t11off + 90 * idx + 19);

                auto t11_xy_xxzz = primBuffer.data(t11off + 90 * idx + 20);

                auto t11_xy_xyyy = primBuffer.data(t11off + 90 * idx + 21);

                auto t11_xy_xyyz = primBuffer.data(t11off + 90 * idx + 22);

                auto t11_xy_xyzz = primBuffer.data(t11off + 90 * idx + 23);

                auto t11_xy_xzzz = primBuffer.data(t11off + 90 * idx + 24);

                auto t11_xy_yyyy = primBuffer.data(t11off + 90 * idx + 25);

                auto t11_xy_yyyz = primBuffer.data(t11off + 90 * idx + 26);

                auto t11_xy_yyzz = primBuffer.data(t11off + 90 * idx + 27);

                auto t11_xy_yzzz = primBuffer.data(t11off + 90 * idx + 28);

                auto t11_xy_zzzz = primBuffer.data(t11off + 90 * idx + 29);

                auto t11_xz_xxxx = primBuffer.data(t11off + 90 * idx + 30);

                auto t11_xz_xxxy = primBuffer.data(t11off + 90 * idx + 31);

                auto t11_xz_xxxz = primBuffer.data(t11off + 90 * idx + 32);

                auto t11_xz_xxyy = primBuffer.data(t11off + 90 * idx + 33);

                auto t11_xz_xxyz = primBuffer.data(t11off + 90 * idx + 34);

                auto t11_xz_xxzz = primBuffer.data(t11off + 90 * idx + 35);

                auto t11_xz_xyyy = primBuffer.data(t11off + 90 * idx + 36);

                auto t11_xz_xyyz = primBuffer.data(t11off + 90 * idx + 37);

                auto t11_xz_xyzz = primBuffer.data(t11off + 90 * idx + 38);

                auto t11_xz_xzzz = primBuffer.data(t11off + 90 * idx + 39);

                auto t11_xz_yyyy = primBuffer.data(t11off + 90 * idx + 40);

                auto t11_xz_yyyz = primBuffer.data(t11off + 90 * idx + 41);

                auto t11_xz_yyzz = primBuffer.data(t11off + 90 * idx + 42);

                auto t11_xz_yzzz = primBuffer.data(t11off + 90 * idx + 43);

                auto t11_xz_zzzz = primBuffer.data(t11off + 90 * idx + 44);

                auto t11_yy_xxxx = primBuffer.data(t11off + 90 * idx + 45);

                auto t11_yy_xxxy = primBuffer.data(t11off + 90 * idx + 46);

                auto t11_yy_xxxz = primBuffer.data(t11off + 90 * idx + 47);

                auto t11_yy_xxyy = primBuffer.data(t11off + 90 * idx + 48);

                auto t11_yy_xxyz = primBuffer.data(t11off + 90 * idx + 49);

                auto t11_yy_xxzz = primBuffer.data(t11off + 90 * idx + 50);

                auto t11_yy_xyyy = primBuffer.data(t11off + 90 * idx + 51);

                auto t11_yy_xyyz = primBuffer.data(t11off + 90 * idx + 52);

                auto t11_yy_xyzz = primBuffer.data(t11off + 90 * idx + 53);

                auto t11_yy_xzzz = primBuffer.data(t11off + 90 * idx + 54);

                auto t11_yy_yyyy = primBuffer.data(t11off + 90 * idx + 55);

                auto t11_yy_yyyz = primBuffer.data(t11off + 90 * idx + 56);

                auto t11_yy_yyzz = primBuffer.data(t11off + 90 * idx + 57);

                auto t11_yy_yzzz = primBuffer.data(t11off + 90 * idx + 58);

                auto t11_yy_zzzz = primBuffer.data(t11off + 90 * idx + 59);

                auto t11_yz_xxxx = primBuffer.data(t11off + 90 * idx + 60);

                auto t11_yz_xxxy = primBuffer.data(t11off + 90 * idx + 61);

                auto t11_yz_xxxz = primBuffer.data(t11off + 90 * idx + 62);

                auto t11_yz_xxyy = primBuffer.data(t11off + 90 * idx + 63);

                auto t11_yz_xxyz = primBuffer.data(t11off + 90 * idx + 64);

                auto t11_yz_xxzz = primBuffer.data(t11off + 90 * idx + 65);

                auto t11_yz_xyyy = primBuffer.data(t11off + 90 * idx + 66);

                auto t11_yz_xyyz = primBuffer.data(t11off + 90 * idx + 67);

                auto t11_yz_xyzz = primBuffer.data(t11off + 90 * idx + 68);

                auto t11_yz_xzzz = primBuffer.data(t11off + 90 * idx + 69);

                auto t11_yz_yyyy = primBuffer.data(t11off + 90 * idx + 70);

                auto t11_yz_yyyz = primBuffer.data(t11off + 90 * idx + 71);

                auto t11_yz_yyzz = primBuffer.data(t11off + 90 * idx + 72);

                auto t11_yz_yzzz = primBuffer.data(t11off + 90 * idx + 73);

                auto t11_yz_zzzz = primBuffer.data(t11off + 90 * idx + 74);

                auto t11_zz_xxxx = primBuffer.data(t11off + 90 * idx + 75);

                auto t11_zz_xxxy = primBuffer.data(t11off + 90 * idx + 76);

                auto t11_zz_xxxz = primBuffer.data(t11off + 90 * idx + 77);

                auto t11_zz_xxyy = primBuffer.data(t11off + 90 * idx + 78);

                auto t11_zz_xxyz = primBuffer.data(t11off + 90 * idx + 79);

                auto t11_zz_xxzz = primBuffer.data(t11off + 90 * idx + 80);

                auto t11_zz_xyyy = primBuffer.data(t11off + 90 * idx + 81);

                auto t11_zz_xyyz = primBuffer.data(t11off + 90 * idx + 82);

                auto t11_zz_xyzz = primBuffer.data(t11off + 90 * idx + 83);

                auto t11_zz_xzzz = primBuffer.data(t11off + 90 * idx + 84);

                auto t11_zz_yyyy = primBuffer.data(t11off + 90 * idx + 85);

                auto t11_zz_yyyz = primBuffer.data(t11off + 90 * idx + 86);

                auto t11_zz_yyzz = primBuffer.data(t11off + 90 * idx + 87);

                auto t11_zz_yzzz = primBuffer.data(t11off + 90 * idx + 88);

                auto t11_zz_zzzz = primBuffer.data(t11off + 90 * idx + 89);

                // set up pointers to (F|A(0)|G)^(m) integrals

                auto t_xxx_xxxx = primBuffer.data(toff + 150 * idx);

                auto t_xxx_xxxy = primBuffer.data(toff + 150 * idx + 1);

                auto t_xxx_xxxz = primBuffer.data(toff + 150 * idx + 2);

                auto t_xxx_xxyy = primBuffer.data(toff + 150 * idx + 3);

                auto t_xxx_xxyz = primBuffer.data(toff + 150 * idx + 4);

                auto t_xxx_xxzz = primBuffer.data(toff + 150 * idx + 5);

                auto t_xxx_xyyy = primBuffer.data(toff + 150 * idx + 6);

                auto t_xxx_xyyz = primBuffer.data(toff + 150 * idx + 7);

                auto t_xxx_xyzz = primBuffer.data(toff + 150 * idx + 8);

                auto t_xxx_xzzz = primBuffer.data(toff + 150 * idx + 9);

                auto t_xxx_yyyy = primBuffer.data(toff + 150 * idx + 10);

                auto t_xxx_yyyz = primBuffer.data(toff + 150 * idx + 11);

                auto t_xxx_yyzz = primBuffer.data(toff + 150 * idx + 12);

                auto t_xxx_yzzz = primBuffer.data(toff + 150 * idx + 13);

                auto t_xxx_zzzz = primBuffer.data(toff + 150 * idx + 14);

                auto t_xxy_xxxx = primBuffer.data(toff + 150 * idx + 15);

                auto t_xxy_xxxy = primBuffer.data(toff + 150 * idx + 16);

                auto t_xxy_xxxz = primBuffer.data(toff + 150 * idx + 17);

                auto t_xxy_xxyy = primBuffer.data(toff + 150 * idx + 18);

                auto t_xxy_xxyz = primBuffer.data(toff + 150 * idx + 19);

                auto t_xxy_xxzz = primBuffer.data(toff + 150 * idx + 20);

                auto t_xxy_xyyy = primBuffer.data(toff + 150 * idx + 21);

                auto t_xxy_xyyz = primBuffer.data(toff + 150 * idx + 22);

                auto t_xxy_xyzz = primBuffer.data(toff + 150 * idx + 23);

                auto t_xxy_xzzz = primBuffer.data(toff + 150 * idx + 24);

                auto t_xxy_yyyy = primBuffer.data(toff + 150 * idx + 25);

                auto t_xxy_yyyz = primBuffer.data(toff + 150 * idx + 26);

                auto t_xxy_yyzz = primBuffer.data(toff + 150 * idx + 27);

                auto t_xxy_yzzz = primBuffer.data(toff + 150 * idx + 28);

                auto t_xxy_zzzz = primBuffer.data(toff + 150 * idx + 29);

                auto t_xxz_xxxx = primBuffer.data(toff + 150 * idx + 30);

                auto t_xxz_xxxy = primBuffer.data(toff + 150 * idx + 31);

                auto t_xxz_xxxz = primBuffer.data(toff + 150 * idx + 32);

                auto t_xxz_xxyy = primBuffer.data(toff + 150 * idx + 33);

                auto t_xxz_xxyz = primBuffer.data(toff + 150 * idx + 34);

                auto t_xxz_xxzz = primBuffer.data(toff + 150 * idx + 35);

                auto t_xxz_xyyy = primBuffer.data(toff + 150 * idx + 36);

                auto t_xxz_xyyz = primBuffer.data(toff + 150 * idx + 37);

                auto t_xxz_xyzz = primBuffer.data(toff + 150 * idx + 38);

                auto t_xxz_xzzz = primBuffer.data(toff + 150 * idx + 39);

                auto t_xxz_yyyy = primBuffer.data(toff + 150 * idx + 40);

                auto t_xxz_yyyz = primBuffer.data(toff + 150 * idx + 41);

                auto t_xxz_yyzz = primBuffer.data(toff + 150 * idx + 42);

                auto t_xxz_yzzz = primBuffer.data(toff + 150 * idx + 43);

                auto t_xxz_zzzz = primBuffer.data(toff + 150 * idx + 44);

                auto t_xyy_xxxx = primBuffer.data(toff + 150 * idx + 45);

                auto t_xyy_xxxy = primBuffer.data(toff + 150 * idx + 46);

                auto t_xyy_xxxz = primBuffer.data(toff + 150 * idx + 47);

                auto t_xyy_xxyy = primBuffer.data(toff + 150 * idx + 48);

                auto t_xyy_xxyz = primBuffer.data(toff + 150 * idx + 49);

                auto t_xyy_xxzz = primBuffer.data(toff + 150 * idx + 50);

                auto t_xyy_xyyy = primBuffer.data(toff + 150 * idx + 51);

                auto t_xyy_xyyz = primBuffer.data(toff + 150 * idx + 52);

                auto t_xyy_xyzz = primBuffer.data(toff + 150 * idx + 53);

                auto t_xyy_xzzz = primBuffer.data(toff + 150 * idx + 54);

                auto t_xyy_yyyy = primBuffer.data(toff + 150 * idx + 55);

                auto t_xyy_yyyz = primBuffer.data(toff + 150 * idx + 56);

                auto t_xyy_yyzz = primBuffer.data(toff + 150 * idx + 57);

                auto t_xyy_yzzz = primBuffer.data(toff + 150 * idx + 58);

                auto t_xyy_zzzz = primBuffer.data(toff + 150 * idx + 59);

                auto t_xyz_xxxx = primBuffer.data(toff + 150 * idx + 60);

                auto t_xyz_xxxy = primBuffer.data(toff + 150 * idx + 61);

                auto t_xyz_xxxz = primBuffer.data(toff + 150 * idx + 62);

                auto t_xyz_xxyy = primBuffer.data(toff + 150 * idx + 63);

                auto t_xyz_xxyz = primBuffer.data(toff + 150 * idx + 64);

                auto t_xyz_xxzz = primBuffer.data(toff + 150 * idx + 65);

                auto t_xyz_xyyy = primBuffer.data(toff + 150 * idx + 66);

                auto t_xyz_xyyz = primBuffer.data(toff + 150 * idx + 67);

                auto t_xyz_xyzz = primBuffer.data(toff + 150 * idx + 68);

                auto t_xyz_xzzz = primBuffer.data(toff + 150 * idx + 69);

                auto t_xyz_yyyy = primBuffer.data(toff + 150 * idx + 70);

                auto t_xyz_yyyz = primBuffer.data(toff + 150 * idx + 71);

                auto t_xyz_yyzz = primBuffer.data(toff + 150 * idx + 72);

                auto t_xyz_yzzz = primBuffer.data(toff + 150 * idx + 73);

                auto t_xyz_zzzz = primBuffer.data(toff + 150 * idx + 74);

                auto t_xzz_xxxx = primBuffer.data(toff + 150 * idx + 75);

                auto t_xzz_xxxy = primBuffer.data(toff + 150 * idx + 76);

                auto t_xzz_xxxz = primBuffer.data(toff + 150 * idx + 77);

                auto t_xzz_xxyy = primBuffer.data(toff + 150 * idx + 78);

                auto t_xzz_xxyz = primBuffer.data(toff + 150 * idx + 79);

                auto t_xzz_xxzz = primBuffer.data(toff + 150 * idx + 80);

                auto t_xzz_xyyy = primBuffer.data(toff + 150 * idx + 81);

                auto t_xzz_xyyz = primBuffer.data(toff + 150 * idx + 82);

                auto t_xzz_xyzz = primBuffer.data(toff + 150 * idx + 83);

                auto t_xzz_xzzz = primBuffer.data(toff + 150 * idx + 84);

                auto t_xzz_yyyy = primBuffer.data(toff + 150 * idx + 85);

                auto t_xzz_yyyz = primBuffer.data(toff + 150 * idx + 86);

                auto t_xzz_yyzz = primBuffer.data(toff + 150 * idx + 87);

                auto t_xzz_yzzz = primBuffer.data(toff + 150 * idx + 88);

                auto t_xzz_zzzz = primBuffer.data(toff + 150 * idx + 89);

                auto t_yyy_xxxx = primBuffer.data(toff + 150 * idx + 90);

                auto t_yyy_xxxy = primBuffer.data(toff + 150 * idx + 91);

                auto t_yyy_xxxz = primBuffer.data(toff + 150 * idx + 92);

                auto t_yyy_xxyy = primBuffer.data(toff + 150 * idx + 93);

                auto t_yyy_xxyz = primBuffer.data(toff + 150 * idx + 94);

                auto t_yyy_xxzz = primBuffer.data(toff + 150 * idx + 95);

                auto t_yyy_xyyy = primBuffer.data(toff + 150 * idx + 96);

                auto t_yyy_xyyz = primBuffer.data(toff + 150 * idx + 97);

                auto t_yyy_xyzz = primBuffer.data(toff + 150 * idx + 98);

                auto t_yyy_xzzz = primBuffer.data(toff + 150 * idx + 99);

                auto t_yyy_yyyy = primBuffer.data(toff + 150 * idx + 100);

                auto t_yyy_yyyz = primBuffer.data(toff + 150 * idx + 101);

                auto t_yyy_yyzz = primBuffer.data(toff + 150 * idx + 102);

                auto t_yyy_yzzz = primBuffer.data(toff + 150 * idx + 103);

                auto t_yyy_zzzz = primBuffer.data(toff + 150 * idx + 104);

                auto t_yyz_xxxx = primBuffer.data(toff + 150 * idx + 105);

                auto t_yyz_xxxy = primBuffer.data(toff + 150 * idx + 106);

                auto t_yyz_xxxz = primBuffer.data(toff + 150 * idx + 107);

                auto t_yyz_xxyy = primBuffer.data(toff + 150 * idx + 108);

                auto t_yyz_xxyz = primBuffer.data(toff + 150 * idx + 109);

                auto t_yyz_xxzz = primBuffer.data(toff + 150 * idx + 110);

                auto t_yyz_xyyy = primBuffer.data(toff + 150 * idx + 111);

                auto t_yyz_xyyz = primBuffer.data(toff + 150 * idx + 112);

                auto t_yyz_xyzz = primBuffer.data(toff + 150 * idx + 113);

                auto t_yyz_xzzz = primBuffer.data(toff + 150 * idx + 114);

                auto t_yyz_yyyy = primBuffer.data(toff + 150 * idx + 115);

                auto t_yyz_yyyz = primBuffer.data(toff + 150 * idx + 116);

                auto t_yyz_yyzz = primBuffer.data(toff + 150 * idx + 117);

                auto t_yyz_yzzz = primBuffer.data(toff + 150 * idx + 118);

                auto t_yyz_zzzz = primBuffer.data(toff + 150 * idx + 119);

                auto t_yzz_xxxx = primBuffer.data(toff + 150 * idx + 120);

                auto t_yzz_xxxy = primBuffer.data(toff + 150 * idx + 121);

                auto t_yzz_xxxz = primBuffer.data(toff + 150 * idx + 122);

                auto t_yzz_xxyy = primBuffer.data(toff + 150 * idx + 123);

                auto t_yzz_xxyz = primBuffer.data(toff + 150 * idx + 124);

                auto t_yzz_xxzz = primBuffer.data(toff + 150 * idx + 125);

                auto t_yzz_xyyy = primBuffer.data(toff + 150 * idx + 126);

                auto t_yzz_xyyz = primBuffer.data(toff + 150 * idx + 127);

                auto t_yzz_xyzz = primBuffer.data(toff + 150 * idx + 128);

                auto t_yzz_xzzz = primBuffer.data(toff + 150 * idx + 129);

                auto t_yzz_yyyy = primBuffer.data(toff + 150 * idx + 130);

                auto t_yzz_yyyz = primBuffer.data(toff + 150 * idx + 131);

                auto t_yzz_yyzz = primBuffer.data(toff + 150 * idx + 132);

                auto t_yzz_yzzz = primBuffer.data(toff + 150 * idx + 133);

                auto t_yzz_zzzz = primBuffer.data(toff + 150 * idx + 134);

                auto t_zzz_xxxx = primBuffer.data(toff + 150 * idx + 135);

                auto t_zzz_xxxy = primBuffer.data(toff + 150 * idx + 136);

                auto t_zzz_xxxz = primBuffer.data(toff + 150 * idx + 137);

                auto t_zzz_xxyy = primBuffer.data(toff + 150 * idx + 138);

                auto t_zzz_xxyz = primBuffer.data(toff + 150 * idx + 139);

                auto t_zzz_xxzz = primBuffer.data(toff + 150 * idx + 140);

                auto t_zzz_xyyy = primBuffer.data(toff + 150 * idx + 141);

                auto t_zzz_xyyz = primBuffer.data(toff + 150 * idx + 142);

                auto t_zzz_xyzz = primBuffer.data(toff + 150 * idx + 143);

                auto t_zzz_xzzz = primBuffer.data(toff + 150 * idx + 144);

                auto t_zzz_yyyy = primBuffer.data(toff + 150 * idx + 145);

                auto t_zzz_yyyz = primBuffer.data(toff + 150 * idx + 146);

                auto t_zzz_yyzz = primBuffer.data(toff + 150 * idx + 147);

                auto t_zzz_yzzz = primBuffer.data(toff + 150 * idx + 148);

                auto t_zzz_zzzz = primBuffer.data(toff + 150 * idx + 149);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xx_xxx, tk0_xx_xxy,\
                                         tk0_xx_xxz, tk0_xx_xyy, tk0_xx_xyz, tk0_xx_xzz,\
                                         tk0_xx_yyy, tk0_xx_yyz, tk0_xx_yzz, tk0_xx_zzz,\
                                         tk0_xy_xxx, tk0_xy_xxy, tk0_xy_xxz, tk0_xy_xyy,\
                                         tk0_xy_xyz, tk0_xy_xzz, tk0_xy_yyy, tk0_xy_yyz,\
                                         tk0_xy_yzz, tk0_xy_zzz, tk0_xz_xxx, tk0_xz_xxy,\
                                         tk0_xz_xxz, tk0_xz_xyy, tk0_xz_xyz, tk0_xz_xzz,\
                                         tk0_xz_yyy, tk0_xz_yyz, tk0_xz_yzz, tk0_xz_zzz,\
                                         tk0_yy_xxx, tk0_yy_xxy, tk0_yy_xxz, tk0_yy_xyy,\
                                         tk0_yy_xyz, tk0_yy_xzz, tk0_yy_yyy, tk0_yy_yyz,\
                                         tk0_yy_yzz, tk0_yy_zzz, tk0_yz_xxx, tk0_yz_xxy,\
                                         tk0_yz_xxz, tk0_yz_xyy, tk0_yz_xyz, tk0_yz_xzz,\
                                         tk0_yz_yyy, tk0_yz_yyz, tk0_yz_yzz, tk0_yz_zzz,\
                                         tk0_zz_xxx, tk0_zz_xxy, tk0_zz_xxz, tk0_zz_xyy,\
                                         tk0_zz_xyz, tk0_zz_xzz, tk0_zz_yyy, tk0_zz_yyz,\
                                         tk0_zz_yzz, tk0_zz_zzz, tk1_xx_xxx, tk1_xx_xxy,\
                                         tk1_xx_xxz, tk1_xx_xyy, tk1_xx_xyz, tk1_xx_xzz,\
                                         tk1_xx_yyy, tk1_xx_yyz, tk1_xx_yzz, tk1_xx_zzz,\
                                         tk1_xy_xxx, tk1_xy_xxy, tk1_xy_xxz, tk1_xy_xyy,\
                                         tk1_xy_xyz, tk1_xy_xzz, tk1_xy_yyy, tk1_xy_yyz,\
                                         tk1_xy_yzz, tk1_xy_zzz, tk1_xz_xxx, tk1_xz_xxy,\
                                         tk1_xz_xxz, tk1_xz_xyy, tk1_xz_xyz, tk1_xz_xzz,\
                                         tk1_xz_yyy, tk1_xz_yyz, tk1_xz_yzz, tk1_xz_zzz,\
                                         tk1_yy_xxx, tk1_yy_xxy, tk1_yy_xxz, tk1_yy_xyy,\
                                         tk1_yy_xyz, tk1_yy_xzz, tk1_yy_yyy, tk1_yy_yyz,\
                                         tk1_yy_yzz, tk1_yy_zzz, tk1_yz_xxx, tk1_yz_xxy,\
                                         tk1_yz_xxz, tk1_yz_xyy, tk1_yz_xyz, tk1_yz_xzz,\
                                         tk1_yz_yyy, tk1_yz_yyz, tk1_yz_yzz, tk1_yz_zzz,\
                                         tk1_zz_xxx, tk1_zz_xxy, tk1_zz_xxz, tk1_zz_xyy,\
                                         tk1_zz_xyz, tk1_zz_xzz, tk1_zz_yyy, tk1_zz_yyz,\
                                         tk1_zz_yzz, tk1_zz_zzz, t20_x_xxxx, t20_x_xxxy,\
                                         t20_x_xxxz, t20_x_xxyy, t20_x_xxyz, t20_x_xxzz,\
                                         t20_x_xyyy, t20_x_xyyz, t20_x_xyzz, t20_x_xzzz,\
                                         t20_x_yyyy, t20_x_yyyz, t20_x_yyzz, t20_x_yzzz,\
                                         t20_x_zzzz, t20_y_xxxx, t20_y_xxxy, t20_y_xxxz,\
                                         t20_y_xxyy, t20_y_xxyz, t20_y_xxzz, t20_y_xyyy,\
                                         t20_y_xyyz, t20_y_xyzz, t20_y_xzzz, t20_y_yyyy,\
                                         t20_y_yyyz, t20_y_yyzz, t20_y_yzzz, t20_y_zzzz,\
                                         t20_z_xxxx, t20_z_xxxy, t20_z_xxxz, t20_z_xxyy,\
                                         t20_z_xxyz, t20_z_xxzz, t20_z_xyyy, t20_z_xyyz,\
                                         t20_z_xyzz, t20_z_xzzz, t20_z_yyyy, t20_z_yyyz,\
                                         t20_z_yyzz, t20_z_yzzz, t20_z_zzzz, t21_x_xxxx,\
                                         t21_x_xxxy, t21_x_xxxz, t21_x_xxyy, t21_x_xxyz,\
                                         t21_x_xxzz, t21_x_xyyy, t21_x_xyyz, t21_x_xyzz,\
                                         t21_x_xzzz, t21_x_yyyy, t21_x_yyyz, t21_x_yyzz,\
                                         t21_x_yzzz, t21_x_zzzz, t21_y_xxxx, t21_y_xxxy,\
                                         t21_y_xxxz, t21_y_xxyy, t21_y_xxyz, t21_y_xxzz,\
                                         t21_y_xyyy, t21_y_xyyz, t21_y_xyzz, t21_y_xzzz,\
                                         t21_y_yyyy, t21_y_yyyz, t21_y_yyzz, t21_y_yzzz,\
                                         t21_y_zzzz, t21_z_xxxx, t21_z_xxxy, t21_z_xxxz,\
                                         t21_z_xxyy, t21_z_xxyz, t21_z_xxzz, t21_z_xyyy,\
                                         t21_z_xyyz, t21_z_xyzz, t21_z_xzzz, t21_z_yyyy,\
                                         t21_z_yyyz, t21_z_yyzz, t21_z_yzzz, t21_z_zzzz,\
                                         t10_xx_xxxx, t10_xx_xxxy, t10_xx_xxxz,\
                                         t10_xx_xxyy, t10_xx_xxyz, t10_xx_xxzz,\
                                         t10_xx_xyyy, t10_xx_xyyz, t10_xx_xyzz,\
                                         t10_xx_xzzz, t10_xx_yyyy, t10_xx_yyyz,\
                                         t10_xx_yyzz, t10_xx_yzzz, t10_xx_zzzz,\
                                         t10_xy_xxxx, t10_xy_xxxy, t10_xy_xxxz,\
                                         t10_xy_xxyy, t10_xy_xxyz, t10_xy_xxzz,\
                                         t10_xy_xyyy, t10_xy_xyyz, t10_xy_xyzz,\
                                         t10_xy_xzzz, t10_xy_yyyy, t10_xy_yyyz,\
                                         t10_xy_yyzz, t10_xy_yzzz, t10_xy_zzzz,\
                                         t10_xz_xxxx, t10_xz_xxxy, t10_xz_xxxz,\
                                         t10_xz_xxyy, t10_xz_xxyz, t10_xz_xxzz,\
                                         t10_xz_xyyy, t10_xz_xyyz, t10_xz_xyzz,\
                                         t10_xz_xzzz, t10_xz_yyyy, t10_xz_yyyz,\
                                         t10_xz_yyzz, t10_xz_yzzz, t10_xz_zzzz,\
                                         t10_yy_xxxx, t10_yy_xxxy, t10_yy_xxxz,\
                                         t10_yy_xxyy, t10_yy_xxyz, t10_yy_xxzz,\
                                         t10_yy_xyyy, t10_yy_xyyz, t10_yy_xyzz,\
                                         t10_yy_xzzz, t10_yy_yyyy, t10_yy_yyyz,\
                                         t10_yy_yyzz, t10_yy_yzzz, t10_yy_zzzz,\
                                         t10_yz_xxxx, t10_yz_xxxy, t10_yz_xxxz,\
                                         t10_yz_xxyy, t10_yz_xxyz, t10_yz_xxzz,\
                                         t10_yz_xyyy, t10_yz_xyyz, t10_yz_xyzz,\
                                         t10_yz_xzzz, t10_yz_yyyy, t10_yz_yyyz,\
                                         t10_yz_yyzz, t10_yz_yzzz, t10_yz_zzzz,\
                                         t10_zz_xxxx, t10_zz_xxxy, t10_zz_xxxz,\
                                         t10_zz_xxyy, t10_zz_xxyz, t10_zz_xxzz,\
                                         t10_zz_xyyy, t10_zz_xyyz, t10_zz_xyzz,\
                                         t10_zz_xzzz, t10_zz_yyyy, t10_zz_yyyz,\
                                         t10_zz_yyzz, t10_zz_yzzz, t10_zz_zzzz,\
                                         t11_xx_xxxx, t11_xx_xxxy, t11_xx_xxxz,\
                                         t11_xx_xxyy, t11_xx_xxyz, t11_xx_xxzz,\
                                         t11_xx_xyyy, t11_xx_xyyz, t11_xx_xyzz,\
                                         t11_xx_xzzz, t11_xx_yyyy, t11_xx_yyyz,\
                                         t11_xx_yyzz, t11_xx_yzzz, t11_xx_zzzz,\
                                         t11_xy_xxxx, t11_xy_xxxy, t11_xy_xxxz,\
                                         t11_xy_xxyy, t11_xy_xxyz, t11_xy_xxzz,\
                                         t11_xy_xyyy, t11_xy_xyyz, t11_xy_xyzz,\
                                         t11_xy_xzzz, t11_xy_yyyy, t11_xy_yyyz,\
                                         t11_xy_yyzz, t11_xy_yzzz, t11_xy_zzzz,\
                                         t11_xz_xxxx, t11_xz_xxxy, t11_xz_xxxz,\
                                         t11_xz_xxyy, t11_xz_xxyz, t11_xz_xxzz,\
                                         t11_xz_xyyy, t11_xz_xyyz, t11_xz_xyzz,\
                                         t11_xz_xzzz, t11_xz_yyyy, t11_xz_yyyz,\
                                         t11_xz_yyzz, t11_xz_yzzz, t11_xz_zzzz,\
                                         t11_yy_xxxx, t11_yy_xxxy, t11_yy_xxxz,\
                                         t11_yy_xxyy, t11_yy_xxyz, t11_yy_xxzz,\
                                         t11_yy_xyyy, t11_yy_xyyz, t11_yy_xyzz,\
                                         t11_yy_xzzz, t11_yy_yyyy, t11_yy_yyyz,\
                                         t11_yy_yyzz, t11_yy_yzzz, t11_yy_zzzz,\
                                         t11_yz_xxxx, t11_yz_xxxy, t11_yz_xxxz,\
                                         t11_yz_xxyy, t11_yz_xxyz, t11_yz_xxzz,\
                                         t11_yz_xyyy, t11_yz_xyyz, t11_yz_xyzz,\
                                         t11_yz_xzzz, t11_yz_yyyy, t11_yz_yyyz,\
                                         t11_yz_yyzz, t11_yz_yzzz, t11_yz_zzzz,\
                                         t11_zz_xxxx, t11_zz_xxxy, t11_zz_xxxz,\
                                         t11_zz_xxyy, t11_zz_xxyz, t11_zz_xxzz,\
                                         t11_zz_xyyy, t11_zz_xyyz, t11_zz_xyzz,\
                                         t11_zz_xzzz, t11_zz_yyyy, t11_zz_yyyz,\
                                         t11_zz_yyzz, t11_zz_yzzz, t11_zz_zzzz,\
                                         t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, t_xxx_xxyy,\
                                         t_xxx_xxyz, t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz,\
                                         t_xxx_xyzz, t_xxx_xzzz, t_xxx_yyyy, t_xxx_yyyz,\
                                         t_xxx_yyzz, t_xxx_yzzz, t_xxx_zzzz, t_xxy_xxxx,\
                                         t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, t_xxy_xxyz,\
                                         t_xxy_xxzz, t_xxy_xyyy, t_xxy_xyyz, t_xxy_xyzz,\
                                         t_xxy_xzzz, t_xxy_yyyy, t_xxy_yyyz, t_xxy_yyzz,\
                                         t_xxy_yzzz, t_xxy_zzzz, t_xxz_xxxx, t_xxz_xxxy,\
                                         t_xxz_xxxz, t_xxz_xxyy, t_xxz_xxyz, t_xxz_xxzz,\
                                         t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz,\
                                         t_xxz_yyyy, t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz,\
                                         t_xxz_zzzz, t_xyy_xxxx, t_xyy_xxxy, t_xyy_xxxz,\
                                         t_xyy_xxyy, t_xyy_xxyz, t_xyy_xxzz, t_xyy_xyyy,\
                                         t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz, t_xyy_yyyy,\
                                         t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz,\
                                         t_xyz_xxxx, t_xyz_xxxy, t_xyz_xxxz, t_xyz_xxyy,\
                                         t_xyz_xxyz, t_xyz_xxzz, t_xyz_xyyy, t_xyz_xyyz,\
                                         t_xyz_xyzz, t_xyz_xzzz, t_xyz_yyyy, t_xyz_yyyz,\
                                         t_xyz_yyzz, t_xyz_yzzz, t_xyz_zzzz, t_xzz_xxxx,\
                                         t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, t_xzz_xxyz,\
                                         t_xzz_xxzz, t_xzz_xyyy, t_xzz_xyyz, t_xzz_xyzz,\
                                         t_xzz_xzzz, t_xzz_yyyy, t_xzz_yyyz, t_xzz_yyzz,\
                                         t_xzz_yzzz, t_xzz_zzzz, t_yyy_xxxx, t_yyy_xxxy,\
                                         t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz, t_yyy_xxzz,\
                                         t_yyy_xyyy, t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz,\
                                         t_yyy_yyyy, t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz,\
                                         t_yyy_zzzz, t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz,\
                                         t_yyz_xxyy, t_yyz_xxyz, t_yyz_xxzz, t_yyz_xyyy,\
                                         t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz, t_yyz_yyyy,\
                                         t_yyz_yyyz, t_yyz_yyzz, t_yyz_yzzz, t_yyz_zzzz,\
                                         t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, t_yzz_xxyy,\
                                         t_yzz_xxyz, t_yzz_xxzz, t_yzz_xyyy, t_yzz_xyyz,\
                                         t_yzz_xyzz, t_yzz_xzzz, t_yzz_yyyy, t_yzz_yyyz,\
                                         t_yzz_yyzz, t_yzz_yzzz, t_yzz_zzzz, t_zzz_xxxx,\
                                         t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz,\
                                         t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz,\
                                         t_zzz_xzzz, t_zzz_yyyy, t_zzz_yyyz, t_zzz_yyzz,\
                                         t_zzz_yzzz, t_zzz_zzzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxx_xxxx[k] = fra * t10_xx_xxxx[k] - frc * t11_xx_xxxx[k] + f2t * (2.0 * t20_x_xxxx[k] - 2.0 * t21_x_xxxx[k] + 4.0 * tk0_xx_xxx[k] - 4.0 * tk1_xx_xxx[k]);

                    t_xxx_xxxy[k] = fra * t10_xx_xxxy[k] - frc * t11_xx_xxxy[k] + f2t * (2.0 * t20_x_xxxy[k] - 2.0 * t21_x_xxxy[k] + 3.0 * tk0_xx_xxy[k] - 3.0 * tk1_xx_xxy[k]);

                    t_xxx_xxxz[k] = fra * t10_xx_xxxz[k] - frc * t11_xx_xxxz[k] + f2t * (2.0 * t20_x_xxxz[k] - 2.0 * t21_x_xxxz[k] + 3.0 * tk0_xx_xxz[k] - 3.0 * tk1_xx_xxz[k]);

                    t_xxx_xxyy[k] = fra * t10_xx_xxyy[k] - frc * t11_xx_xxyy[k] + f2t * (2.0 * t20_x_xxyy[k] - 2.0 * t21_x_xxyy[k] + 2.0 * tk0_xx_xyy[k] - 2.0 * tk1_xx_xyy[k]);

                    t_xxx_xxyz[k] = fra * t10_xx_xxyz[k] - frc * t11_xx_xxyz[k] + f2t * (2.0 * t20_x_xxyz[k] - 2.0 * t21_x_xxyz[k] + 2.0 * tk0_xx_xyz[k] - 2.0 * tk1_xx_xyz[k]);

                    t_xxx_xxzz[k] = fra * t10_xx_xxzz[k] - frc * t11_xx_xxzz[k] + f2t * (2.0 * t20_x_xxzz[k] - 2.0 * t21_x_xxzz[k] + 2.0 * tk0_xx_xzz[k] - 2.0 * tk1_xx_xzz[k]);

                    t_xxx_xyyy[k] = fra * t10_xx_xyyy[k] - frc * t11_xx_xyyy[k] + f2t * (2.0 * t20_x_xyyy[k] - 2.0 * t21_x_xyyy[k] + tk0_xx_yyy[k] - tk1_xx_yyy[k]);

                    t_xxx_xyyz[k] = fra * t10_xx_xyyz[k] - frc * t11_xx_xyyz[k] + f2t * (2.0 * t20_x_xyyz[k] - 2.0 * t21_x_xyyz[k] + tk0_xx_yyz[k] - tk1_xx_yyz[k]);

                    t_xxx_xyzz[k] = fra * t10_xx_xyzz[k] - frc * t11_xx_xyzz[k] + f2t * (2.0 * t20_x_xyzz[k] - 2.0 * t21_x_xyzz[k] + tk0_xx_yzz[k] - tk1_xx_yzz[k]);

                    t_xxx_xzzz[k] = fra * t10_xx_xzzz[k] - frc * t11_xx_xzzz[k] + f2t * (2.0 * t20_x_xzzz[k] - 2.0 * t21_x_xzzz[k] + tk0_xx_zzz[k] - tk1_xx_zzz[k]);

                    t_xxx_yyyy[k] = fra * t10_xx_yyyy[k] - frc * t11_xx_yyyy[k] + f2t * (2.0 * t20_x_yyyy[k] - 2.0 * t21_x_yyyy[k]);

                    t_xxx_yyyz[k] = fra * t10_xx_yyyz[k] - frc * t11_xx_yyyz[k] + f2t * (2.0 * t20_x_yyyz[k] - 2.0 * t21_x_yyyz[k]);

                    t_xxx_yyzz[k] = fra * t10_xx_yyzz[k] - frc * t11_xx_yyzz[k] + f2t * (2.0 * t20_x_yyzz[k] - 2.0 * t21_x_yyzz[k]);

                    t_xxx_yzzz[k] = fra * t10_xx_yzzz[k] - frc * t11_xx_yzzz[k] + f2t * (2.0 * t20_x_yzzz[k] - 2.0 * t21_x_yzzz[k]);

                    t_xxx_zzzz[k] = fra * t10_xx_zzzz[k] - frc * t11_xx_zzzz[k] + f2t * (2.0 * t20_x_zzzz[k] - 2.0 * t21_x_zzzz[k]);

                    t_xxy_xxxx[k] = fra * t10_xy_xxxx[k] - frc * t11_xy_xxxx[k] + f2t * (t20_y_xxxx[k] - t21_y_xxxx[k] + 4.0 * tk0_xy_xxx[k] - 4.0 * tk1_xy_xxx[k]);

                    t_xxy_xxxy[k] = fra * t10_xy_xxxy[k] - frc * t11_xy_xxxy[k] + f2t * (t20_y_xxxy[k] - t21_y_xxxy[k] + 3.0 * tk0_xy_xxy[k] - 3.0 * tk1_xy_xxy[k]);

                    t_xxy_xxxz[k] = fra * t10_xy_xxxz[k] - frc * t11_xy_xxxz[k] + f2t * (t20_y_xxxz[k] - t21_y_xxxz[k] + 3.0 * tk0_xy_xxz[k] - 3.0 * tk1_xy_xxz[k]);

                    t_xxy_xxyy[k] = fra * t10_xy_xxyy[k] - frc * t11_xy_xxyy[k] + f2t * (t20_y_xxyy[k] - t21_y_xxyy[k] + 2.0 * tk0_xy_xyy[k] - 2.0 * tk1_xy_xyy[k]);

                    t_xxy_xxyz[k] = fra * t10_xy_xxyz[k] - frc * t11_xy_xxyz[k] + f2t * (t20_y_xxyz[k] - t21_y_xxyz[k] + 2.0 * tk0_xy_xyz[k] - 2.0 * tk1_xy_xyz[k]);

                    t_xxy_xxzz[k] = fra * t10_xy_xxzz[k] - frc * t11_xy_xxzz[k] + f2t * (t20_y_xxzz[k] - t21_y_xxzz[k] + 2.0 * tk0_xy_xzz[k] - 2.0 * tk1_xy_xzz[k]);

                    t_xxy_xyyy[k] = fra * t10_xy_xyyy[k] - frc * t11_xy_xyyy[k] + f2t * (t20_y_xyyy[k] - t21_y_xyyy[k] + tk0_xy_yyy[k] - tk1_xy_yyy[k]);

                    t_xxy_xyyz[k] = fra * t10_xy_xyyz[k] - frc * t11_xy_xyyz[k] + f2t * (t20_y_xyyz[k] - t21_y_xyyz[k] + tk0_xy_yyz[k] - tk1_xy_yyz[k]);

                    t_xxy_xyzz[k] = fra * t10_xy_xyzz[k] - frc * t11_xy_xyzz[k] + f2t * (t20_y_xyzz[k] - t21_y_xyzz[k] + tk0_xy_yzz[k] - tk1_xy_yzz[k]);

                    t_xxy_xzzz[k] = fra * t10_xy_xzzz[k] - frc * t11_xy_xzzz[k] + f2t * (t20_y_xzzz[k] - t21_y_xzzz[k] + tk0_xy_zzz[k] - tk1_xy_zzz[k]);

                    t_xxy_yyyy[k] = fra * t10_xy_yyyy[k] - frc * t11_xy_yyyy[k] + f2t * (t20_y_yyyy[k] - t21_y_yyyy[k]);

                    t_xxy_yyyz[k] = fra * t10_xy_yyyz[k] - frc * t11_xy_yyyz[k] + f2t * (t20_y_yyyz[k] - t21_y_yyyz[k]);

                    t_xxy_yyzz[k] = fra * t10_xy_yyzz[k] - frc * t11_xy_yyzz[k] + f2t * (t20_y_yyzz[k] - t21_y_yyzz[k]);

                    t_xxy_yzzz[k] = fra * t10_xy_yzzz[k] - frc * t11_xy_yzzz[k] + f2t * (t20_y_yzzz[k] - t21_y_yzzz[k]);

                    t_xxy_zzzz[k] = fra * t10_xy_zzzz[k] - frc * t11_xy_zzzz[k] + f2t * (t20_y_zzzz[k] - t21_y_zzzz[k]);

                    t_xxz_xxxx[k] = fra * t10_xz_xxxx[k] - frc * t11_xz_xxxx[k] + f2t * (t20_z_xxxx[k] - t21_z_xxxx[k] + 4.0 * tk0_xz_xxx[k] - 4.0 * tk1_xz_xxx[k]);

                    t_xxz_xxxy[k] = fra * t10_xz_xxxy[k] - frc * t11_xz_xxxy[k] + f2t * (t20_z_xxxy[k] - t21_z_xxxy[k] + 3.0 * tk0_xz_xxy[k] - 3.0 * tk1_xz_xxy[k]);

                    t_xxz_xxxz[k] = fra * t10_xz_xxxz[k] - frc * t11_xz_xxxz[k] + f2t * (t20_z_xxxz[k] - t21_z_xxxz[k] + 3.0 * tk0_xz_xxz[k] - 3.0 * tk1_xz_xxz[k]);

                    t_xxz_xxyy[k] = fra * t10_xz_xxyy[k] - frc * t11_xz_xxyy[k] + f2t * (t20_z_xxyy[k] - t21_z_xxyy[k] + 2.0 * tk0_xz_xyy[k] - 2.0 * tk1_xz_xyy[k]);

                    t_xxz_xxyz[k] = fra * t10_xz_xxyz[k] - frc * t11_xz_xxyz[k] + f2t * (t20_z_xxyz[k] - t21_z_xxyz[k] + 2.0 * tk0_xz_xyz[k] - 2.0 * tk1_xz_xyz[k]);

                    t_xxz_xxzz[k] = fra * t10_xz_xxzz[k] - frc * t11_xz_xxzz[k] + f2t * (t20_z_xxzz[k] - t21_z_xxzz[k] + 2.0 * tk0_xz_xzz[k] - 2.0 * tk1_xz_xzz[k]);

                    t_xxz_xyyy[k] = fra * t10_xz_xyyy[k] - frc * t11_xz_xyyy[k] + f2t * (t20_z_xyyy[k] - t21_z_xyyy[k] + tk0_xz_yyy[k] - tk1_xz_yyy[k]);

                    t_xxz_xyyz[k] = fra * t10_xz_xyyz[k] - frc * t11_xz_xyyz[k] + f2t * (t20_z_xyyz[k] - t21_z_xyyz[k] + tk0_xz_yyz[k] - tk1_xz_yyz[k]);

                    t_xxz_xyzz[k] = fra * t10_xz_xyzz[k] - frc * t11_xz_xyzz[k] + f2t * (t20_z_xyzz[k] - t21_z_xyzz[k] + tk0_xz_yzz[k] - tk1_xz_yzz[k]);

                    t_xxz_xzzz[k] = fra * t10_xz_xzzz[k] - frc * t11_xz_xzzz[k] + f2t * (t20_z_xzzz[k] - t21_z_xzzz[k] + tk0_xz_zzz[k] - tk1_xz_zzz[k]);

                    t_xxz_yyyy[k] = fra * t10_xz_yyyy[k] - frc * t11_xz_yyyy[k] + f2t * (t20_z_yyyy[k] - t21_z_yyyy[k]);

                    t_xxz_yyyz[k] = fra * t10_xz_yyyz[k] - frc * t11_xz_yyyz[k] + f2t * (t20_z_yyyz[k] - t21_z_yyyz[k]);

                    t_xxz_yyzz[k] = fra * t10_xz_yyzz[k] - frc * t11_xz_yyzz[k] + f2t * (t20_z_yyzz[k] - t21_z_yyzz[k]);

                    t_xxz_yzzz[k] = fra * t10_xz_yzzz[k] - frc * t11_xz_yzzz[k] + f2t * (t20_z_yzzz[k] - t21_z_yzzz[k]);

                    t_xxz_zzzz[k] = fra * t10_xz_zzzz[k] - frc * t11_xz_zzzz[k] + f2t * (t20_z_zzzz[k] - t21_z_zzzz[k]);

                    t_xyy_xxxx[k] = fra * t10_yy_xxxx[k] - frc * t11_yy_xxxx[k] + f2t * (4.0 * tk0_yy_xxx[k] - 4.0 * tk1_yy_xxx[k]);

                    t_xyy_xxxy[k] = fra * t10_yy_xxxy[k] - frc * t11_yy_xxxy[k] + f2t * (3.0 * tk0_yy_xxy[k] - 3.0 * tk1_yy_xxy[k]);

                    t_xyy_xxxz[k] = fra * t10_yy_xxxz[k] - frc * t11_yy_xxxz[k] + f2t * (3.0 * tk0_yy_xxz[k] - 3.0 * tk1_yy_xxz[k]);

                    t_xyy_xxyy[k] = fra * t10_yy_xxyy[k] - frc * t11_yy_xxyy[k] + f2t * (2.0 * tk0_yy_xyy[k] - 2.0 * tk1_yy_xyy[k]);

                    t_xyy_xxyz[k] = fra * t10_yy_xxyz[k] - frc * t11_yy_xxyz[k] + f2t * (2.0 * tk0_yy_xyz[k] - 2.0 * tk1_yy_xyz[k]);

                    t_xyy_xxzz[k] = fra * t10_yy_xxzz[k] - frc * t11_yy_xxzz[k] + f2t * (2.0 * tk0_yy_xzz[k] - 2.0 * tk1_yy_xzz[k]);

                    t_xyy_xyyy[k] = fra * t10_yy_xyyy[k] - frc * t11_yy_xyyy[k] + f2t * (tk0_yy_yyy[k] - tk1_yy_yyy[k]);

                    t_xyy_xyyz[k] = fra * t10_yy_xyyz[k] - frc * t11_yy_xyyz[k] + f2t * (tk0_yy_yyz[k] - tk1_yy_yyz[k]);

                    t_xyy_xyzz[k] = fra * t10_yy_xyzz[k] - frc * t11_yy_xyzz[k] + f2t * (tk0_yy_yzz[k] - tk1_yy_yzz[k]);

                    t_xyy_xzzz[k] = fra * t10_yy_xzzz[k] - frc * t11_yy_xzzz[k] + f2t * (tk0_yy_zzz[k] - tk1_yy_zzz[k]);

                    t_xyy_yyyy[k] = fra * t10_yy_yyyy[k] - frc * t11_yy_yyyy[k];

                    t_xyy_yyyz[k] = fra * t10_yy_yyyz[k] - frc * t11_yy_yyyz[k];

                    t_xyy_yyzz[k] = fra * t10_yy_yyzz[k] - frc * t11_yy_yyzz[k];

                    t_xyy_yzzz[k] = fra * t10_yy_yzzz[k] - frc * t11_yy_yzzz[k];

                    t_xyy_zzzz[k] = fra * t10_yy_zzzz[k] - frc * t11_yy_zzzz[k];

                    t_xyz_xxxx[k] = fra * t10_yz_xxxx[k] - frc * t11_yz_xxxx[k] + f2t * (4.0 * tk0_yz_xxx[k] - 4.0 * tk1_yz_xxx[k]);

                    t_xyz_xxxy[k] = fra * t10_yz_xxxy[k] - frc * t11_yz_xxxy[k] + f2t * (3.0 * tk0_yz_xxy[k] - 3.0 * tk1_yz_xxy[k]);

                    t_xyz_xxxz[k] = fra * t10_yz_xxxz[k] - frc * t11_yz_xxxz[k] + f2t * (3.0 * tk0_yz_xxz[k] - 3.0 * tk1_yz_xxz[k]);

                    t_xyz_xxyy[k] = fra * t10_yz_xxyy[k] - frc * t11_yz_xxyy[k] + f2t * (2.0 * tk0_yz_xyy[k] - 2.0 * tk1_yz_xyy[k]);

                    t_xyz_xxyz[k] = fra * t10_yz_xxyz[k] - frc * t11_yz_xxyz[k] + f2t * (2.0 * tk0_yz_xyz[k] - 2.0 * tk1_yz_xyz[k]);

                    t_xyz_xxzz[k] = fra * t10_yz_xxzz[k] - frc * t11_yz_xxzz[k] + f2t * (2.0 * tk0_yz_xzz[k] - 2.0 * tk1_yz_xzz[k]);

                    t_xyz_xyyy[k] = fra * t10_yz_xyyy[k] - frc * t11_yz_xyyy[k] + f2t * (tk0_yz_yyy[k] - tk1_yz_yyy[k]);

                    t_xyz_xyyz[k] = fra * t10_yz_xyyz[k] - frc * t11_yz_xyyz[k] + f2t * (tk0_yz_yyz[k] - tk1_yz_yyz[k]);

                    t_xyz_xyzz[k] = fra * t10_yz_xyzz[k] - frc * t11_yz_xyzz[k] + f2t * (tk0_yz_yzz[k] - tk1_yz_yzz[k]);

                    t_xyz_xzzz[k] = fra * t10_yz_xzzz[k] - frc * t11_yz_xzzz[k] + f2t * (tk0_yz_zzz[k] - tk1_yz_zzz[k]);

                    t_xyz_yyyy[k] = fra * t10_yz_yyyy[k] - frc * t11_yz_yyyy[k];

                    t_xyz_yyyz[k] = fra * t10_yz_yyyz[k] - frc * t11_yz_yyyz[k];

                    t_xyz_yyzz[k] = fra * t10_yz_yyzz[k] - frc * t11_yz_yyzz[k];

                    t_xyz_yzzz[k] = fra * t10_yz_yzzz[k] - frc * t11_yz_yzzz[k];

                    t_xyz_zzzz[k] = fra * t10_yz_zzzz[k] - frc * t11_yz_zzzz[k];

                    t_xzz_xxxx[k] = fra * t10_zz_xxxx[k] - frc * t11_zz_xxxx[k] + f2t * (4.0 * tk0_zz_xxx[k] - 4.0 * tk1_zz_xxx[k]);

                    t_xzz_xxxy[k] = fra * t10_zz_xxxy[k] - frc * t11_zz_xxxy[k] + f2t * (3.0 * tk0_zz_xxy[k] - 3.0 * tk1_zz_xxy[k]);

                    t_xzz_xxxz[k] = fra * t10_zz_xxxz[k] - frc * t11_zz_xxxz[k] + f2t * (3.0 * tk0_zz_xxz[k] - 3.0 * tk1_zz_xxz[k]);

                    t_xzz_xxyy[k] = fra * t10_zz_xxyy[k] - frc * t11_zz_xxyy[k] + f2t * (2.0 * tk0_zz_xyy[k] - 2.0 * tk1_zz_xyy[k]);

                    t_xzz_xxyz[k] = fra * t10_zz_xxyz[k] - frc * t11_zz_xxyz[k] + f2t * (2.0 * tk0_zz_xyz[k] - 2.0 * tk1_zz_xyz[k]);

                    t_xzz_xxzz[k] = fra * t10_zz_xxzz[k] - frc * t11_zz_xxzz[k] + f2t * (2.0 * tk0_zz_xzz[k] - 2.0 * tk1_zz_xzz[k]);

                    t_xzz_xyyy[k] = fra * t10_zz_xyyy[k] - frc * t11_zz_xyyy[k] + f2t * (tk0_zz_yyy[k] - tk1_zz_yyy[k]);

                    t_xzz_xyyz[k] = fra * t10_zz_xyyz[k] - frc * t11_zz_xyyz[k] + f2t * (tk0_zz_yyz[k] - tk1_zz_yyz[k]);

                    t_xzz_xyzz[k] = fra * t10_zz_xyzz[k] - frc * t11_zz_xyzz[k] + f2t * (tk0_zz_yzz[k] - tk1_zz_yzz[k]);

                    t_xzz_xzzz[k] = fra * t10_zz_xzzz[k] - frc * t11_zz_xzzz[k] + f2t * (tk0_zz_zzz[k] - tk1_zz_zzz[k]);

                    t_xzz_yyyy[k] = fra * t10_zz_yyyy[k] - frc * t11_zz_yyyy[k];

                    t_xzz_yyyz[k] = fra * t10_zz_yyyz[k] - frc * t11_zz_yyyz[k];

                    t_xzz_yyzz[k] = fra * t10_zz_yyzz[k] - frc * t11_zz_yyzz[k];

                    t_xzz_yzzz[k] = fra * t10_zz_yzzz[k] - frc * t11_zz_yzzz[k];

                    t_xzz_zzzz[k] = fra * t10_zz_zzzz[k] - frc * t11_zz_zzzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyy_xxxx[k] = fra * t10_yy_xxxx[k] - frc * t11_yy_xxxx[k] + f2t * (2.0 * t20_y_xxxx[k] - 2.0 * t21_y_xxxx[k]);

                    t_yyy_xxxy[k] = fra * t10_yy_xxxy[k] - frc * t11_yy_xxxy[k] + f2t * (2.0 * t20_y_xxxy[k] - 2.0 * t21_y_xxxy[k] + tk0_yy_xxx[k] - tk1_yy_xxx[k]);

                    t_yyy_xxxz[k] = fra * t10_yy_xxxz[k] - frc * t11_yy_xxxz[k] + f2t * (2.0 * t20_y_xxxz[k] - 2.0 * t21_y_xxxz[k]);

                    t_yyy_xxyy[k] = fra * t10_yy_xxyy[k] - frc * t11_yy_xxyy[k] + f2t * (2.0 * t20_y_xxyy[k] - 2.0 * t21_y_xxyy[k] + 2.0 * tk0_yy_xxy[k] - 2.0 * tk1_yy_xxy[k]);

                    t_yyy_xxyz[k] = fra * t10_yy_xxyz[k] - frc * t11_yy_xxyz[k] + f2t * (2.0 * t20_y_xxyz[k] - 2.0 * t21_y_xxyz[k] + tk0_yy_xxz[k] - tk1_yy_xxz[k]);

                    t_yyy_xxzz[k] = fra * t10_yy_xxzz[k] - frc * t11_yy_xxzz[k] + f2t * (2.0 * t20_y_xxzz[k] - 2.0 * t21_y_xxzz[k]);

                    t_yyy_xyyy[k] = fra * t10_yy_xyyy[k] - frc * t11_yy_xyyy[k] + f2t * (2.0 * t20_y_xyyy[k] - 2.0 * t21_y_xyyy[k] + 3.0 * tk0_yy_xyy[k] - 3.0 * tk1_yy_xyy[k]);

                    t_yyy_xyyz[k] = fra * t10_yy_xyyz[k] - frc * t11_yy_xyyz[k] + f2t * (2.0 * t20_y_xyyz[k] - 2.0 * t21_y_xyyz[k] + 2.0 * tk0_yy_xyz[k] - 2.0 * tk1_yy_xyz[k]);

                    t_yyy_xyzz[k] = fra * t10_yy_xyzz[k] - frc * t11_yy_xyzz[k] + f2t * (2.0 * t20_y_xyzz[k] - 2.0 * t21_y_xyzz[k] + tk0_yy_xzz[k] - tk1_yy_xzz[k]);

                    t_yyy_xzzz[k] = fra * t10_yy_xzzz[k] - frc * t11_yy_xzzz[k] + f2t * (2.0 * t20_y_xzzz[k] - 2.0 * t21_y_xzzz[k]);

                    t_yyy_yyyy[k] = fra * t10_yy_yyyy[k] - frc * t11_yy_yyyy[k] + f2t * (2.0 * t20_y_yyyy[k] - 2.0 * t21_y_yyyy[k] + 4.0 * tk0_yy_yyy[k] - 4.0 * tk1_yy_yyy[k]);

                    t_yyy_yyyz[k] = fra * t10_yy_yyyz[k] - frc * t11_yy_yyyz[k] + f2t * (2.0 * t20_y_yyyz[k] - 2.0 * t21_y_yyyz[k] + 3.0 * tk0_yy_yyz[k] - 3.0 * tk1_yy_yyz[k]);

                    t_yyy_yyzz[k] = fra * t10_yy_yyzz[k] - frc * t11_yy_yyzz[k] + f2t * (2.0 * t20_y_yyzz[k] - 2.0 * t21_y_yyzz[k] + 2.0 * tk0_yy_yzz[k] - 2.0 * tk1_yy_yzz[k]);

                    t_yyy_yzzz[k] = fra * t10_yy_yzzz[k] - frc * t11_yy_yzzz[k] + f2t * (2.0 * t20_y_yzzz[k] - 2.0 * t21_y_yzzz[k] + tk0_yy_zzz[k] - tk1_yy_zzz[k]);

                    t_yyy_zzzz[k] = fra * t10_yy_zzzz[k] - frc * t11_yy_zzzz[k] + f2t * (2.0 * t20_y_zzzz[k] - 2.0 * t21_y_zzzz[k]);

                    t_yyz_xxxx[k] = fra * t10_yz_xxxx[k] - frc * t11_yz_xxxx[k] + f2t * (t20_z_xxxx[k] - t21_z_xxxx[k]);

                    t_yyz_xxxy[k] = fra * t10_yz_xxxy[k] - frc * t11_yz_xxxy[k] + f2t * (t20_z_xxxy[k] - t21_z_xxxy[k] + tk0_yz_xxx[k] - tk1_yz_xxx[k]);

                    t_yyz_xxxz[k] = fra * t10_yz_xxxz[k] - frc * t11_yz_xxxz[k] + f2t * (t20_z_xxxz[k] - t21_z_xxxz[k]);

                    t_yyz_xxyy[k] = fra * t10_yz_xxyy[k] - frc * t11_yz_xxyy[k] + f2t * (t20_z_xxyy[k] - t21_z_xxyy[k] + 2.0 * tk0_yz_xxy[k] - 2.0 * tk1_yz_xxy[k]);

                    t_yyz_xxyz[k] = fra * t10_yz_xxyz[k] - frc * t11_yz_xxyz[k] + f2t * (t20_z_xxyz[k] - t21_z_xxyz[k] + tk0_yz_xxz[k] - tk1_yz_xxz[k]);

                    t_yyz_xxzz[k] = fra * t10_yz_xxzz[k] - frc * t11_yz_xxzz[k] + f2t * (t20_z_xxzz[k] - t21_z_xxzz[k]);

                    t_yyz_xyyy[k] = fra * t10_yz_xyyy[k] - frc * t11_yz_xyyy[k] + f2t * (t20_z_xyyy[k] - t21_z_xyyy[k] + 3.0 * tk0_yz_xyy[k] - 3.0 * tk1_yz_xyy[k]);

                    t_yyz_xyyz[k] = fra * t10_yz_xyyz[k] - frc * t11_yz_xyyz[k] + f2t * (t20_z_xyyz[k] - t21_z_xyyz[k] + 2.0 * tk0_yz_xyz[k] - 2.0 * tk1_yz_xyz[k]);

                    t_yyz_xyzz[k] = fra * t10_yz_xyzz[k] - frc * t11_yz_xyzz[k] + f2t * (t20_z_xyzz[k] - t21_z_xyzz[k] + tk0_yz_xzz[k] - tk1_yz_xzz[k]);

                    t_yyz_xzzz[k] = fra * t10_yz_xzzz[k] - frc * t11_yz_xzzz[k] + f2t * (t20_z_xzzz[k] - t21_z_xzzz[k]);

                    t_yyz_yyyy[k] = fra * t10_yz_yyyy[k] - frc * t11_yz_yyyy[k] + f2t * (t20_z_yyyy[k] - t21_z_yyyy[k] + 4.0 * tk0_yz_yyy[k] - 4.0 * tk1_yz_yyy[k]);

                    t_yyz_yyyz[k] = fra * t10_yz_yyyz[k] - frc * t11_yz_yyyz[k] + f2t * (t20_z_yyyz[k] - t21_z_yyyz[k] + 3.0 * tk0_yz_yyz[k] - 3.0 * tk1_yz_yyz[k]);

                    t_yyz_yyzz[k] = fra * t10_yz_yyzz[k] - frc * t11_yz_yyzz[k] + f2t * (t20_z_yyzz[k] - t21_z_yyzz[k] + 2.0 * tk0_yz_yzz[k] - 2.0 * tk1_yz_yzz[k]);

                    t_yyz_yzzz[k] = fra * t10_yz_yzzz[k] - frc * t11_yz_yzzz[k] + f2t * (t20_z_yzzz[k] - t21_z_yzzz[k] + tk0_yz_zzz[k] - tk1_yz_zzz[k]);

                    t_yyz_zzzz[k] = fra * t10_yz_zzzz[k] - frc * t11_yz_zzzz[k] + f2t * (t20_z_zzzz[k] - t21_z_zzzz[k]);

                    t_yzz_xxxx[k] = fra * t10_zz_xxxx[k] - frc * t11_zz_xxxx[k];

                    t_yzz_xxxy[k] = fra * t10_zz_xxxy[k] - frc * t11_zz_xxxy[k] + f2t * (tk0_zz_xxx[k] - tk1_zz_xxx[k]);

                    t_yzz_xxxz[k] = fra * t10_zz_xxxz[k] - frc * t11_zz_xxxz[k];

                    t_yzz_xxyy[k] = fra * t10_zz_xxyy[k] - frc * t11_zz_xxyy[k] + f2t * (2.0 * tk0_zz_xxy[k] - 2.0 * tk1_zz_xxy[k]);

                    t_yzz_xxyz[k] = fra * t10_zz_xxyz[k] - frc * t11_zz_xxyz[k] + f2t * (tk0_zz_xxz[k] - tk1_zz_xxz[k]);

                    t_yzz_xxzz[k] = fra * t10_zz_xxzz[k] - frc * t11_zz_xxzz[k];

                    t_yzz_xyyy[k] = fra * t10_zz_xyyy[k] - frc * t11_zz_xyyy[k] + f2t * (3.0 * tk0_zz_xyy[k] - 3.0 * tk1_zz_xyy[k]);

                    t_yzz_xyyz[k] = fra * t10_zz_xyyz[k] - frc * t11_zz_xyyz[k] + f2t * (2.0 * tk0_zz_xyz[k] - 2.0 * tk1_zz_xyz[k]);

                    t_yzz_xyzz[k] = fra * t10_zz_xyzz[k] - frc * t11_zz_xyzz[k] + f2t * (tk0_zz_xzz[k] - tk1_zz_xzz[k]);

                    t_yzz_xzzz[k] = fra * t10_zz_xzzz[k] - frc * t11_zz_xzzz[k];

                    t_yzz_yyyy[k] = fra * t10_zz_yyyy[k] - frc * t11_zz_yyyy[k] + f2t * (4.0 * tk0_zz_yyy[k] - 4.0 * tk1_zz_yyy[k]);

                    t_yzz_yyyz[k] = fra * t10_zz_yyyz[k] - frc * t11_zz_yyyz[k] + f2t * (3.0 * tk0_zz_yyz[k] - 3.0 * tk1_zz_yyz[k]);

                    t_yzz_yyzz[k] = fra * t10_zz_yyzz[k] - frc * t11_zz_yyzz[k] + f2t * (2.0 * tk0_zz_yzz[k] - 2.0 * tk1_zz_yzz[k]);

                    t_yzz_yzzz[k] = fra * t10_zz_yzzz[k] - frc * t11_zz_yzzz[k] + f2t * (tk0_zz_zzz[k] - tk1_zz_zzz[k]);

                    t_yzz_zzzz[k] = fra * t10_zz_zzzz[k] - frc * t11_zz_zzzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzz_xxxx[k] = fra * t10_zz_xxxx[k] - frc * t11_zz_xxxx[k] + f2t * (2.0 * t20_z_xxxx[k] - 2.0 * t21_z_xxxx[k]);

                    t_zzz_xxxy[k] = fra * t10_zz_xxxy[k] - frc * t11_zz_xxxy[k] + f2t * (2.0 * t20_z_xxxy[k] - 2.0 * t21_z_xxxy[k]);

                    t_zzz_xxxz[k] = fra * t10_zz_xxxz[k] - frc * t11_zz_xxxz[k] + f2t * (2.0 * t20_z_xxxz[k] - 2.0 * t21_z_xxxz[k] + tk0_zz_xxx[k] - tk1_zz_xxx[k]);

                    t_zzz_xxyy[k] = fra * t10_zz_xxyy[k] - frc * t11_zz_xxyy[k] + f2t * (2.0 * t20_z_xxyy[k] - 2.0 * t21_z_xxyy[k]);

                    t_zzz_xxyz[k] = fra * t10_zz_xxyz[k] - frc * t11_zz_xxyz[k] + f2t * (2.0 * t20_z_xxyz[k] - 2.0 * t21_z_xxyz[k] + tk0_zz_xxy[k] - tk1_zz_xxy[k]);

                    t_zzz_xxzz[k] = fra * t10_zz_xxzz[k] - frc * t11_zz_xxzz[k] + f2t * (2.0 * t20_z_xxzz[k] - 2.0 * t21_z_xxzz[k] + 2.0 * tk0_zz_xxz[k] - 2.0 * tk1_zz_xxz[k]);

                    t_zzz_xyyy[k] = fra * t10_zz_xyyy[k] - frc * t11_zz_xyyy[k] + f2t * (2.0 * t20_z_xyyy[k] - 2.0 * t21_z_xyyy[k]);

                    t_zzz_xyyz[k] = fra * t10_zz_xyyz[k] - frc * t11_zz_xyyz[k] + f2t * (2.0 * t20_z_xyyz[k] - 2.0 * t21_z_xyyz[k] + tk0_zz_xyy[k] - tk1_zz_xyy[k]);

                    t_zzz_xyzz[k] = fra * t10_zz_xyzz[k] - frc * t11_zz_xyzz[k] + f2t * (2.0 * t20_z_xyzz[k] - 2.0 * t21_z_xyzz[k] + 2.0 * tk0_zz_xyz[k] - 2.0 * tk1_zz_xyz[k]);

                    t_zzz_xzzz[k] = fra * t10_zz_xzzz[k] - frc * t11_zz_xzzz[k] + f2t * (2.0 * t20_z_xzzz[k] - 2.0 * t21_z_xzzz[k] + 3.0 * tk0_zz_xzz[k] - 3.0 * tk1_zz_xzz[k]);

                    t_zzz_yyyy[k] = fra * t10_zz_yyyy[k] - frc * t11_zz_yyyy[k] + f2t * (2.0 * t20_z_yyyy[k] - 2.0 * t21_z_yyyy[k]);

                    t_zzz_yyyz[k] = fra * t10_zz_yyyz[k] - frc * t11_zz_yyyz[k] + f2t * (2.0 * t20_z_yyyz[k] - 2.0 * t21_z_yyyz[k] + tk0_zz_yyy[k] - tk1_zz_yyy[k]);

                    t_zzz_yyzz[k] = fra * t10_zz_yyzz[k] - frc * t11_zz_yyzz[k] + f2t * (2.0 * t20_z_yyzz[k] - 2.0 * t21_z_yyzz[k] + 2.0 * tk0_zz_yyz[k] - 2.0 * tk1_zz_yyz[k]);

                    t_zzz_yzzz[k] = fra * t10_zz_yzzz[k] - frc * t11_zz_yzzz[k] + f2t * (2.0 * t20_z_yzzz[k] - 2.0 * t21_z_yzzz[k] + 3.0 * tk0_zz_yzz[k] - 3.0 * tk1_zz_yzz[k]);

                    t_zzz_zzzz[k] = fra * t10_zz_zzzz[k] - frc * t11_zz_zzzz[k] + f2t * (2.0 * t20_z_zzzz[k] - 2.0 * t21_z_zzzz[k] + 4.0 * tk0_zz_zzz[k] - 4.0 * tk1_zz_zzz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForGF(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {4, 3, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 4, 3);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {4, 3, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 3, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 3, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 3, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 2, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 2, i + 1});

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

                // set up pointers to (F|A(0)|D)^(m) integrals

                auto tk0_xxx_xx = primBuffer.data(tk0off + 60 * idx);

                auto tk0_xxx_xy = primBuffer.data(tk0off + 60 * idx + 1);

                auto tk0_xxx_xz = primBuffer.data(tk0off + 60 * idx + 2);

                auto tk0_xxx_yy = primBuffer.data(tk0off + 60 * idx + 3);

                auto tk0_xxx_yz = primBuffer.data(tk0off + 60 * idx + 4);

                auto tk0_xxx_zz = primBuffer.data(tk0off + 60 * idx + 5);

                auto tk0_xxy_xx = primBuffer.data(tk0off + 60 * idx + 6);

                auto tk0_xxy_xy = primBuffer.data(tk0off + 60 * idx + 7);

                auto tk0_xxy_xz = primBuffer.data(tk0off + 60 * idx + 8);

                auto tk0_xxy_yy = primBuffer.data(tk0off + 60 * idx + 9);

                auto tk0_xxy_yz = primBuffer.data(tk0off + 60 * idx + 10);

                auto tk0_xxy_zz = primBuffer.data(tk0off + 60 * idx + 11);

                auto tk0_xxz_xx = primBuffer.data(tk0off + 60 * idx + 12);

                auto tk0_xxz_xy = primBuffer.data(tk0off + 60 * idx + 13);

                auto tk0_xxz_xz = primBuffer.data(tk0off + 60 * idx + 14);

                auto tk0_xxz_yy = primBuffer.data(tk0off + 60 * idx + 15);

                auto tk0_xxz_yz = primBuffer.data(tk0off + 60 * idx + 16);

                auto tk0_xxz_zz = primBuffer.data(tk0off + 60 * idx + 17);

                auto tk0_xyy_xx = primBuffer.data(tk0off + 60 * idx + 18);

                auto tk0_xyy_xy = primBuffer.data(tk0off + 60 * idx + 19);

                auto tk0_xyy_xz = primBuffer.data(tk0off + 60 * idx + 20);

                auto tk0_xyy_yy = primBuffer.data(tk0off + 60 * idx + 21);

                auto tk0_xyy_yz = primBuffer.data(tk0off + 60 * idx + 22);

                auto tk0_xyy_zz = primBuffer.data(tk0off + 60 * idx + 23);

                auto tk0_xyz_xx = primBuffer.data(tk0off + 60 * idx + 24);

                auto tk0_xyz_xy = primBuffer.data(tk0off + 60 * idx + 25);

                auto tk0_xyz_xz = primBuffer.data(tk0off + 60 * idx + 26);

                auto tk0_xyz_yy = primBuffer.data(tk0off + 60 * idx + 27);

                auto tk0_xyz_yz = primBuffer.data(tk0off + 60 * idx + 28);

                auto tk0_xyz_zz = primBuffer.data(tk0off + 60 * idx + 29);

                auto tk0_xzz_xx = primBuffer.data(tk0off + 60 * idx + 30);

                auto tk0_xzz_xy = primBuffer.data(tk0off + 60 * idx + 31);

                auto tk0_xzz_xz = primBuffer.data(tk0off + 60 * idx + 32);

                auto tk0_xzz_yy = primBuffer.data(tk0off + 60 * idx + 33);

                auto tk0_xzz_yz = primBuffer.data(tk0off + 60 * idx + 34);

                auto tk0_xzz_zz = primBuffer.data(tk0off + 60 * idx + 35);

                auto tk0_yyy_xx = primBuffer.data(tk0off + 60 * idx + 36);

                auto tk0_yyy_xy = primBuffer.data(tk0off + 60 * idx + 37);

                auto tk0_yyy_xz = primBuffer.data(tk0off + 60 * idx + 38);

                auto tk0_yyy_yy = primBuffer.data(tk0off + 60 * idx + 39);

                auto tk0_yyy_yz = primBuffer.data(tk0off + 60 * idx + 40);

                auto tk0_yyy_zz = primBuffer.data(tk0off + 60 * idx + 41);

                auto tk0_yyz_xx = primBuffer.data(tk0off + 60 * idx + 42);

                auto tk0_yyz_xy = primBuffer.data(tk0off + 60 * idx + 43);

                auto tk0_yyz_xz = primBuffer.data(tk0off + 60 * idx + 44);

                auto tk0_yyz_yy = primBuffer.data(tk0off + 60 * idx + 45);

                auto tk0_yyz_yz = primBuffer.data(tk0off + 60 * idx + 46);

                auto tk0_yyz_zz = primBuffer.data(tk0off + 60 * idx + 47);

                auto tk0_yzz_xx = primBuffer.data(tk0off + 60 * idx + 48);

                auto tk0_yzz_xy = primBuffer.data(tk0off + 60 * idx + 49);

                auto tk0_yzz_xz = primBuffer.data(tk0off + 60 * idx + 50);

                auto tk0_yzz_yy = primBuffer.data(tk0off + 60 * idx + 51);

                auto tk0_yzz_yz = primBuffer.data(tk0off + 60 * idx + 52);

                auto tk0_yzz_zz = primBuffer.data(tk0off + 60 * idx + 53);

                auto tk0_zzz_xx = primBuffer.data(tk0off + 60 * idx + 54);

                auto tk0_zzz_xy = primBuffer.data(tk0off + 60 * idx + 55);

                auto tk0_zzz_xz = primBuffer.data(tk0off + 60 * idx + 56);

                auto tk0_zzz_yy = primBuffer.data(tk0off + 60 * idx + 57);

                auto tk0_zzz_yz = primBuffer.data(tk0off + 60 * idx + 58);

                auto tk0_zzz_zz = primBuffer.data(tk0off + 60 * idx + 59);

                // set up pointers to (F|A(0)|D)^(m+1) integrals

                auto tk1_xxx_xx = primBuffer.data(tk1off + 60 * idx);

                auto tk1_xxx_xy = primBuffer.data(tk1off + 60 * idx + 1);

                auto tk1_xxx_xz = primBuffer.data(tk1off + 60 * idx + 2);

                auto tk1_xxx_yy = primBuffer.data(tk1off + 60 * idx + 3);

                auto tk1_xxx_yz = primBuffer.data(tk1off + 60 * idx + 4);

                auto tk1_xxx_zz = primBuffer.data(tk1off + 60 * idx + 5);

                auto tk1_xxy_xx = primBuffer.data(tk1off + 60 * idx + 6);

                auto tk1_xxy_xy = primBuffer.data(tk1off + 60 * idx + 7);

                auto tk1_xxy_xz = primBuffer.data(tk1off + 60 * idx + 8);

                auto tk1_xxy_yy = primBuffer.data(tk1off + 60 * idx + 9);

                auto tk1_xxy_yz = primBuffer.data(tk1off + 60 * idx + 10);

                auto tk1_xxy_zz = primBuffer.data(tk1off + 60 * idx + 11);

                auto tk1_xxz_xx = primBuffer.data(tk1off + 60 * idx + 12);

                auto tk1_xxz_xy = primBuffer.data(tk1off + 60 * idx + 13);

                auto tk1_xxz_xz = primBuffer.data(tk1off + 60 * idx + 14);

                auto tk1_xxz_yy = primBuffer.data(tk1off + 60 * idx + 15);

                auto tk1_xxz_yz = primBuffer.data(tk1off + 60 * idx + 16);

                auto tk1_xxz_zz = primBuffer.data(tk1off + 60 * idx + 17);

                auto tk1_xyy_xx = primBuffer.data(tk1off + 60 * idx + 18);

                auto tk1_xyy_xy = primBuffer.data(tk1off + 60 * idx + 19);

                auto tk1_xyy_xz = primBuffer.data(tk1off + 60 * idx + 20);

                auto tk1_xyy_yy = primBuffer.data(tk1off + 60 * idx + 21);

                auto tk1_xyy_yz = primBuffer.data(tk1off + 60 * idx + 22);

                auto tk1_xyy_zz = primBuffer.data(tk1off + 60 * idx + 23);

                auto tk1_xyz_xx = primBuffer.data(tk1off + 60 * idx + 24);

                auto tk1_xyz_xy = primBuffer.data(tk1off + 60 * idx + 25);

                auto tk1_xyz_xz = primBuffer.data(tk1off + 60 * idx + 26);

                auto tk1_xyz_yy = primBuffer.data(tk1off + 60 * idx + 27);

                auto tk1_xyz_yz = primBuffer.data(tk1off + 60 * idx + 28);

                auto tk1_xyz_zz = primBuffer.data(tk1off + 60 * idx + 29);

                auto tk1_xzz_xx = primBuffer.data(tk1off + 60 * idx + 30);

                auto tk1_xzz_xy = primBuffer.data(tk1off + 60 * idx + 31);

                auto tk1_xzz_xz = primBuffer.data(tk1off + 60 * idx + 32);

                auto tk1_xzz_yy = primBuffer.data(tk1off + 60 * idx + 33);

                auto tk1_xzz_yz = primBuffer.data(tk1off + 60 * idx + 34);

                auto tk1_xzz_zz = primBuffer.data(tk1off + 60 * idx + 35);

                auto tk1_yyy_xx = primBuffer.data(tk1off + 60 * idx + 36);

                auto tk1_yyy_xy = primBuffer.data(tk1off + 60 * idx + 37);

                auto tk1_yyy_xz = primBuffer.data(tk1off + 60 * idx + 38);

                auto tk1_yyy_yy = primBuffer.data(tk1off + 60 * idx + 39);

                auto tk1_yyy_yz = primBuffer.data(tk1off + 60 * idx + 40);

                auto tk1_yyy_zz = primBuffer.data(tk1off + 60 * idx + 41);

                auto tk1_yyz_xx = primBuffer.data(tk1off + 60 * idx + 42);

                auto tk1_yyz_xy = primBuffer.data(tk1off + 60 * idx + 43);

                auto tk1_yyz_xz = primBuffer.data(tk1off + 60 * idx + 44);

                auto tk1_yyz_yy = primBuffer.data(tk1off + 60 * idx + 45);

                auto tk1_yyz_yz = primBuffer.data(tk1off + 60 * idx + 46);

                auto tk1_yyz_zz = primBuffer.data(tk1off + 60 * idx + 47);

                auto tk1_yzz_xx = primBuffer.data(tk1off + 60 * idx + 48);

                auto tk1_yzz_xy = primBuffer.data(tk1off + 60 * idx + 49);

                auto tk1_yzz_xz = primBuffer.data(tk1off + 60 * idx + 50);

                auto tk1_yzz_yy = primBuffer.data(tk1off + 60 * idx + 51);

                auto tk1_yzz_yz = primBuffer.data(tk1off + 60 * idx + 52);

                auto tk1_yzz_zz = primBuffer.data(tk1off + 60 * idx + 53);

                auto tk1_zzz_xx = primBuffer.data(tk1off + 60 * idx + 54);

                auto tk1_zzz_xy = primBuffer.data(tk1off + 60 * idx + 55);

                auto tk1_zzz_xz = primBuffer.data(tk1off + 60 * idx + 56);

                auto tk1_zzz_yy = primBuffer.data(tk1off + 60 * idx + 57);

                auto tk1_zzz_yz = primBuffer.data(tk1off + 60 * idx + 58);

                auto tk1_zzz_zz = primBuffer.data(tk1off + 60 * idx + 59);

                // set up pointers to (D|A(0)|F)^(m) integrals

                auto t20_xx_xxx = primBuffer.data(t20off + 60 * idx);

                auto t20_xx_xxy = primBuffer.data(t20off + 60 * idx + 1);

                auto t20_xx_xxz = primBuffer.data(t20off + 60 * idx + 2);

                auto t20_xx_xyy = primBuffer.data(t20off + 60 * idx + 3);

                auto t20_xx_xyz = primBuffer.data(t20off + 60 * idx + 4);

                auto t20_xx_xzz = primBuffer.data(t20off + 60 * idx + 5);

                auto t20_xx_yyy = primBuffer.data(t20off + 60 * idx + 6);

                auto t20_xx_yyz = primBuffer.data(t20off + 60 * idx + 7);

                auto t20_xx_yzz = primBuffer.data(t20off + 60 * idx + 8);

                auto t20_xx_zzz = primBuffer.data(t20off + 60 * idx + 9);

                auto t20_xy_xxx = primBuffer.data(t20off + 60 * idx + 10);

                auto t20_xy_xxy = primBuffer.data(t20off + 60 * idx + 11);

                auto t20_xy_xxz = primBuffer.data(t20off + 60 * idx + 12);

                auto t20_xy_xyy = primBuffer.data(t20off + 60 * idx + 13);

                auto t20_xy_xyz = primBuffer.data(t20off + 60 * idx + 14);

                auto t20_xy_xzz = primBuffer.data(t20off + 60 * idx + 15);

                auto t20_xy_yyy = primBuffer.data(t20off + 60 * idx + 16);

                auto t20_xy_yyz = primBuffer.data(t20off + 60 * idx + 17);

                auto t20_xy_yzz = primBuffer.data(t20off + 60 * idx + 18);

                auto t20_xy_zzz = primBuffer.data(t20off + 60 * idx + 19);

                auto t20_xz_xxx = primBuffer.data(t20off + 60 * idx + 20);

                auto t20_xz_xxy = primBuffer.data(t20off + 60 * idx + 21);

                auto t20_xz_xxz = primBuffer.data(t20off + 60 * idx + 22);

                auto t20_xz_xyy = primBuffer.data(t20off + 60 * idx + 23);

                auto t20_xz_xyz = primBuffer.data(t20off + 60 * idx + 24);

                auto t20_xz_xzz = primBuffer.data(t20off + 60 * idx + 25);

                auto t20_xz_yyy = primBuffer.data(t20off + 60 * idx + 26);

                auto t20_xz_yyz = primBuffer.data(t20off + 60 * idx + 27);

                auto t20_xz_yzz = primBuffer.data(t20off + 60 * idx + 28);

                auto t20_xz_zzz = primBuffer.data(t20off + 60 * idx + 29);

                auto t20_yy_xxx = primBuffer.data(t20off + 60 * idx + 30);

                auto t20_yy_xxy = primBuffer.data(t20off + 60 * idx + 31);

                auto t20_yy_xxz = primBuffer.data(t20off + 60 * idx + 32);

                auto t20_yy_xyy = primBuffer.data(t20off + 60 * idx + 33);

                auto t20_yy_xyz = primBuffer.data(t20off + 60 * idx + 34);

                auto t20_yy_xzz = primBuffer.data(t20off + 60 * idx + 35);

                auto t20_yy_yyy = primBuffer.data(t20off + 60 * idx + 36);

                auto t20_yy_yyz = primBuffer.data(t20off + 60 * idx + 37);

                auto t20_yy_yzz = primBuffer.data(t20off + 60 * idx + 38);

                auto t20_yy_zzz = primBuffer.data(t20off + 60 * idx + 39);

                auto t20_yz_xxx = primBuffer.data(t20off + 60 * idx + 40);

                auto t20_yz_xxy = primBuffer.data(t20off + 60 * idx + 41);

                auto t20_yz_xxz = primBuffer.data(t20off + 60 * idx + 42);

                auto t20_yz_xyy = primBuffer.data(t20off + 60 * idx + 43);

                auto t20_yz_xyz = primBuffer.data(t20off + 60 * idx + 44);

                auto t20_yz_xzz = primBuffer.data(t20off + 60 * idx + 45);

                auto t20_yz_yyy = primBuffer.data(t20off + 60 * idx + 46);

                auto t20_yz_yyz = primBuffer.data(t20off + 60 * idx + 47);

                auto t20_yz_yzz = primBuffer.data(t20off + 60 * idx + 48);

                auto t20_yz_zzz = primBuffer.data(t20off + 60 * idx + 49);

                auto t20_zz_xxx = primBuffer.data(t20off + 60 * idx + 50);

                auto t20_zz_xxy = primBuffer.data(t20off + 60 * idx + 51);

                auto t20_zz_xxz = primBuffer.data(t20off + 60 * idx + 52);

                auto t20_zz_xyy = primBuffer.data(t20off + 60 * idx + 53);

                auto t20_zz_xyz = primBuffer.data(t20off + 60 * idx + 54);

                auto t20_zz_xzz = primBuffer.data(t20off + 60 * idx + 55);

                auto t20_zz_yyy = primBuffer.data(t20off + 60 * idx + 56);

                auto t20_zz_yyz = primBuffer.data(t20off + 60 * idx + 57);

                auto t20_zz_yzz = primBuffer.data(t20off + 60 * idx + 58);

                auto t20_zz_zzz = primBuffer.data(t20off + 60 * idx + 59);

                // set up pointers to (D|A(0)|F)^(m+1) integrals

                auto t21_xx_xxx = primBuffer.data(t21off + 60 * idx);

                auto t21_xx_xxy = primBuffer.data(t21off + 60 * idx + 1);

                auto t21_xx_xxz = primBuffer.data(t21off + 60 * idx + 2);

                auto t21_xx_xyy = primBuffer.data(t21off + 60 * idx + 3);

                auto t21_xx_xyz = primBuffer.data(t21off + 60 * idx + 4);

                auto t21_xx_xzz = primBuffer.data(t21off + 60 * idx + 5);

                auto t21_xx_yyy = primBuffer.data(t21off + 60 * idx + 6);

                auto t21_xx_yyz = primBuffer.data(t21off + 60 * idx + 7);

                auto t21_xx_yzz = primBuffer.data(t21off + 60 * idx + 8);

                auto t21_xx_zzz = primBuffer.data(t21off + 60 * idx + 9);

                auto t21_xy_xxx = primBuffer.data(t21off + 60 * idx + 10);

                auto t21_xy_xxy = primBuffer.data(t21off + 60 * idx + 11);

                auto t21_xy_xxz = primBuffer.data(t21off + 60 * idx + 12);

                auto t21_xy_xyy = primBuffer.data(t21off + 60 * idx + 13);

                auto t21_xy_xyz = primBuffer.data(t21off + 60 * idx + 14);

                auto t21_xy_xzz = primBuffer.data(t21off + 60 * idx + 15);

                auto t21_xy_yyy = primBuffer.data(t21off + 60 * idx + 16);

                auto t21_xy_yyz = primBuffer.data(t21off + 60 * idx + 17);

                auto t21_xy_yzz = primBuffer.data(t21off + 60 * idx + 18);

                auto t21_xy_zzz = primBuffer.data(t21off + 60 * idx + 19);

                auto t21_xz_xxx = primBuffer.data(t21off + 60 * idx + 20);

                auto t21_xz_xxy = primBuffer.data(t21off + 60 * idx + 21);

                auto t21_xz_xxz = primBuffer.data(t21off + 60 * idx + 22);

                auto t21_xz_xyy = primBuffer.data(t21off + 60 * idx + 23);

                auto t21_xz_xyz = primBuffer.data(t21off + 60 * idx + 24);

                auto t21_xz_xzz = primBuffer.data(t21off + 60 * idx + 25);

                auto t21_xz_yyy = primBuffer.data(t21off + 60 * idx + 26);

                auto t21_xz_yyz = primBuffer.data(t21off + 60 * idx + 27);

                auto t21_xz_yzz = primBuffer.data(t21off + 60 * idx + 28);

                auto t21_xz_zzz = primBuffer.data(t21off + 60 * idx + 29);

                auto t21_yy_xxx = primBuffer.data(t21off + 60 * idx + 30);

                auto t21_yy_xxy = primBuffer.data(t21off + 60 * idx + 31);

                auto t21_yy_xxz = primBuffer.data(t21off + 60 * idx + 32);

                auto t21_yy_xyy = primBuffer.data(t21off + 60 * idx + 33);

                auto t21_yy_xyz = primBuffer.data(t21off + 60 * idx + 34);

                auto t21_yy_xzz = primBuffer.data(t21off + 60 * idx + 35);

                auto t21_yy_yyy = primBuffer.data(t21off + 60 * idx + 36);

                auto t21_yy_yyz = primBuffer.data(t21off + 60 * idx + 37);

                auto t21_yy_yzz = primBuffer.data(t21off + 60 * idx + 38);

                auto t21_yy_zzz = primBuffer.data(t21off + 60 * idx + 39);

                auto t21_yz_xxx = primBuffer.data(t21off + 60 * idx + 40);

                auto t21_yz_xxy = primBuffer.data(t21off + 60 * idx + 41);

                auto t21_yz_xxz = primBuffer.data(t21off + 60 * idx + 42);

                auto t21_yz_xyy = primBuffer.data(t21off + 60 * idx + 43);

                auto t21_yz_xyz = primBuffer.data(t21off + 60 * idx + 44);

                auto t21_yz_xzz = primBuffer.data(t21off + 60 * idx + 45);

                auto t21_yz_yyy = primBuffer.data(t21off + 60 * idx + 46);

                auto t21_yz_yyz = primBuffer.data(t21off + 60 * idx + 47);

                auto t21_yz_yzz = primBuffer.data(t21off + 60 * idx + 48);

                auto t21_yz_zzz = primBuffer.data(t21off + 60 * idx + 49);

                auto t21_zz_xxx = primBuffer.data(t21off + 60 * idx + 50);

                auto t21_zz_xxy = primBuffer.data(t21off + 60 * idx + 51);

                auto t21_zz_xxz = primBuffer.data(t21off + 60 * idx + 52);

                auto t21_zz_xyy = primBuffer.data(t21off + 60 * idx + 53);

                auto t21_zz_xyz = primBuffer.data(t21off + 60 * idx + 54);

                auto t21_zz_xzz = primBuffer.data(t21off + 60 * idx + 55);

                auto t21_zz_yyy = primBuffer.data(t21off + 60 * idx + 56);

                auto t21_zz_yyz = primBuffer.data(t21off + 60 * idx + 57);

                auto t21_zz_yzz = primBuffer.data(t21off + 60 * idx + 58);

                auto t21_zz_zzz = primBuffer.data(t21off + 60 * idx + 59);

                // set up pointers to (F|A(0)|F)^(m) integrals

                auto t10_xxx_xxx = primBuffer.data(t10off + 100 * idx);

                auto t10_xxx_xxy = primBuffer.data(t10off + 100 * idx + 1);

                auto t10_xxx_xxz = primBuffer.data(t10off + 100 * idx + 2);

                auto t10_xxx_xyy = primBuffer.data(t10off + 100 * idx + 3);

                auto t10_xxx_xyz = primBuffer.data(t10off + 100 * idx + 4);

                auto t10_xxx_xzz = primBuffer.data(t10off + 100 * idx + 5);

                auto t10_xxx_yyy = primBuffer.data(t10off + 100 * idx + 6);

                auto t10_xxx_yyz = primBuffer.data(t10off + 100 * idx + 7);

                auto t10_xxx_yzz = primBuffer.data(t10off + 100 * idx + 8);

                auto t10_xxx_zzz = primBuffer.data(t10off + 100 * idx + 9);

                auto t10_xxy_xxx = primBuffer.data(t10off + 100 * idx + 10);

                auto t10_xxy_xxy = primBuffer.data(t10off + 100 * idx + 11);

                auto t10_xxy_xxz = primBuffer.data(t10off + 100 * idx + 12);

                auto t10_xxy_xyy = primBuffer.data(t10off + 100 * idx + 13);

                auto t10_xxy_xyz = primBuffer.data(t10off + 100 * idx + 14);

                auto t10_xxy_xzz = primBuffer.data(t10off + 100 * idx + 15);

                auto t10_xxy_yyy = primBuffer.data(t10off + 100 * idx + 16);

                auto t10_xxy_yyz = primBuffer.data(t10off + 100 * idx + 17);

                auto t10_xxy_yzz = primBuffer.data(t10off + 100 * idx + 18);

                auto t10_xxy_zzz = primBuffer.data(t10off + 100 * idx + 19);

                auto t10_xxz_xxx = primBuffer.data(t10off + 100 * idx + 20);

                auto t10_xxz_xxy = primBuffer.data(t10off + 100 * idx + 21);

                auto t10_xxz_xxz = primBuffer.data(t10off + 100 * idx + 22);

                auto t10_xxz_xyy = primBuffer.data(t10off + 100 * idx + 23);

                auto t10_xxz_xyz = primBuffer.data(t10off + 100 * idx + 24);

                auto t10_xxz_xzz = primBuffer.data(t10off + 100 * idx + 25);

                auto t10_xxz_yyy = primBuffer.data(t10off + 100 * idx + 26);

                auto t10_xxz_yyz = primBuffer.data(t10off + 100 * idx + 27);

                auto t10_xxz_yzz = primBuffer.data(t10off + 100 * idx + 28);

                auto t10_xxz_zzz = primBuffer.data(t10off + 100 * idx + 29);

                auto t10_xyy_xxx = primBuffer.data(t10off + 100 * idx + 30);

                auto t10_xyy_xxy = primBuffer.data(t10off + 100 * idx + 31);

                auto t10_xyy_xxz = primBuffer.data(t10off + 100 * idx + 32);

                auto t10_xyy_xyy = primBuffer.data(t10off + 100 * idx + 33);

                auto t10_xyy_xyz = primBuffer.data(t10off + 100 * idx + 34);

                auto t10_xyy_xzz = primBuffer.data(t10off + 100 * idx + 35);

                auto t10_xyy_yyy = primBuffer.data(t10off + 100 * idx + 36);

                auto t10_xyy_yyz = primBuffer.data(t10off + 100 * idx + 37);

                auto t10_xyy_yzz = primBuffer.data(t10off + 100 * idx + 38);

                auto t10_xyy_zzz = primBuffer.data(t10off + 100 * idx + 39);

                auto t10_xyz_xxx = primBuffer.data(t10off + 100 * idx + 40);

                auto t10_xyz_xxy = primBuffer.data(t10off + 100 * idx + 41);

                auto t10_xyz_xxz = primBuffer.data(t10off + 100 * idx + 42);

                auto t10_xyz_xyy = primBuffer.data(t10off + 100 * idx + 43);

                auto t10_xyz_xyz = primBuffer.data(t10off + 100 * idx + 44);

                auto t10_xyz_xzz = primBuffer.data(t10off + 100 * idx + 45);

                auto t10_xyz_yyy = primBuffer.data(t10off + 100 * idx + 46);

                auto t10_xyz_yyz = primBuffer.data(t10off + 100 * idx + 47);

                auto t10_xyz_yzz = primBuffer.data(t10off + 100 * idx + 48);

                auto t10_xyz_zzz = primBuffer.data(t10off + 100 * idx + 49);

                auto t10_xzz_xxx = primBuffer.data(t10off + 100 * idx + 50);

                auto t10_xzz_xxy = primBuffer.data(t10off + 100 * idx + 51);

                auto t10_xzz_xxz = primBuffer.data(t10off + 100 * idx + 52);

                auto t10_xzz_xyy = primBuffer.data(t10off + 100 * idx + 53);

                auto t10_xzz_xyz = primBuffer.data(t10off + 100 * idx + 54);

                auto t10_xzz_xzz = primBuffer.data(t10off + 100 * idx + 55);

                auto t10_xzz_yyy = primBuffer.data(t10off + 100 * idx + 56);

                auto t10_xzz_yyz = primBuffer.data(t10off + 100 * idx + 57);

                auto t10_xzz_yzz = primBuffer.data(t10off + 100 * idx + 58);

                auto t10_xzz_zzz = primBuffer.data(t10off + 100 * idx + 59);

                auto t10_yyy_xxx = primBuffer.data(t10off + 100 * idx + 60);

                auto t10_yyy_xxy = primBuffer.data(t10off + 100 * idx + 61);

                auto t10_yyy_xxz = primBuffer.data(t10off + 100 * idx + 62);

                auto t10_yyy_xyy = primBuffer.data(t10off + 100 * idx + 63);

                auto t10_yyy_xyz = primBuffer.data(t10off + 100 * idx + 64);

                auto t10_yyy_xzz = primBuffer.data(t10off + 100 * idx + 65);

                auto t10_yyy_yyy = primBuffer.data(t10off + 100 * idx + 66);

                auto t10_yyy_yyz = primBuffer.data(t10off + 100 * idx + 67);

                auto t10_yyy_yzz = primBuffer.data(t10off + 100 * idx + 68);

                auto t10_yyy_zzz = primBuffer.data(t10off + 100 * idx + 69);

                auto t10_yyz_xxx = primBuffer.data(t10off + 100 * idx + 70);

                auto t10_yyz_xxy = primBuffer.data(t10off + 100 * idx + 71);

                auto t10_yyz_xxz = primBuffer.data(t10off + 100 * idx + 72);

                auto t10_yyz_xyy = primBuffer.data(t10off + 100 * idx + 73);

                auto t10_yyz_xyz = primBuffer.data(t10off + 100 * idx + 74);

                auto t10_yyz_xzz = primBuffer.data(t10off + 100 * idx + 75);

                auto t10_yyz_yyy = primBuffer.data(t10off + 100 * idx + 76);

                auto t10_yyz_yyz = primBuffer.data(t10off + 100 * idx + 77);

                auto t10_yyz_yzz = primBuffer.data(t10off + 100 * idx + 78);

                auto t10_yyz_zzz = primBuffer.data(t10off + 100 * idx + 79);

                auto t10_yzz_xxx = primBuffer.data(t10off + 100 * idx + 80);

                auto t10_yzz_xxy = primBuffer.data(t10off + 100 * idx + 81);

                auto t10_yzz_xxz = primBuffer.data(t10off + 100 * idx + 82);

                auto t10_yzz_xyy = primBuffer.data(t10off + 100 * idx + 83);

                auto t10_yzz_xyz = primBuffer.data(t10off + 100 * idx + 84);

                auto t10_yzz_xzz = primBuffer.data(t10off + 100 * idx + 85);

                auto t10_yzz_yyy = primBuffer.data(t10off + 100 * idx + 86);

                auto t10_yzz_yyz = primBuffer.data(t10off + 100 * idx + 87);

                auto t10_yzz_yzz = primBuffer.data(t10off + 100 * idx + 88);

                auto t10_yzz_zzz = primBuffer.data(t10off + 100 * idx + 89);

                auto t10_zzz_xxx = primBuffer.data(t10off + 100 * idx + 90);

                auto t10_zzz_xxy = primBuffer.data(t10off + 100 * idx + 91);

                auto t10_zzz_xxz = primBuffer.data(t10off + 100 * idx + 92);

                auto t10_zzz_xyy = primBuffer.data(t10off + 100 * idx + 93);

                auto t10_zzz_xyz = primBuffer.data(t10off + 100 * idx + 94);

                auto t10_zzz_xzz = primBuffer.data(t10off + 100 * idx + 95);

                auto t10_zzz_yyy = primBuffer.data(t10off + 100 * idx + 96);

                auto t10_zzz_yyz = primBuffer.data(t10off + 100 * idx + 97);

                auto t10_zzz_yzz = primBuffer.data(t10off + 100 * idx + 98);

                auto t10_zzz_zzz = primBuffer.data(t10off + 100 * idx + 99);

                // set up pointers to (F|A(0)|F)^(m+1) integrals

                auto t11_xxx_xxx = primBuffer.data(t11off + 100 * idx);

                auto t11_xxx_xxy = primBuffer.data(t11off + 100 * idx + 1);

                auto t11_xxx_xxz = primBuffer.data(t11off + 100 * idx + 2);

                auto t11_xxx_xyy = primBuffer.data(t11off + 100 * idx + 3);

                auto t11_xxx_xyz = primBuffer.data(t11off + 100 * idx + 4);

                auto t11_xxx_xzz = primBuffer.data(t11off + 100 * idx + 5);

                auto t11_xxx_yyy = primBuffer.data(t11off + 100 * idx + 6);

                auto t11_xxx_yyz = primBuffer.data(t11off + 100 * idx + 7);

                auto t11_xxx_yzz = primBuffer.data(t11off + 100 * idx + 8);

                auto t11_xxx_zzz = primBuffer.data(t11off + 100 * idx + 9);

                auto t11_xxy_xxx = primBuffer.data(t11off + 100 * idx + 10);

                auto t11_xxy_xxy = primBuffer.data(t11off + 100 * idx + 11);

                auto t11_xxy_xxz = primBuffer.data(t11off + 100 * idx + 12);

                auto t11_xxy_xyy = primBuffer.data(t11off + 100 * idx + 13);

                auto t11_xxy_xyz = primBuffer.data(t11off + 100 * idx + 14);

                auto t11_xxy_xzz = primBuffer.data(t11off + 100 * idx + 15);

                auto t11_xxy_yyy = primBuffer.data(t11off + 100 * idx + 16);

                auto t11_xxy_yyz = primBuffer.data(t11off + 100 * idx + 17);

                auto t11_xxy_yzz = primBuffer.data(t11off + 100 * idx + 18);

                auto t11_xxy_zzz = primBuffer.data(t11off + 100 * idx + 19);

                auto t11_xxz_xxx = primBuffer.data(t11off + 100 * idx + 20);

                auto t11_xxz_xxy = primBuffer.data(t11off + 100 * idx + 21);

                auto t11_xxz_xxz = primBuffer.data(t11off + 100 * idx + 22);

                auto t11_xxz_xyy = primBuffer.data(t11off + 100 * idx + 23);

                auto t11_xxz_xyz = primBuffer.data(t11off + 100 * idx + 24);

                auto t11_xxz_xzz = primBuffer.data(t11off + 100 * idx + 25);

                auto t11_xxz_yyy = primBuffer.data(t11off + 100 * idx + 26);

                auto t11_xxz_yyz = primBuffer.data(t11off + 100 * idx + 27);

                auto t11_xxz_yzz = primBuffer.data(t11off + 100 * idx + 28);

                auto t11_xxz_zzz = primBuffer.data(t11off + 100 * idx + 29);

                auto t11_xyy_xxx = primBuffer.data(t11off + 100 * idx + 30);

                auto t11_xyy_xxy = primBuffer.data(t11off + 100 * idx + 31);

                auto t11_xyy_xxz = primBuffer.data(t11off + 100 * idx + 32);

                auto t11_xyy_xyy = primBuffer.data(t11off + 100 * idx + 33);

                auto t11_xyy_xyz = primBuffer.data(t11off + 100 * idx + 34);

                auto t11_xyy_xzz = primBuffer.data(t11off + 100 * idx + 35);

                auto t11_xyy_yyy = primBuffer.data(t11off + 100 * idx + 36);

                auto t11_xyy_yyz = primBuffer.data(t11off + 100 * idx + 37);

                auto t11_xyy_yzz = primBuffer.data(t11off + 100 * idx + 38);

                auto t11_xyy_zzz = primBuffer.data(t11off + 100 * idx + 39);

                auto t11_xyz_xxx = primBuffer.data(t11off + 100 * idx + 40);

                auto t11_xyz_xxy = primBuffer.data(t11off + 100 * idx + 41);

                auto t11_xyz_xxz = primBuffer.data(t11off + 100 * idx + 42);

                auto t11_xyz_xyy = primBuffer.data(t11off + 100 * idx + 43);

                auto t11_xyz_xyz = primBuffer.data(t11off + 100 * idx + 44);

                auto t11_xyz_xzz = primBuffer.data(t11off + 100 * idx + 45);

                auto t11_xyz_yyy = primBuffer.data(t11off + 100 * idx + 46);

                auto t11_xyz_yyz = primBuffer.data(t11off + 100 * idx + 47);

                auto t11_xyz_yzz = primBuffer.data(t11off + 100 * idx + 48);

                auto t11_xyz_zzz = primBuffer.data(t11off + 100 * idx + 49);

                auto t11_xzz_xxx = primBuffer.data(t11off + 100 * idx + 50);

                auto t11_xzz_xxy = primBuffer.data(t11off + 100 * idx + 51);

                auto t11_xzz_xxz = primBuffer.data(t11off + 100 * idx + 52);

                auto t11_xzz_xyy = primBuffer.data(t11off + 100 * idx + 53);

                auto t11_xzz_xyz = primBuffer.data(t11off + 100 * idx + 54);

                auto t11_xzz_xzz = primBuffer.data(t11off + 100 * idx + 55);

                auto t11_xzz_yyy = primBuffer.data(t11off + 100 * idx + 56);

                auto t11_xzz_yyz = primBuffer.data(t11off + 100 * idx + 57);

                auto t11_xzz_yzz = primBuffer.data(t11off + 100 * idx + 58);

                auto t11_xzz_zzz = primBuffer.data(t11off + 100 * idx + 59);

                auto t11_yyy_xxx = primBuffer.data(t11off + 100 * idx + 60);

                auto t11_yyy_xxy = primBuffer.data(t11off + 100 * idx + 61);

                auto t11_yyy_xxz = primBuffer.data(t11off + 100 * idx + 62);

                auto t11_yyy_xyy = primBuffer.data(t11off + 100 * idx + 63);

                auto t11_yyy_xyz = primBuffer.data(t11off + 100 * idx + 64);

                auto t11_yyy_xzz = primBuffer.data(t11off + 100 * idx + 65);

                auto t11_yyy_yyy = primBuffer.data(t11off + 100 * idx + 66);

                auto t11_yyy_yyz = primBuffer.data(t11off + 100 * idx + 67);

                auto t11_yyy_yzz = primBuffer.data(t11off + 100 * idx + 68);

                auto t11_yyy_zzz = primBuffer.data(t11off + 100 * idx + 69);

                auto t11_yyz_xxx = primBuffer.data(t11off + 100 * idx + 70);

                auto t11_yyz_xxy = primBuffer.data(t11off + 100 * idx + 71);

                auto t11_yyz_xxz = primBuffer.data(t11off + 100 * idx + 72);

                auto t11_yyz_xyy = primBuffer.data(t11off + 100 * idx + 73);

                auto t11_yyz_xyz = primBuffer.data(t11off + 100 * idx + 74);

                auto t11_yyz_xzz = primBuffer.data(t11off + 100 * idx + 75);

                auto t11_yyz_yyy = primBuffer.data(t11off + 100 * idx + 76);

                auto t11_yyz_yyz = primBuffer.data(t11off + 100 * idx + 77);

                auto t11_yyz_yzz = primBuffer.data(t11off + 100 * idx + 78);

                auto t11_yyz_zzz = primBuffer.data(t11off + 100 * idx + 79);

                auto t11_yzz_xxx = primBuffer.data(t11off + 100 * idx + 80);

                auto t11_yzz_xxy = primBuffer.data(t11off + 100 * idx + 81);

                auto t11_yzz_xxz = primBuffer.data(t11off + 100 * idx + 82);

                auto t11_yzz_xyy = primBuffer.data(t11off + 100 * idx + 83);

                auto t11_yzz_xyz = primBuffer.data(t11off + 100 * idx + 84);

                auto t11_yzz_xzz = primBuffer.data(t11off + 100 * idx + 85);

                auto t11_yzz_yyy = primBuffer.data(t11off + 100 * idx + 86);

                auto t11_yzz_yyz = primBuffer.data(t11off + 100 * idx + 87);

                auto t11_yzz_yzz = primBuffer.data(t11off + 100 * idx + 88);

                auto t11_yzz_zzz = primBuffer.data(t11off + 100 * idx + 89);

                auto t11_zzz_xxx = primBuffer.data(t11off + 100 * idx + 90);

                auto t11_zzz_xxy = primBuffer.data(t11off + 100 * idx + 91);

                auto t11_zzz_xxz = primBuffer.data(t11off + 100 * idx + 92);

                auto t11_zzz_xyy = primBuffer.data(t11off + 100 * idx + 93);

                auto t11_zzz_xyz = primBuffer.data(t11off + 100 * idx + 94);

                auto t11_zzz_xzz = primBuffer.data(t11off + 100 * idx + 95);

                auto t11_zzz_yyy = primBuffer.data(t11off + 100 * idx + 96);

                auto t11_zzz_yyz = primBuffer.data(t11off + 100 * idx + 97);

                auto t11_zzz_yzz = primBuffer.data(t11off + 100 * idx + 98);

                auto t11_zzz_zzz = primBuffer.data(t11off + 100 * idx + 99);

                // set up pointers to (G|A(0)|F)^(m) integrals

                auto t_xxxx_xxx = primBuffer.data(toff + 150 * idx);

                auto t_xxxx_xxy = primBuffer.data(toff + 150 * idx + 1);

                auto t_xxxx_xxz = primBuffer.data(toff + 150 * idx + 2);

                auto t_xxxx_xyy = primBuffer.data(toff + 150 * idx + 3);

                auto t_xxxx_xyz = primBuffer.data(toff + 150 * idx + 4);

                auto t_xxxx_xzz = primBuffer.data(toff + 150 * idx + 5);

                auto t_xxxx_yyy = primBuffer.data(toff + 150 * idx + 6);

                auto t_xxxx_yyz = primBuffer.data(toff + 150 * idx + 7);

                auto t_xxxx_yzz = primBuffer.data(toff + 150 * idx + 8);

                auto t_xxxx_zzz = primBuffer.data(toff + 150 * idx + 9);

                auto t_xxxy_xxx = primBuffer.data(toff + 150 * idx + 10);

                auto t_xxxy_xxy = primBuffer.data(toff + 150 * idx + 11);

                auto t_xxxy_xxz = primBuffer.data(toff + 150 * idx + 12);

                auto t_xxxy_xyy = primBuffer.data(toff + 150 * idx + 13);

                auto t_xxxy_xyz = primBuffer.data(toff + 150 * idx + 14);

                auto t_xxxy_xzz = primBuffer.data(toff + 150 * idx + 15);

                auto t_xxxy_yyy = primBuffer.data(toff + 150 * idx + 16);

                auto t_xxxy_yyz = primBuffer.data(toff + 150 * idx + 17);

                auto t_xxxy_yzz = primBuffer.data(toff + 150 * idx + 18);

                auto t_xxxy_zzz = primBuffer.data(toff + 150 * idx + 19);

                auto t_xxxz_xxx = primBuffer.data(toff + 150 * idx + 20);

                auto t_xxxz_xxy = primBuffer.data(toff + 150 * idx + 21);

                auto t_xxxz_xxz = primBuffer.data(toff + 150 * idx + 22);

                auto t_xxxz_xyy = primBuffer.data(toff + 150 * idx + 23);

                auto t_xxxz_xyz = primBuffer.data(toff + 150 * idx + 24);

                auto t_xxxz_xzz = primBuffer.data(toff + 150 * idx + 25);

                auto t_xxxz_yyy = primBuffer.data(toff + 150 * idx + 26);

                auto t_xxxz_yyz = primBuffer.data(toff + 150 * idx + 27);

                auto t_xxxz_yzz = primBuffer.data(toff + 150 * idx + 28);

                auto t_xxxz_zzz = primBuffer.data(toff + 150 * idx + 29);

                auto t_xxyy_xxx = primBuffer.data(toff + 150 * idx + 30);

                auto t_xxyy_xxy = primBuffer.data(toff + 150 * idx + 31);

                auto t_xxyy_xxz = primBuffer.data(toff + 150 * idx + 32);

                auto t_xxyy_xyy = primBuffer.data(toff + 150 * idx + 33);

                auto t_xxyy_xyz = primBuffer.data(toff + 150 * idx + 34);

                auto t_xxyy_xzz = primBuffer.data(toff + 150 * idx + 35);

                auto t_xxyy_yyy = primBuffer.data(toff + 150 * idx + 36);

                auto t_xxyy_yyz = primBuffer.data(toff + 150 * idx + 37);

                auto t_xxyy_yzz = primBuffer.data(toff + 150 * idx + 38);

                auto t_xxyy_zzz = primBuffer.data(toff + 150 * idx + 39);

                auto t_xxyz_xxx = primBuffer.data(toff + 150 * idx + 40);

                auto t_xxyz_xxy = primBuffer.data(toff + 150 * idx + 41);

                auto t_xxyz_xxz = primBuffer.data(toff + 150 * idx + 42);

                auto t_xxyz_xyy = primBuffer.data(toff + 150 * idx + 43);

                auto t_xxyz_xyz = primBuffer.data(toff + 150 * idx + 44);

                auto t_xxyz_xzz = primBuffer.data(toff + 150 * idx + 45);

                auto t_xxyz_yyy = primBuffer.data(toff + 150 * idx + 46);

                auto t_xxyz_yyz = primBuffer.data(toff + 150 * idx + 47);

                auto t_xxyz_yzz = primBuffer.data(toff + 150 * idx + 48);

                auto t_xxyz_zzz = primBuffer.data(toff + 150 * idx + 49);

                auto t_xxzz_xxx = primBuffer.data(toff + 150 * idx + 50);

                auto t_xxzz_xxy = primBuffer.data(toff + 150 * idx + 51);

                auto t_xxzz_xxz = primBuffer.data(toff + 150 * idx + 52);

                auto t_xxzz_xyy = primBuffer.data(toff + 150 * idx + 53);

                auto t_xxzz_xyz = primBuffer.data(toff + 150 * idx + 54);

                auto t_xxzz_xzz = primBuffer.data(toff + 150 * idx + 55);

                auto t_xxzz_yyy = primBuffer.data(toff + 150 * idx + 56);

                auto t_xxzz_yyz = primBuffer.data(toff + 150 * idx + 57);

                auto t_xxzz_yzz = primBuffer.data(toff + 150 * idx + 58);

                auto t_xxzz_zzz = primBuffer.data(toff + 150 * idx + 59);

                auto t_xyyy_xxx = primBuffer.data(toff + 150 * idx + 60);

                auto t_xyyy_xxy = primBuffer.data(toff + 150 * idx + 61);

                auto t_xyyy_xxz = primBuffer.data(toff + 150 * idx + 62);

                auto t_xyyy_xyy = primBuffer.data(toff + 150 * idx + 63);

                auto t_xyyy_xyz = primBuffer.data(toff + 150 * idx + 64);

                auto t_xyyy_xzz = primBuffer.data(toff + 150 * idx + 65);

                auto t_xyyy_yyy = primBuffer.data(toff + 150 * idx + 66);

                auto t_xyyy_yyz = primBuffer.data(toff + 150 * idx + 67);

                auto t_xyyy_yzz = primBuffer.data(toff + 150 * idx + 68);

                auto t_xyyy_zzz = primBuffer.data(toff + 150 * idx + 69);

                auto t_xyyz_xxx = primBuffer.data(toff + 150 * idx + 70);

                auto t_xyyz_xxy = primBuffer.data(toff + 150 * idx + 71);

                auto t_xyyz_xxz = primBuffer.data(toff + 150 * idx + 72);

                auto t_xyyz_xyy = primBuffer.data(toff + 150 * idx + 73);

                auto t_xyyz_xyz = primBuffer.data(toff + 150 * idx + 74);

                auto t_xyyz_xzz = primBuffer.data(toff + 150 * idx + 75);

                auto t_xyyz_yyy = primBuffer.data(toff + 150 * idx + 76);

                auto t_xyyz_yyz = primBuffer.data(toff + 150 * idx + 77);

                auto t_xyyz_yzz = primBuffer.data(toff + 150 * idx + 78);

                auto t_xyyz_zzz = primBuffer.data(toff + 150 * idx + 79);

                auto t_xyzz_xxx = primBuffer.data(toff + 150 * idx + 80);

                auto t_xyzz_xxy = primBuffer.data(toff + 150 * idx + 81);

                auto t_xyzz_xxz = primBuffer.data(toff + 150 * idx + 82);

                auto t_xyzz_xyy = primBuffer.data(toff + 150 * idx + 83);

                auto t_xyzz_xyz = primBuffer.data(toff + 150 * idx + 84);

                auto t_xyzz_xzz = primBuffer.data(toff + 150 * idx + 85);

                auto t_xyzz_yyy = primBuffer.data(toff + 150 * idx + 86);

                auto t_xyzz_yyz = primBuffer.data(toff + 150 * idx + 87);

                auto t_xyzz_yzz = primBuffer.data(toff + 150 * idx + 88);

                auto t_xyzz_zzz = primBuffer.data(toff + 150 * idx + 89);

                auto t_xzzz_xxx = primBuffer.data(toff + 150 * idx + 90);

                auto t_xzzz_xxy = primBuffer.data(toff + 150 * idx + 91);

                auto t_xzzz_xxz = primBuffer.data(toff + 150 * idx + 92);

                auto t_xzzz_xyy = primBuffer.data(toff + 150 * idx + 93);

                auto t_xzzz_xyz = primBuffer.data(toff + 150 * idx + 94);

                auto t_xzzz_xzz = primBuffer.data(toff + 150 * idx + 95);

                auto t_xzzz_yyy = primBuffer.data(toff + 150 * idx + 96);

                auto t_xzzz_yyz = primBuffer.data(toff + 150 * idx + 97);

                auto t_xzzz_yzz = primBuffer.data(toff + 150 * idx + 98);

                auto t_xzzz_zzz = primBuffer.data(toff + 150 * idx + 99);

                auto t_yyyy_xxx = primBuffer.data(toff + 150 * idx + 100);

                auto t_yyyy_xxy = primBuffer.data(toff + 150 * idx + 101);

                auto t_yyyy_xxz = primBuffer.data(toff + 150 * idx + 102);

                auto t_yyyy_xyy = primBuffer.data(toff + 150 * idx + 103);

                auto t_yyyy_xyz = primBuffer.data(toff + 150 * idx + 104);

                auto t_yyyy_xzz = primBuffer.data(toff + 150 * idx + 105);

                auto t_yyyy_yyy = primBuffer.data(toff + 150 * idx + 106);

                auto t_yyyy_yyz = primBuffer.data(toff + 150 * idx + 107);

                auto t_yyyy_yzz = primBuffer.data(toff + 150 * idx + 108);

                auto t_yyyy_zzz = primBuffer.data(toff + 150 * idx + 109);

                auto t_yyyz_xxx = primBuffer.data(toff + 150 * idx + 110);

                auto t_yyyz_xxy = primBuffer.data(toff + 150 * idx + 111);

                auto t_yyyz_xxz = primBuffer.data(toff + 150 * idx + 112);

                auto t_yyyz_xyy = primBuffer.data(toff + 150 * idx + 113);

                auto t_yyyz_xyz = primBuffer.data(toff + 150 * idx + 114);

                auto t_yyyz_xzz = primBuffer.data(toff + 150 * idx + 115);

                auto t_yyyz_yyy = primBuffer.data(toff + 150 * idx + 116);

                auto t_yyyz_yyz = primBuffer.data(toff + 150 * idx + 117);

                auto t_yyyz_yzz = primBuffer.data(toff + 150 * idx + 118);

                auto t_yyyz_zzz = primBuffer.data(toff + 150 * idx + 119);

                auto t_yyzz_xxx = primBuffer.data(toff + 150 * idx + 120);

                auto t_yyzz_xxy = primBuffer.data(toff + 150 * idx + 121);

                auto t_yyzz_xxz = primBuffer.data(toff + 150 * idx + 122);

                auto t_yyzz_xyy = primBuffer.data(toff + 150 * idx + 123);

                auto t_yyzz_xyz = primBuffer.data(toff + 150 * idx + 124);

                auto t_yyzz_xzz = primBuffer.data(toff + 150 * idx + 125);

                auto t_yyzz_yyy = primBuffer.data(toff + 150 * idx + 126);

                auto t_yyzz_yyz = primBuffer.data(toff + 150 * idx + 127);

                auto t_yyzz_yzz = primBuffer.data(toff + 150 * idx + 128);

                auto t_yyzz_zzz = primBuffer.data(toff + 150 * idx + 129);

                auto t_yzzz_xxx = primBuffer.data(toff + 150 * idx + 130);

                auto t_yzzz_xxy = primBuffer.data(toff + 150 * idx + 131);

                auto t_yzzz_xxz = primBuffer.data(toff + 150 * idx + 132);

                auto t_yzzz_xyy = primBuffer.data(toff + 150 * idx + 133);

                auto t_yzzz_xyz = primBuffer.data(toff + 150 * idx + 134);

                auto t_yzzz_xzz = primBuffer.data(toff + 150 * idx + 135);

                auto t_yzzz_yyy = primBuffer.data(toff + 150 * idx + 136);

                auto t_yzzz_yyz = primBuffer.data(toff + 150 * idx + 137);

                auto t_yzzz_yzz = primBuffer.data(toff + 150 * idx + 138);

                auto t_yzzz_zzz = primBuffer.data(toff + 150 * idx + 139);

                auto t_zzzz_xxx = primBuffer.data(toff + 150 * idx + 140);

                auto t_zzzz_xxy = primBuffer.data(toff + 150 * idx + 141);

                auto t_zzzz_xxz = primBuffer.data(toff + 150 * idx + 142);

                auto t_zzzz_xyy = primBuffer.data(toff + 150 * idx + 143);

                auto t_zzzz_xyz = primBuffer.data(toff + 150 * idx + 144);

                auto t_zzzz_xzz = primBuffer.data(toff + 150 * idx + 145);

                auto t_zzzz_yyy = primBuffer.data(toff + 150 * idx + 146);

                auto t_zzzz_yyz = primBuffer.data(toff + 150 * idx + 147);

                auto t_zzzz_yzz = primBuffer.data(toff + 150 * idx + 148);

                auto t_zzzz_zzz = primBuffer.data(toff + 150 * idx + 149);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xxx_xx, tk0_xxx_xy,\
                                         tk0_xxx_xz, tk0_xxx_yy, tk0_xxx_yz, tk0_xxx_zz,\
                                         tk0_xxy_xx, tk0_xxy_xy, tk0_xxy_xz, tk0_xxy_yy,\
                                         tk0_xxy_yz, tk0_xxy_zz, tk0_xxz_xx, tk0_xxz_xy,\
                                         tk0_xxz_xz, tk0_xxz_yy, tk0_xxz_yz, tk0_xxz_zz,\
                                         tk0_xyy_xx, tk0_xyy_xy, tk0_xyy_xz, tk0_xyy_yy,\
                                         tk0_xyy_yz, tk0_xyy_zz, tk0_xyz_xx, tk0_xyz_xy,\
                                         tk0_xyz_xz, tk0_xyz_yy, tk0_xyz_yz, tk0_xyz_zz,\
                                         tk0_xzz_xx, tk0_xzz_xy, tk0_xzz_xz, tk0_xzz_yy,\
                                         tk0_xzz_yz, tk0_xzz_zz, tk0_yyy_xx, tk0_yyy_xy,\
                                         tk0_yyy_xz, tk0_yyy_yy, tk0_yyy_yz, tk0_yyy_zz,\
                                         tk0_yyz_xx, tk0_yyz_xy, tk0_yyz_xz, tk0_yyz_yy,\
                                         tk0_yyz_yz, tk0_yyz_zz, tk0_yzz_xx, tk0_yzz_xy,\
                                         tk0_yzz_xz, tk0_yzz_yy, tk0_yzz_yz, tk0_yzz_zz,\
                                         tk0_zzz_xx, tk0_zzz_xy, tk0_zzz_xz, tk0_zzz_yy,\
                                         tk0_zzz_yz, tk0_zzz_zz, tk1_xxx_xx, tk1_xxx_xy,\
                                         tk1_xxx_xz, tk1_xxx_yy, tk1_xxx_yz, tk1_xxx_zz,\
                                         tk1_xxy_xx, tk1_xxy_xy, tk1_xxy_xz, tk1_xxy_yy,\
                                         tk1_xxy_yz, tk1_xxy_zz, tk1_xxz_xx, tk1_xxz_xy,\
                                         tk1_xxz_xz, tk1_xxz_yy, tk1_xxz_yz, tk1_xxz_zz,\
                                         tk1_xyy_xx, tk1_xyy_xy, tk1_xyy_xz, tk1_xyy_yy,\
                                         tk1_xyy_yz, tk1_xyy_zz, tk1_xyz_xx, tk1_xyz_xy,\
                                         tk1_xyz_xz, tk1_xyz_yy, tk1_xyz_yz, tk1_xyz_zz,\
                                         tk1_xzz_xx, tk1_xzz_xy, tk1_xzz_xz, tk1_xzz_yy,\
                                         tk1_xzz_yz, tk1_xzz_zz, tk1_yyy_xx, tk1_yyy_xy,\
                                         tk1_yyy_xz, tk1_yyy_yy, tk1_yyy_yz, tk1_yyy_zz,\
                                         tk1_yyz_xx, tk1_yyz_xy, tk1_yyz_xz, tk1_yyz_yy,\
                                         tk1_yyz_yz, tk1_yyz_zz, tk1_yzz_xx, tk1_yzz_xy,\
                                         tk1_yzz_xz, tk1_yzz_yy, tk1_yzz_yz, tk1_yzz_zz,\
                                         tk1_zzz_xx, tk1_zzz_xy, tk1_zzz_xz, tk1_zzz_yy,\
                                         tk1_zzz_yz, tk1_zzz_zz, t20_xx_xxx, t20_xx_xxy,\
                                         t20_xx_xxz, t20_xx_xyy, t20_xx_xyz, t20_xx_xzz,\
                                         t20_xx_yyy, t20_xx_yyz, t20_xx_yzz, t20_xx_zzz,\
                                         t20_xy_xxx, t20_xy_xxy, t20_xy_xxz, t20_xy_xyy,\
                                         t20_xy_xyz, t20_xy_xzz, t20_xy_yyy, t20_xy_yyz,\
                                         t20_xy_yzz, t20_xy_zzz, t20_xz_xxx, t20_xz_xxy,\
                                         t20_xz_xxz, t20_xz_xyy, t20_xz_xyz, t20_xz_xzz,\
                                         t20_xz_yyy, t20_xz_yyz, t20_xz_yzz, t20_xz_zzz,\
                                         t20_yy_xxx, t20_yy_xxy, t20_yy_xxz, t20_yy_xyy,\
                                         t20_yy_xyz, t20_yy_xzz, t20_yy_yyy, t20_yy_yyz,\
                                         t20_yy_yzz, t20_yy_zzz, t20_yz_xxx, t20_yz_xxy,\
                                         t20_yz_xxz, t20_yz_xyy, t20_yz_xyz, t20_yz_xzz,\
                                         t20_yz_yyy, t20_yz_yyz, t20_yz_yzz, t20_yz_zzz,\
                                         t20_zz_xxx, t20_zz_xxy, t20_zz_xxz, t20_zz_xyy,\
                                         t20_zz_xyz, t20_zz_xzz, t20_zz_yyy, t20_zz_yyz,\
                                         t20_zz_yzz, t20_zz_zzz, t21_xx_xxx, t21_xx_xxy,\
                                         t21_xx_xxz, t21_xx_xyy, t21_xx_xyz, t21_xx_xzz,\
                                         t21_xx_yyy, t21_xx_yyz, t21_xx_yzz, t21_xx_zzz,\
                                         t21_xy_xxx, t21_xy_xxy, t21_xy_xxz, t21_xy_xyy,\
                                         t21_xy_xyz, t21_xy_xzz, t21_xy_yyy, t21_xy_yyz,\
                                         t21_xy_yzz, t21_xy_zzz, t21_xz_xxx, t21_xz_xxy,\
                                         t21_xz_xxz, t21_xz_xyy, t21_xz_xyz, t21_xz_xzz,\
                                         t21_xz_yyy, t21_xz_yyz, t21_xz_yzz, t21_xz_zzz,\
                                         t21_yy_xxx, t21_yy_xxy, t21_yy_xxz, t21_yy_xyy,\
                                         t21_yy_xyz, t21_yy_xzz, t21_yy_yyy, t21_yy_yyz,\
                                         t21_yy_yzz, t21_yy_zzz, t21_yz_xxx, t21_yz_xxy,\
                                         t21_yz_xxz, t21_yz_xyy, t21_yz_xyz, t21_yz_xzz,\
                                         t21_yz_yyy, t21_yz_yyz, t21_yz_yzz, t21_yz_zzz,\
                                         t21_zz_xxx, t21_zz_xxy, t21_zz_xxz, t21_zz_xyy,\
                                         t21_zz_xyz, t21_zz_xzz, t21_zz_yyy, t21_zz_yyz,\
                                         t21_zz_yzz, t21_zz_zzz, t10_xxx_xxx, t10_xxx_xxy,\
                                         t10_xxx_xxz, t10_xxx_xyy, t10_xxx_xyz,\
                                         t10_xxx_xzz, t10_xxx_yyy, t10_xxx_yyz,\
                                         t10_xxx_yzz, t10_xxx_zzz, t10_xxy_xxx,\
                                         t10_xxy_xxy, t10_xxy_xxz, t10_xxy_xyy,\
                                         t10_xxy_xyz, t10_xxy_xzz, t10_xxy_yyy,\
                                         t10_xxy_yyz, t10_xxy_yzz, t10_xxy_zzz,\
                                         t10_xxz_xxx, t10_xxz_xxy, t10_xxz_xxz,\
                                         t10_xxz_xyy, t10_xxz_xyz, t10_xxz_xzz,\
                                         t10_xxz_yyy, t10_xxz_yyz, t10_xxz_yzz,\
                                         t10_xxz_zzz, t10_xyy_xxx, t10_xyy_xxy,\
                                         t10_xyy_xxz, t10_xyy_xyy, t10_xyy_xyz,\
                                         t10_xyy_xzz, t10_xyy_yyy, t10_xyy_yyz,\
                                         t10_xyy_yzz, t10_xyy_zzz, t10_xyz_xxx,\
                                         t10_xyz_xxy, t10_xyz_xxz, t10_xyz_xyy,\
                                         t10_xyz_xyz, t10_xyz_xzz, t10_xyz_yyy,\
                                         t10_xyz_yyz, t10_xyz_yzz, t10_xyz_zzz,\
                                         t10_xzz_xxx, t10_xzz_xxy, t10_xzz_xxz,\
                                         t10_xzz_xyy, t10_xzz_xyz, t10_xzz_xzz,\
                                         t10_xzz_yyy, t10_xzz_yyz, t10_xzz_yzz,\
                                         t10_xzz_zzz, t10_yyy_xxx, t10_yyy_xxy,\
                                         t10_yyy_xxz, t10_yyy_xyy, t10_yyy_xyz,\
                                         t10_yyy_xzz, t10_yyy_yyy, t10_yyy_yyz,\
                                         t10_yyy_yzz, t10_yyy_zzz, t10_yyz_xxx,\
                                         t10_yyz_xxy, t10_yyz_xxz, t10_yyz_xyy,\
                                         t10_yyz_xyz, t10_yyz_xzz, t10_yyz_yyy,\
                                         t10_yyz_yyz, t10_yyz_yzz, t10_yyz_zzz,\
                                         t10_yzz_xxx, t10_yzz_xxy, t10_yzz_xxz,\
                                         t10_yzz_xyy, t10_yzz_xyz, t10_yzz_xzz,\
                                         t10_yzz_yyy, t10_yzz_yyz, t10_yzz_yzz,\
                                         t10_yzz_zzz, t10_zzz_xxx, t10_zzz_xxy,\
                                         t10_zzz_xxz, t10_zzz_xyy, t10_zzz_xyz,\
                                         t10_zzz_xzz, t10_zzz_yyy, t10_zzz_yyz,\
                                         t10_zzz_yzz, t10_zzz_zzz, t11_xxx_xxx,\
                                         t11_xxx_xxy, t11_xxx_xxz, t11_xxx_xyy,\
                                         t11_xxx_xyz, t11_xxx_xzz, t11_xxx_yyy,\
                                         t11_xxx_yyz, t11_xxx_yzz, t11_xxx_zzz,\
                                         t11_xxy_xxx, t11_xxy_xxy, t11_xxy_xxz,\
                                         t11_xxy_xyy, t11_xxy_xyz, t11_xxy_xzz,\
                                         t11_xxy_yyy, t11_xxy_yyz, t11_xxy_yzz,\
                                         t11_xxy_zzz, t11_xxz_xxx, t11_xxz_xxy,\
                                         t11_xxz_xxz, t11_xxz_xyy, t11_xxz_xyz,\
                                         t11_xxz_xzz, t11_xxz_yyy, t11_xxz_yyz,\
                                         t11_xxz_yzz, t11_xxz_zzz, t11_xyy_xxx,\
                                         t11_xyy_xxy, t11_xyy_xxz, t11_xyy_xyy,\
                                         t11_xyy_xyz, t11_xyy_xzz, t11_xyy_yyy,\
                                         t11_xyy_yyz, t11_xyy_yzz, t11_xyy_zzz,\
                                         t11_xyz_xxx, t11_xyz_xxy, t11_xyz_xxz,\
                                         t11_xyz_xyy, t11_xyz_xyz, t11_xyz_xzz,\
                                         t11_xyz_yyy, t11_xyz_yyz, t11_xyz_yzz,\
                                         t11_xyz_zzz, t11_xzz_xxx, t11_xzz_xxy,\
                                         t11_xzz_xxz, t11_xzz_xyy, t11_xzz_xyz,\
                                         t11_xzz_xzz, t11_xzz_yyy, t11_xzz_yyz,\
                                         t11_xzz_yzz, t11_xzz_zzz, t11_yyy_xxx,\
                                         t11_yyy_xxy, t11_yyy_xxz, t11_yyy_xyy,\
                                         t11_yyy_xyz, t11_yyy_xzz, t11_yyy_yyy,\
                                         t11_yyy_yyz, t11_yyy_yzz, t11_yyy_zzz,\
                                         t11_yyz_xxx, t11_yyz_xxy, t11_yyz_xxz,\
                                         t11_yyz_xyy, t11_yyz_xyz, t11_yyz_xzz,\
                                         t11_yyz_yyy, t11_yyz_yyz, t11_yyz_yzz,\
                                         t11_yyz_zzz, t11_yzz_xxx, t11_yzz_xxy,\
                                         t11_yzz_xxz, t11_yzz_xyy, t11_yzz_xyz,\
                                         t11_yzz_xzz, t11_yzz_yyy, t11_yzz_yyz,\
                                         t11_yzz_yzz, t11_yzz_zzz, t11_zzz_xxx,\
                                         t11_zzz_xxy, t11_zzz_xxz, t11_zzz_xyy,\
                                         t11_zzz_xyz, t11_zzz_xzz, t11_zzz_yyy,\
                                         t11_zzz_yyz, t11_zzz_yzz, t11_zzz_zzz,\
                                         t_xxxx_xxx, t_xxxx_xxy, t_xxxx_xxz, t_xxxx_xyy,\
                                         t_xxxx_xyz, t_xxxx_xzz, t_xxxx_yyy, t_xxxx_yyz,\
                                         t_xxxx_yzz, t_xxxx_zzz, t_xxxy_xxx, t_xxxy_xxy,\
                                         t_xxxy_xxz, t_xxxy_xyy, t_xxxy_xyz, t_xxxy_xzz,\
                                         t_xxxy_yyy, t_xxxy_yyz, t_xxxy_yzz, t_xxxy_zzz,\
                                         t_xxxz_xxx, t_xxxz_xxy, t_xxxz_xxz, t_xxxz_xyy,\
                                         t_xxxz_xyz, t_xxxz_xzz, t_xxxz_yyy, t_xxxz_yyz,\
                                         t_xxxz_yzz, t_xxxz_zzz, t_xxyy_xxx, t_xxyy_xxy,\
                                         t_xxyy_xxz, t_xxyy_xyy, t_xxyy_xyz, t_xxyy_xzz,\
                                         t_xxyy_yyy, t_xxyy_yyz, t_xxyy_yzz, t_xxyy_zzz,\
                                         t_xxyz_xxx, t_xxyz_xxy, t_xxyz_xxz, t_xxyz_xyy,\
                                         t_xxyz_xyz, t_xxyz_xzz, t_xxyz_yyy, t_xxyz_yyz,\
                                         t_xxyz_yzz, t_xxyz_zzz, t_xxzz_xxx, t_xxzz_xxy,\
                                         t_xxzz_xxz, t_xxzz_xyy, t_xxzz_xyz, t_xxzz_xzz,\
                                         t_xxzz_yyy, t_xxzz_yyz, t_xxzz_yzz, t_xxzz_zzz,\
                                         t_xyyy_xxx, t_xyyy_xxy, t_xyyy_xxz, t_xyyy_xyy,\
                                         t_xyyy_xyz, t_xyyy_xzz, t_xyyy_yyy, t_xyyy_yyz,\
                                         t_xyyy_yzz, t_xyyy_zzz, t_xyyz_xxx, t_xyyz_xxy,\
                                         t_xyyz_xxz, t_xyyz_xyy, t_xyyz_xyz, t_xyyz_xzz,\
                                         t_xyyz_yyy, t_xyyz_yyz, t_xyyz_yzz, t_xyyz_zzz,\
                                         t_xyzz_xxx, t_xyzz_xxy, t_xyzz_xxz, t_xyzz_xyy,\
                                         t_xyzz_xyz, t_xyzz_xzz, t_xyzz_yyy, t_xyzz_yyz,\
                                         t_xyzz_yzz, t_xyzz_zzz, t_xzzz_xxx, t_xzzz_xxy,\
                                         t_xzzz_xxz, t_xzzz_xyy, t_xzzz_xyz, t_xzzz_xzz,\
                                         t_xzzz_yyy, t_xzzz_yyz, t_xzzz_yzz, t_xzzz_zzz,\
                                         t_yyyy_xxx, t_yyyy_xxy, t_yyyy_xxz, t_yyyy_xyy,\
                                         t_yyyy_xyz, t_yyyy_xzz, t_yyyy_yyy, t_yyyy_yyz,\
                                         t_yyyy_yzz, t_yyyy_zzz, t_yyyz_xxx, t_yyyz_xxy,\
                                         t_yyyz_xxz, t_yyyz_xyy, t_yyyz_xyz, t_yyyz_xzz,\
                                         t_yyyz_yyy, t_yyyz_yyz, t_yyyz_yzz, t_yyyz_zzz,\
                                         t_yyzz_xxx, t_yyzz_xxy, t_yyzz_xxz, t_yyzz_xyy,\
                                         t_yyzz_xyz, t_yyzz_xzz, t_yyzz_yyy, t_yyzz_yyz,\
                                         t_yyzz_yzz, t_yyzz_zzz, t_yzzz_xxx, t_yzzz_xxy,\
                                         t_yzzz_xxz, t_yzzz_xyy, t_yzzz_xyz, t_yzzz_xzz,\
                                         t_yzzz_yyy, t_yzzz_yyz, t_yzzz_yzz, t_yzzz_zzz,\
                                         t_zzzz_xxx, t_zzzz_xxy, t_zzzz_xxz, t_zzzz_xyy,\
                                         t_zzzz_xyz, t_zzzz_xzz, t_zzzz_yyy, t_zzzz_yyz,\
                                         t_zzzz_yzz, t_zzzz_zzz, pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxxx_xxx[k] = fra * t10_xxx_xxx[k] - frc * t11_xxx_xxx[k] + f2t * (3.0 * t20_xx_xxx[k] - 3.0 * t21_xx_xxx[k] + 3.0 * tk0_xxx_xx[k] - 3.0 * tk1_xxx_xx[k]);

                    t_xxxx_xxy[k] = fra * t10_xxx_xxy[k] - frc * t11_xxx_xxy[k] + f2t * (3.0 * t20_xx_xxy[k] - 3.0 * t21_xx_xxy[k] + 2.0 * tk0_xxx_xy[k] - 2.0 * tk1_xxx_xy[k]);

                    t_xxxx_xxz[k] = fra * t10_xxx_xxz[k] - frc * t11_xxx_xxz[k] + f2t * (3.0 * t20_xx_xxz[k] - 3.0 * t21_xx_xxz[k] + 2.0 * tk0_xxx_xz[k] - 2.0 * tk1_xxx_xz[k]);

                    t_xxxx_xyy[k] = fra * t10_xxx_xyy[k] - frc * t11_xxx_xyy[k] + f2t * (3.0 * t20_xx_xyy[k] - 3.0 * t21_xx_xyy[k] + tk0_xxx_yy[k] - tk1_xxx_yy[k]);

                    t_xxxx_xyz[k] = fra * t10_xxx_xyz[k] - frc * t11_xxx_xyz[k] + f2t * (3.0 * t20_xx_xyz[k] - 3.0 * t21_xx_xyz[k] + tk0_xxx_yz[k] - tk1_xxx_yz[k]);

                    t_xxxx_xzz[k] = fra * t10_xxx_xzz[k] - frc * t11_xxx_xzz[k] + f2t * (3.0 * t20_xx_xzz[k] - 3.0 * t21_xx_xzz[k] + tk0_xxx_zz[k] - tk1_xxx_zz[k]);

                    t_xxxx_yyy[k] = fra * t10_xxx_yyy[k] - frc * t11_xxx_yyy[k] + f2t * (3.0 * t20_xx_yyy[k] - 3.0 * t21_xx_yyy[k]);

                    t_xxxx_yyz[k] = fra * t10_xxx_yyz[k] - frc * t11_xxx_yyz[k] + f2t * (3.0 * t20_xx_yyz[k] - 3.0 * t21_xx_yyz[k]);

                    t_xxxx_yzz[k] = fra * t10_xxx_yzz[k] - frc * t11_xxx_yzz[k] + f2t * (3.0 * t20_xx_yzz[k] - 3.0 * t21_xx_yzz[k]);

                    t_xxxx_zzz[k] = fra * t10_xxx_zzz[k] - frc * t11_xxx_zzz[k] + f2t * (3.0 * t20_xx_zzz[k] - 3.0 * t21_xx_zzz[k]);

                    t_xxxy_xxx[k] = fra * t10_xxy_xxx[k] - frc * t11_xxy_xxx[k] + f2t * (2.0 * t20_xy_xxx[k] - 2.0 * t21_xy_xxx[k] + 3.0 * tk0_xxy_xx[k] - 3.0 * tk1_xxy_xx[k]);

                    t_xxxy_xxy[k] = fra * t10_xxy_xxy[k] - frc * t11_xxy_xxy[k] + f2t * (2.0 * t20_xy_xxy[k] - 2.0 * t21_xy_xxy[k] + 2.0 * tk0_xxy_xy[k] - 2.0 * tk1_xxy_xy[k]);

                    t_xxxy_xxz[k] = fra * t10_xxy_xxz[k] - frc * t11_xxy_xxz[k] + f2t * (2.0 * t20_xy_xxz[k] - 2.0 * t21_xy_xxz[k] + 2.0 * tk0_xxy_xz[k] - 2.0 * tk1_xxy_xz[k]);

                    t_xxxy_xyy[k] = fra * t10_xxy_xyy[k] - frc * t11_xxy_xyy[k] + f2t * (2.0 * t20_xy_xyy[k] - 2.0 * t21_xy_xyy[k] + tk0_xxy_yy[k] - tk1_xxy_yy[k]);

                    t_xxxy_xyz[k] = fra * t10_xxy_xyz[k] - frc * t11_xxy_xyz[k] + f2t * (2.0 * t20_xy_xyz[k] - 2.0 * t21_xy_xyz[k] + tk0_xxy_yz[k] - tk1_xxy_yz[k]);

                    t_xxxy_xzz[k] = fra * t10_xxy_xzz[k] - frc * t11_xxy_xzz[k] + f2t * (2.0 * t20_xy_xzz[k] - 2.0 * t21_xy_xzz[k] + tk0_xxy_zz[k] - tk1_xxy_zz[k]);

                    t_xxxy_yyy[k] = fra * t10_xxy_yyy[k] - frc * t11_xxy_yyy[k] + f2t * (2.0 * t20_xy_yyy[k] - 2.0 * t21_xy_yyy[k]);

                    t_xxxy_yyz[k] = fra * t10_xxy_yyz[k] - frc * t11_xxy_yyz[k] + f2t * (2.0 * t20_xy_yyz[k] - 2.0 * t21_xy_yyz[k]);

                    t_xxxy_yzz[k] = fra * t10_xxy_yzz[k] - frc * t11_xxy_yzz[k] + f2t * (2.0 * t20_xy_yzz[k] - 2.0 * t21_xy_yzz[k]);

                    t_xxxy_zzz[k] = fra * t10_xxy_zzz[k] - frc * t11_xxy_zzz[k] + f2t * (2.0 * t20_xy_zzz[k] - 2.0 * t21_xy_zzz[k]);

                    t_xxxz_xxx[k] = fra * t10_xxz_xxx[k] - frc * t11_xxz_xxx[k] + f2t * (2.0 * t20_xz_xxx[k] - 2.0 * t21_xz_xxx[k] + 3.0 * tk0_xxz_xx[k] - 3.0 * tk1_xxz_xx[k]);

                    t_xxxz_xxy[k] = fra * t10_xxz_xxy[k] - frc * t11_xxz_xxy[k] + f2t * (2.0 * t20_xz_xxy[k] - 2.0 * t21_xz_xxy[k] + 2.0 * tk0_xxz_xy[k] - 2.0 * tk1_xxz_xy[k]);

                    t_xxxz_xxz[k] = fra * t10_xxz_xxz[k] - frc * t11_xxz_xxz[k] + f2t * (2.0 * t20_xz_xxz[k] - 2.0 * t21_xz_xxz[k] + 2.0 * tk0_xxz_xz[k] - 2.0 * tk1_xxz_xz[k]);

                    t_xxxz_xyy[k] = fra * t10_xxz_xyy[k] - frc * t11_xxz_xyy[k] + f2t * (2.0 * t20_xz_xyy[k] - 2.0 * t21_xz_xyy[k] + tk0_xxz_yy[k] - tk1_xxz_yy[k]);

                    t_xxxz_xyz[k] = fra * t10_xxz_xyz[k] - frc * t11_xxz_xyz[k] + f2t * (2.0 * t20_xz_xyz[k] - 2.0 * t21_xz_xyz[k] + tk0_xxz_yz[k] - tk1_xxz_yz[k]);

                    t_xxxz_xzz[k] = fra * t10_xxz_xzz[k] - frc * t11_xxz_xzz[k] + f2t * (2.0 * t20_xz_xzz[k] - 2.0 * t21_xz_xzz[k] + tk0_xxz_zz[k] - tk1_xxz_zz[k]);

                    t_xxxz_yyy[k] = fra * t10_xxz_yyy[k] - frc * t11_xxz_yyy[k] + f2t * (2.0 * t20_xz_yyy[k] - 2.0 * t21_xz_yyy[k]);

                    t_xxxz_yyz[k] = fra * t10_xxz_yyz[k] - frc * t11_xxz_yyz[k] + f2t * (2.0 * t20_xz_yyz[k] - 2.0 * t21_xz_yyz[k]);

                    t_xxxz_yzz[k] = fra * t10_xxz_yzz[k] - frc * t11_xxz_yzz[k] + f2t * (2.0 * t20_xz_yzz[k] - 2.0 * t21_xz_yzz[k]);

                    t_xxxz_zzz[k] = fra * t10_xxz_zzz[k] - frc * t11_xxz_zzz[k] + f2t * (2.0 * t20_xz_zzz[k] - 2.0 * t21_xz_zzz[k]);

                    t_xxyy_xxx[k] = fra * t10_xyy_xxx[k] - frc * t11_xyy_xxx[k] + f2t * (t20_yy_xxx[k] - t21_yy_xxx[k] + 3.0 * tk0_xyy_xx[k] - 3.0 * tk1_xyy_xx[k]);

                    t_xxyy_xxy[k] = fra * t10_xyy_xxy[k] - frc * t11_xyy_xxy[k] + f2t * (t20_yy_xxy[k] - t21_yy_xxy[k] + 2.0 * tk0_xyy_xy[k] - 2.0 * tk1_xyy_xy[k]);

                    t_xxyy_xxz[k] = fra * t10_xyy_xxz[k] - frc * t11_xyy_xxz[k] + f2t * (t20_yy_xxz[k] - t21_yy_xxz[k] + 2.0 * tk0_xyy_xz[k] - 2.0 * tk1_xyy_xz[k]);

                    t_xxyy_xyy[k] = fra * t10_xyy_xyy[k] - frc * t11_xyy_xyy[k] + f2t * (t20_yy_xyy[k] - t21_yy_xyy[k] + tk0_xyy_yy[k] - tk1_xyy_yy[k]);

                    t_xxyy_xyz[k] = fra * t10_xyy_xyz[k] - frc * t11_xyy_xyz[k] + f2t * (t20_yy_xyz[k] - t21_yy_xyz[k] + tk0_xyy_yz[k] - tk1_xyy_yz[k]);

                    t_xxyy_xzz[k] = fra * t10_xyy_xzz[k] - frc * t11_xyy_xzz[k] + f2t * (t20_yy_xzz[k] - t21_yy_xzz[k] + tk0_xyy_zz[k] - tk1_xyy_zz[k]);

                    t_xxyy_yyy[k] = fra * t10_xyy_yyy[k] - frc * t11_xyy_yyy[k] + f2t * (t20_yy_yyy[k] - t21_yy_yyy[k]);

                    t_xxyy_yyz[k] = fra * t10_xyy_yyz[k] - frc * t11_xyy_yyz[k] + f2t * (t20_yy_yyz[k] - t21_yy_yyz[k]);

                    t_xxyy_yzz[k] = fra * t10_xyy_yzz[k] - frc * t11_xyy_yzz[k] + f2t * (t20_yy_yzz[k] - t21_yy_yzz[k]);

                    t_xxyy_zzz[k] = fra * t10_xyy_zzz[k] - frc * t11_xyy_zzz[k] + f2t * (t20_yy_zzz[k] - t21_yy_zzz[k]);

                    t_xxyz_xxx[k] = fra * t10_xyz_xxx[k] - frc * t11_xyz_xxx[k] + f2t * (t20_yz_xxx[k] - t21_yz_xxx[k] + 3.0 * tk0_xyz_xx[k] - 3.0 * tk1_xyz_xx[k]);

                    t_xxyz_xxy[k] = fra * t10_xyz_xxy[k] - frc * t11_xyz_xxy[k] + f2t * (t20_yz_xxy[k] - t21_yz_xxy[k] + 2.0 * tk0_xyz_xy[k] - 2.0 * tk1_xyz_xy[k]);

                    t_xxyz_xxz[k] = fra * t10_xyz_xxz[k] - frc * t11_xyz_xxz[k] + f2t * (t20_yz_xxz[k] - t21_yz_xxz[k] + 2.0 * tk0_xyz_xz[k] - 2.0 * tk1_xyz_xz[k]);

                    t_xxyz_xyy[k] = fra * t10_xyz_xyy[k] - frc * t11_xyz_xyy[k] + f2t * (t20_yz_xyy[k] - t21_yz_xyy[k] + tk0_xyz_yy[k] - tk1_xyz_yy[k]);

                    t_xxyz_xyz[k] = fra * t10_xyz_xyz[k] - frc * t11_xyz_xyz[k] + f2t * (t20_yz_xyz[k] - t21_yz_xyz[k] + tk0_xyz_yz[k] - tk1_xyz_yz[k]);

                    t_xxyz_xzz[k] = fra * t10_xyz_xzz[k] - frc * t11_xyz_xzz[k] + f2t * (t20_yz_xzz[k] - t21_yz_xzz[k] + tk0_xyz_zz[k] - tk1_xyz_zz[k]);

                    t_xxyz_yyy[k] = fra * t10_xyz_yyy[k] - frc * t11_xyz_yyy[k] + f2t * (t20_yz_yyy[k] - t21_yz_yyy[k]);

                    t_xxyz_yyz[k] = fra * t10_xyz_yyz[k] - frc * t11_xyz_yyz[k] + f2t * (t20_yz_yyz[k] - t21_yz_yyz[k]);

                    t_xxyz_yzz[k] = fra * t10_xyz_yzz[k] - frc * t11_xyz_yzz[k] + f2t * (t20_yz_yzz[k] - t21_yz_yzz[k]);

                    t_xxyz_zzz[k] = fra * t10_xyz_zzz[k] - frc * t11_xyz_zzz[k] + f2t * (t20_yz_zzz[k] - t21_yz_zzz[k]);

                    t_xxzz_xxx[k] = fra * t10_xzz_xxx[k] - frc * t11_xzz_xxx[k] + f2t * (t20_zz_xxx[k] - t21_zz_xxx[k] + 3.0 * tk0_xzz_xx[k] - 3.0 * tk1_xzz_xx[k]);

                    t_xxzz_xxy[k] = fra * t10_xzz_xxy[k] - frc * t11_xzz_xxy[k] + f2t * (t20_zz_xxy[k] - t21_zz_xxy[k] + 2.0 * tk0_xzz_xy[k] - 2.0 * tk1_xzz_xy[k]);

                    t_xxzz_xxz[k] = fra * t10_xzz_xxz[k] - frc * t11_xzz_xxz[k] + f2t * (t20_zz_xxz[k] - t21_zz_xxz[k] + 2.0 * tk0_xzz_xz[k] - 2.0 * tk1_xzz_xz[k]);

                    t_xxzz_xyy[k] = fra * t10_xzz_xyy[k] - frc * t11_xzz_xyy[k] + f2t * (t20_zz_xyy[k] - t21_zz_xyy[k] + tk0_xzz_yy[k] - tk1_xzz_yy[k]);

                    t_xxzz_xyz[k] = fra * t10_xzz_xyz[k] - frc * t11_xzz_xyz[k] + f2t * (t20_zz_xyz[k] - t21_zz_xyz[k] + tk0_xzz_yz[k] - tk1_xzz_yz[k]);

                    t_xxzz_xzz[k] = fra * t10_xzz_xzz[k] - frc * t11_xzz_xzz[k] + f2t * (t20_zz_xzz[k] - t21_zz_xzz[k] + tk0_xzz_zz[k] - tk1_xzz_zz[k]);

                    t_xxzz_yyy[k] = fra * t10_xzz_yyy[k] - frc * t11_xzz_yyy[k] + f2t * (t20_zz_yyy[k] - t21_zz_yyy[k]);

                    t_xxzz_yyz[k] = fra * t10_xzz_yyz[k] - frc * t11_xzz_yyz[k] + f2t * (t20_zz_yyz[k] - t21_zz_yyz[k]);

                    t_xxzz_yzz[k] = fra * t10_xzz_yzz[k] - frc * t11_xzz_yzz[k] + f2t * (t20_zz_yzz[k] - t21_zz_yzz[k]);

                    t_xxzz_zzz[k] = fra * t10_xzz_zzz[k] - frc * t11_xzz_zzz[k] + f2t * (t20_zz_zzz[k] - t21_zz_zzz[k]);

                    t_xyyy_xxx[k] = fra * t10_yyy_xxx[k] - frc * t11_yyy_xxx[k] + f2t * (3.0 * tk0_yyy_xx[k] - 3.0 * tk1_yyy_xx[k]);

                    t_xyyy_xxy[k] = fra * t10_yyy_xxy[k] - frc * t11_yyy_xxy[k] + f2t * (2.0 * tk0_yyy_xy[k] - 2.0 * tk1_yyy_xy[k]);

                    t_xyyy_xxz[k] = fra * t10_yyy_xxz[k] - frc * t11_yyy_xxz[k] + f2t * (2.0 * tk0_yyy_xz[k] - 2.0 * tk1_yyy_xz[k]);

                    t_xyyy_xyy[k] = fra * t10_yyy_xyy[k] - frc * t11_yyy_xyy[k] + f2t * (tk0_yyy_yy[k] - tk1_yyy_yy[k]);

                    t_xyyy_xyz[k] = fra * t10_yyy_xyz[k] - frc * t11_yyy_xyz[k] + f2t * (tk0_yyy_yz[k] - tk1_yyy_yz[k]);

                    t_xyyy_xzz[k] = fra * t10_yyy_xzz[k] - frc * t11_yyy_xzz[k] + f2t * (tk0_yyy_zz[k] - tk1_yyy_zz[k]);

                    t_xyyy_yyy[k] = fra * t10_yyy_yyy[k] - frc * t11_yyy_yyy[k];

                    t_xyyy_yyz[k] = fra * t10_yyy_yyz[k] - frc * t11_yyy_yyz[k];

                    t_xyyy_yzz[k] = fra * t10_yyy_yzz[k] - frc * t11_yyy_yzz[k];

                    t_xyyy_zzz[k] = fra * t10_yyy_zzz[k] - frc * t11_yyy_zzz[k];

                    t_xyyz_xxx[k] = fra * t10_yyz_xxx[k] - frc * t11_yyz_xxx[k] + f2t * (3.0 * tk0_yyz_xx[k] - 3.0 * tk1_yyz_xx[k]);

                    t_xyyz_xxy[k] = fra * t10_yyz_xxy[k] - frc * t11_yyz_xxy[k] + f2t * (2.0 * tk0_yyz_xy[k] - 2.0 * tk1_yyz_xy[k]);

                    t_xyyz_xxz[k] = fra * t10_yyz_xxz[k] - frc * t11_yyz_xxz[k] + f2t * (2.0 * tk0_yyz_xz[k] - 2.0 * tk1_yyz_xz[k]);

                    t_xyyz_xyy[k] = fra * t10_yyz_xyy[k] - frc * t11_yyz_xyy[k] + f2t * (tk0_yyz_yy[k] - tk1_yyz_yy[k]);

                    t_xyyz_xyz[k] = fra * t10_yyz_xyz[k] - frc * t11_yyz_xyz[k] + f2t * (tk0_yyz_yz[k] - tk1_yyz_yz[k]);

                    t_xyyz_xzz[k] = fra * t10_yyz_xzz[k] - frc * t11_yyz_xzz[k] + f2t * (tk0_yyz_zz[k] - tk1_yyz_zz[k]);

                    t_xyyz_yyy[k] = fra * t10_yyz_yyy[k] - frc * t11_yyz_yyy[k];

                    t_xyyz_yyz[k] = fra * t10_yyz_yyz[k] - frc * t11_yyz_yyz[k];

                    t_xyyz_yzz[k] = fra * t10_yyz_yzz[k] - frc * t11_yyz_yzz[k];

                    t_xyyz_zzz[k] = fra * t10_yyz_zzz[k] - frc * t11_yyz_zzz[k];

                    t_xyzz_xxx[k] = fra * t10_yzz_xxx[k] - frc * t11_yzz_xxx[k] + f2t * (3.0 * tk0_yzz_xx[k] - 3.0 * tk1_yzz_xx[k]);

                    t_xyzz_xxy[k] = fra * t10_yzz_xxy[k] - frc * t11_yzz_xxy[k] + f2t * (2.0 * tk0_yzz_xy[k] - 2.0 * tk1_yzz_xy[k]);

                    t_xyzz_xxz[k] = fra * t10_yzz_xxz[k] - frc * t11_yzz_xxz[k] + f2t * (2.0 * tk0_yzz_xz[k] - 2.0 * tk1_yzz_xz[k]);

                    t_xyzz_xyy[k] = fra * t10_yzz_xyy[k] - frc * t11_yzz_xyy[k] + f2t * (tk0_yzz_yy[k] - tk1_yzz_yy[k]);

                    t_xyzz_xyz[k] = fra * t10_yzz_xyz[k] - frc * t11_yzz_xyz[k] + f2t * (tk0_yzz_yz[k] - tk1_yzz_yz[k]);

                    t_xyzz_xzz[k] = fra * t10_yzz_xzz[k] - frc * t11_yzz_xzz[k] + f2t * (tk0_yzz_zz[k] - tk1_yzz_zz[k]);

                    t_xyzz_yyy[k] = fra * t10_yzz_yyy[k] - frc * t11_yzz_yyy[k];

                    t_xyzz_yyz[k] = fra * t10_yzz_yyz[k] - frc * t11_yzz_yyz[k];

                    t_xyzz_yzz[k] = fra * t10_yzz_yzz[k] - frc * t11_yzz_yzz[k];

                    t_xyzz_zzz[k] = fra * t10_yzz_zzz[k] - frc * t11_yzz_zzz[k];

                    t_xzzz_xxx[k] = fra * t10_zzz_xxx[k] - frc * t11_zzz_xxx[k] + f2t * (3.0 * tk0_zzz_xx[k] - 3.0 * tk1_zzz_xx[k]);

                    t_xzzz_xxy[k] = fra * t10_zzz_xxy[k] - frc * t11_zzz_xxy[k] + f2t * (2.0 * tk0_zzz_xy[k] - 2.0 * tk1_zzz_xy[k]);

                    t_xzzz_xxz[k] = fra * t10_zzz_xxz[k] - frc * t11_zzz_xxz[k] + f2t * (2.0 * tk0_zzz_xz[k] - 2.0 * tk1_zzz_xz[k]);

                    t_xzzz_xyy[k] = fra * t10_zzz_xyy[k] - frc * t11_zzz_xyy[k] + f2t * (tk0_zzz_yy[k] - tk1_zzz_yy[k]);

                    t_xzzz_xyz[k] = fra * t10_zzz_xyz[k] - frc * t11_zzz_xyz[k] + f2t * (tk0_zzz_yz[k] - tk1_zzz_yz[k]);

                    t_xzzz_xzz[k] = fra * t10_zzz_xzz[k] - frc * t11_zzz_xzz[k] + f2t * (tk0_zzz_zz[k] - tk1_zzz_zz[k]);

                    t_xzzz_yyy[k] = fra * t10_zzz_yyy[k] - frc * t11_zzz_yyy[k];

                    t_xzzz_yyz[k] = fra * t10_zzz_yyz[k] - frc * t11_zzz_yyz[k];

                    t_xzzz_yzz[k] = fra * t10_zzz_yzz[k] - frc * t11_zzz_yzz[k];

                    t_xzzz_zzz[k] = fra * t10_zzz_zzz[k] - frc * t11_zzz_zzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyyy_xxx[k] = fra * t10_yyy_xxx[k] - frc * t11_yyy_xxx[k] + f2t * (3.0 * t20_yy_xxx[k] - 3.0 * t21_yy_xxx[k]);

                    t_yyyy_xxy[k] = fra * t10_yyy_xxy[k] - frc * t11_yyy_xxy[k] + f2t * (3.0 * t20_yy_xxy[k] - 3.0 * t21_yy_xxy[k] + tk0_yyy_xx[k] - tk1_yyy_xx[k]);

                    t_yyyy_xxz[k] = fra * t10_yyy_xxz[k] - frc * t11_yyy_xxz[k] + f2t * (3.0 * t20_yy_xxz[k] - 3.0 * t21_yy_xxz[k]);

                    t_yyyy_xyy[k] = fra * t10_yyy_xyy[k] - frc * t11_yyy_xyy[k] + f2t * (3.0 * t20_yy_xyy[k] - 3.0 * t21_yy_xyy[k] + 2.0 * tk0_yyy_xy[k] - 2.0 * tk1_yyy_xy[k]);

                    t_yyyy_xyz[k] = fra * t10_yyy_xyz[k] - frc * t11_yyy_xyz[k] + f2t * (3.0 * t20_yy_xyz[k] - 3.0 * t21_yy_xyz[k] + tk0_yyy_xz[k] - tk1_yyy_xz[k]);

                    t_yyyy_xzz[k] = fra * t10_yyy_xzz[k] - frc * t11_yyy_xzz[k] + f2t * (3.0 * t20_yy_xzz[k] - 3.0 * t21_yy_xzz[k]);

                    t_yyyy_yyy[k] = fra * t10_yyy_yyy[k] - frc * t11_yyy_yyy[k] + f2t * (3.0 * t20_yy_yyy[k] - 3.0 * t21_yy_yyy[k] + 3.0 * tk0_yyy_yy[k] - 3.0 * tk1_yyy_yy[k]);

                    t_yyyy_yyz[k] = fra * t10_yyy_yyz[k] - frc * t11_yyy_yyz[k] + f2t * (3.0 * t20_yy_yyz[k] - 3.0 * t21_yy_yyz[k] + 2.0 * tk0_yyy_yz[k] - 2.0 * tk1_yyy_yz[k]);

                    t_yyyy_yzz[k] = fra * t10_yyy_yzz[k] - frc * t11_yyy_yzz[k] + f2t * (3.0 * t20_yy_yzz[k] - 3.0 * t21_yy_yzz[k] + tk0_yyy_zz[k] - tk1_yyy_zz[k]);

                    t_yyyy_zzz[k] = fra * t10_yyy_zzz[k] - frc * t11_yyy_zzz[k] + f2t * (3.0 * t20_yy_zzz[k] - 3.0 * t21_yy_zzz[k]);

                    t_yyyz_xxx[k] = fra * t10_yyz_xxx[k] - frc * t11_yyz_xxx[k] + f2t * (2.0 * t20_yz_xxx[k] - 2.0 * t21_yz_xxx[k]);

                    t_yyyz_xxy[k] = fra * t10_yyz_xxy[k] - frc * t11_yyz_xxy[k] + f2t * (2.0 * t20_yz_xxy[k] - 2.0 * t21_yz_xxy[k] + tk0_yyz_xx[k] - tk1_yyz_xx[k]);

                    t_yyyz_xxz[k] = fra * t10_yyz_xxz[k] - frc * t11_yyz_xxz[k] + f2t * (2.0 * t20_yz_xxz[k] - 2.0 * t21_yz_xxz[k]);

                    t_yyyz_xyy[k] = fra * t10_yyz_xyy[k] - frc * t11_yyz_xyy[k] + f2t * (2.0 * t20_yz_xyy[k] - 2.0 * t21_yz_xyy[k] + 2.0 * tk0_yyz_xy[k] - 2.0 * tk1_yyz_xy[k]);

                    t_yyyz_xyz[k] = fra * t10_yyz_xyz[k] - frc * t11_yyz_xyz[k] + f2t * (2.0 * t20_yz_xyz[k] - 2.0 * t21_yz_xyz[k] + tk0_yyz_xz[k] - tk1_yyz_xz[k]);

                    t_yyyz_xzz[k] = fra * t10_yyz_xzz[k] - frc * t11_yyz_xzz[k] + f2t * (2.0 * t20_yz_xzz[k] - 2.0 * t21_yz_xzz[k]);

                    t_yyyz_yyy[k] = fra * t10_yyz_yyy[k] - frc * t11_yyz_yyy[k] + f2t * (2.0 * t20_yz_yyy[k] - 2.0 * t21_yz_yyy[k] + 3.0 * tk0_yyz_yy[k] - 3.0 * tk1_yyz_yy[k]);

                    t_yyyz_yyz[k] = fra * t10_yyz_yyz[k] - frc * t11_yyz_yyz[k] + f2t * (2.0 * t20_yz_yyz[k] - 2.0 * t21_yz_yyz[k] + 2.0 * tk0_yyz_yz[k] - 2.0 * tk1_yyz_yz[k]);

                    t_yyyz_yzz[k] = fra * t10_yyz_yzz[k] - frc * t11_yyz_yzz[k] + f2t * (2.0 * t20_yz_yzz[k] - 2.0 * t21_yz_yzz[k] + tk0_yyz_zz[k] - tk1_yyz_zz[k]);

                    t_yyyz_zzz[k] = fra * t10_yyz_zzz[k] - frc * t11_yyz_zzz[k] + f2t * (2.0 * t20_yz_zzz[k] - 2.0 * t21_yz_zzz[k]);

                    t_yyzz_xxx[k] = fra * t10_yzz_xxx[k] - frc * t11_yzz_xxx[k] + f2t * (t20_zz_xxx[k] - t21_zz_xxx[k]);

                    t_yyzz_xxy[k] = fra * t10_yzz_xxy[k] - frc * t11_yzz_xxy[k] + f2t * (t20_zz_xxy[k] - t21_zz_xxy[k] + tk0_yzz_xx[k] - tk1_yzz_xx[k]);

                    t_yyzz_xxz[k] = fra * t10_yzz_xxz[k] - frc * t11_yzz_xxz[k] + f2t * (t20_zz_xxz[k] - t21_zz_xxz[k]);

                    t_yyzz_xyy[k] = fra * t10_yzz_xyy[k] - frc * t11_yzz_xyy[k] + f2t * (t20_zz_xyy[k] - t21_zz_xyy[k] + 2.0 * tk0_yzz_xy[k] - 2.0 * tk1_yzz_xy[k]);

                    t_yyzz_xyz[k] = fra * t10_yzz_xyz[k] - frc * t11_yzz_xyz[k] + f2t * (t20_zz_xyz[k] - t21_zz_xyz[k] + tk0_yzz_xz[k] - tk1_yzz_xz[k]);

                    t_yyzz_xzz[k] = fra * t10_yzz_xzz[k] - frc * t11_yzz_xzz[k] + f2t * (t20_zz_xzz[k] - t21_zz_xzz[k]);

                    t_yyzz_yyy[k] = fra * t10_yzz_yyy[k] - frc * t11_yzz_yyy[k] + f2t * (t20_zz_yyy[k] - t21_zz_yyy[k] + 3.0 * tk0_yzz_yy[k] - 3.0 * tk1_yzz_yy[k]);

                    t_yyzz_yyz[k] = fra * t10_yzz_yyz[k] - frc * t11_yzz_yyz[k] + f2t * (t20_zz_yyz[k] - t21_zz_yyz[k] + 2.0 * tk0_yzz_yz[k] - 2.0 * tk1_yzz_yz[k]);

                    t_yyzz_yzz[k] = fra * t10_yzz_yzz[k] - frc * t11_yzz_yzz[k] + f2t * (t20_zz_yzz[k] - t21_zz_yzz[k] + tk0_yzz_zz[k] - tk1_yzz_zz[k]);

                    t_yyzz_zzz[k] = fra * t10_yzz_zzz[k] - frc * t11_yzz_zzz[k] + f2t * (t20_zz_zzz[k] - t21_zz_zzz[k]);

                    t_yzzz_xxx[k] = fra * t10_zzz_xxx[k] - frc * t11_zzz_xxx[k];

                    t_yzzz_xxy[k] = fra * t10_zzz_xxy[k] - frc * t11_zzz_xxy[k] + f2t * (tk0_zzz_xx[k] - tk1_zzz_xx[k]);

                    t_yzzz_xxz[k] = fra * t10_zzz_xxz[k] - frc * t11_zzz_xxz[k];

                    t_yzzz_xyy[k] = fra * t10_zzz_xyy[k] - frc * t11_zzz_xyy[k] + f2t * (2.0 * tk0_zzz_xy[k] - 2.0 * tk1_zzz_xy[k]);

                    t_yzzz_xyz[k] = fra * t10_zzz_xyz[k] - frc * t11_zzz_xyz[k] + f2t * (tk0_zzz_xz[k] - tk1_zzz_xz[k]);

                    t_yzzz_xzz[k] = fra * t10_zzz_xzz[k] - frc * t11_zzz_xzz[k];

                    t_yzzz_yyy[k] = fra * t10_zzz_yyy[k] - frc * t11_zzz_yyy[k] + f2t * (3.0 * tk0_zzz_yy[k] - 3.0 * tk1_zzz_yy[k]);

                    t_yzzz_yyz[k] = fra * t10_zzz_yyz[k] - frc * t11_zzz_yyz[k] + f2t * (2.0 * tk0_zzz_yz[k] - 2.0 * tk1_zzz_yz[k]);

                    t_yzzz_yzz[k] = fra * t10_zzz_yzz[k] - frc * t11_zzz_yzz[k] + f2t * (tk0_zzz_zz[k] - tk1_zzz_zz[k]);

                    t_yzzz_zzz[k] = fra * t10_zzz_zzz[k] - frc * t11_zzz_zzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzzz_xxx[k] = fra * t10_zzz_xxx[k] - frc * t11_zzz_xxx[k] + f2t * (3.0 * t20_zz_xxx[k] - 3.0 * t21_zz_xxx[k]);

                    t_zzzz_xxy[k] = fra * t10_zzz_xxy[k] - frc * t11_zzz_xxy[k] + f2t * (3.0 * t20_zz_xxy[k] - 3.0 * t21_zz_xxy[k]);

                    t_zzzz_xxz[k] = fra * t10_zzz_xxz[k] - frc * t11_zzz_xxz[k] + f2t * (3.0 * t20_zz_xxz[k] - 3.0 * t21_zz_xxz[k] + tk0_zzz_xx[k] - tk1_zzz_xx[k]);

                    t_zzzz_xyy[k] = fra * t10_zzz_xyy[k] - frc * t11_zzz_xyy[k] + f2t * (3.0 * t20_zz_xyy[k] - 3.0 * t21_zz_xyy[k]);

                    t_zzzz_xyz[k] = fra * t10_zzz_xyz[k] - frc * t11_zzz_xyz[k] + f2t * (3.0 * t20_zz_xyz[k] - 3.0 * t21_zz_xyz[k] + tk0_zzz_xy[k] - tk1_zzz_xy[k]);

                    t_zzzz_xzz[k] = fra * t10_zzz_xzz[k] - frc * t11_zzz_xzz[k] + f2t * (3.0 * t20_zz_xzz[k] - 3.0 * t21_zz_xzz[k] + 2.0 * tk0_zzz_xz[k] - 2.0 * tk1_zzz_xz[k]);

                    t_zzzz_yyy[k] = fra * t10_zzz_yyy[k] - frc * t11_zzz_yyy[k] + f2t * (3.0 * t20_zz_yyy[k] - 3.0 * t21_zz_yyy[k]);

                    t_zzzz_yyz[k] = fra * t10_zzz_yyz[k] - frc * t11_zzz_yyz[k] + f2t * (3.0 * t20_zz_yyz[k] - 3.0 * t21_zz_yyz[k] + tk0_zzz_yy[k] - tk1_zzz_yy[k]);

                    t_zzzz_yzz[k] = fra * t10_zzz_yzz[k] - frc * t11_zzz_yzz[k] + f2t * (3.0 * t20_zz_yzz[k] - 3.0 * t21_zz_yzz[k] + 2.0 * tk0_zzz_yz[k] - 2.0 * tk1_zzz_yz[k]);

                    t_zzzz_zzz[k] = fra * t10_zzz_zzz[k] - frc * t11_zzz_zzz[k] + f2t * (3.0 * t20_zz_zzz[k] - 3.0 * t21_zz_zzz[k] + 3.0 * tk0_zzz_zz[k] - 3.0 * tk1_zzz_zz[k]);
                }

                idx++;
            }
        }
    }
    
    void
    compNuclearPotentialForGG(      CMemBlock2D<double>&  primBuffer,
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

        if (!genfunc::isInVector(recPattern, {4, 4, 0})) return;

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute primitive integrals up to required order

        auto bord = genfunc::maxOrderOfPair(recPattern, 4, 4);

        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer

            auto toff   = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {4, 4, i});

            auto t10off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 4, i});

            auto t11off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 4, i + 1});

            auto t20off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 4, i});

            auto t21off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {2, 4, i + 1});

            auto tk0off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 3, i});

            auto tk1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {3, 3, i + 1});

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

                // set up pointers to (F|A(0)|F)^(m) integrals

                auto tk0_xxx_xxx = primBuffer.data(tk0off + 100 * idx);

                auto tk0_xxx_xxy = primBuffer.data(tk0off + 100 * idx + 1);

                auto tk0_xxx_xxz = primBuffer.data(tk0off + 100 * idx + 2);

                auto tk0_xxx_xyy = primBuffer.data(tk0off + 100 * idx + 3);

                auto tk0_xxx_xyz = primBuffer.data(tk0off + 100 * idx + 4);

                auto tk0_xxx_xzz = primBuffer.data(tk0off + 100 * idx + 5);

                auto tk0_xxx_yyy = primBuffer.data(tk0off + 100 * idx + 6);

                auto tk0_xxx_yyz = primBuffer.data(tk0off + 100 * idx + 7);

                auto tk0_xxx_yzz = primBuffer.data(tk0off + 100 * idx + 8);

                auto tk0_xxx_zzz = primBuffer.data(tk0off + 100 * idx + 9);

                auto tk0_xxy_xxx = primBuffer.data(tk0off + 100 * idx + 10);

                auto tk0_xxy_xxy = primBuffer.data(tk0off + 100 * idx + 11);

                auto tk0_xxy_xxz = primBuffer.data(tk0off + 100 * idx + 12);

                auto tk0_xxy_xyy = primBuffer.data(tk0off + 100 * idx + 13);

                auto tk0_xxy_xyz = primBuffer.data(tk0off + 100 * idx + 14);

                auto tk0_xxy_xzz = primBuffer.data(tk0off + 100 * idx + 15);

                auto tk0_xxy_yyy = primBuffer.data(tk0off + 100 * idx + 16);

                auto tk0_xxy_yyz = primBuffer.data(tk0off + 100 * idx + 17);

                auto tk0_xxy_yzz = primBuffer.data(tk0off + 100 * idx + 18);

                auto tk0_xxy_zzz = primBuffer.data(tk0off + 100 * idx + 19);

                auto tk0_xxz_xxx = primBuffer.data(tk0off + 100 * idx + 20);

                auto tk0_xxz_xxy = primBuffer.data(tk0off + 100 * idx + 21);

                auto tk0_xxz_xxz = primBuffer.data(tk0off + 100 * idx + 22);

                auto tk0_xxz_xyy = primBuffer.data(tk0off + 100 * idx + 23);

                auto tk0_xxz_xyz = primBuffer.data(tk0off + 100 * idx + 24);

                auto tk0_xxz_xzz = primBuffer.data(tk0off + 100 * idx + 25);

                auto tk0_xxz_yyy = primBuffer.data(tk0off + 100 * idx + 26);

                auto tk0_xxz_yyz = primBuffer.data(tk0off + 100 * idx + 27);

                auto tk0_xxz_yzz = primBuffer.data(tk0off + 100 * idx + 28);

                auto tk0_xxz_zzz = primBuffer.data(tk0off + 100 * idx + 29);

                auto tk0_xyy_xxx = primBuffer.data(tk0off + 100 * idx + 30);

                auto tk0_xyy_xxy = primBuffer.data(tk0off + 100 * idx + 31);

                auto tk0_xyy_xxz = primBuffer.data(tk0off + 100 * idx + 32);

                auto tk0_xyy_xyy = primBuffer.data(tk0off + 100 * idx + 33);

                auto tk0_xyy_xyz = primBuffer.data(tk0off + 100 * idx + 34);

                auto tk0_xyy_xzz = primBuffer.data(tk0off + 100 * idx + 35);

                auto tk0_xyy_yyy = primBuffer.data(tk0off + 100 * idx + 36);

                auto tk0_xyy_yyz = primBuffer.data(tk0off + 100 * idx + 37);

                auto tk0_xyy_yzz = primBuffer.data(tk0off + 100 * idx + 38);

                auto tk0_xyy_zzz = primBuffer.data(tk0off + 100 * idx + 39);

                auto tk0_xyz_xxx = primBuffer.data(tk0off + 100 * idx + 40);

                auto tk0_xyz_xxy = primBuffer.data(tk0off + 100 * idx + 41);

                auto tk0_xyz_xxz = primBuffer.data(tk0off + 100 * idx + 42);

                auto tk0_xyz_xyy = primBuffer.data(tk0off + 100 * idx + 43);

                auto tk0_xyz_xyz = primBuffer.data(tk0off + 100 * idx + 44);

                auto tk0_xyz_xzz = primBuffer.data(tk0off + 100 * idx + 45);

                auto tk0_xyz_yyy = primBuffer.data(tk0off + 100 * idx + 46);

                auto tk0_xyz_yyz = primBuffer.data(tk0off + 100 * idx + 47);

                auto tk0_xyz_yzz = primBuffer.data(tk0off + 100 * idx + 48);

                auto tk0_xyz_zzz = primBuffer.data(tk0off + 100 * idx + 49);

                auto tk0_xzz_xxx = primBuffer.data(tk0off + 100 * idx + 50);

                auto tk0_xzz_xxy = primBuffer.data(tk0off + 100 * idx + 51);

                auto tk0_xzz_xxz = primBuffer.data(tk0off + 100 * idx + 52);

                auto tk0_xzz_xyy = primBuffer.data(tk0off + 100 * idx + 53);

                auto tk0_xzz_xyz = primBuffer.data(tk0off + 100 * idx + 54);

                auto tk0_xzz_xzz = primBuffer.data(tk0off + 100 * idx + 55);

                auto tk0_xzz_yyy = primBuffer.data(tk0off + 100 * idx + 56);

                auto tk0_xzz_yyz = primBuffer.data(tk0off + 100 * idx + 57);

                auto tk0_xzz_yzz = primBuffer.data(tk0off + 100 * idx + 58);

                auto tk0_xzz_zzz = primBuffer.data(tk0off + 100 * idx + 59);

                auto tk0_yyy_xxx = primBuffer.data(tk0off + 100 * idx + 60);

                auto tk0_yyy_xxy = primBuffer.data(tk0off + 100 * idx + 61);

                auto tk0_yyy_xxz = primBuffer.data(tk0off + 100 * idx + 62);

                auto tk0_yyy_xyy = primBuffer.data(tk0off + 100 * idx + 63);

                auto tk0_yyy_xyz = primBuffer.data(tk0off + 100 * idx + 64);

                auto tk0_yyy_xzz = primBuffer.data(tk0off + 100 * idx + 65);

                auto tk0_yyy_yyy = primBuffer.data(tk0off + 100 * idx + 66);

                auto tk0_yyy_yyz = primBuffer.data(tk0off + 100 * idx + 67);

                auto tk0_yyy_yzz = primBuffer.data(tk0off + 100 * idx + 68);

                auto tk0_yyy_zzz = primBuffer.data(tk0off + 100 * idx + 69);

                auto tk0_yyz_xxx = primBuffer.data(tk0off + 100 * idx + 70);

                auto tk0_yyz_xxy = primBuffer.data(tk0off + 100 * idx + 71);

                auto tk0_yyz_xxz = primBuffer.data(tk0off + 100 * idx + 72);

                auto tk0_yyz_xyy = primBuffer.data(tk0off + 100 * idx + 73);

                auto tk0_yyz_xyz = primBuffer.data(tk0off + 100 * idx + 74);

                auto tk0_yyz_xzz = primBuffer.data(tk0off + 100 * idx + 75);

                auto tk0_yyz_yyy = primBuffer.data(tk0off + 100 * idx + 76);

                auto tk0_yyz_yyz = primBuffer.data(tk0off + 100 * idx + 77);

                auto tk0_yyz_yzz = primBuffer.data(tk0off + 100 * idx + 78);

                auto tk0_yyz_zzz = primBuffer.data(tk0off + 100 * idx + 79);

                auto tk0_yzz_xxx = primBuffer.data(tk0off + 100 * idx + 80);

                auto tk0_yzz_xxy = primBuffer.data(tk0off + 100 * idx + 81);

                auto tk0_yzz_xxz = primBuffer.data(tk0off + 100 * idx + 82);

                auto tk0_yzz_xyy = primBuffer.data(tk0off + 100 * idx + 83);

                auto tk0_yzz_xyz = primBuffer.data(tk0off + 100 * idx + 84);

                auto tk0_yzz_xzz = primBuffer.data(tk0off + 100 * idx + 85);

                auto tk0_yzz_yyy = primBuffer.data(tk0off + 100 * idx + 86);

                auto tk0_yzz_yyz = primBuffer.data(tk0off + 100 * idx + 87);

                auto tk0_yzz_yzz = primBuffer.data(tk0off + 100 * idx + 88);

                auto tk0_yzz_zzz = primBuffer.data(tk0off + 100 * idx + 89);

                auto tk0_zzz_xxx = primBuffer.data(tk0off + 100 * idx + 90);

                auto tk0_zzz_xxy = primBuffer.data(tk0off + 100 * idx + 91);

                auto tk0_zzz_xxz = primBuffer.data(tk0off + 100 * idx + 92);

                auto tk0_zzz_xyy = primBuffer.data(tk0off + 100 * idx + 93);

                auto tk0_zzz_xyz = primBuffer.data(tk0off + 100 * idx + 94);

                auto tk0_zzz_xzz = primBuffer.data(tk0off + 100 * idx + 95);

                auto tk0_zzz_yyy = primBuffer.data(tk0off + 100 * idx + 96);

                auto tk0_zzz_yyz = primBuffer.data(tk0off + 100 * idx + 97);

                auto tk0_zzz_yzz = primBuffer.data(tk0off + 100 * idx + 98);

                auto tk0_zzz_zzz = primBuffer.data(tk0off + 100 * idx + 99);

                // set up pointers to (F|A(0)|F)^(m+1) integrals

                auto tk1_xxx_xxx = primBuffer.data(tk1off + 100 * idx);

                auto tk1_xxx_xxy = primBuffer.data(tk1off + 100 * idx + 1);

                auto tk1_xxx_xxz = primBuffer.data(tk1off + 100 * idx + 2);

                auto tk1_xxx_xyy = primBuffer.data(tk1off + 100 * idx + 3);

                auto tk1_xxx_xyz = primBuffer.data(tk1off + 100 * idx + 4);

                auto tk1_xxx_xzz = primBuffer.data(tk1off + 100 * idx + 5);

                auto tk1_xxx_yyy = primBuffer.data(tk1off + 100 * idx + 6);

                auto tk1_xxx_yyz = primBuffer.data(tk1off + 100 * idx + 7);

                auto tk1_xxx_yzz = primBuffer.data(tk1off + 100 * idx + 8);

                auto tk1_xxx_zzz = primBuffer.data(tk1off + 100 * idx + 9);

                auto tk1_xxy_xxx = primBuffer.data(tk1off + 100 * idx + 10);

                auto tk1_xxy_xxy = primBuffer.data(tk1off + 100 * idx + 11);

                auto tk1_xxy_xxz = primBuffer.data(tk1off + 100 * idx + 12);

                auto tk1_xxy_xyy = primBuffer.data(tk1off + 100 * idx + 13);

                auto tk1_xxy_xyz = primBuffer.data(tk1off + 100 * idx + 14);

                auto tk1_xxy_xzz = primBuffer.data(tk1off + 100 * idx + 15);

                auto tk1_xxy_yyy = primBuffer.data(tk1off + 100 * idx + 16);

                auto tk1_xxy_yyz = primBuffer.data(tk1off + 100 * idx + 17);

                auto tk1_xxy_yzz = primBuffer.data(tk1off + 100 * idx + 18);

                auto tk1_xxy_zzz = primBuffer.data(tk1off + 100 * idx + 19);

                auto tk1_xxz_xxx = primBuffer.data(tk1off + 100 * idx + 20);

                auto tk1_xxz_xxy = primBuffer.data(tk1off + 100 * idx + 21);

                auto tk1_xxz_xxz = primBuffer.data(tk1off + 100 * idx + 22);

                auto tk1_xxz_xyy = primBuffer.data(tk1off + 100 * idx + 23);

                auto tk1_xxz_xyz = primBuffer.data(tk1off + 100 * idx + 24);

                auto tk1_xxz_xzz = primBuffer.data(tk1off + 100 * idx + 25);

                auto tk1_xxz_yyy = primBuffer.data(tk1off + 100 * idx + 26);

                auto tk1_xxz_yyz = primBuffer.data(tk1off + 100 * idx + 27);

                auto tk1_xxz_yzz = primBuffer.data(tk1off + 100 * idx + 28);

                auto tk1_xxz_zzz = primBuffer.data(tk1off + 100 * idx + 29);

                auto tk1_xyy_xxx = primBuffer.data(tk1off + 100 * idx + 30);

                auto tk1_xyy_xxy = primBuffer.data(tk1off + 100 * idx + 31);

                auto tk1_xyy_xxz = primBuffer.data(tk1off + 100 * idx + 32);

                auto tk1_xyy_xyy = primBuffer.data(tk1off + 100 * idx + 33);

                auto tk1_xyy_xyz = primBuffer.data(tk1off + 100 * idx + 34);

                auto tk1_xyy_xzz = primBuffer.data(tk1off + 100 * idx + 35);

                auto tk1_xyy_yyy = primBuffer.data(tk1off + 100 * idx + 36);

                auto tk1_xyy_yyz = primBuffer.data(tk1off + 100 * idx + 37);

                auto tk1_xyy_yzz = primBuffer.data(tk1off + 100 * idx + 38);

                auto tk1_xyy_zzz = primBuffer.data(tk1off + 100 * idx + 39);

                auto tk1_xyz_xxx = primBuffer.data(tk1off + 100 * idx + 40);

                auto tk1_xyz_xxy = primBuffer.data(tk1off + 100 * idx + 41);

                auto tk1_xyz_xxz = primBuffer.data(tk1off + 100 * idx + 42);

                auto tk1_xyz_xyy = primBuffer.data(tk1off + 100 * idx + 43);

                auto tk1_xyz_xyz = primBuffer.data(tk1off + 100 * idx + 44);

                auto tk1_xyz_xzz = primBuffer.data(tk1off + 100 * idx + 45);

                auto tk1_xyz_yyy = primBuffer.data(tk1off + 100 * idx + 46);

                auto tk1_xyz_yyz = primBuffer.data(tk1off + 100 * idx + 47);

                auto tk1_xyz_yzz = primBuffer.data(tk1off + 100 * idx + 48);

                auto tk1_xyz_zzz = primBuffer.data(tk1off + 100 * idx + 49);

                auto tk1_xzz_xxx = primBuffer.data(tk1off + 100 * idx + 50);

                auto tk1_xzz_xxy = primBuffer.data(tk1off + 100 * idx + 51);

                auto tk1_xzz_xxz = primBuffer.data(tk1off + 100 * idx + 52);

                auto tk1_xzz_xyy = primBuffer.data(tk1off + 100 * idx + 53);

                auto tk1_xzz_xyz = primBuffer.data(tk1off + 100 * idx + 54);

                auto tk1_xzz_xzz = primBuffer.data(tk1off + 100 * idx + 55);

                auto tk1_xzz_yyy = primBuffer.data(tk1off + 100 * idx + 56);

                auto tk1_xzz_yyz = primBuffer.data(tk1off + 100 * idx + 57);

                auto tk1_xzz_yzz = primBuffer.data(tk1off + 100 * idx + 58);

                auto tk1_xzz_zzz = primBuffer.data(tk1off + 100 * idx + 59);

                auto tk1_yyy_xxx = primBuffer.data(tk1off + 100 * idx + 60);

                auto tk1_yyy_xxy = primBuffer.data(tk1off + 100 * idx + 61);

                auto tk1_yyy_xxz = primBuffer.data(tk1off + 100 * idx + 62);

                auto tk1_yyy_xyy = primBuffer.data(tk1off + 100 * idx + 63);

                auto tk1_yyy_xyz = primBuffer.data(tk1off + 100 * idx + 64);

                auto tk1_yyy_xzz = primBuffer.data(tk1off + 100 * idx + 65);

                auto tk1_yyy_yyy = primBuffer.data(tk1off + 100 * idx + 66);

                auto tk1_yyy_yyz = primBuffer.data(tk1off + 100 * idx + 67);

                auto tk1_yyy_yzz = primBuffer.data(tk1off + 100 * idx + 68);

                auto tk1_yyy_zzz = primBuffer.data(tk1off + 100 * idx + 69);

                auto tk1_yyz_xxx = primBuffer.data(tk1off + 100 * idx + 70);

                auto tk1_yyz_xxy = primBuffer.data(tk1off + 100 * idx + 71);

                auto tk1_yyz_xxz = primBuffer.data(tk1off + 100 * idx + 72);

                auto tk1_yyz_xyy = primBuffer.data(tk1off + 100 * idx + 73);

                auto tk1_yyz_xyz = primBuffer.data(tk1off + 100 * idx + 74);

                auto tk1_yyz_xzz = primBuffer.data(tk1off + 100 * idx + 75);

                auto tk1_yyz_yyy = primBuffer.data(tk1off + 100 * idx + 76);

                auto tk1_yyz_yyz = primBuffer.data(tk1off + 100 * idx + 77);

                auto tk1_yyz_yzz = primBuffer.data(tk1off + 100 * idx + 78);

                auto tk1_yyz_zzz = primBuffer.data(tk1off + 100 * idx + 79);

                auto tk1_yzz_xxx = primBuffer.data(tk1off + 100 * idx + 80);

                auto tk1_yzz_xxy = primBuffer.data(tk1off + 100 * idx + 81);

                auto tk1_yzz_xxz = primBuffer.data(tk1off + 100 * idx + 82);

                auto tk1_yzz_xyy = primBuffer.data(tk1off + 100 * idx + 83);

                auto tk1_yzz_xyz = primBuffer.data(tk1off + 100 * idx + 84);

                auto tk1_yzz_xzz = primBuffer.data(tk1off + 100 * idx + 85);

                auto tk1_yzz_yyy = primBuffer.data(tk1off + 100 * idx + 86);

                auto tk1_yzz_yyz = primBuffer.data(tk1off + 100 * idx + 87);

                auto tk1_yzz_yzz = primBuffer.data(tk1off + 100 * idx + 88);

                auto tk1_yzz_zzz = primBuffer.data(tk1off + 100 * idx + 89);

                auto tk1_zzz_xxx = primBuffer.data(tk1off + 100 * idx + 90);

                auto tk1_zzz_xxy = primBuffer.data(tk1off + 100 * idx + 91);

                auto tk1_zzz_xxz = primBuffer.data(tk1off + 100 * idx + 92);

                auto tk1_zzz_xyy = primBuffer.data(tk1off + 100 * idx + 93);

                auto tk1_zzz_xyz = primBuffer.data(tk1off + 100 * idx + 94);

                auto tk1_zzz_xzz = primBuffer.data(tk1off + 100 * idx + 95);

                auto tk1_zzz_yyy = primBuffer.data(tk1off + 100 * idx + 96);

                auto tk1_zzz_yyz = primBuffer.data(tk1off + 100 * idx + 97);

                auto tk1_zzz_yzz = primBuffer.data(tk1off + 100 * idx + 98);

                auto tk1_zzz_zzz = primBuffer.data(tk1off + 100 * idx + 99);

                // set up pointers to (D|A(0)|G)^(m) integrals

                auto t20_xx_xxxx = primBuffer.data(t20off + 90 * idx);

                auto t20_xx_xxxy = primBuffer.data(t20off + 90 * idx + 1);

                auto t20_xx_xxxz = primBuffer.data(t20off + 90 * idx + 2);

                auto t20_xx_xxyy = primBuffer.data(t20off + 90 * idx + 3);

                auto t20_xx_xxyz = primBuffer.data(t20off + 90 * idx + 4);

                auto t20_xx_xxzz = primBuffer.data(t20off + 90 * idx + 5);

                auto t20_xx_xyyy = primBuffer.data(t20off + 90 * idx + 6);

                auto t20_xx_xyyz = primBuffer.data(t20off + 90 * idx + 7);

                auto t20_xx_xyzz = primBuffer.data(t20off + 90 * idx + 8);

                auto t20_xx_xzzz = primBuffer.data(t20off + 90 * idx + 9);

                auto t20_xx_yyyy = primBuffer.data(t20off + 90 * idx + 10);

                auto t20_xx_yyyz = primBuffer.data(t20off + 90 * idx + 11);

                auto t20_xx_yyzz = primBuffer.data(t20off + 90 * idx + 12);

                auto t20_xx_yzzz = primBuffer.data(t20off + 90 * idx + 13);

                auto t20_xx_zzzz = primBuffer.data(t20off + 90 * idx + 14);

                auto t20_xy_xxxx = primBuffer.data(t20off + 90 * idx + 15);

                auto t20_xy_xxxy = primBuffer.data(t20off + 90 * idx + 16);

                auto t20_xy_xxxz = primBuffer.data(t20off + 90 * idx + 17);

                auto t20_xy_xxyy = primBuffer.data(t20off + 90 * idx + 18);

                auto t20_xy_xxyz = primBuffer.data(t20off + 90 * idx + 19);

                auto t20_xy_xxzz = primBuffer.data(t20off + 90 * idx + 20);

                auto t20_xy_xyyy = primBuffer.data(t20off + 90 * idx + 21);

                auto t20_xy_xyyz = primBuffer.data(t20off + 90 * idx + 22);

                auto t20_xy_xyzz = primBuffer.data(t20off + 90 * idx + 23);

                auto t20_xy_xzzz = primBuffer.data(t20off + 90 * idx + 24);

                auto t20_xy_yyyy = primBuffer.data(t20off + 90 * idx + 25);

                auto t20_xy_yyyz = primBuffer.data(t20off + 90 * idx + 26);

                auto t20_xy_yyzz = primBuffer.data(t20off + 90 * idx + 27);

                auto t20_xy_yzzz = primBuffer.data(t20off + 90 * idx + 28);

                auto t20_xy_zzzz = primBuffer.data(t20off + 90 * idx + 29);

                auto t20_xz_xxxx = primBuffer.data(t20off + 90 * idx + 30);

                auto t20_xz_xxxy = primBuffer.data(t20off + 90 * idx + 31);

                auto t20_xz_xxxz = primBuffer.data(t20off + 90 * idx + 32);

                auto t20_xz_xxyy = primBuffer.data(t20off + 90 * idx + 33);

                auto t20_xz_xxyz = primBuffer.data(t20off + 90 * idx + 34);

                auto t20_xz_xxzz = primBuffer.data(t20off + 90 * idx + 35);

                auto t20_xz_xyyy = primBuffer.data(t20off + 90 * idx + 36);

                auto t20_xz_xyyz = primBuffer.data(t20off + 90 * idx + 37);

                auto t20_xz_xyzz = primBuffer.data(t20off + 90 * idx + 38);

                auto t20_xz_xzzz = primBuffer.data(t20off + 90 * idx + 39);

                auto t20_xz_yyyy = primBuffer.data(t20off + 90 * idx + 40);

                auto t20_xz_yyyz = primBuffer.data(t20off + 90 * idx + 41);

                auto t20_xz_yyzz = primBuffer.data(t20off + 90 * idx + 42);

                auto t20_xz_yzzz = primBuffer.data(t20off + 90 * idx + 43);

                auto t20_xz_zzzz = primBuffer.data(t20off + 90 * idx + 44);

                auto t20_yy_xxxx = primBuffer.data(t20off + 90 * idx + 45);

                auto t20_yy_xxxy = primBuffer.data(t20off + 90 * idx + 46);

                auto t20_yy_xxxz = primBuffer.data(t20off + 90 * idx + 47);

                auto t20_yy_xxyy = primBuffer.data(t20off + 90 * idx + 48);

                auto t20_yy_xxyz = primBuffer.data(t20off + 90 * idx + 49);

                auto t20_yy_xxzz = primBuffer.data(t20off + 90 * idx + 50);

                auto t20_yy_xyyy = primBuffer.data(t20off + 90 * idx + 51);

                auto t20_yy_xyyz = primBuffer.data(t20off + 90 * idx + 52);

                auto t20_yy_xyzz = primBuffer.data(t20off + 90 * idx + 53);

                auto t20_yy_xzzz = primBuffer.data(t20off + 90 * idx + 54);

                auto t20_yy_yyyy = primBuffer.data(t20off + 90 * idx + 55);

                auto t20_yy_yyyz = primBuffer.data(t20off + 90 * idx + 56);

                auto t20_yy_yyzz = primBuffer.data(t20off + 90 * idx + 57);

                auto t20_yy_yzzz = primBuffer.data(t20off + 90 * idx + 58);

                auto t20_yy_zzzz = primBuffer.data(t20off + 90 * idx + 59);

                auto t20_yz_xxxx = primBuffer.data(t20off + 90 * idx + 60);

                auto t20_yz_xxxy = primBuffer.data(t20off + 90 * idx + 61);

                auto t20_yz_xxxz = primBuffer.data(t20off + 90 * idx + 62);

                auto t20_yz_xxyy = primBuffer.data(t20off + 90 * idx + 63);

                auto t20_yz_xxyz = primBuffer.data(t20off + 90 * idx + 64);

                auto t20_yz_xxzz = primBuffer.data(t20off + 90 * idx + 65);

                auto t20_yz_xyyy = primBuffer.data(t20off + 90 * idx + 66);

                auto t20_yz_xyyz = primBuffer.data(t20off + 90 * idx + 67);

                auto t20_yz_xyzz = primBuffer.data(t20off + 90 * idx + 68);

                auto t20_yz_xzzz = primBuffer.data(t20off + 90 * idx + 69);

                auto t20_yz_yyyy = primBuffer.data(t20off + 90 * idx + 70);

                auto t20_yz_yyyz = primBuffer.data(t20off + 90 * idx + 71);

                auto t20_yz_yyzz = primBuffer.data(t20off + 90 * idx + 72);

                auto t20_yz_yzzz = primBuffer.data(t20off + 90 * idx + 73);

                auto t20_yz_zzzz = primBuffer.data(t20off + 90 * idx + 74);

                auto t20_zz_xxxx = primBuffer.data(t20off + 90 * idx + 75);

                auto t20_zz_xxxy = primBuffer.data(t20off + 90 * idx + 76);

                auto t20_zz_xxxz = primBuffer.data(t20off + 90 * idx + 77);

                auto t20_zz_xxyy = primBuffer.data(t20off + 90 * idx + 78);

                auto t20_zz_xxyz = primBuffer.data(t20off + 90 * idx + 79);

                auto t20_zz_xxzz = primBuffer.data(t20off + 90 * idx + 80);

                auto t20_zz_xyyy = primBuffer.data(t20off + 90 * idx + 81);

                auto t20_zz_xyyz = primBuffer.data(t20off + 90 * idx + 82);

                auto t20_zz_xyzz = primBuffer.data(t20off + 90 * idx + 83);

                auto t20_zz_xzzz = primBuffer.data(t20off + 90 * idx + 84);

                auto t20_zz_yyyy = primBuffer.data(t20off + 90 * idx + 85);

                auto t20_zz_yyyz = primBuffer.data(t20off + 90 * idx + 86);

                auto t20_zz_yyzz = primBuffer.data(t20off + 90 * idx + 87);

                auto t20_zz_yzzz = primBuffer.data(t20off + 90 * idx + 88);

                auto t20_zz_zzzz = primBuffer.data(t20off + 90 * idx + 89);

                // set up pointers to (D|A(0)|G)^(m+1) integrals

                auto t21_xx_xxxx = primBuffer.data(t21off + 90 * idx);

                auto t21_xx_xxxy = primBuffer.data(t21off + 90 * idx + 1);

                auto t21_xx_xxxz = primBuffer.data(t21off + 90 * idx + 2);

                auto t21_xx_xxyy = primBuffer.data(t21off + 90 * idx + 3);

                auto t21_xx_xxyz = primBuffer.data(t21off + 90 * idx + 4);

                auto t21_xx_xxzz = primBuffer.data(t21off + 90 * idx + 5);

                auto t21_xx_xyyy = primBuffer.data(t21off + 90 * idx + 6);

                auto t21_xx_xyyz = primBuffer.data(t21off + 90 * idx + 7);

                auto t21_xx_xyzz = primBuffer.data(t21off + 90 * idx + 8);

                auto t21_xx_xzzz = primBuffer.data(t21off + 90 * idx + 9);

                auto t21_xx_yyyy = primBuffer.data(t21off + 90 * idx + 10);

                auto t21_xx_yyyz = primBuffer.data(t21off + 90 * idx + 11);

                auto t21_xx_yyzz = primBuffer.data(t21off + 90 * idx + 12);

                auto t21_xx_yzzz = primBuffer.data(t21off + 90 * idx + 13);

                auto t21_xx_zzzz = primBuffer.data(t21off + 90 * idx + 14);

                auto t21_xy_xxxx = primBuffer.data(t21off + 90 * idx + 15);

                auto t21_xy_xxxy = primBuffer.data(t21off + 90 * idx + 16);

                auto t21_xy_xxxz = primBuffer.data(t21off + 90 * idx + 17);

                auto t21_xy_xxyy = primBuffer.data(t21off + 90 * idx + 18);

                auto t21_xy_xxyz = primBuffer.data(t21off + 90 * idx + 19);

                auto t21_xy_xxzz = primBuffer.data(t21off + 90 * idx + 20);

                auto t21_xy_xyyy = primBuffer.data(t21off + 90 * idx + 21);

                auto t21_xy_xyyz = primBuffer.data(t21off + 90 * idx + 22);

                auto t21_xy_xyzz = primBuffer.data(t21off + 90 * idx + 23);

                auto t21_xy_xzzz = primBuffer.data(t21off + 90 * idx + 24);

                auto t21_xy_yyyy = primBuffer.data(t21off + 90 * idx + 25);

                auto t21_xy_yyyz = primBuffer.data(t21off + 90 * idx + 26);

                auto t21_xy_yyzz = primBuffer.data(t21off + 90 * idx + 27);

                auto t21_xy_yzzz = primBuffer.data(t21off + 90 * idx + 28);

                auto t21_xy_zzzz = primBuffer.data(t21off + 90 * idx + 29);

                auto t21_xz_xxxx = primBuffer.data(t21off + 90 * idx + 30);

                auto t21_xz_xxxy = primBuffer.data(t21off + 90 * idx + 31);

                auto t21_xz_xxxz = primBuffer.data(t21off + 90 * idx + 32);

                auto t21_xz_xxyy = primBuffer.data(t21off + 90 * idx + 33);

                auto t21_xz_xxyz = primBuffer.data(t21off + 90 * idx + 34);

                auto t21_xz_xxzz = primBuffer.data(t21off + 90 * idx + 35);

                auto t21_xz_xyyy = primBuffer.data(t21off + 90 * idx + 36);

                auto t21_xz_xyyz = primBuffer.data(t21off + 90 * idx + 37);

                auto t21_xz_xyzz = primBuffer.data(t21off + 90 * idx + 38);

                auto t21_xz_xzzz = primBuffer.data(t21off + 90 * idx + 39);

                auto t21_xz_yyyy = primBuffer.data(t21off + 90 * idx + 40);

                auto t21_xz_yyyz = primBuffer.data(t21off + 90 * idx + 41);

                auto t21_xz_yyzz = primBuffer.data(t21off + 90 * idx + 42);

                auto t21_xz_yzzz = primBuffer.data(t21off + 90 * idx + 43);

                auto t21_xz_zzzz = primBuffer.data(t21off + 90 * idx + 44);

                auto t21_yy_xxxx = primBuffer.data(t21off + 90 * idx + 45);

                auto t21_yy_xxxy = primBuffer.data(t21off + 90 * idx + 46);

                auto t21_yy_xxxz = primBuffer.data(t21off + 90 * idx + 47);

                auto t21_yy_xxyy = primBuffer.data(t21off + 90 * idx + 48);

                auto t21_yy_xxyz = primBuffer.data(t21off + 90 * idx + 49);

                auto t21_yy_xxzz = primBuffer.data(t21off + 90 * idx + 50);

                auto t21_yy_xyyy = primBuffer.data(t21off + 90 * idx + 51);

                auto t21_yy_xyyz = primBuffer.data(t21off + 90 * idx + 52);

                auto t21_yy_xyzz = primBuffer.data(t21off + 90 * idx + 53);

                auto t21_yy_xzzz = primBuffer.data(t21off + 90 * idx + 54);

                auto t21_yy_yyyy = primBuffer.data(t21off + 90 * idx + 55);

                auto t21_yy_yyyz = primBuffer.data(t21off + 90 * idx + 56);

                auto t21_yy_yyzz = primBuffer.data(t21off + 90 * idx + 57);

                auto t21_yy_yzzz = primBuffer.data(t21off + 90 * idx + 58);

                auto t21_yy_zzzz = primBuffer.data(t21off + 90 * idx + 59);

                auto t21_yz_xxxx = primBuffer.data(t21off + 90 * idx + 60);

                auto t21_yz_xxxy = primBuffer.data(t21off + 90 * idx + 61);

                auto t21_yz_xxxz = primBuffer.data(t21off + 90 * idx + 62);

                auto t21_yz_xxyy = primBuffer.data(t21off + 90 * idx + 63);

                auto t21_yz_xxyz = primBuffer.data(t21off + 90 * idx + 64);

                auto t21_yz_xxzz = primBuffer.data(t21off + 90 * idx + 65);

                auto t21_yz_xyyy = primBuffer.data(t21off + 90 * idx + 66);

                auto t21_yz_xyyz = primBuffer.data(t21off + 90 * idx + 67);

                auto t21_yz_xyzz = primBuffer.data(t21off + 90 * idx + 68);

                auto t21_yz_xzzz = primBuffer.data(t21off + 90 * idx + 69);

                auto t21_yz_yyyy = primBuffer.data(t21off + 90 * idx + 70);

                auto t21_yz_yyyz = primBuffer.data(t21off + 90 * idx + 71);

                auto t21_yz_yyzz = primBuffer.data(t21off + 90 * idx + 72);

                auto t21_yz_yzzz = primBuffer.data(t21off + 90 * idx + 73);

                auto t21_yz_zzzz = primBuffer.data(t21off + 90 * idx + 74);

                auto t21_zz_xxxx = primBuffer.data(t21off + 90 * idx + 75);

                auto t21_zz_xxxy = primBuffer.data(t21off + 90 * idx + 76);

                auto t21_zz_xxxz = primBuffer.data(t21off + 90 * idx + 77);

                auto t21_zz_xxyy = primBuffer.data(t21off + 90 * idx + 78);

                auto t21_zz_xxyz = primBuffer.data(t21off + 90 * idx + 79);

                auto t21_zz_xxzz = primBuffer.data(t21off + 90 * idx + 80);

                auto t21_zz_xyyy = primBuffer.data(t21off + 90 * idx + 81);

                auto t21_zz_xyyz = primBuffer.data(t21off + 90 * idx + 82);

                auto t21_zz_xyzz = primBuffer.data(t21off + 90 * idx + 83);

                auto t21_zz_xzzz = primBuffer.data(t21off + 90 * idx + 84);

                auto t21_zz_yyyy = primBuffer.data(t21off + 90 * idx + 85);

                auto t21_zz_yyyz = primBuffer.data(t21off + 90 * idx + 86);

                auto t21_zz_yyzz = primBuffer.data(t21off + 90 * idx + 87);

                auto t21_zz_yzzz = primBuffer.data(t21off + 90 * idx + 88);

                auto t21_zz_zzzz = primBuffer.data(t21off + 90 * idx + 89);

                // set up pointers to (F|A(0)|G)^(m) integrals

                auto t10_xxx_xxxx = primBuffer.data(t10off + 150 * idx);

                auto t10_xxx_xxxy = primBuffer.data(t10off + 150 * idx + 1);

                auto t10_xxx_xxxz = primBuffer.data(t10off + 150 * idx + 2);

                auto t10_xxx_xxyy = primBuffer.data(t10off + 150 * idx + 3);

                auto t10_xxx_xxyz = primBuffer.data(t10off + 150 * idx + 4);

                auto t10_xxx_xxzz = primBuffer.data(t10off + 150 * idx + 5);

                auto t10_xxx_xyyy = primBuffer.data(t10off + 150 * idx + 6);

                auto t10_xxx_xyyz = primBuffer.data(t10off + 150 * idx + 7);

                auto t10_xxx_xyzz = primBuffer.data(t10off + 150 * idx + 8);

                auto t10_xxx_xzzz = primBuffer.data(t10off + 150 * idx + 9);

                auto t10_xxx_yyyy = primBuffer.data(t10off + 150 * idx + 10);

                auto t10_xxx_yyyz = primBuffer.data(t10off + 150 * idx + 11);

                auto t10_xxx_yyzz = primBuffer.data(t10off + 150 * idx + 12);

                auto t10_xxx_yzzz = primBuffer.data(t10off + 150 * idx + 13);

                auto t10_xxx_zzzz = primBuffer.data(t10off + 150 * idx + 14);

                auto t10_xxy_xxxx = primBuffer.data(t10off + 150 * idx + 15);

                auto t10_xxy_xxxy = primBuffer.data(t10off + 150 * idx + 16);

                auto t10_xxy_xxxz = primBuffer.data(t10off + 150 * idx + 17);

                auto t10_xxy_xxyy = primBuffer.data(t10off + 150 * idx + 18);

                auto t10_xxy_xxyz = primBuffer.data(t10off + 150 * idx + 19);

                auto t10_xxy_xxzz = primBuffer.data(t10off + 150 * idx + 20);

                auto t10_xxy_xyyy = primBuffer.data(t10off + 150 * idx + 21);

                auto t10_xxy_xyyz = primBuffer.data(t10off + 150 * idx + 22);

                auto t10_xxy_xyzz = primBuffer.data(t10off + 150 * idx + 23);

                auto t10_xxy_xzzz = primBuffer.data(t10off + 150 * idx + 24);

                auto t10_xxy_yyyy = primBuffer.data(t10off + 150 * idx + 25);

                auto t10_xxy_yyyz = primBuffer.data(t10off + 150 * idx + 26);

                auto t10_xxy_yyzz = primBuffer.data(t10off + 150 * idx + 27);

                auto t10_xxy_yzzz = primBuffer.data(t10off + 150 * idx + 28);

                auto t10_xxy_zzzz = primBuffer.data(t10off + 150 * idx + 29);

                auto t10_xxz_xxxx = primBuffer.data(t10off + 150 * idx + 30);

                auto t10_xxz_xxxy = primBuffer.data(t10off + 150 * idx + 31);

                auto t10_xxz_xxxz = primBuffer.data(t10off + 150 * idx + 32);

                auto t10_xxz_xxyy = primBuffer.data(t10off + 150 * idx + 33);

                auto t10_xxz_xxyz = primBuffer.data(t10off + 150 * idx + 34);

                auto t10_xxz_xxzz = primBuffer.data(t10off + 150 * idx + 35);

                auto t10_xxz_xyyy = primBuffer.data(t10off + 150 * idx + 36);

                auto t10_xxz_xyyz = primBuffer.data(t10off + 150 * idx + 37);

                auto t10_xxz_xyzz = primBuffer.data(t10off + 150 * idx + 38);

                auto t10_xxz_xzzz = primBuffer.data(t10off + 150 * idx + 39);

                auto t10_xxz_yyyy = primBuffer.data(t10off + 150 * idx + 40);

                auto t10_xxz_yyyz = primBuffer.data(t10off + 150 * idx + 41);

                auto t10_xxz_yyzz = primBuffer.data(t10off + 150 * idx + 42);

                auto t10_xxz_yzzz = primBuffer.data(t10off + 150 * idx + 43);

                auto t10_xxz_zzzz = primBuffer.data(t10off + 150 * idx + 44);

                auto t10_xyy_xxxx = primBuffer.data(t10off + 150 * idx + 45);

                auto t10_xyy_xxxy = primBuffer.data(t10off + 150 * idx + 46);

                auto t10_xyy_xxxz = primBuffer.data(t10off + 150 * idx + 47);

                auto t10_xyy_xxyy = primBuffer.data(t10off + 150 * idx + 48);

                auto t10_xyy_xxyz = primBuffer.data(t10off + 150 * idx + 49);

                auto t10_xyy_xxzz = primBuffer.data(t10off + 150 * idx + 50);

                auto t10_xyy_xyyy = primBuffer.data(t10off + 150 * idx + 51);

                auto t10_xyy_xyyz = primBuffer.data(t10off + 150 * idx + 52);

                auto t10_xyy_xyzz = primBuffer.data(t10off + 150 * idx + 53);

                auto t10_xyy_xzzz = primBuffer.data(t10off + 150 * idx + 54);

                auto t10_xyy_yyyy = primBuffer.data(t10off + 150 * idx + 55);

                auto t10_xyy_yyyz = primBuffer.data(t10off + 150 * idx + 56);

                auto t10_xyy_yyzz = primBuffer.data(t10off + 150 * idx + 57);

                auto t10_xyy_yzzz = primBuffer.data(t10off + 150 * idx + 58);

                auto t10_xyy_zzzz = primBuffer.data(t10off + 150 * idx + 59);

                auto t10_xyz_xxxx = primBuffer.data(t10off + 150 * idx + 60);

                auto t10_xyz_xxxy = primBuffer.data(t10off + 150 * idx + 61);

                auto t10_xyz_xxxz = primBuffer.data(t10off + 150 * idx + 62);

                auto t10_xyz_xxyy = primBuffer.data(t10off + 150 * idx + 63);

                auto t10_xyz_xxyz = primBuffer.data(t10off + 150 * idx + 64);

                auto t10_xyz_xxzz = primBuffer.data(t10off + 150 * idx + 65);

                auto t10_xyz_xyyy = primBuffer.data(t10off + 150 * idx + 66);

                auto t10_xyz_xyyz = primBuffer.data(t10off + 150 * idx + 67);

                auto t10_xyz_xyzz = primBuffer.data(t10off + 150 * idx + 68);

                auto t10_xyz_xzzz = primBuffer.data(t10off + 150 * idx + 69);

                auto t10_xyz_yyyy = primBuffer.data(t10off + 150 * idx + 70);

                auto t10_xyz_yyyz = primBuffer.data(t10off + 150 * idx + 71);

                auto t10_xyz_yyzz = primBuffer.data(t10off + 150 * idx + 72);

                auto t10_xyz_yzzz = primBuffer.data(t10off + 150 * idx + 73);

                auto t10_xyz_zzzz = primBuffer.data(t10off + 150 * idx + 74);

                auto t10_xzz_xxxx = primBuffer.data(t10off + 150 * idx + 75);

                auto t10_xzz_xxxy = primBuffer.data(t10off + 150 * idx + 76);

                auto t10_xzz_xxxz = primBuffer.data(t10off + 150 * idx + 77);

                auto t10_xzz_xxyy = primBuffer.data(t10off + 150 * idx + 78);

                auto t10_xzz_xxyz = primBuffer.data(t10off + 150 * idx + 79);

                auto t10_xzz_xxzz = primBuffer.data(t10off + 150 * idx + 80);

                auto t10_xzz_xyyy = primBuffer.data(t10off + 150 * idx + 81);

                auto t10_xzz_xyyz = primBuffer.data(t10off + 150 * idx + 82);

                auto t10_xzz_xyzz = primBuffer.data(t10off + 150 * idx + 83);

                auto t10_xzz_xzzz = primBuffer.data(t10off + 150 * idx + 84);

                auto t10_xzz_yyyy = primBuffer.data(t10off + 150 * idx + 85);

                auto t10_xzz_yyyz = primBuffer.data(t10off + 150 * idx + 86);

                auto t10_xzz_yyzz = primBuffer.data(t10off + 150 * idx + 87);

                auto t10_xzz_yzzz = primBuffer.data(t10off + 150 * idx + 88);

                auto t10_xzz_zzzz = primBuffer.data(t10off + 150 * idx + 89);

                auto t10_yyy_xxxx = primBuffer.data(t10off + 150 * idx + 90);

                auto t10_yyy_xxxy = primBuffer.data(t10off + 150 * idx + 91);

                auto t10_yyy_xxxz = primBuffer.data(t10off + 150 * idx + 92);

                auto t10_yyy_xxyy = primBuffer.data(t10off + 150 * idx + 93);

                auto t10_yyy_xxyz = primBuffer.data(t10off + 150 * idx + 94);

                auto t10_yyy_xxzz = primBuffer.data(t10off + 150 * idx + 95);

                auto t10_yyy_xyyy = primBuffer.data(t10off + 150 * idx + 96);

                auto t10_yyy_xyyz = primBuffer.data(t10off + 150 * idx + 97);

                auto t10_yyy_xyzz = primBuffer.data(t10off + 150 * idx + 98);

                auto t10_yyy_xzzz = primBuffer.data(t10off + 150 * idx + 99);

                auto t10_yyy_yyyy = primBuffer.data(t10off + 150 * idx + 100);

                auto t10_yyy_yyyz = primBuffer.data(t10off + 150 * idx + 101);

                auto t10_yyy_yyzz = primBuffer.data(t10off + 150 * idx + 102);

                auto t10_yyy_yzzz = primBuffer.data(t10off + 150 * idx + 103);

                auto t10_yyy_zzzz = primBuffer.data(t10off + 150 * idx + 104);

                auto t10_yyz_xxxx = primBuffer.data(t10off + 150 * idx + 105);

                auto t10_yyz_xxxy = primBuffer.data(t10off + 150 * idx + 106);

                auto t10_yyz_xxxz = primBuffer.data(t10off + 150 * idx + 107);

                auto t10_yyz_xxyy = primBuffer.data(t10off + 150 * idx + 108);

                auto t10_yyz_xxyz = primBuffer.data(t10off + 150 * idx + 109);

                auto t10_yyz_xxzz = primBuffer.data(t10off + 150 * idx + 110);

                auto t10_yyz_xyyy = primBuffer.data(t10off + 150 * idx + 111);

                auto t10_yyz_xyyz = primBuffer.data(t10off + 150 * idx + 112);

                auto t10_yyz_xyzz = primBuffer.data(t10off + 150 * idx + 113);

                auto t10_yyz_xzzz = primBuffer.data(t10off + 150 * idx + 114);

                auto t10_yyz_yyyy = primBuffer.data(t10off + 150 * idx + 115);

                auto t10_yyz_yyyz = primBuffer.data(t10off + 150 * idx + 116);

                auto t10_yyz_yyzz = primBuffer.data(t10off + 150 * idx + 117);

                auto t10_yyz_yzzz = primBuffer.data(t10off + 150 * idx + 118);

                auto t10_yyz_zzzz = primBuffer.data(t10off + 150 * idx + 119);

                auto t10_yzz_xxxx = primBuffer.data(t10off + 150 * idx + 120);

                auto t10_yzz_xxxy = primBuffer.data(t10off + 150 * idx + 121);

                auto t10_yzz_xxxz = primBuffer.data(t10off + 150 * idx + 122);

                auto t10_yzz_xxyy = primBuffer.data(t10off + 150 * idx + 123);

                auto t10_yzz_xxyz = primBuffer.data(t10off + 150 * idx + 124);

                auto t10_yzz_xxzz = primBuffer.data(t10off + 150 * idx + 125);

                auto t10_yzz_xyyy = primBuffer.data(t10off + 150 * idx + 126);

                auto t10_yzz_xyyz = primBuffer.data(t10off + 150 * idx + 127);

                auto t10_yzz_xyzz = primBuffer.data(t10off + 150 * idx + 128);

                auto t10_yzz_xzzz = primBuffer.data(t10off + 150 * idx + 129);

                auto t10_yzz_yyyy = primBuffer.data(t10off + 150 * idx + 130);

                auto t10_yzz_yyyz = primBuffer.data(t10off + 150 * idx + 131);

                auto t10_yzz_yyzz = primBuffer.data(t10off + 150 * idx + 132);

                auto t10_yzz_yzzz = primBuffer.data(t10off + 150 * idx + 133);

                auto t10_yzz_zzzz = primBuffer.data(t10off + 150 * idx + 134);

                auto t10_zzz_xxxx = primBuffer.data(t10off + 150 * idx + 135);

                auto t10_zzz_xxxy = primBuffer.data(t10off + 150 * idx + 136);

                auto t10_zzz_xxxz = primBuffer.data(t10off + 150 * idx + 137);

                auto t10_zzz_xxyy = primBuffer.data(t10off + 150 * idx + 138);

                auto t10_zzz_xxyz = primBuffer.data(t10off + 150 * idx + 139);

                auto t10_zzz_xxzz = primBuffer.data(t10off + 150 * idx + 140);

                auto t10_zzz_xyyy = primBuffer.data(t10off + 150 * idx + 141);

                auto t10_zzz_xyyz = primBuffer.data(t10off + 150 * idx + 142);

                auto t10_zzz_xyzz = primBuffer.data(t10off + 150 * idx + 143);

                auto t10_zzz_xzzz = primBuffer.data(t10off + 150 * idx + 144);

                auto t10_zzz_yyyy = primBuffer.data(t10off + 150 * idx + 145);

                auto t10_zzz_yyyz = primBuffer.data(t10off + 150 * idx + 146);

                auto t10_zzz_yyzz = primBuffer.data(t10off + 150 * idx + 147);

                auto t10_zzz_yzzz = primBuffer.data(t10off + 150 * idx + 148);

                auto t10_zzz_zzzz = primBuffer.data(t10off + 150 * idx + 149);

                // set up pointers to (F|A(0)|G)^(m+1) integrals

                auto t11_xxx_xxxx = primBuffer.data(t11off + 150 * idx);

                auto t11_xxx_xxxy = primBuffer.data(t11off + 150 * idx + 1);

                auto t11_xxx_xxxz = primBuffer.data(t11off + 150 * idx + 2);

                auto t11_xxx_xxyy = primBuffer.data(t11off + 150 * idx + 3);

                auto t11_xxx_xxyz = primBuffer.data(t11off + 150 * idx + 4);

                auto t11_xxx_xxzz = primBuffer.data(t11off + 150 * idx + 5);

                auto t11_xxx_xyyy = primBuffer.data(t11off + 150 * idx + 6);

                auto t11_xxx_xyyz = primBuffer.data(t11off + 150 * idx + 7);

                auto t11_xxx_xyzz = primBuffer.data(t11off + 150 * idx + 8);

                auto t11_xxx_xzzz = primBuffer.data(t11off + 150 * idx + 9);

                auto t11_xxx_yyyy = primBuffer.data(t11off + 150 * idx + 10);

                auto t11_xxx_yyyz = primBuffer.data(t11off + 150 * idx + 11);

                auto t11_xxx_yyzz = primBuffer.data(t11off + 150 * idx + 12);

                auto t11_xxx_yzzz = primBuffer.data(t11off + 150 * idx + 13);

                auto t11_xxx_zzzz = primBuffer.data(t11off + 150 * idx + 14);

                auto t11_xxy_xxxx = primBuffer.data(t11off + 150 * idx + 15);

                auto t11_xxy_xxxy = primBuffer.data(t11off + 150 * idx + 16);

                auto t11_xxy_xxxz = primBuffer.data(t11off + 150 * idx + 17);

                auto t11_xxy_xxyy = primBuffer.data(t11off + 150 * idx + 18);

                auto t11_xxy_xxyz = primBuffer.data(t11off + 150 * idx + 19);

                auto t11_xxy_xxzz = primBuffer.data(t11off + 150 * idx + 20);

                auto t11_xxy_xyyy = primBuffer.data(t11off + 150 * idx + 21);

                auto t11_xxy_xyyz = primBuffer.data(t11off + 150 * idx + 22);

                auto t11_xxy_xyzz = primBuffer.data(t11off + 150 * idx + 23);

                auto t11_xxy_xzzz = primBuffer.data(t11off + 150 * idx + 24);

                auto t11_xxy_yyyy = primBuffer.data(t11off + 150 * idx + 25);

                auto t11_xxy_yyyz = primBuffer.data(t11off + 150 * idx + 26);

                auto t11_xxy_yyzz = primBuffer.data(t11off + 150 * idx + 27);

                auto t11_xxy_yzzz = primBuffer.data(t11off + 150 * idx + 28);

                auto t11_xxy_zzzz = primBuffer.data(t11off + 150 * idx + 29);

                auto t11_xxz_xxxx = primBuffer.data(t11off + 150 * idx + 30);

                auto t11_xxz_xxxy = primBuffer.data(t11off + 150 * idx + 31);

                auto t11_xxz_xxxz = primBuffer.data(t11off + 150 * idx + 32);

                auto t11_xxz_xxyy = primBuffer.data(t11off + 150 * idx + 33);

                auto t11_xxz_xxyz = primBuffer.data(t11off + 150 * idx + 34);

                auto t11_xxz_xxzz = primBuffer.data(t11off + 150 * idx + 35);

                auto t11_xxz_xyyy = primBuffer.data(t11off + 150 * idx + 36);

                auto t11_xxz_xyyz = primBuffer.data(t11off + 150 * idx + 37);

                auto t11_xxz_xyzz = primBuffer.data(t11off + 150 * idx + 38);

                auto t11_xxz_xzzz = primBuffer.data(t11off + 150 * idx + 39);

                auto t11_xxz_yyyy = primBuffer.data(t11off + 150 * idx + 40);

                auto t11_xxz_yyyz = primBuffer.data(t11off + 150 * idx + 41);

                auto t11_xxz_yyzz = primBuffer.data(t11off + 150 * idx + 42);

                auto t11_xxz_yzzz = primBuffer.data(t11off + 150 * idx + 43);

                auto t11_xxz_zzzz = primBuffer.data(t11off + 150 * idx + 44);

                auto t11_xyy_xxxx = primBuffer.data(t11off + 150 * idx + 45);

                auto t11_xyy_xxxy = primBuffer.data(t11off + 150 * idx + 46);

                auto t11_xyy_xxxz = primBuffer.data(t11off + 150 * idx + 47);

                auto t11_xyy_xxyy = primBuffer.data(t11off + 150 * idx + 48);

                auto t11_xyy_xxyz = primBuffer.data(t11off + 150 * idx + 49);

                auto t11_xyy_xxzz = primBuffer.data(t11off + 150 * idx + 50);

                auto t11_xyy_xyyy = primBuffer.data(t11off + 150 * idx + 51);

                auto t11_xyy_xyyz = primBuffer.data(t11off + 150 * idx + 52);

                auto t11_xyy_xyzz = primBuffer.data(t11off + 150 * idx + 53);

                auto t11_xyy_xzzz = primBuffer.data(t11off + 150 * idx + 54);

                auto t11_xyy_yyyy = primBuffer.data(t11off + 150 * idx + 55);

                auto t11_xyy_yyyz = primBuffer.data(t11off + 150 * idx + 56);

                auto t11_xyy_yyzz = primBuffer.data(t11off + 150 * idx + 57);

                auto t11_xyy_yzzz = primBuffer.data(t11off + 150 * idx + 58);

                auto t11_xyy_zzzz = primBuffer.data(t11off + 150 * idx + 59);

                auto t11_xyz_xxxx = primBuffer.data(t11off + 150 * idx + 60);

                auto t11_xyz_xxxy = primBuffer.data(t11off + 150 * idx + 61);

                auto t11_xyz_xxxz = primBuffer.data(t11off + 150 * idx + 62);

                auto t11_xyz_xxyy = primBuffer.data(t11off + 150 * idx + 63);

                auto t11_xyz_xxyz = primBuffer.data(t11off + 150 * idx + 64);

                auto t11_xyz_xxzz = primBuffer.data(t11off + 150 * idx + 65);

                auto t11_xyz_xyyy = primBuffer.data(t11off + 150 * idx + 66);

                auto t11_xyz_xyyz = primBuffer.data(t11off + 150 * idx + 67);

                auto t11_xyz_xyzz = primBuffer.data(t11off + 150 * idx + 68);

                auto t11_xyz_xzzz = primBuffer.data(t11off + 150 * idx + 69);

                auto t11_xyz_yyyy = primBuffer.data(t11off + 150 * idx + 70);

                auto t11_xyz_yyyz = primBuffer.data(t11off + 150 * idx + 71);

                auto t11_xyz_yyzz = primBuffer.data(t11off + 150 * idx + 72);

                auto t11_xyz_yzzz = primBuffer.data(t11off + 150 * idx + 73);

                auto t11_xyz_zzzz = primBuffer.data(t11off + 150 * idx + 74);

                auto t11_xzz_xxxx = primBuffer.data(t11off + 150 * idx + 75);

                auto t11_xzz_xxxy = primBuffer.data(t11off + 150 * idx + 76);

                auto t11_xzz_xxxz = primBuffer.data(t11off + 150 * idx + 77);

                auto t11_xzz_xxyy = primBuffer.data(t11off + 150 * idx + 78);

                auto t11_xzz_xxyz = primBuffer.data(t11off + 150 * idx + 79);

                auto t11_xzz_xxzz = primBuffer.data(t11off + 150 * idx + 80);

                auto t11_xzz_xyyy = primBuffer.data(t11off + 150 * idx + 81);

                auto t11_xzz_xyyz = primBuffer.data(t11off + 150 * idx + 82);

                auto t11_xzz_xyzz = primBuffer.data(t11off + 150 * idx + 83);

                auto t11_xzz_xzzz = primBuffer.data(t11off + 150 * idx + 84);

                auto t11_xzz_yyyy = primBuffer.data(t11off + 150 * idx + 85);

                auto t11_xzz_yyyz = primBuffer.data(t11off + 150 * idx + 86);

                auto t11_xzz_yyzz = primBuffer.data(t11off + 150 * idx + 87);

                auto t11_xzz_yzzz = primBuffer.data(t11off + 150 * idx + 88);

                auto t11_xzz_zzzz = primBuffer.data(t11off + 150 * idx + 89);

                auto t11_yyy_xxxx = primBuffer.data(t11off + 150 * idx + 90);

                auto t11_yyy_xxxy = primBuffer.data(t11off + 150 * idx + 91);

                auto t11_yyy_xxxz = primBuffer.data(t11off + 150 * idx + 92);

                auto t11_yyy_xxyy = primBuffer.data(t11off + 150 * idx + 93);

                auto t11_yyy_xxyz = primBuffer.data(t11off + 150 * idx + 94);

                auto t11_yyy_xxzz = primBuffer.data(t11off + 150 * idx + 95);

                auto t11_yyy_xyyy = primBuffer.data(t11off + 150 * idx + 96);

                auto t11_yyy_xyyz = primBuffer.data(t11off + 150 * idx + 97);

                auto t11_yyy_xyzz = primBuffer.data(t11off + 150 * idx + 98);

                auto t11_yyy_xzzz = primBuffer.data(t11off + 150 * idx + 99);

                auto t11_yyy_yyyy = primBuffer.data(t11off + 150 * idx + 100);

                auto t11_yyy_yyyz = primBuffer.data(t11off + 150 * idx + 101);

                auto t11_yyy_yyzz = primBuffer.data(t11off + 150 * idx + 102);

                auto t11_yyy_yzzz = primBuffer.data(t11off + 150 * idx + 103);

                auto t11_yyy_zzzz = primBuffer.data(t11off + 150 * idx + 104);

                auto t11_yyz_xxxx = primBuffer.data(t11off + 150 * idx + 105);

                auto t11_yyz_xxxy = primBuffer.data(t11off + 150 * idx + 106);

                auto t11_yyz_xxxz = primBuffer.data(t11off + 150 * idx + 107);

                auto t11_yyz_xxyy = primBuffer.data(t11off + 150 * idx + 108);

                auto t11_yyz_xxyz = primBuffer.data(t11off + 150 * idx + 109);

                auto t11_yyz_xxzz = primBuffer.data(t11off + 150 * idx + 110);

                auto t11_yyz_xyyy = primBuffer.data(t11off + 150 * idx + 111);

                auto t11_yyz_xyyz = primBuffer.data(t11off + 150 * idx + 112);

                auto t11_yyz_xyzz = primBuffer.data(t11off + 150 * idx + 113);

                auto t11_yyz_xzzz = primBuffer.data(t11off + 150 * idx + 114);

                auto t11_yyz_yyyy = primBuffer.data(t11off + 150 * idx + 115);

                auto t11_yyz_yyyz = primBuffer.data(t11off + 150 * idx + 116);

                auto t11_yyz_yyzz = primBuffer.data(t11off + 150 * idx + 117);

                auto t11_yyz_yzzz = primBuffer.data(t11off + 150 * idx + 118);

                auto t11_yyz_zzzz = primBuffer.data(t11off + 150 * idx + 119);

                auto t11_yzz_xxxx = primBuffer.data(t11off + 150 * idx + 120);

                auto t11_yzz_xxxy = primBuffer.data(t11off + 150 * idx + 121);

                auto t11_yzz_xxxz = primBuffer.data(t11off + 150 * idx + 122);

                auto t11_yzz_xxyy = primBuffer.data(t11off + 150 * idx + 123);

                auto t11_yzz_xxyz = primBuffer.data(t11off + 150 * idx + 124);

                auto t11_yzz_xxzz = primBuffer.data(t11off + 150 * idx + 125);

                auto t11_yzz_xyyy = primBuffer.data(t11off + 150 * idx + 126);

                auto t11_yzz_xyyz = primBuffer.data(t11off + 150 * idx + 127);

                auto t11_yzz_xyzz = primBuffer.data(t11off + 150 * idx + 128);

                auto t11_yzz_xzzz = primBuffer.data(t11off + 150 * idx + 129);

                auto t11_yzz_yyyy = primBuffer.data(t11off + 150 * idx + 130);

                auto t11_yzz_yyyz = primBuffer.data(t11off + 150 * idx + 131);

                auto t11_yzz_yyzz = primBuffer.data(t11off + 150 * idx + 132);

                auto t11_yzz_yzzz = primBuffer.data(t11off + 150 * idx + 133);

                auto t11_yzz_zzzz = primBuffer.data(t11off + 150 * idx + 134);

                auto t11_zzz_xxxx = primBuffer.data(t11off + 150 * idx + 135);

                auto t11_zzz_xxxy = primBuffer.data(t11off + 150 * idx + 136);

                auto t11_zzz_xxxz = primBuffer.data(t11off + 150 * idx + 137);

                auto t11_zzz_xxyy = primBuffer.data(t11off + 150 * idx + 138);

                auto t11_zzz_xxyz = primBuffer.data(t11off + 150 * idx + 139);

                auto t11_zzz_xxzz = primBuffer.data(t11off + 150 * idx + 140);

                auto t11_zzz_xyyy = primBuffer.data(t11off + 150 * idx + 141);

                auto t11_zzz_xyyz = primBuffer.data(t11off + 150 * idx + 142);

                auto t11_zzz_xyzz = primBuffer.data(t11off + 150 * idx + 143);

                auto t11_zzz_xzzz = primBuffer.data(t11off + 150 * idx + 144);

                auto t11_zzz_yyyy = primBuffer.data(t11off + 150 * idx + 145);

                auto t11_zzz_yyyz = primBuffer.data(t11off + 150 * idx + 146);

                auto t11_zzz_yyzz = primBuffer.data(t11off + 150 * idx + 147);

                auto t11_zzz_yzzz = primBuffer.data(t11off + 150 * idx + 148);

                auto t11_zzz_zzzz = primBuffer.data(t11off + 150 * idx + 149);

                // set up pointers to (G|A(0)|G)^(m) integrals

                auto t_xxxx_xxxx = primBuffer.data(toff + 225 * idx);

                auto t_xxxx_xxxy = primBuffer.data(toff + 225 * idx + 1);

                auto t_xxxx_xxxz = primBuffer.data(toff + 225 * idx + 2);

                auto t_xxxx_xxyy = primBuffer.data(toff + 225 * idx + 3);

                auto t_xxxx_xxyz = primBuffer.data(toff + 225 * idx + 4);

                auto t_xxxx_xxzz = primBuffer.data(toff + 225 * idx + 5);

                auto t_xxxx_xyyy = primBuffer.data(toff + 225 * idx + 6);

                auto t_xxxx_xyyz = primBuffer.data(toff + 225 * idx + 7);

                auto t_xxxx_xyzz = primBuffer.data(toff + 225 * idx + 8);

                auto t_xxxx_xzzz = primBuffer.data(toff + 225 * idx + 9);

                auto t_xxxx_yyyy = primBuffer.data(toff + 225 * idx + 10);

                auto t_xxxx_yyyz = primBuffer.data(toff + 225 * idx + 11);

                auto t_xxxx_yyzz = primBuffer.data(toff + 225 * idx + 12);

                auto t_xxxx_yzzz = primBuffer.data(toff + 225 * idx + 13);

                auto t_xxxx_zzzz = primBuffer.data(toff + 225 * idx + 14);

                auto t_xxxy_xxxx = primBuffer.data(toff + 225 * idx + 15);

                auto t_xxxy_xxxy = primBuffer.data(toff + 225 * idx + 16);

                auto t_xxxy_xxxz = primBuffer.data(toff + 225 * idx + 17);

                auto t_xxxy_xxyy = primBuffer.data(toff + 225 * idx + 18);

                auto t_xxxy_xxyz = primBuffer.data(toff + 225 * idx + 19);

                auto t_xxxy_xxzz = primBuffer.data(toff + 225 * idx + 20);

                auto t_xxxy_xyyy = primBuffer.data(toff + 225 * idx + 21);

                auto t_xxxy_xyyz = primBuffer.data(toff + 225 * idx + 22);

                auto t_xxxy_xyzz = primBuffer.data(toff + 225 * idx + 23);

                auto t_xxxy_xzzz = primBuffer.data(toff + 225 * idx + 24);

                auto t_xxxy_yyyy = primBuffer.data(toff + 225 * idx + 25);

                auto t_xxxy_yyyz = primBuffer.data(toff + 225 * idx + 26);

                auto t_xxxy_yyzz = primBuffer.data(toff + 225 * idx + 27);

                auto t_xxxy_yzzz = primBuffer.data(toff + 225 * idx + 28);

                auto t_xxxy_zzzz = primBuffer.data(toff + 225 * idx + 29);

                auto t_xxxz_xxxx = primBuffer.data(toff + 225 * idx + 30);

                auto t_xxxz_xxxy = primBuffer.data(toff + 225 * idx + 31);

                auto t_xxxz_xxxz = primBuffer.data(toff + 225 * idx + 32);

                auto t_xxxz_xxyy = primBuffer.data(toff + 225 * idx + 33);

                auto t_xxxz_xxyz = primBuffer.data(toff + 225 * idx + 34);

                auto t_xxxz_xxzz = primBuffer.data(toff + 225 * idx + 35);

                auto t_xxxz_xyyy = primBuffer.data(toff + 225 * idx + 36);

                auto t_xxxz_xyyz = primBuffer.data(toff + 225 * idx + 37);

                auto t_xxxz_xyzz = primBuffer.data(toff + 225 * idx + 38);

                auto t_xxxz_xzzz = primBuffer.data(toff + 225 * idx + 39);

                auto t_xxxz_yyyy = primBuffer.data(toff + 225 * idx + 40);

                auto t_xxxz_yyyz = primBuffer.data(toff + 225 * idx + 41);

                auto t_xxxz_yyzz = primBuffer.data(toff + 225 * idx + 42);

                auto t_xxxz_yzzz = primBuffer.data(toff + 225 * idx + 43);

                auto t_xxxz_zzzz = primBuffer.data(toff + 225 * idx + 44);

                auto t_xxyy_xxxx = primBuffer.data(toff + 225 * idx + 45);

                auto t_xxyy_xxxy = primBuffer.data(toff + 225 * idx + 46);

                auto t_xxyy_xxxz = primBuffer.data(toff + 225 * idx + 47);

                auto t_xxyy_xxyy = primBuffer.data(toff + 225 * idx + 48);

                auto t_xxyy_xxyz = primBuffer.data(toff + 225 * idx + 49);

                auto t_xxyy_xxzz = primBuffer.data(toff + 225 * idx + 50);

                auto t_xxyy_xyyy = primBuffer.data(toff + 225 * idx + 51);

                auto t_xxyy_xyyz = primBuffer.data(toff + 225 * idx + 52);

                auto t_xxyy_xyzz = primBuffer.data(toff + 225 * idx + 53);

                auto t_xxyy_xzzz = primBuffer.data(toff + 225 * idx + 54);

                auto t_xxyy_yyyy = primBuffer.data(toff + 225 * idx + 55);

                auto t_xxyy_yyyz = primBuffer.data(toff + 225 * idx + 56);

                auto t_xxyy_yyzz = primBuffer.data(toff + 225 * idx + 57);

                auto t_xxyy_yzzz = primBuffer.data(toff + 225 * idx + 58);

                auto t_xxyy_zzzz = primBuffer.data(toff + 225 * idx + 59);

                auto t_xxyz_xxxx = primBuffer.data(toff + 225 * idx + 60);

                auto t_xxyz_xxxy = primBuffer.data(toff + 225 * idx + 61);

                auto t_xxyz_xxxz = primBuffer.data(toff + 225 * idx + 62);

                auto t_xxyz_xxyy = primBuffer.data(toff + 225 * idx + 63);

                auto t_xxyz_xxyz = primBuffer.data(toff + 225 * idx + 64);

                auto t_xxyz_xxzz = primBuffer.data(toff + 225 * idx + 65);

                auto t_xxyz_xyyy = primBuffer.data(toff + 225 * idx + 66);

                auto t_xxyz_xyyz = primBuffer.data(toff + 225 * idx + 67);

                auto t_xxyz_xyzz = primBuffer.data(toff + 225 * idx + 68);

                auto t_xxyz_xzzz = primBuffer.data(toff + 225 * idx + 69);

                auto t_xxyz_yyyy = primBuffer.data(toff + 225 * idx + 70);

                auto t_xxyz_yyyz = primBuffer.data(toff + 225 * idx + 71);

                auto t_xxyz_yyzz = primBuffer.data(toff + 225 * idx + 72);

                auto t_xxyz_yzzz = primBuffer.data(toff + 225 * idx + 73);

                auto t_xxyz_zzzz = primBuffer.data(toff + 225 * idx + 74);

                auto t_xxzz_xxxx = primBuffer.data(toff + 225 * idx + 75);

                auto t_xxzz_xxxy = primBuffer.data(toff + 225 * idx + 76);

                auto t_xxzz_xxxz = primBuffer.data(toff + 225 * idx + 77);

                auto t_xxzz_xxyy = primBuffer.data(toff + 225 * idx + 78);

                auto t_xxzz_xxyz = primBuffer.data(toff + 225 * idx + 79);

                auto t_xxzz_xxzz = primBuffer.data(toff + 225 * idx + 80);

                auto t_xxzz_xyyy = primBuffer.data(toff + 225 * idx + 81);

                auto t_xxzz_xyyz = primBuffer.data(toff + 225 * idx + 82);

                auto t_xxzz_xyzz = primBuffer.data(toff + 225 * idx + 83);

                auto t_xxzz_xzzz = primBuffer.data(toff + 225 * idx + 84);

                auto t_xxzz_yyyy = primBuffer.data(toff + 225 * idx + 85);

                auto t_xxzz_yyyz = primBuffer.data(toff + 225 * idx + 86);

                auto t_xxzz_yyzz = primBuffer.data(toff + 225 * idx + 87);

                auto t_xxzz_yzzz = primBuffer.data(toff + 225 * idx + 88);

                auto t_xxzz_zzzz = primBuffer.data(toff + 225 * idx + 89);

                auto t_xyyy_xxxx = primBuffer.data(toff + 225 * idx + 90);

                auto t_xyyy_xxxy = primBuffer.data(toff + 225 * idx + 91);

                auto t_xyyy_xxxz = primBuffer.data(toff + 225 * idx + 92);

                auto t_xyyy_xxyy = primBuffer.data(toff + 225 * idx + 93);

                auto t_xyyy_xxyz = primBuffer.data(toff + 225 * idx + 94);

                auto t_xyyy_xxzz = primBuffer.data(toff + 225 * idx + 95);

                auto t_xyyy_xyyy = primBuffer.data(toff + 225 * idx + 96);

                auto t_xyyy_xyyz = primBuffer.data(toff + 225 * idx + 97);

                auto t_xyyy_xyzz = primBuffer.data(toff + 225 * idx + 98);

                auto t_xyyy_xzzz = primBuffer.data(toff + 225 * idx + 99);

                auto t_xyyy_yyyy = primBuffer.data(toff + 225 * idx + 100);

                auto t_xyyy_yyyz = primBuffer.data(toff + 225 * idx + 101);

                auto t_xyyy_yyzz = primBuffer.data(toff + 225 * idx + 102);

                auto t_xyyy_yzzz = primBuffer.data(toff + 225 * idx + 103);

                auto t_xyyy_zzzz = primBuffer.data(toff + 225 * idx + 104);

                auto t_xyyz_xxxx = primBuffer.data(toff + 225 * idx + 105);

                auto t_xyyz_xxxy = primBuffer.data(toff + 225 * idx + 106);

                auto t_xyyz_xxxz = primBuffer.data(toff + 225 * idx + 107);

                auto t_xyyz_xxyy = primBuffer.data(toff + 225 * idx + 108);

                auto t_xyyz_xxyz = primBuffer.data(toff + 225 * idx + 109);

                auto t_xyyz_xxzz = primBuffer.data(toff + 225 * idx + 110);

                auto t_xyyz_xyyy = primBuffer.data(toff + 225 * idx + 111);

                auto t_xyyz_xyyz = primBuffer.data(toff + 225 * idx + 112);

                auto t_xyyz_xyzz = primBuffer.data(toff + 225 * idx + 113);

                auto t_xyyz_xzzz = primBuffer.data(toff + 225 * idx + 114);

                auto t_xyyz_yyyy = primBuffer.data(toff + 225 * idx + 115);

                auto t_xyyz_yyyz = primBuffer.data(toff + 225 * idx + 116);

                auto t_xyyz_yyzz = primBuffer.data(toff + 225 * idx + 117);

                auto t_xyyz_yzzz = primBuffer.data(toff + 225 * idx + 118);

                auto t_xyyz_zzzz = primBuffer.data(toff + 225 * idx + 119);

                auto t_xyzz_xxxx = primBuffer.data(toff + 225 * idx + 120);

                auto t_xyzz_xxxy = primBuffer.data(toff + 225 * idx + 121);

                auto t_xyzz_xxxz = primBuffer.data(toff + 225 * idx + 122);

                auto t_xyzz_xxyy = primBuffer.data(toff + 225 * idx + 123);

                auto t_xyzz_xxyz = primBuffer.data(toff + 225 * idx + 124);

                auto t_xyzz_xxzz = primBuffer.data(toff + 225 * idx + 125);

                auto t_xyzz_xyyy = primBuffer.data(toff + 225 * idx + 126);

                auto t_xyzz_xyyz = primBuffer.data(toff + 225 * idx + 127);

                auto t_xyzz_xyzz = primBuffer.data(toff + 225 * idx + 128);

                auto t_xyzz_xzzz = primBuffer.data(toff + 225 * idx + 129);

                auto t_xyzz_yyyy = primBuffer.data(toff + 225 * idx + 130);

                auto t_xyzz_yyyz = primBuffer.data(toff + 225 * idx + 131);

                auto t_xyzz_yyzz = primBuffer.data(toff + 225 * idx + 132);

                auto t_xyzz_yzzz = primBuffer.data(toff + 225 * idx + 133);

                auto t_xyzz_zzzz = primBuffer.data(toff + 225 * idx + 134);

                auto t_xzzz_xxxx = primBuffer.data(toff + 225 * idx + 135);

                auto t_xzzz_xxxy = primBuffer.data(toff + 225 * idx + 136);

                auto t_xzzz_xxxz = primBuffer.data(toff + 225 * idx + 137);

                auto t_xzzz_xxyy = primBuffer.data(toff + 225 * idx + 138);

                auto t_xzzz_xxyz = primBuffer.data(toff + 225 * idx + 139);

                auto t_xzzz_xxzz = primBuffer.data(toff + 225 * idx + 140);

                auto t_xzzz_xyyy = primBuffer.data(toff + 225 * idx + 141);

                auto t_xzzz_xyyz = primBuffer.data(toff + 225 * idx + 142);

                auto t_xzzz_xyzz = primBuffer.data(toff + 225 * idx + 143);

                auto t_xzzz_xzzz = primBuffer.data(toff + 225 * idx + 144);

                auto t_xzzz_yyyy = primBuffer.data(toff + 225 * idx + 145);

                auto t_xzzz_yyyz = primBuffer.data(toff + 225 * idx + 146);

                auto t_xzzz_yyzz = primBuffer.data(toff + 225 * idx + 147);

                auto t_xzzz_yzzz = primBuffer.data(toff + 225 * idx + 148);

                auto t_xzzz_zzzz = primBuffer.data(toff + 225 * idx + 149);

                auto t_yyyy_xxxx = primBuffer.data(toff + 225 * idx + 150);

                auto t_yyyy_xxxy = primBuffer.data(toff + 225 * idx + 151);

                auto t_yyyy_xxxz = primBuffer.data(toff + 225 * idx + 152);

                auto t_yyyy_xxyy = primBuffer.data(toff + 225 * idx + 153);

                auto t_yyyy_xxyz = primBuffer.data(toff + 225 * idx + 154);

                auto t_yyyy_xxzz = primBuffer.data(toff + 225 * idx + 155);

                auto t_yyyy_xyyy = primBuffer.data(toff + 225 * idx + 156);

                auto t_yyyy_xyyz = primBuffer.data(toff + 225 * idx + 157);

                auto t_yyyy_xyzz = primBuffer.data(toff + 225 * idx + 158);

                auto t_yyyy_xzzz = primBuffer.data(toff + 225 * idx + 159);

                auto t_yyyy_yyyy = primBuffer.data(toff + 225 * idx + 160);

                auto t_yyyy_yyyz = primBuffer.data(toff + 225 * idx + 161);

                auto t_yyyy_yyzz = primBuffer.data(toff + 225 * idx + 162);

                auto t_yyyy_yzzz = primBuffer.data(toff + 225 * idx + 163);

                auto t_yyyy_zzzz = primBuffer.data(toff + 225 * idx + 164);

                auto t_yyyz_xxxx = primBuffer.data(toff + 225 * idx + 165);

                auto t_yyyz_xxxy = primBuffer.data(toff + 225 * idx + 166);

                auto t_yyyz_xxxz = primBuffer.data(toff + 225 * idx + 167);

                auto t_yyyz_xxyy = primBuffer.data(toff + 225 * idx + 168);

                auto t_yyyz_xxyz = primBuffer.data(toff + 225 * idx + 169);

                auto t_yyyz_xxzz = primBuffer.data(toff + 225 * idx + 170);

                auto t_yyyz_xyyy = primBuffer.data(toff + 225 * idx + 171);

                auto t_yyyz_xyyz = primBuffer.data(toff + 225 * idx + 172);

                auto t_yyyz_xyzz = primBuffer.data(toff + 225 * idx + 173);

                auto t_yyyz_xzzz = primBuffer.data(toff + 225 * idx + 174);

                auto t_yyyz_yyyy = primBuffer.data(toff + 225 * idx + 175);

                auto t_yyyz_yyyz = primBuffer.data(toff + 225 * idx + 176);

                auto t_yyyz_yyzz = primBuffer.data(toff + 225 * idx + 177);

                auto t_yyyz_yzzz = primBuffer.data(toff + 225 * idx + 178);

                auto t_yyyz_zzzz = primBuffer.data(toff + 225 * idx + 179);

                auto t_yyzz_xxxx = primBuffer.data(toff + 225 * idx + 180);

                auto t_yyzz_xxxy = primBuffer.data(toff + 225 * idx + 181);

                auto t_yyzz_xxxz = primBuffer.data(toff + 225 * idx + 182);

                auto t_yyzz_xxyy = primBuffer.data(toff + 225 * idx + 183);

                auto t_yyzz_xxyz = primBuffer.data(toff + 225 * idx + 184);

                auto t_yyzz_xxzz = primBuffer.data(toff + 225 * idx + 185);

                auto t_yyzz_xyyy = primBuffer.data(toff + 225 * idx + 186);

                auto t_yyzz_xyyz = primBuffer.data(toff + 225 * idx + 187);

                auto t_yyzz_xyzz = primBuffer.data(toff + 225 * idx + 188);

                auto t_yyzz_xzzz = primBuffer.data(toff + 225 * idx + 189);

                auto t_yyzz_yyyy = primBuffer.data(toff + 225 * idx + 190);

                auto t_yyzz_yyyz = primBuffer.data(toff + 225 * idx + 191);

                auto t_yyzz_yyzz = primBuffer.data(toff + 225 * idx + 192);

                auto t_yyzz_yzzz = primBuffer.data(toff + 225 * idx + 193);

                auto t_yyzz_zzzz = primBuffer.data(toff + 225 * idx + 194);

                auto t_yzzz_xxxx = primBuffer.data(toff + 225 * idx + 195);

                auto t_yzzz_xxxy = primBuffer.data(toff + 225 * idx + 196);

                auto t_yzzz_xxxz = primBuffer.data(toff + 225 * idx + 197);

                auto t_yzzz_xxyy = primBuffer.data(toff + 225 * idx + 198);

                auto t_yzzz_xxyz = primBuffer.data(toff + 225 * idx + 199);

                auto t_yzzz_xxzz = primBuffer.data(toff + 225 * idx + 200);

                auto t_yzzz_xyyy = primBuffer.data(toff + 225 * idx + 201);

                auto t_yzzz_xyyz = primBuffer.data(toff + 225 * idx + 202);

                auto t_yzzz_xyzz = primBuffer.data(toff + 225 * idx + 203);

                auto t_yzzz_xzzz = primBuffer.data(toff + 225 * idx + 204);

                auto t_yzzz_yyyy = primBuffer.data(toff + 225 * idx + 205);

                auto t_yzzz_yyyz = primBuffer.data(toff + 225 * idx + 206);

                auto t_yzzz_yyzz = primBuffer.data(toff + 225 * idx + 207);

                auto t_yzzz_yzzz = primBuffer.data(toff + 225 * idx + 208);

                auto t_yzzz_zzzz = primBuffer.data(toff + 225 * idx + 209);

                auto t_zzzz_xxxx = primBuffer.data(toff + 225 * idx + 210);

                auto t_zzzz_xxxy = primBuffer.data(toff + 225 * idx + 211);

                auto t_zzzz_xxxz = primBuffer.data(toff + 225 * idx + 212);

                auto t_zzzz_xxyy = primBuffer.data(toff + 225 * idx + 213);

                auto t_zzzz_xxyz = primBuffer.data(toff + 225 * idx + 214);

                auto t_zzzz_xxzz = primBuffer.data(toff + 225 * idx + 215);

                auto t_zzzz_xyyy = primBuffer.data(toff + 225 * idx + 216);

                auto t_zzzz_xyyz = primBuffer.data(toff + 225 * idx + 217);

                auto t_zzzz_xyzz = primBuffer.data(toff + 225 * idx + 218);

                auto t_zzzz_xzzz = primBuffer.data(toff + 225 * idx + 219);

                auto t_zzzz_yyyy = primBuffer.data(toff + 225 * idx + 220);

                auto t_zzzz_yyyz = primBuffer.data(toff + 225 * idx + 221);

                auto t_zzzz_yyzz = primBuffer.data(toff + 225 * idx + 222);

                auto t_zzzz_yzzz = primBuffer.data(toff + 225 * idx + 223);

                auto t_zzzz_zzzz = primBuffer.data(toff + 225 * idx + 224);

                #pragma omp simd aligned(fx, pax, pay, paz, tk0_xxx_xxx, tk0_xxx_xxy,\
                                         tk0_xxx_xxz, tk0_xxx_xyy, tk0_xxx_xyz,\
                                         tk0_xxx_xzz, tk0_xxx_yyy, tk0_xxx_yyz,\
                                         tk0_xxx_yzz, tk0_xxx_zzz, tk0_xxy_xxx,\
                                         tk0_xxy_xxy, tk0_xxy_xxz, tk0_xxy_xyy,\
                                         tk0_xxy_xyz, tk0_xxy_xzz, tk0_xxy_yyy,\
                                         tk0_xxy_yyz, tk0_xxy_yzz, tk0_xxy_zzz,\
                                         tk0_xxz_xxx, tk0_xxz_xxy, tk0_xxz_xxz,\
                                         tk0_xxz_xyy, tk0_xxz_xyz, tk0_xxz_xzz,\
                                         tk0_xxz_yyy, tk0_xxz_yyz, tk0_xxz_yzz,\
                                         tk0_xxz_zzz, tk0_xyy_xxx, tk0_xyy_xxy,\
                                         tk0_xyy_xxz, tk0_xyy_xyy, tk0_xyy_xyz,\
                                         tk0_xyy_xzz, tk0_xyy_yyy, tk0_xyy_yyz,\
                                         tk0_xyy_yzz, tk0_xyy_zzz, tk0_xyz_xxx,\
                                         tk0_xyz_xxy, tk0_xyz_xxz, tk0_xyz_xyy,\
                                         tk0_xyz_xyz, tk0_xyz_xzz, tk0_xyz_yyy,\
                                         tk0_xyz_yyz, tk0_xyz_yzz, tk0_xyz_zzz,\
                                         tk0_xzz_xxx, tk0_xzz_xxy, tk0_xzz_xxz,\
                                         tk0_xzz_xyy, tk0_xzz_xyz, tk0_xzz_xzz,\
                                         tk0_xzz_yyy, tk0_xzz_yyz, tk0_xzz_yzz,\
                                         tk0_xzz_zzz, tk0_yyy_xxx, tk0_yyy_xxy,\
                                         tk0_yyy_xxz, tk0_yyy_xyy, tk0_yyy_xyz,\
                                         tk0_yyy_xzz, tk0_yyy_yyy, tk0_yyy_yyz,\
                                         tk0_yyy_yzz, tk0_yyy_zzz, tk0_yyz_xxx,\
                                         tk0_yyz_xxy, tk0_yyz_xxz, tk0_yyz_xyy,\
                                         tk0_yyz_xyz, tk0_yyz_xzz, tk0_yyz_yyy,\
                                         tk0_yyz_yyz, tk0_yyz_yzz, tk0_yyz_zzz,\
                                         tk0_yzz_xxx, tk0_yzz_xxy, tk0_yzz_xxz,\
                                         tk0_yzz_xyy, tk0_yzz_xyz, tk0_yzz_xzz,\
                                         tk0_yzz_yyy, tk0_yzz_yyz, tk0_yzz_yzz,\
                                         tk0_yzz_zzz, tk0_zzz_xxx, tk0_zzz_xxy,\
                                         tk0_zzz_xxz, tk0_zzz_xyy, tk0_zzz_xyz,\
                                         tk0_zzz_xzz, tk0_zzz_yyy, tk0_zzz_yyz,\
                                         tk0_zzz_yzz, tk0_zzz_zzz, tk1_xxx_xxx,\
                                         tk1_xxx_xxy, tk1_xxx_xxz, tk1_xxx_xyy,\
                                         tk1_xxx_xyz, tk1_xxx_xzz, tk1_xxx_yyy,\
                                         tk1_xxx_yyz, tk1_xxx_yzz, tk1_xxx_zzz,\
                                         tk1_xxy_xxx, tk1_xxy_xxy, tk1_xxy_xxz,\
                                         tk1_xxy_xyy, tk1_xxy_xyz, tk1_xxy_xzz,\
                                         tk1_xxy_yyy, tk1_xxy_yyz, tk1_xxy_yzz,\
                                         tk1_xxy_zzz, tk1_xxz_xxx, tk1_xxz_xxy,\
                                         tk1_xxz_xxz, tk1_xxz_xyy, tk1_xxz_xyz,\
                                         tk1_xxz_xzz, tk1_xxz_yyy, tk1_xxz_yyz,\
                                         tk1_xxz_yzz, tk1_xxz_zzz, tk1_xyy_xxx,\
                                         tk1_xyy_xxy, tk1_xyy_xxz, tk1_xyy_xyy,\
                                         tk1_xyy_xyz, tk1_xyy_xzz, tk1_xyy_yyy,\
                                         tk1_xyy_yyz, tk1_xyy_yzz, tk1_xyy_zzz,\
                                         tk1_xyz_xxx, tk1_xyz_xxy, tk1_xyz_xxz,\
                                         tk1_xyz_xyy, tk1_xyz_xyz, tk1_xyz_xzz,\
                                         tk1_xyz_yyy, tk1_xyz_yyz, tk1_xyz_yzz,\
                                         tk1_xyz_zzz, tk1_xzz_xxx, tk1_xzz_xxy,\
                                         tk1_xzz_xxz, tk1_xzz_xyy, tk1_xzz_xyz,\
                                         tk1_xzz_xzz, tk1_xzz_yyy, tk1_xzz_yyz,\
                                         tk1_xzz_yzz, tk1_xzz_zzz, tk1_yyy_xxx,\
                                         tk1_yyy_xxy, tk1_yyy_xxz, tk1_yyy_xyy,\
                                         tk1_yyy_xyz, tk1_yyy_xzz, tk1_yyy_yyy,\
                                         tk1_yyy_yyz, tk1_yyy_yzz, tk1_yyy_zzz,\
                                         tk1_yyz_xxx, tk1_yyz_xxy, tk1_yyz_xxz,\
                                         tk1_yyz_xyy, tk1_yyz_xyz, tk1_yyz_xzz,\
                                         tk1_yyz_yyy, tk1_yyz_yyz, tk1_yyz_yzz,\
                                         tk1_yyz_zzz, tk1_yzz_xxx, tk1_yzz_xxy,\
                                         tk1_yzz_xxz, tk1_yzz_xyy, tk1_yzz_xyz,\
                                         tk1_yzz_xzz, tk1_yzz_yyy, tk1_yzz_yyz,\
                                         tk1_yzz_yzz, tk1_yzz_zzz, tk1_zzz_xxx,\
                                         tk1_zzz_xxy, tk1_zzz_xxz, tk1_zzz_xyy,\
                                         tk1_zzz_xyz, tk1_zzz_xzz, tk1_zzz_yyy,\
                                         tk1_zzz_yyz, tk1_zzz_yzz, tk1_zzz_zzz,\
                                         t20_xx_xxxx, t20_xx_xxxy, t20_xx_xxxz,\
                                         t20_xx_xxyy, t20_xx_xxyz, t20_xx_xxzz,\
                                         t20_xx_xyyy, t20_xx_xyyz, t20_xx_xyzz,\
                                         t20_xx_xzzz, t20_xx_yyyy, t20_xx_yyyz,\
                                         t20_xx_yyzz, t20_xx_yzzz, t20_xx_zzzz,\
                                         t20_xy_xxxx, t20_xy_xxxy, t20_xy_xxxz,\
                                         t20_xy_xxyy, t20_xy_xxyz, t20_xy_xxzz,\
                                         t20_xy_xyyy, t20_xy_xyyz, t20_xy_xyzz,\
                                         t20_xy_xzzz, t20_xy_yyyy, t20_xy_yyyz,\
                                         t20_xy_yyzz, t20_xy_yzzz, t20_xy_zzzz,\
                                         t20_xz_xxxx, t20_xz_xxxy, t20_xz_xxxz,\
                                         t20_xz_xxyy, t20_xz_xxyz, t20_xz_xxzz,\
                                         t20_xz_xyyy, t20_xz_xyyz, t20_xz_xyzz,\
                                         t20_xz_xzzz, t20_xz_yyyy, t20_xz_yyyz,\
                                         t20_xz_yyzz, t20_xz_yzzz, t20_xz_zzzz,\
                                         t20_yy_xxxx, t20_yy_xxxy, t20_yy_xxxz,\
                                         t20_yy_xxyy, t20_yy_xxyz, t20_yy_xxzz,\
                                         t20_yy_xyyy, t20_yy_xyyz, t20_yy_xyzz,\
                                         t20_yy_xzzz, t20_yy_yyyy, t20_yy_yyyz,\
                                         t20_yy_yyzz, t20_yy_yzzz, t20_yy_zzzz,\
                                         t20_yz_xxxx, t20_yz_xxxy, t20_yz_xxxz,\
                                         t20_yz_xxyy, t20_yz_xxyz, t20_yz_xxzz,\
                                         t20_yz_xyyy, t20_yz_xyyz, t20_yz_xyzz,\
                                         t20_yz_xzzz, t20_yz_yyyy, t20_yz_yyyz,\
                                         t20_yz_yyzz, t20_yz_yzzz, t20_yz_zzzz,\
                                         t20_zz_xxxx, t20_zz_xxxy, t20_zz_xxxz,\
                                         t20_zz_xxyy, t20_zz_xxyz, t20_zz_xxzz,\
                                         t20_zz_xyyy, t20_zz_xyyz, t20_zz_xyzz,\
                                         t20_zz_xzzz, t20_zz_yyyy, t20_zz_yyyz,\
                                         t20_zz_yyzz, t20_zz_yzzz, t20_zz_zzzz,\
                                         t21_xx_xxxx, t21_xx_xxxy, t21_xx_xxxz,\
                                         t21_xx_xxyy, t21_xx_xxyz, t21_xx_xxzz,\
                                         t21_xx_xyyy, t21_xx_xyyz, t21_xx_xyzz,\
                                         t21_xx_xzzz, t21_xx_yyyy, t21_xx_yyyz,\
                                         t21_xx_yyzz, t21_xx_yzzz, t21_xx_zzzz,\
                                         t21_xy_xxxx, t21_xy_xxxy, t21_xy_xxxz,\
                                         t21_xy_xxyy, t21_xy_xxyz, t21_xy_xxzz,\
                                         t21_xy_xyyy, t21_xy_xyyz, t21_xy_xyzz,\
                                         t21_xy_xzzz, t21_xy_yyyy, t21_xy_yyyz,\
                                         t21_xy_yyzz, t21_xy_yzzz, t21_xy_zzzz,\
                                         t21_xz_xxxx, t21_xz_xxxy, t21_xz_xxxz,\
                                         t21_xz_xxyy, t21_xz_xxyz, t21_xz_xxzz,\
                                         t21_xz_xyyy, t21_xz_xyyz, t21_xz_xyzz,\
                                         t21_xz_xzzz, t21_xz_yyyy, t21_xz_yyyz,\
                                         t21_xz_yyzz, t21_xz_yzzz, t21_xz_zzzz,\
                                         t21_yy_xxxx, t21_yy_xxxy, t21_yy_xxxz,\
                                         t21_yy_xxyy, t21_yy_xxyz, t21_yy_xxzz,\
                                         t21_yy_xyyy, t21_yy_xyyz, t21_yy_xyzz,\
                                         t21_yy_xzzz, t21_yy_yyyy, t21_yy_yyyz,\
                                         t21_yy_yyzz, t21_yy_yzzz, t21_yy_zzzz,\
                                         t21_yz_xxxx, t21_yz_xxxy, t21_yz_xxxz,\
                                         t21_yz_xxyy, t21_yz_xxyz, t21_yz_xxzz,\
                                         t21_yz_xyyy, t21_yz_xyyz, t21_yz_xyzz,\
                                         t21_yz_xzzz, t21_yz_yyyy, t21_yz_yyyz,\
                                         t21_yz_yyzz, t21_yz_yzzz, t21_yz_zzzz,\
                                         t21_zz_xxxx, t21_zz_xxxy, t21_zz_xxxz,\
                                         t21_zz_xxyy, t21_zz_xxyz, t21_zz_xxzz,\
                                         t21_zz_xyyy, t21_zz_xyyz, t21_zz_xyzz,\
                                         t21_zz_xzzz, t21_zz_yyyy, t21_zz_yyyz,\
                                         t21_zz_yyzz, t21_zz_yzzz, t21_zz_zzzz,\
                                         t10_xxx_xxxx, t10_xxx_xxxy, t10_xxx_xxxz,\
                                         t10_xxx_xxyy, t10_xxx_xxyz, t10_xxx_xxzz,\
                                         t10_xxx_xyyy, t10_xxx_xyyz, t10_xxx_xyzz,\
                                         t10_xxx_xzzz, t10_xxx_yyyy, t10_xxx_yyyz,\
                                         t10_xxx_yyzz, t10_xxx_yzzz, t10_xxx_zzzz,\
                                         t10_xxy_xxxx, t10_xxy_xxxy, t10_xxy_xxxz,\
                                         t10_xxy_xxyy, t10_xxy_xxyz, t10_xxy_xxzz,\
                                         t10_xxy_xyyy, t10_xxy_xyyz, t10_xxy_xyzz,\
                                         t10_xxy_xzzz, t10_xxy_yyyy, t10_xxy_yyyz,\
                                         t10_xxy_yyzz, t10_xxy_yzzz, t10_xxy_zzzz,\
                                         t10_xxz_xxxx, t10_xxz_xxxy, t10_xxz_xxxz,\
                                         t10_xxz_xxyy, t10_xxz_xxyz, t10_xxz_xxzz,\
                                         t10_xxz_xyyy, t10_xxz_xyyz, t10_xxz_xyzz,\
                                         t10_xxz_xzzz, t10_xxz_yyyy, t10_xxz_yyyz,\
                                         t10_xxz_yyzz, t10_xxz_yzzz, t10_xxz_zzzz,\
                                         t10_xyy_xxxx, t10_xyy_xxxy, t10_xyy_xxxz,\
                                         t10_xyy_xxyy, t10_xyy_xxyz, t10_xyy_xxzz,\
                                         t10_xyy_xyyy, t10_xyy_xyyz, t10_xyy_xyzz,\
                                         t10_xyy_xzzz, t10_xyy_yyyy, t10_xyy_yyyz,\
                                         t10_xyy_yyzz, t10_xyy_yzzz, t10_xyy_zzzz,\
                                         t10_xyz_xxxx, t10_xyz_xxxy, t10_xyz_xxxz,\
                                         t10_xyz_xxyy, t10_xyz_xxyz, t10_xyz_xxzz,\
                                         t10_xyz_xyyy, t10_xyz_xyyz, t10_xyz_xyzz,\
                                         t10_xyz_xzzz, t10_xyz_yyyy, t10_xyz_yyyz,\
                                         t10_xyz_yyzz, t10_xyz_yzzz, t10_xyz_zzzz,\
                                         t10_xzz_xxxx, t10_xzz_xxxy, t10_xzz_xxxz,\
                                         t10_xzz_xxyy, t10_xzz_xxyz, t10_xzz_xxzz,\
                                         t10_xzz_xyyy, t10_xzz_xyyz, t10_xzz_xyzz,\
                                         t10_xzz_xzzz, t10_xzz_yyyy, t10_xzz_yyyz,\
                                         t10_xzz_yyzz, t10_xzz_yzzz, t10_xzz_zzzz,\
                                         t10_yyy_xxxx, t10_yyy_xxxy, t10_yyy_xxxz,\
                                         t10_yyy_xxyy, t10_yyy_xxyz, t10_yyy_xxzz,\
                                         t10_yyy_xyyy, t10_yyy_xyyz, t10_yyy_xyzz,\
                                         t10_yyy_xzzz, t10_yyy_yyyy, t10_yyy_yyyz,\
                                         t10_yyy_yyzz, t10_yyy_yzzz, t10_yyy_zzzz,\
                                         t10_yyz_xxxx, t10_yyz_xxxy, t10_yyz_xxxz,\
                                         t10_yyz_xxyy, t10_yyz_xxyz, t10_yyz_xxzz,\
                                         t10_yyz_xyyy, t10_yyz_xyyz, t10_yyz_xyzz,\
                                         t10_yyz_xzzz, t10_yyz_yyyy, t10_yyz_yyyz,\
                                         t10_yyz_yyzz, t10_yyz_yzzz, t10_yyz_zzzz,\
                                         t10_yzz_xxxx, t10_yzz_xxxy, t10_yzz_xxxz,\
                                         t10_yzz_xxyy, t10_yzz_xxyz, t10_yzz_xxzz,\
                                         t10_yzz_xyyy, t10_yzz_xyyz, t10_yzz_xyzz,\
                                         t10_yzz_xzzz, t10_yzz_yyyy, t10_yzz_yyyz,\
                                         t10_yzz_yyzz, t10_yzz_yzzz, t10_yzz_zzzz,\
                                         t10_zzz_xxxx, t10_zzz_xxxy, t10_zzz_xxxz,\
                                         t10_zzz_xxyy, t10_zzz_xxyz, t10_zzz_xxzz,\
                                         t10_zzz_xyyy, t10_zzz_xyyz, t10_zzz_xyzz,\
                                         t10_zzz_xzzz, t10_zzz_yyyy, t10_zzz_yyyz,\
                                         t10_zzz_yyzz, t10_zzz_yzzz, t10_zzz_zzzz,\
                                         t11_xxx_xxxx, t11_xxx_xxxy, t11_xxx_xxxz,\
                                         t11_xxx_xxyy, t11_xxx_xxyz, t11_xxx_xxzz,\
                                         t11_xxx_xyyy, t11_xxx_xyyz, t11_xxx_xyzz,\
                                         t11_xxx_xzzz, t11_xxx_yyyy, t11_xxx_yyyz,\
                                         t11_xxx_yyzz, t11_xxx_yzzz, t11_xxx_zzzz,\
                                         t11_xxy_xxxx, t11_xxy_xxxy, t11_xxy_xxxz,\
                                         t11_xxy_xxyy, t11_xxy_xxyz, t11_xxy_xxzz,\
                                         t11_xxy_xyyy, t11_xxy_xyyz, t11_xxy_xyzz,\
                                         t11_xxy_xzzz, t11_xxy_yyyy, t11_xxy_yyyz,\
                                         t11_xxy_yyzz, t11_xxy_yzzz, t11_xxy_zzzz,\
                                         t11_xxz_xxxx, t11_xxz_xxxy, t11_xxz_xxxz,\
                                         t11_xxz_xxyy, t11_xxz_xxyz, t11_xxz_xxzz,\
                                         t11_xxz_xyyy, t11_xxz_xyyz, t11_xxz_xyzz,\
                                         t11_xxz_xzzz, t11_xxz_yyyy, t11_xxz_yyyz,\
                                         t11_xxz_yyzz, t11_xxz_yzzz, t11_xxz_zzzz,\
                                         t11_xyy_xxxx, t11_xyy_xxxy, t11_xyy_xxxz,\
                                         t11_xyy_xxyy, t11_xyy_xxyz, t11_xyy_xxzz,\
                                         t11_xyy_xyyy, t11_xyy_xyyz, t11_xyy_xyzz,\
                                         t11_xyy_xzzz, t11_xyy_yyyy, t11_xyy_yyyz,\
                                         t11_xyy_yyzz, t11_xyy_yzzz, t11_xyy_zzzz,\
                                         t11_xyz_xxxx, t11_xyz_xxxy, t11_xyz_xxxz,\
                                         t11_xyz_xxyy, t11_xyz_xxyz, t11_xyz_xxzz,\
                                         t11_xyz_xyyy, t11_xyz_xyyz, t11_xyz_xyzz,\
                                         t11_xyz_xzzz, t11_xyz_yyyy, t11_xyz_yyyz,\
                                         t11_xyz_yyzz, t11_xyz_yzzz, t11_xyz_zzzz,\
                                         t11_xzz_xxxx, t11_xzz_xxxy, t11_xzz_xxxz,\
                                         t11_xzz_xxyy, t11_xzz_xxyz, t11_xzz_xxzz,\
                                         t11_xzz_xyyy, t11_xzz_xyyz, t11_xzz_xyzz,\
                                         t11_xzz_xzzz, t11_xzz_yyyy, t11_xzz_yyyz,\
                                         t11_xzz_yyzz, t11_xzz_yzzz, t11_xzz_zzzz,\
                                         t11_yyy_xxxx, t11_yyy_xxxy, t11_yyy_xxxz,\
                                         t11_yyy_xxyy, t11_yyy_xxyz, t11_yyy_xxzz,\
                                         t11_yyy_xyyy, t11_yyy_xyyz, t11_yyy_xyzz,\
                                         t11_yyy_xzzz, t11_yyy_yyyy, t11_yyy_yyyz,\
                                         t11_yyy_yyzz, t11_yyy_yzzz, t11_yyy_zzzz,\
                                         t11_yyz_xxxx, t11_yyz_xxxy, t11_yyz_xxxz,\
                                         t11_yyz_xxyy, t11_yyz_xxyz, t11_yyz_xxzz,\
                                         t11_yyz_xyyy, t11_yyz_xyyz, t11_yyz_xyzz,\
                                         t11_yyz_xzzz, t11_yyz_yyyy, t11_yyz_yyyz,\
                                         t11_yyz_yyzz, t11_yyz_yzzz, t11_yyz_zzzz,\
                                         t11_yzz_xxxx, t11_yzz_xxxy, t11_yzz_xxxz,\
                                         t11_yzz_xxyy, t11_yzz_xxyz, t11_yzz_xxzz,\
                                         t11_yzz_xyyy, t11_yzz_xyyz, t11_yzz_xyzz,\
                                         t11_yzz_xzzz, t11_yzz_yyyy, t11_yzz_yyyz,\
                                         t11_yzz_yyzz, t11_yzz_yzzz, t11_yzz_zzzz,\
                                         t11_zzz_xxxx, t11_zzz_xxxy, t11_zzz_xxxz,\
                                         t11_zzz_xxyy, t11_zzz_xxyz, t11_zzz_xxzz,\
                                         t11_zzz_xyyy, t11_zzz_xyyz, t11_zzz_xyzz,\
                                         t11_zzz_xzzz, t11_zzz_yyyy, t11_zzz_yyyz,\
                                         t11_zzz_yyzz, t11_zzz_yzzz, t11_zzz_zzzz,\
                                         t_xxxx_xxxx, t_xxxx_xxxy, t_xxxx_xxxz,\
                                         t_xxxx_xxyy, t_xxxx_xxyz, t_xxxx_xxzz,\
                                         t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz,\
                                         t_xxxx_xzzz, t_xxxx_yyyy, t_xxxx_yyyz,\
                                         t_xxxx_yyzz, t_xxxx_yzzz, t_xxxx_zzzz,\
                                         t_xxxy_xxxx, t_xxxy_xxxy, t_xxxy_xxxz,\
                                         t_xxxy_xxyy, t_xxxy_xxyz, t_xxxy_xxzz,\
                                         t_xxxy_xyyy, t_xxxy_xyyz, t_xxxy_xyzz,\
                                         t_xxxy_xzzz, t_xxxy_yyyy, t_xxxy_yyyz,\
                                         t_xxxy_yyzz, t_xxxy_yzzz, t_xxxy_zzzz,\
                                         t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz,\
                                         t_xxxz_xxyy, t_xxxz_xxyz, t_xxxz_xxzz,\
                                         t_xxxz_xyyy, t_xxxz_xyyz, t_xxxz_xyzz,\
                                         t_xxxz_xzzz, t_xxxz_yyyy, t_xxxz_yyyz,\
                                         t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz,\
                                         t_xxyy_xxxx, t_xxyy_xxxy, t_xxyy_xxxz,\
                                         t_xxyy_xxyy, t_xxyy_xxyz, t_xxyy_xxzz,\
                                         t_xxyy_xyyy, t_xxyy_xyyz, t_xxyy_xyzz,\
                                         t_xxyy_xzzz, t_xxyy_yyyy, t_xxyy_yyyz,\
                                         t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz,\
                                         t_xxyz_xxxx, t_xxyz_xxxy, t_xxyz_xxxz,\
                                         t_xxyz_xxyy, t_xxyz_xxyz, t_xxyz_xxzz,\
                                         t_xxyz_xyyy, t_xxyz_xyyz, t_xxyz_xyzz,\
                                         t_xxyz_xzzz, t_xxyz_yyyy, t_xxyz_yyyz,\
                                         t_xxyz_yyzz, t_xxyz_yzzz, t_xxyz_zzzz,\
                                         t_xxzz_xxxx, t_xxzz_xxxy, t_xxzz_xxxz,\
                                         t_xxzz_xxyy, t_xxzz_xxyz, t_xxzz_xxzz,\
                                         t_xxzz_xyyy, t_xxzz_xyyz, t_xxzz_xyzz,\
                                         t_xxzz_xzzz, t_xxzz_yyyy, t_xxzz_yyyz,\
                                         t_xxzz_yyzz, t_xxzz_yzzz, t_xxzz_zzzz,\
                                         t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz,\
                                         t_xyyy_xxyy, t_xyyy_xxyz, t_xyyy_xxzz,\
                                         t_xyyy_xyyy, t_xyyy_xyyz, t_xyyy_xyzz,\
                                         t_xyyy_xzzz, t_xyyy_yyyy, t_xyyy_yyyz,\
                                         t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz,\
                                         t_xyyz_xxxx, t_xyyz_xxxy, t_xyyz_xxxz,\
                                         t_xyyz_xxyy, t_xyyz_xxyz, t_xyyz_xxzz,\
                                         t_xyyz_xyyy, t_xyyz_xyyz, t_xyyz_xyzz,\
                                         t_xyyz_xzzz, t_xyyz_yyyy, t_xyyz_yyyz,\
                                         t_xyyz_yyzz, t_xyyz_yzzz, t_xyyz_zzzz,\
                                         t_xyzz_xxxx, t_xyzz_xxxy, t_xyzz_xxxz,\
                                         t_xyzz_xxyy, t_xyzz_xxyz, t_xyzz_xxzz,\
                                         t_xyzz_xyyy, t_xyzz_xyyz, t_xyzz_xyzz,\
                                         t_xyzz_xzzz, t_xyzz_yyyy, t_xyzz_yyyz,\
                                         t_xyzz_yyzz, t_xyzz_yzzz, t_xyzz_zzzz,\
                                         t_xzzz_xxxx, t_xzzz_xxxy, t_xzzz_xxxz,\
                                         t_xzzz_xxyy, t_xzzz_xxyz, t_xzzz_xxzz,\
                                         t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz,\
                                         t_xzzz_xzzz, t_xzzz_yyyy, t_xzzz_yyyz,\
                                         t_xzzz_yyzz, t_xzzz_yzzz, t_xzzz_zzzz,\
                                         t_yyyy_xxxx, t_yyyy_xxxy, t_yyyy_xxxz,\
                                         t_yyyy_xxyy, t_yyyy_xxyz, t_yyyy_xxzz,\
                                         t_yyyy_xyyy, t_yyyy_xyyz, t_yyyy_xyzz,\
                                         t_yyyy_xzzz, t_yyyy_yyyy, t_yyyy_yyyz,\
                                         t_yyyy_yyzz, t_yyyy_yzzz, t_yyyy_zzzz,\
                                         t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz,\
                                         t_yyyz_xxyy, t_yyyz_xxyz, t_yyyz_xxzz,\
                                         t_yyyz_xyyy, t_yyyz_xyyz, t_yyyz_xyzz,\
                                         t_yyyz_xzzz, t_yyyz_yyyy, t_yyyz_yyyz,\
                                         t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz,\
                                         t_yyzz_xxxx, t_yyzz_xxxy, t_yyzz_xxxz,\
                                         t_yyzz_xxyy, t_yyzz_xxyz, t_yyzz_xxzz,\
                                         t_yyzz_xyyy, t_yyzz_xyyz, t_yyzz_xyzz,\
                                         t_yyzz_xzzz, t_yyzz_yyyy, t_yyzz_yyyz,\
                                         t_yyzz_yyzz, t_yyzz_yzzz, t_yyzz_zzzz,\
                                         t_yzzz_xxxx, t_yzzz_xxxy, t_yzzz_xxxz,\
                                         t_yzzz_xxyy, t_yzzz_xxyz, t_yzzz_xxzz,\
                                         t_yzzz_xyyy, t_yzzz_xyyz, t_yzzz_xyzz,\
                                         t_yzzz_xzzz, t_yzzz_yyyy, t_yzzz_yyyz,\
                                         t_yzzz_yyzz, t_yzzz_yzzz, t_yzzz_zzzz,\
                                         t_zzzz_xxxx, t_zzzz_xxxy, t_zzzz_xxxz,\
                                         t_zzzz_xxyy, t_zzzz_xxyz, t_zzzz_xxzz,\
                                         t_zzzz_xyyy, t_zzzz_xyyz, t_zzzz_xyzz,\
                                         t_zzzz_xzzz, t_zzzz_yyyy, t_zzzz_yyyz,\
                                         t_zzzz_yyzz, t_zzzz_yzzz, t_zzzz_zzzz,\
                                         pcx, pcy, pcz: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    // scaled prefactor

                    double f2t = 0.50 * fx[k];

                    // leading x component

                    double fra = pax[k];

                    double frc = pcx[k];

                    t_xxxx_xxxx[k] = fra * t10_xxx_xxxx[k] - frc * t11_xxx_xxxx[k] + f2t * (3.0 * t20_xx_xxxx[k] - 3.0 * t21_xx_xxxx[k] + 4.0 * tk0_xxx_xxx[k] - 4.0 * tk1_xxx_xxx[k]);

                    t_xxxx_xxxy[k] = fra * t10_xxx_xxxy[k] - frc * t11_xxx_xxxy[k] + f2t * (3.0 * t20_xx_xxxy[k] - 3.0 * t21_xx_xxxy[k] + 3.0 * tk0_xxx_xxy[k] - 3.0 * tk1_xxx_xxy[k]);

                    t_xxxx_xxxz[k] = fra * t10_xxx_xxxz[k] - frc * t11_xxx_xxxz[k] + f2t * (3.0 * t20_xx_xxxz[k] - 3.0 * t21_xx_xxxz[k] + 3.0 * tk0_xxx_xxz[k] - 3.0 * tk1_xxx_xxz[k]);

                    t_xxxx_xxyy[k] = fra * t10_xxx_xxyy[k] - frc * t11_xxx_xxyy[k] + f2t * (3.0 * t20_xx_xxyy[k] - 3.0 * t21_xx_xxyy[k] + 2.0 * tk0_xxx_xyy[k] - 2.0 * tk1_xxx_xyy[k]);

                    t_xxxx_xxyz[k] = fra * t10_xxx_xxyz[k] - frc * t11_xxx_xxyz[k] + f2t * (3.0 * t20_xx_xxyz[k] - 3.0 * t21_xx_xxyz[k] + 2.0 * tk0_xxx_xyz[k] - 2.0 * tk1_xxx_xyz[k]);

                    t_xxxx_xxzz[k] = fra * t10_xxx_xxzz[k] - frc * t11_xxx_xxzz[k] + f2t * (3.0 * t20_xx_xxzz[k] - 3.0 * t21_xx_xxzz[k] + 2.0 * tk0_xxx_xzz[k] - 2.0 * tk1_xxx_xzz[k]);

                    t_xxxx_xyyy[k] = fra * t10_xxx_xyyy[k] - frc * t11_xxx_xyyy[k] + f2t * (3.0 * t20_xx_xyyy[k] - 3.0 * t21_xx_xyyy[k] + tk0_xxx_yyy[k] - tk1_xxx_yyy[k]);

                    t_xxxx_xyyz[k] = fra * t10_xxx_xyyz[k] - frc * t11_xxx_xyyz[k] + f2t * (3.0 * t20_xx_xyyz[k] - 3.0 * t21_xx_xyyz[k] + tk0_xxx_yyz[k] - tk1_xxx_yyz[k]);

                    t_xxxx_xyzz[k] = fra * t10_xxx_xyzz[k] - frc * t11_xxx_xyzz[k] + f2t * (3.0 * t20_xx_xyzz[k] - 3.0 * t21_xx_xyzz[k] + tk0_xxx_yzz[k] - tk1_xxx_yzz[k]);

                    t_xxxx_xzzz[k] = fra * t10_xxx_xzzz[k] - frc * t11_xxx_xzzz[k] + f2t * (3.0 * t20_xx_xzzz[k] - 3.0 * t21_xx_xzzz[k] + tk0_xxx_zzz[k] - tk1_xxx_zzz[k]);

                    t_xxxx_yyyy[k] = fra * t10_xxx_yyyy[k] - frc * t11_xxx_yyyy[k] + f2t * (3.0 * t20_xx_yyyy[k] - 3.0 * t21_xx_yyyy[k]);

                    t_xxxx_yyyz[k] = fra * t10_xxx_yyyz[k] - frc * t11_xxx_yyyz[k] + f2t * (3.0 * t20_xx_yyyz[k] - 3.0 * t21_xx_yyyz[k]);

                    t_xxxx_yyzz[k] = fra * t10_xxx_yyzz[k] - frc * t11_xxx_yyzz[k] + f2t * (3.0 * t20_xx_yyzz[k] - 3.0 * t21_xx_yyzz[k]);

                    t_xxxx_yzzz[k] = fra * t10_xxx_yzzz[k] - frc * t11_xxx_yzzz[k] + f2t * (3.0 * t20_xx_yzzz[k] - 3.0 * t21_xx_yzzz[k]);

                    t_xxxx_zzzz[k] = fra * t10_xxx_zzzz[k] - frc * t11_xxx_zzzz[k] + f2t * (3.0 * t20_xx_zzzz[k] - 3.0 * t21_xx_zzzz[k]);

                    t_xxxy_xxxx[k] = fra * t10_xxy_xxxx[k] - frc * t11_xxy_xxxx[k] + f2t * (2.0 * t20_xy_xxxx[k] - 2.0 * t21_xy_xxxx[k] + 4.0 * tk0_xxy_xxx[k] - 4.0 * tk1_xxy_xxx[k]);

                    t_xxxy_xxxy[k] = fra * t10_xxy_xxxy[k] - frc * t11_xxy_xxxy[k] + f2t * (2.0 * t20_xy_xxxy[k] - 2.0 * t21_xy_xxxy[k] + 3.0 * tk0_xxy_xxy[k] - 3.0 * tk1_xxy_xxy[k]);

                    t_xxxy_xxxz[k] = fra * t10_xxy_xxxz[k] - frc * t11_xxy_xxxz[k] + f2t * (2.0 * t20_xy_xxxz[k] - 2.0 * t21_xy_xxxz[k] + 3.0 * tk0_xxy_xxz[k] - 3.0 * tk1_xxy_xxz[k]);

                    t_xxxy_xxyy[k] = fra * t10_xxy_xxyy[k] - frc * t11_xxy_xxyy[k] + f2t * (2.0 * t20_xy_xxyy[k] - 2.0 * t21_xy_xxyy[k] + 2.0 * tk0_xxy_xyy[k] - 2.0 * tk1_xxy_xyy[k]);

                    t_xxxy_xxyz[k] = fra * t10_xxy_xxyz[k] - frc * t11_xxy_xxyz[k] + f2t * (2.0 * t20_xy_xxyz[k] - 2.0 * t21_xy_xxyz[k] + 2.0 * tk0_xxy_xyz[k] - 2.0 * tk1_xxy_xyz[k]);

                    t_xxxy_xxzz[k] = fra * t10_xxy_xxzz[k] - frc * t11_xxy_xxzz[k] + f2t * (2.0 * t20_xy_xxzz[k] - 2.0 * t21_xy_xxzz[k] + 2.0 * tk0_xxy_xzz[k] - 2.0 * tk1_xxy_xzz[k]);

                    t_xxxy_xyyy[k] = fra * t10_xxy_xyyy[k] - frc * t11_xxy_xyyy[k] + f2t * (2.0 * t20_xy_xyyy[k] - 2.0 * t21_xy_xyyy[k] + tk0_xxy_yyy[k] - tk1_xxy_yyy[k]);

                    t_xxxy_xyyz[k] = fra * t10_xxy_xyyz[k] - frc * t11_xxy_xyyz[k] + f2t * (2.0 * t20_xy_xyyz[k] - 2.0 * t21_xy_xyyz[k] + tk0_xxy_yyz[k] - tk1_xxy_yyz[k]);

                    t_xxxy_xyzz[k] = fra * t10_xxy_xyzz[k] - frc * t11_xxy_xyzz[k] + f2t * (2.0 * t20_xy_xyzz[k] - 2.0 * t21_xy_xyzz[k] + tk0_xxy_yzz[k] - tk1_xxy_yzz[k]);

                    t_xxxy_xzzz[k] = fra * t10_xxy_xzzz[k] - frc * t11_xxy_xzzz[k] + f2t * (2.0 * t20_xy_xzzz[k] - 2.0 * t21_xy_xzzz[k] + tk0_xxy_zzz[k] - tk1_xxy_zzz[k]);

                    t_xxxy_yyyy[k] = fra * t10_xxy_yyyy[k] - frc * t11_xxy_yyyy[k] + f2t * (2.0 * t20_xy_yyyy[k] - 2.0 * t21_xy_yyyy[k]);

                    t_xxxy_yyyz[k] = fra * t10_xxy_yyyz[k] - frc * t11_xxy_yyyz[k] + f2t * (2.0 * t20_xy_yyyz[k] - 2.0 * t21_xy_yyyz[k]);

                    t_xxxy_yyzz[k] = fra * t10_xxy_yyzz[k] - frc * t11_xxy_yyzz[k] + f2t * (2.0 * t20_xy_yyzz[k] - 2.0 * t21_xy_yyzz[k]);

                    t_xxxy_yzzz[k] = fra * t10_xxy_yzzz[k] - frc * t11_xxy_yzzz[k] + f2t * (2.0 * t20_xy_yzzz[k] - 2.0 * t21_xy_yzzz[k]);

                    t_xxxy_zzzz[k] = fra * t10_xxy_zzzz[k] - frc * t11_xxy_zzzz[k] + f2t * (2.0 * t20_xy_zzzz[k] - 2.0 * t21_xy_zzzz[k]);

                    t_xxxz_xxxx[k] = fra * t10_xxz_xxxx[k] - frc * t11_xxz_xxxx[k] + f2t * (2.0 * t20_xz_xxxx[k] - 2.0 * t21_xz_xxxx[k] + 4.0 * tk0_xxz_xxx[k] - 4.0 * tk1_xxz_xxx[k]);

                    t_xxxz_xxxy[k] = fra * t10_xxz_xxxy[k] - frc * t11_xxz_xxxy[k] + f2t * (2.0 * t20_xz_xxxy[k] - 2.0 * t21_xz_xxxy[k] + 3.0 * tk0_xxz_xxy[k] - 3.0 * tk1_xxz_xxy[k]);

                    t_xxxz_xxxz[k] = fra * t10_xxz_xxxz[k] - frc * t11_xxz_xxxz[k] + f2t * (2.0 * t20_xz_xxxz[k] - 2.0 * t21_xz_xxxz[k] + 3.0 * tk0_xxz_xxz[k] - 3.0 * tk1_xxz_xxz[k]);

                    t_xxxz_xxyy[k] = fra * t10_xxz_xxyy[k] - frc * t11_xxz_xxyy[k] + f2t * (2.0 * t20_xz_xxyy[k] - 2.0 * t21_xz_xxyy[k] + 2.0 * tk0_xxz_xyy[k] - 2.0 * tk1_xxz_xyy[k]);

                    t_xxxz_xxyz[k] = fra * t10_xxz_xxyz[k] - frc * t11_xxz_xxyz[k] + f2t * (2.0 * t20_xz_xxyz[k] - 2.0 * t21_xz_xxyz[k] + 2.0 * tk0_xxz_xyz[k] - 2.0 * tk1_xxz_xyz[k]);

                    t_xxxz_xxzz[k] = fra * t10_xxz_xxzz[k] - frc * t11_xxz_xxzz[k] + f2t * (2.0 * t20_xz_xxzz[k] - 2.0 * t21_xz_xxzz[k] + 2.0 * tk0_xxz_xzz[k] - 2.0 * tk1_xxz_xzz[k]);

                    t_xxxz_xyyy[k] = fra * t10_xxz_xyyy[k] - frc * t11_xxz_xyyy[k] + f2t * (2.0 * t20_xz_xyyy[k] - 2.0 * t21_xz_xyyy[k] + tk0_xxz_yyy[k] - tk1_xxz_yyy[k]);

                    t_xxxz_xyyz[k] = fra * t10_xxz_xyyz[k] - frc * t11_xxz_xyyz[k] + f2t * (2.0 * t20_xz_xyyz[k] - 2.0 * t21_xz_xyyz[k] + tk0_xxz_yyz[k] - tk1_xxz_yyz[k]);

                    t_xxxz_xyzz[k] = fra * t10_xxz_xyzz[k] - frc * t11_xxz_xyzz[k] + f2t * (2.0 * t20_xz_xyzz[k] - 2.0 * t21_xz_xyzz[k] + tk0_xxz_yzz[k] - tk1_xxz_yzz[k]);

                    t_xxxz_xzzz[k] = fra * t10_xxz_xzzz[k] - frc * t11_xxz_xzzz[k] + f2t * (2.0 * t20_xz_xzzz[k] - 2.0 * t21_xz_xzzz[k] + tk0_xxz_zzz[k] - tk1_xxz_zzz[k]);

                    t_xxxz_yyyy[k] = fra * t10_xxz_yyyy[k] - frc * t11_xxz_yyyy[k] + f2t * (2.0 * t20_xz_yyyy[k] - 2.0 * t21_xz_yyyy[k]);

                    t_xxxz_yyyz[k] = fra * t10_xxz_yyyz[k] - frc * t11_xxz_yyyz[k] + f2t * (2.0 * t20_xz_yyyz[k] - 2.0 * t21_xz_yyyz[k]);

                    t_xxxz_yyzz[k] = fra * t10_xxz_yyzz[k] - frc * t11_xxz_yyzz[k] + f2t * (2.0 * t20_xz_yyzz[k] - 2.0 * t21_xz_yyzz[k]);

                    t_xxxz_yzzz[k] = fra * t10_xxz_yzzz[k] - frc * t11_xxz_yzzz[k] + f2t * (2.0 * t20_xz_yzzz[k] - 2.0 * t21_xz_yzzz[k]);

                    t_xxxz_zzzz[k] = fra * t10_xxz_zzzz[k] - frc * t11_xxz_zzzz[k] + f2t * (2.0 * t20_xz_zzzz[k] - 2.0 * t21_xz_zzzz[k]);

                    t_xxyy_xxxx[k] = fra * t10_xyy_xxxx[k] - frc * t11_xyy_xxxx[k] + f2t * (t20_yy_xxxx[k] - t21_yy_xxxx[k] + 4.0 * tk0_xyy_xxx[k] - 4.0 * tk1_xyy_xxx[k]);

                    t_xxyy_xxxy[k] = fra * t10_xyy_xxxy[k] - frc * t11_xyy_xxxy[k] + f2t * (t20_yy_xxxy[k] - t21_yy_xxxy[k] + 3.0 * tk0_xyy_xxy[k] - 3.0 * tk1_xyy_xxy[k]);

                    t_xxyy_xxxz[k] = fra * t10_xyy_xxxz[k] - frc * t11_xyy_xxxz[k] + f2t * (t20_yy_xxxz[k] - t21_yy_xxxz[k] + 3.0 * tk0_xyy_xxz[k] - 3.0 * tk1_xyy_xxz[k]);

                    t_xxyy_xxyy[k] = fra * t10_xyy_xxyy[k] - frc * t11_xyy_xxyy[k] + f2t * (t20_yy_xxyy[k] - t21_yy_xxyy[k] + 2.0 * tk0_xyy_xyy[k] - 2.0 * tk1_xyy_xyy[k]);

                    t_xxyy_xxyz[k] = fra * t10_xyy_xxyz[k] - frc * t11_xyy_xxyz[k] + f2t * (t20_yy_xxyz[k] - t21_yy_xxyz[k] + 2.0 * tk0_xyy_xyz[k] - 2.0 * tk1_xyy_xyz[k]);

                    t_xxyy_xxzz[k] = fra * t10_xyy_xxzz[k] - frc * t11_xyy_xxzz[k] + f2t * (t20_yy_xxzz[k] - t21_yy_xxzz[k] + 2.0 * tk0_xyy_xzz[k] - 2.0 * tk1_xyy_xzz[k]);

                    t_xxyy_xyyy[k] = fra * t10_xyy_xyyy[k] - frc * t11_xyy_xyyy[k] + f2t * (t20_yy_xyyy[k] - t21_yy_xyyy[k] + tk0_xyy_yyy[k] - tk1_xyy_yyy[k]);

                    t_xxyy_xyyz[k] = fra * t10_xyy_xyyz[k] - frc * t11_xyy_xyyz[k] + f2t * (t20_yy_xyyz[k] - t21_yy_xyyz[k] + tk0_xyy_yyz[k] - tk1_xyy_yyz[k]);

                    t_xxyy_xyzz[k] = fra * t10_xyy_xyzz[k] - frc * t11_xyy_xyzz[k] + f2t * (t20_yy_xyzz[k] - t21_yy_xyzz[k] + tk0_xyy_yzz[k] - tk1_xyy_yzz[k]);

                    t_xxyy_xzzz[k] = fra * t10_xyy_xzzz[k] - frc * t11_xyy_xzzz[k] + f2t * (t20_yy_xzzz[k] - t21_yy_xzzz[k] + tk0_xyy_zzz[k] - tk1_xyy_zzz[k]);

                    t_xxyy_yyyy[k] = fra * t10_xyy_yyyy[k] - frc * t11_xyy_yyyy[k] + f2t * (t20_yy_yyyy[k] - t21_yy_yyyy[k]);

                    t_xxyy_yyyz[k] = fra * t10_xyy_yyyz[k] - frc * t11_xyy_yyyz[k] + f2t * (t20_yy_yyyz[k] - t21_yy_yyyz[k]);

                    t_xxyy_yyzz[k] = fra * t10_xyy_yyzz[k] - frc * t11_xyy_yyzz[k] + f2t * (t20_yy_yyzz[k] - t21_yy_yyzz[k]);

                    t_xxyy_yzzz[k] = fra * t10_xyy_yzzz[k] - frc * t11_xyy_yzzz[k] + f2t * (t20_yy_yzzz[k] - t21_yy_yzzz[k]);

                    t_xxyy_zzzz[k] = fra * t10_xyy_zzzz[k] - frc * t11_xyy_zzzz[k] + f2t * (t20_yy_zzzz[k] - t21_yy_zzzz[k]);

                    t_xxyz_xxxx[k] = fra * t10_xyz_xxxx[k] - frc * t11_xyz_xxxx[k] + f2t * (t20_yz_xxxx[k] - t21_yz_xxxx[k] + 4.0 * tk0_xyz_xxx[k] - 4.0 * tk1_xyz_xxx[k]);

                    t_xxyz_xxxy[k] = fra * t10_xyz_xxxy[k] - frc * t11_xyz_xxxy[k] + f2t * (t20_yz_xxxy[k] - t21_yz_xxxy[k] + 3.0 * tk0_xyz_xxy[k] - 3.0 * tk1_xyz_xxy[k]);

                    t_xxyz_xxxz[k] = fra * t10_xyz_xxxz[k] - frc * t11_xyz_xxxz[k] + f2t * (t20_yz_xxxz[k] - t21_yz_xxxz[k] + 3.0 * tk0_xyz_xxz[k] - 3.0 * tk1_xyz_xxz[k]);

                    t_xxyz_xxyy[k] = fra * t10_xyz_xxyy[k] - frc * t11_xyz_xxyy[k] + f2t * (t20_yz_xxyy[k] - t21_yz_xxyy[k] + 2.0 * tk0_xyz_xyy[k] - 2.0 * tk1_xyz_xyy[k]);

                    t_xxyz_xxyz[k] = fra * t10_xyz_xxyz[k] - frc * t11_xyz_xxyz[k] + f2t * (t20_yz_xxyz[k] - t21_yz_xxyz[k] + 2.0 * tk0_xyz_xyz[k] - 2.0 * tk1_xyz_xyz[k]);

                    t_xxyz_xxzz[k] = fra * t10_xyz_xxzz[k] - frc * t11_xyz_xxzz[k] + f2t * (t20_yz_xxzz[k] - t21_yz_xxzz[k] + 2.0 * tk0_xyz_xzz[k] - 2.0 * tk1_xyz_xzz[k]);

                    t_xxyz_xyyy[k] = fra * t10_xyz_xyyy[k] - frc * t11_xyz_xyyy[k] + f2t * (t20_yz_xyyy[k] - t21_yz_xyyy[k] + tk0_xyz_yyy[k] - tk1_xyz_yyy[k]);

                    t_xxyz_xyyz[k] = fra * t10_xyz_xyyz[k] - frc * t11_xyz_xyyz[k] + f2t * (t20_yz_xyyz[k] - t21_yz_xyyz[k] + tk0_xyz_yyz[k] - tk1_xyz_yyz[k]);

                    t_xxyz_xyzz[k] = fra * t10_xyz_xyzz[k] - frc * t11_xyz_xyzz[k] + f2t * (t20_yz_xyzz[k] - t21_yz_xyzz[k] + tk0_xyz_yzz[k] - tk1_xyz_yzz[k]);

                    t_xxyz_xzzz[k] = fra * t10_xyz_xzzz[k] - frc * t11_xyz_xzzz[k] + f2t * (t20_yz_xzzz[k] - t21_yz_xzzz[k] + tk0_xyz_zzz[k] - tk1_xyz_zzz[k]);

                    t_xxyz_yyyy[k] = fra * t10_xyz_yyyy[k] - frc * t11_xyz_yyyy[k] + f2t * (t20_yz_yyyy[k] - t21_yz_yyyy[k]);

                    t_xxyz_yyyz[k] = fra * t10_xyz_yyyz[k] - frc * t11_xyz_yyyz[k] + f2t * (t20_yz_yyyz[k] - t21_yz_yyyz[k]);

                    t_xxyz_yyzz[k] = fra * t10_xyz_yyzz[k] - frc * t11_xyz_yyzz[k] + f2t * (t20_yz_yyzz[k] - t21_yz_yyzz[k]);

                    t_xxyz_yzzz[k] = fra * t10_xyz_yzzz[k] - frc * t11_xyz_yzzz[k] + f2t * (t20_yz_yzzz[k] - t21_yz_yzzz[k]);

                    t_xxyz_zzzz[k] = fra * t10_xyz_zzzz[k] - frc * t11_xyz_zzzz[k] + f2t * (t20_yz_zzzz[k] - t21_yz_zzzz[k]);

                    t_xxzz_xxxx[k] = fra * t10_xzz_xxxx[k] - frc * t11_xzz_xxxx[k] + f2t * (t20_zz_xxxx[k] - t21_zz_xxxx[k] + 4.0 * tk0_xzz_xxx[k] - 4.0 * tk1_xzz_xxx[k]);

                    t_xxzz_xxxy[k] = fra * t10_xzz_xxxy[k] - frc * t11_xzz_xxxy[k] + f2t * (t20_zz_xxxy[k] - t21_zz_xxxy[k] + 3.0 * tk0_xzz_xxy[k] - 3.0 * tk1_xzz_xxy[k]);

                    t_xxzz_xxxz[k] = fra * t10_xzz_xxxz[k] - frc * t11_xzz_xxxz[k] + f2t * (t20_zz_xxxz[k] - t21_zz_xxxz[k] + 3.0 * tk0_xzz_xxz[k] - 3.0 * tk1_xzz_xxz[k]);

                    t_xxzz_xxyy[k] = fra * t10_xzz_xxyy[k] - frc * t11_xzz_xxyy[k] + f2t * (t20_zz_xxyy[k] - t21_zz_xxyy[k] + 2.0 * tk0_xzz_xyy[k] - 2.0 * tk1_xzz_xyy[k]);

                    t_xxzz_xxyz[k] = fra * t10_xzz_xxyz[k] - frc * t11_xzz_xxyz[k] + f2t * (t20_zz_xxyz[k] - t21_zz_xxyz[k] + 2.0 * tk0_xzz_xyz[k] - 2.0 * tk1_xzz_xyz[k]);

                    t_xxzz_xxzz[k] = fra * t10_xzz_xxzz[k] - frc * t11_xzz_xxzz[k] + f2t * (t20_zz_xxzz[k] - t21_zz_xxzz[k] + 2.0 * tk0_xzz_xzz[k] - 2.0 * tk1_xzz_xzz[k]);

                    t_xxzz_xyyy[k] = fra * t10_xzz_xyyy[k] - frc * t11_xzz_xyyy[k] + f2t * (t20_zz_xyyy[k] - t21_zz_xyyy[k] + tk0_xzz_yyy[k] - tk1_xzz_yyy[k]);

                    t_xxzz_xyyz[k] = fra * t10_xzz_xyyz[k] - frc * t11_xzz_xyyz[k] + f2t * (t20_zz_xyyz[k] - t21_zz_xyyz[k] + tk0_xzz_yyz[k] - tk1_xzz_yyz[k]);

                    t_xxzz_xyzz[k] = fra * t10_xzz_xyzz[k] - frc * t11_xzz_xyzz[k] + f2t * (t20_zz_xyzz[k] - t21_zz_xyzz[k] + tk0_xzz_yzz[k] - tk1_xzz_yzz[k]);

                    t_xxzz_xzzz[k] = fra * t10_xzz_xzzz[k] - frc * t11_xzz_xzzz[k] + f2t * (t20_zz_xzzz[k] - t21_zz_xzzz[k] + tk0_xzz_zzz[k] - tk1_xzz_zzz[k]);

                    t_xxzz_yyyy[k] = fra * t10_xzz_yyyy[k] - frc * t11_xzz_yyyy[k] + f2t * (t20_zz_yyyy[k] - t21_zz_yyyy[k]);

                    t_xxzz_yyyz[k] = fra * t10_xzz_yyyz[k] - frc * t11_xzz_yyyz[k] + f2t * (t20_zz_yyyz[k] - t21_zz_yyyz[k]);

                    t_xxzz_yyzz[k] = fra * t10_xzz_yyzz[k] - frc * t11_xzz_yyzz[k] + f2t * (t20_zz_yyzz[k] - t21_zz_yyzz[k]);

                    t_xxzz_yzzz[k] = fra * t10_xzz_yzzz[k] - frc * t11_xzz_yzzz[k] + f2t * (t20_zz_yzzz[k] - t21_zz_yzzz[k]);

                    t_xxzz_zzzz[k] = fra * t10_xzz_zzzz[k] - frc * t11_xzz_zzzz[k] + f2t * (t20_zz_zzzz[k] - t21_zz_zzzz[k]);

                    t_xyyy_xxxx[k] = fra * t10_yyy_xxxx[k] - frc * t11_yyy_xxxx[k] + f2t * (4.0 * tk0_yyy_xxx[k] - 4.0 * tk1_yyy_xxx[k]);

                    t_xyyy_xxxy[k] = fra * t10_yyy_xxxy[k] - frc * t11_yyy_xxxy[k] + f2t * (3.0 * tk0_yyy_xxy[k] - 3.0 * tk1_yyy_xxy[k]);

                    t_xyyy_xxxz[k] = fra * t10_yyy_xxxz[k] - frc * t11_yyy_xxxz[k] + f2t * (3.0 * tk0_yyy_xxz[k] - 3.0 * tk1_yyy_xxz[k]);

                    t_xyyy_xxyy[k] = fra * t10_yyy_xxyy[k] - frc * t11_yyy_xxyy[k] + f2t * (2.0 * tk0_yyy_xyy[k] - 2.0 * tk1_yyy_xyy[k]);

                    t_xyyy_xxyz[k] = fra * t10_yyy_xxyz[k] - frc * t11_yyy_xxyz[k] + f2t * (2.0 * tk0_yyy_xyz[k] - 2.0 * tk1_yyy_xyz[k]);

                    t_xyyy_xxzz[k] = fra * t10_yyy_xxzz[k] - frc * t11_yyy_xxzz[k] + f2t * (2.0 * tk0_yyy_xzz[k] - 2.0 * tk1_yyy_xzz[k]);

                    t_xyyy_xyyy[k] = fra * t10_yyy_xyyy[k] - frc * t11_yyy_xyyy[k] + f2t * (tk0_yyy_yyy[k] - tk1_yyy_yyy[k]);

                    t_xyyy_xyyz[k] = fra * t10_yyy_xyyz[k] - frc * t11_yyy_xyyz[k] + f2t * (tk0_yyy_yyz[k] - tk1_yyy_yyz[k]);

                    t_xyyy_xyzz[k] = fra * t10_yyy_xyzz[k] - frc * t11_yyy_xyzz[k] + f2t * (tk0_yyy_yzz[k] - tk1_yyy_yzz[k]);

                    t_xyyy_xzzz[k] = fra * t10_yyy_xzzz[k] - frc * t11_yyy_xzzz[k] + f2t * (tk0_yyy_zzz[k] - tk1_yyy_zzz[k]);

                    t_xyyy_yyyy[k] = fra * t10_yyy_yyyy[k] - frc * t11_yyy_yyyy[k];

                    t_xyyy_yyyz[k] = fra * t10_yyy_yyyz[k] - frc * t11_yyy_yyyz[k];

                    t_xyyy_yyzz[k] = fra * t10_yyy_yyzz[k] - frc * t11_yyy_yyzz[k];

                    t_xyyy_yzzz[k] = fra * t10_yyy_yzzz[k] - frc * t11_yyy_yzzz[k];

                    t_xyyy_zzzz[k] = fra * t10_yyy_zzzz[k] - frc * t11_yyy_zzzz[k];

                    t_xyyz_xxxx[k] = fra * t10_yyz_xxxx[k] - frc * t11_yyz_xxxx[k] + f2t * (4.0 * tk0_yyz_xxx[k] - 4.0 * tk1_yyz_xxx[k]);

                    t_xyyz_xxxy[k] = fra * t10_yyz_xxxy[k] - frc * t11_yyz_xxxy[k] + f2t * (3.0 * tk0_yyz_xxy[k] - 3.0 * tk1_yyz_xxy[k]);

                    t_xyyz_xxxz[k] = fra * t10_yyz_xxxz[k] - frc * t11_yyz_xxxz[k] + f2t * (3.0 * tk0_yyz_xxz[k] - 3.0 * tk1_yyz_xxz[k]);

                    t_xyyz_xxyy[k] = fra * t10_yyz_xxyy[k] - frc * t11_yyz_xxyy[k] + f2t * (2.0 * tk0_yyz_xyy[k] - 2.0 * tk1_yyz_xyy[k]);

                    t_xyyz_xxyz[k] = fra * t10_yyz_xxyz[k] - frc * t11_yyz_xxyz[k] + f2t * (2.0 * tk0_yyz_xyz[k] - 2.0 * tk1_yyz_xyz[k]);

                    t_xyyz_xxzz[k] = fra * t10_yyz_xxzz[k] - frc * t11_yyz_xxzz[k] + f2t * (2.0 * tk0_yyz_xzz[k] - 2.0 * tk1_yyz_xzz[k]);

                    t_xyyz_xyyy[k] = fra * t10_yyz_xyyy[k] - frc * t11_yyz_xyyy[k] + f2t * (tk0_yyz_yyy[k] - tk1_yyz_yyy[k]);

                    t_xyyz_xyyz[k] = fra * t10_yyz_xyyz[k] - frc * t11_yyz_xyyz[k] + f2t * (tk0_yyz_yyz[k] - tk1_yyz_yyz[k]);

                    t_xyyz_xyzz[k] = fra * t10_yyz_xyzz[k] - frc * t11_yyz_xyzz[k] + f2t * (tk0_yyz_yzz[k] - tk1_yyz_yzz[k]);

                    t_xyyz_xzzz[k] = fra * t10_yyz_xzzz[k] - frc * t11_yyz_xzzz[k] + f2t * (tk0_yyz_zzz[k] - tk1_yyz_zzz[k]);

                    t_xyyz_yyyy[k] = fra * t10_yyz_yyyy[k] - frc * t11_yyz_yyyy[k];

                    t_xyyz_yyyz[k] = fra * t10_yyz_yyyz[k] - frc * t11_yyz_yyyz[k];

                    t_xyyz_yyzz[k] = fra * t10_yyz_yyzz[k] - frc * t11_yyz_yyzz[k];

                    t_xyyz_yzzz[k] = fra * t10_yyz_yzzz[k] - frc * t11_yyz_yzzz[k];

                    t_xyyz_zzzz[k] = fra * t10_yyz_zzzz[k] - frc * t11_yyz_zzzz[k];

                    t_xyzz_xxxx[k] = fra * t10_yzz_xxxx[k] - frc * t11_yzz_xxxx[k] + f2t * (4.0 * tk0_yzz_xxx[k] - 4.0 * tk1_yzz_xxx[k]);

                    t_xyzz_xxxy[k] = fra * t10_yzz_xxxy[k] - frc * t11_yzz_xxxy[k] + f2t * (3.0 * tk0_yzz_xxy[k] - 3.0 * tk1_yzz_xxy[k]);

                    t_xyzz_xxxz[k] = fra * t10_yzz_xxxz[k] - frc * t11_yzz_xxxz[k] + f2t * (3.0 * tk0_yzz_xxz[k] - 3.0 * tk1_yzz_xxz[k]);

                    t_xyzz_xxyy[k] = fra * t10_yzz_xxyy[k] - frc * t11_yzz_xxyy[k] + f2t * (2.0 * tk0_yzz_xyy[k] - 2.0 * tk1_yzz_xyy[k]);

                    t_xyzz_xxyz[k] = fra * t10_yzz_xxyz[k] - frc * t11_yzz_xxyz[k] + f2t * (2.0 * tk0_yzz_xyz[k] - 2.0 * tk1_yzz_xyz[k]);

                    t_xyzz_xxzz[k] = fra * t10_yzz_xxzz[k] - frc * t11_yzz_xxzz[k] + f2t * (2.0 * tk0_yzz_xzz[k] - 2.0 * tk1_yzz_xzz[k]);

                    t_xyzz_xyyy[k] = fra * t10_yzz_xyyy[k] - frc * t11_yzz_xyyy[k] + f2t * (tk0_yzz_yyy[k] - tk1_yzz_yyy[k]);

                    t_xyzz_xyyz[k] = fra * t10_yzz_xyyz[k] - frc * t11_yzz_xyyz[k] + f2t * (tk0_yzz_yyz[k] - tk1_yzz_yyz[k]);

                    t_xyzz_xyzz[k] = fra * t10_yzz_xyzz[k] - frc * t11_yzz_xyzz[k] + f2t * (tk0_yzz_yzz[k] - tk1_yzz_yzz[k]);

                    t_xyzz_xzzz[k] = fra * t10_yzz_xzzz[k] - frc * t11_yzz_xzzz[k] + f2t * (tk0_yzz_zzz[k] - tk1_yzz_zzz[k]);

                    t_xyzz_yyyy[k] = fra * t10_yzz_yyyy[k] - frc * t11_yzz_yyyy[k];

                    t_xyzz_yyyz[k] = fra * t10_yzz_yyyz[k] - frc * t11_yzz_yyyz[k];

                    t_xyzz_yyzz[k] = fra * t10_yzz_yyzz[k] - frc * t11_yzz_yyzz[k];

                    t_xyzz_yzzz[k] = fra * t10_yzz_yzzz[k] - frc * t11_yzz_yzzz[k];

                    t_xyzz_zzzz[k] = fra * t10_yzz_zzzz[k] - frc * t11_yzz_zzzz[k];

                    t_xzzz_xxxx[k] = fra * t10_zzz_xxxx[k] - frc * t11_zzz_xxxx[k] + f2t * (4.0 * tk0_zzz_xxx[k] - 4.0 * tk1_zzz_xxx[k]);

                    t_xzzz_xxxy[k] = fra * t10_zzz_xxxy[k] - frc * t11_zzz_xxxy[k] + f2t * (3.0 * tk0_zzz_xxy[k] - 3.0 * tk1_zzz_xxy[k]);

                    t_xzzz_xxxz[k] = fra * t10_zzz_xxxz[k] - frc * t11_zzz_xxxz[k] + f2t * (3.0 * tk0_zzz_xxz[k] - 3.0 * tk1_zzz_xxz[k]);

                    t_xzzz_xxyy[k] = fra * t10_zzz_xxyy[k] - frc * t11_zzz_xxyy[k] + f2t * (2.0 * tk0_zzz_xyy[k] - 2.0 * tk1_zzz_xyy[k]);

                    t_xzzz_xxyz[k] = fra * t10_zzz_xxyz[k] - frc * t11_zzz_xxyz[k] + f2t * (2.0 * tk0_zzz_xyz[k] - 2.0 * tk1_zzz_xyz[k]);

                    t_xzzz_xxzz[k] = fra * t10_zzz_xxzz[k] - frc * t11_zzz_xxzz[k] + f2t * (2.0 * tk0_zzz_xzz[k] - 2.0 * tk1_zzz_xzz[k]);

                    t_xzzz_xyyy[k] = fra * t10_zzz_xyyy[k] - frc * t11_zzz_xyyy[k] + f2t * (tk0_zzz_yyy[k] - tk1_zzz_yyy[k]);

                    t_xzzz_xyyz[k] = fra * t10_zzz_xyyz[k] - frc * t11_zzz_xyyz[k] + f2t * (tk0_zzz_yyz[k] - tk1_zzz_yyz[k]);

                    t_xzzz_xyzz[k] = fra * t10_zzz_xyzz[k] - frc * t11_zzz_xyzz[k] + f2t * (tk0_zzz_yzz[k] - tk1_zzz_yzz[k]);

                    t_xzzz_xzzz[k] = fra * t10_zzz_xzzz[k] - frc * t11_zzz_xzzz[k] + f2t * (tk0_zzz_zzz[k] - tk1_zzz_zzz[k]);

                    t_xzzz_yyyy[k] = fra * t10_zzz_yyyy[k] - frc * t11_zzz_yyyy[k];

                    t_xzzz_yyyz[k] = fra * t10_zzz_yyyz[k] - frc * t11_zzz_yyyz[k];

                    t_xzzz_yyzz[k] = fra * t10_zzz_yyzz[k] - frc * t11_zzz_yyzz[k];

                    t_xzzz_yzzz[k] = fra * t10_zzz_yzzz[k] - frc * t11_zzz_yzzz[k];

                    t_xzzz_zzzz[k] = fra * t10_zzz_zzzz[k] - frc * t11_zzz_zzzz[k];

                    // leading y component

                    fra = pay[k];

                    frc = pcy[k];

                    t_yyyy_xxxx[k] = fra * t10_yyy_xxxx[k] - frc * t11_yyy_xxxx[k] + f2t * (3.0 * t20_yy_xxxx[k] - 3.0 * t21_yy_xxxx[k]);

                    t_yyyy_xxxy[k] = fra * t10_yyy_xxxy[k] - frc * t11_yyy_xxxy[k] + f2t * (3.0 * t20_yy_xxxy[k] - 3.0 * t21_yy_xxxy[k] + tk0_yyy_xxx[k] - tk1_yyy_xxx[k]);

                    t_yyyy_xxxz[k] = fra * t10_yyy_xxxz[k] - frc * t11_yyy_xxxz[k] + f2t * (3.0 * t20_yy_xxxz[k] - 3.0 * t21_yy_xxxz[k]);

                    t_yyyy_xxyy[k] = fra * t10_yyy_xxyy[k] - frc * t11_yyy_xxyy[k] + f2t * (3.0 * t20_yy_xxyy[k] - 3.0 * t21_yy_xxyy[k] + 2.0 * tk0_yyy_xxy[k] - 2.0 * tk1_yyy_xxy[k]);

                    t_yyyy_xxyz[k] = fra * t10_yyy_xxyz[k] - frc * t11_yyy_xxyz[k] + f2t * (3.0 * t20_yy_xxyz[k] - 3.0 * t21_yy_xxyz[k] + tk0_yyy_xxz[k] - tk1_yyy_xxz[k]);

                    t_yyyy_xxzz[k] = fra * t10_yyy_xxzz[k] - frc * t11_yyy_xxzz[k] + f2t * (3.0 * t20_yy_xxzz[k] - 3.0 * t21_yy_xxzz[k]);

                    t_yyyy_xyyy[k] = fra * t10_yyy_xyyy[k] - frc * t11_yyy_xyyy[k] + f2t * (3.0 * t20_yy_xyyy[k] - 3.0 * t21_yy_xyyy[k] + 3.0 * tk0_yyy_xyy[k] - 3.0 * tk1_yyy_xyy[k]);

                    t_yyyy_xyyz[k] = fra * t10_yyy_xyyz[k] - frc * t11_yyy_xyyz[k] + f2t * (3.0 * t20_yy_xyyz[k] - 3.0 * t21_yy_xyyz[k] + 2.0 * tk0_yyy_xyz[k] - 2.0 * tk1_yyy_xyz[k]);

                    t_yyyy_xyzz[k] = fra * t10_yyy_xyzz[k] - frc * t11_yyy_xyzz[k] + f2t * (3.0 * t20_yy_xyzz[k] - 3.0 * t21_yy_xyzz[k] + tk0_yyy_xzz[k] - tk1_yyy_xzz[k]);

                    t_yyyy_xzzz[k] = fra * t10_yyy_xzzz[k] - frc * t11_yyy_xzzz[k] + f2t * (3.0 * t20_yy_xzzz[k] - 3.0 * t21_yy_xzzz[k]);

                    t_yyyy_yyyy[k] = fra * t10_yyy_yyyy[k] - frc * t11_yyy_yyyy[k] + f2t * (3.0 * t20_yy_yyyy[k] - 3.0 * t21_yy_yyyy[k] + 4.0 * tk0_yyy_yyy[k] - 4.0 * tk1_yyy_yyy[k]);

                    t_yyyy_yyyz[k] = fra * t10_yyy_yyyz[k] - frc * t11_yyy_yyyz[k] + f2t * (3.0 * t20_yy_yyyz[k] - 3.0 * t21_yy_yyyz[k] + 3.0 * tk0_yyy_yyz[k] - 3.0 * tk1_yyy_yyz[k]);

                    t_yyyy_yyzz[k] = fra * t10_yyy_yyzz[k] - frc * t11_yyy_yyzz[k] + f2t * (3.0 * t20_yy_yyzz[k] - 3.0 * t21_yy_yyzz[k] + 2.0 * tk0_yyy_yzz[k] - 2.0 * tk1_yyy_yzz[k]);

                    t_yyyy_yzzz[k] = fra * t10_yyy_yzzz[k] - frc * t11_yyy_yzzz[k] + f2t * (3.0 * t20_yy_yzzz[k] - 3.0 * t21_yy_yzzz[k] + tk0_yyy_zzz[k] - tk1_yyy_zzz[k]);

                    t_yyyy_zzzz[k] = fra * t10_yyy_zzzz[k] - frc * t11_yyy_zzzz[k] + f2t * (3.0 * t20_yy_zzzz[k] - 3.0 * t21_yy_zzzz[k]);

                    t_yyyz_xxxx[k] = fra * t10_yyz_xxxx[k] - frc * t11_yyz_xxxx[k] + f2t * (2.0 * t20_yz_xxxx[k] - 2.0 * t21_yz_xxxx[k]);

                    t_yyyz_xxxy[k] = fra * t10_yyz_xxxy[k] - frc * t11_yyz_xxxy[k] + f2t * (2.0 * t20_yz_xxxy[k] - 2.0 * t21_yz_xxxy[k] + tk0_yyz_xxx[k] - tk1_yyz_xxx[k]);

                    t_yyyz_xxxz[k] = fra * t10_yyz_xxxz[k] - frc * t11_yyz_xxxz[k] + f2t * (2.0 * t20_yz_xxxz[k] - 2.0 * t21_yz_xxxz[k]);

                    t_yyyz_xxyy[k] = fra * t10_yyz_xxyy[k] - frc * t11_yyz_xxyy[k] + f2t * (2.0 * t20_yz_xxyy[k] - 2.0 * t21_yz_xxyy[k] + 2.0 * tk0_yyz_xxy[k] - 2.0 * tk1_yyz_xxy[k]);

                    t_yyyz_xxyz[k] = fra * t10_yyz_xxyz[k] - frc * t11_yyz_xxyz[k] + f2t * (2.0 * t20_yz_xxyz[k] - 2.0 * t21_yz_xxyz[k] + tk0_yyz_xxz[k] - tk1_yyz_xxz[k]);

                    t_yyyz_xxzz[k] = fra * t10_yyz_xxzz[k] - frc * t11_yyz_xxzz[k] + f2t * (2.0 * t20_yz_xxzz[k] - 2.0 * t21_yz_xxzz[k]);

                    t_yyyz_xyyy[k] = fra * t10_yyz_xyyy[k] - frc * t11_yyz_xyyy[k] + f2t * (2.0 * t20_yz_xyyy[k] - 2.0 * t21_yz_xyyy[k] + 3.0 * tk0_yyz_xyy[k] - 3.0 * tk1_yyz_xyy[k]);

                    t_yyyz_xyyz[k] = fra * t10_yyz_xyyz[k] - frc * t11_yyz_xyyz[k] + f2t * (2.0 * t20_yz_xyyz[k] - 2.0 * t21_yz_xyyz[k] + 2.0 * tk0_yyz_xyz[k] - 2.0 * tk1_yyz_xyz[k]);

                    t_yyyz_xyzz[k] = fra * t10_yyz_xyzz[k] - frc * t11_yyz_xyzz[k] + f2t * (2.0 * t20_yz_xyzz[k] - 2.0 * t21_yz_xyzz[k] + tk0_yyz_xzz[k] - tk1_yyz_xzz[k]);

                    t_yyyz_xzzz[k] = fra * t10_yyz_xzzz[k] - frc * t11_yyz_xzzz[k] + f2t * (2.0 * t20_yz_xzzz[k] - 2.0 * t21_yz_xzzz[k]);

                    t_yyyz_yyyy[k] = fra * t10_yyz_yyyy[k] - frc * t11_yyz_yyyy[k] + f2t * (2.0 * t20_yz_yyyy[k] - 2.0 * t21_yz_yyyy[k] + 4.0 * tk0_yyz_yyy[k] - 4.0 * tk1_yyz_yyy[k]);

                    t_yyyz_yyyz[k] = fra * t10_yyz_yyyz[k] - frc * t11_yyz_yyyz[k] + f2t * (2.0 * t20_yz_yyyz[k] - 2.0 * t21_yz_yyyz[k] + 3.0 * tk0_yyz_yyz[k] - 3.0 * tk1_yyz_yyz[k]);

                    t_yyyz_yyzz[k] = fra * t10_yyz_yyzz[k] - frc * t11_yyz_yyzz[k] + f2t * (2.0 * t20_yz_yyzz[k] - 2.0 * t21_yz_yyzz[k] + 2.0 * tk0_yyz_yzz[k] - 2.0 * tk1_yyz_yzz[k]);

                    t_yyyz_yzzz[k] = fra * t10_yyz_yzzz[k] - frc * t11_yyz_yzzz[k] + f2t * (2.0 * t20_yz_yzzz[k] - 2.0 * t21_yz_yzzz[k] + tk0_yyz_zzz[k] - tk1_yyz_zzz[k]);

                    t_yyyz_zzzz[k] = fra * t10_yyz_zzzz[k] - frc * t11_yyz_zzzz[k] + f2t * (2.0 * t20_yz_zzzz[k] - 2.0 * t21_yz_zzzz[k]);

                    t_yyzz_xxxx[k] = fra * t10_yzz_xxxx[k] - frc * t11_yzz_xxxx[k] + f2t * (t20_zz_xxxx[k] - t21_zz_xxxx[k]);

                    t_yyzz_xxxy[k] = fra * t10_yzz_xxxy[k] - frc * t11_yzz_xxxy[k] + f2t * (t20_zz_xxxy[k] - t21_zz_xxxy[k] + tk0_yzz_xxx[k] - tk1_yzz_xxx[k]);

                    t_yyzz_xxxz[k] = fra * t10_yzz_xxxz[k] - frc * t11_yzz_xxxz[k] + f2t * (t20_zz_xxxz[k] - t21_zz_xxxz[k]);

                    t_yyzz_xxyy[k] = fra * t10_yzz_xxyy[k] - frc * t11_yzz_xxyy[k] + f2t * (t20_zz_xxyy[k] - t21_zz_xxyy[k] + 2.0 * tk0_yzz_xxy[k] - 2.0 * tk1_yzz_xxy[k]);

                    t_yyzz_xxyz[k] = fra * t10_yzz_xxyz[k] - frc * t11_yzz_xxyz[k] + f2t * (t20_zz_xxyz[k] - t21_zz_xxyz[k] + tk0_yzz_xxz[k] - tk1_yzz_xxz[k]);

                    t_yyzz_xxzz[k] = fra * t10_yzz_xxzz[k] - frc * t11_yzz_xxzz[k] + f2t * (t20_zz_xxzz[k] - t21_zz_xxzz[k]);

                    t_yyzz_xyyy[k] = fra * t10_yzz_xyyy[k] - frc * t11_yzz_xyyy[k] + f2t * (t20_zz_xyyy[k] - t21_zz_xyyy[k] + 3.0 * tk0_yzz_xyy[k] - 3.0 * tk1_yzz_xyy[k]);

                    t_yyzz_xyyz[k] = fra * t10_yzz_xyyz[k] - frc * t11_yzz_xyyz[k] + f2t * (t20_zz_xyyz[k] - t21_zz_xyyz[k] + 2.0 * tk0_yzz_xyz[k] - 2.0 * tk1_yzz_xyz[k]);

                    t_yyzz_xyzz[k] = fra * t10_yzz_xyzz[k] - frc * t11_yzz_xyzz[k] + f2t * (t20_zz_xyzz[k] - t21_zz_xyzz[k] + tk0_yzz_xzz[k] - tk1_yzz_xzz[k]);

                    t_yyzz_xzzz[k] = fra * t10_yzz_xzzz[k] - frc * t11_yzz_xzzz[k] + f2t * (t20_zz_xzzz[k] - t21_zz_xzzz[k]);

                    t_yyzz_yyyy[k] = fra * t10_yzz_yyyy[k] - frc * t11_yzz_yyyy[k] + f2t * (t20_zz_yyyy[k] - t21_zz_yyyy[k] + 4.0 * tk0_yzz_yyy[k] - 4.0 * tk1_yzz_yyy[k]);

                    t_yyzz_yyyz[k] = fra * t10_yzz_yyyz[k] - frc * t11_yzz_yyyz[k] + f2t * (t20_zz_yyyz[k] - t21_zz_yyyz[k] + 3.0 * tk0_yzz_yyz[k] - 3.0 * tk1_yzz_yyz[k]);

                    t_yyzz_yyzz[k] = fra * t10_yzz_yyzz[k] - frc * t11_yzz_yyzz[k] + f2t * (t20_zz_yyzz[k] - t21_zz_yyzz[k] + 2.0 * tk0_yzz_yzz[k] - 2.0 * tk1_yzz_yzz[k]);

                    t_yyzz_yzzz[k] = fra * t10_yzz_yzzz[k] - frc * t11_yzz_yzzz[k] + f2t * (t20_zz_yzzz[k] - t21_zz_yzzz[k] + tk0_yzz_zzz[k] - tk1_yzz_zzz[k]);

                    t_yyzz_zzzz[k] = fra * t10_yzz_zzzz[k] - frc * t11_yzz_zzzz[k] + f2t * (t20_zz_zzzz[k] - t21_zz_zzzz[k]);

                    t_yzzz_xxxx[k] = fra * t10_zzz_xxxx[k] - frc * t11_zzz_xxxx[k];

                    t_yzzz_xxxy[k] = fra * t10_zzz_xxxy[k] - frc * t11_zzz_xxxy[k] + f2t * (tk0_zzz_xxx[k] - tk1_zzz_xxx[k]);

                    t_yzzz_xxxz[k] = fra * t10_zzz_xxxz[k] - frc * t11_zzz_xxxz[k];

                    t_yzzz_xxyy[k] = fra * t10_zzz_xxyy[k] - frc * t11_zzz_xxyy[k] + f2t * (2.0 * tk0_zzz_xxy[k] - 2.0 * tk1_zzz_xxy[k]);

                    t_yzzz_xxyz[k] = fra * t10_zzz_xxyz[k] - frc * t11_zzz_xxyz[k] + f2t * (tk0_zzz_xxz[k] - tk1_zzz_xxz[k]);

                    t_yzzz_xxzz[k] = fra * t10_zzz_xxzz[k] - frc * t11_zzz_xxzz[k];

                    t_yzzz_xyyy[k] = fra * t10_zzz_xyyy[k] - frc * t11_zzz_xyyy[k] + f2t * (3.0 * tk0_zzz_xyy[k] - 3.0 * tk1_zzz_xyy[k]);

                    t_yzzz_xyyz[k] = fra * t10_zzz_xyyz[k] - frc * t11_zzz_xyyz[k] + f2t * (2.0 * tk0_zzz_xyz[k] - 2.0 * tk1_zzz_xyz[k]);

                    t_yzzz_xyzz[k] = fra * t10_zzz_xyzz[k] - frc * t11_zzz_xyzz[k] + f2t * (tk0_zzz_xzz[k] - tk1_zzz_xzz[k]);

                    t_yzzz_xzzz[k] = fra * t10_zzz_xzzz[k] - frc * t11_zzz_xzzz[k];

                    t_yzzz_yyyy[k] = fra * t10_zzz_yyyy[k] - frc * t11_zzz_yyyy[k] + f2t * (4.0 * tk0_zzz_yyy[k] - 4.0 * tk1_zzz_yyy[k]);

                    t_yzzz_yyyz[k] = fra * t10_zzz_yyyz[k] - frc * t11_zzz_yyyz[k] + f2t * (3.0 * tk0_zzz_yyz[k] - 3.0 * tk1_zzz_yyz[k]);

                    t_yzzz_yyzz[k] = fra * t10_zzz_yyzz[k] - frc * t11_zzz_yyzz[k] + f2t * (2.0 * tk0_zzz_yzz[k] - 2.0 * tk1_zzz_yzz[k]);

                    t_yzzz_yzzz[k] = fra * t10_zzz_yzzz[k] - frc * t11_zzz_yzzz[k] + f2t * (tk0_zzz_zzz[k] - tk1_zzz_zzz[k]);

                    t_yzzz_zzzz[k] = fra * t10_zzz_zzzz[k] - frc * t11_zzz_zzzz[k];

                    // leading z component

                    fra = paz[k];

                    frc = pcz[k];

                    t_zzzz_xxxx[k] = fra * t10_zzz_xxxx[k] - frc * t11_zzz_xxxx[k] + f2t * (3.0 * t20_zz_xxxx[k] - 3.0 * t21_zz_xxxx[k]);

                    t_zzzz_xxxy[k] = fra * t10_zzz_xxxy[k] - frc * t11_zzz_xxxy[k] + f2t * (3.0 * t20_zz_xxxy[k] - 3.0 * t21_zz_xxxy[k]);

                    t_zzzz_xxxz[k] = fra * t10_zzz_xxxz[k] - frc * t11_zzz_xxxz[k] + f2t * (3.0 * t20_zz_xxxz[k] - 3.0 * t21_zz_xxxz[k] + tk0_zzz_xxx[k] - tk1_zzz_xxx[k]);

                    t_zzzz_xxyy[k] = fra * t10_zzz_xxyy[k] - frc * t11_zzz_xxyy[k] + f2t * (3.0 * t20_zz_xxyy[k] - 3.0 * t21_zz_xxyy[k]);

                    t_zzzz_xxyz[k] = fra * t10_zzz_xxyz[k] - frc * t11_zzz_xxyz[k] + f2t * (3.0 * t20_zz_xxyz[k] - 3.0 * t21_zz_xxyz[k] + tk0_zzz_xxy[k] - tk1_zzz_xxy[k]);

                    t_zzzz_xxzz[k] = fra * t10_zzz_xxzz[k] - frc * t11_zzz_xxzz[k] + f2t * (3.0 * t20_zz_xxzz[k] - 3.0 * t21_zz_xxzz[k] + 2.0 * tk0_zzz_xxz[k] - 2.0 * tk1_zzz_xxz[k]);

                    t_zzzz_xyyy[k] = fra * t10_zzz_xyyy[k] - frc * t11_zzz_xyyy[k] + f2t * (3.0 * t20_zz_xyyy[k] - 3.0 * t21_zz_xyyy[k]);

                    t_zzzz_xyyz[k] = fra * t10_zzz_xyyz[k] - frc * t11_zzz_xyyz[k] + f2t * (3.0 * t20_zz_xyyz[k] - 3.0 * t21_zz_xyyz[k] + tk0_zzz_xyy[k] - tk1_zzz_xyy[k]);

                    t_zzzz_xyzz[k] = fra * t10_zzz_xyzz[k] - frc * t11_zzz_xyzz[k] + f2t * (3.0 * t20_zz_xyzz[k] - 3.0 * t21_zz_xyzz[k] + 2.0 * tk0_zzz_xyz[k] - 2.0 * tk1_zzz_xyz[k]);

                    t_zzzz_xzzz[k] = fra * t10_zzz_xzzz[k] - frc * t11_zzz_xzzz[k] + f2t * (3.0 * t20_zz_xzzz[k] - 3.0 * t21_zz_xzzz[k] + 3.0 * tk0_zzz_xzz[k] - 3.0 * tk1_zzz_xzz[k]);

                    t_zzzz_yyyy[k] = fra * t10_zzz_yyyy[k] - frc * t11_zzz_yyyy[k] + f2t * (3.0 * t20_zz_yyyy[k] - 3.0 * t21_zz_yyyy[k]);

                    t_zzzz_yyyz[k] = fra * t10_zzz_yyyz[k] - frc * t11_zzz_yyyz[k] + f2t * (3.0 * t20_zz_yyyz[k] - 3.0 * t21_zz_yyyz[k] + tk0_zzz_yyy[k] - tk1_zzz_yyy[k]);

                    t_zzzz_yyzz[k] = fra * t10_zzz_yyzz[k] - frc * t11_zzz_yyzz[k] + f2t * (3.0 * t20_zz_yyzz[k] - 3.0 * t21_zz_yyzz[k] + 2.0 * tk0_zzz_yyz[k] - 2.0 * tk1_zzz_yyz[k]);

                    t_zzzz_yzzz[k] = fra * t10_zzz_yzzz[k] - frc * t11_zzz_yzzz[k] + f2t * (3.0 * t20_zz_yzzz[k] - 3.0 * t21_zz_yzzz[k] + 3.0 * tk0_zzz_yzz[k] - 3.0 * tk1_zzz_yzz[k]);

                    t_zzzz_zzzz[k] = fra * t10_zzz_zzzz[k] - frc * t11_zzz_zzzz[k] + f2t * (3.0 * t20_zz_zzzz[k] - 3.0 * t21_zz_zzzz[k] + 4.0 * tk0_zzz_zzz[k] - 4.0 * tk1_zzz_zzz[k]);
                }

                idx++;
            }
        }
    }
    
} // npotrecfunc namespace
