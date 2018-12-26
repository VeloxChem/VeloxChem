//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFunc.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

namespace kinrecfunc { // kinrecfunc namespace
    
    void
    compKineticEnergyForSS(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  abDistances,
                           const CGtoBlock&            braGtoBlock,
                           const CGtoBlock&            ketGtoBlock,
                           const int32_t               iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // fetch up pi values
        
        auto fpi = mathconst::getPiValue();
        
        // get position of integrals in primitves buffer
        
        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});
        
        auto koff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(4 * idx);
            
            auto fz = osFactors.data(4 * idx + 1);
            
            // set up primitives buffer data
            
            auto s_0_0 = primBuffer.data(soff + idx);

            auto t_0_0 = primBuffer.data(koff + idx);
            
            #pragma omp simd aligned(s_0_0, t_0_0, fx, fz, abx, aby, abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double r2ab = abx[j] * abx[j] + aby[j] * aby[j] + abz[j] * abz[j];
                
                s_0_0[j] = std::pow(fpi * fx[j], 1.5) * std::exp(-fz[j] * r2ab);
                
                t_0_0[j] = fz[j] * (3.0 - 2.0 * fz[j] * r2ab) * s_0_0[j];
            }
            
            idx++;
        }
    }
    
    void
    compKineticEnergyForSP(      CMemBlock2D<double>&  primBuffer,
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

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(s1off + idx);

            // set up pointers to (S|T|S) integrals

            auto t_0_0 = primBuffer.data(k1off + idx);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(soff + 3 * idx);

            auto s_0_y = primBuffer.data(soff + 3 * idx + 1);

            auto s_0_z = primBuffer.data(soff + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(koff + 3 * idx);

            auto t_0_y = primBuffer.data(koff + 3 * idx + 1);

            auto t_0_z = primBuffer.data(koff + 3 * idx + 2);

            #pragma omp simd aligned(fz, pbx, pby, pbz, s_0_0, t_0_0, s_0_x,\
                                     s_0_y, s_0_z, t_0_x, t_0_y, t_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pbx[j];

                s_0_x[j] = fr * s_0_0[j];

                t_0_x[j] = fr * t_0_0[j] + f2z * s_0_x[j];

                // leading y component

                fr = pby[j];

                s_0_y[j] = fr * s_0_0[j];

                t_0_y[j] = fr * t_0_0[j] + f2z * s_0_y[j];

                // leading z component

                fr = pbz[j];

                s_0_z[j] = fr * s_0_0[j];

                t_0_z[j] = fr * t_0_0[j] + f2z * s_0_z[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForPS(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(s1off + idx);

            // set up pointers to (S|T|S) integrals

            auto t_0_0 = primBuffer.data(k1off + idx);

            // set up pointers to (P|S) integrals

            auto s_x_0 = primBuffer.data(soff + 3 * idx);

            auto s_y_0 = primBuffer.data(soff + 3 * idx + 1);

            auto s_z_0 = primBuffer.data(soff + 3 * idx + 2);

            // set up pointers to (P|T|S) integrals

            auto t_x_0 = primBuffer.data(koff + 3 * idx);

            auto t_y_0 = primBuffer.data(koff + 3 * idx + 1);

            auto t_z_0 = primBuffer.data(koff + 3 * idx + 2);

            #pragma omp simd aligned(fz, pax, pay, paz, s_0_0, t_0_0, s_x_0,\
                                     s_y_0, s_z_0, t_x_0, t_y_0, t_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pax[j];

                s_x_0[j] = fr * s_0_0[j];

                t_x_0[j] = fr * t_0_0[j] + f2z * s_x_0[j];

                // leading y component

                fr = pay[j];

                s_y_0[j] = fr * s_0_0[j];

                t_y_0[j] = fr * t_0_0[j] + f2z * s_y_0[j];

                // leading z component

                fr = paz[j];

                s_z_0[j] = fr * s_0_0[j];

                t_z_0[j] = fr * t_0_0[j] + f2z * s_z_0[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForPP(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(skoff + idx);

            // set up pointers to (S|T|S) integrals

            auto t_0_0 = primBuffer.data(kkoff + idx);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s1off + 3 * idx);

            auto s_0_y = primBuffer.data(s1off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s1off + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(k1off + 3 * idx);

            auto t_0_y = primBuffer.data(k1off + 3 * idx + 1);

            auto t_0_z = primBuffer.data(k1off + 3 * idx + 2);

            // set up pointers to (P|P) integrals

            auto s_x_x = primBuffer.data(soff + 9 * idx);

            auto s_x_y = primBuffer.data(soff + 9 * idx + 1);

            auto s_x_z = primBuffer.data(soff + 9 * idx + 2);

            auto s_y_x = primBuffer.data(soff + 9 * idx + 3);

            auto s_y_y = primBuffer.data(soff + 9 * idx + 4);

            auto s_y_z = primBuffer.data(soff + 9 * idx + 5);

            auto s_z_x = primBuffer.data(soff + 9 * idx + 6);

            auto s_z_y = primBuffer.data(soff + 9 * idx + 7);

            auto s_z_z = primBuffer.data(soff + 9 * idx + 8);

            // set up pointers to (P|T|P) integrals

            auto t_x_x = primBuffer.data(koff + 9 * idx);

            auto t_x_y = primBuffer.data(koff + 9 * idx + 1);

            auto t_x_z = primBuffer.data(koff + 9 * idx + 2);

            auto t_y_x = primBuffer.data(koff + 9 * idx + 3);

            auto t_y_y = primBuffer.data(koff + 9 * idx + 4);

            auto t_y_z = primBuffer.data(koff + 9 * idx + 5);

            auto t_z_x = primBuffer.data(koff + 9 * idx + 6);

            auto t_z_y = primBuffer.data(koff + 9 * idx + 7);

            auto t_z_z = primBuffer.data(koff + 9 * idx + 8);

            #pragma omp simd aligned(fx, fz, pax, pay, paz, s_0_0, t_0_0, s_0_x,\
                                     s_0_y, s_0_z, t_0_x, t_0_y, t_0_z, s_x_x,\
                                     s_x_y, s_x_z, s_y_x, s_y_y, s_y_z, s_z_x,\
                                     s_z_y, s_z_z, t_x_x, t_x_y, t_x_z, t_y_x,\
                                     t_y_y, t_y_z, t_z_x, t_z_y, t_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pax[j];

                s_x_x[j] = fr * s_0_x[j] + f2t * s_0_0[j];

                t_x_x[j] = fr * t_0_x[j] + f2t * t_0_0[j] + f2z * s_x_x[j];

                s_x_y[j] = fr * s_0_y[j];

                t_x_y[j] = fr * t_0_y[j] + f2z * s_x_y[j];

                s_x_z[j] = fr * s_0_z[j];

                t_x_z[j] = fr * t_0_z[j] + f2z * s_x_z[j];

                // leading y component

                fr = pay[j];

                s_y_x[j] = fr * s_0_x[j];

                t_y_x[j] = fr * t_0_x[j] + f2z * s_y_x[j];

                s_y_y[j] = fr * s_0_y[j] + f2t * s_0_0[j];

                t_y_y[j] = fr * t_0_y[j] + f2t * t_0_0[j] + f2z * s_y_y[j];

                s_y_z[j] = fr * s_0_z[j];

                t_y_z[j] = fr * t_0_z[j] + f2z * s_y_z[j];

                // leading z component

                fr = paz[j];

                s_z_x[j] = fr * s_0_x[j];

                t_z_x[j] = fr * t_0_x[j] + f2z * s_z_x[j];

                s_z_y[j] = fr * s_0_y[j];

                t_z_y[j] = fr * t_0_y[j] + f2z * s_z_y[j];

                s_z_z[j] = fr * s_0_z[j] + f2t * s_0_0[j];

                t_z_z[j] = fr * t_0_z[j] + f2t * t_0_0[j] + f2z * s_z_z[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForSD(      CMemBlock2D<double>&  primBuffer,
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

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(s2off + idx);

            // set up pointers to (S|T|S) integrals

            auto t_0_0 = primBuffer.data(k2off + idx);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s1off + 3 * idx);

            auto s_0_y = primBuffer.data(s1off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s1off + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(k1off + 3 * idx);

            auto t_0_y = primBuffer.data(k1off + 3 * idx + 1);

            auto t_0_z = primBuffer.data(k1off + 3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(soff + 6 * idx);

            auto s_0_xy = primBuffer.data(soff + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(soff + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(soff + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(soff + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(soff + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(koff + 6 * idx);

            auto t_0_xy = primBuffer.data(koff + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(koff + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(koff + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(koff + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(koff + 6 * idx + 5);

            #pragma omp simd aligned(fx, fz, fgb, pbx, pby, pbz, s_0_0, t_0_0,\
                                     s_0_x, s_0_y, s_0_z, t_0_x, t_0_y, t_0_z,\
                                     s_0_xx, s_0_xy, s_0_xz, s_0_yy, s_0_yz, s_0_zz,\
                                     t_0_xx, t_0_xy, t_0_xz, t_0_yy, t_0_yz, t_0_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2b = 0.50 * fgb[j];

                // leading x component

                double fr = pbx[j];

                s_0_xx[j] = fr * s_0_x[j] + f2t * s_0_0[j];

                t_0_xx[j] = fr * t_0_x[j] + f2t * t_0_0[j] + f2z * (s_0_xx[j] - f2b * s_0_0[j]);

                s_0_xy[j] = fr * s_0_y[j];

                t_0_xy[j] = fr * t_0_y[j] + f2z * s_0_xy[j];

                s_0_xz[j] = fr * s_0_z[j];

                t_0_xz[j] = fr * t_0_z[j] + f2z * s_0_xz[j];

                // leading y component

                fr = pby[j];

                s_0_yy[j] = fr * s_0_y[j] + f2t * s_0_0[j];

                t_0_yy[j] = fr * t_0_y[j] + f2t * t_0_0[j] + f2z * (s_0_yy[j] - f2b * s_0_0[j]);

                s_0_yz[j] = fr * s_0_z[j];

                t_0_yz[j] = fr * t_0_z[j] + f2z * s_0_yz[j];

                // leading z component

                fr = pbz[j];

                s_0_zz[j] = fr * s_0_z[j] + f2t * s_0_0[j];

                t_0_zz[j] = fr * t_0_z[j] + f2t * t_0_0[j] + f2z * (s_0_zz[j] - f2b * s_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDS(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|S) integrals

            auto s_0_0 = primBuffer.data(s2off + idx);

            // set up pointers to (S|T|S) integrals

            auto t_0_0 = primBuffer.data(k2off + idx);

            // set up pointers to (P|S) integrals

            auto s_x_0 = primBuffer.data(s1off + 3 * idx);

            auto s_y_0 = primBuffer.data(s1off + 3 * idx + 1);

            auto s_z_0 = primBuffer.data(s1off + 3 * idx + 2);

            // set up pointers to (P|T|S) integrals

            auto t_x_0 = primBuffer.data(k1off + 3 * idx);

            auto t_y_0 = primBuffer.data(k1off + 3 * idx + 1);

            auto t_z_0 = primBuffer.data(k1off + 3 * idx + 2);

            // set up pointers to (D|S) integrals

            auto s_xx_0 = primBuffer.data(soff + 6 * idx);

            auto s_xy_0 = primBuffer.data(soff + 6 * idx + 1);

            auto s_xz_0 = primBuffer.data(soff + 6 * idx + 2);

            auto s_yy_0 = primBuffer.data(soff + 6 * idx + 3);

            auto s_yz_0 = primBuffer.data(soff + 6 * idx + 4);

            auto s_zz_0 = primBuffer.data(soff + 6 * idx + 5);

            // set up pointers to (D|T|S) integrals

            auto t_xx_0 = primBuffer.data(koff + 6 * idx);

            auto t_xy_0 = primBuffer.data(koff + 6 * idx + 1);

            auto t_xz_0 = primBuffer.data(koff + 6 * idx + 2);

            auto t_yy_0 = primBuffer.data(koff + 6 * idx + 3);

            auto t_yz_0 = primBuffer.data(koff + 6 * idx + 4);

            auto t_zz_0 = primBuffer.data(koff + 6 * idx + 5);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_0_0, t_0_0,\
                                     s_x_0, s_y_0, s_z_0, t_x_0, t_y_0, t_z_0,\
                                     s_xx_0, s_xy_0, s_xz_0, s_yy_0, s_yz_0, s_zz_0,\
                                     t_xx_0, t_xy_0, t_xz_0, t_yy_0, t_yz_0, t_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xx_0[j] = fr * s_x_0[j] + f2t * s_0_0[j];

                t_xx_0[j] = fr * t_x_0[j] + f2t * t_0_0[j] + f2z * (s_xx_0[j] - f2a * s_0_0[j]);

                s_xy_0[j] = fr * s_y_0[j];

                t_xy_0[j] = fr * t_y_0[j] + f2z * s_xy_0[j];

                s_xz_0[j] = fr * s_z_0[j];

                t_xz_0[j] = fr * t_z_0[j] + f2z * s_xz_0[j];

                // leading y component

                fr = pay[j];

                s_yy_0[j] = fr * s_y_0[j] + f2t * s_0_0[j];

                t_yy_0[j] = fr * t_y_0[j] + f2t * t_0_0[j] + f2z * (s_yy_0[j] - f2a * s_0_0[j]);

                s_yz_0[j] = fr * s_z_0[j];

                t_yz_0[j] = fr * t_z_0[j] + f2z * s_yz_0[j];

                // leading z component

                fr = paz[j];

                s_zz_0[j] = fr * s_z_0[j] + f2t * s_0_0[j];

                t_zz_0[j] = fr * t_z_0[j] + f2t * t_0_0[j] + f2z * (s_zz_0[j] - f2a * s_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForPD(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(skoff + 3 * idx);

            auto s_0_y = primBuffer.data(skoff + 3 * idx + 1);

            auto s_0_z = primBuffer.data(skoff + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(kkoff + 3 * idx);

            auto t_0_y = primBuffer.data(kkoff + 3 * idx + 1);

            auto t_0_z = primBuffer.data(kkoff + 3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(s1off + 6 * idx);

            auto s_0_xy = primBuffer.data(s1off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(s1off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(s1off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(s1off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(s1off + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(k1off + 6 * idx);

            auto t_0_xy = primBuffer.data(k1off + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(k1off + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(k1off + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(k1off + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(k1off + 6 * idx + 5);

            // set up pointers to (P|D) integrals

            auto s_x_xx = primBuffer.data(soff + 18 * idx);

            auto s_x_xy = primBuffer.data(soff + 18 * idx + 1);

            auto s_x_xz = primBuffer.data(soff + 18 * idx + 2);

            auto s_x_yy = primBuffer.data(soff + 18 * idx + 3);

            auto s_x_yz = primBuffer.data(soff + 18 * idx + 4);

            auto s_x_zz = primBuffer.data(soff + 18 * idx + 5);

            auto s_y_xx = primBuffer.data(soff + 18 * idx + 6);

            auto s_y_xy = primBuffer.data(soff + 18 * idx + 7);

            auto s_y_xz = primBuffer.data(soff + 18 * idx + 8);

            auto s_y_yy = primBuffer.data(soff + 18 * idx + 9);

            auto s_y_yz = primBuffer.data(soff + 18 * idx + 10);

            auto s_y_zz = primBuffer.data(soff + 18 * idx + 11);

            auto s_z_xx = primBuffer.data(soff + 18 * idx + 12);

            auto s_z_xy = primBuffer.data(soff + 18 * idx + 13);

            auto s_z_xz = primBuffer.data(soff + 18 * idx + 14);

            auto s_z_yy = primBuffer.data(soff + 18 * idx + 15);

            auto s_z_yz = primBuffer.data(soff + 18 * idx + 16);

            auto s_z_zz = primBuffer.data(soff + 18 * idx + 17);

            // set up pointers to (P|T|D) integrals

            auto t_x_xx = primBuffer.data(koff + 18 * idx);

            auto t_x_xy = primBuffer.data(koff + 18 * idx + 1);

            auto t_x_xz = primBuffer.data(koff + 18 * idx + 2);

            auto t_x_yy = primBuffer.data(koff + 18 * idx + 3);

            auto t_x_yz = primBuffer.data(koff + 18 * idx + 4);

            auto t_x_zz = primBuffer.data(koff + 18 * idx + 5);

            auto t_y_xx = primBuffer.data(koff + 18 * idx + 6);

            auto t_y_xy = primBuffer.data(koff + 18 * idx + 7);

            auto t_y_xz = primBuffer.data(koff + 18 * idx + 8);

            auto t_y_yy = primBuffer.data(koff + 18 * idx + 9);

            auto t_y_yz = primBuffer.data(koff + 18 * idx + 10);

            auto t_y_zz = primBuffer.data(koff + 18 * idx + 11);

            auto t_z_xx = primBuffer.data(koff + 18 * idx + 12);

            auto t_z_xy = primBuffer.data(koff + 18 * idx + 13);

            auto t_z_xz = primBuffer.data(koff + 18 * idx + 14);

            auto t_z_yy = primBuffer.data(koff + 18 * idx + 15);

            auto t_z_yz = primBuffer.data(koff + 18 * idx + 16);

            auto t_z_zz = primBuffer.data(koff + 18 * idx + 17);

            #pragma omp simd aligned(fx, fz, pax, pay, paz, s_0_x, s_0_y, s_0_z,\
                                     t_0_x, t_0_y, t_0_z, s_0_xx, s_0_xy, s_0_xz,\
                                     s_0_yy, s_0_yz, s_0_zz, t_0_xx, t_0_xy, t_0_xz,\
                                     t_0_yy, t_0_yz, t_0_zz, s_x_xx, s_x_xy, s_x_xz,\
                                     s_x_yy, s_x_yz, s_x_zz, s_y_xx, s_y_xy, s_y_xz,\
                                     s_y_yy, s_y_yz, s_y_zz, s_z_xx, s_z_xy, s_z_xz,\
                                     s_z_yy, s_z_yz, s_z_zz, t_x_xx, t_x_xy, t_x_xz,\
                                     t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy, t_y_xz,\
                                     t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy, t_z_xz,\
                                     t_z_yy, t_z_yz, t_z_zz: VLX_ALIGN)
             for (int32_t j = 0; j < nprim; j++)
             {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pax[j];

                s_x_xx[j] = fr * s_0_xx[j] + f2t * 2.0 * s_0_x[j];

                t_x_xx[j] = fr * t_0_xx[j] + f2t * 2.0 * t_0_x[j] + f2z * s_x_xx[j];

                s_x_xy[j] = fr * s_0_xy[j] + f2t * s_0_y[j];

                t_x_xy[j] = fr * t_0_xy[j] + f2t * t_0_y[j] + f2z * s_x_xy[j];

                s_x_xz[j] = fr * s_0_xz[j] + f2t * s_0_z[j];

                t_x_xz[j] = fr * t_0_xz[j] + f2t * t_0_z[j] + f2z * s_x_xz[j];

                s_x_yy[j] = fr * s_0_yy[j];

                t_x_yy[j] = fr * t_0_yy[j] + f2z * s_x_yy[j];

                s_x_yz[j] = fr * s_0_yz[j];

                t_x_yz[j] = fr * t_0_yz[j] + f2z * s_x_yz[j];

                s_x_zz[j] = fr * s_0_zz[j];

                t_x_zz[j] = fr * t_0_zz[j] + f2z * s_x_zz[j];

                // leading y component

                fr = pay[j];

                s_y_xx[j] = fr * s_0_xx[j];

                t_y_xx[j] = fr * t_0_xx[j] + f2z * s_y_xx[j];

                s_y_xy[j] = fr * s_0_xy[j] + f2t * s_0_x[j];

                t_y_xy[j] = fr * t_0_xy[j] + f2t * t_0_x[j] + f2z * s_y_xy[j];

                s_y_xz[j] = fr * s_0_xz[j];

                t_y_xz[j] = fr * t_0_xz[j] + f2z * s_y_xz[j];

                s_y_yy[j] = fr * s_0_yy[j] + f2t * 2.0 * s_0_y[j];

                t_y_yy[j] = fr * t_0_yy[j] + f2t * 2.0 * t_0_y[j] + f2z * s_y_yy[j];

                s_y_yz[j] = fr * s_0_yz[j] + f2t * s_0_z[j];

                t_y_yz[j] = fr * t_0_yz[j] + f2t * t_0_z[j] + f2z * s_y_yz[j];

                s_y_zz[j] = fr * s_0_zz[j];

                t_y_zz[j] = fr * t_0_zz[j] + f2z * s_y_zz[j];

                // leading z component

                fr = paz[j];

                s_z_xx[j] = fr * s_0_xx[j];

                t_z_xx[j] = fr * t_0_xx[j] + f2z * s_z_xx[j];

                s_z_xy[j] = fr * s_0_xy[j];

                t_z_xy[j] = fr * t_0_xy[j] + f2z * s_z_xy[j];

                s_z_xz[j] = fr * s_0_xz[j] + f2t * s_0_x[j];

                t_z_xz[j] = fr * t_0_xz[j] + f2t * t_0_x[j] + f2z * s_z_xz[j];

                s_z_yy[j] = fr * s_0_yy[j];

                t_z_yy[j] = fr * t_0_yy[j] + f2z * s_z_yy[j];

                s_z_yz[j] = fr * s_0_yz[j] + f2t * s_0_y[j];

                t_z_yz[j] = fr * t_0_yz[j] + f2t * t_0_y[j] + f2z * s_z_yz[j];

                s_z_zz[j] = fr * s_0_zz[j] + f2t * 2.0 * s_0_z[j];

                t_z_zz[j] = fr * t_0_zz[j] + f2t * 2.0 * t_0_z[j] + f2z * s_z_zz[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDP(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|S) integrals

            auto s_x_0 = primBuffer.data(skoff + 3 * idx);

            auto s_y_0 = primBuffer.data(skoff + 3 * idx + 1);

            auto s_z_0 = primBuffer.data(skoff + 3 * idx + 2);

            // set up pointers to (P|T|S) integrals

            auto t_x_0 = primBuffer.data(kkoff + 3 * idx);

            auto t_y_0 = primBuffer.data(kkoff + 3 * idx + 1);

            auto t_z_0 = primBuffer.data(kkoff + 3 * idx + 2);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s2off + 3 * idx);

            auto s_0_y = primBuffer.data(s2off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s2off + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(k2off + 3 * idx);

            auto t_0_y = primBuffer.data(k2off + 3 * idx + 1);

            auto t_0_z = primBuffer.data(k2off + 3 * idx + 2);

            // set up pointers to (P|P) integrals

            auto s_x_x = primBuffer.data(s1off + 9 * idx);

            auto s_x_y = primBuffer.data(s1off + 9 * idx + 1);

            auto s_x_z = primBuffer.data(s1off + 9 * idx + 2);

            auto s_y_x = primBuffer.data(s1off + 9 * idx + 3);

            auto s_y_y = primBuffer.data(s1off + 9 * idx + 4);

            auto s_y_z = primBuffer.data(s1off + 9 * idx + 5);

            auto s_z_x = primBuffer.data(s1off + 9 * idx + 6);

            auto s_z_y = primBuffer.data(s1off + 9 * idx + 7);

            auto s_z_z = primBuffer.data(s1off + 9 * idx + 8);

            // set up pointers to (P|T|P) integrals

            auto t_x_x = primBuffer.data(k1off + 9 * idx);

            auto t_x_y = primBuffer.data(k1off + 9 * idx + 1);

            auto t_x_z = primBuffer.data(k1off + 9 * idx + 2);

            auto t_y_x = primBuffer.data(k1off + 9 * idx + 3);

            auto t_y_y = primBuffer.data(k1off + 9 * idx + 4);

            auto t_y_z = primBuffer.data(k1off + 9 * idx + 5);

            auto t_z_x = primBuffer.data(k1off + 9 * idx + 6);

            auto t_z_y = primBuffer.data(k1off + 9 * idx + 7);

            auto t_z_z = primBuffer.data(k1off + 9 * idx + 8);

            // set up pointers to (D|P) integrals

            auto s_xx_x = primBuffer.data(soff + 18 * idx);

            auto s_xx_y = primBuffer.data(soff + 18 * idx + 1);

            auto s_xx_z = primBuffer.data(soff + 18 * idx + 2);

            auto s_xy_x = primBuffer.data(soff + 18 * idx + 3);

            auto s_xy_y = primBuffer.data(soff + 18 * idx + 4);

            auto s_xy_z = primBuffer.data(soff + 18 * idx + 5);

            auto s_xz_x = primBuffer.data(soff + 18 * idx + 6);

            auto s_xz_y = primBuffer.data(soff + 18 * idx + 7);

            auto s_xz_z = primBuffer.data(soff + 18 * idx + 8);

            auto s_yy_x = primBuffer.data(soff + 18 * idx + 9);

            auto s_yy_y = primBuffer.data(soff + 18 * idx + 10);

            auto s_yy_z = primBuffer.data(soff + 18 * idx + 11);

            auto s_yz_x = primBuffer.data(soff + 18 * idx + 12);

            auto s_yz_y = primBuffer.data(soff + 18 * idx + 13);

            auto s_yz_z = primBuffer.data(soff + 18 * idx + 14);

            auto s_zz_x = primBuffer.data(soff + 18 * idx + 15);

            auto s_zz_y = primBuffer.data(soff + 18 * idx + 16);

            auto s_zz_z = primBuffer.data(soff + 18 * idx + 17);

            // set up pointers to (D|T|P) integrals

            auto t_xx_x = primBuffer.data(koff + 18 * idx);

            auto t_xx_y = primBuffer.data(koff + 18 * idx + 1);

            auto t_xx_z = primBuffer.data(koff + 18 * idx + 2);

            auto t_xy_x = primBuffer.data(koff + 18 * idx + 3);

            auto t_xy_y = primBuffer.data(koff + 18 * idx + 4);

            auto t_xy_z = primBuffer.data(koff + 18 * idx + 5);

            auto t_xz_x = primBuffer.data(koff + 18 * idx + 6);

            auto t_xz_y = primBuffer.data(koff + 18 * idx + 7);

            auto t_xz_z = primBuffer.data(koff + 18 * idx + 8);

            auto t_yy_x = primBuffer.data(koff + 18 * idx + 9);

            auto t_yy_y = primBuffer.data(koff + 18 * idx + 10);

            auto t_yy_z = primBuffer.data(koff + 18 * idx + 11);

            auto t_yz_x = primBuffer.data(koff + 18 * idx + 12);

            auto t_yz_y = primBuffer.data(koff + 18 * idx + 13);

            auto t_yz_z = primBuffer.data(koff + 18 * idx + 14);

            auto t_zz_x = primBuffer.data(koff + 18 * idx + 15);

            auto t_zz_y = primBuffer.data(koff + 18 * idx + 16);

            auto t_zz_z = primBuffer.data(koff + 18 * idx + 17);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_x_0, s_y_0,\
                                     s_z_0, t_x_0, t_y_0, t_z_0, s_0_x, s_0_y,\
                                     s_0_z, t_0_x, t_0_y, t_0_z, s_x_x, s_x_y,\
                                     s_x_z, s_y_x, s_y_y, s_y_z, s_z_x, s_z_y,\
                                     s_z_z, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y,\
                                     t_y_z, t_z_x, t_z_y, t_z_z, s_xx_x, s_xx_y,\
                                     s_xx_z, s_xy_x, s_xy_y, s_xy_z, s_xz_x, s_xz_y,\
                                     s_xz_z, s_yy_x, s_yy_y, s_yy_z, s_yz_x, s_yz_y,\
                                     s_yz_z, s_zz_x, s_zz_y, s_zz_z, t_xx_x, t_xx_y,\
                                     t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y,\
                                     t_xz_z, t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y,\
                                     t_yz_z, t_zz_x, t_zz_y, t_zz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xx_x[j] = fr * s_x_x[j] + f2t * (s_0_x[j] + s_x_0[j]);

                t_xx_x[j] = fr * t_x_x[j] + f2t * (t_0_x[j] + t_x_0[j]) + f2z * (s_xx_x[j] - f2a * s_0_x[j]);

                s_xx_y[j] = fr * s_x_y[j] + f2t * s_0_y[j];

                t_xx_y[j] = fr * t_x_y[j] + f2t * t_0_y[j] + f2z * (s_xx_y[j] - f2a * s_0_y[j]);

                s_xx_z[j] = fr * s_x_z[j] + f2t * s_0_z[j];

                t_xx_z[j] = fr * t_x_z[j] + f2t * t_0_z[j] + f2z * (s_xx_z[j] - f2a * s_0_z[j]);

                s_xy_x[j] = fr * s_y_x[j] + f2t * s_y_0[j];

                t_xy_x[j] = fr * t_y_x[j] + f2t * t_y_0[j] + f2z * s_xy_x[j];

                s_xy_y[j] = fr * s_y_y[j];

                t_xy_y[j] = fr * t_y_y[j] + f2z * s_xy_y[j];

                s_xy_z[j] = fr * s_y_z[j];

                t_xy_z[j] = fr * t_y_z[j] + f2z * s_xy_z[j];

                s_xz_x[j] = fr * s_z_x[j] + f2t * s_z_0[j];

                t_xz_x[j] = fr * t_z_x[j] + f2t * t_z_0[j] + f2z * s_xz_x[j];

                s_xz_y[j] = fr * s_z_y[j];

                t_xz_y[j] = fr * t_z_y[j] + f2z * s_xz_y[j];

                s_xz_z[j] = fr * s_z_z[j];

                t_xz_z[j] = fr * t_z_z[j] + f2z * s_xz_z[j];

                // leading y component

                fr = pay[j];

                s_yy_x[j] = fr * s_y_x[j] + f2t * s_0_x[j];

                t_yy_x[j] = fr * t_y_x[j] + f2t * t_0_x[j] + f2z * (s_yy_x[j] - f2a * s_0_x[j]);

                s_yy_y[j] = fr * s_y_y[j] + f2t * (s_0_y[j] + s_y_0[j]);

                t_yy_y[j] = fr * t_y_y[j] + f2t * (t_0_y[j] + t_y_0[j]) + f2z * (s_yy_y[j] - f2a * s_0_y[j]);

                s_yy_z[j] = fr * s_y_z[j] + f2t * s_0_z[j];

                t_yy_z[j] = fr * t_y_z[j] + f2t * t_0_z[j] + f2z * (s_yy_z[j] - f2a * s_0_z[j]);

                s_yz_x[j] = fr * s_z_x[j];

                t_yz_x[j] = fr * t_z_x[j] + f2z * s_yz_x[j];

                s_yz_y[j] = fr * s_z_y[j] + f2t * s_z_0[j];

                t_yz_y[j] = fr * t_z_y[j] + f2t * t_z_0[j] + f2z * s_yz_y[j];

                s_yz_z[j] = fr * s_z_z[j];

                t_yz_z[j] = fr * t_z_z[j] + f2z * s_yz_z[j];

                // leading z component

                fr = paz[j];

                s_zz_x[j] = fr * s_z_x[j] + f2t * s_0_x[j];

                t_zz_x[j] = fr * t_z_x[j] + f2t * t_0_x[j] + f2z * (s_zz_x[j] - f2a * s_0_x[j]);

                s_zz_y[j] = fr * s_z_y[j] + f2t * s_0_y[j];

                t_zz_y[j] = fr * t_z_y[j] + f2t * t_0_y[j] + f2z * (s_zz_y[j] - f2a * s_0_y[j]);

                s_zz_z[j] = fr * s_z_z[j] + f2t * (s_0_z[j] + s_z_0[j]);

                t_zz_z[j] = fr * t_z_z[j] + f2t * (t_0_z[j] + t_z_0[j]) + f2z * (s_zz_z[j] - f2a * s_0_z[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDD(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|P) integrals

            auto s_x_x = primBuffer.data(skoff + 9 * idx);

            auto s_x_y = primBuffer.data(skoff + 9 * idx + 1);

            auto s_x_z = primBuffer.data(skoff + 9 * idx + 2);

            auto s_y_x = primBuffer.data(skoff + 9 * idx + 3);

            auto s_y_y = primBuffer.data(skoff + 9 * idx + 4);

            auto s_y_z = primBuffer.data(skoff + 9 * idx + 5);

            auto s_z_x = primBuffer.data(skoff + 9 * idx + 6);

            auto s_z_y = primBuffer.data(skoff + 9 * idx + 7);

            auto s_z_z = primBuffer.data(skoff + 9 * idx + 8);

            // set up pointers to (P|T|P) integrals

            auto t_x_x = primBuffer.data(kkoff + 9 * idx);

            auto t_x_y = primBuffer.data(kkoff + 9 * idx + 1);

            auto t_x_z = primBuffer.data(kkoff + 9 * idx + 2);

            auto t_y_x = primBuffer.data(kkoff + 9 * idx + 3);

            auto t_y_y = primBuffer.data(kkoff + 9 * idx + 4);

            auto t_y_z = primBuffer.data(kkoff + 9 * idx + 5);

            auto t_z_x = primBuffer.data(kkoff + 9 * idx + 6);

            auto t_z_y = primBuffer.data(kkoff + 9 * idx + 7);

            auto t_z_z = primBuffer.data(kkoff + 9 * idx + 8);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(s2off + 6 * idx);

            auto s_0_xy = primBuffer.data(s2off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(s2off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(s2off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(s2off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(s2off + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(k2off + 6 * idx);

            auto t_0_xy = primBuffer.data(k2off + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(k2off + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(k2off + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(k2off + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(k2off + 6 * idx + 5);

            // set up pointers to (P|D) integrals

            auto s_x_xx = primBuffer.data(s1off + 18 * idx);

            auto s_x_xy = primBuffer.data(s1off + 18 * idx + 1);

            auto s_x_xz = primBuffer.data(s1off + 18 * idx + 2);

            auto s_x_yy = primBuffer.data(s1off + 18 * idx + 3);

            auto s_x_yz = primBuffer.data(s1off + 18 * idx + 4);

            auto s_x_zz = primBuffer.data(s1off + 18 * idx + 5);

            auto s_y_xx = primBuffer.data(s1off + 18 * idx + 6);

            auto s_y_xy = primBuffer.data(s1off + 18 * idx + 7);

            auto s_y_xz = primBuffer.data(s1off + 18 * idx + 8);

            auto s_y_yy = primBuffer.data(s1off + 18 * idx + 9);

            auto s_y_yz = primBuffer.data(s1off + 18 * idx + 10);

            auto s_y_zz = primBuffer.data(s1off + 18 * idx + 11);

            auto s_z_xx = primBuffer.data(s1off + 18 * idx + 12);

            auto s_z_xy = primBuffer.data(s1off + 18 * idx + 13);

            auto s_z_xz = primBuffer.data(s1off + 18 * idx + 14);

            auto s_z_yy = primBuffer.data(s1off + 18 * idx + 15);

            auto s_z_yz = primBuffer.data(s1off + 18 * idx + 16);

            auto s_z_zz = primBuffer.data(s1off + 18 * idx + 17);

            // set up pointers to (P|T|D) integrals

            auto t_x_xx = primBuffer.data(k1off + 18 * idx);

            auto t_x_xy = primBuffer.data(k1off + 18 * idx + 1);

            auto t_x_xz = primBuffer.data(k1off + 18 * idx + 2);

            auto t_x_yy = primBuffer.data(k1off + 18 * idx + 3);

            auto t_x_yz = primBuffer.data(k1off + 18 * idx + 4);

            auto t_x_zz = primBuffer.data(k1off + 18 * idx + 5);

            auto t_y_xx = primBuffer.data(k1off + 18 * idx + 6);

            auto t_y_xy = primBuffer.data(k1off + 18 * idx + 7);

            auto t_y_xz = primBuffer.data(k1off + 18 * idx + 8);

            auto t_y_yy = primBuffer.data(k1off + 18 * idx + 9);

            auto t_y_yz = primBuffer.data(k1off + 18 * idx + 10);

            auto t_y_zz = primBuffer.data(k1off + 18 * idx + 11);

            auto t_z_xx = primBuffer.data(k1off + 18 * idx + 12);

            auto t_z_xy = primBuffer.data(k1off + 18 * idx + 13);

            auto t_z_xz = primBuffer.data(k1off + 18 * idx + 14);

            auto t_z_yy = primBuffer.data(k1off + 18 * idx + 15);

            auto t_z_yz = primBuffer.data(k1off + 18 * idx + 16);

            auto t_z_zz = primBuffer.data(k1off + 18 * idx + 17);

            // set up pointers to (D|D) integrals

            auto s_xx_xx = primBuffer.data(soff + 36 * idx);

            auto s_xx_xy = primBuffer.data(soff + 36 * idx + 1);

            auto s_xx_xz = primBuffer.data(soff + 36 * idx + 2);

            auto s_xx_yy = primBuffer.data(soff + 36 * idx + 3);

            auto s_xx_yz = primBuffer.data(soff + 36 * idx + 4);

            auto s_xx_zz = primBuffer.data(soff + 36 * idx + 5);

            auto s_xy_xx = primBuffer.data(soff + 36 * idx + 6);

            auto s_xy_xy = primBuffer.data(soff + 36 * idx + 7);

            auto s_xy_xz = primBuffer.data(soff + 36 * idx + 8);

            auto s_xy_yy = primBuffer.data(soff + 36 * idx + 9);

            auto s_xy_yz = primBuffer.data(soff + 36 * idx + 10);

            auto s_xy_zz = primBuffer.data(soff + 36 * idx + 11);

            auto s_xz_xx = primBuffer.data(soff + 36 * idx + 12);

            auto s_xz_xy = primBuffer.data(soff + 36 * idx + 13);

            auto s_xz_xz = primBuffer.data(soff + 36 * idx + 14);

            auto s_xz_yy = primBuffer.data(soff + 36 * idx + 15);

            auto s_xz_yz = primBuffer.data(soff + 36 * idx + 16);

            auto s_xz_zz = primBuffer.data(soff + 36 * idx + 17);

            auto s_yy_xx = primBuffer.data(soff + 36 * idx + 18);

            auto s_yy_xy = primBuffer.data(soff + 36 * idx + 19);

            auto s_yy_xz = primBuffer.data(soff + 36 * idx + 20);

            auto s_yy_yy = primBuffer.data(soff + 36 * idx + 21);

            auto s_yy_yz = primBuffer.data(soff + 36 * idx + 22);

            auto s_yy_zz = primBuffer.data(soff + 36 * idx + 23);

            auto s_yz_xx = primBuffer.data(soff + 36 * idx + 24);

            auto s_yz_xy = primBuffer.data(soff + 36 * idx + 25);

            auto s_yz_xz = primBuffer.data(soff + 36 * idx + 26);

            auto s_yz_yy = primBuffer.data(soff + 36 * idx + 27);

            auto s_yz_yz = primBuffer.data(soff + 36 * idx + 28);

            auto s_yz_zz = primBuffer.data(soff + 36 * idx + 29);

            auto s_zz_xx = primBuffer.data(soff + 36 * idx + 30);

            auto s_zz_xy = primBuffer.data(soff + 36 * idx + 31);

            auto s_zz_xz = primBuffer.data(soff + 36 * idx + 32);

            auto s_zz_yy = primBuffer.data(soff + 36 * idx + 33);

            auto s_zz_yz = primBuffer.data(soff + 36 * idx + 34);

            auto s_zz_zz = primBuffer.data(soff + 36 * idx + 35);

            // set up pointers to (D|T|D) integrals

            auto t_xx_xx = primBuffer.data(koff + 36 * idx);

            auto t_xx_xy = primBuffer.data(koff + 36 * idx + 1);

            auto t_xx_xz = primBuffer.data(koff + 36 * idx + 2);

            auto t_xx_yy = primBuffer.data(koff + 36 * idx + 3);

            auto t_xx_yz = primBuffer.data(koff + 36 * idx + 4);

            auto t_xx_zz = primBuffer.data(koff + 36 * idx + 5);

            auto t_xy_xx = primBuffer.data(koff + 36 * idx + 6);

            auto t_xy_xy = primBuffer.data(koff + 36 * idx + 7);

            auto t_xy_xz = primBuffer.data(koff + 36 * idx + 8);

            auto t_xy_yy = primBuffer.data(koff + 36 * idx + 9);

            auto t_xy_yz = primBuffer.data(koff + 36 * idx + 10);

            auto t_xy_zz = primBuffer.data(koff + 36 * idx + 11);

            auto t_xz_xx = primBuffer.data(koff + 36 * idx + 12);

            auto t_xz_xy = primBuffer.data(koff + 36 * idx + 13);

            auto t_xz_xz = primBuffer.data(koff + 36 * idx + 14);

            auto t_xz_yy = primBuffer.data(koff + 36 * idx + 15);

            auto t_xz_yz = primBuffer.data(koff + 36 * idx + 16);

            auto t_xz_zz = primBuffer.data(koff + 36 * idx + 17);

            auto t_yy_xx = primBuffer.data(koff + 36 * idx + 18);

            auto t_yy_xy = primBuffer.data(koff + 36 * idx + 19);

            auto t_yy_xz = primBuffer.data(koff + 36 * idx + 20);

            auto t_yy_yy = primBuffer.data(koff + 36 * idx + 21);

            auto t_yy_yz = primBuffer.data(koff + 36 * idx + 22);

            auto t_yy_zz = primBuffer.data(koff + 36 * idx + 23);

            auto t_yz_xx = primBuffer.data(koff + 36 * idx + 24);

            auto t_yz_xy = primBuffer.data(koff + 36 * idx + 25);

            auto t_yz_xz = primBuffer.data(koff + 36 * idx + 26);

            auto t_yz_yy = primBuffer.data(koff + 36 * idx + 27);

            auto t_yz_yz = primBuffer.data(koff + 36 * idx + 28);

            auto t_yz_zz = primBuffer.data(koff + 36 * idx + 29);

            auto t_zz_xx = primBuffer.data(koff + 36 * idx + 30);

            auto t_zz_xy = primBuffer.data(koff + 36 * idx + 31);

            auto t_zz_xz = primBuffer.data(koff + 36 * idx + 32);

            auto t_zz_yy = primBuffer.data(koff + 36 * idx + 33);

            auto t_zz_yz = primBuffer.data(koff + 36 * idx + 34);

            auto t_zz_zz = primBuffer.data(koff + 36 * idx + 35);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_x_x, s_x_y,\
                                     s_x_z, s_y_x, s_y_y, s_y_z, s_z_x, s_z_y,\
                                     s_z_z, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y,\
                                     t_y_z, t_z_x, t_z_y, t_z_z, s_0_xx, s_0_xy,\
                                     s_0_xz, s_0_yy, s_0_yz, s_0_zz, t_0_xx, t_0_xy,\
                                     t_0_xz, t_0_yy, t_0_yz, t_0_zz, s_x_xx, s_x_xy,\
                                     s_x_xz, s_x_yy, s_x_yz, s_x_zz, s_y_xx, s_y_xy,\
                                     s_y_xz, s_y_yy, s_y_yz, s_y_zz, s_z_xx, s_z_xy,\
                                     s_z_xz, s_z_yy, s_z_yz, s_z_zz, t_x_xx, t_x_xy,\
                                     t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy,\
                                     t_y_xz, t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy,\
                                     t_z_xz, t_z_yy, t_z_yz, t_z_zz, s_xx_xx, s_xx_xy,\
                                     s_xx_xz, s_xx_yy, s_xx_yz, s_xx_zz, s_xy_xx,\
                                     s_xy_xy, s_xy_xz, s_xy_yy, s_xy_yz, s_xy_zz,\
                                     s_xz_xx, s_xz_xy, s_xz_xz, s_xz_yy, s_xz_yz,\
                                     s_xz_zz, s_yy_xx, s_yy_xy, s_yy_xz, s_yy_yy,\
                                     s_yy_yz, s_yy_zz, s_yz_xx, s_yz_xy, s_yz_xz,\
                                     s_yz_yy, s_yz_yz, s_yz_zz, s_zz_xx, s_zz_xy,\
                                     s_zz_xz, s_zz_yy, s_zz_yz, s_zz_zz, t_xx_xx,\
                                     t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz,\
                                     t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz,\
                                     t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy,\
                                     t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy, t_yy_xz,\
                                     t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, t_yz_xy,\
                                     t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx,\
                                     t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xx_xx[j] = fr * s_x_xx[j] + f2t * (s_0_xx[j] + 2.0 * s_x_x[j]);

                t_xx_xx[j] = fr * t_x_xx[j] + f2t * (t_0_xx[j] + 2.0 * t_x_x[j]) + f2z * (s_xx_xx[j] - f2a * s_0_xx[j]);

                s_xx_xy[j] = fr * s_x_xy[j] + f2t * (s_0_xy[j] + s_x_y[j]);

                t_xx_xy[j] = fr * t_x_xy[j] + f2t * (t_0_xy[j] + t_x_y[j]) + f2z * (s_xx_xy[j] - f2a * s_0_xy[j]);

                s_xx_xz[j] = fr * s_x_xz[j] + f2t * (s_0_xz[j] + s_x_z[j]);

                t_xx_xz[j] = fr * t_x_xz[j] + f2t * (t_0_xz[j] + t_x_z[j]) + f2z * (s_xx_xz[j] - f2a * s_0_xz[j]);

                s_xx_yy[j] = fr * s_x_yy[j] + f2t * s_0_yy[j];

                t_xx_yy[j] = fr * t_x_yy[j] + f2t * t_0_yy[j] + f2z * (s_xx_yy[j] - f2a * s_0_yy[j]);

                s_xx_yz[j] = fr * s_x_yz[j] + f2t * s_0_yz[j];

                t_xx_yz[j] = fr * t_x_yz[j] + f2t * t_0_yz[j] + f2z * (s_xx_yz[j] - f2a * s_0_yz[j]);

                s_xx_zz[j] = fr * s_x_zz[j] + f2t * s_0_zz[j];

                t_xx_zz[j] = fr * t_x_zz[j] + f2t * t_0_zz[j] + f2z * (s_xx_zz[j] - f2a * s_0_zz[j]);

                s_xy_xx[j] = fr * s_y_xx[j] + f2t * 2.0 * s_y_x[j];

                t_xy_xx[j] = fr * t_y_xx[j] + f2t * 2.0 * t_y_x[j] + f2z * s_xy_xx[j];

                s_xy_xy[j] = fr * s_y_xy[j] + f2t * s_y_y[j];

                t_xy_xy[j] = fr * t_y_xy[j] + f2t * t_y_y[j] + f2z * s_xy_xy[j];

                s_xy_xz[j] = fr * s_y_xz[j] + f2t * s_y_z[j];

                t_xy_xz[j] = fr * t_y_xz[j] + f2t * t_y_z[j] + f2z * s_xy_xz[j];

                s_xy_yy[j] = fr * s_y_yy[j];

                t_xy_yy[j] = fr * t_y_yy[j] + f2z * s_xy_yy[j];

                s_xy_yz[j] = fr * s_y_yz[j];

                t_xy_yz[j] = fr * t_y_yz[j] + f2z * s_xy_yz[j];

                s_xy_zz[j] = fr * s_y_zz[j];

                t_xy_zz[j] = fr * t_y_zz[j] + f2z * s_xy_zz[j];

                s_xz_xx[j] = fr * s_z_xx[j] + f2t * 2.0 * s_z_x[j];

                t_xz_xx[j] = fr * t_z_xx[j] + f2t * 2.0 * t_z_x[j] + f2z * s_xz_xx[j];

                s_xz_xy[j] = fr * s_z_xy[j] + f2t * s_z_y[j];

                t_xz_xy[j] = fr * t_z_xy[j] + f2t * t_z_y[j] + f2z * s_xz_xy[j];

                s_xz_xz[j] = fr * s_z_xz[j] + f2t * s_z_z[j];

                t_xz_xz[j] = fr * t_z_xz[j] + f2t * t_z_z[j] + f2z * s_xz_xz[j];

                s_xz_yy[j] = fr * s_z_yy[j];

                t_xz_yy[j] = fr * t_z_yy[j] + f2z * s_xz_yy[j];

                s_xz_yz[j] = fr * s_z_yz[j];

                t_xz_yz[j] = fr * t_z_yz[j] + f2z * s_xz_yz[j];

                s_xz_zz[j] = fr * s_z_zz[j];

                t_xz_zz[j] = fr * t_z_zz[j] + f2z * s_xz_zz[j];

                // leading y component

                fr = pay[j];

                s_yy_xx[j] = fr * s_y_xx[j] + f2t * s_0_xx[j];

                t_yy_xx[j] = fr * t_y_xx[j] + f2t * t_0_xx[j] + f2z * (s_yy_xx[j] - f2a * s_0_xx[j]);

                s_yy_xy[j] = fr * s_y_xy[j] + f2t * (s_0_xy[j] + s_y_x[j]);

                t_yy_xy[j] = fr * t_y_xy[j] + f2t * (t_0_xy[j] + t_y_x[j]) + f2z * (s_yy_xy[j] - f2a * s_0_xy[j]);

                s_yy_xz[j] = fr * s_y_xz[j] + f2t * s_0_xz[j];

                t_yy_xz[j] = fr * t_y_xz[j] + f2t * t_0_xz[j] + f2z * (s_yy_xz[j] - f2a * s_0_xz[j]);

                s_yy_yy[j] = fr * s_y_yy[j] + f2t * (s_0_yy[j] + 2.0 * s_y_y[j]);

                t_yy_yy[j] = fr * t_y_yy[j] + f2t * (t_0_yy[j] + 2.0 * t_y_y[j]) + f2z * (s_yy_yy[j] - f2a * s_0_yy[j]);

                s_yy_yz[j] = fr * s_y_yz[j] + f2t * (s_0_yz[j] + s_y_z[j]);

                t_yy_yz[j] = fr * t_y_yz[j] + f2t * (t_0_yz[j] + t_y_z[j]) + f2z * (s_yy_yz[j] - f2a * s_0_yz[j]);

                s_yy_zz[j] = fr * s_y_zz[j] + f2t * s_0_zz[j];

                t_yy_zz[j] = fr * t_y_zz[j] + f2t * t_0_zz[j] + f2z * (s_yy_zz[j] - f2a * s_0_zz[j]);

                s_yz_xx[j] = fr * s_z_xx[j];

                t_yz_xx[j] = fr * t_z_xx[j] + f2z * s_yz_xx[j];

                s_yz_xy[j] = fr * s_z_xy[j] + f2t * s_z_x[j];

                t_yz_xy[j] = fr * t_z_xy[j] + f2t * t_z_x[j] + f2z * s_yz_xy[j];

                s_yz_xz[j] = fr * s_z_xz[j];

                t_yz_xz[j] = fr * t_z_xz[j] + f2z * s_yz_xz[j];

                s_yz_yy[j] = fr * s_z_yy[j] + f2t * 2.0 * s_z_y[j];

                t_yz_yy[j] = fr * t_z_yy[j] + f2t * 2.0 * t_z_y[j] + f2z * s_yz_yy[j];

                s_yz_yz[j] = fr * s_z_yz[j] + f2t * s_z_z[j];

                t_yz_yz[j] = fr * t_z_yz[j] + f2t * t_z_z[j] + f2z * s_yz_yz[j];

                s_yz_zz[j] = fr * s_z_zz[j];

                t_yz_zz[j] = fr * t_z_zz[j] + f2z * s_yz_zz[j];

                // leading z component

                fr = paz[j];

                s_zz_xx[j] = fr * s_z_xx[j] + f2t * s_0_xx[j];

                t_zz_xx[j] = fr * t_z_xx[j] + f2t * t_0_xx[j] + f2z * (s_zz_xx[j] - f2a * s_0_xx[j]);

                s_zz_xy[j] = fr * s_z_xy[j] + f2t * s_0_xy[j];

                t_zz_xy[j] = fr * t_z_xy[j] + f2t * t_0_xy[j] + f2z * (s_zz_xy[j] - f2a * s_0_xy[j]);

                s_zz_xz[j] = fr * s_z_xz[j] + f2t * (s_0_xz[j] + s_z_x[j]);

                t_zz_xz[j] = fr * t_z_xz[j] + f2t * (t_0_xz[j] + t_z_x[j]) + f2z * (s_zz_xz[j] - f2a * s_0_xz[j]);

                s_zz_yy[j] = fr * s_z_yy[j] + f2t * s_0_yy[j];

                t_zz_yy[j] = fr * t_z_yy[j] + f2t * t_0_yy[j] + f2z * (s_zz_yy[j] - f2a * s_0_yy[j]);

                s_zz_yz[j] = fr * s_z_yz[j] + f2t * (s_0_yz[j] + s_z_y[j]);

                t_zz_yz[j] = fr * t_z_yz[j] + f2t * (t_0_yz[j] + t_z_y[j]) + f2z * (s_zz_yz[j] - f2a * s_0_yz[j]);

                s_zz_zz[j] = fr * s_z_zz[j] + f2t * (s_0_zz[j] + 2.0 * s_z_z[j]);

                t_zz_zz[j] = fr * t_z_zz[j] + f2t * (t_0_zz[j] + 2.0 * t_z_z[j]) + f2z * (s_zz_zz[j] - f2a * s_0_zz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForSF(      CMemBlock2D<double>&  primBuffer,
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

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(s2off + 3 * idx);

            auto s_0_y = primBuffer.data(s2off + 3 * idx + 1);

            auto s_0_z = primBuffer.data(s2off + 3 * idx + 2);

            // set up pointers to (S|T|P) integrals

            auto t_0_x = primBuffer.data(k2off + 3 * idx);

            auto t_0_y = primBuffer.data(k2off + 3 * idx + 1);

            auto t_0_z = primBuffer.data(k2off + 3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(s1off + 6 * idx);

            auto s_0_xy = primBuffer.data(s1off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(s1off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(s1off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(s1off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(s1off + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(k1off + 6 * idx);

            auto t_0_xy = primBuffer.data(k1off + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(k1off + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(k1off + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(k1off + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(k1off + 6 * idx + 5);

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

            // set up pointers to (S|T|F) integrals

            auto t_0_xxx = primBuffer.data(koff + 10 * idx);

            auto t_0_xxy = primBuffer.data(koff + 10 * idx + 1);

            auto t_0_xxz = primBuffer.data(koff + 10 * idx + 2);

            auto t_0_xyy = primBuffer.data(koff + 10 * idx + 3);

            auto t_0_xyz = primBuffer.data(koff + 10 * idx + 4);

            auto t_0_xzz = primBuffer.data(koff + 10 * idx + 5);

            auto t_0_yyy = primBuffer.data(koff + 10 * idx + 6);

            auto t_0_yyz = primBuffer.data(koff + 10 * idx + 7);

            auto t_0_yzz = primBuffer.data(koff + 10 * idx + 8);

            auto t_0_zzz = primBuffer.data(koff + 10 * idx + 9);

            #pragma omp simd aligned(fx, fz, fgb, pbx, pby, pbz, s_0_x, s_0_y,\
                                     s_0_z, t_0_x, t_0_y, t_0_z, s_0_xx, s_0_xy,\
                                     s_0_xz, s_0_yy, s_0_yz, s_0_zz, t_0_xx, t_0_xy,\
                                     t_0_xz, t_0_yy, t_0_yz, t_0_zz, s_0_xxx, s_0_xxy,\
                                     s_0_xxz, s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz, t_0_xxx, t_0_xxy,\
                                     t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy,\
                                     t_0_yyz, t_0_yzz, t_0_zzz: VLX_ALIGN)
             for (int32_t j = 0; j < nprim; j++)
             {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2b = 0.50 * fgb[j];

                // leading x component

                double fr = pbx[j];

                s_0_xxx[j] = fr * s_0_xx[j] + f2t * 2.0 * s_0_x[j];

                t_0_xxx[j] = fr * t_0_xx[j] + f2t * 2.0 * t_0_x[j] + f2z * (s_0_xxx[j] - f2b * 2.0 * s_0_x[j]);

                s_0_xxy[j] = fr * s_0_xy[j] + f2t * s_0_y[j];

                t_0_xxy[j] = fr * t_0_xy[j] + f2t * t_0_y[j] + f2z * (s_0_xxy[j] - f2b * s_0_y[j]);

                s_0_xxz[j] = fr * s_0_xz[j] + f2t * s_0_z[j];

                t_0_xxz[j] = fr * t_0_xz[j] + f2t * t_0_z[j] + f2z * (s_0_xxz[j] - f2b * s_0_z[j]);

                s_0_xyy[j] = fr * s_0_yy[j];

                t_0_xyy[j] = fr * t_0_yy[j] + f2z * s_0_xyy[j];

                s_0_xyz[j] = fr * s_0_yz[j];

                t_0_xyz[j] = fr * t_0_yz[j] + f2z * s_0_xyz[j];

                s_0_xzz[j] = fr * s_0_zz[j];

                t_0_xzz[j] = fr * t_0_zz[j] + f2z * s_0_xzz[j];

                // leading y component

                fr = pby[j];

                s_0_yyy[j] = fr * s_0_yy[j] + f2t * 2.0 * s_0_y[j];

                t_0_yyy[j] = fr * t_0_yy[j] + f2t * 2.0 * t_0_y[j] + f2z * (s_0_yyy[j] - f2b * 2.0 * s_0_y[j]);

                s_0_yyz[j] = fr * s_0_yz[j] + f2t * s_0_z[j];

                t_0_yyz[j] = fr * t_0_yz[j] + f2t * t_0_z[j] + f2z * (s_0_yyz[j] - f2b * s_0_z[j]);

                s_0_yzz[j] = fr * s_0_zz[j];

                t_0_yzz[j] = fr * t_0_zz[j] + f2z * s_0_yzz[j];

                // leading z component

                fr = pbz[j];

                s_0_zzz[j] = fr * s_0_zz[j] + f2t * 2.0 * s_0_z[j];

                t_0_zzz[j] = fr * t_0_zz[j] + f2t * 2.0 * t_0_z[j] + f2z * (s_0_zzz[j] - f2b * 2.0 * s_0_z[j]);
            }

            idx++;
        }
    }
    
    
    void
    compKineticEnergyForFS(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|S) integrals

            auto s_x_0 = primBuffer.data(s2off + 3 * idx);

            auto s_y_0 = primBuffer.data(s2off + 3 * idx + 1);

            auto s_z_0 = primBuffer.data(s2off + 3 * idx + 2);

            // set up pointers to (P|T|S) integrals

            auto t_x_0 = primBuffer.data(k2off + 3 * idx);

            auto t_y_0 = primBuffer.data(k2off + 3 * idx + 1);

            auto t_z_0 = primBuffer.data(k2off + 3 * idx + 2);

            // set up pointers to (D|S) integrals

            auto s_xx_0 = primBuffer.data(s1off + 6 * idx);

            auto s_xy_0 = primBuffer.data(s1off + 6 * idx + 1);

            auto s_xz_0 = primBuffer.data(s1off + 6 * idx + 2);

            auto s_yy_0 = primBuffer.data(s1off + 6 * idx + 3);

            auto s_yz_0 = primBuffer.data(s1off + 6 * idx + 4);

            auto s_zz_0 = primBuffer.data(s1off + 6 * idx + 5);

            // set up pointers to (D|T|S) integrals

            auto t_xx_0 = primBuffer.data(k1off + 6 * idx);

            auto t_xy_0 = primBuffer.data(k1off + 6 * idx + 1);

            auto t_xz_0 = primBuffer.data(k1off + 6 * idx + 2);

            auto t_yy_0 = primBuffer.data(k1off + 6 * idx + 3);

            auto t_yz_0 = primBuffer.data(k1off + 6 * idx + 4);

            auto t_zz_0 = primBuffer.data(k1off + 6 * idx + 5);

            // set up pointers to (F|S) integrals

            auto s_xxx_0 = primBuffer.data(soff + 10 * idx);

            auto s_xxy_0 = primBuffer.data(soff + 10 * idx + 1);

            auto s_xxz_0 = primBuffer.data(soff + 10 * idx + 2);

            auto s_xyy_0 = primBuffer.data(soff + 10 * idx + 3);

            auto s_xyz_0 = primBuffer.data(soff + 10 * idx + 4);

            auto s_xzz_0 = primBuffer.data(soff + 10 * idx + 5);

            auto s_yyy_0 = primBuffer.data(soff + 10 * idx + 6);

            auto s_yyz_0 = primBuffer.data(soff + 10 * idx + 7);

            auto s_yzz_0 = primBuffer.data(soff + 10 * idx + 8);

            auto s_zzz_0 = primBuffer.data(soff + 10 * idx + 9);

            // set up pointers to (F|T|S) integrals

            auto t_xxx_0 = primBuffer.data(koff + 10 * idx);

            auto t_xxy_0 = primBuffer.data(koff + 10 * idx + 1);

            auto t_xxz_0 = primBuffer.data(koff + 10 * idx + 2);

            auto t_xyy_0 = primBuffer.data(koff + 10 * idx + 3);

            auto t_xyz_0 = primBuffer.data(koff + 10 * idx + 4);

            auto t_xzz_0 = primBuffer.data(koff + 10 * idx + 5);

            auto t_yyy_0 = primBuffer.data(koff + 10 * idx + 6);

            auto t_yyz_0 = primBuffer.data(koff + 10 * idx + 7);

            auto t_yzz_0 = primBuffer.data(koff + 10 * idx + 8);

            auto t_zzz_0 = primBuffer.data(koff + 10 * idx + 9);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_x_0, s_y_0,\
                                     s_z_0, t_x_0, t_y_0, t_z_0, s_xx_0, s_xy_0,\
                                     s_xz_0, s_yy_0, s_yz_0, s_zz_0, t_xx_0, t_xy_0,\
                                     t_xz_0, t_yy_0, t_yz_0, t_zz_0, s_xxx_0, s_xxy_0,\
                                     s_xxz_0, s_xyy_0, s_xyz_0, s_xzz_0, s_yyy_0,\
                                     s_yyz_0, s_yzz_0, s_zzz_0, t_xxx_0, t_xxy_0,\
                                     t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, t_yyy_0,\
                                     t_yyz_0, t_yzz_0, t_zzz_0: VLX_ALIGN)
             for (int32_t j = 0; j < nprim; j++)
             {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxx_0[j] = fr * s_xx_0[j] + f2t * 2.0 * s_x_0[j];

                t_xxx_0[j] = fr * t_xx_0[j] + f2t * 2.0 * t_x_0[j] + f2z * (s_xxx_0[j] - f2a * 2.0 * s_x_0[j]);

                s_xxy_0[j] = fr * s_xy_0[j] + f2t * s_y_0[j];

                t_xxy_0[j] = fr * t_xy_0[j] + f2t * t_y_0[j] + f2z * (s_xxy_0[j] - f2a * s_y_0[j]);

                s_xxz_0[j] = fr * s_xz_0[j] + f2t * s_z_0[j];

                t_xxz_0[j] = fr * t_xz_0[j] + f2t * t_z_0[j] + f2z * (s_xxz_0[j] - f2a * s_z_0[j]);

                s_xyy_0[j] = fr * s_yy_0[j];

                t_xyy_0[j] = fr * t_yy_0[j] + f2z * s_xyy_0[j];

                s_xyz_0[j] = fr * s_yz_0[j];

                t_xyz_0[j] = fr * t_yz_0[j] + f2z * s_xyz_0[j];

                s_xzz_0[j] = fr * s_zz_0[j];

                t_xzz_0[j] = fr * t_zz_0[j] + f2z * s_xzz_0[j];

                // leading y component

                fr = pay[j];

                s_yyy_0[j] = fr * s_yy_0[j] + f2t * 2.0 * s_y_0[j];

                t_yyy_0[j] = fr * t_yy_0[j] + f2t * 2.0 * t_y_0[j] + f2z * (s_yyy_0[j] - f2a * 2.0 * s_y_0[j]);

                s_yyz_0[j] = fr * s_yz_0[j] + f2t * s_z_0[j];

                t_yyz_0[j] = fr * t_yz_0[j] + f2t * t_z_0[j] + f2z * (s_yyz_0[j] - f2a * s_z_0[j]);

                s_yzz_0[j] = fr * s_zz_0[j];

                t_yzz_0[j] = fr * t_zz_0[j] + f2z * s_yzz_0[j];

                // leading z component

                fr = paz[j];

                s_zzz_0[j] = fr * s_zz_0[j] + f2t * 2.0 * s_z_0[j];

                t_zzz_0[j] = fr * t_zz_0[j] + f2t * 2.0 * t_z_0[j] + f2z * (s_zzz_0[j] - f2a * 2.0 * s_z_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForPF(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(skoff + 6 * idx);

            auto s_0_xy = primBuffer.data(skoff + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(skoff + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(skoff + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(skoff + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(skoff + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(kkoff + 6 * idx);

            auto t_0_xy = primBuffer.data(kkoff + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(kkoff + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(kkoff + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(kkoff + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(kkoff + 6 * idx + 5);

            // set up pointers to (S|F) integrals

            auto s_0_xxx = primBuffer.data(s1off + 10 * idx);

            auto s_0_xxy = primBuffer.data(s1off + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(s1off + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(s1off + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(s1off + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(s1off + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(s1off + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(s1off + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(s1off + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(s1off + 10 * idx + 9);

            // set up pointers to (S|T|F) integrals

            auto t_0_xxx = primBuffer.data(k1off + 10 * idx);

            auto t_0_xxy = primBuffer.data(k1off + 10 * idx + 1);

            auto t_0_xxz = primBuffer.data(k1off + 10 * idx + 2);

            auto t_0_xyy = primBuffer.data(k1off + 10 * idx + 3);

            auto t_0_xyz = primBuffer.data(k1off + 10 * idx + 4);

            auto t_0_xzz = primBuffer.data(k1off + 10 * idx + 5);

            auto t_0_yyy = primBuffer.data(k1off + 10 * idx + 6);

            auto t_0_yyz = primBuffer.data(k1off + 10 * idx + 7);

            auto t_0_yzz = primBuffer.data(k1off + 10 * idx + 8);

            auto t_0_zzz = primBuffer.data(k1off + 10 * idx + 9);

            // set up pointers to (P|F) integrals

            auto s_x_xxx = primBuffer.data(soff + 30 * idx);

            auto s_x_xxy = primBuffer.data(soff + 30 * idx + 1);

            auto s_x_xxz = primBuffer.data(soff + 30 * idx + 2);

            auto s_x_xyy = primBuffer.data(soff + 30 * idx + 3);

            auto s_x_xyz = primBuffer.data(soff + 30 * idx + 4);

            auto s_x_xzz = primBuffer.data(soff + 30 * idx + 5);

            auto s_x_yyy = primBuffer.data(soff + 30 * idx + 6);

            auto s_x_yyz = primBuffer.data(soff + 30 * idx + 7);

            auto s_x_yzz = primBuffer.data(soff + 30 * idx + 8);

            auto s_x_zzz = primBuffer.data(soff + 30 * idx + 9);

            auto s_y_xxx = primBuffer.data(soff + 30 * idx + 10);

            auto s_y_xxy = primBuffer.data(soff + 30 * idx + 11);

            auto s_y_xxz = primBuffer.data(soff + 30 * idx + 12);

            auto s_y_xyy = primBuffer.data(soff + 30 * idx + 13);

            auto s_y_xyz = primBuffer.data(soff + 30 * idx + 14);

            auto s_y_xzz = primBuffer.data(soff + 30 * idx + 15);

            auto s_y_yyy = primBuffer.data(soff + 30 * idx + 16);

            auto s_y_yyz = primBuffer.data(soff + 30 * idx + 17);

            auto s_y_yzz = primBuffer.data(soff + 30 * idx + 18);

            auto s_y_zzz = primBuffer.data(soff + 30 * idx + 19);

            auto s_z_xxx = primBuffer.data(soff + 30 * idx + 20);

            auto s_z_xxy = primBuffer.data(soff + 30 * idx + 21);

            auto s_z_xxz = primBuffer.data(soff + 30 * idx + 22);

            auto s_z_xyy = primBuffer.data(soff + 30 * idx + 23);

            auto s_z_xyz = primBuffer.data(soff + 30 * idx + 24);

            auto s_z_xzz = primBuffer.data(soff + 30 * idx + 25);

            auto s_z_yyy = primBuffer.data(soff + 30 * idx + 26);

            auto s_z_yyz = primBuffer.data(soff + 30 * idx + 27);

            auto s_z_yzz = primBuffer.data(soff + 30 * idx + 28);

            auto s_z_zzz = primBuffer.data(soff + 30 * idx + 29);

            // set up pointers to (P|T|F) integrals

            auto t_x_xxx = primBuffer.data(koff + 30 * idx);

            auto t_x_xxy = primBuffer.data(koff + 30 * idx + 1);

            auto t_x_xxz = primBuffer.data(koff + 30 * idx + 2);

            auto t_x_xyy = primBuffer.data(koff + 30 * idx + 3);

            auto t_x_xyz = primBuffer.data(koff + 30 * idx + 4);

            auto t_x_xzz = primBuffer.data(koff + 30 * idx + 5);

            auto t_x_yyy = primBuffer.data(koff + 30 * idx + 6);

            auto t_x_yyz = primBuffer.data(koff + 30 * idx + 7);

            auto t_x_yzz = primBuffer.data(koff + 30 * idx + 8);

            auto t_x_zzz = primBuffer.data(koff + 30 * idx + 9);

            auto t_y_xxx = primBuffer.data(koff + 30 * idx + 10);

            auto t_y_xxy = primBuffer.data(koff + 30 * idx + 11);

            auto t_y_xxz = primBuffer.data(koff + 30 * idx + 12);

            auto t_y_xyy = primBuffer.data(koff + 30 * idx + 13);

            auto t_y_xyz = primBuffer.data(koff + 30 * idx + 14);

            auto t_y_xzz = primBuffer.data(koff + 30 * idx + 15);

            auto t_y_yyy = primBuffer.data(koff + 30 * idx + 16);

            auto t_y_yyz = primBuffer.data(koff + 30 * idx + 17);

            auto t_y_yzz = primBuffer.data(koff + 30 * idx + 18);

            auto t_y_zzz = primBuffer.data(koff + 30 * idx + 19);

            auto t_z_xxx = primBuffer.data(koff + 30 * idx + 20);

            auto t_z_xxy = primBuffer.data(koff + 30 * idx + 21);

            auto t_z_xxz = primBuffer.data(koff + 30 * idx + 22);

            auto t_z_xyy = primBuffer.data(koff + 30 * idx + 23);

            auto t_z_xyz = primBuffer.data(koff + 30 * idx + 24);

            auto t_z_xzz = primBuffer.data(koff + 30 * idx + 25);

            auto t_z_yyy = primBuffer.data(koff + 30 * idx + 26);

            auto t_z_yyz = primBuffer.data(koff + 30 * idx + 27);

            auto t_z_yzz = primBuffer.data(koff + 30 * idx + 28);

            auto t_z_zzz = primBuffer.data(koff + 30 * idx + 29);

            #pragma omp simd aligned(fx, fz, pax, pay, paz, s_0_xx, s_0_xy, s_0_xz,\
                                     s_0_yy, s_0_yz, s_0_zz, t_0_xx, t_0_xy, t_0_xz,\
                                     t_0_yy, t_0_yz, t_0_zz, s_0_xxx, s_0_xxy,\
                                     s_0_xxz, s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz, t_0_xxx, t_0_xxy,\
                                     t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy,\
                                     t_0_yyz, t_0_yzz, t_0_zzz, s_x_xxx, s_x_xxy,\
                                     s_x_xxz, s_x_xyy, s_x_xyz, s_x_xzz, s_x_yyy,\
                                     s_x_yyz, s_x_yzz, s_x_zzz, s_y_xxx, s_y_xxy,\
                                     s_y_xxz, s_y_xyy, s_y_xyz, s_y_xzz, s_y_yyy,\
                                     s_y_yyz, s_y_yzz, s_y_zzz, s_z_xxx, s_z_xxy,\
                                     s_z_xxz, s_z_xyy, s_z_xyz, s_z_xzz, s_z_yyy,\
                                     s_z_yyz, s_z_yzz, s_z_zzz, t_x_xxx, t_x_xxy,\
                                     t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy,\
                                     t_x_yyz, t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy,\
                                     t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy,\
                                     t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy,\
                                     t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy,\
                                     t_z_yyz, t_z_yzz, t_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pax[j];

                s_x_xxx[j] = fr * s_0_xxx[j] + f2t * 3.0 * s_0_xx[j];

                t_x_xxx[j] = fr * t_0_xxx[j] + f2t * 3.0 * t_0_xx[j] + f2z * s_x_xxx[j];

                s_x_xxy[j] = fr * s_0_xxy[j] + f2t * 2.0 * s_0_xy[j];

                t_x_xxy[j] = fr * t_0_xxy[j] + f2t * 2.0 * t_0_xy[j] + f2z * s_x_xxy[j];

                s_x_xxz[j] = fr * s_0_xxz[j] + f2t * 2.0 * s_0_xz[j];

                t_x_xxz[j] = fr * t_0_xxz[j] + f2t * 2.0 * t_0_xz[j] + f2z * s_x_xxz[j];

                s_x_xyy[j] = fr * s_0_xyy[j] + f2t * s_0_yy[j];

                t_x_xyy[j] = fr * t_0_xyy[j] + f2t * t_0_yy[j] + f2z * s_x_xyy[j];

                s_x_xyz[j] = fr * s_0_xyz[j] + f2t * s_0_yz[j];

                t_x_xyz[j] = fr * t_0_xyz[j] + f2t * t_0_yz[j] + f2z * s_x_xyz[j];

                s_x_xzz[j] = fr * s_0_xzz[j] + f2t * s_0_zz[j];

                t_x_xzz[j] = fr * t_0_xzz[j] + f2t * t_0_zz[j] + f2z * s_x_xzz[j];

                s_x_yyy[j] = fr * s_0_yyy[j];

                t_x_yyy[j] = fr * t_0_yyy[j] + f2z * s_x_yyy[j];

                s_x_yyz[j] = fr * s_0_yyz[j];

                t_x_yyz[j] = fr * t_0_yyz[j] + f2z * s_x_yyz[j];

                s_x_yzz[j] = fr * s_0_yzz[j];

                t_x_yzz[j] = fr * t_0_yzz[j] + f2z * s_x_yzz[j];

                s_x_zzz[j] = fr * s_0_zzz[j];

                t_x_zzz[j] = fr * t_0_zzz[j] + f2z * s_x_zzz[j];

                // leading y component

                fr = pay[j];

                s_y_xxx[j] = fr * s_0_xxx[j];

                t_y_xxx[j] = fr * t_0_xxx[j] + f2z * s_y_xxx[j];

                s_y_xxy[j] = fr * s_0_xxy[j] + f2t * s_0_xx[j];

                t_y_xxy[j] = fr * t_0_xxy[j] + f2t * t_0_xx[j] + f2z * s_y_xxy[j];

                s_y_xxz[j] = fr * s_0_xxz[j];

                t_y_xxz[j] = fr * t_0_xxz[j] + f2z * s_y_xxz[j];

                s_y_xyy[j] = fr * s_0_xyy[j] + f2t * 2.0 * s_0_xy[j];

                t_y_xyy[j] = fr * t_0_xyy[j] + f2t * 2.0 * t_0_xy[j] + f2z * s_y_xyy[j];

                s_y_xyz[j] = fr * s_0_xyz[j] + f2t * s_0_xz[j];

                t_y_xyz[j] = fr * t_0_xyz[j] + f2t * t_0_xz[j] + f2z * s_y_xyz[j];

                s_y_xzz[j] = fr * s_0_xzz[j];

                t_y_xzz[j] = fr * t_0_xzz[j] + f2z * s_y_xzz[j];

                s_y_yyy[j] = fr * s_0_yyy[j] + f2t * 3.0 * s_0_yy[j];

                t_y_yyy[j] = fr * t_0_yyy[j] + f2t * 3.0 * t_0_yy[j] + f2z * s_y_yyy[j];

                s_y_yyz[j] = fr * s_0_yyz[j] + f2t * 2.0 * s_0_yz[j];

                t_y_yyz[j] = fr * t_0_yyz[j] + f2t * 2.0 * t_0_yz[j] + f2z * s_y_yyz[j];

                s_y_yzz[j] = fr * s_0_yzz[j] + f2t * s_0_zz[j];

                t_y_yzz[j] = fr * t_0_yzz[j] + f2t * t_0_zz[j] + f2z * s_y_yzz[j];

                s_y_zzz[j] = fr * s_0_zzz[j];

                t_y_zzz[j] = fr * t_0_zzz[j] + f2z * s_y_zzz[j];

                // leading z component

                fr = paz[j];

                s_z_xxx[j] = fr * s_0_xxx[j];

                t_z_xxx[j] = fr * t_0_xxx[j] + f2z * s_z_xxx[j];

                s_z_xxy[j] = fr * s_0_xxy[j];

                t_z_xxy[j] = fr * t_0_xxy[j] + f2z * s_z_xxy[j];

                s_z_xxz[j] = fr * s_0_xxz[j] + f2t * s_0_xx[j];

                t_z_xxz[j] = fr * t_0_xxz[j] + f2t * t_0_xx[j] + f2z * s_z_xxz[j];

                s_z_xyy[j] = fr * s_0_xyy[j];

                t_z_xyy[j] = fr * t_0_xyy[j] + f2z * s_z_xyy[j];

                s_z_xyz[j] = fr * s_0_xyz[j] + f2t * s_0_xy[j];

                t_z_xyz[j] = fr * t_0_xyz[j] + f2t * t_0_xy[j] + f2z * s_z_xyz[j];

                s_z_xzz[j] = fr * s_0_xzz[j] + f2t * 2.0 * s_0_xz[j];

                t_z_xzz[j] = fr * t_0_xzz[j] + f2t * 2.0 * t_0_xz[j] + f2z * s_z_xzz[j];

                s_z_yyy[j] = fr * s_0_yyy[j];

                t_z_yyy[j] = fr * t_0_yyy[j] + f2z * s_z_yyy[j];

                s_z_yyz[j] = fr * s_0_yyz[j] + f2t * s_0_yy[j];

                t_z_yyz[j] = fr * t_0_yyz[j] + f2t * t_0_yy[j] + f2z * s_z_yyz[j];

                s_z_yzz[j] = fr * s_0_yzz[j] + f2t * 2.0 * s_0_yz[j];

                t_z_yzz[j] = fr * t_0_yzz[j] + f2t * 2.0 * t_0_yz[j] + f2z * s_z_yzz[j];

                s_z_zzz[j] = fr * s_0_zzz[j] + f2t * 3.0 * s_0_zz[j];

                t_z_zzz[j] = fr * t_0_zzz[j] + f2t * 3.0 * t_0_zz[j] + f2z * s_z_zzz[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFP(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|S) integrals

            auto s_xx_0 = primBuffer.data(skoff + 6 * idx);

            auto s_xy_0 = primBuffer.data(skoff + 6 * idx + 1);

            auto s_xz_0 = primBuffer.data(skoff + 6 * idx + 2);

            auto s_yy_0 = primBuffer.data(skoff + 6 * idx + 3);

            auto s_yz_0 = primBuffer.data(skoff + 6 * idx + 4);

            auto s_zz_0 = primBuffer.data(skoff + 6 * idx + 5);

            // set up pointers to (D|T|S) integrals

            auto t_xx_0 = primBuffer.data(kkoff + 6 * idx);

            auto t_xy_0 = primBuffer.data(kkoff + 6 * idx + 1);

            auto t_xz_0 = primBuffer.data(kkoff + 6 * idx + 2);

            auto t_yy_0 = primBuffer.data(kkoff + 6 * idx + 3);

            auto t_yz_0 = primBuffer.data(kkoff + 6 * idx + 4);

            auto t_zz_0 = primBuffer.data(kkoff + 6 * idx + 5);

            // set up pointers to (P|P) integrals

            auto s_x_x = primBuffer.data(s2off + 9 * idx);

            auto s_x_y = primBuffer.data(s2off + 9 * idx + 1);

            auto s_x_z = primBuffer.data(s2off + 9 * idx + 2);

            auto s_y_x = primBuffer.data(s2off + 9 * idx + 3);

            auto s_y_y = primBuffer.data(s2off + 9 * idx + 4);

            auto s_y_z = primBuffer.data(s2off + 9 * idx + 5);

            auto s_z_x = primBuffer.data(s2off + 9 * idx + 6);

            auto s_z_y = primBuffer.data(s2off + 9 * idx + 7);

            auto s_z_z = primBuffer.data(s2off + 9 * idx + 8);

            // set up pointers to (P|T|P) integrals

            auto t_x_x = primBuffer.data(k2off + 9 * idx);

            auto t_x_y = primBuffer.data(k2off + 9 * idx + 1);

            auto t_x_z = primBuffer.data(k2off + 9 * idx + 2);

            auto t_y_x = primBuffer.data(k2off + 9 * idx + 3);

            auto t_y_y = primBuffer.data(k2off + 9 * idx + 4);

            auto t_y_z = primBuffer.data(k2off + 9 * idx + 5);

            auto t_z_x = primBuffer.data(k2off + 9 * idx + 6);

            auto t_z_y = primBuffer.data(k2off + 9 * idx + 7);

            auto t_z_z = primBuffer.data(k2off + 9 * idx + 8);

            // set up pointers to (D|P) integrals

            auto s_xx_x = primBuffer.data(s1off + 18 * idx);

            auto s_xx_y = primBuffer.data(s1off + 18 * idx + 1);

            auto s_xx_z = primBuffer.data(s1off + 18 * idx + 2);

            auto s_xy_x = primBuffer.data(s1off + 18 * idx + 3);

            auto s_xy_y = primBuffer.data(s1off + 18 * idx + 4);

            auto s_xy_z = primBuffer.data(s1off + 18 * idx + 5);

            auto s_xz_x = primBuffer.data(s1off + 18 * idx + 6);

            auto s_xz_y = primBuffer.data(s1off + 18 * idx + 7);

            auto s_xz_z = primBuffer.data(s1off + 18 * idx + 8);

            auto s_yy_x = primBuffer.data(s1off + 18 * idx + 9);

            auto s_yy_y = primBuffer.data(s1off + 18 * idx + 10);

            auto s_yy_z = primBuffer.data(s1off + 18 * idx + 11);

            auto s_yz_x = primBuffer.data(s1off + 18 * idx + 12);

            auto s_yz_y = primBuffer.data(s1off + 18 * idx + 13);

            auto s_yz_z = primBuffer.data(s1off + 18 * idx + 14);

            auto s_zz_x = primBuffer.data(s1off + 18 * idx + 15);

            auto s_zz_y = primBuffer.data(s1off + 18 * idx + 16);

            auto s_zz_z = primBuffer.data(s1off + 18 * idx + 17);

            // set up pointers to (D|T|P) integrals

            auto t_xx_x = primBuffer.data(k1off + 18 * idx);

            auto t_xx_y = primBuffer.data(k1off + 18 * idx + 1);

            auto t_xx_z = primBuffer.data(k1off + 18 * idx + 2);

            auto t_xy_x = primBuffer.data(k1off + 18 * idx + 3);

            auto t_xy_y = primBuffer.data(k1off + 18 * idx + 4);

            auto t_xy_z = primBuffer.data(k1off + 18 * idx + 5);

            auto t_xz_x = primBuffer.data(k1off + 18 * idx + 6);

            auto t_xz_y = primBuffer.data(k1off + 18 * idx + 7);

            auto t_xz_z = primBuffer.data(k1off + 18 * idx + 8);

            auto t_yy_x = primBuffer.data(k1off + 18 * idx + 9);

            auto t_yy_y = primBuffer.data(k1off + 18 * idx + 10);

            auto t_yy_z = primBuffer.data(k1off + 18 * idx + 11);

            auto t_yz_x = primBuffer.data(k1off + 18 * idx + 12);

            auto t_yz_y = primBuffer.data(k1off + 18 * idx + 13);

            auto t_yz_z = primBuffer.data(k1off + 18 * idx + 14);

            auto t_zz_x = primBuffer.data(k1off + 18 * idx + 15);

            auto t_zz_y = primBuffer.data(k1off + 18 * idx + 16);

            auto t_zz_z = primBuffer.data(k1off + 18 * idx + 17);

            // set up pointers to (F|P) integrals

            auto s_xxx_x = primBuffer.data(soff + 30 * idx);

            auto s_xxx_y = primBuffer.data(soff + 30 * idx + 1);

            auto s_xxx_z = primBuffer.data(soff + 30 * idx + 2);

            auto s_xxy_x = primBuffer.data(soff + 30 * idx + 3);

            auto s_xxy_y = primBuffer.data(soff + 30 * idx + 4);

            auto s_xxy_z = primBuffer.data(soff + 30 * idx + 5);

            auto s_xxz_x = primBuffer.data(soff + 30 * idx + 6);

            auto s_xxz_y = primBuffer.data(soff + 30 * idx + 7);

            auto s_xxz_z = primBuffer.data(soff + 30 * idx + 8);

            auto s_xyy_x = primBuffer.data(soff + 30 * idx + 9);

            auto s_xyy_y = primBuffer.data(soff + 30 * idx + 10);

            auto s_xyy_z = primBuffer.data(soff + 30 * idx + 11);

            auto s_xyz_x = primBuffer.data(soff + 30 * idx + 12);

            auto s_xyz_y = primBuffer.data(soff + 30 * idx + 13);

            auto s_xyz_z = primBuffer.data(soff + 30 * idx + 14);

            auto s_xzz_x = primBuffer.data(soff + 30 * idx + 15);

            auto s_xzz_y = primBuffer.data(soff + 30 * idx + 16);

            auto s_xzz_z = primBuffer.data(soff + 30 * idx + 17);

            auto s_yyy_x = primBuffer.data(soff + 30 * idx + 18);

            auto s_yyy_y = primBuffer.data(soff + 30 * idx + 19);

            auto s_yyy_z = primBuffer.data(soff + 30 * idx + 20);

            auto s_yyz_x = primBuffer.data(soff + 30 * idx + 21);

            auto s_yyz_y = primBuffer.data(soff + 30 * idx + 22);

            auto s_yyz_z = primBuffer.data(soff + 30 * idx + 23);

            auto s_yzz_x = primBuffer.data(soff + 30 * idx + 24);

            auto s_yzz_y = primBuffer.data(soff + 30 * idx + 25);

            auto s_yzz_z = primBuffer.data(soff + 30 * idx + 26);

            auto s_zzz_x = primBuffer.data(soff + 30 * idx + 27);

            auto s_zzz_y = primBuffer.data(soff + 30 * idx + 28);

            auto s_zzz_z = primBuffer.data(soff + 30 * idx + 29);

            // set up pointers to (F|T|P) integrals

            auto t_xxx_x = primBuffer.data(koff + 30 * idx);

            auto t_xxx_y = primBuffer.data(koff + 30 * idx + 1);

            auto t_xxx_z = primBuffer.data(koff + 30 * idx + 2);

            auto t_xxy_x = primBuffer.data(koff + 30 * idx + 3);

            auto t_xxy_y = primBuffer.data(koff + 30 * idx + 4);

            auto t_xxy_z = primBuffer.data(koff + 30 * idx + 5);

            auto t_xxz_x = primBuffer.data(koff + 30 * idx + 6);

            auto t_xxz_y = primBuffer.data(koff + 30 * idx + 7);

            auto t_xxz_z = primBuffer.data(koff + 30 * idx + 8);

            auto t_xyy_x = primBuffer.data(koff + 30 * idx + 9);

            auto t_xyy_y = primBuffer.data(koff + 30 * idx + 10);

            auto t_xyy_z = primBuffer.data(koff + 30 * idx + 11);

            auto t_xyz_x = primBuffer.data(koff + 30 * idx + 12);

            auto t_xyz_y = primBuffer.data(koff + 30 * idx + 13);

            auto t_xyz_z = primBuffer.data(koff + 30 * idx + 14);

            auto t_xzz_x = primBuffer.data(koff + 30 * idx + 15);

            auto t_xzz_y = primBuffer.data(koff + 30 * idx + 16);

            auto t_xzz_z = primBuffer.data(koff + 30 * idx + 17);

            auto t_yyy_x = primBuffer.data(koff + 30 * idx + 18);

            auto t_yyy_y = primBuffer.data(koff + 30 * idx + 19);

            auto t_yyy_z = primBuffer.data(koff + 30 * idx + 20);

            auto t_yyz_x = primBuffer.data(koff + 30 * idx + 21);

            auto t_yyz_y = primBuffer.data(koff + 30 * idx + 22);

            auto t_yyz_z = primBuffer.data(koff + 30 * idx + 23);

            auto t_yzz_x = primBuffer.data(koff + 30 * idx + 24);

            auto t_yzz_y = primBuffer.data(koff + 30 * idx + 25);

            auto t_yzz_z = primBuffer.data(koff + 30 * idx + 26);

            auto t_zzz_x = primBuffer.data(koff + 30 * idx + 27);

            auto t_zzz_y = primBuffer.data(koff + 30 * idx + 28);

            auto t_zzz_z = primBuffer.data(koff + 30 * idx + 29);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xx_0, s_xy_0,\
                                     s_xz_0, s_yy_0, s_yz_0, s_zz_0, t_xx_0, t_xy_0,\
                                     t_xz_0, t_yy_0, t_yz_0, t_zz_0, s_x_x, s_x_y,\
                                     s_x_z, s_y_x, s_y_y, s_y_z, s_z_x, s_z_y,\
                                     s_z_z, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y,\
                                     t_y_z, t_z_x, t_z_y, t_z_z, s_xx_x, s_xx_y,\
                                     s_xx_z, s_xy_x, s_xy_y, s_xy_z, s_xz_x, s_xz_y,\
                                     s_xz_z, s_yy_x, s_yy_y, s_yy_z, s_yz_x, s_yz_y,\
                                     s_yz_z, s_zz_x, s_zz_y, s_zz_z, t_xx_x, t_xx_y,\
                                     t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y,\
                                     t_xz_z, t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y,\
                                     t_yz_z, t_zz_x, t_zz_y, t_zz_z, s_xxx_x, s_xxx_y,\
                                     s_xxx_z, s_xxy_x, s_xxy_y, s_xxy_z, s_xxz_x,\
                                     s_xxz_y, s_xxz_z, s_xyy_x, s_xyy_y, s_xyy_z,\
                                     s_xyz_x, s_xyz_y, s_xyz_z, s_xzz_x, s_xzz_y,\
                                     s_xzz_z, s_yyy_x, s_yyy_y, s_yyy_z, s_yyz_x,\
                                     s_yyz_y, s_yyz_z, s_yzz_x, s_yzz_y, s_yzz_z,\
                                     s_zzz_x, s_zzz_y, s_zzz_z, t_xxx_x, t_xxx_y,\
                                     t_xxx_z, t_xxy_x, t_xxy_y, t_xxy_z, t_xxz_x,\
                                     t_xxz_y, t_xxz_z, t_xyy_x, t_xyy_y, t_xyy_z,\
                                     t_xyz_x, t_xyz_y, t_xyz_z, t_xzz_x, t_xzz_y,\
                                     t_xzz_z, t_yyy_x, t_yyy_y, t_yyy_z, t_yyz_x,\
                                     t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z,\
                                     t_zzz_x, t_zzz_y, t_zzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxx_x[j] = fr * s_xx_x[j] + f2t * (2.0 * s_x_x[j] + s_xx_0[j]);

                t_xxx_x[j] = fr * t_xx_x[j] + f2t * (2.0 * t_x_x[j] + t_xx_0[j]) + f2z * (s_xxx_x[j] - f2a * 2.0 * s_x_x[j]);

                s_xxx_y[j] = fr * s_xx_y[j] + f2t * 2.0 * s_x_y[j];

                t_xxx_y[j] = fr * t_xx_y[j] + f2t * 2.0 * t_x_y[j] + f2z * (s_xxx_y[j] - f2a * 2.0 * s_x_y[j]);

                s_xxx_z[j] = fr * s_xx_z[j] + f2t * 2.0 * s_x_z[j];

                t_xxx_z[j] = fr * t_xx_z[j] + f2t * 2.0 * t_x_z[j] + f2z * (s_xxx_z[j] - f2a * 2.0 * s_x_z[j]);

                s_xxy_x[j] = fr * s_xy_x[j] + f2t * (s_y_x[j] + s_xy_0[j]);

                t_xxy_x[j] = fr * t_xy_x[j] + f2t * (t_y_x[j] + t_xy_0[j]) + f2z * (s_xxy_x[j] - f2a * s_y_x[j]);

                s_xxy_y[j] = fr * s_xy_y[j] + f2t * s_y_y[j];

                t_xxy_y[j] = fr * t_xy_y[j] + f2t * t_y_y[j] + f2z * (s_xxy_y[j] - f2a * s_y_y[j]);

                s_xxy_z[j] = fr * s_xy_z[j] + f2t * s_y_z[j];

                t_xxy_z[j] = fr * t_xy_z[j] + f2t * t_y_z[j] + f2z * (s_xxy_z[j] - f2a * s_y_z[j]);

                s_xxz_x[j] = fr * s_xz_x[j] + f2t * (s_z_x[j] + s_xz_0[j]);

                t_xxz_x[j] = fr * t_xz_x[j] + f2t * (t_z_x[j] + t_xz_0[j]) + f2z * (s_xxz_x[j] - f2a * s_z_x[j]);

                s_xxz_y[j] = fr * s_xz_y[j] + f2t * s_z_y[j];

                t_xxz_y[j] = fr * t_xz_y[j] + f2t * t_z_y[j] + f2z * (s_xxz_y[j] - f2a * s_z_y[j]);

                s_xxz_z[j] = fr * s_xz_z[j] + f2t * s_z_z[j];

                t_xxz_z[j] = fr * t_xz_z[j] + f2t * t_z_z[j] + f2z * (s_xxz_z[j] - f2a * s_z_z[j]);

                s_xyy_x[j] = fr * s_yy_x[j] + f2t * s_yy_0[j];

                t_xyy_x[j] = fr * t_yy_x[j] + f2t * t_yy_0[j] + f2z * s_xyy_x[j];

                s_xyy_y[j] = fr * s_yy_y[j];

                t_xyy_y[j] = fr * t_yy_y[j] + f2z * s_xyy_y[j];

                s_xyy_z[j] = fr * s_yy_z[j];

                t_xyy_z[j] = fr * t_yy_z[j] + f2z * s_xyy_z[j];

                s_xyz_x[j] = fr * s_yz_x[j] + f2t * s_yz_0[j];

                t_xyz_x[j] = fr * t_yz_x[j] + f2t * t_yz_0[j] + f2z * s_xyz_x[j];

                s_xyz_y[j] = fr * s_yz_y[j];

                t_xyz_y[j] = fr * t_yz_y[j] + f2z * s_xyz_y[j];

                s_xyz_z[j] = fr * s_yz_z[j];

                t_xyz_z[j] = fr * t_yz_z[j] + f2z * s_xyz_z[j];

                s_xzz_x[j] = fr * s_zz_x[j] + f2t * s_zz_0[j];

                t_xzz_x[j] = fr * t_zz_x[j] + f2t * t_zz_0[j] + f2z * s_xzz_x[j];

                s_xzz_y[j] = fr * s_zz_y[j];

                t_xzz_y[j] = fr * t_zz_y[j] + f2z * s_xzz_y[j];

                s_xzz_z[j] = fr * s_zz_z[j];

                t_xzz_z[j] = fr * t_zz_z[j] + f2z * s_xzz_z[j];

                // leading y component

                fr = pay[j];

                s_yyy_x[j] = fr * s_yy_x[j] + f2t * 2.0 * s_y_x[j];

                t_yyy_x[j] = fr * t_yy_x[j] + f2t * 2.0 * t_y_x[j] + f2z * (s_yyy_x[j] - f2a * 2.0 * s_y_x[j]);

                s_yyy_y[j] = fr * s_yy_y[j] + f2t * (2.0 * s_y_y[j] + s_yy_0[j]);

                t_yyy_y[j] = fr * t_yy_y[j] + f2t * (2.0 * t_y_y[j] + t_yy_0[j]) + f2z * (s_yyy_y[j] - f2a * 2.0 * s_y_y[j]);

                s_yyy_z[j] = fr * s_yy_z[j] + f2t * 2.0 * s_y_z[j];

                t_yyy_z[j] = fr * t_yy_z[j] + f2t * 2.0 * t_y_z[j] + f2z * (s_yyy_z[j] - f2a * 2.0 * s_y_z[j]);

                s_yyz_x[j] = fr * s_yz_x[j] + f2t * s_z_x[j];

                t_yyz_x[j] = fr * t_yz_x[j] + f2t * t_z_x[j] + f2z * (s_yyz_x[j] - f2a * s_z_x[j]);

                s_yyz_y[j] = fr * s_yz_y[j] + f2t * (s_z_y[j] + s_yz_0[j]);

                t_yyz_y[j] = fr * t_yz_y[j] + f2t * (t_z_y[j] + t_yz_0[j]) + f2z * (s_yyz_y[j] - f2a * s_z_y[j]);

                s_yyz_z[j] = fr * s_yz_z[j] + f2t * s_z_z[j];

                t_yyz_z[j] = fr * t_yz_z[j] + f2t * t_z_z[j] + f2z * (s_yyz_z[j] - f2a * s_z_z[j]);

                s_yzz_x[j] = fr * s_zz_x[j];

                t_yzz_x[j] = fr * t_zz_x[j] + f2z * s_yzz_x[j];

                s_yzz_y[j] = fr * s_zz_y[j] + f2t * s_zz_0[j];

                t_yzz_y[j] = fr * t_zz_y[j] + f2t * t_zz_0[j] + f2z * s_yzz_y[j];

                s_yzz_z[j] = fr * s_zz_z[j];

                t_yzz_z[j] = fr * t_zz_z[j] + f2z * s_yzz_z[j];

                // leading z component

                fr = paz[j];

                s_zzz_x[j] = fr * s_zz_x[j] + f2t * 2.0 * s_z_x[j];

                t_zzz_x[j] = fr * t_zz_x[j] + f2t * 2.0 * t_z_x[j] + f2z * (s_zzz_x[j] - f2a * 2.0 * s_z_x[j]);

                s_zzz_y[j] = fr * s_zz_y[j] + f2t * 2.0 * s_z_y[j];

                t_zzz_y[j] = fr * t_zz_y[j] + f2t * 2.0 * t_z_y[j] + f2z * (s_zzz_y[j] - f2a * 2.0 * s_z_y[j]);

                s_zzz_z[j] = fr * s_zz_z[j] + f2t * (2.0 * s_z_z[j] + s_zz_0[j]);

                t_zzz_z[j] = fr * t_zz_z[j] + f2t * (2.0 * t_z_z[j] + t_zz_0[j]) + f2z * (s_zzz_z[j] - f2a * 2.0 * s_z_z[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDF(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|D) integrals

            auto s_x_xx = primBuffer.data(skoff + 18 * idx);

            auto s_x_xy = primBuffer.data(skoff + 18 * idx + 1);

            auto s_x_xz = primBuffer.data(skoff + 18 * idx + 2);

            auto s_x_yy = primBuffer.data(skoff + 18 * idx + 3);

            auto s_x_yz = primBuffer.data(skoff + 18 * idx + 4);

            auto s_x_zz = primBuffer.data(skoff + 18 * idx + 5);

            auto s_y_xx = primBuffer.data(skoff + 18 * idx + 6);

            auto s_y_xy = primBuffer.data(skoff + 18 * idx + 7);

            auto s_y_xz = primBuffer.data(skoff + 18 * idx + 8);

            auto s_y_yy = primBuffer.data(skoff + 18 * idx + 9);

            auto s_y_yz = primBuffer.data(skoff + 18 * idx + 10);

            auto s_y_zz = primBuffer.data(skoff + 18 * idx + 11);

            auto s_z_xx = primBuffer.data(skoff + 18 * idx + 12);

            auto s_z_xy = primBuffer.data(skoff + 18 * idx + 13);

            auto s_z_xz = primBuffer.data(skoff + 18 * idx + 14);

            auto s_z_yy = primBuffer.data(skoff + 18 * idx + 15);

            auto s_z_yz = primBuffer.data(skoff + 18 * idx + 16);

            auto s_z_zz = primBuffer.data(skoff + 18 * idx + 17);

            // set up pointers to (P|T|D) integrals

            auto t_x_xx = primBuffer.data(kkoff + 18 * idx);

            auto t_x_xy = primBuffer.data(kkoff + 18 * idx + 1);

            auto t_x_xz = primBuffer.data(kkoff + 18 * idx + 2);

            auto t_x_yy = primBuffer.data(kkoff + 18 * idx + 3);

            auto t_x_yz = primBuffer.data(kkoff + 18 * idx + 4);

            auto t_x_zz = primBuffer.data(kkoff + 18 * idx + 5);

            auto t_y_xx = primBuffer.data(kkoff + 18 * idx + 6);

            auto t_y_xy = primBuffer.data(kkoff + 18 * idx + 7);

            auto t_y_xz = primBuffer.data(kkoff + 18 * idx + 8);

            auto t_y_yy = primBuffer.data(kkoff + 18 * idx + 9);

            auto t_y_yz = primBuffer.data(kkoff + 18 * idx + 10);

            auto t_y_zz = primBuffer.data(kkoff + 18 * idx + 11);

            auto t_z_xx = primBuffer.data(kkoff + 18 * idx + 12);

            auto t_z_xy = primBuffer.data(kkoff + 18 * idx + 13);

            auto t_z_xz = primBuffer.data(kkoff + 18 * idx + 14);

            auto t_z_yy = primBuffer.data(kkoff + 18 * idx + 15);

            auto t_z_yz = primBuffer.data(kkoff + 18 * idx + 16);

            auto t_z_zz = primBuffer.data(kkoff + 18 * idx + 17);

            // set up pointers to (S|F) integrals

            auto s_0_xxx = primBuffer.data(s2off + 10 * idx);

            auto s_0_xxy = primBuffer.data(s2off + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(s2off + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(s2off + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(s2off + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(s2off + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(s2off + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(s2off + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(s2off + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(s2off + 10 * idx + 9);

            // set up pointers to (S|T|F) integrals

            auto t_0_xxx = primBuffer.data(k2off + 10 * idx);

            auto t_0_xxy = primBuffer.data(k2off + 10 * idx + 1);

            auto t_0_xxz = primBuffer.data(k2off + 10 * idx + 2);

            auto t_0_xyy = primBuffer.data(k2off + 10 * idx + 3);

            auto t_0_xyz = primBuffer.data(k2off + 10 * idx + 4);

            auto t_0_xzz = primBuffer.data(k2off + 10 * idx + 5);

            auto t_0_yyy = primBuffer.data(k2off + 10 * idx + 6);

            auto t_0_yyz = primBuffer.data(k2off + 10 * idx + 7);

            auto t_0_yzz = primBuffer.data(k2off + 10 * idx + 8);

            auto t_0_zzz = primBuffer.data(k2off + 10 * idx + 9);

            // set up pointers to (P|F) integrals

            auto s_x_xxx = primBuffer.data(s1off + 30 * idx);

            auto s_x_xxy = primBuffer.data(s1off + 30 * idx + 1);

            auto s_x_xxz = primBuffer.data(s1off + 30 * idx + 2);

            auto s_x_xyy = primBuffer.data(s1off + 30 * idx + 3);

            auto s_x_xyz = primBuffer.data(s1off + 30 * idx + 4);

            auto s_x_xzz = primBuffer.data(s1off + 30 * idx + 5);

            auto s_x_yyy = primBuffer.data(s1off + 30 * idx + 6);

            auto s_x_yyz = primBuffer.data(s1off + 30 * idx + 7);

            auto s_x_yzz = primBuffer.data(s1off + 30 * idx + 8);

            auto s_x_zzz = primBuffer.data(s1off + 30 * idx + 9);

            auto s_y_xxx = primBuffer.data(s1off + 30 * idx + 10);

            auto s_y_xxy = primBuffer.data(s1off + 30 * idx + 11);

            auto s_y_xxz = primBuffer.data(s1off + 30 * idx + 12);

            auto s_y_xyy = primBuffer.data(s1off + 30 * idx + 13);

            auto s_y_xyz = primBuffer.data(s1off + 30 * idx + 14);

            auto s_y_xzz = primBuffer.data(s1off + 30 * idx + 15);

            auto s_y_yyy = primBuffer.data(s1off + 30 * idx + 16);

            auto s_y_yyz = primBuffer.data(s1off + 30 * idx + 17);

            auto s_y_yzz = primBuffer.data(s1off + 30 * idx + 18);

            auto s_y_zzz = primBuffer.data(s1off + 30 * idx + 19);

            auto s_z_xxx = primBuffer.data(s1off + 30 * idx + 20);

            auto s_z_xxy = primBuffer.data(s1off + 30 * idx + 21);

            auto s_z_xxz = primBuffer.data(s1off + 30 * idx + 22);

            auto s_z_xyy = primBuffer.data(s1off + 30 * idx + 23);

            auto s_z_xyz = primBuffer.data(s1off + 30 * idx + 24);

            auto s_z_xzz = primBuffer.data(s1off + 30 * idx + 25);

            auto s_z_yyy = primBuffer.data(s1off + 30 * idx + 26);

            auto s_z_yyz = primBuffer.data(s1off + 30 * idx + 27);

            auto s_z_yzz = primBuffer.data(s1off + 30 * idx + 28);

            auto s_z_zzz = primBuffer.data(s1off + 30 * idx + 29);

            // set up pointers to (P|T|F) integrals

            auto t_x_xxx = primBuffer.data(k1off + 30 * idx);

            auto t_x_xxy = primBuffer.data(k1off + 30 * idx + 1);

            auto t_x_xxz = primBuffer.data(k1off + 30 * idx + 2);

            auto t_x_xyy = primBuffer.data(k1off + 30 * idx + 3);

            auto t_x_xyz = primBuffer.data(k1off + 30 * idx + 4);

            auto t_x_xzz = primBuffer.data(k1off + 30 * idx + 5);

            auto t_x_yyy = primBuffer.data(k1off + 30 * idx + 6);

            auto t_x_yyz = primBuffer.data(k1off + 30 * idx + 7);

            auto t_x_yzz = primBuffer.data(k1off + 30 * idx + 8);

            auto t_x_zzz = primBuffer.data(k1off + 30 * idx + 9);

            auto t_y_xxx = primBuffer.data(k1off + 30 * idx + 10);

            auto t_y_xxy = primBuffer.data(k1off + 30 * idx + 11);

            auto t_y_xxz = primBuffer.data(k1off + 30 * idx + 12);

            auto t_y_xyy = primBuffer.data(k1off + 30 * idx + 13);

            auto t_y_xyz = primBuffer.data(k1off + 30 * idx + 14);

            auto t_y_xzz = primBuffer.data(k1off + 30 * idx + 15);

            auto t_y_yyy = primBuffer.data(k1off + 30 * idx + 16);

            auto t_y_yyz = primBuffer.data(k1off + 30 * idx + 17);

            auto t_y_yzz = primBuffer.data(k1off + 30 * idx + 18);

            auto t_y_zzz = primBuffer.data(k1off + 30 * idx + 19);

            auto t_z_xxx = primBuffer.data(k1off + 30 * idx + 20);

            auto t_z_xxy = primBuffer.data(k1off + 30 * idx + 21);

            auto t_z_xxz = primBuffer.data(k1off + 30 * idx + 22);

            auto t_z_xyy = primBuffer.data(k1off + 30 * idx + 23);

            auto t_z_xyz = primBuffer.data(k1off + 30 * idx + 24);

            auto t_z_xzz = primBuffer.data(k1off + 30 * idx + 25);

            auto t_z_yyy = primBuffer.data(k1off + 30 * idx + 26);

            auto t_z_yyz = primBuffer.data(k1off + 30 * idx + 27);

            auto t_z_yzz = primBuffer.data(k1off + 30 * idx + 28);

            auto t_z_zzz = primBuffer.data(k1off + 30 * idx + 29);

            // set up pointers to (D|F) integrals

            auto s_xx_xxx = primBuffer.data(soff + 60 * idx);

            auto s_xx_xxy = primBuffer.data(soff + 60 * idx + 1);

            auto s_xx_xxz = primBuffer.data(soff + 60 * idx + 2);

            auto s_xx_xyy = primBuffer.data(soff + 60 * idx + 3);

            auto s_xx_xyz = primBuffer.data(soff + 60 * idx + 4);

            auto s_xx_xzz = primBuffer.data(soff + 60 * idx + 5);

            auto s_xx_yyy = primBuffer.data(soff + 60 * idx + 6);

            auto s_xx_yyz = primBuffer.data(soff + 60 * idx + 7);

            auto s_xx_yzz = primBuffer.data(soff + 60 * idx + 8);

            auto s_xx_zzz = primBuffer.data(soff + 60 * idx + 9);

            auto s_xy_xxx = primBuffer.data(soff + 60 * idx + 10);

            auto s_xy_xxy = primBuffer.data(soff + 60 * idx + 11);

            auto s_xy_xxz = primBuffer.data(soff + 60 * idx + 12);

            auto s_xy_xyy = primBuffer.data(soff + 60 * idx + 13);

            auto s_xy_xyz = primBuffer.data(soff + 60 * idx + 14);

            auto s_xy_xzz = primBuffer.data(soff + 60 * idx + 15);

            auto s_xy_yyy = primBuffer.data(soff + 60 * idx + 16);

            auto s_xy_yyz = primBuffer.data(soff + 60 * idx + 17);

            auto s_xy_yzz = primBuffer.data(soff + 60 * idx + 18);

            auto s_xy_zzz = primBuffer.data(soff + 60 * idx + 19);

            auto s_xz_xxx = primBuffer.data(soff + 60 * idx + 20);

            auto s_xz_xxy = primBuffer.data(soff + 60 * idx + 21);

            auto s_xz_xxz = primBuffer.data(soff + 60 * idx + 22);

            auto s_xz_xyy = primBuffer.data(soff + 60 * idx + 23);

            auto s_xz_xyz = primBuffer.data(soff + 60 * idx + 24);

            auto s_xz_xzz = primBuffer.data(soff + 60 * idx + 25);

            auto s_xz_yyy = primBuffer.data(soff + 60 * idx + 26);

            auto s_xz_yyz = primBuffer.data(soff + 60 * idx + 27);

            auto s_xz_yzz = primBuffer.data(soff + 60 * idx + 28);

            auto s_xz_zzz = primBuffer.data(soff + 60 * idx + 29);

            auto s_yy_xxx = primBuffer.data(soff + 60 * idx + 30);

            auto s_yy_xxy = primBuffer.data(soff + 60 * idx + 31);

            auto s_yy_xxz = primBuffer.data(soff + 60 * idx + 32);

            auto s_yy_xyy = primBuffer.data(soff + 60 * idx + 33);

            auto s_yy_xyz = primBuffer.data(soff + 60 * idx + 34);

            auto s_yy_xzz = primBuffer.data(soff + 60 * idx + 35);

            auto s_yy_yyy = primBuffer.data(soff + 60 * idx + 36);

            auto s_yy_yyz = primBuffer.data(soff + 60 * idx + 37);

            auto s_yy_yzz = primBuffer.data(soff + 60 * idx + 38);

            auto s_yy_zzz = primBuffer.data(soff + 60 * idx + 39);

            auto s_yz_xxx = primBuffer.data(soff + 60 * idx + 40);

            auto s_yz_xxy = primBuffer.data(soff + 60 * idx + 41);

            auto s_yz_xxz = primBuffer.data(soff + 60 * idx + 42);

            auto s_yz_xyy = primBuffer.data(soff + 60 * idx + 43);

            auto s_yz_xyz = primBuffer.data(soff + 60 * idx + 44);

            auto s_yz_xzz = primBuffer.data(soff + 60 * idx + 45);

            auto s_yz_yyy = primBuffer.data(soff + 60 * idx + 46);

            auto s_yz_yyz = primBuffer.data(soff + 60 * idx + 47);

            auto s_yz_yzz = primBuffer.data(soff + 60 * idx + 48);

            auto s_yz_zzz = primBuffer.data(soff + 60 * idx + 49);

            auto s_zz_xxx = primBuffer.data(soff + 60 * idx + 50);

            auto s_zz_xxy = primBuffer.data(soff + 60 * idx + 51);

            auto s_zz_xxz = primBuffer.data(soff + 60 * idx + 52);

            auto s_zz_xyy = primBuffer.data(soff + 60 * idx + 53);

            auto s_zz_xyz = primBuffer.data(soff + 60 * idx + 54);

            auto s_zz_xzz = primBuffer.data(soff + 60 * idx + 55);

            auto s_zz_yyy = primBuffer.data(soff + 60 * idx + 56);

            auto s_zz_yyz = primBuffer.data(soff + 60 * idx + 57);

            auto s_zz_yzz = primBuffer.data(soff + 60 * idx + 58);

            auto s_zz_zzz = primBuffer.data(soff + 60 * idx + 59);

            // set up pointers to (D|T|F) integrals

            auto t_xx_xxx = primBuffer.data(koff + 60 * idx);

            auto t_xx_xxy = primBuffer.data(koff + 60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(koff + 60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(koff + 60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(koff + 60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(koff + 60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(koff + 60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(koff + 60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(koff + 60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(koff + 60 * idx + 9);

            auto t_xy_xxx = primBuffer.data(koff + 60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(koff + 60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(koff + 60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(koff + 60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(koff + 60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(koff + 60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(koff + 60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(koff + 60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(koff + 60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(koff + 60 * idx + 19);

            auto t_xz_xxx = primBuffer.data(koff + 60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(koff + 60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(koff + 60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(koff + 60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(koff + 60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(koff + 60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(koff + 60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(koff + 60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(koff + 60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(koff + 60 * idx + 29);

            auto t_yy_xxx = primBuffer.data(koff + 60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(koff + 60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(koff + 60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(koff + 60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(koff + 60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(koff + 60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(koff + 60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(koff + 60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(koff + 60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(koff + 60 * idx + 39);

            auto t_yz_xxx = primBuffer.data(koff + 60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(koff + 60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(koff + 60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(koff + 60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(koff + 60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(koff + 60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(koff + 60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(koff + 60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(koff + 60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(koff + 60 * idx + 49);

            auto t_zz_xxx = primBuffer.data(koff + 60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(koff + 60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(koff + 60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(koff + 60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(koff + 60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(koff + 60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(koff + 60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(koff + 60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(koff + 60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(koff + 60 * idx + 59);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_x_xx, s_x_xy,\
                                     s_x_xz, s_x_yy, s_x_yz, s_x_zz, s_y_xx, s_y_xy,\
                                     s_y_xz, s_y_yy, s_y_yz, s_y_zz, s_z_xx, s_z_xy,\
                                     s_z_xz, s_z_yy, s_z_yz, s_z_zz, t_x_xx, t_x_xy,\
                                     t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy,\
                                     t_y_xz, t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy,\
                                     t_z_xz, t_z_yy, t_z_yz, t_z_zz, s_0_xxx, s_0_xxy,\
                                     s_0_xxz, s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz, t_0_xxx, t_0_xxy,\
                                     t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy,\
                                     t_0_yyz, t_0_yzz, t_0_zzz, s_x_xxx, s_x_xxy,\
                                     s_x_xxz, s_x_xyy, s_x_xyz, s_x_xzz, s_x_yyy,\
                                     s_x_yyz, s_x_yzz, s_x_zzz, s_y_xxx, s_y_xxy,\
                                     s_y_xxz, s_y_xyy, s_y_xyz, s_y_xzz, s_y_yyy,\
                                     s_y_yyz, s_y_yzz, s_y_zzz, s_z_xxx, s_z_xxy,\
                                     s_z_xxz, s_z_xyy, s_z_xyz, s_z_xzz, s_z_yyy,\
                                     s_z_yyz, s_z_yzz, s_z_zzz, t_x_xxx, t_x_xxy,\
                                     t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy,\
                                     t_x_yyz, t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy,\
                                     t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy,\
                                     t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy,\
                                     t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy,\
                                     t_z_yyz, t_z_yzz, t_z_zzz, s_xx_xxx, s_xx_xxy,\
                                     s_xx_xxz, s_xx_xyy, s_xx_xyz, s_xx_xzz, s_xx_yyy,\
                                     s_xx_yyz, s_xx_yzz, s_xx_zzz, s_xy_xxx, s_xy_xxy,\
                                     s_xy_xxz, s_xy_xyy, s_xy_xyz, s_xy_xzz, s_xy_yyy,\
                                     s_xy_yyz, s_xy_yzz, s_xy_zzz, s_xz_xxx, s_xz_xxy,\
                                     s_xz_xxz, s_xz_xyy, s_xz_xyz, s_xz_xzz, s_xz_yyy,\
                                     s_xz_yyz, s_xz_yzz, s_xz_zzz, s_yy_xxx, s_yy_xxy,\
                                     s_yy_xxz, s_yy_xyy, s_yy_xyz, s_yy_xzz, s_yy_yyy,\
                                     s_yy_yyz, s_yy_yzz, s_yy_zzz, s_yz_xxx, s_yz_xxy,\
                                     s_yz_xxz, s_yz_xyy, s_yz_xyz, s_yz_xzz, s_yz_yyy,\
                                     s_yz_yyz, s_yz_yzz, s_yz_zzz, s_zz_xxx, s_zz_xxy,\
                                     s_zz_xxz, s_zz_xyy, s_zz_xyz, s_zz_xzz, s_zz_yyy,\
                                     s_zz_yyz, s_zz_yzz, s_zz_zzz, t_xx_xxx, t_xx_xxy,\
                                     t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, t_xx_yyy,\
                                     t_xx_yyz, t_xx_yzz, t_xx_zzz, t_xy_xxx, t_xy_xxy,\
                                     t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy,\
                                     t_xy_yyz, t_xy_yzz, t_xy_zzz, t_xz_xxx, t_xz_xxy,\
                                     t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy,\
                                     t_xz_yyz, t_xz_yzz, t_xz_zzz, t_yy_xxx, t_yy_xxy,\
                                     t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, t_yy_yyy,\
                                     t_yy_yyz, t_yy_yzz, t_yy_zzz, t_yz_xxx, t_yz_xxy,\
                                     t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy,\
                                     t_yz_yyz, t_yz_yzz, t_yz_zzz, t_zz_xxx, t_zz_xxy,\
                                     t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, t_zz_yyy,\
                                     t_zz_yyz, t_zz_yzz, t_zz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xx_xxx[j] = fr * s_x_xxx[j] + f2t * (s_0_xxx[j] + 3.0 * s_x_xx[j]);

                t_xx_xxx[j] = fr * t_x_xxx[j] + f2t * (t_0_xxx[j] + 3.0 * t_x_xx[j]) + f2z * (s_xx_xxx[j] - f2a * s_0_xxx[j]);

                s_xx_xxy[j] = fr * s_x_xxy[j] + f2t * (s_0_xxy[j] + 2.0 * s_x_xy[j]);

                t_xx_xxy[j] = fr * t_x_xxy[j] + f2t * (t_0_xxy[j] + 2.0 * t_x_xy[j]) + f2z * (s_xx_xxy[j] - f2a * s_0_xxy[j]);

                s_xx_xxz[j] = fr * s_x_xxz[j] + f2t * (s_0_xxz[j] + 2.0 * s_x_xz[j]);

                t_xx_xxz[j] = fr * t_x_xxz[j] + f2t * (t_0_xxz[j] + 2.0 * t_x_xz[j]) + f2z * (s_xx_xxz[j] - f2a * s_0_xxz[j]);

                s_xx_xyy[j] = fr * s_x_xyy[j] + f2t * (s_0_xyy[j] + s_x_yy[j]);

                t_xx_xyy[j] = fr * t_x_xyy[j] + f2t * (t_0_xyy[j] + t_x_yy[j]) + f2z * (s_xx_xyy[j] - f2a * s_0_xyy[j]);

                s_xx_xyz[j] = fr * s_x_xyz[j] + f2t * (s_0_xyz[j] + s_x_yz[j]);

                t_xx_xyz[j] = fr * t_x_xyz[j] + f2t * (t_0_xyz[j] + t_x_yz[j]) + f2z * (s_xx_xyz[j] - f2a * s_0_xyz[j]);

                s_xx_xzz[j] = fr * s_x_xzz[j] + f2t * (s_0_xzz[j] + s_x_zz[j]);

                t_xx_xzz[j] = fr * t_x_xzz[j] + f2t * (t_0_xzz[j] + t_x_zz[j]) + f2z * (s_xx_xzz[j] - f2a * s_0_xzz[j]);

                s_xx_yyy[j] = fr * s_x_yyy[j] + f2t * s_0_yyy[j];

                t_xx_yyy[j] = fr * t_x_yyy[j] + f2t * t_0_yyy[j] + f2z * (s_xx_yyy[j] - f2a * s_0_yyy[j]);

                s_xx_yyz[j] = fr * s_x_yyz[j] + f2t * s_0_yyz[j];

                t_xx_yyz[j] = fr * t_x_yyz[j] + f2t * t_0_yyz[j] + f2z * (s_xx_yyz[j] - f2a * s_0_yyz[j]);

                s_xx_yzz[j] = fr * s_x_yzz[j] + f2t * s_0_yzz[j];

                t_xx_yzz[j] = fr * t_x_yzz[j] + f2t * t_0_yzz[j] + f2z * (s_xx_yzz[j] - f2a * s_0_yzz[j]);

                s_xx_zzz[j] = fr * s_x_zzz[j] + f2t * s_0_zzz[j];

                t_xx_zzz[j] = fr * t_x_zzz[j] + f2t * t_0_zzz[j] + f2z * (s_xx_zzz[j] - f2a * s_0_zzz[j]);

                s_xy_xxx[j] = fr * s_y_xxx[j] + f2t * 3.0 * s_y_xx[j];

                t_xy_xxx[j] = fr * t_y_xxx[j] + f2t * 3.0 * t_y_xx[j] + f2z * s_xy_xxx[j];

                s_xy_xxy[j] = fr * s_y_xxy[j] + f2t * 2.0 * s_y_xy[j];

                t_xy_xxy[j] = fr * t_y_xxy[j] + f2t * 2.0 * t_y_xy[j] + f2z * s_xy_xxy[j];

                s_xy_xxz[j] = fr * s_y_xxz[j] + f2t * 2.0 * s_y_xz[j];

                t_xy_xxz[j] = fr * t_y_xxz[j] + f2t * 2.0 * t_y_xz[j] + f2z * s_xy_xxz[j];

                s_xy_xyy[j] = fr * s_y_xyy[j] + f2t * s_y_yy[j];

                t_xy_xyy[j] = fr * t_y_xyy[j] + f2t * t_y_yy[j] + f2z * s_xy_xyy[j];

                s_xy_xyz[j] = fr * s_y_xyz[j] + f2t * s_y_yz[j];

                t_xy_xyz[j] = fr * t_y_xyz[j] + f2t * t_y_yz[j] + f2z * s_xy_xyz[j];

                s_xy_xzz[j] = fr * s_y_xzz[j] + f2t * s_y_zz[j];

                t_xy_xzz[j] = fr * t_y_xzz[j] + f2t * t_y_zz[j] + f2z * s_xy_xzz[j];

                s_xy_yyy[j] = fr * s_y_yyy[j];

                t_xy_yyy[j] = fr * t_y_yyy[j] + f2z * s_xy_yyy[j];

                s_xy_yyz[j] = fr * s_y_yyz[j];

                t_xy_yyz[j] = fr * t_y_yyz[j] + f2z * s_xy_yyz[j];

                s_xy_yzz[j] = fr * s_y_yzz[j];

                t_xy_yzz[j] = fr * t_y_yzz[j] + f2z * s_xy_yzz[j];

                s_xy_zzz[j] = fr * s_y_zzz[j];

                t_xy_zzz[j] = fr * t_y_zzz[j] + f2z * s_xy_zzz[j];

                s_xz_xxx[j] = fr * s_z_xxx[j] + f2t * 3.0 * s_z_xx[j];

                t_xz_xxx[j] = fr * t_z_xxx[j] + f2t * 3.0 * t_z_xx[j] + f2z * s_xz_xxx[j];

                s_xz_xxy[j] = fr * s_z_xxy[j] + f2t * 2.0 * s_z_xy[j];

                t_xz_xxy[j] = fr * t_z_xxy[j] + f2t * 2.0 * t_z_xy[j] + f2z * s_xz_xxy[j];

                s_xz_xxz[j] = fr * s_z_xxz[j] + f2t * 2.0 * s_z_xz[j];

                t_xz_xxz[j] = fr * t_z_xxz[j] + f2t * 2.0 * t_z_xz[j] + f2z * s_xz_xxz[j];

                s_xz_xyy[j] = fr * s_z_xyy[j] + f2t * s_z_yy[j];

                t_xz_xyy[j] = fr * t_z_xyy[j] + f2t * t_z_yy[j] + f2z * s_xz_xyy[j];

                s_xz_xyz[j] = fr * s_z_xyz[j] + f2t * s_z_yz[j];

                t_xz_xyz[j] = fr * t_z_xyz[j] + f2t * t_z_yz[j] + f2z * s_xz_xyz[j];

                s_xz_xzz[j] = fr * s_z_xzz[j] + f2t * s_z_zz[j];

                t_xz_xzz[j] = fr * t_z_xzz[j] + f2t * t_z_zz[j] + f2z * s_xz_xzz[j];

                s_xz_yyy[j] = fr * s_z_yyy[j];

                t_xz_yyy[j] = fr * t_z_yyy[j] + f2z * s_xz_yyy[j];

                s_xz_yyz[j] = fr * s_z_yyz[j];

                t_xz_yyz[j] = fr * t_z_yyz[j] + f2z * s_xz_yyz[j];

                s_xz_yzz[j] = fr * s_z_yzz[j];

                t_xz_yzz[j] = fr * t_z_yzz[j] + f2z * s_xz_yzz[j];

                s_xz_zzz[j] = fr * s_z_zzz[j];

                t_xz_zzz[j] = fr * t_z_zzz[j] + f2z * s_xz_zzz[j];

                // leading y component

                fr = pay[j];

                s_yy_xxx[j] = fr * s_y_xxx[j] + f2t * s_0_xxx[j];

                t_yy_xxx[j] = fr * t_y_xxx[j] + f2t * t_0_xxx[j] + f2z * (s_yy_xxx[j] - f2a * s_0_xxx[j]);

                s_yy_xxy[j] = fr * s_y_xxy[j] + f2t * (s_0_xxy[j] + s_y_xx[j]);

                t_yy_xxy[j] = fr * t_y_xxy[j] + f2t * (t_0_xxy[j] + t_y_xx[j]) + f2z * (s_yy_xxy[j] - f2a * s_0_xxy[j]);

                s_yy_xxz[j] = fr * s_y_xxz[j] + f2t * s_0_xxz[j];

                t_yy_xxz[j] = fr * t_y_xxz[j] + f2t * t_0_xxz[j] + f2z * (s_yy_xxz[j] - f2a * s_0_xxz[j]);

                s_yy_xyy[j] = fr * s_y_xyy[j] + f2t * (s_0_xyy[j] + 2.0 * s_y_xy[j]);

                t_yy_xyy[j] = fr * t_y_xyy[j] + f2t * (t_0_xyy[j] + 2.0 * t_y_xy[j]) + f2z * (s_yy_xyy[j] - f2a * s_0_xyy[j]);

                s_yy_xyz[j] = fr * s_y_xyz[j] + f2t * (s_0_xyz[j] + s_y_xz[j]);

                t_yy_xyz[j] = fr * t_y_xyz[j] + f2t * (t_0_xyz[j] + t_y_xz[j]) + f2z * (s_yy_xyz[j] - f2a * s_0_xyz[j]);

                s_yy_xzz[j] = fr * s_y_xzz[j] + f2t * s_0_xzz[j];

                t_yy_xzz[j] = fr * t_y_xzz[j] + f2t * t_0_xzz[j] + f2z * (s_yy_xzz[j] - f2a * s_0_xzz[j]);

                s_yy_yyy[j] = fr * s_y_yyy[j] + f2t * (s_0_yyy[j] + 3.0 * s_y_yy[j]);

                t_yy_yyy[j] = fr * t_y_yyy[j] + f2t * (t_0_yyy[j] + 3.0 * t_y_yy[j]) + f2z * (s_yy_yyy[j] - f2a * s_0_yyy[j]);

                s_yy_yyz[j] = fr * s_y_yyz[j] + f2t * (s_0_yyz[j] + 2.0 * s_y_yz[j]);

                t_yy_yyz[j] = fr * t_y_yyz[j] + f2t * (t_0_yyz[j] + 2.0 * t_y_yz[j]) + f2z * (s_yy_yyz[j] - f2a * s_0_yyz[j]);

                s_yy_yzz[j] = fr * s_y_yzz[j] + f2t * (s_0_yzz[j] + s_y_zz[j]);

                t_yy_yzz[j] = fr * t_y_yzz[j] + f2t * (t_0_yzz[j] + t_y_zz[j]) + f2z * (s_yy_yzz[j] - f2a * s_0_yzz[j]);

                s_yy_zzz[j] = fr * s_y_zzz[j] + f2t * s_0_zzz[j];

                t_yy_zzz[j] = fr * t_y_zzz[j] + f2t * t_0_zzz[j] + f2z * (s_yy_zzz[j] - f2a * s_0_zzz[j]);

                s_yz_xxx[j] = fr * s_z_xxx[j];

                t_yz_xxx[j] = fr * t_z_xxx[j] + f2z * s_yz_xxx[j];

                s_yz_xxy[j] = fr * s_z_xxy[j] + f2t * s_z_xx[j];

                t_yz_xxy[j] = fr * t_z_xxy[j] + f2t * t_z_xx[j] + f2z * s_yz_xxy[j];

                s_yz_xxz[j] = fr * s_z_xxz[j];

                t_yz_xxz[j] = fr * t_z_xxz[j] + f2z * s_yz_xxz[j];

                s_yz_xyy[j] = fr * s_z_xyy[j] + f2t * 2.0 * s_z_xy[j];

                t_yz_xyy[j] = fr * t_z_xyy[j] + f2t * 2.0 * t_z_xy[j] + f2z * s_yz_xyy[j];

                s_yz_xyz[j] = fr * s_z_xyz[j] + f2t * s_z_xz[j];

                t_yz_xyz[j] = fr * t_z_xyz[j] + f2t * t_z_xz[j] + f2z * s_yz_xyz[j];

                s_yz_xzz[j] = fr * s_z_xzz[j];

                t_yz_xzz[j] = fr * t_z_xzz[j] + f2z * s_yz_xzz[j];

                s_yz_yyy[j] = fr * s_z_yyy[j] + f2t * 3.0 * s_z_yy[j];

                t_yz_yyy[j] = fr * t_z_yyy[j] + f2t * 3.0 * t_z_yy[j] + f2z * s_yz_yyy[j];

                s_yz_yyz[j] = fr * s_z_yyz[j] + f2t * 2.0 * s_z_yz[j];

                t_yz_yyz[j] = fr * t_z_yyz[j] + f2t * 2.0 * t_z_yz[j] + f2z * s_yz_yyz[j];

                s_yz_yzz[j] = fr * s_z_yzz[j] + f2t * s_z_zz[j];

                t_yz_yzz[j] = fr * t_z_yzz[j] + f2t * t_z_zz[j] + f2z * s_yz_yzz[j];

                s_yz_zzz[j] = fr * s_z_zzz[j];

                t_yz_zzz[j] = fr * t_z_zzz[j] + f2z * s_yz_zzz[j];

                // leading z component

                fr = paz[j];

                s_zz_xxx[j] = fr * s_z_xxx[j] + f2t * s_0_xxx[j];

                t_zz_xxx[j] = fr * t_z_xxx[j] + f2t * t_0_xxx[j] + f2z * (s_zz_xxx[j] - f2a * s_0_xxx[j]);

                s_zz_xxy[j] = fr * s_z_xxy[j] + f2t * s_0_xxy[j];

                t_zz_xxy[j] = fr * t_z_xxy[j] + f2t * t_0_xxy[j] + f2z * (s_zz_xxy[j] - f2a * s_0_xxy[j]);

                s_zz_xxz[j] = fr * s_z_xxz[j] + f2t * (s_0_xxz[j] + s_z_xx[j]);

                t_zz_xxz[j] = fr * t_z_xxz[j] + f2t * (t_0_xxz[j] + t_z_xx[j]) + f2z * (s_zz_xxz[j] - f2a * s_0_xxz[j]);

                s_zz_xyy[j] = fr * s_z_xyy[j] + f2t * s_0_xyy[j];

                t_zz_xyy[j] = fr * t_z_xyy[j] + f2t * t_0_xyy[j] + f2z * (s_zz_xyy[j] - f2a * s_0_xyy[j]);

                s_zz_xyz[j] = fr * s_z_xyz[j] + f2t * (s_0_xyz[j] + s_z_xy[j]);

                t_zz_xyz[j] = fr * t_z_xyz[j] + f2t * (t_0_xyz[j] + t_z_xy[j]) + f2z * (s_zz_xyz[j] - f2a * s_0_xyz[j]);

                s_zz_xzz[j] = fr * s_z_xzz[j] + f2t * (s_0_xzz[j] + 2.0 * s_z_xz[j]);

                t_zz_xzz[j] = fr * t_z_xzz[j] + f2t * (t_0_xzz[j] + 2.0 * t_z_xz[j]) + f2z * (s_zz_xzz[j] - f2a * s_0_xzz[j]);

                s_zz_yyy[j] = fr * s_z_yyy[j] + f2t * s_0_yyy[j];

                t_zz_yyy[j] = fr * t_z_yyy[j] + f2t * t_0_yyy[j] + f2z * (s_zz_yyy[j] - f2a * s_0_yyy[j]);

                s_zz_yyz[j] = fr * s_z_yyz[j] + f2t * (s_0_yyz[j] + s_z_yy[j]);

                t_zz_yyz[j] = fr * t_z_yyz[j] + f2t * (t_0_yyz[j] + t_z_yy[j]) + f2z * (s_zz_yyz[j] - f2a * s_0_yyz[j]);

                s_zz_yzz[j] = fr * s_z_yzz[j] + f2t * (s_0_yzz[j] + 2.0 * s_z_yz[j]);

                t_zz_yzz[j] = fr * t_z_yzz[j] + f2t * (t_0_yzz[j] + 2.0 * t_z_yz[j]) + f2z * (s_zz_yzz[j] - f2a * s_0_yzz[j]);

                s_zz_zzz[j] = fr * s_z_zzz[j] + f2t * (s_0_zzz[j] + 3.0 * s_z_zz[j]);

                t_zz_zzz[j] = fr * t_z_zzz[j] + f2t * (t_0_zzz[j] + 3.0 * t_z_zz[j]) + f2z * (s_zz_zzz[j] - f2a * s_0_zzz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFD(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|P) integrals

            auto s_xx_x = primBuffer.data(skoff + 18 * idx);

            auto s_xx_y = primBuffer.data(skoff + 18 * idx + 1);

            auto s_xx_z = primBuffer.data(skoff + 18 * idx + 2);

            auto s_xy_x = primBuffer.data(skoff + 18 * idx + 3);

            auto s_xy_y = primBuffer.data(skoff + 18 * idx + 4);

            auto s_xy_z = primBuffer.data(skoff + 18 * idx + 5);

            auto s_xz_x = primBuffer.data(skoff + 18 * idx + 6);

            auto s_xz_y = primBuffer.data(skoff + 18 * idx + 7);

            auto s_xz_z = primBuffer.data(skoff + 18 * idx + 8);

            auto s_yy_x = primBuffer.data(skoff + 18 * idx + 9);

            auto s_yy_y = primBuffer.data(skoff + 18 * idx + 10);

            auto s_yy_z = primBuffer.data(skoff + 18 * idx + 11);

            auto s_yz_x = primBuffer.data(skoff + 18 * idx + 12);

            auto s_yz_y = primBuffer.data(skoff + 18 * idx + 13);

            auto s_yz_z = primBuffer.data(skoff + 18 * idx + 14);

            auto s_zz_x = primBuffer.data(skoff + 18 * idx + 15);

            auto s_zz_y = primBuffer.data(skoff + 18 * idx + 16);

            auto s_zz_z = primBuffer.data(skoff + 18 * idx + 17);

            // set up pointers to (D|T|P) integrals

            auto t_xx_x = primBuffer.data(kkoff + 18 * idx);

            auto t_xx_y = primBuffer.data(kkoff + 18 * idx + 1);

            auto t_xx_z = primBuffer.data(kkoff + 18 * idx + 2);

            auto t_xy_x = primBuffer.data(kkoff + 18 * idx + 3);

            auto t_xy_y = primBuffer.data(kkoff + 18 * idx + 4);

            auto t_xy_z = primBuffer.data(kkoff + 18 * idx + 5);

            auto t_xz_x = primBuffer.data(kkoff + 18 * idx + 6);

            auto t_xz_y = primBuffer.data(kkoff + 18 * idx + 7);

            auto t_xz_z = primBuffer.data(kkoff + 18 * idx + 8);

            auto t_yy_x = primBuffer.data(kkoff + 18 * idx + 9);

            auto t_yy_y = primBuffer.data(kkoff + 18 * idx + 10);

            auto t_yy_z = primBuffer.data(kkoff + 18 * idx + 11);

            auto t_yz_x = primBuffer.data(kkoff + 18 * idx + 12);

            auto t_yz_y = primBuffer.data(kkoff + 18 * idx + 13);

            auto t_yz_z = primBuffer.data(kkoff + 18 * idx + 14);

            auto t_zz_x = primBuffer.data(kkoff + 18 * idx + 15);

            auto t_zz_y = primBuffer.data(kkoff + 18 * idx + 16);

            auto t_zz_z = primBuffer.data(kkoff + 18 * idx + 17);

            // set up pointers to (P|D) integrals

            auto s_x_xx = primBuffer.data(s2off + 18 * idx);

            auto s_x_xy = primBuffer.data(s2off + 18 * idx + 1);

            auto s_x_xz = primBuffer.data(s2off + 18 * idx + 2);

            auto s_x_yy = primBuffer.data(s2off + 18 * idx + 3);

            auto s_x_yz = primBuffer.data(s2off + 18 * idx + 4);

            auto s_x_zz = primBuffer.data(s2off + 18 * idx + 5);

            auto s_y_xx = primBuffer.data(s2off + 18 * idx + 6);

            auto s_y_xy = primBuffer.data(s2off + 18 * idx + 7);

            auto s_y_xz = primBuffer.data(s2off + 18 * idx + 8);

            auto s_y_yy = primBuffer.data(s2off + 18 * idx + 9);

            auto s_y_yz = primBuffer.data(s2off + 18 * idx + 10);

            auto s_y_zz = primBuffer.data(s2off + 18 * idx + 11);

            auto s_z_xx = primBuffer.data(s2off + 18 * idx + 12);

            auto s_z_xy = primBuffer.data(s2off + 18 * idx + 13);

            auto s_z_xz = primBuffer.data(s2off + 18 * idx + 14);

            auto s_z_yy = primBuffer.data(s2off + 18 * idx + 15);

            auto s_z_yz = primBuffer.data(s2off + 18 * idx + 16);

            auto s_z_zz = primBuffer.data(s2off + 18 * idx + 17);

            // set up pointers to (P|T|D) integrals

            auto t_x_xx = primBuffer.data(k2off + 18 * idx);

            auto t_x_xy = primBuffer.data(k2off + 18 * idx + 1);

            auto t_x_xz = primBuffer.data(k2off + 18 * idx + 2);

            auto t_x_yy = primBuffer.data(k2off + 18 * idx + 3);

            auto t_x_yz = primBuffer.data(k2off + 18 * idx + 4);

            auto t_x_zz = primBuffer.data(k2off + 18 * idx + 5);

            auto t_y_xx = primBuffer.data(k2off + 18 * idx + 6);

            auto t_y_xy = primBuffer.data(k2off + 18 * idx + 7);

            auto t_y_xz = primBuffer.data(k2off + 18 * idx + 8);

            auto t_y_yy = primBuffer.data(k2off + 18 * idx + 9);

            auto t_y_yz = primBuffer.data(k2off + 18 * idx + 10);

            auto t_y_zz = primBuffer.data(k2off + 18 * idx + 11);

            auto t_z_xx = primBuffer.data(k2off + 18 * idx + 12);

            auto t_z_xy = primBuffer.data(k2off + 18 * idx + 13);

            auto t_z_xz = primBuffer.data(k2off + 18 * idx + 14);

            auto t_z_yy = primBuffer.data(k2off + 18 * idx + 15);

            auto t_z_yz = primBuffer.data(k2off + 18 * idx + 16);

            auto t_z_zz = primBuffer.data(k2off + 18 * idx + 17);

            // set up pointers to (D|D) integrals

            auto s_xx_xx = primBuffer.data(s1off + 36 * idx);

            auto s_xx_xy = primBuffer.data(s1off + 36 * idx + 1);

            auto s_xx_xz = primBuffer.data(s1off + 36 * idx + 2);

            auto s_xx_yy = primBuffer.data(s1off + 36 * idx + 3);

            auto s_xx_yz = primBuffer.data(s1off + 36 * idx + 4);

            auto s_xx_zz = primBuffer.data(s1off + 36 * idx + 5);

            auto s_xy_xx = primBuffer.data(s1off + 36 * idx + 6);

            auto s_xy_xy = primBuffer.data(s1off + 36 * idx + 7);

            auto s_xy_xz = primBuffer.data(s1off + 36 * idx + 8);

            auto s_xy_yy = primBuffer.data(s1off + 36 * idx + 9);

            auto s_xy_yz = primBuffer.data(s1off + 36 * idx + 10);

            auto s_xy_zz = primBuffer.data(s1off + 36 * idx + 11);

            auto s_xz_xx = primBuffer.data(s1off + 36 * idx + 12);

            auto s_xz_xy = primBuffer.data(s1off + 36 * idx + 13);

            auto s_xz_xz = primBuffer.data(s1off + 36 * idx + 14);

            auto s_xz_yy = primBuffer.data(s1off + 36 * idx + 15);

            auto s_xz_yz = primBuffer.data(s1off + 36 * idx + 16);

            auto s_xz_zz = primBuffer.data(s1off + 36 * idx + 17);

            auto s_yy_xx = primBuffer.data(s1off + 36 * idx + 18);

            auto s_yy_xy = primBuffer.data(s1off + 36 * idx + 19);

            auto s_yy_xz = primBuffer.data(s1off + 36 * idx + 20);

            auto s_yy_yy = primBuffer.data(s1off + 36 * idx + 21);

            auto s_yy_yz = primBuffer.data(s1off + 36 * idx + 22);

            auto s_yy_zz = primBuffer.data(s1off + 36 * idx + 23);

            auto s_yz_xx = primBuffer.data(s1off + 36 * idx + 24);

            auto s_yz_xy = primBuffer.data(s1off + 36 * idx + 25);

            auto s_yz_xz = primBuffer.data(s1off + 36 * idx + 26);

            auto s_yz_yy = primBuffer.data(s1off + 36 * idx + 27);

            auto s_yz_yz = primBuffer.data(s1off + 36 * idx + 28);

            auto s_yz_zz = primBuffer.data(s1off + 36 * idx + 29);

            auto s_zz_xx = primBuffer.data(s1off + 36 * idx + 30);

            auto s_zz_xy = primBuffer.data(s1off + 36 * idx + 31);

            auto s_zz_xz = primBuffer.data(s1off + 36 * idx + 32);

            auto s_zz_yy = primBuffer.data(s1off + 36 * idx + 33);

            auto s_zz_yz = primBuffer.data(s1off + 36 * idx + 34);

            auto s_zz_zz = primBuffer.data(s1off + 36 * idx + 35);

            // set up pointers to (D|T|D) integrals

            auto t_xx_xx = primBuffer.data(k1off + 36 * idx);

            auto t_xx_xy = primBuffer.data(k1off + 36 * idx + 1);

            auto t_xx_xz = primBuffer.data(k1off + 36 * idx + 2);

            auto t_xx_yy = primBuffer.data(k1off + 36 * idx + 3);

            auto t_xx_yz = primBuffer.data(k1off + 36 * idx + 4);

            auto t_xx_zz = primBuffer.data(k1off + 36 * idx + 5);

            auto t_xy_xx = primBuffer.data(k1off + 36 * idx + 6);

            auto t_xy_xy = primBuffer.data(k1off + 36 * idx + 7);

            auto t_xy_xz = primBuffer.data(k1off + 36 * idx + 8);

            auto t_xy_yy = primBuffer.data(k1off + 36 * idx + 9);

            auto t_xy_yz = primBuffer.data(k1off + 36 * idx + 10);

            auto t_xy_zz = primBuffer.data(k1off + 36 * idx + 11);

            auto t_xz_xx = primBuffer.data(k1off + 36 * idx + 12);

            auto t_xz_xy = primBuffer.data(k1off + 36 * idx + 13);

            auto t_xz_xz = primBuffer.data(k1off + 36 * idx + 14);

            auto t_xz_yy = primBuffer.data(k1off + 36 * idx + 15);

            auto t_xz_yz = primBuffer.data(k1off + 36 * idx + 16);

            auto t_xz_zz = primBuffer.data(k1off + 36 * idx + 17);

            auto t_yy_xx = primBuffer.data(k1off + 36 * idx + 18);

            auto t_yy_xy = primBuffer.data(k1off + 36 * idx + 19);

            auto t_yy_xz = primBuffer.data(k1off + 36 * idx + 20);

            auto t_yy_yy = primBuffer.data(k1off + 36 * idx + 21);

            auto t_yy_yz = primBuffer.data(k1off + 36 * idx + 22);

            auto t_yy_zz = primBuffer.data(k1off + 36 * idx + 23);

            auto t_yz_xx = primBuffer.data(k1off + 36 * idx + 24);

            auto t_yz_xy = primBuffer.data(k1off + 36 * idx + 25);

            auto t_yz_xz = primBuffer.data(k1off + 36 * idx + 26);

            auto t_yz_yy = primBuffer.data(k1off + 36 * idx + 27);

            auto t_yz_yz = primBuffer.data(k1off + 36 * idx + 28);

            auto t_yz_zz = primBuffer.data(k1off + 36 * idx + 29);

            auto t_zz_xx = primBuffer.data(k1off + 36 * idx + 30);

            auto t_zz_xy = primBuffer.data(k1off + 36 * idx + 31);

            auto t_zz_xz = primBuffer.data(k1off + 36 * idx + 32);

            auto t_zz_yy = primBuffer.data(k1off + 36 * idx + 33);

            auto t_zz_yz = primBuffer.data(k1off + 36 * idx + 34);

            auto t_zz_zz = primBuffer.data(k1off + 36 * idx + 35);

            // set up pointers to (F|D) integrals

            auto s_xxx_xx = primBuffer.data(soff + 60 * idx);

            auto s_xxx_xy = primBuffer.data(soff + 60 * idx + 1);

            auto s_xxx_xz = primBuffer.data(soff + 60 * idx + 2);

            auto s_xxx_yy = primBuffer.data(soff + 60 * idx + 3);

            auto s_xxx_yz = primBuffer.data(soff + 60 * idx + 4);

            auto s_xxx_zz = primBuffer.data(soff + 60 * idx + 5);

            auto s_xxy_xx = primBuffer.data(soff + 60 * idx + 6);

            auto s_xxy_xy = primBuffer.data(soff + 60 * idx + 7);

            auto s_xxy_xz = primBuffer.data(soff + 60 * idx + 8);

            auto s_xxy_yy = primBuffer.data(soff + 60 * idx + 9);

            auto s_xxy_yz = primBuffer.data(soff + 60 * idx + 10);

            auto s_xxy_zz = primBuffer.data(soff + 60 * idx + 11);

            auto s_xxz_xx = primBuffer.data(soff + 60 * idx + 12);

            auto s_xxz_xy = primBuffer.data(soff + 60 * idx + 13);

            auto s_xxz_xz = primBuffer.data(soff + 60 * idx + 14);

            auto s_xxz_yy = primBuffer.data(soff + 60 * idx + 15);

            auto s_xxz_yz = primBuffer.data(soff + 60 * idx + 16);

            auto s_xxz_zz = primBuffer.data(soff + 60 * idx + 17);

            auto s_xyy_xx = primBuffer.data(soff + 60 * idx + 18);

            auto s_xyy_xy = primBuffer.data(soff + 60 * idx + 19);

            auto s_xyy_xz = primBuffer.data(soff + 60 * idx + 20);

            auto s_xyy_yy = primBuffer.data(soff + 60 * idx + 21);

            auto s_xyy_yz = primBuffer.data(soff + 60 * idx + 22);

            auto s_xyy_zz = primBuffer.data(soff + 60 * idx + 23);

            auto s_xyz_xx = primBuffer.data(soff + 60 * idx + 24);

            auto s_xyz_xy = primBuffer.data(soff + 60 * idx + 25);

            auto s_xyz_xz = primBuffer.data(soff + 60 * idx + 26);

            auto s_xyz_yy = primBuffer.data(soff + 60 * idx + 27);

            auto s_xyz_yz = primBuffer.data(soff + 60 * idx + 28);

            auto s_xyz_zz = primBuffer.data(soff + 60 * idx + 29);

            auto s_xzz_xx = primBuffer.data(soff + 60 * idx + 30);

            auto s_xzz_xy = primBuffer.data(soff + 60 * idx + 31);

            auto s_xzz_xz = primBuffer.data(soff + 60 * idx + 32);

            auto s_xzz_yy = primBuffer.data(soff + 60 * idx + 33);

            auto s_xzz_yz = primBuffer.data(soff + 60 * idx + 34);

            auto s_xzz_zz = primBuffer.data(soff + 60 * idx + 35);

            auto s_yyy_xx = primBuffer.data(soff + 60 * idx + 36);

            auto s_yyy_xy = primBuffer.data(soff + 60 * idx + 37);

            auto s_yyy_xz = primBuffer.data(soff + 60 * idx + 38);

            auto s_yyy_yy = primBuffer.data(soff + 60 * idx + 39);

            auto s_yyy_yz = primBuffer.data(soff + 60 * idx + 40);

            auto s_yyy_zz = primBuffer.data(soff + 60 * idx + 41);

            auto s_yyz_xx = primBuffer.data(soff + 60 * idx + 42);

            auto s_yyz_xy = primBuffer.data(soff + 60 * idx + 43);

            auto s_yyz_xz = primBuffer.data(soff + 60 * idx + 44);

            auto s_yyz_yy = primBuffer.data(soff + 60 * idx + 45);

            auto s_yyz_yz = primBuffer.data(soff + 60 * idx + 46);

            auto s_yyz_zz = primBuffer.data(soff + 60 * idx + 47);

            auto s_yzz_xx = primBuffer.data(soff + 60 * idx + 48);

            auto s_yzz_xy = primBuffer.data(soff + 60 * idx + 49);

            auto s_yzz_xz = primBuffer.data(soff + 60 * idx + 50);

            auto s_yzz_yy = primBuffer.data(soff + 60 * idx + 51);

            auto s_yzz_yz = primBuffer.data(soff + 60 * idx + 52);

            auto s_yzz_zz = primBuffer.data(soff + 60 * idx + 53);

            auto s_zzz_xx = primBuffer.data(soff + 60 * idx + 54);

            auto s_zzz_xy = primBuffer.data(soff + 60 * idx + 55);

            auto s_zzz_xz = primBuffer.data(soff + 60 * idx + 56);

            auto s_zzz_yy = primBuffer.data(soff + 60 * idx + 57);

            auto s_zzz_yz = primBuffer.data(soff + 60 * idx + 58);

            auto s_zzz_zz = primBuffer.data(soff + 60 * idx + 59);

            // set up pointers to (F|T|D) integrals

            auto t_xxx_xx = primBuffer.data(koff + 60 * idx);

            auto t_xxx_xy = primBuffer.data(koff + 60 * idx + 1);

            auto t_xxx_xz = primBuffer.data(koff + 60 * idx + 2);

            auto t_xxx_yy = primBuffer.data(koff + 60 * idx + 3);

            auto t_xxx_yz = primBuffer.data(koff + 60 * idx + 4);

            auto t_xxx_zz = primBuffer.data(koff + 60 * idx + 5);

            auto t_xxy_xx = primBuffer.data(koff + 60 * idx + 6);

            auto t_xxy_xy = primBuffer.data(koff + 60 * idx + 7);

            auto t_xxy_xz = primBuffer.data(koff + 60 * idx + 8);

            auto t_xxy_yy = primBuffer.data(koff + 60 * idx + 9);

            auto t_xxy_yz = primBuffer.data(koff + 60 * idx + 10);

            auto t_xxy_zz = primBuffer.data(koff + 60 * idx + 11);

            auto t_xxz_xx = primBuffer.data(koff + 60 * idx + 12);

            auto t_xxz_xy = primBuffer.data(koff + 60 * idx + 13);

            auto t_xxz_xz = primBuffer.data(koff + 60 * idx + 14);

            auto t_xxz_yy = primBuffer.data(koff + 60 * idx + 15);

            auto t_xxz_yz = primBuffer.data(koff + 60 * idx + 16);

            auto t_xxz_zz = primBuffer.data(koff + 60 * idx + 17);

            auto t_xyy_xx = primBuffer.data(koff + 60 * idx + 18);

            auto t_xyy_xy = primBuffer.data(koff + 60 * idx + 19);

            auto t_xyy_xz = primBuffer.data(koff + 60 * idx + 20);

            auto t_xyy_yy = primBuffer.data(koff + 60 * idx + 21);

            auto t_xyy_yz = primBuffer.data(koff + 60 * idx + 22);

            auto t_xyy_zz = primBuffer.data(koff + 60 * idx + 23);

            auto t_xyz_xx = primBuffer.data(koff + 60 * idx + 24);

            auto t_xyz_xy = primBuffer.data(koff + 60 * idx + 25);

            auto t_xyz_xz = primBuffer.data(koff + 60 * idx + 26);

            auto t_xyz_yy = primBuffer.data(koff + 60 * idx + 27);

            auto t_xyz_yz = primBuffer.data(koff + 60 * idx + 28);

            auto t_xyz_zz = primBuffer.data(koff + 60 * idx + 29);

            auto t_xzz_xx = primBuffer.data(koff + 60 * idx + 30);

            auto t_xzz_xy = primBuffer.data(koff + 60 * idx + 31);

            auto t_xzz_xz = primBuffer.data(koff + 60 * idx + 32);

            auto t_xzz_yy = primBuffer.data(koff + 60 * idx + 33);

            auto t_xzz_yz = primBuffer.data(koff + 60 * idx + 34);

            auto t_xzz_zz = primBuffer.data(koff + 60 * idx + 35);

            auto t_yyy_xx = primBuffer.data(koff + 60 * idx + 36);

            auto t_yyy_xy = primBuffer.data(koff + 60 * idx + 37);

            auto t_yyy_xz = primBuffer.data(koff + 60 * idx + 38);

            auto t_yyy_yy = primBuffer.data(koff + 60 * idx + 39);

            auto t_yyy_yz = primBuffer.data(koff + 60 * idx + 40);

            auto t_yyy_zz = primBuffer.data(koff + 60 * idx + 41);

            auto t_yyz_xx = primBuffer.data(koff + 60 * idx + 42);

            auto t_yyz_xy = primBuffer.data(koff + 60 * idx + 43);

            auto t_yyz_xz = primBuffer.data(koff + 60 * idx + 44);

            auto t_yyz_yy = primBuffer.data(koff + 60 * idx + 45);

            auto t_yyz_yz = primBuffer.data(koff + 60 * idx + 46);

            auto t_yyz_zz = primBuffer.data(koff + 60 * idx + 47);

            auto t_yzz_xx = primBuffer.data(koff + 60 * idx + 48);

            auto t_yzz_xy = primBuffer.data(koff + 60 * idx + 49);

            auto t_yzz_xz = primBuffer.data(koff + 60 * idx + 50);

            auto t_yzz_yy = primBuffer.data(koff + 60 * idx + 51);

            auto t_yzz_yz = primBuffer.data(koff + 60 * idx + 52);

            auto t_yzz_zz = primBuffer.data(koff + 60 * idx + 53);

            auto t_zzz_xx = primBuffer.data(koff + 60 * idx + 54);

            auto t_zzz_xy = primBuffer.data(koff + 60 * idx + 55);

            auto t_zzz_xz = primBuffer.data(koff + 60 * idx + 56);

            auto t_zzz_yy = primBuffer.data(koff + 60 * idx + 57);

            auto t_zzz_yz = primBuffer.data(koff + 60 * idx + 58);

            auto t_zzz_zz = primBuffer.data(koff + 60 * idx + 59);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xx_x, s_xx_y,\
                                     s_xx_z, s_xy_x, s_xy_y, s_xy_z, s_xz_x, s_xz_y,\
                                     s_xz_z, s_yy_x, s_yy_y, s_yy_z, s_yz_x, s_yz_y,\
                                     s_yz_z, s_zz_x, s_zz_y, s_zz_z, t_xx_x, t_xx_y,\
                                     t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y,\
                                     t_xz_z, t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y,\
                                     t_yz_z, t_zz_x, t_zz_y, t_zz_z, s_x_xx, s_x_xy,\
                                     s_x_xz, s_x_yy, s_x_yz, s_x_zz, s_y_xx, s_y_xy,\
                                     s_y_xz, s_y_yy, s_y_yz, s_y_zz, s_z_xx, s_z_xy,\
                                     s_z_xz, s_z_yy, s_z_yz, s_z_zz, t_x_xx, t_x_xy,\
                                     t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy,\
                                     t_y_xz, t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy,\
                                     t_z_xz, t_z_yy, t_z_yz, t_z_zz, s_xx_xx, s_xx_xy,\
                                     s_xx_xz, s_xx_yy, s_xx_yz, s_xx_zz, s_xy_xx,\
                                     s_xy_xy, s_xy_xz, s_xy_yy, s_xy_yz, s_xy_zz,\
                                     s_xz_xx, s_xz_xy, s_xz_xz, s_xz_yy, s_xz_yz,\
                                     s_xz_zz, s_yy_xx, s_yy_xy, s_yy_xz, s_yy_yy,\
                                     s_yy_yz, s_yy_zz, s_yz_xx, s_yz_xy, s_yz_xz,\
                                     s_yz_yy, s_yz_yz, s_yz_zz, s_zz_xx, s_zz_xy,\
                                     s_zz_xz, s_zz_yy, s_zz_yz, s_zz_zz, t_xx_xx,\
                                     t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz,\
                                     t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz,\
                                     t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy,\
                                     t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy, t_yy_xz,\
                                     t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, t_yz_xy,\
                                     t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx,\
                                     t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz,\
                                     s_xxx_xx, s_xxx_xy, s_xxx_xz, s_xxx_yy, s_xxx_yz,\
                                     s_xxx_zz, s_xxy_xx, s_xxy_xy, s_xxy_xz, s_xxy_yy,\
                                     s_xxy_yz, s_xxy_zz, s_xxz_xx, s_xxz_xy, s_xxz_xz,\
                                     s_xxz_yy, s_xxz_yz, s_xxz_zz, s_xyy_xx, s_xyy_xy,\
                                     s_xyy_xz, s_xyy_yy, s_xyy_yz, s_xyy_zz, s_xyz_xx,\
                                     s_xyz_xy, s_xyz_xz, s_xyz_yy, s_xyz_yz, s_xyz_zz,\
                                     s_xzz_xx, s_xzz_xy, s_xzz_xz, s_xzz_yy, s_xzz_yz,\
                                     s_xzz_zz, s_yyy_xx, s_yyy_xy, s_yyy_xz, s_yyy_yy,\
                                     s_yyy_yz, s_yyy_zz, s_yyz_xx, s_yyz_xy, s_yyz_xz,\
                                     s_yyz_yy, s_yyz_yz, s_yyz_zz, s_yzz_xx, s_yzz_xy,\
                                     s_yzz_xz, s_yzz_yy, s_yzz_yz, s_yzz_zz, s_zzz_xx,\
                                     s_zzz_xy, s_zzz_xz, s_zzz_yy, s_zzz_yz, s_zzz_zz,\
                                     t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_yy, t_xxx_yz,\
                                     t_xxx_zz, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy,\
                                     t_xxy_yz, t_xxy_zz, t_xxz_xx, t_xxz_xy, t_xxz_xz,\
                                     t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, t_xyy_xy,\
                                     t_xyy_xz, t_xyy_yy, t_xyy_yz, t_xyy_zz, t_xyz_xx,\
                                     t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz,\
                                     t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy, t_xzz_yz,\
                                     t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy,\
                                     t_yyy_yz, t_yyy_zz, t_yyz_xx, t_yyz_xy, t_yyz_xz,\
                                     t_yyz_yy, t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy,\
                                     t_yzz_xz, t_yzz_yy, t_yzz_yz, t_yzz_zz, t_zzz_xx,\
                                     t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, t_zzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxx_xx[j] = fr * s_xx_xx[j] + f2t * (2.0 * s_x_xx[j] + 2.0 * s_xx_x[j]);

                t_xxx_xx[j] = fr * t_xx_xx[j] + f2t * (2.0 * t_x_xx[j] + 2.0 * t_xx_x[j]) + f2z * (s_xxx_xx[j] - f2a * 2.0 * s_x_xx[j]);

                s_xxx_xy[j] = fr * s_xx_xy[j] + f2t * (2.0 * s_x_xy[j] + s_xx_y[j]);

                t_xxx_xy[j] = fr * t_xx_xy[j] + f2t * (2.0 * t_x_xy[j] + t_xx_y[j]) + f2z * (s_xxx_xy[j] - f2a * 2.0 * s_x_xy[j]);

                s_xxx_xz[j] = fr * s_xx_xz[j] + f2t * (2.0 * s_x_xz[j] + s_xx_z[j]);

                t_xxx_xz[j] = fr * t_xx_xz[j] + f2t * (2.0 * t_x_xz[j] + t_xx_z[j]) + f2z * (s_xxx_xz[j] - f2a * 2.0 * s_x_xz[j]);

                s_xxx_yy[j] = fr * s_xx_yy[j] + f2t * 2.0 * s_x_yy[j];

                t_xxx_yy[j] = fr * t_xx_yy[j] + f2t * 2.0 * t_x_yy[j] + f2z * (s_xxx_yy[j] - f2a * 2.0 * s_x_yy[j]);

                s_xxx_yz[j] = fr * s_xx_yz[j] + f2t * 2.0 * s_x_yz[j];

                t_xxx_yz[j] = fr * t_xx_yz[j] + f2t * 2.0 * t_x_yz[j] + f2z * (s_xxx_yz[j] - f2a * 2.0 * s_x_yz[j]);

                s_xxx_zz[j] = fr * s_xx_zz[j] + f2t * 2.0 * s_x_zz[j];

                t_xxx_zz[j] = fr * t_xx_zz[j] + f2t * 2.0 * t_x_zz[j] + f2z * (s_xxx_zz[j] - f2a * 2.0 * s_x_zz[j]);

                s_xxy_xx[j] = fr * s_xy_xx[j] + f2t * (s_y_xx[j] + 2.0 * s_xy_x[j]);

                t_xxy_xx[j] = fr * t_xy_xx[j] + f2t * (t_y_xx[j] + 2.0 * t_xy_x[j]) + f2z * (s_xxy_xx[j] - f2a * s_y_xx[j]);

                s_xxy_xy[j] = fr * s_xy_xy[j] + f2t * (s_y_xy[j] + s_xy_y[j]);

                t_xxy_xy[j] = fr * t_xy_xy[j] + f2t * (t_y_xy[j] + t_xy_y[j]) + f2z * (s_xxy_xy[j] - f2a * s_y_xy[j]);

                s_xxy_xz[j] = fr * s_xy_xz[j] + f2t * (s_y_xz[j] + s_xy_z[j]);

                t_xxy_xz[j] = fr * t_xy_xz[j] + f2t * (t_y_xz[j] + t_xy_z[j]) + f2z * (s_xxy_xz[j] - f2a * s_y_xz[j]);

                s_xxy_yy[j] = fr * s_xy_yy[j] + f2t * s_y_yy[j];

                t_xxy_yy[j] = fr * t_xy_yy[j] + f2t * t_y_yy[j] + f2z * (s_xxy_yy[j] - f2a * s_y_yy[j]);

                s_xxy_yz[j] = fr * s_xy_yz[j] + f2t * s_y_yz[j];

                t_xxy_yz[j] = fr * t_xy_yz[j] + f2t * t_y_yz[j] + f2z * (s_xxy_yz[j] - f2a * s_y_yz[j]);

                s_xxy_zz[j] = fr * s_xy_zz[j] + f2t * s_y_zz[j];

                t_xxy_zz[j] = fr * t_xy_zz[j] + f2t * t_y_zz[j] + f2z * (s_xxy_zz[j] - f2a * s_y_zz[j]);

                s_xxz_xx[j] = fr * s_xz_xx[j] + f2t * (s_z_xx[j] + 2.0 * s_xz_x[j]);

                t_xxz_xx[j] = fr * t_xz_xx[j] + f2t * (t_z_xx[j] + 2.0 * t_xz_x[j]) + f2z * (s_xxz_xx[j] - f2a * s_z_xx[j]);

                s_xxz_xy[j] = fr * s_xz_xy[j] + f2t * (s_z_xy[j] + s_xz_y[j]);

                t_xxz_xy[j] = fr * t_xz_xy[j] + f2t * (t_z_xy[j] + t_xz_y[j]) + f2z * (s_xxz_xy[j] - f2a * s_z_xy[j]);

                s_xxz_xz[j] = fr * s_xz_xz[j] + f2t * (s_z_xz[j] + s_xz_z[j]);

                t_xxz_xz[j] = fr * t_xz_xz[j] + f2t * (t_z_xz[j] + t_xz_z[j]) + f2z * (s_xxz_xz[j] - f2a * s_z_xz[j]);

                s_xxz_yy[j] = fr * s_xz_yy[j] + f2t * s_z_yy[j];

                t_xxz_yy[j] = fr * t_xz_yy[j] + f2t * t_z_yy[j] + f2z * (s_xxz_yy[j] - f2a * s_z_yy[j]);

                s_xxz_yz[j] = fr * s_xz_yz[j] + f2t * s_z_yz[j];

                t_xxz_yz[j] = fr * t_xz_yz[j] + f2t * t_z_yz[j] + f2z * (s_xxz_yz[j] - f2a * s_z_yz[j]);

                s_xxz_zz[j] = fr * s_xz_zz[j] + f2t * s_z_zz[j];

                t_xxz_zz[j] = fr * t_xz_zz[j] + f2t * t_z_zz[j] + f2z * (s_xxz_zz[j] - f2a * s_z_zz[j]);

                s_xyy_xx[j] = fr * s_yy_xx[j] + f2t * 2.0 * s_yy_x[j];

                t_xyy_xx[j] = fr * t_yy_xx[j] + f2t * 2.0 * t_yy_x[j] + f2z * s_xyy_xx[j];

                s_xyy_xy[j] = fr * s_yy_xy[j] + f2t * s_yy_y[j];

                t_xyy_xy[j] = fr * t_yy_xy[j] + f2t * t_yy_y[j] + f2z * s_xyy_xy[j];

                s_xyy_xz[j] = fr * s_yy_xz[j] + f2t * s_yy_z[j];

                t_xyy_xz[j] = fr * t_yy_xz[j] + f2t * t_yy_z[j] + f2z * s_xyy_xz[j];

                s_xyy_yy[j] = fr * s_yy_yy[j];

                t_xyy_yy[j] = fr * t_yy_yy[j] + f2z * s_xyy_yy[j];

                s_xyy_yz[j] = fr * s_yy_yz[j];

                t_xyy_yz[j] = fr * t_yy_yz[j] + f2z * s_xyy_yz[j];

                s_xyy_zz[j] = fr * s_yy_zz[j];

                t_xyy_zz[j] = fr * t_yy_zz[j] + f2z * s_xyy_zz[j];

                s_xyz_xx[j] = fr * s_yz_xx[j] + f2t * 2.0 * s_yz_x[j];

                t_xyz_xx[j] = fr * t_yz_xx[j] + f2t * 2.0 * t_yz_x[j] + f2z * s_xyz_xx[j];

                s_xyz_xy[j] = fr * s_yz_xy[j] + f2t * s_yz_y[j];

                t_xyz_xy[j] = fr * t_yz_xy[j] + f2t * t_yz_y[j] + f2z * s_xyz_xy[j];

                s_xyz_xz[j] = fr * s_yz_xz[j] + f2t * s_yz_z[j];

                t_xyz_xz[j] = fr * t_yz_xz[j] + f2t * t_yz_z[j] + f2z * s_xyz_xz[j];

                s_xyz_yy[j] = fr * s_yz_yy[j];

                t_xyz_yy[j] = fr * t_yz_yy[j] + f2z * s_xyz_yy[j];

                s_xyz_yz[j] = fr * s_yz_yz[j];

                t_xyz_yz[j] = fr * t_yz_yz[j] + f2z * s_xyz_yz[j];

                s_xyz_zz[j] = fr * s_yz_zz[j];

                t_xyz_zz[j] = fr * t_yz_zz[j] + f2z * s_xyz_zz[j];

                s_xzz_xx[j] = fr * s_zz_xx[j] + f2t * 2.0 * s_zz_x[j];

                t_xzz_xx[j] = fr * t_zz_xx[j] + f2t * 2.0 * t_zz_x[j] + f2z * s_xzz_xx[j];

                s_xzz_xy[j] = fr * s_zz_xy[j] + f2t * s_zz_y[j];

                t_xzz_xy[j] = fr * t_zz_xy[j] + f2t * t_zz_y[j] + f2z * s_xzz_xy[j];

                s_xzz_xz[j] = fr * s_zz_xz[j] + f2t * s_zz_z[j];

                t_xzz_xz[j] = fr * t_zz_xz[j] + f2t * t_zz_z[j] + f2z * s_xzz_xz[j];

                s_xzz_yy[j] = fr * s_zz_yy[j];

                t_xzz_yy[j] = fr * t_zz_yy[j] + f2z * s_xzz_yy[j];

                s_xzz_yz[j] = fr * s_zz_yz[j];

                t_xzz_yz[j] = fr * t_zz_yz[j] + f2z * s_xzz_yz[j];

                s_xzz_zz[j] = fr * s_zz_zz[j];

                t_xzz_zz[j] = fr * t_zz_zz[j] + f2z * s_xzz_zz[j];

                // leading y component

                fr = pay[j];

                s_yyy_xx[j] = fr * s_yy_xx[j] + f2t * 2.0 * s_y_xx[j];

                t_yyy_xx[j] = fr * t_yy_xx[j] + f2t * 2.0 * t_y_xx[j] + f2z * (s_yyy_xx[j] - f2a * 2.0 * s_y_xx[j]);

                s_yyy_xy[j] = fr * s_yy_xy[j] + f2t * (2.0 * s_y_xy[j] + s_yy_x[j]);

                t_yyy_xy[j] = fr * t_yy_xy[j] + f2t * (2.0 * t_y_xy[j] + t_yy_x[j]) + f2z * (s_yyy_xy[j] - f2a * 2.0 * s_y_xy[j]);

                s_yyy_xz[j] = fr * s_yy_xz[j] + f2t * 2.0 * s_y_xz[j];

                t_yyy_xz[j] = fr * t_yy_xz[j] + f2t * 2.0 * t_y_xz[j] + f2z * (s_yyy_xz[j] - f2a * 2.0 * s_y_xz[j]);

                s_yyy_yy[j] = fr * s_yy_yy[j] + f2t * (2.0 * s_y_yy[j] + 2.0 * s_yy_y[j]);

                t_yyy_yy[j] = fr * t_yy_yy[j] + f2t * (2.0 * t_y_yy[j] + 2.0 * t_yy_y[j]) + f2z * (s_yyy_yy[j] - f2a * 2.0 * s_y_yy[j]);

                s_yyy_yz[j] = fr * s_yy_yz[j] + f2t * (2.0 * s_y_yz[j] + s_yy_z[j]);

                t_yyy_yz[j] = fr * t_yy_yz[j] + f2t * (2.0 * t_y_yz[j] + t_yy_z[j]) + f2z * (s_yyy_yz[j] - f2a * 2.0 * s_y_yz[j]);

                s_yyy_zz[j] = fr * s_yy_zz[j] + f2t * 2.0 * s_y_zz[j];

                t_yyy_zz[j] = fr * t_yy_zz[j] + f2t * 2.0 * t_y_zz[j] + f2z * (s_yyy_zz[j] - f2a * 2.0 * s_y_zz[j]);

                s_yyz_xx[j] = fr * s_yz_xx[j] + f2t * s_z_xx[j];

                t_yyz_xx[j] = fr * t_yz_xx[j] + f2t * t_z_xx[j] + f2z * (s_yyz_xx[j] - f2a * s_z_xx[j]);

                s_yyz_xy[j] = fr * s_yz_xy[j] + f2t * (s_z_xy[j] + s_yz_x[j]);

                t_yyz_xy[j] = fr * t_yz_xy[j] + f2t * (t_z_xy[j] + t_yz_x[j]) + f2z * (s_yyz_xy[j] - f2a * s_z_xy[j]);

                s_yyz_xz[j] = fr * s_yz_xz[j] + f2t * s_z_xz[j];

                t_yyz_xz[j] = fr * t_yz_xz[j] + f2t * t_z_xz[j] + f2z * (s_yyz_xz[j] - f2a * s_z_xz[j]);

                s_yyz_yy[j] = fr * s_yz_yy[j] + f2t * (s_z_yy[j] + 2.0 * s_yz_y[j]);

                t_yyz_yy[j] = fr * t_yz_yy[j] + f2t * (t_z_yy[j] + 2.0 * t_yz_y[j]) + f2z * (s_yyz_yy[j] - f2a * s_z_yy[j]);

                s_yyz_yz[j] = fr * s_yz_yz[j] + f2t * (s_z_yz[j] + s_yz_z[j]);

                t_yyz_yz[j] = fr * t_yz_yz[j] + f2t * (t_z_yz[j] + t_yz_z[j]) + f2z * (s_yyz_yz[j] - f2a * s_z_yz[j]);

                s_yyz_zz[j] = fr * s_yz_zz[j] + f2t * s_z_zz[j];

                t_yyz_zz[j] = fr * t_yz_zz[j] + f2t * t_z_zz[j] + f2z * (s_yyz_zz[j] - f2a * s_z_zz[j]);

                s_yzz_xx[j] = fr * s_zz_xx[j];

                t_yzz_xx[j] = fr * t_zz_xx[j] + f2z * s_yzz_xx[j];

                s_yzz_xy[j] = fr * s_zz_xy[j] + f2t * s_zz_x[j];

                t_yzz_xy[j] = fr * t_zz_xy[j] + f2t * t_zz_x[j] + f2z * s_yzz_xy[j];

                s_yzz_xz[j] = fr * s_zz_xz[j];

                t_yzz_xz[j] = fr * t_zz_xz[j] + f2z * s_yzz_xz[j];

                s_yzz_yy[j] = fr * s_zz_yy[j] + f2t * 2.0 * s_zz_y[j];

                t_yzz_yy[j] = fr * t_zz_yy[j] + f2t * 2.0 * t_zz_y[j] + f2z * s_yzz_yy[j];

                s_yzz_yz[j] = fr * s_zz_yz[j] + f2t * s_zz_z[j];

                t_yzz_yz[j] = fr * t_zz_yz[j] + f2t * t_zz_z[j] + f2z * s_yzz_yz[j];

                s_yzz_zz[j] = fr * s_zz_zz[j];

                t_yzz_zz[j] = fr * t_zz_zz[j] + f2z * s_yzz_zz[j];

                // leading z component

                fr = paz[j];

                s_zzz_xx[j] = fr * s_zz_xx[j] + f2t * 2.0 * s_z_xx[j];

                t_zzz_xx[j] = fr * t_zz_xx[j] + f2t * 2.0 * t_z_xx[j] + f2z * (s_zzz_xx[j] - f2a * 2.0 * s_z_xx[j]);

                s_zzz_xy[j] = fr * s_zz_xy[j] + f2t * 2.0 * s_z_xy[j];

                t_zzz_xy[j] = fr * t_zz_xy[j] + f2t * 2.0 * t_z_xy[j] + f2z * (s_zzz_xy[j] - f2a * 2.0 * s_z_xy[j]);

                s_zzz_xz[j] = fr * s_zz_xz[j] + f2t * (2.0 * s_z_xz[j] + s_zz_x[j]);

                t_zzz_xz[j] = fr * t_zz_xz[j] + f2t * (2.0 * t_z_xz[j] + t_zz_x[j]) + f2z * (s_zzz_xz[j] - f2a * 2.0 * s_z_xz[j]);

                s_zzz_yy[j] = fr * s_zz_yy[j] + f2t * 2.0 * s_z_yy[j];

                t_zzz_yy[j] = fr * t_zz_yy[j] + f2t * 2.0 * t_z_yy[j] + f2z * (s_zzz_yy[j] - f2a * 2.0 * s_z_yy[j]);

                s_zzz_yz[j] = fr * s_zz_yz[j] + f2t * (2.0 * s_z_yz[j] + s_zz_y[j]);

                t_zzz_yz[j] = fr * t_zz_yz[j] + f2t * (2.0 * t_z_yz[j] + t_zz_y[j]) + f2z * (s_zzz_yz[j] - f2a * 2.0 * s_z_yz[j]);

                s_zzz_zz[j] = fr * s_zz_zz[j] + f2t * (2.0 * s_z_zz[j] + 2.0 * s_zz_z[j]);

                t_zzz_zz[j] = fr * t_zz_zz[j] + f2t * (2.0 * t_z_zz[j] + 2.0 * t_zz_z[j]) + f2z * (s_zzz_zz[j] - f2a * 2.0 * s_z_zz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFF(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|D) integrals

            auto s_xx_xx = primBuffer.data(skoff + 36 * idx);

            auto s_xx_xy = primBuffer.data(skoff + 36 * idx + 1);

            auto s_xx_xz = primBuffer.data(skoff + 36 * idx + 2);

            auto s_xx_yy = primBuffer.data(skoff + 36 * idx + 3);

            auto s_xx_yz = primBuffer.data(skoff + 36 * idx + 4);

            auto s_xx_zz = primBuffer.data(skoff + 36 * idx + 5);

            auto s_xy_xx = primBuffer.data(skoff + 36 * idx + 6);

            auto s_xy_xy = primBuffer.data(skoff + 36 * idx + 7);

            auto s_xy_xz = primBuffer.data(skoff + 36 * idx + 8);

            auto s_xy_yy = primBuffer.data(skoff + 36 * idx + 9);

            auto s_xy_yz = primBuffer.data(skoff + 36 * idx + 10);

            auto s_xy_zz = primBuffer.data(skoff + 36 * idx + 11);

            auto s_xz_xx = primBuffer.data(skoff + 36 * idx + 12);

            auto s_xz_xy = primBuffer.data(skoff + 36 * idx + 13);

            auto s_xz_xz = primBuffer.data(skoff + 36 * idx + 14);

            auto s_xz_yy = primBuffer.data(skoff + 36 * idx + 15);

            auto s_xz_yz = primBuffer.data(skoff + 36 * idx + 16);

            auto s_xz_zz = primBuffer.data(skoff + 36 * idx + 17);

            auto s_yy_xx = primBuffer.data(skoff + 36 * idx + 18);

            auto s_yy_xy = primBuffer.data(skoff + 36 * idx + 19);

            auto s_yy_xz = primBuffer.data(skoff + 36 * idx + 20);

            auto s_yy_yy = primBuffer.data(skoff + 36 * idx + 21);

            auto s_yy_yz = primBuffer.data(skoff + 36 * idx + 22);

            auto s_yy_zz = primBuffer.data(skoff + 36 * idx + 23);

            auto s_yz_xx = primBuffer.data(skoff + 36 * idx + 24);

            auto s_yz_xy = primBuffer.data(skoff + 36 * idx + 25);

            auto s_yz_xz = primBuffer.data(skoff + 36 * idx + 26);

            auto s_yz_yy = primBuffer.data(skoff + 36 * idx + 27);

            auto s_yz_yz = primBuffer.data(skoff + 36 * idx + 28);

            auto s_yz_zz = primBuffer.data(skoff + 36 * idx + 29);

            auto s_zz_xx = primBuffer.data(skoff + 36 * idx + 30);

            auto s_zz_xy = primBuffer.data(skoff + 36 * idx + 31);

            auto s_zz_xz = primBuffer.data(skoff + 36 * idx + 32);

            auto s_zz_yy = primBuffer.data(skoff + 36 * idx + 33);

            auto s_zz_yz = primBuffer.data(skoff + 36 * idx + 34);

            auto s_zz_zz = primBuffer.data(skoff + 36 * idx + 35);

            // set up pointers to (D|T|D) integrals

            auto t_xx_xx = primBuffer.data(kkoff + 36 * idx);

            auto t_xx_xy = primBuffer.data(kkoff + 36 * idx + 1);

            auto t_xx_xz = primBuffer.data(kkoff + 36 * idx + 2);

            auto t_xx_yy = primBuffer.data(kkoff + 36 * idx + 3);

            auto t_xx_yz = primBuffer.data(kkoff + 36 * idx + 4);

            auto t_xx_zz = primBuffer.data(kkoff + 36 * idx + 5);

            auto t_xy_xx = primBuffer.data(kkoff + 36 * idx + 6);

            auto t_xy_xy = primBuffer.data(kkoff + 36 * idx + 7);

            auto t_xy_xz = primBuffer.data(kkoff + 36 * idx + 8);

            auto t_xy_yy = primBuffer.data(kkoff + 36 * idx + 9);

            auto t_xy_yz = primBuffer.data(kkoff + 36 * idx + 10);

            auto t_xy_zz = primBuffer.data(kkoff + 36 * idx + 11);

            auto t_xz_xx = primBuffer.data(kkoff + 36 * idx + 12);

            auto t_xz_xy = primBuffer.data(kkoff + 36 * idx + 13);

            auto t_xz_xz = primBuffer.data(kkoff + 36 * idx + 14);

            auto t_xz_yy = primBuffer.data(kkoff + 36 * idx + 15);

            auto t_xz_yz = primBuffer.data(kkoff + 36 * idx + 16);

            auto t_xz_zz = primBuffer.data(kkoff + 36 * idx + 17);

            auto t_yy_xx = primBuffer.data(kkoff + 36 * idx + 18);

            auto t_yy_xy = primBuffer.data(kkoff + 36 * idx + 19);

            auto t_yy_xz = primBuffer.data(kkoff + 36 * idx + 20);

            auto t_yy_yy = primBuffer.data(kkoff + 36 * idx + 21);

            auto t_yy_yz = primBuffer.data(kkoff + 36 * idx + 22);

            auto t_yy_zz = primBuffer.data(kkoff + 36 * idx + 23);

            auto t_yz_xx = primBuffer.data(kkoff + 36 * idx + 24);

            auto t_yz_xy = primBuffer.data(kkoff + 36 * idx + 25);

            auto t_yz_xz = primBuffer.data(kkoff + 36 * idx + 26);

            auto t_yz_yy = primBuffer.data(kkoff + 36 * idx + 27);

            auto t_yz_yz = primBuffer.data(kkoff + 36 * idx + 28);

            auto t_yz_zz = primBuffer.data(kkoff + 36 * idx + 29);

            auto t_zz_xx = primBuffer.data(kkoff + 36 * idx + 30);

            auto t_zz_xy = primBuffer.data(kkoff + 36 * idx + 31);

            auto t_zz_xz = primBuffer.data(kkoff + 36 * idx + 32);

            auto t_zz_yy = primBuffer.data(kkoff + 36 * idx + 33);

            auto t_zz_yz = primBuffer.data(kkoff + 36 * idx + 34);

            auto t_zz_zz = primBuffer.data(kkoff + 36 * idx + 35);

            // set up pointers to (P|F) integrals

            auto s_x_xxx = primBuffer.data(s2off + 30 * idx);

            auto s_x_xxy = primBuffer.data(s2off + 30 * idx + 1);

            auto s_x_xxz = primBuffer.data(s2off + 30 * idx + 2);

            auto s_x_xyy = primBuffer.data(s2off + 30 * idx + 3);

            auto s_x_xyz = primBuffer.data(s2off + 30 * idx + 4);

            auto s_x_xzz = primBuffer.data(s2off + 30 * idx + 5);

            auto s_x_yyy = primBuffer.data(s2off + 30 * idx + 6);

            auto s_x_yyz = primBuffer.data(s2off + 30 * idx + 7);

            auto s_x_yzz = primBuffer.data(s2off + 30 * idx + 8);

            auto s_x_zzz = primBuffer.data(s2off + 30 * idx + 9);

            auto s_y_xxx = primBuffer.data(s2off + 30 * idx + 10);

            auto s_y_xxy = primBuffer.data(s2off + 30 * idx + 11);

            auto s_y_xxz = primBuffer.data(s2off + 30 * idx + 12);

            auto s_y_xyy = primBuffer.data(s2off + 30 * idx + 13);

            auto s_y_xyz = primBuffer.data(s2off + 30 * idx + 14);

            auto s_y_xzz = primBuffer.data(s2off + 30 * idx + 15);

            auto s_y_yyy = primBuffer.data(s2off + 30 * idx + 16);

            auto s_y_yyz = primBuffer.data(s2off + 30 * idx + 17);

            auto s_y_yzz = primBuffer.data(s2off + 30 * idx + 18);

            auto s_y_zzz = primBuffer.data(s2off + 30 * idx + 19);

            auto s_z_xxx = primBuffer.data(s2off + 30 * idx + 20);

            auto s_z_xxy = primBuffer.data(s2off + 30 * idx + 21);

            auto s_z_xxz = primBuffer.data(s2off + 30 * idx + 22);

            auto s_z_xyy = primBuffer.data(s2off + 30 * idx + 23);

            auto s_z_xyz = primBuffer.data(s2off + 30 * idx + 24);

            auto s_z_xzz = primBuffer.data(s2off + 30 * idx + 25);

            auto s_z_yyy = primBuffer.data(s2off + 30 * idx + 26);

            auto s_z_yyz = primBuffer.data(s2off + 30 * idx + 27);

            auto s_z_yzz = primBuffer.data(s2off + 30 * idx + 28);

            auto s_z_zzz = primBuffer.data(s2off + 30 * idx + 29);

            // set up pointers to (P|T|F) integrals

            auto t_x_xxx = primBuffer.data(k2off + 30 * idx);

            auto t_x_xxy = primBuffer.data(k2off + 30 * idx + 1);

            auto t_x_xxz = primBuffer.data(k2off + 30 * idx + 2);

            auto t_x_xyy = primBuffer.data(k2off + 30 * idx + 3);

            auto t_x_xyz = primBuffer.data(k2off + 30 * idx + 4);

            auto t_x_xzz = primBuffer.data(k2off + 30 * idx + 5);

            auto t_x_yyy = primBuffer.data(k2off + 30 * idx + 6);

            auto t_x_yyz = primBuffer.data(k2off + 30 * idx + 7);

            auto t_x_yzz = primBuffer.data(k2off + 30 * idx + 8);

            auto t_x_zzz = primBuffer.data(k2off + 30 * idx + 9);

            auto t_y_xxx = primBuffer.data(k2off + 30 * idx + 10);

            auto t_y_xxy = primBuffer.data(k2off + 30 * idx + 11);

            auto t_y_xxz = primBuffer.data(k2off + 30 * idx + 12);

            auto t_y_xyy = primBuffer.data(k2off + 30 * idx + 13);

            auto t_y_xyz = primBuffer.data(k2off + 30 * idx + 14);

            auto t_y_xzz = primBuffer.data(k2off + 30 * idx + 15);

            auto t_y_yyy = primBuffer.data(k2off + 30 * idx + 16);

            auto t_y_yyz = primBuffer.data(k2off + 30 * idx + 17);

            auto t_y_yzz = primBuffer.data(k2off + 30 * idx + 18);

            auto t_y_zzz = primBuffer.data(k2off + 30 * idx + 19);

            auto t_z_xxx = primBuffer.data(k2off + 30 * idx + 20);

            auto t_z_xxy = primBuffer.data(k2off + 30 * idx + 21);

            auto t_z_xxz = primBuffer.data(k2off + 30 * idx + 22);

            auto t_z_xyy = primBuffer.data(k2off + 30 * idx + 23);

            auto t_z_xyz = primBuffer.data(k2off + 30 * idx + 24);

            auto t_z_xzz = primBuffer.data(k2off + 30 * idx + 25);

            auto t_z_yyy = primBuffer.data(k2off + 30 * idx + 26);

            auto t_z_yyz = primBuffer.data(k2off + 30 * idx + 27);

            auto t_z_yzz = primBuffer.data(k2off + 30 * idx + 28);

            auto t_z_zzz = primBuffer.data(k2off + 30 * idx + 29);

            // set up pointers to (D|F) integrals

            auto s_xx_xxx = primBuffer.data(s1off + 60 * idx);

            auto s_xx_xxy = primBuffer.data(s1off + 60 * idx + 1);

            auto s_xx_xxz = primBuffer.data(s1off + 60 * idx + 2);

            auto s_xx_xyy = primBuffer.data(s1off + 60 * idx + 3);

            auto s_xx_xyz = primBuffer.data(s1off + 60 * idx + 4);

            auto s_xx_xzz = primBuffer.data(s1off + 60 * idx + 5);

            auto s_xx_yyy = primBuffer.data(s1off + 60 * idx + 6);

            auto s_xx_yyz = primBuffer.data(s1off + 60 * idx + 7);

            auto s_xx_yzz = primBuffer.data(s1off + 60 * idx + 8);

            auto s_xx_zzz = primBuffer.data(s1off + 60 * idx + 9);

            auto s_xy_xxx = primBuffer.data(s1off + 60 * idx + 10);

            auto s_xy_xxy = primBuffer.data(s1off + 60 * idx + 11);

            auto s_xy_xxz = primBuffer.data(s1off + 60 * idx + 12);

            auto s_xy_xyy = primBuffer.data(s1off + 60 * idx + 13);

            auto s_xy_xyz = primBuffer.data(s1off + 60 * idx + 14);

            auto s_xy_xzz = primBuffer.data(s1off + 60 * idx + 15);

            auto s_xy_yyy = primBuffer.data(s1off + 60 * idx + 16);

            auto s_xy_yyz = primBuffer.data(s1off + 60 * idx + 17);

            auto s_xy_yzz = primBuffer.data(s1off + 60 * idx + 18);

            auto s_xy_zzz = primBuffer.data(s1off + 60 * idx + 19);

            auto s_xz_xxx = primBuffer.data(s1off + 60 * idx + 20);

            auto s_xz_xxy = primBuffer.data(s1off + 60 * idx + 21);

            auto s_xz_xxz = primBuffer.data(s1off + 60 * idx + 22);

            auto s_xz_xyy = primBuffer.data(s1off + 60 * idx + 23);

            auto s_xz_xyz = primBuffer.data(s1off + 60 * idx + 24);

            auto s_xz_xzz = primBuffer.data(s1off + 60 * idx + 25);

            auto s_xz_yyy = primBuffer.data(s1off + 60 * idx + 26);

            auto s_xz_yyz = primBuffer.data(s1off + 60 * idx + 27);

            auto s_xz_yzz = primBuffer.data(s1off + 60 * idx + 28);

            auto s_xz_zzz = primBuffer.data(s1off + 60 * idx + 29);

            auto s_yy_xxx = primBuffer.data(s1off + 60 * idx + 30);

            auto s_yy_xxy = primBuffer.data(s1off + 60 * idx + 31);

            auto s_yy_xxz = primBuffer.data(s1off + 60 * idx + 32);

            auto s_yy_xyy = primBuffer.data(s1off + 60 * idx + 33);

            auto s_yy_xyz = primBuffer.data(s1off + 60 * idx + 34);

            auto s_yy_xzz = primBuffer.data(s1off + 60 * idx + 35);

            auto s_yy_yyy = primBuffer.data(s1off + 60 * idx + 36);

            auto s_yy_yyz = primBuffer.data(s1off + 60 * idx + 37);

            auto s_yy_yzz = primBuffer.data(s1off + 60 * idx + 38);

            auto s_yy_zzz = primBuffer.data(s1off + 60 * idx + 39);

            auto s_yz_xxx = primBuffer.data(s1off + 60 * idx + 40);

            auto s_yz_xxy = primBuffer.data(s1off + 60 * idx + 41);

            auto s_yz_xxz = primBuffer.data(s1off + 60 * idx + 42);

            auto s_yz_xyy = primBuffer.data(s1off + 60 * idx + 43);

            auto s_yz_xyz = primBuffer.data(s1off + 60 * idx + 44);

            auto s_yz_xzz = primBuffer.data(s1off + 60 * idx + 45);

            auto s_yz_yyy = primBuffer.data(s1off + 60 * idx + 46);

            auto s_yz_yyz = primBuffer.data(s1off + 60 * idx + 47);

            auto s_yz_yzz = primBuffer.data(s1off + 60 * idx + 48);

            auto s_yz_zzz = primBuffer.data(s1off + 60 * idx + 49);

            auto s_zz_xxx = primBuffer.data(s1off + 60 * idx + 50);

            auto s_zz_xxy = primBuffer.data(s1off + 60 * idx + 51);

            auto s_zz_xxz = primBuffer.data(s1off + 60 * idx + 52);

            auto s_zz_xyy = primBuffer.data(s1off + 60 * idx + 53);

            auto s_zz_xyz = primBuffer.data(s1off + 60 * idx + 54);

            auto s_zz_xzz = primBuffer.data(s1off + 60 * idx + 55);

            auto s_zz_yyy = primBuffer.data(s1off + 60 * idx + 56);

            auto s_zz_yyz = primBuffer.data(s1off + 60 * idx + 57);

            auto s_zz_yzz = primBuffer.data(s1off + 60 * idx + 58);

            auto s_zz_zzz = primBuffer.data(s1off + 60 * idx + 59);

            // set up pointers to (D|T|F) integrals

            auto t_xx_xxx = primBuffer.data(k1off + 60 * idx);

            auto t_xx_xxy = primBuffer.data(k1off + 60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(k1off + 60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(k1off + 60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(k1off + 60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(k1off + 60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(k1off + 60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(k1off + 60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(k1off + 60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(k1off + 60 * idx + 9);

            auto t_xy_xxx = primBuffer.data(k1off + 60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(k1off + 60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(k1off + 60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(k1off + 60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(k1off + 60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(k1off + 60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(k1off + 60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(k1off + 60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(k1off + 60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(k1off + 60 * idx + 19);

            auto t_xz_xxx = primBuffer.data(k1off + 60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(k1off + 60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(k1off + 60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(k1off + 60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(k1off + 60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(k1off + 60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(k1off + 60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(k1off + 60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(k1off + 60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(k1off + 60 * idx + 29);

            auto t_yy_xxx = primBuffer.data(k1off + 60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(k1off + 60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(k1off + 60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(k1off + 60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(k1off + 60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(k1off + 60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(k1off + 60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(k1off + 60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(k1off + 60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(k1off + 60 * idx + 39);

            auto t_yz_xxx = primBuffer.data(k1off + 60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(k1off + 60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(k1off + 60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(k1off + 60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(k1off + 60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(k1off + 60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(k1off + 60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(k1off + 60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(k1off + 60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(k1off + 60 * idx + 49);

            auto t_zz_xxx = primBuffer.data(k1off + 60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(k1off + 60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(k1off + 60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(k1off + 60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(k1off + 60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(k1off + 60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(k1off + 60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(k1off + 60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(k1off + 60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(k1off + 60 * idx + 59);

            // set up pointers to (F|F) integrals

            auto s_xxx_xxx = primBuffer.data(soff + 100 * idx);

            auto s_xxx_xxy = primBuffer.data(soff + 100 * idx + 1);

            auto s_xxx_xxz = primBuffer.data(soff + 100 * idx + 2);

            auto s_xxx_xyy = primBuffer.data(soff + 100 * idx + 3);

            auto s_xxx_xyz = primBuffer.data(soff + 100 * idx + 4);

            auto s_xxx_xzz = primBuffer.data(soff + 100 * idx + 5);

            auto s_xxx_yyy = primBuffer.data(soff + 100 * idx + 6);

            auto s_xxx_yyz = primBuffer.data(soff + 100 * idx + 7);

            auto s_xxx_yzz = primBuffer.data(soff + 100 * idx + 8);

            auto s_xxx_zzz = primBuffer.data(soff + 100 * idx + 9);

            auto s_xxy_xxx = primBuffer.data(soff + 100 * idx + 10);

            auto s_xxy_xxy = primBuffer.data(soff + 100 * idx + 11);

            auto s_xxy_xxz = primBuffer.data(soff + 100 * idx + 12);

            auto s_xxy_xyy = primBuffer.data(soff + 100 * idx + 13);

            auto s_xxy_xyz = primBuffer.data(soff + 100 * idx + 14);

            auto s_xxy_xzz = primBuffer.data(soff + 100 * idx + 15);

            auto s_xxy_yyy = primBuffer.data(soff + 100 * idx + 16);

            auto s_xxy_yyz = primBuffer.data(soff + 100 * idx + 17);

            auto s_xxy_yzz = primBuffer.data(soff + 100 * idx + 18);

            auto s_xxy_zzz = primBuffer.data(soff + 100 * idx + 19);

            auto s_xxz_xxx = primBuffer.data(soff + 100 * idx + 20);

            auto s_xxz_xxy = primBuffer.data(soff + 100 * idx + 21);

            auto s_xxz_xxz = primBuffer.data(soff + 100 * idx + 22);

            auto s_xxz_xyy = primBuffer.data(soff + 100 * idx + 23);

            auto s_xxz_xyz = primBuffer.data(soff + 100 * idx + 24);

            auto s_xxz_xzz = primBuffer.data(soff + 100 * idx + 25);

            auto s_xxz_yyy = primBuffer.data(soff + 100 * idx + 26);

            auto s_xxz_yyz = primBuffer.data(soff + 100 * idx + 27);

            auto s_xxz_yzz = primBuffer.data(soff + 100 * idx + 28);

            auto s_xxz_zzz = primBuffer.data(soff + 100 * idx + 29);

            auto s_xyy_xxx = primBuffer.data(soff + 100 * idx + 30);

            auto s_xyy_xxy = primBuffer.data(soff + 100 * idx + 31);

            auto s_xyy_xxz = primBuffer.data(soff + 100 * idx + 32);

            auto s_xyy_xyy = primBuffer.data(soff + 100 * idx + 33);

            auto s_xyy_xyz = primBuffer.data(soff + 100 * idx + 34);

            auto s_xyy_xzz = primBuffer.data(soff + 100 * idx + 35);

            auto s_xyy_yyy = primBuffer.data(soff + 100 * idx + 36);

            auto s_xyy_yyz = primBuffer.data(soff + 100 * idx + 37);

            auto s_xyy_yzz = primBuffer.data(soff + 100 * idx + 38);

            auto s_xyy_zzz = primBuffer.data(soff + 100 * idx + 39);

            auto s_xyz_xxx = primBuffer.data(soff + 100 * idx + 40);

            auto s_xyz_xxy = primBuffer.data(soff + 100 * idx + 41);

            auto s_xyz_xxz = primBuffer.data(soff + 100 * idx + 42);

            auto s_xyz_xyy = primBuffer.data(soff + 100 * idx + 43);

            auto s_xyz_xyz = primBuffer.data(soff + 100 * idx + 44);

            auto s_xyz_xzz = primBuffer.data(soff + 100 * idx + 45);

            auto s_xyz_yyy = primBuffer.data(soff + 100 * idx + 46);

            auto s_xyz_yyz = primBuffer.data(soff + 100 * idx + 47);

            auto s_xyz_yzz = primBuffer.data(soff + 100 * idx + 48);

            auto s_xyz_zzz = primBuffer.data(soff + 100 * idx + 49);

            auto s_xzz_xxx = primBuffer.data(soff + 100 * idx + 50);

            auto s_xzz_xxy = primBuffer.data(soff + 100 * idx + 51);

            auto s_xzz_xxz = primBuffer.data(soff + 100 * idx + 52);

            auto s_xzz_xyy = primBuffer.data(soff + 100 * idx + 53);

            auto s_xzz_xyz = primBuffer.data(soff + 100 * idx + 54);

            auto s_xzz_xzz = primBuffer.data(soff + 100 * idx + 55);

            auto s_xzz_yyy = primBuffer.data(soff + 100 * idx + 56);

            auto s_xzz_yyz = primBuffer.data(soff + 100 * idx + 57);

            auto s_xzz_yzz = primBuffer.data(soff + 100 * idx + 58);

            auto s_xzz_zzz = primBuffer.data(soff + 100 * idx + 59);

            auto s_yyy_xxx = primBuffer.data(soff + 100 * idx + 60);

            auto s_yyy_xxy = primBuffer.data(soff + 100 * idx + 61);

            auto s_yyy_xxz = primBuffer.data(soff + 100 * idx + 62);

            auto s_yyy_xyy = primBuffer.data(soff + 100 * idx + 63);

            auto s_yyy_xyz = primBuffer.data(soff + 100 * idx + 64);

            auto s_yyy_xzz = primBuffer.data(soff + 100 * idx + 65);

            auto s_yyy_yyy = primBuffer.data(soff + 100 * idx + 66);

            auto s_yyy_yyz = primBuffer.data(soff + 100 * idx + 67);

            auto s_yyy_yzz = primBuffer.data(soff + 100 * idx + 68);

            auto s_yyy_zzz = primBuffer.data(soff + 100 * idx + 69);

            auto s_yyz_xxx = primBuffer.data(soff + 100 * idx + 70);

            auto s_yyz_xxy = primBuffer.data(soff + 100 * idx + 71);

            auto s_yyz_xxz = primBuffer.data(soff + 100 * idx + 72);

            auto s_yyz_xyy = primBuffer.data(soff + 100 * idx + 73);

            auto s_yyz_xyz = primBuffer.data(soff + 100 * idx + 74);

            auto s_yyz_xzz = primBuffer.data(soff + 100 * idx + 75);

            auto s_yyz_yyy = primBuffer.data(soff + 100 * idx + 76);

            auto s_yyz_yyz = primBuffer.data(soff + 100 * idx + 77);

            auto s_yyz_yzz = primBuffer.data(soff + 100 * idx + 78);

            auto s_yyz_zzz = primBuffer.data(soff + 100 * idx + 79);

            auto s_yzz_xxx = primBuffer.data(soff + 100 * idx + 80);

            auto s_yzz_xxy = primBuffer.data(soff + 100 * idx + 81);

            auto s_yzz_xxz = primBuffer.data(soff + 100 * idx + 82);

            auto s_yzz_xyy = primBuffer.data(soff + 100 * idx + 83);

            auto s_yzz_xyz = primBuffer.data(soff + 100 * idx + 84);

            auto s_yzz_xzz = primBuffer.data(soff + 100 * idx + 85);

            auto s_yzz_yyy = primBuffer.data(soff + 100 * idx + 86);

            auto s_yzz_yyz = primBuffer.data(soff + 100 * idx + 87);

            auto s_yzz_yzz = primBuffer.data(soff + 100 * idx + 88);

            auto s_yzz_zzz = primBuffer.data(soff + 100 * idx + 89);

            auto s_zzz_xxx = primBuffer.data(soff + 100 * idx + 90);

            auto s_zzz_xxy = primBuffer.data(soff + 100 * idx + 91);

            auto s_zzz_xxz = primBuffer.data(soff + 100 * idx + 92);

            auto s_zzz_xyy = primBuffer.data(soff + 100 * idx + 93);

            auto s_zzz_xyz = primBuffer.data(soff + 100 * idx + 94);

            auto s_zzz_xzz = primBuffer.data(soff + 100 * idx + 95);

            auto s_zzz_yyy = primBuffer.data(soff + 100 * idx + 96);

            auto s_zzz_yyz = primBuffer.data(soff + 100 * idx + 97);

            auto s_zzz_yzz = primBuffer.data(soff + 100 * idx + 98);

            auto s_zzz_zzz = primBuffer.data(soff + 100 * idx + 99);

            // set up pointers to (F|T|F) integrals

            auto t_xxx_xxx = primBuffer.data(koff + 100 * idx);

            auto t_xxx_xxy = primBuffer.data(koff + 100 * idx + 1);

            auto t_xxx_xxz = primBuffer.data(koff + 100 * idx + 2);

            auto t_xxx_xyy = primBuffer.data(koff + 100 * idx + 3);

            auto t_xxx_xyz = primBuffer.data(koff + 100 * idx + 4);

            auto t_xxx_xzz = primBuffer.data(koff + 100 * idx + 5);

            auto t_xxx_yyy = primBuffer.data(koff + 100 * idx + 6);

            auto t_xxx_yyz = primBuffer.data(koff + 100 * idx + 7);

            auto t_xxx_yzz = primBuffer.data(koff + 100 * idx + 8);

            auto t_xxx_zzz = primBuffer.data(koff + 100 * idx + 9);

            auto t_xxy_xxx = primBuffer.data(koff + 100 * idx + 10);

            auto t_xxy_xxy = primBuffer.data(koff + 100 * idx + 11);

            auto t_xxy_xxz = primBuffer.data(koff + 100 * idx + 12);

            auto t_xxy_xyy = primBuffer.data(koff + 100 * idx + 13);

            auto t_xxy_xyz = primBuffer.data(koff + 100 * idx + 14);

            auto t_xxy_xzz = primBuffer.data(koff + 100 * idx + 15);

            auto t_xxy_yyy = primBuffer.data(koff + 100 * idx + 16);

            auto t_xxy_yyz = primBuffer.data(koff + 100 * idx + 17);

            auto t_xxy_yzz = primBuffer.data(koff + 100 * idx + 18);

            auto t_xxy_zzz = primBuffer.data(koff + 100 * idx + 19);

            auto t_xxz_xxx = primBuffer.data(koff + 100 * idx + 20);

            auto t_xxz_xxy = primBuffer.data(koff + 100 * idx + 21);

            auto t_xxz_xxz = primBuffer.data(koff + 100 * idx + 22);

            auto t_xxz_xyy = primBuffer.data(koff + 100 * idx + 23);

            auto t_xxz_xyz = primBuffer.data(koff + 100 * idx + 24);

            auto t_xxz_xzz = primBuffer.data(koff + 100 * idx + 25);

            auto t_xxz_yyy = primBuffer.data(koff + 100 * idx + 26);

            auto t_xxz_yyz = primBuffer.data(koff + 100 * idx + 27);

            auto t_xxz_yzz = primBuffer.data(koff + 100 * idx + 28);

            auto t_xxz_zzz = primBuffer.data(koff + 100 * idx + 29);

            auto t_xyy_xxx = primBuffer.data(koff + 100 * idx + 30);

            auto t_xyy_xxy = primBuffer.data(koff + 100 * idx + 31);

            auto t_xyy_xxz = primBuffer.data(koff + 100 * idx + 32);

            auto t_xyy_xyy = primBuffer.data(koff + 100 * idx + 33);

            auto t_xyy_xyz = primBuffer.data(koff + 100 * idx + 34);

            auto t_xyy_xzz = primBuffer.data(koff + 100 * idx + 35);

            auto t_xyy_yyy = primBuffer.data(koff + 100 * idx + 36);

            auto t_xyy_yyz = primBuffer.data(koff + 100 * idx + 37);

            auto t_xyy_yzz = primBuffer.data(koff + 100 * idx + 38);

            auto t_xyy_zzz = primBuffer.data(koff + 100 * idx + 39);

            auto t_xyz_xxx = primBuffer.data(koff + 100 * idx + 40);

            auto t_xyz_xxy = primBuffer.data(koff + 100 * idx + 41);

            auto t_xyz_xxz = primBuffer.data(koff + 100 * idx + 42);

            auto t_xyz_xyy = primBuffer.data(koff + 100 * idx + 43);

            auto t_xyz_xyz = primBuffer.data(koff + 100 * idx + 44);

            auto t_xyz_xzz = primBuffer.data(koff + 100 * idx + 45);

            auto t_xyz_yyy = primBuffer.data(koff + 100 * idx + 46);

            auto t_xyz_yyz = primBuffer.data(koff + 100 * idx + 47);

            auto t_xyz_yzz = primBuffer.data(koff + 100 * idx + 48);

            auto t_xyz_zzz = primBuffer.data(koff + 100 * idx + 49);

            auto t_xzz_xxx = primBuffer.data(koff + 100 * idx + 50);

            auto t_xzz_xxy = primBuffer.data(koff + 100 * idx + 51);

            auto t_xzz_xxz = primBuffer.data(koff + 100 * idx + 52);

            auto t_xzz_xyy = primBuffer.data(koff + 100 * idx + 53);

            auto t_xzz_xyz = primBuffer.data(koff + 100 * idx + 54);

            auto t_xzz_xzz = primBuffer.data(koff + 100 * idx + 55);

            auto t_xzz_yyy = primBuffer.data(koff + 100 * idx + 56);

            auto t_xzz_yyz = primBuffer.data(koff + 100 * idx + 57);

            auto t_xzz_yzz = primBuffer.data(koff + 100 * idx + 58);

            auto t_xzz_zzz = primBuffer.data(koff + 100 * idx + 59);

            auto t_yyy_xxx = primBuffer.data(koff + 100 * idx + 60);

            auto t_yyy_xxy = primBuffer.data(koff + 100 * idx + 61);

            auto t_yyy_xxz = primBuffer.data(koff + 100 * idx + 62);

            auto t_yyy_xyy = primBuffer.data(koff + 100 * idx + 63);

            auto t_yyy_xyz = primBuffer.data(koff + 100 * idx + 64);

            auto t_yyy_xzz = primBuffer.data(koff + 100 * idx + 65);

            auto t_yyy_yyy = primBuffer.data(koff + 100 * idx + 66);

            auto t_yyy_yyz = primBuffer.data(koff + 100 * idx + 67);

            auto t_yyy_yzz = primBuffer.data(koff + 100 * idx + 68);

            auto t_yyy_zzz = primBuffer.data(koff + 100 * idx + 69);

            auto t_yyz_xxx = primBuffer.data(koff + 100 * idx + 70);

            auto t_yyz_xxy = primBuffer.data(koff + 100 * idx + 71);

            auto t_yyz_xxz = primBuffer.data(koff + 100 * idx + 72);

            auto t_yyz_xyy = primBuffer.data(koff + 100 * idx + 73);

            auto t_yyz_xyz = primBuffer.data(koff + 100 * idx + 74);

            auto t_yyz_xzz = primBuffer.data(koff + 100 * idx + 75);

            auto t_yyz_yyy = primBuffer.data(koff + 100 * idx + 76);

            auto t_yyz_yyz = primBuffer.data(koff + 100 * idx + 77);

            auto t_yyz_yzz = primBuffer.data(koff + 100 * idx + 78);

            auto t_yyz_zzz = primBuffer.data(koff + 100 * idx + 79);

            auto t_yzz_xxx = primBuffer.data(koff + 100 * idx + 80);

            auto t_yzz_xxy = primBuffer.data(koff + 100 * idx + 81);

            auto t_yzz_xxz = primBuffer.data(koff + 100 * idx + 82);

            auto t_yzz_xyy = primBuffer.data(koff + 100 * idx + 83);

            auto t_yzz_xyz = primBuffer.data(koff + 100 * idx + 84);

            auto t_yzz_xzz = primBuffer.data(koff + 100 * idx + 85);

            auto t_yzz_yyy = primBuffer.data(koff + 100 * idx + 86);

            auto t_yzz_yyz = primBuffer.data(koff + 100 * idx + 87);

            auto t_yzz_yzz = primBuffer.data(koff + 100 * idx + 88);

            auto t_yzz_zzz = primBuffer.data(koff + 100 * idx + 89);

            auto t_zzz_xxx = primBuffer.data(koff + 100 * idx + 90);

            auto t_zzz_xxy = primBuffer.data(koff + 100 * idx + 91);

            auto t_zzz_xxz = primBuffer.data(koff + 100 * idx + 92);

            auto t_zzz_xyy = primBuffer.data(koff + 100 * idx + 93);

            auto t_zzz_xyz = primBuffer.data(koff + 100 * idx + 94);

            auto t_zzz_xzz = primBuffer.data(koff + 100 * idx + 95);

            auto t_zzz_yyy = primBuffer.data(koff + 100 * idx + 96);

            auto t_zzz_yyz = primBuffer.data(koff + 100 * idx + 97);

            auto t_zzz_yzz = primBuffer.data(koff + 100 * idx + 98);

            auto t_zzz_zzz = primBuffer.data(koff + 100 * idx + 99);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xx_xx, s_xx_xy,\
                                     s_xx_xz, s_xx_yy, s_xx_yz, s_xx_zz, s_xy_xx,\
                                     s_xy_xy, s_xy_xz, s_xy_yy, s_xy_yz, s_xy_zz,\
                                     s_xz_xx, s_xz_xy, s_xz_xz, s_xz_yy, s_xz_yz,\
                                     s_xz_zz, s_yy_xx, s_yy_xy, s_yy_xz, s_yy_yy,\
                                     s_yy_yz, s_yy_zz, s_yz_xx, s_yz_xy, s_yz_xz,\
                                     s_yz_yy, s_yz_yz, s_yz_zz, s_zz_xx, s_zz_xy,\
                                     s_zz_xz, s_zz_yy, s_zz_yz, s_zz_zz, t_xx_xx,\
                                     t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz,\
                                     t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz,\
                                     t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy,\
                                     t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy, t_yy_xz,\
                                     t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, t_yz_xy,\
                                     t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx,\
                                     t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz,\
                                     s_x_xxx, s_x_xxy, s_x_xxz, s_x_xyy, s_x_xyz,\
                                     s_x_xzz, s_x_yyy, s_x_yyz, s_x_yzz, s_x_zzz,\
                                     s_y_xxx, s_y_xxy, s_y_xxz, s_y_xyy, s_y_xyz,\
                                     s_y_xzz, s_y_yyy, s_y_yyz, s_y_yzz, s_y_zzz,\
                                     s_z_xxx, s_z_xxy, s_z_xxz, s_z_xyy, s_z_xyz,\
                                     s_z_xzz, s_z_yyy, s_z_yyz, s_z_yzz, s_z_zzz,\
                                     t_x_xxx, t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz,\
                                     t_x_xzz, t_x_yyy, t_x_yyz, t_x_yzz, t_x_zzz,\
                                     t_y_xxx, t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz,\
                                     t_y_xzz, t_y_yyy, t_y_yyz, t_y_yzz, t_y_zzz,\
                                     t_z_xxx, t_z_xxy, t_z_xxz, t_z_xyy, t_z_xyz,\
                                     t_z_xzz, t_z_yyy, t_z_yyz, t_z_yzz, t_z_zzz,\
                                     s_xx_xxx, s_xx_xxy, s_xx_xxz, s_xx_xyy, s_xx_xyz,\
                                     s_xx_xzz, s_xx_yyy, s_xx_yyz, s_xx_yzz, s_xx_zzz,\
                                     s_xy_xxx, s_xy_xxy, s_xy_xxz, s_xy_xyy, s_xy_xyz,\
                                     s_xy_xzz, s_xy_yyy, s_xy_yyz, s_xy_yzz, s_xy_zzz,\
                                     s_xz_xxx, s_xz_xxy, s_xz_xxz, s_xz_xyy, s_xz_xyz,\
                                     s_xz_xzz, s_xz_yyy, s_xz_yyz, s_xz_yzz, s_xz_zzz,\
                                     s_yy_xxx, s_yy_xxy, s_yy_xxz, s_yy_xyy, s_yy_xyz,\
                                     s_yy_xzz, s_yy_yyy, s_yy_yyz, s_yy_yzz, s_yy_zzz,\
                                     s_yz_xxx, s_yz_xxy, s_yz_xxz, s_yz_xyy, s_yz_xyz,\
                                     s_yz_xzz, s_yz_yyy, s_yz_yyz, s_yz_yzz, s_yz_zzz,\
                                     s_zz_xxx, s_zz_xxy, s_zz_xxz, s_zz_xyy, s_zz_xyz,\
                                     s_zz_xzz, s_zz_yyy, s_zz_yyz, s_zz_yzz, s_zz_zzz,\
                                     t_xx_xxx, t_xx_xxy, t_xx_xxz, t_xx_xyy, t_xx_xyz,\
                                     t_xx_xzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz,\
                                     t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy, t_xy_xyz,\
                                     t_xy_xzz, t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz,\
                                     t_xz_xxx, t_xz_xxy, t_xz_xxz, t_xz_xyy, t_xz_xyz,\
                                     t_xz_xzz, t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz,\
                                     t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy, t_yy_xyz,\
                                     t_yy_xzz, t_yy_yyy, t_yy_yyz, t_yy_yzz, t_yy_zzz,\
                                     t_yz_xxx, t_yz_xxy, t_yz_xxz, t_yz_xyy, t_yz_xyz,\
                                     t_yz_xzz, t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz,\
                                     t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy, t_zz_xyz,\
                                     t_zz_xzz, t_zz_yyy, t_zz_yyz, t_zz_yzz, t_zz_zzz,\
                                     s_xxx_xxx, s_xxx_xxy, s_xxx_xxz, s_xxx_xyy,\
                                     s_xxx_xyz, s_xxx_xzz, s_xxx_yyy, s_xxx_yyz,\
                                     s_xxx_yzz, s_xxx_zzz, s_xxy_xxx, s_xxy_xxy,\
                                     s_xxy_xxz, s_xxy_xyy, s_xxy_xyz, s_xxy_xzz,\
                                     s_xxy_yyy, s_xxy_yyz, s_xxy_yzz, s_xxy_zzz,\
                                     s_xxz_xxx, s_xxz_xxy, s_xxz_xxz, s_xxz_xyy,\
                                     s_xxz_xyz, s_xxz_xzz, s_xxz_yyy, s_xxz_yyz,\
                                     s_xxz_yzz, s_xxz_zzz, s_xyy_xxx, s_xyy_xxy,\
                                     s_xyy_xxz, s_xyy_xyy, s_xyy_xyz, s_xyy_xzz,\
                                     s_xyy_yyy, s_xyy_yyz, s_xyy_yzz, s_xyy_zzz,\
                                     s_xyz_xxx, s_xyz_xxy, s_xyz_xxz, s_xyz_xyy,\
                                     s_xyz_xyz, s_xyz_xzz, s_xyz_yyy, s_xyz_yyz,\
                                     s_xyz_yzz, s_xyz_zzz, s_xzz_xxx, s_xzz_xxy,\
                                     s_xzz_xxz, s_xzz_xyy, s_xzz_xyz, s_xzz_xzz,\
                                     s_xzz_yyy, s_xzz_yyz, s_xzz_yzz, s_xzz_zzz,\
                                     s_yyy_xxx, s_yyy_xxy, s_yyy_xxz, s_yyy_xyy,\
                                     s_yyy_xyz, s_yyy_xzz, s_yyy_yyy, s_yyy_yyz,\
                                     s_yyy_yzz, s_yyy_zzz, s_yyz_xxx, s_yyz_xxy,\
                                     s_yyz_xxz, s_yyz_xyy, s_yyz_xyz, s_yyz_xzz,\
                                     s_yyz_yyy, s_yyz_yyz, s_yyz_yzz, s_yyz_zzz,\
                                     s_yzz_xxx, s_yzz_xxy, s_yzz_xxz, s_yzz_xyy,\
                                     s_yzz_xyz, s_yzz_xzz, s_yzz_yyy, s_yzz_yyz,\
                                     s_yzz_yzz, s_yzz_zzz, s_zzz_xxx, s_zzz_xxy,\
                                     s_zzz_xxz, s_zzz_xyy, s_zzz_xyz, s_zzz_xzz,\
                                     s_zzz_yyy, s_zzz_yyz, s_zzz_yzz, s_zzz_zzz,\
                                     t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, t_xxx_xyy,\
                                     t_xxx_xyz, t_xxx_xzz, t_xxx_yyy, t_xxx_yyz,\
                                     t_xxx_yzz, t_xxx_zzz, t_xxy_xxx, t_xxy_xxy,\
                                     t_xxy_xxz, t_xxy_xyy, t_xxy_xyz, t_xxy_xzz,\
                                     t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, t_xxy_zzz,\
                                     t_xxz_xxx, t_xxz_xxy, t_xxz_xxz, t_xxz_xyy,\
                                     t_xxz_xyz, t_xxz_xzz, t_xxz_yyy, t_xxz_yyz,\
                                     t_xxz_yzz, t_xxz_zzz, t_xyy_xxx, t_xyy_xxy,\
                                     t_xyy_xxz, t_xyy_xyy, t_xyy_xyz, t_xyy_xzz,\
                                     t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, t_xyy_zzz,\
                                     t_xyz_xxx, t_xyz_xxy, t_xyz_xxz, t_xyz_xyy,\
                                     t_xyz_xyz, t_xyz_xzz, t_xyz_yyy, t_xyz_yyz,\
                                     t_xyz_yzz, t_xyz_zzz, t_xzz_xxx, t_xzz_xxy,\
                                     t_xzz_xxz, t_xzz_xyy, t_xzz_xyz, t_xzz_xzz,\
                                     t_xzz_yyy, t_xzz_yyz, t_xzz_yzz, t_xzz_zzz,\
                                     t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, t_yyy_xyy,\
                                     t_yyy_xyz, t_yyy_xzz, t_yyy_yyy, t_yyy_yyz,\
                                     t_yyy_yzz, t_yyy_zzz, t_yyz_xxx, t_yyz_xxy,\
                                     t_yyz_xxz, t_yyz_xyy, t_yyz_xyz, t_yyz_xzz,\
                                     t_yyz_yyy, t_yyz_yyz, t_yyz_yzz, t_yyz_zzz,\
                                     t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy,\
                                     t_yzz_xyz, t_yzz_xzz, t_yzz_yyy, t_yzz_yyz,\
                                     t_yzz_yzz, t_yzz_zzz, t_zzz_xxx, t_zzz_xxy,\
                                     t_zzz_xxz, t_zzz_xyy, t_zzz_xyz, t_zzz_xzz,\
                                     t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, t_zzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxx_xxx[j] = fr * s_xx_xxx[j] + f2t * (2.0 * s_x_xxx[j] + 3.0 * s_xx_xx[j]);

                t_xxx_xxx[j] = fr * t_xx_xxx[j] + f2t * (2.0 * t_x_xxx[j] + 3.0 * t_xx_xx[j]) + f2z * (s_xxx_xxx[j] - f2a * 2.0 * s_x_xxx[j]);

                s_xxx_xxy[j] = fr * s_xx_xxy[j] + f2t * (2.0 * s_x_xxy[j] + 2.0 * s_xx_xy[j]);

                t_xxx_xxy[j] = fr * t_xx_xxy[j] + f2t * (2.0 * t_x_xxy[j] + 2.0 * t_xx_xy[j]) + f2z * (s_xxx_xxy[j] - f2a * 2.0 * s_x_xxy[j]);

                s_xxx_xxz[j] = fr * s_xx_xxz[j] + f2t * (2.0 * s_x_xxz[j] + 2.0 * s_xx_xz[j]);

                t_xxx_xxz[j] = fr * t_xx_xxz[j] + f2t * (2.0 * t_x_xxz[j] + 2.0 * t_xx_xz[j]) + f2z * (s_xxx_xxz[j] - f2a * 2.0 * s_x_xxz[j]);

                s_xxx_xyy[j] = fr * s_xx_xyy[j] + f2t * (2.0 * s_x_xyy[j] + s_xx_yy[j]);

                t_xxx_xyy[j] = fr * t_xx_xyy[j] + f2t * (2.0 * t_x_xyy[j] + t_xx_yy[j]) + f2z * (s_xxx_xyy[j] - f2a * 2.0 * s_x_xyy[j]);

                s_xxx_xyz[j] = fr * s_xx_xyz[j] + f2t * (2.0 * s_x_xyz[j] + s_xx_yz[j]);

                t_xxx_xyz[j] = fr * t_xx_xyz[j] + f2t * (2.0 * t_x_xyz[j] + t_xx_yz[j]) + f2z * (s_xxx_xyz[j] - f2a * 2.0 * s_x_xyz[j]);

                s_xxx_xzz[j] = fr * s_xx_xzz[j] + f2t * (2.0 * s_x_xzz[j] + s_xx_zz[j]);

                t_xxx_xzz[j] = fr * t_xx_xzz[j] + f2t * (2.0 * t_x_xzz[j] + t_xx_zz[j]) + f2z * (s_xxx_xzz[j] - f2a * 2.0 * s_x_xzz[j]);

                s_xxx_yyy[j] = fr * s_xx_yyy[j] + f2t * 2.0 * s_x_yyy[j];

                t_xxx_yyy[j] = fr * t_xx_yyy[j] + f2t * 2.0 * t_x_yyy[j] + f2z * (s_xxx_yyy[j] - f2a * 2.0 * s_x_yyy[j]);

                s_xxx_yyz[j] = fr * s_xx_yyz[j] + f2t * 2.0 * s_x_yyz[j];

                t_xxx_yyz[j] = fr * t_xx_yyz[j] + f2t * 2.0 * t_x_yyz[j] + f2z * (s_xxx_yyz[j] - f2a * 2.0 * s_x_yyz[j]);

                s_xxx_yzz[j] = fr * s_xx_yzz[j] + f2t * 2.0 * s_x_yzz[j];

                t_xxx_yzz[j] = fr * t_xx_yzz[j] + f2t * 2.0 * t_x_yzz[j] + f2z * (s_xxx_yzz[j] - f2a * 2.0 * s_x_yzz[j]);

                s_xxx_zzz[j] = fr * s_xx_zzz[j] + f2t * 2.0 * s_x_zzz[j];

                t_xxx_zzz[j] = fr * t_xx_zzz[j] + f2t * 2.0 * t_x_zzz[j] + f2z * (s_xxx_zzz[j] - f2a * 2.0 * s_x_zzz[j]);

                s_xxy_xxx[j] = fr * s_xy_xxx[j] + f2t * (s_y_xxx[j] + 3.0 * s_xy_xx[j]);

                t_xxy_xxx[j] = fr * t_xy_xxx[j] + f2t * (t_y_xxx[j] + 3.0 * t_xy_xx[j]) + f2z * (s_xxy_xxx[j] - f2a * s_y_xxx[j]);

                s_xxy_xxy[j] = fr * s_xy_xxy[j] + f2t * (s_y_xxy[j] + 2.0 * s_xy_xy[j]);

                t_xxy_xxy[j] = fr * t_xy_xxy[j] + f2t * (t_y_xxy[j] + 2.0 * t_xy_xy[j]) + f2z * (s_xxy_xxy[j] - f2a * s_y_xxy[j]);

                s_xxy_xxz[j] = fr * s_xy_xxz[j] + f2t * (s_y_xxz[j] + 2.0 * s_xy_xz[j]);

                t_xxy_xxz[j] = fr * t_xy_xxz[j] + f2t * (t_y_xxz[j] + 2.0 * t_xy_xz[j]) + f2z * (s_xxy_xxz[j] - f2a * s_y_xxz[j]);

                s_xxy_xyy[j] = fr * s_xy_xyy[j] + f2t * (s_y_xyy[j] + s_xy_yy[j]);

                t_xxy_xyy[j] = fr * t_xy_xyy[j] + f2t * (t_y_xyy[j] + t_xy_yy[j]) + f2z * (s_xxy_xyy[j] - f2a * s_y_xyy[j]);

                s_xxy_xyz[j] = fr * s_xy_xyz[j] + f2t * (s_y_xyz[j] + s_xy_yz[j]);

                t_xxy_xyz[j] = fr * t_xy_xyz[j] + f2t * (t_y_xyz[j] + t_xy_yz[j]) + f2z * (s_xxy_xyz[j] - f2a * s_y_xyz[j]);

                s_xxy_xzz[j] = fr * s_xy_xzz[j] + f2t * (s_y_xzz[j] + s_xy_zz[j]);

                t_xxy_xzz[j] = fr * t_xy_xzz[j] + f2t * (t_y_xzz[j] + t_xy_zz[j]) + f2z * (s_xxy_xzz[j] - f2a * s_y_xzz[j]);

                s_xxy_yyy[j] = fr * s_xy_yyy[j] + f2t * s_y_yyy[j];

                t_xxy_yyy[j] = fr * t_xy_yyy[j] + f2t * t_y_yyy[j] + f2z * (s_xxy_yyy[j] - f2a * s_y_yyy[j]);

                s_xxy_yyz[j] = fr * s_xy_yyz[j] + f2t * s_y_yyz[j];

                t_xxy_yyz[j] = fr * t_xy_yyz[j] + f2t * t_y_yyz[j] + f2z * (s_xxy_yyz[j] - f2a * s_y_yyz[j]);

                s_xxy_yzz[j] = fr * s_xy_yzz[j] + f2t * s_y_yzz[j];

                t_xxy_yzz[j] = fr * t_xy_yzz[j] + f2t * t_y_yzz[j] + f2z * (s_xxy_yzz[j] - f2a * s_y_yzz[j]);

                s_xxy_zzz[j] = fr * s_xy_zzz[j] + f2t * s_y_zzz[j];

                t_xxy_zzz[j] = fr * t_xy_zzz[j] + f2t * t_y_zzz[j] + f2z * (s_xxy_zzz[j] - f2a * s_y_zzz[j]);

                s_xxz_xxx[j] = fr * s_xz_xxx[j] + f2t * (s_z_xxx[j] + 3.0 * s_xz_xx[j]);

                t_xxz_xxx[j] = fr * t_xz_xxx[j] + f2t * (t_z_xxx[j] + 3.0 * t_xz_xx[j]) + f2z * (s_xxz_xxx[j] - f2a * s_z_xxx[j]);

                s_xxz_xxy[j] = fr * s_xz_xxy[j] + f2t * (s_z_xxy[j] + 2.0 * s_xz_xy[j]);

                t_xxz_xxy[j] = fr * t_xz_xxy[j] + f2t * (t_z_xxy[j] + 2.0 * t_xz_xy[j]) + f2z * (s_xxz_xxy[j] - f2a * s_z_xxy[j]);

                s_xxz_xxz[j] = fr * s_xz_xxz[j] + f2t * (s_z_xxz[j] + 2.0 * s_xz_xz[j]);

                t_xxz_xxz[j] = fr * t_xz_xxz[j] + f2t * (t_z_xxz[j] + 2.0 * t_xz_xz[j]) + f2z * (s_xxz_xxz[j] - f2a * s_z_xxz[j]);

                s_xxz_xyy[j] = fr * s_xz_xyy[j] + f2t * (s_z_xyy[j] + s_xz_yy[j]);

                t_xxz_xyy[j] = fr * t_xz_xyy[j] + f2t * (t_z_xyy[j] + t_xz_yy[j]) + f2z * (s_xxz_xyy[j] - f2a * s_z_xyy[j]);

                s_xxz_xyz[j] = fr * s_xz_xyz[j] + f2t * (s_z_xyz[j] + s_xz_yz[j]);

                t_xxz_xyz[j] = fr * t_xz_xyz[j] + f2t * (t_z_xyz[j] + t_xz_yz[j]) + f2z * (s_xxz_xyz[j] - f2a * s_z_xyz[j]);

                s_xxz_xzz[j] = fr * s_xz_xzz[j] + f2t * (s_z_xzz[j] + s_xz_zz[j]);

                t_xxz_xzz[j] = fr * t_xz_xzz[j] + f2t * (t_z_xzz[j] + t_xz_zz[j]) + f2z * (s_xxz_xzz[j] - f2a * s_z_xzz[j]);

                s_xxz_yyy[j] = fr * s_xz_yyy[j] + f2t * s_z_yyy[j];

                t_xxz_yyy[j] = fr * t_xz_yyy[j] + f2t * t_z_yyy[j] + f2z * (s_xxz_yyy[j] - f2a * s_z_yyy[j]);

                s_xxz_yyz[j] = fr * s_xz_yyz[j] + f2t * s_z_yyz[j];

                t_xxz_yyz[j] = fr * t_xz_yyz[j] + f2t * t_z_yyz[j] + f2z * (s_xxz_yyz[j] - f2a * s_z_yyz[j]);

                s_xxz_yzz[j] = fr * s_xz_yzz[j] + f2t * s_z_yzz[j];

                t_xxz_yzz[j] = fr * t_xz_yzz[j] + f2t * t_z_yzz[j] + f2z * (s_xxz_yzz[j] - f2a * s_z_yzz[j]);

                s_xxz_zzz[j] = fr * s_xz_zzz[j] + f2t * s_z_zzz[j];

                t_xxz_zzz[j] = fr * t_xz_zzz[j] + f2t * t_z_zzz[j] + f2z * (s_xxz_zzz[j] - f2a * s_z_zzz[j]);

                s_xyy_xxx[j] = fr * s_yy_xxx[j] + f2t * 3.0 * s_yy_xx[j];

                t_xyy_xxx[j] = fr * t_yy_xxx[j] + f2t * 3.0 * t_yy_xx[j] + f2z * s_xyy_xxx[j];

                s_xyy_xxy[j] = fr * s_yy_xxy[j] + f2t * 2.0 * s_yy_xy[j];

                t_xyy_xxy[j] = fr * t_yy_xxy[j] + f2t * 2.0 * t_yy_xy[j] + f2z * s_xyy_xxy[j];

                s_xyy_xxz[j] = fr * s_yy_xxz[j] + f2t * 2.0 * s_yy_xz[j];

                t_xyy_xxz[j] = fr * t_yy_xxz[j] + f2t * 2.0 * t_yy_xz[j] + f2z * s_xyy_xxz[j];

                s_xyy_xyy[j] = fr * s_yy_xyy[j] + f2t * s_yy_yy[j];

                t_xyy_xyy[j] = fr * t_yy_xyy[j] + f2t * t_yy_yy[j] + f2z * s_xyy_xyy[j];

                s_xyy_xyz[j] = fr * s_yy_xyz[j] + f2t * s_yy_yz[j];

                t_xyy_xyz[j] = fr * t_yy_xyz[j] + f2t * t_yy_yz[j] + f2z * s_xyy_xyz[j];

                s_xyy_xzz[j] = fr * s_yy_xzz[j] + f2t * s_yy_zz[j];

                t_xyy_xzz[j] = fr * t_yy_xzz[j] + f2t * t_yy_zz[j] + f2z * s_xyy_xzz[j];

                s_xyy_yyy[j] = fr * s_yy_yyy[j];

                t_xyy_yyy[j] = fr * t_yy_yyy[j] + f2z * s_xyy_yyy[j];

                s_xyy_yyz[j] = fr * s_yy_yyz[j];

                t_xyy_yyz[j] = fr * t_yy_yyz[j] + f2z * s_xyy_yyz[j];

                s_xyy_yzz[j] = fr * s_yy_yzz[j];

                t_xyy_yzz[j] = fr * t_yy_yzz[j] + f2z * s_xyy_yzz[j];

                s_xyy_zzz[j] = fr * s_yy_zzz[j];

                t_xyy_zzz[j] = fr * t_yy_zzz[j] + f2z * s_xyy_zzz[j];

                s_xyz_xxx[j] = fr * s_yz_xxx[j] + f2t * 3.0 * s_yz_xx[j];

                t_xyz_xxx[j] = fr * t_yz_xxx[j] + f2t * 3.0 * t_yz_xx[j] + f2z * s_xyz_xxx[j];

                s_xyz_xxy[j] = fr * s_yz_xxy[j] + f2t * 2.0 * s_yz_xy[j];

                t_xyz_xxy[j] = fr * t_yz_xxy[j] + f2t * 2.0 * t_yz_xy[j] + f2z * s_xyz_xxy[j];

                s_xyz_xxz[j] = fr * s_yz_xxz[j] + f2t * 2.0 * s_yz_xz[j];

                t_xyz_xxz[j] = fr * t_yz_xxz[j] + f2t * 2.0 * t_yz_xz[j] + f2z * s_xyz_xxz[j];

                s_xyz_xyy[j] = fr * s_yz_xyy[j] + f2t * s_yz_yy[j];

                t_xyz_xyy[j] = fr * t_yz_xyy[j] + f2t * t_yz_yy[j] + f2z * s_xyz_xyy[j];

                s_xyz_xyz[j] = fr * s_yz_xyz[j] + f2t * s_yz_yz[j];

                t_xyz_xyz[j] = fr * t_yz_xyz[j] + f2t * t_yz_yz[j] + f2z * s_xyz_xyz[j];

                s_xyz_xzz[j] = fr * s_yz_xzz[j] + f2t * s_yz_zz[j];

                t_xyz_xzz[j] = fr * t_yz_xzz[j] + f2t * t_yz_zz[j] + f2z * s_xyz_xzz[j];

                s_xyz_yyy[j] = fr * s_yz_yyy[j];

                t_xyz_yyy[j] = fr * t_yz_yyy[j] + f2z * s_xyz_yyy[j];

                s_xyz_yyz[j] = fr * s_yz_yyz[j];

                t_xyz_yyz[j] = fr * t_yz_yyz[j] + f2z * s_xyz_yyz[j];

                s_xyz_yzz[j] = fr * s_yz_yzz[j];

                t_xyz_yzz[j] = fr * t_yz_yzz[j] + f2z * s_xyz_yzz[j];

                s_xyz_zzz[j] = fr * s_yz_zzz[j];

                t_xyz_zzz[j] = fr * t_yz_zzz[j] + f2z * s_xyz_zzz[j];

                s_xzz_xxx[j] = fr * s_zz_xxx[j] + f2t * 3.0 * s_zz_xx[j];

                t_xzz_xxx[j] = fr * t_zz_xxx[j] + f2t * 3.0 * t_zz_xx[j] + f2z * s_xzz_xxx[j];

                s_xzz_xxy[j] = fr * s_zz_xxy[j] + f2t * 2.0 * s_zz_xy[j];

                t_xzz_xxy[j] = fr * t_zz_xxy[j] + f2t * 2.0 * t_zz_xy[j] + f2z * s_xzz_xxy[j];

                s_xzz_xxz[j] = fr * s_zz_xxz[j] + f2t * 2.0 * s_zz_xz[j];

                t_xzz_xxz[j] = fr * t_zz_xxz[j] + f2t * 2.0 * t_zz_xz[j] + f2z * s_xzz_xxz[j];

                s_xzz_xyy[j] = fr * s_zz_xyy[j] + f2t * s_zz_yy[j];

                t_xzz_xyy[j] = fr * t_zz_xyy[j] + f2t * t_zz_yy[j] + f2z * s_xzz_xyy[j];

                s_xzz_xyz[j] = fr * s_zz_xyz[j] + f2t * s_zz_yz[j];

                t_xzz_xyz[j] = fr * t_zz_xyz[j] + f2t * t_zz_yz[j] + f2z * s_xzz_xyz[j];

                s_xzz_xzz[j] = fr * s_zz_xzz[j] + f2t * s_zz_zz[j];

                t_xzz_xzz[j] = fr * t_zz_xzz[j] + f2t * t_zz_zz[j] + f2z * s_xzz_xzz[j];

                s_xzz_yyy[j] = fr * s_zz_yyy[j];

                t_xzz_yyy[j] = fr * t_zz_yyy[j] + f2z * s_xzz_yyy[j];

                s_xzz_yyz[j] = fr * s_zz_yyz[j];

                t_xzz_yyz[j] = fr * t_zz_yyz[j] + f2z * s_xzz_yyz[j];

                s_xzz_yzz[j] = fr * s_zz_yzz[j];

                t_xzz_yzz[j] = fr * t_zz_yzz[j] + f2z * s_xzz_yzz[j];

                s_xzz_zzz[j] = fr * s_zz_zzz[j];

                t_xzz_zzz[j] = fr * t_zz_zzz[j] + f2z * s_xzz_zzz[j];

                // leading y component

                fr = pay[j];

                s_yyy_xxx[j] = fr * s_yy_xxx[j] + f2t * 2.0 * s_y_xxx[j];

                t_yyy_xxx[j] = fr * t_yy_xxx[j] + f2t * 2.0 * t_y_xxx[j] + f2z * (s_yyy_xxx[j] - f2a * 2.0 * s_y_xxx[j]);

                s_yyy_xxy[j] = fr * s_yy_xxy[j] + f2t * (2.0 * s_y_xxy[j] + s_yy_xx[j]);

                t_yyy_xxy[j] = fr * t_yy_xxy[j] + f2t * (2.0 * t_y_xxy[j] + t_yy_xx[j]) + f2z * (s_yyy_xxy[j] - f2a * 2.0 * s_y_xxy[j]);

                s_yyy_xxz[j] = fr * s_yy_xxz[j] + f2t * 2.0 * s_y_xxz[j];

                t_yyy_xxz[j] = fr * t_yy_xxz[j] + f2t * 2.0 * t_y_xxz[j] + f2z * (s_yyy_xxz[j] - f2a * 2.0 * s_y_xxz[j]);

                s_yyy_xyy[j] = fr * s_yy_xyy[j] + f2t * (2.0 * s_y_xyy[j] + 2.0 * s_yy_xy[j]);

                t_yyy_xyy[j] = fr * t_yy_xyy[j] + f2t * (2.0 * t_y_xyy[j] + 2.0 * t_yy_xy[j]) + f2z * (s_yyy_xyy[j] - f2a * 2.0 * s_y_xyy[j]);

                s_yyy_xyz[j] = fr * s_yy_xyz[j] + f2t * (2.0 * s_y_xyz[j] + s_yy_xz[j]);

                t_yyy_xyz[j] = fr * t_yy_xyz[j] + f2t * (2.0 * t_y_xyz[j] + t_yy_xz[j]) + f2z * (s_yyy_xyz[j] - f2a * 2.0 * s_y_xyz[j]);

                s_yyy_xzz[j] = fr * s_yy_xzz[j] + f2t * 2.0 * s_y_xzz[j];

                t_yyy_xzz[j] = fr * t_yy_xzz[j] + f2t * 2.0 * t_y_xzz[j] + f2z * (s_yyy_xzz[j] - f2a * 2.0 * s_y_xzz[j]);

                s_yyy_yyy[j] = fr * s_yy_yyy[j] + f2t * (2.0 * s_y_yyy[j] + 3.0 * s_yy_yy[j]);

                t_yyy_yyy[j] = fr * t_yy_yyy[j] + f2t * (2.0 * t_y_yyy[j] + 3.0 * t_yy_yy[j]) + f2z * (s_yyy_yyy[j] - f2a * 2.0 * s_y_yyy[j]);

                s_yyy_yyz[j] = fr * s_yy_yyz[j] + f2t * (2.0 * s_y_yyz[j] + 2.0 * s_yy_yz[j]);

                t_yyy_yyz[j] = fr * t_yy_yyz[j] + f2t * (2.0 * t_y_yyz[j] + 2.0 * t_yy_yz[j]) + f2z * (s_yyy_yyz[j] - f2a * 2.0 * s_y_yyz[j]);

                s_yyy_yzz[j] = fr * s_yy_yzz[j] + f2t * (2.0 * s_y_yzz[j] + s_yy_zz[j]);

                t_yyy_yzz[j] = fr * t_yy_yzz[j] + f2t * (2.0 * t_y_yzz[j] + t_yy_zz[j]) + f2z * (s_yyy_yzz[j] - f2a * 2.0 * s_y_yzz[j]);

                s_yyy_zzz[j] = fr * s_yy_zzz[j] + f2t * 2.0 * s_y_zzz[j];

                t_yyy_zzz[j] = fr * t_yy_zzz[j] + f2t * 2.0 * t_y_zzz[j] + f2z * (s_yyy_zzz[j] - f2a * 2.0 * s_y_zzz[j]);

                s_yyz_xxx[j] = fr * s_yz_xxx[j] + f2t * s_z_xxx[j];

                t_yyz_xxx[j] = fr * t_yz_xxx[j] + f2t * t_z_xxx[j] + f2z * (s_yyz_xxx[j] - f2a * s_z_xxx[j]);

                s_yyz_xxy[j] = fr * s_yz_xxy[j] + f2t * (s_z_xxy[j] + s_yz_xx[j]);

                t_yyz_xxy[j] = fr * t_yz_xxy[j] + f2t * (t_z_xxy[j] + t_yz_xx[j]) + f2z * (s_yyz_xxy[j] - f2a * s_z_xxy[j]);

                s_yyz_xxz[j] = fr * s_yz_xxz[j] + f2t * s_z_xxz[j];

                t_yyz_xxz[j] = fr * t_yz_xxz[j] + f2t * t_z_xxz[j] + f2z * (s_yyz_xxz[j] - f2a * s_z_xxz[j]);

                s_yyz_xyy[j] = fr * s_yz_xyy[j] + f2t * (s_z_xyy[j] + 2.0 * s_yz_xy[j]);

                t_yyz_xyy[j] = fr * t_yz_xyy[j] + f2t * (t_z_xyy[j] + 2.0 * t_yz_xy[j]) + f2z * (s_yyz_xyy[j] - f2a * s_z_xyy[j]);

                s_yyz_xyz[j] = fr * s_yz_xyz[j] + f2t * (s_z_xyz[j] + s_yz_xz[j]);

                t_yyz_xyz[j] = fr * t_yz_xyz[j] + f2t * (t_z_xyz[j] + t_yz_xz[j]) + f2z * (s_yyz_xyz[j] - f2a * s_z_xyz[j]);

                s_yyz_xzz[j] = fr * s_yz_xzz[j] + f2t * s_z_xzz[j];

                t_yyz_xzz[j] = fr * t_yz_xzz[j] + f2t * t_z_xzz[j] + f2z * (s_yyz_xzz[j] - f2a * s_z_xzz[j]);

                s_yyz_yyy[j] = fr * s_yz_yyy[j] + f2t * (s_z_yyy[j] + 3.0 * s_yz_yy[j]);

                t_yyz_yyy[j] = fr * t_yz_yyy[j] + f2t * (t_z_yyy[j] + 3.0 * t_yz_yy[j]) + f2z * (s_yyz_yyy[j] - f2a * s_z_yyy[j]);

                s_yyz_yyz[j] = fr * s_yz_yyz[j] + f2t * (s_z_yyz[j] + 2.0 * s_yz_yz[j]);

                t_yyz_yyz[j] = fr * t_yz_yyz[j] + f2t * (t_z_yyz[j] + 2.0 * t_yz_yz[j]) + f2z * (s_yyz_yyz[j] - f2a * s_z_yyz[j]);

                s_yyz_yzz[j] = fr * s_yz_yzz[j] + f2t * (s_z_yzz[j] + s_yz_zz[j]);

                t_yyz_yzz[j] = fr * t_yz_yzz[j] + f2t * (t_z_yzz[j] + t_yz_zz[j]) + f2z * (s_yyz_yzz[j] - f2a * s_z_yzz[j]);

                s_yyz_zzz[j] = fr * s_yz_zzz[j] + f2t * s_z_zzz[j];

                t_yyz_zzz[j] = fr * t_yz_zzz[j] + f2t * t_z_zzz[j] + f2z * (s_yyz_zzz[j] - f2a * s_z_zzz[j]);

                s_yzz_xxx[j] = fr * s_zz_xxx[j];

                t_yzz_xxx[j] = fr * t_zz_xxx[j] + f2z * s_yzz_xxx[j];

                s_yzz_xxy[j] = fr * s_zz_xxy[j] + f2t * s_zz_xx[j];

                t_yzz_xxy[j] = fr * t_zz_xxy[j] + f2t * t_zz_xx[j] + f2z * s_yzz_xxy[j];

                s_yzz_xxz[j] = fr * s_zz_xxz[j];

                t_yzz_xxz[j] = fr * t_zz_xxz[j] + f2z * s_yzz_xxz[j];

                s_yzz_xyy[j] = fr * s_zz_xyy[j] + f2t * 2.0 * s_zz_xy[j];

                t_yzz_xyy[j] = fr * t_zz_xyy[j] + f2t * 2.0 * t_zz_xy[j] + f2z * s_yzz_xyy[j];

                s_yzz_xyz[j] = fr * s_zz_xyz[j] + f2t * s_zz_xz[j];

                t_yzz_xyz[j] = fr * t_zz_xyz[j] + f2t * t_zz_xz[j] + f2z * s_yzz_xyz[j];

                s_yzz_xzz[j] = fr * s_zz_xzz[j];

                t_yzz_xzz[j] = fr * t_zz_xzz[j] + f2z * s_yzz_xzz[j];

                s_yzz_yyy[j] = fr * s_zz_yyy[j] + f2t * 3.0 * s_zz_yy[j];

                t_yzz_yyy[j] = fr * t_zz_yyy[j] + f2t * 3.0 * t_zz_yy[j] + f2z * s_yzz_yyy[j];

                s_yzz_yyz[j] = fr * s_zz_yyz[j] + f2t * 2.0 * s_zz_yz[j];

                t_yzz_yyz[j] = fr * t_zz_yyz[j] + f2t * 2.0 * t_zz_yz[j] + f2z * s_yzz_yyz[j];

                s_yzz_yzz[j] = fr * s_zz_yzz[j] + f2t * s_zz_zz[j];

                t_yzz_yzz[j] = fr * t_zz_yzz[j] + f2t * t_zz_zz[j] + f2z * s_yzz_yzz[j];

                s_yzz_zzz[j] = fr * s_zz_zzz[j];

                t_yzz_zzz[j] = fr * t_zz_zzz[j] + f2z * s_yzz_zzz[j];

                // leading z component

                fr = paz[j];

                s_zzz_xxx[j] = fr * s_zz_xxx[j] + f2t * 2.0 * s_z_xxx[j];

                t_zzz_xxx[j] = fr * t_zz_xxx[j] + f2t * 2.0 * t_z_xxx[j] + f2z * (s_zzz_xxx[j] - f2a * 2.0 * s_z_xxx[j]);

                s_zzz_xxy[j] = fr * s_zz_xxy[j] + f2t * 2.0 * s_z_xxy[j];

                t_zzz_xxy[j] = fr * t_zz_xxy[j] + f2t * 2.0 * t_z_xxy[j] + f2z * (s_zzz_xxy[j] - f2a * 2.0 * s_z_xxy[j]);

                s_zzz_xxz[j] = fr * s_zz_xxz[j] + f2t * (2.0 * s_z_xxz[j] + s_zz_xx[j]);

                t_zzz_xxz[j] = fr * t_zz_xxz[j] + f2t * (2.0 * t_z_xxz[j] + t_zz_xx[j]) + f2z * (s_zzz_xxz[j] - f2a * 2.0 * s_z_xxz[j]);

                s_zzz_xyy[j] = fr * s_zz_xyy[j] + f2t * 2.0 * s_z_xyy[j];

                t_zzz_xyy[j] = fr * t_zz_xyy[j] + f2t * 2.0 * t_z_xyy[j] + f2z * (s_zzz_xyy[j] - f2a * 2.0 * s_z_xyy[j]);

                s_zzz_xyz[j] = fr * s_zz_xyz[j] + f2t * (2.0 * s_z_xyz[j] + s_zz_xy[j]);

                t_zzz_xyz[j] = fr * t_zz_xyz[j] + f2t * (2.0 * t_z_xyz[j] + t_zz_xy[j]) + f2z * (s_zzz_xyz[j] - f2a * 2.0 * s_z_xyz[j]);

                s_zzz_xzz[j] = fr * s_zz_xzz[j] + f2t * (2.0 * s_z_xzz[j] + 2.0 * s_zz_xz[j]);

                t_zzz_xzz[j] = fr * t_zz_xzz[j] + f2t * (2.0 * t_z_xzz[j] + 2.0 * t_zz_xz[j]) + f2z * (s_zzz_xzz[j] - f2a * 2.0 * s_z_xzz[j]);

                s_zzz_yyy[j] = fr * s_zz_yyy[j] + f2t * 2.0 * s_z_yyy[j];

                t_zzz_yyy[j] = fr * t_zz_yyy[j] + f2t * 2.0 * t_z_yyy[j] + f2z * (s_zzz_yyy[j] - f2a * 2.0 * s_z_yyy[j]);

                s_zzz_yyz[j] = fr * s_zz_yyz[j] + f2t * (2.0 * s_z_yyz[j] + s_zz_yy[j]);

                t_zzz_yyz[j] = fr * t_zz_yyz[j] + f2t * (2.0 * t_z_yyz[j] + t_zz_yy[j]) + f2z * (s_zzz_yyz[j] - f2a * 2.0 * s_z_yyz[j]);

                s_zzz_yzz[j] = fr * s_zz_yzz[j] + f2t * (2.0 * s_z_yzz[j] + 2.0 * s_zz_yz[j]);

                t_zzz_yzz[j] = fr * t_zz_yzz[j] + f2t * (2.0 * t_z_yzz[j] + 2.0 * t_zz_yz[j]) + f2z * (s_zzz_yzz[j] - f2a * 2.0 * s_z_yzz[j]);

                s_zzz_zzz[j] = fr * s_zz_zzz[j] + f2t * (2.0 * s_z_zzz[j] + 3.0 * s_zz_zz[j]);

                t_zzz_zzz[j] = fr * t_zz_zzz[j] + f2t * (2.0 * t_z_zzz[j] + 3.0 * t_zz_zz[j]) + f2z * (s_zzz_zzz[j] - f2a * 2.0 * s_z_zzz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForSG(      CMemBlock2D<double>&  primBuffer,
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

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|D) integrals

            auto s_0_xx = primBuffer.data(s2off + 6 * idx);

            auto s_0_xy = primBuffer.data(s2off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(s2off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(s2off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(s2off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(s2off + 6 * idx + 5);

            // set up pointers to (S|T|D) integrals

            auto t_0_xx = primBuffer.data(k2off + 6 * idx);

            auto t_0_xy = primBuffer.data(k2off + 6 * idx + 1);

            auto t_0_xz = primBuffer.data(k2off + 6 * idx + 2);

            auto t_0_yy = primBuffer.data(k2off + 6 * idx + 3);

            auto t_0_yz = primBuffer.data(k2off + 6 * idx + 4);

            auto t_0_zz = primBuffer.data(k2off + 6 * idx + 5);

            // set up pointers to (S|F) integrals

            auto s_0_xxx = primBuffer.data(s1off + 10 * idx);

            auto s_0_xxy = primBuffer.data(s1off + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(s1off + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(s1off + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(s1off + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(s1off + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(s1off + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(s1off + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(s1off + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(s1off + 10 * idx + 9);

            // set up pointers to (S|T|F) integrals

            auto t_0_xxx = primBuffer.data(k1off + 10 * idx);

            auto t_0_xxy = primBuffer.data(k1off + 10 * idx + 1);

            auto t_0_xxz = primBuffer.data(k1off + 10 * idx + 2);

            auto t_0_xyy = primBuffer.data(k1off + 10 * idx + 3);

            auto t_0_xyz = primBuffer.data(k1off + 10 * idx + 4);

            auto t_0_xzz = primBuffer.data(k1off + 10 * idx + 5);

            auto t_0_yyy = primBuffer.data(k1off + 10 * idx + 6);

            auto t_0_yyz = primBuffer.data(k1off + 10 * idx + 7);

            auto t_0_yzz = primBuffer.data(k1off + 10 * idx + 8);

            auto t_0_zzz = primBuffer.data(k1off + 10 * idx + 9);

            // set up pointers to (S|G) integrals

            auto s_0_xxxx = primBuffer.data(soff + 15 * idx);

            auto s_0_xxxy = primBuffer.data(soff + 15 * idx + 1);

            auto s_0_xxxz = primBuffer.data(soff + 15 * idx + 2);

            auto s_0_xxyy = primBuffer.data(soff + 15 * idx + 3);

            auto s_0_xxyz = primBuffer.data(soff + 15 * idx + 4);

            auto s_0_xxzz = primBuffer.data(soff + 15 * idx + 5);

            auto s_0_xyyy = primBuffer.data(soff + 15 * idx + 6);

            auto s_0_xyyz = primBuffer.data(soff + 15 * idx + 7);

            auto s_0_xyzz = primBuffer.data(soff + 15 * idx + 8);

            auto s_0_xzzz = primBuffer.data(soff + 15 * idx + 9);

            auto s_0_yyyy = primBuffer.data(soff + 15 * idx + 10);

            auto s_0_yyyz = primBuffer.data(soff + 15 * idx + 11);

            auto s_0_yyzz = primBuffer.data(soff + 15 * idx + 12);

            auto s_0_yzzz = primBuffer.data(soff + 15 * idx + 13);

            auto s_0_zzzz = primBuffer.data(soff + 15 * idx + 14);

            // set up pointers to (S|T|G) integrals

            auto t_0_xxxx = primBuffer.data(koff + 15 * idx);

            auto t_0_xxxy = primBuffer.data(koff + 15 * idx + 1);

            auto t_0_xxxz = primBuffer.data(koff + 15 * idx + 2);

            auto t_0_xxyy = primBuffer.data(koff + 15 * idx + 3);

            auto t_0_xxyz = primBuffer.data(koff + 15 * idx + 4);

            auto t_0_xxzz = primBuffer.data(koff + 15 * idx + 5);

            auto t_0_xyyy = primBuffer.data(koff + 15 * idx + 6);

            auto t_0_xyyz = primBuffer.data(koff + 15 * idx + 7);

            auto t_0_xyzz = primBuffer.data(koff + 15 * idx + 8);

            auto t_0_xzzz = primBuffer.data(koff + 15 * idx + 9);

            auto t_0_yyyy = primBuffer.data(koff + 15 * idx + 10);

            auto t_0_yyyz = primBuffer.data(koff + 15 * idx + 11);

            auto t_0_yyzz = primBuffer.data(koff + 15 * idx + 12);

            auto t_0_yzzz = primBuffer.data(koff + 15 * idx + 13);

            auto t_0_zzzz = primBuffer.data(koff + 15 * idx + 14);

            #pragma omp simd aligned(fx, fz, fgb, pbx, pby, pbz, s_0_xx, s_0_xy,\
                                     s_0_xz, s_0_yy, s_0_yz, s_0_zz, t_0_xx, t_0_xy,\
                                     t_0_xz, t_0_yy, t_0_yz, t_0_zz, s_0_xxx, s_0_xxy,\
                                     s_0_xxz, s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz, t_0_xxx, t_0_xxy,\
                                     t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy,\
                                     t_0_yyz, t_0_yzz, t_0_zzz, s_0_xxxx, s_0_xxxy,\
                                     s_0_xxxz, s_0_xxyy, s_0_xxyz, s_0_xxzz, s_0_xyyy,\
                                     s_0_xyyz, s_0_xyzz, s_0_xzzz, s_0_yyyy, s_0_yyyz,\
                                     s_0_yyzz, s_0_yzzz, s_0_zzzz, t_0_xxxx, t_0_xxxy,\
                                     t_0_xxxz, t_0_xxyy, t_0_xxyz, t_0_xxzz, t_0_xyyy,\
                                     t_0_xyyz, t_0_xyzz, t_0_xzzz, t_0_yyyy, t_0_yyyz,\
                                     t_0_yyzz, t_0_yzzz, t_0_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2b = 0.50 * fgb[j];

                // leading x component

                double fr = pbx[j];

                s_0_xxxx[j] = fr * s_0_xxx[j] + f2t * 3.0 * s_0_xx[j];

                t_0_xxxx[j] = fr * t_0_xxx[j] + f2t * 3.0 * t_0_xx[j] + f2z * (s_0_xxxx[j] - f2b * 3.0 * s_0_xx[j]);

                s_0_xxxy[j] = fr * s_0_xxy[j] + f2t * 2.0 * s_0_xy[j];

                t_0_xxxy[j] = fr * t_0_xxy[j] + f2t * 2.0 * t_0_xy[j] + f2z * (s_0_xxxy[j] - f2b * 2.0 * s_0_xy[j]);

                s_0_xxxz[j] = fr * s_0_xxz[j] + f2t * 2.0 * s_0_xz[j];

                t_0_xxxz[j] = fr * t_0_xxz[j] + f2t * 2.0 * t_0_xz[j] + f2z * (s_0_xxxz[j] - f2b * 2.0 * s_0_xz[j]);

                s_0_xxyy[j] = fr * s_0_xyy[j] + f2t * s_0_yy[j];

                t_0_xxyy[j] = fr * t_0_xyy[j] + f2t * t_0_yy[j] + f2z * (s_0_xxyy[j] - f2b * s_0_yy[j]);

                s_0_xxyz[j] = fr * s_0_xyz[j] + f2t * s_0_yz[j];

                t_0_xxyz[j] = fr * t_0_xyz[j] + f2t * t_0_yz[j] + f2z * (s_0_xxyz[j] - f2b * s_0_yz[j]);

                s_0_xxzz[j] = fr * s_0_xzz[j] + f2t * s_0_zz[j];

                t_0_xxzz[j] = fr * t_0_xzz[j] + f2t * t_0_zz[j] + f2z * (s_0_xxzz[j] - f2b * s_0_zz[j]);

                s_0_xyyy[j] = fr * s_0_yyy[j];

                t_0_xyyy[j] = fr * t_0_yyy[j] + f2z * s_0_xyyy[j];

                s_0_xyyz[j] = fr * s_0_yyz[j];

                t_0_xyyz[j] = fr * t_0_yyz[j] + f2z * s_0_xyyz[j];

                s_0_xyzz[j] = fr * s_0_yzz[j];

                t_0_xyzz[j] = fr * t_0_yzz[j] + f2z * s_0_xyzz[j];

                s_0_xzzz[j] = fr * s_0_zzz[j];

                t_0_xzzz[j] = fr * t_0_zzz[j] + f2z * s_0_xzzz[j];

                // leading y component

                fr = pby[j];

                s_0_yyyy[j] = fr * s_0_yyy[j] + f2t * 3.0 * s_0_yy[j];

                t_0_yyyy[j] = fr * t_0_yyy[j] + f2t * 3.0 * t_0_yy[j] + f2z * (s_0_yyyy[j] - f2b * 3.0 * s_0_yy[j]);

                s_0_yyyz[j] = fr * s_0_yyz[j] + f2t * 2.0 * s_0_yz[j];

                t_0_yyyz[j] = fr * t_0_yyz[j] + f2t * 2.0 * t_0_yz[j] + f2z * (s_0_yyyz[j] - f2b * 2.0 * s_0_yz[j]);

                s_0_yyzz[j] = fr * s_0_yzz[j] + f2t * s_0_zz[j];

                t_0_yyzz[j] = fr * t_0_yzz[j] + f2t * t_0_zz[j] + f2z * (s_0_yyzz[j] - f2b * s_0_zz[j]);

                s_0_yzzz[j] = fr * s_0_zzz[j];

                t_0_yzzz[j] = fr * t_0_zzz[j] + f2z * s_0_yzzz[j];

                // leading z component

                fr = pbz[j];

                s_0_zzzz[j] = fr * s_0_zzz[j] + f2t * 3.0 * s_0_zz[j];

                t_0_zzzz[j] = fr * t_0_zzz[j] + f2t * 3.0 * t_0_zz[j] + f2z * (s_0_zzzz[j] - f2b * 3.0 * s_0_zz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGS(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 0, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 0, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|S) integrals

            auto s_xx_0 = primBuffer.data(s2off + 6 * idx);

            auto s_xy_0 = primBuffer.data(s2off + 6 * idx + 1);

            auto s_xz_0 = primBuffer.data(s2off + 6 * idx + 2);

            auto s_yy_0 = primBuffer.data(s2off + 6 * idx + 3);

            auto s_yz_0 = primBuffer.data(s2off + 6 * idx + 4);

            auto s_zz_0 = primBuffer.data(s2off + 6 * idx + 5);

            // set up pointers to (D|T|S) integrals

            auto t_xx_0 = primBuffer.data(k2off + 6 * idx);

            auto t_xy_0 = primBuffer.data(k2off + 6 * idx + 1);

            auto t_xz_0 = primBuffer.data(k2off + 6 * idx + 2);

            auto t_yy_0 = primBuffer.data(k2off + 6 * idx + 3);

            auto t_yz_0 = primBuffer.data(k2off + 6 * idx + 4);

            auto t_zz_0 = primBuffer.data(k2off + 6 * idx + 5);

            // set up pointers to (F|S) integrals

            auto s_xxx_0 = primBuffer.data(s1off + 10 * idx);

            auto s_xxy_0 = primBuffer.data(s1off + 10 * idx + 1);

            auto s_xxz_0 = primBuffer.data(s1off + 10 * idx + 2);

            auto s_xyy_0 = primBuffer.data(s1off + 10 * idx + 3);

            auto s_xyz_0 = primBuffer.data(s1off + 10 * idx + 4);

            auto s_xzz_0 = primBuffer.data(s1off + 10 * idx + 5);

            auto s_yyy_0 = primBuffer.data(s1off + 10 * idx + 6);

            auto s_yyz_0 = primBuffer.data(s1off + 10 * idx + 7);

            auto s_yzz_0 = primBuffer.data(s1off + 10 * idx + 8);

            auto s_zzz_0 = primBuffer.data(s1off + 10 * idx + 9);

            // set up pointers to (F|T|S) integrals

            auto t_xxx_0 = primBuffer.data(k1off + 10 * idx);

            auto t_xxy_0 = primBuffer.data(k1off + 10 * idx + 1);

            auto t_xxz_0 = primBuffer.data(k1off + 10 * idx + 2);

            auto t_xyy_0 = primBuffer.data(k1off + 10 * idx + 3);

            auto t_xyz_0 = primBuffer.data(k1off + 10 * idx + 4);

            auto t_xzz_0 = primBuffer.data(k1off + 10 * idx + 5);

            auto t_yyy_0 = primBuffer.data(k1off + 10 * idx + 6);

            auto t_yyz_0 = primBuffer.data(k1off + 10 * idx + 7);

            auto t_yzz_0 = primBuffer.data(k1off + 10 * idx + 8);

            auto t_zzz_0 = primBuffer.data(k1off + 10 * idx + 9);

            // set up pointers to (G|S) integrals

            auto s_xxxx_0 = primBuffer.data(soff + 15 * idx);

            auto s_xxxy_0 = primBuffer.data(soff + 15 * idx + 1);

            auto s_xxxz_0 = primBuffer.data(soff + 15 * idx + 2);

            auto s_xxyy_0 = primBuffer.data(soff + 15 * idx + 3);

            auto s_xxyz_0 = primBuffer.data(soff + 15 * idx + 4);

            auto s_xxzz_0 = primBuffer.data(soff + 15 * idx + 5);

            auto s_xyyy_0 = primBuffer.data(soff + 15 * idx + 6);

            auto s_xyyz_0 = primBuffer.data(soff + 15 * idx + 7);

            auto s_xyzz_0 = primBuffer.data(soff + 15 * idx + 8);

            auto s_xzzz_0 = primBuffer.data(soff + 15 * idx + 9);

            auto s_yyyy_0 = primBuffer.data(soff + 15 * idx + 10);

            auto s_yyyz_0 = primBuffer.data(soff + 15 * idx + 11);

            auto s_yyzz_0 = primBuffer.data(soff + 15 * idx + 12);

            auto s_yzzz_0 = primBuffer.data(soff + 15 * idx + 13);

            auto s_zzzz_0 = primBuffer.data(soff + 15 * idx + 14);

            // set up pointers to (G|T|S) integrals

            auto t_xxxx_0 = primBuffer.data(koff + 15 * idx);

            auto t_xxxy_0 = primBuffer.data(koff + 15 * idx + 1);

            auto t_xxxz_0 = primBuffer.data(koff + 15 * idx + 2);

            auto t_xxyy_0 = primBuffer.data(koff + 15 * idx + 3);

            auto t_xxyz_0 = primBuffer.data(koff + 15 * idx + 4);

            auto t_xxzz_0 = primBuffer.data(koff + 15 * idx + 5);

            auto t_xyyy_0 = primBuffer.data(koff + 15 * idx + 6);

            auto t_xyyz_0 = primBuffer.data(koff + 15 * idx + 7);

            auto t_xyzz_0 = primBuffer.data(koff + 15 * idx + 8);

            auto t_xzzz_0 = primBuffer.data(koff + 15 * idx + 9);

            auto t_yyyy_0 = primBuffer.data(koff + 15 * idx + 10);

            auto t_yyyz_0 = primBuffer.data(koff + 15 * idx + 11);

            auto t_yyzz_0 = primBuffer.data(koff + 15 * idx + 12);

            auto t_yzzz_0 = primBuffer.data(koff + 15 * idx + 13);

            auto t_zzzz_0 = primBuffer.data(koff + 15 * idx + 14);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xx_0, s_xy_0,\
                                     s_xz_0, s_yy_0, s_yz_0, s_zz_0, t_xx_0, t_xy_0,\
                                     t_xz_0, t_yy_0, t_yz_0, t_zz_0, s_xxx_0, s_xxy_0,\
                                     s_xxz_0, s_xyy_0, s_xyz_0, s_xzz_0, s_yyy_0,\
                                     s_yyz_0, s_yzz_0, s_zzz_0, t_xxx_0, t_xxy_0,\
                                     t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, t_yyy_0,\
                                     t_yyz_0, t_yzz_0, t_zzz_0, s_xxxx_0, s_xxxy_0,\
                                     s_xxxz_0, s_xxyy_0, s_xxyz_0, s_xxzz_0, s_xyyy_0,\
                                     s_xyyz_0, s_xyzz_0, s_xzzz_0, s_yyyy_0, s_yyyz_0,\
                                     s_yyzz_0, s_yzzz_0, s_zzzz_0, t_xxxx_0, t_xxxy_0,\
                                     t_xxxz_0, t_xxyy_0, t_xxyz_0, t_xxzz_0, t_xyyy_0,\
                                     t_xyyz_0, t_xyzz_0, t_xzzz_0, t_yyyy_0, t_yyyz_0,\
                                     t_yyzz_0, t_yzzz_0, t_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxxx_0[j] = fr * s_xxx_0[j] + f2t * 3.0 * s_xx_0[j];

                t_xxxx_0[j] = fr * t_xxx_0[j] + f2t * 3.0 * t_xx_0[j] + f2z * (s_xxxx_0[j] - f2a * 3.0 * s_xx_0[j]);

                s_xxxy_0[j] = fr * s_xxy_0[j] + f2t * 2.0 * s_xy_0[j];

                t_xxxy_0[j] = fr * t_xxy_0[j] + f2t * 2.0 * t_xy_0[j] + f2z * (s_xxxy_0[j] - f2a * 2.0 * s_xy_0[j]);

                s_xxxz_0[j] = fr * s_xxz_0[j] + f2t * 2.0 * s_xz_0[j];

                t_xxxz_0[j] = fr * t_xxz_0[j] + f2t * 2.0 * t_xz_0[j] + f2z * (s_xxxz_0[j] - f2a * 2.0 * s_xz_0[j]);

                s_xxyy_0[j] = fr * s_xyy_0[j] + f2t * s_yy_0[j];

                t_xxyy_0[j] = fr * t_xyy_0[j] + f2t * t_yy_0[j] + f2z * (s_xxyy_0[j] - f2a * s_yy_0[j]);

                s_xxyz_0[j] = fr * s_xyz_0[j] + f2t * s_yz_0[j];

                t_xxyz_0[j] = fr * t_xyz_0[j] + f2t * t_yz_0[j] + f2z * (s_xxyz_0[j] - f2a * s_yz_0[j]);

                s_xxzz_0[j] = fr * s_xzz_0[j] + f2t * s_zz_0[j];

                t_xxzz_0[j] = fr * t_xzz_0[j] + f2t * t_zz_0[j] + f2z * (s_xxzz_0[j] - f2a * s_zz_0[j]);

                s_xyyy_0[j] = fr * s_yyy_0[j];

                t_xyyy_0[j] = fr * t_yyy_0[j] + f2z * s_xyyy_0[j];

                s_xyyz_0[j] = fr * s_yyz_0[j];

                t_xyyz_0[j] = fr * t_yyz_0[j] + f2z * s_xyyz_0[j];

                s_xyzz_0[j] = fr * s_yzz_0[j];

                t_xyzz_0[j] = fr * t_yzz_0[j] + f2z * s_xyzz_0[j];

                s_xzzz_0[j] = fr * s_zzz_0[j];

                t_xzzz_0[j] = fr * t_zzz_0[j] + f2z * s_xzzz_0[j];

                // leading y component

                fr = pay[j];

                s_yyyy_0[j] = fr * s_yyy_0[j] + f2t * 3.0 * s_yy_0[j];

                t_yyyy_0[j] = fr * t_yyy_0[j] + f2t * 3.0 * t_yy_0[j] + f2z * (s_yyyy_0[j] - f2a * 3.0 * s_yy_0[j]);

                s_yyyz_0[j] = fr * s_yyz_0[j] + f2t * 2.0 * s_yz_0[j];

                t_yyyz_0[j] = fr * t_yyz_0[j] + f2t * 2.0 * t_yz_0[j] + f2z * (s_yyyz_0[j] - f2a * 2.0 * s_yz_0[j]);

                s_yyzz_0[j] = fr * s_yzz_0[j] + f2t * s_zz_0[j];

                t_yyzz_0[j] = fr * t_yzz_0[j] + f2t * t_zz_0[j] + f2z * (s_yyzz_0[j] - f2a * s_zz_0[j]);

                s_yzzz_0[j] = fr * s_zzz_0[j];

                t_yzzz_0[j] = fr * t_zzz_0[j] + f2z * s_yzzz_0[j];

                // leading z component

                fr = paz[j];

                s_zzzz_0[j] = fr * s_zzz_0[j] + f2t * 3.0 * s_zz_0[j];

                t_zzzz_0[j] = fr * t_zzz_0[j] + f2t * 3.0 * t_zz_0[j] + f2z * (s_zzzz_0[j] - f2a * 3.0 * s_zz_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForPG(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|F) integrals

            auto s_0_xxx = primBuffer.data(skoff + 10 * idx);

            auto s_0_xxy = primBuffer.data(skoff + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(skoff + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(skoff + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(skoff + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(skoff + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(skoff + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(skoff + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(skoff + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(skoff + 10 * idx + 9);

            // set up pointers to (S|T|F) integrals

            auto t_0_xxx = primBuffer.data(kkoff + 10 * idx);

            auto t_0_xxy = primBuffer.data(kkoff + 10 * idx + 1);

            auto t_0_xxz = primBuffer.data(kkoff + 10 * idx + 2);

            auto t_0_xyy = primBuffer.data(kkoff + 10 * idx + 3);

            auto t_0_xyz = primBuffer.data(kkoff + 10 * idx + 4);

            auto t_0_xzz = primBuffer.data(kkoff + 10 * idx + 5);

            auto t_0_yyy = primBuffer.data(kkoff + 10 * idx + 6);

            auto t_0_yyz = primBuffer.data(kkoff + 10 * idx + 7);

            auto t_0_yzz = primBuffer.data(kkoff + 10 * idx + 8);

            auto t_0_zzz = primBuffer.data(kkoff + 10 * idx + 9);

            // set up pointers to (S|G) integrals

            auto s_0_xxxx = primBuffer.data(s1off + 15 * idx);

            auto s_0_xxxy = primBuffer.data(s1off + 15 * idx + 1);

            auto s_0_xxxz = primBuffer.data(s1off + 15 * idx + 2);

            auto s_0_xxyy = primBuffer.data(s1off + 15 * idx + 3);

            auto s_0_xxyz = primBuffer.data(s1off + 15 * idx + 4);

            auto s_0_xxzz = primBuffer.data(s1off + 15 * idx + 5);

            auto s_0_xyyy = primBuffer.data(s1off + 15 * idx + 6);

            auto s_0_xyyz = primBuffer.data(s1off + 15 * idx + 7);

            auto s_0_xyzz = primBuffer.data(s1off + 15 * idx + 8);

            auto s_0_xzzz = primBuffer.data(s1off + 15 * idx + 9);

            auto s_0_yyyy = primBuffer.data(s1off + 15 * idx + 10);

            auto s_0_yyyz = primBuffer.data(s1off + 15 * idx + 11);

            auto s_0_yyzz = primBuffer.data(s1off + 15 * idx + 12);

            auto s_0_yzzz = primBuffer.data(s1off + 15 * idx + 13);

            auto s_0_zzzz = primBuffer.data(s1off + 15 * idx + 14);

            // set up pointers to (S|T|G) integrals

            auto t_0_xxxx = primBuffer.data(k1off + 15 * idx);

            auto t_0_xxxy = primBuffer.data(k1off + 15 * idx + 1);

            auto t_0_xxxz = primBuffer.data(k1off + 15 * idx + 2);

            auto t_0_xxyy = primBuffer.data(k1off + 15 * idx + 3);

            auto t_0_xxyz = primBuffer.data(k1off + 15 * idx + 4);

            auto t_0_xxzz = primBuffer.data(k1off + 15 * idx + 5);

            auto t_0_xyyy = primBuffer.data(k1off + 15 * idx + 6);

            auto t_0_xyyz = primBuffer.data(k1off + 15 * idx + 7);

            auto t_0_xyzz = primBuffer.data(k1off + 15 * idx + 8);

            auto t_0_xzzz = primBuffer.data(k1off + 15 * idx + 9);

            auto t_0_yyyy = primBuffer.data(k1off + 15 * idx + 10);

            auto t_0_yyyz = primBuffer.data(k1off + 15 * idx + 11);

            auto t_0_yyzz = primBuffer.data(k1off + 15 * idx + 12);

            auto t_0_yzzz = primBuffer.data(k1off + 15 * idx + 13);

            auto t_0_zzzz = primBuffer.data(k1off + 15 * idx + 14);

            // set up pointers to (P|G) integrals

            auto s_x_xxxx = primBuffer.data(soff + 45 * idx);

            auto s_x_xxxy = primBuffer.data(soff + 45 * idx + 1);

            auto s_x_xxxz = primBuffer.data(soff + 45 * idx + 2);

            auto s_x_xxyy = primBuffer.data(soff + 45 * idx + 3);

            auto s_x_xxyz = primBuffer.data(soff + 45 * idx + 4);

            auto s_x_xxzz = primBuffer.data(soff + 45 * idx + 5);

            auto s_x_xyyy = primBuffer.data(soff + 45 * idx + 6);

            auto s_x_xyyz = primBuffer.data(soff + 45 * idx + 7);

            auto s_x_xyzz = primBuffer.data(soff + 45 * idx + 8);

            auto s_x_xzzz = primBuffer.data(soff + 45 * idx + 9);

            auto s_x_yyyy = primBuffer.data(soff + 45 * idx + 10);

            auto s_x_yyyz = primBuffer.data(soff + 45 * idx + 11);

            auto s_x_yyzz = primBuffer.data(soff + 45 * idx + 12);

            auto s_x_yzzz = primBuffer.data(soff + 45 * idx + 13);

            auto s_x_zzzz = primBuffer.data(soff + 45 * idx + 14);

            auto s_y_xxxx = primBuffer.data(soff + 45 * idx + 15);

            auto s_y_xxxy = primBuffer.data(soff + 45 * idx + 16);

            auto s_y_xxxz = primBuffer.data(soff + 45 * idx + 17);

            auto s_y_xxyy = primBuffer.data(soff + 45 * idx + 18);

            auto s_y_xxyz = primBuffer.data(soff + 45 * idx + 19);

            auto s_y_xxzz = primBuffer.data(soff + 45 * idx + 20);

            auto s_y_xyyy = primBuffer.data(soff + 45 * idx + 21);

            auto s_y_xyyz = primBuffer.data(soff + 45 * idx + 22);

            auto s_y_xyzz = primBuffer.data(soff + 45 * idx + 23);

            auto s_y_xzzz = primBuffer.data(soff + 45 * idx + 24);

            auto s_y_yyyy = primBuffer.data(soff + 45 * idx + 25);

            auto s_y_yyyz = primBuffer.data(soff + 45 * idx + 26);

            auto s_y_yyzz = primBuffer.data(soff + 45 * idx + 27);

            auto s_y_yzzz = primBuffer.data(soff + 45 * idx + 28);

            auto s_y_zzzz = primBuffer.data(soff + 45 * idx + 29);

            auto s_z_xxxx = primBuffer.data(soff + 45 * idx + 30);

            auto s_z_xxxy = primBuffer.data(soff + 45 * idx + 31);

            auto s_z_xxxz = primBuffer.data(soff + 45 * idx + 32);

            auto s_z_xxyy = primBuffer.data(soff + 45 * idx + 33);

            auto s_z_xxyz = primBuffer.data(soff + 45 * idx + 34);

            auto s_z_xxzz = primBuffer.data(soff + 45 * idx + 35);

            auto s_z_xyyy = primBuffer.data(soff + 45 * idx + 36);

            auto s_z_xyyz = primBuffer.data(soff + 45 * idx + 37);

            auto s_z_xyzz = primBuffer.data(soff + 45 * idx + 38);

            auto s_z_xzzz = primBuffer.data(soff + 45 * idx + 39);

            auto s_z_yyyy = primBuffer.data(soff + 45 * idx + 40);

            auto s_z_yyyz = primBuffer.data(soff + 45 * idx + 41);

            auto s_z_yyzz = primBuffer.data(soff + 45 * idx + 42);

            auto s_z_yzzz = primBuffer.data(soff + 45 * idx + 43);

            auto s_z_zzzz = primBuffer.data(soff + 45 * idx + 44);

            // set up pointers to (P|T|G) integrals

            auto t_x_xxxx = primBuffer.data(koff + 45 * idx);

            auto t_x_xxxy = primBuffer.data(koff + 45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(koff + 45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(koff + 45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(koff + 45 * idx + 4);

            auto t_x_xxzz = primBuffer.data(koff + 45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(koff + 45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(koff + 45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(koff + 45 * idx + 8);

            auto t_x_xzzz = primBuffer.data(koff + 45 * idx + 9);

            auto t_x_yyyy = primBuffer.data(koff + 45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(koff + 45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(koff + 45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(koff + 45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(koff + 45 * idx + 14);

            auto t_y_xxxx = primBuffer.data(koff + 45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(koff + 45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(koff + 45 * idx + 17);

            auto t_y_xxyy = primBuffer.data(koff + 45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(koff + 45 * idx + 19);

            auto t_y_xxzz = primBuffer.data(koff + 45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(koff + 45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(koff + 45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(koff + 45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(koff + 45 * idx + 24);

            auto t_y_yyyy = primBuffer.data(koff + 45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(koff + 45 * idx + 26);

            auto t_y_yyzz = primBuffer.data(koff + 45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(koff + 45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(koff + 45 * idx + 29);

            auto t_z_xxxx = primBuffer.data(koff + 45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(koff + 45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(koff + 45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(koff + 45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(koff + 45 * idx + 34);

            auto t_z_xxzz = primBuffer.data(koff + 45 * idx + 35);

            auto t_z_xyyy = primBuffer.data(koff + 45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(koff + 45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(koff + 45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(koff + 45 * idx + 39);

            auto t_z_yyyy = primBuffer.data(koff + 45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(koff + 45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(koff + 45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(koff + 45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(koff + 45 * idx + 44);

            #pragma omp simd aligned(fx, fz, pax, pay, paz, s_0_xxx, s_0_xxy, s_0_xxz,\
                                     s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy, s_0_yyz,\
                                     s_0_yzz, s_0_zzz, t_0_xxx, t_0_xxy, t_0_xxz,\
                                     t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy, t_0_yyz,\
                                     t_0_yzz, t_0_zzz, s_0_xxxx, s_0_xxxy, s_0_xxxz,\
                                     s_0_xxyy, s_0_xxyz, s_0_xxzz, s_0_xyyy, s_0_xyyz,\
                                     s_0_xyzz, s_0_xzzz, s_0_yyyy, s_0_yyyz, s_0_yyzz,\
                                     s_0_yzzz, s_0_zzzz, t_0_xxxx, t_0_xxxy, t_0_xxxz,\
                                     t_0_xxyy, t_0_xxyz, t_0_xxzz, t_0_xyyy, t_0_xyyz,\
                                     t_0_xyzz, t_0_xzzz, t_0_yyyy, t_0_yyyz, t_0_yyzz,\
                                     t_0_yzzz, t_0_zzzz, s_x_xxxx, s_x_xxxy, s_x_xxxz,\
                                     s_x_xxyy, s_x_xxyz, s_x_xxzz, s_x_xyyy, s_x_xyyz,\
                                     s_x_xyzz, s_x_xzzz, s_x_yyyy, s_x_yyyz, s_x_yyzz,\
                                     s_x_yzzz, s_x_zzzz, s_y_xxxx, s_y_xxxy, s_y_xxxz,\
                                     s_y_xxyy, s_y_xxyz, s_y_xxzz, s_y_xyyy, s_y_xyyz,\
                                     s_y_xyzz, s_y_xzzz, s_y_yyyy, s_y_yyyz, s_y_yyzz,\
                                     s_y_yzzz, s_y_zzzz, s_z_xxxx, s_z_xxxy, s_z_xxxz,\
                                     s_z_xxyy, s_z_xxyz, s_z_xxzz, s_z_xyyy, s_z_xyyz,\
                                     s_z_xyzz, s_z_xzzz, s_z_yyyy, s_z_yyyz, s_z_yyzz,\
                                     s_z_yzzz, s_z_zzzz, t_x_xxxx, t_x_xxxy, t_x_xxxz,\
                                     t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy, t_x_xyyz,\
                                     t_x_xyzz, t_x_xzzz, t_x_yyyy, t_x_yyyz, t_x_yyzz,\
                                     t_x_yzzz, t_x_zzzz, t_y_xxxx, t_y_xxxy, t_y_xxxz,\
                                     t_y_xxyy, t_y_xxyz, t_y_xxzz, t_y_xyyy, t_y_xyyz,\
                                     t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz, t_y_yyzz,\
                                     t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy, t_z_xxxz,\
                                     t_z_xxyy, t_z_xxyz, t_z_xxzz, t_z_xyyy, t_z_xyyz,\
                                     t_z_xyzz, t_z_xzzz, t_z_yyyy, t_z_yyyz, t_z_yyzz,\
                                     t_z_yzzz, t_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                // leading x component

                double fr = pax[j];

                s_x_xxxx[j] = fr * s_0_xxxx[j] + f2t * 4.0 * s_0_xxx[j];

                t_x_xxxx[j] = fr * t_0_xxxx[j] + f2t * 4.0 * t_0_xxx[j] + f2z * s_x_xxxx[j];

                s_x_xxxy[j] = fr * s_0_xxxy[j] + f2t * 3.0 * s_0_xxy[j];

                t_x_xxxy[j] = fr * t_0_xxxy[j] + f2t * 3.0 * t_0_xxy[j] + f2z * s_x_xxxy[j];

                s_x_xxxz[j] = fr * s_0_xxxz[j] + f2t * 3.0 * s_0_xxz[j];

                t_x_xxxz[j] = fr * t_0_xxxz[j] + f2t * 3.0 * t_0_xxz[j] + f2z * s_x_xxxz[j];

                s_x_xxyy[j] = fr * s_0_xxyy[j] + f2t * 2.0 * s_0_xyy[j];

                t_x_xxyy[j] = fr * t_0_xxyy[j] + f2t * 2.0 * t_0_xyy[j] + f2z * s_x_xxyy[j];

                s_x_xxyz[j] = fr * s_0_xxyz[j] + f2t * 2.0 * s_0_xyz[j];

                t_x_xxyz[j] = fr * t_0_xxyz[j] + f2t * 2.0 * t_0_xyz[j] + f2z * s_x_xxyz[j];

                s_x_xxzz[j] = fr * s_0_xxzz[j] + f2t * 2.0 * s_0_xzz[j];

                t_x_xxzz[j] = fr * t_0_xxzz[j] + f2t * 2.0 * t_0_xzz[j] + f2z * s_x_xxzz[j];

                s_x_xyyy[j] = fr * s_0_xyyy[j] + f2t * s_0_yyy[j];

                t_x_xyyy[j] = fr * t_0_xyyy[j] + f2t * t_0_yyy[j] + f2z * s_x_xyyy[j];

                s_x_xyyz[j] = fr * s_0_xyyz[j] + f2t * s_0_yyz[j];

                t_x_xyyz[j] = fr * t_0_xyyz[j] + f2t * t_0_yyz[j] + f2z * s_x_xyyz[j];

                s_x_xyzz[j] = fr * s_0_xyzz[j] + f2t * s_0_yzz[j];

                t_x_xyzz[j] = fr * t_0_xyzz[j] + f2t * t_0_yzz[j] + f2z * s_x_xyzz[j];

                s_x_xzzz[j] = fr * s_0_xzzz[j] + f2t * s_0_zzz[j];

                t_x_xzzz[j] = fr * t_0_xzzz[j] + f2t * t_0_zzz[j] + f2z * s_x_xzzz[j];

                s_x_yyyy[j] = fr * s_0_yyyy[j];

                t_x_yyyy[j] = fr * t_0_yyyy[j] + f2z * s_x_yyyy[j];

                s_x_yyyz[j] = fr * s_0_yyyz[j];

                t_x_yyyz[j] = fr * t_0_yyyz[j] + f2z * s_x_yyyz[j];

                s_x_yyzz[j] = fr * s_0_yyzz[j];

                t_x_yyzz[j] = fr * t_0_yyzz[j] + f2z * s_x_yyzz[j];

                s_x_yzzz[j] = fr * s_0_yzzz[j];

                t_x_yzzz[j] = fr * t_0_yzzz[j] + f2z * s_x_yzzz[j];

                s_x_zzzz[j] = fr * s_0_zzzz[j];

                t_x_zzzz[j] = fr * t_0_zzzz[j] + f2z * s_x_zzzz[j];

                // leading y component

                fr = pay[j];

                s_y_xxxx[j] = fr * s_0_xxxx[j];

                t_y_xxxx[j] = fr * t_0_xxxx[j] + f2z * s_y_xxxx[j];

                s_y_xxxy[j] = fr * s_0_xxxy[j] + f2t * s_0_xxx[j];

                t_y_xxxy[j] = fr * t_0_xxxy[j] + f2t * t_0_xxx[j] + f2z * s_y_xxxy[j];

                s_y_xxxz[j] = fr * s_0_xxxz[j];

                t_y_xxxz[j] = fr * t_0_xxxz[j] + f2z * s_y_xxxz[j];

                s_y_xxyy[j] = fr * s_0_xxyy[j] + f2t * 2.0 * s_0_xxy[j];

                t_y_xxyy[j] = fr * t_0_xxyy[j] + f2t * 2.0 * t_0_xxy[j] + f2z * s_y_xxyy[j];

                s_y_xxyz[j] = fr * s_0_xxyz[j] + f2t * s_0_xxz[j];

                t_y_xxyz[j] = fr * t_0_xxyz[j] + f2t * t_0_xxz[j] + f2z * s_y_xxyz[j];

                s_y_xxzz[j] = fr * s_0_xxzz[j];

                t_y_xxzz[j] = fr * t_0_xxzz[j] + f2z * s_y_xxzz[j];

                s_y_xyyy[j] = fr * s_0_xyyy[j] + f2t * 3.0 * s_0_xyy[j];

                t_y_xyyy[j] = fr * t_0_xyyy[j] + f2t * 3.0 * t_0_xyy[j] + f2z * s_y_xyyy[j];

                s_y_xyyz[j] = fr * s_0_xyyz[j] + f2t * 2.0 * s_0_xyz[j];

                t_y_xyyz[j] = fr * t_0_xyyz[j] + f2t * 2.0 * t_0_xyz[j] + f2z * s_y_xyyz[j];

                s_y_xyzz[j] = fr * s_0_xyzz[j] + f2t * s_0_xzz[j];

                t_y_xyzz[j] = fr * t_0_xyzz[j] + f2t * t_0_xzz[j] + f2z * s_y_xyzz[j];

                s_y_xzzz[j] = fr * s_0_xzzz[j];

                t_y_xzzz[j] = fr * t_0_xzzz[j] + f2z * s_y_xzzz[j];

                s_y_yyyy[j] = fr * s_0_yyyy[j] + f2t * 4.0 * s_0_yyy[j];

                t_y_yyyy[j] = fr * t_0_yyyy[j] + f2t * 4.0 * t_0_yyy[j] + f2z * s_y_yyyy[j];

                s_y_yyyz[j] = fr * s_0_yyyz[j] + f2t * 3.0 * s_0_yyz[j];

                t_y_yyyz[j] = fr * t_0_yyyz[j] + f2t * 3.0 * t_0_yyz[j] + f2z * s_y_yyyz[j];

                s_y_yyzz[j] = fr * s_0_yyzz[j] + f2t * 2.0 * s_0_yzz[j];

                t_y_yyzz[j] = fr * t_0_yyzz[j] + f2t * 2.0 * t_0_yzz[j] + f2z * s_y_yyzz[j];

                s_y_yzzz[j] = fr * s_0_yzzz[j] + f2t * s_0_zzz[j];

                t_y_yzzz[j] = fr * t_0_yzzz[j] + f2t * t_0_zzz[j] + f2z * s_y_yzzz[j];

                s_y_zzzz[j] = fr * s_0_zzzz[j];

                t_y_zzzz[j] = fr * t_0_zzzz[j] + f2z * s_y_zzzz[j];

                // leading z component

                fr = paz[j];

                s_z_xxxx[j] = fr * s_0_xxxx[j];

                t_z_xxxx[j] = fr * t_0_xxxx[j] + f2z * s_z_xxxx[j];

                s_z_xxxy[j] = fr * s_0_xxxy[j];

                t_z_xxxy[j] = fr * t_0_xxxy[j] + f2z * s_z_xxxy[j];

                s_z_xxxz[j] = fr * s_0_xxxz[j] + f2t * s_0_xxx[j];

                t_z_xxxz[j] = fr * t_0_xxxz[j] + f2t * t_0_xxx[j] + f2z * s_z_xxxz[j];

                s_z_xxyy[j] = fr * s_0_xxyy[j];

                t_z_xxyy[j] = fr * t_0_xxyy[j] + f2z * s_z_xxyy[j];

                s_z_xxyz[j] = fr * s_0_xxyz[j] + f2t * s_0_xxy[j];

                t_z_xxyz[j] = fr * t_0_xxyz[j] + f2t * t_0_xxy[j] + f2z * s_z_xxyz[j];

                s_z_xxzz[j] = fr * s_0_xxzz[j] + f2t * 2.0 * s_0_xxz[j];

                t_z_xxzz[j] = fr * t_0_xxzz[j] + f2t * 2.0 * t_0_xxz[j] + f2z * s_z_xxzz[j];

                s_z_xyyy[j] = fr * s_0_xyyy[j];

                t_z_xyyy[j] = fr * t_0_xyyy[j] + f2z * s_z_xyyy[j];

                s_z_xyyz[j] = fr * s_0_xyyz[j] + f2t * s_0_xyy[j];

                t_z_xyyz[j] = fr * t_0_xyyz[j] + f2t * t_0_xyy[j] + f2z * s_z_xyyz[j];

                s_z_xyzz[j] = fr * s_0_xyzz[j] + f2t * 2.0 * s_0_xyz[j];

                t_z_xyzz[j] = fr * t_0_xyzz[j] + f2t * 2.0 * t_0_xyz[j] + f2z * s_z_xyzz[j];

                s_z_xzzz[j] = fr * s_0_xzzz[j] + f2t * 3.0 * s_0_xzz[j];

                t_z_xzzz[j] = fr * t_0_xzzz[j] + f2t * 3.0 * t_0_xzz[j] + f2z * s_z_xzzz[j];

                s_z_yyyy[j] = fr * s_0_yyyy[j];

                t_z_yyyy[j] = fr * t_0_yyyy[j] + f2z * s_z_yyyy[j];

                s_z_yyyz[j] = fr * s_0_yyyz[j] + f2t * s_0_yyy[j];

                t_z_yyyz[j] = fr * t_0_yyyz[j] + f2t * t_0_yyy[j] + f2z * s_z_yyyz[j];

                s_z_yyzz[j] = fr * s_0_yyzz[j] + f2t * 2.0 * s_0_yyz[j];

                t_z_yyzz[j] = fr * t_0_yyzz[j] + f2t * 2.0 * t_0_yyz[j] + f2z * s_z_yyzz[j];

                s_z_yzzz[j] = fr * s_0_yzzz[j] + f2t * 3.0 * s_0_yzz[j];

                t_z_yzzz[j] = fr * t_0_yzzz[j] + f2t * 3.0 * t_0_yzz[j] + f2z * s_z_yzzz[j];

                s_z_zzzz[j] = fr * s_0_zzzz[j] + f2t * 4.0 * s_0_zzz[j];

                t_z_zzzz[j] = fr * t_0_zzzz[j] + f2t * 4.0 * t_0_zzz[j] + f2z * s_z_zzzz[j];
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGP(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 1, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 1, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|S) integrals

            auto s_xxx_0 = primBuffer.data(skoff + 10 * idx);

            auto s_xxy_0 = primBuffer.data(skoff + 10 * idx + 1);

            auto s_xxz_0 = primBuffer.data(skoff + 10 * idx + 2);

            auto s_xyy_0 = primBuffer.data(skoff + 10 * idx + 3);

            auto s_xyz_0 = primBuffer.data(skoff + 10 * idx + 4);

            auto s_xzz_0 = primBuffer.data(skoff + 10 * idx + 5);

            auto s_yyy_0 = primBuffer.data(skoff + 10 * idx + 6);

            auto s_yyz_0 = primBuffer.data(skoff + 10 * idx + 7);

            auto s_yzz_0 = primBuffer.data(skoff + 10 * idx + 8);

            auto s_zzz_0 = primBuffer.data(skoff + 10 * idx + 9);

            // set up pointers to (F|T|S) integrals

            auto t_xxx_0 = primBuffer.data(kkoff + 10 * idx);

            auto t_xxy_0 = primBuffer.data(kkoff + 10 * idx + 1);

            auto t_xxz_0 = primBuffer.data(kkoff + 10 * idx + 2);

            auto t_xyy_0 = primBuffer.data(kkoff + 10 * idx + 3);

            auto t_xyz_0 = primBuffer.data(kkoff + 10 * idx + 4);

            auto t_xzz_0 = primBuffer.data(kkoff + 10 * idx + 5);

            auto t_yyy_0 = primBuffer.data(kkoff + 10 * idx + 6);

            auto t_yyz_0 = primBuffer.data(kkoff + 10 * idx + 7);

            auto t_yzz_0 = primBuffer.data(kkoff + 10 * idx + 8);

            auto t_zzz_0 = primBuffer.data(kkoff + 10 * idx + 9);

            // set up pointers to (D|P) integrals

            auto s_xx_x = primBuffer.data(s2off + 18 * idx);

            auto s_xx_y = primBuffer.data(s2off + 18 * idx + 1);

            auto s_xx_z = primBuffer.data(s2off + 18 * idx + 2);

            auto s_xy_x = primBuffer.data(s2off + 18 * idx + 3);

            auto s_xy_y = primBuffer.data(s2off + 18 * idx + 4);

            auto s_xy_z = primBuffer.data(s2off + 18 * idx + 5);

            auto s_xz_x = primBuffer.data(s2off + 18 * idx + 6);

            auto s_xz_y = primBuffer.data(s2off + 18 * idx + 7);

            auto s_xz_z = primBuffer.data(s2off + 18 * idx + 8);

            auto s_yy_x = primBuffer.data(s2off + 18 * idx + 9);

            auto s_yy_y = primBuffer.data(s2off + 18 * idx + 10);

            auto s_yy_z = primBuffer.data(s2off + 18 * idx + 11);

            auto s_yz_x = primBuffer.data(s2off + 18 * idx + 12);

            auto s_yz_y = primBuffer.data(s2off + 18 * idx + 13);

            auto s_yz_z = primBuffer.data(s2off + 18 * idx + 14);

            auto s_zz_x = primBuffer.data(s2off + 18 * idx + 15);

            auto s_zz_y = primBuffer.data(s2off + 18 * idx + 16);

            auto s_zz_z = primBuffer.data(s2off + 18 * idx + 17);

            // set up pointers to (D|T|P) integrals

            auto t_xx_x = primBuffer.data(k2off + 18 * idx);

            auto t_xx_y = primBuffer.data(k2off + 18 * idx + 1);

            auto t_xx_z = primBuffer.data(k2off + 18 * idx + 2);

            auto t_xy_x = primBuffer.data(k2off + 18 * idx + 3);

            auto t_xy_y = primBuffer.data(k2off + 18 * idx + 4);

            auto t_xy_z = primBuffer.data(k2off + 18 * idx + 5);

            auto t_xz_x = primBuffer.data(k2off + 18 * idx + 6);

            auto t_xz_y = primBuffer.data(k2off + 18 * idx + 7);

            auto t_xz_z = primBuffer.data(k2off + 18 * idx + 8);

            auto t_yy_x = primBuffer.data(k2off + 18 * idx + 9);

            auto t_yy_y = primBuffer.data(k2off + 18 * idx + 10);

            auto t_yy_z = primBuffer.data(k2off + 18 * idx + 11);

            auto t_yz_x = primBuffer.data(k2off + 18 * idx + 12);

            auto t_yz_y = primBuffer.data(k2off + 18 * idx + 13);

            auto t_yz_z = primBuffer.data(k2off + 18 * idx + 14);

            auto t_zz_x = primBuffer.data(k2off + 18 * idx + 15);

            auto t_zz_y = primBuffer.data(k2off + 18 * idx + 16);

            auto t_zz_z = primBuffer.data(k2off + 18 * idx + 17);

            // set up pointers to (F|P) integrals

            auto s_xxx_x = primBuffer.data(s1off + 30 * idx);

            auto s_xxx_y = primBuffer.data(s1off + 30 * idx + 1);

            auto s_xxx_z = primBuffer.data(s1off + 30 * idx + 2);

            auto s_xxy_x = primBuffer.data(s1off + 30 * idx + 3);

            auto s_xxy_y = primBuffer.data(s1off + 30 * idx + 4);

            auto s_xxy_z = primBuffer.data(s1off + 30 * idx + 5);

            auto s_xxz_x = primBuffer.data(s1off + 30 * idx + 6);

            auto s_xxz_y = primBuffer.data(s1off + 30 * idx + 7);

            auto s_xxz_z = primBuffer.data(s1off + 30 * idx + 8);

            auto s_xyy_x = primBuffer.data(s1off + 30 * idx + 9);

            auto s_xyy_y = primBuffer.data(s1off + 30 * idx + 10);

            auto s_xyy_z = primBuffer.data(s1off + 30 * idx + 11);

            auto s_xyz_x = primBuffer.data(s1off + 30 * idx + 12);

            auto s_xyz_y = primBuffer.data(s1off + 30 * idx + 13);

            auto s_xyz_z = primBuffer.data(s1off + 30 * idx + 14);

            auto s_xzz_x = primBuffer.data(s1off + 30 * idx + 15);

            auto s_xzz_y = primBuffer.data(s1off + 30 * idx + 16);

            auto s_xzz_z = primBuffer.data(s1off + 30 * idx + 17);

            auto s_yyy_x = primBuffer.data(s1off + 30 * idx + 18);

            auto s_yyy_y = primBuffer.data(s1off + 30 * idx + 19);

            auto s_yyy_z = primBuffer.data(s1off + 30 * idx + 20);

            auto s_yyz_x = primBuffer.data(s1off + 30 * idx + 21);

            auto s_yyz_y = primBuffer.data(s1off + 30 * idx + 22);

            auto s_yyz_z = primBuffer.data(s1off + 30 * idx + 23);

            auto s_yzz_x = primBuffer.data(s1off + 30 * idx + 24);

            auto s_yzz_y = primBuffer.data(s1off + 30 * idx + 25);

            auto s_yzz_z = primBuffer.data(s1off + 30 * idx + 26);

            auto s_zzz_x = primBuffer.data(s1off + 30 * idx + 27);

            auto s_zzz_y = primBuffer.data(s1off + 30 * idx + 28);

            auto s_zzz_z = primBuffer.data(s1off + 30 * idx + 29);

            // set up pointers to (F|T|P) integrals

            auto t_xxx_x = primBuffer.data(k1off + 30 * idx);

            auto t_xxx_y = primBuffer.data(k1off + 30 * idx + 1);

            auto t_xxx_z = primBuffer.data(k1off + 30 * idx + 2);

            auto t_xxy_x = primBuffer.data(k1off + 30 * idx + 3);

            auto t_xxy_y = primBuffer.data(k1off + 30 * idx + 4);

            auto t_xxy_z = primBuffer.data(k1off + 30 * idx + 5);

            auto t_xxz_x = primBuffer.data(k1off + 30 * idx + 6);

            auto t_xxz_y = primBuffer.data(k1off + 30 * idx + 7);

            auto t_xxz_z = primBuffer.data(k1off + 30 * idx + 8);

            auto t_xyy_x = primBuffer.data(k1off + 30 * idx + 9);

            auto t_xyy_y = primBuffer.data(k1off + 30 * idx + 10);

            auto t_xyy_z = primBuffer.data(k1off + 30 * idx + 11);

            auto t_xyz_x = primBuffer.data(k1off + 30 * idx + 12);

            auto t_xyz_y = primBuffer.data(k1off + 30 * idx + 13);

            auto t_xyz_z = primBuffer.data(k1off + 30 * idx + 14);

            auto t_xzz_x = primBuffer.data(k1off + 30 * idx + 15);

            auto t_xzz_y = primBuffer.data(k1off + 30 * idx + 16);

            auto t_xzz_z = primBuffer.data(k1off + 30 * idx + 17);

            auto t_yyy_x = primBuffer.data(k1off + 30 * idx + 18);

            auto t_yyy_y = primBuffer.data(k1off + 30 * idx + 19);

            auto t_yyy_z = primBuffer.data(k1off + 30 * idx + 20);

            auto t_yyz_x = primBuffer.data(k1off + 30 * idx + 21);

            auto t_yyz_y = primBuffer.data(k1off + 30 * idx + 22);

            auto t_yyz_z = primBuffer.data(k1off + 30 * idx + 23);

            auto t_yzz_x = primBuffer.data(k1off + 30 * idx + 24);

            auto t_yzz_y = primBuffer.data(k1off + 30 * idx + 25);

            auto t_yzz_z = primBuffer.data(k1off + 30 * idx + 26);

            auto t_zzz_x = primBuffer.data(k1off + 30 * idx + 27);

            auto t_zzz_y = primBuffer.data(k1off + 30 * idx + 28);

            auto t_zzz_z = primBuffer.data(k1off + 30 * idx + 29);

            // set up pointers to (G|P) integrals

            auto s_xxxx_x = primBuffer.data(soff + 45 * idx);

            auto s_xxxx_y = primBuffer.data(soff + 45 * idx + 1);

            auto s_xxxx_z = primBuffer.data(soff + 45 * idx + 2);

            auto s_xxxy_x = primBuffer.data(soff + 45 * idx + 3);

            auto s_xxxy_y = primBuffer.data(soff + 45 * idx + 4);

            auto s_xxxy_z = primBuffer.data(soff + 45 * idx + 5);

            auto s_xxxz_x = primBuffer.data(soff + 45 * idx + 6);

            auto s_xxxz_y = primBuffer.data(soff + 45 * idx + 7);

            auto s_xxxz_z = primBuffer.data(soff + 45 * idx + 8);

            auto s_xxyy_x = primBuffer.data(soff + 45 * idx + 9);

            auto s_xxyy_y = primBuffer.data(soff + 45 * idx + 10);

            auto s_xxyy_z = primBuffer.data(soff + 45 * idx + 11);

            auto s_xxyz_x = primBuffer.data(soff + 45 * idx + 12);

            auto s_xxyz_y = primBuffer.data(soff + 45 * idx + 13);

            auto s_xxyz_z = primBuffer.data(soff + 45 * idx + 14);

            auto s_xxzz_x = primBuffer.data(soff + 45 * idx + 15);

            auto s_xxzz_y = primBuffer.data(soff + 45 * idx + 16);

            auto s_xxzz_z = primBuffer.data(soff + 45 * idx + 17);

            auto s_xyyy_x = primBuffer.data(soff + 45 * idx + 18);

            auto s_xyyy_y = primBuffer.data(soff + 45 * idx + 19);

            auto s_xyyy_z = primBuffer.data(soff + 45 * idx + 20);

            auto s_xyyz_x = primBuffer.data(soff + 45 * idx + 21);

            auto s_xyyz_y = primBuffer.data(soff + 45 * idx + 22);

            auto s_xyyz_z = primBuffer.data(soff + 45 * idx + 23);

            auto s_xyzz_x = primBuffer.data(soff + 45 * idx + 24);

            auto s_xyzz_y = primBuffer.data(soff + 45 * idx + 25);

            auto s_xyzz_z = primBuffer.data(soff + 45 * idx + 26);

            auto s_xzzz_x = primBuffer.data(soff + 45 * idx + 27);

            auto s_xzzz_y = primBuffer.data(soff + 45 * idx + 28);

            auto s_xzzz_z = primBuffer.data(soff + 45 * idx + 29);

            auto s_yyyy_x = primBuffer.data(soff + 45 * idx + 30);

            auto s_yyyy_y = primBuffer.data(soff + 45 * idx + 31);

            auto s_yyyy_z = primBuffer.data(soff + 45 * idx + 32);

            auto s_yyyz_x = primBuffer.data(soff + 45 * idx + 33);

            auto s_yyyz_y = primBuffer.data(soff + 45 * idx + 34);

            auto s_yyyz_z = primBuffer.data(soff + 45 * idx + 35);

            auto s_yyzz_x = primBuffer.data(soff + 45 * idx + 36);

            auto s_yyzz_y = primBuffer.data(soff + 45 * idx + 37);

            auto s_yyzz_z = primBuffer.data(soff + 45 * idx + 38);

            auto s_yzzz_x = primBuffer.data(soff + 45 * idx + 39);

            auto s_yzzz_y = primBuffer.data(soff + 45 * idx + 40);

            auto s_yzzz_z = primBuffer.data(soff + 45 * idx + 41);

            auto s_zzzz_x = primBuffer.data(soff + 45 * idx + 42);

            auto s_zzzz_y = primBuffer.data(soff + 45 * idx + 43);

            auto s_zzzz_z = primBuffer.data(soff + 45 * idx + 44);

            // set up pointers to (G|T|P) integrals

            auto t_xxxx_x = primBuffer.data(koff + 45 * idx);

            auto t_xxxx_y = primBuffer.data(koff + 45 * idx + 1);

            auto t_xxxx_z = primBuffer.data(koff + 45 * idx + 2);

            auto t_xxxy_x = primBuffer.data(koff + 45 * idx + 3);

            auto t_xxxy_y = primBuffer.data(koff + 45 * idx + 4);

            auto t_xxxy_z = primBuffer.data(koff + 45 * idx + 5);

            auto t_xxxz_x = primBuffer.data(koff + 45 * idx + 6);

            auto t_xxxz_y = primBuffer.data(koff + 45 * idx + 7);

            auto t_xxxz_z = primBuffer.data(koff + 45 * idx + 8);

            auto t_xxyy_x = primBuffer.data(koff + 45 * idx + 9);

            auto t_xxyy_y = primBuffer.data(koff + 45 * idx + 10);

            auto t_xxyy_z = primBuffer.data(koff + 45 * idx + 11);

            auto t_xxyz_x = primBuffer.data(koff + 45 * idx + 12);

            auto t_xxyz_y = primBuffer.data(koff + 45 * idx + 13);

            auto t_xxyz_z = primBuffer.data(koff + 45 * idx + 14);

            auto t_xxzz_x = primBuffer.data(koff + 45 * idx + 15);

            auto t_xxzz_y = primBuffer.data(koff + 45 * idx + 16);

            auto t_xxzz_z = primBuffer.data(koff + 45 * idx + 17);

            auto t_xyyy_x = primBuffer.data(koff + 45 * idx + 18);

            auto t_xyyy_y = primBuffer.data(koff + 45 * idx + 19);

            auto t_xyyy_z = primBuffer.data(koff + 45 * idx + 20);

            auto t_xyyz_x = primBuffer.data(koff + 45 * idx + 21);

            auto t_xyyz_y = primBuffer.data(koff + 45 * idx + 22);

            auto t_xyyz_z = primBuffer.data(koff + 45 * idx + 23);

            auto t_xyzz_x = primBuffer.data(koff + 45 * idx + 24);

            auto t_xyzz_y = primBuffer.data(koff + 45 * idx + 25);

            auto t_xyzz_z = primBuffer.data(koff + 45 * idx + 26);

            auto t_xzzz_x = primBuffer.data(koff + 45 * idx + 27);

            auto t_xzzz_y = primBuffer.data(koff + 45 * idx + 28);

            auto t_xzzz_z = primBuffer.data(koff + 45 * idx + 29);

            auto t_yyyy_x = primBuffer.data(koff + 45 * idx + 30);

            auto t_yyyy_y = primBuffer.data(koff + 45 * idx + 31);

            auto t_yyyy_z = primBuffer.data(koff + 45 * idx + 32);

            auto t_yyyz_x = primBuffer.data(koff + 45 * idx + 33);

            auto t_yyyz_y = primBuffer.data(koff + 45 * idx + 34);

            auto t_yyyz_z = primBuffer.data(koff + 45 * idx + 35);

            auto t_yyzz_x = primBuffer.data(koff + 45 * idx + 36);

            auto t_yyzz_y = primBuffer.data(koff + 45 * idx + 37);

            auto t_yyzz_z = primBuffer.data(koff + 45 * idx + 38);

            auto t_yzzz_x = primBuffer.data(koff + 45 * idx + 39);

            auto t_yzzz_y = primBuffer.data(koff + 45 * idx + 40);

            auto t_yzzz_z = primBuffer.data(koff + 45 * idx + 41);

            auto t_zzzz_x = primBuffer.data(koff + 45 * idx + 42);

            auto t_zzzz_y = primBuffer.data(koff + 45 * idx + 43);

            auto t_zzzz_z = primBuffer.data(koff + 45 * idx + 44);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xxx_0, s_xxy_0,\
                                     s_xxz_0, s_xyy_0, s_xyz_0, s_xzz_0, s_yyy_0,\
                                     s_yyz_0, s_yzz_0, s_zzz_0, t_xxx_0, t_xxy_0,\
                                     t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, t_yyy_0,\
                                     t_yyz_0, t_yzz_0, t_zzz_0, s_xx_x, s_xx_y,\
                                     s_xx_z, s_xy_x, s_xy_y, s_xy_z, s_xz_x, s_xz_y,\
                                     s_xz_z, s_yy_x, s_yy_y, s_yy_z, s_yz_x, s_yz_y,\
                                     s_yz_z, s_zz_x, s_zz_y, s_zz_z, t_xx_x, t_xx_y,\
                                     t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y,\
                                     t_xz_z, t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y,\
                                     t_yz_z, t_zz_x, t_zz_y, t_zz_z, s_xxx_x, s_xxx_y,\
                                     s_xxx_z, s_xxy_x, s_xxy_y, s_xxy_z, s_xxz_x,\
                                     s_xxz_y, s_xxz_z, s_xyy_x, s_xyy_y, s_xyy_z,\
                                     s_xyz_x, s_xyz_y, s_xyz_z, s_xzz_x, s_xzz_y,\
                                     s_xzz_z, s_yyy_x, s_yyy_y, s_yyy_z, s_yyz_x,\
                                     s_yyz_y, s_yyz_z, s_yzz_x, s_yzz_y, s_yzz_z,\
                                     s_zzz_x, s_zzz_y, s_zzz_z, t_xxx_x, t_xxx_y,\
                                     t_xxx_z, t_xxy_x, t_xxy_y, t_xxy_z, t_xxz_x,\
                                     t_xxz_y, t_xxz_z, t_xyy_x, t_xyy_y, t_xyy_z,\
                                     t_xyz_x, t_xyz_y, t_xyz_z, t_xzz_x, t_xzz_y,\
                                     t_xzz_z, t_yyy_x, t_yyy_y, t_yyy_z, t_yyz_x,\
                                     t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z,\
                                     t_zzz_x, t_zzz_y, t_zzz_z, s_xxxx_x, s_xxxx_y,\
                                     s_xxxx_z, s_xxxy_x, s_xxxy_y, s_xxxy_z, s_xxxz_x,\
                                     s_xxxz_y, s_xxxz_z, s_xxyy_x, s_xxyy_y, s_xxyy_z,\
                                     s_xxyz_x, s_xxyz_y, s_xxyz_z, s_xxzz_x, s_xxzz_y,\
                                     s_xxzz_z, s_xyyy_x, s_xyyy_y, s_xyyy_z, s_xyyz_x,\
                                     s_xyyz_y, s_xyyz_z, s_xyzz_x, s_xyzz_y, s_xyzz_z,\
                                     s_xzzz_x, s_xzzz_y, s_xzzz_z, s_yyyy_x, s_yyyy_y,\
                                     s_yyyy_z, s_yyyz_x, s_yyyz_y, s_yyyz_z, s_yyzz_x,\
                                     s_yyzz_y, s_yyzz_z, s_yzzz_x, s_yzzz_y, s_yzzz_z,\
                                     s_zzzz_x, s_zzzz_y, s_zzzz_z, t_xxxx_x, t_xxxx_y,\
                                     t_xxxx_z, t_xxxy_x, t_xxxy_y, t_xxxy_z, t_xxxz_x,\
                                     t_xxxz_y, t_xxxz_z, t_xxyy_x, t_xxyy_y, t_xxyy_z,\
                                     t_xxyz_x, t_xxyz_y, t_xxyz_z, t_xxzz_x, t_xxzz_y,\
                                     t_xxzz_z, t_xyyy_x, t_xyyy_y, t_xyyy_z, t_xyyz_x,\
                                     t_xyyz_y, t_xyyz_z, t_xyzz_x, t_xyzz_y, t_xyzz_z,\
                                     t_xzzz_x, t_xzzz_y, t_xzzz_z, t_yyyy_x, t_yyyy_y,\
                                     t_yyyy_z, t_yyyz_x, t_yyyz_y, t_yyyz_z, t_yyzz_x,\
                                     t_yyzz_y, t_yyzz_z, t_yzzz_x, t_yzzz_y, t_yzzz_z,\
                                     t_zzzz_x, t_zzzz_y, t_zzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxxx_x[j] = fr * s_xxx_x[j] + f2t * (3.0 * s_xx_x[j] + s_xxx_0[j]);

                t_xxxx_x[j] = fr * t_xxx_x[j] + f2t * (3.0 * t_xx_x[j] + t_xxx_0[j]) + f2z * (s_xxxx_x[j] - f2a * 3.0 * s_xx_x[j]);

                s_xxxx_y[j] = fr * s_xxx_y[j] + f2t * 3.0 * s_xx_y[j];

                t_xxxx_y[j] = fr * t_xxx_y[j] + f2t * 3.0 * t_xx_y[j] + f2z * (s_xxxx_y[j] - f2a * 3.0 * s_xx_y[j]);

                s_xxxx_z[j] = fr * s_xxx_z[j] + f2t * 3.0 * s_xx_z[j];

                t_xxxx_z[j] = fr * t_xxx_z[j] + f2t * 3.0 * t_xx_z[j] + f2z * (s_xxxx_z[j] - f2a * 3.0 * s_xx_z[j]);

                s_xxxy_x[j] = fr * s_xxy_x[j] + f2t * (2.0 * s_xy_x[j] + s_xxy_0[j]);

                t_xxxy_x[j] = fr * t_xxy_x[j] + f2t * (2.0 * t_xy_x[j] + t_xxy_0[j]) + f2z * (s_xxxy_x[j] - f2a * 2.0 * s_xy_x[j]);

                s_xxxy_y[j] = fr * s_xxy_y[j] + f2t * 2.0 * s_xy_y[j];

                t_xxxy_y[j] = fr * t_xxy_y[j] + f2t * 2.0 * t_xy_y[j] + f2z * (s_xxxy_y[j] - f2a * 2.0 * s_xy_y[j]);

                s_xxxy_z[j] = fr * s_xxy_z[j] + f2t * 2.0 * s_xy_z[j];

                t_xxxy_z[j] = fr * t_xxy_z[j] + f2t * 2.0 * t_xy_z[j] + f2z * (s_xxxy_z[j] - f2a * 2.0 * s_xy_z[j]);

                s_xxxz_x[j] = fr * s_xxz_x[j] + f2t * (2.0 * s_xz_x[j] + s_xxz_0[j]);

                t_xxxz_x[j] = fr * t_xxz_x[j] + f2t * (2.0 * t_xz_x[j] + t_xxz_0[j]) + f2z * (s_xxxz_x[j] - f2a * 2.0 * s_xz_x[j]);

                s_xxxz_y[j] = fr * s_xxz_y[j] + f2t * 2.0 * s_xz_y[j];

                t_xxxz_y[j] = fr * t_xxz_y[j] + f2t * 2.0 * t_xz_y[j] + f2z * (s_xxxz_y[j] - f2a * 2.0 * s_xz_y[j]);

                s_xxxz_z[j] = fr * s_xxz_z[j] + f2t * 2.0 * s_xz_z[j];

                t_xxxz_z[j] = fr * t_xxz_z[j] + f2t * 2.0 * t_xz_z[j] + f2z * (s_xxxz_z[j] - f2a * 2.0 * s_xz_z[j]);

                s_xxyy_x[j] = fr * s_xyy_x[j] + f2t * (s_yy_x[j] + s_xyy_0[j]);

                t_xxyy_x[j] = fr * t_xyy_x[j] + f2t * (t_yy_x[j] + t_xyy_0[j]) + f2z * (s_xxyy_x[j] - f2a * s_yy_x[j]);

                s_xxyy_y[j] = fr * s_xyy_y[j] + f2t * s_yy_y[j];

                t_xxyy_y[j] = fr * t_xyy_y[j] + f2t * t_yy_y[j] + f2z * (s_xxyy_y[j] - f2a * s_yy_y[j]);

                s_xxyy_z[j] = fr * s_xyy_z[j] + f2t * s_yy_z[j];

                t_xxyy_z[j] = fr * t_xyy_z[j] + f2t * t_yy_z[j] + f2z * (s_xxyy_z[j] - f2a * s_yy_z[j]);

                s_xxyz_x[j] = fr * s_xyz_x[j] + f2t * (s_yz_x[j] + s_xyz_0[j]);

                t_xxyz_x[j] = fr * t_xyz_x[j] + f2t * (t_yz_x[j] + t_xyz_0[j]) + f2z * (s_xxyz_x[j] - f2a * s_yz_x[j]);

                s_xxyz_y[j] = fr * s_xyz_y[j] + f2t * s_yz_y[j];

                t_xxyz_y[j] = fr * t_xyz_y[j] + f2t * t_yz_y[j] + f2z * (s_xxyz_y[j] - f2a * s_yz_y[j]);

                s_xxyz_z[j] = fr * s_xyz_z[j] + f2t * s_yz_z[j];

                t_xxyz_z[j] = fr * t_xyz_z[j] + f2t * t_yz_z[j] + f2z * (s_xxyz_z[j] - f2a * s_yz_z[j]);

                s_xxzz_x[j] = fr * s_xzz_x[j] + f2t * (s_zz_x[j] + s_xzz_0[j]);

                t_xxzz_x[j] = fr * t_xzz_x[j] + f2t * (t_zz_x[j] + t_xzz_0[j]) + f2z * (s_xxzz_x[j] - f2a * s_zz_x[j]);

                s_xxzz_y[j] = fr * s_xzz_y[j] + f2t * s_zz_y[j];

                t_xxzz_y[j] = fr * t_xzz_y[j] + f2t * t_zz_y[j] + f2z * (s_xxzz_y[j] - f2a * s_zz_y[j]);

                s_xxzz_z[j] = fr * s_xzz_z[j] + f2t * s_zz_z[j];

                t_xxzz_z[j] = fr * t_xzz_z[j] + f2t * t_zz_z[j] + f2z * (s_xxzz_z[j] - f2a * s_zz_z[j]);

                s_xyyy_x[j] = fr * s_yyy_x[j] + f2t * s_yyy_0[j];

                t_xyyy_x[j] = fr * t_yyy_x[j] + f2t * t_yyy_0[j] + f2z * s_xyyy_x[j];

                s_xyyy_y[j] = fr * s_yyy_y[j];

                t_xyyy_y[j] = fr * t_yyy_y[j] + f2z * s_xyyy_y[j];

                s_xyyy_z[j] = fr * s_yyy_z[j];

                t_xyyy_z[j] = fr * t_yyy_z[j] + f2z * s_xyyy_z[j];

                s_xyyz_x[j] = fr * s_yyz_x[j] + f2t * s_yyz_0[j];

                t_xyyz_x[j] = fr * t_yyz_x[j] + f2t * t_yyz_0[j] + f2z * s_xyyz_x[j];

                s_xyyz_y[j] = fr * s_yyz_y[j];

                t_xyyz_y[j] = fr * t_yyz_y[j] + f2z * s_xyyz_y[j];

                s_xyyz_z[j] = fr * s_yyz_z[j];

                t_xyyz_z[j] = fr * t_yyz_z[j] + f2z * s_xyyz_z[j];

                s_xyzz_x[j] = fr * s_yzz_x[j] + f2t * s_yzz_0[j];

                t_xyzz_x[j] = fr * t_yzz_x[j] + f2t * t_yzz_0[j] + f2z * s_xyzz_x[j];

                s_xyzz_y[j] = fr * s_yzz_y[j];

                t_xyzz_y[j] = fr * t_yzz_y[j] + f2z * s_xyzz_y[j];

                s_xyzz_z[j] = fr * s_yzz_z[j];

                t_xyzz_z[j] = fr * t_yzz_z[j] + f2z * s_xyzz_z[j];

                s_xzzz_x[j] = fr * s_zzz_x[j] + f2t * s_zzz_0[j];

                t_xzzz_x[j] = fr * t_zzz_x[j] + f2t * t_zzz_0[j] + f2z * s_xzzz_x[j];

                s_xzzz_y[j] = fr * s_zzz_y[j];

                t_xzzz_y[j] = fr * t_zzz_y[j] + f2z * s_xzzz_y[j];

                s_xzzz_z[j] = fr * s_zzz_z[j];

                t_xzzz_z[j] = fr * t_zzz_z[j] + f2z * s_xzzz_z[j];

                // leading y component

                fr = pay[j];

                s_yyyy_x[j] = fr * s_yyy_x[j] + f2t * 3.0 * s_yy_x[j];

                t_yyyy_x[j] = fr * t_yyy_x[j] + f2t * 3.0 * t_yy_x[j] + f2z * (s_yyyy_x[j] - f2a * 3.0 * s_yy_x[j]);

                s_yyyy_y[j] = fr * s_yyy_y[j] + f2t * (3.0 * s_yy_y[j] + s_yyy_0[j]);

                t_yyyy_y[j] = fr * t_yyy_y[j] + f2t * (3.0 * t_yy_y[j] + t_yyy_0[j]) + f2z * (s_yyyy_y[j] - f2a * 3.0 * s_yy_y[j]);

                s_yyyy_z[j] = fr * s_yyy_z[j] + f2t * 3.0 * s_yy_z[j];

                t_yyyy_z[j] = fr * t_yyy_z[j] + f2t * 3.0 * t_yy_z[j] + f2z * (s_yyyy_z[j] - f2a * 3.0 * s_yy_z[j]);

                s_yyyz_x[j] = fr * s_yyz_x[j] + f2t * 2.0 * s_yz_x[j];

                t_yyyz_x[j] = fr * t_yyz_x[j] + f2t * 2.0 * t_yz_x[j] + f2z * (s_yyyz_x[j] - f2a * 2.0 * s_yz_x[j]);

                s_yyyz_y[j] = fr * s_yyz_y[j] + f2t * (2.0 * s_yz_y[j] + s_yyz_0[j]);

                t_yyyz_y[j] = fr * t_yyz_y[j] + f2t * (2.0 * t_yz_y[j] + t_yyz_0[j]) + f2z * (s_yyyz_y[j] - f2a * 2.0 * s_yz_y[j]);

                s_yyyz_z[j] = fr * s_yyz_z[j] + f2t * 2.0 * s_yz_z[j];

                t_yyyz_z[j] = fr * t_yyz_z[j] + f2t * 2.0 * t_yz_z[j] + f2z * (s_yyyz_z[j] - f2a * 2.0 * s_yz_z[j]);

                s_yyzz_x[j] = fr * s_yzz_x[j] + f2t * s_zz_x[j];

                t_yyzz_x[j] = fr * t_yzz_x[j] + f2t * t_zz_x[j] + f2z * (s_yyzz_x[j] - f2a * s_zz_x[j]);

                s_yyzz_y[j] = fr * s_yzz_y[j] + f2t * (s_zz_y[j] + s_yzz_0[j]);

                t_yyzz_y[j] = fr * t_yzz_y[j] + f2t * (t_zz_y[j] + t_yzz_0[j]) + f2z * (s_yyzz_y[j] - f2a * s_zz_y[j]);

                s_yyzz_z[j] = fr * s_yzz_z[j] + f2t * s_zz_z[j];

                t_yyzz_z[j] = fr * t_yzz_z[j] + f2t * t_zz_z[j] + f2z * (s_yyzz_z[j] - f2a * s_zz_z[j]);

                s_yzzz_x[j] = fr * s_zzz_x[j];

                t_yzzz_x[j] = fr * t_zzz_x[j] + f2z * s_yzzz_x[j];

                s_yzzz_y[j] = fr * s_zzz_y[j] + f2t * s_zzz_0[j];

                t_yzzz_y[j] = fr * t_zzz_y[j] + f2t * t_zzz_0[j] + f2z * s_yzzz_y[j];

                s_yzzz_z[j] = fr * s_zzz_z[j];

                t_yzzz_z[j] = fr * t_zzz_z[j] + f2z * s_yzzz_z[j];

                // leading z component

                fr = paz[j];

                s_zzzz_x[j] = fr * s_zzz_x[j] + f2t * 3.0 * s_zz_x[j];

                t_zzzz_x[j] = fr * t_zzz_x[j] + f2t * 3.0 * t_zz_x[j] + f2z * (s_zzzz_x[j] - f2a * 3.0 * s_zz_x[j]);

                s_zzzz_y[j] = fr * s_zzz_y[j] + f2t * 3.0 * s_zz_y[j];

                t_zzzz_y[j] = fr * t_zzz_y[j] + f2t * 3.0 * t_zz_y[j] + f2z * (s_zzzz_y[j] - f2a * 3.0 * s_zz_y[j]);

                s_zzzz_z[j] = fr * s_zzz_z[j] + f2t * (3.0 * s_zz_z[j] + s_zzz_0[j]);

                t_zzzz_z[j] = fr * t_zzz_z[j] + f2t * (3.0 * t_zz_z[j] + t_zzz_0[j]) + f2z * (s_zzzz_z[j] - f2a * 3.0 * s_zz_z[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDG(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|F) integrals

            auto s_x_xxx = primBuffer.data(skoff + 30 * idx);

            auto s_x_xxy = primBuffer.data(skoff + 30 * idx + 1);

            auto s_x_xxz = primBuffer.data(skoff + 30 * idx + 2);

            auto s_x_xyy = primBuffer.data(skoff + 30 * idx + 3);

            auto s_x_xyz = primBuffer.data(skoff + 30 * idx + 4);

            auto s_x_xzz = primBuffer.data(skoff + 30 * idx + 5);

            auto s_x_yyy = primBuffer.data(skoff + 30 * idx + 6);

            auto s_x_yyz = primBuffer.data(skoff + 30 * idx + 7);

            auto s_x_yzz = primBuffer.data(skoff + 30 * idx + 8);

            auto s_x_zzz = primBuffer.data(skoff + 30 * idx + 9);

            auto s_y_xxx = primBuffer.data(skoff + 30 * idx + 10);

            auto s_y_xxy = primBuffer.data(skoff + 30 * idx + 11);

            auto s_y_xxz = primBuffer.data(skoff + 30 * idx + 12);

            auto s_y_xyy = primBuffer.data(skoff + 30 * idx + 13);

            auto s_y_xyz = primBuffer.data(skoff + 30 * idx + 14);

            auto s_y_xzz = primBuffer.data(skoff + 30 * idx + 15);

            auto s_y_yyy = primBuffer.data(skoff + 30 * idx + 16);

            auto s_y_yyz = primBuffer.data(skoff + 30 * idx + 17);

            auto s_y_yzz = primBuffer.data(skoff + 30 * idx + 18);

            auto s_y_zzz = primBuffer.data(skoff + 30 * idx + 19);

            auto s_z_xxx = primBuffer.data(skoff + 30 * idx + 20);

            auto s_z_xxy = primBuffer.data(skoff + 30 * idx + 21);

            auto s_z_xxz = primBuffer.data(skoff + 30 * idx + 22);

            auto s_z_xyy = primBuffer.data(skoff + 30 * idx + 23);

            auto s_z_xyz = primBuffer.data(skoff + 30 * idx + 24);

            auto s_z_xzz = primBuffer.data(skoff + 30 * idx + 25);

            auto s_z_yyy = primBuffer.data(skoff + 30 * idx + 26);

            auto s_z_yyz = primBuffer.data(skoff + 30 * idx + 27);

            auto s_z_yzz = primBuffer.data(skoff + 30 * idx + 28);

            auto s_z_zzz = primBuffer.data(skoff + 30 * idx + 29);

            // set up pointers to (P|T|F) integrals

            auto t_x_xxx = primBuffer.data(kkoff + 30 * idx);

            auto t_x_xxy = primBuffer.data(kkoff + 30 * idx + 1);

            auto t_x_xxz = primBuffer.data(kkoff + 30 * idx + 2);

            auto t_x_xyy = primBuffer.data(kkoff + 30 * idx + 3);

            auto t_x_xyz = primBuffer.data(kkoff + 30 * idx + 4);

            auto t_x_xzz = primBuffer.data(kkoff + 30 * idx + 5);

            auto t_x_yyy = primBuffer.data(kkoff + 30 * idx + 6);

            auto t_x_yyz = primBuffer.data(kkoff + 30 * idx + 7);

            auto t_x_yzz = primBuffer.data(kkoff + 30 * idx + 8);

            auto t_x_zzz = primBuffer.data(kkoff + 30 * idx + 9);

            auto t_y_xxx = primBuffer.data(kkoff + 30 * idx + 10);

            auto t_y_xxy = primBuffer.data(kkoff + 30 * idx + 11);

            auto t_y_xxz = primBuffer.data(kkoff + 30 * idx + 12);

            auto t_y_xyy = primBuffer.data(kkoff + 30 * idx + 13);

            auto t_y_xyz = primBuffer.data(kkoff + 30 * idx + 14);

            auto t_y_xzz = primBuffer.data(kkoff + 30 * idx + 15);

            auto t_y_yyy = primBuffer.data(kkoff + 30 * idx + 16);

            auto t_y_yyz = primBuffer.data(kkoff + 30 * idx + 17);

            auto t_y_yzz = primBuffer.data(kkoff + 30 * idx + 18);

            auto t_y_zzz = primBuffer.data(kkoff + 30 * idx + 19);

            auto t_z_xxx = primBuffer.data(kkoff + 30 * idx + 20);

            auto t_z_xxy = primBuffer.data(kkoff + 30 * idx + 21);

            auto t_z_xxz = primBuffer.data(kkoff + 30 * idx + 22);

            auto t_z_xyy = primBuffer.data(kkoff + 30 * idx + 23);

            auto t_z_xyz = primBuffer.data(kkoff + 30 * idx + 24);

            auto t_z_xzz = primBuffer.data(kkoff + 30 * idx + 25);

            auto t_z_yyy = primBuffer.data(kkoff + 30 * idx + 26);

            auto t_z_yyz = primBuffer.data(kkoff + 30 * idx + 27);

            auto t_z_yzz = primBuffer.data(kkoff + 30 * idx + 28);

            auto t_z_zzz = primBuffer.data(kkoff + 30 * idx + 29);

            // set up pointers to (S|G) integrals

            auto s_0_xxxx = primBuffer.data(s2off + 15 * idx);

            auto s_0_xxxy = primBuffer.data(s2off + 15 * idx + 1);

            auto s_0_xxxz = primBuffer.data(s2off + 15 * idx + 2);

            auto s_0_xxyy = primBuffer.data(s2off + 15 * idx + 3);

            auto s_0_xxyz = primBuffer.data(s2off + 15 * idx + 4);

            auto s_0_xxzz = primBuffer.data(s2off + 15 * idx + 5);

            auto s_0_xyyy = primBuffer.data(s2off + 15 * idx + 6);

            auto s_0_xyyz = primBuffer.data(s2off + 15 * idx + 7);

            auto s_0_xyzz = primBuffer.data(s2off + 15 * idx + 8);

            auto s_0_xzzz = primBuffer.data(s2off + 15 * idx + 9);

            auto s_0_yyyy = primBuffer.data(s2off + 15 * idx + 10);

            auto s_0_yyyz = primBuffer.data(s2off + 15 * idx + 11);

            auto s_0_yyzz = primBuffer.data(s2off + 15 * idx + 12);

            auto s_0_yzzz = primBuffer.data(s2off + 15 * idx + 13);

            auto s_0_zzzz = primBuffer.data(s2off + 15 * idx + 14);

            // set up pointers to (S|T|G) integrals

            auto t_0_xxxx = primBuffer.data(k2off + 15 * idx);

            auto t_0_xxxy = primBuffer.data(k2off + 15 * idx + 1);

            auto t_0_xxxz = primBuffer.data(k2off + 15 * idx + 2);

            auto t_0_xxyy = primBuffer.data(k2off + 15 * idx + 3);

            auto t_0_xxyz = primBuffer.data(k2off + 15 * idx + 4);

            auto t_0_xxzz = primBuffer.data(k2off + 15 * idx + 5);

            auto t_0_xyyy = primBuffer.data(k2off + 15 * idx + 6);

            auto t_0_xyyz = primBuffer.data(k2off + 15 * idx + 7);

            auto t_0_xyzz = primBuffer.data(k2off + 15 * idx + 8);

            auto t_0_xzzz = primBuffer.data(k2off + 15 * idx + 9);

            auto t_0_yyyy = primBuffer.data(k2off + 15 * idx + 10);

            auto t_0_yyyz = primBuffer.data(k2off + 15 * idx + 11);

            auto t_0_yyzz = primBuffer.data(k2off + 15 * idx + 12);

            auto t_0_yzzz = primBuffer.data(k2off + 15 * idx + 13);

            auto t_0_zzzz = primBuffer.data(k2off + 15 * idx + 14);

            // set up pointers to (P|G) integrals

            auto s_x_xxxx = primBuffer.data(s1off + 45 * idx);

            auto s_x_xxxy = primBuffer.data(s1off + 45 * idx + 1);

            auto s_x_xxxz = primBuffer.data(s1off + 45 * idx + 2);

            auto s_x_xxyy = primBuffer.data(s1off + 45 * idx + 3);

            auto s_x_xxyz = primBuffer.data(s1off + 45 * idx + 4);

            auto s_x_xxzz = primBuffer.data(s1off + 45 * idx + 5);

            auto s_x_xyyy = primBuffer.data(s1off + 45 * idx + 6);

            auto s_x_xyyz = primBuffer.data(s1off + 45 * idx + 7);

            auto s_x_xyzz = primBuffer.data(s1off + 45 * idx + 8);

            auto s_x_xzzz = primBuffer.data(s1off + 45 * idx + 9);

            auto s_x_yyyy = primBuffer.data(s1off + 45 * idx + 10);

            auto s_x_yyyz = primBuffer.data(s1off + 45 * idx + 11);

            auto s_x_yyzz = primBuffer.data(s1off + 45 * idx + 12);

            auto s_x_yzzz = primBuffer.data(s1off + 45 * idx + 13);

            auto s_x_zzzz = primBuffer.data(s1off + 45 * idx + 14);

            auto s_y_xxxx = primBuffer.data(s1off + 45 * idx + 15);

            auto s_y_xxxy = primBuffer.data(s1off + 45 * idx + 16);

            auto s_y_xxxz = primBuffer.data(s1off + 45 * idx + 17);

            auto s_y_xxyy = primBuffer.data(s1off + 45 * idx + 18);

            auto s_y_xxyz = primBuffer.data(s1off + 45 * idx + 19);

            auto s_y_xxzz = primBuffer.data(s1off + 45 * idx + 20);

            auto s_y_xyyy = primBuffer.data(s1off + 45 * idx + 21);

            auto s_y_xyyz = primBuffer.data(s1off + 45 * idx + 22);

            auto s_y_xyzz = primBuffer.data(s1off + 45 * idx + 23);

            auto s_y_xzzz = primBuffer.data(s1off + 45 * idx + 24);

            auto s_y_yyyy = primBuffer.data(s1off + 45 * idx + 25);

            auto s_y_yyyz = primBuffer.data(s1off + 45 * idx + 26);

            auto s_y_yyzz = primBuffer.data(s1off + 45 * idx + 27);

            auto s_y_yzzz = primBuffer.data(s1off + 45 * idx + 28);

            auto s_y_zzzz = primBuffer.data(s1off + 45 * idx + 29);

            auto s_z_xxxx = primBuffer.data(s1off + 45 * idx + 30);

            auto s_z_xxxy = primBuffer.data(s1off + 45 * idx + 31);

            auto s_z_xxxz = primBuffer.data(s1off + 45 * idx + 32);

            auto s_z_xxyy = primBuffer.data(s1off + 45 * idx + 33);

            auto s_z_xxyz = primBuffer.data(s1off + 45 * idx + 34);

            auto s_z_xxzz = primBuffer.data(s1off + 45 * idx + 35);

            auto s_z_xyyy = primBuffer.data(s1off + 45 * idx + 36);

            auto s_z_xyyz = primBuffer.data(s1off + 45 * idx + 37);

            auto s_z_xyzz = primBuffer.data(s1off + 45 * idx + 38);

            auto s_z_xzzz = primBuffer.data(s1off + 45 * idx + 39);

            auto s_z_yyyy = primBuffer.data(s1off + 45 * idx + 40);

            auto s_z_yyyz = primBuffer.data(s1off + 45 * idx + 41);

            auto s_z_yyzz = primBuffer.data(s1off + 45 * idx + 42);

            auto s_z_yzzz = primBuffer.data(s1off + 45 * idx + 43);

            auto s_z_zzzz = primBuffer.data(s1off + 45 * idx + 44);

            // set up pointers to (P|T|G) integrals

            auto t_x_xxxx = primBuffer.data(k1off + 45 * idx);

            auto t_x_xxxy = primBuffer.data(k1off + 45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(k1off + 45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(k1off + 45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(k1off + 45 * idx + 4);

            auto t_x_xxzz = primBuffer.data(k1off + 45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(k1off + 45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(k1off + 45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(k1off + 45 * idx + 8);

            auto t_x_xzzz = primBuffer.data(k1off + 45 * idx + 9);

            auto t_x_yyyy = primBuffer.data(k1off + 45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(k1off + 45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(k1off + 45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(k1off + 45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(k1off + 45 * idx + 14);

            auto t_y_xxxx = primBuffer.data(k1off + 45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(k1off + 45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(k1off + 45 * idx + 17);

            auto t_y_xxyy = primBuffer.data(k1off + 45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(k1off + 45 * idx + 19);

            auto t_y_xxzz = primBuffer.data(k1off + 45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(k1off + 45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(k1off + 45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(k1off + 45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(k1off + 45 * idx + 24);

            auto t_y_yyyy = primBuffer.data(k1off + 45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(k1off + 45 * idx + 26);

            auto t_y_yyzz = primBuffer.data(k1off + 45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(k1off + 45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(k1off + 45 * idx + 29);

            auto t_z_xxxx = primBuffer.data(k1off + 45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(k1off + 45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(k1off + 45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(k1off + 45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(k1off + 45 * idx + 34);

            auto t_z_xxzz = primBuffer.data(k1off + 45 * idx + 35);

            auto t_z_xyyy = primBuffer.data(k1off + 45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(k1off + 45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(k1off + 45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(k1off + 45 * idx + 39);

            auto t_z_yyyy = primBuffer.data(k1off + 45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(k1off + 45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(k1off + 45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(k1off + 45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(k1off + 45 * idx + 44);

            // set up pointers to (D|G) integrals

            auto s_xx_xxxx = primBuffer.data(soff + 90 * idx);

            auto s_xx_xxxy = primBuffer.data(soff + 90 * idx + 1);

            auto s_xx_xxxz = primBuffer.data(soff + 90 * idx + 2);

            auto s_xx_xxyy = primBuffer.data(soff + 90 * idx + 3);

            auto s_xx_xxyz = primBuffer.data(soff + 90 * idx + 4);

            auto s_xx_xxzz = primBuffer.data(soff + 90 * idx + 5);

            auto s_xx_xyyy = primBuffer.data(soff + 90 * idx + 6);

            auto s_xx_xyyz = primBuffer.data(soff + 90 * idx + 7);

            auto s_xx_xyzz = primBuffer.data(soff + 90 * idx + 8);

            auto s_xx_xzzz = primBuffer.data(soff + 90 * idx + 9);

            auto s_xx_yyyy = primBuffer.data(soff + 90 * idx + 10);

            auto s_xx_yyyz = primBuffer.data(soff + 90 * idx + 11);

            auto s_xx_yyzz = primBuffer.data(soff + 90 * idx + 12);

            auto s_xx_yzzz = primBuffer.data(soff + 90 * idx + 13);

            auto s_xx_zzzz = primBuffer.data(soff + 90 * idx + 14);

            auto s_xy_xxxx = primBuffer.data(soff + 90 * idx + 15);

            auto s_xy_xxxy = primBuffer.data(soff + 90 * idx + 16);

            auto s_xy_xxxz = primBuffer.data(soff + 90 * idx + 17);

            auto s_xy_xxyy = primBuffer.data(soff + 90 * idx + 18);

            auto s_xy_xxyz = primBuffer.data(soff + 90 * idx + 19);

            auto s_xy_xxzz = primBuffer.data(soff + 90 * idx + 20);

            auto s_xy_xyyy = primBuffer.data(soff + 90 * idx + 21);

            auto s_xy_xyyz = primBuffer.data(soff + 90 * idx + 22);

            auto s_xy_xyzz = primBuffer.data(soff + 90 * idx + 23);

            auto s_xy_xzzz = primBuffer.data(soff + 90 * idx + 24);

            auto s_xy_yyyy = primBuffer.data(soff + 90 * idx + 25);

            auto s_xy_yyyz = primBuffer.data(soff + 90 * idx + 26);

            auto s_xy_yyzz = primBuffer.data(soff + 90 * idx + 27);

            auto s_xy_yzzz = primBuffer.data(soff + 90 * idx + 28);

            auto s_xy_zzzz = primBuffer.data(soff + 90 * idx + 29);

            auto s_xz_xxxx = primBuffer.data(soff + 90 * idx + 30);

            auto s_xz_xxxy = primBuffer.data(soff + 90 * idx + 31);

            auto s_xz_xxxz = primBuffer.data(soff + 90 * idx + 32);

            auto s_xz_xxyy = primBuffer.data(soff + 90 * idx + 33);

            auto s_xz_xxyz = primBuffer.data(soff + 90 * idx + 34);

            auto s_xz_xxzz = primBuffer.data(soff + 90 * idx + 35);

            auto s_xz_xyyy = primBuffer.data(soff + 90 * idx + 36);

            auto s_xz_xyyz = primBuffer.data(soff + 90 * idx + 37);

            auto s_xz_xyzz = primBuffer.data(soff + 90 * idx + 38);

            auto s_xz_xzzz = primBuffer.data(soff + 90 * idx + 39);

            auto s_xz_yyyy = primBuffer.data(soff + 90 * idx + 40);

            auto s_xz_yyyz = primBuffer.data(soff + 90 * idx + 41);

            auto s_xz_yyzz = primBuffer.data(soff + 90 * idx + 42);

            auto s_xz_yzzz = primBuffer.data(soff + 90 * idx + 43);

            auto s_xz_zzzz = primBuffer.data(soff + 90 * idx + 44);

            auto s_yy_xxxx = primBuffer.data(soff + 90 * idx + 45);

            auto s_yy_xxxy = primBuffer.data(soff + 90 * idx + 46);

            auto s_yy_xxxz = primBuffer.data(soff + 90 * idx + 47);

            auto s_yy_xxyy = primBuffer.data(soff + 90 * idx + 48);

            auto s_yy_xxyz = primBuffer.data(soff + 90 * idx + 49);

            auto s_yy_xxzz = primBuffer.data(soff + 90 * idx + 50);

            auto s_yy_xyyy = primBuffer.data(soff + 90 * idx + 51);

            auto s_yy_xyyz = primBuffer.data(soff + 90 * idx + 52);

            auto s_yy_xyzz = primBuffer.data(soff + 90 * idx + 53);

            auto s_yy_xzzz = primBuffer.data(soff + 90 * idx + 54);

            auto s_yy_yyyy = primBuffer.data(soff + 90 * idx + 55);

            auto s_yy_yyyz = primBuffer.data(soff + 90 * idx + 56);

            auto s_yy_yyzz = primBuffer.data(soff + 90 * idx + 57);

            auto s_yy_yzzz = primBuffer.data(soff + 90 * idx + 58);

            auto s_yy_zzzz = primBuffer.data(soff + 90 * idx + 59);

            auto s_yz_xxxx = primBuffer.data(soff + 90 * idx + 60);

            auto s_yz_xxxy = primBuffer.data(soff + 90 * idx + 61);

            auto s_yz_xxxz = primBuffer.data(soff + 90 * idx + 62);

            auto s_yz_xxyy = primBuffer.data(soff + 90 * idx + 63);

            auto s_yz_xxyz = primBuffer.data(soff + 90 * idx + 64);

            auto s_yz_xxzz = primBuffer.data(soff + 90 * idx + 65);

            auto s_yz_xyyy = primBuffer.data(soff + 90 * idx + 66);

            auto s_yz_xyyz = primBuffer.data(soff + 90 * idx + 67);

            auto s_yz_xyzz = primBuffer.data(soff + 90 * idx + 68);

            auto s_yz_xzzz = primBuffer.data(soff + 90 * idx + 69);

            auto s_yz_yyyy = primBuffer.data(soff + 90 * idx + 70);

            auto s_yz_yyyz = primBuffer.data(soff + 90 * idx + 71);

            auto s_yz_yyzz = primBuffer.data(soff + 90 * idx + 72);

            auto s_yz_yzzz = primBuffer.data(soff + 90 * idx + 73);

            auto s_yz_zzzz = primBuffer.data(soff + 90 * idx + 74);

            auto s_zz_xxxx = primBuffer.data(soff + 90 * idx + 75);

            auto s_zz_xxxy = primBuffer.data(soff + 90 * idx + 76);

            auto s_zz_xxxz = primBuffer.data(soff + 90 * idx + 77);

            auto s_zz_xxyy = primBuffer.data(soff + 90 * idx + 78);

            auto s_zz_xxyz = primBuffer.data(soff + 90 * idx + 79);

            auto s_zz_xxzz = primBuffer.data(soff + 90 * idx + 80);

            auto s_zz_xyyy = primBuffer.data(soff + 90 * idx + 81);

            auto s_zz_xyyz = primBuffer.data(soff + 90 * idx + 82);

            auto s_zz_xyzz = primBuffer.data(soff + 90 * idx + 83);

            auto s_zz_xzzz = primBuffer.data(soff + 90 * idx + 84);

            auto s_zz_yyyy = primBuffer.data(soff + 90 * idx + 85);

            auto s_zz_yyyz = primBuffer.data(soff + 90 * idx + 86);

            auto s_zz_yyzz = primBuffer.data(soff + 90 * idx + 87);

            auto s_zz_yzzz = primBuffer.data(soff + 90 * idx + 88);

            auto s_zz_zzzz = primBuffer.data(soff + 90 * idx + 89);

            // set up pointers to (D|T|G) integrals

            auto t_xx_xxxx = primBuffer.data(koff + 90 * idx);

            auto t_xx_xxxy = primBuffer.data(koff + 90 * idx + 1);

            auto t_xx_xxxz = primBuffer.data(koff + 90 * idx + 2);

            auto t_xx_xxyy = primBuffer.data(koff + 90 * idx + 3);

            auto t_xx_xxyz = primBuffer.data(koff + 90 * idx + 4);

            auto t_xx_xxzz = primBuffer.data(koff + 90 * idx + 5);

            auto t_xx_xyyy = primBuffer.data(koff + 90 * idx + 6);

            auto t_xx_xyyz = primBuffer.data(koff + 90 * idx + 7);

            auto t_xx_xyzz = primBuffer.data(koff + 90 * idx + 8);

            auto t_xx_xzzz = primBuffer.data(koff + 90 * idx + 9);

            auto t_xx_yyyy = primBuffer.data(koff + 90 * idx + 10);

            auto t_xx_yyyz = primBuffer.data(koff + 90 * idx + 11);

            auto t_xx_yyzz = primBuffer.data(koff + 90 * idx + 12);

            auto t_xx_yzzz = primBuffer.data(koff + 90 * idx + 13);

            auto t_xx_zzzz = primBuffer.data(koff + 90 * idx + 14);

            auto t_xy_xxxx = primBuffer.data(koff + 90 * idx + 15);

            auto t_xy_xxxy = primBuffer.data(koff + 90 * idx + 16);

            auto t_xy_xxxz = primBuffer.data(koff + 90 * idx + 17);

            auto t_xy_xxyy = primBuffer.data(koff + 90 * idx + 18);

            auto t_xy_xxyz = primBuffer.data(koff + 90 * idx + 19);

            auto t_xy_xxzz = primBuffer.data(koff + 90 * idx + 20);

            auto t_xy_xyyy = primBuffer.data(koff + 90 * idx + 21);

            auto t_xy_xyyz = primBuffer.data(koff + 90 * idx + 22);

            auto t_xy_xyzz = primBuffer.data(koff + 90 * idx + 23);

            auto t_xy_xzzz = primBuffer.data(koff + 90 * idx + 24);

            auto t_xy_yyyy = primBuffer.data(koff + 90 * idx + 25);

            auto t_xy_yyyz = primBuffer.data(koff + 90 * idx + 26);

            auto t_xy_yyzz = primBuffer.data(koff + 90 * idx + 27);

            auto t_xy_yzzz = primBuffer.data(koff + 90 * idx + 28);

            auto t_xy_zzzz = primBuffer.data(koff + 90 * idx + 29);

            auto t_xz_xxxx = primBuffer.data(koff + 90 * idx + 30);

            auto t_xz_xxxy = primBuffer.data(koff + 90 * idx + 31);

            auto t_xz_xxxz = primBuffer.data(koff + 90 * idx + 32);

            auto t_xz_xxyy = primBuffer.data(koff + 90 * idx + 33);

            auto t_xz_xxyz = primBuffer.data(koff + 90 * idx + 34);

            auto t_xz_xxzz = primBuffer.data(koff + 90 * idx + 35);

            auto t_xz_xyyy = primBuffer.data(koff + 90 * idx + 36);

            auto t_xz_xyyz = primBuffer.data(koff + 90 * idx + 37);

            auto t_xz_xyzz = primBuffer.data(koff + 90 * idx + 38);

            auto t_xz_xzzz = primBuffer.data(koff + 90 * idx + 39);

            auto t_xz_yyyy = primBuffer.data(koff + 90 * idx + 40);

            auto t_xz_yyyz = primBuffer.data(koff + 90 * idx + 41);

            auto t_xz_yyzz = primBuffer.data(koff + 90 * idx + 42);

            auto t_xz_yzzz = primBuffer.data(koff + 90 * idx + 43);

            auto t_xz_zzzz = primBuffer.data(koff + 90 * idx + 44);

            auto t_yy_xxxx = primBuffer.data(koff + 90 * idx + 45);

            auto t_yy_xxxy = primBuffer.data(koff + 90 * idx + 46);

            auto t_yy_xxxz = primBuffer.data(koff + 90 * idx + 47);

            auto t_yy_xxyy = primBuffer.data(koff + 90 * idx + 48);

            auto t_yy_xxyz = primBuffer.data(koff + 90 * idx + 49);

            auto t_yy_xxzz = primBuffer.data(koff + 90 * idx + 50);

            auto t_yy_xyyy = primBuffer.data(koff + 90 * idx + 51);

            auto t_yy_xyyz = primBuffer.data(koff + 90 * idx + 52);

            auto t_yy_xyzz = primBuffer.data(koff + 90 * idx + 53);

            auto t_yy_xzzz = primBuffer.data(koff + 90 * idx + 54);

            auto t_yy_yyyy = primBuffer.data(koff + 90 * idx + 55);

            auto t_yy_yyyz = primBuffer.data(koff + 90 * idx + 56);

            auto t_yy_yyzz = primBuffer.data(koff + 90 * idx + 57);

            auto t_yy_yzzz = primBuffer.data(koff + 90 * idx + 58);

            auto t_yy_zzzz = primBuffer.data(koff + 90 * idx + 59);

            auto t_yz_xxxx = primBuffer.data(koff + 90 * idx + 60);

            auto t_yz_xxxy = primBuffer.data(koff + 90 * idx + 61);

            auto t_yz_xxxz = primBuffer.data(koff + 90 * idx + 62);

            auto t_yz_xxyy = primBuffer.data(koff + 90 * idx + 63);

            auto t_yz_xxyz = primBuffer.data(koff + 90 * idx + 64);

            auto t_yz_xxzz = primBuffer.data(koff + 90 * idx + 65);

            auto t_yz_xyyy = primBuffer.data(koff + 90 * idx + 66);

            auto t_yz_xyyz = primBuffer.data(koff + 90 * idx + 67);

            auto t_yz_xyzz = primBuffer.data(koff + 90 * idx + 68);

            auto t_yz_xzzz = primBuffer.data(koff + 90 * idx + 69);

            auto t_yz_yyyy = primBuffer.data(koff + 90 * idx + 70);

            auto t_yz_yyyz = primBuffer.data(koff + 90 * idx + 71);

            auto t_yz_yyzz = primBuffer.data(koff + 90 * idx + 72);

            auto t_yz_yzzz = primBuffer.data(koff + 90 * idx + 73);

            auto t_yz_zzzz = primBuffer.data(koff + 90 * idx + 74);

            auto t_zz_xxxx = primBuffer.data(koff + 90 * idx + 75);

            auto t_zz_xxxy = primBuffer.data(koff + 90 * idx + 76);

            auto t_zz_xxxz = primBuffer.data(koff + 90 * idx + 77);

            auto t_zz_xxyy = primBuffer.data(koff + 90 * idx + 78);

            auto t_zz_xxyz = primBuffer.data(koff + 90 * idx + 79);

            auto t_zz_xxzz = primBuffer.data(koff + 90 * idx + 80);

            auto t_zz_xyyy = primBuffer.data(koff + 90 * idx + 81);

            auto t_zz_xyyz = primBuffer.data(koff + 90 * idx + 82);

            auto t_zz_xyzz = primBuffer.data(koff + 90 * idx + 83);

            auto t_zz_xzzz = primBuffer.data(koff + 90 * idx + 84);

            auto t_zz_yyyy = primBuffer.data(koff + 90 * idx + 85);

            auto t_zz_yyyz = primBuffer.data(koff + 90 * idx + 86);

            auto t_zz_yyzz = primBuffer.data(koff + 90 * idx + 87);

            auto t_zz_yzzz = primBuffer.data(koff + 90 * idx + 88);

            auto t_zz_zzzz = primBuffer.data(koff + 90 * idx + 89);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_x_xxx, s_x_xxy,\
                                     s_x_xxz, s_x_xyy, s_x_xyz, s_x_xzz, s_x_yyy,\
                                     s_x_yyz, s_x_yzz, s_x_zzz, s_y_xxx, s_y_xxy,\
                                     s_y_xxz, s_y_xyy, s_y_xyz, s_y_xzz, s_y_yyy,\
                                     s_y_yyz, s_y_yzz, s_y_zzz, s_z_xxx, s_z_xxy,\
                                     s_z_xxz, s_z_xyy, s_z_xyz, s_z_xzz, s_z_yyy,\
                                     s_z_yyz, s_z_yzz, s_z_zzz, t_x_xxx, t_x_xxy,\
                                     t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy,\
                                     t_x_yyz, t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy,\
                                     t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy,\
                                     t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy,\
                                     t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy,\
                                     t_z_yyz, t_z_yzz, t_z_zzz, s_0_xxxx, s_0_xxxy,\
                                     s_0_xxxz, s_0_xxyy, s_0_xxyz, s_0_xxzz, s_0_xyyy,\
                                     s_0_xyyz, s_0_xyzz, s_0_xzzz, s_0_yyyy, s_0_yyyz,\
                                     s_0_yyzz, s_0_yzzz, s_0_zzzz, t_0_xxxx, t_0_xxxy,\
                                     t_0_xxxz, t_0_xxyy, t_0_xxyz, t_0_xxzz, t_0_xyyy,\
                                     t_0_xyyz, t_0_xyzz, t_0_xzzz, t_0_yyyy, t_0_yyyz,\
                                     t_0_yyzz, t_0_yzzz, t_0_zzzz, s_x_xxxx, s_x_xxxy,\
                                     s_x_xxxz, s_x_xxyy, s_x_xxyz, s_x_xxzz, s_x_xyyy,\
                                     s_x_xyyz, s_x_xyzz, s_x_xzzz, s_x_yyyy, s_x_yyyz,\
                                     s_x_yyzz, s_x_yzzz, s_x_zzzz, s_y_xxxx, s_y_xxxy,\
                                     s_y_xxxz, s_y_xxyy, s_y_xxyz, s_y_xxzz, s_y_xyyy,\
                                     s_y_xyyz, s_y_xyzz, s_y_xzzz, s_y_yyyy, s_y_yyyz,\
                                     s_y_yyzz, s_y_yzzz, s_y_zzzz, s_z_xxxx, s_z_xxxy,\
                                     s_z_xxxz, s_z_xxyy, s_z_xxyz, s_z_xxzz, s_z_xyyy,\
                                     s_z_xyyz, s_z_xyzz, s_z_xzzz, s_z_yyyy, s_z_yyyz,\
                                     s_z_yyzz, s_z_yzzz, s_z_zzzz, t_x_xxxx, t_x_xxxy,\
                                     t_x_xxxz, t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy,\
                                     t_x_xyyz, t_x_xyzz, t_x_xzzz, t_x_yyyy, t_x_yyyz,\
                                     t_x_yyzz, t_x_yzzz, t_x_zzzz, t_y_xxxx, t_y_xxxy,\
                                     t_y_xxxz, t_y_xxyy, t_y_xxyz, t_y_xxzz, t_y_xyyy,\
                                     t_y_xyyz, t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz,\
                                     t_y_yyzz, t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy,\
                                     t_z_xxxz, t_z_xxyy, t_z_xxyz, t_z_xxzz, t_z_xyyy,\
                                     t_z_xyyz, t_z_xyzz, t_z_xzzz, t_z_yyyy, t_z_yyyz,\
                                     t_z_yyzz, t_z_yzzz, t_z_zzzz, s_xx_xxxx, s_xx_xxxy,\
                                     s_xx_xxxz, s_xx_xxyy, s_xx_xxyz, s_xx_xxzz,\
                                     s_xx_xyyy, s_xx_xyyz, s_xx_xyzz, s_xx_xzzz,\
                                     s_xx_yyyy, s_xx_yyyz, s_xx_yyzz, s_xx_yzzz,\
                                     s_xx_zzzz, s_xy_xxxx, s_xy_xxxy, s_xy_xxxz,\
                                     s_xy_xxyy, s_xy_xxyz, s_xy_xxzz, s_xy_xyyy,\
                                     s_xy_xyyz, s_xy_xyzz, s_xy_xzzz, s_xy_yyyy,\
                                     s_xy_yyyz, s_xy_yyzz, s_xy_yzzz, s_xy_zzzz,\
                                     s_xz_xxxx, s_xz_xxxy, s_xz_xxxz, s_xz_xxyy,\
                                     s_xz_xxyz, s_xz_xxzz, s_xz_xyyy, s_xz_xyyz,\
                                     s_xz_xyzz, s_xz_xzzz, s_xz_yyyy, s_xz_yyyz,\
                                     s_xz_yyzz, s_xz_yzzz, s_xz_zzzz, s_yy_xxxx,\
                                     s_yy_xxxy, s_yy_xxxz, s_yy_xxyy, s_yy_xxyz,\
                                     s_yy_xxzz, s_yy_xyyy, s_yy_xyyz, s_yy_xyzz,\
                                     s_yy_xzzz, s_yy_yyyy, s_yy_yyyz, s_yy_yyzz,\
                                     s_yy_yzzz, s_yy_zzzz, s_yz_xxxx, s_yz_xxxy,\
                                     s_yz_xxxz, s_yz_xxyy, s_yz_xxyz, s_yz_xxzz,\
                                     s_yz_xyyy, s_yz_xyyz, s_yz_xyzz, s_yz_xzzz,\
                                     s_yz_yyyy, s_yz_yyyz, s_yz_yyzz, s_yz_yzzz,\
                                     s_yz_zzzz, s_zz_xxxx, s_zz_xxxy, s_zz_xxxz,\
                                     s_zz_xxyy, s_zz_xxyz, s_zz_xxzz, s_zz_xyyy,\
                                     s_zz_xyyz, s_zz_xyzz, s_zz_xzzz, s_zz_yyyy,\
                                     s_zz_yyyz, s_zz_yyzz, s_zz_yzzz, s_zz_zzzz,\
                                     t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy,\
                                     t_xx_xxyz, t_xx_xxzz, t_xx_xyyy, t_xx_xyyz,\
                                     t_xx_xyzz, t_xx_xzzz, t_xx_yyyy, t_xx_yyyz,\
                                     t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx,\
                                     t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz,\
                                     t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz,\
                                     t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz,\
                                     t_xy_yzzz, t_xy_zzzz, t_xz_xxxx, t_xz_xxxy,\
                                     t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, t_xz_xxzz,\
                                     t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz,\
                                     t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz,\
                                     t_xz_zzzz, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz,\
                                     t_yy_xxyy, t_yy_xxyz, t_yy_xxzz, t_yy_xyyy,\
                                     t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy,\
                                     t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, t_yy_zzzz,\
                                     t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, t_yz_xxyy,\
                                     t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz,\
                                     t_yz_xyzz, t_yz_xzzz, t_yz_yyyy, t_yz_yyyz,\
                                     t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_zz_xxxx,\
                                     t_zz_xxxy, t_zz_xxxz, t_zz_xxyy, t_zz_xxyz,\
                                     t_zz_xxzz, t_zz_xyyy, t_zz_xyyz, t_zz_xyzz,\
                                     t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz,\
                                     t_zz_yzzz, t_zz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xx_xxxx[j] = fr * s_x_xxxx[j] + f2t * (s_0_xxxx[j] + 4.0 * s_x_xxx[j]);

                t_xx_xxxx[j] = fr * t_x_xxxx[j] + f2t * (t_0_xxxx[j] + 4.0 * t_x_xxx[j]) + f2z * (s_xx_xxxx[j] - f2a * s_0_xxxx[j]);

                s_xx_xxxy[j] = fr * s_x_xxxy[j] + f2t * (s_0_xxxy[j] + 3.0 * s_x_xxy[j]);

                t_xx_xxxy[j] = fr * t_x_xxxy[j] + f2t * (t_0_xxxy[j] + 3.0 * t_x_xxy[j]) + f2z * (s_xx_xxxy[j] - f2a * s_0_xxxy[j]);

                s_xx_xxxz[j] = fr * s_x_xxxz[j] + f2t * (s_0_xxxz[j] + 3.0 * s_x_xxz[j]);

                t_xx_xxxz[j] = fr * t_x_xxxz[j] + f2t * (t_0_xxxz[j] + 3.0 * t_x_xxz[j]) + f2z * (s_xx_xxxz[j] - f2a * s_0_xxxz[j]);

                s_xx_xxyy[j] = fr * s_x_xxyy[j] + f2t * (s_0_xxyy[j] + 2.0 * s_x_xyy[j]);

                t_xx_xxyy[j] = fr * t_x_xxyy[j] + f2t * (t_0_xxyy[j] + 2.0 * t_x_xyy[j]) + f2z * (s_xx_xxyy[j] - f2a * s_0_xxyy[j]);

                s_xx_xxyz[j] = fr * s_x_xxyz[j] + f2t * (s_0_xxyz[j] + 2.0 * s_x_xyz[j]);

                t_xx_xxyz[j] = fr * t_x_xxyz[j] + f2t * (t_0_xxyz[j] + 2.0 * t_x_xyz[j]) + f2z * (s_xx_xxyz[j] - f2a * s_0_xxyz[j]);

                s_xx_xxzz[j] = fr * s_x_xxzz[j] + f2t * (s_0_xxzz[j] + 2.0 * s_x_xzz[j]);

                t_xx_xxzz[j] = fr * t_x_xxzz[j] + f2t * (t_0_xxzz[j] + 2.0 * t_x_xzz[j]) + f2z * (s_xx_xxzz[j] - f2a * s_0_xxzz[j]);

                s_xx_xyyy[j] = fr * s_x_xyyy[j] + f2t * (s_0_xyyy[j] + s_x_yyy[j]);

                t_xx_xyyy[j] = fr * t_x_xyyy[j] + f2t * (t_0_xyyy[j] + t_x_yyy[j]) + f2z * (s_xx_xyyy[j] - f2a * s_0_xyyy[j]);

                s_xx_xyyz[j] = fr * s_x_xyyz[j] + f2t * (s_0_xyyz[j] + s_x_yyz[j]);

                t_xx_xyyz[j] = fr * t_x_xyyz[j] + f2t * (t_0_xyyz[j] + t_x_yyz[j]) + f2z * (s_xx_xyyz[j] - f2a * s_0_xyyz[j]);

                s_xx_xyzz[j] = fr * s_x_xyzz[j] + f2t * (s_0_xyzz[j] + s_x_yzz[j]);

                t_xx_xyzz[j] = fr * t_x_xyzz[j] + f2t * (t_0_xyzz[j] + t_x_yzz[j]) + f2z * (s_xx_xyzz[j] - f2a * s_0_xyzz[j]);

                s_xx_xzzz[j] = fr * s_x_xzzz[j] + f2t * (s_0_xzzz[j] + s_x_zzz[j]);

                t_xx_xzzz[j] = fr * t_x_xzzz[j] + f2t * (t_0_xzzz[j] + t_x_zzz[j]) + f2z * (s_xx_xzzz[j] - f2a * s_0_xzzz[j]);

                s_xx_yyyy[j] = fr * s_x_yyyy[j] + f2t * s_0_yyyy[j];

                t_xx_yyyy[j] = fr * t_x_yyyy[j] + f2t * t_0_yyyy[j] + f2z * (s_xx_yyyy[j] - f2a * s_0_yyyy[j]);

                s_xx_yyyz[j] = fr * s_x_yyyz[j] + f2t * s_0_yyyz[j];

                t_xx_yyyz[j] = fr * t_x_yyyz[j] + f2t * t_0_yyyz[j] + f2z * (s_xx_yyyz[j] - f2a * s_0_yyyz[j]);

                s_xx_yyzz[j] = fr * s_x_yyzz[j] + f2t * s_0_yyzz[j];

                t_xx_yyzz[j] = fr * t_x_yyzz[j] + f2t * t_0_yyzz[j] + f2z * (s_xx_yyzz[j] - f2a * s_0_yyzz[j]);

                s_xx_yzzz[j] = fr * s_x_yzzz[j] + f2t * s_0_yzzz[j];

                t_xx_yzzz[j] = fr * t_x_yzzz[j] + f2t * t_0_yzzz[j] + f2z * (s_xx_yzzz[j] - f2a * s_0_yzzz[j]);

                s_xx_zzzz[j] = fr * s_x_zzzz[j] + f2t * s_0_zzzz[j];

                t_xx_zzzz[j] = fr * t_x_zzzz[j] + f2t * t_0_zzzz[j] + f2z * (s_xx_zzzz[j] - f2a * s_0_zzzz[j]);

                s_xy_xxxx[j] = fr * s_y_xxxx[j] + f2t * 4.0 * s_y_xxx[j];

                t_xy_xxxx[j] = fr * t_y_xxxx[j] + f2t * 4.0 * t_y_xxx[j] + f2z * s_xy_xxxx[j];

                s_xy_xxxy[j] = fr * s_y_xxxy[j] + f2t * 3.0 * s_y_xxy[j];

                t_xy_xxxy[j] = fr * t_y_xxxy[j] + f2t * 3.0 * t_y_xxy[j] + f2z * s_xy_xxxy[j];

                s_xy_xxxz[j] = fr * s_y_xxxz[j] + f2t * 3.0 * s_y_xxz[j];

                t_xy_xxxz[j] = fr * t_y_xxxz[j] + f2t * 3.0 * t_y_xxz[j] + f2z * s_xy_xxxz[j];

                s_xy_xxyy[j] = fr * s_y_xxyy[j] + f2t * 2.0 * s_y_xyy[j];

                t_xy_xxyy[j] = fr * t_y_xxyy[j] + f2t * 2.0 * t_y_xyy[j] + f2z * s_xy_xxyy[j];

                s_xy_xxyz[j] = fr * s_y_xxyz[j] + f2t * 2.0 * s_y_xyz[j];

                t_xy_xxyz[j] = fr * t_y_xxyz[j] + f2t * 2.0 * t_y_xyz[j] + f2z * s_xy_xxyz[j];

                s_xy_xxzz[j] = fr * s_y_xxzz[j] + f2t * 2.0 * s_y_xzz[j];

                t_xy_xxzz[j] = fr * t_y_xxzz[j] + f2t * 2.0 * t_y_xzz[j] + f2z * s_xy_xxzz[j];

                s_xy_xyyy[j] = fr * s_y_xyyy[j] + f2t * s_y_yyy[j];

                t_xy_xyyy[j] = fr * t_y_xyyy[j] + f2t * t_y_yyy[j] + f2z * s_xy_xyyy[j];

                s_xy_xyyz[j] = fr * s_y_xyyz[j] + f2t * s_y_yyz[j];

                t_xy_xyyz[j] = fr * t_y_xyyz[j] + f2t * t_y_yyz[j] + f2z * s_xy_xyyz[j];

                s_xy_xyzz[j] = fr * s_y_xyzz[j] + f2t * s_y_yzz[j];

                t_xy_xyzz[j] = fr * t_y_xyzz[j] + f2t * t_y_yzz[j] + f2z * s_xy_xyzz[j];

                s_xy_xzzz[j] = fr * s_y_xzzz[j] + f2t * s_y_zzz[j];

                t_xy_xzzz[j] = fr * t_y_xzzz[j] + f2t * t_y_zzz[j] + f2z * s_xy_xzzz[j];

                s_xy_yyyy[j] = fr * s_y_yyyy[j];

                t_xy_yyyy[j] = fr * t_y_yyyy[j] + f2z * s_xy_yyyy[j];

                s_xy_yyyz[j] = fr * s_y_yyyz[j];

                t_xy_yyyz[j] = fr * t_y_yyyz[j] + f2z * s_xy_yyyz[j];

                s_xy_yyzz[j] = fr * s_y_yyzz[j];

                t_xy_yyzz[j] = fr * t_y_yyzz[j] + f2z * s_xy_yyzz[j];

                s_xy_yzzz[j] = fr * s_y_yzzz[j];

                t_xy_yzzz[j] = fr * t_y_yzzz[j] + f2z * s_xy_yzzz[j];

                s_xy_zzzz[j] = fr * s_y_zzzz[j];

                t_xy_zzzz[j] = fr * t_y_zzzz[j] + f2z * s_xy_zzzz[j];

                s_xz_xxxx[j] = fr * s_z_xxxx[j] + f2t * 4.0 * s_z_xxx[j];

                t_xz_xxxx[j] = fr * t_z_xxxx[j] + f2t * 4.0 * t_z_xxx[j] + f2z * s_xz_xxxx[j];

                s_xz_xxxy[j] = fr * s_z_xxxy[j] + f2t * 3.0 * s_z_xxy[j];

                t_xz_xxxy[j] = fr * t_z_xxxy[j] + f2t * 3.0 * t_z_xxy[j] + f2z * s_xz_xxxy[j];

                s_xz_xxxz[j] = fr * s_z_xxxz[j] + f2t * 3.0 * s_z_xxz[j];

                t_xz_xxxz[j] = fr * t_z_xxxz[j] + f2t * 3.0 * t_z_xxz[j] + f2z * s_xz_xxxz[j];

                s_xz_xxyy[j] = fr * s_z_xxyy[j] + f2t * 2.0 * s_z_xyy[j];

                t_xz_xxyy[j] = fr * t_z_xxyy[j] + f2t * 2.0 * t_z_xyy[j] + f2z * s_xz_xxyy[j];

                s_xz_xxyz[j] = fr * s_z_xxyz[j] + f2t * 2.0 * s_z_xyz[j];

                t_xz_xxyz[j] = fr * t_z_xxyz[j] + f2t * 2.0 * t_z_xyz[j] + f2z * s_xz_xxyz[j];

                s_xz_xxzz[j] = fr * s_z_xxzz[j] + f2t * 2.0 * s_z_xzz[j];

                t_xz_xxzz[j] = fr * t_z_xxzz[j] + f2t * 2.0 * t_z_xzz[j] + f2z * s_xz_xxzz[j];

                s_xz_xyyy[j] = fr * s_z_xyyy[j] + f2t * s_z_yyy[j];

                t_xz_xyyy[j] = fr * t_z_xyyy[j] + f2t * t_z_yyy[j] + f2z * s_xz_xyyy[j];

                s_xz_xyyz[j] = fr * s_z_xyyz[j] + f2t * s_z_yyz[j];

                t_xz_xyyz[j] = fr * t_z_xyyz[j] + f2t * t_z_yyz[j] + f2z * s_xz_xyyz[j];

                s_xz_xyzz[j] = fr * s_z_xyzz[j] + f2t * s_z_yzz[j];

                t_xz_xyzz[j] = fr * t_z_xyzz[j] + f2t * t_z_yzz[j] + f2z * s_xz_xyzz[j];

                s_xz_xzzz[j] = fr * s_z_xzzz[j] + f2t * s_z_zzz[j];

                t_xz_xzzz[j] = fr * t_z_xzzz[j] + f2t * t_z_zzz[j] + f2z * s_xz_xzzz[j];

                s_xz_yyyy[j] = fr * s_z_yyyy[j];

                t_xz_yyyy[j] = fr * t_z_yyyy[j] + f2z * s_xz_yyyy[j];

                s_xz_yyyz[j] = fr * s_z_yyyz[j];

                t_xz_yyyz[j] = fr * t_z_yyyz[j] + f2z * s_xz_yyyz[j];

                s_xz_yyzz[j] = fr * s_z_yyzz[j];

                t_xz_yyzz[j] = fr * t_z_yyzz[j] + f2z * s_xz_yyzz[j];

                s_xz_yzzz[j] = fr * s_z_yzzz[j];

                t_xz_yzzz[j] = fr * t_z_yzzz[j] + f2z * s_xz_yzzz[j];

                s_xz_zzzz[j] = fr * s_z_zzzz[j];

                t_xz_zzzz[j] = fr * t_z_zzzz[j] + f2z * s_xz_zzzz[j];

                // leading y component

                fr = pay[j];

                s_yy_xxxx[j] = fr * s_y_xxxx[j] + f2t * s_0_xxxx[j];

                t_yy_xxxx[j] = fr * t_y_xxxx[j] + f2t * t_0_xxxx[j] + f2z * (s_yy_xxxx[j] - f2a * s_0_xxxx[j]);

                s_yy_xxxy[j] = fr * s_y_xxxy[j] + f2t * (s_0_xxxy[j] + s_y_xxx[j]);

                t_yy_xxxy[j] = fr * t_y_xxxy[j] + f2t * (t_0_xxxy[j] + t_y_xxx[j]) + f2z * (s_yy_xxxy[j] - f2a * s_0_xxxy[j]);

                s_yy_xxxz[j] = fr * s_y_xxxz[j] + f2t * s_0_xxxz[j];

                t_yy_xxxz[j] = fr * t_y_xxxz[j] + f2t * t_0_xxxz[j] + f2z * (s_yy_xxxz[j] - f2a * s_0_xxxz[j]);

                s_yy_xxyy[j] = fr * s_y_xxyy[j] + f2t * (s_0_xxyy[j] + 2.0 * s_y_xxy[j]);

                t_yy_xxyy[j] = fr * t_y_xxyy[j] + f2t * (t_0_xxyy[j] + 2.0 * t_y_xxy[j]) + f2z * (s_yy_xxyy[j] - f2a * s_0_xxyy[j]);

                s_yy_xxyz[j] = fr * s_y_xxyz[j] + f2t * (s_0_xxyz[j] + s_y_xxz[j]);

                t_yy_xxyz[j] = fr * t_y_xxyz[j] + f2t * (t_0_xxyz[j] + t_y_xxz[j]) + f2z * (s_yy_xxyz[j] - f2a * s_0_xxyz[j]);

                s_yy_xxzz[j] = fr * s_y_xxzz[j] + f2t * s_0_xxzz[j];

                t_yy_xxzz[j] = fr * t_y_xxzz[j] + f2t * t_0_xxzz[j] + f2z * (s_yy_xxzz[j] - f2a * s_0_xxzz[j]);

                s_yy_xyyy[j] = fr * s_y_xyyy[j] + f2t * (s_0_xyyy[j] + 3.0 * s_y_xyy[j]);

                t_yy_xyyy[j] = fr * t_y_xyyy[j] + f2t * (t_0_xyyy[j] + 3.0 * t_y_xyy[j]) + f2z * (s_yy_xyyy[j] - f2a * s_0_xyyy[j]);

                s_yy_xyyz[j] = fr * s_y_xyyz[j] + f2t * (s_0_xyyz[j] + 2.0 * s_y_xyz[j]);

                t_yy_xyyz[j] = fr * t_y_xyyz[j] + f2t * (t_0_xyyz[j] + 2.0 * t_y_xyz[j]) + f2z * (s_yy_xyyz[j] - f2a * s_0_xyyz[j]);

                s_yy_xyzz[j] = fr * s_y_xyzz[j] + f2t * (s_0_xyzz[j] + s_y_xzz[j]);

                t_yy_xyzz[j] = fr * t_y_xyzz[j] + f2t * (t_0_xyzz[j] + t_y_xzz[j]) + f2z * (s_yy_xyzz[j] - f2a * s_0_xyzz[j]);

                s_yy_xzzz[j] = fr * s_y_xzzz[j] + f2t * s_0_xzzz[j];

                t_yy_xzzz[j] = fr * t_y_xzzz[j] + f2t * t_0_xzzz[j] + f2z * (s_yy_xzzz[j] - f2a * s_0_xzzz[j]);

                s_yy_yyyy[j] = fr * s_y_yyyy[j] + f2t * (s_0_yyyy[j] + 4.0 * s_y_yyy[j]);

                t_yy_yyyy[j] = fr * t_y_yyyy[j] + f2t * (t_0_yyyy[j] + 4.0 * t_y_yyy[j]) + f2z * (s_yy_yyyy[j] - f2a * s_0_yyyy[j]);

                s_yy_yyyz[j] = fr * s_y_yyyz[j] + f2t * (s_0_yyyz[j] + 3.0 * s_y_yyz[j]);

                t_yy_yyyz[j] = fr * t_y_yyyz[j] + f2t * (t_0_yyyz[j] + 3.0 * t_y_yyz[j]) + f2z * (s_yy_yyyz[j] - f2a * s_0_yyyz[j]);

                s_yy_yyzz[j] = fr * s_y_yyzz[j] + f2t * (s_0_yyzz[j] + 2.0 * s_y_yzz[j]);

                t_yy_yyzz[j] = fr * t_y_yyzz[j] + f2t * (t_0_yyzz[j] + 2.0 * t_y_yzz[j]) + f2z * (s_yy_yyzz[j] - f2a * s_0_yyzz[j]);

                s_yy_yzzz[j] = fr * s_y_yzzz[j] + f2t * (s_0_yzzz[j] + s_y_zzz[j]);

                t_yy_yzzz[j] = fr * t_y_yzzz[j] + f2t * (t_0_yzzz[j] + t_y_zzz[j]) + f2z * (s_yy_yzzz[j] - f2a * s_0_yzzz[j]);

                s_yy_zzzz[j] = fr * s_y_zzzz[j] + f2t * s_0_zzzz[j];

                t_yy_zzzz[j] = fr * t_y_zzzz[j] + f2t * t_0_zzzz[j] + f2z * (s_yy_zzzz[j] - f2a * s_0_zzzz[j]);

                s_yz_xxxx[j] = fr * s_z_xxxx[j];

                t_yz_xxxx[j] = fr * t_z_xxxx[j] + f2z * s_yz_xxxx[j];

                s_yz_xxxy[j] = fr * s_z_xxxy[j] + f2t * s_z_xxx[j];

                t_yz_xxxy[j] = fr * t_z_xxxy[j] + f2t * t_z_xxx[j] + f2z * s_yz_xxxy[j];

                s_yz_xxxz[j] = fr * s_z_xxxz[j];

                t_yz_xxxz[j] = fr * t_z_xxxz[j] + f2z * s_yz_xxxz[j];

                s_yz_xxyy[j] = fr * s_z_xxyy[j] + f2t * 2.0 * s_z_xxy[j];

                t_yz_xxyy[j] = fr * t_z_xxyy[j] + f2t * 2.0 * t_z_xxy[j] + f2z * s_yz_xxyy[j];

                s_yz_xxyz[j] = fr * s_z_xxyz[j] + f2t * s_z_xxz[j];

                t_yz_xxyz[j] = fr * t_z_xxyz[j] + f2t * t_z_xxz[j] + f2z * s_yz_xxyz[j];

                s_yz_xxzz[j] = fr * s_z_xxzz[j];

                t_yz_xxzz[j] = fr * t_z_xxzz[j] + f2z * s_yz_xxzz[j];

                s_yz_xyyy[j] = fr * s_z_xyyy[j] + f2t * 3.0 * s_z_xyy[j];

                t_yz_xyyy[j] = fr * t_z_xyyy[j] + f2t * 3.0 * t_z_xyy[j] + f2z * s_yz_xyyy[j];

                s_yz_xyyz[j] = fr * s_z_xyyz[j] + f2t * 2.0 * s_z_xyz[j];

                t_yz_xyyz[j] = fr * t_z_xyyz[j] + f2t * 2.0 * t_z_xyz[j] + f2z * s_yz_xyyz[j];

                s_yz_xyzz[j] = fr * s_z_xyzz[j] + f2t * s_z_xzz[j];

                t_yz_xyzz[j] = fr * t_z_xyzz[j] + f2t * t_z_xzz[j] + f2z * s_yz_xyzz[j];

                s_yz_xzzz[j] = fr * s_z_xzzz[j];

                t_yz_xzzz[j] = fr * t_z_xzzz[j] + f2z * s_yz_xzzz[j];

                s_yz_yyyy[j] = fr * s_z_yyyy[j] + f2t * 4.0 * s_z_yyy[j];

                t_yz_yyyy[j] = fr * t_z_yyyy[j] + f2t * 4.0 * t_z_yyy[j] + f2z * s_yz_yyyy[j];

                s_yz_yyyz[j] = fr * s_z_yyyz[j] + f2t * 3.0 * s_z_yyz[j];

                t_yz_yyyz[j] = fr * t_z_yyyz[j] + f2t * 3.0 * t_z_yyz[j] + f2z * s_yz_yyyz[j];

                s_yz_yyzz[j] = fr * s_z_yyzz[j] + f2t * 2.0 * s_z_yzz[j];

                t_yz_yyzz[j] = fr * t_z_yyzz[j] + f2t * 2.0 * t_z_yzz[j] + f2z * s_yz_yyzz[j];

                s_yz_yzzz[j] = fr * s_z_yzzz[j] + f2t * s_z_zzz[j];

                t_yz_yzzz[j] = fr * t_z_yzzz[j] + f2t * t_z_zzz[j] + f2z * s_yz_yzzz[j];

                s_yz_zzzz[j] = fr * s_z_zzzz[j];

                t_yz_zzzz[j] = fr * t_z_zzzz[j] + f2z * s_yz_zzzz[j];

                // leading z component

                fr = paz[j];

                s_zz_xxxx[j] = fr * s_z_xxxx[j] + f2t * s_0_xxxx[j];

                t_zz_xxxx[j] = fr * t_z_xxxx[j] + f2t * t_0_xxxx[j] + f2z * (s_zz_xxxx[j] - f2a * s_0_xxxx[j]);

                s_zz_xxxy[j] = fr * s_z_xxxy[j] + f2t * s_0_xxxy[j];

                t_zz_xxxy[j] = fr * t_z_xxxy[j] + f2t * t_0_xxxy[j] + f2z * (s_zz_xxxy[j] - f2a * s_0_xxxy[j]);

                s_zz_xxxz[j] = fr * s_z_xxxz[j] + f2t * (s_0_xxxz[j] + s_z_xxx[j]);

                t_zz_xxxz[j] = fr * t_z_xxxz[j] + f2t * (t_0_xxxz[j] + t_z_xxx[j]) + f2z * (s_zz_xxxz[j] - f2a * s_0_xxxz[j]);

                s_zz_xxyy[j] = fr * s_z_xxyy[j] + f2t * s_0_xxyy[j];

                t_zz_xxyy[j] = fr * t_z_xxyy[j] + f2t * t_0_xxyy[j] + f2z * (s_zz_xxyy[j] - f2a * s_0_xxyy[j]);

                s_zz_xxyz[j] = fr * s_z_xxyz[j] + f2t * (s_0_xxyz[j] + s_z_xxy[j]);

                t_zz_xxyz[j] = fr * t_z_xxyz[j] + f2t * (t_0_xxyz[j] + t_z_xxy[j]) + f2z * (s_zz_xxyz[j] - f2a * s_0_xxyz[j]);

                s_zz_xxzz[j] = fr * s_z_xxzz[j] + f2t * (s_0_xxzz[j] + 2.0 * s_z_xxz[j]);

                t_zz_xxzz[j] = fr * t_z_xxzz[j] + f2t * (t_0_xxzz[j] + 2.0 * t_z_xxz[j]) + f2z * (s_zz_xxzz[j] - f2a * s_0_xxzz[j]);

                s_zz_xyyy[j] = fr * s_z_xyyy[j] + f2t * s_0_xyyy[j];

                t_zz_xyyy[j] = fr * t_z_xyyy[j] + f2t * t_0_xyyy[j] + f2z * (s_zz_xyyy[j] - f2a * s_0_xyyy[j]);

                s_zz_xyyz[j] = fr * s_z_xyyz[j] + f2t * (s_0_xyyz[j] + s_z_xyy[j]);

                t_zz_xyyz[j] = fr * t_z_xyyz[j] + f2t * (t_0_xyyz[j] + t_z_xyy[j]) + f2z * (s_zz_xyyz[j] - f2a * s_0_xyyz[j]);

                s_zz_xyzz[j] = fr * s_z_xyzz[j] + f2t * (s_0_xyzz[j] + 2.0 * s_z_xyz[j]);

                t_zz_xyzz[j] = fr * t_z_xyzz[j] + f2t * (t_0_xyzz[j] + 2.0 * t_z_xyz[j]) + f2z * (s_zz_xyzz[j] - f2a * s_0_xyzz[j]);

                s_zz_xzzz[j] = fr * s_z_xzzz[j] + f2t * (s_0_xzzz[j] + 3.0 * s_z_xzz[j]);

                t_zz_xzzz[j] = fr * t_z_xzzz[j] + f2t * (t_0_xzzz[j] + 3.0 * t_z_xzz[j]) + f2z * (s_zz_xzzz[j] - f2a * s_0_xzzz[j]);

                s_zz_yyyy[j] = fr * s_z_yyyy[j] + f2t * s_0_yyyy[j];

                t_zz_yyyy[j] = fr * t_z_yyyy[j] + f2t * t_0_yyyy[j] + f2z * (s_zz_yyyy[j] - f2a * s_0_yyyy[j]);

                s_zz_yyyz[j] = fr * s_z_yyyz[j] + f2t * (s_0_yyyz[j] + s_z_yyy[j]);

                t_zz_yyyz[j] = fr * t_z_yyyz[j] + f2t * (t_0_yyyz[j] + t_z_yyy[j]) + f2z * (s_zz_yyyz[j] - f2a * s_0_yyyz[j]);

                s_zz_yyzz[j] = fr * s_z_yyzz[j] + f2t * (s_0_yyzz[j] + 2.0 * s_z_yyz[j]);

                t_zz_yyzz[j] = fr * t_z_yyzz[j] + f2t * (t_0_yyzz[j] + 2.0 * t_z_yyz[j]) + f2z * (s_zz_yyzz[j] - f2a * s_0_yyzz[j]);

                s_zz_yzzz[j] = fr * s_z_yzzz[j] + f2t * (s_0_yzzz[j] + 3.0 * s_z_yzz[j]);

                t_zz_yzzz[j] = fr * t_z_yzzz[j] + f2t * (t_0_yzzz[j] + 3.0 * t_z_yzz[j]) + f2z * (s_zz_yzzz[j] - f2a * s_0_yzzz[j]);

                s_zz_zzzz[j] = fr * s_z_zzzz[j] + f2t * (s_0_zzzz[j] + 4.0 * s_z_zzz[j]);

                t_zz_zzzz[j] = fr * t_z_zzzz[j] + f2t * (t_0_zzzz[j] + 4.0 * t_z_zzz[j]) + f2z * (s_zz_zzzz[j] - f2a * s_0_zzzz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGD(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 2, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 2, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|P) integrals

            auto s_xxx_x = primBuffer.data(skoff + 30 * idx);

            auto s_xxx_y = primBuffer.data(skoff + 30 * idx + 1);

            auto s_xxx_z = primBuffer.data(skoff + 30 * idx + 2);

            auto s_xxy_x = primBuffer.data(skoff + 30 * idx + 3);

            auto s_xxy_y = primBuffer.data(skoff + 30 * idx + 4);

            auto s_xxy_z = primBuffer.data(skoff + 30 * idx + 5);

            auto s_xxz_x = primBuffer.data(skoff + 30 * idx + 6);

            auto s_xxz_y = primBuffer.data(skoff + 30 * idx + 7);

            auto s_xxz_z = primBuffer.data(skoff + 30 * idx + 8);

            auto s_xyy_x = primBuffer.data(skoff + 30 * idx + 9);

            auto s_xyy_y = primBuffer.data(skoff + 30 * idx + 10);

            auto s_xyy_z = primBuffer.data(skoff + 30 * idx + 11);

            auto s_xyz_x = primBuffer.data(skoff + 30 * idx + 12);

            auto s_xyz_y = primBuffer.data(skoff + 30 * idx + 13);

            auto s_xyz_z = primBuffer.data(skoff + 30 * idx + 14);

            auto s_xzz_x = primBuffer.data(skoff + 30 * idx + 15);

            auto s_xzz_y = primBuffer.data(skoff + 30 * idx + 16);

            auto s_xzz_z = primBuffer.data(skoff + 30 * idx + 17);

            auto s_yyy_x = primBuffer.data(skoff + 30 * idx + 18);

            auto s_yyy_y = primBuffer.data(skoff + 30 * idx + 19);

            auto s_yyy_z = primBuffer.data(skoff + 30 * idx + 20);

            auto s_yyz_x = primBuffer.data(skoff + 30 * idx + 21);

            auto s_yyz_y = primBuffer.data(skoff + 30 * idx + 22);

            auto s_yyz_z = primBuffer.data(skoff + 30 * idx + 23);

            auto s_yzz_x = primBuffer.data(skoff + 30 * idx + 24);

            auto s_yzz_y = primBuffer.data(skoff + 30 * idx + 25);

            auto s_yzz_z = primBuffer.data(skoff + 30 * idx + 26);

            auto s_zzz_x = primBuffer.data(skoff + 30 * idx + 27);

            auto s_zzz_y = primBuffer.data(skoff + 30 * idx + 28);

            auto s_zzz_z = primBuffer.data(skoff + 30 * idx + 29);

            // set up pointers to (F|T|P) integrals

            auto t_xxx_x = primBuffer.data(kkoff + 30 * idx);

            auto t_xxx_y = primBuffer.data(kkoff + 30 * idx + 1);

            auto t_xxx_z = primBuffer.data(kkoff + 30 * idx + 2);

            auto t_xxy_x = primBuffer.data(kkoff + 30 * idx + 3);

            auto t_xxy_y = primBuffer.data(kkoff + 30 * idx + 4);

            auto t_xxy_z = primBuffer.data(kkoff + 30 * idx + 5);

            auto t_xxz_x = primBuffer.data(kkoff + 30 * idx + 6);

            auto t_xxz_y = primBuffer.data(kkoff + 30 * idx + 7);

            auto t_xxz_z = primBuffer.data(kkoff + 30 * idx + 8);

            auto t_xyy_x = primBuffer.data(kkoff + 30 * idx + 9);

            auto t_xyy_y = primBuffer.data(kkoff + 30 * idx + 10);

            auto t_xyy_z = primBuffer.data(kkoff + 30 * idx + 11);

            auto t_xyz_x = primBuffer.data(kkoff + 30 * idx + 12);

            auto t_xyz_y = primBuffer.data(kkoff + 30 * idx + 13);

            auto t_xyz_z = primBuffer.data(kkoff + 30 * idx + 14);

            auto t_xzz_x = primBuffer.data(kkoff + 30 * idx + 15);

            auto t_xzz_y = primBuffer.data(kkoff + 30 * idx + 16);

            auto t_xzz_z = primBuffer.data(kkoff + 30 * idx + 17);

            auto t_yyy_x = primBuffer.data(kkoff + 30 * idx + 18);

            auto t_yyy_y = primBuffer.data(kkoff + 30 * idx + 19);

            auto t_yyy_z = primBuffer.data(kkoff + 30 * idx + 20);

            auto t_yyz_x = primBuffer.data(kkoff + 30 * idx + 21);

            auto t_yyz_y = primBuffer.data(kkoff + 30 * idx + 22);

            auto t_yyz_z = primBuffer.data(kkoff + 30 * idx + 23);

            auto t_yzz_x = primBuffer.data(kkoff + 30 * idx + 24);

            auto t_yzz_y = primBuffer.data(kkoff + 30 * idx + 25);

            auto t_yzz_z = primBuffer.data(kkoff + 30 * idx + 26);

            auto t_zzz_x = primBuffer.data(kkoff + 30 * idx + 27);

            auto t_zzz_y = primBuffer.data(kkoff + 30 * idx + 28);

            auto t_zzz_z = primBuffer.data(kkoff + 30 * idx + 29);

            // set up pointers to (D|D) integrals

            auto s_xx_xx = primBuffer.data(s2off + 36 * idx);

            auto s_xx_xy = primBuffer.data(s2off + 36 * idx + 1);

            auto s_xx_xz = primBuffer.data(s2off + 36 * idx + 2);

            auto s_xx_yy = primBuffer.data(s2off + 36 * idx + 3);

            auto s_xx_yz = primBuffer.data(s2off + 36 * idx + 4);

            auto s_xx_zz = primBuffer.data(s2off + 36 * idx + 5);

            auto s_xy_xx = primBuffer.data(s2off + 36 * idx + 6);

            auto s_xy_xy = primBuffer.data(s2off + 36 * idx + 7);

            auto s_xy_xz = primBuffer.data(s2off + 36 * idx + 8);

            auto s_xy_yy = primBuffer.data(s2off + 36 * idx + 9);

            auto s_xy_yz = primBuffer.data(s2off + 36 * idx + 10);

            auto s_xy_zz = primBuffer.data(s2off + 36 * idx + 11);

            auto s_xz_xx = primBuffer.data(s2off + 36 * idx + 12);

            auto s_xz_xy = primBuffer.data(s2off + 36 * idx + 13);

            auto s_xz_xz = primBuffer.data(s2off + 36 * idx + 14);

            auto s_xz_yy = primBuffer.data(s2off + 36 * idx + 15);

            auto s_xz_yz = primBuffer.data(s2off + 36 * idx + 16);

            auto s_xz_zz = primBuffer.data(s2off + 36 * idx + 17);

            auto s_yy_xx = primBuffer.data(s2off + 36 * idx + 18);

            auto s_yy_xy = primBuffer.data(s2off + 36 * idx + 19);

            auto s_yy_xz = primBuffer.data(s2off + 36 * idx + 20);

            auto s_yy_yy = primBuffer.data(s2off + 36 * idx + 21);

            auto s_yy_yz = primBuffer.data(s2off + 36 * idx + 22);

            auto s_yy_zz = primBuffer.data(s2off + 36 * idx + 23);

            auto s_yz_xx = primBuffer.data(s2off + 36 * idx + 24);

            auto s_yz_xy = primBuffer.data(s2off + 36 * idx + 25);

            auto s_yz_xz = primBuffer.data(s2off + 36 * idx + 26);

            auto s_yz_yy = primBuffer.data(s2off + 36 * idx + 27);

            auto s_yz_yz = primBuffer.data(s2off + 36 * idx + 28);

            auto s_yz_zz = primBuffer.data(s2off + 36 * idx + 29);

            auto s_zz_xx = primBuffer.data(s2off + 36 * idx + 30);

            auto s_zz_xy = primBuffer.data(s2off + 36 * idx + 31);

            auto s_zz_xz = primBuffer.data(s2off + 36 * idx + 32);

            auto s_zz_yy = primBuffer.data(s2off + 36 * idx + 33);

            auto s_zz_yz = primBuffer.data(s2off + 36 * idx + 34);

            auto s_zz_zz = primBuffer.data(s2off + 36 * idx + 35);

            // set up pointers to (D|T|D) integrals

            auto t_xx_xx = primBuffer.data(k2off + 36 * idx);

            auto t_xx_xy = primBuffer.data(k2off + 36 * idx + 1);

            auto t_xx_xz = primBuffer.data(k2off + 36 * idx + 2);

            auto t_xx_yy = primBuffer.data(k2off + 36 * idx + 3);

            auto t_xx_yz = primBuffer.data(k2off + 36 * idx + 4);

            auto t_xx_zz = primBuffer.data(k2off + 36 * idx + 5);

            auto t_xy_xx = primBuffer.data(k2off + 36 * idx + 6);

            auto t_xy_xy = primBuffer.data(k2off + 36 * idx + 7);

            auto t_xy_xz = primBuffer.data(k2off + 36 * idx + 8);

            auto t_xy_yy = primBuffer.data(k2off + 36 * idx + 9);

            auto t_xy_yz = primBuffer.data(k2off + 36 * idx + 10);

            auto t_xy_zz = primBuffer.data(k2off + 36 * idx + 11);

            auto t_xz_xx = primBuffer.data(k2off + 36 * idx + 12);

            auto t_xz_xy = primBuffer.data(k2off + 36 * idx + 13);

            auto t_xz_xz = primBuffer.data(k2off + 36 * idx + 14);

            auto t_xz_yy = primBuffer.data(k2off + 36 * idx + 15);

            auto t_xz_yz = primBuffer.data(k2off + 36 * idx + 16);

            auto t_xz_zz = primBuffer.data(k2off + 36 * idx + 17);

            auto t_yy_xx = primBuffer.data(k2off + 36 * idx + 18);

            auto t_yy_xy = primBuffer.data(k2off + 36 * idx + 19);

            auto t_yy_xz = primBuffer.data(k2off + 36 * idx + 20);

            auto t_yy_yy = primBuffer.data(k2off + 36 * idx + 21);

            auto t_yy_yz = primBuffer.data(k2off + 36 * idx + 22);

            auto t_yy_zz = primBuffer.data(k2off + 36 * idx + 23);

            auto t_yz_xx = primBuffer.data(k2off + 36 * idx + 24);

            auto t_yz_xy = primBuffer.data(k2off + 36 * idx + 25);

            auto t_yz_xz = primBuffer.data(k2off + 36 * idx + 26);

            auto t_yz_yy = primBuffer.data(k2off + 36 * idx + 27);

            auto t_yz_yz = primBuffer.data(k2off + 36 * idx + 28);

            auto t_yz_zz = primBuffer.data(k2off + 36 * idx + 29);

            auto t_zz_xx = primBuffer.data(k2off + 36 * idx + 30);

            auto t_zz_xy = primBuffer.data(k2off + 36 * idx + 31);

            auto t_zz_xz = primBuffer.data(k2off + 36 * idx + 32);

            auto t_zz_yy = primBuffer.data(k2off + 36 * idx + 33);

            auto t_zz_yz = primBuffer.data(k2off + 36 * idx + 34);

            auto t_zz_zz = primBuffer.data(k2off + 36 * idx + 35);

            // set up pointers to (F|D) integrals

            auto s_xxx_xx = primBuffer.data(s1off + 60 * idx);

            auto s_xxx_xy = primBuffer.data(s1off + 60 * idx + 1);

            auto s_xxx_xz = primBuffer.data(s1off + 60 * idx + 2);

            auto s_xxx_yy = primBuffer.data(s1off + 60 * idx + 3);

            auto s_xxx_yz = primBuffer.data(s1off + 60 * idx + 4);

            auto s_xxx_zz = primBuffer.data(s1off + 60 * idx + 5);

            auto s_xxy_xx = primBuffer.data(s1off + 60 * idx + 6);

            auto s_xxy_xy = primBuffer.data(s1off + 60 * idx + 7);

            auto s_xxy_xz = primBuffer.data(s1off + 60 * idx + 8);

            auto s_xxy_yy = primBuffer.data(s1off + 60 * idx + 9);

            auto s_xxy_yz = primBuffer.data(s1off + 60 * idx + 10);

            auto s_xxy_zz = primBuffer.data(s1off + 60 * idx + 11);

            auto s_xxz_xx = primBuffer.data(s1off + 60 * idx + 12);

            auto s_xxz_xy = primBuffer.data(s1off + 60 * idx + 13);

            auto s_xxz_xz = primBuffer.data(s1off + 60 * idx + 14);

            auto s_xxz_yy = primBuffer.data(s1off + 60 * idx + 15);

            auto s_xxz_yz = primBuffer.data(s1off + 60 * idx + 16);

            auto s_xxz_zz = primBuffer.data(s1off + 60 * idx + 17);

            auto s_xyy_xx = primBuffer.data(s1off + 60 * idx + 18);

            auto s_xyy_xy = primBuffer.data(s1off + 60 * idx + 19);

            auto s_xyy_xz = primBuffer.data(s1off + 60 * idx + 20);

            auto s_xyy_yy = primBuffer.data(s1off + 60 * idx + 21);

            auto s_xyy_yz = primBuffer.data(s1off + 60 * idx + 22);

            auto s_xyy_zz = primBuffer.data(s1off + 60 * idx + 23);

            auto s_xyz_xx = primBuffer.data(s1off + 60 * idx + 24);

            auto s_xyz_xy = primBuffer.data(s1off + 60 * idx + 25);

            auto s_xyz_xz = primBuffer.data(s1off + 60 * idx + 26);

            auto s_xyz_yy = primBuffer.data(s1off + 60 * idx + 27);

            auto s_xyz_yz = primBuffer.data(s1off + 60 * idx + 28);

            auto s_xyz_zz = primBuffer.data(s1off + 60 * idx + 29);

            auto s_xzz_xx = primBuffer.data(s1off + 60 * idx + 30);

            auto s_xzz_xy = primBuffer.data(s1off + 60 * idx + 31);

            auto s_xzz_xz = primBuffer.data(s1off + 60 * idx + 32);

            auto s_xzz_yy = primBuffer.data(s1off + 60 * idx + 33);

            auto s_xzz_yz = primBuffer.data(s1off + 60 * idx + 34);

            auto s_xzz_zz = primBuffer.data(s1off + 60 * idx + 35);

            auto s_yyy_xx = primBuffer.data(s1off + 60 * idx + 36);

            auto s_yyy_xy = primBuffer.data(s1off + 60 * idx + 37);

            auto s_yyy_xz = primBuffer.data(s1off + 60 * idx + 38);

            auto s_yyy_yy = primBuffer.data(s1off + 60 * idx + 39);

            auto s_yyy_yz = primBuffer.data(s1off + 60 * idx + 40);

            auto s_yyy_zz = primBuffer.data(s1off + 60 * idx + 41);

            auto s_yyz_xx = primBuffer.data(s1off + 60 * idx + 42);

            auto s_yyz_xy = primBuffer.data(s1off + 60 * idx + 43);

            auto s_yyz_xz = primBuffer.data(s1off + 60 * idx + 44);

            auto s_yyz_yy = primBuffer.data(s1off + 60 * idx + 45);

            auto s_yyz_yz = primBuffer.data(s1off + 60 * idx + 46);

            auto s_yyz_zz = primBuffer.data(s1off + 60 * idx + 47);

            auto s_yzz_xx = primBuffer.data(s1off + 60 * idx + 48);

            auto s_yzz_xy = primBuffer.data(s1off + 60 * idx + 49);

            auto s_yzz_xz = primBuffer.data(s1off + 60 * idx + 50);

            auto s_yzz_yy = primBuffer.data(s1off + 60 * idx + 51);

            auto s_yzz_yz = primBuffer.data(s1off + 60 * idx + 52);

            auto s_yzz_zz = primBuffer.data(s1off + 60 * idx + 53);

            auto s_zzz_xx = primBuffer.data(s1off + 60 * idx + 54);

            auto s_zzz_xy = primBuffer.data(s1off + 60 * idx + 55);

            auto s_zzz_xz = primBuffer.data(s1off + 60 * idx + 56);

            auto s_zzz_yy = primBuffer.data(s1off + 60 * idx + 57);

            auto s_zzz_yz = primBuffer.data(s1off + 60 * idx + 58);

            auto s_zzz_zz = primBuffer.data(s1off + 60 * idx + 59);

            // set up pointers to (F|T|D) integrals

            auto t_xxx_xx = primBuffer.data(k1off + 60 * idx);

            auto t_xxx_xy = primBuffer.data(k1off + 60 * idx + 1);

            auto t_xxx_xz = primBuffer.data(k1off + 60 * idx + 2);

            auto t_xxx_yy = primBuffer.data(k1off + 60 * idx + 3);

            auto t_xxx_yz = primBuffer.data(k1off + 60 * idx + 4);

            auto t_xxx_zz = primBuffer.data(k1off + 60 * idx + 5);

            auto t_xxy_xx = primBuffer.data(k1off + 60 * idx + 6);

            auto t_xxy_xy = primBuffer.data(k1off + 60 * idx + 7);

            auto t_xxy_xz = primBuffer.data(k1off + 60 * idx + 8);

            auto t_xxy_yy = primBuffer.data(k1off + 60 * idx + 9);

            auto t_xxy_yz = primBuffer.data(k1off + 60 * idx + 10);

            auto t_xxy_zz = primBuffer.data(k1off + 60 * idx + 11);

            auto t_xxz_xx = primBuffer.data(k1off + 60 * idx + 12);

            auto t_xxz_xy = primBuffer.data(k1off + 60 * idx + 13);

            auto t_xxz_xz = primBuffer.data(k1off + 60 * idx + 14);

            auto t_xxz_yy = primBuffer.data(k1off + 60 * idx + 15);

            auto t_xxz_yz = primBuffer.data(k1off + 60 * idx + 16);

            auto t_xxz_zz = primBuffer.data(k1off + 60 * idx + 17);

            auto t_xyy_xx = primBuffer.data(k1off + 60 * idx + 18);

            auto t_xyy_xy = primBuffer.data(k1off + 60 * idx + 19);

            auto t_xyy_xz = primBuffer.data(k1off + 60 * idx + 20);

            auto t_xyy_yy = primBuffer.data(k1off + 60 * idx + 21);

            auto t_xyy_yz = primBuffer.data(k1off + 60 * idx + 22);

            auto t_xyy_zz = primBuffer.data(k1off + 60 * idx + 23);

            auto t_xyz_xx = primBuffer.data(k1off + 60 * idx + 24);

            auto t_xyz_xy = primBuffer.data(k1off + 60 * idx + 25);

            auto t_xyz_xz = primBuffer.data(k1off + 60 * idx + 26);

            auto t_xyz_yy = primBuffer.data(k1off + 60 * idx + 27);

            auto t_xyz_yz = primBuffer.data(k1off + 60 * idx + 28);

            auto t_xyz_zz = primBuffer.data(k1off + 60 * idx + 29);

            auto t_xzz_xx = primBuffer.data(k1off + 60 * idx + 30);

            auto t_xzz_xy = primBuffer.data(k1off + 60 * idx + 31);

            auto t_xzz_xz = primBuffer.data(k1off + 60 * idx + 32);

            auto t_xzz_yy = primBuffer.data(k1off + 60 * idx + 33);

            auto t_xzz_yz = primBuffer.data(k1off + 60 * idx + 34);

            auto t_xzz_zz = primBuffer.data(k1off + 60 * idx + 35);

            auto t_yyy_xx = primBuffer.data(k1off + 60 * idx + 36);

            auto t_yyy_xy = primBuffer.data(k1off + 60 * idx + 37);

            auto t_yyy_xz = primBuffer.data(k1off + 60 * idx + 38);

            auto t_yyy_yy = primBuffer.data(k1off + 60 * idx + 39);

            auto t_yyy_yz = primBuffer.data(k1off + 60 * idx + 40);

            auto t_yyy_zz = primBuffer.data(k1off + 60 * idx + 41);

            auto t_yyz_xx = primBuffer.data(k1off + 60 * idx + 42);

            auto t_yyz_xy = primBuffer.data(k1off + 60 * idx + 43);

            auto t_yyz_xz = primBuffer.data(k1off + 60 * idx + 44);

            auto t_yyz_yy = primBuffer.data(k1off + 60 * idx + 45);

            auto t_yyz_yz = primBuffer.data(k1off + 60 * idx + 46);

            auto t_yyz_zz = primBuffer.data(k1off + 60 * idx + 47);

            auto t_yzz_xx = primBuffer.data(k1off + 60 * idx + 48);

            auto t_yzz_xy = primBuffer.data(k1off + 60 * idx + 49);

            auto t_yzz_xz = primBuffer.data(k1off + 60 * idx + 50);

            auto t_yzz_yy = primBuffer.data(k1off + 60 * idx + 51);

            auto t_yzz_yz = primBuffer.data(k1off + 60 * idx + 52);

            auto t_yzz_zz = primBuffer.data(k1off + 60 * idx + 53);

            auto t_zzz_xx = primBuffer.data(k1off + 60 * idx + 54);

            auto t_zzz_xy = primBuffer.data(k1off + 60 * idx + 55);

            auto t_zzz_xz = primBuffer.data(k1off + 60 * idx + 56);

            auto t_zzz_yy = primBuffer.data(k1off + 60 * idx + 57);

            auto t_zzz_yz = primBuffer.data(k1off + 60 * idx + 58);

            auto t_zzz_zz = primBuffer.data(k1off + 60 * idx + 59);

            // set up pointers to (G|D) integrals

            auto s_xxxx_xx = primBuffer.data(soff + 90 * idx);

            auto s_xxxx_xy = primBuffer.data(soff + 90 * idx + 1);

            auto s_xxxx_xz = primBuffer.data(soff + 90 * idx + 2);

            auto s_xxxx_yy = primBuffer.data(soff + 90 * idx + 3);

            auto s_xxxx_yz = primBuffer.data(soff + 90 * idx + 4);

            auto s_xxxx_zz = primBuffer.data(soff + 90 * idx + 5);

            auto s_xxxy_xx = primBuffer.data(soff + 90 * idx + 6);

            auto s_xxxy_xy = primBuffer.data(soff + 90 * idx + 7);

            auto s_xxxy_xz = primBuffer.data(soff + 90 * idx + 8);

            auto s_xxxy_yy = primBuffer.data(soff + 90 * idx + 9);

            auto s_xxxy_yz = primBuffer.data(soff + 90 * idx + 10);

            auto s_xxxy_zz = primBuffer.data(soff + 90 * idx + 11);

            auto s_xxxz_xx = primBuffer.data(soff + 90 * idx + 12);

            auto s_xxxz_xy = primBuffer.data(soff + 90 * idx + 13);

            auto s_xxxz_xz = primBuffer.data(soff + 90 * idx + 14);

            auto s_xxxz_yy = primBuffer.data(soff + 90 * idx + 15);

            auto s_xxxz_yz = primBuffer.data(soff + 90 * idx + 16);

            auto s_xxxz_zz = primBuffer.data(soff + 90 * idx + 17);

            auto s_xxyy_xx = primBuffer.data(soff + 90 * idx + 18);

            auto s_xxyy_xy = primBuffer.data(soff + 90 * idx + 19);

            auto s_xxyy_xz = primBuffer.data(soff + 90 * idx + 20);

            auto s_xxyy_yy = primBuffer.data(soff + 90 * idx + 21);

            auto s_xxyy_yz = primBuffer.data(soff + 90 * idx + 22);

            auto s_xxyy_zz = primBuffer.data(soff + 90 * idx + 23);

            auto s_xxyz_xx = primBuffer.data(soff + 90 * idx + 24);

            auto s_xxyz_xy = primBuffer.data(soff + 90 * idx + 25);

            auto s_xxyz_xz = primBuffer.data(soff + 90 * idx + 26);

            auto s_xxyz_yy = primBuffer.data(soff + 90 * idx + 27);

            auto s_xxyz_yz = primBuffer.data(soff + 90 * idx + 28);

            auto s_xxyz_zz = primBuffer.data(soff + 90 * idx + 29);

            auto s_xxzz_xx = primBuffer.data(soff + 90 * idx + 30);

            auto s_xxzz_xy = primBuffer.data(soff + 90 * idx + 31);

            auto s_xxzz_xz = primBuffer.data(soff + 90 * idx + 32);

            auto s_xxzz_yy = primBuffer.data(soff + 90 * idx + 33);

            auto s_xxzz_yz = primBuffer.data(soff + 90 * idx + 34);

            auto s_xxzz_zz = primBuffer.data(soff + 90 * idx + 35);

            auto s_xyyy_xx = primBuffer.data(soff + 90 * idx + 36);

            auto s_xyyy_xy = primBuffer.data(soff + 90 * idx + 37);

            auto s_xyyy_xz = primBuffer.data(soff + 90 * idx + 38);

            auto s_xyyy_yy = primBuffer.data(soff + 90 * idx + 39);

            auto s_xyyy_yz = primBuffer.data(soff + 90 * idx + 40);

            auto s_xyyy_zz = primBuffer.data(soff + 90 * idx + 41);

            auto s_xyyz_xx = primBuffer.data(soff + 90 * idx + 42);

            auto s_xyyz_xy = primBuffer.data(soff + 90 * idx + 43);

            auto s_xyyz_xz = primBuffer.data(soff + 90 * idx + 44);

            auto s_xyyz_yy = primBuffer.data(soff + 90 * idx + 45);

            auto s_xyyz_yz = primBuffer.data(soff + 90 * idx + 46);

            auto s_xyyz_zz = primBuffer.data(soff + 90 * idx + 47);

            auto s_xyzz_xx = primBuffer.data(soff + 90 * idx + 48);

            auto s_xyzz_xy = primBuffer.data(soff + 90 * idx + 49);

            auto s_xyzz_xz = primBuffer.data(soff + 90 * idx + 50);

            auto s_xyzz_yy = primBuffer.data(soff + 90 * idx + 51);

            auto s_xyzz_yz = primBuffer.data(soff + 90 * idx + 52);

            auto s_xyzz_zz = primBuffer.data(soff + 90 * idx + 53);

            auto s_xzzz_xx = primBuffer.data(soff + 90 * idx + 54);

            auto s_xzzz_xy = primBuffer.data(soff + 90 * idx + 55);

            auto s_xzzz_xz = primBuffer.data(soff + 90 * idx + 56);

            auto s_xzzz_yy = primBuffer.data(soff + 90 * idx + 57);

            auto s_xzzz_yz = primBuffer.data(soff + 90 * idx + 58);

            auto s_xzzz_zz = primBuffer.data(soff + 90 * idx + 59);

            auto s_yyyy_xx = primBuffer.data(soff + 90 * idx + 60);

            auto s_yyyy_xy = primBuffer.data(soff + 90 * idx + 61);

            auto s_yyyy_xz = primBuffer.data(soff + 90 * idx + 62);

            auto s_yyyy_yy = primBuffer.data(soff + 90 * idx + 63);

            auto s_yyyy_yz = primBuffer.data(soff + 90 * idx + 64);

            auto s_yyyy_zz = primBuffer.data(soff + 90 * idx + 65);

            auto s_yyyz_xx = primBuffer.data(soff + 90 * idx + 66);

            auto s_yyyz_xy = primBuffer.data(soff + 90 * idx + 67);

            auto s_yyyz_xz = primBuffer.data(soff + 90 * idx + 68);

            auto s_yyyz_yy = primBuffer.data(soff + 90 * idx + 69);

            auto s_yyyz_yz = primBuffer.data(soff + 90 * idx + 70);

            auto s_yyyz_zz = primBuffer.data(soff + 90 * idx + 71);

            auto s_yyzz_xx = primBuffer.data(soff + 90 * idx + 72);

            auto s_yyzz_xy = primBuffer.data(soff + 90 * idx + 73);

            auto s_yyzz_xz = primBuffer.data(soff + 90 * idx + 74);

            auto s_yyzz_yy = primBuffer.data(soff + 90 * idx + 75);

            auto s_yyzz_yz = primBuffer.data(soff + 90 * idx + 76);

            auto s_yyzz_zz = primBuffer.data(soff + 90 * idx + 77);

            auto s_yzzz_xx = primBuffer.data(soff + 90 * idx + 78);

            auto s_yzzz_xy = primBuffer.data(soff + 90 * idx + 79);

            auto s_yzzz_xz = primBuffer.data(soff + 90 * idx + 80);

            auto s_yzzz_yy = primBuffer.data(soff + 90 * idx + 81);

            auto s_yzzz_yz = primBuffer.data(soff + 90 * idx + 82);

            auto s_yzzz_zz = primBuffer.data(soff + 90 * idx + 83);

            auto s_zzzz_xx = primBuffer.data(soff + 90 * idx + 84);

            auto s_zzzz_xy = primBuffer.data(soff + 90 * idx + 85);

            auto s_zzzz_xz = primBuffer.data(soff + 90 * idx + 86);

            auto s_zzzz_yy = primBuffer.data(soff + 90 * idx + 87);

            auto s_zzzz_yz = primBuffer.data(soff + 90 * idx + 88);

            auto s_zzzz_zz = primBuffer.data(soff + 90 * idx + 89);

            // set up pointers to (G|T|D) integrals

            auto t_xxxx_xx = primBuffer.data(koff + 90 * idx);

            auto t_xxxx_xy = primBuffer.data(koff + 90 * idx + 1);

            auto t_xxxx_xz = primBuffer.data(koff + 90 * idx + 2);

            auto t_xxxx_yy = primBuffer.data(koff + 90 * idx + 3);

            auto t_xxxx_yz = primBuffer.data(koff + 90 * idx + 4);

            auto t_xxxx_zz = primBuffer.data(koff + 90 * idx + 5);

            auto t_xxxy_xx = primBuffer.data(koff + 90 * idx + 6);

            auto t_xxxy_xy = primBuffer.data(koff + 90 * idx + 7);

            auto t_xxxy_xz = primBuffer.data(koff + 90 * idx + 8);

            auto t_xxxy_yy = primBuffer.data(koff + 90 * idx + 9);

            auto t_xxxy_yz = primBuffer.data(koff + 90 * idx + 10);

            auto t_xxxy_zz = primBuffer.data(koff + 90 * idx + 11);

            auto t_xxxz_xx = primBuffer.data(koff + 90 * idx + 12);

            auto t_xxxz_xy = primBuffer.data(koff + 90 * idx + 13);

            auto t_xxxz_xz = primBuffer.data(koff + 90 * idx + 14);

            auto t_xxxz_yy = primBuffer.data(koff + 90 * idx + 15);

            auto t_xxxz_yz = primBuffer.data(koff + 90 * idx + 16);

            auto t_xxxz_zz = primBuffer.data(koff + 90 * idx + 17);

            auto t_xxyy_xx = primBuffer.data(koff + 90 * idx + 18);

            auto t_xxyy_xy = primBuffer.data(koff + 90 * idx + 19);

            auto t_xxyy_xz = primBuffer.data(koff + 90 * idx + 20);

            auto t_xxyy_yy = primBuffer.data(koff + 90 * idx + 21);

            auto t_xxyy_yz = primBuffer.data(koff + 90 * idx + 22);

            auto t_xxyy_zz = primBuffer.data(koff + 90 * idx + 23);

            auto t_xxyz_xx = primBuffer.data(koff + 90 * idx + 24);

            auto t_xxyz_xy = primBuffer.data(koff + 90 * idx + 25);

            auto t_xxyz_xz = primBuffer.data(koff + 90 * idx + 26);

            auto t_xxyz_yy = primBuffer.data(koff + 90 * idx + 27);

            auto t_xxyz_yz = primBuffer.data(koff + 90 * idx + 28);

            auto t_xxyz_zz = primBuffer.data(koff + 90 * idx + 29);

            auto t_xxzz_xx = primBuffer.data(koff + 90 * idx + 30);

            auto t_xxzz_xy = primBuffer.data(koff + 90 * idx + 31);

            auto t_xxzz_xz = primBuffer.data(koff + 90 * idx + 32);

            auto t_xxzz_yy = primBuffer.data(koff + 90 * idx + 33);

            auto t_xxzz_yz = primBuffer.data(koff + 90 * idx + 34);

            auto t_xxzz_zz = primBuffer.data(koff + 90 * idx + 35);

            auto t_xyyy_xx = primBuffer.data(koff + 90 * idx + 36);

            auto t_xyyy_xy = primBuffer.data(koff + 90 * idx + 37);

            auto t_xyyy_xz = primBuffer.data(koff + 90 * idx + 38);

            auto t_xyyy_yy = primBuffer.data(koff + 90 * idx + 39);

            auto t_xyyy_yz = primBuffer.data(koff + 90 * idx + 40);

            auto t_xyyy_zz = primBuffer.data(koff + 90 * idx + 41);

            auto t_xyyz_xx = primBuffer.data(koff + 90 * idx + 42);

            auto t_xyyz_xy = primBuffer.data(koff + 90 * idx + 43);

            auto t_xyyz_xz = primBuffer.data(koff + 90 * idx + 44);

            auto t_xyyz_yy = primBuffer.data(koff + 90 * idx + 45);

            auto t_xyyz_yz = primBuffer.data(koff + 90 * idx + 46);

            auto t_xyyz_zz = primBuffer.data(koff + 90 * idx + 47);

            auto t_xyzz_xx = primBuffer.data(koff + 90 * idx + 48);

            auto t_xyzz_xy = primBuffer.data(koff + 90 * idx + 49);

            auto t_xyzz_xz = primBuffer.data(koff + 90 * idx + 50);

            auto t_xyzz_yy = primBuffer.data(koff + 90 * idx + 51);

            auto t_xyzz_yz = primBuffer.data(koff + 90 * idx + 52);

            auto t_xyzz_zz = primBuffer.data(koff + 90 * idx + 53);

            auto t_xzzz_xx = primBuffer.data(koff + 90 * idx + 54);

            auto t_xzzz_xy = primBuffer.data(koff + 90 * idx + 55);

            auto t_xzzz_xz = primBuffer.data(koff + 90 * idx + 56);

            auto t_xzzz_yy = primBuffer.data(koff + 90 * idx + 57);

            auto t_xzzz_yz = primBuffer.data(koff + 90 * idx + 58);

            auto t_xzzz_zz = primBuffer.data(koff + 90 * idx + 59);

            auto t_yyyy_xx = primBuffer.data(koff + 90 * idx + 60);

            auto t_yyyy_xy = primBuffer.data(koff + 90 * idx + 61);

            auto t_yyyy_xz = primBuffer.data(koff + 90 * idx + 62);

            auto t_yyyy_yy = primBuffer.data(koff + 90 * idx + 63);

            auto t_yyyy_yz = primBuffer.data(koff + 90 * idx + 64);

            auto t_yyyy_zz = primBuffer.data(koff + 90 * idx + 65);

            auto t_yyyz_xx = primBuffer.data(koff + 90 * idx + 66);

            auto t_yyyz_xy = primBuffer.data(koff + 90 * idx + 67);

            auto t_yyyz_xz = primBuffer.data(koff + 90 * idx + 68);

            auto t_yyyz_yy = primBuffer.data(koff + 90 * idx + 69);

            auto t_yyyz_yz = primBuffer.data(koff + 90 * idx + 70);

            auto t_yyyz_zz = primBuffer.data(koff + 90 * idx + 71);

            auto t_yyzz_xx = primBuffer.data(koff + 90 * idx + 72);

            auto t_yyzz_xy = primBuffer.data(koff + 90 * idx + 73);

            auto t_yyzz_xz = primBuffer.data(koff + 90 * idx + 74);

            auto t_yyzz_yy = primBuffer.data(koff + 90 * idx + 75);

            auto t_yyzz_yz = primBuffer.data(koff + 90 * idx + 76);

            auto t_yyzz_zz = primBuffer.data(koff + 90 * idx + 77);

            auto t_yzzz_xx = primBuffer.data(koff + 90 * idx + 78);

            auto t_yzzz_xy = primBuffer.data(koff + 90 * idx + 79);

            auto t_yzzz_xz = primBuffer.data(koff + 90 * idx + 80);

            auto t_yzzz_yy = primBuffer.data(koff + 90 * idx + 81);

            auto t_yzzz_yz = primBuffer.data(koff + 90 * idx + 82);

            auto t_yzzz_zz = primBuffer.data(koff + 90 * idx + 83);

            auto t_zzzz_xx = primBuffer.data(koff + 90 * idx + 84);

            auto t_zzzz_xy = primBuffer.data(koff + 90 * idx + 85);

            auto t_zzzz_xz = primBuffer.data(koff + 90 * idx + 86);

            auto t_zzzz_yy = primBuffer.data(koff + 90 * idx + 87);

            auto t_zzzz_yz = primBuffer.data(koff + 90 * idx + 88);

            auto t_zzzz_zz = primBuffer.data(koff + 90 * idx + 89);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xxx_x, s_xxx_y,\
                                     s_xxx_z, s_xxy_x, s_xxy_y, s_xxy_z, s_xxz_x,\
                                     s_xxz_y, s_xxz_z, s_xyy_x, s_xyy_y, s_xyy_z,\
                                     s_xyz_x, s_xyz_y, s_xyz_z, s_xzz_x, s_xzz_y,\
                                     s_xzz_z, s_yyy_x, s_yyy_y, s_yyy_z, s_yyz_x,\
                                     s_yyz_y, s_yyz_z, s_yzz_x, s_yzz_y, s_yzz_z,\
                                     s_zzz_x, s_zzz_y, s_zzz_z, t_xxx_x, t_xxx_y,\
                                     t_xxx_z, t_xxy_x, t_xxy_y, t_xxy_z, t_xxz_x,\
                                     t_xxz_y, t_xxz_z, t_xyy_x, t_xyy_y, t_xyy_z,\
                                     t_xyz_x, t_xyz_y, t_xyz_z, t_xzz_x, t_xzz_y,\
                                     t_xzz_z, t_yyy_x, t_yyy_y, t_yyy_z, t_yyz_x,\
                                     t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z,\
                                     t_zzz_x, t_zzz_y, t_zzz_z, s_xx_xx, s_xx_xy,\
                                     s_xx_xz, s_xx_yy, s_xx_yz, s_xx_zz, s_xy_xx,\
                                     s_xy_xy, s_xy_xz, s_xy_yy, s_xy_yz, s_xy_zz,\
                                     s_xz_xx, s_xz_xy, s_xz_xz, s_xz_yy, s_xz_yz,\
                                     s_xz_zz, s_yy_xx, s_yy_xy, s_yy_xz, s_yy_yy,\
                                     s_yy_yz, s_yy_zz, s_yz_xx, s_yz_xy, s_yz_xz,\
                                     s_yz_yy, s_yz_yz, s_yz_zz, s_zz_xx, s_zz_xy,\
                                     s_zz_xz, s_zz_yy, s_zz_yz, s_zz_zz, t_xx_xx,\
                                     t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz,\
                                     t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz,\
                                     t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy,\
                                     t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy, t_yy_xz,\
                                     t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, t_yz_xy,\
                                     t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx,\
                                     t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz,\
                                     s_xxx_xx, s_xxx_xy, s_xxx_xz, s_xxx_yy, s_xxx_yz,\
                                     s_xxx_zz, s_xxy_xx, s_xxy_xy, s_xxy_xz, s_xxy_yy,\
                                     s_xxy_yz, s_xxy_zz, s_xxz_xx, s_xxz_xy, s_xxz_xz,\
                                     s_xxz_yy, s_xxz_yz, s_xxz_zz, s_xyy_xx, s_xyy_xy,\
                                     s_xyy_xz, s_xyy_yy, s_xyy_yz, s_xyy_zz, s_xyz_xx,\
                                     s_xyz_xy, s_xyz_xz, s_xyz_yy, s_xyz_yz, s_xyz_zz,\
                                     s_xzz_xx, s_xzz_xy, s_xzz_xz, s_xzz_yy, s_xzz_yz,\
                                     s_xzz_zz, s_yyy_xx, s_yyy_xy, s_yyy_xz, s_yyy_yy,\
                                     s_yyy_yz, s_yyy_zz, s_yyz_xx, s_yyz_xy, s_yyz_xz,\
                                     s_yyz_yy, s_yyz_yz, s_yyz_zz, s_yzz_xx, s_yzz_xy,\
                                     s_yzz_xz, s_yzz_yy, s_yzz_yz, s_yzz_zz, s_zzz_xx,\
                                     s_zzz_xy, s_zzz_xz, s_zzz_yy, s_zzz_yz, s_zzz_zz,\
                                     t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_yy, t_xxx_yz,\
                                     t_xxx_zz, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy,\
                                     t_xxy_yz, t_xxy_zz, t_xxz_xx, t_xxz_xy, t_xxz_xz,\
                                     t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, t_xyy_xy,\
                                     t_xyy_xz, t_xyy_yy, t_xyy_yz, t_xyy_zz, t_xyz_xx,\
                                     t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz,\
                                     t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy, t_xzz_yz,\
                                     t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy,\
                                     t_yyy_yz, t_yyy_zz, t_yyz_xx, t_yyz_xy, t_yyz_xz,\
                                     t_yyz_yy, t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy,\
                                     t_yzz_xz, t_yzz_yy, t_yzz_yz, t_yzz_zz, t_zzz_xx,\
                                     t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, t_zzz_zz,\
                                     s_xxxx_xx, s_xxxx_xy, s_xxxx_xz, s_xxxx_yy,\
                                     s_xxxx_yz, s_xxxx_zz, s_xxxy_xx, s_xxxy_xy,\
                                     s_xxxy_xz, s_xxxy_yy, s_xxxy_yz, s_xxxy_zz,\
                                     s_xxxz_xx, s_xxxz_xy, s_xxxz_xz, s_xxxz_yy,\
                                     s_xxxz_yz, s_xxxz_zz, s_xxyy_xx, s_xxyy_xy,\
                                     s_xxyy_xz, s_xxyy_yy, s_xxyy_yz, s_xxyy_zz,\
                                     s_xxyz_xx, s_xxyz_xy, s_xxyz_xz, s_xxyz_yy,\
                                     s_xxyz_yz, s_xxyz_zz, s_xxzz_xx, s_xxzz_xy,\
                                     s_xxzz_xz, s_xxzz_yy, s_xxzz_yz, s_xxzz_zz,\
                                     s_xyyy_xx, s_xyyy_xy, s_xyyy_xz, s_xyyy_yy,\
                                     s_xyyy_yz, s_xyyy_zz, s_xyyz_xx, s_xyyz_xy,\
                                     s_xyyz_xz, s_xyyz_yy, s_xyyz_yz, s_xyyz_zz,\
                                     s_xyzz_xx, s_xyzz_xy, s_xyzz_xz, s_xyzz_yy,\
                                     s_xyzz_yz, s_xyzz_zz, s_xzzz_xx, s_xzzz_xy,\
                                     s_xzzz_xz, s_xzzz_yy, s_xzzz_yz, s_xzzz_zz,\
                                     s_yyyy_xx, s_yyyy_xy, s_yyyy_xz, s_yyyy_yy,\
                                     s_yyyy_yz, s_yyyy_zz, s_yyyz_xx, s_yyyz_xy,\
                                     s_yyyz_xz, s_yyyz_yy, s_yyyz_yz, s_yyyz_zz,\
                                     s_yyzz_xx, s_yyzz_xy, s_yyzz_xz, s_yyzz_yy,\
                                     s_yyzz_yz, s_yyzz_zz, s_yzzz_xx, s_yzzz_xy,\
                                     s_yzzz_xz, s_yzzz_yy, s_yzzz_yz, s_yzzz_zz,\
                                     s_zzzz_xx, s_zzzz_xy, s_zzzz_xz, s_zzzz_yy,\
                                     s_zzzz_yz, s_zzzz_zz, t_xxxx_xx, t_xxxx_xy,\
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
                                     t_zzzz_xz, t_zzzz_yy, t_zzzz_yz, t_zzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxxx_xx[j] = fr * s_xxx_xx[j] + f2t * (3.0 * s_xx_xx[j] + 2.0 * s_xxx_x[j]);

                t_xxxx_xx[j] = fr * t_xxx_xx[j] + f2t * (3.0 * t_xx_xx[j] + 2.0 * t_xxx_x[j]) + f2z * (s_xxxx_xx[j] - f2a * 3.0 * s_xx_xx[j]);

                s_xxxx_xy[j] = fr * s_xxx_xy[j] + f2t * (3.0 * s_xx_xy[j] + s_xxx_y[j]);

                t_xxxx_xy[j] = fr * t_xxx_xy[j] + f2t * (3.0 * t_xx_xy[j] + t_xxx_y[j]) + f2z * (s_xxxx_xy[j] - f2a * 3.0 * s_xx_xy[j]);

                s_xxxx_xz[j] = fr * s_xxx_xz[j] + f2t * (3.0 * s_xx_xz[j] + s_xxx_z[j]);

                t_xxxx_xz[j] = fr * t_xxx_xz[j] + f2t * (3.0 * t_xx_xz[j] + t_xxx_z[j]) + f2z * (s_xxxx_xz[j] - f2a * 3.0 * s_xx_xz[j]);

                s_xxxx_yy[j] = fr * s_xxx_yy[j] + f2t * 3.0 * s_xx_yy[j];

                t_xxxx_yy[j] = fr * t_xxx_yy[j] + f2t * 3.0 * t_xx_yy[j] + f2z * (s_xxxx_yy[j] - f2a * 3.0 * s_xx_yy[j]);

                s_xxxx_yz[j] = fr * s_xxx_yz[j] + f2t * 3.0 * s_xx_yz[j];

                t_xxxx_yz[j] = fr * t_xxx_yz[j] + f2t * 3.0 * t_xx_yz[j] + f2z * (s_xxxx_yz[j] - f2a * 3.0 * s_xx_yz[j]);

                s_xxxx_zz[j] = fr * s_xxx_zz[j] + f2t * 3.0 * s_xx_zz[j];

                t_xxxx_zz[j] = fr * t_xxx_zz[j] + f2t * 3.0 * t_xx_zz[j] + f2z * (s_xxxx_zz[j] - f2a * 3.0 * s_xx_zz[j]);

                s_xxxy_xx[j] = fr * s_xxy_xx[j] + f2t * (2.0 * s_xy_xx[j] + 2.0 * s_xxy_x[j]);

                t_xxxy_xx[j] = fr * t_xxy_xx[j] + f2t * (2.0 * t_xy_xx[j] + 2.0 * t_xxy_x[j]) + f2z * (s_xxxy_xx[j] - f2a * 2.0 * s_xy_xx[j]);

                s_xxxy_xy[j] = fr * s_xxy_xy[j] + f2t * (2.0 * s_xy_xy[j] + s_xxy_y[j]);

                t_xxxy_xy[j] = fr * t_xxy_xy[j] + f2t * (2.0 * t_xy_xy[j] + t_xxy_y[j]) + f2z * (s_xxxy_xy[j] - f2a * 2.0 * s_xy_xy[j]);

                s_xxxy_xz[j] = fr * s_xxy_xz[j] + f2t * (2.0 * s_xy_xz[j] + s_xxy_z[j]);

                t_xxxy_xz[j] = fr * t_xxy_xz[j] + f2t * (2.0 * t_xy_xz[j] + t_xxy_z[j]) + f2z * (s_xxxy_xz[j] - f2a * 2.0 * s_xy_xz[j]);

                s_xxxy_yy[j] = fr * s_xxy_yy[j] + f2t * 2.0 * s_xy_yy[j];

                t_xxxy_yy[j] = fr * t_xxy_yy[j] + f2t * 2.0 * t_xy_yy[j] + f2z * (s_xxxy_yy[j] - f2a * 2.0 * s_xy_yy[j]);

                s_xxxy_yz[j] = fr * s_xxy_yz[j] + f2t * 2.0 * s_xy_yz[j];

                t_xxxy_yz[j] = fr * t_xxy_yz[j] + f2t * 2.0 * t_xy_yz[j] + f2z * (s_xxxy_yz[j] - f2a * 2.0 * s_xy_yz[j]);

                s_xxxy_zz[j] = fr * s_xxy_zz[j] + f2t * 2.0 * s_xy_zz[j];

                t_xxxy_zz[j] = fr * t_xxy_zz[j] + f2t * 2.0 * t_xy_zz[j] + f2z * (s_xxxy_zz[j] - f2a * 2.0 * s_xy_zz[j]);

                s_xxxz_xx[j] = fr * s_xxz_xx[j] + f2t * (2.0 * s_xz_xx[j] + 2.0 * s_xxz_x[j]);

                t_xxxz_xx[j] = fr * t_xxz_xx[j] + f2t * (2.0 * t_xz_xx[j] + 2.0 * t_xxz_x[j]) + f2z * (s_xxxz_xx[j] - f2a * 2.0 * s_xz_xx[j]);

                s_xxxz_xy[j] = fr * s_xxz_xy[j] + f2t * (2.0 * s_xz_xy[j] + s_xxz_y[j]);

                t_xxxz_xy[j] = fr * t_xxz_xy[j] + f2t * (2.0 * t_xz_xy[j] + t_xxz_y[j]) + f2z * (s_xxxz_xy[j] - f2a * 2.0 * s_xz_xy[j]);

                s_xxxz_xz[j] = fr * s_xxz_xz[j] + f2t * (2.0 * s_xz_xz[j] + s_xxz_z[j]);

                t_xxxz_xz[j] = fr * t_xxz_xz[j] + f2t * (2.0 * t_xz_xz[j] + t_xxz_z[j]) + f2z * (s_xxxz_xz[j] - f2a * 2.0 * s_xz_xz[j]);

                s_xxxz_yy[j] = fr * s_xxz_yy[j] + f2t * 2.0 * s_xz_yy[j];

                t_xxxz_yy[j] = fr * t_xxz_yy[j] + f2t * 2.0 * t_xz_yy[j] + f2z * (s_xxxz_yy[j] - f2a * 2.0 * s_xz_yy[j]);

                s_xxxz_yz[j] = fr * s_xxz_yz[j] + f2t * 2.0 * s_xz_yz[j];

                t_xxxz_yz[j] = fr * t_xxz_yz[j] + f2t * 2.0 * t_xz_yz[j] + f2z * (s_xxxz_yz[j] - f2a * 2.0 * s_xz_yz[j]);

                s_xxxz_zz[j] = fr * s_xxz_zz[j] + f2t * 2.0 * s_xz_zz[j];

                t_xxxz_zz[j] = fr * t_xxz_zz[j] + f2t * 2.0 * t_xz_zz[j] + f2z * (s_xxxz_zz[j] - f2a * 2.0 * s_xz_zz[j]);

                s_xxyy_xx[j] = fr * s_xyy_xx[j] + f2t * (s_yy_xx[j] + 2.0 * s_xyy_x[j]);

                t_xxyy_xx[j] = fr * t_xyy_xx[j] + f2t * (t_yy_xx[j] + 2.0 * t_xyy_x[j]) + f2z * (s_xxyy_xx[j] - f2a * s_yy_xx[j]);

                s_xxyy_xy[j] = fr * s_xyy_xy[j] + f2t * (s_yy_xy[j] + s_xyy_y[j]);

                t_xxyy_xy[j] = fr * t_xyy_xy[j] + f2t * (t_yy_xy[j] + t_xyy_y[j]) + f2z * (s_xxyy_xy[j] - f2a * s_yy_xy[j]);

                s_xxyy_xz[j] = fr * s_xyy_xz[j] + f2t * (s_yy_xz[j] + s_xyy_z[j]);

                t_xxyy_xz[j] = fr * t_xyy_xz[j] + f2t * (t_yy_xz[j] + t_xyy_z[j]) + f2z * (s_xxyy_xz[j] - f2a * s_yy_xz[j]);

                s_xxyy_yy[j] = fr * s_xyy_yy[j] + f2t * s_yy_yy[j];

                t_xxyy_yy[j] = fr * t_xyy_yy[j] + f2t * t_yy_yy[j] + f2z * (s_xxyy_yy[j] - f2a * s_yy_yy[j]);

                s_xxyy_yz[j] = fr * s_xyy_yz[j] + f2t * s_yy_yz[j];

                t_xxyy_yz[j] = fr * t_xyy_yz[j] + f2t * t_yy_yz[j] + f2z * (s_xxyy_yz[j] - f2a * s_yy_yz[j]);

                s_xxyy_zz[j] = fr * s_xyy_zz[j] + f2t * s_yy_zz[j];

                t_xxyy_zz[j] = fr * t_xyy_zz[j] + f2t * t_yy_zz[j] + f2z * (s_xxyy_zz[j] - f2a * s_yy_zz[j]);

                s_xxyz_xx[j] = fr * s_xyz_xx[j] + f2t * (s_yz_xx[j] + 2.0 * s_xyz_x[j]);

                t_xxyz_xx[j] = fr * t_xyz_xx[j] + f2t * (t_yz_xx[j] + 2.0 * t_xyz_x[j]) + f2z * (s_xxyz_xx[j] - f2a * s_yz_xx[j]);

                s_xxyz_xy[j] = fr * s_xyz_xy[j] + f2t * (s_yz_xy[j] + s_xyz_y[j]);

                t_xxyz_xy[j] = fr * t_xyz_xy[j] + f2t * (t_yz_xy[j] + t_xyz_y[j]) + f2z * (s_xxyz_xy[j] - f2a * s_yz_xy[j]);

                s_xxyz_xz[j] = fr * s_xyz_xz[j] + f2t * (s_yz_xz[j] + s_xyz_z[j]);

                t_xxyz_xz[j] = fr * t_xyz_xz[j] + f2t * (t_yz_xz[j] + t_xyz_z[j]) + f2z * (s_xxyz_xz[j] - f2a * s_yz_xz[j]);

                s_xxyz_yy[j] = fr * s_xyz_yy[j] + f2t * s_yz_yy[j];

                t_xxyz_yy[j] = fr * t_xyz_yy[j] + f2t * t_yz_yy[j] + f2z * (s_xxyz_yy[j] - f2a * s_yz_yy[j]);

                s_xxyz_yz[j] = fr * s_xyz_yz[j] + f2t * s_yz_yz[j];

                t_xxyz_yz[j] = fr * t_xyz_yz[j] + f2t * t_yz_yz[j] + f2z * (s_xxyz_yz[j] - f2a * s_yz_yz[j]);

                s_xxyz_zz[j] = fr * s_xyz_zz[j] + f2t * s_yz_zz[j];

                t_xxyz_zz[j] = fr * t_xyz_zz[j] + f2t * t_yz_zz[j] + f2z * (s_xxyz_zz[j] - f2a * s_yz_zz[j]);

                s_xxzz_xx[j] = fr * s_xzz_xx[j] + f2t * (s_zz_xx[j] + 2.0 * s_xzz_x[j]);

                t_xxzz_xx[j] = fr * t_xzz_xx[j] + f2t * (t_zz_xx[j] + 2.0 * t_xzz_x[j]) + f2z * (s_xxzz_xx[j] - f2a * s_zz_xx[j]);

                s_xxzz_xy[j] = fr * s_xzz_xy[j] + f2t * (s_zz_xy[j] + s_xzz_y[j]);

                t_xxzz_xy[j] = fr * t_xzz_xy[j] + f2t * (t_zz_xy[j] + t_xzz_y[j]) + f2z * (s_xxzz_xy[j] - f2a * s_zz_xy[j]);

                s_xxzz_xz[j] = fr * s_xzz_xz[j] + f2t * (s_zz_xz[j] + s_xzz_z[j]);

                t_xxzz_xz[j] = fr * t_xzz_xz[j] + f2t * (t_zz_xz[j] + t_xzz_z[j]) + f2z * (s_xxzz_xz[j] - f2a * s_zz_xz[j]);

                s_xxzz_yy[j] = fr * s_xzz_yy[j] + f2t * s_zz_yy[j];

                t_xxzz_yy[j] = fr * t_xzz_yy[j] + f2t * t_zz_yy[j] + f2z * (s_xxzz_yy[j] - f2a * s_zz_yy[j]);

                s_xxzz_yz[j] = fr * s_xzz_yz[j] + f2t * s_zz_yz[j];

                t_xxzz_yz[j] = fr * t_xzz_yz[j] + f2t * t_zz_yz[j] + f2z * (s_xxzz_yz[j] - f2a * s_zz_yz[j]);

                s_xxzz_zz[j] = fr * s_xzz_zz[j] + f2t * s_zz_zz[j];

                t_xxzz_zz[j] = fr * t_xzz_zz[j] + f2t * t_zz_zz[j] + f2z * (s_xxzz_zz[j] - f2a * s_zz_zz[j]);

                s_xyyy_xx[j] = fr * s_yyy_xx[j] + f2t * 2.0 * s_yyy_x[j];

                t_xyyy_xx[j] = fr * t_yyy_xx[j] + f2t * 2.0 * t_yyy_x[j] + f2z * s_xyyy_xx[j];

                s_xyyy_xy[j] = fr * s_yyy_xy[j] + f2t * s_yyy_y[j];

                t_xyyy_xy[j] = fr * t_yyy_xy[j] + f2t * t_yyy_y[j] + f2z * s_xyyy_xy[j];

                s_xyyy_xz[j] = fr * s_yyy_xz[j] + f2t * s_yyy_z[j];

                t_xyyy_xz[j] = fr * t_yyy_xz[j] + f2t * t_yyy_z[j] + f2z * s_xyyy_xz[j];

                s_xyyy_yy[j] = fr * s_yyy_yy[j];

                t_xyyy_yy[j] = fr * t_yyy_yy[j] + f2z * s_xyyy_yy[j];

                s_xyyy_yz[j] = fr * s_yyy_yz[j];

                t_xyyy_yz[j] = fr * t_yyy_yz[j] + f2z * s_xyyy_yz[j];

                s_xyyy_zz[j] = fr * s_yyy_zz[j];

                t_xyyy_zz[j] = fr * t_yyy_zz[j] + f2z * s_xyyy_zz[j];

                s_xyyz_xx[j] = fr * s_yyz_xx[j] + f2t * 2.0 * s_yyz_x[j];

                t_xyyz_xx[j] = fr * t_yyz_xx[j] + f2t * 2.0 * t_yyz_x[j] + f2z * s_xyyz_xx[j];

                s_xyyz_xy[j] = fr * s_yyz_xy[j] + f2t * s_yyz_y[j];

                t_xyyz_xy[j] = fr * t_yyz_xy[j] + f2t * t_yyz_y[j] + f2z * s_xyyz_xy[j];

                s_xyyz_xz[j] = fr * s_yyz_xz[j] + f2t * s_yyz_z[j];

                t_xyyz_xz[j] = fr * t_yyz_xz[j] + f2t * t_yyz_z[j] + f2z * s_xyyz_xz[j];

                s_xyyz_yy[j] = fr * s_yyz_yy[j];

                t_xyyz_yy[j] = fr * t_yyz_yy[j] + f2z * s_xyyz_yy[j];

                s_xyyz_yz[j] = fr * s_yyz_yz[j];

                t_xyyz_yz[j] = fr * t_yyz_yz[j] + f2z * s_xyyz_yz[j];

                s_xyyz_zz[j] = fr * s_yyz_zz[j];

                t_xyyz_zz[j] = fr * t_yyz_zz[j] + f2z * s_xyyz_zz[j];

                s_xyzz_xx[j] = fr * s_yzz_xx[j] + f2t * 2.0 * s_yzz_x[j];

                t_xyzz_xx[j] = fr * t_yzz_xx[j] + f2t * 2.0 * t_yzz_x[j] + f2z * s_xyzz_xx[j];

                s_xyzz_xy[j] = fr * s_yzz_xy[j] + f2t * s_yzz_y[j];

                t_xyzz_xy[j] = fr * t_yzz_xy[j] + f2t * t_yzz_y[j] + f2z * s_xyzz_xy[j];

                s_xyzz_xz[j] = fr * s_yzz_xz[j] + f2t * s_yzz_z[j];

                t_xyzz_xz[j] = fr * t_yzz_xz[j] + f2t * t_yzz_z[j] + f2z * s_xyzz_xz[j];

                s_xyzz_yy[j] = fr * s_yzz_yy[j];

                t_xyzz_yy[j] = fr * t_yzz_yy[j] + f2z * s_xyzz_yy[j];

                s_xyzz_yz[j] = fr * s_yzz_yz[j];

                t_xyzz_yz[j] = fr * t_yzz_yz[j] + f2z * s_xyzz_yz[j];

                s_xyzz_zz[j] = fr * s_yzz_zz[j];

                t_xyzz_zz[j] = fr * t_yzz_zz[j] + f2z * s_xyzz_zz[j];

                s_xzzz_xx[j] = fr * s_zzz_xx[j] + f2t * 2.0 * s_zzz_x[j];

                t_xzzz_xx[j] = fr * t_zzz_xx[j] + f2t * 2.0 * t_zzz_x[j] + f2z * s_xzzz_xx[j];

                s_xzzz_xy[j] = fr * s_zzz_xy[j] + f2t * s_zzz_y[j];

                t_xzzz_xy[j] = fr * t_zzz_xy[j] + f2t * t_zzz_y[j] + f2z * s_xzzz_xy[j];

                s_xzzz_xz[j] = fr * s_zzz_xz[j] + f2t * s_zzz_z[j];

                t_xzzz_xz[j] = fr * t_zzz_xz[j] + f2t * t_zzz_z[j] + f2z * s_xzzz_xz[j];

                s_xzzz_yy[j] = fr * s_zzz_yy[j];

                t_xzzz_yy[j] = fr * t_zzz_yy[j] + f2z * s_xzzz_yy[j];

                s_xzzz_yz[j] = fr * s_zzz_yz[j];

                t_xzzz_yz[j] = fr * t_zzz_yz[j] + f2z * s_xzzz_yz[j];

                s_xzzz_zz[j] = fr * s_zzz_zz[j];

                t_xzzz_zz[j] = fr * t_zzz_zz[j] + f2z * s_xzzz_zz[j];

                // leading y component

                fr = pay[j];

                s_yyyy_xx[j] = fr * s_yyy_xx[j] + f2t * 3.0 * s_yy_xx[j];

                t_yyyy_xx[j] = fr * t_yyy_xx[j] + f2t * 3.0 * t_yy_xx[j] + f2z * (s_yyyy_xx[j] - f2a * 3.0 * s_yy_xx[j]);

                s_yyyy_xy[j] = fr * s_yyy_xy[j] + f2t * (3.0 * s_yy_xy[j] + s_yyy_x[j]);

                t_yyyy_xy[j] = fr * t_yyy_xy[j] + f2t * (3.0 * t_yy_xy[j] + t_yyy_x[j]) + f2z * (s_yyyy_xy[j] - f2a * 3.0 * s_yy_xy[j]);

                s_yyyy_xz[j] = fr * s_yyy_xz[j] + f2t * 3.0 * s_yy_xz[j];

                t_yyyy_xz[j] = fr * t_yyy_xz[j] + f2t * 3.0 * t_yy_xz[j] + f2z * (s_yyyy_xz[j] - f2a * 3.0 * s_yy_xz[j]);

                s_yyyy_yy[j] = fr * s_yyy_yy[j] + f2t * (3.0 * s_yy_yy[j] + 2.0 * s_yyy_y[j]);

                t_yyyy_yy[j] = fr * t_yyy_yy[j] + f2t * (3.0 * t_yy_yy[j] + 2.0 * t_yyy_y[j]) + f2z * (s_yyyy_yy[j] - f2a * 3.0 * s_yy_yy[j]);

                s_yyyy_yz[j] = fr * s_yyy_yz[j] + f2t * (3.0 * s_yy_yz[j] + s_yyy_z[j]);

                t_yyyy_yz[j] = fr * t_yyy_yz[j] + f2t * (3.0 * t_yy_yz[j] + t_yyy_z[j]) + f2z * (s_yyyy_yz[j] - f2a * 3.0 * s_yy_yz[j]);

                s_yyyy_zz[j] = fr * s_yyy_zz[j] + f2t * 3.0 * s_yy_zz[j];

                t_yyyy_zz[j] = fr * t_yyy_zz[j] + f2t * 3.0 * t_yy_zz[j] + f2z * (s_yyyy_zz[j] - f2a * 3.0 * s_yy_zz[j]);

                s_yyyz_xx[j] = fr * s_yyz_xx[j] + f2t * 2.0 * s_yz_xx[j];

                t_yyyz_xx[j] = fr * t_yyz_xx[j] + f2t * 2.0 * t_yz_xx[j] + f2z * (s_yyyz_xx[j] - f2a * 2.0 * s_yz_xx[j]);

                s_yyyz_xy[j] = fr * s_yyz_xy[j] + f2t * (2.0 * s_yz_xy[j] + s_yyz_x[j]);

                t_yyyz_xy[j] = fr * t_yyz_xy[j] + f2t * (2.0 * t_yz_xy[j] + t_yyz_x[j]) + f2z * (s_yyyz_xy[j] - f2a * 2.0 * s_yz_xy[j]);

                s_yyyz_xz[j] = fr * s_yyz_xz[j] + f2t * 2.0 * s_yz_xz[j];

                t_yyyz_xz[j] = fr * t_yyz_xz[j] + f2t * 2.0 * t_yz_xz[j] + f2z * (s_yyyz_xz[j] - f2a * 2.0 * s_yz_xz[j]);

                s_yyyz_yy[j] = fr * s_yyz_yy[j] + f2t * (2.0 * s_yz_yy[j] + 2.0 * s_yyz_y[j]);

                t_yyyz_yy[j] = fr * t_yyz_yy[j] + f2t * (2.0 * t_yz_yy[j] + 2.0 * t_yyz_y[j]) + f2z * (s_yyyz_yy[j] - f2a * 2.0 * s_yz_yy[j]);

                s_yyyz_yz[j] = fr * s_yyz_yz[j] + f2t * (2.0 * s_yz_yz[j] + s_yyz_z[j]);

                t_yyyz_yz[j] = fr * t_yyz_yz[j] + f2t * (2.0 * t_yz_yz[j] + t_yyz_z[j]) + f2z * (s_yyyz_yz[j] - f2a * 2.0 * s_yz_yz[j]);

                s_yyyz_zz[j] = fr * s_yyz_zz[j] + f2t * 2.0 * s_yz_zz[j];

                t_yyyz_zz[j] = fr * t_yyz_zz[j] + f2t * 2.0 * t_yz_zz[j] + f2z * (s_yyyz_zz[j] - f2a * 2.0 * s_yz_zz[j]);

                s_yyzz_xx[j] = fr * s_yzz_xx[j] + f2t * s_zz_xx[j];

                t_yyzz_xx[j] = fr * t_yzz_xx[j] + f2t * t_zz_xx[j] + f2z * (s_yyzz_xx[j] - f2a * s_zz_xx[j]);

                s_yyzz_xy[j] = fr * s_yzz_xy[j] + f2t * (s_zz_xy[j] + s_yzz_x[j]);

                t_yyzz_xy[j] = fr * t_yzz_xy[j] + f2t * (t_zz_xy[j] + t_yzz_x[j]) + f2z * (s_yyzz_xy[j] - f2a * s_zz_xy[j]);

                s_yyzz_xz[j] = fr * s_yzz_xz[j] + f2t * s_zz_xz[j];

                t_yyzz_xz[j] = fr * t_yzz_xz[j] + f2t * t_zz_xz[j] + f2z * (s_yyzz_xz[j] - f2a * s_zz_xz[j]);

                s_yyzz_yy[j] = fr * s_yzz_yy[j] + f2t * (s_zz_yy[j] + 2.0 * s_yzz_y[j]);

                t_yyzz_yy[j] = fr * t_yzz_yy[j] + f2t * (t_zz_yy[j] + 2.0 * t_yzz_y[j]) + f2z * (s_yyzz_yy[j] - f2a * s_zz_yy[j]);

                s_yyzz_yz[j] = fr * s_yzz_yz[j] + f2t * (s_zz_yz[j] + s_yzz_z[j]);

                t_yyzz_yz[j] = fr * t_yzz_yz[j] + f2t * (t_zz_yz[j] + t_yzz_z[j]) + f2z * (s_yyzz_yz[j] - f2a * s_zz_yz[j]);

                s_yyzz_zz[j] = fr * s_yzz_zz[j] + f2t * s_zz_zz[j];

                t_yyzz_zz[j] = fr * t_yzz_zz[j] + f2t * t_zz_zz[j] + f2z * (s_yyzz_zz[j] - f2a * s_zz_zz[j]);

                s_yzzz_xx[j] = fr * s_zzz_xx[j];

                t_yzzz_xx[j] = fr * t_zzz_xx[j] + f2z * s_yzzz_xx[j];

                s_yzzz_xy[j] = fr * s_zzz_xy[j] + f2t * s_zzz_x[j];

                t_yzzz_xy[j] = fr * t_zzz_xy[j] + f2t * t_zzz_x[j] + f2z * s_yzzz_xy[j];

                s_yzzz_xz[j] = fr * s_zzz_xz[j];

                t_yzzz_xz[j] = fr * t_zzz_xz[j] + f2z * s_yzzz_xz[j];

                s_yzzz_yy[j] = fr * s_zzz_yy[j] + f2t * 2.0 * s_zzz_y[j];

                t_yzzz_yy[j] = fr * t_zzz_yy[j] + f2t * 2.0 * t_zzz_y[j] + f2z * s_yzzz_yy[j];

                s_yzzz_yz[j] = fr * s_zzz_yz[j] + f2t * s_zzz_z[j];

                t_yzzz_yz[j] = fr * t_zzz_yz[j] + f2t * t_zzz_z[j] + f2z * s_yzzz_yz[j];

                s_yzzz_zz[j] = fr * s_zzz_zz[j];

                t_yzzz_zz[j] = fr * t_zzz_zz[j] + f2z * s_yzzz_zz[j];

                // leading z component

                fr = paz[j];

                s_zzzz_xx[j] = fr * s_zzz_xx[j] + f2t * 3.0 * s_zz_xx[j];

                t_zzzz_xx[j] = fr * t_zzz_xx[j] + f2t * 3.0 * t_zz_xx[j] + f2z * (s_zzzz_xx[j] - f2a * 3.0 * s_zz_xx[j]);

                s_zzzz_xy[j] = fr * s_zzz_xy[j] + f2t * 3.0 * s_zz_xy[j];

                t_zzzz_xy[j] = fr * t_zzz_xy[j] + f2t * 3.0 * t_zz_xy[j] + f2z * (s_zzzz_xy[j] - f2a * 3.0 * s_zz_xy[j]);

                s_zzzz_xz[j] = fr * s_zzz_xz[j] + f2t * (3.0 * s_zz_xz[j] + s_zzz_x[j]);

                t_zzzz_xz[j] = fr * t_zzz_xz[j] + f2t * (3.0 * t_zz_xz[j] + t_zzz_x[j]) + f2z * (s_zzzz_xz[j] - f2a * 3.0 * s_zz_xz[j]);

                s_zzzz_yy[j] = fr * s_zzz_yy[j] + f2t * 3.0 * s_zz_yy[j];

                t_zzzz_yy[j] = fr * t_zzz_yy[j] + f2t * 3.0 * t_zz_yy[j] + f2z * (s_zzzz_yy[j] - f2a * 3.0 * s_zz_yy[j]);

                s_zzzz_yz[j] = fr * s_zzz_yz[j] + f2t * (3.0 * s_zz_yz[j] + s_zzz_y[j]);

                t_zzzz_yz[j] = fr * t_zzz_yz[j] + f2t * (3.0 * t_zz_yz[j] + t_zzz_y[j]) + f2z * (s_zzzz_yz[j] - f2a * 3.0 * s_zz_yz[j]);

                s_zzzz_zz[j] = fr * s_zzz_zz[j] + f2t * (3.0 * s_zz_zz[j] + 2.0 * s_zzz_z[j]);

                t_zzzz_zz[j] = fr * t_zzz_zz[j] + f2t * (3.0 * t_zz_zz[j] + 2.0 * t_zzz_z[j]) + f2z * (s_zzzz_zz[j] - f2a * 3.0 * s_zz_zz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFG(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|F) integrals

            auto s_xx_xxx = primBuffer.data(skoff + 60 * idx);

            auto s_xx_xxy = primBuffer.data(skoff + 60 * idx + 1);

            auto s_xx_xxz = primBuffer.data(skoff + 60 * idx + 2);

            auto s_xx_xyy = primBuffer.data(skoff + 60 * idx + 3);

            auto s_xx_xyz = primBuffer.data(skoff + 60 * idx + 4);

            auto s_xx_xzz = primBuffer.data(skoff + 60 * idx + 5);

            auto s_xx_yyy = primBuffer.data(skoff + 60 * idx + 6);

            auto s_xx_yyz = primBuffer.data(skoff + 60 * idx + 7);

            auto s_xx_yzz = primBuffer.data(skoff + 60 * idx + 8);

            auto s_xx_zzz = primBuffer.data(skoff + 60 * idx + 9);

            auto s_xy_xxx = primBuffer.data(skoff + 60 * idx + 10);

            auto s_xy_xxy = primBuffer.data(skoff + 60 * idx + 11);

            auto s_xy_xxz = primBuffer.data(skoff + 60 * idx + 12);

            auto s_xy_xyy = primBuffer.data(skoff + 60 * idx + 13);

            auto s_xy_xyz = primBuffer.data(skoff + 60 * idx + 14);

            auto s_xy_xzz = primBuffer.data(skoff + 60 * idx + 15);

            auto s_xy_yyy = primBuffer.data(skoff + 60 * idx + 16);

            auto s_xy_yyz = primBuffer.data(skoff + 60 * idx + 17);

            auto s_xy_yzz = primBuffer.data(skoff + 60 * idx + 18);

            auto s_xy_zzz = primBuffer.data(skoff + 60 * idx + 19);

            auto s_xz_xxx = primBuffer.data(skoff + 60 * idx + 20);

            auto s_xz_xxy = primBuffer.data(skoff + 60 * idx + 21);

            auto s_xz_xxz = primBuffer.data(skoff + 60 * idx + 22);

            auto s_xz_xyy = primBuffer.data(skoff + 60 * idx + 23);

            auto s_xz_xyz = primBuffer.data(skoff + 60 * idx + 24);

            auto s_xz_xzz = primBuffer.data(skoff + 60 * idx + 25);

            auto s_xz_yyy = primBuffer.data(skoff + 60 * idx + 26);

            auto s_xz_yyz = primBuffer.data(skoff + 60 * idx + 27);

            auto s_xz_yzz = primBuffer.data(skoff + 60 * idx + 28);

            auto s_xz_zzz = primBuffer.data(skoff + 60 * idx + 29);

            auto s_yy_xxx = primBuffer.data(skoff + 60 * idx + 30);

            auto s_yy_xxy = primBuffer.data(skoff + 60 * idx + 31);

            auto s_yy_xxz = primBuffer.data(skoff + 60 * idx + 32);

            auto s_yy_xyy = primBuffer.data(skoff + 60 * idx + 33);

            auto s_yy_xyz = primBuffer.data(skoff + 60 * idx + 34);

            auto s_yy_xzz = primBuffer.data(skoff + 60 * idx + 35);

            auto s_yy_yyy = primBuffer.data(skoff + 60 * idx + 36);

            auto s_yy_yyz = primBuffer.data(skoff + 60 * idx + 37);

            auto s_yy_yzz = primBuffer.data(skoff + 60 * idx + 38);

            auto s_yy_zzz = primBuffer.data(skoff + 60 * idx + 39);

            auto s_yz_xxx = primBuffer.data(skoff + 60 * idx + 40);

            auto s_yz_xxy = primBuffer.data(skoff + 60 * idx + 41);

            auto s_yz_xxz = primBuffer.data(skoff + 60 * idx + 42);

            auto s_yz_xyy = primBuffer.data(skoff + 60 * idx + 43);

            auto s_yz_xyz = primBuffer.data(skoff + 60 * idx + 44);

            auto s_yz_xzz = primBuffer.data(skoff + 60 * idx + 45);

            auto s_yz_yyy = primBuffer.data(skoff + 60 * idx + 46);

            auto s_yz_yyz = primBuffer.data(skoff + 60 * idx + 47);

            auto s_yz_yzz = primBuffer.data(skoff + 60 * idx + 48);

            auto s_yz_zzz = primBuffer.data(skoff + 60 * idx + 49);

            auto s_zz_xxx = primBuffer.data(skoff + 60 * idx + 50);

            auto s_zz_xxy = primBuffer.data(skoff + 60 * idx + 51);

            auto s_zz_xxz = primBuffer.data(skoff + 60 * idx + 52);

            auto s_zz_xyy = primBuffer.data(skoff + 60 * idx + 53);

            auto s_zz_xyz = primBuffer.data(skoff + 60 * idx + 54);

            auto s_zz_xzz = primBuffer.data(skoff + 60 * idx + 55);

            auto s_zz_yyy = primBuffer.data(skoff + 60 * idx + 56);

            auto s_zz_yyz = primBuffer.data(skoff + 60 * idx + 57);

            auto s_zz_yzz = primBuffer.data(skoff + 60 * idx + 58);

            auto s_zz_zzz = primBuffer.data(skoff + 60 * idx + 59);

            // set up pointers to (D|T|F) integrals

            auto t_xx_xxx = primBuffer.data(kkoff + 60 * idx);

            auto t_xx_xxy = primBuffer.data(kkoff + 60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(kkoff + 60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(kkoff + 60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(kkoff + 60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(kkoff + 60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(kkoff + 60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(kkoff + 60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(kkoff + 60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(kkoff + 60 * idx + 9);

            auto t_xy_xxx = primBuffer.data(kkoff + 60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(kkoff + 60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(kkoff + 60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(kkoff + 60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(kkoff + 60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(kkoff + 60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(kkoff + 60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(kkoff + 60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(kkoff + 60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(kkoff + 60 * idx + 19);

            auto t_xz_xxx = primBuffer.data(kkoff + 60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(kkoff + 60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(kkoff + 60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(kkoff + 60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(kkoff + 60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(kkoff + 60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(kkoff + 60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(kkoff + 60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(kkoff + 60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(kkoff + 60 * idx + 29);

            auto t_yy_xxx = primBuffer.data(kkoff + 60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(kkoff + 60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(kkoff + 60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(kkoff + 60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(kkoff + 60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(kkoff + 60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(kkoff + 60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(kkoff + 60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(kkoff + 60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(kkoff + 60 * idx + 39);

            auto t_yz_xxx = primBuffer.data(kkoff + 60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(kkoff + 60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(kkoff + 60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(kkoff + 60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(kkoff + 60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(kkoff + 60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(kkoff + 60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(kkoff + 60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(kkoff + 60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(kkoff + 60 * idx + 49);

            auto t_zz_xxx = primBuffer.data(kkoff + 60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(kkoff + 60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(kkoff + 60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(kkoff + 60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(kkoff + 60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(kkoff + 60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(kkoff + 60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(kkoff + 60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(kkoff + 60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(kkoff + 60 * idx + 59);

            // set up pointers to (P|G) integrals

            auto s_x_xxxx = primBuffer.data(s2off + 45 * idx);

            auto s_x_xxxy = primBuffer.data(s2off + 45 * idx + 1);

            auto s_x_xxxz = primBuffer.data(s2off + 45 * idx + 2);

            auto s_x_xxyy = primBuffer.data(s2off + 45 * idx + 3);

            auto s_x_xxyz = primBuffer.data(s2off + 45 * idx + 4);

            auto s_x_xxzz = primBuffer.data(s2off + 45 * idx + 5);

            auto s_x_xyyy = primBuffer.data(s2off + 45 * idx + 6);

            auto s_x_xyyz = primBuffer.data(s2off + 45 * idx + 7);

            auto s_x_xyzz = primBuffer.data(s2off + 45 * idx + 8);

            auto s_x_xzzz = primBuffer.data(s2off + 45 * idx + 9);

            auto s_x_yyyy = primBuffer.data(s2off + 45 * idx + 10);

            auto s_x_yyyz = primBuffer.data(s2off + 45 * idx + 11);

            auto s_x_yyzz = primBuffer.data(s2off + 45 * idx + 12);

            auto s_x_yzzz = primBuffer.data(s2off + 45 * idx + 13);

            auto s_x_zzzz = primBuffer.data(s2off + 45 * idx + 14);

            auto s_y_xxxx = primBuffer.data(s2off + 45 * idx + 15);

            auto s_y_xxxy = primBuffer.data(s2off + 45 * idx + 16);

            auto s_y_xxxz = primBuffer.data(s2off + 45 * idx + 17);

            auto s_y_xxyy = primBuffer.data(s2off + 45 * idx + 18);

            auto s_y_xxyz = primBuffer.data(s2off + 45 * idx + 19);

            auto s_y_xxzz = primBuffer.data(s2off + 45 * idx + 20);

            auto s_y_xyyy = primBuffer.data(s2off + 45 * idx + 21);

            auto s_y_xyyz = primBuffer.data(s2off + 45 * idx + 22);

            auto s_y_xyzz = primBuffer.data(s2off + 45 * idx + 23);

            auto s_y_xzzz = primBuffer.data(s2off + 45 * idx + 24);

            auto s_y_yyyy = primBuffer.data(s2off + 45 * idx + 25);

            auto s_y_yyyz = primBuffer.data(s2off + 45 * idx + 26);

            auto s_y_yyzz = primBuffer.data(s2off + 45 * idx + 27);

            auto s_y_yzzz = primBuffer.data(s2off + 45 * idx + 28);

            auto s_y_zzzz = primBuffer.data(s2off + 45 * idx + 29);

            auto s_z_xxxx = primBuffer.data(s2off + 45 * idx + 30);

            auto s_z_xxxy = primBuffer.data(s2off + 45 * idx + 31);

            auto s_z_xxxz = primBuffer.data(s2off + 45 * idx + 32);

            auto s_z_xxyy = primBuffer.data(s2off + 45 * idx + 33);

            auto s_z_xxyz = primBuffer.data(s2off + 45 * idx + 34);

            auto s_z_xxzz = primBuffer.data(s2off + 45 * idx + 35);

            auto s_z_xyyy = primBuffer.data(s2off + 45 * idx + 36);

            auto s_z_xyyz = primBuffer.data(s2off + 45 * idx + 37);

            auto s_z_xyzz = primBuffer.data(s2off + 45 * idx + 38);

            auto s_z_xzzz = primBuffer.data(s2off + 45 * idx + 39);

            auto s_z_yyyy = primBuffer.data(s2off + 45 * idx + 40);

            auto s_z_yyyz = primBuffer.data(s2off + 45 * idx + 41);

            auto s_z_yyzz = primBuffer.data(s2off + 45 * idx + 42);

            auto s_z_yzzz = primBuffer.data(s2off + 45 * idx + 43);

            auto s_z_zzzz = primBuffer.data(s2off + 45 * idx + 44);

            // set up pointers to (P|T|G) integrals

            auto t_x_xxxx = primBuffer.data(k2off + 45 * idx);

            auto t_x_xxxy = primBuffer.data(k2off + 45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(k2off + 45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(k2off + 45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(k2off + 45 * idx + 4);

            auto t_x_xxzz = primBuffer.data(k2off + 45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(k2off + 45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(k2off + 45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(k2off + 45 * idx + 8);

            auto t_x_xzzz = primBuffer.data(k2off + 45 * idx + 9);

            auto t_x_yyyy = primBuffer.data(k2off + 45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(k2off + 45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(k2off + 45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(k2off + 45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(k2off + 45 * idx + 14);

            auto t_y_xxxx = primBuffer.data(k2off + 45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(k2off + 45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(k2off + 45 * idx + 17);

            auto t_y_xxyy = primBuffer.data(k2off + 45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(k2off + 45 * idx + 19);

            auto t_y_xxzz = primBuffer.data(k2off + 45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(k2off + 45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(k2off + 45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(k2off + 45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(k2off + 45 * idx + 24);

            auto t_y_yyyy = primBuffer.data(k2off + 45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(k2off + 45 * idx + 26);

            auto t_y_yyzz = primBuffer.data(k2off + 45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(k2off + 45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(k2off + 45 * idx + 29);

            auto t_z_xxxx = primBuffer.data(k2off + 45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(k2off + 45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(k2off + 45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(k2off + 45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(k2off + 45 * idx + 34);

            auto t_z_xxzz = primBuffer.data(k2off + 45 * idx + 35);

            auto t_z_xyyy = primBuffer.data(k2off + 45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(k2off + 45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(k2off + 45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(k2off + 45 * idx + 39);

            auto t_z_yyyy = primBuffer.data(k2off + 45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(k2off + 45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(k2off + 45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(k2off + 45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(k2off + 45 * idx + 44);

            // set up pointers to (D|G) integrals

            auto s_xx_xxxx = primBuffer.data(s1off + 90 * idx);

            auto s_xx_xxxy = primBuffer.data(s1off + 90 * idx + 1);

            auto s_xx_xxxz = primBuffer.data(s1off + 90 * idx + 2);

            auto s_xx_xxyy = primBuffer.data(s1off + 90 * idx + 3);

            auto s_xx_xxyz = primBuffer.data(s1off + 90 * idx + 4);

            auto s_xx_xxzz = primBuffer.data(s1off + 90 * idx + 5);

            auto s_xx_xyyy = primBuffer.data(s1off + 90 * idx + 6);

            auto s_xx_xyyz = primBuffer.data(s1off + 90 * idx + 7);

            auto s_xx_xyzz = primBuffer.data(s1off + 90 * idx + 8);

            auto s_xx_xzzz = primBuffer.data(s1off + 90 * idx + 9);

            auto s_xx_yyyy = primBuffer.data(s1off + 90 * idx + 10);

            auto s_xx_yyyz = primBuffer.data(s1off + 90 * idx + 11);

            auto s_xx_yyzz = primBuffer.data(s1off + 90 * idx + 12);

            auto s_xx_yzzz = primBuffer.data(s1off + 90 * idx + 13);

            auto s_xx_zzzz = primBuffer.data(s1off + 90 * idx + 14);

            auto s_xy_xxxx = primBuffer.data(s1off + 90 * idx + 15);

            auto s_xy_xxxy = primBuffer.data(s1off + 90 * idx + 16);

            auto s_xy_xxxz = primBuffer.data(s1off + 90 * idx + 17);

            auto s_xy_xxyy = primBuffer.data(s1off + 90 * idx + 18);

            auto s_xy_xxyz = primBuffer.data(s1off + 90 * idx + 19);

            auto s_xy_xxzz = primBuffer.data(s1off + 90 * idx + 20);

            auto s_xy_xyyy = primBuffer.data(s1off + 90 * idx + 21);

            auto s_xy_xyyz = primBuffer.data(s1off + 90 * idx + 22);

            auto s_xy_xyzz = primBuffer.data(s1off + 90 * idx + 23);

            auto s_xy_xzzz = primBuffer.data(s1off + 90 * idx + 24);

            auto s_xy_yyyy = primBuffer.data(s1off + 90 * idx + 25);

            auto s_xy_yyyz = primBuffer.data(s1off + 90 * idx + 26);

            auto s_xy_yyzz = primBuffer.data(s1off + 90 * idx + 27);

            auto s_xy_yzzz = primBuffer.data(s1off + 90 * idx + 28);

            auto s_xy_zzzz = primBuffer.data(s1off + 90 * idx + 29);

            auto s_xz_xxxx = primBuffer.data(s1off + 90 * idx + 30);

            auto s_xz_xxxy = primBuffer.data(s1off + 90 * idx + 31);

            auto s_xz_xxxz = primBuffer.data(s1off + 90 * idx + 32);

            auto s_xz_xxyy = primBuffer.data(s1off + 90 * idx + 33);

            auto s_xz_xxyz = primBuffer.data(s1off + 90 * idx + 34);

            auto s_xz_xxzz = primBuffer.data(s1off + 90 * idx + 35);

            auto s_xz_xyyy = primBuffer.data(s1off + 90 * idx + 36);

            auto s_xz_xyyz = primBuffer.data(s1off + 90 * idx + 37);

            auto s_xz_xyzz = primBuffer.data(s1off + 90 * idx + 38);

            auto s_xz_xzzz = primBuffer.data(s1off + 90 * idx + 39);

            auto s_xz_yyyy = primBuffer.data(s1off + 90 * idx + 40);

            auto s_xz_yyyz = primBuffer.data(s1off + 90 * idx + 41);

            auto s_xz_yyzz = primBuffer.data(s1off + 90 * idx + 42);

            auto s_xz_yzzz = primBuffer.data(s1off + 90 * idx + 43);

            auto s_xz_zzzz = primBuffer.data(s1off + 90 * idx + 44);

            auto s_yy_xxxx = primBuffer.data(s1off + 90 * idx + 45);

            auto s_yy_xxxy = primBuffer.data(s1off + 90 * idx + 46);

            auto s_yy_xxxz = primBuffer.data(s1off + 90 * idx + 47);

            auto s_yy_xxyy = primBuffer.data(s1off + 90 * idx + 48);

            auto s_yy_xxyz = primBuffer.data(s1off + 90 * idx + 49);

            auto s_yy_xxzz = primBuffer.data(s1off + 90 * idx + 50);

            auto s_yy_xyyy = primBuffer.data(s1off + 90 * idx + 51);

            auto s_yy_xyyz = primBuffer.data(s1off + 90 * idx + 52);

            auto s_yy_xyzz = primBuffer.data(s1off + 90 * idx + 53);

            auto s_yy_xzzz = primBuffer.data(s1off + 90 * idx + 54);

            auto s_yy_yyyy = primBuffer.data(s1off + 90 * idx + 55);

            auto s_yy_yyyz = primBuffer.data(s1off + 90 * idx + 56);

            auto s_yy_yyzz = primBuffer.data(s1off + 90 * idx + 57);

            auto s_yy_yzzz = primBuffer.data(s1off + 90 * idx + 58);

            auto s_yy_zzzz = primBuffer.data(s1off + 90 * idx + 59);

            auto s_yz_xxxx = primBuffer.data(s1off + 90 * idx + 60);

            auto s_yz_xxxy = primBuffer.data(s1off + 90 * idx + 61);

            auto s_yz_xxxz = primBuffer.data(s1off + 90 * idx + 62);

            auto s_yz_xxyy = primBuffer.data(s1off + 90 * idx + 63);

            auto s_yz_xxyz = primBuffer.data(s1off + 90 * idx + 64);

            auto s_yz_xxzz = primBuffer.data(s1off + 90 * idx + 65);

            auto s_yz_xyyy = primBuffer.data(s1off + 90 * idx + 66);

            auto s_yz_xyyz = primBuffer.data(s1off + 90 * idx + 67);

            auto s_yz_xyzz = primBuffer.data(s1off + 90 * idx + 68);

            auto s_yz_xzzz = primBuffer.data(s1off + 90 * idx + 69);

            auto s_yz_yyyy = primBuffer.data(s1off + 90 * idx + 70);

            auto s_yz_yyyz = primBuffer.data(s1off + 90 * idx + 71);

            auto s_yz_yyzz = primBuffer.data(s1off + 90 * idx + 72);

            auto s_yz_yzzz = primBuffer.data(s1off + 90 * idx + 73);

            auto s_yz_zzzz = primBuffer.data(s1off + 90 * idx + 74);

            auto s_zz_xxxx = primBuffer.data(s1off + 90 * idx + 75);

            auto s_zz_xxxy = primBuffer.data(s1off + 90 * idx + 76);

            auto s_zz_xxxz = primBuffer.data(s1off + 90 * idx + 77);

            auto s_zz_xxyy = primBuffer.data(s1off + 90 * idx + 78);

            auto s_zz_xxyz = primBuffer.data(s1off + 90 * idx + 79);

            auto s_zz_xxzz = primBuffer.data(s1off + 90 * idx + 80);

            auto s_zz_xyyy = primBuffer.data(s1off + 90 * idx + 81);

            auto s_zz_xyyz = primBuffer.data(s1off + 90 * idx + 82);

            auto s_zz_xyzz = primBuffer.data(s1off + 90 * idx + 83);

            auto s_zz_xzzz = primBuffer.data(s1off + 90 * idx + 84);

            auto s_zz_yyyy = primBuffer.data(s1off + 90 * idx + 85);

            auto s_zz_yyyz = primBuffer.data(s1off + 90 * idx + 86);

            auto s_zz_yyzz = primBuffer.data(s1off + 90 * idx + 87);

            auto s_zz_yzzz = primBuffer.data(s1off + 90 * idx + 88);

            auto s_zz_zzzz = primBuffer.data(s1off + 90 * idx + 89);

            // set up pointers to (D|T|G) integrals

            auto t_xx_xxxx = primBuffer.data(k1off + 90 * idx);

            auto t_xx_xxxy = primBuffer.data(k1off + 90 * idx + 1);

            auto t_xx_xxxz = primBuffer.data(k1off + 90 * idx + 2);

            auto t_xx_xxyy = primBuffer.data(k1off + 90 * idx + 3);

            auto t_xx_xxyz = primBuffer.data(k1off + 90 * idx + 4);

            auto t_xx_xxzz = primBuffer.data(k1off + 90 * idx + 5);

            auto t_xx_xyyy = primBuffer.data(k1off + 90 * idx + 6);

            auto t_xx_xyyz = primBuffer.data(k1off + 90 * idx + 7);

            auto t_xx_xyzz = primBuffer.data(k1off + 90 * idx + 8);

            auto t_xx_xzzz = primBuffer.data(k1off + 90 * idx + 9);

            auto t_xx_yyyy = primBuffer.data(k1off + 90 * idx + 10);

            auto t_xx_yyyz = primBuffer.data(k1off + 90 * idx + 11);

            auto t_xx_yyzz = primBuffer.data(k1off + 90 * idx + 12);

            auto t_xx_yzzz = primBuffer.data(k1off + 90 * idx + 13);

            auto t_xx_zzzz = primBuffer.data(k1off + 90 * idx + 14);

            auto t_xy_xxxx = primBuffer.data(k1off + 90 * idx + 15);

            auto t_xy_xxxy = primBuffer.data(k1off + 90 * idx + 16);

            auto t_xy_xxxz = primBuffer.data(k1off + 90 * idx + 17);

            auto t_xy_xxyy = primBuffer.data(k1off + 90 * idx + 18);

            auto t_xy_xxyz = primBuffer.data(k1off + 90 * idx + 19);

            auto t_xy_xxzz = primBuffer.data(k1off + 90 * idx + 20);

            auto t_xy_xyyy = primBuffer.data(k1off + 90 * idx + 21);

            auto t_xy_xyyz = primBuffer.data(k1off + 90 * idx + 22);

            auto t_xy_xyzz = primBuffer.data(k1off + 90 * idx + 23);

            auto t_xy_xzzz = primBuffer.data(k1off + 90 * idx + 24);

            auto t_xy_yyyy = primBuffer.data(k1off + 90 * idx + 25);

            auto t_xy_yyyz = primBuffer.data(k1off + 90 * idx + 26);

            auto t_xy_yyzz = primBuffer.data(k1off + 90 * idx + 27);

            auto t_xy_yzzz = primBuffer.data(k1off + 90 * idx + 28);

            auto t_xy_zzzz = primBuffer.data(k1off + 90 * idx + 29);

            auto t_xz_xxxx = primBuffer.data(k1off + 90 * idx + 30);

            auto t_xz_xxxy = primBuffer.data(k1off + 90 * idx + 31);

            auto t_xz_xxxz = primBuffer.data(k1off + 90 * idx + 32);

            auto t_xz_xxyy = primBuffer.data(k1off + 90 * idx + 33);

            auto t_xz_xxyz = primBuffer.data(k1off + 90 * idx + 34);

            auto t_xz_xxzz = primBuffer.data(k1off + 90 * idx + 35);

            auto t_xz_xyyy = primBuffer.data(k1off + 90 * idx + 36);

            auto t_xz_xyyz = primBuffer.data(k1off + 90 * idx + 37);

            auto t_xz_xyzz = primBuffer.data(k1off + 90 * idx + 38);

            auto t_xz_xzzz = primBuffer.data(k1off + 90 * idx + 39);

            auto t_xz_yyyy = primBuffer.data(k1off + 90 * idx + 40);

            auto t_xz_yyyz = primBuffer.data(k1off + 90 * idx + 41);

            auto t_xz_yyzz = primBuffer.data(k1off + 90 * idx + 42);

            auto t_xz_yzzz = primBuffer.data(k1off + 90 * idx + 43);

            auto t_xz_zzzz = primBuffer.data(k1off + 90 * idx + 44);

            auto t_yy_xxxx = primBuffer.data(k1off + 90 * idx + 45);

            auto t_yy_xxxy = primBuffer.data(k1off + 90 * idx + 46);

            auto t_yy_xxxz = primBuffer.data(k1off + 90 * idx + 47);

            auto t_yy_xxyy = primBuffer.data(k1off + 90 * idx + 48);

            auto t_yy_xxyz = primBuffer.data(k1off + 90 * idx + 49);

            auto t_yy_xxzz = primBuffer.data(k1off + 90 * idx + 50);

            auto t_yy_xyyy = primBuffer.data(k1off + 90 * idx + 51);

            auto t_yy_xyyz = primBuffer.data(k1off + 90 * idx + 52);

            auto t_yy_xyzz = primBuffer.data(k1off + 90 * idx + 53);

            auto t_yy_xzzz = primBuffer.data(k1off + 90 * idx + 54);

            auto t_yy_yyyy = primBuffer.data(k1off + 90 * idx + 55);

            auto t_yy_yyyz = primBuffer.data(k1off + 90 * idx + 56);

            auto t_yy_yyzz = primBuffer.data(k1off + 90 * idx + 57);

            auto t_yy_yzzz = primBuffer.data(k1off + 90 * idx + 58);

            auto t_yy_zzzz = primBuffer.data(k1off + 90 * idx + 59);

            auto t_yz_xxxx = primBuffer.data(k1off + 90 * idx + 60);

            auto t_yz_xxxy = primBuffer.data(k1off + 90 * idx + 61);

            auto t_yz_xxxz = primBuffer.data(k1off + 90 * idx + 62);

            auto t_yz_xxyy = primBuffer.data(k1off + 90 * idx + 63);

            auto t_yz_xxyz = primBuffer.data(k1off + 90 * idx + 64);

            auto t_yz_xxzz = primBuffer.data(k1off + 90 * idx + 65);

            auto t_yz_xyyy = primBuffer.data(k1off + 90 * idx + 66);

            auto t_yz_xyyz = primBuffer.data(k1off + 90 * idx + 67);

            auto t_yz_xyzz = primBuffer.data(k1off + 90 * idx + 68);

            auto t_yz_xzzz = primBuffer.data(k1off + 90 * idx + 69);

            auto t_yz_yyyy = primBuffer.data(k1off + 90 * idx + 70);

            auto t_yz_yyyz = primBuffer.data(k1off + 90 * idx + 71);

            auto t_yz_yyzz = primBuffer.data(k1off + 90 * idx + 72);

            auto t_yz_yzzz = primBuffer.data(k1off + 90 * idx + 73);

            auto t_yz_zzzz = primBuffer.data(k1off + 90 * idx + 74);

            auto t_zz_xxxx = primBuffer.data(k1off + 90 * idx + 75);

            auto t_zz_xxxy = primBuffer.data(k1off + 90 * idx + 76);

            auto t_zz_xxxz = primBuffer.data(k1off + 90 * idx + 77);

            auto t_zz_xxyy = primBuffer.data(k1off + 90 * idx + 78);

            auto t_zz_xxyz = primBuffer.data(k1off + 90 * idx + 79);

            auto t_zz_xxzz = primBuffer.data(k1off + 90 * idx + 80);

            auto t_zz_xyyy = primBuffer.data(k1off + 90 * idx + 81);

            auto t_zz_xyyz = primBuffer.data(k1off + 90 * idx + 82);

            auto t_zz_xyzz = primBuffer.data(k1off + 90 * idx + 83);

            auto t_zz_xzzz = primBuffer.data(k1off + 90 * idx + 84);

            auto t_zz_yyyy = primBuffer.data(k1off + 90 * idx + 85);

            auto t_zz_yyyz = primBuffer.data(k1off + 90 * idx + 86);

            auto t_zz_yyzz = primBuffer.data(k1off + 90 * idx + 87);

            auto t_zz_yzzz = primBuffer.data(k1off + 90 * idx + 88);

            auto t_zz_zzzz = primBuffer.data(k1off + 90 * idx + 89);

            // set up pointers to (F|G) integrals

            auto s_xxx_xxxx = primBuffer.data(soff + 150 * idx);

            auto s_xxx_xxxy = primBuffer.data(soff + 150 * idx + 1);

            auto s_xxx_xxxz = primBuffer.data(soff + 150 * idx + 2);

            auto s_xxx_xxyy = primBuffer.data(soff + 150 * idx + 3);

            auto s_xxx_xxyz = primBuffer.data(soff + 150 * idx + 4);

            auto s_xxx_xxzz = primBuffer.data(soff + 150 * idx + 5);

            auto s_xxx_xyyy = primBuffer.data(soff + 150 * idx + 6);

            auto s_xxx_xyyz = primBuffer.data(soff + 150 * idx + 7);

            auto s_xxx_xyzz = primBuffer.data(soff + 150 * idx + 8);

            auto s_xxx_xzzz = primBuffer.data(soff + 150 * idx + 9);

            auto s_xxx_yyyy = primBuffer.data(soff + 150 * idx + 10);

            auto s_xxx_yyyz = primBuffer.data(soff + 150 * idx + 11);

            auto s_xxx_yyzz = primBuffer.data(soff + 150 * idx + 12);

            auto s_xxx_yzzz = primBuffer.data(soff + 150 * idx + 13);

            auto s_xxx_zzzz = primBuffer.data(soff + 150 * idx + 14);

            auto s_xxy_xxxx = primBuffer.data(soff + 150 * idx + 15);

            auto s_xxy_xxxy = primBuffer.data(soff + 150 * idx + 16);

            auto s_xxy_xxxz = primBuffer.data(soff + 150 * idx + 17);

            auto s_xxy_xxyy = primBuffer.data(soff + 150 * idx + 18);

            auto s_xxy_xxyz = primBuffer.data(soff + 150 * idx + 19);

            auto s_xxy_xxzz = primBuffer.data(soff + 150 * idx + 20);

            auto s_xxy_xyyy = primBuffer.data(soff + 150 * idx + 21);

            auto s_xxy_xyyz = primBuffer.data(soff + 150 * idx + 22);

            auto s_xxy_xyzz = primBuffer.data(soff + 150 * idx + 23);

            auto s_xxy_xzzz = primBuffer.data(soff + 150 * idx + 24);

            auto s_xxy_yyyy = primBuffer.data(soff + 150 * idx + 25);

            auto s_xxy_yyyz = primBuffer.data(soff + 150 * idx + 26);

            auto s_xxy_yyzz = primBuffer.data(soff + 150 * idx + 27);

            auto s_xxy_yzzz = primBuffer.data(soff + 150 * idx + 28);

            auto s_xxy_zzzz = primBuffer.data(soff + 150 * idx + 29);

            auto s_xxz_xxxx = primBuffer.data(soff + 150 * idx + 30);

            auto s_xxz_xxxy = primBuffer.data(soff + 150 * idx + 31);

            auto s_xxz_xxxz = primBuffer.data(soff + 150 * idx + 32);

            auto s_xxz_xxyy = primBuffer.data(soff + 150 * idx + 33);

            auto s_xxz_xxyz = primBuffer.data(soff + 150 * idx + 34);

            auto s_xxz_xxzz = primBuffer.data(soff + 150 * idx + 35);

            auto s_xxz_xyyy = primBuffer.data(soff + 150 * idx + 36);

            auto s_xxz_xyyz = primBuffer.data(soff + 150 * idx + 37);

            auto s_xxz_xyzz = primBuffer.data(soff + 150 * idx + 38);

            auto s_xxz_xzzz = primBuffer.data(soff + 150 * idx + 39);

            auto s_xxz_yyyy = primBuffer.data(soff + 150 * idx + 40);

            auto s_xxz_yyyz = primBuffer.data(soff + 150 * idx + 41);

            auto s_xxz_yyzz = primBuffer.data(soff + 150 * idx + 42);

            auto s_xxz_yzzz = primBuffer.data(soff + 150 * idx + 43);

            auto s_xxz_zzzz = primBuffer.data(soff + 150 * idx + 44);

            auto s_xyy_xxxx = primBuffer.data(soff + 150 * idx + 45);

            auto s_xyy_xxxy = primBuffer.data(soff + 150 * idx + 46);

            auto s_xyy_xxxz = primBuffer.data(soff + 150 * idx + 47);

            auto s_xyy_xxyy = primBuffer.data(soff + 150 * idx + 48);

            auto s_xyy_xxyz = primBuffer.data(soff + 150 * idx + 49);

            auto s_xyy_xxzz = primBuffer.data(soff + 150 * idx + 50);

            auto s_xyy_xyyy = primBuffer.data(soff + 150 * idx + 51);

            auto s_xyy_xyyz = primBuffer.data(soff + 150 * idx + 52);

            auto s_xyy_xyzz = primBuffer.data(soff + 150 * idx + 53);

            auto s_xyy_xzzz = primBuffer.data(soff + 150 * idx + 54);

            auto s_xyy_yyyy = primBuffer.data(soff + 150 * idx + 55);

            auto s_xyy_yyyz = primBuffer.data(soff + 150 * idx + 56);

            auto s_xyy_yyzz = primBuffer.data(soff + 150 * idx + 57);

            auto s_xyy_yzzz = primBuffer.data(soff + 150 * idx + 58);

            auto s_xyy_zzzz = primBuffer.data(soff + 150 * idx + 59);

            auto s_xyz_xxxx = primBuffer.data(soff + 150 * idx + 60);

            auto s_xyz_xxxy = primBuffer.data(soff + 150 * idx + 61);

            auto s_xyz_xxxz = primBuffer.data(soff + 150 * idx + 62);

            auto s_xyz_xxyy = primBuffer.data(soff + 150 * idx + 63);

            auto s_xyz_xxyz = primBuffer.data(soff + 150 * idx + 64);

            auto s_xyz_xxzz = primBuffer.data(soff + 150 * idx + 65);

            auto s_xyz_xyyy = primBuffer.data(soff + 150 * idx + 66);

            auto s_xyz_xyyz = primBuffer.data(soff + 150 * idx + 67);

            auto s_xyz_xyzz = primBuffer.data(soff + 150 * idx + 68);

            auto s_xyz_xzzz = primBuffer.data(soff + 150 * idx + 69);

            auto s_xyz_yyyy = primBuffer.data(soff + 150 * idx + 70);

            auto s_xyz_yyyz = primBuffer.data(soff + 150 * idx + 71);

            auto s_xyz_yyzz = primBuffer.data(soff + 150 * idx + 72);

            auto s_xyz_yzzz = primBuffer.data(soff + 150 * idx + 73);

            auto s_xyz_zzzz = primBuffer.data(soff + 150 * idx + 74);

            auto s_xzz_xxxx = primBuffer.data(soff + 150 * idx + 75);

            auto s_xzz_xxxy = primBuffer.data(soff + 150 * idx + 76);

            auto s_xzz_xxxz = primBuffer.data(soff + 150 * idx + 77);

            auto s_xzz_xxyy = primBuffer.data(soff + 150 * idx + 78);

            auto s_xzz_xxyz = primBuffer.data(soff + 150 * idx + 79);

            auto s_xzz_xxzz = primBuffer.data(soff + 150 * idx + 80);

            auto s_xzz_xyyy = primBuffer.data(soff + 150 * idx + 81);

            auto s_xzz_xyyz = primBuffer.data(soff + 150 * idx + 82);

            auto s_xzz_xyzz = primBuffer.data(soff + 150 * idx + 83);

            auto s_xzz_xzzz = primBuffer.data(soff + 150 * idx + 84);

            auto s_xzz_yyyy = primBuffer.data(soff + 150 * idx + 85);

            auto s_xzz_yyyz = primBuffer.data(soff + 150 * idx + 86);

            auto s_xzz_yyzz = primBuffer.data(soff + 150 * idx + 87);

            auto s_xzz_yzzz = primBuffer.data(soff + 150 * idx + 88);

            auto s_xzz_zzzz = primBuffer.data(soff + 150 * idx + 89);

            auto s_yyy_xxxx = primBuffer.data(soff + 150 * idx + 90);

            auto s_yyy_xxxy = primBuffer.data(soff + 150 * idx + 91);

            auto s_yyy_xxxz = primBuffer.data(soff + 150 * idx + 92);

            auto s_yyy_xxyy = primBuffer.data(soff + 150 * idx + 93);

            auto s_yyy_xxyz = primBuffer.data(soff + 150 * idx + 94);

            auto s_yyy_xxzz = primBuffer.data(soff + 150 * idx + 95);

            auto s_yyy_xyyy = primBuffer.data(soff + 150 * idx + 96);

            auto s_yyy_xyyz = primBuffer.data(soff + 150 * idx + 97);

            auto s_yyy_xyzz = primBuffer.data(soff + 150 * idx + 98);

            auto s_yyy_xzzz = primBuffer.data(soff + 150 * idx + 99);

            auto s_yyy_yyyy = primBuffer.data(soff + 150 * idx + 100);

            auto s_yyy_yyyz = primBuffer.data(soff + 150 * idx + 101);

            auto s_yyy_yyzz = primBuffer.data(soff + 150 * idx + 102);

            auto s_yyy_yzzz = primBuffer.data(soff + 150 * idx + 103);

            auto s_yyy_zzzz = primBuffer.data(soff + 150 * idx + 104);

            auto s_yyz_xxxx = primBuffer.data(soff + 150 * idx + 105);

            auto s_yyz_xxxy = primBuffer.data(soff + 150 * idx + 106);

            auto s_yyz_xxxz = primBuffer.data(soff + 150 * idx + 107);

            auto s_yyz_xxyy = primBuffer.data(soff + 150 * idx + 108);

            auto s_yyz_xxyz = primBuffer.data(soff + 150 * idx + 109);

            auto s_yyz_xxzz = primBuffer.data(soff + 150 * idx + 110);

            auto s_yyz_xyyy = primBuffer.data(soff + 150 * idx + 111);

            auto s_yyz_xyyz = primBuffer.data(soff + 150 * idx + 112);

            auto s_yyz_xyzz = primBuffer.data(soff + 150 * idx + 113);

            auto s_yyz_xzzz = primBuffer.data(soff + 150 * idx + 114);

            auto s_yyz_yyyy = primBuffer.data(soff + 150 * idx + 115);

            auto s_yyz_yyyz = primBuffer.data(soff + 150 * idx + 116);

            auto s_yyz_yyzz = primBuffer.data(soff + 150 * idx + 117);

            auto s_yyz_yzzz = primBuffer.data(soff + 150 * idx + 118);

            auto s_yyz_zzzz = primBuffer.data(soff + 150 * idx + 119);

            auto s_yzz_xxxx = primBuffer.data(soff + 150 * idx + 120);

            auto s_yzz_xxxy = primBuffer.data(soff + 150 * idx + 121);

            auto s_yzz_xxxz = primBuffer.data(soff + 150 * idx + 122);

            auto s_yzz_xxyy = primBuffer.data(soff + 150 * idx + 123);

            auto s_yzz_xxyz = primBuffer.data(soff + 150 * idx + 124);

            auto s_yzz_xxzz = primBuffer.data(soff + 150 * idx + 125);

            auto s_yzz_xyyy = primBuffer.data(soff + 150 * idx + 126);

            auto s_yzz_xyyz = primBuffer.data(soff + 150 * idx + 127);

            auto s_yzz_xyzz = primBuffer.data(soff + 150 * idx + 128);

            auto s_yzz_xzzz = primBuffer.data(soff + 150 * idx + 129);

            auto s_yzz_yyyy = primBuffer.data(soff + 150 * idx + 130);

            auto s_yzz_yyyz = primBuffer.data(soff + 150 * idx + 131);

            auto s_yzz_yyzz = primBuffer.data(soff + 150 * idx + 132);

            auto s_yzz_yzzz = primBuffer.data(soff + 150 * idx + 133);

            auto s_yzz_zzzz = primBuffer.data(soff + 150 * idx + 134);

            auto s_zzz_xxxx = primBuffer.data(soff + 150 * idx + 135);

            auto s_zzz_xxxy = primBuffer.data(soff + 150 * idx + 136);

            auto s_zzz_xxxz = primBuffer.data(soff + 150 * idx + 137);

            auto s_zzz_xxyy = primBuffer.data(soff + 150 * idx + 138);

            auto s_zzz_xxyz = primBuffer.data(soff + 150 * idx + 139);

            auto s_zzz_xxzz = primBuffer.data(soff + 150 * idx + 140);

            auto s_zzz_xyyy = primBuffer.data(soff + 150 * idx + 141);

            auto s_zzz_xyyz = primBuffer.data(soff + 150 * idx + 142);

            auto s_zzz_xyzz = primBuffer.data(soff + 150 * idx + 143);

            auto s_zzz_xzzz = primBuffer.data(soff + 150 * idx + 144);

            auto s_zzz_yyyy = primBuffer.data(soff + 150 * idx + 145);

            auto s_zzz_yyyz = primBuffer.data(soff + 150 * idx + 146);

            auto s_zzz_yyzz = primBuffer.data(soff + 150 * idx + 147);

            auto s_zzz_yzzz = primBuffer.data(soff + 150 * idx + 148);

            auto s_zzz_zzzz = primBuffer.data(soff + 150 * idx + 149);

            // set up pointers to (F|T|G) integrals

            auto t_xxx_xxxx = primBuffer.data(koff + 150 * idx);

            auto t_xxx_xxxy = primBuffer.data(koff + 150 * idx + 1);

            auto t_xxx_xxxz = primBuffer.data(koff + 150 * idx + 2);

            auto t_xxx_xxyy = primBuffer.data(koff + 150 * idx + 3);

            auto t_xxx_xxyz = primBuffer.data(koff + 150 * idx + 4);

            auto t_xxx_xxzz = primBuffer.data(koff + 150 * idx + 5);

            auto t_xxx_xyyy = primBuffer.data(koff + 150 * idx + 6);

            auto t_xxx_xyyz = primBuffer.data(koff + 150 * idx + 7);

            auto t_xxx_xyzz = primBuffer.data(koff + 150 * idx + 8);

            auto t_xxx_xzzz = primBuffer.data(koff + 150 * idx + 9);

            auto t_xxx_yyyy = primBuffer.data(koff + 150 * idx + 10);

            auto t_xxx_yyyz = primBuffer.data(koff + 150 * idx + 11);

            auto t_xxx_yyzz = primBuffer.data(koff + 150 * idx + 12);

            auto t_xxx_yzzz = primBuffer.data(koff + 150 * idx + 13);

            auto t_xxx_zzzz = primBuffer.data(koff + 150 * idx + 14);

            auto t_xxy_xxxx = primBuffer.data(koff + 150 * idx + 15);

            auto t_xxy_xxxy = primBuffer.data(koff + 150 * idx + 16);

            auto t_xxy_xxxz = primBuffer.data(koff + 150 * idx + 17);

            auto t_xxy_xxyy = primBuffer.data(koff + 150 * idx + 18);

            auto t_xxy_xxyz = primBuffer.data(koff + 150 * idx + 19);

            auto t_xxy_xxzz = primBuffer.data(koff + 150 * idx + 20);

            auto t_xxy_xyyy = primBuffer.data(koff + 150 * idx + 21);

            auto t_xxy_xyyz = primBuffer.data(koff + 150 * idx + 22);

            auto t_xxy_xyzz = primBuffer.data(koff + 150 * idx + 23);

            auto t_xxy_xzzz = primBuffer.data(koff + 150 * idx + 24);

            auto t_xxy_yyyy = primBuffer.data(koff + 150 * idx + 25);

            auto t_xxy_yyyz = primBuffer.data(koff + 150 * idx + 26);

            auto t_xxy_yyzz = primBuffer.data(koff + 150 * idx + 27);

            auto t_xxy_yzzz = primBuffer.data(koff + 150 * idx + 28);

            auto t_xxy_zzzz = primBuffer.data(koff + 150 * idx + 29);

            auto t_xxz_xxxx = primBuffer.data(koff + 150 * idx + 30);

            auto t_xxz_xxxy = primBuffer.data(koff + 150 * idx + 31);

            auto t_xxz_xxxz = primBuffer.data(koff + 150 * idx + 32);

            auto t_xxz_xxyy = primBuffer.data(koff + 150 * idx + 33);

            auto t_xxz_xxyz = primBuffer.data(koff + 150 * idx + 34);

            auto t_xxz_xxzz = primBuffer.data(koff + 150 * idx + 35);

            auto t_xxz_xyyy = primBuffer.data(koff + 150 * idx + 36);

            auto t_xxz_xyyz = primBuffer.data(koff + 150 * idx + 37);

            auto t_xxz_xyzz = primBuffer.data(koff + 150 * idx + 38);

            auto t_xxz_xzzz = primBuffer.data(koff + 150 * idx + 39);

            auto t_xxz_yyyy = primBuffer.data(koff + 150 * idx + 40);

            auto t_xxz_yyyz = primBuffer.data(koff + 150 * idx + 41);

            auto t_xxz_yyzz = primBuffer.data(koff + 150 * idx + 42);

            auto t_xxz_yzzz = primBuffer.data(koff + 150 * idx + 43);

            auto t_xxz_zzzz = primBuffer.data(koff + 150 * idx + 44);

            auto t_xyy_xxxx = primBuffer.data(koff + 150 * idx + 45);

            auto t_xyy_xxxy = primBuffer.data(koff + 150 * idx + 46);

            auto t_xyy_xxxz = primBuffer.data(koff + 150 * idx + 47);

            auto t_xyy_xxyy = primBuffer.data(koff + 150 * idx + 48);

            auto t_xyy_xxyz = primBuffer.data(koff + 150 * idx + 49);

            auto t_xyy_xxzz = primBuffer.data(koff + 150 * idx + 50);

            auto t_xyy_xyyy = primBuffer.data(koff + 150 * idx + 51);

            auto t_xyy_xyyz = primBuffer.data(koff + 150 * idx + 52);

            auto t_xyy_xyzz = primBuffer.data(koff + 150 * idx + 53);

            auto t_xyy_xzzz = primBuffer.data(koff + 150 * idx + 54);

            auto t_xyy_yyyy = primBuffer.data(koff + 150 * idx + 55);

            auto t_xyy_yyyz = primBuffer.data(koff + 150 * idx + 56);

            auto t_xyy_yyzz = primBuffer.data(koff + 150 * idx + 57);

            auto t_xyy_yzzz = primBuffer.data(koff + 150 * idx + 58);

            auto t_xyy_zzzz = primBuffer.data(koff + 150 * idx + 59);

            auto t_xyz_xxxx = primBuffer.data(koff + 150 * idx + 60);

            auto t_xyz_xxxy = primBuffer.data(koff + 150 * idx + 61);

            auto t_xyz_xxxz = primBuffer.data(koff + 150 * idx + 62);

            auto t_xyz_xxyy = primBuffer.data(koff + 150 * idx + 63);

            auto t_xyz_xxyz = primBuffer.data(koff + 150 * idx + 64);

            auto t_xyz_xxzz = primBuffer.data(koff + 150 * idx + 65);

            auto t_xyz_xyyy = primBuffer.data(koff + 150 * idx + 66);

            auto t_xyz_xyyz = primBuffer.data(koff + 150 * idx + 67);

            auto t_xyz_xyzz = primBuffer.data(koff + 150 * idx + 68);

            auto t_xyz_xzzz = primBuffer.data(koff + 150 * idx + 69);

            auto t_xyz_yyyy = primBuffer.data(koff + 150 * idx + 70);

            auto t_xyz_yyyz = primBuffer.data(koff + 150 * idx + 71);

            auto t_xyz_yyzz = primBuffer.data(koff + 150 * idx + 72);

            auto t_xyz_yzzz = primBuffer.data(koff + 150 * idx + 73);

            auto t_xyz_zzzz = primBuffer.data(koff + 150 * idx + 74);

            auto t_xzz_xxxx = primBuffer.data(koff + 150 * idx + 75);

            auto t_xzz_xxxy = primBuffer.data(koff + 150 * idx + 76);

            auto t_xzz_xxxz = primBuffer.data(koff + 150 * idx + 77);

            auto t_xzz_xxyy = primBuffer.data(koff + 150 * idx + 78);

            auto t_xzz_xxyz = primBuffer.data(koff + 150 * idx + 79);

            auto t_xzz_xxzz = primBuffer.data(koff + 150 * idx + 80);

            auto t_xzz_xyyy = primBuffer.data(koff + 150 * idx + 81);

            auto t_xzz_xyyz = primBuffer.data(koff + 150 * idx + 82);

            auto t_xzz_xyzz = primBuffer.data(koff + 150 * idx + 83);

            auto t_xzz_xzzz = primBuffer.data(koff + 150 * idx + 84);

            auto t_xzz_yyyy = primBuffer.data(koff + 150 * idx + 85);

            auto t_xzz_yyyz = primBuffer.data(koff + 150 * idx + 86);

            auto t_xzz_yyzz = primBuffer.data(koff + 150 * idx + 87);

            auto t_xzz_yzzz = primBuffer.data(koff + 150 * idx + 88);

            auto t_xzz_zzzz = primBuffer.data(koff + 150 * idx + 89);

            auto t_yyy_xxxx = primBuffer.data(koff + 150 * idx + 90);

            auto t_yyy_xxxy = primBuffer.data(koff + 150 * idx + 91);

            auto t_yyy_xxxz = primBuffer.data(koff + 150 * idx + 92);

            auto t_yyy_xxyy = primBuffer.data(koff + 150 * idx + 93);

            auto t_yyy_xxyz = primBuffer.data(koff + 150 * idx + 94);

            auto t_yyy_xxzz = primBuffer.data(koff + 150 * idx + 95);

            auto t_yyy_xyyy = primBuffer.data(koff + 150 * idx + 96);

            auto t_yyy_xyyz = primBuffer.data(koff + 150 * idx + 97);

            auto t_yyy_xyzz = primBuffer.data(koff + 150 * idx + 98);

            auto t_yyy_xzzz = primBuffer.data(koff + 150 * idx + 99);

            auto t_yyy_yyyy = primBuffer.data(koff + 150 * idx + 100);

            auto t_yyy_yyyz = primBuffer.data(koff + 150 * idx + 101);

            auto t_yyy_yyzz = primBuffer.data(koff + 150 * idx + 102);

            auto t_yyy_yzzz = primBuffer.data(koff + 150 * idx + 103);

            auto t_yyy_zzzz = primBuffer.data(koff + 150 * idx + 104);

            auto t_yyz_xxxx = primBuffer.data(koff + 150 * idx + 105);

            auto t_yyz_xxxy = primBuffer.data(koff + 150 * idx + 106);

            auto t_yyz_xxxz = primBuffer.data(koff + 150 * idx + 107);

            auto t_yyz_xxyy = primBuffer.data(koff + 150 * idx + 108);

            auto t_yyz_xxyz = primBuffer.data(koff + 150 * idx + 109);

            auto t_yyz_xxzz = primBuffer.data(koff + 150 * idx + 110);

            auto t_yyz_xyyy = primBuffer.data(koff + 150 * idx + 111);

            auto t_yyz_xyyz = primBuffer.data(koff + 150 * idx + 112);

            auto t_yyz_xyzz = primBuffer.data(koff + 150 * idx + 113);

            auto t_yyz_xzzz = primBuffer.data(koff + 150 * idx + 114);

            auto t_yyz_yyyy = primBuffer.data(koff + 150 * idx + 115);

            auto t_yyz_yyyz = primBuffer.data(koff + 150 * idx + 116);

            auto t_yyz_yyzz = primBuffer.data(koff + 150 * idx + 117);

            auto t_yyz_yzzz = primBuffer.data(koff + 150 * idx + 118);

            auto t_yyz_zzzz = primBuffer.data(koff + 150 * idx + 119);

            auto t_yzz_xxxx = primBuffer.data(koff + 150 * idx + 120);

            auto t_yzz_xxxy = primBuffer.data(koff + 150 * idx + 121);

            auto t_yzz_xxxz = primBuffer.data(koff + 150 * idx + 122);

            auto t_yzz_xxyy = primBuffer.data(koff + 150 * idx + 123);

            auto t_yzz_xxyz = primBuffer.data(koff + 150 * idx + 124);

            auto t_yzz_xxzz = primBuffer.data(koff + 150 * idx + 125);

            auto t_yzz_xyyy = primBuffer.data(koff + 150 * idx + 126);

            auto t_yzz_xyyz = primBuffer.data(koff + 150 * idx + 127);

            auto t_yzz_xyzz = primBuffer.data(koff + 150 * idx + 128);

            auto t_yzz_xzzz = primBuffer.data(koff + 150 * idx + 129);

            auto t_yzz_yyyy = primBuffer.data(koff + 150 * idx + 130);

            auto t_yzz_yyyz = primBuffer.data(koff + 150 * idx + 131);

            auto t_yzz_yyzz = primBuffer.data(koff + 150 * idx + 132);

            auto t_yzz_yzzz = primBuffer.data(koff + 150 * idx + 133);

            auto t_yzz_zzzz = primBuffer.data(koff + 150 * idx + 134);

            auto t_zzz_xxxx = primBuffer.data(koff + 150 * idx + 135);

            auto t_zzz_xxxy = primBuffer.data(koff + 150 * idx + 136);

            auto t_zzz_xxxz = primBuffer.data(koff + 150 * idx + 137);

            auto t_zzz_xxyy = primBuffer.data(koff + 150 * idx + 138);

            auto t_zzz_xxyz = primBuffer.data(koff + 150 * idx + 139);

            auto t_zzz_xxzz = primBuffer.data(koff + 150 * idx + 140);

            auto t_zzz_xyyy = primBuffer.data(koff + 150 * idx + 141);

            auto t_zzz_xyyz = primBuffer.data(koff + 150 * idx + 142);

            auto t_zzz_xyzz = primBuffer.data(koff + 150 * idx + 143);

            auto t_zzz_xzzz = primBuffer.data(koff + 150 * idx + 144);

            auto t_zzz_yyyy = primBuffer.data(koff + 150 * idx + 145);

            auto t_zzz_yyyz = primBuffer.data(koff + 150 * idx + 146);

            auto t_zzz_yyzz = primBuffer.data(koff + 150 * idx + 147);

            auto t_zzz_yzzz = primBuffer.data(koff + 150 * idx + 148);

            auto t_zzz_zzzz = primBuffer.data(koff + 150 * idx + 149);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xx_xxx, s_xx_xxy,\
                                     s_xx_xxz, s_xx_xyy, s_xx_xyz, s_xx_xzz, s_xx_yyy,\
                                     s_xx_yyz, s_xx_yzz, s_xx_zzz, s_xy_xxx, s_xy_xxy,\
                                     s_xy_xxz, s_xy_xyy, s_xy_xyz, s_xy_xzz, s_xy_yyy,\
                                     s_xy_yyz, s_xy_yzz, s_xy_zzz, s_xz_xxx, s_xz_xxy,\
                                     s_xz_xxz, s_xz_xyy, s_xz_xyz, s_xz_xzz, s_xz_yyy,\
                                     s_xz_yyz, s_xz_yzz, s_xz_zzz, s_yy_xxx, s_yy_xxy,\
                                     s_yy_xxz, s_yy_xyy, s_yy_xyz, s_yy_xzz, s_yy_yyy,\
                                     s_yy_yyz, s_yy_yzz, s_yy_zzz, s_yz_xxx, s_yz_xxy,\
                                     s_yz_xxz, s_yz_xyy, s_yz_xyz, s_yz_xzz, s_yz_yyy,\
                                     s_yz_yyz, s_yz_yzz, s_yz_zzz, s_zz_xxx, s_zz_xxy,\
                                     s_zz_xxz, s_zz_xyy, s_zz_xyz, s_zz_xzz, s_zz_yyy,\
                                     s_zz_yyz, s_zz_yzz, s_zz_zzz, t_xx_xxx, t_xx_xxy,\
                                     t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, t_xx_yyy,\
                                     t_xx_yyz, t_xx_yzz, t_xx_zzz, t_xy_xxx, t_xy_xxy,\
                                     t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy,\
                                     t_xy_yyz, t_xy_yzz, t_xy_zzz, t_xz_xxx, t_xz_xxy,\
                                     t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy,\
                                     t_xz_yyz, t_xz_yzz, t_xz_zzz, t_yy_xxx, t_yy_xxy,\
                                     t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, t_yy_yyy,\
                                     t_yy_yyz, t_yy_yzz, t_yy_zzz, t_yz_xxx, t_yz_xxy,\
                                     t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy,\
                                     t_yz_yyz, t_yz_yzz, t_yz_zzz, t_zz_xxx, t_zz_xxy,\
                                     t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, t_zz_yyy,\
                                     t_zz_yyz, t_zz_yzz, t_zz_zzz, s_x_xxxx, s_x_xxxy,\
                                     s_x_xxxz, s_x_xxyy, s_x_xxyz, s_x_xxzz, s_x_xyyy,\
                                     s_x_xyyz, s_x_xyzz, s_x_xzzz, s_x_yyyy, s_x_yyyz,\
                                     s_x_yyzz, s_x_yzzz, s_x_zzzz, s_y_xxxx, s_y_xxxy,\
                                     s_y_xxxz, s_y_xxyy, s_y_xxyz, s_y_xxzz, s_y_xyyy,\
                                     s_y_xyyz, s_y_xyzz, s_y_xzzz, s_y_yyyy, s_y_yyyz,\
                                     s_y_yyzz, s_y_yzzz, s_y_zzzz, s_z_xxxx, s_z_xxxy,\
                                     s_z_xxxz, s_z_xxyy, s_z_xxyz, s_z_xxzz, s_z_xyyy,\
                                     s_z_xyyz, s_z_xyzz, s_z_xzzz, s_z_yyyy, s_z_yyyz,\
                                     s_z_yyzz, s_z_yzzz, s_z_zzzz, t_x_xxxx, t_x_xxxy,\
                                     t_x_xxxz, t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy,\
                                     t_x_xyyz, t_x_xyzz, t_x_xzzz, t_x_yyyy, t_x_yyyz,\
                                     t_x_yyzz, t_x_yzzz, t_x_zzzz, t_y_xxxx, t_y_xxxy,\
                                     t_y_xxxz, t_y_xxyy, t_y_xxyz, t_y_xxzz, t_y_xyyy,\
                                     t_y_xyyz, t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz,\
                                     t_y_yyzz, t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy,\
                                     t_z_xxxz, t_z_xxyy, t_z_xxyz, t_z_xxzz, t_z_xyyy,\
                                     t_z_xyyz, t_z_xyzz, t_z_xzzz, t_z_yyyy, t_z_yyyz,\
                                     t_z_yyzz, t_z_yzzz, t_z_zzzz, s_xx_xxxx, s_xx_xxxy,\
                                     s_xx_xxxz, s_xx_xxyy, s_xx_xxyz, s_xx_xxzz,\
                                     s_xx_xyyy, s_xx_xyyz, s_xx_xyzz, s_xx_xzzz,\
                                     s_xx_yyyy, s_xx_yyyz, s_xx_yyzz, s_xx_yzzz,\
                                     s_xx_zzzz, s_xy_xxxx, s_xy_xxxy, s_xy_xxxz,\
                                     s_xy_xxyy, s_xy_xxyz, s_xy_xxzz, s_xy_xyyy,\
                                     s_xy_xyyz, s_xy_xyzz, s_xy_xzzz, s_xy_yyyy,\
                                     s_xy_yyyz, s_xy_yyzz, s_xy_yzzz, s_xy_zzzz,\
                                     s_xz_xxxx, s_xz_xxxy, s_xz_xxxz, s_xz_xxyy,\
                                     s_xz_xxyz, s_xz_xxzz, s_xz_xyyy, s_xz_xyyz,\
                                     s_xz_xyzz, s_xz_xzzz, s_xz_yyyy, s_xz_yyyz,\
                                     s_xz_yyzz, s_xz_yzzz, s_xz_zzzz, s_yy_xxxx,\
                                     s_yy_xxxy, s_yy_xxxz, s_yy_xxyy, s_yy_xxyz,\
                                     s_yy_xxzz, s_yy_xyyy, s_yy_xyyz, s_yy_xyzz,\
                                     s_yy_xzzz, s_yy_yyyy, s_yy_yyyz, s_yy_yyzz,\
                                     s_yy_yzzz, s_yy_zzzz, s_yz_xxxx, s_yz_xxxy,\
                                     s_yz_xxxz, s_yz_xxyy, s_yz_xxyz, s_yz_xxzz,\
                                     s_yz_xyyy, s_yz_xyyz, s_yz_xyzz, s_yz_xzzz,\
                                     s_yz_yyyy, s_yz_yyyz, s_yz_yyzz, s_yz_yzzz,\
                                     s_yz_zzzz, s_zz_xxxx, s_zz_xxxy, s_zz_xxxz,\
                                     s_zz_xxyy, s_zz_xxyz, s_zz_xxzz, s_zz_xyyy,\
                                     s_zz_xyyz, s_zz_xyzz, s_zz_xzzz, s_zz_yyyy,\
                                     s_zz_yyyz, s_zz_yyzz, s_zz_yzzz, s_zz_zzzz,\
                                     t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy,\
                                     t_xx_xxyz, t_xx_xxzz, t_xx_xyyy, t_xx_xyyz,\
                                     t_xx_xyzz, t_xx_xzzz, t_xx_yyyy, t_xx_yyyz,\
                                     t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx,\
                                     t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz,\
                                     t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz,\
                                     t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz,\
                                     t_xy_yzzz, t_xy_zzzz, t_xz_xxxx, t_xz_xxxy,\
                                     t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, t_xz_xxzz,\
                                     t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz,\
                                     t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz,\
                                     t_xz_zzzz, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz,\
                                     t_yy_xxyy, t_yy_xxyz, t_yy_xxzz, t_yy_xyyy,\
                                     t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy,\
                                     t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, t_yy_zzzz,\
                                     t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, t_yz_xxyy,\
                                     t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz,\
                                     t_yz_xyzz, t_yz_xzzz, t_yz_yyyy, t_yz_yyyz,\
                                     t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_zz_xxxx,\
                                     t_zz_xxxy, t_zz_xxxz, t_zz_xxyy, t_zz_xxyz,\
                                     t_zz_xxzz, t_zz_xyyy, t_zz_xyyz, t_zz_xyzz,\
                                     t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz,\
                                     t_zz_yzzz, t_zz_zzzz, s_xxx_xxxx, s_xxx_xxxy,\
                                     s_xxx_xxxz, s_xxx_xxyy, s_xxx_xxyz, s_xxx_xxzz,\
                                     s_xxx_xyyy, s_xxx_xyyz, s_xxx_xyzz, s_xxx_xzzz,\
                                     s_xxx_yyyy, s_xxx_yyyz, s_xxx_yyzz, s_xxx_yzzz,\
                                     s_xxx_zzzz, s_xxy_xxxx, s_xxy_xxxy, s_xxy_xxxz,\
                                     s_xxy_xxyy, s_xxy_xxyz, s_xxy_xxzz, s_xxy_xyyy,\
                                     s_xxy_xyyz, s_xxy_xyzz, s_xxy_xzzz, s_xxy_yyyy,\
                                     s_xxy_yyyz, s_xxy_yyzz, s_xxy_yzzz, s_xxy_zzzz,\
                                     s_xxz_xxxx, s_xxz_xxxy, s_xxz_xxxz, s_xxz_xxyy,\
                                     s_xxz_xxyz, s_xxz_xxzz, s_xxz_xyyy, s_xxz_xyyz,\
                                     s_xxz_xyzz, s_xxz_xzzz, s_xxz_yyyy, s_xxz_yyyz,\
                                     s_xxz_yyzz, s_xxz_yzzz, s_xxz_zzzz, s_xyy_xxxx,\
                                     s_xyy_xxxy, s_xyy_xxxz, s_xyy_xxyy, s_xyy_xxyz,\
                                     s_xyy_xxzz, s_xyy_xyyy, s_xyy_xyyz, s_xyy_xyzz,\
                                     s_xyy_xzzz, s_xyy_yyyy, s_xyy_yyyz, s_xyy_yyzz,\
                                     s_xyy_yzzz, s_xyy_zzzz, s_xyz_xxxx, s_xyz_xxxy,\
                                     s_xyz_xxxz, s_xyz_xxyy, s_xyz_xxyz, s_xyz_xxzz,\
                                     s_xyz_xyyy, s_xyz_xyyz, s_xyz_xyzz, s_xyz_xzzz,\
                                     s_xyz_yyyy, s_xyz_yyyz, s_xyz_yyzz, s_xyz_yzzz,\
                                     s_xyz_zzzz, s_xzz_xxxx, s_xzz_xxxy, s_xzz_xxxz,\
                                     s_xzz_xxyy, s_xzz_xxyz, s_xzz_xxzz, s_xzz_xyyy,\
                                     s_xzz_xyyz, s_xzz_xyzz, s_xzz_xzzz, s_xzz_yyyy,\
                                     s_xzz_yyyz, s_xzz_yyzz, s_xzz_yzzz, s_xzz_zzzz,\
                                     s_yyy_xxxx, s_yyy_xxxy, s_yyy_xxxz, s_yyy_xxyy,\
                                     s_yyy_xxyz, s_yyy_xxzz, s_yyy_xyyy, s_yyy_xyyz,\
                                     s_yyy_xyzz, s_yyy_xzzz, s_yyy_yyyy, s_yyy_yyyz,\
                                     s_yyy_yyzz, s_yyy_yzzz, s_yyy_zzzz, s_yyz_xxxx,\
                                     s_yyz_xxxy, s_yyz_xxxz, s_yyz_xxyy, s_yyz_xxyz,\
                                     s_yyz_xxzz, s_yyz_xyyy, s_yyz_xyyz, s_yyz_xyzz,\
                                     s_yyz_xzzz, s_yyz_yyyy, s_yyz_yyyz, s_yyz_yyzz,\
                                     s_yyz_yzzz, s_yyz_zzzz, s_yzz_xxxx, s_yzz_xxxy,\
                                     s_yzz_xxxz, s_yzz_xxyy, s_yzz_xxyz, s_yzz_xxzz,\
                                     s_yzz_xyyy, s_yzz_xyyz, s_yzz_xyzz, s_yzz_xzzz,\
                                     s_yzz_yyyy, s_yzz_yyyz, s_yzz_yyzz, s_yzz_yzzz,\
                                     s_yzz_zzzz, s_zzz_xxxx, s_zzz_xxxy, s_zzz_xxxz,\
                                     s_zzz_xxyy, s_zzz_xxyz, s_zzz_xxzz, s_zzz_xyyy,\
                                     s_zzz_xyyz, s_zzz_xyzz, s_zzz_xzzz, s_zzz_yyyy,\
                                     s_zzz_yyyz, s_zzz_yyzz, s_zzz_yzzz, s_zzz_zzzz,\
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
                                     t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxx_xxxx[j] = fr * s_xx_xxxx[j] + f2t * (2.0 * s_x_xxxx[j] + 4.0 * s_xx_xxx[j]);

                t_xxx_xxxx[j] = fr * t_xx_xxxx[j] + f2t * (2.0 * t_x_xxxx[j] + 4.0 * t_xx_xxx[j]) + f2z * (s_xxx_xxxx[j] - f2a * 2.0 * s_x_xxxx[j]);

                s_xxx_xxxy[j] = fr * s_xx_xxxy[j] + f2t * (2.0 * s_x_xxxy[j] + 3.0 * s_xx_xxy[j]);

                t_xxx_xxxy[j] = fr * t_xx_xxxy[j] + f2t * (2.0 * t_x_xxxy[j] + 3.0 * t_xx_xxy[j]) + f2z * (s_xxx_xxxy[j] - f2a * 2.0 * s_x_xxxy[j]);

                s_xxx_xxxz[j] = fr * s_xx_xxxz[j] + f2t * (2.0 * s_x_xxxz[j] + 3.0 * s_xx_xxz[j]);

                t_xxx_xxxz[j] = fr * t_xx_xxxz[j] + f2t * (2.0 * t_x_xxxz[j] + 3.0 * t_xx_xxz[j]) + f2z * (s_xxx_xxxz[j] - f2a * 2.0 * s_x_xxxz[j]);

                s_xxx_xxyy[j] = fr * s_xx_xxyy[j] + f2t * (2.0 * s_x_xxyy[j] + 2.0 * s_xx_xyy[j]);

                t_xxx_xxyy[j] = fr * t_xx_xxyy[j] + f2t * (2.0 * t_x_xxyy[j] + 2.0 * t_xx_xyy[j]) + f2z * (s_xxx_xxyy[j] - f2a * 2.0 * s_x_xxyy[j]);

                s_xxx_xxyz[j] = fr * s_xx_xxyz[j] + f2t * (2.0 * s_x_xxyz[j] + 2.0 * s_xx_xyz[j]);

                t_xxx_xxyz[j] = fr * t_xx_xxyz[j] + f2t * (2.0 * t_x_xxyz[j] + 2.0 * t_xx_xyz[j]) + f2z * (s_xxx_xxyz[j] - f2a * 2.0 * s_x_xxyz[j]);

                s_xxx_xxzz[j] = fr * s_xx_xxzz[j] + f2t * (2.0 * s_x_xxzz[j] + 2.0 * s_xx_xzz[j]);

                t_xxx_xxzz[j] = fr * t_xx_xxzz[j] + f2t * (2.0 * t_x_xxzz[j] + 2.0 * t_xx_xzz[j]) + f2z * (s_xxx_xxzz[j] - f2a * 2.0 * s_x_xxzz[j]);

                s_xxx_xyyy[j] = fr * s_xx_xyyy[j] + f2t * (2.0 * s_x_xyyy[j] + s_xx_yyy[j]);

                t_xxx_xyyy[j] = fr * t_xx_xyyy[j] + f2t * (2.0 * t_x_xyyy[j] + t_xx_yyy[j]) + f2z * (s_xxx_xyyy[j] - f2a * 2.0 * s_x_xyyy[j]);

                s_xxx_xyyz[j] = fr * s_xx_xyyz[j] + f2t * (2.0 * s_x_xyyz[j] + s_xx_yyz[j]);

                t_xxx_xyyz[j] = fr * t_xx_xyyz[j] + f2t * (2.0 * t_x_xyyz[j] + t_xx_yyz[j]) + f2z * (s_xxx_xyyz[j] - f2a * 2.0 * s_x_xyyz[j]);

                s_xxx_xyzz[j] = fr * s_xx_xyzz[j] + f2t * (2.0 * s_x_xyzz[j] + s_xx_yzz[j]);

                t_xxx_xyzz[j] = fr * t_xx_xyzz[j] + f2t * (2.0 * t_x_xyzz[j] + t_xx_yzz[j]) + f2z * (s_xxx_xyzz[j] - f2a * 2.0 * s_x_xyzz[j]);

                s_xxx_xzzz[j] = fr * s_xx_xzzz[j] + f2t * (2.0 * s_x_xzzz[j] + s_xx_zzz[j]);

                t_xxx_xzzz[j] = fr * t_xx_xzzz[j] + f2t * (2.0 * t_x_xzzz[j] + t_xx_zzz[j]) + f2z * (s_xxx_xzzz[j] - f2a * 2.0 * s_x_xzzz[j]);

                s_xxx_yyyy[j] = fr * s_xx_yyyy[j] + f2t * 2.0 * s_x_yyyy[j];

                t_xxx_yyyy[j] = fr * t_xx_yyyy[j] + f2t * 2.0 * t_x_yyyy[j] + f2z * (s_xxx_yyyy[j] - f2a * 2.0 * s_x_yyyy[j]);

                s_xxx_yyyz[j] = fr * s_xx_yyyz[j] + f2t * 2.0 * s_x_yyyz[j];

                t_xxx_yyyz[j] = fr * t_xx_yyyz[j] + f2t * 2.0 * t_x_yyyz[j] + f2z * (s_xxx_yyyz[j] - f2a * 2.0 * s_x_yyyz[j]);

                s_xxx_yyzz[j] = fr * s_xx_yyzz[j] + f2t * 2.0 * s_x_yyzz[j];

                t_xxx_yyzz[j] = fr * t_xx_yyzz[j] + f2t * 2.0 * t_x_yyzz[j] + f2z * (s_xxx_yyzz[j] - f2a * 2.0 * s_x_yyzz[j]);

                s_xxx_yzzz[j] = fr * s_xx_yzzz[j] + f2t * 2.0 * s_x_yzzz[j];

                t_xxx_yzzz[j] = fr * t_xx_yzzz[j] + f2t * 2.0 * t_x_yzzz[j] + f2z * (s_xxx_yzzz[j] - f2a * 2.0 * s_x_yzzz[j]);

                s_xxx_zzzz[j] = fr * s_xx_zzzz[j] + f2t * 2.0 * s_x_zzzz[j];

                t_xxx_zzzz[j] = fr * t_xx_zzzz[j] + f2t * 2.0 * t_x_zzzz[j] + f2z * (s_xxx_zzzz[j] - f2a * 2.0 * s_x_zzzz[j]);

                s_xxy_xxxx[j] = fr * s_xy_xxxx[j] + f2t * (s_y_xxxx[j] + 4.0 * s_xy_xxx[j]);

                t_xxy_xxxx[j] = fr * t_xy_xxxx[j] + f2t * (t_y_xxxx[j] + 4.0 * t_xy_xxx[j]) + f2z * (s_xxy_xxxx[j] - f2a * s_y_xxxx[j]);

                s_xxy_xxxy[j] = fr * s_xy_xxxy[j] + f2t * (s_y_xxxy[j] + 3.0 * s_xy_xxy[j]);

                t_xxy_xxxy[j] = fr * t_xy_xxxy[j] + f2t * (t_y_xxxy[j] + 3.0 * t_xy_xxy[j]) + f2z * (s_xxy_xxxy[j] - f2a * s_y_xxxy[j]);

                s_xxy_xxxz[j] = fr * s_xy_xxxz[j] + f2t * (s_y_xxxz[j] + 3.0 * s_xy_xxz[j]);

                t_xxy_xxxz[j] = fr * t_xy_xxxz[j] + f2t * (t_y_xxxz[j] + 3.0 * t_xy_xxz[j]) + f2z * (s_xxy_xxxz[j] - f2a * s_y_xxxz[j]);

                s_xxy_xxyy[j] = fr * s_xy_xxyy[j] + f2t * (s_y_xxyy[j] + 2.0 * s_xy_xyy[j]);

                t_xxy_xxyy[j] = fr * t_xy_xxyy[j] + f2t * (t_y_xxyy[j] + 2.0 * t_xy_xyy[j]) + f2z * (s_xxy_xxyy[j] - f2a * s_y_xxyy[j]);

                s_xxy_xxyz[j] = fr * s_xy_xxyz[j] + f2t * (s_y_xxyz[j] + 2.0 * s_xy_xyz[j]);

                t_xxy_xxyz[j] = fr * t_xy_xxyz[j] + f2t * (t_y_xxyz[j] + 2.0 * t_xy_xyz[j]) + f2z * (s_xxy_xxyz[j] - f2a * s_y_xxyz[j]);

                s_xxy_xxzz[j] = fr * s_xy_xxzz[j] + f2t * (s_y_xxzz[j] + 2.0 * s_xy_xzz[j]);

                t_xxy_xxzz[j] = fr * t_xy_xxzz[j] + f2t * (t_y_xxzz[j] + 2.0 * t_xy_xzz[j]) + f2z * (s_xxy_xxzz[j] - f2a * s_y_xxzz[j]);

                s_xxy_xyyy[j] = fr * s_xy_xyyy[j] + f2t * (s_y_xyyy[j] + s_xy_yyy[j]);

                t_xxy_xyyy[j] = fr * t_xy_xyyy[j] + f2t * (t_y_xyyy[j] + t_xy_yyy[j]) + f2z * (s_xxy_xyyy[j] - f2a * s_y_xyyy[j]);

                s_xxy_xyyz[j] = fr * s_xy_xyyz[j] + f2t * (s_y_xyyz[j] + s_xy_yyz[j]);

                t_xxy_xyyz[j] = fr * t_xy_xyyz[j] + f2t * (t_y_xyyz[j] + t_xy_yyz[j]) + f2z * (s_xxy_xyyz[j] - f2a * s_y_xyyz[j]);

                s_xxy_xyzz[j] = fr * s_xy_xyzz[j] + f2t * (s_y_xyzz[j] + s_xy_yzz[j]);

                t_xxy_xyzz[j] = fr * t_xy_xyzz[j] + f2t * (t_y_xyzz[j] + t_xy_yzz[j]) + f2z * (s_xxy_xyzz[j] - f2a * s_y_xyzz[j]);

                s_xxy_xzzz[j] = fr * s_xy_xzzz[j] + f2t * (s_y_xzzz[j] + s_xy_zzz[j]);

                t_xxy_xzzz[j] = fr * t_xy_xzzz[j] + f2t * (t_y_xzzz[j] + t_xy_zzz[j]) + f2z * (s_xxy_xzzz[j] - f2a * s_y_xzzz[j]);

                s_xxy_yyyy[j] = fr * s_xy_yyyy[j] + f2t * s_y_yyyy[j];

                t_xxy_yyyy[j] = fr * t_xy_yyyy[j] + f2t * t_y_yyyy[j] + f2z * (s_xxy_yyyy[j] - f2a * s_y_yyyy[j]);

                s_xxy_yyyz[j] = fr * s_xy_yyyz[j] + f2t * s_y_yyyz[j];

                t_xxy_yyyz[j] = fr * t_xy_yyyz[j] + f2t * t_y_yyyz[j] + f2z * (s_xxy_yyyz[j] - f2a * s_y_yyyz[j]);

                s_xxy_yyzz[j] = fr * s_xy_yyzz[j] + f2t * s_y_yyzz[j];

                t_xxy_yyzz[j] = fr * t_xy_yyzz[j] + f2t * t_y_yyzz[j] + f2z * (s_xxy_yyzz[j] - f2a * s_y_yyzz[j]);

                s_xxy_yzzz[j] = fr * s_xy_yzzz[j] + f2t * s_y_yzzz[j];

                t_xxy_yzzz[j] = fr * t_xy_yzzz[j] + f2t * t_y_yzzz[j] + f2z * (s_xxy_yzzz[j] - f2a * s_y_yzzz[j]);

                s_xxy_zzzz[j] = fr * s_xy_zzzz[j] + f2t * s_y_zzzz[j];

                t_xxy_zzzz[j] = fr * t_xy_zzzz[j] + f2t * t_y_zzzz[j] + f2z * (s_xxy_zzzz[j] - f2a * s_y_zzzz[j]);

                s_xxz_xxxx[j] = fr * s_xz_xxxx[j] + f2t * (s_z_xxxx[j] + 4.0 * s_xz_xxx[j]);

                t_xxz_xxxx[j] = fr * t_xz_xxxx[j] + f2t * (t_z_xxxx[j] + 4.0 * t_xz_xxx[j]) + f2z * (s_xxz_xxxx[j] - f2a * s_z_xxxx[j]);

                s_xxz_xxxy[j] = fr * s_xz_xxxy[j] + f2t * (s_z_xxxy[j] + 3.0 * s_xz_xxy[j]);

                t_xxz_xxxy[j] = fr * t_xz_xxxy[j] + f2t * (t_z_xxxy[j] + 3.0 * t_xz_xxy[j]) + f2z * (s_xxz_xxxy[j] - f2a * s_z_xxxy[j]);

                s_xxz_xxxz[j] = fr * s_xz_xxxz[j] + f2t * (s_z_xxxz[j] + 3.0 * s_xz_xxz[j]);

                t_xxz_xxxz[j] = fr * t_xz_xxxz[j] + f2t * (t_z_xxxz[j] + 3.0 * t_xz_xxz[j]) + f2z * (s_xxz_xxxz[j] - f2a * s_z_xxxz[j]);

                s_xxz_xxyy[j] = fr * s_xz_xxyy[j] + f2t * (s_z_xxyy[j] + 2.0 * s_xz_xyy[j]);

                t_xxz_xxyy[j] = fr * t_xz_xxyy[j] + f2t * (t_z_xxyy[j] + 2.0 * t_xz_xyy[j]) + f2z * (s_xxz_xxyy[j] - f2a * s_z_xxyy[j]);

                s_xxz_xxyz[j] = fr * s_xz_xxyz[j] + f2t * (s_z_xxyz[j] + 2.0 * s_xz_xyz[j]);

                t_xxz_xxyz[j] = fr * t_xz_xxyz[j] + f2t * (t_z_xxyz[j] + 2.0 * t_xz_xyz[j]) + f2z * (s_xxz_xxyz[j] - f2a * s_z_xxyz[j]);

                s_xxz_xxzz[j] = fr * s_xz_xxzz[j] + f2t * (s_z_xxzz[j] + 2.0 * s_xz_xzz[j]);

                t_xxz_xxzz[j] = fr * t_xz_xxzz[j] + f2t * (t_z_xxzz[j] + 2.0 * t_xz_xzz[j]) + f2z * (s_xxz_xxzz[j] - f2a * s_z_xxzz[j]);

                s_xxz_xyyy[j] = fr * s_xz_xyyy[j] + f2t * (s_z_xyyy[j] + s_xz_yyy[j]);

                t_xxz_xyyy[j] = fr * t_xz_xyyy[j] + f2t * (t_z_xyyy[j] + t_xz_yyy[j]) + f2z * (s_xxz_xyyy[j] - f2a * s_z_xyyy[j]);

                s_xxz_xyyz[j] = fr * s_xz_xyyz[j] + f2t * (s_z_xyyz[j] + s_xz_yyz[j]);

                t_xxz_xyyz[j] = fr * t_xz_xyyz[j] + f2t * (t_z_xyyz[j] + t_xz_yyz[j]) + f2z * (s_xxz_xyyz[j] - f2a * s_z_xyyz[j]);

                s_xxz_xyzz[j] = fr * s_xz_xyzz[j] + f2t * (s_z_xyzz[j] + s_xz_yzz[j]);

                t_xxz_xyzz[j] = fr * t_xz_xyzz[j] + f2t * (t_z_xyzz[j] + t_xz_yzz[j]) + f2z * (s_xxz_xyzz[j] - f2a * s_z_xyzz[j]);

                s_xxz_xzzz[j] = fr * s_xz_xzzz[j] + f2t * (s_z_xzzz[j] + s_xz_zzz[j]);

                t_xxz_xzzz[j] = fr * t_xz_xzzz[j] + f2t * (t_z_xzzz[j] + t_xz_zzz[j]) + f2z * (s_xxz_xzzz[j] - f2a * s_z_xzzz[j]);

                s_xxz_yyyy[j] = fr * s_xz_yyyy[j] + f2t * s_z_yyyy[j];

                t_xxz_yyyy[j] = fr * t_xz_yyyy[j] + f2t * t_z_yyyy[j] + f2z * (s_xxz_yyyy[j] - f2a * s_z_yyyy[j]);

                s_xxz_yyyz[j] = fr * s_xz_yyyz[j] + f2t * s_z_yyyz[j];

                t_xxz_yyyz[j] = fr * t_xz_yyyz[j] + f2t * t_z_yyyz[j] + f2z * (s_xxz_yyyz[j] - f2a * s_z_yyyz[j]);

                s_xxz_yyzz[j] = fr * s_xz_yyzz[j] + f2t * s_z_yyzz[j];

                t_xxz_yyzz[j] = fr * t_xz_yyzz[j] + f2t * t_z_yyzz[j] + f2z * (s_xxz_yyzz[j] - f2a * s_z_yyzz[j]);

                s_xxz_yzzz[j] = fr * s_xz_yzzz[j] + f2t * s_z_yzzz[j];

                t_xxz_yzzz[j] = fr * t_xz_yzzz[j] + f2t * t_z_yzzz[j] + f2z * (s_xxz_yzzz[j] - f2a * s_z_yzzz[j]);

                s_xxz_zzzz[j] = fr * s_xz_zzzz[j] + f2t * s_z_zzzz[j];

                t_xxz_zzzz[j] = fr * t_xz_zzzz[j] + f2t * t_z_zzzz[j] + f2z * (s_xxz_zzzz[j] - f2a * s_z_zzzz[j]);

                s_xyy_xxxx[j] = fr * s_yy_xxxx[j] + f2t * 4.0 * s_yy_xxx[j];

                t_xyy_xxxx[j] = fr * t_yy_xxxx[j] + f2t * 4.0 * t_yy_xxx[j] + f2z * s_xyy_xxxx[j];

                s_xyy_xxxy[j] = fr * s_yy_xxxy[j] + f2t * 3.0 * s_yy_xxy[j];

                t_xyy_xxxy[j] = fr * t_yy_xxxy[j] + f2t * 3.0 * t_yy_xxy[j] + f2z * s_xyy_xxxy[j];

                s_xyy_xxxz[j] = fr * s_yy_xxxz[j] + f2t * 3.0 * s_yy_xxz[j];

                t_xyy_xxxz[j] = fr * t_yy_xxxz[j] + f2t * 3.0 * t_yy_xxz[j] + f2z * s_xyy_xxxz[j];

                s_xyy_xxyy[j] = fr * s_yy_xxyy[j] + f2t * 2.0 * s_yy_xyy[j];

                t_xyy_xxyy[j] = fr * t_yy_xxyy[j] + f2t * 2.0 * t_yy_xyy[j] + f2z * s_xyy_xxyy[j];

                s_xyy_xxyz[j] = fr * s_yy_xxyz[j] + f2t * 2.0 * s_yy_xyz[j];

                t_xyy_xxyz[j] = fr * t_yy_xxyz[j] + f2t * 2.0 * t_yy_xyz[j] + f2z * s_xyy_xxyz[j];

                s_xyy_xxzz[j] = fr * s_yy_xxzz[j] + f2t * 2.0 * s_yy_xzz[j];

                t_xyy_xxzz[j] = fr * t_yy_xxzz[j] + f2t * 2.0 * t_yy_xzz[j] + f2z * s_xyy_xxzz[j];

                s_xyy_xyyy[j] = fr * s_yy_xyyy[j] + f2t * s_yy_yyy[j];

                t_xyy_xyyy[j] = fr * t_yy_xyyy[j] + f2t * t_yy_yyy[j] + f2z * s_xyy_xyyy[j];

                s_xyy_xyyz[j] = fr * s_yy_xyyz[j] + f2t * s_yy_yyz[j];

                t_xyy_xyyz[j] = fr * t_yy_xyyz[j] + f2t * t_yy_yyz[j] + f2z * s_xyy_xyyz[j];

                s_xyy_xyzz[j] = fr * s_yy_xyzz[j] + f2t * s_yy_yzz[j];

                t_xyy_xyzz[j] = fr * t_yy_xyzz[j] + f2t * t_yy_yzz[j] + f2z * s_xyy_xyzz[j];

                s_xyy_xzzz[j] = fr * s_yy_xzzz[j] + f2t * s_yy_zzz[j];

                t_xyy_xzzz[j] = fr * t_yy_xzzz[j] + f2t * t_yy_zzz[j] + f2z * s_xyy_xzzz[j];

                s_xyy_yyyy[j] = fr * s_yy_yyyy[j];

                t_xyy_yyyy[j] = fr * t_yy_yyyy[j] + f2z * s_xyy_yyyy[j];

                s_xyy_yyyz[j] = fr * s_yy_yyyz[j];

                t_xyy_yyyz[j] = fr * t_yy_yyyz[j] + f2z * s_xyy_yyyz[j];

                s_xyy_yyzz[j] = fr * s_yy_yyzz[j];

                t_xyy_yyzz[j] = fr * t_yy_yyzz[j] + f2z * s_xyy_yyzz[j];

                s_xyy_yzzz[j] = fr * s_yy_yzzz[j];

                t_xyy_yzzz[j] = fr * t_yy_yzzz[j] + f2z * s_xyy_yzzz[j];

                s_xyy_zzzz[j] = fr * s_yy_zzzz[j];

                t_xyy_zzzz[j] = fr * t_yy_zzzz[j] + f2z * s_xyy_zzzz[j];

                s_xyz_xxxx[j] = fr * s_yz_xxxx[j] + f2t * 4.0 * s_yz_xxx[j];

                t_xyz_xxxx[j] = fr * t_yz_xxxx[j] + f2t * 4.0 * t_yz_xxx[j] + f2z * s_xyz_xxxx[j];

                s_xyz_xxxy[j] = fr * s_yz_xxxy[j] + f2t * 3.0 * s_yz_xxy[j];

                t_xyz_xxxy[j] = fr * t_yz_xxxy[j] + f2t * 3.0 * t_yz_xxy[j] + f2z * s_xyz_xxxy[j];

                s_xyz_xxxz[j] = fr * s_yz_xxxz[j] + f2t * 3.0 * s_yz_xxz[j];

                t_xyz_xxxz[j] = fr * t_yz_xxxz[j] + f2t * 3.0 * t_yz_xxz[j] + f2z * s_xyz_xxxz[j];

                s_xyz_xxyy[j] = fr * s_yz_xxyy[j] + f2t * 2.0 * s_yz_xyy[j];

                t_xyz_xxyy[j] = fr * t_yz_xxyy[j] + f2t * 2.0 * t_yz_xyy[j] + f2z * s_xyz_xxyy[j];

                s_xyz_xxyz[j] = fr * s_yz_xxyz[j] + f2t * 2.0 * s_yz_xyz[j];

                t_xyz_xxyz[j] = fr * t_yz_xxyz[j] + f2t * 2.0 * t_yz_xyz[j] + f2z * s_xyz_xxyz[j];

                s_xyz_xxzz[j] = fr * s_yz_xxzz[j] + f2t * 2.0 * s_yz_xzz[j];

                t_xyz_xxzz[j] = fr * t_yz_xxzz[j] + f2t * 2.0 * t_yz_xzz[j] + f2z * s_xyz_xxzz[j];

                s_xyz_xyyy[j] = fr * s_yz_xyyy[j] + f2t * s_yz_yyy[j];

                t_xyz_xyyy[j] = fr * t_yz_xyyy[j] + f2t * t_yz_yyy[j] + f2z * s_xyz_xyyy[j];

                s_xyz_xyyz[j] = fr * s_yz_xyyz[j] + f2t * s_yz_yyz[j];

                t_xyz_xyyz[j] = fr * t_yz_xyyz[j] + f2t * t_yz_yyz[j] + f2z * s_xyz_xyyz[j];

                s_xyz_xyzz[j] = fr * s_yz_xyzz[j] + f2t * s_yz_yzz[j];

                t_xyz_xyzz[j] = fr * t_yz_xyzz[j] + f2t * t_yz_yzz[j] + f2z * s_xyz_xyzz[j];

                s_xyz_xzzz[j] = fr * s_yz_xzzz[j] + f2t * s_yz_zzz[j];

                t_xyz_xzzz[j] = fr * t_yz_xzzz[j] + f2t * t_yz_zzz[j] + f2z * s_xyz_xzzz[j];

                s_xyz_yyyy[j] = fr * s_yz_yyyy[j];

                t_xyz_yyyy[j] = fr * t_yz_yyyy[j] + f2z * s_xyz_yyyy[j];

                s_xyz_yyyz[j] = fr * s_yz_yyyz[j];

                t_xyz_yyyz[j] = fr * t_yz_yyyz[j] + f2z * s_xyz_yyyz[j];

                s_xyz_yyzz[j] = fr * s_yz_yyzz[j];

                t_xyz_yyzz[j] = fr * t_yz_yyzz[j] + f2z * s_xyz_yyzz[j];

                s_xyz_yzzz[j] = fr * s_yz_yzzz[j];

                t_xyz_yzzz[j] = fr * t_yz_yzzz[j] + f2z * s_xyz_yzzz[j];

                s_xyz_zzzz[j] = fr * s_yz_zzzz[j];

                t_xyz_zzzz[j] = fr * t_yz_zzzz[j] + f2z * s_xyz_zzzz[j];

                s_xzz_xxxx[j] = fr * s_zz_xxxx[j] + f2t * 4.0 * s_zz_xxx[j];

                t_xzz_xxxx[j] = fr * t_zz_xxxx[j] + f2t * 4.0 * t_zz_xxx[j] + f2z * s_xzz_xxxx[j];

                s_xzz_xxxy[j] = fr * s_zz_xxxy[j] + f2t * 3.0 * s_zz_xxy[j];

                t_xzz_xxxy[j] = fr * t_zz_xxxy[j] + f2t * 3.0 * t_zz_xxy[j] + f2z * s_xzz_xxxy[j];

                s_xzz_xxxz[j] = fr * s_zz_xxxz[j] + f2t * 3.0 * s_zz_xxz[j];

                t_xzz_xxxz[j] = fr * t_zz_xxxz[j] + f2t * 3.0 * t_zz_xxz[j] + f2z * s_xzz_xxxz[j];

                s_xzz_xxyy[j] = fr * s_zz_xxyy[j] + f2t * 2.0 * s_zz_xyy[j];

                t_xzz_xxyy[j] = fr * t_zz_xxyy[j] + f2t * 2.0 * t_zz_xyy[j] + f2z * s_xzz_xxyy[j];

                s_xzz_xxyz[j] = fr * s_zz_xxyz[j] + f2t * 2.0 * s_zz_xyz[j];

                t_xzz_xxyz[j] = fr * t_zz_xxyz[j] + f2t * 2.0 * t_zz_xyz[j] + f2z * s_xzz_xxyz[j];

                s_xzz_xxzz[j] = fr * s_zz_xxzz[j] + f2t * 2.0 * s_zz_xzz[j];

                t_xzz_xxzz[j] = fr * t_zz_xxzz[j] + f2t * 2.0 * t_zz_xzz[j] + f2z * s_xzz_xxzz[j];

                s_xzz_xyyy[j] = fr * s_zz_xyyy[j] + f2t * s_zz_yyy[j];

                t_xzz_xyyy[j] = fr * t_zz_xyyy[j] + f2t * t_zz_yyy[j] + f2z * s_xzz_xyyy[j];

                s_xzz_xyyz[j] = fr * s_zz_xyyz[j] + f2t * s_zz_yyz[j];

                t_xzz_xyyz[j] = fr * t_zz_xyyz[j] + f2t * t_zz_yyz[j] + f2z * s_xzz_xyyz[j];

                s_xzz_xyzz[j] = fr * s_zz_xyzz[j] + f2t * s_zz_yzz[j];

                t_xzz_xyzz[j] = fr * t_zz_xyzz[j] + f2t * t_zz_yzz[j] + f2z * s_xzz_xyzz[j];

                s_xzz_xzzz[j] = fr * s_zz_xzzz[j] + f2t * s_zz_zzz[j];

                t_xzz_xzzz[j] = fr * t_zz_xzzz[j] + f2t * t_zz_zzz[j] + f2z * s_xzz_xzzz[j];

                s_xzz_yyyy[j] = fr * s_zz_yyyy[j];

                t_xzz_yyyy[j] = fr * t_zz_yyyy[j] + f2z * s_xzz_yyyy[j];

                s_xzz_yyyz[j] = fr * s_zz_yyyz[j];

                t_xzz_yyyz[j] = fr * t_zz_yyyz[j] + f2z * s_xzz_yyyz[j];

                s_xzz_yyzz[j] = fr * s_zz_yyzz[j];

                t_xzz_yyzz[j] = fr * t_zz_yyzz[j] + f2z * s_xzz_yyzz[j];

                s_xzz_yzzz[j] = fr * s_zz_yzzz[j];

                t_xzz_yzzz[j] = fr * t_zz_yzzz[j] + f2z * s_xzz_yzzz[j];

                s_xzz_zzzz[j] = fr * s_zz_zzzz[j];

                t_xzz_zzzz[j] = fr * t_zz_zzzz[j] + f2z * s_xzz_zzzz[j];

                // leading y component

                fr = pay[j];

                s_yyy_xxxx[j] = fr * s_yy_xxxx[j] + f2t * 2.0 * s_y_xxxx[j];

                t_yyy_xxxx[j] = fr * t_yy_xxxx[j] + f2t * 2.0 * t_y_xxxx[j] + f2z * (s_yyy_xxxx[j] - f2a * 2.0 * s_y_xxxx[j]);

                s_yyy_xxxy[j] = fr * s_yy_xxxy[j] + f2t * (2.0 * s_y_xxxy[j] + s_yy_xxx[j]);

                t_yyy_xxxy[j] = fr * t_yy_xxxy[j] + f2t * (2.0 * t_y_xxxy[j] + t_yy_xxx[j]) + f2z * (s_yyy_xxxy[j] - f2a * 2.0 * s_y_xxxy[j]);

                s_yyy_xxxz[j] = fr * s_yy_xxxz[j] + f2t * 2.0 * s_y_xxxz[j];

                t_yyy_xxxz[j] = fr * t_yy_xxxz[j] + f2t * 2.0 * t_y_xxxz[j] + f2z * (s_yyy_xxxz[j] - f2a * 2.0 * s_y_xxxz[j]);

                s_yyy_xxyy[j] = fr * s_yy_xxyy[j] + f2t * (2.0 * s_y_xxyy[j] + 2.0 * s_yy_xxy[j]);

                t_yyy_xxyy[j] = fr * t_yy_xxyy[j] + f2t * (2.0 * t_y_xxyy[j] + 2.0 * t_yy_xxy[j]) + f2z * (s_yyy_xxyy[j] - f2a * 2.0 * s_y_xxyy[j]);

                s_yyy_xxyz[j] = fr * s_yy_xxyz[j] + f2t * (2.0 * s_y_xxyz[j] + s_yy_xxz[j]);

                t_yyy_xxyz[j] = fr * t_yy_xxyz[j] + f2t * (2.0 * t_y_xxyz[j] + t_yy_xxz[j]) + f2z * (s_yyy_xxyz[j] - f2a * 2.0 * s_y_xxyz[j]);

                s_yyy_xxzz[j] = fr * s_yy_xxzz[j] + f2t * 2.0 * s_y_xxzz[j];

                t_yyy_xxzz[j] = fr * t_yy_xxzz[j] + f2t * 2.0 * t_y_xxzz[j] + f2z * (s_yyy_xxzz[j] - f2a * 2.0 * s_y_xxzz[j]);

                s_yyy_xyyy[j] = fr * s_yy_xyyy[j] + f2t * (2.0 * s_y_xyyy[j] + 3.0 * s_yy_xyy[j]);

                t_yyy_xyyy[j] = fr * t_yy_xyyy[j] + f2t * (2.0 * t_y_xyyy[j] + 3.0 * t_yy_xyy[j]) + f2z * (s_yyy_xyyy[j] - f2a * 2.0 * s_y_xyyy[j]);

                s_yyy_xyyz[j] = fr * s_yy_xyyz[j] + f2t * (2.0 * s_y_xyyz[j] + 2.0 * s_yy_xyz[j]);

                t_yyy_xyyz[j] = fr * t_yy_xyyz[j] + f2t * (2.0 * t_y_xyyz[j] + 2.0 * t_yy_xyz[j]) + f2z * (s_yyy_xyyz[j] - f2a * 2.0 * s_y_xyyz[j]);

                s_yyy_xyzz[j] = fr * s_yy_xyzz[j] + f2t * (2.0 * s_y_xyzz[j] + s_yy_xzz[j]);

                t_yyy_xyzz[j] = fr * t_yy_xyzz[j] + f2t * (2.0 * t_y_xyzz[j] + t_yy_xzz[j]) + f2z * (s_yyy_xyzz[j] - f2a * 2.0 * s_y_xyzz[j]);

                s_yyy_xzzz[j] = fr * s_yy_xzzz[j] + f2t * 2.0 * s_y_xzzz[j];

                t_yyy_xzzz[j] = fr * t_yy_xzzz[j] + f2t * 2.0 * t_y_xzzz[j] + f2z * (s_yyy_xzzz[j] - f2a * 2.0 * s_y_xzzz[j]);

                s_yyy_yyyy[j] = fr * s_yy_yyyy[j] + f2t * (2.0 * s_y_yyyy[j] + 4.0 * s_yy_yyy[j]);

                t_yyy_yyyy[j] = fr * t_yy_yyyy[j] + f2t * (2.0 * t_y_yyyy[j] + 4.0 * t_yy_yyy[j]) + f2z * (s_yyy_yyyy[j] - f2a * 2.0 * s_y_yyyy[j]);

                s_yyy_yyyz[j] = fr * s_yy_yyyz[j] + f2t * (2.0 * s_y_yyyz[j] + 3.0 * s_yy_yyz[j]);

                t_yyy_yyyz[j] = fr * t_yy_yyyz[j] + f2t * (2.0 * t_y_yyyz[j] + 3.0 * t_yy_yyz[j]) + f2z * (s_yyy_yyyz[j] - f2a * 2.0 * s_y_yyyz[j]);

                s_yyy_yyzz[j] = fr * s_yy_yyzz[j] + f2t * (2.0 * s_y_yyzz[j] + 2.0 * s_yy_yzz[j]);

                t_yyy_yyzz[j] = fr * t_yy_yyzz[j] + f2t * (2.0 * t_y_yyzz[j] + 2.0 * t_yy_yzz[j]) + f2z * (s_yyy_yyzz[j] - f2a * 2.0 * s_y_yyzz[j]);

                s_yyy_yzzz[j] = fr * s_yy_yzzz[j] + f2t * (2.0 * s_y_yzzz[j] + s_yy_zzz[j]);

                t_yyy_yzzz[j] = fr * t_yy_yzzz[j] + f2t * (2.0 * t_y_yzzz[j] + t_yy_zzz[j]) + f2z * (s_yyy_yzzz[j] - f2a * 2.0 * s_y_yzzz[j]);

                s_yyy_zzzz[j] = fr * s_yy_zzzz[j] + f2t * 2.0 * s_y_zzzz[j];

                t_yyy_zzzz[j] = fr * t_yy_zzzz[j] + f2t * 2.0 * t_y_zzzz[j] + f2z * (s_yyy_zzzz[j] - f2a * 2.0 * s_y_zzzz[j]);

                s_yyz_xxxx[j] = fr * s_yz_xxxx[j] + f2t * s_z_xxxx[j];

                t_yyz_xxxx[j] = fr * t_yz_xxxx[j] + f2t * t_z_xxxx[j] + f2z * (s_yyz_xxxx[j] - f2a * s_z_xxxx[j]);

                s_yyz_xxxy[j] = fr * s_yz_xxxy[j] + f2t * (s_z_xxxy[j] + s_yz_xxx[j]);

                t_yyz_xxxy[j] = fr * t_yz_xxxy[j] + f2t * (t_z_xxxy[j] + t_yz_xxx[j]) + f2z * (s_yyz_xxxy[j] - f2a * s_z_xxxy[j]);

                s_yyz_xxxz[j] = fr * s_yz_xxxz[j] + f2t * s_z_xxxz[j];

                t_yyz_xxxz[j] = fr * t_yz_xxxz[j] + f2t * t_z_xxxz[j] + f2z * (s_yyz_xxxz[j] - f2a * s_z_xxxz[j]);

                s_yyz_xxyy[j] = fr * s_yz_xxyy[j] + f2t * (s_z_xxyy[j] + 2.0 * s_yz_xxy[j]);

                t_yyz_xxyy[j] = fr * t_yz_xxyy[j] + f2t * (t_z_xxyy[j] + 2.0 * t_yz_xxy[j]) + f2z * (s_yyz_xxyy[j] - f2a * s_z_xxyy[j]);

                s_yyz_xxyz[j] = fr * s_yz_xxyz[j] + f2t * (s_z_xxyz[j] + s_yz_xxz[j]);

                t_yyz_xxyz[j] = fr * t_yz_xxyz[j] + f2t * (t_z_xxyz[j] + t_yz_xxz[j]) + f2z * (s_yyz_xxyz[j] - f2a * s_z_xxyz[j]);

                s_yyz_xxzz[j] = fr * s_yz_xxzz[j] + f2t * s_z_xxzz[j];

                t_yyz_xxzz[j] = fr * t_yz_xxzz[j] + f2t * t_z_xxzz[j] + f2z * (s_yyz_xxzz[j] - f2a * s_z_xxzz[j]);

                s_yyz_xyyy[j] = fr * s_yz_xyyy[j] + f2t * (s_z_xyyy[j] + 3.0 * s_yz_xyy[j]);

                t_yyz_xyyy[j] = fr * t_yz_xyyy[j] + f2t * (t_z_xyyy[j] + 3.0 * t_yz_xyy[j]) + f2z * (s_yyz_xyyy[j] - f2a * s_z_xyyy[j]);

                s_yyz_xyyz[j] = fr * s_yz_xyyz[j] + f2t * (s_z_xyyz[j] + 2.0 * s_yz_xyz[j]);

                t_yyz_xyyz[j] = fr * t_yz_xyyz[j] + f2t * (t_z_xyyz[j] + 2.0 * t_yz_xyz[j]) + f2z * (s_yyz_xyyz[j] - f2a * s_z_xyyz[j]);

                s_yyz_xyzz[j] = fr * s_yz_xyzz[j] + f2t * (s_z_xyzz[j] + s_yz_xzz[j]);

                t_yyz_xyzz[j] = fr * t_yz_xyzz[j] + f2t * (t_z_xyzz[j] + t_yz_xzz[j]) + f2z * (s_yyz_xyzz[j] - f2a * s_z_xyzz[j]);

                s_yyz_xzzz[j] = fr * s_yz_xzzz[j] + f2t * s_z_xzzz[j];

                t_yyz_xzzz[j] = fr * t_yz_xzzz[j] + f2t * t_z_xzzz[j] + f2z * (s_yyz_xzzz[j] - f2a * s_z_xzzz[j]);

                s_yyz_yyyy[j] = fr * s_yz_yyyy[j] + f2t * (s_z_yyyy[j] + 4.0 * s_yz_yyy[j]);

                t_yyz_yyyy[j] = fr * t_yz_yyyy[j] + f2t * (t_z_yyyy[j] + 4.0 * t_yz_yyy[j]) + f2z * (s_yyz_yyyy[j] - f2a * s_z_yyyy[j]);

                s_yyz_yyyz[j] = fr * s_yz_yyyz[j] + f2t * (s_z_yyyz[j] + 3.0 * s_yz_yyz[j]);

                t_yyz_yyyz[j] = fr * t_yz_yyyz[j] + f2t * (t_z_yyyz[j] + 3.0 * t_yz_yyz[j]) + f2z * (s_yyz_yyyz[j] - f2a * s_z_yyyz[j]);

                s_yyz_yyzz[j] = fr * s_yz_yyzz[j] + f2t * (s_z_yyzz[j] + 2.0 * s_yz_yzz[j]);

                t_yyz_yyzz[j] = fr * t_yz_yyzz[j] + f2t * (t_z_yyzz[j] + 2.0 * t_yz_yzz[j]) + f2z * (s_yyz_yyzz[j] - f2a * s_z_yyzz[j]);

                s_yyz_yzzz[j] = fr * s_yz_yzzz[j] + f2t * (s_z_yzzz[j] + s_yz_zzz[j]);

                t_yyz_yzzz[j] = fr * t_yz_yzzz[j] + f2t * (t_z_yzzz[j] + t_yz_zzz[j]) + f2z * (s_yyz_yzzz[j] - f2a * s_z_yzzz[j]);

                s_yyz_zzzz[j] = fr * s_yz_zzzz[j] + f2t * s_z_zzzz[j];

                t_yyz_zzzz[j] = fr * t_yz_zzzz[j] + f2t * t_z_zzzz[j] + f2z * (s_yyz_zzzz[j] - f2a * s_z_zzzz[j]);

                s_yzz_xxxx[j] = fr * s_zz_xxxx[j];

                t_yzz_xxxx[j] = fr * t_zz_xxxx[j] + f2z * s_yzz_xxxx[j];

                s_yzz_xxxy[j] = fr * s_zz_xxxy[j] + f2t * s_zz_xxx[j];

                t_yzz_xxxy[j] = fr * t_zz_xxxy[j] + f2t * t_zz_xxx[j] + f2z * s_yzz_xxxy[j];

                s_yzz_xxxz[j] = fr * s_zz_xxxz[j];

                t_yzz_xxxz[j] = fr * t_zz_xxxz[j] + f2z * s_yzz_xxxz[j];

                s_yzz_xxyy[j] = fr * s_zz_xxyy[j] + f2t * 2.0 * s_zz_xxy[j];

                t_yzz_xxyy[j] = fr * t_zz_xxyy[j] + f2t * 2.0 * t_zz_xxy[j] + f2z * s_yzz_xxyy[j];

                s_yzz_xxyz[j] = fr * s_zz_xxyz[j] + f2t * s_zz_xxz[j];

                t_yzz_xxyz[j] = fr * t_zz_xxyz[j] + f2t * t_zz_xxz[j] + f2z * s_yzz_xxyz[j];

                s_yzz_xxzz[j] = fr * s_zz_xxzz[j];

                t_yzz_xxzz[j] = fr * t_zz_xxzz[j] + f2z * s_yzz_xxzz[j];

                s_yzz_xyyy[j] = fr * s_zz_xyyy[j] + f2t * 3.0 * s_zz_xyy[j];

                t_yzz_xyyy[j] = fr * t_zz_xyyy[j] + f2t * 3.0 * t_zz_xyy[j] + f2z * s_yzz_xyyy[j];

                s_yzz_xyyz[j] = fr * s_zz_xyyz[j] + f2t * 2.0 * s_zz_xyz[j];

                t_yzz_xyyz[j] = fr * t_zz_xyyz[j] + f2t * 2.0 * t_zz_xyz[j] + f2z * s_yzz_xyyz[j];

                s_yzz_xyzz[j] = fr * s_zz_xyzz[j] + f2t * s_zz_xzz[j];

                t_yzz_xyzz[j] = fr * t_zz_xyzz[j] + f2t * t_zz_xzz[j] + f2z * s_yzz_xyzz[j];

                s_yzz_xzzz[j] = fr * s_zz_xzzz[j];

                t_yzz_xzzz[j] = fr * t_zz_xzzz[j] + f2z * s_yzz_xzzz[j];

                s_yzz_yyyy[j] = fr * s_zz_yyyy[j] + f2t * 4.0 * s_zz_yyy[j];

                t_yzz_yyyy[j] = fr * t_zz_yyyy[j] + f2t * 4.0 * t_zz_yyy[j] + f2z * s_yzz_yyyy[j];

                s_yzz_yyyz[j] = fr * s_zz_yyyz[j] + f2t * 3.0 * s_zz_yyz[j];

                t_yzz_yyyz[j] = fr * t_zz_yyyz[j] + f2t * 3.0 * t_zz_yyz[j] + f2z * s_yzz_yyyz[j];

                s_yzz_yyzz[j] = fr * s_zz_yyzz[j] + f2t * 2.0 * s_zz_yzz[j];

                t_yzz_yyzz[j] = fr * t_zz_yyzz[j] + f2t * 2.0 * t_zz_yzz[j] + f2z * s_yzz_yyzz[j];

                s_yzz_yzzz[j] = fr * s_zz_yzzz[j] + f2t * s_zz_zzz[j];

                t_yzz_yzzz[j] = fr * t_zz_yzzz[j] + f2t * t_zz_zzz[j] + f2z * s_yzz_yzzz[j];

                s_yzz_zzzz[j] = fr * s_zz_zzzz[j];

                t_yzz_zzzz[j] = fr * t_zz_zzzz[j] + f2z * s_yzz_zzzz[j];

                // leading z component

                fr = paz[j];

                s_zzz_xxxx[j] = fr * s_zz_xxxx[j] + f2t * 2.0 * s_z_xxxx[j];

                t_zzz_xxxx[j] = fr * t_zz_xxxx[j] + f2t * 2.0 * t_z_xxxx[j] + f2z * (s_zzz_xxxx[j] - f2a * 2.0 * s_z_xxxx[j]);

                s_zzz_xxxy[j] = fr * s_zz_xxxy[j] + f2t * 2.0 * s_z_xxxy[j];

                t_zzz_xxxy[j] = fr * t_zz_xxxy[j] + f2t * 2.0 * t_z_xxxy[j] + f2z * (s_zzz_xxxy[j] - f2a * 2.0 * s_z_xxxy[j]);

                s_zzz_xxxz[j] = fr * s_zz_xxxz[j] + f2t * (2.0 * s_z_xxxz[j] + s_zz_xxx[j]);

                t_zzz_xxxz[j] = fr * t_zz_xxxz[j] + f2t * (2.0 * t_z_xxxz[j] + t_zz_xxx[j]) + f2z * (s_zzz_xxxz[j] - f2a * 2.0 * s_z_xxxz[j]);

                s_zzz_xxyy[j] = fr * s_zz_xxyy[j] + f2t * 2.0 * s_z_xxyy[j];

                t_zzz_xxyy[j] = fr * t_zz_xxyy[j] + f2t * 2.0 * t_z_xxyy[j] + f2z * (s_zzz_xxyy[j] - f2a * 2.0 * s_z_xxyy[j]);

                s_zzz_xxyz[j] = fr * s_zz_xxyz[j] + f2t * (2.0 * s_z_xxyz[j] + s_zz_xxy[j]);

                t_zzz_xxyz[j] = fr * t_zz_xxyz[j] + f2t * (2.0 * t_z_xxyz[j] + t_zz_xxy[j]) + f2z * (s_zzz_xxyz[j] - f2a * 2.0 * s_z_xxyz[j]);

                s_zzz_xxzz[j] = fr * s_zz_xxzz[j] + f2t * (2.0 * s_z_xxzz[j] + 2.0 * s_zz_xxz[j]);

                t_zzz_xxzz[j] = fr * t_zz_xxzz[j] + f2t * (2.0 * t_z_xxzz[j] + 2.0 * t_zz_xxz[j]) + f2z * (s_zzz_xxzz[j] - f2a * 2.0 * s_z_xxzz[j]);

                s_zzz_xyyy[j] = fr * s_zz_xyyy[j] + f2t * 2.0 * s_z_xyyy[j];

                t_zzz_xyyy[j] = fr * t_zz_xyyy[j] + f2t * 2.0 * t_z_xyyy[j] + f2z * (s_zzz_xyyy[j] - f2a * 2.0 * s_z_xyyy[j]);

                s_zzz_xyyz[j] = fr * s_zz_xyyz[j] + f2t * (2.0 * s_z_xyyz[j] + s_zz_xyy[j]);

                t_zzz_xyyz[j] = fr * t_zz_xyyz[j] + f2t * (2.0 * t_z_xyyz[j] + t_zz_xyy[j]) + f2z * (s_zzz_xyyz[j] - f2a * 2.0 * s_z_xyyz[j]);

                s_zzz_xyzz[j] = fr * s_zz_xyzz[j] + f2t * (2.0 * s_z_xyzz[j] + 2.0 * s_zz_xyz[j]);

                t_zzz_xyzz[j] = fr * t_zz_xyzz[j] + f2t * (2.0 * t_z_xyzz[j] + 2.0 * t_zz_xyz[j]) + f2z * (s_zzz_xyzz[j] - f2a * 2.0 * s_z_xyzz[j]);

                s_zzz_xzzz[j] = fr * s_zz_xzzz[j] + f2t * (2.0 * s_z_xzzz[j] + 3.0 * s_zz_xzz[j]);

                t_zzz_xzzz[j] = fr * t_zz_xzzz[j] + f2t * (2.0 * t_z_xzzz[j] + 3.0 * t_zz_xzz[j]) + f2z * (s_zzz_xzzz[j] - f2a * 2.0 * s_z_xzzz[j]);

                s_zzz_yyyy[j] = fr * s_zz_yyyy[j] + f2t * 2.0 * s_z_yyyy[j];

                t_zzz_yyyy[j] = fr * t_zz_yyyy[j] + f2t * 2.0 * t_z_yyyy[j] + f2z * (s_zzz_yyyy[j] - f2a * 2.0 * s_z_yyyy[j]);

                s_zzz_yyyz[j] = fr * s_zz_yyyz[j] + f2t * (2.0 * s_z_yyyz[j] + s_zz_yyy[j]);

                t_zzz_yyyz[j] = fr * t_zz_yyyz[j] + f2t * (2.0 * t_z_yyyz[j] + t_zz_yyy[j]) + f2z * (s_zzz_yyyz[j] - f2a * 2.0 * s_z_yyyz[j]);

                s_zzz_yyzz[j] = fr * s_zz_yyzz[j] + f2t * (2.0 * s_z_yyzz[j] + 2.0 * s_zz_yyz[j]);

                t_zzz_yyzz[j] = fr * t_zz_yyzz[j] + f2t * (2.0 * t_z_yyzz[j] + 2.0 * t_zz_yyz[j]) + f2z * (s_zzz_yyzz[j] - f2a * 2.0 * s_z_yyzz[j]);

                s_zzz_yzzz[j] = fr * s_zz_yzzz[j] + f2t * (2.0 * s_z_yzzz[j] + 3.0 * s_zz_yzz[j]);

                t_zzz_yzzz[j] = fr * t_zz_yzzz[j] + f2t * (2.0 * t_z_yzzz[j] + 3.0 * t_zz_yzz[j]) + f2z * (s_zzz_yzzz[j] - f2a * 2.0 * s_z_yzzz[j]);

                s_zzz_zzzz[j] = fr * s_zz_zzzz[j] + f2t * (2.0 * s_z_zzzz[j] + 4.0 * s_zz_zzz[j]);

                t_zzz_zzzz[j] = fr * t_zz_zzzz[j] + f2t * (2.0 * t_z_zzzz[j] + 4.0 * t_zz_zzz[j]) + f2z * (s_zzz_zzzz[j] - f2a * 2.0 * s_z_zzzz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGF(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 3, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 3, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|D) integrals

            auto s_xxx_xx = primBuffer.data(skoff + 60 * idx);

            auto s_xxx_xy = primBuffer.data(skoff + 60 * idx + 1);

            auto s_xxx_xz = primBuffer.data(skoff + 60 * idx + 2);

            auto s_xxx_yy = primBuffer.data(skoff + 60 * idx + 3);

            auto s_xxx_yz = primBuffer.data(skoff + 60 * idx + 4);

            auto s_xxx_zz = primBuffer.data(skoff + 60 * idx + 5);

            auto s_xxy_xx = primBuffer.data(skoff + 60 * idx + 6);

            auto s_xxy_xy = primBuffer.data(skoff + 60 * idx + 7);

            auto s_xxy_xz = primBuffer.data(skoff + 60 * idx + 8);

            auto s_xxy_yy = primBuffer.data(skoff + 60 * idx + 9);

            auto s_xxy_yz = primBuffer.data(skoff + 60 * idx + 10);

            auto s_xxy_zz = primBuffer.data(skoff + 60 * idx + 11);

            auto s_xxz_xx = primBuffer.data(skoff + 60 * idx + 12);

            auto s_xxz_xy = primBuffer.data(skoff + 60 * idx + 13);

            auto s_xxz_xz = primBuffer.data(skoff + 60 * idx + 14);

            auto s_xxz_yy = primBuffer.data(skoff + 60 * idx + 15);

            auto s_xxz_yz = primBuffer.data(skoff + 60 * idx + 16);

            auto s_xxz_zz = primBuffer.data(skoff + 60 * idx + 17);

            auto s_xyy_xx = primBuffer.data(skoff + 60 * idx + 18);

            auto s_xyy_xy = primBuffer.data(skoff + 60 * idx + 19);

            auto s_xyy_xz = primBuffer.data(skoff + 60 * idx + 20);

            auto s_xyy_yy = primBuffer.data(skoff + 60 * idx + 21);

            auto s_xyy_yz = primBuffer.data(skoff + 60 * idx + 22);

            auto s_xyy_zz = primBuffer.data(skoff + 60 * idx + 23);

            auto s_xyz_xx = primBuffer.data(skoff + 60 * idx + 24);

            auto s_xyz_xy = primBuffer.data(skoff + 60 * idx + 25);

            auto s_xyz_xz = primBuffer.data(skoff + 60 * idx + 26);

            auto s_xyz_yy = primBuffer.data(skoff + 60 * idx + 27);

            auto s_xyz_yz = primBuffer.data(skoff + 60 * idx + 28);

            auto s_xyz_zz = primBuffer.data(skoff + 60 * idx + 29);

            auto s_xzz_xx = primBuffer.data(skoff + 60 * idx + 30);

            auto s_xzz_xy = primBuffer.data(skoff + 60 * idx + 31);

            auto s_xzz_xz = primBuffer.data(skoff + 60 * idx + 32);

            auto s_xzz_yy = primBuffer.data(skoff + 60 * idx + 33);

            auto s_xzz_yz = primBuffer.data(skoff + 60 * idx + 34);

            auto s_xzz_zz = primBuffer.data(skoff + 60 * idx + 35);

            auto s_yyy_xx = primBuffer.data(skoff + 60 * idx + 36);

            auto s_yyy_xy = primBuffer.data(skoff + 60 * idx + 37);

            auto s_yyy_xz = primBuffer.data(skoff + 60 * idx + 38);

            auto s_yyy_yy = primBuffer.data(skoff + 60 * idx + 39);

            auto s_yyy_yz = primBuffer.data(skoff + 60 * idx + 40);

            auto s_yyy_zz = primBuffer.data(skoff + 60 * idx + 41);

            auto s_yyz_xx = primBuffer.data(skoff + 60 * idx + 42);

            auto s_yyz_xy = primBuffer.data(skoff + 60 * idx + 43);

            auto s_yyz_xz = primBuffer.data(skoff + 60 * idx + 44);

            auto s_yyz_yy = primBuffer.data(skoff + 60 * idx + 45);

            auto s_yyz_yz = primBuffer.data(skoff + 60 * idx + 46);

            auto s_yyz_zz = primBuffer.data(skoff + 60 * idx + 47);

            auto s_yzz_xx = primBuffer.data(skoff + 60 * idx + 48);

            auto s_yzz_xy = primBuffer.data(skoff + 60 * idx + 49);

            auto s_yzz_xz = primBuffer.data(skoff + 60 * idx + 50);

            auto s_yzz_yy = primBuffer.data(skoff + 60 * idx + 51);

            auto s_yzz_yz = primBuffer.data(skoff + 60 * idx + 52);

            auto s_yzz_zz = primBuffer.data(skoff + 60 * idx + 53);

            auto s_zzz_xx = primBuffer.data(skoff + 60 * idx + 54);

            auto s_zzz_xy = primBuffer.data(skoff + 60 * idx + 55);

            auto s_zzz_xz = primBuffer.data(skoff + 60 * idx + 56);

            auto s_zzz_yy = primBuffer.data(skoff + 60 * idx + 57);

            auto s_zzz_yz = primBuffer.data(skoff + 60 * idx + 58);

            auto s_zzz_zz = primBuffer.data(skoff + 60 * idx + 59);

            // set up pointers to (F|T|D) integrals

            auto t_xxx_xx = primBuffer.data(kkoff + 60 * idx);

            auto t_xxx_xy = primBuffer.data(kkoff + 60 * idx + 1);

            auto t_xxx_xz = primBuffer.data(kkoff + 60 * idx + 2);

            auto t_xxx_yy = primBuffer.data(kkoff + 60 * idx + 3);

            auto t_xxx_yz = primBuffer.data(kkoff + 60 * idx + 4);

            auto t_xxx_zz = primBuffer.data(kkoff + 60 * idx + 5);

            auto t_xxy_xx = primBuffer.data(kkoff + 60 * idx + 6);

            auto t_xxy_xy = primBuffer.data(kkoff + 60 * idx + 7);

            auto t_xxy_xz = primBuffer.data(kkoff + 60 * idx + 8);

            auto t_xxy_yy = primBuffer.data(kkoff + 60 * idx + 9);

            auto t_xxy_yz = primBuffer.data(kkoff + 60 * idx + 10);

            auto t_xxy_zz = primBuffer.data(kkoff + 60 * idx + 11);

            auto t_xxz_xx = primBuffer.data(kkoff + 60 * idx + 12);

            auto t_xxz_xy = primBuffer.data(kkoff + 60 * idx + 13);

            auto t_xxz_xz = primBuffer.data(kkoff + 60 * idx + 14);

            auto t_xxz_yy = primBuffer.data(kkoff + 60 * idx + 15);

            auto t_xxz_yz = primBuffer.data(kkoff + 60 * idx + 16);

            auto t_xxz_zz = primBuffer.data(kkoff + 60 * idx + 17);

            auto t_xyy_xx = primBuffer.data(kkoff + 60 * idx + 18);

            auto t_xyy_xy = primBuffer.data(kkoff + 60 * idx + 19);

            auto t_xyy_xz = primBuffer.data(kkoff + 60 * idx + 20);

            auto t_xyy_yy = primBuffer.data(kkoff + 60 * idx + 21);

            auto t_xyy_yz = primBuffer.data(kkoff + 60 * idx + 22);

            auto t_xyy_zz = primBuffer.data(kkoff + 60 * idx + 23);

            auto t_xyz_xx = primBuffer.data(kkoff + 60 * idx + 24);

            auto t_xyz_xy = primBuffer.data(kkoff + 60 * idx + 25);

            auto t_xyz_xz = primBuffer.data(kkoff + 60 * idx + 26);

            auto t_xyz_yy = primBuffer.data(kkoff + 60 * idx + 27);

            auto t_xyz_yz = primBuffer.data(kkoff + 60 * idx + 28);

            auto t_xyz_zz = primBuffer.data(kkoff + 60 * idx + 29);

            auto t_xzz_xx = primBuffer.data(kkoff + 60 * idx + 30);

            auto t_xzz_xy = primBuffer.data(kkoff + 60 * idx + 31);

            auto t_xzz_xz = primBuffer.data(kkoff + 60 * idx + 32);

            auto t_xzz_yy = primBuffer.data(kkoff + 60 * idx + 33);

            auto t_xzz_yz = primBuffer.data(kkoff + 60 * idx + 34);

            auto t_xzz_zz = primBuffer.data(kkoff + 60 * idx + 35);

            auto t_yyy_xx = primBuffer.data(kkoff + 60 * idx + 36);

            auto t_yyy_xy = primBuffer.data(kkoff + 60 * idx + 37);

            auto t_yyy_xz = primBuffer.data(kkoff + 60 * idx + 38);

            auto t_yyy_yy = primBuffer.data(kkoff + 60 * idx + 39);

            auto t_yyy_yz = primBuffer.data(kkoff + 60 * idx + 40);

            auto t_yyy_zz = primBuffer.data(kkoff + 60 * idx + 41);

            auto t_yyz_xx = primBuffer.data(kkoff + 60 * idx + 42);

            auto t_yyz_xy = primBuffer.data(kkoff + 60 * idx + 43);

            auto t_yyz_xz = primBuffer.data(kkoff + 60 * idx + 44);

            auto t_yyz_yy = primBuffer.data(kkoff + 60 * idx + 45);

            auto t_yyz_yz = primBuffer.data(kkoff + 60 * idx + 46);

            auto t_yyz_zz = primBuffer.data(kkoff + 60 * idx + 47);

            auto t_yzz_xx = primBuffer.data(kkoff + 60 * idx + 48);

            auto t_yzz_xy = primBuffer.data(kkoff + 60 * idx + 49);

            auto t_yzz_xz = primBuffer.data(kkoff + 60 * idx + 50);

            auto t_yzz_yy = primBuffer.data(kkoff + 60 * idx + 51);

            auto t_yzz_yz = primBuffer.data(kkoff + 60 * idx + 52);

            auto t_yzz_zz = primBuffer.data(kkoff + 60 * idx + 53);

            auto t_zzz_xx = primBuffer.data(kkoff + 60 * idx + 54);

            auto t_zzz_xy = primBuffer.data(kkoff + 60 * idx + 55);

            auto t_zzz_xz = primBuffer.data(kkoff + 60 * idx + 56);

            auto t_zzz_yy = primBuffer.data(kkoff + 60 * idx + 57);

            auto t_zzz_yz = primBuffer.data(kkoff + 60 * idx + 58);

            auto t_zzz_zz = primBuffer.data(kkoff + 60 * idx + 59);

            // set up pointers to (D|F) integrals

            auto s_xx_xxx = primBuffer.data(s2off + 60 * idx);

            auto s_xx_xxy = primBuffer.data(s2off + 60 * idx + 1);

            auto s_xx_xxz = primBuffer.data(s2off + 60 * idx + 2);

            auto s_xx_xyy = primBuffer.data(s2off + 60 * idx + 3);

            auto s_xx_xyz = primBuffer.data(s2off + 60 * idx + 4);

            auto s_xx_xzz = primBuffer.data(s2off + 60 * idx + 5);

            auto s_xx_yyy = primBuffer.data(s2off + 60 * idx + 6);

            auto s_xx_yyz = primBuffer.data(s2off + 60 * idx + 7);

            auto s_xx_yzz = primBuffer.data(s2off + 60 * idx + 8);

            auto s_xx_zzz = primBuffer.data(s2off + 60 * idx + 9);

            auto s_xy_xxx = primBuffer.data(s2off + 60 * idx + 10);

            auto s_xy_xxy = primBuffer.data(s2off + 60 * idx + 11);

            auto s_xy_xxz = primBuffer.data(s2off + 60 * idx + 12);

            auto s_xy_xyy = primBuffer.data(s2off + 60 * idx + 13);

            auto s_xy_xyz = primBuffer.data(s2off + 60 * idx + 14);

            auto s_xy_xzz = primBuffer.data(s2off + 60 * idx + 15);

            auto s_xy_yyy = primBuffer.data(s2off + 60 * idx + 16);

            auto s_xy_yyz = primBuffer.data(s2off + 60 * idx + 17);

            auto s_xy_yzz = primBuffer.data(s2off + 60 * idx + 18);

            auto s_xy_zzz = primBuffer.data(s2off + 60 * idx + 19);

            auto s_xz_xxx = primBuffer.data(s2off + 60 * idx + 20);

            auto s_xz_xxy = primBuffer.data(s2off + 60 * idx + 21);

            auto s_xz_xxz = primBuffer.data(s2off + 60 * idx + 22);

            auto s_xz_xyy = primBuffer.data(s2off + 60 * idx + 23);

            auto s_xz_xyz = primBuffer.data(s2off + 60 * idx + 24);

            auto s_xz_xzz = primBuffer.data(s2off + 60 * idx + 25);

            auto s_xz_yyy = primBuffer.data(s2off + 60 * idx + 26);

            auto s_xz_yyz = primBuffer.data(s2off + 60 * idx + 27);

            auto s_xz_yzz = primBuffer.data(s2off + 60 * idx + 28);

            auto s_xz_zzz = primBuffer.data(s2off + 60 * idx + 29);

            auto s_yy_xxx = primBuffer.data(s2off + 60 * idx + 30);

            auto s_yy_xxy = primBuffer.data(s2off + 60 * idx + 31);

            auto s_yy_xxz = primBuffer.data(s2off + 60 * idx + 32);

            auto s_yy_xyy = primBuffer.data(s2off + 60 * idx + 33);

            auto s_yy_xyz = primBuffer.data(s2off + 60 * idx + 34);

            auto s_yy_xzz = primBuffer.data(s2off + 60 * idx + 35);

            auto s_yy_yyy = primBuffer.data(s2off + 60 * idx + 36);

            auto s_yy_yyz = primBuffer.data(s2off + 60 * idx + 37);

            auto s_yy_yzz = primBuffer.data(s2off + 60 * idx + 38);

            auto s_yy_zzz = primBuffer.data(s2off + 60 * idx + 39);

            auto s_yz_xxx = primBuffer.data(s2off + 60 * idx + 40);

            auto s_yz_xxy = primBuffer.data(s2off + 60 * idx + 41);

            auto s_yz_xxz = primBuffer.data(s2off + 60 * idx + 42);

            auto s_yz_xyy = primBuffer.data(s2off + 60 * idx + 43);

            auto s_yz_xyz = primBuffer.data(s2off + 60 * idx + 44);

            auto s_yz_xzz = primBuffer.data(s2off + 60 * idx + 45);

            auto s_yz_yyy = primBuffer.data(s2off + 60 * idx + 46);

            auto s_yz_yyz = primBuffer.data(s2off + 60 * idx + 47);

            auto s_yz_yzz = primBuffer.data(s2off + 60 * idx + 48);

            auto s_yz_zzz = primBuffer.data(s2off + 60 * idx + 49);

            auto s_zz_xxx = primBuffer.data(s2off + 60 * idx + 50);

            auto s_zz_xxy = primBuffer.data(s2off + 60 * idx + 51);

            auto s_zz_xxz = primBuffer.data(s2off + 60 * idx + 52);

            auto s_zz_xyy = primBuffer.data(s2off + 60 * idx + 53);

            auto s_zz_xyz = primBuffer.data(s2off + 60 * idx + 54);

            auto s_zz_xzz = primBuffer.data(s2off + 60 * idx + 55);

            auto s_zz_yyy = primBuffer.data(s2off + 60 * idx + 56);

            auto s_zz_yyz = primBuffer.data(s2off + 60 * idx + 57);

            auto s_zz_yzz = primBuffer.data(s2off + 60 * idx + 58);

            auto s_zz_zzz = primBuffer.data(s2off + 60 * idx + 59);

            // set up pointers to (D|T|F) integrals

            auto t_xx_xxx = primBuffer.data(k2off + 60 * idx);

            auto t_xx_xxy = primBuffer.data(k2off + 60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(k2off + 60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(k2off + 60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(k2off + 60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(k2off + 60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(k2off + 60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(k2off + 60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(k2off + 60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(k2off + 60 * idx + 9);

            auto t_xy_xxx = primBuffer.data(k2off + 60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(k2off + 60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(k2off + 60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(k2off + 60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(k2off + 60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(k2off + 60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(k2off + 60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(k2off + 60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(k2off + 60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(k2off + 60 * idx + 19);

            auto t_xz_xxx = primBuffer.data(k2off + 60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(k2off + 60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(k2off + 60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(k2off + 60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(k2off + 60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(k2off + 60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(k2off + 60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(k2off + 60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(k2off + 60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(k2off + 60 * idx + 29);

            auto t_yy_xxx = primBuffer.data(k2off + 60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(k2off + 60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(k2off + 60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(k2off + 60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(k2off + 60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(k2off + 60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(k2off + 60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(k2off + 60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(k2off + 60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(k2off + 60 * idx + 39);

            auto t_yz_xxx = primBuffer.data(k2off + 60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(k2off + 60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(k2off + 60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(k2off + 60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(k2off + 60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(k2off + 60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(k2off + 60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(k2off + 60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(k2off + 60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(k2off + 60 * idx + 49);

            auto t_zz_xxx = primBuffer.data(k2off + 60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(k2off + 60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(k2off + 60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(k2off + 60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(k2off + 60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(k2off + 60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(k2off + 60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(k2off + 60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(k2off + 60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(k2off + 60 * idx + 59);

            // set up pointers to (F|F) integrals

            auto s_xxx_xxx = primBuffer.data(s1off + 100 * idx);

            auto s_xxx_xxy = primBuffer.data(s1off + 100 * idx + 1);

            auto s_xxx_xxz = primBuffer.data(s1off + 100 * idx + 2);

            auto s_xxx_xyy = primBuffer.data(s1off + 100 * idx + 3);

            auto s_xxx_xyz = primBuffer.data(s1off + 100 * idx + 4);

            auto s_xxx_xzz = primBuffer.data(s1off + 100 * idx + 5);

            auto s_xxx_yyy = primBuffer.data(s1off + 100 * idx + 6);

            auto s_xxx_yyz = primBuffer.data(s1off + 100 * idx + 7);

            auto s_xxx_yzz = primBuffer.data(s1off + 100 * idx + 8);

            auto s_xxx_zzz = primBuffer.data(s1off + 100 * idx + 9);

            auto s_xxy_xxx = primBuffer.data(s1off + 100 * idx + 10);

            auto s_xxy_xxy = primBuffer.data(s1off + 100 * idx + 11);

            auto s_xxy_xxz = primBuffer.data(s1off + 100 * idx + 12);

            auto s_xxy_xyy = primBuffer.data(s1off + 100 * idx + 13);

            auto s_xxy_xyz = primBuffer.data(s1off + 100 * idx + 14);

            auto s_xxy_xzz = primBuffer.data(s1off + 100 * idx + 15);

            auto s_xxy_yyy = primBuffer.data(s1off + 100 * idx + 16);

            auto s_xxy_yyz = primBuffer.data(s1off + 100 * idx + 17);

            auto s_xxy_yzz = primBuffer.data(s1off + 100 * idx + 18);

            auto s_xxy_zzz = primBuffer.data(s1off + 100 * idx + 19);

            auto s_xxz_xxx = primBuffer.data(s1off + 100 * idx + 20);

            auto s_xxz_xxy = primBuffer.data(s1off + 100 * idx + 21);

            auto s_xxz_xxz = primBuffer.data(s1off + 100 * idx + 22);

            auto s_xxz_xyy = primBuffer.data(s1off + 100 * idx + 23);

            auto s_xxz_xyz = primBuffer.data(s1off + 100 * idx + 24);

            auto s_xxz_xzz = primBuffer.data(s1off + 100 * idx + 25);

            auto s_xxz_yyy = primBuffer.data(s1off + 100 * idx + 26);

            auto s_xxz_yyz = primBuffer.data(s1off + 100 * idx + 27);

            auto s_xxz_yzz = primBuffer.data(s1off + 100 * idx + 28);

            auto s_xxz_zzz = primBuffer.data(s1off + 100 * idx + 29);

            auto s_xyy_xxx = primBuffer.data(s1off + 100 * idx + 30);

            auto s_xyy_xxy = primBuffer.data(s1off + 100 * idx + 31);

            auto s_xyy_xxz = primBuffer.data(s1off + 100 * idx + 32);

            auto s_xyy_xyy = primBuffer.data(s1off + 100 * idx + 33);

            auto s_xyy_xyz = primBuffer.data(s1off + 100 * idx + 34);

            auto s_xyy_xzz = primBuffer.data(s1off + 100 * idx + 35);

            auto s_xyy_yyy = primBuffer.data(s1off + 100 * idx + 36);

            auto s_xyy_yyz = primBuffer.data(s1off + 100 * idx + 37);

            auto s_xyy_yzz = primBuffer.data(s1off + 100 * idx + 38);

            auto s_xyy_zzz = primBuffer.data(s1off + 100 * idx + 39);

            auto s_xyz_xxx = primBuffer.data(s1off + 100 * idx + 40);

            auto s_xyz_xxy = primBuffer.data(s1off + 100 * idx + 41);

            auto s_xyz_xxz = primBuffer.data(s1off + 100 * idx + 42);

            auto s_xyz_xyy = primBuffer.data(s1off + 100 * idx + 43);

            auto s_xyz_xyz = primBuffer.data(s1off + 100 * idx + 44);

            auto s_xyz_xzz = primBuffer.data(s1off + 100 * idx + 45);

            auto s_xyz_yyy = primBuffer.data(s1off + 100 * idx + 46);

            auto s_xyz_yyz = primBuffer.data(s1off + 100 * idx + 47);

            auto s_xyz_yzz = primBuffer.data(s1off + 100 * idx + 48);

            auto s_xyz_zzz = primBuffer.data(s1off + 100 * idx + 49);

            auto s_xzz_xxx = primBuffer.data(s1off + 100 * idx + 50);

            auto s_xzz_xxy = primBuffer.data(s1off + 100 * idx + 51);

            auto s_xzz_xxz = primBuffer.data(s1off + 100 * idx + 52);

            auto s_xzz_xyy = primBuffer.data(s1off + 100 * idx + 53);

            auto s_xzz_xyz = primBuffer.data(s1off + 100 * idx + 54);

            auto s_xzz_xzz = primBuffer.data(s1off + 100 * idx + 55);

            auto s_xzz_yyy = primBuffer.data(s1off + 100 * idx + 56);

            auto s_xzz_yyz = primBuffer.data(s1off + 100 * idx + 57);

            auto s_xzz_yzz = primBuffer.data(s1off + 100 * idx + 58);

            auto s_xzz_zzz = primBuffer.data(s1off + 100 * idx + 59);

            auto s_yyy_xxx = primBuffer.data(s1off + 100 * idx + 60);

            auto s_yyy_xxy = primBuffer.data(s1off + 100 * idx + 61);

            auto s_yyy_xxz = primBuffer.data(s1off + 100 * idx + 62);

            auto s_yyy_xyy = primBuffer.data(s1off + 100 * idx + 63);

            auto s_yyy_xyz = primBuffer.data(s1off + 100 * idx + 64);

            auto s_yyy_xzz = primBuffer.data(s1off + 100 * idx + 65);

            auto s_yyy_yyy = primBuffer.data(s1off + 100 * idx + 66);

            auto s_yyy_yyz = primBuffer.data(s1off + 100 * idx + 67);

            auto s_yyy_yzz = primBuffer.data(s1off + 100 * idx + 68);

            auto s_yyy_zzz = primBuffer.data(s1off + 100 * idx + 69);

            auto s_yyz_xxx = primBuffer.data(s1off + 100 * idx + 70);

            auto s_yyz_xxy = primBuffer.data(s1off + 100 * idx + 71);

            auto s_yyz_xxz = primBuffer.data(s1off + 100 * idx + 72);

            auto s_yyz_xyy = primBuffer.data(s1off + 100 * idx + 73);

            auto s_yyz_xyz = primBuffer.data(s1off + 100 * idx + 74);

            auto s_yyz_xzz = primBuffer.data(s1off + 100 * idx + 75);

            auto s_yyz_yyy = primBuffer.data(s1off + 100 * idx + 76);

            auto s_yyz_yyz = primBuffer.data(s1off + 100 * idx + 77);

            auto s_yyz_yzz = primBuffer.data(s1off + 100 * idx + 78);

            auto s_yyz_zzz = primBuffer.data(s1off + 100 * idx + 79);

            auto s_yzz_xxx = primBuffer.data(s1off + 100 * idx + 80);

            auto s_yzz_xxy = primBuffer.data(s1off + 100 * idx + 81);

            auto s_yzz_xxz = primBuffer.data(s1off + 100 * idx + 82);

            auto s_yzz_xyy = primBuffer.data(s1off + 100 * idx + 83);

            auto s_yzz_xyz = primBuffer.data(s1off + 100 * idx + 84);

            auto s_yzz_xzz = primBuffer.data(s1off + 100 * idx + 85);

            auto s_yzz_yyy = primBuffer.data(s1off + 100 * idx + 86);

            auto s_yzz_yyz = primBuffer.data(s1off + 100 * idx + 87);

            auto s_yzz_yzz = primBuffer.data(s1off + 100 * idx + 88);

            auto s_yzz_zzz = primBuffer.data(s1off + 100 * idx + 89);

            auto s_zzz_xxx = primBuffer.data(s1off + 100 * idx + 90);

            auto s_zzz_xxy = primBuffer.data(s1off + 100 * idx + 91);

            auto s_zzz_xxz = primBuffer.data(s1off + 100 * idx + 92);

            auto s_zzz_xyy = primBuffer.data(s1off + 100 * idx + 93);

            auto s_zzz_xyz = primBuffer.data(s1off + 100 * idx + 94);

            auto s_zzz_xzz = primBuffer.data(s1off + 100 * idx + 95);

            auto s_zzz_yyy = primBuffer.data(s1off + 100 * idx + 96);

            auto s_zzz_yyz = primBuffer.data(s1off + 100 * idx + 97);

            auto s_zzz_yzz = primBuffer.data(s1off + 100 * idx + 98);

            auto s_zzz_zzz = primBuffer.data(s1off + 100 * idx + 99);

            // set up pointers to (F|T|F) integrals

            auto t_xxx_xxx = primBuffer.data(k1off + 100 * idx);

            auto t_xxx_xxy = primBuffer.data(k1off + 100 * idx + 1);

            auto t_xxx_xxz = primBuffer.data(k1off + 100 * idx + 2);

            auto t_xxx_xyy = primBuffer.data(k1off + 100 * idx + 3);

            auto t_xxx_xyz = primBuffer.data(k1off + 100 * idx + 4);

            auto t_xxx_xzz = primBuffer.data(k1off + 100 * idx + 5);

            auto t_xxx_yyy = primBuffer.data(k1off + 100 * idx + 6);

            auto t_xxx_yyz = primBuffer.data(k1off + 100 * idx + 7);

            auto t_xxx_yzz = primBuffer.data(k1off + 100 * idx + 8);

            auto t_xxx_zzz = primBuffer.data(k1off + 100 * idx + 9);

            auto t_xxy_xxx = primBuffer.data(k1off + 100 * idx + 10);

            auto t_xxy_xxy = primBuffer.data(k1off + 100 * idx + 11);

            auto t_xxy_xxz = primBuffer.data(k1off + 100 * idx + 12);

            auto t_xxy_xyy = primBuffer.data(k1off + 100 * idx + 13);

            auto t_xxy_xyz = primBuffer.data(k1off + 100 * idx + 14);

            auto t_xxy_xzz = primBuffer.data(k1off + 100 * idx + 15);

            auto t_xxy_yyy = primBuffer.data(k1off + 100 * idx + 16);

            auto t_xxy_yyz = primBuffer.data(k1off + 100 * idx + 17);

            auto t_xxy_yzz = primBuffer.data(k1off + 100 * idx + 18);

            auto t_xxy_zzz = primBuffer.data(k1off + 100 * idx + 19);

            auto t_xxz_xxx = primBuffer.data(k1off + 100 * idx + 20);

            auto t_xxz_xxy = primBuffer.data(k1off + 100 * idx + 21);

            auto t_xxz_xxz = primBuffer.data(k1off + 100 * idx + 22);

            auto t_xxz_xyy = primBuffer.data(k1off + 100 * idx + 23);

            auto t_xxz_xyz = primBuffer.data(k1off + 100 * idx + 24);

            auto t_xxz_xzz = primBuffer.data(k1off + 100 * idx + 25);

            auto t_xxz_yyy = primBuffer.data(k1off + 100 * idx + 26);

            auto t_xxz_yyz = primBuffer.data(k1off + 100 * idx + 27);

            auto t_xxz_yzz = primBuffer.data(k1off + 100 * idx + 28);

            auto t_xxz_zzz = primBuffer.data(k1off + 100 * idx + 29);

            auto t_xyy_xxx = primBuffer.data(k1off + 100 * idx + 30);

            auto t_xyy_xxy = primBuffer.data(k1off + 100 * idx + 31);

            auto t_xyy_xxz = primBuffer.data(k1off + 100 * idx + 32);

            auto t_xyy_xyy = primBuffer.data(k1off + 100 * idx + 33);

            auto t_xyy_xyz = primBuffer.data(k1off + 100 * idx + 34);

            auto t_xyy_xzz = primBuffer.data(k1off + 100 * idx + 35);

            auto t_xyy_yyy = primBuffer.data(k1off + 100 * idx + 36);

            auto t_xyy_yyz = primBuffer.data(k1off + 100 * idx + 37);

            auto t_xyy_yzz = primBuffer.data(k1off + 100 * idx + 38);

            auto t_xyy_zzz = primBuffer.data(k1off + 100 * idx + 39);

            auto t_xyz_xxx = primBuffer.data(k1off + 100 * idx + 40);

            auto t_xyz_xxy = primBuffer.data(k1off + 100 * idx + 41);

            auto t_xyz_xxz = primBuffer.data(k1off + 100 * idx + 42);

            auto t_xyz_xyy = primBuffer.data(k1off + 100 * idx + 43);

            auto t_xyz_xyz = primBuffer.data(k1off + 100 * idx + 44);

            auto t_xyz_xzz = primBuffer.data(k1off + 100 * idx + 45);

            auto t_xyz_yyy = primBuffer.data(k1off + 100 * idx + 46);

            auto t_xyz_yyz = primBuffer.data(k1off + 100 * idx + 47);

            auto t_xyz_yzz = primBuffer.data(k1off + 100 * idx + 48);

            auto t_xyz_zzz = primBuffer.data(k1off + 100 * idx + 49);

            auto t_xzz_xxx = primBuffer.data(k1off + 100 * idx + 50);

            auto t_xzz_xxy = primBuffer.data(k1off + 100 * idx + 51);

            auto t_xzz_xxz = primBuffer.data(k1off + 100 * idx + 52);

            auto t_xzz_xyy = primBuffer.data(k1off + 100 * idx + 53);

            auto t_xzz_xyz = primBuffer.data(k1off + 100 * idx + 54);

            auto t_xzz_xzz = primBuffer.data(k1off + 100 * idx + 55);

            auto t_xzz_yyy = primBuffer.data(k1off + 100 * idx + 56);

            auto t_xzz_yyz = primBuffer.data(k1off + 100 * idx + 57);

            auto t_xzz_yzz = primBuffer.data(k1off + 100 * idx + 58);

            auto t_xzz_zzz = primBuffer.data(k1off + 100 * idx + 59);

            auto t_yyy_xxx = primBuffer.data(k1off + 100 * idx + 60);

            auto t_yyy_xxy = primBuffer.data(k1off + 100 * idx + 61);

            auto t_yyy_xxz = primBuffer.data(k1off + 100 * idx + 62);

            auto t_yyy_xyy = primBuffer.data(k1off + 100 * idx + 63);

            auto t_yyy_xyz = primBuffer.data(k1off + 100 * idx + 64);

            auto t_yyy_xzz = primBuffer.data(k1off + 100 * idx + 65);

            auto t_yyy_yyy = primBuffer.data(k1off + 100 * idx + 66);

            auto t_yyy_yyz = primBuffer.data(k1off + 100 * idx + 67);

            auto t_yyy_yzz = primBuffer.data(k1off + 100 * idx + 68);

            auto t_yyy_zzz = primBuffer.data(k1off + 100 * idx + 69);

            auto t_yyz_xxx = primBuffer.data(k1off + 100 * idx + 70);

            auto t_yyz_xxy = primBuffer.data(k1off + 100 * idx + 71);

            auto t_yyz_xxz = primBuffer.data(k1off + 100 * idx + 72);

            auto t_yyz_xyy = primBuffer.data(k1off + 100 * idx + 73);

            auto t_yyz_xyz = primBuffer.data(k1off + 100 * idx + 74);

            auto t_yyz_xzz = primBuffer.data(k1off + 100 * idx + 75);

            auto t_yyz_yyy = primBuffer.data(k1off + 100 * idx + 76);

            auto t_yyz_yyz = primBuffer.data(k1off + 100 * idx + 77);

            auto t_yyz_yzz = primBuffer.data(k1off + 100 * idx + 78);

            auto t_yyz_zzz = primBuffer.data(k1off + 100 * idx + 79);

            auto t_yzz_xxx = primBuffer.data(k1off + 100 * idx + 80);

            auto t_yzz_xxy = primBuffer.data(k1off + 100 * idx + 81);

            auto t_yzz_xxz = primBuffer.data(k1off + 100 * idx + 82);

            auto t_yzz_xyy = primBuffer.data(k1off + 100 * idx + 83);

            auto t_yzz_xyz = primBuffer.data(k1off + 100 * idx + 84);

            auto t_yzz_xzz = primBuffer.data(k1off + 100 * idx + 85);

            auto t_yzz_yyy = primBuffer.data(k1off + 100 * idx + 86);

            auto t_yzz_yyz = primBuffer.data(k1off + 100 * idx + 87);

            auto t_yzz_yzz = primBuffer.data(k1off + 100 * idx + 88);

            auto t_yzz_zzz = primBuffer.data(k1off + 100 * idx + 89);

            auto t_zzz_xxx = primBuffer.data(k1off + 100 * idx + 90);

            auto t_zzz_xxy = primBuffer.data(k1off + 100 * idx + 91);

            auto t_zzz_xxz = primBuffer.data(k1off + 100 * idx + 92);

            auto t_zzz_xyy = primBuffer.data(k1off + 100 * idx + 93);

            auto t_zzz_xyz = primBuffer.data(k1off + 100 * idx + 94);

            auto t_zzz_xzz = primBuffer.data(k1off + 100 * idx + 95);

            auto t_zzz_yyy = primBuffer.data(k1off + 100 * idx + 96);

            auto t_zzz_yyz = primBuffer.data(k1off + 100 * idx + 97);

            auto t_zzz_yzz = primBuffer.data(k1off + 100 * idx + 98);

            auto t_zzz_zzz = primBuffer.data(k1off + 100 * idx + 99);

            // set up pointers to (G|F) integrals

            auto s_xxxx_xxx = primBuffer.data(soff + 150 * idx);

            auto s_xxxx_xxy = primBuffer.data(soff + 150 * idx + 1);

            auto s_xxxx_xxz = primBuffer.data(soff + 150 * idx + 2);

            auto s_xxxx_xyy = primBuffer.data(soff + 150 * idx + 3);

            auto s_xxxx_xyz = primBuffer.data(soff + 150 * idx + 4);

            auto s_xxxx_xzz = primBuffer.data(soff + 150 * idx + 5);

            auto s_xxxx_yyy = primBuffer.data(soff + 150 * idx + 6);

            auto s_xxxx_yyz = primBuffer.data(soff + 150 * idx + 7);

            auto s_xxxx_yzz = primBuffer.data(soff + 150 * idx + 8);

            auto s_xxxx_zzz = primBuffer.data(soff + 150 * idx + 9);

            auto s_xxxy_xxx = primBuffer.data(soff + 150 * idx + 10);

            auto s_xxxy_xxy = primBuffer.data(soff + 150 * idx + 11);

            auto s_xxxy_xxz = primBuffer.data(soff + 150 * idx + 12);

            auto s_xxxy_xyy = primBuffer.data(soff + 150 * idx + 13);

            auto s_xxxy_xyz = primBuffer.data(soff + 150 * idx + 14);

            auto s_xxxy_xzz = primBuffer.data(soff + 150 * idx + 15);

            auto s_xxxy_yyy = primBuffer.data(soff + 150 * idx + 16);

            auto s_xxxy_yyz = primBuffer.data(soff + 150 * idx + 17);

            auto s_xxxy_yzz = primBuffer.data(soff + 150 * idx + 18);

            auto s_xxxy_zzz = primBuffer.data(soff + 150 * idx + 19);

            auto s_xxxz_xxx = primBuffer.data(soff + 150 * idx + 20);

            auto s_xxxz_xxy = primBuffer.data(soff + 150 * idx + 21);

            auto s_xxxz_xxz = primBuffer.data(soff + 150 * idx + 22);

            auto s_xxxz_xyy = primBuffer.data(soff + 150 * idx + 23);

            auto s_xxxz_xyz = primBuffer.data(soff + 150 * idx + 24);

            auto s_xxxz_xzz = primBuffer.data(soff + 150 * idx + 25);

            auto s_xxxz_yyy = primBuffer.data(soff + 150 * idx + 26);

            auto s_xxxz_yyz = primBuffer.data(soff + 150 * idx + 27);

            auto s_xxxz_yzz = primBuffer.data(soff + 150 * idx + 28);

            auto s_xxxz_zzz = primBuffer.data(soff + 150 * idx + 29);

            auto s_xxyy_xxx = primBuffer.data(soff + 150 * idx + 30);

            auto s_xxyy_xxy = primBuffer.data(soff + 150 * idx + 31);

            auto s_xxyy_xxz = primBuffer.data(soff + 150 * idx + 32);

            auto s_xxyy_xyy = primBuffer.data(soff + 150 * idx + 33);

            auto s_xxyy_xyz = primBuffer.data(soff + 150 * idx + 34);

            auto s_xxyy_xzz = primBuffer.data(soff + 150 * idx + 35);

            auto s_xxyy_yyy = primBuffer.data(soff + 150 * idx + 36);

            auto s_xxyy_yyz = primBuffer.data(soff + 150 * idx + 37);

            auto s_xxyy_yzz = primBuffer.data(soff + 150 * idx + 38);

            auto s_xxyy_zzz = primBuffer.data(soff + 150 * idx + 39);

            auto s_xxyz_xxx = primBuffer.data(soff + 150 * idx + 40);

            auto s_xxyz_xxy = primBuffer.data(soff + 150 * idx + 41);

            auto s_xxyz_xxz = primBuffer.data(soff + 150 * idx + 42);

            auto s_xxyz_xyy = primBuffer.data(soff + 150 * idx + 43);

            auto s_xxyz_xyz = primBuffer.data(soff + 150 * idx + 44);

            auto s_xxyz_xzz = primBuffer.data(soff + 150 * idx + 45);

            auto s_xxyz_yyy = primBuffer.data(soff + 150 * idx + 46);

            auto s_xxyz_yyz = primBuffer.data(soff + 150 * idx + 47);

            auto s_xxyz_yzz = primBuffer.data(soff + 150 * idx + 48);

            auto s_xxyz_zzz = primBuffer.data(soff + 150 * idx + 49);

            auto s_xxzz_xxx = primBuffer.data(soff + 150 * idx + 50);

            auto s_xxzz_xxy = primBuffer.data(soff + 150 * idx + 51);

            auto s_xxzz_xxz = primBuffer.data(soff + 150 * idx + 52);

            auto s_xxzz_xyy = primBuffer.data(soff + 150 * idx + 53);

            auto s_xxzz_xyz = primBuffer.data(soff + 150 * idx + 54);

            auto s_xxzz_xzz = primBuffer.data(soff + 150 * idx + 55);

            auto s_xxzz_yyy = primBuffer.data(soff + 150 * idx + 56);

            auto s_xxzz_yyz = primBuffer.data(soff + 150 * idx + 57);

            auto s_xxzz_yzz = primBuffer.data(soff + 150 * idx + 58);

            auto s_xxzz_zzz = primBuffer.data(soff + 150 * idx + 59);

            auto s_xyyy_xxx = primBuffer.data(soff + 150 * idx + 60);

            auto s_xyyy_xxy = primBuffer.data(soff + 150 * idx + 61);

            auto s_xyyy_xxz = primBuffer.data(soff + 150 * idx + 62);

            auto s_xyyy_xyy = primBuffer.data(soff + 150 * idx + 63);

            auto s_xyyy_xyz = primBuffer.data(soff + 150 * idx + 64);

            auto s_xyyy_xzz = primBuffer.data(soff + 150 * idx + 65);

            auto s_xyyy_yyy = primBuffer.data(soff + 150 * idx + 66);

            auto s_xyyy_yyz = primBuffer.data(soff + 150 * idx + 67);

            auto s_xyyy_yzz = primBuffer.data(soff + 150 * idx + 68);

            auto s_xyyy_zzz = primBuffer.data(soff + 150 * idx + 69);

            auto s_xyyz_xxx = primBuffer.data(soff + 150 * idx + 70);

            auto s_xyyz_xxy = primBuffer.data(soff + 150 * idx + 71);

            auto s_xyyz_xxz = primBuffer.data(soff + 150 * idx + 72);

            auto s_xyyz_xyy = primBuffer.data(soff + 150 * idx + 73);

            auto s_xyyz_xyz = primBuffer.data(soff + 150 * idx + 74);

            auto s_xyyz_xzz = primBuffer.data(soff + 150 * idx + 75);

            auto s_xyyz_yyy = primBuffer.data(soff + 150 * idx + 76);

            auto s_xyyz_yyz = primBuffer.data(soff + 150 * idx + 77);

            auto s_xyyz_yzz = primBuffer.data(soff + 150 * idx + 78);

            auto s_xyyz_zzz = primBuffer.data(soff + 150 * idx + 79);

            auto s_xyzz_xxx = primBuffer.data(soff + 150 * idx + 80);

            auto s_xyzz_xxy = primBuffer.data(soff + 150 * idx + 81);

            auto s_xyzz_xxz = primBuffer.data(soff + 150 * idx + 82);

            auto s_xyzz_xyy = primBuffer.data(soff + 150 * idx + 83);

            auto s_xyzz_xyz = primBuffer.data(soff + 150 * idx + 84);

            auto s_xyzz_xzz = primBuffer.data(soff + 150 * idx + 85);

            auto s_xyzz_yyy = primBuffer.data(soff + 150 * idx + 86);

            auto s_xyzz_yyz = primBuffer.data(soff + 150 * idx + 87);

            auto s_xyzz_yzz = primBuffer.data(soff + 150 * idx + 88);

            auto s_xyzz_zzz = primBuffer.data(soff + 150 * idx + 89);

            auto s_xzzz_xxx = primBuffer.data(soff + 150 * idx + 90);

            auto s_xzzz_xxy = primBuffer.data(soff + 150 * idx + 91);

            auto s_xzzz_xxz = primBuffer.data(soff + 150 * idx + 92);

            auto s_xzzz_xyy = primBuffer.data(soff + 150 * idx + 93);

            auto s_xzzz_xyz = primBuffer.data(soff + 150 * idx + 94);

            auto s_xzzz_xzz = primBuffer.data(soff + 150 * idx + 95);

            auto s_xzzz_yyy = primBuffer.data(soff + 150 * idx + 96);

            auto s_xzzz_yyz = primBuffer.data(soff + 150 * idx + 97);

            auto s_xzzz_yzz = primBuffer.data(soff + 150 * idx + 98);

            auto s_xzzz_zzz = primBuffer.data(soff + 150 * idx + 99);

            auto s_yyyy_xxx = primBuffer.data(soff + 150 * idx + 100);

            auto s_yyyy_xxy = primBuffer.data(soff + 150 * idx + 101);

            auto s_yyyy_xxz = primBuffer.data(soff + 150 * idx + 102);

            auto s_yyyy_xyy = primBuffer.data(soff + 150 * idx + 103);

            auto s_yyyy_xyz = primBuffer.data(soff + 150 * idx + 104);

            auto s_yyyy_xzz = primBuffer.data(soff + 150 * idx + 105);

            auto s_yyyy_yyy = primBuffer.data(soff + 150 * idx + 106);

            auto s_yyyy_yyz = primBuffer.data(soff + 150 * idx + 107);

            auto s_yyyy_yzz = primBuffer.data(soff + 150 * idx + 108);

            auto s_yyyy_zzz = primBuffer.data(soff + 150 * idx + 109);

            auto s_yyyz_xxx = primBuffer.data(soff + 150 * idx + 110);

            auto s_yyyz_xxy = primBuffer.data(soff + 150 * idx + 111);

            auto s_yyyz_xxz = primBuffer.data(soff + 150 * idx + 112);

            auto s_yyyz_xyy = primBuffer.data(soff + 150 * idx + 113);

            auto s_yyyz_xyz = primBuffer.data(soff + 150 * idx + 114);

            auto s_yyyz_xzz = primBuffer.data(soff + 150 * idx + 115);

            auto s_yyyz_yyy = primBuffer.data(soff + 150 * idx + 116);

            auto s_yyyz_yyz = primBuffer.data(soff + 150 * idx + 117);

            auto s_yyyz_yzz = primBuffer.data(soff + 150 * idx + 118);

            auto s_yyyz_zzz = primBuffer.data(soff + 150 * idx + 119);

            auto s_yyzz_xxx = primBuffer.data(soff + 150 * idx + 120);

            auto s_yyzz_xxy = primBuffer.data(soff + 150 * idx + 121);

            auto s_yyzz_xxz = primBuffer.data(soff + 150 * idx + 122);

            auto s_yyzz_xyy = primBuffer.data(soff + 150 * idx + 123);

            auto s_yyzz_xyz = primBuffer.data(soff + 150 * idx + 124);

            auto s_yyzz_xzz = primBuffer.data(soff + 150 * idx + 125);

            auto s_yyzz_yyy = primBuffer.data(soff + 150 * idx + 126);

            auto s_yyzz_yyz = primBuffer.data(soff + 150 * idx + 127);

            auto s_yyzz_yzz = primBuffer.data(soff + 150 * idx + 128);

            auto s_yyzz_zzz = primBuffer.data(soff + 150 * idx + 129);

            auto s_yzzz_xxx = primBuffer.data(soff + 150 * idx + 130);

            auto s_yzzz_xxy = primBuffer.data(soff + 150 * idx + 131);

            auto s_yzzz_xxz = primBuffer.data(soff + 150 * idx + 132);

            auto s_yzzz_xyy = primBuffer.data(soff + 150 * idx + 133);

            auto s_yzzz_xyz = primBuffer.data(soff + 150 * idx + 134);

            auto s_yzzz_xzz = primBuffer.data(soff + 150 * idx + 135);

            auto s_yzzz_yyy = primBuffer.data(soff + 150 * idx + 136);

            auto s_yzzz_yyz = primBuffer.data(soff + 150 * idx + 137);

            auto s_yzzz_yzz = primBuffer.data(soff + 150 * idx + 138);

            auto s_yzzz_zzz = primBuffer.data(soff + 150 * idx + 139);

            auto s_zzzz_xxx = primBuffer.data(soff + 150 * idx + 140);

            auto s_zzzz_xxy = primBuffer.data(soff + 150 * idx + 141);

            auto s_zzzz_xxz = primBuffer.data(soff + 150 * idx + 142);

            auto s_zzzz_xyy = primBuffer.data(soff + 150 * idx + 143);

            auto s_zzzz_xyz = primBuffer.data(soff + 150 * idx + 144);

            auto s_zzzz_xzz = primBuffer.data(soff + 150 * idx + 145);

            auto s_zzzz_yyy = primBuffer.data(soff + 150 * idx + 146);

            auto s_zzzz_yyz = primBuffer.data(soff + 150 * idx + 147);

            auto s_zzzz_yzz = primBuffer.data(soff + 150 * idx + 148);

            auto s_zzzz_zzz = primBuffer.data(soff + 150 * idx + 149);

            // set up pointers to (G|T|F) integrals

            auto t_xxxx_xxx = primBuffer.data(koff + 150 * idx);

            auto t_xxxx_xxy = primBuffer.data(koff + 150 * idx + 1);

            auto t_xxxx_xxz = primBuffer.data(koff + 150 * idx + 2);

            auto t_xxxx_xyy = primBuffer.data(koff + 150 * idx + 3);

            auto t_xxxx_xyz = primBuffer.data(koff + 150 * idx + 4);

            auto t_xxxx_xzz = primBuffer.data(koff + 150 * idx + 5);

            auto t_xxxx_yyy = primBuffer.data(koff + 150 * idx + 6);

            auto t_xxxx_yyz = primBuffer.data(koff + 150 * idx + 7);

            auto t_xxxx_yzz = primBuffer.data(koff + 150 * idx + 8);

            auto t_xxxx_zzz = primBuffer.data(koff + 150 * idx + 9);

            auto t_xxxy_xxx = primBuffer.data(koff + 150 * idx + 10);

            auto t_xxxy_xxy = primBuffer.data(koff + 150 * idx + 11);

            auto t_xxxy_xxz = primBuffer.data(koff + 150 * idx + 12);

            auto t_xxxy_xyy = primBuffer.data(koff + 150 * idx + 13);

            auto t_xxxy_xyz = primBuffer.data(koff + 150 * idx + 14);

            auto t_xxxy_xzz = primBuffer.data(koff + 150 * idx + 15);

            auto t_xxxy_yyy = primBuffer.data(koff + 150 * idx + 16);

            auto t_xxxy_yyz = primBuffer.data(koff + 150 * idx + 17);

            auto t_xxxy_yzz = primBuffer.data(koff + 150 * idx + 18);

            auto t_xxxy_zzz = primBuffer.data(koff + 150 * idx + 19);

            auto t_xxxz_xxx = primBuffer.data(koff + 150 * idx + 20);

            auto t_xxxz_xxy = primBuffer.data(koff + 150 * idx + 21);

            auto t_xxxz_xxz = primBuffer.data(koff + 150 * idx + 22);

            auto t_xxxz_xyy = primBuffer.data(koff + 150 * idx + 23);

            auto t_xxxz_xyz = primBuffer.data(koff + 150 * idx + 24);

            auto t_xxxz_xzz = primBuffer.data(koff + 150 * idx + 25);

            auto t_xxxz_yyy = primBuffer.data(koff + 150 * idx + 26);

            auto t_xxxz_yyz = primBuffer.data(koff + 150 * idx + 27);

            auto t_xxxz_yzz = primBuffer.data(koff + 150 * idx + 28);

            auto t_xxxz_zzz = primBuffer.data(koff + 150 * idx + 29);

            auto t_xxyy_xxx = primBuffer.data(koff + 150 * idx + 30);

            auto t_xxyy_xxy = primBuffer.data(koff + 150 * idx + 31);

            auto t_xxyy_xxz = primBuffer.data(koff + 150 * idx + 32);

            auto t_xxyy_xyy = primBuffer.data(koff + 150 * idx + 33);

            auto t_xxyy_xyz = primBuffer.data(koff + 150 * idx + 34);

            auto t_xxyy_xzz = primBuffer.data(koff + 150 * idx + 35);

            auto t_xxyy_yyy = primBuffer.data(koff + 150 * idx + 36);

            auto t_xxyy_yyz = primBuffer.data(koff + 150 * idx + 37);

            auto t_xxyy_yzz = primBuffer.data(koff + 150 * idx + 38);

            auto t_xxyy_zzz = primBuffer.data(koff + 150 * idx + 39);

            auto t_xxyz_xxx = primBuffer.data(koff + 150 * idx + 40);

            auto t_xxyz_xxy = primBuffer.data(koff + 150 * idx + 41);

            auto t_xxyz_xxz = primBuffer.data(koff + 150 * idx + 42);

            auto t_xxyz_xyy = primBuffer.data(koff + 150 * idx + 43);

            auto t_xxyz_xyz = primBuffer.data(koff + 150 * idx + 44);

            auto t_xxyz_xzz = primBuffer.data(koff + 150 * idx + 45);

            auto t_xxyz_yyy = primBuffer.data(koff + 150 * idx + 46);

            auto t_xxyz_yyz = primBuffer.data(koff + 150 * idx + 47);

            auto t_xxyz_yzz = primBuffer.data(koff + 150 * idx + 48);

            auto t_xxyz_zzz = primBuffer.data(koff + 150 * idx + 49);

            auto t_xxzz_xxx = primBuffer.data(koff + 150 * idx + 50);

            auto t_xxzz_xxy = primBuffer.data(koff + 150 * idx + 51);

            auto t_xxzz_xxz = primBuffer.data(koff + 150 * idx + 52);

            auto t_xxzz_xyy = primBuffer.data(koff + 150 * idx + 53);

            auto t_xxzz_xyz = primBuffer.data(koff + 150 * idx + 54);

            auto t_xxzz_xzz = primBuffer.data(koff + 150 * idx + 55);

            auto t_xxzz_yyy = primBuffer.data(koff + 150 * idx + 56);

            auto t_xxzz_yyz = primBuffer.data(koff + 150 * idx + 57);

            auto t_xxzz_yzz = primBuffer.data(koff + 150 * idx + 58);

            auto t_xxzz_zzz = primBuffer.data(koff + 150 * idx + 59);

            auto t_xyyy_xxx = primBuffer.data(koff + 150 * idx + 60);

            auto t_xyyy_xxy = primBuffer.data(koff + 150 * idx + 61);

            auto t_xyyy_xxz = primBuffer.data(koff + 150 * idx + 62);

            auto t_xyyy_xyy = primBuffer.data(koff + 150 * idx + 63);

            auto t_xyyy_xyz = primBuffer.data(koff + 150 * idx + 64);

            auto t_xyyy_xzz = primBuffer.data(koff + 150 * idx + 65);

            auto t_xyyy_yyy = primBuffer.data(koff + 150 * idx + 66);

            auto t_xyyy_yyz = primBuffer.data(koff + 150 * idx + 67);

            auto t_xyyy_yzz = primBuffer.data(koff + 150 * idx + 68);

            auto t_xyyy_zzz = primBuffer.data(koff + 150 * idx + 69);

            auto t_xyyz_xxx = primBuffer.data(koff + 150 * idx + 70);

            auto t_xyyz_xxy = primBuffer.data(koff + 150 * idx + 71);

            auto t_xyyz_xxz = primBuffer.data(koff + 150 * idx + 72);

            auto t_xyyz_xyy = primBuffer.data(koff + 150 * idx + 73);

            auto t_xyyz_xyz = primBuffer.data(koff + 150 * idx + 74);

            auto t_xyyz_xzz = primBuffer.data(koff + 150 * idx + 75);

            auto t_xyyz_yyy = primBuffer.data(koff + 150 * idx + 76);

            auto t_xyyz_yyz = primBuffer.data(koff + 150 * idx + 77);

            auto t_xyyz_yzz = primBuffer.data(koff + 150 * idx + 78);

            auto t_xyyz_zzz = primBuffer.data(koff + 150 * idx + 79);

            auto t_xyzz_xxx = primBuffer.data(koff + 150 * idx + 80);

            auto t_xyzz_xxy = primBuffer.data(koff + 150 * idx + 81);

            auto t_xyzz_xxz = primBuffer.data(koff + 150 * idx + 82);

            auto t_xyzz_xyy = primBuffer.data(koff + 150 * idx + 83);

            auto t_xyzz_xyz = primBuffer.data(koff + 150 * idx + 84);

            auto t_xyzz_xzz = primBuffer.data(koff + 150 * idx + 85);

            auto t_xyzz_yyy = primBuffer.data(koff + 150 * idx + 86);

            auto t_xyzz_yyz = primBuffer.data(koff + 150 * idx + 87);

            auto t_xyzz_yzz = primBuffer.data(koff + 150 * idx + 88);

            auto t_xyzz_zzz = primBuffer.data(koff + 150 * idx + 89);

            auto t_xzzz_xxx = primBuffer.data(koff + 150 * idx + 90);

            auto t_xzzz_xxy = primBuffer.data(koff + 150 * idx + 91);

            auto t_xzzz_xxz = primBuffer.data(koff + 150 * idx + 92);

            auto t_xzzz_xyy = primBuffer.data(koff + 150 * idx + 93);

            auto t_xzzz_xyz = primBuffer.data(koff + 150 * idx + 94);

            auto t_xzzz_xzz = primBuffer.data(koff + 150 * idx + 95);

            auto t_xzzz_yyy = primBuffer.data(koff + 150 * idx + 96);

            auto t_xzzz_yyz = primBuffer.data(koff + 150 * idx + 97);

            auto t_xzzz_yzz = primBuffer.data(koff + 150 * idx + 98);

            auto t_xzzz_zzz = primBuffer.data(koff + 150 * idx + 99);

            auto t_yyyy_xxx = primBuffer.data(koff + 150 * idx + 100);

            auto t_yyyy_xxy = primBuffer.data(koff + 150 * idx + 101);

            auto t_yyyy_xxz = primBuffer.data(koff + 150 * idx + 102);

            auto t_yyyy_xyy = primBuffer.data(koff + 150 * idx + 103);

            auto t_yyyy_xyz = primBuffer.data(koff + 150 * idx + 104);

            auto t_yyyy_xzz = primBuffer.data(koff + 150 * idx + 105);

            auto t_yyyy_yyy = primBuffer.data(koff + 150 * idx + 106);

            auto t_yyyy_yyz = primBuffer.data(koff + 150 * idx + 107);

            auto t_yyyy_yzz = primBuffer.data(koff + 150 * idx + 108);

            auto t_yyyy_zzz = primBuffer.data(koff + 150 * idx + 109);

            auto t_yyyz_xxx = primBuffer.data(koff + 150 * idx + 110);

            auto t_yyyz_xxy = primBuffer.data(koff + 150 * idx + 111);

            auto t_yyyz_xxz = primBuffer.data(koff + 150 * idx + 112);

            auto t_yyyz_xyy = primBuffer.data(koff + 150 * idx + 113);

            auto t_yyyz_xyz = primBuffer.data(koff + 150 * idx + 114);

            auto t_yyyz_xzz = primBuffer.data(koff + 150 * idx + 115);

            auto t_yyyz_yyy = primBuffer.data(koff + 150 * idx + 116);

            auto t_yyyz_yyz = primBuffer.data(koff + 150 * idx + 117);

            auto t_yyyz_yzz = primBuffer.data(koff + 150 * idx + 118);

            auto t_yyyz_zzz = primBuffer.data(koff + 150 * idx + 119);

            auto t_yyzz_xxx = primBuffer.data(koff + 150 * idx + 120);

            auto t_yyzz_xxy = primBuffer.data(koff + 150 * idx + 121);

            auto t_yyzz_xxz = primBuffer.data(koff + 150 * idx + 122);

            auto t_yyzz_xyy = primBuffer.data(koff + 150 * idx + 123);

            auto t_yyzz_xyz = primBuffer.data(koff + 150 * idx + 124);

            auto t_yyzz_xzz = primBuffer.data(koff + 150 * idx + 125);

            auto t_yyzz_yyy = primBuffer.data(koff + 150 * idx + 126);

            auto t_yyzz_yyz = primBuffer.data(koff + 150 * idx + 127);

            auto t_yyzz_yzz = primBuffer.data(koff + 150 * idx + 128);

            auto t_yyzz_zzz = primBuffer.data(koff + 150 * idx + 129);

            auto t_yzzz_xxx = primBuffer.data(koff + 150 * idx + 130);

            auto t_yzzz_xxy = primBuffer.data(koff + 150 * idx + 131);

            auto t_yzzz_xxz = primBuffer.data(koff + 150 * idx + 132);

            auto t_yzzz_xyy = primBuffer.data(koff + 150 * idx + 133);

            auto t_yzzz_xyz = primBuffer.data(koff + 150 * idx + 134);

            auto t_yzzz_xzz = primBuffer.data(koff + 150 * idx + 135);

            auto t_yzzz_yyy = primBuffer.data(koff + 150 * idx + 136);

            auto t_yzzz_yyz = primBuffer.data(koff + 150 * idx + 137);

            auto t_yzzz_yzz = primBuffer.data(koff + 150 * idx + 138);

            auto t_yzzz_zzz = primBuffer.data(koff + 150 * idx + 139);

            auto t_zzzz_xxx = primBuffer.data(koff + 150 * idx + 140);

            auto t_zzzz_xxy = primBuffer.data(koff + 150 * idx + 141);

            auto t_zzzz_xxz = primBuffer.data(koff + 150 * idx + 142);

            auto t_zzzz_xyy = primBuffer.data(koff + 150 * idx + 143);

            auto t_zzzz_xyz = primBuffer.data(koff + 150 * idx + 144);

            auto t_zzzz_xzz = primBuffer.data(koff + 150 * idx + 145);

            auto t_zzzz_yyy = primBuffer.data(koff + 150 * idx + 146);

            auto t_zzzz_yyz = primBuffer.data(koff + 150 * idx + 147);

            auto t_zzzz_yzz = primBuffer.data(koff + 150 * idx + 148);

            auto t_zzzz_zzz = primBuffer.data(koff + 150 * idx + 149);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xxx_xx, s_xxx_xy,\
                                     s_xxx_xz, s_xxx_yy, s_xxx_yz, s_xxx_zz, s_xxy_xx,\
                                     s_xxy_xy, s_xxy_xz, s_xxy_yy, s_xxy_yz, s_xxy_zz,\
                                     s_xxz_xx, s_xxz_xy, s_xxz_xz, s_xxz_yy, s_xxz_yz,\
                                     s_xxz_zz, s_xyy_xx, s_xyy_xy, s_xyy_xz, s_xyy_yy,\
                                     s_xyy_yz, s_xyy_zz, s_xyz_xx, s_xyz_xy, s_xyz_xz,\
                                     s_xyz_yy, s_xyz_yz, s_xyz_zz, s_xzz_xx, s_xzz_xy,\
                                     s_xzz_xz, s_xzz_yy, s_xzz_yz, s_xzz_zz, s_yyy_xx,\
                                     s_yyy_xy, s_yyy_xz, s_yyy_yy, s_yyy_yz, s_yyy_zz,\
                                     s_yyz_xx, s_yyz_xy, s_yyz_xz, s_yyz_yy, s_yyz_yz,\
                                     s_yyz_zz, s_yzz_xx, s_yzz_xy, s_yzz_xz, s_yzz_yy,\
                                     s_yzz_yz, s_yzz_zz, s_zzz_xx, s_zzz_xy, s_zzz_xz,\
                                     s_zzz_yy, s_zzz_yz, s_zzz_zz, t_xxx_xx, t_xxx_xy,\
                                     t_xxx_xz, t_xxx_yy, t_xxx_yz, t_xxx_zz, t_xxy_xx,\
                                     t_xxy_xy, t_xxy_xz, t_xxy_yy, t_xxy_yz, t_xxy_zz,\
                                     t_xxz_xx, t_xxz_xy, t_xxz_xz, t_xxz_yy, t_xxz_yz,\
                                     t_xxz_zz, t_xyy_xx, t_xyy_xy, t_xyy_xz, t_xyy_yy,\
                                     t_xyy_yz, t_xyy_zz, t_xyz_xx, t_xyz_xy, t_xyz_xz,\
                                     t_xyz_yy, t_xyz_yz, t_xyz_zz, t_xzz_xx, t_xzz_xy,\
                                     t_xzz_xz, t_xzz_yy, t_xzz_yz, t_xzz_zz, t_yyy_xx,\
                                     t_yyy_xy, t_yyy_xz, t_yyy_yy, t_yyy_yz, t_yyy_zz,\
                                     t_yyz_xx, t_yyz_xy, t_yyz_xz, t_yyz_yy, t_yyz_yz,\
                                     t_yyz_zz, t_yzz_xx, t_yzz_xy, t_yzz_xz, t_yzz_yy,\
                                     t_yzz_yz, t_yzz_zz, t_zzz_xx, t_zzz_xy, t_zzz_xz,\
                                     t_zzz_yy, t_zzz_yz, t_zzz_zz, s_xx_xxx, s_xx_xxy,\
                                     s_xx_xxz, s_xx_xyy, s_xx_xyz, s_xx_xzz, s_xx_yyy,\
                                     s_xx_yyz, s_xx_yzz, s_xx_zzz, s_xy_xxx, s_xy_xxy,\
                                     s_xy_xxz, s_xy_xyy, s_xy_xyz, s_xy_xzz, s_xy_yyy,\
                                     s_xy_yyz, s_xy_yzz, s_xy_zzz, s_xz_xxx, s_xz_xxy,\
                                     s_xz_xxz, s_xz_xyy, s_xz_xyz, s_xz_xzz, s_xz_yyy,\
                                     s_xz_yyz, s_xz_yzz, s_xz_zzz, s_yy_xxx, s_yy_xxy,\
                                     s_yy_xxz, s_yy_xyy, s_yy_xyz, s_yy_xzz, s_yy_yyy,\
                                     s_yy_yyz, s_yy_yzz, s_yy_zzz, s_yz_xxx, s_yz_xxy,\
                                     s_yz_xxz, s_yz_xyy, s_yz_xyz, s_yz_xzz, s_yz_yyy,\
                                     s_yz_yyz, s_yz_yzz, s_yz_zzz, s_zz_xxx, s_zz_xxy,\
                                     s_zz_xxz, s_zz_xyy, s_zz_xyz, s_zz_xzz, s_zz_yyy,\
                                     s_zz_yyz, s_zz_yzz, s_zz_zzz, t_xx_xxx, t_xx_xxy,\
                                     t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, t_xx_yyy,\
                                     t_xx_yyz, t_xx_yzz, t_xx_zzz, t_xy_xxx, t_xy_xxy,\
                                     t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy,\
                                     t_xy_yyz, t_xy_yzz, t_xy_zzz, t_xz_xxx, t_xz_xxy,\
                                     t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy,\
                                     t_xz_yyz, t_xz_yzz, t_xz_zzz, t_yy_xxx, t_yy_xxy,\
                                     t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, t_yy_yyy,\
                                     t_yy_yyz, t_yy_yzz, t_yy_zzz, t_yz_xxx, t_yz_xxy,\
                                     t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy,\
                                     t_yz_yyz, t_yz_yzz, t_yz_zzz, t_zz_xxx, t_zz_xxy,\
                                     t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, t_zz_yyy,\
                                     t_zz_yyz, t_zz_yzz, t_zz_zzz, s_xxx_xxx, s_xxx_xxy,\
                                     s_xxx_xxz, s_xxx_xyy, s_xxx_xyz, s_xxx_xzz,\
                                     s_xxx_yyy, s_xxx_yyz, s_xxx_yzz, s_xxx_zzz,\
                                     s_xxy_xxx, s_xxy_xxy, s_xxy_xxz, s_xxy_xyy,\
                                     s_xxy_xyz, s_xxy_xzz, s_xxy_yyy, s_xxy_yyz,\
                                     s_xxy_yzz, s_xxy_zzz, s_xxz_xxx, s_xxz_xxy,\
                                     s_xxz_xxz, s_xxz_xyy, s_xxz_xyz, s_xxz_xzz,\
                                     s_xxz_yyy, s_xxz_yyz, s_xxz_yzz, s_xxz_zzz,\
                                     s_xyy_xxx, s_xyy_xxy, s_xyy_xxz, s_xyy_xyy,\
                                     s_xyy_xyz, s_xyy_xzz, s_xyy_yyy, s_xyy_yyz,\
                                     s_xyy_yzz, s_xyy_zzz, s_xyz_xxx, s_xyz_xxy,\
                                     s_xyz_xxz, s_xyz_xyy, s_xyz_xyz, s_xyz_xzz,\
                                     s_xyz_yyy, s_xyz_yyz, s_xyz_yzz, s_xyz_zzz,\
                                     s_xzz_xxx, s_xzz_xxy, s_xzz_xxz, s_xzz_xyy,\
                                     s_xzz_xyz, s_xzz_xzz, s_xzz_yyy, s_xzz_yyz,\
                                     s_xzz_yzz, s_xzz_zzz, s_yyy_xxx, s_yyy_xxy,\
                                     s_yyy_xxz, s_yyy_xyy, s_yyy_xyz, s_yyy_xzz,\
                                     s_yyy_yyy, s_yyy_yyz, s_yyy_yzz, s_yyy_zzz,\
                                     s_yyz_xxx, s_yyz_xxy, s_yyz_xxz, s_yyz_xyy,\
                                     s_yyz_xyz, s_yyz_xzz, s_yyz_yyy, s_yyz_yyz,\
                                     s_yyz_yzz, s_yyz_zzz, s_yzz_xxx, s_yzz_xxy,\
                                     s_yzz_xxz, s_yzz_xyy, s_yzz_xyz, s_yzz_xzz,\
                                     s_yzz_yyy, s_yzz_yyz, s_yzz_yzz, s_yzz_zzz,\
                                     s_zzz_xxx, s_zzz_xxy, s_zzz_xxz, s_zzz_xyy,\
                                     s_zzz_xyz, s_zzz_xzz, s_zzz_yyy, s_zzz_yyz,\
                                     s_zzz_yzz, s_zzz_zzz, t_xxx_xxx, t_xxx_xxy,\
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
                                     t_zzz_yzz, t_zzz_zzz, s_xxxx_xxx, s_xxxx_xxy,\
                                     s_xxxx_xxz, s_xxxx_xyy, s_xxxx_xyz, s_xxxx_xzz,\
                                     s_xxxx_yyy, s_xxxx_yyz, s_xxxx_yzz, s_xxxx_zzz,\
                                     s_xxxy_xxx, s_xxxy_xxy, s_xxxy_xxz, s_xxxy_xyy,\
                                     s_xxxy_xyz, s_xxxy_xzz, s_xxxy_yyy, s_xxxy_yyz,\
                                     s_xxxy_yzz, s_xxxy_zzz, s_xxxz_xxx, s_xxxz_xxy,\
                                     s_xxxz_xxz, s_xxxz_xyy, s_xxxz_xyz, s_xxxz_xzz,\
                                     s_xxxz_yyy, s_xxxz_yyz, s_xxxz_yzz, s_xxxz_zzz,\
                                     s_xxyy_xxx, s_xxyy_xxy, s_xxyy_xxz, s_xxyy_xyy,\
                                     s_xxyy_xyz, s_xxyy_xzz, s_xxyy_yyy, s_xxyy_yyz,\
                                     s_xxyy_yzz, s_xxyy_zzz, s_xxyz_xxx, s_xxyz_xxy,\
                                     s_xxyz_xxz, s_xxyz_xyy, s_xxyz_xyz, s_xxyz_xzz,\
                                     s_xxyz_yyy, s_xxyz_yyz, s_xxyz_yzz, s_xxyz_zzz,\
                                     s_xxzz_xxx, s_xxzz_xxy, s_xxzz_xxz, s_xxzz_xyy,\
                                     s_xxzz_xyz, s_xxzz_xzz, s_xxzz_yyy, s_xxzz_yyz,\
                                     s_xxzz_yzz, s_xxzz_zzz, s_xyyy_xxx, s_xyyy_xxy,\
                                     s_xyyy_xxz, s_xyyy_xyy, s_xyyy_xyz, s_xyyy_xzz,\
                                     s_xyyy_yyy, s_xyyy_yyz, s_xyyy_yzz, s_xyyy_zzz,\
                                     s_xyyz_xxx, s_xyyz_xxy, s_xyyz_xxz, s_xyyz_xyy,\
                                     s_xyyz_xyz, s_xyyz_xzz, s_xyyz_yyy, s_xyyz_yyz,\
                                     s_xyyz_yzz, s_xyyz_zzz, s_xyzz_xxx, s_xyzz_xxy,\
                                     s_xyzz_xxz, s_xyzz_xyy, s_xyzz_xyz, s_xyzz_xzz,\
                                     s_xyzz_yyy, s_xyzz_yyz, s_xyzz_yzz, s_xyzz_zzz,\
                                     s_xzzz_xxx, s_xzzz_xxy, s_xzzz_xxz, s_xzzz_xyy,\
                                     s_xzzz_xyz, s_xzzz_xzz, s_xzzz_yyy, s_xzzz_yyz,\
                                     s_xzzz_yzz, s_xzzz_zzz, s_yyyy_xxx, s_yyyy_xxy,\
                                     s_yyyy_xxz, s_yyyy_xyy, s_yyyy_xyz, s_yyyy_xzz,\
                                     s_yyyy_yyy, s_yyyy_yyz, s_yyyy_yzz, s_yyyy_zzz,\
                                     s_yyyz_xxx, s_yyyz_xxy, s_yyyz_xxz, s_yyyz_xyy,\
                                     s_yyyz_xyz, s_yyyz_xzz, s_yyyz_yyy, s_yyyz_yyz,\
                                     s_yyyz_yzz, s_yyyz_zzz, s_yyzz_xxx, s_yyzz_xxy,\
                                     s_yyzz_xxz, s_yyzz_xyy, s_yyzz_xyz, s_yyzz_xzz,\
                                     s_yyzz_yyy, s_yyzz_yyz, s_yyzz_yzz, s_yyzz_zzz,\
                                     s_yzzz_xxx, s_yzzz_xxy, s_yzzz_xxz, s_yzzz_xyy,\
                                     s_yzzz_xyz, s_yzzz_xzz, s_yzzz_yyy, s_yzzz_yyz,\
                                     s_yzzz_yzz, s_yzzz_zzz, s_zzzz_xxx, s_zzzz_xxy,\
                                     s_zzzz_xxz, s_zzzz_xyy, s_zzzz_xyz, s_zzzz_xzz,\
                                     s_zzzz_yyy, s_zzzz_yyz, s_zzzz_yzz, s_zzzz_zzz,\
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
                                     t_zzzz_yzz, t_zzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxxx_xxx[j] = fr * s_xxx_xxx[j] + f2t * (3.0 * s_xx_xxx[j] + 3.0 * s_xxx_xx[j]);

                t_xxxx_xxx[j] = fr * t_xxx_xxx[j] + f2t * (3.0 * t_xx_xxx[j] + 3.0 * t_xxx_xx[j]) + f2z * (s_xxxx_xxx[j] - f2a * 3.0 * s_xx_xxx[j]);

                s_xxxx_xxy[j] = fr * s_xxx_xxy[j] + f2t * (3.0 * s_xx_xxy[j] + 2.0 * s_xxx_xy[j]);

                t_xxxx_xxy[j] = fr * t_xxx_xxy[j] + f2t * (3.0 * t_xx_xxy[j] + 2.0 * t_xxx_xy[j]) + f2z * (s_xxxx_xxy[j] - f2a * 3.0 * s_xx_xxy[j]);

                s_xxxx_xxz[j] = fr * s_xxx_xxz[j] + f2t * (3.0 * s_xx_xxz[j] + 2.0 * s_xxx_xz[j]);

                t_xxxx_xxz[j] = fr * t_xxx_xxz[j] + f2t * (3.0 * t_xx_xxz[j] + 2.0 * t_xxx_xz[j]) + f2z * (s_xxxx_xxz[j] - f2a * 3.0 * s_xx_xxz[j]);

                s_xxxx_xyy[j] = fr * s_xxx_xyy[j] + f2t * (3.0 * s_xx_xyy[j] + s_xxx_yy[j]);

                t_xxxx_xyy[j] = fr * t_xxx_xyy[j] + f2t * (3.0 * t_xx_xyy[j] + t_xxx_yy[j]) + f2z * (s_xxxx_xyy[j] - f2a * 3.0 * s_xx_xyy[j]);

                s_xxxx_xyz[j] = fr * s_xxx_xyz[j] + f2t * (3.0 * s_xx_xyz[j] + s_xxx_yz[j]);

                t_xxxx_xyz[j] = fr * t_xxx_xyz[j] + f2t * (3.0 * t_xx_xyz[j] + t_xxx_yz[j]) + f2z * (s_xxxx_xyz[j] - f2a * 3.0 * s_xx_xyz[j]);

                s_xxxx_xzz[j] = fr * s_xxx_xzz[j] + f2t * (3.0 * s_xx_xzz[j] + s_xxx_zz[j]);

                t_xxxx_xzz[j] = fr * t_xxx_xzz[j] + f2t * (3.0 * t_xx_xzz[j] + t_xxx_zz[j]) + f2z * (s_xxxx_xzz[j] - f2a * 3.0 * s_xx_xzz[j]);

                s_xxxx_yyy[j] = fr * s_xxx_yyy[j] + f2t * 3.0 * s_xx_yyy[j];

                t_xxxx_yyy[j] = fr * t_xxx_yyy[j] + f2t * 3.0 * t_xx_yyy[j] + f2z * (s_xxxx_yyy[j] - f2a * 3.0 * s_xx_yyy[j]);

                s_xxxx_yyz[j] = fr * s_xxx_yyz[j] + f2t * 3.0 * s_xx_yyz[j];

                t_xxxx_yyz[j] = fr * t_xxx_yyz[j] + f2t * 3.0 * t_xx_yyz[j] + f2z * (s_xxxx_yyz[j] - f2a * 3.0 * s_xx_yyz[j]);

                s_xxxx_yzz[j] = fr * s_xxx_yzz[j] + f2t * 3.0 * s_xx_yzz[j];

                t_xxxx_yzz[j] = fr * t_xxx_yzz[j] + f2t * 3.0 * t_xx_yzz[j] + f2z * (s_xxxx_yzz[j] - f2a * 3.0 * s_xx_yzz[j]);

                s_xxxx_zzz[j] = fr * s_xxx_zzz[j] + f2t * 3.0 * s_xx_zzz[j];

                t_xxxx_zzz[j] = fr * t_xxx_zzz[j] + f2t * 3.0 * t_xx_zzz[j] + f2z * (s_xxxx_zzz[j] - f2a * 3.0 * s_xx_zzz[j]);

                s_xxxy_xxx[j] = fr * s_xxy_xxx[j] + f2t * (2.0 * s_xy_xxx[j] + 3.0 * s_xxy_xx[j]);

                t_xxxy_xxx[j] = fr * t_xxy_xxx[j] + f2t * (2.0 * t_xy_xxx[j] + 3.0 * t_xxy_xx[j]) + f2z * (s_xxxy_xxx[j] - f2a * 2.0 * s_xy_xxx[j]);

                s_xxxy_xxy[j] = fr * s_xxy_xxy[j] + f2t * (2.0 * s_xy_xxy[j] + 2.0 * s_xxy_xy[j]);

                t_xxxy_xxy[j] = fr * t_xxy_xxy[j] + f2t * (2.0 * t_xy_xxy[j] + 2.0 * t_xxy_xy[j]) + f2z * (s_xxxy_xxy[j] - f2a * 2.0 * s_xy_xxy[j]);

                s_xxxy_xxz[j] = fr * s_xxy_xxz[j] + f2t * (2.0 * s_xy_xxz[j] + 2.0 * s_xxy_xz[j]);

                t_xxxy_xxz[j] = fr * t_xxy_xxz[j] + f2t * (2.0 * t_xy_xxz[j] + 2.0 * t_xxy_xz[j]) + f2z * (s_xxxy_xxz[j] - f2a * 2.0 * s_xy_xxz[j]);

                s_xxxy_xyy[j] = fr * s_xxy_xyy[j] + f2t * (2.0 * s_xy_xyy[j] + s_xxy_yy[j]);

                t_xxxy_xyy[j] = fr * t_xxy_xyy[j] + f2t * (2.0 * t_xy_xyy[j] + t_xxy_yy[j]) + f2z * (s_xxxy_xyy[j] - f2a * 2.0 * s_xy_xyy[j]);

                s_xxxy_xyz[j] = fr * s_xxy_xyz[j] + f2t * (2.0 * s_xy_xyz[j] + s_xxy_yz[j]);

                t_xxxy_xyz[j] = fr * t_xxy_xyz[j] + f2t * (2.0 * t_xy_xyz[j] + t_xxy_yz[j]) + f2z * (s_xxxy_xyz[j] - f2a * 2.0 * s_xy_xyz[j]);

                s_xxxy_xzz[j] = fr * s_xxy_xzz[j] + f2t * (2.0 * s_xy_xzz[j] + s_xxy_zz[j]);

                t_xxxy_xzz[j] = fr * t_xxy_xzz[j] + f2t * (2.0 * t_xy_xzz[j] + t_xxy_zz[j]) + f2z * (s_xxxy_xzz[j] - f2a * 2.0 * s_xy_xzz[j]);

                s_xxxy_yyy[j] = fr * s_xxy_yyy[j] + f2t * 2.0 * s_xy_yyy[j];

                t_xxxy_yyy[j] = fr * t_xxy_yyy[j] + f2t * 2.0 * t_xy_yyy[j] + f2z * (s_xxxy_yyy[j] - f2a * 2.0 * s_xy_yyy[j]);

                s_xxxy_yyz[j] = fr * s_xxy_yyz[j] + f2t * 2.0 * s_xy_yyz[j];

                t_xxxy_yyz[j] = fr * t_xxy_yyz[j] + f2t * 2.0 * t_xy_yyz[j] + f2z * (s_xxxy_yyz[j] - f2a * 2.0 * s_xy_yyz[j]);

                s_xxxy_yzz[j] = fr * s_xxy_yzz[j] + f2t * 2.0 * s_xy_yzz[j];

                t_xxxy_yzz[j] = fr * t_xxy_yzz[j] + f2t * 2.0 * t_xy_yzz[j] + f2z * (s_xxxy_yzz[j] - f2a * 2.0 * s_xy_yzz[j]);

                s_xxxy_zzz[j] = fr * s_xxy_zzz[j] + f2t * 2.0 * s_xy_zzz[j];

                t_xxxy_zzz[j] = fr * t_xxy_zzz[j] + f2t * 2.0 * t_xy_zzz[j] + f2z * (s_xxxy_zzz[j] - f2a * 2.0 * s_xy_zzz[j]);

                s_xxxz_xxx[j] = fr * s_xxz_xxx[j] + f2t * (2.0 * s_xz_xxx[j] + 3.0 * s_xxz_xx[j]);

                t_xxxz_xxx[j] = fr * t_xxz_xxx[j] + f2t * (2.0 * t_xz_xxx[j] + 3.0 * t_xxz_xx[j]) + f2z * (s_xxxz_xxx[j] - f2a * 2.0 * s_xz_xxx[j]);

                s_xxxz_xxy[j] = fr * s_xxz_xxy[j] + f2t * (2.0 * s_xz_xxy[j] + 2.0 * s_xxz_xy[j]);

                t_xxxz_xxy[j] = fr * t_xxz_xxy[j] + f2t * (2.0 * t_xz_xxy[j] + 2.0 * t_xxz_xy[j]) + f2z * (s_xxxz_xxy[j] - f2a * 2.0 * s_xz_xxy[j]);

                s_xxxz_xxz[j] = fr * s_xxz_xxz[j] + f2t * (2.0 * s_xz_xxz[j] + 2.0 * s_xxz_xz[j]);

                t_xxxz_xxz[j] = fr * t_xxz_xxz[j] + f2t * (2.0 * t_xz_xxz[j] + 2.0 * t_xxz_xz[j]) + f2z * (s_xxxz_xxz[j] - f2a * 2.0 * s_xz_xxz[j]);

                s_xxxz_xyy[j] = fr * s_xxz_xyy[j] + f2t * (2.0 * s_xz_xyy[j] + s_xxz_yy[j]);

                t_xxxz_xyy[j] = fr * t_xxz_xyy[j] + f2t * (2.0 * t_xz_xyy[j] + t_xxz_yy[j]) + f2z * (s_xxxz_xyy[j] - f2a * 2.0 * s_xz_xyy[j]);

                s_xxxz_xyz[j] = fr * s_xxz_xyz[j] + f2t * (2.0 * s_xz_xyz[j] + s_xxz_yz[j]);

                t_xxxz_xyz[j] = fr * t_xxz_xyz[j] + f2t * (2.0 * t_xz_xyz[j] + t_xxz_yz[j]) + f2z * (s_xxxz_xyz[j] - f2a * 2.0 * s_xz_xyz[j]);

                s_xxxz_xzz[j] = fr * s_xxz_xzz[j] + f2t * (2.0 * s_xz_xzz[j] + s_xxz_zz[j]);

                t_xxxz_xzz[j] = fr * t_xxz_xzz[j] + f2t * (2.0 * t_xz_xzz[j] + t_xxz_zz[j]) + f2z * (s_xxxz_xzz[j] - f2a * 2.0 * s_xz_xzz[j]);

                s_xxxz_yyy[j] = fr * s_xxz_yyy[j] + f2t * 2.0 * s_xz_yyy[j];

                t_xxxz_yyy[j] = fr * t_xxz_yyy[j] + f2t * 2.0 * t_xz_yyy[j] + f2z * (s_xxxz_yyy[j] - f2a * 2.0 * s_xz_yyy[j]);

                s_xxxz_yyz[j] = fr * s_xxz_yyz[j] + f2t * 2.0 * s_xz_yyz[j];

                t_xxxz_yyz[j] = fr * t_xxz_yyz[j] + f2t * 2.0 * t_xz_yyz[j] + f2z * (s_xxxz_yyz[j] - f2a * 2.0 * s_xz_yyz[j]);

                s_xxxz_yzz[j] = fr * s_xxz_yzz[j] + f2t * 2.0 * s_xz_yzz[j];

                t_xxxz_yzz[j] = fr * t_xxz_yzz[j] + f2t * 2.0 * t_xz_yzz[j] + f2z * (s_xxxz_yzz[j] - f2a * 2.0 * s_xz_yzz[j]);

                s_xxxz_zzz[j] = fr * s_xxz_zzz[j] + f2t * 2.0 * s_xz_zzz[j];

                t_xxxz_zzz[j] = fr * t_xxz_zzz[j] + f2t * 2.0 * t_xz_zzz[j] + f2z * (s_xxxz_zzz[j] - f2a * 2.0 * s_xz_zzz[j]);

                s_xxyy_xxx[j] = fr * s_xyy_xxx[j] + f2t * (s_yy_xxx[j] + 3.0 * s_xyy_xx[j]);

                t_xxyy_xxx[j] = fr * t_xyy_xxx[j] + f2t * (t_yy_xxx[j] + 3.0 * t_xyy_xx[j]) + f2z * (s_xxyy_xxx[j] - f2a * s_yy_xxx[j]);

                s_xxyy_xxy[j] = fr * s_xyy_xxy[j] + f2t * (s_yy_xxy[j] + 2.0 * s_xyy_xy[j]);

                t_xxyy_xxy[j] = fr * t_xyy_xxy[j] + f2t * (t_yy_xxy[j] + 2.0 * t_xyy_xy[j]) + f2z * (s_xxyy_xxy[j] - f2a * s_yy_xxy[j]);

                s_xxyy_xxz[j] = fr * s_xyy_xxz[j] + f2t * (s_yy_xxz[j] + 2.0 * s_xyy_xz[j]);

                t_xxyy_xxz[j] = fr * t_xyy_xxz[j] + f2t * (t_yy_xxz[j] + 2.0 * t_xyy_xz[j]) + f2z * (s_xxyy_xxz[j] - f2a * s_yy_xxz[j]);

                s_xxyy_xyy[j] = fr * s_xyy_xyy[j] + f2t * (s_yy_xyy[j] + s_xyy_yy[j]);

                t_xxyy_xyy[j] = fr * t_xyy_xyy[j] + f2t * (t_yy_xyy[j] + t_xyy_yy[j]) + f2z * (s_xxyy_xyy[j] - f2a * s_yy_xyy[j]);

                s_xxyy_xyz[j] = fr * s_xyy_xyz[j] + f2t * (s_yy_xyz[j] + s_xyy_yz[j]);

                t_xxyy_xyz[j] = fr * t_xyy_xyz[j] + f2t * (t_yy_xyz[j] + t_xyy_yz[j]) + f2z * (s_xxyy_xyz[j] - f2a * s_yy_xyz[j]);

                s_xxyy_xzz[j] = fr * s_xyy_xzz[j] + f2t * (s_yy_xzz[j] + s_xyy_zz[j]);

                t_xxyy_xzz[j] = fr * t_xyy_xzz[j] + f2t * (t_yy_xzz[j] + t_xyy_zz[j]) + f2z * (s_xxyy_xzz[j] - f2a * s_yy_xzz[j]);

                s_xxyy_yyy[j] = fr * s_xyy_yyy[j] + f2t * s_yy_yyy[j];

                t_xxyy_yyy[j] = fr * t_xyy_yyy[j] + f2t * t_yy_yyy[j] + f2z * (s_xxyy_yyy[j] - f2a * s_yy_yyy[j]);

                s_xxyy_yyz[j] = fr * s_xyy_yyz[j] + f2t * s_yy_yyz[j];

                t_xxyy_yyz[j] = fr * t_xyy_yyz[j] + f2t * t_yy_yyz[j] + f2z * (s_xxyy_yyz[j] - f2a * s_yy_yyz[j]);

                s_xxyy_yzz[j] = fr * s_xyy_yzz[j] + f2t * s_yy_yzz[j];

                t_xxyy_yzz[j] = fr * t_xyy_yzz[j] + f2t * t_yy_yzz[j] + f2z * (s_xxyy_yzz[j] - f2a * s_yy_yzz[j]);

                s_xxyy_zzz[j] = fr * s_xyy_zzz[j] + f2t * s_yy_zzz[j];

                t_xxyy_zzz[j] = fr * t_xyy_zzz[j] + f2t * t_yy_zzz[j] + f2z * (s_xxyy_zzz[j] - f2a * s_yy_zzz[j]);

                s_xxyz_xxx[j] = fr * s_xyz_xxx[j] + f2t * (s_yz_xxx[j] + 3.0 * s_xyz_xx[j]);

                t_xxyz_xxx[j] = fr * t_xyz_xxx[j] + f2t * (t_yz_xxx[j] + 3.0 * t_xyz_xx[j]) + f2z * (s_xxyz_xxx[j] - f2a * s_yz_xxx[j]);

                s_xxyz_xxy[j] = fr * s_xyz_xxy[j] + f2t * (s_yz_xxy[j] + 2.0 * s_xyz_xy[j]);

                t_xxyz_xxy[j] = fr * t_xyz_xxy[j] + f2t * (t_yz_xxy[j] + 2.0 * t_xyz_xy[j]) + f2z * (s_xxyz_xxy[j] - f2a * s_yz_xxy[j]);

                s_xxyz_xxz[j] = fr * s_xyz_xxz[j] + f2t * (s_yz_xxz[j] + 2.0 * s_xyz_xz[j]);

                t_xxyz_xxz[j] = fr * t_xyz_xxz[j] + f2t * (t_yz_xxz[j] + 2.0 * t_xyz_xz[j]) + f2z * (s_xxyz_xxz[j] - f2a * s_yz_xxz[j]);

                s_xxyz_xyy[j] = fr * s_xyz_xyy[j] + f2t * (s_yz_xyy[j] + s_xyz_yy[j]);

                t_xxyz_xyy[j] = fr * t_xyz_xyy[j] + f2t * (t_yz_xyy[j] + t_xyz_yy[j]) + f2z * (s_xxyz_xyy[j] - f2a * s_yz_xyy[j]);

                s_xxyz_xyz[j] = fr * s_xyz_xyz[j] + f2t * (s_yz_xyz[j] + s_xyz_yz[j]);

                t_xxyz_xyz[j] = fr * t_xyz_xyz[j] + f2t * (t_yz_xyz[j] + t_xyz_yz[j]) + f2z * (s_xxyz_xyz[j] - f2a * s_yz_xyz[j]);

                s_xxyz_xzz[j] = fr * s_xyz_xzz[j] + f2t * (s_yz_xzz[j] + s_xyz_zz[j]);

                t_xxyz_xzz[j] = fr * t_xyz_xzz[j] + f2t * (t_yz_xzz[j] + t_xyz_zz[j]) + f2z * (s_xxyz_xzz[j] - f2a * s_yz_xzz[j]);

                s_xxyz_yyy[j] = fr * s_xyz_yyy[j] + f2t * s_yz_yyy[j];

                t_xxyz_yyy[j] = fr * t_xyz_yyy[j] + f2t * t_yz_yyy[j] + f2z * (s_xxyz_yyy[j] - f2a * s_yz_yyy[j]);

                s_xxyz_yyz[j] = fr * s_xyz_yyz[j] + f2t * s_yz_yyz[j];

                t_xxyz_yyz[j] = fr * t_xyz_yyz[j] + f2t * t_yz_yyz[j] + f2z * (s_xxyz_yyz[j] - f2a * s_yz_yyz[j]);

                s_xxyz_yzz[j] = fr * s_xyz_yzz[j] + f2t * s_yz_yzz[j];

                t_xxyz_yzz[j] = fr * t_xyz_yzz[j] + f2t * t_yz_yzz[j] + f2z * (s_xxyz_yzz[j] - f2a * s_yz_yzz[j]);

                s_xxyz_zzz[j] = fr * s_xyz_zzz[j] + f2t * s_yz_zzz[j];

                t_xxyz_zzz[j] = fr * t_xyz_zzz[j] + f2t * t_yz_zzz[j] + f2z * (s_xxyz_zzz[j] - f2a * s_yz_zzz[j]);

                s_xxzz_xxx[j] = fr * s_xzz_xxx[j] + f2t * (s_zz_xxx[j] + 3.0 * s_xzz_xx[j]);

                t_xxzz_xxx[j] = fr * t_xzz_xxx[j] + f2t * (t_zz_xxx[j] + 3.0 * t_xzz_xx[j]) + f2z * (s_xxzz_xxx[j] - f2a * s_zz_xxx[j]);

                s_xxzz_xxy[j] = fr * s_xzz_xxy[j] + f2t * (s_zz_xxy[j] + 2.0 * s_xzz_xy[j]);

                t_xxzz_xxy[j] = fr * t_xzz_xxy[j] + f2t * (t_zz_xxy[j] + 2.0 * t_xzz_xy[j]) + f2z * (s_xxzz_xxy[j] - f2a * s_zz_xxy[j]);

                s_xxzz_xxz[j] = fr * s_xzz_xxz[j] + f2t * (s_zz_xxz[j] + 2.0 * s_xzz_xz[j]);

                t_xxzz_xxz[j] = fr * t_xzz_xxz[j] + f2t * (t_zz_xxz[j] + 2.0 * t_xzz_xz[j]) + f2z * (s_xxzz_xxz[j] - f2a * s_zz_xxz[j]);

                s_xxzz_xyy[j] = fr * s_xzz_xyy[j] + f2t * (s_zz_xyy[j] + s_xzz_yy[j]);

                t_xxzz_xyy[j] = fr * t_xzz_xyy[j] + f2t * (t_zz_xyy[j] + t_xzz_yy[j]) + f2z * (s_xxzz_xyy[j] - f2a * s_zz_xyy[j]);

                s_xxzz_xyz[j] = fr * s_xzz_xyz[j] + f2t * (s_zz_xyz[j] + s_xzz_yz[j]);

                t_xxzz_xyz[j] = fr * t_xzz_xyz[j] + f2t * (t_zz_xyz[j] + t_xzz_yz[j]) + f2z * (s_xxzz_xyz[j] - f2a * s_zz_xyz[j]);

                s_xxzz_xzz[j] = fr * s_xzz_xzz[j] + f2t * (s_zz_xzz[j] + s_xzz_zz[j]);

                t_xxzz_xzz[j] = fr * t_xzz_xzz[j] + f2t * (t_zz_xzz[j] + t_xzz_zz[j]) + f2z * (s_xxzz_xzz[j] - f2a * s_zz_xzz[j]);

                s_xxzz_yyy[j] = fr * s_xzz_yyy[j] + f2t * s_zz_yyy[j];

                t_xxzz_yyy[j] = fr * t_xzz_yyy[j] + f2t * t_zz_yyy[j] + f2z * (s_xxzz_yyy[j] - f2a * s_zz_yyy[j]);

                s_xxzz_yyz[j] = fr * s_xzz_yyz[j] + f2t * s_zz_yyz[j];

                t_xxzz_yyz[j] = fr * t_xzz_yyz[j] + f2t * t_zz_yyz[j] + f2z * (s_xxzz_yyz[j] - f2a * s_zz_yyz[j]);

                s_xxzz_yzz[j] = fr * s_xzz_yzz[j] + f2t * s_zz_yzz[j];

                t_xxzz_yzz[j] = fr * t_xzz_yzz[j] + f2t * t_zz_yzz[j] + f2z * (s_xxzz_yzz[j] - f2a * s_zz_yzz[j]);

                s_xxzz_zzz[j] = fr * s_xzz_zzz[j] + f2t * s_zz_zzz[j];

                t_xxzz_zzz[j] = fr * t_xzz_zzz[j] + f2t * t_zz_zzz[j] + f2z * (s_xxzz_zzz[j] - f2a * s_zz_zzz[j]);

                s_xyyy_xxx[j] = fr * s_yyy_xxx[j] + f2t * 3.0 * s_yyy_xx[j];

                t_xyyy_xxx[j] = fr * t_yyy_xxx[j] + f2t * 3.0 * t_yyy_xx[j] + f2z * s_xyyy_xxx[j];

                s_xyyy_xxy[j] = fr * s_yyy_xxy[j] + f2t * 2.0 * s_yyy_xy[j];

                t_xyyy_xxy[j] = fr * t_yyy_xxy[j] + f2t * 2.0 * t_yyy_xy[j] + f2z * s_xyyy_xxy[j];

                s_xyyy_xxz[j] = fr * s_yyy_xxz[j] + f2t * 2.0 * s_yyy_xz[j];

                t_xyyy_xxz[j] = fr * t_yyy_xxz[j] + f2t * 2.0 * t_yyy_xz[j] + f2z * s_xyyy_xxz[j];

                s_xyyy_xyy[j] = fr * s_yyy_xyy[j] + f2t * s_yyy_yy[j];

                t_xyyy_xyy[j] = fr * t_yyy_xyy[j] + f2t * t_yyy_yy[j] + f2z * s_xyyy_xyy[j];

                s_xyyy_xyz[j] = fr * s_yyy_xyz[j] + f2t * s_yyy_yz[j];

                t_xyyy_xyz[j] = fr * t_yyy_xyz[j] + f2t * t_yyy_yz[j] + f2z * s_xyyy_xyz[j];

                s_xyyy_xzz[j] = fr * s_yyy_xzz[j] + f2t * s_yyy_zz[j];

                t_xyyy_xzz[j] = fr * t_yyy_xzz[j] + f2t * t_yyy_zz[j] + f2z * s_xyyy_xzz[j];

                s_xyyy_yyy[j] = fr * s_yyy_yyy[j];

                t_xyyy_yyy[j] = fr * t_yyy_yyy[j] + f2z * s_xyyy_yyy[j];

                s_xyyy_yyz[j] = fr * s_yyy_yyz[j];

                t_xyyy_yyz[j] = fr * t_yyy_yyz[j] + f2z * s_xyyy_yyz[j];

                s_xyyy_yzz[j] = fr * s_yyy_yzz[j];

                t_xyyy_yzz[j] = fr * t_yyy_yzz[j] + f2z * s_xyyy_yzz[j];

                s_xyyy_zzz[j] = fr * s_yyy_zzz[j];

                t_xyyy_zzz[j] = fr * t_yyy_zzz[j] + f2z * s_xyyy_zzz[j];

                s_xyyz_xxx[j] = fr * s_yyz_xxx[j] + f2t * 3.0 * s_yyz_xx[j];

                t_xyyz_xxx[j] = fr * t_yyz_xxx[j] + f2t * 3.0 * t_yyz_xx[j] + f2z * s_xyyz_xxx[j];

                s_xyyz_xxy[j] = fr * s_yyz_xxy[j] + f2t * 2.0 * s_yyz_xy[j];

                t_xyyz_xxy[j] = fr * t_yyz_xxy[j] + f2t * 2.0 * t_yyz_xy[j] + f2z * s_xyyz_xxy[j];

                s_xyyz_xxz[j] = fr * s_yyz_xxz[j] + f2t * 2.0 * s_yyz_xz[j];

                t_xyyz_xxz[j] = fr * t_yyz_xxz[j] + f2t * 2.0 * t_yyz_xz[j] + f2z * s_xyyz_xxz[j];

                s_xyyz_xyy[j] = fr * s_yyz_xyy[j] + f2t * s_yyz_yy[j];

                t_xyyz_xyy[j] = fr * t_yyz_xyy[j] + f2t * t_yyz_yy[j] + f2z * s_xyyz_xyy[j];

                s_xyyz_xyz[j] = fr * s_yyz_xyz[j] + f2t * s_yyz_yz[j];

                t_xyyz_xyz[j] = fr * t_yyz_xyz[j] + f2t * t_yyz_yz[j] + f2z * s_xyyz_xyz[j];

                s_xyyz_xzz[j] = fr * s_yyz_xzz[j] + f2t * s_yyz_zz[j];

                t_xyyz_xzz[j] = fr * t_yyz_xzz[j] + f2t * t_yyz_zz[j] + f2z * s_xyyz_xzz[j];

                s_xyyz_yyy[j] = fr * s_yyz_yyy[j];

                t_xyyz_yyy[j] = fr * t_yyz_yyy[j] + f2z * s_xyyz_yyy[j];

                s_xyyz_yyz[j] = fr * s_yyz_yyz[j];

                t_xyyz_yyz[j] = fr * t_yyz_yyz[j] + f2z * s_xyyz_yyz[j];

                s_xyyz_yzz[j] = fr * s_yyz_yzz[j];

                t_xyyz_yzz[j] = fr * t_yyz_yzz[j] + f2z * s_xyyz_yzz[j];

                s_xyyz_zzz[j] = fr * s_yyz_zzz[j];

                t_xyyz_zzz[j] = fr * t_yyz_zzz[j] + f2z * s_xyyz_zzz[j];

                s_xyzz_xxx[j] = fr * s_yzz_xxx[j] + f2t * 3.0 * s_yzz_xx[j];

                t_xyzz_xxx[j] = fr * t_yzz_xxx[j] + f2t * 3.0 * t_yzz_xx[j] + f2z * s_xyzz_xxx[j];

                s_xyzz_xxy[j] = fr * s_yzz_xxy[j] + f2t * 2.0 * s_yzz_xy[j];

                t_xyzz_xxy[j] = fr * t_yzz_xxy[j] + f2t * 2.0 * t_yzz_xy[j] + f2z * s_xyzz_xxy[j];

                s_xyzz_xxz[j] = fr * s_yzz_xxz[j] + f2t * 2.0 * s_yzz_xz[j];

                t_xyzz_xxz[j] = fr * t_yzz_xxz[j] + f2t * 2.0 * t_yzz_xz[j] + f2z * s_xyzz_xxz[j];

                s_xyzz_xyy[j] = fr * s_yzz_xyy[j] + f2t * s_yzz_yy[j];

                t_xyzz_xyy[j] = fr * t_yzz_xyy[j] + f2t * t_yzz_yy[j] + f2z * s_xyzz_xyy[j];

                s_xyzz_xyz[j] = fr * s_yzz_xyz[j] + f2t * s_yzz_yz[j];

                t_xyzz_xyz[j] = fr * t_yzz_xyz[j] + f2t * t_yzz_yz[j] + f2z * s_xyzz_xyz[j];

                s_xyzz_xzz[j] = fr * s_yzz_xzz[j] + f2t * s_yzz_zz[j];

                t_xyzz_xzz[j] = fr * t_yzz_xzz[j] + f2t * t_yzz_zz[j] + f2z * s_xyzz_xzz[j];

                s_xyzz_yyy[j] = fr * s_yzz_yyy[j];

                t_xyzz_yyy[j] = fr * t_yzz_yyy[j] + f2z * s_xyzz_yyy[j];

                s_xyzz_yyz[j] = fr * s_yzz_yyz[j];

                t_xyzz_yyz[j] = fr * t_yzz_yyz[j] + f2z * s_xyzz_yyz[j];

                s_xyzz_yzz[j] = fr * s_yzz_yzz[j];

                t_xyzz_yzz[j] = fr * t_yzz_yzz[j] + f2z * s_xyzz_yzz[j];

                s_xyzz_zzz[j] = fr * s_yzz_zzz[j];

                t_xyzz_zzz[j] = fr * t_yzz_zzz[j] + f2z * s_xyzz_zzz[j];

                s_xzzz_xxx[j] = fr * s_zzz_xxx[j] + f2t * 3.0 * s_zzz_xx[j];

                t_xzzz_xxx[j] = fr * t_zzz_xxx[j] + f2t * 3.0 * t_zzz_xx[j] + f2z * s_xzzz_xxx[j];

                s_xzzz_xxy[j] = fr * s_zzz_xxy[j] + f2t * 2.0 * s_zzz_xy[j];

                t_xzzz_xxy[j] = fr * t_zzz_xxy[j] + f2t * 2.0 * t_zzz_xy[j] + f2z * s_xzzz_xxy[j];

                s_xzzz_xxz[j] = fr * s_zzz_xxz[j] + f2t * 2.0 * s_zzz_xz[j];

                t_xzzz_xxz[j] = fr * t_zzz_xxz[j] + f2t * 2.0 * t_zzz_xz[j] + f2z * s_xzzz_xxz[j];

                s_xzzz_xyy[j] = fr * s_zzz_xyy[j] + f2t * s_zzz_yy[j];

                t_xzzz_xyy[j] = fr * t_zzz_xyy[j] + f2t * t_zzz_yy[j] + f2z * s_xzzz_xyy[j];

                s_xzzz_xyz[j] = fr * s_zzz_xyz[j] + f2t * s_zzz_yz[j];

                t_xzzz_xyz[j] = fr * t_zzz_xyz[j] + f2t * t_zzz_yz[j] + f2z * s_xzzz_xyz[j];

                s_xzzz_xzz[j] = fr * s_zzz_xzz[j] + f2t * s_zzz_zz[j];

                t_xzzz_xzz[j] = fr * t_zzz_xzz[j] + f2t * t_zzz_zz[j] + f2z * s_xzzz_xzz[j];

                s_xzzz_yyy[j] = fr * s_zzz_yyy[j];

                t_xzzz_yyy[j] = fr * t_zzz_yyy[j] + f2z * s_xzzz_yyy[j];

                s_xzzz_yyz[j] = fr * s_zzz_yyz[j];

                t_xzzz_yyz[j] = fr * t_zzz_yyz[j] + f2z * s_xzzz_yyz[j];

                s_xzzz_yzz[j] = fr * s_zzz_yzz[j];

                t_xzzz_yzz[j] = fr * t_zzz_yzz[j] + f2z * s_xzzz_yzz[j];

                s_xzzz_zzz[j] = fr * s_zzz_zzz[j];

                t_xzzz_zzz[j] = fr * t_zzz_zzz[j] + f2z * s_xzzz_zzz[j];

                // leading y component

                fr = pay[j];

                s_yyyy_xxx[j] = fr * s_yyy_xxx[j] + f2t * 3.0 * s_yy_xxx[j];

                t_yyyy_xxx[j] = fr * t_yyy_xxx[j] + f2t * 3.0 * t_yy_xxx[j] + f2z * (s_yyyy_xxx[j] - f2a * 3.0 * s_yy_xxx[j]);

                s_yyyy_xxy[j] = fr * s_yyy_xxy[j] + f2t * (3.0 * s_yy_xxy[j] + s_yyy_xx[j]);

                t_yyyy_xxy[j] = fr * t_yyy_xxy[j] + f2t * (3.0 * t_yy_xxy[j] + t_yyy_xx[j]) + f2z * (s_yyyy_xxy[j] - f2a * 3.0 * s_yy_xxy[j]);

                s_yyyy_xxz[j] = fr * s_yyy_xxz[j] + f2t * 3.0 * s_yy_xxz[j];

                t_yyyy_xxz[j] = fr * t_yyy_xxz[j] + f2t * 3.0 * t_yy_xxz[j] + f2z * (s_yyyy_xxz[j] - f2a * 3.0 * s_yy_xxz[j]);

                s_yyyy_xyy[j] = fr * s_yyy_xyy[j] + f2t * (3.0 * s_yy_xyy[j] + 2.0 * s_yyy_xy[j]);

                t_yyyy_xyy[j] = fr * t_yyy_xyy[j] + f2t * (3.0 * t_yy_xyy[j] + 2.0 * t_yyy_xy[j]) + f2z * (s_yyyy_xyy[j] - f2a * 3.0 * s_yy_xyy[j]);

                s_yyyy_xyz[j] = fr * s_yyy_xyz[j] + f2t * (3.0 * s_yy_xyz[j] + s_yyy_xz[j]);

                t_yyyy_xyz[j] = fr * t_yyy_xyz[j] + f2t * (3.0 * t_yy_xyz[j] + t_yyy_xz[j]) + f2z * (s_yyyy_xyz[j] - f2a * 3.0 * s_yy_xyz[j]);

                s_yyyy_xzz[j] = fr * s_yyy_xzz[j] + f2t * 3.0 * s_yy_xzz[j];

                t_yyyy_xzz[j] = fr * t_yyy_xzz[j] + f2t * 3.0 * t_yy_xzz[j] + f2z * (s_yyyy_xzz[j] - f2a * 3.0 * s_yy_xzz[j]);

                s_yyyy_yyy[j] = fr * s_yyy_yyy[j] + f2t * (3.0 * s_yy_yyy[j] + 3.0 * s_yyy_yy[j]);

                t_yyyy_yyy[j] = fr * t_yyy_yyy[j] + f2t * (3.0 * t_yy_yyy[j] + 3.0 * t_yyy_yy[j]) + f2z * (s_yyyy_yyy[j] - f2a * 3.0 * s_yy_yyy[j]);

                s_yyyy_yyz[j] = fr * s_yyy_yyz[j] + f2t * (3.0 * s_yy_yyz[j] + 2.0 * s_yyy_yz[j]);

                t_yyyy_yyz[j] = fr * t_yyy_yyz[j] + f2t * (3.0 * t_yy_yyz[j] + 2.0 * t_yyy_yz[j]) + f2z * (s_yyyy_yyz[j] - f2a * 3.0 * s_yy_yyz[j]);

                s_yyyy_yzz[j] = fr * s_yyy_yzz[j] + f2t * (3.0 * s_yy_yzz[j] + s_yyy_zz[j]);

                t_yyyy_yzz[j] = fr * t_yyy_yzz[j] + f2t * (3.0 * t_yy_yzz[j] + t_yyy_zz[j]) + f2z * (s_yyyy_yzz[j] - f2a * 3.0 * s_yy_yzz[j]);

                s_yyyy_zzz[j] = fr * s_yyy_zzz[j] + f2t * 3.0 * s_yy_zzz[j];

                t_yyyy_zzz[j] = fr * t_yyy_zzz[j] + f2t * 3.0 * t_yy_zzz[j] + f2z * (s_yyyy_zzz[j] - f2a * 3.0 * s_yy_zzz[j]);

                s_yyyz_xxx[j] = fr * s_yyz_xxx[j] + f2t * 2.0 * s_yz_xxx[j];

                t_yyyz_xxx[j] = fr * t_yyz_xxx[j] + f2t * 2.0 * t_yz_xxx[j] + f2z * (s_yyyz_xxx[j] - f2a * 2.0 * s_yz_xxx[j]);

                s_yyyz_xxy[j] = fr * s_yyz_xxy[j] + f2t * (2.0 * s_yz_xxy[j] + s_yyz_xx[j]);

                t_yyyz_xxy[j] = fr * t_yyz_xxy[j] + f2t * (2.0 * t_yz_xxy[j] + t_yyz_xx[j]) + f2z * (s_yyyz_xxy[j] - f2a * 2.0 * s_yz_xxy[j]);

                s_yyyz_xxz[j] = fr * s_yyz_xxz[j] + f2t * 2.0 * s_yz_xxz[j];

                t_yyyz_xxz[j] = fr * t_yyz_xxz[j] + f2t * 2.0 * t_yz_xxz[j] + f2z * (s_yyyz_xxz[j] - f2a * 2.0 * s_yz_xxz[j]);

                s_yyyz_xyy[j] = fr * s_yyz_xyy[j] + f2t * (2.0 * s_yz_xyy[j] + 2.0 * s_yyz_xy[j]);

                t_yyyz_xyy[j] = fr * t_yyz_xyy[j] + f2t * (2.0 * t_yz_xyy[j] + 2.0 * t_yyz_xy[j]) + f2z * (s_yyyz_xyy[j] - f2a * 2.0 * s_yz_xyy[j]);

                s_yyyz_xyz[j] = fr * s_yyz_xyz[j] + f2t * (2.0 * s_yz_xyz[j] + s_yyz_xz[j]);

                t_yyyz_xyz[j] = fr * t_yyz_xyz[j] + f2t * (2.0 * t_yz_xyz[j] + t_yyz_xz[j]) + f2z * (s_yyyz_xyz[j] - f2a * 2.0 * s_yz_xyz[j]);

                s_yyyz_xzz[j] = fr * s_yyz_xzz[j] + f2t * 2.0 * s_yz_xzz[j];

                t_yyyz_xzz[j] = fr * t_yyz_xzz[j] + f2t * 2.0 * t_yz_xzz[j] + f2z * (s_yyyz_xzz[j] - f2a * 2.0 * s_yz_xzz[j]);

                s_yyyz_yyy[j] = fr * s_yyz_yyy[j] + f2t * (2.0 * s_yz_yyy[j] + 3.0 * s_yyz_yy[j]);

                t_yyyz_yyy[j] = fr * t_yyz_yyy[j] + f2t * (2.0 * t_yz_yyy[j] + 3.0 * t_yyz_yy[j]) + f2z * (s_yyyz_yyy[j] - f2a * 2.0 * s_yz_yyy[j]);

                s_yyyz_yyz[j] = fr * s_yyz_yyz[j] + f2t * (2.0 * s_yz_yyz[j] + 2.0 * s_yyz_yz[j]);

                t_yyyz_yyz[j] = fr * t_yyz_yyz[j] + f2t * (2.0 * t_yz_yyz[j] + 2.0 * t_yyz_yz[j]) + f2z * (s_yyyz_yyz[j] - f2a * 2.0 * s_yz_yyz[j]);

                s_yyyz_yzz[j] = fr * s_yyz_yzz[j] + f2t * (2.0 * s_yz_yzz[j] + s_yyz_zz[j]);

                t_yyyz_yzz[j] = fr * t_yyz_yzz[j] + f2t * (2.0 * t_yz_yzz[j] + t_yyz_zz[j]) + f2z * (s_yyyz_yzz[j] - f2a * 2.0 * s_yz_yzz[j]);

                s_yyyz_zzz[j] = fr * s_yyz_zzz[j] + f2t * 2.0 * s_yz_zzz[j];

                t_yyyz_zzz[j] = fr * t_yyz_zzz[j] + f2t * 2.0 * t_yz_zzz[j] + f2z * (s_yyyz_zzz[j] - f2a * 2.0 * s_yz_zzz[j]);

                s_yyzz_xxx[j] = fr * s_yzz_xxx[j] + f2t * s_zz_xxx[j];

                t_yyzz_xxx[j] = fr * t_yzz_xxx[j] + f2t * t_zz_xxx[j] + f2z * (s_yyzz_xxx[j] - f2a * s_zz_xxx[j]);

                s_yyzz_xxy[j] = fr * s_yzz_xxy[j] + f2t * (s_zz_xxy[j] + s_yzz_xx[j]);

                t_yyzz_xxy[j] = fr * t_yzz_xxy[j] + f2t * (t_zz_xxy[j] + t_yzz_xx[j]) + f2z * (s_yyzz_xxy[j] - f2a * s_zz_xxy[j]);

                s_yyzz_xxz[j] = fr * s_yzz_xxz[j] + f2t * s_zz_xxz[j];

                t_yyzz_xxz[j] = fr * t_yzz_xxz[j] + f2t * t_zz_xxz[j] + f2z * (s_yyzz_xxz[j] - f2a * s_zz_xxz[j]);

                s_yyzz_xyy[j] = fr * s_yzz_xyy[j] + f2t * (s_zz_xyy[j] + 2.0 * s_yzz_xy[j]);

                t_yyzz_xyy[j] = fr * t_yzz_xyy[j] + f2t * (t_zz_xyy[j] + 2.0 * t_yzz_xy[j]) + f2z * (s_yyzz_xyy[j] - f2a * s_zz_xyy[j]);

                s_yyzz_xyz[j] = fr * s_yzz_xyz[j] + f2t * (s_zz_xyz[j] + s_yzz_xz[j]);

                t_yyzz_xyz[j] = fr * t_yzz_xyz[j] + f2t * (t_zz_xyz[j] + t_yzz_xz[j]) + f2z * (s_yyzz_xyz[j] - f2a * s_zz_xyz[j]);

                s_yyzz_xzz[j] = fr * s_yzz_xzz[j] + f2t * s_zz_xzz[j];

                t_yyzz_xzz[j] = fr * t_yzz_xzz[j] + f2t * t_zz_xzz[j] + f2z * (s_yyzz_xzz[j] - f2a * s_zz_xzz[j]);

                s_yyzz_yyy[j] = fr * s_yzz_yyy[j] + f2t * (s_zz_yyy[j] + 3.0 * s_yzz_yy[j]);

                t_yyzz_yyy[j] = fr * t_yzz_yyy[j] + f2t * (t_zz_yyy[j] + 3.0 * t_yzz_yy[j]) + f2z * (s_yyzz_yyy[j] - f2a * s_zz_yyy[j]);

                s_yyzz_yyz[j] = fr * s_yzz_yyz[j] + f2t * (s_zz_yyz[j] + 2.0 * s_yzz_yz[j]);

                t_yyzz_yyz[j] = fr * t_yzz_yyz[j] + f2t * (t_zz_yyz[j] + 2.0 * t_yzz_yz[j]) + f2z * (s_yyzz_yyz[j] - f2a * s_zz_yyz[j]);

                s_yyzz_yzz[j] = fr * s_yzz_yzz[j] + f2t * (s_zz_yzz[j] + s_yzz_zz[j]);

                t_yyzz_yzz[j] = fr * t_yzz_yzz[j] + f2t * (t_zz_yzz[j] + t_yzz_zz[j]) + f2z * (s_yyzz_yzz[j] - f2a * s_zz_yzz[j]);

                s_yyzz_zzz[j] = fr * s_yzz_zzz[j] + f2t * s_zz_zzz[j];

                t_yyzz_zzz[j] = fr * t_yzz_zzz[j] + f2t * t_zz_zzz[j] + f2z * (s_yyzz_zzz[j] - f2a * s_zz_zzz[j]);

                s_yzzz_xxx[j] = fr * s_zzz_xxx[j];

                t_yzzz_xxx[j] = fr * t_zzz_xxx[j] + f2z * s_yzzz_xxx[j];

                s_yzzz_xxy[j] = fr * s_zzz_xxy[j] + f2t * s_zzz_xx[j];

                t_yzzz_xxy[j] = fr * t_zzz_xxy[j] + f2t * t_zzz_xx[j] + f2z * s_yzzz_xxy[j];

                s_yzzz_xxz[j] = fr * s_zzz_xxz[j];

                t_yzzz_xxz[j] = fr * t_zzz_xxz[j] + f2z * s_yzzz_xxz[j];

                s_yzzz_xyy[j] = fr * s_zzz_xyy[j] + f2t * 2.0 * s_zzz_xy[j];

                t_yzzz_xyy[j] = fr * t_zzz_xyy[j] + f2t * 2.0 * t_zzz_xy[j] + f2z * s_yzzz_xyy[j];

                s_yzzz_xyz[j] = fr * s_zzz_xyz[j] + f2t * s_zzz_xz[j];

                t_yzzz_xyz[j] = fr * t_zzz_xyz[j] + f2t * t_zzz_xz[j] + f2z * s_yzzz_xyz[j];

                s_yzzz_xzz[j] = fr * s_zzz_xzz[j];

                t_yzzz_xzz[j] = fr * t_zzz_xzz[j] + f2z * s_yzzz_xzz[j];

                s_yzzz_yyy[j] = fr * s_zzz_yyy[j] + f2t * 3.0 * s_zzz_yy[j];

                t_yzzz_yyy[j] = fr * t_zzz_yyy[j] + f2t * 3.0 * t_zzz_yy[j] + f2z * s_yzzz_yyy[j];

                s_yzzz_yyz[j] = fr * s_zzz_yyz[j] + f2t * 2.0 * s_zzz_yz[j];

                t_yzzz_yyz[j] = fr * t_zzz_yyz[j] + f2t * 2.0 * t_zzz_yz[j] + f2z * s_yzzz_yyz[j];

                s_yzzz_yzz[j] = fr * s_zzz_yzz[j] + f2t * s_zzz_zz[j];

                t_yzzz_yzz[j] = fr * t_zzz_yzz[j] + f2t * t_zzz_zz[j] + f2z * s_yzzz_yzz[j];

                s_yzzz_zzz[j] = fr * s_zzz_zzz[j];

                t_yzzz_zzz[j] = fr * t_zzz_zzz[j] + f2z * s_yzzz_zzz[j];

                // leading z component

                fr = paz[j];

                s_zzzz_xxx[j] = fr * s_zzz_xxx[j] + f2t * 3.0 * s_zz_xxx[j];

                t_zzzz_xxx[j] = fr * t_zzz_xxx[j] + f2t * 3.0 * t_zz_xxx[j] + f2z * (s_zzzz_xxx[j] - f2a * 3.0 * s_zz_xxx[j]);

                s_zzzz_xxy[j] = fr * s_zzz_xxy[j] + f2t * 3.0 * s_zz_xxy[j];

                t_zzzz_xxy[j] = fr * t_zzz_xxy[j] + f2t * 3.0 * t_zz_xxy[j] + f2z * (s_zzzz_xxy[j] - f2a * 3.0 * s_zz_xxy[j]);

                s_zzzz_xxz[j] = fr * s_zzz_xxz[j] + f2t * (3.0 * s_zz_xxz[j] + s_zzz_xx[j]);

                t_zzzz_xxz[j] = fr * t_zzz_xxz[j] + f2t * (3.0 * t_zz_xxz[j] + t_zzz_xx[j]) + f2z * (s_zzzz_xxz[j] - f2a * 3.0 * s_zz_xxz[j]);

                s_zzzz_xyy[j] = fr * s_zzz_xyy[j] + f2t * 3.0 * s_zz_xyy[j];

                t_zzzz_xyy[j] = fr * t_zzz_xyy[j] + f2t * 3.0 * t_zz_xyy[j] + f2z * (s_zzzz_xyy[j] - f2a * 3.0 * s_zz_xyy[j]);

                s_zzzz_xyz[j] = fr * s_zzz_xyz[j] + f2t * (3.0 * s_zz_xyz[j] + s_zzz_xy[j]);

                t_zzzz_xyz[j] = fr * t_zzz_xyz[j] + f2t * (3.0 * t_zz_xyz[j] + t_zzz_xy[j]) + f2z * (s_zzzz_xyz[j] - f2a * 3.0 * s_zz_xyz[j]);

                s_zzzz_xzz[j] = fr * s_zzz_xzz[j] + f2t * (3.0 * s_zz_xzz[j] + 2.0 * s_zzz_xz[j]);

                t_zzzz_xzz[j] = fr * t_zzz_xzz[j] + f2t * (3.0 * t_zz_xzz[j] + 2.0 * t_zzz_xz[j]) + f2z * (s_zzzz_xzz[j] - f2a * 3.0 * s_zz_xzz[j]);

                s_zzzz_yyy[j] = fr * s_zzz_yyy[j] + f2t * 3.0 * s_zz_yyy[j];

                t_zzzz_yyy[j] = fr * t_zzz_yyy[j] + f2t * 3.0 * t_zz_yyy[j] + f2z * (s_zzzz_yyy[j] - f2a * 3.0 * s_zz_yyy[j]);

                s_zzzz_yyz[j] = fr * s_zzz_yyz[j] + f2t * (3.0 * s_zz_yyz[j] + s_zzz_yy[j]);

                t_zzzz_yyz[j] = fr * t_zzz_yyz[j] + f2t * (3.0 * t_zz_yyz[j] + t_zzz_yy[j]) + f2z * (s_zzzz_yyz[j] - f2a * 3.0 * s_zz_yyz[j]);

                s_zzzz_yzz[j] = fr * s_zzz_yzz[j] + f2t * (3.0 * s_zz_yzz[j] + 2.0 * s_zzz_yz[j]);

                t_zzzz_yzz[j] = fr * t_zzz_yzz[j] + f2t * (3.0 * t_zz_yzz[j] + 2.0 * t_zzz_yz[j]) + f2z * (s_zzzz_yzz[j] - f2a * 3.0 * s_zz_yzz[j]);

                s_zzzz_zzz[j] = fr * s_zzz_zzz[j] + f2t * (3.0 * s_zz_zzz[j] + 3.0 * s_zzz_zz[j]);

                t_zzzz_zzz[j] = fr * t_zzz_zzz[j] + f2t * (3.0 * t_zz_zzz[j] + 3.0 * t_zzz_zz[j]) + f2z * (s_zzzz_zzz[j] - f2a * 3.0 * s_zz_zzz[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGG(      CMemBlock2D<double>&  primBuffer,
                           const CVecThreeIndexes&     recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  paDistances,
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

        // get position of integrals in primitves buffer

        auto soff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 4, 1});

        auto koff  = genfunc::findTripleIndex(recIndexes, recPattern, {4, 4, 0});

        auto s1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, 1});

        auto k1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, 0});

        auto s2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 1});

        auto k2off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, 0});

        auto skoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 1});

        auto kkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, 0});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|F) integrals

            auto s_xxx_xxx = primBuffer.data(skoff + 100 * idx);

            auto s_xxx_xxy = primBuffer.data(skoff + 100 * idx + 1);

            auto s_xxx_xxz = primBuffer.data(skoff + 100 * idx + 2);

            auto s_xxx_xyy = primBuffer.data(skoff + 100 * idx + 3);

            auto s_xxx_xyz = primBuffer.data(skoff + 100 * idx + 4);

            auto s_xxx_xzz = primBuffer.data(skoff + 100 * idx + 5);

            auto s_xxx_yyy = primBuffer.data(skoff + 100 * idx + 6);

            auto s_xxx_yyz = primBuffer.data(skoff + 100 * idx + 7);

            auto s_xxx_yzz = primBuffer.data(skoff + 100 * idx + 8);

            auto s_xxx_zzz = primBuffer.data(skoff + 100 * idx + 9);

            auto s_xxy_xxx = primBuffer.data(skoff + 100 * idx + 10);

            auto s_xxy_xxy = primBuffer.data(skoff + 100 * idx + 11);

            auto s_xxy_xxz = primBuffer.data(skoff + 100 * idx + 12);

            auto s_xxy_xyy = primBuffer.data(skoff + 100 * idx + 13);

            auto s_xxy_xyz = primBuffer.data(skoff + 100 * idx + 14);

            auto s_xxy_xzz = primBuffer.data(skoff + 100 * idx + 15);

            auto s_xxy_yyy = primBuffer.data(skoff + 100 * idx + 16);

            auto s_xxy_yyz = primBuffer.data(skoff + 100 * idx + 17);

            auto s_xxy_yzz = primBuffer.data(skoff + 100 * idx + 18);

            auto s_xxy_zzz = primBuffer.data(skoff + 100 * idx + 19);

            auto s_xxz_xxx = primBuffer.data(skoff + 100 * idx + 20);

            auto s_xxz_xxy = primBuffer.data(skoff + 100 * idx + 21);

            auto s_xxz_xxz = primBuffer.data(skoff + 100 * idx + 22);

            auto s_xxz_xyy = primBuffer.data(skoff + 100 * idx + 23);

            auto s_xxz_xyz = primBuffer.data(skoff + 100 * idx + 24);

            auto s_xxz_xzz = primBuffer.data(skoff + 100 * idx + 25);

            auto s_xxz_yyy = primBuffer.data(skoff + 100 * idx + 26);

            auto s_xxz_yyz = primBuffer.data(skoff + 100 * idx + 27);

            auto s_xxz_yzz = primBuffer.data(skoff + 100 * idx + 28);

            auto s_xxz_zzz = primBuffer.data(skoff + 100 * idx + 29);

            auto s_xyy_xxx = primBuffer.data(skoff + 100 * idx + 30);

            auto s_xyy_xxy = primBuffer.data(skoff + 100 * idx + 31);

            auto s_xyy_xxz = primBuffer.data(skoff + 100 * idx + 32);

            auto s_xyy_xyy = primBuffer.data(skoff + 100 * idx + 33);

            auto s_xyy_xyz = primBuffer.data(skoff + 100 * idx + 34);

            auto s_xyy_xzz = primBuffer.data(skoff + 100 * idx + 35);

            auto s_xyy_yyy = primBuffer.data(skoff + 100 * idx + 36);

            auto s_xyy_yyz = primBuffer.data(skoff + 100 * idx + 37);

            auto s_xyy_yzz = primBuffer.data(skoff + 100 * idx + 38);

            auto s_xyy_zzz = primBuffer.data(skoff + 100 * idx + 39);

            auto s_xyz_xxx = primBuffer.data(skoff + 100 * idx + 40);

            auto s_xyz_xxy = primBuffer.data(skoff + 100 * idx + 41);

            auto s_xyz_xxz = primBuffer.data(skoff + 100 * idx + 42);

            auto s_xyz_xyy = primBuffer.data(skoff + 100 * idx + 43);

            auto s_xyz_xyz = primBuffer.data(skoff + 100 * idx + 44);

            auto s_xyz_xzz = primBuffer.data(skoff + 100 * idx + 45);

            auto s_xyz_yyy = primBuffer.data(skoff + 100 * idx + 46);

            auto s_xyz_yyz = primBuffer.data(skoff + 100 * idx + 47);

            auto s_xyz_yzz = primBuffer.data(skoff + 100 * idx + 48);

            auto s_xyz_zzz = primBuffer.data(skoff + 100 * idx + 49);

            auto s_xzz_xxx = primBuffer.data(skoff + 100 * idx + 50);

            auto s_xzz_xxy = primBuffer.data(skoff + 100 * idx + 51);

            auto s_xzz_xxz = primBuffer.data(skoff + 100 * idx + 52);

            auto s_xzz_xyy = primBuffer.data(skoff + 100 * idx + 53);

            auto s_xzz_xyz = primBuffer.data(skoff + 100 * idx + 54);

            auto s_xzz_xzz = primBuffer.data(skoff + 100 * idx + 55);

            auto s_xzz_yyy = primBuffer.data(skoff + 100 * idx + 56);

            auto s_xzz_yyz = primBuffer.data(skoff + 100 * idx + 57);

            auto s_xzz_yzz = primBuffer.data(skoff + 100 * idx + 58);

            auto s_xzz_zzz = primBuffer.data(skoff + 100 * idx + 59);

            auto s_yyy_xxx = primBuffer.data(skoff + 100 * idx + 60);

            auto s_yyy_xxy = primBuffer.data(skoff + 100 * idx + 61);

            auto s_yyy_xxz = primBuffer.data(skoff + 100 * idx + 62);

            auto s_yyy_xyy = primBuffer.data(skoff + 100 * idx + 63);

            auto s_yyy_xyz = primBuffer.data(skoff + 100 * idx + 64);

            auto s_yyy_xzz = primBuffer.data(skoff + 100 * idx + 65);

            auto s_yyy_yyy = primBuffer.data(skoff + 100 * idx + 66);

            auto s_yyy_yyz = primBuffer.data(skoff + 100 * idx + 67);

            auto s_yyy_yzz = primBuffer.data(skoff + 100 * idx + 68);

            auto s_yyy_zzz = primBuffer.data(skoff + 100 * idx + 69);

            auto s_yyz_xxx = primBuffer.data(skoff + 100 * idx + 70);

            auto s_yyz_xxy = primBuffer.data(skoff + 100 * idx + 71);

            auto s_yyz_xxz = primBuffer.data(skoff + 100 * idx + 72);

            auto s_yyz_xyy = primBuffer.data(skoff + 100 * idx + 73);

            auto s_yyz_xyz = primBuffer.data(skoff + 100 * idx + 74);

            auto s_yyz_xzz = primBuffer.data(skoff + 100 * idx + 75);

            auto s_yyz_yyy = primBuffer.data(skoff + 100 * idx + 76);

            auto s_yyz_yyz = primBuffer.data(skoff + 100 * idx + 77);

            auto s_yyz_yzz = primBuffer.data(skoff + 100 * idx + 78);

            auto s_yyz_zzz = primBuffer.data(skoff + 100 * idx + 79);

            auto s_yzz_xxx = primBuffer.data(skoff + 100 * idx + 80);

            auto s_yzz_xxy = primBuffer.data(skoff + 100 * idx + 81);

            auto s_yzz_xxz = primBuffer.data(skoff + 100 * idx + 82);

            auto s_yzz_xyy = primBuffer.data(skoff + 100 * idx + 83);

            auto s_yzz_xyz = primBuffer.data(skoff + 100 * idx + 84);

            auto s_yzz_xzz = primBuffer.data(skoff + 100 * idx + 85);

            auto s_yzz_yyy = primBuffer.data(skoff + 100 * idx + 86);

            auto s_yzz_yyz = primBuffer.data(skoff + 100 * idx + 87);

            auto s_yzz_yzz = primBuffer.data(skoff + 100 * idx + 88);

            auto s_yzz_zzz = primBuffer.data(skoff + 100 * idx + 89);

            auto s_zzz_xxx = primBuffer.data(skoff + 100 * idx + 90);

            auto s_zzz_xxy = primBuffer.data(skoff + 100 * idx + 91);

            auto s_zzz_xxz = primBuffer.data(skoff + 100 * idx + 92);

            auto s_zzz_xyy = primBuffer.data(skoff + 100 * idx + 93);

            auto s_zzz_xyz = primBuffer.data(skoff + 100 * idx + 94);

            auto s_zzz_xzz = primBuffer.data(skoff + 100 * idx + 95);

            auto s_zzz_yyy = primBuffer.data(skoff + 100 * idx + 96);

            auto s_zzz_yyz = primBuffer.data(skoff + 100 * idx + 97);

            auto s_zzz_yzz = primBuffer.data(skoff + 100 * idx + 98);

            auto s_zzz_zzz = primBuffer.data(skoff + 100 * idx + 99);

            // set up pointers to (F|T|F) integrals

            auto t_xxx_xxx = primBuffer.data(kkoff + 100 * idx);

            auto t_xxx_xxy = primBuffer.data(kkoff + 100 * idx + 1);

            auto t_xxx_xxz = primBuffer.data(kkoff + 100 * idx + 2);

            auto t_xxx_xyy = primBuffer.data(kkoff + 100 * idx + 3);

            auto t_xxx_xyz = primBuffer.data(kkoff + 100 * idx + 4);

            auto t_xxx_xzz = primBuffer.data(kkoff + 100 * idx + 5);

            auto t_xxx_yyy = primBuffer.data(kkoff + 100 * idx + 6);

            auto t_xxx_yyz = primBuffer.data(kkoff + 100 * idx + 7);

            auto t_xxx_yzz = primBuffer.data(kkoff + 100 * idx + 8);

            auto t_xxx_zzz = primBuffer.data(kkoff + 100 * idx + 9);

            auto t_xxy_xxx = primBuffer.data(kkoff + 100 * idx + 10);

            auto t_xxy_xxy = primBuffer.data(kkoff + 100 * idx + 11);

            auto t_xxy_xxz = primBuffer.data(kkoff + 100 * idx + 12);

            auto t_xxy_xyy = primBuffer.data(kkoff + 100 * idx + 13);

            auto t_xxy_xyz = primBuffer.data(kkoff + 100 * idx + 14);

            auto t_xxy_xzz = primBuffer.data(kkoff + 100 * idx + 15);

            auto t_xxy_yyy = primBuffer.data(kkoff + 100 * idx + 16);

            auto t_xxy_yyz = primBuffer.data(kkoff + 100 * idx + 17);

            auto t_xxy_yzz = primBuffer.data(kkoff + 100 * idx + 18);

            auto t_xxy_zzz = primBuffer.data(kkoff + 100 * idx + 19);

            auto t_xxz_xxx = primBuffer.data(kkoff + 100 * idx + 20);

            auto t_xxz_xxy = primBuffer.data(kkoff + 100 * idx + 21);

            auto t_xxz_xxz = primBuffer.data(kkoff + 100 * idx + 22);

            auto t_xxz_xyy = primBuffer.data(kkoff + 100 * idx + 23);

            auto t_xxz_xyz = primBuffer.data(kkoff + 100 * idx + 24);

            auto t_xxz_xzz = primBuffer.data(kkoff + 100 * idx + 25);

            auto t_xxz_yyy = primBuffer.data(kkoff + 100 * idx + 26);

            auto t_xxz_yyz = primBuffer.data(kkoff + 100 * idx + 27);

            auto t_xxz_yzz = primBuffer.data(kkoff + 100 * idx + 28);

            auto t_xxz_zzz = primBuffer.data(kkoff + 100 * idx + 29);

            auto t_xyy_xxx = primBuffer.data(kkoff + 100 * idx + 30);

            auto t_xyy_xxy = primBuffer.data(kkoff + 100 * idx + 31);

            auto t_xyy_xxz = primBuffer.data(kkoff + 100 * idx + 32);

            auto t_xyy_xyy = primBuffer.data(kkoff + 100 * idx + 33);

            auto t_xyy_xyz = primBuffer.data(kkoff + 100 * idx + 34);

            auto t_xyy_xzz = primBuffer.data(kkoff + 100 * idx + 35);

            auto t_xyy_yyy = primBuffer.data(kkoff + 100 * idx + 36);

            auto t_xyy_yyz = primBuffer.data(kkoff + 100 * idx + 37);

            auto t_xyy_yzz = primBuffer.data(kkoff + 100 * idx + 38);

            auto t_xyy_zzz = primBuffer.data(kkoff + 100 * idx + 39);

            auto t_xyz_xxx = primBuffer.data(kkoff + 100 * idx + 40);

            auto t_xyz_xxy = primBuffer.data(kkoff + 100 * idx + 41);

            auto t_xyz_xxz = primBuffer.data(kkoff + 100 * idx + 42);

            auto t_xyz_xyy = primBuffer.data(kkoff + 100 * idx + 43);

            auto t_xyz_xyz = primBuffer.data(kkoff + 100 * idx + 44);

            auto t_xyz_xzz = primBuffer.data(kkoff + 100 * idx + 45);

            auto t_xyz_yyy = primBuffer.data(kkoff + 100 * idx + 46);

            auto t_xyz_yyz = primBuffer.data(kkoff + 100 * idx + 47);

            auto t_xyz_yzz = primBuffer.data(kkoff + 100 * idx + 48);

            auto t_xyz_zzz = primBuffer.data(kkoff + 100 * idx + 49);

            auto t_xzz_xxx = primBuffer.data(kkoff + 100 * idx + 50);

            auto t_xzz_xxy = primBuffer.data(kkoff + 100 * idx + 51);

            auto t_xzz_xxz = primBuffer.data(kkoff + 100 * idx + 52);

            auto t_xzz_xyy = primBuffer.data(kkoff + 100 * idx + 53);

            auto t_xzz_xyz = primBuffer.data(kkoff + 100 * idx + 54);

            auto t_xzz_xzz = primBuffer.data(kkoff + 100 * idx + 55);

            auto t_xzz_yyy = primBuffer.data(kkoff + 100 * idx + 56);

            auto t_xzz_yyz = primBuffer.data(kkoff + 100 * idx + 57);

            auto t_xzz_yzz = primBuffer.data(kkoff + 100 * idx + 58);

            auto t_xzz_zzz = primBuffer.data(kkoff + 100 * idx + 59);

            auto t_yyy_xxx = primBuffer.data(kkoff + 100 * idx + 60);

            auto t_yyy_xxy = primBuffer.data(kkoff + 100 * idx + 61);

            auto t_yyy_xxz = primBuffer.data(kkoff + 100 * idx + 62);

            auto t_yyy_xyy = primBuffer.data(kkoff + 100 * idx + 63);

            auto t_yyy_xyz = primBuffer.data(kkoff + 100 * idx + 64);

            auto t_yyy_xzz = primBuffer.data(kkoff + 100 * idx + 65);

            auto t_yyy_yyy = primBuffer.data(kkoff + 100 * idx + 66);

            auto t_yyy_yyz = primBuffer.data(kkoff + 100 * idx + 67);

            auto t_yyy_yzz = primBuffer.data(kkoff + 100 * idx + 68);

            auto t_yyy_zzz = primBuffer.data(kkoff + 100 * idx + 69);

            auto t_yyz_xxx = primBuffer.data(kkoff + 100 * idx + 70);

            auto t_yyz_xxy = primBuffer.data(kkoff + 100 * idx + 71);

            auto t_yyz_xxz = primBuffer.data(kkoff + 100 * idx + 72);

            auto t_yyz_xyy = primBuffer.data(kkoff + 100 * idx + 73);

            auto t_yyz_xyz = primBuffer.data(kkoff + 100 * idx + 74);

            auto t_yyz_xzz = primBuffer.data(kkoff + 100 * idx + 75);

            auto t_yyz_yyy = primBuffer.data(kkoff + 100 * idx + 76);

            auto t_yyz_yyz = primBuffer.data(kkoff + 100 * idx + 77);

            auto t_yyz_yzz = primBuffer.data(kkoff + 100 * idx + 78);

            auto t_yyz_zzz = primBuffer.data(kkoff + 100 * idx + 79);

            auto t_yzz_xxx = primBuffer.data(kkoff + 100 * idx + 80);

            auto t_yzz_xxy = primBuffer.data(kkoff + 100 * idx + 81);

            auto t_yzz_xxz = primBuffer.data(kkoff + 100 * idx + 82);

            auto t_yzz_xyy = primBuffer.data(kkoff + 100 * idx + 83);

            auto t_yzz_xyz = primBuffer.data(kkoff + 100 * idx + 84);

            auto t_yzz_xzz = primBuffer.data(kkoff + 100 * idx + 85);

            auto t_yzz_yyy = primBuffer.data(kkoff + 100 * idx + 86);

            auto t_yzz_yyz = primBuffer.data(kkoff + 100 * idx + 87);

            auto t_yzz_yzz = primBuffer.data(kkoff + 100 * idx + 88);

            auto t_yzz_zzz = primBuffer.data(kkoff + 100 * idx + 89);

            auto t_zzz_xxx = primBuffer.data(kkoff + 100 * idx + 90);

            auto t_zzz_xxy = primBuffer.data(kkoff + 100 * idx + 91);

            auto t_zzz_xxz = primBuffer.data(kkoff + 100 * idx + 92);

            auto t_zzz_xyy = primBuffer.data(kkoff + 100 * idx + 93);

            auto t_zzz_xyz = primBuffer.data(kkoff + 100 * idx + 94);

            auto t_zzz_xzz = primBuffer.data(kkoff + 100 * idx + 95);

            auto t_zzz_yyy = primBuffer.data(kkoff + 100 * idx + 96);

            auto t_zzz_yyz = primBuffer.data(kkoff + 100 * idx + 97);

            auto t_zzz_yzz = primBuffer.data(kkoff + 100 * idx + 98);

            auto t_zzz_zzz = primBuffer.data(kkoff + 100 * idx + 99);

            // set up pointers to (D|G) integrals

            auto s_xx_xxxx = primBuffer.data(s2off + 90 * idx);

            auto s_xx_xxxy = primBuffer.data(s2off + 90 * idx + 1);

            auto s_xx_xxxz = primBuffer.data(s2off + 90 * idx + 2);

            auto s_xx_xxyy = primBuffer.data(s2off + 90 * idx + 3);

            auto s_xx_xxyz = primBuffer.data(s2off + 90 * idx + 4);

            auto s_xx_xxzz = primBuffer.data(s2off + 90 * idx + 5);

            auto s_xx_xyyy = primBuffer.data(s2off + 90 * idx + 6);

            auto s_xx_xyyz = primBuffer.data(s2off + 90 * idx + 7);

            auto s_xx_xyzz = primBuffer.data(s2off + 90 * idx + 8);

            auto s_xx_xzzz = primBuffer.data(s2off + 90 * idx + 9);

            auto s_xx_yyyy = primBuffer.data(s2off + 90 * idx + 10);

            auto s_xx_yyyz = primBuffer.data(s2off + 90 * idx + 11);

            auto s_xx_yyzz = primBuffer.data(s2off + 90 * idx + 12);

            auto s_xx_yzzz = primBuffer.data(s2off + 90 * idx + 13);

            auto s_xx_zzzz = primBuffer.data(s2off + 90 * idx + 14);

            auto s_xy_xxxx = primBuffer.data(s2off + 90 * idx + 15);

            auto s_xy_xxxy = primBuffer.data(s2off + 90 * idx + 16);

            auto s_xy_xxxz = primBuffer.data(s2off + 90 * idx + 17);

            auto s_xy_xxyy = primBuffer.data(s2off + 90 * idx + 18);

            auto s_xy_xxyz = primBuffer.data(s2off + 90 * idx + 19);

            auto s_xy_xxzz = primBuffer.data(s2off + 90 * idx + 20);

            auto s_xy_xyyy = primBuffer.data(s2off + 90 * idx + 21);

            auto s_xy_xyyz = primBuffer.data(s2off + 90 * idx + 22);

            auto s_xy_xyzz = primBuffer.data(s2off + 90 * idx + 23);

            auto s_xy_xzzz = primBuffer.data(s2off + 90 * idx + 24);

            auto s_xy_yyyy = primBuffer.data(s2off + 90 * idx + 25);

            auto s_xy_yyyz = primBuffer.data(s2off + 90 * idx + 26);

            auto s_xy_yyzz = primBuffer.data(s2off + 90 * idx + 27);

            auto s_xy_yzzz = primBuffer.data(s2off + 90 * idx + 28);

            auto s_xy_zzzz = primBuffer.data(s2off + 90 * idx + 29);

            auto s_xz_xxxx = primBuffer.data(s2off + 90 * idx + 30);

            auto s_xz_xxxy = primBuffer.data(s2off + 90 * idx + 31);

            auto s_xz_xxxz = primBuffer.data(s2off + 90 * idx + 32);

            auto s_xz_xxyy = primBuffer.data(s2off + 90 * idx + 33);

            auto s_xz_xxyz = primBuffer.data(s2off + 90 * idx + 34);

            auto s_xz_xxzz = primBuffer.data(s2off + 90 * idx + 35);

            auto s_xz_xyyy = primBuffer.data(s2off + 90 * idx + 36);

            auto s_xz_xyyz = primBuffer.data(s2off + 90 * idx + 37);

            auto s_xz_xyzz = primBuffer.data(s2off + 90 * idx + 38);

            auto s_xz_xzzz = primBuffer.data(s2off + 90 * idx + 39);

            auto s_xz_yyyy = primBuffer.data(s2off + 90 * idx + 40);

            auto s_xz_yyyz = primBuffer.data(s2off + 90 * idx + 41);

            auto s_xz_yyzz = primBuffer.data(s2off + 90 * idx + 42);

            auto s_xz_yzzz = primBuffer.data(s2off + 90 * idx + 43);

            auto s_xz_zzzz = primBuffer.data(s2off + 90 * idx + 44);

            auto s_yy_xxxx = primBuffer.data(s2off + 90 * idx + 45);

            auto s_yy_xxxy = primBuffer.data(s2off + 90 * idx + 46);

            auto s_yy_xxxz = primBuffer.data(s2off + 90 * idx + 47);

            auto s_yy_xxyy = primBuffer.data(s2off + 90 * idx + 48);

            auto s_yy_xxyz = primBuffer.data(s2off + 90 * idx + 49);

            auto s_yy_xxzz = primBuffer.data(s2off + 90 * idx + 50);

            auto s_yy_xyyy = primBuffer.data(s2off + 90 * idx + 51);

            auto s_yy_xyyz = primBuffer.data(s2off + 90 * idx + 52);

            auto s_yy_xyzz = primBuffer.data(s2off + 90 * idx + 53);

            auto s_yy_xzzz = primBuffer.data(s2off + 90 * idx + 54);

            auto s_yy_yyyy = primBuffer.data(s2off + 90 * idx + 55);

            auto s_yy_yyyz = primBuffer.data(s2off + 90 * idx + 56);

            auto s_yy_yyzz = primBuffer.data(s2off + 90 * idx + 57);

            auto s_yy_yzzz = primBuffer.data(s2off + 90 * idx + 58);

            auto s_yy_zzzz = primBuffer.data(s2off + 90 * idx + 59);

            auto s_yz_xxxx = primBuffer.data(s2off + 90 * idx + 60);

            auto s_yz_xxxy = primBuffer.data(s2off + 90 * idx + 61);

            auto s_yz_xxxz = primBuffer.data(s2off + 90 * idx + 62);

            auto s_yz_xxyy = primBuffer.data(s2off + 90 * idx + 63);

            auto s_yz_xxyz = primBuffer.data(s2off + 90 * idx + 64);

            auto s_yz_xxzz = primBuffer.data(s2off + 90 * idx + 65);

            auto s_yz_xyyy = primBuffer.data(s2off + 90 * idx + 66);

            auto s_yz_xyyz = primBuffer.data(s2off + 90 * idx + 67);

            auto s_yz_xyzz = primBuffer.data(s2off + 90 * idx + 68);

            auto s_yz_xzzz = primBuffer.data(s2off + 90 * idx + 69);

            auto s_yz_yyyy = primBuffer.data(s2off + 90 * idx + 70);

            auto s_yz_yyyz = primBuffer.data(s2off + 90 * idx + 71);

            auto s_yz_yyzz = primBuffer.data(s2off + 90 * idx + 72);

            auto s_yz_yzzz = primBuffer.data(s2off + 90 * idx + 73);

            auto s_yz_zzzz = primBuffer.data(s2off + 90 * idx + 74);

            auto s_zz_xxxx = primBuffer.data(s2off + 90 * idx + 75);

            auto s_zz_xxxy = primBuffer.data(s2off + 90 * idx + 76);

            auto s_zz_xxxz = primBuffer.data(s2off + 90 * idx + 77);

            auto s_zz_xxyy = primBuffer.data(s2off + 90 * idx + 78);

            auto s_zz_xxyz = primBuffer.data(s2off + 90 * idx + 79);

            auto s_zz_xxzz = primBuffer.data(s2off + 90 * idx + 80);

            auto s_zz_xyyy = primBuffer.data(s2off + 90 * idx + 81);

            auto s_zz_xyyz = primBuffer.data(s2off + 90 * idx + 82);

            auto s_zz_xyzz = primBuffer.data(s2off + 90 * idx + 83);

            auto s_zz_xzzz = primBuffer.data(s2off + 90 * idx + 84);

            auto s_zz_yyyy = primBuffer.data(s2off + 90 * idx + 85);

            auto s_zz_yyyz = primBuffer.data(s2off + 90 * idx + 86);

            auto s_zz_yyzz = primBuffer.data(s2off + 90 * idx + 87);

            auto s_zz_yzzz = primBuffer.data(s2off + 90 * idx + 88);

            auto s_zz_zzzz = primBuffer.data(s2off + 90 * idx + 89);

            // set up pointers to (D|T|G) integrals

            auto t_xx_xxxx = primBuffer.data(k2off + 90 * idx);

            auto t_xx_xxxy = primBuffer.data(k2off + 90 * idx + 1);

            auto t_xx_xxxz = primBuffer.data(k2off + 90 * idx + 2);

            auto t_xx_xxyy = primBuffer.data(k2off + 90 * idx + 3);

            auto t_xx_xxyz = primBuffer.data(k2off + 90 * idx + 4);

            auto t_xx_xxzz = primBuffer.data(k2off + 90 * idx + 5);

            auto t_xx_xyyy = primBuffer.data(k2off + 90 * idx + 6);

            auto t_xx_xyyz = primBuffer.data(k2off + 90 * idx + 7);

            auto t_xx_xyzz = primBuffer.data(k2off + 90 * idx + 8);

            auto t_xx_xzzz = primBuffer.data(k2off + 90 * idx + 9);

            auto t_xx_yyyy = primBuffer.data(k2off + 90 * idx + 10);

            auto t_xx_yyyz = primBuffer.data(k2off + 90 * idx + 11);

            auto t_xx_yyzz = primBuffer.data(k2off + 90 * idx + 12);

            auto t_xx_yzzz = primBuffer.data(k2off + 90 * idx + 13);

            auto t_xx_zzzz = primBuffer.data(k2off + 90 * idx + 14);

            auto t_xy_xxxx = primBuffer.data(k2off + 90 * idx + 15);

            auto t_xy_xxxy = primBuffer.data(k2off + 90 * idx + 16);

            auto t_xy_xxxz = primBuffer.data(k2off + 90 * idx + 17);

            auto t_xy_xxyy = primBuffer.data(k2off + 90 * idx + 18);

            auto t_xy_xxyz = primBuffer.data(k2off + 90 * idx + 19);

            auto t_xy_xxzz = primBuffer.data(k2off + 90 * idx + 20);

            auto t_xy_xyyy = primBuffer.data(k2off + 90 * idx + 21);

            auto t_xy_xyyz = primBuffer.data(k2off + 90 * idx + 22);

            auto t_xy_xyzz = primBuffer.data(k2off + 90 * idx + 23);

            auto t_xy_xzzz = primBuffer.data(k2off + 90 * idx + 24);

            auto t_xy_yyyy = primBuffer.data(k2off + 90 * idx + 25);

            auto t_xy_yyyz = primBuffer.data(k2off + 90 * idx + 26);

            auto t_xy_yyzz = primBuffer.data(k2off + 90 * idx + 27);

            auto t_xy_yzzz = primBuffer.data(k2off + 90 * idx + 28);

            auto t_xy_zzzz = primBuffer.data(k2off + 90 * idx + 29);

            auto t_xz_xxxx = primBuffer.data(k2off + 90 * idx + 30);

            auto t_xz_xxxy = primBuffer.data(k2off + 90 * idx + 31);

            auto t_xz_xxxz = primBuffer.data(k2off + 90 * idx + 32);

            auto t_xz_xxyy = primBuffer.data(k2off + 90 * idx + 33);

            auto t_xz_xxyz = primBuffer.data(k2off + 90 * idx + 34);

            auto t_xz_xxzz = primBuffer.data(k2off + 90 * idx + 35);

            auto t_xz_xyyy = primBuffer.data(k2off + 90 * idx + 36);

            auto t_xz_xyyz = primBuffer.data(k2off + 90 * idx + 37);

            auto t_xz_xyzz = primBuffer.data(k2off + 90 * idx + 38);

            auto t_xz_xzzz = primBuffer.data(k2off + 90 * idx + 39);

            auto t_xz_yyyy = primBuffer.data(k2off + 90 * idx + 40);

            auto t_xz_yyyz = primBuffer.data(k2off + 90 * idx + 41);

            auto t_xz_yyzz = primBuffer.data(k2off + 90 * idx + 42);

            auto t_xz_yzzz = primBuffer.data(k2off + 90 * idx + 43);

            auto t_xz_zzzz = primBuffer.data(k2off + 90 * idx + 44);

            auto t_yy_xxxx = primBuffer.data(k2off + 90 * idx + 45);

            auto t_yy_xxxy = primBuffer.data(k2off + 90 * idx + 46);

            auto t_yy_xxxz = primBuffer.data(k2off + 90 * idx + 47);

            auto t_yy_xxyy = primBuffer.data(k2off + 90 * idx + 48);

            auto t_yy_xxyz = primBuffer.data(k2off + 90 * idx + 49);

            auto t_yy_xxzz = primBuffer.data(k2off + 90 * idx + 50);

            auto t_yy_xyyy = primBuffer.data(k2off + 90 * idx + 51);

            auto t_yy_xyyz = primBuffer.data(k2off + 90 * idx + 52);

            auto t_yy_xyzz = primBuffer.data(k2off + 90 * idx + 53);

            auto t_yy_xzzz = primBuffer.data(k2off + 90 * idx + 54);

            auto t_yy_yyyy = primBuffer.data(k2off + 90 * idx + 55);

            auto t_yy_yyyz = primBuffer.data(k2off + 90 * idx + 56);

            auto t_yy_yyzz = primBuffer.data(k2off + 90 * idx + 57);

            auto t_yy_yzzz = primBuffer.data(k2off + 90 * idx + 58);

            auto t_yy_zzzz = primBuffer.data(k2off + 90 * idx + 59);

            auto t_yz_xxxx = primBuffer.data(k2off + 90 * idx + 60);

            auto t_yz_xxxy = primBuffer.data(k2off + 90 * idx + 61);

            auto t_yz_xxxz = primBuffer.data(k2off + 90 * idx + 62);

            auto t_yz_xxyy = primBuffer.data(k2off + 90 * idx + 63);

            auto t_yz_xxyz = primBuffer.data(k2off + 90 * idx + 64);

            auto t_yz_xxzz = primBuffer.data(k2off + 90 * idx + 65);

            auto t_yz_xyyy = primBuffer.data(k2off + 90 * idx + 66);

            auto t_yz_xyyz = primBuffer.data(k2off + 90 * idx + 67);

            auto t_yz_xyzz = primBuffer.data(k2off + 90 * idx + 68);

            auto t_yz_xzzz = primBuffer.data(k2off + 90 * idx + 69);

            auto t_yz_yyyy = primBuffer.data(k2off + 90 * idx + 70);

            auto t_yz_yyyz = primBuffer.data(k2off + 90 * idx + 71);

            auto t_yz_yyzz = primBuffer.data(k2off + 90 * idx + 72);

            auto t_yz_yzzz = primBuffer.data(k2off + 90 * idx + 73);

            auto t_yz_zzzz = primBuffer.data(k2off + 90 * idx + 74);

            auto t_zz_xxxx = primBuffer.data(k2off + 90 * idx + 75);

            auto t_zz_xxxy = primBuffer.data(k2off + 90 * idx + 76);

            auto t_zz_xxxz = primBuffer.data(k2off + 90 * idx + 77);

            auto t_zz_xxyy = primBuffer.data(k2off + 90 * idx + 78);

            auto t_zz_xxyz = primBuffer.data(k2off + 90 * idx + 79);

            auto t_zz_xxzz = primBuffer.data(k2off + 90 * idx + 80);

            auto t_zz_xyyy = primBuffer.data(k2off + 90 * idx + 81);

            auto t_zz_xyyz = primBuffer.data(k2off + 90 * idx + 82);

            auto t_zz_xyzz = primBuffer.data(k2off + 90 * idx + 83);

            auto t_zz_xzzz = primBuffer.data(k2off + 90 * idx + 84);

            auto t_zz_yyyy = primBuffer.data(k2off + 90 * idx + 85);

            auto t_zz_yyyz = primBuffer.data(k2off + 90 * idx + 86);

            auto t_zz_yyzz = primBuffer.data(k2off + 90 * idx + 87);

            auto t_zz_yzzz = primBuffer.data(k2off + 90 * idx + 88);

            auto t_zz_zzzz = primBuffer.data(k2off + 90 * idx + 89);

            // set up pointers to (F|G) integrals

            auto s_xxx_xxxx = primBuffer.data(s1off + 150 * idx);

            auto s_xxx_xxxy = primBuffer.data(s1off + 150 * idx + 1);

            auto s_xxx_xxxz = primBuffer.data(s1off + 150 * idx + 2);

            auto s_xxx_xxyy = primBuffer.data(s1off + 150 * idx + 3);

            auto s_xxx_xxyz = primBuffer.data(s1off + 150 * idx + 4);

            auto s_xxx_xxzz = primBuffer.data(s1off + 150 * idx + 5);

            auto s_xxx_xyyy = primBuffer.data(s1off + 150 * idx + 6);

            auto s_xxx_xyyz = primBuffer.data(s1off + 150 * idx + 7);

            auto s_xxx_xyzz = primBuffer.data(s1off + 150 * idx + 8);

            auto s_xxx_xzzz = primBuffer.data(s1off + 150 * idx + 9);

            auto s_xxx_yyyy = primBuffer.data(s1off + 150 * idx + 10);

            auto s_xxx_yyyz = primBuffer.data(s1off + 150 * idx + 11);

            auto s_xxx_yyzz = primBuffer.data(s1off + 150 * idx + 12);

            auto s_xxx_yzzz = primBuffer.data(s1off + 150 * idx + 13);

            auto s_xxx_zzzz = primBuffer.data(s1off + 150 * idx + 14);

            auto s_xxy_xxxx = primBuffer.data(s1off + 150 * idx + 15);

            auto s_xxy_xxxy = primBuffer.data(s1off + 150 * idx + 16);

            auto s_xxy_xxxz = primBuffer.data(s1off + 150 * idx + 17);

            auto s_xxy_xxyy = primBuffer.data(s1off + 150 * idx + 18);

            auto s_xxy_xxyz = primBuffer.data(s1off + 150 * idx + 19);

            auto s_xxy_xxzz = primBuffer.data(s1off + 150 * idx + 20);

            auto s_xxy_xyyy = primBuffer.data(s1off + 150 * idx + 21);

            auto s_xxy_xyyz = primBuffer.data(s1off + 150 * idx + 22);

            auto s_xxy_xyzz = primBuffer.data(s1off + 150 * idx + 23);

            auto s_xxy_xzzz = primBuffer.data(s1off + 150 * idx + 24);

            auto s_xxy_yyyy = primBuffer.data(s1off + 150 * idx + 25);

            auto s_xxy_yyyz = primBuffer.data(s1off + 150 * idx + 26);

            auto s_xxy_yyzz = primBuffer.data(s1off + 150 * idx + 27);

            auto s_xxy_yzzz = primBuffer.data(s1off + 150 * idx + 28);

            auto s_xxy_zzzz = primBuffer.data(s1off + 150 * idx + 29);

            auto s_xxz_xxxx = primBuffer.data(s1off + 150 * idx + 30);

            auto s_xxz_xxxy = primBuffer.data(s1off + 150 * idx + 31);

            auto s_xxz_xxxz = primBuffer.data(s1off + 150 * idx + 32);

            auto s_xxz_xxyy = primBuffer.data(s1off + 150 * idx + 33);

            auto s_xxz_xxyz = primBuffer.data(s1off + 150 * idx + 34);

            auto s_xxz_xxzz = primBuffer.data(s1off + 150 * idx + 35);

            auto s_xxz_xyyy = primBuffer.data(s1off + 150 * idx + 36);

            auto s_xxz_xyyz = primBuffer.data(s1off + 150 * idx + 37);

            auto s_xxz_xyzz = primBuffer.data(s1off + 150 * idx + 38);

            auto s_xxz_xzzz = primBuffer.data(s1off + 150 * idx + 39);

            auto s_xxz_yyyy = primBuffer.data(s1off + 150 * idx + 40);

            auto s_xxz_yyyz = primBuffer.data(s1off + 150 * idx + 41);

            auto s_xxz_yyzz = primBuffer.data(s1off + 150 * idx + 42);

            auto s_xxz_yzzz = primBuffer.data(s1off + 150 * idx + 43);

            auto s_xxz_zzzz = primBuffer.data(s1off + 150 * idx + 44);

            auto s_xyy_xxxx = primBuffer.data(s1off + 150 * idx + 45);

            auto s_xyy_xxxy = primBuffer.data(s1off + 150 * idx + 46);

            auto s_xyy_xxxz = primBuffer.data(s1off + 150 * idx + 47);

            auto s_xyy_xxyy = primBuffer.data(s1off + 150 * idx + 48);

            auto s_xyy_xxyz = primBuffer.data(s1off + 150 * idx + 49);

            auto s_xyy_xxzz = primBuffer.data(s1off + 150 * idx + 50);

            auto s_xyy_xyyy = primBuffer.data(s1off + 150 * idx + 51);

            auto s_xyy_xyyz = primBuffer.data(s1off + 150 * idx + 52);

            auto s_xyy_xyzz = primBuffer.data(s1off + 150 * idx + 53);

            auto s_xyy_xzzz = primBuffer.data(s1off + 150 * idx + 54);

            auto s_xyy_yyyy = primBuffer.data(s1off + 150 * idx + 55);

            auto s_xyy_yyyz = primBuffer.data(s1off + 150 * idx + 56);

            auto s_xyy_yyzz = primBuffer.data(s1off + 150 * idx + 57);

            auto s_xyy_yzzz = primBuffer.data(s1off + 150 * idx + 58);

            auto s_xyy_zzzz = primBuffer.data(s1off + 150 * idx + 59);

            auto s_xyz_xxxx = primBuffer.data(s1off + 150 * idx + 60);

            auto s_xyz_xxxy = primBuffer.data(s1off + 150 * idx + 61);

            auto s_xyz_xxxz = primBuffer.data(s1off + 150 * idx + 62);

            auto s_xyz_xxyy = primBuffer.data(s1off + 150 * idx + 63);

            auto s_xyz_xxyz = primBuffer.data(s1off + 150 * idx + 64);

            auto s_xyz_xxzz = primBuffer.data(s1off + 150 * idx + 65);

            auto s_xyz_xyyy = primBuffer.data(s1off + 150 * idx + 66);

            auto s_xyz_xyyz = primBuffer.data(s1off + 150 * idx + 67);

            auto s_xyz_xyzz = primBuffer.data(s1off + 150 * idx + 68);

            auto s_xyz_xzzz = primBuffer.data(s1off + 150 * idx + 69);

            auto s_xyz_yyyy = primBuffer.data(s1off + 150 * idx + 70);

            auto s_xyz_yyyz = primBuffer.data(s1off + 150 * idx + 71);

            auto s_xyz_yyzz = primBuffer.data(s1off + 150 * idx + 72);

            auto s_xyz_yzzz = primBuffer.data(s1off + 150 * idx + 73);

            auto s_xyz_zzzz = primBuffer.data(s1off + 150 * idx + 74);

            auto s_xzz_xxxx = primBuffer.data(s1off + 150 * idx + 75);

            auto s_xzz_xxxy = primBuffer.data(s1off + 150 * idx + 76);

            auto s_xzz_xxxz = primBuffer.data(s1off + 150 * idx + 77);

            auto s_xzz_xxyy = primBuffer.data(s1off + 150 * idx + 78);

            auto s_xzz_xxyz = primBuffer.data(s1off + 150 * idx + 79);

            auto s_xzz_xxzz = primBuffer.data(s1off + 150 * idx + 80);

            auto s_xzz_xyyy = primBuffer.data(s1off + 150 * idx + 81);

            auto s_xzz_xyyz = primBuffer.data(s1off + 150 * idx + 82);

            auto s_xzz_xyzz = primBuffer.data(s1off + 150 * idx + 83);

            auto s_xzz_xzzz = primBuffer.data(s1off + 150 * idx + 84);

            auto s_xzz_yyyy = primBuffer.data(s1off + 150 * idx + 85);

            auto s_xzz_yyyz = primBuffer.data(s1off + 150 * idx + 86);

            auto s_xzz_yyzz = primBuffer.data(s1off + 150 * idx + 87);

            auto s_xzz_yzzz = primBuffer.data(s1off + 150 * idx + 88);

            auto s_xzz_zzzz = primBuffer.data(s1off + 150 * idx + 89);

            auto s_yyy_xxxx = primBuffer.data(s1off + 150 * idx + 90);

            auto s_yyy_xxxy = primBuffer.data(s1off + 150 * idx + 91);

            auto s_yyy_xxxz = primBuffer.data(s1off + 150 * idx + 92);

            auto s_yyy_xxyy = primBuffer.data(s1off + 150 * idx + 93);

            auto s_yyy_xxyz = primBuffer.data(s1off + 150 * idx + 94);

            auto s_yyy_xxzz = primBuffer.data(s1off + 150 * idx + 95);

            auto s_yyy_xyyy = primBuffer.data(s1off + 150 * idx + 96);

            auto s_yyy_xyyz = primBuffer.data(s1off + 150 * idx + 97);

            auto s_yyy_xyzz = primBuffer.data(s1off + 150 * idx + 98);

            auto s_yyy_xzzz = primBuffer.data(s1off + 150 * idx + 99);

            auto s_yyy_yyyy = primBuffer.data(s1off + 150 * idx + 100);

            auto s_yyy_yyyz = primBuffer.data(s1off + 150 * idx + 101);

            auto s_yyy_yyzz = primBuffer.data(s1off + 150 * idx + 102);

            auto s_yyy_yzzz = primBuffer.data(s1off + 150 * idx + 103);

            auto s_yyy_zzzz = primBuffer.data(s1off + 150 * idx + 104);

            auto s_yyz_xxxx = primBuffer.data(s1off + 150 * idx + 105);

            auto s_yyz_xxxy = primBuffer.data(s1off + 150 * idx + 106);

            auto s_yyz_xxxz = primBuffer.data(s1off + 150 * idx + 107);

            auto s_yyz_xxyy = primBuffer.data(s1off + 150 * idx + 108);

            auto s_yyz_xxyz = primBuffer.data(s1off + 150 * idx + 109);

            auto s_yyz_xxzz = primBuffer.data(s1off + 150 * idx + 110);

            auto s_yyz_xyyy = primBuffer.data(s1off + 150 * idx + 111);

            auto s_yyz_xyyz = primBuffer.data(s1off + 150 * idx + 112);

            auto s_yyz_xyzz = primBuffer.data(s1off + 150 * idx + 113);

            auto s_yyz_xzzz = primBuffer.data(s1off + 150 * idx + 114);

            auto s_yyz_yyyy = primBuffer.data(s1off + 150 * idx + 115);

            auto s_yyz_yyyz = primBuffer.data(s1off + 150 * idx + 116);

            auto s_yyz_yyzz = primBuffer.data(s1off + 150 * idx + 117);

            auto s_yyz_yzzz = primBuffer.data(s1off + 150 * idx + 118);

            auto s_yyz_zzzz = primBuffer.data(s1off + 150 * idx + 119);

            auto s_yzz_xxxx = primBuffer.data(s1off + 150 * idx + 120);

            auto s_yzz_xxxy = primBuffer.data(s1off + 150 * idx + 121);

            auto s_yzz_xxxz = primBuffer.data(s1off + 150 * idx + 122);

            auto s_yzz_xxyy = primBuffer.data(s1off + 150 * idx + 123);

            auto s_yzz_xxyz = primBuffer.data(s1off + 150 * idx + 124);

            auto s_yzz_xxzz = primBuffer.data(s1off + 150 * idx + 125);

            auto s_yzz_xyyy = primBuffer.data(s1off + 150 * idx + 126);

            auto s_yzz_xyyz = primBuffer.data(s1off + 150 * idx + 127);

            auto s_yzz_xyzz = primBuffer.data(s1off + 150 * idx + 128);

            auto s_yzz_xzzz = primBuffer.data(s1off + 150 * idx + 129);

            auto s_yzz_yyyy = primBuffer.data(s1off + 150 * idx + 130);

            auto s_yzz_yyyz = primBuffer.data(s1off + 150 * idx + 131);

            auto s_yzz_yyzz = primBuffer.data(s1off + 150 * idx + 132);

            auto s_yzz_yzzz = primBuffer.data(s1off + 150 * idx + 133);

            auto s_yzz_zzzz = primBuffer.data(s1off + 150 * idx + 134);

            auto s_zzz_xxxx = primBuffer.data(s1off + 150 * idx + 135);

            auto s_zzz_xxxy = primBuffer.data(s1off + 150 * idx + 136);

            auto s_zzz_xxxz = primBuffer.data(s1off + 150 * idx + 137);

            auto s_zzz_xxyy = primBuffer.data(s1off + 150 * idx + 138);

            auto s_zzz_xxyz = primBuffer.data(s1off + 150 * idx + 139);

            auto s_zzz_xxzz = primBuffer.data(s1off + 150 * idx + 140);

            auto s_zzz_xyyy = primBuffer.data(s1off + 150 * idx + 141);

            auto s_zzz_xyyz = primBuffer.data(s1off + 150 * idx + 142);

            auto s_zzz_xyzz = primBuffer.data(s1off + 150 * idx + 143);

            auto s_zzz_xzzz = primBuffer.data(s1off + 150 * idx + 144);

            auto s_zzz_yyyy = primBuffer.data(s1off + 150 * idx + 145);

            auto s_zzz_yyyz = primBuffer.data(s1off + 150 * idx + 146);

            auto s_zzz_yyzz = primBuffer.data(s1off + 150 * idx + 147);

            auto s_zzz_yzzz = primBuffer.data(s1off + 150 * idx + 148);

            auto s_zzz_zzzz = primBuffer.data(s1off + 150 * idx + 149);

            // set up pointers to (F|T|G) integrals

            auto t_xxx_xxxx = primBuffer.data(k1off + 150 * idx);

            auto t_xxx_xxxy = primBuffer.data(k1off + 150 * idx + 1);

            auto t_xxx_xxxz = primBuffer.data(k1off + 150 * idx + 2);

            auto t_xxx_xxyy = primBuffer.data(k1off + 150 * idx + 3);

            auto t_xxx_xxyz = primBuffer.data(k1off + 150 * idx + 4);

            auto t_xxx_xxzz = primBuffer.data(k1off + 150 * idx + 5);

            auto t_xxx_xyyy = primBuffer.data(k1off + 150 * idx + 6);

            auto t_xxx_xyyz = primBuffer.data(k1off + 150 * idx + 7);

            auto t_xxx_xyzz = primBuffer.data(k1off + 150 * idx + 8);

            auto t_xxx_xzzz = primBuffer.data(k1off + 150 * idx + 9);

            auto t_xxx_yyyy = primBuffer.data(k1off + 150 * idx + 10);

            auto t_xxx_yyyz = primBuffer.data(k1off + 150 * idx + 11);

            auto t_xxx_yyzz = primBuffer.data(k1off + 150 * idx + 12);

            auto t_xxx_yzzz = primBuffer.data(k1off + 150 * idx + 13);

            auto t_xxx_zzzz = primBuffer.data(k1off + 150 * idx + 14);

            auto t_xxy_xxxx = primBuffer.data(k1off + 150 * idx + 15);

            auto t_xxy_xxxy = primBuffer.data(k1off + 150 * idx + 16);

            auto t_xxy_xxxz = primBuffer.data(k1off + 150 * idx + 17);

            auto t_xxy_xxyy = primBuffer.data(k1off + 150 * idx + 18);

            auto t_xxy_xxyz = primBuffer.data(k1off + 150 * idx + 19);

            auto t_xxy_xxzz = primBuffer.data(k1off + 150 * idx + 20);

            auto t_xxy_xyyy = primBuffer.data(k1off + 150 * idx + 21);

            auto t_xxy_xyyz = primBuffer.data(k1off + 150 * idx + 22);

            auto t_xxy_xyzz = primBuffer.data(k1off + 150 * idx + 23);

            auto t_xxy_xzzz = primBuffer.data(k1off + 150 * idx + 24);

            auto t_xxy_yyyy = primBuffer.data(k1off + 150 * idx + 25);

            auto t_xxy_yyyz = primBuffer.data(k1off + 150 * idx + 26);

            auto t_xxy_yyzz = primBuffer.data(k1off + 150 * idx + 27);

            auto t_xxy_yzzz = primBuffer.data(k1off + 150 * idx + 28);

            auto t_xxy_zzzz = primBuffer.data(k1off + 150 * idx + 29);

            auto t_xxz_xxxx = primBuffer.data(k1off + 150 * idx + 30);

            auto t_xxz_xxxy = primBuffer.data(k1off + 150 * idx + 31);

            auto t_xxz_xxxz = primBuffer.data(k1off + 150 * idx + 32);

            auto t_xxz_xxyy = primBuffer.data(k1off + 150 * idx + 33);

            auto t_xxz_xxyz = primBuffer.data(k1off + 150 * idx + 34);

            auto t_xxz_xxzz = primBuffer.data(k1off + 150 * idx + 35);

            auto t_xxz_xyyy = primBuffer.data(k1off + 150 * idx + 36);

            auto t_xxz_xyyz = primBuffer.data(k1off + 150 * idx + 37);

            auto t_xxz_xyzz = primBuffer.data(k1off + 150 * idx + 38);

            auto t_xxz_xzzz = primBuffer.data(k1off + 150 * idx + 39);

            auto t_xxz_yyyy = primBuffer.data(k1off + 150 * idx + 40);

            auto t_xxz_yyyz = primBuffer.data(k1off + 150 * idx + 41);

            auto t_xxz_yyzz = primBuffer.data(k1off + 150 * idx + 42);

            auto t_xxz_yzzz = primBuffer.data(k1off + 150 * idx + 43);

            auto t_xxz_zzzz = primBuffer.data(k1off + 150 * idx + 44);

            auto t_xyy_xxxx = primBuffer.data(k1off + 150 * idx + 45);

            auto t_xyy_xxxy = primBuffer.data(k1off + 150 * idx + 46);

            auto t_xyy_xxxz = primBuffer.data(k1off + 150 * idx + 47);

            auto t_xyy_xxyy = primBuffer.data(k1off + 150 * idx + 48);

            auto t_xyy_xxyz = primBuffer.data(k1off + 150 * idx + 49);

            auto t_xyy_xxzz = primBuffer.data(k1off + 150 * idx + 50);

            auto t_xyy_xyyy = primBuffer.data(k1off + 150 * idx + 51);

            auto t_xyy_xyyz = primBuffer.data(k1off + 150 * idx + 52);

            auto t_xyy_xyzz = primBuffer.data(k1off + 150 * idx + 53);

            auto t_xyy_xzzz = primBuffer.data(k1off + 150 * idx + 54);

            auto t_xyy_yyyy = primBuffer.data(k1off + 150 * idx + 55);

            auto t_xyy_yyyz = primBuffer.data(k1off + 150 * idx + 56);

            auto t_xyy_yyzz = primBuffer.data(k1off + 150 * idx + 57);

            auto t_xyy_yzzz = primBuffer.data(k1off + 150 * idx + 58);

            auto t_xyy_zzzz = primBuffer.data(k1off + 150 * idx + 59);

            auto t_xyz_xxxx = primBuffer.data(k1off + 150 * idx + 60);

            auto t_xyz_xxxy = primBuffer.data(k1off + 150 * idx + 61);

            auto t_xyz_xxxz = primBuffer.data(k1off + 150 * idx + 62);

            auto t_xyz_xxyy = primBuffer.data(k1off + 150 * idx + 63);

            auto t_xyz_xxyz = primBuffer.data(k1off + 150 * idx + 64);

            auto t_xyz_xxzz = primBuffer.data(k1off + 150 * idx + 65);

            auto t_xyz_xyyy = primBuffer.data(k1off + 150 * idx + 66);

            auto t_xyz_xyyz = primBuffer.data(k1off + 150 * idx + 67);

            auto t_xyz_xyzz = primBuffer.data(k1off + 150 * idx + 68);

            auto t_xyz_xzzz = primBuffer.data(k1off + 150 * idx + 69);

            auto t_xyz_yyyy = primBuffer.data(k1off + 150 * idx + 70);

            auto t_xyz_yyyz = primBuffer.data(k1off + 150 * idx + 71);

            auto t_xyz_yyzz = primBuffer.data(k1off + 150 * idx + 72);

            auto t_xyz_yzzz = primBuffer.data(k1off + 150 * idx + 73);

            auto t_xyz_zzzz = primBuffer.data(k1off + 150 * idx + 74);

            auto t_xzz_xxxx = primBuffer.data(k1off + 150 * idx + 75);

            auto t_xzz_xxxy = primBuffer.data(k1off + 150 * idx + 76);

            auto t_xzz_xxxz = primBuffer.data(k1off + 150 * idx + 77);

            auto t_xzz_xxyy = primBuffer.data(k1off + 150 * idx + 78);

            auto t_xzz_xxyz = primBuffer.data(k1off + 150 * idx + 79);

            auto t_xzz_xxzz = primBuffer.data(k1off + 150 * idx + 80);

            auto t_xzz_xyyy = primBuffer.data(k1off + 150 * idx + 81);

            auto t_xzz_xyyz = primBuffer.data(k1off + 150 * idx + 82);

            auto t_xzz_xyzz = primBuffer.data(k1off + 150 * idx + 83);

            auto t_xzz_xzzz = primBuffer.data(k1off + 150 * idx + 84);

            auto t_xzz_yyyy = primBuffer.data(k1off + 150 * idx + 85);

            auto t_xzz_yyyz = primBuffer.data(k1off + 150 * idx + 86);

            auto t_xzz_yyzz = primBuffer.data(k1off + 150 * idx + 87);

            auto t_xzz_yzzz = primBuffer.data(k1off + 150 * idx + 88);

            auto t_xzz_zzzz = primBuffer.data(k1off + 150 * idx + 89);

            auto t_yyy_xxxx = primBuffer.data(k1off + 150 * idx + 90);

            auto t_yyy_xxxy = primBuffer.data(k1off + 150 * idx + 91);

            auto t_yyy_xxxz = primBuffer.data(k1off + 150 * idx + 92);

            auto t_yyy_xxyy = primBuffer.data(k1off + 150 * idx + 93);

            auto t_yyy_xxyz = primBuffer.data(k1off + 150 * idx + 94);

            auto t_yyy_xxzz = primBuffer.data(k1off + 150 * idx + 95);

            auto t_yyy_xyyy = primBuffer.data(k1off + 150 * idx + 96);

            auto t_yyy_xyyz = primBuffer.data(k1off + 150 * idx + 97);

            auto t_yyy_xyzz = primBuffer.data(k1off + 150 * idx + 98);

            auto t_yyy_xzzz = primBuffer.data(k1off + 150 * idx + 99);

            auto t_yyy_yyyy = primBuffer.data(k1off + 150 * idx + 100);

            auto t_yyy_yyyz = primBuffer.data(k1off + 150 * idx + 101);

            auto t_yyy_yyzz = primBuffer.data(k1off + 150 * idx + 102);

            auto t_yyy_yzzz = primBuffer.data(k1off + 150 * idx + 103);

            auto t_yyy_zzzz = primBuffer.data(k1off + 150 * idx + 104);

            auto t_yyz_xxxx = primBuffer.data(k1off + 150 * idx + 105);

            auto t_yyz_xxxy = primBuffer.data(k1off + 150 * idx + 106);

            auto t_yyz_xxxz = primBuffer.data(k1off + 150 * idx + 107);

            auto t_yyz_xxyy = primBuffer.data(k1off + 150 * idx + 108);

            auto t_yyz_xxyz = primBuffer.data(k1off + 150 * idx + 109);

            auto t_yyz_xxzz = primBuffer.data(k1off + 150 * idx + 110);

            auto t_yyz_xyyy = primBuffer.data(k1off + 150 * idx + 111);

            auto t_yyz_xyyz = primBuffer.data(k1off + 150 * idx + 112);

            auto t_yyz_xyzz = primBuffer.data(k1off + 150 * idx + 113);

            auto t_yyz_xzzz = primBuffer.data(k1off + 150 * idx + 114);

            auto t_yyz_yyyy = primBuffer.data(k1off + 150 * idx + 115);

            auto t_yyz_yyyz = primBuffer.data(k1off + 150 * idx + 116);

            auto t_yyz_yyzz = primBuffer.data(k1off + 150 * idx + 117);

            auto t_yyz_yzzz = primBuffer.data(k1off + 150 * idx + 118);

            auto t_yyz_zzzz = primBuffer.data(k1off + 150 * idx + 119);

            auto t_yzz_xxxx = primBuffer.data(k1off + 150 * idx + 120);

            auto t_yzz_xxxy = primBuffer.data(k1off + 150 * idx + 121);

            auto t_yzz_xxxz = primBuffer.data(k1off + 150 * idx + 122);

            auto t_yzz_xxyy = primBuffer.data(k1off + 150 * idx + 123);

            auto t_yzz_xxyz = primBuffer.data(k1off + 150 * idx + 124);

            auto t_yzz_xxzz = primBuffer.data(k1off + 150 * idx + 125);

            auto t_yzz_xyyy = primBuffer.data(k1off + 150 * idx + 126);

            auto t_yzz_xyyz = primBuffer.data(k1off + 150 * idx + 127);

            auto t_yzz_xyzz = primBuffer.data(k1off + 150 * idx + 128);

            auto t_yzz_xzzz = primBuffer.data(k1off + 150 * idx + 129);

            auto t_yzz_yyyy = primBuffer.data(k1off + 150 * idx + 130);

            auto t_yzz_yyyz = primBuffer.data(k1off + 150 * idx + 131);

            auto t_yzz_yyzz = primBuffer.data(k1off + 150 * idx + 132);

            auto t_yzz_yzzz = primBuffer.data(k1off + 150 * idx + 133);

            auto t_yzz_zzzz = primBuffer.data(k1off + 150 * idx + 134);

            auto t_zzz_xxxx = primBuffer.data(k1off + 150 * idx + 135);

            auto t_zzz_xxxy = primBuffer.data(k1off + 150 * idx + 136);

            auto t_zzz_xxxz = primBuffer.data(k1off + 150 * idx + 137);

            auto t_zzz_xxyy = primBuffer.data(k1off + 150 * idx + 138);

            auto t_zzz_xxyz = primBuffer.data(k1off + 150 * idx + 139);

            auto t_zzz_xxzz = primBuffer.data(k1off + 150 * idx + 140);

            auto t_zzz_xyyy = primBuffer.data(k1off + 150 * idx + 141);

            auto t_zzz_xyyz = primBuffer.data(k1off + 150 * idx + 142);

            auto t_zzz_xyzz = primBuffer.data(k1off + 150 * idx + 143);

            auto t_zzz_xzzz = primBuffer.data(k1off + 150 * idx + 144);

            auto t_zzz_yyyy = primBuffer.data(k1off + 150 * idx + 145);

            auto t_zzz_yyyz = primBuffer.data(k1off + 150 * idx + 146);

            auto t_zzz_yyzz = primBuffer.data(k1off + 150 * idx + 147);

            auto t_zzz_yzzz = primBuffer.data(k1off + 150 * idx + 148);

            auto t_zzz_zzzz = primBuffer.data(k1off + 150 * idx + 149);

            // set up pointers to (G|G) integrals

            auto s_xxxx_xxxx = primBuffer.data(soff + 225 * idx);

            auto s_xxxx_xxxy = primBuffer.data(soff + 225 * idx + 1);

            auto s_xxxx_xxxz = primBuffer.data(soff + 225 * idx + 2);

            auto s_xxxx_xxyy = primBuffer.data(soff + 225 * idx + 3);

            auto s_xxxx_xxyz = primBuffer.data(soff + 225 * idx + 4);

            auto s_xxxx_xxzz = primBuffer.data(soff + 225 * idx + 5);

            auto s_xxxx_xyyy = primBuffer.data(soff + 225 * idx + 6);

            auto s_xxxx_xyyz = primBuffer.data(soff + 225 * idx + 7);

            auto s_xxxx_xyzz = primBuffer.data(soff + 225 * idx + 8);

            auto s_xxxx_xzzz = primBuffer.data(soff + 225 * idx + 9);

            auto s_xxxx_yyyy = primBuffer.data(soff + 225 * idx + 10);

            auto s_xxxx_yyyz = primBuffer.data(soff + 225 * idx + 11);

            auto s_xxxx_yyzz = primBuffer.data(soff + 225 * idx + 12);

            auto s_xxxx_yzzz = primBuffer.data(soff + 225 * idx + 13);

            auto s_xxxx_zzzz = primBuffer.data(soff + 225 * idx + 14);

            auto s_xxxy_xxxx = primBuffer.data(soff + 225 * idx + 15);

            auto s_xxxy_xxxy = primBuffer.data(soff + 225 * idx + 16);

            auto s_xxxy_xxxz = primBuffer.data(soff + 225 * idx + 17);

            auto s_xxxy_xxyy = primBuffer.data(soff + 225 * idx + 18);

            auto s_xxxy_xxyz = primBuffer.data(soff + 225 * idx + 19);

            auto s_xxxy_xxzz = primBuffer.data(soff + 225 * idx + 20);

            auto s_xxxy_xyyy = primBuffer.data(soff + 225 * idx + 21);

            auto s_xxxy_xyyz = primBuffer.data(soff + 225 * idx + 22);

            auto s_xxxy_xyzz = primBuffer.data(soff + 225 * idx + 23);

            auto s_xxxy_xzzz = primBuffer.data(soff + 225 * idx + 24);

            auto s_xxxy_yyyy = primBuffer.data(soff + 225 * idx + 25);

            auto s_xxxy_yyyz = primBuffer.data(soff + 225 * idx + 26);

            auto s_xxxy_yyzz = primBuffer.data(soff + 225 * idx + 27);

            auto s_xxxy_yzzz = primBuffer.data(soff + 225 * idx + 28);

            auto s_xxxy_zzzz = primBuffer.data(soff + 225 * idx + 29);

            auto s_xxxz_xxxx = primBuffer.data(soff + 225 * idx + 30);

            auto s_xxxz_xxxy = primBuffer.data(soff + 225 * idx + 31);

            auto s_xxxz_xxxz = primBuffer.data(soff + 225 * idx + 32);

            auto s_xxxz_xxyy = primBuffer.data(soff + 225 * idx + 33);

            auto s_xxxz_xxyz = primBuffer.data(soff + 225 * idx + 34);

            auto s_xxxz_xxzz = primBuffer.data(soff + 225 * idx + 35);

            auto s_xxxz_xyyy = primBuffer.data(soff + 225 * idx + 36);

            auto s_xxxz_xyyz = primBuffer.data(soff + 225 * idx + 37);

            auto s_xxxz_xyzz = primBuffer.data(soff + 225 * idx + 38);

            auto s_xxxz_xzzz = primBuffer.data(soff + 225 * idx + 39);

            auto s_xxxz_yyyy = primBuffer.data(soff + 225 * idx + 40);

            auto s_xxxz_yyyz = primBuffer.data(soff + 225 * idx + 41);

            auto s_xxxz_yyzz = primBuffer.data(soff + 225 * idx + 42);

            auto s_xxxz_yzzz = primBuffer.data(soff + 225 * idx + 43);

            auto s_xxxz_zzzz = primBuffer.data(soff + 225 * idx + 44);

            auto s_xxyy_xxxx = primBuffer.data(soff + 225 * idx + 45);

            auto s_xxyy_xxxy = primBuffer.data(soff + 225 * idx + 46);

            auto s_xxyy_xxxz = primBuffer.data(soff + 225 * idx + 47);

            auto s_xxyy_xxyy = primBuffer.data(soff + 225 * idx + 48);

            auto s_xxyy_xxyz = primBuffer.data(soff + 225 * idx + 49);

            auto s_xxyy_xxzz = primBuffer.data(soff + 225 * idx + 50);

            auto s_xxyy_xyyy = primBuffer.data(soff + 225 * idx + 51);

            auto s_xxyy_xyyz = primBuffer.data(soff + 225 * idx + 52);

            auto s_xxyy_xyzz = primBuffer.data(soff + 225 * idx + 53);

            auto s_xxyy_xzzz = primBuffer.data(soff + 225 * idx + 54);

            auto s_xxyy_yyyy = primBuffer.data(soff + 225 * idx + 55);

            auto s_xxyy_yyyz = primBuffer.data(soff + 225 * idx + 56);

            auto s_xxyy_yyzz = primBuffer.data(soff + 225 * idx + 57);

            auto s_xxyy_yzzz = primBuffer.data(soff + 225 * idx + 58);

            auto s_xxyy_zzzz = primBuffer.data(soff + 225 * idx + 59);

            auto s_xxyz_xxxx = primBuffer.data(soff + 225 * idx + 60);

            auto s_xxyz_xxxy = primBuffer.data(soff + 225 * idx + 61);

            auto s_xxyz_xxxz = primBuffer.data(soff + 225 * idx + 62);

            auto s_xxyz_xxyy = primBuffer.data(soff + 225 * idx + 63);

            auto s_xxyz_xxyz = primBuffer.data(soff + 225 * idx + 64);

            auto s_xxyz_xxzz = primBuffer.data(soff + 225 * idx + 65);

            auto s_xxyz_xyyy = primBuffer.data(soff + 225 * idx + 66);

            auto s_xxyz_xyyz = primBuffer.data(soff + 225 * idx + 67);

            auto s_xxyz_xyzz = primBuffer.data(soff + 225 * idx + 68);

            auto s_xxyz_xzzz = primBuffer.data(soff + 225 * idx + 69);

            auto s_xxyz_yyyy = primBuffer.data(soff + 225 * idx + 70);

            auto s_xxyz_yyyz = primBuffer.data(soff + 225 * idx + 71);

            auto s_xxyz_yyzz = primBuffer.data(soff + 225 * idx + 72);

            auto s_xxyz_yzzz = primBuffer.data(soff + 225 * idx + 73);

            auto s_xxyz_zzzz = primBuffer.data(soff + 225 * idx + 74);

            auto s_xxzz_xxxx = primBuffer.data(soff + 225 * idx + 75);

            auto s_xxzz_xxxy = primBuffer.data(soff + 225 * idx + 76);

            auto s_xxzz_xxxz = primBuffer.data(soff + 225 * idx + 77);

            auto s_xxzz_xxyy = primBuffer.data(soff + 225 * idx + 78);

            auto s_xxzz_xxyz = primBuffer.data(soff + 225 * idx + 79);

            auto s_xxzz_xxzz = primBuffer.data(soff + 225 * idx + 80);

            auto s_xxzz_xyyy = primBuffer.data(soff + 225 * idx + 81);

            auto s_xxzz_xyyz = primBuffer.data(soff + 225 * idx + 82);

            auto s_xxzz_xyzz = primBuffer.data(soff + 225 * idx + 83);

            auto s_xxzz_xzzz = primBuffer.data(soff + 225 * idx + 84);

            auto s_xxzz_yyyy = primBuffer.data(soff + 225 * idx + 85);

            auto s_xxzz_yyyz = primBuffer.data(soff + 225 * idx + 86);

            auto s_xxzz_yyzz = primBuffer.data(soff + 225 * idx + 87);

            auto s_xxzz_yzzz = primBuffer.data(soff + 225 * idx + 88);

            auto s_xxzz_zzzz = primBuffer.data(soff + 225 * idx + 89);

            auto s_xyyy_xxxx = primBuffer.data(soff + 225 * idx + 90);

            auto s_xyyy_xxxy = primBuffer.data(soff + 225 * idx + 91);

            auto s_xyyy_xxxz = primBuffer.data(soff + 225 * idx + 92);

            auto s_xyyy_xxyy = primBuffer.data(soff + 225 * idx + 93);

            auto s_xyyy_xxyz = primBuffer.data(soff + 225 * idx + 94);

            auto s_xyyy_xxzz = primBuffer.data(soff + 225 * idx + 95);

            auto s_xyyy_xyyy = primBuffer.data(soff + 225 * idx + 96);

            auto s_xyyy_xyyz = primBuffer.data(soff + 225 * idx + 97);

            auto s_xyyy_xyzz = primBuffer.data(soff + 225 * idx + 98);

            auto s_xyyy_xzzz = primBuffer.data(soff + 225 * idx + 99);

            auto s_xyyy_yyyy = primBuffer.data(soff + 225 * idx + 100);

            auto s_xyyy_yyyz = primBuffer.data(soff + 225 * idx + 101);

            auto s_xyyy_yyzz = primBuffer.data(soff + 225 * idx + 102);

            auto s_xyyy_yzzz = primBuffer.data(soff + 225 * idx + 103);

            auto s_xyyy_zzzz = primBuffer.data(soff + 225 * idx + 104);

            auto s_xyyz_xxxx = primBuffer.data(soff + 225 * idx + 105);

            auto s_xyyz_xxxy = primBuffer.data(soff + 225 * idx + 106);

            auto s_xyyz_xxxz = primBuffer.data(soff + 225 * idx + 107);

            auto s_xyyz_xxyy = primBuffer.data(soff + 225 * idx + 108);

            auto s_xyyz_xxyz = primBuffer.data(soff + 225 * idx + 109);

            auto s_xyyz_xxzz = primBuffer.data(soff + 225 * idx + 110);

            auto s_xyyz_xyyy = primBuffer.data(soff + 225 * idx + 111);

            auto s_xyyz_xyyz = primBuffer.data(soff + 225 * idx + 112);

            auto s_xyyz_xyzz = primBuffer.data(soff + 225 * idx + 113);

            auto s_xyyz_xzzz = primBuffer.data(soff + 225 * idx + 114);

            auto s_xyyz_yyyy = primBuffer.data(soff + 225 * idx + 115);

            auto s_xyyz_yyyz = primBuffer.data(soff + 225 * idx + 116);

            auto s_xyyz_yyzz = primBuffer.data(soff + 225 * idx + 117);

            auto s_xyyz_yzzz = primBuffer.data(soff + 225 * idx + 118);

            auto s_xyyz_zzzz = primBuffer.data(soff + 225 * idx + 119);

            auto s_xyzz_xxxx = primBuffer.data(soff + 225 * idx + 120);

            auto s_xyzz_xxxy = primBuffer.data(soff + 225 * idx + 121);

            auto s_xyzz_xxxz = primBuffer.data(soff + 225 * idx + 122);

            auto s_xyzz_xxyy = primBuffer.data(soff + 225 * idx + 123);

            auto s_xyzz_xxyz = primBuffer.data(soff + 225 * idx + 124);

            auto s_xyzz_xxzz = primBuffer.data(soff + 225 * idx + 125);

            auto s_xyzz_xyyy = primBuffer.data(soff + 225 * idx + 126);

            auto s_xyzz_xyyz = primBuffer.data(soff + 225 * idx + 127);

            auto s_xyzz_xyzz = primBuffer.data(soff + 225 * idx + 128);

            auto s_xyzz_xzzz = primBuffer.data(soff + 225 * idx + 129);

            auto s_xyzz_yyyy = primBuffer.data(soff + 225 * idx + 130);

            auto s_xyzz_yyyz = primBuffer.data(soff + 225 * idx + 131);

            auto s_xyzz_yyzz = primBuffer.data(soff + 225 * idx + 132);

            auto s_xyzz_yzzz = primBuffer.data(soff + 225 * idx + 133);

            auto s_xyzz_zzzz = primBuffer.data(soff + 225 * idx + 134);

            auto s_xzzz_xxxx = primBuffer.data(soff + 225 * idx + 135);

            auto s_xzzz_xxxy = primBuffer.data(soff + 225 * idx + 136);

            auto s_xzzz_xxxz = primBuffer.data(soff + 225 * idx + 137);

            auto s_xzzz_xxyy = primBuffer.data(soff + 225 * idx + 138);

            auto s_xzzz_xxyz = primBuffer.data(soff + 225 * idx + 139);

            auto s_xzzz_xxzz = primBuffer.data(soff + 225 * idx + 140);

            auto s_xzzz_xyyy = primBuffer.data(soff + 225 * idx + 141);

            auto s_xzzz_xyyz = primBuffer.data(soff + 225 * idx + 142);

            auto s_xzzz_xyzz = primBuffer.data(soff + 225 * idx + 143);

            auto s_xzzz_xzzz = primBuffer.data(soff + 225 * idx + 144);

            auto s_xzzz_yyyy = primBuffer.data(soff + 225 * idx + 145);

            auto s_xzzz_yyyz = primBuffer.data(soff + 225 * idx + 146);

            auto s_xzzz_yyzz = primBuffer.data(soff + 225 * idx + 147);

            auto s_xzzz_yzzz = primBuffer.data(soff + 225 * idx + 148);

            auto s_xzzz_zzzz = primBuffer.data(soff + 225 * idx + 149);

            auto s_yyyy_xxxx = primBuffer.data(soff + 225 * idx + 150);

            auto s_yyyy_xxxy = primBuffer.data(soff + 225 * idx + 151);

            auto s_yyyy_xxxz = primBuffer.data(soff + 225 * idx + 152);

            auto s_yyyy_xxyy = primBuffer.data(soff + 225 * idx + 153);

            auto s_yyyy_xxyz = primBuffer.data(soff + 225 * idx + 154);

            auto s_yyyy_xxzz = primBuffer.data(soff + 225 * idx + 155);

            auto s_yyyy_xyyy = primBuffer.data(soff + 225 * idx + 156);

            auto s_yyyy_xyyz = primBuffer.data(soff + 225 * idx + 157);

            auto s_yyyy_xyzz = primBuffer.data(soff + 225 * idx + 158);

            auto s_yyyy_xzzz = primBuffer.data(soff + 225 * idx + 159);

            auto s_yyyy_yyyy = primBuffer.data(soff + 225 * idx + 160);

            auto s_yyyy_yyyz = primBuffer.data(soff + 225 * idx + 161);

            auto s_yyyy_yyzz = primBuffer.data(soff + 225 * idx + 162);

            auto s_yyyy_yzzz = primBuffer.data(soff + 225 * idx + 163);

            auto s_yyyy_zzzz = primBuffer.data(soff + 225 * idx + 164);

            auto s_yyyz_xxxx = primBuffer.data(soff + 225 * idx + 165);

            auto s_yyyz_xxxy = primBuffer.data(soff + 225 * idx + 166);

            auto s_yyyz_xxxz = primBuffer.data(soff + 225 * idx + 167);

            auto s_yyyz_xxyy = primBuffer.data(soff + 225 * idx + 168);

            auto s_yyyz_xxyz = primBuffer.data(soff + 225 * idx + 169);

            auto s_yyyz_xxzz = primBuffer.data(soff + 225 * idx + 170);

            auto s_yyyz_xyyy = primBuffer.data(soff + 225 * idx + 171);

            auto s_yyyz_xyyz = primBuffer.data(soff + 225 * idx + 172);

            auto s_yyyz_xyzz = primBuffer.data(soff + 225 * idx + 173);

            auto s_yyyz_xzzz = primBuffer.data(soff + 225 * idx + 174);

            auto s_yyyz_yyyy = primBuffer.data(soff + 225 * idx + 175);

            auto s_yyyz_yyyz = primBuffer.data(soff + 225 * idx + 176);

            auto s_yyyz_yyzz = primBuffer.data(soff + 225 * idx + 177);

            auto s_yyyz_yzzz = primBuffer.data(soff + 225 * idx + 178);

            auto s_yyyz_zzzz = primBuffer.data(soff + 225 * idx + 179);

            auto s_yyzz_xxxx = primBuffer.data(soff + 225 * idx + 180);

            auto s_yyzz_xxxy = primBuffer.data(soff + 225 * idx + 181);

            auto s_yyzz_xxxz = primBuffer.data(soff + 225 * idx + 182);

            auto s_yyzz_xxyy = primBuffer.data(soff + 225 * idx + 183);

            auto s_yyzz_xxyz = primBuffer.data(soff + 225 * idx + 184);

            auto s_yyzz_xxzz = primBuffer.data(soff + 225 * idx + 185);

            auto s_yyzz_xyyy = primBuffer.data(soff + 225 * idx + 186);

            auto s_yyzz_xyyz = primBuffer.data(soff + 225 * idx + 187);

            auto s_yyzz_xyzz = primBuffer.data(soff + 225 * idx + 188);

            auto s_yyzz_xzzz = primBuffer.data(soff + 225 * idx + 189);

            auto s_yyzz_yyyy = primBuffer.data(soff + 225 * idx + 190);

            auto s_yyzz_yyyz = primBuffer.data(soff + 225 * idx + 191);

            auto s_yyzz_yyzz = primBuffer.data(soff + 225 * idx + 192);

            auto s_yyzz_yzzz = primBuffer.data(soff + 225 * idx + 193);

            auto s_yyzz_zzzz = primBuffer.data(soff + 225 * idx + 194);

            auto s_yzzz_xxxx = primBuffer.data(soff + 225 * idx + 195);

            auto s_yzzz_xxxy = primBuffer.data(soff + 225 * idx + 196);

            auto s_yzzz_xxxz = primBuffer.data(soff + 225 * idx + 197);

            auto s_yzzz_xxyy = primBuffer.data(soff + 225 * idx + 198);

            auto s_yzzz_xxyz = primBuffer.data(soff + 225 * idx + 199);

            auto s_yzzz_xxzz = primBuffer.data(soff + 225 * idx + 200);

            auto s_yzzz_xyyy = primBuffer.data(soff + 225 * idx + 201);

            auto s_yzzz_xyyz = primBuffer.data(soff + 225 * idx + 202);

            auto s_yzzz_xyzz = primBuffer.data(soff + 225 * idx + 203);

            auto s_yzzz_xzzz = primBuffer.data(soff + 225 * idx + 204);

            auto s_yzzz_yyyy = primBuffer.data(soff + 225 * idx + 205);

            auto s_yzzz_yyyz = primBuffer.data(soff + 225 * idx + 206);

            auto s_yzzz_yyzz = primBuffer.data(soff + 225 * idx + 207);

            auto s_yzzz_yzzz = primBuffer.data(soff + 225 * idx + 208);

            auto s_yzzz_zzzz = primBuffer.data(soff + 225 * idx + 209);

            auto s_zzzz_xxxx = primBuffer.data(soff + 225 * idx + 210);

            auto s_zzzz_xxxy = primBuffer.data(soff + 225 * idx + 211);

            auto s_zzzz_xxxz = primBuffer.data(soff + 225 * idx + 212);

            auto s_zzzz_xxyy = primBuffer.data(soff + 225 * idx + 213);

            auto s_zzzz_xxyz = primBuffer.data(soff + 225 * idx + 214);

            auto s_zzzz_xxzz = primBuffer.data(soff + 225 * idx + 215);

            auto s_zzzz_xyyy = primBuffer.data(soff + 225 * idx + 216);

            auto s_zzzz_xyyz = primBuffer.data(soff + 225 * idx + 217);

            auto s_zzzz_xyzz = primBuffer.data(soff + 225 * idx + 218);

            auto s_zzzz_xzzz = primBuffer.data(soff + 225 * idx + 219);

            auto s_zzzz_yyyy = primBuffer.data(soff + 225 * idx + 220);

            auto s_zzzz_yyyz = primBuffer.data(soff + 225 * idx + 221);

            auto s_zzzz_yyzz = primBuffer.data(soff + 225 * idx + 222);

            auto s_zzzz_yzzz = primBuffer.data(soff + 225 * idx + 223);

            auto s_zzzz_zzzz = primBuffer.data(soff + 225 * idx + 224);

            // set up pointers to (G|T|G) integrals

            auto t_xxxx_xxxx = primBuffer.data(koff + 225 * idx);

            auto t_xxxx_xxxy = primBuffer.data(koff + 225 * idx + 1);

            auto t_xxxx_xxxz = primBuffer.data(koff + 225 * idx + 2);

            auto t_xxxx_xxyy = primBuffer.data(koff + 225 * idx + 3);

            auto t_xxxx_xxyz = primBuffer.data(koff + 225 * idx + 4);

            auto t_xxxx_xxzz = primBuffer.data(koff + 225 * idx + 5);

            auto t_xxxx_xyyy = primBuffer.data(koff + 225 * idx + 6);

            auto t_xxxx_xyyz = primBuffer.data(koff + 225 * idx + 7);

            auto t_xxxx_xyzz = primBuffer.data(koff + 225 * idx + 8);

            auto t_xxxx_xzzz = primBuffer.data(koff + 225 * idx + 9);

            auto t_xxxx_yyyy = primBuffer.data(koff + 225 * idx + 10);

            auto t_xxxx_yyyz = primBuffer.data(koff + 225 * idx + 11);

            auto t_xxxx_yyzz = primBuffer.data(koff + 225 * idx + 12);

            auto t_xxxx_yzzz = primBuffer.data(koff + 225 * idx + 13);

            auto t_xxxx_zzzz = primBuffer.data(koff + 225 * idx + 14);

            auto t_xxxy_xxxx = primBuffer.data(koff + 225 * idx + 15);

            auto t_xxxy_xxxy = primBuffer.data(koff + 225 * idx + 16);

            auto t_xxxy_xxxz = primBuffer.data(koff + 225 * idx + 17);

            auto t_xxxy_xxyy = primBuffer.data(koff + 225 * idx + 18);

            auto t_xxxy_xxyz = primBuffer.data(koff + 225 * idx + 19);

            auto t_xxxy_xxzz = primBuffer.data(koff + 225 * idx + 20);

            auto t_xxxy_xyyy = primBuffer.data(koff + 225 * idx + 21);

            auto t_xxxy_xyyz = primBuffer.data(koff + 225 * idx + 22);

            auto t_xxxy_xyzz = primBuffer.data(koff + 225 * idx + 23);

            auto t_xxxy_xzzz = primBuffer.data(koff + 225 * idx + 24);

            auto t_xxxy_yyyy = primBuffer.data(koff + 225 * idx + 25);

            auto t_xxxy_yyyz = primBuffer.data(koff + 225 * idx + 26);

            auto t_xxxy_yyzz = primBuffer.data(koff + 225 * idx + 27);

            auto t_xxxy_yzzz = primBuffer.data(koff + 225 * idx + 28);

            auto t_xxxy_zzzz = primBuffer.data(koff + 225 * idx + 29);

            auto t_xxxz_xxxx = primBuffer.data(koff + 225 * idx + 30);

            auto t_xxxz_xxxy = primBuffer.data(koff + 225 * idx + 31);

            auto t_xxxz_xxxz = primBuffer.data(koff + 225 * idx + 32);

            auto t_xxxz_xxyy = primBuffer.data(koff + 225 * idx + 33);

            auto t_xxxz_xxyz = primBuffer.data(koff + 225 * idx + 34);

            auto t_xxxz_xxzz = primBuffer.data(koff + 225 * idx + 35);

            auto t_xxxz_xyyy = primBuffer.data(koff + 225 * idx + 36);

            auto t_xxxz_xyyz = primBuffer.data(koff + 225 * idx + 37);

            auto t_xxxz_xyzz = primBuffer.data(koff + 225 * idx + 38);

            auto t_xxxz_xzzz = primBuffer.data(koff + 225 * idx + 39);

            auto t_xxxz_yyyy = primBuffer.data(koff + 225 * idx + 40);

            auto t_xxxz_yyyz = primBuffer.data(koff + 225 * idx + 41);

            auto t_xxxz_yyzz = primBuffer.data(koff + 225 * idx + 42);

            auto t_xxxz_yzzz = primBuffer.data(koff + 225 * idx + 43);

            auto t_xxxz_zzzz = primBuffer.data(koff + 225 * idx + 44);

            auto t_xxyy_xxxx = primBuffer.data(koff + 225 * idx + 45);

            auto t_xxyy_xxxy = primBuffer.data(koff + 225 * idx + 46);

            auto t_xxyy_xxxz = primBuffer.data(koff + 225 * idx + 47);

            auto t_xxyy_xxyy = primBuffer.data(koff + 225 * idx + 48);

            auto t_xxyy_xxyz = primBuffer.data(koff + 225 * idx + 49);

            auto t_xxyy_xxzz = primBuffer.data(koff + 225 * idx + 50);

            auto t_xxyy_xyyy = primBuffer.data(koff + 225 * idx + 51);

            auto t_xxyy_xyyz = primBuffer.data(koff + 225 * idx + 52);

            auto t_xxyy_xyzz = primBuffer.data(koff + 225 * idx + 53);

            auto t_xxyy_xzzz = primBuffer.data(koff + 225 * idx + 54);

            auto t_xxyy_yyyy = primBuffer.data(koff + 225 * idx + 55);

            auto t_xxyy_yyyz = primBuffer.data(koff + 225 * idx + 56);

            auto t_xxyy_yyzz = primBuffer.data(koff + 225 * idx + 57);

            auto t_xxyy_yzzz = primBuffer.data(koff + 225 * idx + 58);

            auto t_xxyy_zzzz = primBuffer.data(koff + 225 * idx + 59);

            auto t_xxyz_xxxx = primBuffer.data(koff + 225 * idx + 60);

            auto t_xxyz_xxxy = primBuffer.data(koff + 225 * idx + 61);

            auto t_xxyz_xxxz = primBuffer.data(koff + 225 * idx + 62);

            auto t_xxyz_xxyy = primBuffer.data(koff + 225 * idx + 63);

            auto t_xxyz_xxyz = primBuffer.data(koff + 225 * idx + 64);

            auto t_xxyz_xxzz = primBuffer.data(koff + 225 * idx + 65);

            auto t_xxyz_xyyy = primBuffer.data(koff + 225 * idx + 66);

            auto t_xxyz_xyyz = primBuffer.data(koff + 225 * idx + 67);

            auto t_xxyz_xyzz = primBuffer.data(koff + 225 * idx + 68);

            auto t_xxyz_xzzz = primBuffer.data(koff + 225 * idx + 69);

            auto t_xxyz_yyyy = primBuffer.data(koff + 225 * idx + 70);

            auto t_xxyz_yyyz = primBuffer.data(koff + 225 * idx + 71);

            auto t_xxyz_yyzz = primBuffer.data(koff + 225 * idx + 72);

            auto t_xxyz_yzzz = primBuffer.data(koff + 225 * idx + 73);

            auto t_xxyz_zzzz = primBuffer.data(koff + 225 * idx + 74);

            auto t_xxzz_xxxx = primBuffer.data(koff + 225 * idx + 75);

            auto t_xxzz_xxxy = primBuffer.data(koff + 225 * idx + 76);

            auto t_xxzz_xxxz = primBuffer.data(koff + 225 * idx + 77);

            auto t_xxzz_xxyy = primBuffer.data(koff + 225 * idx + 78);

            auto t_xxzz_xxyz = primBuffer.data(koff + 225 * idx + 79);

            auto t_xxzz_xxzz = primBuffer.data(koff + 225 * idx + 80);

            auto t_xxzz_xyyy = primBuffer.data(koff + 225 * idx + 81);

            auto t_xxzz_xyyz = primBuffer.data(koff + 225 * idx + 82);

            auto t_xxzz_xyzz = primBuffer.data(koff + 225 * idx + 83);

            auto t_xxzz_xzzz = primBuffer.data(koff + 225 * idx + 84);

            auto t_xxzz_yyyy = primBuffer.data(koff + 225 * idx + 85);

            auto t_xxzz_yyyz = primBuffer.data(koff + 225 * idx + 86);

            auto t_xxzz_yyzz = primBuffer.data(koff + 225 * idx + 87);

            auto t_xxzz_yzzz = primBuffer.data(koff + 225 * idx + 88);

            auto t_xxzz_zzzz = primBuffer.data(koff + 225 * idx + 89);

            auto t_xyyy_xxxx = primBuffer.data(koff + 225 * idx + 90);

            auto t_xyyy_xxxy = primBuffer.data(koff + 225 * idx + 91);

            auto t_xyyy_xxxz = primBuffer.data(koff + 225 * idx + 92);

            auto t_xyyy_xxyy = primBuffer.data(koff + 225 * idx + 93);

            auto t_xyyy_xxyz = primBuffer.data(koff + 225 * idx + 94);

            auto t_xyyy_xxzz = primBuffer.data(koff + 225 * idx + 95);

            auto t_xyyy_xyyy = primBuffer.data(koff + 225 * idx + 96);

            auto t_xyyy_xyyz = primBuffer.data(koff + 225 * idx + 97);

            auto t_xyyy_xyzz = primBuffer.data(koff + 225 * idx + 98);

            auto t_xyyy_xzzz = primBuffer.data(koff + 225 * idx + 99);

            auto t_xyyy_yyyy = primBuffer.data(koff + 225 * idx + 100);

            auto t_xyyy_yyyz = primBuffer.data(koff + 225 * idx + 101);

            auto t_xyyy_yyzz = primBuffer.data(koff + 225 * idx + 102);

            auto t_xyyy_yzzz = primBuffer.data(koff + 225 * idx + 103);

            auto t_xyyy_zzzz = primBuffer.data(koff + 225 * idx + 104);

            auto t_xyyz_xxxx = primBuffer.data(koff + 225 * idx + 105);

            auto t_xyyz_xxxy = primBuffer.data(koff + 225 * idx + 106);

            auto t_xyyz_xxxz = primBuffer.data(koff + 225 * idx + 107);

            auto t_xyyz_xxyy = primBuffer.data(koff + 225 * idx + 108);

            auto t_xyyz_xxyz = primBuffer.data(koff + 225 * idx + 109);

            auto t_xyyz_xxzz = primBuffer.data(koff + 225 * idx + 110);

            auto t_xyyz_xyyy = primBuffer.data(koff + 225 * idx + 111);

            auto t_xyyz_xyyz = primBuffer.data(koff + 225 * idx + 112);

            auto t_xyyz_xyzz = primBuffer.data(koff + 225 * idx + 113);

            auto t_xyyz_xzzz = primBuffer.data(koff + 225 * idx + 114);

            auto t_xyyz_yyyy = primBuffer.data(koff + 225 * idx + 115);

            auto t_xyyz_yyyz = primBuffer.data(koff + 225 * idx + 116);

            auto t_xyyz_yyzz = primBuffer.data(koff + 225 * idx + 117);

            auto t_xyyz_yzzz = primBuffer.data(koff + 225 * idx + 118);

            auto t_xyyz_zzzz = primBuffer.data(koff + 225 * idx + 119);

            auto t_xyzz_xxxx = primBuffer.data(koff + 225 * idx + 120);

            auto t_xyzz_xxxy = primBuffer.data(koff + 225 * idx + 121);

            auto t_xyzz_xxxz = primBuffer.data(koff + 225 * idx + 122);

            auto t_xyzz_xxyy = primBuffer.data(koff + 225 * idx + 123);

            auto t_xyzz_xxyz = primBuffer.data(koff + 225 * idx + 124);

            auto t_xyzz_xxzz = primBuffer.data(koff + 225 * idx + 125);

            auto t_xyzz_xyyy = primBuffer.data(koff + 225 * idx + 126);

            auto t_xyzz_xyyz = primBuffer.data(koff + 225 * idx + 127);

            auto t_xyzz_xyzz = primBuffer.data(koff + 225 * idx + 128);

            auto t_xyzz_xzzz = primBuffer.data(koff + 225 * idx + 129);

            auto t_xyzz_yyyy = primBuffer.data(koff + 225 * idx + 130);

            auto t_xyzz_yyyz = primBuffer.data(koff + 225 * idx + 131);

            auto t_xyzz_yyzz = primBuffer.data(koff + 225 * idx + 132);

            auto t_xyzz_yzzz = primBuffer.data(koff + 225 * idx + 133);

            auto t_xyzz_zzzz = primBuffer.data(koff + 225 * idx + 134);

            auto t_xzzz_xxxx = primBuffer.data(koff + 225 * idx + 135);

            auto t_xzzz_xxxy = primBuffer.data(koff + 225 * idx + 136);

            auto t_xzzz_xxxz = primBuffer.data(koff + 225 * idx + 137);

            auto t_xzzz_xxyy = primBuffer.data(koff + 225 * idx + 138);

            auto t_xzzz_xxyz = primBuffer.data(koff + 225 * idx + 139);

            auto t_xzzz_xxzz = primBuffer.data(koff + 225 * idx + 140);

            auto t_xzzz_xyyy = primBuffer.data(koff + 225 * idx + 141);

            auto t_xzzz_xyyz = primBuffer.data(koff + 225 * idx + 142);

            auto t_xzzz_xyzz = primBuffer.data(koff + 225 * idx + 143);

            auto t_xzzz_xzzz = primBuffer.data(koff + 225 * idx + 144);

            auto t_xzzz_yyyy = primBuffer.data(koff + 225 * idx + 145);

            auto t_xzzz_yyyz = primBuffer.data(koff + 225 * idx + 146);

            auto t_xzzz_yyzz = primBuffer.data(koff + 225 * idx + 147);

            auto t_xzzz_yzzz = primBuffer.data(koff + 225 * idx + 148);

            auto t_xzzz_zzzz = primBuffer.data(koff + 225 * idx + 149);

            auto t_yyyy_xxxx = primBuffer.data(koff + 225 * idx + 150);

            auto t_yyyy_xxxy = primBuffer.data(koff + 225 * idx + 151);

            auto t_yyyy_xxxz = primBuffer.data(koff + 225 * idx + 152);

            auto t_yyyy_xxyy = primBuffer.data(koff + 225 * idx + 153);

            auto t_yyyy_xxyz = primBuffer.data(koff + 225 * idx + 154);

            auto t_yyyy_xxzz = primBuffer.data(koff + 225 * idx + 155);

            auto t_yyyy_xyyy = primBuffer.data(koff + 225 * idx + 156);

            auto t_yyyy_xyyz = primBuffer.data(koff + 225 * idx + 157);

            auto t_yyyy_xyzz = primBuffer.data(koff + 225 * idx + 158);

            auto t_yyyy_xzzz = primBuffer.data(koff + 225 * idx + 159);

            auto t_yyyy_yyyy = primBuffer.data(koff + 225 * idx + 160);

            auto t_yyyy_yyyz = primBuffer.data(koff + 225 * idx + 161);

            auto t_yyyy_yyzz = primBuffer.data(koff + 225 * idx + 162);

            auto t_yyyy_yzzz = primBuffer.data(koff + 225 * idx + 163);

            auto t_yyyy_zzzz = primBuffer.data(koff + 225 * idx + 164);

            auto t_yyyz_xxxx = primBuffer.data(koff + 225 * idx + 165);

            auto t_yyyz_xxxy = primBuffer.data(koff + 225 * idx + 166);

            auto t_yyyz_xxxz = primBuffer.data(koff + 225 * idx + 167);

            auto t_yyyz_xxyy = primBuffer.data(koff + 225 * idx + 168);

            auto t_yyyz_xxyz = primBuffer.data(koff + 225 * idx + 169);

            auto t_yyyz_xxzz = primBuffer.data(koff + 225 * idx + 170);

            auto t_yyyz_xyyy = primBuffer.data(koff + 225 * idx + 171);

            auto t_yyyz_xyyz = primBuffer.data(koff + 225 * idx + 172);

            auto t_yyyz_xyzz = primBuffer.data(koff + 225 * idx + 173);

            auto t_yyyz_xzzz = primBuffer.data(koff + 225 * idx + 174);

            auto t_yyyz_yyyy = primBuffer.data(koff + 225 * idx + 175);

            auto t_yyyz_yyyz = primBuffer.data(koff + 225 * idx + 176);

            auto t_yyyz_yyzz = primBuffer.data(koff + 225 * idx + 177);

            auto t_yyyz_yzzz = primBuffer.data(koff + 225 * idx + 178);

            auto t_yyyz_zzzz = primBuffer.data(koff + 225 * idx + 179);

            auto t_yyzz_xxxx = primBuffer.data(koff + 225 * idx + 180);

            auto t_yyzz_xxxy = primBuffer.data(koff + 225 * idx + 181);

            auto t_yyzz_xxxz = primBuffer.data(koff + 225 * idx + 182);

            auto t_yyzz_xxyy = primBuffer.data(koff + 225 * idx + 183);

            auto t_yyzz_xxyz = primBuffer.data(koff + 225 * idx + 184);

            auto t_yyzz_xxzz = primBuffer.data(koff + 225 * idx + 185);

            auto t_yyzz_xyyy = primBuffer.data(koff + 225 * idx + 186);

            auto t_yyzz_xyyz = primBuffer.data(koff + 225 * idx + 187);

            auto t_yyzz_xyzz = primBuffer.data(koff + 225 * idx + 188);

            auto t_yyzz_xzzz = primBuffer.data(koff + 225 * idx + 189);

            auto t_yyzz_yyyy = primBuffer.data(koff + 225 * idx + 190);

            auto t_yyzz_yyyz = primBuffer.data(koff + 225 * idx + 191);

            auto t_yyzz_yyzz = primBuffer.data(koff + 225 * idx + 192);

            auto t_yyzz_yzzz = primBuffer.data(koff + 225 * idx + 193);

            auto t_yyzz_zzzz = primBuffer.data(koff + 225 * idx + 194);

            auto t_yzzz_xxxx = primBuffer.data(koff + 225 * idx + 195);

            auto t_yzzz_xxxy = primBuffer.data(koff + 225 * idx + 196);

            auto t_yzzz_xxxz = primBuffer.data(koff + 225 * idx + 197);

            auto t_yzzz_xxyy = primBuffer.data(koff + 225 * idx + 198);

            auto t_yzzz_xxyz = primBuffer.data(koff + 225 * idx + 199);

            auto t_yzzz_xxzz = primBuffer.data(koff + 225 * idx + 200);

            auto t_yzzz_xyyy = primBuffer.data(koff + 225 * idx + 201);

            auto t_yzzz_xyyz = primBuffer.data(koff + 225 * idx + 202);

            auto t_yzzz_xyzz = primBuffer.data(koff + 225 * idx + 203);

            auto t_yzzz_xzzz = primBuffer.data(koff + 225 * idx + 204);

            auto t_yzzz_yyyy = primBuffer.data(koff + 225 * idx + 205);

            auto t_yzzz_yyyz = primBuffer.data(koff + 225 * idx + 206);

            auto t_yzzz_yyzz = primBuffer.data(koff + 225 * idx + 207);

            auto t_yzzz_yzzz = primBuffer.data(koff + 225 * idx + 208);

            auto t_yzzz_zzzz = primBuffer.data(koff + 225 * idx + 209);

            auto t_zzzz_xxxx = primBuffer.data(koff + 225 * idx + 210);

            auto t_zzzz_xxxy = primBuffer.data(koff + 225 * idx + 211);

            auto t_zzzz_xxxz = primBuffer.data(koff + 225 * idx + 212);

            auto t_zzzz_xxyy = primBuffer.data(koff + 225 * idx + 213);

            auto t_zzzz_xxyz = primBuffer.data(koff + 225 * idx + 214);

            auto t_zzzz_xxzz = primBuffer.data(koff + 225 * idx + 215);

            auto t_zzzz_xyyy = primBuffer.data(koff + 225 * idx + 216);

            auto t_zzzz_xyyz = primBuffer.data(koff + 225 * idx + 217);

            auto t_zzzz_xyzz = primBuffer.data(koff + 225 * idx + 218);

            auto t_zzzz_xzzz = primBuffer.data(koff + 225 * idx + 219);

            auto t_zzzz_yyyy = primBuffer.data(koff + 225 * idx + 220);

            auto t_zzzz_yyyz = primBuffer.data(koff + 225 * idx + 221);

            auto t_zzzz_yyzz = primBuffer.data(koff + 225 * idx + 222);

            auto t_zzzz_yzzz = primBuffer.data(koff + 225 * idx + 223);

            auto t_zzzz_zzzz = primBuffer.data(koff + 225 * idx + 224);

            #pragma omp simd aligned(fx, fz, fga, pax, pay, paz, s_xxx_xxx, s_xxx_xxy,\
                                     s_xxx_xxz, s_xxx_xyy, s_xxx_xyz, s_xxx_xzz,\
                                     s_xxx_yyy, s_xxx_yyz, s_xxx_yzz, s_xxx_zzz,\
                                     s_xxy_xxx, s_xxy_xxy, s_xxy_xxz, s_xxy_xyy,\
                                     s_xxy_xyz, s_xxy_xzz, s_xxy_yyy, s_xxy_yyz,\
                                     s_xxy_yzz, s_xxy_zzz, s_xxz_xxx, s_xxz_xxy,\
                                     s_xxz_xxz, s_xxz_xyy, s_xxz_xyz, s_xxz_xzz,\
                                     s_xxz_yyy, s_xxz_yyz, s_xxz_yzz, s_xxz_zzz,\
                                     s_xyy_xxx, s_xyy_xxy, s_xyy_xxz, s_xyy_xyy,\
                                     s_xyy_xyz, s_xyy_xzz, s_xyy_yyy, s_xyy_yyz,\
                                     s_xyy_yzz, s_xyy_zzz, s_xyz_xxx, s_xyz_xxy,\
                                     s_xyz_xxz, s_xyz_xyy, s_xyz_xyz, s_xyz_xzz,\
                                     s_xyz_yyy, s_xyz_yyz, s_xyz_yzz, s_xyz_zzz,\
                                     s_xzz_xxx, s_xzz_xxy, s_xzz_xxz, s_xzz_xyy,\
                                     s_xzz_xyz, s_xzz_xzz, s_xzz_yyy, s_xzz_yyz,\
                                     s_xzz_yzz, s_xzz_zzz, s_yyy_xxx, s_yyy_xxy,\
                                     s_yyy_xxz, s_yyy_xyy, s_yyy_xyz, s_yyy_xzz,\
                                     s_yyy_yyy, s_yyy_yyz, s_yyy_yzz, s_yyy_zzz,\
                                     s_yyz_xxx, s_yyz_xxy, s_yyz_xxz, s_yyz_xyy,\
                                     s_yyz_xyz, s_yyz_xzz, s_yyz_yyy, s_yyz_yyz,\
                                     s_yyz_yzz, s_yyz_zzz, s_yzz_xxx, s_yzz_xxy,\
                                     s_yzz_xxz, s_yzz_xyy, s_yzz_xyz, s_yzz_xzz,\
                                     s_yzz_yyy, s_yzz_yyz, s_yzz_yzz, s_yzz_zzz,\
                                     s_zzz_xxx, s_zzz_xxy, s_zzz_xxz, s_zzz_xyy,\
                                     s_zzz_xyz, s_zzz_xzz, s_zzz_yyy, s_zzz_yyz,\
                                     s_zzz_yzz, s_zzz_zzz, t_xxx_xxx, t_xxx_xxy,\
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
                                     t_zzz_yzz, t_zzz_zzz, s_xx_xxxx, s_xx_xxxy,\
                                     s_xx_xxxz, s_xx_xxyy, s_xx_xxyz, s_xx_xxzz,\
                                     s_xx_xyyy, s_xx_xyyz, s_xx_xyzz, s_xx_xzzz,\
                                     s_xx_yyyy, s_xx_yyyz, s_xx_yyzz, s_xx_yzzz,\
                                     s_xx_zzzz, s_xy_xxxx, s_xy_xxxy, s_xy_xxxz,\
                                     s_xy_xxyy, s_xy_xxyz, s_xy_xxzz, s_xy_xyyy,\
                                     s_xy_xyyz, s_xy_xyzz, s_xy_xzzz, s_xy_yyyy,\
                                     s_xy_yyyz, s_xy_yyzz, s_xy_yzzz, s_xy_zzzz,\
                                     s_xz_xxxx, s_xz_xxxy, s_xz_xxxz, s_xz_xxyy,\
                                     s_xz_xxyz, s_xz_xxzz, s_xz_xyyy, s_xz_xyyz,\
                                     s_xz_xyzz, s_xz_xzzz, s_xz_yyyy, s_xz_yyyz,\
                                     s_xz_yyzz, s_xz_yzzz, s_xz_zzzz, s_yy_xxxx,\
                                     s_yy_xxxy, s_yy_xxxz, s_yy_xxyy, s_yy_xxyz,\
                                     s_yy_xxzz, s_yy_xyyy, s_yy_xyyz, s_yy_xyzz,\
                                     s_yy_xzzz, s_yy_yyyy, s_yy_yyyz, s_yy_yyzz,\
                                     s_yy_yzzz, s_yy_zzzz, s_yz_xxxx, s_yz_xxxy,\
                                     s_yz_xxxz, s_yz_xxyy, s_yz_xxyz, s_yz_xxzz,\
                                     s_yz_xyyy, s_yz_xyyz, s_yz_xyzz, s_yz_xzzz,\
                                     s_yz_yyyy, s_yz_yyyz, s_yz_yyzz, s_yz_yzzz,\
                                     s_yz_zzzz, s_zz_xxxx, s_zz_xxxy, s_zz_xxxz,\
                                     s_zz_xxyy, s_zz_xxyz, s_zz_xxzz, s_zz_xyyy,\
                                     s_zz_xyyz, s_zz_xyzz, s_zz_xzzz, s_zz_yyyy,\
                                     s_zz_yyyz, s_zz_yyzz, s_zz_yzzz, s_zz_zzzz,\
                                     t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy,\
                                     t_xx_xxyz, t_xx_xxzz, t_xx_xyyy, t_xx_xyyz,\
                                     t_xx_xyzz, t_xx_xzzz, t_xx_yyyy, t_xx_yyyz,\
                                     t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx,\
                                     t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz,\
                                     t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz,\
                                     t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz,\
                                     t_xy_yzzz, t_xy_zzzz, t_xz_xxxx, t_xz_xxxy,\
                                     t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, t_xz_xxzz,\
                                     t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz,\
                                     t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz,\
                                     t_xz_zzzz, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz,\
                                     t_yy_xxyy, t_yy_xxyz, t_yy_xxzz, t_yy_xyyy,\
                                     t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy,\
                                     t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, t_yy_zzzz,\
                                     t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, t_yz_xxyy,\
                                     t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz,\
                                     t_yz_xyzz, t_yz_xzzz, t_yz_yyyy, t_yz_yyyz,\
                                     t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_zz_xxxx,\
                                     t_zz_xxxy, t_zz_xxxz, t_zz_xxyy, t_zz_xxyz,\
                                     t_zz_xxzz, t_zz_xyyy, t_zz_xyyz, t_zz_xyzz,\
                                     t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz,\
                                     t_zz_yzzz, t_zz_zzzz, s_xxx_xxxx, s_xxx_xxxy,\
                                     s_xxx_xxxz, s_xxx_xxyy, s_xxx_xxyz, s_xxx_xxzz,\
                                     s_xxx_xyyy, s_xxx_xyyz, s_xxx_xyzz, s_xxx_xzzz,\
                                     s_xxx_yyyy, s_xxx_yyyz, s_xxx_yyzz, s_xxx_yzzz,\
                                     s_xxx_zzzz, s_xxy_xxxx, s_xxy_xxxy, s_xxy_xxxz,\
                                     s_xxy_xxyy, s_xxy_xxyz, s_xxy_xxzz, s_xxy_xyyy,\
                                     s_xxy_xyyz, s_xxy_xyzz, s_xxy_xzzz, s_xxy_yyyy,\
                                     s_xxy_yyyz, s_xxy_yyzz, s_xxy_yzzz, s_xxy_zzzz,\
                                     s_xxz_xxxx, s_xxz_xxxy, s_xxz_xxxz, s_xxz_xxyy,\
                                     s_xxz_xxyz, s_xxz_xxzz, s_xxz_xyyy, s_xxz_xyyz,\
                                     s_xxz_xyzz, s_xxz_xzzz, s_xxz_yyyy, s_xxz_yyyz,\
                                     s_xxz_yyzz, s_xxz_yzzz, s_xxz_zzzz, s_xyy_xxxx,\
                                     s_xyy_xxxy, s_xyy_xxxz, s_xyy_xxyy, s_xyy_xxyz,\
                                     s_xyy_xxzz, s_xyy_xyyy, s_xyy_xyyz, s_xyy_xyzz,\
                                     s_xyy_xzzz, s_xyy_yyyy, s_xyy_yyyz, s_xyy_yyzz,\
                                     s_xyy_yzzz, s_xyy_zzzz, s_xyz_xxxx, s_xyz_xxxy,\
                                     s_xyz_xxxz, s_xyz_xxyy, s_xyz_xxyz, s_xyz_xxzz,\
                                     s_xyz_xyyy, s_xyz_xyyz, s_xyz_xyzz, s_xyz_xzzz,\
                                     s_xyz_yyyy, s_xyz_yyyz, s_xyz_yyzz, s_xyz_yzzz,\
                                     s_xyz_zzzz, s_xzz_xxxx, s_xzz_xxxy, s_xzz_xxxz,\
                                     s_xzz_xxyy, s_xzz_xxyz, s_xzz_xxzz, s_xzz_xyyy,\
                                     s_xzz_xyyz, s_xzz_xyzz, s_xzz_xzzz, s_xzz_yyyy,\
                                     s_xzz_yyyz, s_xzz_yyzz, s_xzz_yzzz, s_xzz_zzzz,\
                                     s_yyy_xxxx, s_yyy_xxxy, s_yyy_xxxz, s_yyy_xxyy,\
                                     s_yyy_xxyz, s_yyy_xxzz, s_yyy_xyyy, s_yyy_xyyz,\
                                     s_yyy_xyzz, s_yyy_xzzz, s_yyy_yyyy, s_yyy_yyyz,\
                                     s_yyy_yyzz, s_yyy_yzzz, s_yyy_zzzz, s_yyz_xxxx,\
                                     s_yyz_xxxy, s_yyz_xxxz, s_yyz_xxyy, s_yyz_xxyz,\
                                     s_yyz_xxzz, s_yyz_xyyy, s_yyz_xyyz, s_yyz_xyzz,\
                                     s_yyz_xzzz, s_yyz_yyyy, s_yyz_yyyz, s_yyz_yyzz,\
                                     s_yyz_yzzz, s_yyz_zzzz, s_yzz_xxxx, s_yzz_xxxy,\
                                     s_yzz_xxxz, s_yzz_xxyy, s_yzz_xxyz, s_yzz_xxzz,\
                                     s_yzz_xyyy, s_yzz_xyyz, s_yzz_xyzz, s_yzz_xzzz,\
                                     s_yzz_yyyy, s_yzz_yyyz, s_yzz_yyzz, s_yzz_yzzz,\
                                     s_yzz_zzzz, s_zzz_xxxx, s_zzz_xxxy, s_zzz_xxxz,\
                                     s_zzz_xxyy, s_zzz_xxyz, s_zzz_xxzz, s_zzz_xyyy,\
                                     s_zzz_xyyz, s_zzz_xyzz, s_zzz_xzzz, s_zzz_yyyy,\
                                     s_zzz_yyyz, s_zzz_yyzz, s_zzz_yzzz, s_zzz_zzzz,\
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
                                     t_zzz_yzzz, t_zzz_zzzz, s_xxxx_xxxx, s_xxxx_xxxy,\
                                     s_xxxx_xxxz, s_xxxx_xxyy, s_xxxx_xxyz, s_xxxx_xxzz,\
                                     s_xxxx_xyyy, s_xxxx_xyyz, s_xxxx_xyzz, s_xxxx_xzzz,\
                                     s_xxxx_yyyy, s_xxxx_yyyz, s_xxxx_yyzz, s_xxxx_yzzz,\
                                     s_xxxx_zzzz, s_xxxy_xxxx, s_xxxy_xxxy, s_xxxy_xxxz,\
                                     s_xxxy_xxyy, s_xxxy_xxyz, s_xxxy_xxzz, s_xxxy_xyyy,\
                                     s_xxxy_xyyz, s_xxxy_xyzz, s_xxxy_xzzz, s_xxxy_yyyy,\
                                     s_xxxy_yyyz, s_xxxy_yyzz, s_xxxy_yzzz, s_xxxy_zzzz,\
                                     s_xxxz_xxxx, s_xxxz_xxxy, s_xxxz_xxxz, s_xxxz_xxyy,\
                                     s_xxxz_xxyz, s_xxxz_xxzz, s_xxxz_xyyy, s_xxxz_xyyz,\
                                     s_xxxz_xyzz, s_xxxz_xzzz, s_xxxz_yyyy, s_xxxz_yyyz,\
                                     s_xxxz_yyzz, s_xxxz_yzzz, s_xxxz_zzzz, s_xxyy_xxxx,\
                                     s_xxyy_xxxy, s_xxyy_xxxz, s_xxyy_xxyy, s_xxyy_xxyz,\
                                     s_xxyy_xxzz, s_xxyy_xyyy, s_xxyy_xyyz, s_xxyy_xyzz,\
                                     s_xxyy_xzzz, s_xxyy_yyyy, s_xxyy_yyyz, s_xxyy_yyzz,\
                                     s_xxyy_yzzz, s_xxyy_zzzz, s_xxyz_xxxx, s_xxyz_xxxy,\
                                     s_xxyz_xxxz, s_xxyz_xxyy, s_xxyz_xxyz, s_xxyz_xxzz,\
                                     s_xxyz_xyyy, s_xxyz_xyyz, s_xxyz_xyzz, s_xxyz_xzzz,\
                                     s_xxyz_yyyy, s_xxyz_yyyz, s_xxyz_yyzz, s_xxyz_yzzz,\
                                     s_xxyz_zzzz, s_xxzz_xxxx, s_xxzz_xxxy, s_xxzz_xxxz,\
                                     s_xxzz_xxyy, s_xxzz_xxyz, s_xxzz_xxzz, s_xxzz_xyyy,\
                                     s_xxzz_xyyz, s_xxzz_xyzz, s_xxzz_xzzz, s_xxzz_yyyy,\
                                     s_xxzz_yyyz, s_xxzz_yyzz, s_xxzz_yzzz, s_xxzz_zzzz,\
                                     s_xyyy_xxxx, s_xyyy_xxxy, s_xyyy_xxxz, s_xyyy_xxyy,\
                                     s_xyyy_xxyz, s_xyyy_xxzz, s_xyyy_xyyy, s_xyyy_xyyz,\
                                     s_xyyy_xyzz, s_xyyy_xzzz, s_xyyy_yyyy, s_xyyy_yyyz,\
                                     s_xyyy_yyzz, s_xyyy_yzzz, s_xyyy_zzzz, s_xyyz_xxxx,\
                                     s_xyyz_xxxy, s_xyyz_xxxz, s_xyyz_xxyy, s_xyyz_xxyz,\
                                     s_xyyz_xxzz, s_xyyz_xyyy, s_xyyz_xyyz, s_xyyz_xyzz,\
                                     s_xyyz_xzzz, s_xyyz_yyyy, s_xyyz_yyyz, s_xyyz_yyzz,\
                                     s_xyyz_yzzz, s_xyyz_zzzz, s_xyzz_xxxx, s_xyzz_xxxy,\
                                     s_xyzz_xxxz, s_xyzz_xxyy, s_xyzz_xxyz, s_xyzz_xxzz,\
                                     s_xyzz_xyyy, s_xyzz_xyyz, s_xyzz_xyzz, s_xyzz_xzzz,\
                                     s_xyzz_yyyy, s_xyzz_yyyz, s_xyzz_yyzz, s_xyzz_yzzz,\
                                     s_xyzz_zzzz, s_xzzz_xxxx, s_xzzz_xxxy, s_xzzz_xxxz,\
                                     s_xzzz_xxyy, s_xzzz_xxyz, s_xzzz_xxzz, s_xzzz_xyyy,\
                                     s_xzzz_xyyz, s_xzzz_xyzz, s_xzzz_xzzz, s_xzzz_yyyy,\
                                     s_xzzz_yyyz, s_xzzz_yyzz, s_xzzz_yzzz, s_xzzz_zzzz,\
                                     s_yyyy_xxxx, s_yyyy_xxxy, s_yyyy_xxxz, s_yyyy_xxyy,\
                                     s_yyyy_xxyz, s_yyyy_xxzz, s_yyyy_xyyy, s_yyyy_xyyz,\
                                     s_yyyy_xyzz, s_yyyy_xzzz, s_yyyy_yyyy, s_yyyy_yyyz,\
                                     s_yyyy_yyzz, s_yyyy_yzzz, s_yyyy_zzzz, s_yyyz_xxxx,\
                                     s_yyyz_xxxy, s_yyyz_xxxz, s_yyyz_xxyy, s_yyyz_xxyz,\
                                     s_yyyz_xxzz, s_yyyz_xyyy, s_yyyz_xyyz, s_yyyz_xyzz,\
                                     s_yyyz_xzzz, s_yyyz_yyyy, s_yyyz_yyyz, s_yyyz_yyzz,\
                                     s_yyyz_yzzz, s_yyyz_zzzz, s_yyzz_xxxx, s_yyzz_xxxy,\
                                     s_yyzz_xxxz, s_yyzz_xxyy, s_yyzz_xxyz, s_yyzz_xxzz,\
                                     s_yyzz_xyyy, s_yyzz_xyyz, s_yyzz_xyzz, s_yyzz_xzzz,\
                                     s_yyzz_yyyy, s_yyzz_yyyz, s_yyzz_yyzz, s_yyzz_yzzz,\
                                     s_yyzz_zzzz, s_yzzz_xxxx, s_yzzz_xxxy, s_yzzz_xxxz,\
                                     s_yzzz_xxyy, s_yzzz_xxyz, s_yzzz_xxzz, s_yzzz_xyyy,\
                                     s_yzzz_xyyz, s_yzzz_xyzz, s_yzzz_xzzz, s_yzzz_yyyy,\
                                     s_yzzz_yyyz, s_yzzz_yyzz, s_yzzz_yzzz, s_yzzz_zzzz,\
                                     s_zzzz_xxxx, s_zzzz_xxxy, s_zzzz_xxxz, s_zzzz_xxyy,\
                                     s_zzzz_xxyz, s_zzzz_xxzz, s_zzzz_xyyy, s_zzzz_xyyz,\
                                     s_zzzz_xyzz, s_zzzz_xzzz, s_zzzz_yyyy, s_zzzz_yyyz,\
                                     s_zzzz_yyzz, s_zzzz_yzzz, s_zzzz_zzzz, t_xxxx_xxxx,\
                                     t_xxxx_xxxy, t_xxxx_xxxz, t_xxxx_xxyy, t_xxxx_xxyz,\
                                     t_xxxx_xxzz, t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz,\
                                     t_xxxx_xzzz, t_xxxx_yyyy, t_xxxx_yyyz, t_xxxx_yyzz,\
                                     t_xxxx_yzzz, t_xxxx_zzzz, t_xxxy_xxxx, t_xxxy_xxxy,\
                                     t_xxxy_xxxz, t_xxxy_xxyy, t_xxxy_xxyz, t_xxxy_xxzz,\
                                     t_xxxy_xyyy, t_xxxy_xyyz, t_xxxy_xyzz, t_xxxy_xzzz,\
                                     t_xxxy_yyyy, t_xxxy_yyyz, t_xxxy_yyzz, t_xxxy_yzzz,\
                                     t_xxxy_zzzz, t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz,\
                                     t_xxxz_xxyy, t_xxxz_xxyz, t_xxxz_xxzz, t_xxxz_xyyy,\
                                     t_xxxz_xyyz, t_xxxz_xyzz, t_xxxz_xzzz, t_xxxz_yyyy,\
                                     t_xxxz_yyyz, t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz,\
                                     t_xxyy_xxxx, t_xxyy_xxxy, t_xxyy_xxxz, t_xxyy_xxyy,\
                                     t_xxyy_xxyz, t_xxyy_xxzz, t_xxyy_xyyy, t_xxyy_xyyz,\
                                     t_xxyy_xyzz, t_xxyy_xzzz, t_xxyy_yyyy, t_xxyy_yyyz,\
                                     t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz, t_xxyz_xxxx,\
                                     t_xxyz_xxxy, t_xxyz_xxxz, t_xxyz_xxyy, t_xxyz_xxyz,\
                                     t_xxyz_xxzz, t_xxyz_xyyy, t_xxyz_xyyz, t_xxyz_xyzz,\
                                     t_xxyz_xzzz, t_xxyz_yyyy, t_xxyz_yyyz, t_xxyz_yyzz,\
                                     t_xxyz_yzzz, t_xxyz_zzzz, t_xxzz_xxxx, t_xxzz_xxxy,\
                                     t_xxzz_xxxz, t_xxzz_xxyy, t_xxzz_xxyz, t_xxzz_xxzz,\
                                     t_xxzz_xyyy, t_xxzz_xyyz, t_xxzz_xyzz, t_xxzz_xzzz,\
                                     t_xxzz_yyyy, t_xxzz_yyyz, t_xxzz_yyzz, t_xxzz_yzzz,\
                                     t_xxzz_zzzz, t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz,\
                                     t_xyyy_xxyy, t_xyyy_xxyz, t_xyyy_xxzz, t_xyyy_xyyy,\
                                     t_xyyy_xyyz, t_xyyy_xyzz, t_xyyy_xzzz, t_xyyy_yyyy,\
                                     t_xyyy_yyyz, t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz,\
                                     t_xyyz_xxxx, t_xyyz_xxxy, t_xyyz_xxxz, t_xyyz_xxyy,\
                                     t_xyyz_xxyz, t_xyyz_xxzz, t_xyyz_xyyy, t_xyyz_xyyz,\
                                     t_xyyz_xyzz, t_xyyz_xzzz, t_xyyz_yyyy, t_xyyz_yyyz,\
                                     t_xyyz_yyzz, t_xyyz_yzzz, t_xyyz_zzzz, t_xyzz_xxxx,\
                                     t_xyzz_xxxy, t_xyzz_xxxz, t_xyzz_xxyy, t_xyzz_xxyz,\
                                     t_xyzz_xxzz, t_xyzz_xyyy, t_xyzz_xyyz, t_xyzz_xyzz,\
                                     t_xyzz_xzzz, t_xyzz_yyyy, t_xyzz_yyyz, t_xyzz_yyzz,\
                                     t_xyzz_yzzz, t_xyzz_zzzz, t_xzzz_xxxx, t_xzzz_xxxy,\
                                     t_xzzz_xxxz, t_xzzz_xxyy, t_xzzz_xxyz, t_xzzz_xxzz,\
                                     t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz, t_xzzz_xzzz,\
                                     t_xzzz_yyyy, t_xzzz_yyyz, t_xzzz_yyzz, t_xzzz_yzzz,\
                                     t_xzzz_zzzz, t_yyyy_xxxx, t_yyyy_xxxy, t_yyyy_xxxz,\
                                     t_yyyy_xxyy, t_yyyy_xxyz, t_yyyy_xxzz, t_yyyy_xyyy,\
                                     t_yyyy_xyyz, t_yyyy_xyzz, t_yyyy_xzzz, t_yyyy_yyyy,\
                                     t_yyyy_yyyz, t_yyyy_yyzz, t_yyyy_yzzz, t_yyyy_zzzz,\
                                     t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz, t_yyyz_xxyy,\
                                     t_yyyz_xxyz, t_yyyz_xxzz, t_yyyz_xyyy, t_yyyz_xyyz,\
                                     t_yyyz_xyzz, t_yyyz_xzzz, t_yyyz_yyyy, t_yyyz_yyyz,\
                                     t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz, t_yyzz_xxxx,\
                                     t_yyzz_xxxy, t_yyzz_xxxz, t_yyzz_xxyy, t_yyzz_xxyz,\
                                     t_yyzz_xxzz, t_yyzz_xyyy, t_yyzz_xyyz, t_yyzz_xyzz,\
                                     t_yyzz_xzzz, t_yyzz_yyyy, t_yyzz_yyyz, t_yyzz_yyzz,\
                                     t_yyzz_yzzz, t_yyzz_zzzz, t_yzzz_xxxx, t_yzzz_xxxy,\
                                     t_yzzz_xxxz, t_yzzz_xxyy, t_yzzz_xxyz, t_yzzz_xxzz,\
                                     t_yzzz_xyyy, t_yzzz_xyyz, t_yzzz_xyzz, t_yzzz_xzzz,\
                                     t_yzzz_yyyy, t_yzzz_yyyz, t_yzzz_yyzz, t_yzzz_yzzz,\
                                     t_yzzz_zzzz, t_zzzz_xxxx, t_zzzz_xxxy, t_zzzz_xxxz,\
                                     t_zzzz_xxyy, t_zzzz_xxyz, t_zzzz_xxzz, t_zzzz_xyyy,\
                                     t_zzzz_xyyz, t_zzzz_xyzz, t_zzzz_xzzz, t_zzzz_yyyy,\
                                     t_zzzz_yyyz, t_zzzz_yyzz, t_zzzz_yzzz, t_zzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                double f2z = 2.00 * fz[j];

                double f2a = 0.50 * fga[j];

                // leading x component

                double fr = pax[j];

                s_xxxx_xxxx[j] = fr * s_xxx_xxxx[j] + f2t * (3.0 * s_xx_xxxx[j] + 4.0 * s_xxx_xxx[j]);

                t_xxxx_xxxx[j] = fr * t_xxx_xxxx[j] + f2t * (3.0 * t_xx_xxxx[j] + 4.0 * t_xxx_xxx[j]) + f2z * (s_xxxx_xxxx[j] - f2a * 3.0 * s_xx_xxxx[j]);

                s_xxxx_xxxy[j] = fr * s_xxx_xxxy[j] + f2t * (3.0 * s_xx_xxxy[j] + 3.0 * s_xxx_xxy[j]);

                t_xxxx_xxxy[j] = fr * t_xxx_xxxy[j] + f2t * (3.0 * t_xx_xxxy[j] + 3.0 * t_xxx_xxy[j]) + f2z * (s_xxxx_xxxy[j] - f2a * 3.0 * s_xx_xxxy[j]);

                s_xxxx_xxxz[j] = fr * s_xxx_xxxz[j] + f2t * (3.0 * s_xx_xxxz[j] + 3.0 * s_xxx_xxz[j]);

                t_xxxx_xxxz[j] = fr * t_xxx_xxxz[j] + f2t * (3.0 * t_xx_xxxz[j] + 3.0 * t_xxx_xxz[j]) + f2z * (s_xxxx_xxxz[j] - f2a * 3.0 * s_xx_xxxz[j]);

                s_xxxx_xxyy[j] = fr * s_xxx_xxyy[j] + f2t * (3.0 * s_xx_xxyy[j] + 2.0 * s_xxx_xyy[j]);

                t_xxxx_xxyy[j] = fr * t_xxx_xxyy[j] + f2t * (3.0 * t_xx_xxyy[j] + 2.0 * t_xxx_xyy[j]) + f2z * (s_xxxx_xxyy[j] - f2a * 3.0 * s_xx_xxyy[j]);

                s_xxxx_xxyz[j] = fr * s_xxx_xxyz[j] + f2t * (3.0 * s_xx_xxyz[j] + 2.0 * s_xxx_xyz[j]);

                t_xxxx_xxyz[j] = fr * t_xxx_xxyz[j] + f2t * (3.0 * t_xx_xxyz[j] + 2.0 * t_xxx_xyz[j]) + f2z * (s_xxxx_xxyz[j] - f2a * 3.0 * s_xx_xxyz[j]);

                s_xxxx_xxzz[j] = fr * s_xxx_xxzz[j] + f2t * (3.0 * s_xx_xxzz[j] + 2.0 * s_xxx_xzz[j]);

                t_xxxx_xxzz[j] = fr * t_xxx_xxzz[j] + f2t * (3.0 * t_xx_xxzz[j] + 2.0 * t_xxx_xzz[j]) + f2z * (s_xxxx_xxzz[j] - f2a * 3.0 * s_xx_xxzz[j]);

                s_xxxx_xyyy[j] = fr * s_xxx_xyyy[j] + f2t * (3.0 * s_xx_xyyy[j] + s_xxx_yyy[j]);

                t_xxxx_xyyy[j] = fr * t_xxx_xyyy[j] + f2t * (3.0 * t_xx_xyyy[j] + t_xxx_yyy[j]) + f2z * (s_xxxx_xyyy[j] - f2a * 3.0 * s_xx_xyyy[j]);

                s_xxxx_xyyz[j] = fr * s_xxx_xyyz[j] + f2t * (3.0 * s_xx_xyyz[j] + s_xxx_yyz[j]);

                t_xxxx_xyyz[j] = fr * t_xxx_xyyz[j] + f2t * (3.0 * t_xx_xyyz[j] + t_xxx_yyz[j]) + f2z * (s_xxxx_xyyz[j] - f2a * 3.0 * s_xx_xyyz[j]);

                s_xxxx_xyzz[j] = fr * s_xxx_xyzz[j] + f2t * (3.0 * s_xx_xyzz[j] + s_xxx_yzz[j]);

                t_xxxx_xyzz[j] = fr * t_xxx_xyzz[j] + f2t * (3.0 * t_xx_xyzz[j] + t_xxx_yzz[j]) + f2z * (s_xxxx_xyzz[j] - f2a * 3.0 * s_xx_xyzz[j]);

                s_xxxx_xzzz[j] = fr * s_xxx_xzzz[j] + f2t * (3.0 * s_xx_xzzz[j] + s_xxx_zzz[j]);

                t_xxxx_xzzz[j] = fr * t_xxx_xzzz[j] + f2t * (3.0 * t_xx_xzzz[j] + t_xxx_zzz[j]) + f2z * (s_xxxx_xzzz[j] - f2a * 3.0 * s_xx_xzzz[j]);

                s_xxxx_yyyy[j] = fr * s_xxx_yyyy[j] + f2t * 3.0 * s_xx_yyyy[j];

                t_xxxx_yyyy[j] = fr * t_xxx_yyyy[j] + f2t * 3.0 * t_xx_yyyy[j] + f2z * (s_xxxx_yyyy[j] - f2a * 3.0 * s_xx_yyyy[j]);

                s_xxxx_yyyz[j] = fr * s_xxx_yyyz[j] + f2t * 3.0 * s_xx_yyyz[j];

                t_xxxx_yyyz[j] = fr * t_xxx_yyyz[j] + f2t * 3.0 * t_xx_yyyz[j] + f2z * (s_xxxx_yyyz[j] - f2a * 3.0 * s_xx_yyyz[j]);

                s_xxxx_yyzz[j] = fr * s_xxx_yyzz[j] + f2t * 3.0 * s_xx_yyzz[j];

                t_xxxx_yyzz[j] = fr * t_xxx_yyzz[j] + f2t * 3.0 * t_xx_yyzz[j] + f2z * (s_xxxx_yyzz[j] - f2a * 3.0 * s_xx_yyzz[j]);

                s_xxxx_yzzz[j] = fr * s_xxx_yzzz[j] + f2t * 3.0 * s_xx_yzzz[j];

                t_xxxx_yzzz[j] = fr * t_xxx_yzzz[j] + f2t * 3.0 * t_xx_yzzz[j] + f2z * (s_xxxx_yzzz[j] - f2a * 3.0 * s_xx_yzzz[j]);

                s_xxxx_zzzz[j] = fr * s_xxx_zzzz[j] + f2t * 3.0 * s_xx_zzzz[j];

                t_xxxx_zzzz[j] = fr * t_xxx_zzzz[j] + f2t * 3.0 * t_xx_zzzz[j] + f2z * (s_xxxx_zzzz[j] - f2a * 3.0 * s_xx_zzzz[j]);

                s_xxxy_xxxx[j] = fr * s_xxy_xxxx[j] + f2t * (2.0 * s_xy_xxxx[j] + 4.0 * s_xxy_xxx[j]);

                t_xxxy_xxxx[j] = fr * t_xxy_xxxx[j] + f2t * (2.0 * t_xy_xxxx[j] + 4.0 * t_xxy_xxx[j]) + f2z * (s_xxxy_xxxx[j] - f2a * 2.0 * s_xy_xxxx[j]);

                s_xxxy_xxxy[j] = fr * s_xxy_xxxy[j] + f2t * (2.0 * s_xy_xxxy[j] + 3.0 * s_xxy_xxy[j]);

                t_xxxy_xxxy[j] = fr * t_xxy_xxxy[j] + f2t * (2.0 * t_xy_xxxy[j] + 3.0 * t_xxy_xxy[j]) + f2z * (s_xxxy_xxxy[j] - f2a * 2.0 * s_xy_xxxy[j]);

                s_xxxy_xxxz[j] = fr * s_xxy_xxxz[j] + f2t * (2.0 * s_xy_xxxz[j] + 3.0 * s_xxy_xxz[j]);

                t_xxxy_xxxz[j] = fr * t_xxy_xxxz[j] + f2t * (2.0 * t_xy_xxxz[j] + 3.0 * t_xxy_xxz[j]) + f2z * (s_xxxy_xxxz[j] - f2a * 2.0 * s_xy_xxxz[j]);

                s_xxxy_xxyy[j] = fr * s_xxy_xxyy[j] + f2t * (2.0 * s_xy_xxyy[j] + 2.0 * s_xxy_xyy[j]);

                t_xxxy_xxyy[j] = fr * t_xxy_xxyy[j] + f2t * (2.0 * t_xy_xxyy[j] + 2.0 * t_xxy_xyy[j]) + f2z * (s_xxxy_xxyy[j] - f2a * 2.0 * s_xy_xxyy[j]);

                s_xxxy_xxyz[j] = fr * s_xxy_xxyz[j] + f2t * (2.0 * s_xy_xxyz[j] + 2.0 * s_xxy_xyz[j]);

                t_xxxy_xxyz[j] = fr * t_xxy_xxyz[j] + f2t * (2.0 * t_xy_xxyz[j] + 2.0 * t_xxy_xyz[j]) + f2z * (s_xxxy_xxyz[j] - f2a * 2.0 * s_xy_xxyz[j]);

                s_xxxy_xxzz[j] = fr * s_xxy_xxzz[j] + f2t * (2.0 * s_xy_xxzz[j] + 2.0 * s_xxy_xzz[j]);

                t_xxxy_xxzz[j] = fr * t_xxy_xxzz[j] + f2t * (2.0 * t_xy_xxzz[j] + 2.0 * t_xxy_xzz[j]) + f2z * (s_xxxy_xxzz[j] - f2a * 2.0 * s_xy_xxzz[j]);

                s_xxxy_xyyy[j] = fr * s_xxy_xyyy[j] + f2t * (2.0 * s_xy_xyyy[j] + s_xxy_yyy[j]);

                t_xxxy_xyyy[j] = fr * t_xxy_xyyy[j] + f2t * (2.0 * t_xy_xyyy[j] + t_xxy_yyy[j]) + f2z * (s_xxxy_xyyy[j] - f2a * 2.0 * s_xy_xyyy[j]);

                s_xxxy_xyyz[j] = fr * s_xxy_xyyz[j] + f2t * (2.0 * s_xy_xyyz[j] + s_xxy_yyz[j]);

                t_xxxy_xyyz[j] = fr * t_xxy_xyyz[j] + f2t * (2.0 * t_xy_xyyz[j] + t_xxy_yyz[j]) + f2z * (s_xxxy_xyyz[j] - f2a * 2.0 * s_xy_xyyz[j]);

                s_xxxy_xyzz[j] = fr * s_xxy_xyzz[j] + f2t * (2.0 * s_xy_xyzz[j] + s_xxy_yzz[j]);

                t_xxxy_xyzz[j] = fr * t_xxy_xyzz[j] + f2t * (2.0 * t_xy_xyzz[j] + t_xxy_yzz[j]) + f2z * (s_xxxy_xyzz[j] - f2a * 2.0 * s_xy_xyzz[j]);

                s_xxxy_xzzz[j] = fr * s_xxy_xzzz[j] + f2t * (2.0 * s_xy_xzzz[j] + s_xxy_zzz[j]);

                t_xxxy_xzzz[j] = fr * t_xxy_xzzz[j] + f2t * (2.0 * t_xy_xzzz[j] + t_xxy_zzz[j]) + f2z * (s_xxxy_xzzz[j] - f2a * 2.0 * s_xy_xzzz[j]);

                s_xxxy_yyyy[j] = fr * s_xxy_yyyy[j] + f2t * 2.0 * s_xy_yyyy[j];

                t_xxxy_yyyy[j] = fr * t_xxy_yyyy[j] + f2t * 2.0 * t_xy_yyyy[j] + f2z * (s_xxxy_yyyy[j] - f2a * 2.0 * s_xy_yyyy[j]);

                s_xxxy_yyyz[j] = fr * s_xxy_yyyz[j] + f2t * 2.0 * s_xy_yyyz[j];

                t_xxxy_yyyz[j] = fr * t_xxy_yyyz[j] + f2t * 2.0 * t_xy_yyyz[j] + f2z * (s_xxxy_yyyz[j] - f2a * 2.0 * s_xy_yyyz[j]);

                s_xxxy_yyzz[j] = fr * s_xxy_yyzz[j] + f2t * 2.0 * s_xy_yyzz[j];

                t_xxxy_yyzz[j] = fr * t_xxy_yyzz[j] + f2t * 2.0 * t_xy_yyzz[j] + f2z * (s_xxxy_yyzz[j] - f2a * 2.0 * s_xy_yyzz[j]);

                s_xxxy_yzzz[j] = fr * s_xxy_yzzz[j] + f2t * 2.0 * s_xy_yzzz[j];

                t_xxxy_yzzz[j] = fr * t_xxy_yzzz[j] + f2t * 2.0 * t_xy_yzzz[j] + f2z * (s_xxxy_yzzz[j] - f2a * 2.0 * s_xy_yzzz[j]);

                s_xxxy_zzzz[j] = fr * s_xxy_zzzz[j] + f2t * 2.0 * s_xy_zzzz[j];

                t_xxxy_zzzz[j] = fr * t_xxy_zzzz[j] + f2t * 2.0 * t_xy_zzzz[j] + f2z * (s_xxxy_zzzz[j] - f2a * 2.0 * s_xy_zzzz[j]);

                s_xxxz_xxxx[j] = fr * s_xxz_xxxx[j] + f2t * (2.0 * s_xz_xxxx[j] + 4.0 * s_xxz_xxx[j]);

                t_xxxz_xxxx[j] = fr * t_xxz_xxxx[j] + f2t * (2.0 * t_xz_xxxx[j] + 4.0 * t_xxz_xxx[j]) + f2z * (s_xxxz_xxxx[j] - f2a * 2.0 * s_xz_xxxx[j]);

                s_xxxz_xxxy[j] = fr * s_xxz_xxxy[j] + f2t * (2.0 * s_xz_xxxy[j] + 3.0 * s_xxz_xxy[j]);

                t_xxxz_xxxy[j] = fr * t_xxz_xxxy[j] + f2t * (2.0 * t_xz_xxxy[j] + 3.0 * t_xxz_xxy[j]) + f2z * (s_xxxz_xxxy[j] - f2a * 2.0 * s_xz_xxxy[j]);

                s_xxxz_xxxz[j] = fr * s_xxz_xxxz[j] + f2t * (2.0 * s_xz_xxxz[j] + 3.0 * s_xxz_xxz[j]);

                t_xxxz_xxxz[j] = fr * t_xxz_xxxz[j] + f2t * (2.0 * t_xz_xxxz[j] + 3.0 * t_xxz_xxz[j]) + f2z * (s_xxxz_xxxz[j] - f2a * 2.0 * s_xz_xxxz[j]);

                s_xxxz_xxyy[j] = fr * s_xxz_xxyy[j] + f2t * (2.0 * s_xz_xxyy[j] + 2.0 * s_xxz_xyy[j]);

                t_xxxz_xxyy[j] = fr * t_xxz_xxyy[j] + f2t * (2.0 * t_xz_xxyy[j] + 2.0 * t_xxz_xyy[j]) + f2z * (s_xxxz_xxyy[j] - f2a * 2.0 * s_xz_xxyy[j]);

                s_xxxz_xxyz[j] = fr * s_xxz_xxyz[j] + f2t * (2.0 * s_xz_xxyz[j] + 2.0 * s_xxz_xyz[j]);

                t_xxxz_xxyz[j] = fr * t_xxz_xxyz[j] + f2t * (2.0 * t_xz_xxyz[j] + 2.0 * t_xxz_xyz[j]) + f2z * (s_xxxz_xxyz[j] - f2a * 2.0 * s_xz_xxyz[j]);

                s_xxxz_xxzz[j] = fr * s_xxz_xxzz[j] + f2t * (2.0 * s_xz_xxzz[j] + 2.0 * s_xxz_xzz[j]);

                t_xxxz_xxzz[j] = fr * t_xxz_xxzz[j] + f2t * (2.0 * t_xz_xxzz[j] + 2.0 * t_xxz_xzz[j]) + f2z * (s_xxxz_xxzz[j] - f2a * 2.0 * s_xz_xxzz[j]);

                s_xxxz_xyyy[j] = fr * s_xxz_xyyy[j] + f2t * (2.0 * s_xz_xyyy[j] + s_xxz_yyy[j]);

                t_xxxz_xyyy[j] = fr * t_xxz_xyyy[j] + f2t * (2.0 * t_xz_xyyy[j] + t_xxz_yyy[j]) + f2z * (s_xxxz_xyyy[j] - f2a * 2.0 * s_xz_xyyy[j]);

                s_xxxz_xyyz[j] = fr * s_xxz_xyyz[j] + f2t * (2.0 * s_xz_xyyz[j] + s_xxz_yyz[j]);

                t_xxxz_xyyz[j] = fr * t_xxz_xyyz[j] + f2t * (2.0 * t_xz_xyyz[j] + t_xxz_yyz[j]) + f2z * (s_xxxz_xyyz[j] - f2a * 2.0 * s_xz_xyyz[j]);

                s_xxxz_xyzz[j] = fr * s_xxz_xyzz[j] + f2t * (2.0 * s_xz_xyzz[j] + s_xxz_yzz[j]);

                t_xxxz_xyzz[j] = fr * t_xxz_xyzz[j] + f2t * (2.0 * t_xz_xyzz[j] + t_xxz_yzz[j]) + f2z * (s_xxxz_xyzz[j] - f2a * 2.0 * s_xz_xyzz[j]);

                s_xxxz_xzzz[j] = fr * s_xxz_xzzz[j] + f2t * (2.0 * s_xz_xzzz[j] + s_xxz_zzz[j]);

                t_xxxz_xzzz[j] = fr * t_xxz_xzzz[j] + f2t * (2.0 * t_xz_xzzz[j] + t_xxz_zzz[j]) + f2z * (s_xxxz_xzzz[j] - f2a * 2.0 * s_xz_xzzz[j]);

                s_xxxz_yyyy[j] = fr * s_xxz_yyyy[j] + f2t * 2.0 * s_xz_yyyy[j];

                t_xxxz_yyyy[j] = fr * t_xxz_yyyy[j] + f2t * 2.0 * t_xz_yyyy[j] + f2z * (s_xxxz_yyyy[j] - f2a * 2.0 * s_xz_yyyy[j]);

                s_xxxz_yyyz[j] = fr * s_xxz_yyyz[j] + f2t * 2.0 * s_xz_yyyz[j];

                t_xxxz_yyyz[j] = fr * t_xxz_yyyz[j] + f2t * 2.0 * t_xz_yyyz[j] + f2z * (s_xxxz_yyyz[j] - f2a * 2.0 * s_xz_yyyz[j]);

                s_xxxz_yyzz[j] = fr * s_xxz_yyzz[j] + f2t * 2.0 * s_xz_yyzz[j];

                t_xxxz_yyzz[j] = fr * t_xxz_yyzz[j] + f2t * 2.0 * t_xz_yyzz[j] + f2z * (s_xxxz_yyzz[j] - f2a * 2.0 * s_xz_yyzz[j]);

                s_xxxz_yzzz[j] = fr * s_xxz_yzzz[j] + f2t * 2.0 * s_xz_yzzz[j];

                t_xxxz_yzzz[j] = fr * t_xxz_yzzz[j] + f2t * 2.0 * t_xz_yzzz[j] + f2z * (s_xxxz_yzzz[j] - f2a * 2.0 * s_xz_yzzz[j]);

                s_xxxz_zzzz[j] = fr * s_xxz_zzzz[j] + f2t * 2.0 * s_xz_zzzz[j];

                t_xxxz_zzzz[j] = fr * t_xxz_zzzz[j] + f2t * 2.0 * t_xz_zzzz[j] + f2z * (s_xxxz_zzzz[j] - f2a * 2.0 * s_xz_zzzz[j]);

                s_xxyy_xxxx[j] = fr * s_xyy_xxxx[j] + f2t * (s_yy_xxxx[j] + 4.0 * s_xyy_xxx[j]);

                t_xxyy_xxxx[j] = fr * t_xyy_xxxx[j] + f2t * (t_yy_xxxx[j] + 4.0 * t_xyy_xxx[j]) + f2z * (s_xxyy_xxxx[j] - f2a * s_yy_xxxx[j]);

                s_xxyy_xxxy[j] = fr * s_xyy_xxxy[j] + f2t * (s_yy_xxxy[j] + 3.0 * s_xyy_xxy[j]);

                t_xxyy_xxxy[j] = fr * t_xyy_xxxy[j] + f2t * (t_yy_xxxy[j] + 3.0 * t_xyy_xxy[j]) + f2z * (s_xxyy_xxxy[j] - f2a * s_yy_xxxy[j]);

                s_xxyy_xxxz[j] = fr * s_xyy_xxxz[j] + f2t * (s_yy_xxxz[j] + 3.0 * s_xyy_xxz[j]);

                t_xxyy_xxxz[j] = fr * t_xyy_xxxz[j] + f2t * (t_yy_xxxz[j] + 3.0 * t_xyy_xxz[j]) + f2z * (s_xxyy_xxxz[j] - f2a * s_yy_xxxz[j]);

                s_xxyy_xxyy[j] = fr * s_xyy_xxyy[j] + f2t * (s_yy_xxyy[j] + 2.0 * s_xyy_xyy[j]);

                t_xxyy_xxyy[j] = fr * t_xyy_xxyy[j] + f2t * (t_yy_xxyy[j] + 2.0 * t_xyy_xyy[j]) + f2z * (s_xxyy_xxyy[j] - f2a * s_yy_xxyy[j]);

                s_xxyy_xxyz[j] = fr * s_xyy_xxyz[j] + f2t * (s_yy_xxyz[j] + 2.0 * s_xyy_xyz[j]);

                t_xxyy_xxyz[j] = fr * t_xyy_xxyz[j] + f2t * (t_yy_xxyz[j] + 2.0 * t_xyy_xyz[j]) + f2z * (s_xxyy_xxyz[j] - f2a * s_yy_xxyz[j]);

                s_xxyy_xxzz[j] = fr * s_xyy_xxzz[j] + f2t * (s_yy_xxzz[j] + 2.0 * s_xyy_xzz[j]);

                t_xxyy_xxzz[j] = fr * t_xyy_xxzz[j] + f2t * (t_yy_xxzz[j] + 2.0 * t_xyy_xzz[j]) + f2z * (s_xxyy_xxzz[j] - f2a * s_yy_xxzz[j]);

                s_xxyy_xyyy[j] = fr * s_xyy_xyyy[j] + f2t * (s_yy_xyyy[j] + s_xyy_yyy[j]);

                t_xxyy_xyyy[j] = fr * t_xyy_xyyy[j] + f2t * (t_yy_xyyy[j] + t_xyy_yyy[j]) + f2z * (s_xxyy_xyyy[j] - f2a * s_yy_xyyy[j]);

                s_xxyy_xyyz[j] = fr * s_xyy_xyyz[j] + f2t * (s_yy_xyyz[j] + s_xyy_yyz[j]);

                t_xxyy_xyyz[j] = fr * t_xyy_xyyz[j] + f2t * (t_yy_xyyz[j] + t_xyy_yyz[j]) + f2z * (s_xxyy_xyyz[j] - f2a * s_yy_xyyz[j]);

                s_xxyy_xyzz[j] = fr * s_xyy_xyzz[j] + f2t * (s_yy_xyzz[j] + s_xyy_yzz[j]);

                t_xxyy_xyzz[j] = fr * t_xyy_xyzz[j] + f2t * (t_yy_xyzz[j] + t_xyy_yzz[j]) + f2z * (s_xxyy_xyzz[j] - f2a * s_yy_xyzz[j]);

                s_xxyy_xzzz[j] = fr * s_xyy_xzzz[j] + f2t * (s_yy_xzzz[j] + s_xyy_zzz[j]);

                t_xxyy_xzzz[j] = fr * t_xyy_xzzz[j] + f2t * (t_yy_xzzz[j] + t_xyy_zzz[j]) + f2z * (s_xxyy_xzzz[j] - f2a * s_yy_xzzz[j]);

                s_xxyy_yyyy[j] = fr * s_xyy_yyyy[j] + f2t * s_yy_yyyy[j];

                t_xxyy_yyyy[j] = fr * t_xyy_yyyy[j] + f2t * t_yy_yyyy[j] + f2z * (s_xxyy_yyyy[j] - f2a * s_yy_yyyy[j]);

                s_xxyy_yyyz[j] = fr * s_xyy_yyyz[j] + f2t * s_yy_yyyz[j];

                t_xxyy_yyyz[j] = fr * t_xyy_yyyz[j] + f2t * t_yy_yyyz[j] + f2z * (s_xxyy_yyyz[j] - f2a * s_yy_yyyz[j]);

                s_xxyy_yyzz[j] = fr * s_xyy_yyzz[j] + f2t * s_yy_yyzz[j];

                t_xxyy_yyzz[j] = fr * t_xyy_yyzz[j] + f2t * t_yy_yyzz[j] + f2z * (s_xxyy_yyzz[j] - f2a * s_yy_yyzz[j]);

                s_xxyy_yzzz[j] = fr * s_xyy_yzzz[j] + f2t * s_yy_yzzz[j];

                t_xxyy_yzzz[j] = fr * t_xyy_yzzz[j] + f2t * t_yy_yzzz[j] + f2z * (s_xxyy_yzzz[j] - f2a * s_yy_yzzz[j]);

                s_xxyy_zzzz[j] = fr * s_xyy_zzzz[j] + f2t * s_yy_zzzz[j];

                t_xxyy_zzzz[j] = fr * t_xyy_zzzz[j] + f2t * t_yy_zzzz[j] + f2z * (s_xxyy_zzzz[j] - f2a * s_yy_zzzz[j]);

                s_xxyz_xxxx[j] = fr * s_xyz_xxxx[j] + f2t * (s_yz_xxxx[j] + 4.0 * s_xyz_xxx[j]);

                t_xxyz_xxxx[j] = fr * t_xyz_xxxx[j] + f2t * (t_yz_xxxx[j] + 4.0 * t_xyz_xxx[j]) + f2z * (s_xxyz_xxxx[j] - f2a * s_yz_xxxx[j]);

                s_xxyz_xxxy[j] = fr * s_xyz_xxxy[j] + f2t * (s_yz_xxxy[j] + 3.0 * s_xyz_xxy[j]);

                t_xxyz_xxxy[j] = fr * t_xyz_xxxy[j] + f2t * (t_yz_xxxy[j] + 3.0 * t_xyz_xxy[j]) + f2z * (s_xxyz_xxxy[j] - f2a * s_yz_xxxy[j]);

                s_xxyz_xxxz[j] = fr * s_xyz_xxxz[j] + f2t * (s_yz_xxxz[j] + 3.0 * s_xyz_xxz[j]);

                t_xxyz_xxxz[j] = fr * t_xyz_xxxz[j] + f2t * (t_yz_xxxz[j] + 3.0 * t_xyz_xxz[j]) + f2z * (s_xxyz_xxxz[j] - f2a * s_yz_xxxz[j]);

                s_xxyz_xxyy[j] = fr * s_xyz_xxyy[j] + f2t * (s_yz_xxyy[j] + 2.0 * s_xyz_xyy[j]);

                t_xxyz_xxyy[j] = fr * t_xyz_xxyy[j] + f2t * (t_yz_xxyy[j] + 2.0 * t_xyz_xyy[j]) + f2z * (s_xxyz_xxyy[j] - f2a * s_yz_xxyy[j]);

                s_xxyz_xxyz[j] = fr * s_xyz_xxyz[j] + f2t * (s_yz_xxyz[j] + 2.0 * s_xyz_xyz[j]);

                t_xxyz_xxyz[j] = fr * t_xyz_xxyz[j] + f2t * (t_yz_xxyz[j] + 2.0 * t_xyz_xyz[j]) + f2z * (s_xxyz_xxyz[j] - f2a * s_yz_xxyz[j]);

                s_xxyz_xxzz[j] = fr * s_xyz_xxzz[j] + f2t * (s_yz_xxzz[j] + 2.0 * s_xyz_xzz[j]);

                t_xxyz_xxzz[j] = fr * t_xyz_xxzz[j] + f2t * (t_yz_xxzz[j] + 2.0 * t_xyz_xzz[j]) + f2z * (s_xxyz_xxzz[j] - f2a * s_yz_xxzz[j]);

                s_xxyz_xyyy[j] = fr * s_xyz_xyyy[j] + f2t * (s_yz_xyyy[j] + s_xyz_yyy[j]);

                t_xxyz_xyyy[j] = fr * t_xyz_xyyy[j] + f2t * (t_yz_xyyy[j] + t_xyz_yyy[j]) + f2z * (s_xxyz_xyyy[j] - f2a * s_yz_xyyy[j]);

                s_xxyz_xyyz[j] = fr * s_xyz_xyyz[j] + f2t * (s_yz_xyyz[j] + s_xyz_yyz[j]);

                t_xxyz_xyyz[j] = fr * t_xyz_xyyz[j] + f2t * (t_yz_xyyz[j] + t_xyz_yyz[j]) + f2z * (s_xxyz_xyyz[j] - f2a * s_yz_xyyz[j]);

                s_xxyz_xyzz[j] = fr * s_xyz_xyzz[j] + f2t * (s_yz_xyzz[j] + s_xyz_yzz[j]);

                t_xxyz_xyzz[j] = fr * t_xyz_xyzz[j] + f2t * (t_yz_xyzz[j] + t_xyz_yzz[j]) + f2z * (s_xxyz_xyzz[j] - f2a * s_yz_xyzz[j]);

                s_xxyz_xzzz[j] = fr * s_xyz_xzzz[j] + f2t * (s_yz_xzzz[j] + s_xyz_zzz[j]);

                t_xxyz_xzzz[j] = fr * t_xyz_xzzz[j] + f2t * (t_yz_xzzz[j] + t_xyz_zzz[j]) + f2z * (s_xxyz_xzzz[j] - f2a * s_yz_xzzz[j]);

                s_xxyz_yyyy[j] = fr * s_xyz_yyyy[j] + f2t * s_yz_yyyy[j];

                t_xxyz_yyyy[j] = fr * t_xyz_yyyy[j] + f2t * t_yz_yyyy[j] + f2z * (s_xxyz_yyyy[j] - f2a * s_yz_yyyy[j]);

                s_xxyz_yyyz[j] = fr * s_xyz_yyyz[j] + f2t * s_yz_yyyz[j];

                t_xxyz_yyyz[j] = fr * t_xyz_yyyz[j] + f2t * t_yz_yyyz[j] + f2z * (s_xxyz_yyyz[j] - f2a * s_yz_yyyz[j]);

                s_xxyz_yyzz[j] = fr * s_xyz_yyzz[j] + f2t * s_yz_yyzz[j];

                t_xxyz_yyzz[j] = fr * t_xyz_yyzz[j] + f2t * t_yz_yyzz[j] + f2z * (s_xxyz_yyzz[j] - f2a * s_yz_yyzz[j]);

                s_xxyz_yzzz[j] = fr * s_xyz_yzzz[j] + f2t * s_yz_yzzz[j];

                t_xxyz_yzzz[j] = fr * t_xyz_yzzz[j] + f2t * t_yz_yzzz[j] + f2z * (s_xxyz_yzzz[j] - f2a * s_yz_yzzz[j]);

                s_xxyz_zzzz[j] = fr * s_xyz_zzzz[j] + f2t * s_yz_zzzz[j];

                t_xxyz_zzzz[j] = fr * t_xyz_zzzz[j] + f2t * t_yz_zzzz[j] + f2z * (s_xxyz_zzzz[j] - f2a * s_yz_zzzz[j]);

                s_xxzz_xxxx[j] = fr * s_xzz_xxxx[j] + f2t * (s_zz_xxxx[j] + 4.0 * s_xzz_xxx[j]);

                t_xxzz_xxxx[j] = fr * t_xzz_xxxx[j] + f2t * (t_zz_xxxx[j] + 4.0 * t_xzz_xxx[j]) + f2z * (s_xxzz_xxxx[j] - f2a * s_zz_xxxx[j]);

                s_xxzz_xxxy[j] = fr * s_xzz_xxxy[j] + f2t * (s_zz_xxxy[j] + 3.0 * s_xzz_xxy[j]);

                t_xxzz_xxxy[j] = fr * t_xzz_xxxy[j] + f2t * (t_zz_xxxy[j] + 3.0 * t_xzz_xxy[j]) + f2z * (s_xxzz_xxxy[j] - f2a * s_zz_xxxy[j]);

                s_xxzz_xxxz[j] = fr * s_xzz_xxxz[j] + f2t * (s_zz_xxxz[j] + 3.0 * s_xzz_xxz[j]);

                t_xxzz_xxxz[j] = fr * t_xzz_xxxz[j] + f2t * (t_zz_xxxz[j] + 3.0 * t_xzz_xxz[j]) + f2z * (s_xxzz_xxxz[j] - f2a * s_zz_xxxz[j]);

                s_xxzz_xxyy[j] = fr * s_xzz_xxyy[j] + f2t * (s_zz_xxyy[j] + 2.0 * s_xzz_xyy[j]);

                t_xxzz_xxyy[j] = fr * t_xzz_xxyy[j] + f2t * (t_zz_xxyy[j] + 2.0 * t_xzz_xyy[j]) + f2z * (s_xxzz_xxyy[j] - f2a * s_zz_xxyy[j]);

                s_xxzz_xxyz[j] = fr * s_xzz_xxyz[j] + f2t * (s_zz_xxyz[j] + 2.0 * s_xzz_xyz[j]);

                t_xxzz_xxyz[j] = fr * t_xzz_xxyz[j] + f2t * (t_zz_xxyz[j] + 2.0 * t_xzz_xyz[j]) + f2z * (s_xxzz_xxyz[j] - f2a * s_zz_xxyz[j]);

                s_xxzz_xxzz[j] = fr * s_xzz_xxzz[j] + f2t * (s_zz_xxzz[j] + 2.0 * s_xzz_xzz[j]);

                t_xxzz_xxzz[j] = fr * t_xzz_xxzz[j] + f2t * (t_zz_xxzz[j] + 2.0 * t_xzz_xzz[j]) + f2z * (s_xxzz_xxzz[j] - f2a * s_zz_xxzz[j]);

                s_xxzz_xyyy[j] = fr * s_xzz_xyyy[j] + f2t * (s_zz_xyyy[j] + s_xzz_yyy[j]);

                t_xxzz_xyyy[j] = fr * t_xzz_xyyy[j] + f2t * (t_zz_xyyy[j] + t_xzz_yyy[j]) + f2z * (s_xxzz_xyyy[j] - f2a * s_zz_xyyy[j]);

                s_xxzz_xyyz[j] = fr * s_xzz_xyyz[j] + f2t * (s_zz_xyyz[j] + s_xzz_yyz[j]);

                t_xxzz_xyyz[j] = fr * t_xzz_xyyz[j] + f2t * (t_zz_xyyz[j] + t_xzz_yyz[j]) + f2z * (s_xxzz_xyyz[j] - f2a * s_zz_xyyz[j]);

                s_xxzz_xyzz[j] = fr * s_xzz_xyzz[j] + f2t * (s_zz_xyzz[j] + s_xzz_yzz[j]);

                t_xxzz_xyzz[j] = fr * t_xzz_xyzz[j] + f2t * (t_zz_xyzz[j] + t_xzz_yzz[j]) + f2z * (s_xxzz_xyzz[j] - f2a * s_zz_xyzz[j]);

                s_xxzz_xzzz[j] = fr * s_xzz_xzzz[j] + f2t * (s_zz_xzzz[j] + s_xzz_zzz[j]);

                t_xxzz_xzzz[j] = fr * t_xzz_xzzz[j] + f2t * (t_zz_xzzz[j] + t_xzz_zzz[j]) + f2z * (s_xxzz_xzzz[j] - f2a * s_zz_xzzz[j]);

                s_xxzz_yyyy[j] = fr * s_xzz_yyyy[j] + f2t * s_zz_yyyy[j];

                t_xxzz_yyyy[j] = fr * t_xzz_yyyy[j] + f2t * t_zz_yyyy[j] + f2z * (s_xxzz_yyyy[j] - f2a * s_zz_yyyy[j]);

                s_xxzz_yyyz[j] = fr * s_xzz_yyyz[j] + f2t * s_zz_yyyz[j];

                t_xxzz_yyyz[j] = fr * t_xzz_yyyz[j] + f2t * t_zz_yyyz[j] + f2z * (s_xxzz_yyyz[j] - f2a * s_zz_yyyz[j]);

                s_xxzz_yyzz[j] = fr * s_xzz_yyzz[j] + f2t * s_zz_yyzz[j];

                t_xxzz_yyzz[j] = fr * t_xzz_yyzz[j] + f2t * t_zz_yyzz[j] + f2z * (s_xxzz_yyzz[j] - f2a * s_zz_yyzz[j]);

                s_xxzz_yzzz[j] = fr * s_xzz_yzzz[j] + f2t * s_zz_yzzz[j];

                t_xxzz_yzzz[j] = fr * t_xzz_yzzz[j] + f2t * t_zz_yzzz[j] + f2z * (s_xxzz_yzzz[j] - f2a * s_zz_yzzz[j]);

                s_xxzz_zzzz[j] = fr * s_xzz_zzzz[j] + f2t * s_zz_zzzz[j];

                t_xxzz_zzzz[j] = fr * t_xzz_zzzz[j] + f2t * t_zz_zzzz[j] + f2z * (s_xxzz_zzzz[j] - f2a * s_zz_zzzz[j]);

                s_xyyy_xxxx[j] = fr * s_yyy_xxxx[j] + f2t * 4.0 * s_yyy_xxx[j];

                t_xyyy_xxxx[j] = fr * t_yyy_xxxx[j] + f2t * 4.0 * t_yyy_xxx[j] + f2z * s_xyyy_xxxx[j];

                s_xyyy_xxxy[j] = fr * s_yyy_xxxy[j] + f2t * 3.0 * s_yyy_xxy[j];

                t_xyyy_xxxy[j] = fr * t_yyy_xxxy[j] + f2t * 3.0 * t_yyy_xxy[j] + f2z * s_xyyy_xxxy[j];

                s_xyyy_xxxz[j] = fr * s_yyy_xxxz[j] + f2t * 3.0 * s_yyy_xxz[j];

                t_xyyy_xxxz[j] = fr * t_yyy_xxxz[j] + f2t * 3.0 * t_yyy_xxz[j] + f2z * s_xyyy_xxxz[j];

                s_xyyy_xxyy[j] = fr * s_yyy_xxyy[j] + f2t * 2.0 * s_yyy_xyy[j];

                t_xyyy_xxyy[j] = fr * t_yyy_xxyy[j] + f2t * 2.0 * t_yyy_xyy[j] + f2z * s_xyyy_xxyy[j];

                s_xyyy_xxyz[j] = fr * s_yyy_xxyz[j] + f2t * 2.0 * s_yyy_xyz[j];

                t_xyyy_xxyz[j] = fr * t_yyy_xxyz[j] + f2t * 2.0 * t_yyy_xyz[j] + f2z * s_xyyy_xxyz[j];

                s_xyyy_xxzz[j] = fr * s_yyy_xxzz[j] + f2t * 2.0 * s_yyy_xzz[j];

                t_xyyy_xxzz[j] = fr * t_yyy_xxzz[j] + f2t * 2.0 * t_yyy_xzz[j] + f2z * s_xyyy_xxzz[j];

                s_xyyy_xyyy[j] = fr * s_yyy_xyyy[j] + f2t * s_yyy_yyy[j];

                t_xyyy_xyyy[j] = fr * t_yyy_xyyy[j] + f2t * t_yyy_yyy[j] + f2z * s_xyyy_xyyy[j];

                s_xyyy_xyyz[j] = fr * s_yyy_xyyz[j] + f2t * s_yyy_yyz[j];

                t_xyyy_xyyz[j] = fr * t_yyy_xyyz[j] + f2t * t_yyy_yyz[j] + f2z * s_xyyy_xyyz[j];

                s_xyyy_xyzz[j] = fr * s_yyy_xyzz[j] + f2t * s_yyy_yzz[j];

                t_xyyy_xyzz[j] = fr * t_yyy_xyzz[j] + f2t * t_yyy_yzz[j] + f2z * s_xyyy_xyzz[j];

                s_xyyy_xzzz[j] = fr * s_yyy_xzzz[j] + f2t * s_yyy_zzz[j];

                t_xyyy_xzzz[j] = fr * t_yyy_xzzz[j] + f2t * t_yyy_zzz[j] + f2z * s_xyyy_xzzz[j];

                s_xyyy_yyyy[j] = fr * s_yyy_yyyy[j];

                t_xyyy_yyyy[j] = fr * t_yyy_yyyy[j] + f2z * s_xyyy_yyyy[j];

                s_xyyy_yyyz[j] = fr * s_yyy_yyyz[j];

                t_xyyy_yyyz[j] = fr * t_yyy_yyyz[j] + f2z * s_xyyy_yyyz[j];

                s_xyyy_yyzz[j] = fr * s_yyy_yyzz[j];

                t_xyyy_yyzz[j] = fr * t_yyy_yyzz[j] + f2z * s_xyyy_yyzz[j];

                s_xyyy_yzzz[j] = fr * s_yyy_yzzz[j];

                t_xyyy_yzzz[j] = fr * t_yyy_yzzz[j] + f2z * s_xyyy_yzzz[j];

                s_xyyy_zzzz[j] = fr * s_yyy_zzzz[j];

                t_xyyy_zzzz[j] = fr * t_yyy_zzzz[j] + f2z * s_xyyy_zzzz[j];

                s_xyyz_xxxx[j] = fr * s_yyz_xxxx[j] + f2t * 4.0 * s_yyz_xxx[j];

                t_xyyz_xxxx[j] = fr * t_yyz_xxxx[j] + f2t * 4.0 * t_yyz_xxx[j] + f2z * s_xyyz_xxxx[j];

                s_xyyz_xxxy[j] = fr * s_yyz_xxxy[j] + f2t * 3.0 * s_yyz_xxy[j];

                t_xyyz_xxxy[j] = fr * t_yyz_xxxy[j] + f2t * 3.0 * t_yyz_xxy[j] + f2z * s_xyyz_xxxy[j];

                s_xyyz_xxxz[j] = fr * s_yyz_xxxz[j] + f2t * 3.0 * s_yyz_xxz[j];

                t_xyyz_xxxz[j] = fr * t_yyz_xxxz[j] + f2t * 3.0 * t_yyz_xxz[j] + f2z * s_xyyz_xxxz[j];

                s_xyyz_xxyy[j] = fr * s_yyz_xxyy[j] + f2t * 2.0 * s_yyz_xyy[j];

                t_xyyz_xxyy[j] = fr * t_yyz_xxyy[j] + f2t * 2.0 * t_yyz_xyy[j] + f2z * s_xyyz_xxyy[j];

                s_xyyz_xxyz[j] = fr * s_yyz_xxyz[j] + f2t * 2.0 * s_yyz_xyz[j];

                t_xyyz_xxyz[j] = fr * t_yyz_xxyz[j] + f2t * 2.0 * t_yyz_xyz[j] + f2z * s_xyyz_xxyz[j];

                s_xyyz_xxzz[j] = fr * s_yyz_xxzz[j] + f2t * 2.0 * s_yyz_xzz[j];

                t_xyyz_xxzz[j] = fr * t_yyz_xxzz[j] + f2t * 2.0 * t_yyz_xzz[j] + f2z * s_xyyz_xxzz[j];

                s_xyyz_xyyy[j] = fr * s_yyz_xyyy[j] + f2t * s_yyz_yyy[j];

                t_xyyz_xyyy[j] = fr * t_yyz_xyyy[j] + f2t * t_yyz_yyy[j] + f2z * s_xyyz_xyyy[j];

                s_xyyz_xyyz[j] = fr * s_yyz_xyyz[j] + f2t * s_yyz_yyz[j];

                t_xyyz_xyyz[j] = fr * t_yyz_xyyz[j] + f2t * t_yyz_yyz[j] + f2z * s_xyyz_xyyz[j];

                s_xyyz_xyzz[j] = fr * s_yyz_xyzz[j] + f2t * s_yyz_yzz[j];

                t_xyyz_xyzz[j] = fr * t_yyz_xyzz[j] + f2t * t_yyz_yzz[j] + f2z * s_xyyz_xyzz[j];

                s_xyyz_xzzz[j] = fr * s_yyz_xzzz[j] + f2t * s_yyz_zzz[j];

                t_xyyz_xzzz[j] = fr * t_yyz_xzzz[j] + f2t * t_yyz_zzz[j] + f2z * s_xyyz_xzzz[j];

                s_xyyz_yyyy[j] = fr * s_yyz_yyyy[j];

                t_xyyz_yyyy[j] = fr * t_yyz_yyyy[j] + f2z * s_xyyz_yyyy[j];

                s_xyyz_yyyz[j] = fr * s_yyz_yyyz[j];

                t_xyyz_yyyz[j] = fr * t_yyz_yyyz[j] + f2z * s_xyyz_yyyz[j];

                s_xyyz_yyzz[j] = fr * s_yyz_yyzz[j];

                t_xyyz_yyzz[j] = fr * t_yyz_yyzz[j] + f2z * s_xyyz_yyzz[j];

                s_xyyz_yzzz[j] = fr * s_yyz_yzzz[j];

                t_xyyz_yzzz[j] = fr * t_yyz_yzzz[j] + f2z * s_xyyz_yzzz[j];

                s_xyyz_zzzz[j] = fr * s_yyz_zzzz[j];

                t_xyyz_zzzz[j] = fr * t_yyz_zzzz[j] + f2z * s_xyyz_zzzz[j];

                s_xyzz_xxxx[j] = fr * s_yzz_xxxx[j] + f2t * 4.0 * s_yzz_xxx[j];

                t_xyzz_xxxx[j] = fr * t_yzz_xxxx[j] + f2t * 4.0 * t_yzz_xxx[j] + f2z * s_xyzz_xxxx[j];

                s_xyzz_xxxy[j] = fr * s_yzz_xxxy[j] + f2t * 3.0 * s_yzz_xxy[j];

                t_xyzz_xxxy[j] = fr * t_yzz_xxxy[j] + f2t * 3.0 * t_yzz_xxy[j] + f2z * s_xyzz_xxxy[j];

                s_xyzz_xxxz[j] = fr * s_yzz_xxxz[j] + f2t * 3.0 * s_yzz_xxz[j];

                t_xyzz_xxxz[j] = fr * t_yzz_xxxz[j] + f2t * 3.0 * t_yzz_xxz[j] + f2z * s_xyzz_xxxz[j];

                s_xyzz_xxyy[j] = fr * s_yzz_xxyy[j] + f2t * 2.0 * s_yzz_xyy[j];

                t_xyzz_xxyy[j] = fr * t_yzz_xxyy[j] + f2t * 2.0 * t_yzz_xyy[j] + f2z * s_xyzz_xxyy[j];

                s_xyzz_xxyz[j] = fr * s_yzz_xxyz[j] + f2t * 2.0 * s_yzz_xyz[j];

                t_xyzz_xxyz[j] = fr * t_yzz_xxyz[j] + f2t * 2.0 * t_yzz_xyz[j] + f2z * s_xyzz_xxyz[j];

                s_xyzz_xxzz[j] = fr * s_yzz_xxzz[j] + f2t * 2.0 * s_yzz_xzz[j];

                t_xyzz_xxzz[j] = fr * t_yzz_xxzz[j] + f2t * 2.0 * t_yzz_xzz[j] + f2z * s_xyzz_xxzz[j];

                s_xyzz_xyyy[j] = fr * s_yzz_xyyy[j] + f2t * s_yzz_yyy[j];

                t_xyzz_xyyy[j] = fr * t_yzz_xyyy[j] + f2t * t_yzz_yyy[j] + f2z * s_xyzz_xyyy[j];

                s_xyzz_xyyz[j] = fr * s_yzz_xyyz[j] + f2t * s_yzz_yyz[j];

                t_xyzz_xyyz[j] = fr * t_yzz_xyyz[j] + f2t * t_yzz_yyz[j] + f2z * s_xyzz_xyyz[j];

                s_xyzz_xyzz[j] = fr * s_yzz_xyzz[j] + f2t * s_yzz_yzz[j];

                t_xyzz_xyzz[j] = fr * t_yzz_xyzz[j] + f2t * t_yzz_yzz[j] + f2z * s_xyzz_xyzz[j];

                s_xyzz_xzzz[j] = fr * s_yzz_xzzz[j] + f2t * s_yzz_zzz[j];

                t_xyzz_xzzz[j] = fr * t_yzz_xzzz[j] + f2t * t_yzz_zzz[j] + f2z * s_xyzz_xzzz[j];

                s_xyzz_yyyy[j] = fr * s_yzz_yyyy[j];

                t_xyzz_yyyy[j] = fr * t_yzz_yyyy[j] + f2z * s_xyzz_yyyy[j];

                s_xyzz_yyyz[j] = fr * s_yzz_yyyz[j];

                t_xyzz_yyyz[j] = fr * t_yzz_yyyz[j] + f2z * s_xyzz_yyyz[j];

                s_xyzz_yyzz[j] = fr * s_yzz_yyzz[j];

                t_xyzz_yyzz[j] = fr * t_yzz_yyzz[j] + f2z * s_xyzz_yyzz[j];

                s_xyzz_yzzz[j] = fr * s_yzz_yzzz[j];

                t_xyzz_yzzz[j] = fr * t_yzz_yzzz[j] + f2z * s_xyzz_yzzz[j];

                s_xyzz_zzzz[j] = fr * s_yzz_zzzz[j];

                t_xyzz_zzzz[j] = fr * t_yzz_zzzz[j] + f2z * s_xyzz_zzzz[j];

                s_xzzz_xxxx[j] = fr * s_zzz_xxxx[j] + f2t * 4.0 * s_zzz_xxx[j];

                t_xzzz_xxxx[j] = fr * t_zzz_xxxx[j] + f2t * 4.0 * t_zzz_xxx[j] + f2z * s_xzzz_xxxx[j];

                s_xzzz_xxxy[j] = fr * s_zzz_xxxy[j] + f2t * 3.0 * s_zzz_xxy[j];

                t_xzzz_xxxy[j] = fr * t_zzz_xxxy[j] + f2t * 3.0 * t_zzz_xxy[j] + f2z * s_xzzz_xxxy[j];

                s_xzzz_xxxz[j] = fr * s_zzz_xxxz[j] + f2t * 3.0 * s_zzz_xxz[j];

                t_xzzz_xxxz[j] = fr * t_zzz_xxxz[j] + f2t * 3.0 * t_zzz_xxz[j] + f2z * s_xzzz_xxxz[j];

                s_xzzz_xxyy[j] = fr * s_zzz_xxyy[j] + f2t * 2.0 * s_zzz_xyy[j];

                t_xzzz_xxyy[j] = fr * t_zzz_xxyy[j] + f2t * 2.0 * t_zzz_xyy[j] + f2z * s_xzzz_xxyy[j];

                s_xzzz_xxyz[j] = fr * s_zzz_xxyz[j] + f2t * 2.0 * s_zzz_xyz[j];

                t_xzzz_xxyz[j] = fr * t_zzz_xxyz[j] + f2t * 2.0 * t_zzz_xyz[j] + f2z * s_xzzz_xxyz[j];

                s_xzzz_xxzz[j] = fr * s_zzz_xxzz[j] + f2t * 2.0 * s_zzz_xzz[j];

                t_xzzz_xxzz[j] = fr * t_zzz_xxzz[j] + f2t * 2.0 * t_zzz_xzz[j] + f2z * s_xzzz_xxzz[j];

                s_xzzz_xyyy[j] = fr * s_zzz_xyyy[j] + f2t * s_zzz_yyy[j];

                t_xzzz_xyyy[j] = fr * t_zzz_xyyy[j] + f2t * t_zzz_yyy[j] + f2z * s_xzzz_xyyy[j];

                s_xzzz_xyyz[j] = fr * s_zzz_xyyz[j] + f2t * s_zzz_yyz[j];

                t_xzzz_xyyz[j] = fr * t_zzz_xyyz[j] + f2t * t_zzz_yyz[j] + f2z * s_xzzz_xyyz[j];

                s_xzzz_xyzz[j] = fr * s_zzz_xyzz[j] + f2t * s_zzz_yzz[j];

                t_xzzz_xyzz[j] = fr * t_zzz_xyzz[j] + f2t * t_zzz_yzz[j] + f2z * s_xzzz_xyzz[j];

                s_xzzz_xzzz[j] = fr * s_zzz_xzzz[j] + f2t * s_zzz_zzz[j];

                t_xzzz_xzzz[j] = fr * t_zzz_xzzz[j] + f2t * t_zzz_zzz[j] + f2z * s_xzzz_xzzz[j];

                s_xzzz_yyyy[j] = fr * s_zzz_yyyy[j];

                t_xzzz_yyyy[j] = fr * t_zzz_yyyy[j] + f2z * s_xzzz_yyyy[j];

                s_xzzz_yyyz[j] = fr * s_zzz_yyyz[j];

                t_xzzz_yyyz[j] = fr * t_zzz_yyyz[j] + f2z * s_xzzz_yyyz[j];

                s_xzzz_yyzz[j] = fr * s_zzz_yyzz[j];

                t_xzzz_yyzz[j] = fr * t_zzz_yyzz[j] + f2z * s_xzzz_yyzz[j];

                s_xzzz_yzzz[j] = fr * s_zzz_yzzz[j];

                t_xzzz_yzzz[j] = fr * t_zzz_yzzz[j] + f2z * s_xzzz_yzzz[j];

                s_xzzz_zzzz[j] = fr * s_zzz_zzzz[j];

                t_xzzz_zzzz[j] = fr * t_zzz_zzzz[j] + f2z * s_xzzz_zzzz[j];

                // leading y component

                fr = pay[j];

                s_yyyy_xxxx[j] = fr * s_yyy_xxxx[j] + f2t * 3.0 * s_yy_xxxx[j];

                t_yyyy_xxxx[j] = fr * t_yyy_xxxx[j] + f2t * 3.0 * t_yy_xxxx[j] + f2z * (s_yyyy_xxxx[j] - f2a * 3.0 * s_yy_xxxx[j]);

                s_yyyy_xxxy[j] = fr * s_yyy_xxxy[j] + f2t * (3.0 * s_yy_xxxy[j] + s_yyy_xxx[j]);

                t_yyyy_xxxy[j] = fr * t_yyy_xxxy[j] + f2t * (3.0 * t_yy_xxxy[j] + t_yyy_xxx[j]) + f2z * (s_yyyy_xxxy[j] - f2a * 3.0 * s_yy_xxxy[j]);

                s_yyyy_xxxz[j] = fr * s_yyy_xxxz[j] + f2t * 3.0 * s_yy_xxxz[j];

                t_yyyy_xxxz[j] = fr * t_yyy_xxxz[j] + f2t * 3.0 * t_yy_xxxz[j] + f2z * (s_yyyy_xxxz[j] - f2a * 3.0 * s_yy_xxxz[j]);

                s_yyyy_xxyy[j] = fr * s_yyy_xxyy[j] + f2t * (3.0 * s_yy_xxyy[j] + 2.0 * s_yyy_xxy[j]);

                t_yyyy_xxyy[j] = fr * t_yyy_xxyy[j] + f2t * (3.0 * t_yy_xxyy[j] + 2.0 * t_yyy_xxy[j]) + f2z * (s_yyyy_xxyy[j] - f2a * 3.0 * s_yy_xxyy[j]);

                s_yyyy_xxyz[j] = fr * s_yyy_xxyz[j] + f2t * (3.0 * s_yy_xxyz[j] + s_yyy_xxz[j]);

                t_yyyy_xxyz[j] = fr * t_yyy_xxyz[j] + f2t * (3.0 * t_yy_xxyz[j] + t_yyy_xxz[j]) + f2z * (s_yyyy_xxyz[j] - f2a * 3.0 * s_yy_xxyz[j]);

                s_yyyy_xxzz[j] = fr * s_yyy_xxzz[j] + f2t * 3.0 * s_yy_xxzz[j];

                t_yyyy_xxzz[j] = fr * t_yyy_xxzz[j] + f2t * 3.0 * t_yy_xxzz[j] + f2z * (s_yyyy_xxzz[j] - f2a * 3.0 * s_yy_xxzz[j]);

                s_yyyy_xyyy[j] = fr * s_yyy_xyyy[j] + f2t * (3.0 * s_yy_xyyy[j] + 3.0 * s_yyy_xyy[j]);

                t_yyyy_xyyy[j] = fr * t_yyy_xyyy[j] + f2t * (3.0 * t_yy_xyyy[j] + 3.0 * t_yyy_xyy[j]) + f2z * (s_yyyy_xyyy[j] - f2a * 3.0 * s_yy_xyyy[j]);

                s_yyyy_xyyz[j] = fr * s_yyy_xyyz[j] + f2t * (3.0 * s_yy_xyyz[j] + 2.0 * s_yyy_xyz[j]);

                t_yyyy_xyyz[j] = fr * t_yyy_xyyz[j] + f2t * (3.0 * t_yy_xyyz[j] + 2.0 * t_yyy_xyz[j]) + f2z * (s_yyyy_xyyz[j] - f2a * 3.0 * s_yy_xyyz[j]);

                s_yyyy_xyzz[j] = fr * s_yyy_xyzz[j] + f2t * (3.0 * s_yy_xyzz[j] + s_yyy_xzz[j]);

                t_yyyy_xyzz[j] = fr * t_yyy_xyzz[j] + f2t * (3.0 * t_yy_xyzz[j] + t_yyy_xzz[j]) + f2z * (s_yyyy_xyzz[j] - f2a * 3.0 * s_yy_xyzz[j]);

                s_yyyy_xzzz[j] = fr * s_yyy_xzzz[j] + f2t * 3.0 * s_yy_xzzz[j];

                t_yyyy_xzzz[j] = fr * t_yyy_xzzz[j] + f2t * 3.0 * t_yy_xzzz[j] + f2z * (s_yyyy_xzzz[j] - f2a * 3.0 * s_yy_xzzz[j]);

                s_yyyy_yyyy[j] = fr * s_yyy_yyyy[j] + f2t * (3.0 * s_yy_yyyy[j] + 4.0 * s_yyy_yyy[j]);

                t_yyyy_yyyy[j] = fr * t_yyy_yyyy[j] + f2t * (3.0 * t_yy_yyyy[j] + 4.0 * t_yyy_yyy[j]) + f2z * (s_yyyy_yyyy[j] - f2a * 3.0 * s_yy_yyyy[j]);

                s_yyyy_yyyz[j] = fr * s_yyy_yyyz[j] + f2t * (3.0 * s_yy_yyyz[j] + 3.0 * s_yyy_yyz[j]);

                t_yyyy_yyyz[j] = fr * t_yyy_yyyz[j] + f2t * (3.0 * t_yy_yyyz[j] + 3.0 * t_yyy_yyz[j]) + f2z * (s_yyyy_yyyz[j] - f2a * 3.0 * s_yy_yyyz[j]);

                s_yyyy_yyzz[j] = fr * s_yyy_yyzz[j] + f2t * (3.0 * s_yy_yyzz[j] + 2.0 * s_yyy_yzz[j]);

                t_yyyy_yyzz[j] = fr * t_yyy_yyzz[j] + f2t * (3.0 * t_yy_yyzz[j] + 2.0 * t_yyy_yzz[j]) + f2z * (s_yyyy_yyzz[j] - f2a * 3.0 * s_yy_yyzz[j]);

                s_yyyy_yzzz[j] = fr * s_yyy_yzzz[j] + f2t * (3.0 * s_yy_yzzz[j] + s_yyy_zzz[j]);

                t_yyyy_yzzz[j] = fr * t_yyy_yzzz[j] + f2t * (3.0 * t_yy_yzzz[j] + t_yyy_zzz[j]) + f2z * (s_yyyy_yzzz[j] - f2a * 3.0 * s_yy_yzzz[j]);

                s_yyyy_zzzz[j] = fr * s_yyy_zzzz[j] + f2t * 3.0 * s_yy_zzzz[j];

                t_yyyy_zzzz[j] = fr * t_yyy_zzzz[j] + f2t * 3.0 * t_yy_zzzz[j] + f2z * (s_yyyy_zzzz[j] - f2a * 3.0 * s_yy_zzzz[j]);

                s_yyyz_xxxx[j] = fr * s_yyz_xxxx[j] + f2t * 2.0 * s_yz_xxxx[j];

                t_yyyz_xxxx[j] = fr * t_yyz_xxxx[j] + f2t * 2.0 * t_yz_xxxx[j] + f2z * (s_yyyz_xxxx[j] - f2a * 2.0 * s_yz_xxxx[j]);

                s_yyyz_xxxy[j] = fr * s_yyz_xxxy[j] + f2t * (2.0 * s_yz_xxxy[j] + s_yyz_xxx[j]);

                t_yyyz_xxxy[j] = fr * t_yyz_xxxy[j] + f2t * (2.0 * t_yz_xxxy[j] + t_yyz_xxx[j]) + f2z * (s_yyyz_xxxy[j] - f2a * 2.0 * s_yz_xxxy[j]);

                s_yyyz_xxxz[j] = fr * s_yyz_xxxz[j] + f2t * 2.0 * s_yz_xxxz[j];

                t_yyyz_xxxz[j] = fr * t_yyz_xxxz[j] + f2t * 2.0 * t_yz_xxxz[j] + f2z * (s_yyyz_xxxz[j] - f2a * 2.0 * s_yz_xxxz[j]);

                s_yyyz_xxyy[j] = fr * s_yyz_xxyy[j] + f2t * (2.0 * s_yz_xxyy[j] + 2.0 * s_yyz_xxy[j]);

                t_yyyz_xxyy[j] = fr * t_yyz_xxyy[j] + f2t * (2.0 * t_yz_xxyy[j] + 2.0 * t_yyz_xxy[j]) + f2z * (s_yyyz_xxyy[j] - f2a * 2.0 * s_yz_xxyy[j]);

                s_yyyz_xxyz[j] = fr * s_yyz_xxyz[j] + f2t * (2.0 * s_yz_xxyz[j] + s_yyz_xxz[j]);

                t_yyyz_xxyz[j] = fr * t_yyz_xxyz[j] + f2t * (2.0 * t_yz_xxyz[j] + t_yyz_xxz[j]) + f2z * (s_yyyz_xxyz[j] - f2a * 2.0 * s_yz_xxyz[j]);

                s_yyyz_xxzz[j] = fr * s_yyz_xxzz[j] + f2t * 2.0 * s_yz_xxzz[j];

                t_yyyz_xxzz[j] = fr * t_yyz_xxzz[j] + f2t * 2.0 * t_yz_xxzz[j] + f2z * (s_yyyz_xxzz[j] - f2a * 2.0 * s_yz_xxzz[j]);

                s_yyyz_xyyy[j] = fr * s_yyz_xyyy[j] + f2t * (2.0 * s_yz_xyyy[j] + 3.0 * s_yyz_xyy[j]);

                t_yyyz_xyyy[j] = fr * t_yyz_xyyy[j] + f2t * (2.0 * t_yz_xyyy[j] + 3.0 * t_yyz_xyy[j]) + f2z * (s_yyyz_xyyy[j] - f2a * 2.0 * s_yz_xyyy[j]);

                s_yyyz_xyyz[j] = fr * s_yyz_xyyz[j] + f2t * (2.0 * s_yz_xyyz[j] + 2.0 * s_yyz_xyz[j]);

                t_yyyz_xyyz[j] = fr * t_yyz_xyyz[j] + f2t * (2.0 * t_yz_xyyz[j] + 2.0 * t_yyz_xyz[j]) + f2z * (s_yyyz_xyyz[j] - f2a * 2.0 * s_yz_xyyz[j]);

                s_yyyz_xyzz[j] = fr * s_yyz_xyzz[j] + f2t * (2.0 * s_yz_xyzz[j] + s_yyz_xzz[j]);

                t_yyyz_xyzz[j] = fr * t_yyz_xyzz[j] + f2t * (2.0 * t_yz_xyzz[j] + t_yyz_xzz[j]) + f2z * (s_yyyz_xyzz[j] - f2a * 2.0 * s_yz_xyzz[j]);

                s_yyyz_xzzz[j] = fr * s_yyz_xzzz[j] + f2t * 2.0 * s_yz_xzzz[j];

                t_yyyz_xzzz[j] = fr * t_yyz_xzzz[j] + f2t * 2.0 * t_yz_xzzz[j] + f2z * (s_yyyz_xzzz[j] - f2a * 2.0 * s_yz_xzzz[j]);

                s_yyyz_yyyy[j] = fr * s_yyz_yyyy[j] + f2t * (2.0 * s_yz_yyyy[j] + 4.0 * s_yyz_yyy[j]);

                t_yyyz_yyyy[j] = fr * t_yyz_yyyy[j] + f2t * (2.0 * t_yz_yyyy[j] + 4.0 * t_yyz_yyy[j]) + f2z * (s_yyyz_yyyy[j] - f2a * 2.0 * s_yz_yyyy[j]);

                s_yyyz_yyyz[j] = fr * s_yyz_yyyz[j] + f2t * (2.0 * s_yz_yyyz[j] + 3.0 * s_yyz_yyz[j]);

                t_yyyz_yyyz[j] = fr * t_yyz_yyyz[j] + f2t * (2.0 * t_yz_yyyz[j] + 3.0 * t_yyz_yyz[j]) + f2z * (s_yyyz_yyyz[j] - f2a * 2.0 * s_yz_yyyz[j]);

                s_yyyz_yyzz[j] = fr * s_yyz_yyzz[j] + f2t * (2.0 * s_yz_yyzz[j] + 2.0 * s_yyz_yzz[j]);

                t_yyyz_yyzz[j] = fr * t_yyz_yyzz[j] + f2t * (2.0 * t_yz_yyzz[j] + 2.0 * t_yyz_yzz[j]) + f2z * (s_yyyz_yyzz[j] - f2a * 2.0 * s_yz_yyzz[j]);

                s_yyyz_yzzz[j] = fr * s_yyz_yzzz[j] + f2t * (2.0 * s_yz_yzzz[j] + s_yyz_zzz[j]);

                t_yyyz_yzzz[j] = fr * t_yyz_yzzz[j] + f2t * (2.0 * t_yz_yzzz[j] + t_yyz_zzz[j]) + f2z * (s_yyyz_yzzz[j] - f2a * 2.0 * s_yz_yzzz[j]);

                s_yyyz_zzzz[j] = fr * s_yyz_zzzz[j] + f2t * 2.0 * s_yz_zzzz[j];

                t_yyyz_zzzz[j] = fr * t_yyz_zzzz[j] + f2t * 2.0 * t_yz_zzzz[j] + f2z * (s_yyyz_zzzz[j] - f2a * 2.0 * s_yz_zzzz[j]);

                s_yyzz_xxxx[j] = fr * s_yzz_xxxx[j] + f2t * s_zz_xxxx[j];

                t_yyzz_xxxx[j] = fr * t_yzz_xxxx[j] + f2t * t_zz_xxxx[j] + f2z * (s_yyzz_xxxx[j] - f2a * s_zz_xxxx[j]);

                s_yyzz_xxxy[j] = fr * s_yzz_xxxy[j] + f2t * (s_zz_xxxy[j] + s_yzz_xxx[j]);

                t_yyzz_xxxy[j] = fr * t_yzz_xxxy[j] + f2t * (t_zz_xxxy[j] + t_yzz_xxx[j]) + f2z * (s_yyzz_xxxy[j] - f2a * s_zz_xxxy[j]);

                s_yyzz_xxxz[j] = fr * s_yzz_xxxz[j] + f2t * s_zz_xxxz[j];

                t_yyzz_xxxz[j] = fr * t_yzz_xxxz[j] + f2t * t_zz_xxxz[j] + f2z * (s_yyzz_xxxz[j] - f2a * s_zz_xxxz[j]);

                s_yyzz_xxyy[j] = fr * s_yzz_xxyy[j] + f2t * (s_zz_xxyy[j] + 2.0 * s_yzz_xxy[j]);

                t_yyzz_xxyy[j] = fr * t_yzz_xxyy[j] + f2t * (t_zz_xxyy[j] + 2.0 * t_yzz_xxy[j]) + f2z * (s_yyzz_xxyy[j] - f2a * s_zz_xxyy[j]);

                s_yyzz_xxyz[j] = fr * s_yzz_xxyz[j] + f2t * (s_zz_xxyz[j] + s_yzz_xxz[j]);

                t_yyzz_xxyz[j] = fr * t_yzz_xxyz[j] + f2t * (t_zz_xxyz[j] + t_yzz_xxz[j]) + f2z * (s_yyzz_xxyz[j] - f2a * s_zz_xxyz[j]);

                s_yyzz_xxzz[j] = fr * s_yzz_xxzz[j] + f2t * s_zz_xxzz[j];

                t_yyzz_xxzz[j] = fr * t_yzz_xxzz[j] + f2t * t_zz_xxzz[j] + f2z * (s_yyzz_xxzz[j] - f2a * s_zz_xxzz[j]);

                s_yyzz_xyyy[j] = fr * s_yzz_xyyy[j] + f2t * (s_zz_xyyy[j] + 3.0 * s_yzz_xyy[j]);

                t_yyzz_xyyy[j] = fr * t_yzz_xyyy[j] + f2t * (t_zz_xyyy[j] + 3.0 * t_yzz_xyy[j]) + f2z * (s_yyzz_xyyy[j] - f2a * s_zz_xyyy[j]);

                s_yyzz_xyyz[j] = fr * s_yzz_xyyz[j] + f2t * (s_zz_xyyz[j] + 2.0 * s_yzz_xyz[j]);

                t_yyzz_xyyz[j] = fr * t_yzz_xyyz[j] + f2t * (t_zz_xyyz[j] + 2.0 * t_yzz_xyz[j]) + f2z * (s_yyzz_xyyz[j] - f2a * s_zz_xyyz[j]);

                s_yyzz_xyzz[j] = fr * s_yzz_xyzz[j] + f2t * (s_zz_xyzz[j] + s_yzz_xzz[j]);

                t_yyzz_xyzz[j] = fr * t_yzz_xyzz[j] + f2t * (t_zz_xyzz[j] + t_yzz_xzz[j]) + f2z * (s_yyzz_xyzz[j] - f2a * s_zz_xyzz[j]);

                s_yyzz_xzzz[j] = fr * s_yzz_xzzz[j] + f2t * s_zz_xzzz[j];

                t_yyzz_xzzz[j] = fr * t_yzz_xzzz[j] + f2t * t_zz_xzzz[j] + f2z * (s_yyzz_xzzz[j] - f2a * s_zz_xzzz[j]);

                s_yyzz_yyyy[j] = fr * s_yzz_yyyy[j] + f2t * (s_zz_yyyy[j] + 4.0 * s_yzz_yyy[j]);

                t_yyzz_yyyy[j] = fr * t_yzz_yyyy[j] + f2t * (t_zz_yyyy[j] + 4.0 * t_yzz_yyy[j]) + f2z * (s_yyzz_yyyy[j] - f2a * s_zz_yyyy[j]);

                s_yyzz_yyyz[j] = fr * s_yzz_yyyz[j] + f2t * (s_zz_yyyz[j] + 3.0 * s_yzz_yyz[j]);

                t_yyzz_yyyz[j] = fr * t_yzz_yyyz[j] + f2t * (t_zz_yyyz[j] + 3.0 * t_yzz_yyz[j]) + f2z * (s_yyzz_yyyz[j] - f2a * s_zz_yyyz[j]);

                s_yyzz_yyzz[j] = fr * s_yzz_yyzz[j] + f2t * (s_zz_yyzz[j] + 2.0 * s_yzz_yzz[j]);

                t_yyzz_yyzz[j] = fr * t_yzz_yyzz[j] + f2t * (t_zz_yyzz[j] + 2.0 * t_yzz_yzz[j]) + f2z * (s_yyzz_yyzz[j] - f2a * s_zz_yyzz[j]);

                s_yyzz_yzzz[j] = fr * s_yzz_yzzz[j] + f2t * (s_zz_yzzz[j] + s_yzz_zzz[j]);

                t_yyzz_yzzz[j] = fr * t_yzz_yzzz[j] + f2t * (t_zz_yzzz[j] + t_yzz_zzz[j]) + f2z * (s_yyzz_yzzz[j] - f2a * s_zz_yzzz[j]);

                s_yyzz_zzzz[j] = fr * s_yzz_zzzz[j] + f2t * s_zz_zzzz[j];

                t_yyzz_zzzz[j] = fr * t_yzz_zzzz[j] + f2t * t_zz_zzzz[j] + f2z * (s_yyzz_zzzz[j] - f2a * s_zz_zzzz[j]);

                s_yzzz_xxxx[j] = fr * s_zzz_xxxx[j];

                t_yzzz_xxxx[j] = fr * t_zzz_xxxx[j] + f2z * s_yzzz_xxxx[j];

                s_yzzz_xxxy[j] = fr * s_zzz_xxxy[j] + f2t * s_zzz_xxx[j];

                t_yzzz_xxxy[j] = fr * t_zzz_xxxy[j] + f2t * t_zzz_xxx[j] + f2z * s_yzzz_xxxy[j];

                s_yzzz_xxxz[j] = fr * s_zzz_xxxz[j];

                t_yzzz_xxxz[j] = fr * t_zzz_xxxz[j] + f2z * s_yzzz_xxxz[j];

                s_yzzz_xxyy[j] = fr * s_zzz_xxyy[j] + f2t * 2.0 * s_zzz_xxy[j];

                t_yzzz_xxyy[j] = fr * t_zzz_xxyy[j] + f2t * 2.0 * t_zzz_xxy[j] + f2z * s_yzzz_xxyy[j];

                s_yzzz_xxyz[j] = fr * s_zzz_xxyz[j] + f2t * s_zzz_xxz[j];

                t_yzzz_xxyz[j] = fr * t_zzz_xxyz[j] + f2t * t_zzz_xxz[j] + f2z * s_yzzz_xxyz[j];

                s_yzzz_xxzz[j] = fr * s_zzz_xxzz[j];

                t_yzzz_xxzz[j] = fr * t_zzz_xxzz[j] + f2z * s_yzzz_xxzz[j];

                s_yzzz_xyyy[j] = fr * s_zzz_xyyy[j] + f2t * 3.0 * s_zzz_xyy[j];

                t_yzzz_xyyy[j] = fr * t_zzz_xyyy[j] + f2t * 3.0 * t_zzz_xyy[j] + f2z * s_yzzz_xyyy[j];

                s_yzzz_xyyz[j] = fr * s_zzz_xyyz[j] + f2t * 2.0 * s_zzz_xyz[j];

                t_yzzz_xyyz[j] = fr * t_zzz_xyyz[j] + f2t * 2.0 * t_zzz_xyz[j] + f2z * s_yzzz_xyyz[j];

                s_yzzz_xyzz[j] = fr * s_zzz_xyzz[j] + f2t * s_zzz_xzz[j];

                t_yzzz_xyzz[j] = fr * t_zzz_xyzz[j] + f2t * t_zzz_xzz[j] + f2z * s_yzzz_xyzz[j];

                s_yzzz_xzzz[j] = fr * s_zzz_xzzz[j];

                t_yzzz_xzzz[j] = fr * t_zzz_xzzz[j] + f2z * s_yzzz_xzzz[j];

                s_yzzz_yyyy[j] = fr * s_zzz_yyyy[j] + f2t * 4.0 * s_zzz_yyy[j];

                t_yzzz_yyyy[j] = fr * t_zzz_yyyy[j] + f2t * 4.0 * t_zzz_yyy[j] + f2z * s_yzzz_yyyy[j];

                s_yzzz_yyyz[j] = fr * s_zzz_yyyz[j] + f2t * 3.0 * s_zzz_yyz[j];

                t_yzzz_yyyz[j] = fr * t_zzz_yyyz[j] + f2t * 3.0 * t_zzz_yyz[j] + f2z * s_yzzz_yyyz[j];

                s_yzzz_yyzz[j] = fr * s_zzz_yyzz[j] + f2t * 2.0 * s_zzz_yzz[j];

                t_yzzz_yyzz[j] = fr * t_zzz_yyzz[j] + f2t * 2.0 * t_zzz_yzz[j] + f2z * s_yzzz_yyzz[j];

                s_yzzz_yzzz[j] = fr * s_zzz_yzzz[j] + f2t * s_zzz_zzz[j];

                t_yzzz_yzzz[j] = fr * t_zzz_yzzz[j] + f2t * t_zzz_zzz[j] + f2z * s_yzzz_yzzz[j];

                s_yzzz_zzzz[j] = fr * s_zzz_zzzz[j];

                t_yzzz_zzzz[j] = fr * t_zzz_zzzz[j] + f2z * s_yzzz_zzzz[j];

                // leading z component

                fr = paz[j];

                s_zzzz_xxxx[j] = fr * s_zzz_xxxx[j] + f2t * 3.0 * s_zz_xxxx[j];

                t_zzzz_xxxx[j] = fr * t_zzz_xxxx[j] + f2t * 3.0 * t_zz_xxxx[j] + f2z * (s_zzzz_xxxx[j] - f2a * 3.0 * s_zz_xxxx[j]);

                s_zzzz_xxxy[j] = fr * s_zzz_xxxy[j] + f2t * 3.0 * s_zz_xxxy[j];

                t_zzzz_xxxy[j] = fr * t_zzz_xxxy[j] + f2t * 3.0 * t_zz_xxxy[j] + f2z * (s_zzzz_xxxy[j] - f2a * 3.0 * s_zz_xxxy[j]);

                s_zzzz_xxxz[j] = fr * s_zzz_xxxz[j] + f2t * (3.0 * s_zz_xxxz[j] + s_zzz_xxx[j]);

                t_zzzz_xxxz[j] = fr * t_zzz_xxxz[j] + f2t * (3.0 * t_zz_xxxz[j] + t_zzz_xxx[j]) + f2z * (s_zzzz_xxxz[j] - f2a * 3.0 * s_zz_xxxz[j]);

                s_zzzz_xxyy[j] = fr * s_zzz_xxyy[j] + f2t * 3.0 * s_zz_xxyy[j];

                t_zzzz_xxyy[j] = fr * t_zzz_xxyy[j] + f2t * 3.0 * t_zz_xxyy[j] + f2z * (s_zzzz_xxyy[j] - f2a * 3.0 * s_zz_xxyy[j]);

                s_zzzz_xxyz[j] = fr * s_zzz_xxyz[j] + f2t * (3.0 * s_zz_xxyz[j] + s_zzz_xxy[j]);

                t_zzzz_xxyz[j] = fr * t_zzz_xxyz[j] + f2t * (3.0 * t_zz_xxyz[j] + t_zzz_xxy[j]) + f2z * (s_zzzz_xxyz[j] - f2a * 3.0 * s_zz_xxyz[j]);

                s_zzzz_xxzz[j] = fr * s_zzz_xxzz[j] + f2t * (3.0 * s_zz_xxzz[j] + 2.0 * s_zzz_xxz[j]);

                t_zzzz_xxzz[j] = fr * t_zzz_xxzz[j] + f2t * (3.0 * t_zz_xxzz[j] + 2.0 * t_zzz_xxz[j]) + f2z * (s_zzzz_xxzz[j] - f2a * 3.0 * s_zz_xxzz[j]);

                s_zzzz_xyyy[j] = fr * s_zzz_xyyy[j] + f2t * 3.0 * s_zz_xyyy[j];

                t_zzzz_xyyy[j] = fr * t_zzz_xyyy[j] + f2t * 3.0 * t_zz_xyyy[j] + f2z * (s_zzzz_xyyy[j] - f2a * 3.0 * s_zz_xyyy[j]);

                s_zzzz_xyyz[j] = fr * s_zzz_xyyz[j] + f2t * (3.0 * s_zz_xyyz[j] + s_zzz_xyy[j]);

                t_zzzz_xyyz[j] = fr * t_zzz_xyyz[j] + f2t * (3.0 * t_zz_xyyz[j] + t_zzz_xyy[j]) + f2z * (s_zzzz_xyyz[j] - f2a * 3.0 * s_zz_xyyz[j]);

                s_zzzz_xyzz[j] = fr * s_zzz_xyzz[j] + f2t * (3.0 * s_zz_xyzz[j] + 2.0 * s_zzz_xyz[j]);

                t_zzzz_xyzz[j] = fr * t_zzz_xyzz[j] + f2t * (3.0 * t_zz_xyzz[j] + 2.0 * t_zzz_xyz[j]) + f2z * (s_zzzz_xyzz[j] - f2a * 3.0 * s_zz_xyzz[j]);

                s_zzzz_xzzz[j] = fr * s_zzz_xzzz[j] + f2t * (3.0 * s_zz_xzzz[j] + 3.0 * s_zzz_xzz[j]);

                t_zzzz_xzzz[j] = fr * t_zzz_xzzz[j] + f2t * (3.0 * t_zz_xzzz[j] + 3.0 * t_zzz_xzz[j]) + f2z * (s_zzzz_xzzz[j] - f2a * 3.0 * s_zz_xzzz[j]);

                s_zzzz_yyyy[j] = fr * s_zzz_yyyy[j] + f2t * 3.0 * s_zz_yyyy[j];

                t_zzzz_yyyy[j] = fr * t_zzz_yyyy[j] + f2t * 3.0 * t_zz_yyyy[j] + f2z * (s_zzzz_yyyy[j] - f2a * 3.0 * s_zz_yyyy[j]);

                s_zzzz_yyyz[j] = fr * s_zzz_yyyz[j] + f2t * (3.0 * s_zz_yyyz[j] + s_zzz_yyy[j]);

                t_zzzz_yyyz[j] = fr * t_zzz_yyyz[j] + f2t * (3.0 * t_zz_yyyz[j] + t_zzz_yyy[j]) + f2z * (s_zzzz_yyyz[j] - f2a * 3.0 * s_zz_yyyz[j]);

                s_zzzz_yyzz[j] = fr * s_zzz_yyzz[j] + f2t * (3.0 * s_zz_yyzz[j] + 2.0 * s_zzz_yyz[j]);

                t_zzzz_yyzz[j] = fr * t_zzz_yyzz[j] + f2t * (3.0 * t_zz_yyzz[j] + 2.0 * t_zzz_yyz[j]) + f2z * (s_zzzz_yyzz[j] - f2a * 3.0 * s_zz_yyzz[j]);

                s_zzzz_yzzz[j] = fr * s_zzz_yzzz[j] + f2t * (3.0 * s_zz_yzzz[j] + 3.0 * s_zzz_yzz[j]);

                t_zzzz_yzzz[j] = fr * t_zzz_yzzz[j] + f2t * (3.0 * t_zz_yzzz[j] + 3.0 * t_zzz_yzz[j]) + f2z * (s_zzzz_yzzz[j] - f2a * 3.0 * s_zz_yzzz[j]);

                s_zzzz_zzzz[j] = fr * s_zzz_zzzz[j] + f2t * (3.0 * s_zz_zzzz[j] + 4.0 * s_zzz_zzz[j]);

                t_zzzz_zzzz[j] = fr * t_zzz_zzzz[j] + f2t * (3.0 * t_zz_zzzz[j] + 4.0 * t_zzz_zzz[j]) + f2z * (s_zzzz_zzzz[j] - f2a * 3.0 * s_zz_zzzz[j]);
            }

            idx++;
        }
    }
    
} // kinrecfunc namespace
