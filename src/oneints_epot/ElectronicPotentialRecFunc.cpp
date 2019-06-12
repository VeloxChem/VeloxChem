//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronicPotentialRecFunc.hpp"

#include <cmath>

#include "GenFunc.hpp"
#include "MathConst.hpp"

namespace epotrecfunc {  // epotrecfunc namespace

void
compElectronicPotentialForSS(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CBoysFunction&        bfTable,
                             CMemBlock<double>&          bfArguments,
                             CMemBlock2D<double>&        bfValues,
                             const int32_t               bfOrder,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  abDistances,
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

    // set up up pi prefactor

    auto fpi = 2.0 * std::pow(mathconst::getPiValue(), 2.50);

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up Obara-Saika prefactors

        auto fz = osFactors.data(6 * idx + 1);

        // compute Boys function argument

        auto fargs = bfArguments.data();

        #pragma omp simd aligned(fargs, fz, abx, aby, abz: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            fargs[j] = fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +

                                abz[j] * abz[j]);
        }

        // evaluate Boys function values

        bfTable.compute(bfValues, bfArguments, bfOrder);

        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(6 * idx);

        auto fia = osFactors.data(6 * idx + 2);

        auto fib = osFactors.data(6 * idx + 4);

        auto fb = bnorm[i];

        // compute overlap scaling factor

        #pragma omp simd aligned(fx, fia, fib, knorm, fargs: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            fargs[j] = fb * knorm[j] * fpi * std::sqrt(fx[j]) * fia[j] * fib[j];
        }

        // distribute (S|g(r,r')|S) integrals

        for (int32_t j = 0; j <= bfOrder; j++)
        {
            auto pidx = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, j});

            if (pidx != -1)
            {
                auto g_0_0 = primBuffer.data(pidx + idx);

                auto bvals = bfValues.data(j);

                #pragma omp simd aligned(g_0_0, bvals, fargs: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    g_0_0[k] = bvals[k] * fargs[k];
                }
            }
        }

        idx++;
    }
}

void
compElectronicPotentialForSP(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  pbDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 0, 1);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        // get position of integrals in primitves buffer

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|S)^(m+1) integrals

            auto g_0_0 = primBuffer.data(g1off + idx);

            // set up pointers to (S|g(r,r')|P)^(m) integrals

            auto g_0_x = primBuffer.data(goff + 3 * idx);

            auto g_0_y = primBuffer.data(goff + 3 * idx + 1);

            auto g_0_z = primBuffer.data(goff + 3 * idx + 2);

            #pragma omp simd aligned(pbx, pby, pbz, g_0_x, g_0_y, g_0_z,\
                                         g_0_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                double fact = g_0_0[k];

                g_0_x[k] = fact * pbx[k];

                g_0_y[k] = fact * pby[k];

                g_0_z[k] = fact * pbz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForPS(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 1, 0);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        // get position of integrals in primitves buffer

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|S)^(m+1) integrals

            auto g_0_0 = primBuffer.data(g1off + idx);

            // set up pointers to (P|g(r,r')|S)^(m) integrals

            auto g_x_0 = primBuffer.data(goff + 3 * idx);

            auto g_y_0 = primBuffer.data(goff + 3 * idx + 1);

            auto g_z_0 = primBuffer.data(goff + 3 * idx + 2);

            #pragma omp simd aligned(pax, pay, paz, g_x_0, g_y_0, g_z_0,\
                                         g_0_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                double fact = g_0_0[k];

                g_x_0[k] = fact * pax[k];

                g_y_0[k] = fact * pay[k];

                g_z_0[k] = fact * paz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForPP(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 1, 1);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|S)^(m+1) integrals

            auto gk_0_0 = primBuffer.data(gkoff + idx);

            // set up pointers to (S|g(r,r')|P)^(m+1) integrals

            auto g1_0_x = primBuffer.data(g1off + 3 * idx);

            auto g1_0_y = primBuffer.data(g1off + 3 * idx + 1);

            auto g1_0_z = primBuffer.data(g1off + 3 * idx + 2);

            // set up pointers to (P|g(r,r')|P)^(m) integrals

            auto g_x_x = primBuffer.data(goff + 9 * idx);

            auto g_x_y = primBuffer.data(goff + 9 * idx + 1);

            auto g_x_z = primBuffer.data(goff + 9 * idx + 2);

            auto g_y_x = primBuffer.data(goff + 9 * idx + 3);

            auto g_y_y = primBuffer.data(goff + 9 * idx + 4);

            auto g_y_z = primBuffer.data(goff + 9 * idx + 5);

            auto g_z_x = primBuffer.data(goff + 9 * idx + 6);

            auto g_z_y = primBuffer.data(goff + 9 * idx + 7);

            auto g_z_z = primBuffer.data(goff + 9 * idx + 8);

            #pragma omp simd aligned(fx, pax, pay, paz, gk_0_0, g1_0_x, g1_0_y,\
                                         g1_0_z, g_x_x, g_x_y, g_x_z, g_y_x, g_y_y,\
                                         g_y_z, g_z_x, g_z_y, g_z_z: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k] * gk_0_0[k];

                // leading x component

                double fr = pax[k];

                g_x_x[k] = fr * g1_0_x[k] + f2t;

                g_x_y[k] = fr * g1_0_y[k];

                g_x_z[k] = fr * g1_0_z[k];

                // leading y component

                fr = pay[k];

                g_y_x[k] = fr * g1_0_x[k];

                g_y_y[k] = fr * g1_0_y[k] + f2t;

                g_y_z[k] = fr * g1_0_z[k];

                // leading z component

                fr = paz[k];

                g_z_x[k] = fr * g1_0_x[k];

                g_z_y[k] = fr * g1_0_y[k];

                g_z_z[k] = fr * g1_0_z[k] + f2t;
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForSD(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  pbDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 0, 2);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 4);

            auto fz = osFactors.data(6 * idx + 5);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|S)^(m) integrals

            auto g20_0_0 = primBuffer.data(g20off + idx);

            // set up pointers to (S|g(r,r')|S)^(m+1) integrals

            auto g21_0_0 = primBuffer.data(g21off + idx);

            // set up pointers to (S|g(r,r')|P)^(m+1) integrals

            auto g1_0_x = primBuffer.data(g1off + 3 * idx);

            auto g1_0_y = primBuffer.data(g1off + 3 * idx + 1);

            auto g1_0_z = primBuffer.data(g1off + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|D)^(m) integrals

            auto g_0_xx = primBuffer.data(goff + 6 * idx);

            auto g_0_xy = primBuffer.data(goff + 6 * idx + 1);

            auto g_0_xz = primBuffer.data(goff + 6 * idx + 2);

            auto g_0_yy = primBuffer.data(goff + 6 * idx + 3);

            auto g_0_yz = primBuffer.data(goff + 6 * idx + 4);

            auto g_0_zz = primBuffer.data(goff + 6 * idx + 5);

            #pragma omp simd aligned(fg, fz, pbx, pby, pbz, g20_0_0, g21_0_0,\
                                         g1_0_x, g1_0_y, g1_0_z, g_0_xx, g_0_xy,\
                                         g_0_xz, g_0_yy, g_0_yz, g_0_zz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors

                double f2g = 0.50 * fg[k] * (g20_0_0[k] - fz[k] * g21_0_0[k]);

                // leading x component

                double fr = pbx[k];

                g_0_xx[k] = fr * g1_0_x[k] + f2g;

                g_0_xy[k] = fr * g1_0_y[k];

                g_0_xz[k] = fr * g1_0_z[k];

                // leading y component

                fr = pby[k];

                g_0_yy[k] = fr * g1_0_y[k] + f2g;

                g_0_yz[k] = fr * g1_0_z[k];

                // leading z component

                g_0_zz[k] = pbz[k] * g1_0_z[k] + f2g;
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForDS(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 2, 0);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|S)^(m) integrals

            auto g20_0_0 = primBuffer.data(g20off + idx);

            // set up pointers to (S|g(r,r')|S)^(m+1) integrals

            auto g21_0_0 = primBuffer.data(g21off + idx);

            // set up pointers to (P|g(r,r')|S)^(m+1) integrals

            auto g1_x_0 = primBuffer.data(g1off + 3 * idx);

            auto g1_y_0 = primBuffer.data(g1off + 3 * idx + 1);

            auto g1_z_0 = primBuffer.data(g1off + 3 * idx + 2);

            // set up pointers to (D|g(r,r')|S)^(m) integrals

            auto g_xx_0 = primBuffer.data(goff + 6 * idx);

            auto g_xy_0 = primBuffer.data(goff + 6 * idx + 1);

            auto g_xz_0 = primBuffer.data(goff + 6 * idx + 2);

            auto g_yy_0 = primBuffer.data(goff + 6 * idx + 3);

            auto g_yz_0 = primBuffer.data(goff + 6 * idx + 4);

            auto g_zz_0 = primBuffer.data(goff + 6 * idx + 5);

            #pragma omp simd aligned(fg, fz, pax, pay, paz, g20_0_0, g21_0_0,\
                                         g1_x_0, g1_y_0, g1_z_0, g_xx_0, g_xy_0,\
                                         g_xz_0, g_yy_0, g_yz_0, g_zz_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors for bra

                double f2g = 0.50 * fg[k] * (g20_0_0[k] - fz[k] * g21_0_0[k]);

                // leading x component

                double fr = pax[k];

                g_xx_0[k] = fr * g1_x_0[k] + f2g;

                g_xy_0[k] = fr * g1_y_0[k];

                g_xz_0[k] = fr * g1_z_0[k];

                // leading y component

                fr = pay[k];

                g_yy_0[k] = fr * g1_y_0[k] + f2g;

                g_yz_0[k] = fr * g1_z_0[k];

                // leading z component

                fr = paz[k];

                g_zz_0[k] = fr * g1_z_0[k] + f2g;
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForPD(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 1, 2);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|P)^(m+1) integrals

            auto gk_0_x = primBuffer.data(gkoff + 3 * idx);

            auto gk_0_y = primBuffer.data(gkoff + 3 * idx + 1);

            auto gk_0_z = primBuffer.data(gkoff + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|D)^(m+1) integrals

            auto g1_0_xx = primBuffer.data(g1off + 6 * idx);

            auto g1_0_xy = primBuffer.data(g1off + 6 * idx + 1);

            auto g1_0_xz = primBuffer.data(g1off + 6 * idx + 2);

            auto g1_0_yy = primBuffer.data(g1off + 6 * idx + 3);

            auto g1_0_yz = primBuffer.data(g1off + 6 * idx + 4);

            auto g1_0_zz = primBuffer.data(g1off + 6 * idx + 5);

            // set up pointers to (P|g(r,r')|D)^(m) integrals

            auto g_x_xx = primBuffer.data(goff + 18 * idx);

            auto g_x_xy = primBuffer.data(goff + 18 * idx + 1);

            auto g_x_xz = primBuffer.data(goff + 18 * idx + 2);

            auto g_x_yy = primBuffer.data(goff + 18 * idx + 3);

            auto g_x_yz = primBuffer.data(goff + 18 * idx + 4);

            auto g_x_zz = primBuffer.data(goff + 18 * idx + 5);

            auto g_y_xx = primBuffer.data(goff + 18 * idx + 6);

            auto g_y_xy = primBuffer.data(goff + 18 * idx + 7);

            auto g_y_xz = primBuffer.data(goff + 18 * idx + 8);

            auto g_y_yy = primBuffer.data(goff + 18 * idx + 9);

            auto g_y_yz = primBuffer.data(goff + 18 * idx + 10);

            auto g_y_zz = primBuffer.data(goff + 18 * idx + 11);

            auto g_z_xx = primBuffer.data(goff + 18 * idx + 12);

            auto g_z_xy = primBuffer.data(goff + 18 * idx + 13);

            auto g_z_xz = primBuffer.data(goff + 18 * idx + 14);

            auto g_z_yy = primBuffer.data(goff + 18 * idx + 15);

            auto g_z_yz = primBuffer.data(goff + 18 * idx + 16);

            auto g_z_zz = primBuffer.data(goff + 18 * idx + 17);

            #pragma omp simd aligned(fx, pax, pay, paz, gk_0_x, gk_0_y, gk_0_z,\
                                         g1_0_xx, g1_0_xy, g1_0_xz, g1_0_yy, g1_0_yz,\
                                         g1_0_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy,\
                                         g_x_yz, g_x_zz, g_y_xx, g_y_xy, g_y_xz,\
                                         g_y_yy, g_y_yz, g_y_zz, g_z_xx, g_z_xy,\
                                         g_z_xz, g_z_yy, g_z_yz, g_z_zz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // leading x component

                double fr = pax[k];

                g_x_xx[k] = fr * g1_0_xx[k] + 2.0 * f2t * gk_0_x[k];

                g_x_xy[k] = fr * g1_0_xy[k] + f2t * gk_0_y[k];

                g_x_xz[k] = fr * g1_0_xz[k] + f2t * gk_0_z[k];

                g_x_yy[k] = fr * g1_0_yy[k];

                g_x_yz[k] = fr * g1_0_yz[k];

                g_x_zz[k] = fr * g1_0_zz[k];

                // leading y component

                fr = pay[k];

                g_y_xx[k] = fr * g1_0_xx[k];

                g_y_xy[k] = fr * g1_0_xy[k] + f2t * gk_0_x[k];

                g_y_xz[k] = fr * g1_0_xz[k];

                g_y_yy[k] = fr * g1_0_yy[k] + 2.0 * f2t * gk_0_y[k];

                g_y_yz[k] = fr * g1_0_yz[k] + f2t * gk_0_z[k];

                g_y_zz[k] = fr * g1_0_zz[k];

                // leading z component

                fr = paz[k];

                g_z_xx[k] = fr * g1_0_xx[k];

                g_z_xy[k] = fr * g1_0_xy[k];

                g_z_xz[k] = fr * g1_0_xz[k] + f2t * gk_0_x[k];

                g_z_yy[k] = fr * g1_0_yy[k];

                g_z_yz[k] = fr * g1_0_yz[k] + f2t * gk_0_y[k];

                g_z_zz[k] = fr * g1_0_zz[k] + 2.0 * f2t * gk_0_z[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForDP(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 2, 1);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|g(r,r')|S)^(m+1) integrals

            auto gk_x_0 = primBuffer.data(gkoff + 3 * idx);

            auto gk_y_0 = primBuffer.data(gkoff + 3 * idx + 1);

            auto gk_z_0 = primBuffer.data(gkoff + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|P)^(m) integrals

            auto g20_0_x = primBuffer.data(g20off + 3 * idx);

            auto g20_0_y = primBuffer.data(g20off + 3 * idx + 1);

            auto g20_0_z = primBuffer.data(g20off + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|P)^(m+1) integrals

            auto g21_0_x = primBuffer.data(g21off + 3 * idx);

            auto g21_0_y = primBuffer.data(g21off + 3 * idx + 1);

            auto g21_0_z = primBuffer.data(g21off + 3 * idx + 2);

            // set up pointers to (P|g(r,r')|P)^(m+1) integrals

            auto g1_x_x = primBuffer.data(g1off + 9 * idx);

            auto g1_x_y = primBuffer.data(g1off + 9 * idx + 1);

            auto g1_x_z = primBuffer.data(g1off + 9 * idx + 2);

            auto g1_y_x = primBuffer.data(g1off + 9 * idx + 3);

            auto g1_y_y = primBuffer.data(g1off + 9 * idx + 4);

            auto g1_y_z = primBuffer.data(g1off + 9 * idx + 5);

            auto g1_z_x = primBuffer.data(g1off + 9 * idx + 6);

            auto g1_z_y = primBuffer.data(g1off + 9 * idx + 7);

            auto g1_z_z = primBuffer.data(g1off + 9 * idx + 8);

            // set up pointers to (D|g(r,r')|P)^(m) integrals

            auto g_xx_x = primBuffer.data(goff + 18 * idx);

            auto g_xx_y = primBuffer.data(goff + 18 * idx + 1);

            auto g_xx_z = primBuffer.data(goff + 18 * idx + 2);

            auto g_xy_x = primBuffer.data(goff + 18 * idx + 3);

            auto g_xy_y = primBuffer.data(goff + 18 * idx + 4);

            auto g_xy_z = primBuffer.data(goff + 18 * idx + 5);

            auto g_xz_x = primBuffer.data(goff + 18 * idx + 6);

            auto g_xz_y = primBuffer.data(goff + 18 * idx + 7);

            auto g_xz_z = primBuffer.data(goff + 18 * idx + 8);

            auto g_yy_x = primBuffer.data(goff + 18 * idx + 9);

            auto g_yy_y = primBuffer.data(goff + 18 * idx + 10);

            auto g_yy_z = primBuffer.data(goff + 18 * idx + 11);

            auto g_yz_x = primBuffer.data(goff + 18 * idx + 12);

            auto g_yz_y = primBuffer.data(goff + 18 * idx + 13);

            auto g_yz_z = primBuffer.data(goff + 18 * idx + 14);

            auto g_zz_x = primBuffer.data(goff + 18 * idx + 15);

            auto g_zz_y = primBuffer.data(goff + 18 * idx + 16);

            auto g_zz_z = primBuffer.data(goff + 18 * idx + 17);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_x_0, gk_y_0,\
                                         gk_z_0, g20_0_x, g20_0_y, g20_0_z, g21_0_x,\
                                         g21_0_y, g21_0_z, g1_x_x, g1_x_y, g1_x_z,\
                                         g1_y_x, g1_y_y, g1_y_z, g1_z_x, g1_z_y,\
                                         g1_z_z, g_xx_x, g_xx_y, g_xx_z, g_xy_x,\
                                         g_xy_y, g_xy_z, g_xz_x, g_xz_y, g_xz_z,\
                                         g_yy_x, g_yy_y, g_yy_z, g_yz_x, g_yz_y,\
                                         g_yz_z, g_zz_x, g_zz_y, g_zz_z: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xx_x[k] = fr * g1_x_x[k] + f2g * (g20_0_x[k] - fgz * g21_0_x[k]) + f2t * gk_x_0[k];

                g_xx_y[k] = fr * g1_x_y[k] + f2g * (g20_0_y[k] - fgz * g21_0_y[k]);

                g_xx_z[k] = fr * g1_x_z[k] + f2g * (g20_0_z[k] - fgz * g21_0_z[k]);

                g_xy_x[k] = fr * g1_y_x[k] + f2t * gk_y_0[k];

                g_xy_y[k] = fr * g1_y_y[k];

                g_xy_z[k] = fr * g1_y_z[k];

                g_xz_x[k] = fr * g1_z_x[k] + f2t * gk_z_0[k];

                g_xz_y[k] = fr * g1_z_y[k];

                g_xz_z[k] = fr * g1_z_z[k];

                // leading y component

                fr = pay[k];

                g_yy_x[k] = fr * g1_y_x[k] + f2g * (g20_0_x[k] - fgz * g21_0_x[k]);

                g_yy_y[k] = fr * g1_y_y[k] + f2g * (g20_0_y[k] - fgz * g21_0_y[k]) + f2t * gk_y_0[k];

                g_yy_z[k] = fr * g1_y_z[k] + f2g * (g20_0_z[k] - fgz * g21_0_z[k]);

                g_yz_x[k] = fr * g1_z_x[k];

                g_yz_y[k] = fr * g1_z_y[k] + f2t * gk_z_0[k];

                g_yz_z[k] = fr * g1_z_z[k];

                // leading z component

                fr = paz[k];

                g_zz_x[k] = fr * g1_z_x[k] + f2g * (g20_0_x[k] - fgz * g21_0_x[k]);

                g_zz_y[k] = fr * g1_z_y[k] + f2g * (g20_0_y[k] - fgz * g21_0_y[k]);

                g_zz_z[k] = fr * g1_z_z[k] + f2g * (g20_0_z[k] - fgz * g21_0_z[k]) + f2t * gk_z_0[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForDD(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 2, 2);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|g(r,r')|P)^(m+1) integrals

            auto gk_x_x = primBuffer.data(gkoff + 9 * idx);

            auto gk_x_y = primBuffer.data(gkoff + 9 * idx + 1);

            auto gk_x_z = primBuffer.data(gkoff + 9 * idx + 2);

            auto gk_y_x = primBuffer.data(gkoff + 9 * idx + 3);

            auto gk_y_y = primBuffer.data(gkoff + 9 * idx + 4);

            auto gk_y_z = primBuffer.data(gkoff + 9 * idx + 5);

            auto gk_z_x = primBuffer.data(gkoff + 9 * idx + 6);

            auto gk_z_y = primBuffer.data(gkoff + 9 * idx + 7);

            auto gk_z_z = primBuffer.data(gkoff + 9 * idx + 8);

            // set up pointers to (S|g(r,r')|D)^(m) integrals

            auto g20_0_xx = primBuffer.data(g20off + 6 * idx);

            auto g20_0_xy = primBuffer.data(g20off + 6 * idx + 1);

            auto g20_0_xz = primBuffer.data(g20off + 6 * idx + 2);

            auto g20_0_yy = primBuffer.data(g20off + 6 * idx + 3);

            auto g20_0_yz = primBuffer.data(g20off + 6 * idx + 4);

            auto g20_0_zz = primBuffer.data(g20off + 6 * idx + 5);

            // set up pointers to (S|g(r,r')|D)^(m+1) integrals

            auto g21_0_xx = primBuffer.data(g21off + 6 * idx);

            auto g21_0_xy = primBuffer.data(g21off + 6 * idx + 1);

            auto g21_0_xz = primBuffer.data(g21off + 6 * idx + 2);

            auto g21_0_yy = primBuffer.data(g21off + 6 * idx + 3);

            auto g21_0_yz = primBuffer.data(g21off + 6 * idx + 4);

            auto g21_0_zz = primBuffer.data(g21off + 6 * idx + 5);

            // set up pointers to (P|g(r,r')|D)^(m+1) integrals

            auto g1_x_xx = primBuffer.data(g1off + 18 * idx);

            auto g1_x_xy = primBuffer.data(g1off + 18 * idx + 1);

            auto g1_x_xz = primBuffer.data(g1off + 18 * idx + 2);

            auto g1_x_yy = primBuffer.data(g1off + 18 * idx + 3);

            auto g1_x_yz = primBuffer.data(g1off + 18 * idx + 4);

            auto g1_x_zz = primBuffer.data(g1off + 18 * idx + 5);

            auto g1_y_xx = primBuffer.data(g1off + 18 * idx + 6);

            auto g1_y_xy = primBuffer.data(g1off + 18 * idx + 7);

            auto g1_y_xz = primBuffer.data(g1off + 18 * idx + 8);

            auto g1_y_yy = primBuffer.data(g1off + 18 * idx + 9);

            auto g1_y_yz = primBuffer.data(g1off + 18 * idx + 10);

            auto g1_y_zz = primBuffer.data(g1off + 18 * idx + 11);

            auto g1_z_xx = primBuffer.data(g1off + 18 * idx + 12);

            auto g1_z_xy = primBuffer.data(g1off + 18 * idx + 13);

            auto g1_z_xz = primBuffer.data(g1off + 18 * idx + 14);

            auto g1_z_yy = primBuffer.data(g1off + 18 * idx + 15);

            auto g1_z_yz = primBuffer.data(g1off + 18 * idx + 16);

            auto g1_z_zz = primBuffer.data(g1off + 18 * idx + 17);

            // set up pointers to (D|g(r,r')|D)^(m) integrals

            auto g_xx_xx = primBuffer.data(goff + 36 * idx);

            auto g_xx_xy = primBuffer.data(goff + 36 * idx + 1);

            auto g_xx_xz = primBuffer.data(goff + 36 * idx + 2);

            auto g_xx_yy = primBuffer.data(goff + 36 * idx + 3);

            auto g_xx_yz = primBuffer.data(goff + 36 * idx + 4);

            auto g_xx_zz = primBuffer.data(goff + 36 * idx + 5);

            auto g_xy_xx = primBuffer.data(goff + 36 * idx + 6);

            auto g_xy_xy = primBuffer.data(goff + 36 * idx + 7);

            auto g_xy_xz = primBuffer.data(goff + 36 * idx + 8);

            auto g_xy_yy = primBuffer.data(goff + 36 * idx + 9);

            auto g_xy_yz = primBuffer.data(goff + 36 * idx + 10);

            auto g_xy_zz = primBuffer.data(goff + 36 * idx + 11);

            auto g_xz_xx = primBuffer.data(goff + 36 * idx + 12);

            auto g_xz_xy = primBuffer.data(goff + 36 * idx + 13);

            auto g_xz_xz = primBuffer.data(goff + 36 * idx + 14);

            auto g_xz_yy = primBuffer.data(goff + 36 * idx + 15);

            auto g_xz_yz = primBuffer.data(goff + 36 * idx + 16);

            auto g_xz_zz = primBuffer.data(goff + 36 * idx + 17);

            auto g_yy_xx = primBuffer.data(goff + 36 * idx + 18);

            auto g_yy_xy = primBuffer.data(goff + 36 * idx + 19);

            auto g_yy_xz = primBuffer.data(goff + 36 * idx + 20);

            auto g_yy_yy = primBuffer.data(goff + 36 * idx + 21);

            auto g_yy_yz = primBuffer.data(goff + 36 * idx + 22);

            auto g_yy_zz = primBuffer.data(goff + 36 * idx + 23);

            auto g_yz_xx = primBuffer.data(goff + 36 * idx + 24);

            auto g_yz_xy = primBuffer.data(goff + 36 * idx + 25);

            auto g_yz_xz = primBuffer.data(goff + 36 * idx + 26);

            auto g_yz_yy = primBuffer.data(goff + 36 * idx + 27);

            auto g_yz_yz = primBuffer.data(goff + 36 * idx + 28);

            auto g_yz_zz = primBuffer.data(goff + 36 * idx + 29);

            auto g_zz_xx = primBuffer.data(goff + 36 * idx + 30);

            auto g_zz_xy = primBuffer.data(goff + 36 * idx + 31);

            auto g_zz_xz = primBuffer.data(goff + 36 * idx + 32);

            auto g_zz_yy = primBuffer.data(goff + 36 * idx + 33);

            auto g_zz_yz = primBuffer.data(goff + 36 * idx + 34);

            auto g_zz_zz = primBuffer.data(goff + 36 * idx + 35);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_x_x, gk_x_y,\
                                         gk_x_z, gk_y_x, gk_y_y, gk_y_z, gk_z_x,\
                                         gk_z_y, gk_z_z, g20_0_xx, g20_0_xy, g20_0_xz,\
                                         g20_0_yy, g20_0_yz, g20_0_zz, g21_0_xx,\
                                         g21_0_xy, g21_0_xz, g21_0_yy, g21_0_yz,\
                                         g21_0_zz, g1_x_xx, g1_x_xy, g1_x_xz, g1_x_yy,\
                                         g1_x_yz, g1_x_zz, g1_y_xx, g1_y_xy, g1_y_xz,\
                                         g1_y_yy, g1_y_yz, g1_y_zz, g1_z_xx, g1_z_xy,\
                                         g1_z_xz, g1_z_yy, g1_z_yz, g1_z_zz, g_xx_xx,\
                                         g_xx_xy, g_xx_xz, g_xx_yy, g_xx_yz, g_xx_zz,\
                                         g_xy_xx, g_xy_xy, g_xy_xz, g_xy_yy, g_xy_yz,\
                                         g_xy_zz, g_xz_xx, g_xz_xy, g_xz_xz, g_xz_yy,\
                                         g_xz_yz, g_xz_zz, g_yy_xx, g_yy_xy, g_yy_xz,\
                                         g_yy_yy, g_yy_yz, g_yy_zz, g_yz_xx, g_yz_xy,\
                                         g_yz_xz, g_yz_yy, g_yz_yz, g_yz_zz, g_zz_xx,\
                                         g_zz_xy, g_zz_xz, g_zz_yy, g_zz_yz, g_zz_zz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xx_xx[k] = fr * g1_x_xx[k] + f2g * (g20_0_xx[k] - fgz * g21_0_xx[k]) + 2.0 * f2t * gk_x_x[k];

                g_xx_xy[k] = fr * g1_x_xy[k] + f2g * (g20_0_xy[k] - fgz * g21_0_xy[k]) + f2t * gk_x_y[k];

                g_xx_xz[k] = fr * g1_x_xz[k] + f2g * (g20_0_xz[k] - fgz * g21_0_xz[k]) + f2t * gk_x_z[k];

                g_xx_yy[k] = fr * g1_x_yy[k] + f2g * (g20_0_yy[k] - fgz * g21_0_yy[k]);

                g_xx_yz[k] = fr * g1_x_yz[k] + f2g * (g20_0_yz[k] - fgz * g21_0_yz[k]);

                g_xx_zz[k] = fr * g1_x_zz[k] + f2g * (g20_0_zz[k] - fgz * g21_0_zz[k]);

                g_xy_xx[k] = fr * g1_y_xx[k] + 2.0 * f2t * gk_y_x[k];

                g_xy_xy[k] = fr * g1_y_xy[k] + f2t * gk_y_y[k];

                g_xy_xz[k] = fr * g1_y_xz[k] + f2t * gk_y_z[k];

                g_xy_yy[k] = fr * g1_y_yy[k];

                g_xy_yz[k] = fr * g1_y_yz[k];

                g_xy_zz[k] = fr * g1_y_zz[k];

                g_xz_xx[k] = fr * g1_z_xx[k] + 2.0 * f2t * gk_z_x[k];

                g_xz_xy[k] = fr * g1_z_xy[k] + f2t * gk_z_y[k];

                g_xz_xz[k] = fr * g1_z_xz[k] + f2t * gk_z_z[k];

                g_xz_yy[k] = fr * g1_z_yy[k];

                g_xz_yz[k] = fr * g1_z_yz[k];

                g_xz_zz[k] = fr * g1_z_zz[k];

                // leading y component

                fr = pay[k];

                g_yy_xx[k] = fr * g1_y_xx[k] + f2g * (g20_0_xx[k] - fgz * g21_0_xx[k]);

                g_yy_xy[k] = fr * g1_y_xy[k] + f2g * (g20_0_xy[k] - fgz * g21_0_xy[k]) + f2t * gk_y_x[k];

                g_yy_xz[k] = fr * g1_y_xz[k] + f2g * (g20_0_xz[k] - fgz * g21_0_xz[k]);

                g_yy_yy[k] = fr * g1_y_yy[k] + f2g * (g20_0_yy[k] - fgz * g21_0_yy[k]) + 2.0 * f2t * gk_y_y[k];

                g_yy_yz[k] = fr * g1_y_yz[k] + f2g * (g20_0_yz[k] - fgz * g21_0_yz[k]) + f2t * gk_y_z[k];

                g_yy_zz[k] = fr * g1_y_zz[k] + f2g * (g20_0_zz[k] - fgz * g21_0_zz[k]);

                g_yz_xx[k] = fr * g1_z_xx[k];

                g_yz_xy[k] = fr * g1_z_xy[k] + f2t * gk_z_x[k];

                g_yz_xz[k] = fr * g1_z_xz[k];

                g_yz_yy[k] = fr * g1_z_yy[k] + 2.0 * f2t * gk_z_y[k];

                g_yz_yz[k] = fr * g1_z_yz[k] + f2t * gk_z_z[k];

                g_yz_zz[k] = fr * g1_z_zz[k];

                // leading z component

                fr = paz[k];

                g_zz_xx[k] = fr * g1_z_xx[k] + f2g * (g20_0_xx[k] - fgz * g21_0_xx[k]);

                g_zz_xy[k] = fr * g1_z_xy[k] + f2g * (g20_0_xy[k] - fgz * g21_0_xy[k]);

                g_zz_xz[k] = fr * g1_z_xz[k] + f2g * (g20_0_xz[k] - fgz * g21_0_xz[k]) + f2t * gk_z_x[k];

                g_zz_yy[k] = fr * g1_z_yy[k] + f2g * (g20_0_yy[k] - fgz * g21_0_yy[k]);

                g_zz_yz[k] = fr * g1_z_yz[k] + f2g * (g20_0_yz[k] - fgz * g21_0_yz[k]) + f2t * gk_z_y[k];

                g_zz_zz[k] = fr * g1_z_zz[k] + f2g * (g20_0_zz[k] - fgz * g21_0_zz[k]) + 2.0 * f2t * gk_z_z[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForSF(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  pbDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 0, 3);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 1, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 4);

            auto fz = osFactors.data(6 * idx + 5);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|P)^(m) integrals

            auto g20_0_x = primBuffer.data(g20off + 3 * idx);

            auto g20_0_y = primBuffer.data(g20off + 3 * idx + 1);

            auto g20_0_z = primBuffer.data(g20off + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|P)^(m+1) integrals

            auto g21_0_x = primBuffer.data(g21off + 3 * idx);

            auto g21_0_y = primBuffer.data(g21off + 3 * idx + 1);

            auto g21_0_z = primBuffer.data(g21off + 3 * idx + 2);

            // set up pointers to (S|g(r,r')|D)^(m+1) integrals

            auto g1_0_xx = primBuffer.data(g1off + 6 * idx);

            auto g1_0_xy = primBuffer.data(g1off + 6 * idx + 1);

            auto g1_0_xz = primBuffer.data(g1off + 6 * idx + 2);

            auto g1_0_yy = primBuffer.data(g1off + 6 * idx + 3);

            auto g1_0_yz = primBuffer.data(g1off + 6 * idx + 4);

            auto g1_0_zz = primBuffer.data(g1off + 6 * idx + 5);

            // set up pointers to (S|g(r,r')|F)^(m) integrals

            auto g_0_xxx = primBuffer.data(goff + 10 * idx);

            auto g_0_xxy = primBuffer.data(goff + 10 * idx + 1);

            auto g_0_xxz = primBuffer.data(goff + 10 * idx + 2);

            auto g_0_xyy = primBuffer.data(goff + 10 * idx + 3);

            auto g_0_xyz = primBuffer.data(goff + 10 * idx + 4);

            auto g_0_xzz = primBuffer.data(goff + 10 * idx + 5);

            auto g_0_yyy = primBuffer.data(goff + 10 * idx + 6);

            auto g_0_yyz = primBuffer.data(goff + 10 * idx + 7);

            auto g_0_yzz = primBuffer.data(goff + 10 * idx + 8);

            auto g_0_zzz = primBuffer.data(goff + 10 * idx + 9);

            #pragma omp simd aligned(fg, fz, pbx, pby, pbz, g20_0_x, g20_0_y,\
                                         g20_0_z, g21_0_x, g21_0_y, g21_0_z, g1_0_xx,\
                                         g1_0_xy, g1_0_xz, g1_0_yy, g1_0_yz, g1_0_zz,\
                                         g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz,\
                                         g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pbx[k];

                g_0_xxx[k] = fr * g1_0_xx[k] + f2g * (2.0 * g20_0_x[k] - 2.0 * fgz * g21_0_x[k]);

                g_0_xxy[k] = fr * g1_0_xy[k] + f2g * (g20_0_y[k] - fgz * g21_0_y[k]);

                g_0_xxz[k] = fr * g1_0_xz[k] + f2g * (g20_0_z[k] - fgz * g21_0_z[k]);

                g_0_xyy[k] = fr * g1_0_yy[k];

                g_0_xyz[k] = fr * g1_0_yz[k];

                g_0_xzz[k] = fr * g1_0_zz[k];

                // leading y component

                fr = pby[k];

                g_0_yyy[k] = fr * g1_0_yy[k] + f2g * (2.0 * g20_0_y[k] - 2.0 * fgz * g21_0_y[k]);

                g_0_yyz[k] = fr * g1_0_yz[k] + f2g * (g20_0_z[k] - fgz * g21_0_z[k]);

                g_0_yzz[k] = fr * g1_0_zz[k];

                // leading z component

                g_0_zzz[k] = pbz[k] * g1_0_zz[k] + f2g * (2.0 * g20_0_z[k] - 2.0 * fgz * g21_0_z[k]);
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForFS(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 3, 0);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|g(r,r')|S)^(m) integrals

            auto g20_x_0 = primBuffer.data(g20off + 3 * idx);

            auto g20_y_0 = primBuffer.data(g20off + 3 * idx + 1);

            auto g20_z_0 = primBuffer.data(g20off + 3 * idx + 2);

            // set up pointers to (P|g(r,r')|S)^(m+1) integrals

            auto g21_x_0 = primBuffer.data(g21off + 3 * idx);

            auto g21_y_0 = primBuffer.data(g21off + 3 * idx + 1);

            auto g21_z_0 = primBuffer.data(g21off + 3 * idx + 2);

            // set up pointers to (D|g(r,r')|S)^(m+1) integrals

            auto g1_xx_0 = primBuffer.data(g1off + 6 * idx);

            auto g1_xy_0 = primBuffer.data(g1off + 6 * idx + 1);

            auto g1_xz_0 = primBuffer.data(g1off + 6 * idx + 2);

            auto g1_yy_0 = primBuffer.data(g1off + 6 * idx + 3);

            auto g1_yz_0 = primBuffer.data(g1off + 6 * idx + 4);

            auto g1_zz_0 = primBuffer.data(g1off + 6 * idx + 5);

            // set up pointers to (F|g(r,r')|S)^(m) integrals

            auto g_xxx_0 = primBuffer.data(goff + 10 * idx);

            auto g_xxy_0 = primBuffer.data(goff + 10 * idx + 1);

            auto g_xxz_0 = primBuffer.data(goff + 10 * idx + 2);

            auto g_xyy_0 = primBuffer.data(goff + 10 * idx + 3);

            auto g_xyz_0 = primBuffer.data(goff + 10 * idx + 4);

            auto g_xzz_0 = primBuffer.data(goff + 10 * idx + 5);

            auto g_yyy_0 = primBuffer.data(goff + 10 * idx + 6);

            auto g_yyz_0 = primBuffer.data(goff + 10 * idx + 7);

            auto g_yzz_0 = primBuffer.data(goff + 10 * idx + 8);

            auto g_zzz_0 = primBuffer.data(goff + 10 * idx + 9);

            #pragma omp simd aligned(fg, fz, pax, pay, paz, g20_x_0, g20_y_0,\
                                         g20_z_0, g21_x_0, g21_y_0, g21_z_0, g1_xx_0,\
                                         g1_xy_0, g1_xz_0, g1_yy_0, g1_yz_0, g1_zz_0,\
                                         g_xxx_0, g_xxy_0, g_xxz_0, g_xyy_0, g_xyz_0,\
                                         g_xzz_0, g_yyy_0, g_yyz_0, g_yzz_0, g_zzz_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxx_0[k] = fr * g1_xx_0[k] + f2g * (2.0 * g20_x_0[k] - 2.0 * fgz * g21_x_0[k]);

                g_xxy_0[k] = fr * g1_xy_0[k] + f2g * (g20_y_0[k] - fgz * g21_y_0[k]);

                g_xxz_0[k] = fr * g1_xz_0[k] + f2g * (g20_z_0[k] - fgz * g21_z_0[k]);

                g_xyy_0[k] = fr * g1_yy_0[k];

                g_xyz_0[k] = fr * g1_yz_0[k];

                g_xzz_0[k] = fr * g1_zz_0[k];

                // leading y component

                fr = pay[k];

                g_yyy_0[k] = fr * g1_yy_0[k] + f2g * (2.0 * g20_y_0[k] - 2.0 * fgz * g21_y_0[k]);

                g_yyz_0[k] = fr * g1_yz_0[k] + f2g * (g20_z_0[k] - fgz * g21_z_0[k]);

                g_yzz_0[k] = fr * g1_zz_0[k];

                // leading z component

                fr = paz[k];

                g_zzz_0[k] = fr * g1_zz_0[k] + f2g * (2.0 * g20_z_0[k] - 2.0 * fgz * g21_z_0[k]);
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForPF(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 1, 3);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|D)^(m+1) integrals

            auto gk_0_xx = primBuffer.data(gkoff + 6 * idx);

            auto gk_0_xy = primBuffer.data(gkoff + 6 * idx + 1);

            auto gk_0_xz = primBuffer.data(gkoff + 6 * idx + 2);

            auto gk_0_yy = primBuffer.data(gkoff + 6 * idx + 3);

            auto gk_0_yz = primBuffer.data(gkoff + 6 * idx + 4);

            auto gk_0_zz = primBuffer.data(gkoff + 6 * idx + 5);

            // set up pointers to (S|g(r,r')|F)^(m+1) integrals

            auto g1_0_xxx = primBuffer.data(g1off + 10 * idx);

            auto g1_0_xxy = primBuffer.data(g1off + 10 * idx + 1);

            auto g1_0_xxz = primBuffer.data(g1off + 10 * idx + 2);

            auto g1_0_xyy = primBuffer.data(g1off + 10 * idx + 3);

            auto g1_0_xyz = primBuffer.data(g1off + 10 * idx + 4);

            auto g1_0_xzz = primBuffer.data(g1off + 10 * idx + 5);

            auto g1_0_yyy = primBuffer.data(g1off + 10 * idx + 6);

            auto g1_0_yyz = primBuffer.data(g1off + 10 * idx + 7);

            auto g1_0_yzz = primBuffer.data(g1off + 10 * idx + 8);

            auto g1_0_zzz = primBuffer.data(g1off + 10 * idx + 9);

            // set up pointers to (P|g(r,r')|F)^(m) integrals

            auto g_x_xxx = primBuffer.data(goff + 30 * idx);

            auto g_x_xxy = primBuffer.data(goff + 30 * idx + 1);

            auto g_x_xxz = primBuffer.data(goff + 30 * idx + 2);

            auto g_x_xyy = primBuffer.data(goff + 30 * idx + 3);

            auto g_x_xyz = primBuffer.data(goff + 30 * idx + 4);

            auto g_x_xzz = primBuffer.data(goff + 30 * idx + 5);

            auto g_x_yyy = primBuffer.data(goff + 30 * idx + 6);

            auto g_x_yyz = primBuffer.data(goff + 30 * idx + 7);

            auto g_x_yzz = primBuffer.data(goff + 30 * idx + 8);

            auto g_x_zzz = primBuffer.data(goff + 30 * idx + 9);

            auto g_y_xxx = primBuffer.data(goff + 30 * idx + 10);

            auto g_y_xxy = primBuffer.data(goff + 30 * idx + 11);

            auto g_y_xxz = primBuffer.data(goff + 30 * idx + 12);

            auto g_y_xyy = primBuffer.data(goff + 30 * idx + 13);

            auto g_y_xyz = primBuffer.data(goff + 30 * idx + 14);

            auto g_y_xzz = primBuffer.data(goff + 30 * idx + 15);

            auto g_y_yyy = primBuffer.data(goff + 30 * idx + 16);

            auto g_y_yyz = primBuffer.data(goff + 30 * idx + 17);

            auto g_y_yzz = primBuffer.data(goff + 30 * idx + 18);

            auto g_y_zzz = primBuffer.data(goff + 30 * idx + 19);

            auto g_z_xxx = primBuffer.data(goff + 30 * idx + 20);

            auto g_z_xxy = primBuffer.data(goff + 30 * idx + 21);

            auto g_z_xxz = primBuffer.data(goff + 30 * idx + 22);

            auto g_z_xyy = primBuffer.data(goff + 30 * idx + 23);

            auto g_z_xyz = primBuffer.data(goff + 30 * idx + 24);

            auto g_z_xzz = primBuffer.data(goff + 30 * idx + 25);

            auto g_z_yyy = primBuffer.data(goff + 30 * idx + 26);

            auto g_z_yyz = primBuffer.data(goff + 30 * idx + 27);

            auto g_z_yzz = primBuffer.data(goff + 30 * idx + 28);

            auto g_z_zzz = primBuffer.data(goff + 30 * idx + 29);

            #pragma omp simd aligned(fx, pax, pay, paz, gk_0_xx, gk_0_xy, gk_0_xz,\
                                         gk_0_yy, gk_0_yz, gk_0_zz, g1_0_xxx, g1_0_xxy,\
                                         g1_0_xxz, g1_0_xyy, g1_0_xyz, g1_0_xzz,\
                                         g1_0_yyy, g1_0_yyz, g1_0_yzz, g1_0_zzz,\
                                         g_x_xxx, g_x_xxy, g_x_xxz, g_x_xyy, g_x_xyz,\
                                         g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz, g_x_zzz,\
                                         g_y_xxx, g_y_xxy, g_y_xxz, g_y_xyy, g_y_xyz,\
                                         g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz, g_y_zzz,\
                                         g_z_xxx, g_z_xxy, g_z_xxz, g_z_xyy, g_z_xyz,\
                                         g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz, g_z_zzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // leading x component

                double fr = pax[k];

                g_x_xxx[k] = fr * g1_0_xxx[k] + 3.0 * f2t * gk_0_xx[k];

                g_x_xxy[k] = fr * g1_0_xxy[k] + 2.0 * f2t * gk_0_xy[k];

                g_x_xxz[k] = fr * g1_0_xxz[k] + 2.0 * f2t * gk_0_xz[k];

                g_x_xyy[k] = fr * g1_0_xyy[k] + f2t * gk_0_yy[k];

                g_x_xyz[k] = fr * g1_0_xyz[k] + f2t * gk_0_yz[k];

                g_x_xzz[k] = fr * g1_0_xzz[k] + f2t * gk_0_zz[k];

                g_x_yyy[k] = fr * g1_0_yyy[k];

                g_x_yyz[k] = fr * g1_0_yyz[k];

                g_x_yzz[k] = fr * g1_0_yzz[k];

                g_x_zzz[k] = fr * g1_0_zzz[k];

                // leading y component

                fr = pay[k];

                g_y_xxx[k] = fr * g1_0_xxx[k];

                g_y_xxy[k] = fr * g1_0_xxy[k] + f2t * gk_0_xx[k];

                g_y_xxz[k] = fr * g1_0_xxz[k];

                g_y_xyy[k] = fr * g1_0_xyy[k] + 2.0 * f2t * gk_0_xy[k];

                g_y_xyz[k] = fr * g1_0_xyz[k] + f2t * gk_0_xz[k];

                g_y_xzz[k] = fr * g1_0_xzz[k];

                g_y_yyy[k] = fr * g1_0_yyy[k] + 3.0 * f2t * gk_0_yy[k];

                g_y_yyz[k] = fr * g1_0_yyz[k] + 2.0 * f2t * gk_0_yz[k];

                g_y_yzz[k] = fr * g1_0_yzz[k] + f2t * gk_0_zz[k];

                g_y_zzz[k] = fr * g1_0_zzz[k];

                // leading z component

                fr = paz[k];

                g_z_xxx[k] = fr * g1_0_xxx[k];

                g_z_xxy[k] = fr * g1_0_xxy[k];

                g_z_xxz[k] = fr * g1_0_xxz[k] + f2t * gk_0_xx[k];

                g_z_xyy[k] = fr * g1_0_xyy[k];

                g_z_xyz[k] = fr * g1_0_xyz[k] + f2t * gk_0_xy[k];

                g_z_xzz[k] = fr * g1_0_xzz[k] + 2.0 * f2t * gk_0_xz[k];

                g_z_yyy[k] = fr * g1_0_yyy[k];

                g_z_yyz[k] = fr * g1_0_yyz[k] + f2t * gk_0_yy[k];

                g_z_yzz[k] = fr * g1_0_yzz[k] + 2.0 * f2t * gk_0_yz[k];

                g_z_zzz[k] = fr * g1_0_zzz[k] + 3.0 * f2t * gk_0_zz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForFP(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 3, 1);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 1, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|g(r,r')|S)^(m+1) integrals

            auto gk_xx_0 = primBuffer.data(gkoff + 6 * idx);

            auto gk_xy_0 = primBuffer.data(gkoff + 6 * idx + 1);

            auto gk_xz_0 = primBuffer.data(gkoff + 6 * idx + 2);

            auto gk_yy_0 = primBuffer.data(gkoff + 6 * idx + 3);

            auto gk_yz_0 = primBuffer.data(gkoff + 6 * idx + 4);

            auto gk_zz_0 = primBuffer.data(gkoff + 6 * idx + 5);

            // set up pointers to (P|g(r,r')|P)^(m) integrals

            auto g20_x_x = primBuffer.data(g20off + 9 * idx);

            auto g20_x_y = primBuffer.data(g20off + 9 * idx + 1);

            auto g20_x_z = primBuffer.data(g20off + 9 * idx + 2);

            auto g20_y_x = primBuffer.data(g20off + 9 * idx + 3);

            auto g20_y_y = primBuffer.data(g20off + 9 * idx + 4);

            auto g20_y_z = primBuffer.data(g20off + 9 * idx + 5);

            auto g20_z_x = primBuffer.data(g20off + 9 * idx + 6);

            auto g20_z_y = primBuffer.data(g20off + 9 * idx + 7);

            auto g20_z_z = primBuffer.data(g20off + 9 * idx + 8);

            // set up pointers to (P|g(r,r')|P)^(m+1) integrals

            auto g21_x_x = primBuffer.data(g21off + 9 * idx);

            auto g21_x_y = primBuffer.data(g21off + 9 * idx + 1);

            auto g21_x_z = primBuffer.data(g21off + 9 * idx + 2);

            auto g21_y_x = primBuffer.data(g21off + 9 * idx + 3);

            auto g21_y_y = primBuffer.data(g21off + 9 * idx + 4);

            auto g21_y_z = primBuffer.data(g21off + 9 * idx + 5);

            auto g21_z_x = primBuffer.data(g21off + 9 * idx + 6);

            auto g21_z_y = primBuffer.data(g21off + 9 * idx + 7);

            auto g21_z_z = primBuffer.data(g21off + 9 * idx + 8);

            // set up pointers to (D|g(r,r')|P)^(m+1) integrals

            auto g1_xx_x = primBuffer.data(g1off + 18 * idx);

            auto g1_xx_y = primBuffer.data(g1off + 18 * idx + 1);

            auto g1_xx_z = primBuffer.data(g1off + 18 * idx + 2);

            auto g1_xy_x = primBuffer.data(g1off + 18 * idx + 3);

            auto g1_xy_y = primBuffer.data(g1off + 18 * idx + 4);

            auto g1_xy_z = primBuffer.data(g1off + 18 * idx + 5);

            auto g1_xz_x = primBuffer.data(g1off + 18 * idx + 6);

            auto g1_xz_y = primBuffer.data(g1off + 18 * idx + 7);

            auto g1_xz_z = primBuffer.data(g1off + 18 * idx + 8);

            auto g1_yy_x = primBuffer.data(g1off + 18 * idx + 9);

            auto g1_yy_y = primBuffer.data(g1off + 18 * idx + 10);

            auto g1_yy_z = primBuffer.data(g1off + 18 * idx + 11);

            auto g1_yz_x = primBuffer.data(g1off + 18 * idx + 12);

            auto g1_yz_y = primBuffer.data(g1off + 18 * idx + 13);

            auto g1_yz_z = primBuffer.data(g1off + 18 * idx + 14);

            auto g1_zz_x = primBuffer.data(g1off + 18 * idx + 15);

            auto g1_zz_y = primBuffer.data(g1off + 18 * idx + 16);

            auto g1_zz_z = primBuffer.data(g1off + 18 * idx + 17);

            // set up pointers to (F|g(r,r')|P)^(m) integrals

            auto g_xxx_x = primBuffer.data(goff + 30 * idx);

            auto g_xxx_y = primBuffer.data(goff + 30 * idx + 1);

            auto g_xxx_z = primBuffer.data(goff + 30 * idx + 2);

            auto g_xxy_x = primBuffer.data(goff + 30 * idx + 3);

            auto g_xxy_y = primBuffer.data(goff + 30 * idx + 4);

            auto g_xxy_z = primBuffer.data(goff + 30 * idx + 5);

            auto g_xxz_x = primBuffer.data(goff + 30 * idx + 6);

            auto g_xxz_y = primBuffer.data(goff + 30 * idx + 7);

            auto g_xxz_z = primBuffer.data(goff + 30 * idx + 8);

            auto g_xyy_x = primBuffer.data(goff + 30 * idx + 9);

            auto g_xyy_y = primBuffer.data(goff + 30 * idx + 10);

            auto g_xyy_z = primBuffer.data(goff + 30 * idx + 11);

            auto g_xyz_x = primBuffer.data(goff + 30 * idx + 12);

            auto g_xyz_y = primBuffer.data(goff + 30 * idx + 13);

            auto g_xyz_z = primBuffer.data(goff + 30 * idx + 14);

            auto g_xzz_x = primBuffer.data(goff + 30 * idx + 15);

            auto g_xzz_y = primBuffer.data(goff + 30 * idx + 16);

            auto g_xzz_z = primBuffer.data(goff + 30 * idx + 17);

            auto g_yyy_x = primBuffer.data(goff + 30 * idx + 18);

            auto g_yyy_y = primBuffer.data(goff + 30 * idx + 19);

            auto g_yyy_z = primBuffer.data(goff + 30 * idx + 20);

            auto g_yyz_x = primBuffer.data(goff + 30 * idx + 21);

            auto g_yyz_y = primBuffer.data(goff + 30 * idx + 22);

            auto g_yyz_z = primBuffer.data(goff + 30 * idx + 23);

            auto g_yzz_x = primBuffer.data(goff + 30 * idx + 24);

            auto g_yzz_y = primBuffer.data(goff + 30 * idx + 25);

            auto g_yzz_z = primBuffer.data(goff + 30 * idx + 26);

            auto g_zzz_x = primBuffer.data(goff + 30 * idx + 27);

            auto g_zzz_y = primBuffer.data(goff + 30 * idx + 28);

            auto g_zzz_z = primBuffer.data(goff + 30 * idx + 29);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xx_0, gk_xy_0,\
                                         gk_xz_0, gk_yy_0, gk_yz_0, gk_zz_0, g20_x_x,\
                                         g20_x_y, g20_x_z, g20_y_x, g20_y_y, g20_y_z,\
                                         g20_z_x, g20_z_y, g20_z_z, g21_x_x, g21_x_y,\
                                         g21_x_z, g21_y_x, g21_y_y, g21_y_z, g21_z_x,\
                                         g21_z_y, g21_z_z, g1_xx_x, g1_xx_y, g1_xx_z,\
                                         g1_xy_x, g1_xy_y, g1_xy_z, g1_xz_x, g1_xz_y,\
                                         g1_xz_z, g1_yy_x, g1_yy_y, g1_yy_z, g1_yz_x,\
                                         g1_yz_y, g1_yz_z, g1_zz_x, g1_zz_y, g1_zz_z,\
                                         g_xxx_x, g_xxx_y, g_xxx_z, g_xxy_x, g_xxy_y,\
                                         g_xxy_z, g_xxz_x, g_xxz_y, g_xxz_z, g_xyy_x,\
                                         g_xyy_y, g_xyy_z, g_xyz_x, g_xyz_y, g_xyz_z,\
                                         g_xzz_x, g_xzz_y, g_xzz_z, g_yyy_x, g_yyy_y,\
                                         g_yyy_z, g_yyz_x, g_yyz_y, g_yyz_z, g_yzz_x,\
                                         g_yzz_y, g_yzz_z, g_zzz_x, g_zzz_y, g_zzz_z: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxx_x[k] = fr * g1_xx_x[k] + f2g * (2.0 * g20_x_x[k] - 2.0 * fgz * g21_x_x[k]) + f2t * gk_xx_0[k];

                g_xxx_y[k] = fr * g1_xx_y[k] + f2g * (2.0 * g20_x_y[k] - 2.0 * fgz * g21_x_y[k]);

                g_xxx_z[k] = fr * g1_xx_z[k] + f2g * (2.0 * g20_x_z[k] - 2.0 * fgz * g21_x_z[k]);

                g_xxy_x[k] = fr * g1_xy_x[k] + f2g * (g20_y_x[k] - fgz * g21_y_x[k]) + f2t * gk_xy_0[k];

                g_xxy_y[k] = fr * g1_xy_y[k] + f2g * (g20_y_y[k] - fgz * g21_y_y[k]);

                g_xxy_z[k] = fr * g1_xy_z[k] + f2g * (g20_y_z[k] - fgz * g21_y_z[k]);

                g_xxz_x[k] = fr * g1_xz_x[k] + f2g * (g20_z_x[k] - fgz * g21_z_x[k]) + f2t * gk_xz_0[k];

                g_xxz_y[k] = fr * g1_xz_y[k] + f2g * (g20_z_y[k] - fgz * g21_z_y[k]);

                g_xxz_z[k] = fr * g1_xz_z[k] + f2g * (g20_z_z[k] - fgz * g21_z_z[k]);

                g_xyy_x[k] = fr * g1_yy_x[k] + f2t * gk_yy_0[k];

                g_xyy_y[k] = fr * g1_yy_y[k];

                g_xyy_z[k] = fr * g1_yy_z[k];

                g_xyz_x[k] = fr * g1_yz_x[k] + f2t * gk_yz_0[k];

                g_xyz_y[k] = fr * g1_yz_y[k];

                g_xyz_z[k] = fr * g1_yz_z[k];

                g_xzz_x[k] = fr * g1_zz_x[k] + f2t * gk_zz_0[k];

                g_xzz_y[k] = fr * g1_zz_y[k];

                g_xzz_z[k] = fr * g1_zz_z[k];

                // leading y component

                fr = pay[k];

                g_yyy_x[k] = fr * g1_yy_x[k] + f2g * (2.0 * g20_y_x[k] - 2.0 * fgz * g21_y_x[k]);

                g_yyy_y[k] = fr * g1_yy_y[k] + f2g * (2.0 * g20_y_y[k] - 2.0 * fgz * g21_y_y[k]) + f2t * gk_yy_0[k];

                g_yyy_z[k] = fr * g1_yy_z[k] + f2g * (2.0 * g20_y_z[k] - 2.0 * fgz * g21_y_z[k]);

                g_yyz_x[k] = fr * g1_yz_x[k] + f2g * (g20_z_x[k] - fgz * g21_z_x[k]);

                g_yyz_y[k] = fr * g1_yz_y[k] + f2g * (g20_z_y[k] - fgz * g21_z_y[k]) + f2t * gk_yz_0[k];

                g_yyz_z[k] = fr * g1_yz_z[k] + f2g * (g20_z_z[k] - fgz * g21_z_z[k]);

                g_yzz_x[k] = fr * g1_zz_x[k];

                g_yzz_y[k] = fr * g1_zz_y[k] + f2t * gk_zz_0[k];

                g_yzz_z[k] = fr * g1_zz_z[k];

                // leading z component

                fr = paz[k];

                g_zzz_x[k] = fr * g1_zz_x[k] + f2g * (2.0 * g20_z_x[k] - 2.0 * fgz * g21_z_x[k]);

                g_zzz_y[k] = fr * g1_zz_y[k] + f2g * (2.0 * g20_z_y[k] - 2.0 * fgz * g21_z_y[k]);

                g_zzz_z[k] = fr * g1_zz_z[k] + f2g * (2.0 * g20_z_z[k] - 2.0 * fgz * g21_z_z[k]) + f2t * gk_zz_0[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForDF(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 2, 3);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|g(r,r')|D)^(m+1) integrals

            auto gk_x_xx = primBuffer.data(gkoff + 18 * idx);

            auto gk_x_xy = primBuffer.data(gkoff + 18 * idx + 1);

            auto gk_x_xz = primBuffer.data(gkoff + 18 * idx + 2);

            auto gk_x_yy = primBuffer.data(gkoff + 18 * idx + 3);

            auto gk_x_yz = primBuffer.data(gkoff + 18 * idx + 4);

            auto gk_x_zz = primBuffer.data(gkoff + 18 * idx + 5);

            auto gk_y_xx = primBuffer.data(gkoff + 18 * idx + 6);

            auto gk_y_xy = primBuffer.data(gkoff + 18 * idx + 7);

            auto gk_y_xz = primBuffer.data(gkoff + 18 * idx + 8);

            auto gk_y_yy = primBuffer.data(gkoff + 18 * idx + 9);

            auto gk_y_yz = primBuffer.data(gkoff + 18 * idx + 10);

            auto gk_y_zz = primBuffer.data(gkoff + 18 * idx + 11);

            auto gk_z_xx = primBuffer.data(gkoff + 18 * idx + 12);

            auto gk_z_xy = primBuffer.data(gkoff + 18 * idx + 13);

            auto gk_z_xz = primBuffer.data(gkoff + 18 * idx + 14);

            auto gk_z_yy = primBuffer.data(gkoff + 18 * idx + 15);

            auto gk_z_yz = primBuffer.data(gkoff + 18 * idx + 16);

            auto gk_z_zz = primBuffer.data(gkoff + 18 * idx + 17);

            // set up pointers to (S|g(r,r')|F)^(m) integrals

            auto g20_0_xxx = primBuffer.data(g20off + 10 * idx);

            auto g20_0_xxy = primBuffer.data(g20off + 10 * idx + 1);

            auto g20_0_xxz = primBuffer.data(g20off + 10 * idx + 2);

            auto g20_0_xyy = primBuffer.data(g20off + 10 * idx + 3);

            auto g20_0_xyz = primBuffer.data(g20off + 10 * idx + 4);

            auto g20_0_xzz = primBuffer.data(g20off + 10 * idx + 5);

            auto g20_0_yyy = primBuffer.data(g20off + 10 * idx + 6);

            auto g20_0_yyz = primBuffer.data(g20off + 10 * idx + 7);

            auto g20_0_yzz = primBuffer.data(g20off + 10 * idx + 8);

            auto g20_0_zzz = primBuffer.data(g20off + 10 * idx + 9);

            // set up pointers to (S|g(r,r')|F)^(m+1) integrals

            auto g21_0_xxx = primBuffer.data(g21off + 10 * idx);

            auto g21_0_xxy = primBuffer.data(g21off + 10 * idx + 1);

            auto g21_0_xxz = primBuffer.data(g21off + 10 * idx + 2);

            auto g21_0_xyy = primBuffer.data(g21off + 10 * idx + 3);

            auto g21_0_xyz = primBuffer.data(g21off + 10 * idx + 4);

            auto g21_0_xzz = primBuffer.data(g21off + 10 * idx + 5);

            auto g21_0_yyy = primBuffer.data(g21off + 10 * idx + 6);

            auto g21_0_yyz = primBuffer.data(g21off + 10 * idx + 7);

            auto g21_0_yzz = primBuffer.data(g21off + 10 * idx + 8);

            auto g21_0_zzz = primBuffer.data(g21off + 10 * idx + 9);

            // set up pointers to (P|g(r,r')|F)^(m+1) integrals

            auto g1_x_xxx = primBuffer.data(g1off + 30 * idx);

            auto g1_x_xxy = primBuffer.data(g1off + 30 * idx + 1);

            auto g1_x_xxz = primBuffer.data(g1off + 30 * idx + 2);

            auto g1_x_xyy = primBuffer.data(g1off + 30 * idx + 3);

            auto g1_x_xyz = primBuffer.data(g1off + 30 * idx + 4);

            auto g1_x_xzz = primBuffer.data(g1off + 30 * idx + 5);

            auto g1_x_yyy = primBuffer.data(g1off + 30 * idx + 6);

            auto g1_x_yyz = primBuffer.data(g1off + 30 * idx + 7);

            auto g1_x_yzz = primBuffer.data(g1off + 30 * idx + 8);

            auto g1_x_zzz = primBuffer.data(g1off + 30 * idx + 9);

            auto g1_y_xxx = primBuffer.data(g1off + 30 * idx + 10);

            auto g1_y_xxy = primBuffer.data(g1off + 30 * idx + 11);

            auto g1_y_xxz = primBuffer.data(g1off + 30 * idx + 12);

            auto g1_y_xyy = primBuffer.data(g1off + 30 * idx + 13);

            auto g1_y_xyz = primBuffer.data(g1off + 30 * idx + 14);

            auto g1_y_xzz = primBuffer.data(g1off + 30 * idx + 15);

            auto g1_y_yyy = primBuffer.data(g1off + 30 * idx + 16);

            auto g1_y_yyz = primBuffer.data(g1off + 30 * idx + 17);

            auto g1_y_yzz = primBuffer.data(g1off + 30 * idx + 18);

            auto g1_y_zzz = primBuffer.data(g1off + 30 * idx + 19);

            auto g1_z_xxx = primBuffer.data(g1off + 30 * idx + 20);

            auto g1_z_xxy = primBuffer.data(g1off + 30 * idx + 21);

            auto g1_z_xxz = primBuffer.data(g1off + 30 * idx + 22);

            auto g1_z_xyy = primBuffer.data(g1off + 30 * idx + 23);

            auto g1_z_xyz = primBuffer.data(g1off + 30 * idx + 24);

            auto g1_z_xzz = primBuffer.data(g1off + 30 * idx + 25);

            auto g1_z_yyy = primBuffer.data(g1off + 30 * idx + 26);

            auto g1_z_yyz = primBuffer.data(g1off + 30 * idx + 27);

            auto g1_z_yzz = primBuffer.data(g1off + 30 * idx + 28);

            auto g1_z_zzz = primBuffer.data(g1off + 30 * idx + 29);

            // set up pointers to (D|g(r,r')|F)^(m) integrals

            auto g_xx_xxx = primBuffer.data(goff + 60 * idx);

            auto g_xx_xxy = primBuffer.data(goff + 60 * idx + 1);

            auto g_xx_xxz = primBuffer.data(goff + 60 * idx + 2);

            auto g_xx_xyy = primBuffer.data(goff + 60 * idx + 3);

            auto g_xx_xyz = primBuffer.data(goff + 60 * idx + 4);

            auto g_xx_xzz = primBuffer.data(goff + 60 * idx + 5);

            auto g_xx_yyy = primBuffer.data(goff + 60 * idx + 6);

            auto g_xx_yyz = primBuffer.data(goff + 60 * idx + 7);

            auto g_xx_yzz = primBuffer.data(goff + 60 * idx + 8);

            auto g_xx_zzz = primBuffer.data(goff + 60 * idx + 9);

            auto g_xy_xxx = primBuffer.data(goff + 60 * idx + 10);

            auto g_xy_xxy = primBuffer.data(goff + 60 * idx + 11);

            auto g_xy_xxz = primBuffer.data(goff + 60 * idx + 12);

            auto g_xy_xyy = primBuffer.data(goff + 60 * idx + 13);

            auto g_xy_xyz = primBuffer.data(goff + 60 * idx + 14);

            auto g_xy_xzz = primBuffer.data(goff + 60 * idx + 15);

            auto g_xy_yyy = primBuffer.data(goff + 60 * idx + 16);

            auto g_xy_yyz = primBuffer.data(goff + 60 * idx + 17);

            auto g_xy_yzz = primBuffer.data(goff + 60 * idx + 18);

            auto g_xy_zzz = primBuffer.data(goff + 60 * idx + 19);

            auto g_xz_xxx = primBuffer.data(goff + 60 * idx + 20);

            auto g_xz_xxy = primBuffer.data(goff + 60 * idx + 21);

            auto g_xz_xxz = primBuffer.data(goff + 60 * idx + 22);

            auto g_xz_xyy = primBuffer.data(goff + 60 * idx + 23);

            auto g_xz_xyz = primBuffer.data(goff + 60 * idx + 24);

            auto g_xz_xzz = primBuffer.data(goff + 60 * idx + 25);

            auto g_xz_yyy = primBuffer.data(goff + 60 * idx + 26);

            auto g_xz_yyz = primBuffer.data(goff + 60 * idx + 27);

            auto g_xz_yzz = primBuffer.data(goff + 60 * idx + 28);

            auto g_xz_zzz = primBuffer.data(goff + 60 * idx + 29);

            auto g_yy_xxx = primBuffer.data(goff + 60 * idx + 30);

            auto g_yy_xxy = primBuffer.data(goff + 60 * idx + 31);

            auto g_yy_xxz = primBuffer.data(goff + 60 * idx + 32);

            auto g_yy_xyy = primBuffer.data(goff + 60 * idx + 33);

            auto g_yy_xyz = primBuffer.data(goff + 60 * idx + 34);

            auto g_yy_xzz = primBuffer.data(goff + 60 * idx + 35);

            auto g_yy_yyy = primBuffer.data(goff + 60 * idx + 36);

            auto g_yy_yyz = primBuffer.data(goff + 60 * idx + 37);

            auto g_yy_yzz = primBuffer.data(goff + 60 * idx + 38);

            auto g_yy_zzz = primBuffer.data(goff + 60 * idx + 39);

            auto g_yz_xxx = primBuffer.data(goff + 60 * idx + 40);

            auto g_yz_xxy = primBuffer.data(goff + 60 * idx + 41);

            auto g_yz_xxz = primBuffer.data(goff + 60 * idx + 42);

            auto g_yz_xyy = primBuffer.data(goff + 60 * idx + 43);

            auto g_yz_xyz = primBuffer.data(goff + 60 * idx + 44);

            auto g_yz_xzz = primBuffer.data(goff + 60 * idx + 45);

            auto g_yz_yyy = primBuffer.data(goff + 60 * idx + 46);

            auto g_yz_yyz = primBuffer.data(goff + 60 * idx + 47);

            auto g_yz_yzz = primBuffer.data(goff + 60 * idx + 48);

            auto g_yz_zzz = primBuffer.data(goff + 60 * idx + 49);

            auto g_zz_xxx = primBuffer.data(goff + 60 * idx + 50);

            auto g_zz_xxy = primBuffer.data(goff + 60 * idx + 51);

            auto g_zz_xxz = primBuffer.data(goff + 60 * idx + 52);

            auto g_zz_xyy = primBuffer.data(goff + 60 * idx + 53);

            auto g_zz_xyz = primBuffer.data(goff + 60 * idx + 54);

            auto g_zz_xzz = primBuffer.data(goff + 60 * idx + 55);

            auto g_zz_yyy = primBuffer.data(goff + 60 * idx + 56);

            auto g_zz_yyz = primBuffer.data(goff + 60 * idx + 57);

            auto g_zz_yzz = primBuffer.data(goff + 60 * idx + 58);

            auto g_zz_zzz = primBuffer.data(goff + 60 * idx + 59);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_x_xx, gk_x_xy,\
                                         gk_x_xz, gk_x_yy, gk_x_yz, gk_x_zz, gk_y_xx,\
                                         gk_y_xy, gk_y_xz, gk_y_yy, gk_y_yz, gk_y_zz,\
                                         gk_z_xx, gk_z_xy, gk_z_xz, gk_z_yy, gk_z_yz,\
                                         gk_z_zz, g20_0_xxx, g20_0_xxy, g20_0_xxz,\
                                         g20_0_xyy, g20_0_xyz, g20_0_xzz, g20_0_yyy,\
                                         g20_0_yyz, g20_0_yzz, g20_0_zzz, g21_0_xxx,\
                                         g21_0_xxy, g21_0_xxz, g21_0_xyy, g21_0_xyz,\
                                         g21_0_xzz, g21_0_yyy, g21_0_yyz, g21_0_yzz,\
                                         g21_0_zzz, g1_x_xxx, g1_x_xxy, g1_x_xxz,\
                                         g1_x_xyy, g1_x_xyz, g1_x_xzz, g1_x_yyy,\
                                         g1_x_yyz, g1_x_yzz, g1_x_zzz, g1_y_xxx,\
                                         g1_y_xxy, g1_y_xxz, g1_y_xyy, g1_y_xyz,\
                                         g1_y_xzz, g1_y_yyy, g1_y_yyz, g1_y_yzz,\
                                         g1_y_zzz, g1_z_xxx, g1_z_xxy, g1_z_xxz,\
                                         g1_z_xyy, g1_z_xyz, g1_z_xzz, g1_z_yyy,\
                                         g1_z_yyz, g1_z_yzz, g1_z_zzz, g_xx_xxx,\
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
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xx_xxx[k] = fr * g1_x_xxx[k] + f2g * (g20_0_xxx[k] - fgz * g21_0_xxx[k]) + 3.0 * f2t * gk_x_xx[k];

                g_xx_xxy[k] = fr * g1_x_xxy[k] + f2g * (g20_0_xxy[k] - fgz * g21_0_xxy[k]) + 2.0 * f2t * gk_x_xy[k];

                g_xx_xxz[k] = fr * g1_x_xxz[k] + f2g * (g20_0_xxz[k] - fgz * g21_0_xxz[k]) + 2.0 * f2t * gk_x_xz[k];

                g_xx_xyy[k] = fr * g1_x_xyy[k] + f2g * (g20_0_xyy[k] - fgz * g21_0_xyy[k]) + f2t * gk_x_yy[k];

                g_xx_xyz[k] = fr * g1_x_xyz[k] + f2g * (g20_0_xyz[k] - fgz * g21_0_xyz[k]) + f2t * gk_x_yz[k];

                g_xx_xzz[k] = fr * g1_x_xzz[k] + f2g * (g20_0_xzz[k] - fgz * g21_0_xzz[k]) + f2t * gk_x_zz[k];

                g_xx_yyy[k] = fr * g1_x_yyy[k] + f2g * (g20_0_yyy[k] - fgz * g21_0_yyy[k]);

                g_xx_yyz[k] = fr * g1_x_yyz[k] + f2g * (g20_0_yyz[k] - fgz * g21_0_yyz[k]);

                g_xx_yzz[k] = fr * g1_x_yzz[k] + f2g * (g20_0_yzz[k] - fgz * g21_0_yzz[k]);

                g_xx_zzz[k] = fr * g1_x_zzz[k] + f2g * (g20_0_zzz[k] - fgz * g21_0_zzz[k]);

                g_xy_xxx[k] = fr * g1_y_xxx[k] + 3.0 * f2t * gk_y_xx[k];

                g_xy_xxy[k] = fr * g1_y_xxy[k] + 2.0 * f2t * gk_y_xy[k];

                g_xy_xxz[k] = fr * g1_y_xxz[k] + 2.0 * f2t * gk_y_xz[k];

                g_xy_xyy[k] = fr * g1_y_xyy[k] + f2t * gk_y_yy[k];

                g_xy_xyz[k] = fr * g1_y_xyz[k] + f2t * gk_y_yz[k];

                g_xy_xzz[k] = fr * g1_y_xzz[k] + f2t * gk_y_zz[k];

                g_xy_yyy[k] = fr * g1_y_yyy[k];

                g_xy_yyz[k] = fr * g1_y_yyz[k];

                g_xy_yzz[k] = fr * g1_y_yzz[k];

                g_xy_zzz[k] = fr * g1_y_zzz[k];

                g_xz_xxx[k] = fr * g1_z_xxx[k] + 3.0 * f2t * gk_z_xx[k];

                g_xz_xxy[k] = fr * g1_z_xxy[k] + 2.0 * f2t * gk_z_xy[k];

                g_xz_xxz[k] = fr * g1_z_xxz[k] + 2.0 * f2t * gk_z_xz[k];

                g_xz_xyy[k] = fr * g1_z_xyy[k] + f2t * gk_z_yy[k];

                g_xz_xyz[k] = fr * g1_z_xyz[k] + f2t * gk_z_yz[k];

                g_xz_xzz[k] = fr * g1_z_xzz[k] + f2t * gk_z_zz[k];

                g_xz_yyy[k] = fr * g1_z_yyy[k];

                g_xz_yyz[k] = fr * g1_z_yyz[k];

                g_xz_yzz[k] = fr * g1_z_yzz[k];

                g_xz_zzz[k] = fr * g1_z_zzz[k];

                // leading y component

                fr = pay[k];

                g_yy_xxx[k] = fr * g1_y_xxx[k] + f2g * (g20_0_xxx[k] - fgz * g21_0_xxx[k]);

                g_yy_xxy[k] = fr * g1_y_xxy[k] + f2g * (g20_0_xxy[k] - fgz * g21_0_xxy[k]) + f2t * gk_y_xx[k];

                g_yy_xxz[k] = fr * g1_y_xxz[k] + f2g * (g20_0_xxz[k] - fgz * g21_0_xxz[k]);

                g_yy_xyy[k] = fr * g1_y_xyy[k] + f2g * (g20_0_xyy[k] - fgz * g21_0_xyy[k]) + 2.0 * f2t * gk_y_xy[k];

                g_yy_xyz[k] = fr * g1_y_xyz[k] + f2g * (g20_0_xyz[k] - fgz * g21_0_xyz[k]) + f2t * gk_y_xz[k];

                g_yy_xzz[k] = fr * g1_y_xzz[k] + f2g * (g20_0_xzz[k] - fgz * g21_0_xzz[k]);

                g_yy_yyy[k] = fr * g1_y_yyy[k] + f2g * (g20_0_yyy[k] - fgz * g21_0_yyy[k]) + 3.0 * f2t * gk_y_yy[k];

                g_yy_yyz[k] = fr * g1_y_yyz[k] + f2g * (g20_0_yyz[k] - fgz * g21_0_yyz[k]) + 2.0 * f2t * gk_y_yz[k];

                g_yy_yzz[k] = fr * g1_y_yzz[k] + f2g * (g20_0_yzz[k] - fgz * g21_0_yzz[k]) + f2t * gk_y_zz[k];

                g_yy_zzz[k] = fr * g1_y_zzz[k] + f2g * (g20_0_zzz[k] - fgz * g21_0_zzz[k]);

                g_yz_xxx[k] = fr * g1_z_xxx[k];

                g_yz_xxy[k] = fr * g1_z_xxy[k] + f2t * gk_z_xx[k];

                g_yz_xxz[k] = fr * g1_z_xxz[k];

                g_yz_xyy[k] = fr * g1_z_xyy[k] + 2.0 * f2t * gk_z_xy[k];

                g_yz_xyz[k] = fr * g1_z_xyz[k] + f2t * gk_z_xz[k];

                g_yz_xzz[k] = fr * g1_z_xzz[k];

                g_yz_yyy[k] = fr * g1_z_yyy[k] + 3.0 * f2t * gk_z_yy[k];

                g_yz_yyz[k] = fr * g1_z_yyz[k] + 2.0 * f2t * gk_z_yz[k];

                g_yz_yzz[k] = fr * g1_z_yzz[k] + f2t * gk_z_zz[k];

                g_yz_zzz[k] = fr * g1_z_zzz[k];

                // leading z component

                fr = paz[k];

                g_zz_xxx[k] = fr * g1_z_xxx[k] + f2g * (g20_0_xxx[k] - fgz * g21_0_xxx[k]);

                g_zz_xxy[k] = fr * g1_z_xxy[k] + f2g * (g20_0_xxy[k] - fgz * g21_0_xxy[k]);

                g_zz_xxz[k] = fr * g1_z_xxz[k] + f2g * (g20_0_xxz[k] - fgz * g21_0_xxz[k]) + f2t * gk_z_xx[k];

                g_zz_xyy[k] = fr * g1_z_xyy[k] + f2g * (g20_0_xyy[k] - fgz * g21_0_xyy[k]);

                g_zz_xyz[k] = fr * g1_z_xyz[k] + f2g * (g20_0_xyz[k] - fgz * g21_0_xyz[k]) + f2t * gk_z_xy[k];

                g_zz_xzz[k] = fr * g1_z_xzz[k] + f2g * (g20_0_xzz[k] - fgz * g21_0_xzz[k]) + 2.0 * f2t * gk_z_xz[k];

                g_zz_yyy[k] = fr * g1_z_yyy[k] + f2g * (g20_0_yyy[k] - fgz * g21_0_yyy[k]);

                g_zz_yyz[k] = fr * g1_z_yyz[k] + f2g * (g20_0_yyz[k] - fgz * g21_0_yyz[k]) + f2t * gk_z_yy[k];

                g_zz_yzz[k] = fr * g1_z_yzz[k] + f2g * (g20_0_yzz[k] - fgz * g21_0_yzz[k]) + 2.0 * f2t * gk_z_yz[k];

                g_zz_zzz[k] = fr * g1_z_zzz[k] + f2g * (g20_0_zzz[k] - fgz * g21_0_zzz[k]) + 3.0 * f2t * gk_z_zz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForFD(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 3, 2);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 2, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|g(r,r')|P)^(m+1) integrals

            auto gk_xx_x = primBuffer.data(gkoff + 18 * idx);

            auto gk_xx_y = primBuffer.data(gkoff + 18 * idx + 1);

            auto gk_xx_z = primBuffer.data(gkoff + 18 * idx + 2);

            auto gk_xy_x = primBuffer.data(gkoff + 18 * idx + 3);

            auto gk_xy_y = primBuffer.data(gkoff + 18 * idx + 4);

            auto gk_xy_z = primBuffer.data(gkoff + 18 * idx + 5);

            auto gk_xz_x = primBuffer.data(gkoff + 18 * idx + 6);

            auto gk_xz_y = primBuffer.data(gkoff + 18 * idx + 7);

            auto gk_xz_z = primBuffer.data(gkoff + 18 * idx + 8);

            auto gk_yy_x = primBuffer.data(gkoff + 18 * idx + 9);

            auto gk_yy_y = primBuffer.data(gkoff + 18 * idx + 10);

            auto gk_yy_z = primBuffer.data(gkoff + 18 * idx + 11);

            auto gk_yz_x = primBuffer.data(gkoff + 18 * idx + 12);

            auto gk_yz_y = primBuffer.data(gkoff + 18 * idx + 13);

            auto gk_yz_z = primBuffer.data(gkoff + 18 * idx + 14);

            auto gk_zz_x = primBuffer.data(gkoff + 18 * idx + 15);

            auto gk_zz_y = primBuffer.data(gkoff + 18 * idx + 16);

            auto gk_zz_z = primBuffer.data(gkoff + 18 * idx + 17);

            // set up pointers to (P|g(r,r')|D)^(m) integrals

            auto g20_x_xx = primBuffer.data(g20off + 18 * idx);

            auto g20_x_xy = primBuffer.data(g20off + 18 * idx + 1);

            auto g20_x_xz = primBuffer.data(g20off + 18 * idx + 2);

            auto g20_x_yy = primBuffer.data(g20off + 18 * idx + 3);

            auto g20_x_yz = primBuffer.data(g20off + 18 * idx + 4);

            auto g20_x_zz = primBuffer.data(g20off + 18 * idx + 5);

            auto g20_y_xx = primBuffer.data(g20off + 18 * idx + 6);

            auto g20_y_xy = primBuffer.data(g20off + 18 * idx + 7);

            auto g20_y_xz = primBuffer.data(g20off + 18 * idx + 8);

            auto g20_y_yy = primBuffer.data(g20off + 18 * idx + 9);

            auto g20_y_yz = primBuffer.data(g20off + 18 * idx + 10);

            auto g20_y_zz = primBuffer.data(g20off + 18 * idx + 11);

            auto g20_z_xx = primBuffer.data(g20off + 18 * idx + 12);

            auto g20_z_xy = primBuffer.data(g20off + 18 * idx + 13);

            auto g20_z_xz = primBuffer.data(g20off + 18 * idx + 14);

            auto g20_z_yy = primBuffer.data(g20off + 18 * idx + 15);

            auto g20_z_yz = primBuffer.data(g20off + 18 * idx + 16);

            auto g20_z_zz = primBuffer.data(g20off + 18 * idx + 17);

            // set up pointers to (P|g(r,r')|D)^(m+1) integrals

            auto g21_x_xx = primBuffer.data(g21off + 18 * idx);

            auto g21_x_xy = primBuffer.data(g21off + 18 * idx + 1);

            auto g21_x_xz = primBuffer.data(g21off + 18 * idx + 2);

            auto g21_x_yy = primBuffer.data(g21off + 18 * idx + 3);

            auto g21_x_yz = primBuffer.data(g21off + 18 * idx + 4);

            auto g21_x_zz = primBuffer.data(g21off + 18 * idx + 5);

            auto g21_y_xx = primBuffer.data(g21off + 18 * idx + 6);

            auto g21_y_xy = primBuffer.data(g21off + 18 * idx + 7);

            auto g21_y_xz = primBuffer.data(g21off + 18 * idx + 8);

            auto g21_y_yy = primBuffer.data(g21off + 18 * idx + 9);

            auto g21_y_yz = primBuffer.data(g21off + 18 * idx + 10);

            auto g21_y_zz = primBuffer.data(g21off + 18 * idx + 11);

            auto g21_z_xx = primBuffer.data(g21off + 18 * idx + 12);

            auto g21_z_xy = primBuffer.data(g21off + 18 * idx + 13);

            auto g21_z_xz = primBuffer.data(g21off + 18 * idx + 14);

            auto g21_z_yy = primBuffer.data(g21off + 18 * idx + 15);

            auto g21_z_yz = primBuffer.data(g21off + 18 * idx + 16);

            auto g21_z_zz = primBuffer.data(g21off + 18 * idx + 17);

            // set up pointers to (D|g(r,r')|D)^(m+1) integrals

            auto g1_xx_xx = primBuffer.data(g1off + 36 * idx);

            auto g1_xx_xy = primBuffer.data(g1off + 36 * idx + 1);

            auto g1_xx_xz = primBuffer.data(g1off + 36 * idx + 2);

            auto g1_xx_yy = primBuffer.data(g1off + 36 * idx + 3);

            auto g1_xx_yz = primBuffer.data(g1off + 36 * idx + 4);

            auto g1_xx_zz = primBuffer.data(g1off + 36 * idx + 5);

            auto g1_xy_xx = primBuffer.data(g1off + 36 * idx + 6);

            auto g1_xy_xy = primBuffer.data(g1off + 36 * idx + 7);

            auto g1_xy_xz = primBuffer.data(g1off + 36 * idx + 8);

            auto g1_xy_yy = primBuffer.data(g1off + 36 * idx + 9);

            auto g1_xy_yz = primBuffer.data(g1off + 36 * idx + 10);

            auto g1_xy_zz = primBuffer.data(g1off + 36 * idx + 11);

            auto g1_xz_xx = primBuffer.data(g1off + 36 * idx + 12);

            auto g1_xz_xy = primBuffer.data(g1off + 36 * idx + 13);

            auto g1_xz_xz = primBuffer.data(g1off + 36 * idx + 14);

            auto g1_xz_yy = primBuffer.data(g1off + 36 * idx + 15);

            auto g1_xz_yz = primBuffer.data(g1off + 36 * idx + 16);

            auto g1_xz_zz = primBuffer.data(g1off + 36 * idx + 17);

            auto g1_yy_xx = primBuffer.data(g1off + 36 * idx + 18);

            auto g1_yy_xy = primBuffer.data(g1off + 36 * idx + 19);

            auto g1_yy_xz = primBuffer.data(g1off + 36 * idx + 20);

            auto g1_yy_yy = primBuffer.data(g1off + 36 * idx + 21);

            auto g1_yy_yz = primBuffer.data(g1off + 36 * idx + 22);

            auto g1_yy_zz = primBuffer.data(g1off + 36 * idx + 23);

            auto g1_yz_xx = primBuffer.data(g1off + 36 * idx + 24);

            auto g1_yz_xy = primBuffer.data(g1off + 36 * idx + 25);

            auto g1_yz_xz = primBuffer.data(g1off + 36 * idx + 26);

            auto g1_yz_yy = primBuffer.data(g1off + 36 * idx + 27);

            auto g1_yz_yz = primBuffer.data(g1off + 36 * idx + 28);

            auto g1_yz_zz = primBuffer.data(g1off + 36 * idx + 29);

            auto g1_zz_xx = primBuffer.data(g1off + 36 * idx + 30);

            auto g1_zz_xy = primBuffer.data(g1off + 36 * idx + 31);

            auto g1_zz_xz = primBuffer.data(g1off + 36 * idx + 32);

            auto g1_zz_yy = primBuffer.data(g1off + 36 * idx + 33);

            auto g1_zz_yz = primBuffer.data(g1off + 36 * idx + 34);

            auto g1_zz_zz = primBuffer.data(g1off + 36 * idx + 35);

            // set up pointers to (F|g(r,r')|D)^(m) integrals

            auto g_xxx_xx = primBuffer.data(goff + 60 * idx);

            auto g_xxx_xy = primBuffer.data(goff + 60 * idx + 1);

            auto g_xxx_xz = primBuffer.data(goff + 60 * idx + 2);

            auto g_xxx_yy = primBuffer.data(goff + 60 * idx + 3);

            auto g_xxx_yz = primBuffer.data(goff + 60 * idx + 4);

            auto g_xxx_zz = primBuffer.data(goff + 60 * idx + 5);

            auto g_xxy_xx = primBuffer.data(goff + 60 * idx + 6);

            auto g_xxy_xy = primBuffer.data(goff + 60 * idx + 7);

            auto g_xxy_xz = primBuffer.data(goff + 60 * idx + 8);

            auto g_xxy_yy = primBuffer.data(goff + 60 * idx + 9);

            auto g_xxy_yz = primBuffer.data(goff + 60 * idx + 10);

            auto g_xxy_zz = primBuffer.data(goff + 60 * idx + 11);

            auto g_xxz_xx = primBuffer.data(goff + 60 * idx + 12);

            auto g_xxz_xy = primBuffer.data(goff + 60 * idx + 13);

            auto g_xxz_xz = primBuffer.data(goff + 60 * idx + 14);

            auto g_xxz_yy = primBuffer.data(goff + 60 * idx + 15);

            auto g_xxz_yz = primBuffer.data(goff + 60 * idx + 16);

            auto g_xxz_zz = primBuffer.data(goff + 60 * idx + 17);

            auto g_xyy_xx = primBuffer.data(goff + 60 * idx + 18);

            auto g_xyy_xy = primBuffer.data(goff + 60 * idx + 19);

            auto g_xyy_xz = primBuffer.data(goff + 60 * idx + 20);

            auto g_xyy_yy = primBuffer.data(goff + 60 * idx + 21);

            auto g_xyy_yz = primBuffer.data(goff + 60 * idx + 22);

            auto g_xyy_zz = primBuffer.data(goff + 60 * idx + 23);

            auto g_xyz_xx = primBuffer.data(goff + 60 * idx + 24);

            auto g_xyz_xy = primBuffer.data(goff + 60 * idx + 25);

            auto g_xyz_xz = primBuffer.data(goff + 60 * idx + 26);

            auto g_xyz_yy = primBuffer.data(goff + 60 * idx + 27);

            auto g_xyz_yz = primBuffer.data(goff + 60 * idx + 28);

            auto g_xyz_zz = primBuffer.data(goff + 60 * idx + 29);

            auto g_xzz_xx = primBuffer.data(goff + 60 * idx + 30);

            auto g_xzz_xy = primBuffer.data(goff + 60 * idx + 31);

            auto g_xzz_xz = primBuffer.data(goff + 60 * idx + 32);

            auto g_xzz_yy = primBuffer.data(goff + 60 * idx + 33);

            auto g_xzz_yz = primBuffer.data(goff + 60 * idx + 34);

            auto g_xzz_zz = primBuffer.data(goff + 60 * idx + 35);

            auto g_yyy_xx = primBuffer.data(goff + 60 * idx + 36);

            auto g_yyy_xy = primBuffer.data(goff + 60 * idx + 37);

            auto g_yyy_xz = primBuffer.data(goff + 60 * idx + 38);

            auto g_yyy_yy = primBuffer.data(goff + 60 * idx + 39);

            auto g_yyy_yz = primBuffer.data(goff + 60 * idx + 40);

            auto g_yyy_zz = primBuffer.data(goff + 60 * idx + 41);

            auto g_yyz_xx = primBuffer.data(goff + 60 * idx + 42);

            auto g_yyz_xy = primBuffer.data(goff + 60 * idx + 43);

            auto g_yyz_xz = primBuffer.data(goff + 60 * idx + 44);

            auto g_yyz_yy = primBuffer.data(goff + 60 * idx + 45);

            auto g_yyz_yz = primBuffer.data(goff + 60 * idx + 46);

            auto g_yyz_zz = primBuffer.data(goff + 60 * idx + 47);

            auto g_yzz_xx = primBuffer.data(goff + 60 * idx + 48);

            auto g_yzz_xy = primBuffer.data(goff + 60 * idx + 49);

            auto g_yzz_xz = primBuffer.data(goff + 60 * idx + 50);

            auto g_yzz_yy = primBuffer.data(goff + 60 * idx + 51);

            auto g_yzz_yz = primBuffer.data(goff + 60 * idx + 52);

            auto g_yzz_zz = primBuffer.data(goff + 60 * idx + 53);

            auto g_zzz_xx = primBuffer.data(goff + 60 * idx + 54);

            auto g_zzz_xy = primBuffer.data(goff + 60 * idx + 55);

            auto g_zzz_xz = primBuffer.data(goff + 60 * idx + 56);

            auto g_zzz_yy = primBuffer.data(goff + 60 * idx + 57);

            auto g_zzz_yz = primBuffer.data(goff + 60 * idx + 58);

            auto g_zzz_zz = primBuffer.data(goff + 60 * idx + 59);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xx_x, gk_xx_y,\
                                         gk_xx_z, gk_xy_x, gk_xy_y, gk_xy_z, gk_xz_x,\
                                         gk_xz_y, gk_xz_z, gk_yy_x, gk_yy_y, gk_yy_z,\
                                         gk_yz_x, gk_yz_y, gk_yz_z, gk_zz_x, gk_zz_y,\
                                         gk_zz_z, g20_x_xx, g20_x_xy, g20_x_xz,\
                                         g20_x_yy, g20_x_yz, g20_x_zz, g20_y_xx,\
                                         g20_y_xy, g20_y_xz, g20_y_yy, g20_y_yz,\
                                         g20_y_zz, g20_z_xx, g20_z_xy, g20_z_xz,\
                                         g20_z_yy, g20_z_yz, g20_z_zz, g21_x_xx,\
                                         g21_x_xy, g21_x_xz, g21_x_yy, g21_x_yz,\
                                         g21_x_zz, g21_y_xx, g21_y_xy, g21_y_xz,\
                                         g21_y_yy, g21_y_yz, g21_y_zz, g21_z_xx,\
                                         g21_z_xy, g21_z_xz, g21_z_yy, g21_z_yz,\
                                         g21_z_zz, g1_xx_xx, g1_xx_xy, g1_xx_xz,\
                                         g1_xx_yy, g1_xx_yz, g1_xx_zz, g1_xy_xx,\
                                         g1_xy_xy, g1_xy_xz, g1_xy_yy, g1_xy_yz,\
                                         g1_xy_zz, g1_xz_xx, g1_xz_xy, g1_xz_xz,\
                                         g1_xz_yy, g1_xz_yz, g1_xz_zz, g1_yy_xx,\
                                         g1_yy_xy, g1_yy_xz, g1_yy_yy, g1_yy_yz,\
                                         g1_yy_zz, g1_yz_xx, g1_yz_xy, g1_yz_xz,\
                                         g1_yz_yy, g1_yz_yz, g1_yz_zz, g1_zz_xx,\
                                         g1_zz_xy, g1_zz_xz, g1_zz_yy, g1_zz_yz,\
                                         g1_zz_zz, g_xxx_xx, g_xxx_xy, g_xxx_xz,\
                                         g_xxx_yy, g_xxx_yz, g_xxx_zz, g_xxy_xx,\
                                         g_xxy_xy, g_xxy_xz, g_xxy_yy, g_xxy_yz,\
                                         g_xxy_zz, g_xxz_xx, g_xxz_xy, g_xxz_xz,\
                                         g_xxz_yy, g_xxz_yz, g_xxz_zz, g_xyy_xx,\
                                         g_xyy_xy, g_xyy_xz, g_xyy_yy, g_xyy_yz,\
                                         g_xyy_zz, g_xyz_xx, g_xyz_xy, g_xyz_xz,\
                                         g_xyz_yy, g_xyz_yz, g_xyz_zz, g_xzz_xx,\
                                         g_xzz_xy, g_xzz_xz, g_xzz_yy, g_xzz_yz,\
                                         g_xzz_zz, g_yyy_xx, g_yyy_xy, g_yyy_xz,\
                                         g_yyy_yy, g_yyy_yz, g_yyy_zz, g_yyz_xx,\
                                         g_yyz_xy, g_yyz_xz, g_yyz_yy, g_yyz_yz,\
                                         g_yyz_zz, g_yzz_xx, g_yzz_xy, g_yzz_xz,\
                                         g_yzz_yy, g_yzz_yz, g_yzz_zz, g_zzz_xx,\
                                         g_zzz_xy, g_zzz_xz, g_zzz_yy, g_zzz_yz,\
                                         g_zzz_zz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxx_xx[k] = fr * g1_xx_xx[k] + f2g * (2.0 * g20_x_xx[k] - 2.0 * fgz * g21_x_xx[k]) + 2.0 * f2t * gk_xx_x[k];

                g_xxx_xy[k] = fr * g1_xx_xy[k] + f2g * (2.0 * g20_x_xy[k] - 2.0 * fgz * g21_x_xy[k]) + f2t * gk_xx_y[k];

                g_xxx_xz[k] = fr * g1_xx_xz[k] + f2g * (2.0 * g20_x_xz[k] - 2.0 * fgz * g21_x_xz[k]) + f2t * gk_xx_z[k];

                g_xxx_yy[k] = fr * g1_xx_yy[k] + f2g * (2.0 * g20_x_yy[k] - 2.0 * fgz * g21_x_yy[k]);

                g_xxx_yz[k] = fr * g1_xx_yz[k] + f2g * (2.0 * g20_x_yz[k] - 2.0 * fgz * g21_x_yz[k]);

                g_xxx_zz[k] = fr * g1_xx_zz[k] + f2g * (2.0 * g20_x_zz[k] - 2.0 * fgz * g21_x_zz[k]);

                g_xxy_xx[k] = fr * g1_xy_xx[k] + f2g * (g20_y_xx[k] - fgz * g21_y_xx[k]) + 2.0 * f2t * gk_xy_x[k];

                g_xxy_xy[k] = fr * g1_xy_xy[k] + f2g * (g20_y_xy[k] - fgz * g21_y_xy[k]) + f2t * gk_xy_y[k];

                g_xxy_xz[k] = fr * g1_xy_xz[k] + f2g * (g20_y_xz[k] - fgz * g21_y_xz[k]) + f2t * gk_xy_z[k];

                g_xxy_yy[k] = fr * g1_xy_yy[k] + f2g * (g20_y_yy[k] - fgz * g21_y_yy[k]);

                g_xxy_yz[k] = fr * g1_xy_yz[k] + f2g * (g20_y_yz[k] - fgz * g21_y_yz[k]);

                g_xxy_zz[k] = fr * g1_xy_zz[k] + f2g * (g20_y_zz[k] - fgz * g21_y_zz[k]);

                g_xxz_xx[k] = fr * g1_xz_xx[k] + f2g * (g20_z_xx[k] - fgz * g21_z_xx[k]) + 2.0 * f2t * gk_xz_x[k];

                g_xxz_xy[k] = fr * g1_xz_xy[k] + f2g * (g20_z_xy[k] - fgz * g21_z_xy[k]) + f2t * gk_xz_y[k];

                g_xxz_xz[k] = fr * g1_xz_xz[k] + f2g * (g20_z_xz[k] - fgz * g21_z_xz[k]) + f2t * gk_xz_z[k];

                g_xxz_yy[k] = fr * g1_xz_yy[k] + f2g * (g20_z_yy[k] - fgz * g21_z_yy[k]);

                g_xxz_yz[k] = fr * g1_xz_yz[k] + f2g * (g20_z_yz[k] - fgz * g21_z_yz[k]);

                g_xxz_zz[k] = fr * g1_xz_zz[k] + f2g * (g20_z_zz[k] - fgz * g21_z_zz[k]);

                g_xyy_xx[k] = fr * g1_yy_xx[k] + 2.0 * f2t * gk_yy_x[k];

                g_xyy_xy[k] = fr * g1_yy_xy[k] + f2t * gk_yy_y[k];

                g_xyy_xz[k] = fr * g1_yy_xz[k] + f2t * gk_yy_z[k];

                g_xyy_yy[k] = fr * g1_yy_yy[k];

                g_xyy_yz[k] = fr * g1_yy_yz[k];

                g_xyy_zz[k] = fr * g1_yy_zz[k];

                g_xyz_xx[k] = fr * g1_yz_xx[k] + 2.0 * f2t * gk_yz_x[k];

                g_xyz_xy[k] = fr * g1_yz_xy[k] + f2t * gk_yz_y[k];

                g_xyz_xz[k] = fr * g1_yz_xz[k] + f2t * gk_yz_z[k];

                g_xyz_yy[k] = fr * g1_yz_yy[k];

                g_xyz_yz[k] = fr * g1_yz_yz[k];

                g_xyz_zz[k] = fr * g1_yz_zz[k];

                g_xzz_xx[k] = fr * g1_zz_xx[k] + 2.0 * f2t * gk_zz_x[k];

                g_xzz_xy[k] = fr * g1_zz_xy[k] + f2t * gk_zz_y[k];

                g_xzz_xz[k] = fr * g1_zz_xz[k] + f2t * gk_zz_z[k];

                g_xzz_yy[k] = fr * g1_zz_yy[k];

                g_xzz_yz[k] = fr * g1_zz_yz[k];

                g_xzz_zz[k] = fr * g1_zz_zz[k];

                // leading y component

                fr = pay[k];

                g_yyy_xx[k] = fr * g1_yy_xx[k] + f2g * (2.0 * g20_y_xx[k] - 2.0 * fgz * g21_y_xx[k]);

                g_yyy_xy[k] = fr * g1_yy_xy[k] + f2g * (2.0 * g20_y_xy[k] - 2.0 * fgz * g21_y_xy[k]) + f2t * gk_yy_x[k];

                g_yyy_xz[k] = fr * g1_yy_xz[k] + f2g * (2.0 * g20_y_xz[k] - 2.0 * fgz * g21_y_xz[k]);

                g_yyy_yy[k] = fr * g1_yy_yy[k] + f2g * (2.0 * g20_y_yy[k] - 2.0 * fgz * g21_y_yy[k]) + 2.0 * f2t * gk_yy_y[k];

                g_yyy_yz[k] = fr * g1_yy_yz[k] + f2g * (2.0 * g20_y_yz[k] - 2.0 * fgz * g21_y_yz[k]) + f2t * gk_yy_z[k];

                g_yyy_zz[k] = fr * g1_yy_zz[k] + f2g * (2.0 * g20_y_zz[k] - 2.0 * fgz * g21_y_zz[k]);

                g_yyz_xx[k] = fr * g1_yz_xx[k] + f2g * (g20_z_xx[k] - fgz * g21_z_xx[k]);

                g_yyz_xy[k] = fr * g1_yz_xy[k] + f2g * (g20_z_xy[k] - fgz * g21_z_xy[k]) + f2t * gk_yz_x[k];

                g_yyz_xz[k] = fr * g1_yz_xz[k] + f2g * (g20_z_xz[k] - fgz * g21_z_xz[k]);

                g_yyz_yy[k] = fr * g1_yz_yy[k] + f2g * (g20_z_yy[k] - fgz * g21_z_yy[k]) + 2.0 * f2t * gk_yz_y[k];

                g_yyz_yz[k] = fr * g1_yz_yz[k] + f2g * (g20_z_yz[k] - fgz * g21_z_yz[k]) + f2t * gk_yz_z[k];

                g_yyz_zz[k] = fr * g1_yz_zz[k] + f2g * (g20_z_zz[k] - fgz * g21_z_zz[k]);

                g_yzz_xx[k] = fr * g1_zz_xx[k];

                g_yzz_xy[k] = fr * g1_zz_xy[k] + f2t * gk_zz_x[k];

                g_yzz_xz[k] = fr * g1_zz_xz[k];

                g_yzz_yy[k] = fr * g1_zz_yy[k] + 2.0 * f2t * gk_zz_y[k];

                g_yzz_yz[k] = fr * g1_zz_yz[k] + f2t * gk_zz_z[k];

                g_yzz_zz[k] = fr * g1_zz_zz[k];

                // leading z component

                fr = paz[k];

                g_zzz_xx[k] = fr * g1_zz_xx[k] + f2g * (2.0 * g20_z_xx[k] - 2.0 * fgz * g21_z_xx[k]);

                g_zzz_xy[k] = fr * g1_zz_xy[k] + f2g * (2.0 * g20_z_xy[k] - 2.0 * fgz * g21_z_xy[k]);

                g_zzz_xz[k] = fr * g1_zz_xz[k] + f2g * (2.0 * g20_z_xz[k] - 2.0 * fgz * g21_z_xz[k]) + f2t * gk_zz_x[k];

                g_zzz_yy[k] = fr * g1_zz_yy[k] + f2g * (2.0 * g20_z_yy[k] - 2.0 * fgz * g21_z_yy[k]);

                g_zzz_yz[k] = fr * g1_zz_yz[k] + f2g * (2.0 * g20_z_yz[k] - 2.0 * fgz * g21_z_yz[k]) + f2t * gk_zz_y[k];

                g_zzz_zz[k] = fr * g1_zz_zz[k] + f2g * (2.0 * g20_z_zz[k] - 2.0 * fgz * g21_z_zz[k]) + 2.0 * f2t * gk_zz_z[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForFF(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 3, 3);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|g(r,r')|D)^(m+1) integrals

            auto gk_xx_xx = primBuffer.data(gkoff + 36 * idx);

            auto gk_xx_xy = primBuffer.data(gkoff + 36 * idx + 1);

            auto gk_xx_xz = primBuffer.data(gkoff + 36 * idx + 2);

            auto gk_xx_yy = primBuffer.data(gkoff + 36 * idx + 3);

            auto gk_xx_yz = primBuffer.data(gkoff + 36 * idx + 4);

            auto gk_xx_zz = primBuffer.data(gkoff + 36 * idx + 5);

            auto gk_xy_xx = primBuffer.data(gkoff + 36 * idx + 6);

            auto gk_xy_xy = primBuffer.data(gkoff + 36 * idx + 7);

            auto gk_xy_xz = primBuffer.data(gkoff + 36 * idx + 8);

            auto gk_xy_yy = primBuffer.data(gkoff + 36 * idx + 9);

            auto gk_xy_yz = primBuffer.data(gkoff + 36 * idx + 10);

            auto gk_xy_zz = primBuffer.data(gkoff + 36 * idx + 11);

            auto gk_xz_xx = primBuffer.data(gkoff + 36 * idx + 12);

            auto gk_xz_xy = primBuffer.data(gkoff + 36 * idx + 13);

            auto gk_xz_xz = primBuffer.data(gkoff + 36 * idx + 14);

            auto gk_xz_yy = primBuffer.data(gkoff + 36 * idx + 15);

            auto gk_xz_yz = primBuffer.data(gkoff + 36 * idx + 16);

            auto gk_xz_zz = primBuffer.data(gkoff + 36 * idx + 17);

            auto gk_yy_xx = primBuffer.data(gkoff + 36 * idx + 18);

            auto gk_yy_xy = primBuffer.data(gkoff + 36 * idx + 19);

            auto gk_yy_xz = primBuffer.data(gkoff + 36 * idx + 20);

            auto gk_yy_yy = primBuffer.data(gkoff + 36 * idx + 21);

            auto gk_yy_yz = primBuffer.data(gkoff + 36 * idx + 22);

            auto gk_yy_zz = primBuffer.data(gkoff + 36 * idx + 23);

            auto gk_yz_xx = primBuffer.data(gkoff + 36 * idx + 24);

            auto gk_yz_xy = primBuffer.data(gkoff + 36 * idx + 25);

            auto gk_yz_xz = primBuffer.data(gkoff + 36 * idx + 26);

            auto gk_yz_yy = primBuffer.data(gkoff + 36 * idx + 27);

            auto gk_yz_yz = primBuffer.data(gkoff + 36 * idx + 28);

            auto gk_yz_zz = primBuffer.data(gkoff + 36 * idx + 29);

            auto gk_zz_xx = primBuffer.data(gkoff + 36 * idx + 30);

            auto gk_zz_xy = primBuffer.data(gkoff + 36 * idx + 31);

            auto gk_zz_xz = primBuffer.data(gkoff + 36 * idx + 32);

            auto gk_zz_yy = primBuffer.data(gkoff + 36 * idx + 33);

            auto gk_zz_yz = primBuffer.data(gkoff + 36 * idx + 34);

            auto gk_zz_zz = primBuffer.data(gkoff + 36 * idx + 35);

            // set up pointers to (P|g(r,r')|F)^(m) integrals

            auto g20_x_xxx = primBuffer.data(g20off + 30 * idx);

            auto g20_x_xxy = primBuffer.data(g20off + 30 * idx + 1);

            auto g20_x_xxz = primBuffer.data(g20off + 30 * idx + 2);

            auto g20_x_xyy = primBuffer.data(g20off + 30 * idx + 3);

            auto g20_x_xyz = primBuffer.data(g20off + 30 * idx + 4);

            auto g20_x_xzz = primBuffer.data(g20off + 30 * idx + 5);

            auto g20_x_yyy = primBuffer.data(g20off + 30 * idx + 6);

            auto g20_x_yyz = primBuffer.data(g20off + 30 * idx + 7);

            auto g20_x_yzz = primBuffer.data(g20off + 30 * idx + 8);

            auto g20_x_zzz = primBuffer.data(g20off + 30 * idx + 9);

            auto g20_y_xxx = primBuffer.data(g20off + 30 * idx + 10);

            auto g20_y_xxy = primBuffer.data(g20off + 30 * idx + 11);

            auto g20_y_xxz = primBuffer.data(g20off + 30 * idx + 12);

            auto g20_y_xyy = primBuffer.data(g20off + 30 * idx + 13);

            auto g20_y_xyz = primBuffer.data(g20off + 30 * idx + 14);

            auto g20_y_xzz = primBuffer.data(g20off + 30 * idx + 15);

            auto g20_y_yyy = primBuffer.data(g20off + 30 * idx + 16);

            auto g20_y_yyz = primBuffer.data(g20off + 30 * idx + 17);

            auto g20_y_yzz = primBuffer.data(g20off + 30 * idx + 18);

            auto g20_y_zzz = primBuffer.data(g20off + 30 * idx + 19);

            auto g20_z_xxx = primBuffer.data(g20off + 30 * idx + 20);

            auto g20_z_xxy = primBuffer.data(g20off + 30 * idx + 21);

            auto g20_z_xxz = primBuffer.data(g20off + 30 * idx + 22);

            auto g20_z_xyy = primBuffer.data(g20off + 30 * idx + 23);

            auto g20_z_xyz = primBuffer.data(g20off + 30 * idx + 24);

            auto g20_z_xzz = primBuffer.data(g20off + 30 * idx + 25);

            auto g20_z_yyy = primBuffer.data(g20off + 30 * idx + 26);

            auto g20_z_yyz = primBuffer.data(g20off + 30 * idx + 27);

            auto g20_z_yzz = primBuffer.data(g20off + 30 * idx + 28);

            auto g20_z_zzz = primBuffer.data(g20off + 30 * idx + 29);

            // set up pointers to (P|g(r,r')|F)^(m+1) integrals

            auto g21_x_xxx = primBuffer.data(g21off + 30 * idx);

            auto g21_x_xxy = primBuffer.data(g21off + 30 * idx + 1);

            auto g21_x_xxz = primBuffer.data(g21off + 30 * idx + 2);

            auto g21_x_xyy = primBuffer.data(g21off + 30 * idx + 3);

            auto g21_x_xyz = primBuffer.data(g21off + 30 * idx + 4);

            auto g21_x_xzz = primBuffer.data(g21off + 30 * idx + 5);

            auto g21_x_yyy = primBuffer.data(g21off + 30 * idx + 6);

            auto g21_x_yyz = primBuffer.data(g21off + 30 * idx + 7);

            auto g21_x_yzz = primBuffer.data(g21off + 30 * idx + 8);

            auto g21_x_zzz = primBuffer.data(g21off + 30 * idx + 9);

            auto g21_y_xxx = primBuffer.data(g21off + 30 * idx + 10);

            auto g21_y_xxy = primBuffer.data(g21off + 30 * idx + 11);

            auto g21_y_xxz = primBuffer.data(g21off + 30 * idx + 12);

            auto g21_y_xyy = primBuffer.data(g21off + 30 * idx + 13);

            auto g21_y_xyz = primBuffer.data(g21off + 30 * idx + 14);

            auto g21_y_xzz = primBuffer.data(g21off + 30 * idx + 15);

            auto g21_y_yyy = primBuffer.data(g21off + 30 * idx + 16);

            auto g21_y_yyz = primBuffer.data(g21off + 30 * idx + 17);

            auto g21_y_yzz = primBuffer.data(g21off + 30 * idx + 18);

            auto g21_y_zzz = primBuffer.data(g21off + 30 * idx + 19);

            auto g21_z_xxx = primBuffer.data(g21off + 30 * idx + 20);

            auto g21_z_xxy = primBuffer.data(g21off + 30 * idx + 21);

            auto g21_z_xxz = primBuffer.data(g21off + 30 * idx + 22);

            auto g21_z_xyy = primBuffer.data(g21off + 30 * idx + 23);

            auto g21_z_xyz = primBuffer.data(g21off + 30 * idx + 24);

            auto g21_z_xzz = primBuffer.data(g21off + 30 * idx + 25);

            auto g21_z_yyy = primBuffer.data(g21off + 30 * idx + 26);

            auto g21_z_yyz = primBuffer.data(g21off + 30 * idx + 27);

            auto g21_z_yzz = primBuffer.data(g21off + 30 * idx + 28);

            auto g21_z_zzz = primBuffer.data(g21off + 30 * idx + 29);

            // set up pointers to (D|g(r,r')|F)^(m+1) integrals

            auto g1_xx_xxx = primBuffer.data(g1off + 60 * idx);

            auto g1_xx_xxy = primBuffer.data(g1off + 60 * idx + 1);

            auto g1_xx_xxz = primBuffer.data(g1off + 60 * idx + 2);

            auto g1_xx_xyy = primBuffer.data(g1off + 60 * idx + 3);

            auto g1_xx_xyz = primBuffer.data(g1off + 60 * idx + 4);

            auto g1_xx_xzz = primBuffer.data(g1off + 60 * idx + 5);

            auto g1_xx_yyy = primBuffer.data(g1off + 60 * idx + 6);

            auto g1_xx_yyz = primBuffer.data(g1off + 60 * idx + 7);

            auto g1_xx_yzz = primBuffer.data(g1off + 60 * idx + 8);

            auto g1_xx_zzz = primBuffer.data(g1off + 60 * idx + 9);

            auto g1_xy_xxx = primBuffer.data(g1off + 60 * idx + 10);

            auto g1_xy_xxy = primBuffer.data(g1off + 60 * idx + 11);

            auto g1_xy_xxz = primBuffer.data(g1off + 60 * idx + 12);

            auto g1_xy_xyy = primBuffer.data(g1off + 60 * idx + 13);

            auto g1_xy_xyz = primBuffer.data(g1off + 60 * idx + 14);

            auto g1_xy_xzz = primBuffer.data(g1off + 60 * idx + 15);

            auto g1_xy_yyy = primBuffer.data(g1off + 60 * idx + 16);

            auto g1_xy_yyz = primBuffer.data(g1off + 60 * idx + 17);

            auto g1_xy_yzz = primBuffer.data(g1off + 60 * idx + 18);

            auto g1_xy_zzz = primBuffer.data(g1off + 60 * idx + 19);

            auto g1_xz_xxx = primBuffer.data(g1off + 60 * idx + 20);

            auto g1_xz_xxy = primBuffer.data(g1off + 60 * idx + 21);

            auto g1_xz_xxz = primBuffer.data(g1off + 60 * idx + 22);

            auto g1_xz_xyy = primBuffer.data(g1off + 60 * idx + 23);

            auto g1_xz_xyz = primBuffer.data(g1off + 60 * idx + 24);

            auto g1_xz_xzz = primBuffer.data(g1off + 60 * idx + 25);

            auto g1_xz_yyy = primBuffer.data(g1off + 60 * idx + 26);

            auto g1_xz_yyz = primBuffer.data(g1off + 60 * idx + 27);

            auto g1_xz_yzz = primBuffer.data(g1off + 60 * idx + 28);

            auto g1_xz_zzz = primBuffer.data(g1off + 60 * idx + 29);

            auto g1_yy_xxx = primBuffer.data(g1off + 60 * idx + 30);

            auto g1_yy_xxy = primBuffer.data(g1off + 60 * idx + 31);

            auto g1_yy_xxz = primBuffer.data(g1off + 60 * idx + 32);

            auto g1_yy_xyy = primBuffer.data(g1off + 60 * idx + 33);

            auto g1_yy_xyz = primBuffer.data(g1off + 60 * idx + 34);

            auto g1_yy_xzz = primBuffer.data(g1off + 60 * idx + 35);

            auto g1_yy_yyy = primBuffer.data(g1off + 60 * idx + 36);

            auto g1_yy_yyz = primBuffer.data(g1off + 60 * idx + 37);

            auto g1_yy_yzz = primBuffer.data(g1off + 60 * idx + 38);

            auto g1_yy_zzz = primBuffer.data(g1off + 60 * idx + 39);

            auto g1_yz_xxx = primBuffer.data(g1off + 60 * idx + 40);

            auto g1_yz_xxy = primBuffer.data(g1off + 60 * idx + 41);

            auto g1_yz_xxz = primBuffer.data(g1off + 60 * idx + 42);

            auto g1_yz_xyy = primBuffer.data(g1off + 60 * idx + 43);

            auto g1_yz_xyz = primBuffer.data(g1off + 60 * idx + 44);

            auto g1_yz_xzz = primBuffer.data(g1off + 60 * idx + 45);

            auto g1_yz_yyy = primBuffer.data(g1off + 60 * idx + 46);

            auto g1_yz_yyz = primBuffer.data(g1off + 60 * idx + 47);

            auto g1_yz_yzz = primBuffer.data(g1off + 60 * idx + 48);

            auto g1_yz_zzz = primBuffer.data(g1off + 60 * idx + 49);

            auto g1_zz_xxx = primBuffer.data(g1off + 60 * idx + 50);

            auto g1_zz_xxy = primBuffer.data(g1off + 60 * idx + 51);

            auto g1_zz_xxz = primBuffer.data(g1off + 60 * idx + 52);

            auto g1_zz_xyy = primBuffer.data(g1off + 60 * idx + 53);

            auto g1_zz_xyz = primBuffer.data(g1off + 60 * idx + 54);

            auto g1_zz_xzz = primBuffer.data(g1off + 60 * idx + 55);

            auto g1_zz_yyy = primBuffer.data(g1off + 60 * idx + 56);

            auto g1_zz_yyz = primBuffer.data(g1off + 60 * idx + 57);

            auto g1_zz_yzz = primBuffer.data(g1off + 60 * idx + 58);

            auto g1_zz_zzz = primBuffer.data(g1off + 60 * idx + 59);

            // set up pointers to (F|g(r,r')|F)^(m) integrals

            auto g_xxx_xxx = primBuffer.data(goff + 100 * idx);

            auto g_xxx_xxy = primBuffer.data(goff + 100 * idx + 1);

            auto g_xxx_xxz = primBuffer.data(goff + 100 * idx + 2);

            auto g_xxx_xyy = primBuffer.data(goff + 100 * idx + 3);

            auto g_xxx_xyz = primBuffer.data(goff + 100 * idx + 4);

            auto g_xxx_xzz = primBuffer.data(goff + 100 * idx + 5);

            auto g_xxx_yyy = primBuffer.data(goff + 100 * idx + 6);

            auto g_xxx_yyz = primBuffer.data(goff + 100 * idx + 7);

            auto g_xxx_yzz = primBuffer.data(goff + 100 * idx + 8);

            auto g_xxx_zzz = primBuffer.data(goff + 100 * idx + 9);

            auto g_xxy_xxx = primBuffer.data(goff + 100 * idx + 10);

            auto g_xxy_xxy = primBuffer.data(goff + 100 * idx + 11);

            auto g_xxy_xxz = primBuffer.data(goff + 100 * idx + 12);

            auto g_xxy_xyy = primBuffer.data(goff + 100 * idx + 13);

            auto g_xxy_xyz = primBuffer.data(goff + 100 * idx + 14);

            auto g_xxy_xzz = primBuffer.data(goff + 100 * idx + 15);

            auto g_xxy_yyy = primBuffer.data(goff + 100 * idx + 16);

            auto g_xxy_yyz = primBuffer.data(goff + 100 * idx + 17);

            auto g_xxy_yzz = primBuffer.data(goff + 100 * idx + 18);

            auto g_xxy_zzz = primBuffer.data(goff + 100 * idx + 19);

            auto g_xxz_xxx = primBuffer.data(goff + 100 * idx + 20);

            auto g_xxz_xxy = primBuffer.data(goff + 100 * idx + 21);

            auto g_xxz_xxz = primBuffer.data(goff + 100 * idx + 22);

            auto g_xxz_xyy = primBuffer.data(goff + 100 * idx + 23);

            auto g_xxz_xyz = primBuffer.data(goff + 100 * idx + 24);

            auto g_xxz_xzz = primBuffer.data(goff + 100 * idx + 25);

            auto g_xxz_yyy = primBuffer.data(goff + 100 * idx + 26);

            auto g_xxz_yyz = primBuffer.data(goff + 100 * idx + 27);

            auto g_xxz_yzz = primBuffer.data(goff + 100 * idx + 28);

            auto g_xxz_zzz = primBuffer.data(goff + 100 * idx + 29);

            auto g_xyy_xxx = primBuffer.data(goff + 100 * idx + 30);

            auto g_xyy_xxy = primBuffer.data(goff + 100 * idx + 31);

            auto g_xyy_xxz = primBuffer.data(goff + 100 * idx + 32);

            auto g_xyy_xyy = primBuffer.data(goff + 100 * idx + 33);

            auto g_xyy_xyz = primBuffer.data(goff + 100 * idx + 34);

            auto g_xyy_xzz = primBuffer.data(goff + 100 * idx + 35);

            auto g_xyy_yyy = primBuffer.data(goff + 100 * idx + 36);

            auto g_xyy_yyz = primBuffer.data(goff + 100 * idx + 37);

            auto g_xyy_yzz = primBuffer.data(goff + 100 * idx + 38);

            auto g_xyy_zzz = primBuffer.data(goff + 100 * idx + 39);

            auto g_xyz_xxx = primBuffer.data(goff + 100 * idx + 40);

            auto g_xyz_xxy = primBuffer.data(goff + 100 * idx + 41);

            auto g_xyz_xxz = primBuffer.data(goff + 100 * idx + 42);

            auto g_xyz_xyy = primBuffer.data(goff + 100 * idx + 43);

            auto g_xyz_xyz = primBuffer.data(goff + 100 * idx + 44);

            auto g_xyz_xzz = primBuffer.data(goff + 100 * idx + 45);

            auto g_xyz_yyy = primBuffer.data(goff + 100 * idx + 46);

            auto g_xyz_yyz = primBuffer.data(goff + 100 * idx + 47);

            auto g_xyz_yzz = primBuffer.data(goff + 100 * idx + 48);

            auto g_xyz_zzz = primBuffer.data(goff + 100 * idx + 49);

            auto g_xzz_xxx = primBuffer.data(goff + 100 * idx + 50);

            auto g_xzz_xxy = primBuffer.data(goff + 100 * idx + 51);

            auto g_xzz_xxz = primBuffer.data(goff + 100 * idx + 52);

            auto g_xzz_xyy = primBuffer.data(goff + 100 * idx + 53);

            auto g_xzz_xyz = primBuffer.data(goff + 100 * idx + 54);

            auto g_xzz_xzz = primBuffer.data(goff + 100 * idx + 55);

            auto g_xzz_yyy = primBuffer.data(goff + 100 * idx + 56);

            auto g_xzz_yyz = primBuffer.data(goff + 100 * idx + 57);

            auto g_xzz_yzz = primBuffer.data(goff + 100 * idx + 58);

            auto g_xzz_zzz = primBuffer.data(goff + 100 * idx + 59);

            auto g_yyy_xxx = primBuffer.data(goff + 100 * idx + 60);

            auto g_yyy_xxy = primBuffer.data(goff + 100 * idx + 61);

            auto g_yyy_xxz = primBuffer.data(goff + 100 * idx + 62);

            auto g_yyy_xyy = primBuffer.data(goff + 100 * idx + 63);

            auto g_yyy_xyz = primBuffer.data(goff + 100 * idx + 64);

            auto g_yyy_xzz = primBuffer.data(goff + 100 * idx + 65);

            auto g_yyy_yyy = primBuffer.data(goff + 100 * idx + 66);

            auto g_yyy_yyz = primBuffer.data(goff + 100 * idx + 67);

            auto g_yyy_yzz = primBuffer.data(goff + 100 * idx + 68);

            auto g_yyy_zzz = primBuffer.data(goff + 100 * idx + 69);

            auto g_yyz_xxx = primBuffer.data(goff + 100 * idx + 70);

            auto g_yyz_xxy = primBuffer.data(goff + 100 * idx + 71);

            auto g_yyz_xxz = primBuffer.data(goff + 100 * idx + 72);

            auto g_yyz_xyy = primBuffer.data(goff + 100 * idx + 73);

            auto g_yyz_xyz = primBuffer.data(goff + 100 * idx + 74);

            auto g_yyz_xzz = primBuffer.data(goff + 100 * idx + 75);

            auto g_yyz_yyy = primBuffer.data(goff + 100 * idx + 76);

            auto g_yyz_yyz = primBuffer.data(goff + 100 * idx + 77);

            auto g_yyz_yzz = primBuffer.data(goff + 100 * idx + 78);

            auto g_yyz_zzz = primBuffer.data(goff + 100 * idx + 79);

            auto g_yzz_xxx = primBuffer.data(goff + 100 * idx + 80);

            auto g_yzz_xxy = primBuffer.data(goff + 100 * idx + 81);

            auto g_yzz_xxz = primBuffer.data(goff + 100 * idx + 82);

            auto g_yzz_xyy = primBuffer.data(goff + 100 * idx + 83);

            auto g_yzz_xyz = primBuffer.data(goff + 100 * idx + 84);

            auto g_yzz_xzz = primBuffer.data(goff + 100 * idx + 85);

            auto g_yzz_yyy = primBuffer.data(goff + 100 * idx + 86);

            auto g_yzz_yyz = primBuffer.data(goff + 100 * idx + 87);

            auto g_yzz_yzz = primBuffer.data(goff + 100 * idx + 88);

            auto g_yzz_zzz = primBuffer.data(goff + 100 * idx + 89);

            auto g_zzz_xxx = primBuffer.data(goff + 100 * idx + 90);

            auto g_zzz_xxy = primBuffer.data(goff + 100 * idx + 91);

            auto g_zzz_xxz = primBuffer.data(goff + 100 * idx + 92);

            auto g_zzz_xyy = primBuffer.data(goff + 100 * idx + 93);

            auto g_zzz_xyz = primBuffer.data(goff + 100 * idx + 94);

            auto g_zzz_xzz = primBuffer.data(goff + 100 * idx + 95);

            auto g_zzz_yyy = primBuffer.data(goff + 100 * idx + 96);

            auto g_zzz_yyz = primBuffer.data(goff + 100 * idx + 97);

            auto g_zzz_yzz = primBuffer.data(goff + 100 * idx + 98);

            auto g_zzz_zzz = primBuffer.data(goff + 100 * idx + 99);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xx_xx, gk_xx_xy,\
                                         gk_xx_xz, gk_xx_yy, gk_xx_yz, gk_xx_zz,\
                                         gk_xy_xx, gk_xy_xy, gk_xy_xz, gk_xy_yy,\
                                         gk_xy_yz, gk_xy_zz, gk_xz_xx, gk_xz_xy,\
                                         gk_xz_xz, gk_xz_yy, gk_xz_yz, gk_xz_zz,\
                                         gk_yy_xx, gk_yy_xy, gk_yy_xz, gk_yy_yy,\
                                         gk_yy_yz, gk_yy_zz, gk_yz_xx, gk_yz_xy,\
                                         gk_yz_xz, gk_yz_yy, gk_yz_yz, gk_yz_zz,\
                                         gk_zz_xx, gk_zz_xy, gk_zz_xz, gk_zz_yy,\
                                         gk_zz_yz, gk_zz_zz, g20_x_xxx, g20_x_xxy,\
                                         g20_x_xxz, g20_x_xyy, g20_x_xyz, g20_x_xzz,\
                                         g20_x_yyy, g20_x_yyz, g20_x_yzz, g20_x_zzz,\
                                         g20_y_xxx, g20_y_xxy, g20_y_xxz, g20_y_xyy,\
                                         g20_y_xyz, g20_y_xzz, g20_y_yyy, g20_y_yyz,\
                                         g20_y_yzz, g20_y_zzz, g20_z_xxx, g20_z_xxy,\
                                         g20_z_xxz, g20_z_xyy, g20_z_xyz, g20_z_xzz,\
                                         g20_z_yyy, g20_z_yyz, g20_z_yzz, g20_z_zzz,\
                                         g21_x_xxx, g21_x_xxy, g21_x_xxz, g21_x_xyy,\
                                         g21_x_xyz, g21_x_xzz, g21_x_yyy, g21_x_yyz,\
                                         g21_x_yzz, g21_x_zzz, g21_y_xxx, g21_y_xxy,\
                                         g21_y_xxz, g21_y_xyy, g21_y_xyz, g21_y_xzz,\
                                         g21_y_yyy, g21_y_yyz, g21_y_yzz, g21_y_zzz,\
                                         g21_z_xxx, g21_z_xxy, g21_z_xxz, g21_z_xyy,\
                                         g21_z_xyz, g21_z_xzz, g21_z_yyy, g21_z_yyz,\
                                         g21_z_yzz, g21_z_zzz, g1_xx_xxx, g1_xx_xxy,\
                                         g1_xx_xxz, g1_xx_xyy, g1_xx_xyz, g1_xx_xzz,\
                                         g1_xx_yyy, g1_xx_yyz, g1_xx_yzz, g1_xx_zzz,\
                                         g1_xy_xxx, g1_xy_xxy, g1_xy_xxz, g1_xy_xyy,\
                                         g1_xy_xyz, g1_xy_xzz, g1_xy_yyy, g1_xy_yyz,\
                                         g1_xy_yzz, g1_xy_zzz, g1_xz_xxx, g1_xz_xxy,\
                                         g1_xz_xxz, g1_xz_xyy, g1_xz_xyz, g1_xz_xzz,\
                                         g1_xz_yyy, g1_xz_yyz, g1_xz_yzz, g1_xz_zzz,\
                                         g1_yy_xxx, g1_yy_xxy, g1_yy_xxz, g1_yy_xyy,\
                                         g1_yy_xyz, g1_yy_xzz, g1_yy_yyy, g1_yy_yyz,\
                                         g1_yy_yzz, g1_yy_zzz, g1_yz_xxx, g1_yz_xxy,\
                                         g1_yz_xxz, g1_yz_xyy, g1_yz_xyz, g1_yz_xzz,\
                                         g1_yz_yyy, g1_yz_yyz, g1_yz_yzz, g1_yz_zzz,\
                                         g1_zz_xxx, g1_zz_xxy, g1_zz_xxz, g1_zz_xyy,\
                                         g1_zz_xyz, g1_zz_xzz, g1_zz_yyy, g1_zz_yyz,\
                                         g1_zz_yzz, g1_zz_zzz, g_xxx_xxx, g_xxx_xxy,\
                                         g_xxx_xxz, g_xxx_xyy, g_xxx_xyz, g_xxx_xzz,\
                                         g_xxx_yyy, g_xxx_yyz, g_xxx_yzz, g_xxx_zzz,\
                                         g_xxy_xxx, g_xxy_xxy, g_xxy_xxz, g_xxy_xyy,\
                                         g_xxy_xyz, g_xxy_xzz, g_xxy_yyy, g_xxy_yyz,\
                                         g_xxy_yzz, g_xxy_zzz, g_xxz_xxx, g_xxz_xxy,\
                                         g_xxz_xxz, g_xxz_xyy, g_xxz_xyz, g_xxz_xzz,\
                                         g_xxz_yyy, g_xxz_yyz, g_xxz_yzz, g_xxz_zzz,\
                                         g_xyy_xxx, g_xyy_xxy, g_xyy_xxz, g_xyy_xyy,\
                                         g_xyy_xyz, g_xyy_xzz, g_xyy_yyy, g_xyy_yyz,\
                                         g_xyy_yzz, g_xyy_zzz, g_xyz_xxx, g_xyz_xxy,\
                                         g_xyz_xxz, g_xyz_xyy, g_xyz_xyz, g_xyz_xzz,\
                                         g_xyz_yyy, g_xyz_yyz, g_xyz_yzz, g_xyz_zzz,\
                                         g_xzz_xxx, g_xzz_xxy, g_xzz_xxz, g_xzz_xyy,\
                                         g_xzz_xyz, g_xzz_xzz, g_xzz_yyy, g_xzz_yyz,\
                                         g_xzz_yzz, g_xzz_zzz, g_yyy_xxx, g_yyy_xxy,\
                                         g_yyy_xxz, g_yyy_xyy, g_yyy_xyz, g_yyy_xzz,\
                                         g_yyy_yyy, g_yyy_yyz, g_yyy_yzz, g_yyy_zzz,\
                                         g_yyz_xxx, g_yyz_xxy, g_yyz_xxz, g_yyz_xyy,\
                                         g_yyz_xyz, g_yyz_xzz, g_yyz_yyy, g_yyz_yyz,\
                                         g_yyz_yzz, g_yyz_zzz, g_yzz_xxx, g_yzz_xxy,\
                                         g_yzz_xxz, g_yzz_xyy, g_yzz_xyz, g_yzz_xzz,\
                                         g_yzz_yyy, g_yzz_yyz, g_yzz_yzz, g_yzz_zzz,\
                                         g_zzz_xxx, g_zzz_xxy, g_zzz_xxz, g_zzz_xyy,\
                                         g_zzz_xyz, g_zzz_xzz, g_zzz_yyy, g_zzz_yyz,\
                                         g_zzz_yzz, g_zzz_zzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxx_xxx[k] = fr * g1_xx_xxx[k] + f2g * (2.0 * g20_x_xxx[k] - 2.0 * fgz * g21_x_xxx[k]) + 3.0 * f2t * gk_xx_xx[k];

                g_xxx_xxy[k] = fr * g1_xx_xxy[k] + f2g * (2.0 * g20_x_xxy[k] - 2.0 * fgz * g21_x_xxy[k]) + 2.0 * f2t * gk_xx_xy[k];

                g_xxx_xxz[k] = fr * g1_xx_xxz[k] + f2g * (2.0 * g20_x_xxz[k] - 2.0 * fgz * g21_x_xxz[k]) + 2.0 * f2t * gk_xx_xz[k];

                g_xxx_xyy[k] = fr * g1_xx_xyy[k] + f2g * (2.0 * g20_x_xyy[k] - 2.0 * fgz * g21_x_xyy[k]) + f2t * gk_xx_yy[k];

                g_xxx_xyz[k] = fr * g1_xx_xyz[k] + f2g * (2.0 * g20_x_xyz[k] - 2.0 * fgz * g21_x_xyz[k]) + f2t * gk_xx_yz[k];

                g_xxx_xzz[k] = fr * g1_xx_xzz[k] + f2g * (2.0 * g20_x_xzz[k] - 2.0 * fgz * g21_x_xzz[k]) + f2t * gk_xx_zz[k];

                g_xxx_yyy[k] = fr * g1_xx_yyy[k] + f2g * (2.0 * g20_x_yyy[k] - 2.0 * fgz * g21_x_yyy[k]);

                g_xxx_yyz[k] = fr * g1_xx_yyz[k] + f2g * (2.0 * g20_x_yyz[k] - 2.0 * fgz * g21_x_yyz[k]);

                g_xxx_yzz[k] = fr * g1_xx_yzz[k] + f2g * (2.0 * g20_x_yzz[k] - 2.0 * fgz * g21_x_yzz[k]);

                g_xxx_zzz[k] = fr * g1_xx_zzz[k] + f2g * (2.0 * g20_x_zzz[k] - 2.0 * fgz * g21_x_zzz[k]);

                g_xxy_xxx[k] = fr * g1_xy_xxx[k] + f2g * (g20_y_xxx[k] - fgz * g21_y_xxx[k]) + 3.0 * f2t * gk_xy_xx[k];

                g_xxy_xxy[k] = fr * g1_xy_xxy[k] + f2g * (g20_y_xxy[k] - fgz * g21_y_xxy[k]) + 2.0 * f2t * gk_xy_xy[k];

                g_xxy_xxz[k] = fr * g1_xy_xxz[k] + f2g * (g20_y_xxz[k] - fgz * g21_y_xxz[k]) + 2.0 * f2t * gk_xy_xz[k];

                g_xxy_xyy[k] = fr * g1_xy_xyy[k] + f2g * (g20_y_xyy[k] - fgz * g21_y_xyy[k]) + f2t * gk_xy_yy[k];

                g_xxy_xyz[k] = fr * g1_xy_xyz[k] + f2g * (g20_y_xyz[k] - fgz * g21_y_xyz[k]) + f2t * gk_xy_yz[k];

                g_xxy_xzz[k] = fr * g1_xy_xzz[k] + f2g * (g20_y_xzz[k] - fgz * g21_y_xzz[k]) + f2t * gk_xy_zz[k];

                g_xxy_yyy[k] = fr * g1_xy_yyy[k] + f2g * (g20_y_yyy[k] - fgz * g21_y_yyy[k]);

                g_xxy_yyz[k] = fr * g1_xy_yyz[k] + f2g * (g20_y_yyz[k] - fgz * g21_y_yyz[k]);

                g_xxy_yzz[k] = fr * g1_xy_yzz[k] + f2g * (g20_y_yzz[k] - fgz * g21_y_yzz[k]);

                g_xxy_zzz[k] = fr * g1_xy_zzz[k] + f2g * (g20_y_zzz[k] - fgz * g21_y_zzz[k]);

                g_xxz_xxx[k] = fr * g1_xz_xxx[k] + f2g * (g20_z_xxx[k] - fgz * g21_z_xxx[k]) + 3.0 * f2t * gk_xz_xx[k];

                g_xxz_xxy[k] = fr * g1_xz_xxy[k] + f2g * (g20_z_xxy[k] - fgz * g21_z_xxy[k]) + 2.0 * f2t * gk_xz_xy[k];

                g_xxz_xxz[k] = fr * g1_xz_xxz[k] + f2g * (g20_z_xxz[k] - fgz * g21_z_xxz[k]) + 2.0 * f2t * gk_xz_xz[k];

                g_xxz_xyy[k] = fr * g1_xz_xyy[k] + f2g * (g20_z_xyy[k] - fgz * g21_z_xyy[k]) + f2t * gk_xz_yy[k];

                g_xxz_xyz[k] = fr * g1_xz_xyz[k] + f2g * (g20_z_xyz[k] - fgz * g21_z_xyz[k]) + f2t * gk_xz_yz[k];

                g_xxz_xzz[k] = fr * g1_xz_xzz[k] + f2g * (g20_z_xzz[k] - fgz * g21_z_xzz[k]) + f2t * gk_xz_zz[k];

                g_xxz_yyy[k] = fr * g1_xz_yyy[k] + f2g * (g20_z_yyy[k] - fgz * g21_z_yyy[k]);

                g_xxz_yyz[k] = fr * g1_xz_yyz[k] + f2g * (g20_z_yyz[k] - fgz * g21_z_yyz[k]);

                g_xxz_yzz[k] = fr * g1_xz_yzz[k] + f2g * (g20_z_yzz[k] - fgz * g21_z_yzz[k]);

                g_xxz_zzz[k] = fr * g1_xz_zzz[k] + f2g * (g20_z_zzz[k] - fgz * g21_z_zzz[k]);

                g_xyy_xxx[k] = fr * g1_yy_xxx[k] + 3.0 * f2t * gk_yy_xx[k];

                g_xyy_xxy[k] = fr * g1_yy_xxy[k] + 2.0 * f2t * gk_yy_xy[k];

                g_xyy_xxz[k] = fr * g1_yy_xxz[k] + 2.0 * f2t * gk_yy_xz[k];

                g_xyy_xyy[k] = fr * g1_yy_xyy[k] + f2t * gk_yy_yy[k];

                g_xyy_xyz[k] = fr * g1_yy_xyz[k] + f2t * gk_yy_yz[k];

                g_xyy_xzz[k] = fr * g1_yy_xzz[k] + f2t * gk_yy_zz[k];

                g_xyy_yyy[k] = fr * g1_yy_yyy[k];

                g_xyy_yyz[k] = fr * g1_yy_yyz[k];

                g_xyy_yzz[k] = fr * g1_yy_yzz[k];

                g_xyy_zzz[k] = fr * g1_yy_zzz[k];

                g_xyz_xxx[k] = fr * g1_yz_xxx[k] + 3.0 * f2t * gk_yz_xx[k];

                g_xyz_xxy[k] = fr * g1_yz_xxy[k] + 2.0 * f2t * gk_yz_xy[k];

                g_xyz_xxz[k] = fr * g1_yz_xxz[k] + 2.0 * f2t * gk_yz_xz[k];

                g_xyz_xyy[k] = fr * g1_yz_xyy[k] + f2t * gk_yz_yy[k];

                g_xyz_xyz[k] = fr * g1_yz_xyz[k] + f2t * gk_yz_yz[k];

                g_xyz_xzz[k] = fr * g1_yz_xzz[k] + f2t * gk_yz_zz[k];

                g_xyz_yyy[k] = fr * g1_yz_yyy[k];

                g_xyz_yyz[k] = fr * g1_yz_yyz[k];

                g_xyz_yzz[k] = fr * g1_yz_yzz[k];

                g_xyz_zzz[k] = fr * g1_yz_zzz[k];

                g_xzz_xxx[k] = fr * g1_zz_xxx[k] + 3.0 * f2t * gk_zz_xx[k];

                g_xzz_xxy[k] = fr * g1_zz_xxy[k] + 2.0 * f2t * gk_zz_xy[k];

                g_xzz_xxz[k] = fr * g1_zz_xxz[k] + 2.0 * f2t * gk_zz_xz[k];

                g_xzz_xyy[k] = fr * g1_zz_xyy[k] + f2t * gk_zz_yy[k];

                g_xzz_xyz[k] = fr * g1_zz_xyz[k] + f2t * gk_zz_yz[k];

                g_xzz_xzz[k] = fr * g1_zz_xzz[k] + f2t * gk_zz_zz[k];

                g_xzz_yyy[k] = fr * g1_zz_yyy[k];

                g_xzz_yyz[k] = fr * g1_zz_yyz[k];

                g_xzz_yzz[k] = fr * g1_zz_yzz[k];

                g_xzz_zzz[k] = fr * g1_zz_zzz[k];

                // leading y component

                fr = pay[k];

                g_yyy_xxx[k] = fr * g1_yy_xxx[k] + f2g * (2.0 * g20_y_xxx[k] - 2.0 * fgz * g21_y_xxx[k]);

                g_yyy_xxy[k] = fr * g1_yy_xxy[k] + f2g * (2.0 * g20_y_xxy[k] - 2.0 * fgz * g21_y_xxy[k]) + f2t * gk_yy_xx[k];

                g_yyy_xxz[k] = fr * g1_yy_xxz[k] + f2g * (2.0 * g20_y_xxz[k] - 2.0 * fgz * g21_y_xxz[k]);

                g_yyy_xyy[k] = fr * g1_yy_xyy[k] + f2g * (2.0 * g20_y_xyy[k] - 2.0 * fgz * g21_y_xyy[k]) + 2.0 * f2t * gk_yy_xy[k];

                g_yyy_xyz[k] = fr * g1_yy_xyz[k] + f2g * (2.0 * g20_y_xyz[k] - 2.0 * fgz * g21_y_xyz[k]) + f2t * gk_yy_xz[k];

                g_yyy_xzz[k] = fr * g1_yy_xzz[k] + f2g * (2.0 * g20_y_xzz[k] - 2.0 * fgz * g21_y_xzz[k]);

                g_yyy_yyy[k] = fr * g1_yy_yyy[k] + f2g * (2.0 * g20_y_yyy[k] - 2.0 * fgz * g21_y_yyy[k]) + 3.0 * f2t * gk_yy_yy[k];

                g_yyy_yyz[k] = fr * g1_yy_yyz[k] + f2g * (2.0 * g20_y_yyz[k] - 2.0 * fgz * g21_y_yyz[k]) + 2.0 * f2t * gk_yy_yz[k];

                g_yyy_yzz[k] = fr * g1_yy_yzz[k] + f2g * (2.0 * g20_y_yzz[k] - 2.0 * fgz * g21_y_yzz[k]) + f2t * gk_yy_zz[k];

                g_yyy_zzz[k] = fr * g1_yy_zzz[k] + f2g * (2.0 * g20_y_zzz[k] - 2.0 * fgz * g21_y_zzz[k]);

                g_yyz_xxx[k] = fr * g1_yz_xxx[k] + f2g * (g20_z_xxx[k] - fgz * g21_z_xxx[k]);

                g_yyz_xxy[k] = fr * g1_yz_xxy[k] + f2g * (g20_z_xxy[k] - fgz * g21_z_xxy[k]) + f2t * gk_yz_xx[k];

                g_yyz_xxz[k] = fr * g1_yz_xxz[k] + f2g * (g20_z_xxz[k] - fgz * g21_z_xxz[k]);

                g_yyz_xyy[k] = fr * g1_yz_xyy[k] + f2g * (g20_z_xyy[k] - fgz * g21_z_xyy[k]) + 2.0 * f2t * gk_yz_xy[k];

                g_yyz_xyz[k] = fr * g1_yz_xyz[k] + f2g * (g20_z_xyz[k] - fgz * g21_z_xyz[k]) + f2t * gk_yz_xz[k];

                g_yyz_xzz[k] = fr * g1_yz_xzz[k] + f2g * (g20_z_xzz[k] - fgz * g21_z_xzz[k]);

                g_yyz_yyy[k] = fr * g1_yz_yyy[k] + f2g * (g20_z_yyy[k] - fgz * g21_z_yyy[k]) + 3.0 * f2t * gk_yz_yy[k];

                g_yyz_yyz[k] = fr * g1_yz_yyz[k] + f2g * (g20_z_yyz[k] - fgz * g21_z_yyz[k]) + 2.0 * f2t * gk_yz_yz[k];

                g_yyz_yzz[k] = fr * g1_yz_yzz[k] + f2g * (g20_z_yzz[k] - fgz * g21_z_yzz[k]) + f2t * gk_yz_zz[k];

                g_yyz_zzz[k] = fr * g1_yz_zzz[k] + f2g * (g20_z_zzz[k] - fgz * g21_z_zzz[k]);

                g_yzz_xxx[k] = fr * g1_zz_xxx[k];

                g_yzz_xxy[k] = fr * g1_zz_xxy[k] + f2t * gk_zz_xx[k];

                g_yzz_xxz[k] = fr * g1_zz_xxz[k];

                g_yzz_xyy[k] = fr * g1_zz_xyy[k] + 2.0 * f2t * gk_zz_xy[k];

                g_yzz_xyz[k] = fr * g1_zz_xyz[k] + f2t * gk_zz_xz[k];

                g_yzz_xzz[k] = fr * g1_zz_xzz[k];

                g_yzz_yyy[k] = fr * g1_zz_yyy[k] + 3.0 * f2t * gk_zz_yy[k];

                g_yzz_yyz[k] = fr * g1_zz_yyz[k] + 2.0 * f2t * gk_zz_yz[k];

                g_yzz_yzz[k] = fr * g1_zz_yzz[k] + f2t * gk_zz_zz[k];

                g_yzz_zzz[k] = fr * g1_zz_zzz[k];

                // leading z component

                fr = paz[k];

                g_zzz_xxx[k] = fr * g1_zz_xxx[k] + f2g * (2.0 * g20_z_xxx[k] - 2.0 * fgz * g21_z_xxx[k]);

                g_zzz_xxy[k] = fr * g1_zz_xxy[k] + f2g * (2.0 * g20_z_xxy[k] - 2.0 * fgz * g21_z_xxy[k]);

                g_zzz_xxz[k] = fr * g1_zz_xxz[k] + f2g * (2.0 * g20_z_xxz[k] - 2.0 * fgz * g21_z_xxz[k]) + f2t * gk_zz_xx[k];

                g_zzz_xyy[k] = fr * g1_zz_xyy[k] + f2g * (2.0 * g20_z_xyy[k] - 2.0 * fgz * g21_z_xyy[k]);

                g_zzz_xyz[k] = fr * g1_zz_xyz[k] + f2g * (2.0 * g20_z_xyz[k] - 2.0 * fgz * g21_z_xyz[k]) + f2t * gk_zz_xy[k];

                g_zzz_xzz[k] = fr * g1_zz_xzz[k] + f2g * (2.0 * g20_z_xzz[k] - 2.0 * fgz * g21_z_xzz[k]) + 2.0 * f2t * gk_zz_xz[k];

                g_zzz_yyy[k] = fr * g1_zz_yyy[k] + f2g * (2.0 * g20_z_yyy[k] - 2.0 * fgz * g21_z_yyy[k]);

                g_zzz_yyz[k] = fr * g1_zz_yyz[k] + f2g * (2.0 * g20_z_yyz[k] - 2.0 * fgz * g21_z_yyz[k]) + f2t * gk_zz_yy[k];

                g_zzz_yzz[k] = fr * g1_zz_yzz[k] + f2g * (2.0 * g20_z_yzz[k] - 2.0 * fgz * g21_z_yzz[k]) + 2.0 * f2t * gk_zz_yz[k];

                g_zzz_zzz[k] = fr * g1_zz_zzz[k] + f2g * (2.0 * g20_z_zzz[k] - 2.0 * fgz * g21_z_zzz[k]) + 3.0 * f2t * gk_zz_zz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForSG(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  pbDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 0, 4);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 2, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 4);

            auto fz = osFactors.data(6 * idx + 5);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|D)^(m) integrals

            auto g20_0_xx = primBuffer.data(g20off + 6 * idx);

            auto g20_0_xy = primBuffer.data(g20off + 6 * idx + 1);

            auto g20_0_xz = primBuffer.data(g20off + 6 * idx + 2);

            auto g20_0_yy = primBuffer.data(g20off + 6 * idx + 3);

            auto g20_0_yz = primBuffer.data(g20off + 6 * idx + 4);

            auto g20_0_zz = primBuffer.data(g20off + 6 * idx + 5);

            // set up pointers to (S|g(r,r')|D)^(m+1) integrals

            auto g21_0_xx = primBuffer.data(g21off + 6 * idx);

            auto g21_0_xy = primBuffer.data(g21off + 6 * idx + 1);

            auto g21_0_xz = primBuffer.data(g21off + 6 * idx + 2);

            auto g21_0_yy = primBuffer.data(g21off + 6 * idx + 3);

            auto g21_0_yz = primBuffer.data(g21off + 6 * idx + 4);

            auto g21_0_zz = primBuffer.data(g21off + 6 * idx + 5);

            // set up pointers to (S|g(r,r')|F)^(m+1) integrals

            auto g1_0_xxx = primBuffer.data(g1off + 10 * idx);

            auto g1_0_xxy = primBuffer.data(g1off + 10 * idx + 1);

            auto g1_0_xxz = primBuffer.data(g1off + 10 * idx + 2);

            auto g1_0_xyy = primBuffer.data(g1off + 10 * idx + 3);

            auto g1_0_xyz = primBuffer.data(g1off + 10 * idx + 4);

            auto g1_0_xzz = primBuffer.data(g1off + 10 * idx + 5);

            auto g1_0_yyy = primBuffer.data(g1off + 10 * idx + 6);

            auto g1_0_yyz = primBuffer.data(g1off + 10 * idx + 7);

            auto g1_0_yzz = primBuffer.data(g1off + 10 * idx + 8);

            auto g1_0_zzz = primBuffer.data(g1off + 10 * idx + 9);

            // set up pointers to (S|g(r,r')|G)^(m) integrals

            auto g_0_xxxx = primBuffer.data(goff + 15 * idx);

            auto g_0_xxxy = primBuffer.data(goff + 15 * idx + 1);

            auto g_0_xxxz = primBuffer.data(goff + 15 * idx + 2);

            auto g_0_xxyy = primBuffer.data(goff + 15 * idx + 3);

            auto g_0_xxyz = primBuffer.data(goff + 15 * idx + 4);

            auto g_0_xxzz = primBuffer.data(goff + 15 * idx + 5);

            auto g_0_xyyy = primBuffer.data(goff + 15 * idx + 6);

            auto g_0_xyyz = primBuffer.data(goff + 15 * idx + 7);

            auto g_0_xyzz = primBuffer.data(goff + 15 * idx + 8);

            auto g_0_xzzz = primBuffer.data(goff + 15 * idx + 9);

            auto g_0_yyyy = primBuffer.data(goff + 15 * idx + 10);

            auto g_0_yyyz = primBuffer.data(goff + 15 * idx + 11);

            auto g_0_yyzz = primBuffer.data(goff + 15 * idx + 12);

            auto g_0_yzzz = primBuffer.data(goff + 15 * idx + 13);

            auto g_0_zzzz = primBuffer.data(goff + 15 * idx + 14);

            #pragma omp simd aligned(fg, fz, pbx, pby, pbz, g20_0_xx, g20_0_xy,\
                                         g20_0_xz, g20_0_yy, g20_0_yz, g20_0_zz,\
                                         g21_0_xx, g21_0_xy, g21_0_xz, g21_0_yy,\
                                         g21_0_yz, g21_0_zz, g1_0_xxx, g1_0_xxy,\
                                         g1_0_xxz, g1_0_xyy, g1_0_xyz, g1_0_xzz,\
                                         g1_0_yyy, g1_0_yyz, g1_0_yzz, g1_0_zzz,\
                                         g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy,\
                                         g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz,\
                                         g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz,\
                                         g_0_yyzz, g_0_yzzz, g_0_zzzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pbx[k];

                g_0_xxxx[k] = fr * g1_0_xxx[k] + f2g * (3.0 * g20_0_xx[k] - 3.0 * fgz * g21_0_xx[k]);

                g_0_xxxy[k] = fr * g1_0_xxy[k] + f2g * (2.0 * g20_0_xy[k] - 2.0 * fgz * g21_0_xy[k]);

                g_0_xxxz[k] = fr * g1_0_xxz[k] + f2g * (2.0 * g20_0_xz[k] - 2.0 * fgz * g21_0_xz[k]);

                g_0_xxyy[k] = fr * g1_0_xyy[k] + f2g * (g20_0_yy[k] - fgz * g21_0_yy[k]);

                g_0_xxyz[k] = fr * g1_0_xyz[k] + f2g * (g20_0_yz[k] - fgz * g21_0_yz[k]);

                g_0_xxzz[k] = fr * g1_0_xzz[k] + f2g * (g20_0_zz[k] - fgz * g21_0_zz[k]);

                g_0_xyyy[k] = fr * g1_0_yyy[k];

                g_0_xyyz[k] = fr * g1_0_yyz[k];

                g_0_xyzz[k] = fr * g1_0_yzz[k];

                g_0_xzzz[k] = fr * g1_0_zzz[k];

                // leading y component

                fr = pby[k];

                g_0_yyyy[k] = fr * g1_0_yyy[k] + f2g * (3.0 * g20_0_yy[k] - 3.0 * fgz * g21_0_yy[k]);

                g_0_yyyz[k] = fr * g1_0_yyz[k] + f2g * (2.0 * g20_0_yz[k] - 2.0 * fgz * g21_0_yz[k]);

                g_0_yyzz[k] = fr * g1_0_yzz[k] + f2g * (g20_0_zz[k] - fgz * g21_0_zz[k]);

                g_0_yzzz[k] = fr * g1_0_zzz[k];

                // leading z component

                g_0_zzzz[k] = pbz[k] * g1_0_zzz[k] + f2g * (3.0 * g20_0_zz[k] - 3.0 * fgz * g21_0_zz[k]);
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForGS(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 4, 0);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {4, 0, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|g(r,r')|S)^(m) integrals

            auto g20_xx_0 = primBuffer.data(g20off + 6 * idx);

            auto g20_xy_0 = primBuffer.data(g20off + 6 * idx + 1);

            auto g20_xz_0 = primBuffer.data(g20off + 6 * idx + 2);

            auto g20_yy_0 = primBuffer.data(g20off + 6 * idx + 3);

            auto g20_yz_0 = primBuffer.data(g20off + 6 * idx + 4);

            auto g20_zz_0 = primBuffer.data(g20off + 6 * idx + 5);

            // set up pointers to (D|g(r,r')|S)^(m+1) integrals

            auto g21_xx_0 = primBuffer.data(g21off + 6 * idx);

            auto g21_xy_0 = primBuffer.data(g21off + 6 * idx + 1);

            auto g21_xz_0 = primBuffer.data(g21off + 6 * idx + 2);

            auto g21_yy_0 = primBuffer.data(g21off + 6 * idx + 3);

            auto g21_yz_0 = primBuffer.data(g21off + 6 * idx + 4);

            auto g21_zz_0 = primBuffer.data(g21off + 6 * idx + 5);

            // set up pointers to (F|g(r,r')|S)^(m+1) integrals

            auto g1_xxx_0 = primBuffer.data(g1off + 10 * idx);

            auto g1_xxy_0 = primBuffer.data(g1off + 10 * idx + 1);

            auto g1_xxz_0 = primBuffer.data(g1off + 10 * idx + 2);

            auto g1_xyy_0 = primBuffer.data(g1off + 10 * idx + 3);

            auto g1_xyz_0 = primBuffer.data(g1off + 10 * idx + 4);

            auto g1_xzz_0 = primBuffer.data(g1off + 10 * idx + 5);

            auto g1_yyy_0 = primBuffer.data(g1off + 10 * idx + 6);

            auto g1_yyz_0 = primBuffer.data(g1off + 10 * idx + 7);

            auto g1_yzz_0 = primBuffer.data(g1off + 10 * idx + 8);

            auto g1_zzz_0 = primBuffer.data(g1off + 10 * idx + 9);

            // set up pointers to (G|g(r,r')|S)^(m) integrals

            auto g_xxxx_0 = primBuffer.data(goff + 15 * idx);

            auto g_xxxy_0 = primBuffer.data(goff + 15 * idx + 1);

            auto g_xxxz_0 = primBuffer.data(goff + 15 * idx + 2);

            auto g_xxyy_0 = primBuffer.data(goff + 15 * idx + 3);

            auto g_xxyz_0 = primBuffer.data(goff + 15 * idx + 4);

            auto g_xxzz_0 = primBuffer.data(goff + 15 * idx + 5);

            auto g_xyyy_0 = primBuffer.data(goff + 15 * idx + 6);

            auto g_xyyz_0 = primBuffer.data(goff + 15 * idx + 7);

            auto g_xyzz_0 = primBuffer.data(goff + 15 * idx + 8);

            auto g_xzzz_0 = primBuffer.data(goff + 15 * idx + 9);

            auto g_yyyy_0 = primBuffer.data(goff + 15 * idx + 10);

            auto g_yyyz_0 = primBuffer.data(goff + 15 * idx + 11);

            auto g_yyzz_0 = primBuffer.data(goff + 15 * idx + 12);

            auto g_yzzz_0 = primBuffer.data(goff + 15 * idx + 13);

            auto g_zzzz_0 = primBuffer.data(goff + 15 * idx + 14);

            #pragma omp simd aligned(fg, fz, pax, pay, paz, g20_xx_0, g20_xy_0,\
                                         g20_xz_0, g20_yy_0, g20_yz_0, g20_zz_0,\
                                         g21_xx_0, g21_xy_0, g21_xz_0, g21_yy_0,\
                                         g21_yz_0, g21_zz_0, g1_xxx_0, g1_xxy_0,\
                                         g1_xxz_0, g1_xyy_0, g1_xyz_0, g1_xzz_0,\
                                         g1_yyy_0, g1_yyz_0, g1_yzz_0, g1_zzz_0,\
                                         g_xxxx_0, g_xxxy_0, g_xxxz_0, g_xxyy_0,\
                                         g_xxyz_0, g_xxzz_0, g_xyyy_0, g_xyyz_0,\
                                         g_xyzz_0, g_xzzz_0, g_yyyy_0, g_yyyz_0,\
                                         g_yyzz_0, g_yzzz_0, g_zzzz_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxxx_0[k] = fr * g1_xxx_0[k] + f2g * (3.0 * g20_xx_0[k] - 3.0 * fgz * g21_xx_0[k]);

                g_xxxy_0[k] = fr * g1_xxy_0[k] + f2g * (2.0 * g20_xy_0[k] - 2.0 * fgz * g21_xy_0[k]);

                g_xxxz_0[k] = fr * g1_xxz_0[k] + f2g * (2.0 * g20_xz_0[k] - 2.0 * fgz * g21_xz_0[k]);

                g_xxyy_0[k] = fr * g1_xyy_0[k] + f2g * (g20_yy_0[k] - fgz * g21_yy_0[k]);

                g_xxyz_0[k] = fr * g1_xyz_0[k] + f2g * (g20_yz_0[k] - fgz * g21_yz_0[k]);

                g_xxzz_0[k] = fr * g1_xzz_0[k] + f2g * (g20_zz_0[k] - fgz * g21_zz_0[k]);

                g_xyyy_0[k] = fr * g1_yyy_0[k];

                g_xyyz_0[k] = fr * g1_yyz_0[k];

                g_xyzz_0[k] = fr * g1_yzz_0[k];

                g_xzzz_0[k] = fr * g1_zzz_0[k];

                // leading y component

                fr = pay[k];

                g_yyyy_0[k] = fr * g1_yyy_0[k] + f2g * (3.0 * g20_yy_0[k] - 3.0 * fgz * g21_yy_0[k]);

                g_yyyz_0[k] = fr * g1_yyz_0[k] + f2g * (2.0 * g20_yz_0[k] - 2.0 * fgz * g21_yz_0[k]);

                g_yyzz_0[k] = fr * g1_yzz_0[k] + f2g * (g20_zz_0[k] - fgz * g21_zz_0[k]);

                g_yzzz_0[k] = fr * g1_zzz_0[k];

                // leading z component

                g_zzzz_0[k] = paz[k] * g1_zzz_0[k] + f2g * (3.0 * g20_zz_0[k] - 3.0 * fgz * g21_zz_0[k]);
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForPG(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 1, 4);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 3, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (S|g(r,r')|F)^(m+1) integrals

            auto gk_0_xxx = primBuffer.data(gkoff + 10 * idx);

            auto gk_0_xxy = primBuffer.data(gkoff + 10 * idx + 1);

            auto gk_0_xxz = primBuffer.data(gkoff + 10 * idx + 2);

            auto gk_0_xyy = primBuffer.data(gkoff + 10 * idx + 3);

            auto gk_0_xyz = primBuffer.data(gkoff + 10 * idx + 4);

            auto gk_0_xzz = primBuffer.data(gkoff + 10 * idx + 5);

            auto gk_0_yyy = primBuffer.data(gkoff + 10 * idx + 6);

            auto gk_0_yyz = primBuffer.data(gkoff + 10 * idx + 7);

            auto gk_0_yzz = primBuffer.data(gkoff + 10 * idx + 8);

            auto gk_0_zzz = primBuffer.data(gkoff + 10 * idx + 9);

            // set up pointers to (S|g(r,r')|G)^(m+1) integrals

            auto g1_0_xxxx = primBuffer.data(g1off + 15 * idx);

            auto g1_0_xxxy = primBuffer.data(g1off + 15 * idx + 1);

            auto g1_0_xxxz = primBuffer.data(g1off + 15 * idx + 2);

            auto g1_0_xxyy = primBuffer.data(g1off + 15 * idx + 3);

            auto g1_0_xxyz = primBuffer.data(g1off + 15 * idx + 4);

            auto g1_0_xxzz = primBuffer.data(g1off + 15 * idx + 5);

            auto g1_0_xyyy = primBuffer.data(g1off + 15 * idx + 6);

            auto g1_0_xyyz = primBuffer.data(g1off + 15 * idx + 7);

            auto g1_0_xyzz = primBuffer.data(g1off + 15 * idx + 8);

            auto g1_0_xzzz = primBuffer.data(g1off + 15 * idx + 9);

            auto g1_0_yyyy = primBuffer.data(g1off + 15 * idx + 10);

            auto g1_0_yyyz = primBuffer.data(g1off + 15 * idx + 11);

            auto g1_0_yyzz = primBuffer.data(g1off + 15 * idx + 12);

            auto g1_0_yzzz = primBuffer.data(g1off + 15 * idx + 13);

            auto g1_0_zzzz = primBuffer.data(g1off + 15 * idx + 14);

            // set up pointers to (P|g(r,r')|G)^(m) integrals

            auto g_x_xxxx = primBuffer.data(goff + 45 * idx);

            auto g_x_xxxy = primBuffer.data(goff + 45 * idx + 1);

            auto g_x_xxxz = primBuffer.data(goff + 45 * idx + 2);

            auto g_x_xxyy = primBuffer.data(goff + 45 * idx + 3);

            auto g_x_xxyz = primBuffer.data(goff + 45 * idx + 4);

            auto g_x_xxzz = primBuffer.data(goff + 45 * idx + 5);

            auto g_x_xyyy = primBuffer.data(goff + 45 * idx + 6);

            auto g_x_xyyz = primBuffer.data(goff + 45 * idx + 7);

            auto g_x_xyzz = primBuffer.data(goff + 45 * idx + 8);

            auto g_x_xzzz = primBuffer.data(goff + 45 * idx + 9);

            auto g_x_yyyy = primBuffer.data(goff + 45 * idx + 10);

            auto g_x_yyyz = primBuffer.data(goff + 45 * idx + 11);

            auto g_x_yyzz = primBuffer.data(goff + 45 * idx + 12);

            auto g_x_yzzz = primBuffer.data(goff + 45 * idx + 13);

            auto g_x_zzzz = primBuffer.data(goff + 45 * idx + 14);

            auto g_y_xxxx = primBuffer.data(goff + 45 * idx + 15);

            auto g_y_xxxy = primBuffer.data(goff + 45 * idx + 16);

            auto g_y_xxxz = primBuffer.data(goff + 45 * idx + 17);

            auto g_y_xxyy = primBuffer.data(goff + 45 * idx + 18);

            auto g_y_xxyz = primBuffer.data(goff + 45 * idx + 19);

            auto g_y_xxzz = primBuffer.data(goff + 45 * idx + 20);

            auto g_y_xyyy = primBuffer.data(goff + 45 * idx + 21);

            auto g_y_xyyz = primBuffer.data(goff + 45 * idx + 22);

            auto g_y_xyzz = primBuffer.data(goff + 45 * idx + 23);

            auto g_y_xzzz = primBuffer.data(goff + 45 * idx + 24);

            auto g_y_yyyy = primBuffer.data(goff + 45 * idx + 25);

            auto g_y_yyyz = primBuffer.data(goff + 45 * idx + 26);

            auto g_y_yyzz = primBuffer.data(goff + 45 * idx + 27);

            auto g_y_yzzz = primBuffer.data(goff + 45 * idx + 28);

            auto g_y_zzzz = primBuffer.data(goff + 45 * idx + 29);

            auto g_z_xxxx = primBuffer.data(goff + 45 * idx + 30);

            auto g_z_xxxy = primBuffer.data(goff + 45 * idx + 31);

            auto g_z_xxxz = primBuffer.data(goff + 45 * idx + 32);

            auto g_z_xxyy = primBuffer.data(goff + 45 * idx + 33);

            auto g_z_xxyz = primBuffer.data(goff + 45 * idx + 34);

            auto g_z_xxzz = primBuffer.data(goff + 45 * idx + 35);

            auto g_z_xyyy = primBuffer.data(goff + 45 * idx + 36);

            auto g_z_xyyz = primBuffer.data(goff + 45 * idx + 37);

            auto g_z_xyzz = primBuffer.data(goff + 45 * idx + 38);

            auto g_z_xzzz = primBuffer.data(goff + 45 * idx + 39);

            auto g_z_yyyy = primBuffer.data(goff + 45 * idx + 40);

            auto g_z_yyyz = primBuffer.data(goff + 45 * idx + 41);

            auto g_z_yyzz = primBuffer.data(goff + 45 * idx + 42);

            auto g_z_yzzz = primBuffer.data(goff + 45 * idx + 43);

            auto g_z_zzzz = primBuffer.data(goff + 45 * idx + 44);

            #pragma omp simd aligned(fx, pax, pay, paz, gk_0_xxx, gk_0_xxy,\
                                         gk_0_xxz, gk_0_xyy, gk_0_xyz, gk_0_xzz,\
                                         gk_0_yyy, gk_0_yyz, gk_0_yzz, gk_0_zzz,\
                                         g1_0_xxxx, g1_0_xxxy, g1_0_xxxz, g1_0_xxyy,\
                                         g1_0_xxyz, g1_0_xxzz, g1_0_xyyy, g1_0_xyyz,\
                                         g1_0_xyzz, g1_0_xzzz, g1_0_yyyy, g1_0_yyyz,\
                                         g1_0_yyzz, g1_0_yzzz, g1_0_zzzz, g_x_xxxx,\
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
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // leading x component

                double fr = pax[k];

                g_x_xxxx[k] = fr * g1_0_xxxx[k] + 4.0 * f2t * gk_0_xxx[k];

                g_x_xxxy[k] = fr * g1_0_xxxy[k] + 3.0 * f2t * gk_0_xxy[k];

                g_x_xxxz[k] = fr * g1_0_xxxz[k] + 3.0 * f2t * gk_0_xxz[k];

                g_x_xxyy[k] = fr * g1_0_xxyy[k] + 2.0 * f2t * gk_0_xyy[k];

                g_x_xxyz[k] = fr * g1_0_xxyz[k] + 2.0 * f2t * gk_0_xyz[k];

                g_x_xxzz[k] = fr * g1_0_xxzz[k] + 2.0 * f2t * gk_0_xzz[k];

                g_x_xyyy[k] = fr * g1_0_xyyy[k] + f2t * gk_0_yyy[k];

                g_x_xyyz[k] = fr * g1_0_xyyz[k] + f2t * gk_0_yyz[k];

                g_x_xyzz[k] = fr * g1_0_xyzz[k] + f2t * gk_0_yzz[k];

                g_x_xzzz[k] = fr * g1_0_xzzz[k] + f2t * gk_0_zzz[k];

                g_x_yyyy[k] = fr * g1_0_yyyy[k];

                g_x_yyyz[k] = fr * g1_0_yyyz[k];

                g_x_yyzz[k] = fr * g1_0_yyzz[k];

                g_x_yzzz[k] = fr * g1_0_yzzz[k];

                g_x_zzzz[k] = fr * g1_0_zzzz[k];

                // leading y component

                fr = pay[k];

                g_y_xxxx[k] = fr * g1_0_xxxx[k];

                g_y_xxxy[k] = fr * g1_0_xxxy[k] + f2t * gk_0_xxx[k];

                g_y_xxxz[k] = fr * g1_0_xxxz[k];

                g_y_xxyy[k] = fr * g1_0_xxyy[k] + 2.0 * f2t * gk_0_xxy[k];

                g_y_xxyz[k] = fr * g1_0_xxyz[k] + f2t * gk_0_xxz[k];

                g_y_xxzz[k] = fr * g1_0_xxzz[k];

                g_y_xyyy[k] = fr * g1_0_xyyy[k] + 3.0 * f2t * gk_0_xyy[k];

                g_y_xyyz[k] = fr * g1_0_xyyz[k] + 2.0 * f2t * gk_0_xyz[k];

                g_y_xyzz[k] = fr * g1_0_xyzz[k] + f2t * gk_0_xzz[k];

                g_y_xzzz[k] = fr * g1_0_xzzz[k];

                g_y_yyyy[k] = fr * g1_0_yyyy[k] + 4.0 * f2t * gk_0_yyy[k];

                g_y_yyyz[k] = fr * g1_0_yyyz[k] + 3.0 * f2t * gk_0_yyz[k];

                g_y_yyzz[k] = fr * g1_0_yyzz[k] + 2.0 * f2t * gk_0_yzz[k];

                g_y_yzzz[k] = fr * g1_0_yzzz[k] + f2t * gk_0_zzz[k];

                g_y_zzzz[k] = fr * g1_0_zzzz[k];

                // leading z component

                fr = paz[k];

                g_z_xxxx[k] = fr * g1_0_xxxx[k];

                g_z_xxxy[k] = fr * g1_0_xxxy[k];

                g_z_xxxz[k] = fr * g1_0_xxxz[k] + f2t * gk_0_xxx[k];

                g_z_xxyy[k] = fr * g1_0_xxyy[k];

                g_z_xxyz[k] = fr * g1_0_xxyz[k] + f2t * gk_0_xxy[k];

                g_z_xxzz[k] = fr * g1_0_xxzz[k] + 2.0 * f2t * gk_0_xxz[k];

                g_z_xyyy[k] = fr * g1_0_xyyy[k];

                g_z_xyyz[k] = fr * g1_0_xyyz[k] + f2t * gk_0_xyy[k];

                g_z_xyzz[k] = fr * g1_0_xyzz[k] + 2.0 * f2t * gk_0_xyz[k];

                g_z_xzzz[k] = fr * g1_0_xzzz[k] + 3.0 * f2t * gk_0_xzz[k];

                g_z_yyyy[k] = fr * g1_0_yyyy[k];

                g_z_yyyz[k] = fr * g1_0_yyyz[k] + f2t * gk_0_yyy[k];

                g_z_yyzz[k] = fr * g1_0_yyzz[k] + 2.0 * f2t * gk_0_yyz[k];

                g_z_yzzz[k] = fr * g1_0_yzzz[k] + 3.0 * f2t * gk_0_yzz[k];

                g_z_zzzz[k] = fr * g1_0_zzzz[k] + 4.0 * f2t * gk_0_zzz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForGP(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 4, 1);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {4, 1, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 1, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 0, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|g(r,r')|S)^(m+1) integrals

            auto gk_xxx_0 = primBuffer.data(gkoff + 10 * idx);

            auto gk_xxy_0 = primBuffer.data(gkoff + 10 * idx + 1);

            auto gk_xxz_0 = primBuffer.data(gkoff + 10 * idx + 2);

            auto gk_xyy_0 = primBuffer.data(gkoff + 10 * idx + 3);

            auto gk_xyz_0 = primBuffer.data(gkoff + 10 * idx + 4);

            auto gk_xzz_0 = primBuffer.data(gkoff + 10 * idx + 5);

            auto gk_yyy_0 = primBuffer.data(gkoff + 10 * idx + 6);

            auto gk_yyz_0 = primBuffer.data(gkoff + 10 * idx + 7);

            auto gk_yzz_0 = primBuffer.data(gkoff + 10 * idx + 8);

            auto gk_zzz_0 = primBuffer.data(gkoff + 10 * idx + 9);

            // set up pointers to (D|g(r,r')|P)^(m) integrals

            auto g20_xx_x = primBuffer.data(g20off + 18 * idx);

            auto g20_xx_y = primBuffer.data(g20off + 18 * idx + 1);

            auto g20_xx_z = primBuffer.data(g20off + 18 * idx + 2);

            auto g20_xy_x = primBuffer.data(g20off + 18 * idx + 3);

            auto g20_xy_y = primBuffer.data(g20off + 18 * idx + 4);

            auto g20_xy_z = primBuffer.data(g20off + 18 * idx + 5);

            auto g20_xz_x = primBuffer.data(g20off + 18 * idx + 6);

            auto g20_xz_y = primBuffer.data(g20off + 18 * idx + 7);

            auto g20_xz_z = primBuffer.data(g20off + 18 * idx + 8);

            auto g20_yy_x = primBuffer.data(g20off + 18 * idx + 9);

            auto g20_yy_y = primBuffer.data(g20off + 18 * idx + 10);

            auto g20_yy_z = primBuffer.data(g20off + 18 * idx + 11);

            auto g20_yz_x = primBuffer.data(g20off + 18 * idx + 12);

            auto g20_yz_y = primBuffer.data(g20off + 18 * idx + 13);

            auto g20_yz_z = primBuffer.data(g20off + 18 * idx + 14);

            auto g20_zz_x = primBuffer.data(g20off + 18 * idx + 15);

            auto g20_zz_y = primBuffer.data(g20off + 18 * idx + 16);

            auto g20_zz_z = primBuffer.data(g20off + 18 * idx + 17);

            // set up pointers to (D|g(r,r')|P)^(m+1) integrals

            auto g21_xx_x = primBuffer.data(g21off + 18 * idx);

            auto g21_xx_y = primBuffer.data(g21off + 18 * idx + 1);

            auto g21_xx_z = primBuffer.data(g21off + 18 * idx + 2);

            auto g21_xy_x = primBuffer.data(g21off + 18 * idx + 3);

            auto g21_xy_y = primBuffer.data(g21off + 18 * idx + 4);

            auto g21_xy_z = primBuffer.data(g21off + 18 * idx + 5);

            auto g21_xz_x = primBuffer.data(g21off + 18 * idx + 6);

            auto g21_xz_y = primBuffer.data(g21off + 18 * idx + 7);

            auto g21_xz_z = primBuffer.data(g21off + 18 * idx + 8);

            auto g21_yy_x = primBuffer.data(g21off + 18 * idx + 9);

            auto g21_yy_y = primBuffer.data(g21off + 18 * idx + 10);

            auto g21_yy_z = primBuffer.data(g21off + 18 * idx + 11);

            auto g21_yz_x = primBuffer.data(g21off + 18 * idx + 12);

            auto g21_yz_y = primBuffer.data(g21off + 18 * idx + 13);

            auto g21_yz_z = primBuffer.data(g21off + 18 * idx + 14);

            auto g21_zz_x = primBuffer.data(g21off + 18 * idx + 15);

            auto g21_zz_y = primBuffer.data(g21off + 18 * idx + 16);

            auto g21_zz_z = primBuffer.data(g21off + 18 * idx + 17);

            // set up pointers to (F|g(r,r')|P)^(m+1) integrals

            auto g1_xxx_x = primBuffer.data(g1off + 30 * idx);

            auto g1_xxx_y = primBuffer.data(g1off + 30 * idx + 1);

            auto g1_xxx_z = primBuffer.data(g1off + 30 * idx + 2);

            auto g1_xxy_x = primBuffer.data(g1off + 30 * idx + 3);

            auto g1_xxy_y = primBuffer.data(g1off + 30 * idx + 4);

            auto g1_xxy_z = primBuffer.data(g1off + 30 * idx + 5);

            auto g1_xxz_x = primBuffer.data(g1off + 30 * idx + 6);

            auto g1_xxz_y = primBuffer.data(g1off + 30 * idx + 7);

            auto g1_xxz_z = primBuffer.data(g1off + 30 * idx + 8);

            auto g1_xyy_x = primBuffer.data(g1off + 30 * idx + 9);

            auto g1_xyy_y = primBuffer.data(g1off + 30 * idx + 10);

            auto g1_xyy_z = primBuffer.data(g1off + 30 * idx + 11);

            auto g1_xyz_x = primBuffer.data(g1off + 30 * idx + 12);

            auto g1_xyz_y = primBuffer.data(g1off + 30 * idx + 13);

            auto g1_xyz_z = primBuffer.data(g1off + 30 * idx + 14);

            auto g1_xzz_x = primBuffer.data(g1off + 30 * idx + 15);

            auto g1_xzz_y = primBuffer.data(g1off + 30 * idx + 16);

            auto g1_xzz_z = primBuffer.data(g1off + 30 * idx + 17);

            auto g1_yyy_x = primBuffer.data(g1off + 30 * idx + 18);

            auto g1_yyy_y = primBuffer.data(g1off + 30 * idx + 19);

            auto g1_yyy_z = primBuffer.data(g1off + 30 * idx + 20);

            auto g1_yyz_x = primBuffer.data(g1off + 30 * idx + 21);

            auto g1_yyz_y = primBuffer.data(g1off + 30 * idx + 22);

            auto g1_yyz_z = primBuffer.data(g1off + 30 * idx + 23);

            auto g1_yzz_x = primBuffer.data(g1off + 30 * idx + 24);

            auto g1_yzz_y = primBuffer.data(g1off + 30 * idx + 25);

            auto g1_yzz_z = primBuffer.data(g1off + 30 * idx + 26);

            auto g1_zzz_x = primBuffer.data(g1off + 30 * idx + 27);

            auto g1_zzz_y = primBuffer.data(g1off + 30 * idx + 28);

            auto g1_zzz_z = primBuffer.data(g1off + 30 * idx + 29);

            // set up pointers to (G|g(r,r')|P)^(m) integrals

            auto g_xxxx_x = primBuffer.data(goff + 45 * idx);

            auto g_xxxx_y = primBuffer.data(goff + 45 * idx + 1);

            auto g_xxxx_z = primBuffer.data(goff + 45 * idx + 2);

            auto g_xxxy_x = primBuffer.data(goff + 45 * idx + 3);

            auto g_xxxy_y = primBuffer.data(goff + 45 * idx + 4);

            auto g_xxxy_z = primBuffer.data(goff + 45 * idx + 5);

            auto g_xxxz_x = primBuffer.data(goff + 45 * idx + 6);

            auto g_xxxz_y = primBuffer.data(goff + 45 * idx + 7);

            auto g_xxxz_z = primBuffer.data(goff + 45 * idx + 8);

            auto g_xxyy_x = primBuffer.data(goff + 45 * idx + 9);

            auto g_xxyy_y = primBuffer.data(goff + 45 * idx + 10);

            auto g_xxyy_z = primBuffer.data(goff + 45 * idx + 11);

            auto g_xxyz_x = primBuffer.data(goff + 45 * idx + 12);

            auto g_xxyz_y = primBuffer.data(goff + 45 * idx + 13);

            auto g_xxyz_z = primBuffer.data(goff + 45 * idx + 14);

            auto g_xxzz_x = primBuffer.data(goff + 45 * idx + 15);

            auto g_xxzz_y = primBuffer.data(goff + 45 * idx + 16);

            auto g_xxzz_z = primBuffer.data(goff + 45 * idx + 17);

            auto g_xyyy_x = primBuffer.data(goff + 45 * idx + 18);

            auto g_xyyy_y = primBuffer.data(goff + 45 * idx + 19);

            auto g_xyyy_z = primBuffer.data(goff + 45 * idx + 20);

            auto g_xyyz_x = primBuffer.data(goff + 45 * idx + 21);

            auto g_xyyz_y = primBuffer.data(goff + 45 * idx + 22);

            auto g_xyyz_z = primBuffer.data(goff + 45 * idx + 23);

            auto g_xyzz_x = primBuffer.data(goff + 45 * idx + 24);

            auto g_xyzz_y = primBuffer.data(goff + 45 * idx + 25);

            auto g_xyzz_z = primBuffer.data(goff + 45 * idx + 26);

            auto g_xzzz_x = primBuffer.data(goff + 45 * idx + 27);

            auto g_xzzz_y = primBuffer.data(goff + 45 * idx + 28);

            auto g_xzzz_z = primBuffer.data(goff + 45 * idx + 29);

            auto g_yyyy_x = primBuffer.data(goff + 45 * idx + 30);

            auto g_yyyy_y = primBuffer.data(goff + 45 * idx + 31);

            auto g_yyyy_z = primBuffer.data(goff + 45 * idx + 32);

            auto g_yyyz_x = primBuffer.data(goff + 45 * idx + 33);

            auto g_yyyz_y = primBuffer.data(goff + 45 * idx + 34);

            auto g_yyyz_z = primBuffer.data(goff + 45 * idx + 35);

            auto g_yyzz_x = primBuffer.data(goff + 45 * idx + 36);

            auto g_yyzz_y = primBuffer.data(goff + 45 * idx + 37);

            auto g_yyzz_z = primBuffer.data(goff + 45 * idx + 38);

            auto g_yzzz_x = primBuffer.data(goff + 45 * idx + 39);

            auto g_yzzz_y = primBuffer.data(goff + 45 * idx + 40);

            auto g_yzzz_z = primBuffer.data(goff + 45 * idx + 41);

            auto g_zzzz_x = primBuffer.data(goff + 45 * idx + 42);

            auto g_zzzz_y = primBuffer.data(goff + 45 * idx + 43);

            auto g_zzzz_z = primBuffer.data(goff + 45 * idx + 44);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xxx_0, gk_xxy_0,\
                                         gk_xxz_0, gk_xyy_0, gk_xyz_0, gk_xzz_0,\
                                         gk_yyy_0, gk_yyz_0, gk_yzz_0, gk_zzz_0,\
                                         g20_xx_x, g20_xx_y, g20_xx_z, g20_xy_x,\
                                         g20_xy_y, g20_xy_z, g20_xz_x, g20_xz_y,\
                                         g20_xz_z, g20_yy_x, g20_yy_y, g20_yy_z,\
                                         g20_yz_x, g20_yz_y, g20_yz_z, g20_zz_x,\
                                         g20_zz_y, g20_zz_z, g21_xx_x, g21_xx_y,\
                                         g21_xx_z, g21_xy_x, g21_xy_y, g21_xy_z,\
                                         g21_xz_x, g21_xz_y, g21_xz_z, g21_yy_x,\
                                         g21_yy_y, g21_yy_z, g21_yz_x, g21_yz_y,\
                                         g21_yz_z, g21_zz_x, g21_zz_y, g21_zz_z,\
                                         g1_xxx_x, g1_xxx_y, g1_xxx_z, g1_xxy_x,\
                                         g1_xxy_y, g1_xxy_z, g1_xxz_x, g1_xxz_y,\
                                         g1_xxz_z, g1_xyy_x, g1_xyy_y, g1_xyy_z,\
                                         g1_xyz_x, g1_xyz_y, g1_xyz_z, g1_xzz_x,\
                                         g1_xzz_y, g1_xzz_z, g1_yyy_x, g1_yyy_y,\
                                         g1_yyy_z, g1_yyz_x, g1_yyz_y, g1_yyz_z,\
                                         g1_yzz_x, g1_yzz_y, g1_yzz_z, g1_zzz_x,\
                                         g1_zzz_y, g1_zzz_z, g_xxxx_x, g_xxxx_y,\
                                         g_xxxx_z, g_xxxy_x, g_xxxy_y, g_xxxy_z,\
                                         g_xxxz_x, g_xxxz_y, g_xxxz_z, g_xxyy_x,\
                                         g_xxyy_y, g_xxyy_z, g_xxyz_x, g_xxyz_y,\
                                         g_xxyz_z, g_xxzz_x, g_xxzz_y, g_xxzz_z,\
                                         g_xyyy_x, g_xyyy_y, g_xyyy_z, g_xyyz_x,\
                                         g_xyyz_y, g_xyyz_z, g_xyzz_x, g_xyzz_y,\
                                         g_xyzz_z, g_xzzz_x, g_xzzz_y, g_xzzz_z,\
                                         g_yyyy_x, g_yyyy_y, g_yyyy_z, g_yyyz_x,\
                                         g_yyyz_y, g_yyyz_z, g_yyzz_x, g_yyzz_y,\
                                         g_yyzz_z, g_yzzz_x, g_yzzz_y, g_yzzz_z,\
                                         g_zzzz_x, g_zzzz_y, g_zzzz_z: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxxx_x[k] = fr * g1_xxx_x[k] + f2g * (3.0 * g20_xx_x[k] - 3.0 * fgz * g21_xx_x[k]) + f2t * gk_xxx_0[k];

                g_xxxx_y[k] = fr * g1_xxx_y[k] + f2g * (3.0 * g20_xx_y[k] - 3.0 * fgz * g21_xx_y[k]);

                g_xxxx_z[k] = fr * g1_xxx_z[k] + f2g * (3.0 * g20_xx_z[k] - 3.0 * fgz * g21_xx_z[k]);

                g_xxxy_x[k] = fr * g1_xxy_x[k] + f2g * (2.0 * g20_xy_x[k] - 2.0 * fgz * g21_xy_x[k]) + f2t * gk_xxy_0[k];

                g_xxxy_y[k] = fr * g1_xxy_y[k] + f2g * (2.0 * g20_xy_y[k] - 2.0 * fgz * g21_xy_y[k]);

                g_xxxy_z[k] = fr * g1_xxy_z[k] + f2g * (2.0 * g20_xy_z[k] - 2.0 * fgz * g21_xy_z[k]);

                g_xxxz_x[k] = fr * g1_xxz_x[k] + f2g * (2.0 * g20_xz_x[k] - 2.0 * fgz * g21_xz_x[k]) + f2t * gk_xxz_0[k];

                g_xxxz_y[k] = fr * g1_xxz_y[k] + f2g * (2.0 * g20_xz_y[k] - 2.0 * fgz * g21_xz_y[k]);

                g_xxxz_z[k] = fr * g1_xxz_z[k] + f2g * (2.0 * g20_xz_z[k] - 2.0 * fgz * g21_xz_z[k]);

                g_xxyy_x[k] = fr * g1_xyy_x[k] + f2g * (g20_yy_x[k] - fgz * g21_yy_x[k]) + f2t * gk_xyy_0[k];

                g_xxyy_y[k] = fr * g1_xyy_y[k] + f2g * (g20_yy_y[k] - fgz * g21_yy_y[k]);

                g_xxyy_z[k] = fr * g1_xyy_z[k] + f2g * (g20_yy_z[k] - fgz * g21_yy_z[k]);

                g_xxyz_x[k] = fr * g1_xyz_x[k] + f2g * (g20_yz_x[k] - fgz * g21_yz_x[k]) + f2t * gk_xyz_0[k];

                g_xxyz_y[k] = fr * g1_xyz_y[k] + f2g * (g20_yz_y[k] - fgz * g21_yz_y[k]);

                g_xxyz_z[k] = fr * g1_xyz_z[k] + f2g * (g20_yz_z[k] - fgz * g21_yz_z[k]);

                g_xxzz_x[k] = fr * g1_xzz_x[k] + f2g * (g20_zz_x[k] - fgz * g21_zz_x[k]) + f2t * gk_xzz_0[k];

                g_xxzz_y[k] = fr * g1_xzz_y[k] + f2g * (g20_zz_y[k] - fgz * g21_zz_y[k]);

                g_xxzz_z[k] = fr * g1_xzz_z[k] + f2g * (g20_zz_z[k] - fgz * g21_zz_z[k]);

                g_xyyy_x[k] = fr * g1_yyy_x[k] + f2t * gk_yyy_0[k];

                g_xyyy_y[k] = fr * g1_yyy_y[k];

                g_xyyy_z[k] = fr * g1_yyy_z[k];

                g_xyyz_x[k] = fr * g1_yyz_x[k] + f2t * gk_yyz_0[k];

                g_xyyz_y[k] = fr * g1_yyz_y[k];

                g_xyyz_z[k] = fr * g1_yyz_z[k];

                g_xyzz_x[k] = fr * g1_yzz_x[k] + f2t * gk_yzz_0[k];

                g_xyzz_y[k] = fr * g1_yzz_y[k];

                g_xyzz_z[k] = fr * g1_yzz_z[k];

                g_xzzz_x[k] = fr * g1_zzz_x[k] + f2t * gk_zzz_0[k];

                g_xzzz_y[k] = fr * g1_zzz_y[k];

                g_xzzz_z[k] = fr * g1_zzz_z[k];

                // leading y component

                fr = pay[k];

                g_yyyy_x[k] = fr * g1_yyy_x[k] + f2g * (3.0 * g20_yy_x[k] - 3.0 * fgz * g21_yy_x[k]);

                g_yyyy_y[k] = fr * g1_yyy_y[k] + f2g * (3.0 * g20_yy_y[k] - 3.0 * fgz * g21_yy_y[k]) + f2t * gk_yyy_0[k];

                g_yyyy_z[k] = fr * g1_yyy_z[k] + f2g * (3.0 * g20_yy_z[k] - 3.0 * fgz * g21_yy_z[k]);

                g_yyyz_x[k] = fr * g1_yyz_x[k] + f2g * (2.0 * g20_yz_x[k] - 2.0 * fgz * g21_yz_x[k]);

                g_yyyz_y[k] = fr * g1_yyz_y[k] + f2g * (2.0 * g20_yz_y[k] - 2.0 * fgz * g21_yz_y[k]) + f2t * gk_yyz_0[k];

                g_yyyz_z[k] = fr * g1_yyz_z[k] + f2g * (2.0 * g20_yz_z[k] - 2.0 * fgz * g21_yz_z[k]);

                g_yyzz_x[k] = fr * g1_yzz_x[k] + f2g * (g20_zz_x[k] - fgz * g21_zz_x[k]);

                g_yyzz_y[k] = fr * g1_yzz_y[k] + f2g * (g20_zz_y[k] - fgz * g21_zz_y[k]) + f2t * gk_yzz_0[k];

                g_yyzz_z[k] = fr * g1_yzz_z[k] + f2g * (g20_zz_z[k] - fgz * g21_zz_z[k]);

                g_yzzz_x[k] = fr * g1_zzz_x[k];

                g_yzzz_y[k] = fr * g1_zzz_y[k] + f2t * gk_zzz_0[k];

                g_yzzz_z[k] = fr * g1_zzz_z[k];

                // leading z component

                fr = paz[k];

                g_zzzz_x[k] = fr * g1_zzz_x[k] + f2g * (3.0 * g20_zz_x[k] - 3.0 * fgz * g21_zz_x[k]);

                g_zzzz_y[k] = fr * g1_zzz_y[k] + f2g * (3.0 * g20_zz_y[k] - 3.0 * fgz * g21_zz_y[k]);

                g_zzzz_z[k] = fr * g1_zzz_z[k] + f2g * (3.0 * g20_zz_z[k] - 3.0 * fgz * g21_zz_z[k]) + f2t * gk_zzz_0[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForDG(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 2, 4);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {0, 4, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {1, 3, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (P|g(r,r')|F)^(m+1) integrals

            auto gk_x_xxx = primBuffer.data(gkoff + 30 * idx);

            auto gk_x_xxy = primBuffer.data(gkoff + 30 * idx + 1);

            auto gk_x_xxz = primBuffer.data(gkoff + 30 * idx + 2);

            auto gk_x_xyy = primBuffer.data(gkoff + 30 * idx + 3);

            auto gk_x_xyz = primBuffer.data(gkoff + 30 * idx + 4);

            auto gk_x_xzz = primBuffer.data(gkoff + 30 * idx + 5);

            auto gk_x_yyy = primBuffer.data(gkoff + 30 * idx + 6);

            auto gk_x_yyz = primBuffer.data(gkoff + 30 * idx + 7);

            auto gk_x_yzz = primBuffer.data(gkoff + 30 * idx + 8);

            auto gk_x_zzz = primBuffer.data(gkoff + 30 * idx + 9);

            auto gk_y_xxx = primBuffer.data(gkoff + 30 * idx + 10);

            auto gk_y_xxy = primBuffer.data(gkoff + 30 * idx + 11);

            auto gk_y_xxz = primBuffer.data(gkoff + 30 * idx + 12);

            auto gk_y_xyy = primBuffer.data(gkoff + 30 * idx + 13);

            auto gk_y_xyz = primBuffer.data(gkoff + 30 * idx + 14);

            auto gk_y_xzz = primBuffer.data(gkoff + 30 * idx + 15);

            auto gk_y_yyy = primBuffer.data(gkoff + 30 * idx + 16);

            auto gk_y_yyz = primBuffer.data(gkoff + 30 * idx + 17);

            auto gk_y_yzz = primBuffer.data(gkoff + 30 * idx + 18);

            auto gk_y_zzz = primBuffer.data(gkoff + 30 * idx + 19);

            auto gk_z_xxx = primBuffer.data(gkoff + 30 * idx + 20);

            auto gk_z_xxy = primBuffer.data(gkoff + 30 * idx + 21);

            auto gk_z_xxz = primBuffer.data(gkoff + 30 * idx + 22);

            auto gk_z_xyy = primBuffer.data(gkoff + 30 * idx + 23);

            auto gk_z_xyz = primBuffer.data(gkoff + 30 * idx + 24);

            auto gk_z_xzz = primBuffer.data(gkoff + 30 * idx + 25);

            auto gk_z_yyy = primBuffer.data(gkoff + 30 * idx + 26);

            auto gk_z_yyz = primBuffer.data(gkoff + 30 * idx + 27);

            auto gk_z_yzz = primBuffer.data(gkoff + 30 * idx + 28);

            auto gk_z_zzz = primBuffer.data(gkoff + 30 * idx + 29);

            // set up pointers to (S|g(r,r')|G)^(m) integrals

            auto g20_0_xxxx = primBuffer.data(g20off + 15 * idx);

            auto g20_0_xxxy = primBuffer.data(g20off + 15 * idx + 1);

            auto g20_0_xxxz = primBuffer.data(g20off + 15 * idx + 2);

            auto g20_0_xxyy = primBuffer.data(g20off + 15 * idx + 3);

            auto g20_0_xxyz = primBuffer.data(g20off + 15 * idx + 4);

            auto g20_0_xxzz = primBuffer.data(g20off + 15 * idx + 5);

            auto g20_0_xyyy = primBuffer.data(g20off + 15 * idx + 6);

            auto g20_0_xyyz = primBuffer.data(g20off + 15 * idx + 7);

            auto g20_0_xyzz = primBuffer.data(g20off + 15 * idx + 8);

            auto g20_0_xzzz = primBuffer.data(g20off + 15 * idx + 9);

            auto g20_0_yyyy = primBuffer.data(g20off + 15 * idx + 10);

            auto g20_0_yyyz = primBuffer.data(g20off + 15 * idx + 11);

            auto g20_0_yyzz = primBuffer.data(g20off + 15 * idx + 12);

            auto g20_0_yzzz = primBuffer.data(g20off + 15 * idx + 13);

            auto g20_0_zzzz = primBuffer.data(g20off + 15 * idx + 14);

            // set up pointers to (S|g(r,r')|G)^(m+1) integrals

            auto g21_0_xxxx = primBuffer.data(g21off + 15 * idx);

            auto g21_0_xxxy = primBuffer.data(g21off + 15 * idx + 1);

            auto g21_0_xxxz = primBuffer.data(g21off + 15 * idx + 2);

            auto g21_0_xxyy = primBuffer.data(g21off + 15 * idx + 3);

            auto g21_0_xxyz = primBuffer.data(g21off + 15 * idx + 4);

            auto g21_0_xxzz = primBuffer.data(g21off + 15 * idx + 5);

            auto g21_0_xyyy = primBuffer.data(g21off + 15 * idx + 6);

            auto g21_0_xyyz = primBuffer.data(g21off + 15 * idx + 7);

            auto g21_0_xyzz = primBuffer.data(g21off + 15 * idx + 8);

            auto g21_0_xzzz = primBuffer.data(g21off + 15 * idx + 9);

            auto g21_0_yyyy = primBuffer.data(g21off + 15 * idx + 10);

            auto g21_0_yyyz = primBuffer.data(g21off + 15 * idx + 11);

            auto g21_0_yyzz = primBuffer.data(g21off + 15 * idx + 12);

            auto g21_0_yzzz = primBuffer.data(g21off + 15 * idx + 13);

            auto g21_0_zzzz = primBuffer.data(g21off + 15 * idx + 14);

            // set up pointers to (P|g(r,r')|G)^(m+1) integrals

            auto g1_x_xxxx = primBuffer.data(g1off + 45 * idx);

            auto g1_x_xxxy = primBuffer.data(g1off + 45 * idx + 1);

            auto g1_x_xxxz = primBuffer.data(g1off + 45 * idx + 2);

            auto g1_x_xxyy = primBuffer.data(g1off + 45 * idx + 3);

            auto g1_x_xxyz = primBuffer.data(g1off + 45 * idx + 4);

            auto g1_x_xxzz = primBuffer.data(g1off + 45 * idx + 5);

            auto g1_x_xyyy = primBuffer.data(g1off + 45 * idx + 6);

            auto g1_x_xyyz = primBuffer.data(g1off + 45 * idx + 7);

            auto g1_x_xyzz = primBuffer.data(g1off + 45 * idx + 8);

            auto g1_x_xzzz = primBuffer.data(g1off + 45 * idx + 9);

            auto g1_x_yyyy = primBuffer.data(g1off + 45 * idx + 10);

            auto g1_x_yyyz = primBuffer.data(g1off + 45 * idx + 11);

            auto g1_x_yyzz = primBuffer.data(g1off + 45 * idx + 12);

            auto g1_x_yzzz = primBuffer.data(g1off + 45 * idx + 13);

            auto g1_x_zzzz = primBuffer.data(g1off + 45 * idx + 14);

            auto g1_y_xxxx = primBuffer.data(g1off + 45 * idx + 15);

            auto g1_y_xxxy = primBuffer.data(g1off + 45 * idx + 16);

            auto g1_y_xxxz = primBuffer.data(g1off + 45 * idx + 17);

            auto g1_y_xxyy = primBuffer.data(g1off + 45 * idx + 18);

            auto g1_y_xxyz = primBuffer.data(g1off + 45 * idx + 19);

            auto g1_y_xxzz = primBuffer.data(g1off + 45 * idx + 20);

            auto g1_y_xyyy = primBuffer.data(g1off + 45 * idx + 21);

            auto g1_y_xyyz = primBuffer.data(g1off + 45 * idx + 22);

            auto g1_y_xyzz = primBuffer.data(g1off + 45 * idx + 23);

            auto g1_y_xzzz = primBuffer.data(g1off + 45 * idx + 24);

            auto g1_y_yyyy = primBuffer.data(g1off + 45 * idx + 25);

            auto g1_y_yyyz = primBuffer.data(g1off + 45 * idx + 26);

            auto g1_y_yyzz = primBuffer.data(g1off + 45 * idx + 27);

            auto g1_y_yzzz = primBuffer.data(g1off + 45 * idx + 28);

            auto g1_y_zzzz = primBuffer.data(g1off + 45 * idx + 29);

            auto g1_z_xxxx = primBuffer.data(g1off + 45 * idx + 30);

            auto g1_z_xxxy = primBuffer.data(g1off + 45 * idx + 31);

            auto g1_z_xxxz = primBuffer.data(g1off + 45 * idx + 32);

            auto g1_z_xxyy = primBuffer.data(g1off + 45 * idx + 33);

            auto g1_z_xxyz = primBuffer.data(g1off + 45 * idx + 34);

            auto g1_z_xxzz = primBuffer.data(g1off + 45 * idx + 35);

            auto g1_z_xyyy = primBuffer.data(g1off + 45 * idx + 36);

            auto g1_z_xyyz = primBuffer.data(g1off + 45 * idx + 37);

            auto g1_z_xyzz = primBuffer.data(g1off + 45 * idx + 38);

            auto g1_z_xzzz = primBuffer.data(g1off + 45 * idx + 39);

            auto g1_z_yyyy = primBuffer.data(g1off + 45 * idx + 40);

            auto g1_z_yyyz = primBuffer.data(g1off + 45 * idx + 41);

            auto g1_z_yyzz = primBuffer.data(g1off + 45 * idx + 42);

            auto g1_z_yzzz = primBuffer.data(g1off + 45 * idx + 43);

            auto g1_z_zzzz = primBuffer.data(g1off + 45 * idx + 44);

            // set up pointers to (D|g(r,r')|G)^(m) integrals

            auto g_xx_xxxx = primBuffer.data(goff + 90 * idx);

            auto g_xx_xxxy = primBuffer.data(goff + 90 * idx + 1);

            auto g_xx_xxxz = primBuffer.data(goff + 90 * idx + 2);

            auto g_xx_xxyy = primBuffer.data(goff + 90 * idx + 3);

            auto g_xx_xxyz = primBuffer.data(goff + 90 * idx + 4);

            auto g_xx_xxzz = primBuffer.data(goff + 90 * idx + 5);

            auto g_xx_xyyy = primBuffer.data(goff + 90 * idx + 6);

            auto g_xx_xyyz = primBuffer.data(goff + 90 * idx + 7);

            auto g_xx_xyzz = primBuffer.data(goff + 90 * idx + 8);

            auto g_xx_xzzz = primBuffer.data(goff + 90 * idx + 9);

            auto g_xx_yyyy = primBuffer.data(goff + 90 * idx + 10);

            auto g_xx_yyyz = primBuffer.data(goff + 90 * idx + 11);

            auto g_xx_yyzz = primBuffer.data(goff + 90 * idx + 12);

            auto g_xx_yzzz = primBuffer.data(goff + 90 * idx + 13);

            auto g_xx_zzzz = primBuffer.data(goff + 90 * idx + 14);

            auto g_xy_xxxx = primBuffer.data(goff + 90 * idx + 15);

            auto g_xy_xxxy = primBuffer.data(goff + 90 * idx + 16);

            auto g_xy_xxxz = primBuffer.data(goff + 90 * idx + 17);

            auto g_xy_xxyy = primBuffer.data(goff + 90 * idx + 18);

            auto g_xy_xxyz = primBuffer.data(goff + 90 * idx + 19);

            auto g_xy_xxzz = primBuffer.data(goff + 90 * idx + 20);

            auto g_xy_xyyy = primBuffer.data(goff + 90 * idx + 21);

            auto g_xy_xyyz = primBuffer.data(goff + 90 * idx + 22);

            auto g_xy_xyzz = primBuffer.data(goff + 90 * idx + 23);

            auto g_xy_xzzz = primBuffer.data(goff + 90 * idx + 24);

            auto g_xy_yyyy = primBuffer.data(goff + 90 * idx + 25);

            auto g_xy_yyyz = primBuffer.data(goff + 90 * idx + 26);

            auto g_xy_yyzz = primBuffer.data(goff + 90 * idx + 27);

            auto g_xy_yzzz = primBuffer.data(goff + 90 * idx + 28);

            auto g_xy_zzzz = primBuffer.data(goff + 90 * idx + 29);

            auto g_xz_xxxx = primBuffer.data(goff + 90 * idx + 30);

            auto g_xz_xxxy = primBuffer.data(goff + 90 * idx + 31);

            auto g_xz_xxxz = primBuffer.data(goff + 90 * idx + 32);

            auto g_xz_xxyy = primBuffer.data(goff + 90 * idx + 33);

            auto g_xz_xxyz = primBuffer.data(goff + 90 * idx + 34);

            auto g_xz_xxzz = primBuffer.data(goff + 90 * idx + 35);

            auto g_xz_xyyy = primBuffer.data(goff + 90 * idx + 36);

            auto g_xz_xyyz = primBuffer.data(goff + 90 * idx + 37);

            auto g_xz_xyzz = primBuffer.data(goff + 90 * idx + 38);

            auto g_xz_xzzz = primBuffer.data(goff + 90 * idx + 39);

            auto g_xz_yyyy = primBuffer.data(goff + 90 * idx + 40);

            auto g_xz_yyyz = primBuffer.data(goff + 90 * idx + 41);

            auto g_xz_yyzz = primBuffer.data(goff + 90 * idx + 42);

            auto g_xz_yzzz = primBuffer.data(goff + 90 * idx + 43);

            auto g_xz_zzzz = primBuffer.data(goff + 90 * idx + 44);

            auto g_yy_xxxx = primBuffer.data(goff + 90 * idx + 45);

            auto g_yy_xxxy = primBuffer.data(goff + 90 * idx + 46);

            auto g_yy_xxxz = primBuffer.data(goff + 90 * idx + 47);

            auto g_yy_xxyy = primBuffer.data(goff + 90 * idx + 48);

            auto g_yy_xxyz = primBuffer.data(goff + 90 * idx + 49);

            auto g_yy_xxzz = primBuffer.data(goff + 90 * idx + 50);

            auto g_yy_xyyy = primBuffer.data(goff + 90 * idx + 51);

            auto g_yy_xyyz = primBuffer.data(goff + 90 * idx + 52);

            auto g_yy_xyzz = primBuffer.data(goff + 90 * idx + 53);

            auto g_yy_xzzz = primBuffer.data(goff + 90 * idx + 54);

            auto g_yy_yyyy = primBuffer.data(goff + 90 * idx + 55);

            auto g_yy_yyyz = primBuffer.data(goff + 90 * idx + 56);

            auto g_yy_yyzz = primBuffer.data(goff + 90 * idx + 57);

            auto g_yy_yzzz = primBuffer.data(goff + 90 * idx + 58);

            auto g_yy_zzzz = primBuffer.data(goff + 90 * idx + 59);

            auto g_yz_xxxx = primBuffer.data(goff + 90 * idx + 60);

            auto g_yz_xxxy = primBuffer.data(goff + 90 * idx + 61);

            auto g_yz_xxxz = primBuffer.data(goff + 90 * idx + 62);

            auto g_yz_xxyy = primBuffer.data(goff + 90 * idx + 63);

            auto g_yz_xxyz = primBuffer.data(goff + 90 * idx + 64);

            auto g_yz_xxzz = primBuffer.data(goff + 90 * idx + 65);

            auto g_yz_xyyy = primBuffer.data(goff + 90 * idx + 66);

            auto g_yz_xyyz = primBuffer.data(goff + 90 * idx + 67);

            auto g_yz_xyzz = primBuffer.data(goff + 90 * idx + 68);

            auto g_yz_xzzz = primBuffer.data(goff + 90 * idx + 69);

            auto g_yz_yyyy = primBuffer.data(goff + 90 * idx + 70);

            auto g_yz_yyyz = primBuffer.data(goff + 90 * idx + 71);

            auto g_yz_yyzz = primBuffer.data(goff + 90 * idx + 72);

            auto g_yz_yzzz = primBuffer.data(goff + 90 * idx + 73);

            auto g_yz_zzzz = primBuffer.data(goff + 90 * idx + 74);

            auto g_zz_xxxx = primBuffer.data(goff + 90 * idx + 75);

            auto g_zz_xxxy = primBuffer.data(goff + 90 * idx + 76);

            auto g_zz_xxxz = primBuffer.data(goff + 90 * idx + 77);

            auto g_zz_xxyy = primBuffer.data(goff + 90 * idx + 78);

            auto g_zz_xxyz = primBuffer.data(goff + 90 * idx + 79);

            auto g_zz_xxzz = primBuffer.data(goff + 90 * idx + 80);

            auto g_zz_xyyy = primBuffer.data(goff + 90 * idx + 81);

            auto g_zz_xyyz = primBuffer.data(goff + 90 * idx + 82);

            auto g_zz_xyzz = primBuffer.data(goff + 90 * idx + 83);

            auto g_zz_xzzz = primBuffer.data(goff + 90 * idx + 84);

            auto g_zz_yyyy = primBuffer.data(goff + 90 * idx + 85);

            auto g_zz_yyyz = primBuffer.data(goff + 90 * idx + 86);

            auto g_zz_yyzz = primBuffer.data(goff + 90 * idx + 87);

            auto g_zz_yzzz = primBuffer.data(goff + 90 * idx + 88);

            auto g_zz_zzzz = primBuffer.data(goff + 90 * idx + 89);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_x_xxx, gk_x_xxy,\
                                         gk_x_xxz, gk_x_xyy, gk_x_xyz, gk_x_xzz,\
                                         gk_x_yyy, gk_x_yyz, gk_x_yzz, gk_x_zzz,\
                                         gk_y_xxx, gk_y_xxy, gk_y_xxz, gk_y_xyy,\
                                         gk_y_xyz, gk_y_xzz, gk_y_yyy, gk_y_yyz,\
                                         gk_y_yzz, gk_y_zzz, gk_z_xxx, gk_z_xxy,\
                                         gk_z_xxz, gk_z_xyy, gk_z_xyz, gk_z_xzz,\
                                         gk_z_yyy, gk_z_yyz, gk_z_yzz, gk_z_zzz,\
                                         g20_0_xxxx, g20_0_xxxy, g20_0_xxxz, g20_0_xxyy,\
                                         g20_0_xxyz, g20_0_xxzz, g20_0_xyyy, g20_0_xyyz,\
                                         g20_0_xyzz, g20_0_xzzz, g20_0_yyyy, g20_0_yyyz,\
                                         g20_0_yyzz, g20_0_yzzz, g20_0_zzzz, g21_0_xxxx,\
                                         g21_0_xxxy, g21_0_xxxz, g21_0_xxyy, g21_0_xxyz,\
                                         g21_0_xxzz, g21_0_xyyy, g21_0_xyyz, g21_0_xyzz,\
                                         g21_0_xzzz, g21_0_yyyy, g21_0_yyyz, g21_0_yyzz,\
                                         g21_0_yzzz, g21_0_zzzz, g1_x_xxxx, g1_x_xxxy,\
                                         g1_x_xxxz, g1_x_xxyy, g1_x_xxyz, g1_x_xxzz,\
                                         g1_x_xyyy, g1_x_xyyz, g1_x_xyzz, g1_x_xzzz,\
                                         g1_x_yyyy, g1_x_yyyz, g1_x_yyzz, g1_x_yzzz,\
                                         g1_x_zzzz, g1_y_xxxx, g1_y_xxxy, g1_y_xxxz,\
                                         g1_y_xxyy, g1_y_xxyz, g1_y_xxzz, g1_y_xyyy,\
                                         g1_y_xyyz, g1_y_xyzz, g1_y_xzzz, g1_y_yyyy,\
                                         g1_y_yyyz, g1_y_yyzz, g1_y_yzzz, g1_y_zzzz,\
                                         g1_z_xxxx, g1_z_xxxy, g1_z_xxxz, g1_z_xxyy,\
                                         g1_z_xxyz, g1_z_xxzz, g1_z_xyyy, g1_z_xyyz,\
                                         g1_z_xyzz, g1_z_xzzz, g1_z_yyyy, g1_z_yyyz,\
                                         g1_z_yyzz, g1_z_yzzz, g1_z_zzzz, g_xx_xxxx,\
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
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xx_xxxx[k] = fr * g1_x_xxxx[k] + f2g * (g20_0_xxxx[k] - fgz * g21_0_xxxx[k]) + 4.0 * f2t * gk_x_xxx[k];

                g_xx_xxxy[k] = fr * g1_x_xxxy[k] + f2g * (g20_0_xxxy[k] - fgz * g21_0_xxxy[k]) + 3.0 * f2t * gk_x_xxy[k];

                g_xx_xxxz[k] = fr * g1_x_xxxz[k] + f2g * (g20_0_xxxz[k] - fgz * g21_0_xxxz[k]) + 3.0 * f2t * gk_x_xxz[k];

                g_xx_xxyy[k] = fr * g1_x_xxyy[k] + f2g * (g20_0_xxyy[k] - fgz * g21_0_xxyy[k]) + 2.0 * f2t * gk_x_xyy[k];

                g_xx_xxyz[k] = fr * g1_x_xxyz[k] + f2g * (g20_0_xxyz[k] - fgz * g21_0_xxyz[k]) + 2.0 * f2t * gk_x_xyz[k];

                g_xx_xxzz[k] = fr * g1_x_xxzz[k] + f2g * (g20_0_xxzz[k] - fgz * g21_0_xxzz[k]) + 2.0 * f2t * gk_x_xzz[k];

                g_xx_xyyy[k] = fr * g1_x_xyyy[k] + f2g * (g20_0_xyyy[k] - fgz * g21_0_xyyy[k]) + f2t * gk_x_yyy[k];

                g_xx_xyyz[k] = fr * g1_x_xyyz[k] + f2g * (g20_0_xyyz[k] - fgz * g21_0_xyyz[k]) + f2t * gk_x_yyz[k];

                g_xx_xyzz[k] = fr * g1_x_xyzz[k] + f2g * (g20_0_xyzz[k] - fgz * g21_0_xyzz[k]) + f2t * gk_x_yzz[k];

                g_xx_xzzz[k] = fr * g1_x_xzzz[k] + f2g * (g20_0_xzzz[k] - fgz * g21_0_xzzz[k]) + f2t * gk_x_zzz[k];

                g_xx_yyyy[k] = fr * g1_x_yyyy[k] + f2g * (g20_0_yyyy[k] - fgz * g21_0_yyyy[k]);

                g_xx_yyyz[k] = fr * g1_x_yyyz[k] + f2g * (g20_0_yyyz[k] - fgz * g21_0_yyyz[k]);

                g_xx_yyzz[k] = fr * g1_x_yyzz[k] + f2g * (g20_0_yyzz[k] - fgz * g21_0_yyzz[k]);

                g_xx_yzzz[k] = fr * g1_x_yzzz[k] + f2g * (g20_0_yzzz[k] - fgz * g21_0_yzzz[k]);

                g_xx_zzzz[k] = fr * g1_x_zzzz[k] + f2g * (g20_0_zzzz[k] - fgz * g21_0_zzzz[k]);

                g_xy_xxxx[k] = fr * g1_y_xxxx[k] + 4.0 * f2t * gk_y_xxx[k];

                g_xy_xxxy[k] = fr * g1_y_xxxy[k] + 3.0 * f2t * gk_y_xxy[k];

                g_xy_xxxz[k] = fr * g1_y_xxxz[k] + 3.0 * f2t * gk_y_xxz[k];

                g_xy_xxyy[k] = fr * g1_y_xxyy[k] + 2.0 * f2t * gk_y_xyy[k];

                g_xy_xxyz[k] = fr * g1_y_xxyz[k] + 2.0 * f2t * gk_y_xyz[k];

                g_xy_xxzz[k] = fr * g1_y_xxzz[k] + 2.0 * f2t * gk_y_xzz[k];

                g_xy_xyyy[k] = fr * g1_y_xyyy[k] + f2t * gk_y_yyy[k];

                g_xy_xyyz[k] = fr * g1_y_xyyz[k] + f2t * gk_y_yyz[k];

                g_xy_xyzz[k] = fr * g1_y_xyzz[k] + f2t * gk_y_yzz[k];

                g_xy_xzzz[k] = fr * g1_y_xzzz[k] + f2t * gk_y_zzz[k];

                g_xy_yyyy[k] = fr * g1_y_yyyy[k];

                g_xy_yyyz[k] = fr * g1_y_yyyz[k];

                g_xy_yyzz[k] = fr * g1_y_yyzz[k];

                g_xy_yzzz[k] = fr * g1_y_yzzz[k];

                g_xy_zzzz[k] = fr * g1_y_zzzz[k];

                g_xz_xxxx[k] = fr * g1_z_xxxx[k] + 4.0 * f2t * gk_z_xxx[k];

                g_xz_xxxy[k] = fr * g1_z_xxxy[k] + 3.0 * f2t * gk_z_xxy[k];

                g_xz_xxxz[k] = fr * g1_z_xxxz[k] + 3.0 * f2t * gk_z_xxz[k];

                g_xz_xxyy[k] = fr * g1_z_xxyy[k] + 2.0 * f2t * gk_z_xyy[k];

                g_xz_xxyz[k] = fr * g1_z_xxyz[k] + 2.0 * f2t * gk_z_xyz[k];

                g_xz_xxzz[k] = fr * g1_z_xxzz[k] + 2.0 * f2t * gk_z_xzz[k];

                g_xz_xyyy[k] = fr * g1_z_xyyy[k] + f2t * gk_z_yyy[k];

                g_xz_xyyz[k] = fr * g1_z_xyyz[k] + f2t * gk_z_yyz[k];

                g_xz_xyzz[k] = fr * g1_z_xyzz[k] + f2t * gk_z_yzz[k];

                g_xz_xzzz[k] = fr * g1_z_xzzz[k] + f2t * gk_z_zzz[k];

                g_xz_yyyy[k] = fr * g1_z_yyyy[k];

                g_xz_yyyz[k] = fr * g1_z_yyyz[k];

                g_xz_yyzz[k] = fr * g1_z_yyzz[k];

                g_xz_yzzz[k] = fr * g1_z_yzzz[k];

                g_xz_zzzz[k] = fr * g1_z_zzzz[k];

                // leading y component

                fr = pay[k];

                g_yy_xxxx[k] = fr * g1_y_xxxx[k] + f2g * (g20_0_xxxx[k] - fgz * g21_0_xxxx[k]);

                g_yy_xxxy[k] = fr * g1_y_xxxy[k] + f2g * (g20_0_xxxy[k] - fgz * g21_0_xxxy[k]) + f2t * gk_y_xxx[k];

                g_yy_xxxz[k] = fr * g1_y_xxxz[k] + f2g * (g20_0_xxxz[k] - fgz * g21_0_xxxz[k]);

                g_yy_xxyy[k] = fr * g1_y_xxyy[k] + f2g * (g20_0_xxyy[k] - fgz * g21_0_xxyy[k]) + 2.0 * f2t * gk_y_xxy[k];

                g_yy_xxyz[k] = fr * g1_y_xxyz[k] + f2g * (g20_0_xxyz[k] - fgz * g21_0_xxyz[k]) + f2t * gk_y_xxz[k];

                g_yy_xxzz[k] = fr * g1_y_xxzz[k] + f2g * (g20_0_xxzz[k] - fgz * g21_0_xxzz[k]);

                g_yy_xyyy[k] = fr * g1_y_xyyy[k] + f2g * (g20_0_xyyy[k] - fgz * g21_0_xyyy[k]) + 3.0 * f2t * gk_y_xyy[k];

                g_yy_xyyz[k] = fr * g1_y_xyyz[k] + f2g * (g20_0_xyyz[k] - fgz * g21_0_xyyz[k]) + 2.0 * f2t * gk_y_xyz[k];

                g_yy_xyzz[k] = fr * g1_y_xyzz[k] + f2g * (g20_0_xyzz[k] - fgz * g21_0_xyzz[k]) + f2t * gk_y_xzz[k];

                g_yy_xzzz[k] = fr * g1_y_xzzz[k] + f2g * (g20_0_xzzz[k] - fgz * g21_0_xzzz[k]);

                g_yy_yyyy[k] = fr * g1_y_yyyy[k] + f2g * (g20_0_yyyy[k] - fgz * g21_0_yyyy[k]) + 4.0 * f2t * gk_y_yyy[k];

                g_yy_yyyz[k] = fr * g1_y_yyyz[k] + f2g * (g20_0_yyyz[k] - fgz * g21_0_yyyz[k]) + 3.0 * f2t * gk_y_yyz[k];

                g_yy_yyzz[k] = fr * g1_y_yyzz[k] + f2g * (g20_0_yyzz[k] - fgz * g21_0_yyzz[k]) + 2.0 * f2t * gk_y_yzz[k];

                g_yy_yzzz[k] = fr * g1_y_yzzz[k] + f2g * (g20_0_yzzz[k] - fgz * g21_0_yzzz[k]) + f2t * gk_y_zzz[k];

                g_yy_zzzz[k] = fr * g1_y_zzzz[k] + f2g * (g20_0_zzzz[k] - fgz * g21_0_zzzz[k]);

                g_yz_xxxx[k] = fr * g1_z_xxxx[k];

                g_yz_xxxy[k] = fr * g1_z_xxxy[k] + f2t * gk_z_xxx[k];

                g_yz_xxxz[k] = fr * g1_z_xxxz[k];

                g_yz_xxyy[k] = fr * g1_z_xxyy[k] + 2.0 * f2t * gk_z_xxy[k];

                g_yz_xxyz[k] = fr * g1_z_xxyz[k] + f2t * gk_z_xxz[k];

                g_yz_xxzz[k] = fr * g1_z_xxzz[k];

                g_yz_xyyy[k] = fr * g1_z_xyyy[k] + 3.0 * f2t * gk_z_xyy[k];

                g_yz_xyyz[k] = fr * g1_z_xyyz[k] + 2.0 * f2t * gk_z_xyz[k];

                g_yz_xyzz[k] = fr * g1_z_xyzz[k] + f2t * gk_z_xzz[k];

                g_yz_xzzz[k] = fr * g1_z_xzzz[k];

                g_yz_yyyy[k] = fr * g1_z_yyyy[k] + 4.0 * f2t * gk_z_yyy[k];

                g_yz_yyyz[k] = fr * g1_z_yyyz[k] + 3.0 * f2t * gk_z_yyz[k];

                g_yz_yyzz[k] = fr * g1_z_yyzz[k] + 2.0 * f2t * gk_z_yzz[k];

                g_yz_yzzz[k] = fr * g1_z_yzzz[k] + f2t * gk_z_zzz[k];

                g_yz_zzzz[k] = fr * g1_z_zzzz[k];

                // leading z component

                fr = paz[k];

                g_zz_xxxx[k] = fr * g1_z_xxxx[k] + f2g * (g20_0_xxxx[k] - fgz * g21_0_xxxx[k]);

                g_zz_xxxy[k] = fr * g1_z_xxxy[k] + f2g * (g20_0_xxxy[k] - fgz * g21_0_xxxy[k]);

                g_zz_xxxz[k] = fr * g1_z_xxxz[k] + f2g * (g20_0_xxxz[k] - fgz * g21_0_xxxz[k]) + f2t * gk_z_xxx[k];

                g_zz_xxyy[k] = fr * g1_z_xxyy[k] + f2g * (g20_0_xxyy[k] - fgz * g21_0_xxyy[k]);

                g_zz_xxyz[k] = fr * g1_z_xxyz[k] + f2g * (g20_0_xxyz[k] - fgz * g21_0_xxyz[k]) + f2t * gk_z_xxy[k];

                g_zz_xxzz[k] = fr * g1_z_xxzz[k] + f2g * (g20_0_xxzz[k] - fgz * g21_0_xxzz[k]) + 2.0 * f2t * gk_z_xxz[k];

                g_zz_xyyy[k] = fr * g1_z_xyyy[k] + f2g * (g20_0_xyyy[k] - fgz * g21_0_xyyy[k]);

                g_zz_xyyz[k] = fr * g1_z_xyyz[k] + f2g * (g20_0_xyyz[k] - fgz * g21_0_xyyz[k]) + f2t * gk_z_xyy[k];

                g_zz_xyzz[k] = fr * g1_z_xyzz[k] + f2g * (g20_0_xyzz[k] - fgz * g21_0_xyzz[k]) + 2.0 * f2t * gk_z_xyz[k];

                g_zz_xzzz[k] = fr * g1_z_xzzz[k] + f2g * (g20_0_xzzz[k] - fgz * g21_0_xzzz[k]) + 3.0 * f2t * gk_z_xzz[k];

                g_zz_yyyy[k] = fr * g1_z_yyyy[k] + f2g * (g20_0_yyyy[k] - fgz * g21_0_yyyy[k]);

                g_zz_yyyz[k] = fr * g1_z_yyyz[k] + f2g * (g20_0_yyyz[k] - fgz * g21_0_yyyz[k]) + f2t * gk_z_yyy[k];

                g_zz_yyzz[k] = fr * g1_z_yyzz[k] + f2g * (g20_0_yyzz[k] - fgz * g21_0_yyzz[k]) + 2.0 * f2t * gk_z_yyz[k];

                g_zz_yzzz[k] = fr * g1_z_yzzz[k] + f2g * (g20_0_yzzz[k] - fgz * g21_0_yzzz[k]) + 3.0 * f2t * gk_z_yzz[k];

                g_zz_zzzz[k] = fr * g1_z_zzzz[k] + f2g * (g20_0_zzzz[k] - fgz * g21_0_zzzz[k]) + 4.0 * f2t * gk_z_zzz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForGD(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 4, 2);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {4, 2, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 2, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 1, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|g(r,r')|P)^(m+1) integrals

            auto gk_xxx_x = primBuffer.data(gkoff + 30 * idx);

            auto gk_xxx_y = primBuffer.data(gkoff + 30 * idx + 1);

            auto gk_xxx_z = primBuffer.data(gkoff + 30 * idx + 2);

            auto gk_xxy_x = primBuffer.data(gkoff + 30 * idx + 3);

            auto gk_xxy_y = primBuffer.data(gkoff + 30 * idx + 4);

            auto gk_xxy_z = primBuffer.data(gkoff + 30 * idx + 5);

            auto gk_xxz_x = primBuffer.data(gkoff + 30 * idx + 6);

            auto gk_xxz_y = primBuffer.data(gkoff + 30 * idx + 7);

            auto gk_xxz_z = primBuffer.data(gkoff + 30 * idx + 8);

            auto gk_xyy_x = primBuffer.data(gkoff + 30 * idx + 9);

            auto gk_xyy_y = primBuffer.data(gkoff + 30 * idx + 10);

            auto gk_xyy_z = primBuffer.data(gkoff + 30 * idx + 11);

            auto gk_xyz_x = primBuffer.data(gkoff + 30 * idx + 12);

            auto gk_xyz_y = primBuffer.data(gkoff + 30 * idx + 13);

            auto gk_xyz_z = primBuffer.data(gkoff + 30 * idx + 14);

            auto gk_xzz_x = primBuffer.data(gkoff + 30 * idx + 15);

            auto gk_xzz_y = primBuffer.data(gkoff + 30 * idx + 16);

            auto gk_xzz_z = primBuffer.data(gkoff + 30 * idx + 17);

            auto gk_yyy_x = primBuffer.data(gkoff + 30 * idx + 18);

            auto gk_yyy_y = primBuffer.data(gkoff + 30 * idx + 19);

            auto gk_yyy_z = primBuffer.data(gkoff + 30 * idx + 20);

            auto gk_yyz_x = primBuffer.data(gkoff + 30 * idx + 21);

            auto gk_yyz_y = primBuffer.data(gkoff + 30 * idx + 22);

            auto gk_yyz_z = primBuffer.data(gkoff + 30 * idx + 23);

            auto gk_yzz_x = primBuffer.data(gkoff + 30 * idx + 24);

            auto gk_yzz_y = primBuffer.data(gkoff + 30 * idx + 25);

            auto gk_yzz_z = primBuffer.data(gkoff + 30 * idx + 26);

            auto gk_zzz_x = primBuffer.data(gkoff + 30 * idx + 27);

            auto gk_zzz_y = primBuffer.data(gkoff + 30 * idx + 28);

            auto gk_zzz_z = primBuffer.data(gkoff + 30 * idx + 29);

            // set up pointers to (D|g(r,r')|D)^(m) integrals

            auto g20_xx_xx = primBuffer.data(g20off + 36 * idx);

            auto g20_xx_xy = primBuffer.data(g20off + 36 * idx + 1);

            auto g20_xx_xz = primBuffer.data(g20off + 36 * idx + 2);

            auto g20_xx_yy = primBuffer.data(g20off + 36 * idx + 3);

            auto g20_xx_yz = primBuffer.data(g20off + 36 * idx + 4);

            auto g20_xx_zz = primBuffer.data(g20off + 36 * idx + 5);

            auto g20_xy_xx = primBuffer.data(g20off + 36 * idx + 6);

            auto g20_xy_xy = primBuffer.data(g20off + 36 * idx + 7);

            auto g20_xy_xz = primBuffer.data(g20off + 36 * idx + 8);

            auto g20_xy_yy = primBuffer.data(g20off + 36 * idx + 9);

            auto g20_xy_yz = primBuffer.data(g20off + 36 * idx + 10);

            auto g20_xy_zz = primBuffer.data(g20off + 36 * idx + 11);

            auto g20_xz_xx = primBuffer.data(g20off + 36 * idx + 12);

            auto g20_xz_xy = primBuffer.data(g20off + 36 * idx + 13);

            auto g20_xz_xz = primBuffer.data(g20off + 36 * idx + 14);

            auto g20_xz_yy = primBuffer.data(g20off + 36 * idx + 15);

            auto g20_xz_yz = primBuffer.data(g20off + 36 * idx + 16);

            auto g20_xz_zz = primBuffer.data(g20off + 36 * idx + 17);

            auto g20_yy_xx = primBuffer.data(g20off + 36 * idx + 18);

            auto g20_yy_xy = primBuffer.data(g20off + 36 * idx + 19);

            auto g20_yy_xz = primBuffer.data(g20off + 36 * idx + 20);

            auto g20_yy_yy = primBuffer.data(g20off + 36 * idx + 21);

            auto g20_yy_yz = primBuffer.data(g20off + 36 * idx + 22);

            auto g20_yy_zz = primBuffer.data(g20off + 36 * idx + 23);

            auto g20_yz_xx = primBuffer.data(g20off + 36 * idx + 24);

            auto g20_yz_xy = primBuffer.data(g20off + 36 * idx + 25);

            auto g20_yz_xz = primBuffer.data(g20off + 36 * idx + 26);

            auto g20_yz_yy = primBuffer.data(g20off + 36 * idx + 27);

            auto g20_yz_yz = primBuffer.data(g20off + 36 * idx + 28);

            auto g20_yz_zz = primBuffer.data(g20off + 36 * idx + 29);

            auto g20_zz_xx = primBuffer.data(g20off + 36 * idx + 30);

            auto g20_zz_xy = primBuffer.data(g20off + 36 * idx + 31);

            auto g20_zz_xz = primBuffer.data(g20off + 36 * idx + 32);

            auto g20_zz_yy = primBuffer.data(g20off + 36 * idx + 33);

            auto g20_zz_yz = primBuffer.data(g20off + 36 * idx + 34);

            auto g20_zz_zz = primBuffer.data(g20off + 36 * idx + 35);

            // set up pointers to (D|g(r,r')|D)^(m+1) integrals

            auto g21_xx_xx = primBuffer.data(g21off + 36 * idx);

            auto g21_xx_xy = primBuffer.data(g21off + 36 * idx + 1);

            auto g21_xx_xz = primBuffer.data(g21off + 36 * idx + 2);

            auto g21_xx_yy = primBuffer.data(g21off + 36 * idx + 3);

            auto g21_xx_yz = primBuffer.data(g21off + 36 * idx + 4);

            auto g21_xx_zz = primBuffer.data(g21off + 36 * idx + 5);

            auto g21_xy_xx = primBuffer.data(g21off + 36 * idx + 6);

            auto g21_xy_xy = primBuffer.data(g21off + 36 * idx + 7);

            auto g21_xy_xz = primBuffer.data(g21off + 36 * idx + 8);

            auto g21_xy_yy = primBuffer.data(g21off + 36 * idx + 9);

            auto g21_xy_yz = primBuffer.data(g21off + 36 * idx + 10);

            auto g21_xy_zz = primBuffer.data(g21off + 36 * idx + 11);

            auto g21_xz_xx = primBuffer.data(g21off + 36 * idx + 12);

            auto g21_xz_xy = primBuffer.data(g21off + 36 * idx + 13);

            auto g21_xz_xz = primBuffer.data(g21off + 36 * idx + 14);

            auto g21_xz_yy = primBuffer.data(g21off + 36 * idx + 15);

            auto g21_xz_yz = primBuffer.data(g21off + 36 * idx + 16);

            auto g21_xz_zz = primBuffer.data(g21off + 36 * idx + 17);

            auto g21_yy_xx = primBuffer.data(g21off + 36 * idx + 18);

            auto g21_yy_xy = primBuffer.data(g21off + 36 * idx + 19);

            auto g21_yy_xz = primBuffer.data(g21off + 36 * idx + 20);

            auto g21_yy_yy = primBuffer.data(g21off + 36 * idx + 21);

            auto g21_yy_yz = primBuffer.data(g21off + 36 * idx + 22);

            auto g21_yy_zz = primBuffer.data(g21off + 36 * idx + 23);

            auto g21_yz_xx = primBuffer.data(g21off + 36 * idx + 24);

            auto g21_yz_xy = primBuffer.data(g21off + 36 * idx + 25);

            auto g21_yz_xz = primBuffer.data(g21off + 36 * idx + 26);

            auto g21_yz_yy = primBuffer.data(g21off + 36 * idx + 27);

            auto g21_yz_yz = primBuffer.data(g21off + 36 * idx + 28);

            auto g21_yz_zz = primBuffer.data(g21off + 36 * idx + 29);

            auto g21_zz_xx = primBuffer.data(g21off + 36 * idx + 30);

            auto g21_zz_xy = primBuffer.data(g21off + 36 * idx + 31);

            auto g21_zz_xz = primBuffer.data(g21off + 36 * idx + 32);

            auto g21_zz_yy = primBuffer.data(g21off + 36 * idx + 33);

            auto g21_zz_yz = primBuffer.data(g21off + 36 * idx + 34);

            auto g21_zz_zz = primBuffer.data(g21off + 36 * idx + 35);

            // set up pointers to (F|g(r,r')|D)^(m+1) integrals

            auto g1_xxx_xx = primBuffer.data(g1off + 60 * idx);

            auto g1_xxx_xy = primBuffer.data(g1off + 60 * idx + 1);

            auto g1_xxx_xz = primBuffer.data(g1off + 60 * idx + 2);

            auto g1_xxx_yy = primBuffer.data(g1off + 60 * idx + 3);

            auto g1_xxx_yz = primBuffer.data(g1off + 60 * idx + 4);

            auto g1_xxx_zz = primBuffer.data(g1off + 60 * idx + 5);

            auto g1_xxy_xx = primBuffer.data(g1off + 60 * idx + 6);

            auto g1_xxy_xy = primBuffer.data(g1off + 60 * idx + 7);

            auto g1_xxy_xz = primBuffer.data(g1off + 60 * idx + 8);

            auto g1_xxy_yy = primBuffer.data(g1off + 60 * idx + 9);

            auto g1_xxy_yz = primBuffer.data(g1off + 60 * idx + 10);

            auto g1_xxy_zz = primBuffer.data(g1off + 60 * idx + 11);

            auto g1_xxz_xx = primBuffer.data(g1off + 60 * idx + 12);

            auto g1_xxz_xy = primBuffer.data(g1off + 60 * idx + 13);

            auto g1_xxz_xz = primBuffer.data(g1off + 60 * idx + 14);

            auto g1_xxz_yy = primBuffer.data(g1off + 60 * idx + 15);

            auto g1_xxz_yz = primBuffer.data(g1off + 60 * idx + 16);

            auto g1_xxz_zz = primBuffer.data(g1off + 60 * idx + 17);

            auto g1_xyy_xx = primBuffer.data(g1off + 60 * idx + 18);

            auto g1_xyy_xy = primBuffer.data(g1off + 60 * idx + 19);

            auto g1_xyy_xz = primBuffer.data(g1off + 60 * idx + 20);

            auto g1_xyy_yy = primBuffer.data(g1off + 60 * idx + 21);

            auto g1_xyy_yz = primBuffer.data(g1off + 60 * idx + 22);

            auto g1_xyy_zz = primBuffer.data(g1off + 60 * idx + 23);

            auto g1_xyz_xx = primBuffer.data(g1off + 60 * idx + 24);

            auto g1_xyz_xy = primBuffer.data(g1off + 60 * idx + 25);

            auto g1_xyz_xz = primBuffer.data(g1off + 60 * idx + 26);

            auto g1_xyz_yy = primBuffer.data(g1off + 60 * idx + 27);

            auto g1_xyz_yz = primBuffer.data(g1off + 60 * idx + 28);

            auto g1_xyz_zz = primBuffer.data(g1off + 60 * idx + 29);

            auto g1_xzz_xx = primBuffer.data(g1off + 60 * idx + 30);

            auto g1_xzz_xy = primBuffer.data(g1off + 60 * idx + 31);

            auto g1_xzz_xz = primBuffer.data(g1off + 60 * idx + 32);

            auto g1_xzz_yy = primBuffer.data(g1off + 60 * idx + 33);

            auto g1_xzz_yz = primBuffer.data(g1off + 60 * idx + 34);

            auto g1_xzz_zz = primBuffer.data(g1off + 60 * idx + 35);

            auto g1_yyy_xx = primBuffer.data(g1off + 60 * idx + 36);

            auto g1_yyy_xy = primBuffer.data(g1off + 60 * idx + 37);

            auto g1_yyy_xz = primBuffer.data(g1off + 60 * idx + 38);

            auto g1_yyy_yy = primBuffer.data(g1off + 60 * idx + 39);

            auto g1_yyy_yz = primBuffer.data(g1off + 60 * idx + 40);

            auto g1_yyy_zz = primBuffer.data(g1off + 60 * idx + 41);

            auto g1_yyz_xx = primBuffer.data(g1off + 60 * idx + 42);

            auto g1_yyz_xy = primBuffer.data(g1off + 60 * idx + 43);

            auto g1_yyz_xz = primBuffer.data(g1off + 60 * idx + 44);

            auto g1_yyz_yy = primBuffer.data(g1off + 60 * idx + 45);

            auto g1_yyz_yz = primBuffer.data(g1off + 60 * idx + 46);

            auto g1_yyz_zz = primBuffer.data(g1off + 60 * idx + 47);

            auto g1_yzz_xx = primBuffer.data(g1off + 60 * idx + 48);

            auto g1_yzz_xy = primBuffer.data(g1off + 60 * idx + 49);

            auto g1_yzz_xz = primBuffer.data(g1off + 60 * idx + 50);

            auto g1_yzz_yy = primBuffer.data(g1off + 60 * idx + 51);

            auto g1_yzz_yz = primBuffer.data(g1off + 60 * idx + 52);

            auto g1_yzz_zz = primBuffer.data(g1off + 60 * idx + 53);

            auto g1_zzz_xx = primBuffer.data(g1off + 60 * idx + 54);

            auto g1_zzz_xy = primBuffer.data(g1off + 60 * idx + 55);

            auto g1_zzz_xz = primBuffer.data(g1off + 60 * idx + 56);

            auto g1_zzz_yy = primBuffer.data(g1off + 60 * idx + 57);

            auto g1_zzz_yz = primBuffer.data(g1off + 60 * idx + 58);

            auto g1_zzz_zz = primBuffer.data(g1off + 60 * idx + 59);

            // set up pointers to (G|g(r,r')|D)^(m) integrals

            auto g_xxxx_xx = primBuffer.data(goff + 90 * idx);

            auto g_xxxx_xy = primBuffer.data(goff + 90 * idx + 1);

            auto g_xxxx_xz = primBuffer.data(goff + 90 * idx + 2);

            auto g_xxxx_yy = primBuffer.data(goff + 90 * idx + 3);

            auto g_xxxx_yz = primBuffer.data(goff + 90 * idx + 4);

            auto g_xxxx_zz = primBuffer.data(goff + 90 * idx + 5);

            auto g_xxxy_xx = primBuffer.data(goff + 90 * idx + 6);

            auto g_xxxy_xy = primBuffer.data(goff + 90 * idx + 7);

            auto g_xxxy_xz = primBuffer.data(goff + 90 * idx + 8);

            auto g_xxxy_yy = primBuffer.data(goff + 90 * idx + 9);

            auto g_xxxy_yz = primBuffer.data(goff + 90 * idx + 10);

            auto g_xxxy_zz = primBuffer.data(goff + 90 * idx + 11);

            auto g_xxxz_xx = primBuffer.data(goff + 90 * idx + 12);

            auto g_xxxz_xy = primBuffer.data(goff + 90 * idx + 13);

            auto g_xxxz_xz = primBuffer.data(goff + 90 * idx + 14);

            auto g_xxxz_yy = primBuffer.data(goff + 90 * idx + 15);

            auto g_xxxz_yz = primBuffer.data(goff + 90 * idx + 16);

            auto g_xxxz_zz = primBuffer.data(goff + 90 * idx + 17);

            auto g_xxyy_xx = primBuffer.data(goff + 90 * idx + 18);

            auto g_xxyy_xy = primBuffer.data(goff + 90 * idx + 19);

            auto g_xxyy_xz = primBuffer.data(goff + 90 * idx + 20);

            auto g_xxyy_yy = primBuffer.data(goff + 90 * idx + 21);

            auto g_xxyy_yz = primBuffer.data(goff + 90 * idx + 22);

            auto g_xxyy_zz = primBuffer.data(goff + 90 * idx + 23);

            auto g_xxyz_xx = primBuffer.data(goff + 90 * idx + 24);

            auto g_xxyz_xy = primBuffer.data(goff + 90 * idx + 25);

            auto g_xxyz_xz = primBuffer.data(goff + 90 * idx + 26);

            auto g_xxyz_yy = primBuffer.data(goff + 90 * idx + 27);

            auto g_xxyz_yz = primBuffer.data(goff + 90 * idx + 28);

            auto g_xxyz_zz = primBuffer.data(goff + 90 * idx + 29);

            auto g_xxzz_xx = primBuffer.data(goff + 90 * idx + 30);

            auto g_xxzz_xy = primBuffer.data(goff + 90 * idx + 31);

            auto g_xxzz_xz = primBuffer.data(goff + 90 * idx + 32);

            auto g_xxzz_yy = primBuffer.data(goff + 90 * idx + 33);

            auto g_xxzz_yz = primBuffer.data(goff + 90 * idx + 34);

            auto g_xxzz_zz = primBuffer.data(goff + 90 * idx + 35);

            auto g_xyyy_xx = primBuffer.data(goff + 90 * idx + 36);

            auto g_xyyy_xy = primBuffer.data(goff + 90 * idx + 37);

            auto g_xyyy_xz = primBuffer.data(goff + 90 * idx + 38);

            auto g_xyyy_yy = primBuffer.data(goff + 90 * idx + 39);

            auto g_xyyy_yz = primBuffer.data(goff + 90 * idx + 40);

            auto g_xyyy_zz = primBuffer.data(goff + 90 * idx + 41);

            auto g_xyyz_xx = primBuffer.data(goff + 90 * idx + 42);

            auto g_xyyz_xy = primBuffer.data(goff + 90 * idx + 43);

            auto g_xyyz_xz = primBuffer.data(goff + 90 * idx + 44);

            auto g_xyyz_yy = primBuffer.data(goff + 90 * idx + 45);

            auto g_xyyz_yz = primBuffer.data(goff + 90 * idx + 46);

            auto g_xyyz_zz = primBuffer.data(goff + 90 * idx + 47);

            auto g_xyzz_xx = primBuffer.data(goff + 90 * idx + 48);

            auto g_xyzz_xy = primBuffer.data(goff + 90 * idx + 49);

            auto g_xyzz_xz = primBuffer.data(goff + 90 * idx + 50);

            auto g_xyzz_yy = primBuffer.data(goff + 90 * idx + 51);

            auto g_xyzz_yz = primBuffer.data(goff + 90 * idx + 52);

            auto g_xyzz_zz = primBuffer.data(goff + 90 * idx + 53);

            auto g_xzzz_xx = primBuffer.data(goff + 90 * idx + 54);

            auto g_xzzz_xy = primBuffer.data(goff + 90 * idx + 55);

            auto g_xzzz_xz = primBuffer.data(goff + 90 * idx + 56);

            auto g_xzzz_yy = primBuffer.data(goff + 90 * idx + 57);

            auto g_xzzz_yz = primBuffer.data(goff + 90 * idx + 58);

            auto g_xzzz_zz = primBuffer.data(goff + 90 * idx + 59);

            auto g_yyyy_xx = primBuffer.data(goff + 90 * idx + 60);

            auto g_yyyy_xy = primBuffer.data(goff + 90 * idx + 61);

            auto g_yyyy_xz = primBuffer.data(goff + 90 * idx + 62);

            auto g_yyyy_yy = primBuffer.data(goff + 90 * idx + 63);

            auto g_yyyy_yz = primBuffer.data(goff + 90 * idx + 64);

            auto g_yyyy_zz = primBuffer.data(goff + 90 * idx + 65);

            auto g_yyyz_xx = primBuffer.data(goff + 90 * idx + 66);

            auto g_yyyz_xy = primBuffer.data(goff + 90 * idx + 67);

            auto g_yyyz_xz = primBuffer.data(goff + 90 * idx + 68);

            auto g_yyyz_yy = primBuffer.data(goff + 90 * idx + 69);

            auto g_yyyz_yz = primBuffer.data(goff + 90 * idx + 70);

            auto g_yyyz_zz = primBuffer.data(goff + 90 * idx + 71);

            auto g_yyzz_xx = primBuffer.data(goff + 90 * idx + 72);

            auto g_yyzz_xy = primBuffer.data(goff + 90 * idx + 73);

            auto g_yyzz_xz = primBuffer.data(goff + 90 * idx + 74);

            auto g_yyzz_yy = primBuffer.data(goff + 90 * idx + 75);

            auto g_yyzz_yz = primBuffer.data(goff + 90 * idx + 76);

            auto g_yyzz_zz = primBuffer.data(goff + 90 * idx + 77);

            auto g_yzzz_xx = primBuffer.data(goff + 90 * idx + 78);

            auto g_yzzz_xy = primBuffer.data(goff + 90 * idx + 79);

            auto g_yzzz_xz = primBuffer.data(goff + 90 * idx + 80);

            auto g_yzzz_yy = primBuffer.data(goff + 90 * idx + 81);

            auto g_yzzz_yz = primBuffer.data(goff + 90 * idx + 82);

            auto g_yzzz_zz = primBuffer.data(goff + 90 * idx + 83);

            auto g_zzzz_xx = primBuffer.data(goff + 90 * idx + 84);

            auto g_zzzz_xy = primBuffer.data(goff + 90 * idx + 85);

            auto g_zzzz_xz = primBuffer.data(goff + 90 * idx + 86);

            auto g_zzzz_yy = primBuffer.data(goff + 90 * idx + 87);

            auto g_zzzz_yz = primBuffer.data(goff + 90 * idx + 88);

            auto g_zzzz_zz = primBuffer.data(goff + 90 * idx + 89);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xxx_x, gk_xxx_y,\
                                         gk_xxx_z, gk_xxy_x, gk_xxy_y, gk_xxy_z,\
                                         gk_xxz_x, gk_xxz_y, gk_xxz_z, gk_xyy_x,\
                                         gk_xyy_y, gk_xyy_z, gk_xyz_x, gk_xyz_y,\
                                         gk_xyz_z, gk_xzz_x, gk_xzz_y, gk_xzz_z,\
                                         gk_yyy_x, gk_yyy_y, gk_yyy_z, gk_yyz_x,\
                                         gk_yyz_y, gk_yyz_z, gk_yzz_x, gk_yzz_y,\
                                         gk_yzz_z, gk_zzz_x, gk_zzz_y, gk_zzz_z,\
                                         g20_xx_xx, g20_xx_xy, g20_xx_xz, g20_xx_yy,\
                                         g20_xx_yz, g20_xx_zz, g20_xy_xx, g20_xy_xy,\
                                         g20_xy_xz, g20_xy_yy, g20_xy_yz, g20_xy_zz,\
                                         g20_xz_xx, g20_xz_xy, g20_xz_xz, g20_xz_yy,\
                                         g20_xz_yz, g20_xz_zz, g20_yy_xx, g20_yy_xy,\
                                         g20_yy_xz, g20_yy_yy, g20_yy_yz, g20_yy_zz,\
                                         g20_yz_xx, g20_yz_xy, g20_yz_xz, g20_yz_yy,\
                                         g20_yz_yz, g20_yz_zz, g20_zz_xx, g20_zz_xy,\
                                         g20_zz_xz, g20_zz_yy, g20_zz_yz, g20_zz_zz,\
                                         g21_xx_xx, g21_xx_xy, g21_xx_xz, g21_xx_yy,\
                                         g21_xx_yz, g21_xx_zz, g21_xy_xx, g21_xy_xy,\
                                         g21_xy_xz, g21_xy_yy, g21_xy_yz, g21_xy_zz,\
                                         g21_xz_xx, g21_xz_xy, g21_xz_xz, g21_xz_yy,\
                                         g21_xz_yz, g21_xz_zz, g21_yy_xx, g21_yy_xy,\
                                         g21_yy_xz, g21_yy_yy, g21_yy_yz, g21_yy_zz,\
                                         g21_yz_xx, g21_yz_xy, g21_yz_xz, g21_yz_yy,\
                                         g21_yz_yz, g21_yz_zz, g21_zz_xx, g21_zz_xy,\
                                         g21_zz_xz, g21_zz_yy, g21_zz_yz, g21_zz_zz,\
                                         g1_xxx_xx, g1_xxx_xy, g1_xxx_xz, g1_xxx_yy,\
                                         g1_xxx_yz, g1_xxx_zz, g1_xxy_xx, g1_xxy_xy,\
                                         g1_xxy_xz, g1_xxy_yy, g1_xxy_yz, g1_xxy_zz,\
                                         g1_xxz_xx, g1_xxz_xy, g1_xxz_xz, g1_xxz_yy,\
                                         g1_xxz_yz, g1_xxz_zz, g1_xyy_xx, g1_xyy_xy,\
                                         g1_xyy_xz, g1_xyy_yy, g1_xyy_yz, g1_xyy_zz,\
                                         g1_xyz_xx, g1_xyz_xy, g1_xyz_xz, g1_xyz_yy,\
                                         g1_xyz_yz, g1_xyz_zz, g1_xzz_xx, g1_xzz_xy,\
                                         g1_xzz_xz, g1_xzz_yy, g1_xzz_yz, g1_xzz_zz,\
                                         g1_yyy_xx, g1_yyy_xy, g1_yyy_xz, g1_yyy_yy,\
                                         g1_yyy_yz, g1_yyy_zz, g1_yyz_xx, g1_yyz_xy,\
                                         g1_yyz_xz, g1_yyz_yy, g1_yyz_yz, g1_yyz_zz,\
                                         g1_yzz_xx, g1_yzz_xy, g1_yzz_xz, g1_yzz_yy,\
                                         g1_yzz_yz, g1_yzz_zz, g1_zzz_xx, g1_zzz_xy,\
                                         g1_zzz_xz, g1_zzz_yy, g1_zzz_yz, g1_zzz_zz,\
                                         g_xxxx_xx, g_xxxx_xy, g_xxxx_xz, g_xxxx_yy,\
                                         g_xxxx_yz, g_xxxx_zz, g_xxxy_xx, g_xxxy_xy,\
                                         g_xxxy_xz, g_xxxy_yy, g_xxxy_yz, g_xxxy_zz,\
                                         g_xxxz_xx, g_xxxz_xy, g_xxxz_xz, g_xxxz_yy,\
                                         g_xxxz_yz, g_xxxz_zz, g_xxyy_xx, g_xxyy_xy,\
                                         g_xxyy_xz, g_xxyy_yy, g_xxyy_yz, g_xxyy_zz,\
                                         g_xxyz_xx, g_xxyz_xy, g_xxyz_xz, g_xxyz_yy,\
                                         g_xxyz_yz, g_xxyz_zz, g_xxzz_xx, g_xxzz_xy,\
                                         g_xxzz_xz, g_xxzz_yy, g_xxzz_yz, g_xxzz_zz,\
                                         g_xyyy_xx, g_xyyy_xy, g_xyyy_xz, g_xyyy_yy,\
                                         g_xyyy_yz, g_xyyy_zz, g_xyyz_xx, g_xyyz_xy,\
                                         g_xyyz_xz, g_xyyz_yy, g_xyyz_yz, g_xyyz_zz,\
                                         g_xyzz_xx, g_xyzz_xy, g_xyzz_xz, g_xyzz_yy,\
                                         g_xyzz_yz, g_xyzz_zz, g_xzzz_xx, g_xzzz_xy,\
                                         g_xzzz_xz, g_xzzz_yy, g_xzzz_yz, g_xzzz_zz,\
                                         g_yyyy_xx, g_yyyy_xy, g_yyyy_xz, g_yyyy_yy,\
                                         g_yyyy_yz, g_yyyy_zz, g_yyyz_xx, g_yyyz_xy,\
                                         g_yyyz_xz, g_yyyz_yy, g_yyyz_yz, g_yyyz_zz,\
                                         g_yyzz_xx, g_yyzz_xy, g_yyzz_xz, g_yyzz_yy,\
                                         g_yyzz_yz, g_yyzz_zz, g_yzzz_xx, g_yzzz_xy,\
                                         g_yzzz_xz, g_yzzz_yy, g_yzzz_yz, g_yzzz_zz,\
                                         g_zzzz_xx, g_zzzz_xy, g_zzzz_xz, g_zzzz_yy,\
                                         g_zzzz_yz, g_zzzz_zz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxxx_xx[k] = fr * g1_xxx_xx[k] + f2g * (3.0 * g20_xx_xx[k] - 3.0 * fgz * g21_xx_xx[k]) + 2.0 * f2t * gk_xxx_x[k];

                g_xxxx_xy[k] = fr * g1_xxx_xy[k] + f2g * (3.0 * g20_xx_xy[k] - 3.0 * fgz * g21_xx_xy[k]) + f2t * gk_xxx_y[k];

                g_xxxx_xz[k] = fr * g1_xxx_xz[k] + f2g * (3.0 * g20_xx_xz[k] - 3.0 * fgz * g21_xx_xz[k]) + f2t * gk_xxx_z[k];

                g_xxxx_yy[k] = fr * g1_xxx_yy[k] + f2g * (3.0 * g20_xx_yy[k] - 3.0 * fgz * g21_xx_yy[k]);

                g_xxxx_yz[k] = fr * g1_xxx_yz[k] + f2g * (3.0 * g20_xx_yz[k] - 3.0 * fgz * g21_xx_yz[k]);

                g_xxxx_zz[k] = fr * g1_xxx_zz[k] + f2g * (3.0 * g20_xx_zz[k] - 3.0 * fgz * g21_xx_zz[k]);

                g_xxxy_xx[k] = fr * g1_xxy_xx[k] + f2g * (2.0 * g20_xy_xx[k] - 2.0 * fgz * g21_xy_xx[k]) + 2.0 * f2t * gk_xxy_x[k];

                g_xxxy_xy[k] = fr * g1_xxy_xy[k] + f2g * (2.0 * g20_xy_xy[k] - 2.0 * fgz * g21_xy_xy[k]) + f2t * gk_xxy_y[k];

                g_xxxy_xz[k] = fr * g1_xxy_xz[k] + f2g * (2.0 * g20_xy_xz[k] - 2.0 * fgz * g21_xy_xz[k]) + f2t * gk_xxy_z[k];

                g_xxxy_yy[k] = fr * g1_xxy_yy[k] + f2g * (2.0 * g20_xy_yy[k] - 2.0 * fgz * g21_xy_yy[k]);

                g_xxxy_yz[k] = fr * g1_xxy_yz[k] + f2g * (2.0 * g20_xy_yz[k] - 2.0 * fgz * g21_xy_yz[k]);

                g_xxxy_zz[k] = fr * g1_xxy_zz[k] + f2g * (2.0 * g20_xy_zz[k] - 2.0 * fgz * g21_xy_zz[k]);

                g_xxxz_xx[k] = fr * g1_xxz_xx[k] + f2g * (2.0 * g20_xz_xx[k] - 2.0 * fgz * g21_xz_xx[k]) + 2.0 * f2t * gk_xxz_x[k];

                g_xxxz_xy[k] = fr * g1_xxz_xy[k] + f2g * (2.0 * g20_xz_xy[k] - 2.0 * fgz * g21_xz_xy[k]) + f2t * gk_xxz_y[k];

                g_xxxz_xz[k] = fr * g1_xxz_xz[k] + f2g * (2.0 * g20_xz_xz[k] - 2.0 * fgz * g21_xz_xz[k]) + f2t * gk_xxz_z[k];

                g_xxxz_yy[k] = fr * g1_xxz_yy[k] + f2g * (2.0 * g20_xz_yy[k] - 2.0 * fgz * g21_xz_yy[k]);

                g_xxxz_yz[k] = fr * g1_xxz_yz[k] + f2g * (2.0 * g20_xz_yz[k] - 2.0 * fgz * g21_xz_yz[k]);

                g_xxxz_zz[k] = fr * g1_xxz_zz[k] + f2g * (2.0 * g20_xz_zz[k] - 2.0 * fgz * g21_xz_zz[k]);

                g_xxyy_xx[k] = fr * g1_xyy_xx[k] + f2g * (g20_yy_xx[k] - fgz * g21_yy_xx[k]) + 2.0 * f2t * gk_xyy_x[k];

                g_xxyy_xy[k] = fr * g1_xyy_xy[k] + f2g * (g20_yy_xy[k] - fgz * g21_yy_xy[k]) + f2t * gk_xyy_y[k];

                g_xxyy_xz[k] = fr * g1_xyy_xz[k] + f2g * (g20_yy_xz[k] - fgz * g21_yy_xz[k]) + f2t * gk_xyy_z[k];

                g_xxyy_yy[k] = fr * g1_xyy_yy[k] + f2g * (g20_yy_yy[k] - fgz * g21_yy_yy[k]);

                g_xxyy_yz[k] = fr * g1_xyy_yz[k] + f2g * (g20_yy_yz[k] - fgz * g21_yy_yz[k]);

                g_xxyy_zz[k] = fr * g1_xyy_zz[k] + f2g * (g20_yy_zz[k] - fgz * g21_yy_zz[k]);

                g_xxyz_xx[k] = fr * g1_xyz_xx[k] + f2g * (g20_yz_xx[k] - fgz * g21_yz_xx[k]) + 2.0 * f2t * gk_xyz_x[k];

                g_xxyz_xy[k] = fr * g1_xyz_xy[k] + f2g * (g20_yz_xy[k] - fgz * g21_yz_xy[k]) + f2t * gk_xyz_y[k];

                g_xxyz_xz[k] = fr * g1_xyz_xz[k] + f2g * (g20_yz_xz[k] - fgz * g21_yz_xz[k]) + f2t * gk_xyz_z[k];

                g_xxyz_yy[k] = fr * g1_xyz_yy[k] + f2g * (g20_yz_yy[k] - fgz * g21_yz_yy[k]);

                g_xxyz_yz[k] = fr * g1_xyz_yz[k] + f2g * (g20_yz_yz[k] - fgz * g21_yz_yz[k]);

                g_xxyz_zz[k] = fr * g1_xyz_zz[k] + f2g * (g20_yz_zz[k] - fgz * g21_yz_zz[k]);

                g_xxzz_xx[k] = fr * g1_xzz_xx[k] + f2g * (g20_zz_xx[k] - fgz * g21_zz_xx[k]) + 2.0 * f2t * gk_xzz_x[k];

                g_xxzz_xy[k] = fr * g1_xzz_xy[k] + f2g * (g20_zz_xy[k] - fgz * g21_zz_xy[k]) + f2t * gk_xzz_y[k];

                g_xxzz_xz[k] = fr * g1_xzz_xz[k] + f2g * (g20_zz_xz[k] - fgz * g21_zz_xz[k]) + f2t * gk_xzz_z[k];

                g_xxzz_yy[k] = fr * g1_xzz_yy[k] + f2g * (g20_zz_yy[k] - fgz * g21_zz_yy[k]);

                g_xxzz_yz[k] = fr * g1_xzz_yz[k] + f2g * (g20_zz_yz[k] - fgz * g21_zz_yz[k]);

                g_xxzz_zz[k] = fr * g1_xzz_zz[k] + f2g * (g20_zz_zz[k] - fgz * g21_zz_zz[k]);

                g_xyyy_xx[k] = fr * g1_yyy_xx[k] + 2.0 * f2t * gk_yyy_x[k];

                g_xyyy_xy[k] = fr * g1_yyy_xy[k] + f2t * gk_yyy_y[k];

                g_xyyy_xz[k] = fr * g1_yyy_xz[k] + f2t * gk_yyy_z[k];

                g_xyyy_yy[k] = fr * g1_yyy_yy[k];

                g_xyyy_yz[k] = fr * g1_yyy_yz[k];

                g_xyyy_zz[k] = fr * g1_yyy_zz[k];

                g_xyyz_xx[k] = fr * g1_yyz_xx[k] + 2.0 * f2t * gk_yyz_x[k];

                g_xyyz_xy[k] = fr * g1_yyz_xy[k] + f2t * gk_yyz_y[k];

                g_xyyz_xz[k] = fr * g1_yyz_xz[k] + f2t * gk_yyz_z[k];

                g_xyyz_yy[k] = fr * g1_yyz_yy[k];

                g_xyyz_yz[k] = fr * g1_yyz_yz[k];

                g_xyyz_zz[k] = fr * g1_yyz_zz[k];

                g_xyzz_xx[k] = fr * g1_yzz_xx[k] + 2.0 * f2t * gk_yzz_x[k];

                g_xyzz_xy[k] = fr * g1_yzz_xy[k] + f2t * gk_yzz_y[k];

                g_xyzz_xz[k] = fr * g1_yzz_xz[k] + f2t * gk_yzz_z[k];

                g_xyzz_yy[k] = fr * g1_yzz_yy[k];

                g_xyzz_yz[k] = fr * g1_yzz_yz[k];

                g_xyzz_zz[k] = fr * g1_yzz_zz[k];

                g_xzzz_xx[k] = fr * g1_zzz_xx[k] + 2.0 * f2t * gk_zzz_x[k];

                g_xzzz_xy[k] = fr * g1_zzz_xy[k] + f2t * gk_zzz_y[k];

                g_xzzz_xz[k] = fr * g1_zzz_xz[k] + f2t * gk_zzz_z[k];

                g_xzzz_yy[k] = fr * g1_zzz_yy[k];

                g_xzzz_yz[k] = fr * g1_zzz_yz[k];

                g_xzzz_zz[k] = fr * g1_zzz_zz[k];

                // leading y component

                fr = pay[k];

                g_yyyy_xx[k] = fr * g1_yyy_xx[k] + f2g * (3.0 * g20_yy_xx[k] - 3.0 * fgz * g21_yy_xx[k]);

                g_yyyy_xy[k] = fr * g1_yyy_xy[k] + f2g * (3.0 * g20_yy_xy[k] - 3.0 * fgz * g21_yy_xy[k]) + f2t * gk_yyy_x[k];

                g_yyyy_xz[k] = fr * g1_yyy_xz[k] + f2g * (3.0 * g20_yy_xz[k] - 3.0 * fgz * g21_yy_xz[k]);

                g_yyyy_yy[k] = fr * g1_yyy_yy[k] + f2g * (3.0 * g20_yy_yy[k] - 3.0 * fgz * g21_yy_yy[k]) + 2.0 * f2t * gk_yyy_y[k];

                g_yyyy_yz[k] = fr * g1_yyy_yz[k] + f2g * (3.0 * g20_yy_yz[k] - 3.0 * fgz * g21_yy_yz[k]) + f2t * gk_yyy_z[k];

                g_yyyy_zz[k] = fr * g1_yyy_zz[k] + f2g * (3.0 * g20_yy_zz[k] - 3.0 * fgz * g21_yy_zz[k]);

                g_yyyz_xx[k] = fr * g1_yyz_xx[k] + f2g * (2.0 * g20_yz_xx[k] - 2.0 * fgz * g21_yz_xx[k]);

                g_yyyz_xy[k] = fr * g1_yyz_xy[k] + f2g * (2.0 * g20_yz_xy[k] - 2.0 * fgz * g21_yz_xy[k]) + f2t * gk_yyz_x[k];

                g_yyyz_xz[k] = fr * g1_yyz_xz[k] + f2g * (2.0 * g20_yz_xz[k] - 2.0 * fgz * g21_yz_xz[k]);

                g_yyyz_yy[k] = fr * g1_yyz_yy[k] + f2g * (2.0 * g20_yz_yy[k] - 2.0 * fgz * g21_yz_yy[k]) + 2.0 * f2t * gk_yyz_y[k];

                g_yyyz_yz[k] = fr * g1_yyz_yz[k] + f2g * (2.0 * g20_yz_yz[k] - 2.0 * fgz * g21_yz_yz[k]) + f2t * gk_yyz_z[k];

                g_yyyz_zz[k] = fr * g1_yyz_zz[k] + f2g * (2.0 * g20_yz_zz[k] - 2.0 * fgz * g21_yz_zz[k]);

                g_yyzz_xx[k] = fr * g1_yzz_xx[k] + f2g * (g20_zz_xx[k] - fgz * g21_zz_xx[k]);

                g_yyzz_xy[k] = fr * g1_yzz_xy[k] + f2g * (g20_zz_xy[k] - fgz * g21_zz_xy[k]) + f2t * gk_yzz_x[k];

                g_yyzz_xz[k] = fr * g1_yzz_xz[k] + f2g * (g20_zz_xz[k] - fgz * g21_zz_xz[k]);

                g_yyzz_yy[k] = fr * g1_yzz_yy[k] + f2g * (g20_zz_yy[k] - fgz * g21_zz_yy[k]) + 2.0 * f2t * gk_yzz_y[k];

                g_yyzz_yz[k] = fr * g1_yzz_yz[k] + f2g * (g20_zz_yz[k] - fgz * g21_zz_yz[k]) + f2t * gk_yzz_z[k];

                g_yyzz_zz[k] = fr * g1_yzz_zz[k] + f2g * (g20_zz_zz[k] - fgz * g21_zz_zz[k]);

                g_yzzz_xx[k] = fr * g1_zzz_xx[k];

                g_yzzz_xy[k] = fr * g1_zzz_xy[k] + f2t * gk_zzz_x[k];

                g_yzzz_xz[k] = fr * g1_zzz_xz[k];

                g_yzzz_yy[k] = fr * g1_zzz_yy[k] + 2.0 * f2t * gk_zzz_y[k];

                g_yzzz_yz[k] = fr * g1_zzz_yz[k] + f2t * gk_zzz_z[k];

                g_yzzz_zz[k] = fr * g1_zzz_zz[k];

                // leading z component

                fr = paz[k];

                g_zzzz_xx[k] = fr * g1_zzz_xx[k] + f2g * (3.0 * g20_zz_xx[k] - 3.0 * fgz * g21_zz_xx[k]);

                g_zzzz_xy[k] = fr * g1_zzz_xy[k] + f2g * (3.0 * g20_zz_xy[k] - 3.0 * fgz * g21_zz_xy[k]);

                g_zzzz_xz[k] = fr * g1_zzz_xz[k] + f2g * (3.0 * g20_zz_xz[k] - 3.0 * fgz * g21_zz_xz[k]) + f2t * gk_zzz_x[k];

                g_zzzz_yy[k] = fr * g1_zzz_yy[k] + f2g * (3.0 * g20_zz_yy[k] - 3.0 * fgz * g21_zz_yy[k]);

                g_zzzz_yz[k] = fr * g1_zzz_yz[k] + f2g * (3.0 * g20_zz_yz[k] - 3.0 * fgz * g21_zz_yz[k]) + f2t * gk_zzz_y[k];

                g_zzzz_zz[k] = fr * g1_zzz_zz[k] + f2g * (3.0 * g20_zz_zz[k] - 3.0 * fgz * g21_zz_zz[k]) + 2.0 * f2t * gk_zzz_z[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForFG(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 3, 4);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {1, 4, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (D|g(r,r')|F)^(m+1) integrals

            auto gk_xx_xxx = primBuffer.data(gkoff + 60 * idx);

            auto gk_xx_xxy = primBuffer.data(gkoff + 60 * idx + 1);

            auto gk_xx_xxz = primBuffer.data(gkoff + 60 * idx + 2);

            auto gk_xx_xyy = primBuffer.data(gkoff + 60 * idx + 3);

            auto gk_xx_xyz = primBuffer.data(gkoff + 60 * idx + 4);

            auto gk_xx_xzz = primBuffer.data(gkoff + 60 * idx + 5);

            auto gk_xx_yyy = primBuffer.data(gkoff + 60 * idx + 6);

            auto gk_xx_yyz = primBuffer.data(gkoff + 60 * idx + 7);

            auto gk_xx_yzz = primBuffer.data(gkoff + 60 * idx + 8);

            auto gk_xx_zzz = primBuffer.data(gkoff + 60 * idx + 9);

            auto gk_xy_xxx = primBuffer.data(gkoff + 60 * idx + 10);

            auto gk_xy_xxy = primBuffer.data(gkoff + 60 * idx + 11);

            auto gk_xy_xxz = primBuffer.data(gkoff + 60 * idx + 12);

            auto gk_xy_xyy = primBuffer.data(gkoff + 60 * idx + 13);

            auto gk_xy_xyz = primBuffer.data(gkoff + 60 * idx + 14);

            auto gk_xy_xzz = primBuffer.data(gkoff + 60 * idx + 15);

            auto gk_xy_yyy = primBuffer.data(gkoff + 60 * idx + 16);

            auto gk_xy_yyz = primBuffer.data(gkoff + 60 * idx + 17);

            auto gk_xy_yzz = primBuffer.data(gkoff + 60 * idx + 18);

            auto gk_xy_zzz = primBuffer.data(gkoff + 60 * idx + 19);

            auto gk_xz_xxx = primBuffer.data(gkoff + 60 * idx + 20);

            auto gk_xz_xxy = primBuffer.data(gkoff + 60 * idx + 21);

            auto gk_xz_xxz = primBuffer.data(gkoff + 60 * idx + 22);

            auto gk_xz_xyy = primBuffer.data(gkoff + 60 * idx + 23);

            auto gk_xz_xyz = primBuffer.data(gkoff + 60 * idx + 24);

            auto gk_xz_xzz = primBuffer.data(gkoff + 60 * idx + 25);

            auto gk_xz_yyy = primBuffer.data(gkoff + 60 * idx + 26);

            auto gk_xz_yyz = primBuffer.data(gkoff + 60 * idx + 27);

            auto gk_xz_yzz = primBuffer.data(gkoff + 60 * idx + 28);

            auto gk_xz_zzz = primBuffer.data(gkoff + 60 * idx + 29);

            auto gk_yy_xxx = primBuffer.data(gkoff + 60 * idx + 30);

            auto gk_yy_xxy = primBuffer.data(gkoff + 60 * idx + 31);

            auto gk_yy_xxz = primBuffer.data(gkoff + 60 * idx + 32);

            auto gk_yy_xyy = primBuffer.data(gkoff + 60 * idx + 33);

            auto gk_yy_xyz = primBuffer.data(gkoff + 60 * idx + 34);

            auto gk_yy_xzz = primBuffer.data(gkoff + 60 * idx + 35);

            auto gk_yy_yyy = primBuffer.data(gkoff + 60 * idx + 36);

            auto gk_yy_yyz = primBuffer.data(gkoff + 60 * idx + 37);

            auto gk_yy_yzz = primBuffer.data(gkoff + 60 * idx + 38);

            auto gk_yy_zzz = primBuffer.data(gkoff + 60 * idx + 39);

            auto gk_yz_xxx = primBuffer.data(gkoff + 60 * idx + 40);

            auto gk_yz_xxy = primBuffer.data(gkoff + 60 * idx + 41);

            auto gk_yz_xxz = primBuffer.data(gkoff + 60 * idx + 42);

            auto gk_yz_xyy = primBuffer.data(gkoff + 60 * idx + 43);

            auto gk_yz_xyz = primBuffer.data(gkoff + 60 * idx + 44);

            auto gk_yz_xzz = primBuffer.data(gkoff + 60 * idx + 45);

            auto gk_yz_yyy = primBuffer.data(gkoff + 60 * idx + 46);

            auto gk_yz_yyz = primBuffer.data(gkoff + 60 * idx + 47);

            auto gk_yz_yzz = primBuffer.data(gkoff + 60 * idx + 48);

            auto gk_yz_zzz = primBuffer.data(gkoff + 60 * idx + 49);

            auto gk_zz_xxx = primBuffer.data(gkoff + 60 * idx + 50);

            auto gk_zz_xxy = primBuffer.data(gkoff + 60 * idx + 51);

            auto gk_zz_xxz = primBuffer.data(gkoff + 60 * idx + 52);

            auto gk_zz_xyy = primBuffer.data(gkoff + 60 * idx + 53);

            auto gk_zz_xyz = primBuffer.data(gkoff + 60 * idx + 54);

            auto gk_zz_xzz = primBuffer.data(gkoff + 60 * idx + 55);

            auto gk_zz_yyy = primBuffer.data(gkoff + 60 * idx + 56);

            auto gk_zz_yyz = primBuffer.data(gkoff + 60 * idx + 57);

            auto gk_zz_yzz = primBuffer.data(gkoff + 60 * idx + 58);

            auto gk_zz_zzz = primBuffer.data(gkoff + 60 * idx + 59);

            // set up pointers to (P|g(r,r')|G)^(m) integrals

            auto g20_x_xxxx = primBuffer.data(g20off + 45 * idx);

            auto g20_x_xxxy = primBuffer.data(g20off + 45 * idx + 1);

            auto g20_x_xxxz = primBuffer.data(g20off + 45 * idx + 2);

            auto g20_x_xxyy = primBuffer.data(g20off + 45 * idx + 3);

            auto g20_x_xxyz = primBuffer.data(g20off + 45 * idx + 4);

            auto g20_x_xxzz = primBuffer.data(g20off + 45 * idx + 5);

            auto g20_x_xyyy = primBuffer.data(g20off + 45 * idx + 6);

            auto g20_x_xyyz = primBuffer.data(g20off + 45 * idx + 7);

            auto g20_x_xyzz = primBuffer.data(g20off + 45 * idx + 8);

            auto g20_x_xzzz = primBuffer.data(g20off + 45 * idx + 9);

            auto g20_x_yyyy = primBuffer.data(g20off + 45 * idx + 10);

            auto g20_x_yyyz = primBuffer.data(g20off + 45 * idx + 11);

            auto g20_x_yyzz = primBuffer.data(g20off + 45 * idx + 12);

            auto g20_x_yzzz = primBuffer.data(g20off + 45 * idx + 13);

            auto g20_x_zzzz = primBuffer.data(g20off + 45 * idx + 14);

            auto g20_y_xxxx = primBuffer.data(g20off + 45 * idx + 15);

            auto g20_y_xxxy = primBuffer.data(g20off + 45 * idx + 16);

            auto g20_y_xxxz = primBuffer.data(g20off + 45 * idx + 17);

            auto g20_y_xxyy = primBuffer.data(g20off + 45 * idx + 18);

            auto g20_y_xxyz = primBuffer.data(g20off + 45 * idx + 19);

            auto g20_y_xxzz = primBuffer.data(g20off + 45 * idx + 20);

            auto g20_y_xyyy = primBuffer.data(g20off + 45 * idx + 21);

            auto g20_y_xyyz = primBuffer.data(g20off + 45 * idx + 22);

            auto g20_y_xyzz = primBuffer.data(g20off + 45 * idx + 23);

            auto g20_y_xzzz = primBuffer.data(g20off + 45 * idx + 24);

            auto g20_y_yyyy = primBuffer.data(g20off + 45 * idx + 25);

            auto g20_y_yyyz = primBuffer.data(g20off + 45 * idx + 26);

            auto g20_y_yyzz = primBuffer.data(g20off + 45 * idx + 27);

            auto g20_y_yzzz = primBuffer.data(g20off + 45 * idx + 28);

            auto g20_y_zzzz = primBuffer.data(g20off + 45 * idx + 29);

            auto g20_z_xxxx = primBuffer.data(g20off + 45 * idx + 30);

            auto g20_z_xxxy = primBuffer.data(g20off + 45 * idx + 31);

            auto g20_z_xxxz = primBuffer.data(g20off + 45 * idx + 32);

            auto g20_z_xxyy = primBuffer.data(g20off + 45 * idx + 33);

            auto g20_z_xxyz = primBuffer.data(g20off + 45 * idx + 34);

            auto g20_z_xxzz = primBuffer.data(g20off + 45 * idx + 35);

            auto g20_z_xyyy = primBuffer.data(g20off + 45 * idx + 36);

            auto g20_z_xyyz = primBuffer.data(g20off + 45 * idx + 37);

            auto g20_z_xyzz = primBuffer.data(g20off + 45 * idx + 38);

            auto g20_z_xzzz = primBuffer.data(g20off + 45 * idx + 39);

            auto g20_z_yyyy = primBuffer.data(g20off + 45 * idx + 40);

            auto g20_z_yyyz = primBuffer.data(g20off + 45 * idx + 41);

            auto g20_z_yyzz = primBuffer.data(g20off + 45 * idx + 42);

            auto g20_z_yzzz = primBuffer.data(g20off + 45 * idx + 43);

            auto g20_z_zzzz = primBuffer.data(g20off + 45 * idx + 44);

            // set up pointers to (P|g(r,r')|G)^(m+1) integrals

            auto g21_x_xxxx = primBuffer.data(g21off + 45 * idx);

            auto g21_x_xxxy = primBuffer.data(g21off + 45 * idx + 1);

            auto g21_x_xxxz = primBuffer.data(g21off + 45 * idx + 2);

            auto g21_x_xxyy = primBuffer.data(g21off + 45 * idx + 3);

            auto g21_x_xxyz = primBuffer.data(g21off + 45 * idx + 4);

            auto g21_x_xxzz = primBuffer.data(g21off + 45 * idx + 5);

            auto g21_x_xyyy = primBuffer.data(g21off + 45 * idx + 6);

            auto g21_x_xyyz = primBuffer.data(g21off + 45 * idx + 7);

            auto g21_x_xyzz = primBuffer.data(g21off + 45 * idx + 8);

            auto g21_x_xzzz = primBuffer.data(g21off + 45 * idx + 9);

            auto g21_x_yyyy = primBuffer.data(g21off + 45 * idx + 10);

            auto g21_x_yyyz = primBuffer.data(g21off + 45 * idx + 11);

            auto g21_x_yyzz = primBuffer.data(g21off + 45 * idx + 12);

            auto g21_x_yzzz = primBuffer.data(g21off + 45 * idx + 13);

            auto g21_x_zzzz = primBuffer.data(g21off + 45 * idx + 14);

            auto g21_y_xxxx = primBuffer.data(g21off + 45 * idx + 15);

            auto g21_y_xxxy = primBuffer.data(g21off + 45 * idx + 16);

            auto g21_y_xxxz = primBuffer.data(g21off + 45 * idx + 17);

            auto g21_y_xxyy = primBuffer.data(g21off + 45 * idx + 18);

            auto g21_y_xxyz = primBuffer.data(g21off + 45 * idx + 19);

            auto g21_y_xxzz = primBuffer.data(g21off + 45 * idx + 20);

            auto g21_y_xyyy = primBuffer.data(g21off + 45 * idx + 21);

            auto g21_y_xyyz = primBuffer.data(g21off + 45 * idx + 22);

            auto g21_y_xyzz = primBuffer.data(g21off + 45 * idx + 23);

            auto g21_y_xzzz = primBuffer.data(g21off + 45 * idx + 24);

            auto g21_y_yyyy = primBuffer.data(g21off + 45 * idx + 25);

            auto g21_y_yyyz = primBuffer.data(g21off + 45 * idx + 26);

            auto g21_y_yyzz = primBuffer.data(g21off + 45 * idx + 27);

            auto g21_y_yzzz = primBuffer.data(g21off + 45 * idx + 28);

            auto g21_y_zzzz = primBuffer.data(g21off + 45 * idx + 29);

            auto g21_z_xxxx = primBuffer.data(g21off + 45 * idx + 30);

            auto g21_z_xxxy = primBuffer.data(g21off + 45 * idx + 31);

            auto g21_z_xxxz = primBuffer.data(g21off + 45 * idx + 32);

            auto g21_z_xxyy = primBuffer.data(g21off + 45 * idx + 33);

            auto g21_z_xxyz = primBuffer.data(g21off + 45 * idx + 34);

            auto g21_z_xxzz = primBuffer.data(g21off + 45 * idx + 35);

            auto g21_z_xyyy = primBuffer.data(g21off + 45 * idx + 36);

            auto g21_z_xyyz = primBuffer.data(g21off + 45 * idx + 37);

            auto g21_z_xyzz = primBuffer.data(g21off + 45 * idx + 38);

            auto g21_z_xzzz = primBuffer.data(g21off + 45 * idx + 39);

            auto g21_z_yyyy = primBuffer.data(g21off + 45 * idx + 40);

            auto g21_z_yyyz = primBuffer.data(g21off + 45 * idx + 41);

            auto g21_z_yyzz = primBuffer.data(g21off + 45 * idx + 42);

            auto g21_z_yzzz = primBuffer.data(g21off + 45 * idx + 43);

            auto g21_z_zzzz = primBuffer.data(g21off + 45 * idx + 44);

            // set up pointers to (D|g(r,r')|G)^(m+1) integrals

            auto g1_xx_xxxx = primBuffer.data(g1off + 90 * idx);

            auto g1_xx_xxxy = primBuffer.data(g1off + 90 * idx + 1);

            auto g1_xx_xxxz = primBuffer.data(g1off + 90 * idx + 2);

            auto g1_xx_xxyy = primBuffer.data(g1off + 90 * idx + 3);

            auto g1_xx_xxyz = primBuffer.data(g1off + 90 * idx + 4);

            auto g1_xx_xxzz = primBuffer.data(g1off + 90 * idx + 5);

            auto g1_xx_xyyy = primBuffer.data(g1off + 90 * idx + 6);

            auto g1_xx_xyyz = primBuffer.data(g1off + 90 * idx + 7);

            auto g1_xx_xyzz = primBuffer.data(g1off + 90 * idx + 8);

            auto g1_xx_xzzz = primBuffer.data(g1off + 90 * idx + 9);

            auto g1_xx_yyyy = primBuffer.data(g1off + 90 * idx + 10);

            auto g1_xx_yyyz = primBuffer.data(g1off + 90 * idx + 11);

            auto g1_xx_yyzz = primBuffer.data(g1off + 90 * idx + 12);

            auto g1_xx_yzzz = primBuffer.data(g1off + 90 * idx + 13);

            auto g1_xx_zzzz = primBuffer.data(g1off + 90 * idx + 14);

            auto g1_xy_xxxx = primBuffer.data(g1off + 90 * idx + 15);

            auto g1_xy_xxxy = primBuffer.data(g1off + 90 * idx + 16);

            auto g1_xy_xxxz = primBuffer.data(g1off + 90 * idx + 17);

            auto g1_xy_xxyy = primBuffer.data(g1off + 90 * idx + 18);

            auto g1_xy_xxyz = primBuffer.data(g1off + 90 * idx + 19);

            auto g1_xy_xxzz = primBuffer.data(g1off + 90 * idx + 20);

            auto g1_xy_xyyy = primBuffer.data(g1off + 90 * idx + 21);

            auto g1_xy_xyyz = primBuffer.data(g1off + 90 * idx + 22);

            auto g1_xy_xyzz = primBuffer.data(g1off + 90 * idx + 23);

            auto g1_xy_xzzz = primBuffer.data(g1off + 90 * idx + 24);

            auto g1_xy_yyyy = primBuffer.data(g1off + 90 * idx + 25);

            auto g1_xy_yyyz = primBuffer.data(g1off + 90 * idx + 26);

            auto g1_xy_yyzz = primBuffer.data(g1off + 90 * idx + 27);

            auto g1_xy_yzzz = primBuffer.data(g1off + 90 * idx + 28);

            auto g1_xy_zzzz = primBuffer.data(g1off + 90 * idx + 29);

            auto g1_xz_xxxx = primBuffer.data(g1off + 90 * idx + 30);

            auto g1_xz_xxxy = primBuffer.data(g1off + 90 * idx + 31);

            auto g1_xz_xxxz = primBuffer.data(g1off + 90 * idx + 32);

            auto g1_xz_xxyy = primBuffer.data(g1off + 90 * idx + 33);

            auto g1_xz_xxyz = primBuffer.data(g1off + 90 * idx + 34);

            auto g1_xz_xxzz = primBuffer.data(g1off + 90 * idx + 35);

            auto g1_xz_xyyy = primBuffer.data(g1off + 90 * idx + 36);

            auto g1_xz_xyyz = primBuffer.data(g1off + 90 * idx + 37);

            auto g1_xz_xyzz = primBuffer.data(g1off + 90 * idx + 38);

            auto g1_xz_xzzz = primBuffer.data(g1off + 90 * idx + 39);

            auto g1_xz_yyyy = primBuffer.data(g1off + 90 * idx + 40);

            auto g1_xz_yyyz = primBuffer.data(g1off + 90 * idx + 41);

            auto g1_xz_yyzz = primBuffer.data(g1off + 90 * idx + 42);

            auto g1_xz_yzzz = primBuffer.data(g1off + 90 * idx + 43);

            auto g1_xz_zzzz = primBuffer.data(g1off + 90 * idx + 44);

            auto g1_yy_xxxx = primBuffer.data(g1off + 90 * idx + 45);

            auto g1_yy_xxxy = primBuffer.data(g1off + 90 * idx + 46);

            auto g1_yy_xxxz = primBuffer.data(g1off + 90 * idx + 47);

            auto g1_yy_xxyy = primBuffer.data(g1off + 90 * idx + 48);

            auto g1_yy_xxyz = primBuffer.data(g1off + 90 * idx + 49);

            auto g1_yy_xxzz = primBuffer.data(g1off + 90 * idx + 50);

            auto g1_yy_xyyy = primBuffer.data(g1off + 90 * idx + 51);

            auto g1_yy_xyyz = primBuffer.data(g1off + 90 * idx + 52);

            auto g1_yy_xyzz = primBuffer.data(g1off + 90 * idx + 53);

            auto g1_yy_xzzz = primBuffer.data(g1off + 90 * idx + 54);

            auto g1_yy_yyyy = primBuffer.data(g1off + 90 * idx + 55);

            auto g1_yy_yyyz = primBuffer.data(g1off + 90 * idx + 56);

            auto g1_yy_yyzz = primBuffer.data(g1off + 90 * idx + 57);

            auto g1_yy_yzzz = primBuffer.data(g1off + 90 * idx + 58);

            auto g1_yy_zzzz = primBuffer.data(g1off + 90 * idx + 59);

            auto g1_yz_xxxx = primBuffer.data(g1off + 90 * idx + 60);

            auto g1_yz_xxxy = primBuffer.data(g1off + 90 * idx + 61);

            auto g1_yz_xxxz = primBuffer.data(g1off + 90 * idx + 62);

            auto g1_yz_xxyy = primBuffer.data(g1off + 90 * idx + 63);

            auto g1_yz_xxyz = primBuffer.data(g1off + 90 * idx + 64);

            auto g1_yz_xxzz = primBuffer.data(g1off + 90 * idx + 65);

            auto g1_yz_xyyy = primBuffer.data(g1off + 90 * idx + 66);

            auto g1_yz_xyyz = primBuffer.data(g1off + 90 * idx + 67);

            auto g1_yz_xyzz = primBuffer.data(g1off + 90 * idx + 68);

            auto g1_yz_xzzz = primBuffer.data(g1off + 90 * idx + 69);

            auto g1_yz_yyyy = primBuffer.data(g1off + 90 * idx + 70);

            auto g1_yz_yyyz = primBuffer.data(g1off + 90 * idx + 71);

            auto g1_yz_yyzz = primBuffer.data(g1off + 90 * idx + 72);

            auto g1_yz_yzzz = primBuffer.data(g1off + 90 * idx + 73);

            auto g1_yz_zzzz = primBuffer.data(g1off + 90 * idx + 74);

            auto g1_zz_xxxx = primBuffer.data(g1off + 90 * idx + 75);

            auto g1_zz_xxxy = primBuffer.data(g1off + 90 * idx + 76);

            auto g1_zz_xxxz = primBuffer.data(g1off + 90 * idx + 77);

            auto g1_zz_xxyy = primBuffer.data(g1off + 90 * idx + 78);

            auto g1_zz_xxyz = primBuffer.data(g1off + 90 * idx + 79);

            auto g1_zz_xxzz = primBuffer.data(g1off + 90 * idx + 80);

            auto g1_zz_xyyy = primBuffer.data(g1off + 90 * idx + 81);

            auto g1_zz_xyyz = primBuffer.data(g1off + 90 * idx + 82);

            auto g1_zz_xyzz = primBuffer.data(g1off + 90 * idx + 83);

            auto g1_zz_xzzz = primBuffer.data(g1off + 90 * idx + 84);

            auto g1_zz_yyyy = primBuffer.data(g1off + 90 * idx + 85);

            auto g1_zz_yyyz = primBuffer.data(g1off + 90 * idx + 86);

            auto g1_zz_yyzz = primBuffer.data(g1off + 90 * idx + 87);

            auto g1_zz_yzzz = primBuffer.data(g1off + 90 * idx + 88);

            auto g1_zz_zzzz = primBuffer.data(g1off + 90 * idx + 89);

            // set up pointers to (F|g(r,r')|G)^(m) integrals

            auto g_xxx_xxxx = primBuffer.data(goff + 150 * idx);

            auto g_xxx_xxxy = primBuffer.data(goff + 150 * idx + 1);

            auto g_xxx_xxxz = primBuffer.data(goff + 150 * idx + 2);

            auto g_xxx_xxyy = primBuffer.data(goff + 150 * idx + 3);

            auto g_xxx_xxyz = primBuffer.data(goff + 150 * idx + 4);

            auto g_xxx_xxzz = primBuffer.data(goff + 150 * idx + 5);

            auto g_xxx_xyyy = primBuffer.data(goff + 150 * idx + 6);

            auto g_xxx_xyyz = primBuffer.data(goff + 150 * idx + 7);

            auto g_xxx_xyzz = primBuffer.data(goff + 150 * idx + 8);

            auto g_xxx_xzzz = primBuffer.data(goff + 150 * idx + 9);

            auto g_xxx_yyyy = primBuffer.data(goff + 150 * idx + 10);

            auto g_xxx_yyyz = primBuffer.data(goff + 150 * idx + 11);

            auto g_xxx_yyzz = primBuffer.data(goff + 150 * idx + 12);

            auto g_xxx_yzzz = primBuffer.data(goff + 150 * idx + 13);

            auto g_xxx_zzzz = primBuffer.data(goff + 150 * idx + 14);

            auto g_xxy_xxxx = primBuffer.data(goff + 150 * idx + 15);

            auto g_xxy_xxxy = primBuffer.data(goff + 150 * idx + 16);

            auto g_xxy_xxxz = primBuffer.data(goff + 150 * idx + 17);

            auto g_xxy_xxyy = primBuffer.data(goff + 150 * idx + 18);

            auto g_xxy_xxyz = primBuffer.data(goff + 150 * idx + 19);

            auto g_xxy_xxzz = primBuffer.data(goff + 150 * idx + 20);

            auto g_xxy_xyyy = primBuffer.data(goff + 150 * idx + 21);

            auto g_xxy_xyyz = primBuffer.data(goff + 150 * idx + 22);

            auto g_xxy_xyzz = primBuffer.data(goff + 150 * idx + 23);

            auto g_xxy_xzzz = primBuffer.data(goff + 150 * idx + 24);

            auto g_xxy_yyyy = primBuffer.data(goff + 150 * idx + 25);

            auto g_xxy_yyyz = primBuffer.data(goff + 150 * idx + 26);

            auto g_xxy_yyzz = primBuffer.data(goff + 150 * idx + 27);

            auto g_xxy_yzzz = primBuffer.data(goff + 150 * idx + 28);

            auto g_xxy_zzzz = primBuffer.data(goff + 150 * idx + 29);

            auto g_xxz_xxxx = primBuffer.data(goff + 150 * idx + 30);

            auto g_xxz_xxxy = primBuffer.data(goff + 150 * idx + 31);

            auto g_xxz_xxxz = primBuffer.data(goff + 150 * idx + 32);

            auto g_xxz_xxyy = primBuffer.data(goff + 150 * idx + 33);

            auto g_xxz_xxyz = primBuffer.data(goff + 150 * idx + 34);

            auto g_xxz_xxzz = primBuffer.data(goff + 150 * idx + 35);

            auto g_xxz_xyyy = primBuffer.data(goff + 150 * idx + 36);

            auto g_xxz_xyyz = primBuffer.data(goff + 150 * idx + 37);

            auto g_xxz_xyzz = primBuffer.data(goff + 150 * idx + 38);

            auto g_xxz_xzzz = primBuffer.data(goff + 150 * idx + 39);

            auto g_xxz_yyyy = primBuffer.data(goff + 150 * idx + 40);

            auto g_xxz_yyyz = primBuffer.data(goff + 150 * idx + 41);

            auto g_xxz_yyzz = primBuffer.data(goff + 150 * idx + 42);

            auto g_xxz_yzzz = primBuffer.data(goff + 150 * idx + 43);

            auto g_xxz_zzzz = primBuffer.data(goff + 150 * idx + 44);

            auto g_xyy_xxxx = primBuffer.data(goff + 150 * idx + 45);

            auto g_xyy_xxxy = primBuffer.data(goff + 150 * idx + 46);

            auto g_xyy_xxxz = primBuffer.data(goff + 150 * idx + 47);

            auto g_xyy_xxyy = primBuffer.data(goff + 150 * idx + 48);

            auto g_xyy_xxyz = primBuffer.data(goff + 150 * idx + 49);

            auto g_xyy_xxzz = primBuffer.data(goff + 150 * idx + 50);

            auto g_xyy_xyyy = primBuffer.data(goff + 150 * idx + 51);

            auto g_xyy_xyyz = primBuffer.data(goff + 150 * idx + 52);

            auto g_xyy_xyzz = primBuffer.data(goff + 150 * idx + 53);

            auto g_xyy_xzzz = primBuffer.data(goff + 150 * idx + 54);

            auto g_xyy_yyyy = primBuffer.data(goff + 150 * idx + 55);

            auto g_xyy_yyyz = primBuffer.data(goff + 150 * idx + 56);

            auto g_xyy_yyzz = primBuffer.data(goff + 150 * idx + 57);

            auto g_xyy_yzzz = primBuffer.data(goff + 150 * idx + 58);

            auto g_xyy_zzzz = primBuffer.data(goff + 150 * idx + 59);

            auto g_xyz_xxxx = primBuffer.data(goff + 150 * idx + 60);

            auto g_xyz_xxxy = primBuffer.data(goff + 150 * idx + 61);

            auto g_xyz_xxxz = primBuffer.data(goff + 150 * idx + 62);

            auto g_xyz_xxyy = primBuffer.data(goff + 150 * idx + 63);

            auto g_xyz_xxyz = primBuffer.data(goff + 150 * idx + 64);

            auto g_xyz_xxzz = primBuffer.data(goff + 150 * idx + 65);

            auto g_xyz_xyyy = primBuffer.data(goff + 150 * idx + 66);

            auto g_xyz_xyyz = primBuffer.data(goff + 150 * idx + 67);

            auto g_xyz_xyzz = primBuffer.data(goff + 150 * idx + 68);

            auto g_xyz_xzzz = primBuffer.data(goff + 150 * idx + 69);

            auto g_xyz_yyyy = primBuffer.data(goff + 150 * idx + 70);

            auto g_xyz_yyyz = primBuffer.data(goff + 150 * idx + 71);

            auto g_xyz_yyzz = primBuffer.data(goff + 150 * idx + 72);

            auto g_xyz_yzzz = primBuffer.data(goff + 150 * idx + 73);

            auto g_xyz_zzzz = primBuffer.data(goff + 150 * idx + 74);

            auto g_xzz_xxxx = primBuffer.data(goff + 150 * idx + 75);

            auto g_xzz_xxxy = primBuffer.data(goff + 150 * idx + 76);

            auto g_xzz_xxxz = primBuffer.data(goff + 150 * idx + 77);

            auto g_xzz_xxyy = primBuffer.data(goff + 150 * idx + 78);

            auto g_xzz_xxyz = primBuffer.data(goff + 150 * idx + 79);

            auto g_xzz_xxzz = primBuffer.data(goff + 150 * idx + 80);

            auto g_xzz_xyyy = primBuffer.data(goff + 150 * idx + 81);

            auto g_xzz_xyyz = primBuffer.data(goff + 150 * idx + 82);

            auto g_xzz_xyzz = primBuffer.data(goff + 150 * idx + 83);

            auto g_xzz_xzzz = primBuffer.data(goff + 150 * idx + 84);

            auto g_xzz_yyyy = primBuffer.data(goff + 150 * idx + 85);

            auto g_xzz_yyyz = primBuffer.data(goff + 150 * idx + 86);

            auto g_xzz_yyzz = primBuffer.data(goff + 150 * idx + 87);

            auto g_xzz_yzzz = primBuffer.data(goff + 150 * idx + 88);

            auto g_xzz_zzzz = primBuffer.data(goff + 150 * idx + 89);

            auto g_yyy_xxxx = primBuffer.data(goff + 150 * idx + 90);

            auto g_yyy_xxxy = primBuffer.data(goff + 150 * idx + 91);

            auto g_yyy_xxxz = primBuffer.data(goff + 150 * idx + 92);

            auto g_yyy_xxyy = primBuffer.data(goff + 150 * idx + 93);

            auto g_yyy_xxyz = primBuffer.data(goff + 150 * idx + 94);

            auto g_yyy_xxzz = primBuffer.data(goff + 150 * idx + 95);

            auto g_yyy_xyyy = primBuffer.data(goff + 150 * idx + 96);

            auto g_yyy_xyyz = primBuffer.data(goff + 150 * idx + 97);

            auto g_yyy_xyzz = primBuffer.data(goff + 150 * idx + 98);

            auto g_yyy_xzzz = primBuffer.data(goff + 150 * idx + 99);

            auto g_yyy_yyyy = primBuffer.data(goff + 150 * idx + 100);

            auto g_yyy_yyyz = primBuffer.data(goff + 150 * idx + 101);

            auto g_yyy_yyzz = primBuffer.data(goff + 150 * idx + 102);

            auto g_yyy_yzzz = primBuffer.data(goff + 150 * idx + 103);

            auto g_yyy_zzzz = primBuffer.data(goff + 150 * idx + 104);

            auto g_yyz_xxxx = primBuffer.data(goff + 150 * idx + 105);

            auto g_yyz_xxxy = primBuffer.data(goff + 150 * idx + 106);

            auto g_yyz_xxxz = primBuffer.data(goff + 150 * idx + 107);

            auto g_yyz_xxyy = primBuffer.data(goff + 150 * idx + 108);

            auto g_yyz_xxyz = primBuffer.data(goff + 150 * idx + 109);

            auto g_yyz_xxzz = primBuffer.data(goff + 150 * idx + 110);

            auto g_yyz_xyyy = primBuffer.data(goff + 150 * idx + 111);

            auto g_yyz_xyyz = primBuffer.data(goff + 150 * idx + 112);

            auto g_yyz_xyzz = primBuffer.data(goff + 150 * idx + 113);

            auto g_yyz_xzzz = primBuffer.data(goff + 150 * idx + 114);

            auto g_yyz_yyyy = primBuffer.data(goff + 150 * idx + 115);

            auto g_yyz_yyyz = primBuffer.data(goff + 150 * idx + 116);

            auto g_yyz_yyzz = primBuffer.data(goff + 150 * idx + 117);

            auto g_yyz_yzzz = primBuffer.data(goff + 150 * idx + 118);

            auto g_yyz_zzzz = primBuffer.data(goff + 150 * idx + 119);

            auto g_yzz_xxxx = primBuffer.data(goff + 150 * idx + 120);

            auto g_yzz_xxxy = primBuffer.data(goff + 150 * idx + 121);

            auto g_yzz_xxxz = primBuffer.data(goff + 150 * idx + 122);

            auto g_yzz_xxyy = primBuffer.data(goff + 150 * idx + 123);

            auto g_yzz_xxyz = primBuffer.data(goff + 150 * idx + 124);

            auto g_yzz_xxzz = primBuffer.data(goff + 150 * idx + 125);

            auto g_yzz_xyyy = primBuffer.data(goff + 150 * idx + 126);

            auto g_yzz_xyyz = primBuffer.data(goff + 150 * idx + 127);

            auto g_yzz_xyzz = primBuffer.data(goff + 150 * idx + 128);

            auto g_yzz_xzzz = primBuffer.data(goff + 150 * idx + 129);

            auto g_yzz_yyyy = primBuffer.data(goff + 150 * idx + 130);

            auto g_yzz_yyyz = primBuffer.data(goff + 150 * idx + 131);

            auto g_yzz_yyzz = primBuffer.data(goff + 150 * idx + 132);

            auto g_yzz_yzzz = primBuffer.data(goff + 150 * idx + 133);

            auto g_yzz_zzzz = primBuffer.data(goff + 150 * idx + 134);

            auto g_zzz_xxxx = primBuffer.data(goff + 150 * idx + 135);

            auto g_zzz_xxxy = primBuffer.data(goff + 150 * idx + 136);

            auto g_zzz_xxxz = primBuffer.data(goff + 150 * idx + 137);

            auto g_zzz_xxyy = primBuffer.data(goff + 150 * idx + 138);

            auto g_zzz_xxyz = primBuffer.data(goff + 150 * idx + 139);

            auto g_zzz_xxzz = primBuffer.data(goff + 150 * idx + 140);

            auto g_zzz_xyyy = primBuffer.data(goff + 150 * idx + 141);

            auto g_zzz_xyyz = primBuffer.data(goff + 150 * idx + 142);

            auto g_zzz_xyzz = primBuffer.data(goff + 150 * idx + 143);

            auto g_zzz_xzzz = primBuffer.data(goff + 150 * idx + 144);

            auto g_zzz_yyyy = primBuffer.data(goff + 150 * idx + 145);

            auto g_zzz_yyyz = primBuffer.data(goff + 150 * idx + 146);

            auto g_zzz_yyzz = primBuffer.data(goff + 150 * idx + 147);

            auto g_zzz_yzzz = primBuffer.data(goff + 150 * idx + 148);

            auto g_zzz_zzzz = primBuffer.data(goff + 150 * idx + 149);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xx_xxx,\
                                         gk_xx_xxy, gk_xx_xxz, gk_xx_xyy, gk_xx_xyz,\
                                         gk_xx_xzz, gk_xx_yyy, gk_xx_yyz, gk_xx_yzz,\
                                         gk_xx_zzz, gk_xy_xxx, gk_xy_xxy, gk_xy_xxz,\
                                         gk_xy_xyy, gk_xy_xyz, gk_xy_xzz, gk_xy_yyy,\
                                         gk_xy_yyz, gk_xy_yzz, gk_xy_zzz, gk_xz_xxx,\
                                         gk_xz_xxy, gk_xz_xxz, gk_xz_xyy, gk_xz_xyz,\
                                         gk_xz_xzz, gk_xz_yyy, gk_xz_yyz, gk_xz_yzz,\
                                         gk_xz_zzz, gk_yy_xxx, gk_yy_xxy, gk_yy_xxz,\
                                         gk_yy_xyy, gk_yy_xyz, gk_yy_xzz, gk_yy_yyy,\
                                         gk_yy_yyz, gk_yy_yzz, gk_yy_zzz, gk_yz_xxx,\
                                         gk_yz_xxy, gk_yz_xxz, gk_yz_xyy, gk_yz_xyz,\
                                         gk_yz_xzz, gk_yz_yyy, gk_yz_yyz, gk_yz_yzz,\
                                         gk_yz_zzz, gk_zz_xxx, gk_zz_xxy, gk_zz_xxz,\
                                         gk_zz_xyy, gk_zz_xyz, gk_zz_xzz, gk_zz_yyy,\
                                         gk_zz_yyz, gk_zz_yzz, gk_zz_zzz, g20_x_xxxx,\
                                         g20_x_xxxy, g20_x_xxxz, g20_x_xxyy, g20_x_xxyz,\
                                         g20_x_xxzz, g20_x_xyyy, g20_x_xyyz, g20_x_xyzz,\
                                         g20_x_xzzz, g20_x_yyyy, g20_x_yyyz, g20_x_yyzz,\
                                         g20_x_yzzz, g20_x_zzzz, g20_y_xxxx, g20_y_xxxy,\
                                         g20_y_xxxz, g20_y_xxyy, g20_y_xxyz, g20_y_xxzz,\
                                         g20_y_xyyy, g20_y_xyyz, g20_y_xyzz, g20_y_xzzz,\
                                         g20_y_yyyy, g20_y_yyyz, g20_y_yyzz, g20_y_yzzz,\
                                         g20_y_zzzz, g20_z_xxxx, g20_z_xxxy, g20_z_xxxz,\
                                         g20_z_xxyy, g20_z_xxyz, g20_z_xxzz, g20_z_xyyy,\
                                         g20_z_xyyz, g20_z_xyzz, g20_z_xzzz, g20_z_yyyy,\
                                         g20_z_yyyz, g20_z_yyzz, g20_z_yzzz, g20_z_zzzz,\
                                         g21_x_xxxx, g21_x_xxxy, g21_x_xxxz, g21_x_xxyy,\
                                         g21_x_xxyz, g21_x_xxzz, g21_x_xyyy, g21_x_xyyz,\
                                         g21_x_xyzz, g21_x_xzzz, g21_x_yyyy, g21_x_yyyz,\
                                         g21_x_yyzz, g21_x_yzzz, g21_x_zzzz, g21_y_xxxx,\
                                         g21_y_xxxy, g21_y_xxxz, g21_y_xxyy, g21_y_xxyz,\
                                         g21_y_xxzz, g21_y_xyyy, g21_y_xyyz, g21_y_xyzz,\
                                         g21_y_xzzz, g21_y_yyyy, g21_y_yyyz, g21_y_yyzz,\
                                         g21_y_yzzz, g21_y_zzzz, g21_z_xxxx, g21_z_xxxy,\
                                         g21_z_xxxz, g21_z_xxyy, g21_z_xxyz, g21_z_xxzz,\
                                         g21_z_xyyy, g21_z_xyyz, g21_z_xyzz, g21_z_xzzz,\
                                         g21_z_yyyy, g21_z_yyyz, g21_z_yyzz, g21_z_yzzz,\
                                         g21_z_zzzz, g1_xx_xxxx, g1_xx_xxxy, g1_xx_xxxz,\
                                         g1_xx_xxyy, g1_xx_xxyz, g1_xx_xxzz, g1_xx_xyyy,\
                                         g1_xx_xyyz, g1_xx_xyzz, g1_xx_xzzz, g1_xx_yyyy,\
                                         g1_xx_yyyz, g1_xx_yyzz, g1_xx_yzzz, g1_xx_zzzz,\
                                         g1_xy_xxxx, g1_xy_xxxy, g1_xy_xxxz, g1_xy_xxyy,\
                                         g1_xy_xxyz, g1_xy_xxzz, g1_xy_xyyy, g1_xy_xyyz,\
                                         g1_xy_xyzz, g1_xy_xzzz, g1_xy_yyyy, g1_xy_yyyz,\
                                         g1_xy_yyzz, g1_xy_yzzz, g1_xy_zzzz, g1_xz_xxxx,\
                                         g1_xz_xxxy, g1_xz_xxxz, g1_xz_xxyy, g1_xz_xxyz,\
                                         g1_xz_xxzz, g1_xz_xyyy, g1_xz_xyyz, g1_xz_xyzz,\
                                         g1_xz_xzzz, g1_xz_yyyy, g1_xz_yyyz, g1_xz_yyzz,\
                                         g1_xz_yzzz, g1_xz_zzzz, g1_yy_xxxx, g1_yy_xxxy,\
                                         g1_yy_xxxz, g1_yy_xxyy, g1_yy_xxyz, g1_yy_xxzz,\
                                         g1_yy_xyyy, g1_yy_xyyz, g1_yy_xyzz, g1_yy_xzzz,\
                                         g1_yy_yyyy, g1_yy_yyyz, g1_yy_yyzz, g1_yy_yzzz,\
                                         g1_yy_zzzz, g1_yz_xxxx, g1_yz_xxxy, g1_yz_xxxz,\
                                         g1_yz_xxyy, g1_yz_xxyz, g1_yz_xxzz, g1_yz_xyyy,\
                                         g1_yz_xyyz, g1_yz_xyzz, g1_yz_xzzz, g1_yz_yyyy,\
                                         g1_yz_yyyz, g1_yz_yyzz, g1_yz_yzzz, g1_yz_zzzz,\
                                         g1_zz_xxxx, g1_zz_xxxy, g1_zz_xxxz, g1_zz_xxyy,\
                                         g1_zz_xxyz, g1_zz_xxzz, g1_zz_xyyy, g1_zz_xyyz,\
                                         g1_zz_xyzz, g1_zz_xzzz, g1_zz_yyyy, g1_zz_yyyz,\
                                         g1_zz_yyzz, g1_zz_yzzz, g1_zz_zzzz, g_xxx_xxxx,\
                                         g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz,\
                                         g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz,\
                                         g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz,\
                                         g_xxx_yzzz, g_xxx_zzzz, g_xxy_xxxx, g_xxy_xxxy,\
                                         g_xxy_xxxz, g_xxy_xxyy, g_xxy_xxyz, g_xxy_xxzz,\
                                         g_xxy_xyyy, g_xxy_xyyz, g_xxy_xyzz, g_xxy_xzzz,\
                                         g_xxy_yyyy, g_xxy_yyyz, g_xxy_yyzz, g_xxy_yzzz,\
                                         g_xxy_zzzz, g_xxz_xxxx, g_xxz_xxxy, g_xxz_xxxz,\
                                         g_xxz_xxyy, g_xxz_xxyz, g_xxz_xxzz, g_xxz_xyyy,\
                                         g_xxz_xyyz, g_xxz_xyzz, g_xxz_xzzz, g_xxz_yyyy,\
                                         g_xxz_yyyz, g_xxz_yyzz, g_xxz_yzzz, g_xxz_zzzz,\
                                         g_xyy_xxxx, g_xyy_xxxy, g_xyy_xxxz, g_xyy_xxyy,\
                                         g_xyy_xxyz, g_xyy_xxzz, g_xyy_xyyy, g_xyy_xyyz,\
                                         g_xyy_xyzz, g_xyy_xzzz, g_xyy_yyyy, g_xyy_yyyz,\
                                         g_xyy_yyzz, g_xyy_yzzz, g_xyy_zzzz, g_xyz_xxxx,\
                                         g_xyz_xxxy, g_xyz_xxxz, g_xyz_xxyy, g_xyz_xxyz,\
                                         g_xyz_xxzz, g_xyz_xyyy, g_xyz_xyyz, g_xyz_xyzz,\
                                         g_xyz_xzzz, g_xyz_yyyy, g_xyz_yyyz, g_xyz_yyzz,\
                                         g_xyz_yzzz, g_xyz_zzzz, g_xzz_xxxx, g_xzz_xxxy,\
                                         g_xzz_xxxz, g_xzz_xxyy, g_xzz_xxyz, g_xzz_xxzz,\
                                         g_xzz_xyyy, g_xzz_xyyz, g_xzz_xyzz, g_xzz_xzzz,\
                                         g_xzz_yyyy, g_xzz_yyyz, g_xzz_yyzz, g_xzz_yzzz,\
                                         g_xzz_zzzz, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz,\
                                         g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy,\
                                         g_yyy_xyyz, g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy,\
                                         g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz,\
                                         g_yyz_xxxx, g_yyz_xxxy, g_yyz_xxxz, g_yyz_xxyy,\
                                         g_yyz_xxyz, g_yyz_xxzz, g_yyz_xyyy, g_yyz_xyyz,\
                                         g_yyz_xyzz, g_yyz_xzzz, g_yyz_yyyy, g_yyz_yyyz,\
                                         g_yyz_yyzz, g_yyz_yzzz, g_yyz_zzzz, g_yzz_xxxx,\
                                         g_yzz_xxxy, g_yzz_xxxz, g_yzz_xxyy, g_yzz_xxyz,\
                                         g_yzz_xxzz, g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyzz,\
                                         g_yzz_xzzz, g_yzz_yyyy, g_yzz_yyyz, g_yzz_yyzz,\
                                         g_yzz_yzzz, g_yzz_zzzz, g_zzz_xxxx, g_zzz_xxxy,\
                                         g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz,\
                                         g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz,\
                                         g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz,\
                                         g_zzz_zzzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxx_xxxx[k] = fr * g1_xx_xxxx[k] + f2g * (2.0 * g20_x_xxxx[k] - 2.0 * fgz * g21_x_xxxx[k]) + 4.0 * f2t * gk_xx_xxx[k];

                g_xxx_xxxy[k] = fr * g1_xx_xxxy[k] + f2g * (2.0 * g20_x_xxxy[k] - 2.0 * fgz * g21_x_xxxy[k]) + 3.0 * f2t * gk_xx_xxy[k];

                g_xxx_xxxz[k] = fr * g1_xx_xxxz[k] + f2g * (2.0 * g20_x_xxxz[k] - 2.0 * fgz * g21_x_xxxz[k]) + 3.0 * f2t * gk_xx_xxz[k];

                g_xxx_xxyy[k] = fr * g1_xx_xxyy[k] + f2g * (2.0 * g20_x_xxyy[k] - 2.0 * fgz * g21_x_xxyy[k]) + 2.0 * f2t * gk_xx_xyy[k];

                g_xxx_xxyz[k] = fr * g1_xx_xxyz[k] + f2g * (2.0 * g20_x_xxyz[k] - 2.0 * fgz * g21_x_xxyz[k]) + 2.0 * f2t * gk_xx_xyz[k];

                g_xxx_xxzz[k] = fr * g1_xx_xxzz[k] + f2g * (2.0 * g20_x_xxzz[k] - 2.0 * fgz * g21_x_xxzz[k]) + 2.0 * f2t * gk_xx_xzz[k];

                g_xxx_xyyy[k] = fr * g1_xx_xyyy[k] + f2g * (2.0 * g20_x_xyyy[k] - 2.0 * fgz * g21_x_xyyy[k]) + f2t * gk_xx_yyy[k];

                g_xxx_xyyz[k] = fr * g1_xx_xyyz[k] + f2g * (2.0 * g20_x_xyyz[k] - 2.0 * fgz * g21_x_xyyz[k]) + f2t * gk_xx_yyz[k];

                g_xxx_xyzz[k] = fr * g1_xx_xyzz[k] + f2g * (2.0 * g20_x_xyzz[k] - 2.0 * fgz * g21_x_xyzz[k]) + f2t * gk_xx_yzz[k];

                g_xxx_xzzz[k] = fr * g1_xx_xzzz[k] + f2g * (2.0 * g20_x_xzzz[k] - 2.0 * fgz * g21_x_xzzz[k]) + f2t * gk_xx_zzz[k];

                g_xxx_yyyy[k] = fr * g1_xx_yyyy[k] + f2g * (2.0 * g20_x_yyyy[k] - 2.0 * fgz * g21_x_yyyy[k]);

                g_xxx_yyyz[k] = fr * g1_xx_yyyz[k] + f2g * (2.0 * g20_x_yyyz[k] - 2.0 * fgz * g21_x_yyyz[k]);

                g_xxx_yyzz[k] = fr * g1_xx_yyzz[k] + f2g * (2.0 * g20_x_yyzz[k] - 2.0 * fgz * g21_x_yyzz[k]);

                g_xxx_yzzz[k] = fr * g1_xx_yzzz[k] + f2g * (2.0 * g20_x_yzzz[k] - 2.0 * fgz * g21_x_yzzz[k]);

                g_xxx_zzzz[k] = fr * g1_xx_zzzz[k] + f2g * (2.0 * g20_x_zzzz[k] - 2.0 * fgz * g21_x_zzzz[k]);

                g_xxy_xxxx[k] = fr * g1_xy_xxxx[k] + f2g * (g20_y_xxxx[k] - fgz * g21_y_xxxx[k]) + 4.0 * f2t * gk_xy_xxx[k];

                g_xxy_xxxy[k] = fr * g1_xy_xxxy[k] + f2g * (g20_y_xxxy[k] - fgz * g21_y_xxxy[k]) + 3.0 * f2t * gk_xy_xxy[k];

                g_xxy_xxxz[k] = fr * g1_xy_xxxz[k] + f2g * (g20_y_xxxz[k] - fgz * g21_y_xxxz[k]) + 3.0 * f2t * gk_xy_xxz[k];

                g_xxy_xxyy[k] = fr * g1_xy_xxyy[k] + f2g * (g20_y_xxyy[k] - fgz * g21_y_xxyy[k]) + 2.0 * f2t * gk_xy_xyy[k];

                g_xxy_xxyz[k] = fr * g1_xy_xxyz[k] + f2g * (g20_y_xxyz[k] - fgz * g21_y_xxyz[k]) + 2.0 * f2t * gk_xy_xyz[k];

                g_xxy_xxzz[k] = fr * g1_xy_xxzz[k] + f2g * (g20_y_xxzz[k] - fgz * g21_y_xxzz[k]) + 2.0 * f2t * gk_xy_xzz[k];

                g_xxy_xyyy[k] = fr * g1_xy_xyyy[k] + f2g * (g20_y_xyyy[k] - fgz * g21_y_xyyy[k]) + f2t * gk_xy_yyy[k];

                g_xxy_xyyz[k] = fr * g1_xy_xyyz[k] + f2g * (g20_y_xyyz[k] - fgz * g21_y_xyyz[k]) + f2t * gk_xy_yyz[k];

                g_xxy_xyzz[k] = fr * g1_xy_xyzz[k] + f2g * (g20_y_xyzz[k] - fgz * g21_y_xyzz[k]) + f2t * gk_xy_yzz[k];

                g_xxy_xzzz[k] = fr * g1_xy_xzzz[k] + f2g * (g20_y_xzzz[k] - fgz * g21_y_xzzz[k]) + f2t * gk_xy_zzz[k];

                g_xxy_yyyy[k] = fr * g1_xy_yyyy[k] + f2g * (g20_y_yyyy[k] - fgz * g21_y_yyyy[k]);

                g_xxy_yyyz[k] = fr * g1_xy_yyyz[k] + f2g * (g20_y_yyyz[k] - fgz * g21_y_yyyz[k]);

                g_xxy_yyzz[k] = fr * g1_xy_yyzz[k] + f2g * (g20_y_yyzz[k] - fgz * g21_y_yyzz[k]);

                g_xxy_yzzz[k] = fr * g1_xy_yzzz[k] + f2g * (g20_y_yzzz[k] - fgz * g21_y_yzzz[k]);

                g_xxy_zzzz[k] = fr * g1_xy_zzzz[k] + f2g * (g20_y_zzzz[k] - fgz * g21_y_zzzz[k]);

                g_xxz_xxxx[k] = fr * g1_xz_xxxx[k] + f2g * (g20_z_xxxx[k] - fgz * g21_z_xxxx[k]) + 4.0 * f2t * gk_xz_xxx[k];

                g_xxz_xxxy[k] = fr * g1_xz_xxxy[k] + f2g * (g20_z_xxxy[k] - fgz * g21_z_xxxy[k]) + 3.0 * f2t * gk_xz_xxy[k];

                g_xxz_xxxz[k] = fr * g1_xz_xxxz[k] + f2g * (g20_z_xxxz[k] - fgz * g21_z_xxxz[k]) + 3.0 * f2t * gk_xz_xxz[k];

                g_xxz_xxyy[k] = fr * g1_xz_xxyy[k] + f2g * (g20_z_xxyy[k] - fgz * g21_z_xxyy[k]) + 2.0 * f2t * gk_xz_xyy[k];

                g_xxz_xxyz[k] = fr * g1_xz_xxyz[k] + f2g * (g20_z_xxyz[k] - fgz * g21_z_xxyz[k]) + 2.0 * f2t * gk_xz_xyz[k];

                g_xxz_xxzz[k] = fr * g1_xz_xxzz[k] + f2g * (g20_z_xxzz[k] - fgz * g21_z_xxzz[k]) + 2.0 * f2t * gk_xz_xzz[k];

                g_xxz_xyyy[k] = fr * g1_xz_xyyy[k] + f2g * (g20_z_xyyy[k] - fgz * g21_z_xyyy[k]) + f2t * gk_xz_yyy[k];

                g_xxz_xyyz[k] = fr * g1_xz_xyyz[k] + f2g * (g20_z_xyyz[k] - fgz * g21_z_xyyz[k]) + f2t * gk_xz_yyz[k];

                g_xxz_xyzz[k] = fr * g1_xz_xyzz[k] + f2g * (g20_z_xyzz[k] - fgz * g21_z_xyzz[k]) + f2t * gk_xz_yzz[k];

                g_xxz_xzzz[k] = fr * g1_xz_xzzz[k] + f2g * (g20_z_xzzz[k] - fgz * g21_z_xzzz[k]) + f2t * gk_xz_zzz[k];

                g_xxz_yyyy[k] = fr * g1_xz_yyyy[k] + f2g * (g20_z_yyyy[k] - fgz * g21_z_yyyy[k]);

                g_xxz_yyyz[k] = fr * g1_xz_yyyz[k] + f2g * (g20_z_yyyz[k] - fgz * g21_z_yyyz[k]);

                g_xxz_yyzz[k] = fr * g1_xz_yyzz[k] + f2g * (g20_z_yyzz[k] - fgz * g21_z_yyzz[k]);

                g_xxz_yzzz[k] = fr * g1_xz_yzzz[k] + f2g * (g20_z_yzzz[k] - fgz * g21_z_yzzz[k]);

                g_xxz_zzzz[k] = fr * g1_xz_zzzz[k] + f2g * (g20_z_zzzz[k] - fgz * g21_z_zzzz[k]);

                g_xyy_xxxx[k] = fr * g1_yy_xxxx[k] + 4.0 * f2t * gk_yy_xxx[k];

                g_xyy_xxxy[k] = fr * g1_yy_xxxy[k] + 3.0 * f2t * gk_yy_xxy[k];

                g_xyy_xxxz[k] = fr * g1_yy_xxxz[k] + 3.0 * f2t * gk_yy_xxz[k];

                g_xyy_xxyy[k] = fr * g1_yy_xxyy[k] + 2.0 * f2t * gk_yy_xyy[k];

                g_xyy_xxyz[k] = fr * g1_yy_xxyz[k] + 2.0 * f2t * gk_yy_xyz[k];

                g_xyy_xxzz[k] = fr * g1_yy_xxzz[k] + 2.0 * f2t * gk_yy_xzz[k];

                g_xyy_xyyy[k] = fr * g1_yy_xyyy[k] + f2t * gk_yy_yyy[k];

                g_xyy_xyyz[k] = fr * g1_yy_xyyz[k] + f2t * gk_yy_yyz[k];

                g_xyy_xyzz[k] = fr * g1_yy_xyzz[k] + f2t * gk_yy_yzz[k];

                g_xyy_xzzz[k] = fr * g1_yy_xzzz[k] + f2t * gk_yy_zzz[k];

                g_xyy_yyyy[k] = fr * g1_yy_yyyy[k];

                g_xyy_yyyz[k] = fr * g1_yy_yyyz[k];

                g_xyy_yyzz[k] = fr * g1_yy_yyzz[k];

                g_xyy_yzzz[k] = fr * g1_yy_yzzz[k];

                g_xyy_zzzz[k] = fr * g1_yy_zzzz[k];

                g_xyz_xxxx[k] = fr * g1_yz_xxxx[k] + 4.0 * f2t * gk_yz_xxx[k];

                g_xyz_xxxy[k] = fr * g1_yz_xxxy[k] + 3.0 * f2t * gk_yz_xxy[k];

                g_xyz_xxxz[k] = fr * g1_yz_xxxz[k] + 3.0 * f2t * gk_yz_xxz[k];

                g_xyz_xxyy[k] = fr * g1_yz_xxyy[k] + 2.0 * f2t * gk_yz_xyy[k];

                g_xyz_xxyz[k] = fr * g1_yz_xxyz[k] + 2.0 * f2t * gk_yz_xyz[k];

                g_xyz_xxzz[k] = fr * g1_yz_xxzz[k] + 2.0 * f2t * gk_yz_xzz[k];

                g_xyz_xyyy[k] = fr * g1_yz_xyyy[k] + f2t * gk_yz_yyy[k];

                g_xyz_xyyz[k] = fr * g1_yz_xyyz[k] + f2t * gk_yz_yyz[k];

                g_xyz_xyzz[k] = fr * g1_yz_xyzz[k] + f2t * gk_yz_yzz[k];

                g_xyz_xzzz[k] = fr * g1_yz_xzzz[k] + f2t * gk_yz_zzz[k];

                g_xyz_yyyy[k] = fr * g1_yz_yyyy[k];

                g_xyz_yyyz[k] = fr * g1_yz_yyyz[k];

                g_xyz_yyzz[k] = fr * g1_yz_yyzz[k];

                g_xyz_yzzz[k] = fr * g1_yz_yzzz[k];

                g_xyz_zzzz[k] = fr * g1_yz_zzzz[k];

                g_xzz_xxxx[k] = fr * g1_zz_xxxx[k] + 4.0 * f2t * gk_zz_xxx[k];

                g_xzz_xxxy[k] = fr * g1_zz_xxxy[k] + 3.0 * f2t * gk_zz_xxy[k];

                g_xzz_xxxz[k] = fr * g1_zz_xxxz[k] + 3.0 * f2t * gk_zz_xxz[k];

                g_xzz_xxyy[k] = fr * g1_zz_xxyy[k] + 2.0 * f2t * gk_zz_xyy[k];

                g_xzz_xxyz[k] = fr * g1_zz_xxyz[k] + 2.0 * f2t * gk_zz_xyz[k];

                g_xzz_xxzz[k] = fr * g1_zz_xxzz[k] + 2.0 * f2t * gk_zz_xzz[k];

                g_xzz_xyyy[k] = fr * g1_zz_xyyy[k] + f2t * gk_zz_yyy[k];

                g_xzz_xyyz[k] = fr * g1_zz_xyyz[k] + f2t * gk_zz_yyz[k];

                g_xzz_xyzz[k] = fr * g1_zz_xyzz[k] + f2t * gk_zz_yzz[k];

                g_xzz_xzzz[k] = fr * g1_zz_xzzz[k] + f2t * gk_zz_zzz[k];

                g_xzz_yyyy[k] = fr * g1_zz_yyyy[k];

                g_xzz_yyyz[k] = fr * g1_zz_yyyz[k];

                g_xzz_yyzz[k] = fr * g1_zz_yyzz[k];

                g_xzz_yzzz[k] = fr * g1_zz_yzzz[k];

                g_xzz_zzzz[k] = fr * g1_zz_zzzz[k];

                // leading y component

                fr = pay[k];

                g_yyy_xxxx[k] = fr * g1_yy_xxxx[k] + f2g * (2.0 * g20_y_xxxx[k] - 2.0 * fgz * g21_y_xxxx[k]);

                g_yyy_xxxy[k] = fr * g1_yy_xxxy[k] + f2g * (2.0 * g20_y_xxxy[k] - 2.0 * fgz * g21_y_xxxy[k]) + f2t * gk_yy_xxx[k];

                g_yyy_xxxz[k] = fr * g1_yy_xxxz[k] + f2g * (2.0 * g20_y_xxxz[k] - 2.0 * fgz * g21_y_xxxz[k]);

                g_yyy_xxyy[k] = fr * g1_yy_xxyy[k] + f2g * (2.0 * g20_y_xxyy[k] - 2.0 * fgz * g21_y_xxyy[k]) + 2.0 * f2t * gk_yy_xxy[k];

                g_yyy_xxyz[k] = fr * g1_yy_xxyz[k] + f2g * (2.0 * g20_y_xxyz[k] - 2.0 * fgz * g21_y_xxyz[k]) + f2t * gk_yy_xxz[k];

                g_yyy_xxzz[k] = fr * g1_yy_xxzz[k] + f2g * (2.0 * g20_y_xxzz[k] - 2.0 * fgz * g21_y_xxzz[k]);

                g_yyy_xyyy[k] = fr * g1_yy_xyyy[k] + f2g * (2.0 * g20_y_xyyy[k] - 2.0 * fgz * g21_y_xyyy[k]) + 3.0 * f2t * gk_yy_xyy[k];

                g_yyy_xyyz[k] = fr * g1_yy_xyyz[k] + f2g * (2.0 * g20_y_xyyz[k] - 2.0 * fgz * g21_y_xyyz[k]) + 2.0 * f2t * gk_yy_xyz[k];

                g_yyy_xyzz[k] = fr * g1_yy_xyzz[k] + f2g * (2.0 * g20_y_xyzz[k] - 2.0 * fgz * g21_y_xyzz[k]) + f2t * gk_yy_xzz[k];

                g_yyy_xzzz[k] = fr * g1_yy_xzzz[k] + f2g * (2.0 * g20_y_xzzz[k] - 2.0 * fgz * g21_y_xzzz[k]);

                g_yyy_yyyy[k] = fr * g1_yy_yyyy[k] + f2g * (2.0 * g20_y_yyyy[k] - 2.0 * fgz * g21_y_yyyy[k]) + 4.0 * f2t * gk_yy_yyy[k];

                g_yyy_yyyz[k] = fr * g1_yy_yyyz[k] + f2g * (2.0 * g20_y_yyyz[k] - 2.0 * fgz * g21_y_yyyz[k]) + 3.0 * f2t * gk_yy_yyz[k];

                g_yyy_yyzz[k] = fr * g1_yy_yyzz[k] + f2g * (2.0 * g20_y_yyzz[k] - 2.0 * fgz * g21_y_yyzz[k]) + 2.0 * f2t * gk_yy_yzz[k];

                g_yyy_yzzz[k] = fr * g1_yy_yzzz[k] + f2g * (2.0 * g20_y_yzzz[k] - 2.0 * fgz * g21_y_yzzz[k]) + f2t * gk_yy_zzz[k];

                g_yyy_zzzz[k] = fr * g1_yy_zzzz[k] + f2g * (2.0 * g20_y_zzzz[k] - 2.0 * fgz * g21_y_zzzz[k]);

                g_yyz_xxxx[k] = fr * g1_yz_xxxx[k] + f2g * (g20_z_xxxx[k] - fgz * g21_z_xxxx[k]);

                g_yyz_xxxy[k] = fr * g1_yz_xxxy[k] + f2g * (g20_z_xxxy[k] - fgz * g21_z_xxxy[k]) + f2t * gk_yz_xxx[k];

                g_yyz_xxxz[k] = fr * g1_yz_xxxz[k] + f2g * (g20_z_xxxz[k] - fgz * g21_z_xxxz[k]);

                g_yyz_xxyy[k] = fr * g1_yz_xxyy[k] + f2g * (g20_z_xxyy[k] - fgz * g21_z_xxyy[k]) + 2.0 * f2t * gk_yz_xxy[k];

                g_yyz_xxyz[k] = fr * g1_yz_xxyz[k] + f2g * (g20_z_xxyz[k] - fgz * g21_z_xxyz[k]) + f2t * gk_yz_xxz[k];

                g_yyz_xxzz[k] = fr * g1_yz_xxzz[k] + f2g * (g20_z_xxzz[k] - fgz * g21_z_xxzz[k]);

                g_yyz_xyyy[k] = fr * g1_yz_xyyy[k] + f2g * (g20_z_xyyy[k] - fgz * g21_z_xyyy[k]) + 3.0 * f2t * gk_yz_xyy[k];

                g_yyz_xyyz[k] = fr * g1_yz_xyyz[k] + f2g * (g20_z_xyyz[k] - fgz * g21_z_xyyz[k]) + 2.0 * f2t * gk_yz_xyz[k];

                g_yyz_xyzz[k] = fr * g1_yz_xyzz[k] + f2g * (g20_z_xyzz[k] - fgz * g21_z_xyzz[k]) + f2t * gk_yz_xzz[k];

                g_yyz_xzzz[k] = fr * g1_yz_xzzz[k] + f2g * (g20_z_xzzz[k] - fgz * g21_z_xzzz[k]);

                g_yyz_yyyy[k] = fr * g1_yz_yyyy[k] + f2g * (g20_z_yyyy[k] - fgz * g21_z_yyyy[k]) + 4.0 * f2t * gk_yz_yyy[k];

                g_yyz_yyyz[k] = fr * g1_yz_yyyz[k] + f2g * (g20_z_yyyz[k] - fgz * g21_z_yyyz[k]) + 3.0 * f2t * gk_yz_yyz[k];

                g_yyz_yyzz[k] = fr * g1_yz_yyzz[k] + f2g * (g20_z_yyzz[k] - fgz * g21_z_yyzz[k]) + 2.0 * f2t * gk_yz_yzz[k];

                g_yyz_yzzz[k] = fr * g1_yz_yzzz[k] + f2g * (g20_z_yzzz[k] - fgz * g21_z_yzzz[k]) + f2t * gk_yz_zzz[k];

                g_yyz_zzzz[k] = fr * g1_yz_zzzz[k] + f2g * (g20_z_zzzz[k] - fgz * g21_z_zzzz[k]);

                g_yzz_xxxx[k] = fr * g1_zz_xxxx[k];

                g_yzz_xxxy[k] = fr * g1_zz_xxxy[k] + f2t * gk_zz_xxx[k];

                g_yzz_xxxz[k] = fr * g1_zz_xxxz[k];

                g_yzz_xxyy[k] = fr * g1_zz_xxyy[k] + 2.0 * f2t * gk_zz_xxy[k];

                g_yzz_xxyz[k] = fr * g1_zz_xxyz[k] + f2t * gk_zz_xxz[k];

                g_yzz_xxzz[k] = fr * g1_zz_xxzz[k];

                g_yzz_xyyy[k] = fr * g1_zz_xyyy[k] + 3.0 * f2t * gk_zz_xyy[k];

                g_yzz_xyyz[k] = fr * g1_zz_xyyz[k] + 2.0 * f2t * gk_zz_xyz[k];

                g_yzz_xyzz[k] = fr * g1_zz_xyzz[k] + f2t * gk_zz_xzz[k];

                g_yzz_xzzz[k] = fr * g1_zz_xzzz[k];

                g_yzz_yyyy[k] = fr * g1_zz_yyyy[k] + 4.0 * f2t * gk_zz_yyy[k];

                g_yzz_yyyz[k] = fr * g1_zz_yyyz[k] + 3.0 * f2t * gk_zz_yyz[k];

                g_yzz_yyzz[k] = fr * g1_zz_yyzz[k] + 2.0 * f2t * gk_zz_yzz[k];

                g_yzz_yzzz[k] = fr * g1_zz_yzzz[k] + f2t * gk_zz_zzz[k];

                g_yzz_zzzz[k] = fr * g1_zz_zzzz[k];

                // leading z component

                fr = paz[k];

                g_zzz_xxxx[k] = fr * g1_zz_xxxx[k] + f2g * (2.0 * g20_z_xxxx[k] - 2.0 * fgz * g21_z_xxxx[k]);

                g_zzz_xxxy[k] = fr * g1_zz_xxxy[k] + f2g * (2.0 * g20_z_xxxy[k] - 2.0 * fgz * g21_z_xxxy[k]);

                g_zzz_xxxz[k] = fr * g1_zz_xxxz[k] + f2g * (2.0 * g20_z_xxxz[k] - 2.0 * fgz * g21_z_xxxz[k]) + f2t * gk_zz_xxx[k];

                g_zzz_xxyy[k] = fr * g1_zz_xxyy[k] + f2g * (2.0 * g20_z_xxyy[k] - 2.0 * fgz * g21_z_xxyy[k]);

                g_zzz_xxyz[k] = fr * g1_zz_xxyz[k] + f2g * (2.0 * g20_z_xxyz[k] - 2.0 * fgz * g21_z_xxyz[k]) + f2t * gk_zz_xxy[k];

                g_zzz_xxzz[k] = fr * g1_zz_xxzz[k] + f2g * (2.0 * g20_z_xxzz[k] - 2.0 * fgz * g21_z_xxzz[k]) + 2.0 * f2t * gk_zz_xxz[k];

                g_zzz_xyyy[k] = fr * g1_zz_xyyy[k] + f2g * (2.0 * g20_z_xyyy[k] - 2.0 * fgz * g21_z_xyyy[k]);

                g_zzz_xyyz[k] = fr * g1_zz_xyyz[k] + f2g * (2.0 * g20_z_xyyz[k] - 2.0 * fgz * g21_z_xyyz[k]) + f2t * gk_zz_xyy[k];

                g_zzz_xyzz[k] = fr * g1_zz_xyzz[k] + f2g * (2.0 * g20_z_xyzz[k] - 2.0 * fgz * g21_z_xyzz[k]) + 2.0 * f2t * gk_zz_xyz[k];

                g_zzz_xzzz[k] = fr * g1_zz_xzzz[k] + f2g * (2.0 * g20_z_xzzz[k] - 2.0 * fgz * g21_z_xzzz[k]) + 3.0 * f2t * gk_zz_xzz[k];

                g_zzz_yyyy[k] = fr * g1_zz_yyyy[k] + f2g * (2.0 * g20_z_yyyy[k] - 2.0 * fgz * g21_z_yyyy[k]);

                g_zzz_yyyz[k] = fr * g1_zz_yyyz[k] + f2g * (2.0 * g20_z_yyyz[k] - 2.0 * fgz * g21_z_yyyz[k]) + f2t * gk_zz_yyy[k];

                g_zzz_yyzz[k] = fr * g1_zz_yyzz[k] + f2g * (2.0 * g20_z_yyzz[k] - 2.0 * fgz * g21_z_yyzz[k]) + 2.0 * f2t * gk_zz_yyz[k];

                g_zzz_yzzz[k] = fr * g1_zz_yzzz[k] + f2g * (2.0 * g20_z_yzzz[k] - 2.0 * fgz * g21_z_yzzz[k]) + 3.0 * f2t * gk_zz_yzz[k];

                g_zzz_zzzz[k] = fr * g1_zz_zzzz[k] + f2g * (2.0 * g20_z_zzzz[k] - 2.0 * fgz * g21_z_zzzz[k]) + 4.0 * f2t * gk_zz_zzz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForGF(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 4, 3);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {4, 3, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 3, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 2, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|g(r,r')|D)^(m+1) integrals

            auto gk_xxx_xx = primBuffer.data(gkoff + 60 * idx);

            auto gk_xxx_xy = primBuffer.data(gkoff + 60 * idx + 1);

            auto gk_xxx_xz = primBuffer.data(gkoff + 60 * idx + 2);

            auto gk_xxx_yy = primBuffer.data(gkoff + 60 * idx + 3);

            auto gk_xxx_yz = primBuffer.data(gkoff + 60 * idx + 4);

            auto gk_xxx_zz = primBuffer.data(gkoff + 60 * idx + 5);

            auto gk_xxy_xx = primBuffer.data(gkoff + 60 * idx + 6);

            auto gk_xxy_xy = primBuffer.data(gkoff + 60 * idx + 7);

            auto gk_xxy_xz = primBuffer.data(gkoff + 60 * idx + 8);

            auto gk_xxy_yy = primBuffer.data(gkoff + 60 * idx + 9);

            auto gk_xxy_yz = primBuffer.data(gkoff + 60 * idx + 10);

            auto gk_xxy_zz = primBuffer.data(gkoff + 60 * idx + 11);

            auto gk_xxz_xx = primBuffer.data(gkoff + 60 * idx + 12);

            auto gk_xxz_xy = primBuffer.data(gkoff + 60 * idx + 13);

            auto gk_xxz_xz = primBuffer.data(gkoff + 60 * idx + 14);

            auto gk_xxz_yy = primBuffer.data(gkoff + 60 * idx + 15);

            auto gk_xxz_yz = primBuffer.data(gkoff + 60 * idx + 16);

            auto gk_xxz_zz = primBuffer.data(gkoff + 60 * idx + 17);

            auto gk_xyy_xx = primBuffer.data(gkoff + 60 * idx + 18);

            auto gk_xyy_xy = primBuffer.data(gkoff + 60 * idx + 19);

            auto gk_xyy_xz = primBuffer.data(gkoff + 60 * idx + 20);

            auto gk_xyy_yy = primBuffer.data(gkoff + 60 * idx + 21);

            auto gk_xyy_yz = primBuffer.data(gkoff + 60 * idx + 22);

            auto gk_xyy_zz = primBuffer.data(gkoff + 60 * idx + 23);

            auto gk_xyz_xx = primBuffer.data(gkoff + 60 * idx + 24);

            auto gk_xyz_xy = primBuffer.data(gkoff + 60 * idx + 25);

            auto gk_xyz_xz = primBuffer.data(gkoff + 60 * idx + 26);

            auto gk_xyz_yy = primBuffer.data(gkoff + 60 * idx + 27);

            auto gk_xyz_yz = primBuffer.data(gkoff + 60 * idx + 28);

            auto gk_xyz_zz = primBuffer.data(gkoff + 60 * idx + 29);

            auto gk_xzz_xx = primBuffer.data(gkoff + 60 * idx + 30);

            auto gk_xzz_xy = primBuffer.data(gkoff + 60 * idx + 31);

            auto gk_xzz_xz = primBuffer.data(gkoff + 60 * idx + 32);

            auto gk_xzz_yy = primBuffer.data(gkoff + 60 * idx + 33);

            auto gk_xzz_yz = primBuffer.data(gkoff + 60 * idx + 34);

            auto gk_xzz_zz = primBuffer.data(gkoff + 60 * idx + 35);

            auto gk_yyy_xx = primBuffer.data(gkoff + 60 * idx + 36);

            auto gk_yyy_xy = primBuffer.data(gkoff + 60 * idx + 37);

            auto gk_yyy_xz = primBuffer.data(gkoff + 60 * idx + 38);

            auto gk_yyy_yy = primBuffer.data(gkoff + 60 * idx + 39);

            auto gk_yyy_yz = primBuffer.data(gkoff + 60 * idx + 40);

            auto gk_yyy_zz = primBuffer.data(gkoff + 60 * idx + 41);

            auto gk_yyz_xx = primBuffer.data(gkoff + 60 * idx + 42);

            auto gk_yyz_xy = primBuffer.data(gkoff + 60 * idx + 43);

            auto gk_yyz_xz = primBuffer.data(gkoff + 60 * idx + 44);

            auto gk_yyz_yy = primBuffer.data(gkoff + 60 * idx + 45);

            auto gk_yyz_yz = primBuffer.data(gkoff + 60 * idx + 46);

            auto gk_yyz_zz = primBuffer.data(gkoff + 60 * idx + 47);

            auto gk_yzz_xx = primBuffer.data(gkoff + 60 * idx + 48);

            auto gk_yzz_xy = primBuffer.data(gkoff + 60 * idx + 49);

            auto gk_yzz_xz = primBuffer.data(gkoff + 60 * idx + 50);

            auto gk_yzz_yy = primBuffer.data(gkoff + 60 * idx + 51);

            auto gk_yzz_yz = primBuffer.data(gkoff + 60 * idx + 52);

            auto gk_yzz_zz = primBuffer.data(gkoff + 60 * idx + 53);

            auto gk_zzz_xx = primBuffer.data(gkoff + 60 * idx + 54);

            auto gk_zzz_xy = primBuffer.data(gkoff + 60 * idx + 55);

            auto gk_zzz_xz = primBuffer.data(gkoff + 60 * idx + 56);

            auto gk_zzz_yy = primBuffer.data(gkoff + 60 * idx + 57);

            auto gk_zzz_yz = primBuffer.data(gkoff + 60 * idx + 58);

            auto gk_zzz_zz = primBuffer.data(gkoff + 60 * idx + 59);

            // set up pointers to (D|g(r,r')|F)^(m) integrals

            auto g20_xx_xxx = primBuffer.data(g20off + 60 * idx);

            auto g20_xx_xxy = primBuffer.data(g20off + 60 * idx + 1);

            auto g20_xx_xxz = primBuffer.data(g20off + 60 * idx + 2);

            auto g20_xx_xyy = primBuffer.data(g20off + 60 * idx + 3);

            auto g20_xx_xyz = primBuffer.data(g20off + 60 * idx + 4);

            auto g20_xx_xzz = primBuffer.data(g20off + 60 * idx + 5);

            auto g20_xx_yyy = primBuffer.data(g20off + 60 * idx + 6);

            auto g20_xx_yyz = primBuffer.data(g20off + 60 * idx + 7);

            auto g20_xx_yzz = primBuffer.data(g20off + 60 * idx + 8);

            auto g20_xx_zzz = primBuffer.data(g20off + 60 * idx + 9);

            auto g20_xy_xxx = primBuffer.data(g20off + 60 * idx + 10);

            auto g20_xy_xxy = primBuffer.data(g20off + 60 * idx + 11);

            auto g20_xy_xxz = primBuffer.data(g20off + 60 * idx + 12);

            auto g20_xy_xyy = primBuffer.data(g20off + 60 * idx + 13);

            auto g20_xy_xyz = primBuffer.data(g20off + 60 * idx + 14);

            auto g20_xy_xzz = primBuffer.data(g20off + 60 * idx + 15);

            auto g20_xy_yyy = primBuffer.data(g20off + 60 * idx + 16);

            auto g20_xy_yyz = primBuffer.data(g20off + 60 * idx + 17);

            auto g20_xy_yzz = primBuffer.data(g20off + 60 * idx + 18);

            auto g20_xy_zzz = primBuffer.data(g20off + 60 * idx + 19);

            auto g20_xz_xxx = primBuffer.data(g20off + 60 * idx + 20);

            auto g20_xz_xxy = primBuffer.data(g20off + 60 * idx + 21);

            auto g20_xz_xxz = primBuffer.data(g20off + 60 * idx + 22);

            auto g20_xz_xyy = primBuffer.data(g20off + 60 * idx + 23);

            auto g20_xz_xyz = primBuffer.data(g20off + 60 * idx + 24);

            auto g20_xz_xzz = primBuffer.data(g20off + 60 * idx + 25);

            auto g20_xz_yyy = primBuffer.data(g20off + 60 * idx + 26);

            auto g20_xz_yyz = primBuffer.data(g20off + 60 * idx + 27);

            auto g20_xz_yzz = primBuffer.data(g20off + 60 * idx + 28);

            auto g20_xz_zzz = primBuffer.data(g20off + 60 * idx + 29);

            auto g20_yy_xxx = primBuffer.data(g20off + 60 * idx + 30);

            auto g20_yy_xxy = primBuffer.data(g20off + 60 * idx + 31);

            auto g20_yy_xxz = primBuffer.data(g20off + 60 * idx + 32);

            auto g20_yy_xyy = primBuffer.data(g20off + 60 * idx + 33);

            auto g20_yy_xyz = primBuffer.data(g20off + 60 * idx + 34);

            auto g20_yy_xzz = primBuffer.data(g20off + 60 * idx + 35);

            auto g20_yy_yyy = primBuffer.data(g20off + 60 * idx + 36);

            auto g20_yy_yyz = primBuffer.data(g20off + 60 * idx + 37);

            auto g20_yy_yzz = primBuffer.data(g20off + 60 * idx + 38);

            auto g20_yy_zzz = primBuffer.data(g20off + 60 * idx + 39);

            auto g20_yz_xxx = primBuffer.data(g20off + 60 * idx + 40);

            auto g20_yz_xxy = primBuffer.data(g20off + 60 * idx + 41);

            auto g20_yz_xxz = primBuffer.data(g20off + 60 * idx + 42);

            auto g20_yz_xyy = primBuffer.data(g20off + 60 * idx + 43);

            auto g20_yz_xyz = primBuffer.data(g20off + 60 * idx + 44);

            auto g20_yz_xzz = primBuffer.data(g20off + 60 * idx + 45);

            auto g20_yz_yyy = primBuffer.data(g20off + 60 * idx + 46);

            auto g20_yz_yyz = primBuffer.data(g20off + 60 * idx + 47);

            auto g20_yz_yzz = primBuffer.data(g20off + 60 * idx + 48);

            auto g20_yz_zzz = primBuffer.data(g20off + 60 * idx + 49);

            auto g20_zz_xxx = primBuffer.data(g20off + 60 * idx + 50);

            auto g20_zz_xxy = primBuffer.data(g20off + 60 * idx + 51);

            auto g20_zz_xxz = primBuffer.data(g20off + 60 * idx + 52);

            auto g20_zz_xyy = primBuffer.data(g20off + 60 * idx + 53);

            auto g20_zz_xyz = primBuffer.data(g20off + 60 * idx + 54);

            auto g20_zz_xzz = primBuffer.data(g20off + 60 * idx + 55);

            auto g20_zz_yyy = primBuffer.data(g20off + 60 * idx + 56);

            auto g20_zz_yyz = primBuffer.data(g20off + 60 * idx + 57);

            auto g20_zz_yzz = primBuffer.data(g20off + 60 * idx + 58);

            auto g20_zz_zzz = primBuffer.data(g20off + 60 * idx + 59);

            // set up pointers to (D|g(r,r')|F)^(m+1) integrals

            auto g21_xx_xxx = primBuffer.data(g21off + 60 * idx);

            auto g21_xx_xxy = primBuffer.data(g21off + 60 * idx + 1);

            auto g21_xx_xxz = primBuffer.data(g21off + 60 * idx + 2);

            auto g21_xx_xyy = primBuffer.data(g21off + 60 * idx + 3);

            auto g21_xx_xyz = primBuffer.data(g21off + 60 * idx + 4);

            auto g21_xx_xzz = primBuffer.data(g21off + 60 * idx + 5);

            auto g21_xx_yyy = primBuffer.data(g21off + 60 * idx + 6);

            auto g21_xx_yyz = primBuffer.data(g21off + 60 * idx + 7);

            auto g21_xx_yzz = primBuffer.data(g21off + 60 * idx + 8);

            auto g21_xx_zzz = primBuffer.data(g21off + 60 * idx + 9);

            auto g21_xy_xxx = primBuffer.data(g21off + 60 * idx + 10);

            auto g21_xy_xxy = primBuffer.data(g21off + 60 * idx + 11);

            auto g21_xy_xxz = primBuffer.data(g21off + 60 * idx + 12);

            auto g21_xy_xyy = primBuffer.data(g21off + 60 * idx + 13);

            auto g21_xy_xyz = primBuffer.data(g21off + 60 * idx + 14);

            auto g21_xy_xzz = primBuffer.data(g21off + 60 * idx + 15);

            auto g21_xy_yyy = primBuffer.data(g21off + 60 * idx + 16);

            auto g21_xy_yyz = primBuffer.data(g21off + 60 * idx + 17);

            auto g21_xy_yzz = primBuffer.data(g21off + 60 * idx + 18);

            auto g21_xy_zzz = primBuffer.data(g21off + 60 * idx + 19);

            auto g21_xz_xxx = primBuffer.data(g21off + 60 * idx + 20);

            auto g21_xz_xxy = primBuffer.data(g21off + 60 * idx + 21);

            auto g21_xz_xxz = primBuffer.data(g21off + 60 * idx + 22);

            auto g21_xz_xyy = primBuffer.data(g21off + 60 * idx + 23);

            auto g21_xz_xyz = primBuffer.data(g21off + 60 * idx + 24);

            auto g21_xz_xzz = primBuffer.data(g21off + 60 * idx + 25);

            auto g21_xz_yyy = primBuffer.data(g21off + 60 * idx + 26);

            auto g21_xz_yyz = primBuffer.data(g21off + 60 * idx + 27);

            auto g21_xz_yzz = primBuffer.data(g21off + 60 * idx + 28);

            auto g21_xz_zzz = primBuffer.data(g21off + 60 * idx + 29);

            auto g21_yy_xxx = primBuffer.data(g21off + 60 * idx + 30);

            auto g21_yy_xxy = primBuffer.data(g21off + 60 * idx + 31);

            auto g21_yy_xxz = primBuffer.data(g21off + 60 * idx + 32);

            auto g21_yy_xyy = primBuffer.data(g21off + 60 * idx + 33);

            auto g21_yy_xyz = primBuffer.data(g21off + 60 * idx + 34);

            auto g21_yy_xzz = primBuffer.data(g21off + 60 * idx + 35);

            auto g21_yy_yyy = primBuffer.data(g21off + 60 * idx + 36);

            auto g21_yy_yyz = primBuffer.data(g21off + 60 * idx + 37);

            auto g21_yy_yzz = primBuffer.data(g21off + 60 * idx + 38);

            auto g21_yy_zzz = primBuffer.data(g21off + 60 * idx + 39);

            auto g21_yz_xxx = primBuffer.data(g21off + 60 * idx + 40);

            auto g21_yz_xxy = primBuffer.data(g21off + 60 * idx + 41);

            auto g21_yz_xxz = primBuffer.data(g21off + 60 * idx + 42);

            auto g21_yz_xyy = primBuffer.data(g21off + 60 * idx + 43);

            auto g21_yz_xyz = primBuffer.data(g21off + 60 * idx + 44);

            auto g21_yz_xzz = primBuffer.data(g21off + 60 * idx + 45);

            auto g21_yz_yyy = primBuffer.data(g21off + 60 * idx + 46);

            auto g21_yz_yyz = primBuffer.data(g21off + 60 * idx + 47);

            auto g21_yz_yzz = primBuffer.data(g21off + 60 * idx + 48);

            auto g21_yz_zzz = primBuffer.data(g21off + 60 * idx + 49);

            auto g21_zz_xxx = primBuffer.data(g21off + 60 * idx + 50);

            auto g21_zz_xxy = primBuffer.data(g21off + 60 * idx + 51);

            auto g21_zz_xxz = primBuffer.data(g21off + 60 * idx + 52);

            auto g21_zz_xyy = primBuffer.data(g21off + 60 * idx + 53);

            auto g21_zz_xyz = primBuffer.data(g21off + 60 * idx + 54);

            auto g21_zz_xzz = primBuffer.data(g21off + 60 * idx + 55);

            auto g21_zz_yyy = primBuffer.data(g21off + 60 * idx + 56);

            auto g21_zz_yyz = primBuffer.data(g21off + 60 * idx + 57);

            auto g21_zz_yzz = primBuffer.data(g21off + 60 * idx + 58);

            auto g21_zz_zzz = primBuffer.data(g21off + 60 * idx + 59);

            // set up pointers to (F|g(r,r')|F)^(m+1) integrals

            auto g1_xxx_xxx = primBuffer.data(g1off + 100 * idx);

            auto g1_xxx_xxy = primBuffer.data(g1off + 100 * idx + 1);

            auto g1_xxx_xxz = primBuffer.data(g1off + 100 * idx + 2);

            auto g1_xxx_xyy = primBuffer.data(g1off + 100 * idx + 3);

            auto g1_xxx_xyz = primBuffer.data(g1off + 100 * idx + 4);

            auto g1_xxx_xzz = primBuffer.data(g1off + 100 * idx + 5);

            auto g1_xxx_yyy = primBuffer.data(g1off + 100 * idx + 6);

            auto g1_xxx_yyz = primBuffer.data(g1off + 100 * idx + 7);

            auto g1_xxx_yzz = primBuffer.data(g1off + 100 * idx + 8);

            auto g1_xxx_zzz = primBuffer.data(g1off + 100 * idx + 9);

            auto g1_xxy_xxx = primBuffer.data(g1off + 100 * idx + 10);

            auto g1_xxy_xxy = primBuffer.data(g1off + 100 * idx + 11);

            auto g1_xxy_xxz = primBuffer.data(g1off + 100 * idx + 12);

            auto g1_xxy_xyy = primBuffer.data(g1off + 100 * idx + 13);

            auto g1_xxy_xyz = primBuffer.data(g1off + 100 * idx + 14);

            auto g1_xxy_xzz = primBuffer.data(g1off + 100 * idx + 15);

            auto g1_xxy_yyy = primBuffer.data(g1off + 100 * idx + 16);

            auto g1_xxy_yyz = primBuffer.data(g1off + 100 * idx + 17);

            auto g1_xxy_yzz = primBuffer.data(g1off + 100 * idx + 18);

            auto g1_xxy_zzz = primBuffer.data(g1off + 100 * idx + 19);

            auto g1_xxz_xxx = primBuffer.data(g1off + 100 * idx + 20);

            auto g1_xxz_xxy = primBuffer.data(g1off + 100 * idx + 21);

            auto g1_xxz_xxz = primBuffer.data(g1off + 100 * idx + 22);

            auto g1_xxz_xyy = primBuffer.data(g1off + 100 * idx + 23);

            auto g1_xxz_xyz = primBuffer.data(g1off + 100 * idx + 24);

            auto g1_xxz_xzz = primBuffer.data(g1off + 100 * idx + 25);

            auto g1_xxz_yyy = primBuffer.data(g1off + 100 * idx + 26);

            auto g1_xxz_yyz = primBuffer.data(g1off + 100 * idx + 27);

            auto g1_xxz_yzz = primBuffer.data(g1off + 100 * idx + 28);

            auto g1_xxz_zzz = primBuffer.data(g1off + 100 * idx + 29);

            auto g1_xyy_xxx = primBuffer.data(g1off + 100 * idx + 30);

            auto g1_xyy_xxy = primBuffer.data(g1off + 100 * idx + 31);

            auto g1_xyy_xxz = primBuffer.data(g1off + 100 * idx + 32);

            auto g1_xyy_xyy = primBuffer.data(g1off + 100 * idx + 33);

            auto g1_xyy_xyz = primBuffer.data(g1off + 100 * idx + 34);

            auto g1_xyy_xzz = primBuffer.data(g1off + 100 * idx + 35);

            auto g1_xyy_yyy = primBuffer.data(g1off + 100 * idx + 36);

            auto g1_xyy_yyz = primBuffer.data(g1off + 100 * idx + 37);

            auto g1_xyy_yzz = primBuffer.data(g1off + 100 * idx + 38);

            auto g1_xyy_zzz = primBuffer.data(g1off + 100 * idx + 39);

            auto g1_xyz_xxx = primBuffer.data(g1off + 100 * idx + 40);

            auto g1_xyz_xxy = primBuffer.data(g1off + 100 * idx + 41);

            auto g1_xyz_xxz = primBuffer.data(g1off + 100 * idx + 42);

            auto g1_xyz_xyy = primBuffer.data(g1off + 100 * idx + 43);

            auto g1_xyz_xyz = primBuffer.data(g1off + 100 * idx + 44);

            auto g1_xyz_xzz = primBuffer.data(g1off + 100 * idx + 45);

            auto g1_xyz_yyy = primBuffer.data(g1off + 100 * idx + 46);

            auto g1_xyz_yyz = primBuffer.data(g1off + 100 * idx + 47);

            auto g1_xyz_yzz = primBuffer.data(g1off + 100 * idx + 48);

            auto g1_xyz_zzz = primBuffer.data(g1off + 100 * idx + 49);

            auto g1_xzz_xxx = primBuffer.data(g1off + 100 * idx + 50);

            auto g1_xzz_xxy = primBuffer.data(g1off + 100 * idx + 51);

            auto g1_xzz_xxz = primBuffer.data(g1off + 100 * idx + 52);

            auto g1_xzz_xyy = primBuffer.data(g1off + 100 * idx + 53);

            auto g1_xzz_xyz = primBuffer.data(g1off + 100 * idx + 54);

            auto g1_xzz_xzz = primBuffer.data(g1off + 100 * idx + 55);

            auto g1_xzz_yyy = primBuffer.data(g1off + 100 * idx + 56);

            auto g1_xzz_yyz = primBuffer.data(g1off + 100 * idx + 57);

            auto g1_xzz_yzz = primBuffer.data(g1off + 100 * idx + 58);

            auto g1_xzz_zzz = primBuffer.data(g1off + 100 * idx + 59);

            auto g1_yyy_xxx = primBuffer.data(g1off + 100 * idx + 60);

            auto g1_yyy_xxy = primBuffer.data(g1off + 100 * idx + 61);

            auto g1_yyy_xxz = primBuffer.data(g1off + 100 * idx + 62);

            auto g1_yyy_xyy = primBuffer.data(g1off + 100 * idx + 63);

            auto g1_yyy_xyz = primBuffer.data(g1off + 100 * idx + 64);

            auto g1_yyy_xzz = primBuffer.data(g1off + 100 * idx + 65);

            auto g1_yyy_yyy = primBuffer.data(g1off + 100 * idx + 66);

            auto g1_yyy_yyz = primBuffer.data(g1off + 100 * idx + 67);

            auto g1_yyy_yzz = primBuffer.data(g1off + 100 * idx + 68);

            auto g1_yyy_zzz = primBuffer.data(g1off + 100 * idx + 69);

            auto g1_yyz_xxx = primBuffer.data(g1off + 100 * idx + 70);

            auto g1_yyz_xxy = primBuffer.data(g1off + 100 * idx + 71);

            auto g1_yyz_xxz = primBuffer.data(g1off + 100 * idx + 72);

            auto g1_yyz_xyy = primBuffer.data(g1off + 100 * idx + 73);

            auto g1_yyz_xyz = primBuffer.data(g1off + 100 * idx + 74);

            auto g1_yyz_xzz = primBuffer.data(g1off + 100 * idx + 75);

            auto g1_yyz_yyy = primBuffer.data(g1off + 100 * idx + 76);

            auto g1_yyz_yyz = primBuffer.data(g1off + 100 * idx + 77);

            auto g1_yyz_yzz = primBuffer.data(g1off + 100 * idx + 78);

            auto g1_yyz_zzz = primBuffer.data(g1off + 100 * idx + 79);

            auto g1_yzz_xxx = primBuffer.data(g1off + 100 * idx + 80);

            auto g1_yzz_xxy = primBuffer.data(g1off + 100 * idx + 81);

            auto g1_yzz_xxz = primBuffer.data(g1off + 100 * idx + 82);

            auto g1_yzz_xyy = primBuffer.data(g1off + 100 * idx + 83);

            auto g1_yzz_xyz = primBuffer.data(g1off + 100 * idx + 84);

            auto g1_yzz_xzz = primBuffer.data(g1off + 100 * idx + 85);

            auto g1_yzz_yyy = primBuffer.data(g1off + 100 * idx + 86);

            auto g1_yzz_yyz = primBuffer.data(g1off + 100 * idx + 87);

            auto g1_yzz_yzz = primBuffer.data(g1off + 100 * idx + 88);

            auto g1_yzz_zzz = primBuffer.data(g1off + 100 * idx + 89);

            auto g1_zzz_xxx = primBuffer.data(g1off + 100 * idx + 90);

            auto g1_zzz_xxy = primBuffer.data(g1off + 100 * idx + 91);

            auto g1_zzz_xxz = primBuffer.data(g1off + 100 * idx + 92);

            auto g1_zzz_xyy = primBuffer.data(g1off + 100 * idx + 93);

            auto g1_zzz_xyz = primBuffer.data(g1off + 100 * idx + 94);

            auto g1_zzz_xzz = primBuffer.data(g1off + 100 * idx + 95);

            auto g1_zzz_yyy = primBuffer.data(g1off + 100 * idx + 96);

            auto g1_zzz_yyz = primBuffer.data(g1off + 100 * idx + 97);

            auto g1_zzz_yzz = primBuffer.data(g1off + 100 * idx + 98);

            auto g1_zzz_zzz = primBuffer.data(g1off + 100 * idx + 99);

            // set up pointers to (G|g(r,r')|F)^(m) integrals

            auto g_xxxx_xxx = primBuffer.data(goff + 150 * idx);

            auto g_xxxx_xxy = primBuffer.data(goff + 150 * idx + 1);

            auto g_xxxx_xxz = primBuffer.data(goff + 150 * idx + 2);

            auto g_xxxx_xyy = primBuffer.data(goff + 150 * idx + 3);

            auto g_xxxx_xyz = primBuffer.data(goff + 150 * idx + 4);

            auto g_xxxx_xzz = primBuffer.data(goff + 150 * idx + 5);

            auto g_xxxx_yyy = primBuffer.data(goff + 150 * idx + 6);

            auto g_xxxx_yyz = primBuffer.data(goff + 150 * idx + 7);

            auto g_xxxx_yzz = primBuffer.data(goff + 150 * idx + 8);

            auto g_xxxx_zzz = primBuffer.data(goff + 150 * idx + 9);

            auto g_xxxy_xxx = primBuffer.data(goff + 150 * idx + 10);

            auto g_xxxy_xxy = primBuffer.data(goff + 150 * idx + 11);

            auto g_xxxy_xxz = primBuffer.data(goff + 150 * idx + 12);

            auto g_xxxy_xyy = primBuffer.data(goff + 150 * idx + 13);

            auto g_xxxy_xyz = primBuffer.data(goff + 150 * idx + 14);

            auto g_xxxy_xzz = primBuffer.data(goff + 150 * idx + 15);

            auto g_xxxy_yyy = primBuffer.data(goff + 150 * idx + 16);

            auto g_xxxy_yyz = primBuffer.data(goff + 150 * idx + 17);

            auto g_xxxy_yzz = primBuffer.data(goff + 150 * idx + 18);

            auto g_xxxy_zzz = primBuffer.data(goff + 150 * idx + 19);

            auto g_xxxz_xxx = primBuffer.data(goff + 150 * idx + 20);

            auto g_xxxz_xxy = primBuffer.data(goff + 150 * idx + 21);

            auto g_xxxz_xxz = primBuffer.data(goff + 150 * idx + 22);

            auto g_xxxz_xyy = primBuffer.data(goff + 150 * idx + 23);

            auto g_xxxz_xyz = primBuffer.data(goff + 150 * idx + 24);

            auto g_xxxz_xzz = primBuffer.data(goff + 150 * idx + 25);

            auto g_xxxz_yyy = primBuffer.data(goff + 150 * idx + 26);

            auto g_xxxz_yyz = primBuffer.data(goff + 150 * idx + 27);

            auto g_xxxz_yzz = primBuffer.data(goff + 150 * idx + 28);

            auto g_xxxz_zzz = primBuffer.data(goff + 150 * idx + 29);

            auto g_xxyy_xxx = primBuffer.data(goff + 150 * idx + 30);

            auto g_xxyy_xxy = primBuffer.data(goff + 150 * idx + 31);

            auto g_xxyy_xxz = primBuffer.data(goff + 150 * idx + 32);

            auto g_xxyy_xyy = primBuffer.data(goff + 150 * idx + 33);

            auto g_xxyy_xyz = primBuffer.data(goff + 150 * idx + 34);

            auto g_xxyy_xzz = primBuffer.data(goff + 150 * idx + 35);

            auto g_xxyy_yyy = primBuffer.data(goff + 150 * idx + 36);

            auto g_xxyy_yyz = primBuffer.data(goff + 150 * idx + 37);

            auto g_xxyy_yzz = primBuffer.data(goff + 150 * idx + 38);

            auto g_xxyy_zzz = primBuffer.data(goff + 150 * idx + 39);

            auto g_xxyz_xxx = primBuffer.data(goff + 150 * idx + 40);

            auto g_xxyz_xxy = primBuffer.data(goff + 150 * idx + 41);

            auto g_xxyz_xxz = primBuffer.data(goff + 150 * idx + 42);

            auto g_xxyz_xyy = primBuffer.data(goff + 150 * idx + 43);

            auto g_xxyz_xyz = primBuffer.data(goff + 150 * idx + 44);

            auto g_xxyz_xzz = primBuffer.data(goff + 150 * idx + 45);

            auto g_xxyz_yyy = primBuffer.data(goff + 150 * idx + 46);

            auto g_xxyz_yyz = primBuffer.data(goff + 150 * idx + 47);

            auto g_xxyz_yzz = primBuffer.data(goff + 150 * idx + 48);

            auto g_xxyz_zzz = primBuffer.data(goff + 150 * idx + 49);

            auto g_xxzz_xxx = primBuffer.data(goff + 150 * idx + 50);

            auto g_xxzz_xxy = primBuffer.data(goff + 150 * idx + 51);

            auto g_xxzz_xxz = primBuffer.data(goff + 150 * idx + 52);

            auto g_xxzz_xyy = primBuffer.data(goff + 150 * idx + 53);

            auto g_xxzz_xyz = primBuffer.data(goff + 150 * idx + 54);

            auto g_xxzz_xzz = primBuffer.data(goff + 150 * idx + 55);

            auto g_xxzz_yyy = primBuffer.data(goff + 150 * idx + 56);

            auto g_xxzz_yyz = primBuffer.data(goff + 150 * idx + 57);

            auto g_xxzz_yzz = primBuffer.data(goff + 150 * idx + 58);

            auto g_xxzz_zzz = primBuffer.data(goff + 150 * idx + 59);

            auto g_xyyy_xxx = primBuffer.data(goff + 150 * idx + 60);

            auto g_xyyy_xxy = primBuffer.data(goff + 150 * idx + 61);

            auto g_xyyy_xxz = primBuffer.data(goff + 150 * idx + 62);

            auto g_xyyy_xyy = primBuffer.data(goff + 150 * idx + 63);

            auto g_xyyy_xyz = primBuffer.data(goff + 150 * idx + 64);

            auto g_xyyy_xzz = primBuffer.data(goff + 150 * idx + 65);

            auto g_xyyy_yyy = primBuffer.data(goff + 150 * idx + 66);

            auto g_xyyy_yyz = primBuffer.data(goff + 150 * idx + 67);

            auto g_xyyy_yzz = primBuffer.data(goff + 150 * idx + 68);

            auto g_xyyy_zzz = primBuffer.data(goff + 150 * idx + 69);

            auto g_xyyz_xxx = primBuffer.data(goff + 150 * idx + 70);

            auto g_xyyz_xxy = primBuffer.data(goff + 150 * idx + 71);

            auto g_xyyz_xxz = primBuffer.data(goff + 150 * idx + 72);

            auto g_xyyz_xyy = primBuffer.data(goff + 150 * idx + 73);

            auto g_xyyz_xyz = primBuffer.data(goff + 150 * idx + 74);

            auto g_xyyz_xzz = primBuffer.data(goff + 150 * idx + 75);

            auto g_xyyz_yyy = primBuffer.data(goff + 150 * idx + 76);

            auto g_xyyz_yyz = primBuffer.data(goff + 150 * idx + 77);

            auto g_xyyz_yzz = primBuffer.data(goff + 150 * idx + 78);

            auto g_xyyz_zzz = primBuffer.data(goff + 150 * idx + 79);

            auto g_xyzz_xxx = primBuffer.data(goff + 150 * idx + 80);

            auto g_xyzz_xxy = primBuffer.data(goff + 150 * idx + 81);

            auto g_xyzz_xxz = primBuffer.data(goff + 150 * idx + 82);

            auto g_xyzz_xyy = primBuffer.data(goff + 150 * idx + 83);

            auto g_xyzz_xyz = primBuffer.data(goff + 150 * idx + 84);

            auto g_xyzz_xzz = primBuffer.data(goff + 150 * idx + 85);

            auto g_xyzz_yyy = primBuffer.data(goff + 150 * idx + 86);

            auto g_xyzz_yyz = primBuffer.data(goff + 150 * idx + 87);

            auto g_xyzz_yzz = primBuffer.data(goff + 150 * idx + 88);

            auto g_xyzz_zzz = primBuffer.data(goff + 150 * idx + 89);

            auto g_xzzz_xxx = primBuffer.data(goff + 150 * idx + 90);

            auto g_xzzz_xxy = primBuffer.data(goff + 150 * idx + 91);

            auto g_xzzz_xxz = primBuffer.data(goff + 150 * idx + 92);

            auto g_xzzz_xyy = primBuffer.data(goff + 150 * idx + 93);

            auto g_xzzz_xyz = primBuffer.data(goff + 150 * idx + 94);

            auto g_xzzz_xzz = primBuffer.data(goff + 150 * idx + 95);

            auto g_xzzz_yyy = primBuffer.data(goff + 150 * idx + 96);

            auto g_xzzz_yyz = primBuffer.data(goff + 150 * idx + 97);

            auto g_xzzz_yzz = primBuffer.data(goff + 150 * idx + 98);

            auto g_xzzz_zzz = primBuffer.data(goff + 150 * idx + 99);

            auto g_yyyy_xxx = primBuffer.data(goff + 150 * idx + 100);

            auto g_yyyy_xxy = primBuffer.data(goff + 150 * idx + 101);

            auto g_yyyy_xxz = primBuffer.data(goff + 150 * idx + 102);

            auto g_yyyy_xyy = primBuffer.data(goff + 150 * idx + 103);

            auto g_yyyy_xyz = primBuffer.data(goff + 150 * idx + 104);

            auto g_yyyy_xzz = primBuffer.data(goff + 150 * idx + 105);

            auto g_yyyy_yyy = primBuffer.data(goff + 150 * idx + 106);

            auto g_yyyy_yyz = primBuffer.data(goff + 150 * idx + 107);

            auto g_yyyy_yzz = primBuffer.data(goff + 150 * idx + 108);

            auto g_yyyy_zzz = primBuffer.data(goff + 150 * idx + 109);

            auto g_yyyz_xxx = primBuffer.data(goff + 150 * idx + 110);

            auto g_yyyz_xxy = primBuffer.data(goff + 150 * idx + 111);

            auto g_yyyz_xxz = primBuffer.data(goff + 150 * idx + 112);

            auto g_yyyz_xyy = primBuffer.data(goff + 150 * idx + 113);

            auto g_yyyz_xyz = primBuffer.data(goff + 150 * idx + 114);

            auto g_yyyz_xzz = primBuffer.data(goff + 150 * idx + 115);

            auto g_yyyz_yyy = primBuffer.data(goff + 150 * idx + 116);

            auto g_yyyz_yyz = primBuffer.data(goff + 150 * idx + 117);

            auto g_yyyz_yzz = primBuffer.data(goff + 150 * idx + 118);

            auto g_yyyz_zzz = primBuffer.data(goff + 150 * idx + 119);

            auto g_yyzz_xxx = primBuffer.data(goff + 150 * idx + 120);

            auto g_yyzz_xxy = primBuffer.data(goff + 150 * idx + 121);

            auto g_yyzz_xxz = primBuffer.data(goff + 150 * idx + 122);

            auto g_yyzz_xyy = primBuffer.data(goff + 150 * idx + 123);

            auto g_yyzz_xyz = primBuffer.data(goff + 150 * idx + 124);

            auto g_yyzz_xzz = primBuffer.data(goff + 150 * idx + 125);

            auto g_yyzz_yyy = primBuffer.data(goff + 150 * idx + 126);

            auto g_yyzz_yyz = primBuffer.data(goff + 150 * idx + 127);

            auto g_yyzz_yzz = primBuffer.data(goff + 150 * idx + 128);

            auto g_yyzz_zzz = primBuffer.data(goff + 150 * idx + 129);

            auto g_yzzz_xxx = primBuffer.data(goff + 150 * idx + 130);

            auto g_yzzz_xxy = primBuffer.data(goff + 150 * idx + 131);

            auto g_yzzz_xxz = primBuffer.data(goff + 150 * idx + 132);

            auto g_yzzz_xyy = primBuffer.data(goff + 150 * idx + 133);

            auto g_yzzz_xyz = primBuffer.data(goff + 150 * idx + 134);

            auto g_yzzz_xzz = primBuffer.data(goff + 150 * idx + 135);

            auto g_yzzz_yyy = primBuffer.data(goff + 150 * idx + 136);

            auto g_yzzz_yyz = primBuffer.data(goff + 150 * idx + 137);

            auto g_yzzz_yzz = primBuffer.data(goff + 150 * idx + 138);

            auto g_yzzz_zzz = primBuffer.data(goff + 150 * idx + 139);

            auto g_zzzz_xxx = primBuffer.data(goff + 150 * idx + 140);

            auto g_zzzz_xxy = primBuffer.data(goff + 150 * idx + 141);

            auto g_zzzz_xxz = primBuffer.data(goff + 150 * idx + 142);

            auto g_zzzz_xyy = primBuffer.data(goff + 150 * idx + 143);

            auto g_zzzz_xyz = primBuffer.data(goff + 150 * idx + 144);

            auto g_zzzz_xzz = primBuffer.data(goff + 150 * idx + 145);

            auto g_zzzz_yyy = primBuffer.data(goff + 150 * idx + 146);

            auto g_zzzz_yyz = primBuffer.data(goff + 150 * idx + 147);

            auto g_zzzz_yzz = primBuffer.data(goff + 150 * idx + 148);

            auto g_zzzz_zzz = primBuffer.data(goff + 150 * idx + 149);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xxx_xx,\
                                         gk_xxx_xy, gk_xxx_xz, gk_xxx_yy, gk_xxx_yz,\
                                         gk_xxx_zz, gk_xxy_xx, gk_xxy_xy, gk_xxy_xz,\
                                         gk_xxy_yy, gk_xxy_yz, gk_xxy_zz, gk_xxz_xx,\
                                         gk_xxz_xy, gk_xxz_xz, gk_xxz_yy, gk_xxz_yz,\
                                         gk_xxz_zz, gk_xyy_xx, gk_xyy_xy, gk_xyy_xz,\
                                         gk_xyy_yy, gk_xyy_yz, gk_xyy_zz, gk_xyz_xx,\
                                         gk_xyz_xy, gk_xyz_xz, gk_xyz_yy, gk_xyz_yz,\
                                         gk_xyz_zz, gk_xzz_xx, gk_xzz_xy, gk_xzz_xz,\
                                         gk_xzz_yy, gk_xzz_yz, gk_xzz_zz, gk_yyy_xx,\
                                         gk_yyy_xy, gk_yyy_xz, gk_yyy_yy, gk_yyy_yz,\
                                         gk_yyy_zz, gk_yyz_xx, gk_yyz_xy, gk_yyz_xz,\
                                         gk_yyz_yy, gk_yyz_yz, gk_yyz_zz, gk_yzz_xx,\
                                         gk_yzz_xy, gk_yzz_xz, gk_yzz_yy, gk_yzz_yz,\
                                         gk_yzz_zz, gk_zzz_xx, gk_zzz_xy, gk_zzz_xz,\
                                         gk_zzz_yy, gk_zzz_yz, gk_zzz_zz, g20_xx_xxx,\
                                         g20_xx_xxy, g20_xx_xxz, g20_xx_xyy, g20_xx_xyz,\
                                         g20_xx_xzz, g20_xx_yyy, g20_xx_yyz, g20_xx_yzz,\
                                         g20_xx_zzz, g20_xy_xxx, g20_xy_xxy, g20_xy_xxz,\
                                         g20_xy_xyy, g20_xy_xyz, g20_xy_xzz, g20_xy_yyy,\
                                         g20_xy_yyz, g20_xy_yzz, g20_xy_zzz, g20_xz_xxx,\
                                         g20_xz_xxy, g20_xz_xxz, g20_xz_xyy, g20_xz_xyz,\
                                         g20_xz_xzz, g20_xz_yyy, g20_xz_yyz, g20_xz_yzz,\
                                         g20_xz_zzz, g20_yy_xxx, g20_yy_xxy, g20_yy_xxz,\
                                         g20_yy_xyy, g20_yy_xyz, g20_yy_xzz, g20_yy_yyy,\
                                         g20_yy_yyz, g20_yy_yzz, g20_yy_zzz, g20_yz_xxx,\
                                         g20_yz_xxy, g20_yz_xxz, g20_yz_xyy, g20_yz_xyz,\
                                         g20_yz_xzz, g20_yz_yyy, g20_yz_yyz, g20_yz_yzz,\
                                         g20_yz_zzz, g20_zz_xxx, g20_zz_xxy, g20_zz_xxz,\
                                         g20_zz_xyy, g20_zz_xyz, g20_zz_xzz, g20_zz_yyy,\
                                         g20_zz_yyz, g20_zz_yzz, g20_zz_zzz, g21_xx_xxx,\
                                         g21_xx_xxy, g21_xx_xxz, g21_xx_xyy, g21_xx_xyz,\
                                         g21_xx_xzz, g21_xx_yyy, g21_xx_yyz, g21_xx_yzz,\
                                         g21_xx_zzz, g21_xy_xxx, g21_xy_xxy, g21_xy_xxz,\
                                         g21_xy_xyy, g21_xy_xyz, g21_xy_xzz, g21_xy_yyy,\
                                         g21_xy_yyz, g21_xy_yzz, g21_xy_zzz, g21_xz_xxx,\
                                         g21_xz_xxy, g21_xz_xxz, g21_xz_xyy, g21_xz_xyz,\
                                         g21_xz_xzz, g21_xz_yyy, g21_xz_yyz, g21_xz_yzz,\
                                         g21_xz_zzz, g21_yy_xxx, g21_yy_xxy, g21_yy_xxz,\
                                         g21_yy_xyy, g21_yy_xyz, g21_yy_xzz, g21_yy_yyy,\
                                         g21_yy_yyz, g21_yy_yzz, g21_yy_zzz, g21_yz_xxx,\
                                         g21_yz_xxy, g21_yz_xxz, g21_yz_xyy, g21_yz_xyz,\
                                         g21_yz_xzz, g21_yz_yyy, g21_yz_yyz, g21_yz_yzz,\
                                         g21_yz_zzz, g21_zz_xxx, g21_zz_xxy, g21_zz_xxz,\
                                         g21_zz_xyy, g21_zz_xyz, g21_zz_xzz, g21_zz_yyy,\
                                         g21_zz_yyz, g21_zz_yzz, g21_zz_zzz, g1_xxx_xxx,\
                                         g1_xxx_xxy, g1_xxx_xxz, g1_xxx_xyy, g1_xxx_xyz,\
                                         g1_xxx_xzz, g1_xxx_yyy, g1_xxx_yyz, g1_xxx_yzz,\
                                         g1_xxx_zzz, g1_xxy_xxx, g1_xxy_xxy, g1_xxy_xxz,\
                                         g1_xxy_xyy, g1_xxy_xyz, g1_xxy_xzz, g1_xxy_yyy,\
                                         g1_xxy_yyz, g1_xxy_yzz, g1_xxy_zzz, g1_xxz_xxx,\
                                         g1_xxz_xxy, g1_xxz_xxz, g1_xxz_xyy, g1_xxz_xyz,\
                                         g1_xxz_xzz, g1_xxz_yyy, g1_xxz_yyz, g1_xxz_yzz,\
                                         g1_xxz_zzz, g1_xyy_xxx, g1_xyy_xxy, g1_xyy_xxz,\
                                         g1_xyy_xyy, g1_xyy_xyz, g1_xyy_xzz, g1_xyy_yyy,\
                                         g1_xyy_yyz, g1_xyy_yzz, g1_xyy_zzz, g1_xyz_xxx,\
                                         g1_xyz_xxy, g1_xyz_xxz, g1_xyz_xyy, g1_xyz_xyz,\
                                         g1_xyz_xzz, g1_xyz_yyy, g1_xyz_yyz, g1_xyz_yzz,\
                                         g1_xyz_zzz, g1_xzz_xxx, g1_xzz_xxy, g1_xzz_xxz,\
                                         g1_xzz_xyy, g1_xzz_xyz, g1_xzz_xzz, g1_xzz_yyy,\
                                         g1_xzz_yyz, g1_xzz_yzz, g1_xzz_zzz, g1_yyy_xxx,\
                                         g1_yyy_xxy, g1_yyy_xxz, g1_yyy_xyy, g1_yyy_xyz,\
                                         g1_yyy_xzz, g1_yyy_yyy, g1_yyy_yyz, g1_yyy_yzz,\
                                         g1_yyy_zzz, g1_yyz_xxx, g1_yyz_xxy, g1_yyz_xxz,\
                                         g1_yyz_xyy, g1_yyz_xyz, g1_yyz_xzz, g1_yyz_yyy,\
                                         g1_yyz_yyz, g1_yyz_yzz, g1_yyz_zzz, g1_yzz_xxx,\
                                         g1_yzz_xxy, g1_yzz_xxz, g1_yzz_xyy, g1_yzz_xyz,\
                                         g1_yzz_xzz, g1_yzz_yyy, g1_yzz_yyz, g1_yzz_yzz,\
                                         g1_yzz_zzz, g1_zzz_xxx, g1_zzz_xxy, g1_zzz_xxz,\
                                         g1_zzz_xyy, g1_zzz_xyz, g1_zzz_xzz, g1_zzz_yyy,\
                                         g1_zzz_yyz, g1_zzz_yzz, g1_zzz_zzz, g_xxxx_xxx,\
                                         g_xxxx_xxy, g_xxxx_xxz, g_xxxx_xyy, g_xxxx_xyz,\
                                         g_xxxx_xzz, g_xxxx_yyy, g_xxxx_yyz, g_xxxx_yzz,\
                                         g_xxxx_zzz, g_xxxy_xxx, g_xxxy_xxy, g_xxxy_xxz,\
                                         g_xxxy_xyy, g_xxxy_xyz, g_xxxy_xzz, g_xxxy_yyy,\
                                         g_xxxy_yyz, g_xxxy_yzz, g_xxxy_zzz, g_xxxz_xxx,\
                                         g_xxxz_xxy, g_xxxz_xxz, g_xxxz_xyy, g_xxxz_xyz,\
                                         g_xxxz_xzz, g_xxxz_yyy, g_xxxz_yyz, g_xxxz_yzz,\
                                         g_xxxz_zzz, g_xxyy_xxx, g_xxyy_xxy, g_xxyy_xxz,\
                                         g_xxyy_xyy, g_xxyy_xyz, g_xxyy_xzz, g_xxyy_yyy,\
                                         g_xxyy_yyz, g_xxyy_yzz, g_xxyy_zzz, g_xxyz_xxx,\
                                         g_xxyz_xxy, g_xxyz_xxz, g_xxyz_xyy, g_xxyz_xyz,\
                                         g_xxyz_xzz, g_xxyz_yyy, g_xxyz_yyz, g_xxyz_yzz,\
                                         g_xxyz_zzz, g_xxzz_xxx, g_xxzz_xxy, g_xxzz_xxz,\
                                         g_xxzz_xyy, g_xxzz_xyz, g_xxzz_xzz, g_xxzz_yyy,\
                                         g_xxzz_yyz, g_xxzz_yzz, g_xxzz_zzz, g_xyyy_xxx,\
                                         g_xyyy_xxy, g_xyyy_xxz, g_xyyy_xyy, g_xyyy_xyz,\
                                         g_xyyy_xzz, g_xyyy_yyy, g_xyyy_yyz, g_xyyy_yzz,\
                                         g_xyyy_zzz, g_xyyz_xxx, g_xyyz_xxy, g_xyyz_xxz,\
                                         g_xyyz_xyy, g_xyyz_xyz, g_xyyz_xzz, g_xyyz_yyy,\
                                         g_xyyz_yyz, g_xyyz_yzz, g_xyyz_zzz, g_xyzz_xxx,\
                                         g_xyzz_xxy, g_xyzz_xxz, g_xyzz_xyy, g_xyzz_xyz,\
                                         g_xyzz_xzz, g_xyzz_yyy, g_xyzz_yyz, g_xyzz_yzz,\
                                         g_xyzz_zzz, g_xzzz_xxx, g_xzzz_xxy, g_xzzz_xxz,\
                                         g_xzzz_xyy, g_xzzz_xyz, g_xzzz_xzz, g_xzzz_yyy,\
                                         g_xzzz_yyz, g_xzzz_yzz, g_xzzz_zzz, g_yyyy_xxx,\
                                         g_yyyy_xxy, g_yyyy_xxz, g_yyyy_xyy, g_yyyy_xyz,\
                                         g_yyyy_xzz, g_yyyy_yyy, g_yyyy_yyz, g_yyyy_yzz,\
                                         g_yyyy_zzz, g_yyyz_xxx, g_yyyz_xxy, g_yyyz_xxz,\
                                         g_yyyz_xyy, g_yyyz_xyz, g_yyyz_xzz, g_yyyz_yyy,\
                                         g_yyyz_yyz, g_yyyz_yzz, g_yyyz_zzz, g_yyzz_xxx,\
                                         g_yyzz_xxy, g_yyzz_xxz, g_yyzz_xyy, g_yyzz_xyz,\
                                         g_yyzz_xzz, g_yyzz_yyy, g_yyzz_yyz, g_yyzz_yzz,\
                                         g_yyzz_zzz, g_yzzz_xxx, g_yzzz_xxy, g_yzzz_xxz,\
                                         g_yzzz_xyy, g_yzzz_xyz, g_yzzz_xzz, g_yzzz_yyy,\
                                         g_yzzz_yyz, g_yzzz_yzz, g_yzzz_zzz, g_zzzz_xxx,\
                                         g_zzzz_xxy, g_zzzz_xxz, g_zzzz_xyy, g_zzzz_xyz,\
                                         g_zzzz_xzz, g_zzzz_yyy, g_zzzz_yyz, g_zzzz_yzz,\
                                         g_zzzz_zzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxxx_xxx[k] = fr * g1_xxx_xxx[k] + f2g * (3.0 * g20_xx_xxx[k] - 3.0 * fgz * g21_xx_xxx[k]) + 3.0 * f2t * gk_xxx_xx[k];

                g_xxxx_xxy[k] = fr * g1_xxx_xxy[k] + f2g * (3.0 * g20_xx_xxy[k] - 3.0 * fgz * g21_xx_xxy[k]) + 2.0 * f2t * gk_xxx_xy[k];

                g_xxxx_xxz[k] = fr * g1_xxx_xxz[k] + f2g * (3.0 * g20_xx_xxz[k] - 3.0 * fgz * g21_xx_xxz[k]) + 2.0 * f2t * gk_xxx_xz[k];

                g_xxxx_xyy[k] = fr * g1_xxx_xyy[k] + f2g * (3.0 * g20_xx_xyy[k] - 3.0 * fgz * g21_xx_xyy[k]) + f2t * gk_xxx_yy[k];

                g_xxxx_xyz[k] = fr * g1_xxx_xyz[k] + f2g * (3.0 * g20_xx_xyz[k] - 3.0 * fgz * g21_xx_xyz[k]) + f2t * gk_xxx_yz[k];

                g_xxxx_xzz[k] = fr * g1_xxx_xzz[k] + f2g * (3.0 * g20_xx_xzz[k] - 3.0 * fgz * g21_xx_xzz[k]) + f2t * gk_xxx_zz[k];

                g_xxxx_yyy[k] = fr * g1_xxx_yyy[k] + f2g * (3.0 * g20_xx_yyy[k] - 3.0 * fgz * g21_xx_yyy[k]);

                g_xxxx_yyz[k] = fr * g1_xxx_yyz[k] + f2g * (3.0 * g20_xx_yyz[k] - 3.0 * fgz * g21_xx_yyz[k]);

                g_xxxx_yzz[k] = fr * g1_xxx_yzz[k] + f2g * (3.0 * g20_xx_yzz[k] - 3.0 * fgz * g21_xx_yzz[k]);

                g_xxxx_zzz[k] = fr * g1_xxx_zzz[k] + f2g * (3.0 * g20_xx_zzz[k] - 3.0 * fgz * g21_xx_zzz[k]);

                g_xxxy_xxx[k] = fr * g1_xxy_xxx[k] + f2g * (2.0 * g20_xy_xxx[k] - 2.0 * fgz * g21_xy_xxx[k]) + 3.0 * f2t * gk_xxy_xx[k];

                g_xxxy_xxy[k] = fr * g1_xxy_xxy[k] + f2g * (2.0 * g20_xy_xxy[k] - 2.0 * fgz * g21_xy_xxy[k]) + 2.0 * f2t * gk_xxy_xy[k];

                g_xxxy_xxz[k] = fr * g1_xxy_xxz[k] + f2g * (2.0 * g20_xy_xxz[k] - 2.0 * fgz * g21_xy_xxz[k]) + 2.0 * f2t * gk_xxy_xz[k];

                g_xxxy_xyy[k] = fr * g1_xxy_xyy[k] + f2g * (2.0 * g20_xy_xyy[k] - 2.0 * fgz * g21_xy_xyy[k]) + f2t * gk_xxy_yy[k];

                g_xxxy_xyz[k] = fr * g1_xxy_xyz[k] + f2g * (2.0 * g20_xy_xyz[k] - 2.0 * fgz * g21_xy_xyz[k]) + f2t * gk_xxy_yz[k];

                g_xxxy_xzz[k] = fr * g1_xxy_xzz[k] + f2g * (2.0 * g20_xy_xzz[k] - 2.0 * fgz * g21_xy_xzz[k]) + f2t * gk_xxy_zz[k];

                g_xxxy_yyy[k] = fr * g1_xxy_yyy[k] + f2g * (2.0 * g20_xy_yyy[k] - 2.0 * fgz * g21_xy_yyy[k]);

                g_xxxy_yyz[k] = fr * g1_xxy_yyz[k] + f2g * (2.0 * g20_xy_yyz[k] - 2.0 * fgz * g21_xy_yyz[k]);

                g_xxxy_yzz[k] = fr * g1_xxy_yzz[k] + f2g * (2.0 * g20_xy_yzz[k] - 2.0 * fgz * g21_xy_yzz[k]);

                g_xxxy_zzz[k] = fr * g1_xxy_zzz[k] + f2g * (2.0 * g20_xy_zzz[k] - 2.0 * fgz * g21_xy_zzz[k]);

                g_xxxz_xxx[k] = fr * g1_xxz_xxx[k] + f2g * (2.0 * g20_xz_xxx[k] - 2.0 * fgz * g21_xz_xxx[k]) + 3.0 * f2t * gk_xxz_xx[k];

                g_xxxz_xxy[k] = fr * g1_xxz_xxy[k] + f2g * (2.0 * g20_xz_xxy[k] - 2.0 * fgz * g21_xz_xxy[k]) + 2.0 * f2t * gk_xxz_xy[k];

                g_xxxz_xxz[k] = fr * g1_xxz_xxz[k] + f2g * (2.0 * g20_xz_xxz[k] - 2.0 * fgz * g21_xz_xxz[k]) + 2.0 * f2t * gk_xxz_xz[k];

                g_xxxz_xyy[k] = fr * g1_xxz_xyy[k] + f2g * (2.0 * g20_xz_xyy[k] - 2.0 * fgz * g21_xz_xyy[k]) + f2t * gk_xxz_yy[k];

                g_xxxz_xyz[k] = fr * g1_xxz_xyz[k] + f2g * (2.0 * g20_xz_xyz[k] - 2.0 * fgz * g21_xz_xyz[k]) + f2t * gk_xxz_yz[k];

                g_xxxz_xzz[k] = fr * g1_xxz_xzz[k] + f2g * (2.0 * g20_xz_xzz[k] - 2.0 * fgz * g21_xz_xzz[k]) + f2t * gk_xxz_zz[k];

                g_xxxz_yyy[k] = fr * g1_xxz_yyy[k] + f2g * (2.0 * g20_xz_yyy[k] - 2.0 * fgz * g21_xz_yyy[k]);

                g_xxxz_yyz[k] = fr * g1_xxz_yyz[k] + f2g * (2.0 * g20_xz_yyz[k] - 2.0 * fgz * g21_xz_yyz[k]);

                g_xxxz_yzz[k] = fr * g1_xxz_yzz[k] + f2g * (2.0 * g20_xz_yzz[k] - 2.0 * fgz * g21_xz_yzz[k]);

                g_xxxz_zzz[k] = fr * g1_xxz_zzz[k] + f2g * (2.0 * g20_xz_zzz[k] - 2.0 * fgz * g21_xz_zzz[k]);

                g_xxyy_xxx[k] = fr * g1_xyy_xxx[k] + f2g * (g20_yy_xxx[k] - fgz * g21_yy_xxx[k]) + 3.0 * f2t * gk_xyy_xx[k];

                g_xxyy_xxy[k] = fr * g1_xyy_xxy[k] + f2g * (g20_yy_xxy[k] - fgz * g21_yy_xxy[k]) + 2.0 * f2t * gk_xyy_xy[k];

                g_xxyy_xxz[k] = fr * g1_xyy_xxz[k] + f2g * (g20_yy_xxz[k] - fgz * g21_yy_xxz[k]) + 2.0 * f2t * gk_xyy_xz[k];

                g_xxyy_xyy[k] = fr * g1_xyy_xyy[k] + f2g * (g20_yy_xyy[k] - fgz * g21_yy_xyy[k]) + f2t * gk_xyy_yy[k];

                g_xxyy_xyz[k] = fr * g1_xyy_xyz[k] + f2g * (g20_yy_xyz[k] - fgz * g21_yy_xyz[k]) + f2t * gk_xyy_yz[k];

                g_xxyy_xzz[k] = fr * g1_xyy_xzz[k] + f2g * (g20_yy_xzz[k] - fgz * g21_yy_xzz[k]) + f2t * gk_xyy_zz[k];

                g_xxyy_yyy[k] = fr * g1_xyy_yyy[k] + f2g * (g20_yy_yyy[k] - fgz * g21_yy_yyy[k]);

                g_xxyy_yyz[k] = fr * g1_xyy_yyz[k] + f2g * (g20_yy_yyz[k] - fgz * g21_yy_yyz[k]);

                g_xxyy_yzz[k] = fr * g1_xyy_yzz[k] + f2g * (g20_yy_yzz[k] - fgz * g21_yy_yzz[k]);

                g_xxyy_zzz[k] = fr * g1_xyy_zzz[k] + f2g * (g20_yy_zzz[k] - fgz * g21_yy_zzz[k]);

                g_xxyz_xxx[k] = fr * g1_xyz_xxx[k] + f2g * (g20_yz_xxx[k] - fgz * g21_yz_xxx[k]) + 3.0 * f2t * gk_xyz_xx[k];

                g_xxyz_xxy[k] = fr * g1_xyz_xxy[k] + f2g * (g20_yz_xxy[k] - fgz * g21_yz_xxy[k]) + 2.0 * f2t * gk_xyz_xy[k];

                g_xxyz_xxz[k] = fr * g1_xyz_xxz[k] + f2g * (g20_yz_xxz[k] - fgz * g21_yz_xxz[k]) + 2.0 * f2t * gk_xyz_xz[k];

                g_xxyz_xyy[k] = fr * g1_xyz_xyy[k] + f2g * (g20_yz_xyy[k] - fgz * g21_yz_xyy[k]) + f2t * gk_xyz_yy[k];

                g_xxyz_xyz[k] = fr * g1_xyz_xyz[k] + f2g * (g20_yz_xyz[k] - fgz * g21_yz_xyz[k]) + f2t * gk_xyz_yz[k];

                g_xxyz_xzz[k] = fr * g1_xyz_xzz[k] + f2g * (g20_yz_xzz[k] - fgz * g21_yz_xzz[k]) + f2t * gk_xyz_zz[k];

                g_xxyz_yyy[k] = fr * g1_xyz_yyy[k] + f2g * (g20_yz_yyy[k] - fgz * g21_yz_yyy[k]);

                g_xxyz_yyz[k] = fr * g1_xyz_yyz[k] + f2g * (g20_yz_yyz[k] - fgz * g21_yz_yyz[k]);

                g_xxyz_yzz[k] = fr * g1_xyz_yzz[k] + f2g * (g20_yz_yzz[k] - fgz * g21_yz_yzz[k]);

                g_xxyz_zzz[k] = fr * g1_xyz_zzz[k] + f2g * (g20_yz_zzz[k] - fgz * g21_yz_zzz[k]);

                g_xxzz_xxx[k] = fr * g1_xzz_xxx[k] + f2g * (g20_zz_xxx[k] - fgz * g21_zz_xxx[k]) + 3.0 * f2t * gk_xzz_xx[k];

                g_xxzz_xxy[k] = fr * g1_xzz_xxy[k] + f2g * (g20_zz_xxy[k] - fgz * g21_zz_xxy[k]) + 2.0 * f2t * gk_xzz_xy[k];

                g_xxzz_xxz[k] = fr * g1_xzz_xxz[k] + f2g * (g20_zz_xxz[k] - fgz * g21_zz_xxz[k]) + 2.0 * f2t * gk_xzz_xz[k];

                g_xxzz_xyy[k] = fr * g1_xzz_xyy[k] + f2g * (g20_zz_xyy[k] - fgz * g21_zz_xyy[k]) + f2t * gk_xzz_yy[k];

                g_xxzz_xyz[k] = fr * g1_xzz_xyz[k] + f2g * (g20_zz_xyz[k] - fgz * g21_zz_xyz[k]) + f2t * gk_xzz_yz[k];

                g_xxzz_xzz[k] = fr * g1_xzz_xzz[k] + f2g * (g20_zz_xzz[k] - fgz * g21_zz_xzz[k]) + f2t * gk_xzz_zz[k];

                g_xxzz_yyy[k] = fr * g1_xzz_yyy[k] + f2g * (g20_zz_yyy[k] - fgz * g21_zz_yyy[k]);

                g_xxzz_yyz[k] = fr * g1_xzz_yyz[k] + f2g * (g20_zz_yyz[k] - fgz * g21_zz_yyz[k]);

                g_xxzz_yzz[k] = fr * g1_xzz_yzz[k] + f2g * (g20_zz_yzz[k] - fgz * g21_zz_yzz[k]);

                g_xxzz_zzz[k] = fr * g1_xzz_zzz[k] + f2g * (g20_zz_zzz[k] - fgz * g21_zz_zzz[k]);

                g_xyyy_xxx[k] = fr * g1_yyy_xxx[k] + 3.0 * f2t * gk_yyy_xx[k];

                g_xyyy_xxy[k] = fr * g1_yyy_xxy[k] + 2.0 * f2t * gk_yyy_xy[k];

                g_xyyy_xxz[k] = fr * g1_yyy_xxz[k] + 2.0 * f2t * gk_yyy_xz[k];

                g_xyyy_xyy[k] = fr * g1_yyy_xyy[k] + f2t * gk_yyy_yy[k];

                g_xyyy_xyz[k] = fr * g1_yyy_xyz[k] + f2t * gk_yyy_yz[k];

                g_xyyy_xzz[k] = fr * g1_yyy_xzz[k] + f2t * gk_yyy_zz[k];

                g_xyyy_yyy[k] = fr * g1_yyy_yyy[k];

                g_xyyy_yyz[k] = fr * g1_yyy_yyz[k];

                g_xyyy_yzz[k] = fr * g1_yyy_yzz[k];

                g_xyyy_zzz[k] = fr * g1_yyy_zzz[k];

                g_xyyz_xxx[k] = fr * g1_yyz_xxx[k] + 3.0 * f2t * gk_yyz_xx[k];

                g_xyyz_xxy[k] = fr * g1_yyz_xxy[k] + 2.0 * f2t * gk_yyz_xy[k];

                g_xyyz_xxz[k] = fr * g1_yyz_xxz[k] + 2.0 * f2t * gk_yyz_xz[k];

                g_xyyz_xyy[k] = fr * g1_yyz_xyy[k] + f2t * gk_yyz_yy[k];

                g_xyyz_xyz[k] = fr * g1_yyz_xyz[k] + f2t * gk_yyz_yz[k];

                g_xyyz_xzz[k] = fr * g1_yyz_xzz[k] + f2t * gk_yyz_zz[k];

                g_xyyz_yyy[k] = fr * g1_yyz_yyy[k];

                g_xyyz_yyz[k] = fr * g1_yyz_yyz[k];

                g_xyyz_yzz[k] = fr * g1_yyz_yzz[k];

                g_xyyz_zzz[k] = fr * g1_yyz_zzz[k];

                g_xyzz_xxx[k] = fr * g1_yzz_xxx[k] + 3.0 * f2t * gk_yzz_xx[k];

                g_xyzz_xxy[k] = fr * g1_yzz_xxy[k] + 2.0 * f2t * gk_yzz_xy[k];

                g_xyzz_xxz[k] = fr * g1_yzz_xxz[k] + 2.0 * f2t * gk_yzz_xz[k];

                g_xyzz_xyy[k] = fr * g1_yzz_xyy[k] + f2t * gk_yzz_yy[k];

                g_xyzz_xyz[k] = fr * g1_yzz_xyz[k] + f2t * gk_yzz_yz[k];

                g_xyzz_xzz[k] = fr * g1_yzz_xzz[k] + f2t * gk_yzz_zz[k];

                g_xyzz_yyy[k] = fr * g1_yzz_yyy[k];

                g_xyzz_yyz[k] = fr * g1_yzz_yyz[k];

                g_xyzz_yzz[k] = fr * g1_yzz_yzz[k];

                g_xyzz_zzz[k] = fr * g1_yzz_zzz[k];

                g_xzzz_xxx[k] = fr * g1_zzz_xxx[k] + 3.0 * f2t * gk_zzz_xx[k];

                g_xzzz_xxy[k] = fr * g1_zzz_xxy[k] + 2.0 * f2t * gk_zzz_xy[k];

                g_xzzz_xxz[k] = fr * g1_zzz_xxz[k] + 2.0 * f2t * gk_zzz_xz[k];

                g_xzzz_xyy[k] = fr * g1_zzz_xyy[k] + f2t * gk_zzz_yy[k];

                g_xzzz_xyz[k] = fr * g1_zzz_xyz[k] + f2t * gk_zzz_yz[k];

                g_xzzz_xzz[k] = fr * g1_zzz_xzz[k] + f2t * gk_zzz_zz[k];

                g_xzzz_yyy[k] = fr * g1_zzz_yyy[k];

                g_xzzz_yyz[k] = fr * g1_zzz_yyz[k];

                g_xzzz_yzz[k] = fr * g1_zzz_yzz[k];

                g_xzzz_zzz[k] = fr * g1_zzz_zzz[k];

                // leading y component

                fr = pay[k];

                g_yyyy_xxx[k] = fr * g1_yyy_xxx[k] + f2g * (3.0 * g20_yy_xxx[k] - 3.0 * fgz * g21_yy_xxx[k]);

                g_yyyy_xxy[k] = fr * g1_yyy_xxy[k] + f2g * (3.0 * g20_yy_xxy[k] - 3.0 * fgz * g21_yy_xxy[k]) + f2t * gk_yyy_xx[k];

                g_yyyy_xxz[k] = fr * g1_yyy_xxz[k] + f2g * (3.0 * g20_yy_xxz[k] - 3.0 * fgz * g21_yy_xxz[k]);

                g_yyyy_xyy[k] = fr * g1_yyy_xyy[k] + f2g * (3.0 * g20_yy_xyy[k] - 3.0 * fgz * g21_yy_xyy[k]) + 2.0 * f2t * gk_yyy_xy[k];

                g_yyyy_xyz[k] = fr * g1_yyy_xyz[k] + f2g * (3.0 * g20_yy_xyz[k] - 3.0 * fgz * g21_yy_xyz[k]) + f2t * gk_yyy_xz[k];

                g_yyyy_xzz[k] = fr * g1_yyy_xzz[k] + f2g * (3.0 * g20_yy_xzz[k] - 3.0 * fgz * g21_yy_xzz[k]);

                g_yyyy_yyy[k] = fr * g1_yyy_yyy[k] + f2g * (3.0 * g20_yy_yyy[k] - 3.0 * fgz * g21_yy_yyy[k]) + 3.0 * f2t * gk_yyy_yy[k];

                g_yyyy_yyz[k] = fr * g1_yyy_yyz[k] + f2g * (3.0 * g20_yy_yyz[k] - 3.0 * fgz * g21_yy_yyz[k]) + 2.0 * f2t * gk_yyy_yz[k];

                g_yyyy_yzz[k] = fr * g1_yyy_yzz[k] + f2g * (3.0 * g20_yy_yzz[k] - 3.0 * fgz * g21_yy_yzz[k]) + f2t * gk_yyy_zz[k];

                g_yyyy_zzz[k] = fr * g1_yyy_zzz[k] + f2g * (3.0 * g20_yy_zzz[k] - 3.0 * fgz * g21_yy_zzz[k]);

                g_yyyz_xxx[k] = fr * g1_yyz_xxx[k] + f2g * (2.0 * g20_yz_xxx[k] - 2.0 * fgz * g21_yz_xxx[k]);

                g_yyyz_xxy[k] = fr * g1_yyz_xxy[k] + f2g * (2.0 * g20_yz_xxy[k] - 2.0 * fgz * g21_yz_xxy[k]) + f2t * gk_yyz_xx[k];

                g_yyyz_xxz[k] = fr * g1_yyz_xxz[k] + f2g * (2.0 * g20_yz_xxz[k] - 2.0 * fgz * g21_yz_xxz[k]);

                g_yyyz_xyy[k] = fr * g1_yyz_xyy[k] + f2g * (2.0 * g20_yz_xyy[k] - 2.0 * fgz * g21_yz_xyy[k]) + 2.0 * f2t * gk_yyz_xy[k];

                g_yyyz_xyz[k] = fr * g1_yyz_xyz[k] + f2g * (2.0 * g20_yz_xyz[k] - 2.0 * fgz * g21_yz_xyz[k]) + f2t * gk_yyz_xz[k];

                g_yyyz_xzz[k] = fr * g1_yyz_xzz[k] + f2g * (2.0 * g20_yz_xzz[k] - 2.0 * fgz * g21_yz_xzz[k]);

                g_yyyz_yyy[k] = fr * g1_yyz_yyy[k] + f2g * (2.0 * g20_yz_yyy[k] - 2.0 * fgz * g21_yz_yyy[k]) + 3.0 * f2t * gk_yyz_yy[k];

                g_yyyz_yyz[k] = fr * g1_yyz_yyz[k] + f2g * (2.0 * g20_yz_yyz[k] - 2.0 * fgz * g21_yz_yyz[k]) + 2.0 * f2t * gk_yyz_yz[k];

                g_yyyz_yzz[k] = fr * g1_yyz_yzz[k] + f2g * (2.0 * g20_yz_yzz[k] - 2.0 * fgz * g21_yz_yzz[k]) + f2t * gk_yyz_zz[k];

                g_yyyz_zzz[k] = fr * g1_yyz_zzz[k] + f2g * (2.0 * g20_yz_zzz[k] - 2.0 * fgz * g21_yz_zzz[k]);

                g_yyzz_xxx[k] = fr * g1_yzz_xxx[k] + f2g * (g20_zz_xxx[k] - fgz * g21_zz_xxx[k]);

                g_yyzz_xxy[k] = fr * g1_yzz_xxy[k] + f2g * (g20_zz_xxy[k] - fgz * g21_zz_xxy[k]) + f2t * gk_yzz_xx[k];

                g_yyzz_xxz[k] = fr * g1_yzz_xxz[k] + f2g * (g20_zz_xxz[k] - fgz * g21_zz_xxz[k]);

                g_yyzz_xyy[k] = fr * g1_yzz_xyy[k] + f2g * (g20_zz_xyy[k] - fgz * g21_zz_xyy[k]) + 2.0 * f2t * gk_yzz_xy[k];

                g_yyzz_xyz[k] = fr * g1_yzz_xyz[k] + f2g * (g20_zz_xyz[k] - fgz * g21_zz_xyz[k]) + f2t * gk_yzz_xz[k];

                g_yyzz_xzz[k] = fr * g1_yzz_xzz[k] + f2g * (g20_zz_xzz[k] - fgz * g21_zz_xzz[k]);

                g_yyzz_yyy[k] = fr * g1_yzz_yyy[k] + f2g * (g20_zz_yyy[k] - fgz * g21_zz_yyy[k]) + 3.0 * f2t * gk_yzz_yy[k];

                g_yyzz_yyz[k] = fr * g1_yzz_yyz[k] + f2g * (g20_zz_yyz[k] - fgz * g21_zz_yyz[k]) + 2.0 * f2t * gk_yzz_yz[k];

                g_yyzz_yzz[k] = fr * g1_yzz_yzz[k] + f2g * (g20_zz_yzz[k] - fgz * g21_zz_yzz[k]) + f2t * gk_yzz_zz[k];

                g_yyzz_zzz[k] = fr * g1_yzz_zzz[k] + f2g * (g20_zz_zzz[k] - fgz * g21_zz_zzz[k]);

                g_yzzz_xxx[k] = fr * g1_zzz_xxx[k];

                g_yzzz_xxy[k] = fr * g1_zzz_xxy[k] + f2t * gk_zzz_xx[k];

                g_yzzz_xxz[k] = fr * g1_zzz_xxz[k];

                g_yzzz_xyy[k] = fr * g1_zzz_xyy[k] + 2.0 * f2t * gk_zzz_xy[k];

                g_yzzz_xyz[k] = fr * g1_zzz_xyz[k] + f2t * gk_zzz_xz[k];

                g_yzzz_xzz[k] = fr * g1_zzz_xzz[k];

                g_yzzz_yyy[k] = fr * g1_zzz_yyy[k] + 3.0 * f2t * gk_zzz_yy[k];

                g_yzzz_yyz[k] = fr * g1_zzz_yyz[k] + 2.0 * f2t * gk_zzz_yz[k];

                g_yzzz_yzz[k] = fr * g1_zzz_yzz[k] + f2t * gk_zzz_zz[k];

                g_yzzz_zzz[k] = fr * g1_zzz_zzz[k];

                // leading z component

                fr = paz[k];

                g_zzzz_xxx[k] = fr * g1_zzz_xxx[k] + f2g * (3.0 * g20_zz_xxx[k] - 3.0 * fgz * g21_zz_xxx[k]);

                g_zzzz_xxy[k] = fr * g1_zzz_xxy[k] + f2g * (3.0 * g20_zz_xxy[k] - 3.0 * fgz * g21_zz_xxy[k]);

                g_zzzz_xxz[k] = fr * g1_zzz_xxz[k] + f2g * (3.0 * g20_zz_xxz[k] - 3.0 * fgz * g21_zz_xxz[k]) + f2t * gk_zzz_xx[k];

                g_zzzz_xyy[k] = fr * g1_zzz_xyy[k] + f2g * (3.0 * g20_zz_xyy[k] - 3.0 * fgz * g21_zz_xyy[k]);

                g_zzzz_xyz[k] = fr * g1_zzz_xyz[k] + f2g * (3.0 * g20_zz_xyz[k] - 3.0 * fgz * g21_zz_xyz[k]) + f2t * gk_zzz_xy[k];

                g_zzzz_xzz[k] = fr * g1_zzz_xzz[k] + f2g * (3.0 * g20_zz_xzz[k] - 3.0 * fgz * g21_zz_xzz[k]) + 2.0 * f2t * gk_zzz_xz[k];

                g_zzzz_yyy[k] = fr * g1_zzz_yyy[k] + f2g * (3.0 * g20_zz_yyy[k] - 3.0 * fgz * g21_zz_yyy[k]);

                g_zzzz_yyz[k] = fr * g1_zzz_yyz[k] + f2g * (3.0 * g20_zz_yyz[k] - 3.0 * fgz * g21_zz_yyz[k]) + f2t * gk_zzz_yy[k];

                g_zzzz_yzz[k] = fr * g1_zzz_yzz[k] + f2g * (3.0 * g20_zz_yzz[k] - 3.0 * fgz * g21_zz_yzz[k]) + 2.0 * f2t * gk_zzz_yz[k];

                g_zzzz_zzz[k] = fr * g1_zzz_zzz[k] + f2g * (3.0 * g20_zz_zzz[k] - 3.0 * fgz * g21_zz_zzz[k]) + 3.0 * f2t * gk_zzz_zz[k];
            }

            idx++;
        }
    }
}

void
compElectronicPotentialForGG(CMemBlock2D<double>&        primBuffer,
                             const CVecThreeIndexes&     recPattern,
                             const std::vector<int32_t>& recIndexes,
                             const CMemBlock2D<double>&  osFactors,
                             const CMemBlock2D<double>&  paDistances,
                             const CGtoBlock&            braGtoBlock,
                             const CGtoBlock&            ketGtoBlock,
                             const int32_t               iContrGto)
{
    auto bord = genfunc::maxOrderOfPair(recPattern, 4, 4);

    // skip integrals if not included in recursion pattern

    if (bord < 0) return;

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // compute primitive integrals up to required order

    for (int32_t i = 0; i <= bord; i++)
    {
        // get position of integrals in primitves buffer

        auto goff = genfunc::findTripleIndex(recIndexes, recPattern, {4, 4, i});

        // skip integrals if this order is not required

        if (goff == -1) continue;

        auto g1off = genfunc::findTripleIndex(recIndexes, recPattern, {3, 4, i + 1});

        auto g20off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, i});

        auto g21off = genfunc::findTripleIndex(recIndexes, recPattern, {2, 4, i + 1});

        auto gkoff = genfunc::findTripleIndex(recIndexes, recPattern, {3, 3, i + 1});

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(6 * idx);

            auto fg = osFactors.data(6 * idx + 2);

            auto fz = osFactors.data(6 * idx + 3);

            // set up pointers to distances R(PA)

            auto pax = paDistances.data(3 * idx);

            auto pay = paDistances.data(3 * idx + 1);

            auto paz = paDistances.data(3 * idx + 2);

            // set up pointers to (F|g(r,r')|F)^(m+1) integrals

            auto gk_xxx_xxx = primBuffer.data(gkoff + 100 * idx);

            auto gk_xxx_xxy = primBuffer.data(gkoff + 100 * idx + 1);

            auto gk_xxx_xxz = primBuffer.data(gkoff + 100 * idx + 2);

            auto gk_xxx_xyy = primBuffer.data(gkoff + 100 * idx + 3);

            auto gk_xxx_xyz = primBuffer.data(gkoff + 100 * idx + 4);

            auto gk_xxx_xzz = primBuffer.data(gkoff + 100 * idx + 5);

            auto gk_xxx_yyy = primBuffer.data(gkoff + 100 * idx + 6);

            auto gk_xxx_yyz = primBuffer.data(gkoff + 100 * idx + 7);

            auto gk_xxx_yzz = primBuffer.data(gkoff + 100 * idx + 8);

            auto gk_xxx_zzz = primBuffer.data(gkoff + 100 * idx + 9);

            auto gk_xxy_xxx = primBuffer.data(gkoff + 100 * idx + 10);

            auto gk_xxy_xxy = primBuffer.data(gkoff + 100 * idx + 11);

            auto gk_xxy_xxz = primBuffer.data(gkoff + 100 * idx + 12);

            auto gk_xxy_xyy = primBuffer.data(gkoff + 100 * idx + 13);

            auto gk_xxy_xyz = primBuffer.data(gkoff + 100 * idx + 14);

            auto gk_xxy_xzz = primBuffer.data(gkoff + 100 * idx + 15);

            auto gk_xxy_yyy = primBuffer.data(gkoff + 100 * idx + 16);

            auto gk_xxy_yyz = primBuffer.data(gkoff + 100 * idx + 17);

            auto gk_xxy_yzz = primBuffer.data(gkoff + 100 * idx + 18);

            auto gk_xxy_zzz = primBuffer.data(gkoff + 100 * idx + 19);

            auto gk_xxz_xxx = primBuffer.data(gkoff + 100 * idx + 20);

            auto gk_xxz_xxy = primBuffer.data(gkoff + 100 * idx + 21);

            auto gk_xxz_xxz = primBuffer.data(gkoff + 100 * idx + 22);

            auto gk_xxz_xyy = primBuffer.data(gkoff + 100 * idx + 23);

            auto gk_xxz_xyz = primBuffer.data(gkoff + 100 * idx + 24);

            auto gk_xxz_xzz = primBuffer.data(gkoff + 100 * idx + 25);

            auto gk_xxz_yyy = primBuffer.data(gkoff + 100 * idx + 26);

            auto gk_xxz_yyz = primBuffer.data(gkoff + 100 * idx + 27);

            auto gk_xxz_yzz = primBuffer.data(gkoff + 100 * idx + 28);

            auto gk_xxz_zzz = primBuffer.data(gkoff + 100 * idx + 29);

            auto gk_xyy_xxx = primBuffer.data(gkoff + 100 * idx + 30);

            auto gk_xyy_xxy = primBuffer.data(gkoff + 100 * idx + 31);

            auto gk_xyy_xxz = primBuffer.data(gkoff + 100 * idx + 32);

            auto gk_xyy_xyy = primBuffer.data(gkoff + 100 * idx + 33);

            auto gk_xyy_xyz = primBuffer.data(gkoff + 100 * idx + 34);

            auto gk_xyy_xzz = primBuffer.data(gkoff + 100 * idx + 35);

            auto gk_xyy_yyy = primBuffer.data(gkoff + 100 * idx + 36);

            auto gk_xyy_yyz = primBuffer.data(gkoff + 100 * idx + 37);

            auto gk_xyy_yzz = primBuffer.data(gkoff + 100 * idx + 38);

            auto gk_xyy_zzz = primBuffer.data(gkoff + 100 * idx + 39);

            auto gk_xyz_xxx = primBuffer.data(gkoff + 100 * idx + 40);

            auto gk_xyz_xxy = primBuffer.data(gkoff + 100 * idx + 41);

            auto gk_xyz_xxz = primBuffer.data(gkoff + 100 * idx + 42);

            auto gk_xyz_xyy = primBuffer.data(gkoff + 100 * idx + 43);

            auto gk_xyz_xyz = primBuffer.data(gkoff + 100 * idx + 44);

            auto gk_xyz_xzz = primBuffer.data(gkoff + 100 * idx + 45);

            auto gk_xyz_yyy = primBuffer.data(gkoff + 100 * idx + 46);

            auto gk_xyz_yyz = primBuffer.data(gkoff + 100 * idx + 47);

            auto gk_xyz_yzz = primBuffer.data(gkoff + 100 * idx + 48);

            auto gk_xyz_zzz = primBuffer.data(gkoff + 100 * idx + 49);

            auto gk_xzz_xxx = primBuffer.data(gkoff + 100 * idx + 50);

            auto gk_xzz_xxy = primBuffer.data(gkoff + 100 * idx + 51);

            auto gk_xzz_xxz = primBuffer.data(gkoff + 100 * idx + 52);

            auto gk_xzz_xyy = primBuffer.data(gkoff + 100 * idx + 53);

            auto gk_xzz_xyz = primBuffer.data(gkoff + 100 * idx + 54);

            auto gk_xzz_xzz = primBuffer.data(gkoff + 100 * idx + 55);

            auto gk_xzz_yyy = primBuffer.data(gkoff + 100 * idx + 56);

            auto gk_xzz_yyz = primBuffer.data(gkoff + 100 * idx + 57);

            auto gk_xzz_yzz = primBuffer.data(gkoff + 100 * idx + 58);

            auto gk_xzz_zzz = primBuffer.data(gkoff + 100 * idx + 59);

            auto gk_yyy_xxx = primBuffer.data(gkoff + 100 * idx + 60);

            auto gk_yyy_xxy = primBuffer.data(gkoff + 100 * idx + 61);

            auto gk_yyy_xxz = primBuffer.data(gkoff + 100 * idx + 62);

            auto gk_yyy_xyy = primBuffer.data(gkoff + 100 * idx + 63);

            auto gk_yyy_xyz = primBuffer.data(gkoff + 100 * idx + 64);

            auto gk_yyy_xzz = primBuffer.data(gkoff + 100 * idx + 65);

            auto gk_yyy_yyy = primBuffer.data(gkoff + 100 * idx + 66);

            auto gk_yyy_yyz = primBuffer.data(gkoff + 100 * idx + 67);

            auto gk_yyy_yzz = primBuffer.data(gkoff + 100 * idx + 68);

            auto gk_yyy_zzz = primBuffer.data(gkoff + 100 * idx + 69);

            auto gk_yyz_xxx = primBuffer.data(gkoff + 100 * idx + 70);

            auto gk_yyz_xxy = primBuffer.data(gkoff + 100 * idx + 71);

            auto gk_yyz_xxz = primBuffer.data(gkoff + 100 * idx + 72);

            auto gk_yyz_xyy = primBuffer.data(gkoff + 100 * idx + 73);

            auto gk_yyz_xyz = primBuffer.data(gkoff + 100 * idx + 74);

            auto gk_yyz_xzz = primBuffer.data(gkoff + 100 * idx + 75);

            auto gk_yyz_yyy = primBuffer.data(gkoff + 100 * idx + 76);

            auto gk_yyz_yyz = primBuffer.data(gkoff + 100 * idx + 77);

            auto gk_yyz_yzz = primBuffer.data(gkoff + 100 * idx + 78);

            auto gk_yyz_zzz = primBuffer.data(gkoff + 100 * idx + 79);

            auto gk_yzz_xxx = primBuffer.data(gkoff + 100 * idx + 80);

            auto gk_yzz_xxy = primBuffer.data(gkoff + 100 * idx + 81);

            auto gk_yzz_xxz = primBuffer.data(gkoff + 100 * idx + 82);

            auto gk_yzz_xyy = primBuffer.data(gkoff + 100 * idx + 83);

            auto gk_yzz_xyz = primBuffer.data(gkoff + 100 * idx + 84);

            auto gk_yzz_xzz = primBuffer.data(gkoff + 100 * idx + 85);

            auto gk_yzz_yyy = primBuffer.data(gkoff + 100 * idx + 86);

            auto gk_yzz_yyz = primBuffer.data(gkoff + 100 * idx + 87);

            auto gk_yzz_yzz = primBuffer.data(gkoff + 100 * idx + 88);

            auto gk_yzz_zzz = primBuffer.data(gkoff + 100 * idx + 89);

            auto gk_zzz_xxx = primBuffer.data(gkoff + 100 * idx + 90);

            auto gk_zzz_xxy = primBuffer.data(gkoff + 100 * idx + 91);

            auto gk_zzz_xxz = primBuffer.data(gkoff + 100 * idx + 92);

            auto gk_zzz_xyy = primBuffer.data(gkoff + 100 * idx + 93);

            auto gk_zzz_xyz = primBuffer.data(gkoff + 100 * idx + 94);

            auto gk_zzz_xzz = primBuffer.data(gkoff + 100 * idx + 95);

            auto gk_zzz_yyy = primBuffer.data(gkoff + 100 * idx + 96);

            auto gk_zzz_yyz = primBuffer.data(gkoff + 100 * idx + 97);

            auto gk_zzz_yzz = primBuffer.data(gkoff + 100 * idx + 98);

            auto gk_zzz_zzz = primBuffer.data(gkoff + 100 * idx + 99);

            // set up pointers to (D|g(r,r')|G)^(m) integrals

            auto g20_xx_xxxx = primBuffer.data(g20off + 90 * idx);

            auto g20_xx_xxxy = primBuffer.data(g20off + 90 * idx + 1);

            auto g20_xx_xxxz = primBuffer.data(g20off + 90 * idx + 2);

            auto g20_xx_xxyy = primBuffer.data(g20off + 90 * idx + 3);

            auto g20_xx_xxyz = primBuffer.data(g20off + 90 * idx + 4);

            auto g20_xx_xxzz = primBuffer.data(g20off + 90 * idx + 5);

            auto g20_xx_xyyy = primBuffer.data(g20off + 90 * idx + 6);

            auto g20_xx_xyyz = primBuffer.data(g20off + 90 * idx + 7);

            auto g20_xx_xyzz = primBuffer.data(g20off + 90 * idx + 8);

            auto g20_xx_xzzz = primBuffer.data(g20off + 90 * idx + 9);

            auto g20_xx_yyyy = primBuffer.data(g20off + 90 * idx + 10);

            auto g20_xx_yyyz = primBuffer.data(g20off + 90 * idx + 11);

            auto g20_xx_yyzz = primBuffer.data(g20off + 90 * idx + 12);

            auto g20_xx_yzzz = primBuffer.data(g20off + 90 * idx + 13);

            auto g20_xx_zzzz = primBuffer.data(g20off + 90 * idx + 14);

            auto g20_xy_xxxx = primBuffer.data(g20off + 90 * idx + 15);

            auto g20_xy_xxxy = primBuffer.data(g20off + 90 * idx + 16);

            auto g20_xy_xxxz = primBuffer.data(g20off + 90 * idx + 17);

            auto g20_xy_xxyy = primBuffer.data(g20off + 90 * idx + 18);

            auto g20_xy_xxyz = primBuffer.data(g20off + 90 * idx + 19);

            auto g20_xy_xxzz = primBuffer.data(g20off + 90 * idx + 20);

            auto g20_xy_xyyy = primBuffer.data(g20off + 90 * idx + 21);

            auto g20_xy_xyyz = primBuffer.data(g20off + 90 * idx + 22);

            auto g20_xy_xyzz = primBuffer.data(g20off + 90 * idx + 23);

            auto g20_xy_xzzz = primBuffer.data(g20off + 90 * idx + 24);

            auto g20_xy_yyyy = primBuffer.data(g20off + 90 * idx + 25);

            auto g20_xy_yyyz = primBuffer.data(g20off + 90 * idx + 26);

            auto g20_xy_yyzz = primBuffer.data(g20off + 90 * idx + 27);

            auto g20_xy_yzzz = primBuffer.data(g20off + 90 * idx + 28);

            auto g20_xy_zzzz = primBuffer.data(g20off + 90 * idx + 29);

            auto g20_xz_xxxx = primBuffer.data(g20off + 90 * idx + 30);

            auto g20_xz_xxxy = primBuffer.data(g20off + 90 * idx + 31);

            auto g20_xz_xxxz = primBuffer.data(g20off + 90 * idx + 32);

            auto g20_xz_xxyy = primBuffer.data(g20off + 90 * idx + 33);

            auto g20_xz_xxyz = primBuffer.data(g20off + 90 * idx + 34);

            auto g20_xz_xxzz = primBuffer.data(g20off + 90 * idx + 35);

            auto g20_xz_xyyy = primBuffer.data(g20off + 90 * idx + 36);

            auto g20_xz_xyyz = primBuffer.data(g20off + 90 * idx + 37);

            auto g20_xz_xyzz = primBuffer.data(g20off + 90 * idx + 38);

            auto g20_xz_xzzz = primBuffer.data(g20off + 90 * idx + 39);

            auto g20_xz_yyyy = primBuffer.data(g20off + 90 * idx + 40);

            auto g20_xz_yyyz = primBuffer.data(g20off + 90 * idx + 41);

            auto g20_xz_yyzz = primBuffer.data(g20off + 90 * idx + 42);

            auto g20_xz_yzzz = primBuffer.data(g20off + 90 * idx + 43);

            auto g20_xz_zzzz = primBuffer.data(g20off + 90 * idx + 44);

            auto g20_yy_xxxx = primBuffer.data(g20off + 90 * idx + 45);

            auto g20_yy_xxxy = primBuffer.data(g20off + 90 * idx + 46);

            auto g20_yy_xxxz = primBuffer.data(g20off + 90 * idx + 47);

            auto g20_yy_xxyy = primBuffer.data(g20off + 90 * idx + 48);

            auto g20_yy_xxyz = primBuffer.data(g20off + 90 * idx + 49);

            auto g20_yy_xxzz = primBuffer.data(g20off + 90 * idx + 50);

            auto g20_yy_xyyy = primBuffer.data(g20off + 90 * idx + 51);

            auto g20_yy_xyyz = primBuffer.data(g20off + 90 * idx + 52);

            auto g20_yy_xyzz = primBuffer.data(g20off + 90 * idx + 53);

            auto g20_yy_xzzz = primBuffer.data(g20off + 90 * idx + 54);

            auto g20_yy_yyyy = primBuffer.data(g20off + 90 * idx + 55);

            auto g20_yy_yyyz = primBuffer.data(g20off + 90 * idx + 56);

            auto g20_yy_yyzz = primBuffer.data(g20off + 90 * idx + 57);

            auto g20_yy_yzzz = primBuffer.data(g20off + 90 * idx + 58);

            auto g20_yy_zzzz = primBuffer.data(g20off + 90 * idx + 59);

            auto g20_yz_xxxx = primBuffer.data(g20off + 90 * idx + 60);

            auto g20_yz_xxxy = primBuffer.data(g20off + 90 * idx + 61);

            auto g20_yz_xxxz = primBuffer.data(g20off + 90 * idx + 62);

            auto g20_yz_xxyy = primBuffer.data(g20off + 90 * idx + 63);

            auto g20_yz_xxyz = primBuffer.data(g20off + 90 * idx + 64);

            auto g20_yz_xxzz = primBuffer.data(g20off + 90 * idx + 65);

            auto g20_yz_xyyy = primBuffer.data(g20off + 90 * idx + 66);

            auto g20_yz_xyyz = primBuffer.data(g20off + 90 * idx + 67);

            auto g20_yz_xyzz = primBuffer.data(g20off + 90 * idx + 68);

            auto g20_yz_xzzz = primBuffer.data(g20off + 90 * idx + 69);

            auto g20_yz_yyyy = primBuffer.data(g20off + 90 * idx + 70);

            auto g20_yz_yyyz = primBuffer.data(g20off + 90 * idx + 71);

            auto g20_yz_yyzz = primBuffer.data(g20off + 90 * idx + 72);

            auto g20_yz_yzzz = primBuffer.data(g20off + 90 * idx + 73);

            auto g20_yz_zzzz = primBuffer.data(g20off + 90 * idx + 74);

            auto g20_zz_xxxx = primBuffer.data(g20off + 90 * idx + 75);

            auto g20_zz_xxxy = primBuffer.data(g20off + 90 * idx + 76);

            auto g20_zz_xxxz = primBuffer.data(g20off + 90 * idx + 77);

            auto g20_zz_xxyy = primBuffer.data(g20off + 90 * idx + 78);

            auto g20_zz_xxyz = primBuffer.data(g20off + 90 * idx + 79);

            auto g20_zz_xxzz = primBuffer.data(g20off + 90 * idx + 80);

            auto g20_zz_xyyy = primBuffer.data(g20off + 90 * idx + 81);

            auto g20_zz_xyyz = primBuffer.data(g20off + 90 * idx + 82);

            auto g20_zz_xyzz = primBuffer.data(g20off + 90 * idx + 83);

            auto g20_zz_xzzz = primBuffer.data(g20off + 90 * idx + 84);

            auto g20_zz_yyyy = primBuffer.data(g20off + 90 * idx + 85);

            auto g20_zz_yyyz = primBuffer.data(g20off + 90 * idx + 86);

            auto g20_zz_yyzz = primBuffer.data(g20off + 90 * idx + 87);

            auto g20_zz_yzzz = primBuffer.data(g20off + 90 * idx + 88);

            auto g20_zz_zzzz = primBuffer.data(g20off + 90 * idx + 89);

            // set up pointers to (D|g(r,r')|G)^(m+1) integrals

            auto g21_xx_xxxx = primBuffer.data(g21off + 90 * idx);

            auto g21_xx_xxxy = primBuffer.data(g21off + 90 * idx + 1);

            auto g21_xx_xxxz = primBuffer.data(g21off + 90 * idx + 2);

            auto g21_xx_xxyy = primBuffer.data(g21off + 90 * idx + 3);

            auto g21_xx_xxyz = primBuffer.data(g21off + 90 * idx + 4);

            auto g21_xx_xxzz = primBuffer.data(g21off + 90 * idx + 5);

            auto g21_xx_xyyy = primBuffer.data(g21off + 90 * idx + 6);

            auto g21_xx_xyyz = primBuffer.data(g21off + 90 * idx + 7);

            auto g21_xx_xyzz = primBuffer.data(g21off + 90 * idx + 8);

            auto g21_xx_xzzz = primBuffer.data(g21off + 90 * idx + 9);

            auto g21_xx_yyyy = primBuffer.data(g21off + 90 * idx + 10);

            auto g21_xx_yyyz = primBuffer.data(g21off + 90 * idx + 11);

            auto g21_xx_yyzz = primBuffer.data(g21off + 90 * idx + 12);

            auto g21_xx_yzzz = primBuffer.data(g21off + 90 * idx + 13);

            auto g21_xx_zzzz = primBuffer.data(g21off + 90 * idx + 14);

            auto g21_xy_xxxx = primBuffer.data(g21off + 90 * idx + 15);

            auto g21_xy_xxxy = primBuffer.data(g21off + 90 * idx + 16);

            auto g21_xy_xxxz = primBuffer.data(g21off + 90 * idx + 17);

            auto g21_xy_xxyy = primBuffer.data(g21off + 90 * idx + 18);

            auto g21_xy_xxyz = primBuffer.data(g21off + 90 * idx + 19);

            auto g21_xy_xxzz = primBuffer.data(g21off + 90 * idx + 20);

            auto g21_xy_xyyy = primBuffer.data(g21off + 90 * idx + 21);

            auto g21_xy_xyyz = primBuffer.data(g21off + 90 * idx + 22);

            auto g21_xy_xyzz = primBuffer.data(g21off + 90 * idx + 23);

            auto g21_xy_xzzz = primBuffer.data(g21off + 90 * idx + 24);

            auto g21_xy_yyyy = primBuffer.data(g21off + 90 * idx + 25);

            auto g21_xy_yyyz = primBuffer.data(g21off + 90 * idx + 26);

            auto g21_xy_yyzz = primBuffer.data(g21off + 90 * idx + 27);

            auto g21_xy_yzzz = primBuffer.data(g21off + 90 * idx + 28);

            auto g21_xy_zzzz = primBuffer.data(g21off + 90 * idx + 29);

            auto g21_xz_xxxx = primBuffer.data(g21off + 90 * idx + 30);

            auto g21_xz_xxxy = primBuffer.data(g21off + 90 * idx + 31);

            auto g21_xz_xxxz = primBuffer.data(g21off + 90 * idx + 32);

            auto g21_xz_xxyy = primBuffer.data(g21off + 90 * idx + 33);

            auto g21_xz_xxyz = primBuffer.data(g21off + 90 * idx + 34);

            auto g21_xz_xxzz = primBuffer.data(g21off + 90 * idx + 35);

            auto g21_xz_xyyy = primBuffer.data(g21off + 90 * idx + 36);

            auto g21_xz_xyyz = primBuffer.data(g21off + 90 * idx + 37);

            auto g21_xz_xyzz = primBuffer.data(g21off + 90 * idx + 38);

            auto g21_xz_xzzz = primBuffer.data(g21off + 90 * idx + 39);

            auto g21_xz_yyyy = primBuffer.data(g21off + 90 * idx + 40);

            auto g21_xz_yyyz = primBuffer.data(g21off + 90 * idx + 41);

            auto g21_xz_yyzz = primBuffer.data(g21off + 90 * idx + 42);

            auto g21_xz_yzzz = primBuffer.data(g21off + 90 * idx + 43);

            auto g21_xz_zzzz = primBuffer.data(g21off + 90 * idx + 44);

            auto g21_yy_xxxx = primBuffer.data(g21off + 90 * idx + 45);

            auto g21_yy_xxxy = primBuffer.data(g21off + 90 * idx + 46);

            auto g21_yy_xxxz = primBuffer.data(g21off + 90 * idx + 47);

            auto g21_yy_xxyy = primBuffer.data(g21off + 90 * idx + 48);

            auto g21_yy_xxyz = primBuffer.data(g21off + 90 * idx + 49);

            auto g21_yy_xxzz = primBuffer.data(g21off + 90 * idx + 50);

            auto g21_yy_xyyy = primBuffer.data(g21off + 90 * idx + 51);

            auto g21_yy_xyyz = primBuffer.data(g21off + 90 * idx + 52);

            auto g21_yy_xyzz = primBuffer.data(g21off + 90 * idx + 53);

            auto g21_yy_xzzz = primBuffer.data(g21off + 90 * idx + 54);

            auto g21_yy_yyyy = primBuffer.data(g21off + 90 * idx + 55);

            auto g21_yy_yyyz = primBuffer.data(g21off + 90 * idx + 56);

            auto g21_yy_yyzz = primBuffer.data(g21off + 90 * idx + 57);

            auto g21_yy_yzzz = primBuffer.data(g21off + 90 * idx + 58);

            auto g21_yy_zzzz = primBuffer.data(g21off + 90 * idx + 59);

            auto g21_yz_xxxx = primBuffer.data(g21off + 90 * idx + 60);

            auto g21_yz_xxxy = primBuffer.data(g21off + 90 * idx + 61);

            auto g21_yz_xxxz = primBuffer.data(g21off + 90 * idx + 62);

            auto g21_yz_xxyy = primBuffer.data(g21off + 90 * idx + 63);

            auto g21_yz_xxyz = primBuffer.data(g21off + 90 * idx + 64);

            auto g21_yz_xxzz = primBuffer.data(g21off + 90 * idx + 65);

            auto g21_yz_xyyy = primBuffer.data(g21off + 90 * idx + 66);

            auto g21_yz_xyyz = primBuffer.data(g21off + 90 * idx + 67);

            auto g21_yz_xyzz = primBuffer.data(g21off + 90 * idx + 68);

            auto g21_yz_xzzz = primBuffer.data(g21off + 90 * idx + 69);

            auto g21_yz_yyyy = primBuffer.data(g21off + 90 * idx + 70);

            auto g21_yz_yyyz = primBuffer.data(g21off + 90 * idx + 71);

            auto g21_yz_yyzz = primBuffer.data(g21off + 90 * idx + 72);

            auto g21_yz_yzzz = primBuffer.data(g21off + 90 * idx + 73);

            auto g21_yz_zzzz = primBuffer.data(g21off + 90 * idx + 74);

            auto g21_zz_xxxx = primBuffer.data(g21off + 90 * idx + 75);

            auto g21_zz_xxxy = primBuffer.data(g21off + 90 * idx + 76);

            auto g21_zz_xxxz = primBuffer.data(g21off + 90 * idx + 77);

            auto g21_zz_xxyy = primBuffer.data(g21off + 90 * idx + 78);

            auto g21_zz_xxyz = primBuffer.data(g21off + 90 * idx + 79);

            auto g21_zz_xxzz = primBuffer.data(g21off + 90 * idx + 80);

            auto g21_zz_xyyy = primBuffer.data(g21off + 90 * idx + 81);

            auto g21_zz_xyyz = primBuffer.data(g21off + 90 * idx + 82);

            auto g21_zz_xyzz = primBuffer.data(g21off + 90 * idx + 83);

            auto g21_zz_xzzz = primBuffer.data(g21off + 90 * idx + 84);

            auto g21_zz_yyyy = primBuffer.data(g21off + 90 * idx + 85);

            auto g21_zz_yyyz = primBuffer.data(g21off + 90 * idx + 86);

            auto g21_zz_yyzz = primBuffer.data(g21off + 90 * idx + 87);

            auto g21_zz_yzzz = primBuffer.data(g21off + 90 * idx + 88);

            auto g21_zz_zzzz = primBuffer.data(g21off + 90 * idx + 89);

            // set up pointers to (F|g(r,r')|G)^(m+1) integrals

            auto g1_xxx_xxxx = primBuffer.data(g1off + 150 * idx);

            auto g1_xxx_xxxy = primBuffer.data(g1off + 150 * idx + 1);

            auto g1_xxx_xxxz = primBuffer.data(g1off + 150 * idx + 2);

            auto g1_xxx_xxyy = primBuffer.data(g1off + 150 * idx + 3);

            auto g1_xxx_xxyz = primBuffer.data(g1off + 150 * idx + 4);

            auto g1_xxx_xxzz = primBuffer.data(g1off + 150 * idx + 5);

            auto g1_xxx_xyyy = primBuffer.data(g1off + 150 * idx + 6);

            auto g1_xxx_xyyz = primBuffer.data(g1off + 150 * idx + 7);

            auto g1_xxx_xyzz = primBuffer.data(g1off + 150 * idx + 8);

            auto g1_xxx_xzzz = primBuffer.data(g1off + 150 * idx + 9);

            auto g1_xxx_yyyy = primBuffer.data(g1off + 150 * idx + 10);

            auto g1_xxx_yyyz = primBuffer.data(g1off + 150 * idx + 11);

            auto g1_xxx_yyzz = primBuffer.data(g1off + 150 * idx + 12);

            auto g1_xxx_yzzz = primBuffer.data(g1off + 150 * idx + 13);

            auto g1_xxx_zzzz = primBuffer.data(g1off + 150 * idx + 14);

            auto g1_xxy_xxxx = primBuffer.data(g1off + 150 * idx + 15);

            auto g1_xxy_xxxy = primBuffer.data(g1off + 150 * idx + 16);

            auto g1_xxy_xxxz = primBuffer.data(g1off + 150 * idx + 17);

            auto g1_xxy_xxyy = primBuffer.data(g1off + 150 * idx + 18);

            auto g1_xxy_xxyz = primBuffer.data(g1off + 150 * idx + 19);

            auto g1_xxy_xxzz = primBuffer.data(g1off + 150 * idx + 20);

            auto g1_xxy_xyyy = primBuffer.data(g1off + 150 * idx + 21);

            auto g1_xxy_xyyz = primBuffer.data(g1off + 150 * idx + 22);

            auto g1_xxy_xyzz = primBuffer.data(g1off + 150 * idx + 23);

            auto g1_xxy_xzzz = primBuffer.data(g1off + 150 * idx + 24);

            auto g1_xxy_yyyy = primBuffer.data(g1off + 150 * idx + 25);

            auto g1_xxy_yyyz = primBuffer.data(g1off + 150 * idx + 26);

            auto g1_xxy_yyzz = primBuffer.data(g1off + 150 * idx + 27);

            auto g1_xxy_yzzz = primBuffer.data(g1off + 150 * idx + 28);

            auto g1_xxy_zzzz = primBuffer.data(g1off + 150 * idx + 29);

            auto g1_xxz_xxxx = primBuffer.data(g1off + 150 * idx + 30);

            auto g1_xxz_xxxy = primBuffer.data(g1off + 150 * idx + 31);

            auto g1_xxz_xxxz = primBuffer.data(g1off + 150 * idx + 32);

            auto g1_xxz_xxyy = primBuffer.data(g1off + 150 * idx + 33);

            auto g1_xxz_xxyz = primBuffer.data(g1off + 150 * idx + 34);

            auto g1_xxz_xxzz = primBuffer.data(g1off + 150 * idx + 35);

            auto g1_xxz_xyyy = primBuffer.data(g1off + 150 * idx + 36);

            auto g1_xxz_xyyz = primBuffer.data(g1off + 150 * idx + 37);

            auto g1_xxz_xyzz = primBuffer.data(g1off + 150 * idx + 38);

            auto g1_xxz_xzzz = primBuffer.data(g1off + 150 * idx + 39);

            auto g1_xxz_yyyy = primBuffer.data(g1off + 150 * idx + 40);

            auto g1_xxz_yyyz = primBuffer.data(g1off + 150 * idx + 41);

            auto g1_xxz_yyzz = primBuffer.data(g1off + 150 * idx + 42);

            auto g1_xxz_yzzz = primBuffer.data(g1off + 150 * idx + 43);

            auto g1_xxz_zzzz = primBuffer.data(g1off + 150 * idx + 44);

            auto g1_xyy_xxxx = primBuffer.data(g1off + 150 * idx + 45);

            auto g1_xyy_xxxy = primBuffer.data(g1off + 150 * idx + 46);

            auto g1_xyy_xxxz = primBuffer.data(g1off + 150 * idx + 47);

            auto g1_xyy_xxyy = primBuffer.data(g1off + 150 * idx + 48);

            auto g1_xyy_xxyz = primBuffer.data(g1off + 150 * idx + 49);

            auto g1_xyy_xxzz = primBuffer.data(g1off + 150 * idx + 50);

            auto g1_xyy_xyyy = primBuffer.data(g1off + 150 * idx + 51);

            auto g1_xyy_xyyz = primBuffer.data(g1off + 150 * idx + 52);

            auto g1_xyy_xyzz = primBuffer.data(g1off + 150 * idx + 53);

            auto g1_xyy_xzzz = primBuffer.data(g1off + 150 * idx + 54);

            auto g1_xyy_yyyy = primBuffer.data(g1off + 150 * idx + 55);

            auto g1_xyy_yyyz = primBuffer.data(g1off + 150 * idx + 56);

            auto g1_xyy_yyzz = primBuffer.data(g1off + 150 * idx + 57);

            auto g1_xyy_yzzz = primBuffer.data(g1off + 150 * idx + 58);

            auto g1_xyy_zzzz = primBuffer.data(g1off + 150 * idx + 59);

            auto g1_xyz_xxxx = primBuffer.data(g1off + 150 * idx + 60);

            auto g1_xyz_xxxy = primBuffer.data(g1off + 150 * idx + 61);

            auto g1_xyz_xxxz = primBuffer.data(g1off + 150 * idx + 62);

            auto g1_xyz_xxyy = primBuffer.data(g1off + 150 * idx + 63);

            auto g1_xyz_xxyz = primBuffer.data(g1off + 150 * idx + 64);

            auto g1_xyz_xxzz = primBuffer.data(g1off + 150 * idx + 65);

            auto g1_xyz_xyyy = primBuffer.data(g1off + 150 * idx + 66);

            auto g1_xyz_xyyz = primBuffer.data(g1off + 150 * idx + 67);

            auto g1_xyz_xyzz = primBuffer.data(g1off + 150 * idx + 68);

            auto g1_xyz_xzzz = primBuffer.data(g1off + 150 * idx + 69);

            auto g1_xyz_yyyy = primBuffer.data(g1off + 150 * idx + 70);

            auto g1_xyz_yyyz = primBuffer.data(g1off + 150 * idx + 71);

            auto g1_xyz_yyzz = primBuffer.data(g1off + 150 * idx + 72);

            auto g1_xyz_yzzz = primBuffer.data(g1off + 150 * idx + 73);

            auto g1_xyz_zzzz = primBuffer.data(g1off + 150 * idx + 74);

            auto g1_xzz_xxxx = primBuffer.data(g1off + 150 * idx + 75);

            auto g1_xzz_xxxy = primBuffer.data(g1off + 150 * idx + 76);

            auto g1_xzz_xxxz = primBuffer.data(g1off + 150 * idx + 77);

            auto g1_xzz_xxyy = primBuffer.data(g1off + 150 * idx + 78);

            auto g1_xzz_xxyz = primBuffer.data(g1off + 150 * idx + 79);

            auto g1_xzz_xxzz = primBuffer.data(g1off + 150 * idx + 80);

            auto g1_xzz_xyyy = primBuffer.data(g1off + 150 * idx + 81);

            auto g1_xzz_xyyz = primBuffer.data(g1off + 150 * idx + 82);

            auto g1_xzz_xyzz = primBuffer.data(g1off + 150 * idx + 83);

            auto g1_xzz_xzzz = primBuffer.data(g1off + 150 * idx + 84);

            auto g1_xzz_yyyy = primBuffer.data(g1off + 150 * idx + 85);

            auto g1_xzz_yyyz = primBuffer.data(g1off + 150 * idx + 86);

            auto g1_xzz_yyzz = primBuffer.data(g1off + 150 * idx + 87);

            auto g1_xzz_yzzz = primBuffer.data(g1off + 150 * idx + 88);

            auto g1_xzz_zzzz = primBuffer.data(g1off + 150 * idx + 89);

            auto g1_yyy_xxxx = primBuffer.data(g1off + 150 * idx + 90);

            auto g1_yyy_xxxy = primBuffer.data(g1off + 150 * idx + 91);

            auto g1_yyy_xxxz = primBuffer.data(g1off + 150 * idx + 92);

            auto g1_yyy_xxyy = primBuffer.data(g1off + 150 * idx + 93);

            auto g1_yyy_xxyz = primBuffer.data(g1off + 150 * idx + 94);

            auto g1_yyy_xxzz = primBuffer.data(g1off + 150 * idx + 95);

            auto g1_yyy_xyyy = primBuffer.data(g1off + 150 * idx + 96);

            auto g1_yyy_xyyz = primBuffer.data(g1off + 150 * idx + 97);

            auto g1_yyy_xyzz = primBuffer.data(g1off + 150 * idx + 98);

            auto g1_yyy_xzzz = primBuffer.data(g1off + 150 * idx + 99);

            auto g1_yyy_yyyy = primBuffer.data(g1off + 150 * idx + 100);

            auto g1_yyy_yyyz = primBuffer.data(g1off + 150 * idx + 101);

            auto g1_yyy_yyzz = primBuffer.data(g1off + 150 * idx + 102);

            auto g1_yyy_yzzz = primBuffer.data(g1off + 150 * idx + 103);

            auto g1_yyy_zzzz = primBuffer.data(g1off + 150 * idx + 104);

            auto g1_yyz_xxxx = primBuffer.data(g1off + 150 * idx + 105);

            auto g1_yyz_xxxy = primBuffer.data(g1off + 150 * idx + 106);

            auto g1_yyz_xxxz = primBuffer.data(g1off + 150 * idx + 107);

            auto g1_yyz_xxyy = primBuffer.data(g1off + 150 * idx + 108);

            auto g1_yyz_xxyz = primBuffer.data(g1off + 150 * idx + 109);

            auto g1_yyz_xxzz = primBuffer.data(g1off + 150 * idx + 110);

            auto g1_yyz_xyyy = primBuffer.data(g1off + 150 * idx + 111);

            auto g1_yyz_xyyz = primBuffer.data(g1off + 150 * idx + 112);

            auto g1_yyz_xyzz = primBuffer.data(g1off + 150 * idx + 113);

            auto g1_yyz_xzzz = primBuffer.data(g1off + 150 * idx + 114);

            auto g1_yyz_yyyy = primBuffer.data(g1off + 150 * idx + 115);

            auto g1_yyz_yyyz = primBuffer.data(g1off + 150 * idx + 116);

            auto g1_yyz_yyzz = primBuffer.data(g1off + 150 * idx + 117);

            auto g1_yyz_yzzz = primBuffer.data(g1off + 150 * idx + 118);

            auto g1_yyz_zzzz = primBuffer.data(g1off + 150 * idx + 119);

            auto g1_yzz_xxxx = primBuffer.data(g1off + 150 * idx + 120);

            auto g1_yzz_xxxy = primBuffer.data(g1off + 150 * idx + 121);

            auto g1_yzz_xxxz = primBuffer.data(g1off + 150 * idx + 122);

            auto g1_yzz_xxyy = primBuffer.data(g1off + 150 * idx + 123);

            auto g1_yzz_xxyz = primBuffer.data(g1off + 150 * idx + 124);

            auto g1_yzz_xxzz = primBuffer.data(g1off + 150 * idx + 125);

            auto g1_yzz_xyyy = primBuffer.data(g1off + 150 * idx + 126);

            auto g1_yzz_xyyz = primBuffer.data(g1off + 150 * idx + 127);

            auto g1_yzz_xyzz = primBuffer.data(g1off + 150 * idx + 128);

            auto g1_yzz_xzzz = primBuffer.data(g1off + 150 * idx + 129);

            auto g1_yzz_yyyy = primBuffer.data(g1off + 150 * idx + 130);

            auto g1_yzz_yyyz = primBuffer.data(g1off + 150 * idx + 131);

            auto g1_yzz_yyzz = primBuffer.data(g1off + 150 * idx + 132);

            auto g1_yzz_yzzz = primBuffer.data(g1off + 150 * idx + 133);

            auto g1_yzz_zzzz = primBuffer.data(g1off + 150 * idx + 134);

            auto g1_zzz_xxxx = primBuffer.data(g1off + 150 * idx + 135);

            auto g1_zzz_xxxy = primBuffer.data(g1off + 150 * idx + 136);

            auto g1_zzz_xxxz = primBuffer.data(g1off + 150 * idx + 137);

            auto g1_zzz_xxyy = primBuffer.data(g1off + 150 * idx + 138);

            auto g1_zzz_xxyz = primBuffer.data(g1off + 150 * idx + 139);

            auto g1_zzz_xxzz = primBuffer.data(g1off + 150 * idx + 140);

            auto g1_zzz_xyyy = primBuffer.data(g1off + 150 * idx + 141);

            auto g1_zzz_xyyz = primBuffer.data(g1off + 150 * idx + 142);

            auto g1_zzz_xyzz = primBuffer.data(g1off + 150 * idx + 143);

            auto g1_zzz_xzzz = primBuffer.data(g1off + 150 * idx + 144);

            auto g1_zzz_yyyy = primBuffer.data(g1off + 150 * idx + 145);

            auto g1_zzz_yyyz = primBuffer.data(g1off + 150 * idx + 146);

            auto g1_zzz_yyzz = primBuffer.data(g1off + 150 * idx + 147);

            auto g1_zzz_yzzz = primBuffer.data(g1off + 150 * idx + 148);

            auto g1_zzz_zzzz = primBuffer.data(g1off + 150 * idx + 149);

            // set up pointers to (G|g(r,r')|G)^(m) integrals

            auto g_xxxx_xxxx = primBuffer.data(goff + 225 * idx);

            auto g_xxxx_xxxy = primBuffer.data(goff + 225 * idx + 1);

            auto g_xxxx_xxxz = primBuffer.data(goff + 225 * idx + 2);

            auto g_xxxx_xxyy = primBuffer.data(goff + 225 * idx + 3);

            auto g_xxxx_xxyz = primBuffer.data(goff + 225 * idx + 4);

            auto g_xxxx_xxzz = primBuffer.data(goff + 225 * idx + 5);

            auto g_xxxx_xyyy = primBuffer.data(goff + 225 * idx + 6);

            auto g_xxxx_xyyz = primBuffer.data(goff + 225 * idx + 7);

            auto g_xxxx_xyzz = primBuffer.data(goff + 225 * idx + 8);

            auto g_xxxx_xzzz = primBuffer.data(goff + 225 * idx + 9);

            auto g_xxxx_yyyy = primBuffer.data(goff + 225 * idx + 10);

            auto g_xxxx_yyyz = primBuffer.data(goff + 225 * idx + 11);

            auto g_xxxx_yyzz = primBuffer.data(goff + 225 * idx + 12);

            auto g_xxxx_yzzz = primBuffer.data(goff + 225 * idx + 13);

            auto g_xxxx_zzzz = primBuffer.data(goff + 225 * idx + 14);

            auto g_xxxy_xxxx = primBuffer.data(goff + 225 * idx + 15);

            auto g_xxxy_xxxy = primBuffer.data(goff + 225 * idx + 16);

            auto g_xxxy_xxxz = primBuffer.data(goff + 225 * idx + 17);

            auto g_xxxy_xxyy = primBuffer.data(goff + 225 * idx + 18);

            auto g_xxxy_xxyz = primBuffer.data(goff + 225 * idx + 19);

            auto g_xxxy_xxzz = primBuffer.data(goff + 225 * idx + 20);

            auto g_xxxy_xyyy = primBuffer.data(goff + 225 * idx + 21);

            auto g_xxxy_xyyz = primBuffer.data(goff + 225 * idx + 22);

            auto g_xxxy_xyzz = primBuffer.data(goff + 225 * idx + 23);

            auto g_xxxy_xzzz = primBuffer.data(goff + 225 * idx + 24);

            auto g_xxxy_yyyy = primBuffer.data(goff + 225 * idx + 25);

            auto g_xxxy_yyyz = primBuffer.data(goff + 225 * idx + 26);

            auto g_xxxy_yyzz = primBuffer.data(goff + 225 * idx + 27);

            auto g_xxxy_yzzz = primBuffer.data(goff + 225 * idx + 28);

            auto g_xxxy_zzzz = primBuffer.data(goff + 225 * idx + 29);

            auto g_xxxz_xxxx = primBuffer.data(goff + 225 * idx + 30);

            auto g_xxxz_xxxy = primBuffer.data(goff + 225 * idx + 31);

            auto g_xxxz_xxxz = primBuffer.data(goff + 225 * idx + 32);

            auto g_xxxz_xxyy = primBuffer.data(goff + 225 * idx + 33);

            auto g_xxxz_xxyz = primBuffer.data(goff + 225 * idx + 34);

            auto g_xxxz_xxzz = primBuffer.data(goff + 225 * idx + 35);

            auto g_xxxz_xyyy = primBuffer.data(goff + 225 * idx + 36);

            auto g_xxxz_xyyz = primBuffer.data(goff + 225 * idx + 37);

            auto g_xxxz_xyzz = primBuffer.data(goff + 225 * idx + 38);

            auto g_xxxz_xzzz = primBuffer.data(goff + 225 * idx + 39);

            auto g_xxxz_yyyy = primBuffer.data(goff + 225 * idx + 40);

            auto g_xxxz_yyyz = primBuffer.data(goff + 225 * idx + 41);

            auto g_xxxz_yyzz = primBuffer.data(goff + 225 * idx + 42);

            auto g_xxxz_yzzz = primBuffer.data(goff + 225 * idx + 43);

            auto g_xxxz_zzzz = primBuffer.data(goff + 225 * idx + 44);

            auto g_xxyy_xxxx = primBuffer.data(goff + 225 * idx + 45);

            auto g_xxyy_xxxy = primBuffer.data(goff + 225 * idx + 46);

            auto g_xxyy_xxxz = primBuffer.data(goff + 225 * idx + 47);

            auto g_xxyy_xxyy = primBuffer.data(goff + 225 * idx + 48);

            auto g_xxyy_xxyz = primBuffer.data(goff + 225 * idx + 49);

            auto g_xxyy_xxzz = primBuffer.data(goff + 225 * idx + 50);

            auto g_xxyy_xyyy = primBuffer.data(goff + 225 * idx + 51);

            auto g_xxyy_xyyz = primBuffer.data(goff + 225 * idx + 52);

            auto g_xxyy_xyzz = primBuffer.data(goff + 225 * idx + 53);

            auto g_xxyy_xzzz = primBuffer.data(goff + 225 * idx + 54);

            auto g_xxyy_yyyy = primBuffer.data(goff + 225 * idx + 55);

            auto g_xxyy_yyyz = primBuffer.data(goff + 225 * idx + 56);

            auto g_xxyy_yyzz = primBuffer.data(goff + 225 * idx + 57);

            auto g_xxyy_yzzz = primBuffer.data(goff + 225 * idx + 58);

            auto g_xxyy_zzzz = primBuffer.data(goff + 225 * idx + 59);

            auto g_xxyz_xxxx = primBuffer.data(goff + 225 * idx + 60);

            auto g_xxyz_xxxy = primBuffer.data(goff + 225 * idx + 61);

            auto g_xxyz_xxxz = primBuffer.data(goff + 225 * idx + 62);

            auto g_xxyz_xxyy = primBuffer.data(goff + 225 * idx + 63);

            auto g_xxyz_xxyz = primBuffer.data(goff + 225 * idx + 64);

            auto g_xxyz_xxzz = primBuffer.data(goff + 225 * idx + 65);

            auto g_xxyz_xyyy = primBuffer.data(goff + 225 * idx + 66);

            auto g_xxyz_xyyz = primBuffer.data(goff + 225 * idx + 67);

            auto g_xxyz_xyzz = primBuffer.data(goff + 225 * idx + 68);

            auto g_xxyz_xzzz = primBuffer.data(goff + 225 * idx + 69);

            auto g_xxyz_yyyy = primBuffer.data(goff + 225 * idx + 70);

            auto g_xxyz_yyyz = primBuffer.data(goff + 225 * idx + 71);

            auto g_xxyz_yyzz = primBuffer.data(goff + 225 * idx + 72);

            auto g_xxyz_yzzz = primBuffer.data(goff + 225 * idx + 73);

            auto g_xxyz_zzzz = primBuffer.data(goff + 225 * idx + 74);

            auto g_xxzz_xxxx = primBuffer.data(goff + 225 * idx + 75);

            auto g_xxzz_xxxy = primBuffer.data(goff + 225 * idx + 76);

            auto g_xxzz_xxxz = primBuffer.data(goff + 225 * idx + 77);

            auto g_xxzz_xxyy = primBuffer.data(goff + 225 * idx + 78);

            auto g_xxzz_xxyz = primBuffer.data(goff + 225 * idx + 79);

            auto g_xxzz_xxzz = primBuffer.data(goff + 225 * idx + 80);

            auto g_xxzz_xyyy = primBuffer.data(goff + 225 * idx + 81);

            auto g_xxzz_xyyz = primBuffer.data(goff + 225 * idx + 82);

            auto g_xxzz_xyzz = primBuffer.data(goff + 225 * idx + 83);

            auto g_xxzz_xzzz = primBuffer.data(goff + 225 * idx + 84);

            auto g_xxzz_yyyy = primBuffer.data(goff + 225 * idx + 85);

            auto g_xxzz_yyyz = primBuffer.data(goff + 225 * idx + 86);

            auto g_xxzz_yyzz = primBuffer.data(goff + 225 * idx + 87);

            auto g_xxzz_yzzz = primBuffer.data(goff + 225 * idx + 88);

            auto g_xxzz_zzzz = primBuffer.data(goff + 225 * idx + 89);

            auto g_xyyy_xxxx = primBuffer.data(goff + 225 * idx + 90);

            auto g_xyyy_xxxy = primBuffer.data(goff + 225 * idx + 91);

            auto g_xyyy_xxxz = primBuffer.data(goff + 225 * idx + 92);

            auto g_xyyy_xxyy = primBuffer.data(goff + 225 * idx + 93);

            auto g_xyyy_xxyz = primBuffer.data(goff + 225 * idx + 94);

            auto g_xyyy_xxzz = primBuffer.data(goff + 225 * idx + 95);

            auto g_xyyy_xyyy = primBuffer.data(goff + 225 * idx + 96);

            auto g_xyyy_xyyz = primBuffer.data(goff + 225 * idx + 97);

            auto g_xyyy_xyzz = primBuffer.data(goff + 225 * idx + 98);

            auto g_xyyy_xzzz = primBuffer.data(goff + 225 * idx + 99);

            auto g_xyyy_yyyy = primBuffer.data(goff + 225 * idx + 100);

            auto g_xyyy_yyyz = primBuffer.data(goff + 225 * idx + 101);

            auto g_xyyy_yyzz = primBuffer.data(goff + 225 * idx + 102);

            auto g_xyyy_yzzz = primBuffer.data(goff + 225 * idx + 103);

            auto g_xyyy_zzzz = primBuffer.data(goff + 225 * idx + 104);

            auto g_xyyz_xxxx = primBuffer.data(goff + 225 * idx + 105);

            auto g_xyyz_xxxy = primBuffer.data(goff + 225 * idx + 106);

            auto g_xyyz_xxxz = primBuffer.data(goff + 225 * idx + 107);

            auto g_xyyz_xxyy = primBuffer.data(goff + 225 * idx + 108);

            auto g_xyyz_xxyz = primBuffer.data(goff + 225 * idx + 109);

            auto g_xyyz_xxzz = primBuffer.data(goff + 225 * idx + 110);

            auto g_xyyz_xyyy = primBuffer.data(goff + 225 * idx + 111);

            auto g_xyyz_xyyz = primBuffer.data(goff + 225 * idx + 112);

            auto g_xyyz_xyzz = primBuffer.data(goff + 225 * idx + 113);

            auto g_xyyz_xzzz = primBuffer.data(goff + 225 * idx + 114);

            auto g_xyyz_yyyy = primBuffer.data(goff + 225 * idx + 115);

            auto g_xyyz_yyyz = primBuffer.data(goff + 225 * idx + 116);

            auto g_xyyz_yyzz = primBuffer.data(goff + 225 * idx + 117);

            auto g_xyyz_yzzz = primBuffer.data(goff + 225 * idx + 118);

            auto g_xyyz_zzzz = primBuffer.data(goff + 225 * idx + 119);

            auto g_xyzz_xxxx = primBuffer.data(goff + 225 * idx + 120);

            auto g_xyzz_xxxy = primBuffer.data(goff + 225 * idx + 121);

            auto g_xyzz_xxxz = primBuffer.data(goff + 225 * idx + 122);

            auto g_xyzz_xxyy = primBuffer.data(goff + 225 * idx + 123);

            auto g_xyzz_xxyz = primBuffer.data(goff + 225 * idx + 124);

            auto g_xyzz_xxzz = primBuffer.data(goff + 225 * idx + 125);

            auto g_xyzz_xyyy = primBuffer.data(goff + 225 * idx + 126);

            auto g_xyzz_xyyz = primBuffer.data(goff + 225 * idx + 127);

            auto g_xyzz_xyzz = primBuffer.data(goff + 225 * idx + 128);

            auto g_xyzz_xzzz = primBuffer.data(goff + 225 * idx + 129);

            auto g_xyzz_yyyy = primBuffer.data(goff + 225 * idx + 130);

            auto g_xyzz_yyyz = primBuffer.data(goff + 225 * idx + 131);

            auto g_xyzz_yyzz = primBuffer.data(goff + 225 * idx + 132);

            auto g_xyzz_yzzz = primBuffer.data(goff + 225 * idx + 133);

            auto g_xyzz_zzzz = primBuffer.data(goff + 225 * idx + 134);

            auto g_xzzz_xxxx = primBuffer.data(goff + 225 * idx + 135);

            auto g_xzzz_xxxy = primBuffer.data(goff + 225 * idx + 136);

            auto g_xzzz_xxxz = primBuffer.data(goff + 225 * idx + 137);

            auto g_xzzz_xxyy = primBuffer.data(goff + 225 * idx + 138);

            auto g_xzzz_xxyz = primBuffer.data(goff + 225 * idx + 139);

            auto g_xzzz_xxzz = primBuffer.data(goff + 225 * idx + 140);

            auto g_xzzz_xyyy = primBuffer.data(goff + 225 * idx + 141);

            auto g_xzzz_xyyz = primBuffer.data(goff + 225 * idx + 142);

            auto g_xzzz_xyzz = primBuffer.data(goff + 225 * idx + 143);

            auto g_xzzz_xzzz = primBuffer.data(goff + 225 * idx + 144);

            auto g_xzzz_yyyy = primBuffer.data(goff + 225 * idx + 145);

            auto g_xzzz_yyyz = primBuffer.data(goff + 225 * idx + 146);

            auto g_xzzz_yyzz = primBuffer.data(goff + 225 * idx + 147);

            auto g_xzzz_yzzz = primBuffer.data(goff + 225 * idx + 148);

            auto g_xzzz_zzzz = primBuffer.data(goff + 225 * idx + 149);

            auto g_yyyy_xxxx = primBuffer.data(goff + 225 * idx + 150);

            auto g_yyyy_xxxy = primBuffer.data(goff + 225 * idx + 151);

            auto g_yyyy_xxxz = primBuffer.data(goff + 225 * idx + 152);

            auto g_yyyy_xxyy = primBuffer.data(goff + 225 * idx + 153);

            auto g_yyyy_xxyz = primBuffer.data(goff + 225 * idx + 154);

            auto g_yyyy_xxzz = primBuffer.data(goff + 225 * idx + 155);

            auto g_yyyy_xyyy = primBuffer.data(goff + 225 * idx + 156);

            auto g_yyyy_xyyz = primBuffer.data(goff + 225 * idx + 157);

            auto g_yyyy_xyzz = primBuffer.data(goff + 225 * idx + 158);

            auto g_yyyy_xzzz = primBuffer.data(goff + 225 * idx + 159);

            auto g_yyyy_yyyy = primBuffer.data(goff + 225 * idx + 160);

            auto g_yyyy_yyyz = primBuffer.data(goff + 225 * idx + 161);

            auto g_yyyy_yyzz = primBuffer.data(goff + 225 * idx + 162);

            auto g_yyyy_yzzz = primBuffer.data(goff + 225 * idx + 163);

            auto g_yyyy_zzzz = primBuffer.data(goff + 225 * idx + 164);

            auto g_yyyz_xxxx = primBuffer.data(goff + 225 * idx + 165);

            auto g_yyyz_xxxy = primBuffer.data(goff + 225 * idx + 166);

            auto g_yyyz_xxxz = primBuffer.data(goff + 225 * idx + 167);

            auto g_yyyz_xxyy = primBuffer.data(goff + 225 * idx + 168);

            auto g_yyyz_xxyz = primBuffer.data(goff + 225 * idx + 169);

            auto g_yyyz_xxzz = primBuffer.data(goff + 225 * idx + 170);

            auto g_yyyz_xyyy = primBuffer.data(goff + 225 * idx + 171);

            auto g_yyyz_xyyz = primBuffer.data(goff + 225 * idx + 172);

            auto g_yyyz_xyzz = primBuffer.data(goff + 225 * idx + 173);

            auto g_yyyz_xzzz = primBuffer.data(goff + 225 * idx + 174);

            auto g_yyyz_yyyy = primBuffer.data(goff + 225 * idx + 175);

            auto g_yyyz_yyyz = primBuffer.data(goff + 225 * idx + 176);

            auto g_yyyz_yyzz = primBuffer.data(goff + 225 * idx + 177);

            auto g_yyyz_yzzz = primBuffer.data(goff + 225 * idx + 178);

            auto g_yyyz_zzzz = primBuffer.data(goff + 225 * idx + 179);

            auto g_yyzz_xxxx = primBuffer.data(goff + 225 * idx + 180);

            auto g_yyzz_xxxy = primBuffer.data(goff + 225 * idx + 181);

            auto g_yyzz_xxxz = primBuffer.data(goff + 225 * idx + 182);

            auto g_yyzz_xxyy = primBuffer.data(goff + 225 * idx + 183);

            auto g_yyzz_xxyz = primBuffer.data(goff + 225 * idx + 184);

            auto g_yyzz_xxzz = primBuffer.data(goff + 225 * idx + 185);

            auto g_yyzz_xyyy = primBuffer.data(goff + 225 * idx + 186);

            auto g_yyzz_xyyz = primBuffer.data(goff + 225 * idx + 187);

            auto g_yyzz_xyzz = primBuffer.data(goff + 225 * idx + 188);

            auto g_yyzz_xzzz = primBuffer.data(goff + 225 * idx + 189);

            auto g_yyzz_yyyy = primBuffer.data(goff + 225 * idx + 190);

            auto g_yyzz_yyyz = primBuffer.data(goff + 225 * idx + 191);

            auto g_yyzz_yyzz = primBuffer.data(goff + 225 * idx + 192);

            auto g_yyzz_yzzz = primBuffer.data(goff + 225 * idx + 193);

            auto g_yyzz_zzzz = primBuffer.data(goff + 225 * idx + 194);

            auto g_yzzz_xxxx = primBuffer.data(goff + 225 * idx + 195);

            auto g_yzzz_xxxy = primBuffer.data(goff + 225 * idx + 196);

            auto g_yzzz_xxxz = primBuffer.data(goff + 225 * idx + 197);

            auto g_yzzz_xxyy = primBuffer.data(goff + 225 * idx + 198);

            auto g_yzzz_xxyz = primBuffer.data(goff + 225 * idx + 199);

            auto g_yzzz_xxzz = primBuffer.data(goff + 225 * idx + 200);

            auto g_yzzz_xyyy = primBuffer.data(goff + 225 * idx + 201);

            auto g_yzzz_xyyz = primBuffer.data(goff + 225 * idx + 202);

            auto g_yzzz_xyzz = primBuffer.data(goff + 225 * idx + 203);

            auto g_yzzz_xzzz = primBuffer.data(goff + 225 * idx + 204);

            auto g_yzzz_yyyy = primBuffer.data(goff + 225 * idx + 205);

            auto g_yzzz_yyyz = primBuffer.data(goff + 225 * idx + 206);

            auto g_yzzz_yyzz = primBuffer.data(goff + 225 * idx + 207);

            auto g_yzzz_yzzz = primBuffer.data(goff + 225 * idx + 208);

            auto g_yzzz_zzzz = primBuffer.data(goff + 225 * idx + 209);

            auto g_zzzz_xxxx = primBuffer.data(goff + 225 * idx + 210);

            auto g_zzzz_xxxy = primBuffer.data(goff + 225 * idx + 211);

            auto g_zzzz_xxxz = primBuffer.data(goff + 225 * idx + 212);

            auto g_zzzz_xxyy = primBuffer.data(goff + 225 * idx + 213);

            auto g_zzzz_xxyz = primBuffer.data(goff + 225 * idx + 214);

            auto g_zzzz_xxzz = primBuffer.data(goff + 225 * idx + 215);

            auto g_zzzz_xyyy = primBuffer.data(goff + 225 * idx + 216);

            auto g_zzzz_xyyz = primBuffer.data(goff + 225 * idx + 217);

            auto g_zzzz_xyzz = primBuffer.data(goff + 225 * idx + 218);

            auto g_zzzz_xzzz = primBuffer.data(goff + 225 * idx + 219);

            auto g_zzzz_yyyy = primBuffer.data(goff + 225 * idx + 220);

            auto g_zzzz_yyyz = primBuffer.data(goff + 225 * idx + 221);

            auto g_zzzz_yyzz = primBuffer.data(goff + 225 * idx + 222);

            auto g_zzzz_yzzz = primBuffer.data(goff + 225 * idx + 223);

            auto g_zzzz_zzzz = primBuffer.data(goff + 225 * idx + 224);

            #pragma omp simd aligned(fx, fg, fz, pax, pay, paz, gk_xxx_xxx,\
                                         gk_xxx_xxy, gk_xxx_xxz, gk_xxx_xyy, gk_xxx_xyz,\
                                         gk_xxx_xzz, gk_xxx_yyy, gk_xxx_yyz, gk_xxx_yzz,\
                                         gk_xxx_zzz, gk_xxy_xxx, gk_xxy_xxy, gk_xxy_xxz,\
                                         gk_xxy_xyy, gk_xxy_xyz, gk_xxy_xzz, gk_xxy_yyy,\
                                         gk_xxy_yyz, gk_xxy_yzz, gk_xxy_zzz, gk_xxz_xxx,\
                                         gk_xxz_xxy, gk_xxz_xxz, gk_xxz_xyy, gk_xxz_xyz,\
                                         gk_xxz_xzz, gk_xxz_yyy, gk_xxz_yyz, gk_xxz_yzz,\
                                         gk_xxz_zzz, gk_xyy_xxx, gk_xyy_xxy, gk_xyy_xxz,\
                                         gk_xyy_xyy, gk_xyy_xyz, gk_xyy_xzz, gk_xyy_yyy,\
                                         gk_xyy_yyz, gk_xyy_yzz, gk_xyy_zzz, gk_xyz_xxx,\
                                         gk_xyz_xxy, gk_xyz_xxz, gk_xyz_xyy, gk_xyz_xyz,\
                                         gk_xyz_xzz, gk_xyz_yyy, gk_xyz_yyz, gk_xyz_yzz,\
                                         gk_xyz_zzz, gk_xzz_xxx, gk_xzz_xxy, gk_xzz_xxz,\
                                         gk_xzz_xyy, gk_xzz_xyz, gk_xzz_xzz, gk_xzz_yyy,\
                                         gk_xzz_yyz, gk_xzz_yzz, gk_xzz_zzz, gk_yyy_xxx,\
                                         gk_yyy_xxy, gk_yyy_xxz, gk_yyy_xyy, gk_yyy_xyz,\
                                         gk_yyy_xzz, gk_yyy_yyy, gk_yyy_yyz, gk_yyy_yzz,\
                                         gk_yyy_zzz, gk_yyz_xxx, gk_yyz_xxy, gk_yyz_xxz,\
                                         gk_yyz_xyy, gk_yyz_xyz, gk_yyz_xzz, gk_yyz_yyy,\
                                         gk_yyz_yyz, gk_yyz_yzz, gk_yyz_zzz, gk_yzz_xxx,\
                                         gk_yzz_xxy, gk_yzz_xxz, gk_yzz_xyy, gk_yzz_xyz,\
                                         gk_yzz_xzz, gk_yzz_yyy, gk_yzz_yyz, gk_yzz_yzz,\
                                         gk_yzz_zzz, gk_zzz_xxx, gk_zzz_xxy, gk_zzz_xxz,\
                                         gk_zzz_xyy, gk_zzz_xyz, gk_zzz_xzz, gk_zzz_yyy,\
                                         gk_zzz_yyz, gk_zzz_yzz, gk_zzz_zzz, g20_xx_xxxx,\
                                         g20_xx_xxxy, g20_xx_xxxz, g20_xx_xxyy,\
                                         g20_xx_xxyz, g20_xx_xxzz, g20_xx_xyyy,\
                                         g20_xx_xyyz, g20_xx_xyzz, g20_xx_xzzz,\
                                         g20_xx_yyyy, g20_xx_yyyz, g20_xx_yyzz,\
                                         g20_xx_yzzz, g20_xx_zzzz, g20_xy_xxxx,\
                                         g20_xy_xxxy, g20_xy_xxxz, g20_xy_xxyy,\
                                         g20_xy_xxyz, g20_xy_xxzz, g20_xy_xyyy,\
                                         g20_xy_xyyz, g20_xy_xyzz, g20_xy_xzzz,\
                                         g20_xy_yyyy, g20_xy_yyyz, g20_xy_yyzz,\
                                         g20_xy_yzzz, g20_xy_zzzz, g20_xz_xxxx,\
                                         g20_xz_xxxy, g20_xz_xxxz, g20_xz_xxyy,\
                                         g20_xz_xxyz, g20_xz_xxzz, g20_xz_xyyy,\
                                         g20_xz_xyyz, g20_xz_xyzz, g20_xz_xzzz,\
                                         g20_xz_yyyy, g20_xz_yyyz, g20_xz_yyzz,\
                                         g20_xz_yzzz, g20_xz_zzzz, g20_yy_xxxx,\
                                         g20_yy_xxxy, g20_yy_xxxz, g20_yy_xxyy,\
                                         g20_yy_xxyz, g20_yy_xxzz, g20_yy_xyyy,\
                                         g20_yy_xyyz, g20_yy_xyzz, g20_yy_xzzz,\
                                         g20_yy_yyyy, g20_yy_yyyz, g20_yy_yyzz,\
                                         g20_yy_yzzz, g20_yy_zzzz, g20_yz_xxxx,\
                                         g20_yz_xxxy, g20_yz_xxxz, g20_yz_xxyy,\
                                         g20_yz_xxyz, g20_yz_xxzz, g20_yz_xyyy,\
                                         g20_yz_xyyz, g20_yz_xyzz, g20_yz_xzzz,\
                                         g20_yz_yyyy, g20_yz_yyyz, g20_yz_yyzz,\
                                         g20_yz_yzzz, g20_yz_zzzz, g20_zz_xxxx,\
                                         g20_zz_xxxy, g20_zz_xxxz, g20_zz_xxyy,\
                                         g20_zz_xxyz, g20_zz_xxzz, g20_zz_xyyy,\
                                         g20_zz_xyyz, g20_zz_xyzz, g20_zz_xzzz,\
                                         g20_zz_yyyy, g20_zz_yyyz, g20_zz_yyzz,\
                                         g20_zz_yzzz, g20_zz_zzzz, g21_xx_xxxx,\
                                         g21_xx_xxxy, g21_xx_xxxz, g21_xx_xxyy,\
                                         g21_xx_xxyz, g21_xx_xxzz, g21_xx_xyyy,\
                                         g21_xx_xyyz, g21_xx_xyzz, g21_xx_xzzz,\
                                         g21_xx_yyyy, g21_xx_yyyz, g21_xx_yyzz,\
                                         g21_xx_yzzz, g21_xx_zzzz, g21_xy_xxxx,\
                                         g21_xy_xxxy, g21_xy_xxxz, g21_xy_xxyy,\
                                         g21_xy_xxyz, g21_xy_xxzz, g21_xy_xyyy,\
                                         g21_xy_xyyz, g21_xy_xyzz, g21_xy_xzzz,\
                                         g21_xy_yyyy, g21_xy_yyyz, g21_xy_yyzz,\
                                         g21_xy_yzzz, g21_xy_zzzz, g21_xz_xxxx,\
                                         g21_xz_xxxy, g21_xz_xxxz, g21_xz_xxyy,\
                                         g21_xz_xxyz, g21_xz_xxzz, g21_xz_xyyy,\
                                         g21_xz_xyyz, g21_xz_xyzz, g21_xz_xzzz,\
                                         g21_xz_yyyy, g21_xz_yyyz, g21_xz_yyzz,\
                                         g21_xz_yzzz, g21_xz_zzzz, g21_yy_xxxx,\
                                         g21_yy_xxxy, g21_yy_xxxz, g21_yy_xxyy,\
                                         g21_yy_xxyz, g21_yy_xxzz, g21_yy_xyyy,\
                                         g21_yy_xyyz, g21_yy_xyzz, g21_yy_xzzz,\
                                         g21_yy_yyyy, g21_yy_yyyz, g21_yy_yyzz,\
                                         g21_yy_yzzz, g21_yy_zzzz, g21_yz_xxxx,\
                                         g21_yz_xxxy, g21_yz_xxxz, g21_yz_xxyy,\
                                         g21_yz_xxyz, g21_yz_xxzz, g21_yz_xyyy,\
                                         g21_yz_xyyz, g21_yz_xyzz, g21_yz_xzzz,\
                                         g21_yz_yyyy, g21_yz_yyyz, g21_yz_yyzz,\
                                         g21_yz_yzzz, g21_yz_zzzz, g21_zz_xxxx,\
                                         g21_zz_xxxy, g21_zz_xxxz, g21_zz_xxyy,\
                                         g21_zz_xxyz, g21_zz_xxzz, g21_zz_xyyy,\
                                         g21_zz_xyyz, g21_zz_xyzz, g21_zz_xzzz,\
                                         g21_zz_yyyy, g21_zz_yyyz, g21_zz_yyzz,\
                                         g21_zz_yzzz, g21_zz_zzzz, g1_xxx_xxxx,\
                                         g1_xxx_xxxy, g1_xxx_xxxz, g1_xxx_xxyy,\
                                         g1_xxx_xxyz, g1_xxx_xxzz, g1_xxx_xyyy,\
                                         g1_xxx_xyyz, g1_xxx_xyzz, g1_xxx_xzzz,\
                                         g1_xxx_yyyy, g1_xxx_yyyz, g1_xxx_yyzz,\
                                         g1_xxx_yzzz, g1_xxx_zzzz, g1_xxy_xxxx,\
                                         g1_xxy_xxxy, g1_xxy_xxxz, g1_xxy_xxyy,\
                                         g1_xxy_xxyz, g1_xxy_xxzz, g1_xxy_xyyy,\
                                         g1_xxy_xyyz, g1_xxy_xyzz, g1_xxy_xzzz,\
                                         g1_xxy_yyyy, g1_xxy_yyyz, g1_xxy_yyzz,\
                                         g1_xxy_yzzz, g1_xxy_zzzz, g1_xxz_xxxx,\
                                         g1_xxz_xxxy, g1_xxz_xxxz, g1_xxz_xxyy,\
                                         g1_xxz_xxyz, g1_xxz_xxzz, g1_xxz_xyyy,\
                                         g1_xxz_xyyz, g1_xxz_xyzz, g1_xxz_xzzz,\
                                         g1_xxz_yyyy, g1_xxz_yyyz, g1_xxz_yyzz,\
                                         g1_xxz_yzzz, g1_xxz_zzzz, g1_xyy_xxxx,\
                                         g1_xyy_xxxy, g1_xyy_xxxz, g1_xyy_xxyy,\
                                         g1_xyy_xxyz, g1_xyy_xxzz, g1_xyy_xyyy,\
                                         g1_xyy_xyyz, g1_xyy_xyzz, g1_xyy_xzzz,\
                                         g1_xyy_yyyy, g1_xyy_yyyz, g1_xyy_yyzz,\
                                         g1_xyy_yzzz, g1_xyy_zzzz, g1_xyz_xxxx,\
                                         g1_xyz_xxxy, g1_xyz_xxxz, g1_xyz_xxyy,\
                                         g1_xyz_xxyz, g1_xyz_xxzz, g1_xyz_xyyy,\
                                         g1_xyz_xyyz, g1_xyz_xyzz, g1_xyz_xzzz,\
                                         g1_xyz_yyyy, g1_xyz_yyyz, g1_xyz_yyzz,\
                                         g1_xyz_yzzz, g1_xyz_zzzz, g1_xzz_xxxx,\
                                         g1_xzz_xxxy, g1_xzz_xxxz, g1_xzz_xxyy,\
                                         g1_xzz_xxyz, g1_xzz_xxzz, g1_xzz_xyyy,\
                                         g1_xzz_xyyz, g1_xzz_xyzz, g1_xzz_xzzz,\
                                         g1_xzz_yyyy, g1_xzz_yyyz, g1_xzz_yyzz,\
                                         g1_xzz_yzzz, g1_xzz_zzzz, g1_yyy_xxxx,\
                                         g1_yyy_xxxy, g1_yyy_xxxz, g1_yyy_xxyy,\
                                         g1_yyy_xxyz, g1_yyy_xxzz, g1_yyy_xyyy,\
                                         g1_yyy_xyyz, g1_yyy_xyzz, g1_yyy_xzzz,\
                                         g1_yyy_yyyy, g1_yyy_yyyz, g1_yyy_yyzz,\
                                         g1_yyy_yzzz, g1_yyy_zzzz, g1_yyz_xxxx,\
                                         g1_yyz_xxxy, g1_yyz_xxxz, g1_yyz_xxyy,\
                                         g1_yyz_xxyz, g1_yyz_xxzz, g1_yyz_xyyy,\
                                         g1_yyz_xyyz, g1_yyz_xyzz, g1_yyz_xzzz,\
                                         g1_yyz_yyyy, g1_yyz_yyyz, g1_yyz_yyzz,\
                                         g1_yyz_yzzz, g1_yyz_zzzz, g1_yzz_xxxx,\
                                         g1_yzz_xxxy, g1_yzz_xxxz, g1_yzz_xxyy,\
                                         g1_yzz_xxyz, g1_yzz_xxzz, g1_yzz_xyyy,\
                                         g1_yzz_xyyz, g1_yzz_xyzz, g1_yzz_xzzz,\
                                         g1_yzz_yyyy, g1_yzz_yyyz, g1_yzz_yyzz,\
                                         g1_yzz_yzzz, g1_yzz_zzzz, g1_zzz_xxxx,\
                                         g1_zzz_xxxy, g1_zzz_xxxz, g1_zzz_xxyy,\
                                         g1_zzz_xxyz, g1_zzz_xxzz, g1_zzz_xyyy,\
                                         g1_zzz_xyyz, g1_zzz_xyzz, g1_zzz_xzzz,\
                                         g1_zzz_yyyy, g1_zzz_yyyz, g1_zzz_yyzz,\
                                         g1_zzz_yzzz, g1_zzz_zzzz, g_xxxx_xxxx,\
                                         g_xxxx_xxxy, g_xxxx_xxxz, g_xxxx_xxyy,\
                                         g_xxxx_xxyz, g_xxxx_xxzz, g_xxxx_xyyy,\
                                         g_xxxx_xyyz, g_xxxx_xyzz, g_xxxx_xzzz,\
                                         g_xxxx_yyyy, g_xxxx_yyyz, g_xxxx_yyzz,\
                                         g_xxxx_yzzz, g_xxxx_zzzz, g_xxxy_xxxx,\
                                         g_xxxy_xxxy, g_xxxy_xxxz, g_xxxy_xxyy,\
                                         g_xxxy_xxyz, g_xxxy_xxzz, g_xxxy_xyyy,\
                                         g_xxxy_xyyz, g_xxxy_xyzz, g_xxxy_xzzz,\
                                         g_xxxy_yyyy, g_xxxy_yyyz, g_xxxy_yyzz,\
                                         g_xxxy_yzzz, g_xxxy_zzzz, g_xxxz_xxxx,\
                                         g_xxxz_xxxy, g_xxxz_xxxz, g_xxxz_xxyy,\
                                         g_xxxz_xxyz, g_xxxz_xxzz, g_xxxz_xyyy,\
                                         g_xxxz_xyyz, g_xxxz_xyzz, g_xxxz_xzzz,\
                                         g_xxxz_yyyy, g_xxxz_yyyz, g_xxxz_yyzz,\
                                         g_xxxz_yzzz, g_xxxz_zzzz, g_xxyy_xxxx,\
                                         g_xxyy_xxxy, g_xxyy_xxxz, g_xxyy_xxyy,\
                                         g_xxyy_xxyz, g_xxyy_xxzz, g_xxyy_xyyy,\
                                         g_xxyy_xyyz, g_xxyy_xyzz, g_xxyy_xzzz,\
                                         g_xxyy_yyyy, g_xxyy_yyyz, g_xxyy_yyzz,\
                                         g_xxyy_yzzz, g_xxyy_zzzz, g_xxyz_xxxx,\
                                         g_xxyz_xxxy, g_xxyz_xxxz, g_xxyz_xxyy,\
                                         g_xxyz_xxyz, g_xxyz_xxzz, g_xxyz_xyyy,\
                                         g_xxyz_xyyz, g_xxyz_xyzz, g_xxyz_xzzz,\
                                         g_xxyz_yyyy, g_xxyz_yyyz, g_xxyz_yyzz,\
                                         g_xxyz_yzzz, g_xxyz_zzzz, g_xxzz_xxxx,\
                                         g_xxzz_xxxy, g_xxzz_xxxz, g_xxzz_xxyy,\
                                         g_xxzz_xxyz, g_xxzz_xxzz, g_xxzz_xyyy,\
                                         g_xxzz_xyyz, g_xxzz_xyzz, g_xxzz_xzzz,\
                                         g_xxzz_yyyy, g_xxzz_yyyz, g_xxzz_yyzz,\
                                         g_xxzz_yzzz, g_xxzz_zzzz, g_xyyy_xxxx,\
                                         g_xyyy_xxxy, g_xyyy_xxxz, g_xyyy_xxyy,\
                                         g_xyyy_xxyz, g_xyyy_xxzz, g_xyyy_xyyy,\
                                         g_xyyy_xyyz, g_xyyy_xyzz, g_xyyy_xzzz,\
                                         g_xyyy_yyyy, g_xyyy_yyyz, g_xyyy_yyzz,\
                                         g_xyyy_yzzz, g_xyyy_zzzz, g_xyyz_xxxx,\
                                         g_xyyz_xxxy, g_xyyz_xxxz, g_xyyz_xxyy,\
                                         g_xyyz_xxyz, g_xyyz_xxzz, g_xyyz_xyyy,\
                                         g_xyyz_xyyz, g_xyyz_xyzz, g_xyyz_xzzz,\
                                         g_xyyz_yyyy, g_xyyz_yyyz, g_xyyz_yyzz,\
                                         g_xyyz_yzzz, g_xyyz_zzzz, g_xyzz_xxxx,\
                                         g_xyzz_xxxy, g_xyzz_xxxz, g_xyzz_xxyy,\
                                         g_xyzz_xxyz, g_xyzz_xxzz, g_xyzz_xyyy,\
                                         g_xyzz_xyyz, g_xyzz_xyzz, g_xyzz_xzzz,\
                                         g_xyzz_yyyy, g_xyzz_yyyz, g_xyzz_yyzz,\
                                         g_xyzz_yzzz, g_xyzz_zzzz, g_xzzz_xxxx,\
                                         g_xzzz_xxxy, g_xzzz_xxxz, g_xzzz_xxyy,\
                                         g_xzzz_xxyz, g_xzzz_xxzz, g_xzzz_xyyy,\
                                         g_xzzz_xyyz, g_xzzz_xyzz, g_xzzz_xzzz,\
                                         g_xzzz_yyyy, g_xzzz_yyyz, g_xzzz_yyzz,\
                                         g_xzzz_yzzz, g_xzzz_zzzz, g_yyyy_xxxx,\
                                         g_yyyy_xxxy, g_yyyy_xxxz, g_yyyy_xxyy,\
                                         g_yyyy_xxyz, g_yyyy_xxzz, g_yyyy_xyyy,\
                                         g_yyyy_xyyz, g_yyyy_xyzz, g_yyyy_xzzz,\
                                         g_yyyy_yyyy, g_yyyy_yyyz, g_yyyy_yyzz,\
                                         g_yyyy_yzzz, g_yyyy_zzzz, g_yyyz_xxxx,\
                                         g_yyyz_xxxy, g_yyyz_xxxz, g_yyyz_xxyy,\
                                         g_yyyz_xxyz, g_yyyz_xxzz, g_yyyz_xyyy,\
                                         g_yyyz_xyyz, g_yyyz_xyzz, g_yyyz_xzzz,\
                                         g_yyyz_yyyy, g_yyyz_yyyz, g_yyyz_yyzz,\
                                         g_yyyz_yzzz, g_yyyz_zzzz, g_yyzz_xxxx,\
                                         g_yyzz_xxxy, g_yyzz_xxxz, g_yyzz_xxyy,\
                                         g_yyzz_xxyz, g_yyzz_xxzz, g_yyzz_xyyy,\
                                         g_yyzz_xyyz, g_yyzz_xyzz, g_yyzz_xzzz,\
                                         g_yyzz_yyyy, g_yyzz_yyyz, g_yyzz_yyzz,\
                                         g_yyzz_yzzz, g_yyzz_zzzz, g_yzzz_xxxx,\
                                         g_yzzz_xxxy, g_yzzz_xxxz, g_yzzz_xxyy,\
                                         g_yzzz_xxyz, g_yzzz_xxzz, g_yzzz_xyyy,\
                                         g_yzzz_xyyz, g_yzzz_xyzz, g_yzzz_xzzz,\
                                         g_yzzz_yyyy, g_yzzz_yyyz, g_yzzz_yyzz,\
                                         g_yzzz_yzzz, g_yzzz_zzzz, g_zzzz_xxxx,\
                                         g_zzzz_xxxy, g_zzzz_xxxz, g_zzzz_xxyy,\
                                         g_zzzz_xxyz, g_zzzz_xxzz, g_zzzz_xyyy,\
                                         g_zzzz_xyyz, g_zzzz_xyzz, g_zzzz_xzzz,\
                                         g_zzzz_yyyy, g_zzzz_yyyz, g_zzzz_yyzz,\
                                         g_zzzz_yzzz, g_zzzz_zzzz: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                // scaled prefactor for ket

                double f2t = 0.50 * fx[k];

                // scaled prefactors for bra

                double f2g = 0.50 * fg[k];

                double fgz = fz[k];

                // leading x component

                double fr = pax[k];

                g_xxxx_xxxx[k] = fr * g1_xxx_xxxx[k] + f2g * (3.0 * g20_xx_xxxx[k] - 3.0 * fgz * g21_xx_xxxx[k]) + 4.0 * f2t * gk_xxx_xxx[k];

                g_xxxx_xxxy[k] = fr * g1_xxx_xxxy[k] + f2g * (3.0 * g20_xx_xxxy[k] - 3.0 * fgz * g21_xx_xxxy[k]) + 3.0 * f2t * gk_xxx_xxy[k];

                g_xxxx_xxxz[k] = fr * g1_xxx_xxxz[k] + f2g * (3.0 * g20_xx_xxxz[k] - 3.0 * fgz * g21_xx_xxxz[k]) + 3.0 * f2t * gk_xxx_xxz[k];

                g_xxxx_xxyy[k] = fr * g1_xxx_xxyy[k] + f2g * (3.0 * g20_xx_xxyy[k] - 3.0 * fgz * g21_xx_xxyy[k]) + 2.0 * f2t * gk_xxx_xyy[k];

                g_xxxx_xxyz[k] = fr * g1_xxx_xxyz[k] + f2g * (3.0 * g20_xx_xxyz[k] - 3.0 * fgz * g21_xx_xxyz[k]) + 2.0 * f2t * gk_xxx_xyz[k];

                g_xxxx_xxzz[k] = fr * g1_xxx_xxzz[k] + f2g * (3.0 * g20_xx_xxzz[k] - 3.0 * fgz * g21_xx_xxzz[k]) + 2.0 * f2t * gk_xxx_xzz[k];

                g_xxxx_xyyy[k] = fr * g1_xxx_xyyy[k] + f2g * (3.0 * g20_xx_xyyy[k] - 3.0 * fgz * g21_xx_xyyy[k]) + f2t * gk_xxx_yyy[k];

                g_xxxx_xyyz[k] = fr * g1_xxx_xyyz[k] + f2g * (3.0 * g20_xx_xyyz[k] - 3.0 * fgz * g21_xx_xyyz[k]) + f2t * gk_xxx_yyz[k];

                g_xxxx_xyzz[k] = fr * g1_xxx_xyzz[k] + f2g * (3.0 * g20_xx_xyzz[k] - 3.0 * fgz * g21_xx_xyzz[k]) + f2t * gk_xxx_yzz[k];

                g_xxxx_xzzz[k] = fr * g1_xxx_xzzz[k] + f2g * (3.0 * g20_xx_xzzz[k] - 3.0 * fgz * g21_xx_xzzz[k]) + f2t * gk_xxx_zzz[k];

                g_xxxx_yyyy[k] = fr * g1_xxx_yyyy[k] + f2g * (3.0 * g20_xx_yyyy[k] - 3.0 * fgz * g21_xx_yyyy[k]);

                g_xxxx_yyyz[k] = fr * g1_xxx_yyyz[k] + f2g * (3.0 * g20_xx_yyyz[k] - 3.0 * fgz * g21_xx_yyyz[k]);

                g_xxxx_yyzz[k] = fr * g1_xxx_yyzz[k] + f2g * (3.0 * g20_xx_yyzz[k] - 3.0 * fgz * g21_xx_yyzz[k]);

                g_xxxx_yzzz[k] = fr * g1_xxx_yzzz[k] + f2g * (3.0 * g20_xx_yzzz[k] - 3.0 * fgz * g21_xx_yzzz[k]);

                g_xxxx_zzzz[k] = fr * g1_xxx_zzzz[k] + f2g * (3.0 * g20_xx_zzzz[k] - 3.0 * fgz * g21_xx_zzzz[k]);

                g_xxxy_xxxx[k] = fr * g1_xxy_xxxx[k] + f2g * (2.0 * g20_xy_xxxx[k] - 2.0 * fgz * g21_xy_xxxx[k]) + 4.0 * f2t * gk_xxy_xxx[k];

                g_xxxy_xxxy[k] = fr * g1_xxy_xxxy[k] + f2g * (2.0 * g20_xy_xxxy[k] - 2.0 * fgz * g21_xy_xxxy[k]) + 3.0 * f2t * gk_xxy_xxy[k];

                g_xxxy_xxxz[k] = fr * g1_xxy_xxxz[k] + f2g * (2.0 * g20_xy_xxxz[k] - 2.0 * fgz * g21_xy_xxxz[k]) + 3.0 * f2t * gk_xxy_xxz[k];

                g_xxxy_xxyy[k] = fr * g1_xxy_xxyy[k] + f2g * (2.0 * g20_xy_xxyy[k] - 2.0 * fgz * g21_xy_xxyy[k]) + 2.0 * f2t * gk_xxy_xyy[k];

                g_xxxy_xxyz[k] = fr * g1_xxy_xxyz[k] + f2g * (2.0 * g20_xy_xxyz[k] - 2.0 * fgz * g21_xy_xxyz[k]) + 2.0 * f2t * gk_xxy_xyz[k];

                g_xxxy_xxzz[k] = fr * g1_xxy_xxzz[k] + f2g * (2.0 * g20_xy_xxzz[k] - 2.0 * fgz * g21_xy_xxzz[k]) + 2.0 * f2t * gk_xxy_xzz[k];

                g_xxxy_xyyy[k] = fr * g1_xxy_xyyy[k] + f2g * (2.0 * g20_xy_xyyy[k] - 2.0 * fgz * g21_xy_xyyy[k]) + f2t * gk_xxy_yyy[k];

                g_xxxy_xyyz[k] = fr * g1_xxy_xyyz[k] + f2g * (2.0 * g20_xy_xyyz[k] - 2.0 * fgz * g21_xy_xyyz[k]) + f2t * gk_xxy_yyz[k];

                g_xxxy_xyzz[k] = fr * g1_xxy_xyzz[k] + f2g * (2.0 * g20_xy_xyzz[k] - 2.0 * fgz * g21_xy_xyzz[k]) + f2t * gk_xxy_yzz[k];

                g_xxxy_xzzz[k] = fr * g1_xxy_xzzz[k] + f2g * (2.0 * g20_xy_xzzz[k] - 2.0 * fgz * g21_xy_xzzz[k]) + f2t * gk_xxy_zzz[k];

                g_xxxy_yyyy[k] = fr * g1_xxy_yyyy[k] + f2g * (2.0 * g20_xy_yyyy[k] - 2.0 * fgz * g21_xy_yyyy[k]);

                g_xxxy_yyyz[k] = fr * g1_xxy_yyyz[k] + f2g * (2.0 * g20_xy_yyyz[k] - 2.0 * fgz * g21_xy_yyyz[k]);

                g_xxxy_yyzz[k] = fr * g1_xxy_yyzz[k] + f2g * (2.0 * g20_xy_yyzz[k] - 2.0 * fgz * g21_xy_yyzz[k]);

                g_xxxy_yzzz[k] = fr * g1_xxy_yzzz[k] + f2g * (2.0 * g20_xy_yzzz[k] - 2.0 * fgz * g21_xy_yzzz[k]);

                g_xxxy_zzzz[k] = fr * g1_xxy_zzzz[k] + f2g * (2.0 * g20_xy_zzzz[k] - 2.0 * fgz * g21_xy_zzzz[k]);

                g_xxxz_xxxx[k] = fr * g1_xxz_xxxx[k] + f2g * (2.0 * g20_xz_xxxx[k] - 2.0 * fgz * g21_xz_xxxx[k]) + 4.0 * f2t * gk_xxz_xxx[k];

                g_xxxz_xxxy[k] = fr * g1_xxz_xxxy[k] + f2g * (2.0 * g20_xz_xxxy[k] - 2.0 * fgz * g21_xz_xxxy[k]) + 3.0 * f2t * gk_xxz_xxy[k];

                g_xxxz_xxxz[k] = fr * g1_xxz_xxxz[k] + f2g * (2.0 * g20_xz_xxxz[k] - 2.0 * fgz * g21_xz_xxxz[k]) + 3.0 * f2t * gk_xxz_xxz[k];

                g_xxxz_xxyy[k] = fr * g1_xxz_xxyy[k] + f2g * (2.0 * g20_xz_xxyy[k] - 2.0 * fgz * g21_xz_xxyy[k]) + 2.0 * f2t * gk_xxz_xyy[k];

                g_xxxz_xxyz[k] = fr * g1_xxz_xxyz[k] + f2g * (2.0 * g20_xz_xxyz[k] - 2.0 * fgz * g21_xz_xxyz[k]) + 2.0 * f2t * gk_xxz_xyz[k];

                g_xxxz_xxzz[k] = fr * g1_xxz_xxzz[k] + f2g * (2.0 * g20_xz_xxzz[k] - 2.0 * fgz * g21_xz_xxzz[k]) + 2.0 * f2t * gk_xxz_xzz[k];

                g_xxxz_xyyy[k] = fr * g1_xxz_xyyy[k] + f2g * (2.0 * g20_xz_xyyy[k] - 2.0 * fgz * g21_xz_xyyy[k]) + f2t * gk_xxz_yyy[k];

                g_xxxz_xyyz[k] = fr * g1_xxz_xyyz[k] + f2g * (2.0 * g20_xz_xyyz[k] - 2.0 * fgz * g21_xz_xyyz[k]) + f2t * gk_xxz_yyz[k];

                g_xxxz_xyzz[k] = fr * g1_xxz_xyzz[k] + f2g * (2.0 * g20_xz_xyzz[k] - 2.0 * fgz * g21_xz_xyzz[k]) + f2t * gk_xxz_yzz[k];

                g_xxxz_xzzz[k] = fr * g1_xxz_xzzz[k] + f2g * (2.0 * g20_xz_xzzz[k] - 2.0 * fgz * g21_xz_xzzz[k]) + f2t * gk_xxz_zzz[k];

                g_xxxz_yyyy[k] = fr * g1_xxz_yyyy[k] + f2g * (2.0 * g20_xz_yyyy[k] - 2.0 * fgz * g21_xz_yyyy[k]);

                g_xxxz_yyyz[k] = fr * g1_xxz_yyyz[k] + f2g * (2.0 * g20_xz_yyyz[k] - 2.0 * fgz * g21_xz_yyyz[k]);

                g_xxxz_yyzz[k] = fr * g1_xxz_yyzz[k] + f2g * (2.0 * g20_xz_yyzz[k] - 2.0 * fgz * g21_xz_yyzz[k]);

                g_xxxz_yzzz[k] = fr * g1_xxz_yzzz[k] + f2g * (2.0 * g20_xz_yzzz[k] - 2.0 * fgz * g21_xz_yzzz[k]);

                g_xxxz_zzzz[k] = fr * g1_xxz_zzzz[k] + f2g * (2.0 * g20_xz_zzzz[k] - 2.0 * fgz * g21_xz_zzzz[k]);

                g_xxyy_xxxx[k] = fr * g1_xyy_xxxx[k] + f2g * (g20_yy_xxxx[k] - fgz * g21_yy_xxxx[k]) + 4.0 * f2t * gk_xyy_xxx[k];

                g_xxyy_xxxy[k] = fr * g1_xyy_xxxy[k] + f2g * (g20_yy_xxxy[k] - fgz * g21_yy_xxxy[k]) + 3.0 * f2t * gk_xyy_xxy[k];

                g_xxyy_xxxz[k] = fr * g1_xyy_xxxz[k] + f2g * (g20_yy_xxxz[k] - fgz * g21_yy_xxxz[k]) + 3.0 * f2t * gk_xyy_xxz[k];

                g_xxyy_xxyy[k] = fr * g1_xyy_xxyy[k] + f2g * (g20_yy_xxyy[k] - fgz * g21_yy_xxyy[k]) + 2.0 * f2t * gk_xyy_xyy[k];

                g_xxyy_xxyz[k] = fr * g1_xyy_xxyz[k] + f2g * (g20_yy_xxyz[k] - fgz * g21_yy_xxyz[k]) + 2.0 * f2t * gk_xyy_xyz[k];

                g_xxyy_xxzz[k] = fr * g1_xyy_xxzz[k] + f2g * (g20_yy_xxzz[k] - fgz * g21_yy_xxzz[k]) + 2.0 * f2t * gk_xyy_xzz[k];

                g_xxyy_xyyy[k] = fr * g1_xyy_xyyy[k] + f2g * (g20_yy_xyyy[k] - fgz * g21_yy_xyyy[k]) + f2t * gk_xyy_yyy[k];

                g_xxyy_xyyz[k] = fr * g1_xyy_xyyz[k] + f2g * (g20_yy_xyyz[k] - fgz * g21_yy_xyyz[k]) + f2t * gk_xyy_yyz[k];

                g_xxyy_xyzz[k] = fr * g1_xyy_xyzz[k] + f2g * (g20_yy_xyzz[k] - fgz * g21_yy_xyzz[k]) + f2t * gk_xyy_yzz[k];

                g_xxyy_xzzz[k] = fr * g1_xyy_xzzz[k] + f2g * (g20_yy_xzzz[k] - fgz * g21_yy_xzzz[k]) + f2t * gk_xyy_zzz[k];

                g_xxyy_yyyy[k] = fr * g1_xyy_yyyy[k] + f2g * (g20_yy_yyyy[k] - fgz * g21_yy_yyyy[k]);

                g_xxyy_yyyz[k] = fr * g1_xyy_yyyz[k] + f2g * (g20_yy_yyyz[k] - fgz * g21_yy_yyyz[k]);

                g_xxyy_yyzz[k] = fr * g1_xyy_yyzz[k] + f2g * (g20_yy_yyzz[k] - fgz * g21_yy_yyzz[k]);

                g_xxyy_yzzz[k] = fr * g1_xyy_yzzz[k] + f2g * (g20_yy_yzzz[k] - fgz * g21_yy_yzzz[k]);

                g_xxyy_zzzz[k] = fr * g1_xyy_zzzz[k] + f2g * (g20_yy_zzzz[k] - fgz * g21_yy_zzzz[k]);

                g_xxyz_xxxx[k] = fr * g1_xyz_xxxx[k] + f2g * (g20_yz_xxxx[k] - fgz * g21_yz_xxxx[k]) + 4.0 * f2t * gk_xyz_xxx[k];

                g_xxyz_xxxy[k] = fr * g1_xyz_xxxy[k] + f2g * (g20_yz_xxxy[k] - fgz * g21_yz_xxxy[k]) + 3.0 * f2t * gk_xyz_xxy[k];

                g_xxyz_xxxz[k] = fr * g1_xyz_xxxz[k] + f2g * (g20_yz_xxxz[k] - fgz * g21_yz_xxxz[k]) + 3.0 * f2t * gk_xyz_xxz[k];

                g_xxyz_xxyy[k] = fr * g1_xyz_xxyy[k] + f2g * (g20_yz_xxyy[k] - fgz * g21_yz_xxyy[k]) + 2.0 * f2t * gk_xyz_xyy[k];

                g_xxyz_xxyz[k] = fr * g1_xyz_xxyz[k] + f2g * (g20_yz_xxyz[k] - fgz * g21_yz_xxyz[k]) + 2.0 * f2t * gk_xyz_xyz[k];

                g_xxyz_xxzz[k] = fr * g1_xyz_xxzz[k] + f2g * (g20_yz_xxzz[k] - fgz * g21_yz_xxzz[k]) + 2.0 * f2t * gk_xyz_xzz[k];

                g_xxyz_xyyy[k] = fr * g1_xyz_xyyy[k] + f2g * (g20_yz_xyyy[k] - fgz * g21_yz_xyyy[k]) + f2t * gk_xyz_yyy[k];

                g_xxyz_xyyz[k] = fr * g1_xyz_xyyz[k] + f2g * (g20_yz_xyyz[k] - fgz * g21_yz_xyyz[k]) + f2t * gk_xyz_yyz[k];

                g_xxyz_xyzz[k] = fr * g1_xyz_xyzz[k] + f2g * (g20_yz_xyzz[k] - fgz * g21_yz_xyzz[k]) + f2t * gk_xyz_yzz[k];

                g_xxyz_xzzz[k] = fr * g1_xyz_xzzz[k] + f2g * (g20_yz_xzzz[k] - fgz * g21_yz_xzzz[k]) + f2t * gk_xyz_zzz[k];

                g_xxyz_yyyy[k] = fr * g1_xyz_yyyy[k] + f2g * (g20_yz_yyyy[k] - fgz * g21_yz_yyyy[k]);

                g_xxyz_yyyz[k] = fr * g1_xyz_yyyz[k] + f2g * (g20_yz_yyyz[k] - fgz * g21_yz_yyyz[k]);

                g_xxyz_yyzz[k] = fr * g1_xyz_yyzz[k] + f2g * (g20_yz_yyzz[k] - fgz * g21_yz_yyzz[k]);

                g_xxyz_yzzz[k] = fr * g1_xyz_yzzz[k] + f2g * (g20_yz_yzzz[k] - fgz * g21_yz_yzzz[k]);

                g_xxyz_zzzz[k] = fr * g1_xyz_zzzz[k] + f2g * (g20_yz_zzzz[k] - fgz * g21_yz_zzzz[k]);

                g_xxzz_xxxx[k] = fr * g1_xzz_xxxx[k] + f2g * (g20_zz_xxxx[k] - fgz * g21_zz_xxxx[k]) + 4.0 * f2t * gk_xzz_xxx[k];

                g_xxzz_xxxy[k] = fr * g1_xzz_xxxy[k] + f2g * (g20_zz_xxxy[k] - fgz * g21_zz_xxxy[k]) + 3.0 * f2t * gk_xzz_xxy[k];

                g_xxzz_xxxz[k] = fr * g1_xzz_xxxz[k] + f2g * (g20_zz_xxxz[k] - fgz * g21_zz_xxxz[k]) + 3.0 * f2t * gk_xzz_xxz[k];

                g_xxzz_xxyy[k] = fr * g1_xzz_xxyy[k] + f2g * (g20_zz_xxyy[k] - fgz * g21_zz_xxyy[k]) + 2.0 * f2t * gk_xzz_xyy[k];

                g_xxzz_xxyz[k] = fr * g1_xzz_xxyz[k] + f2g * (g20_zz_xxyz[k] - fgz * g21_zz_xxyz[k]) + 2.0 * f2t * gk_xzz_xyz[k];

                g_xxzz_xxzz[k] = fr * g1_xzz_xxzz[k] + f2g * (g20_zz_xxzz[k] - fgz * g21_zz_xxzz[k]) + 2.0 * f2t * gk_xzz_xzz[k];

                g_xxzz_xyyy[k] = fr * g1_xzz_xyyy[k] + f2g * (g20_zz_xyyy[k] - fgz * g21_zz_xyyy[k]) + f2t * gk_xzz_yyy[k];

                g_xxzz_xyyz[k] = fr * g1_xzz_xyyz[k] + f2g * (g20_zz_xyyz[k] - fgz * g21_zz_xyyz[k]) + f2t * gk_xzz_yyz[k];

                g_xxzz_xyzz[k] = fr * g1_xzz_xyzz[k] + f2g * (g20_zz_xyzz[k] - fgz * g21_zz_xyzz[k]) + f2t * gk_xzz_yzz[k];

                g_xxzz_xzzz[k] = fr * g1_xzz_xzzz[k] + f2g * (g20_zz_xzzz[k] - fgz * g21_zz_xzzz[k]) + f2t * gk_xzz_zzz[k];

                g_xxzz_yyyy[k] = fr * g1_xzz_yyyy[k] + f2g * (g20_zz_yyyy[k] - fgz * g21_zz_yyyy[k]);

                g_xxzz_yyyz[k] = fr * g1_xzz_yyyz[k] + f2g * (g20_zz_yyyz[k] - fgz * g21_zz_yyyz[k]);

                g_xxzz_yyzz[k] = fr * g1_xzz_yyzz[k] + f2g * (g20_zz_yyzz[k] - fgz * g21_zz_yyzz[k]);

                g_xxzz_yzzz[k] = fr * g1_xzz_yzzz[k] + f2g * (g20_zz_yzzz[k] - fgz * g21_zz_yzzz[k]);

                g_xxzz_zzzz[k] = fr * g1_xzz_zzzz[k] + f2g * (g20_zz_zzzz[k] - fgz * g21_zz_zzzz[k]);

                g_xyyy_xxxx[k] = fr * g1_yyy_xxxx[k] + 4.0 * f2t * gk_yyy_xxx[k];

                g_xyyy_xxxy[k] = fr * g1_yyy_xxxy[k] + 3.0 * f2t * gk_yyy_xxy[k];

                g_xyyy_xxxz[k] = fr * g1_yyy_xxxz[k] + 3.0 * f2t * gk_yyy_xxz[k];

                g_xyyy_xxyy[k] = fr * g1_yyy_xxyy[k] + 2.0 * f2t * gk_yyy_xyy[k];

                g_xyyy_xxyz[k] = fr * g1_yyy_xxyz[k] + 2.0 * f2t * gk_yyy_xyz[k];

                g_xyyy_xxzz[k] = fr * g1_yyy_xxzz[k] + 2.0 * f2t * gk_yyy_xzz[k];

                g_xyyy_xyyy[k] = fr * g1_yyy_xyyy[k] + f2t * gk_yyy_yyy[k];

                g_xyyy_xyyz[k] = fr * g1_yyy_xyyz[k] + f2t * gk_yyy_yyz[k];

                g_xyyy_xyzz[k] = fr * g1_yyy_xyzz[k] + f2t * gk_yyy_yzz[k];

                g_xyyy_xzzz[k] = fr * g1_yyy_xzzz[k] + f2t * gk_yyy_zzz[k];

                g_xyyy_yyyy[k] = fr * g1_yyy_yyyy[k];

                g_xyyy_yyyz[k] = fr * g1_yyy_yyyz[k];

                g_xyyy_yyzz[k] = fr * g1_yyy_yyzz[k];

                g_xyyy_yzzz[k] = fr * g1_yyy_yzzz[k];

                g_xyyy_zzzz[k] = fr * g1_yyy_zzzz[k];

                g_xyyz_xxxx[k] = fr * g1_yyz_xxxx[k] + 4.0 * f2t * gk_yyz_xxx[k];

                g_xyyz_xxxy[k] = fr * g1_yyz_xxxy[k] + 3.0 * f2t * gk_yyz_xxy[k];

                g_xyyz_xxxz[k] = fr * g1_yyz_xxxz[k] + 3.0 * f2t * gk_yyz_xxz[k];

                g_xyyz_xxyy[k] = fr * g1_yyz_xxyy[k] + 2.0 * f2t * gk_yyz_xyy[k];

                g_xyyz_xxyz[k] = fr * g1_yyz_xxyz[k] + 2.0 * f2t * gk_yyz_xyz[k];

                g_xyyz_xxzz[k] = fr * g1_yyz_xxzz[k] + 2.0 * f2t * gk_yyz_xzz[k];

                g_xyyz_xyyy[k] = fr * g1_yyz_xyyy[k] + f2t * gk_yyz_yyy[k];

                g_xyyz_xyyz[k] = fr * g1_yyz_xyyz[k] + f2t * gk_yyz_yyz[k];

                g_xyyz_xyzz[k] = fr * g1_yyz_xyzz[k] + f2t * gk_yyz_yzz[k];

                g_xyyz_xzzz[k] = fr * g1_yyz_xzzz[k] + f2t * gk_yyz_zzz[k];

                g_xyyz_yyyy[k] = fr * g1_yyz_yyyy[k];

                g_xyyz_yyyz[k] = fr * g1_yyz_yyyz[k];

                g_xyyz_yyzz[k] = fr * g1_yyz_yyzz[k];

                g_xyyz_yzzz[k] = fr * g1_yyz_yzzz[k];

                g_xyyz_zzzz[k] = fr * g1_yyz_zzzz[k];

                g_xyzz_xxxx[k] = fr * g1_yzz_xxxx[k] + 4.0 * f2t * gk_yzz_xxx[k];

                g_xyzz_xxxy[k] = fr * g1_yzz_xxxy[k] + 3.0 * f2t * gk_yzz_xxy[k];

                g_xyzz_xxxz[k] = fr * g1_yzz_xxxz[k] + 3.0 * f2t * gk_yzz_xxz[k];

                g_xyzz_xxyy[k] = fr * g1_yzz_xxyy[k] + 2.0 * f2t * gk_yzz_xyy[k];

                g_xyzz_xxyz[k] = fr * g1_yzz_xxyz[k] + 2.0 * f2t * gk_yzz_xyz[k];

                g_xyzz_xxzz[k] = fr * g1_yzz_xxzz[k] + 2.0 * f2t * gk_yzz_xzz[k];

                g_xyzz_xyyy[k] = fr * g1_yzz_xyyy[k] + f2t * gk_yzz_yyy[k];

                g_xyzz_xyyz[k] = fr * g1_yzz_xyyz[k] + f2t * gk_yzz_yyz[k];

                g_xyzz_xyzz[k] = fr * g1_yzz_xyzz[k] + f2t * gk_yzz_yzz[k];

                g_xyzz_xzzz[k] = fr * g1_yzz_xzzz[k] + f2t * gk_yzz_zzz[k];

                g_xyzz_yyyy[k] = fr * g1_yzz_yyyy[k];

                g_xyzz_yyyz[k] = fr * g1_yzz_yyyz[k];

                g_xyzz_yyzz[k] = fr * g1_yzz_yyzz[k];

                g_xyzz_yzzz[k] = fr * g1_yzz_yzzz[k];

                g_xyzz_zzzz[k] = fr * g1_yzz_zzzz[k];

                g_xzzz_xxxx[k] = fr * g1_zzz_xxxx[k] + 4.0 * f2t * gk_zzz_xxx[k];

                g_xzzz_xxxy[k] = fr * g1_zzz_xxxy[k] + 3.0 * f2t * gk_zzz_xxy[k];

                g_xzzz_xxxz[k] = fr * g1_zzz_xxxz[k] + 3.0 * f2t * gk_zzz_xxz[k];

                g_xzzz_xxyy[k] = fr * g1_zzz_xxyy[k] + 2.0 * f2t * gk_zzz_xyy[k];

                g_xzzz_xxyz[k] = fr * g1_zzz_xxyz[k] + 2.0 * f2t * gk_zzz_xyz[k];

                g_xzzz_xxzz[k] = fr * g1_zzz_xxzz[k] + 2.0 * f2t * gk_zzz_xzz[k];

                g_xzzz_xyyy[k] = fr * g1_zzz_xyyy[k] + f2t * gk_zzz_yyy[k];

                g_xzzz_xyyz[k] = fr * g1_zzz_xyyz[k] + f2t * gk_zzz_yyz[k];

                g_xzzz_xyzz[k] = fr * g1_zzz_xyzz[k] + f2t * gk_zzz_yzz[k];

                g_xzzz_xzzz[k] = fr * g1_zzz_xzzz[k] + f2t * gk_zzz_zzz[k];

                g_xzzz_yyyy[k] = fr * g1_zzz_yyyy[k];

                g_xzzz_yyyz[k] = fr * g1_zzz_yyyz[k];

                g_xzzz_yyzz[k] = fr * g1_zzz_yyzz[k];

                g_xzzz_yzzz[k] = fr * g1_zzz_yzzz[k];

                g_xzzz_zzzz[k] = fr * g1_zzz_zzzz[k];

                // leading y component

                fr = pay[k];

                g_yyyy_xxxx[k] = fr * g1_yyy_xxxx[k] + f2g * (3.0 * g20_yy_xxxx[k] - 3.0 * fgz * g21_yy_xxxx[k]);

                g_yyyy_xxxy[k] = fr * g1_yyy_xxxy[k] + f2g * (3.0 * g20_yy_xxxy[k] - 3.0 * fgz * g21_yy_xxxy[k]) + f2t * gk_yyy_xxx[k];

                g_yyyy_xxxz[k] = fr * g1_yyy_xxxz[k] + f2g * (3.0 * g20_yy_xxxz[k] - 3.0 * fgz * g21_yy_xxxz[k]);

                g_yyyy_xxyy[k] = fr * g1_yyy_xxyy[k] + f2g * (3.0 * g20_yy_xxyy[k] - 3.0 * fgz * g21_yy_xxyy[k]) + 2.0 * f2t * gk_yyy_xxy[k];

                g_yyyy_xxyz[k] = fr * g1_yyy_xxyz[k] + f2g * (3.0 * g20_yy_xxyz[k] - 3.0 * fgz * g21_yy_xxyz[k]) + f2t * gk_yyy_xxz[k];

                g_yyyy_xxzz[k] = fr * g1_yyy_xxzz[k] + f2g * (3.0 * g20_yy_xxzz[k] - 3.0 * fgz * g21_yy_xxzz[k]);

                g_yyyy_xyyy[k] = fr * g1_yyy_xyyy[k] + f2g * (3.0 * g20_yy_xyyy[k] - 3.0 * fgz * g21_yy_xyyy[k]) + 3.0 * f2t * gk_yyy_xyy[k];

                g_yyyy_xyyz[k] = fr * g1_yyy_xyyz[k] + f2g * (3.0 * g20_yy_xyyz[k] - 3.0 * fgz * g21_yy_xyyz[k]) + 2.0 * f2t * gk_yyy_xyz[k];

                g_yyyy_xyzz[k] = fr * g1_yyy_xyzz[k] + f2g * (3.0 * g20_yy_xyzz[k] - 3.0 * fgz * g21_yy_xyzz[k]) + f2t * gk_yyy_xzz[k];

                g_yyyy_xzzz[k] = fr * g1_yyy_xzzz[k] + f2g * (3.0 * g20_yy_xzzz[k] - 3.0 * fgz * g21_yy_xzzz[k]);

                g_yyyy_yyyy[k] = fr * g1_yyy_yyyy[k] + f2g * (3.0 * g20_yy_yyyy[k] - 3.0 * fgz * g21_yy_yyyy[k]) + 4.0 * f2t * gk_yyy_yyy[k];

                g_yyyy_yyyz[k] = fr * g1_yyy_yyyz[k] + f2g * (3.0 * g20_yy_yyyz[k] - 3.0 * fgz * g21_yy_yyyz[k]) + 3.0 * f2t * gk_yyy_yyz[k];

                g_yyyy_yyzz[k] = fr * g1_yyy_yyzz[k] + f2g * (3.0 * g20_yy_yyzz[k] - 3.0 * fgz * g21_yy_yyzz[k]) + 2.0 * f2t * gk_yyy_yzz[k];

                g_yyyy_yzzz[k] = fr * g1_yyy_yzzz[k] + f2g * (3.0 * g20_yy_yzzz[k] - 3.0 * fgz * g21_yy_yzzz[k]) + f2t * gk_yyy_zzz[k];

                g_yyyy_zzzz[k] = fr * g1_yyy_zzzz[k] + f2g * (3.0 * g20_yy_zzzz[k] - 3.0 * fgz * g21_yy_zzzz[k]);

                g_yyyz_xxxx[k] = fr * g1_yyz_xxxx[k] + f2g * (2.0 * g20_yz_xxxx[k] - 2.0 * fgz * g21_yz_xxxx[k]);

                g_yyyz_xxxy[k] = fr * g1_yyz_xxxy[k] + f2g * (2.0 * g20_yz_xxxy[k] - 2.0 * fgz * g21_yz_xxxy[k]) + f2t * gk_yyz_xxx[k];

                g_yyyz_xxxz[k] = fr * g1_yyz_xxxz[k] + f2g * (2.0 * g20_yz_xxxz[k] - 2.0 * fgz * g21_yz_xxxz[k]);

                g_yyyz_xxyy[k] = fr * g1_yyz_xxyy[k] + f2g * (2.0 * g20_yz_xxyy[k] - 2.0 * fgz * g21_yz_xxyy[k]) + 2.0 * f2t * gk_yyz_xxy[k];

                g_yyyz_xxyz[k] = fr * g1_yyz_xxyz[k] + f2g * (2.0 * g20_yz_xxyz[k] - 2.0 * fgz * g21_yz_xxyz[k]) + f2t * gk_yyz_xxz[k];

                g_yyyz_xxzz[k] = fr * g1_yyz_xxzz[k] + f2g * (2.0 * g20_yz_xxzz[k] - 2.0 * fgz * g21_yz_xxzz[k]);

                g_yyyz_xyyy[k] = fr * g1_yyz_xyyy[k] + f2g * (2.0 * g20_yz_xyyy[k] - 2.0 * fgz * g21_yz_xyyy[k]) + 3.0 * f2t * gk_yyz_xyy[k];

                g_yyyz_xyyz[k] = fr * g1_yyz_xyyz[k] + f2g * (2.0 * g20_yz_xyyz[k] - 2.0 * fgz * g21_yz_xyyz[k]) + 2.0 * f2t * gk_yyz_xyz[k];

                g_yyyz_xyzz[k] = fr * g1_yyz_xyzz[k] + f2g * (2.0 * g20_yz_xyzz[k] - 2.0 * fgz * g21_yz_xyzz[k]) + f2t * gk_yyz_xzz[k];

                g_yyyz_xzzz[k] = fr * g1_yyz_xzzz[k] + f2g * (2.0 * g20_yz_xzzz[k] - 2.0 * fgz * g21_yz_xzzz[k]);

                g_yyyz_yyyy[k] = fr * g1_yyz_yyyy[k] + f2g * (2.0 * g20_yz_yyyy[k] - 2.0 * fgz * g21_yz_yyyy[k]) + 4.0 * f2t * gk_yyz_yyy[k];

                g_yyyz_yyyz[k] = fr * g1_yyz_yyyz[k] + f2g * (2.0 * g20_yz_yyyz[k] - 2.0 * fgz * g21_yz_yyyz[k]) + 3.0 * f2t * gk_yyz_yyz[k];

                g_yyyz_yyzz[k] = fr * g1_yyz_yyzz[k] + f2g * (2.0 * g20_yz_yyzz[k] - 2.0 * fgz * g21_yz_yyzz[k]) + 2.0 * f2t * gk_yyz_yzz[k];

                g_yyyz_yzzz[k] = fr * g1_yyz_yzzz[k] + f2g * (2.0 * g20_yz_yzzz[k] - 2.0 * fgz * g21_yz_yzzz[k]) + f2t * gk_yyz_zzz[k];

                g_yyyz_zzzz[k] = fr * g1_yyz_zzzz[k] + f2g * (2.0 * g20_yz_zzzz[k] - 2.0 * fgz * g21_yz_zzzz[k]);

                g_yyzz_xxxx[k] = fr * g1_yzz_xxxx[k] + f2g * (g20_zz_xxxx[k] - fgz * g21_zz_xxxx[k]);

                g_yyzz_xxxy[k] = fr * g1_yzz_xxxy[k] + f2g * (g20_zz_xxxy[k] - fgz * g21_zz_xxxy[k]) + f2t * gk_yzz_xxx[k];

                g_yyzz_xxxz[k] = fr * g1_yzz_xxxz[k] + f2g * (g20_zz_xxxz[k] - fgz * g21_zz_xxxz[k]);

                g_yyzz_xxyy[k] = fr * g1_yzz_xxyy[k] + f2g * (g20_zz_xxyy[k] - fgz * g21_zz_xxyy[k]) + 2.0 * f2t * gk_yzz_xxy[k];

                g_yyzz_xxyz[k] = fr * g1_yzz_xxyz[k] + f2g * (g20_zz_xxyz[k] - fgz * g21_zz_xxyz[k]) + f2t * gk_yzz_xxz[k];

                g_yyzz_xxzz[k] = fr * g1_yzz_xxzz[k] + f2g * (g20_zz_xxzz[k] - fgz * g21_zz_xxzz[k]);

                g_yyzz_xyyy[k] = fr * g1_yzz_xyyy[k] + f2g * (g20_zz_xyyy[k] - fgz * g21_zz_xyyy[k]) + 3.0 * f2t * gk_yzz_xyy[k];

                g_yyzz_xyyz[k] = fr * g1_yzz_xyyz[k] + f2g * (g20_zz_xyyz[k] - fgz * g21_zz_xyyz[k]) + 2.0 * f2t * gk_yzz_xyz[k];

                g_yyzz_xyzz[k] = fr * g1_yzz_xyzz[k] + f2g * (g20_zz_xyzz[k] - fgz * g21_zz_xyzz[k]) + f2t * gk_yzz_xzz[k];

                g_yyzz_xzzz[k] = fr * g1_yzz_xzzz[k] + f2g * (g20_zz_xzzz[k] - fgz * g21_zz_xzzz[k]);

                g_yyzz_yyyy[k] = fr * g1_yzz_yyyy[k] + f2g * (g20_zz_yyyy[k] - fgz * g21_zz_yyyy[k]) + 4.0 * f2t * gk_yzz_yyy[k];

                g_yyzz_yyyz[k] = fr * g1_yzz_yyyz[k] + f2g * (g20_zz_yyyz[k] - fgz * g21_zz_yyyz[k]) + 3.0 * f2t * gk_yzz_yyz[k];

                g_yyzz_yyzz[k] = fr * g1_yzz_yyzz[k] + f2g * (g20_zz_yyzz[k] - fgz * g21_zz_yyzz[k]) + 2.0 * f2t * gk_yzz_yzz[k];

                g_yyzz_yzzz[k] = fr * g1_yzz_yzzz[k] + f2g * (g20_zz_yzzz[k] - fgz * g21_zz_yzzz[k]) + f2t * gk_yzz_zzz[k];

                g_yyzz_zzzz[k] = fr * g1_yzz_zzzz[k] + f2g * (g20_zz_zzzz[k] - fgz * g21_zz_zzzz[k]);

                g_yzzz_xxxx[k] = fr * g1_zzz_xxxx[k];

                g_yzzz_xxxy[k] = fr * g1_zzz_xxxy[k] + f2t * gk_zzz_xxx[k];

                g_yzzz_xxxz[k] = fr * g1_zzz_xxxz[k];

                g_yzzz_xxyy[k] = fr * g1_zzz_xxyy[k] + 2.0 * f2t * gk_zzz_xxy[k];

                g_yzzz_xxyz[k] = fr * g1_zzz_xxyz[k] + f2t * gk_zzz_xxz[k];

                g_yzzz_xxzz[k] = fr * g1_zzz_xxzz[k];

                g_yzzz_xyyy[k] = fr * g1_zzz_xyyy[k] + 3.0 * f2t * gk_zzz_xyy[k];

                g_yzzz_xyyz[k] = fr * g1_zzz_xyyz[k] + 2.0 * f2t * gk_zzz_xyz[k];

                g_yzzz_xyzz[k] = fr * g1_zzz_xyzz[k] + f2t * gk_zzz_xzz[k];

                g_yzzz_xzzz[k] = fr * g1_zzz_xzzz[k];

                g_yzzz_yyyy[k] = fr * g1_zzz_yyyy[k] + 4.0 * f2t * gk_zzz_yyy[k];

                g_yzzz_yyyz[k] = fr * g1_zzz_yyyz[k] + 3.0 * f2t * gk_zzz_yyz[k];

                g_yzzz_yyzz[k] = fr * g1_zzz_yyzz[k] + 2.0 * f2t * gk_zzz_yzz[k];

                g_yzzz_yzzz[k] = fr * g1_zzz_yzzz[k] + f2t * gk_zzz_zzz[k];

                g_yzzz_zzzz[k] = fr * g1_zzz_zzzz[k];

                // leading z component

                fr = paz[k];

                g_zzzz_xxxx[k] = fr * g1_zzz_xxxx[k] + f2g * (3.0 * g20_zz_xxxx[k] - 3.0 * fgz * g21_zz_xxxx[k]);

                g_zzzz_xxxy[k] = fr * g1_zzz_xxxy[k] + f2g * (3.0 * g20_zz_xxxy[k] - 3.0 * fgz * g21_zz_xxxy[k]);

                g_zzzz_xxxz[k] = fr * g1_zzz_xxxz[k] + f2g * (3.0 * g20_zz_xxxz[k] - 3.0 * fgz * g21_zz_xxxz[k]) + f2t * gk_zzz_xxx[k];

                g_zzzz_xxyy[k] = fr * g1_zzz_xxyy[k] + f2g * (3.0 * g20_zz_xxyy[k] - 3.0 * fgz * g21_zz_xxyy[k]);

                g_zzzz_xxyz[k] = fr * g1_zzz_xxyz[k] + f2g * (3.0 * g20_zz_xxyz[k] - 3.0 * fgz * g21_zz_xxyz[k]) + f2t * gk_zzz_xxy[k];

                g_zzzz_xxzz[k] = fr * g1_zzz_xxzz[k] + f2g * (3.0 * g20_zz_xxzz[k] - 3.0 * fgz * g21_zz_xxzz[k]) + 2.0 * f2t * gk_zzz_xxz[k];

                g_zzzz_xyyy[k] = fr * g1_zzz_xyyy[k] + f2g * (3.0 * g20_zz_xyyy[k] - 3.0 * fgz * g21_zz_xyyy[k]);

                g_zzzz_xyyz[k] = fr * g1_zzz_xyyz[k] + f2g * (3.0 * g20_zz_xyyz[k] - 3.0 * fgz * g21_zz_xyyz[k]) + f2t * gk_zzz_xyy[k];

                g_zzzz_xyzz[k] = fr * g1_zzz_xyzz[k] + f2g * (3.0 * g20_zz_xyzz[k] - 3.0 * fgz * g21_zz_xyzz[k]) + 2.0 * f2t * gk_zzz_xyz[k];

                g_zzzz_xzzz[k] = fr * g1_zzz_xzzz[k] + f2g * (3.0 * g20_zz_xzzz[k] - 3.0 * fgz * g21_zz_xzzz[k]) + 3.0 * f2t * gk_zzz_xzz[k];

                g_zzzz_yyyy[k] = fr * g1_zzz_yyyy[k] + f2g * (3.0 * g20_zz_yyyy[k] - 3.0 * fgz * g21_zz_yyyy[k]);

                g_zzzz_yyyz[k] = fr * g1_zzz_yyyz[k] + f2g * (3.0 * g20_zz_yyyz[k] - 3.0 * fgz * g21_zz_yyyz[k]) + f2t * gk_zzz_yyy[k];

                g_zzzz_yyzz[k] = fr * g1_zzz_yyzz[k] + f2g * (3.0 * g20_zz_yyzz[k] - 3.0 * fgz * g21_zz_yyzz[k]) + 2.0 * f2t * gk_zzz_yyz[k];

                g_zzzz_yzzz[k] = fr * g1_zzz_yzzz[k] + f2g * (3.0 * g20_zz_yzzz[k] - 3.0 * fgz * g21_zz_yzzz[k]) + 3.0 * f2t * gk_zzz_yzz[k];

                g_zzzz_zzzz[k] = fr * g1_zzz_zzzz[k] + f2g * (3.0 * g20_zz_zzzz[k] - 3.0 * fgz * g21_zz_zzzz[k]) + 4.0 * f2t * gk_zzz_zzz[k];
            }

            idx++;
        }
    }
}

}  // namespace epotrecfunc
