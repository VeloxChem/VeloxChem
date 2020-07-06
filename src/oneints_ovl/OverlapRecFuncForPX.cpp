//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "OverlapRecFuncForPX.hpp"

namespace ovlrecfunc {  // ovlrecfunc namespace

void
compOverlapForPP(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_1_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx);

        auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1);

        auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2);

        auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx);

        // set up pointers to integrals

        auto ts_x_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx);

        auto ts_x_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 1);

        auto ts_x_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 2);

        auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3);

        auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4);

        auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5);

        auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6);

        auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7);

        auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_0_0, ts_0_x_0, ts_0_y_0, ts_0_z_0, ts_x_x_0, \
                                     ts_x_y_0, ts_x_z_0, ts_y_x_0, ts_y_y_0, ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_x_x_0[j] = pa_x[j] * ts_0_x_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

            ts_x_y_0[j] = pa_x[j] * ts_0_y_0[j];

            ts_x_z_0[j] = pa_x[j] * ts_0_z_0[j];

            ts_y_x_0[j] = pa_y[j] * ts_0_x_0[j];

            ts_y_y_0[j] = pa_y[j] * ts_0_y_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

            ts_y_z_0[j] = pa_y[j] * ts_0_z_0[j];

            ts_z_x_0[j] = pa_z[j] * ts_0_x_0[j];

            ts_z_y_0[j] = pa_z[j] * ts_0_y_0[j];

            ts_z_z_0[j] = pa_z[j] * ts_0_z_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];
        }

        idx++;
    }
}

void
compOverlapForPD(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_1_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx);

        auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1);

        auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2);

        // set up pointers to integrals

        auto ts_x_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx);

        auto ts_x_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 1);

        auto ts_x_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 2);

        auto ts_x_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 3);

        auto ts_x_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 4);

        auto ts_x_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 5);

        auto ts_y_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 6);

        auto ts_y_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 7);

        auto ts_y_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 8);

        auto ts_y_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 9);

        auto ts_y_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 10);

        auto ts_y_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 11);

        auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12);

        auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13);

        auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14);

        auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15);

        auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16);

        auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_x_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_y_0, \
                                     ts_0_yy_0, ts_0_yz_0, ts_0_z_0, ts_0_zz_0, ts_x_xx_0, ts_x_xy_0, ts_x_xz_0, \
                                     ts_x_yy_0, ts_x_yz_0, ts_x_zz_0, ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, \
                                     ts_y_yz_0, ts_y_zz_0, ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, \
                                     ts_z_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_x_xx_0[j] = pa_x[j] * ts_0_xx_0[j] + fl1_fx * ts_0_x_0[j];

            ts_x_xy_0[j] = pa_x[j] * ts_0_xy_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

            ts_x_xz_0[j] = pa_x[j] * ts_0_xz_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

            ts_x_yy_0[j] = pa_x[j] * ts_0_yy_0[j];

            ts_x_yz_0[j] = pa_x[j] * ts_0_yz_0[j];

            ts_x_zz_0[j] = pa_x[j] * ts_0_zz_0[j];

            ts_y_xx_0[j] = pa_y[j] * ts_0_xx_0[j];

            ts_y_xy_0[j] = pa_y[j] * ts_0_xy_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

            ts_y_xz_0[j] = pa_y[j] * ts_0_xz_0[j];

            ts_y_yy_0[j] = pa_y[j] * ts_0_yy_0[j] + fl1_fx * ts_0_y_0[j];

            ts_y_yz_0[j] = pa_y[j] * ts_0_yz_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

            ts_y_zz_0[j] = pa_y[j] * ts_0_zz_0[j];

            ts_z_xx_0[j] = pa_z[j] * ts_0_xx_0[j];

            ts_z_xy_0[j] = pa_z[j] * ts_0_xy_0[j];

            ts_z_xz_0[j] = pa_z[j] * ts_0_xz_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

            ts_z_yy_0[j] = pa_z[j] * ts_0_yy_0[j];

            ts_z_yz_0[j] = pa_z[j] * ts_0_yz_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

            ts_z_zz_0[j] = pa_z[j] * ts_0_zz_0[j] + fl1_fx * ts_0_z_0[j];
        }

        idx++;
    }
}

void
compOverlapForDP(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_2_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_x_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx);

        auto ts_x_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 1);

        auto ts_x_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 2);

        auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3);

        auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4);

        auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5);

        auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6);

        auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7);

        auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8);

        auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx);

        auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1);

        auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2);

        auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx);

        auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1);

        auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2);

        // set up pointers to integrals

        auto ts_xx_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx);

        auto ts_xx_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 1);

        auto ts_xx_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 2);

        auto ts_xy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 3);

        auto ts_xy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 4);

        auto ts_xy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 5);

        auto ts_xz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 6);

        auto ts_xz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 7);

        auto ts_xz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 8);

        auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9);

        auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10);

        auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11);

        auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12);

        auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13);

        auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14);

        auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15);

        auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16);

        auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_x_0, ts_0_y_0, ts_0_z_0, ts_x_0_0, ts_x_x_0, \
                                     ts_x_y_0, ts_x_z_0, ts_xx_x_0, ts_xx_y_0, ts_xx_z_0, ts_xy_x_0, ts_xy_y_0, \
                                     ts_xy_z_0, ts_xz_x_0, ts_xz_y_0, ts_xz_z_0, ts_y_0_0, ts_y_x_0, ts_y_y_0, ts_y_z_0, \
                                     ts_yy_x_0, ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, ts_z_0_0, \
                                     ts_z_x_0, ts_z_y_0, ts_z_z_0, ts_zz_x_0, ts_zz_y_0, ts_zz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_xx_x_0[j] = pa_x[j] * ts_x_x_0[j] + 0.5 * fl1_fx * ts_0_x_0[j] + 0.5 * fl1_fx * ts_x_0_0[j];

            ts_xx_y_0[j] = pa_x[j] * ts_x_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

            ts_xx_z_0[j] = pa_x[j] * ts_x_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

            ts_xy_x_0[j] = pa_x[j] * ts_y_x_0[j] + 0.5 * fl1_fx * ts_y_0_0[j];

            ts_xy_y_0[j] = pa_x[j] * ts_y_y_0[j];

            ts_xy_z_0[j] = pa_x[j] * ts_y_z_0[j];

            ts_xz_x_0[j] = pa_x[j] * ts_z_x_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

            ts_xz_y_0[j] = pa_x[j] * ts_z_y_0[j];

            ts_xz_z_0[j] = pa_x[j] * ts_z_z_0[j];

            ts_yy_x_0[j] = pa_y[j] * ts_y_x_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

            ts_yy_y_0[j] = pa_y[j] * ts_y_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j] + 0.5 * fl1_fx * ts_y_0_0[j];

            ts_yy_z_0[j] = pa_y[j] * ts_y_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

            ts_yz_x_0[j] = pa_y[j] * ts_z_x_0[j];

            ts_yz_y_0[j] = pa_y[j] * ts_z_y_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

            ts_yz_z_0[j] = pa_y[j] * ts_z_z_0[j];

            ts_zz_x_0[j] = pa_z[j] * ts_z_x_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

            ts_zz_y_0[j] = pa_z[j] * ts_z_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

            ts_zz_z_0[j] = pa_z[j] * ts_z_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];
        }

        idx++;
    }
}

void
compOverlapForPF(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_1_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx);

        auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1);

        auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2);

        auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3);

        auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4);

        auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5);

        auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6);

        auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7);

        auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8);

        auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9);

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto ts_x_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx);

        auto ts_x_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 1);

        auto ts_x_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 2);

        auto ts_x_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 3);

        auto ts_x_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 4);

        auto ts_x_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 5);

        auto ts_x_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 6);

        auto ts_x_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 7);

        auto ts_x_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 8);

        auto ts_x_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 9);

        auto ts_y_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 10);

        auto ts_y_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 11);

        auto ts_y_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 12);

        auto ts_y_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 13);

        auto ts_y_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 14);

        auto ts_y_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 15);

        auto ts_y_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 16);

        auto ts_y_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 17);

        auto ts_y_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 18);

        auto ts_y_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 19);

        auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20);

        auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21);

        auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22);

        auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23);

        auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24);

        auto ts_z_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 25);

        auto ts_z_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 26);

        auto ts_z_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 27);

        auto ts_z_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 28);

        auto ts_z_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 29);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_xx_0, ts_0_xxx_0, ts_0_xxy_0, ts_0_xxz_0, \
                                     ts_0_xy_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xz_0, ts_0_xzz_0, ts_0_yy_0, ts_0_yyy_0, \
                                     ts_0_yyz_0, ts_0_yz_0, ts_0_yzz_0, ts_0_zz_0, ts_0_zzz_0, ts_x_xxx_0, ts_x_xxy_0, \
                                     ts_x_xxz_0, ts_x_xyy_0, ts_x_xyz_0, ts_x_xzz_0, ts_x_yyy_0, ts_x_yyz_0, ts_x_yzz_0, \
                                     ts_x_zzz_0, ts_y_xxx_0, ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, \
                                     ts_y_yyy_0, ts_y_yyz_0, ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, \
                                     ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_x_xxx_0[j] = pa_x[j] * ts_0_xxx_0[j] + 1.5 * fl1_fx * ts_0_xx_0[j];

            ts_x_xxy_0[j] = pa_x[j] * ts_0_xxy_0[j] + fl1_fx * ts_0_xy_0[j];

            ts_x_xxz_0[j] = pa_x[j] * ts_0_xxz_0[j] + fl1_fx * ts_0_xz_0[j];

            ts_x_xyy_0[j] = pa_x[j] * ts_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

            ts_x_xyz_0[j] = pa_x[j] * ts_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

            ts_x_xzz_0[j] = pa_x[j] * ts_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

            ts_x_yyy_0[j] = pa_x[j] * ts_0_yyy_0[j];

            ts_x_yyz_0[j] = pa_x[j] * ts_0_yyz_0[j];

            ts_x_yzz_0[j] = pa_x[j] * ts_0_yzz_0[j];

            ts_x_zzz_0[j] = pa_x[j] * ts_0_zzz_0[j];

            ts_y_xxx_0[j] = pa_y[j] * ts_0_xxx_0[j];

            ts_y_xxy_0[j] = pa_y[j] * ts_0_xxy_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

            ts_y_xxz_0[j] = pa_y[j] * ts_0_xxz_0[j];

            ts_y_xyy_0[j] = pa_y[j] * ts_0_xyy_0[j] + fl1_fx * ts_0_xy_0[j];

            ts_y_xyz_0[j] = pa_y[j] * ts_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_xz_0[j];

            ts_y_xzz_0[j] = pa_y[j] * ts_0_xzz_0[j];

            ts_y_yyy_0[j] = pa_y[j] * ts_0_yyy_0[j] + 1.5 * fl1_fx * ts_0_yy_0[j];

            ts_y_yyz_0[j] = pa_y[j] * ts_0_yyz_0[j] + fl1_fx * ts_0_yz_0[j];

            ts_y_yzz_0[j] = pa_y[j] * ts_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

            ts_y_zzz_0[j] = pa_y[j] * ts_0_zzz_0[j];

            ts_z_xxx_0[j] = pa_z[j] * ts_0_xxx_0[j];

            ts_z_xxy_0[j] = pa_z[j] * ts_0_xxy_0[j];

            ts_z_xxz_0[j] = pa_z[j] * ts_0_xxz_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

            ts_z_xyy_0[j] = pa_z[j] * ts_0_xyy_0[j];

            ts_z_xyz_0[j] = pa_z[j] * ts_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_xy_0[j];

            ts_z_xzz_0[j] = pa_z[j] * ts_0_xzz_0[j] + fl1_fx * ts_0_xz_0[j];

            ts_z_yyy_0[j] = pa_z[j] * ts_0_yyy_0[j];

            ts_z_yyz_0[j] = pa_z[j] * ts_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

            ts_z_yzz_0[j] = pa_z[j] * ts_0_yzz_0[j] + fl1_fx * ts_0_yz_0[j];

            ts_z_zzz_0[j] = pa_z[j] * ts_0_zzz_0[j] + 1.5 * fl1_fx * ts_0_zz_0[j];
        }

        idx++;
    }
}

void
compOverlapForFP(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_3_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_xx_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx);

        auto ts_xx_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 1);

        auto ts_xx_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 2);

        auto ts_xy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 3);

        auto ts_xy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 4);

        auto ts_xy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 5);

        auto ts_xz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 6);

        auto ts_xz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 7);

        auto ts_xz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 8);

        auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9);

        auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10);

        auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11);

        auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12);

        auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13);

        auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14);

        auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15);

        auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16);

        auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17);

        auto ts_x_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx);

        auto ts_x_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 1);

        auto ts_x_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 2);

        auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3);

        auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4);

        auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5);

        auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6);

        auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7);

        auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8);

        auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx);

        auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1);

        auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2);

        auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3);

        auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4);

        auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto ts_xxx_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx);

        auto ts_xxx_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 1);

        auto ts_xxx_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 2);

        auto ts_xxy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 3);

        auto ts_xxy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 4);

        auto ts_xxy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 5);

        auto ts_xxz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 6);

        auto ts_xxz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 7);

        auto ts_xxz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 8);

        auto ts_xyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 9);

        auto ts_xyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 10);

        auto ts_xyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 11);

        auto ts_xyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 12);

        auto ts_xyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 13);

        auto ts_xyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 14);

        auto ts_xzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 15);

        auto ts_xzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 16);

        auto ts_xzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 17);

        auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18);

        auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19);

        auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20);

        auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21);

        auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22);

        auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23);

        auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24);

        auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25);

        auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26);

        auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27);

        auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28);

        auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_x_x_0, ts_x_y_0, ts_x_z_0, ts_xx_0_0, ts_xx_x_0, \
                                     ts_xx_y_0, ts_xx_z_0, ts_xxx_x_0, ts_xxx_y_0, ts_xxx_z_0, ts_xxy_x_0, ts_xxy_y_0, \
                                     ts_xxy_z_0, ts_xxz_x_0, ts_xxz_y_0, ts_xxz_z_0, ts_xy_0_0, ts_xy_x_0, ts_xy_y_0, \
                                     ts_xy_z_0, ts_xyy_x_0, ts_xyy_y_0, ts_xyy_z_0, ts_xyz_x_0, ts_xyz_y_0, ts_xyz_z_0, \
                                     ts_xz_0_0, ts_xz_x_0, ts_xz_y_0, ts_xz_z_0, ts_xzz_x_0, ts_xzz_y_0, ts_xzz_z_0, \
                                     ts_y_x_0, ts_y_y_0, ts_y_z_0, ts_yy_0_0, ts_yy_x_0, ts_yy_y_0, ts_yy_z_0, \
                                     ts_yyy_x_0, ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, ts_yyz_y_0, ts_yyz_z_0, ts_yz_0_0, \
                                     ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, ts_yzz_x_0, ts_yzz_y_0, ts_yzz_z_0, ts_z_x_0, \
                                     ts_z_y_0, ts_z_z_0, ts_zz_0_0, ts_zz_x_0, ts_zz_y_0, ts_zz_z_0, ts_zzz_x_0, \
                                     ts_zzz_y_0, ts_zzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_xxx_x_0[j] = pa_x[j] * ts_xx_x_0[j] + fl1_fx * ts_x_x_0[j] + 0.5 * fl1_fx * ts_xx_0_0[j];

            ts_xxx_y_0[j] = pa_x[j] * ts_xx_y_0[j] + fl1_fx * ts_x_y_0[j];

            ts_xxx_z_0[j] = pa_x[j] * ts_xx_z_0[j] + fl1_fx * ts_x_z_0[j];

            ts_xxy_x_0[j] = pa_x[j] * ts_xy_x_0[j] + 0.5 * fl1_fx * ts_y_x_0[j] + 0.5 * fl1_fx * ts_xy_0_0[j];

            ts_xxy_y_0[j] = pa_x[j] * ts_xy_y_0[j] + 0.5 * fl1_fx * ts_y_y_0[j];

            ts_xxy_z_0[j] = pa_x[j] * ts_xy_z_0[j] + 0.5 * fl1_fx * ts_y_z_0[j];

            ts_xxz_x_0[j] = pa_x[j] * ts_xz_x_0[j] + 0.5 * fl1_fx * ts_z_x_0[j] + 0.5 * fl1_fx * ts_xz_0_0[j];

            ts_xxz_y_0[j] = pa_x[j] * ts_xz_y_0[j] + 0.5 * fl1_fx * ts_z_y_0[j];

            ts_xxz_z_0[j] = pa_x[j] * ts_xz_z_0[j] + 0.5 * fl1_fx * ts_z_z_0[j];

            ts_xyy_x_0[j] = pa_x[j] * ts_yy_x_0[j] + 0.5 * fl1_fx * ts_yy_0_0[j];

            ts_xyy_y_0[j] = pa_x[j] * ts_yy_y_0[j];

            ts_xyy_z_0[j] = pa_x[j] * ts_yy_z_0[j];

            ts_xyz_x_0[j] = pa_x[j] * ts_yz_x_0[j] + 0.5 * fl1_fx * ts_yz_0_0[j];

            ts_xyz_y_0[j] = pa_x[j] * ts_yz_y_0[j];

            ts_xyz_z_0[j] = pa_x[j] * ts_yz_z_0[j];

            ts_xzz_x_0[j] = pa_x[j] * ts_zz_x_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

            ts_xzz_y_0[j] = pa_x[j] * ts_zz_y_0[j];

            ts_xzz_z_0[j] = pa_x[j] * ts_zz_z_0[j];

            ts_yyy_x_0[j] = pa_y[j] * ts_yy_x_0[j] + fl1_fx * ts_y_x_0[j];

            ts_yyy_y_0[j] = pa_y[j] * ts_yy_y_0[j] + fl1_fx * ts_y_y_0[j] + 0.5 * fl1_fx * ts_yy_0_0[j];

            ts_yyy_z_0[j] = pa_y[j] * ts_yy_z_0[j] + fl1_fx * ts_y_z_0[j];

            ts_yyz_x_0[j] = pa_y[j] * ts_yz_x_0[j] + 0.5 * fl1_fx * ts_z_x_0[j];

            ts_yyz_y_0[j] = pa_y[j] * ts_yz_y_0[j] + 0.5 * fl1_fx * ts_z_y_0[j] + 0.5 * fl1_fx * ts_yz_0_0[j];

            ts_yyz_z_0[j] = pa_y[j] * ts_yz_z_0[j] + 0.5 * fl1_fx * ts_z_z_0[j];

            ts_yzz_x_0[j] = pa_y[j] * ts_zz_x_0[j];

            ts_yzz_y_0[j] = pa_y[j] * ts_zz_y_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

            ts_yzz_z_0[j] = pa_y[j] * ts_zz_z_0[j];

            ts_zzz_x_0[j] = pa_z[j] * ts_zz_x_0[j] + fl1_fx * ts_z_x_0[j];

            ts_zzz_y_0[j] = pa_z[j] * ts_zz_y_0[j] + fl1_fx * ts_z_y_0[j];

            ts_zzz_z_0[j] = pa_z[j] * ts_zz_z_0[j] + fl1_fx * ts_z_z_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];
        }

        idx++;
    }
}

void
compOverlapForPG(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx);

        auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1);

        auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2);

        auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3);

        auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4);

        auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5);

        auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6);

        auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7);

        auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8);

        auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9);

        auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10);

        auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11);

        auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12);

        auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13);

        auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14);

        auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx);

        auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1);

        auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2);

        auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3);

        auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4);

        auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5);

        auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6);

        auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7);

        auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8);

        auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9);

        // set up pointers to integrals

        auto ts_x_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx);

        auto ts_x_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 1);

        auto ts_x_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 2);

        auto ts_x_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 3);

        auto ts_x_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 4);

        auto ts_x_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 5);

        auto ts_x_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 6);

        auto ts_x_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 7);

        auto ts_x_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 8);

        auto ts_x_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 9);

        auto ts_x_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 10);

        auto ts_x_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 11);

        auto ts_x_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 12);

        auto ts_x_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 13);

        auto ts_x_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 14);

        auto ts_y_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 15);

        auto ts_y_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 16);

        auto ts_y_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 17);

        auto ts_y_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 18);

        auto ts_y_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 19);

        auto ts_y_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 20);

        auto ts_y_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 21);

        auto ts_y_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 22);

        auto ts_y_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 23);

        auto ts_y_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 24);

        auto ts_y_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 25);

        auto ts_y_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 26);

        auto ts_y_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 27);

        auto ts_y_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 28);

        auto ts_y_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 29);

        auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30);

        auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31);

        auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32);

        auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33);

        auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34);

        auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35);

        auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36);

        auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37);

        auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38);

        auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39);

        auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40);

        auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41);

        auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42);

        auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43);

        auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_xxx_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, \
                                     ts_0_xxy_0, ts_0_xxyy_0, ts_0_xxyz_0, ts_0_xxz_0, ts_0_xxzz_0, ts_0_xyy_0, \
                                     ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyz_0, ts_0_xyzz_0, ts_0_xzz_0, ts_0_xzzz_0, \
                                     ts_0_yyy_0, ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyz_0, ts_0_yyzz_0, ts_0_yzz_0, \
                                     ts_0_yzzz_0, ts_0_zzz_0, ts_0_zzzz_0, ts_x_xxxx_0, ts_x_xxxy_0, ts_x_xxxz_0, \
                                     ts_x_xxyy_0, ts_x_xxyz_0, ts_x_xxzz_0, ts_x_xyyy_0, ts_x_xyyz_0, ts_x_xyzz_0, \
                                     ts_x_xzzz_0, ts_x_yyyy_0, ts_x_yyyz_0, ts_x_yyzz_0, ts_x_yzzz_0, ts_x_zzzz_0, \
                                     ts_y_xxxx_0, ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, \
                                     ts_y_xyyy_0, ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, \
                                     ts_y_yyzz_0, ts_y_yzzz_0, ts_y_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, \
                                     ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, \
                                     ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_x_xxxx_0[j] = pa_x[j] * ts_0_xxxx_0[j] + 2.0 * fl1_fx * ts_0_xxx_0[j];

            ts_x_xxxy_0[j] = pa_x[j] * ts_0_xxxy_0[j] + 1.5 * fl1_fx * ts_0_xxy_0[j];

            ts_x_xxxz_0[j] = pa_x[j] * ts_0_xxxz_0[j] + 1.5 * fl1_fx * ts_0_xxz_0[j];

            ts_x_xxyy_0[j] = pa_x[j] * ts_0_xxyy_0[j] + fl1_fx * ts_0_xyy_0[j];

            ts_x_xxyz_0[j] = pa_x[j] * ts_0_xxyz_0[j] + fl1_fx * ts_0_xyz_0[j];

            ts_x_xxzz_0[j] = pa_x[j] * ts_0_xxzz_0[j] + fl1_fx * ts_0_xzz_0[j];

            ts_x_xyyy_0[j] = pa_x[j] * ts_0_xyyy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

            ts_x_xyyz_0[j] = pa_x[j] * ts_0_xyyz_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

            ts_x_xyzz_0[j] = pa_x[j] * ts_0_xyzz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

            ts_x_xzzz_0[j] = pa_x[j] * ts_0_xzzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

            ts_x_yyyy_0[j] = pa_x[j] * ts_0_yyyy_0[j];

            ts_x_yyyz_0[j] = pa_x[j] * ts_0_yyyz_0[j];

            ts_x_yyzz_0[j] = pa_x[j] * ts_0_yyzz_0[j];

            ts_x_yzzz_0[j] = pa_x[j] * ts_0_yzzz_0[j];

            ts_x_zzzz_0[j] = pa_x[j] * ts_0_zzzz_0[j];

            ts_y_xxxx_0[j] = pa_y[j] * ts_0_xxxx_0[j];

            ts_y_xxxy_0[j] = pa_y[j] * ts_0_xxxy_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

            ts_y_xxxz_0[j] = pa_y[j] * ts_0_xxxz_0[j];

            ts_y_xxyy_0[j] = pa_y[j] * ts_0_xxyy_0[j] + fl1_fx * ts_0_xxy_0[j];

            ts_y_xxyz_0[j] = pa_y[j] * ts_0_xxyz_0[j] + 0.5 * fl1_fx * ts_0_xxz_0[j];

            ts_y_xxzz_0[j] = pa_y[j] * ts_0_xxzz_0[j];

            ts_y_xyyy_0[j] = pa_y[j] * ts_0_xyyy_0[j] + 1.5 * fl1_fx * ts_0_xyy_0[j];

            ts_y_xyyz_0[j] = pa_y[j] * ts_0_xyyz_0[j] + fl1_fx * ts_0_xyz_0[j];

            ts_y_xyzz_0[j] = pa_y[j] * ts_0_xyzz_0[j] + 0.5 * fl1_fx * ts_0_xzz_0[j];

            ts_y_xzzz_0[j] = pa_y[j] * ts_0_xzzz_0[j];

            ts_y_yyyy_0[j] = pa_y[j] * ts_0_yyyy_0[j] + 2.0 * fl1_fx * ts_0_yyy_0[j];

            ts_y_yyyz_0[j] = pa_y[j] * ts_0_yyyz_0[j] + 1.5 * fl1_fx * ts_0_yyz_0[j];

            ts_y_yyzz_0[j] = pa_y[j] * ts_0_yyzz_0[j] + fl1_fx * ts_0_yzz_0[j];

            ts_y_yzzz_0[j] = pa_y[j] * ts_0_yzzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

            ts_y_zzzz_0[j] = pa_y[j] * ts_0_zzzz_0[j];

            ts_z_xxxx_0[j] = pa_z[j] * ts_0_xxxx_0[j];

            ts_z_xxxy_0[j] = pa_z[j] * ts_0_xxxy_0[j];

            ts_z_xxxz_0[j] = pa_z[j] * ts_0_xxxz_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

            ts_z_xxyy_0[j] = pa_z[j] * ts_0_xxyy_0[j];

            ts_z_xxyz_0[j] = pa_z[j] * ts_0_xxyz_0[j] + 0.5 * fl1_fx * ts_0_xxy_0[j];

            ts_z_xxzz_0[j] = pa_z[j] * ts_0_xxzz_0[j] + fl1_fx * ts_0_xxz_0[j];

            ts_z_xyyy_0[j] = pa_z[j] * ts_0_xyyy_0[j];

            ts_z_xyyz_0[j] = pa_z[j] * ts_0_xyyz_0[j] + 0.5 * fl1_fx * ts_0_xyy_0[j];

            ts_z_xyzz_0[j] = pa_z[j] * ts_0_xyzz_0[j] + fl1_fx * ts_0_xyz_0[j];

            ts_z_xzzz_0[j] = pa_z[j] * ts_0_xzzz_0[j] + 1.5 * fl1_fx * ts_0_xzz_0[j];

            ts_z_yyyy_0[j] = pa_z[j] * ts_0_yyyy_0[j];

            ts_z_yyyz_0[j] = pa_z[j] * ts_0_yyyz_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

            ts_z_yyzz_0[j] = pa_z[j] * ts_0_yyzz_0[j] + fl1_fx * ts_0_yyz_0[j];

            ts_z_yzzz_0[j] = pa_z[j] * ts_0_yzzz_0[j] + 1.5 * fl1_fx * ts_0_yzz_0[j];

            ts_z_zzzz_0[j] = pa_z[j] * ts_0_zzzz_0[j] + 2.0 * fl1_fx * ts_0_zzz_0[j];
        }

        idx++;
    }
}

void
compOverlapForGP(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_xxx_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx);

        auto ts_xxx_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 1);

        auto ts_xxx_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 2);

        auto ts_xxy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 3);

        auto ts_xxy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 4);

        auto ts_xxy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 5);

        auto ts_xxz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 6);

        auto ts_xxz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 7);

        auto ts_xxz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 8);

        auto ts_xyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 9);

        auto ts_xyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 10);

        auto ts_xyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 11);

        auto ts_xyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 12);

        auto ts_xyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 13);

        auto ts_xyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 14);

        auto ts_xzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 15);

        auto ts_xzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 16);

        auto ts_xzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 17);

        auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18);

        auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19);

        auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20);

        auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21);

        auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22);

        auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23);

        auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24);

        auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25);

        auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26);

        auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27);

        auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28);

        auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29);

        auto ts_xx_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx);

        auto ts_xx_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 1);

        auto ts_xx_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 2);

        auto ts_xy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 3);

        auto ts_xy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 4);

        auto ts_xy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 5);

        auto ts_xz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 6);

        auto ts_xz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 7);

        auto ts_xz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 8);

        auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9);

        auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10);

        auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11);

        auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12);

        auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13);

        auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14);

        auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15);

        auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16);

        auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17);

        auto ts_xxx_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx);

        auto ts_xxy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 1);

        auto ts_xxz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 2);

        auto ts_xyy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 3);

        auto ts_xyz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 4);

        auto ts_xzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 5);

        auto ts_yyy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 6);

        auto ts_yyz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 7);

        auto ts_yzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 8);

        auto ts_zzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 9);

        // set up pointers to integrals

        auto ts_xxxx_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx);

        auto ts_xxxx_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 1);

        auto ts_xxxx_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 2);

        auto ts_xxxy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 3);

        auto ts_xxxy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 4);

        auto ts_xxxy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 5);

        auto ts_xxxz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 6);

        auto ts_xxxz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 7);

        auto ts_xxxz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 8);

        auto ts_xxyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 9);

        auto ts_xxyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 10);

        auto ts_xxyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 11);

        auto ts_xxyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 12);

        auto ts_xxyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 13);

        auto ts_xxyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 14);

        auto ts_xxzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 15);

        auto ts_xxzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 16);

        auto ts_xxzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 17);

        auto ts_xyyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 18);

        auto ts_xyyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 19);

        auto ts_xyyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 20);

        auto ts_xyyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 21);

        auto ts_xyyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 22);

        auto ts_xyyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 23);

        auto ts_xyzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 24);

        auto ts_xyzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 25);

        auto ts_xyzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 26);

        auto ts_xzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 27);

        auto ts_xzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 28);

        auto ts_xzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 29);

        auto ts_yyyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 30);

        auto ts_yyyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 31);

        auto ts_yyyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 32);

        auto ts_yyyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 33);

        auto ts_yyyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 34);

        auto ts_yyyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 35);

        auto ts_yyzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 36);

        auto ts_yyzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 37);

        auto ts_yyzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 38);

        auto ts_yzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 39);

        auto ts_yzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 40);

        auto ts_yzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 41);

        auto ts_zzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 42);

        auto ts_zzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 43);

        auto ts_zzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 44);

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_xx_x_0, ts_xx_y_0, ts_xx_z_0, ts_xxx_0_0, \
                                     ts_xxx_x_0, ts_xxx_y_0, ts_xxx_z_0, ts_xxxx_x_0, ts_xxxx_y_0, ts_xxxx_z_0, \
                                     ts_xxxy_x_0, ts_xxxy_y_0, ts_xxxy_z_0, ts_xxxz_x_0, ts_xxxz_y_0, ts_xxxz_z_0, \
                                     ts_xxy_0_0, ts_xxy_x_0, ts_xxy_y_0, ts_xxy_z_0, ts_xxyy_x_0, ts_xxyy_y_0, \
                                     ts_xxyy_z_0, ts_xxyz_x_0, ts_xxyz_y_0, ts_xxyz_z_0, ts_xxz_0_0, ts_xxz_x_0, \
                                     ts_xxz_y_0, ts_xxz_z_0, ts_xxzz_x_0, ts_xxzz_y_0, ts_xxzz_z_0, ts_xy_x_0, \
                                     ts_xy_y_0, ts_xy_z_0, ts_xyy_0_0, ts_xyy_x_0, ts_xyy_y_0, ts_xyy_z_0, ts_xyyy_x_0, \
                                     ts_xyyy_y_0, ts_xyyy_z_0, ts_xyyz_x_0, ts_xyyz_y_0, ts_xyyz_z_0, ts_xyz_0_0, \
                                     ts_xyz_x_0, ts_xyz_y_0, ts_xyz_z_0, ts_xyzz_x_0, ts_xyzz_y_0, ts_xyzz_z_0, \
                                     ts_xz_x_0, ts_xz_y_0, ts_xz_z_0, ts_xzz_0_0, ts_xzz_x_0, ts_xzz_y_0, ts_xzz_z_0, \
                                     ts_xzzz_x_0, ts_xzzz_y_0, ts_xzzz_z_0, ts_yy_x_0, ts_yy_y_0, ts_yy_z_0, ts_yyy_0_0, \
                                     ts_yyy_x_0, ts_yyy_y_0, ts_yyy_z_0, ts_yyyy_x_0, ts_yyyy_y_0, ts_yyyy_z_0, \
                                     ts_yyyz_x_0, ts_yyyz_y_0, ts_yyyz_z_0, ts_yyz_0_0, ts_yyz_x_0, ts_yyz_y_0, \
                                     ts_yyz_z_0, ts_yyzz_x_0, ts_yyzz_y_0, ts_yyzz_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, \
                                     ts_yzz_0_0, ts_yzz_x_0, ts_yzz_y_0, ts_yzz_z_0, ts_yzzz_x_0, ts_yzzz_y_0, \
                                     ts_yzzz_z_0, ts_zz_x_0, ts_zz_y_0, ts_zz_z_0, ts_zzz_0_0, ts_zzz_x_0, ts_zzz_y_0, \
                                     ts_zzz_z_0, ts_zzzz_x_0, ts_zzzz_y_0, ts_zzzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_xxxx_x_0[j] = pa_x[j] * ts_xxx_x_0[j] + 1.5 * fl1_fx * ts_xx_x_0[j] + 0.5 * fl1_fx * ts_xxx_0_0[j];

            ts_xxxx_y_0[j] = pa_x[j] * ts_xxx_y_0[j] + 1.5 * fl1_fx * ts_xx_y_0[j];

            ts_xxxx_z_0[j] = pa_x[j] * ts_xxx_z_0[j] + 1.5 * fl1_fx * ts_xx_z_0[j];

            ts_xxxy_x_0[j] = pa_x[j] * ts_xxy_x_0[j] + fl1_fx * ts_xy_x_0[j] + 0.5 * fl1_fx * ts_xxy_0_0[j];

            ts_xxxy_y_0[j] = pa_x[j] * ts_xxy_y_0[j] + fl1_fx * ts_xy_y_0[j];

            ts_xxxy_z_0[j] = pa_x[j] * ts_xxy_z_0[j] + fl1_fx * ts_xy_z_0[j];

            ts_xxxz_x_0[j] = pa_x[j] * ts_xxz_x_0[j] + fl1_fx * ts_xz_x_0[j] + 0.5 * fl1_fx * ts_xxz_0_0[j];

            ts_xxxz_y_0[j] = pa_x[j] * ts_xxz_y_0[j] + fl1_fx * ts_xz_y_0[j];

            ts_xxxz_z_0[j] = pa_x[j] * ts_xxz_z_0[j] + fl1_fx * ts_xz_z_0[j];

            ts_xxyy_x_0[j] = pa_x[j] * ts_xyy_x_0[j] + 0.5 * fl1_fx * ts_yy_x_0[j] + 0.5 * fl1_fx * ts_xyy_0_0[j];

            ts_xxyy_y_0[j] = pa_x[j] * ts_xyy_y_0[j] + 0.5 * fl1_fx * ts_yy_y_0[j];

            ts_xxyy_z_0[j] = pa_x[j] * ts_xyy_z_0[j] + 0.5 * fl1_fx * ts_yy_z_0[j];

            ts_xxyz_x_0[j] = pa_x[j] * ts_xyz_x_0[j] + 0.5 * fl1_fx * ts_yz_x_0[j] + 0.5 * fl1_fx * ts_xyz_0_0[j];

            ts_xxyz_y_0[j] = pa_x[j] * ts_xyz_y_0[j] + 0.5 * fl1_fx * ts_yz_y_0[j];

            ts_xxyz_z_0[j] = pa_x[j] * ts_xyz_z_0[j] + 0.5 * fl1_fx * ts_yz_z_0[j];

            ts_xxzz_x_0[j] = pa_x[j] * ts_xzz_x_0[j] + 0.5 * fl1_fx * ts_zz_x_0[j] + 0.5 * fl1_fx * ts_xzz_0_0[j];

            ts_xxzz_y_0[j] = pa_x[j] * ts_xzz_y_0[j] + 0.5 * fl1_fx * ts_zz_y_0[j];

            ts_xxzz_z_0[j] = pa_x[j] * ts_xzz_z_0[j] + 0.5 * fl1_fx * ts_zz_z_0[j];

            ts_xyyy_x_0[j] = pa_x[j] * ts_yyy_x_0[j] + 0.5 * fl1_fx * ts_yyy_0_0[j];

            ts_xyyy_y_0[j] = pa_x[j] * ts_yyy_y_0[j];

            ts_xyyy_z_0[j] = pa_x[j] * ts_yyy_z_0[j];

            ts_xyyz_x_0[j] = pa_x[j] * ts_yyz_x_0[j] + 0.5 * fl1_fx * ts_yyz_0_0[j];

            ts_xyyz_y_0[j] = pa_x[j] * ts_yyz_y_0[j];

            ts_xyyz_z_0[j] = pa_x[j] * ts_yyz_z_0[j];

            ts_xyzz_x_0[j] = pa_x[j] * ts_yzz_x_0[j] + 0.5 * fl1_fx * ts_yzz_0_0[j];

            ts_xyzz_y_0[j] = pa_x[j] * ts_yzz_y_0[j];

            ts_xyzz_z_0[j] = pa_x[j] * ts_yzz_z_0[j];

            ts_xzzz_x_0[j] = pa_x[j] * ts_zzz_x_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];

            ts_xzzz_y_0[j] = pa_x[j] * ts_zzz_y_0[j];

            ts_xzzz_z_0[j] = pa_x[j] * ts_zzz_z_0[j];

            ts_yyyy_x_0[j] = pa_y[j] * ts_yyy_x_0[j] + 1.5 * fl1_fx * ts_yy_x_0[j];

            ts_yyyy_y_0[j] = pa_y[j] * ts_yyy_y_0[j] + 1.5 * fl1_fx * ts_yy_y_0[j] + 0.5 * fl1_fx * ts_yyy_0_0[j];

            ts_yyyy_z_0[j] = pa_y[j] * ts_yyy_z_0[j] + 1.5 * fl1_fx * ts_yy_z_0[j];

            ts_yyyz_x_0[j] = pa_y[j] * ts_yyz_x_0[j] + fl1_fx * ts_yz_x_0[j];

            ts_yyyz_y_0[j] = pa_y[j] * ts_yyz_y_0[j] + fl1_fx * ts_yz_y_0[j] + 0.5 * fl1_fx * ts_yyz_0_0[j];

            ts_yyyz_z_0[j] = pa_y[j] * ts_yyz_z_0[j] + fl1_fx * ts_yz_z_0[j];

            ts_yyzz_x_0[j] = pa_y[j] * ts_yzz_x_0[j] + 0.5 * fl1_fx * ts_zz_x_0[j];

            ts_yyzz_y_0[j] = pa_y[j] * ts_yzz_y_0[j] + 0.5 * fl1_fx * ts_zz_y_0[j] + 0.5 * fl1_fx * ts_yzz_0_0[j];

            ts_yyzz_z_0[j] = pa_y[j] * ts_yzz_z_0[j] + 0.5 * fl1_fx * ts_zz_z_0[j];

            ts_yzzz_x_0[j] = pa_y[j] * ts_zzz_x_0[j];

            ts_yzzz_y_0[j] = pa_y[j] * ts_zzz_y_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];

            ts_yzzz_z_0[j] = pa_y[j] * ts_zzz_z_0[j];

            ts_zzzz_x_0[j] = pa_z[j] * ts_zzz_x_0[j] + 1.5 * fl1_fx * ts_zz_x_0[j];

            ts_zzzz_y_0[j] = pa_z[j] * ts_zzz_y_0[j] + 1.5 * fl1_fx * ts_zz_y_0[j];

            ts_zzzz_z_0[j] = pa_z[j] * ts_zzz_z_0[j] + 1.5 * fl1_fx * ts_zz_z_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];
        }

        idx++;
    }
}

}  // namespace ovlrecfunc
