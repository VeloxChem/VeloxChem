//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AngularMomentumRecFuncForPX.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForPP(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tlx_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + idx);

        auto tly_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + bdim + idx);

        auto tlz_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + 2 * bdim + idx);

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx);

        auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx);

        auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx);

        auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1);

        auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2);

        auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2);

        // set up pointers to integrals

        auto tlx_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx);

        auto tly_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx);

        auto tlz_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx);

        auto tlx_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 1);

        auto tly_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tlz_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tlx_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 2);

        auto tly_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tlz_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3);

        auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4);

        auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5);

        auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6);

        auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7);

        auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8);

        auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_0_x_0, tdx_0_y_0, tdx_0_z_0, tdy_0_x_0, \
                                     tdy_0_y_0, tdy_0_z_0, tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, tlx_0_0_0, tlx_0_x_0, \
                                     tlx_0_y_0, tlx_0_z_0, tlx_x_x_0, tlx_x_y_0, tlx_x_z_0, tlx_y_x_0, tlx_y_y_0, \
                                     tlx_y_z_0, tlx_z_x_0, tlx_z_y_0, tlx_z_z_0, tly_0_0_0, tly_0_x_0, tly_0_y_0, \
                                     tly_0_z_0, tly_x_x_0, tly_x_y_0, tly_x_z_0, tly_y_x_0, tly_y_y_0, tly_y_z_0, \
                                     tly_z_x_0, tly_z_y_0, tly_z_z_0, tlz_0_0_0, tlz_0_x_0, tlz_0_y_0, tlz_0_z_0, \
                                     tlz_x_x_0, tlz_x_y_0, tlz_x_z_0, tlz_y_x_0, tlz_y_y_0, tlz_y_z_0, tlz_z_x_0, \
                                     tlz_z_y_0, tlz_z_z_0, tpx_0_x_0, tpx_0_y_0, tpx_0_z_0, tpy_0_x_0, tpy_0_y_0, \
                                     tpy_0_z_0, tpz_0_x_0, tpz_0_y_0, tpz_0_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_x_x_0[j] = pa_x[j] * tlx_0_x_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j];

            tly_x_x_0[j] = pa_x[j] * tly_0_x_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j] + fl1_fx * fl1_fgb * tdz_0_x_0[j];

            tlz_x_x_0[j] = pa_x[j] * tlz_0_x_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] - 0.5 * fl1_fx * tpy_0_x_0[j] - fl1_fx * fl1_fgb * tdy_0_x_0[j];

            tlx_x_y_0[j] = pa_x[j] * tlx_0_y_0[j];

            tly_x_y_0[j] = pa_x[j] * tly_0_y_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j] + fl1_fx * fl1_fgb * tdz_0_y_0[j];

            tlz_x_y_0[j] = pa_x[j] * tlz_0_y_0[j] - 0.5 * fl1_fx * tpy_0_y_0[j] - fl1_fx * fl1_fgb * tdy_0_y_0[j];

            tlx_x_z_0[j] = pa_x[j] * tlx_0_z_0[j];

            tly_x_z_0[j] = pa_x[j] * tly_0_z_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j] + fl1_fx * fl1_fgb * tdz_0_z_0[j];

            tlz_x_z_0[j] = pa_x[j] * tlz_0_z_0[j] - 0.5 * fl1_fx * tpy_0_z_0[j] - fl1_fx * fl1_fgb * tdy_0_z_0[j];

            tlx_y_x_0[j] = pa_y[j] * tlx_0_x_0[j] - 0.5 * fl1_fx * tpz_0_x_0[j] - fl1_fx * fl1_fgb * tdz_0_x_0[j];

            tly_y_x_0[j] = pa_y[j] * tly_0_x_0[j];

            tlz_y_x_0[j] = pa_y[j] * tlz_0_x_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j] + fl1_fx * fl1_fgb * tdx_0_x_0[j];

            tlx_y_y_0[j] = pa_y[j] * tlx_0_y_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] - 0.5 * fl1_fx * tpz_0_y_0[j] - fl1_fx * fl1_fgb * tdz_0_y_0[j];

            tly_y_y_0[j] = pa_y[j] * tly_0_y_0[j] + 0.5 * fl1_fx * tly_0_0_0[j];

            tlz_y_y_0[j] = pa_y[j] * tlz_0_y_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] + fl1_fx * fl1_fgb * tdx_0_y_0[j];

            tlx_y_z_0[j] = pa_y[j] * tlx_0_z_0[j] - 0.5 * fl1_fx * tpz_0_z_0[j] - fl1_fx * fl1_fgb * tdz_0_z_0[j];

            tly_y_z_0[j] = pa_y[j] * tly_0_z_0[j];

            tlz_y_z_0[j] = pa_y[j] * tlz_0_z_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] + fl1_fx * fl1_fgb * tdx_0_z_0[j];

            tlx_z_x_0[j] = pa_z[j] * tlx_0_x_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j] + fl1_fx * fl1_fgb * tdy_0_x_0[j];

            tly_z_x_0[j] = pa_z[j] * tly_0_x_0[j] - 0.5 * fl1_fx * tpx_0_x_0[j] - fl1_fx * fl1_fgb * tdx_0_x_0[j];

            tlz_z_x_0[j] = pa_z[j] * tlz_0_x_0[j];

            tlx_z_y_0[j] = pa_z[j] * tlx_0_y_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j] + fl1_fx * fl1_fgb * tdy_0_y_0[j];

            tly_z_y_0[j] = pa_z[j] * tly_0_y_0[j] - 0.5 * fl1_fx * tpx_0_y_0[j] - fl1_fx * fl1_fgb * tdx_0_y_0[j];

            tlz_z_y_0[j] = pa_z[j] * tlz_0_y_0[j];

            tlx_z_z_0[j] = pa_z[j] * tlx_0_z_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] + fl1_fx * fl1_fgb * tdy_0_z_0[j];

            tly_z_z_0[j] = pa_z[j] * tly_0_z_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] - 0.5 * fl1_fx * tpx_0_z_0[j] - fl1_fx * fl1_fgb * tdx_0_z_0[j];

            tlz_z_z_0[j] = pa_z[j] * tlz_0_z_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPD(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForPD_0_27(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForPD_27_54(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForPD_0_27(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,27)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tdx_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx);

        auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx);

        auto tdz_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx);

        auto tdx_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 1);

        auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tdz_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tdx_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 2);

        auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tdz_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5);

        // set up pointers to integrals

        auto tlx_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx);

        auto tly_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx);

        auto tlz_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx);

        auto tlx_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 1);

        auto tly_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 1);

        auto tlz_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 1);

        auto tlx_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 2);

        auto tly_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 2);

        auto tlz_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 2);

        auto tlx_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 3);

        auto tly_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 3);

        auto tlz_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 3);

        auto tlx_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 4);

        auto tly_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 4);

        auto tlz_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 4);

        auto tlx_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 5);

        auto tly_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 5);

        auto tlz_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 5);

        auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6);

        auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7);

        auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8);

        auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8);

        // Batch of Integrals (0,27)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdy_0_xx_0, \
                                     tdy_0_xy_0, tdy_0_xz_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_zz_0, tdz_0_xx_0, tdz_0_xy_0, \
                                     tdz_0_xz_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_zz_0, tlx_0_x_0, tlx_0_xx_0, tlx_0_xy_0, \
                                     tlx_0_xz_0, tlx_0_y_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_z_0, tlx_0_zz_0, tlx_x_xx_0, \
                                     tlx_x_xy_0, tlx_x_xz_0, tlx_x_yy_0, tlx_x_yz_0, tlx_x_zz_0, tlx_y_xx_0, tlx_y_xy_0, \
                                     tlx_y_xz_0, tly_0_x_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_y_0, tly_0_yy_0, \
                                     tly_0_yz_0, tly_0_z_0, tly_0_zz_0, tly_x_xx_0, tly_x_xy_0, tly_x_xz_0, tly_x_yy_0, \
                                     tly_x_yz_0, tly_x_zz_0, tly_y_xx_0, tly_y_xy_0, tly_y_xz_0, tlz_0_x_0, tlz_0_xx_0, \
                                     tlz_0_xy_0, tlz_0_xz_0, tlz_0_y_0, tlz_0_yy_0, tlz_0_yz_0, tlz_0_z_0, tlz_0_zz_0, \
                                     tlz_x_xx_0, tlz_x_xy_0, tlz_x_xz_0, tlz_x_yy_0, tlz_x_yz_0, tlz_x_zz_0, tlz_y_xx_0, \
                                     tlz_y_xy_0, tlz_y_xz_0, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, tpy_0_xx_0, tpy_0_xy_0, \
                                     tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, tpz_0_xx_0, tpz_0_xy_0, tpz_0_xz_0, \
                                     tpz_0_yy_0, tpz_0_yz_0, tpz_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_x_xx_0[j] = pa_x[j] * tlx_0_xx_0[j] + fl1_fx * tlx_0_x_0[j];

            tly_x_xx_0[j] = pa_x[j] * tly_0_xx_0[j] + fl1_fx * tly_0_x_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j] + fl1_fx * fl1_fgb * tdz_0_xx_0[j];

            tlz_x_xx_0[j] = pa_x[j] * tlz_0_xx_0[j] + fl1_fx * tlz_0_x_0[j] - 0.5 * fl1_fx * tpy_0_xx_0[j] - fl1_fx * fl1_fgb * tdy_0_xx_0[j];

            tlx_x_xy_0[j] = pa_x[j] * tlx_0_xy_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j];

            tly_x_xy_0[j] = pa_x[j] * tly_0_xy_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] + fl1_fx * fl1_fgb * tdz_0_xy_0[j];

            tlz_x_xy_0[j] = pa_x[j] * tlz_0_xy_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j] - 0.5 * fl1_fx * tpy_0_xy_0[j] - fl1_fx * fl1_fgb * tdy_0_xy_0[j];

            tlx_x_xz_0[j] = pa_x[j] * tlx_0_xz_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j];

            tly_x_xz_0[j] = pa_x[j] * tly_0_xz_0[j] + 0.5 * fl1_fx * tly_0_z_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j] + fl1_fx * fl1_fgb * tdz_0_xz_0[j];

            tlz_x_xz_0[j] = pa_x[j] * tlz_0_xz_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] - 0.5 * fl1_fx * tpy_0_xz_0[j] - fl1_fx * fl1_fgb * tdy_0_xz_0[j];

            tlx_x_yy_0[j] = pa_x[j] * tlx_0_yy_0[j];

            tly_x_yy_0[j] = pa_x[j] * tly_0_yy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j] + fl1_fx * fl1_fgb * tdz_0_yy_0[j];

            tlz_x_yy_0[j] = pa_x[j] * tlz_0_yy_0[j] - 0.5 * fl1_fx * tpy_0_yy_0[j] - fl1_fx * fl1_fgb * tdy_0_yy_0[j];

            tlx_x_yz_0[j] = pa_x[j] * tlx_0_yz_0[j];

            tly_x_yz_0[j] = pa_x[j] * tly_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j] + fl1_fx * fl1_fgb * tdz_0_yz_0[j];

            tlz_x_yz_0[j] = pa_x[j] * tlz_0_yz_0[j] - 0.5 * fl1_fx * tpy_0_yz_0[j] - fl1_fx * fl1_fgb * tdy_0_yz_0[j];

            tlx_x_zz_0[j] = pa_x[j] * tlx_0_zz_0[j];

            tly_x_zz_0[j] = pa_x[j] * tly_0_zz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j] + fl1_fx * fl1_fgb * tdz_0_zz_0[j];

            tlz_x_zz_0[j] = pa_x[j] * tlz_0_zz_0[j] - 0.5 * fl1_fx * tpy_0_zz_0[j] - fl1_fx * fl1_fgb * tdy_0_zz_0[j];

            tlx_y_xx_0[j] = pa_y[j] * tlx_0_xx_0[j] - 0.5 * fl1_fx * tpz_0_xx_0[j] - fl1_fx * fl1_fgb * tdz_0_xx_0[j];

            tly_y_xx_0[j] = pa_y[j] * tly_0_xx_0[j];

            tlz_y_xx_0[j] = pa_y[j] * tlz_0_xx_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j] + fl1_fx * fl1_fgb * tdx_0_xx_0[j];

            tlx_y_xy_0[j] = pa_y[j] * tlx_0_xy_0[j] + 0.5 * fl1_fx * tlx_0_x_0[j] - 0.5 * fl1_fx * tpz_0_xy_0[j] - fl1_fx * fl1_fgb * tdz_0_xy_0[j];

            tly_y_xy_0[j] = pa_y[j] * tly_0_xy_0[j] + 0.5 * fl1_fx * tly_0_x_0[j];

            tlz_y_xy_0[j] = pa_y[j] * tlz_0_xy_0[j] + 0.5 * fl1_fx * tlz_0_x_0[j] + 0.5 * fl1_fx * tpx_0_xy_0[j] + fl1_fx * fl1_fgb * tdx_0_xy_0[j];

            tlx_y_xz_0[j] = pa_y[j] * tlx_0_xz_0[j] - 0.5 * fl1_fx * tpz_0_xz_0[j] - fl1_fx * fl1_fgb * tdz_0_xz_0[j];

            tly_y_xz_0[j] = pa_y[j] * tly_0_xz_0[j];

            tlz_y_xz_0[j] = pa_y[j] * tlz_0_xz_0[j] + 0.5 * fl1_fx * tpx_0_xz_0[j] + fl1_fx * fl1_fgb * tdx_0_xz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPD_27_54(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (27,54)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpx_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 3);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 4);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 5);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tdx_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx);

        auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx);

        auto tdx_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 1);

        auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tdx_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 2);

        auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tdx_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 3);

        auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tdx_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 4);

        auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tdx_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 5);

        auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5);

        // set up pointers to integrals

        auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9);

        auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10);

        auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11);

        auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12);

        auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13);

        auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14);

        auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15);

        auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16);

        auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17);

        auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17);

        // Batch of Integrals (27,54)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_yy_0, \
                                     tdx_0_yz_0, tdx_0_zz_0, tdy_0_xx_0, tdy_0_xy_0, tdy_0_xz_0, tdy_0_yy_0, tdy_0_yz_0, \
                                     tdy_0_zz_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_zz_0, tlx_0_x_0, tlx_0_xx_0, tlx_0_xy_0, \
                                     tlx_0_xz_0, tlx_0_y_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_z_0, tlx_0_zz_0, tlx_y_yy_0, \
                                     tlx_y_yz_0, tlx_y_zz_0, tlx_z_xx_0, tlx_z_xy_0, tlx_z_xz_0, tlx_z_yy_0, tlx_z_yz_0, \
                                     tlx_z_zz_0, tly_0_x_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_y_0, tly_0_yy_0, \
                                     tly_0_yz_0, tly_0_z_0, tly_0_zz_0, tly_y_yy_0, tly_y_yz_0, tly_y_zz_0, tly_z_xx_0, \
                                     tly_z_xy_0, tly_z_xz_0, tly_z_yy_0, tly_z_yz_0, tly_z_zz_0, tlz_0_x_0, tlz_0_xx_0, \
                                     tlz_0_xy_0, tlz_0_xz_0, tlz_0_y_0, tlz_0_yy_0, tlz_0_yz_0, tlz_0_z_0, tlz_0_zz_0, \
                                     tlz_y_yy_0, tlz_y_yz_0, tlz_y_zz_0, tlz_z_xx_0, tlz_z_xy_0, tlz_z_xz_0, tlz_z_yy_0, \
                                     tlz_z_yz_0, tlz_z_zz_0, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, tpx_0_yy_0, tpx_0_yz_0, \
                                     tpx_0_zz_0, tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, \
                                     tpz_0_yy_0, tpz_0_yz_0, tpz_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_y_yy_0[j] = pa_y[j] * tlx_0_yy_0[j] + fl1_fx * tlx_0_y_0[j] - 0.5 * fl1_fx * tpz_0_yy_0[j] - fl1_fx * fl1_fgb * tdz_0_yy_0[j];

            tly_y_yy_0[j] = pa_y[j] * tly_0_yy_0[j] + fl1_fx * tly_0_y_0[j];

            tlz_y_yy_0[j] = pa_y[j] * tlz_0_yy_0[j] + fl1_fx * tlz_0_y_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] + fl1_fx * fl1_fgb * tdx_0_yy_0[j];

            tlx_y_yz_0[j] = pa_y[j] * tlx_0_yz_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j] - 0.5 * fl1_fx * tpz_0_yz_0[j] - fl1_fx * fl1_fgb * tdz_0_yz_0[j];

            tly_y_yz_0[j] = pa_y[j] * tly_0_yz_0[j] + 0.5 * fl1_fx * tly_0_z_0[j];

            tlz_y_yz_0[j] = pa_y[j] * tlz_0_yz_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] + fl1_fx * fl1_fgb * tdx_0_yz_0[j];

            tlx_y_zz_0[j] = pa_y[j] * tlx_0_zz_0[j] - 0.5 * fl1_fx * tpz_0_zz_0[j] - fl1_fx * fl1_fgb * tdz_0_zz_0[j];

            tly_y_zz_0[j] = pa_y[j] * tly_0_zz_0[j];

            tlz_y_zz_0[j] = pa_y[j] * tlz_0_zz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] + fl1_fx * fl1_fgb * tdx_0_zz_0[j];

            tlx_z_xx_0[j] = pa_z[j] * tlx_0_xx_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j] + fl1_fx * fl1_fgb * tdy_0_xx_0[j];

            tly_z_xx_0[j] = pa_z[j] * tly_0_xx_0[j] - 0.5 * fl1_fx * tpx_0_xx_0[j] - fl1_fx * fl1_fgb * tdx_0_xx_0[j];

            tlz_z_xx_0[j] = pa_z[j] * tlz_0_xx_0[j];

            tlx_z_xy_0[j] = pa_z[j] * tlx_0_xy_0[j] + 0.5 * fl1_fx * tpy_0_xy_0[j] + fl1_fx * fl1_fgb * tdy_0_xy_0[j];

            tly_z_xy_0[j] = pa_z[j] * tly_0_xy_0[j] - 0.5 * fl1_fx * tpx_0_xy_0[j] - fl1_fx * fl1_fgb * tdx_0_xy_0[j];

            tlz_z_xy_0[j] = pa_z[j] * tlz_0_xy_0[j];

            tlx_z_xz_0[j] = pa_z[j] * tlx_0_xz_0[j] + 0.5 * fl1_fx * tlx_0_x_0[j] + 0.5 * fl1_fx * tpy_0_xz_0[j] + fl1_fx * fl1_fgb * tdy_0_xz_0[j];

            tly_z_xz_0[j] = pa_z[j] * tly_0_xz_0[j] + 0.5 * fl1_fx * tly_0_x_0[j] - 0.5 * fl1_fx * tpx_0_xz_0[j] - fl1_fx * fl1_fgb * tdx_0_xz_0[j];

            tlz_z_xz_0[j] = pa_z[j] * tlz_0_xz_0[j] + 0.5 * fl1_fx * tlz_0_x_0[j];

            tlx_z_yy_0[j] = pa_z[j] * tlx_0_yy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j] + fl1_fx * fl1_fgb * tdy_0_yy_0[j];

            tly_z_yy_0[j] = pa_z[j] * tly_0_yy_0[j] - 0.5 * fl1_fx * tpx_0_yy_0[j] - fl1_fx * fl1_fgb * tdx_0_yy_0[j];

            tlz_z_yy_0[j] = pa_z[j] * tlz_0_yy_0[j];

            tlx_z_yz_0[j] = pa_z[j] * tlx_0_yz_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j] + fl1_fx * fl1_fgb * tdy_0_yz_0[j];

            tly_z_yz_0[j] = pa_z[j] * tly_0_yz_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] - 0.5 * fl1_fx * tpx_0_yz_0[j] - fl1_fx * fl1_fgb * tdx_0_yz_0[j];

            tlz_z_yz_0[j] = pa_z[j] * tlz_0_yz_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j];

            tlx_z_zz_0[j] = pa_z[j] * tlx_0_zz_0[j] + fl1_fx * tlx_0_z_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] + fl1_fx * fl1_fgb * tdy_0_zz_0[j];

            tly_z_zz_0[j] = pa_z[j] * tly_0_zz_0[j] + fl1_fx * tly_0_z_0[j] - 0.5 * fl1_fx * tpx_0_zz_0[j] - fl1_fx * fl1_fgb * tdx_0_zz_0[j];

            tlz_z_zz_0[j] = pa_z[j] * tlz_0_zz_0[j] + fl1_fx * tlz_0_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForDP(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForDP_0_27(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForDP_27_54(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForDP_0_27(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,27)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_2_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tlx_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx);

        auto tly_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx);

        auto tlz_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx);

        auto tlx_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 1);

        auto tly_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tlz_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tlx_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 2);

        auto tly_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tlz_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3);

        auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4);

        auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5);

        auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6);

        auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7);

        auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8);

        auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tlx_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx);

        auto tly_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx);

        auto tlz_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx);

        auto tlx_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 1);

        auto tly_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 2);

        auto tly_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tpy_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx);

        auto tpz_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx);

        auto tpy_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tpz_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tpy_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tpz_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tdy_x_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx);

        auto tdz_x_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx);

        auto tdy_x_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tdz_x_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tdy_x_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tdz_x_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tdy_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tdz_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tdy_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tdz_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tdy_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tdz_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8);

        // set up pointers to integrals

        auto tlx_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx);

        auto tly_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx);

        auto tlz_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx);

        auto tlx_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 1);

        auto tly_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tlz_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tlx_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 2);

        auto tly_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tlz_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tlx_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 3);

        auto tly_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tlz_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tlx_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 4);

        auto tly_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tlz_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tlx_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 5);

        auto tly_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tlz_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tlx_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 6);

        auto tly_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tlz_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tlx_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 7);

        auto tly_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tlz_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tlx_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 8);

        auto tly_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tlz_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 8);

        // Batch of Integrals (0,27)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_x_x_0, tdy_x_y_0, tdy_x_z_0, tdy_y_x_0, tdy_y_y_0, \
                                     tdy_y_z_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, tdz_x_x_0, tdz_x_y_0, tdz_x_z_0, \
                                     tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, tlx_0_x_0, \
                                     tlx_0_y_0, tlx_0_z_0, tlx_x_0_0, tlx_x_x_0, tlx_x_y_0, tlx_x_z_0, tlx_xx_x_0, \
                                     tlx_xx_y_0, tlx_xx_z_0, tlx_xy_x_0, tlx_xy_y_0, tlx_xy_z_0, tlx_xz_x_0, tlx_xz_y_0, \
                                     tlx_xz_z_0, tlx_y_0_0, tlx_y_x_0, tlx_y_y_0, tlx_y_z_0, tlx_z_0_0, tlx_z_x_0, \
                                     tlx_z_y_0, tlx_z_z_0, tly_0_x_0, tly_0_y_0, tly_0_z_0, tly_x_0_0, tly_x_x_0, \
                                     tly_x_y_0, tly_x_z_0, tly_xx_x_0, tly_xx_y_0, tly_xx_z_0, tly_xy_x_0, tly_xy_y_0, \
                                     tly_xy_z_0, tly_xz_x_0, tly_xz_y_0, tly_xz_z_0, tly_y_0_0, tly_y_x_0, tly_y_y_0, \
                                     tly_y_z_0, tly_z_0_0, tly_z_x_0, tly_z_y_0, tly_z_z_0, tlz_0_x_0, tlz_0_y_0, \
                                     tlz_0_z_0, tlz_x_0_0, tlz_x_x_0, tlz_x_y_0, tlz_x_z_0, tlz_xx_x_0, tlz_xx_y_0, \
                                     tlz_xx_z_0, tlz_xy_x_0, tlz_xy_y_0, tlz_xy_z_0, tlz_xz_x_0, tlz_xz_y_0, tlz_xz_z_0, \
                                     tlz_y_0_0, tlz_y_x_0, tlz_y_y_0, tlz_y_z_0, tlz_z_0_0, tlz_z_x_0, tlz_z_y_0, \
                                     tlz_z_z_0, tpy_x_x_0, tpy_x_y_0, tpy_x_z_0, tpy_y_x_0, tpy_y_y_0, tpy_y_z_0, \
                                     tpy_z_x_0, tpy_z_y_0, tpy_z_z_0, tpz_x_x_0, tpz_x_y_0, tpz_x_z_0, tpz_y_x_0, \
                                     tpz_y_y_0, tpz_y_z_0, tpz_z_x_0, tpz_z_y_0, tpz_z_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xx_x_0[j] = pa_x[j] * tlx_x_x_0[j] + 0.5 * fl1_fx * tlx_0_x_0[j] + 0.5 * fl1_fx * tlx_x_0_0[j];

            tly_xx_x_0[j] = pa_x[j] * tly_x_x_0[j] + 0.5 * fl1_fx * tly_0_x_0[j] + 0.5 * fl1_fx * tly_x_0_0[j] + 0.5 * fl1_fx * tpz_x_x_0[j] +
                            fl1_fx * fl1_fgb * tdz_x_x_0[j];

            tlz_xx_x_0[j] = pa_x[j] * tlz_x_x_0[j] + 0.5 * fl1_fx * tlz_0_x_0[j] + 0.5 * fl1_fx * tlz_x_0_0[j] - 0.5 * fl1_fx * tpy_x_x_0[j] -
                            fl1_fx * fl1_fgb * tdy_x_x_0[j];

            tlx_xx_y_0[j] = pa_x[j] * tlx_x_y_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j];

            tly_xx_y_0[j] = pa_x[j] * tly_x_y_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] + 0.5 * fl1_fx * tpz_x_y_0[j] + fl1_fx * fl1_fgb * tdz_x_y_0[j];

            tlz_xx_y_0[j] = pa_x[j] * tlz_x_y_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j] - 0.5 * fl1_fx * tpy_x_y_0[j] - fl1_fx * fl1_fgb * tdy_x_y_0[j];

            tlx_xx_z_0[j] = pa_x[j] * tlx_x_z_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j];

            tly_xx_z_0[j] = pa_x[j] * tly_x_z_0[j] + 0.5 * fl1_fx * tly_0_z_0[j] + 0.5 * fl1_fx * tpz_x_z_0[j] + fl1_fx * fl1_fgb * tdz_x_z_0[j];

            tlz_xx_z_0[j] = pa_x[j] * tlz_x_z_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] - 0.5 * fl1_fx * tpy_x_z_0[j] - fl1_fx * fl1_fgb * tdy_x_z_0[j];

            tlx_xy_x_0[j] = pa_x[j] * tlx_y_x_0[j] + 0.5 * fl1_fx * tlx_y_0_0[j];

            tly_xy_x_0[j] = pa_x[j] * tly_y_x_0[j] + 0.5 * fl1_fx * tly_y_0_0[j] + 0.5 * fl1_fx * tpz_y_x_0[j] + fl1_fx * fl1_fgb * tdz_y_x_0[j];

            tlz_xy_x_0[j] = pa_x[j] * tlz_y_x_0[j] + 0.5 * fl1_fx * tlz_y_0_0[j] - 0.5 * fl1_fx * tpy_y_x_0[j] - fl1_fx * fl1_fgb * tdy_y_x_0[j];

            tlx_xy_y_0[j] = pa_x[j] * tlx_y_y_0[j];

            tly_xy_y_0[j] = pa_x[j] * tly_y_y_0[j] + 0.5 * fl1_fx * tpz_y_y_0[j] + fl1_fx * fl1_fgb * tdz_y_y_0[j];

            tlz_xy_y_0[j] = pa_x[j] * tlz_y_y_0[j] - 0.5 * fl1_fx * tpy_y_y_0[j] - fl1_fx * fl1_fgb * tdy_y_y_0[j];

            tlx_xy_z_0[j] = pa_x[j] * tlx_y_z_0[j];

            tly_xy_z_0[j] = pa_x[j] * tly_y_z_0[j] + 0.5 * fl1_fx * tpz_y_z_0[j] + fl1_fx * fl1_fgb * tdz_y_z_0[j];

            tlz_xy_z_0[j] = pa_x[j] * tlz_y_z_0[j] - 0.5 * fl1_fx * tpy_y_z_0[j] - fl1_fx * fl1_fgb * tdy_y_z_0[j];

            tlx_xz_x_0[j] = pa_x[j] * tlx_z_x_0[j] + 0.5 * fl1_fx * tlx_z_0_0[j];

            tly_xz_x_0[j] = pa_x[j] * tly_z_x_0[j] + 0.5 * fl1_fx * tly_z_0_0[j] + 0.5 * fl1_fx * tpz_z_x_0[j] + fl1_fx * fl1_fgb * tdz_z_x_0[j];

            tlz_xz_x_0[j] = pa_x[j] * tlz_z_x_0[j] + 0.5 * fl1_fx * tlz_z_0_0[j] - 0.5 * fl1_fx * tpy_z_x_0[j] - fl1_fx * fl1_fgb * tdy_z_x_0[j];

            tlx_xz_y_0[j] = pa_x[j] * tlx_z_y_0[j];

            tly_xz_y_0[j] = pa_x[j] * tly_z_y_0[j] + 0.5 * fl1_fx * tpz_z_y_0[j] + fl1_fx * fl1_fgb * tdz_z_y_0[j];

            tlz_xz_y_0[j] = pa_x[j] * tlz_z_y_0[j] - 0.5 * fl1_fx * tpy_z_y_0[j] - fl1_fx * fl1_fgb * tdy_z_y_0[j];

            tlx_xz_z_0[j] = pa_x[j] * tlx_z_z_0[j];

            tly_xz_z_0[j] = pa_x[j] * tly_z_z_0[j] + 0.5 * fl1_fx * tpz_z_z_0[j] + fl1_fx * fl1_fgb * tdz_z_z_0[j];

            tlz_xz_z_0[j] = pa_x[j] * tlz_z_z_0[j] - 0.5 * fl1_fx * tpy_z_z_0[j] - fl1_fx * fl1_fgb * tdy_z_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForDP_27_54(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (27,54)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_2_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3);

        auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4);

        auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5);

        auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6);

        auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7);

        auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8);

        auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tlx_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 1);

        auto tly_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 2);

        auto tly_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tdx_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 3);

        auto tdz_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tdx_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 4);

        auto tdz_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tdx_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 5);

        auto tdz_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6);

        auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7);

        auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8);

        auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8);

        // set up pointers to integrals

        auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9);

        auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10);

        auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11);

        auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12);

        auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13);

        auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14);

        auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15);

        auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16);

        auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17);

        auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17);

        // Batch of Integrals (27,54)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_y_x_0, tdx_y_y_0, tdx_y_z_0, tdx_z_x_0, tdx_z_y_0, \
                                     tdx_z_z_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, \
                                     tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, tlx_0_x_0, tlx_0_y_0, tlx_0_z_0, tlx_y_0_0, \
                                     tlx_y_x_0, tlx_y_y_0, tlx_y_z_0, tlx_yy_x_0, tlx_yy_y_0, tlx_yy_z_0, tlx_yz_x_0, \
                                     tlx_yz_y_0, tlx_yz_z_0, tlx_z_0_0, tlx_z_x_0, tlx_z_y_0, tlx_z_z_0, tlx_zz_x_0, \
                                     tlx_zz_y_0, tlx_zz_z_0, tly_0_x_0, tly_0_y_0, tly_0_z_0, tly_y_0_0, tly_y_x_0, \
                                     tly_y_y_0, tly_y_z_0, tly_yy_x_0, tly_yy_y_0, tly_yy_z_0, tly_yz_x_0, tly_yz_y_0, \
                                     tly_yz_z_0, tly_z_0_0, tly_z_x_0, tly_z_y_0, tly_z_z_0, tly_zz_x_0, tly_zz_y_0, \
                                     tly_zz_z_0, tlz_0_x_0, tlz_0_y_0, tlz_0_z_0, tlz_y_0_0, tlz_y_x_0, tlz_y_y_0, \
                                     tlz_y_z_0, tlz_yy_x_0, tlz_yy_y_0, tlz_yy_z_0, tlz_yz_x_0, tlz_yz_y_0, tlz_yz_z_0, \
                                     tlz_z_0_0, tlz_z_x_0, tlz_z_y_0, tlz_z_z_0, tlz_zz_x_0, tlz_zz_y_0, tlz_zz_z_0, \
                                     tpx_y_x_0, tpx_y_y_0, tpx_y_z_0, tpx_z_x_0, tpx_z_y_0, tpx_z_z_0, tpy_z_x_0, \
                                     tpy_z_y_0, tpy_z_z_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, tpz_z_x_0, tpz_z_y_0, \
                                     tpz_z_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yy_x_0[j] = pa_y[j] * tlx_y_x_0[j] + 0.5 * fl1_fx * tlx_0_x_0[j] - 0.5 * fl1_fx * tpz_y_x_0[j] - fl1_fx * fl1_fgb * tdz_y_x_0[j];

            tly_yy_x_0[j] = pa_y[j] * tly_y_x_0[j] + 0.5 * fl1_fx * tly_0_x_0[j];

            tlz_yy_x_0[j] = pa_y[j] * tlz_y_x_0[j] + 0.5 * fl1_fx * tlz_0_x_0[j] + 0.5 * fl1_fx * tpx_y_x_0[j] + fl1_fx * fl1_fgb * tdx_y_x_0[j];

            tlx_yy_y_0[j] = pa_y[j] * tlx_y_y_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j] + 0.5 * fl1_fx * tlx_y_0_0[j] - 0.5 * fl1_fx * tpz_y_y_0[j] -
                            fl1_fx * fl1_fgb * tdz_y_y_0[j];

            tly_yy_y_0[j] = pa_y[j] * tly_y_y_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] + 0.5 * fl1_fx * tly_y_0_0[j];

            tlz_yy_y_0[j] = pa_y[j] * tlz_y_y_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j] + 0.5 * fl1_fx * tlz_y_0_0[j] + 0.5 * fl1_fx * tpx_y_y_0[j] +
                            fl1_fx * fl1_fgb * tdx_y_y_0[j];

            tlx_yy_z_0[j] = pa_y[j] * tlx_y_z_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j] - 0.5 * fl1_fx * tpz_y_z_0[j] - fl1_fx * fl1_fgb * tdz_y_z_0[j];

            tly_yy_z_0[j] = pa_y[j] * tly_y_z_0[j] + 0.5 * fl1_fx * tly_0_z_0[j];

            tlz_yy_z_0[j] = pa_y[j] * tlz_y_z_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] + 0.5 * fl1_fx * tpx_y_z_0[j] + fl1_fx * fl1_fgb * tdx_y_z_0[j];

            tlx_yz_x_0[j] = pa_y[j] * tlx_z_x_0[j] - 0.5 * fl1_fx * tpz_z_x_0[j] - fl1_fx * fl1_fgb * tdz_z_x_0[j];

            tly_yz_x_0[j] = pa_y[j] * tly_z_x_0[j];

            tlz_yz_x_0[j] = pa_y[j] * tlz_z_x_0[j] + 0.5 * fl1_fx * tpx_z_x_0[j] + fl1_fx * fl1_fgb * tdx_z_x_0[j];

            tlx_yz_y_0[j] = pa_y[j] * tlx_z_y_0[j] + 0.5 * fl1_fx * tlx_z_0_0[j] - 0.5 * fl1_fx * tpz_z_y_0[j] - fl1_fx * fl1_fgb * tdz_z_y_0[j];

            tly_yz_y_0[j] = pa_y[j] * tly_z_y_0[j] + 0.5 * fl1_fx * tly_z_0_0[j];

            tlz_yz_y_0[j] = pa_y[j] * tlz_z_y_0[j] + 0.5 * fl1_fx * tlz_z_0_0[j] + 0.5 * fl1_fx * tpx_z_y_0[j] + fl1_fx * fl1_fgb * tdx_z_y_0[j];

            tlx_yz_z_0[j] = pa_y[j] * tlx_z_z_0[j] - 0.5 * fl1_fx * tpz_z_z_0[j] - fl1_fx * fl1_fgb * tdz_z_z_0[j];

            tly_yz_z_0[j] = pa_y[j] * tly_z_z_0[j];

            tlz_yz_z_0[j] = pa_y[j] * tlz_z_z_0[j] + 0.5 * fl1_fx * tpx_z_z_0[j] + fl1_fx * fl1_fgb * tdx_z_z_0[j];

            tlx_zz_x_0[j] = pa_z[j] * tlx_z_x_0[j] + 0.5 * fl1_fx * tlx_0_x_0[j] + 0.5 * fl1_fx * tpy_z_x_0[j] + fl1_fx * fl1_fgb * tdy_z_x_0[j];

            tly_zz_x_0[j] = pa_z[j] * tly_z_x_0[j] + 0.5 * fl1_fx * tly_0_x_0[j] - 0.5 * fl1_fx * tpx_z_x_0[j] - fl1_fx * fl1_fgb * tdx_z_x_0[j];

            tlz_zz_x_0[j] = pa_z[j] * tlz_z_x_0[j] + 0.5 * fl1_fx * tlz_0_x_0[j];

            tlx_zz_y_0[j] = pa_z[j] * tlx_z_y_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j] + 0.5 * fl1_fx * tpy_z_y_0[j] + fl1_fx * fl1_fgb * tdy_z_y_0[j];

            tly_zz_y_0[j] = pa_z[j] * tly_z_y_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] - 0.5 * fl1_fx * tpx_z_y_0[j] - fl1_fx * fl1_fgb * tdx_z_y_0[j];

            tlz_zz_y_0[j] = pa_z[j] * tlz_z_y_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j];

            tlx_zz_z_0[j] = pa_z[j] * tlx_z_z_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j] + 0.5 * fl1_fx * tlx_z_0_0[j] + 0.5 * fl1_fx * tpy_z_z_0[j] +
                            fl1_fx * fl1_fgb * tdy_z_z_0[j];

            tly_zz_z_0[j] = pa_z[j] * tly_z_z_0[j] + 0.5 * fl1_fx * tly_0_z_0[j] + 0.5 * fl1_fx * tly_z_0_0[j] - 0.5 * fl1_fx * tpx_z_z_0[j] -
                            fl1_fx * fl1_fgb * tdx_z_z_0[j];

            tlz_zz_z_0[j] = pa_z[j] * tlz_z_z_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] + 0.5 * fl1_fx * tlz_z_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPF(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForPF_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForPF_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForPF_0_45(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpy_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tpy_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tpy_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tpy_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tpy_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tdx_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx);

        auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx);

        auto tdz_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx);

        auto tdx_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 1);

        auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tdz_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tdx_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 2);

        auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tdz_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tdx_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 3);

        auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tdz_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tdx_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 4);

        auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tdz_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9);

        // set up pointers to integrals

        auto tlx_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx);

        auto tly_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx);

        auto tlz_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx);

        auto tlx_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 1);

        auto tly_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tlz_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tlx_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 2);

        auto tly_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tlz_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tlx_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 3);

        auto tly_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tlz_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tlx_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 4);

        auto tly_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tlz_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tlx_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 5);

        auto tly_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tlz_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tlx_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 6);

        auto tly_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tlz_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tlx_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 7);

        auto tly_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tlz_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tlx_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 8);

        auto tly_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tlz_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tlx_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 9);

        auto tly_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tlz_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 9);

        auto tlx_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 10);

        auto tly_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tlz_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tlx_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 11);

        auto tly_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tlz_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tlx_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 12);

        auto tly_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tlz_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tlx_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 13);

        auto tly_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tlz_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tlx_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 14);

        auto tly_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tlz_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, tdx_0_xyy_0, \
                                     tdx_0_xyz_0, tdy_0_xxx_0, tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xyy_0, tdy_0_xyz_0, \
                                     tdy_0_xzz_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yzz_0, tdy_0_zzz_0, tdz_0_xxx_0, \
                                     tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xzz_0, tdz_0_yyy_0, \
                                     tdz_0_yyz_0, tdz_0_yzz_0, tdz_0_zzz_0, tlx_0_xx_0, tlx_0_xxx_0, tlx_0_xxy_0, \
                                     tlx_0_xxz_0, tlx_0_xy_0, tlx_0_xyy_0, tlx_0_xyz_0, tlx_0_xz_0, tlx_0_xzz_0, \
                                     tlx_0_yy_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yz_0, tlx_0_yzz_0, tlx_0_zz_0, \
                                     tlx_0_zzz_0, tlx_x_xxx_0, tlx_x_xxy_0, tlx_x_xxz_0, tlx_x_xyy_0, tlx_x_xyz_0, \
                                     tlx_x_xzz_0, tlx_x_yyy_0, tlx_x_yyz_0, tlx_x_yzz_0, tlx_x_zzz_0, tlx_y_xxx_0, \
                                     tlx_y_xxy_0, tlx_y_xxz_0, tlx_y_xyy_0, tlx_y_xyz_0, tly_0_xx_0, tly_0_xxx_0, \
                                     tly_0_xxy_0, tly_0_xxz_0, tly_0_xy_0, tly_0_xyy_0, tly_0_xyz_0, tly_0_xz_0, \
                                     tly_0_xzz_0, tly_0_yy_0, tly_0_yyy_0, tly_0_yyz_0, tly_0_yz_0, tly_0_yzz_0, \
                                     tly_0_zz_0, tly_0_zzz_0, tly_x_xxx_0, tly_x_xxy_0, tly_x_xxz_0, tly_x_xyy_0, \
                                     tly_x_xyz_0, tly_x_xzz_0, tly_x_yyy_0, tly_x_yyz_0, tly_x_yzz_0, tly_x_zzz_0, \
                                     tly_y_xxx_0, tly_y_xxy_0, tly_y_xxz_0, tly_y_xyy_0, tly_y_xyz_0, tlz_0_xx_0, \
                                     tlz_0_xxx_0, tlz_0_xxy_0, tlz_0_xxz_0, tlz_0_xy_0, tlz_0_xyy_0, tlz_0_xyz_0, \
                                     tlz_0_xz_0, tlz_0_xzz_0, tlz_0_yy_0, tlz_0_yyy_0, tlz_0_yyz_0, tlz_0_yz_0, \
                                     tlz_0_yzz_0, tlz_0_zz_0, tlz_0_zzz_0, tlz_x_xxx_0, tlz_x_xxy_0, tlz_x_xxz_0, \
                                     tlz_x_xyy_0, tlz_x_xyz_0, tlz_x_xzz_0, tlz_x_yyy_0, tlz_x_yyz_0, tlz_x_yzz_0, \
                                     tlz_x_zzz_0, tlz_y_xxx_0, tlz_y_xxy_0, tlz_y_xxz_0, tlz_y_xyy_0, tlz_y_xyz_0, \
                                     tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, tpx_0_xyy_0, tpx_0_xyz_0, tpy_0_xxx_0, \
                                     tpy_0_xxy_0, tpy_0_xxz_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xzz_0, tpy_0_yyy_0, \
                                     tpy_0_yyz_0, tpy_0_yzz_0, tpy_0_zzz_0, tpz_0_xxx_0, tpz_0_xxy_0, tpz_0_xxz_0, \
                                     tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xzz_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yzz_0, \
                                     tpz_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_x_xxx_0[j] = pa_x[j] * tlx_0_xxx_0[j] + 1.5 * fl1_fx * tlx_0_xx_0[j];

            tly_x_xxx_0[j] =
                pa_x[j] * tly_0_xxx_0[j] + 1.5 * fl1_fx * tly_0_xx_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j] + fl1_fx * fl1_fgb * tdz_0_xxx_0[j];

            tlz_x_xxx_0[j] =
                pa_x[j] * tlz_0_xxx_0[j] + 1.5 * fl1_fx * tlz_0_xx_0[j] - 0.5 * fl1_fx * tpy_0_xxx_0[j] - fl1_fx * fl1_fgb * tdy_0_xxx_0[j];

            tlx_x_xxy_0[j] = pa_x[j] * tlx_0_xxy_0[j] + fl1_fx * tlx_0_xy_0[j];

            tly_x_xxy_0[j] = pa_x[j] * tly_0_xxy_0[j] + fl1_fx * tly_0_xy_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] + fl1_fx * fl1_fgb * tdz_0_xxy_0[j];

            tlz_x_xxy_0[j] = pa_x[j] * tlz_0_xxy_0[j] + fl1_fx * tlz_0_xy_0[j] - 0.5 * fl1_fx * tpy_0_xxy_0[j] - fl1_fx * fl1_fgb * tdy_0_xxy_0[j];

            tlx_x_xxz_0[j] = pa_x[j] * tlx_0_xxz_0[j] + fl1_fx * tlx_0_xz_0[j];

            tly_x_xxz_0[j] = pa_x[j] * tly_0_xxz_0[j] + fl1_fx * tly_0_xz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j] + fl1_fx * fl1_fgb * tdz_0_xxz_0[j];

            tlz_x_xxz_0[j] = pa_x[j] * tlz_0_xxz_0[j] + fl1_fx * tlz_0_xz_0[j] - 0.5 * fl1_fx * tpy_0_xxz_0[j] - fl1_fx * fl1_fgb * tdy_0_xxz_0[j];

            tlx_x_xyy_0[j] = pa_x[j] * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j];

            tly_x_xyy_0[j] =
                pa_x[j] * tly_0_xyy_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] + fl1_fx * fl1_fgb * tdz_0_xyy_0[j];

            tlz_x_xyy_0[j] =
                pa_x[j] * tlz_0_xyy_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j] - 0.5 * fl1_fx * tpy_0_xyy_0[j] - fl1_fx * fl1_fgb * tdy_0_xyy_0[j];

            tlx_x_xyz_0[j] = pa_x[j] * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_0_yz_0[j];

            tly_x_xyz_0[j] =
                pa_x[j] * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_xyz_0[j] + fl1_fx * fl1_fgb * tdz_0_xyz_0[j];

            tlz_x_xyz_0[j] =
                pa_x[j] * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_0_yz_0[j] - 0.5 * fl1_fx * tpy_0_xyz_0[j] - fl1_fx * fl1_fgb * tdy_0_xyz_0[j];

            tlx_x_xzz_0[j] = pa_x[j] * tlx_0_xzz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j];

            tly_x_xzz_0[j] =
                pa_x[j] * tly_0_xzz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j] + fl1_fx * fl1_fgb * tdz_0_xzz_0[j];

            tlz_x_xzz_0[j] =
                pa_x[j] * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] - 0.5 * fl1_fx * tpy_0_xzz_0[j] - fl1_fx * fl1_fgb * tdy_0_xzz_0[j];

            tlx_x_yyy_0[j] = pa_x[j] * tlx_0_yyy_0[j];

            tly_x_yyy_0[j] = pa_x[j] * tly_0_yyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j] + fl1_fx * fl1_fgb * tdz_0_yyy_0[j];

            tlz_x_yyy_0[j] = pa_x[j] * tlz_0_yyy_0[j] - 0.5 * fl1_fx * tpy_0_yyy_0[j] - fl1_fx * fl1_fgb * tdy_0_yyy_0[j];

            tlx_x_yyz_0[j] = pa_x[j] * tlx_0_yyz_0[j];

            tly_x_yyz_0[j] = pa_x[j] * tly_0_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j] + fl1_fx * fl1_fgb * tdz_0_yyz_0[j];

            tlz_x_yyz_0[j] = pa_x[j] * tlz_0_yyz_0[j] - 0.5 * fl1_fx * tpy_0_yyz_0[j] - fl1_fx * fl1_fgb * tdy_0_yyz_0[j];

            tlx_x_yzz_0[j] = pa_x[j] * tlx_0_yzz_0[j];

            tly_x_yzz_0[j] = pa_x[j] * tly_0_yzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j] + fl1_fx * fl1_fgb * tdz_0_yzz_0[j];

            tlz_x_yzz_0[j] = pa_x[j] * tlz_0_yzz_0[j] - 0.5 * fl1_fx * tpy_0_yzz_0[j] - fl1_fx * fl1_fgb * tdy_0_yzz_0[j];

            tlx_x_zzz_0[j] = pa_x[j] * tlx_0_zzz_0[j];

            tly_x_zzz_0[j] = pa_x[j] * tly_0_zzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j] + fl1_fx * fl1_fgb * tdz_0_zzz_0[j];

            tlz_x_zzz_0[j] = pa_x[j] * tlz_0_zzz_0[j] - 0.5 * fl1_fx * tpy_0_zzz_0[j] - fl1_fx * fl1_fgb * tdy_0_zzz_0[j];

            tlx_y_xxx_0[j] = pa_y[j] * tlx_0_xxx_0[j] - 0.5 * fl1_fx * tpz_0_xxx_0[j] - fl1_fx * fl1_fgb * tdz_0_xxx_0[j];

            tly_y_xxx_0[j] = pa_y[j] * tly_0_xxx_0[j];

            tlz_y_xxx_0[j] = pa_y[j] * tlz_0_xxx_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j] + fl1_fx * fl1_fgb * tdx_0_xxx_0[j];

            tlx_y_xxy_0[j] =
                pa_y[j] * tlx_0_xxy_0[j] + 0.5 * fl1_fx * tlx_0_xx_0[j] - 0.5 * fl1_fx * tpz_0_xxy_0[j] - fl1_fx * fl1_fgb * tdz_0_xxy_0[j];

            tly_y_xxy_0[j] = pa_y[j] * tly_0_xxy_0[j] + 0.5 * fl1_fx * tly_0_xx_0[j];

            tlz_y_xxy_0[j] =
                pa_y[j] * tlz_0_xxy_0[j] + 0.5 * fl1_fx * tlz_0_xx_0[j] + 0.5 * fl1_fx * tpx_0_xxy_0[j] + fl1_fx * fl1_fgb * tdx_0_xxy_0[j];

            tlx_y_xxz_0[j] = pa_y[j] * tlx_0_xxz_0[j] - 0.5 * fl1_fx * tpz_0_xxz_0[j] - fl1_fx * fl1_fgb * tdz_0_xxz_0[j];

            tly_y_xxz_0[j] = pa_y[j] * tly_0_xxz_0[j];

            tlz_y_xxz_0[j] = pa_y[j] * tlz_0_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xxz_0[j] + fl1_fx * fl1_fgb * tdx_0_xxz_0[j];

            tlx_y_xyy_0[j] = pa_y[j] * tlx_0_xyy_0[j] + fl1_fx * tlx_0_xy_0[j] - 0.5 * fl1_fx * tpz_0_xyy_0[j] - fl1_fx * fl1_fgb * tdz_0_xyy_0[j];

            tly_y_xyy_0[j] = pa_y[j] * tly_0_xyy_0[j] + fl1_fx * tly_0_xy_0[j];

            tlz_y_xyy_0[j] = pa_y[j] * tlz_0_xyy_0[j] + fl1_fx * tlz_0_xy_0[j] + 0.5 * fl1_fx * tpx_0_xyy_0[j] + fl1_fx * fl1_fgb * tdx_0_xyy_0[j];

            tlx_y_xyz_0[j] =
                pa_y[j] * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_0_xz_0[j] - 0.5 * fl1_fx * tpz_0_xyz_0[j] - fl1_fx * fl1_fgb * tdz_0_xyz_0[j];

            tly_y_xyz_0[j] = pa_y[j] * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_0_xz_0[j];

            tlz_y_xyz_0[j] =
                pa_y[j] * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_0_xz_0[j] + 0.5 * fl1_fx * tpx_0_xyz_0[j] + fl1_fx * fl1_fgb * tdx_0_xyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPF_45_90(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

        auto tpy_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tpx_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 6);

        auto tpy_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tpx_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 7);

        auto tpy_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tpx_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 8);

        auto tpy_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tpx_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 9);

        auto tpy_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tdx_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx);

        auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx);

        auto tdx_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 1);

        auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tdx_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 2);

        auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tdx_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 3);

        auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tdx_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 4);

        auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tdx_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 5);

        auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tdx_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 6);

        auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tdx_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 7);

        auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tdx_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 8);

        auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tdx_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 9);

        auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9);

        // set up pointers to integrals

        auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15);

        auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16);

        auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tlz_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tlx_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 17);

        auto tly_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tlz_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tlx_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 18);

        auto tly_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tlz_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tlx_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 19);

        auto tly_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tlz_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20);

        auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21);

        auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22);

        auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23);

        auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24);

        auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25);

        auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26);

        auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27);

        auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28);

        auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29);

        auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, tdx_0_xyy_0, \
                                     tdx_0_xyz_0, tdx_0_xzz_0, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yzz_0, tdx_0_zzz_0, \
                                     tdy_0_xxx_0, tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xzz_0, \
                                     tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yzz_0, tdy_0_zzz_0, tdz_0_xzz_0, tdz_0_yyy_0, \
                                     tdz_0_yyz_0, tdz_0_yzz_0, tdz_0_zzz_0, tlx_0_xx_0, tlx_0_xxx_0, tlx_0_xxy_0, \
                                     tlx_0_xxz_0, tlx_0_xy_0, tlx_0_xyy_0, tlx_0_xyz_0, tlx_0_xz_0, tlx_0_xzz_0, \
                                     tlx_0_yy_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yz_0, tlx_0_yzz_0, tlx_0_zz_0, \
                                     tlx_0_zzz_0, tlx_y_xzz_0, tlx_y_yyy_0, tlx_y_yyz_0, tlx_y_yzz_0, tlx_y_zzz_0, \
                                     tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, tlx_z_xyy_0, tlx_z_xyz_0, tlx_z_xzz_0, \
                                     tlx_z_yyy_0, tlx_z_yyz_0, tlx_z_yzz_0, tlx_z_zzz_0, tly_0_xx_0, tly_0_xxx_0, \
                                     tly_0_xxy_0, tly_0_xxz_0, tly_0_xy_0, tly_0_xyy_0, tly_0_xyz_0, tly_0_xz_0, \
                                     tly_0_xzz_0, tly_0_yy_0, tly_0_yyy_0, tly_0_yyz_0, tly_0_yz_0, tly_0_yzz_0, \
                                     tly_0_zz_0, tly_0_zzz_0, tly_y_xzz_0, tly_y_yyy_0, tly_y_yyz_0, tly_y_yzz_0, \
                                     tly_y_zzz_0, tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xyy_0, tly_z_xyz_0, \
                                     tly_z_xzz_0, tly_z_yyy_0, tly_z_yyz_0, tly_z_yzz_0, tly_z_zzz_0, tlz_0_xx_0, \
                                     tlz_0_xxx_0, tlz_0_xxy_0, tlz_0_xxz_0, tlz_0_xy_0, tlz_0_xyy_0, tlz_0_xyz_0, \
                                     tlz_0_xz_0, tlz_0_xzz_0, tlz_0_yy_0, tlz_0_yyy_0, tlz_0_yyz_0, tlz_0_yz_0, \
                                     tlz_0_yzz_0, tlz_0_zz_0, tlz_0_zzz_0, tlz_y_xzz_0, tlz_y_yyy_0, tlz_y_yyz_0, \
                                     tlz_y_yzz_0, tlz_y_zzz_0, tlz_z_xxx_0, tlz_z_xxy_0, tlz_z_xxz_0, tlz_z_xyy_0, \
                                     tlz_z_xyz_0, tlz_z_xzz_0, tlz_z_yyy_0, tlz_z_yyz_0, tlz_z_yzz_0, tlz_z_zzz_0, \
                                     tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, tpx_0_xyy_0, tpx_0_xyz_0, tpx_0_xzz_0, \
                                     tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yzz_0, tpx_0_zzz_0, tpy_0_xxx_0, tpy_0_xxy_0, \
                                     tpy_0_xxz_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xzz_0, tpy_0_yyy_0, tpy_0_yyz_0, \
                                     tpy_0_yzz_0, tpy_0_zzz_0, tpz_0_xzz_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yzz_0, \
                                     tpz_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_y_xzz_0[j] = pa_y[j] * tlx_0_xzz_0[j] - 0.5 * fl1_fx * tpz_0_xzz_0[j] - fl1_fx * fl1_fgb * tdz_0_xzz_0[j];

            tly_y_xzz_0[j] = pa_y[j] * tly_0_xzz_0[j];

            tlz_y_xzz_0[j] = pa_y[j] * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tpx_0_xzz_0[j] + fl1_fx * fl1_fgb * tdx_0_xzz_0[j];

            tlx_y_yyy_0[j] =
                pa_y[j] * tlx_0_yyy_0[j] + 1.5 * fl1_fx * tlx_0_yy_0[j] - 0.5 * fl1_fx * tpz_0_yyy_0[j] - fl1_fx * fl1_fgb * tdz_0_yyy_0[j];

            tly_y_yyy_0[j] = pa_y[j] * tly_0_yyy_0[j] + 1.5 * fl1_fx * tly_0_yy_0[j];

            tlz_y_yyy_0[j] =
                pa_y[j] * tlz_0_yyy_0[j] + 1.5 * fl1_fx * tlz_0_yy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j] + fl1_fx * fl1_fgb * tdx_0_yyy_0[j];

            tlx_y_yyz_0[j] = pa_y[j] * tlx_0_yyz_0[j] + fl1_fx * tlx_0_yz_0[j] - 0.5 * fl1_fx * tpz_0_yyz_0[j] - fl1_fx * fl1_fgb * tdz_0_yyz_0[j];

            tly_y_yyz_0[j] = pa_y[j] * tly_0_yyz_0[j] + fl1_fx * tly_0_yz_0[j];

            tlz_y_yyz_0[j] = pa_y[j] * tlz_0_yyz_0[j] + fl1_fx * tlz_0_yz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] + fl1_fx * fl1_fgb * tdx_0_yyz_0[j];

            tlx_y_yzz_0[j] =
                pa_y[j] * tlx_0_yzz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j] - 0.5 * fl1_fx * tpz_0_yzz_0[j] - fl1_fx * fl1_fgb * tdz_0_yzz_0[j];

            tly_y_yzz_0[j] = pa_y[j] * tly_0_yzz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j];

            tlz_y_yzz_0[j] =
                pa_y[j] * tlz_0_yzz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] + fl1_fx * fl1_fgb * tdx_0_yzz_0[j];

            tlx_y_zzz_0[j] = pa_y[j] * tlx_0_zzz_0[j] - 0.5 * fl1_fx * tpz_0_zzz_0[j] - fl1_fx * fl1_fgb * tdz_0_zzz_0[j];

            tly_y_zzz_0[j] = pa_y[j] * tly_0_zzz_0[j];

            tlz_y_zzz_0[j] = pa_y[j] * tlz_0_zzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j] + fl1_fx * fl1_fgb * tdx_0_zzz_0[j];

            tlx_z_xxx_0[j] = pa_z[j] * tlx_0_xxx_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j] + fl1_fx * fl1_fgb * tdy_0_xxx_0[j];

            tly_z_xxx_0[j] = pa_z[j] * tly_0_xxx_0[j] - 0.5 * fl1_fx * tpx_0_xxx_0[j] - fl1_fx * fl1_fgb * tdx_0_xxx_0[j];

            tlz_z_xxx_0[j] = pa_z[j] * tlz_0_xxx_0[j];

            tlx_z_xxy_0[j] = pa_z[j] * tlx_0_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xxy_0[j] + fl1_fx * fl1_fgb * tdy_0_xxy_0[j];

            tly_z_xxy_0[j] = pa_z[j] * tly_0_xxy_0[j] - 0.5 * fl1_fx * tpx_0_xxy_0[j] - fl1_fx * fl1_fgb * tdx_0_xxy_0[j];

            tlz_z_xxy_0[j] = pa_z[j] * tlz_0_xxy_0[j];

            tlx_z_xxz_0[j] =
                pa_z[j] * tlx_0_xxz_0[j] + 0.5 * fl1_fx * tlx_0_xx_0[j] + 0.5 * fl1_fx * tpy_0_xxz_0[j] + fl1_fx * fl1_fgb * tdy_0_xxz_0[j];

            tly_z_xxz_0[j] =
                pa_z[j] * tly_0_xxz_0[j] + 0.5 * fl1_fx * tly_0_xx_0[j] - 0.5 * fl1_fx * tpx_0_xxz_0[j] - fl1_fx * fl1_fgb * tdx_0_xxz_0[j];

            tlz_z_xxz_0[j] = pa_z[j] * tlz_0_xxz_0[j] + 0.5 * fl1_fx * tlz_0_xx_0[j];

            tlx_z_xyy_0[j] = pa_z[j] * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tpy_0_xyy_0[j] + fl1_fx * fl1_fgb * tdy_0_xyy_0[j];

            tly_z_xyy_0[j] = pa_z[j] * tly_0_xyy_0[j] - 0.5 * fl1_fx * tpx_0_xyy_0[j] - fl1_fx * fl1_fgb * tdx_0_xyy_0[j];

            tlz_z_xyy_0[j] = pa_z[j] * tlz_0_xyy_0[j];

            tlx_z_xyz_0[j] =
                pa_z[j] * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_0_xy_0[j] + 0.5 * fl1_fx * tpy_0_xyz_0[j] + fl1_fx * fl1_fgb * tdy_0_xyz_0[j];

            tly_z_xyz_0[j] =
                pa_z[j] * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_0_xy_0[j] - 0.5 * fl1_fx * tpx_0_xyz_0[j] - fl1_fx * fl1_fgb * tdx_0_xyz_0[j];

            tlz_z_xyz_0[j] = pa_z[j] * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_0_xy_0[j];

            tlx_z_xzz_0[j] = pa_z[j] * tlx_0_xzz_0[j] + fl1_fx * tlx_0_xz_0[j] + 0.5 * fl1_fx * tpy_0_xzz_0[j] + fl1_fx * fl1_fgb * tdy_0_xzz_0[j];

            tly_z_xzz_0[j] = pa_z[j] * tly_0_xzz_0[j] + fl1_fx * tly_0_xz_0[j] - 0.5 * fl1_fx * tpx_0_xzz_0[j] - fl1_fx * fl1_fgb * tdx_0_xzz_0[j];

            tlz_z_xzz_0[j] = pa_z[j] * tlz_0_xzz_0[j] + fl1_fx * tlz_0_xz_0[j];

            tlx_z_yyy_0[j] = pa_z[j] * tlx_0_yyy_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j] + fl1_fx * fl1_fgb * tdy_0_yyy_0[j];

            tly_z_yyy_0[j] = pa_z[j] * tly_0_yyy_0[j] - 0.5 * fl1_fx * tpx_0_yyy_0[j] - fl1_fx * fl1_fgb * tdx_0_yyy_0[j];

            tlz_z_yyy_0[j] = pa_z[j] * tlz_0_yyy_0[j];

            tlx_z_yyz_0[j] =
                pa_z[j] * tlx_0_yyz_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j] + 0.5 * fl1_fx * tpy_0_yyz_0[j] + fl1_fx * fl1_fgb * tdy_0_yyz_0[j];

            tly_z_yyz_0[j] =
                pa_z[j] * tly_0_yyz_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] - 0.5 * fl1_fx * tpx_0_yyz_0[j] - fl1_fx * fl1_fgb * tdx_0_yyz_0[j];

            tlz_z_yyz_0[j] = pa_z[j] * tlz_0_yyz_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j];

            tlx_z_yzz_0[j] = pa_z[j] * tlx_0_yzz_0[j] + fl1_fx * tlx_0_yz_0[j] + 0.5 * fl1_fx * tpy_0_yzz_0[j] + fl1_fx * fl1_fgb * tdy_0_yzz_0[j];

            tly_z_yzz_0[j] = pa_z[j] * tly_0_yzz_0[j] + fl1_fx * tly_0_yz_0[j] - 0.5 * fl1_fx * tpx_0_yzz_0[j] - fl1_fx * fl1_fgb * tdx_0_yzz_0[j];

            tlz_z_yzz_0[j] = pa_z[j] * tlz_0_yzz_0[j] + fl1_fx * tlz_0_yz_0[j];

            tlx_z_zzz_0[j] =
                pa_z[j] * tlx_0_zzz_0[j] + 1.5 * fl1_fx * tlx_0_zz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j] + fl1_fx * fl1_fgb * tdy_0_zzz_0[j];

            tly_z_zzz_0[j] =
                pa_z[j] * tly_0_zzz_0[j] + 1.5 * fl1_fx * tly_0_zz_0[j] - 0.5 * fl1_fx * tpx_0_zzz_0[j] - fl1_fx * fl1_fgb * tdx_0_zzz_0[j];

            tlz_z_zzz_0[j] = pa_z[j] * tlz_0_zzz_0[j] + 1.5 * fl1_fx * tlz_0_zz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFP(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForFP_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFP_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForFP_0_45(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tlx_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx);

        auto tly_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx);

        auto tlz_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx);

        auto tlx_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 1);

        auto tly_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tlz_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tlx_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 2);

        auto tly_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tlz_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tlx_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 3);

        auto tly_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tlz_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tlx_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 4);

        auto tly_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tlz_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tlx_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 5);

        auto tly_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tlz_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tlx_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 6);

        auto tly_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tlz_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tlx_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 7);

        auto tly_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tlz_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tlx_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 8);

        auto tly_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tlz_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9);

        auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10);

        auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11);

        auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12);

        auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13);

        auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14);

        auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx);

        auto tly_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx);

        auto tlz_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx);

        auto tlx_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 1);

        auto tly_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tlz_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tlx_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 2);

        auto tly_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tlz_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3);

        auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4);

        auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5);

        auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6);

        auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7);

        auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8);

        auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tlx_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx);

        auto tly_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx);

        auto tlz_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx);

        auto tlx_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 1);

        auto tly_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 2);

        auto tly_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 3);

        auto tly_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 4);

        auto tly_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpy_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx);

        auto tpz_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx);

        auto tpy_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tpy_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tpy_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tpy_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tpy_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tpy_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tpy_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tpy_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tdy_xx_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx);

        auto tdz_xx_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx);

        auto tdy_xx_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tdz_xx_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tdy_xx_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tdz_xx_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tdy_xy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tdz_xy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tdy_xy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tdz_xy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tdy_xy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tdz_xy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tdy_xz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tdz_xz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tdy_xz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tdz_xz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tdy_xz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tdz_xz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tdy_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tdz_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tdy_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tdz_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tdy_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tdz_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tdy_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tdz_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tdy_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tdz_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tdy_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tdz_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 14);

        // set up pointers to integrals

        auto tlx_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx);

        auto tly_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx);

        auto tlz_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx);

        auto tlx_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 1);

        auto tly_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tlz_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tlx_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 2);

        auto tly_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tlz_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tlx_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 3);

        auto tly_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tlz_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tlx_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 4);

        auto tly_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tlz_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tlx_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 5);

        auto tly_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tlz_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tlx_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 6);

        auto tly_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tlz_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tlx_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 7);

        auto tly_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tlz_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tlx_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 8);

        auto tly_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tlz_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tlx_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 9);

        auto tly_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tlz_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tlx_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 10);

        auto tly_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tlz_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tlx_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 11);

        auto tly_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tlz_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tlx_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 12);

        auto tly_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tlz_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tlx_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 13);

        auto tly_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tlz_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tlx_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 14);

        auto tly_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tlz_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xx_x_0, tdy_xx_y_0, tdy_xx_z_0, tdy_xy_x_0, tdy_xy_y_0, \
                                     tdy_xy_z_0, tdy_xz_x_0, tdy_xz_y_0, tdy_xz_z_0, tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, \
                                     tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, tdz_xx_x_0, tdz_xx_y_0, tdz_xx_z_0, tdz_xy_x_0, \
                                     tdz_xy_y_0, tdz_xy_z_0, tdz_xz_x_0, tdz_xz_y_0, tdz_xz_z_0, tdz_yy_x_0, tdz_yy_y_0, \
                                     tdz_yy_z_0, tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, tlx_x_x_0, tlx_x_y_0, tlx_x_z_0, \
                                     tlx_xx_0_0, tlx_xx_x_0, tlx_xx_y_0, tlx_xx_z_0, tlx_xxx_x_0, tlx_xxx_y_0, \
                                     tlx_xxx_z_0, tlx_xxy_x_0, tlx_xxy_y_0, tlx_xxy_z_0, tlx_xxz_x_0, tlx_xxz_y_0, \
                                     tlx_xxz_z_0, tlx_xy_0_0, tlx_xy_x_0, tlx_xy_y_0, tlx_xy_z_0, tlx_xyy_x_0, \
                                     tlx_xyy_y_0, tlx_xyy_z_0, tlx_xyz_x_0, tlx_xyz_y_0, tlx_xyz_z_0, tlx_xz_0_0, \
                                     tlx_xz_x_0, tlx_xz_y_0, tlx_xz_z_0, tlx_y_x_0, tlx_y_y_0, tlx_y_z_0, tlx_yy_0_0, \
                                     tlx_yy_x_0, tlx_yy_y_0, tlx_yy_z_0, tlx_yz_0_0, tlx_yz_x_0, tlx_yz_y_0, tlx_yz_z_0, \
                                     tlx_z_x_0, tlx_z_y_0, tlx_z_z_0, tly_x_x_0, tly_x_y_0, tly_x_z_0, tly_xx_0_0, \
                                     tly_xx_x_0, tly_xx_y_0, tly_xx_z_0, tly_xxx_x_0, tly_xxx_y_0, tly_xxx_z_0, \
                                     tly_xxy_x_0, tly_xxy_y_0, tly_xxy_z_0, tly_xxz_x_0, tly_xxz_y_0, tly_xxz_z_0, \
                                     tly_xy_0_0, tly_xy_x_0, tly_xy_y_0, tly_xy_z_0, tly_xyy_x_0, tly_xyy_y_0, \
                                     tly_xyy_z_0, tly_xyz_x_0, tly_xyz_y_0, tly_xyz_z_0, tly_xz_0_0, tly_xz_x_0, \
                                     tly_xz_y_0, tly_xz_z_0, tly_y_x_0, tly_y_y_0, tly_y_z_0, tly_yy_0_0, tly_yy_x_0, \
                                     tly_yy_y_0, tly_yy_z_0, tly_yz_0_0, tly_yz_x_0, tly_yz_y_0, tly_yz_z_0, tly_z_x_0, \
                                     tly_z_y_0, tly_z_z_0, tlz_x_x_0, tlz_x_y_0, tlz_x_z_0, tlz_xx_0_0, tlz_xx_x_0, \
                                     tlz_xx_y_0, tlz_xx_z_0, tlz_xxx_x_0, tlz_xxx_y_0, tlz_xxx_z_0, tlz_xxy_x_0, \
                                     tlz_xxy_y_0, tlz_xxy_z_0, tlz_xxz_x_0, tlz_xxz_y_0, tlz_xxz_z_0, tlz_xy_0_0, \
                                     tlz_xy_x_0, tlz_xy_y_0, tlz_xy_z_0, tlz_xyy_x_0, tlz_xyy_y_0, tlz_xyy_z_0, \
                                     tlz_xyz_x_0, tlz_xyz_y_0, tlz_xyz_z_0, tlz_xz_0_0, tlz_xz_x_0, tlz_xz_y_0, \
                                     tlz_xz_z_0, tlz_y_x_0, tlz_y_y_0, tlz_y_z_0, tlz_yy_0_0, tlz_yy_x_0, tlz_yy_y_0, \
                                     tlz_yy_z_0, tlz_yz_0_0, tlz_yz_x_0, tlz_yz_y_0, tlz_yz_z_0, tlz_z_x_0, tlz_z_y_0, \
                                     tlz_z_z_0, tpy_xx_x_0, tpy_xx_y_0, tpy_xx_z_0, tpy_xy_x_0, tpy_xy_y_0, tpy_xy_z_0, \
                                     tpy_xz_x_0, tpy_xz_y_0, tpy_xz_z_0, tpy_yy_x_0, tpy_yy_y_0, tpy_yy_z_0, tpy_yz_x_0, \
                                     tpy_yz_y_0, tpy_yz_z_0, tpz_xx_x_0, tpz_xx_y_0, tpz_xx_z_0, tpz_xy_x_0, tpz_xy_y_0, \
                                     tpz_xy_z_0, tpz_xz_x_0, tpz_xz_y_0, tpz_xz_z_0, tpz_yy_x_0, tpz_yy_y_0, tpz_yy_z_0, \
                                     tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxx_x_0[j] = pa_x[j] * tlx_xx_x_0[j] + fl1_fx * tlx_x_x_0[j] + 0.5 * fl1_fx * tlx_xx_0_0[j];

            tly_xxx_x_0[j] = pa_x[j] * tly_xx_x_0[j] + fl1_fx * tly_x_x_0[j] + 0.5 * fl1_fx * tly_xx_0_0[j] + 0.5 * fl1_fx * tpz_xx_x_0[j] +
                             fl1_fx * fl1_fgb * tdz_xx_x_0[j];

            tlz_xxx_x_0[j] = pa_x[j] * tlz_xx_x_0[j] + fl1_fx * tlz_x_x_0[j] + 0.5 * fl1_fx * tlz_xx_0_0[j] - 0.5 * fl1_fx * tpy_xx_x_0[j] -
                             fl1_fx * fl1_fgb * tdy_xx_x_0[j];

            tlx_xxx_y_0[j] = pa_x[j] * tlx_xx_y_0[j] + fl1_fx * tlx_x_y_0[j];

            tly_xxx_y_0[j] = pa_x[j] * tly_xx_y_0[j] + fl1_fx * tly_x_y_0[j] + 0.5 * fl1_fx * tpz_xx_y_0[j] + fl1_fx * fl1_fgb * tdz_xx_y_0[j];

            tlz_xxx_y_0[j] = pa_x[j] * tlz_xx_y_0[j] + fl1_fx * tlz_x_y_0[j] - 0.5 * fl1_fx * tpy_xx_y_0[j] - fl1_fx * fl1_fgb * tdy_xx_y_0[j];

            tlx_xxx_z_0[j] = pa_x[j] * tlx_xx_z_0[j] + fl1_fx * tlx_x_z_0[j];

            tly_xxx_z_0[j] = pa_x[j] * tly_xx_z_0[j] + fl1_fx * tly_x_z_0[j] + 0.5 * fl1_fx * tpz_xx_z_0[j] + fl1_fx * fl1_fgb * tdz_xx_z_0[j];

            tlz_xxx_z_0[j] = pa_x[j] * tlz_xx_z_0[j] + fl1_fx * tlz_x_z_0[j] - 0.5 * fl1_fx * tpy_xx_z_0[j] - fl1_fx * fl1_fgb * tdy_xx_z_0[j];

            tlx_xxy_x_0[j] = pa_x[j] * tlx_xy_x_0[j] + 0.5 * fl1_fx * tlx_y_x_0[j] + 0.5 * fl1_fx * tlx_xy_0_0[j];

            tly_xxy_x_0[j] = pa_x[j] * tly_xy_x_0[j] + 0.5 * fl1_fx * tly_y_x_0[j] + 0.5 * fl1_fx * tly_xy_0_0[j] + 0.5 * fl1_fx * tpz_xy_x_0[j] +
                             fl1_fx * fl1_fgb * tdz_xy_x_0[j];

            tlz_xxy_x_0[j] = pa_x[j] * tlz_xy_x_0[j] + 0.5 * fl1_fx * tlz_y_x_0[j] + 0.5 * fl1_fx * tlz_xy_0_0[j] - 0.5 * fl1_fx * tpy_xy_x_0[j] -
                             fl1_fx * fl1_fgb * tdy_xy_x_0[j];

            tlx_xxy_y_0[j] = pa_x[j] * tlx_xy_y_0[j] + 0.5 * fl1_fx * tlx_y_y_0[j];

            tly_xxy_y_0[j] = pa_x[j] * tly_xy_y_0[j] + 0.5 * fl1_fx * tly_y_y_0[j] + 0.5 * fl1_fx * tpz_xy_y_0[j] + fl1_fx * fl1_fgb * tdz_xy_y_0[j];

            tlz_xxy_y_0[j] = pa_x[j] * tlz_xy_y_0[j] + 0.5 * fl1_fx * tlz_y_y_0[j] - 0.5 * fl1_fx * tpy_xy_y_0[j] - fl1_fx * fl1_fgb * tdy_xy_y_0[j];

            tlx_xxy_z_0[j] = pa_x[j] * tlx_xy_z_0[j] + 0.5 * fl1_fx * tlx_y_z_0[j];

            tly_xxy_z_0[j] = pa_x[j] * tly_xy_z_0[j] + 0.5 * fl1_fx * tly_y_z_0[j] + 0.5 * fl1_fx * tpz_xy_z_0[j] + fl1_fx * fl1_fgb * tdz_xy_z_0[j];

            tlz_xxy_z_0[j] = pa_x[j] * tlz_xy_z_0[j] + 0.5 * fl1_fx * tlz_y_z_0[j] - 0.5 * fl1_fx * tpy_xy_z_0[j] - fl1_fx * fl1_fgb * tdy_xy_z_0[j];

            tlx_xxz_x_0[j] = pa_x[j] * tlx_xz_x_0[j] + 0.5 * fl1_fx * tlx_z_x_0[j] + 0.5 * fl1_fx * tlx_xz_0_0[j];

            tly_xxz_x_0[j] = pa_x[j] * tly_xz_x_0[j] + 0.5 * fl1_fx * tly_z_x_0[j] + 0.5 * fl1_fx * tly_xz_0_0[j] + 0.5 * fl1_fx * tpz_xz_x_0[j] +
                             fl1_fx * fl1_fgb * tdz_xz_x_0[j];

            tlz_xxz_x_0[j] = pa_x[j] * tlz_xz_x_0[j] + 0.5 * fl1_fx * tlz_z_x_0[j] + 0.5 * fl1_fx * tlz_xz_0_0[j] - 0.5 * fl1_fx * tpy_xz_x_0[j] -
                             fl1_fx * fl1_fgb * tdy_xz_x_0[j];

            tlx_xxz_y_0[j] = pa_x[j] * tlx_xz_y_0[j] + 0.5 * fl1_fx * tlx_z_y_0[j];

            tly_xxz_y_0[j] = pa_x[j] * tly_xz_y_0[j] + 0.5 * fl1_fx * tly_z_y_0[j] + 0.5 * fl1_fx * tpz_xz_y_0[j] + fl1_fx * fl1_fgb * tdz_xz_y_0[j];

            tlz_xxz_y_0[j] = pa_x[j] * tlz_xz_y_0[j] + 0.5 * fl1_fx * tlz_z_y_0[j] - 0.5 * fl1_fx * tpy_xz_y_0[j] - fl1_fx * fl1_fgb * tdy_xz_y_0[j];

            tlx_xxz_z_0[j] = pa_x[j] * tlx_xz_z_0[j] + 0.5 * fl1_fx * tlx_z_z_0[j];

            tly_xxz_z_0[j] = pa_x[j] * tly_xz_z_0[j] + 0.5 * fl1_fx * tly_z_z_0[j] + 0.5 * fl1_fx * tpz_xz_z_0[j] + fl1_fx * fl1_fgb * tdz_xz_z_0[j];

            tlz_xxz_z_0[j] = pa_x[j] * tlz_xz_z_0[j] + 0.5 * fl1_fx * tlz_z_z_0[j] - 0.5 * fl1_fx * tpy_xz_z_0[j] - fl1_fx * fl1_fgb * tdy_xz_z_0[j];

            tlx_xyy_x_0[j] = pa_x[j] * tlx_yy_x_0[j] + 0.5 * fl1_fx * tlx_yy_0_0[j];

            tly_xyy_x_0[j] = pa_x[j] * tly_yy_x_0[j] + 0.5 * fl1_fx * tly_yy_0_0[j] + 0.5 * fl1_fx * tpz_yy_x_0[j] + fl1_fx * fl1_fgb * tdz_yy_x_0[j];

            tlz_xyy_x_0[j] = pa_x[j] * tlz_yy_x_0[j] + 0.5 * fl1_fx * tlz_yy_0_0[j] - 0.5 * fl1_fx * tpy_yy_x_0[j] - fl1_fx * fl1_fgb * tdy_yy_x_0[j];

            tlx_xyy_y_0[j] = pa_x[j] * tlx_yy_y_0[j];

            tly_xyy_y_0[j] = pa_x[j] * tly_yy_y_0[j] + 0.5 * fl1_fx * tpz_yy_y_0[j] + fl1_fx * fl1_fgb * tdz_yy_y_0[j];

            tlz_xyy_y_0[j] = pa_x[j] * tlz_yy_y_0[j] - 0.5 * fl1_fx * tpy_yy_y_0[j] - fl1_fx * fl1_fgb * tdy_yy_y_0[j];

            tlx_xyy_z_0[j] = pa_x[j] * tlx_yy_z_0[j];

            tly_xyy_z_0[j] = pa_x[j] * tly_yy_z_0[j] + 0.5 * fl1_fx * tpz_yy_z_0[j] + fl1_fx * fl1_fgb * tdz_yy_z_0[j];

            tlz_xyy_z_0[j] = pa_x[j] * tlz_yy_z_0[j] - 0.5 * fl1_fx * tpy_yy_z_0[j] - fl1_fx * fl1_fgb * tdy_yy_z_0[j];

            tlx_xyz_x_0[j] = pa_x[j] * tlx_yz_x_0[j] + 0.5 * fl1_fx * tlx_yz_0_0[j];

            tly_xyz_x_0[j] = pa_x[j] * tly_yz_x_0[j] + 0.5 * fl1_fx * tly_yz_0_0[j] + 0.5 * fl1_fx * tpz_yz_x_0[j] + fl1_fx * fl1_fgb * tdz_yz_x_0[j];

            tlz_xyz_x_0[j] = pa_x[j] * tlz_yz_x_0[j] + 0.5 * fl1_fx * tlz_yz_0_0[j] - 0.5 * fl1_fx * tpy_yz_x_0[j] - fl1_fx * fl1_fgb * tdy_yz_x_0[j];

            tlx_xyz_y_0[j] = pa_x[j] * tlx_yz_y_0[j];

            tly_xyz_y_0[j] = pa_x[j] * tly_yz_y_0[j] + 0.5 * fl1_fx * tpz_yz_y_0[j] + fl1_fx * fl1_fgb * tdz_yz_y_0[j];

            tlz_xyz_y_0[j] = pa_x[j] * tlz_yz_y_0[j] - 0.5 * fl1_fx * tpy_yz_y_0[j] - fl1_fx * fl1_fgb * tdy_yz_y_0[j];

            tlx_xyz_z_0[j] = pa_x[j] * tlx_yz_z_0[j];

            tly_xyz_z_0[j] = pa_x[j] * tly_yz_z_0[j] + 0.5 * fl1_fx * tpz_yz_z_0[j] + fl1_fx * fl1_fgb * tdz_yz_z_0[j];

            tlz_xyz_z_0[j] = pa_x[j] * tlz_yz_z_0[j] - 0.5 * fl1_fx * tpy_yz_z_0[j] - fl1_fx * fl1_fgb * tdy_yz_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFP_45_90(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9);

        auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10);

        auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11);

        auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12);

        auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13);

        auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14);

        auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15);

        auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16);

        auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17);

        auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3);

        auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4);

        auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5);

        auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6);

        auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7);

        auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8);

        auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tlx_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 3);

        auto tly_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 4);

        auto tly_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 5);

        auto tly_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tdx_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 9);

        auto tdz_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tdx_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 10);

        auto tdz_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tdx_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 11);

        auto tdz_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tdx_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 12);

        auto tdz_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tdx_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 13);

        auto tdz_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tdx_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 14);

        auto tdz_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15);

        auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16);

        auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17);

        auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17);

        // set up pointers to integrals

        auto tlx_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 15);

        auto tly_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tlz_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tlx_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 16);

        auto tly_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tlz_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tlx_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 17);

        auto tly_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tlz_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tlx_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 18);

        auto tly_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tlz_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tlx_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 19);

        auto tly_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tlz_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tlx_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 20);

        auto tly_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tlz_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tlx_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 21);

        auto tly_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tlz_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tlx_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 22);

        auto tly_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tlz_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tlx_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 23);

        auto tly_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tlz_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tlx_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 24);

        auto tly_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tlz_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tlx_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 25);

        auto tly_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tlz_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tlx_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 26);

        auto tly_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tlz_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tlx_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 27);

        auto tly_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tlz_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tlx_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 28);

        auto tly_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tlz_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tlx_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 29);

        auto tly_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tlz_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_yy_x_0, tdx_yy_y_0, tdx_yy_z_0, tdx_yz_x_0, \
                                     tdx_yz_y_0, tdx_yz_z_0, tdx_zz_x_0, tdx_zz_y_0, tdx_zz_z_0, tdy_zz_x_0, tdy_zz_y_0, \
                                     tdy_zz_z_0, tdz_yy_x_0, tdz_yy_y_0, tdz_yy_z_0, tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, \
                                     tdz_zz_x_0, tdz_zz_y_0, tdz_zz_z_0, tlx_xzz_x_0, tlx_xzz_y_0, tlx_xzz_z_0, \
                                     tlx_y_x_0, tlx_y_y_0, tlx_y_z_0, tlx_yy_0_0, tlx_yy_x_0, tlx_yy_y_0, tlx_yy_z_0, \
                                     tlx_yyy_x_0, tlx_yyy_y_0, tlx_yyy_z_0, tlx_yyz_x_0, tlx_yyz_y_0, tlx_yyz_z_0, \
                                     tlx_yz_0_0, tlx_yz_x_0, tlx_yz_y_0, tlx_yz_z_0, tlx_yzz_x_0, tlx_yzz_y_0, \
                                     tlx_yzz_z_0, tlx_z_x_0, tlx_z_y_0, tlx_z_z_0, tlx_zz_0_0, tlx_zz_x_0, tlx_zz_y_0, \
                                     tlx_zz_z_0, tlx_zzz_x_0, tlx_zzz_y_0, tlx_zzz_z_0, tly_xzz_x_0, tly_xzz_y_0, \
                                     tly_xzz_z_0, tly_y_x_0, tly_y_y_0, tly_y_z_0, tly_yy_0_0, tly_yy_x_0, tly_yy_y_0, \
                                     tly_yy_z_0, tly_yyy_x_0, tly_yyy_y_0, tly_yyy_z_0, tly_yyz_x_0, tly_yyz_y_0, \
                                     tly_yyz_z_0, tly_yz_0_0, tly_yz_x_0, tly_yz_y_0, tly_yz_z_0, tly_yzz_x_0, \
                                     tly_yzz_y_0, tly_yzz_z_0, tly_z_x_0, tly_z_y_0, tly_z_z_0, tly_zz_0_0, tly_zz_x_0, \
                                     tly_zz_y_0, tly_zz_z_0, tly_zzz_x_0, tly_zzz_y_0, tly_zzz_z_0, tlz_xzz_x_0, \
                                     tlz_xzz_y_0, tlz_xzz_z_0, tlz_y_x_0, tlz_y_y_0, tlz_y_z_0, tlz_yy_0_0, tlz_yy_x_0, \
                                     tlz_yy_y_0, tlz_yy_z_0, tlz_yyy_x_0, tlz_yyy_y_0, tlz_yyy_z_0, tlz_yyz_x_0, \
                                     tlz_yyz_y_0, tlz_yyz_z_0, tlz_yz_0_0, tlz_yz_x_0, tlz_yz_y_0, tlz_yz_z_0, \
                                     tlz_yzz_x_0, tlz_yzz_y_0, tlz_yzz_z_0, tlz_z_x_0, tlz_z_y_0, tlz_z_z_0, tlz_zz_0_0, \
                                     tlz_zz_x_0, tlz_zz_y_0, tlz_zz_z_0, tlz_zzz_x_0, tlz_zzz_y_0, tlz_zzz_z_0, \
                                     tpx_yy_x_0, tpx_yy_y_0, tpx_yy_z_0, tpx_yz_x_0, tpx_yz_y_0, tpx_yz_z_0, tpx_zz_x_0, \
                                     tpx_zz_y_0, tpx_zz_z_0, tpy_zz_x_0, tpy_zz_y_0, tpy_zz_z_0, tpz_yy_x_0, tpz_yy_y_0, \
                                     tpz_yy_z_0, tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0, tpz_zz_x_0, tpz_zz_y_0, tpz_zz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xzz_x_0[j] = pa_x[j] * tlx_zz_x_0[j] + 0.5 * fl1_fx * tlx_zz_0_0[j];

            tly_xzz_x_0[j] = pa_x[j] * tly_zz_x_0[j] + 0.5 * fl1_fx * tly_zz_0_0[j] + 0.5 * fl1_fx * tpz_zz_x_0[j] + fl1_fx * fl1_fgb * tdz_zz_x_0[j];

            tlz_xzz_x_0[j] = pa_x[j] * tlz_zz_x_0[j] + 0.5 * fl1_fx * tlz_zz_0_0[j] - 0.5 * fl1_fx * tpy_zz_x_0[j] - fl1_fx * fl1_fgb * tdy_zz_x_0[j];

            tlx_xzz_y_0[j] = pa_x[j] * tlx_zz_y_0[j];

            tly_xzz_y_0[j] = pa_x[j] * tly_zz_y_0[j] + 0.5 * fl1_fx * tpz_zz_y_0[j] + fl1_fx * fl1_fgb * tdz_zz_y_0[j];

            tlz_xzz_y_0[j] = pa_x[j] * tlz_zz_y_0[j] - 0.5 * fl1_fx * tpy_zz_y_0[j] - fl1_fx * fl1_fgb * tdy_zz_y_0[j];

            tlx_xzz_z_0[j] = pa_x[j] * tlx_zz_z_0[j];

            tly_xzz_z_0[j] = pa_x[j] * tly_zz_z_0[j] + 0.5 * fl1_fx * tpz_zz_z_0[j] + fl1_fx * fl1_fgb * tdz_zz_z_0[j];

            tlz_xzz_z_0[j] = pa_x[j] * tlz_zz_z_0[j] - 0.5 * fl1_fx * tpy_zz_z_0[j] - fl1_fx * fl1_fgb * tdy_zz_z_0[j];

            tlx_yyy_x_0[j] = pa_y[j] * tlx_yy_x_0[j] + fl1_fx * tlx_y_x_0[j] - 0.5 * fl1_fx * tpz_yy_x_0[j] - fl1_fx * fl1_fgb * tdz_yy_x_0[j];

            tly_yyy_x_0[j] = pa_y[j] * tly_yy_x_0[j] + fl1_fx * tly_y_x_0[j];

            tlz_yyy_x_0[j] = pa_y[j] * tlz_yy_x_0[j] + fl1_fx * tlz_y_x_0[j] + 0.5 * fl1_fx * tpx_yy_x_0[j] + fl1_fx * fl1_fgb * tdx_yy_x_0[j];

            tlx_yyy_y_0[j] = pa_y[j] * tlx_yy_y_0[j] + fl1_fx * tlx_y_y_0[j] + 0.5 * fl1_fx * tlx_yy_0_0[j] - 0.5 * fl1_fx * tpz_yy_y_0[j] -
                             fl1_fx * fl1_fgb * tdz_yy_y_0[j];

            tly_yyy_y_0[j] = pa_y[j] * tly_yy_y_0[j] + fl1_fx * tly_y_y_0[j] + 0.5 * fl1_fx * tly_yy_0_0[j];

            tlz_yyy_y_0[j] = pa_y[j] * tlz_yy_y_0[j] + fl1_fx * tlz_y_y_0[j] + 0.5 * fl1_fx * tlz_yy_0_0[j] + 0.5 * fl1_fx * tpx_yy_y_0[j] +
                             fl1_fx * fl1_fgb * tdx_yy_y_0[j];

            tlx_yyy_z_0[j] = pa_y[j] * tlx_yy_z_0[j] + fl1_fx * tlx_y_z_0[j] - 0.5 * fl1_fx * tpz_yy_z_0[j] - fl1_fx * fl1_fgb * tdz_yy_z_0[j];

            tly_yyy_z_0[j] = pa_y[j] * tly_yy_z_0[j] + fl1_fx * tly_y_z_0[j];

            tlz_yyy_z_0[j] = pa_y[j] * tlz_yy_z_0[j] + fl1_fx * tlz_y_z_0[j] + 0.5 * fl1_fx * tpx_yy_z_0[j] + fl1_fx * fl1_fgb * tdx_yy_z_0[j];

            tlx_yyz_x_0[j] = pa_y[j] * tlx_yz_x_0[j] + 0.5 * fl1_fx * tlx_z_x_0[j] - 0.5 * fl1_fx * tpz_yz_x_0[j] - fl1_fx * fl1_fgb * tdz_yz_x_0[j];

            tly_yyz_x_0[j] = pa_y[j] * tly_yz_x_0[j] + 0.5 * fl1_fx * tly_z_x_0[j];

            tlz_yyz_x_0[j] = pa_y[j] * tlz_yz_x_0[j] + 0.5 * fl1_fx * tlz_z_x_0[j] + 0.5 * fl1_fx * tpx_yz_x_0[j] + fl1_fx * fl1_fgb * tdx_yz_x_0[j];

            tlx_yyz_y_0[j] = pa_y[j] * tlx_yz_y_0[j] + 0.5 * fl1_fx * tlx_z_y_0[j] + 0.5 * fl1_fx * tlx_yz_0_0[j] - 0.5 * fl1_fx * tpz_yz_y_0[j] -
                             fl1_fx * fl1_fgb * tdz_yz_y_0[j];

            tly_yyz_y_0[j] = pa_y[j] * tly_yz_y_0[j] + 0.5 * fl1_fx * tly_z_y_0[j] + 0.5 * fl1_fx * tly_yz_0_0[j];

            tlz_yyz_y_0[j] = pa_y[j] * tlz_yz_y_0[j] + 0.5 * fl1_fx * tlz_z_y_0[j] + 0.5 * fl1_fx * tlz_yz_0_0[j] + 0.5 * fl1_fx * tpx_yz_y_0[j] +
                             fl1_fx * fl1_fgb * tdx_yz_y_0[j];

            tlx_yyz_z_0[j] = pa_y[j] * tlx_yz_z_0[j] + 0.5 * fl1_fx * tlx_z_z_0[j] - 0.5 * fl1_fx * tpz_yz_z_0[j] - fl1_fx * fl1_fgb * tdz_yz_z_0[j];

            tly_yyz_z_0[j] = pa_y[j] * tly_yz_z_0[j] + 0.5 * fl1_fx * tly_z_z_0[j];

            tlz_yyz_z_0[j] = pa_y[j] * tlz_yz_z_0[j] + 0.5 * fl1_fx * tlz_z_z_0[j] + 0.5 * fl1_fx * tpx_yz_z_0[j] + fl1_fx * fl1_fgb * tdx_yz_z_0[j];

            tlx_yzz_x_0[j] = pa_y[j] * tlx_zz_x_0[j] - 0.5 * fl1_fx * tpz_zz_x_0[j] - fl1_fx * fl1_fgb * tdz_zz_x_0[j];

            tly_yzz_x_0[j] = pa_y[j] * tly_zz_x_0[j];

            tlz_yzz_x_0[j] = pa_y[j] * tlz_zz_x_0[j] + 0.5 * fl1_fx * tpx_zz_x_0[j] + fl1_fx * fl1_fgb * tdx_zz_x_0[j];

            tlx_yzz_y_0[j] = pa_y[j] * tlx_zz_y_0[j] + 0.5 * fl1_fx * tlx_zz_0_0[j] - 0.5 * fl1_fx * tpz_zz_y_0[j] - fl1_fx * fl1_fgb * tdz_zz_y_0[j];

            tly_yzz_y_0[j] = pa_y[j] * tly_zz_y_0[j] + 0.5 * fl1_fx * tly_zz_0_0[j];

            tlz_yzz_y_0[j] = pa_y[j] * tlz_zz_y_0[j] + 0.5 * fl1_fx * tlz_zz_0_0[j] + 0.5 * fl1_fx * tpx_zz_y_0[j] + fl1_fx * fl1_fgb * tdx_zz_y_0[j];

            tlx_yzz_z_0[j] = pa_y[j] * tlx_zz_z_0[j] - 0.5 * fl1_fx * tpz_zz_z_0[j] - fl1_fx * fl1_fgb * tdz_zz_z_0[j];

            tly_yzz_z_0[j] = pa_y[j] * tly_zz_z_0[j];

            tlz_yzz_z_0[j] = pa_y[j] * tlz_zz_z_0[j] + 0.5 * fl1_fx * tpx_zz_z_0[j] + fl1_fx * fl1_fgb * tdx_zz_z_0[j];

            tlx_zzz_x_0[j] = pa_z[j] * tlx_zz_x_0[j] + fl1_fx * tlx_z_x_0[j] + 0.5 * fl1_fx * tpy_zz_x_0[j] + fl1_fx * fl1_fgb * tdy_zz_x_0[j];

            tly_zzz_x_0[j] = pa_z[j] * tly_zz_x_0[j] + fl1_fx * tly_z_x_0[j] - 0.5 * fl1_fx * tpx_zz_x_0[j] - fl1_fx * fl1_fgb * tdx_zz_x_0[j];

            tlz_zzz_x_0[j] = pa_z[j] * tlz_zz_x_0[j] + fl1_fx * tlz_z_x_0[j];

            tlx_zzz_y_0[j] = pa_z[j] * tlx_zz_y_0[j] + fl1_fx * tlx_z_y_0[j] + 0.5 * fl1_fx * tpy_zz_y_0[j] + fl1_fx * fl1_fgb * tdy_zz_y_0[j];

            tly_zzz_y_0[j] = pa_z[j] * tly_zz_y_0[j] + fl1_fx * tly_z_y_0[j] - 0.5 * fl1_fx * tpx_zz_y_0[j] - fl1_fx * fl1_fgb * tdx_zz_y_0[j];

            tlz_zzz_y_0[j] = pa_z[j] * tlz_zz_y_0[j] + fl1_fx * tlz_z_y_0[j];

            tlx_zzz_z_0[j] = pa_z[j] * tlx_zz_z_0[j] + fl1_fx * tlx_z_z_0[j] + 0.5 * fl1_fx * tlx_zz_0_0[j] + 0.5 * fl1_fx * tpy_zz_z_0[j] +
                             fl1_fx * fl1_fgb * tdy_zz_z_0[j];

            tly_zzz_z_0[j] = pa_z[j] * tly_zz_z_0[j] + fl1_fx * tly_z_z_0[j] + 0.5 * fl1_fx * tly_zz_0_0[j] - 0.5 * fl1_fx * tpx_zz_z_0[j] -
                             fl1_fx * fl1_fgb * tdx_zz_z_0[j];

            tlz_zzz_z_0[j] = pa_z[j] * tlz_zz_z_0[j] + fl1_fx * tlz_z_z_0[j] + 0.5 * fl1_fx * tlz_zz_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPG(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForPG_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForPG_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForPG_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForPG_0_45(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx);

        auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx);

        auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx);

        auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1);

        auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2);

        auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3);

        auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4);

        auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5);

        auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6);

        auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7);

        auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8);

        auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9);

        auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10);

        auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11);

        auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12);

        auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13);

        auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14);

        auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tdy_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx);

        auto tdz_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx);

        auto tdy_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tdz_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tdy_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tdz_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tdy_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tdz_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tdy_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tdz_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tdy_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tdz_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tdy_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tdz_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tdy_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tdz_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tdy_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tdz_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tdy_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tdz_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tdy_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tdz_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tdy_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tdz_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tdy_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tdz_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tdy_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tdz_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tdy_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tdz_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 14);

        // set up pointers to integrals

        auto tlx_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx);

        auto tly_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx);

        auto tlz_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx);

        auto tlx_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 1);

        auto tly_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 1);

        auto tlz_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 1);

        auto tlx_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 2);

        auto tly_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 2);

        auto tlz_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 2);

        auto tlx_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 3);

        auto tly_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 3);

        auto tlz_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 3);

        auto tlx_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 4);

        auto tly_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 4);

        auto tlz_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 4);

        auto tlx_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 5);

        auto tly_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 5);

        auto tlz_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 5);

        auto tlx_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 6);

        auto tly_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 6);

        auto tlz_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 6);

        auto tlx_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 7);

        auto tly_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 7);

        auto tlz_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 7);

        auto tlx_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 8);

        auto tly_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 8);

        auto tlz_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 8);

        auto tlx_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 9);

        auto tly_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 9);

        auto tlz_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 9);

        auto tlx_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 10);

        auto tly_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 10);

        auto tlz_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 10);

        auto tlx_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 11);

        auto tly_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 11);

        auto tlz_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 11);

        auto tlx_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 12);

        auto tly_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 12);

        auto tlz_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 12);

        auto tlx_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 13);

        auto tly_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 13);

        auto tlz_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 13);

        auto tlx_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 14);

        auto tly_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 14);

        auto tlz_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxyy_0, \
                                     tdy_0_xxyz_0, tdy_0_xxzz_0, tdy_0_xyyy_0, tdy_0_xyyz_0, tdy_0_xyzz_0, tdy_0_xzzz_0, \
                                     tdy_0_yyyy_0, tdy_0_yyyz_0, tdy_0_yyzz_0, tdy_0_yzzz_0, tdy_0_zzzz_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxzz_0, tdz_0_xyyy_0, \
                                     tdz_0_xyyz_0, tdz_0_xyzz_0, tdz_0_xzzz_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyzz_0, \
                                     tdz_0_yzzz_0, tdz_0_zzzz_0, tlx_0_xxx_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, \
                                     tlx_0_xxy_0, tlx_0_xxyy_0, tlx_0_xxyz_0, tlx_0_xxz_0, tlx_0_xxzz_0, tlx_0_xyy_0, \
                                     tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyz_0, tlx_0_xyzz_0, tlx_0_xzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyy_0, tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyz_0, tlx_0_yyzz_0, tlx_0_yzz_0, \
                                     tlx_0_yzzz_0, tlx_0_zzz_0, tlx_0_zzzz_0, tlx_x_xxxx_0, tlx_x_xxxy_0, tlx_x_xxxz_0, \
                                     tlx_x_xxyy_0, tlx_x_xxyz_0, tlx_x_xxzz_0, tlx_x_xyyy_0, tlx_x_xyyz_0, tlx_x_xyzz_0, \
                                     tlx_x_xzzz_0, tlx_x_yyyy_0, tlx_x_yyyz_0, tlx_x_yyzz_0, tlx_x_yzzz_0, tlx_x_zzzz_0, \
                                     tly_0_xxx_0, tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxy_0, tly_0_xxyy_0, \
                                     tly_0_xxyz_0, tly_0_xxz_0, tly_0_xxzz_0, tly_0_xyy_0, tly_0_xyyy_0, tly_0_xyyz_0, \
                                     tly_0_xyz_0, tly_0_xyzz_0, tly_0_xzz_0, tly_0_xzzz_0, tly_0_yyy_0, tly_0_yyyy_0, \
                                     tly_0_yyyz_0, tly_0_yyz_0, tly_0_yyzz_0, tly_0_yzz_0, tly_0_yzzz_0, tly_0_zzz_0, \
                                     tly_0_zzzz_0, tly_x_xxxx_0, tly_x_xxxy_0, tly_x_xxxz_0, tly_x_xxyy_0, tly_x_xxyz_0, \
                                     tly_x_xxzz_0, tly_x_xyyy_0, tly_x_xyyz_0, tly_x_xyzz_0, tly_x_xzzz_0, tly_x_yyyy_0, \
                                     tly_x_yyyz_0, tly_x_yyzz_0, tly_x_yzzz_0, tly_x_zzzz_0, tlz_0_xxx_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxy_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxz_0, \
                                     tlz_0_xxzz_0, tlz_0_xyy_0, tlz_0_xyyy_0, tlz_0_xyyz_0, tlz_0_xyz_0, tlz_0_xyzz_0, \
                                     tlz_0_xzz_0, tlz_0_xzzz_0, tlz_0_yyy_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyz_0, \
                                     tlz_0_yyzz_0, tlz_0_yzz_0, tlz_0_yzzz_0, tlz_0_zzz_0, tlz_0_zzzz_0, tlz_x_xxxx_0, \
                                     tlz_x_xxxy_0, tlz_x_xxxz_0, tlz_x_xxyy_0, tlz_x_xxyz_0, tlz_x_xxzz_0, tlz_x_xyyy_0, \
                                     tlz_x_xyyz_0, tlz_x_xyzz_0, tlz_x_xzzz_0, tlz_x_yyyy_0, tlz_x_yyyz_0, tlz_x_yyzz_0, \
                                     tlz_x_yzzz_0, tlz_x_zzzz_0, tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxyy_0, \
                                     tpy_0_xxyz_0, tpy_0_xxzz_0, tpy_0_xyyy_0, tpy_0_xyyz_0, tpy_0_xyzz_0, tpy_0_xzzz_0, \
                                     tpy_0_yyyy_0, tpy_0_yyyz_0, tpy_0_yyzz_0, tpy_0_yzzz_0, tpy_0_zzzz_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxzz_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyzz_0, tpz_0_xzzz_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yzzz_0, tpz_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_x_xxxx_0[j] = pa_x[j] * tlx_0_xxxx_0[j] + 2.0 * fl1_fx * tlx_0_xxx_0[j];

            tly_x_xxxx_0[j] =
                pa_x[j] * tly_0_xxxx_0[j] + 2.0 * fl1_fx * tly_0_xxx_0[j] + 0.5 * fl1_fx * tpz_0_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_0_xxxx_0[j];

            tlz_x_xxxx_0[j] =
                pa_x[j] * tlz_0_xxxx_0[j] + 2.0 * fl1_fx * tlz_0_xxx_0[j] - 0.5 * fl1_fx * tpy_0_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_0_xxxx_0[j];

            tlx_x_xxxy_0[j] = pa_x[j] * tlx_0_xxxy_0[j] + 1.5 * fl1_fx * tlx_0_xxy_0[j];

            tly_x_xxxy_0[j] =
                pa_x[j] * tly_0_xxxy_0[j] + 1.5 * fl1_fx * tly_0_xxy_0[j] + 0.5 * fl1_fx * tpz_0_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_0_xxxy_0[j];

            tlz_x_xxxy_0[j] =
                pa_x[j] * tlz_0_xxxy_0[j] + 1.5 * fl1_fx * tlz_0_xxy_0[j] - 0.5 * fl1_fx * tpy_0_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_0_xxxy_0[j];

            tlx_x_xxxz_0[j] = pa_x[j] * tlx_0_xxxz_0[j] + 1.5 * fl1_fx * tlx_0_xxz_0[j];

            tly_x_xxxz_0[j] =
                pa_x[j] * tly_0_xxxz_0[j] + 1.5 * fl1_fx * tly_0_xxz_0[j] + 0.5 * fl1_fx * tpz_0_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_0_xxxz_0[j];

            tlz_x_xxxz_0[j] =
                pa_x[j] * tlz_0_xxxz_0[j] + 1.5 * fl1_fx * tlz_0_xxz_0[j] - 0.5 * fl1_fx * tpy_0_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_0_xxxz_0[j];

            tlx_x_xxyy_0[j] = pa_x[j] * tlx_0_xxyy_0[j] + fl1_fx * tlx_0_xyy_0[j];

            tly_x_xxyy_0[j] =
                pa_x[j] * tly_0_xxyy_0[j] + fl1_fx * tly_0_xyy_0[j] + 0.5 * fl1_fx * tpz_0_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_0_xxyy_0[j];

            tlz_x_xxyy_0[j] =
                pa_x[j] * tlz_0_xxyy_0[j] + fl1_fx * tlz_0_xyy_0[j] - 0.5 * fl1_fx * tpy_0_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_0_xxyy_0[j];

            tlx_x_xxyz_0[j] = pa_x[j] * tlx_0_xxyz_0[j] + fl1_fx * tlx_0_xyz_0[j];

            tly_x_xxyz_0[j] =
                pa_x[j] * tly_0_xxyz_0[j] + fl1_fx * tly_0_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_0_xxyz_0[j];

            tlz_x_xxyz_0[j] =
                pa_x[j] * tlz_0_xxyz_0[j] + fl1_fx * tlz_0_xyz_0[j] - 0.5 * fl1_fx * tpy_0_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_0_xxyz_0[j];

            tlx_x_xxzz_0[j] = pa_x[j] * tlx_0_xxzz_0[j] + fl1_fx * tlx_0_xzz_0[j];

            tly_x_xxzz_0[j] =
                pa_x[j] * tly_0_xxzz_0[j] + fl1_fx * tly_0_xzz_0[j] + 0.5 * fl1_fx * tpz_0_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_0_xxzz_0[j];

            tlz_x_xxzz_0[j] =
                pa_x[j] * tlz_0_xxzz_0[j] + fl1_fx * tlz_0_xzz_0[j] - 0.5 * fl1_fx * tpy_0_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_0_xxzz_0[j];

            tlx_x_xyyy_0[j] = pa_x[j] * tlx_0_xyyy_0[j] + 0.5 * fl1_fx * tlx_0_yyy_0[j];

            tly_x_xyyy_0[j] =
                pa_x[j] * tly_0_xyyy_0[j] + 0.5 * fl1_fx * tly_0_yyy_0[j] + 0.5 * fl1_fx * tpz_0_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_0_xyyy_0[j];

            tlz_x_xyyy_0[j] =
                pa_x[j] * tlz_0_xyyy_0[j] + 0.5 * fl1_fx * tlz_0_yyy_0[j] - 0.5 * fl1_fx * tpy_0_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_0_xyyy_0[j];

            tlx_x_xyyz_0[j] = pa_x[j] * tlx_0_xyyz_0[j] + 0.5 * fl1_fx * tlx_0_yyz_0[j];

            tly_x_xyyz_0[j] =
                pa_x[j] * tly_0_xyyz_0[j] + 0.5 * fl1_fx * tly_0_yyz_0[j] + 0.5 * fl1_fx * tpz_0_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_0_xyyz_0[j];

            tlz_x_xyyz_0[j] =
                pa_x[j] * tlz_0_xyyz_0[j] + 0.5 * fl1_fx * tlz_0_yyz_0[j] - 0.5 * fl1_fx * tpy_0_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_0_xyyz_0[j];

            tlx_x_xyzz_0[j] = pa_x[j] * tlx_0_xyzz_0[j] + 0.5 * fl1_fx * tlx_0_yzz_0[j];

            tly_x_xyzz_0[j] =
                pa_x[j] * tly_0_xyzz_0[j] + 0.5 * fl1_fx * tly_0_yzz_0[j] + 0.5 * fl1_fx * tpz_0_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_0_xyzz_0[j];

            tlz_x_xyzz_0[j] =
                pa_x[j] * tlz_0_xyzz_0[j] + 0.5 * fl1_fx * tlz_0_yzz_0[j] - 0.5 * fl1_fx * tpy_0_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_0_xyzz_0[j];

            tlx_x_xzzz_0[j] = pa_x[j] * tlx_0_xzzz_0[j] + 0.5 * fl1_fx * tlx_0_zzz_0[j];

            tly_x_xzzz_0[j] =
                pa_x[j] * tly_0_xzzz_0[j] + 0.5 * fl1_fx * tly_0_zzz_0[j] + 0.5 * fl1_fx * tpz_0_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_0_xzzz_0[j];

            tlz_x_xzzz_0[j] =
                pa_x[j] * tlz_0_xzzz_0[j] + 0.5 * fl1_fx * tlz_0_zzz_0[j] - 0.5 * fl1_fx * tpy_0_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_0_xzzz_0[j];

            tlx_x_yyyy_0[j] = pa_x[j] * tlx_0_yyyy_0[j];

            tly_x_yyyy_0[j] = pa_x[j] * tly_0_yyyy_0[j] + 0.5 * fl1_fx * tpz_0_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_0_yyyy_0[j];

            tlz_x_yyyy_0[j] = pa_x[j] * tlz_0_yyyy_0[j] - 0.5 * fl1_fx * tpy_0_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_0_yyyy_0[j];

            tlx_x_yyyz_0[j] = pa_x[j] * tlx_0_yyyz_0[j];

            tly_x_yyyz_0[j] = pa_x[j] * tly_0_yyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_0_yyyz_0[j];

            tlz_x_yyyz_0[j] = pa_x[j] * tlz_0_yyyz_0[j] - 0.5 * fl1_fx * tpy_0_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_0_yyyz_0[j];

            tlx_x_yyzz_0[j] = pa_x[j] * tlx_0_yyzz_0[j];

            tly_x_yyzz_0[j] = pa_x[j] * tly_0_yyzz_0[j] + 0.5 * fl1_fx * tpz_0_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_0_yyzz_0[j];

            tlz_x_yyzz_0[j] = pa_x[j] * tlz_0_yyzz_0[j] - 0.5 * fl1_fx * tpy_0_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_0_yyzz_0[j];

            tlx_x_yzzz_0[j] = pa_x[j] * tlx_0_yzzz_0[j];

            tly_x_yzzz_0[j] = pa_x[j] * tly_0_yzzz_0[j] + 0.5 * fl1_fx * tpz_0_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_0_yzzz_0[j];

            tlz_x_yzzz_0[j] = pa_x[j] * tlz_0_yzzz_0[j] - 0.5 * fl1_fx * tpy_0_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_0_yzzz_0[j];

            tlx_x_zzzz_0[j] = pa_x[j] * tlx_0_zzzz_0[j];

            tly_x_zzzz_0[j] = pa_x[j] * tly_0_zzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_0_zzzz_0[j];

            tlz_x_zzzz_0[j] = pa_x[j] * tlz_0_zzzz_0[j] - 0.5 * fl1_fx * tpy_0_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPG_45_90(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx);

        auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx);

        auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx);

        auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1);

        auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2);

        auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3);

        auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4);

        auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5);

        auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6);

        auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7);

        auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8);

        auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9);

        auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10);

        auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11);

        auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12);

        auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13);

        auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14);

        auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tdx_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx);

        auto tdz_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx);

        auto tdx_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 1);

        auto tdz_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tdx_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 2);

        auto tdz_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tdx_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 3);

        auto tdz_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tdx_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 4);

        auto tdz_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tdx_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 5);

        auto tdz_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tdx_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 6);

        auto tdz_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tdx_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 7);

        auto tdz_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tdx_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 8);

        auto tdz_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tdx_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 9);

        auto tdz_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tdx_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 10);

        auto tdz_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tdx_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 11);

        auto tdz_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tdx_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 12);

        auto tdz_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tdx_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 13);

        auto tdz_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tdx_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 14);

        auto tdz_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 14);

        // set up pointers to integrals

        auto tlx_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 15);

        auto tly_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tlz_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tlx_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 16);

        auto tly_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tlz_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tlx_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 17);

        auto tly_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tlz_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tlx_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 18);

        auto tly_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tlz_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tlx_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 19);

        auto tly_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tlz_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tlx_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 20);

        auto tly_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tlz_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tlx_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 21);

        auto tly_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tlz_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tlx_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 22);

        auto tly_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tlz_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tlx_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 23);

        auto tly_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tlz_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tlx_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 24);

        auto tly_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tlz_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tlx_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 25);

        auto tly_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tlz_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tlx_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 26);

        auto tly_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tlz_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tlx_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 27);

        auto tly_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tlz_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tlx_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 28);

        auto tly_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tlz_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tlx_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 29);

        auto tly_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tlz_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, tdx_0_xxyy_0, \
                                     tdx_0_xxyz_0, tdx_0_xxzz_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyzz_0, tdx_0_yzzz_0, tdx_0_zzzz_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxzz_0, tdz_0_xyyy_0, \
                                     tdz_0_xyyz_0, tdz_0_xyzz_0, tdz_0_xzzz_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyzz_0, \
                                     tdz_0_yzzz_0, tdz_0_zzzz_0, tlx_0_xxx_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, \
                                     tlx_0_xxy_0, tlx_0_xxyy_0, tlx_0_xxyz_0, tlx_0_xxz_0, tlx_0_xxzz_0, tlx_0_xyy_0, \
                                     tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyz_0, tlx_0_xyzz_0, tlx_0_xzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyy_0, tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyz_0, tlx_0_yyzz_0, tlx_0_yzz_0, \
                                     tlx_0_yzzz_0, tlx_0_zzz_0, tlx_0_zzzz_0, tlx_y_xxxx_0, tlx_y_xxxy_0, tlx_y_xxxz_0, \
                                     tlx_y_xxyy_0, tlx_y_xxyz_0, tlx_y_xxzz_0, tlx_y_xyyy_0, tlx_y_xyyz_0, tlx_y_xyzz_0, \
                                     tlx_y_xzzz_0, tlx_y_yyyy_0, tlx_y_yyyz_0, tlx_y_yyzz_0, tlx_y_yzzz_0, tlx_y_zzzz_0, \
                                     tly_0_xxx_0, tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxy_0, tly_0_xxyy_0, \
                                     tly_0_xxyz_0, tly_0_xxz_0, tly_0_xxzz_0, tly_0_xyy_0, tly_0_xyyy_0, tly_0_xyyz_0, \
                                     tly_0_xyz_0, tly_0_xyzz_0, tly_0_xzz_0, tly_0_xzzz_0, tly_0_yyy_0, tly_0_yyyy_0, \
                                     tly_0_yyyz_0, tly_0_yyz_0, tly_0_yyzz_0, tly_0_yzz_0, tly_0_yzzz_0, tly_0_zzz_0, \
                                     tly_0_zzzz_0, tly_y_xxxx_0, tly_y_xxxy_0, tly_y_xxxz_0, tly_y_xxyy_0, tly_y_xxyz_0, \
                                     tly_y_xxzz_0, tly_y_xyyy_0, tly_y_xyyz_0, tly_y_xyzz_0, tly_y_xzzz_0, tly_y_yyyy_0, \
                                     tly_y_yyyz_0, tly_y_yyzz_0, tly_y_yzzz_0, tly_y_zzzz_0, tlz_0_xxx_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxy_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxz_0, \
                                     tlz_0_xxzz_0, tlz_0_xyy_0, tlz_0_xyyy_0, tlz_0_xyyz_0, tlz_0_xyz_0, tlz_0_xyzz_0, \
                                     tlz_0_xzz_0, tlz_0_xzzz_0, tlz_0_yyy_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyz_0, \
                                     tlz_0_yyzz_0, tlz_0_yzz_0, tlz_0_yzzz_0, tlz_0_zzz_0, tlz_0_zzzz_0, tlz_y_xxxx_0, \
                                     tlz_y_xxxy_0, tlz_y_xxxz_0, tlz_y_xxyy_0, tlz_y_xxyz_0, tlz_y_xxzz_0, tlz_y_xyyy_0, \
                                     tlz_y_xyyz_0, tlz_y_xyzz_0, tlz_y_xzzz_0, tlz_y_yyyy_0, tlz_y_yyyz_0, tlz_y_yyzz_0, \
                                     tlz_y_yzzz_0, tlz_y_zzzz_0, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxyy_0, \
                                     tpx_0_xxyz_0, tpx_0_xxzz_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyzz_0, tpx_0_yzzz_0, tpx_0_zzzz_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxzz_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyzz_0, tpz_0_xzzz_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yzzz_0, tpz_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_y_xxxx_0[j] = pa_y[j] * tlx_0_xxxx_0[j] - 0.5 * fl1_fx * tpz_0_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_0_xxxx_0[j];

            tly_y_xxxx_0[j] = pa_y[j] * tly_0_xxxx_0[j];

            tlz_y_xxxx_0[j] = pa_y[j] * tlz_0_xxxx_0[j] + 0.5 * fl1_fx * tpx_0_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_0_xxxx_0[j];

            tlx_y_xxxy_0[j] =
                pa_y[j] * tlx_0_xxxy_0[j] + 0.5 * fl1_fx * tlx_0_xxx_0[j] - 0.5 * fl1_fx * tpz_0_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_0_xxxy_0[j];

            tly_y_xxxy_0[j] = pa_y[j] * tly_0_xxxy_0[j] + 0.5 * fl1_fx * tly_0_xxx_0[j];

            tlz_y_xxxy_0[j] =
                pa_y[j] * tlz_0_xxxy_0[j] + 0.5 * fl1_fx * tlz_0_xxx_0[j] + 0.5 * fl1_fx * tpx_0_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_0_xxxy_0[j];

            tlx_y_xxxz_0[j] = pa_y[j] * tlx_0_xxxz_0[j] - 0.5 * fl1_fx * tpz_0_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_0_xxxz_0[j];

            tly_y_xxxz_0[j] = pa_y[j] * tly_0_xxxz_0[j];

            tlz_y_xxxz_0[j] = pa_y[j] * tlz_0_xxxz_0[j] + 0.5 * fl1_fx * tpx_0_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_0_xxxz_0[j];

            tlx_y_xxyy_0[j] =
                pa_y[j] * tlx_0_xxyy_0[j] + fl1_fx * tlx_0_xxy_0[j] - 0.5 * fl1_fx * tpz_0_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_0_xxyy_0[j];

            tly_y_xxyy_0[j] = pa_y[j] * tly_0_xxyy_0[j] + fl1_fx * tly_0_xxy_0[j];

            tlz_y_xxyy_0[j] =
                pa_y[j] * tlz_0_xxyy_0[j] + fl1_fx * tlz_0_xxy_0[j] + 0.5 * fl1_fx * tpx_0_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_0_xxyy_0[j];

            tlx_y_xxyz_0[j] =
                pa_y[j] * tlx_0_xxyz_0[j] + 0.5 * fl1_fx * tlx_0_xxz_0[j] - 0.5 * fl1_fx * tpz_0_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_0_xxyz_0[j];

            tly_y_xxyz_0[j] = pa_y[j] * tly_0_xxyz_0[j] + 0.5 * fl1_fx * tly_0_xxz_0[j];

            tlz_y_xxyz_0[j] =
                pa_y[j] * tlz_0_xxyz_0[j] + 0.5 * fl1_fx * tlz_0_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_0_xxyz_0[j];

            tlx_y_xxzz_0[j] = pa_y[j] * tlx_0_xxzz_0[j] - 0.5 * fl1_fx * tpz_0_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_0_xxzz_0[j];

            tly_y_xxzz_0[j] = pa_y[j] * tly_0_xxzz_0[j];

            tlz_y_xxzz_0[j] = pa_y[j] * tlz_0_xxzz_0[j] + 0.5 * fl1_fx * tpx_0_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_0_xxzz_0[j];

            tlx_y_xyyy_0[j] =
                pa_y[j] * tlx_0_xyyy_0[j] + 1.5 * fl1_fx * tlx_0_xyy_0[j] - 0.5 * fl1_fx * tpz_0_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_0_xyyy_0[j];

            tly_y_xyyy_0[j] = pa_y[j] * tly_0_xyyy_0[j] + 1.5 * fl1_fx * tly_0_xyy_0[j];

            tlz_y_xyyy_0[j] =
                pa_y[j] * tlz_0_xyyy_0[j] + 1.5 * fl1_fx * tlz_0_xyy_0[j] + 0.5 * fl1_fx * tpx_0_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_0_xyyy_0[j];

            tlx_y_xyyz_0[j] =
                pa_y[j] * tlx_0_xyyz_0[j] + fl1_fx * tlx_0_xyz_0[j] - 0.5 * fl1_fx * tpz_0_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_0_xyyz_0[j];

            tly_y_xyyz_0[j] = pa_y[j] * tly_0_xyyz_0[j] + fl1_fx * tly_0_xyz_0[j];

            tlz_y_xyyz_0[j] =
                pa_y[j] * tlz_0_xyyz_0[j] + fl1_fx * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_0_xyyz_0[j];

            tlx_y_xyzz_0[j] =
                pa_y[j] * tlx_0_xyzz_0[j] + 0.5 * fl1_fx * tlx_0_xzz_0[j] - 0.5 * fl1_fx * tpz_0_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_0_xyzz_0[j];

            tly_y_xyzz_0[j] = pa_y[j] * tly_0_xyzz_0[j] + 0.5 * fl1_fx * tly_0_xzz_0[j];

            tlz_y_xyzz_0[j] =
                pa_y[j] * tlz_0_xyzz_0[j] + 0.5 * fl1_fx * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tpx_0_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_0_xyzz_0[j];

            tlx_y_xzzz_0[j] = pa_y[j] * tlx_0_xzzz_0[j] - 0.5 * fl1_fx * tpz_0_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_0_xzzz_0[j];

            tly_y_xzzz_0[j] = pa_y[j] * tly_0_xzzz_0[j];

            tlz_y_xzzz_0[j] = pa_y[j] * tlz_0_xzzz_0[j] + 0.5 * fl1_fx * tpx_0_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_0_xzzz_0[j];

            tlx_y_yyyy_0[j] =
                pa_y[j] * tlx_0_yyyy_0[j] + 2.0 * fl1_fx * tlx_0_yyy_0[j] - 0.5 * fl1_fx * tpz_0_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_0_yyyy_0[j];

            tly_y_yyyy_0[j] = pa_y[j] * tly_0_yyyy_0[j] + 2.0 * fl1_fx * tly_0_yyy_0[j];

            tlz_y_yyyy_0[j] =
                pa_y[j] * tlz_0_yyyy_0[j] + 2.0 * fl1_fx * tlz_0_yyy_0[j] + 0.5 * fl1_fx * tpx_0_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_0_yyyy_0[j];

            tlx_y_yyyz_0[j] =
                pa_y[j] * tlx_0_yyyz_0[j] + 1.5 * fl1_fx * tlx_0_yyz_0[j] - 0.5 * fl1_fx * tpz_0_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_0_yyyz_0[j];

            tly_y_yyyz_0[j] = pa_y[j] * tly_0_yyyz_0[j] + 1.5 * fl1_fx * tly_0_yyz_0[j];

            tlz_y_yyyz_0[j] =
                pa_y[j] * tlz_0_yyyz_0[j] + 1.5 * fl1_fx * tlz_0_yyz_0[j] + 0.5 * fl1_fx * tpx_0_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_0_yyyz_0[j];

            tlx_y_yyzz_0[j] =
                pa_y[j] * tlx_0_yyzz_0[j] + fl1_fx * tlx_0_yzz_0[j] - 0.5 * fl1_fx * tpz_0_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_0_yyzz_0[j];

            tly_y_yyzz_0[j] = pa_y[j] * tly_0_yyzz_0[j] + fl1_fx * tly_0_yzz_0[j];

            tlz_y_yyzz_0[j] =
                pa_y[j] * tlz_0_yyzz_0[j] + fl1_fx * tlz_0_yzz_0[j] + 0.5 * fl1_fx * tpx_0_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_0_yyzz_0[j];

            tlx_y_yzzz_0[j] =
                pa_y[j] * tlx_0_yzzz_0[j] + 0.5 * fl1_fx * tlx_0_zzz_0[j] - 0.5 * fl1_fx * tpz_0_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_0_yzzz_0[j];

            tly_y_yzzz_0[j] = pa_y[j] * tly_0_yzzz_0[j] + 0.5 * fl1_fx * tly_0_zzz_0[j];

            tlz_y_yzzz_0[j] =
                pa_y[j] * tlz_0_yzzz_0[j] + 0.5 * fl1_fx * tlz_0_zzz_0[j] + 0.5 * fl1_fx * tpx_0_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_0_yzzz_0[j];

            tlx_y_zzzz_0[j] = pa_y[j] * tlx_0_zzzz_0[j] - 0.5 * fl1_fx * tpz_0_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_0_zzzz_0[j];

            tly_y_zzzz_0[j] = pa_y[j] * tly_0_zzzz_0[j];

            tlz_y_zzzz_0[j] = pa_y[j] * tlz_0_zzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPG_90_135(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (90,135)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx);

        auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx);

        auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx);

        auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1);

        auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2);

        auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3);

        auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4);

        auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5);

        auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6);

        auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7);

        auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8);

        auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9);

        auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10);

        auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11);

        auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12);

        auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13);

        auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14);

        auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tdx_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx);

        auto tdy_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx);

        auto tdx_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 1);

        auto tdy_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tdx_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 2);

        auto tdy_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tdx_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 3);

        auto tdy_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tdx_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 4);

        auto tdy_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tdx_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 5);

        auto tdy_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tdx_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 6);

        auto tdy_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tdx_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 7);

        auto tdy_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tdx_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 8);

        auto tdy_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tdx_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 9);

        auto tdy_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tdx_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 10);

        auto tdy_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tdx_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 11);

        auto tdy_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tdx_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 12);

        auto tdy_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tdx_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 13);

        auto tdy_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tdx_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 14);

        auto tdy_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 14);

        // set up pointers to integrals

        auto tlx_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 30);

        auto tly_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tlz_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tlx_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 31);

        auto tly_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tlz_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tlx_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 32);

        auto tly_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tlz_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tlx_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 33);

        auto tly_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tlz_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tlx_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 34);

        auto tly_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tlz_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tlx_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 35);

        auto tly_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tlz_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tlx_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 36);

        auto tly_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tlz_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tlx_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 37);

        auto tly_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tlz_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tlx_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 38);

        auto tly_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tlz_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tlx_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 39);

        auto tly_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tlz_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tlx_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 40);

        auto tly_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tlz_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tlx_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 41);

        auto tly_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tlz_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tlx_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 42);

        auto tly_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tlz_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tlx_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 43);

        auto tly_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tlz_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tlx_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 44);

        auto tly_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tlz_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_z, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, tdx_0_xxyy_0, \
                                     tdx_0_xxyz_0, tdx_0_xxzz_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyzz_0, tdx_0_yzzz_0, tdx_0_zzzz_0, tdy_0_xxxx_0, \
                                     tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxyy_0, tdy_0_xxyz_0, tdy_0_xxzz_0, tdy_0_xyyy_0, \
                                     tdy_0_xyyz_0, tdy_0_xyzz_0, tdy_0_xzzz_0, tdy_0_yyyy_0, tdy_0_yyyz_0, tdy_0_yyzz_0, \
                                     tdy_0_yzzz_0, tdy_0_zzzz_0, tlx_0_xxx_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, \
                                     tlx_0_xxy_0, tlx_0_xxyy_0, tlx_0_xxyz_0, tlx_0_xxz_0, tlx_0_xxzz_0, tlx_0_xyy_0, \
                                     tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyz_0, tlx_0_xyzz_0, tlx_0_xzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyy_0, tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyz_0, tlx_0_yyzz_0, tlx_0_yzz_0, \
                                     tlx_0_yzzz_0, tlx_0_zzz_0, tlx_0_zzzz_0, tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, \
                                     tlx_z_xxyy_0, tlx_z_xxyz_0, tlx_z_xxzz_0, tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyzz_0, \
                                     tlx_z_xzzz_0, tlx_z_yyyy_0, tlx_z_yyyz_0, tlx_z_yyzz_0, tlx_z_yzzz_0, tlx_z_zzzz_0, \
                                     tly_0_xxx_0, tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxy_0, tly_0_xxyy_0, \
                                     tly_0_xxyz_0, tly_0_xxz_0, tly_0_xxzz_0, tly_0_xyy_0, tly_0_xyyy_0, tly_0_xyyz_0, \
                                     tly_0_xyz_0, tly_0_xyzz_0, tly_0_xzz_0, tly_0_xzzz_0, tly_0_yyy_0, tly_0_yyyy_0, \
                                     tly_0_yyyz_0, tly_0_yyz_0, tly_0_yyzz_0, tly_0_yzz_0, tly_0_yzzz_0, tly_0_zzz_0, \
                                     tly_0_zzzz_0, tly_z_xxxx_0, tly_z_xxxy_0, tly_z_xxxz_0, tly_z_xxyy_0, tly_z_xxyz_0, \
                                     tly_z_xxzz_0, tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyzz_0, tly_z_xzzz_0, tly_z_yyyy_0, \
                                     tly_z_yyyz_0, tly_z_yyzz_0, tly_z_yzzz_0, tly_z_zzzz_0, tlz_0_xxx_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxy_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxz_0, \
                                     tlz_0_xxzz_0, tlz_0_xyy_0, tlz_0_xyyy_0, tlz_0_xyyz_0, tlz_0_xyz_0, tlz_0_xyzz_0, \
                                     tlz_0_xzz_0, tlz_0_xzzz_0, tlz_0_yyy_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyz_0, \
                                     tlz_0_yyzz_0, tlz_0_yzz_0, tlz_0_yzzz_0, tlz_0_zzz_0, tlz_0_zzzz_0, tlz_z_xxxx_0, \
                                     tlz_z_xxxy_0, tlz_z_xxxz_0, tlz_z_xxyy_0, tlz_z_xxyz_0, tlz_z_xxzz_0, tlz_z_xyyy_0, \
                                     tlz_z_xyyz_0, tlz_z_xyzz_0, tlz_z_xzzz_0, tlz_z_yyyy_0, tlz_z_yyyz_0, tlz_z_yyzz_0, \
                                     tlz_z_yzzz_0, tlz_z_zzzz_0, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxyy_0, \
                                     tpx_0_xxyz_0, tpx_0_xxzz_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyzz_0, tpx_0_yzzz_0, tpx_0_zzzz_0, tpy_0_xxxx_0, \
                                     tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxyy_0, tpy_0_xxyz_0, tpy_0_xxzz_0, tpy_0_xyyy_0, \
                                     tpy_0_xyyz_0, tpy_0_xyzz_0, tpy_0_xzzz_0, tpy_0_yyyy_0, tpy_0_yyyz_0, tpy_0_yyzz_0, \
                                     tpy_0_yzzz_0, tpy_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_z_xxxx_0[j] = pa_z[j] * tlx_0_xxxx_0[j] + 0.5 * fl1_fx * tpy_0_xxxx_0[j] + fl1_fx * fl1_fgb * tdy_0_xxxx_0[j];

            tly_z_xxxx_0[j] = pa_z[j] * tly_0_xxxx_0[j] - 0.5 * fl1_fx * tpx_0_xxxx_0[j] - fl1_fx * fl1_fgb * tdx_0_xxxx_0[j];

            tlz_z_xxxx_0[j] = pa_z[j] * tlz_0_xxxx_0[j];

            tlx_z_xxxy_0[j] = pa_z[j] * tlx_0_xxxy_0[j] + 0.5 * fl1_fx * tpy_0_xxxy_0[j] + fl1_fx * fl1_fgb * tdy_0_xxxy_0[j];

            tly_z_xxxy_0[j] = pa_z[j] * tly_0_xxxy_0[j] - 0.5 * fl1_fx * tpx_0_xxxy_0[j] - fl1_fx * fl1_fgb * tdx_0_xxxy_0[j];

            tlz_z_xxxy_0[j] = pa_z[j] * tlz_0_xxxy_0[j];

            tlx_z_xxxz_0[j] =
                pa_z[j] * tlx_0_xxxz_0[j] + 0.5 * fl1_fx * tlx_0_xxx_0[j] + 0.5 * fl1_fx * tpy_0_xxxz_0[j] + fl1_fx * fl1_fgb * tdy_0_xxxz_0[j];

            tly_z_xxxz_0[j] =
                pa_z[j] * tly_0_xxxz_0[j] + 0.5 * fl1_fx * tly_0_xxx_0[j] - 0.5 * fl1_fx * tpx_0_xxxz_0[j] - fl1_fx * fl1_fgb * tdx_0_xxxz_0[j];

            tlz_z_xxxz_0[j] = pa_z[j] * tlz_0_xxxz_0[j] + 0.5 * fl1_fx * tlz_0_xxx_0[j];

            tlx_z_xxyy_0[j] = pa_z[j] * tlx_0_xxyy_0[j] + 0.5 * fl1_fx * tpy_0_xxyy_0[j] + fl1_fx * fl1_fgb * tdy_0_xxyy_0[j];

            tly_z_xxyy_0[j] = pa_z[j] * tly_0_xxyy_0[j] - 0.5 * fl1_fx * tpx_0_xxyy_0[j] - fl1_fx * fl1_fgb * tdx_0_xxyy_0[j];

            tlz_z_xxyy_0[j] = pa_z[j] * tlz_0_xxyy_0[j];

            tlx_z_xxyz_0[j] =
                pa_z[j] * tlx_0_xxyz_0[j] + 0.5 * fl1_fx * tlx_0_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xxyz_0[j] + fl1_fx * fl1_fgb * tdy_0_xxyz_0[j];

            tly_z_xxyz_0[j] =
                pa_z[j] * tly_0_xxyz_0[j] + 0.5 * fl1_fx * tly_0_xxy_0[j] - 0.5 * fl1_fx * tpx_0_xxyz_0[j] - fl1_fx * fl1_fgb * tdx_0_xxyz_0[j];

            tlz_z_xxyz_0[j] = pa_z[j] * tlz_0_xxyz_0[j] + 0.5 * fl1_fx * tlz_0_xxy_0[j];

            tlx_z_xxzz_0[j] =
                pa_z[j] * tlx_0_xxzz_0[j] + fl1_fx * tlx_0_xxz_0[j] + 0.5 * fl1_fx * tpy_0_xxzz_0[j] + fl1_fx * fl1_fgb * tdy_0_xxzz_0[j];

            tly_z_xxzz_0[j] =
                pa_z[j] * tly_0_xxzz_0[j] + fl1_fx * tly_0_xxz_0[j] - 0.5 * fl1_fx * tpx_0_xxzz_0[j] - fl1_fx * fl1_fgb * tdx_0_xxzz_0[j];

            tlz_z_xxzz_0[j] = pa_z[j] * tlz_0_xxzz_0[j] + fl1_fx * tlz_0_xxz_0[j];

            tlx_z_xyyy_0[j] = pa_z[j] * tlx_0_xyyy_0[j] + 0.5 * fl1_fx * tpy_0_xyyy_0[j] + fl1_fx * fl1_fgb * tdy_0_xyyy_0[j];

            tly_z_xyyy_0[j] = pa_z[j] * tly_0_xyyy_0[j] - 0.5 * fl1_fx * tpx_0_xyyy_0[j] - fl1_fx * fl1_fgb * tdx_0_xyyy_0[j];

            tlz_z_xyyy_0[j] = pa_z[j] * tlz_0_xyyy_0[j];

            tlx_z_xyyz_0[j] =
                pa_z[j] * tlx_0_xyyz_0[j] + 0.5 * fl1_fx * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tpy_0_xyyz_0[j] + fl1_fx * fl1_fgb * tdy_0_xyyz_0[j];

            tly_z_xyyz_0[j] =
                pa_z[j] * tly_0_xyyz_0[j] + 0.5 * fl1_fx * tly_0_xyy_0[j] - 0.5 * fl1_fx * tpx_0_xyyz_0[j] - fl1_fx * fl1_fgb * tdx_0_xyyz_0[j];

            tlz_z_xyyz_0[j] = pa_z[j] * tlz_0_xyyz_0[j] + 0.5 * fl1_fx * tlz_0_xyy_0[j];

            tlx_z_xyzz_0[j] =
                pa_z[j] * tlx_0_xyzz_0[j] + fl1_fx * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xyzz_0[j] + fl1_fx * fl1_fgb * tdy_0_xyzz_0[j];

            tly_z_xyzz_0[j] =
                pa_z[j] * tly_0_xyzz_0[j] + fl1_fx * tly_0_xyz_0[j] - 0.5 * fl1_fx * tpx_0_xyzz_0[j] - fl1_fx * fl1_fgb * tdx_0_xyzz_0[j];

            tlz_z_xyzz_0[j] = pa_z[j] * tlz_0_xyzz_0[j] + fl1_fx * tlz_0_xyz_0[j];

            tlx_z_xzzz_0[j] =
                pa_z[j] * tlx_0_xzzz_0[j] + 1.5 * fl1_fx * tlx_0_xzz_0[j] + 0.5 * fl1_fx * tpy_0_xzzz_0[j] + fl1_fx * fl1_fgb * tdy_0_xzzz_0[j];

            tly_z_xzzz_0[j] =
                pa_z[j] * tly_0_xzzz_0[j] + 1.5 * fl1_fx * tly_0_xzz_0[j] - 0.5 * fl1_fx * tpx_0_xzzz_0[j] - fl1_fx * fl1_fgb * tdx_0_xzzz_0[j];

            tlz_z_xzzz_0[j] = pa_z[j] * tlz_0_xzzz_0[j] + 1.5 * fl1_fx * tlz_0_xzz_0[j];

            tlx_z_yyyy_0[j] = pa_z[j] * tlx_0_yyyy_0[j] + 0.5 * fl1_fx * tpy_0_yyyy_0[j] + fl1_fx * fl1_fgb * tdy_0_yyyy_0[j];

            tly_z_yyyy_0[j] = pa_z[j] * tly_0_yyyy_0[j] - 0.5 * fl1_fx * tpx_0_yyyy_0[j] - fl1_fx * fl1_fgb * tdx_0_yyyy_0[j];

            tlz_z_yyyy_0[j] = pa_z[j] * tlz_0_yyyy_0[j];

            tlx_z_yyyz_0[j] =
                pa_z[j] * tlx_0_yyyz_0[j] + 0.5 * fl1_fx * tlx_0_yyy_0[j] + 0.5 * fl1_fx * tpy_0_yyyz_0[j] + fl1_fx * fl1_fgb * tdy_0_yyyz_0[j];

            tly_z_yyyz_0[j] =
                pa_z[j] * tly_0_yyyz_0[j] + 0.5 * fl1_fx * tly_0_yyy_0[j] - 0.5 * fl1_fx * tpx_0_yyyz_0[j] - fl1_fx * fl1_fgb * tdx_0_yyyz_0[j];

            tlz_z_yyyz_0[j] = pa_z[j] * tlz_0_yyyz_0[j] + 0.5 * fl1_fx * tlz_0_yyy_0[j];

            tlx_z_yyzz_0[j] =
                pa_z[j] * tlx_0_yyzz_0[j] + fl1_fx * tlx_0_yyz_0[j] + 0.5 * fl1_fx * tpy_0_yyzz_0[j] + fl1_fx * fl1_fgb * tdy_0_yyzz_0[j];

            tly_z_yyzz_0[j] =
                pa_z[j] * tly_0_yyzz_0[j] + fl1_fx * tly_0_yyz_0[j] - 0.5 * fl1_fx * tpx_0_yyzz_0[j] - fl1_fx * fl1_fgb * tdx_0_yyzz_0[j];

            tlz_z_yyzz_0[j] = pa_z[j] * tlz_0_yyzz_0[j] + fl1_fx * tlz_0_yyz_0[j];

            tlx_z_yzzz_0[j] =
                pa_z[j] * tlx_0_yzzz_0[j] + 1.5 * fl1_fx * tlx_0_yzz_0[j] + 0.5 * fl1_fx * tpy_0_yzzz_0[j] + fl1_fx * fl1_fgb * tdy_0_yzzz_0[j];

            tly_z_yzzz_0[j] =
                pa_z[j] * tly_0_yzzz_0[j] + 1.5 * fl1_fx * tly_0_yzz_0[j] - 0.5 * fl1_fx * tpx_0_yzzz_0[j] - fl1_fx * fl1_fgb * tdx_0_yzzz_0[j];

            tlz_z_yzzz_0[j] = pa_z[j] * tlz_0_yzzz_0[j] + 1.5 * fl1_fx * tlz_0_yzz_0[j];

            tlx_z_zzzz_0[j] =
                pa_z[j] * tlx_0_zzzz_0[j] + 2.0 * fl1_fx * tlx_0_zzz_0[j] + 0.5 * fl1_fx * tpy_0_zzzz_0[j] + fl1_fx * fl1_fgb * tdy_0_zzzz_0[j];

            tly_z_zzzz_0[j] =
                pa_z[j] * tly_0_zzzz_0[j] + 2.0 * fl1_fx * tly_0_zzz_0[j] - 0.5 * fl1_fx * tpx_0_zzzz_0[j] - fl1_fx * fl1_fgb * tdx_0_zzzz_0[j];

            tlz_z_zzzz_0[j] = pa_z[j] * tlz_0_zzzz_0[j] + 2.0 * fl1_fx * tlz_0_zzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGP(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForGP_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGP_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGP_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForGP_0_45(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tlx_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx);

        auto tly_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx);

        auto tlz_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx);

        auto tlx_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 1);

        auto tly_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tlz_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tlx_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 2);

        auto tly_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tlz_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tlx_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 3);

        auto tly_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tlz_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tlx_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 4);

        auto tly_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tlz_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tlx_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 5);

        auto tly_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tlz_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tlx_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 6);

        auto tly_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tlz_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tlx_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 7);

        auto tly_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tlz_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tlx_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 8);

        auto tly_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tlz_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tlx_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 9);

        auto tly_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tlz_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tlx_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 10);

        auto tly_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tlz_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tlx_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 11);

        auto tly_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tlz_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tlx_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 12);

        auto tly_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tlz_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tlx_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 13);

        auto tly_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tlz_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tlx_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 14);

        auto tly_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tlz_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 14);

        auto tlx_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx);

        auto tly_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx);

        auto tlz_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx);

        auto tlx_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 1);

        auto tly_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tlz_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tlx_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 2);

        auto tly_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tlz_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tlx_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 3);

        auto tly_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tlz_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tlx_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 4);

        auto tly_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tlz_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tlx_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 5);

        auto tly_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tlz_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tlx_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 6);

        auto tly_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tlz_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tlx_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 7);

        auto tly_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tlz_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tlx_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 8);

        auto tly_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tlz_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9);

        auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10);

        auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11);

        auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12);

        auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13);

        auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14);

        auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx);

        auto tly_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx);

        auto tlz_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx);

        auto tlx_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 1);

        auto tly_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 2);

        auto tly_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 3);

        auto tly_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 4);

        auto tly_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tpy_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx);

        auto tpz_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx);

        auto tpy_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tpy_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tpy_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tpy_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tpy_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tpy_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tpy_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tpy_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tpy_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tpy_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tpy_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tpy_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tpy_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tpy_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 14);

        auto tdy_xxx_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx);

        auto tdz_xxx_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx);

        auto tdy_xxx_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tdz_xxx_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tdy_xxx_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tdz_xxx_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tdy_xxy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tdz_xxy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tdy_xxy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tdz_xxy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tdy_xxy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tdz_xxy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tdy_xxz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tdz_xxz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tdy_xxz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tdz_xxz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tdy_xxz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tdz_xxz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tdy_xyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tdz_xyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tdy_xyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tdz_xyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tdy_xyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tdz_xyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tdy_xyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tdz_xyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tdy_xyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tdz_xyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tdy_xyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tdz_xyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 14);

        // set up pointers to integrals

        auto tlx_xxxx_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx);

        auto tly_xxxx_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx);

        auto tlz_xxxx_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx);

        auto tlx_xxxx_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 1);

        auto tly_xxxx_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 1);

        auto tlz_xxxx_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 1);

        auto tlx_xxxx_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 2);

        auto tly_xxxx_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 2);

        auto tlz_xxxx_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 2);

        auto tlx_xxxy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 3);

        auto tly_xxxy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 3);

        auto tlz_xxxy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 3);

        auto tlx_xxxy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 4);

        auto tly_xxxy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 4);

        auto tlz_xxxy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 4);

        auto tlx_xxxy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 5);

        auto tly_xxxy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 5);

        auto tlz_xxxy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 5);

        auto tlx_xxxz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 6);

        auto tly_xxxz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 6);

        auto tlz_xxxz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 6);

        auto tlx_xxxz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 7);

        auto tly_xxxz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 7);

        auto tlz_xxxz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 7);

        auto tlx_xxxz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 8);

        auto tly_xxxz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 8);

        auto tlz_xxxz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 8);

        auto tlx_xxyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 9);

        auto tly_xxyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 9);

        auto tlz_xxyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 9);

        auto tlx_xxyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 10);

        auto tly_xxyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 10);

        auto tlz_xxyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 10);

        auto tlx_xxyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 11);

        auto tly_xxyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 11);

        auto tlz_xxyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 11);

        auto tlx_xxyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 12);

        auto tly_xxyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 12);

        auto tlz_xxyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 12);

        auto tlx_xxyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 13);

        auto tly_xxyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 13);

        auto tlz_xxyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 13);

        auto tlx_xxyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 14);

        auto tly_xxyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 14);

        auto tlz_xxyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxx_x_0, tdy_xxx_y_0, tdy_xxx_z_0, tdy_xxy_x_0, \
                                     tdy_xxy_y_0, tdy_xxy_z_0, tdy_xxz_x_0, tdy_xxz_y_0, tdy_xxz_z_0, tdy_xyy_x_0, \
                                     tdy_xyy_y_0, tdy_xyy_z_0, tdy_xyz_x_0, tdy_xyz_y_0, tdy_xyz_z_0, tdz_xxx_x_0, \
                                     tdz_xxx_y_0, tdz_xxx_z_0, tdz_xxy_x_0, tdz_xxy_y_0, tdz_xxy_z_0, tdz_xxz_x_0, \
                                     tdz_xxz_y_0, tdz_xxz_z_0, tdz_xyy_x_0, tdz_xyy_y_0, tdz_xyy_z_0, tdz_xyz_x_0, \
                                     tdz_xyz_y_0, tdz_xyz_z_0, tlx_xx_x_0, tlx_xx_y_0, tlx_xx_z_0, tlx_xxx_0_0, \
                                     tlx_xxx_x_0, tlx_xxx_y_0, tlx_xxx_z_0, tlx_xxxx_x_0, tlx_xxxx_y_0, tlx_xxxx_z_0, \
                                     tlx_xxxy_x_0, tlx_xxxy_y_0, tlx_xxxy_z_0, tlx_xxxz_x_0, tlx_xxxz_y_0, tlx_xxxz_z_0, \
                                     tlx_xxy_0_0, tlx_xxy_x_0, tlx_xxy_y_0, tlx_xxy_z_0, tlx_xxyy_x_0, tlx_xxyy_y_0, \
                                     tlx_xxyy_z_0, tlx_xxyz_x_0, tlx_xxyz_y_0, tlx_xxyz_z_0, tlx_xxz_0_0, tlx_xxz_x_0, \
                                     tlx_xxz_y_0, tlx_xxz_z_0, tlx_xy_x_0, tlx_xy_y_0, tlx_xy_z_0, tlx_xyy_0_0, \
                                     tlx_xyy_x_0, tlx_xyy_y_0, tlx_xyy_z_0, tlx_xyz_0_0, tlx_xyz_x_0, tlx_xyz_y_0, \
                                     tlx_xyz_z_0, tlx_xz_x_0, tlx_xz_y_0, tlx_xz_z_0, tlx_yy_x_0, tlx_yy_y_0, tlx_yy_z_0, \
                                     tlx_yz_x_0, tlx_yz_y_0, tlx_yz_z_0, tly_xx_x_0, tly_xx_y_0, tly_xx_z_0, \
                                     tly_xxx_0_0, tly_xxx_x_0, tly_xxx_y_0, tly_xxx_z_0, tly_xxxx_x_0, tly_xxxx_y_0, \
                                     tly_xxxx_z_0, tly_xxxy_x_0, tly_xxxy_y_0, tly_xxxy_z_0, tly_xxxz_x_0, tly_xxxz_y_0, \
                                     tly_xxxz_z_0, tly_xxy_0_0, tly_xxy_x_0, tly_xxy_y_0, tly_xxy_z_0, tly_xxyy_x_0, \
                                     tly_xxyy_y_0, tly_xxyy_z_0, tly_xxyz_x_0, tly_xxyz_y_0, tly_xxyz_z_0, tly_xxz_0_0, \
                                     tly_xxz_x_0, tly_xxz_y_0, tly_xxz_z_0, tly_xy_x_0, tly_xy_y_0, tly_xy_z_0, \
                                     tly_xyy_0_0, tly_xyy_x_0, tly_xyy_y_0, tly_xyy_z_0, tly_xyz_0_0, tly_xyz_x_0, \
                                     tly_xyz_y_0, tly_xyz_z_0, tly_xz_x_0, tly_xz_y_0, tly_xz_z_0, tly_yy_x_0, \
                                     tly_yy_y_0, tly_yy_z_0, tly_yz_x_0, tly_yz_y_0, tly_yz_z_0, tlz_xx_x_0, tlz_xx_y_0, \
                                     tlz_xx_z_0, tlz_xxx_0_0, tlz_xxx_x_0, tlz_xxx_y_0, tlz_xxx_z_0, tlz_xxxx_x_0, \
                                     tlz_xxxx_y_0, tlz_xxxx_z_0, tlz_xxxy_x_0, tlz_xxxy_y_0, tlz_xxxy_z_0, tlz_xxxz_x_0, \
                                     tlz_xxxz_y_0, tlz_xxxz_z_0, tlz_xxy_0_0, tlz_xxy_x_0, tlz_xxy_y_0, tlz_xxy_z_0, \
                                     tlz_xxyy_x_0, tlz_xxyy_y_0, tlz_xxyy_z_0, tlz_xxyz_x_0, tlz_xxyz_y_0, tlz_xxyz_z_0, \
                                     tlz_xxz_0_0, tlz_xxz_x_0, tlz_xxz_y_0, tlz_xxz_z_0, tlz_xy_x_0, tlz_xy_y_0, \
                                     tlz_xy_z_0, tlz_xyy_0_0, tlz_xyy_x_0, tlz_xyy_y_0, tlz_xyy_z_0, tlz_xyz_0_0, \
                                     tlz_xyz_x_0, tlz_xyz_y_0, tlz_xyz_z_0, tlz_xz_x_0, tlz_xz_y_0, tlz_xz_z_0, \
                                     tlz_yy_x_0, tlz_yy_y_0, tlz_yy_z_0, tlz_yz_x_0, tlz_yz_y_0, tlz_yz_z_0, \
                                     tpy_xxx_x_0, tpy_xxx_y_0, tpy_xxx_z_0, tpy_xxy_x_0, tpy_xxy_y_0, tpy_xxy_z_0, \
                                     tpy_xxz_x_0, tpy_xxz_y_0, tpy_xxz_z_0, tpy_xyy_x_0, tpy_xyy_y_0, tpy_xyy_z_0, \
                                     tpy_xyz_x_0, tpy_xyz_y_0, tpy_xyz_z_0, tpz_xxx_x_0, tpz_xxx_y_0, tpz_xxx_z_0, \
                                     tpz_xxy_x_0, tpz_xxy_y_0, tpz_xxy_z_0, tpz_xxz_x_0, tpz_xxz_y_0, tpz_xxz_z_0, \
                                     tpz_xyy_x_0, tpz_xyy_y_0, tpz_xyy_z_0, tpz_xyz_x_0, tpz_xyz_y_0, tpz_xyz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxxx_x_0[j] = pa_x[j] * tlx_xxx_x_0[j] + 1.5 * fl1_fx * tlx_xx_x_0[j] + 0.5 * fl1_fx * tlx_xxx_0_0[j];

            tly_xxxx_x_0[j] = pa_x[j] * tly_xxx_x_0[j] + 1.5 * fl1_fx * tly_xx_x_0[j] + 0.5 * fl1_fx * tly_xxx_0_0[j] +
                              0.5 * fl1_fx * tpz_xxx_x_0[j] + fl1_fx * fl1_fgb * tdz_xxx_x_0[j];

            tlz_xxxx_x_0[j] = pa_x[j] * tlz_xxx_x_0[j] + 1.5 * fl1_fx * tlz_xx_x_0[j] + 0.5 * fl1_fx * tlz_xxx_0_0[j] -
                              0.5 * fl1_fx * tpy_xxx_x_0[j] - fl1_fx * fl1_fgb * tdy_xxx_x_0[j];

            tlx_xxxx_y_0[j] = pa_x[j] * tlx_xxx_y_0[j] + 1.5 * fl1_fx * tlx_xx_y_0[j];

            tly_xxxx_y_0[j] =
                pa_x[j] * tly_xxx_y_0[j] + 1.5 * fl1_fx * tly_xx_y_0[j] + 0.5 * fl1_fx * tpz_xxx_y_0[j] + fl1_fx * fl1_fgb * tdz_xxx_y_0[j];

            tlz_xxxx_y_0[j] =
                pa_x[j] * tlz_xxx_y_0[j] + 1.5 * fl1_fx * tlz_xx_y_0[j] - 0.5 * fl1_fx * tpy_xxx_y_0[j] - fl1_fx * fl1_fgb * tdy_xxx_y_0[j];

            tlx_xxxx_z_0[j] = pa_x[j] * tlx_xxx_z_0[j] + 1.5 * fl1_fx * tlx_xx_z_0[j];

            tly_xxxx_z_0[j] =
                pa_x[j] * tly_xxx_z_0[j] + 1.5 * fl1_fx * tly_xx_z_0[j] + 0.5 * fl1_fx * tpz_xxx_z_0[j] + fl1_fx * fl1_fgb * tdz_xxx_z_0[j];

            tlz_xxxx_z_0[j] =
                pa_x[j] * tlz_xxx_z_0[j] + 1.5 * fl1_fx * tlz_xx_z_0[j] - 0.5 * fl1_fx * tpy_xxx_z_0[j] - fl1_fx * fl1_fgb * tdy_xxx_z_0[j];

            tlx_xxxy_x_0[j] = pa_x[j] * tlx_xxy_x_0[j] + fl1_fx * tlx_xy_x_0[j] + 0.5 * fl1_fx * tlx_xxy_0_0[j];

            tly_xxxy_x_0[j] = pa_x[j] * tly_xxy_x_0[j] + fl1_fx * tly_xy_x_0[j] + 0.5 * fl1_fx * tly_xxy_0_0[j] + 0.5 * fl1_fx * tpz_xxy_x_0[j] +
                              fl1_fx * fl1_fgb * tdz_xxy_x_0[j];

            tlz_xxxy_x_0[j] = pa_x[j] * tlz_xxy_x_0[j] + fl1_fx * tlz_xy_x_0[j] + 0.5 * fl1_fx * tlz_xxy_0_0[j] - 0.5 * fl1_fx * tpy_xxy_x_0[j] -
                              fl1_fx * fl1_fgb * tdy_xxy_x_0[j];

            tlx_xxxy_y_0[j] = pa_x[j] * tlx_xxy_y_0[j] + fl1_fx * tlx_xy_y_0[j];

            tly_xxxy_y_0[j] = pa_x[j] * tly_xxy_y_0[j] + fl1_fx * tly_xy_y_0[j] + 0.5 * fl1_fx * tpz_xxy_y_0[j] + fl1_fx * fl1_fgb * tdz_xxy_y_0[j];

            tlz_xxxy_y_0[j] = pa_x[j] * tlz_xxy_y_0[j] + fl1_fx * tlz_xy_y_0[j] - 0.5 * fl1_fx * tpy_xxy_y_0[j] - fl1_fx * fl1_fgb * tdy_xxy_y_0[j];

            tlx_xxxy_z_0[j] = pa_x[j] * tlx_xxy_z_0[j] + fl1_fx * tlx_xy_z_0[j];

            tly_xxxy_z_0[j] = pa_x[j] * tly_xxy_z_0[j] + fl1_fx * tly_xy_z_0[j] + 0.5 * fl1_fx * tpz_xxy_z_0[j] + fl1_fx * fl1_fgb * tdz_xxy_z_0[j];

            tlz_xxxy_z_0[j] = pa_x[j] * tlz_xxy_z_0[j] + fl1_fx * tlz_xy_z_0[j] - 0.5 * fl1_fx * tpy_xxy_z_0[j] - fl1_fx * fl1_fgb * tdy_xxy_z_0[j];

            tlx_xxxz_x_0[j] = pa_x[j] * tlx_xxz_x_0[j] + fl1_fx * tlx_xz_x_0[j] + 0.5 * fl1_fx * tlx_xxz_0_0[j];

            tly_xxxz_x_0[j] = pa_x[j] * tly_xxz_x_0[j] + fl1_fx * tly_xz_x_0[j] + 0.5 * fl1_fx * tly_xxz_0_0[j] + 0.5 * fl1_fx * tpz_xxz_x_0[j] +
                              fl1_fx * fl1_fgb * tdz_xxz_x_0[j];

            tlz_xxxz_x_0[j] = pa_x[j] * tlz_xxz_x_0[j] + fl1_fx * tlz_xz_x_0[j] + 0.5 * fl1_fx * tlz_xxz_0_0[j] - 0.5 * fl1_fx * tpy_xxz_x_0[j] -
                              fl1_fx * fl1_fgb * tdy_xxz_x_0[j];

            tlx_xxxz_y_0[j] = pa_x[j] * tlx_xxz_y_0[j] + fl1_fx * tlx_xz_y_0[j];

            tly_xxxz_y_0[j] = pa_x[j] * tly_xxz_y_0[j] + fl1_fx * tly_xz_y_0[j] + 0.5 * fl1_fx * tpz_xxz_y_0[j] + fl1_fx * fl1_fgb * tdz_xxz_y_0[j];

            tlz_xxxz_y_0[j] = pa_x[j] * tlz_xxz_y_0[j] + fl1_fx * tlz_xz_y_0[j] - 0.5 * fl1_fx * tpy_xxz_y_0[j] - fl1_fx * fl1_fgb * tdy_xxz_y_0[j];

            tlx_xxxz_z_0[j] = pa_x[j] * tlx_xxz_z_0[j] + fl1_fx * tlx_xz_z_0[j];

            tly_xxxz_z_0[j] = pa_x[j] * tly_xxz_z_0[j] + fl1_fx * tly_xz_z_0[j] + 0.5 * fl1_fx * tpz_xxz_z_0[j] + fl1_fx * fl1_fgb * tdz_xxz_z_0[j];

            tlz_xxxz_z_0[j] = pa_x[j] * tlz_xxz_z_0[j] + fl1_fx * tlz_xz_z_0[j] - 0.5 * fl1_fx * tpy_xxz_z_0[j] - fl1_fx * fl1_fgb * tdy_xxz_z_0[j];

            tlx_xxyy_x_0[j] = pa_x[j] * tlx_xyy_x_0[j] + 0.5 * fl1_fx * tlx_yy_x_0[j] + 0.5 * fl1_fx * tlx_xyy_0_0[j];

            tly_xxyy_x_0[j] = pa_x[j] * tly_xyy_x_0[j] + 0.5 * fl1_fx * tly_yy_x_0[j] + 0.5 * fl1_fx * tly_xyy_0_0[j] +
                              0.5 * fl1_fx * tpz_xyy_x_0[j] + fl1_fx * fl1_fgb * tdz_xyy_x_0[j];

            tlz_xxyy_x_0[j] = pa_x[j] * tlz_xyy_x_0[j] + 0.5 * fl1_fx * tlz_yy_x_0[j] + 0.5 * fl1_fx * tlz_xyy_0_0[j] -
                              0.5 * fl1_fx * tpy_xyy_x_0[j] - fl1_fx * fl1_fgb * tdy_xyy_x_0[j];

            tlx_xxyy_y_0[j] = pa_x[j] * tlx_xyy_y_0[j] + 0.5 * fl1_fx * tlx_yy_y_0[j];

            tly_xxyy_y_0[j] =
                pa_x[j] * tly_xyy_y_0[j] + 0.5 * fl1_fx * tly_yy_y_0[j] + 0.5 * fl1_fx * tpz_xyy_y_0[j] + fl1_fx * fl1_fgb * tdz_xyy_y_0[j];

            tlz_xxyy_y_0[j] =
                pa_x[j] * tlz_xyy_y_0[j] + 0.5 * fl1_fx * tlz_yy_y_0[j] - 0.5 * fl1_fx * tpy_xyy_y_0[j] - fl1_fx * fl1_fgb * tdy_xyy_y_0[j];

            tlx_xxyy_z_0[j] = pa_x[j] * tlx_xyy_z_0[j] + 0.5 * fl1_fx * tlx_yy_z_0[j];

            tly_xxyy_z_0[j] =
                pa_x[j] * tly_xyy_z_0[j] + 0.5 * fl1_fx * tly_yy_z_0[j] + 0.5 * fl1_fx * tpz_xyy_z_0[j] + fl1_fx * fl1_fgb * tdz_xyy_z_0[j];

            tlz_xxyy_z_0[j] =
                pa_x[j] * tlz_xyy_z_0[j] + 0.5 * fl1_fx * tlz_yy_z_0[j] - 0.5 * fl1_fx * tpy_xyy_z_0[j] - fl1_fx * fl1_fgb * tdy_xyy_z_0[j];

            tlx_xxyz_x_0[j] = pa_x[j] * tlx_xyz_x_0[j] + 0.5 * fl1_fx * tlx_yz_x_0[j] + 0.5 * fl1_fx * tlx_xyz_0_0[j];

            tly_xxyz_x_0[j] = pa_x[j] * tly_xyz_x_0[j] + 0.5 * fl1_fx * tly_yz_x_0[j] + 0.5 * fl1_fx * tly_xyz_0_0[j] +
                              0.5 * fl1_fx * tpz_xyz_x_0[j] + fl1_fx * fl1_fgb * tdz_xyz_x_0[j];

            tlz_xxyz_x_0[j] = pa_x[j] * tlz_xyz_x_0[j] + 0.5 * fl1_fx * tlz_yz_x_0[j] + 0.5 * fl1_fx * tlz_xyz_0_0[j] -
                              0.5 * fl1_fx * tpy_xyz_x_0[j] - fl1_fx * fl1_fgb * tdy_xyz_x_0[j];

            tlx_xxyz_y_0[j] = pa_x[j] * tlx_xyz_y_0[j] + 0.5 * fl1_fx * tlx_yz_y_0[j];

            tly_xxyz_y_0[j] =
                pa_x[j] * tly_xyz_y_0[j] + 0.5 * fl1_fx * tly_yz_y_0[j] + 0.5 * fl1_fx * tpz_xyz_y_0[j] + fl1_fx * fl1_fgb * tdz_xyz_y_0[j];

            tlz_xxyz_y_0[j] =
                pa_x[j] * tlz_xyz_y_0[j] + 0.5 * fl1_fx * tlz_yz_y_0[j] - 0.5 * fl1_fx * tpy_xyz_y_0[j] - fl1_fx * fl1_fgb * tdy_xyz_y_0[j];

            tlx_xxyz_z_0[j] = pa_x[j] * tlx_xyz_z_0[j] + 0.5 * fl1_fx * tlx_yz_z_0[j];

            tly_xxyz_z_0[j] =
                pa_x[j] * tly_xyz_z_0[j] + 0.5 * fl1_fx * tly_yz_z_0[j] + 0.5 * fl1_fx * tpz_xyz_z_0[j] + fl1_fx * fl1_fgb * tdz_xyz_z_0[j];

            tlz_xxyz_z_0[j] =
                pa_x[j] * tlz_xyz_z_0[j] + 0.5 * fl1_fx * tlz_yz_z_0[j] - 0.5 * fl1_fx * tpy_xyz_z_0[j] - fl1_fx * fl1_fgb * tdy_xyz_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGP_45_90(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tlx_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 15);

        auto tly_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tlz_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tlx_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 16);

        auto tly_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tlz_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tlx_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 17);

        auto tly_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tlz_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tlx_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 18);

        auto tly_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tlz_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tlx_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 19);

        auto tly_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tlz_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tlx_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 20);

        auto tly_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tlz_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tlx_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 21);

        auto tly_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tlz_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tlx_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 22);

        auto tly_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tlz_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tlx_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 23);

        auto tly_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tlz_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tlx_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 24);

        auto tly_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tlz_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tlx_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 25);

        auto tly_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tlz_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tlx_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 26);

        auto tly_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tlz_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tlx_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 27);

        auto tly_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tlz_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tlx_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 28);

        auto tly_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tlz_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tlx_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 29);

        auto tly_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tlz_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15);

        auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16);

        auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17);

        auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tlx_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 5);

        auto tly_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 6);

        auto tly_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 7);

        auto tly_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 8);

        auto tly_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 9);

        auto tly_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 9);

        auto tpy_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tpy_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tpy_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tpy_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tpy_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tpy_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tdy_xzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tdz_xzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tdy_xzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tdz_xzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tdy_xzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tdz_xzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tdy_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tdz_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tdy_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tdz_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tdy_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tdz_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tdy_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tdz_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tdy_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tdz_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tdy_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tdz_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tdy_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tdz_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tdy_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tdz_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tdy_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tdz_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tdy_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tdz_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tdy_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tdz_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tdy_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tdz_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 29);

        // set up pointers to integrals

        auto tlx_xxzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 15);

        auto tly_xxzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 15);

        auto tlz_xxzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 15);

        auto tlx_xxzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 16);

        auto tly_xxzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 16);

        auto tlz_xxzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 16);

        auto tlx_xxzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 17);

        auto tly_xxzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 17);

        auto tlz_xxzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 17);

        auto tlx_xyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 18);

        auto tly_xyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 18);

        auto tlz_xyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 18);

        auto tlx_xyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 19);

        auto tly_xyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 19);

        auto tlz_xyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 19);

        auto tlx_xyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 20);

        auto tly_xyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 20);

        auto tlz_xyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 20);

        auto tlx_xyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 21);

        auto tly_xyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 21);

        auto tlz_xyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 21);

        auto tlx_xyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 22);

        auto tly_xyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 22);

        auto tlz_xyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 22);

        auto tlx_xyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 23);

        auto tly_xyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 23);

        auto tlz_xyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 23);

        auto tlx_xyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 24);

        auto tly_xyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 24);

        auto tlz_xyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 24);

        auto tlx_xyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 25);

        auto tly_xyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 25);

        auto tlz_xyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 25);

        auto tlx_xyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 26);

        auto tly_xyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 26);

        auto tlz_xyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 26);

        auto tlx_xzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 27);

        auto tly_xzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 27);

        auto tlz_xzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 27);

        auto tlx_xzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 28);

        auto tly_xzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 28);

        auto tlz_xzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 28);

        auto tlx_xzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 29);

        auto tly_xzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 29);

        auto tlz_xzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xzz_x_0, tdy_xzz_y_0, tdy_xzz_z_0, tdy_yyy_x_0, \
                                     tdy_yyy_y_0, tdy_yyy_z_0, tdy_yyz_x_0, tdy_yyz_y_0, tdy_yyz_z_0, tdy_yzz_x_0, \
                                     tdy_yzz_y_0, tdy_yzz_z_0, tdy_zzz_x_0, tdy_zzz_y_0, tdy_zzz_z_0, tdz_xzz_x_0, \
                                     tdz_xzz_y_0, tdz_xzz_z_0, tdz_yyy_x_0, tdz_yyy_y_0, tdz_yyy_z_0, tdz_yyz_x_0, \
                                     tdz_yyz_y_0, tdz_yyz_z_0, tdz_yzz_x_0, tdz_yzz_y_0, tdz_yzz_z_0, tdz_zzz_x_0, \
                                     tdz_zzz_y_0, tdz_zzz_z_0, tlx_xxzz_x_0, tlx_xxzz_y_0, tlx_xxzz_z_0, tlx_xyyy_x_0, \
                                     tlx_xyyy_y_0, tlx_xyyy_z_0, tlx_xyyz_x_0, tlx_xyyz_y_0, tlx_xyyz_z_0, tlx_xyzz_x_0, \
                                     tlx_xyzz_y_0, tlx_xyzz_z_0, tlx_xzz_0_0, tlx_xzz_x_0, tlx_xzz_y_0, tlx_xzz_z_0, \
                                     tlx_xzzz_x_0, tlx_xzzz_y_0, tlx_xzzz_z_0, tlx_yyy_0_0, tlx_yyy_x_0, tlx_yyy_y_0, \
                                     tlx_yyy_z_0, tlx_yyz_0_0, tlx_yyz_x_0, tlx_yyz_y_0, tlx_yyz_z_0, tlx_yzz_0_0, \
                                     tlx_yzz_x_0, tlx_yzz_y_0, tlx_yzz_z_0, tlx_zz_x_0, tlx_zz_y_0, tlx_zz_z_0, \
                                     tlx_zzz_0_0, tlx_zzz_x_0, tlx_zzz_y_0, tlx_zzz_z_0, tly_xxzz_x_0, tly_xxzz_y_0, \
                                     tly_xxzz_z_0, tly_xyyy_x_0, tly_xyyy_y_0, tly_xyyy_z_0, tly_xyyz_x_0, tly_xyyz_y_0, \
                                     tly_xyyz_z_0, tly_xyzz_x_0, tly_xyzz_y_0, tly_xyzz_z_0, tly_xzz_0_0, tly_xzz_x_0, \
                                     tly_xzz_y_0, tly_xzz_z_0, tly_xzzz_x_0, tly_xzzz_y_0, tly_xzzz_z_0, tly_yyy_0_0, \
                                     tly_yyy_x_0, tly_yyy_y_0, tly_yyy_z_0, tly_yyz_0_0, tly_yyz_x_0, tly_yyz_y_0, \
                                     tly_yyz_z_0, tly_yzz_0_0, tly_yzz_x_0, tly_yzz_y_0, tly_yzz_z_0, tly_zz_x_0, \
                                     tly_zz_y_0, tly_zz_z_0, tly_zzz_0_0, tly_zzz_x_0, tly_zzz_y_0, tly_zzz_z_0, \
                                     tlz_xxzz_x_0, tlz_xxzz_y_0, tlz_xxzz_z_0, tlz_xyyy_x_0, tlz_xyyy_y_0, tlz_xyyy_z_0, \
                                     tlz_xyyz_x_0, tlz_xyyz_y_0, tlz_xyyz_z_0, tlz_xyzz_x_0, tlz_xyzz_y_0, tlz_xyzz_z_0, \
                                     tlz_xzz_0_0, tlz_xzz_x_0, tlz_xzz_y_0, tlz_xzz_z_0, tlz_xzzz_x_0, tlz_xzzz_y_0, \
                                     tlz_xzzz_z_0, tlz_yyy_0_0, tlz_yyy_x_0, tlz_yyy_y_0, tlz_yyy_z_0, tlz_yyz_0_0, \
                                     tlz_yyz_x_0, tlz_yyz_y_0, tlz_yyz_z_0, tlz_yzz_0_0, tlz_yzz_x_0, tlz_yzz_y_0, \
                                     tlz_yzz_z_0, tlz_zz_x_0, tlz_zz_y_0, tlz_zz_z_0, tlz_zzz_0_0, tlz_zzz_x_0, \
                                     tlz_zzz_y_0, tlz_zzz_z_0, tpy_xzz_x_0, tpy_xzz_y_0, tpy_xzz_z_0, tpy_yyy_x_0, \
                                     tpy_yyy_y_0, tpy_yyy_z_0, tpy_yyz_x_0, tpy_yyz_y_0, tpy_yyz_z_0, tpy_yzz_x_0, \
                                     tpy_yzz_y_0, tpy_yzz_z_0, tpy_zzz_x_0, tpy_zzz_y_0, tpy_zzz_z_0, tpz_xzz_x_0, \
                                     tpz_xzz_y_0, tpz_xzz_z_0, tpz_yyy_x_0, tpz_yyy_y_0, tpz_yyy_z_0, tpz_yyz_x_0, \
                                     tpz_yyz_y_0, tpz_yyz_z_0, tpz_yzz_x_0, tpz_yzz_y_0, tpz_yzz_z_0, tpz_zzz_x_0, \
                                     tpz_zzz_y_0, tpz_zzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxzz_x_0[j] = pa_x[j] * tlx_xzz_x_0[j] + 0.5 * fl1_fx * tlx_zz_x_0[j] + 0.5 * fl1_fx * tlx_xzz_0_0[j];

            tly_xxzz_x_0[j] = pa_x[j] * tly_xzz_x_0[j] + 0.5 * fl1_fx * tly_zz_x_0[j] + 0.5 * fl1_fx * tly_xzz_0_0[j] +
                              0.5 * fl1_fx * tpz_xzz_x_0[j] + fl1_fx * fl1_fgb * tdz_xzz_x_0[j];

            tlz_xxzz_x_0[j] = pa_x[j] * tlz_xzz_x_0[j] + 0.5 * fl1_fx * tlz_zz_x_0[j] + 0.5 * fl1_fx * tlz_xzz_0_0[j] -
                              0.5 * fl1_fx * tpy_xzz_x_0[j] - fl1_fx * fl1_fgb * tdy_xzz_x_0[j];

            tlx_xxzz_y_0[j] = pa_x[j] * tlx_xzz_y_0[j] + 0.5 * fl1_fx * tlx_zz_y_0[j];

            tly_xxzz_y_0[j] =
                pa_x[j] * tly_xzz_y_0[j] + 0.5 * fl1_fx * tly_zz_y_0[j] + 0.5 * fl1_fx * tpz_xzz_y_0[j] + fl1_fx * fl1_fgb * tdz_xzz_y_0[j];

            tlz_xxzz_y_0[j] =
                pa_x[j] * tlz_xzz_y_0[j] + 0.5 * fl1_fx * tlz_zz_y_0[j] - 0.5 * fl1_fx * tpy_xzz_y_0[j] - fl1_fx * fl1_fgb * tdy_xzz_y_0[j];

            tlx_xxzz_z_0[j] = pa_x[j] * tlx_xzz_z_0[j] + 0.5 * fl1_fx * tlx_zz_z_0[j];

            tly_xxzz_z_0[j] =
                pa_x[j] * tly_xzz_z_0[j] + 0.5 * fl1_fx * tly_zz_z_0[j] + 0.5 * fl1_fx * tpz_xzz_z_0[j] + fl1_fx * fl1_fgb * tdz_xzz_z_0[j];

            tlz_xxzz_z_0[j] =
                pa_x[j] * tlz_xzz_z_0[j] + 0.5 * fl1_fx * tlz_zz_z_0[j] - 0.5 * fl1_fx * tpy_xzz_z_0[j] - fl1_fx * fl1_fgb * tdy_xzz_z_0[j];

            tlx_xyyy_x_0[j] = pa_x[j] * tlx_yyy_x_0[j] + 0.5 * fl1_fx * tlx_yyy_0_0[j];

            tly_xyyy_x_0[j] =
                pa_x[j] * tly_yyy_x_0[j] + 0.5 * fl1_fx * tly_yyy_0_0[j] + 0.5 * fl1_fx * tpz_yyy_x_0[j] + fl1_fx * fl1_fgb * tdz_yyy_x_0[j];

            tlz_xyyy_x_0[j] =
                pa_x[j] * tlz_yyy_x_0[j] + 0.5 * fl1_fx * tlz_yyy_0_0[j] - 0.5 * fl1_fx * tpy_yyy_x_0[j] - fl1_fx * fl1_fgb * tdy_yyy_x_0[j];

            tlx_xyyy_y_0[j] = pa_x[j] * tlx_yyy_y_0[j];

            tly_xyyy_y_0[j] = pa_x[j] * tly_yyy_y_0[j] + 0.5 * fl1_fx * tpz_yyy_y_0[j] + fl1_fx * fl1_fgb * tdz_yyy_y_0[j];

            tlz_xyyy_y_0[j] = pa_x[j] * tlz_yyy_y_0[j] - 0.5 * fl1_fx * tpy_yyy_y_0[j] - fl1_fx * fl1_fgb * tdy_yyy_y_0[j];

            tlx_xyyy_z_0[j] = pa_x[j] * tlx_yyy_z_0[j];

            tly_xyyy_z_0[j] = pa_x[j] * tly_yyy_z_0[j] + 0.5 * fl1_fx * tpz_yyy_z_0[j] + fl1_fx * fl1_fgb * tdz_yyy_z_0[j];

            tlz_xyyy_z_0[j] = pa_x[j] * tlz_yyy_z_0[j] - 0.5 * fl1_fx * tpy_yyy_z_0[j] - fl1_fx * fl1_fgb * tdy_yyy_z_0[j];

            tlx_xyyz_x_0[j] = pa_x[j] * tlx_yyz_x_0[j] + 0.5 * fl1_fx * tlx_yyz_0_0[j];

            tly_xyyz_x_0[j] =
                pa_x[j] * tly_yyz_x_0[j] + 0.5 * fl1_fx * tly_yyz_0_0[j] + 0.5 * fl1_fx * tpz_yyz_x_0[j] + fl1_fx * fl1_fgb * tdz_yyz_x_0[j];

            tlz_xyyz_x_0[j] =
                pa_x[j] * tlz_yyz_x_0[j] + 0.5 * fl1_fx * tlz_yyz_0_0[j] - 0.5 * fl1_fx * tpy_yyz_x_0[j] - fl1_fx * fl1_fgb * tdy_yyz_x_0[j];

            tlx_xyyz_y_0[j] = pa_x[j] * tlx_yyz_y_0[j];

            tly_xyyz_y_0[j] = pa_x[j] * tly_yyz_y_0[j] + 0.5 * fl1_fx * tpz_yyz_y_0[j] + fl1_fx * fl1_fgb * tdz_yyz_y_0[j];

            tlz_xyyz_y_0[j] = pa_x[j] * tlz_yyz_y_0[j] - 0.5 * fl1_fx * tpy_yyz_y_0[j] - fl1_fx * fl1_fgb * tdy_yyz_y_0[j];

            tlx_xyyz_z_0[j] = pa_x[j] * tlx_yyz_z_0[j];

            tly_xyyz_z_0[j] = pa_x[j] * tly_yyz_z_0[j] + 0.5 * fl1_fx * tpz_yyz_z_0[j] + fl1_fx * fl1_fgb * tdz_yyz_z_0[j];

            tlz_xyyz_z_0[j] = pa_x[j] * tlz_yyz_z_0[j] - 0.5 * fl1_fx * tpy_yyz_z_0[j] - fl1_fx * fl1_fgb * tdy_yyz_z_0[j];

            tlx_xyzz_x_0[j] = pa_x[j] * tlx_yzz_x_0[j] + 0.5 * fl1_fx * tlx_yzz_0_0[j];

            tly_xyzz_x_0[j] =
                pa_x[j] * tly_yzz_x_0[j] + 0.5 * fl1_fx * tly_yzz_0_0[j] + 0.5 * fl1_fx * tpz_yzz_x_0[j] + fl1_fx * fl1_fgb * tdz_yzz_x_0[j];

            tlz_xyzz_x_0[j] =
                pa_x[j] * tlz_yzz_x_0[j] + 0.5 * fl1_fx * tlz_yzz_0_0[j] - 0.5 * fl1_fx * tpy_yzz_x_0[j] - fl1_fx * fl1_fgb * tdy_yzz_x_0[j];

            tlx_xyzz_y_0[j] = pa_x[j] * tlx_yzz_y_0[j];

            tly_xyzz_y_0[j] = pa_x[j] * tly_yzz_y_0[j] + 0.5 * fl1_fx * tpz_yzz_y_0[j] + fl1_fx * fl1_fgb * tdz_yzz_y_0[j];

            tlz_xyzz_y_0[j] = pa_x[j] * tlz_yzz_y_0[j] - 0.5 * fl1_fx * tpy_yzz_y_0[j] - fl1_fx * fl1_fgb * tdy_yzz_y_0[j];

            tlx_xyzz_z_0[j] = pa_x[j] * tlx_yzz_z_0[j];

            tly_xyzz_z_0[j] = pa_x[j] * tly_yzz_z_0[j] + 0.5 * fl1_fx * tpz_yzz_z_0[j] + fl1_fx * fl1_fgb * tdz_yzz_z_0[j];

            tlz_xyzz_z_0[j] = pa_x[j] * tlz_yzz_z_0[j] - 0.5 * fl1_fx * tpy_yzz_z_0[j] - fl1_fx * fl1_fgb * tdy_yzz_z_0[j];

            tlx_xzzz_x_0[j] = pa_x[j] * tlx_zzz_x_0[j] + 0.5 * fl1_fx * tlx_zzz_0_0[j];

            tly_xzzz_x_0[j] =
                pa_x[j] * tly_zzz_x_0[j] + 0.5 * fl1_fx * tly_zzz_0_0[j] + 0.5 * fl1_fx * tpz_zzz_x_0[j] + fl1_fx * fl1_fgb * tdz_zzz_x_0[j];

            tlz_xzzz_x_0[j] =
                pa_x[j] * tlz_zzz_x_0[j] + 0.5 * fl1_fx * tlz_zzz_0_0[j] - 0.5 * fl1_fx * tpy_zzz_x_0[j] - fl1_fx * fl1_fgb * tdy_zzz_x_0[j];

            tlx_xzzz_y_0[j] = pa_x[j] * tlx_zzz_y_0[j];

            tly_xzzz_y_0[j] = pa_x[j] * tly_zzz_y_0[j] + 0.5 * fl1_fx * tpz_zzz_y_0[j] + fl1_fx * fl1_fgb * tdz_zzz_y_0[j];

            tlz_xzzz_y_0[j] = pa_x[j] * tlz_zzz_y_0[j] - 0.5 * fl1_fx * tpy_zzz_y_0[j] - fl1_fx * fl1_fgb * tdy_zzz_y_0[j];

            tlx_xzzz_z_0[j] = pa_x[j] * tlx_zzz_z_0[j];

            tly_xzzz_z_0[j] = pa_x[j] * tly_zzz_z_0[j] + 0.5 * fl1_fx * tpz_zzz_z_0[j] + fl1_fx * fl1_fgb * tdz_zzz_z_0[j];

            tlz_xzzz_z_0[j] = pa_x[j] * tlz_zzz_z_0[j] - 0.5 * fl1_fx * tpy_zzz_z_0[j] - fl1_fx * fl1_fgb * tdy_zzz_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGP_90_135(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (90,135)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 18);

        auto tly_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tlz_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tlx_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 19);

        auto tly_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tlz_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tlx_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 20);

        auto tly_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tlz_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tlx_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 21);

        auto tly_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tlz_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tlx_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 22);

        auto tly_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tlz_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tlx_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 23);

        auto tly_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tlz_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tlx_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 24);

        auto tly_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tlz_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tlx_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 25);

        auto tly_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tlz_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tlx_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 26);

        auto tly_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tlz_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tlx_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 27);

        auto tly_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tlz_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tlx_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 28);

        auto tly_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tlz_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tlx_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 29);

        auto tly_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tlz_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9);

        auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10);

        auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11);

        auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12);

        auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13);

        auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14);

        auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15);

        auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16);

        auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17);

        auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tlx_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 6);

        auto tly_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 7);

        auto tly_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 8);

        auto tly_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 9);

        auto tly_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 9);

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpz_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 27);

        auto tpy_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 28);

        auto tpy_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 29);

        auto tpy_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tdx_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 18);

        auto tdz_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tdx_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 19);

        auto tdz_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tdx_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 20);

        auto tdz_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tdx_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 21);

        auto tdz_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tdx_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 22);

        auto tdz_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tdx_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 23);

        auto tdz_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tdx_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 24);

        auto tdz_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tdx_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 25);

        auto tdz_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tdx_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 26);

        auto tdz_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tdx_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 27);

        auto tdy_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tdz_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tdx_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 28);

        auto tdy_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tdz_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tdx_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 29);

        auto tdy_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tdz_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 29);

        // set up pointers to integrals

        auto tlx_yyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 30);

        auto tly_yyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 30);

        auto tlz_yyyy_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 30);

        auto tlx_yyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 31);

        auto tly_yyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 31);

        auto tlz_yyyy_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 31);

        auto tlx_yyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 32);

        auto tly_yyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 32);

        auto tlz_yyyy_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 32);

        auto tlx_yyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 33);

        auto tly_yyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 33);

        auto tlz_yyyz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 33);

        auto tlx_yyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 34);

        auto tly_yyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 34);

        auto tlz_yyyz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 34);

        auto tlx_yyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 35);

        auto tly_yyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 35);

        auto tlz_yyyz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 35);

        auto tlx_yyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 36);

        auto tly_yyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 36);

        auto tlz_yyzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 36);

        auto tlx_yyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 37);

        auto tly_yyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 37);

        auto tlz_yyzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 37);

        auto tlx_yyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 38);

        auto tly_yyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 38);

        auto tlz_yyzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 38);

        auto tlx_yzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 39);

        auto tly_yzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 39);

        auto tlz_yzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 39);

        auto tlx_yzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 40);

        auto tly_yzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 40);

        auto tlz_yzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 40);

        auto tlx_yzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 41);

        auto tly_yzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 41);

        auto tlz_yzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 41);

        auto tlx_zzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 42);

        auto tly_zzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 42);

        auto tlz_zzzz_x_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 42);

        auto tlx_zzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 43);

        auto tly_zzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 43);

        auto tlz_zzzz_y_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 43);

        auto tlx_zzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * idx + 44);

        auto tly_zzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 45 * bdim + 45 * idx + 44);

        auto tlz_zzzz_z_0 = primBuffer.data(pidx_l_4_1_m0 + 90 * bdim + 45 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_yyy_x_0, tdx_yyy_y_0, tdx_yyy_z_0, tdx_yyz_x_0, \
                                     tdx_yyz_y_0, tdx_yyz_z_0, tdx_yzz_x_0, tdx_yzz_y_0, tdx_yzz_z_0, tdx_zzz_x_0, \
                                     tdx_zzz_y_0, tdx_zzz_z_0, tdy_zzz_x_0, tdy_zzz_y_0, tdy_zzz_z_0, tdz_yyy_x_0, \
                                     tdz_yyy_y_0, tdz_yyy_z_0, tdz_yyz_x_0, tdz_yyz_y_0, tdz_yyz_z_0, tdz_yzz_x_0, \
                                     tdz_yzz_y_0, tdz_yzz_z_0, tdz_zzz_x_0, tdz_zzz_y_0, tdz_zzz_z_0, tlx_yy_x_0, \
                                     tlx_yy_y_0, tlx_yy_z_0, tlx_yyy_0_0, tlx_yyy_x_0, tlx_yyy_y_0, tlx_yyy_z_0, \
                                     tlx_yyyy_x_0, tlx_yyyy_y_0, tlx_yyyy_z_0, tlx_yyyz_x_0, tlx_yyyz_y_0, tlx_yyyz_z_0, \
                                     tlx_yyz_0_0, tlx_yyz_x_0, tlx_yyz_y_0, tlx_yyz_z_0, tlx_yyzz_x_0, tlx_yyzz_y_0, \
                                     tlx_yyzz_z_0, tlx_yz_x_0, tlx_yz_y_0, tlx_yz_z_0, tlx_yzz_0_0, tlx_yzz_x_0, \
                                     tlx_yzz_y_0, tlx_yzz_z_0, tlx_yzzz_x_0, tlx_yzzz_y_0, tlx_yzzz_z_0, tlx_zz_x_0, \
                                     tlx_zz_y_0, tlx_zz_z_0, tlx_zzz_0_0, tlx_zzz_x_0, tlx_zzz_y_0, tlx_zzz_z_0, \
                                     tlx_zzzz_x_0, tlx_zzzz_y_0, tlx_zzzz_z_0, tly_yy_x_0, tly_yy_y_0, tly_yy_z_0, \
                                     tly_yyy_0_0, tly_yyy_x_0, tly_yyy_y_0, tly_yyy_z_0, tly_yyyy_x_0, tly_yyyy_y_0, \
                                     tly_yyyy_z_0, tly_yyyz_x_0, tly_yyyz_y_0, tly_yyyz_z_0, tly_yyz_0_0, tly_yyz_x_0, \
                                     tly_yyz_y_0, tly_yyz_z_0, tly_yyzz_x_0, tly_yyzz_y_0, tly_yyzz_z_0, tly_yz_x_0, \
                                     tly_yz_y_0, tly_yz_z_0, tly_yzz_0_0, tly_yzz_x_0, tly_yzz_y_0, tly_yzz_z_0, \
                                     tly_yzzz_x_0, tly_yzzz_y_0, tly_yzzz_z_0, tly_zz_x_0, tly_zz_y_0, tly_zz_z_0, \
                                     tly_zzz_0_0, tly_zzz_x_0, tly_zzz_y_0, tly_zzz_z_0, tly_zzzz_x_0, tly_zzzz_y_0, \
                                     tly_zzzz_z_0, tlz_yy_x_0, tlz_yy_y_0, tlz_yy_z_0, tlz_yyy_0_0, tlz_yyy_x_0, \
                                     tlz_yyy_y_0, tlz_yyy_z_0, tlz_yyyy_x_0, tlz_yyyy_y_0, tlz_yyyy_z_0, tlz_yyyz_x_0, \
                                     tlz_yyyz_y_0, tlz_yyyz_z_0, tlz_yyz_0_0, tlz_yyz_x_0, tlz_yyz_y_0, tlz_yyz_z_0, \
                                     tlz_yyzz_x_0, tlz_yyzz_y_0, tlz_yyzz_z_0, tlz_yz_x_0, tlz_yz_y_0, tlz_yz_z_0, \
                                     tlz_yzz_0_0, tlz_yzz_x_0, tlz_yzz_y_0, tlz_yzz_z_0, tlz_yzzz_x_0, tlz_yzzz_y_0, \
                                     tlz_yzzz_z_0, tlz_zz_x_0, tlz_zz_y_0, tlz_zz_z_0, tlz_zzz_0_0, tlz_zzz_x_0, \
                                     tlz_zzz_y_0, tlz_zzz_z_0, tlz_zzzz_x_0, tlz_zzzz_y_0, tlz_zzzz_z_0, tpx_yyy_x_0, \
                                     tpx_yyy_y_0, tpx_yyy_z_0, tpx_yyz_x_0, tpx_yyz_y_0, tpx_yyz_z_0, tpx_yzz_x_0, \
                                     tpx_yzz_y_0, tpx_yzz_z_0, tpx_zzz_x_0, tpx_zzz_y_0, tpx_zzz_z_0, tpy_zzz_x_0, \
                                     tpy_zzz_y_0, tpy_zzz_z_0, tpz_yyy_x_0, tpz_yyy_y_0, tpz_yyy_z_0, tpz_yyz_x_0, \
                                     tpz_yyz_y_0, tpz_yyz_z_0, tpz_yzz_x_0, tpz_yzz_y_0, tpz_yzz_z_0, tpz_zzz_x_0, \
                                     tpz_zzz_y_0, tpz_zzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyyy_x_0[j] =
                pa_y[j] * tlx_yyy_x_0[j] + 1.5 * fl1_fx * tlx_yy_x_0[j] - 0.5 * fl1_fx * tpz_yyy_x_0[j] - fl1_fx * fl1_fgb * tdz_yyy_x_0[j];

            tly_yyyy_x_0[j] = pa_y[j] * tly_yyy_x_0[j] + 1.5 * fl1_fx * tly_yy_x_0[j];

            tlz_yyyy_x_0[j] =
                pa_y[j] * tlz_yyy_x_0[j] + 1.5 * fl1_fx * tlz_yy_x_0[j] + 0.5 * fl1_fx * tpx_yyy_x_0[j] + fl1_fx * fl1_fgb * tdx_yyy_x_0[j];

            tlx_yyyy_y_0[j] = pa_y[j] * tlx_yyy_y_0[j] + 1.5 * fl1_fx * tlx_yy_y_0[j] + 0.5 * fl1_fx * tlx_yyy_0_0[j] -
                              0.5 * fl1_fx * tpz_yyy_y_0[j] - fl1_fx * fl1_fgb * tdz_yyy_y_0[j];

            tly_yyyy_y_0[j] = pa_y[j] * tly_yyy_y_0[j] + 1.5 * fl1_fx * tly_yy_y_0[j] + 0.5 * fl1_fx * tly_yyy_0_0[j];

            tlz_yyyy_y_0[j] = pa_y[j] * tlz_yyy_y_0[j] + 1.5 * fl1_fx * tlz_yy_y_0[j] + 0.5 * fl1_fx * tlz_yyy_0_0[j] +
                              0.5 * fl1_fx * tpx_yyy_y_0[j] + fl1_fx * fl1_fgb * tdx_yyy_y_0[j];

            tlx_yyyy_z_0[j] =
                pa_y[j] * tlx_yyy_z_0[j] + 1.5 * fl1_fx * tlx_yy_z_0[j] - 0.5 * fl1_fx * tpz_yyy_z_0[j] - fl1_fx * fl1_fgb * tdz_yyy_z_0[j];

            tly_yyyy_z_0[j] = pa_y[j] * tly_yyy_z_0[j] + 1.5 * fl1_fx * tly_yy_z_0[j];

            tlz_yyyy_z_0[j] =
                pa_y[j] * tlz_yyy_z_0[j] + 1.5 * fl1_fx * tlz_yy_z_0[j] + 0.5 * fl1_fx * tpx_yyy_z_0[j] + fl1_fx * fl1_fgb * tdx_yyy_z_0[j];

            tlx_yyyz_x_0[j] = pa_y[j] * tlx_yyz_x_0[j] + fl1_fx * tlx_yz_x_0[j] - 0.5 * fl1_fx * tpz_yyz_x_0[j] - fl1_fx * fl1_fgb * tdz_yyz_x_0[j];

            tly_yyyz_x_0[j] = pa_y[j] * tly_yyz_x_0[j] + fl1_fx * tly_yz_x_0[j];

            tlz_yyyz_x_0[j] = pa_y[j] * tlz_yyz_x_0[j] + fl1_fx * tlz_yz_x_0[j] + 0.5 * fl1_fx * tpx_yyz_x_0[j] + fl1_fx * fl1_fgb * tdx_yyz_x_0[j];

            tlx_yyyz_y_0[j] = pa_y[j] * tlx_yyz_y_0[j] + fl1_fx * tlx_yz_y_0[j] + 0.5 * fl1_fx * tlx_yyz_0_0[j] - 0.5 * fl1_fx * tpz_yyz_y_0[j] -
                              fl1_fx * fl1_fgb * tdz_yyz_y_0[j];

            tly_yyyz_y_0[j] = pa_y[j] * tly_yyz_y_0[j] + fl1_fx * tly_yz_y_0[j] + 0.5 * fl1_fx * tly_yyz_0_0[j];

            tlz_yyyz_y_0[j] = pa_y[j] * tlz_yyz_y_0[j] + fl1_fx * tlz_yz_y_0[j] + 0.5 * fl1_fx * tlz_yyz_0_0[j] + 0.5 * fl1_fx * tpx_yyz_y_0[j] +
                              fl1_fx * fl1_fgb * tdx_yyz_y_0[j];

            tlx_yyyz_z_0[j] = pa_y[j] * tlx_yyz_z_0[j] + fl1_fx * tlx_yz_z_0[j] - 0.5 * fl1_fx * tpz_yyz_z_0[j] - fl1_fx * fl1_fgb * tdz_yyz_z_0[j];

            tly_yyyz_z_0[j] = pa_y[j] * tly_yyz_z_0[j] + fl1_fx * tly_yz_z_0[j];

            tlz_yyyz_z_0[j] = pa_y[j] * tlz_yyz_z_0[j] + fl1_fx * tlz_yz_z_0[j] + 0.5 * fl1_fx * tpx_yyz_z_0[j] + fl1_fx * fl1_fgb * tdx_yyz_z_0[j];

            tlx_yyzz_x_0[j] =
                pa_y[j] * tlx_yzz_x_0[j] + 0.5 * fl1_fx * tlx_zz_x_0[j] - 0.5 * fl1_fx * tpz_yzz_x_0[j] - fl1_fx * fl1_fgb * tdz_yzz_x_0[j];

            tly_yyzz_x_0[j] = pa_y[j] * tly_yzz_x_0[j] + 0.5 * fl1_fx * tly_zz_x_0[j];

            tlz_yyzz_x_0[j] =
                pa_y[j] * tlz_yzz_x_0[j] + 0.5 * fl1_fx * tlz_zz_x_0[j] + 0.5 * fl1_fx * tpx_yzz_x_0[j] + fl1_fx * fl1_fgb * tdx_yzz_x_0[j];

            tlx_yyzz_y_0[j] = pa_y[j] * tlx_yzz_y_0[j] + 0.5 * fl1_fx * tlx_zz_y_0[j] + 0.5 * fl1_fx * tlx_yzz_0_0[j] -
                              0.5 * fl1_fx * tpz_yzz_y_0[j] - fl1_fx * fl1_fgb * tdz_yzz_y_0[j];

            tly_yyzz_y_0[j] = pa_y[j] * tly_yzz_y_0[j] + 0.5 * fl1_fx * tly_zz_y_0[j] + 0.5 * fl1_fx * tly_yzz_0_0[j];

            tlz_yyzz_y_0[j] = pa_y[j] * tlz_yzz_y_0[j] + 0.5 * fl1_fx * tlz_zz_y_0[j] + 0.5 * fl1_fx * tlz_yzz_0_0[j] +
                              0.5 * fl1_fx * tpx_yzz_y_0[j] + fl1_fx * fl1_fgb * tdx_yzz_y_0[j];

            tlx_yyzz_z_0[j] =
                pa_y[j] * tlx_yzz_z_0[j] + 0.5 * fl1_fx * tlx_zz_z_0[j] - 0.5 * fl1_fx * tpz_yzz_z_0[j] - fl1_fx * fl1_fgb * tdz_yzz_z_0[j];

            tly_yyzz_z_0[j] = pa_y[j] * tly_yzz_z_0[j] + 0.5 * fl1_fx * tly_zz_z_0[j];

            tlz_yyzz_z_0[j] =
                pa_y[j] * tlz_yzz_z_0[j] + 0.5 * fl1_fx * tlz_zz_z_0[j] + 0.5 * fl1_fx * tpx_yzz_z_0[j] + fl1_fx * fl1_fgb * tdx_yzz_z_0[j];

            tlx_yzzz_x_0[j] = pa_y[j] * tlx_zzz_x_0[j] - 0.5 * fl1_fx * tpz_zzz_x_0[j] - fl1_fx * fl1_fgb * tdz_zzz_x_0[j];

            tly_yzzz_x_0[j] = pa_y[j] * tly_zzz_x_0[j];

            tlz_yzzz_x_0[j] = pa_y[j] * tlz_zzz_x_0[j] + 0.5 * fl1_fx * tpx_zzz_x_0[j] + fl1_fx * fl1_fgb * tdx_zzz_x_0[j];

            tlx_yzzz_y_0[j] =
                pa_y[j] * tlx_zzz_y_0[j] + 0.5 * fl1_fx * tlx_zzz_0_0[j] - 0.5 * fl1_fx * tpz_zzz_y_0[j] - fl1_fx * fl1_fgb * tdz_zzz_y_0[j];

            tly_yzzz_y_0[j] = pa_y[j] * tly_zzz_y_0[j] + 0.5 * fl1_fx * tly_zzz_0_0[j];

            tlz_yzzz_y_0[j] =
                pa_y[j] * tlz_zzz_y_0[j] + 0.5 * fl1_fx * tlz_zzz_0_0[j] + 0.5 * fl1_fx * tpx_zzz_y_0[j] + fl1_fx * fl1_fgb * tdx_zzz_y_0[j];

            tlx_yzzz_z_0[j] = pa_y[j] * tlx_zzz_z_0[j] - 0.5 * fl1_fx * tpz_zzz_z_0[j] - fl1_fx * fl1_fgb * tdz_zzz_z_0[j];

            tly_yzzz_z_0[j] = pa_y[j] * tly_zzz_z_0[j];

            tlz_yzzz_z_0[j] = pa_y[j] * tlz_zzz_z_0[j] + 0.5 * fl1_fx * tpx_zzz_z_0[j] + fl1_fx * fl1_fgb * tdx_zzz_z_0[j];

            tlx_zzzz_x_0[j] =
                pa_z[j] * tlx_zzz_x_0[j] + 1.5 * fl1_fx * tlx_zz_x_0[j] + 0.5 * fl1_fx * tpy_zzz_x_0[j] + fl1_fx * fl1_fgb * tdy_zzz_x_0[j];

            tly_zzzz_x_0[j] =
                pa_z[j] * tly_zzz_x_0[j] + 1.5 * fl1_fx * tly_zz_x_0[j] - 0.5 * fl1_fx * tpx_zzz_x_0[j] - fl1_fx * fl1_fgb * tdx_zzz_x_0[j];

            tlz_zzzz_x_0[j] = pa_z[j] * tlz_zzz_x_0[j] + 1.5 * fl1_fx * tlz_zz_x_0[j];

            tlx_zzzz_y_0[j] =
                pa_z[j] * tlx_zzz_y_0[j] + 1.5 * fl1_fx * tlx_zz_y_0[j] + 0.5 * fl1_fx * tpy_zzz_y_0[j] + fl1_fx * fl1_fgb * tdy_zzz_y_0[j];

            tly_zzzz_y_0[j] =
                pa_z[j] * tly_zzz_y_0[j] + 1.5 * fl1_fx * tly_zz_y_0[j] - 0.5 * fl1_fx * tpx_zzz_y_0[j] - fl1_fx * fl1_fgb * tdx_zzz_y_0[j];

            tlz_zzzz_y_0[j] = pa_z[j] * tlz_zzz_y_0[j] + 1.5 * fl1_fx * tlz_zz_y_0[j];

            tlx_zzzz_z_0[j] = pa_z[j] * tlx_zzz_z_0[j] + 1.5 * fl1_fx * tlx_zz_z_0[j] + 0.5 * fl1_fx * tlx_zzz_0_0[j] +
                              0.5 * fl1_fx * tpy_zzz_z_0[j] + fl1_fx * fl1_fgb * tdy_zzz_z_0[j];

            tly_zzzz_z_0[j] = pa_z[j] * tly_zzz_z_0[j] + 1.5 * fl1_fx * tly_zz_z_0[j] + 0.5 * fl1_fx * tly_zzz_0_0[j] -
                              0.5 * fl1_fx * tpx_zzz_z_0[j] - fl1_fx * fl1_fgb * tdx_zzz_z_0[j];

            tlz_zzzz_z_0[j] = pa_z[j] * tlz_zzz_z_0[j] + 1.5 * fl1_fx * tlz_zz_z_0[j] + 0.5 * fl1_fx * tlz_zzz_0_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
