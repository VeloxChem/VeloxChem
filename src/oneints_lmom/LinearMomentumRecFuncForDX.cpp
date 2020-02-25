//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "LinearMomentumRecFuncForDX.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForDD(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForDD_0_36(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDD_36_72(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDD_72_108(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForDD_0_36(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (0,36)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx);

        auto tpy_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx);

        auto tpz_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx);

        auto tpx_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 1);

        auto tpy_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 2);

        auto tpy_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 3);

        auto tpy_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 4);

        auto tpy_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 5);

        auto tpy_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 3);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 4);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 5);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx);

        auto tpy_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx);

        auto tpz_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx);

        auto tpx_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 1);

        auto tpy_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tpz_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tpx_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 2);

        auto tpy_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tpz_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

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

        // set up pointers to integrals

        auto tpx_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx);

        auto tpy_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx);

        auto tpz_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx);

        auto tpx_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 1);

        auto tpy_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tpz_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tpx_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 2);

        auto tpy_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tpz_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tpx_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 3);

        auto tpy_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tpz_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tpx_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 4);

        auto tpy_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tpz_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tpx_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 5);

        auto tpy_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tpz_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tpx_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 6);

        auto tpy_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tpz_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tpx_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 7);

        auto tpy_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tpz_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tpx_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 8);

        auto tpy_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tpz_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tpx_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 9);

        auto tpy_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tpz_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tpx_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 10);

        auto tpy_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tpz_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tpx_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 11);

        auto tpy_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tpz_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 11);

        // Batch of Integrals (0,36)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, tpx_0_yy_0, tpx_0_yz_0, \
                                     tpx_0_zz_0, tpx_x_x_0, tpx_x_xx_0, tpx_x_xy_0, tpx_x_xz_0, tpx_x_y_0, tpx_x_yy_0, \
                                     tpx_x_yz_0, tpx_x_z_0, tpx_x_zz_0, tpx_xx_xx_0, tpx_xx_xy_0, tpx_xx_xz_0, \
                                     tpx_xx_yy_0, tpx_xx_yz_0, tpx_xx_zz_0, tpx_xy_xx_0, tpx_xy_xy_0, tpx_xy_xz_0, \
                                     tpx_xy_yy_0, tpx_xy_yz_0, tpx_xy_zz_0, tpx_y_x_0, tpx_y_xx_0, tpx_y_xy_0, \
                                     tpx_y_xz_0, tpx_y_y_0, tpx_y_yy_0, tpx_y_yz_0, tpx_y_z_0, tpx_y_zz_0, tpy_0_xx_0, \
                                     tpy_0_xy_0, tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, tpy_x_x_0, tpy_x_xx_0, \
                                     tpy_x_xy_0, tpy_x_xz_0, tpy_x_y_0, tpy_x_yy_0, tpy_x_yz_0, tpy_x_z_0, tpy_x_zz_0, \
                                     tpy_xx_xx_0, tpy_xx_xy_0, tpy_xx_xz_0, tpy_xx_yy_0, tpy_xx_yz_0, tpy_xx_zz_0, \
                                     tpy_xy_xx_0, tpy_xy_xy_0, tpy_xy_xz_0, tpy_xy_yy_0, tpy_xy_yz_0, tpy_xy_zz_0, \
                                     tpy_y_x_0, tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpy_y_y_0, tpy_y_yy_0, tpy_y_yz_0, \
                                     tpy_y_z_0, tpy_y_zz_0, tpz_0_xx_0, tpz_0_xy_0, tpz_0_xz_0, tpz_0_yy_0, tpz_0_yz_0, \
                                     tpz_0_zz_0, tpz_x_x_0, tpz_x_xx_0, tpz_x_xy_0, tpz_x_xz_0, tpz_x_y_0, tpz_x_yy_0, \
                                     tpz_x_yz_0, tpz_x_z_0, tpz_x_zz_0, tpz_xx_xx_0, tpz_xx_xy_0, tpz_xx_xz_0, \
                                     tpz_xx_yy_0, tpz_xx_yz_0, tpz_xx_zz_0, tpz_xy_xx_0, tpz_xy_xy_0, tpz_xy_xz_0, \
                                     tpz_xy_yy_0, tpz_xy_yz_0, tpz_xy_zz_0, tpz_y_x_0, tpz_y_xx_0, tpz_y_xy_0, \
                                     tpz_y_xz_0, tpz_y_y_0, tpz_y_yy_0, tpz_y_yz_0, tpz_y_z_0, tpz_y_zz_0, ts_x_xx_0, \
                                     ts_x_xy_0, ts_x_xz_0, ts_x_yy_0, ts_x_yz_0, ts_x_zz_0, ts_y_xx_0, ts_y_xy_0, \
                                     ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, ts_y_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xx_xx_0[j] = pa_x[j] * tpx_x_xx_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j] + fl1_fx * tpx_x_x_0[j] - fl1_fgb * fl1_fx * ts_x_xx_0[j];

            tpy_xx_xx_0[j] = pa_x[j] * tpy_x_xx_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j] + fl1_fx * tpy_x_x_0[j];

            tpz_xx_xx_0[j] = pa_x[j] * tpz_x_xx_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j] + fl1_fx * tpz_x_x_0[j];

            tpx_xx_xy_0[j] = pa_x[j] * tpx_x_xy_0[j] + 0.5 * fl1_fx * tpx_0_xy_0[j] + 0.5 * fl1_fx * tpx_x_y_0[j] - fl1_fgb * fl1_fx * ts_x_xy_0[j];

            tpy_xx_xy_0[j] = pa_x[j] * tpy_x_xy_0[j] + 0.5 * fl1_fx * tpy_0_xy_0[j] + 0.5 * fl1_fx * tpy_x_y_0[j];

            tpz_xx_xy_0[j] = pa_x[j] * tpz_x_xy_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] + 0.5 * fl1_fx * tpz_x_y_0[j];

            tpx_xx_xz_0[j] = pa_x[j] * tpx_x_xz_0[j] + 0.5 * fl1_fx * tpx_0_xz_0[j] + 0.5 * fl1_fx * tpx_x_z_0[j] - fl1_fgb * fl1_fx * ts_x_xz_0[j];

            tpy_xx_xz_0[j] = pa_x[j] * tpy_x_xz_0[j] + 0.5 * fl1_fx * tpy_0_xz_0[j] + 0.5 * fl1_fx * tpy_x_z_0[j];

            tpz_xx_xz_0[j] = pa_x[j] * tpz_x_xz_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j] + 0.5 * fl1_fx * tpz_x_z_0[j];

            tpx_xx_yy_0[j] = pa_x[j] * tpx_x_yy_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] - fl1_fgb * fl1_fx * ts_x_yy_0[j];

            tpy_xx_yy_0[j] = pa_x[j] * tpy_x_yy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j];

            tpz_xx_yy_0[j] = pa_x[j] * tpz_x_yy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j];

            tpx_xx_yz_0[j] = pa_x[j] * tpx_x_yz_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] - fl1_fgb * fl1_fx * ts_x_yz_0[j];

            tpy_xx_yz_0[j] = pa_x[j] * tpy_x_yz_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j];

            tpz_xx_yz_0[j] = pa_x[j] * tpz_x_yz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j];

            tpx_xx_zz_0[j] = pa_x[j] * tpx_x_zz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] - fl1_fgb * fl1_fx * ts_x_zz_0[j];

            tpy_xx_zz_0[j] = pa_x[j] * tpy_x_zz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j];

            tpz_xx_zz_0[j] = pa_x[j] * tpz_x_zz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];

            tpx_xy_xx_0[j] = pa_x[j] * tpx_y_xx_0[j] + fl1_fx * tpx_y_x_0[j] - fl1_fgb * fl1_fx * ts_y_xx_0[j];

            tpy_xy_xx_0[j] = pa_x[j] * tpy_y_xx_0[j] + fl1_fx * tpy_y_x_0[j];

            tpz_xy_xx_0[j] = pa_x[j] * tpz_y_xx_0[j] + fl1_fx * tpz_y_x_0[j];

            tpx_xy_xy_0[j] = pa_x[j] * tpx_y_xy_0[j] + 0.5 * fl1_fx * tpx_y_y_0[j] - fl1_fgb * fl1_fx * ts_y_xy_0[j];

            tpy_xy_xy_0[j] = pa_x[j] * tpy_y_xy_0[j] + 0.5 * fl1_fx * tpy_y_y_0[j];

            tpz_xy_xy_0[j] = pa_x[j] * tpz_y_xy_0[j] + 0.5 * fl1_fx * tpz_y_y_0[j];

            tpx_xy_xz_0[j] = pa_x[j] * tpx_y_xz_0[j] + 0.5 * fl1_fx * tpx_y_z_0[j] - fl1_fgb * fl1_fx * ts_y_xz_0[j];

            tpy_xy_xz_0[j] = pa_x[j] * tpy_y_xz_0[j] + 0.5 * fl1_fx * tpy_y_z_0[j];

            tpz_xy_xz_0[j] = pa_x[j] * tpz_y_xz_0[j] + 0.5 * fl1_fx * tpz_y_z_0[j];

            tpx_xy_yy_0[j] = pa_x[j] * tpx_y_yy_0[j] - fl1_fgb * fl1_fx * ts_y_yy_0[j];

            tpy_xy_yy_0[j] = pa_x[j] * tpy_y_yy_0[j];

            tpz_xy_yy_0[j] = pa_x[j] * tpz_y_yy_0[j];

            tpx_xy_yz_0[j] = pa_x[j] * tpx_y_yz_0[j] - fl1_fgb * fl1_fx * ts_y_yz_0[j];

            tpy_xy_yz_0[j] = pa_x[j] * tpy_y_yz_0[j];

            tpz_xy_yz_0[j] = pa_x[j] * tpz_y_yz_0[j];

            tpx_xy_zz_0[j] = pa_x[j] * tpx_y_zz_0[j] - fl1_fgb * fl1_fx * ts_y_zz_0[j];

            tpy_xy_zz_0[j] = pa_x[j] * tpy_y_zz_0[j];

            tpz_xy_zz_0[j] = pa_x[j] * tpz_y_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDD_36_72(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (36,72)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 3);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 4);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 5);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

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

        // set up pointers to integrals

        auto tpx_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 12);

        auto tpy_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tpz_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tpx_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 13);

        auto tpy_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tpz_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tpx_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 14);

        auto tpy_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tpz_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tpx_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 15);

        auto tpy_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tpz_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tpx_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 16);

        auto tpy_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tpz_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tpx_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 17);

        auto tpy_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tpz_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        // Batch of Integrals (36,72)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, tpx_0_yy_0, \
                                     tpx_0_yz_0, tpx_0_zz_0, tpx_xz_xx_0, tpx_xz_xy_0, tpx_xz_xz_0, tpx_xz_yy_0, \
                                     tpx_xz_yz_0, tpx_xz_zz_0, tpx_y_x_0, tpx_y_xx_0, tpx_y_xy_0, tpx_y_xz_0, tpx_y_y_0, \
                                     tpx_y_yy_0, tpx_y_yz_0, tpx_y_z_0, tpx_y_zz_0, tpx_yy_xx_0, tpx_yy_xy_0, \
                                     tpx_yy_xz_0, tpx_yy_yy_0, tpx_yy_yz_0, tpx_yy_zz_0, tpx_z_x_0, tpx_z_xx_0, \
                                     tpx_z_xy_0, tpx_z_xz_0, tpx_z_y_0, tpx_z_yy_0, tpx_z_yz_0, tpx_z_z_0, tpx_z_zz_0, \
                                     tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, \
                                     tpy_xz_xx_0, tpy_xz_xy_0, tpy_xz_xz_0, tpy_xz_yy_0, tpy_xz_yz_0, tpy_xz_zz_0, \
                                     tpy_y_x_0, tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpy_y_y_0, tpy_y_yy_0, tpy_y_yz_0, \
                                     tpy_y_z_0, tpy_y_zz_0, tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_yy_yy_0, \
                                     tpy_yy_yz_0, tpy_yy_zz_0, tpy_z_x_0, tpy_z_xx_0, tpy_z_xy_0, tpy_z_xz_0, tpy_z_y_0, \
                                     tpy_z_yy_0, tpy_z_yz_0, tpy_z_z_0, tpy_z_zz_0, tpz_0_xx_0, tpz_0_xy_0, tpz_0_xz_0, \
                                     tpz_0_yy_0, tpz_0_yz_0, tpz_0_zz_0, tpz_xz_xx_0, tpz_xz_xy_0, tpz_xz_xz_0, \
                                     tpz_xz_yy_0, tpz_xz_yz_0, tpz_xz_zz_0, tpz_y_x_0, tpz_y_xx_0, tpz_y_xy_0, \
                                     tpz_y_xz_0, tpz_y_y_0, tpz_y_yy_0, tpz_y_yz_0, tpz_y_z_0, tpz_y_zz_0, tpz_yy_xx_0, \
                                     tpz_yy_xy_0, tpz_yy_xz_0, tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_zz_0, tpz_z_x_0, \
                                     tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_z_y_0, tpz_z_yy_0, tpz_z_yz_0, tpz_z_z_0, \
                                     tpz_z_zz_0, ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, ts_y_zz_0, \
                                     ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, ts_z_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xz_xx_0[j] = pa_x[j] * tpx_z_xx_0[j] + fl1_fx * tpx_z_x_0[j] - fl1_fgb * fl1_fx * ts_z_xx_0[j];

            tpy_xz_xx_0[j] = pa_x[j] * tpy_z_xx_0[j] + fl1_fx * tpy_z_x_0[j];

            tpz_xz_xx_0[j] = pa_x[j] * tpz_z_xx_0[j] + fl1_fx * tpz_z_x_0[j];

            tpx_xz_xy_0[j] = pa_x[j] * tpx_z_xy_0[j] + 0.5 * fl1_fx * tpx_z_y_0[j] - fl1_fgb * fl1_fx * ts_z_xy_0[j];

            tpy_xz_xy_0[j] = pa_x[j] * tpy_z_xy_0[j] + 0.5 * fl1_fx * tpy_z_y_0[j];

            tpz_xz_xy_0[j] = pa_x[j] * tpz_z_xy_0[j] + 0.5 * fl1_fx * tpz_z_y_0[j];

            tpx_xz_xz_0[j] = pa_x[j] * tpx_z_xz_0[j] + 0.5 * fl1_fx * tpx_z_z_0[j] - fl1_fgb * fl1_fx * ts_z_xz_0[j];

            tpy_xz_xz_0[j] = pa_x[j] * tpy_z_xz_0[j] + 0.5 * fl1_fx * tpy_z_z_0[j];

            tpz_xz_xz_0[j] = pa_x[j] * tpz_z_xz_0[j] + 0.5 * fl1_fx * tpz_z_z_0[j];

            tpx_xz_yy_0[j] = pa_x[j] * tpx_z_yy_0[j] - fl1_fgb * fl1_fx * ts_z_yy_0[j];

            tpy_xz_yy_0[j] = pa_x[j] * tpy_z_yy_0[j];

            tpz_xz_yy_0[j] = pa_x[j] * tpz_z_yy_0[j];

            tpx_xz_yz_0[j] = pa_x[j] * tpx_z_yz_0[j] - fl1_fgb * fl1_fx * ts_z_yz_0[j];

            tpy_xz_yz_0[j] = pa_x[j] * tpy_z_yz_0[j];

            tpz_xz_yz_0[j] = pa_x[j] * tpz_z_yz_0[j];

            tpx_xz_zz_0[j] = pa_x[j] * tpx_z_zz_0[j] - fl1_fgb * fl1_fx * ts_z_zz_0[j];

            tpy_xz_zz_0[j] = pa_x[j] * tpy_z_zz_0[j];

            tpz_xz_zz_0[j] = pa_x[j] * tpz_z_zz_0[j];

            tpx_yy_xx_0[j] = pa_y[j] * tpx_y_xx_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j];

            tpy_yy_xx_0[j] = pa_y[j] * tpy_y_xx_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j] - fl1_fgb * fl1_fx * ts_y_xx_0[j];

            tpz_yy_xx_0[j] = pa_y[j] * tpz_y_xx_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j];

            tpx_yy_xy_0[j] = pa_y[j] * tpx_y_xy_0[j] + 0.5 * fl1_fx * tpx_0_xy_0[j] + 0.5 * fl1_fx * tpx_y_x_0[j];

            tpy_yy_xy_0[j] = pa_y[j] * tpy_y_xy_0[j] + 0.5 * fl1_fx * tpy_0_xy_0[j] + 0.5 * fl1_fx * tpy_y_x_0[j] - fl1_fgb * fl1_fx * ts_y_xy_0[j];

            tpz_yy_xy_0[j] = pa_y[j] * tpz_y_xy_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] + 0.5 * fl1_fx * tpz_y_x_0[j];

            tpx_yy_xz_0[j] = pa_y[j] * tpx_y_xz_0[j] + 0.5 * fl1_fx * tpx_0_xz_0[j];

            tpy_yy_xz_0[j] = pa_y[j] * tpy_y_xz_0[j] + 0.5 * fl1_fx * tpy_0_xz_0[j] - fl1_fgb * fl1_fx * ts_y_xz_0[j];

            tpz_yy_xz_0[j] = pa_y[j] * tpz_y_xz_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j];

            tpx_yy_yy_0[j] = pa_y[j] * tpx_y_yy_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] + fl1_fx * tpx_y_y_0[j];

            tpy_yy_yy_0[j] = pa_y[j] * tpy_y_yy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j] + fl1_fx * tpy_y_y_0[j] - fl1_fgb * fl1_fx * ts_y_yy_0[j];

            tpz_yy_yy_0[j] = pa_y[j] * tpz_y_yy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j] + fl1_fx * tpz_y_y_0[j];

            tpx_yy_yz_0[j] = pa_y[j] * tpx_y_yz_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] + 0.5 * fl1_fx * tpx_y_z_0[j];

            tpy_yy_yz_0[j] = pa_y[j] * tpy_y_yz_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j] + 0.5 * fl1_fx * tpy_y_z_0[j] - fl1_fgb * fl1_fx * ts_y_yz_0[j];

            tpz_yy_yz_0[j] = pa_y[j] * tpz_y_yz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j] + 0.5 * fl1_fx * tpz_y_z_0[j];

            tpx_yy_zz_0[j] = pa_y[j] * tpx_y_zz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j];

            tpy_yy_zz_0[j] = pa_y[j] * tpy_y_zz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] - fl1_fgb * fl1_fx * ts_y_zz_0[j];

            tpz_yy_zz_0[j] = pa_y[j] * tpz_y_zz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDD_72_108(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (72,108)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 3);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 4);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 5);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12);

        auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13);

        auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14);

        auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15);

        auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16);

        auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17);

        // set up pointers to integrals

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        // Batch of Integrals (72,108)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, tpx_0_yy_0, \
                                     tpx_0_yz_0, tpx_0_zz_0, tpx_yz_xx_0, tpx_yz_xy_0, tpx_yz_xz_0, tpx_yz_yy_0, \
                                     tpx_yz_yz_0, tpx_yz_zz_0, tpx_z_x_0, tpx_z_xx_0, tpx_z_xy_0, tpx_z_xz_0, tpx_z_y_0, \
                                     tpx_z_yy_0, tpx_z_yz_0, tpx_z_z_0, tpx_z_zz_0, tpx_zz_xx_0, tpx_zz_xy_0, \
                                     tpx_zz_xz_0, tpx_zz_yy_0, tpx_zz_yz_0, tpx_zz_zz_0, tpy_0_xx_0, tpy_0_xy_0, \
                                     tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, tpy_yz_xx_0, tpy_yz_xy_0, \
                                     tpy_yz_xz_0, tpy_yz_yy_0, tpy_yz_yz_0, tpy_yz_zz_0, tpy_z_x_0, tpy_z_xx_0, \
                                     tpy_z_xy_0, tpy_z_xz_0, tpy_z_y_0, tpy_z_yy_0, tpy_z_yz_0, tpy_z_z_0, tpy_z_zz_0, \
                                     tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpy_zz_yy_0, tpy_zz_yz_0, tpy_zz_zz_0, \
                                     tpz_0_xx_0, tpz_0_xy_0, tpz_0_xz_0, tpz_0_yy_0, tpz_0_yz_0, tpz_0_zz_0, \
                                     tpz_yz_xx_0, tpz_yz_xy_0, tpz_yz_xz_0, tpz_yz_yy_0, tpz_yz_yz_0, tpz_yz_zz_0, \
                                     tpz_z_x_0, tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_z_y_0, tpz_z_yy_0, tpz_z_yz_0, \
                                     tpz_z_z_0, tpz_z_zz_0, tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_yy_0, \
                                     tpz_zz_yz_0, tpz_zz_zz_0, ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, \
                                     ts_z_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yz_xx_0[j] = pa_y[j] * tpx_z_xx_0[j];

            tpy_yz_xx_0[j] = pa_y[j] * tpy_z_xx_0[j] - fl1_fgb * fl1_fx * ts_z_xx_0[j];

            tpz_yz_xx_0[j] = pa_y[j] * tpz_z_xx_0[j];

            tpx_yz_xy_0[j] = pa_y[j] * tpx_z_xy_0[j] + 0.5 * fl1_fx * tpx_z_x_0[j];

            tpy_yz_xy_0[j] = pa_y[j] * tpy_z_xy_0[j] + 0.5 * fl1_fx * tpy_z_x_0[j] - fl1_fgb * fl1_fx * ts_z_xy_0[j];

            tpz_yz_xy_0[j] = pa_y[j] * tpz_z_xy_0[j] + 0.5 * fl1_fx * tpz_z_x_0[j];

            tpx_yz_xz_0[j] = pa_y[j] * tpx_z_xz_0[j];

            tpy_yz_xz_0[j] = pa_y[j] * tpy_z_xz_0[j] - fl1_fgb * fl1_fx * ts_z_xz_0[j];

            tpz_yz_xz_0[j] = pa_y[j] * tpz_z_xz_0[j];

            tpx_yz_yy_0[j] = pa_y[j] * tpx_z_yy_0[j] + fl1_fx * tpx_z_y_0[j];

            tpy_yz_yy_0[j] = pa_y[j] * tpy_z_yy_0[j] + fl1_fx * tpy_z_y_0[j] - fl1_fgb * fl1_fx * ts_z_yy_0[j];

            tpz_yz_yy_0[j] = pa_y[j] * tpz_z_yy_0[j] + fl1_fx * tpz_z_y_0[j];

            tpx_yz_yz_0[j] = pa_y[j] * tpx_z_yz_0[j] + 0.5 * fl1_fx * tpx_z_z_0[j];

            tpy_yz_yz_0[j] = pa_y[j] * tpy_z_yz_0[j] + 0.5 * fl1_fx * tpy_z_z_0[j] - fl1_fgb * fl1_fx * ts_z_yz_0[j];

            tpz_yz_yz_0[j] = pa_y[j] * tpz_z_yz_0[j] + 0.5 * fl1_fx * tpz_z_z_0[j];

            tpx_yz_zz_0[j] = pa_y[j] * tpx_z_zz_0[j];

            tpy_yz_zz_0[j] = pa_y[j] * tpy_z_zz_0[j] - fl1_fgb * fl1_fx * ts_z_zz_0[j];

            tpz_yz_zz_0[j] = pa_y[j] * tpz_z_zz_0[j];

            tpx_zz_xx_0[j] = pa_z[j] * tpx_z_xx_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j];

            tpy_zz_xx_0[j] = pa_z[j] * tpy_z_xx_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j];

            tpz_zz_xx_0[j] = pa_z[j] * tpz_z_xx_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j] - fl1_fgb * fl1_fx * ts_z_xx_0[j];

            tpx_zz_xy_0[j] = pa_z[j] * tpx_z_xy_0[j] + 0.5 * fl1_fx * tpx_0_xy_0[j];

            tpy_zz_xy_0[j] = pa_z[j] * tpy_z_xy_0[j] + 0.5 * fl1_fx * tpy_0_xy_0[j];

            tpz_zz_xy_0[j] = pa_z[j] * tpz_z_xy_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] - fl1_fgb * fl1_fx * ts_z_xy_0[j];

            tpx_zz_xz_0[j] = pa_z[j] * tpx_z_xz_0[j] + 0.5 * fl1_fx * tpx_0_xz_0[j] + 0.5 * fl1_fx * tpx_z_x_0[j];

            tpy_zz_xz_0[j] = pa_z[j] * tpy_z_xz_0[j] + 0.5 * fl1_fx * tpy_0_xz_0[j] + 0.5 * fl1_fx * tpy_z_x_0[j];

            tpz_zz_xz_0[j] = pa_z[j] * tpz_z_xz_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j] + 0.5 * fl1_fx * tpz_z_x_0[j] - fl1_fgb * fl1_fx * ts_z_xz_0[j];

            tpx_zz_yy_0[j] = pa_z[j] * tpx_z_yy_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j];

            tpy_zz_yy_0[j] = pa_z[j] * tpy_z_yy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j];

            tpz_zz_yy_0[j] = pa_z[j] * tpz_z_yy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j] - fl1_fgb * fl1_fx * ts_z_yy_0[j];

            tpx_zz_yz_0[j] = pa_z[j] * tpx_z_yz_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] + 0.5 * fl1_fx * tpx_z_y_0[j];

            tpy_zz_yz_0[j] = pa_z[j] * tpy_z_yz_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j] + 0.5 * fl1_fx * tpy_z_y_0[j];

            tpz_zz_yz_0[j] = pa_z[j] * tpz_z_yz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j] + 0.5 * fl1_fx * tpz_z_y_0[j] - fl1_fgb * fl1_fx * ts_z_yz_0[j];

            tpx_zz_zz_0[j] = pa_z[j] * tpx_z_zz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] + fl1_fx * tpx_z_z_0[j];

            tpy_zz_zz_0[j] = pa_z[j] * tpy_z_zz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] + fl1_fx * tpy_z_z_0[j];

            tpz_zz_zz_0[j] = pa_z[j] * tpz_z_zz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j] + fl1_fx * tpz_z_z_0[j] - fl1_fgb * fl1_fx * ts_z_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForDF_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDF_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDF_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDF_135_180(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForDF_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx);

        auto tpy_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx);

        auto tpz_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx);

        auto tpx_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 1);

        auto tpy_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 2);

        auto tpy_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 3);

        auto tpy_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 4);

        auto tpy_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 5);

        auto tpy_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 6);

        auto tpy_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 7);

        auto tpy_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 8);

        auto tpy_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 9);

        auto tpy_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

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

        auto tpx_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx);

        auto tpy_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx);

        auto tpz_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx);

        auto tpx_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 1);

        auto tpy_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 2);

        auto tpy_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 3);

        auto tpy_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 4);

        auto tpy_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 5);

        auto tpy_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

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

        // set up pointers to integrals

        auto tpx_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx);

        auto tpy_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx);

        auto tpz_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx);

        auto tpx_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 1);

        auto tpy_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tpx_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 2);

        auto tpy_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tpx_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 3);

        auto tpy_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tpx_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 4);

        auto tpy_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tpx_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 5);

        auto tpy_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tpx_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 6);

        auto tpy_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tpx_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 7);

        auto tpy_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tpx_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 8);

        auto tpy_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tpx_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 9);

        auto tpy_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tpx_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 10);

        auto tpy_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tpx_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 11);

        auto tpy_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto tpz_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tpx_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 12);

        auto tpy_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tpx_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 13);

        auto tpy_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tpx_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 14);

        auto tpy_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, tpx_0_xyy_0, \
                                     tpx_0_xyz_0, tpx_0_xzz_0, tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yzz_0, tpx_0_zzz_0, \
                                     tpx_x_xx_0, tpx_x_xxx_0, tpx_x_xxy_0, tpx_x_xxz_0, tpx_x_xy_0, tpx_x_xyy_0, \
                                     tpx_x_xyz_0, tpx_x_xz_0, tpx_x_xzz_0, tpx_x_yy_0, tpx_x_yyy_0, tpx_x_yyz_0, \
                                     tpx_x_yz_0, tpx_x_yzz_0, tpx_x_zz_0, tpx_x_zzz_0, tpx_xx_xxx_0, tpx_xx_xxy_0, \
                                     tpx_xx_xxz_0, tpx_xx_xyy_0, tpx_xx_xyz_0, tpx_xx_xzz_0, tpx_xx_yyy_0, tpx_xx_yyz_0, \
                                     tpx_xx_yzz_0, tpx_xx_zzz_0, tpx_xy_xxx_0, tpx_xy_xxy_0, tpx_xy_xxz_0, tpx_xy_xyy_0, \
                                     tpx_xy_xyz_0, tpx_y_xx_0, tpx_y_xxx_0, tpx_y_xxy_0, tpx_y_xxz_0, tpx_y_xy_0, \
                                     tpx_y_xyy_0, tpx_y_xyz_0, tpx_y_xz_0, tpx_y_yy_0, tpx_y_yz_0, tpy_0_xxx_0, \
                                     tpy_0_xxy_0, tpy_0_xxz_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xzz_0, tpy_0_yyy_0, \
                                     tpy_0_yyz_0, tpy_0_yzz_0, tpy_0_zzz_0, tpy_x_xx_0, tpy_x_xxx_0, tpy_x_xxy_0, \
                                     tpy_x_xxz_0, tpy_x_xy_0, tpy_x_xyy_0, tpy_x_xyz_0, tpy_x_xz_0, tpy_x_xzz_0, \
                                     tpy_x_yy_0, tpy_x_yyy_0, tpy_x_yyz_0, tpy_x_yz_0, tpy_x_yzz_0, tpy_x_zz_0, \
                                     tpy_x_zzz_0, tpy_xx_xxx_0, tpy_xx_xxy_0, tpy_xx_xxz_0, tpy_xx_xyy_0, tpy_xx_xyz_0, \
                                     tpy_xx_xzz_0, tpy_xx_yyy_0, tpy_xx_yyz_0, tpy_xx_yzz_0, tpy_xx_zzz_0, tpy_xy_xxx_0, \
                                     tpy_xy_xxy_0, tpy_xy_xxz_0, tpy_xy_xyy_0, tpy_xy_xyz_0, tpy_y_xx_0, tpy_y_xxx_0, \
                                     tpy_y_xxy_0, tpy_y_xxz_0, tpy_y_xy_0, tpy_y_xyy_0, tpy_y_xyz_0, tpy_y_xz_0, \
                                     tpy_y_yy_0, tpy_y_yz_0, tpz_0_xxx_0, tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xyy_0, \
                                     tpz_0_xyz_0, tpz_0_xzz_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yzz_0, tpz_0_zzz_0, \
                                     tpz_x_xx_0, tpz_x_xxx_0, tpz_x_xxy_0, tpz_x_xxz_0, tpz_x_xy_0, tpz_x_xyy_0, \
                                     tpz_x_xyz_0, tpz_x_xz_0, tpz_x_xzz_0, tpz_x_yy_0, tpz_x_yyy_0, tpz_x_yyz_0, \
                                     tpz_x_yz_0, tpz_x_yzz_0, tpz_x_zz_0, tpz_x_zzz_0, tpz_xx_xxx_0, tpz_xx_xxy_0, \
                                     tpz_xx_xxz_0, tpz_xx_xyy_0, tpz_xx_xyz_0, tpz_xx_xzz_0, tpz_xx_yyy_0, tpz_xx_yyz_0, \
                                     tpz_xx_yzz_0, tpz_xx_zzz_0, tpz_xy_xxx_0, tpz_xy_xxy_0, tpz_xy_xxz_0, tpz_xy_xyy_0, \
                                     tpz_xy_xyz_0, tpz_y_xx_0, tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xy_0, \
                                     tpz_y_xyy_0, tpz_y_xyz_0, tpz_y_xz_0, tpz_y_yy_0, tpz_y_yz_0, ts_x_xxx_0, \
                                     ts_x_xxy_0, ts_x_xxz_0, ts_x_xyy_0, ts_x_xyz_0, ts_x_xzz_0, ts_x_yyy_0, ts_x_yyz_0, \
                                     ts_x_yzz_0, ts_x_zzz_0, ts_y_xxx_0, ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xx_xxx_0[j] =
                pa_x[j] * tpx_x_xxx_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j] + 1.5 * fl1_fx * tpx_x_xx_0[j] - fl1_fgb * fl1_fx * ts_x_xxx_0[j];

            tpy_xx_xxx_0[j] = pa_x[j] * tpy_x_xxx_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j] + 1.5 * fl1_fx * tpy_x_xx_0[j];

            tpz_xx_xxx_0[j] = pa_x[j] * tpz_x_xxx_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j] + 1.5 * fl1_fx * tpz_x_xx_0[j];

            tpx_xx_xxy_0[j] = pa_x[j] * tpx_x_xxy_0[j] + 0.5 * fl1_fx * tpx_0_xxy_0[j] + fl1_fx * tpx_x_xy_0[j] - fl1_fgb * fl1_fx * ts_x_xxy_0[j];

            tpy_xx_xxy_0[j] = pa_x[j] * tpy_x_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xxy_0[j] + fl1_fx * tpy_x_xy_0[j];

            tpz_xx_xxy_0[j] = pa_x[j] * tpz_x_xxy_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] + fl1_fx * tpz_x_xy_0[j];

            tpx_xx_xxz_0[j] = pa_x[j] * tpx_x_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xxz_0[j] + fl1_fx * tpx_x_xz_0[j] - fl1_fgb * fl1_fx * ts_x_xxz_0[j];

            tpy_xx_xxz_0[j] = pa_x[j] * tpy_x_xxz_0[j] + 0.5 * fl1_fx * tpy_0_xxz_0[j] + fl1_fx * tpy_x_xz_0[j];

            tpz_xx_xxz_0[j] = pa_x[j] * tpz_x_xxz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j] + fl1_fx * tpz_x_xz_0[j];

            tpx_xx_xyy_0[j] =
                pa_x[j] * tpx_x_xyy_0[j] + 0.5 * fl1_fx * tpx_0_xyy_0[j] + 0.5 * fl1_fx * tpx_x_yy_0[j] - fl1_fgb * fl1_fx * ts_x_xyy_0[j];

            tpy_xx_xyy_0[j] = pa_x[j] * tpy_x_xyy_0[j] + 0.5 * fl1_fx * tpy_0_xyy_0[j] + 0.5 * fl1_fx * tpy_x_yy_0[j];

            tpz_xx_xyy_0[j] = pa_x[j] * tpz_x_xyy_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] + 0.5 * fl1_fx * tpz_x_yy_0[j];

            tpx_xx_xyz_0[j] =
                pa_x[j] * tpx_x_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_x_yz_0[j] - fl1_fgb * fl1_fx * ts_x_xyz_0[j];

            tpy_xx_xyz_0[j] = pa_x[j] * tpy_x_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_x_yz_0[j];

            tpz_xx_xyz_0[j] = pa_x[j] * tpz_x_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_x_yz_0[j];

            tpx_xx_xzz_0[j] =
                pa_x[j] * tpx_x_xzz_0[j] + 0.5 * fl1_fx * tpx_0_xzz_0[j] + 0.5 * fl1_fx * tpx_x_zz_0[j] - fl1_fgb * fl1_fx * ts_x_xzz_0[j];

            tpy_xx_xzz_0[j] = pa_x[j] * tpy_x_xzz_0[j] + 0.5 * fl1_fx * tpy_0_xzz_0[j] + 0.5 * fl1_fx * tpy_x_zz_0[j];

            tpz_xx_xzz_0[j] = pa_x[j] * tpz_x_xzz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j] + 0.5 * fl1_fx * tpz_x_zz_0[j];

            tpx_xx_yyy_0[j] = pa_x[j] * tpx_x_yyy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_x_yyy_0[j];

            tpy_xx_yyy_0[j] = pa_x[j] * tpy_x_yyy_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j];

            tpz_xx_yyy_0[j] = pa_x[j] * tpz_x_yyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j];

            tpx_xx_yyz_0[j] = pa_x[j] * tpx_x_yyz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] - fl1_fgb * fl1_fx * ts_x_yyz_0[j];

            tpy_xx_yyz_0[j] = pa_x[j] * tpy_x_yyz_0[j] + 0.5 * fl1_fx * tpy_0_yyz_0[j];

            tpz_xx_yyz_0[j] = pa_x[j] * tpz_x_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j];

            tpx_xx_yzz_0[j] = pa_x[j] * tpx_x_yzz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] - fl1_fgb * fl1_fx * ts_x_yzz_0[j];

            tpy_xx_yzz_0[j] = pa_x[j] * tpy_x_yzz_0[j] + 0.5 * fl1_fx * tpy_0_yzz_0[j];

            tpz_xx_yzz_0[j] = pa_x[j] * tpz_x_yzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j];

            tpx_xx_zzz_0[j] = pa_x[j] * tpx_x_zzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_x_zzz_0[j];

            tpy_xx_zzz_0[j] = pa_x[j] * tpy_x_zzz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j];

            tpz_xx_zzz_0[j] = pa_x[j] * tpz_x_zzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j];

            tpx_xy_xxx_0[j] = pa_x[j] * tpx_y_xxx_0[j] + 1.5 * fl1_fx * tpx_y_xx_0[j] - fl1_fgb * fl1_fx * ts_y_xxx_0[j];

            tpy_xy_xxx_0[j] = pa_x[j] * tpy_y_xxx_0[j] + 1.5 * fl1_fx * tpy_y_xx_0[j];

            tpz_xy_xxx_0[j] = pa_x[j] * tpz_y_xxx_0[j] + 1.5 * fl1_fx * tpz_y_xx_0[j];

            tpx_xy_xxy_0[j] = pa_x[j] * tpx_y_xxy_0[j] + fl1_fx * tpx_y_xy_0[j] - fl1_fgb * fl1_fx * ts_y_xxy_0[j];

            tpy_xy_xxy_0[j] = pa_x[j] * tpy_y_xxy_0[j] + fl1_fx * tpy_y_xy_0[j];

            tpz_xy_xxy_0[j] = pa_x[j] * tpz_y_xxy_0[j] + fl1_fx * tpz_y_xy_0[j];

            tpx_xy_xxz_0[j] = pa_x[j] * tpx_y_xxz_0[j] + fl1_fx * tpx_y_xz_0[j] - fl1_fgb * fl1_fx * ts_y_xxz_0[j];

            tpy_xy_xxz_0[j] = pa_x[j] * tpy_y_xxz_0[j] + fl1_fx * tpy_y_xz_0[j];

            tpz_xy_xxz_0[j] = pa_x[j] * tpz_y_xxz_0[j] + fl1_fx * tpz_y_xz_0[j];

            tpx_xy_xyy_0[j] = pa_x[j] * tpx_y_xyy_0[j] + 0.5 * fl1_fx * tpx_y_yy_0[j] - fl1_fgb * fl1_fx * ts_y_xyy_0[j];

            tpy_xy_xyy_0[j] = pa_x[j] * tpy_y_xyy_0[j] + 0.5 * fl1_fx * tpy_y_yy_0[j];

            tpz_xy_xyy_0[j] = pa_x[j] * tpz_y_xyy_0[j] + 0.5 * fl1_fx * tpz_y_yy_0[j];

            tpx_xy_xyz_0[j] = pa_x[j] * tpx_y_xyz_0[j] + 0.5 * fl1_fx * tpx_y_yz_0[j] - fl1_fgb * fl1_fx * ts_y_xyz_0[j];

            tpy_xy_xyz_0[j] = pa_x[j] * tpy_y_xyz_0[j] + 0.5 * fl1_fx * tpy_y_yz_0[j];

            tpz_xy_xyz_0[j] = pa_x[j] * tpz_y_xyz_0[j] + 0.5 * fl1_fx * tpz_y_yz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDF_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

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

        // set up pointers to integrals

        auto tpx_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 15);

        auto tpy_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tpx_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 16);

        auto tpy_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tpz_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 16);

        auto tpx_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 17);

        auto tpy_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tpx_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 18);

        auto tpy_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tpx_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 19);

        auto tpy_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tpx_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 20);

        auto tpy_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tpx_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 21);

        auto tpy_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tpx_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 22);

        auto tpy_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tpx_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 23);

        auto tpy_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tpx_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 24);

        auto tpy_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tpx_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 25);

        auto tpy_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tpx_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 26);

        auto tpy_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tpx_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 27);

        auto tpy_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tpx_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 28);

        auto tpy_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tpx_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 29);

        auto tpy_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xy_xzz_0, tpx_xy_yyy_0, tpx_xy_yyz_0, tpx_xy_yzz_0, \
                                     tpx_xy_zzz_0, tpx_xz_xxx_0, tpx_xz_xxy_0, tpx_xz_xxz_0, tpx_xz_xyy_0, tpx_xz_xyz_0, \
                                     tpx_xz_xzz_0, tpx_xz_yyy_0, tpx_xz_yyz_0, tpx_xz_yzz_0, tpx_xz_zzz_0, tpx_y_xzz_0, \
                                     tpx_y_yyy_0, tpx_y_yyz_0, tpx_y_yzz_0, tpx_y_zz_0, tpx_y_zzz_0, tpx_z_xx_0, \
                                     tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xy_0, tpx_z_xyy_0, tpx_z_xyz_0, \
                                     tpx_z_xz_0, tpx_z_xzz_0, tpx_z_yy_0, tpx_z_yyy_0, tpx_z_yyz_0, tpx_z_yz_0, \
                                     tpx_z_yzz_0, tpx_z_zz_0, tpx_z_zzz_0, tpy_xy_xzz_0, tpy_xy_yyy_0, tpy_xy_yyz_0, \
                                     tpy_xy_yzz_0, tpy_xy_zzz_0, tpy_xz_xxx_0, tpy_xz_xxy_0, tpy_xz_xxz_0, tpy_xz_xyy_0, \
                                     tpy_xz_xyz_0, tpy_xz_xzz_0, tpy_xz_yyy_0, tpy_xz_yyz_0, tpy_xz_yzz_0, tpy_xz_zzz_0, \
                                     tpy_y_xzz_0, tpy_y_yyy_0, tpy_y_yyz_0, tpy_y_yzz_0, tpy_y_zz_0, tpy_y_zzz_0, \
                                     tpy_z_xx_0, tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, tpy_z_xy_0, tpy_z_xyy_0, \
                                     tpy_z_xyz_0, tpy_z_xz_0, tpy_z_xzz_0, tpy_z_yy_0, tpy_z_yyy_0, tpy_z_yyz_0, \
                                     tpy_z_yz_0, tpy_z_yzz_0, tpy_z_zz_0, tpy_z_zzz_0, tpz_xy_xzz_0, tpz_xy_yyy_0, \
                                     tpz_xy_yyz_0, tpz_xy_yzz_0, tpz_xy_zzz_0, tpz_xz_xxx_0, tpz_xz_xxy_0, tpz_xz_xxz_0, \
                                     tpz_xz_xyy_0, tpz_xz_xyz_0, tpz_xz_xzz_0, tpz_xz_yyy_0, tpz_xz_yyz_0, tpz_xz_yzz_0, \
                                     tpz_xz_zzz_0, tpz_y_xzz_0, tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yzz_0, tpz_y_zz_0, \
                                     tpz_y_zzz_0, tpz_z_xx_0, tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xy_0, \
                                     tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xz_0, tpz_z_xzz_0, tpz_z_yy_0, tpz_z_yyy_0, \
                                     tpz_z_yyz_0, tpz_z_yz_0, tpz_z_yzz_0, tpz_z_zz_0, tpz_z_zzz_0, ts_y_xzz_0, \
                                     ts_y_yyy_0, ts_y_yyz_0, ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, \
                                     ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xy_xzz_0[j] = pa_x[j] * tpx_y_xzz_0[j] + 0.5 * fl1_fx * tpx_y_zz_0[j] - fl1_fgb * fl1_fx * ts_y_xzz_0[j];

            tpy_xy_xzz_0[j] = pa_x[j] * tpy_y_xzz_0[j] + 0.5 * fl1_fx * tpy_y_zz_0[j];

            tpz_xy_xzz_0[j] = pa_x[j] * tpz_y_xzz_0[j] + 0.5 * fl1_fx * tpz_y_zz_0[j];

            tpx_xy_yyy_0[j] = pa_x[j] * tpx_y_yyy_0[j] - fl1_fgb * fl1_fx * ts_y_yyy_0[j];

            tpy_xy_yyy_0[j] = pa_x[j] * tpy_y_yyy_0[j];

            tpz_xy_yyy_0[j] = pa_x[j] * tpz_y_yyy_0[j];

            tpx_xy_yyz_0[j] = pa_x[j] * tpx_y_yyz_0[j] - fl1_fgb * fl1_fx * ts_y_yyz_0[j];

            tpy_xy_yyz_0[j] = pa_x[j] * tpy_y_yyz_0[j];

            tpz_xy_yyz_0[j] = pa_x[j] * tpz_y_yyz_0[j];

            tpx_xy_yzz_0[j] = pa_x[j] * tpx_y_yzz_0[j] - fl1_fgb * fl1_fx * ts_y_yzz_0[j];

            tpy_xy_yzz_0[j] = pa_x[j] * tpy_y_yzz_0[j];

            tpz_xy_yzz_0[j] = pa_x[j] * tpz_y_yzz_0[j];

            tpx_xy_zzz_0[j] = pa_x[j] * tpx_y_zzz_0[j] - fl1_fgb * fl1_fx * ts_y_zzz_0[j];

            tpy_xy_zzz_0[j] = pa_x[j] * tpy_y_zzz_0[j];

            tpz_xy_zzz_0[j] = pa_x[j] * tpz_y_zzz_0[j];

            tpx_xz_xxx_0[j] = pa_x[j] * tpx_z_xxx_0[j] + 1.5 * fl1_fx * tpx_z_xx_0[j] - fl1_fgb * fl1_fx * ts_z_xxx_0[j];

            tpy_xz_xxx_0[j] = pa_x[j] * tpy_z_xxx_0[j] + 1.5 * fl1_fx * tpy_z_xx_0[j];

            tpz_xz_xxx_0[j] = pa_x[j] * tpz_z_xxx_0[j] + 1.5 * fl1_fx * tpz_z_xx_0[j];

            tpx_xz_xxy_0[j] = pa_x[j] * tpx_z_xxy_0[j] + fl1_fx * tpx_z_xy_0[j] - fl1_fgb * fl1_fx * ts_z_xxy_0[j];

            tpy_xz_xxy_0[j] = pa_x[j] * tpy_z_xxy_0[j] + fl1_fx * tpy_z_xy_0[j];

            tpz_xz_xxy_0[j] = pa_x[j] * tpz_z_xxy_0[j] + fl1_fx * tpz_z_xy_0[j];

            tpx_xz_xxz_0[j] = pa_x[j] * tpx_z_xxz_0[j] + fl1_fx * tpx_z_xz_0[j] - fl1_fgb * fl1_fx * ts_z_xxz_0[j];

            tpy_xz_xxz_0[j] = pa_x[j] * tpy_z_xxz_0[j] + fl1_fx * tpy_z_xz_0[j];

            tpz_xz_xxz_0[j] = pa_x[j] * tpz_z_xxz_0[j] + fl1_fx * tpz_z_xz_0[j];

            tpx_xz_xyy_0[j] = pa_x[j] * tpx_z_xyy_0[j] + 0.5 * fl1_fx * tpx_z_yy_0[j] - fl1_fgb * fl1_fx * ts_z_xyy_0[j];

            tpy_xz_xyy_0[j] = pa_x[j] * tpy_z_xyy_0[j] + 0.5 * fl1_fx * tpy_z_yy_0[j];

            tpz_xz_xyy_0[j] = pa_x[j] * tpz_z_xyy_0[j] + 0.5 * fl1_fx * tpz_z_yy_0[j];

            tpx_xz_xyz_0[j] = pa_x[j] * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_z_yz_0[j] - fl1_fgb * fl1_fx * ts_z_xyz_0[j];

            tpy_xz_xyz_0[j] = pa_x[j] * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_z_yz_0[j];

            tpz_xz_xyz_0[j] = pa_x[j] * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_z_yz_0[j];

            tpx_xz_xzz_0[j] = pa_x[j] * tpx_z_xzz_0[j] + 0.5 * fl1_fx * tpx_z_zz_0[j] - fl1_fgb * fl1_fx * ts_z_xzz_0[j];

            tpy_xz_xzz_0[j] = pa_x[j] * tpy_z_xzz_0[j] + 0.5 * fl1_fx * tpy_z_zz_0[j];

            tpz_xz_xzz_0[j] = pa_x[j] * tpz_z_xzz_0[j] + 0.5 * fl1_fx * tpz_z_zz_0[j];

            tpx_xz_yyy_0[j] = pa_x[j] * tpx_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyy_0[j];

            tpy_xz_yyy_0[j] = pa_x[j] * tpy_z_yyy_0[j];

            tpz_xz_yyy_0[j] = pa_x[j] * tpz_z_yyy_0[j];

            tpx_xz_yyz_0[j] = pa_x[j] * tpx_z_yyz_0[j] - fl1_fgb * fl1_fx * ts_z_yyz_0[j];

            tpy_xz_yyz_0[j] = pa_x[j] * tpy_z_yyz_0[j];

            tpz_xz_yyz_0[j] = pa_x[j] * tpz_z_yyz_0[j];

            tpx_xz_yzz_0[j] = pa_x[j] * tpx_z_yzz_0[j] - fl1_fgb * fl1_fx * ts_z_yzz_0[j];

            tpy_xz_yzz_0[j] = pa_x[j] * tpy_z_yzz_0[j];

            tpz_xz_yzz_0[j] = pa_x[j] * tpz_z_yzz_0[j];

            tpx_xz_zzz_0[j] = pa_x[j] * tpx_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_z_zzz_0[j];

            tpy_xz_zzz_0[j] = pa_x[j] * tpy_z_zzz_0[j];

            tpz_xz_zzz_0[j] = pa_x[j] * tpz_z_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDF_90_135(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

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

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

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

        // set up pointers to integrals

        auto tpx_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 30);

        auto tpy_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 31);

        auto tpy_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 32);

        auto tpy_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 33);

        auto tpy_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 34);

        auto tpy_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 35);

        auto tpy_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tpx_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 36);

        auto tpy_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 37);

        auto tpy_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 38);

        auto tpy_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 39);

        auto tpy_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 40);

        auto tpy_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 41);

        auto tpy_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 42);

        auto tpy_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 43);

        auto tpy_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 44);

        auto tpy_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, tpx_0_xyy_0, \
                                     tpx_0_xyz_0, tpx_0_xzz_0, tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yzz_0, tpx_0_zzz_0, \
                                     tpx_y_xx_0, tpx_y_xxx_0, tpx_y_xxy_0, tpx_y_xxz_0, tpx_y_xy_0, tpx_y_xyy_0, \
                                     tpx_y_xyz_0, tpx_y_xz_0, tpx_y_xzz_0, tpx_y_yy_0, tpx_y_yyy_0, tpx_y_yyz_0, \
                                     tpx_y_yz_0, tpx_y_yzz_0, tpx_y_zz_0, tpx_y_zzz_0, tpx_yy_xxx_0, tpx_yy_xxy_0, \
                                     tpx_yy_xxz_0, tpx_yy_xyy_0, tpx_yy_xyz_0, tpx_yy_xzz_0, tpx_yy_yyy_0, tpx_yy_yyz_0, \
                                     tpx_yy_yzz_0, tpx_yy_zzz_0, tpx_yz_xxx_0, tpx_yz_xxy_0, tpx_yz_xxz_0, tpx_yz_xyy_0, \
                                     tpx_yz_xyz_0, tpx_z_xx_0, tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xy_0, \
                                     tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xz_0, tpy_0_xxx_0, tpy_0_xxy_0, tpy_0_xxz_0, \
                                     tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xzz_0, tpy_0_yyy_0, tpy_0_yyz_0, tpy_0_yzz_0, \
                                     tpy_0_zzz_0, tpy_y_xx_0, tpy_y_xxx_0, tpy_y_xxy_0, tpy_y_xxz_0, tpy_y_xy_0, \
                                     tpy_y_xyy_0, tpy_y_xyz_0, tpy_y_xz_0, tpy_y_xzz_0, tpy_y_yy_0, tpy_y_yyy_0, \
                                     tpy_y_yyz_0, tpy_y_yz_0, tpy_y_yzz_0, tpy_y_zz_0, tpy_y_zzz_0, tpy_yy_xxx_0, \
                                     tpy_yy_xxy_0, tpy_yy_xxz_0, tpy_yy_xyy_0, tpy_yy_xyz_0, tpy_yy_xzz_0, tpy_yy_yyy_0, \
                                     tpy_yy_yyz_0, tpy_yy_yzz_0, tpy_yy_zzz_0, tpy_yz_xxx_0, tpy_yz_xxy_0, tpy_yz_xxz_0, \
                                     tpy_yz_xyy_0, tpy_yz_xyz_0, tpy_z_xx_0, tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, \
                                     tpy_z_xy_0, tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xz_0, tpz_0_xxx_0, tpz_0_xxy_0, \
                                     tpz_0_xxz_0, tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xzz_0, tpz_0_yyy_0, tpz_0_yyz_0, \
                                     tpz_0_yzz_0, tpz_0_zzz_0, tpz_y_xx_0, tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, \
                                     tpz_y_xy_0, tpz_y_xyy_0, tpz_y_xyz_0, tpz_y_xz_0, tpz_y_xzz_0, tpz_y_yy_0, \
                                     tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yz_0, tpz_y_yzz_0, tpz_y_zz_0, tpz_y_zzz_0, \
                                     tpz_yy_xxx_0, tpz_yy_xxy_0, tpz_yy_xxz_0, tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_yy_xzz_0, \
                                     tpz_yy_yyy_0, tpz_yy_yyz_0, tpz_yy_yzz_0, tpz_yy_zzz_0, tpz_yz_xxx_0, tpz_yz_xxy_0, \
                                     tpz_yz_xxz_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_z_xx_0, tpz_z_xxx_0, tpz_z_xxy_0, \
                                     tpz_z_xxz_0, tpz_z_xy_0, tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xz_0, ts_y_xxx_0, \
                                     ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, ts_y_yyy_0, ts_y_yyz_0, \
                                     ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yy_xxx_0[j] = pa_y[j] * tpx_y_xxx_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j];

            tpy_yy_xxx_0[j] = pa_y[j] * tpy_y_xxx_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_y_xxx_0[j];

            tpz_yy_xxx_0[j] = pa_y[j] * tpz_y_xxx_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j];

            tpx_yy_xxy_0[j] = pa_y[j] * tpx_y_xxy_0[j] + 0.5 * fl1_fx * tpx_0_xxy_0[j] + 0.5 * fl1_fx * tpx_y_xx_0[j];

            tpy_yy_xxy_0[j] =
                pa_y[j] * tpy_y_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xxy_0[j] + 0.5 * fl1_fx * tpy_y_xx_0[j] - fl1_fgb * fl1_fx * ts_y_xxy_0[j];

            tpz_yy_xxy_0[j] = pa_y[j] * tpz_y_xxy_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] + 0.5 * fl1_fx * tpz_y_xx_0[j];

            tpx_yy_xxz_0[j] = pa_y[j] * tpx_y_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xxz_0[j];

            tpy_yy_xxz_0[j] = pa_y[j] * tpy_y_xxz_0[j] + 0.5 * fl1_fx * tpy_0_xxz_0[j] - fl1_fgb * fl1_fx * ts_y_xxz_0[j];

            tpz_yy_xxz_0[j] = pa_y[j] * tpz_y_xxz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j];

            tpx_yy_xyy_0[j] = pa_y[j] * tpx_y_xyy_0[j] + 0.5 * fl1_fx * tpx_0_xyy_0[j] + fl1_fx * tpx_y_xy_0[j];

            tpy_yy_xyy_0[j] = pa_y[j] * tpy_y_xyy_0[j] + 0.5 * fl1_fx * tpy_0_xyy_0[j] + fl1_fx * tpy_y_xy_0[j] - fl1_fgb * fl1_fx * ts_y_xyy_0[j];

            tpz_yy_xyy_0[j] = pa_y[j] * tpz_y_xyy_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] + fl1_fx * tpz_y_xy_0[j];

            tpx_yy_xyz_0[j] = pa_y[j] * tpx_y_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_y_xz_0[j];

            tpy_yy_xyz_0[j] =
                pa_y[j] * tpy_y_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_y_xz_0[j] - fl1_fgb * fl1_fx * ts_y_xyz_0[j];

            tpz_yy_xyz_0[j] = pa_y[j] * tpz_y_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_y_xz_0[j];

            tpx_yy_xzz_0[j] = pa_y[j] * tpx_y_xzz_0[j] + 0.5 * fl1_fx * tpx_0_xzz_0[j];

            tpy_yy_xzz_0[j] = pa_y[j] * tpy_y_xzz_0[j] + 0.5 * fl1_fx * tpy_0_xzz_0[j] - fl1_fgb * fl1_fx * ts_y_xzz_0[j];

            tpz_yy_xzz_0[j] = pa_y[j] * tpz_y_xzz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j];

            tpx_yy_yyy_0[j] = pa_y[j] * tpx_y_yyy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j] + 1.5 * fl1_fx * tpx_y_yy_0[j];

            tpy_yy_yyy_0[j] =
                pa_y[j] * tpy_y_yyy_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j] + 1.5 * fl1_fx * tpy_y_yy_0[j] - fl1_fgb * fl1_fx * ts_y_yyy_0[j];

            tpz_yy_yyy_0[j] = pa_y[j] * tpz_y_yyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j] + 1.5 * fl1_fx * tpz_y_yy_0[j];

            tpx_yy_yyz_0[j] = pa_y[j] * tpx_y_yyz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] + fl1_fx * tpx_y_yz_0[j];

            tpy_yy_yyz_0[j] = pa_y[j] * tpy_y_yyz_0[j] + 0.5 * fl1_fx * tpy_0_yyz_0[j] + fl1_fx * tpy_y_yz_0[j] - fl1_fgb * fl1_fx * ts_y_yyz_0[j];

            tpz_yy_yyz_0[j] = pa_y[j] * tpz_y_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j] + fl1_fx * tpz_y_yz_0[j];

            tpx_yy_yzz_0[j] = pa_y[j] * tpx_y_yzz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] + 0.5 * fl1_fx * tpx_y_zz_0[j];

            tpy_yy_yzz_0[j] =
                pa_y[j] * tpy_y_yzz_0[j] + 0.5 * fl1_fx * tpy_0_yzz_0[j] + 0.5 * fl1_fx * tpy_y_zz_0[j] - fl1_fgb * fl1_fx * ts_y_yzz_0[j];

            tpz_yy_yzz_0[j] = pa_y[j] * tpz_y_yzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j] + 0.5 * fl1_fx * tpz_y_zz_0[j];

            tpx_yy_zzz_0[j] = pa_y[j] * tpx_y_zzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j];

            tpy_yy_zzz_0[j] = pa_y[j] * tpy_y_zzz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_y_zzz_0[j];

            tpz_yy_zzz_0[j] = pa_y[j] * tpz_y_zzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j];

            tpx_yz_xxx_0[j] = pa_y[j] * tpx_z_xxx_0[j];

            tpy_yz_xxx_0[j] = pa_y[j] * tpy_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxx_0[j];

            tpz_yz_xxx_0[j] = pa_y[j] * tpz_z_xxx_0[j];

            tpx_yz_xxy_0[j] = pa_y[j] * tpx_z_xxy_0[j] + 0.5 * fl1_fx * tpx_z_xx_0[j];

            tpy_yz_xxy_0[j] = pa_y[j] * tpy_z_xxy_0[j] + 0.5 * fl1_fx * tpy_z_xx_0[j] - fl1_fgb * fl1_fx * ts_z_xxy_0[j];

            tpz_yz_xxy_0[j] = pa_y[j] * tpz_z_xxy_0[j] + 0.5 * fl1_fx * tpz_z_xx_0[j];

            tpx_yz_xxz_0[j] = pa_y[j] * tpx_z_xxz_0[j];

            tpy_yz_xxz_0[j] = pa_y[j] * tpy_z_xxz_0[j] - fl1_fgb * fl1_fx * ts_z_xxz_0[j];

            tpz_yz_xxz_0[j] = pa_y[j] * tpz_z_xxz_0[j];

            tpx_yz_xyy_0[j] = pa_y[j] * tpx_z_xyy_0[j] + fl1_fx * tpx_z_xy_0[j];

            tpy_yz_xyy_0[j] = pa_y[j] * tpy_z_xyy_0[j] + fl1_fx * tpy_z_xy_0[j] - fl1_fgb * fl1_fx * ts_z_xyy_0[j];

            tpz_yz_xyy_0[j] = pa_y[j] * tpz_z_xyy_0[j] + fl1_fx * tpz_z_xy_0[j];

            tpx_yz_xyz_0[j] = pa_y[j] * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_z_xz_0[j];

            tpy_yz_xyz_0[j] = pa_y[j] * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_z_xz_0[j] - fl1_fgb * fl1_fx * ts_z_xyz_0[j];

            tpz_yz_xyz_0[j] = pa_y[j] * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_z_xz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDF_135_180(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (135,180)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

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

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

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

        // set up pointers to integrals

        auto tpx_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 45);

        auto tpy_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 46);

        auto tpy_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 47);

        auto tpy_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 48);

        auto tpy_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 49);

        auto tpy_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 53);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 54);

        auto tpy_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 55);

        auto tpy_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 56);

        auto tpy_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 57);

        auto tpy_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 58);

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 59);

        auto tpy_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 59);

        // Batch of Integrals (135,180)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, tpx_0_xyy_0, \
                                     tpx_0_xyz_0, tpx_0_xzz_0, tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yzz_0, tpx_0_zzz_0, \
                                     tpx_yz_xzz_0, tpx_yz_yyy_0, tpx_yz_yyz_0, tpx_yz_yzz_0, tpx_yz_zzz_0, tpx_z_xx_0, \
                                     tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xy_0, tpx_z_xyy_0, tpx_z_xyz_0, \
                                     tpx_z_xz_0, tpx_z_xzz_0, tpx_z_yy_0, tpx_z_yyy_0, tpx_z_yyz_0, tpx_z_yz_0, \
                                     tpx_z_yzz_0, tpx_z_zz_0, tpx_z_zzz_0, tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, \
                                     tpx_zz_xyy_0, tpx_zz_xyz_0, tpx_zz_xzz_0, tpx_zz_yyy_0, tpx_zz_yyz_0, tpx_zz_yzz_0, \
                                     tpx_zz_zzz_0, tpy_0_xxx_0, tpy_0_xxy_0, tpy_0_xxz_0, tpy_0_xyy_0, tpy_0_xyz_0, \
                                     tpy_0_xzz_0, tpy_0_yyy_0, tpy_0_yyz_0, tpy_0_yzz_0, tpy_0_zzz_0, tpy_yz_xzz_0, \
                                     tpy_yz_yyy_0, tpy_yz_yyz_0, tpy_yz_yzz_0, tpy_yz_zzz_0, tpy_z_xx_0, tpy_z_xxx_0, \
                                     tpy_z_xxy_0, tpy_z_xxz_0, tpy_z_xy_0, tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xz_0, \
                                     tpy_z_xzz_0, tpy_z_yy_0, tpy_z_yyy_0, tpy_z_yyz_0, tpy_z_yz_0, tpy_z_yzz_0, \
                                     tpy_z_zz_0, tpy_z_zzz_0, tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xyy_0, \
                                     tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, tpy_zz_yyz_0, tpy_zz_yzz_0, tpy_zz_zzz_0, \
                                     tpz_0_xxx_0, tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xzz_0, \
                                     tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yzz_0, tpz_0_zzz_0, tpz_yz_xzz_0, tpz_yz_yyy_0, \
                                     tpz_yz_yyz_0, tpz_yz_yzz_0, tpz_yz_zzz_0, tpz_z_xx_0, tpz_z_xxx_0, tpz_z_xxy_0, \
                                     tpz_z_xxz_0, tpz_z_xy_0, tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xz_0, tpz_z_xzz_0, \
                                     tpz_z_yy_0, tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yz_0, tpz_z_yzz_0, tpz_z_zz_0, \
                                     tpz_z_zzz_0, tpz_zz_xxx_0, tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xyy_0, tpz_zz_xyz_0, \
                                     tpz_zz_xzz_0, tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yzz_0, tpz_zz_zzz_0, ts_z_xxx_0, \
                                     ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, \
                                     ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yz_xzz_0[j] = pa_y[j] * tpx_z_xzz_0[j];

            tpy_yz_xzz_0[j] = pa_y[j] * tpy_z_xzz_0[j] - fl1_fgb * fl1_fx * ts_z_xzz_0[j];

            tpz_yz_xzz_0[j] = pa_y[j] * tpz_z_xzz_0[j];

            tpx_yz_yyy_0[j] = pa_y[j] * tpx_z_yyy_0[j] + 1.5 * fl1_fx * tpx_z_yy_0[j];

            tpy_yz_yyy_0[j] = pa_y[j] * tpy_z_yyy_0[j] + 1.5 * fl1_fx * tpy_z_yy_0[j] - fl1_fgb * fl1_fx * ts_z_yyy_0[j];

            tpz_yz_yyy_0[j] = pa_y[j] * tpz_z_yyy_0[j] + 1.5 * fl1_fx * tpz_z_yy_0[j];

            tpx_yz_yyz_0[j] = pa_y[j] * tpx_z_yyz_0[j] + fl1_fx * tpx_z_yz_0[j];

            tpy_yz_yyz_0[j] = pa_y[j] * tpy_z_yyz_0[j] + fl1_fx * tpy_z_yz_0[j] - fl1_fgb * fl1_fx * ts_z_yyz_0[j];

            tpz_yz_yyz_0[j] = pa_y[j] * tpz_z_yyz_0[j] + fl1_fx * tpz_z_yz_0[j];

            tpx_yz_yzz_0[j] = pa_y[j] * tpx_z_yzz_0[j] + 0.5 * fl1_fx * tpx_z_zz_0[j];

            tpy_yz_yzz_0[j] = pa_y[j] * tpy_z_yzz_0[j] + 0.5 * fl1_fx * tpy_z_zz_0[j] - fl1_fgb * fl1_fx * ts_z_yzz_0[j];

            tpz_yz_yzz_0[j] = pa_y[j] * tpz_z_yzz_0[j] + 0.5 * fl1_fx * tpz_z_zz_0[j];

            tpx_yz_zzz_0[j] = pa_y[j] * tpx_z_zzz_0[j];

            tpy_yz_zzz_0[j] = pa_y[j] * tpy_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_z_zzz_0[j];

            tpz_yz_zzz_0[j] = pa_y[j] * tpz_z_zzz_0[j];

            tpx_zz_xxx_0[j] = pa_z[j] * tpx_z_xxx_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j];

            tpy_zz_xxx_0[j] = pa_z[j] * tpy_z_xxx_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j];

            tpz_zz_xxx_0[j] = pa_z[j] * tpz_z_xxx_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxx_0[j];

            tpx_zz_xxy_0[j] = pa_z[j] * tpx_z_xxy_0[j] + 0.5 * fl1_fx * tpx_0_xxy_0[j];

            tpy_zz_xxy_0[j] = pa_z[j] * tpy_z_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xxy_0[j];

            tpz_zz_xxy_0[j] = pa_z[j] * tpz_z_xxy_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] - fl1_fgb * fl1_fx * ts_z_xxy_0[j];

            tpx_zz_xxz_0[j] = pa_z[j] * tpx_z_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xxz_0[j] + 0.5 * fl1_fx * tpx_z_xx_0[j];

            tpy_zz_xxz_0[j] = pa_z[j] * tpy_z_xxz_0[j] + 0.5 * fl1_fx * tpy_0_xxz_0[j] + 0.5 * fl1_fx * tpy_z_xx_0[j];

            tpz_zz_xxz_0[j] =
                pa_z[j] * tpz_z_xxz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j] + 0.5 * fl1_fx * tpz_z_xx_0[j] - fl1_fgb * fl1_fx * ts_z_xxz_0[j];

            tpx_zz_xyy_0[j] = pa_z[j] * tpx_z_xyy_0[j] + 0.5 * fl1_fx * tpx_0_xyy_0[j];

            tpy_zz_xyy_0[j] = pa_z[j] * tpy_z_xyy_0[j] + 0.5 * fl1_fx * tpy_0_xyy_0[j];

            tpz_zz_xyy_0[j] = pa_z[j] * tpz_z_xyy_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] - fl1_fgb * fl1_fx * ts_z_xyy_0[j];

            tpx_zz_xyz_0[j] = pa_z[j] * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_z_xy_0[j];

            tpy_zz_xyz_0[j] = pa_z[j] * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_z_xy_0[j];

            tpz_zz_xyz_0[j] =
                pa_z[j] * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_z_xy_0[j] - fl1_fgb * fl1_fx * ts_z_xyz_0[j];

            tpx_zz_xzz_0[j] = pa_z[j] * tpx_z_xzz_0[j] + 0.5 * fl1_fx * tpx_0_xzz_0[j] + fl1_fx * tpx_z_xz_0[j];

            tpy_zz_xzz_0[j] = pa_z[j] * tpy_z_xzz_0[j] + 0.5 * fl1_fx * tpy_0_xzz_0[j] + fl1_fx * tpy_z_xz_0[j];

            tpz_zz_xzz_0[j] = pa_z[j] * tpz_z_xzz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j] + fl1_fx * tpz_z_xz_0[j] - fl1_fgb * fl1_fx * ts_z_xzz_0[j];

            tpx_zz_yyy_0[j] = pa_z[j] * tpx_z_yyy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j];

            tpy_zz_yyy_0[j] = pa_z[j] * tpy_z_yyy_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j];

            tpz_zz_yyy_0[j] = pa_z[j] * tpz_z_yyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyy_0[j];

            tpx_zz_yyz_0[j] = pa_z[j] * tpx_z_yyz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] + 0.5 * fl1_fx * tpx_z_yy_0[j];

            tpy_zz_yyz_0[j] = pa_z[j] * tpy_z_yyz_0[j] + 0.5 * fl1_fx * tpy_0_yyz_0[j] + 0.5 * fl1_fx * tpy_z_yy_0[j];

            tpz_zz_yyz_0[j] =
                pa_z[j] * tpz_z_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j] + 0.5 * fl1_fx * tpz_z_yy_0[j] - fl1_fgb * fl1_fx * ts_z_yyz_0[j];

            tpx_zz_yzz_0[j] = pa_z[j] * tpx_z_yzz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] + fl1_fx * tpx_z_yz_0[j];

            tpy_zz_yzz_0[j] = pa_z[j] * tpy_z_yzz_0[j] + 0.5 * fl1_fx * tpy_0_yzz_0[j] + fl1_fx * tpy_z_yz_0[j];

            tpz_zz_yzz_0[j] = pa_z[j] * tpz_z_yzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j] + fl1_fx * tpz_z_yz_0[j] - fl1_fgb * fl1_fx * ts_z_yzz_0[j];

            tpx_zz_zzz_0[j] = pa_z[j] * tpx_z_zzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j] + 1.5 * fl1_fx * tpx_z_zz_0[j];

            tpy_zz_zzz_0[j] = pa_z[j] * tpy_z_zzz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j] + 1.5 * fl1_fx * tpy_z_zz_0[j];

            tpz_zz_zzz_0[j] =
                pa_z[j] * tpz_z_zzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j] + 1.5 * fl1_fx * tpz_z_zz_0[j] - fl1_fgb * fl1_fx * ts_z_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFD(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForFD_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFD_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFD_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFD_135_180(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForFD_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx);

        auto tpy_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx);

        auto tpz_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx);

        auto tpx_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 1);

        auto tpy_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tpz_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tpx_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 2);

        auto tpy_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tpz_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tpx_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 3);

        auto tpy_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tpz_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tpx_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 4);

        auto tpy_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tpz_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tpx_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 5);

        auto tpy_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tpz_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tpx_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 6);

        auto tpy_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tpz_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tpx_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 7);

        auto tpy_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tpz_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tpx_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 8);

        auto tpy_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tpz_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tpx_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 9);

        auto tpy_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tpz_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tpx_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 10);

        auto tpy_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tpz_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tpx_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 11);

        auto tpy_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tpz_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 11);

        auto tpx_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 12);

        auto tpy_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tpz_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tpx_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 13);

        auto tpy_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tpz_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tpx_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 14);

        auto tpy_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tpz_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tpx_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx);

        auto tpy_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx);

        auto tpz_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx);

        auto tpx_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 1);

        auto tpy_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 2);

        auto tpy_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 3);

        auto tpy_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 4);

        auto tpy_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 5);

        auto tpy_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx);

        auto tpy_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx);

        auto tpz_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx);

        auto tpx_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 1);

        auto tpy_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 2);

        auto tpy_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 3);

        auto tpy_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 4);

        auto tpy_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 5);

        auto tpy_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 6);

        auto tpy_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 7);

        auto tpy_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 8);

        auto tpy_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto ts_xx_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx);

        auto ts_xx_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 1);

        auto ts_xx_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 2);

        auto ts_xx_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 3);

        auto ts_xx_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 4);

        auto ts_xx_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 5);

        auto ts_xy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 6);

        auto ts_xy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 7);

        auto ts_xy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 8);

        auto ts_xy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 9);

        auto ts_xy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 10);

        auto ts_xy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 11);

        auto ts_xz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 12);

        auto ts_xz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 13);

        auto ts_xz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 14);

        // set up pointers to integrals

        auto tpx_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx);

        auto tpy_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx);

        auto tpz_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx);

        auto tpx_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 1);

        auto tpy_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 1);

        auto tpx_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 2);

        auto tpy_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 2);

        auto tpx_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 3);

        auto tpy_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 3);

        auto tpx_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 4);

        auto tpy_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 4);

        auto tpx_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 5);

        auto tpy_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 5);

        auto tpx_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 6);

        auto tpy_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 6);

        auto tpx_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 7);

        auto tpy_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 7);

        auto tpx_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 8);

        auto tpy_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 8);

        auto tpx_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 9);

        auto tpy_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 9);

        auto tpx_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 10);

        auto tpy_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 10);

        auto tpx_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 11);

        auto tpy_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 11);

        auto tpz_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 11);

        auto tpx_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 12);

        auto tpy_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 12);

        auto tpx_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 13);

        auto tpy_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 13);

        auto tpx_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 14);

        auto tpy_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_x_xx_0, tpx_x_xy_0, tpx_x_xz_0, tpx_x_yy_0, tpx_x_yz_0, \
                                     tpx_x_zz_0, tpx_xx_x_0, tpx_xx_xx_0, tpx_xx_xy_0, tpx_xx_xz_0, tpx_xx_y_0, \
                                     tpx_xx_yy_0, tpx_xx_yz_0, tpx_xx_z_0, tpx_xx_zz_0, tpx_xxx_xx_0, tpx_xxx_xy_0, \
                                     tpx_xxx_xz_0, tpx_xxx_yy_0, tpx_xxx_yz_0, tpx_xxx_zz_0, tpx_xxy_xx_0, tpx_xxy_xy_0, \
                                     tpx_xxy_xz_0, tpx_xxy_yy_0, tpx_xxy_yz_0, tpx_xxy_zz_0, tpx_xxz_xx_0, tpx_xxz_xy_0, \
                                     tpx_xxz_xz_0, tpx_xy_x_0, tpx_xy_xx_0, tpx_xy_xy_0, tpx_xy_xz_0, tpx_xy_y_0, \
                                     tpx_xy_yy_0, tpx_xy_yz_0, tpx_xy_z_0, tpx_xy_zz_0, tpx_xz_x_0, tpx_xz_xx_0, \
                                     tpx_xz_xy_0, tpx_xz_xz_0, tpx_xz_y_0, tpx_xz_z_0, tpx_y_xx_0, tpx_y_xy_0, \
                                     tpx_y_xz_0, tpx_y_yy_0, tpx_y_yz_0, tpx_y_zz_0, tpx_z_xx_0, tpx_z_xy_0, tpx_z_xz_0, \
                                     tpy_x_xx_0, tpy_x_xy_0, tpy_x_xz_0, tpy_x_yy_0, tpy_x_yz_0, tpy_x_zz_0, tpy_xx_x_0, \
                                     tpy_xx_xx_0, tpy_xx_xy_0, tpy_xx_xz_0, tpy_xx_y_0, tpy_xx_yy_0, tpy_xx_yz_0, \
                                     tpy_xx_z_0, tpy_xx_zz_0, tpy_xxx_xx_0, tpy_xxx_xy_0, tpy_xxx_xz_0, tpy_xxx_yy_0, \
                                     tpy_xxx_yz_0, tpy_xxx_zz_0, tpy_xxy_xx_0, tpy_xxy_xy_0, tpy_xxy_xz_0, tpy_xxy_yy_0, \
                                     tpy_xxy_yz_0, tpy_xxy_zz_0, tpy_xxz_xx_0, tpy_xxz_xy_0, tpy_xxz_xz_0, tpy_xy_x_0, \
                                     tpy_xy_xx_0, tpy_xy_xy_0, tpy_xy_xz_0, tpy_xy_y_0, tpy_xy_yy_0, tpy_xy_yz_0, \
                                     tpy_xy_z_0, tpy_xy_zz_0, tpy_xz_x_0, tpy_xz_xx_0, tpy_xz_xy_0, tpy_xz_xz_0, \
                                     tpy_xz_y_0, tpy_xz_z_0, tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpy_y_yy_0, tpy_y_yz_0, \
                                     tpy_y_zz_0, tpy_z_xx_0, tpy_z_xy_0, tpy_z_xz_0, tpz_x_xx_0, tpz_x_xy_0, tpz_x_xz_0, \
                                     tpz_x_yy_0, tpz_x_yz_0, tpz_x_zz_0, tpz_xx_x_0, tpz_xx_xx_0, tpz_xx_xy_0, \
                                     tpz_xx_xz_0, tpz_xx_y_0, tpz_xx_yy_0, tpz_xx_yz_0, tpz_xx_z_0, tpz_xx_zz_0, \
                                     tpz_xxx_xx_0, tpz_xxx_xy_0, tpz_xxx_xz_0, tpz_xxx_yy_0, tpz_xxx_yz_0, tpz_xxx_zz_0, \
                                     tpz_xxy_xx_0, tpz_xxy_xy_0, tpz_xxy_xz_0, tpz_xxy_yy_0, tpz_xxy_yz_0, tpz_xxy_zz_0, \
                                     tpz_xxz_xx_0, tpz_xxz_xy_0, tpz_xxz_xz_0, tpz_xy_x_0, tpz_xy_xx_0, tpz_xy_xy_0, \
                                     tpz_xy_xz_0, tpz_xy_y_0, tpz_xy_yy_0, tpz_xy_yz_0, tpz_xy_z_0, tpz_xy_zz_0, \
                                     tpz_xz_x_0, tpz_xz_xx_0, tpz_xz_xy_0, tpz_xz_xz_0, tpz_xz_y_0, tpz_xz_z_0, \
                                     tpz_y_xx_0, tpz_y_xy_0, tpz_y_xz_0, tpz_y_yy_0, tpz_y_yz_0, tpz_y_zz_0, tpz_z_xx_0, \
                                     tpz_z_xy_0, tpz_z_xz_0, ts_xx_xx_0, ts_xx_xy_0, ts_xx_xz_0, ts_xx_yy_0, ts_xx_yz_0, \
                                     ts_xx_zz_0, ts_xy_xx_0, ts_xy_xy_0, ts_xy_xz_0, ts_xy_yy_0, ts_xy_yz_0, ts_xy_zz_0, \
                                     ts_xz_xx_0, ts_xz_xy_0, ts_xz_xz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxx_xx_0[j] = pa_x[j] * tpx_xx_xx_0[j] + fl1_fx * tpx_x_xx_0[j] + fl1_fx * tpx_xx_x_0[j] - fl1_fgb * fl1_fx * ts_xx_xx_0[j];

            tpy_xxx_xx_0[j] = pa_x[j] * tpy_xx_xx_0[j] + fl1_fx * tpy_x_xx_0[j] + fl1_fx * tpy_xx_x_0[j];

            tpz_xxx_xx_0[j] = pa_x[j] * tpz_xx_xx_0[j] + fl1_fx * tpz_x_xx_0[j] + fl1_fx * tpz_xx_x_0[j];

            tpx_xxx_xy_0[j] = pa_x[j] * tpx_xx_xy_0[j] + fl1_fx * tpx_x_xy_0[j] + 0.5 * fl1_fx * tpx_xx_y_0[j] - fl1_fgb * fl1_fx * ts_xx_xy_0[j];

            tpy_xxx_xy_0[j] = pa_x[j] * tpy_xx_xy_0[j] + fl1_fx * tpy_x_xy_0[j] + 0.5 * fl1_fx * tpy_xx_y_0[j];

            tpz_xxx_xy_0[j] = pa_x[j] * tpz_xx_xy_0[j] + fl1_fx * tpz_x_xy_0[j] + 0.5 * fl1_fx * tpz_xx_y_0[j];

            tpx_xxx_xz_0[j] = pa_x[j] * tpx_xx_xz_0[j] + fl1_fx * tpx_x_xz_0[j] + 0.5 * fl1_fx * tpx_xx_z_0[j] - fl1_fgb * fl1_fx * ts_xx_xz_0[j];

            tpy_xxx_xz_0[j] = pa_x[j] * tpy_xx_xz_0[j] + fl1_fx * tpy_x_xz_0[j] + 0.5 * fl1_fx * tpy_xx_z_0[j];

            tpz_xxx_xz_0[j] = pa_x[j] * tpz_xx_xz_0[j] + fl1_fx * tpz_x_xz_0[j] + 0.5 * fl1_fx * tpz_xx_z_0[j];

            tpx_xxx_yy_0[j] = pa_x[j] * tpx_xx_yy_0[j] + fl1_fx * tpx_x_yy_0[j] - fl1_fgb * fl1_fx * ts_xx_yy_0[j];

            tpy_xxx_yy_0[j] = pa_x[j] * tpy_xx_yy_0[j] + fl1_fx * tpy_x_yy_0[j];

            tpz_xxx_yy_0[j] = pa_x[j] * tpz_xx_yy_0[j] + fl1_fx * tpz_x_yy_0[j];

            tpx_xxx_yz_0[j] = pa_x[j] * tpx_xx_yz_0[j] + fl1_fx * tpx_x_yz_0[j] - fl1_fgb * fl1_fx * ts_xx_yz_0[j];

            tpy_xxx_yz_0[j] = pa_x[j] * tpy_xx_yz_0[j] + fl1_fx * tpy_x_yz_0[j];

            tpz_xxx_yz_0[j] = pa_x[j] * tpz_xx_yz_0[j] + fl1_fx * tpz_x_yz_0[j];

            tpx_xxx_zz_0[j] = pa_x[j] * tpx_xx_zz_0[j] + fl1_fx * tpx_x_zz_0[j] - fl1_fgb * fl1_fx * ts_xx_zz_0[j];

            tpy_xxx_zz_0[j] = pa_x[j] * tpy_xx_zz_0[j] + fl1_fx * tpy_x_zz_0[j];

            tpz_xxx_zz_0[j] = pa_x[j] * tpz_xx_zz_0[j] + fl1_fx * tpz_x_zz_0[j];

            tpx_xxy_xx_0[j] = pa_x[j] * tpx_xy_xx_0[j] + 0.5 * fl1_fx * tpx_y_xx_0[j] + fl1_fx * tpx_xy_x_0[j] - fl1_fgb * fl1_fx * ts_xy_xx_0[j];

            tpy_xxy_xx_0[j] = pa_x[j] * tpy_xy_xx_0[j] + 0.5 * fl1_fx * tpy_y_xx_0[j] + fl1_fx * tpy_xy_x_0[j];

            tpz_xxy_xx_0[j] = pa_x[j] * tpz_xy_xx_0[j] + 0.5 * fl1_fx * tpz_y_xx_0[j] + fl1_fx * tpz_xy_x_0[j];

            tpx_xxy_xy_0[j] =
                pa_x[j] * tpx_xy_xy_0[j] + 0.5 * fl1_fx * tpx_y_xy_0[j] + 0.5 * fl1_fx * tpx_xy_y_0[j] - fl1_fgb * fl1_fx * ts_xy_xy_0[j];

            tpy_xxy_xy_0[j] = pa_x[j] * tpy_xy_xy_0[j] + 0.5 * fl1_fx * tpy_y_xy_0[j] + 0.5 * fl1_fx * tpy_xy_y_0[j];

            tpz_xxy_xy_0[j] = pa_x[j] * tpz_xy_xy_0[j] + 0.5 * fl1_fx * tpz_y_xy_0[j] + 0.5 * fl1_fx * tpz_xy_y_0[j];

            tpx_xxy_xz_0[j] =
                pa_x[j] * tpx_xy_xz_0[j] + 0.5 * fl1_fx * tpx_y_xz_0[j] + 0.5 * fl1_fx * tpx_xy_z_0[j] - fl1_fgb * fl1_fx * ts_xy_xz_0[j];

            tpy_xxy_xz_0[j] = pa_x[j] * tpy_xy_xz_0[j] + 0.5 * fl1_fx * tpy_y_xz_0[j] + 0.5 * fl1_fx * tpy_xy_z_0[j];

            tpz_xxy_xz_0[j] = pa_x[j] * tpz_xy_xz_0[j] + 0.5 * fl1_fx * tpz_y_xz_0[j] + 0.5 * fl1_fx * tpz_xy_z_0[j];

            tpx_xxy_yy_0[j] = pa_x[j] * tpx_xy_yy_0[j] + 0.5 * fl1_fx * tpx_y_yy_0[j] - fl1_fgb * fl1_fx * ts_xy_yy_0[j];

            tpy_xxy_yy_0[j] = pa_x[j] * tpy_xy_yy_0[j] + 0.5 * fl1_fx * tpy_y_yy_0[j];

            tpz_xxy_yy_0[j] = pa_x[j] * tpz_xy_yy_0[j] + 0.5 * fl1_fx * tpz_y_yy_0[j];

            tpx_xxy_yz_0[j] = pa_x[j] * tpx_xy_yz_0[j] + 0.5 * fl1_fx * tpx_y_yz_0[j] - fl1_fgb * fl1_fx * ts_xy_yz_0[j];

            tpy_xxy_yz_0[j] = pa_x[j] * tpy_xy_yz_0[j] + 0.5 * fl1_fx * tpy_y_yz_0[j];

            tpz_xxy_yz_0[j] = pa_x[j] * tpz_xy_yz_0[j] + 0.5 * fl1_fx * tpz_y_yz_0[j];

            tpx_xxy_zz_0[j] = pa_x[j] * tpx_xy_zz_0[j] + 0.5 * fl1_fx * tpx_y_zz_0[j] - fl1_fgb * fl1_fx * ts_xy_zz_0[j];

            tpy_xxy_zz_0[j] = pa_x[j] * tpy_xy_zz_0[j] + 0.5 * fl1_fx * tpy_y_zz_0[j];

            tpz_xxy_zz_0[j] = pa_x[j] * tpz_xy_zz_0[j] + 0.5 * fl1_fx * tpz_y_zz_0[j];

            tpx_xxz_xx_0[j] = pa_x[j] * tpx_xz_xx_0[j] + 0.5 * fl1_fx * tpx_z_xx_0[j] + fl1_fx * tpx_xz_x_0[j] - fl1_fgb * fl1_fx * ts_xz_xx_0[j];

            tpy_xxz_xx_0[j] = pa_x[j] * tpy_xz_xx_0[j] + 0.5 * fl1_fx * tpy_z_xx_0[j] + fl1_fx * tpy_xz_x_0[j];

            tpz_xxz_xx_0[j] = pa_x[j] * tpz_xz_xx_0[j] + 0.5 * fl1_fx * tpz_z_xx_0[j] + fl1_fx * tpz_xz_x_0[j];

            tpx_xxz_xy_0[j] =
                pa_x[j] * tpx_xz_xy_0[j] + 0.5 * fl1_fx * tpx_z_xy_0[j] + 0.5 * fl1_fx * tpx_xz_y_0[j] - fl1_fgb * fl1_fx * ts_xz_xy_0[j];

            tpy_xxz_xy_0[j] = pa_x[j] * tpy_xz_xy_0[j] + 0.5 * fl1_fx * tpy_z_xy_0[j] + 0.5 * fl1_fx * tpy_xz_y_0[j];

            tpz_xxz_xy_0[j] = pa_x[j] * tpz_xz_xy_0[j] + 0.5 * fl1_fx * tpz_z_xy_0[j] + 0.5 * fl1_fx * tpz_xz_y_0[j];

            tpx_xxz_xz_0[j] =
                pa_x[j] * tpx_xz_xz_0[j] + 0.5 * fl1_fx * tpx_z_xz_0[j] + 0.5 * fl1_fx * tpx_xz_z_0[j] - fl1_fgb * fl1_fx * ts_xz_xz_0[j];

            tpy_xxz_xz_0[j] = pa_x[j] * tpy_xz_xz_0[j] + 0.5 * fl1_fx * tpy_z_xz_0[j] + 0.5 * fl1_fx * tpy_xz_z_0[j];

            tpz_xxz_xz_0[j] = pa_x[j] * tpz_xz_xz_0[j] + 0.5 * fl1_fx * tpz_z_xz_0[j] + 0.5 * fl1_fx * tpz_xz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFD_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 15);

        auto tpy_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tpz_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tpx_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 16);

        auto tpy_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tpz_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tpx_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 17);

        auto tpy_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tpz_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto ts_xz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 15);

        auto ts_xz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 16);

        auto ts_xz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 17);

        auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18);

        auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19);

        auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20);

        auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21);

        auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22);

        auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23);

        auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24);

        auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25);

        auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26);

        auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27);

        auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28);

        auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29);

        // set up pointers to integrals

        auto tpx_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 15);

        auto tpy_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 15);

        auto tpx_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 16);

        auto tpy_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 16);

        auto tpz_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 16);

        auto tpx_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 17);

        auto tpy_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 17);

        auto tpx_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 18);

        auto tpy_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 18);

        auto tpx_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 19);

        auto tpy_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 19);

        auto tpx_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 20);

        auto tpy_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 20);

        auto tpx_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 21);

        auto tpy_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 21);

        auto tpx_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 22);

        auto tpy_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 22);

        auto tpx_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 23);

        auto tpy_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 23);

        auto tpx_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 24);

        auto tpy_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 24);

        auto tpx_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 25);

        auto tpy_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 25);

        auto tpx_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 26);

        auto tpy_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 26);

        auto tpx_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 27);

        auto tpy_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 27);

        auto tpx_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 28);

        auto tpy_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 28);

        auto tpx_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 29);

        auto tpy_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxz_yy_0, tpx_xxz_yz_0, tpx_xxz_zz_0, tpx_xyy_xx_0, \
                                     tpx_xyy_xy_0, tpx_xyy_xz_0, tpx_xyy_yy_0, tpx_xyy_yz_0, tpx_xyy_zz_0, tpx_xyz_xx_0, \
                                     tpx_xyz_xy_0, tpx_xyz_xz_0, tpx_xyz_yy_0, tpx_xyz_yz_0, tpx_xyz_zz_0, tpx_xz_yy_0, \
                                     tpx_xz_yz_0, tpx_xz_zz_0, tpx_yy_x_0, tpx_yy_xx_0, tpx_yy_xy_0, tpx_yy_xz_0, \
                                     tpx_yy_y_0, tpx_yy_yy_0, tpx_yy_yz_0, tpx_yy_z_0, tpx_yy_zz_0, tpx_yz_x_0, \
                                     tpx_yz_xx_0, tpx_yz_xy_0, tpx_yz_xz_0, tpx_yz_y_0, tpx_yz_yy_0, tpx_yz_yz_0, \
                                     tpx_yz_z_0, tpx_yz_zz_0, tpx_z_yy_0, tpx_z_yz_0, tpx_z_zz_0, tpy_xxz_yy_0, \
                                     tpy_xxz_yz_0, tpy_xxz_zz_0, tpy_xyy_xx_0, tpy_xyy_xy_0, tpy_xyy_xz_0, tpy_xyy_yy_0, \
                                     tpy_xyy_yz_0, tpy_xyy_zz_0, tpy_xyz_xx_0, tpy_xyz_xy_0, tpy_xyz_xz_0, tpy_xyz_yy_0, \
                                     tpy_xyz_yz_0, tpy_xyz_zz_0, tpy_xz_yy_0, tpy_xz_yz_0, tpy_xz_zz_0, tpy_yy_x_0, \
                                     tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_yy_y_0, tpy_yy_yy_0, tpy_yy_yz_0, \
                                     tpy_yy_z_0, tpy_yy_zz_0, tpy_yz_x_0, tpy_yz_xx_0, tpy_yz_xy_0, tpy_yz_xz_0, \
                                     tpy_yz_y_0, tpy_yz_yy_0, tpy_yz_yz_0, tpy_yz_z_0, tpy_yz_zz_0, tpy_z_yy_0, \
                                     tpy_z_yz_0, tpy_z_zz_0, tpz_xxz_yy_0, tpz_xxz_yz_0, tpz_xxz_zz_0, tpz_xyy_xx_0, \
                                     tpz_xyy_xy_0, tpz_xyy_xz_0, tpz_xyy_yy_0, tpz_xyy_yz_0, tpz_xyy_zz_0, tpz_xyz_xx_0, \
                                     tpz_xyz_xy_0, tpz_xyz_xz_0, tpz_xyz_yy_0, tpz_xyz_yz_0, tpz_xyz_zz_0, tpz_xz_yy_0, \
                                     tpz_xz_yz_0, tpz_xz_zz_0, tpz_yy_x_0, tpz_yy_xx_0, tpz_yy_xy_0, tpz_yy_xz_0, \
                                     tpz_yy_y_0, tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_z_0, tpz_yy_zz_0, tpz_yz_x_0, \
                                     tpz_yz_xx_0, tpz_yz_xy_0, tpz_yz_xz_0, tpz_yz_y_0, tpz_yz_yy_0, tpz_yz_yz_0, \
                                     tpz_yz_z_0, tpz_yz_zz_0, tpz_z_yy_0, tpz_z_yz_0, tpz_z_zz_0, ts_xz_yy_0, \
                                     ts_xz_yz_0, ts_xz_zz_0, ts_yy_xx_0, ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, \
                                     ts_yy_zz_0, ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, ts_yz_yy_0, ts_yz_yz_0, ts_yz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxz_yy_0[j] = pa_x[j] * tpx_xz_yy_0[j] + 0.5 * fl1_fx * tpx_z_yy_0[j] - fl1_fgb * fl1_fx * ts_xz_yy_0[j];

            tpy_xxz_yy_0[j] = pa_x[j] * tpy_xz_yy_0[j] + 0.5 * fl1_fx * tpy_z_yy_0[j];

            tpz_xxz_yy_0[j] = pa_x[j] * tpz_xz_yy_0[j] + 0.5 * fl1_fx * tpz_z_yy_0[j];

            tpx_xxz_yz_0[j] = pa_x[j] * tpx_xz_yz_0[j] + 0.5 * fl1_fx * tpx_z_yz_0[j] - fl1_fgb * fl1_fx * ts_xz_yz_0[j];

            tpy_xxz_yz_0[j] = pa_x[j] * tpy_xz_yz_0[j] + 0.5 * fl1_fx * tpy_z_yz_0[j];

            tpz_xxz_yz_0[j] = pa_x[j] * tpz_xz_yz_0[j] + 0.5 * fl1_fx * tpz_z_yz_0[j];

            tpx_xxz_zz_0[j] = pa_x[j] * tpx_xz_zz_0[j] + 0.5 * fl1_fx * tpx_z_zz_0[j] - fl1_fgb * fl1_fx * ts_xz_zz_0[j];

            tpy_xxz_zz_0[j] = pa_x[j] * tpy_xz_zz_0[j] + 0.5 * fl1_fx * tpy_z_zz_0[j];

            tpz_xxz_zz_0[j] = pa_x[j] * tpz_xz_zz_0[j] + 0.5 * fl1_fx * tpz_z_zz_0[j];

            tpx_xyy_xx_0[j] = pa_x[j] * tpx_yy_xx_0[j] + fl1_fx * tpx_yy_x_0[j] - fl1_fgb * fl1_fx * ts_yy_xx_0[j];

            tpy_xyy_xx_0[j] = pa_x[j] * tpy_yy_xx_0[j] + fl1_fx * tpy_yy_x_0[j];

            tpz_xyy_xx_0[j] = pa_x[j] * tpz_yy_xx_0[j] + fl1_fx * tpz_yy_x_0[j];

            tpx_xyy_xy_0[j] = pa_x[j] * tpx_yy_xy_0[j] + 0.5 * fl1_fx * tpx_yy_y_0[j] - fl1_fgb * fl1_fx * ts_yy_xy_0[j];

            tpy_xyy_xy_0[j] = pa_x[j] * tpy_yy_xy_0[j] + 0.5 * fl1_fx * tpy_yy_y_0[j];

            tpz_xyy_xy_0[j] = pa_x[j] * tpz_yy_xy_0[j] + 0.5 * fl1_fx * tpz_yy_y_0[j];

            tpx_xyy_xz_0[j] = pa_x[j] * tpx_yy_xz_0[j] + 0.5 * fl1_fx * tpx_yy_z_0[j] - fl1_fgb * fl1_fx * ts_yy_xz_0[j];

            tpy_xyy_xz_0[j] = pa_x[j] * tpy_yy_xz_0[j] + 0.5 * fl1_fx * tpy_yy_z_0[j];

            tpz_xyy_xz_0[j] = pa_x[j] * tpz_yy_xz_0[j] + 0.5 * fl1_fx * tpz_yy_z_0[j];

            tpx_xyy_yy_0[j] = pa_x[j] * tpx_yy_yy_0[j] - fl1_fgb * fl1_fx * ts_yy_yy_0[j];

            tpy_xyy_yy_0[j] = pa_x[j] * tpy_yy_yy_0[j];

            tpz_xyy_yy_0[j] = pa_x[j] * tpz_yy_yy_0[j];

            tpx_xyy_yz_0[j] = pa_x[j] * tpx_yy_yz_0[j] - fl1_fgb * fl1_fx * ts_yy_yz_0[j];

            tpy_xyy_yz_0[j] = pa_x[j] * tpy_yy_yz_0[j];

            tpz_xyy_yz_0[j] = pa_x[j] * tpz_yy_yz_0[j];

            tpx_xyy_zz_0[j] = pa_x[j] * tpx_yy_zz_0[j] - fl1_fgb * fl1_fx * ts_yy_zz_0[j];

            tpy_xyy_zz_0[j] = pa_x[j] * tpy_yy_zz_0[j];

            tpz_xyy_zz_0[j] = pa_x[j] * tpz_yy_zz_0[j];

            tpx_xyz_xx_0[j] = pa_x[j] * tpx_yz_xx_0[j] + fl1_fx * tpx_yz_x_0[j] - fl1_fgb * fl1_fx * ts_yz_xx_0[j];

            tpy_xyz_xx_0[j] = pa_x[j] * tpy_yz_xx_0[j] + fl1_fx * tpy_yz_x_0[j];

            tpz_xyz_xx_0[j] = pa_x[j] * tpz_yz_xx_0[j] + fl1_fx * tpz_yz_x_0[j];

            tpx_xyz_xy_0[j] = pa_x[j] * tpx_yz_xy_0[j] + 0.5 * fl1_fx * tpx_yz_y_0[j] - fl1_fgb * fl1_fx * ts_yz_xy_0[j];

            tpy_xyz_xy_0[j] = pa_x[j] * tpy_yz_xy_0[j] + 0.5 * fl1_fx * tpy_yz_y_0[j];

            tpz_xyz_xy_0[j] = pa_x[j] * tpz_yz_xy_0[j] + 0.5 * fl1_fx * tpz_yz_y_0[j];

            tpx_xyz_xz_0[j] = pa_x[j] * tpx_yz_xz_0[j] + 0.5 * fl1_fx * tpx_yz_z_0[j] - fl1_fgb * fl1_fx * ts_yz_xz_0[j];

            tpy_xyz_xz_0[j] = pa_x[j] * tpy_yz_xz_0[j] + 0.5 * fl1_fx * tpy_yz_z_0[j];

            tpz_xyz_xz_0[j] = pa_x[j] * tpz_yz_xz_0[j] + 0.5 * fl1_fx * tpz_yz_z_0[j];

            tpx_xyz_yy_0[j] = pa_x[j] * tpx_yz_yy_0[j] - fl1_fgb * fl1_fx * ts_yz_yy_0[j];

            tpy_xyz_yy_0[j] = pa_x[j] * tpy_yz_yy_0[j];

            tpz_xyz_yy_0[j] = pa_x[j] * tpz_yz_yy_0[j];

            tpx_xyz_yz_0[j] = pa_x[j] * tpx_yz_yz_0[j] - fl1_fgb * fl1_fx * ts_yz_yz_0[j];

            tpy_xyz_yz_0[j] = pa_x[j] * tpy_yz_yz_0[j];

            tpz_xyz_yz_0[j] = pa_x[j] * tpz_yz_yz_0[j];

            tpx_xyz_zz_0[j] = pa_x[j] * tpx_yz_zz_0[j] - fl1_fgb * fl1_fx * ts_yz_zz_0[j];

            tpy_xyz_zz_0[j] = pa_x[j] * tpy_yz_zz_0[j];

            tpz_xyz_zz_0[j] = pa_x[j] * tpz_yz_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFD_90_135(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18);

        auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19);

        auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20);

        auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21);

        auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22);

        auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23);

        auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24);

        auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25);

        auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26);

        auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30);

        auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31);

        auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32);

        auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33);

        auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34);

        auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35);

        // set up pointers to integrals

        auto tpx_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 30);

        auto tpy_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 31);

        auto tpy_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 32);

        auto tpy_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 33);

        auto tpy_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 34);

        auto tpy_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 35);

        auto tpy_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 35);

        auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36);

        auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37);

        auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38);

        auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39);

        auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40);

        auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41);

        auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42);

        auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43);

        auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44);

        auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_xzz_xx_0, tpx_xzz_xy_0, tpx_xzz_xz_0, \
                                     tpx_xzz_yy_0, tpx_xzz_yz_0, tpx_xzz_zz_0, tpx_y_xx_0, tpx_y_xy_0, tpx_y_xz_0, \
                                     tpx_y_yy_0, tpx_y_yz_0, tpx_y_zz_0, tpx_yy_x_0, tpx_yy_xx_0, tpx_yy_xy_0, \
                                     tpx_yy_xz_0, tpx_yy_y_0, tpx_yy_yy_0, tpx_yy_yz_0, tpx_yy_z_0, tpx_yy_zz_0, \
                                     tpx_yyy_xx_0, tpx_yyy_xy_0, tpx_yyy_xz_0, tpx_yyy_yy_0, tpx_yyy_yz_0, tpx_yyy_zz_0, \
                                     tpx_yyz_xx_0, tpx_yyz_xy_0, tpx_yyz_xz_0, tpx_yz_x_0, tpx_yz_xx_0, tpx_yz_xy_0, \
                                     tpx_yz_xz_0, tpx_z_xx_0, tpx_z_xy_0, tpx_z_xz_0, tpx_zz_x_0, tpx_zz_xx_0, \
                                     tpx_zz_xy_0, tpx_zz_xz_0, tpx_zz_y_0, tpx_zz_yy_0, tpx_zz_yz_0, tpx_zz_z_0, \
                                     tpx_zz_zz_0, tpy_xzz_xx_0, tpy_xzz_xy_0, tpy_xzz_xz_0, tpy_xzz_yy_0, tpy_xzz_yz_0, \
                                     tpy_xzz_zz_0, tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpy_y_yy_0, tpy_y_yz_0, tpy_y_zz_0, \
                                     tpy_yy_x_0, tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_yy_y_0, tpy_yy_yy_0, \
                                     tpy_yy_yz_0, tpy_yy_z_0, tpy_yy_zz_0, tpy_yyy_xx_0, tpy_yyy_xy_0, tpy_yyy_xz_0, \
                                     tpy_yyy_yy_0, tpy_yyy_yz_0, tpy_yyy_zz_0, tpy_yyz_xx_0, tpy_yyz_xy_0, tpy_yyz_xz_0, \
                                     tpy_yz_x_0, tpy_yz_xx_0, tpy_yz_xy_0, tpy_yz_xz_0, tpy_z_xx_0, tpy_z_xy_0, \
                                     tpy_z_xz_0, tpy_zz_x_0, tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpy_zz_y_0, \
                                     tpy_zz_yy_0, tpy_zz_yz_0, tpy_zz_z_0, tpy_zz_zz_0, tpz_xzz_xx_0, tpz_xzz_xy_0, \
                                     tpz_xzz_xz_0, tpz_xzz_yy_0, tpz_xzz_yz_0, tpz_xzz_zz_0, tpz_y_xx_0, tpz_y_xy_0, \
                                     tpz_y_xz_0, tpz_y_yy_0, tpz_y_yz_0, tpz_y_zz_0, tpz_yy_x_0, tpz_yy_xx_0, \
                                     tpz_yy_xy_0, tpz_yy_xz_0, tpz_yy_y_0, tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_z_0, \
                                     tpz_yy_zz_0, tpz_yyy_xx_0, tpz_yyy_xy_0, tpz_yyy_xz_0, tpz_yyy_yy_0, tpz_yyy_yz_0, \
                                     tpz_yyy_zz_0, tpz_yyz_xx_0, tpz_yyz_xy_0, tpz_yyz_xz_0, tpz_yz_x_0, tpz_yz_xx_0, \
                                     tpz_yz_xy_0, tpz_yz_xz_0, tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_zz_x_0, \
                                     tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_y_0, tpz_zz_yy_0, tpz_zz_yz_0, \
                                     tpz_zz_z_0, tpz_zz_zz_0, ts_yy_xx_0, ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, \
                                     ts_yy_yz_0, ts_yy_zz_0, ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, ts_zz_xx_0, ts_zz_xy_0, \
                                     ts_zz_xz_0, ts_zz_yy_0, ts_zz_yz_0, ts_zz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xzz_xx_0[j] = pa_x[j] * tpx_zz_xx_0[j] + fl1_fx * tpx_zz_x_0[j] - fl1_fgb * fl1_fx * ts_zz_xx_0[j];

            tpy_xzz_xx_0[j] = pa_x[j] * tpy_zz_xx_0[j] + fl1_fx * tpy_zz_x_0[j];

            tpz_xzz_xx_0[j] = pa_x[j] * tpz_zz_xx_0[j] + fl1_fx * tpz_zz_x_0[j];

            tpx_xzz_xy_0[j] = pa_x[j] * tpx_zz_xy_0[j] + 0.5 * fl1_fx * tpx_zz_y_0[j] - fl1_fgb * fl1_fx * ts_zz_xy_0[j];

            tpy_xzz_xy_0[j] = pa_x[j] * tpy_zz_xy_0[j] + 0.5 * fl1_fx * tpy_zz_y_0[j];

            tpz_xzz_xy_0[j] = pa_x[j] * tpz_zz_xy_0[j] + 0.5 * fl1_fx * tpz_zz_y_0[j];

            tpx_xzz_xz_0[j] = pa_x[j] * tpx_zz_xz_0[j] + 0.5 * fl1_fx * tpx_zz_z_0[j] - fl1_fgb * fl1_fx * ts_zz_xz_0[j];

            tpy_xzz_xz_0[j] = pa_x[j] * tpy_zz_xz_0[j] + 0.5 * fl1_fx * tpy_zz_z_0[j];

            tpz_xzz_xz_0[j] = pa_x[j] * tpz_zz_xz_0[j] + 0.5 * fl1_fx * tpz_zz_z_0[j];

            tpx_xzz_yy_0[j] = pa_x[j] * tpx_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_zz_yy_0[j];

            tpy_xzz_yy_0[j] = pa_x[j] * tpy_zz_yy_0[j];

            tpz_xzz_yy_0[j] = pa_x[j] * tpz_zz_yy_0[j];

            tpx_xzz_yz_0[j] = pa_x[j] * tpx_zz_yz_0[j] - fl1_fgb * fl1_fx * ts_zz_yz_0[j];

            tpy_xzz_yz_0[j] = pa_x[j] * tpy_zz_yz_0[j];

            tpz_xzz_yz_0[j] = pa_x[j] * tpz_zz_yz_0[j];

            tpx_xzz_zz_0[j] = pa_x[j] * tpx_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_zz_zz_0[j];

            tpy_xzz_zz_0[j] = pa_x[j] * tpy_zz_zz_0[j];

            tpz_xzz_zz_0[j] = pa_x[j] * tpz_zz_zz_0[j];

            tpx_yyy_xx_0[j] = pa_y[j] * tpx_yy_xx_0[j] + fl1_fx * tpx_y_xx_0[j];

            tpy_yyy_xx_0[j] = pa_y[j] * tpy_yy_xx_0[j] + fl1_fx * tpy_y_xx_0[j] - fl1_fgb * fl1_fx * ts_yy_xx_0[j];

            tpz_yyy_xx_0[j] = pa_y[j] * tpz_yy_xx_0[j] + fl1_fx * tpz_y_xx_0[j];

            tpx_yyy_xy_0[j] = pa_y[j] * tpx_yy_xy_0[j] + fl1_fx * tpx_y_xy_0[j] + 0.5 * fl1_fx * tpx_yy_x_0[j];

            tpy_yyy_xy_0[j] = pa_y[j] * tpy_yy_xy_0[j] + fl1_fx * tpy_y_xy_0[j] + 0.5 * fl1_fx * tpy_yy_x_0[j] - fl1_fgb * fl1_fx * ts_yy_xy_0[j];

            tpz_yyy_xy_0[j] = pa_y[j] * tpz_yy_xy_0[j] + fl1_fx * tpz_y_xy_0[j] + 0.5 * fl1_fx * tpz_yy_x_0[j];

            tpx_yyy_xz_0[j] = pa_y[j] * tpx_yy_xz_0[j] + fl1_fx * tpx_y_xz_0[j];

            tpy_yyy_xz_0[j] = pa_y[j] * tpy_yy_xz_0[j] + fl1_fx * tpy_y_xz_0[j] - fl1_fgb * fl1_fx * ts_yy_xz_0[j];

            tpz_yyy_xz_0[j] = pa_y[j] * tpz_yy_xz_0[j] + fl1_fx * tpz_y_xz_0[j];

            tpx_yyy_yy_0[j] = pa_y[j] * tpx_yy_yy_0[j] + fl1_fx * tpx_y_yy_0[j] + fl1_fx * tpx_yy_y_0[j];

            tpy_yyy_yy_0[j] = pa_y[j] * tpy_yy_yy_0[j] + fl1_fx * tpy_y_yy_0[j] + fl1_fx * tpy_yy_y_0[j] - fl1_fgb * fl1_fx * ts_yy_yy_0[j];

            tpz_yyy_yy_0[j] = pa_y[j] * tpz_yy_yy_0[j] + fl1_fx * tpz_y_yy_0[j] + fl1_fx * tpz_yy_y_0[j];

            tpx_yyy_yz_0[j] = pa_y[j] * tpx_yy_yz_0[j] + fl1_fx * tpx_y_yz_0[j] + 0.5 * fl1_fx * tpx_yy_z_0[j];

            tpy_yyy_yz_0[j] = pa_y[j] * tpy_yy_yz_0[j] + fl1_fx * tpy_y_yz_0[j] + 0.5 * fl1_fx * tpy_yy_z_0[j] - fl1_fgb * fl1_fx * ts_yy_yz_0[j];

            tpz_yyy_yz_0[j] = pa_y[j] * tpz_yy_yz_0[j] + fl1_fx * tpz_y_yz_0[j] + 0.5 * fl1_fx * tpz_yy_z_0[j];

            tpx_yyy_zz_0[j] = pa_y[j] * tpx_yy_zz_0[j] + fl1_fx * tpx_y_zz_0[j];

            tpy_yyy_zz_0[j] = pa_y[j] * tpy_yy_zz_0[j] + fl1_fx * tpy_y_zz_0[j] - fl1_fgb * fl1_fx * ts_yy_zz_0[j];

            tpz_yyy_zz_0[j] = pa_y[j] * tpz_yy_zz_0[j] + fl1_fx * tpz_y_zz_0[j];

            tpx_yyz_xx_0[j] = pa_y[j] * tpx_yz_xx_0[j] + 0.5 * fl1_fx * tpx_z_xx_0[j];

            tpy_yyz_xx_0[j] = pa_y[j] * tpy_yz_xx_0[j] + 0.5 * fl1_fx * tpy_z_xx_0[j] - fl1_fgb * fl1_fx * ts_yz_xx_0[j];

            tpz_yyz_xx_0[j] = pa_y[j] * tpz_yz_xx_0[j] + 0.5 * fl1_fx * tpz_z_xx_0[j];

            tpx_yyz_xy_0[j] = pa_y[j] * tpx_yz_xy_0[j] + 0.5 * fl1_fx * tpx_z_xy_0[j] + 0.5 * fl1_fx * tpx_yz_x_0[j];

            tpy_yyz_xy_0[j] =
                pa_y[j] * tpy_yz_xy_0[j] + 0.5 * fl1_fx * tpy_z_xy_0[j] + 0.5 * fl1_fx * tpy_yz_x_0[j] - fl1_fgb * fl1_fx * ts_yz_xy_0[j];

            tpz_yyz_xy_0[j] = pa_y[j] * tpz_yz_xy_0[j] + 0.5 * fl1_fx * tpz_z_xy_0[j] + 0.5 * fl1_fx * tpz_yz_x_0[j];

            tpx_yyz_xz_0[j] = pa_y[j] * tpx_yz_xz_0[j] + 0.5 * fl1_fx * tpx_z_xz_0[j];

            tpy_yyz_xz_0[j] = pa_y[j] * tpy_yz_xz_0[j] + 0.5 * fl1_fx * tpy_z_xz_0[j] - fl1_fgb * fl1_fx * ts_yz_xz_0[j];

            tpz_yyz_xz_0[j] = pa_y[j] * tpz_yz_xz_0[j] + 0.5 * fl1_fx * tpz_z_xz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFD_135_180(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (135,180)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

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

        auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27);

        auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28);

        auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29);

        auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30);

        auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31);

        auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32);

        auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33);

        auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34);

        auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35);

        // set up pointers to integrals

        auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45);

        auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46);

        auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47);

        auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48);

        auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49);

        auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50);

        auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51);

        auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52);

        auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53);

        auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56);

        auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57);

        auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58);

        auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59);

        auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59);

        // Batch of Integrals (135,180)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yyz_yy_0, tpx_yyz_yz_0, tpx_yyz_zz_0, tpx_yz_y_0, \
                                     tpx_yz_yy_0, tpx_yz_yz_0, tpx_yz_z_0, tpx_yz_zz_0, tpx_yzz_xx_0, tpx_yzz_xy_0, \
                                     tpx_yzz_xz_0, tpx_yzz_yy_0, tpx_yzz_yz_0, tpx_yzz_zz_0, tpx_z_xx_0, tpx_z_xy_0, \
                                     tpx_z_xz_0, tpx_z_yy_0, tpx_z_yz_0, tpx_z_zz_0, tpx_zz_x_0, tpx_zz_xx_0, \
                                     tpx_zz_xy_0, tpx_zz_xz_0, tpx_zz_y_0, tpx_zz_yy_0, tpx_zz_yz_0, tpx_zz_z_0, \
                                     tpx_zz_zz_0, tpx_zzz_xx_0, tpx_zzz_xy_0, tpx_zzz_xz_0, tpx_zzz_yy_0, tpx_zzz_yz_0, \
                                     tpx_zzz_zz_0, tpy_yyz_yy_0, tpy_yyz_yz_0, tpy_yyz_zz_0, tpy_yz_y_0, tpy_yz_yy_0, \
                                     tpy_yz_yz_0, tpy_yz_z_0, tpy_yz_zz_0, tpy_yzz_xx_0, tpy_yzz_xy_0, tpy_yzz_xz_0, \
                                     tpy_yzz_yy_0, tpy_yzz_yz_0, tpy_yzz_zz_0, tpy_z_xx_0, tpy_z_xy_0, tpy_z_xz_0, \
                                     tpy_z_yy_0, tpy_z_yz_0, tpy_z_zz_0, tpy_zz_x_0, tpy_zz_xx_0, tpy_zz_xy_0, \
                                     tpy_zz_xz_0, tpy_zz_y_0, tpy_zz_yy_0, tpy_zz_yz_0, tpy_zz_z_0, tpy_zz_zz_0, \
                                     tpy_zzz_xx_0, tpy_zzz_xy_0, tpy_zzz_xz_0, tpy_zzz_yy_0, tpy_zzz_yz_0, tpy_zzz_zz_0, \
                                     tpz_yyz_yy_0, tpz_yyz_yz_0, tpz_yyz_zz_0, tpz_yz_y_0, tpz_yz_yy_0, tpz_yz_yz_0, \
                                     tpz_yz_z_0, tpz_yz_zz_0, tpz_yzz_xx_0, tpz_yzz_xy_0, tpz_yzz_xz_0, tpz_yzz_yy_0, \
                                     tpz_yzz_yz_0, tpz_yzz_zz_0, tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_z_yy_0, \
                                     tpz_z_yz_0, tpz_z_zz_0, tpz_zz_x_0, tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, \
                                     tpz_zz_y_0, tpz_zz_yy_0, tpz_zz_yz_0, tpz_zz_z_0, tpz_zz_zz_0, tpz_zzz_xx_0, \
                                     tpz_zzz_xy_0, tpz_zzz_xz_0, tpz_zzz_yy_0, tpz_zzz_yz_0, tpz_zzz_zz_0, ts_yz_yy_0, \
                                     ts_yz_yz_0, ts_yz_zz_0, ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, ts_zz_yy_0, ts_zz_yz_0, \
                                     ts_zz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyz_yy_0[j] = pa_y[j] * tpx_yz_yy_0[j] + 0.5 * fl1_fx * tpx_z_yy_0[j] + fl1_fx * tpx_yz_y_0[j];

            tpy_yyz_yy_0[j] = pa_y[j] * tpy_yz_yy_0[j] + 0.5 * fl1_fx * tpy_z_yy_0[j] + fl1_fx * tpy_yz_y_0[j] - fl1_fgb * fl1_fx * ts_yz_yy_0[j];

            tpz_yyz_yy_0[j] = pa_y[j] * tpz_yz_yy_0[j] + 0.5 * fl1_fx * tpz_z_yy_0[j] + fl1_fx * tpz_yz_y_0[j];

            tpx_yyz_yz_0[j] = pa_y[j] * tpx_yz_yz_0[j] + 0.5 * fl1_fx * tpx_z_yz_0[j] + 0.5 * fl1_fx * tpx_yz_z_0[j];

            tpy_yyz_yz_0[j] =
                pa_y[j] * tpy_yz_yz_0[j] + 0.5 * fl1_fx * tpy_z_yz_0[j] + 0.5 * fl1_fx * tpy_yz_z_0[j] - fl1_fgb * fl1_fx * ts_yz_yz_0[j];

            tpz_yyz_yz_0[j] = pa_y[j] * tpz_yz_yz_0[j] + 0.5 * fl1_fx * tpz_z_yz_0[j] + 0.5 * fl1_fx * tpz_yz_z_0[j];

            tpx_yyz_zz_0[j] = pa_y[j] * tpx_yz_zz_0[j] + 0.5 * fl1_fx * tpx_z_zz_0[j];

            tpy_yyz_zz_0[j] = pa_y[j] * tpy_yz_zz_0[j] + 0.5 * fl1_fx * tpy_z_zz_0[j] - fl1_fgb * fl1_fx * ts_yz_zz_0[j];

            tpz_yyz_zz_0[j] = pa_y[j] * tpz_yz_zz_0[j] + 0.5 * fl1_fx * tpz_z_zz_0[j];

            tpx_yzz_xx_0[j] = pa_y[j] * tpx_zz_xx_0[j];

            tpy_yzz_xx_0[j] = pa_y[j] * tpy_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_zz_xx_0[j];

            tpz_yzz_xx_0[j] = pa_y[j] * tpz_zz_xx_0[j];

            tpx_yzz_xy_0[j] = pa_y[j] * tpx_zz_xy_0[j] + 0.5 * fl1_fx * tpx_zz_x_0[j];

            tpy_yzz_xy_0[j] = pa_y[j] * tpy_zz_xy_0[j] + 0.5 * fl1_fx * tpy_zz_x_0[j] - fl1_fgb * fl1_fx * ts_zz_xy_0[j];

            tpz_yzz_xy_0[j] = pa_y[j] * tpz_zz_xy_0[j] + 0.5 * fl1_fx * tpz_zz_x_0[j];

            tpx_yzz_xz_0[j] = pa_y[j] * tpx_zz_xz_0[j];

            tpy_yzz_xz_0[j] = pa_y[j] * tpy_zz_xz_0[j] - fl1_fgb * fl1_fx * ts_zz_xz_0[j];

            tpz_yzz_xz_0[j] = pa_y[j] * tpz_zz_xz_0[j];

            tpx_yzz_yy_0[j] = pa_y[j] * tpx_zz_yy_0[j] + fl1_fx * tpx_zz_y_0[j];

            tpy_yzz_yy_0[j] = pa_y[j] * tpy_zz_yy_0[j] + fl1_fx * tpy_zz_y_0[j] - fl1_fgb * fl1_fx * ts_zz_yy_0[j];

            tpz_yzz_yy_0[j] = pa_y[j] * tpz_zz_yy_0[j] + fl1_fx * tpz_zz_y_0[j];

            tpx_yzz_yz_0[j] = pa_y[j] * tpx_zz_yz_0[j] + 0.5 * fl1_fx * tpx_zz_z_0[j];

            tpy_yzz_yz_0[j] = pa_y[j] * tpy_zz_yz_0[j] + 0.5 * fl1_fx * tpy_zz_z_0[j] - fl1_fgb * fl1_fx * ts_zz_yz_0[j];

            tpz_yzz_yz_0[j] = pa_y[j] * tpz_zz_yz_0[j] + 0.5 * fl1_fx * tpz_zz_z_0[j];

            tpx_yzz_zz_0[j] = pa_y[j] * tpx_zz_zz_0[j];

            tpy_yzz_zz_0[j] = pa_y[j] * tpy_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_zz_zz_0[j];

            tpz_yzz_zz_0[j] = pa_y[j] * tpz_zz_zz_0[j];

            tpx_zzz_xx_0[j] = pa_z[j] * tpx_zz_xx_0[j] + fl1_fx * tpx_z_xx_0[j];

            tpy_zzz_xx_0[j] = pa_z[j] * tpy_zz_xx_0[j] + fl1_fx * tpy_z_xx_0[j];

            tpz_zzz_xx_0[j] = pa_z[j] * tpz_zz_xx_0[j] + fl1_fx * tpz_z_xx_0[j] - fl1_fgb * fl1_fx * ts_zz_xx_0[j];

            tpx_zzz_xy_0[j] = pa_z[j] * tpx_zz_xy_0[j] + fl1_fx * tpx_z_xy_0[j];

            tpy_zzz_xy_0[j] = pa_z[j] * tpy_zz_xy_0[j] + fl1_fx * tpy_z_xy_0[j];

            tpz_zzz_xy_0[j] = pa_z[j] * tpz_zz_xy_0[j] + fl1_fx * tpz_z_xy_0[j] - fl1_fgb * fl1_fx * ts_zz_xy_0[j];

            tpx_zzz_xz_0[j] = pa_z[j] * tpx_zz_xz_0[j] + fl1_fx * tpx_z_xz_0[j] + 0.5 * fl1_fx * tpx_zz_x_0[j];

            tpy_zzz_xz_0[j] = pa_z[j] * tpy_zz_xz_0[j] + fl1_fx * tpy_z_xz_0[j] + 0.5 * fl1_fx * tpy_zz_x_0[j];

            tpz_zzz_xz_0[j] = pa_z[j] * tpz_zz_xz_0[j] + fl1_fx * tpz_z_xz_0[j] + 0.5 * fl1_fx * tpz_zz_x_0[j] - fl1_fgb * fl1_fx * ts_zz_xz_0[j];

            tpx_zzz_yy_0[j] = pa_z[j] * tpx_zz_yy_0[j] + fl1_fx * tpx_z_yy_0[j];

            tpy_zzz_yy_0[j] = pa_z[j] * tpy_zz_yy_0[j] + fl1_fx * tpy_z_yy_0[j];

            tpz_zzz_yy_0[j] = pa_z[j] * tpz_zz_yy_0[j] + fl1_fx * tpz_z_yy_0[j] - fl1_fgb * fl1_fx * ts_zz_yy_0[j];

            tpx_zzz_yz_0[j] = pa_z[j] * tpx_zz_yz_0[j] + fl1_fx * tpx_z_yz_0[j] + 0.5 * fl1_fx * tpx_zz_y_0[j];

            tpy_zzz_yz_0[j] = pa_z[j] * tpy_zz_yz_0[j] + fl1_fx * tpy_z_yz_0[j] + 0.5 * fl1_fx * tpy_zz_y_0[j];

            tpz_zzz_yz_0[j] = pa_z[j] * tpz_zz_yz_0[j] + fl1_fx * tpz_z_yz_0[j] + 0.5 * fl1_fx * tpz_zz_y_0[j] - fl1_fgb * fl1_fx * ts_zz_yz_0[j];

            tpx_zzz_zz_0[j] = pa_z[j] * tpx_zz_zz_0[j] + fl1_fx * tpx_z_zz_0[j] + fl1_fx * tpx_zz_z_0[j];

            tpy_zzz_zz_0[j] = pa_z[j] * tpy_zz_zz_0[j] + fl1_fx * tpy_z_zz_0[j] + fl1_fx * tpy_zz_z_0[j];

            tpz_zzz_zz_0[j] = pa_z[j] * tpz_zz_zz_0[j] + fl1_fx * tpz_z_zz_0[j] + fl1_fx * tpz_zz_z_0[j] - fl1_fgb * fl1_fx * ts_zz_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForDG_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDG_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDG_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDG_135_180(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDG_180_225(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDG_225_270(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForDG_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx);

        auto tpy_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx);

        auto tpz_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx);

        auto tpx_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 1);

        auto tpy_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 1);

        auto tpz_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 1);

        auto tpx_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 2);

        auto tpy_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 2);

        auto tpz_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 2);

        auto tpx_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 3);

        auto tpy_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 3);

        auto tpz_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 3);

        auto tpx_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 4);

        auto tpy_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 4);

        auto tpz_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 4);

        auto tpx_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 5);

        auto tpy_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 5);

        auto tpz_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 5);

        auto tpx_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 6);

        auto tpy_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 6);

        auto tpz_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 6);

        auto tpx_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 7);

        auto tpy_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 7);

        auto tpz_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 7);

        auto tpx_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 8);

        auto tpy_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 8);

        auto tpz_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 8);

        auto tpx_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 9);

        auto tpy_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 9);

        auto tpz_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 9);

        auto tpx_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 10);

        auto tpy_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 10);

        auto tpz_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 10);

        auto tpx_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 11);

        auto tpy_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 11);

        auto tpz_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 11);

        auto tpx_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 12);

        auto tpy_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 12);

        auto tpz_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 12);

        auto tpx_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 13);

        auto tpy_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 13);

        auto tpz_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 13);

        auto tpx_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 14);

        auto tpy_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 14);

        auto tpz_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 14);

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx);

        auto tpy_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx);

        auto tpz_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx);

        auto tpx_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 1);

        auto tpy_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 2);

        auto tpy_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 3);

        auto tpy_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 4);

        auto tpy_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 5);

        auto tpy_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 6);

        auto tpy_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 7);

        auto tpy_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 8);

        auto tpy_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 9);

        auto tpy_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 9);

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

        // set up pointers to integrals

        auto tpx_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx);

        auto tpy_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx);

        auto tpz_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx);

        auto tpx_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 1);

        auto tpy_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tpz_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tpx_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 2);

        auto tpy_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tpz_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tpx_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 3);

        auto tpy_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tpz_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tpx_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 4);

        auto tpy_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tpz_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tpx_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 5);

        auto tpy_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tpz_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tpx_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 6);

        auto tpy_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tpz_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tpx_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 7);

        auto tpy_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tpz_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tpx_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 8);

        auto tpy_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tpz_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tpx_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 9);

        auto tpy_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tpz_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tpx_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 10);

        auto tpy_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tpz_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tpx_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 11);

        auto tpy_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tpz_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tpx_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 12);

        auto tpy_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tpz_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tpx_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 13);

        auto tpy_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tpz_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tpx_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 14);

        auto tpy_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tpz_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxyy_0, \
                                     tpx_0_xxyz_0, tpx_0_xxzz_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyzz_0, tpx_0_yzzz_0, tpx_0_zzzz_0, tpx_x_xxx_0, \
                                     tpx_x_xxxx_0, tpx_x_xxxy_0, tpx_x_xxxz_0, tpx_x_xxy_0, tpx_x_xxyy_0, tpx_x_xxyz_0, \
                                     tpx_x_xxz_0, tpx_x_xxzz_0, tpx_x_xyy_0, tpx_x_xyyy_0, tpx_x_xyyz_0, tpx_x_xyz_0, \
                                     tpx_x_xyzz_0, tpx_x_xzz_0, tpx_x_xzzz_0, tpx_x_yyy_0, tpx_x_yyyy_0, tpx_x_yyyz_0, \
                                     tpx_x_yyz_0, tpx_x_yyzz_0, tpx_x_yzz_0, tpx_x_yzzz_0, tpx_x_zzz_0, tpx_x_zzzz_0, \
                                     tpx_xx_xxxx_0, tpx_xx_xxxy_0, tpx_xx_xxxz_0, tpx_xx_xxyy_0, tpx_xx_xxyz_0, \
                                     tpx_xx_xxzz_0, tpx_xx_xyyy_0, tpx_xx_xyyz_0, tpx_xx_xyzz_0, tpx_xx_xzzz_0, \
                                     tpx_xx_yyyy_0, tpx_xx_yyyz_0, tpx_xx_yyzz_0, tpx_xx_yzzz_0, tpx_xx_zzzz_0, \
                                     tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxyy_0, tpy_0_xxyz_0, tpy_0_xxzz_0, \
                                     tpy_0_xyyy_0, tpy_0_xyyz_0, tpy_0_xyzz_0, tpy_0_xzzz_0, tpy_0_yyyy_0, tpy_0_yyyz_0, \
                                     tpy_0_yyzz_0, tpy_0_yzzz_0, tpy_0_zzzz_0, tpy_x_xxx_0, tpy_x_xxxx_0, tpy_x_xxxy_0, \
                                     tpy_x_xxxz_0, tpy_x_xxy_0, tpy_x_xxyy_0, tpy_x_xxyz_0, tpy_x_xxz_0, tpy_x_xxzz_0, \
                                     tpy_x_xyy_0, tpy_x_xyyy_0, tpy_x_xyyz_0, tpy_x_xyz_0, tpy_x_xyzz_0, tpy_x_xzz_0, \
                                     tpy_x_xzzz_0, tpy_x_yyy_0, tpy_x_yyyy_0, tpy_x_yyyz_0, tpy_x_yyz_0, tpy_x_yyzz_0, \
                                     tpy_x_yzz_0, tpy_x_yzzz_0, tpy_x_zzz_0, tpy_x_zzzz_0, tpy_xx_xxxx_0, \
                                     tpy_xx_xxxy_0, tpy_xx_xxxz_0, tpy_xx_xxyy_0, tpy_xx_xxyz_0, tpy_xx_xxzz_0, \
                                     tpy_xx_xyyy_0, tpy_xx_xyyz_0, tpy_xx_xyzz_0, tpy_xx_xzzz_0, tpy_xx_yyyy_0, \
                                     tpy_xx_yyyz_0, tpy_xx_yyzz_0, tpy_xx_yzzz_0, tpy_xx_zzzz_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxzz_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyzz_0, tpz_0_xzzz_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yzzz_0, tpz_0_zzzz_0, tpz_x_xxx_0, tpz_x_xxxx_0, tpz_x_xxxy_0, tpz_x_xxxz_0, \
                                     tpz_x_xxy_0, tpz_x_xxyy_0, tpz_x_xxyz_0, tpz_x_xxz_0, tpz_x_xxzz_0, tpz_x_xyy_0, \
                                     tpz_x_xyyy_0, tpz_x_xyyz_0, tpz_x_xyz_0, tpz_x_xyzz_0, tpz_x_xzz_0, tpz_x_xzzz_0, \
                                     tpz_x_yyy_0, tpz_x_yyyy_0, tpz_x_yyyz_0, tpz_x_yyz_0, tpz_x_yyzz_0, tpz_x_yzz_0, \
                                     tpz_x_yzzz_0, tpz_x_zzz_0, tpz_x_zzzz_0, tpz_xx_xxxx_0, tpz_xx_xxxy_0, \
                                     tpz_xx_xxxz_0, tpz_xx_xxyy_0, tpz_xx_xxyz_0, tpz_xx_xxzz_0, tpz_xx_xyyy_0, \
                                     tpz_xx_xyyz_0, tpz_xx_xyzz_0, tpz_xx_xzzz_0, tpz_xx_yyyy_0, tpz_xx_yyyz_0, \
                                     tpz_xx_yyzz_0, tpz_xx_yzzz_0, tpz_xx_zzzz_0, ts_x_xxxx_0, ts_x_xxxy_0, ts_x_xxxz_0, \
                                     ts_x_xxyy_0, ts_x_xxyz_0, ts_x_xxzz_0, ts_x_xyyy_0, ts_x_xyyz_0, ts_x_xyzz_0, \
                                     ts_x_xzzz_0, ts_x_yyyy_0, ts_x_yyyz_0, ts_x_yyzz_0, ts_x_yzzz_0, ts_x_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xx_xxxx_0[j] =
                pa_x[j] * tpx_x_xxxx_0[j] + 0.5 * fl1_fx * tpx_0_xxxx_0[j] + 2.0 * fl1_fx * tpx_x_xxx_0[j] - fl1_fgb * fl1_fx * ts_x_xxxx_0[j];

            tpy_xx_xxxx_0[j] = pa_x[j] * tpy_x_xxxx_0[j] + 0.5 * fl1_fx * tpy_0_xxxx_0[j] + 2.0 * fl1_fx * tpy_x_xxx_0[j];

            tpz_xx_xxxx_0[j] = pa_x[j] * tpz_x_xxxx_0[j] + 0.5 * fl1_fx * tpz_0_xxxx_0[j] + 2.0 * fl1_fx * tpz_x_xxx_0[j];

            tpx_xx_xxxy_0[j] =
                pa_x[j] * tpx_x_xxxy_0[j] + 0.5 * fl1_fx * tpx_0_xxxy_0[j] + 1.5 * fl1_fx * tpx_x_xxy_0[j] - fl1_fgb * fl1_fx * ts_x_xxxy_0[j];

            tpy_xx_xxxy_0[j] = pa_x[j] * tpy_x_xxxy_0[j] + 0.5 * fl1_fx * tpy_0_xxxy_0[j] + 1.5 * fl1_fx * tpy_x_xxy_0[j];

            tpz_xx_xxxy_0[j] = pa_x[j] * tpz_x_xxxy_0[j] + 0.5 * fl1_fx * tpz_0_xxxy_0[j] + 1.5 * fl1_fx * tpz_x_xxy_0[j];

            tpx_xx_xxxz_0[j] =
                pa_x[j] * tpx_x_xxxz_0[j] + 0.5 * fl1_fx * tpx_0_xxxz_0[j] + 1.5 * fl1_fx * tpx_x_xxz_0[j] - fl1_fgb * fl1_fx * ts_x_xxxz_0[j];

            tpy_xx_xxxz_0[j] = pa_x[j] * tpy_x_xxxz_0[j] + 0.5 * fl1_fx * tpy_0_xxxz_0[j] + 1.5 * fl1_fx * tpy_x_xxz_0[j];

            tpz_xx_xxxz_0[j] = pa_x[j] * tpz_x_xxxz_0[j] + 0.5 * fl1_fx * tpz_0_xxxz_0[j] + 1.5 * fl1_fx * tpz_x_xxz_0[j];

            tpx_xx_xxyy_0[j] =
                pa_x[j] * tpx_x_xxyy_0[j] + 0.5 * fl1_fx * tpx_0_xxyy_0[j] + fl1_fx * tpx_x_xyy_0[j] - fl1_fgb * fl1_fx * ts_x_xxyy_0[j];

            tpy_xx_xxyy_0[j] = pa_x[j] * tpy_x_xxyy_0[j] + 0.5 * fl1_fx * tpy_0_xxyy_0[j] + fl1_fx * tpy_x_xyy_0[j];

            tpz_xx_xxyy_0[j] = pa_x[j] * tpz_x_xxyy_0[j] + 0.5 * fl1_fx * tpz_0_xxyy_0[j] + fl1_fx * tpz_x_xyy_0[j];

            tpx_xx_xxyz_0[j] =
                pa_x[j] * tpx_x_xxyz_0[j] + 0.5 * fl1_fx * tpx_0_xxyz_0[j] + fl1_fx * tpx_x_xyz_0[j] - fl1_fgb * fl1_fx * ts_x_xxyz_0[j];

            tpy_xx_xxyz_0[j] = pa_x[j] * tpy_x_xxyz_0[j] + 0.5 * fl1_fx * tpy_0_xxyz_0[j] + fl1_fx * tpy_x_xyz_0[j];

            tpz_xx_xxyz_0[j] = pa_x[j] * tpz_x_xxyz_0[j] + 0.5 * fl1_fx * tpz_0_xxyz_0[j] + fl1_fx * tpz_x_xyz_0[j];

            tpx_xx_xxzz_0[j] =
                pa_x[j] * tpx_x_xxzz_0[j] + 0.5 * fl1_fx * tpx_0_xxzz_0[j] + fl1_fx * tpx_x_xzz_0[j] - fl1_fgb * fl1_fx * ts_x_xxzz_0[j];

            tpy_xx_xxzz_0[j] = pa_x[j] * tpy_x_xxzz_0[j] + 0.5 * fl1_fx * tpy_0_xxzz_0[j] + fl1_fx * tpy_x_xzz_0[j];

            tpz_xx_xxzz_0[j] = pa_x[j] * tpz_x_xxzz_0[j] + 0.5 * fl1_fx * tpz_0_xxzz_0[j] + fl1_fx * tpz_x_xzz_0[j];

            tpx_xx_xyyy_0[j] =
                pa_x[j] * tpx_x_xyyy_0[j] + 0.5 * fl1_fx * tpx_0_xyyy_0[j] + 0.5 * fl1_fx * tpx_x_yyy_0[j] - fl1_fgb * fl1_fx * ts_x_xyyy_0[j];

            tpy_xx_xyyy_0[j] = pa_x[j] * tpy_x_xyyy_0[j] + 0.5 * fl1_fx * tpy_0_xyyy_0[j] + 0.5 * fl1_fx * tpy_x_yyy_0[j];

            tpz_xx_xyyy_0[j] = pa_x[j] * tpz_x_xyyy_0[j] + 0.5 * fl1_fx * tpz_0_xyyy_0[j] + 0.5 * fl1_fx * tpz_x_yyy_0[j];

            tpx_xx_xyyz_0[j] =
                pa_x[j] * tpx_x_xyyz_0[j] + 0.5 * fl1_fx * tpx_0_xyyz_0[j] + 0.5 * fl1_fx * tpx_x_yyz_0[j] - fl1_fgb * fl1_fx * ts_x_xyyz_0[j];

            tpy_xx_xyyz_0[j] = pa_x[j] * tpy_x_xyyz_0[j] + 0.5 * fl1_fx * tpy_0_xyyz_0[j] + 0.5 * fl1_fx * tpy_x_yyz_0[j];

            tpz_xx_xyyz_0[j] = pa_x[j] * tpz_x_xyyz_0[j] + 0.5 * fl1_fx * tpz_0_xyyz_0[j] + 0.5 * fl1_fx * tpz_x_yyz_0[j];

            tpx_xx_xyzz_0[j] =
                pa_x[j] * tpx_x_xyzz_0[j] + 0.5 * fl1_fx * tpx_0_xyzz_0[j] + 0.5 * fl1_fx * tpx_x_yzz_0[j] - fl1_fgb * fl1_fx * ts_x_xyzz_0[j];

            tpy_xx_xyzz_0[j] = pa_x[j] * tpy_x_xyzz_0[j] + 0.5 * fl1_fx * tpy_0_xyzz_0[j] + 0.5 * fl1_fx * tpy_x_yzz_0[j];

            tpz_xx_xyzz_0[j] = pa_x[j] * tpz_x_xyzz_0[j] + 0.5 * fl1_fx * tpz_0_xyzz_0[j] + 0.5 * fl1_fx * tpz_x_yzz_0[j];

            tpx_xx_xzzz_0[j] =
                pa_x[j] * tpx_x_xzzz_0[j] + 0.5 * fl1_fx * tpx_0_xzzz_0[j] + 0.5 * fl1_fx * tpx_x_zzz_0[j] - fl1_fgb * fl1_fx * ts_x_xzzz_0[j];

            tpy_xx_xzzz_0[j] = pa_x[j] * tpy_x_xzzz_0[j] + 0.5 * fl1_fx * tpy_0_xzzz_0[j] + 0.5 * fl1_fx * tpy_x_zzz_0[j];

            tpz_xx_xzzz_0[j] = pa_x[j] * tpz_x_xzzz_0[j] + 0.5 * fl1_fx * tpz_0_xzzz_0[j] + 0.5 * fl1_fx * tpz_x_zzz_0[j];

            tpx_xx_yyyy_0[j] = pa_x[j] * tpx_x_yyyy_0[j] + 0.5 * fl1_fx * tpx_0_yyyy_0[j] - fl1_fgb * fl1_fx * ts_x_yyyy_0[j];

            tpy_xx_yyyy_0[j] = pa_x[j] * tpy_x_yyyy_0[j] + 0.5 * fl1_fx * tpy_0_yyyy_0[j];

            tpz_xx_yyyy_0[j] = pa_x[j] * tpz_x_yyyy_0[j] + 0.5 * fl1_fx * tpz_0_yyyy_0[j];

            tpx_xx_yyyz_0[j] = pa_x[j] * tpx_x_yyyz_0[j] + 0.5 * fl1_fx * tpx_0_yyyz_0[j] - fl1_fgb * fl1_fx * ts_x_yyyz_0[j];

            tpy_xx_yyyz_0[j] = pa_x[j] * tpy_x_yyyz_0[j] + 0.5 * fl1_fx * tpy_0_yyyz_0[j];

            tpz_xx_yyyz_0[j] = pa_x[j] * tpz_x_yyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyyz_0[j];

            tpx_xx_yyzz_0[j] = pa_x[j] * tpx_x_yyzz_0[j] + 0.5 * fl1_fx * tpx_0_yyzz_0[j] - fl1_fgb * fl1_fx * ts_x_yyzz_0[j];

            tpy_xx_yyzz_0[j] = pa_x[j] * tpy_x_yyzz_0[j] + 0.5 * fl1_fx * tpy_0_yyzz_0[j];

            tpz_xx_yyzz_0[j] = pa_x[j] * tpz_x_yyzz_0[j] + 0.5 * fl1_fx * tpz_0_yyzz_0[j];

            tpx_xx_yzzz_0[j] = pa_x[j] * tpx_x_yzzz_0[j] + 0.5 * fl1_fx * tpx_0_yzzz_0[j] - fl1_fgb * fl1_fx * ts_x_yzzz_0[j];

            tpy_xx_yzzz_0[j] = pa_x[j] * tpy_x_yzzz_0[j] + 0.5 * fl1_fx * tpy_0_yzzz_0[j];

            tpz_xx_yzzz_0[j] = pa_x[j] * tpz_x_yzzz_0[j] + 0.5 * fl1_fx * tpz_0_yzzz_0[j];

            tpx_xx_zzzz_0[j] = pa_x[j] * tpx_x_zzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzzz_0[j] - fl1_fgb * fl1_fx * ts_x_zzzz_0[j];

            tpy_xx_zzzz_0[j] = pa_x[j] * tpy_x_zzzz_0[j] + 0.5 * fl1_fx * tpy_0_zzzz_0[j];

            tpz_xx_zzzz_0[j] = pa_x[j] * tpz_x_zzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15);

        auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16);

        auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17);

        auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18);

        auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19);

        auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20);

        auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21);

        auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22);

        auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23);

        auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24);

        auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25);

        auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26);

        auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27);

        auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28);

        auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29);

        auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

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

        // set up pointers to integrals

        auto tpx_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 15);

        auto tpy_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tpz_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tpx_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 16);

        auto tpy_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tpz_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tpx_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 17);

        auto tpy_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tpz_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tpx_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 18);

        auto tpy_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tpz_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tpx_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 19);

        auto tpy_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tpz_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tpx_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 20);

        auto tpy_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tpz_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tpx_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 21);

        auto tpy_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tpz_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tpx_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 22);

        auto tpy_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tpz_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tpx_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 23);

        auto tpy_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tpz_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tpx_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 24);

        auto tpy_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tpz_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tpx_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 25);

        auto tpy_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tpz_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tpx_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 26);

        auto tpy_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tpz_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tpx_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 27);

        auto tpy_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tpz_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tpx_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 28);

        auto tpy_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tpz_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tpx_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 29);

        auto tpy_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tpz_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xy_xxxx_0, tpx_xy_xxxy_0, tpx_xy_xxxz_0, \
                                     tpx_xy_xxyy_0, tpx_xy_xxyz_0, tpx_xy_xxzz_0, tpx_xy_xyyy_0, tpx_xy_xyyz_0, \
                                     tpx_xy_xyzz_0, tpx_xy_xzzz_0, tpx_xy_yyyy_0, tpx_xy_yyyz_0, tpx_xy_yyzz_0, \
                                     tpx_xy_yzzz_0, tpx_xy_zzzz_0, tpx_y_xxx_0, tpx_y_xxxx_0, tpx_y_xxxy_0, tpx_y_xxxz_0, \
                                     tpx_y_xxy_0, tpx_y_xxyy_0, tpx_y_xxyz_0, tpx_y_xxz_0, tpx_y_xxzz_0, tpx_y_xyy_0, \
                                     tpx_y_xyyy_0, tpx_y_xyyz_0, tpx_y_xyz_0, tpx_y_xyzz_0, tpx_y_xzz_0, tpx_y_xzzz_0, \
                                     tpx_y_yyy_0, tpx_y_yyyy_0, tpx_y_yyyz_0, tpx_y_yyz_0, tpx_y_yyzz_0, tpx_y_yzz_0, \
                                     tpx_y_yzzz_0, tpx_y_zzz_0, tpx_y_zzzz_0, tpy_xy_xxxx_0, tpy_xy_xxxy_0, \
                                     tpy_xy_xxxz_0, tpy_xy_xxyy_0, tpy_xy_xxyz_0, tpy_xy_xxzz_0, tpy_xy_xyyy_0, \
                                     tpy_xy_xyyz_0, tpy_xy_xyzz_0, tpy_xy_xzzz_0, tpy_xy_yyyy_0, tpy_xy_yyyz_0, \
                                     tpy_xy_yyzz_0, tpy_xy_yzzz_0, tpy_xy_zzzz_0, tpy_y_xxx_0, tpy_y_xxxx_0, \
                                     tpy_y_xxxy_0, tpy_y_xxxz_0, tpy_y_xxy_0, tpy_y_xxyy_0, tpy_y_xxyz_0, tpy_y_xxz_0, \
                                     tpy_y_xxzz_0, tpy_y_xyy_0, tpy_y_xyyy_0, tpy_y_xyyz_0, tpy_y_xyz_0, tpy_y_xyzz_0, \
                                     tpy_y_xzz_0, tpy_y_xzzz_0, tpy_y_yyy_0, tpy_y_yyyy_0, tpy_y_yyyz_0, tpy_y_yyz_0, \
                                     tpy_y_yyzz_0, tpy_y_yzz_0, tpy_y_yzzz_0, tpy_y_zzz_0, tpy_y_zzzz_0, tpz_xy_xxxx_0, \
                                     tpz_xy_xxxy_0, tpz_xy_xxxz_0, tpz_xy_xxyy_0, tpz_xy_xxyz_0, tpz_xy_xxzz_0, \
                                     tpz_xy_xyyy_0, tpz_xy_xyyz_0, tpz_xy_xyzz_0, tpz_xy_xzzz_0, tpz_xy_yyyy_0, \
                                     tpz_xy_yyyz_0, tpz_xy_yyzz_0, tpz_xy_yzzz_0, tpz_xy_zzzz_0, tpz_y_xxx_0, \
                                     tpz_y_xxxx_0, tpz_y_xxxy_0, tpz_y_xxxz_0, tpz_y_xxy_0, tpz_y_xxyy_0, tpz_y_xxyz_0, \
                                     tpz_y_xxz_0, tpz_y_xxzz_0, tpz_y_xyy_0, tpz_y_xyyy_0, tpz_y_xyyz_0, tpz_y_xyz_0, \
                                     tpz_y_xyzz_0, tpz_y_xzz_0, tpz_y_xzzz_0, tpz_y_yyy_0, tpz_y_yyyy_0, tpz_y_yyyz_0, \
                                     tpz_y_yyz_0, tpz_y_yyzz_0, tpz_y_yzz_0, tpz_y_yzzz_0, tpz_y_zzz_0, tpz_y_zzzz_0, \
                                     ts_y_xxxx_0, ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, \
                                     ts_y_xyyy_0, ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, \
                                     ts_y_yyzz_0, ts_y_yzzz_0, ts_y_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xy_xxxx_0[j] = pa_x[j] * tpx_y_xxxx_0[j] + 2.0 * fl1_fx * tpx_y_xxx_0[j] - fl1_fgb * fl1_fx * ts_y_xxxx_0[j];

            tpy_xy_xxxx_0[j] = pa_x[j] * tpy_y_xxxx_0[j] + 2.0 * fl1_fx * tpy_y_xxx_0[j];

            tpz_xy_xxxx_0[j] = pa_x[j] * tpz_y_xxxx_0[j] + 2.0 * fl1_fx * tpz_y_xxx_0[j];

            tpx_xy_xxxy_0[j] = pa_x[j] * tpx_y_xxxy_0[j] + 1.5 * fl1_fx * tpx_y_xxy_0[j] - fl1_fgb * fl1_fx * ts_y_xxxy_0[j];

            tpy_xy_xxxy_0[j] = pa_x[j] * tpy_y_xxxy_0[j] + 1.5 * fl1_fx * tpy_y_xxy_0[j];

            tpz_xy_xxxy_0[j] = pa_x[j] * tpz_y_xxxy_0[j] + 1.5 * fl1_fx * tpz_y_xxy_0[j];

            tpx_xy_xxxz_0[j] = pa_x[j] * tpx_y_xxxz_0[j] + 1.5 * fl1_fx * tpx_y_xxz_0[j] - fl1_fgb * fl1_fx * ts_y_xxxz_0[j];

            tpy_xy_xxxz_0[j] = pa_x[j] * tpy_y_xxxz_0[j] + 1.5 * fl1_fx * tpy_y_xxz_0[j];

            tpz_xy_xxxz_0[j] = pa_x[j] * tpz_y_xxxz_0[j] + 1.5 * fl1_fx * tpz_y_xxz_0[j];

            tpx_xy_xxyy_0[j] = pa_x[j] * tpx_y_xxyy_0[j] + fl1_fx * tpx_y_xyy_0[j] - fl1_fgb * fl1_fx * ts_y_xxyy_0[j];

            tpy_xy_xxyy_0[j] = pa_x[j] * tpy_y_xxyy_0[j] + fl1_fx * tpy_y_xyy_0[j];

            tpz_xy_xxyy_0[j] = pa_x[j] * tpz_y_xxyy_0[j] + fl1_fx * tpz_y_xyy_0[j];

            tpx_xy_xxyz_0[j] = pa_x[j] * tpx_y_xxyz_0[j] + fl1_fx * tpx_y_xyz_0[j] - fl1_fgb * fl1_fx * ts_y_xxyz_0[j];

            tpy_xy_xxyz_0[j] = pa_x[j] * tpy_y_xxyz_0[j] + fl1_fx * tpy_y_xyz_0[j];

            tpz_xy_xxyz_0[j] = pa_x[j] * tpz_y_xxyz_0[j] + fl1_fx * tpz_y_xyz_0[j];

            tpx_xy_xxzz_0[j] = pa_x[j] * tpx_y_xxzz_0[j] + fl1_fx * tpx_y_xzz_0[j] - fl1_fgb * fl1_fx * ts_y_xxzz_0[j];

            tpy_xy_xxzz_0[j] = pa_x[j] * tpy_y_xxzz_0[j] + fl1_fx * tpy_y_xzz_0[j];

            tpz_xy_xxzz_0[j] = pa_x[j] * tpz_y_xxzz_0[j] + fl1_fx * tpz_y_xzz_0[j];

            tpx_xy_xyyy_0[j] = pa_x[j] * tpx_y_xyyy_0[j] + 0.5 * fl1_fx * tpx_y_yyy_0[j] - fl1_fgb * fl1_fx * ts_y_xyyy_0[j];

            tpy_xy_xyyy_0[j] = pa_x[j] * tpy_y_xyyy_0[j] + 0.5 * fl1_fx * tpy_y_yyy_0[j];

            tpz_xy_xyyy_0[j] = pa_x[j] * tpz_y_xyyy_0[j] + 0.5 * fl1_fx * tpz_y_yyy_0[j];

            tpx_xy_xyyz_0[j] = pa_x[j] * tpx_y_xyyz_0[j] + 0.5 * fl1_fx * tpx_y_yyz_0[j] - fl1_fgb * fl1_fx * ts_y_xyyz_0[j];

            tpy_xy_xyyz_0[j] = pa_x[j] * tpy_y_xyyz_0[j] + 0.5 * fl1_fx * tpy_y_yyz_0[j];

            tpz_xy_xyyz_0[j] = pa_x[j] * tpz_y_xyyz_0[j] + 0.5 * fl1_fx * tpz_y_yyz_0[j];

            tpx_xy_xyzz_0[j] = pa_x[j] * tpx_y_xyzz_0[j] + 0.5 * fl1_fx * tpx_y_yzz_0[j] - fl1_fgb * fl1_fx * ts_y_xyzz_0[j];

            tpy_xy_xyzz_0[j] = pa_x[j] * tpy_y_xyzz_0[j] + 0.5 * fl1_fx * tpy_y_yzz_0[j];

            tpz_xy_xyzz_0[j] = pa_x[j] * tpz_y_xyzz_0[j] + 0.5 * fl1_fx * tpz_y_yzz_0[j];

            tpx_xy_xzzz_0[j] = pa_x[j] * tpx_y_xzzz_0[j] + 0.5 * fl1_fx * tpx_y_zzz_0[j] - fl1_fgb * fl1_fx * ts_y_xzzz_0[j];

            tpy_xy_xzzz_0[j] = pa_x[j] * tpy_y_xzzz_0[j] + 0.5 * fl1_fx * tpy_y_zzz_0[j];

            tpz_xy_xzzz_0[j] = pa_x[j] * tpz_y_xzzz_0[j] + 0.5 * fl1_fx * tpz_y_zzz_0[j];

            tpx_xy_yyyy_0[j] = pa_x[j] * tpx_y_yyyy_0[j] - fl1_fgb * fl1_fx * ts_y_yyyy_0[j];

            tpy_xy_yyyy_0[j] = pa_x[j] * tpy_y_yyyy_0[j];

            tpz_xy_yyyy_0[j] = pa_x[j] * tpz_y_yyyy_0[j];

            tpx_xy_yyyz_0[j] = pa_x[j] * tpx_y_yyyz_0[j] - fl1_fgb * fl1_fx * ts_y_yyyz_0[j];

            tpy_xy_yyyz_0[j] = pa_x[j] * tpy_y_yyyz_0[j];

            tpz_xy_yyyz_0[j] = pa_x[j] * tpz_y_yyyz_0[j];

            tpx_xy_yyzz_0[j] = pa_x[j] * tpx_y_yyzz_0[j] - fl1_fgb * fl1_fx * ts_y_yyzz_0[j];

            tpy_xy_yyzz_0[j] = pa_x[j] * tpy_y_yyzz_0[j];

            tpz_xy_yyzz_0[j] = pa_x[j] * tpz_y_yyzz_0[j];

            tpx_xy_yzzz_0[j] = pa_x[j] * tpx_y_yzzz_0[j] - fl1_fgb * fl1_fx * ts_y_yzzz_0[j];

            tpy_xy_yzzz_0[j] = pa_x[j] * tpy_y_yzzz_0[j];

            tpz_xy_yzzz_0[j] = pa_x[j] * tpz_y_yzzz_0[j];

            tpx_xy_zzzz_0[j] = pa_x[j] * tpx_y_zzzz_0[j] - fl1_fgb * fl1_fx * ts_y_zzzz_0[j];

            tpy_xy_zzzz_0[j] = pa_x[j] * tpy_y_zzzz_0[j];

            tpz_xy_zzzz_0[j] = pa_x[j] * tpz_y_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG_90_135(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

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

        // set up pointers to integrals

        auto tpx_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 30);

        auto tpy_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tpz_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tpx_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 31);

        auto tpy_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tpz_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tpx_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 32);

        auto tpy_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tpz_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tpx_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 33);

        auto tpy_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tpz_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tpx_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 34);

        auto tpy_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tpz_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tpx_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 35);

        auto tpy_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tpz_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tpx_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 36);

        auto tpy_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tpz_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tpx_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 37);

        auto tpy_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tpz_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tpx_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 38);

        auto tpy_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tpz_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tpx_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 39);

        auto tpy_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tpz_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tpx_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 40);

        auto tpy_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tpz_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tpx_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 41);

        auto tpy_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tpz_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tpx_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 42);

        auto tpy_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tpz_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tpx_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 43);

        auto tpy_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tpz_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tpx_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 44);

        auto tpy_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tpz_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xz_xxxx_0, tpx_xz_xxxy_0, tpx_xz_xxxz_0, \
                                     tpx_xz_xxyy_0, tpx_xz_xxyz_0, tpx_xz_xxzz_0, tpx_xz_xyyy_0, tpx_xz_xyyz_0, \
                                     tpx_xz_xyzz_0, tpx_xz_xzzz_0, tpx_xz_yyyy_0, tpx_xz_yyyz_0, tpx_xz_yyzz_0, \
                                     tpx_xz_yzzz_0, tpx_xz_zzzz_0, tpx_z_xxx_0, tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, \
                                     tpx_z_xxy_0, tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxz_0, tpx_z_xxzz_0, tpx_z_xyy_0, \
                                     tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyz_0, tpx_z_xyzz_0, tpx_z_xzz_0, tpx_z_xzzz_0, \
                                     tpx_z_yyy_0, tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyz_0, tpx_z_yyzz_0, tpx_z_yzz_0, \
                                     tpx_z_yzzz_0, tpx_z_zzz_0, tpx_z_zzzz_0, tpy_xz_xxxx_0, tpy_xz_xxxy_0, \
                                     tpy_xz_xxxz_0, tpy_xz_xxyy_0, tpy_xz_xxyz_0, tpy_xz_xxzz_0, tpy_xz_xyyy_0, \
                                     tpy_xz_xyyz_0, tpy_xz_xyzz_0, tpy_xz_xzzz_0, tpy_xz_yyyy_0, tpy_xz_yyyz_0, \
                                     tpy_xz_yyzz_0, tpy_xz_yzzz_0, tpy_xz_zzzz_0, tpy_z_xxx_0, tpy_z_xxxx_0, \
                                     tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxy_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxz_0, \
                                     tpy_z_xxzz_0, tpy_z_xyy_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyz_0, tpy_z_xyzz_0, \
                                     tpy_z_xzz_0, tpy_z_xzzz_0, tpy_z_yyy_0, tpy_z_yyyy_0, tpy_z_yyyz_0, tpy_z_yyz_0, \
                                     tpy_z_yyzz_0, tpy_z_yzz_0, tpy_z_yzzz_0, tpy_z_zzz_0, tpy_z_zzzz_0, tpz_xz_xxxx_0, \
                                     tpz_xz_xxxy_0, tpz_xz_xxxz_0, tpz_xz_xxyy_0, tpz_xz_xxyz_0, tpz_xz_xxzz_0, \
                                     tpz_xz_xyyy_0, tpz_xz_xyyz_0, tpz_xz_xyzz_0, tpz_xz_xzzz_0, tpz_xz_yyyy_0, \
                                     tpz_xz_yyyz_0, tpz_xz_yyzz_0, tpz_xz_yzzz_0, tpz_xz_zzzz_0, tpz_z_xxx_0, \
                                     tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, tpz_z_xxy_0, tpz_z_xxyy_0, tpz_z_xxyz_0, \
                                     tpz_z_xxz_0, tpz_z_xxzz_0, tpz_z_xyy_0, tpz_z_xyyy_0, tpz_z_xyyz_0, tpz_z_xyz_0, \
                                     tpz_z_xyzz_0, tpz_z_xzz_0, tpz_z_xzzz_0, tpz_z_yyy_0, tpz_z_yyyy_0, tpz_z_yyyz_0, \
                                     tpz_z_yyz_0, tpz_z_yyzz_0, tpz_z_yzz_0, tpz_z_yzzz_0, tpz_z_zzz_0, tpz_z_zzzz_0, \
                                     ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, \
                                     ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, \
                                     ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xz_xxxx_0[j] = pa_x[j] * tpx_z_xxxx_0[j] + 2.0 * fl1_fx * tpx_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxxx_0[j];

            tpy_xz_xxxx_0[j] = pa_x[j] * tpy_z_xxxx_0[j] + 2.0 * fl1_fx * tpy_z_xxx_0[j];

            tpz_xz_xxxx_0[j] = pa_x[j] * tpz_z_xxxx_0[j] + 2.0 * fl1_fx * tpz_z_xxx_0[j];

            tpx_xz_xxxy_0[j] = pa_x[j] * tpx_z_xxxy_0[j] + 1.5 * fl1_fx * tpx_z_xxy_0[j] - fl1_fgb * fl1_fx * ts_z_xxxy_0[j];

            tpy_xz_xxxy_0[j] = pa_x[j] * tpy_z_xxxy_0[j] + 1.5 * fl1_fx * tpy_z_xxy_0[j];

            tpz_xz_xxxy_0[j] = pa_x[j] * tpz_z_xxxy_0[j] + 1.5 * fl1_fx * tpz_z_xxy_0[j];

            tpx_xz_xxxz_0[j] = pa_x[j] * tpx_z_xxxz_0[j] + 1.5 * fl1_fx * tpx_z_xxz_0[j] - fl1_fgb * fl1_fx * ts_z_xxxz_0[j];

            tpy_xz_xxxz_0[j] = pa_x[j] * tpy_z_xxxz_0[j] + 1.5 * fl1_fx * tpy_z_xxz_0[j];

            tpz_xz_xxxz_0[j] = pa_x[j] * tpz_z_xxxz_0[j] + 1.5 * fl1_fx * tpz_z_xxz_0[j];

            tpx_xz_xxyy_0[j] = pa_x[j] * tpx_z_xxyy_0[j] + fl1_fx * tpx_z_xyy_0[j] - fl1_fgb * fl1_fx * ts_z_xxyy_0[j];

            tpy_xz_xxyy_0[j] = pa_x[j] * tpy_z_xxyy_0[j] + fl1_fx * tpy_z_xyy_0[j];

            tpz_xz_xxyy_0[j] = pa_x[j] * tpz_z_xxyy_0[j] + fl1_fx * tpz_z_xyy_0[j];

            tpx_xz_xxyz_0[j] = pa_x[j] * tpx_z_xxyz_0[j] + fl1_fx * tpx_z_xyz_0[j] - fl1_fgb * fl1_fx * ts_z_xxyz_0[j];

            tpy_xz_xxyz_0[j] = pa_x[j] * tpy_z_xxyz_0[j] + fl1_fx * tpy_z_xyz_0[j];

            tpz_xz_xxyz_0[j] = pa_x[j] * tpz_z_xxyz_0[j] + fl1_fx * tpz_z_xyz_0[j];

            tpx_xz_xxzz_0[j] = pa_x[j] * tpx_z_xxzz_0[j] + fl1_fx * tpx_z_xzz_0[j] - fl1_fgb * fl1_fx * ts_z_xxzz_0[j];

            tpy_xz_xxzz_0[j] = pa_x[j] * tpy_z_xxzz_0[j] + fl1_fx * tpy_z_xzz_0[j];

            tpz_xz_xxzz_0[j] = pa_x[j] * tpz_z_xxzz_0[j] + fl1_fx * tpz_z_xzz_0[j];

            tpx_xz_xyyy_0[j] = pa_x[j] * tpx_z_xyyy_0[j] + 0.5 * fl1_fx * tpx_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_z_xyyy_0[j];

            tpy_xz_xyyy_0[j] = pa_x[j] * tpy_z_xyyy_0[j] + 0.5 * fl1_fx * tpy_z_yyy_0[j];

            tpz_xz_xyyy_0[j] = pa_x[j] * tpz_z_xyyy_0[j] + 0.5 * fl1_fx * tpz_z_yyy_0[j];

            tpx_xz_xyyz_0[j] = pa_x[j] * tpx_z_xyyz_0[j] + 0.5 * fl1_fx * tpx_z_yyz_0[j] - fl1_fgb * fl1_fx * ts_z_xyyz_0[j];

            tpy_xz_xyyz_0[j] = pa_x[j] * tpy_z_xyyz_0[j] + 0.5 * fl1_fx * tpy_z_yyz_0[j];

            tpz_xz_xyyz_0[j] = pa_x[j] * tpz_z_xyyz_0[j] + 0.5 * fl1_fx * tpz_z_yyz_0[j];

            tpx_xz_xyzz_0[j] = pa_x[j] * tpx_z_xyzz_0[j] + 0.5 * fl1_fx * tpx_z_yzz_0[j] - fl1_fgb * fl1_fx * ts_z_xyzz_0[j];

            tpy_xz_xyzz_0[j] = pa_x[j] * tpy_z_xyzz_0[j] + 0.5 * fl1_fx * tpy_z_yzz_0[j];

            tpz_xz_xyzz_0[j] = pa_x[j] * tpz_z_xyzz_0[j] + 0.5 * fl1_fx * tpz_z_yzz_0[j];

            tpx_xz_xzzz_0[j] = pa_x[j] * tpx_z_xzzz_0[j] + 0.5 * fl1_fx * tpx_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_z_xzzz_0[j];

            tpy_xz_xzzz_0[j] = pa_x[j] * tpy_z_xzzz_0[j] + 0.5 * fl1_fx * tpy_z_zzz_0[j];

            tpz_xz_xzzz_0[j] = pa_x[j] * tpz_z_xzzz_0[j] + 0.5 * fl1_fx * tpz_z_zzz_0[j];

            tpx_xz_yyyy_0[j] = pa_x[j] * tpx_z_yyyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyyy_0[j];

            tpy_xz_yyyy_0[j] = pa_x[j] * tpy_z_yyyy_0[j];

            tpz_xz_yyyy_0[j] = pa_x[j] * tpz_z_yyyy_0[j];

            tpx_xz_yyyz_0[j] = pa_x[j] * tpx_z_yyyz_0[j] - fl1_fgb * fl1_fx * ts_z_yyyz_0[j];

            tpy_xz_yyyz_0[j] = pa_x[j] * tpy_z_yyyz_0[j];

            tpz_xz_yyyz_0[j] = pa_x[j] * tpz_z_yyyz_0[j];

            tpx_xz_yyzz_0[j] = pa_x[j] * tpx_z_yyzz_0[j] - fl1_fgb * fl1_fx * ts_z_yyzz_0[j];

            tpy_xz_yyzz_0[j] = pa_x[j] * tpy_z_yyzz_0[j];

            tpz_xz_yyzz_0[j] = pa_x[j] * tpz_z_yyzz_0[j];

            tpx_xz_yzzz_0[j] = pa_x[j] * tpx_z_yzzz_0[j] - fl1_fgb * fl1_fx * ts_z_yzzz_0[j];

            tpy_xz_yzzz_0[j] = pa_x[j] * tpy_z_yzzz_0[j];

            tpz_xz_yzzz_0[j] = pa_x[j] * tpz_z_yzzz_0[j];

            tpx_xz_zzzz_0[j] = pa_x[j] * tpx_z_zzzz_0[j] - fl1_fgb * fl1_fx * ts_z_zzzz_0[j];

            tpy_xz_zzzz_0[j] = pa_x[j] * tpy_z_zzzz_0[j];

            tpz_xz_zzzz_0[j] = pa_x[j] * tpz_z_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG_135_180(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (135,180)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15);

        auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16);

        auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17);

        auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18);

        auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19);

        auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20);

        auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21);

        auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22);

        auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23);

        auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24);

        auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25);

        auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26);

        auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27);

        auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28);

        auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29);

        auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

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

        // set up pointers to integrals

        auto tpx_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 45);

        auto tpy_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tpz_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tpx_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 46);

        auto tpy_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tpz_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tpx_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 47);

        auto tpy_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tpz_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tpx_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 48);

        auto tpy_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tpz_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tpx_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 49);

        auto tpy_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tpx_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 50);

        auto tpy_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tpz_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tpx_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 51);

        auto tpy_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tpz_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tpx_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 52);

        auto tpy_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tpz_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tpx_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 53);

        auto tpy_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tpz_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tpx_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 54);

        auto tpy_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tpz_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tpx_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 55);

        auto tpy_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tpz_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tpx_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 56);

        auto tpy_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tpz_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tpx_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 57);

        auto tpy_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tpz_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tpx_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 58);

        auto tpy_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tpz_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tpx_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 59);

        auto tpy_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tpz_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 59);

        // Batch of Integrals (135,180)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxyy_0, \
                                     tpx_0_xxyz_0, tpx_0_xxzz_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyzz_0, tpx_0_yzzz_0, tpx_0_zzzz_0, tpx_y_xxx_0, \
                                     tpx_y_xxxx_0, tpx_y_xxxy_0, tpx_y_xxxz_0, tpx_y_xxy_0, tpx_y_xxyy_0, tpx_y_xxyz_0, \
                                     tpx_y_xxz_0, tpx_y_xxzz_0, tpx_y_xyy_0, tpx_y_xyyy_0, tpx_y_xyyz_0, tpx_y_xyz_0, \
                                     tpx_y_xyzz_0, tpx_y_xzz_0, tpx_y_xzzz_0, tpx_y_yyy_0, tpx_y_yyyy_0, tpx_y_yyyz_0, \
                                     tpx_y_yyz_0, tpx_y_yyzz_0, tpx_y_yzz_0, tpx_y_yzzz_0, tpx_y_zzz_0, tpx_y_zzzz_0, \
                                     tpx_yy_xxxx_0, tpx_yy_xxxy_0, tpx_yy_xxxz_0, tpx_yy_xxyy_0, tpx_yy_xxyz_0, \
                                     tpx_yy_xxzz_0, tpx_yy_xyyy_0, tpx_yy_xyyz_0, tpx_yy_xyzz_0, tpx_yy_xzzz_0, \
                                     tpx_yy_yyyy_0, tpx_yy_yyyz_0, tpx_yy_yyzz_0, tpx_yy_yzzz_0, tpx_yy_zzzz_0, \
                                     tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxyy_0, tpy_0_xxyz_0, tpy_0_xxzz_0, \
                                     tpy_0_xyyy_0, tpy_0_xyyz_0, tpy_0_xyzz_0, tpy_0_xzzz_0, tpy_0_yyyy_0, tpy_0_yyyz_0, \
                                     tpy_0_yyzz_0, tpy_0_yzzz_0, tpy_0_zzzz_0, tpy_y_xxx_0, tpy_y_xxxx_0, tpy_y_xxxy_0, \
                                     tpy_y_xxxz_0, tpy_y_xxy_0, tpy_y_xxyy_0, tpy_y_xxyz_0, tpy_y_xxz_0, tpy_y_xxzz_0, \
                                     tpy_y_xyy_0, tpy_y_xyyy_0, tpy_y_xyyz_0, tpy_y_xyz_0, tpy_y_xyzz_0, tpy_y_xzz_0, \
                                     tpy_y_xzzz_0, tpy_y_yyy_0, tpy_y_yyyy_0, tpy_y_yyyz_0, tpy_y_yyz_0, tpy_y_yyzz_0, \
                                     tpy_y_yzz_0, tpy_y_yzzz_0, tpy_y_zzz_0, tpy_y_zzzz_0, tpy_yy_xxxx_0, \
                                     tpy_yy_xxxy_0, tpy_yy_xxxz_0, tpy_yy_xxyy_0, tpy_yy_xxyz_0, tpy_yy_xxzz_0, \
                                     tpy_yy_xyyy_0, tpy_yy_xyyz_0, tpy_yy_xyzz_0, tpy_yy_xzzz_0, tpy_yy_yyyy_0, \
                                     tpy_yy_yyyz_0, tpy_yy_yyzz_0, tpy_yy_yzzz_0, tpy_yy_zzzz_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxzz_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyzz_0, tpz_0_xzzz_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yzzz_0, tpz_0_zzzz_0, tpz_y_xxx_0, tpz_y_xxxx_0, tpz_y_xxxy_0, tpz_y_xxxz_0, \
                                     tpz_y_xxy_0, tpz_y_xxyy_0, tpz_y_xxyz_0, tpz_y_xxz_0, tpz_y_xxzz_0, tpz_y_xyy_0, \
                                     tpz_y_xyyy_0, tpz_y_xyyz_0, tpz_y_xyz_0, tpz_y_xyzz_0, tpz_y_xzz_0, tpz_y_xzzz_0, \
                                     tpz_y_yyy_0, tpz_y_yyyy_0, tpz_y_yyyz_0, tpz_y_yyz_0, tpz_y_yyzz_0, tpz_y_yzz_0, \
                                     tpz_y_yzzz_0, tpz_y_zzz_0, tpz_y_zzzz_0, tpz_yy_xxxx_0, tpz_yy_xxxy_0, \
                                     tpz_yy_xxxz_0, tpz_yy_xxyy_0, tpz_yy_xxyz_0, tpz_yy_xxzz_0, tpz_yy_xyyy_0, \
                                     tpz_yy_xyyz_0, tpz_yy_xyzz_0, tpz_yy_xzzz_0, tpz_yy_yyyy_0, tpz_yy_yyyz_0, \
                                     tpz_yy_yyzz_0, tpz_yy_yzzz_0, tpz_yy_zzzz_0, ts_y_xxxx_0, ts_y_xxxy_0, ts_y_xxxz_0, \
                                     ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, ts_y_xyyz_0, ts_y_xyzz_0, \
                                     ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, ts_y_yzzz_0, ts_y_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yy_xxxx_0[j] = pa_y[j] * tpx_y_xxxx_0[j] + 0.5 * fl1_fx * tpx_0_xxxx_0[j];

            tpy_yy_xxxx_0[j] = pa_y[j] * tpy_y_xxxx_0[j] + 0.5 * fl1_fx * tpy_0_xxxx_0[j] - fl1_fgb * fl1_fx * ts_y_xxxx_0[j];

            tpz_yy_xxxx_0[j] = pa_y[j] * tpz_y_xxxx_0[j] + 0.5 * fl1_fx * tpz_0_xxxx_0[j];

            tpx_yy_xxxy_0[j] = pa_y[j] * tpx_y_xxxy_0[j] + 0.5 * fl1_fx * tpx_0_xxxy_0[j] + 0.5 * fl1_fx * tpx_y_xxx_0[j];

            tpy_yy_xxxy_0[j] =
                pa_y[j] * tpy_y_xxxy_0[j] + 0.5 * fl1_fx * tpy_0_xxxy_0[j] + 0.5 * fl1_fx * tpy_y_xxx_0[j] - fl1_fgb * fl1_fx * ts_y_xxxy_0[j];

            tpz_yy_xxxy_0[j] = pa_y[j] * tpz_y_xxxy_0[j] + 0.5 * fl1_fx * tpz_0_xxxy_0[j] + 0.5 * fl1_fx * tpz_y_xxx_0[j];

            tpx_yy_xxxz_0[j] = pa_y[j] * tpx_y_xxxz_0[j] + 0.5 * fl1_fx * tpx_0_xxxz_0[j];

            tpy_yy_xxxz_0[j] = pa_y[j] * tpy_y_xxxz_0[j] + 0.5 * fl1_fx * tpy_0_xxxz_0[j] - fl1_fgb * fl1_fx * ts_y_xxxz_0[j];

            tpz_yy_xxxz_0[j] = pa_y[j] * tpz_y_xxxz_0[j] + 0.5 * fl1_fx * tpz_0_xxxz_0[j];

            tpx_yy_xxyy_0[j] = pa_y[j] * tpx_y_xxyy_0[j] + 0.5 * fl1_fx * tpx_0_xxyy_0[j] + fl1_fx * tpx_y_xxy_0[j];

            tpy_yy_xxyy_0[j] =
                pa_y[j] * tpy_y_xxyy_0[j] + 0.5 * fl1_fx * tpy_0_xxyy_0[j] + fl1_fx * tpy_y_xxy_0[j] - fl1_fgb * fl1_fx * ts_y_xxyy_0[j];

            tpz_yy_xxyy_0[j] = pa_y[j] * tpz_y_xxyy_0[j] + 0.5 * fl1_fx * tpz_0_xxyy_0[j] + fl1_fx * tpz_y_xxy_0[j];

            tpx_yy_xxyz_0[j] = pa_y[j] * tpx_y_xxyz_0[j] + 0.5 * fl1_fx * tpx_0_xxyz_0[j] + 0.5 * fl1_fx * tpx_y_xxz_0[j];

            tpy_yy_xxyz_0[j] =
                pa_y[j] * tpy_y_xxyz_0[j] + 0.5 * fl1_fx * tpy_0_xxyz_0[j] + 0.5 * fl1_fx * tpy_y_xxz_0[j] - fl1_fgb * fl1_fx * ts_y_xxyz_0[j];

            tpz_yy_xxyz_0[j] = pa_y[j] * tpz_y_xxyz_0[j] + 0.5 * fl1_fx * tpz_0_xxyz_0[j] + 0.5 * fl1_fx * tpz_y_xxz_0[j];

            tpx_yy_xxzz_0[j] = pa_y[j] * tpx_y_xxzz_0[j] + 0.5 * fl1_fx * tpx_0_xxzz_0[j];

            tpy_yy_xxzz_0[j] = pa_y[j] * tpy_y_xxzz_0[j] + 0.5 * fl1_fx * tpy_0_xxzz_0[j] - fl1_fgb * fl1_fx * ts_y_xxzz_0[j];

            tpz_yy_xxzz_0[j] = pa_y[j] * tpz_y_xxzz_0[j] + 0.5 * fl1_fx * tpz_0_xxzz_0[j];

            tpx_yy_xyyy_0[j] = pa_y[j] * tpx_y_xyyy_0[j] + 0.5 * fl1_fx * tpx_0_xyyy_0[j] + 1.5 * fl1_fx * tpx_y_xyy_0[j];

            tpy_yy_xyyy_0[j] =
                pa_y[j] * tpy_y_xyyy_0[j] + 0.5 * fl1_fx * tpy_0_xyyy_0[j] + 1.5 * fl1_fx * tpy_y_xyy_0[j] - fl1_fgb * fl1_fx * ts_y_xyyy_0[j];

            tpz_yy_xyyy_0[j] = pa_y[j] * tpz_y_xyyy_0[j] + 0.5 * fl1_fx * tpz_0_xyyy_0[j] + 1.5 * fl1_fx * tpz_y_xyy_0[j];

            tpx_yy_xyyz_0[j] = pa_y[j] * tpx_y_xyyz_0[j] + 0.5 * fl1_fx * tpx_0_xyyz_0[j] + fl1_fx * tpx_y_xyz_0[j];

            tpy_yy_xyyz_0[j] =
                pa_y[j] * tpy_y_xyyz_0[j] + 0.5 * fl1_fx * tpy_0_xyyz_0[j] + fl1_fx * tpy_y_xyz_0[j] - fl1_fgb * fl1_fx * ts_y_xyyz_0[j];

            tpz_yy_xyyz_0[j] = pa_y[j] * tpz_y_xyyz_0[j] + 0.5 * fl1_fx * tpz_0_xyyz_0[j] + fl1_fx * tpz_y_xyz_0[j];

            tpx_yy_xyzz_0[j] = pa_y[j] * tpx_y_xyzz_0[j] + 0.5 * fl1_fx * tpx_0_xyzz_0[j] + 0.5 * fl1_fx * tpx_y_xzz_0[j];

            tpy_yy_xyzz_0[j] =
                pa_y[j] * tpy_y_xyzz_0[j] + 0.5 * fl1_fx * tpy_0_xyzz_0[j] + 0.5 * fl1_fx * tpy_y_xzz_0[j] - fl1_fgb * fl1_fx * ts_y_xyzz_0[j];

            tpz_yy_xyzz_0[j] = pa_y[j] * tpz_y_xyzz_0[j] + 0.5 * fl1_fx * tpz_0_xyzz_0[j] + 0.5 * fl1_fx * tpz_y_xzz_0[j];

            tpx_yy_xzzz_0[j] = pa_y[j] * tpx_y_xzzz_0[j] + 0.5 * fl1_fx * tpx_0_xzzz_0[j];

            tpy_yy_xzzz_0[j] = pa_y[j] * tpy_y_xzzz_0[j] + 0.5 * fl1_fx * tpy_0_xzzz_0[j] - fl1_fgb * fl1_fx * ts_y_xzzz_0[j];

            tpz_yy_xzzz_0[j] = pa_y[j] * tpz_y_xzzz_0[j] + 0.5 * fl1_fx * tpz_0_xzzz_0[j];

            tpx_yy_yyyy_0[j] = pa_y[j] * tpx_y_yyyy_0[j] + 0.5 * fl1_fx * tpx_0_yyyy_0[j] + 2.0 * fl1_fx * tpx_y_yyy_0[j];

            tpy_yy_yyyy_0[j] =
                pa_y[j] * tpy_y_yyyy_0[j] + 0.5 * fl1_fx * tpy_0_yyyy_0[j] + 2.0 * fl1_fx * tpy_y_yyy_0[j] - fl1_fgb * fl1_fx * ts_y_yyyy_0[j];

            tpz_yy_yyyy_0[j] = pa_y[j] * tpz_y_yyyy_0[j] + 0.5 * fl1_fx * tpz_0_yyyy_0[j] + 2.0 * fl1_fx * tpz_y_yyy_0[j];

            tpx_yy_yyyz_0[j] = pa_y[j] * tpx_y_yyyz_0[j] + 0.5 * fl1_fx * tpx_0_yyyz_0[j] + 1.5 * fl1_fx * tpx_y_yyz_0[j];

            tpy_yy_yyyz_0[j] =
                pa_y[j] * tpy_y_yyyz_0[j] + 0.5 * fl1_fx * tpy_0_yyyz_0[j] + 1.5 * fl1_fx * tpy_y_yyz_0[j] - fl1_fgb * fl1_fx * ts_y_yyyz_0[j];

            tpz_yy_yyyz_0[j] = pa_y[j] * tpz_y_yyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyyz_0[j] + 1.5 * fl1_fx * tpz_y_yyz_0[j];

            tpx_yy_yyzz_0[j] = pa_y[j] * tpx_y_yyzz_0[j] + 0.5 * fl1_fx * tpx_0_yyzz_0[j] + fl1_fx * tpx_y_yzz_0[j];

            tpy_yy_yyzz_0[j] =
                pa_y[j] * tpy_y_yyzz_0[j] + 0.5 * fl1_fx * tpy_0_yyzz_0[j] + fl1_fx * tpy_y_yzz_0[j] - fl1_fgb * fl1_fx * ts_y_yyzz_0[j];

            tpz_yy_yyzz_0[j] = pa_y[j] * tpz_y_yyzz_0[j] + 0.5 * fl1_fx * tpz_0_yyzz_0[j] + fl1_fx * tpz_y_yzz_0[j];

            tpx_yy_yzzz_0[j] = pa_y[j] * tpx_y_yzzz_0[j] + 0.5 * fl1_fx * tpx_0_yzzz_0[j] + 0.5 * fl1_fx * tpx_y_zzz_0[j];

            tpy_yy_yzzz_0[j] =
                pa_y[j] * tpy_y_yzzz_0[j] + 0.5 * fl1_fx * tpy_0_yzzz_0[j] + 0.5 * fl1_fx * tpy_y_zzz_0[j] - fl1_fgb * fl1_fx * ts_y_yzzz_0[j];

            tpz_yy_yzzz_0[j] = pa_y[j] * tpz_y_yzzz_0[j] + 0.5 * fl1_fx * tpz_0_yzzz_0[j] + 0.5 * fl1_fx * tpz_y_zzz_0[j];

            tpx_yy_zzzz_0[j] = pa_y[j] * tpx_y_zzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzzz_0[j];

            tpy_yy_zzzz_0[j] = pa_y[j] * tpy_y_zzzz_0[j] + 0.5 * fl1_fx * tpy_0_zzzz_0[j] - fl1_fgb * fl1_fx * ts_y_zzzz_0[j];

            tpz_yy_zzzz_0[j] = pa_y[j] * tpz_y_zzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG_180_225(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (180,225)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

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

        // set up pointers to integrals

        auto tpx_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 60);

        auto tpy_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tpz_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tpx_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 61);

        auto tpy_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tpz_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tpx_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 62);

        auto tpy_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tpz_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tpx_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 63);

        auto tpy_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tpz_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tpx_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 64);

        auto tpy_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tpz_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tpx_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 65);

        auto tpy_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 66);

        auto tpy_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tpz_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tpx_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 67);

        auto tpy_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tpz_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tpx_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 68);

        auto tpy_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tpz_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tpx_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 69);

        auto tpy_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tpz_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tpx_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 70);

        auto tpy_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tpz_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tpx_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 71);

        auto tpy_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tpz_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tpx_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 72);

        auto tpy_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tpz_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tpx_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 73);

        auto tpy_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tpz_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tpx_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 74);

        auto tpy_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tpz_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 74);

        // Batch of Integrals (180,225)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yz_xxxx_0, tpx_yz_xxxy_0, tpx_yz_xxxz_0, \
                                     tpx_yz_xxyy_0, tpx_yz_xxyz_0, tpx_yz_xxzz_0, tpx_yz_xyyy_0, tpx_yz_xyyz_0, \
                                     tpx_yz_xyzz_0, tpx_yz_xzzz_0, tpx_yz_yyyy_0, tpx_yz_yyyz_0, tpx_yz_yyzz_0, \
                                     tpx_yz_yzzz_0, tpx_yz_zzzz_0, tpx_z_xxx_0, tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, \
                                     tpx_z_xxy_0, tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxz_0, tpx_z_xxzz_0, tpx_z_xyy_0, \
                                     tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyz_0, tpx_z_xyzz_0, tpx_z_xzz_0, tpx_z_xzzz_0, \
                                     tpx_z_yyy_0, tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyz_0, tpx_z_yyzz_0, tpx_z_yzz_0, \
                                     tpx_z_yzzz_0, tpx_z_zzz_0, tpx_z_zzzz_0, tpy_yz_xxxx_0, tpy_yz_xxxy_0, \
                                     tpy_yz_xxxz_0, tpy_yz_xxyy_0, tpy_yz_xxyz_0, tpy_yz_xxzz_0, tpy_yz_xyyy_0, \
                                     tpy_yz_xyyz_0, tpy_yz_xyzz_0, tpy_yz_xzzz_0, tpy_yz_yyyy_0, tpy_yz_yyyz_0, \
                                     tpy_yz_yyzz_0, tpy_yz_yzzz_0, tpy_yz_zzzz_0, tpy_z_xxx_0, tpy_z_xxxx_0, \
                                     tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxy_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxz_0, \
                                     tpy_z_xxzz_0, tpy_z_xyy_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyz_0, tpy_z_xyzz_0, \
                                     tpy_z_xzz_0, tpy_z_xzzz_0, tpy_z_yyy_0, tpy_z_yyyy_0, tpy_z_yyyz_0, tpy_z_yyz_0, \
                                     tpy_z_yyzz_0, tpy_z_yzz_0, tpy_z_yzzz_0, tpy_z_zzz_0, tpy_z_zzzz_0, tpz_yz_xxxx_0, \
                                     tpz_yz_xxxy_0, tpz_yz_xxxz_0, tpz_yz_xxyy_0, tpz_yz_xxyz_0, tpz_yz_xxzz_0, \
                                     tpz_yz_xyyy_0, tpz_yz_xyyz_0, tpz_yz_xyzz_0, tpz_yz_xzzz_0, tpz_yz_yyyy_0, \
                                     tpz_yz_yyyz_0, tpz_yz_yyzz_0, tpz_yz_yzzz_0, tpz_yz_zzzz_0, tpz_z_xxx_0, \
                                     tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, tpz_z_xxy_0, tpz_z_xxyy_0, tpz_z_xxyz_0, \
                                     tpz_z_xxz_0, tpz_z_xxzz_0, tpz_z_xyy_0, tpz_z_xyyy_0, tpz_z_xyyz_0, tpz_z_xyz_0, \
                                     tpz_z_xyzz_0, tpz_z_xzz_0, tpz_z_xzzz_0, tpz_z_yyy_0, tpz_z_yyyy_0, tpz_z_yyyz_0, \
                                     tpz_z_yyz_0, tpz_z_yyzz_0, tpz_z_yzz_0, tpz_z_yzzz_0, tpz_z_zzz_0, tpz_z_zzzz_0, \
                                     ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, \
                                     ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, \
                                     ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yz_xxxx_0[j] = pa_y[j] * tpx_z_xxxx_0[j];

            tpy_yz_xxxx_0[j] = pa_y[j] * tpy_z_xxxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxxx_0[j];

            tpz_yz_xxxx_0[j] = pa_y[j] * tpz_z_xxxx_0[j];

            tpx_yz_xxxy_0[j] = pa_y[j] * tpx_z_xxxy_0[j] + 0.5 * fl1_fx * tpx_z_xxx_0[j];

            tpy_yz_xxxy_0[j] = pa_y[j] * tpy_z_xxxy_0[j] + 0.5 * fl1_fx * tpy_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxxy_0[j];

            tpz_yz_xxxy_0[j] = pa_y[j] * tpz_z_xxxy_0[j] + 0.5 * fl1_fx * tpz_z_xxx_0[j];

            tpx_yz_xxxz_0[j] = pa_y[j] * tpx_z_xxxz_0[j];

            tpy_yz_xxxz_0[j] = pa_y[j] * tpy_z_xxxz_0[j] - fl1_fgb * fl1_fx * ts_z_xxxz_0[j];

            tpz_yz_xxxz_0[j] = pa_y[j] * tpz_z_xxxz_0[j];

            tpx_yz_xxyy_0[j] = pa_y[j] * tpx_z_xxyy_0[j] + fl1_fx * tpx_z_xxy_0[j];

            tpy_yz_xxyy_0[j] = pa_y[j] * tpy_z_xxyy_0[j] + fl1_fx * tpy_z_xxy_0[j] - fl1_fgb * fl1_fx * ts_z_xxyy_0[j];

            tpz_yz_xxyy_0[j] = pa_y[j] * tpz_z_xxyy_0[j] + fl1_fx * tpz_z_xxy_0[j];

            tpx_yz_xxyz_0[j] = pa_y[j] * tpx_z_xxyz_0[j] + 0.5 * fl1_fx * tpx_z_xxz_0[j];

            tpy_yz_xxyz_0[j] = pa_y[j] * tpy_z_xxyz_0[j] + 0.5 * fl1_fx * tpy_z_xxz_0[j] - fl1_fgb * fl1_fx * ts_z_xxyz_0[j];

            tpz_yz_xxyz_0[j] = pa_y[j] * tpz_z_xxyz_0[j] + 0.5 * fl1_fx * tpz_z_xxz_0[j];

            tpx_yz_xxzz_0[j] = pa_y[j] * tpx_z_xxzz_0[j];

            tpy_yz_xxzz_0[j] = pa_y[j] * tpy_z_xxzz_0[j] - fl1_fgb * fl1_fx * ts_z_xxzz_0[j];

            tpz_yz_xxzz_0[j] = pa_y[j] * tpz_z_xxzz_0[j];

            tpx_yz_xyyy_0[j] = pa_y[j] * tpx_z_xyyy_0[j] + 1.5 * fl1_fx * tpx_z_xyy_0[j];

            tpy_yz_xyyy_0[j] = pa_y[j] * tpy_z_xyyy_0[j] + 1.5 * fl1_fx * tpy_z_xyy_0[j] - fl1_fgb * fl1_fx * ts_z_xyyy_0[j];

            tpz_yz_xyyy_0[j] = pa_y[j] * tpz_z_xyyy_0[j] + 1.5 * fl1_fx * tpz_z_xyy_0[j];

            tpx_yz_xyyz_0[j] = pa_y[j] * tpx_z_xyyz_0[j] + fl1_fx * tpx_z_xyz_0[j];

            tpy_yz_xyyz_0[j] = pa_y[j] * tpy_z_xyyz_0[j] + fl1_fx * tpy_z_xyz_0[j] - fl1_fgb * fl1_fx * ts_z_xyyz_0[j];

            tpz_yz_xyyz_0[j] = pa_y[j] * tpz_z_xyyz_0[j] + fl1_fx * tpz_z_xyz_0[j];

            tpx_yz_xyzz_0[j] = pa_y[j] * tpx_z_xyzz_0[j] + 0.5 * fl1_fx * tpx_z_xzz_0[j];

            tpy_yz_xyzz_0[j] = pa_y[j] * tpy_z_xyzz_0[j] + 0.5 * fl1_fx * tpy_z_xzz_0[j] - fl1_fgb * fl1_fx * ts_z_xyzz_0[j];

            tpz_yz_xyzz_0[j] = pa_y[j] * tpz_z_xyzz_0[j] + 0.5 * fl1_fx * tpz_z_xzz_0[j];

            tpx_yz_xzzz_0[j] = pa_y[j] * tpx_z_xzzz_0[j];

            tpy_yz_xzzz_0[j] = pa_y[j] * tpy_z_xzzz_0[j] - fl1_fgb * fl1_fx * ts_z_xzzz_0[j];

            tpz_yz_xzzz_0[j] = pa_y[j] * tpz_z_xzzz_0[j];

            tpx_yz_yyyy_0[j] = pa_y[j] * tpx_z_yyyy_0[j] + 2.0 * fl1_fx * tpx_z_yyy_0[j];

            tpy_yz_yyyy_0[j] = pa_y[j] * tpy_z_yyyy_0[j] + 2.0 * fl1_fx * tpy_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyyy_0[j];

            tpz_yz_yyyy_0[j] = pa_y[j] * tpz_z_yyyy_0[j] + 2.0 * fl1_fx * tpz_z_yyy_0[j];

            tpx_yz_yyyz_0[j] = pa_y[j] * tpx_z_yyyz_0[j] + 1.5 * fl1_fx * tpx_z_yyz_0[j];

            tpy_yz_yyyz_0[j] = pa_y[j] * tpy_z_yyyz_0[j] + 1.5 * fl1_fx * tpy_z_yyz_0[j] - fl1_fgb * fl1_fx * ts_z_yyyz_0[j];

            tpz_yz_yyyz_0[j] = pa_y[j] * tpz_z_yyyz_0[j] + 1.5 * fl1_fx * tpz_z_yyz_0[j];

            tpx_yz_yyzz_0[j] = pa_y[j] * tpx_z_yyzz_0[j] + fl1_fx * tpx_z_yzz_0[j];

            tpy_yz_yyzz_0[j] = pa_y[j] * tpy_z_yyzz_0[j] + fl1_fx * tpy_z_yzz_0[j] - fl1_fgb * fl1_fx * ts_z_yyzz_0[j];

            tpz_yz_yyzz_0[j] = pa_y[j] * tpz_z_yyzz_0[j] + fl1_fx * tpz_z_yzz_0[j];

            tpx_yz_yzzz_0[j] = pa_y[j] * tpx_z_yzzz_0[j] + 0.5 * fl1_fx * tpx_z_zzz_0[j];

            tpy_yz_yzzz_0[j] = pa_y[j] * tpy_z_yzzz_0[j] + 0.5 * fl1_fx * tpy_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_z_yzzz_0[j];

            tpz_yz_yzzz_0[j] = pa_y[j] * tpz_z_yzzz_0[j] + 0.5 * fl1_fx * tpz_z_zzz_0[j];

            tpx_yz_zzzz_0[j] = pa_y[j] * tpx_z_zzzz_0[j];

            tpy_yz_zzzz_0[j] = pa_y[j] * tpy_z_zzzz_0[j] - fl1_fgb * fl1_fx * ts_z_zzzz_0[j];

            tpz_yz_zzzz_0[j] = pa_y[j] * tpz_z_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDG_225_270(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (225,270)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

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

        // set up pointers to integrals

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpx_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 88);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tpx_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 89);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        // Batch of Integrals (225,270)

        #pragma omp simd aligned(fgb, fx, pa_z, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxyy_0, \
                                     tpx_0_xxyz_0, tpx_0_xxzz_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyzz_0, tpx_0_yzzz_0, tpx_0_zzzz_0, tpx_z_xxx_0, \
                                     tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, tpx_z_xxy_0, tpx_z_xxyy_0, tpx_z_xxyz_0, \
                                     tpx_z_xxz_0, tpx_z_xxzz_0, tpx_z_xyy_0, tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyz_0, \
                                     tpx_z_xyzz_0, tpx_z_xzz_0, tpx_z_xzzz_0, tpx_z_yyy_0, tpx_z_yyyy_0, tpx_z_yyyz_0, \
                                     tpx_z_yyz_0, tpx_z_yyzz_0, tpx_z_yzz_0, tpx_z_yzzz_0, tpx_z_zzz_0, tpx_z_zzzz_0, \
                                     tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, \
                                     tpx_zz_xxzz_0, tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyzz_0, tpx_zz_xzzz_0, \
                                     tpx_zz_yyyy_0, tpx_zz_yyyz_0, tpx_zz_yyzz_0, tpx_zz_yzzz_0, tpx_zz_zzzz_0, \
                                     tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxyy_0, tpy_0_xxyz_0, tpy_0_xxzz_0, \
                                     tpy_0_xyyy_0, tpy_0_xyyz_0, tpy_0_xyzz_0, tpy_0_xzzz_0, tpy_0_yyyy_0, tpy_0_yyyz_0, \
                                     tpy_0_yyzz_0, tpy_0_yzzz_0, tpy_0_zzzz_0, tpy_z_xxx_0, tpy_z_xxxx_0, tpy_z_xxxy_0, \
                                     tpy_z_xxxz_0, tpy_z_xxy_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxz_0, tpy_z_xxzz_0, \
                                     tpy_z_xyy_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyz_0, tpy_z_xyzz_0, tpy_z_xzz_0, \
                                     tpy_z_xzzz_0, tpy_z_yyy_0, tpy_z_yyyy_0, tpy_z_yyyz_0, tpy_z_yyz_0, tpy_z_yyzz_0, \
                                     tpy_z_yzz_0, tpy_z_yzzz_0, tpy_z_zzz_0, tpy_z_zzzz_0, tpy_zz_xxxx_0, \
                                     tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, \
                                     tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyzz_0, tpy_zz_xzzz_0, tpy_zz_yyyy_0, \
                                     tpy_zz_yyyz_0, tpy_zz_yyzz_0, tpy_zz_yzzz_0, tpy_zz_zzzz_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxzz_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyzz_0, tpz_0_xzzz_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yzzz_0, tpz_0_zzzz_0, tpz_z_xxx_0, tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, \
                                     tpz_z_xxy_0, tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxz_0, tpz_z_xxzz_0, tpz_z_xyy_0, \
                                     tpz_z_xyyy_0, tpz_z_xyyz_0, tpz_z_xyz_0, tpz_z_xyzz_0, tpz_z_xzz_0, tpz_z_xzzz_0, \
                                     tpz_z_yyy_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyz_0, tpz_z_yyzz_0, tpz_z_yzz_0, \
                                     tpz_z_yzzz_0, tpz_z_zzz_0, tpz_z_zzzz_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, \
                                     tpz_zz_xxxz_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, tpz_zz_xxzz_0, tpz_zz_xyyy_0, \
                                     tpz_zz_xyyz_0, tpz_zz_xyzz_0, tpz_zz_xzzz_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, \
                                     tpz_zz_yyzz_0, tpz_zz_yzzz_0, tpz_zz_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, \
                                     ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, \
                                     ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_zz_xxxx_0[j] = pa_z[j] * tpx_z_xxxx_0[j] + 0.5 * fl1_fx * tpx_0_xxxx_0[j];

            tpy_zz_xxxx_0[j] = pa_z[j] * tpy_z_xxxx_0[j] + 0.5 * fl1_fx * tpy_0_xxxx_0[j];

            tpz_zz_xxxx_0[j] = pa_z[j] * tpz_z_xxxx_0[j] + 0.5 * fl1_fx * tpz_0_xxxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxxx_0[j];

            tpx_zz_xxxy_0[j] = pa_z[j] * tpx_z_xxxy_0[j] + 0.5 * fl1_fx * tpx_0_xxxy_0[j];

            tpy_zz_xxxy_0[j] = pa_z[j] * tpy_z_xxxy_0[j] + 0.5 * fl1_fx * tpy_0_xxxy_0[j];

            tpz_zz_xxxy_0[j] = pa_z[j] * tpz_z_xxxy_0[j] + 0.5 * fl1_fx * tpz_0_xxxy_0[j] - fl1_fgb * fl1_fx * ts_z_xxxy_0[j];

            tpx_zz_xxxz_0[j] = pa_z[j] * tpx_z_xxxz_0[j] + 0.5 * fl1_fx * tpx_0_xxxz_0[j] + 0.5 * fl1_fx * tpx_z_xxx_0[j];

            tpy_zz_xxxz_0[j] = pa_z[j] * tpy_z_xxxz_0[j] + 0.5 * fl1_fx * tpy_0_xxxz_0[j] + 0.5 * fl1_fx * tpy_z_xxx_0[j];

            tpz_zz_xxxz_0[j] =
                pa_z[j] * tpz_z_xxxz_0[j] + 0.5 * fl1_fx * tpz_0_xxxz_0[j] + 0.5 * fl1_fx * tpz_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_z_xxxz_0[j];

            tpx_zz_xxyy_0[j] = pa_z[j] * tpx_z_xxyy_0[j] + 0.5 * fl1_fx * tpx_0_xxyy_0[j];

            tpy_zz_xxyy_0[j] = pa_z[j] * tpy_z_xxyy_0[j] + 0.5 * fl1_fx * tpy_0_xxyy_0[j];

            tpz_zz_xxyy_0[j] = pa_z[j] * tpz_z_xxyy_0[j] + 0.5 * fl1_fx * tpz_0_xxyy_0[j] - fl1_fgb * fl1_fx * ts_z_xxyy_0[j];

            tpx_zz_xxyz_0[j] = pa_z[j] * tpx_z_xxyz_0[j] + 0.5 * fl1_fx * tpx_0_xxyz_0[j] + 0.5 * fl1_fx * tpx_z_xxy_0[j];

            tpy_zz_xxyz_0[j] = pa_z[j] * tpy_z_xxyz_0[j] + 0.5 * fl1_fx * tpy_0_xxyz_0[j] + 0.5 * fl1_fx * tpy_z_xxy_0[j];

            tpz_zz_xxyz_0[j] =
                pa_z[j] * tpz_z_xxyz_0[j] + 0.5 * fl1_fx * tpz_0_xxyz_0[j] + 0.5 * fl1_fx * tpz_z_xxy_0[j] - fl1_fgb * fl1_fx * ts_z_xxyz_0[j];

            tpx_zz_xxzz_0[j] = pa_z[j] * tpx_z_xxzz_0[j] + 0.5 * fl1_fx * tpx_0_xxzz_0[j] + fl1_fx * tpx_z_xxz_0[j];

            tpy_zz_xxzz_0[j] = pa_z[j] * tpy_z_xxzz_0[j] + 0.5 * fl1_fx * tpy_0_xxzz_0[j] + fl1_fx * tpy_z_xxz_0[j];

            tpz_zz_xxzz_0[j] =
                pa_z[j] * tpz_z_xxzz_0[j] + 0.5 * fl1_fx * tpz_0_xxzz_0[j] + fl1_fx * tpz_z_xxz_0[j] - fl1_fgb * fl1_fx * ts_z_xxzz_0[j];

            tpx_zz_xyyy_0[j] = pa_z[j] * tpx_z_xyyy_0[j] + 0.5 * fl1_fx * tpx_0_xyyy_0[j];

            tpy_zz_xyyy_0[j] = pa_z[j] * tpy_z_xyyy_0[j] + 0.5 * fl1_fx * tpy_0_xyyy_0[j];

            tpz_zz_xyyy_0[j] = pa_z[j] * tpz_z_xyyy_0[j] + 0.5 * fl1_fx * tpz_0_xyyy_0[j] - fl1_fgb * fl1_fx * ts_z_xyyy_0[j];

            tpx_zz_xyyz_0[j] = pa_z[j] * tpx_z_xyyz_0[j] + 0.5 * fl1_fx * tpx_0_xyyz_0[j] + 0.5 * fl1_fx * tpx_z_xyy_0[j];

            tpy_zz_xyyz_0[j] = pa_z[j] * tpy_z_xyyz_0[j] + 0.5 * fl1_fx * tpy_0_xyyz_0[j] + 0.5 * fl1_fx * tpy_z_xyy_0[j];

            tpz_zz_xyyz_0[j] =
                pa_z[j] * tpz_z_xyyz_0[j] + 0.5 * fl1_fx * tpz_0_xyyz_0[j] + 0.5 * fl1_fx * tpz_z_xyy_0[j] - fl1_fgb * fl1_fx * ts_z_xyyz_0[j];

            tpx_zz_xyzz_0[j] = pa_z[j] * tpx_z_xyzz_0[j] + 0.5 * fl1_fx * tpx_0_xyzz_0[j] + fl1_fx * tpx_z_xyz_0[j];

            tpy_zz_xyzz_0[j] = pa_z[j] * tpy_z_xyzz_0[j] + 0.5 * fl1_fx * tpy_0_xyzz_0[j] + fl1_fx * tpy_z_xyz_0[j];

            tpz_zz_xyzz_0[j] =
                pa_z[j] * tpz_z_xyzz_0[j] + 0.5 * fl1_fx * tpz_0_xyzz_0[j] + fl1_fx * tpz_z_xyz_0[j] - fl1_fgb * fl1_fx * ts_z_xyzz_0[j];

            tpx_zz_xzzz_0[j] = pa_z[j] * tpx_z_xzzz_0[j] + 0.5 * fl1_fx * tpx_0_xzzz_0[j] + 1.5 * fl1_fx * tpx_z_xzz_0[j];

            tpy_zz_xzzz_0[j] = pa_z[j] * tpy_z_xzzz_0[j] + 0.5 * fl1_fx * tpy_0_xzzz_0[j] + 1.5 * fl1_fx * tpy_z_xzz_0[j];

            tpz_zz_xzzz_0[j] =
                pa_z[j] * tpz_z_xzzz_0[j] + 0.5 * fl1_fx * tpz_0_xzzz_0[j] + 1.5 * fl1_fx * tpz_z_xzz_0[j] - fl1_fgb * fl1_fx * ts_z_xzzz_0[j];

            tpx_zz_yyyy_0[j] = pa_z[j] * tpx_z_yyyy_0[j] + 0.5 * fl1_fx * tpx_0_yyyy_0[j];

            tpy_zz_yyyy_0[j] = pa_z[j] * tpy_z_yyyy_0[j] + 0.5 * fl1_fx * tpy_0_yyyy_0[j];

            tpz_zz_yyyy_0[j] = pa_z[j] * tpz_z_yyyy_0[j] + 0.5 * fl1_fx * tpz_0_yyyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyyy_0[j];

            tpx_zz_yyyz_0[j] = pa_z[j] * tpx_z_yyyz_0[j] + 0.5 * fl1_fx * tpx_0_yyyz_0[j] + 0.5 * fl1_fx * tpx_z_yyy_0[j];

            tpy_zz_yyyz_0[j] = pa_z[j] * tpy_z_yyyz_0[j] + 0.5 * fl1_fx * tpy_0_yyyz_0[j] + 0.5 * fl1_fx * tpy_z_yyy_0[j];

            tpz_zz_yyyz_0[j] =
                pa_z[j] * tpz_z_yyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyyz_0[j] + 0.5 * fl1_fx * tpz_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_z_yyyz_0[j];

            tpx_zz_yyzz_0[j] = pa_z[j] * tpx_z_yyzz_0[j] + 0.5 * fl1_fx * tpx_0_yyzz_0[j] + fl1_fx * tpx_z_yyz_0[j];

            tpy_zz_yyzz_0[j] = pa_z[j] * tpy_z_yyzz_0[j] + 0.5 * fl1_fx * tpy_0_yyzz_0[j] + fl1_fx * tpy_z_yyz_0[j];

            tpz_zz_yyzz_0[j] =
                pa_z[j] * tpz_z_yyzz_0[j] + 0.5 * fl1_fx * tpz_0_yyzz_0[j] + fl1_fx * tpz_z_yyz_0[j] - fl1_fgb * fl1_fx * ts_z_yyzz_0[j];

            tpx_zz_yzzz_0[j] = pa_z[j] * tpx_z_yzzz_0[j] + 0.5 * fl1_fx * tpx_0_yzzz_0[j] + 1.5 * fl1_fx * tpx_z_yzz_0[j];

            tpy_zz_yzzz_0[j] = pa_z[j] * tpy_z_yzzz_0[j] + 0.5 * fl1_fx * tpy_0_yzzz_0[j] + 1.5 * fl1_fx * tpy_z_yzz_0[j];

            tpz_zz_yzzz_0[j] =
                pa_z[j] * tpz_z_yzzz_0[j] + 0.5 * fl1_fx * tpz_0_yzzz_0[j] + 1.5 * fl1_fx * tpz_z_yzz_0[j] - fl1_fgb * fl1_fx * ts_z_yzzz_0[j];

            tpx_zz_zzzz_0[j] = pa_z[j] * tpx_z_zzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzzz_0[j] + 2.0 * fl1_fx * tpx_z_zzz_0[j];

            tpy_zz_zzzz_0[j] = pa_z[j] * tpy_z_zzzz_0[j] + 0.5 * fl1_fx * tpy_0_zzzz_0[j] + 2.0 * fl1_fx * tpy_z_zzz_0[j];

            tpz_zz_zzzz_0[j] =
                pa_z[j] * tpz_z_zzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzzz_0[j] + 2.0 * fl1_fx * tpz_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_z_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForGD_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGD_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGD_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGD_135_180(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGD_180_225(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGD_225_270(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForGD_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx);

        auto tpy_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx);

        auto tpz_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx);

        auto tpx_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 1);

        auto tpy_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 1);

        auto tpx_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 2);

        auto tpy_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 2);

        auto tpx_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 3);

        auto tpy_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 3);

        auto tpx_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 4);

        auto tpy_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 4);

        auto tpx_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 5);

        auto tpy_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 5);

        auto tpx_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 6);

        auto tpy_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 6);

        auto tpx_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 7);

        auto tpy_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 7);

        auto tpx_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 8);

        auto tpy_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 8);

        auto tpx_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 9);

        auto tpy_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 9);

        auto tpx_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 10);

        auto tpy_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 10);

        auto tpx_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 11);

        auto tpy_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 11);

        auto tpz_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 11);

        auto tpx_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 12);

        auto tpy_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 12);

        auto tpx_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 13);

        auto tpy_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 13);

        auto tpx_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 14);

        auto tpy_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 14);

        auto tpx_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx);

        auto tpy_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx);

        auto tpz_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx);

        auto tpx_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 1);

        auto tpy_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tpz_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tpx_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 2);

        auto tpy_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tpz_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tpx_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 3);

        auto tpy_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tpz_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tpx_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 4);

        auto tpy_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tpz_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tpx_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 5);

        auto tpy_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tpz_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tpx_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 6);

        auto tpy_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tpz_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tpx_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 7);

        auto tpy_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tpz_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tpx_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 8);

        auto tpy_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tpz_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tpx_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 9);

        auto tpy_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tpz_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tpx_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 10);

        auto tpy_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tpz_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tpx_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 11);

        auto tpy_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tpz_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 11);

        auto tpx_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 12);

        auto tpy_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tpz_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tpx_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 13);

        auto tpy_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tpz_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tpx_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 14);

        auto tpy_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tpz_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tpx_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx);

        auto tpy_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx);

        auto tpz_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx);

        auto tpx_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 1);

        auto tpy_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 2);

        auto tpy_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 3);

        auto tpy_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 4);

        auto tpy_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 5);

        auto tpy_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 6);

        auto tpy_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 7);

        auto tpy_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 8);

        auto tpy_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto ts_xxx_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx);

        auto ts_xxx_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 1);

        auto ts_xxx_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 2);

        auto ts_xxx_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 3);

        auto ts_xxx_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 4);

        auto ts_xxx_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 5);

        auto ts_xxy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 6);

        auto ts_xxy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 7);

        auto ts_xxy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 8);

        auto ts_xxy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 9);

        auto ts_xxy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 10);

        auto ts_xxy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 11);

        auto ts_xxz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 12);

        auto ts_xxz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 13);

        auto ts_xxz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 14);

        // set up pointers to integrals

        auto tpx_xxxx_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx);

        auto tpy_xxxx_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx);

        auto tpz_xxxx_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx);

        auto tpx_xxxx_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 1);

        auto tpy_xxxx_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 1);

        auto tpz_xxxx_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 1);

        auto tpx_xxxx_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 2);

        auto tpy_xxxx_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 2);

        auto tpz_xxxx_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 2);

        auto tpx_xxxx_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 3);

        auto tpy_xxxx_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 3);

        auto tpz_xxxx_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 3);

        auto tpx_xxxx_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 4);

        auto tpy_xxxx_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 4);

        auto tpz_xxxx_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 4);

        auto tpx_xxxx_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 5);

        auto tpy_xxxx_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 5);

        auto tpz_xxxx_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 5);

        auto tpx_xxxy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 6);

        auto tpy_xxxy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 6);

        auto tpz_xxxy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 6);

        auto tpx_xxxy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 7);

        auto tpy_xxxy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 7);

        auto tpz_xxxy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 7);

        auto tpx_xxxy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 8);

        auto tpy_xxxy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 8);

        auto tpz_xxxy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 8);

        auto tpx_xxxy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 9);

        auto tpy_xxxy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 9);

        auto tpz_xxxy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 9);

        auto tpx_xxxy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 10);

        auto tpy_xxxy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 10);

        auto tpz_xxxy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 10);

        auto tpx_xxxy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 11);

        auto tpy_xxxy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 11);

        auto tpz_xxxy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 11);

        auto tpx_xxxz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 12);

        auto tpy_xxxz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 12);

        auto tpz_xxxz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 12);

        auto tpx_xxxz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 13);

        auto tpy_xxxz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 13);

        auto tpz_xxxz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 13);

        auto tpx_xxxz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 14);

        auto tpy_xxxz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 14);

        auto tpz_xxxz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xx_xx_0, tpx_xx_xy_0, tpx_xx_xz_0, tpx_xx_yy_0, \
                                     tpx_xx_yz_0, tpx_xx_zz_0, tpx_xxx_x_0, tpx_xxx_xx_0, tpx_xxx_xy_0, tpx_xxx_xz_0, \
                                     tpx_xxx_y_0, tpx_xxx_yy_0, tpx_xxx_yz_0, tpx_xxx_z_0, tpx_xxx_zz_0, tpx_xxxx_xx_0, \
                                     tpx_xxxx_xy_0, tpx_xxxx_xz_0, tpx_xxxx_yy_0, tpx_xxxx_yz_0, tpx_xxxx_zz_0, \
                                     tpx_xxxy_xx_0, tpx_xxxy_xy_0, tpx_xxxy_xz_0, tpx_xxxy_yy_0, tpx_xxxy_yz_0, \
                                     tpx_xxxy_zz_0, tpx_xxxz_xx_0, tpx_xxxz_xy_0, tpx_xxxz_xz_0, tpx_xxy_x_0, \
                                     tpx_xxy_xx_0, tpx_xxy_xy_0, tpx_xxy_xz_0, tpx_xxy_y_0, tpx_xxy_yy_0, tpx_xxy_yz_0, \
                                     tpx_xxy_z_0, tpx_xxy_zz_0, tpx_xxz_x_0, tpx_xxz_xx_0, tpx_xxz_xy_0, tpx_xxz_xz_0, \
                                     tpx_xxz_y_0, tpx_xxz_z_0, tpx_xy_xx_0, tpx_xy_xy_0, tpx_xy_xz_0, tpx_xy_yy_0, \
                                     tpx_xy_yz_0, tpx_xy_zz_0, tpx_xz_xx_0, tpx_xz_xy_0, tpx_xz_xz_0, tpy_xx_xx_0, \
                                     tpy_xx_xy_0, tpy_xx_xz_0, tpy_xx_yy_0, tpy_xx_yz_0, tpy_xx_zz_0, tpy_xxx_x_0, \
                                     tpy_xxx_xx_0, tpy_xxx_xy_0, tpy_xxx_xz_0, tpy_xxx_y_0, tpy_xxx_yy_0, tpy_xxx_yz_0, \
                                     tpy_xxx_z_0, tpy_xxx_zz_0, tpy_xxxx_xx_0, tpy_xxxx_xy_0, tpy_xxxx_xz_0, \
                                     tpy_xxxx_yy_0, tpy_xxxx_yz_0, tpy_xxxx_zz_0, tpy_xxxy_xx_0, tpy_xxxy_xy_0, \
                                     tpy_xxxy_xz_0, tpy_xxxy_yy_0, tpy_xxxy_yz_0, tpy_xxxy_zz_0, tpy_xxxz_xx_0, \
                                     tpy_xxxz_xy_0, tpy_xxxz_xz_0, tpy_xxy_x_0, tpy_xxy_xx_0, tpy_xxy_xy_0, tpy_xxy_xz_0, \
                                     tpy_xxy_y_0, tpy_xxy_yy_0, tpy_xxy_yz_0, tpy_xxy_z_0, tpy_xxy_zz_0, tpy_xxz_x_0, \
                                     tpy_xxz_xx_0, tpy_xxz_xy_0, tpy_xxz_xz_0, tpy_xxz_y_0, tpy_xxz_z_0, tpy_xy_xx_0, \
                                     tpy_xy_xy_0, tpy_xy_xz_0, tpy_xy_yy_0, tpy_xy_yz_0, tpy_xy_zz_0, tpy_xz_xx_0, \
                                     tpy_xz_xy_0, tpy_xz_xz_0, tpz_xx_xx_0, tpz_xx_xy_0, tpz_xx_xz_0, tpz_xx_yy_0, \
                                     tpz_xx_yz_0, tpz_xx_zz_0, tpz_xxx_x_0, tpz_xxx_xx_0, tpz_xxx_xy_0, tpz_xxx_xz_0, \
                                     tpz_xxx_y_0, tpz_xxx_yy_0, tpz_xxx_yz_0, tpz_xxx_z_0, tpz_xxx_zz_0, tpz_xxxx_xx_0, \
                                     tpz_xxxx_xy_0, tpz_xxxx_xz_0, tpz_xxxx_yy_0, tpz_xxxx_yz_0, tpz_xxxx_zz_0, \
                                     tpz_xxxy_xx_0, tpz_xxxy_xy_0, tpz_xxxy_xz_0, tpz_xxxy_yy_0, tpz_xxxy_yz_0, \
                                     tpz_xxxy_zz_0, tpz_xxxz_xx_0, tpz_xxxz_xy_0, tpz_xxxz_xz_0, tpz_xxy_x_0, \
                                     tpz_xxy_xx_0, tpz_xxy_xy_0, tpz_xxy_xz_0, tpz_xxy_y_0, tpz_xxy_yy_0, tpz_xxy_yz_0, \
                                     tpz_xxy_z_0, tpz_xxy_zz_0, tpz_xxz_x_0, tpz_xxz_xx_0, tpz_xxz_xy_0, tpz_xxz_xz_0, \
                                     tpz_xxz_y_0, tpz_xxz_z_0, tpz_xy_xx_0, tpz_xy_xy_0, tpz_xy_xz_0, tpz_xy_yy_0, \
                                     tpz_xy_yz_0, tpz_xy_zz_0, tpz_xz_xx_0, tpz_xz_xy_0, tpz_xz_xz_0, ts_xxx_xx_0, \
                                     ts_xxx_xy_0, ts_xxx_xz_0, ts_xxx_yy_0, ts_xxx_yz_0, ts_xxx_zz_0, ts_xxy_xx_0, \
                                     ts_xxy_xy_0, ts_xxy_xz_0, ts_xxy_yy_0, ts_xxy_yz_0, ts_xxy_zz_0, ts_xxz_xx_0, \
                                     ts_xxz_xy_0, ts_xxz_xz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxx_xx_0[j] =
                pa_x[j] * tpx_xxx_xx_0[j] + 1.5 * fl1_fx * tpx_xx_xx_0[j] + fl1_fx * tpx_xxx_x_0[j] - fl1_fgb * fl1_fx * ts_xxx_xx_0[j];

            tpy_xxxx_xx_0[j] = pa_x[j] * tpy_xxx_xx_0[j] + 1.5 * fl1_fx * tpy_xx_xx_0[j] + fl1_fx * tpy_xxx_x_0[j];

            tpz_xxxx_xx_0[j] = pa_x[j] * tpz_xxx_xx_0[j] + 1.5 * fl1_fx * tpz_xx_xx_0[j] + fl1_fx * tpz_xxx_x_0[j];

            tpx_xxxx_xy_0[j] =
                pa_x[j] * tpx_xxx_xy_0[j] + 1.5 * fl1_fx * tpx_xx_xy_0[j] + 0.5 * fl1_fx * tpx_xxx_y_0[j] - fl1_fgb * fl1_fx * ts_xxx_xy_0[j];

            tpy_xxxx_xy_0[j] = pa_x[j] * tpy_xxx_xy_0[j] + 1.5 * fl1_fx * tpy_xx_xy_0[j] + 0.5 * fl1_fx * tpy_xxx_y_0[j];

            tpz_xxxx_xy_0[j] = pa_x[j] * tpz_xxx_xy_0[j] + 1.5 * fl1_fx * tpz_xx_xy_0[j] + 0.5 * fl1_fx * tpz_xxx_y_0[j];

            tpx_xxxx_xz_0[j] =
                pa_x[j] * tpx_xxx_xz_0[j] + 1.5 * fl1_fx * tpx_xx_xz_0[j] + 0.5 * fl1_fx * tpx_xxx_z_0[j] - fl1_fgb * fl1_fx * ts_xxx_xz_0[j];

            tpy_xxxx_xz_0[j] = pa_x[j] * tpy_xxx_xz_0[j] + 1.5 * fl1_fx * tpy_xx_xz_0[j] + 0.5 * fl1_fx * tpy_xxx_z_0[j];

            tpz_xxxx_xz_0[j] = pa_x[j] * tpz_xxx_xz_0[j] + 1.5 * fl1_fx * tpz_xx_xz_0[j] + 0.5 * fl1_fx * tpz_xxx_z_0[j];

            tpx_xxxx_yy_0[j] = pa_x[j] * tpx_xxx_yy_0[j] + 1.5 * fl1_fx * tpx_xx_yy_0[j] - fl1_fgb * fl1_fx * ts_xxx_yy_0[j];

            tpy_xxxx_yy_0[j] = pa_x[j] * tpy_xxx_yy_0[j] + 1.5 * fl1_fx * tpy_xx_yy_0[j];

            tpz_xxxx_yy_0[j] = pa_x[j] * tpz_xxx_yy_0[j] + 1.5 * fl1_fx * tpz_xx_yy_0[j];

            tpx_xxxx_yz_0[j] = pa_x[j] * tpx_xxx_yz_0[j] + 1.5 * fl1_fx * tpx_xx_yz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yz_0[j];

            tpy_xxxx_yz_0[j] = pa_x[j] * tpy_xxx_yz_0[j] + 1.5 * fl1_fx * tpy_xx_yz_0[j];

            tpz_xxxx_yz_0[j] = pa_x[j] * tpz_xxx_yz_0[j] + 1.5 * fl1_fx * tpz_xx_yz_0[j];

            tpx_xxxx_zz_0[j] = pa_x[j] * tpx_xxx_zz_0[j] + 1.5 * fl1_fx * tpx_xx_zz_0[j] - fl1_fgb * fl1_fx * ts_xxx_zz_0[j];

            tpy_xxxx_zz_0[j] = pa_x[j] * tpy_xxx_zz_0[j] + 1.5 * fl1_fx * tpy_xx_zz_0[j];

            tpz_xxxx_zz_0[j] = pa_x[j] * tpz_xxx_zz_0[j] + 1.5 * fl1_fx * tpz_xx_zz_0[j];

            tpx_xxxy_xx_0[j] = pa_x[j] * tpx_xxy_xx_0[j] + fl1_fx * tpx_xy_xx_0[j] + fl1_fx * tpx_xxy_x_0[j] - fl1_fgb * fl1_fx * ts_xxy_xx_0[j];

            tpy_xxxy_xx_0[j] = pa_x[j] * tpy_xxy_xx_0[j] + fl1_fx * tpy_xy_xx_0[j] + fl1_fx * tpy_xxy_x_0[j];

            tpz_xxxy_xx_0[j] = pa_x[j] * tpz_xxy_xx_0[j] + fl1_fx * tpz_xy_xx_0[j] + fl1_fx * tpz_xxy_x_0[j];

            tpx_xxxy_xy_0[j] =
                pa_x[j] * tpx_xxy_xy_0[j] + fl1_fx * tpx_xy_xy_0[j] + 0.5 * fl1_fx * tpx_xxy_y_0[j] - fl1_fgb * fl1_fx * ts_xxy_xy_0[j];

            tpy_xxxy_xy_0[j] = pa_x[j] * tpy_xxy_xy_0[j] + fl1_fx * tpy_xy_xy_0[j] + 0.5 * fl1_fx * tpy_xxy_y_0[j];

            tpz_xxxy_xy_0[j] = pa_x[j] * tpz_xxy_xy_0[j] + fl1_fx * tpz_xy_xy_0[j] + 0.5 * fl1_fx * tpz_xxy_y_0[j];

            tpx_xxxy_xz_0[j] =
                pa_x[j] * tpx_xxy_xz_0[j] + fl1_fx * tpx_xy_xz_0[j] + 0.5 * fl1_fx * tpx_xxy_z_0[j] - fl1_fgb * fl1_fx * ts_xxy_xz_0[j];

            tpy_xxxy_xz_0[j] = pa_x[j] * tpy_xxy_xz_0[j] + fl1_fx * tpy_xy_xz_0[j] + 0.5 * fl1_fx * tpy_xxy_z_0[j];

            tpz_xxxy_xz_0[j] = pa_x[j] * tpz_xxy_xz_0[j] + fl1_fx * tpz_xy_xz_0[j] + 0.5 * fl1_fx * tpz_xxy_z_0[j];

            tpx_xxxy_yy_0[j] = pa_x[j] * tpx_xxy_yy_0[j] + fl1_fx * tpx_xy_yy_0[j] - fl1_fgb * fl1_fx * ts_xxy_yy_0[j];

            tpy_xxxy_yy_0[j] = pa_x[j] * tpy_xxy_yy_0[j] + fl1_fx * tpy_xy_yy_0[j];

            tpz_xxxy_yy_0[j] = pa_x[j] * tpz_xxy_yy_0[j] + fl1_fx * tpz_xy_yy_0[j];

            tpx_xxxy_yz_0[j] = pa_x[j] * tpx_xxy_yz_0[j] + fl1_fx * tpx_xy_yz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yz_0[j];

            tpy_xxxy_yz_0[j] = pa_x[j] * tpy_xxy_yz_0[j] + fl1_fx * tpy_xy_yz_0[j];

            tpz_xxxy_yz_0[j] = pa_x[j] * tpz_xxy_yz_0[j] + fl1_fx * tpz_xy_yz_0[j];

            tpx_xxxy_zz_0[j] = pa_x[j] * tpx_xxy_zz_0[j] + fl1_fx * tpx_xy_zz_0[j] - fl1_fgb * fl1_fx * ts_xxy_zz_0[j];

            tpy_xxxy_zz_0[j] = pa_x[j] * tpy_xxy_zz_0[j] + fl1_fx * tpy_xy_zz_0[j];

            tpz_xxxy_zz_0[j] = pa_x[j] * tpz_xxy_zz_0[j] + fl1_fx * tpz_xy_zz_0[j];

            tpx_xxxz_xx_0[j] = pa_x[j] * tpx_xxz_xx_0[j] + fl1_fx * tpx_xz_xx_0[j] + fl1_fx * tpx_xxz_x_0[j] - fl1_fgb * fl1_fx * ts_xxz_xx_0[j];

            tpy_xxxz_xx_0[j] = pa_x[j] * tpy_xxz_xx_0[j] + fl1_fx * tpy_xz_xx_0[j] + fl1_fx * tpy_xxz_x_0[j];

            tpz_xxxz_xx_0[j] = pa_x[j] * tpz_xxz_xx_0[j] + fl1_fx * tpz_xz_xx_0[j] + fl1_fx * tpz_xxz_x_0[j];

            tpx_xxxz_xy_0[j] =
                pa_x[j] * tpx_xxz_xy_0[j] + fl1_fx * tpx_xz_xy_0[j] + 0.5 * fl1_fx * tpx_xxz_y_0[j] - fl1_fgb * fl1_fx * ts_xxz_xy_0[j];

            tpy_xxxz_xy_0[j] = pa_x[j] * tpy_xxz_xy_0[j] + fl1_fx * tpy_xz_xy_0[j] + 0.5 * fl1_fx * tpy_xxz_y_0[j];

            tpz_xxxz_xy_0[j] = pa_x[j] * tpz_xxz_xy_0[j] + fl1_fx * tpz_xz_xy_0[j] + 0.5 * fl1_fx * tpz_xxz_y_0[j];

            tpx_xxxz_xz_0[j] =
                pa_x[j] * tpx_xxz_xz_0[j] + fl1_fx * tpx_xz_xz_0[j] + 0.5 * fl1_fx * tpx_xxz_z_0[j] - fl1_fgb * fl1_fx * ts_xxz_xz_0[j];

            tpy_xxxz_xz_0[j] = pa_x[j] * tpy_xxz_xz_0[j] + fl1_fx * tpy_xz_xz_0[j] + 0.5 * fl1_fx * tpy_xxz_z_0[j];

            tpz_xxxz_xz_0[j] = pa_x[j] * tpz_xxz_xz_0[j] + fl1_fx * tpz_xz_xz_0[j] + 0.5 * fl1_fx * tpz_xxz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 15);

        auto tpy_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 15);

        auto tpx_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 16);

        auto tpy_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 16);

        auto tpz_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 16);

        auto tpx_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 17);

        auto tpy_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 17);

        auto tpx_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 18);

        auto tpy_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 18);

        auto tpx_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 19);

        auto tpy_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 19);

        auto tpx_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 20);

        auto tpy_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 20);

        auto tpx_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 21);

        auto tpy_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 21);

        auto tpx_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 22);

        auto tpy_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 22);

        auto tpx_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 23);

        auto tpy_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 23);

        auto tpx_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 24);

        auto tpy_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 24);

        auto tpx_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 25);

        auto tpy_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 25);

        auto tpx_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 26);

        auto tpy_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 26);

        auto tpx_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 27);

        auto tpy_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 27);

        auto tpx_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 28);

        auto tpy_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 28);

        auto tpx_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 29);

        auto tpy_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 29);

        auto tpx_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 15);

        auto tpy_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tpz_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tpx_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 16);

        auto tpy_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tpz_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tpx_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 17);

        auto tpy_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tpz_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 9);

        auto tpy_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 10);

        auto tpy_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 11);

        auto tpy_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 12);

        auto tpy_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 13);

        auto tpy_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 14);

        auto tpy_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 14);

        auto ts_xxz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 15);

        auto ts_xxz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 16);

        auto ts_xxz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 17);

        auto ts_xyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 18);

        auto ts_xyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 19);

        auto ts_xyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 20);

        auto ts_xyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 21);

        auto ts_xyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 22);

        auto ts_xyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 23);

        auto ts_xyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 24);

        auto ts_xyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 25);

        auto ts_xyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 26);

        auto ts_xyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 27);

        auto ts_xyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 28);

        auto ts_xyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 29);

        // set up pointers to integrals

        auto tpx_xxxz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 15);

        auto tpy_xxxz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 15);

        auto tpz_xxxz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 15);

        auto tpx_xxxz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 16);

        auto tpy_xxxz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 16);

        auto tpz_xxxz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 16);

        auto tpx_xxxz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 17);

        auto tpy_xxxz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 17);

        auto tpz_xxxz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 17);

        auto tpx_xxyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 18);

        auto tpy_xxyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 18);

        auto tpz_xxyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 18);

        auto tpx_xxyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 19);

        auto tpy_xxyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 19);

        auto tpz_xxyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 19);

        auto tpx_xxyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 20);

        auto tpy_xxyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 20);

        auto tpz_xxyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 20);

        auto tpx_xxyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 21);

        auto tpy_xxyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 21);

        auto tpz_xxyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 21);

        auto tpx_xxyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 22);

        auto tpy_xxyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 22);

        auto tpz_xxyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 22);

        auto tpx_xxyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 23);

        auto tpy_xxyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 23);

        auto tpz_xxyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 23);

        auto tpx_xxyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 24);

        auto tpy_xxyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 24);

        auto tpz_xxyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 24);

        auto tpx_xxyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 25);

        auto tpy_xxyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 25);

        auto tpz_xxyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 25);

        auto tpx_xxyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 26);

        auto tpy_xxyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 26);

        auto tpz_xxyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 26);

        auto tpx_xxyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 27);

        auto tpy_xxyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 27);

        auto tpz_xxyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 27);

        auto tpx_xxyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 28);

        auto tpy_xxyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 28);

        auto tpz_xxyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 28);

        auto tpx_xxyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 29);

        auto tpy_xxyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 29);

        auto tpz_xxyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxxz_yy_0, tpx_xxxz_yz_0, tpx_xxxz_zz_0, \
                                     tpx_xxyy_xx_0, tpx_xxyy_xy_0, tpx_xxyy_xz_0, tpx_xxyy_yy_0, tpx_xxyy_yz_0, \
                                     tpx_xxyy_zz_0, tpx_xxyz_xx_0, tpx_xxyz_xy_0, tpx_xxyz_xz_0, tpx_xxyz_yy_0, \
                                     tpx_xxyz_yz_0, tpx_xxyz_zz_0, tpx_xxz_yy_0, tpx_xxz_yz_0, tpx_xxz_zz_0, tpx_xyy_x_0, \
                                     tpx_xyy_xx_0, tpx_xyy_xy_0, tpx_xyy_xz_0, tpx_xyy_y_0, tpx_xyy_yy_0, tpx_xyy_yz_0, \
                                     tpx_xyy_z_0, tpx_xyy_zz_0, tpx_xyz_x_0, tpx_xyz_xx_0, tpx_xyz_xy_0, tpx_xyz_xz_0, \
                                     tpx_xyz_y_0, tpx_xyz_yy_0, tpx_xyz_yz_0, tpx_xyz_z_0, tpx_xyz_zz_0, tpx_xz_yy_0, \
                                     tpx_xz_yz_0, tpx_xz_zz_0, tpx_yy_xx_0, tpx_yy_xy_0, tpx_yy_xz_0, tpx_yy_yy_0, \
                                     tpx_yy_yz_0, tpx_yy_zz_0, tpx_yz_xx_0, tpx_yz_xy_0, tpx_yz_xz_0, tpx_yz_yy_0, \
                                     tpx_yz_yz_0, tpx_yz_zz_0, tpy_xxxz_yy_0, tpy_xxxz_yz_0, tpy_xxxz_zz_0, \
                                     tpy_xxyy_xx_0, tpy_xxyy_xy_0, tpy_xxyy_xz_0, tpy_xxyy_yy_0, tpy_xxyy_yz_0, \
                                     tpy_xxyy_zz_0, tpy_xxyz_xx_0, tpy_xxyz_xy_0, tpy_xxyz_xz_0, tpy_xxyz_yy_0, \
                                     tpy_xxyz_yz_0, tpy_xxyz_zz_0, tpy_xxz_yy_0, tpy_xxz_yz_0, tpy_xxz_zz_0, tpy_xyy_x_0, \
                                     tpy_xyy_xx_0, tpy_xyy_xy_0, tpy_xyy_xz_0, tpy_xyy_y_0, tpy_xyy_yy_0, tpy_xyy_yz_0, \
                                     tpy_xyy_z_0, tpy_xyy_zz_0, tpy_xyz_x_0, tpy_xyz_xx_0, tpy_xyz_xy_0, tpy_xyz_xz_0, \
                                     tpy_xyz_y_0, tpy_xyz_yy_0, tpy_xyz_yz_0, tpy_xyz_z_0, tpy_xyz_zz_0, tpy_xz_yy_0, \
                                     tpy_xz_yz_0, tpy_xz_zz_0, tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_yy_yy_0, \
                                     tpy_yy_yz_0, tpy_yy_zz_0, tpy_yz_xx_0, tpy_yz_xy_0, tpy_yz_xz_0, tpy_yz_yy_0, \
                                     tpy_yz_yz_0, tpy_yz_zz_0, tpz_xxxz_yy_0, tpz_xxxz_yz_0, tpz_xxxz_zz_0, \
                                     tpz_xxyy_xx_0, tpz_xxyy_xy_0, tpz_xxyy_xz_0, tpz_xxyy_yy_0, tpz_xxyy_yz_0, \
                                     tpz_xxyy_zz_0, tpz_xxyz_xx_0, tpz_xxyz_xy_0, tpz_xxyz_xz_0, tpz_xxyz_yy_0, \
                                     tpz_xxyz_yz_0, tpz_xxyz_zz_0, tpz_xxz_yy_0, tpz_xxz_yz_0, tpz_xxz_zz_0, tpz_xyy_x_0, \
                                     tpz_xyy_xx_0, tpz_xyy_xy_0, tpz_xyy_xz_0, tpz_xyy_y_0, tpz_xyy_yy_0, tpz_xyy_yz_0, \
                                     tpz_xyy_z_0, tpz_xyy_zz_0, tpz_xyz_x_0, tpz_xyz_xx_0, tpz_xyz_xy_0, tpz_xyz_xz_0, \
                                     tpz_xyz_y_0, tpz_xyz_yy_0, tpz_xyz_yz_0, tpz_xyz_z_0, tpz_xyz_zz_0, tpz_xz_yy_0, \
                                     tpz_xz_yz_0, tpz_xz_zz_0, tpz_yy_xx_0, tpz_yy_xy_0, tpz_yy_xz_0, tpz_yy_yy_0, \
                                     tpz_yy_yz_0, tpz_yy_zz_0, tpz_yz_xx_0, tpz_yz_xy_0, tpz_yz_xz_0, tpz_yz_yy_0, \
                                     tpz_yz_yz_0, tpz_yz_zz_0, ts_xxz_yy_0, ts_xxz_yz_0, ts_xxz_zz_0, ts_xyy_xx_0, \
                                     ts_xyy_xy_0, ts_xyy_xz_0, ts_xyy_yy_0, ts_xyy_yz_0, ts_xyy_zz_0, ts_xyz_xx_0, \
                                     ts_xyz_xy_0, ts_xyz_xz_0, ts_xyz_yy_0, ts_xyz_yz_0, ts_xyz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxz_yy_0[j] = pa_x[j] * tpx_xxz_yy_0[j] + fl1_fx * tpx_xz_yy_0[j] - fl1_fgb * fl1_fx * ts_xxz_yy_0[j];

            tpy_xxxz_yy_0[j] = pa_x[j] * tpy_xxz_yy_0[j] + fl1_fx * tpy_xz_yy_0[j];

            tpz_xxxz_yy_0[j] = pa_x[j] * tpz_xxz_yy_0[j] + fl1_fx * tpz_xz_yy_0[j];

            tpx_xxxz_yz_0[j] = pa_x[j] * tpx_xxz_yz_0[j] + fl1_fx * tpx_xz_yz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yz_0[j];

            tpy_xxxz_yz_0[j] = pa_x[j] * tpy_xxz_yz_0[j] + fl1_fx * tpy_xz_yz_0[j];

            tpz_xxxz_yz_0[j] = pa_x[j] * tpz_xxz_yz_0[j] + fl1_fx * tpz_xz_yz_0[j];

            tpx_xxxz_zz_0[j] = pa_x[j] * tpx_xxz_zz_0[j] + fl1_fx * tpx_xz_zz_0[j] - fl1_fgb * fl1_fx * ts_xxz_zz_0[j];

            tpy_xxxz_zz_0[j] = pa_x[j] * tpy_xxz_zz_0[j] + fl1_fx * tpy_xz_zz_0[j];

            tpz_xxxz_zz_0[j] = pa_x[j] * tpz_xxz_zz_0[j] + fl1_fx * tpz_xz_zz_0[j];

            tpx_xxyy_xx_0[j] =
                pa_x[j] * tpx_xyy_xx_0[j] + 0.5 * fl1_fx * tpx_yy_xx_0[j] + fl1_fx * tpx_xyy_x_0[j] - fl1_fgb * fl1_fx * ts_xyy_xx_0[j];

            tpy_xxyy_xx_0[j] = pa_x[j] * tpy_xyy_xx_0[j] + 0.5 * fl1_fx * tpy_yy_xx_0[j] + fl1_fx * tpy_xyy_x_0[j];

            tpz_xxyy_xx_0[j] = pa_x[j] * tpz_xyy_xx_0[j] + 0.5 * fl1_fx * tpz_yy_xx_0[j] + fl1_fx * tpz_xyy_x_0[j];

            tpx_xxyy_xy_0[j] =
                pa_x[j] * tpx_xyy_xy_0[j] + 0.5 * fl1_fx * tpx_yy_xy_0[j] + 0.5 * fl1_fx * tpx_xyy_y_0[j] - fl1_fgb * fl1_fx * ts_xyy_xy_0[j];

            tpy_xxyy_xy_0[j] = pa_x[j] * tpy_xyy_xy_0[j] + 0.5 * fl1_fx * tpy_yy_xy_0[j] + 0.5 * fl1_fx * tpy_xyy_y_0[j];

            tpz_xxyy_xy_0[j] = pa_x[j] * tpz_xyy_xy_0[j] + 0.5 * fl1_fx * tpz_yy_xy_0[j] + 0.5 * fl1_fx * tpz_xyy_y_0[j];

            tpx_xxyy_xz_0[j] =
                pa_x[j] * tpx_xyy_xz_0[j] + 0.5 * fl1_fx * tpx_yy_xz_0[j] + 0.5 * fl1_fx * tpx_xyy_z_0[j] - fl1_fgb * fl1_fx * ts_xyy_xz_0[j];

            tpy_xxyy_xz_0[j] = pa_x[j] * tpy_xyy_xz_0[j] + 0.5 * fl1_fx * tpy_yy_xz_0[j] + 0.5 * fl1_fx * tpy_xyy_z_0[j];

            tpz_xxyy_xz_0[j] = pa_x[j] * tpz_xyy_xz_0[j] + 0.5 * fl1_fx * tpz_yy_xz_0[j] + 0.5 * fl1_fx * tpz_xyy_z_0[j];

            tpx_xxyy_yy_0[j] = pa_x[j] * tpx_xyy_yy_0[j] + 0.5 * fl1_fx * tpx_yy_yy_0[j] - fl1_fgb * fl1_fx * ts_xyy_yy_0[j];

            tpy_xxyy_yy_0[j] = pa_x[j] * tpy_xyy_yy_0[j] + 0.5 * fl1_fx * tpy_yy_yy_0[j];

            tpz_xxyy_yy_0[j] = pa_x[j] * tpz_xyy_yy_0[j] + 0.5 * fl1_fx * tpz_yy_yy_0[j];

            tpx_xxyy_yz_0[j] = pa_x[j] * tpx_xyy_yz_0[j] + 0.5 * fl1_fx * tpx_yy_yz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yz_0[j];

            tpy_xxyy_yz_0[j] = pa_x[j] * tpy_xyy_yz_0[j] + 0.5 * fl1_fx * tpy_yy_yz_0[j];

            tpz_xxyy_yz_0[j] = pa_x[j] * tpz_xyy_yz_0[j] + 0.5 * fl1_fx * tpz_yy_yz_0[j];

            tpx_xxyy_zz_0[j] = pa_x[j] * tpx_xyy_zz_0[j] + 0.5 * fl1_fx * tpx_yy_zz_0[j] - fl1_fgb * fl1_fx * ts_xyy_zz_0[j];

            tpy_xxyy_zz_0[j] = pa_x[j] * tpy_xyy_zz_0[j] + 0.5 * fl1_fx * tpy_yy_zz_0[j];

            tpz_xxyy_zz_0[j] = pa_x[j] * tpz_xyy_zz_0[j] + 0.5 * fl1_fx * tpz_yy_zz_0[j];

            tpx_xxyz_xx_0[j] =
                pa_x[j] * tpx_xyz_xx_0[j] + 0.5 * fl1_fx * tpx_yz_xx_0[j] + fl1_fx * tpx_xyz_x_0[j] - fl1_fgb * fl1_fx * ts_xyz_xx_0[j];

            tpy_xxyz_xx_0[j] = pa_x[j] * tpy_xyz_xx_0[j] + 0.5 * fl1_fx * tpy_yz_xx_0[j] + fl1_fx * tpy_xyz_x_0[j];

            tpz_xxyz_xx_0[j] = pa_x[j] * tpz_xyz_xx_0[j] + 0.5 * fl1_fx * tpz_yz_xx_0[j] + fl1_fx * tpz_xyz_x_0[j];

            tpx_xxyz_xy_0[j] =
                pa_x[j] * tpx_xyz_xy_0[j] + 0.5 * fl1_fx * tpx_yz_xy_0[j] + 0.5 * fl1_fx * tpx_xyz_y_0[j] - fl1_fgb * fl1_fx * ts_xyz_xy_0[j];

            tpy_xxyz_xy_0[j] = pa_x[j] * tpy_xyz_xy_0[j] + 0.5 * fl1_fx * tpy_yz_xy_0[j] + 0.5 * fl1_fx * tpy_xyz_y_0[j];

            tpz_xxyz_xy_0[j] = pa_x[j] * tpz_xyz_xy_0[j] + 0.5 * fl1_fx * tpz_yz_xy_0[j] + 0.5 * fl1_fx * tpz_xyz_y_0[j];

            tpx_xxyz_xz_0[j] =
                pa_x[j] * tpx_xyz_xz_0[j] + 0.5 * fl1_fx * tpx_yz_xz_0[j] + 0.5 * fl1_fx * tpx_xyz_z_0[j] - fl1_fgb * fl1_fx * ts_xyz_xz_0[j];

            tpy_xxyz_xz_0[j] = pa_x[j] * tpy_xyz_xz_0[j] + 0.5 * fl1_fx * tpy_yz_xz_0[j] + 0.5 * fl1_fx * tpy_xyz_z_0[j];

            tpz_xxyz_xz_0[j] = pa_x[j] * tpz_xyz_xz_0[j] + 0.5 * fl1_fx * tpz_yz_xz_0[j] + 0.5 * fl1_fx * tpz_xyz_z_0[j];

            tpx_xxyz_yy_0[j] = pa_x[j] * tpx_xyz_yy_0[j] + 0.5 * fl1_fx * tpx_yz_yy_0[j] - fl1_fgb * fl1_fx * ts_xyz_yy_0[j];

            tpy_xxyz_yy_0[j] = pa_x[j] * tpy_xyz_yy_0[j] + 0.5 * fl1_fx * tpy_yz_yy_0[j];

            tpz_xxyz_yy_0[j] = pa_x[j] * tpz_xyz_yy_0[j] + 0.5 * fl1_fx * tpz_yz_yy_0[j];

            tpx_xxyz_yz_0[j] = pa_x[j] * tpx_xyz_yz_0[j] + 0.5 * fl1_fx * tpx_yz_yz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yz_0[j];

            tpy_xxyz_yz_0[j] = pa_x[j] * tpy_xyz_yz_0[j] + 0.5 * fl1_fx * tpy_yz_yz_0[j];

            tpz_xxyz_yz_0[j] = pa_x[j] * tpz_xyz_yz_0[j] + 0.5 * fl1_fx * tpz_yz_yz_0[j];

            tpx_xxyz_zz_0[j] = pa_x[j] * tpx_xyz_zz_0[j] + 0.5 * fl1_fx * tpx_yz_zz_0[j] - fl1_fgb * fl1_fx * ts_xyz_zz_0[j];

            tpy_xxyz_zz_0[j] = pa_x[j] * tpy_xyz_zz_0[j] + 0.5 * fl1_fx * tpy_yz_zz_0[j];

            tpz_xxyz_zz_0[j] = pa_x[j] * tpz_xyz_zz_0[j] + 0.5 * fl1_fx * tpz_yz_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD_90_135(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 30);

        auto tpy_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 31);

        auto tpy_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 32);

        auto tpy_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 33);

        auto tpy_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 34);

        auto tpy_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 35);

        auto tpy_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 35);

        auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36);

        auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37);

        auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38);

        auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39);

        auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40);

        auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41);

        auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42);

        auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43);

        auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44);

        auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 15);

        auto tpy_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 16);

        auto tpy_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 17);

        auto tpy_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto ts_xzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 30);

        auto ts_xzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 31);

        auto ts_xzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 32);

        auto ts_xzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 33);

        auto ts_xzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 34);

        auto ts_xzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 35);

        auto ts_yyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 36);

        auto ts_yyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 37);

        auto ts_yyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 38);

        auto ts_yyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 39);

        auto ts_yyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 40);

        auto ts_yyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 41);

        auto ts_yyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 42);

        auto ts_yyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 43);

        auto ts_yyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 44);

        // set up pointers to integrals

        auto tpx_xxzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 30);

        auto tpy_xxzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 30);

        auto tpz_xxzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 30);

        auto tpx_xxzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 31);

        auto tpy_xxzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 31);

        auto tpz_xxzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 31);

        auto tpx_xxzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 32);

        auto tpy_xxzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 32);

        auto tpz_xxzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 32);

        auto tpx_xxzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 33);

        auto tpy_xxzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 33);

        auto tpz_xxzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 33);

        auto tpx_xxzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 34);

        auto tpy_xxzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 34);

        auto tpz_xxzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 34);

        auto tpx_xxzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 35);

        auto tpy_xxzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 35);

        auto tpz_xxzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 35);

        auto tpx_xyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 36);

        auto tpy_xyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 36);

        auto tpz_xyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 36);

        auto tpx_xyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 37);

        auto tpy_xyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 37);

        auto tpz_xyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 37);

        auto tpx_xyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 38);

        auto tpy_xyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 38);

        auto tpz_xyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 38);

        auto tpx_xyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 39);

        auto tpy_xyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 39);

        auto tpz_xyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 39);

        auto tpx_xyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 40);

        auto tpy_xyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 40);

        auto tpz_xyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 40);

        auto tpx_xyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 41);

        auto tpy_xyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 41);

        auto tpz_xyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 41);

        auto tpx_xyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 42);

        auto tpy_xyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 42);

        auto tpz_xyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 42);

        auto tpx_xyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 43);

        auto tpy_xyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 43);

        auto tpz_xyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 43);

        auto tpx_xyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 44);

        auto tpy_xyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 44);

        auto tpz_xyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxzz_xx_0, tpx_xxzz_xy_0, tpx_xxzz_xz_0, \
                                     tpx_xxzz_yy_0, tpx_xxzz_yz_0, tpx_xxzz_zz_0, tpx_xyyy_xx_0, tpx_xyyy_xy_0, \
                                     tpx_xyyy_xz_0, tpx_xyyy_yy_0, tpx_xyyy_yz_0, tpx_xyyy_zz_0, tpx_xyyz_xx_0, \
                                     tpx_xyyz_xy_0, tpx_xyyz_xz_0, tpx_xzz_x_0, tpx_xzz_xx_0, tpx_xzz_xy_0, tpx_xzz_xz_0, \
                                     tpx_xzz_y_0, tpx_xzz_yy_0, tpx_xzz_yz_0, tpx_xzz_z_0, tpx_xzz_zz_0, tpx_yyy_x_0, \
                                     tpx_yyy_xx_0, tpx_yyy_xy_0, tpx_yyy_xz_0, tpx_yyy_y_0, tpx_yyy_yy_0, tpx_yyy_yz_0, \
                                     tpx_yyy_z_0, tpx_yyy_zz_0, tpx_yyz_x_0, tpx_yyz_xx_0, tpx_yyz_xy_0, tpx_yyz_xz_0, \
                                     tpx_yyz_y_0, tpx_yyz_z_0, tpx_zz_xx_0, tpx_zz_xy_0, tpx_zz_xz_0, tpx_zz_yy_0, \
                                     tpx_zz_yz_0, tpx_zz_zz_0, tpy_xxzz_xx_0, tpy_xxzz_xy_0, tpy_xxzz_xz_0, \
                                     tpy_xxzz_yy_0, tpy_xxzz_yz_0, tpy_xxzz_zz_0, tpy_xyyy_xx_0, tpy_xyyy_xy_0, \
                                     tpy_xyyy_xz_0, tpy_xyyy_yy_0, tpy_xyyy_yz_0, tpy_xyyy_zz_0, tpy_xyyz_xx_0, \
                                     tpy_xyyz_xy_0, tpy_xyyz_xz_0, tpy_xzz_x_0, tpy_xzz_xx_0, tpy_xzz_xy_0, tpy_xzz_xz_0, \
                                     tpy_xzz_y_0, tpy_xzz_yy_0, tpy_xzz_yz_0, tpy_xzz_z_0, tpy_xzz_zz_0, tpy_yyy_x_0, \
                                     tpy_yyy_xx_0, tpy_yyy_xy_0, tpy_yyy_xz_0, tpy_yyy_y_0, tpy_yyy_yy_0, tpy_yyy_yz_0, \
                                     tpy_yyy_z_0, tpy_yyy_zz_0, tpy_yyz_x_0, tpy_yyz_xx_0, tpy_yyz_xy_0, tpy_yyz_xz_0, \
                                     tpy_yyz_y_0, tpy_yyz_z_0, tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpy_zz_yy_0, \
                                     tpy_zz_yz_0, tpy_zz_zz_0, tpz_xxzz_xx_0, tpz_xxzz_xy_0, tpz_xxzz_xz_0, \
                                     tpz_xxzz_yy_0, tpz_xxzz_yz_0, tpz_xxzz_zz_0, tpz_xyyy_xx_0, tpz_xyyy_xy_0, \
                                     tpz_xyyy_xz_0, tpz_xyyy_yy_0, tpz_xyyy_yz_0, tpz_xyyy_zz_0, tpz_xyyz_xx_0, \
                                     tpz_xyyz_xy_0, tpz_xyyz_xz_0, tpz_xzz_x_0, tpz_xzz_xx_0, tpz_xzz_xy_0, tpz_xzz_xz_0, \
                                     tpz_xzz_y_0, tpz_xzz_yy_0, tpz_xzz_yz_0, tpz_xzz_z_0, tpz_xzz_zz_0, tpz_yyy_x_0, \
                                     tpz_yyy_xx_0, tpz_yyy_xy_0, tpz_yyy_xz_0, tpz_yyy_y_0, tpz_yyy_yy_0, tpz_yyy_yz_0, \
                                     tpz_yyy_z_0, tpz_yyy_zz_0, tpz_yyz_x_0, tpz_yyz_xx_0, tpz_yyz_xy_0, tpz_yyz_xz_0, \
                                     tpz_yyz_y_0, tpz_yyz_z_0, tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_yy_0, \
                                     tpz_zz_yz_0, tpz_zz_zz_0, ts_xzz_xx_0, ts_xzz_xy_0, ts_xzz_xz_0, ts_xzz_yy_0, \
                                     ts_xzz_yz_0, ts_xzz_zz_0, ts_yyy_xx_0, ts_yyy_xy_0, ts_yyy_xz_0, ts_yyy_yy_0, \
                                     ts_yyy_yz_0, ts_yyy_zz_0, ts_yyz_xx_0, ts_yyz_xy_0, ts_yyz_xz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxzz_xx_0[j] =
                pa_x[j] * tpx_xzz_xx_0[j] + 0.5 * fl1_fx * tpx_zz_xx_0[j] + fl1_fx * tpx_xzz_x_0[j] - fl1_fgb * fl1_fx * ts_xzz_xx_0[j];

            tpy_xxzz_xx_0[j] = pa_x[j] * tpy_xzz_xx_0[j] + 0.5 * fl1_fx * tpy_zz_xx_0[j] + fl1_fx * tpy_xzz_x_0[j];

            tpz_xxzz_xx_0[j] = pa_x[j] * tpz_xzz_xx_0[j] + 0.5 * fl1_fx * tpz_zz_xx_0[j] + fl1_fx * tpz_xzz_x_0[j];

            tpx_xxzz_xy_0[j] =
                pa_x[j] * tpx_xzz_xy_0[j] + 0.5 * fl1_fx * tpx_zz_xy_0[j] + 0.5 * fl1_fx * tpx_xzz_y_0[j] - fl1_fgb * fl1_fx * ts_xzz_xy_0[j];

            tpy_xxzz_xy_0[j] = pa_x[j] * tpy_xzz_xy_0[j] + 0.5 * fl1_fx * tpy_zz_xy_0[j] + 0.5 * fl1_fx * tpy_xzz_y_0[j];

            tpz_xxzz_xy_0[j] = pa_x[j] * tpz_xzz_xy_0[j] + 0.5 * fl1_fx * tpz_zz_xy_0[j] + 0.5 * fl1_fx * tpz_xzz_y_0[j];

            tpx_xxzz_xz_0[j] =
                pa_x[j] * tpx_xzz_xz_0[j] + 0.5 * fl1_fx * tpx_zz_xz_0[j] + 0.5 * fl1_fx * tpx_xzz_z_0[j] - fl1_fgb * fl1_fx * ts_xzz_xz_0[j];

            tpy_xxzz_xz_0[j] = pa_x[j] * tpy_xzz_xz_0[j] + 0.5 * fl1_fx * tpy_zz_xz_0[j] + 0.5 * fl1_fx * tpy_xzz_z_0[j];

            tpz_xxzz_xz_0[j] = pa_x[j] * tpz_xzz_xz_0[j] + 0.5 * fl1_fx * tpz_zz_xz_0[j] + 0.5 * fl1_fx * tpz_xzz_z_0[j];

            tpx_xxzz_yy_0[j] = pa_x[j] * tpx_xzz_yy_0[j] + 0.5 * fl1_fx * tpx_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_xzz_yy_0[j];

            tpy_xxzz_yy_0[j] = pa_x[j] * tpy_xzz_yy_0[j] + 0.5 * fl1_fx * tpy_zz_yy_0[j];

            tpz_xxzz_yy_0[j] = pa_x[j] * tpz_xzz_yy_0[j] + 0.5 * fl1_fx * tpz_zz_yy_0[j];

            tpx_xxzz_yz_0[j] = pa_x[j] * tpx_xzz_yz_0[j] + 0.5 * fl1_fx * tpx_zz_yz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yz_0[j];

            tpy_xxzz_yz_0[j] = pa_x[j] * tpy_xzz_yz_0[j] + 0.5 * fl1_fx * tpy_zz_yz_0[j];

            tpz_xxzz_yz_0[j] = pa_x[j] * tpz_xzz_yz_0[j] + 0.5 * fl1_fx * tpz_zz_yz_0[j];

            tpx_xxzz_zz_0[j] = pa_x[j] * tpx_xzz_zz_0[j] + 0.5 * fl1_fx * tpx_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_xzz_zz_0[j];

            tpy_xxzz_zz_0[j] = pa_x[j] * tpy_xzz_zz_0[j] + 0.5 * fl1_fx * tpy_zz_zz_0[j];

            tpz_xxzz_zz_0[j] = pa_x[j] * tpz_xzz_zz_0[j] + 0.5 * fl1_fx * tpz_zz_zz_0[j];

            tpx_xyyy_xx_0[j] = pa_x[j] * tpx_yyy_xx_0[j] + fl1_fx * tpx_yyy_x_0[j] - fl1_fgb * fl1_fx * ts_yyy_xx_0[j];

            tpy_xyyy_xx_0[j] = pa_x[j] * tpy_yyy_xx_0[j] + fl1_fx * tpy_yyy_x_0[j];

            tpz_xyyy_xx_0[j] = pa_x[j] * tpz_yyy_xx_0[j] + fl1_fx * tpz_yyy_x_0[j];

            tpx_xyyy_xy_0[j] = pa_x[j] * tpx_yyy_xy_0[j] + 0.5 * fl1_fx * tpx_yyy_y_0[j] - fl1_fgb * fl1_fx * ts_yyy_xy_0[j];

            tpy_xyyy_xy_0[j] = pa_x[j] * tpy_yyy_xy_0[j] + 0.5 * fl1_fx * tpy_yyy_y_0[j];

            tpz_xyyy_xy_0[j] = pa_x[j] * tpz_yyy_xy_0[j] + 0.5 * fl1_fx * tpz_yyy_y_0[j];

            tpx_xyyy_xz_0[j] = pa_x[j] * tpx_yyy_xz_0[j] + 0.5 * fl1_fx * tpx_yyy_z_0[j] - fl1_fgb * fl1_fx * ts_yyy_xz_0[j];

            tpy_xyyy_xz_0[j] = pa_x[j] * tpy_yyy_xz_0[j] + 0.5 * fl1_fx * tpy_yyy_z_0[j];

            tpz_xyyy_xz_0[j] = pa_x[j] * tpz_yyy_xz_0[j] + 0.5 * fl1_fx * tpz_yyy_z_0[j];

            tpx_xyyy_yy_0[j] = pa_x[j] * tpx_yyy_yy_0[j] - fl1_fgb * fl1_fx * ts_yyy_yy_0[j];

            tpy_xyyy_yy_0[j] = pa_x[j] * tpy_yyy_yy_0[j];

            tpz_xyyy_yy_0[j] = pa_x[j] * tpz_yyy_yy_0[j];

            tpx_xyyy_yz_0[j] = pa_x[j] * tpx_yyy_yz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yz_0[j];

            tpy_xyyy_yz_0[j] = pa_x[j] * tpy_yyy_yz_0[j];

            tpz_xyyy_yz_0[j] = pa_x[j] * tpz_yyy_yz_0[j];

            tpx_xyyy_zz_0[j] = pa_x[j] * tpx_yyy_zz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zz_0[j];

            tpy_xyyy_zz_0[j] = pa_x[j] * tpy_yyy_zz_0[j];

            tpz_xyyy_zz_0[j] = pa_x[j] * tpz_yyy_zz_0[j];

            tpx_xyyz_xx_0[j] = pa_x[j] * tpx_yyz_xx_0[j] + fl1_fx * tpx_yyz_x_0[j] - fl1_fgb * fl1_fx * ts_yyz_xx_0[j];

            tpy_xyyz_xx_0[j] = pa_x[j] * tpy_yyz_xx_0[j] + fl1_fx * tpy_yyz_x_0[j];

            tpz_xyyz_xx_0[j] = pa_x[j] * tpz_yyz_xx_0[j] + fl1_fx * tpz_yyz_x_0[j];

            tpx_xyyz_xy_0[j] = pa_x[j] * tpx_yyz_xy_0[j] + 0.5 * fl1_fx * tpx_yyz_y_0[j] - fl1_fgb * fl1_fx * ts_yyz_xy_0[j];

            tpy_xyyz_xy_0[j] = pa_x[j] * tpy_yyz_xy_0[j] + 0.5 * fl1_fx * tpy_yyz_y_0[j];

            tpz_xyyz_xy_0[j] = pa_x[j] * tpz_yyz_xy_0[j] + 0.5 * fl1_fx * tpz_yyz_y_0[j];

            tpx_xyyz_xz_0[j] = pa_x[j] * tpx_yyz_xz_0[j] + 0.5 * fl1_fx * tpx_yyz_z_0[j] - fl1_fgb * fl1_fx * ts_yyz_xz_0[j];

            tpy_xyyz_xz_0[j] = pa_x[j] * tpy_yyz_xz_0[j] + 0.5 * fl1_fx * tpy_yyz_z_0[j];

            tpz_xyyz_xz_0[j] = pa_x[j] * tpz_yyz_xz_0[j] + 0.5 * fl1_fx * tpz_yyz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD_135_180(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (135,180)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45);

        auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46);

        auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47);

        auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48);

        auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49);

        auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50);

        auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51);

        auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52);

        auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53);

        auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56);

        auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57);

        auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58);

        auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59);

        auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

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

        auto ts_yyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 45);

        auto ts_yyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 46);

        auto ts_yyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 47);

        auto ts_yzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 48);

        auto ts_yzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 49);

        auto ts_yzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 50);

        auto ts_yzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 51);

        auto ts_yzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 52);

        auto ts_yzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 53);

        auto ts_zzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 54);

        auto ts_zzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 55);

        auto ts_zzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 56);

        auto ts_zzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 57);

        auto ts_zzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 58);

        auto ts_zzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 59);

        // set up pointers to integrals

        auto tpx_xyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 45);

        auto tpy_xyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 45);

        auto tpz_xyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 45);

        auto tpx_xyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 46);

        auto tpy_xyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 46);

        auto tpz_xyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 46);

        auto tpx_xyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 47);

        auto tpy_xyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 47);

        auto tpz_xyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 47);

        auto tpx_xyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 48);

        auto tpy_xyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 48);

        auto tpz_xyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 48);

        auto tpx_xyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 49);

        auto tpy_xyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_xyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 49);

        auto tpx_xyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 50);

        auto tpy_xyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 50);

        auto tpz_xyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 50);

        auto tpx_xyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 51);

        auto tpy_xyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 51);

        auto tpz_xyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 51);

        auto tpx_xyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 52);

        auto tpy_xyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 52);

        auto tpz_xyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 52);

        auto tpx_xyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 53);

        auto tpy_xyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 53);

        auto tpz_xyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 53);

        auto tpx_xzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 54);

        auto tpy_xzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 54);

        auto tpz_xzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 54);

        auto tpx_xzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 55);

        auto tpy_xzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 55);

        auto tpz_xzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 55);

        auto tpx_xzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 56);

        auto tpy_xzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 56);

        auto tpz_xzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 56);

        auto tpx_xzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 57);

        auto tpy_xzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 57);

        auto tpz_xzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 57);

        auto tpx_xzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 58);

        auto tpy_xzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 58);

        auto tpz_xzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 58);

        auto tpx_xzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 59);

        auto tpy_xzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 59);

        auto tpz_xzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 59);

        // Batch of Integrals (135,180)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyyz_yy_0, tpx_xyyz_yz_0, tpx_xyyz_zz_0, \
                                     tpx_xyzz_xx_0, tpx_xyzz_xy_0, tpx_xyzz_xz_0, tpx_xyzz_yy_0, tpx_xyzz_yz_0, \
                                     tpx_xyzz_zz_0, tpx_xzzz_xx_0, tpx_xzzz_xy_0, tpx_xzzz_xz_0, tpx_xzzz_yy_0, \
                                     tpx_xzzz_yz_0, tpx_xzzz_zz_0, tpx_yyz_yy_0, tpx_yyz_yz_0, tpx_yyz_zz_0, tpx_yzz_x_0, \
                                     tpx_yzz_xx_0, tpx_yzz_xy_0, tpx_yzz_xz_0, tpx_yzz_y_0, tpx_yzz_yy_0, tpx_yzz_yz_0, \
                                     tpx_yzz_z_0, tpx_yzz_zz_0, tpx_zzz_x_0, tpx_zzz_xx_0, tpx_zzz_xy_0, tpx_zzz_xz_0, \
                                     tpx_zzz_y_0, tpx_zzz_yy_0, tpx_zzz_yz_0, tpx_zzz_z_0, tpx_zzz_zz_0, tpy_xyyz_yy_0, \
                                     tpy_xyyz_yz_0, tpy_xyyz_zz_0, tpy_xyzz_xx_0, tpy_xyzz_xy_0, tpy_xyzz_xz_0, \
                                     tpy_xyzz_yy_0, tpy_xyzz_yz_0, tpy_xyzz_zz_0, tpy_xzzz_xx_0, tpy_xzzz_xy_0, \
                                     tpy_xzzz_xz_0, tpy_xzzz_yy_0, tpy_xzzz_yz_0, tpy_xzzz_zz_0, tpy_yyz_yy_0, \
                                     tpy_yyz_yz_0, tpy_yyz_zz_0, tpy_yzz_x_0, tpy_yzz_xx_0, tpy_yzz_xy_0, tpy_yzz_xz_0, \
                                     tpy_yzz_y_0, tpy_yzz_yy_0, tpy_yzz_yz_0, tpy_yzz_z_0, tpy_yzz_zz_0, tpy_zzz_x_0, \
                                     tpy_zzz_xx_0, tpy_zzz_xy_0, tpy_zzz_xz_0, tpy_zzz_y_0, tpy_zzz_yy_0, tpy_zzz_yz_0, \
                                     tpy_zzz_z_0, tpy_zzz_zz_0, tpz_xyyz_yy_0, tpz_xyyz_yz_0, tpz_xyyz_zz_0, \
                                     tpz_xyzz_xx_0, tpz_xyzz_xy_0, tpz_xyzz_xz_0, tpz_xyzz_yy_0, tpz_xyzz_yz_0, \
                                     tpz_xyzz_zz_0, tpz_xzzz_xx_0, tpz_xzzz_xy_0, tpz_xzzz_xz_0, tpz_xzzz_yy_0, \
                                     tpz_xzzz_yz_0, tpz_xzzz_zz_0, tpz_yyz_yy_0, tpz_yyz_yz_0, tpz_yyz_zz_0, tpz_yzz_x_0, \
                                     tpz_yzz_xx_0, tpz_yzz_xy_0, tpz_yzz_xz_0, tpz_yzz_y_0, tpz_yzz_yy_0, tpz_yzz_yz_0, \
                                     tpz_yzz_z_0, tpz_yzz_zz_0, tpz_zzz_x_0, tpz_zzz_xx_0, tpz_zzz_xy_0, tpz_zzz_xz_0, \
                                     tpz_zzz_y_0, tpz_zzz_yy_0, tpz_zzz_yz_0, tpz_zzz_z_0, tpz_zzz_zz_0, ts_yyz_yy_0, \
                                     ts_yyz_yz_0, ts_yyz_zz_0, ts_yzz_xx_0, ts_yzz_xy_0, ts_yzz_xz_0, ts_yzz_yy_0, \
                                     ts_yzz_yz_0, ts_yzz_zz_0, ts_zzz_xx_0, ts_zzz_xy_0, ts_zzz_xz_0, ts_zzz_yy_0, \
                                     ts_zzz_yz_0, ts_zzz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xyyz_yy_0[j] = pa_x[j] * tpx_yyz_yy_0[j] - fl1_fgb * fl1_fx * ts_yyz_yy_0[j];

            tpy_xyyz_yy_0[j] = pa_x[j] * tpy_yyz_yy_0[j];

            tpz_xyyz_yy_0[j] = pa_x[j] * tpz_yyz_yy_0[j];

            tpx_xyyz_yz_0[j] = pa_x[j] * tpx_yyz_yz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yz_0[j];

            tpy_xyyz_yz_0[j] = pa_x[j] * tpy_yyz_yz_0[j];

            tpz_xyyz_yz_0[j] = pa_x[j] * tpz_yyz_yz_0[j];

            tpx_xyyz_zz_0[j] = pa_x[j] * tpx_yyz_zz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zz_0[j];

            tpy_xyyz_zz_0[j] = pa_x[j] * tpy_yyz_zz_0[j];

            tpz_xyyz_zz_0[j] = pa_x[j] * tpz_yyz_zz_0[j];

            tpx_xyzz_xx_0[j] = pa_x[j] * tpx_yzz_xx_0[j] + fl1_fx * tpx_yzz_x_0[j] - fl1_fgb * fl1_fx * ts_yzz_xx_0[j];

            tpy_xyzz_xx_0[j] = pa_x[j] * tpy_yzz_xx_0[j] + fl1_fx * tpy_yzz_x_0[j];

            tpz_xyzz_xx_0[j] = pa_x[j] * tpz_yzz_xx_0[j] + fl1_fx * tpz_yzz_x_0[j];

            tpx_xyzz_xy_0[j] = pa_x[j] * tpx_yzz_xy_0[j] + 0.5 * fl1_fx * tpx_yzz_y_0[j] - fl1_fgb * fl1_fx * ts_yzz_xy_0[j];

            tpy_xyzz_xy_0[j] = pa_x[j] * tpy_yzz_xy_0[j] + 0.5 * fl1_fx * tpy_yzz_y_0[j];

            tpz_xyzz_xy_0[j] = pa_x[j] * tpz_yzz_xy_0[j] + 0.5 * fl1_fx * tpz_yzz_y_0[j];

            tpx_xyzz_xz_0[j] = pa_x[j] * tpx_yzz_xz_0[j] + 0.5 * fl1_fx * tpx_yzz_z_0[j] - fl1_fgb * fl1_fx * ts_yzz_xz_0[j];

            tpy_xyzz_xz_0[j] = pa_x[j] * tpy_yzz_xz_0[j] + 0.5 * fl1_fx * tpy_yzz_z_0[j];

            tpz_xyzz_xz_0[j] = pa_x[j] * tpz_yzz_xz_0[j] + 0.5 * fl1_fx * tpz_yzz_z_0[j];

            tpx_xyzz_yy_0[j] = pa_x[j] * tpx_yzz_yy_0[j] - fl1_fgb * fl1_fx * ts_yzz_yy_0[j];

            tpy_xyzz_yy_0[j] = pa_x[j] * tpy_yzz_yy_0[j];

            tpz_xyzz_yy_0[j] = pa_x[j] * tpz_yzz_yy_0[j];

            tpx_xyzz_yz_0[j] = pa_x[j] * tpx_yzz_yz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yz_0[j];

            tpy_xyzz_yz_0[j] = pa_x[j] * tpy_yzz_yz_0[j];

            tpz_xyzz_yz_0[j] = pa_x[j] * tpz_yzz_yz_0[j];

            tpx_xyzz_zz_0[j] = pa_x[j] * tpx_yzz_zz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zz_0[j];

            tpy_xyzz_zz_0[j] = pa_x[j] * tpy_yzz_zz_0[j];

            tpz_xyzz_zz_0[j] = pa_x[j] * tpz_yzz_zz_0[j];

            tpx_xzzz_xx_0[j] = pa_x[j] * tpx_zzz_xx_0[j] + fl1_fx * tpx_zzz_x_0[j] - fl1_fgb * fl1_fx * ts_zzz_xx_0[j];

            tpy_xzzz_xx_0[j] = pa_x[j] * tpy_zzz_xx_0[j] + fl1_fx * tpy_zzz_x_0[j];

            tpz_xzzz_xx_0[j] = pa_x[j] * tpz_zzz_xx_0[j] + fl1_fx * tpz_zzz_x_0[j];

            tpx_xzzz_xy_0[j] = pa_x[j] * tpx_zzz_xy_0[j] + 0.5 * fl1_fx * tpx_zzz_y_0[j] - fl1_fgb * fl1_fx * ts_zzz_xy_0[j];

            tpy_xzzz_xy_0[j] = pa_x[j] * tpy_zzz_xy_0[j] + 0.5 * fl1_fx * tpy_zzz_y_0[j];

            tpz_xzzz_xy_0[j] = pa_x[j] * tpz_zzz_xy_0[j] + 0.5 * fl1_fx * tpz_zzz_y_0[j];

            tpx_xzzz_xz_0[j] = pa_x[j] * tpx_zzz_xz_0[j] + 0.5 * fl1_fx * tpx_zzz_z_0[j] - fl1_fgb * fl1_fx * ts_zzz_xz_0[j];

            tpy_xzzz_xz_0[j] = pa_x[j] * tpy_zzz_xz_0[j] + 0.5 * fl1_fx * tpy_zzz_z_0[j];

            tpz_xzzz_xz_0[j] = pa_x[j] * tpz_zzz_xz_0[j] + 0.5 * fl1_fx * tpz_zzz_z_0[j];

            tpx_xzzz_yy_0[j] = pa_x[j] * tpx_zzz_yy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yy_0[j];

            tpy_xzzz_yy_0[j] = pa_x[j] * tpy_zzz_yy_0[j];

            tpz_xzzz_yy_0[j] = pa_x[j] * tpz_zzz_yy_0[j];

            tpx_xzzz_yz_0[j] = pa_x[j] * tpx_zzz_yz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yz_0[j];

            tpy_xzzz_yz_0[j] = pa_x[j] * tpy_zzz_yz_0[j];

            tpz_xzzz_yz_0[j] = pa_x[j] * tpz_zzz_yz_0[j];

            tpx_xzzz_zz_0[j] = pa_x[j] * tpx_zzz_zz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zz_0[j];

            tpy_xzzz_zz_0[j] = pa_x[j] * tpy_zzz_zz_0[j];

            tpz_xzzz_zz_0[j] = pa_x[j] * tpz_zzz_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD_180_225(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (180,225)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36);

        auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37);

        auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38);

        auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39);

        auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40);

        auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41);

        auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42);

        auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43);

        auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44);

        auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45);

        auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46);

        auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47);

        auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48);

        auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49);

        auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50);

        auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto ts_yyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 36);

        auto ts_yyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 37);

        auto ts_yyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 38);

        auto ts_yyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 39);

        auto ts_yyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 40);

        auto ts_yyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 41);

        auto ts_yyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 42);

        auto ts_yyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 43);

        auto ts_yyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 44);

        auto ts_yyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 45);

        auto ts_yyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 46);

        auto ts_yyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 47);

        auto ts_yzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 48);

        auto ts_yzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 49);

        auto ts_yzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 50);

        // set up pointers to integrals

        auto tpx_yyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 60);

        auto tpy_yyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 60);

        auto tpz_yyyy_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 60);

        auto tpx_yyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 61);

        auto tpy_yyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 61);

        auto tpz_yyyy_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 61);

        auto tpx_yyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 62);

        auto tpy_yyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 62);

        auto tpz_yyyy_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 62);

        auto tpx_yyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 63);

        auto tpy_yyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 63);

        auto tpz_yyyy_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 63);

        auto tpx_yyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 64);

        auto tpy_yyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 64);

        auto tpz_yyyy_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 64);

        auto tpx_yyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 65);

        auto tpy_yyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yyyy_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 66);

        auto tpy_yyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 66);

        auto tpz_yyyz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 66);

        auto tpx_yyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 67);

        auto tpy_yyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 67);

        auto tpz_yyyz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 67);

        auto tpx_yyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 68);

        auto tpy_yyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 68);

        auto tpz_yyyz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 68);

        auto tpx_yyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 69);

        auto tpy_yyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 69);

        auto tpz_yyyz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 69);

        auto tpx_yyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 70);

        auto tpy_yyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 70);

        auto tpz_yyyz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 70);

        auto tpx_yyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 71);

        auto tpy_yyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 71);

        auto tpz_yyyz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 71);

        auto tpx_yyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 72);

        auto tpy_yyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 72);

        auto tpz_yyzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 72);

        auto tpx_yyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 73);

        auto tpy_yyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 73);

        auto tpz_yyzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 73);

        auto tpx_yyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 74);

        auto tpy_yyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 74);

        auto tpz_yyzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 74);

        // Batch of Integrals (180,225)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yy_xx_0, tpx_yy_xy_0, tpx_yy_xz_0, tpx_yy_yy_0, \
                                     tpx_yy_yz_0, tpx_yy_zz_0, tpx_yyy_x_0, tpx_yyy_xx_0, tpx_yyy_xy_0, tpx_yyy_xz_0, \
                                     tpx_yyy_y_0, tpx_yyy_yy_0, tpx_yyy_yz_0, tpx_yyy_z_0, tpx_yyy_zz_0, tpx_yyyy_xx_0, \
                                     tpx_yyyy_xy_0, tpx_yyyy_xz_0, tpx_yyyy_yy_0, tpx_yyyy_yz_0, tpx_yyyy_zz_0, \
                                     tpx_yyyz_xx_0, tpx_yyyz_xy_0, tpx_yyyz_xz_0, tpx_yyyz_yy_0, tpx_yyyz_yz_0, \
                                     tpx_yyyz_zz_0, tpx_yyz_x_0, tpx_yyz_xx_0, tpx_yyz_xy_0, tpx_yyz_xz_0, tpx_yyz_y_0, \
                                     tpx_yyz_yy_0, tpx_yyz_yz_0, tpx_yyz_z_0, tpx_yyz_zz_0, tpx_yyzz_xx_0, \
                                     tpx_yyzz_xy_0, tpx_yyzz_xz_0, tpx_yz_xx_0, tpx_yz_xy_0, tpx_yz_xz_0, tpx_yz_yy_0, \
                                     tpx_yz_yz_0, tpx_yz_zz_0, tpx_yzz_x_0, tpx_yzz_xx_0, tpx_yzz_xy_0, tpx_yzz_xz_0, \
                                     tpx_zz_xx_0, tpx_zz_xy_0, tpx_zz_xz_0, tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, \
                                     tpy_yy_yy_0, tpy_yy_yz_0, tpy_yy_zz_0, tpy_yyy_x_0, tpy_yyy_xx_0, tpy_yyy_xy_0, \
                                     tpy_yyy_xz_0, tpy_yyy_y_0, tpy_yyy_yy_0, tpy_yyy_yz_0, tpy_yyy_z_0, tpy_yyy_zz_0, \
                                     tpy_yyyy_xx_0, tpy_yyyy_xy_0, tpy_yyyy_xz_0, tpy_yyyy_yy_0, tpy_yyyy_yz_0, \
                                     tpy_yyyy_zz_0, tpy_yyyz_xx_0, tpy_yyyz_xy_0, tpy_yyyz_xz_0, tpy_yyyz_yy_0, \
                                     tpy_yyyz_yz_0, tpy_yyyz_zz_0, tpy_yyz_x_0, tpy_yyz_xx_0, tpy_yyz_xy_0, tpy_yyz_xz_0, \
                                     tpy_yyz_y_0, tpy_yyz_yy_0, tpy_yyz_yz_0, tpy_yyz_z_0, tpy_yyz_zz_0, tpy_yyzz_xx_0, \
                                     tpy_yyzz_xy_0, tpy_yyzz_xz_0, tpy_yz_xx_0, tpy_yz_xy_0, tpy_yz_xz_0, tpy_yz_yy_0, \
                                     tpy_yz_yz_0, tpy_yz_zz_0, tpy_yzz_x_0, tpy_yzz_xx_0, tpy_yzz_xy_0, tpy_yzz_xz_0, \
                                     tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpz_yy_xx_0, tpz_yy_xy_0, tpz_yy_xz_0, \
                                     tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_zz_0, tpz_yyy_x_0, tpz_yyy_xx_0, tpz_yyy_xy_0, \
                                     tpz_yyy_xz_0, tpz_yyy_y_0, tpz_yyy_yy_0, tpz_yyy_yz_0, tpz_yyy_z_0, tpz_yyy_zz_0, \
                                     tpz_yyyy_xx_0, tpz_yyyy_xy_0, tpz_yyyy_xz_0, tpz_yyyy_yy_0, tpz_yyyy_yz_0, \
                                     tpz_yyyy_zz_0, tpz_yyyz_xx_0, tpz_yyyz_xy_0, tpz_yyyz_xz_0, tpz_yyyz_yy_0, \
                                     tpz_yyyz_yz_0, tpz_yyyz_zz_0, tpz_yyz_x_0, tpz_yyz_xx_0, tpz_yyz_xy_0, tpz_yyz_xz_0, \
                                     tpz_yyz_y_0, tpz_yyz_yy_0, tpz_yyz_yz_0, tpz_yyz_z_0, tpz_yyz_zz_0, tpz_yyzz_xx_0, \
                                     tpz_yyzz_xy_0, tpz_yyzz_xz_0, tpz_yz_xx_0, tpz_yz_xy_0, tpz_yz_xz_0, tpz_yz_yy_0, \
                                     tpz_yz_yz_0, tpz_yz_zz_0, tpz_yzz_x_0, tpz_yzz_xx_0, tpz_yzz_xy_0, tpz_yzz_xz_0, \
                                     tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, ts_yyy_xx_0, ts_yyy_xy_0, ts_yyy_xz_0, \
                                     ts_yyy_yy_0, ts_yyy_yz_0, ts_yyy_zz_0, ts_yyz_xx_0, ts_yyz_xy_0, ts_yyz_xz_0, \
                                     ts_yyz_yy_0, ts_yyz_yz_0, ts_yyz_zz_0, ts_yzz_xx_0, ts_yzz_xy_0, ts_yzz_xz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyyy_xx_0[j] = pa_y[j] * tpx_yyy_xx_0[j] + 1.5 * fl1_fx * tpx_yy_xx_0[j];

            tpy_yyyy_xx_0[j] = pa_y[j] * tpy_yyy_xx_0[j] + 1.5 * fl1_fx * tpy_yy_xx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xx_0[j];

            tpz_yyyy_xx_0[j] = pa_y[j] * tpz_yyy_xx_0[j] + 1.5 * fl1_fx * tpz_yy_xx_0[j];

            tpx_yyyy_xy_0[j] = pa_y[j] * tpx_yyy_xy_0[j] + 1.5 * fl1_fx * tpx_yy_xy_0[j] + 0.5 * fl1_fx * tpx_yyy_x_0[j];

            tpy_yyyy_xy_0[j] =
                pa_y[j] * tpy_yyy_xy_0[j] + 1.5 * fl1_fx * tpy_yy_xy_0[j] + 0.5 * fl1_fx * tpy_yyy_x_0[j] - fl1_fgb * fl1_fx * ts_yyy_xy_0[j];

            tpz_yyyy_xy_0[j] = pa_y[j] * tpz_yyy_xy_0[j] + 1.5 * fl1_fx * tpz_yy_xy_0[j] + 0.5 * fl1_fx * tpz_yyy_x_0[j];

            tpx_yyyy_xz_0[j] = pa_y[j] * tpx_yyy_xz_0[j] + 1.5 * fl1_fx * tpx_yy_xz_0[j];

            tpy_yyyy_xz_0[j] = pa_y[j] * tpy_yyy_xz_0[j] + 1.5 * fl1_fx * tpy_yy_xz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xz_0[j];

            tpz_yyyy_xz_0[j] = pa_y[j] * tpz_yyy_xz_0[j] + 1.5 * fl1_fx * tpz_yy_xz_0[j];

            tpx_yyyy_yy_0[j] = pa_y[j] * tpx_yyy_yy_0[j] + 1.5 * fl1_fx * tpx_yy_yy_0[j] + fl1_fx * tpx_yyy_y_0[j];

            tpy_yyyy_yy_0[j] =
                pa_y[j] * tpy_yyy_yy_0[j] + 1.5 * fl1_fx * tpy_yy_yy_0[j] + fl1_fx * tpy_yyy_y_0[j] - fl1_fgb * fl1_fx * ts_yyy_yy_0[j];

            tpz_yyyy_yy_0[j] = pa_y[j] * tpz_yyy_yy_0[j] + 1.5 * fl1_fx * tpz_yy_yy_0[j] + fl1_fx * tpz_yyy_y_0[j];

            tpx_yyyy_yz_0[j] = pa_y[j] * tpx_yyy_yz_0[j] + 1.5 * fl1_fx * tpx_yy_yz_0[j] + 0.5 * fl1_fx * tpx_yyy_z_0[j];

            tpy_yyyy_yz_0[j] =
                pa_y[j] * tpy_yyy_yz_0[j] + 1.5 * fl1_fx * tpy_yy_yz_0[j] + 0.5 * fl1_fx * tpy_yyy_z_0[j] - fl1_fgb * fl1_fx * ts_yyy_yz_0[j];

            tpz_yyyy_yz_0[j] = pa_y[j] * tpz_yyy_yz_0[j] + 1.5 * fl1_fx * tpz_yy_yz_0[j] + 0.5 * fl1_fx * tpz_yyy_z_0[j];

            tpx_yyyy_zz_0[j] = pa_y[j] * tpx_yyy_zz_0[j] + 1.5 * fl1_fx * tpx_yy_zz_0[j];

            tpy_yyyy_zz_0[j] = pa_y[j] * tpy_yyy_zz_0[j] + 1.5 * fl1_fx * tpy_yy_zz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zz_0[j];

            tpz_yyyy_zz_0[j] = pa_y[j] * tpz_yyy_zz_0[j] + 1.5 * fl1_fx * tpz_yy_zz_0[j];

            tpx_yyyz_xx_0[j] = pa_y[j] * tpx_yyz_xx_0[j] + fl1_fx * tpx_yz_xx_0[j];

            tpy_yyyz_xx_0[j] = pa_y[j] * tpy_yyz_xx_0[j] + fl1_fx * tpy_yz_xx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xx_0[j];

            tpz_yyyz_xx_0[j] = pa_y[j] * tpz_yyz_xx_0[j] + fl1_fx * tpz_yz_xx_0[j];

            tpx_yyyz_xy_0[j] = pa_y[j] * tpx_yyz_xy_0[j] + fl1_fx * tpx_yz_xy_0[j] + 0.5 * fl1_fx * tpx_yyz_x_0[j];

            tpy_yyyz_xy_0[j] =
                pa_y[j] * tpy_yyz_xy_0[j] + fl1_fx * tpy_yz_xy_0[j] + 0.5 * fl1_fx * tpy_yyz_x_0[j] - fl1_fgb * fl1_fx * ts_yyz_xy_0[j];

            tpz_yyyz_xy_0[j] = pa_y[j] * tpz_yyz_xy_0[j] + fl1_fx * tpz_yz_xy_0[j] + 0.5 * fl1_fx * tpz_yyz_x_0[j];

            tpx_yyyz_xz_0[j] = pa_y[j] * tpx_yyz_xz_0[j] + fl1_fx * tpx_yz_xz_0[j];

            tpy_yyyz_xz_0[j] = pa_y[j] * tpy_yyz_xz_0[j] + fl1_fx * tpy_yz_xz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xz_0[j];

            tpz_yyyz_xz_0[j] = pa_y[j] * tpz_yyz_xz_0[j] + fl1_fx * tpz_yz_xz_0[j];

            tpx_yyyz_yy_0[j] = pa_y[j] * tpx_yyz_yy_0[j] + fl1_fx * tpx_yz_yy_0[j] + fl1_fx * tpx_yyz_y_0[j];

            tpy_yyyz_yy_0[j] = pa_y[j] * tpy_yyz_yy_0[j] + fl1_fx * tpy_yz_yy_0[j] + fl1_fx * tpy_yyz_y_0[j] - fl1_fgb * fl1_fx * ts_yyz_yy_0[j];

            tpz_yyyz_yy_0[j] = pa_y[j] * tpz_yyz_yy_0[j] + fl1_fx * tpz_yz_yy_0[j] + fl1_fx * tpz_yyz_y_0[j];

            tpx_yyyz_yz_0[j] = pa_y[j] * tpx_yyz_yz_0[j] + fl1_fx * tpx_yz_yz_0[j] + 0.5 * fl1_fx * tpx_yyz_z_0[j];

            tpy_yyyz_yz_0[j] =
                pa_y[j] * tpy_yyz_yz_0[j] + fl1_fx * tpy_yz_yz_0[j] + 0.5 * fl1_fx * tpy_yyz_z_0[j] - fl1_fgb * fl1_fx * ts_yyz_yz_0[j];

            tpz_yyyz_yz_0[j] = pa_y[j] * tpz_yyz_yz_0[j] + fl1_fx * tpz_yz_yz_0[j] + 0.5 * fl1_fx * tpz_yyz_z_0[j];

            tpx_yyyz_zz_0[j] = pa_y[j] * tpx_yyz_zz_0[j] + fl1_fx * tpx_yz_zz_0[j];

            tpy_yyyz_zz_0[j] = pa_y[j] * tpy_yyz_zz_0[j] + fl1_fx * tpy_yz_zz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zz_0[j];

            tpz_yyyz_zz_0[j] = pa_y[j] * tpz_yyz_zz_0[j] + fl1_fx * tpz_yz_zz_0[j];

            tpx_yyzz_xx_0[j] = pa_y[j] * tpx_yzz_xx_0[j] + 0.5 * fl1_fx * tpx_zz_xx_0[j];

            tpy_yyzz_xx_0[j] = pa_y[j] * tpy_yzz_xx_0[j] + 0.5 * fl1_fx * tpy_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xx_0[j];

            tpz_yyzz_xx_0[j] = pa_y[j] * tpz_yzz_xx_0[j] + 0.5 * fl1_fx * tpz_zz_xx_0[j];

            tpx_yyzz_xy_0[j] = pa_y[j] * tpx_yzz_xy_0[j] + 0.5 * fl1_fx * tpx_zz_xy_0[j] + 0.5 * fl1_fx * tpx_yzz_x_0[j];

            tpy_yyzz_xy_0[j] =
                pa_y[j] * tpy_yzz_xy_0[j] + 0.5 * fl1_fx * tpy_zz_xy_0[j] + 0.5 * fl1_fx * tpy_yzz_x_0[j] - fl1_fgb * fl1_fx * ts_yzz_xy_0[j];

            tpz_yyzz_xy_0[j] = pa_y[j] * tpz_yzz_xy_0[j] + 0.5 * fl1_fx * tpz_zz_xy_0[j] + 0.5 * fl1_fx * tpz_yzz_x_0[j];

            tpx_yyzz_xz_0[j] = pa_y[j] * tpx_yzz_xz_0[j] + 0.5 * fl1_fx * tpx_zz_xz_0[j];

            tpy_yyzz_xz_0[j] = pa_y[j] * tpy_yzz_xz_0[j] + 0.5 * fl1_fx * tpy_zz_xz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xz_0[j];

            tpz_yyzz_xz_0[j] = pa_y[j] * tpz_yzz_xz_0[j] + 0.5 * fl1_fx * tpz_zz_xz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGD_225_270(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (225,270)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51);

        auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52);

        auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53);

        auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56);

        auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57);

        auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58);

        auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59);

        auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

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

        auto ts_yzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 51);

        auto ts_yzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 52);

        auto ts_yzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 53);

        auto ts_zzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 54);

        auto ts_zzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 55);

        auto ts_zzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 56);

        auto ts_zzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 57);

        auto ts_zzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 58);

        auto ts_zzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 59);

        // set up pointers to integrals

        auto tpx_yyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 75);

        auto tpy_yyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_yyzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_yyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 76);

        auto tpy_yyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_yyzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_yyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 77);

        auto tpy_yyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_yyzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_yzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 78);

        auto tpy_yzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_yzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_yzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 79);

        auto tpy_yzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_yzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_yzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 80);

        auto tpy_yzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_yzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_yzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 81);

        auto tpy_yzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_yzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_yzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 82);

        auto tpy_yzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_yzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_yzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 83);

        auto tpy_yzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_yzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 84);

        auto tpy_zzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zzzz_xx_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 85);

        auto tpy_zzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zzzz_xy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 86);

        auto tpy_zzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zzzz_xz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 87);

        auto tpy_zzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zzzz_yy_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 87);

        auto tpx_zzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 88);

        auto tpy_zzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zzzz_yz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 88);

        auto tpx_zzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * idx + 89);

        auto tpy_zzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zzzz_zz_0 = primBuffer.data(pidx_p_4_2_m0 + 180 * bdim + 90 * idx + 89);

        // Batch of Integrals (225,270)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yyzz_yy_0, tpx_yyzz_yz_0, tpx_yyzz_zz_0, \
                                     tpx_yzz_y_0, tpx_yzz_yy_0, tpx_yzz_yz_0, tpx_yzz_z_0, tpx_yzz_zz_0, tpx_yzzz_xx_0, \
                                     tpx_yzzz_xy_0, tpx_yzzz_xz_0, tpx_yzzz_yy_0, tpx_yzzz_yz_0, tpx_yzzz_zz_0, \
                                     tpx_zz_xx_0, tpx_zz_xy_0, tpx_zz_xz_0, tpx_zz_yy_0, tpx_zz_yz_0, tpx_zz_zz_0, \
                                     tpx_zzz_x_0, tpx_zzz_xx_0, tpx_zzz_xy_0, tpx_zzz_xz_0, tpx_zzz_y_0, tpx_zzz_yy_0, \
                                     tpx_zzz_yz_0, tpx_zzz_z_0, tpx_zzz_zz_0, tpx_zzzz_xx_0, tpx_zzzz_xy_0, \
                                     tpx_zzzz_xz_0, tpx_zzzz_yy_0, tpx_zzzz_yz_0, tpx_zzzz_zz_0, tpy_yyzz_yy_0, \
                                     tpy_yyzz_yz_0, tpy_yyzz_zz_0, tpy_yzz_y_0, tpy_yzz_yy_0, tpy_yzz_yz_0, tpy_yzz_z_0, \
                                     tpy_yzz_zz_0, tpy_yzzz_xx_0, tpy_yzzz_xy_0, tpy_yzzz_xz_0, tpy_yzzz_yy_0, \
                                     tpy_yzzz_yz_0, tpy_yzzz_zz_0, tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpy_zz_yy_0, \
                                     tpy_zz_yz_0, tpy_zz_zz_0, tpy_zzz_x_0, tpy_zzz_xx_0, tpy_zzz_xy_0, tpy_zzz_xz_0, \
                                     tpy_zzz_y_0, tpy_zzz_yy_0, tpy_zzz_yz_0, tpy_zzz_z_0, tpy_zzz_zz_0, tpy_zzzz_xx_0, \
                                     tpy_zzzz_xy_0, tpy_zzzz_xz_0, tpy_zzzz_yy_0, tpy_zzzz_yz_0, tpy_zzzz_zz_0, \
                                     tpz_yyzz_yy_0, tpz_yyzz_yz_0, tpz_yyzz_zz_0, tpz_yzz_y_0, tpz_yzz_yy_0, \
                                     tpz_yzz_yz_0, tpz_yzz_z_0, tpz_yzz_zz_0, tpz_yzzz_xx_0, tpz_yzzz_xy_0, \
                                     tpz_yzzz_xz_0, tpz_yzzz_yy_0, tpz_yzzz_yz_0, tpz_yzzz_zz_0, tpz_zz_xx_0, \
                                     tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_yy_0, tpz_zz_yz_0, tpz_zz_zz_0, tpz_zzz_x_0, \
                                     tpz_zzz_xx_0, tpz_zzz_xy_0, tpz_zzz_xz_0, tpz_zzz_y_0, tpz_zzz_yy_0, tpz_zzz_yz_0, \
                                     tpz_zzz_z_0, tpz_zzz_zz_0, tpz_zzzz_xx_0, tpz_zzzz_xy_0, tpz_zzzz_xz_0, \
                                     tpz_zzzz_yy_0, tpz_zzzz_yz_0, tpz_zzzz_zz_0, ts_yzz_yy_0, ts_yzz_yz_0, ts_yzz_zz_0, \
                                     ts_zzz_xx_0, ts_zzz_xy_0, ts_zzz_xz_0, ts_zzz_yy_0, ts_zzz_yz_0, ts_zzz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyzz_yy_0[j] = pa_y[j] * tpx_yzz_yy_0[j] + 0.5 * fl1_fx * tpx_zz_yy_0[j] + fl1_fx * tpx_yzz_y_0[j];

            tpy_yyzz_yy_0[j] =
                pa_y[j] * tpy_yzz_yy_0[j] + 0.5 * fl1_fx * tpy_zz_yy_0[j] + fl1_fx * tpy_yzz_y_0[j] - fl1_fgb * fl1_fx * ts_yzz_yy_0[j];

            tpz_yyzz_yy_0[j] = pa_y[j] * tpz_yzz_yy_0[j] + 0.5 * fl1_fx * tpz_zz_yy_0[j] + fl1_fx * tpz_yzz_y_0[j];

            tpx_yyzz_yz_0[j] = pa_y[j] * tpx_yzz_yz_0[j] + 0.5 * fl1_fx * tpx_zz_yz_0[j] + 0.5 * fl1_fx * tpx_yzz_z_0[j];

            tpy_yyzz_yz_0[j] =
                pa_y[j] * tpy_yzz_yz_0[j] + 0.5 * fl1_fx * tpy_zz_yz_0[j] + 0.5 * fl1_fx * tpy_yzz_z_0[j] - fl1_fgb * fl1_fx * ts_yzz_yz_0[j];

            tpz_yyzz_yz_0[j] = pa_y[j] * tpz_yzz_yz_0[j] + 0.5 * fl1_fx * tpz_zz_yz_0[j] + 0.5 * fl1_fx * tpz_yzz_z_0[j];

            tpx_yyzz_zz_0[j] = pa_y[j] * tpx_yzz_zz_0[j] + 0.5 * fl1_fx * tpx_zz_zz_0[j];

            tpy_yyzz_zz_0[j] = pa_y[j] * tpy_yzz_zz_0[j] + 0.5 * fl1_fx * tpy_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zz_0[j];

            tpz_yyzz_zz_0[j] = pa_y[j] * tpz_yzz_zz_0[j] + 0.5 * fl1_fx * tpz_zz_zz_0[j];

            tpx_yzzz_xx_0[j] = pa_y[j] * tpx_zzz_xx_0[j];

            tpy_yzzz_xx_0[j] = pa_y[j] * tpy_zzz_xx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xx_0[j];

            tpz_yzzz_xx_0[j] = pa_y[j] * tpz_zzz_xx_0[j];

            tpx_yzzz_xy_0[j] = pa_y[j] * tpx_zzz_xy_0[j] + 0.5 * fl1_fx * tpx_zzz_x_0[j];

            tpy_yzzz_xy_0[j] = pa_y[j] * tpy_zzz_xy_0[j] + 0.5 * fl1_fx * tpy_zzz_x_0[j] - fl1_fgb * fl1_fx * ts_zzz_xy_0[j];

            tpz_yzzz_xy_0[j] = pa_y[j] * tpz_zzz_xy_0[j] + 0.5 * fl1_fx * tpz_zzz_x_0[j];

            tpx_yzzz_xz_0[j] = pa_y[j] * tpx_zzz_xz_0[j];

            tpy_yzzz_xz_0[j] = pa_y[j] * tpy_zzz_xz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xz_0[j];

            tpz_yzzz_xz_0[j] = pa_y[j] * tpz_zzz_xz_0[j];

            tpx_yzzz_yy_0[j] = pa_y[j] * tpx_zzz_yy_0[j] + fl1_fx * tpx_zzz_y_0[j];

            tpy_yzzz_yy_0[j] = pa_y[j] * tpy_zzz_yy_0[j] + fl1_fx * tpy_zzz_y_0[j] - fl1_fgb * fl1_fx * ts_zzz_yy_0[j];

            tpz_yzzz_yy_0[j] = pa_y[j] * tpz_zzz_yy_0[j] + fl1_fx * tpz_zzz_y_0[j];

            tpx_yzzz_yz_0[j] = pa_y[j] * tpx_zzz_yz_0[j] + 0.5 * fl1_fx * tpx_zzz_z_0[j];

            tpy_yzzz_yz_0[j] = pa_y[j] * tpy_zzz_yz_0[j] + 0.5 * fl1_fx * tpy_zzz_z_0[j] - fl1_fgb * fl1_fx * ts_zzz_yz_0[j];

            tpz_yzzz_yz_0[j] = pa_y[j] * tpz_zzz_yz_0[j] + 0.5 * fl1_fx * tpz_zzz_z_0[j];

            tpx_yzzz_zz_0[j] = pa_y[j] * tpx_zzz_zz_0[j];

            tpy_yzzz_zz_0[j] = pa_y[j] * tpy_zzz_zz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zz_0[j];

            tpz_yzzz_zz_0[j] = pa_y[j] * tpz_zzz_zz_0[j];

            tpx_zzzz_xx_0[j] = pa_z[j] * tpx_zzz_xx_0[j] + 1.5 * fl1_fx * tpx_zz_xx_0[j];

            tpy_zzzz_xx_0[j] = pa_z[j] * tpy_zzz_xx_0[j] + 1.5 * fl1_fx * tpy_zz_xx_0[j];

            tpz_zzzz_xx_0[j] = pa_z[j] * tpz_zzz_xx_0[j] + 1.5 * fl1_fx * tpz_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xx_0[j];

            tpx_zzzz_xy_0[j] = pa_z[j] * tpx_zzz_xy_0[j] + 1.5 * fl1_fx * tpx_zz_xy_0[j];

            tpy_zzzz_xy_0[j] = pa_z[j] * tpy_zzz_xy_0[j] + 1.5 * fl1_fx * tpy_zz_xy_0[j];

            tpz_zzzz_xy_0[j] = pa_z[j] * tpz_zzz_xy_0[j] + 1.5 * fl1_fx * tpz_zz_xy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xy_0[j];

            tpx_zzzz_xz_0[j] = pa_z[j] * tpx_zzz_xz_0[j] + 1.5 * fl1_fx * tpx_zz_xz_0[j] + 0.5 * fl1_fx * tpx_zzz_x_0[j];

            tpy_zzzz_xz_0[j] = pa_z[j] * tpy_zzz_xz_0[j] + 1.5 * fl1_fx * tpy_zz_xz_0[j] + 0.5 * fl1_fx * tpy_zzz_x_0[j];

            tpz_zzzz_xz_0[j] =
                pa_z[j] * tpz_zzz_xz_0[j] + 1.5 * fl1_fx * tpz_zz_xz_0[j] + 0.5 * fl1_fx * tpz_zzz_x_0[j] - fl1_fgb * fl1_fx * ts_zzz_xz_0[j];

            tpx_zzzz_yy_0[j] = pa_z[j] * tpx_zzz_yy_0[j] + 1.5 * fl1_fx * tpx_zz_yy_0[j];

            tpy_zzzz_yy_0[j] = pa_z[j] * tpy_zzz_yy_0[j] + 1.5 * fl1_fx * tpy_zz_yy_0[j];

            tpz_zzzz_yy_0[j] = pa_z[j] * tpz_zzz_yy_0[j] + 1.5 * fl1_fx * tpz_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yy_0[j];

            tpx_zzzz_yz_0[j] = pa_z[j] * tpx_zzz_yz_0[j] + 1.5 * fl1_fx * tpx_zz_yz_0[j] + 0.5 * fl1_fx * tpx_zzz_y_0[j];

            tpy_zzzz_yz_0[j] = pa_z[j] * tpy_zzz_yz_0[j] + 1.5 * fl1_fx * tpy_zz_yz_0[j] + 0.5 * fl1_fx * tpy_zzz_y_0[j];

            tpz_zzzz_yz_0[j] =
                pa_z[j] * tpz_zzz_yz_0[j] + 1.5 * fl1_fx * tpz_zz_yz_0[j] + 0.5 * fl1_fx * tpz_zzz_y_0[j] - fl1_fgb * fl1_fx * ts_zzz_yz_0[j];

            tpx_zzzz_zz_0[j] = pa_z[j] * tpx_zzz_zz_0[j] + 1.5 * fl1_fx * tpx_zz_zz_0[j] + fl1_fx * tpx_zzz_z_0[j];

            tpy_zzzz_zz_0[j] = pa_z[j] * tpy_zzz_zz_0[j] + 1.5 * fl1_fx * tpy_zz_zz_0[j] + fl1_fx * tpy_zzz_z_0[j];

            tpz_zzzz_zz_0[j] =
                pa_z[j] * tpz_zzz_zz_0[j] + 1.5 * fl1_fx * tpz_zz_zz_0[j] + fl1_fx * tpz_zzz_z_0[j] - fl1_fgb * fl1_fx * ts_zzz_zz_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
