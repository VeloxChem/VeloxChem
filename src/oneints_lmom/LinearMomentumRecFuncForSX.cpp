//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "LinearMomentumRecFuncForSX.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForSS(CMemBlock2D<double>&       primBuffer,
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

    // set up pointer to overlap integrals

    auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // set up pointer to linear momentum integrals

    auto pidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    if (pidx == -1) return;

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to ditances R(PA)

        auto pax = paDistances.data(3 * idx);

        auto pay = paDistances.data(3 * idx + 1);

        auto paz = paDistances.data(3 * idx + 2);

        // set up primitives buffer data

        auto fovl = primBuffer.data(sidx + idx);

        auto fmomx = primBuffer.data(pidx + idx);

        auto fmomy = primBuffer.data(pidx + bdim + idx);

        auto fmomz = primBuffer.data(pidx + 2 * bdim + idx);

        #pragma omp simd aligned(fovl, fmomx, fmomy, fmomz, fga: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fx = 2.0 * fga[j] * fovl[j];

            fmomx[j] = fx * pax[j];

            fmomy[j] = fx * pay[j];

            fmomz[j] = fx * paz[j];
        }

        idx++;
    }
}

void
compLinearMomentumForSP(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& pbDistances,
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

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_0_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx);

        // set up pointers to integrals

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tpx_0_0_0, tpx_0_x_0, tpx_0_y_0, tpx_0_z_0, \
                                     tpy_0_0_0, tpy_0_x_0, tpy_0_y_0, tpy_0_z_0, tpz_0_0_0, tpz_0_x_0, tpz_0_y_0, \
                                     tpz_0_z_0, ts_0_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tpx_0_x_0[j] = pb_x[j] * tpx_0_0_0[j] + fl1_fga * fl1_fx * ts_0_0_0[j];

            tpy_0_x_0[j] = pb_x[j] * tpy_0_0_0[j];

            tpz_0_x_0[j] = pb_x[j] * tpz_0_0_0[j];

            tpx_0_y_0[j] = pb_y[j] * tpx_0_0_0[j];

            tpy_0_y_0[j] = pb_y[j] * tpy_0_0_0[j] + fl1_fga * fl1_fx * ts_0_0_0[j];

            tpz_0_y_0[j] = pb_y[j] * tpz_0_0_0[j];

            tpx_0_z_0[j] = pb_z[j] * tpx_0_0_0[j];

            tpy_0_z_0[j] = pb_z[j] * tpy_0_0_0[j];

            tpz_0_z_0[j] = pb_z[j] * tpz_0_0_0[j] + fl1_fga * fl1_fx * ts_0_0_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPS(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx);

        // set up pointers to integrals

        auto tpx_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx);

        auto tpy_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx);

        auto tpz_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_0_0_0, tpx_x_0_0, tpx_y_0_0, tpx_z_0_0, \
                                     tpy_0_0_0, tpy_x_0_0, tpy_y_0_0, tpy_z_0_0, tpz_0_0_0, tpz_x_0_0, tpz_y_0_0, \
                                     tpz_z_0_0, ts_0_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_x_0_0[j] = pa_x[j] * tpx_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_0_0[j];

            tpy_x_0_0[j] = pa_x[j] * tpy_0_0_0[j];

            tpz_x_0_0[j] = pa_x[j] * tpz_0_0_0[j];

            tpx_y_0_0[j] = pa_y[j] * tpx_0_0_0[j];

            tpy_y_0_0[j] = pa_y[j] * tpy_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_0_0[j];

            tpz_y_0_0[j] = pa_y[j] * tpz_0_0_0[j];

            tpx_z_0_0[j] = pa_z[j] * tpx_0_0_0[j];

            tpy_z_0_0[j] = pa_z[j] * tpy_0_0_0[j];

            tpz_z_0_0[j] = pa_z[j] * tpz_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_0_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForSD(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& pbDistances,
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

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_0_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx);

        auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1);

        auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2);

        // set up pointers to integrals

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

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tpx_0_0_0, tpx_0_x_0, tpx_0_xx_0, tpx_0_xy_0, \
                                     tpx_0_xz_0, tpx_0_y_0, tpx_0_yy_0, tpx_0_yz_0, tpx_0_z_0, tpx_0_zz_0, tpy_0_0_0, \
                                     tpy_0_x_0, tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_y_0, tpy_0_yy_0, tpy_0_yz_0, \
                                     tpy_0_z_0, tpy_0_zz_0, tpz_0_0_0, tpz_0_x_0, tpz_0_xx_0, tpz_0_xy_0, tpz_0_xz_0, \
                                     tpz_0_y_0, tpz_0_yy_0, tpz_0_yz_0, tpz_0_z_0, tpz_0_zz_0, ts_0_x_0, ts_0_y_0, \
                                     ts_0_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tpx_0_xx_0[j] = pb_x[j] * tpx_0_x_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j] + fl1_fga * fl1_fx * ts_0_x_0[j];

            tpy_0_xx_0[j] = pb_x[j] * tpy_0_x_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_0_xx_0[j] = pb_x[j] * tpz_0_x_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_0_xy_0[j] = pb_x[j] * tpx_0_y_0[j] + fl1_fga * fl1_fx * ts_0_y_0[j];

            tpy_0_xy_0[j] = pb_x[j] * tpy_0_y_0[j];

            tpz_0_xy_0[j] = pb_x[j] * tpz_0_y_0[j];

            tpx_0_xz_0[j] = pb_x[j] * tpx_0_z_0[j] + fl1_fga * fl1_fx * ts_0_z_0[j];

            tpy_0_xz_0[j] = pb_x[j] * tpy_0_z_0[j];

            tpz_0_xz_0[j] = pb_x[j] * tpz_0_z_0[j];

            tpx_0_yy_0[j] = pb_y[j] * tpx_0_y_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_0_yy_0[j] = pb_y[j] * tpy_0_y_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j] + fl1_fga * fl1_fx * ts_0_y_0[j];

            tpz_0_yy_0[j] = pb_y[j] * tpz_0_y_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_0_yz_0[j] = pb_y[j] * tpx_0_z_0[j];

            tpy_0_yz_0[j] = pb_y[j] * tpy_0_z_0[j] + fl1_fga * fl1_fx * ts_0_z_0[j];

            tpz_0_yz_0[j] = pb_y[j] * tpz_0_z_0[j];

            tpx_0_zz_0[j] = pb_z[j] * tpx_0_z_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_0_zz_0[j] = pb_z[j] * tpy_0_z_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_0_zz_0[j] = pb_z[j] * tpz_0_z_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j] + fl1_fga * fl1_fx * ts_0_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDS(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx);

        auto tpy_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx);

        auto tpz_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx);

        auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1);

        auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2);

        // set up pointers to integrals

        auto tpx_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx);

        auto tpy_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx);

        auto tpz_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx);

        auto tpx_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 1);

        auto tpy_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 2);

        auto tpy_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 5);

        auto tpy_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 5);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_0_0_0, tpx_x_0_0, tpx_xx_0_0, tpx_xy_0_0, \
                                     tpx_xz_0_0, tpx_y_0_0, tpx_yy_0_0, tpx_yz_0_0, tpx_z_0_0, tpx_zz_0_0, tpy_0_0_0, \
                                     tpy_x_0_0, tpy_xx_0_0, tpy_xy_0_0, tpy_xz_0_0, tpy_y_0_0, tpy_yy_0_0, tpy_yz_0_0, \
                                     tpy_z_0_0, tpy_zz_0_0, tpz_0_0_0, tpz_x_0_0, tpz_xx_0_0, tpz_xy_0_0, tpz_xz_0_0, \
                                     tpz_y_0_0, tpz_yy_0_0, tpz_yz_0_0, tpz_z_0_0, tpz_zz_0_0, ts_x_0_0, ts_y_0_0, \
                                     ts_z_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xx_0_0[j] = pa_x[j] * tpx_x_0_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j] - fl1_fgb * fl1_fx * ts_x_0_0[j];

            tpy_xx_0_0[j] = pa_x[j] * tpy_x_0_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_xx_0_0[j] = pa_x[j] * tpz_x_0_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_xy_0_0[j] = pa_x[j] * tpx_y_0_0[j] - fl1_fgb * fl1_fx * ts_y_0_0[j];

            tpy_xy_0_0[j] = pa_x[j] * tpy_y_0_0[j];

            tpz_xy_0_0[j] = pa_x[j] * tpz_y_0_0[j];

            tpx_xz_0_0[j] = pa_x[j] * tpx_z_0_0[j] - fl1_fgb * fl1_fx * ts_z_0_0[j];

            tpy_xz_0_0[j] = pa_x[j] * tpy_z_0_0[j];

            tpz_xz_0_0[j] = pa_x[j] * tpz_z_0_0[j];

            tpx_yy_0_0[j] = pa_y[j] * tpx_y_0_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_yy_0_0[j] = pa_y[j] * tpy_y_0_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j] - fl1_fgb * fl1_fx * ts_y_0_0[j];

            tpz_yy_0_0[j] = pa_y[j] * tpz_y_0_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_yz_0_0[j] = pa_y[j] * tpx_z_0_0[j];

            tpy_yz_0_0[j] = pa_y[j] * tpy_z_0_0[j] - fl1_fgb * fl1_fx * ts_z_0_0[j];

            tpz_yz_0_0[j] = pa_y[j] * tpz_z_0_0[j];

            tpx_zz_0_0[j] = pa_z[j] * tpx_z_0_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_zz_0_0[j] = pa_z[j] * tpy_z_0_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_zz_0_0[j] = pa_z[j] * tpz_z_0_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j] - fl1_fgb * fl1_fx * ts_z_0_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForSF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& pbDistances,
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

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_0_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

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

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        // set up pointers to integrals

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

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tpx_0_x_0, tpx_0_xx_0, tpx_0_xxx_0, tpx_0_xxy_0, \
                                     tpx_0_xxz_0, tpx_0_xy_0, tpx_0_xyy_0, tpx_0_xyz_0, tpx_0_xz_0, tpx_0_xzz_0, \
                                     tpx_0_y_0, tpx_0_yy_0, tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yz_0, tpx_0_yzz_0, \
                                     tpx_0_z_0, tpx_0_zz_0, tpx_0_zzz_0, tpy_0_x_0, tpy_0_xx_0, tpy_0_xxx_0, \
                                     tpy_0_xxy_0, tpy_0_xxz_0, tpy_0_xy_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xz_0, \
                                     tpy_0_xzz_0, tpy_0_y_0, tpy_0_yy_0, tpy_0_yyy_0, tpy_0_yyz_0, tpy_0_yz_0, \
                                     tpy_0_yzz_0, tpy_0_z_0, tpy_0_zz_0, tpy_0_zzz_0, tpz_0_x_0, tpz_0_xx_0, tpz_0_xxx_0, \
                                     tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xy_0, tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xz_0, \
                                     tpz_0_xzz_0, tpz_0_y_0, tpz_0_yy_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yz_0, \
                                     tpz_0_yzz_0, tpz_0_z_0, tpz_0_zz_0, tpz_0_zzz_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, \
                                     ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tpx_0_xxx_0[j] = pb_x[j] * tpx_0_xx_0[j] + fl1_fx * tpx_0_x_0[j] + fl1_fga * fl1_fx * ts_0_xx_0[j];

            tpy_0_xxx_0[j] = pb_x[j] * tpy_0_xx_0[j] + fl1_fx * tpy_0_x_0[j];

            tpz_0_xxx_0[j] = pb_x[j] * tpz_0_xx_0[j] + fl1_fx * tpz_0_x_0[j];

            tpx_0_xxy_0[j] = pb_x[j] * tpx_0_xy_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] + fl1_fga * fl1_fx * ts_0_xy_0[j];

            tpy_0_xxy_0[j] = pb_x[j] * tpy_0_xy_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j];

            tpz_0_xxy_0[j] = pb_x[j] * tpz_0_xy_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j];

            tpx_0_xxz_0[j] = pb_x[j] * tpx_0_xz_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] + fl1_fga * fl1_fx * ts_0_xz_0[j];

            tpy_0_xxz_0[j] = pb_x[j] * tpy_0_xz_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j];

            tpz_0_xxz_0[j] = pb_x[j] * tpz_0_xz_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_0_xyy_0[j] = pb_x[j] * tpx_0_yy_0[j] + fl1_fga * fl1_fx * ts_0_yy_0[j];

            tpy_0_xyy_0[j] = pb_x[j] * tpy_0_yy_0[j];

            tpz_0_xyy_0[j] = pb_x[j] * tpz_0_yy_0[j];

            tpx_0_xyz_0[j] = pb_x[j] * tpx_0_yz_0[j] + fl1_fga * fl1_fx * ts_0_yz_0[j];

            tpy_0_xyz_0[j] = pb_x[j] * tpy_0_yz_0[j];

            tpz_0_xyz_0[j] = pb_x[j] * tpz_0_yz_0[j];

            tpx_0_xzz_0[j] = pb_x[j] * tpx_0_zz_0[j] + fl1_fga * fl1_fx * ts_0_zz_0[j];

            tpy_0_xzz_0[j] = pb_x[j] * tpy_0_zz_0[j];

            tpz_0_xzz_0[j] = pb_x[j] * tpz_0_zz_0[j];

            tpx_0_yyy_0[j] = pb_y[j] * tpx_0_yy_0[j] + fl1_fx * tpx_0_y_0[j];

            tpy_0_yyy_0[j] = pb_y[j] * tpy_0_yy_0[j] + fl1_fx * tpy_0_y_0[j] + fl1_fga * fl1_fx * ts_0_yy_0[j];

            tpz_0_yyy_0[j] = pb_y[j] * tpz_0_yy_0[j] + fl1_fx * tpz_0_y_0[j];

            tpx_0_yyz_0[j] = pb_y[j] * tpx_0_yz_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j];

            tpy_0_yyz_0[j] = pb_y[j] * tpy_0_yz_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] + fl1_fga * fl1_fx * ts_0_yz_0[j];

            tpz_0_yyz_0[j] = pb_y[j] * tpz_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_0_yzz_0[j] = pb_y[j] * tpx_0_zz_0[j];

            tpy_0_yzz_0[j] = pb_y[j] * tpy_0_zz_0[j] + fl1_fga * fl1_fx * ts_0_zz_0[j];

            tpz_0_yzz_0[j] = pb_y[j] * tpz_0_zz_0[j];

            tpx_0_zzz_0[j] = pb_z[j] * tpx_0_zz_0[j] + fl1_fx * tpx_0_z_0[j];

            tpy_0_zzz_0[j] = pb_z[j] * tpy_0_zz_0[j] + fl1_fx * tpy_0_z_0[j];

            tpz_0_zzz_0[j] = pb_z[j] * tpz_0_zz_0[j] + fl1_fx * tpz_0_z_0[j] + fl1_fga * fl1_fx * ts_0_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFS(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx);

        auto tpy_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx);

        auto tpz_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx);

        auto tpx_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 1);

        auto tpy_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 2);

        auto tpy_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 5);

        auto tpy_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto tpx_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx);

        auto tpy_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx);

        auto tpz_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx);

        auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1);

        auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2);

        auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3);

        auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4);

        auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto tpx_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx);

        auto tpy_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx);

        auto tpz_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx);

        auto tpx_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 1);

        auto tpy_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 2);

        auto tpy_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 3);

        auto tpy_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 4);

        auto tpy_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 5);

        auto tpy_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tpx_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 6);

        auto tpy_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tpx_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 7);

        auto tpy_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tpx_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 8);

        auto tpy_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tpx_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 9);

        auto tpy_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 9);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_x_0_0, tpx_xx_0_0, tpx_xxx_0_0, tpx_xxy_0_0, \
                                     tpx_xxz_0_0, tpx_xy_0_0, tpx_xyy_0_0, tpx_xyz_0_0, tpx_xz_0_0, tpx_xzz_0_0, \
                                     tpx_y_0_0, tpx_yy_0_0, tpx_yyy_0_0, tpx_yyz_0_0, tpx_yz_0_0, tpx_yzz_0_0, \
                                     tpx_z_0_0, tpx_zz_0_0, tpx_zzz_0_0, tpy_x_0_0, tpy_xx_0_0, tpy_xxx_0_0, \
                                     tpy_xxy_0_0, tpy_xxz_0_0, tpy_xy_0_0, tpy_xyy_0_0, tpy_xyz_0_0, tpy_xz_0_0, \
                                     tpy_xzz_0_0, tpy_y_0_0, tpy_yy_0_0, tpy_yyy_0_0, tpy_yyz_0_0, tpy_yz_0_0, \
                                     tpy_yzz_0_0, tpy_z_0_0, tpy_zz_0_0, tpy_zzz_0_0, tpz_x_0_0, tpz_xx_0_0, tpz_xxx_0_0, \
                                     tpz_xxy_0_0, tpz_xxz_0_0, tpz_xy_0_0, tpz_xyy_0_0, tpz_xyz_0_0, tpz_xz_0_0, \
                                     tpz_xzz_0_0, tpz_y_0_0, tpz_yy_0_0, tpz_yyy_0_0, tpz_yyz_0_0, tpz_yz_0_0, \
                                     tpz_yzz_0_0, tpz_z_0_0, tpz_zz_0_0, tpz_zzz_0_0, ts_xx_0_0, ts_xy_0_0, ts_xz_0_0, \
                                     ts_yy_0_0, ts_yz_0_0, ts_zz_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxx_0_0[j] = pa_x[j] * tpx_xx_0_0[j] + fl1_fx * tpx_x_0_0[j] - fl1_fgb * fl1_fx * ts_xx_0_0[j];

            tpy_xxx_0_0[j] = pa_x[j] * tpy_xx_0_0[j] + fl1_fx * tpy_x_0_0[j];

            tpz_xxx_0_0[j] = pa_x[j] * tpz_xx_0_0[j] + fl1_fx * tpz_x_0_0[j];

            tpx_xxy_0_0[j] = pa_x[j] * tpx_xy_0_0[j] + 0.5 * fl1_fx * tpx_y_0_0[j] - fl1_fgb * fl1_fx * ts_xy_0_0[j];

            tpy_xxy_0_0[j] = pa_x[j] * tpy_xy_0_0[j] + 0.5 * fl1_fx * tpy_y_0_0[j];

            tpz_xxy_0_0[j] = pa_x[j] * tpz_xy_0_0[j] + 0.5 * fl1_fx * tpz_y_0_0[j];

            tpx_xxz_0_0[j] = pa_x[j] * tpx_xz_0_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j] - fl1_fgb * fl1_fx * ts_xz_0_0[j];

            tpy_xxz_0_0[j] = pa_x[j] * tpy_xz_0_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j];

            tpz_xxz_0_0[j] = pa_x[j] * tpz_xz_0_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j];

            tpx_xyy_0_0[j] = pa_x[j] * tpx_yy_0_0[j] - fl1_fgb * fl1_fx * ts_yy_0_0[j];

            tpy_xyy_0_0[j] = pa_x[j] * tpy_yy_0_0[j];

            tpz_xyy_0_0[j] = pa_x[j] * tpz_yy_0_0[j];

            tpx_xyz_0_0[j] = pa_x[j] * tpx_yz_0_0[j] - fl1_fgb * fl1_fx * ts_yz_0_0[j];

            tpy_xyz_0_0[j] = pa_x[j] * tpy_yz_0_0[j];

            tpz_xyz_0_0[j] = pa_x[j] * tpz_yz_0_0[j];

            tpx_xzz_0_0[j] = pa_x[j] * tpx_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zz_0_0[j];

            tpy_xzz_0_0[j] = pa_x[j] * tpy_zz_0_0[j];

            tpz_xzz_0_0[j] = pa_x[j] * tpz_zz_0_0[j];

            tpx_yyy_0_0[j] = pa_y[j] * tpx_yy_0_0[j] + fl1_fx * tpx_y_0_0[j];

            tpy_yyy_0_0[j] = pa_y[j] * tpy_yy_0_0[j] + fl1_fx * tpy_y_0_0[j] - fl1_fgb * fl1_fx * ts_yy_0_0[j];

            tpz_yyy_0_0[j] = pa_y[j] * tpz_yy_0_0[j] + fl1_fx * tpz_y_0_0[j];

            tpx_yyz_0_0[j] = pa_y[j] * tpx_yz_0_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j];

            tpy_yyz_0_0[j] = pa_y[j] * tpy_yz_0_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j] - fl1_fgb * fl1_fx * ts_yz_0_0[j];

            tpz_yyz_0_0[j] = pa_y[j] * tpz_yz_0_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j];

            tpx_yzz_0_0[j] = pa_y[j] * tpx_zz_0_0[j];

            tpy_yzz_0_0[j] = pa_y[j] * tpy_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zz_0_0[j];

            tpz_yzz_0_0[j] = pa_y[j] * tpz_zz_0_0[j];

            tpx_zzz_0_0[j] = pa_z[j] * tpx_zz_0_0[j] + fl1_fx * tpx_z_0_0[j];

            tpy_zzz_0_0[j] = pa_z[j] * tpy_zz_0_0[j] + fl1_fx * tpy_z_0_0[j];

            tpz_zzz_0_0[j] = pa_z[j] * tpz_zz_0_0[j] + fl1_fx * tpz_z_0_0[j] - fl1_fgb * fl1_fx * ts_zz_0_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForSG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& pbDistances,
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

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_0_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

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

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tpx_0_xx_0, tpx_0_xxx_0, tpx_0_xxxx_0, \
                                     tpx_0_xxxy_0, tpx_0_xxxz_0, tpx_0_xxy_0, tpx_0_xxyy_0, tpx_0_xxyz_0, tpx_0_xxz_0, \
                                     tpx_0_xxzz_0, tpx_0_xy_0, tpx_0_xyy_0, tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyz_0, \
                                     tpx_0_xyzz_0, tpx_0_xz_0, tpx_0_xzz_0, tpx_0_xzzz_0, tpx_0_yy_0, tpx_0_yyy_0, \
                                     tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyz_0, tpx_0_yyzz_0, tpx_0_yz_0, tpx_0_yzz_0, \
                                     tpx_0_yzzz_0, tpx_0_zz_0, tpx_0_zzz_0, tpx_0_zzzz_0, tpy_0_xx_0, tpy_0_xxx_0, \
                                     tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxy_0, tpy_0_xxyy_0, tpy_0_xxyz_0, \
                                     tpy_0_xxz_0, tpy_0_xxzz_0, tpy_0_xy_0, tpy_0_xyy_0, tpy_0_xyyy_0, tpy_0_xyyz_0, \
                                     tpy_0_xyz_0, tpy_0_xyzz_0, tpy_0_xz_0, tpy_0_xzz_0, tpy_0_xzzz_0, tpy_0_yy_0, \
                                     tpy_0_yyy_0, tpy_0_yyyy_0, tpy_0_yyyz_0, tpy_0_yyz_0, tpy_0_yyzz_0, tpy_0_yz_0, \
                                     tpy_0_yzz_0, tpy_0_yzzz_0, tpy_0_zz_0, tpy_0_zzz_0, tpy_0_zzzz_0, tpz_0_xx_0, \
                                     tpz_0_xxx_0, tpz_0_xxxx_0, tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxy_0, tpz_0_xxyy_0, \
                                     tpz_0_xxyz_0, tpz_0_xxz_0, tpz_0_xxzz_0, tpz_0_xy_0, tpz_0_xyy_0, tpz_0_xyyy_0, \
                                     tpz_0_xyyz_0, tpz_0_xyz_0, tpz_0_xyzz_0, tpz_0_xz_0, tpz_0_xzz_0, tpz_0_xzzz_0, \
                                     tpz_0_yy_0, tpz_0_yyy_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyz_0, tpz_0_yyzz_0, \
                                     tpz_0_yz_0, tpz_0_yzz_0, tpz_0_yzzz_0, tpz_0_zz_0, tpz_0_zzz_0, tpz_0_zzzz_0, \
                                     ts_0_xxx_0, ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, \
                                     ts_0_yyz_0, ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tpx_0_xxxx_0[j] = pb_x[j] * tpx_0_xxx_0[j] + 1.5 * fl1_fx * tpx_0_xx_0[j] + fl1_fga * fl1_fx * ts_0_xxx_0[j];

            tpy_0_xxxx_0[j] = pb_x[j] * tpy_0_xxx_0[j] + 1.5 * fl1_fx * tpy_0_xx_0[j];

            tpz_0_xxxx_0[j] = pb_x[j] * tpz_0_xxx_0[j] + 1.5 * fl1_fx * tpz_0_xx_0[j];

            tpx_0_xxxy_0[j] = pb_x[j] * tpx_0_xxy_0[j] + fl1_fx * tpx_0_xy_0[j] + fl1_fga * fl1_fx * ts_0_xxy_0[j];

            tpy_0_xxxy_0[j] = pb_x[j] * tpy_0_xxy_0[j] + fl1_fx * tpy_0_xy_0[j];

            tpz_0_xxxy_0[j] = pb_x[j] * tpz_0_xxy_0[j] + fl1_fx * tpz_0_xy_0[j];

            tpx_0_xxxz_0[j] = pb_x[j] * tpx_0_xxz_0[j] + fl1_fx * tpx_0_xz_0[j] + fl1_fga * fl1_fx * ts_0_xxz_0[j];

            tpy_0_xxxz_0[j] = pb_x[j] * tpy_0_xxz_0[j] + fl1_fx * tpy_0_xz_0[j];

            tpz_0_xxxz_0[j] = pb_x[j] * tpz_0_xxz_0[j] + fl1_fx * tpz_0_xz_0[j];

            tpx_0_xxyy_0[j] = pb_x[j] * tpx_0_xyy_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] + fl1_fga * fl1_fx * ts_0_xyy_0[j];

            tpy_0_xxyy_0[j] = pb_x[j] * tpy_0_xyy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j];

            tpz_0_xxyy_0[j] = pb_x[j] * tpz_0_xyy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j];

            tpx_0_xxyz_0[j] = pb_x[j] * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] + fl1_fga * fl1_fx * ts_0_xyz_0[j];

            tpy_0_xxyz_0[j] = pb_x[j] * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j];

            tpz_0_xxyz_0[j] = pb_x[j] * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j];

            tpx_0_xxzz_0[j] = pb_x[j] * tpx_0_xzz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] + fl1_fga * fl1_fx * ts_0_xzz_0[j];

            tpy_0_xxzz_0[j] = pb_x[j] * tpy_0_xzz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j];

            tpz_0_xxzz_0[j] = pb_x[j] * tpz_0_xzz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];

            tpx_0_xyyy_0[j] = pb_x[j] * tpx_0_yyy_0[j] + fl1_fga * fl1_fx * ts_0_yyy_0[j];

            tpy_0_xyyy_0[j] = pb_x[j] * tpy_0_yyy_0[j];

            tpz_0_xyyy_0[j] = pb_x[j] * tpz_0_yyy_0[j];

            tpx_0_xyyz_0[j] = pb_x[j] * tpx_0_yyz_0[j] + fl1_fga * fl1_fx * ts_0_yyz_0[j];

            tpy_0_xyyz_0[j] = pb_x[j] * tpy_0_yyz_0[j];

            tpz_0_xyyz_0[j] = pb_x[j] * tpz_0_yyz_0[j];

            tpx_0_xyzz_0[j] = pb_x[j] * tpx_0_yzz_0[j] + fl1_fga * fl1_fx * ts_0_yzz_0[j];

            tpy_0_xyzz_0[j] = pb_x[j] * tpy_0_yzz_0[j];

            tpz_0_xyzz_0[j] = pb_x[j] * tpz_0_yzz_0[j];

            tpx_0_xzzz_0[j] = pb_x[j] * tpx_0_zzz_0[j] + fl1_fga * fl1_fx * ts_0_zzz_0[j];

            tpy_0_xzzz_0[j] = pb_x[j] * tpy_0_zzz_0[j];

            tpz_0_xzzz_0[j] = pb_x[j] * tpz_0_zzz_0[j];

            tpx_0_yyyy_0[j] = pb_y[j] * tpx_0_yyy_0[j] + 1.5 * fl1_fx * tpx_0_yy_0[j];

            tpy_0_yyyy_0[j] = pb_y[j] * tpy_0_yyy_0[j] + 1.5 * fl1_fx * tpy_0_yy_0[j] + fl1_fga * fl1_fx * ts_0_yyy_0[j];

            tpz_0_yyyy_0[j] = pb_y[j] * tpz_0_yyy_0[j] + 1.5 * fl1_fx * tpz_0_yy_0[j];

            tpx_0_yyyz_0[j] = pb_y[j] * tpx_0_yyz_0[j] + fl1_fx * tpx_0_yz_0[j];

            tpy_0_yyyz_0[j] = pb_y[j] * tpy_0_yyz_0[j] + fl1_fx * tpy_0_yz_0[j] + fl1_fga * fl1_fx * ts_0_yyz_0[j];

            tpz_0_yyyz_0[j] = pb_y[j] * tpz_0_yyz_0[j] + fl1_fx * tpz_0_yz_0[j];

            tpx_0_yyzz_0[j] = pb_y[j] * tpx_0_yzz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j];

            tpy_0_yyzz_0[j] = pb_y[j] * tpy_0_yzz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] + fl1_fga * fl1_fx * ts_0_yzz_0[j];

            tpz_0_yyzz_0[j] = pb_y[j] * tpz_0_yzz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];

            tpx_0_yzzz_0[j] = pb_y[j] * tpx_0_zzz_0[j];

            tpy_0_yzzz_0[j] = pb_y[j] * tpy_0_zzz_0[j] + fl1_fga * fl1_fx * ts_0_zzz_0[j];

            tpz_0_yzzz_0[j] = pb_y[j] * tpz_0_zzz_0[j];

            tpx_0_zzzz_0[j] = pb_z[j] * tpx_0_zzz_0[j] + 1.5 * fl1_fx * tpx_0_zz_0[j];

            tpy_0_zzzz_0[j] = pb_z[j] * tpy_0_zzz_0[j] + 1.5 * fl1_fx * tpy_0_zz_0[j];

            tpz_0_zzzz_0[j] = pb_z[j] * tpz_0_zzz_0[j] + 1.5 * fl1_fx * tpz_0_zz_0[j] + fl1_fga * fl1_fx * ts_0_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGS(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx);

        auto tpy_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx);

        auto tpz_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx);

        auto tpx_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 1);

        auto tpy_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 2);

        auto tpy_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 3);

        auto tpy_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 4);

        auto tpy_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 5);

        auto tpy_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tpx_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 6);

        auto tpy_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tpx_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 7);

        auto tpy_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tpx_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 8);

        auto tpy_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tpx_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 9);

        auto tpy_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 9);

        auto tpx_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx);

        auto tpy_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx);

        auto tpz_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx);

        auto tpx_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 1);

        auto tpy_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 2);

        auto tpy_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 5);

        auto tpy_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 5);

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

        auto tpx_xxxx_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx);

        auto tpy_xxxx_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx);

        auto tpz_xxxx_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx);

        auto tpx_xxxy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 1);

        auto tpy_xxxy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_xxxy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_xxxz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 2);

        auto tpy_xxxz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_xxxz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_xxyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 3);

        auto tpy_xxyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_xxyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_xxyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 4);

        auto tpy_xxyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_xxyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_xxzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 5);

        auto tpy_xxzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_xxzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_xyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 6);

        auto tpy_xyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_xyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_xyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 7);

        auto tpy_xyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_xyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_xyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 8);

        auto tpy_xyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_xyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_xzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 9);

        auto tpy_xzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_xzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_yyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 10);

        auto tpy_yyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_yyyy_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_yyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 11);

        auto tpy_yyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_yyyz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_yyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 12);

        auto tpy_yyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_yyzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_yzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 13);

        auto tpy_yzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_yzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_zzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * idx + 14);

        auto tpy_zzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_zzzz_0_0 = primBuffer.data(pidx_p_4_0_m0 + 30 * bdim + 15 * idx + 14);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_xx_0_0, tpx_xxx_0_0, tpx_xxxx_0_0, \
                                     tpx_xxxy_0_0, tpx_xxxz_0_0, tpx_xxy_0_0, tpx_xxyy_0_0, tpx_xxyz_0_0, tpx_xxz_0_0, \
                                     tpx_xxzz_0_0, tpx_xy_0_0, tpx_xyy_0_0, tpx_xyyy_0_0, tpx_xyyz_0_0, tpx_xyz_0_0, \
                                     tpx_xyzz_0_0, tpx_xz_0_0, tpx_xzz_0_0, tpx_xzzz_0_0, tpx_yy_0_0, tpx_yyy_0_0, \
                                     tpx_yyyy_0_0, tpx_yyyz_0_0, tpx_yyz_0_0, tpx_yyzz_0_0, tpx_yz_0_0, tpx_yzz_0_0, \
                                     tpx_yzzz_0_0, tpx_zz_0_0, tpx_zzz_0_0, tpx_zzzz_0_0, tpy_xx_0_0, tpy_xxx_0_0, \
                                     tpy_xxxx_0_0, tpy_xxxy_0_0, tpy_xxxz_0_0, tpy_xxy_0_0, tpy_xxyy_0_0, tpy_xxyz_0_0, \
                                     tpy_xxz_0_0, tpy_xxzz_0_0, tpy_xy_0_0, tpy_xyy_0_0, tpy_xyyy_0_0, tpy_xyyz_0_0, \
                                     tpy_xyz_0_0, tpy_xyzz_0_0, tpy_xz_0_0, tpy_xzz_0_0, tpy_xzzz_0_0, tpy_yy_0_0, \
                                     tpy_yyy_0_0, tpy_yyyy_0_0, tpy_yyyz_0_0, tpy_yyz_0_0, tpy_yyzz_0_0, tpy_yz_0_0, \
                                     tpy_yzz_0_0, tpy_yzzz_0_0, tpy_zz_0_0, tpy_zzz_0_0, tpy_zzzz_0_0, tpz_xx_0_0, \
                                     tpz_xxx_0_0, tpz_xxxx_0_0, tpz_xxxy_0_0, tpz_xxxz_0_0, tpz_xxy_0_0, tpz_xxyy_0_0, \
                                     tpz_xxyz_0_0, tpz_xxz_0_0, tpz_xxzz_0_0, tpz_xy_0_0, tpz_xyy_0_0, tpz_xyyy_0_0, \
                                     tpz_xyyz_0_0, tpz_xyz_0_0, tpz_xyzz_0_0, tpz_xz_0_0, tpz_xzz_0_0, tpz_xzzz_0_0, \
                                     tpz_yy_0_0, tpz_yyy_0_0, tpz_yyyy_0_0, tpz_yyyz_0_0, tpz_yyz_0_0, tpz_yyzz_0_0, \
                                     tpz_yz_0_0, tpz_yzz_0_0, tpz_yzzz_0_0, tpz_zz_0_0, tpz_zzz_0_0, tpz_zzzz_0_0, \
                                     ts_xxx_0_0, ts_xxy_0_0, ts_xxz_0_0, ts_xyy_0_0, ts_xyz_0_0, ts_xzz_0_0, ts_yyy_0_0, \
                                     ts_yyz_0_0, ts_yzz_0_0, ts_zzz_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxx_0_0[j] = pa_x[j] * tpx_xxx_0_0[j] + 1.5 * fl1_fx * tpx_xx_0_0[j] - fl1_fgb * fl1_fx * ts_xxx_0_0[j];

            tpy_xxxx_0_0[j] = pa_x[j] * tpy_xxx_0_0[j] + 1.5 * fl1_fx * tpy_xx_0_0[j];

            tpz_xxxx_0_0[j] = pa_x[j] * tpz_xxx_0_0[j] + 1.5 * fl1_fx * tpz_xx_0_0[j];

            tpx_xxxy_0_0[j] = pa_x[j] * tpx_xxy_0_0[j] + fl1_fx * tpx_xy_0_0[j] - fl1_fgb * fl1_fx * ts_xxy_0_0[j];

            tpy_xxxy_0_0[j] = pa_x[j] * tpy_xxy_0_0[j] + fl1_fx * tpy_xy_0_0[j];

            tpz_xxxy_0_0[j] = pa_x[j] * tpz_xxy_0_0[j] + fl1_fx * tpz_xy_0_0[j];

            tpx_xxxz_0_0[j] = pa_x[j] * tpx_xxz_0_0[j] + fl1_fx * tpx_xz_0_0[j] - fl1_fgb * fl1_fx * ts_xxz_0_0[j];

            tpy_xxxz_0_0[j] = pa_x[j] * tpy_xxz_0_0[j] + fl1_fx * tpy_xz_0_0[j];

            tpz_xxxz_0_0[j] = pa_x[j] * tpz_xxz_0_0[j] + fl1_fx * tpz_xz_0_0[j];

            tpx_xxyy_0_0[j] = pa_x[j] * tpx_xyy_0_0[j] + 0.5 * fl1_fx * tpx_yy_0_0[j] - fl1_fgb * fl1_fx * ts_xyy_0_0[j];

            tpy_xxyy_0_0[j] = pa_x[j] * tpy_xyy_0_0[j] + 0.5 * fl1_fx * tpy_yy_0_0[j];

            tpz_xxyy_0_0[j] = pa_x[j] * tpz_xyy_0_0[j] + 0.5 * fl1_fx * tpz_yy_0_0[j];

            tpx_xxyz_0_0[j] = pa_x[j] * tpx_xyz_0_0[j] + 0.5 * fl1_fx * tpx_yz_0_0[j] - fl1_fgb * fl1_fx * ts_xyz_0_0[j];

            tpy_xxyz_0_0[j] = pa_x[j] * tpy_xyz_0_0[j] + 0.5 * fl1_fx * tpy_yz_0_0[j];

            tpz_xxyz_0_0[j] = pa_x[j] * tpz_xyz_0_0[j] + 0.5 * fl1_fx * tpz_yz_0_0[j];

            tpx_xxzz_0_0[j] = pa_x[j] * tpx_xzz_0_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j] - fl1_fgb * fl1_fx * ts_xzz_0_0[j];

            tpy_xxzz_0_0[j] = pa_x[j] * tpy_xzz_0_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j];

            tpz_xxzz_0_0[j] = pa_x[j] * tpz_xzz_0_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j];

            tpx_xyyy_0_0[j] = pa_x[j] * tpx_yyy_0_0[j] - fl1_fgb * fl1_fx * ts_yyy_0_0[j];

            tpy_xyyy_0_0[j] = pa_x[j] * tpy_yyy_0_0[j];

            tpz_xyyy_0_0[j] = pa_x[j] * tpz_yyy_0_0[j];

            tpx_xyyz_0_0[j] = pa_x[j] * tpx_yyz_0_0[j] - fl1_fgb * fl1_fx * ts_yyz_0_0[j];

            tpy_xyyz_0_0[j] = pa_x[j] * tpy_yyz_0_0[j];

            tpz_xyyz_0_0[j] = pa_x[j] * tpz_yyz_0_0[j];

            tpx_xyzz_0_0[j] = pa_x[j] * tpx_yzz_0_0[j] - fl1_fgb * fl1_fx * ts_yzz_0_0[j];

            tpy_xyzz_0_0[j] = pa_x[j] * tpy_yzz_0_0[j];

            tpz_xyzz_0_0[j] = pa_x[j] * tpz_yzz_0_0[j];

            tpx_xzzz_0_0[j] = pa_x[j] * tpx_zzz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_0_0[j];

            tpy_xzzz_0_0[j] = pa_x[j] * tpy_zzz_0_0[j];

            tpz_xzzz_0_0[j] = pa_x[j] * tpz_zzz_0_0[j];

            tpx_yyyy_0_0[j] = pa_y[j] * tpx_yyy_0_0[j] + 1.5 * fl1_fx * tpx_yy_0_0[j];

            tpy_yyyy_0_0[j] = pa_y[j] * tpy_yyy_0_0[j] + 1.5 * fl1_fx * tpy_yy_0_0[j] - fl1_fgb * fl1_fx * ts_yyy_0_0[j];

            tpz_yyyy_0_0[j] = pa_y[j] * tpz_yyy_0_0[j] + 1.5 * fl1_fx * tpz_yy_0_0[j];

            tpx_yyyz_0_0[j] = pa_y[j] * tpx_yyz_0_0[j] + fl1_fx * tpx_yz_0_0[j];

            tpy_yyyz_0_0[j] = pa_y[j] * tpy_yyz_0_0[j] + fl1_fx * tpy_yz_0_0[j] - fl1_fgb * fl1_fx * ts_yyz_0_0[j];

            tpz_yyyz_0_0[j] = pa_y[j] * tpz_yyz_0_0[j] + fl1_fx * tpz_yz_0_0[j];

            tpx_yyzz_0_0[j] = pa_y[j] * tpx_yzz_0_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j];

            tpy_yyzz_0_0[j] = pa_y[j] * tpy_yzz_0_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j] - fl1_fgb * fl1_fx * ts_yzz_0_0[j];

            tpz_yyzz_0_0[j] = pa_y[j] * tpz_yzz_0_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j];

            tpx_yzzz_0_0[j] = pa_y[j] * tpx_zzz_0_0[j];

            tpy_yzzz_0_0[j] = pa_y[j] * tpy_zzz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_0_0[j];

            tpz_yzzz_0_0[j] = pa_y[j] * tpz_zzz_0_0[j];

            tpx_zzzz_0_0[j] = pa_z[j] * tpx_zzz_0_0[j] + 1.5 * fl1_fx * tpx_zz_0_0[j];

            tpy_zzzz_0_0[j] = pa_z[j] * tpy_zzz_0_0[j] + 1.5 * fl1_fx * tpy_zz_0_0[j];

            tpz_zzzz_0_0[j] = pa_z[j] * tpz_zzz_0_0[j] + 1.5 * fl1_fx * tpz_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_0_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
