//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AngularMomentumRecFuncForFG.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForFG(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForFG_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_300_350(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_350_400(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFG_400_450(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForFG_0_50(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (0,50)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx);

        auto tly_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx);

        auto tlz_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx);

        auto tlx_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 1);

        auto tly_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tlz_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tlx_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 2);

        auto tly_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tlz_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tlx_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 3);

        auto tly_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tlz_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tlx_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 4);

        auto tly_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tlz_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tlx_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 5);

        auto tly_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tlz_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tlx_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 6);

        auto tly_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tlz_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tlx_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 7);

        auto tly_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tlz_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tlx_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 8);

        auto tly_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tlz_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tlx_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 9);

        auto tly_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tlz_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tlx_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 10);

        auto tly_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tlz_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tlx_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 11);

        auto tly_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tlz_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tlx_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 12);

        auto tly_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tlz_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tlx_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 13);

        auto tly_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tlz_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tlx_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 14);

        auto tly_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tlz_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tlx_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 15);

        auto tly_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tlz_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tlx_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 16);

        auto tly_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 16);

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

        auto tlx_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 15);

        auto tly_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tlz_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tlx_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 16);

        auto tly_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tlx_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx);

        auto tly_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx);

        auto tlz_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx);

        auto tlx_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 1);

        auto tly_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tlz_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tlx_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 2);

        auto tly_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tlz_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tlx_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 3);

        auto tly_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tlz_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tlx_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 4);

        auto tly_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tlz_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tlx_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 5);

        auto tly_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tlz_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tlx_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 6);

        auto tly_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tlz_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tlx_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 7);

        auto tly_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tlz_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tlx_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 8);

        auto tly_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tlz_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tlx_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 9);

        auto tly_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tlz_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tlx_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 10);

        auto tly_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tlz_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tlx_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 11);

        auto tly_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto tpy_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx);

        auto tpz_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx);

        auto tpy_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tpz_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tpy_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tpz_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tpy_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tpz_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tpy_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tpz_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tpy_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tpz_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tpy_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tpz_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tpy_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tpz_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tpy_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tpz_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tpy_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tpz_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tpy_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tpz_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tpy_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tpz_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tpy_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tpz_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tpy_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tpz_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tpy_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tpz_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tpy_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tpz_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tpz_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tdy_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx);

        auto tdz_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx);

        auto tdy_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tdz_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tdy_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tdz_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tdy_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tdz_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tdy_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tdz_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tdy_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tdz_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tdy_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tdz_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tdy_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tdz_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tdy_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tdz_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tdy_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tdz_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tdy_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tdz_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tdy_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tdz_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tdy_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tdz_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tdy_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tdz_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tdy_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tdz_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tdy_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tdz_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tdz_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 16);

        // set up pointers to integrals

        auto tlx_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx);

        auto tly_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx);

        auto tlz_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx);

        auto tlx_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 1);

        auto tly_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tlz_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tlx_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 2);

        auto tly_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tlz_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tlx_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 3);

        auto tly_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tlz_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tlx_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 4);

        auto tly_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tlz_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tlx_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 5);

        auto tly_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tlz_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tlx_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 6);

        auto tly_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tlz_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tlx_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 7);

        auto tly_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tlz_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tlx_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 8);

        auto tly_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tlz_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tlx_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 9);

        auto tly_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tlz_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tlx_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 10);

        auto tly_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tlz_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tlx_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 11);

        auto tly_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tlz_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tlx_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 12);

        auto tly_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tlz_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tlx_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 13);

        auto tly_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tlz_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tlx_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 14);

        auto tly_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tlz_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tlx_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 15);

        auto tly_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tlz_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 15);

        auto tlx_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 16);

        auto tly_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xx_xxxx_0, tdy_xx_xxxy_0, tdy_xx_xxxz_0, \
                                     tdy_xx_xxyy_0, tdy_xx_xxyz_0, tdy_xx_xxzz_0, tdy_xx_xyyy_0, tdy_xx_xyyz_0, \
                                     tdy_xx_xyzz_0, tdy_xx_xzzz_0, tdy_xx_yyyy_0, tdy_xx_yyyz_0, tdy_xx_yyzz_0, \
                                     tdy_xx_yzzz_0, tdy_xx_zzzz_0, tdy_xy_xxxx_0, tdz_xx_xxxx_0, tdz_xx_xxxy_0, \
                                     tdz_xx_xxxz_0, tdz_xx_xxyy_0, tdz_xx_xxyz_0, tdz_xx_xxzz_0, tdz_xx_xyyy_0, \
                                     tdz_xx_xyyz_0, tdz_xx_xyzz_0, tdz_xx_xzzz_0, tdz_xx_yyyy_0, tdz_xx_yyyz_0, \
                                     tdz_xx_yyzz_0, tdz_xx_yzzz_0, tdz_xx_zzzz_0, tdz_xy_xxxx_0, tdz_xy_xxxy_0, \
                                     tlx_x_xxxx_0, tlx_x_xxxy_0, tlx_x_xxxz_0, tlx_x_xxyy_0, tlx_x_xxyz_0, tlx_x_xxzz_0, \
                                     tlx_x_xyyy_0, tlx_x_xyyz_0, tlx_x_xyzz_0, tlx_x_xzzz_0, tlx_x_yyyy_0, tlx_x_yyyz_0, \
                                     tlx_x_yyzz_0, tlx_x_yzzz_0, tlx_x_zzzz_0, tlx_xx_xxx_0, tlx_xx_xxxx_0, \
                                     tlx_xx_xxxy_0, tlx_xx_xxxz_0, tlx_xx_xxy_0, tlx_xx_xxyy_0, tlx_xx_xxyz_0, \
                                     tlx_xx_xxz_0, tlx_xx_xxzz_0, tlx_xx_xyy_0, tlx_xx_xyyy_0, tlx_xx_xyyz_0, \
                                     tlx_xx_xyz_0, tlx_xx_xyzz_0, tlx_xx_xzz_0, tlx_xx_xzzz_0, tlx_xx_yyy_0, \
                                     tlx_xx_yyyy_0, tlx_xx_yyyz_0, tlx_xx_yyz_0, tlx_xx_yyzz_0, tlx_xx_yzz_0, \
                                     tlx_xx_yzzz_0, tlx_xx_zzz_0, tlx_xx_zzzz_0, tlx_xxx_xxxx_0, tlx_xxx_xxxy_0, \
                                     tlx_xxx_xxxz_0, tlx_xxx_xxyy_0, tlx_xxx_xxyz_0, tlx_xxx_xxzz_0, tlx_xxx_xyyy_0, \
                                     tlx_xxx_xyyz_0, tlx_xxx_xyzz_0, tlx_xxx_xzzz_0, tlx_xxx_yyyy_0, tlx_xxx_yyyz_0, \
                                     tlx_xxx_yyzz_0, tlx_xxx_yzzz_0, tlx_xxx_zzzz_0, tlx_xxy_xxxx_0, tlx_xxy_xxxy_0, \
                                     tlx_xy_xxx_0, tlx_xy_xxxx_0, tlx_xy_xxxy_0, tlx_xy_xxy_0, tlx_y_xxxx_0, \
                                     tlx_y_xxxy_0, tly_x_xxxx_0, tly_x_xxxy_0, tly_x_xxxz_0, tly_x_xxyy_0, tly_x_xxyz_0, \
                                     tly_x_xxzz_0, tly_x_xyyy_0, tly_x_xyyz_0, tly_x_xyzz_0, tly_x_xzzz_0, tly_x_yyyy_0, \
                                     tly_x_yyyz_0, tly_x_yyzz_0, tly_x_yzzz_0, tly_x_zzzz_0, tly_xx_xxx_0, \
                                     tly_xx_xxxx_0, tly_xx_xxxy_0, tly_xx_xxxz_0, tly_xx_xxy_0, tly_xx_xxyy_0, \
                                     tly_xx_xxyz_0, tly_xx_xxz_0, tly_xx_xxzz_0, tly_xx_xyy_0, tly_xx_xyyy_0, \
                                     tly_xx_xyyz_0, tly_xx_xyz_0, tly_xx_xyzz_0, tly_xx_xzz_0, tly_xx_xzzz_0, \
                                     tly_xx_yyy_0, tly_xx_yyyy_0, tly_xx_yyyz_0, tly_xx_yyz_0, tly_xx_yyzz_0, \
                                     tly_xx_yzz_0, tly_xx_yzzz_0, tly_xx_zzz_0, tly_xx_zzzz_0, tly_xxx_xxxx_0, \
                                     tly_xxx_xxxy_0, tly_xxx_xxxz_0, tly_xxx_xxyy_0, tly_xxx_xxyz_0, tly_xxx_xxzz_0, \
                                     tly_xxx_xyyy_0, tly_xxx_xyyz_0, tly_xxx_xyzz_0, tly_xxx_xzzz_0, tly_xxx_yyyy_0, \
                                     tly_xxx_yyyz_0, tly_xxx_yyzz_0, tly_xxx_yzzz_0, tly_xxx_zzzz_0, tly_xxy_xxxx_0, \
                                     tly_xxy_xxxy_0, tly_xy_xxx_0, tly_xy_xxxx_0, tly_xy_xxxy_0, tly_xy_xxy_0, \
                                     tly_y_xxxx_0, tly_y_xxxy_0, tlz_x_xxxx_0, tlz_x_xxxy_0, tlz_x_xxxz_0, tlz_x_xxyy_0, \
                                     tlz_x_xxyz_0, tlz_x_xxzz_0, tlz_x_xyyy_0, tlz_x_xyyz_0, tlz_x_xyzz_0, tlz_x_xzzz_0, \
                                     tlz_x_yyyy_0, tlz_x_yyyz_0, tlz_x_yyzz_0, tlz_x_yzzz_0, tlz_x_zzzz_0, tlz_xx_xxx_0, \
                                     tlz_xx_xxxx_0, tlz_xx_xxxy_0, tlz_xx_xxxz_0, tlz_xx_xxy_0, tlz_xx_xxyy_0, \
                                     tlz_xx_xxyz_0, tlz_xx_xxz_0, tlz_xx_xxzz_0, tlz_xx_xyy_0, tlz_xx_xyyy_0, \
                                     tlz_xx_xyyz_0, tlz_xx_xyz_0, tlz_xx_xyzz_0, tlz_xx_xzz_0, tlz_xx_xzzz_0, \
                                     tlz_xx_yyy_0, tlz_xx_yyyy_0, tlz_xx_yyyz_0, tlz_xx_yyz_0, tlz_xx_yyzz_0, \
                                     tlz_xx_yzz_0, tlz_xx_yzzz_0, tlz_xx_zzz_0, tlz_xx_zzzz_0, tlz_xxx_xxxx_0, \
                                     tlz_xxx_xxxy_0, tlz_xxx_xxxz_0, tlz_xxx_xxyy_0, tlz_xxx_xxyz_0, tlz_xxx_xxzz_0, \
                                     tlz_xxx_xyyy_0, tlz_xxx_xyyz_0, tlz_xxx_xyzz_0, tlz_xxx_xzzz_0, tlz_xxx_yyyy_0, \
                                     tlz_xxx_yyyz_0, tlz_xxx_yyzz_0, tlz_xxx_yzzz_0, tlz_xxx_zzzz_0, tlz_xxy_xxxx_0, \
                                     tlz_xy_xxx_0, tlz_xy_xxxx_0, tlz_y_xxxx_0, tpy_xx_xxxx_0, tpy_xx_xxxy_0, \
                                     tpy_xx_xxxz_0, tpy_xx_xxyy_0, tpy_xx_xxyz_0, tpy_xx_xxzz_0, tpy_xx_xyyy_0, \
                                     tpy_xx_xyyz_0, tpy_xx_xyzz_0, tpy_xx_xzzz_0, tpy_xx_yyyy_0, tpy_xx_yyyz_0, \
                                     tpy_xx_yyzz_0, tpy_xx_yzzz_0, tpy_xx_zzzz_0, tpy_xy_xxxx_0, tpz_xx_xxxx_0, \
                                     tpz_xx_xxxy_0, tpz_xx_xxxz_0, tpz_xx_xxyy_0, tpz_xx_xxyz_0, tpz_xx_xxzz_0, \
                                     tpz_xx_xyyy_0, tpz_xx_xyyz_0, tpz_xx_xyzz_0, tpz_xx_xzzz_0, tpz_xx_yyyy_0, \
                                     tpz_xx_yyyz_0, tpz_xx_yyzz_0, tpz_xx_yzzz_0, tpz_xx_zzzz_0, tpz_xy_xxxx_0, \
                                     tpz_xy_xxxy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxx_xxxx_0[j] = pa_x[j] * tlx_xx_xxxx_0[j] + fl1_fx * tlx_x_xxxx_0[j] + 2.0 * fl1_fx * tlx_xx_xxx_0[j];

            tly_xxx_xxxx_0[j] = pa_x[j] * tly_xx_xxxx_0[j] + fl1_fx * tly_x_xxxx_0[j] + 2.0 * fl1_fx * tly_xx_xxx_0[j] +
                                0.5 * fl1_fx * tpz_xx_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xx_xxxx_0[j];

            tlz_xxx_xxxx_0[j] = pa_x[j] * tlz_xx_xxxx_0[j] + fl1_fx * tlz_x_xxxx_0[j] + 2.0 * fl1_fx * tlz_xx_xxx_0[j] -
                                0.5 * fl1_fx * tpy_xx_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xx_xxxx_0[j];

            tlx_xxx_xxxy_0[j] = pa_x[j] * tlx_xx_xxxy_0[j] + fl1_fx * tlx_x_xxxy_0[j] + 1.5 * fl1_fx * tlx_xx_xxy_0[j];

            tly_xxx_xxxy_0[j] = pa_x[j] * tly_xx_xxxy_0[j] + fl1_fx * tly_x_xxxy_0[j] + 1.5 * fl1_fx * tly_xx_xxy_0[j] +
                                0.5 * fl1_fx * tpz_xx_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xx_xxxy_0[j];

            tlz_xxx_xxxy_0[j] = pa_x[j] * tlz_xx_xxxy_0[j] + fl1_fx * tlz_x_xxxy_0[j] + 1.5 * fl1_fx * tlz_xx_xxy_0[j] -
                                0.5 * fl1_fx * tpy_xx_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xx_xxxy_0[j];

            tlx_xxx_xxxz_0[j] = pa_x[j] * tlx_xx_xxxz_0[j] + fl1_fx * tlx_x_xxxz_0[j] + 1.5 * fl1_fx * tlx_xx_xxz_0[j];

            tly_xxx_xxxz_0[j] = pa_x[j] * tly_xx_xxxz_0[j] + fl1_fx * tly_x_xxxz_0[j] + 1.5 * fl1_fx * tly_xx_xxz_0[j] +
                                0.5 * fl1_fx * tpz_xx_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xx_xxxz_0[j];

            tlz_xxx_xxxz_0[j] = pa_x[j] * tlz_xx_xxxz_0[j] + fl1_fx * tlz_x_xxxz_0[j] + 1.5 * fl1_fx * tlz_xx_xxz_0[j] -
                                0.5 * fl1_fx * tpy_xx_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xx_xxxz_0[j];

            tlx_xxx_xxyy_0[j] = pa_x[j] * tlx_xx_xxyy_0[j] + fl1_fx * tlx_x_xxyy_0[j] + fl1_fx * tlx_xx_xyy_0[j];

            tly_xxx_xxyy_0[j] = pa_x[j] * tly_xx_xxyy_0[j] + fl1_fx * tly_x_xxyy_0[j] + fl1_fx * tly_xx_xyy_0[j] + 0.5 * fl1_fx * tpz_xx_xxyy_0[j] +
                                fl1_fx * fl1_fgb * tdz_xx_xxyy_0[j];

            tlz_xxx_xxyy_0[j] = pa_x[j] * tlz_xx_xxyy_0[j] + fl1_fx * tlz_x_xxyy_0[j] + fl1_fx * tlz_xx_xyy_0[j] - 0.5 * fl1_fx * tpy_xx_xxyy_0[j] -
                                fl1_fx * fl1_fgb * tdy_xx_xxyy_0[j];

            tlx_xxx_xxyz_0[j] = pa_x[j] * tlx_xx_xxyz_0[j] + fl1_fx * tlx_x_xxyz_0[j] + fl1_fx * tlx_xx_xyz_0[j];

            tly_xxx_xxyz_0[j] = pa_x[j] * tly_xx_xxyz_0[j] + fl1_fx * tly_x_xxyz_0[j] + fl1_fx * tly_xx_xyz_0[j] + 0.5 * fl1_fx * tpz_xx_xxyz_0[j] +
                                fl1_fx * fl1_fgb * tdz_xx_xxyz_0[j];

            tlz_xxx_xxyz_0[j] = pa_x[j] * tlz_xx_xxyz_0[j] + fl1_fx * tlz_x_xxyz_0[j] + fl1_fx * tlz_xx_xyz_0[j] - 0.5 * fl1_fx * tpy_xx_xxyz_0[j] -
                                fl1_fx * fl1_fgb * tdy_xx_xxyz_0[j];

            tlx_xxx_xxzz_0[j] = pa_x[j] * tlx_xx_xxzz_0[j] + fl1_fx * tlx_x_xxzz_0[j] + fl1_fx * tlx_xx_xzz_0[j];

            tly_xxx_xxzz_0[j] = pa_x[j] * tly_xx_xxzz_0[j] + fl1_fx * tly_x_xxzz_0[j] + fl1_fx * tly_xx_xzz_0[j] + 0.5 * fl1_fx * tpz_xx_xxzz_0[j] +
                                fl1_fx * fl1_fgb * tdz_xx_xxzz_0[j];

            tlz_xxx_xxzz_0[j] = pa_x[j] * tlz_xx_xxzz_0[j] + fl1_fx * tlz_x_xxzz_0[j] + fl1_fx * tlz_xx_xzz_0[j] - 0.5 * fl1_fx * tpy_xx_xxzz_0[j] -
                                fl1_fx * fl1_fgb * tdy_xx_xxzz_0[j];

            tlx_xxx_xyyy_0[j] = pa_x[j] * tlx_xx_xyyy_0[j] + fl1_fx * tlx_x_xyyy_0[j] + 0.5 * fl1_fx * tlx_xx_yyy_0[j];

            tly_xxx_xyyy_0[j] = pa_x[j] * tly_xx_xyyy_0[j] + fl1_fx * tly_x_xyyy_0[j] + 0.5 * fl1_fx * tly_xx_yyy_0[j] +
                                0.5 * fl1_fx * tpz_xx_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xx_xyyy_0[j];

            tlz_xxx_xyyy_0[j] = pa_x[j] * tlz_xx_xyyy_0[j] + fl1_fx * tlz_x_xyyy_0[j] + 0.5 * fl1_fx * tlz_xx_yyy_0[j] -
                                0.5 * fl1_fx * tpy_xx_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xx_xyyy_0[j];

            tlx_xxx_xyyz_0[j] = pa_x[j] * tlx_xx_xyyz_0[j] + fl1_fx * tlx_x_xyyz_0[j] + 0.5 * fl1_fx * tlx_xx_yyz_0[j];

            tly_xxx_xyyz_0[j] = pa_x[j] * tly_xx_xyyz_0[j] + fl1_fx * tly_x_xyyz_0[j] + 0.5 * fl1_fx * tly_xx_yyz_0[j] +
                                0.5 * fl1_fx * tpz_xx_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xx_xyyz_0[j];

            tlz_xxx_xyyz_0[j] = pa_x[j] * tlz_xx_xyyz_0[j] + fl1_fx * tlz_x_xyyz_0[j] + 0.5 * fl1_fx * tlz_xx_yyz_0[j] -
                                0.5 * fl1_fx * tpy_xx_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xx_xyyz_0[j];

            tlx_xxx_xyzz_0[j] = pa_x[j] * tlx_xx_xyzz_0[j] + fl1_fx * tlx_x_xyzz_0[j] + 0.5 * fl1_fx * tlx_xx_yzz_0[j];

            tly_xxx_xyzz_0[j] = pa_x[j] * tly_xx_xyzz_0[j] + fl1_fx * tly_x_xyzz_0[j] + 0.5 * fl1_fx * tly_xx_yzz_0[j] +
                                0.5 * fl1_fx * tpz_xx_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_xyzz_0[j];

            tlz_xxx_xyzz_0[j] = pa_x[j] * tlz_xx_xyzz_0[j] + fl1_fx * tlz_x_xyzz_0[j] + 0.5 * fl1_fx * tlz_xx_yzz_0[j] -
                                0.5 * fl1_fx * tpy_xx_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_xyzz_0[j];

            tlx_xxx_xzzz_0[j] = pa_x[j] * tlx_xx_xzzz_0[j] + fl1_fx * tlx_x_xzzz_0[j] + 0.5 * fl1_fx * tlx_xx_zzz_0[j];

            tly_xxx_xzzz_0[j] = pa_x[j] * tly_xx_xzzz_0[j] + fl1_fx * tly_x_xzzz_0[j] + 0.5 * fl1_fx * tly_xx_zzz_0[j] +
                                0.5 * fl1_fx * tpz_xx_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_xzzz_0[j];

            tlz_xxx_xzzz_0[j] = pa_x[j] * tlz_xx_xzzz_0[j] + fl1_fx * tlz_x_xzzz_0[j] + 0.5 * fl1_fx * tlz_xx_zzz_0[j] -
                                0.5 * fl1_fx * tpy_xx_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_xzzz_0[j];

            tlx_xxx_yyyy_0[j] = pa_x[j] * tlx_xx_yyyy_0[j] + fl1_fx * tlx_x_yyyy_0[j];

            tly_xxx_yyyy_0[j] =
                pa_x[j] * tly_xx_yyyy_0[j] + fl1_fx * tly_x_yyyy_0[j] + 0.5 * fl1_fx * tpz_xx_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_xx_yyyy_0[j];

            tlz_xxx_yyyy_0[j] =
                pa_x[j] * tlz_xx_yyyy_0[j] + fl1_fx * tlz_x_yyyy_0[j] - 0.5 * fl1_fx * tpy_xx_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_xx_yyyy_0[j];

            tlx_xxx_yyyz_0[j] = pa_x[j] * tlx_xx_yyyz_0[j] + fl1_fx * tlx_x_yyyz_0[j];

            tly_xxx_yyyz_0[j] =
                pa_x[j] * tly_xx_yyyz_0[j] + fl1_fx * tly_x_yyyz_0[j] + 0.5 * fl1_fx * tpz_xx_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yyyz_0[j];

            tlz_xxx_yyyz_0[j] =
                pa_x[j] * tlz_xx_yyyz_0[j] + fl1_fx * tlz_x_yyyz_0[j] - 0.5 * fl1_fx * tpy_xx_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yyyz_0[j];

            tlx_xxx_yyzz_0[j] = pa_x[j] * tlx_xx_yyzz_0[j] + fl1_fx * tlx_x_yyzz_0[j];

            tly_xxx_yyzz_0[j] =
                pa_x[j] * tly_xx_yyzz_0[j] + fl1_fx * tly_x_yyzz_0[j] + 0.5 * fl1_fx * tpz_xx_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yyzz_0[j];

            tlz_xxx_yyzz_0[j] =
                pa_x[j] * tlz_xx_yyzz_0[j] + fl1_fx * tlz_x_yyzz_0[j] - 0.5 * fl1_fx * tpy_xx_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yyzz_0[j];

            tlx_xxx_yzzz_0[j] = pa_x[j] * tlx_xx_yzzz_0[j] + fl1_fx * tlx_x_yzzz_0[j];

            tly_xxx_yzzz_0[j] =
                pa_x[j] * tly_xx_yzzz_0[j] + fl1_fx * tly_x_yzzz_0[j] + 0.5 * fl1_fx * tpz_xx_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yzzz_0[j];

            tlz_xxx_yzzz_0[j] =
                pa_x[j] * tlz_xx_yzzz_0[j] + fl1_fx * tlz_x_yzzz_0[j] - 0.5 * fl1_fx * tpy_xx_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yzzz_0[j];

            tlx_xxx_zzzz_0[j] = pa_x[j] * tlx_xx_zzzz_0[j] + fl1_fx * tlx_x_zzzz_0[j];

            tly_xxx_zzzz_0[j] =
                pa_x[j] * tly_xx_zzzz_0[j] + fl1_fx * tly_x_zzzz_0[j] + 0.5 * fl1_fx * tpz_xx_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_zzzz_0[j];

            tlz_xxx_zzzz_0[j] =
                pa_x[j] * tlz_xx_zzzz_0[j] + fl1_fx * tlz_x_zzzz_0[j] - 0.5 * fl1_fx * tpy_xx_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_zzzz_0[j];

            tlx_xxy_xxxx_0[j] = pa_x[j] * tlx_xy_xxxx_0[j] + 0.5 * fl1_fx * tlx_y_xxxx_0[j] + 2.0 * fl1_fx * tlx_xy_xxx_0[j];

            tly_xxy_xxxx_0[j] = pa_x[j] * tly_xy_xxxx_0[j] + 0.5 * fl1_fx * tly_y_xxxx_0[j] + 2.0 * fl1_fx * tly_xy_xxx_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxxx_0[j];

            tlz_xxy_xxxx_0[j] = pa_x[j] * tlz_xy_xxxx_0[j] + 0.5 * fl1_fx * tlz_y_xxxx_0[j] + 2.0 * fl1_fx * tlz_xy_xxx_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxxx_0[j];

            tlx_xxy_xxxy_0[j] = pa_x[j] * tlx_xy_xxxy_0[j] + 0.5 * fl1_fx * tlx_y_xxxy_0[j] + 1.5 * fl1_fx * tlx_xy_xxy_0[j];

            tly_xxy_xxxy_0[j] = pa_x[j] * tly_xy_xxxy_0[j] + 0.5 * fl1_fx * tly_y_xxxy_0[j] + 1.5 * fl1_fx * tly_xy_xxy_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxxy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_50_100(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (50,100)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlz_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tlx_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 17);

        auto tly_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tlz_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tlx_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 18);

        auto tly_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tlz_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tlx_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 19);

        auto tly_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tlz_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tlx_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 20);

        auto tly_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tlz_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tlx_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 21);

        auto tly_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tlz_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tlx_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 22);

        auto tly_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tlz_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tlx_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 23);

        auto tly_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tlz_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tlx_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 24);

        auto tly_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tlz_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tlx_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 25);

        auto tly_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tlz_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tlx_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 26);

        auto tly_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tlz_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tlx_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 27);

        auto tly_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tlz_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tlx_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 28);

        auto tly_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tlz_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tlx_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 29);

        auto tly_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tlz_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tlx_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 30);

        auto tly_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tlz_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tlx_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 31);

        auto tly_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tlz_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tlx_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 32);

        auto tly_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tlz_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tlx_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 33);

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

        auto tlz_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tlx_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 12);

        auto tly_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tlz_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tlx_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 13);

        auto tly_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tlz_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tlx_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 14);

        auto tly_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tlz_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 14);

        auto tlx_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 15);

        auto tly_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tlz_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tlx_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 16);

        auto tly_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tlz_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 16);

        auto tlx_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 17);

        auto tly_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tlz_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tlx_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 18);

        auto tly_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tlz_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tlx_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 19);

        auto tly_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tlz_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tlx_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 20);

        auto tly_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tlz_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tlx_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 21);

        auto tly_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tlz_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tlx_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 22);

        auto tly_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tlz_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tlx_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 23);

        auto tpy_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tpy_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tpz_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tpy_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tpz_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tpy_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tpz_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tpy_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tpz_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tpy_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tpz_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tpy_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tpz_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tpy_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tpz_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tpy_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tpz_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tpy_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tpz_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tpy_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tpz_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tpy_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tpz_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tpy_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tpz_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tpy_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tpz_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tpy_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tpz_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tpy_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tpz_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tpy_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tpz_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tdy_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tdy_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tdz_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tdy_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tdz_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tdy_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tdz_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tdy_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tdz_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tdy_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tdz_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tdy_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tdz_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tdy_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tdz_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tdy_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tdz_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tdy_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tdz_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tdy_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tdz_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tdy_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tdz_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tdy_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tdz_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tdy_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tdz_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tdy_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tdz_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tdy_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tdz_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tdy_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tdz_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 32);

        // set up pointers to integrals

        auto tlz_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tlx_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 17);

        auto tly_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tlz_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tlx_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 18);

        auto tly_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tlz_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tlx_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 19);

        auto tly_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tlz_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tlx_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 20);

        auto tly_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tlz_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tlx_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 21);

        auto tly_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tlz_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tlx_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 22);

        auto tly_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tlz_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tlx_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 23);

        auto tly_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tlz_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tlx_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 24);

        auto tly_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tlz_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tlx_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 25);

        auto tly_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tlz_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tlx_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 26);

        auto tly_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tlz_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tlx_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 27);

        auto tly_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tlz_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tlx_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 28);

        auto tly_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tlz_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tlx_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 29);

        auto tly_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tlz_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tlx_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 30);

        auto tly_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tlz_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tlx_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 31);

        auto tly_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tlz_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tlx_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 32);

        auto tly_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tlz_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tlx_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xy_xxxy_0, tdy_xy_xxxz_0, tdy_xy_xxyy_0, \
                                     tdy_xy_xxyz_0, tdy_xy_xxzz_0, tdy_xy_xyyy_0, tdy_xy_xyyz_0, tdy_xy_xyzz_0, \
                                     tdy_xy_xzzz_0, tdy_xy_yyyy_0, tdy_xy_yyyz_0, tdy_xy_yyzz_0, tdy_xy_yzzz_0, \
                                     tdy_xy_zzzz_0, tdy_xz_xxxx_0, tdy_xz_xxxy_0, tdy_xz_xxxz_0, tdz_xy_xxxz_0, \
                                     tdz_xy_xxyy_0, tdz_xy_xxyz_0, tdz_xy_xxzz_0, tdz_xy_xyyy_0, tdz_xy_xyyz_0, \
                                     tdz_xy_xyzz_0, tdz_xy_xzzz_0, tdz_xy_yyyy_0, tdz_xy_yyyz_0, tdz_xy_yyzz_0, \
                                     tdz_xy_yzzz_0, tdz_xy_zzzz_0, tdz_xz_xxxx_0, tdz_xz_xxxy_0, tdz_xz_xxxz_0, \
                                     tlx_xxy_xxxz_0, tlx_xxy_xxyy_0, tlx_xxy_xxyz_0, tlx_xxy_xxzz_0, tlx_xxy_xyyy_0, \
                                     tlx_xxy_xyyz_0, tlx_xxy_xyzz_0, tlx_xxy_xzzz_0, tlx_xxy_yyyy_0, tlx_xxy_yyyz_0, \
                                     tlx_xxy_yyzz_0, tlx_xxy_yzzz_0, tlx_xxy_zzzz_0, tlx_xxz_xxxx_0, tlx_xxz_xxxy_0, \
                                     tlx_xxz_xxxz_0, tlx_xxz_xxyy_0, tlx_xy_xxxz_0, tlx_xy_xxyy_0, tlx_xy_xxyz_0, \
                                     tlx_xy_xxz_0, tlx_xy_xxzz_0, tlx_xy_xyy_0, tlx_xy_xyyy_0, tlx_xy_xyyz_0, \
                                     tlx_xy_xyz_0, tlx_xy_xyzz_0, tlx_xy_xzz_0, tlx_xy_xzzz_0, tlx_xy_yyy_0, \
                                     tlx_xy_yyyy_0, tlx_xy_yyyz_0, tlx_xy_yyz_0, tlx_xy_yyzz_0, tlx_xy_yzz_0, \
                                     tlx_xy_yzzz_0, tlx_xy_zzz_0, tlx_xy_zzzz_0, tlx_xz_xxx_0, tlx_xz_xxxx_0, \
                                     tlx_xz_xxxy_0, tlx_xz_xxxz_0, tlx_xz_xxy_0, tlx_xz_xxyy_0, tlx_xz_xxz_0, \
                                     tlx_xz_xyy_0, tlx_y_xxxz_0, tlx_y_xxyy_0, tlx_y_xxyz_0, tlx_y_xxzz_0, tlx_y_xyyy_0, \
                                     tlx_y_xyyz_0, tlx_y_xyzz_0, tlx_y_xzzz_0, tlx_y_yyyy_0, tlx_y_yyyz_0, tlx_y_yyzz_0, \
                                     tlx_y_yzzz_0, tlx_y_zzzz_0, tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, tlx_z_xxyy_0, \
                                     tly_xxy_xxxz_0, tly_xxy_xxyy_0, tly_xxy_xxyz_0, tly_xxy_xxzz_0, tly_xxy_xyyy_0, \
                                     tly_xxy_xyyz_0, tly_xxy_xyzz_0, tly_xxy_xzzz_0, tly_xxy_yyyy_0, tly_xxy_yyyz_0, \
                                     tly_xxy_yyzz_0, tly_xxy_yzzz_0, tly_xxy_zzzz_0, tly_xxz_xxxx_0, tly_xxz_xxxy_0, \
                                     tly_xxz_xxxz_0, tly_xy_xxxz_0, tly_xy_xxyy_0, tly_xy_xxyz_0, tly_xy_xxz_0, \
                                     tly_xy_xxzz_0, tly_xy_xyy_0, tly_xy_xyyy_0, tly_xy_xyyz_0, tly_xy_xyz_0, \
                                     tly_xy_xyzz_0, tly_xy_xzz_0, tly_xy_xzzz_0, tly_xy_yyy_0, tly_xy_yyyy_0, \
                                     tly_xy_yyyz_0, tly_xy_yyz_0, tly_xy_yyzz_0, tly_xy_yzz_0, tly_xy_yzzz_0, \
                                     tly_xy_zzz_0, tly_xy_zzzz_0, tly_xz_xxx_0, tly_xz_xxxx_0, tly_xz_xxxy_0, \
                                     tly_xz_xxxz_0, tly_xz_xxy_0, tly_xz_xxz_0, tly_y_xxxz_0, tly_y_xxyy_0, tly_y_xxyz_0, \
                                     tly_y_xxzz_0, tly_y_xyyy_0, tly_y_xyyz_0, tly_y_xyzz_0, tly_y_xzzz_0, tly_y_yyyy_0, \
                                     tly_y_yyyz_0, tly_y_yyzz_0, tly_y_yzzz_0, tly_y_zzzz_0, tly_z_xxxx_0, tly_z_xxxy_0, \
                                     tly_z_xxxz_0, tlz_xxy_xxxy_0, tlz_xxy_xxxz_0, tlz_xxy_xxyy_0, tlz_xxy_xxyz_0, \
                                     tlz_xxy_xxzz_0, tlz_xxy_xyyy_0, tlz_xxy_xyyz_0, tlz_xxy_xyzz_0, tlz_xxy_xzzz_0, \
                                     tlz_xxy_yyyy_0, tlz_xxy_yyyz_0, tlz_xxy_yyzz_0, tlz_xxy_yzzz_0, tlz_xxy_zzzz_0, \
                                     tlz_xxz_xxxx_0, tlz_xxz_xxxy_0, tlz_xxz_xxxz_0, tlz_xy_xxxy_0, tlz_xy_xxxz_0, \
                                     tlz_xy_xxy_0, tlz_xy_xxyy_0, tlz_xy_xxyz_0, tlz_xy_xxz_0, tlz_xy_xxzz_0, \
                                     tlz_xy_xyy_0, tlz_xy_xyyy_0, tlz_xy_xyyz_0, tlz_xy_xyz_0, tlz_xy_xyzz_0, \
                                     tlz_xy_xzz_0, tlz_xy_xzzz_0, tlz_xy_yyy_0, tlz_xy_yyyy_0, tlz_xy_yyyz_0, \
                                     tlz_xy_yyz_0, tlz_xy_yyzz_0, tlz_xy_yzz_0, tlz_xy_yzzz_0, tlz_xy_zzz_0, \
                                     tlz_xy_zzzz_0, tlz_xz_xxx_0, tlz_xz_xxxx_0, tlz_xz_xxxy_0, tlz_xz_xxxz_0, \
                                     tlz_xz_xxy_0, tlz_xz_xxz_0, tlz_y_xxxy_0, tlz_y_xxxz_0, tlz_y_xxyy_0, tlz_y_xxyz_0, \
                                     tlz_y_xxzz_0, tlz_y_xyyy_0, tlz_y_xyyz_0, tlz_y_xyzz_0, tlz_y_xzzz_0, tlz_y_yyyy_0, \
                                     tlz_y_yyyz_0, tlz_y_yyzz_0, tlz_y_yzzz_0, tlz_y_zzzz_0, tlz_z_xxxx_0, tlz_z_xxxy_0, \
                                     tlz_z_xxxz_0, tpy_xy_xxxy_0, tpy_xy_xxxz_0, tpy_xy_xxyy_0, tpy_xy_xxyz_0, \
                                     tpy_xy_xxzz_0, tpy_xy_xyyy_0, tpy_xy_xyyz_0, tpy_xy_xyzz_0, tpy_xy_xzzz_0, \
                                     tpy_xy_yyyy_0, tpy_xy_yyyz_0, tpy_xy_yyzz_0, tpy_xy_yzzz_0, tpy_xy_zzzz_0, \
                                     tpy_xz_xxxx_0, tpy_xz_xxxy_0, tpy_xz_xxxz_0, tpz_xy_xxxz_0, tpz_xy_xxyy_0, \
                                     tpz_xy_xxyz_0, tpz_xy_xxzz_0, tpz_xy_xyyy_0, tpz_xy_xyyz_0, tpz_xy_xyzz_0, \
                                     tpz_xy_xzzz_0, tpz_xy_yyyy_0, tpz_xy_yyyz_0, tpz_xy_yyzz_0, tpz_xy_yzzz_0, \
                                     tpz_xy_zzzz_0, tpz_xz_xxxx_0, tpz_xz_xxxy_0, tpz_xz_xxxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xxy_xxxy_0[j] = pa_x[j] * tlz_xy_xxxy_0[j] + 0.5 * fl1_fx * tlz_y_xxxy_0[j] + 1.5 * fl1_fx * tlz_xy_xxy_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxxy_0[j];

            tlx_xxy_xxxz_0[j] = pa_x[j] * tlx_xy_xxxz_0[j] + 0.5 * fl1_fx * tlx_y_xxxz_0[j] + 1.5 * fl1_fx * tlx_xy_xxz_0[j];

            tly_xxy_xxxz_0[j] = pa_x[j] * tly_xy_xxxz_0[j] + 0.5 * fl1_fx * tly_y_xxxz_0[j] + 1.5 * fl1_fx * tly_xy_xxz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxxz_0[j];

            tlz_xxy_xxxz_0[j] = pa_x[j] * tlz_xy_xxxz_0[j] + 0.5 * fl1_fx * tlz_y_xxxz_0[j] + 1.5 * fl1_fx * tlz_xy_xxz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxxz_0[j];

            tlx_xxy_xxyy_0[j] = pa_x[j] * tlx_xy_xxyy_0[j] + 0.5 * fl1_fx * tlx_y_xxyy_0[j] + fl1_fx * tlx_xy_xyy_0[j];

            tly_xxy_xxyy_0[j] = pa_x[j] * tly_xy_xxyy_0[j] + 0.5 * fl1_fx * tly_y_xxyy_0[j] + fl1_fx * tly_xy_xyy_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxyy_0[j];

            tlz_xxy_xxyy_0[j] = pa_x[j] * tlz_xy_xxyy_0[j] + 0.5 * fl1_fx * tlz_y_xxyy_0[j] + fl1_fx * tlz_xy_xyy_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxyy_0[j];

            tlx_xxy_xxyz_0[j] = pa_x[j] * tlx_xy_xxyz_0[j] + 0.5 * fl1_fx * tlx_y_xxyz_0[j] + fl1_fx * tlx_xy_xyz_0[j];

            tly_xxy_xxyz_0[j] = pa_x[j] * tly_xy_xxyz_0[j] + 0.5 * fl1_fx * tly_y_xxyz_0[j] + fl1_fx * tly_xy_xyz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxyz_0[j];

            tlz_xxy_xxyz_0[j] = pa_x[j] * tlz_xy_xxyz_0[j] + 0.5 * fl1_fx * tlz_y_xxyz_0[j] + fl1_fx * tlz_xy_xyz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxyz_0[j];

            tlx_xxy_xxzz_0[j] = pa_x[j] * tlx_xy_xxzz_0[j] + 0.5 * fl1_fx * tlx_y_xxzz_0[j] + fl1_fx * tlx_xy_xzz_0[j];

            tly_xxy_xxzz_0[j] = pa_x[j] * tly_xy_xxzz_0[j] + 0.5 * fl1_fx * tly_y_xxzz_0[j] + fl1_fx * tly_xy_xzz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxzz_0[j];

            tlz_xxy_xxzz_0[j] = pa_x[j] * tlz_xy_xxzz_0[j] + 0.5 * fl1_fx * tlz_y_xxzz_0[j] + fl1_fx * tlz_xy_xzz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxzz_0[j];

            tlx_xxy_xyyy_0[j] = pa_x[j] * tlx_xy_xyyy_0[j] + 0.5 * fl1_fx * tlx_y_xyyy_0[j] + 0.5 * fl1_fx * tlx_xy_yyy_0[j];

            tly_xxy_xyyy_0[j] = pa_x[j] * tly_xy_xyyy_0[j] + 0.5 * fl1_fx * tly_y_xyyy_0[j] + 0.5 * fl1_fx * tly_xy_yyy_0[j] +
                                0.5 * fl1_fx * tpz_xy_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xy_xyyy_0[j];

            tlz_xxy_xyyy_0[j] = pa_x[j] * tlz_xy_xyyy_0[j] + 0.5 * fl1_fx * tlz_y_xyyy_0[j] + 0.5 * fl1_fx * tlz_xy_yyy_0[j] -
                                0.5 * fl1_fx * tpy_xy_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xy_xyyy_0[j];

            tlx_xxy_xyyz_0[j] = pa_x[j] * tlx_xy_xyyz_0[j] + 0.5 * fl1_fx * tlx_y_xyyz_0[j] + 0.5 * fl1_fx * tlx_xy_yyz_0[j];

            tly_xxy_xyyz_0[j] = pa_x[j] * tly_xy_xyyz_0[j] + 0.5 * fl1_fx * tly_y_xyyz_0[j] + 0.5 * fl1_fx * tly_xy_yyz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xyyz_0[j];

            tlz_xxy_xyyz_0[j] = pa_x[j] * tlz_xy_xyyz_0[j] + 0.5 * fl1_fx * tlz_y_xyyz_0[j] + 0.5 * fl1_fx * tlz_xy_yyz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xyyz_0[j];

            tlx_xxy_xyzz_0[j] = pa_x[j] * tlx_xy_xyzz_0[j] + 0.5 * fl1_fx * tlx_y_xyzz_0[j] + 0.5 * fl1_fx * tlx_xy_yzz_0[j];

            tly_xxy_xyzz_0[j] = pa_x[j] * tly_xy_xyzz_0[j] + 0.5 * fl1_fx * tly_y_xyzz_0[j] + 0.5 * fl1_fx * tly_xy_yzz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xyzz_0[j];

            tlz_xxy_xyzz_0[j] = pa_x[j] * tlz_xy_xyzz_0[j] + 0.5 * fl1_fx * tlz_y_xyzz_0[j] + 0.5 * fl1_fx * tlz_xy_yzz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xyzz_0[j];

            tlx_xxy_xzzz_0[j] = pa_x[j] * tlx_xy_xzzz_0[j] + 0.5 * fl1_fx * tlx_y_xzzz_0[j] + 0.5 * fl1_fx * tlx_xy_zzz_0[j];

            tly_xxy_xzzz_0[j] = pa_x[j] * tly_xy_xzzz_0[j] + 0.5 * fl1_fx * tly_y_xzzz_0[j] + 0.5 * fl1_fx * tly_xy_zzz_0[j] +
                                0.5 * fl1_fx * tpz_xy_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xzzz_0[j];

            tlz_xxy_xzzz_0[j] = pa_x[j] * tlz_xy_xzzz_0[j] + 0.5 * fl1_fx * tlz_y_xzzz_0[j] + 0.5 * fl1_fx * tlz_xy_zzz_0[j] -
                                0.5 * fl1_fx * tpy_xy_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xzzz_0[j];

            tlx_xxy_yyyy_0[j] = pa_x[j] * tlx_xy_yyyy_0[j] + 0.5 * fl1_fx * tlx_y_yyyy_0[j];

            tly_xxy_yyyy_0[j] =
                pa_x[j] * tly_xy_yyyy_0[j] + 0.5 * fl1_fx * tly_y_yyyy_0[j] + 0.5 * fl1_fx * tpz_xy_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_xy_yyyy_0[j];

            tlz_xxy_yyyy_0[j] =
                pa_x[j] * tlz_xy_yyyy_0[j] + 0.5 * fl1_fx * tlz_y_yyyy_0[j] - 0.5 * fl1_fx * tpy_xy_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_xy_yyyy_0[j];

            tlx_xxy_yyyz_0[j] = pa_x[j] * tlx_xy_yyyz_0[j] + 0.5 * fl1_fx * tlx_y_yyyz_0[j];

            tly_xxy_yyyz_0[j] =
                pa_x[j] * tly_xy_yyyz_0[j] + 0.5 * fl1_fx * tly_y_yyyz_0[j] + 0.5 * fl1_fx * tpz_xy_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yyyz_0[j];

            tlz_xxy_yyyz_0[j] =
                pa_x[j] * tlz_xy_yyyz_0[j] + 0.5 * fl1_fx * tlz_y_yyyz_0[j] - 0.5 * fl1_fx * tpy_xy_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yyyz_0[j];

            tlx_xxy_yyzz_0[j] = pa_x[j] * tlx_xy_yyzz_0[j] + 0.5 * fl1_fx * tlx_y_yyzz_0[j];

            tly_xxy_yyzz_0[j] =
                pa_x[j] * tly_xy_yyzz_0[j] + 0.5 * fl1_fx * tly_y_yyzz_0[j] + 0.5 * fl1_fx * tpz_xy_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yyzz_0[j];

            tlz_xxy_yyzz_0[j] =
                pa_x[j] * tlz_xy_yyzz_0[j] + 0.5 * fl1_fx * tlz_y_yyzz_0[j] - 0.5 * fl1_fx * tpy_xy_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yyzz_0[j];

            tlx_xxy_yzzz_0[j] = pa_x[j] * tlx_xy_yzzz_0[j] + 0.5 * fl1_fx * tlx_y_yzzz_0[j];

            tly_xxy_yzzz_0[j] =
                pa_x[j] * tly_xy_yzzz_0[j] + 0.5 * fl1_fx * tly_y_yzzz_0[j] + 0.5 * fl1_fx * tpz_xy_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yzzz_0[j];

            tlz_xxy_yzzz_0[j] =
                pa_x[j] * tlz_xy_yzzz_0[j] + 0.5 * fl1_fx * tlz_y_yzzz_0[j] - 0.5 * fl1_fx * tpy_xy_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yzzz_0[j];

            tlx_xxy_zzzz_0[j] = pa_x[j] * tlx_xy_zzzz_0[j] + 0.5 * fl1_fx * tlx_y_zzzz_0[j];

            tly_xxy_zzzz_0[j] =
                pa_x[j] * tly_xy_zzzz_0[j] + 0.5 * fl1_fx * tly_y_zzzz_0[j] + 0.5 * fl1_fx * tpz_xy_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_zzzz_0[j];

            tlz_xxy_zzzz_0[j] =
                pa_x[j] * tlz_xy_zzzz_0[j] + 0.5 * fl1_fx * tlz_y_zzzz_0[j] - 0.5 * fl1_fx * tpy_xy_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_zzzz_0[j];

            tlx_xxz_xxxx_0[j] = pa_x[j] * tlx_xz_xxxx_0[j] + 0.5 * fl1_fx * tlx_z_xxxx_0[j] + 2.0 * fl1_fx * tlx_xz_xxx_0[j];

            tly_xxz_xxxx_0[j] = pa_x[j] * tly_xz_xxxx_0[j] + 0.5 * fl1_fx * tly_z_xxxx_0[j] + 2.0 * fl1_fx * tly_xz_xxx_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxxx_0[j];

            tlz_xxz_xxxx_0[j] = pa_x[j] * tlz_xz_xxxx_0[j] + 0.5 * fl1_fx * tlz_z_xxxx_0[j] + 2.0 * fl1_fx * tlz_xz_xxx_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxxx_0[j];

            tlx_xxz_xxxy_0[j] = pa_x[j] * tlx_xz_xxxy_0[j] + 0.5 * fl1_fx * tlx_z_xxxy_0[j] + 1.5 * fl1_fx * tlx_xz_xxy_0[j];

            tly_xxz_xxxy_0[j] = pa_x[j] * tly_xz_xxxy_0[j] + 0.5 * fl1_fx * tly_z_xxxy_0[j] + 1.5 * fl1_fx * tly_xz_xxy_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxxy_0[j];

            tlz_xxz_xxxy_0[j] = pa_x[j] * tlz_xz_xxxy_0[j] + 0.5 * fl1_fx * tlz_z_xxxy_0[j] + 1.5 * fl1_fx * tlz_xz_xxy_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxxy_0[j];

            tlx_xxz_xxxz_0[j] = pa_x[j] * tlx_xz_xxxz_0[j] + 0.5 * fl1_fx * tlx_z_xxxz_0[j] + 1.5 * fl1_fx * tlx_xz_xxz_0[j];

            tly_xxz_xxxz_0[j] = pa_x[j] * tly_xz_xxxz_0[j] + 0.5 * fl1_fx * tly_z_xxxz_0[j] + 1.5 * fl1_fx * tly_xz_xxz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxxz_0[j];

            tlz_xxz_xxxz_0[j] = pa_x[j] * tlz_xz_xxxz_0[j] + 0.5 * fl1_fx * tlz_z_xxxz_0[j] + 1.5 * fl1_fx * tlz_xz_xxz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxxz_0[j];

            tlx_xxz_xxyy_0[j] = pa_x[j] * tlx_xz_xxyy_0[j] + 0.5 * fl1_fx * tlx_z_xxyy_0[j] + fl1_fx * tlx_xz_xyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_100_150(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (100,150)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tly_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tlz_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tlx_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 34);

        auto tly_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tlz_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tlx_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 35);

        auto tly_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tlz_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tlx_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 36);

        auto tly_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tlz_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tlx_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 37);

        auto tly_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tlz_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tlx_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 38);

        auto tly_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tlz_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tlx_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 39);

        auto tly_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tlz_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tlx_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 40);

        auto tly_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tlz_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tlx_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 41);

        auto tly_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tlz_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tlx_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 42);

        auto tly_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tlz_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tlx_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 43);

        auto tly_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tlz_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tlx_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 44);

        auto tly_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tlz_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tlx_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 45);

        auto tly_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tlz_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tlx_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 46);

        auto tly_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tlz_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tlx_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 47);

        auto tly_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tlz_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tlx_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 48);

        auto tly_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tlz_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tlx_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 49);

        auto tly_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tlz_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 49);

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

        auto tly_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tlz_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tlx_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 24);

        auto tly_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tlz_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tlx_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 25);

        auto tly_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tlz_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tlx_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 26);

        auto tly_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tlz_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tlx_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 27);

        auto tly_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tlz_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tlx_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 28);

        auto tly_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tlz_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tlx_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 29);

        auto tly_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tlz_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 29);

        auto tlx_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 30);

        auto tly_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tlz_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tlx_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 31);

        auto tly_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tlz_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tlx_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 32);

        auto tly_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tlz_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tlx_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 33);

        auto tly_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tlz_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tlx_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 34);

        auto tly_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tlz_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tpy_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tpz_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tpy_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tpz_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tpy_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tpz_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tpy_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tpz_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tpy_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tpz_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tpy_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tpz_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tpy_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tpz_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tpy_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tpz_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tpy_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tpz_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tpy_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tpz_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tpy_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tpz_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tpy_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tpz_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tpy_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tpz_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tpy_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tpz_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tpy_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tpz_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tpy_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tpz_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tpy_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tdy_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tdz_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tdy_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tdz_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tdy_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tdz_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tdy_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tdz_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tdy_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tdz_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tdy_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tdz_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tdy_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tdz_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tdy_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tdz_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tdy_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tdz_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tdy_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tdz_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tdy_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tdz_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tdy_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tdz_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tdy_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tdz_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tdy_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tdz_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tdy_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tdz_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tdy_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tdz_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tdy_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tdz_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 49);

        // set up pointers to integrals

        auto tly_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tlz_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tlx_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 34);

        auto tly_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tlz_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tlx_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 35);

        auto tly_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tlz_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tlx_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 36);

        auto tly_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tlz_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tlx_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 37);

        auto tly_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tlz_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tlx_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 38);

        auto tly_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tlz_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tlx_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 39);

        auto tly_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tlz_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tlx_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 40);

        auto tly_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tlz_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tlx_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 41);

        auto tly_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tlz_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tlx_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 42);

        auto tly_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tlz_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tlx_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 43);

        auto tly_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tlz_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tlx_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 44);

        auto tly_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tlz_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tlx_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 45);

        auto tly_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tlz_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tlx_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 46);

        auto tly_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tlz_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tlx_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 47);

        auto tly_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tlz_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tlx_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 48);

        auto tly_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tlz_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 48);

        auto tlx_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 49);

        auto tly_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tlz_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xz_xxyy_0, tdy_xz_xxyz_0, tdy_xz_xxzz_0, \
                                     tdy_xz_xyyy_0, tdy_xz_xyyz_0, tdy_xz_xyzz_0, tdy_xz_xzzz_0, tdy_xz_yyyy_0, \
                                     tdy_xz_yyyz_0, tdy_xz_yyzz_0, tdy_xz_yzzz_0, tdy_xz_zzzz_0, tdy_yy_xxxx_0, \
                                     tdy_yy_xxxy_0, tdy_yy_xxxz_0, tdy_yy_xxyy_0, tdy_yy_xxyz_0, tdz_xz_xxyy_0, \
                                     tdz_xz_xxyz_0, tdz_xz_xxzz_0, tdz_xz_xyyy_0, tdz_xz_xyyz_0, tdz_xz_xyzz_0, \
                                     tdz_xz_xzzz_0, tdz_xz_yyyy_0, tdz_xz_yyyz_0, tdz_xz_yyzz_0, tdz_xz_yzzz_0, \
                                     tdz_xz_zzzz_0, tdz_yy_xxxx_0, tdz_yy_xxxy_0, tdz_yy_xxxz_0, tdz_yy_xxyy_0, \
                                     tdz_yy_xxyz_0, tlx_xxz_xxyz_0, tlx_xxz_xxzz_0, tlx_xxz_xyyy_0, tlx_xxz_xyyz_0, \
                                     tlx_xxz_xyzz_0, tlx_xxz_xzzz_0, tlx_xxz_yyyy_0, tlx_xxz_yyyz_0, tlx_xxz_yyzz_0, \
                                     tlx_xxz_yzzz_0, tlx_xxz_zzzz_0, tlx_xyy_xxxx_0, tlx_xyy_xxxy_0, tlx_xyy_xxxz_0, \
                                     tlx_xyy_xxyy_0, tlx_xyy_xxyz_0, tlx_xz_xxyz_0, tlx_xz_xxzz_0, tlx_xz_xyyy_0, \
                                     tlx_xz_xyyz_0, tlx_xz_xyz_0, tlx_xz_xyzz_0, tlx_xz_xzz_0, tlx_xz_xzzz_0, \
                                     tlx_xz_yyy_0, tlx_xz_yyyy_0, tlx_xz_yyyz_0, tlx_xz_yyz_0, tlx_xz_yyzz_0, \
                                     tlx_xz_yzz_0, tlx_xz_yzzz_0, tlx_xz_zzz_0, tlx_xz_zzzz_0, tlx_yy_xxx_0, \
                                     tlx_yy_xxxx_0, tlx_yy_xxxy_0, tlx_yy_xxxz_0, tlx_yy_xxy_0, tlx_yy_xxyy_0, \
                                     tlx_yy_xxyz_0, tlx_yy_xxz_0, tlx_yy_xyy_0, tlx_yy_xyz_0, tlx_z_xxyz_0, tlx_z_xxzz_0, \
                                     tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyzz_0, tlx_z_xzzz_0, tlx_z_yyyy_0, tlx_z_yyyz_0, \
                                     tlx_z_yyzz_0, tlx_z_yzzz_0, tlx_z_zzzz_0, tly_xxz_xxyy_0, tly_xxz_xxyz_0, \
                                     tly_xxz_xxzz_0, tly_xxz_xyyy_0, tly_xxz_xyyz_0, tly_xxz_xyzz_0, tly_xxz_xzzz_0, \
                                     tly_xxz_yyyy_0, tly_xxz_yyyz_0, tly_xxz_yyzz_0, tly_xxz_yzzz_0, tly_xxz_zzzz_0, \
                                     tly_xyy_xxxx_0, tly_xyy_xxxy_0, tly_xyy_xxxz_0, tly_xyy_xxyy_0, tly_xyy_xxyz_0, \
                                     tly_xz_xxyy_0, tly_xz_xxyz_0, tly_xz_xxzz_0, tly_xz_xyy_0, tly_xz_xyyy_0, \
                                     tly_xz_xyyz_0, tly_xz_xyz_0, tly_xz_xyzz_0, tly_xz_xzz_0, tly_xz_xzzz_0, \
                                     tly_xz_yyy_0, tly_xz_yyyy_0, tly_xz_yyyz_0, tly_xz_yyz_0, tly_xz_yyzz_0, \
                                     tly_xz_yzz_0, tly_xz_yzzz_0, tly_xz_zzz_0, tly_xz_zzzz_0, tly_yy_xxx_0, \
                                     tly_yy_xxxx_0, tly_yy_xxxy_0, tly_yy_xxxz_0, tly_yy_xxy_0, tly_yy_xxyy_0, \
                                     tly_yy_xxyz_0, tly_yy_xxz_0, tly_yy_xyy_0, tly_yy_xyz_0, tly_z_xxyy_0, tly_z_xxyz_0, \
                                     tly_z_xxzz_0, tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyzz_0, tly_z_xzzz_0, tly_z_yyyy_0, \
                                     tly_z_yyyz_0, tly_z_yyzz_0, tly_z_yzzz_0, tly_z_zzzz_0, tlz_xxz_xxyy_0, \
                                     tlz_xxz_xxyz_0, tlz_xxz_xxzz_0, tlz_xxz_xyyy_0, tlz_xxz_xyyz_0, tlz_xxz_xyzz_0, \
                                     tlz_xxz_xzzz_0, tlz_xxz_yyyy_0, tlz_xxz_yyyz_0, tlz_xxz_yyzz_0, tlz_xxz_yzzz_0, \
                                     tlz_xxz_zzzz_0, tlz_xyy_xxxx_0, tlz_xyy_xxxy_0, tlz_xyy_xxxz_0, tlz_xyy_xxyy_0, \
                                     tlz_xyy_xxyz_0, tlz_xz_xxyy_0, tlz_xz_xxyz_0, tlz_xz_xxzz_0, tlz_xz_xyy_0, \
                                     tlz_xz_xyyy_0, tlz_xz_xyyz_0, tlz_xz_xyz_0, tlz_xz_xyzz_0, tlz_xz_xzz_0, \
                                     tlz_xz_xzzz_0, tlz_xz_yyy_0, tlz_xz_yyyy_0, tlz_xz_yyyz_0, tlz_xz_yyz_0, \
                                     tlz_xz_yyzz_0, tlz_xz_yzz_0, tlz_xz_yzzz_0, tlz_xz_zzz_0, tlz_xz_zzzz_0, \
                                     tlz_yy_xxx_0, tlz_yy_xxxx_0, tlz_yy_xxxy_0, tlz_yy_xxxz_0, tlz_yy_xxy_0, \
                                     tlz_yy_xxyy_0, tlz_yy_xxyz_0, tlz_yy_xxz_0, tlz_yy_xyy_0, tlz_yy_xyz_0, \
                                     tlz_z_xxyy_0, tlz_z_xxyz_0, tlz_z_xxzz_0, tlz_z_xyyy_0, tlz_z_xyyz_0, tlz_z_xyzz_0, \
                                     tlz_z_xzzz_0, tlz_z_yyyy_0, tlz_z_yyyz_0, tlz_z_yyzz_0, tlz_z_yzzz_0, tlz_z_zzzz_0, \
                                     tpy_xz_xxyy_0, tpy_xz_xxyz_0, tpy_xz_xxzz_0, tpy_xz_xyyy_0, tpy_xz_xyyz_0, \
                                     tpy_xz_xyzz_0, tpy_xz_xzzz_0, tpy_xz_yyyy_0, tpy_xz_yyyz_0, tpy_xz_yyzz_0, \
                                     tpy_xz_yzzz_0, tpy_xz_zzzz_0, tpy_yy_xxxx_0, tpy_yy_xxxy_0, tpy_yy_xxxz_0, \
                                     tpy_yy_xxyy_0, tpy_yy_xxyz_0, tpz_xz_xxyy_0, tpz_xz_xxyz_0, tpz_xz_xxzz_0, \
                                     tpz_xz_xyyy_0, tpz_xz_xyyz_0, tpz_xz_xyzz_0, tpz_xz_xzzz_0, tpz_xz_yyyy_0, \
                                     tpz_xz_yyyz_0, tpz_xz_yyzz_0, tpz_xz_yzzz_0, tpz_xz_zzzz_0, tpz_yy_xxxx_0, \
                                     tpz_yy_xxxy_0, tpz_yy_xxxz_0, tpz_yy_xxyy_0, tpz_yy_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xxz_xxyy_0[j] = pa_x[j] * tly_xz_xxyy_0[j] + 0.5 * fl1_fx * tly_z_xxyy_0[j] + fl1_fx * tly_xz_xyy_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxyy_0[j];

            tlz_xxz_xxyy_0[j] = pa_x[j] * tlz_xz_xxyy_0[j] + 0.5 * fl1_fx * tlz_z_xxyy_0[j] + fl1_fx * tlz_xz_xyy_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxyy_0[j];

            tlx_xxz_xxyz_0[j] = pa_x[j] * tlx_xz_xxyz_0[j] + 0.5 * fl1_fx * tlx_z_xxyz_0[j] + fl1_fx * tlx_xz_xyz_0[j];

            tly_xxz_xxyz_0[j] = pa_x[j] * tly_xz_xxyz_0[j] + 0.5 * fl1_fx * tly_z_xxyz_0[j] + fl1_fx * tly_xz_xyz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxyz_0[j];

            tlz_xxz_xxyz_0[j] = pa_x[j] * tlz_xz_xxyz_0[j] + 0.5 * fl1_fx * tlz_z_xxyz_0[j] + fl1_fx * tlz_xz_xyz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxyz_0[j];

            tlx_xxz_xxzz_0[j] = pa_x[j] * tlx_xz_xxzz_0[j] + 0.5 * fl1_fx * tlx_z_xxzz_0[j] + fl1_fx * tlx_xz_xzz_0[j];

            tly_xxz_xxzz_0[j] = pa_x[j] * tly_xz_xxzz_0[j] + 0.5 * fl1_fx * tly_z_xxzz_0[j] + fl1_fx * tly_xz_xzz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxzz_0[j];

            tlz_xxz_xxzz_0[j] = pa_x[j] * tlz_xz_xxzz_0[j] + 0.5 * fl1_fx * tlz_z_xxzz_0[j] + fl1_fx * tlz_xz_xzz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxzz_0[j];

            tlx_xxz_xyyy_0[j] = pa_x[j] * tlx_xz_xyyy_0[j] + 0.5 * fl1_fx * tlx_z_xyyy_0[j] + 0.5 * fl1_fx * tlx_xz_yyy_0[j];

            tly_xxz_xyyy_0[j] = pa_x[j] * tly_xz_xyyy_0[j] + 0.5 * fl1_fx * tly_z_xyyy_0[j] + 0.5 * fl1_fx * tly_xz_yyy_0[j] +
                                0.5 * fl1_fx * tpz_xz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xz_xyyy_0[j];

            tlz_xxz_xyyy_0[j] = pa_x[j] * tlz_xz_xyyy_0[j] + 0.5 * fl1_fx * tlz_z_xyyy_0[j] + 0.5 * fl1_fx * tlz_xz_yyy_0[j] -
                                0.5 * fl1_fx * tpy_xz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xz_xyyy_0[j];

            tlx_xxz_xyyz_0[j] = pa_x[j] * tlx_xz_xyyz_0[j] + 0.5 * fl1_fx * tlx_z_xyyz_0[j] + 0.5 * fl1_fx * tlx_xz_yyz_0[j];

            tly_xxz_xyyz_0[j] = pa_x[j] * tly_xz_xyyz_0[j] + 0.5 * fl1_fx * tly_z_xyyz_0[j] + 0.5 * fl1_fx * tly_xz_yyz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xyyz_0[j];

            tlz_xxz_xyyz_0[j] = pa_x[j] * tlz_xz_xyyz_0[j] + 0.5 * fl1_fx * tlz_z_xyyz_0[j] + 0.5 * fl1_fx * tlz_xz_yyz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xyyz_0[j];

            tlx_xxz_xyzz_0[j] = pa_x[j] * tlx_xz_xyzz_0[j] + 0.5 * fl1_fx * tlx_z_xyzz_0[j] + 0.5 * fl1_fx * tlx_xz_yzz_0[j];

            tly_xxz_xyzz_0[j] = pa_x[j] * tly_xz_xyzz_0[j] + 0.5 * fl1_fx * tly_z_xyzz_0[j] + 0.5 * fl1_fx * tly_xz_yzz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xyzz_0[j];

            tlz_xxz_xyzz_0[j] = pa_x[j] * tlz_xz_xyzz_0[j] + 0.5 * fl1_fx * tlz_z_xyzz_0[j] + 0.5 * fl1_fx * tlz_xz_yzz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xyzz_0[j];

            tlx_xxz_xzzz_0[j] = pa_x[j] * tlx_xz_xzzz_0[j] + 0.5 * fl1_fx * tlx_z_xzzz_0[j] + 0.5 * fl1_fx * tlx_xz_zzz_0[j];

            tly_xxz_xzzz_0[j] = pa_x[j] * tly_xz_xzzz_0[j] + 0.5 * fl1_fx * tly_z_xzzz_0[j] + 0.5 * fl1_fx * tly_xz_zzz_0[j] +
                                0.5 * fl1_fx * tpz_xz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xzzz_0[j];

            tlz_xxz_xzzz_0[j] = pa_x[j] * tlz_xz_xzzz_0[j] + 0.5 * fl1_fx * tlz_z_xzzz_0[j] + 0.5 * fl1_fx * tlz_xz_zzz_0[j] -
                                0.5 * fl1_fx * tpy_xz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xzzz_0[j];

            tlx_xxz_yyyy_0[j] = pa_x[j] * tlx_xz_yyyy_0[j] + 0.5 * fl1_fx * tlx_z_yyyy_0[j];

            tly_xxz_yyyy_0[j] =
                pa_x[j] * tly_xz_yyyy_0[j] + 0.5 * fl1_fx * tly_z_yyyy_0[j] + 0.5 * fl1_fx * tpz_xz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_xz_yyyy_0[j];

            tlz_xxz_yyyy_0[j] =
                pa_x[j] * tlz_xz_yyyy_0[j] + 0.5 * fl1_fx * tlz_z_yyyy_0[j] - 0.5 * fl1_fx * tpy_xz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_xz_yyyy_0[j];

            tlx_xxz_yyyz_0[j] = pa_x[j] * tlx_xz_yyyz_0[j] + 0.5 * fl1_fx * tlx_z_yyyz_0[j];

            tly_xxz_yyyz_0[j] =
                pa_x[j] * tly_xz_yyyz_0[j] + 0.5 * fl1_fx * tly_z_yyyz_0[j] + 0.5 * fl1_fx * tpz_xz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yyyz_0[j];

            tlz_xxz_yyyz_0[j] =
                pa_x[j] * tlz_xz_yyyz_0[j] + 0.5 * fl1_fx * tlz_z_yyyz_0[j] - 0.5 * fl1_fx * tpy_xz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yyyz_0[j];

            tlx_xxz_yyzz_0[j] = pa_x[j] * tlx_xz_yyzz_0[j] + 0.5 * fl1_fx * tlx_z_yyzz_0[j];

            tly_xxz_yyzz_0[j] =
                pa_x[j] * tly_xz_yyzz_0[j] + 0.5 * fl1_fx * tly_z_yyzz_0[j] + 0.5 * fl1_fx * tpz_xz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yyzz_0[j];

            tlz_xxz_yyzz_0[j] =
                pa_x[j] * tlz_xz_yyzz_0[j] + 0.5 * fl1_fx * tlz_z_yyzz_0[j] - 0.5 * fl1_fx * tpy_xz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yyzz_0[j];

            tlx_xxz_yzzz_0[j] = pa_x[j] * tlx_xz_yzzz_0[j] + 0.5 * fl1_fx * tlx_z_yzzz_0[j];

            tly_xxz_yzzz_0[j] =
                pa_x[j] * tly_xz_yzzz_0[j] + 0.5 * fl1_fx * tly_z_yzzz_0[j] + 0.5 * fl1_fx * tpz_xz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yzzz_0[j];

            tlz_xxz_yzzz_0[j] =
                pa_x[j] * tlz_xz_yzzz_0[j] + 0.5 * fl1_fx * tlz_z_yzzz_0[j] - 0.5 * fl1_fx * tpy_xz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yzzz_0[j];

            tlx_xxz_zzzz_0[j] = pa_x[j] * tlx_xz_zzzz_0[j] + 0.5 * fl1_fx * tlx_z_zzzz_0[j];

            tly_xxz_zzzz_0[j] =
                pa_x[j] * tly_xz_zzzz_0[j] + 0.5 * fl1_fx * tly_z_zzzz_0[j] + 0.5 * fl1_fx * tpz_xz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_zzzz_0[j];

            tlz_xxz_zzzz_0[j] =
                pa_x[j] * tlz_xz_zzzz_0[j] + 0.5 * fl1_fx * tlz_z_zzzz_0[j] - 0.5 * fl1_fx * tpy_xz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_zzzz_0[j];

            tlx_xyy_xxxx_0[j] = pa_x[j] * tlx_yy_xxxx_0[j] + 2.0 * fl1_fx * tlx_yy_xxx_0[j];

            tly_xyy_xxxx_0[j] =
                pa_x[j] * tly_yy_xxxx_0[j] + 2.0 * fl1_fx * tly_yy_xxx_0[j] + 0.5 * fl1_fx * tpz_yy_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxxx_0[j];

            tlz_xyy_xxxx_0[j] =
                pa_x[j] * tlz_yy_xxxx_0[j] + 2.0 * fl1_fx * tlz_yy_xxx_0[j] - 0.5 * fl1_fx * tpy_yy_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxxx_0[j];

            tlx_xyy_xxxy_0[j] = pa_x[j] * tlx_yy_xxxy_0[j] + 1.5 * fl1_fx * tlx_yy_xxy_0[j];

            tly_xyy_xxxy_0[j] =
                pa_x[j] * tly_yy_xxxy_0[j] + 1.5 * fl1_fx * tly_yy_xxy_0[j] + 0.5 * fl1_fx * tpz_yy_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxxy_0[j];

            tlz_xyy_xxxy_0[j] =
                pa_x[j] * tlz_yy_xxxy_0[j] + 1.5 * fl1_fx * tlz_yy_xxy_0[j] - 0.5 * fl1_fx * tpy_yy_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxxy_0[j];

            tlx_xyy_xxxz_0[j] = pa_x[j] * tlx_yy_xxxz_0[j] + 1.5 * fl1_fx * tlx_yy_xxz_0[j];

            tly_xyy_xxxz_0[j] =
                pa_x[j] * tly_yy_xxxz_0[j] + 1.5 * fl1_fx * tly_yy_xxz_0[j] + 0.5 * fl1_fx * tpz_yy_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxxz_0[j];

            tlz_xyy_xxxz_0[j] =
                pa_x[j] * tlz_yy_xxxz_0[j] + 1.5 * fl1_fx * tlz_yy_xxz_0[j] - 0.5 * fl1_fx * tpy_yy_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxxz_0[j];

            tlx_xyy_xxyy_0[j] = pa_x[j] * tlx_yy_xxyy_0[j] + fl1_fx * tlx_yy_xyy_0[j];

            tly_xyy_xxyy_0[j] =
                pa_x[j] * tly_yy_xxyy_0[j] + fl1_fx * tly_yy_xyy_0[j] + 0.5 * fl1_fx * tpz_yy_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxyy_0[j];

            tlz_xyy_xxyy_0[j] =
                pa_x[j] * tlz_yy_xxyy_0[j] + fl1_fx * tlz_yy_xyy_0[j] - 0.5 * fl1_fx * tpy_yy_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxyy_0[j];

            tlx_xyy_xxyz_0[j] = pa_x[j] * tlx_yy_xxyz_0[j] + fl1_fx * tlx_yy_xyz_0[j];

            tly_xyy_xxyz_0[j] =
                pa_x[j] * tly_yy_xxyz_0[j] + fl1_fx * tly_yy_xyz_0[j] + 0.5 * fl1_fx * tpz_yy_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxyz_0[j];

            tlz_xyy_xxyz_0[j] =
                pa_x[j] * tlz_yy_xxyz_0[j] + fl1_fx * tlz_yy_xyz_0[j] - 0.5 * fl1_fx * tpy_yy_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_150_200(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (150,200)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 50);

        auto tly_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tlz_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tlx_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 51);

        auto tly_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tlz_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tlx_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 52);

        auto tly_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tlz_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tlx_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 53);

        auto tly_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tlz_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tlx_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 54);

        auto tly_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tlz_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tlx_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 55);

        auto tly_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tlz_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tlx_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 56);

        auto tly_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tlz_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tlx_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 57);

        auto tly_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tlz_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tlx_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 58);

        auto tly_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tlz_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tlx_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 59);

        auto tly_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tlz_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tlx_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 60);

        auto tly_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tlz_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tlx_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 61);

        auto tly_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tlz_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tlx_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 62);

        auto tly_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tlz_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tlx_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 63);

        auto tly_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tlz_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tlx_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 64);

        auto tly_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tlz_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tlx_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 65);

        auto tly_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tlz_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tlx_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 66);

        auto tly_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tlx_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 35);

        auto tly_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tlz_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tlx_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 36);

        auto tly_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tlz_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tlx_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 37);

        auto tly_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tlz_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tlx_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 38);

        auto tly_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tlz_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tlx_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 39);

        auto tly_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tlz_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tlx_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 40);

        auto tly_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tlz_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tlx_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 41);

        auto tly_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tlz_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tlx_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 42);

        auto tly_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tlz_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tlx_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 43);

        auto tly_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tlz_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tlx_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 44);

        auto tly_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tlz_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tlx_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 45);

        auto tly_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tlz_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tlx_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 46);

        auto tly_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tpy_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tpz_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tpy_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tpz_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tpy_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tpz_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tpy_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tpz_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tpy_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tpz_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tpy_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tpz_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tpy_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tpz_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tpy_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tpz_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tpy_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tpz_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tpy_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tpz_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tpy_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tpz_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tpy_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tpz_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tpy_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tpz_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tpy_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tpz_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tpy_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tpz_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tpy_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpz_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tdy_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tdz_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tdy_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tdz_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tdy_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tdz_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tdy_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tdz_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tdy_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tdz_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tdy_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tdy_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tdz_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tdy_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tdz_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tdy_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tdz_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tdy_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tdz_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tdy_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tdz_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tdy_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tdz_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tdy_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tdz_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tdy_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tdz_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tdy_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tdz_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tdy_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tdz_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tdz_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 66);

        // set up pointers to integrals

        auto tlx_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 50);

        auto tly_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tlz_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tlx_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 51);

        auto tly_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tlz_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tlx_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 52);

        auto tly_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tlz_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tlx_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 53);

        auto tly_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tlz_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tlx_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 54);

        auto tly_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tlz_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tlx_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 55);

        auto tly_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tlz_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tlx_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 56);

        auto tly_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tlz_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tlx_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 57);

        auto tly_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tlz_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tlx_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 58);

        auto tly_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tlz_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tlx_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 59);

        auto tly_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tlz_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tlx_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 60);

        auto tly_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tlz_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tlx_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 61);

        auto tly_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tlz_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tlx_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 62);

        auto tly_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tlz_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tlx_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 63);

        auto tly_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tlz_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tlx_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 64);

        auto tly_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tlz_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 64);

        auto tlx_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 65);

        auto tly_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tlz_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tlx_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 66);

        auto tly_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yy_xxzz_0, tdy_yy_xyyy_0, tdy_yy_xyyz_0, \
                                     tdy_yy_xyzz_0, tdy_yy_xzzz_0, tdy_yy_yyyy_0, tdy_yy_yyyz_0, tdy_yy_yyzz_0, \
                                     tdy_yy_yzzz_0, tdy_yy_zzzz_0, tdy_yz_xxxx_0, tdy_yz_xxxy_0, tdy_yz_xxxz_0, \
                                     tdy_yz_xxyy_0, tdy_yz_xxyz_0, tdy_yz_xxzz_0, tdz_yy_xxzz_0, tdz_yy_xyyy_0, \
                                     tdz_yy_xyyz_0, tdz_yy_xyzz_0, tdz_yy_xzzz_0, tdz_yy_yyyy_0, tdz_yy_yyyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzzz_0, tdz_yy_zzzz_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, \
                                     tdz_yz_xxxz_0, tdz_yz_xxyy_0, tdz_yz_xxyz_0, tdz_yz_xxzz_0, tdz_yz_xyyy_0, \
                                     tlx_xyy_xxzz_0, tlx_xyy_xyyy_0, tlx_xyy_xyyz_0, tlx_xyy_xyzz_0, tlx_xyy_xzzz_0, \
                                     tlx_xyy_yyyy_0, tlx_xyy_yyyz_0, tlx_xyy_yyzz_0, tlx_xyy_yzzz_0, tlx_xyy_zzzz_0, \
                                     tlx_xyz_xxxx_0, tlx_xyz_xxxy_0, tlx_xyz_xxxz_0, tlx_xyz_xxyy_0, tlx_xyz_xxyz_0, \
                                     tlx_xyz_xxzz_0, tlx_xyz_xyyy_0, tlx_yy_xxzz_0, tlx_yy_xyyy_0, tlx_yy_xyyz_0, \
                                     tlx_yy_xyzz_0, tlx_yy_xzz_0, tlx_yy_xzzz_0, tlx_yy_yyy_0, tlx_yy_yyyy_0, \
                                     tlx_yy_yyyz_0, tlx_yy_yyz_0, tlx_yy_yyzz_0, tlx_yy_yzz_0, tlx_yy_yzzz_0, \
                                     tlx_yy_zzz_0, tlx_yy_zzzz_0, tlx_yz_xxx_0, tlx_yz_xxxx_0, tlx_yz_xxxy_0, \
                                     tlx_yz_xxxz_0, tlx_yz_xxy_0, tlx_yz_xxyy_0, tlx_yz_xxyz_0, tlx_yz_xxz_0, \
                                     tlx_yz_xxzz_0, tlx_yz_xyy_0, tlx_yz_xyyy_0, tlx_yz_xyz_0, tlx_yz_xzz_0, \
                                     tlx_yz_yyy_0, tly_xyy_xxzz_0, tly_xyy_xyyy_0, tly_xyy_xyyz_0, tly_xyy_xyzz_0, \
                                     tly_xyy_xzzz_0, tly_xyy_yyyy_0, tly_xyy_yyyz_0, tly_xyy_yyzz_0, tly_xyy_yzzz_0, \
                                     tly_xyy_zzzz_0, tly_xyz_xxxx_0, tly_xyz_xxxy_0, tly_xyz_xxxz_0, tly_xyz_xxyy_0, \
                                     tly_xyz_xxyz_0, tly_xyz_xxzz_0, tly_xyz_xyyy_0, tly_yy_xxzz_0, tly_yy_xyyy_0, \
                                     tly_yy_xyyz_0, tly_yy_xyzz_0, tly_yy_xzz_0, tly_yy_xzzz_0, tly_yy_yyy_0, \
                                     tly_yy_yyyy_0, tly_yy_yyyz_0, tly_yy_yyz_0, tly_yy_yyzz_0, tly_yy_yzz_0, \
                                     tly_yy_yzzz_0, tly_yy_zzz_0, tly_yy_zzzz_0, tly_yz_xxx_0, tly_yz_xxxx_0, \
                                     tly_yz_xxxy_0, tly_yz_xxxz_0, tly_yz_xxy_0, tly_yz_xxyy_0, tly_yz_xxyz_0, \
                                     tly_yz_xxz_0, tly_yz_xxzz_0, tly_yz_xyy_0, tly_yz_xyyy_0, tly_yz_xyz_0, \
                                     tly_yz_xzz_0, tly_yz_yyy_0, tlz_xyy_xxzz_0, tlz_xyy_xyyy_0, tlz_xyy_xyyz_0, \
                                     tlz_xyy_xyzz_0, tlz_xyy_xzzz_0, tlz_xyy_yyyy_0, tlz_xyy_yyyz_0, tlz_xyy_yyzz_0, \
                                     tlz_xyy_yzzz_0, tlz_xyy_zzzz_0, tlz_xyz_xxxx_0, tlz_xyz_xxxy_0, tlz_xyz_xxxz_0, \
                                     tlz_xyz_xxyy_0, tlz_xyz_xxyz_0, tlz_xyz_xxzz_0, tlz_yy_xxzz_0, tlz_yy_xyyy_0, \
                                     tlz_yy_xyyz_0, tlz_yy_xyzz_0, tlz_yy_xzz_0, tlz_yy_xzzz_0, tlz_yy_yyy_0, \
                                     tlz_yy_yyyy_0, tlz_yy_yyyz_0, tlz_yy_yyz_0, tlz_yy_yyzz_0, tlz_yy_yzz_0, \
                                     tlz_yy_yzzz_0, tlz_yy_zzz_0, tlz_yy_zzzz_0, tlz_yz_xxx_0, tlz_yz_xxxx_0, \
                                     tlz_yz_xxxy_0, tlz_yz_xxxz_0, tlz_yz_xxy_0, tlz_yz_xxyy_0, tlz_yz_xxyz_0, \
                                     tlz_yz_xxz_0, tlz_yz_xxzz_0, tlz_yz_xyy_0, tlz_yz_xyz_0, tlz_yz_xzz_0, \
                                     tpy_yy_xxzz_0, tpy_yy_xyyy_0, tpy_yy_xyyz_0, tpy_yy_xyzz_0, tpy_yy_xzzz_0, \
                                     tpy_yy_yyyy_0, tpy_yy_yyyz_0, tpy_yy_yyzz_0, tpy_yy_yzzz_0, tpy_yy_zzzz_0, \
                                     tpy_yz_xxxx_0, tpy_yz_xxxy_0, tpy_yz_xxxz_0, tpy_yz_xxyy_0, tpy_yz_xxyz_0, \
                                     tpy_yz_xxzz_0, tpz_yy_xxzz_0, tpz_yy_xyyy_0, tpz_yy_xyyz_0, tpz_yy_xyzz_0, \
                                     tpz_yy_xzzz_0, tpz_yy_yyyy_0, tpz_yy_yyyz_0, tpz_yy_yyzz_0, tpz_yy_yzzz_0, \
                                     tpz_yy_zzzz_0, tpz_yz_xxxx_0, tpz_yz_xxxy_0, tpz_yz_xxxz_0, tpz_yz_xxyy_0, \
                                     tpz_yz_xxyz_0, tpz_yz_xxzz_0, tpz_yz_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xyy_xxzz_0[j] = pa_x[j] * tlx_yy_xxzz_0[j] + fl1_fx * tlx_yy_xzz_0[j];

            tly_xyy_xxzz_0[j] =
                pa_x[j] * tly_yy_xxzz_0[j] + fl1_fx * tly_yy_xzz_0[j] + 0.5 * fl1_fx * tpz_yy_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxzz_0[j];

            tlz_xyy_xxzz_0[j] =
                pa_x[j] * tlz_yy_xxzz_0[j] + fl1_fx * tlz_yy_xzz_0[j] - 0.5 * fl1_fx * tpy_yy_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxzz_0[j];

            tlx_xyy_xyyy_0[j] = pa_x[j] * tlx_yy_xyyy_0[j] + 0.5 * fl1_fx * tlx_yy_yyy_0[j];

            tly_xyy_xyyy_0[j] =
                pa_x[j] * tly_yy_xyyy_0[j] + 0.5 * fl1_fx * tly_yy_yyy_0[j] + 0.5 * fl1_fx * tpz_yy_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xyyy_0[j];

            tlz_xyy_xyyy_0[j] =
                pa_x[j] * tlz_yy_xyyy_0[j] + 0.5 * fl1_fx * tlz_yy_yyy_0[j] - 0.5 * fl1_fx * tpy_yy_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xyyy_0[j];

            tlx_xyy_xyyz_0[j] = pa_x[j] * tlx_yy_xyyz_0[j] + 0.5 * fl1_fx * tlx_yy_yyz_0[j];

            tly_xyy_xyyz_0[j] =
                pa_x[j] * tly_yy_xyyz_0[j] + 0.5 * fl1_fx * tly_yy_yyz_0[j] + 0.5 * fl1_fx * tpz_yy_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xyyz_0[j];

            tlz_xyy_xyyz_0[j] =
                pa_x[j] * tlz_yy_xyyz_0[j] + 0.5 * fl1_fx * tlz_yy_yyz_0[j] - 0.5 * fl1_fx * tpy_yy_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xyyz_0[j];

            tlx_xyy_xyzz_0[j] = pa_x[j] * tlx_yy_xyzz_0[j] + 0.5 * fl1_fx * tlx_yy_yzz_0[j];

            tly_xyy_xyzz_0[j] =
                pa_x[j] * tly_yy_xyzz_0[j] + 0.5 * fl1_fx * tly_yy_yzz_0[j] + 0.5 * fl1_fx * tpz_yy_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xyzz_0[j];

            tlz_xyy_xyzz_0[j] =
                pa_x[j] * tlz_yy_xyzz_0[j] + 0.5 * fl1_fx * tlz_yy_yzz_0[j] - 0.5 * fl1_fx * tpy_yy_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xyzz_0[j];

            tlx_xyy_xzzz_0[j] = pa_x[j] * tlx_yy_xzzz_0[j] + 0.5 * fl1_fx * tlx_yy_zzz_0[j];

            tly_xyy_xzzz_0[j] =
                pa_x[j] * tly_yy_xzzz_0[j] + 0.5 * fl1_fx * tly_yy_zzz_0[j] + 0.5 * fl1_fx * tpz_yy_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xzzz_0[j];

            tlz_xyy_xzzz_0[j] =
                pa_x[j] * tlz_yy_xzzz_0[j] + 0.5 * fl1_fx * tlz_yy_zzz_0[j] - 0.5 * fl1_fx * tpy_yy_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xzzz_0[j];

            tlx_xyy_yyyy_0[j] = pa_x[j] * tlx_yy_yyyy_0[j];

            tly_xyy_yyyy_0[j] = pa_x[j] * tly_yy_yyyy_0[j] + 0.5 * fl1_fx * tpz_yy_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_yy_yyyy_0[j];

            tlz_xyy_yyyy_0[j] = pa_x[j] * tlz_yy_yyyy_0[j] - 0.5 * fl1_fx * tpy_yy_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_yy_yyyy_0[j];

            tlx_xyy_yyyz_0[j] = pa_x[j] * tlx_yy_yyyz_0[j];

            tly_xyy_yyyz_0[j] = pa_x[j] * tly_yy_yyyz_0[j] + 0.5 * fl1_fx * tpz_yy_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yyyz_0[j];

            tlz_xyy_yyyz_0[j] = pa_x[j] * tlz_yy_yyyz_0[j] - 0.5 * fl1_fx * tpy_yy_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yyyz_0[j];

            tlx_xyy_yyzz_0[j] = pa_x[j] * tlx_yy_yyzz_0[j];

            tly_xyy_yyzz_0[j] = pa_x[j] * tly_yy_yyzz_0[j] + 0.5 * fl1_fx * tpz_yy_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yyzz_0[j];

            tlz_xyy_yyzz_0[j] = pa_x[j] * tlz_yy_yyzz_0[j] - 0.5 * fl1_fx * tpy_yy_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yyzz_0[j];

            tlx_xyy_yzzz_0[j] = pa_x[j] * tlx_yy_yzzz_0[j];

            tly_xyy_yzzz_0[j] = pa_x[j] * tly_yy_yzzz_0[j] + 0.5 * fl1_fx * tpz_yy_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yzzz_0[j];

            tlz_xyy_yzzz_0[j] = pa_x[j] * tlz_yy_yzzz_0[j] - 0.5 * fl1_fx * tpy_yy_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yzzz_0[j];

            tlx_xyy_zzzz_0[j] = pa_x[j] * tlx_yy_zzzz_0[j];

            tly_xyy_zzzz_0[j] = pa_x[j] * tly_yy_zzzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_zzzz_0[j];

            tlz_xyy_zzzz_0[j] = pa_x[j] * tlz_yy_zzzz_0[j] - 0.5 * fl1_fx * tpy_yy_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_zzzz_0[j];

            tlx_xyz_xxxx_0[j] = pa_x[j] * tlx_yz_xxxx_0[j] + 2.0 * fl1_fx * tlx_yz_xxx_0[j];

            tly_xyz_xxxx_0[j] =
                pa_x[j] * tly_yz_xxxx_0[j] + 2.0 * fl1_fx * tly_yz_xxx_0[j] + 0.5 * fl1_fx * tpz_yz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxxx_0[j];

            tlz_xyz_xxxx_0[j] =
                pa_x[j] * tlz_yz_xxxx_0[j] + 2.0 * fl1_fx * tlz_yz_xxx_0[j] - 0.5 * fl1_fx * tpy_yz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxxx_0[j];

            tlx_xyz_xxxy_0[j] = pa_x[j] * tlx_yz_xxxy_0[j] + 1.5 * fl1_fx * tlx_yz_xxy_0[j];

            tly_xyz_xxxy_0[j] =
                pa_x[j] * tly_yz_xxxy_0[j] + 1.5 * fl1_fx * tly_yz_xxy_0[j] + 0.5 * fl1_fx * tpz_yz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxxy_0[j];

            tlz_xyz_xxxy_0[j] =
                pa_x[j] * tlz_yz_xxxy_0[j] + 1.5 * fl1_fx * tlz_yz_xxy_0[j] - 0.5 * fl1_fx * tpy_yz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxxy_0[j];

            tlx_xyz_xxxz_0[j] = pa_x[j] * tlx_yz_xxxz_0[j] + 1.5 * fl1_fx * tlx_yz_xxz_0[j];

            tly_xyz_xxxz_0[j] =
                pa_x[j] * tly_yz_xxxz_0[j] + 1.5 * fl1_fx * tly_yz_xxz_0[j] + 0.5 * fl1_fx * tpz_yz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxxz_0[j];

            tlz_xyz_xxxz_0[j] =
                pa_x[j] * tlz_yz_xxxz_0[j] + 1.5 * fl1_fx * tlz_yz_xxz_0[j] - 0.5 * fl1_fx * tpy_yz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxxz_0[j];

            tlx_xyz_xxyy_0[j] = pa_x[j] * tlx_yz_xxyy_0[j] + fl1_fx * tlx_yz_xyy_0[j];

            tly_xyz_xxyy_0[j] =
                pa_x[j] * tly_yz_xxyy_0[j] + fl1_fx * tly_yz_xyy_0[j] + 0.5 * fl1_fx * tpz_yz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxyy_0[j];

            tlz_xyz_xxyy_0[j] =
                pa_x[j] * tlz_yz_xxyy_0[j] + fl1_fx * tlz_yz_xyy_0[j] - 0.5 * fl1_fx * tpy_yz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxyy_0[j];

            tlx_xyz_xxyz_0[j] = pa_x[j] * tlx_yz_xxyz_0[j] + fl1_fx * tlx_yz_xyz_0[j];

            tly_xyz_xxyz_0[j] =
                pa_x[j] * tly_yz_xxyz_0[j] + fl1_fx * tly_yz_xyz_0[j] + 0.5 * fl1_fx * tpz_yz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxyz_0[j];

            tlz_xyz_xxyz_0[j] =
                pa_x[j] * tlz_yz_xxyz_0[j] + fl1_fx * tlz_yz_xyz_0[j] - 0.5 * fl1_fx * tpy_yz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxyz_0[j];

            tlx_xyz_xxzz_0[j] = pa_x[j] * tlx_yz_xxzz_0[j] + fl1_fx * tlx_yz_xzz_0[j];

            tly_xyz_xxzz_0[j] =
                pa_x[j] * tly_yz_xxzz_0[j] + fl1_fx * tly_yz_xzz_0[j] + 0.5 * fl1_fx * tpz_yz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxzz_0[j];

            tlz_xyz_xxzz_0[j] =
                pa_x[j] * tlz_yz_xxzz_0[j] + fl1_fx * tlz_yz_xzz_0[j] - 0.5 * fl1_fx * tpy_yz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxzz_0[j];

            tlx_xyz_xyyy_0[j] = pa_x[j] * tlx_yz_xyyy_0[j] + 0.5 * fl1_fx * tlx_yz_yyy_0[j];

            tly_xyz_xyyy_0[j] =
                pa_x[j] * tly_yz_xyyy_0[j] + 0.5 * fl1_fx * tly_yz_yyy_0[j] + 0.5 * fl1_fx * tpz_yz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xyyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_200_250(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (200,250)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlz_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tlx_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 67);

        auto tly_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tlz_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tlx_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 68);

        auto tly_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tlz_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tlx_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 69);

        auto tly_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tlz_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tlx_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 70);

        auto tly_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tlz_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tlx_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 71);

        auto tly_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tlz_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tlx_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 72);

        auto tly_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tlz_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tlx_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 73);

        auto tly_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tlz_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tlx_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 74);

        auto tly_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tlz_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tlz_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tlx_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 47);

        auto tly_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tlz_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tlx_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 48);

        auto tly_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tlz_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tlx_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 49);

        auto tly_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tlz_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tlx_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 50);

        auto tly_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tlz_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tlx_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 51);

        auto tly_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tlz_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tlx_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 52);

        auto tly_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tlz_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tlx_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 53);

        auto tly_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tlz_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tlx_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 54);

        auto tly_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 55);

        auto tly_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tlz_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tlx_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 56);

        auto tly_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tlz_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tlx_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 57);

        auto tly_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tlz_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tlx_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 58);

        auto tpy_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tpy_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tpz_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tpy_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tpz_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tpy_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tpz_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tpy_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tpz_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tpy_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tpz_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tpy_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tpz_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tpy_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tpz_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tpy_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tpz_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tdy_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tdy_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tdz_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tdy_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tdz_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tdy_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tdz_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tdy_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tdz_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tdy_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tdy_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tdz_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tdy_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tdz_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tdy_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tdz_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82);

        // set up pointers to integrals

        auto tlz_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tlx_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 67);

        auto tly_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tlz_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tlx_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 68);

        auto tly_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tlz_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tlx_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 69);

        auto tly_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tlz_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tlx_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 70);

        auto tly_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tlz_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tlx_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 71);

        auto tly_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tlz_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tlx_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 72);

        auto tly_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tlz_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tlx_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 73);

        auto tly_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tlz_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tlx_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 74);

        auto tly_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tlz_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tlx_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 75);

        auto tly_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tlz_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tlx_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 76);

        auto tly_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tlz_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tlx_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 77);

        auto tly_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tlz_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tlx_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 78);

        auto tly_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tlz_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tlx_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 79);

        auto tly_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tlz_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tlx_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 80);

        auto tly_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tlz_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 80);

        auto tlx_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 81);

        auto tly_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tlz_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tlx_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 82);

        auto tly_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tlz_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tlx_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yz_xyyy_0, tdy_yz_xyyz_0, tdy_yz_xyzz_0, \
                                     tdy_yz_xzzz_0, tdy_yz_yyyy_0, tdy_yz_yyyz_0, tdy_yz_yyzz_0, tdy_yz_yzzz_0, \
                                     tdy_yz_zzzz_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxyy_0, \
                                     tdy_zz_xxyz_0, tdy_zz_xxzz_0, tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdz_yz_xyyz_0, \
                                     tdz_yz_xyzz_0, tdz_yz_xzzz_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, tdz_yz_yyzz_0, \
                                     tdz_yz_yzzz_0, tdz_yz_zzzz_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, \
                                     tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxzz_0, tdz_zz_xyyy_0, tdz_zz_xyyz_0, \
                                     tlx_xyz_xyyz_0, tlx_xyz_xyzz_0, tlx_xyz_xzzz_0, tlx_xyz_yyyy_0, tlx_xyz_yyyz_0, \
                                     tlx_xyz_yyzz_0, tlx_xyz_yzzz_0, tlx_xyz_zzzz_0, tlx_xzz_xxxx_0, tlx_xzz_xxxy_0, \
                                     tlx_xzz_xxxz_0, tlx_xzz_xxyy_0, tlx_xzz_xxyz_0, tlx_xzz_xxzz_0, tlx_xzz_xyyy_0, \
                                     tlx_xzz_xyyz_0, tlx_xzz_xyzz_0, tlx_yz_xyyz_0, tlx_yz_xyzz_0, tlx_yz_xzzz_0, \
                                     tlx_yz_yyyy_0, tlx_yz_yyyz_0, tlx_yz_yyz_0, tlx_yz_yyzz_0, tlx_yz_yzz_0, \
                                     tlx_yz_yzzz_0, tlx_yz_zzz_0, tlx_yz_zzzz_0, tlx_zz_xxx_0, tlx_zz_xxxx_0, \
                                     tlx_zz_xxxy_0, tlx_zz_xxxz_0, tlx_zz_xxy_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, \
                                     tlx_zz_xxz_0, tlx_zz_xxzz_0, tlx_zz_xyy_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, \
                                     tlx_zz_xyz_0, tlx_zz_xyzz_0, tlx_zz_xzz_0, tlx_zz_yyy_0, tlx_zz_yyz_0, \
                                     tlx_zz_yzz_0, tly_xyz_xyyz_0, tly_xyz_xyzz_0, tly_xyz_xzzz_0, tly_xyz_yyyy_0, \
                                     tly_xyz_yyyz_0, tly_xyz_yyzz_0, tly_xyz_yzzz_0, tly_xyz_zzzz_0, tly_xzz_xxxx_0, \
                                     tly_xzz_xxxy_0, tly_xzz_xxxz_0, tly_xzz_xxyy_0, tly_xzz_xxyz_0, tly_xzz_xxzz_0, \
                                     tly_xzz_xyyy_0, tly_xzz_xyyz_0, tly_yz_xyyz_0, tly_yz_xyzz_0, tly_yz_xzzz_0, \
                                     tly_yz_yyyy_0, tly_yz_yyyz_0, tly_yz_yyz_0, tly_yz_yyzz_0, tly_yz_yzz_0, \
                                     tly_yz_yzzz_0, tly_yz_zzz_0, tly_yz_zzzz_0, tly_zz_xxx_0, tly_zz_xxxx_0, \
                                     tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxy_0, tly_zz_xxyy_0, tly_zz_xxyz_0, \
                                     tly_zz_xxz_0, tly_zz_xxzz_0, tly_zz_xyy_0, tly_zz_xyyy_0, tly_zz_xyyz_0, \
                                     tly_zz_xyz_0, tly_zz_xzz_0, tly_zz_yyy_0, tly_zz_yyz_0, tlz_xyz_xyyy_0, \
                                     tlz_xyz_xyyz_0, tlz_xyz_xyzz_0, tlz_xyz_xzzz_0, tlz_xyz_yyyy_0, tlz_xyz_yyyz_0, \
                                     tlz_xyz_yyzz_0, tlz_xyz_yzzz_0, tlz_xyz_zzzz_0, tlz_xzz_xxxx_0, tlz_xzz_xxxy_0, \
                                     tlz_xzz_xxxz_0, tlz_xzz_xxyy_0, tlz_xzz_xxyz_0, tlz_xzz_xxzz_0, tlz_xzz_xyyy_0, \
                                     tlz_xzz_xyyz_0, tlz_yz_xyyy_0, tlz_yz_xyyz_0, tlz_yz_xyzz_0, tlz_yz_xzzz_0, \
                                     tlz_yz_yyy_0, tlz_yz_yyyy_0, tlz_yz_yyyz_0, tlz_yz_yyz_0, tlz_yz_yyzz_0, \
                                     tlz_yz_yzz_0, tlz_yz_yzzz_0, tlz_yz_zzz_0, tlz_yz_zzzz_0, tlz_zz_xxx_0, \
                                     tlz_zz_xxxx_0, tlz_zz_xxxy_0, tlz_zz_xxxz_0, tlz_zz_xxy_0, tlz_zz_xxyy_0, \
                                     tlz_zz_xxyz_0, tlz_zz_xxz_0, tlz_zz_xxzz_0, tlz_zz_xyy_0, tlz_zz_xyyy_0, \
                                     tlz_zz_xyyz_0, tlz_zz_xyz_0, tlz_zz_xzz_0, tlz_zz_yyy_0, tlz_zz_yyz_0, \
                                     tpy_yz_xyyy_0, tpy_yz_xyyz_0, tpy_yz_xyzz_0, tpy_yz_xzzz_0, tpy_yz_yyyy_0, \
                                     tpy_yz_yyyz_0, tpy_yz_yyzz_0, tpy_yz_yzzz_0, tpy_yz_zzzz_0, tpy_zz_xxxx_0, \
                                     tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, \
                                     tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpz_yz_xyyz_0, tpz_yz_xyzz_0, tpz_yz_xzzz_0, \
                                     tpz_yz_yyyy_0, tpz_yz_yyyz_0, tpz_yz_yyzz_0, tpz_yz_yzzz_0, tpz_yz_zzzz_0, \
                                     tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, \
                                     tpz_zz_xxzz_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xyz_xyyy_0[j] =
                pa_x[j] * tlz_yz_xyyy_0[j] + 0.5 * fl1_fx * tlz_yz_yyy_0[j] - 0.5 * fl1_fx * tpy_yz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xyyy_0[j];

            tlx_xyz_xyyz_0[j] = pa_x[j] * tlx_yz_xyyz_0[j] + 0.5 * fl1_fx * tlx_yz_yyz_0[j];

            tly_xyz_xyyz_0[j] =
                pa_x[j] * tly_yz_xyyz_0[j] + 0.5 * fl1_fx * tly_yz_yyz_0[j] + 0.5 * fl1_fx * tpz_yz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xyyz_0[j];

            tlz_xyz_xyyz_0[j] =
                pa_x[j] * tlz_yz_xyyz_0[j] + 0.5 * fl1_fx * tlz_yz_yyz_0[j] - 0.5 * fl1_fx * tpy_yz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xyyz_0[j];

            tlx_xyz_xyzz_0[j] = pa_x[j] * tlx_yz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yz_yzz_0[j];

            tly_xyz_xyzz_0[j] =
                pa_x[j] * tly_yz_xyzz_0[j] + 0.5 * fl1_fx * tly_yz_yzz_0[j] + 0.5 * fl1_fx * tpz_yz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xyzz_0[j];

            tlz_xyz_xyzz_0[j] =
                pa_x[j] * tlz_yz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yz_yzz_0[j] - 0.5 * fl1_fx * tpy_yz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xyzz_0[j];

            tlx_xyz_xzzz_0[j] = pa_x[j] * tlx_yz_xzzz_0[j] + 0.5 * fl1_fx * tlx_yz_zzz_0[j];

            tly_xyz_xzzz_0[j] =
                pa_x[j] * tly_yz_xzzz_0[j] + 0.5 * fl1_fx * tly_yz_zzz_0[j] + 0.5 * fl1_fx * tpz_yz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xzzz_0[j];

            tlz_xyz_xzzz_0[j] =
                pa_x[j] * tlz_yz_xzzz_0[j] + 0.5 * fl1_fx * tlz_yz_zzz_0[j] - 0.5 * fl1_fx * tpy_yz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xzzz_0[j];

            tlx_xyz_yyyy_0[j] = pa_x[j] * tlx_yz_yyyy_0[j];

            tly_xyz_yyyy_0[j] = pa_x[j] * tly_yz_yyyy_0[j] + 0.5 * fl1_fx * tpz_yz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_yz_yyyy_0[j];

            tlz_xyz_yyyy_0[j] = pa_x[j] * tlz_yz_yyyy_0[j] - 0.5 * fl1_fx * tpy_yz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_yz_yyyy_0[j];

            tlx_xyz_yyyz_0[j] = pa_x[j] * tlx_yz_yyyz_0[j];

            tly_xyz_yyyz_0[j] = pa_x[j] * tly_yz_yyyz_0[j] + 0.5 * fl1_fx * tpz_yz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yyyz_0[j];

            tlz_xyz_yyyz_0[j] = pa_x[j] * tlz_yz_yyyz_0[j] - 0.5 * fl1_fx * tpy_yz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yyyz_0[j];

            tlx_xyz_yyzz_0[j] = pa_x[j] * tlx_yz_yyzz_0[j];

            tly_xyz_yyzz_0[j] = pa_x[j] * tly_yz_yyzz_0[j] + 0.5 * fl1_fx * tpz_yz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yyzz_0[j];

            tlz_xyz_yyzz_0[j] = pa_x[j] * tlz_yz_yyzz_0[j] - 0.5 * fl1_fx * tpy_yz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yyzz_0[j];

            tlx_xyz_yzzz_0[j] = pa_x[j] * tlx_yz_yzzz_0[j];

            tly_xyz_yzzz_0[j] = pa_x[j] * tly_yz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yzzz_0[j];

            tlz_xyz_yzzz_0[j] = pa_x[j] * tlz_yz_yzzz_0[j] - 0.5 * fl1_fx * tpy_yz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yzzz_0[j];

            tlx_xyz_zzzz_0[j] = pa_x[j] * tlx_yz_zzzz_0[j];

            tly_xyz_zzzz_0[j] = pa_x[j] * tly_yz_zzzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_zzzz_0[j];

            tlz_xyz_zzzz_0[j] = pa_x[j] * tlz_yz_zzzz_0[j] - 0.5 * fl1_fx * tpy_yz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_zzzz_0[j];

            tlx_xzz_xxxx_0[j] = pa_x[j] * tlx_zz_xxxx_0[j] + 2.0 * fl1_fx * tlx_zz_xxx_0[j];

            tly_xzz_xxxx_0[j] =
                pa_x[j] * tly_zz_xxxx_0[j] + 2.0 * fl1_fx * tly_zz_xxx_0[j] + 0.5 * fl1_fx * tpz_zz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxxx_0[j];

            tlz_xzz_xxxx_0[j] =
                pa_x[j] * tlz_zz_xxxx_0[j] + 2.0 * fl1_fx * tlz_zz_xxx_0[j] - 0.5 * fl1_fx * tpy_zz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxxx_0[j];

            tlx_xzz_xxxy_0[j] = pa_x[j] * tlx_zz_xxxy_0[j] + 1.5 * fl1_fx * tlx_zz_xxy_0[j];

            tly_xzz_xxxy_0[j] =
                pa_x[j] * tly_zz_xxxy_0[j] + 1.5 * fl1_fx * tly_zz_xxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxxy_0[j];

            tlz_xzz_xxxy_0[j] =
                pa_x[j] * tlz_zz_xxxy_0[j] + 1.5 * fl1_fx * tlz_zz_xxy_0[j] - 0.5 * fl1_fx * tpy_zz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxxy_0[j];

            tlx_xzz_xxxz_0[j] = pa_x[j] * tlx_zz_xxxz_0[j] + 1.5 * fl1_fx * tlx_zz_xxz_0[j];

            tly_xzz_xxxz_0[j] =
                pa_x[j] * tly_zz_xxxz_0[j] + 1.5 * fl1_fx * tly_zz_xxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxxz_0[j];

            tlz_xzz_xxxz_0[j] =
                pa_x[j] * tlz_zz_xxxz_0[j] + 1.5 * fl1_fx * tlz_zz_xxz_0[j] - 0.5 * fl1_fx * tpy_zz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxxz_0[j];

            tlx_xzz_xxyy_0[j] = pa_x[j] * tlx_zz_xxyy_0[j] + fl1_fx * tlx_zz_xyy_0[j];

            tly_xzz_xxyy_0[j] =
                pa_x[j] * tly_zz_xxyy_0[j] + fl1_fx * tly_zz_xyy_0[j] + 0.5 * fl1_fx * tpz_zz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxyy_0[j];

            tlz_xzz_xxyy_0[j] =
                pa_x[j] * tlz_zz_xxyy_0[j] + fl1_fx * tlz_zz_xyy_0[j] - 0.5 * fl1_fx * tpy_zz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxyy_0[j];

            tlx_xzz_xxyz_0[j] = pa_x[j] * tlx_zz_xxyz_0[j] + fl1_fx * tlx_zz_xyz_0[j];

            tly_xzz_xxyz_0[j] =
                pa_x[j] * tly_zz_xxyz_0[j] + fl1_fx * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxyz_0[j];

            tlz_xzz_xxyz_0[j] =
                pa_x[j] * tlz_zz_xxyz_0[j] + fl1_fx * tlz_zz_xyz_0[j] - 0.5 * fl1_fx * tpy_zz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxyz_0[j];

            tlx_xzz_xxzz_0[j] = pa_x[j] * tlx_zz_xxzz_0[j] + fl1_fx * tlx_zz_xzz_0[j];

            tly_xzz_xxzz_0[j] =
                pa_x[j] * tly_zz_xxzz_0[j] + fl1_fx * tly_zz_xzz_0[j] + 0.5 * fl1_fx * tpz_zz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxzz_0[j];

            tlz_xzz_xxzz_0[j] =
                pa_x[j] * tlz_zz_xxzz_0[j] + fl1_fx * tlz_zz_xzz_0[j] - 0.5 * fl1_fx * tpy_zz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxzz_0[j];

            tlx_xzz_xyyy_0[j] = pa_x[j] * tlx_zz_xyyy_0[j] + 0.5 * fl1_fx * tlx_zz_yyy_0[j];

            tly_xzz_xyyy_0[j] =
                pa_x[j] * tly_zz_xyyy_0[j] + 0.5 * fl1_fx * tly_zz_yyy_0[j] + 0.5 * fl1_fx * tpz_zz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xyyy_0[j];

            tlz_xzz_xyyy_0[j] =
                pa_x[j] * tlz_zz_xyyy_0[j] + 0.5 * fl1_fx * tlz_zz_yyy_0[j] - 0.5 * fl1_fx * tpy_zz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xyyy_0[j];

            tlx_xzz_xyyz_0[j] = pa_x[j] * tlx_zz_xyyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyz_0[j];

            tly_xzz_xyyz_0[j] =
                pa_x[j] * tly_zz_xyyz_0[j] + 0.5 * fl1_fx * tly_zz_yyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xyyz_0[j];

            tlz_xzz_xyyz_0[j] =
                pa_x[j] * tlz_zz_xyyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyz_0[j] - 0.5 * fl1_fx * tpy_zz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xyyz_0[j];

            tlx_xzz_xyzz_0[j] = pa_x[j] * tlx_zz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zz_yzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_250_300(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (250,300)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 45);

        auto tly_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tlz_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tlx_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 46);

        auto tly_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tlz_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tlx_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 47);

        auto tly_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tlz_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tlx_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 48);

        auto tly_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tlz_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tlx_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 49);

        auto tly_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tlz_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tlx_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 50);

        auto tly_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tlz_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tlx_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 51);

        auto tly_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tlz_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tlx_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 52);

        auto tly_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tlz_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tlx_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 53);

        auto tly_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tlz_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tlx_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 54);

        auto tly_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tlz_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

        auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89);

        auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89);

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

        auto tlx_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 30);

        auto tly_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tlz_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tlx_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 31);

        auto tly_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tlz_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tlx_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 32);

        auto tly_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tlz_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tlx_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 33);

        auto tly_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tlz_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tlx_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 34);

        auto tly_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tlz_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tlx_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 35);

        auto tly_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tlz_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tly_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 59);

        auto tly_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tlz_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto tpx_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 45);

        auto tpz_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tpx_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 46);

        auto tpz_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tpx_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 47);

        auto tpz_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tpx_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 48);

        auto tpz_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tpx_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tpx_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 50);

        auto tpz_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tpx_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 51);

        auto tpz_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tpx_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 52);

        auto tpz_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tpx_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 53);

        auto tpz_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tpx_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 54);

        auto tpz_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tdx_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 45);

        auto tdz_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tdx_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 46);

        auto tdz_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tdx_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 47);

        auto tdz_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tdx_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 48);

        auto tdz_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tdx_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 49);

        auto tdz_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tdx_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 50);

        auto tdz_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tdx_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 51);

        auto tdz_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tdx_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 52);

        auto tdz_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tdx_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 53);

        auto tdz_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tdx_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 54);

        auto tdz_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89);

        // set up pointers to integrals

        auto tly_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tlz_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tlx_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 84);

        auto tly_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tlz_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tlx_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 85);

        auto tly_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tlz_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tlx_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 86);

        auto tly_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tlz_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tlx_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 87);

        auto tly_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tlz_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tlx_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 88);

        auto tly_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tlz_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tlx_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 89);

        auto tly_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tlz_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tlx_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 90);

        auto tly_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tlz_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tlx_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 91);

        auto tly_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tlz_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tlx_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 92);

        auto tly_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tlz_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tlx_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 93);

        auto tly_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tlz_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tlx_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 94);

        auto tly_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tlz_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tlx_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 95);

        auto tly_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tlz_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tlx_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 96);

        auto tly_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tlz_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tlx_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 97);

        auto tly_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tlz_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tlx_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 98);

        auto tly_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tlz_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tlx_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 99);

        auto tly_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tlz_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_yy_xxxx_0, tdx_yy_xxxy_0, tdx_yy_xxxz_0, \
                                     tdx_yy_xxyy_0, tdx_yy_xxyz_0, tdx_yy_xxzz_0, tdx_yy_xyyy_0, tdx_yy_xyyz_0, \
                                     tdx_yy_xyzz_0, tdx_yy_xzzz_0, tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, \
                                     tdy_zz_yyyz_0, tdy_zz_yyzz_0, tdy_zz_yzzz_0, tdy_zz_zzzz_0, tdz_yy_xxxx_0, \
                                     tdz_yy_xxxy_0, tdz_yy_xxxz_0, tdz_yy_xxyy_0, tdz_yy_xxyz_0, tdz_yy_xxzz_0, \
                                     tdz_yy_xyyy_0, tdz_yy_xyyz_0, tdz_yy_xyzz_0, tdz_yy_xzzz_0, tdz_zz_xyzz_0, \
                                     tdz_zz_xzzz_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyzz_0, tdz_zz_yzzz_0, \
                                     tdz_zz_zzzz_0, tlx_xzz_xzzz_0, tlx_xzz_yyyy_0, tlx_xzz_yyyz_0, tlx_xzz_yyzz_0, \
                                     tlx_xzz_yzzz_0, tlx_xzz_zzzz_0, tlx_y_xxxx_0, tlx_y_xxxy_0, tlx_y_xxxz_0, \
                                     tlx_y_xxyy_0, tlx_y_xxyz_0, tlx_y_xxzz_0, tlx_y_xyyy_0, tlx_y_xyyz_0, tlx_y_xyzz_0, \
                                     tlx_y_xzzz_0, tlx_yy_xxx_0, tlx_yy_xxxx_0, tlx_yy_xxxy_0, tlx_yy_xxxz_0, \
                                     tlx_yy_xxy_0, tlx_yy_xxyy_0, tlx_yy_xxyz_0, tlx_yy_xxz_0, tlx_yy_xxzz_0, \
                                     tlx_yy_xyy_0, tlx_yy_xyyy_0, tlx_yy_xyyz_0, tlx_yy_xyz_0, tlx_yy_xyzz_0, \
                                     tlx_yy_xzz_0, tlx_yy_xzzz_0, tlx_yyy_xxxx_0, tlx_yyy_xxxy_0, tlx_yyy_xxxz_0, \
                                     tlx_yyy_xxyy_0, tlx_yyy_xxyz_0, tlx_yyy_xxzz_0, tlx_yyy_xyyy_0, tlx_yyy_xyyz_0, \
                                     tlx_yyy_xyzz_0, tlx_yyy_xzzz_0, tlx_zz_xzzz_0, tlx_zz_yyyy_0, tlx_zz_yyyz_0, \
                                     tlx_zz_yyzz_0, tlx_zz_yzzz_0, tlx_zz_zzz_0, tlx_zz_zzzz_0, tly_xzz_xyzz_0, \
                                     tly_xzz_xzzz_0, tly_xzz_yyyy_0, tly_xzz_yyyz_0, tly_xzz_yyzz_0, tly_xzz_yzzz_0, \
                                     tly_xzz_zzzz_0, tly_y_xxxx_0, tly_y_xxxy_0, tly_y_xxxz_0, tly_y_xxyy_0, tly_y_xxyz_0, \
                                     tly_y_xxzz_0, tly_y_xyyy_0, tly_y_xyyz_0, tly_y_xyzz_0, tly_y_xzzz_0, tly_yy_xxx_0, \
                                     tly_yy_xxxx_0, tly_yy_xxxy_0, tly_yy_xxxz_0, tly_yy_xxy_0, tly_yy_xxyy_0, \
                                     tly_yy_xxyz_0, tly_yy_xxz_0, tly_yy_xxzz_0, tly_yy_xyy_0, tly_yy_xyyy_0, \
                                     tly_yy_xyyz_0, tly_yy_xyz_0, tly_yy_xyzz_0, tly_yy_xzz_0, tly_yy_xzzz_0, \
                                     tly_yyy_xxxx_0, tly_yyy_xxxy_0, tly_yyy_xxxz_0, tly_yyy_xxyy_0, tly_yyy_xxyz_0, \
                                     tly_yyy_xxzz_0, tly_yyy_xyyy_0, tly_yyy_xyyz_0, tly_yyy_xyzz_0, tly_yyy_xzzz_0, \
                                     tly_zz_xyzz_0, tly_zz_xzzz_0, tly_zz_yyyy_0, tly_zz_yyyz_0, tly_zz_yyzz_0, \
                                     tly_zz_yzz_0, tly_zz_yzzz_0, tly_zz_zzz_0, tly_zz_zzzz_0, tlz_xzz_xyzz_0, \
                                     tlz_xzz_xzzz_0, tlz_xzz_yyyy_0, tlz_xzz_yyyz_0, tlz_xzz_yyzz_0, tlz_xzz_yzzz_0, \
                                     tlz_xzz_zzzz_0, tlz_y_xxxx_0, tlz_y_xxxy_0, tlz_y_xxxz_0, tlz_y_xxyy_0, tlz_y_xxyz_0, \
                                     tlz_y_xxzz_0, tlz_y_xyyy_0, tlz_y_xyyz_0, tlz_y_xyzz_0, tlz_y_xzzz_0, tlz_yy_xxx_0, \
                                     tlz_yy_xxxx_0, tlz_yy_xxxy_0, tlz_yy_xxxz_0, tlz_yy_xxy_0, tlz_yy_xxyy_0, \
                                     tlz_yy_xxyz_0, tlz_yy_xxz_0, tlz_yy_xxzz_0, tlz_yy_xyy_0, tlz_yy_xyyy_0, \
                                     tlz_yy_xyyz_0, tlz_yy_xyz_0, tlz_yy_xyzz_0, tlz_yy_xzz_0, tlz_yy_xzzz_0, \
                                     tlz_yyy_xxxx_0, tlz_yyy_xxxy_0, tlz_yyy_xxxz_0, tlz_yyy_xxyy_0, tlz_yyy_xxyz_0, \
                                     tlz_yyy_xxzz_0, tlz_yyy_xyyy_0, tlz_yyy_xyyz_0, tlz_yyy_xyzz_0, tlz_yyy_xzzz_0, \
                                     tlz_zz_xyzz_0, tlz_zz_xzzz_0, tlz_zz_yyyy_0, tlz_zz_yyyz_0, tlz_zz_yyzz_0, \
                                     tlz_zz_yzz_0, tlz_zz_yzzz_0, tlz_zz_zzz_0, tlz_zz_zzzz_0, tpx_yy_xxxx_0, \
                                     tpx_yy_xxxy_0, tpx_yy_xxxz_0, tpx_yy_xxyy_0, tpx_yy_xxyz_0, tpx_yy_xxzz_0, \
                                     tpx_yy_xyyy_0, tpx_yy_xyyz_0, tpx_yy_xyzz_0, tpx_yy_xzzz_0, tpy_zz_xyzz_0, \
                                     tpy_zz_xzzz_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, tpy_zz_yyzz_0, tpy_zz_yzzz_0, \
                                     tpy_zz_zzzz_0, tpz_yy_xxxx_0, tpz_yy_xxxy_0, tpz_yy_xxxz_0, tpz_yy_xxyy_0, \
                                     tpz_yy_xxyz_0, tpz_yy_xxzz_0, tpz_yy_xyyy_0, tpz_yy_xyyz_0, tpz_yy_xyzz_0, \
                                     tpz_yy_xzzz_0, tpz_zz_xyzz_0, tpz_zz_xzzz_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, \
                                     tpz_zz_yyzz_0, tpz_zz_yzzz_0, tpz_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xzz_xyzz_0[j] =
                pa_x[j] * tly_zz_xyzz_0[j] + 0.5 * fl1_fx * tly_zz_yzz_0[j] + 0.5 * fl1_fx * tpz_zz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xyzz_0[j];

            tlz_xzz_xyzz_0[j] =
                pa_x[j] * tlz_zz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zz_yzz_0[j] - 0.5 * fl1_fx * tpy_zz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xyzz_0[j];

            tlx_xzz_xzzz_0[j] = pa_x[j] * tlx_zz_xzzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzz_0[j];

            tly_xzz_xzzz_0[j] =
                pa_x[j] * tly_zz_xzzz_0[j] + 0.5 * fl1_fx * tly_zz_zzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xzzz_0[j];

            tlz_xzz_xzzz_0[j] =
                pa_x[j] * tlz_zz_xzzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzz_0[j] - 0.5 * fl1_fx * tpy_zz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xzzz_0[j];

            tlx_xzz_yyyy_0[j] = pa_x[j] * tlx_zz_yyyy_0[j];

            tly_xzz_yyyy_0[j] = pa_x[j] * tly_zz_yyyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_zz_yyyy_0[j];

            tlz_xzz_yyyy_0[j] = pa_x[j] * tlz_zz_yyyy_0[j] - 0.5 * fl1_fx * tpy_zz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_zz_yyyy_0[j];

            tlx_xzz_yyyz_0[j] = pa_x[j] * tlx_zz_yyyz_0[j];

            tly_xzz_yyyz_0[j] = pa_x[j] * tly_zz_yyyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yyyz_0[j];

            tlz_xzz_yyyz_0[j] = pa_x[j] * tlz_zz_yyyz_0[j] - 0.5 * fl1_fx * tpy_zz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yyyz_0[j];

            tlx_xzz_yyzz_0[j] = pa_x[j] * tlx_zz_yyzz_0[j];

            tly_xzz_yyzz_0[j] = pa_x[j] * tly_zz_yyzz_0[j] + 0.5 * fl1_fx * tpz_zz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yyzz_0[j];

            tlz_xzz_yyzz_0[j] = pa_x[j] * tlz_zz_yyzz_0[j] - 0.5 * fl1_fx * tpy_zz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yyzz_0[j];

            tlx_xzz_yzzz_0[j] = pa_x[j] * tlx_zz_yzzz_0[j];

            tly_xzz_yzzz_0[j] = pa_x[j] * tly_zz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yzzz_0[j];

            tlz_xzz_yzzz_0[j] = pa_x[j] * tlz_zz_yzzz_0[j] - 0.5 * fl1_fx * tpy_zz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yzzz_0[j];

            tlx_xzz_zzzz_0[j] = pa_x[j] * tlx_zz_zzzz_0[j];

            tly_xzz_zzzz_0[j] = pa_x[j] * tly_zz_zzzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_zzzz_0[j];

            tlz_xzz_zzzz_0[j] = pa_x[j] * tlz_zz_zzzz_0[j] - 0.5 * fl1_fx * tpy_zz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_zzzz_0[j];

            tlx_yyy_xxxx_0[j] =
                pa_y[j] * tlx_yy_xxxx_0[j] + fl1_fx * tlx_y_xxxx_0[j] - 0.5 * fl1_fx * tpz_yy_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxxx_0[j];

            tly_yyy_xxxx_0[j] = pa_y[j] * tly_yy_xxxx_0[j] + fl1_fx * tly_y_xxxx_0[j];

            tlz_yyy_xxxx_0[j] =
                pa_y[j] * tlz_yy_xxxx_0[j] + fl1_fx * tlz_y_xxxx_0[j] + 0.5 * fl1_fx * tpx_yy_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxxx_0[j];

            tlx_yyy_xxxy_0[j] = pa_y[j] * tlx_yy_xxxy_0[j] + fl1_fx * tlx_y_xxxy_0[j] + 0.5 * fl1_fx * tlx_yy_xxx_0[j] -
                                0.5 * fl1_fx * tpz_yy_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxxy_0[j];

            tly_yyy_xxxy_0[j] = pa_y[j] * tly_yy_xxxy_0[j] + fl1_fx * tly_y_xxxy_0[j] + 0.5 * fl1_fx * tly_yy_xxx_0[j];

            tlz_yyy_xxxy_0[j] = pa_y[j] * tlz_yy_xxxy_0[j] + fl1_fx * tlz_y_xxxy_0[j] + 0.5 * fl1_fx * tlz_yy_xxx_0[j] +
                                0.5 * fl1_fx * tpx_yy_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxxy_0[j];

            tlx_yyy_xxxz_0[j] =
                pa_y[j] * tlx_yy_xxxz_0[j] + fl1_fx * tlx_y_xxxz_0[j] - 0.5 * fl1_fx * tpz_yy_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxxz_0[j];

            tly_yyy_xxxz_0[j] = pa_y[j] * tly_yy_xxxz_0[j] + fl1_fx * tly_y_xxxz_0[j];

            tlz_yyy_xxxz_0[j] =
                pa_y[j] * tlz_yy_xxxz_0[j] + fl1_fx * tlz_y_xxxz_0[j] + 0.5 * fl1_fx * tpx_yy_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxxz_0[j];

            tlx_yyy_xxyy_0[j] = pa_y[j] * tlx_yy_xxyy_0[j] + fl1_fx * tlx_y_xxyy_0[j] + fl1_fx * tlx_yy_xxy_0[j] - 0.5 * fl1_fx * tpz_yy_xxyy_0[j] -
                                fl1_fx * fl1_fgb * tdz_yy_xxyy_0[j];

            tly_yyy_xxyy_0[j] = pa_y[j] * tly_yy_xxyy_0[j] + fl1_fx * tly_y_xxyy_0[j] + fl1_fx * tly_yy_xxy_0[j];

            tlz_yyy_xxyy_0[j] = pa_y[j] * tlz_yy_xxyy_0[j] + fl1_fx * tlz_y_xxyy_0[j] + fl1_fx * tlz_yy_xxy_0[j] + 0.5 * fl1_fx * tpx_yy_xxyy_0[j] +
                                fl1_fx * fl1_fgb * tdx_yy_xxyy_0[j];

            tlx_yyy_xxyz_0[j] = pa_y[j] * tlx_yy_xxyz_0[j] + fl1_fx * tlx_y_xxyz_0[j] + 0.5 * fl1_fx * tlx_yy_xxz_0[j] -
                                0.5 * fl1_fx * tpz_yy_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxyz_0[j];

            tly_yyy_xxyz_0[j] = pa_y[j] * tly_yy_xxyz_0[j] + fl1_fx * tly_y_xxyz_0[j] + 0.5 * fl1_fx * tly_yy_xxz_0[j];

            tlz_yyy_xxyz_0[j] = pa_y[j] * tlz_yy_xxyz_0[j] + fl1_fx * tlz_y_xxyz_0[j] + 0.5 * fl1_fx * tlz_yy_xxz_0[j] +
                                0.5 * fl1_fx * tpx_yy_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxyz_0[j];

            tlx_yyy_xxzz_0[j] =
                pa_y[j] * tlx_yy_xxzz_0[j] + fl1_fx * tlx_y_xxzz_0[j] - 0.5 * fl1_fx * tpz_yy_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxzz_0[j];

            tly_yyy_xxzz_0[j] = pa_y[j] * tly_yy_xxzz_0[j] + fl1_fx * tly_y_xxzz_0[j];

            tlz_yyy_xxzz_0[j] =
                pa_y[j] * tlz_yy_xxzz_0[j] + fl1_fx * tlz_y_xxzz_0[j] + 0.5 * fl1_fx * tpx_yy_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxzz_0[j];

            tlx_yyy_xyyy_0[j] = pa_y[j] * tlx_yy_xyyy_0[j] + fl1_fx * tlx_y_xyyy_0[j] + 1.5 * fl1_fx * tlx_yy_xyy_0[j] -
                                0.5 * fl1_fx * tpz_yy_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_yy_xyyy_0[j];

            tly_yyy_xyyy_0[j] = pa_y[j] * tly_yy_xyyy_0[j] + fl1_fx * tly_y_xyyy_0[j] + 1.5 * fl1_fx * tly_yy_xyy_0[j];

            tlz_yyy_xyyy_0[j] = pa_y[j] * tlz_yy_xyyy_0[j] + fl1_fx * tlz_y_xyyy_0[j] + 1.5 * fl1_fx * tlz_yy_xyy_0[j] +
                                0.5 * fl1_fx * tpx_yy_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_yy_xyyy_0[j];

            tlx_yyy_xyyz_0[j] = pa_y[j] * tlx_yy_xyyz_0[j] + fl1_fx * tlx_y_xyyz_0[j] + fl1_fx * tlx_yy_xyz_0[j] - 0.5 * fl1_fx * tpz_yy_xyyz_0[j] -
                                fl1_fx * fl1_fgb * tdz_yy_xyyz_0[j];

            tly_yyy_xyyz_0[j] = pa_y[j] * tly_yy_xyyz_0[j] + fl1_fx * tly_y_xyyz_0[j] + fl1_fx * tly_yy_xyz_0[j];

            tlz_yyy_xyyz_0[j] = pa_y[j] * tlz_yy_xyyz_0[j] + fl1_fx * tlz_y_xyyz_0[j] + fl1_fx * tlz_yy_xyz_0[j] + 0.5 * fl1_fx * tpx_yy_xyyz_0[j] +
                                fl1_fx * fl1_fgb * tdx_yy_xyyz_0[j];

            tlx_yyy_xyzz_0[j] = pa_y[j] * tlx_yy_xyzz_0[j] + fl1_fx * tlx_y_xyzz_0[j] + 0.5 * fl1_fx * tlx_yy_xzz_0[j] -
                                0.5 * fl1_fx * tpz_yy_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xyzz_0[j];

            tly_yyy_xyzz_0[j] = pa_y[j] * tly_yy_xyzz_0[j] + fl1_fx * tly_y_xyzz_0[j] + 0.5 * fl1_fx * tly_yy_xzz_0[j];

            tlz_yyy_xyzz_0[j] = pa_y[j] * tlz_yy_xyzz_0[j] + fl1_fx * tlz_y_xyzz_0[j] + 0.5 * fl1_fx * tlz_yy_xzz_0[j] +
                                0.5 * fl1_fx * tpx_yy_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xyzz_0[j];

            tlx_yyy_xzzz_0[j] =
                pa_y[j] * tlx_yy_xzzz_0[j] + fl1_fx * tlx_y_xzzz_0[j] - 0.5 * fl1_fx * tpz_yy_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xzzz_0[j];

            tly_yyy_xzzz_0[j] = pa_y[j] * tly_yy_xzzz_0[j] + fl1_fx * tly_y_xzzz_0[j];

            tlz_yyy_xzzz_0[j] =
                pa_y[j] * tlz_yy_xzzz_0[j] + fl1_fx * tlz_y_xzzz_0[j] + 0.5 * fl1_fx * tpx_yy_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_300_350(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (300,350)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 55);

        auto tly_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tlz_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tlx_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 56);

        auto tly_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tlz_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tlx_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 57);

        auto tly_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tlz_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tlx_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 58);

        auto tly_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tlz_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tlx_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 59);

        auto tly_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tlz_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tlx_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 60);

        auto tly_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tlz_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tlx_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 61);

        auto tly_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tlz_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tlx_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 62);

        auto tly_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tlz_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tlx_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 63);

        auto tly_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tlz_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tlx_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 64);

        auto tly_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tlz_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tlx_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 65);

        auto tly_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tlz_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tlx_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 66);

        auto tly_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tlz_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tlx_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 67);

        auto tly_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tlz_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tlx_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 68);

        auto tly_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tlz_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tlx_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 69);

        auto tly_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tlz_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tlx_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 70);

        auto tly_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tlz_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tlx_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 71);

        auto tly_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 71);

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

        auto tlx_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 36);

        auto tly_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tlz_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tlx_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 37);

        auto tly_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tlz_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tlx_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 38);

        auto tly_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tlz_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tlx_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 39);

        auto tly_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tlz_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tlx_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 40);

        auto tly_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tlz_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tlx_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 41);

        auto tly_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tlz_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tlx_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 42);

        auto tly_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tlz_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tlx_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 43);

        auto tly_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tlz_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tlx_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 44);

        auto tly_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tlz_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tlx_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 45);

        auto tly_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tlz_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tlx_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 46);

        auto tly_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tlz_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tlx_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 47);

        auto tly_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tpx_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 55);

        auto tpz_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tpx_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 56);

        auto tpz_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tpx_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 57);

        auto tpz_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tpx_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 58);

        auto tpz_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tpx_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 59);

        auto tpz_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tpx_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 60);

        auto tpz_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tpx_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 61);

        auto tpz_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tpx_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 62);

        auto tpz_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tpx_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 63);

        auto tpz_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tpx_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 64);

        auto tpz_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tpx_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 66);

        auto tpz_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tpx_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 67);

        auto tpz_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tpx_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 68);

        auto tpz_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tpx_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 69);

        auto tpz_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tpx_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 70);

        auto tpz_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tpz_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tdx_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 55);

        auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tdx_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 56);

        auto tdz_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tdx_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 57);

        auto tdz_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tdx_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 58);

        auto tdz_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tdx_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 59);

        auto tdz_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tdx_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 60);

        auto tdz_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tdx_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 61);

        auto tdz_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tdx_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 62);

        auto tdz_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tdx_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 63);

        auto tdz_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tdx_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 64);

        auto tdz_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tdx_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 65);

        auto tdz_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tdx_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 66);

        auto tdz_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tdx_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 67);

        auto tdz_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tdx_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 68);

        auto tdz_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tdx_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 69);

        auto tdz_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tdx_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 70);

        auto tdz_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71);

        // set up pointers to integrals

        auto tlx_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 100);

        auto tly_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tlz_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tlx_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 101);

        auto tly_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tlz_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tlx_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 102);

        auto tly_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tlz_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tlx_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 103);

        auto tly_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tlz_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tlx_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 104);

        auto tly_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tlz_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tlx_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 105);

        auto tly_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tlz_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tlx_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 106);

        auto tly_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tlz_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tlx_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 107);

        auto tly_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tlz_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tlx_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 108);

        auto tly_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tlz_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tlx_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 109);

        auto tly_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tlz_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tlx_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 110);

        auto tly_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tlz_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tlx_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 111);

        auto tly_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tlz_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tlx_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 112);

        auto tly_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tlz_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tlx_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 113);

        auto tly_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tlz_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tlx_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 114);

        auto tly_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tlz_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tlx_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 115);

        auto tly_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tlz_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tlx_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 116);

        auto tly_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 116);

        // Batch of Integrals (300,350)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yy_yyyy_0, tdx_yy_yyyz_0, tdx_yy_yyzz_0, \
                                     tdx_yy_yzzz_0, tdx_yy_zzzz_0, tdx_yz_xxxx_0, tdx_yz_xxxy_0, tdx_yz_xxxz_0, \
                                     tdx_yz_xxyy_0, tdx_yz_xxyz_0, tdx_yz_xxzz_0, tdx_yz_xyyy_0, tdx_yz_xyyz_0, \
                                     tdx_yz_xyzz_0, tdx_yz_xzzz_0, tdx_yz_yyyy_0, tdz_yy_yyyy_0, tdz_yy_yyyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzzz_0, tdz_yy_zzzz_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, \
                                     tdz_yz_xxxz_0, tdz_yz_xxyy_0, tdz_yz_xxyz_0, tdz_yz_xxzz_0, tdz_yz_xyyy_0, \
                                     tdz_yz_xyyz_0, tdz_yz_xyzz_0, tdz_yz_xzzz_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, \
                                     tlx_y_yyyy_0, tlx_y_yyyz_0, tlx_y_yyzz_0, tlx_y_yzzz_0, tlx_y_zzzz_0, tlx_yy_yyy_0, \
                                     tlx_yy_yyyy_0, tlx_yy_yyyz_0, tlx_yy_yyz_0, tlx_yy_yyzz_0, tlx_yy_yzz_0, \
                                     tlx_yy_yzzz_0, tlx_yy_zzz_0, tlx_yy_zzzz_0, tlx_yyy_yyyy_0, tlx_yyy_yyyz_0, \
                                     tlx_yyy_yyzz_0, tlx_yyy_yzzz_0, tlx_yyy_zzzz_0, tlx_yyz_xxxx_0, tlx_yyz_xxxy_0, \
                                     tlx_yyz_xxxz_0, tlx_yyz_xxyy_0, tlx_yyz_xxyz_0, tlx_yyz_xxzz_0, tlx_yyz_xyyy_0, \
                                     tlx_yyz_xyyz_0, tlx_yyz_xyzz_0, tlx_yyz_xzzz_0, tlx_yyz_yyyy_0, tlx_yyz_yyyz_0, \
                                     tlx_yz_xxx_0, tlx_yz_xxxx_0, tlx_yz_xxxy_0, tlx_yz_xxxz_0, tlx_yz_xxy_0, \
                                     tlx_yz_xxyy_0, tlx_yz_xxyz_0, tlx_yz_xxz_0, tlx_yz_xxzz_0, tlx_yz_xyy_0, \
                                     tlx_yz_xyyy_0, tlx_yz_xyyz_0, tlx_yz_xyz_0, tlx_yz_xyzz_0, tlx_yz_xzz_0, \
                                     tlx_yz_xzzz_0, tlx_yz_yyy_0, tlx_yz_yyyy_0, tlx_yz_yyyz_0, tlx_yz_yyz_0, \
                                     tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, tlx_z_xxyy_0, tlx_z_xxyz_0, tlx_z_xxzz_0, \
                                     tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyzz_0, tlx_z_xzzz_0, tlx_z_yyyy_0, tlx_z_yyyz_0, \
                                     tly_y_yyyy_0, tly_y_yyyz_0, tly_y_yyzz_0, tly_y_yzzz_0, tly_y_zzzz_0, tly_yy_yyy_0, \
                                     tly_yy_yyyy_0, tly_yy_yyyz_0, tly_yy_yyz_0, tly_yy_yyzz_0, tly_yy_yzz_0, \
                                     tly_yy_yzzz_0, tly_yy_zzz_0, tly_yy_zzzz_0, tly_yyy_yyyy_0, tly_yyy_yyyz_0, \
                                     tly_yyy_yyzz_0, tly_yyy_yzzz_0, tly_yyy_zzzz_0, tly_yyz_xxxx_0, tly_yyz_xxxy_0, \
                                     tly_yyz_xxxz_0, tly_yyz_xxyy_0, tly_yyz_xxyz_0, tly_yyz_xxzz_0, tly_yyz_xyyy_0, \
                                     tly_yyz_xyyz_0, tly_yyz_xyzz_0, tly_yyz_xzzz_0, tly_yyz_yyyy_0, tly_yyz_yyyz_0, \
                                     tly_yz_xxx_0, tly_yz_xxxx_0, tly_yz_xxxy_0, tly_yz_xxxz_0, tly_yz_xxy_0, \
                                     tly_yz_xxyy_0, tly_yz_xxyz_0, tly_yz_xxz_0, tly_yz_xxzz_0, tly_yz_xyy_0, \
                                     tly_yz_xyyy_0, tly_yz_xyyz_0, tly_yz_xyz_0, tly_yz_xyzz_0, tly_yz_xzz_0, \
                                     tly_yz_xzzz_0, tly_yz_yyy_0, tly_yz_yyyy_0, tly_yz_yyyz_0, tly_yz_yyz_0, \
                                     tly_z_xxxx_0, tly_z_xxxy_0, tly_z_xxxz_0, tly_z_xxyy_0, tly_z_xxyz_0, tly_z_xxzz_0, \
                                     tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyzz_0, tly_z_xzzz_0, tly_z_yyyy_0, tly_z_yyyz_0, \
                                     tlz_y_yyyy_0, tlz_y_yyyz_0, tlz_y_yyzz_0, tlz_y_yzzz_0, tlz_y_zzzz_0, tlz_yy_yyy_0, \
                                     tlz_yy_yyyy_0, tlz_yy_yyyz_0, tlz_yy_yyz_0, tlz_yy_yyzz_0, tlz_yy_yzz_0, \
                                     tlz_yy_yzzz_0, tlz_yy_zzz_0, tlz_yy_zzzz_0, tlz_yyy_yyyy_0, tlz_yyy_yyyz_0, \
                                     tlz_yyy_yyzz_0, tlz_yyy_yzzz_0, tlz_yyy_zzzz_0, tlz_yyz_xxxx_0, tlz_yyz_xxxy_0, \
                                     tlz_yyz_xxxz_0, tlz_yyz_xxyy_0, tlz_yyz_xxyz_0, tlz_yyz_xxzz_0, tlz_yyz_xyyy_0, \
                                     tlz_yyz_xyyz_0, tlz_yyz_xyzz_0, tlz_yyz_xzzz_0, tlz_yyz_yyyy_0, tlz_yz_xxx_0, \
                                     tlz_yz_xxxx_0, tlz_yz_xxxy_0, tlz_yz_xxxz_0, tlz_yz_xxy_0, tlz_yz_xxyy_0, \
                                     tlz_yz_xxyz_0, tlz_yz_xxz_0, tlz_yz_xxzz_0, tlz_yz_xyy_0, tlz_yz_xyyy_0, \
                                     tlz_yz_xyyz_0, tlz_yz_xyz_0, tlz_yz_xyzz_0, tlz_yz_xzz_0, tlz_yz_xzzz_0, \
                                     tlz_yz_yyy_0, tlz_yz_yyyy_0, tlz_z_xxxx_0, tlz_z_xxxy_0, tlz_z_xxxz_0, \
                                     tlz_z_xxyy_0, tlz_z_xxyz_0, tlz_z_xxzz_0, tlz_z_xyyy_0, tlz_z_xyyz_0, tlz_z_xyzz_0, \
                                     tlz_z_xzzz_0, tlz_z_yyyy_0, tpx_yy_yyyy_0, tpx_yy_yyyz_0, tpx_yy_yyzz_0, \
                                     tpx_yy_yzzz_0, tpx_yy_zzzz_0, tpx_yz_xxxx_0, tpx_yz_xxxy_0, tpx_yz_xxxz_0, \
                                     tpx_yz_xxyy_0, tpx_yz_xxyz_0, tpx_yz_xxzz_0, tpx_yz_xyyy_0, tpx_yz_xyyz_0, \
                                     tpx_yz_xyzz_0, tpx_yz_xzzz_0, tpx_yz_yyyy_0, tpz_yy_yyyy_0, tpz_yy_yyyz_0, \
                                     tpz_yy_yyzz_0, tpz_yy_yzzz_0, tpz_yy_zzzz_0, tpz_yz_xxxx_0, tpz_yz_xxxy_0, \
                                     tpz_yz_xxxz_0, tpz_yz_xxyy_0, tpz_yz_xxyz_0, tpz_yz_xxzz_0, tpz_yz_xyyy_0, \
                                     tpz_yz_xyyz_0, tpz_yz_xyzz_0, tpz_yz_xzzz_0, tpz_yz_yyyy_0, tpz_yz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyy_yyyy_0[j] = pa_y[j] * tlx_yy_yyyy_0[j] + fl1_fx * tlx_y_yyyy_0[j] + 2.0 * fl1_fx * tlx_yy_yyy_0[j] -
                                0.5 * fl1_fx * tpz_yy_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_yy_yyyy_0[j];

            tly_yyy_yyyy_0[j] = pa_y[j] * tly_yy_yyyy_0[j] + fl1_fx * tly_y_yyyy_0[j] + 2.0 * fl1_fx * tly_yy_yyy_0[j];

            tlz_yyy_yyyy_0[j] = pa_y[j] * tlz_yy_yyyy_0[j] + fl1_fx * tlz_y_yyyy_0[j] + 2.0 * fl1_fx * tlz_yy_yyy_0[j] +
                                0.5 * fl1_fx * tpx_yy_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_yy_yyyy_0[j];

            tlx_yyy_yyyz_0[j] = pa_y[j] * tlx_yy_yyyz_0[j] + fl1_fx * tlx_y_yyyz_0[j] + 1.5 * fl1_fx * tlx_yy_yyz_0[j] -
                                0.5 * fl1_fx * tpz_yy_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_yy_yyyz_0[j];

            tly_yyy_yyyz_0[j] = pa_y[j] * tly_yy_yyyz_0[j] + fl1_fx * tly_y_yyyz_0[j] + 1.5 * fl1_fx * tly_yy_yyz_0[j];

            tlz_yyy_yyyz_0[j] = pa_y[j] * tlz_yy_yyyz_0[j] + fl1_fx * tlz_y_yyyz_0[j] + 1.5 * fl1_fx * tlz_yy_yyz_0[j] +
                                0.5 * fl1_fx * tpx_yy_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_yy_yyyz_0[j];

            tlx_yyy_yyzz_0[j] = pa_y[j] * tlx_yy_yyzz_0[j] + fl1_fx * tlx_y_yyzz_0[j] + fl1_fx * tlx_yy_yzz_0[j] - 0.5 * fl1_fx * tpz_yy_yyzz_0[j] -
                                fl1_fx * fl1_fgb * tdz_yy_yyzz_0[j];

            tly_yyy_yyzz_0[j] = pa_y[j] * tly_yy_yyzz_0[j] + fl1_fx * tly_y_yyzz_0[j] + fl1_fx * tly_yy_yzz_0[j];

            tlz_yyy_yyzz_0[j] = pa_y[j] * tlz_yy_yyzz_0[j] + fl1_fx * tlz_y_yyzz_0[j] + fl1_fx * tlz_yy_yzz_0[j] + 0.5 * fl1_fx * tpx_yy_yyzz_0[j] +
                                fl1_fx * fl1_fgb * tdx_yy_yyzz_0[j];

            tlx_yyy_yzzz_0[j] = pa_y[j] * tlx_yy_yzzz_0[j] + fl1_fx * tlx_y_yzzz_0[j] + 0.5 * fl1_fx * tlx_yy_zzz_0[j] -
                                0.5 * fl1_fx * tpz_yy_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_yzzz_0[j];

            tly_yyy_yzzz_0[j] = pa_y[j] * tly_yy_yzzz_0[j] + fl1_fx * tly_y_yzzz_0[j] + 0.5 * fl1_fx * tly_yy_zzz_0[j];

            tlz_yyy_yzzz_0[j] = pa_y[j] * tlz_yy_yzzz_0[j] + fl1_fx * tlz_y_yzzz_0[j] + 0.5 * fl1_fx * tlz_yy_zzz_0[j] +
                                0.5 * fl1_fx * tpx_yy_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_yzzz_0[j];

            tlx_yyy_zzzz_0[j] =
                pa_y[j] * tlx_yy_zzzz_0[j] + fl1_fx * tlx_y_zzzz_0[j] - 0.5 * fl1_fx * tpz_yy_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_zzzz_0[j];

            tly_yyy_zzzz_0[j] = pa_y[j] * tly_yy_zzzz_0[j] + fl1_fx * tly_y_zzzz_0[j];

            tlz_yyy_zzzz_0[j] =
                pa_y[j] * tlz_yy_zzzz_0[j] + fl1_fx * tlz_y_zzzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_zzzz_0[j];

            tlx_yyz_xxxx_0[j] =
                pa_y[j] * tlx_yz_xxxx_0[j] + 0.5 * fl1_fx * tlx_z_xxxx_0[j] - 0.5 * fl1_fx * tpz_yz_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxxx_0[j];

            tly_yyz_xxxx_0[j] = pa_y[j] * tly_yz_xxxx_0[j] + 0.5 * fl1_fx * tly_z_xxxx_0[j];

            tlz_yyz_xxxx_0[j] =
                pa_y[j] * tlz_yz_xxxx_0[j] + 0.5 * fl1_fx * tlz_z_xxxx_0[j] + 0.5 * fl1_fx * tpx_yz_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxxx_0[j];

            tlx_yyz_xxxy_0[j] = pa_y[j] * tlx_yz_xxxy_0[j] + 0.5 * fl1_fx * tlx_z_xxxy_0[j] + 0.5 * fl1_fx * tlx_yz_xxx_0[j] -
                                0.5 * fl1_fx * tpz_yz_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxxy_0[j];

            tly_yyz_xxxy_0[j] = pa_y[j] * tly_yz_xxxy_0[j] + 0.5 * fl1_fx * tly_z_xxxy_0[j] + 0.5 * fl1_fx * tly_yz_xxx_0[j];

            tlz_yyz_xxxy_0[j] = pa_y[j] * tlz_yz_xxxy_0[j] + 0.5 * fl1_fx * tlz_z_xxxy_0[j] + 0.5 * fl1_fx * tlz_yz_xxx_0[j] +
                                0.5 * fl1_fx * tpx_yz_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxxy_0[j];

            tlx_yyz_xxxz_0[j] =
                pa_y[j] * tlx_yz_xxxz_0[j] + 0.5 * fl1_fx * tlx_z_xxxz_0[j] - 0.5 * fl1_fx * tpz_yz_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxxz_0[j];

            tly_yyz_xxxz_0[j] = pa_y[j] * tly_yz_xxxz_0[j] + 0.5 * fl1_fx * tly_z_xxxz_0[j];

            tlz_yyz_xxxz_0[j] =
                pa_y[j] * tlz_yz_xxxz_0[j] + 0.5 * fl1_fx * tlz_z_xxxz_0[j] + 0.5 * fl1_fx * tpx_yz_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxxz_0[j];

            tlx_yyz_xxyy_0[j] = pa_y[j] * tlx_yz_xxyy_0[j] + 0.5 * fl1_fx * tlx_z_xxyy_0[j] + fl1_fx * tlx_yz_xxy_0[j] -
                                0.5 * fl1_fx * tpz_yz_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxyy_0[j];

            tly_yyz_xxyy_0[j] = pa_y[j] * tly_yz_xxyy_0[j] + 0.5 * fl1_fx * tly_z_xxyy_0[j] + fl1_fx * tly_yz_xxy_0[j];

            tlz_yyz_xxyy_0[j] = pa_y[j] * tlz_yz_xxyy_0[j] + 0.5 * fl1_fx * tlz_z_xxyy_0[j] + fl1_fx * tlz_yz_xxy_0[j] +
                                0.5 * fl1_fx * tpx_yz_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxyy_0[j];

            tlx_yyz_xxyz_0[j] = pa_y[j] * tlx_yz_xxyz_0[j] + 0.5 * fl1_fx * tlx_z_xxyz_0[j] + 0.5 * fl1_fx * tlx_yz_xxz_0[j] -
                                0.5 * fl1_fx * tpz_yz_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxyz_0[j];

            tly_yyz_xxyz_0[j] = pa_y[j] * tly_yz_xxyz_0[j] + 0.5 * fl1_fx * tly_z_xxyz_0[j] + 0.5 * fl1_fx * tly_yz_xxz_0[j];

            tlz_yyz_xxyz_0[j] = pa_y[j] * tlz_yz_xxyz_0[j] + 0.5 * fl1_fx * tlz_z_xxyz_0[j] + 0.5 * fl1_fx * tlz_yz_xxz_0[j] +
                                0.5 * fl1_fx * tpx_yz_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxyz_0[j];

            tlx_yyz_xxzz_0[j] =
                pa_y[j] * tlx_yz_xxzz_0[j] + 0.5 * fl1_fx * tlx_z_xxzz_0[j] - 0.5 * fl1_fx * tpz_yz_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxzz_0[j];

            tly_yyz_xxzz_0[j] = pa_y[j] * tly_yz_xxzz_0[j] + 0.5 * fl1_fx * tly_z_xxzz_0[j];

            tlz_yyz_xxzz_0[j] =
                pa_y[j] * tlz_yz_xxzz_0[j] + 0.5 * fl1_fx * tlz_z_xxzz_0[j] + 0.5 * fl1_fx * tpx_yz_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxzz_0[j];

            tlx_yyz_xyyy_0[j] = pa_y[j] * tlx_yz_xyyy_0[j] + 0.5 * fl1_fx * tlx_z_xyyy_0[j] + 1.5 * fl1_fx * tlx_yz_xyy_0[j] -
                                0.5 * fl1_fx * tpz_yz_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_yz_xyyy_0[j];

            tly_yyz_xyyy_0[j] = pa_y[j] * tly_yz_xyyy_0[j] + 0.5 * fl1_fx * tly_z_xyyy_0[j] + 1.5 * fl1_fx * tly_yz_xyy_0[j];

            tlz_yyz_xyyy_0[j] = pa_y[j] * tlz_yz_xyyy_0[j] + 0.5 * fl1_fx * tlz_z_xyyy_0[j] + 1.5 * fl1_fx * tlz_yz_xyy_0[j] +
                                0.5 * fl1_fx * tpx_yz_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_yz_xyyy_0[j];

            tlx_yyz_xyyz_0[j] = pa_y[j] * tlx_yz_xyyz_0[j] + 0.5 * fl1_fx * tlx_z_xyyz_0[j] + fl1_fx * tlx_yz_xyz_0[j] -
                                0.5 * fl1_fx * tpz_yz_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xyyz_0[j];

            tly_yyz_xyyz_0[j] = pa_y[j] * tly_yz_xyyz_0[j] + 0.5 * fl1_fx * tly_z_xyyz_0[j] + fl1_fx * tly_yz_xyz_0[j];

            tlz_yyz_xyyz_0[j] = pa_y[j] * tlz_yz_xyyz_0[j] + 0.5 * fl1_fx * tlz_z_xyyz_0[j] + fl1_fx * tlz_yz_xyz_0[j] +
                                0.5 * fl1_fx * tpx_yz_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xyyz_0[j];

            tlx_yyz_xyzz_0[j] = pa_y[j] * tlx_yz_xyzz_0[j] + 0.5 * fl1_fx * tlx_z_xyzz_0[j] + 0.5 * fl1_fx * tlx_yz_xzz_0[j] -
                                0.5 * fl1_fx * tpz_yz_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xyzz_0[j];

            tly_yyz_xyzz_0[j] = pa_y[j] * tly_yz_xyzz_0[j] + 0.5 * fl1_fx * tly_z_xyzz_0[j] + 0.5 * fl1_fx * tly_yz_xzz_0[j];

            tlz_yyz_xyzz_0[j] = pa_y[j] * tlz_yz_xyzz_0[j] + 0.5 * fl1_fx * tlz_z_xyzz_0[j] + 0.5 * fl1_fx * tlz_yz_xzz_0[j] +
                                0.5 * fl1_fx * tpx_yz_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xyzz_0[j];

            tlx_yyz_xzzz_0[j] =
                pa_y[j] * tlx_yz_xzzz_0[j] + 0.5 * fl1_fx * tlx_z_xzzz_0[j] - 0.5 * fl1_fx * tpz_yz_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xzzz_0[j];

            tly_yyz_xzzz_0[j] = pa_y[j] * tly_yz_xzzz_0[j] + 0.5 * fl1_fx * tly_z_xzzz_0[j];

            tlz_yyz_xzzz_0[j] =
                pa_y[j] * tlz_yz_xzzz_0[j] + 0.5 * fl1_fx * tlz_z_xzzz_0[j] + 0.5 * fl1_fx * tpx_yz_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xzzz_0[j];

            tlx_yyz_yyyy_0[j] = pa_y[j] * tlx_yz_yyyy_0[j] + 0.5 * fl1_fx * tlx_z_yyyy_0[j] + 2.0 * fl1_fx * tlx_yz_yyy_0[j] -
                                0.5 * fl1_fx * tpz_yz_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_yz_yyyy_0[j];

            tly_yyz_yyyy_0[j] = pa_y[j] * tly_yz_yyyy_0[j] + 0.5 * fl1_fx * tly_z_yyyy_0[j] + 2.0 * fl1_fx * tly_yz_yyy_0[j];

            tlz_yyz_yyyy_0[j] = pa_y[j] * tlz_yz_yyyy_0[j] + 0.5 * fl1_fx * tlz_z_yyyy_0[j] + 2.0 * fl1_fx * tlz_yz_yyy_0[j] +
                                0.5 * fl1_fx * tpx_yz_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_yz_yyyy_0[j];

            tlx_yyz_yyyz_0[j] = pa_y[j] * tlx_yz_yyyz_0[j] + 0.5 * fl1_fx * tlx_z_yyyz_0[j] + 1.5 * fl1_fx * tlx_yz_yyz_0[j] -
                                0.5 * fl1_fx * tpz_yz_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_yz_yyyz_0[j];

            tly_yyz_yyyz_0[j] = pa_y[j] * tly_yz_yyyz_0[j] + 0.5 * fl1_fx * tly_z_yyyz_0[j] + 1.5 * fl1_fx * tly_yz_yyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_350_400(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (350,400)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlz_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tlx_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 72);

        auto tly_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tlz_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tlx_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 73);

        auto tly_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tlz_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tlx_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 74);

        auto tly_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tlz_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

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

        auto tlz_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tlx_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 48);

        auto tly_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tlz_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tlx_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 49);

        auto tly_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tlz_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tlx_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 50);

        auto tly_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tlz_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tlx_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 51);

        auto tly_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tlz_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tlx_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 52);

        auto tly_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tlz_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tlx_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 53);

        auto tly_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tlz_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tlx_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 54);

        auto tly_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 55);

        auto tly_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tlz_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tlx_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 56);

        auto tly_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tlz_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tlx_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 57);

        auto tly_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tlz_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tlx_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 58);

        auto tly_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 59);

        auto tpx_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 71);

        auto tpx_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 72);

        auto tpz_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tpx_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 73);

        auto tpz_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tpx_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 74);

        auto tpz_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tdx_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 71);

        auto tdx_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 72);

        auto tdz_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tdx_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 73);

        auto tdz_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tdx_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 74);

        auto tdz_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75);

        auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76);

        auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77);

        auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78);

        auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79);

        auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80);

        auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81);

        auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82);

        auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83);

        auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84);

        auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85);

        auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86);

        auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87);

        auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88);

        // set up pointers to integrals

        auto tlz_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tlx_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 117);

        auto tly_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tlz_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tlx_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 118);

        auto tly_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tlz_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tlx_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 119);

        auto tly_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tlz_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tlx_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 120);

        auto tly_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tlz_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tlx_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 121);

        auto tly_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tlz_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tlx_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 122);

        auto tly_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tlz_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tlx_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 123);

        auto tly_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tlz_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tlx_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 124);

        auto tly_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tlz_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tlx_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 125);

        auto tly_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tlz_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tlx_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 126);

        auto tly_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tlz_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tlx_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 127);

        auto tly_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tlz_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tlx_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 128);

        auto tly_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tlz_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tlx_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 129);

        auto tly_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tlz_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tlx_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 130);

        auto tly_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tlz_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tlx_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 131);

        auto tly_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tlz_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tlx_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 132);

        auto tly_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tlz_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tlx_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 133);

        // Batch of Integrals (350,400)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yz_yyyz_0, tdx_yz_yyzz_0, tdx_yz_yzzz_0, \
                                     tdx_yz_zzzz_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, tdx_zz_xxyy_0, \
                                     tdx_zz_xxyz_0, tdx_zz_xxzz_0, tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyzz_0, \
                                     tdx_zz_xzzz_0, tdx_zz_yyyy_0, tdx_zz_yyyz_0, tdx_zz_yyzz_0, tdz_yz_yyzz_0, \
                                     tdz_yz_yzzz_0, tdz_yz_zzzz_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, \
                                     tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxzz_0, tdz_zz_xyyy_0, tdz_zz_xyyz_0, \
                                     tdz_zz_xyzz_0, tdz_zz_xzzz_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyzz_0, \
                                     tdz_zz_yzzz_0, tlx_yyz_yyzz_0, tlx_yyz_yzzz_0, tlx_yyz_zzzz_0, tlx_yz_yyzz_0, \
                                     tlx_yz_yzz_0, tlx_yz_yzzz_0, tlx_yz_zzz_0, tlx_yz_zzzz_0, tlx_yzz_xxxx_0, \
                                     tlx_yzz_xxxy_0, tlx_yzz_xxxz_0, tlx_yzz_xxyy_0, tlx_yzz_xxyz_0, tlx_yzz_xxzz_0, \
                                     tlx_yzz_xyyy_0, tlx_yzz_xyyz_0, tlx_yzz_xyzz_0, tlx_yzz_xzzz_0, tlx_yzz_yyyy_0, \
                                     tlx_yzz_yyyz_0, tlx_yzz_yyzz_0, tlx_yzz_yzzz_0, tlx_z_yyzz_0, tlx_z_yzzz_0, \
                                     tlx_z_zzzz_0, tlx_zz_xxx_0, tlx_zz_xxxx_0, tlx_zz_xxxy_0, tlx_zz_xxxz_0, \
                                     tlx_zz_xxy_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, tlx_zz_xxz_0, tlx_zz_xxzz_0, \
                                     tlx_zz_xyy_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyz_0, tlx_zz_xyzz_0, \
                                     tlx_zz_xzz_0, tlx_zz_xzzz_0, tlx_zz_yyy_0, tlx_zz_yyyy_0, tlx_zz_yyyz_0, \
                                     tlx_zz_yyz_0, tlx_zz_yyzz_0, tlx_zz_yzz_0, tlx_zz_yzzz_0, tlx_zz_zzz_0, \
                                     tly_yyz_yyzz_0, tly_yyz_yzzz_0, tly_yyz_zzzz_0, tly_yz_yyzz_0, tly_yz_yzz_0, \
                                     tly_yz_yzzz_0, tly_yz_zzz_0, tly_yz_zzzz_0, tly_yzz_xxxx_0, tly_yzz_xxxy_0, \
                                     tly_yzz_xxxz_0, tly_yzz_xxyy_0, tly_yzz_xxyz_0, tly_yzz_xxzz_0, tly_yzz_xyyy_0, \
                                     tly_yzz_xyyz_0, tly_yzz_xyzz_0, tly_yzz_xzzz_0, tly_yzz_yyyy_0, tly_yzz_yyyz_0, \
                                     tly_yzz_yyzz_0, tly_z_yyzz_0, tly_z_yzzz_0, tly_z_zzzz_0, tly_zz_xxx_0, \
                                     tly_zz_xxxx_0, tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxy_0, tly_zz_xxyy_0, \
                                     tly_zz_xxyz_0, tly_zz_xxz_0, tly_zz_xxzz_0, tly_zz_xyy_0, tly_zz_xyyy_0, \
                                     tly_zz_xyyz_0, tly_zz_xyz_0, tly_zz_xyzz_0, tly_zz_xzz_0, tly_zz_xzzz_0, \
                                     tly_zz_yyy_0, tly_zz_yyyy_0, tly_zz_yyyz_0, tly_zz_yyz_0, tly_zz_yyzz_0, \
                                     tly_zz_yzz_0, tlz_yyz_yyyz_0, tlz_yyz_yyzz_0, tlz_yyz_yzzz_0, tlz_yyz_zzzz_0, \
                                     tlz_yz_yyyz_0, tlz_yz_yyz_0, tlz_yz_yyzz_0, tlz_yz_yzz_0, tlz_yz_yzzz_0, \
                                     tlz_yz_zzz_0, tlz_yz_zzzz_0, tlz_yzz_xxxx_0, tlz_yzz_xxxy_0, tlz_yzz_xxxz_0, \
                                     tlz_yzz_xxyy_0, tlz_yzz_xxyz_0, tlz_yzz_xxzz_0, tlz_yzz_xyyy_0, tlz_yzz_xyyz_0, \
                                     tlz_yzz_xyzz_0, tlz_yzz_xzzz_0, tlz_yzz_yyyy_0, tlz_yzz_yyyz_0, tlz_yzz_yyzz_0, \
                                     tlz_z_yyyz_0, tlz_z_yyzz_0, tlz_z_yzzz_0, tlz_z_zzzz_0, tlz_zz_xxx_0, \
                                     tlz_zz_xxxx_0, tlz_zz_xxxy_0, tlz_zz_xxxz_0, tlz_zz_xxy_0, tlz_zz_xxyy_0, \
                                     tlz_zz_xxyz_0, tlz_zz_xxz_0, tlz_zz_xxzz_0, tlz_zz_xyy_0, tlz_zz_xyyy_0, \
                                     tlz_zz_xyyz_0, tlz_zz_xyz_0, tlz_zz_xyzz_0, tlz_zz_xzz_0, tlz_zz_xzzz_0, \
                                     tlz_zz_yyy_0, tlz_zz_yyyy_0, tlz_zz_yyyz_0, tlz_zz_yyz_0, tlz_zz_yyzz_0, \
                                     tlz_zz_yzz_0, tpx_yz_yyyz_0, tpx_yz_yyzz_0, tpx_yz_yzzz_0, tpx_yz_zzzz_0, \
                                     tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, \
                                     tpx_zz_xxzz_0, tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyzz_0, tpx_zz_xzzz_0, \
                                     tpx_zz_yyyy_0, tpx_zz_yyyz_0, tpx_zz_yyzz_0, tpz_yz_yyzz_0, tpz_yz_yzzz_0, \
                                     tpz_yz_zzzz_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, tpz_zz_xxyy_0, \
                                     tpz_zz_xxyz_0, tpz_zz_xxzz_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0, tpz_zz_xyzz_0, \
                                     tpz_zz_xzzz_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, tpz_zz_yyzz_0, tpz_zz_yzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_yyz_yyyz_0[j] = pa_y[j] * tlz_yz_yyyz_0[j] + 0.5 * fl1_fx * tlz_z_yyyz_0[j] + 1.5 * fl1_fx * tlz_yz_yyz_0[j] +
                                0.5 * fl1_fx * tpx_yz_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_yz_yyyz_0[j];

            tlx_yyz_yyzz_0[j] = pa_y[j] * tlx_yz_yyzz_0[j] + 0.5 * fl1_fx * tlx_z_yyzz_0[j] + fl1_fx * tlx_yz_yzz_0[j] -
                                0.5 * fl1_fx * tpz_yz_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_yyzz_0[j];

            tly_yyz_yyzz_0[j] = pa_y[j] * tly_yz_yyzz_0[j] + 0.5 * fl1_fx * tly_z_yyzz_0[j] + fl1_fx * tly_yz_yzz_0[j];

            tlz_yyz_yyzz_0[j] = pa_y[j] * tlz_yz_yyzz_0[j] + 0.5 * fl1_fx * tlz_z_yyzz_0[j] + fl1_fx * tlz_yz_yzz_0[j] +
                                0.5 * fl1_fx * tpx_yz_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_yyzz_0[j];

            tlx_yyz_yzzz_0[j] = pa_y[j] * tlx_yz_yzzz_0[j] + 0.5 * fl1_fx * tlx_z_yzzz_0[j] + 0.5 * fl1_fx * tlx_yz_zzz_0[j] -
                                0.5 * fl1_fx * tpz_yz_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_yzzz_0[j];

            tly_yyz_yzzz_0[j] = pa_y[j] * tly_yz_yzzz_0[j] + 0.5 * fl1_fx * tly_z_yzzz_0[j] + 0.5 * fl1_fx * tly_yz_zzz_0[j];

            tlz_yyz_yzzz_0[j] = pa_y[j] * tlz_yz_yzzz_0[j] + 0.5 * fl1_fx * tlz_z_yzzz_0[j] + 0.5 * fl1_fx * tlz_yz_zzz_0[j] +
                                0.5 * fl1_fx * tpx_yz_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_yzzz_0[j];

            tlx_yyz_zzzz_0[j] =
                pa_y[j] * tlx_yz_zzzz_0[j] + 0.5 * fl1_fx * tlx_z_zzzz_0[j] - 0.5 * fl1_fx * tpz_yz_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_zzzz_0[j];

            tly_yyz_zzzz_0[j] = pa_y[j] * tly_yz_zzzz_0[j] + 0.5 * fl1_fx * tly_z_zzzz_0[j];

            tlz_yyz_zzzz_0[j] =
                pa_y[j] * tlz_yz_zzzz_0[j] + 0.5 * fl1_fx * tlz_z_zzzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_zzzz_0[j];

            tlx_yzz_xxxx_0[j] = pa_y[j] * tlx_zz_xxxx_0[j] - 0.5 * fl1_fx * tpz_zz_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxxx_0[j];

            tly_yzz_xxxx_0[j] = pa_y[j] * tly_zz_xxxx_0[j];

            tlz_yzz_xxxx_0[j] = pa_y[j] * tlz_zz_xxxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxxx_0[j];

            tlx_yzz_xxxy_0[j] =
                pa_y[j] * tlx_zz_xxxy_0[j] + 0.5 * fl1_fx * tlx_zz_xxx_0[j] - 0.5 * fl1_fx * tpz_zz_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxxy_0[j];

            tly_yzz_xxxy_0[j] = pa_y[j] * tly_zz_xxxy_0[j] + 0.5 * fl1_fx * tly_zz_xxx_0[j];

            tlz_yzz_xxxy_0[j] =
                pa_y[j] * tlz_zz_xxxy_0[j] + 0.5 * fl1_fx * tlz_zz_xxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxxy_0[j];

            tlx_yzz_xxxz_0[j] = pa_y[j] * tlx_zz_xxxz_0[j] - 0.5 * fl1_fx * tpz_zz_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxxz_0[j];

            tly_yzz_xxxz_0[j] = pa_y[j] * tly_zz_xxxz_0[j];

            tlz_yzz_xxxz_0[j] = pa_y[j] * tlz_zz_xxxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxxz_0[j];

            tlx_yzz_xxyy_0[j] =
                pa_y[j] * tlx_zz_xxyy_0[j] + fl1_fx * tlx_zz_xxy_0[j] - 0.5 * fl1_fx * tpz_zz_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxyy_0[j];

            tly_yzz_xxyy_0[j] = pa_y[j] * tly_zz_xxyy_0[j] + fl1_fx * tly_zz_xxy_0[j];

            tlz_yzz_xxyy_0[j] =
                pa_y[j] * tlz_zz_xxyy_0[j] + fl1_fx * tlz_zz_xxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxyy_0[j];

            tlx_yzz_xxyz_0[j] =
                pa_y[j] * tlx_zz_xxyz_0[j] + 0.5 * fl1_fx * tlx_zz_xxz_0[j] - 0.5 * fl1_fx * tpz_zz_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxyz_0[j];

            tly_yzz_xxyz_0[j] = pa_y[j] * tly_zz_xxyz_0[j] + 0.5 * fl1_fx * tly_zz_xxz_0[j];

            tlz_yzz_xxyz_0[j] =
                pa_y[j] * tlz_zz_xxyz_0[j] + 0.5 * fl1_fx * tlz_zz_xxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxyz_0[j];

            tlx_yzz_xxzz_0[j] = pa_y[j] * tlx_zz_xxzz_0[j] - 0.5 * fl1_fx * tpz_zz_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxzz_0[j];

            tly_yzz_xxzz_0[j] = pa_y[j] * tly_zz_xxzz_0[j];

            tlz_yzz_xxzz_0[j] = pa_y[j] * tlz_zz_xxzz_0[j] + 0.5 * fl1_fx * tpx_zz_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxzz_0[j];

            tlx_yzz_xyyy_0[j] =
                pa_y[j] * tlx_zz_xyyy_0[j] + 1.5 * fl1_fx * tlx_zz_xyy_0[j] - 0.5 * fl1_fx * tpz_zz_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xyyy_0[j];

            tly_yzz_xyyy_0[j] = pa_y[j] * tly_zz_xyyy_0[j] + 1.5 * fl1_fx * tly_zz_xyy_0[j];

            tlz_yzz_xyyy_0[j] =
                pa_y[j] * tlz_zz_xyyy_0[j] + 1.5 * fl1_fx * tlz_zz_xyy_0[j] + 0.5 * fl1_fx * tpx_zz_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xyyy_0[j];

            tlx_yzz_xyyz_0[j] =
                pa_y[j] * tlx_zz_xyyz_0[j] + fl1_fx * tlx_zz_xyz_0[j] - 0.5 * fl1_fx * tpz_zz_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xyyz_0[j];

            tly_yzz_xyyz_0[j] = pa_y[j] * tly_zz_xyyz_0[j] + fl1_fx * tly_zz_xyz_0[j];

            tlz_yzz_xyyz_0[j] =
                pa_y[j] * tlz_zz_xyyz_0[j] + fl1_fx * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xyyz_0[j];

            tlx_yzz_xyzz_0[j] =
                pa_y[j] * tlx_zz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zz_xzz_0[j] - 0.5 * fl1_fx * tpz_zz_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xyzz_0[j];

            tly_yzz_xyzz_0[j] = pa_y[j] * tly_zz_xyzz_0[j] + 0.5 * fl1_fx * tly_zz_xzz_0[j];

            tlz_yzz_xyzz_0[j] =
                pa_y[j] * tlz_zz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zz_xzz_0[j] + 0.5 * fl1_fx * tpx_zz_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xyzz_0[j];

            tlx_yzz_xzzz_0[j] = pa_y[j] * tlx_zz_xzzz_0[j] - 0.5 * fl1_fx * tpz_zz_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xzzz_0[j];

            tly_yzz_xzzz_0[j] = pa_y[j] * tly_zz_xzzz_0[j];

            tlz_yzz_xzzz_0[j] = pa_y[j] * tlz_zz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xzzz_0[j];

            tlx_yzz_yyyy_0[j] =
                pa_y[j] * tlx_zz_yyyy_0[j] + 2.0 * fl1_fx * tlx_zz_yyy_0[j] - 0.5 * fl1_fx * tpz_zz_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_zz_yyyy_0[j];

            tly_yzz_yyyy_0[j] = pa_y[j] * tly_zz_yyyy_0[j] + 2.0 * fl1_fx * tly_zz_yyy_0[j];

            tlz_yzz_yyyy_0[j] =
                pa_y[j] * tlz_zz_yyyy_0[j] + 2.0 * fl1_fx * tlz_zz_yyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_zz_yyyy_0[j];

            tlx_yzz_yyyz_0[j] =
                pa_y[j] * tlx_zz_yyyz_0[j] + 1.5 * fl1_fx * tlx_zz_yyz_0[j] - 0.5 * fl1_fx * tpz_zz_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yyyz_0[j];

            tly_yzz_yyyz_0[j] = pa_y[j] * tly_zz_yyyz_0[j] + 1.5 * fl1_fx * tly_zz_yyz_0[j];

            tlz_yzz_yyyz_0[j] =
                pa_y[j] * tlz_zz_yyyz_0[j] + 1.5 * fl1_fx * tlz_zz_yyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yyyz_0[j];

            tlx_yzz_yyzz_0[j] =
                pa_y[j] * tlx_zz_yyzz_0[j] + fl1_fx * tlx_zz_yzz_0[j] - 0.5 * fl1_fx * tpz_zz_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yyzz_0[j];

            tly_yzz_yyzz_0[j] = pa_y[j] * tly_zz_yyzz_0[j] + fl1_fx * tly_zz_yzz_0[j];

            tlz_yzz_yyzz_0[j] =
                pa_y[j] * tlz_zz_yyzz_0[j] + fl1_fx * tlz_zz_yzz_0[j] + 0.5 * fl1_fx * tpx_zz_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yyzz_0[j];

            tlx_yzz_yzzz_0[j] =
                pa_y[j] * tlx_zz_yzzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzz_0[j] - 0.5 * fl1_fx * tpz_zz_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFG_400_450(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (400,450)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

        auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89);

        auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89);

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

        auto tlx_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 50);

        auto tly_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tlz_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tlx_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 51);

        auto tly_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tlz_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tlx_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 52);

        auto tly_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tlz_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tlx_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 53);

        auto tly_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tlz_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tlx_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 54);

        auto tly_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 55);

        auto tly_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tlz_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tlx_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 56);

        auto tly_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tlz_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tlx_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 57);

        auto tly_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tlz_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tlx_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 58);

        auto tly_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 59);

        auto tly_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tlz_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpx_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 88);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpx_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 89);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75);

        auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76);

        auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77);

        auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78);

        auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79);

        auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80);

        auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81);

        auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82);

        auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83);

        auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84);

        auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85);

        auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86);

        auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87);

        auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88);

        auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tdx_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 89);

        auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89);

        // set up pointers to integrals

        auto tly_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tlz_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tlx_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 134);

        auto tly_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tlz_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tlx_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 135);

        auto tly_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tlz_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tlx_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 136);

        auto tly_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tlz_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tlx_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 137);

        auto tly_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tlz_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tlx_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 138);

        auto tly_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tlz_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tlx_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 139);

        auto tly_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tlz_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tlx_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 140);

        auto tly_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tlz_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tlx_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 141);

        auto tly_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tlz_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tlx_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 142);

        auto tly_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tlz_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tlx_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 143);

        auto tly_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tlz_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tlx_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 144);

        auto tly_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tlz_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tlx_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 145);

        auto tly_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tlz_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tlx_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 146);

        auto tly_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tlz_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tlx_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 147);

        auto tly_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tlz_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tlx_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 148);

        auto tly_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tlz_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tlx_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 149);

        auto tly_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tlz_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 149);

        // Batch of Integrals (400,450)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, \
                                     tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxzz_0, tdx_zz_xyyy_0, tdx_zz_xyyz_0, \
                                     tdx_zz_xyzz_0, tdx_zz_xzzz_0, tdx_zz_yyyy_0, tdx_zz_yyyz_0, tdx_zz_yyzz_0, \
                                     tdx_zz_yzzz_0, tdx_zz_zzzz_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, \
                                     tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxzz_0, tdy_zz_xyyy_0, tdy_zz_xyyz_0, \
                                     tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, tdy_zz_yyyz_0, tdy_zz_yyzz_0, \
                                     tdy_zz_yzzz_0, tdy_zz_zzzz_0, tdz_zz_zzzz_0, tlx_yzz_zzzz_0, tlx_z_xxxx_0, \
                                     tlx_z_xxxy_0, tlx_z_xxxz_0, tlx_z_xxyy_0, tlx_z_xxyz_0, tlx_z_xxzz_0, tlx_z_xyyy_0, \
                                     tlx_z_xyyz_0, tlx_z_xyzz_0, tlx_z_xzzz_0, tlx_z_yyyy_0, tlx_z_yyyz_0, tlx_z_yyzz_0, \
                                     tlx_z_yzzz_0, tlx_z_zzzz_0, tlx_zz_xxx_0, tlx_zz_xxxx_0, tlx_zz_xxxy_0, \
                                     tlx_zz_xxxz_0, tlx_zz_xxy_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, tlx_zz_xxz_0, \
                                     tlx_zz_xxzz_0, tlx_zz_xyy_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyz_0, \
                                     tlx_zz_xyzz_0, tlx_zz_xzz_0, tlx_zz_xzzz_0, tlx_zz_yyy_0, tlx_zz_yyyy_0, \
                                     tlx_zz_yyyz_0, tlx_zz_yyz_0, tlx_zz_yyzz_0, tlx_zz_yzz_0, tlx_zz_yzzz_0, \
                                     tlx_zz_zzz_0, tlx_zz_zzzz_0, tlx_zzz_xxxx_0, tlx_zzz_xxxy_0, tlx_zzz_xxxz_0, \
                                     tlx_zzz_xxyy_0, tlx_zzz_xxyz_0, tlx_zzz_xxzz_0, tlx_zzz_xyyy_0, tlx_zzz_xyyz_0, \
                                     tlx_zzz_xyzz_0, tlx_zzz_xzzz_0, tlx_zzz_yyyy_0, tlx_zzz_yyyz_0, tlx_zzz_yyzz_0, \
                                     tlx_zzz_yzzz_0, tlx_zzz_zzzz_0, tly_yzz_yzzz_0, tly_yzz_zzzz_0, tly_z_xxxx_0, \
                                     tly_z_xxxy_0, tly_z_xxxz_0, tly_z_xxyy_0, tly_z_xxyz_0, tly_z_xxzz_0, tly_z_xyyy_0, \
                                     tly_z_xyyz_0, tly_z_xyzz_0, tly_z_xzzz_0, tly_z_yyyy_0, tly_z_yyyz_0, tly_z_yyzz_0, \
                                     tly_z_yzzz_0, tly_z_zzzz_0, tly_zz_xxx_0, tly_zz_xxxx_0, tly_zz_xxxy_0, \
                                     tly_zz_xxxz_0, tly_zz_xxy_0, tly_zz_xxyy_0, tly_zz_xxyz_0, tly_zz_xxz_0, \
                                     tly_zz_xxzz_0, tly_zz_xyy_0, tly_zz_xyyy_0, tly_zz_xyyz_0, tly_zz_xyz_0, \
                                     tly_zz_xyzz_0, tly_zz_xzz_0, tly_zz_xzzz_0, tly_zz_yyy_0, tly_zz_yyyy_0, \
                                     tly_zz_yyyz_0, tly_zz_yyz_0, tly_zz_yyzz_0, tly_zz_yzz_0, tly_zz_yzzz_0, \
                                     tly_zz_zzz_0, tly_zz_zzzz_0, tly_zzz_xxxx_0, tly_zzz_xxxy_0, tly_zzz_xxxz_0, \
                                     tly_zzz_xxyy_0, tly_zzz_xxyz_0, tly_zzz_xxzz_0, tly_zzz_xyyy_0, tly_zzz_xyyz_0, \
                                     tly_zzz_xyzz_0, tly_zzz_xzzz_0, tly_zzz_yyyy_0, tly_zzz_yyyz_0, tly_zzz_yyzz_0, \
                                     tly_zzz_yzzz_0, tly_zzz_zzzz_0, tlz_yzz_yzzz_0, tlz_yzz_zzzz_0, tlz_z_xxxx_0, \
                                     tlz_z_xxxy_0, tlz_z_xxxz_0, tlz_z_xxyy_0, tlz_z_xxyz_0, tlz_z_xxzz_0, tlz_z_xyyy_0, \
                                     tlz_z_xyyz_0, tlz_z_xyzz_0, tlz_z_xzzz_0, tlz_z_yyyy_0, tlz_z_yyyz_0, tlz_z_yyzz_0, \
                                     tlz_z_yzzz_0, tlz_z_zzzz_0, tlz_zz_xxx_0, tlz_zz_xxxx_0, tlz_zz_xxxy_0, \
                                     tlz_zz_xxxz_0, tlz_zz_xxy_0, tlz_zz_xxyy_0, tlz_zz_xxyz_0, tlz_zz_xxz_0, \
                                     tlz_zz_xxzz_0, tlz_zz_xyy_0, tlz_zz_xyyy_0, tlz_zz_xyyz_0, tlz_zz_xyz_0, \
                                     tlz_zz_xyzz_0, tlz_zz_xzz_0, tlz_zz_xzzz_0, tlz_zz_yyy_0, tlz_zz_yyyy_0, \
                                     tlz_zz_yyyz_0, tlz_zz_yyz_0, tlz_zz_yyzz_0, tlz_zz_yzz_0, tlz_zz_yzzz_0, \
                                     tlz_zz_zzz_0, tlz_zz_zzzz_0, tlz_zzz_xxxx_0, tlz_zzz_xxxy_0, tlz_zzz_xxxz_0, \
                                     tlz_zzz_xxyy_0, tlz_zzz_xxyz_0, tlz_zzz_xxzz_0, tlz_zzz_xyyy_0, tlz_zzz_xyyz_0, \
                                     tlz_zzz_xyzz_0, tlz_zzz_xzzz_0, tlz_zzz_yyyy_0, tlz_zzz_yyyz_0, tlz_zzz_yyzz_0, \
                                     tlz_zzz_yzzz_0, tlz_zzz_zzzz_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, \
                                     tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxzz_0, tpx_zz_xyyy_0, tpx_zz_xyyz_0, \
                                     tpx_zz_xyzz_0, tpx_zz_xzzz_0, tpx_zz_yyyy_0, tpx_zz_yyyz_0, tpx_zz_yyzz_0, \
                                     tpx_zz_yzzz_0, tpx_zz_zzzz_0, tpy_zz_xxxx_0, tpy_zz_xxxy_0, tpy_zz_xxxz_0, \
                                     tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, tpy_zz_xyyy_0, tpy_zz_xyyz_0, \
                                     tpy_zz_xyzz_0, tpy_zz_xzzz_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, tpy_zz_yyzz_0, \
                                     tpy_zz_yzzz_0, tpy_zz_zzzz_0, tpz_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_yzz_yzzz_0[j] = pa_y[j] * tly_zz_yzzz_0[j] + 0.5 * fl1_fx * tly_zz_zzz_0[j];

            tlz_yzz_yzzz_0[j] =
                pa_y[j] * tlz_zz_yzzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yzzz_0[j];

            tlx_yzz_zzzz_0[j] = pa_y[j] * tlx_zz_zzzz_0[j] - 0.5 * fl1_fx * tpz_zz_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_zzzz_0[j];

            tly_yzz_zzzz_0[j] = pa_y[j] * tly_zz_zzzz_0[j];

            tlz_yzz_zzzz_0[j] = pa_y[j] * tlz_zz_zzzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_zzzz_0[j];

            tlx_zzz_xxxx_0[j] =
                pa_z[j] * tlx_zz_xxxx_0[j] + fl1_fx * tlx_z_xxxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxxx_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxxx_0[j];

            tly_zzz_xxxx_0[j] =
                pa_z[j] * tly_zz_xxxx_0[j] + fl1_fx * tly_z_xxxx_0[j] - 0.5 * fl1_fx * tpx_zz_xxxx_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxxx_0[j];

            tlz_zzz_xxxx_0[j] = pa_z[j] * tlz_zz_xxxx_0[j] + fl1_fx * tlz_z_xxxx_0[j];

            tlx_zzz_xxxy_0[j] =
                pa_z[j] * tlx_zz_xxxy_0[j] + fl1_fx * tlx_z_xxxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxxy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxxy_0[j];

            tly_zzz_xxxy_0[j] =
                pa_z[j] * tly_zz_xxxy_0[j] + fl1_fx * tly_z_xxxy_0[j] - 0.5 * fl1_fx * tpx_zz_xxxy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxxy_0[j];

            tlz_zzz_xxxy_0[j] = pa_z[j] * tlz_zz_xxxy_0[j] + fl1_fx * tlz_z_xxxy_0[j];

            tlx_zzz_xxxz_0[j] = pa_z[j] * tlx_zz_xxxz_0[j] + fl1_fx * tlx_z_xxxz_0[j] + 0.5 * fl1_fx * tlx_zz_xxx_0[j] +
                                0.5 * fl1_fx * tpy_zz_xxxz_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxxz_0[j];

            tly_zzz_xxxz_0[j] = pa_z[j] * tly_zz_xxxz_0[j] + fl1_fx * tly_z_xxxz_0[j] + 0.5 * fl1_fx * tly_zz_xxx_0[j] -
                                0.5 * fl1_fx * tpx_zz_xxxz_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxxz_0[j];

            tlz_zzz_xxxz_0[j] = pa_z[j] * tlz_zz_xxxz_0[j] + fl1_fx * tlz_z_xxxz_0[j] + 0.5 * fl1_fx * tlz_zz_xxx_0[j];

            tlx_zzz_xxyy_0[j] =
                pa_z[j] * tlx_zz_xxyy_0[j] + fl1_fx * tlx_z_xxyy_0[j] + 0.5 * fl1_fx * tpy_zz_xxyy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxyy_0[j];

            tly_zzz_xxyy_0[j] =
                pa_z[j] * tly_zz_xxyy_0[j] + fl1_fx * tly_z_xxyy_0[j] - 0.5 * fl1_fx * tpx_zz_xxyy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxyy_0[j];

            tlz_zzz_xxyy_0[j] = pa_z[j] * tlz_zz_xxyy_0[j] + fl1_fx * tlz_z_xxyy_0[j];

            tlx_zzz_xxyz_0[j] = pa_z[j] * tlx_zz_xxyz_0[j] + fl1_fx * tlx_z_xxyz_0[j] + 0.5 * fl1_fx * tlx_zz_xxy_0[j] +
                                0.5 * fl1_fx * tpy_zz_xxyz_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxyz_0[j];

            tly_zzz_xxyz_0[j] = pa_z[j] * tly_zz_xxyz_0[j] + fl1_fx * tly_z_xxyz_0[j] + 0.5 * fl1_fx * tly_zz_xxy_0[j] -
                                0.5 * fl1_fx * tpx_zz_xxyz_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxyz_0[j];

            tlz_zzz_xxyz_0[j] = pa_z[j] * tlz_zz_xxyz_0[j] + fl1_fx * tlz_z_xxyz_0[j] + 0.5 * fl1_fx * tlz_zz_xxy_0[j];

            tlx_zzz_xxzz_0[j] = pa_z[j] * tlx_zz_xxzz_0[j] + fl1_fx * tlx_z_xxzz_0[j] + fl1_fx * tlx_zz_xxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxzz_0[j] +
                                fl1_fx * fl1_fgb * tdy_zz_xxzz_0[j];

            tly_zzz_xxzz_0[j] = pa_z[j] * tly_zz_xxzz_0[j] + fl1_fx * tly_z_xxzz_0[j] + fl1_fx * tly_zz_xxz_0[j] - 0.5 * fl1_fx * tpx_zz_xxzz_0[j] -
                                fl1_fx * fl1_fgb * tdx_zz_xxzz_0[j];

            tlz_zzz_xxzz_0[j] = pa_z[j] * tlz_zz_xxzz_0[j] + fl1_fx * tlz_z_xxzz_0[j] + fl1_fx * tlz_zz_xxz_0[j];

            tlx_zzz_xyyy_0[j] =
                pa_z[j] * tlx_zz_xyyy_0[j] + fl1_fx * tlx_z_xyyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyyy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xyyy_0[j];

            tly_zzz_xyyy_0[j] =
                pa_z[j] * tly_zz_xyyy_0[j] + fl1_fx * tly_z_xyyy_0[j] - 0.5 * fl1_fx * tpx_zz_xyyy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xyyy_0[j];

            tlz_zzz_xyyy_0[j] = pa_z[j] * tlz_zz_xyyy_0[j] + fl1_fx * tlz_z_xyyy_0[j];

            tlx_zzz_xyyz_0[j] = pa_z[j] * tlx_zz_xyyz_0[j] + fl1_fx * tlx_z_xyyz_0[j] + 0.5 * fl1_fx * tlx_zz_xyy_0[j] +
                                0.5 * fl1_fx * tpy_zz_xyyz_0[j] + fl1_fx * fl1_fgb * tdy_zz_xyyz_0[j];

            tly_zzz_xyyz_0[j] = pa_z[j] * tly_zz_xyyz_0[j] + fl1_fx * tly_z_xyyz_0[j] + 0.5 * fl1_fx * tly_zz_xyy_0[j] -
                                0.5 * fl1_fx * tpx_zz_xyyz_0[j] - fl1_fx * fl1_fgb * tdx_zz_xyyz_0[j];

            tlz_zzz_xyyz_0[j] = pa_z[j] * tlz_zz_xyyz_0[j] + fl1_fx * tlz_z_xyyz_0[j] + 0.5 * fl1_fx * tlz_zz_xyy_0[j];

            tlx_zzz_xyzz_0[j] = pa_z[j] * tlx_zz_xyzz_0[j] + fl1_fx * tlx_z_xyzz_0[j] + fl1_fx * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyzz_0[j] +
                                fl1_fx * fl1_fgb * tdy_zz_xyzz_0[j];

            tly_zzz_xyzz_0[j] = pa_z[j] * tly_zz_xyzz_0[j] + fl1_fx * tly_z_xyzz_0[j] + fl1_fx * tly_zz_xyz_0[j] - 0.5 * fl1_fx * tpx_zz_xyzz_0[j] -
                                fl1_fx * fl1_fgb * tdx_zz_xyzz_0[j];

            tlz_zzz_xyzz_0[j] = pa_z[j] * tlz_zz_xyzz_0[j] + fl1_fx * tlz_z_xyzz_0[j] + fl1_fx * tlz_zz_xyz_0[j];

            tlx_zzz_xzzz_0[j] = pa_z[j] * tlx_zz_xzzz_0[j] + fl1_fx * tlx_z_xzzz_0[j] + 1.5 * fl1_fx * tlx_zz_xzz_0[j] +
                                0.5 * fl1_fx * tpy_zz_xzzz_0[j] + fl1_fx * fl1_fgb * tdy_zz_xzzz_0[j];

            tly_zzz_xzzz_0[j] = pa_z[j] * tly_zz_xzzz_0[j] + fl1_fx * tly_z_xzzz_0[j] + 1.5 * fl1_fx * tly_zz_xzz_0[j] -
                                0.5 * fl1_fx * tpx_zz_xzzz_0[j] - fl1_fx * fl1_fgb * tdx_zz_xzzz_0[j];

            tlz_zzz_xzzz_0[j] = pa_z[j] * tlz_zz_xzzz_0[j] + fl1_fx * tlz_z_xzzz_0[j] + 1.5 * fl1_fx * tlz_zz_xzz_0[j];

            tlx_zzz_yyyy_0[j] =
                pa_z[j] * tlx_zz_yyyy_0[j] + fl1_fx * tlx_z_yyyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyyy_0[j] + fl1_fx * fl1_fgb * tdy_zz_yyyy_0[j];

            tly_zzz_yyyy_0[j] =
                pa_z[j] * tly_zz_yyyy_0[j] + fl1_fx * tly_z_yyyy_0[j] - 0.5 * fl1_fx * tpx_zz_yyyy_0[j] - fl1_fx * fl1_fgb * tdx_zz_yyyy_0[j];

            tlz_zzz_yyyy_0[j] = pa_z[j] * tlz_zz_yyyy_0[j] + fl1_fx * tlz_z_yyyy_0[j];

            tlx_zzz_yyyz_0[j] = pa_z[j] * tlx_zz_yyyz_0[j] + fl1_fx * tlx_z_yyyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyy_0[j] +
                                0.5 * fl1_fx * tpy_zz_yyyz_0[j] + fl1_fx * fl1_fgb * tdy_zz_yyyz_0[j];

            tly_zzz_yyyz_0[j] = pa_z[j] * tly_zz_yyyz_0[j] + fl1_fx * tly_z_yyyz_0[j] + 0.5 * fl1_fx * tly_zz_yyy_0[j] -
                                0.5 * fl1_fx * tpx_zz_yyyz_0[j] - fl1_fx * fl1_fgb * tdx_zz_yyyz_0[j];

            tlz_zzz_yyyz_0[j] = pa_z[j] * tlz_zz_yyyz_0[j] + fl1_fx * tlz_z_yyyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyy_0[j];

            tlx_zzz_yyzz_0[j] = pa_z[j] * tlx_zz_yyzz_0[j] + fl1_fx * tlx_z_yyzz_0[j] + fl1_fx * tlx_zz_yyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyzz_0[j] +
                                fl1_fx * fl1_fgb * tdy_zz_yyzz_0[j];

            tly_zzz_yyzz_0[j] = pa_z[j] * tly_zz_yyzz_0[j] + fl1_fx * tly_z_yyzz_0[j] + fl1_fx * tly_zz_yyz_0[j] - 0.5 * fl1_fx * tpx_zz_yyzz_0[j] -
                                fl1_fx * fl1_fgb * tdx_zz_yyzz_0[j];

            tlz_zzz_yyzz_0[j] = pa_z[j] * tlz_zz_yyzz_0[j] + fl1_fx * tlz_z_yyzz_0[j] + fl1_fx * tlz_zz_yyz_0[j];

            tlx_zzz_yzzz_0[j] = pa_z[j] * tlx_zz_yzzz_0[j] + fl1_fx * tlx_z_yzzz_0[j] + 1.5 * fl1_fx * tlx_zz_yzz_0[j] +
                                0.5 * fl1_fx * tpy_zz_yzzz_0[j] + fl1_fx * fl1_fgb * tdy_zz_yzzz_0[j];

            tly_zzz_yzzz_0[j] = pa_z[j] * tly_zz_yzzz_0[j] + fl1_fx * tly_z_yzzz_0[j] + 1.5 * fl1_fx * tly_zz_yzz_0[j] -
                                0.5 * fl1_fx * tpx_zz_yzzz_0[j] - fl1_fx * fl1_fgb * tdx_zz_yzzz_0[j];

            tlz_zzz_yzzz_0[j] = pa_z[j] * tlz_zz_yzzz_0[j] + fl1_fx * tlz_z_yzzz_0[j] + 1.5 * fl1_fx * tlz_zz_yzz_0[j];

            tlx_zzz_zzzz_0[j] = pa_z[j] * tlx_zz_zzzz_0[j] + fl1_fx * tlx_z_zzzz_0[j] + 2.0 * fl1_fx * tlx_zz_zzz_0[j] +
                                0.5 * fl1_fx * tpy_zz_zzzz_0[j] + fl1_fx * fl1_fgb * tdy_zz_zzzz_0[j];

            tly_zzz_zzzz_0[j] = pa_z[j] * tly_zz_zzzz_0[j] + fl1_fx * tly_z_zzzz_0[j] + 2.0 * fl1_fx * tly_zz_zzz_0[j] -
                                0.5 * fl1_fx * tpx_zz_zzzz_0[j] - fl1_fx * fl1_fgb * tdx_zz_zzzz_0[j];

            tlz_zzz_zzzz_0[j] = pa_z[j] * tlz_zz_zzzz_0[j] + fl1_fx * tlz_z_zzzz_0[j] + 2.0 * fl1_fx * tlz_zz_zzz_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
