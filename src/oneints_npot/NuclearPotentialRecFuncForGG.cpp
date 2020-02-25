//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "NuclearPotentialRecFuncForGG.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForGG(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForGG_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForGG_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForGG_90_135(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForGG_135_180(
        primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForGG_180_225(
        primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForGG_0_45(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto ta_xxx_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx);

            auto ta_xxx_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 1);

            auto ta_xxx_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 2);

            auto ta_xxx_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 3);

            auto ta_xxx_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 4);

            auto ta_xxx_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 5);

            auto ta_xxx_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 6);

            auto ta_xxx_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 7);

            auto ta_xxx_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 8);

            auto ta_xxx_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 9);

            auto ta_xxx_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 10);

            auto ta_xxx_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 11);

            auto ta_xxx_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 12);

            auto ta_xxx_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 13);

            auto ta_xxx_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 14);

            auto ta_xxy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 15);

            auto ta_xxy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 16);

            auto ta_xxy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 17);

            auto ta_xxy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 18);

            auto ta_xxy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 19);

            auto ta_xxy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 20);

            auto ta_xxy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 21);

            auto ta_xxy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 22);

            auto ta_xxy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 23);

            auto ta_xxy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 24);

            auto ta_xxy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 25);

            auto ta_xxy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 26);

            auto ta_xxy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 27);

            auto ta_xxy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 28);

            auto ta_xxy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 29);

            auto ta_xxz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 30);

            auto ta_xxz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 31);

            auto ta_xxz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 32);

            auto ta_xxz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 33);

            auto ta_xxz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 34);

            auto ta_xxz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 35);

            auto ta_xxz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 36);

            auto ta_xxz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 37);

            auto ta_xxz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 38);

            auto ta_xxz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 39);

            auto ta_xxz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 40);

            auto ta_xxz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 41);

            auto ta_xxz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 42);

            auto ta_xxz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 43);

            auto ta_xxz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 44);

            auto ta_xxx_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx);

            auto ta_xxx_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 1);

            auto ta_xxx_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 2);

            auto ta_xxx_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 3);

            auto ta_xxx_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 4);

            auto ta_xxx_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 5);

            auto ta_xxx_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 6);

            auto ta_xxx_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 7);

            auto ta_xxx_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 8);

            auto ta_xxx_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 9);

            auto ta_xxx_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 10);

            auto ta_xxx_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 11);

            auto ta_xxx_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 12);

            auto ta_xxx_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 13);

            auto ta_xxx_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 14);

            auto ta_xxy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 15);

            auto ta_xxy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 16);

            auto ta_xxy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 17);

            auto ta_xxy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 18);

            auto ta_xxy_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 19);

            auto ta_xxy_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 20);

            auto ta_xxy_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 21);

            auto ta_xxy_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 22);

            auto ta_xxy_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 23);

            auto ta_xxy_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 24);

            auto ta_xxy_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 25);

            auto ta_xxy_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 26);

            auto ta_xxy_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 27);

            auto ta_xxy_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 28);

            auto ta_xxy_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 29);

            auto ta_xxz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 30);

            auto ta_xxz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 31);

            auto ta_xxz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 32);

            auto ta_xxz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 33);

            auto ta_xxz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 34);

            auto ta_xxz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 35);

            auto ta_xxz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 36);

            auto ta_xxz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 37);

            auto ta_xxz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 38);

            auto ta_xxz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 39);

            auto ta_xxz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 40);

            auto ta_xxz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 41);

            auto ta_xxz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 42);

            auto ta_xxz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 43);

            auto ta_xxz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 44);

            auto ta_xx_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx);

            auto ta_xx_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 1);

            auto ta_xx_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 2);

            auto ta_xx_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 3);

            auto ta_xx_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 4);

            auto ta_xx_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 5);

            auto ta_xx_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 6);

            auto ta_xx_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 7);

            auto ta_xx_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 8);

            auto ta_xx_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 9);

            auto ta_xx_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 10);

            auto ta_xx_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 11);

            auto ta_xx_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 12);

            auto ta_xx_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 13);

            auto ta_xx_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 14);

            auto ta_xy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 15);

            auto ta_xy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 16);

            auto ta_xy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 17);

            auto ta_xy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 18);

            auto ta_xy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 19);

            auto ta_xy_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 20);

            auto ta_xy_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 21);

            auto ta_xy_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 22);

            auto ta_xy_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 23);

            auto ta_xy_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 24);

            auto ta_xy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 25);

            auto ta_xy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 26);

            auto ta_xy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 27);

            auto ta_xy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 28);

            auto ta_xy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 29);

            auto ta_xz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 30);

            auto ta_xz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 31);

            auto ta_xz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 32);

            auto ta_xz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 33);

            auto ta_xz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 34);

            auto ta_xz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 35);

            auto ta_xz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 36);

            auto ta_xz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 37);

            auto ta_xz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 38);

            auto ta_xz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 39);

            auto ta_xz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 40);

            auto ta_xz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 41);

            auto ta_xz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 42);

            auto ta_xz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 43);

            auto ta_xz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 44);

            auto ta_xx_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx);

            auto ta_xx_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 1);

            auto ta_xx_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 2);

            auto ta_xx_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 3);

            auto ta_xx_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 4);

            auto ta_xx_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 5);

            auto ta_xx_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 6);

            auto ta_xx_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 7);

            auto ta_xx_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 8);

            auto ta_xx_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 9);

            auto ta_xx_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 10);

            auto ta_xx_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 11);

            auto ta_xx_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 12);

            auto ta_xx_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 13);

            auto ta_xx_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 14);

            auto ta_xy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 15);

            auto ta_xy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 16);

            auto ta_xy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 17);

            auto ta_xy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 18);

            auto ta_xy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 19);

            auto ta_xy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 20);

            auto ta_xy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 21);

            auto ta_xy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 22);

            auto ta_xy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 23);

            auto ta_xy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 24);

            auto ta_xy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 25);

            auto ta_xy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 26);

            auto ta_xy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 27);

            auto ta_xy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 28);

            auto ta_xy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 29);

            auto ta_xz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 30);

            auto ta_xz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 31);

            auto ta_xz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 32);

            auto ta_xz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 33);

            auto ta_xz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 34);

            auto ta_xz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 35);

            auto ta_xz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 36);

            auto ta_xz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 37);

            auto ta_xz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 38);

            auto ta_xz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 39);

            auto ta_xz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 40);

            auto ta_xz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 41);

            auto ta_xz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 42);

            auto ta_xz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 43);

            auto ta_xz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 44);

            auto ta_xxx_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx);

            auto ta_xxx_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 1);

            auto ta_xxx_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 2);

            auto ta_xxx_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 3);

            auto ta_xxx_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 4);

            auto ta_xxx_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 5);

            auto ta_xxx_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 6);

            auto ta_xxx_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 7);

            auto ta_xxx_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 8);

            auto ta_xxx_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 9);

            auto ta_xxy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 10);

            auto ta_xxy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 11);

            auto ta_xxy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 12);

            auto ta_xxy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 13);

            auto ta_xxy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 14);

            auto ta_xxy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 15);

            auto ta_xxy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 16);

            auto ta_xxy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 17);

            auto ta_xxy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 18);

            auto ta_xxy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 19);

            auto ta_xxz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 20);

            auto ta_xxz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 21);

            auto ta_xxz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 22);

            auto ta_xxz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 23);

            auto ta_xxz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 24);

            auto ta_xxz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 25);

            auto ta_xxz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 26);

            auto ta_xxz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 27);

            auto ta_xxz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 28);

            auto ta_xxz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 29);

            auto ta_xxx_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx);

            auto ta_xxx_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 1);

            auto ta_xxx_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 2);

            auto ta_xxx_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 3);

            auto ta_xxx_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 4);

            auto ta_xxx_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 5);

            auto ta_xxx_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 6);

            auto ta_xxx_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 7);

            auto ta_xxx_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 8);

            auto ta_xxx_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 9);

            auto ta_xxy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 10);

            auto ta_xxy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 11);

            auto ta_xxy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 12);

            auto ta_xxy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 13);

            auto ta_xxy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 14);

            auto ta_xxy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 15);

            auto ta_xxy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 16);

            auto ta_xxy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 17);

            auto ta_xxy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 18);

            auto ta_xxy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 19);

            auto ta_xxz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 20);

            auto ta_xxz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 21);

            auto ta_xxz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 22);

            auto ta_xxz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 23);

            auto ta_xxz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 24);

            auto ta_xxz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 25);

            auto ta_xxz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 26);

            auto ta_xxz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 27);

            auto ta_xxz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 28);

            auto ta_xxz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 29);

            // set up pointers to integrals

            auto ta_xxxx_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx);

            auto ta_xxxx_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 1);

            auto ta_xxxx_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 2);

            auto ta_xxxx_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 3);

            auto ta_xxxx_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 4);

            auto ta_xxxx_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 5);

            auto ta_xxxx_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 6);

            auto ta_xxxx_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 7);

            auto ta_xxxx_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 8);

            auto ta_xxxx_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 9);

            auto ta_xxxx_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 10);

            auto ta_xxxx_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 11);

            auto ta_xxxx_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 12);

            auto ta_xxxx_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 13);

            auto ta_xxxx_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 14);

            auto ta_xxxy_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 15);

            auto ta_xxxy_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 16);

            auto ta_xxxy_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 17);

            auto ta_xxxy_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 18);

            auto ta_xxxy_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 19);

            auto ta_xxxy_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 20);

            auto ta_xxxy_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 21);

            auto ta_xxxy_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 22);

            auto ta_xxxy_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 23);

            auto ta_xxxy_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 24);

            auto ta_xxxy_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 25);

            auto ta_xxxy_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 26);

            auto ta_xxxy_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 27);

            auto ta_xxxy_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 28);

            auto ta_xxxy_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 29);

            auto ta_xxxz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 30);

            auto ta_xxxz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 31);

            auto ta_xxxz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 32);

            auto ta_xxxz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 33);

            auto ta_xxxz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 34);

            auto ta_xxxz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 35);

            auto ta_xxxz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 36);

            auto ta_xxxz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 37);

            auto ta_xxxz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 38);

            auto ta_xxxz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 39);

            auto ta_xxxz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 40);

            auto ta_xxxz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 41);

            auto ta_xxxz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 42);

            auto ta_xxxz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 43);

            auto ta_xxxz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 44);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_xxxx_0, ta_xx_xxxx_1, ta_xx_xxxy_0, ta_xx_xxxy_1, \
                                         ta_xx_xxxz_0, ta_xx_xxxz_1, ta_xx_xxyy_0, ta_xx_xxyy_1, ta_xx_xxyz_0, ta_xx_xxyz_1, \
                                         ta_xx_xxzz_0, ta_xx_xxzz_1, ta_xx_xyyy_0, ta_xx_xyyy_1, ta_xx_xyyz_0, ta_xx_xyyz_1, \
                                         ta_xx_xyzz_0, ta_xx_xyzz_1, ta_xx_xzzz_0, ta_xx_xzzz_1, ta_xx_yyyy_0, ta_xx_yyyy_1, \
                                         ta_xx_yyyz_0, ta_xx_yyyz_1, ta_xx_yyzz_0, ta_xx_yyzz_1, ta_xx_yzzz_0, ta_xx_yzzz_1, \
                                         ta_xx_zzzz_0, ta_xx_zzzz_1, ta_xxx_xxx_0, ta_xxx_xxx_1, ta_xxx_xxxx_0, \
                                         ta_xxx_xxxx_1, ta_xxx_xxxy_0, ta_xxx_xxxy_1, ta_xxx_xxxz_0, ta_xxx_xxxz_1, \
                                         ta_xxx_xxy_0, ta_xxx_xxy_1, ta_xxx_xxyy_0, ta_xxx_xxyy_1, ta_xxx_xxyz_0, \
                                         ta_xxx_xxyz_1, ta_xxx_xxz_0, ta_xxx_xxz_1, ta_xxx_xxzz_0, ta_xxx_xxzz_1, \
                                         ta_xxx_xyy_0, ta_xxx_xyy_1, ta_xxx_xyyy_0, ta_xxx_xyyy_1, ta_xxx_xyyz_0, \
                                         ta_xxx_xyyz_1, ta_xxx_xyz_0, ta_xxx_xyz_1, ta_xxx_xyzz_0, ta_xxx_xyzz_1, \
                                         ta_xxx_xzz_0, ta_xxx_xzz_1, ta_xxx_xzzz_0, ta_xxx_xzzz_1, ta_xxx_yyy_0, \
                                         ta_xxx_yyy_1, ta_xxx_yyyy_0, ta_xxx_yyyy_1, ta_xxx_yyyz_0, ta_xxx_yyyz_1, \
                                         ta_xxx_yyz_0, ta_xxx_yyz_1, ta_xxx_yyzz_0, ta_xxx_yyzz_1, ta_xxx_yzz_0, \
                                         ta_xxx_yzz_1, ta_xxx_yzzz_0, ta_xxx_yzzz_1, ta_xxx_zzz_0, ta_xxx_zzz_1, \
                                         ta_xxx_zzzz_0, ta_xxx_zzzz_1, ta_xxxx_xxxx_0, ta_xxxx_xxxy_0, ta_xxxx_xxxz_0, \
                                         ta_xxxx_xxyy_0, ta_xxxx_xxyz_0, ta_xxxx_xxzz_0, ta_xxxx_xyyy_0, ta_xxxx_xyyz_0, \
                                         ta_xxxx_xyzz_0, ta_xxxx_xzzz_0, ta_xxxx_yyyy_0, ta_xxxx_yyyz_0, ta_xxxx_yyzz_0, \
                                         ta_xxxx_yzzz_0, ta_xxxx_zzzz_0, ta_xxxy_xxxx_0, ta_xxxy_xxxy_0, ta_xxxy_xxxz_0, \
                                         ta_xxxy_xxyy_0, ta_xxxy_xxyz_0, ta_xxxy_xxzz_0, ta_xxxy_xyyy_0, ta_xxxy_xyyz_0, \
                                         ta_xxxy_xyzz_0, ta_xxxy_xzzz_0, ta_xxxy_yyyy_0, ta_xxxy_yyyz_0, ta_xxxy_yyzz_0, \
                                         ta_xxxy_yzzz_0, ta_xxxy_zzzz_0, ta_xxxz_xxxx_0, ta_xxxz_xxxy_0, ta_xxxz_xxxz_0, \
                                         ta_xxxz_xxyy_0, ta_xxxz_xxyz_0, ta_xxxz_xxzz_0, ta_xxxz_xyyy_0, ta_xxxz_xyyz_0, \
                                         ta_xxxz_xyzz_0, ta_xxxz_xzzz_0, ta_xxxz_yyyy_0, ta_xxxz_yyyz_0, ta_xxxz_yyzz_0, \
                                         ta_xxxz_yzzz_0, ta_xxxz_zzzz_0, ta_xxy_xxx_0, ta_xxy_xxx_1, ta_xxy_xxxx_0, \
                                         ta_xxy_xxxx_1, ta_xxy_xxxy_0, ta_xxy_xxxy_1, ta_xxy_xxxz_0, ta_xxy_xxxz_1, \
                                         ta_xxy_xxy_0, ta_xxy_xxy_1, ta_xxy_xxyy_0, ta_xxy_xxyy_1, ta_xxy_xxyz_0, \
                                         ta_xxy_xxyz_1, ta_xxy_xxz_0, ta_xxy_xxz_1, ta_xxy_xxzz_0, ta_xxy_xxzz_1, \
                                         ta_xxy_xyy_0, ta_xxy_xyy_1, ta_xxy_xyyy_0, ta_xxy_xyyy_1, ta_xxy_xyyz_0, \
                                         ta_xxy_xyyz_1, ta_xxy_xyz_0, ta_xxy_xyz_1, ta_xxy_xyzz_0, ta_xxy_xyzz_1, \
                                         ta_xxy_xzz_0, ta_xxy_xzz_1, ta_xxy_xzzz_0, ta_xxy_xzzz_1, ta_xxy_yyy_0, \
                                         ta_xxy_yyy_1, ta_xxy_yyyy_0, ta_xxy_yyyy_1, ta_xxy_yyyz_0, ta_xxy_yyyz_1, \
                                         ta_xxy_yyz_0, ta_xxy_yyz_1, ta_xxy_yyzz_0, ta_xxy_yyzz_1, ta_xxy_yzz_0, \
                                         ta_xxy_yzz_1, ta_xxy_yzzz_0, ta_xxy_yzzz_1, ta_xxy_zzz_0, ta_xxy_zzz_1, \
                                         ta_xxy_zzzz_0, ta_xxy_zzzz_1, ta_xxz_xxx_0, ta_xxz_xxx_1, ta_xxz_xxxx_0, \
                                         ta_xxz_xxxx_1, ta_xxz_xxxy_0, ta_xxz_xxxy_1, ta_xxz_xxxz_0, ta_xxz_xxxz_1, \
                                         ta_xxz_xxy_0, ta_xxz_xxy_1, ta_xxz_xxyy_0, ta_xxz_xxyy_1, ta_xxz_xxyz_0, \
                                         ta_xxz_xxyz_1, ta_xxz_xxz_0, ta_xxz_xxz_1, ta_xxz_xxzz_0, ta_xxz_xxzz_1, \
                                         ta_xxz_xyy_0, ta_xxz_xyy_1, ta_xxz_xyyy_0, ta_xxz_xyyy_1, ta_xxz_xyyz_0, \
                                         ta_xxz_xyyz_1, ta_xxz_xyz_0, ta_xxz_xyz_1, ta_xxz_xyzz_0, ta_xxz_xyzz_1, \
                                         ta_xxz_xzz_0, ta_xxz_xzz_1, ta_xxz_xzzz_0, ta_xxz_xzzz_1, ta_xxz_yyy_0, \
                                         ta_xxz_yyy_1, ta_xxz_yyyy_0, ta_xxz_yyyy_1, ta_xxz_yyyz_0, ta_xxz_yyyz_1, \
                                         ta_xxz_yyz_0, ta_xxz_yyz_1, ta_xxz_yyzz_0, ta_xxz_yyzz_1, ta_xxz_yzz_0, \
                                         ta_xxz_yzz_1, ta_xxz_yzzz_0, ta_xxz_yzzz_1, ta_xxz_zzz_0, ta_xxz_zzz_1, \
                                         ta_xxz_zzzz_0, ta_xxz_zzzz_1, ta_xy_xxxx_0, ta_xy_xxxx_1, ta_xy_xxxy_0, \
                                         ta_xy_xxxy_1, ta_xy_xxxz_0, ta_xy_xxxz_1, ta_xy_xxyy_0, ta_xy_xxyy_1, ta_xy_xxyz_0, \
                                         ta_xy_xxyz_1, ta_xy_xxzz_0, ta_xy_xxzz_1, ta_xy_xyyy_0, ta_xy_xyyy_1, ta_xy_xyyz_0, \
                                         ta_xy_xyyz_1, ta_xy_xyzz_0, ta_xy_xyzz_1, ta_xy_xzzz_0, ta_xy_xzzz_1, ta_xy_yyyy_0, \
                                         ta_xy_yyyy_1, ta_xy_yyyz_0, ta_xy_yyyz_1, ta_xy_yyzz_0, ta_xy_yyzz_1, ta_xy_yzzz_0, \
                                         ta_xy_yzzz_1, ta_xy_zzzz_0, ta_xy_zzzz_1, ta_xz_xxxx_0, ta_xz_xxxx_1, ta_xz_xxxy_0, \
                                         ta_xz_xxxy_1, ta_xz_xxxz_0, ta_xz_xxxz_1, ta_xz_xxyy_0, ta_xz_xxyy_1, ta_xz_xxyz_0, \
                                         ta_xz_xxyz_1, ta_xz_xxzz_0, ta_xz_xxzz_1, ta_xz_xyyy_0, ta_xz_xyyy_1, ta_xz_xyyz_0, \
                                         ta_xz_xyyz_1, ta_xz_xyzz_0, ta_xz_xyzz_1, ta_xz_xzzz_0, ta_xz_xzzz_1, ta_xz_yyyy_0, \
                                         ta_xz_yyyy_1, ta_xz_yyyz_0, ta_xz_yyyz_1, ta_xz_yyzz_0, ta_xz_yyzz_1, ta_xz_yzzz_0, \
                                         ta_xz_yzzz_1, ta_xz_zzzz_0, ta_xz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxxx_xxxx_0[j] = pa_x[j] * ta_xxx_xxxx_0[j] - pc_x[j] * ta_xxx_xxxx_1[j] + 1.5 * fl1_fx * ta_xx_xxxx_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxxx_1[j] + 2.0 * fl1_fx * ta_xxx_xxx_0[j] - 2.0 * fl1_fx * ta_xxx_xxx_1[j];

                ta_xxxx_xxxy_0[j] = pa_x[j] * ta_xxx_xxxy_0[j] - pc_x[j] * ta_xxx_xxxy_1[j] + 1.5 * fl1_fx * ta_xx_xxxy_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxxy_1[j] + 1.5 * fl1_fx * ta_xxx_xxy_0[j] - 1.5 * fl1_fx * ta_xxx_xxy_1[j];

                ta_xxxx_xxxz_0[j] = pa_x[j] * ta_xxx_xxxz_0[j] - pc_x[j] * ta_xxx_xxxz_1[j] + 1.5 * fl1_fx * ta_xx_xxxz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxxz_1[j] + 1.5 * fl1_fx * ta_xxx_xxz_0[j] - 1.5 * fl1_fx * ta_xxx_xxz_1[j];

                ta_xxxx_xxyy_0[j] = pa_x[j] * ta_xxx_xxyy_0[j] - pc_x[j] * ta_xxx_xxyy_1[j] + 1.5 * fl1_fx * ta_xx_xxyy_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxyy_1[j] + fl1_fx * ta_xxx_xyy_0[j] - fl1_fx * ta_xxx_xyy_1[j];

                ta_xxxx_xxyz_0[j] = pa_x[j] * ta_xxx_xxyz_0[j] - pc_x[j] * ta_xxx_xxyz_1[j] + 1.5 * fl1_fx * ta_xx_xxyz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxyz_1[j] + fl1_fx * ta_xxx_xyz_0[j] - fl1_fx * ta_xxx_xyz_1[j];

                ta_xxxx_xxzz_0[j] = pa_x[j] * ta_xxx_xxzz_0[j] - pc_x[j] * ta_xxx_xxzz_1[j] + 1.5 * fl1_fx * ta_xx_xxzz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xxzz_1[j] + fl1_fx * ta_xxx_xzz_0[j] - fl1_fx * ta_xxx_xzz_1[j];

                ta_xxxx_xyyy_0[j] = pa_x[j] * ta_xxx_xyyy_0[j] - pc_x[j] * ta_xxx_xyyy_1[j] + 1.5 * fl1_fx * ta_xx_xyyy_0[j] -
                                    1.5 * fl1_fx * ta_xx_xyyy_1[j] + 0.5 * fl1_fx * ta_xxx_yyy_0[j] - 0.5 * fl1_fx * ta_xxx_yyy_1[j];

                ta_xxxx_xyyz_0[j] = pa_x[j] * ta_xxx_xyyz_0[j] - pc_x[j] * ta_xxx_xyyz_1[j] + 1.5 * fl1_fx * ta_xx_xyyz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xyyz_1[j] + 0.5 * fl1_fx * ta_xxx_yyz_0[j] - 0.5 * fl1_fx * ta_xxx_yyz_1[j];

                ta_xxxx_xyzz_0[j] = pa_x[j] * ta_xxx_xyzz_0[j] - pc_x[j] * ta_xxx_xyzz_1[j] + 1.5 * fl1_fx * ta_xx_xyzz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xyzz_1[j] + 0.5 * fl1_fx * ta_xxx_yzz_0[j] - 0.5 * fl1_fx * ta_xxx_yzz_1[j];

                ta_xxxx_xzzz_0[j] = pa_x[j] * ta_xxx_xzzz_0[j] - pc_x[j] * ta_xxx_xzzz_1[j] + 1.5 * fl1_fx * ta_xx_xzzz_0[j] -
                                    1.5 * fl1_fx * ta_xx_xzzz_1[j] + 0.5 * fl1_fx * ta_xxx_zzz_0[j] - 0.5 * fl1_fx * ta_xxx_zzz_1[j];

                ta_xxxx_yyyy_0[j] =
                    pa_x[j] * ta_xxx_yyyy_0[j] - pc_x[j] * ta_xxx_yyyy_1[j] + 1.5 * fl1_fx * ta_xx_yyyy_0[j] - 1.5 * fl1_fx * ta_xx_yyyy_1[j];

                ta_xxxx_yyyz_0[j] =
                    pa_x[j] * ta_xxx_yyyz_0[j] - pc_x[j] * ta_xxx_yyyz_1[j] + 1.5 * fl1_fx * ta_xx_yyyz_0[j] - 1.5 * fl1_fx * ta_xx_yyyz_1[j];

                ta_xxxx_yyzz_0[j] =
                    pa_x[j] * ta_xxx_yyzz_0[j] - pc_x[j] * ta_xxx_yyzz_1[j] + 1.5 * fl1_fx * ta_xx_yyzz_0[j] - 1.5 * fl1_fx * ta_xx_yyzz_1[j];

                ta_xxxx_yzzz_0[j] =
                    pa_x[j] * ta_xxx_yzzz_0[j] - pc_x[j] * ta_xxx_yzzz_1[j] + 1.5 * fl1_fx * ta_xx_yzzz_0[j] - 1.5 * fl1_fx * ta_xx_yzzz_1[j];

                ta_xxxx_zzzz_0[j] =
                    pa_x[j] * ta_xxx_zzzz_0[j] - pc_x[j] * ta_xxx_zzzz_1[j] + 1.5 * fl1_fx * ta_xx_zzzz_0[j] - 1.5 * fl1_fx * ta_xx_zzzz_1[j];

                ta_xxxy_xxxx_0[j] = pa_x[j] * ta_xxy_xxxx_0[j] - pc_x[j] * ta_xxy_xxxx_1[j] + fl1_fx * ta_xy_xxxx_0[j] - fl1_fx * ta_xy_xxxx_1[j] +
                                    2.0 * fl1_fx * ta_xxy_xxx_0[j] - 2.0 * fl1_fx * ta_xxy_xxx_1[j];

                ta_xxxy_xxxy_0[j] = pa_x[j] * ta_xxy_xxxy_0[j] - pc_x[j] * ta_xxy_xxxy_1[j] + fl1_fx * ta_xy_xxxy_0[j] - fl1_fx * ta_xy_xxxy_1[j] +
                                    1.5 * fl1_fx * ta_xxy_xxy_0[j] - 1.5 * fl1_fx * ta_xxy_xxy_1[j];

                ta_xxxy_xxxz_0[j] = pa_x[j] * ta_xxy_xxxz_0[j] - pc_x[j] * ta_xxy_xxxz_1[j] + fl1_fx * ta_xy_xxxz_0[j] - fl1_fx * ta_xy_xxxz_1[j] +
                                    1.5 * fl1_fx * ta_xxy_xxz_0[j] - 1.5 * fl1_fx * ta_xxy_xxz_1[j];

                ta_xxxy_xxyy_0[j] = pa_x[j] * ta_xxy_xxyy_0[j] - pc_x[j] * ta_xxy_xxyy_1[j] + fl1_fx * ta_xy_xxyy_0[j] - fl1_fx * ta_xy_xxyy_1[j] +
                                    fl1_fx * ta_xxy_xyy_0[j] - fl1_fx * ta_xxy_xyy_1[j];

                ta_xxxy_xxyz_0[j] = pa_x[j] * ta_xxy_xxyz_0[j] - pc_x[j] * ta_xxy_xxyz_1[j] + fl1_fx * ta_xy_xxyz_0[j] - fl1_fx * ta_xy_xxyz_1[j] +
                                    fl1_fx * ta_xxy_xyz_0[j] - fl1_fx * ta_xxy_xyz_1[j];

                ta_xxxy_xxzz_0[j] = pa_x[j] * ta_xxy_xxzz_0[j] - pc_x[j] * ta_xxy_xxzz_1[j] + fl1_fx * ta_xy_xxzz_0[j] - fl1_fx * ta_xy_xxzz_1[j] +
                                    fl1_fx * ta_xxy_xzz_0[j] - fl1_fx * ta_xxy_xzz_1[j];

                ta_xxxy_xyyy_0[j] = pa_x[j] * ta_xxy_xyyy_0[j] - pc_x[j] * ta_xxy_xyyy_1[j] + fl1_fx * ta_xy_xyyy_0[j] - fl1_fx * ta_xy_xyyy_1[j] +
                                    0.5 * fl1_fx * ta_xxy_yyy_0[j] - 0.5 * fl1_fx * ta_xxy_yyy_1[j];

                ta_xxxy_xyyz_0[j] = pa_x[j] * ta_xxy_xyyz_0[j] - pc_x[j] * ta_xxy_xyyz_1[j] + fl1_fx * ta_xy_xyyz_0[j] - fl1_fx * ta_xy_xyyz_1[j] +
                                    0.5 * fl1_fx * ta_xxy_yyz_0[j] - 0.5 * fl1_fx * ta_xxy_yyz_1[j];

                ta_xxxy_xyzz_0[j] = pa_x[j] * ta_xxy_xyzz_0[j] - pc_x[j] * ta_xxy_xyzz_1[j] + fl1_fx * ta_xy_xyzz_0[j] - fl1_fx * ta_xy_xyzz_1[j] +
                                    0.5 * fl1_fx * ta_xxy_yzz_0[j] - 0.5 * fl1_fx * ta_xxy_yzz_1[j];

                ta_xxxy_xzzz_0[j] = pa_x[j] * ta_xxy_xzzz_0[j] - pc_x[j] * ta_xxy_xzzz_1[j] + fl1_fx * ta_xy_xzzz_0[j] - fl1_fx * ta_xy_xzzz_1[j] +
                                    0.5 * fl1_fx * ta_xxy_zzz_0[j] - 0.5 * fl1_fx * ta_xxy_zzz_1[j];

                ta_xxxy_yyyy_0[j] = pa_x[j] * ta_xxy_yyyy_0[j] - pc_x[j] * ta_xxy_yyyy_1[j] + fl1_fx * ta_xy_yyyy_0[j] - fl1_fx * ta_xy_yyyy_1[j];

                ta_xxxy_yyyz_0[j] = pa_x[j] * ta_xxy_yyyz_0[j] - pc_x[j] * ta_xxy_yyyz_1[j] + fl1_fx * ta_xy_yyyz_0[j] - fl1_fx * ta_xy_yyyz_1[j];

                ta_xxxy_yyzz_0[j] = pa_x[j] * ta_xxy_yyzz_0[j] - pc_x[j] * ta_xxy_yyzz_1[j] + fl1_fx * ta_xy_yyzz_0[j] - fl1_fx * ta_xy_yyzz_1[j];

                ta_xxxy_yzzz_0[j] = pa_x[j] * ta_xxy_yzzz_0[j] - pc_x[j] * ta_xxy_yzzz_1[j] + fl1_fx * ta_xy_yzzz_0[j] - fl1_fx * ta_xy_yzzz_1[j];

                ta_xxxy_zzzz_0[j] = pa_x[j] * ta_xxy_zzzz_0[j] - pc_x[j] * ta_xxy_zzzz_1[j] + fl1_fx * ta_xy_zzzz_0[j] - fl1_fx * ta_xy_zzzz_1[j];

                ta_xxxz_xxxx_0[j] = pa_x[j] * ta_xxz_xxxx_0[j] - pc_x[j] * ta_xxz_xxxx_1[j] + fl1_fx * ta_xz_xxxx_0[j] - fl1_fx * ta_xz_xxxx_1[j] +
                                    2.0 * fl1_fx * ta_xxz_xxx_0[j] - 2.0 * fl1_fx * ta_xxz_xxx_1[j];

                ta_xxxz_xxxy_0[j] = pa_x[j] * ta_xxz_xxxy_0[j] - pc_x[j] * ta_xxz_xxxy_1[j] + fl1_fx * ta_xz_xxxy_0[j] - fl1_fx * ta_xz_xxxy_1[j] +
                                    1.5 * fl1_fx * ta_xxz_xxy_0[j] - 1.5 * fl1_fx * ta_xxz_xxy_1[j];

                ta_xxxz_xxxz_0[j] = pa_x[j] * ta_xxz_xxxz_0[j] - pc_x[j] * ta_xxz_xxxz_1[j] + fl1_fx * ta_xz_xxxz_0[j] - fl1_fx * ta_xz_xxxz_1[j] +
                                    1.5 * fl1_fx * ta_xxz_xxz_0[j] - 1.5 * fl1_fx * ta_xxz_xxz_1[j];

                ta_xxxz_xxyy_0[j] = pa_x[j] * ta_xxz_xxyy_0[j] - pc_x[j] * ta_xxz_xxyy_1[j] + fl1_fx * ta_xz_xxyy_0[j] - fl1_fx * ta_xz_xxyy_1[j] +
                                    fl1_fx * ta_xxz_xyy_0[j] - fl1_fx * ta_xxz_xyy_1[j];

                ta_xxxz_xxyz_0[j] = pa_x[j] * ta_xxz_xxyz_0[j] - pc_x[j] * ta_xxz_xxyz_1[j] + fl1_fx * ta_xz_xxyz_0[j] - fl1_fx * ta_xz_xxyz_1[j] +
                                    fl1_fx * ta_xxz_xyz_0[j] - fl1_fx * ta_xxz_xyz_1[j];

                ta_xxxz_xxzz_0[j] = pa_x[j] * ta_xxz_xxzz_0[j] - pc_x[j] * ta_xxz_xxzz_1[j] + fl1_fx * ta_xz_xxzz_0[j] - fl1_fx * ta_xz_xxzz_1[j] +
                                    fl1_fx * ta_xxz_xzz_0[j] - fl1_fx * ta_xxz_xzz_1[j];

                ta_xxxz_xyyy_0[j] = pa_x[j] * ta_xxz_xyyy_0[j] - pc_x[j] * ta_xxz_xyyy_1[j] + fl1_fx * ta_xz_xyyy_0[j] - fl1_fx * ta_xz_xyyy_1[j] +
                                    0.5 * fl1_fx * ta_xxz_yyy_0[j] - 0.5 * fl1_fx * ta_xxz_yyy_1[j];

                ta_xxxz_xyyz_0[j] = pa_x[j] * ta_xxz_xyyz_0[j] - pc_x[j] * ta_xxz_xyyz_1[j] + fl1_fx * ta_xz_xyyz_0[j] - fl1_fx * ta_xz_xyyz_1[j] +
                                    0.5 * fl1_fx * ta_xxz_yyz_0[j] - 0.5 * fl1_fx * ta_xxz_yyz_1[j];

                ta_xxxz_xyzz_0[j] = pa_x[j] * ta_xxz_xyzz_0[j] - pc_x[j] * ta_xxz_xyzz_1[j] + fl1_fx * ta_xz_xyzz_0[j] - fl1_fx * ta_xz_xyzz_1[j] +
                                    0.5 * fl1_fx * ta_xxz_yzz_0[j] - 0.5 * fl1_fx * ta_xxz_yzz_1[j];

                ta_xxxz_xzzz_0[j] = pa_x[j] * ta_xxz_xzzz_0[j] - pc_x[j] * ta_xxz_xzzz_1[j] + fl1_fx * ta_xz_xzzz_0[j] - fl1_fx * ta_xz_xzzz_1[j] +
                                    0.5 * fl1_fx * ta_xxz_zzz_0[j] - 0.5 * fl1_fx * ta_xxz_zzz_1[j];

                ta_xxxz_yyyy_0[j] = pa_x[j] * ta_xxz_yyyy_0[j] - pc_x[j] * ta_xxz_yyyy_1[j] + fl1_fx * ta_xz_yyyy_0[j] - fl1_fx * ta_xz_yyyy_1[j];

                ta_xxxz_yyyz_0[j] = pa_x[j] * ta_xxz_yyyz_0[j] - pc_x[j] * ta_xxz_yyyz_1[j] + fl1_fx * ta_xz_yyyz_0[j] - fl1_fx * ta_xz_yyyz_1[j];

                ta_xxxz_yyzz_0[j] = pa_x[j] * ta_xxz_yyzz_0[j] - pc_x[j] * ta_xxz_yyzz_1[j] + fl1_fx * ta_xz_yyzz_0[j] - fl1_fx * ta_xz_yyzz_1[j];

                ta_xxxz_yzzz_0[j] = pa_x[j] * ta_xxz_yzzz_0[j] - pc_x[j] * ta_xxz_yzzz_1[j] + fl1_fx * ta_xz_yzzz_0[j] - fl1_fx * ta_xz_yzzz_1[j];

                ta_xxxz_zzzz_0[j] = pa_x[j] * ta_xxz_zzzz_0[j] - pc_x[j] * ta_xxz_zzzz_1[j] + fl1_fx * ta_xz_zzzz_0[j] - fl1_fx * ta_xz_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGG_45_90(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto ta_xyy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 45);

            auto ta_xyy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 46);

            auto ta_xyy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 47);

            auto ta_xyy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 48);

            auto ta_xyy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 49);

            auto ta_xyy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 50);

            auto ta_xyy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 51);

            auto ta_xyy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 52);

            auto ta_xyy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 53);

            auto ta_xyy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 54);

            auto ta_xyy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 55);

            auto ta_xyy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 56);

            auto ta_xyy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 57);

            auto ta_xyy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 58);

            auto ta_xyy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 59);

            auto ta_xyz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 60);

            auto ta_xyz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 61);

            auto ta_xyz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 62);

            auto ta_xyz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 63);

            auto ta_xyz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 64);

            auto ta_xyz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 65);

            auto ta_xyz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 66);

            auto ta_xyz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 67);

            auto ta_xyz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 68);

            auto ta_xyz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 69);

            auto ta_xyz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 70);

            auto ta_xyz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 71);

            auto ta_xyz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 72);

            auto ta_xyz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 73);

            auto ta_xyz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 74);

            auto ta_xzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 75);

            auto ta_xzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 76);

            auto ta_xzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 77);

            auto ta_xzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 78);

            auto ta_xzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 79);

            auto ta_xzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 80);

            auto ta_xzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 81);

            auto ta_xzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 82);

            auto ta_xzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 83);

            auto ta_xzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 84);

            auto ta_xzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 85);

            auto ta_xzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 86);

            auto ta_xzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 87);

            auto ta_xzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 88);

            auto ta_xzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 89);

            auto ta_xyy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 45);

            auto ta_xyy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 46);

            auto ta_xyy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 47);

            auto ta_xyy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 48);

            auto ta_xyy_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 49);

            auto ta_xyy_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 50);

            auto ta_xyy_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 51);

            auto ta_xyy_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 52);

            auto ta_xyy_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 53);

            auto ta_xyy_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 54);

            auto ta_xyy_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 55);

            auto ta_xyy_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 56);

            auto ta_xyy_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 57);

            auto ta_xyy_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 58);

            auto ta_xyy_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 59);

            auto ta_xyz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 60);

            auto ta_xyz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 61);

            auto ta_xyz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 62);

            auto ta_xyz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 63);

            auto ta_xyz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 64);

            auto ta_xyz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 65);

            auto ta_xyz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 66);

            auto ta_xyz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 67);

            auto ta_xyz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 68);

            auto ta_xyz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 69);

            auto ta_xyz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 70);

            auto ta_xyz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 71);

            auto ta_xyz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 72);

            auto ta_xyz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 73);

            auto ta_xyz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 74);

            auto ta_xzz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 75);

            auto ta_xzz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 76);

            auto ta_xzz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 77);

            auto ta_xzz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 78);

            auto ta_xzz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 79);

            auto ta_xzz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 80);

            auto ta_xzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 81);

            auto ta_xzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 82);

            auto ta_xzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 83);

            auto ta_xzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 84);

            auto ta_xzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 85);

            auto ta_xzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 86);

            auto ta_xzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 87);

            auto ta_xzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 88);

            auto ta_xzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 89);

            auto ta_yy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 45);

            auto ta_yy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 46);

            auto ta_yy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 47);

            auto ta_yy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 48);

            auto ta_yy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 49);

            auto ta_yy_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 50);

            auto ta_yy_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 51);

            auto ta_yy_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 52);

            auto ta_yy_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 53);

            auto ta_yy_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 54);

            auto ta_yy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 55);

            auto ta_yy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 56);

            auto ta_yy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 57);

            auto ta_yy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 58);

            auto ta_yy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 59);

            auto ta_yz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 60);

            auto ta_yz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 61);

            auto ta_yz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 62);

            auto ta_yz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 63);

            auto ta_yz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 64);

            auto ta_yz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 65);

            auto ta_yz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 66);

            auto ta_yz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 67);

            auto ta_yz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 68);

            auto ta_yz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 69);

            auto ta_yz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 70);

            auto ta_yz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 71);

            auto ta_yz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 72);

            auto ta_yz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 73);

            auto ta_yz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 74);

            auto ta_zz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 75);

            auto ta_zz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 76);

            auto ta_zz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 77);

            auto ta_zz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 78);

            auto ta_zz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 79);

            auto ta_zz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 80);

            auto ta_zz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 81);

            auto ta_zz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 82);

            auto ta_zz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 83);

            auto ta_zz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 84);

            auto ta_zz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 85);

            auto ta_zz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 86);

            auto ta_zz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 87);

            auto ta_zz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 88);

            auto ta_zz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 89);

            auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45);

            auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46);

            auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47);

            auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48);

            auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49);

            auto ta_yy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 50);

            auto ta_yy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 51);

            auto ta_yy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 52);

            auto ta_yy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 53);

            auto ta_yy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 54);

            auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55);

            auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56);

            auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57);

            auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58);

            auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59);

            auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60);

            auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61);

            auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62);

            auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63);

            auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64);

            auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65);

            auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66);

            auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67);

            auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68);

            auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69);

            auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70);

            auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71);

            auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72);

            auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73);

            auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74);

            auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75);

            auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76);

            auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77);

            auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78);

            auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79);

            auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80);

            auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81);

            auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82);

            auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83);

            auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84);

            auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85);

            auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86);

            auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87);

            auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88);

            auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89);

            auto ta_xyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 30);

            auto ta_xyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 31);

            auto ta_xyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 32);

            auto ta_xyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 33);

            auto ta_xyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 34);

            auto ta_xyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 35);

            auto ta_xyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 36);

            auto ta_xyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 37);

            auto ta_xyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 38);

            auto ta_xyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 39);

            auto ta_xyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 40);

            auto ta_xyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 41);

            auto ta_xyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 42);

            auto ta_xyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 43);

            auto ta_xyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 44);

            auto ta_xyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 45);

            auto ta_xyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 46);

            auto ta_xyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 47);

            auto ta_xyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 48);

            auto ta_xyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 49);

            auto ta_xzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 50);

            auto ta_xzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 51);

            auto ta_xzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 52);

            auto ta_xzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 53);

            auto ta_xzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 54);

            auto ta_xzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 55);

            auto ta_xzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 56);

            auto ta_xzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 57);

            auto ta_xzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 58);

            auto ta_xzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 59);

            auto ta_xyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 30);

            auto ta_xyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 31);

            auto ta_xyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 32);

            auto ta_xyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 33);

            auto ta_xyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 34);

            auto ta_xyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 35);

            auto ta_xyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 36);

            auto ta_xyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 37);

            auto ta_xyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 38);

            auto ta_xyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 39);

            auto ta_xyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 40);

            auto ta_xyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 41);

            auto ta_xyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 42);

            auto ta_xyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 43);

            auto ta_xyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 44);

            auto ta_xyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 45);

            auto ta_xyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 46);

            auto ta_xyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 47);

            auto ta_xyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 48);

            auto ta_xyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 49);

            auto ta_xzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 50);

            auto ta_xzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 51);

            auto ta_xzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 52);

            auto ta_xzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 53);

            auto ta_xzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 54);

            auto ta_xzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 55);

            auto ta_xzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 56);

            auto ta_xzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 57);

            auto ta_xzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 58);

            auto ta_xzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 59);

            // set up pointers to integrals

            auto ta_xxyy_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 45);

            auto ta_xxyy_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 46);

            auto ta_xxyy_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 47);

            auto ta_xxyy_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 48);

            auto ta_xxyy_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 49);

            auto ta_xxyy_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 50);

            auto ta_xxyy_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 51);

            auto ta_xxyy_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 52);

            auto ta_xxyy_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 53);

            auto ta_xxyy_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 54);

            auto ta_xxyy_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 55);

            auto ta_xxyy_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 56);

            auto ta_xxyy_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 57);

            auto ta_xxyy_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 58);

            auto ta_xxyy_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 59);

            auto ta_xxyz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 60);

            auto ta_xxyz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 61);

            auto ta_xxyz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 62);

            auto ta_xxyz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 63);

            auto ta_xxyz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 64);

            auto ta_xxyz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 65);

            auto ta_xxyz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 66);

            auto ta_xxyz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 67);

            auto ta_xxyz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 68);

            auto ta_xxyz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 69);

            auto ta_xxyz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 70);

            auto ta_xxyz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 71);

            auto ta_xxyz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 72);

            auto ta_xxyz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 73);

            auto ta_xxyz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 74);

            auto ta_xxzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 75);

            auto ta_xxzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 76);

            auto ta_xxzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 77);

            auto ta_xxzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 78);

            auto ta_xxzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 79);

            auto ta_xxzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 80);

            auto ta_xxzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 81);

            auto ta_xxzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 82);

            auto ta_xxzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 83);

            auto ta_xxzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 84);

            auto ta_xxzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 85);

            auto ta_xxzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 86);

            auto ta_xxzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 87);

            auto ta_xxzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 88);

            auto ta_xxzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 89);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxyy_xxxx_0, ta_xxyy_xxxy_0, ta_xxyy_xxxz_0, \
                                         ta_xxyy_xxyy_0, ta_xxyy_xxyz_0, ta_xxyy_xxzz_0, ta_xxyy_xyyy_0, ta_xxyy_xyyz_0, \
                                         ta_xxyy_xyzz_0, ta_xxyy_xzzz_0, ta_xxyy_yyyy_0, ta_xxyy_yyyz_0, ta_xxyy_yyzz_0, \
                                         ta_xxyy_yzzz_0, ta_xxyy_zzzz_0, ta_xxyz_xxxx_0, ta_xxyz_xxxy_0, ta_xxyz_xxxz_0, \
                                         ta_xxyz_xxyy_0, ta_xxyz_xxyz_0, ta_xxyz_xxzz_0, ta_xxyz_xyyy_0, ta_xxyz_xyyz_0, \
                                         ta_xxyz_xyzz_0, ta_xxyz_xzzz_0, ta_xxyz_yyyy_0, ta_xxyz_yyyz_0, ta_xxyz_yyzz_0, \
                                         ta_xxyz_yzzz_0, ta_xxyz_zzzz_0, ta_xxzz_xxxx_0, ta_xxzz_xxxy_0, ta_xxzz_xxxz_0, \
                                         ta_xxzz_xxyy_0, ta_xxzz_xxyz_0, ta_xxzz_xxzz_0, ta_xxzz_xyyy_0, ta_xxzz_xyyz_0, \
                                         ta_xxzz_xyzz_0, ta_xxzz_xzzz_0, ta_xxzz_yyyy_0, ta_xxzz_yyyz_0, ta_xxzz_yyzz_0, \
                                         ta_xxzz_yzzz_0, ta_xxzz_zzzz_0, ta_xyy_xxx_0, ta_xyy_xxx_1, ta_xyy_xxxx_0, \
                                         ta_xyy_xxxx_1, ta_xyy_xxxy_0, ta_xyy_xxxy_1, ta_xyy_xxxz_0, ta_xyy_xxxz_1, \
                                         ta_xyy_xxy_0, ta_xyy_xxy_1, ta_xyy_xxyy_0, ta_xyy_xxyy_1, ta_xyy_xxyz_0, \
                                         ta_xyy_xxyz_1, ta_xyy_xxz_0, ta_xyy_xxz_1, ta_xyy_xxzz_0, ta_xyy_xxzz_1, \
                                         ta_xyy_xyy_0, ta_xyy_xyy_1, ta_xyy_xyyy_0, ta_xyy_xyyy_1, ta_xyy_xyyz_0, \
                                         ta_xyy_xyyz_1, ta_xyy_xyz_0, ta_xyy_xyz_1, ta_xyy_xyzz_0, ta_xyy_xyzz_1, \
                                         ta_xyy_xzz_0, ta_xyy_xzz_1, ta_xyy_xzzz_0, ta_xyy_xzzz_1, ta_xyy_yyy_0, \
                                         ta_xyy_yyy_1, ta_xyy_yyyy_0, ta_xyy_yyyy_1, ta_xyy_yyyz_0, ta_xyy_yyyz_1, \
                                         ta_xyy_yyz_0, ta_xyy_yyz_1, ta_xyy_yyzz_0, ta_xyy_yyzz_1, ta_xyy_yzz_0, \
                                         ta_xyy_yzz_1, ta_xyy_yzzz_0, ta_xyy_yzzz_1, ta_xyy_zzz_0, ta_xyy_zzz_1, \
                                         ta_xyy_zzzz_0, ta_xyy_zzzz_1, ta_xyz_xxx_0, ta_xyz_xxx_1, ta_xyz_xxxx_0, \
                                         ta_xyz_xxxx_1, ta_xyz_xxxy_0, ta_xyz_xxxy_1, ta_xyz_xxxz_0, ta_xyz_xxxz_1, \
                                         ta_xyz_xxy_0, ta_xyz_xxy_1, ta_xyz_xxyy_0, ta_xyz_xxyy_1, ta_xyz_xxyz_0, \
                                         ta_xyz_xxyz_1, ta_xyz_xxz_0, ta_xyz_xxz_1, ta_xyz_xxzz_0, ta_xyz_xxzz_1, \
                                         ta_xyz_xyy_0, ta_xyz_xyy_1, ta_xyz_xyyy_0, ta_xyz_xyyy_1, ta_xyz_xyyz_0, \
                                         ta_xyz_xyyz_1, ta_xyz_xyz_0, ta_xyz_xyz_1, ta_xyz_xyzz_0, ta_xyz_xyzz_1, \
                                         ta_xyz_xzz_0, ta_xyz_xzz_1, ta_xyz_xzzz_0, ta_xyz_xzzz_1, ta_xyz_yyy_0, \
                                         ta_xyz_yyy_1, ta_xyz_yyyy_0, ta_xyz_yyyy_1, ta_xyz_yyyz_0, ta_xyz_yyyz_1, \
                                         ta_xyz_yyz_0, ta_xyz_yyz_1, ta_xyz_yyzz_0, ta_xyz_yyzz_1, ta_xyz_yzz_0, \
                                         ta_xyz_yzz_1, ta_xyz_yzzz_0, ta_xyz_yzzz_1, ta_xyz_zzz_0, ta_xyz_zzz_1, \
                                         ta_xyz_zzzz_0, ta_xyz_zzzz_1, ta_xzz_xxx_0, ta_xzz_xxx_1, ta_xzz_xxxx_0, \
                                         ta_xzz_xxxx_1, ta_xzz_xxxy_0, ta_xzz_xxxy_1, ta_xzz_xxxz_0, ta_xzz_xxxz_1, \
                                         ta_xzz_xxy_0, ta_xzz_xxy_1, ta_xzz_xxyy_0, ta_xzz_xxyy_1, ta_xzz_xxyz_0, \
                                         ta_xzz_xxyz_1, ta_xzz_xxz_0, ta_xzz_xxz_1, ta_xzz_xxzz_0, ta_xzz_xxzz_1, \
                                         ta_xzz_xyy_0, ta_xzz_xyy_1, ta_xzz_xyyy_0, ta_xzz_xyyy_1, ta_xzz_xyyz_0, \
                                         ta_xzz_xyyz_1, ta_xzz_xyz_0, ta_xzz_xyz_1, ta_xzz_xyzz_0, ta_xzz_xyzz_1, \
                                         ta_xzz_xzz_0, ta_xzz_xzz_1, ta_xzz_xzzz_0, ta_xzz_xzzz_1, ta_xzz_yyy_0, \
                                         ta_xzz_yyy_1, ta_xzz_yyyy_0, ta_xzz_yyyy_1, ta_xzz_yyyz_0, ta_xzz_yyyz_1, \
                                         ta_xzz_yyz_0, ta_xzz_yyz_1, ta_xzz_yyzz_0, ta_xzz_yyzz_1, ta_xzz_yzz_0, \
                                         ta_xzz_yzz_1, ta_xzz_yzzz_0, ta_xzz_yzzz_1, ta_xzz_zzz_0, ta_xzz_zzz_1, \
                                         ta_xzz_zzzz_0, ta_xzz_zzzz_1, ta_yy_xxxx_0, ta_yy_xxxx_1, ta_yy_xxxy_0, \
                                         ta_yy_xxxy_1, ta_yy_xxxz_0, ta_yy_xxxz_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyz_0, \
                                         ta_yy_xxyz_1, ta_yy_xxzz_0, ta_yy_xxzz_1, ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyz_0, \
                                         ta_yy_xyyz_1, ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xzzz_0, ta_yy_xzzz_1, ta_yy_yyyy_0, \
                                         ta_yy_yyyy_1, ta_yy_yyyz_0, ta_yy_yyyz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yzzz_0, \
                                         ta_yy_yzzz_1, ta_yy_zzzz_0, ta_yy_zzzz_1, ta_yz_xxxx_0, ta_yz_xxxx_1, ta_yz_xxxy_0, \
                                         ta_yz_xxxy_1, ta_yz_xxxz_0, ta_yz_xxxz_1, ta_yz_xxyy_0, ta_yz_xxyy_1, ta_yz_xxyz_0, \
                                         ta_yz_xxyz_1, ta_yz_xxzz_0, ta_yz_xxzz_1, ta_yz_xyyy_0, ta_yz_xyyy_1, ta_yz_xyyz_0, \
                                         ta_yz_xyyz_1, ta_yz_xyzz_0, ta_yz_xyzz_1, ta_yz_xzzz_0, ta_yz_xzzz_1, ta_yz_yyyy_0, \
                                         ta_yz_yyyy_1, ta_yz_yyyz_0, ta_yz_yyyz_1, ta_yz_yyzz_0, ta_yz_yyzz_1, ta_yz_yzzz_0, \
                                         ta_yz_yzzz_1, ta_yz_zzzz_0, ta_yz_zzzz_1, ta_zz_xxxx_0, ta_zz_xxxx_1, ta_zz_xxxy_0, \
                                         ta_zz_xxxy_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyz_0, \
                                         ta_zz_xxyz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xyyy_0, ta_zz_xyyy_1, ta_zz_xyyz_0, \
                                         ta_zz_xyyz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_yyyy_0, \
                                         ta_zz_yyyy_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yzzz_0, \
                                         ta_zz_yzzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxyy_xxxx_0[j] = pa_x[j] * ta_xyy_xxxx_0[j] - pc_x[j] * ta_xyy_xxxx_1[j] + 0.5 * fl1_fx * ta_yy_xxxx_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxxx_1[j] + 2.0 * fl1_fx * ta_xyy_xxx_0[j] - 2.0 * fl1_fx * ta_xyy_xxx_1[j];

                ta_xxyy_xxxy_0[j] = pa_x[j] * ta_xyy_xxxy_0[j] - pc_x[j] * ta_xyy_xxxy_1[j] + 0.5 * fl1_fx * ta_yy_xxxy_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxxy_1[j] + 1.5 * fl1_fx * ta_xyy_xxy_0[j] - 1.5 * fl1_fx * ta_xyy_xxy_1[j];

                ta_xxyy_xxxz_0[j] = pa_x[j] * ta_xyy_xxxz_0[j] - pc_x[j] * ta_xyy_xxxz_1[j] + 0.5 * fl1_fx * ta_yy_xxxz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxxz_1[j] + 1.5 * fl1_fx * ta_xyy_xxz_0[j] - 1.5 * fl1_fx * ta_xyy_xxz_1[j];

                ta_xxyy_xxyy_0[j] = pa_x[j] * ta_xyy_xxyy_0[j] - pc_x[j] * ta_xyy_xxyy_1[j] + 0.5 * fl1_fx * ta_yy_xxyy_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxyy_1[j] + fl1_fx * ta_xyy_xyy_0[j] - fl1_fx * ta_xyy_xyy_1[j];

                ta_xxyy_xxyz_0[j] = pa_x[j] * ta_xyy_xxyz_0[j] - pc_x[j] * ta_xyy_xxyz_1[j] + 0.5 * fl1_fx * ta_yy_xxyz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxyz_1[j] + fl1_fx * ta_xyy_xyz_0[j] - fl1_fx * ta_xyy_xyz_1[j];

                ta_xxyy_xxzz_0[j] = pa_x[j] * ta_xyy_xxzz_0[j] - pc_x[j] * ta_xyy_xxzz_1[j] + 0.5 * fl1_fx * ta_yy_xxzz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xxzz_1[j] + fl1_fx * ta_xyy_xzz_0[j] - fl1_fx * ta_xyy_xzz_1[j];

                ta_xxyy_xyyy_0[j] = pa_x[j] * ta_xyy_xyyy_0[j] - pc_x[j] * ta_xyy_xyyy_1[j] + 0.5 * fl1_fx * ta_yy_xyyy_0[j] -
                                    0.5 * fl1_fx * ta_yy_xyyy_1[j] + 0.5 * fl1_fx * ta_xyy_yyy_0[j] - 0.5 * fl1_fx * ta_xyy_yyy_1[j];

                ta_xxyy_xyyz_0[j] = pa_x[j] * ta_xyy_xyyz_0[j] - pc_x[j] * ta_xyy_xyyz_1[j] + 0.5 * fl1_fx * ta_yy_xyyz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xyyz_1[j] + 0.5 * fl1_fx * ta_xyy_yyz_0[j] - 0.5 * fl1_fx * ta_xyy_yyz_1[j];

                ta_xxyy_xyzz_0[j] = pa_x[j] * ta_xyy_xyzz_0[j] - pc_x[j] * ta_xyy_xyzz_1[j] + 0.5 * fl1_fx * ta_yy_xyzz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xyzz_1[j] + 0.5 * fl1_fx * ta_xyy_yzz_0[j] - 0.5 * fl1_fx * ta_xyy_yzz_1[j];

                ta_xxyy_xzzz_0[j] = pa_x[j] * ta_xyy_xzzz_0[j] - pc_x[j] * ta_xyy_xzzz_1[j] + 0.5 * fl1_fx * ta_yy_xzzz_0[j] -
                                    0.5 * fl1_fx * ta_yy_xzzz_1[j] + 0.5 * fl1_fx * ta_xyy_zzz_0[j] - 0.5 * fl1_fx * ta_xyy_zzz_1[j];

                ta_xxyy_yyyy_0[j] =
                    pa_x[j] * ta_xyy_yyyy_0[j] - pc_x[j] * ta_xyy_yyyy_1[j] + 0.5 * fl1_fx * ta_yy_yyyy_0[j] - 0.5 * fl1_fx * ta_yy_yyyy_1[j];

                ta_xxyy_yyyz_0[j] =
                    pa_x[j] * ta_xyy_yyyz_0[j] - pc_x[j] * ta_xyy_yyyz_1[j] + 0.5 * fl1_fx * ta_yy_yyyz_0[j] - 0.5 * fl1_fx * ta_yy_yyyz_1[j];

                ta_xxyy_yyzz_0[j] =
                    pa_x[j] * ta_xyy_yyzz_0[j] - pc_x[j] * ta_xyy_yyzz_1[j] + 0.5 * fl1_fx * ta_yy_yyzz_0[j] - 0.5 * fl1_fx * ta_yy_yyzz_1[j];

                ta_xxyy_yzzz_0[j] =
                    pa_x[j] * ta_xyy_yzzz_0[j] - pc_x[j] * ta_xyy_yzzz_1[j] + 0.5 * fl1_fx * ta_yy_yzzz_0[j] - 0.5 * fl1_fx * ta_yy_yzzz_1[j];

                ta_xxyy_zzzz_0[j] =
                    pa_x[j] * ta_xyy_zzzz_0[j] - pc_x[j] * ta_xyy_zzzz_1[j] + 0.5 * fl1_fx * ta_yy_zzzz_0[j] - 0.5 * fl1_fx * ta_yy_zzzz_1[j];

                ta_xxyz_xxxx_0[j] = pa_x[j] * ta_xyz_xxxx_0[j] - pc_x[j] * ta_xyz_xxxx_1[j] + 0.5 * fl1_fx * ta_yz_xxxx_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxxx_1[j] + 2.0 * fl1_fx * ta_xyz_xxx_0[j] - 2.0 * fl1_fx * ta_xyz_xxx_1[j];

                ta_xxyz_xxxy_0[j] = pa_x[j] * ta_xyz_xxxy_0[j] - pc_x[j] * ta_xyz_xxxy_1[j] + 0.5 * fl1_fx * ta_yz_xxxy_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxxy_1[j] + 1.5 * fl1_fx * ta_xyz_xxy_0[j] - 1.5 * fl1_fx * ta_xyz_xxy_1[j];

                ta_xxyz_xxxz_0[j] = pa_x[j] * ta_xyz_xxxz_0[j] - pc_x[j] * ta_xyz_xxxz_1[j] + 0.5 * fl1_fx * ta_yz_xxxz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxxz_1[j] + 1.5 * fl1_fx * ta_xyz_xxz_0[j] - 1.5 * fl1_fx * ta_xyz_xxz_1[j];

                ta_xxyz_xxyy_0[j] = pa_x[j] * ta_xyz_xxyy_0[j] - pc_x[j] * ta_xyz_xxyy_1[j] + 0.5 * fl1_fx * ta_yz_xxyy_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxyy_1[j] + fl1_fx * ta_xyz_xyy_0[j] - fl1_fx * ta_xyz_xyy_1[j];

                ta_xxyz_xxyz_0[j] = pa_x[j] * ta_xyz_xxyz_0[j] - pc_x[j] * ta_xyz_xxyz_1[j] + 0.5 * fl1_fx * ta_yz_xxyz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxyz_1[j] + fl1_fx * ta_xyz_xyz_0[j] - fl1_fx * ta_xyz_xyz_1[j];

                ta_xxyz_xxzz_0[j] = pa_x[j] * ta_xyz_xxzz_0[j] - pc_x[j] * ta_xyz_xxzz_1[j] + 0.5 * fl1_fx * ta_yz_xxzz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xxzz_1[j] + fl1_fx * ta_xyz_xzz_0[j] - fl1_fx * ta_xyz_xzz_1[j];

                ta_xxyz_xyyy_0[j] = pa_x[j] * ta_xyz_xyyy_0[j] - pc_x[j] * ta_xyz_xyyy_1[j] + 0.5 * fl1_fx * ta_yz_xyyy_0[j] -
                                    0.5 * fl1_fx * ta_yz_xyyy_1[j] + 0.5 * fl1_fx * ta_xyz_yyy_0[j] - 0.5 * fl1_fx * ta_xyz_yyy_1[j];

                ta_xxyz_xyyz_0[j] = pa_x[j] * ta_xyz_xyyz_0[j] - pc_x[j] * ta_xyz_xyyz_1[j] + 0.5 * fl1_fx * ta_yz_xyyz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xyyz_1[j] + 0.5 * fl1_fx * ta_xyz_yyz_0[j] - 0.5 * fl1_fx * ta_xyz_yyz_1[j];

                ta_xxyz_xyzz_0[j] = pa_x[j] * ta_xyz_xyzz_0[j] - pc_x[j] * ta_xyz_xyzz_1[j] + 0.5 * fl1_fx * ta_yz_xyzz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xyzz_1[j] + 0.5 * fl1_fx * ta_xyz_yzz_0[j] - 0.5 * fl1_fx * ta_xyz_yzz_1[j];

                ta_xxyz_xzzz_0[j] = pa_x[j] * ta_xyz_xzzz_0[j] - pc_x[j] * ta_xyz_xzzz_1[j] + 0.5 * fl1_fx * ta_yz_xzzz_0[j] -
                                    0.5 * fl1_fx * ta_yz_xzzz_1[j] + 0.5 * fl1_fx * ta_xyz_zzz_0[j] - 0.5 * fl1_fx * ta_xyz_zzz_1[j];

                ta_xxyz_yyyy_0[j] =
                    pa_x[j] * ta_xyz_yyyy_0[j] - pc_x[j] * ta_xyz_yyyy_1[j] + 0.5 * fl1_fx * ta_yz_yyyy_0[j] - 0.5 * fl1_fx * ta_yz_yyyy_1[j];

                ta_xxyz_yyyz_0[j] =
                    pa_x[j] * ta_xyz_yyyz_0[j] - pc_x[j] * ta_xyz_yyyz_1[j] + 0.5 * fl1_fx * ta_yz_yyyz_0[j] - 0.5 * fl1_fx * ta_yz_yyyz_1[j];

                ta_xxyz_yyzz_0[j] =
                    pa_x[j] * ta_xyz_yyzz_0[j] - pc_x[j] * ta_xyz_yyzz_1[j] + 0.5 * fl1_fx * ta_yz_yyzz_0[j] - 0.5 * fl1_fx * ta_yz_yyzz_1[j];

                ta_xxyz_yzzz_0[j] =
                    pa_x[j] * ta_xyz_yzzz_0[j] - pc_x[j] * ta_xyz_yzzz_1[j] + 0.5 * fl1_fx * ta_yz_yzzz_0[j] - 0.5 * fl1_fx * ta_yz_yzzz_1[j];

                ta_xxyz_zzzz_0[j] =
                    pa_x[j] * ta_xyz_zzzz_0[j] - pc_x[j] * ta_xyz_zzzz_1[j] + 0.5 * fl1_fx * ta_yz_zzzz_0[j] - 0.5 * fl1_fx * ta_yz_zzzz_1[j];

                ta_xxzz_xxxx_0[j] = pa_x[j] * ta_xzz_xxxx_0[j] - pc_x[j] * ta_xzz_xxxx_1[j] + 0.5 * fl1_fx * ta_zz_xxxx_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxxx_1[j] + 2.0 * fl1_fx * ta_xzz_xxx_0[j] - 2.0 * fl1_fx * ta_xzz_xxx_1[j];

                ta_xxzz_xxxy_0[j] = pa_x[j] * ta_xzz_xxxy_0[j] - pc_x[j] * ta_xzz_xxxy_1[j] + 0.5 * fl1_fx * ta_zz_xxxy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxxy_1[j] + 1.5 * fl1_fx * ta_xzz_xxy_0[j] - 1.5 * fl1_fx * ta_xzz_xxy_1[j];

                ta_xxzz_xxxz_0[j] = pa_x[j] * ta_xzz_xxxz_0[j] - pc_x[j] * ta_xzz_xxxz_1[j] + 0.5 * fl1_fx * ta_zz_xxxz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxxz_1[j] + 1.5 * fl1_fx * ta_xzz_xxz_0[j] - 1.5 * fl1_fx * ta_xzz_xxz_1[j];

                ta_xxzz_xxyy_0[j] = pa_x[j] * ta_xzz_xxyy_0[j] - pc_x[j] * ta_xzz_xxyy_1[j] + 0.5 * fl1_fx * ta_zz_xxyy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxyy_1[j] + fl1_fx * ta_xzz_xyy_0[j] - fl1_fx * ta_xzz_xyy_1[j];

                ta_xxzz_xxyz_0[j] = pa_x[j] * ta_xzz_xxyz_0[j] - pc_x[j] * ta_xzz_xxyz_1[j] + 0.5 * fl1_fx * ta_zz_xxyz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxyz_1[j] + fl1_fx * ta_xzz_xyz_0[j] - fl1_fx * ta_xzz_xyz_1[j];

                ta_xxzz_xxzz_0[j] = pa_x[j] * ta_xzz_xxzz_0[j] - pc_x[j] * ta_xzz_xxzz_1[j] + 0.5 * fl1_fx * ta_zz_xxzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxzz_1[j] + fl1_fx * ta_xzz_xzz_0[j] - fl1_fx * ta_xzz_xzz_1[j];

                ta_xxzz_xyyy_0[j] = pa_x[j] * ta_xzz_xyyy_0[j] - pc_x[j] * ta_xzz_xyyy_1[j] + 0.5 * fl1_fx * ta_zz_xyyy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyyy_1[j] + 0.5 * fl1_fx * ta_xzz_yyy_0[j] - 0.5 * fl1_fx * ta_xzz_yyy_1[j];

                ta_xxzz_xyyz_0[j] = pa_x[j] * ta_xzz_xyyz_0[j] - pc_x[j] * ta_xzz_xyyz_1[j] + 0.5 * fl1_fx * ta_zz_xyyz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyyz_1[j] + 0.5 * fl1_fx * ta_xzz_yyz_0[j] - 0.5 * fl1_fx * ta_xzz_yyz_1[j];

                ta_xxzz_xyzz_0[j] = pa_x[j] * ta_xzz_xyzz_0[j] - pc_x[j] * ta_xzz_xyzz_1[j] + 0.5 * fl1_fx * ta_zz_xyzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyzz_1[j] + 0.5 * fl1_fx * ta_xzz_yzz_0[j] - 0.5 * fl1_fx * ta_xzz_yzz_1[j];

                ta_xxzz_xzzz_0[j] = pa_x[j] * ta_xzz_xzzz_0[j] - pc_x[j] * ta_xzz_xzzz_1[j] + 0.5 * fl1_fx * ta_zz_xzzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xzzz_1[j] + 0.5 * fl1_fx * ta_xzz_zzz_0[j] - 0.5 * fl1_fx * ta_xzz_zzz_1[j];

                ta_xxzz_yyyy_0[j] =
                    pa_x[j] * ta_xzz_yyyy_0[j] - pc_x[j] * ta_xzz_yyyy_1[j] + 0.5 * fl1_fx * ta_zz_yyyy_0[j] - 0.5 * fl1_fx * ta_zz_yyyy_1[j];

                ta_xxzz_yyyz_0[j] =
                    pa_x[j] * ta_xzz_yyyz_0[j] - pc_x[j] * ta_xzz_yyyz_1[j] + 0.5 * fl1_fx * ta_zz_yyyz_0[j] - 0.5 * fl1_fx * ta_zz_yyyz_1[j];

                ta_xxzz_yyzz_0[j] =
                    pa_x[j] * ta_xzz_yyzz_0[j] - pc_x[j] * ta_xzz_yyzz_1[j] + 0.5 * fl1_fx * ta_zz_yyzz_0[j] - 0.5 * fl1_fx * ta_zz_yyzz_1[j];

                ta_xxzz_yzzz_0[j] =
                    pa_x[j] * ta_xzz_yzzz_0[j] - pc_x[j] * ta_xzz_yzzz_1[j] + 0.5 * fl1_fx * ta_zz_yzzz_0[j] - 0.5 * fl1_fx * ta_zz_yzzz_1[j];

                ta_xxzz_zzzz_0[j] =
                    pa_x[j] * ta_xzz_zzzz_0[j] - pc_x[j] * ta_xzz_zzzz_1[j] + 0.5 * fl1_fx * ta_zz_zzzz_0[j] - 0.5 * fl1_fx * ta_zz_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGG_90_135(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
{
    // Batch of Integrals (90,135)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto ta_yyy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 90);

            auto ta_yyy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 91);

            auto ta_yyy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 92);

            auto ta_yyy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 93);

            auto ta_yyy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 94);

            auto ta_yyy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 95);

            auto ta_yyy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 96);

            auto ta_yyy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 97);

            auto ta_yyy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 98);

            auto ta_yyy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 99);

            auto ta_yyy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 100);

            auto ta_yyy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 101);

            auto ta_yyy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 102);

            auto ta_yyy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 103);

            auto ta_yyy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 104);

            auto ta_yyz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 105);

            auto ta_yyz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 106);

            auto ta_yyz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 107);

            auto ta_yyz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 108);

            auto ta_yyz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 109);

            auto ta_yyz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 110);

            auto ta_yyz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 111);

            auto ta_yyz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 112);

            auto ta_yyz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 113);

            auto ta_yyz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 114);

            auto ta_yyz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 115);

            auto ta_yyz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 116);

            auto ta_yyz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 117);

            auto ta_yyz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 118);

            auto ta_yyz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 119);

            auto ta_yzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 120);

            auto ta_yzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 121);

            auto ta_yzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 122);

            auto ta_yzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 123);

            auto ta_yzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 124);

            auto ta_yzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 125);

            auto ta_yzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 126);

            auto ta_yzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 127);

            auto ta_yzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 128);

            auto ta_yzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 129);

            auto ta_yzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 130);

            auto ta_yzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 131);

            auto ta_yzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 132);

            auto ta_yzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 133);

            auto ta_yzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 134);

            auto ta_yyy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 90);

            auto ta_yyy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 91);

            auto ta_yyy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 92);

            auto ta_yyy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 93);

            auto ta_yyy_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 94);

            auto ta_yyy_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 95);

            auto ta_yyy_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 96);

            auto ta_yyy_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 97);

            auto ta_yyy_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 98);

            auto ta_yyy_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 99);

            auto ta_yyy_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 100);

            auto ta_yyy_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 101);

            auto ta_yyy_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 102);

            auto ta_yyy_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 103);

            auto ta_yyy_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 104);

            auto ta_yyz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 105);

            auto ta_yyz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 106);

            auto ta_yyz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 107);

            auto ta_yyz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 108);

            auto ta_yyz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 109);

            auto ta_yyz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 110);

            auto ta_yyz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 111);

            auto ta_yyz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 112);

            auto ta_yyz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 113);

            auto ta_yyz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 114);

            auto ta_yyz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 115);

            auto ta_yyz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 116);

            auto ta_yyz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 117);

            auto ta_yyz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 118);

            auto ta_yyz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 119);

            auto ta_yzz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 120);

            auto ta_yzz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 121);

            auto ta_yzz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 122);

            auto ta_yzz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 123);

            auto ta_yzz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 124);

            auto ta_yzz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 125);

            auto ta_yzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 126);

            auto ta_yzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 127);

            auto ta_yzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 128);

            auto ta_yzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 129);

            auto ta_yzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 130);

            auto ta_yzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 131);

            auto ta_yzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 132);

            auto ta_yzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 133);

            auto ta_yzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 134);

            auto ta_yyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 60);

            auto ta_yyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 61);

            auto ta_yyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 62);

            auto ta_yyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 63);

            auto ta_yyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 64);

            auto ta_yyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 65);

            auto ta_yyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 66);

            auto ta_yyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 67);

            auto ta_yyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 68);

            auto ta_yyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 69);

            auto ta_yyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 70);

            auto ta_yyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 71);

            auto ta_yyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 72);

            auto ta_yyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 73);

            auto ta_yyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 74);

            auto ta_yyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 75);

            auto ta_yyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 76);

            auto ta_yyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 77);

            auto ta_yyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 78);

            auto ta_yyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 79);

            auto ta_yzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 80);

            auto ta_yzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 81);

            auto ta_yzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 82);

            auto ta_yzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 83);

            auto ta_yzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 84);

            auto ta_yzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 85);

            auto ta_yzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 86);

            auto ta_yzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 87);

            auto ta_yzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 88);

            auto ta_yzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 89);

            auto ta_yyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 60);

            auto ta_yyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 61);

            auto ta_yyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 62);

            auto ta_yyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 63);

            auto ta_yyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 64);

            auto ta_yyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 65);

            auto ta_yyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 66);

            auto ta_yyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 67);

            auto ta_yyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 68);

            auto ta_yyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 69);

            auto ta_yyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 70);

            auto ta_yyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 71);

            auto ta_yyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 72);

            auto ta_yyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 73);

            auto ta_yyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 74);

            auto ta_yyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 75);

            auto ta_yyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 76);

            auto ta_yyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 77);

            auto ta_yyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 78);

            auto ta_yyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 79);

            auto ta_yzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 80);

            auto ta_yzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 81);

            auto ta_yzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 82);

            auto ta_yzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 83);

            auto ta_yzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 84);

            auto ta_yzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 85);

            auto ta_yzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 86);

            auto ta_yzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 87);

            auto ta_yzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 88);

            auto ta_yzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 89);

            // set up pointers to integrals

            auto ta_xyyy_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 90);

            auto ta_xyyy_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 91);

            auto ta_xyyy_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 92);

            auto ta_xyyy_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 93);

            auto ta_xyyy_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 94);

            auto ta_xyyy_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 95);

            auto ta_xyyy_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 96);

            auto ta_xyyy_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 97);

            auto ta_xyyy_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 98);

            auto ta_xyyy_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 99);

            auto ta_xyyy_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 100);

            auto ta_xyyy_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 101);

            auto ta_xyyy_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 102);

            auto ta_xyyy_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 103);

            auto ta_xyyy_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 104);

            auto ta_xyyz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 105);

            auto ta_xyyz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 106);

            auto ta_xyyz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 107);

            auto ta_xyyz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 108);

            auto ta_xyyz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 109);

            auto ta_xyyz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 110);

            auto ta_xyyz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 111);

            auto ta_xyyz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 112);

            auto ta_xyyz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 113);

            auto ta_xyyz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 114);

            auto ta_xyyz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 115);

            auto ta_xyyz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 116);

            auto ta_xyyz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 117);

            auto ta_xyyz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 118);

            auto ta_xyyz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 119);

            auto ta_xyzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 120);

            auto ta_xyzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 121);

            auto ta_xyzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 122);

            auto ta_xyzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 123);

            auto ta_xyzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 124);

            auto ta_xyzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 125);

            auto ta_xyzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 126);

            auto ta_xyzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 127);

            auto ta_xyzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 128);

            auto ta_xyzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 129);

            auto ta_xyzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 130);

            auto ta_xyzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 131);

            auto ta_xyzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 132);

            auto ta_xyzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 133);

            auto ta_xyzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 134);

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xyyy_xxxx_0, ta_xyyy_xxxy_0, ta_xyyy_xxxz_0, \
                                         ta_xyyy_xxyy_0, ta_xyyy_xxyz_0, ta_xyyy_xxzz_0, ta_xyyy_xyyy_0, ta_xyyy_xyyz_0, \
                                         ta_xyyy_xyzz_0, ta_xyyy_xzzz_0, ta_xyyy_yyyy_0, ta_xyyy_yyyz_0, ta_xyyy_yyzz_0, \
                                         ta_xyyy_yzzz_0, ta_xyyy_zzzz_0, ta_xyyz_xxxx_0, ta_xyyz_xxxy_0, ta_xyyz_xxxz_0, \
                                         ta_xyyz_xxyy_0, ta_xyyz_xxyz_0, ta_xyyz_xxzz_0, ta_xyyz_xyyy_0, ta_xyyz_xyyz_0, \
                                         ta_xyyz_xyzz_0, ta_xyyz_xzzz_0, ta_xyyz_yyyy_0, ta_xyyz_yyyz_0, ta_xyyz_yyzz_0, \
                                         ta_xyyz_yzzz_0, ta_xyyz_zzzz_0, ta_xyzz_xxxx_0, ta_xyzz_xxxy_0, ta_xyzz_xxxz_0, \
                                         ta_xyzz_xxyy_0, ta_xyzz_xxyz_0, ta_xyzz_xxzz_0, ta_xyzz_xyyy_0, ta_xyzz_xyyz_0, \
                                         ta_xyzz_xyzz_0, ta_xyzz_xzzz_0, ta_xyzz_yyyy_0, ta_xyzz_yyyz_0, ta_xyzz_yyzz_0, \
                                         ta_xyzz_yzzz_0, ta_xyzz_zzzz_0, ta_yyy_xxx_0, ta_yyy_xxx_1, ta_yyy_xxxx_0, \
                                         ta_yyy_xxxx_1, ta_yyy_xxxy_0, ta_yyy_xxxy_1, ta_yyy_xxxz_0, ta_yyy_xxxz_1, \
                                         ta_yyy_xxy_0, ta_yyy_xxy_1, ta_yyy_xxyy_0, ta_yyy_xxyy_1, ta_yyy_xxyz_0, \
                                         ta_yyy_xxyz_1, ta_yyy_xxz_0, ta_yyy_xxz_1, ta_yyy_xxzz_0, ta_yyy_xxzz_1, \
                                         ta_yyy_xyy_0, ta_yyy_xyy_1, ta_yyy_xyyy_0, ta_yyy_xyyy_1, ta_yyy_xyyz_0, \
                                         ta_yyy_xyyz_1, ta_yyy_xyz_0, ta_yyy_xyz_1, ta_yyy_xyzz_0, ta_yyy_xyzz_1, \
                                         ta_yyy_xzz_0, ta_yyy_xzz_1, ta_yyy_xzzz_0, ta_yyy_xzzz_1, ta_yyy_yyy_0, \
                                         ta_yyy_yyy_1, ta_yyy_yyyy_0, ta_yyy_yyyy_1, ta_yyy_yyyz_0, ta_yyy_yyyz_1, \
                                         ta_yyy_yyz_0, ta_yyy_yyz_1, ta_yyy_yyzz_0, ta_yyy_yyzz_1, ta_yyy_yzz_0, \
                                         ta_yyy_yzz_1, ta_yyy_yzzz_0, ta_yyy_yzzz_1, ta_yyy_zzz_0, ta_yyy_zzz_1, \
                                         ta_yyy_zzzz_0, ta_yyy_zzzz_1, ta_yyz_xxx_0, ta_yyz_xxx_1, ta_yyz_xxxx_0, \
                                         ta_yyz_xxxx_1, ta_yyz_xxxy_0, ta_yyz_xxxy_1, ta_yyz_xxxz_0, ta_yyz_xxxz_1, \
                                         ta_yyz_xxy_0, ta_yyz_xxy_1, ta_yyz_xxyy_0, ta_yyz_xxyy_1, ta_yyz_xxyz_0, \
                                         ta_yyz_xxyz_1, ta_yyz_xxz_0, ta_yyz_xxz_1, ta_yyz_xxzz_0, ta_yyz_xxzz_1, \
                                         ta_yyz_xyy_0, ta_yyz_xyy_1, ta_yyz_xyyy_0, ta_yyz_xyyy_1, ta_yyz_xyyz_0, \
                                         ta_yyz_xyyz_1, ta_yyz_xyz_0, ta_yyz_xyz_1, ta_yyz_xyzz_0, ta_yyz_xyzz_1, \
                                         ta_yyz_xzz_0, ta_yyz_xzz_1, ta_yyz_xzzz_0, ta_yyz_xzzz_1, ta_yyz_yyy_0, \
                                         ta_yyz_yyy_1, ta_yyz_yyyy_0, ta_yyz_yyyy_1, ta_yyz_yyyz_0, ta_yyz_yyyz_1, \
                                         ta_yyz_yyz_0, ta_yyz_yyz_1, ta_yyz_yyzz_0, ta_yyz_yyzz_1, ta_yyz_yzz_0, \
                                         ta_yyz_yzz_1, ta_yyz_yzzz_0, ta_yyz_yzzz_1, ta_yyz_zzz_0, ta_yyz_zzz_1, \
                                         ta_yyz_zzzz_0, ta_yyz_zzzz_1, ta_yzz_xxx_0, ta_yzz_xxx_1, ta_yzz_xxxx_0, \
                                         ta_yzz_xxxx_1, ta_yzz_xxxy_0, ta_yzz_xxxy_1, ta_yzz_xxxz_0, ta_yzz_xxxz_1, \
                                         ta_yzz_xxy_0, ta_yzz_xxy_1, ta_yzz_xxyy_0, ta_yzz_xxyy_1, ta_yzz_xxyz_0, \
                                         ta_yzz_xxyz_1, ta_yzz_xxz_0, ta_yzz_xxz_1, ta_yzz_xxzz_0, ta_yzz_xxzz_1, \
                                         ta_yzz_xyy_0, ta_yzz_xyy_1, ta_yzz_xyyy_0, ta_yzz_xyyy_1, ta_yzz_xyyz_0, \
                                         ta_yzz_xyyz_1, ta_yzz_xyz_0, ta_yzz_xyz_1, ta_yzz_xyzz_0, ta_yzz_xyzz_1, \
                                         ta_yzz_xzz_0, ta_yzz_xzz_1, ta_yzz_xzzz_0, ta_yzz_xzzz_1, ta_yzz_yyy_0, \
                                         ta_yzz_yyy_1, ta_yzz_yyyy_0, ta_yzz_yyyy_1, ta_yzz_yyyz_0, ta_yzz_yyyz_1, \
                                         ta_yzz_yyz_0, ta_yzz_yyz_1, ta_yzz_yyzz_0, ta_yzz_yyzz_1, ta_yzz_yzz_0, \
                                         ta_yzz_yzz_1, ta_yzz_yzzz_0, ta_yzz_yzzz_1, ta_yzz_zzz_0, ta_yzz_zzz_1, \
                                         ta_yzz_zzzz_0, ta_yzz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xyyy_xxxx_0[j] =
                    pa_x[j] * ta_yyy_xxxx_0[j] - pc_x[j] * ta_yyy_xxxx_1[j] + 2.0 * fl1_fx * ta_yyy_xxx_0[j] - 2.0 * fl1_fx * ta_yyy_xxx_1[j];

                ta_xyyy_xxxy_0[j] =
                    pa_x[j] * ta_yyy_xxxy_0[j] - pc_x[j] * ta_yyy_xxxy_1[j] + 1.5 * fl1_fx * ta_yyy_xxy_0[j] - 1.5 * fl1_fx * ta_yyy_xxy_1[j];

                ta_xyyy_xxxz_0[j] =
                    pa_x[j] * ta_yyy_xxxz_0[j] - pc_x[j] * ta_yyy_xxxz_1[j] + 1.5 * fl1_fx * ta_yyy_xxz_0[j] - 1.5 * fl1_fx * ta_yyy_xxz_1[j];

                ta_xyyy_xxyy_0[j] = pa_x[j] * ta_yyy_xxyy_0[j] - pc_x[j] * ta_yyy_xxyy_1[j] + fl1_fx * ta_yyy_xyy_0[j] - fl1_fx * ta_yyy_xyy_1[j];

                ta_xyyy_xxyz_0[j] = pa_x[j] * ta_yyy_xxyz_0[j] - pc_x[j] * ta_yyy_xxyz_1[j] + fl1_fx * ta_yyy_xyz_0[j] - fl1_fx * ta_yyy_xyz_1[j];

                ta_xyyy_xxzz_0[j] = pa_x[j] * ta_yyy_xxzz_0[j] - pc_x[j] * ta_yyy_xxzz_1[j] + fl1_fx * ta_yyy_xzz_0[j] - fl1_fx * ta_yyy_xzz_1[j];

                ta_xyyy_xyyy_0[j] =
                    pa_x[j] * ta_yyy_xyyy_0[j] - pc_x[j] * ta_yyy_xyyy_1[j] + 0.5 * fl1_fx * ta_yyy_yyy_0[j] - 0.5 * fl1_fx * ta_yyy_yyy_1[j];

                ta_xyyy_xyyz_0[j] =
                    pa_x[j] * ta_yyy_xyyz_0[j] - pc_x[j] * ta_yyy_xyyz_1[j] + 0.5 * fl1_fx * ta_yyy_yyz_0[j] - 0.5 * fl1_fx * ta_yyy_yyz_1[j];

                ta_xyyy_xyzz_0[j] =
                    pa_x[j] * ta_yyy_xyzz_0[j] - pc_x[j] * ta_yyy_xyzz_1[j] + 0.5 * fl1_fx * ta_yyy_yzz_0[j] - 0.5 * fl1_fx * ta_yyy_yzz_1[j];

                ta_xyyy_xzzz_0[j] =
                    pa_x[j] * ta_yyy_xzzz_0[j] - pc_x[j] * ta_yyy_xzzz_1[j] + 0.5 * fl1_fx * ta_yyy_zzz_0[j] - 0.5 * fl1_fx * ta_yyy_zzz_1[j];

                ta_xyyy_yyyy_0[j] = pa_x[j] * ta_yyy_yyyy_0[j] - pc_x[j] * ta_yyy_yyyy_1[j];

                ta_xyyy_yyyz_0[j] = pa_x[j] * ta_yyy_yyyz_0[j] - pc_x[j] * ta_yyy_yyyz_1[j];

                ta_xyyy_yyzz_0[j] = pa_x[j] * ta_yyy_yyzz_0[j] - pc_x[j] * ta_yyy_yyzz_1[j];

                ta_xyyy_yzzz_0[j] = pa_x[j] * ta_yyy_yzzz_0[j] - pc_x[j] * ta_yyy_yzzz_1[j];

                ta_xyyy_zzzz_0[j] = pa_x[j] * ta_yyy_zzzz_0[j] - pc_x[j] * ta_yyy_zzzz_1[j];

                ta_xyyz_xxxx_0[j] =
                    pa_x[j] * ta_yyz_xxxx_0[j] - pc_x[j] * ta_yyz_xxxx_1[j] + 2.0 * fl1_fx * ta_yyz_xxx_0[j] - 2.0 * fl1_fx * ta_yyz_xxx_1[j];

                ta_xyyz_xxxy_0[j] =
                    pa_x[j] * ta_yyz_xxxy_0[j] - pc_x[j] * ta_yyz_xxxy_1[j] + 1.5 * fl1_fx * ta_yyz_xxy_0[j] - 1.5 * fl1_fx * ta_yyz_xxy_1[j];

                ta_xyyz_xxxz_0[j] =
                    pa_x[j] * ta_yyz_xxxz_0[j] - pc_x[j] * ta_yyz_xxxz_1[j] + 1.5 * fl1_fx * ta_yyz_xxz_0[j] - 1.5 * fl1_fx * ta_yyz_xxz_1[j];

                ta_xyyz_xxyy_0[j] = pa_x[j] * ta_yyz_xxyy_0[j] - pc_x[j] * ta_yyz_xxyy_1[j] + fl1_fx * ta_yyz_xyy_0[j] - fl1_fx * ta_yyz_xyy_1[j];

                ta_xyyz_xxyz_0[j] = pa_x[j] * ta_yyz_xxyz_0[j] - pc_x[j] * ta_yyz_xxyz_1[j] + fl1_fx * ta_yyz_xyz_0[j] - fl1_fx * ta_yyz_xyz_1[j];

                ta_xyyz_xxzz_0[j] = pa_x[j] * ta_yyz_xxzz_0[j] - pc_x[j] * ta_yyz_xxzz_1[j] + fl1_fx * ta_yyz_xzz_0[j] - fl1_fx * ta_yyz_xzz_1[j];

                ta_xyyz_xyyy_0[j] =
                    pa_x[j] * ta_yyz_xyyy_0[j] - pc_x[j] * ta_yyz_xyyy_1[j] + 0.5 * fl1_fx * ta_yyz_yyy_0[j] - 0.5 * fl1_fx * ta_yyz_yyy_1[j];

                ta_xyyz_xyyz_0[j] =
                    pa_x[j] * ta_yyz_xyyz_0[j] - pc_x[j] * ta_yyz_xyyz_1[j] + 0.5 * fl1_fx * ta_yyz_yyz_0[j] - 0.5 * fl1_fx * ta_yyz_yyz_1[j];

                ta_xyyz_xyzz_0[j] =
                    pa_x[j] * ta_yyz_xyzz_0[j] - pc_x[j] * ta_yyz_xyzz_1[j] + 0.5 * fl1_fx * ta_yyz_yzz_0[j] - 0.5 * fl1_fx * ta_yyz_yzz_1[j];

                ta_xyyz_xzzz_0[j] =
                    pa_x[j] * ta_yyz_xzzz_0[j] - pc_x[j] * ta_yyz_xzzz_1[j] + 0.5 * fl1_fx * ta_yyz_zzz_0[j] - 0.5 * fl1_fx * ta_yyz_zzz_1[j];

                ta_xyyz_yyyy_0[j] = pa_x[j] * ta_yyz_yyyy_0[j] - pc_x[j] * ta_yyz_yyyy_1[j];

                ta_xyyz_yyyz_0[j] = pa_x[j] * ta_yyz_yyyz_0[j] - pc_x[j] * ta_yyz_yyyz_1[j];

                ta_xyyz_yyzz_0[j] = pa_x[j] * ta_yyz_yyzz_0[j] - pc_x[j] * ta_yyz_yyzz_1[j];

                ta_xyyz_yzzz_0[j] = pa_x[j] * ta_yyz_yzzz_0[j] - pc_x[j] * ta_yyz_yzzz_1[j];

                ta_xyyz_zzzz_0[j] = pa_x[j] * ta_yyz_zzzz_0[j] - pc_x[j] * ta_yyz_zzzz_1[j];

                ta_xyzz_xxxx_0[j] =
                    pa_x[j] * ta_yzz_xxxx_0[j] - pc_x[j] * ta_yzz_xxxx_1[j] + 2.0 * fl1_fx * ta_yzz_xxx_0[j] - 2.0 * fl1_fx * ta_yzz_xxx_1[j];

                ta_xyzz_xxxy_0[j] =
                    pa_x[j] * ta_yzz_xxxy_0[j] - pc_x[j] * ta_yzz_xxxy_1[j] + 1.5 * fl1_fx * ta_yzz_xxy_0[j] - 1.5 * fl1_fx * ta_yzz_xxy_1[j];

                ta_xyzz_xxxz_0[j] =
                    pa_x[j] * ta_yzz_xxxz_0[j] - pc_x[j] * ta_yzz_xxxz_1[j] + 1.5 * fl1_fx * ta_yzz_xxz_0[j] - 1.5 * fl1_fx * ta_yzz_xxz_1[j];

                ta_xyzz_xxyy_0[j] = pa_x[j] * ta_yzz_xxyy_0[j] - pc_x[j] * ta_yzz_xxyy_1[j] + fl1_fx * ta_yzz_xyy_0[j] - fl1_fx * ta_yzz_xyy_1[j];

                ta_xyzz_xxyz_0[j] = pa_x[j] * ta_yzz_xxyz_0[j] - pc_x[j] * ta_yzz_xxyz_1[j] + fl1_fx * ta_yzz_xyz_0[j] - fl1_fx * ta_yzz_xyz_1[j];

                ta_xyzz_xxzz_0[j] = pa_x[j] * ta_yzz_xxzz_0[j] - pc_x[j] * ta_yzz_xxzz_1[j] + fl1_fx * ta_yzz_xzz_0[j] - fl1_fx * ta_yzz_xzz_1[j];

                ta_xyzz_xyyy_0[j] =
                    pa_x[j] * ta_yzz_xyyy_0[j] - pc_x[j] * ta_yzz_xyyy_1[j] + 0.5 * fl1_fx * ta_yzz_yyy_0[j] - 0.5 * fl1_fx * ta_yzz_yyy_1[j];

                ta_xyzz_xyyz_0[j] =
                    pa_x[j] * ta_yzz_xyyz_0[j] - pc_x[j] * ta_yzz_xyyz_1[j] + 0.5 * fl1_fx * ta_yzz_yyz_0[j] - 0.5 * fl1_fx * ta_yzz_yyz_1[j];

                ta_xyzz_xyzz_0[j] =
                    pa_x[j] * ta_yzz_xyzz_0[j] - pc_x[j] * ta_yzz_xyzz_1[j] + 0.5 * fl1_fx * ta_yzz_yzz_0[j] - 0.5 * fl1_fx * ta_yzz_yzz_1[j];

                ta_xyzz_xzzz_0[j] =
                    pa_x[j] * ta_yzz_xzzz_0[j] - pc_x[j] * ta_yzz_xzzz_1[j] + 0.5 * fl1_fx * ta_yzz_zzz_0[j] - 0.5 * fl1_fx * ta_yzz_zzz_1[j];

                ta_xyzz_yyyy_0[j] = pa_x[j] * ta_yzz_yyyy_0[j] - pc_x[j] * ta_yzz_yyyy_1[j];

                ta_xyzz_yyyz_0[j] = pa_x[j] * ta_yzz_yyyz_0[j] - pc_x[j] * ta_yzz_yyyz_1[j];

                ta_xyzz_yyzz_0[j] = pa_x[j] * ta_yzz_yyzz_0[j] - pc_x[j] * ta_yzz_yyzz_1[j];

                ta_xyzz_yzzz_0[j] = pa_x[j] * ta_yzz_yzzz_0[j] - pc_x[j] * ta_yzz_yzzz_1[j];

                ta_xyzz_zzzz_0[j] = pa_x[j] * ta_yzz_zzzz_0[j] - pc_x[j] * ta_yzz_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGG_135_180(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
{
    // Batch of Integrals (135,180)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto ta_yyy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 90);

            auto ta_yyy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 91);

            auto ta_yyy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 92);

            auto ta_yyy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 93);

            auto ta_yyy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 94);

            auto ta_yyy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 95);

            auto ta_yyy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 96);

            auto ta_yyy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 97);

            auto ta_yyy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 98);

            auto ta_yyy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 99);

            auto ta_yyy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 100);

            auto ta_yyy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 101);

            auto ta_yyy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 102);

            auto ta_yyy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 103);

            auto ta_yyy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 104);

            auto ta_yyz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 105);

            auto ta_yyz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 106);

            auto ta_yyz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 107);

            auto ta_yyz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 108);

            auto ta_yyz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 109);

            auto ta_yyz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 110);

            auto ta_yyz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 111);

            auto ta_yyz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 112);

            auto ta_yyz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 113);

            auto ta_yyz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 114);

            auto ta_yyz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 115);

            auto ta_yyz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 116);

            auto ta_yyz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 117);

            auto ta_yyz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 118);

            auto ta_yyz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 119);

            auto ta_zzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 135);

            auto ta_zzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 136);

            auto ta_zzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 137);

            auto ta_zzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 138);

            auto ta_zzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 139);

            auto ta_zzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 140);

            auto ta_zzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 141);

            auto ta_zzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 142);

            auto ta_zzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 143);

            auto ta_zzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 144);

            auto ta_zzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 145);

            auto ta_zzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 146);

            auto ta_zzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 147);

            auto ta_zzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 148);

            auto ta_zzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 149);

            auto ta_yyy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 90);

            auto ta_yyy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 91);

            auto ta_yyy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 92);

            auto ta_yyy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 93);

            auto ta_yyy_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 94);

            auto ta_yyy_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 95);

            auto ta_yyy_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 96);

            auto ta_yyy_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 97);

            auto ta_yyy_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 98);

            auto ta_yyy_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 99);

            auto ta_yyy_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 100);

            auto ta_yyy_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 101);

            auto ta_yyy_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 102);

            auto ta_yyy_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 103);

            auto ta_yyy_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 104);

            auto ta_yyz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 105);

            auto ta_yyz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 106);

            auto ta_yyz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 107);

            auto ta_yyz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 108);

            auto ta_yyz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 109);

            auto ta_yyz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 110);

            auto ta_yyz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 111);

            auto ta_yyz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 112);

            auto ta_yyz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 113);

            auto ta_yyz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 114);

            auto ta_yyz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 115);

            auto ta_yyz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 116);

            auto ta_yyz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 117);

            auto ta_yyz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 118);

            auto ta_yyz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 119);

            auto ta_zzz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 135);

            auto ta_zzz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 136);

            auto ta_zzz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 137);

            auto ta_zzz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 138);

            auto ta_zzz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 139);

            auto ta_zzz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 140);

            auto ta_zzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 141);

            auto ta_zzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 142);

            auto ta_zzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 143);

            auto ta_zzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 144);

            auto ta_zzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 145);

            auto ta_zzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 146);

            auto ta_zzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 147);

            auto ta_zzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 148);

            auto ta_zzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 149);

            auto ta_yy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 45);

            auto ta_yy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 46);

            auto ta_yy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 47);

            auto ta_yy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 48);

            auto ta_yy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 49);

            auto ta_yy_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 50);

            auto ta_yy_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 51);

            auto ta_yy_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 52);

            auto ta_yy_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 53);

            auto ta_yy_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 54);

            auto ta_yy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 55);

            auto ta_yy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 56);

            auto ta_yy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 57);

            auto ta_yy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 58);

            auto ta_yy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 59);

            auto ta_yz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 60);

            auto ta_yz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 61);

            auto ta_yz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 62);

            auto ta_yz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 63);

            auto ta_yz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 64);

            auto ta_yz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 65);

            auto ta_yz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 66);

            auto ta_yz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 67);

            auto ta_yz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 68);

            auto ta_yz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 69);

            auto ta_yz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 70);

            auto ta_yz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 71);

            auto ta_yz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 72);

            auto ta_yz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 73);

            auto ta_yz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 74);

            auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45);

            auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46);

            auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47);

            auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48);

            auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49);

            auto ta_yy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 50);

            auto ta_yy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 51);

            auto ta_yy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 52);

            auto ta_yy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 53);

            auto ta_yy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 54);

            auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55);

            auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56);

            auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57);

            auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58);

            auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59);

            auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60);

            auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61);

            auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62);

            auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63);

            auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64);

            auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65);

            auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66);

            auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67);

            auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68);

            auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69);

            auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70);

            auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71);

            auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72);

            auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73);

            auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74);

            auto ta_yyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 60);

            auto ta_yyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 61);

            auto ta_yyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 62);

            auto ta_yyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 63);

            auto ta_yyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 64);

            auto ta_yyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 65);

            auto ta_yyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 66);

            auto ta_yyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 67);

            auto ta_yyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 68);

            auto ta_yyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 69);

            auto ta_yyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 70);

            auto ta_yyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 71);

            auto ta_yyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 72);

            auto ta_yyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 73);

            auto ta_yyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 74);

            auto ta_yyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 75);

            auto ta_yyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 76);

            auto ta_yyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 77);

            auto ta_yyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 78);

            auto ta_yyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 79);

            auto ta_zzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 90);

            auto ta_zzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 91);

            auto ta_zzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 92);

            auto ta_zzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 93);

            auto ta_zzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 94);

            auto ta_zzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 95);

            auto ta_zzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 96);

            auto ta_zzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 97);

            auto ta_zzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 98);

            auto ta_zzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 99);

            auto ta_yyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 60);

            auto ta_yyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 61);

            auto ta_yyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 62);

            auto ta_yyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 63);

            auto ta_yyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 64);

            auto ta_yyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 65);

            auto ta_yyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 66);

            auto ta_yyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 67);

            auto ta_yyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 68);

            auto ta_yyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 69);

            auto ta_yyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 70);

            auto ta_yyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 71);

            auto ta_yyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 72);

            auto ta_yyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 73);

            auto ta_yyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 74);

            auto ta_yyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 75);

            auto ta_yyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 76);

            auto ta_yyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 77);

            auto ta_yyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 78);

            auto ta_yyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 79);

            auto ta_zzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 90);

            auto ta_zzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 91);

            auto ta_zzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 92);

            auto ta_zzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 93);

            auto ta_zzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 94);

            auto ta_zzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 95);

            auto ta_zzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 96);

            auto ta_zzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 97);

            auto ta_zzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 98);

            auto ta_zzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 99);

            // set up pointers to integrals

            auto ta_xzzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 135);

            auto ta_xzzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 136);

            auto ta_xzzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 137);

            auto ta_xzzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 138);

            auto ta_xzzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 139);

            auto ta_xzzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 140);

            auto ta_xzzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 141);

            auto ta_xzzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 142);

            auto ta_xzzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 143);

            auto ta_xzzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 144);

            auto ta_xzzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 145);

            auto ta_xzzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 146);

            auto ta_xzzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 147);

            auto ta_xzzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 148);

            auto ta_xzzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 149);

            auto ta_yyyy_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 150);

            auto ta_yyyy_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 151);

            auto ta_yyyy_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 152);

            auto ta_yyyy_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 153);

            auto ta_yyyy_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 154);

            auto ta_yyyy_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 155);

            auto ta_yyyy_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 156);

            auto ta_yyyy_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 157);

            auto ta_yyyy_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 158);

            auto ta_yyyy_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 159);

            auto ta_yyyy_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 160);

            auto ta_yyyy_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 161);

            auto ta_yyyy_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 162);

            auto ta_yyyy_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 163);

            auto ta_yyyy_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 164);

            auto ta_yyyz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 165);

            auto ta_yyyz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 166);

            auto ta_yyyz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 167);

            auto ta_yyyz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 168);

            auto ta_yyyz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 169);

            auto ta_yyyz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 170);

            auto ta_yyyz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 171);

            auto ta_yyyz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 172);

            auto ta_yyyz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 173);

            auto ta_yyyz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 174);

            auto ta_yyyz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 175);

            auto ta_yyyz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 176);

            auto ta_yyyz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 177);

            auto ta_yyyz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 178);

            auto ta_yyyz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 179);

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_xzzz_xxxx_0, ta_xzzz_xxxy_0, ta_xzzz_xxxz_0, \
                                         ta_xzzz_xxyy_0, ta_xzzz_xxyz_0, ta_xzzz_xxzz_0, ta_xzzz_xyyy_0, ta_xzzz_xyyz_0, \
                                         ta_xzzz_xyzz_0, ta_xzzz_xzzz_0, ta_xzzz_yyyy_0, ta_xzzz_yyyz_0, ta_xzzz_yyzz_0, \
                                         ta_xzzz_yzzz_0, ta_xzzz_zzzz_0, ta_yy_xxxx_0, ta_yy_xxxx_1, ta_yy_xxxy_0, \
                                         ta_yy_xxxy_1, ta_yy_xxxz_0, ta_yy_xxxz_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyz_0, \
                                         ta_yy_xxyz_1, ta_yy_xxzz_0, ta_yy_xxzz_1, ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyz_0, \
                                         ta_yy_xyyz_1, ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xzzz_0, ta_yy_xzzz_1, ta_yy_yyyy_0, \
                                         ta_yy_yyyy_1, ta_yy_yyyz_0, ta_yy_yyyz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yzzz_0, \
                                         ta_yy_yzzz_1, ta_yy_zzzz_0, ta_yy_zzzz_1, ta_yyy_xxx_0, ta_yyy_xxx_1, \
                                         ta_yyy_xxxx_0, ta_yyy_xxxx_1, ta_yyy_xxxy_0, ta_yyy_xxxy_1, ta_yyy_xxxz_0, \
                                         ta_yyy_xxxz_1, ta_yyy_xxy_0, ta_yyy_xxy_1, ta_yyy_xxyy_0, ta_yyy_xxyy_1, \
                                         ta_yyy_xxyz_0, ta_yyy_xxyz_1, ta_yyy_xxz_0, ta_yyy_xxz_1, ta_yyy_xxzz_0, \
                                         ta_yyy_xxzz_1, ta_yyy_xyy_0, ta_yyy_xyy_1, ta_yyy_xyyy_0, ta_yyy_xyyy_1, \
                                         ta_yyy_xyyz_0, ta_yyy_xyyz_1, ta_yyy_xyz_0, ta_yyy_xyz_1, ta_yyy_xyzz_0, \
                                         ta_yyy_xyzz_1, ta_yyy_xzz_0, ta_yyy_xzz_1, ta_yyy_xzzz_0, ta_yyy_xzzz_1, \
                                         ta_yyy_yyy_0, ta_yyy_yyy_1, ta_yyy_yyyy_0, ta_yyy_yyyy_1, ta_yyy_yyyz_0, \
                                         ta_yyy_yyyz_1, ta_yyy_yyz_0, ta_yyy_yyz_1, ta_yyy_yyzz_0, ta_yyy_yyzz_1, \
                                         ta_yyy_yzz_0, ta_yyy_yzz_1, ta_yyy_yzzz_0, ta_yyy_yzzz_1, ta_yyy_zzz_0, \
                                         ta_yyy_zzz_1, ta_yyy_zzzz_0, ta_yyy_zzzz_1, ta_yyyy_xxxx_0, ta_yyyy_xxxy_0, \
                                         ta_yyyy_xxxz_0, ta_yyyy_xxyy_0, ta_yyyy_xxyz_0, ta_yyyy_xxzz_0, ta_yyyy_xyyy_0, \
                                         ta_yyyy_xyyz_0, ta_yyyy_xyzz_0, ta_yyyy_xzzz_0, ta_yyyy_yyyy_0, ta_yyyy_yyyz_0, \
                                         ta_yyyy_yyzz_0, ta_yyyy_yzzz_0, ta_yyyy_zzzz_0, ta_yyyz_xxxx_0, ta_yyyz_xxxy_0, \
                                         ta_yyyz_xxxz_0, ta_yyyz_xxyy_0, ta_yyyz_xxyz_0, ta_yyyz_xxzz_0, ta_yyyz_xyyy_0, \
                                         ta_yyyz_xyyz_0, ta_yyyz_xyzz_0, ta_yyyz_xzzz_0, ta_yyyz_yyyy_0, ta_yyyz_yyyz_0, \
                                         ta_yyyz_yyzz_0, ta_yyyz_yzzz_0, ta_yyyz_zzzz_0, ta_yyz_xxx_0, ta_yyz_xxx_1, \
                                         ta_yyz_xxxx_0, ta_yyz_xxxx_1, ta_yyz_xxxy_0, ta_yyz_xxxy_1, ta_yyz_xxxz_0, \
                                         ta_yyz_xxxz_1, ta_yyz_xxy_0, ta_yyz_xxy_1, ta_yyz_xxyy_0, ta_yyz_xxyy_1, \
                                         ta_yyz_xxyz_0, ta_yyz_xxyz_1, ta_yyz_xxz_0, ta_yyz_xxz_1, ta_yyz_xxzz_0, \
                                         ta_yyz_xxzz_1, ta_yyz_xyy_0, ta_yyz_xyy_1, ta_yyz_xyyy_0, ta_yyz_xyyy_1, \
                                         ta_yyz_xyyz_0, ta_yyz_xyyz_1, ta_yyz_xyz_0, ta_yyz_xyz_1, ta_yyz_xyzz_0, \
                                         ta_yyz_xyzz_1, ta_yyz_xzz_0, ta_yyz_xzz_1, ta_yyz_xzzz_0, ta_yyz_xzzz_1, \
                                         ta_yyz_yyy_0, ta_yyz_yyy_1, ta_yyz_yyyy_0, ta_yyz_yyyy_1, ta_yyz_yyyz_0, \
                                         ta_yyz_yyyz_1, ta_yyz_yyz_0, ta_yyz_yyz_1, ta_yyz_yyzz_0, ta_yyz_yyzz_1, \
                                         ta_yyz_yzz_0, ta_yyz_yzz_1, ta_yyz_yzzz_0, ta_yyz_yzzz_1, ta_yyz_zzz_0, \
                                         ta_yyz_zzz_1, ta_yyz_zzzz_0, ta_yyz_zzzz_1, ta_yz_xxxx_0, ta_yz_xxxx_1, \
                                         ta_yz_xxxy_0, ta_yz_xxxy_1, ta_yz_xxxz_0, ta_yz_xxxz_1, ta_yz_xxyy_0, ta_yz_xxyy_1, \
                                         ta_yz_xxyz_0, ta_yz_xxyz_1, ta_yz_xxzz_0, ta_yz_xxzz_1, ta_yz_xyyy_0, ta_yz_xyyy_1, \
                                         ta_yz_xyyz_0, ta_yz_xyyz_1, ta_yz_xyzz_0, ta_yz_xyzz_1, ta_yz_xzzz_0, ta_yz_xzzz_1, \
                                         ta_yz_yyyy_0, ta_yz_yyyy_1, ta_yz_yyyz_0, ta_yz_yyyz_1, ta_yz_yyzz_0, ta_yz_yyzz_1, \
                                         ta_yz_yzzz_0, ta_yz_yzzz_1, ta_yz_zzzz_0, ta_yz_zzzz_1, ta_zzz_xxx_0, ta_zzz_xxx_1, \
                                         ta_zzz_xxxx_0, ta_zzz_xxxx_1, ta_zzz_xxxy_0, ta_zzz_xxxy_1, ta_zzz_xxxz_0, \
                                         ta_zzz_xxxz_1, ta_zzz_xxy_0, ta_zzz_xxy_1, ta_zzz_xxyy_0, ta_zzz_xxyy_1, \
                                         ta_zzz_xxyz_0, ta_zzz_xxyz_1, ta_zzz_xxz_0, ta_zzz_xxz_1, ta_zzz_xxzz_0, \
                                         ta_zzz_xxzz_1, ta_zzz_xyy_0, ta_zzz_xyy_1, ta_zzz_xyyy_0, ta_zzz_xyyy_1, \
                                         ta_zzz_xyyz_0, ta_zzz_xyyz_1, ta_zzz_xyz_0, ta_zzz_xyz_1, ta_zzz_xyzz_0, \
                                         ta_zzz_xyzz_1, ta_zzz_xzz_0, ta_zzz_xzz_1, ta_zzz_xzzz_0, ta_zzz_xzzz_1, \
                                         ta_zzz_yyy_0, ta_zzz_yyy_1, ta_zzz_yyyy_0, ta_zzz_yyyy_1, ta_zzz_yyyz_0, \
                                         ta_zzz_yyyz_1, ta_zzz_yyz_0, ta_zzz_yyz_1, ta_zzz_yyzz_0, ta_zzz_yyzz_1, \
                                         ta_zzz_yzz_0, ta_zzz_yzz_1, ta_zzz_yzzz_0, ta_zzz_yzzz_1, ta_zzz_zzz_0, \
                                         ta_zzz_zzz_1, ta_zzz_zzzz_0, ta_zzz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xzzz_xxxx_0[j] =
                    pa_x[j] * ta_zzz_xxxx_0[j] - pc_x[j] * ta_zzz_xxxx_1[j] + 2.0 * fl1_fx * ta_zzz_xxx_0[j] - 2.0 * fl1_fx * ta_zzz_xxx_1[j];

                ta_xzzz_xxxy_0[j] =
                    pa_x[j] * ta_zzz_xxxy_0[j] - pc_x[j] * ta_zzz_xxxy_1[j] + 1.5 * fl1_fx * ta_zzz_xxy_0[j] - 1.5 * fl1_fx * ta_zzz_xxy_1[j];

                ta_xzzz_xxxz_0[j] =
                    pa_x[j] * ta_zzz_xxxz_0[j] - pc_x[j] * ta_zzz_xxxz_1[j] + 1.5 * fl1_fx * ta_zzz_xxz_0[j] - 1.5 * fl1_fx * ta_zzz_xxz_1[j];

                ta_xzzz_xxyy_0[j] = pa_x[j] * ta_zzz_xxyy_0[j] - pc_x[j] * ta_zzz_xxyy_1[j] + fl1_fx * ta_zzz_xyy_0[j] - fl1_fx * ta_zzz_xyy_1[j];

                ta_xzzz_xxyz_0[j] = pa_x[j] * ta_zzz_xxyz_0[j] - pc_x[j] * ta_zzz_xxyz_1[j] + fl1_fx * ta_zzz_xyz_0[j] - fl1_fx * ta_zzz_xyz_1[j];

                ta_xzzz_xxzz_0[j] = pa_x[j] * ta_zzz_xxzz_0[j] - pc_x[j] * ta_zzz_xxzz_1[j] + fl1_fx * ta_zzz_xzz_0[j] - fl1_fx * ta_zzz_xzz_1[j];

                ta_xzzz_xyyy_0[j] =
                    pa_x[j] * ta_zzz_xyyy_0[j] - pc_x[j] * ta_zzz_xyyy_1[j] + 0.5 * fl1_fx * ta_zzz_yyy_0[j] - 0.5 * fl1_fx * ta_zzz_yyy_1[j];

                ta_xzzz_xyyz_0[j] =
                    pa_x[j] * ta_zzz_xyyz_0[j] - pc_x[j] * ta_zzz_xyyz_1[j] + 0.5 * fl1_fx * ta_zzz_yyz_0[j] - 0.5 * fl1_fx * ta_zzz_yyz_1[j];

                ta_xzzz_xyzz_0[j] =
                    pa_x[j] * ta_zzz_xyzz_0[j] - pc_x[j] * ta_zzz_xyzz_1[j] + 0.5 * fl1_fx * ta_zzz_yzz_0[j] - 0.5 * fl1_fx * ta_zzz_yzz_1[j];

                ta_xzzz_xzzz_0[j] =
                    pa_x[j] * ta_zzz_xzzz_0[j] - pc_x[j] * ta_zzz_xzzz_1[j] + 0.5 * fl1_fx * ta_zzz_zzz_0[j] - 0.5 * fl1_fx * ta_zzz_zzz_1[j];

                ta_xzzz_yyyy_0[j] = pa_x[j] * ta_zzz_yyyy_0[j] - pc_x[j] * ta_zzz_yyyy_1[j];

                ta_xzzz_yyyz_0[j] = pa_x[j] * ta_zzz_yyyz_0[j] - pc_x[j] * ta_zzz_yyyz_1[j];

                ta_xzzz_yyzz_0[j] = pa_x[j] * ta_zzz_yyzz_0[j] - pc_x[j] * ta_zzz_yyzz_1[j];

                ta_xzzz_yzzz_0[j] = pa_x[j] * ta_zzz_yzzz_0[j] - pc_x[j] * ta_zzz_yzzz_1[j];

                ta_xzzz_zzzz_0[j] = pa_x[j] * ta_zzz_zzzz_0[j] - pc_x[j] * ta_zzz_zzzz_1[j];

                ta_yyyy_xxxx_0[j] =
                    pa_y[j] * ta_yyy_xxxx_0[j] - pc_y[j] * ta_yyy_xxxx_1[j] + 1.5 * fl1_fx * ta_yy_xxxx_0[j] - 1.5 * fl1_fx * ta_yy_xxxx_1[j];

                ta_yyyy_xxxy_0[j] = pa_y[j] * ta_yyy_xxxy_0[j] - pc_y[j] * ta_yyy_xxxy_1[j] + 1.5 * fl1_fx * ta_yy_xxxy_0[j] -
                                    1.5 * fl1_fx * ta_yy_xxxy_1[j] + 0.5 * fl1_fx * ta_yyy_xxx_0[j] - 0.5 * fl1_fx * ta_yyy_xxx_1[j];

                ta_yyyy_xxxz_0[j] =
                    pa_y[j] * ta_yyy_xxxz_0[j] - pc_y[j] * ta_yyy_xxxz_1[j] + 1.5 * fl1_fx * ta_yy_xxxz_0[j] - 1.5 * fl1_fx * ta_yy_xxxz_1[j];

                ta_yyyy_xxyy_0[j] = pa_y[j] * ta_yyy_xxyy_0[j] - pc_y[j] * ta_yyy_xxyy_1[j] + 1.5 * fl1_fx * ta_yy_xxyy_0[j] -
                                    1.5 * fl1_fx * ta_yy_xxyy_1[j] + fl1_fx * ta_yyy_xxy_0[j] - fl1_fx * ta_yyy_xxy_1[j];

                ta_yyyy_xxyz_0[j] = pa_y[j] * ta_yyy_xxyz_0[j] - pc_y[j] * ta_yyy_xxyz_1[j] + 1.5 * fl1_fx * ta_yy_xxyz_0[j] -
                                    1.5 * fl1_fx * ta_yy_xxyz_1[j] + 0.5 * fl1_fx * ta_yyy_xxz_0[j] - 0.5 * fl1_fx * ta_yyy_xxz_1[j];

                ta_yyyy_xxzz_0[j] =
                    pa_y[j] * ta_yyy_xxzz_0[j] - pc_y[j] * ta_yyy_xxzz_1[j] + 1.5 * fl1_fx * ta_yy_xxzz_0[j] - 1.5 * fl1_fx * ta_yy_xxzz_1[j];

                ta_yyyy_xyyy_0[j] = pa_y[j] * ta_yyy_xyyy_0[j] - pc_y[j] * ta_yyy_xyyy_1[j] + 1.5 * fl1_fx * ta_yy_xyyy_0[j] -
                                    1.5 * fl1_fx * ta_yy_xyyy_1[j] + 1.5 * fl1_fx * ta_yyy_xyy_0[j] - 1.5 * fl1_fx * ta_yyy_xyy_1[j];

                ta_yyyy_xyyz_0[j] = pa_y[j] * ta_yyy_xyyz_0[j] - pc_y[j] * ta_yyy_xyyz_1[j] + 1.5 * fl1_fx * ta_yy_xyyz_0[j] -
                                    1.5 * fl1_fx * ta_yy_xyyz_1[j] + fl1_fx * ta_yyy_xyz_0[j] - fl1_fx * ta_yyy_xyz_1[j];

                ta_yyyy_xyzz_0[j] = pa_y[j] * ta_yyy_xyzz_0[j] - pc_y[j] * ta_yyy_xyzz_1[j] + 1.5 * fl1_fx * ta_yy_xyzz_0[j] -
                                    1.5 * fl1_fx * ta_yy_xyzz_1[j] + 0.5 * fl1_fx * ta_yyy_xzz_0[j] - 0.5 * fl1_fx * ta_yyy_xzz_1[j];

                ta_yyyy_xzzz_0[j] =
                    pa_y[j] * ta_yyy_xzzz_0[j] - pc_y[j] * ta_yyy_xzzz_1[j] + 1.5 * fl1_fx * ta_yy_xzzz_0[j] - 1.5 * fl1_fx * ta_yy_xzzz_1[j];

                ta_yyyy_yyyy_0[j] = pa_y[j] * ta_yyy_yyyy_0[j] - pc_y[j] * ta_yyy_yyyy_1[j] + 1.5 * fl1_fx * ta_yy_yyyy_0[j] -
                                    1.5 * fl1_fx * ta_yy_yyyy_1[j] + 2.0 * fl1_fx * ta_yyy_yyy_0[j] - 2.0 * fl1_fx * ta_yyy_yyy_1[j];

                ta_yyyy_yyyz_0[j] = pa_y[j] * ta_yyy_yyyz_0[j] - pc_y[j] * ta_yyy_yyyz_1[j] + 1.5 * fl1_fx * ta_yy_yyyz_0[j] -
                                    1.5 * fl1_fx * ta_yy_yyyz_1[j] + 1.5 * fl1_fx * ta_yyy_yyz_0[j] - 1.5 * fl1_fx * ta_yyy_yyz_1[j];

                ta_yyyy_yyzz_0[j] = pa_y[j] * ta_yyy_yyzz_0[j] - pc_y[j] * ta_yyy_yyzz_1[j] + 1.5 * fl1_fx * ta_yy_yyzz_0[j] -
                                    1.5 * fl1_fx * ta_yy_yyzz_1[j] + fl1_fx * ta_yyy_yzz_0[j] - fl1_fx * ta_yyy_yzz_1[j];

                ta_yyyy_yzzz_0[j] = pa_y[j] * ta_yyy_yzzz_0[j] - pc_y[j] * ta_yyy_yzzz_1[j] + 1.5 * fl1_fx * ta_yy_yzzz_0[j] -
                                    1.5 * fl1_fx * ta_yy_yzzz_1[j] + 0.5 * fl1_fx * ta_yyy_zzz_0[j] - 0.5 * fl1_fx * ta_yyy_zzz_1[j];

                ta_yyyy_zzzz_0[j] =
                    pa_y[j] * ta_yyy_zzzz_0[j] - pc_y[j] * ta_yyy_zzzz_1[j] + 1.5 * fl1_fx * ta_yy_zzzz_0[j] - 1.5 * fl1_fx * ta_yy_zzzz_1[j];

                ta_yyyz_xxxx_0[j] = pa_y[j] * ta_yyz_xxxx_0[j] - pc_y[j] * ta_yyz_xxxx_1[j] + fl1_fx * ta_yz_xxxx_0[j] - fl1_fx * ta_yz_xxxx_1[j];

                ta_yyyz_xxxy_0[j] = pa_y[j] * ta_yyz_xxxy_0[j] - pc_y[j] * ta_yyz_xxxy_1[j] + fl1_fx * ta_yz_xxxy_0[j] - fl1_fx * ta_yz_xxxy_1[j] +
                                    0.5 * fl1_fx * ta_yyz_xxx_0[j] - 0.5 * fl1_fx * ta_yyz_xxx_1[j];

                ta_yyyz_xxxz_0[j] = pa_y[j] * ta_yyz_xxxz_0[j] - pc_y[j] * ta_yyz_xxxz_1[j] + fl1_fx * ta_yz_xxxz_0[j] - fl1_fx * ta_yz_xxxz_1[j];

                ta_yyyz_xxyy_0[j] = pa_y[j] * ta_yyz_xxyy_0[j] - pc_y[j] * ta_yyz_xxyy_1[j] + fl1_fx * ta_yz_xxyy_0[j] - fl1_fx * ta_yz_xxyy_1[j] +
                                    fl1_fx * ta_yyz_xxy_0[j] - fl1_fx * ta_yyz_xxy_1[j];

                ta_yyyz_xxyz_0[j] = pa_y[j] * ta_yyz_xxyz_0[j] - pc_y[j] * ta_yyz_xxyz_1[j] + fl1_fx * ta_yz_xxyz_0[j] - fl1_fx * ta_yz_xxyz_1[j] +
                                    0.5 * fl1_fx * ta_yyz_xxz_0[j] - 0.5 * fl1_fx * ta_yyz_xxz_1[j];

                ta_yyyz_xxzz_0[j] = pa_y[j] * ta_yyz_xxzz_0[j] - pc_y[j] * ta_yyz_xxzz_1[j] + fl1_fx * ta_yz_xxzz_0[j] - fl1_fx * ta_yz_xxzz_1[j];

                ta_yyyz_xyyy_0[j] = pa_y[j] * ta_yyz_xyyy_0[j] - pc_y[j] * ta_yyz_xyyy_1[j] + fl1_fx * ta_yz_xyyy_0[j] - fl1_fx * ta_yz_xyyy_1[j] +
                                    1.5 * fl1_fx * ta_yyz_xyy_0[j] - 1.5 * fl1_fx * ta_yyz_xyy_1[j];

                ta_yyyz_xyyz_0[j] = pa_y[j] * ta_yyz_xyyz_0[j] - pc_y[j] * ta_yyz_xyyz_1[j] + fl1_fx * ta_yz_xyyz_0[j] - fl1_fx * ta_yz_xyyz_1[j] +
                                    fl1_fx * ta_yyz_xyz_0[j] - fl1_fx * ta_yyz_xyz_1[j];

                ta_yyyz_xyzz_0[j] = pa_y[j] * ta_yyz_xyzz_0[j] - pc_y[j] * ta_yyz_xyzz_1[j] + fl1_fx * ta_yz_xyzz_0[j] - fl1_fx * ta_yz_xyzz_1[j] +
                                    0.5 * fl1_fx * ta_yyz_xzz_0[j] - 0.5 * fl1_fx * ta_yyz_xzz_1[j];

                ta_yyyz_xzzz_0[j] = pa_y[j] * ta_yyz_xzzz_0[j] - pc_y[j] * ta_yyz_xzzz_1[j] + fl1_fx * ta_yz_xzzz_0[j] - fl1_fx * ta_yz_xzzz_1[j];

                ta_yyyz_yyyy_0[j] = pa_y[j] * ta_yyz_yyyy_0[j] - pc_y[j] * ta_yyz_yyyy_1[j] + fl1_fx * ta_yz_yyyy_0[j] - fl1_fx * ta_yz_yyyy_1[j] +
                                    2.0 * fl1_fx * ta_yyz_yyy_0[j] - 2.0 * fl1_fx * ta_yyz_yyy_1[j];

                ta_yyyz_yyyz_0[j] = pa_y[j] * ta_yyz_yyyz_0[j] - pc_y[j] * ta_yyz_yyyz_1[j] + fl1_fx * ta_yz_yyyz_0[j] - fl1_fx * ta_yz_yyyz_1[j] +
                                    1.5 * fl1_fx * ta_yyz_yyz_0[j] - 1.5 * fl1_fx * ta_yyz_yyz_1[j];

                ta_yyyz_yyzz_0[j] = pa_y[j] * ta_yyz_yyzz_0[j] - pc_y[j] * ta_yyz_yyzz_1[j] + fl1_fx * ta_yz_yyzz_0[j] - fl1_fx * ta_yz_yyzz_1[j] +
                                    fl1_fx * ta_yyz_yzz_0[j] - fl1_fx * ta_yyz_yzz_1[j];

                ta_yyyz_yzzz_0[j] = pa_y[j] * ta_yyz_yzzz_0[j] - pc_y[j] * ta_yyz_yzzz_1[j] + fl1_fx * ta_yz_yzzz_0[j] - fl1_fx * ta_yz_yzzz_1[j] +
                                    0.5 * fl1_fx * ta_yyz_zzz_0[j] - 0.5 * fl1_fx * ta_yyz_zzz_1[j];

                ta_yyyz_zzzz_0[j] = pa_y[j] * ta_yyz_zzzz_0[j] - pc_y[j] * ta_yyz_zzzz_1[j] + fl1_fx * ta_yz_zzzz_0[j] - fl1_fx * ta_yz_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGG_180_225(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
{
    // Batch of Integrals (180,225)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_yzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 120);

            auto ta_yzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 121);

            auto ta_yzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 122);

            auto ta_yzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 123);

            auto ta_yzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 124);

            auto ta_yzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 125);

            auto ta_yzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 126);

            auto ta_yzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 127);

            auto ta_yzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 128);

            auto ta_yzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 129);

            auto ta_yzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 130);

            auto ta_yzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 131);

            auto ta_yzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 132);

            auto ta_yzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 133);

            auto ta_yzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 134);

            auto ta_zzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 135);

            auto ta_zzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 136);

            auto ta_zzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 137);

            auto ta_zzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 138);

            auto ta_zzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 139);

            auto ta_zzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 140);

            auto ta_zzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 141);

            auto ta_zzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 142);

            auto ta_zzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 143);

            auto ta_zzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 144);

            auto ta_zzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 145);

            auto ta_zzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 146);

            auto ta_zzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 147);

            auto ta_zzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 148);

            auto ta_zzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 149);

            auto ta_yzz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 120);

            auto ta_yzz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 121);

            auto ta_yzz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 122);

            auto ta_yzz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 123);

            auto ta_yzz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 124);

            auto ta_yzz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 125);

            auto ta_yzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 126);

            auto ta_yzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 127);

            auto ta_yzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 128);

            auto ta_yzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 129);

            auto ta_yzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 130);

            auto ta_yzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 131);

            auto ta_yzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 132);

            auto ta_yzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 133);

            auto ta_yzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 134);

            auto ta_zzz_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 135);

            auto ta_zzz_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 136);

            auto ta_zzz_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 137);

            auto ta_zzz_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 138);

            auto ta_zzz_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 139);

            auto ta_zzz_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 140);

            auto ta_zzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 141);

            auto ta_zzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 142);

            auto ta_zzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 143);

            auto ta_zzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 144);

            auto ta_zzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 145);

            auto ta_zzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 146);

            auto ta_zzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 147);

            auto ta_zzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 148);

            auto ta_zzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 149);

            auto ta_zz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 75);

            auto ta_zz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 76);

            auto ta_zz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 77);

            auto ta_zz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 78);

            auto ta_zz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 79);

            auto ta_zz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 80);

            auto ta_zz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 81);

            auto ta_zz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 82);

            auto ta_zz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 83);

            auto ta_zz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 84);

            auto ta_zz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 85);

            auto ta_zz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 86);

            auto ta_zz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 87);

            auto ta_zz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 88);

            auto ta_zz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 89);

            auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75);

            auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76);

            auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77);

            auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78);

            auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79);

            auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80);

            auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81);

            auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82);

            auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83);

            auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84);

            auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85);

            auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86);

            auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87);

            auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88);

            auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89);

            auto ta_yzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 80);

            auto ta_yzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 81);

            auto ta_yzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 82);

            auto ta_yzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 83);

            auto ta_yzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 84);

            auto ta_yzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 85);

            auto ta_yzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 86);

            auto ta_yzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 87);

            auto ta_yzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 88);

            auto ta_yzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 89);

            auto ta_zzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 90);

            auto ta_zzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 91);

            auto ta_zzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 92);

            auto ta_zzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 93);

            auto ta_zzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 94);

            auto ta_zzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 95);

            auto ta_zzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 96);

            auto ta_zzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 97);

            auto ta_zzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 98);

            auto ta_zzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 99);

            auto ta_yzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 80);

            auto ta_yzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 81);

            auto ta_yzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 82);

            auto ta_yzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 83);

            auto ta_yzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 84);

            auto ta_yzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 85);

            auto ta_yzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 86);

            auto ta_yzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 87);

            auto ta_yzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 88);

            auto ta_yzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 89);

            auto ta_zzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 90);

            auto ta_zzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 91);

            auto ta_zzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 92);

            auto ta_zzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 93);

            auto ta_zzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 94);

            auto ta_zzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 95);

            auto ta_zzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 96);

            auto ta_zzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 97);

            auto ta_zzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 98);

            auto ta_zzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 99);

            // set up pointers to integrals

            auto ta_yyzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 180);

            auto ta_yyzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 181);

            auto ta_yyzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 182);

            auto ta_yyzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 183);

            auto ta_yyzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 184);

            auto ta_yyzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 185);

            auto ta_yyzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 186);

            auto ta_yyzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 187);

            auto ta_yyzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 188);

            auto ta_yyzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 189);

            auto ta_yyzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 190);

            auto ta_yyzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 191);

            auto ta_yyzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 192);

            auto ta_yyzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 193);

            auto ta_yyzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 194);

            auto ta_yzzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 195);

            auto ta_yzzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 196);

            auto ta_yzzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 197);

            auto ta_yzzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 198);

            auto ta_yzzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 199);

            auto ta_yzzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 200);

            auto ta_yzzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 201);

            auto ta_yzzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 202);

            auto ta_yzzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 203);

            auto ta_yzzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 204);

            auto ta_yzzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 205);

            auto ta_yzzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 206);

            auto ta_yzzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 207);

            auto ta_yzzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 208);

            auto ta_yzzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 209);

            auto ta_zzzz_xxxx_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 210);

            auto ta_zzzz_xxxy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 211);

            auto ta_zzzz_xxxz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 212);

            auto ta_zzzz_xxyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 213);

            auto ta_zzzz_xxyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 214);

            auto ta_zzzz_xxzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 215);

            auto ta_zzzz_xyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 216);

            auto ta_zzzz_xyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 217);

            auto ta_zzzz_xyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 218);

            auto ta_zzzz_xzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 219);

            auto ta_zzzz_yyyy_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 220);

            auto ta_zzzz_yyyz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 221);

            auto ta_zzzz_yyzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 222);

            auto ta_zzzz_yzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 223);

            auto ta_zzzz_zzzz_0 = primBuffer.data(pidx_a_4_4_m0 + 225 * idx + 224);

            // Batch of Integrals (180,225)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_yyzz_xxxx_0, ta_yyzz_xxxy_0, ta_yyzz_xxxz_0, \
                                         ta_yyzz_xxyy_0, ta_yyzz_xxyz_0, ta_yyzz_xxzz_0, ta_yyzz_xyyy_0, ta_yyzz_xyyz_0, \
                                         ta_yyzz_xyzz_0, ta_yyzz_xzzz_0, ta_yyzz_yyyy_0, ta_yyzz_yyyz_0, ta_yyzz_yyzz_0, \
                                         ta_yyzz_yzzz_0, ta_yyzz_zzzz_0, ta_yzz_xxx_0, ta_yzz_xxx_1, ta_yzz_xxxx_0, \
                                         ta_yzz_xxxx_1, ta_yzz_xxxy_0, ta_yzz_xxxy_1, ta_yzz_xxxz_0, ta_yzz_xxxz_1, \
                                         ta_yzz_xxy_0, ta_yzz_xxy_1, ta_yzz_xxyy_0, ta_yzz_xxyy_1, ta_yzz_xxyz_0, \
                                         ta_yzz_xxyz_1, ta_yzz_xxz_0, ta_yzz_xxz_1, ta_yzz_xxzz_0, ta_yzz_xxzz_1, \
                                         ta_yzz_xyy_0, ta_yzz_xyy_1, ta_yzz_xyyy_0, ta_yzz_xyyy_1, ta_yzz_xyyz_0, \
                                         ta_yzz_xyyz_1, ta_yzz_xyz_0, ta_yzz_xyz_1, ta_yzz_xyzz_0, ta_yzz_xyzz_1, \
                                         ta_yzz_xzz_0, ta_yzz_xzz_1, ta_yzz_xzzz_0, ta_yzz_xzzz_1, ta_yzz_yyy_0, \
                                         ta_yzz_yyy_1, ta_yzz_yyyy_0, ta_yzz_yyyy_1, ta_yzz_yyyz_0, ta_yzz_yyyz_1, \
                                         ta_yzz_yyz_0, ta_yzz_yyz_1, ta_yzz_yyzz_0, ta_yzz_yyzz_1, ta_yzz_yzz_0, \
                                         ta_yzz_yzz_1, ta_yzz_yzzz_0, ta_yzz_yzzz_1, ta_yzz_zzz_0, ta_yzz_zzz_1, \
                                         ta_yzz_zzzz_0, ta_yzz_zzzz_1, ta_yzzz_xxxx_0, ta_yzzz_xxxy_0, ta_yzzz_xxxz_0, \
                                         ta_yzzz_xxyy_0, ta_yzzz_xxyz_0, ta_yzzz_xxzz_0, ta_yzzz_xyyy_0, ta_yzzz_xyyz_0, \
                                         ta_yzzz_xyzz_0, ta_yzzz_xzzz_0, ta_yzzz_yyyy_0, ta_yzzz_yyyz_0, ta_yzzz_yyzz_0, \
                                         ta_yzzz_yzzz_0, ta_yzzz_zzzz_0, ta_zz_xxxx_0, ta_zz_xxxx_1, ta_zz_xxxy_0, \
                                         ta_zz_xxxy_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyz_0, \
                                         ta_zz_xxyz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xyyy_0, ta_zz_xyyy_1, ta_zz_xyyz_0, \
                                         ta_zz_xyyz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_yyyy_0, \
                                         ta_zz_yyyy_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yzzz_0, \
                                         ta_zz_yzzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1, ta_zzz_xxx_0, ta_zzz_xxx_1, \
                                         ta_zzz_xxxx_0, ta_zzz_xxxx_1, ta_zzz_xxxy_0, ta_zzz_xxxy_1, ta_zzz_xxxz_0, \
                                         ta_zzz_xxxz_1, ta_zzz_xxy_0, ta_zzz_xxy_1, ta_zzz_xxyy_0, ta_zzz_xxyy_1, \
                                         ta_zzz_xxyz_0, ta_zzz_xxyz_1, ta_zzz_xxz_0, ta_zzz_xxz_1, ta_zzz_xxzz_0, \
                                         ta_zzz_xxzz_1, ta_zzz_xyy_0, ta_zzz_xyy_1, ta_zzz_xyyy_0, ta_zzz_xyyy_1, \
                                         ta_zzz_xyyz_0, ta_zzz_xyyz_1, ta_zzz_xyz_0, ta_zzz_xyz_1, ta_zzz_xyzz_0, \
                                         ta_zzz_xyzz_1, ta_zzz_xzz_0, ta_zzz_xzz_1, ta_zzz_xzzz_0, ta_zzz_xzzz_1, \
                                         ta_zzz_yyy_0, ta_zzz_yyy_1, ta_zzz_yyyy_0, ta_zzz_yyyy_1, ta_zzz_yyyz_0, \
                                         ta_zzz_yyyz_1, ta_zzz_yyz_0, ta_zzz_yyz_1, ta_zzz_yyzz_0, ta_zzz_yyzz_1, \
                                         ta_zzz_yzz_0, ta_zzz_yzz_1, ta_zzz_yzzz_0, ta_zzz_yzzz_1, ta_zzz_zzz_0, \
                                         ta_zzz_zzz_1, ta_zzz_zzzz_0, ta_zzz_zzzz_1, ta_zzzz_xxxx_0, ta_zzzz_xxxy_0, \
                                         ta_zzzz_xxxz_0, ta_zzzz_xxyy_0, ta_zzzz_xxyz_0, ta_zzzz_xxzz_0, ta_zzzz_xyyy_0, \
                                         ta_zzzz_xyyz_0, ta_zzzz_xyzz_0, ta_zzzz_xzzz_0, ta_zzzz_yyyy_0, ta_zzzz_yyyz_0, \
                                         ta_zzzz_yyzz_0, ta_zzzz_yzzz_0, ta_zzzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_yyzz_xxxx_0[j] =
                    pa_y[j] * ta_yzz_xxxx_0[j] - pc_y[j] * ta_yzz_xxxx_1[j] + 0.5 * fl1_fx * ta_zz_xxxx_0[j] - 0.5 * fl1_fx * ta_zz_xxxx_1[j];

                ta_yyzz_xxxy_0[j] = pa_y[j] * ta_yzz_xxxy_0[j] - pc_y[j] * ta_yzz_xxxy_1[j] + 0.5 * fl1_fx * ta_zz_xxxy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxxy_1[j] + 0.5 * fl1_fx * ta_yzz_xxx_0[j] - 0.5 * fl1_fx * ta_yzz_xxx_1[j];

                ta_yyzz_xxxz_0[j] =
                    pa_y[j] * ta_yzz_xxxz_0[j] - pc_y[j] * ta_yzz_xxxz_1[j] + 0.5 * fl1_fx * ta_zz_xxxz_0[j] - 0.5 * fl1_fx * ta_zz_xxxz_1[j];

                ta_yyzz_xxyy_0[j] = pa_y[j] * ta_yzz_xxyy_0[j] - pc_y[j] * ta_yzz_xxyy_1[j] + 0.5 * fl1_fx * ta_zz_xxyy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxyy_1[j] + fl1_fx * ta_yzz_xxy_0[j] - fl1_fx * ta_yzz_xxy_1[j];

                ta_yyzz_xxyz_0[j] = pa_y[j] * ta_yzz_xxyz_0[j] - pc_y[j] * ta_yzz_xxyz_1[j] + 0.5 * fl1_fx * ta_zz_xxyz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xxyz_1[j] + 0.5 * fl1_fx * ta_yzz_xxz_0[j] - 0.5 * fl1_fx * ta_yzz_xxz_1[j];

                ta_yyzz_xxzz_0[j] =
                    pa_y[j] * ta_yzz_xxzz_0[j] - pc_y[j] * ta_yzz_xxzz_1[j] + 0.5 * fl1_fx * ta_zz_xxzz_0[j] - 0.5 * fl1_fx * ta_zz_xxzz_1[j];

                ta_yyzz_xyyy_0[j] = pa_y[j] * ta_yzz_xyyy_0[j] - pc_y[j] * ta_yzz_xyyy_1[j] + 0.5 * fl1_fx * ta_zz_xyyy_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyyy_1[j] + 1.5 * fl1_fx * ta_yzz_xyy_0[j] - 1.5 * fl1_fx * ta_yzz_xyy_1[j];

                ta_yyzz_xyyz_0[j] = pa_y[j] * ta_yzz_xyyz_0[j] - pc_y[j] * ta_yzz_xyyz_1[j] + 0.5 * fl1_fx * ta_zz_xyyz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyyz_1[j] + fl1_fx * ta_yzz_xyz_0[j] - fl1_fx * ta_yzz_xyz_1[j];

                ta_yyzz_xyzz_0[j] = pa_y[j] * ta_yzz_xyzz_0[j] - pc_y[j] * ta_yzz_xyzz_1[j] + 0.5 * fl1_fx * ta_zz_xyzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_xyzz_1[j] + 0.5 * fl1_fx * ta_yzz_xzz_0[j] - 0.5 * fl1_fx * ta_yzz_xzz_1[j];

                ta_yyzz_xzzz_0[j] =
                    pa_y[j] * ta_yzz_xzzz_0[j] - pc_y[j] * ta_yzz_xzzz_1[j] + 0.5 * fl1_fx * ta_zz_xzzz_0[j] - 0.5 * fl1_fx * ta_zz_xzzz_1[j];

                ta_yyzz_yyyy_0[j] = pa_y[j] * ta_yzz_yyyy_0[j] - pc_y[j] * ta_yzz_yyyy_1[j] + 0.5 * fl1_fx * ta_zz_yyyy_0[j] -
                                    0.5 * fl1_fx * ta_zz_yyyy_1[j] + 2.0 * fl1_fx * ta_yzz_yyy_0[j] - 2.0 * fl1_fx * ta_yzz_yyy_1[j];

                ta_yyzz_yyyz_0[j] = pa_y[j] * ta_yzz_yyyz_0[j] - pc_y[j] * ta_yzz_yyyz_1[j] + 0.5 * fl1_fx * ta_zz_yyyz_0[j] -
                                    0.5 * fl1_fx * ta_zz_yyyz_1[j] + 1.5 * fl1_fx * ta_yzz_yyz_0[j] - 1.5 * fl1_fx * ta_yzz_yyz_1[j];

                ta_yyzz_yyzz_0[j] = pa_y[j] * ta_yzz_yyzz_0[j] - pc_y[j] * ta_yzz_yyzz_1[j] + 0.5 * fl1_fx * ta_zz_yyzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_yyzz_1[j] + fl1_fx * ta_yzz_yzz_0[j] - fl1_fx * ta_yzz_yzz_1[j];

                ta_yyzz_yzzz_0[j] = pa_y[j] * ta_yzz_yzzz_0[j] - pc_y[j] * ta_yzz_yzzz_1[j] + 0.5 * fl1_fx * ta_zz_yzzz_0[j] -
                                    0.5 * fl1_fx * ta_zz_yzzz_1[j] + 0.5 * fl1_fx * ta_yzz_zzz_0[j] - 0.5 * fl1_fx * ta_yzz_zzz_1[j];

                ta_yyzz_zzzz_0[j] =
                    pa_y[j] * ta_yzz_zzzz_0[j] - pc_y[j] * ta_yzz_zzzz_1[j] + 0.5 * fl1_fx * ta_zz_zzzz_0[j] - 0.5 * fl1_fx * ta_zz_zzzz_1[j];

                ta_yzzz_xxxx_0[j] = pa_y[j] * ta_zzz_xxxx_0[j] - pc_y[j] * ta_zzz_xxxx_1[j];

                ta_yzzz_xxxy_0[j] =
                    pa_y[j] * ta_zzz_xxxy_0[j] - pc_y[j] * ta_zzz_xxxy_1[j] + 0.5 * fl1_fx * ta_zzz_xxx_0[j] - 0.5 * fl1_fx * ta_zzz_xxx_1[j];

                ta_yzzz_xxxz_0[j] = pa_y[j] * ta_zzz_xxxz_0[j] - pc_y[j] * ta_zzz_xxxz_1[j];

                ta_yzzz_xxyy_0[j] = pa_y[j] * ta_zzz_xxyy_0[j] - pc_y[j] * ta_zzz_xxyy_1[j] + fl1_fx * ta_zzz_xxy_0[j] - fl1_fx * ta_zzz_xxy_1[j];

                ta_yzzz_xxyz_0[j] =
                    pa_y[j] * ta_zzz_xxyz_0[j] - pc_y[j] * ta_zzz_xxyz_1[j] + 0.5 * fl1_fx * ta_zzz_xxz_0[j] - 0.5 * fl1_fx * ta_zzz_xxz_1[j];

                ta_yzzz_xxzz_0[j] = pa_y[j] * ta_zzz_xxzz_0[j] - pc_y[j] * ta_zzz_xxzz_1[j];

                ta_yzzz_xyyy_0[j] =
                    pa_y[j] * ta_zzz_xyyy_0[j] - pc_y[j] * ta_zzz_xyyy_1[j] + 1.5 * fl1_fx * ta_zzz_xyy_0[j] - 1.5 * fl1_fx * ta_zzz_xyy_1[j];

                ta_yzzz_xyyz_0[j] = pa_y[j] * ta_zzz_xyyz_0[j] - pc_y[j] * ta_zzz_xyyz_1[j] + fl1_fx * ta_zzz_xyz_0[j] - fl1_fx * ta_zzz_xyz_1[j];

                ta_yzzz_xyzz_0[j] =
                    pa_y[j] * ta_zzz_xyzz_0[j] - pc_y[j] * ta_zzz_xyzz_1[j] + 0.5 * fl1_fx * ta_zzz_xzz_0[j] - 0.5 * fl1_fx * ta_zzz_xzz_1[j];

                ta_yzzz_xzzz_0[j] = pa_y[j] * ta_zzz_xzzz_0[j] - pc_y[j] * ta_zzz_xzzz_1[j];

                ta_yzzz_yyyy_0[j] =
                    pa_y[j] * ta_zzz_yyyy_0[j] - pc_y[j] * ta_zzz_yyyy_1[j] + 2.0 * fl1_fx * ta_zzz_yyy_0[j] - 2.0 * fl1_fx * ta_zzz_yyy_1[j];

                ta_yzzz_yyyz_0[j] =
                    pa_y[j] * ta_zzz_yyyz_0[j] - pc_y[j] * ta_zzz_yyyz_1[j] + 1.5 * fl1_fx * ta_zzz_yyz_0[j] - 1.5 * fl1_fx * ta_zzz_yyz_1[j];

                ta_yzzz_yyzz_0[j] = pa_y[j] * ta_zzz_yyzz_0[j] - pc_y[j] * ta_zzz_yyzz_1[j] + fl1_fx * ta_zzz_yzz_0[j] - fl1_fx * ta_zzz_yzz_1[j];

                ta_yzzz_yzzz_0[j] =
                    pa_y[j] * ta_zzz_yzzz_0[j] - pc_y[j] * ta_zzz_yzzz_1[j] + 0.5 * fl1_fx * ta_zzz_zzz_0[j] - 0.5 * fl1_fx * ta_zzz_zzz_1[j];

                ta_yzzz_zzzz_0[j] = pa_y[j] * ta_zzz_zzzz_0[j] - pc_y[j] * ta_zzz_zzzz_1[j];

                ta_zzzz_xxxx_0[j] =
                    pa_z[j] * ta_zzz_xxxx_0[j] - pc_z[j] * ta_zzz_xxxx_1[j] + 1.5 * fl1_fx * ta_zz_xxxx_0[j] - 1.5 * fl1_fx * ta_zz_xxxx_1[j];

                ta_zzzz_xxxy_0[j] =
                    pa_z[j] * ta_zzz_xxxy_0[j] - pc_z[j] * ta_zzz_xxxy_1[j] + 1.5 * fl1_fx * ta_zz_xxxy_0[j] - 1.5 * fl1_fx * ta_zz_xxxy_1[j];

                ta_zzzz_xxxz_0[j] = pa_z[j] * ta_zzz_xxxz_0[j] - pc_z[j] * ta_zzz_xxxz_1[j] + 1.5 * fl1_fx * ta_zz_xxxz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xxxz_1[j] + 0.5 * fl1_fx * ta_zzz_xxx_0[j] - 0.5 * fl1_fx * ta_zzz_xxx_1[j];

                ta_zzzz_xxyy_0[j] =
                    pa_z[j] * ta_zzz_xxyy_0[j] - pc_z[j] * ta_zzz_xxyy_1[j] + 1.5 * fl1_fx * ta_zz_xxyy_0[j] - 1.5 * fl1_fx * ta_zz_xxyy_1[j];

                ta_zzzz_xxyz_0[j] = pa_z[j] * ta_zzz_xxyz_0[j] - pc_z[j] * ta_zzz_xxyz_1[j] + 1.5 * fl1_fx * ta_zz_xxyz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xxyz_1[j] + 0.5 * fl1_fx * ta_zzz_xxy_0[j] - 0.5 * fl1_fx * ta_zzz_xxy_1[j];

                ta_zzzz_xxzz_0[j] = pa_z[j] * ta_zzz_xxzz_0[j] - pc_z[j] * ta_zzz_xxzz_1[j] + 1.5 * fl1_fx * ta_zz_xxzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xxzz_1[j] + fl1_fx * ta_zzz_xxz_0[j] - fl1_fx * ta_zzz_xxz_1[j];

                ta_zzzz_xyyy_0[j] =
                    pa_z[j] * ta_zzz_xyyy_0[j] - pc_z[j] * ta_zzz_xyyy_1[j] + 1.5 * fl1_fx * ta_zz_xyyy_0[j] - 1.5 * fl1_fx * ta_zz_xyyy_1[j];

                ta_zzzz_xyyz_0[j] = pa_z[j] * ta_zzz_xyyz_0[j] - pc_z[j] * ta_zzz_xyyz_1[j] + 1.5 * fl1_fx * ta_zz_xyyz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xyyz_1[j] + 0.5 * fl1_fx * ta_zzz_xyy_0[j] - 0.5 * fl1_fx * ta_zzz_xyy_1[j];

                ta_zzzz_xyzz_0[j] = pa_z[j] * ta_zzz_xyzz_0[j] - pc_z[j] * ta_zzz_xyzz_1[j] + 1.5 * fl1_fx * ta_zz_xyzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xyzz_1[j] + fl1_fx * ta_zzz_xyz_0[j] - fl1_fx * ta_zzz_xyz_1[j];

                ta_zzzz_xzzz_0[j] = pa_z[j] * ta_zzz_xzzz_0[j] - pc_z[j] * ta_zzz_xzzz_1[j] + 1.5 * fl1_fx * ta_zz_xzzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_xzzz_1[j] + 1.5 * fl1_fx * ta_zzz_xzz_0[j] - 1.5 * fl1_fx * ta_zzz_xzz_1[j];

                ta_zzzz_yyyy_0[j] =
                    pa_z[j] * ta_zzz_yyyy_0[j] - pc_z[j] * ta_zzz_yyyy_1[j] + 1.5 * fl1_fx * ta_zz_yyyy_0[j] - 1.5 * fl1_fx * ta_zz_yyyy_1[j];

                ta_zzzz_yyyz_0[j] = pa_z[j] * ta_zzz_yyyz_0[j] - pc_z[j] * ta_zzz_yyyz_1[j] + 1.5 * fl1_fx * ta_zz_yyyz_0[j] -
                                    1.5 * fl1_fx * ta_zz_yyyz_1[j] + 0.5 * fl1_fx * ta_zzz_yyy_0[j] - 0.5 * fl1_fx * ta_zzz_yyy_1[j];

                ta_zzzz_yyzz_0[j] = pa_z[j] * ta_zzz_yyzz_0[j] - pc_z[j] * ta_zzz_yyzz_1[j] + 1.5 * fl1_fx * ta_zz_yyzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_yyzz_1[j] + fl1_fx * ta_zzz_yyz_0[j] - fl1_fx * ta_zzz_yyz_1[j];

                ta_zzzz_yzzz_0[j] = pa_z[j] * ta_zzz_yzzz_0[j] - pc_z[j] * ta_zzz_yzzz_1[j] + 1.5 * fl1_fx * ta_zz_yzzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_yzzz_1[j] + 1.5 * fl1_fx * ta_zzz_yzz_0[j] - 1.5 * fl1_fx * ta_zzz_yzz_1[j];

                ta_zzzz_zzzz_0[j] = pa_z[j] * ta_zzz_zzzz_0[j] - pc_z[j] * ta_zzz_zzzz_1[j] + 1.5 * fl1_fx * ta_zz_zzzz_0[j] -
                                    1.5 * fl1_fx * ta_zz_zzzz_1[j] + 2.0 * fl1_fx * ta_zzz_zzz_0[j] - 2.0 * fl1_fx * ta_zzz_zzz_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
